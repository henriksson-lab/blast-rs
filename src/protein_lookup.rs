//! Protein word lookup table and scan for blastp/blastx/tblastn.
//!
//! Replaces the O(n^2) brute-force approach with the standard BLAST
//! neighborhood-word lookup table: for each w-mer in the query, all
//! amino-acid words scoring >= threshold (against the scoring matrix)
//! are hashed into a table keyed by word index. During scanning, each
//! subject w-mer is hashed once, checked against a presence vector for
//! fast rejection, and only on a PV hit are backbone entries examined
//! and ungapped extensions triggered.

use crate::matrix::AA_SIZE;
use crate::protein::{protein_ungapped_extend, protein_gapped_align, ncbistdaa_to_char};

/// Result of a protein hit after extension.
#[derive(Debug, Clone)]
pub struct ProteinHit {
    pub query_start: usize,
    pub query_end: usize,
    pub subject_start: usize,
    pub subject_end: usize,
    pub score: i32,
    pub num_ident: i32,
    pub align_length: i32,
    pub mismatches: i32,
    pub gap_opens: i32,
    pub qseq: Option<String>,
    pub sseq: Option<String>,
}

/// Protein word lookup table.
///
/// `backbone[hash]` holds the list of query offsets whose neighborhood
/// contains the word that hashes to `hash`.  `pv` is a presence-vector
/// bit array for fast rejection: if the bit for a hash is clear, the
/// backbone entry is guaranteed empty.
#[allow(dead_code)]
pub struct ProteinLookupTable {
    word_size: usize,
    alphabet_size: usize,
    backbone: Vec<Vec<i32>>,
    pv: Vec<u64>,
}

impl ProteinLookupTable {
    /// Build the lookup table from a query sequence.
    ///
    /// * `query` - NCBIstdaa-encoded query
    /// * `word_size` - typically 3 for blastp
    /// * `matrix` - scoring matrix (AA_SIZE x AA_SIZE)
    /// * `threshold` - minimum neighborhood word score (typically 11)
    pub fn build(
        query: &[u8],
        word_size: usize,
        matrix: &[[i32; AA_SIZE]; AA_SIZE],
        threshold: f64,
    ) -> Self {
        let alphabet_size = AA_SIZE;
        let table_size = alphabet_size.pow(word_size as u32);
        let mut backbone: Vec<Vec<i32>> = vec![Vec::new(); table_size];

        // Precompute per-row maximums for branch-and-bound pruning.
        // row_max[aa] = max score achievable when the query letter is `aa`.
        let mut row_max = [0i32; AA_SIZE];
        for q in 0..AA_SIZE {
            let mut mx = i32::MIN;
            for s in 0..AA_SIZE {
                if matrix[q][s] > mx {
                    mx = matrix[q][s];
                }
            }
            row_max[q] = mx;
        }

        // For each query position, generate all neighboring words and insert.
        let thresh_i = threshold as i32; // integer threshold for comparison
        if query.len() >= word_size {
            for i in 0..=(query.len() - word_size) {
                let query_word = &query[i..i + word_size];

                // Compute max possible score at each suffix position for pruning.
                // suffix_max[k] = sum of row_max for positions k..word_size-1
                let mut suffix_max = vec![0i32; word_size + 1];
                for k in (0..word_size).rev() {
                    suffix_max[k] =
                        suffix_max[k + 1] + row_max[query_word[k] as usize];
                }

                // Recursive enumeration with pruning.
                let mut word_buf = vec![0u8; word_size];
                enumerate_neighbors(
                    query_word,
                    matrix,
                    alphabet_size,
                    word_size,
                    thresh_i,
                    &suffix_max,
                    &mut word_buf,
                    0,
                    0,
                    i as i32,
                    &mut backbone,
                );
            }
        }

        // Build presence vector.
        let pv_len = (table_size + 63) / 64;
        let mut pv = vec![0u64; pv_len];
        for (idx, offsets) in backbone.iter().enumerate() {
            if !offsets.is_empty() {
                pv[idx >> 6] |= 1u64 << (idx & 63);
            }
        }

        ProteinLookupTable {
            word_size,
            alphabet_size,
            backbone,
            pv,
        }
    }
}

/// Recursively enumerate neighboring words, pruning branches where the
/// maximum attainable score falls below `threshold`.
fn enumerate_neighbors(
    query_word: &[u8],
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
    alphabet_size: usize,
    word_size: usize,
    threshold: i32,
    suffix_max: &[i32],
    word_buf: &mut [u8],
    pos: usize,
    score_so_far: i32,
    query_offset: i32,
    backbone: &mut [Vec<i32>],
) {
    if pos == word_size {
        // Compute hash and insert.
        let hash = word_hash(word_buf, alphabet_size);
        backbone[hash].push(query_offset);
        return;
    }

    let q_letter = query_word[pos] as usize;

    for aa in 0..alphabet_size {
        let s = score_so_far + matrix[q_letter][aa];
        // Prune: best possible score from remaining positions.
        if s + suffix_max[pos + 1] < threshold {
            continue;
        }
        word_buf[pos] = aa as u8;
        enumerate_neighbors(
            query_word,
            matrix,
            alphabet_size,
            word_size,
            threshold,
            suffix_max,
            word_buf,
            pos + 1,
            s,
            query_offset,
            backbone,
        );
    }
}

/// Hash a word of length `word_size` with alphabet of size `alphabet_size`.
/// hash = w[0]*alphabet_size^(n-1) + w[1]*alphabet_size^(n-2) + ... + w[n-1]
#[inline]
fn word_hash(word: &[u8], alphabet_size: usize) -> usize {
    let mut h: usize = 0;
    for &b in word {
        h = h * alphabet_size + b as usize;
    }
    h
}

/// Scan a subject sequence against a query using the protein lookup table,
/// performing ungapped extensions for each hit.
///
/// Returns a list of `ProteinHit` sorted by descending score.
pub fn protein_scan(
    query: &[u8],
    subject: &[u8],
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
    word_size: usize,
    threshold: f64,
    x_dropoff: i32,
) -> Vec<ProteinHit> {
    if query.len() < word_size || subject.len() < word_size {
        return Vec::new();
    }

    let table = ProteinLookupTable::build(query, word_size, matrix, threshold);

    // Diagonal tracking to avoid redundant extensions.
    // diagonal = s_pos - q_pos (shifted by query.len() to keep non-negative)
    let diag_count = query.len() + subject.len();
    let mut last_diag_hit: Vec<i32> = vec![-1; diag_count];

    let mut hits: Vec<ProteinHit> = Vec::new();

    // Slide over the subject.
    for s_pos in 0..=(subject.len() - word_size) {
        let s_word = &subject[s_pos..s_pos + word_size];

        // Hash the subject word.
        let hash = word_hash(s_word, table.alphabet_size);

        // PV check.
        if table.pv[hash >> 6] & (1u64 << (hash & 63)) == 0 {
            continue;
        }

        // Check backbone entries.
        for &q_off in &table.backbone[hash] {
            let q_pos = q_off as usize;

            // Diagonal tracking: skip if we already extended on this diagonal
            // at or past this subject position.
            let diag = s_pos + query.len() - q_pos;
            let last = last_diag_hit[diag];
            if last >= 0 && (s_pos as i32) < last + word_size as i32 {
                continue;
            }

            // Ungapped extension.
            if let Some((qs, qe, ss, se, score)) =
                protein_ungapped_extend(query, subject, q_pos, s_pos, matrix, x_dropoff)
            {
                last_diag_hit[diag] = se as i32;
                let alen = (qe - qs) as i32;
                let mut ident = 0i32;
                for k in 0..alen as usize {
                    if qs + k < query.len() && ss + k < subject.len()
                        && query[qs + k] == subject[ss + k] { ident += 1; }
                }
                let qseq: String = query[qs..qe].iter().map(|&b| ncbistdaa_to_char(b)).collect();
                let sseq: String = subject[ss..se].iter().map(|&b| ncbistdaa_to_char(b)).collect();
                hits.push(ProteinHit {
                    query_start: qs,
                    query_end: qe,
                    subject_start: ss,
                    subject_end: se,
                    score,
                    num_ident: ident,
                    align_length: alen,
                    mismatches: alen - ident,
                    gap_opens: 0,
                    qseq: Some(qseq),
                    sseq: Some(sseq),
                });
            }
        }
    }

    // Sort by descending score.
    hits.sort_by(|a, b| b.score.cmp(&a.score));
    hits
}

/// Scan + gapped extension: find ungapped seeds, then perform gapped DP on top hits.
///
/// This is the standard two-phase BLAST approach:
/// 1. Find seeds via lookup table + ungapped extension
/// 2. For seeds above a cutoff, perform gapped alignment with X-dropoff DP
pub fn protein_gapped_scan(
    query: &[u8],
    subject: &[u8],
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
    word_size: usize,
    threshold: f64,
    ungap_x_dropoff: i32,
    gap_open: i32,
    gap_extend: i32,
    gap_x_dropoff: i32,
    ungap_cutoff: i32,
) -> Vec<ProteinHit> {
    // Phase 1: ungapped seeds
    let ungapped = protein_scan(query, subject, matrix, word_size, threshold, ungap_x_dropoff);

    // Phase 2: gapped extension on seeds above cutoff
    let mut gapped_hits = Vec::new();
    // Track which diagonals we've already done gapped extension on
    let diag_count = query.len() + subject.len();
    let mut gapped_diag: Vec<bool> = vec![false; diag_count];

    for uh in &ungapped {
        if uh.score < ungap_cutoff { continue; }

        // Use center of ungapped hit as seed
        let seed_q = (uh.query_start + uh.query_end) / 2;
        let seed_s = (uh.subject_start + uh.subject_end) / 2;
        let diag = seed_s + query.len() - seed_q;
        if diag < diag_count && gapped_diag[diag] { continue; }
        if diag < diag_count { gapped_diag[diag] = true; }

        if let Some(gr) = protein_gapped_align(
            query, subject, seed_q, seed_s,
            matrix, gap_open, gap_extend, gap_x_dropoff,
        ) {
            let q_slice = &query[gr.query_start..gr.query_end];
            let s_slice = &subject[gr.subject_start..gr.subject_end];
            let (qseq, sseq) = gr.edit_script.render_alignment(
                q_slice, s_slice, ncbistdaa_to_char,
            );
            gapped_hits.push(ProteinHit {
                query_start: gr.query_start,
                query_end: gr.query_end,
                subject_start: gr.subject_start,
                subject_end: gr.subject_end,
                score: gr.score,
                num_ident: gr.num_ident,
                align_length: gr.align_length,
                mismatches: gr.mismatches,
                gap_opens: gr.gap_opens,
                qseq: Some(qseq),
                sseq: Some(sseq),
            });
        }
    }

    // If no seeds passed the cutoff, fall back to ungapped hits
    if gapped_hits.is_empty() {
        return ungapped;
    }

    gapped_hits.sort_by(|a, b| b.score.cmp(&a.score));
    gapped_hits
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrix::AA_SIZE;

    fn simple_matrix() -> [[i32; AA_SIZE]; AA_SIZE] {
        let mut m = [[0i32; AA_SIZE]; AA_SIZE];
        for i in 1..21 {
            m[i][i] = 4;
        }
        for i in 1..21 {
            for j in 1..21 {
                if i != j {
                    m[i][j] = -1;
                }
            }
        }
        m
    }

    #[test]
    fn test_word_hash() {
        assert_eq!(word_hash(&[1, 2, 3], 28), 1 * 28 * 28 + 2 * 28 + 3);
        assert_eq!(word_hash(&[0, 0, 0], 28), 0);
        assert_eq!(word_hash(&[0, 0, 1], 28), 1);
    }

    #[test]
    fn test_lookup_table_build() {
        let m = simple_matrix();
        // Query: 3 amino acids, word_size=3, threshold=12 (exact match only with score 4*3=12)
        let query = vec![1u8, 2, 3];
        let table = ProteinLookupTable::build(&query, 3, &m, 12.0);
        let hash = word_hash(&[1, 2, 3], AA_SIZE);
        assert!(
            table.backbone[hash].contains(&0),
            "Exact match word should have query offset 0"
        );
        // PV bit should be set.
        assert_ne!(table.pv[hash >> 6] & (1u64 << (hash & 63)), 0);
    }

    #[test]
    fn test_protein_scan_identical() {
        let m = simple_matrix();
        let query = vec![1u8, 2, 3, 4, 5, 6, 7, 8];
        let subject = query.clone();
        let hits = protein_scan(&query, &subject, &m, 3, 11.0, 20);
        assert!(!hits.is_empty(), "Should find hits for identical sequences");
        // Best hit should cover full length.
        let best = &hits[0];
        assert_eq!(best.score, 32); // 8 * 4
    }

    #[test]
    fn test_protein_scan_no_match() {
        let m = simple_matrix();
        // Query all 1s, subject all 2s — word score = -1*3 = -3 < threshold 11
        let query = vec![1u8, 1, 1, 1, 1];
        let subject = vec![2u8, 2, 2, 2, 2];
        let hits = protein_scan(&query, &subject, &m, 3, 11.0, 20);
        assert!(hits.is_empty(), "Should find no hits for unrelated sequences");
    }

    #[test]
    fn test_protein_scan_short_sequences() {
        let m = simple_matrix();
        let query = vec![1u8, 2];
        let subject = vec![1u8, 2, 3];
        let hits = protein_scan(&query, &subject, &m, 3, 11.0, 20);
        assert!(hits.is_empty(), "Sequences shorter than word_size yield no hits");
    }
}
