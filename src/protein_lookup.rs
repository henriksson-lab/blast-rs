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
/// Uses a CSR (compressed sparse row) layout for cache-friendly access:
/// `data[offsets[hash]..offsets[hash+1]]` holds query offsets whose
/// neighborhood contains the word that hashes to `hash`.
/// `pv` is a presence-vector bit array for fast rejection.
#[allow(dead_code)]
/// Maximum hits stored inline per backbone cell (matches NCBI AA_HITS_PER_CELL).
const HITS_PER_CELL: usize = 3;

/// One cell of the thick backbone. Stores up to HITS_PER_CELL query offsets
/// inline (no pointer chase). Overflows to a separate array.
#[derive(Clone, Copy)]
struct BackboneCell {
    num_used: u16,
    /// If num_used <= HITS_PER_CELL: inline entries.
    /// If num_used > HITS_PER_CELL: entries[0] is overflow cursor into `overflow`.
    entries: [i32; HITS_PER_CELL],
}

pub struct ProteinLookupTable {
    word_size: usize,
    alphabet_size: usize,
    /// Thick backbone: each cell stores up to 3 hits inline.
    backbone: Vec<BackboneCell>,
    /// Overflow array for cells with > HITS_PER_CELL hits.
    overflow: Vec<i32>,
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
        // Table size uses power-of-2 charsize (NCBI approach) for fast shift-based hashing
        let table_size = 1usize << (word_size * CHARSIZE);
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
            // Pre-allocate buffers reused across all query positions.
            let mut suffix_max = vec![0i32; word_size + 1];
            let mut word_buf = vec![0u8; word_size];

            for i in 0..=(query.len() - word_size) {
                let query_word = &query[i..i + word_size];

                // Compute max possible score at each suffix position for pruning.
                // suffix_max[k] = sum of row_max for positions k..word_size-1
                suffix_max[word_size] = 0;
                for k in (0..word_size).rev() {
                    suffix_max[k] =
                        suffix_max[k + 1] + row_max[query_word[k] as usize];
                }

                // Recursive enumeration with pruning.
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
        for (idx, entries) in backbone.iter().enumerate() {
            if !entries.is_empty() {
                pv[idx >> 6] |= 1u64 << (idx & 63);
            }
        }

        // Convert to thick backbone (inline ≤3 hits, overflow for rest).
        // Matches NCBI BLAST+ AaLookupBackboneCell layout.
        let empty_cell = BackboneCell { num_used: 0, entries: [0; HITS_PER_CELL] };
        let mut thick: Vec<BackboneCell> = vec![empty_cell; table_size];
        let mut overflow: Vec<i32> = Vec::new();

        for (idx, entries) in backbone.iter().enumerate() {
            let n = entries.len();
            if n == 0 {
                continue;
            }
            thick[idx].num_used = n as u16;
            if n <= HITS_PER_CELL {
                // Store inline — no pointer chase needed at scan time
                for (i, &val) in entries.iter().enumerate() {
                    thick[idx].entries[i] = val;
                }
            } else {
                // Store overflow cursor in entries[0]
                thick[idx].entries[0] = overflow.len() as i32;
                overflow.extend_from_slice(entries);
            }
        }

        ProteinLookupTable {
            word_size,
            alphabet_size,
            backbone: thick,
            overflow,
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

/// Bits per residue for hashing (ceil(log2(alphabet_size))).
/// AA_SIZE=28, charsize=5 (rounds up to 32). Matches NCBI BLAST+ `charsize`.
const CHARSIZE: usize = 5;
/// Mask for the hash: (1 << (word_size * CHARSIZE)) - 1.
/// For word_size=3: mask = (1 << 15) - 1 = 32767.
const HASH_MASK_W3: usize = (1 << (3 * CHARSIZE)) - 1;

/// Hash a word using NCBI-style shift+or (matching ComputeTableIndex).
/// hash = (w[0] << (n-1)*charsize) | (w[1] << (n-2)*charsize) | ... | w[n-1]
#[inline]
fn word_hash(word: &[u8], _alphabet_size: usize) -> usize {
    let mut h: usize = 0;
    for &b in word {
        h = (h << CHARSIZE) | b as usize;
    }
    h
}

/// Incremental hash update (NCBI ComputeTableIndexIncremental).
#[inline]
fn word_hash_incremental(prev_hash: usize, new_char: u8, mask: usize) -> usize {
    ((prev_hash << CHARSIZE) | new_char as usize) & mask
}

/// Build a merged presence vector from multiple lookup tables.
///
/// The merged PV is the bitwise OR of all individual PVs. A subject word
/// only needs to be checked against individual queries if the merged PV bit
/// is set — this skips ~90% of positions when queries are diverse.
pub fn merge_pv(tables: &[&ProteinLookupTable]) -> Vec<u64> {
    if tables.is_empty() {
        return Vec::new();
    }
    let pv_len = tables[0].pv.len();
    let mut merged = vec![0u64; pv_len];
    for table in tables {
        for (i, &bits) in table.pv.iter().enumerate() {
            merged[i] |= bits;
        }
    }
    merged
}

/// Batch scan: scan ONE subject against MULTIPLE queries using a merged PV.
///
/// For each subject position, the merged PV is checked first (one branch).
/// Only on a hit, individual query tables are checked. This reduces the
/// inner loop from O(positions × queries) to O(positions + hits × queries).
///
/// Returns a Vec of (query_index, Vec<ProteinHit>) for queries that had hits.
pub fn batch_scan_subject(
    queries: &[&[u8]],
    tables: &[&ProteinLookupTable],
    merged_pv: &[u64],
    subject: &[u8],
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
    x_dropoff: i32,
) -> Vec<(usize, Vec<ProteinHit>)> {
    if tables.is_empty() || subject.is_empty() {
        return Vec::new();
    }
    let word_size = tables[0].word_size;
    let alphabet_size = tables[0].alphabet_size;

    if subject.len() < word_size {
        return Vec::new();
    }

    // Per-query results and diagonal tracking
    let num_queries = queries.len();
    let mut results: Vec<Vec<ProteinHit>> = vec![Vec::new(); num_queries];

    // Diagonal tracking per query — lazily initialized
    let diag_count_max = subject.len() + queries.iter().map(|q| q.len()).max().unwrap_or(0);
    let mut diag_bufs: Vec<Vec<i32>> = (0..num_queries)
        .map(|_| vec![-1i32; diag_count_max])
        .collect();

    // Slide over the subject — ONE pass
    for s_pos in 0..=(subject.len() - word_size) {
        let s_word = &subject[s_pos..s_pos + word_size];
        let hash = word_hash(s_word, alphabet_size);

        // Merged PV check — rejects ~90% of positions in one branch
        if merged_pv[hash >> 6] & (1u64 << (hash & 63)) == 0 {
            continue;
        }

        // This position matched the merged PV — check individual queries
        for qi in 0..num_queries {
            let table = tables[qi];
            let query = queries[qi];

            // Individual PV check
            if table.pv[hash >> 6] & (1u64 << (hash & 63)) == 0 {
                continue;
            }

            // Thick backbone lookup
            let cell = &table.backbone[hash];
            let num = cell.num_used as usize;
            if num == 0 { continue; }
            let hit_slice: &[i32] = if num <= HITS_PER_CELL {
                &cell.entries[..num]
            } else {
                let cursor = cell.entries[0] as usize;
                &table.overflow[cursor..cursor + num]
            };

            for &q_off in hit_slice {
                let q_pos = q_off as usize;

                // Diagonal tracking
                let diag = s_pos + query.len() - q_pos;
                if diag >= diag_bufs[qi].len() { continue; }
                let last = diag_bufs[qi][diag];
                if last >= 0 && (s_pos as i32) < last + word_size as i32 {
                    continue;
                }

                // Ungapped extension
                if let Some((qs, qe, ss, se, score, ident)) =
                    protein_ungapped_extend(query, subject, q_pos, s_pos, matrix, x_dropoff)
                {
                    diag_bufs[qi][diag] = se as i32;
                    let alen = (qe - qs) as i32;
                    results[qi].push(ProteinHit {
                        query_start: qs, query_end: qe,
                        subject_start: ss, subject_end: se,
                        score, num_ident: ident,
                        align_length: alen, mismatches: alen - ident,
                        gap_opens: 0, qseq: None, sseq: None,
                    });
                }
            }
        }
    }

    // Collect non-empty results with query index
    results.into_iter().enumerate()
        .filter(|(_, hits)| !hits.is_empty())
        .map(|(qi, mut hits)| {
            hits.sort_by(|a, b| b.score.cmp(&a.score));
            (qi, hits)
        })
        .collect()
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
    protein_scan_with_table(query, subject, matrix, &table, x_dropoff)
}

/// Scan a subject using a pre-built lookup table (avoids rebuilding per subject).
pub fn protein_scan_with_table(
    query: &[u8],
    subject: &[u8],
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
    table: &ProteinLookupTable,
    x_dropoff: i32,
) -> Vec<ProteinHit> {
    let mut diag_buf = Vec::new();
    protein_scan_with_table_reuse(query, subject, matrix, table, x_dropoff, &mut diag_buf)
}

/// Two-hit window size (matches NCBI BLAST+ default for protein).
const TWO_HIT_WINDOW: i32 = 40;

/// Extend left from position (q_start-1, s_start-1) with x-dropoff.
/// Returns (best_score, left_displacement, num_identities).
/// Matches NCBI s_BlastAaExtendLeft.
#[inline]
fn extend_left(
    query: &[u8], subject: &[u8],
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
    q_start: usize, s_start: usize,
    x_dropoff: i32,
) -> (i32, i32, i32) {
    let mut score = 0i32;
    let mut best = 0i32;
    let mut best_d = 0i32;
    let mut ident = 0i32;
    let mut best_ident = 0i32;
    if q_start == 0 || s_start == 0 { return (0, 0, 0); }
    let mut qi = q_start - 1;
    let mut si = s_start - 1;
    let mut d = 1i32;
    loop {
        unsafe {
            let q = *query.as_ptr().add(qi);
            let s = *subject.as_ptr().add(si);
            score += *matrix.get_unchecked(q as usize).get_unchecked(s as usize);
            if q == s { ident += 1; }
        }
        if score > best { best = score; best_d = d; best_ident = ident; }
        if best - score > x_dropoff { break; }
        if qi == 0 || si == 0 { break; }
        qi -= 1;
        si -= 1;
        d += 1;
    }
    (best, best_d, best_ident)
}

/// Extend right from position (q_start, s_start) with x-dropoff.
/// `init_score` is the cumulative score from left extension.
/// Returns (best_score, right_displacement, num_identities).
/// Matches NCBI s_BlastAaExtendRight.
#[inline]
fn extend_right(
    query: &[u8], subject: &[u8],
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
    q_start: usize, s_start: usize,
    x_dropoff: i32, init_score: i32,
) -> (i32, i32, i32) {
    let mut score = init_score;
    let mut best = init_score;
    let mut best_d = 0i32;
    let mut ident = 0i32;
    let mut best_ident = 0i32;
    let mut qi = q_start;
    let mut si = s_start;
    let qlen = query.len();
    let slen = subject.len();
    while qi < qlen && si < slen {
        unsafe {
            let q = *query.as_ptr().add(qi);
            let s = *subject.as_ptr().add(si);
            score += *matrix.get_unchecked(q as usize).get_unchecked(s as usize);
            if q == s { ident += 1; }
        }
        if score > best { best = score; best_d = (qi - q_start + 1) as i32; best_ident = ident; }
        if best - score > x_dropoff { break; }
        qi += 1;
        si += 1;
    }
    (best, best_d, best_ident)
}

/// Like `protein_scan_with_table` but reuses a diagonal tracking buffer.
///
/// Uses the NCBI BLAST+ two-hit algorithm: a word hit on a diagonal only
/// triggers ungapped extension if a second hit was seen on the same diagonal
/// within `TWO_HIT_WINDOW` positions. This dramatically reduces extensions.
pub fn protein_scan_with_table_reuse(
    query: &[u8],
    subject: &[u8],
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
    table: &ProteinLookupTable,
    x_dropoff: i32,
    diag_buf: &mut Vec<i32>,
) -> Vec<ProteinHit> {
    let word_size = table.word_size;
    if query.len() < word_size || subject.len() < word_size {
        return Vec::new();
    }

    let diag_count = query.len() + subject.len();
    diag_buf.clear();
    diag_buf.resize(diag_count, i32::MIN);
    let diag_array = diag_buf;

    let mut hits: Vec<ProteinHit> = Vec::new();
    let ws = word_size as i32;
    let qlen = query.len();
    let mask = (1usize << (word_size * CHARSIZE)) - 1;

    let last_pos = subject.len() - word_size;
    let mut hash = word_hash(&subject[0..word_size], 0);

    // SAFETY: All indices are bounds-checked by loop structure:
    // - s_pos in 0..=last_pos where last_pos = subject.len() - word_size
    // - hash < table_size (= alphabet_size^word_size), PV/backbone sized to table_size
    // - diag < diag_count = query.len() + subject.len()
    // Removing bounds checks in this hot loop matches NCBI C approach.
    let pv = table.pv.as_ptr();
    let backbone = table.backbone.as_ptr();
    let overflow = table.overflow.as_ptr();
    let subj = subject.as_ptr();
    let diag_ptr = diag_array.as_mut_ptr();

    for s_pos in 0..=last_pos {
        if s_pos > 0 {
            unsafe {
                // NCBI ComputeTableIndexIncremental: shift + or + mask (no multiply)
                let new = *subj.add(s_pos + word_size - 1) as usize;
                hash = ((hash << CHARSIZE) | new) & mask;
            }
        }

        unsafe {
            if *pv.add(hash >> 6) & (1u64 << (hash & 63)) == 0 {
                continue;
            }

            let cell = &*backbone.add(hash);
            let num = cell.num_used as usize;
            if num == 0 { continue; }

            let (hit_ptr, hit_len) = if num <= HITS_PER_CELL {
                (cell.entries.as_ptr(), num)
            } else {
                let cursor = cell.entries[0] as usize;
                (overflow.add(cursor), num)
            };

            let s_off = s_pos as i32;

            for i in 0..hit_len {
                let q_pos = *hit_ptr.add(i) as usize;
                let diag = s_pos + qlen - q_pos;
                let last = *diag_ptr.add(diag);

                if last == i32::MIN {
                    *diag_ptr.add(diag) = s_off;
                    continue;
                }
                if last < 0 {
                    if s_off < -(last + 1) { continue; }
                    *diag_ptr.add(diag) = s_off;
                    continue;
                }
                let diff = s_off - last;
                if diff >= TWO_HIT_WINDOW { *diag_ptr.add(diag) = s_off; continue; }
                if diff < ws { continue; }

                // NCBI two-hit extension (s_BlastAaExtendTwoHit):
                // 1. Find best starting point within the word at s_pos
                // 2. Extend LEFT — must reach the first hit (at `last`)
                // 3. Extend RIGHT only if left reached far enough
                let s_left_off = (last + ws) as usize; // end of first hit
                let s_right_off = s_pos;
                let q_right_off = q_pos;

                // Find best start within the word
                let mut wscore = 0i32;
                let mut best_wscore = 0i32;
                let mut right_d = 0usize;
                for k in 0..word_size {
                    let qi = q_right_off + k;
                    let si = s_right_off + k;
                    if qi < qlen && si < subject.len() {
                        wscore += *matrix.get_unchecked(
                            *query.as_ptr().add(qi) as usize
                        ).get_unchecked(
                            *subject.as_ptr().add(si) as usize
                        );
                    }
                    if wscore > best_wscore {
                        best_wscore = wscore;
                        right_d = k + 1;
                    }
                }
                let ext_q = q_right_off + right_d;
                let ext_s = s_right_off + right_d;

                // Extend left from ext_s-1 — must reach s_left_off
                let (left_score, left_d, left_ident) = extend_left(
                    query, subject, matrix, ext_q, ext_s, x_dropoff);

                let reached_first = left_d >= (ext_s as i32 - s_left_off as i32);

                if reached_first {
                    // Extend right with cumulative score
                    let (right_score, right_d_r, right_ident) = extend_right(
                        query, subject, matrix, ext_q, ext_s, x_dropoff, left_score);

                    let total_score = left_score.max(right_score);
                    if total_score > 0 {
                        let qs = ext_q - left_d as usize;
                        let qe = ext_q + right_d_r as usize;
                        let ss = ext_s - left_d as usize;
                        let se = ext_s + right_d_r as usize;
                        *diag_ptr.add(diag) = -(se as i32 + 1);
                        let alen = (qe - qs) as i32;
                        let ident = left_ident + right_ident;
                        hits.push(ProteinHit {
                            query_start: qs, query_end: qe,
                            subject_start: ss, subject_end: se,
                            score: total_score, num_ident: ident,
                            align_length: alen, mismatches: alen - ident,
                            gap_opens: 0, qseq: None, sseq: None,
                        });
                        continue;
                    }
                }
                *diag_ptr.add(diag) = s_off;
            }
        }
    }

    if hits.len() > 1 {
        hits.sort_unstable_by(|a, b| b.score.cmp(&a.score));
    }
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
    let table = ProteinLookupTable::build(query, word_size, matrix, threshold);
    protein_gapped_scan_with_table(
        query, subject, matrix, &table, ungap_x_dropoff,
        gap_open, gap_extend, gap_x_dropoff, ungap_cutoff,
    )
}

/// Gapped scan using a pre-built lookup table.
pub fn protein_gapped_scan_with_table(
    query: &[u8],
    subject: &[u8],
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
    table: &ProteinLookupTable,
    ungap_x_dropoff: i32,
    gap_open: i32,
    gap_extend: i32,
    gap_x_dropoff: i32,
    ungap_cutoff: i32,
) -> Vec<ProteinHit> {
    // Phase 1: ungapped seeds
    let ungapped = protein_scan_with_table(query, subject, matrix, table, ungap_x_dropoff);

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
        // Shift-based hash: h = (h << 5) | b for each byte
        // [1,2,3] = (1 << 10) | (2 << 5) | 3 = 1024 + 64 + 3 = 1091
        assert_eq!(word_hash(&[1, 2, 3], 28), (1 << 10) | (2 << 5) | 3);
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
        let cell = &table.backbone[hash];
        assert!(cell.num_used > 0, "Cell should have hits");
        let hits: &[i32] = if (cell.num_used as usize) <= HITS_PER_CELL {
            &cell.entries[..cell.num_used as usize]
        } else {
            let c = cell.entries[0] as usize;
            &table.overflow[c..c + cell.num_used as usize]
        };
        assert!(
            hits.contains(&0),
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

    #[test]
    fn test_lookup_table_pv_array_exact_match() {
        let m = simple_matrix();
        // Query: [1, 2, 3, 4, 5], word_size=3, threshold=12 (exact match only)
        let query = vec![1u8, 2, 3, 4, 5];
        let table = ProteinLookupTable::build(&query, 3, &m, 12.0);
        // Should have entries for words [1,2,3], [2,3,4], [3,4,5]
        let hashes = [
            word_hash(&[1, 2, 3], AA_SIZE),
            word_hash(&[2, 3, 4], AA_SIZE),
            word_hash(&[3, 4, 5], AA_SIZE),
        ];
        for (i, &hash) in hashes.iter().enumerate() {
            let cell = &table.backbone[hash];
            assert!(cell.num_used > 0, "Word at position {} should have entries", i);
            let hits: &[i32] = if (cell.num_used as usize) <= HITS_PER_CELL {
                &cell.entries[..cell.num_used as usize]
            } else {
                let c = cell.entries[0] as usize;
                &table.overflow[c..c + cell.num_used as usize]
            };
            assert!(hits.contains(&(i as i32)),
                "Word at position {} should map to query offset {}", i, i);
            // PV bit should be set
            assert_ne!(table.pv[hash >> 6] & (1u64 << (hash & 63)), 0,
                "PV bit should be set for word at position {}", i);
        }
    }

    #[test]
    fn test_lookup_table_no_neighbors_high_threshold() {
        let m = simple_matrix();
        // threshold=13 > max possible score (4*3=12), so NO words should match
        let query = vec![1u8, 2, 3, 4, 5];
        let table = ProteinLookupTable::build(&query, 3, &m, 13.0);
        let total_entries: usize = table.backbone.iter()
            .map(|c| c.num_used as usize)
            .sum();
        assert_eq!(total_entries, 0, "No words should meet threshold > max score");
    }

    #[test]
    fn test_lookup_table_neighbors_low_threshold() {
        let m = simple_matrix();
        // threshold=7 allows words with 2 matches + 1 mismatch (4+4-1=7)
        let query = vec![1u8, 2, 3];
        let table = ProteinLookupTable::build(&query, 3, &m, 7.0);
        let total_entries: usize = table.backbone.iter()
            .map(|c| c.num_used as usize)
            .sum();
        assert!(total_entries > 1,
            "Low threshold should produce neighboring words (got {})", total_entries);
    }

    #[test]
    fn test_protein_scan_with_table_reuse() {
        let m = crate::matrix::BLOSUM62;
        let query = vec![1u8, 5, 8, 12, 17, 1, 5, 8]; // AGIMV AGI in NCBIstdaa
        let table = ProteinLookupTable::build(&query, 3, &m, 11.0);
        // Scan two different subjects with same table
        let subject1 = query.clone();
        let subject2 = vec![2u8, 3, 4, 5, 6, 7]; // completely different
        let hits1 = protein_scan_with_table(&query, &subject1, &m, &table, 40);
        let hits2 = protein_scan_with_table(&query, &subject2, &m, &table, 40);
        assert!(!hits1.is_empty(), "Should find hits for identical sequences");
        // hits2 may or may not be empty depending on BLOSUM62 neighborhood
    }

    #[test]
    fn test_word_hash_deterministic() {
        // Same input should always produce same hash
        let h1 = word_hash(&[5, 10, 15], AA_SIZE);
        let h2 = word_hash(&[5, 10, 15], AA_SIZE);
        assert_eq!(h1, h2);
        // Different input should (generally) produce different hash
        let h3 = word_hash(&[5, 10, 16], AA_SIZE);
        assert_ne!(h1, h3);
    }

    #[test]
    fn test_word_hash_range() {
        // All hashes should be in [0, 2^(word_size*CHARSIZE))
        let max_hash = 1usize << (3 * CHARSIZE); // 32768
        for a in 0..AA_SIZE {
            for b in 0..AA_SIZE {
                for c in 0..AA_SIZE {
                    let h = word_hash(&[a as u8, b as u8, c as u8], AA_SIZE);
                    assert!(h < max_hash, "Hash {} out of range for [{},{},{}]", h, a, b, c);
                }
            }
        }
    }

    #[test]
    fn test_protein_gapped_scan_with_table_finds_alignment() {
        let m = crate::matrix::BLOSUM62;
        // Two related protein sequences
        let query = vec![1u8, 5, 8, 12, 17, 1, 5, 8, 12, 17]; // repeat of AGIMV
        let subject = query.clone();
        let table = ProteinLookupTable::build(&query, 3, &m, 11.0);
        let hits = protein_gapped_scan_with_table(
            &query, &subject, &m, &table,
            40,     // ungap_x_dropoff
            11, 1,  // gap_open, gap_extend
            260,    // gap_x_dropoff
            40,     // ungap_cutoff
        );
        assert!(!hits.is_empty(), "Gapped scan should find alignment for identical sequences");
        let best = &hits[0];
        assert!(best.score > 0);
        assert_eq!(best.query_start, 0);
    }
}
