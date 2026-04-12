//! Protein BLAST (blastp) support.
//! Implements scoring matrices and protein word finding.

use crate::encoding::NCBISTDAA_TO_AMINOACID;
use crate::gapinfo::{GapEditScript, GapAlignOpType};
use crate::matrix::AA_SIZE;

/// Convert an NCBIstdaa-encoded amino acid to its single-letter IUPAC character.
pub fn ncbistdaa_to_char(b: u8) -> char {
    if (b as usize) < NCBISTDAA_TO_AMINOACID.len() {
        NCBISTDAA_TO_AMINOACID[b as usize] as u8 as char
    } else {
        'X'
    }
}

/// Score two amino acids using a scoring matrix.
#[inline]
pub fn score_aa(matrix: &[[i32; AA_SIZE]; AA_SIZE], aa1: u8, aa2: u8) -> i32 {
    if (aa1 as usize) < AA_SIZE && (aa2 as usize) < AA_SIZE {
        matrix[aa1 as usize][aa2 as usize]
    } else {
        -4 // default mismatch for unknown residues
    }
}

/// Perform ungapped protein extension from a seed position.
pub fn protein_ungapped_extend(
    query: &[u8],    // NCBIstdaa encoded
    subject: &[u8],  // NCBIstdaa encoded
    q_seed: usize,
    s_seed: usize,
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
    x_dropoff: i32,
) -> Option<(usize, usize, usize, usize, i32)> {
    // Extend right
    let mut score = 0i32;
    let mut best = 0i32;
    let mut best_r = 0usize;
    let (mut qi, mut si) = (q_seed, s_seed);
    while qi < query.len() && si < subject.len() {
        score += score_aa(matrix, query[qi], subject[si]);
        if score > best { best = score; best_r = qi - q_seed + 1; }
        if best - score > x_dropoff { break; }
        qi += 1;
        si += 1;
    }

    // Extend left
    let mut sl = 0i32;
    let mut best_l = 0i32;
    let mut best_left = 0usize;
    if q_seed > 0 && s_seed > 0 {
        qi = q_seed - 1;
        si = s_seed - 1;
        loop {
            sl += score_aa(matrix, query[qi], subject[si]);
            if sl > best_l { best_l = sl; best_left = q_seed - qi; }
            if best_l - sl > x_dropoff { break; }
            if qi == 0 || si == 0 { break; }
            qi -= 1;
            si -= 1;
        }
    }

    let total = best + best_l;
    if total <= 0 { return None; }

    Some((
        q_seed - best_left,
        q_seed + best_r,
        s_seed - best_left,
        s_seed + best_r,
        total,
    ))
}

/// Find neighboring words for a protein word using a scoring matrix.
/// Returns all words that score above the threshold when compared to the query word.
pub fn find_neighboring_words(
    query_word: &[u8],
    word_size: usize,
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
    threshold: f64,
) -> Vec<Vec<u8>> {
    // For word_size=3 (standard blastp), enumerate all 20^3 = 8000 possible words
    // and keep those scoring above threshold
    let mut neighbors = Vec::new();
    let aa_count = 20u8; // standard amino acids

    if word_size == 3 {
        for a in 0..aa_count {
            for b in 0..aa_count {
                for c in 0..aa_count {
                    let s = score_aa(matrix, query_word[0], a + 1)
                        + score_aa(matrix, query_word[1], b + 1)
                        + score_aa(matrix, query_word[2], c + 1);
                    if s as f64 >= threshold {
                        neighbors.push(vec![a + 1, b + 1, c + 1]);
                    }
                }
            }
        }
    }
    neighbors
}

/// Result of a protein gapped alignment.
#[derive(Debug, Clone)]
pub struct ProteinGappedResult {
    pub query_start: usize,
    pub query_end: usize,
    pub subject_start: usize,
    pub subject_end: usize,
    pub score: i32,
    pub num_ident: i32,
    pub align_length: i32,
    pub mismatches: i32,
    pub gap_opens: i32,
    pub edit_script: GapEditScript,
}

const MININT: i32 = i32::MIN / 2;

struct GapDP {
    best: i32,
    best_gap: i32,
}

/// Score-only gapped extension in one direction (left or right from seed).
/// Returns the best score found.
fn protein_gapped_score_one_dir(
    query: &[u8], subject: &[u8],
    m: usize, n: usize,
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
    gap_oe: i32, gap_extend: i32,
    mut x_dropoff: i32,
    reverse: bool,
) -> (i32, usize, usize) {
    if x_dropoff < gap_oe { x_dropoff = gap_oe; }
    if m == 0 || n == 0 { return (0, 0, 0); }

    let max_band = ((x_dropoff / gap_extend.max(1)) as usize + 10).min(n + 1);
    let mut sa = Vec::with_capacity(max_band);
    sa.push(GapDP { best: 0, best_gap: -gap_oe });

    let mut score = -gap_oe;
    while sa.len() < max_band && score >= -x_dropoff {
        sa.push(GapDP { best: score, best_gap: score - gap_oe });
        score -= gap_extend;
    }
    let mut b_size = sa.len();

    let mut best_score = 0i32;
    let mut best_q = 0usize;
    let mut best_s = 0usize;
    let mut first_b = 0usize;

    for ai in 1..=m {
        let a_idx = if reverse { m - ai } else { ai };
        if a_idx >= query.len() { break; }
        let a_letter = query[a_idx] as usize;

        let mut sc = MININT;
        let mut sgr = MININT;
        let mut last_b = first_b;

        for bi in first_b..b_size {
            let b_idx = if reverse {
                n.checked_sub(1 + bi).unwrap_or(usize::MAX)
            } else { bi + 1 };
            if b_idx >= subject.len() { break; }
            let b_letter = subject[b_idx] as usize;

            let sgc = sa[bi].best_gap;
            let mat_score = if a_letter < AA_SIZE && b_letter < AA_SIZE {
                matrix[a_letter][b_letter]
            } else { -4 };
            let next_sc = sa[bi].best + mat_score;

            if sc < sgc { sc = sgc; }
            if sc < sgr { sc = sgr; }

            if best_score - sc > x_dropoff {
                if first_b == bi { first_b += 1; }
                sa[bi].best = MININT;
            } else {
                last_b = bi;
                if sc > best_score {
                    best_score = sc;
                    best_q = ai;
                    best_s = bi + 1;
                }
                sa[bi].best_gap = sc.max(sgc) - gap_oe;
                sgr = sc.max(sgr) - gap_oe;
                sa[bi].best = sc;
                sc = next_sc;
                sgr -= gap_extend;
                continue;
            }

            sa[bi].best_gap = MININT;
            sgr = MININT;
            sa[bi].best = MININT;
            sc = next_sc;
        }

        // Handle possible extension past current band
        if sc >= best_score - x_dropoff && last_b + 1 < max_band {
            let new_bi = last_b + 1;
            if new_bi >= sa.len() {
                sa.push(GapDP { best: MININT, best_gap: MININT });
            }
            if new_bi < b_size || b_size < max_band {
                sa[new_bi].best = sc;
                sa[new_bi].best_gap = sc - gap_oe;
                if sc > best_score {
                    best_score = sc;
                    best_q = ai;
                    best_s = new_bi + 1;
                }
                last_b = new_bi;
            }
        }

        if last_b < first_b { break; }
        b_size = last_b + 1;
    }

    (best_score, best_q, best_s)
}

/// Needleman-Wunsch affine-gap traceback on a sub-region.
/// Returns (edit_script, num_ident, mismatches, gap_opens, align_length).
/// Local alignment (Smith-Waterman) with traceback.
/// Returns (edit_script, score, q_start_in_slice, s_start_in_slice, q_end_in_slice, s_end_in_slice).
fn protein_nw_traceback(
    query: &[u8],
    subject: &[u8],
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
    gap_open: i32,
    gap_extend: i32,
) -> (GapEditScript, i32, usize, usize, usize, usize) {
    let m = query.len();
    let n = subject.len();
    if m == 0 || n == 0 {
        let mut esp = GapEditScript::new();
        if m > 0 { esp.push(GapAlignOpType::Ins, m as i32); }
        if n > 0 { esp.push(GapAlignOpType::Del, n as i32); }
        return (esp, 0, 0, 0, 0, 0);
    }

    let gap_oe = gap_open + gap_extend;

    // DP matrices: H = best score, E = gap-in-query (Del), F = gap-in-subject (Ins)
    // Traceback: 0 = Sub, 1 = Del (gap in query), 2 = Ins (gap in subject)
    let sz = (m + 1) * (n + 1);
    let mut h = vec![i32::MIN / 2; sz];
    let mut e = vec![i32::MIN / 2; sz]; // best score ending with gap in query
    let mut f = vec![i32::MIN / 2; sz]; // best score ending with gap in subject
    let mut tb = vec![0u8; sz];

    let idx = |i: usize, j: usize| -> usize { i * (n + 1) + j };

    // Free end gaps: first row/column initialized to 0 (semi-global)
    h[idx(0, 0)] = 0;
    // Leave h[0][j] = 0 and h[i][0] = 0 (no penalty for leading gaps)

    let mut best_score = 0i32;
    let mut best_i = 0usize;
    let mut best_j = 0usize;

    for i in 1..=m {
        let qi = query[i - 1] as usize;
        for j in 1..=n {
            let sj = subject[j - 1] as usize;

            // Del: gap in query (consume subject)
            let e1 = e[idx(i, j - 1)] - gap_extend;
            let e2 = h[idx(i, j - 1)] - gap_oe;
            e[idx(i, j)] = e1.max(e2);

            // Ins: gap in subject (consume query)
            let f1 = f[idx(i - 1, j)] - gap_extend;
            let f2 = h[idx(i - 1, j)] - gap_oe;
            f[idx(i, j)] = f1.max(f2);

            // Sub: aligned pair
            let mat = if qi < AA_SIZE && sj < AA_SIZE { matrix[qi][sj] } else { -4 };
            let diag = h[idx(i - 1, j - 1)] + mat;

            // Smith-Waterman: allow restart from 0
            let best = diag.max(e[idx(i, j)]).max(f[idx(i, j)]).max(0);
            h[idx(i, j)] = best;

            if best == 0 {
                tb[idx(i, j)] = 3; // local restart
            } else if best == diag {
                tb[idx(i, j)] = 0; // Sub
            } else if best == e[idx(i, j)] {
                tb[idx(i, j)] = 1; // Del
            } else {
                tb[idx(i, j)] = 2; // Ins
            }

            if best > best_score {
                best_score = best;
                best_i = i;
                best_j = j;
            }
        }
    }

    let nw_score = best_score;

    // Traceback from best score position (Smith-Waterman style)
    let mut ops: Vec<(GapAlignOpType, i32)> = Vec::new();
    let mut i = best_i;
    let mut j = best_j;
    while i > 0 && j > 0 && tb[idx(i, j)] != 3 {
        let op = match tb[idx(i, j)] {
            0 => { i -= 1; j -= 1; GapAlignOpType::Sub }
            1 => { j -= 1; GapAlignOpType::Del }
            _ => { i -= 1; GapAlignOpType::Ins }
        };
        if let Some(last) = ops.last_mut() {
            if last.0 == op {
                last.1 += 1;
                continue;
            }
        }
        ops.push((op, 1));
    }
    ops.reverse();

    // i, j are now the start of the local alignment (0-based in the slice)
    let sw_q_start = i;
    let sw_s_start = j;
    let sw_q_end = best_i;  // 1-based → end exclusive
    let sw_s_end = best_j;

    let mut esp = GapEditScript::new();
    for (op, count) in ops {
        esp.push(op, count);
    }
    (esp, nw_score, sw_q_start, sw_s_start, sw_q_end, sw_s_end)
}

/// Perform gapped protein alignment from a seed position using X-dropoff DP.
///
/// This is the protein equivalent of `blast_gapped_score_only` from traceback.rs,
/// using a substitution matrix (e.g., BLOSUM62) instead of match/mismatch rewards.
pub fn protein_gapped_align(
    query: &[u8], subject: &[u8],
    seed_q: usize, seed_s: usize,
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
    gap_open: i32, gap_extend: i32,
    x_dropoff: i32,
) -> Option<ProteinGappedResult> {
    let gap_oe = gap_open + gap_extend;

    // Left extension
    let (score_l, ext_q_l, ext_s_l) = protein_gapped_score_one_dir(
        &query[..seed_q + 1], &subject[..seed_s + 1],
        seed_q + 1, seed_s + 1,
        matrix, gap_oe, gap_extend, x_dropoff, true,
    );

    // Right extension
    let (score_r, ext_q_r, ext_s_r) = if seed_q < query.len() - 1 && seed_s < subject.len() - 1 {
        protein_gapped_score_one_dir(
            &query[seed_q..], &subject[seed_s..],
            query.len() - seed_q - 1, subject.len() - seed_s - 1,
            matrix, gap_oe, gap_extend, x_dropoff, false,
        )
    } else { (0, 0, 0) };

    let total_score = score_l + score_r;
    if total_score <= 0 { return None; }

    // X-drop boundaries with a small margin to let SW find the optimal local
    // alignment. Keep margin small to avoid O(m*n) blowup on long sequences.
    let margin = 20;
    let q_start = (seed_q + 1).saturating_sub(ext_q_l).saturating_sub(margin);
    let q_end = (seed_q + ext_q_r + margin).min(query.len());
    let s_start = (seed_s + 1).saturating_sub(ext_s_l).saturating_sub(margin);
    let s_end = (seed_s + ext_s_r + margin).min(subject.len());
    if q_start >= q_end || s_start >= s_end { return None; }

    // Smith-Waterman local alignment on the x-drop region — produces the
    // optimal local score regardless of seed position.
    let q_slice = &query[q_start..q_end];
    let s_slice = &subject[s_start..s_end];
    let (edit_script, sw_score, sw_qs, sw_ss, sw_qe, sw_se) =
        protein_nw_traceback(q_slice, s_slice, matrix, gap_open, gap_extend);
    let final_score = total_score.max(sw_score);
    if final_score <= 0 { return None; }

    // Adjust boundaries to the SW local alignment region
    let final_q_start = q_start + sw_qs;
    let final_q_end = q_start + sw_qe;
    let final_s_start = s_start + sw_ss;
    let final_s_end = s_start + sw_se;
    if final_q_start >= final_q_end || final_s_start >= final_s_end { return None; }

    let local_q = &query[final_q_start..final_q_end];
    let local_s = &subject[final_s_start..final_s_end];
    let (align_length, num_ident, gap_opens) = edit_script.count_identities(local_q, local_s);
    let mismatches = (align_length - num_ident - gap_opens).max(0);

    Some(ProteinGappedResult {
        query_start: final_q_start,
        query_end: final_q_end,
        subject_start: final_s_start,
        subject_end: final_s_end,
        score: final_score,
        num_ident,
        align_length,
        mismatches,
        gap_opens,
        edit_script,
    })
}

/// Combined ungapped extend + gapped alignment for a seed hit.
/// Returns a ProteinGappedResult if the hit passes extension.
pub fn protein_search_hit(
    query: &[u8], subject: &[u8],
    q_seed: usize, s_seed: usize,
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
    gap_open: i32, gap_extend: i32,
    x_dropoff: i32,
) -> Option<ProteinGappedResult> {
    // First do ungapped extension to filter
    let ug = protein_ungapped_extend(query, subject, q_seed, s_seed, matrix, x_dropoff)?;
    // If ungapped score is reasonable, do full gapped alignment
    let mid_q = (ug.0 + ug.1) / 2;
    let mid_s = (ug.2 + ug.3) / 2;
    protein_gapped_align(query, subject, mid_q, mid_s, matrix, gap_open, gap_extend, x_dropoff)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn simple_blosum62() -> [[i32; AA_SIZE]; AA_SIZE] {
        let mut m = [[0i32; AA_SIZE]; AA_SIZE];
        // Fill diagonal with positive scores (match)
        for i in 1..21 {
            m[i][i] = 4;
        }
        // Fill off-diagonal with negative (mismatch)
        for i in 1..21 {
            for j in 1..21 {
                if i != j { m[i][j] = -1; }
            }
        }
        m
    }

    #[test]
    fn test_score_aa() {
        let m = simple_blosum62();
        assert_eq!(score_aa(&m, 1, 1), 4); // A-A match
        assert_eq!(score_aa(&m, 1, 2), -1); // A-B mismatch
    }

    #[test]
    fn test_protein_extend() {
        let m = simple_blosum62();
        let query = vec![1u8, 2, 3, 4, 5, 6, 7, 8]; // 8 amino acids
        let subject = vec![1u8, 2, 3, 4, 5, 6, 7, 8]; // identical
        let result = protein_ungapped_extend(&query, &subject, 4, 4, &m, 20);
        assert!(result.is_some());
        let (qs, qe, _ss, _se, score) = result.unwrap();
        assert_eq!(score, 32); // 8 matches * 4
        assert_eq!(qe - qs, 8);
    }

    #[test]
    fn test_neighboring_words() {
        let m = simple_blosum62();
        let word = vec![1u8, 2, 3]; // A, B, C
        let neighbors = find_neighboring_words(&word, 3, &m, 11.0);
        // The exact match (1,2,3) scores 4+4+4=12 >= 11, so it should be included
        assert!(neighbors.iter().any(|w| w == &vec![1u8, 2, 3]),
            "Exact match should be a neighbor");
    }

    #[test]
    fn test_protein_gapped_align_produces_edit_script() {
        let m = simple_blosum62();
        let query   = vec![1u8, 2, 3, 4, 5, 6, 7, 8];
        let subject = vec![1u8, 2, 3, 4, 5, 6, 7, 8];
        let result = protein_gapped_align(&query, &subject, 4, 4, &m, 11, 1, 50);
        assert!(result.is_some());
        let r = result.unwrap();
        assert!(r.edit_script.alignment_length() > 0);
        let (qseq, sseq) = r.edit_script.render_alignment(
            &query[r.query_start..r.query_end],
            &subject[r.subject_start..r.subject_end],
            ncbistdaa_to_char,
        );
        assert!(!qseq.is_empty());
        assert!(!sseq.is_empty());
        // Aligned strings should have the same length (alignment property)
        assert_eq!(qseq.len(), sseq.len(),
            "aligned strings must have equal length: qseq={}, sseq={}", qseq, sseq);
    }

    #[test]
    fn test_protein_gapped_align_with_gap() {
        let m = simple_blosum62();
        // Test that NW traceback correctly handles different-length regions
        // (which naturally arise from X-dropoff boundary asymmetry)
        let query   = vec![1u8, 2, 3, 4, 5, 6, 7, 8, 9, 10];
        let subject = vec![1u8, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12];
        let result = protein_gapped_align(&query, &subject, 5, 5, &m, 11, 1, 50);
        assert!(result.is_some());
        let r = result.unwrap();
        let (qseq, sseq) = r.edit_script.render_alignment(
            &query[r.query_start..r.query_end],
            &subject[r.subject_start..r.subject_end],
            ncbistdaa_to_char,
        );
        assert_eq!(qseq.len(), sseq.len(),
            "aligned strings must have equal length: qseq={}, sseq={}", qseq, sseq);
        assert!(!qseq.is_empty());
    }

    #[test]
    fn test_nw_traceback_identical() {
        let m = simple_blosum62();
        let seq = vec![1u8, 2, 3, 4, 5];
        let (esp, _, _, _, _, _) = protein_nw_traceback(&seq, &seq, &m, 11, 1);
        let (qseq, sseq) = esp.render_alignment(&seq, &seq, ncbistdaa_to_char);
        assert_eq!(qseq, sseq);
        assert!(!qseq.contains('-'));
    }

    #[test]
    fn test_sw_traceback_local_alignment() {
        let m = simple_blosum62();
        // Short sequences: SW finds best local alignment
        let query   = vec![1u8, 2, 3, 4, 5, 6];
        let subject = vec![1u8, 2, 4, 5, 6]; // missing 3
        let (esp, score, qs, ss, qe, se) = protein_nw_traceback(&query, &subject, &m, 11, 1);
        assert!(score > 0, "Should have positive score");
        let q_slice = &query[qs..qe];
        let s_slice = &subject[ss..se];
        let (qseq, sseq) = esp.render_alignment(q_slice, s_slice, ncbistdaa_to_char);
        assert!(!qseq.is_empty(), "alignment should not be empty: qseq={}, sseq={}", qseq, sseq);
    }

    #[test]
    fn test_protein_ungapped_extend_xdrop() {
        // Use real BLOSUM62 matrix. Build sequences that match for a stretch
        // then diverge so the X-dropoff terminates extension.
        let matrix = crate::matrix::BLOSUM62;
        // NCBIstdaa: A=1, R=2, N=3, D=4, C=5
        // Query:   ARNDCARND (9 residues, all match subject for first 5, then differ)
        let query   = vec![1u8, 2, 3, 4, 5, 1, 2, 3, 4];
        // Subject: ARNDC + XXXXX (first 5 match, last 4 are very different: W=17)
        let subject = vec![1u8, 2, 3, 4, 5, 17, 17, 17, 17];

        let result = protein_ungapped_extend(&query, &subject, 2, 2, &matrix, 15);
        assert!(result.is_some());
        let (_qs, qe, _ss, _se, score) = result.unwrap();
        // The best region should be the first 5 matching residues.
        // The xdrop should prevent extending far into the mismatching tail.
        assert!(score > 0);
        // End of aligned region should not extend much past position 5.
        assert!(qe <= 7, "X-dropoff should terminate extension, qe={}", qe);
    }

    #[test]
    fn test_protein_ungapped_extend_short_sequence() {
        let matrix = crate::matrix::BLOSUM62;
        // Single residue sequences — should not crash.
        let query   = vec![1u8]; // A
        let subject = vec![1u8]; // A
        let result = protein_ungapped_extend(&query, &subject, 0, 0, &matrix, 20);
        assert!(result.is_some());
        let (_qs, _qe, _ss, _se, score) = result.unwrap();
        assert!(score > 0);

        // Two residues, seed at 0.
        let q2 = vec![1u8, 2]; // A, R
        let s2 = vec![1u8, 2]; // A, R
        let result2 = protein_ungapped_extend(&q2, &s2, 0, 0, &matrix, 20);
        assert!(result2.is_some());

        // Empty-like: query length 1, subject length 1, mismatching.
        let q3 = vec![5u8]; // C
        let s3 = vec![17u8]; // W
        let result3 = protein_ungapped_extend(&q3, &s3, 0, 0, &matrix, 5);
        // This should return None since the score is negative.
        assert!(result3.is_none());
    }

    #[test]
    fn test_protein_gapped_align_identical() {
        let matrix = crate::matrix::BLOSUM62;
        // Identical sequences should produce a perfect alignment with no gaps.
        let seq = vec![1u8, 2, 3, 4, 5, 6, 7, 8, 9, 10]; // A R N D C Q E G H I
        let result = protein_gapped_align(&seq, &seq, 5, 5, &matrix, 11, 1, 100);
        assert!(result.is_some());
        let r = result.unwrap();
        assert!(r.score > 0);
        assert_eq!(r.gap_opens, 0, "Identical sequences should have no gaps");
        assert_eq!(r.mismatches, 0, "Identical sequences should have no mismatches");
        // All residues should be identities.
        assert!(r.num_ident > 0);
    }

    #[test]
    fn test_protein_gapped_align_with_multiple_gaps() {
        use crate::encoding::AMINOACID_TO_NCBISTDAA;
        let matrix = crate::matrix::BLOSUM62;
        // Query has two insertions relative to subject.
        // Query:   ARNDCQ--EGHIKLM--NQRST
        // Subject: ARNDCQWWEGHIKLMPPNQRST
        // We encode so that alignment requires two gaps.
        let q_raw = b"ARNDCQEGHIKLMNQRST";
        let s_raw = b"ARNDCQWWEGHIKLMPPNQRST";
        let query: Vec<u8> = q_raw.iter().map(|&b| AMINOACID_TO_NCBISTDAA[b as usize & 0x7F]).collect();
        let subject: Vec<u8> = s_raw.iter().map(|&b| AMINOACID_TO_NCBISTDAA[b as usize & 0x7F]).collect();

        let result = protein_gapped_align(&query, &subject, 3, 3, &matrix, 11, 1, 100);
        assert!(result.is_some());
        let r = result.unwrap();
        assert!(r.score > 0);
        assert!(r.gap_opens >= 1,
            "Alignment with inserted residues should have gaps, got gap_opens={}", r.gap_opens);
    }

    #[test]
    fn test_protein_gapped_align_score_vs_ungapped() {
        use crate::encoding::AMINOACID_TO_NCBISTDAA;
        let matrix = crate::matrix::BLOSUM62;
        // Two sequences that are mostly similar but have a gap.
        // Gapped alignment should score >= ungapped for the same pair.
        let q_raw = b"ARNDCQEGHIKLMNPQRSTVWY";
        let s_raw = b"ARNDCQXXEGHIKLMNPQRSTVWY"; // two insertions after Q
        let query: Vec<u8> = q_raw.iter().map(|&b| AMINOACID_TO_NCBISTDAA[b as usize & 0x7F]).collect();
        let subject: Vec<u8> = s_raw.iter().map(|&b| AMINOACID_TO_NCBISTDAA[b as usize & 0x7F]).collect();

        // Ungapped extension
        let ug = protein_ungapped_extend(&query, &subject, 3, 3, &matrix, 100);
        let ug_score = ug.map_or(0, |(_,_,_,_,s)| s);

        // Gapped alignment
        let gapped = protein_gapped_align(&query, &subject, 3, 3, &matrix, 11, 1, 100);
        assert!(gapped.is_some());
        let g_score = gapped.unwrap().score;

        assert!(g_score >= ug_score,
            "Gapped score ({}) should be >= ungapped score ({})", g_score, ug_score);
    }

    #[test]
    fn test_gapped_align_srta_vs_p0dpq5() {
        use crate::encoding::AMINOACID_TO_NCBISTDAA;
        let matrix = crate::matrix::BLOSUM62;

        // srtA query vs P0DPQ5 sortase A (NCBI BLAST+ gets score 183)
        let q_raw = b"MIIRHPKKKRIMGKWIIAFWLLSAVGVLLLMPAEASVAKYQQNQQIAAIDRTGTAAETDSSLDVAKIELGDPVGILTIPSISLKLPIYDGTSDKILENGVGITEGTGDITGGNGKNPLIAGHSGLYKDNLFDDLPSVKKGEKFYIKVDGEQHAYQIDRIEEVQKDELQRNFVTYLEPNPNEDRVTLMTCTPKGINTHRFLVYGKRVTFTKSELKDEENKKQKLSWKWLLGSTVFLSVMIIGSLFVYKKKK";
        let s_raw = b"MNKQRIYSIVAILLFVVGGVLIGKPFYDGYQAEKKQTENVQAVQKMDYEKHETEFVDASKIDQPDLAEVANASLDKKQVIGRISIPSVSLELPVLKSSTEKNLLSGAATVKENQVMGKGNYALAGHNMSKKGVLFSDIASLKKGDKIYLYDNENEYEYAVTGVSEVTPDKWEVVEDHGKDEITLITCVSVKDNSKRYVVAGDLVGTKAKK";
        let query: Vec<u8> = q_raw.iter().map(|&b| AMINOACID_TO_NCBISTDAA[b as usize & 0x7F]).collect();
        let subject: Vec<u8> = s_raw.iter().map(|&b| AMINOACID_TO_NCBISTDAA[b as usize & 0x7F]).collect();

        // Gapped alignment from a seed in the middle of the alignment
        let result = protein_gapped_align(&query, &subject, 50, 45, &matrix, 11, 1, 260);
        assert!(result.is_some());
        let r = result.unwrap();
        assert!(r.score >= 150,
            "Gapped alignment should get score >= 150, got {}", r.score);
    }
}
