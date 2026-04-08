//! Protein BLAST (blastp) support.
//! Implements scoring matrices and protein word finding.

use crate::matrix::AA_SIZE;

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

    let q_start = seed_q + 1 - ext_q_l;
    let q_end = seed_q + ext_q_r;
    let s_start = seed_s + 1 - ext_s_l;
    let s_end = seed_s + ext_s_r;

    // Compute identity statistics from the aligned region
    let alen = (q_end - q_start).max(s_end - s_start) as i32;
    let mut ident = 0i32;
    let mut mismatches = 0i32;
    let match_len = (q_end - q_start).min(s_end - s_start);
    for k in 0..match_len {
        if q_start + k < query.len() && s_start + k < subject.len() {
            if query[q_start + k] == subject[s_start + k] {
                ident += 1;
            } else {
                mismatches += 1;
            }
        }
    }
    let gap_opens = ((q_end - q_start) as i32 - (s_end - s_start) as i32).unsigned_abs() as i32;

    Some(ProteinGappedResult {
        query_start: q_start,
        query_end: q_end,
        subject_start: s_start,
        subject_end: s_end,
        score: total_score,
        num_ident: ident,
        align_length: alen.max(match_len as i32),
        mismatches,
        gap_opens,
    })
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
}
