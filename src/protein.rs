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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrix::nucleotide_matrix;

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
        let (qs, qe, ss, se, score) = result.unwrap();
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
