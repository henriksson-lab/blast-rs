//! Position-Specific Scoring Matrix (PSSM) for PSI-BLAST.
//! A PSSM represents the amino acid preferences at each position
//! of a multiple sequence alignment, used for iterative search.

use crate::matrix::AA_SIZE;

/// A Position-Specific Scoring Matrix.
#[derive(Debug, Clone)]
pub struct Pssm {
    /// Scoring matrix: `pssm[position][amino_acid]` = score
    pub scores: Vec<[i32; AA_SIZE]>,
    /// Query length (number of positions)
    pub length: usize,
    /// Information content per position
    pub info_content: Vec<f64>,
}

impl Pssm {
    /// Create a PSSM from a query sequence and a standard scoring matrix.
    /// This is the initial PSSM before any iteration (equivalent to the matrix itself).
    pub fn from_sequence(query: &[u8], matrix: &[[i32; AA_SIZE]; AA_SIZE]) -> Self {
        let length = query.len();
        let mut scores = Vec::with_capacity(length);
        let info_content = vec![0.0; length];

        for &aa in query {
            let aa = aa as usize;
            if aa < AA_SIZE {
                scores.push(matrix[aa]);
            } else {
                scores.push([0; AA_SIZE]);
            }
        }

        Pssm { scores, length, info_content }
    }

    /// Score a subject amino acid at a given position.
    #[inline]
    pub fn score_at(&self, position: usize, aa: u8) -> i32 {
        if position < self.length && (aa as usize) < AA_SIZE {
            self.scores[position][aa as usize]
        } else {
            -4 // default for unknown
        }
    }

    /// Update the PSSM from a set of aligned sequences (PSI-BLAST iteration).
    /// Takes aligned sequences and recomputes position-specific scores.
    pub fn update_from_alignment(
        &mut self,
        aligned_seqs: &[Vec<u8>],
        bg_freqs: &[f64; 20],
        pseudocount_weight: f64,
    ) {
        if aligned_seqs.is_empty() { return; }

        for pos in 0..self.length {
            // Count amino acid frequencies at this position
            let mut counts = [0.0f64; 20];
            let mut total = 0.0;
            for seq in aligned_seqs {
                if pos < seq.len() {
                    let aa = seq[pos] as usize;
                    if aa > 0 && aa <= 20 {
                        counts[aa - 1] += 1.0;
                        total += 1.0;
                    }
                }
            }

            if total == 0.0 { continue; }

            // Compute position-specific frequencies with pseudocounts
            let alpha = total / (total + pseudocount_weight);
            for aa in 0..20 {
                let observed = counts[aa] / total;
                let freq = alpha * observed + (1.0 - alpha) * bg_freqs[aa];
                // Convert to log-odds score
                if freq > 0.0 && bg_freqs[aa] > 0.0 {
                    let log_odds = (freq / bg_freqs[aa]).ln() / 0.3176; // scale by lambda
                    self.scores[pos][aa + 1] = (log_odds * 2.0).round() as i32; // scale factor
                }
            }
        }
    }
}

/// Run one iteration of PSI-BLAST.
/// Takes a query, a set of subject sequences, and the current PSSM.
/// Returns hits that pass the inclusion threshold.
pub fn psi_blast_iteration(
    pssm: &Pssm,
    subjects: &[(String, Vec<u8>)], // (id, NCBIstdaa sequence)
    inclusion_evalue: f64,
    search_space: f64,
    kbp_lambda: f64,
    kbp_k: f64,
) -> Vec<(String, i32, f64)> { // (subject_id, score, evalue)
    let mut results = Vec::new();

    for (subj_id, subj_seq) in subjects {
        if subj_seq.len() < 3 || pssm.length < 3 { continue; }

        let mut best_score = 0i32;
        // Simple scan: try each subject position
        for si in 0..=(subj_seq.len().saturating_sub(pssm.length)) {
            let mut score = 0i32;
            let len = pssm.length.min(subj_seq.len() - si);
            for k in 0..len {
                score += pssm.score_at(k, subj_seq[si + k]);
            }
            best_score = best_score.max(score);
        }

        if best_score > 0 {
            let evalue = kbp_k * search_space * (-kbp_lambda * best_score as f64).exp();
            if evalue <= inclusion_evalue {
                results.push((subj_id.clone(), best_score, evalue));
            }
        }
    }

    results.sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal));
    results
}

#[cfg(test)]
mod tests {
    use super::*;

    fn simple_matrix() -> [[i32; AA_SIZE]; AA_SIZE] {
        let mut m = [[0i32; AA_SIZE]; AA_SIZE];
        for i in 1..21 { m[i][i] = 5; }
        for i in 1..21 { for j in 1..21 { if i != j { m[i][j] = -2; } } }
        m
    }

    #[test]
    fn test_pssm_from_sequence() {
        let query = vec![1u8, 2, 3, 4, 5]; // A, B, C, D, E
        let matrix = simple_matrix();
        let pssm = Pssm::from_sequence(&query, &matrix);
        assert_eq!(pssm.length, 5);
        assert_eq!(pssm.score_at(0, 1), 5); // A at position 0 = match
        assert_eq!(pssm.score_at(0, 2), -2); // B at position 0 = mismatch
    }

    #[test]
    fn test_psi_blast_iteration() {
        let query = vec![1u8, 2, 3, 4, 5];
        let matrix = simple_matrix();
        let pssm = Pssm::from_sequence(&query, &matrix);

        let subjects = vec![
            ("match".to_string(), vec![1u8, 2, 3, 4, 5]),
            ("mismatch".to_string(), vec![10u8, 11, 12, 13, 14]),
        ];

        let results = psi_blast_iteration(&pssm, &subjects, 1.0, 1000.0, 0.3176, 0.134);
        assert!(!results.is_empty(), "Should find matching subject");
        assert_eq!(results[0].0, "match");
    }
}
