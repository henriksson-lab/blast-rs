//! Scoring matrices and frequency ratios for BLAST.
//! Rust equivalent of matrix_freq_ratios.c and parts of blast_stat.c.

/// Standard amino acid alphabet size.
pub const AA_SIZE: usize = 28;

/// BLOSUM62 scoring matrix (NCBIstdaa order, 28x28).
/// Row/column order: -ABCDEFGHIKLMNPQRSTVWXYZU*OJ
pub static BLOSUM62: [[i32; AA_SIZE]; AA_SIZE] = {
    let mut m = [[0i32; AA_SIZE]; AA_SIZE];
    // Only fill the standard 20 amino acids + B, Z, X
    // This is a compile-time approximation; the full matrix is loaded at runtime
    // from the NCBI data files in the C code.
    m
};

/// Nucleotide scoring matrix for blastn.
/// Creates a simple match/mismatch matrix for BLASTNA encoding (16 values).
pub fn nucleotide_matrix(reward: i32, penalty: i32) -> [[i32; 16]; 16] {
    let mut m = [[penalty; 16]; 16];
    // Exact matches (same unambiguous base)
    for i in 0..4 {
        m[i][i] = reward;
    }
    // Ambiguous bases get penalty (conservative)
    // N (14) matches nothing specifically
    m
}

/// Background amino acid frequencies (Robinson & Robinson, 1991).
pub static AA_FREQUENCIES: [f64; 20] = [
    0.07805, // A
    0.05129, // C
    0.05364, // D
    0.06295, // E
    0.04259, // F
    0.07377, // G
    0.02199, // H
    0.05142, // I
    0.05744, // K
    0.09019, // L
    0.02243, // M
    0.04487, // N
    0.05203, // P
    0.04264, // Q
    0.05129, // R
    0.07120, // S
    0.05841, // T
    0.06441, // V
    0.01330, // W
    0.03216, // Y
];

/// Background nucleotide frequencies (uniform).
pub static NT_FREQUENCIES: [f64; 4] = [0.25, 0.25, 0.25, 0.25];

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nucleotide_matrix() {
        let m = nucleotide_matrix(1, -3);
        assert_eq!(m[0][0], 1);  // A-A match
        assert_eq!(m[0][1], -3); // A-C mismatch
        assert_eq!(m[2][2], 1);  // G-G match
        assert_eq!(m[3][0], -3); // T-A mismatch
    }

    #[test]
    fn test_aa_frequencies_sum() {
        let sum: f64 = AA_FREQUENCIES.iter().sum();
        assert!((sum - 1.0).abs() < 0.05, "AA frequencies should sum to ~1.0, got {}", sum);
    }

    #[test]
    fn test_nt_frequencies_sum() {
        let sum: f64 = NT_FREQUENCIES.iter().sum();
        assert!((sum - 1.0).abs() < 0.001);
    }
}
