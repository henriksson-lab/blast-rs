//! Rust equivalent of blast_stat.c — Karlin-Altschul statistics.
//! This is the mathematical core for computing E-values and bit scores.

use crate::math;

/// Karlin-Altschul statistical parameters for one context.
#[derive(Debug, Clone)]
pub struct KarlinBlk {
    pub lambda: f64,  // Lambda parameter
    pub k: f64,       // K parameter
    pub log_k: f64,   // ln(K)
    pub h: f64,       // H (relative entropy)
}

impl KarlinBlk {
    pub fn is_valid(&self) -> bool {
        self.lambda > 0.0 && self.k > 0.0 && self.h > 0.0
    }

    /// Convert a raw score to a bit score.
    pub fn raw_to_bit(&self, raw_score: i32) -> f64 {
        (self.lambda * raw_score as f64 - self.log_k) / std::f64::consts::LN_2
    }

    /// Convert a raw score and search space to an E-value.
    pub fn raw_to_evalue(&self, raw_score: i32, search_space: f64) -> f64 {
        search_space * self.k * (-self.lambda * raw_score as f64).exp()
    }

    /// Convert an E-value to the minimum raw score needed.
    pub fn evalue_to_raw(&self, evalue: f64, search_space: f64) -> i32 {
        let score = -(evalue / (search_space * self.k)).ln() / self.lambda;
        score.ceil() as i32
    }
}

/// Score frequency distribution.
#[derive(Debug, Clone)]
pub struct ScoreFreq {
    pub score_min: i32,
    pub score_max: i32,
    pub obs_min: i32,
    pub obs_max: i32,
    pub score_avg: f64,
    pub sprob: Vec<f64>, // probability for each score value
}

/// Gapped Karlin-Altschul parameters (precomputed tables).
#[derive(Debug, Clone)]
pub struct GappedParams {
    pub gap_open: i32,
    pub gap_extend: i32,
    pub lambda: f64,
    pub k: f64,
    pub h: f64,
    pub alpha: f64,
    pub beta: f64,
}

/// Precomputed gap parameter tables for nucleotide scoring.
/// Format: (gap_open, gap_extend, lambda, K, H, alpha, beta)
pub const BLASTN_PARAMS_1_3: &[(i32, i32, f64, f64, f64, f64, f64)] = &[
    // reward=1, penalty=-3
    (5, 2, 0.208, 0.049, 0.14, 1.24, -0.70),
    (2, 2, 0.278, 0.075, 0.27, 0.96, -0.52),
    (1, 2, 0.308, 0.084, 0.34, 0.86, -0.42),
    (0, 2, 0.368, 0.11, 0.56, 0.72, -0.25),
    (4, 1, 0.228, 0.058, 0.18, 1.14, -0.62),
    (3, 1, 0.248, 0.065, 0.21, 1.06, -0.56),
    (2, 1, 0.278, 0.075, 0.28, 0.95, -0.48),
];

/// Precomputed gap parameter tables for nucleotide scoring.
/// Format: (gap_open, gap_extend, lambda, K, H, alpha, beta)
pub const BLASTN_PARAMS_1_2: &[(i32, i32, f64, f64, f64, f64, f64)] = &[
    // reward=1, penalty=-2
    (0, 0, 0.549, 0.205, 1.02, 0.45, -0.19),  // ungapped/linear
    (2, 2, 0.439, 0.14, 0.64, 0.60, -0.26),
    (1, 2, 0.477, 0.16, 0.78, 0.53, -0.22),
    (0, 2, 0.549, 0.205, 1.02, 0.45, -0.19),
    (3, 1, 0.429, 0.13, 0.60, 0.62, -0.27),
    (2, 1, 0.456, 0.15, 0.71, 0.56, -0.24),
    (1, 1, 0.499, 0.18, 0.88, 0.49, -0.21),
];

/// Look up gapped Karlin-Altschul parameters for nucleotide scoring.
pub fn lookup_nucleotide_params(
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
) -> Option<GappedParams> {
    let table = match (reward, penalty) {
        (1, -3) => BLASTN_PARAMS_1_3,
        (1, -2) => BLASTN_PARAMS_1_2,
        _ => return None,
    };

    for &(go, ge, lambda, k, h, alpha, beta) in table {
        if go == gap_open && ge == gap_extend {
            return Some(GappedParams {
                gap_open: go,
                gap_extend: ge,
                lambda,
                k,
                h,
                alpha,
                beta,
            });
        }
    }
    None
}

/// Compute effective search space.
/// eff_length = db_length - num_seqs * length_adjustment
/// eff_query_length = query_length - length_adjustment
/// search_space = eff_length * eff_query_length
pub fn compute_search_space(
    query_length: i64,
    db_length: i64,
    num_seqs: i32,
    length_adjustment: i32,
) -> f64 {
    let eff_query = (query_length - length_adjustment as i64).max(1);
    let eff_db = (db_length - num_seqs as i64 * length_adjustment as i64).max(1);
    eff_query as f64 * eff_db as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_karlin_blk_valid() {
        let kbp = KarlinBlk {
            lambda: 0.208,
            k: 0.049,
            log_k: 0.049_f64.ln(),
            h: 0.14,
        };
        assert!(kbp.is_valid());
    }

    #[test]
    fn test_raw_to_bit() {
        let kbp = KarlinBlk {
            lambda: 0.208,
            k: 0.049,
            log_k: 0.049_f64.ln(),
            h: 0.14,
        };
        let bit = kbp.raw_to_bit(50);
        assert!(bit > 0.0);
        // lambda * 50 / ln(2) - ln(K)/ln(2) ≈ 0.208*50/0.693 - ln(0.049)/0.693
        // ≈ 15.0 + 4.35 ≈ 19.35
        assert!((bit - 19.35).abs() < 0.5);
    }

    #[test]
    fn test_lookup_nucleotide_params() {
        let params = lookup_nucleotide_params(1, -3, 5, 2).unwrap();
        assert_eq!(params.gap_open, 5);
        assert_eq!(params.gap_extend, 2);
        assert!((params.lambda - 0.208).abs() < 0.001);
    }

    #[test]
    fn test_search_space() {
        let ss = compute_search_space(100, 1000000, 2000, 20);
        assert!(ss > 0.0);
        // (100-20) * (1000000 - 2000*20) = 80 * 960000 = 76800000
        assert_eq!(ss, 76800000.0);
    }
}
