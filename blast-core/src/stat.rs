//! Rust equivalent of blast_stat.c — Karlin-Altschul statistics.
//! This is the mathematical core for computing E-values and bit scores.


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
    (5, 2, 0.208, 0.049, 0.14, 1.24, -0.70),
    (2, 2, 0.278, 0.075, 0.27, 0.96, -0.52),
    (1, 2, 0.308, 0.084, 0.34, 0.86, -0.42),
    (0, 2, 0.368, 0.11, 0.56, 0.72, -0.25),
    (4, 1, 0.228, 0.058, 0.18, 1.14, -0.62),
    (3, 1, 0.248, 0.065, 0.21, 1.06, -0.56),
    (2, 1, 0.278, 0.075, 0.28, 0.95, -0.48),
];

pub const BLASTN_PARAMS_1_2: &[(i32, i32, f64, f64, f64, f64, f64)] = &[
    (0, 0, 0.549, 0.205, 1.02, 0.45, -0.19),
    (2, 2, 0.439, 0.14, 0.64, 0.60, -0.26),
    (1, 2, 0.477, 0.16, 0.78, 0.53, -0.22),
    (0, 2, 0.549, 0.205, 1.02, 0.45, -0.19),
    (3, 1, 0.429, 0.13, 0.60, 0.62, -0.27),
    (2, 1, 0.456, 0.15, 0.71, 0.56, -0.24),
    (1, 1, 0.499, 0.18, 0.88, 0.49, -0.21),
];

pub const BLASTN_PARAMS_2_3: &[(i32, i32, f64, f64, f64, f64, f64)] = &[
    (4, 4, 0.22, 0.065, 0.35, 1.00, -0.43),
    (2, 4, 0.29, 0.086, 0.52, 0.82, -0.31),
    (0, 4, 0.38, 0.13, 0.92, 0.62, -0.17),
    (3, 3, 0.23, 0.068, 0.38, 0.96, -0.40),
    (6, 2, 0.21, 0.057, 0.30, 1.06, -0.48),
    (5, 2, 0.23, 0.068, 0.38, 0.95, -0.39),
    (4, 2, 0.25, 0.078, 0.45, 0.87, -0.34),
];

pub const BLASTN_PARAMS_1_1: &[(i32, i32, f64, f64, f64, f64, f64)] = &[
    (3, 2, 0.549, 0.205, 1.02, 0.45, -0.19),
    (2, 2, 0.659, 0.26, 1.50, 0.37, -0.14),
    (1, 2, 0.730, 0.29, 1.83, 0.33, -0.12),
    (0, 2, 0.854, 0.37, 2.69, 0.27, -0.08),
    (4, 1, 0.549, 0.205, 1.02, 0.45, -0.19),
    (3, 1, 0.624, 0.24, 1.30, 0.39, -0.16),
    (2, 1, 0.693, 0.27, 1.64, 0.35, -0.13),
];

/// Compute ungapped Karlin-Altschul lambda parameter for nucleotide scoring.
/// Solves: 0.25 * exp(lambda * reward) + 0.75 * exp(lambda * penalty) = 1
/// using Newton's method.
pub fn compute_ungapped_lambda(reward: i32, penalty: i32) -> f64 {
    let r = reward as f64;
    let p = penalty as f64;
    // Initial guess
    let mut lambda = 1.0;
    for _ in 0..100 {
        let er = (lambda * r).exp();
        let ep = (lambda * p).exp();
        let f = 0.25 * er + 0.75 * ep - 1.0;
        let fp = 0.25 * r * er + 0.75 * p * ep;
        if fp.abs() < 1e-30 { break; }
        let delta = f / fp;
        lambda -= delta;
        if delta.abs() < 1e-12 { break; }
    }
    lambda
}

/// Compute ungapped Karlin-Altschul K parameter for nucleotide scoring.
/// K = (1 - exp(-lambda)) * (1 - exp(-lambda)) / (sum of prob * exp(lambda*score) * score^2 / 2)
/// Simplified for uniform base frequencies.
pub fn compute_ungapped_k(lambda: f64, reward: i32, penalty: i32) -> f64 {
    let r = reward as f64;
    let p = penalty as f64;
    let er = (lambda * r).exp();
    let ep = (lambda * p).exp();
    // H = lambda * (0.25 * r * er + 0.75 * p * ep)
    let h = lambda * (0.25 * r * er + 0.75 * p * ep);
    // K approximation for simple scoring
    // K = H / (lambda * lambda * variance)
    // For nucleotide with uniform freqs, use the standard approximation
    let variance = 0.25 * r * r * er + 0.75 * p * p * ep;
    if variance.abs() < 1e-30 { return 0.1; }
    h / (lambda * variance)
}

/// Compute ungapped KarlinBlk for nucleotide scoring with given reward/penalty.
pub fn compute_ungapped_kbp(reward: i32, penalty: i32) -> KarlinBlk {
    let lambda = compute_ungapped_lambda(reward, penalty);
    let k = compute_ungapped_k(lambda, reward, penalty);
    let h = lambda * (0.25 * (reward as f64) * (lambda * reward as f64).exp()
        + 0.75 * (penalty as f64) * (lambda * penalty as f64).exp());
    KarlinBlk {
        lambda,
        k,
        log_k: k.ln(),
        h,
    }
}

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
        (2, -3) => BLASTN_PARAMS_2_3,
        (1, -1) => BLASTN_PARAMS_1_1,
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
    fn test_compute_ungapped_lambda() {
        // For reward=1, penalty=-3: known lambda ≈ 1.374
        let lambda = compute_ungapped_lambda(1, -3);
        assert!((lambda - 1.374).abs() < 0.01,
            "lambda for 1/-3 should be ~1.374, got {}", lambda);

        // For reward=2, penalty=-3: different lambda
        let lambda2 = compute_ungapped_lambda(2, -3);
        assert!(lambda2 > 0.0 && lambda2 < lambda,
            "lambda for 2/-3 should be positive and < 1/-3 lambda");
    }

    #[test]
    fn test_compute_ungapped_kbp() {
        let kbp = compute_ungapped_kbp(1, -3);
        assert!(kbp.lambda > 0.0);
        assert!(kbp.k > 0.0);
        assert!(kbp.h > 0.0);
        assert!(kbp.is_valid());
    }

    #[test]
    fn test_search_space() {
        let ss = compute_search_space(100, 1000000, 2000, 20);
        assert!(ss > 0.0);
        // (100-20) * (1000000 - 2000*20) = 80 * 960000 = 76800000
        assert_eq!(ss, 76800000.0);
    }
}
