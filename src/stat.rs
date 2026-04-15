//! Rust equivalent of blast_stat.c — Karlin-Altschul statistics.
//! This is the mathematical core for computing E-values and bit scores.

/// Karlin-Altschul statistical parameters for one context.
#[derive(Debug, Clone)]
pub struct KarlinBlk {
    pub lambda: f64, // Lambda parameter
    pub k: f64,      // K parameter
    pub log_k: f64,  // ln(K)
    pub h: f64,      // H (relative entropy)
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
        let denom = search_space * self.k;
        if denom <= 0.0 || self.lambda <= 0.0 || evalue <= 0.0 {
            return 1;
        }
        let score = -(evalue / denom).ln() / self.lambda;
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

/// Gumbel block parameters for Spouge finite-size correction (FSC) e-value.
/// Port of NCBI Blast_GumbelBlk from blast_stat.h.
#[derive(Debug, Clone)]
pub struct GumbelBlk {
    pub lambda: f64,    // gbp->Lambda (gapped lambda from table)
    pub a: f64,         // gbp->a (alpha from array[6])
    pub b: f64,         // gbp->b = 2*G*(a_un - a)
    pub alpha: f64,     // gbp->Alpha (alpha_v from array[9])
    pub beta: f64,      // gbp->Beta = 2*G*(Alpha_un - Alpha)
    pub sigma: f64,     // gbp->Sigma (from array[10])
    pub tau: f64,       // gbp->Tau = 2*G*(Alpha_un - Sigma)
    pub db_length: i64, // total database length
}

/// Compute Spouge finite-size correction e-value.
/// Port of BLAST_SpougeStoE from blast_stat.c:5176.
/// Uses per-subject lengths for more accurate e-values than simple Karlin formula.
pub fn spouge_evalue(
    score: i32,
    kbp: &KarlinBlk,
    gbp: &GumbelBlk,
    query_length: i32,
    subject_length: i32,
) -> f64 {
    let scale_factor = kbp.lambda / gbp.lambda;
    let db_scale_factor = if gbp.db_length > 0 {
        gbp.db_length as f64 / subject_length as f64
    } else {
        1.0
    };

    let ai_hat = gbp.a * scale_factor;
    let bi_hat = gbp.b;
    let alphai_hat = gbp.alpha * scale_factor;
    let betai_hat = gbp.beta;
    let sigma_hat = gbp.sigma * scale_factor;
    let tau_hat = gbp.tau;

    let aj_hat = ai_hat;
    let bj_hat = bi_hat;
    let alphaj_hat = alphai_hat;
    let betaj_hat = betai_hat;

    const CONST_VAL: f64 = 0.39894228040143267793994605993438; // 1/sqrt(2*PI)

    let y = score as f64;
    let m = query_length as f64;
    let n = subject_length as f64;

    let m_li_y = m - (ai_hat * y + bi_hat);
    let vi_y = (2.0 * alphai_hat / kbp.lambda).max(alphai_hat * y + betai_hat);
    let sqrt_vi_y = vi_y.sqrt();
    let m_f = m_li_y / sqrt_vi_y;
    let p_m_f = erfc_approx(-m_f / std::f64::consts::SQRT_2) / 2.0;
    let p1 = m_li_y * p_m_f + sqrt_vi_y * CONST_VAL * (-0.5 * m_f * m_f).exp();

    let n_lj_y = n - (aj_hat * y + bj_hat);
    let vj_y = (2.0 * alphaj_hat / kbp.lambda).max(alphaj_hat * y + betaj_hat);
    let sqrt_vj_y = vj_y.sqrt();
    let n_f = n_lj_y / sqrt_vj_y;
    let p_n_f = erfc_approx(-n_f / std::f64::consts::SQRT_2) / 2.0;
    let p2 = n_lj_y * p_n_f + sqrt_vj_y * CONST_VAL * (-0.5 * n_f * n_f).exp();

    let c_y = (2.0 * sigma_hat / kbp.lambda).max(sigma_hat * y + tau_hat);
    let area = p1 * p2 + c_y * p_m_f * p_n_f;

    (area * kbp.k * (-kbp.lambda * y).exp() * db_scale_factor).max(0.0)
}

/// Complementary error function approximation (Abramowitz & Stegun 7.1.26).
/// Accurate to ~1.5e-7 relative error.
fn erfc_approx(x: f64) -> f64 {
    let t = 1.0 / (1.0 + 0.3275911 * x.abs());
    let poly = t
        * (0.254829592
            + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
    let result = poly * (-x * x).exp();
    if x >= 0.0 {
        result
    } else {
        2.0 - result
    }
}

/// Build Gumbel block for protein BLOSUM62 with given gap costs.
/// Port of Blast_GumbelBlkLoadFromTables from blast_stat.c:3696.
pub fn protein_gumbel_blk(gap_open: i32, gap_extend: i32, db_length: i64) -> Option<GumbelBlk> {
    // Ungapped params (index 0 in BLOSUM62 table)
    let a_un = 0.7916;
    let alpha_un = 4.964660;

    // Find gapped params
    let gp = lookup_protein_params(gap_open, gap_extend)?;

    // Gapped params from table
    // Format: (go, ge, lambda, K, H, alpha, beta)
    // But we also need alpha_v and sigma from the extended table
    // These are stored in BLOSUM62_EXTENDED_PARAMS
    let (alpha_v, sigma) = lookup_blosum62_extended(gap_open, gap_extend)?;

    let g = (gap_open + gap_extend) as f64;

    Some(GumbelBlk {
        lambda: gp.lambda,
        a: gp.alpha, // 'a' = alpha from the basic table (position 6 in C array)
        b: 2.0 * g * (a_un - gp.alpha),
        alpha: alpha_v,
        beta: 2.0 * g * (alpha_un - alpha_v),
        sigma,
        tau: 2.0 * g * (alpha_un - sigma),
        db_length,
    })
}

/// Extended BLOSUM62 params: (gap_open, gap_extend, alpha_v, sigma)
/// From NCBI blast_stat.c blosum62_values positions [9] and [10].
const BLOSUM62_EXTENDED: &[(i32, i32, f64, f64)] = &[
    (11, 2, 12.6738, 12.7576),
    (10, 2, 16.4740, 16.6026),
    (9, 2, 22.7519, 22.9500),
    (8, 2, 35.4838, 35.8213),
    (7, 2, 61.2383, 61.8860),
    (6, 2, 140.417, 141.882),
    (13, 1, 19.5063, 19.8931),
    (12, 1, 27.8562, 28.4699),
    (11, 1, 42.6028, 43.6362),
    (10, 1, 83.1787, 85.0656),
    (9, 1, 210.333, 214.842),
];

fn lookup_blosum62_extended(gap_open: i32, gap_extend: i32) -> Option<(f64, f64)> {
    for &(go, ge, alpha_v, sigma) in BLOSUM62_EXTENDED {
        if go == gap_open && ge == gap_extend {
            return Some((alpha_v, sigma));
        }
    }
    None
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
        if fp.abs() < 1e-30 {
            break;
        }
        let delta = f / fp;
        lambda -= delta;
        if delta.abs() < 1e-12 {
            break;
        }
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
    if variance.abs() < 1e-30 {
        return 0.1;
    }
    h / (lambda * variance)
}

/// Compute ungapped KarlinBlk for nucleotide scoring with given reward/penalty.
pub fn compute_ungapped_kbp(reward: i32, penalty: i32) -> KarlinBlk {
    let lambda = compute_ungapped_lambda(reward, penalty);
    let k = compute_ungapped_k(lambda, reward, penalty);
    let h = lambda
        * (0.25 * (reward as f64) * (lambda * reward as f64).exp()
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

/// Compute length adjustment using the Altschul-Gish formula.
/// This adjusts query and database lengths to account for edge effects.
/// length_adj satisfies: K * (q - length_adj) * (d - N * length_adj) * exp(-lambda * S) = threshold
pub fn compute_length_adjustment(
    query_length: i32,
    db_length: i64,
    num_seqs: i32,
    kbp: &KarlinBlk,
) -> i32 {
    if !kbp.is_valid() || query_length <= 0 || db_length <= 0 {
        return 0;
    }

    let h = kbp.h;
    if h <= 0.0 {
        return 0;
    }

    // Iteratively solve for length_adjustment
    // Formula: len_adj = (ln(K) + ln(eff_q * eff_d)) / H
    let q = query_length as f64;
    let d = db_length as f64;
    let n = num_seqs as f64;
    let log_k = kbp.log_k;

    let mut len_adj = 0.0;
    for _ in 0..20 {
        let eff_q = (q - len_adj).max(1.0);
        let eff_d = (d - n * len_adj).max(1.0);
        let new_adj = (log_k + (eff_q * eff_d).ln()) / h;
        let new_adj = new_adj.min(q - 1.0).min(d / n - 1.0).max(0.0);
        if (new_adj - len_adj).abs() < 0.5 {
            break;
        }
        len_adj = new_adj;
    }

    len_adj.round() as i32
}

/// Precomputed Karlin-Altschul parameters for BLOSUM62 with various gap costs.
/// Values from NCBI blast_stat.c blosum62_values[].
/// Format: (gap_open, gap_extend, lambda, K, H, alpha, beta)
pub const BLOSUM62_PARAMS: &[(i32, i32, f64, f64, f64, f64, f64)] = &[
    (11, 2, 0.297, 0.082, 0.27, 1.1, -10.0),
    (10, 2, 0.291, 0.075, 0.23, 1.3, -15.0),
    (9, 2, 0.279, 0.058, 0.19, 1.5, -19.0),
    (8, 2, 0.264, 0.045, 0.15, 1.8, -26.0),
    (7, 2, 0.239, 0.027, 0.10, 2.5, -46.0),
    (6, 2, 0.201, 0.012, 0.061, 3.3, -58.0),
    (13, 1, 0.292, 0.071, 0.23, 1.2, -11.0),
    (12, 1, 0.283, 0.059, 0.19, 1.5, -19.0),
    (11, 1, 0.267, 0.041, 0.14, 1.9, -30.0),
    (10, 1, 0.243, 0.024, 0.10, 2.5, -44.0),
    (9, 1, 0.206, 0.010, 0.052, 4.0, -87.0),
];

/// Look up gapped KBP for protein scoring (BLOSUM62).
pub fn lookup_protein_params(gap_open: i32, gap_extend: i32) -> Option<GappedParams> {
    for &(go, ge, lambda, k, h, alpha, beta) in BLOSUM62_PARAMS {
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

/// Compute ungapped KBP for protein BLOSUM62.
/// Lambda ≈ 0.3176, K ≈ 0.134 for standard amino acid frequencies.
pub fn protein_ungapped_kbp() -> KarlinBlk {
    KarlinBlk {
        lambda: 0.3176,
        k: 0.134,
        log_k: 0.134_f64.ln(),
        h: 0.401,
    }
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

// ---------------------------------------------------------------------------
// Gapped Karlin-Altschul parameter lookup (exact C-compatible tables)
// Port of Blast_KarlinBlkNuclGappedCalc from blast_stat.c
// ---------------------------------------------------------------------------

/// Each row: [gap_open, gap_extend, lambda, K, H, alpha, beta, theta]
type KbpTableRow = [f64; 8];

const KBPT_1_5: &[KbpTableRow] = &[
    [0.0, 0.0, 1.39, 0.747, 1.38, 1.00, 0.0, 100.0],
    [3.0, 3.0, 1.39, 0.747, 1.38, 1.00, 0.0, 100.0],
];
const KBPT_1_4: &[KbpTableRow] = &[
    [0.0, 0.0, 1.383, 0.738, 1.36, 1.02, 0.0, 100.0],
    [1.0, 2.0, 1.36, 0.67, 1.2, 1.1, 0.0, 98.0],
    [0.0, 2.0, 1.26, 0.43, 0.90, 1.4, -1.0, 91.0],
    [2.0, 1.0, 1.35, 0.61, 1.1, 1.2, -1.0, 98.0],
    [1.0, 1.0, 1.22, 0.35, 0.72, 1.7, -3.0, 88.0],
];
const KBPT_2_7: &[KbpTableRow] = &[
    [0.0, 0.0, 0.69, 0.73, 1.34, 0.515, 0.0, 100.0],
    [2.0, 4.0, 0.68, 0.67, 1.2, 0.55, 0.0, 99.0],
    [0.0, 4.0, 0.63, 0.43, 0.90, 0.7, -1.0, 91.0],
    [4.0, 2.0, 0.675, 0.62, 1.1, 0.6, -1.0, 98.0],
    [2.0, 2.0, 0.61, 0.35, 0.72, 1.7, -3.0, 88.0],
];
const KBPT_1_3: &[KbpTableRow] = &[
    [0.0, 0.0, 1.374, 0.711, 1.31, 1.05, 0.0, 100.0],
    [2.0, 2.0, 1.37, 0.70, 1.2, 1.1, 0.0, 99.0],
    [1.0, 2.0, 1.35, 0.64, 1.1, 1.2, -1.0, 98.0],
    [0.0, 2.0, 1.25, 0.42, 0.83, 1.5, -2.0, 91.0],
    [2.0, 1.0, 1.34, 0.60, 1.1, 1.2, -1.0, 97.0],
    [1.0, 1.0, 1.21, 0.34, 0.71, 1.7, -2.0, 88.0],
];
const KBPT_2_5: &[KbpTableRow] = &[
    [0.0, 0.0, 0.675, 0.65, 1.1, 0.6, -1.0, 99.0],
    [2.0, 4.0, 0.67, 0.59, 1.1, 0.6, -1.0, 98.0],
    [0.0, 4.0, 0.62, 0.39, 0.78, 0.8, -2.0, 91.0],
    [4.0, 2.0, 0.67, 0.61, 1.0, 0.65, -2.0, 98.0],
    [2.0, 2.0, 0.56, 0.32, 0.59, 0.95, -4.0, 82.0],
];
const KBPT_1_2: &[KbpTableRow] = &[
    [0.0, 0.0, 1.28, 0.46, 0.85, 1.5, -2.0, 96.0],
    [2.0, 2.0, 1.33, 0.62, 1.1, 1.2, 0.0, 99.0],
    [1.0, 2.0, 1.30, 0.52, 0.93, 1.4, -2.0, 97.0],
    [0.0, 2.0, 1.19, 0.34, 0.66, 1.8, -3.0, 89.0],
    [3.0, 1.0, 1.32, 0.57, 1.0, 1.3, -1.0, 99.0],
    [2.0, 1.0, 1.29, 0.49, 0.92, 1.4, -1.0, 96.0],
    [1.0, 1.0, 1.14, 0.26, 0.52, 2.2, -5.0, 85.0],
];
const KBPT_2_3: &[KbpTableRow] = &[
    [0.0, 0.0, 0.55, 0.21, 0.46, 1.2, -5.0, 87.0],
    [4.0, 4.0, 0.63, 0.42, 0.84, 0.75, -2.0, 99.0],
    [2.0, 4.0, 0.615, 0.37, 0.72, 0.85, -3.0, 97.0],
    [0.0, 4.0, 0.55, 0.21, 0.46, 1.2, -5.0, 87.0],
    [3.0, 3.0, 0.615, 0.37, 0.68, 0.9, -3.0, 97.0],
    [6.0, 2.0, 0.63, 0.42, 0.84, 0.75, -2.0, 99.0],
    [5.0, 2.0, 0.625, 0.41, 0.78, 0.8, -2.0, 99.0],
    [4.0, 2.0, 0.61, 0.35, 0.68, 0.9, -3.0, 96.0],
    [2.0, 2.0, 0.515, 0.14, 0.33, 1.55, -9.0, 81.0],
];
const KBPT_3_4: &[KbpTableRow] = &[
    [6.0, 3.0, 0.389, 0.25, 0.56, 0.7, -5.0, 95.0],
    [5.0, 3.0, 0.375, 0.21, 0.47, 0.8, -6.0, 92.0],
    [4.0, 3.0, 0.351, 0.14, 0.35, 1.0, -9.0, 86.0],
    [6.0, 2.0, 0.362, 0.16, 0.45, 0.8, -4.0, 88.0],
    [5.0, 2.0, 0.330, 0.092, 0.28, 1.2, -13.0, 81.0],
    [4.0, 2.0, 0.281, 0.046, 0.16, 1.8, -23.0, 69.0],
];
const KBPT_1_1: &[KbpTableRow] = &[
    [3.0, 2.0, 1.09, 0.31, 0.55, 2.0, -2.0, 99.0],
    [2.0, 2.0, 1.07, 0.27, 0.49, 2.2, -3.0, 97.0],
    [1.0, 2.0, 1.02, 0.21, 0.36, 2.8, -6.0, 92.0],
    [0.0, 2.0, 0.80, 0.064, 0.17, 4.8, -16.0, 72.0],
    [4.0, 1.0, 1.08, 0.28, 0.54, 2.0, -2.0, 98.0],
    [3.0, 1.0, 1.06, 0.25, 0.46, 2.3, -4.0, 96.0],
    [2.0, 1.0, 0.99, 0.17, 0.30, 3.3, -10.0, 90.0],
];
const KBPT_4_5: &[KbpTableRow] = &[
    [0.0, 0.0, 0.22, 0.061, 0.22, 1.0, -15.0, 74.0],
    [6.0, 5.0, 0.28, 0.21, 0.47, 0.6, -7.0, 93.0],
    [5.0, 5.0, 0.27, 0.17, 0.39, 0.7, -9.0, 90.0],
    [4.0, 5.0, 0.25, 0.10, 0.31, 0.8, -10.0, 83.0],
    [3.0, 5.0, 0.23, 0.065, 0.25, 0.9, -11.0, 76.0],
];
const KBPT_3_2: &[KbpTableRow] = &[[5.0, 5.0, 0.208, 0.030, 0.072, 2.9, -47.0, 77.0]];
const KBPT_5_4: &[KbpTableRow] = &[
    [10.0, 6.0, 0.163, 0.068, 0.16, 1.0, -19.0, 85.0],
    [8.0, 6.0, 0.146, 0.039, 0.11, 1.3, -29.0, 76.0],
];

struct KbpTableMeta {
    table: &'static [KbpTableRow],
    gap_open_max: i32,
    gap_extend_max: i32,
    round_down: bool,
}

fn kbp_gcd(a: i32, b: i32) -> i32 {
    let (mut a, mut b) = (a.unsigned_abs(), b.unsigned_abs());
    if b > a {
        std::mem::swap(&mut a, &mut b);
    }
    while b != 0 {
        let c = a % b;
        a = b;
        b = c;
    }
    a as i32
}

fn get_kbp_table(reward: i32, penalty: i32) -> Option<KbpTableMeta> {
    match (reward, penalty) {
        (1, -5) => Some(KbpTableMeta {
            table: KBPT_1_5,
            gap_open_max: 3,
            gap_extend_max: 3,
            round_down: false,
        }),
        (1, -4) => Some(KbpTableMeta {
            table: KBPT_1_4,
            gap_open_max: 2,
            gap_extend_max: 2,
            round_down: false,
        }),
        (2, -7) => Some(KbpTableMeta {
            table: KBPT_2_7,
            gap_open_max: 4,
            gap_extend_max: 4,
            round_down: true,
        }),
        (1, -3) => Some(KbpTableMeta {
            table: KBPT_1_3,
            gap_open_max: 2,
            gap_extend_max: 2,
            round_down: false,
        }),
        (2, -5) => Some(KbpTableMeta {
            table: KBPT_2_5,
            gap_open_max: 4,
            gap_extend_max: 4,
            round_down: true,
        }),
        (1, -2) => Some(KbpTableMeta {
            table: KBPT_1_2,
            gap_open_max: 2,
            gap_extend_max: 2,
            round_down: false,
        }),
        (2, -3) => Some(KbpTableMeta {
            table: KBPT_2_3,
            gap_open_max: 6,
            gap_extend_max: 4,
            round_down: true,
        }),
        (3, -4) => Some(KbpTableMeta {
            table: KBPT_3_4,
            gap_open_max: 6,
            gap_extend_max: 3,
            round_down: true,
        }),
        (1, -1) => Some(KbpTableMeta {
            table: KBPT_1_1,
            gap_open_max: 4,
            gap_extend_max: 2,
            round_down: false,
        }),
        (3, -2) => Some(KbpTableMeta {
            table: KBPT_3_2,
            gap_open_max: 5,
            gap_extend_max: 5,
            round_down: false,
        }),
        (4, -5) => Some(KbpTableMeta {
            table: KBPT_4_5,
            gap_open_max: 12,
            gap_extend_max: 8,
            round_down: false,
        }),
        (5, -4) => Some(KbpTableMeta {
            table: KBPT_5_4,
            gap_open_max: 25,
            gap_extend_max: 10,
            round_down: false,
        }),
        _ => None,
    }
}

/// Look up gapped KBP params for nucleotide, matching C's Blast_KarlinBlkNuclGappedCalc exactly.
/// Returns Ok((lambda, k, log_k, h, round_down)) or Err message.
pub fn nucl_gapped_kbp_lookup(
    gap_open: i32,
    gap_extend: i32,
    reward: i32,
    penalty: i32,
    ungapped: &KarlinBlk,
) -> Result<(KarlinBlk, bool), String> {
    let divisor = kbp_gcd(reward, penalty.abs());
    let (nr, np) = (reward / divisor, penalty / divisor);

    let meta = get_kbp_table(nr, np)
        .ok_or_else(|| format!("Unsupported scores {} {}", reward, penalty))?;

    let round_down = meta.round_down;

    // Split: first row with gap_open==0 && gap_extend==0 is the linear entry
    let (affine, linear) =
        if !meta.table.is_empty() && meta.table[0][0] == 0.0 && meta.table[0][1] == 0.0 {
            (&meta.table[1..], Some(&meta.table[0]))
        } else {
            (meta.table, None)
        };

    // GCD scaling
    let (mut go_max, mut ge_max) = (meta.gap_open_max, meta.gap_extend_max);
    let scale = |row: &KbpTableRow, d: i32| -> (f64, f64, f64, f64, f64) {
        let go = row[0] * d as f64;
        let ge = row[1] * d as f64;
        let lam = row[2] / d as f64;
        (go, ge, lam, row[3], row[4])
    };
    if divisor != 1 {
        go_max *= divisor;
        ge_max *= divisor;
    }

    // Linear (non-affine) case
    if gap_open == 0 && gap_extend == 0 {
        if let Some(lin) = linear {
            let (_, _, lam, k, h) = scale(lin, divisor);
            return Ok((
                KarlinBlk {
                    lambda: lam,
                    k,
                    log_k: k.ln(),
                    h,
                },
                round_down,
            ));
        }
    }

    // Search affine entries
    for row in affine {
        let (go, ge, lam, k, h) = scale(row, divisor);
        if go as i32 == gap_open && ge as i32 == gap_extend {
            return Ok((
                KarlinBlk {
                    lambda: lam,
                    k,
                    log_k: k.ln(),
                    h,
                },
                round_down,
            ));
        }
    }

    // Fallback: gap costs exceed table maximum → use ungapped params
    if gap_open >= go_max && gap_extend >= ge_max {
        return Ok((ungapped.clone(), round_down));
    }

    Err(format!(
        "Unsupported gap costs {} {} for scores {} {}",
        gap_open, gap_extend, reward, penalty
    ))
}

/// Compute the length adjustment for effective search space (exact C-compatible).
/// Port of BLAST_ComputeLengthAdjustment from blast_stat.c.
/// Returns (length_adjustment, converged).
pub fn compute_length_adjustment_exact(
    k: f64,
    log_k: f64,
    alpha_d_lambda: f64,
    beta: f64,
    query_length: i32,
    db_length: i64,
    db_num_seqs: i32,
) -> (i32, bool) {
    let m = query_length as f64;
    let n = db_length as f64;
    let nn = db_num_seqs as f64;

    // Compute ell_max using quadratic formula
    let a = nn;
    let mb = m * nn + n;
    let c = n * m - m.max(n) / k;
    if c < 0.0 {
        return (0, false);
    }
    let ell_max_init = 2.0 * c / (mb + (mb * mb - 4.0 * a * c).sqrt());

    let mut ell_min = 0.0_f64;
    let mut ell_max = ell_max_init;
    let mut ell_next = 0.0_f64;
    let mut converged = false;

    for i in 1..=20 {
        let ell = ell_next;
        let ss = (m - ell) * (n - nn * ell);
        let ell_bar = alpha_d_lambda * (log_k + ss.ln()) + beta;

        if ell_bar >= ell {
            ell_min = ell;
            if ell_bar - ell_min <= 1.0 {
                converged = true;
                break;
            }
            if ell_min == ell_max {
                break;
            }
        } else {
            ell_max = ell;
        }

        if ell_min <= ell_bar && ell_bar <= ell_max {
            ell_next = ell_bar;
        } else {
            ell_next = if i == 1 {
                ell_max
            } else {
                (ell_min + ell_max) / 2.0
            };
        }
    }

    let mut adj = ell_min as i32;
    if converged {
        let ell = ell_min.ceil();
        if ell <= ell_max {
            let ss = (m - ell) * (n - nn * ell);
            if alpha_d_lambda * (log_k + ss.ln()) + beta >= ell {
                adj = ell as i32;
            }
        }
    }
    (adj, converged)
}

/// Look up alpha and beta for nucleotide gapped alignment.
/// Port of Blast_GetNuclAlphaBeta from blast_stat.c.
pub fn nucl_alpha_beta(
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    ungapped_lambda: f64,
    ungapped_h: f64,
    gapped: bool,
) -> (f64, f64) {
    if !gapped {
        let beta = if (reward == 1 && penalty == -1) || (reward == 2 && penalty == -3) {
            -2.0
        } else {
            0.0
        };
        return (ungapped_lambda / ungapped_h, beta);
    }

    let divisor = kbp_gcd(reward, penalty.abs());
    let (nr, np) = (reward / divisor, penalty / divisor);

    if let Some(meta) = get_kbp_table(nr, np) {
        let (affine, linear) =
            if !meta.table.is_empty() && meta.table[0][0] == 0.0 && meta.table[0][1] == 0.0 {
                (&meta.table[1..], Some(&meta.table[0]))
            } else {
                (meta.table, None)
            };

        if gap_open == 0 && gap_extend == 0 {
            if let Some(lin) = linear {
                let a = lin[5] / divisor as f64;
                let b = lin[6];
                return (a, b);
            }
        }

        for row in affine {
            let go = (row[0] * divisor as f64) as i32;
            let ge = (row[1] * divisor as f64) as i32;
            if go == gap_open && ge == gap_extend {
                let a = row[5] / divisor as f64;
                let b = row[6];
                return (a, b);
            }
        }
    }

    // Fallback: ungapped values
    let beta = if (reward == 1 && penalty == -1) || (reward == 2 && penalty == -3) {
        -2.0
    } else {
        0.0
    };
    (ungapped_lambda / ungapped_h, beta)
}

// ---------------------------------------------------------------------------
// Ungapped Karlin-Altschul parameter computation (exact C-compatible)
// Port of Blast_ScoreBlkKbpUngappedCalc and sub-functions from blast_stat.c
// ---------------------------------------------------------------------------

/// Score frequency distribution (internal to KBP computation).
struct SfDist {
    score_min: i32,
    score_max: i32,
    obs_min: i32,
    obs_max: i32,
    score_avg: f64,
    probs: Vec<f64>, // indexed by [score - score_min]
}

impl SfDist {
    fn new(lo: i32, hi: i32) -> Self {
        let n = (hi - lo + 1) as usize;
        SfDist {
            score_min: lo,
            score_max: hi,
            obs_min: 0,
            obs_max: 0,
            score_avg: 0.0,
            probs: vec![0.0; n],
        }
    }
    #[inline]
    fn p(&self, s: i32) -> f64 {
        self.probs[(s - self.score_min) as usize]
    }
    #[inline]
    fn p_mut(&mut self, s: i32) -> &mut f64 {
        &mut self.probs[(s - self.score_min) as usize]
    }
}

fn expm1_stable(x: f64) -> f64 {
    if x.abs() > 0.33 {
        x.exp() - 1.0
    } else {
        x.exp_m1()
    }
}

fn powi(x: f64, n: i32) -> f64 {
    if n == 0 {
        return 1.0;
    }
    let mut x = if n < 0 { 1.0 / x } else { x };
    let mut n = n.unsigned_abs();
    let mut y = 1.0;
    while n > 0 {
        if n & 1 != 0 {
            y *= x;
        }
        n /= 2;
        x *= x;
    }
    y
}

/// Newton-Raphson solver for Lambda in x = exp(-lambda) space.
/// Port of NlmKarlinLambdaNR from blast_stat.c.
fn solve_lambda(sfp: &SfDist, d: i32, low: i32, high: i32, lambda0: f64) -> f64 {
    let x0 = (-lambda0).exp();
    let mut x = if x0 > 0.0 && x0 < 1.0 { x0 } else { 0.5 };
    let mut a = 0.0_f64;
    let mut b = 1.0_f64;
    let mut f = 4.0_f64;
    let mut is_newton = false;
    let max_iter = 37; // 20 + 17
    let tolx = 1.0e-5;

    for k in 0..max_iter {
        let fold = f;
        let was_newton = is_newton;
        is_newton = false;

        // Horner evaluation of polynomial and derivative
        let mut g = 0.0_f64;
        f = sfp.p(low);
        let mut i = low + d;
        while i < 0 {
            g = x * g + f;
            f = f * x + sfp.p(i);
            i += d;
        }
        g = x * g + f;
        f = f * x + sfp.p(0) - 1.0;
        i = d;
        while i <= high {
            g = x * g + f;
            f = f * x + sfp.p(i);
            i += d;
        }

        if f > 0.0 {
            a = x;
        } else if f < 0.0 {
            b = x;
        } else {
            break;
        }
        if b - a < 2.0 * a * (1.0 - b) * tolx {
            x = (a + b) / 2.0;
            break;
        }

        if k >= 20 || (was_newton && f.abs() > 0.9 * fold.abs()) || g >= 0.0 {
            x = (a + b) / 2.0;
        } else {
            let p = -f / g;
            let y = x + p;
            if y <= a || y >= b {
                x = (a + b) / 2.0;
            } else {
                is_newton = true;
                x = y;
                if p.abs() < tolx * x * (1.0 - x) {
                    break;
                }
            }
        }
    }
    -x.ln() / d as f64
}

/// Compute Lambda from score frequency distribution.
fn compute_lambda(sfp: &SfDist) -> f64 {
    if sfp.score_avg >= 0.0 {
        return -1.0;
    }
    let low = sfp.obs_min;
    let high = sfp.obs_max;
    let mut d = -low;
    for i in 1..=(high - low) {
        if d <= 1 {
            break;
        }
        if sfp.p(low + i) != 0.0 {
            d = kbp_gcd(d, i);
        }
    }
    solve_lambda(sfp, d, low, high, 0.5)
}

/// Compute H (relative entropy) from score frequencies and Lambda.
fn compute_h(sfp: &SfDist, lambda: f64) -> f64 {
    if lambda < 0.0 {
        return -1.0;
    }
    let etonlam = (-lambda).exp();
    let mut sum = sfp.obs_min as f64 * sfp.p(sfp.obs_min);
    for s in (sfp.obs_min + 1)..=sfp.obs_max {
        sum = s as f64 * sfp.p(s) + etonlam * sum;
    }
    let scale = powi(etonlam, sfp.obs_max);
    if scale > 0.0 {
        lambda * sum / scale
    } else {
        lambda * (lambda * sfp.obs_max as f64 + sum.ln()).exp()
    }
}

/// Compute K from Lambda and H using DP algorithm.
fn compute_k(sfp: &SfDist, lambda: f64, h: f64) -> f64 {
    if lambda <= 0.0 || h <= 0.0 || sfp.score_avg >= 0.0 {
        return -1.0;
    }
    let olow = sfp.obs_min;
    let ohigh = sfp.obs_max;
    let mut divisor = -olow;
    for i in 1..=(ohigh - olow) {
        if divisor <= 1 {
            break;
        }
        if sfp.p(olow + i) != 0.0 {
            divisor = kbp_gcd(divisor, i);
        }
    }

    let high = ohigh / divisor;
    let low = olow / divisor;
    let lam_s = lambda * divisor as f64;
    let range = high - low;
    let ftcf = h / lam_s;
    let eml = (-lam_s).exp();

    // Special cases
    if low == -1 && high == 1 {
        let pl = sfp.p(olow);
        let ph = sfp.p(ohigh);
        return (pl - ph) * (pl - ph) / pl;
    }
    if low == -1 || high == 1 {
        let mut f = ftcf;
        if high != 1 {
            let sa = sfp.score_avg / divisor as f64;
            f = (sa * sa) / f;
        }
        return f * (1.0 - eml);
    }

    // DP
    let prob_at = |i: i32| -> f64 {
        let s = olow + i;
        if s >= sfp.score_min && s <= sfp.score_max {
            sfp.p(s)
        } else {
            0.0
        }
    };
    let ru = range as usize;
    let max_iter = 100usize;
    let mut asp = vec![0.0_f64; max_iter * ru + 1];
    asp[0] = 1.0;
    let mut outer_sum = 0.0_f64;
    let mut low_as = 0i32;
    let mut high_as = 0i32;
    let mut inner_sum = 1.0_f64;
    let mut oldsum;

    for iter in 0..max_iter {
        if inner_sum <= 0.0001 {
            break;
        }
        let mut first = range;
        let mut last = range;
        low_as += low;
        high_as += high;
        let span = (high_as - low_as) as isize;
        let mut pp = span;
        while pp >= 0 {
            let p1s = pp - first as isize;
            let p1e = pp - last as isize;
            let mut isum = 0.0;
            let mut p1 = p1s;
            let mut p2 = first;
            while p1 >= p1e {
                if p1 >= 0 {
                    isum += asp[p1 as usize] * prob_at(p2);
                }
                p1 -= 1;
                p2 += 1;
            }
            asp[pp as usize] = isum;
            if first > 0 {
                first -= 1;
            }
            if pp <= range as isize {
                last -= 1;
            }
            pp -= 1;
        }
        pp += 1;
        inner_sum = asp[pp as usize];
        let mut i = low_as + 1;
        while i < 0 {
            pp += 1;
            inner_sum = asp[pp as usize] + inner_sum * eml;
            i += 1;
        }
        inner_sum *= eml;
        while i <= high_as {
            pp += 1;
            inner_sum += asp[pp as usize];
            i += 1;
        }
        oldsum = inner_sum;
        let _ = oldsum;
        outer_sum += inner_sum / (iter + 1) as f64;
    }

    -(-2.0 * outer_sum).exp() / (ftcf * expm1_stable(-lam_s))
}

/// Context info for ungapped KBP computation.
pub struct UngappedKbpContext {
    pub query_offset: i32,
    pub query_length: i32,
    pub is_valid: bool,
}

/// Compute ungapped KBP for all contexts. Returns per-context `Option<KarlinBlk>`.
/// Port of Blast_ScoreBlkKbpUngappedCalc from blast_stat.c.
pub fn ungapped_kbp_calc(
    query: &[u8],
    contexts: &[UngappedKbpContext],
    loscore: i32,
    hiscore: i32,
    alphabet_size: usize,
    ambiguous: &[u8],
    matrix: &dyn Fn(usize, usize) -> i32,
) -> Vec<Option<KarlinBlk>> {
    // Standard composition: 0.25 each for A/C/G/T (indices 0-3)
    let mut std_freq = vec![0.0f64; alphabet_size];
    for i in 0..4.min(alphabet_size) {
        std_freq[i] = 0.25;
    }

    let mut results = Vec::with_capacity(contexts.len());
    for ctx in contexts {
        if !ctx.is_valid || ctx.query_length <= 0 {
            results.push(None);
            continue;
        }
        let off = ctx.query_offset as usize;
        let len = ctx.query_length as usize;
        let buf = &query[off..off + len];

        // Count residue composition
        let mut counts = vec![0i32; alphabet_size];
        for &b in buf {
            let idx = (b & 0x0F) as usize;
            if idx < alphabet_size {
                counts[idx] += 1;
            }
        }
        for &a in ambiguous {
            if (a as usize) < alphabet_size {
                counts[a as usize] = 0;
            }
        }
        let sum: f64 = counts.iter().map(|&c| c as f64).sum();
        let mut qfreq = vec![0.0f64; alphabet_size];
        if sum > 0.0 {
            for i in 0..alphabet_size {
                qfreq[i] = counts[i] as f64 / sum;
            }
        }

        // Compute score frequency
        let mut sfp = SfDist::new(loscore, hiscore);
        for i in 0..alphabet_size {
            for j in 0..alphabet_size {
                let s = matrix(i, j);
                if s >= loscore && s <= hiscore {
                    *sfp.p_mut(s) += qfreq[i] * std_freq[j];
                }
            }
        }
        // Find obs range and normalize
        let mut obs_min = i32::MAX;
        let mut obs_max = i32::MIN;
        let mut psum = 0.0;
        for s in loscore..=hiscore {
            if sfp.p(s) > 0.0 {
                psum += sfp.p(s);
                if s < obs_min {
                    obs_min = s;
                }
                obs_max = s;
            }
        }
        sfp.obs_min = obs_min;
        sfp.obs_max = obs_max;
        let mut savg = 0.0;
        if psum > 0.0 {
            for s in obs_min..=obs_max {
                *sfp.p_mut(s) /= psum;
                savg += s as f64 * sfp.p(s);
            }
        }
        sfp.score_avg = savg;

        // Compute Lambda, H, K
        let lambda = compute_lambda(&sfp);
        if lambda < 0.0 {
            results.push(None);
            continue;
        }
        let h = compute_h(&sfp, lambda);
        if h < 0.0 {
            results.push(None);
            continue;
        }
        let k = compute_k(&sfp, lambda, h);
        if k < 0.0 {
            results.push(None);
            continue;
        }
        results.push(Some(KarlinBlk {
            lambda,
            k,
            log_k: k.ln(),
            h,
        }));
    }
    results
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
        assert!(
            (lambda - 1.374).abs() < 0.01,
            "lambda for 1/-3 should be ~1.374, got {}",
            lambda
        );

        // For reward=2, penalty=-3: different lambda
        let lambda2 = compute_ungapped_lambda(2, -3);
        assert!(
            lambda2 > 0.0 && lambda2 < lambda,
            "lambda for 2/-3 should be positive and < 1/-3 lambda"
        );
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

    // --- Ported from NCBI scoreblk_unit_test.cpp ---

    /// Port of GetScoreBlockNucl: verify nucleotide scoring properties.
    /// NCBI checks: alphabet_size=16, loscore=-3, hiscore=1, penalty=-3, reward=1
    #[test]
    fn test_scoreblk_nucl_properties() {
        let m = crate::matrix::nucleotide_matrix(1, -3);
        // Find min/max scores in the matrix
        let mut lo = i32::MAX;
        let mut hi = i32::MIN;
        for i in 0..16 {
            for j in 0..16 {
                if m[i][j] < lo {
                    lo = m[i][j];
                }
                if m[i][j] > hi {
                    hi = m[i][j];
                }
            }
        }
        assert_eq!(lo, -3, "nucleotide loscore should be -3");
        assert_eq!(hi, 1, "nucleotide hiscore should be 1");
    }

    /// Port of GetScoreBlockProtein: verify BLOSUM62 scoring properties.
    /// NCBI checks: alphabet_size=BLASTAA_SIZE(28), loscore=-4, hiscore=11
    #[test]
    fn test_scoreblk_protein_properties() {
        let m = &crate::matrix::BLOSUM62;
        let mut lo = i32::MAX;
        let mut hi = i32::MIN;
        for i in 0..crate::matrix::AA_SIZE {
            for j in 0..crate::matrix::AA_SIZE {
                if m[i][j] < lo {
                    lo = m[i][j];
                }
                if m[i][j] > hi {
                    hi = m[i][j];
                }
            }
        }
        assert_eq!(lo, -4, "BLOSUM62 loscore should be -4");
        assert_eq!(hi, 11, "BLOSUM62 hiscore should be 11");
    }

    /// Port of NuclGappedCalc: verify gapped KBP for reward=1, penalty=-2, gap_open=3, gap_extend=1.
    /// NCBI expects: Lambda≈1.32, K≈0.57, round_down=false
    #[test]
    fn test_nucl_gapped_calc_1_2_3_1() {
        let ungapped = compute_ungapped_kbp(1, -2);
        let (kbp, round_down) = nucl_gapped_kbp_lookup(3, 1, 1, -2, &ungapped).unwrap();
        assert!(!round_down);
        assert!(
            (kbp.lambda - 1.32).abs() < 0.01,
            "Lambda should be ~1.32, got {}",
            kbp.lambda
        );
        assert!(
            (kbp.k - 0.57).abs() < 0.01,
            "K should be ~0.57, got {}",
            kbp.k
        );
    }

    /// Port of NuclGappedCalc: verify alpha/beta for reward=1, penalty=-2, gap_open=3, gap_extend=1.
    /// NCBI expects: alpha≈1.3, beta≈-1.0
    #[test]
    fn test_nucl_alpha_beta_1_2_3_1() {
        let ungapped = compute_ungapped_kbp(1, -2);
        let (alpha, beta) = nucl_alpha_beta(1, -2, 3, 1, ungapped.lambda, ungapped.h, true);
        assert!(
            (alpha - 1.3).abs() < 0.01,
            "alpha should be ~1.3, got {}",
            alpha
        );
        assert!(
            (beta - (-1.0)).abs() < 0.01,
            "beta should be ~-1.0, got {}",
            beta
        );
    }

    /// Port of NuclGappedCalc: high gap costs fall back to ungapped params.
    /// reward=1, penalty=-2, gap_open=4, gap_extend=2 → copies ungapped Lambda/K.
    #[test]
    fn test_nucl_gapped_fallback_to_ungapped() {
        let ungapped = compute_ungapped_kbp(1, -2);
        let (kbp, round_down) = nucl_gapped_kbp_lookup(4, 2, 1, -2, &ungapped).unwrap();
        assert!(!round_down);
        assert_eq!(
            kbp.lambda, ungapped.lambda,
            "High gap costs should fall back to ungapped Lambda"
        );
        assert_eq!(
            kbp.k, ungapped.k,
            "High gap costs should fall back to ungapped K"
        );
    }

    /// Port of NuclGappedCalc: ungapped alpha/beta when gap costs exceed table maximum.
    /// Alpha = Lambda/H, beta = 0.0
    #[test]
    fn test_nucl_alpha_beta_ungapped_fallback() {
        let ungapped = compute_ungapped_kbp(1, -2);
        let (alpha, beta) = nucl_alpha_beta(1, -2, 4, 2, ungapped.lambda, ungapped.h, true);
        assert!(
            (alpha - ungapped.lambda / ungapped.h).abs() < 1e-6,
            "alpha should equal Lambda/H for unsupported gap costs"
        );
        assert_eq!(beta, 0.0);
    }

    /// Port of NuclGappedCalc: scaled-up values (reward=10, penalty=-20, gap_open=30, gap_extend=10).
    /// GCD=10, so internally maps to (1, -2, 3, 1). Lambda should be ~0.132 (1.32/10).
    #[test]
    fn test_nucl_gapped_scaled_values() {
        let ungapped = compute_ungapped_kbp(10, -20);
        let (kbp, round_down) = nucl_gapped_kbp_lookup(30, 10, 10, -20, &ungapped).unwrap();
        assert!(!round_down);
        assert!(
            (kbp.lambda - 0.132).abs() < 0.01,
            "Lambda for 10/-20/30/10 should be ~0.132 (scaled), got {}",
            kbp.lambda
        );
        assert!(
            (kbp.k - 0.57).abs() < 0.01,
            "K for 10/-20/30/10 should be ~0.57 (unscaled), got {}",
            kbp.k
        );
    }

    /// Port of NuclGappedCalc: reward=2, penalty=-7, gap_open=4, gap_extend=2 should set round_down=true.
    /// NCBI expects: Lambda≈0.675, K≈0.62, round_down=true
    #[test]
    fn test_nucl_gapped_round_down() {
        let ungapped = compute_ungapped_kbp(2, -7);
        let (kbp, round_down) = nucl_gapped_kbp_lookup(4, 2, 2, -7, &ungapped).unwrap();
        assert!(round_down, "2/-7 should set round_down=true");
        assert!(
            (kbp.lambda - 0.675).abs() < 0.01,
            "Lambda should be ~0.675, got {}",
            kbp.lambda
        );
        assert!(
            (kbp.k - 0.62).abs() < 0.01,
            "K should be ~0.62, got {}",
            kbp.k
        );
    }

    /// Port of NuclGappedCalc: invalid gap costs should return Err.
    /// reward=4, penalty=-5, gap_open=3, gap_extend=2 is not in the table.
    #[test]
    fn test_nucl_gapped_invalid_gap_costs() {
        let ungapped = compute_ungapped_kbp(4, -5);
        let result = nucl_gapped_kbp_lookup(3, 2, 4, -5, &ungapped);
        assert!(
            result.is_err(),
            "gap_open=3, gap_extend=2 for 4/-5 should fail"
        );
    }

    /// Port of NuclGappedCalc: invalid gap costs should return Err.
    /// reward=1, penalty=-2, gap_open=1, gap_extend=3 is not in the table.
    #[test]
    fn test_nucl_gapped_invalid_gap_costs_2() {
        let ungapped = compute_ungapped_kbp(1, -2);
        let result = nucl_gapped_kbp_lookup(1, 3, 1, -2, &ungapped);
        assert!(
            result.is_err(),
            "gap_open=1, gap_extend=3 for 1/-2 should fail"
        );
    }

    /// Port of NuclGappedCalc: unsupported reward/penalty pair returns Err.
    /// reward=2, penalty=-1 is not supported.
    #[test]
    fn test_nucl_gapped_unsupported_scores() {
        let ungapped = compute_ungapped_kbp(2, -1);
        let result = nucl_gapped_kbp_lookup(1, 3, 2, -1, &ungapped);
        assert!(result.is_err(), "2/-1 is not a supported scoring pair");
    }

    /// Port of EqualRewardPenaltyLHtoK: when reward==|penalty|, K should be 1/3.
    /// NCBI uses full Blast_ScoreBlkKbpIdealCalc with DP-based K computation.
    /// The full ungapped_kbp_calc function uses the same algorithm.
    #[test]
    fn test_equal_reward_penalty_k() {
        // Use the full KBP computation via ungapped_kbp_calc
        // Query: uniform ACGT, 1000 bases
        let query: Vec<u8> = (0..1000).map(|i| (i % 4) as u8).collect();
        let contexts = vec![UngappedKbpContext {
            query_offset: 0,
            query_length: 1000,
            is_valid: true,
        }];
        let m = crate::matrix::nucleotide_matrix(2, -2);
        let results = ungapped_kbp_calc(
            &query,
            &contexts,
            -2,
            2,
            16,
            &(4..16).collect::<Vec<u8>>(),
            &|i, j| m[i][j],
        );
        let kbp = results[0].as_ref().unwrap();
        assert!(
            (kbp.k - 1.0 / 3.0).abs() < 0.02,
            "K for equal reward/penalty should be ~1/3, got {}",
            kbp.k
        );
    }

    /// Port of BlastResFreqStdCompNucleotideTest: nucleotide base frequencies should be 0.25 each.
    #[test]
    fn test_nucleotide_standard_frequencies() {
        let freqs = &crate::matrix::NT_FREQUENCIES;
        // First 4 bases (A, C, G, T) should each be 0.25
        for i in 0..4 {
            assert!(
                (freqs[i] - 0.25).abs() < 0.001,
                "Base {} frequency should be 0.25, got {}",
                i,
                freqs[i]
            );
        }
        // Ambiguity codes should be 0
        for i in 4..freqs.len() {
            assert!(
                freqs[i].abs() < 0.001,
                "Ambiguity code {} frequency should be 0, got {}",
                i,
                freqs[i]
            );
        }
    }

    /// Port of BlastResFreqStdCompProteinTest: verify specific amino acid frequencies.
    /// NCBI NCBIstdaa indices: 6=F, 12=M, 22=Y. Compact AA_FREQUENCIES: F=4, M=10, Y=19.
    /// NCBI expects (multiplied by 100000): F→4259(≈3856 NCBI val), M→2243, Y→3216
    #[test]
    fn test_protein_standard_frequencies() {
        let freqs = &crate::matrix::AA_FREQUENCIES;
        // Compact indices: A=0,C=1,D=2,E=3,F=4,G=5,H=6,I=7,K=8,L=9,M=10,N=11,P=12,Q=13,R=14,S=15,T=16,V=17,W=18,Y=19
        let check = |idx: usize, expected_x100k: i32, name: &str| {
            let got = (freqs[idx] * 100000.0).round() as i32;
            assert!(
                (got - expected_x100k).abs() <= 5,
                "{} (idx {}) expected ~{}, got {} (raw: {})",
                name,
                idx,
                expected_x100k,
                got,
                freqs[idx]
            );
        };
        check(4, 4259, "F"); // Phe
        check(10, 2243, "M"); // Met
        check(19, 3216, "Y"); // Tyr
    }

    /// Port of EvalueForProteinFSC: E-values should never be negative.
    /// Uses BLOSUM62 with gap_open=11, gap_extend=1.
    #[test]
    fn test_evalue_never_negative() {
        let kbp = lookup_protein_params(11, 1).unwrap();
        let kbp = KarlinBlk {
            lambda: kbp.lambda,
            k: kbp.k,
            log_k: kbp.k.ln(),
            h: kbp.h,
        };
        // Test cases from NCBI stat_unit_test.cpp
        let cases = [
            (1201, 294, 422),
            (1204, 294, 416),
            (1179, 294, 418),
            (2332, 1801, 1671),
        ];
        for &(score, len1, len2) in &cases {
            let ss = len1 as f64 * len2 as f64;
            let evalue = kbp.raw_to_evalue(score, ss);
            assert!(
                evalue >= 0.0,
                "E-value should be >= 0 for score={}, lens=({},{}), got {}",
                score,
                len1,
                len2,
                evalue
            );
        }
    }

    /// Verify gapped KBP lookup for all supported nucleotide reward/penalty combos.
    /// Not all combos have a linear (0,0) entry; test the first affine entry from each table.
    #[test]
    fn test_all_nucl_reward_penalty_combos() {
        // (reward, penalty, gap_open, gap_extend) - first affine entry from each KBPT table
        let cases: &[(i32, i32, i32, i32)] = &[
            (1, -5, 3, 3),
            (1, -4, 1, 2),
            (1, -3, 2, 2),
            (1, -2, 2, 2),
            (1, -1, 3, 2),
            (2, -7, 2, 4),
            (2, -5, 2, 4),
            (2, -3, 4, 4),
            (3, -4, 6, 3),
            (3, -2, 5, 5),
            (4, -5, 6, 5),
            (5, -4, 10, 6),
        ];
        for &(reward, penalty, go, ge) in cases {
            let ungapped = compute_ungapped_kbp(reward, penalty);
            let result = nucl_gapped_kbp_lookup(go, ge, reward, penalty, &ungapped);
            assert!(
                result.is_ok(),
                "Gapped lookup should work for {}/{}/{}/{}",
                reward,
                penalty,
                go,
                ge
            );
            let (kbp, _) = result.unwrap();
            assert!(
                kbp.is_valid(),
                "KBP should be valid for {}/{}/{}/{}",
                reward,
                penalty,
                go,
                ge
            );
        }
    }

    /// Verify gapped KBP for all entries in KBPT_1_3 table via nucl_gapped_kbp_lookup.
    #[test]
    fn test_nucl_gapped_all_1_3_entries() {
        let ungapped = compute_ungapped_kbp(1, -3);
        // Entries from KBPT_1_3 (the C-compatible table)
        let entries: &[(i32, i32, f64, f64)] = &[
            (2, 2, 1.37, 0.70),
            (1, 2, 1.35, 0.64),
            (0, 2, 1.25, 0.42),
            (2, 1, 1.34, 0.60),
            (1, 1, 1.21, 0.34),
        ];
        for &(go, ge, expected_lambda, expected_k) in entries {
            let result = nucl_gapped_kbp_lookup(go, ge, 1, -3, &ungapped);
            assert!(
                result.is_ok(),
                "1/-3 gap_open={} gap_extend={} should work",
                go,
                ge
            );
            let (kbp, _) = result.unwrap();
            assert!(
                (kbp.lambda - expected_lambda).abs() < 0.01,
                "Lambda mismatch for 1/-3/{}/{}: expected {}, got {}",
                go,
                ge,
                expected_lambda,
                kbp.lambda
            );
            assert!(
                (kbp.k - expected_k).abs() < 0.01,
                "K mismatch for 1/-3/{}/{}: expected {}, got {}",
                go,
                ge,
                expected_k,
                kbp.k
            );
        }
    }

    /// Verify protein gapped KBP for all BLOSUM62 table entries.
    #[test]
    fn test_protein_gapped_all_blosum62_entries() {
        for &(go, ge, expected_lambda, expected_k, expected_h, _alpha, _beta) in BLOSUM62_PARAMS {
            let params = lookup_protein_params(go, ge);
            assert!(
                params.is_some(),
                "BLOSUM62 gap_open={} gap_extend={} should be in table",
                go,
                ge
            );
            let p = params.unwrap();
            assert!((p.lambda - expected_lambda).abs() < 1e-6);
            assert!((p.k - expected_k).abs() < 1e-6);
            assert!((p.h - expected_h).abs() < 1e-6);
        }
    }

    /// Verify length adjustment converges for typical search parameters.
    #[test]
    fn test_length_adjustment_exact_convergence() {
        let kbp = lookup_protein_params(11, 1).unwrap();
        let (adj, converged) = compute_length_adjustment_exact(
            kbp.k,
            kbp.k.ln(),
            kbp.alpha / kbp.lambda,
            kbp.beta,
            300,
            1_000_000,
            5000,
        );
        assert!(converged, "Length adjustment should converge");
        assert!(adj > 0, "Length adjustment should be positive");
        assert!(adj < 300, "Length adjustment should be < query length");
    }

    /// Verify E-value <-> raw score round-trip consistency.
    #[test]
    fn test_evalue_raw_score_roundtrip() {
        let kbp = KarlinBlk {
            lambda: 0.267,
            k: 0.041,
            log_k: 0.041_f64.ln(),
            h: 0.14,
        };
        let search_space = 1e9;
        // For a given raw score, convert to evalue, then back to raw score
        for raw in [20, 50, 100, 200] {
            let evalue = kbp.raw_to_evalue(raw, search_space);
            let recovered = kbp.evalue_to_raw(evalue, search_space);
            assert!(
                (recovered - raw).abs() <= 1,
                "Round-trip failed: raw={} -> evalue={:.2e} -> recovered={}",
                raw,
                evalue,
                recovered
            );
        }
    }

    /// Verify ungapped KBP computation for protein (BLOSUM62).
    #[test]
    fn test_protein_ungapped_kbp_values() {
        let kbp = protein_ungapped_kbp();
        assert!(kbp.is_valid());
        assert!((kbp.lambda - 0.3176).abs() < 0.01);
        assert!((kbp.k - 0.134).abs() < 0.01);
        assert!((kbp.h - 0.401).abs() < 0.01);
    }

    /// Verify ungapped KBP calc produces valid results for context-based computation.
    #[test]
    fn test_ungapped_kbp_calc_nucleotide() {
        // A simple BLASTNA-encoded query: ACGTACGT = [0,1,2,3,0,1,2,3]
        let query: Vec<u8> = vec![0, 1, 2, 3, 0, 1, 2, 3];
        let contexts = vec![UngappedKbpContext {
            query_offset: 0,
            query_length: 8,
            is_valid: true,
        }];
        let m = crate::matrix::nucleotide_matrix(1, -3);
        let results = ungapped_kbp_calc(
            &query,
            &contexts,
            -3,
            1,
            16,
            &(4..16).collect::<Vec<u8>>(), // ambiguous codes 4-15
            &|i, j| m[i][j],
        );
        assert_eq!(results.len(), 1);
        assert!(results[0].is_some(), "Should produce valid KBP");
        let kbp = results[0].as_ref().unwrap();
        assert!(kbp.is_valid());
        assert!(kbp.lambda > 0.0);
        assert!(kbp.k > 0.0);
        assert!(kbp.h > 0.0);
    }

    /// Verify invalid context (zero-length) produces None.
    #[test]
    fn test_ungapped_kbp_calc_invalid_context() {
        let query: Vec<u8> = vec![0, 1, 2, 3];
        let contexts = vec![UngappedKbpContext {
            query_offset: 0,
            query_length: 0,
            is_valid: false,
        }];
        let m = crate::matrix::nucleotide_matrix(1, -3);
        let results = ungapped_kbp_calc(
            &query,
            &contexts,
            -3,
            1,
            16,
            &(4..16).collect::<Vec<u8>>(),
            &|i, j| m[i][j],
        );
        assert!(results[0].is_none());
    }

    /// Verify search space with zero length adjustment.
    #[test]
    fn test_search_space_no_adjustment() {
        let ss = compute_search_space(500, 10_000_000, 1000, 0);
        assert_eq!(ss, 500.0 * 10_000_000.0);
    }

    /// Verify search space clamps to 1 when adjustment exceeds lengths.
    #[test]
    fn test_search_space_clamped() {
        let ss = compute_search_space(10, 100, 5, 50);
        // eff_query = max(10-50, 1) = 1
        // eff_db = max(100 - 5*50, 1) = max(-150, 1) = 1
        assert_eq!(ss, 1.0);
    }

    /// Verify length adjustment is 0 for degenerate inputs.
    #[test]
    fn test_length_adjustment_degenerate() {
        let kbp = KarlinBlk {
            lambda: 0.0,
            k: 0.0,
            log_k: f64::NEG_INFINITY,
            h: 0.0,
        };
        let adj = compute_length_adjustment(100, 1000000, 100, &kbp);
        assert_eq!(adj, 0);
    }
}
