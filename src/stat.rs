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
/// Format: (gap_open, gap_extend, lambda, K, H, alpha, beta)
pub const BLOSUM62_PARAMS: &[(i32, i32, f64, f64, f64, f64, f64)] = &[
    (11, 1, 0.267, 0.041, 0.14, 1.24, -0.70),
    (10, 1, 0.291, 0.075, 0.23, 1.00, -0.55),
    (9, 1, 0.305, 0.10, 0.30, 0.89, -0.46),
    (8, 2, 0.270, 0.047, 0.15, 1.18, -0.67),
    (7, 2, 0.290, 0.070, 0.22, 1.02, -0.57),
    (6, 2, 0.305, 0.097, 0.29, 0.90, -0.48),
    (11, 2, 0.235, 0.020, 0.073, 1.55, -0.91),
    (10, 2, 0.256, 0.032, 0.11, 1.30, -0.76),
    (9, 2, 0.271, 0.044, 0.15, 1.16, -0.67),
];

/// Look up gapped KBP for protein scoring (BLOSUM62).
pub fn lookup_protein_params(gap_open: i32, gap_extend: i32) -> Option<GappedParams> {
    for &(go, ge, lambda, k, h, alpha, beta) in BLOSUM62_PARAMS {
        if go == gap_open && ge == gap_extend {
            return Some(GappedParams { gap_open: go, gap_extend: ge, lambda, k, h, alpha, beta });
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
const KBPT_3_2: &[KbpTableRow] = &[
    [5.0, 5.0, 0.208, 0.030, 0.072, 2.9, -47.0, 77.0],
];
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
    if b > a { std::mem::swap(&mut a, &mut b); }
    while b != 0 { let c = a % b; a = b; b = c; }
    a as i32
}

fn get_kbp_table(reward: i32, penalty: i32) -> Option<KbpTableMeta> {
    match (reward, penalty) {
        (1, -5) => Some(KbpTableMeta { table: KBPT_1_5, gap_open_max: 3, gap_extend_max: 3, round_down: false }),
        (1, -4) => Some(KbpTableMeta { table: KBPT_1_4, gap_open_max: 2, gap_extend_max: 2, round_down: false }),
        (2, -7) => Some(KbpTableMeta { table: KBPT_2_7, gap_open_max: 4, gap_extend_max: 4, round_down: true }),
        (1, -3) => Some(KbpTableMeta { table: KBPT_1_3, gap_open_max: 2, gap_extend_max: 2, round_down: false }),
        (2, -5) => Some(KbpTableMeta { table: KBPT_2_5, gap_open_max: 4, gap_extend_max: 4, round_down: true }),
        (1, -2) => Some(KbpTableMeta { table: KBPT_1_2, gap_open_max: 2, gap_extend_max: 2, round_down: false }),
        (2, -3) => Some(KbpTableMeta { table: KBPT_2_3, gap_open_max: 6, gap_extend_max: 4, round_down: true }),
        (3, -4) => Some(KbpTableMeta { table: KBPT_3_4, gap_open_max: 6, gap_extend_max: 3, round_down: true }),
        (1, -1) => Some(KbpTableMeta { table: KBPT_1_1, gap_open_max: 4, gap_extend_max: 2, round_down: false }),
        (3, -2) => Some(KbpTableMeta { table: KBPT_3_2, gap_open_max: 5, gap_extend_max: 5, round_down: false }),
        (4, -5) => Some(KbpTableMeta { table: KBPT_4_5, gap_open_max: 12, gap_extend_max: 8, round_down: false }),
        (5, -4) => Some(KbpTableMeta { table: KBPT_5_4, gap_open_max: 25, gap_extend_max: 10, round_down: false }),
        _ => None,
    }
}

/// Look up gapped KBP params for nucleotide, matching C's Blast_KarlinBlkNuclGappedCalc exactly.
/// Returns Ok((lambda, k, log_k, h, round_down)) or Err message.
pub fn nucl_gapped_kbp_lookup(
    gap_open: i32, gap_extend: i32,
    reward: i32, penalty: i32,
    ungapped: &KarlinBlk,
) -> Result<(KarlinBlk, bool), String> {
    let divisor = kbp_gcd(reward, penalty.abs());
    let (nr, np) = (reward / divisor, penalty / divisor);

    let meta = get_kbp_table(nr, np)
        .ok_or_else(|| format!("Unsupported scores {} {}", reward, penalty))?;

    let round_down = meta.round_down;

    // Split: first row with gap_open==0 && gap_extend==0 is the linear entry
    let (affine, linear) = if !meta.table.is_empty()
        && meta.table[0][0] == 0.0 && meta.table[0][1] == 0.0 {
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
            return Ok((KarlinBlk { lambda: lam, k, log_k: k.ln(), h }, round_down));
        }
    }

    // Search affine entries
    for row in affine {
        let (go, ge, lam, k, h) = scale(row, divisor);
        if go as i32 == gap_open && ge as i32 == gap_extend {
            return Ok((KarlinBlk { lambda: lam, k, log_k: k.ln(), h }, round_down));
        }
    }

    // Fallback: gap costs exceed table maximum → use ungapped params
    if gap_open >= go_max && gap_extend >= ge_max {
        return Ok((ungapped.clone(), round_down));
    }

    Err(format!("Unsupported gap costs {} {} for scores {} {}", gap_open, gap_extend, reward, penalty))
}

/// Compute the length adjustment for effective search space (exact C-compatible).
/// Port of BLAST_ComputeLengthAdjustment from blast_stat.c.
/// Returns (length_adjustment, converged).
pub fn compute_length_adjustment_exact(
    k: f64, log_k: f64,
    alpha_d_lambda: f64, beta: f64,
    query_length: i32, db_length: i64, db_num_seqs: i32,
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
            if ell_min == ell_max { break; }
        } else {
            ell_max = ell;
        }

        if ell_min <= ell_bar && ell_bar <= ell_max {
            ell_next = ell_bar;
        } else {
            ell_next = if i == 1 { ell_max } else { (ell_min + ell_max) / 2.0 };
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
    reward: i32, penalty: i32,
    gap_open: i32, gap_extend: i32,
    ungapped_lambda: f64, ungapped_h: f64,
    gapped: bool,
) -> (f64, f64) {
    if !gapped {
        let beta = if (reward == 1 && penalty == -1) || (reward == 2 && penalty == -3) { -2.0 } else { 0.0 };
        return (ungapped_lambda / ungapped_h, beta);
    }

    let divisor = kbp_gcd(reward, penalty.abs());
    let (nr, np) = (reward / divisor, penalty / divisor);

    if let Some(meta) = get_kbp_table(nr, np) {
        let (affine, linear) = if !meta.table.is_empty()
            && meta.table[0][0] == 0.0 && meta.table[0][1] == 0.0 {
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
    let beta = if (reward == 1 && penalty == -1) || (reward == 2 && penalty == -3) { -2.0 } else { 0.0 };
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
        SfDist { score_min: lo, score_max: hi, obs_min: 0, obs_max: 0, score_avg: 0.0, probs: vec![0.0; n] }
    }
    #[inline] fn p(&self, s: i32) -> f64 { self.probs[(s - self.score_min) as usize] }
    #[inline] fn p_mut(&mut self, s: i32) -> &mut f64 { &mut self.probs[(s - self.score_min) as usize] }
}

fn expm1_stable(x: f64) -> f64 {
    if x.abs() > 0.33 { x.exp() - 1.0 } else { x.exp_m1() }
}

fn powi(x: f64, n: i32) -> f64 {
    if n == 0 { return 1.0; }
    let mut x = if n < 0 { 1.0 / x } else { x };
    let mut n = n.unsigned_abs();
    let mut y = 1.0;
    while n > 0 { if n & 1 != 0 { y *= x; } n /= 2; x *= x; }
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
        while i < 0 { g = x * g + f; f = f * x + sfp.p(i); i += d; }
        g = x * g + f;
        f = f * x + sfp.p(0) - 1.0;
        i = d;
        while i <= high { g = x * g + f; f = f * x + sfp.p(i); i += d; }

        if f > 0.0 { a = x; } else if f < 0.0 { b = x; } else { break; }
        if b - a < 2.0 * a * (1.0 - b) * tolx { x = (a + b) / 2.0; break; }

        if k >= 20 || (was_newton && f.abs() > 0.9 * fold.abs()) || g >= 0.0 {
            x = (a + b) / 2.0;
        } else {
            let p = -f / g;
            let y = x + p;
            if y <= a || y >= b { x = (a + b) / 2.0; }
            else { is_newton = true; x = y; if p.abs() < tolx * x * (1.0 - x) { break; } }
        }
    }
    -x.ln() / d as f64
}

/// Compute Lambda from score frequency distribution.
fn compute_lambda(sfp: &SfDist) -> f64 {
    if sfp.score_avg >= 0.0 { return -1.0; }
    let low = sfp.obs_min;
    let high = sfp.obs_max;
    let mut d = -low;
    for i in 1..=(high - low) { if d <= 1 { break; } if sfp.p(low + i) != 0.0 { d = kbp_gcd(d, i); } }
    solve_lambda(sfp, d, low, high, 0.5)
}

/// Compute H (relative entropy) from score frequencies and Lambda.
fn compute_h(sfp: &SfDist, lambda: f64) -> f64 {
    if lambda < 0.0 { return -1.0; }
    let etonlam = (-lambda).exp();
    let mut sum = sfp.obs_min as f64 * sfp.p(sfp.obs_min);
    for s in (sfp.obs_min + 1)..=sfp.obs_max {
        sum = s as f64 * sfp.p(s) + etonlam * sum;
    }
    let scale = powi(etonlam, sfp.obs_max);
    if scale > 0.0 { lambda * sum / scale } else { lambda * (lambda * sfp.obs_max as f64 + sum.ln()).exp() }
}

/// Compute K from Lambda and H using DP algorithm.
fn compute_k(sfp: &SfDist, lambda: f64, h: f64) -> f64 {
    if lambda <= 0.0 || h <= 0.0 || sfp.score_avg >= 0.0 { return -1.0; }
    let olow = sfp.obs_min;
    let ohigh = sfp.obs_max;
    let mut divisor = -olow;
    for i in 1..=(ohigh - olow) { if divisor <= 1 { break; } if sfp.p(olow + i) != 0.0 { divisor = kbp_gcd(divisor, i); } }

    let high = ohigh / divisor;
    let low = olow / divisor;
    let lam_s = lambda * divisor as f64;
    let range = high - low;
    let ftcf = h / lam_s;
    let eml = (-lam_s).exp();

    // Special cases
    if low == -1 && high == 1 {
        let pl = sfp.p(olow); let ph = sfp.p(ohigh);
        return (pl - ph) * (pl - ph) / pl;
    }
    if low == -1 || high == 1 {
        let mut f = ftcf;
        if high != 1 { let sa = sfp.score_avg / divisor as f64; f = (sa * sa) / f; }
        return f * (1.0 - eml);
    }

    // DP
    let prob_at = |i: i32| -> f64 {
        let s = olow + i;
        if s >= sfp.score_min && s <= sfp.score_max { sfp.p(s) } else { 0.0 }
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
        if inner_sum <= 0.0001 { break; }
        let mut first = range;
        let mut last = range;
        low_as += low; high_as += high;
        let span = (high_as - low_as) as isize;
        let mut pp = span;
        while pp >= 0 {
            let p1s = pp - first as isize;
            let p1e = pp - last as isize;
            let mut isum = 0.0;
            let mut p1 = p1s;
            let mut p2 = first;
            while p1 >= p1e {
                if p1 >= 0 { isum += asp[p1 as usize] * prob_at(p2); }
                p1 -= 1; p2 += 1;
            }
            asp[pp as usize] = isum;
            if first > 0 { first -= 1; }
            if pp <= range as isize { last -= 1; }
            pp -= 1;
        }
        pp += 1;
        inner_sum = asp[pp as usize];
        let mut i = low_as + 1;
        while i < 0 { pp += 1; inner_sum = asp[pp as usize] + inner_sum * eml; i += 1; }
        inner_sum *= eml;
        while i <= high_as { pp += 1; inner_sum += asp[pp as usize]; i += 1; }
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
    loscore: i32, hiscore: i32,
    alphabet_size: usize,
    ambiguous: &[u8],
    matrix: &dyn Fn(usize, usize) -> i32,
) -> Vec<Option<KarlinBlk>> {
    // Standard composition: 0.25 each for A/C/G/T (indices 0-3)
    let mut std_freq = vec![0.0f64; alphabet_size];
    for i in 0..4.min(alphabet_size) { std_freq[i] = 0.25; }

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
        for &b in buf { let idx = (b & 0x0F) as usize; if idx < alphabet_size { counts[idx] += 1; } }
        for &a in ambiguous { if (a as usize) < alphabet_size { counts[a as usize] = 0; } }
        let sum: f64 = counts.iter().map(|&c| c as f64).sum();
        let mut qfreq = vec![0.0f64; alphabet_size];
        if sum > 0.0 { for i in 0..alphabet_size { qfreq[i] = counts[i] as f64 / sum; } }

        // Compute score frequency
        let mut sfp = SfDist::new(loscore, hiscore);
        for i in 0..alphabet_size {
            for j in 0..alphabet_size {
                let s = matrix(i, j);
                if s >= loscore && s <= hiscore { *sfp.p_mut(s) += qfreq[i] * std_freq[j]; }
            }
        }
        // Find obs range and normalize
        let mut obs_min = i32::MAX; let mut obs_max = i32::MIN;
        let mut psum = 0.0;
        for s in loscore..=hiscore {
            if sfp.p(s) > 0.0 { psum += sfp.p(s); if s < obs_min { obs_min = s; } obs_max = s; }
        }
        sfp.obs_min = obs_min; sfp.obs_max = obs_max;
        let mut savg = 0.0;
        if psum > 0.0 {
            for s in obs_min..=obs_max { *sfp.p_mut(s) /= psum; savg += s as f64 * sfp.p(s); }
        }
        sfp.score_avg = savg;

        // Compute Lambda, H, K
        let lambda = compute_lambda(&sfp);
        if lambda < 0.0 { results.push(None); continue; }
        let h = compute_h(&sfp, lambda);
        if h < 0.0 { results.push(None); continue; }
        let k = compute_k(&sfp, lambda, h);
        if k < 0.0 { results.push(None); continue; }
        results.push(Some(KarlinBlk { lambda, k, log_k: k.ln(), h }));
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
