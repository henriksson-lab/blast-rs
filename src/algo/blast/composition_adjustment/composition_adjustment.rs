//! Composition-based statistics adjustment.
//! Verbatim port of NCBI blast_kappa.c + composition_adjustment.c.
//!
//! Implements comp_based_stats modes:
//!   0 = no adjustment
//!   1 = unconditional composition-based lambda adjustment
//!   2 = conditional (only if composition diverges significantly)

// NCBI score/joint-probability tables (BLOSUM62 etc.) are declared with
// >16 significant digits. Keep them byte-identical to the C source; f64
// rounds the same as C so the extra digits are harmless.
#![allow(clippy::excessive_precision)]

use crate::matrix::AA_SIZE;

/// Port of NCBI `COMPO_NUM_TRUE_AA` (`composition_constants.h:46`):
/// count of the 20 standard amino acids (excludes ambiguity codes).
pub const COMPO_NUM_TRUE_AA: usize = 20;

/// Port of NCBI `COMPO_LARGEST_ALPHABET` (`composition_constants.h:51`):
/// 28-symbol NCBIstdaa alphabet size used by composition-adjustment
/// buffers (matches `encoding::BLASTAA_SIZE`).
pub const COMPO_LARGEST_ALPHABET: usize = crate::encoding::BLASTAA_SIZE;
/// NCBI `kLambdaErrorTolerance` (`composition_adjustment.c:66`).
pub const LAMBDA_ERROR_TOLERANCE: f64 = 0.000_000_1;
/// NCBI `kLambdaFunctionTolerance` (`composition_adjustment.c:69`).
pub const LAMBDA_FUNCTION_TOLERANCE: f64 = 0.000_01;
/// NCBI `kLambdaIterationLimit` (`composition_adjustment.c:72`).
pub const LAMBDA_ITERATION_LIMIT: usize = 100;
/// NCBI `COMPO_MIN_LAMBDA` (`unified_pvalues.h:47`).
pub const COMPO_MIN_LAMBDA: f64 = 0.034;
pub const LAMBDA_RATIO_LOWER_BOUND: f64 = 0.5;
/// NCBI `kCompositionMargin` (`composition_adjustment.c:61`).
pub const COMPOSITION_MARGIN: i32 = 20;
/// Rust-specific float sentinel for "impossible score" in composition
/// adjustment's float matrix. NCBI's `COMPO_SCORE_MIN`
/// (`composition_constants.h:43`) is `INT2_MIN = -32768`, but its
/// callers are in integer space. Rust uses a more aggressive value
/// (−100 000) so a single float-scale accumulator comparison reaches
/// it regardless of per-cell additions.
pub const COMPO_SCORE_MIN: f64 = -100000.0;

/// NCBIstdaa indices of the 20 true amino acids.
/// Matches NCBI trueCharPositions[].
pub const TRUE_CHAR_POSITIONS: [usize; COMPO_NUM_TRUE_AA] = {
    let residues = crate::encoding::NCBISTDAA_STANDARD_RESIDUES;
    [
        residues[0] as usize,
        residues[1] as usize,
        residues[2] as usize,
        residues[3] as usize,
        residues[4] as usize,
        residues[5] as usize,
        residues[6] as usize,
        residues[7] as usize,
        residues[8] as usize,
        residues[9] as usize,
        residues[10] as usize,
        residues[11] as usize,
        residues[12] as usize,
        residues[13] as usize,
        residues[14] as usize,
        residues[15] as usize,
        residues[16] as usize,
        residues[17] as usize,
        residues[18] as usize,
        residues[19] as usize,
    ]
};

/// NCBI: `ECompoAdjustModes` (`composition_constants.h`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum ECompoAdjustModes {
    NoCompositionBasedStats = 0,
    CompositionBasedStats = 1,
    CompositionMatrixAdjust = 2,
    CompoForceFullMatrixAdjust = 3,
    NumCompoAdjustModes,
}

/// NCBI: `Blast_AminoAcidComposition` (`composition_adjustment.h:54`).
#[derive(Debug, Clone)]
pub struct BlastAminoAcidComposition {
    pub prob: [f64; COMPO_LARGEST_ALPHABET],
    pub num_true_amino_acids: i32,
}

/// NCBI: `Blast_MatrixInfo` (`composition_adjustment.h:83`).
#[derive(Debug, Clone)]
pub struct BlastMatrixInfo {
    pub matrix_name: Option<String>,
    pub start_matrix: Vec<Vec<i32>>,
    pub start_freq_ratios: Vec<Vec<f64>>,
    pub rows: i32,
    pub cols: i32,
    pub position_based: i32,
    pub ungapped_lambda: f64,
}

/// NCBI: `Blast_CompositionWorkspace` (`composition_adjustment.h:112`).
#[derive(Debug, Clone)]
pub struct BlastCompositionWorkspace {
    pub mat_b: [[f64; COMPO_NUM_TRUE_AA]; COMPO_NUM_TRUE_AA],
    pub mat_final: [[f64; COMPO_NUM_TRUE_AA]; COMPO_NUM_TRUE_AA],
    pub first_standard_freq: [f64; COMPO_NUM_TRUE_AA],
    pub second_standard_freq: [f64; COMPO_NUM_TRUE_AA],
}

/// Read amino acid composition from NCBIstdaa sequence.
/// NCBI: Blast_ReadAaComposition (composition_adjustment.c).
pub fn blast_read_aa_composition(sequence: &[u8], alphsize: usize) -> (Vec<f64>, usize) {
    let mut composition = BlastAminoAcidComposition {
        prob: [0.0; COMPO_LARGEST_ALPHABET],
        num_true_amino_acids: 0,
    };
    let alphsize = alphsize.min(COMPO_LARGEST_ALPHABET);
    let mut num_true_amino_acids = 0i32;

    for i in 0..alphsize {
        composition.prob[i] = 0.0;
    }

    for &residue in sequence {
        let idx = residue as usize;
        if idx < COMPO_LARGEST_ALPHABET && (ALPHA_CONVERT[idx] >= 0 || idx == E_SELENOCYSTEINE) {
            composition.prob[idx] += 1.0;
            num_true_amino_acids += 1;
        }
    }

    if composition.prob[E_SELENOCYSTEINE] > 0.0 {
        composition.prob[E_CCHAR] += composition.prob[E_SELENOCYSTEINE];
        composition.prob[E_SELENOCYSTEINE] = 0.0;
    }

    composition.num_true_amino_acids = num_true_amino_acids;
    if num_true_amino_acids > 0 {
        for i in 0..alphsize {
            composition.prob[i] /= num_true_amino_acids as f64;
        }
    }
    (
        composition.prob[..alphsize].to_vec(),
        composition.num_true_amino_acids.max(0) as usize,
    )
}

/// NCBI: Blast_GetCompositionRange (composition_adjustment.c:1204).
pub fn blast_get_composition_range(subject_data: &[u8], start: i32, finish: i32) -> (i32, i32) {
    let length = subject_data.len() as i32;
    let mut left = start.clamp(0, length);
    let finish = finish.clamp(0, length);

    let mut i = left;
    while i > 0 {
        if subject_data[(i - 1) as usize] == E_STOP_CHAR as u8 {
            if i + COMPOSITION_MARGIN < left {
                left = i + COMPOSITION_MARGIN;
            }
            break;
        }
        i -= 1;
    }
    if i == 0 {
        left = 0;
    }

    let mut right = finish;
    i = right;
    while i < length {
        if subject_data[i as usize] == E_STOP_CHAR as u8 {
            if i - COMPOSITION_MARGIN > right {
                right = i - COMPOSITION_MARGIN;
            }
            break;
        }
        i += 1;
    }
    if i == length {
        right = length;
    }
    (left, right)
}

/// NCBI: s_GetScoreRange (composition_adjustment.c).
pub fn s_get_score_range(matrix: &[[i32; AA_SIZE]; AA_SIZE], rows: usize) -> (i32, i32) {
    let mut obs_min = 0i32;
    let mut obs_max = 0i32;
    for irow in 0..rows {
        for &jcol in &TRUE_CHAR_POSITIONS {
            let s = matrix[irow][jcol];
            if s < obs_min && s > COMPO_SCORE_MIN as i32 {
                obs_min = s;
            }
            if s > obs_max {
                obs_max = s;
            }
        }
    }
    (obs_min, obs_max)
}

/// Compute score probability distribution from matrix and compositions.
/// NCBI: s_GetMatrixScoreProbs (composition_adjustment.c).
/// Computes score probabilities from integer matrix and composition probabilities.
/// Only considers true amino acid columns (trueCharPositions), matching NCBI exactly.
pub fn s_get_matrix_score_probs(
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
    alphsize: usize,
    query_prob: &[f64],
    subject_prob: &[f64],
) -> (Vec<f64>, i32, i32) {
    let (obs_min, obs_max) = s_get_score_range(matrix, alphsize);

    let range = (obs_max - obs_min + 1) as usize;
    let mut score_prob = vec![0.0f64; range];

    // Port of s_GetMatrixScoreProbs: accumulate probabilities
    for irow in 0..alphsize {
        for &jcol in &TRUE_CHAR_POSITIONS {
            if jcol >= alphsize {
                continue;
            }
            let s = matrix[irow][jcol];
            if s >= obs_min {
                let idx = (s - obs_min) as usize;
                score_prob[idx] += query_prob[irow] * subject_prob[jcol];
            }
        }
    }
    (score_prob, obs_min, obs_max)
}

/// Compute lambda for a specific pair of compositions.
/// NCBI: Blast_CalcLambdaFullPrecision (composition_adjustment.c).
/// Uses Newton iteration on f(x) where x = exp(-lambda).
pub fn blast_calc_lambda_full_precision(
    score_matrix: &[[f64; COMPO_NUM_TRUE_AA]; COMPO_NUM_TRUE_AA],
    row_prob: &[f64],
    col_prob: &[f64],
    alphsize: usize,
) -> (f64, usize) {
    let mut f = 4.0f64;
    let mut left = 0.0f64;
    let mut right = 1.0f64;
    let mut x = 0.367879441171f64; // exp(-1)
    let mut is_newton = false;

    // Find max score and average score
    let mut max_score = COMPO_SCORE_MIN;
    let mut avg_score = 0.0f64;
    for i in 0..alphsize.min(COMPO_NUM_TRUE_AA) {
        if row_prob[i] == 0.0 {
            continue;
        }
        for j in 0..alphsize.min(COMPO_NUM_TRUE_AA) {
            if col_prob[j] == 0.0 {
                continue;
            }
            if max_score < score_matrix[i][j] {
                max_score = score_matrix[i][j];
            }
            avg_score += row_prob[i] * col_prob[j] * score_matrix[i][j];
        }
    }
    if max_score <= 0.0 || avg_score >= 0.0 {
        return (-1.0, LAMBDA_ITERATION_LIMIT);
    }

    for k in 0..LAMBDA_ITERATION_LIMIT {
        let fold = f;
        let lambda = -x.ln();
        let was_newton = is_newton;
        is_newton = false;

        let x_pow_max = (-max_score * lambda).exp();
        f = -x_pow_max;
        let mut slope = max_score * f / x;

        for i in 0..alphsize.min(COMPO_NUM_TRUE_AA) {
            if row_prob[i] == 0.0 {
                continue;
            }
            for j in 0..alphsize.min(COMPO_NUM_TRUE_AA) {
                if col_prob[j] == 0.0 {
                    continue;
                }
                let ff;
                if max_score != score_matrix[i][j] {
                    let diff = max_score - score_matrix[i][j];
                    ff = row_prob[i] * col_prob[j] * (-lambda * diff).exp();
                    slope += diff * ff / x;
                } else {
                    ff = row_prob[i] * col_prob[j];
                }
                f += ff;
            }
        }

        if f > 0.0 {
            left = x;
        } else if f < 0.0 {
            right = x;
        } else {
            return (-x.ln(), k);
        }

        if right - left <= 2.0 * left * (1.0 - right) * LAMBDA_ERROR_TOLERANCE
            && (f / x_pow_max).abs() <= LAMBDA_FUNCTION_TOLERANCE
        {
            x = (left + right) / 2.0;
            return (-x.ln(), k);
        }

        if (was_newton && f.abs() > 0.9 * fold.abs()) || slope >= 0.0 {
            x = (left + right) / 2.0;
        } else {
            let p = -f / slope;
            let y = x + p;
            if y <= left || y >= right {
                x = (left + right) / 2.0;
            } else {
                is_newton = true;
                x = y;
                if p.abs() <= LAMBDA_ERROR_TOLERANCE * x * (1.0 - x)
                    && (f / x_pow_max).abs() <= LAMBDA_FUNCTION_TOLERANCE
                {
                    return (-x.ln(), k);
                }
            }
        }
    }
    (-1.0, LAMBDA_ITERATION_LIMIT)
}

/// NCBI: Blast_KarlinLambdaNR (blast_stat.c).
/// Compute lambda from a score probability distribution using
/// Newton-Raphson with Horner's polynomial evaluation.
pub fn blast_karlin_lambda_nr(sprob: &[f64], obs_min: i32, obs_max: i32, lambda0: f64) -> f64 {
    let low = obs_min;
    let high = obs_max;
    let range = (high - low + 1).max(0) as usize;
    let mut score_avg = 0.0f64;
    for i in 0..range {
        score_avg += (low + i as i32) as f64 * sprob.get(i).copied().unwrap_or(0.0);
    }
    if score_avg >= 0.0 || crate::stat::blast_score_chk(low, high) != 0 {
        return -1.0;
    }

    // Find greatest common divisor of all scores with nonzero probability
    let mut d = -low;
    for i in 1..=(high - low) {
        if d <= 1 {
            break;
        }
        if sprob[i as usize] != 0.0 {
            d = crate::math::blast_gcd(d, i);
        }
    }
    if d <= 0 {
        d = 1;
    }

    let x0 = (-lambda0).exp();
    let mut x: f64 = if x0 > 0.0 && x0 < 1.0 { x0 } else { 0.5 };
    let mut a = 0.0f64;
    let mut b = 1.0f64;
    let mut f = 4.0f64;
    let mut is_newton = 0i32;

    // NCBI `Blast_KarlinLambdaNR` (`blast_stat.c:2594`):
    //   `NlmKarlinLambdaNR(sprob, d, low, high, lambda0,
    //                      BLAST_KARLIN_LAMBDA_ACCURACY_DEFAULT,
    //                      20, 20 + BLAST_KARLIN_LAMBDA_ITER_DEFAULT, ...)`.
    // Per NCBI's `NlmKarlinLambdaNR` signature
    // (`blast_stat.c:2491`: `..., Int4 itmax, Int4 maxNewton, ...`),
    // argument 7 is `itmax = 20` (total loop iterations) and argument 8
    // is `maxNewton = 20 + BLAST_KARLIN_LAMBDA_ITER_DEFAULT = 37` (Newton
    // budget). With `maxNewton > itmax`, the `k >= maxNewton` test never
    // fires; bisection only happens when Newton convergence stalls or
    // `g >= 0`. We had these swapped previously — a real bug.
    let tolx = 1.0e-5;
    let itmax = 20;
    let max_newton = itmax + 17; // BLAST_KARLIN_LAMBDA_ITER_DEFAULT = 17

    // sprob is indexed relative to sprob[0] = P(obs_min)
    // sprob[-obs_min] = P(0)
    let sprob_centered = |s: i32| -> f64 {
        let idx = (s - low) as usize;
        if idx < sprob.len() {
            sprob[idx]
        } else {
            0.0
        }
    };

    let mut k = 0;
    while k < itmax {
        let fold = f;
        let was_newton = is_newton;
        is_newton = 0;

        // Horner's rule for evaluating polynomial and derivative
        let mut g = 0.0f64;
        f = sprob_centered(low);
        let mut i = low + d;
        while i < 0 {
            g = x * g + f;
            f = f * x + sprob_centered(i);
            i += d;
        }
        g = x * g + f;
        f = f * x + sprob_centered(0) - 1.0;
        i = d;
        while i <= high {
            g = x * g + f;
            f = f * x + sprob_centered(i);
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

        if k >= max_newton || (was_newton != 0 && f.abs() > 0.9 * fold.abs()) || g >= 0.0 {
            x = (a + b) / 2.0;
        } else {
            let p = -f / g;
            let y = x + p;
            if y <= a || y >= b {
                x = (a + b) / 2.0;
            } else {
                is_newton = 1;
                x = y;
                if p.abs() < tolx * x * (1.0 - x) {
                    break;
                }
            }
        }
        k += 1;
    }
    -x.ln() / d as f64
}

/// Build an integer matrix directly from frequency ratios at a given lambda.
/// NCBI: Blast_Int4MatrixFromFreq (composition_adjustment.c).
pub fn blast_int4_matrix_from_freq(
    lambda: f64,
    freq_ratios: &[[f64; AA_SIZE]; AA_SIZE],
) -> [[i32; AA_SIZE]; AA_SIZE] {
    let mut result = [[0i32; AA_SIZE]; AA_SIZE];
    for i in 0..AA_SIZE {
        for j in 0..AA_SIZE {
            let score = if freq_ratios[i][j] == 0.0 {
                COMPO_SCORE_MIN
            } else {
                freq_ratios[i][j].ln() / lambda
            };
            result[i][j] = if score < i32::MIN as f64 {
                i32::MIN
            } else {
                crate::math::blast_nint(score) as i32
            };
        }
    }
    result
}

/// NCBI: BLAST_KarlinEtoP (blast_stat.c).
pub fn blast_karlin_eto_p(x: f64) -> f64 {
    // NCBI: `return -BLAST_Expm1(-x);` (= 1 - exp(-x) stably).
    -crate::math::expm1(-x)
}

/// NCBI: BLAST_KarlinPtoE (blast_stat.c).
///
/// Mirrors NCBI `BLAST_KarlinPtoE` (`blast_stat.c:4175`) exactly:
/// - `p < 0` or `p > 1` → returns `INT4_MIN` (cast to f64).
/// - `p == 1` → returns `INT4_MAX` (cast to f64).
/// - otherwise → `-log1p(-p)`.
pub fn blast_karlin_pto_e(p: f64) -> f64 {
    if !(0.0..=1.0).contains(&p) {
        return i32::MIN as f64;
    }
    if p == 1.0 {
        return i32::MAX as f64;
    }
    -crate::math::blast_log1p(-p)
}

/// Combine composition P-value with alignment P-value using Fisher's method.
/// NCBI: Blast_Overall_P_Value (unified_pvalues.c).
pub fn blast_overall_p_value(p_comp: f64, p_alignment: f64) -> f64 {
    let product = p_comp * p_alignment;
    if product > 0.0 {
        product * (1.0 - product.ln())
    } else {
        product
    }
}

/// Composition P-value from lambda using the precomputed table.
/// NCBI: Blast_CompositionPvalue (unified_pvalues.c).
pub fn blast_composition_pvalue(lambda: f64) -> f64 {
    /// NCBI `LAMBDA_BIN_SIZE` (`unified_pvalues.c:47`).
    const LAMBDA_BIN_SIZE: f64 = 0.001;
    /// NCBI `LOW_LAMBDA_BIN_CUT` (`unified_pvalues.c:53`).
    const LOW_LAMBDA_BIN_CUT: usize = 35;

    let score = (lambda - COMPO_MIN_LAMBDA) / LAMBDA_BIN_SIZE;
    if score < LOW_LAMBDA_BIN_CUT as f64 {
        return P_LAMBDA_TABLE[LOW_LAMBDA_BIN_CUT];
    }
    if score > P_LAMBDA_TABLE.len() as f64 - 1.0 {
        return 1.0;
    }
    let bin = score as usize;
    if bin >= P_LAMBDA_TABLE.len() - 1 {
        return P_LAMBDA_TABLE[P_LAMBDA_TABLE.len() - 1];
    }
    let delta = score - bin as f64;
    (1.0 - delta) * P_LAMBDA_TABLE[bin] + delta * P_LAMBDA_TABLE[bin + 1]
}

/// NCBI P_lambda_table — precomputed composition P-values.
/// 565 entries, bin size 0.001, starting at COMPO_MIN_LAMBDA=0.034.
#[allow(dead_code)]
pub static P_LAMBDA_TABLE: [f64; 565] = [
    1.247750e-07,
    1.247750e-07,
    1.247750e-07,
    1.247750e-07,
    1.247750e-07,
    1.247750e-07,
    1.247750e-07,
    1.247750e-07,
    1.247750e-07,
    1.247750e-07,
    1.247750e-07,
    1.247750e-07,
    2.495500e-07,
    2.495500e-07,
    2.495500e-07,
    2.495500e-07,
    2.495500e-07,
    2.495500e-07,
    2.495500e-07,
    2.495500e-07,
    2.495500e-07,
    2.495500e-07,
    3.743250e-07,
    4.991000e-07,
    4.991000e-07,
    4.991000e-07,
    4.991000e-07,
    7.486490e-07,
    7.486490e-07,
    7.486490e-07,
    7.486490e-07,
    8.734240e-07,
    9.981990e-07,
    9.981990e-07,
    1.122974e-06,
    1.122974e-06,
    1.372523e-06,
    1.497298e-06,
    1.622073e-06,
    1.871622e-06,
    1.871622e-06,
    1.996397e-06,
    2.370721e-06,
    2.620270e-06,
    2.620270e-06,
    2.620270e-06,
    2.745045e-06,
    2.745045e-06,
    2.745045e-06,
    2.745045e-06,
    2.869820e-06,
    3.493695e-06,
    3.493695e-06,
    3.493695e-06,
    3.868019e-06,
    4.117568e-06,
    4.367117e-06,
    4.866216e-06,
    4.866216e-06,
    5.240540e-06,
    5.240540e-06,
    5.864415e-06,
    6.363514e-06,
    6.737838e-06,
    7.236937e-06,
    7.860812e-06,
    8.110361e-06,
    8.484685e-06,
    9.233334e-06,
    9.857209e-06,
    1.035631e-05,
    1.060586e-05,
    1.085541e-05,
    1.210316e-05,
    1.272703e-05,
    1.322613e-05,
    1.409955e-05,
    1.509775e-05,
    1.647027e-05,
    1.784279e-05,
    1.859144e-05,
    1.921532e-05,
    2.008874e-05,
    2.133649e-05,
    2.245946e-05,
    2.333288e-05,
    2.445585e-05,
    2.557882e-05,
    2.732567e-05,
    2.832387e-05,
    3.019549e-05,
    3.231666e-05,
    3.468738e-05,
    3.568558e-05,
    3.718288e-05,
    3.942882e-05,
    4.292251e-05,
    4.479413e-05,
    4.666575e-05,
    4.903647e-05,
    5.140719e-05,
    5.352836e-05,
    5.539998e-05,
    5.714683e-05,
    5.964232e-05,
    6.263691e-05,
    6.550673e-05,
    6.825177e-05,
    7.124636e-05,
    7.411618e-05,
    7.611258e-05,
    7.910717e-05,
    8.235131e-05,
    8.534590e-05,
    8.971302e-05,
    9.457923e-05,
    9.907112e-05,
    1.028144e-04,
    1.068072e-04,
    1.124220e-04,
    1.150423e-04,
    1.202828e-04,
    1.241509e-04,
    1.297657e-04,
    1.352558e-04,
    1.408707e-04,
    1.474838e-04,
    1.535977e-04,
    1.605851e-04,
    1.679469e-04,
    1.766811e-04,
    1.849162e-04,
    1.943991e-04,
    2.023847e-04,
    2.111190e-04,
    2.208514e-04,
    2.298352e-04,
    2.390685e-04,
    2.500487e-04,
    2.599059e-04,
    2.705118e-04,
    2.826149e-04,
    2.953419e-04,
    3.058230e-04,
    3.209207e-04,
    3.336477e-04,
    3.469986e-04,
    3.607238e-04,
    3.758215e-04,
    3.915431e-04,
    4.100098e-04,
    4.283517e-04,
    4.454458e-04,
    4.661584e-04,
    4.845003e-04,
    5.032165e-04,
    5.194372e-04,
    5.399003e-04,
    5.642314e-04,
    5.905589e-04,
    6.173855e-04,
    6.407184e-04,
    6.646751e-04,
    6.896300e-04,
    7.199503e-04,
    7.487733e-04,
    7.780954e-04,
    8.119093e-04,
    8.442260e-04,
    8.782895e-04,
    9.142246e-04,
    9.507836e-04,
    9.890894e-04,
    1.025898e-03,
    1.060710e-03,
    1.102010e-03,
    1.141065e-03,
    1.179870e-03,
    1.225413e-03,
    1.268086e-03,
    1.313254e-03,
    1.362665e-03,
    1.413199e-03,
    1.465479e-03,
    1.514516e-03,
    1.568668e-03,
    1.620325e-03,
    1.678595e-03,
    1.738861e-03,
    1.800375e-03,
    1.863511e-03,
    1.933760e-03,
    1.999641e-03,
    2.066271e-03,
    2.136269e-03,
    2.208639e-03,
    2.282256e-03,
    2.363235e-03,
    2.451451e-03,
    2.535050e-03,
    2.626385e-03,
    2.715973e-03,
    2.811676e-03,
    2.906879e-03,
    3.008695e-03,
    3.108390e-03,
    3.220687e-03,
    3.337351e-03,
    3.458507e-03,
    3.587275e-03,
    3.717165e-03,
    3.853419e-03,
    3.998407e-03,
    4.150882e-03,
    4.307349e-03,
    4.470305e-03,
    4.650105e-03,
    4.835520e-03,
    5.028047e-03,
    5.239291e-03,
    5.459394e-03,
    5.699835e-03,
    5.944019e-03,
    6.201679e-03,
    6.472939e-03,
    6.767532e-03,
    7.067490e-03,
    7.402760e-03,
    7.779704e-03,
    8.163261e-03,
    8.581007e-03,
    9.043672e-03,
    9.538778e-03,
    1.006046e-02,
    1.061571e-02,
    1.121351e-02,
    1.187918e-02,
    1.259427e-02,
    1.336213e-02,
    1.420611e-02,
    1.510486e-02,
    1.609220e-02,
    1.716439e-02,
    1.831930e-02,
    1.961309e-02,
    2.100770e-02,
    2.254767e-02,
    2.423026e-02,
    2.604985e-02,
    2.803414e-02,
    3.023304e-02,
    3.261836e-02,
    3.527531e-02,
    3.818206e-02,
    4.134298e-02,
    4.483305e-02,
    4.870456e-02,
    5.296474e-02,
    5.765964e-02,
    6.281384e-02,
    6.850944e-02,
    7.481781e-02,
    8.173895e-02,
    8.937667e-02,
    9.774769e-02,
    1.069522e-01,
    1.170782e-01,
    1.282548e-01,
    1.405027e-01,
    1.539673e-01,
    1.686191e-01,
    1.846574e-01,
    2.020771e-01,
    2.210459e-01,
    2.415374e-01,
    2.637095e-01,
    2.872893e-01,
    3.124426e-01,
    3.391945e-01,
    3.673232e-01,
    3.966229e-01,
    4.269345e-01,
    4.580085e-01,
    4.895437e-01,
    5.211873e-01,
    5.525534e-01,
    5.832051e-01,
    6.129409e-01,
    6.414552e-01,
    6.685366e-01,
    6.941538e-01,
    7.180993e-01,
    7.403272e-01,
    7.608830e-01,
    7.798100e-01,
    7.972520e-01,
    8.132446e-01,
    8.277589e-01,
    8.410031e-01,
    8.531792e-01,
    8.642644e-01,
    8.743145e-01,
    8.835475e-01,
    8.919562e-01,
    8.996102e-01,
    9.066574e-01,
    9.130520e-01,
    9.189736e-01,
    9.244128e-01,
    9.293928e-01,
    9.339722e-01,
    9.381913e-01,
    9.421451e-01,
    9.457657e-01,
    9.491230e-01,
    9.522394e-01,
    9.551589e-01,
    9.578458e-01,
    9.603525e-01,
    9.626876e-01,
    9.648771e-01,
    9.669361e-01,
    9.688645e-01,
    9.706878e-01,
    9.723798e-01,
    9.739466e-01,
    9.754310e-01,
    9.768135e-01,
    9.781311e-01,
    9.793539e-01,
    9.805091e-01,
    9.816204e-01,
    9.826387e-01,
    9.835868e-01,
    9.844858e-01,
    9.853414e-01,
    9.861397e-01,
    9.869005e-01,
    9.876234e-01,
    9.882862e-01,
    9.889344e-01,
    9.895357e-01,
    9.901063e-01,
    9.906435e-01,
    9.911373e-01,
    9.916119e-01,
    9.920576e-01,
    9.924838e-01,
    9.928730e-01,
    9.932470e-01,
    9.935931e-01,
    9.939352e-01,
    9.942430e-01,
    9.945419e-01,
    9.948185e-01,
    9.950818e-01,
    9.953442e-01,
    9.955881e-01,
    9.958246e-01,
    9.960354e-01,
    9.962523e-01,
    9.964441e-01,
    9.966252e-01,
    9.967933e-01,
    9.969595e-01,
    9.971159e-01,
    9.972685e-01,
    9.974102e-01,
    9.975482e-01,
    9.976761e-01,
    9.977899e-01,
    9.978979e-01,
    9.979996e-01,
    9.980933e-01,
    9.981938e-01,
    9.982869e-01,
    9.983711e-01,
    9.984526e-01,
    9.985349e-01,
    9.986089e-01,
    9.986793e-01,
    9.987437e-01,
    9.988027e-01,
    9.988587e-01,
    9.989131e-01,
    9.989656e-01,
    9.990185e-01,
    9.990667e-01,
    9.991120e-01,
    9.991549e-01,
    9.991993e-01,
    9.992385e-01,
    9.992718e-01,
    9.993088e-01,
    9.993393e-01,
    9.993723e-01,
    9.993982e-01,
    9.994278e-01,
    9.994549e-01,
    9.994825e-01,
    9.995078e-01,
    9.995336e-01,
    9.995536e-01,
    9.995732e-01,
    9.995946e-01,
    9.996155e-01,
    9.996358e-01,
    9.996535e-01,
    9.996697e-01,
    9.996842e-01,
    9.996983e-01,
    9.997115e-01,
    9.997263e-01,
    9.997404e-01,
    9.997520e-01,
    9.997612e-01,
    9.997721e-01,
    9.997834e-01,
    9.997917e-01,
    9.998010e-01,
    9.998097e-01,
    9.998197e-01,
    9.998272e-01,
    9.998347e-01,
    9.998427e-01,
    9.998504e-01,
    9.998564e-01,
    9.998629e-01,
    9.998688e-01,
    9.998745e-01,
    9.998795e-01,
    9.998834e-01,
    9.998882e-01,
    9.998937e-01,
    9.998978e-01,
    9.999021e-01,
    9.999067e-01,
    9.999094e-01,
    9.999134e-01,
    9.999174e-01,
    9.999209e-01,
    9.999240e-01,
    9.999275e-01,
    9.999312e-01,
    9.999336e-01,
    9.999374e-01,
    9.999390e-01,
    9.999415e-01,
    9.999435e-01,
    9.999452e-01,
    9.999480e-01,
    9.999506e-01,
    9.999517e-01,
    9.999541e-01,
    9.999555e-01,
    9.999570e-01,
    9.999581e-01,
    9.999595e-01,
    9.999608e-01,
    9.999623e-01,
    9.999636e-01,
    9.999655e-01,
    9.999665e-01,
    9.999680e-01,
    9.999691e-01,
    9.999702e-01,
    9.999712e-01,
    9.999718e-01,
    9.999733e-01,
    9.999739e-01,
    9.999751e-01,
    9.999761e-01,
    9.999771e-01,
    9.999782e-01,
    9.999792e-01,
    9.999801e-01,
    9.999807e-01,
    9.999813e-01,
    9.999816e-01,
    9.999826e-01,
    9.999831e-01,
    9.999834e-01,
    9.999836e-01,
    9.999838e-01,
    9.999841e-01,
    9.999846e-01,
    9.999854e-01,
    9.999856e-01,
    9.999859e-01,
    9.999870e-01,
    9.999874e-01,
    9.999882e-01,
    9.999885e-01,
    9.999890e-01,
    9.999892e-01,
    9.999892e-01,
    9.999893e-01,
    9.999895e-01,
    9.999900e-01,
    9.999904e-01,
    9.999907e-01,
    9.999909e-01,
    9.999914e-01,
    9.999917e-01,
    9.999923e-01,
    9.999923e-01,
    9.999925e-01,
    9.999928e-01,
    9.999932e-01,
    9.999934e-01,
    9.999937e-01,
    9.999939e-01,
    9.999942e-01,
    9.999948e-01,
    9.999950e-01,
    9.999950e-01,
    9.999952e-01,
    9.999954e-01,
    9.999955e-01,
    9.999957e-01,
    9.999957e-01,
    9.999959e-01,
    9.999959e-01,
    9.999960e-01,
    9.999960e-01,
    9.999962e-01,
    9.999963e-01,
    9.999965e-01,
    9.999965e-01,
    9.999967e-01,
    9.999967e-01,
    9.999968e-01,
    9.999973e-01,
    9.999973e-01,
    9.999974e-01,
    9.999975e-01,
    9.999977e-01,
    9.999977e-01,
    9.999977e-01,
    9.999977e-01,
    9.999979e-01,
    9.999979e-01,
    9.999979e-01,
    9.999982e-01,
    9.999983e-01,
    9.999984e-01,
    9.999985e-01,
    9.999985e-01,
    9.999985e-01,
    9.999985e-01,
    9.999985e-01,
    9.999985e-01,
    9.999985e-01,
    9.999985e-01,
    9.999987e-01,
    9.999987e-01,
    9.999988e-01,
    9.999988e-01,
    9.999988e-01,
    9.999988e-01,
    9.999989e-01,
    9.999989e-01,
    9.999989e-01,
    9.999989e-01,
    9.999992e-01,
];

// ── Full composition matrix adjustment (Blast_CompositionMatrixAdj) ──────

/// NCBIstdaa alphabet to 20-letter ARND... conversion.
/// Port of NCBI alphaConvert[COMPO_LARGEST_ALPHABET].
pub static ALPHA_CONVERT: [i32; 28] = [
    -1, 0, -1, 4, 3, 6, 13, 7, 8, 9, 11, 10, 12, 2, 14, 5, 1, 15, 16, 19, 17, -1, 18, -1, -1, -1,
    -1, -1,
];

/// BLOSUM62 integer scores in ARND... 20-letter alphabet.
/// Port of NCBI `BLOS62[20][20]` from composition_adjustment.c.
pub static BLOS62: [[f64; COMPO_NUM_TRUE_AA]; COMPO_NUM_TRUE_AA] = [
    [
        4.0, -1.0, -2.0, -2.0, 0.0, -1.0, -1.0, 0.0, -2.0, -1.0, -1.0, -1.0, -1.0, -2.0, -1.0, 1.0,
        0.0, -3.0, -2.0, 0.0,
    ],
    [
        -1.0, 5.0, 0.0, -2.0, -3.0, 1.0, 0.0, -2.0, 0.0, -3.0, -2.0, 2.0, -1.0, -3.0, -2.0, -1.0,
        -1.0, -3.0, -2.0, -3.0,
    ],
    [
        -2.0, 0.0, 6.0, 1.0, -3.0, 0.0, 0.0, 0.0, 1.0, -3.0, -3.0, 0.0, -2.0, -3.0, -2.0, 1.0, 0.0,
        -4.0, -2.0, -3.0,
    ],
    [
        -2.0, -2.0, 1.0, 6.0, -3.0, 0.0, 2.0, -1.0, -1.0, -3.0, -4.0, -1.0, -3.0, -3.0, -1.0, 0.0,
        -1.0, -4.0, -3.0, -3.0,
    ],
    [
        0.0, -3.0, -3.0, -3.0, 9.0, -3.0, -4.0, -3.0, -3.0, -1.0, -1.0, -3.0, -1.0, -2.0, -3.0,
        -1.0, -1.0, -2.0, -2.0, -1.0,
    ],
    [
        -1.0, 1.0, 0.0, 0.0, -3.0, 5.0, 2.0, -2.0, 0.0, -3.0, -2.0, 1.0, 0.0, -3.0, -1.0, 0.0,
        -1.0, -2.0, -1.0, -2.0,
    ],
    [
        -1.0, 0.0, 0.0, 2.0, -4.0, 2.0, 5.0, -2.0, 0.0, -3.0, -3.0, 1.0, -2.0, -3.0, -1.0, 0.0,
        -1.0, -3.0, -2.0, -2.0,
    ],
    [
        0.0, -2.0, 0.0, -1.0, -3.0, -2.0, -2.0, 6.0, -2.0, -4.0, -4.0, -2.0, -3.0, -3.0, -2.0, 0.0,
        -2.0, -2.0, -3.0, -3.0,
    ],
    [
        -2.0, 0.0, 1.0, -1.0, -3.0, 0.0, 0.0, -2.0, 8.0, -3.0, -3.0, -1.0, -2.0, -1.0, -2.0, -1.0,
        -2.0, -2.0, 2.0, -3.0,
    ],
    [
        -1.0, -3.0, -3.0, -3.0, -1.0, -3.0, -3.0, -4.0, -3.0, 4.0, 2.0, -3.0, 1.0, 0.0, -3.0, -2.0,
        -1.0, -3.0, -1.0, 3.0,
    ],
    [
        -1.0, -2.0, -3.0, -4.0, -1.0, -2.0, -3.0, -4.0, -3.0, 2.0, 4.0, -2.0, 2.0, 0.0, -3.0, -2.0,
        -1.0, -2.0, -1.0, 1.0,
    ],
    [
        -1.0, 2.0, 0.0, -1.0, -3.0, 1.0, 1.0, -2.0, -1.0, -3.0, -2.0, 5.0, -1.0, -3.0, -1.0, 0.0,
        -1.0, -3.0, -2.0, -2.0,
    ],
    [
        -1.0, -1.0, -2.0, -3.0, -1.0, 0.0, -2.0, -3.0, -2.0, 1.0, 2.0, -1.0, 5.0, 0.0, -2.0, -1.0,
        -1.0, -1.0, -1.0, 1.0,
    ],
    [
        -2.0, -3.0, -3.0, -3.0, -2.0, -3.0, -3.0, -3.0, -1.0, 0.0, 0.0, -3.0, 0.0, 6.0, -4.0, -2.0,
        -2.0, 1.0, 3.0, -1.0,
    ],
    [
        -1.0, -2.0, -2.0, -1.0, -3.0, -1.0, -1.0, -2.0, -2.0, -3.0, -3.0, -1.0, -2.0, -4.0, 7.0,
        -1.0, -1.0, -4.0, -3.0, -2.0,
    ],
    [
        1.0, -1.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0, -2.0, -2.0, 0.0, -1.0, -2.0, -1.0, 4.0,
        1.0, -3.0, -2.0, -2.0,
    ],
    [
        0.0, -1.0, 0.0, -1.0, -1.0, -1.0, -1.0, -2.0, -2.0, -1.0, -1.0, -1.0, -1.0, -2.0, -1.0,
        1.0, 5.0, -2.0, -2.0, 0.0,
    ],
    [
        -3.0, -3.0, -4.0, -4.0, -2.0, -2.0, -3.0, -2.0, -2.0, -3.0, -2.0, -3.0, -1.0, 1.0, -4.0,
        -3.0, -2.0, 11.0, 2.0, -3.0,
    ],
    [
        -2.0, -2.0, -2.0, -3.0, -2.0, -1.0, -2.0, -3.0, 2.0, -1.0, -1.0, -2.0, -1.0, 3.0, -3.0,
        -2.0, -2.0, 2.0, 7.0, -1.0,
    ],
    [
        0.0, -3.0, -3.0, -3.0, -1.0, -2.0, -2.0, -3.0, -3.0, 3.0, 1.0, -2.0, 1.0, -1.0, -2.0, -2.0,
        0.0, -3.0, -1.0, 4.0,
    ],
];

/// BLOSUM62 joint probabilities in ARND... 20-letter alphabet.
/// Port of NCBI BLOSUM62_JOINT_PROBS from matrix_frequency_data.c.
pub static BLOSUM62_JOINT_PROBS: [[f64; COMPO_NUM_TRUE_AA]; COMPO_NUM_TRUE_AA] = [
    [
        2.1497573378347484e-02,
        2.3470224274721213e-03,
        1.9493235258876179e-03,
        2.1674844853066858e-03,
        1.5903351423026848e-03,
        1.9242657898716525e-03,
        2.9879059292799641e-03,
        5.8158526388051033e-03,
        1.1076584657559144e-03,
        3.1880644746334580e-03,
        4.4186245468471547e-03,
        3.3466571942021082e-03,
        1.3412107617355408e-03,
        1.6360627863999076e-03,
        2.1568959784943114e-03,
        6.2524987419815400e-03,
        3.7180506975672363e-03,
        4.0281679108936688e-04,
        1.2999956675626666e-03,
        5.0679056444508912e-03,
    ],
    [
        2.3470224274721213e-03,
        1.7757465118386322e-02,
        1.9786027128591904e-03,
        1.5865480081162602e-03,
        3.9365984789376245e-04,
        2.4858611089731411e-03,
        2.6933867548771758e-03,
        1.7221140903704937e-03,
        1.2407382229440791e-03,
        1.2435878276496955e-03,
        2.4193952633248727e-03,
        6.2339060289407083e-03,
        8.0309461712520876e-04,
        9.3181986323789834e-04,
        9.5783034332718700e-04,
        2.2660898636037261e-03,
        1.7802796534180537e-03,
        2.6571979312581875e-04,
        9.2634607111251918e-04,
        1.5810185245264004e-03,
    ],
    [
        1.9493235258876179e-03,
        1.9786027128591904e-03,
        1.4140291972553610e-02,
        3.7201973506001745e-03,
        4.3845466068066216e-04,
        1.5304436972610567e-03,
        2.2097156829738759e-03,
        2.8591871815612977e-03,
        1.4301072616183181e-03,
        9.9437221166923172e-04,
        1.3690958423974782e-03,
        2.4402105140841090e-03,
        5.2943633069226512e-04,
        7.5004227978192801e-04,
        8.6016459857770028e-04,
        3.1466019144814608e-03,
        2.2360795375444384e-03,
        1.6159545671597605e-04,
        7.0048422794024819e-04,
        1.2014015528772706e-03,
    ],
    [
        2.1674844853066858e-03,
        1.5865480081162602e-03,
        3.7201973506001745e-03,
        2.1274574617480089e-02,
        3.9909227141697264e-04,
        1.6481246723433428e-03,
        4.9158017471929655e-03,
        2.5221102126636373e-03,
        9.5384849402143984e-04,
        1.2347404942429857e-03,
        1.5202051791453383e-03,
        2.4453087721980561e-03,
        4.6429229320514104e-04,
        7.6023722413111566e-04,
        1.2373315413524663e-03,
        2.8035127901697272e-03,
        1.8961512776990257e-03,
        1.6218020183662784e-04,
        5.9842263937853702e-04,
        1.3158365660538270e-03,
    ],
    [
        1.5903351423026848e-03,
        3.9365984789376245e-04,
        4.3845466068066216e-04,
        3.9909227141697264e-04,
        1.1931352277704348e-02,
        3.0937204045913537e-04,
        3.8338775043186374e-04,
        7.6951976030099293e-04,
        2.2976387481074697e-04,
        1.0956590131781735e-03,
        1.5682982157153873e-03,
        5.0124929379033781e-04,
        3.7717165634097634e-04,
        5.1389991547056834e-04,
        3.6111795849154795e-04,
        1.0432626586831986e-03,
        9.3041313726939057e-04,
        1.4474923964368156e-04,
        3.4603772624580643e-04,
        1.3606607271146112e-03,
    ],
    [
        1.9242657898716525e-03,
        2.4858611089731411e-03,
        1.5304436972610567e-03,
        1.6481246723433428e-03,
        3.0937204045913537e-04,
        7.3292255467189687e-03,
        3.5385780499965817e-03,
        1.3683038039160171e-03,
        1.0489026828741754e-03,
        8.9102936026571569e-04,
        1.6174411456311808e-03,
        3.0968229715707327e-03,
        7.3993258722701268e-04,
        5.4255147972143906e-04,
        8.4668181752066874e-04,
        1.8931125300036275e-03,
        1.3796838284921874e-03,
        2.2737931366728891e-04,
        6.7584155312457842e-04,
        1.1660966117775285e-03,
    ],
    [
        2.9879059292799641e-03,
        2.6933867548771758e-03,
        2.2097156829738759e-03,
        4.9158017471929655e-03,
        3.8338775043186374e-04,
        3.5385780499965817e-03,
        1.6133927472163669e-02,
        1.9380952488713059e-03,
        1.3667885452189439e-03,
        1.2192061706431622e-03,
        2.0030316026648431e-03,
        4.1322603720305197e-03,
        6.7909745467514783e-04,
        8.5179405867513139e-04,
        1.4216207127018586e-03,
        2.9539180653600089e-03,
        2.0493063257644955e-03,
        2.6488552587183780e-04,
        8.7044186256788659e-04,
        1.6987763526262680e-03,
    ],
    [
        5.8158526388051033e-03,
        1.7221140903704937e-03,
        2.8591871815612977e-03,
        2.5221102126636373e-03,
        7.6951976030099293e-04,
        1.3683038039160171e-03,
        1.9380952488713059e-03,
        3.7804346453413303e-02,
        9.5813607255887238e-04,
        1.3849118546156933e-03,
        2.0864716056392773e-03,
        2.5392537741810947e-03,
        7.3281559749652399e-04,
        1.1976708695723554e-03,
        1.3641171883713547e-03,
        3.8342830901664762e-03,
        2.1858459940987062e-03,
        4.0740829083805248e-04,
        8.3467413018106177e-04,
        1.8218235950233687e-03,
    ],
    [
        1.1076584657559144e-03,
        1.2407382229440791e-03,
        1.4301072616183181e-03,
        9.5384849402143984e-04,
        2.2976387481074697e-04,
        1.0489026828741754e-03,
        1.3667885452189439e-03,
        9.5813607255887238e-04,
        9.2802502369336622e-03,
        5.8089627083019206e-04,
        9.8696608463236094e-04,
        1.1873625842258938e-03,
        3.8264639620910225e-04,
        8.1041076335565583e-04,
        4.7770135861914477e-04,
        1.1052034635193162e-03,
        7.4371746073077327e-04,
        1.5168037757411286e-04,
        1.5213771111755425e-03,
        6.4882907765797669e-04,
    ],
    [
        3.1880644746334580e-03,
        1.2435878276496955e-03,
        9.9437221166923172e-04,
        1.2347404942429857e-03,
        1.0956590131781735e-03,
        8.9102936026571569e-04,
        1.2192061706431622e-03,
        1.3849118546156933e-03,
        5.8089627083019206e-04,
        1.8441526588740136e-02,
        1.1382470627796603e-02,
        1.5655862274689192e-03,
        2.5081290988482057e-03,
        3.0458868657559346e-03,
        1.0068164685944146e-03,
        1.7225081689171561e-03,
        2.6953622613315018e-03,
        3.6183761166072852e-04,
        1.3821121844492116e-03,
        1.1972663837662637e-02,
    ],
    [
        4.4186245468471547e-03,
        2.4193952633248727e-03,
        1.3690958423974782e-03,
        1.5202051791453383e-03,
        1.5682982157153873e-03,
        1.6174411456311808e-03,
        2.0030316026648431e-03,
        2.0864716056392773e-03,
        9.8696608463236094e-04,
        1.1382470627796603e-02,
        3.7141460156350926e-02,
        2.4634345023228079e-03,
        4.9293545515183088e-03,
        5.4151301166464015e-03,
        1.4146090399381900e-03,
        2.4277107072013821e-03,
        3.3238031308707055e-03,
        7.3206640617832933e-04,
        2.2096734692836624e-03,
        9.4786263030457313e-03,
    ],
    [
        3.3466571942021082e-03,
        6.2339060289407083e-03,
        2.4402105140841090e-03,
        2.4453087721980561e-03,
        5.0124929379033781e-04,
        3.0968229715707327e-03,
        4.1322603720305197e-03,
        2.5392537741810947e-03,
        1.1873625842258938e-03,
        1.5655862274689192e-03,
        2.4634345023228079e-03,
        1.6113385590544604e-02,
        9.0876633395557617e-04,
        9.4875149773685364e-04,
        1.5773020912564391e-03,
        3.1016069999481111e-03,
        2.3467014804084987e-03,
        2.7198500003555514e-04,
        9.9908866586876396e-04,
        1.9360424083099779e-03,
    ],
    [
        1.3412107617355408e-03,
        8.0309461712520876e-04,
        5.2943633069226512e-04,
        4.6429229320514104e-04,
        3.7717165634097634e-04,
        7.3993258722701268e-04,
        6.7909745467514783e-04,
        7.3281559749652399e-04,
        3.8264639620910225e-04,
        2.5081290988482057e-03,
        4.9293545515183088e-03,
        9.0876633395557617e-04,
        4.0477309321969848e-03,
        1.1901770463553603e-03,
        4.0824445213456919e-04,
        8.5603787638552766e-04,
        1.0095451907679563e-03,
        1.9872537223131380e-04,
        5.7145288352831449e-04,
        2.3123361470140736e-03,
    ],
    [
        1.6360627863999076e-03,
        9.3181986323789834e-04,
        7.5004227978192801e-04,
        7.6023722413111566e-04,
        5.1389991547056834e-04,
        5.4255147972143906e-04,
        8.5179405867513139e-04,
        1.1976708695723554e-03,
        8.1041076335565583e-04,
        3.0458868657559346e-03,
        5.4151301166464015e-03,
        9.4875149773685364e-04,
        1.1901770463553603e-03,
        1.8277684015431908e-02,
        5.2528021756783813e-04,
        1.1939618185901600e-03,
        1.1624184369750680e-03,
        8.4917468952377874e-04,
        4.2392005745634370e-03,
        2.5763052227920180e-03,
    ],
    [
        2.1568959784943114e-03,
        9.5783034332718700e-04,
        8.6016459857770028e-04,
        1.2373315413524663e-03,
        3.6111795849154795e-04,
        8.4668181752066874e-04,
        1.4216207127018586e-03,
        1.3641171883713547e-03,
        4.7770135861914477e-04,
        1.0068164685944146e-03,
        1.4146090399381900e-03,
        1.5773020912564391e-03,
        4.0824445213456919e-04,
        5.2528021756783813e-04,
        1.9066033679132538e-02,
        1.6662567934883051e-03,
        1.3511005665728870e-03,
        1.4152209821874487e-04,
        4.5224391125285910e-04,
        1.2451325046931832e-03,
    ],
    [
        6.2524987419815400e-03,
        2.2660898636037261e-03,
        3.1466019144814608e-03,
        2.8035127901697272e-03,
        1.0432626586831986e-03,
        1.8931125300036275e-03,
        2.9539180653600089e-03,
        3.8342830901664762e-03,
        1.1052034635193162e-03,
        1.7225081689171561e-03,
        2.4277107072013821e-03,
        3.1016069999481111e-03,
        8.5603787638552766e-04,
        1.1939618185901600e-03,
        1.6662567934883051e-03,
        1.2585947097159817e-02,
        4.7004857686835334e-03,
        2.8731729176487776e-04,
        1.0299846310599138e-03,
        2.3587292053265561e-03,
    ],
    [
        3.7180506975672363e-03,
        1.7802796534180537e-03,
        2.2360795375444384e-03,
        1.8961512776990257e-03,
        9.3041313726939057e-04,
        1.3796838284921874e-03,
        2.0493063257644955e-03,
        2.1858459940987062e-03,
        7.4371746073077327e-04,
        2.6953622613315018e-03,
        3.3238031308707055e-03,
        2.3467014804084987e-03,
        1.0095451907679563e-03,
        1.1624184369750680e-03,
        1.3511005665728870e-03,
        4.7004857686835334e-03,
        1.2514818886617953e-02,
        2.8575770858467209e-04,
        9.4161039895612720e-04,
        3.6402328079338207e-03,
    ],
    [
        4.0281679108936688e-04,
        2.6571979312581875e-04,
        1.6159545671597605e-04,
        1.6218020183662784e-04,
        1.4474923964368156e-04,
        2.2737931366728891e-04,
        2.6488552587183780e-04,
        4.0740829083805248e-04,
        1.5168037757411286e-04,
        3.6183761166072852e-04,
        7.3206640617832933e-04,
        2.7198500003555514e-04,
        1.9872537223131380e-04,
        8.4917468952377874e-04,
        1.4152209821874487e-04,
        2.8731729176487776e-04,
        2.8575770858467209e-04,
        6.4699301717154852e-03,
        8.8744160259272527e-04,
        3.5578318710317554e-04,
    ],
    [
        1.2999956675626666e-03,
        9.2634607111251918e-04,
        7.0048422794024819e-04,
        5.9842263937853702e-04,
        3.4603772624580643e-04,
        6.7584155312457842e-04,
        8.7044186256788659e-04,
        8.3467413018106177e-04,
        1.5213771111755425e-03,
        1.3821121844492116e-03,
        2.2096734692836624e-03,
        9.9908866586876396e-04,
        5.7145288352831449e-04,
        4.2392005745634370e-03,
        4.5224391125285910e-04,
        1.0299846310599138e-03,
        9.4161039895612720e-04,
        8.8744160259272527e-04,
        1.0246100213822419e-02,
        1.5489827890922993e-03,
    ],
    [
        5.0679056444508912e-03,
        1.5810185245264004e-03,
        1.2014015528772706e-03,
        1.3158365660538270e-03,
        1.3606607271146112e-03,
        1.1660966117775285e-03,
        1.6987763526262680e-03,
        1.8218235950233687e-03,
        6.4882907765797669e-04,
        1.1972663837662637e-02,
        9.4786263030457313e-03,
        1.9360424083099779e-03,
        2.3123361470140736e-03,
        2.5763052227920180e-03,
        1.2451325046931832e-03,
        2.3587292053265561e-03,
        3.6402328079338207e-03,
        3.5578318710317554e-04,
        1.5489827890922993e-03,
        1.9631915140537640e-02,
    ],
];

/// BLOSUM62 background frequencies (20 true amino acids, ARND... order).
/// Port of NCBI BLOSUM62_bg from matrix_frequency_data.c.
pub static BLOSUM62_BG: [f64; COMPO_NUM_TRUE_AA] = [
    7.4216205067993410e-02,
    5.1614486141284638e-02,
    4.4645808512757915e-02,
    5.3626000838554413e-02,
    2.4687457167944848e-02,
    3.4259650591416023e-02,
    5.4311925684587502e-02,
    7.4146941452644999e-02,
    2.6212984805266227e-02,
    6.7917367618953756e-02,
    9.8907868497150955e-02,
    5.8155682303079680e-02,
    2.4990197579643110e-02,
    4.7418459742284751e-02,
    3.8538003320306206e-02,
    5.7229029476494421e-02,
    5.0891364550287033e-02,
    1.3029956129972148e-02,
    3.2281512313758580e-02,
    7.2919098205619245e-02,
];

// Constants from composition_adjustment.c
/// NCBI `kCompoAdjustErrTolerance` (`composition_adjustment.c:73`).
pub const COMPO_ADJUST_ERR_TOLERANCE: f64 = 0.000_000_01;
/// NCBI `kCompoAdjustIterationLimit` (`composition_adjustment.c:75`).
pub const COMPO_ADJUST_ITERATION_LIMIT: usize = 2000;
/// NCBI `kFixedReBlosum62` (`composition_adjustment.c:77`): target
/// relative entropy used by the rel-entropy-constrained matrix
/// adjustment path for BLOSUM62.
#[allow(dead_code)]
pub const FIXED_RE_BLOSUM62: f64 = 0.44;
/// NCBI `kReMatrixAdjustmentPseudocounts` (`redo_alignment.c:83`).
#[allow(dead_code)]
pub const RE_PSEUDOCOUNTS: i32 = 20;

// NCBIstdaa character indices — verbatim port of NCBI's anonymous
// enum in `composition_adjustment.h:46-49`.
#[allow(dead_code)]
pub const E_GAP_CHAR: usize = crate::encoding::NCBISTDAA_GAP as usize; // NCBI `eGapChar`
pub const E_BCHAR: usize = crate::encoding::NCBISTDAA_B as usize; // B = D or N
pub const E_CCHAR: usize = crate::encoding::NCBISTDAA_C as usize;
pub const E_DCHAR: usize = crate::encoding::NCBISTDAA_D as usize;
pub const E_ECHAR: usize = crate::encoding::NCBISTDAA_E as usize;
pub const E_ICHAR: usize = crate::encoding::NCBISTDAA_I as usize;
pub const E_LCHAR: usize = crate::encoding::NCBISTDAA_L as usize;
pub const E_NCHAR: usize = crate::encoding::NCBISTDAA_N as usize;
pub const E_QCHAR: usize = crate::encoding::NCBISTDAA_Q as usize;
pub const E_XCHAR: usize = crate::encoding::NCBISTDAA_X as usize; // X = unknown
pub const E_ZCHAR: usize = crate::encoding::NCBISTDAA_Z as usize; // Z = E or Q
pub const E_SELENOCYSTEINE: usize = crate::encoding::NCBISTDAA_U as usize; // U
pub const E_STOP_CHAR: usize = crate::encoding::NCBISTDAA_STOP as usize; // NCBI `eStopChar` ('*')
pub const E_OCHAR: usize = crate::encoding::NCBISTDAA_O as usize; // O = pyrrolysine
pub const E_JCHAR: usize = crate::encoding::NCBISTDAA_J as usize; // J = I or L
pub const MAXIMUM_X_SCORE: f64 = -1.0;

/// NCBI: s_GatherLetterProbs (composition_adjustment.c).
/// Convert NCBIstdaa probabilities to 20-letter ARND... alphabet.
pub fn s_gather_letter_probs(
    output: &mut [f64; COMPO_NUM_TRUE_AA],
    input: &[f64],
    alphsize: usize,
) {
    for c in 0..alphsize.min(28) {
        if ALPHA_CONVERT[c] >= 0 {
            output[ALPHA_CONVERT[c] as usize] = input[c];
        }
    }
}

/// NCBI: s_UnpackLetterProbs (composition_adjustment.c).
/// Convert 20-letter ARND... probs back to NCBIstdaa alphabet.
pub fn s_unpack_letter_probs(
    std_probs: &mut [f64],
    alphsize: usize,
    probs: &[f64; COMPO_NUM_TRUE_AA],
) {
    for c in 0..alphsize.min(28) {
        if ALPHA_CONVERT[c] >= 0 {
            std_probs[c] = probs[ALPHA_CONVERT[c] as usize];
        } else {
            std_probs[c] = 0.0;
        }
    }
}

/// NCBI: Blast_ApplyPseudocounts (composition_adjustment.c).
pub fn blast_apply_pseudocounts(
    probs20: &mut [f64; COMPO_NUM_TRUE_AA],
    number_of_observations: usize,
    background_probs20: &[f64; COMPO_NUM_TRUE_AA],
    pseudocounts: i32,
) {
    let mut sum = 0.0f64;
    for i in 0..COMPO_NUM_TRUE_AA {
        sum += probs20[i];
    }
    if sum == 0.0 {
        sum = 1.0;
    }
    let dpseudocounts = pseudocounts as f64;
    let weight = dpseudocounts / (number_of_observations as f64 + dpseudocounts);
    for i in 0..COMPO_NUM_TRUE_AA {
        probs20[i] = (1.0 - weight) * probs20[i] / sum + weight * background_probs20[i];
    }
}

/// NCBI: Blast_CalcFreqRatios (composition_adjustment.c).
/// Converts joint probabilities to frequency ratios in place.
pub fn blast_calc_freq_ratios(
    ratios: &mut [[f64; COMPO_NUM_TRUE_AA]; COMPO_NUM_TRUE_AA],
    alphsize: usize,
    row_prob: &[f64],
    col_prob: &[f64],
) {
    for i in 0..alphsize {
        if row_prob[i] > 0.0 {
            for j in 0..alphsize {
                if col_prob[j] > 0.0 {
                    ratios[i][j] /= row_prob[i] * col_prob[j];
                }
            }
        }
    }
}

/// NCBI: Blast_FreqRatioToScore (composition_adjustment.c).
pub fn blast_freq_ratio_to_score(
    matrix: &mut [[f64; COMPO_NUM_TRUE_AA]; COMPO_NUM_TRUE_AA],
    rows: usize,
    cols: usize,
    lambda: f64,
) {
    for i in 0..rows {
        for j in 0..cols {
            if matrix[i][j] == 0.0 {
                matrix[i][j] = COMPO_SCORE_MIN;
            } else {
                matrix[i][j] = matrix[i][j].ln() / lambda;
            }
        }
    }
}

/// NCBI: Blast_TargetFreqEntropy (composition_adjustment.c).
/// Compute relative entropy of target frequencies.
pub fn blast_target_freq_entropy(
    target_freq: &[[f64; COMPO_NUM_TRUE_AA]; COMPO_NUM_TRUE_AA],
) -> f64 {
    let mut row_prob = [0.0f64; COMPO_NUM_TRUE_AA];
    let mut col_prob = [0.0f64; COMPO_NUM_TRUE_AA];
    for i in 0..COMPO_NUM_TRUE_AA {
        for j in 0..COMPO_NUM_TRUE_AA {
            row_prob[i] += target_freq[i][j];
            col_prob[j] += target_freq[i][j];
        }
    }
    let mut entropy = 0.0f64;
    for i in 0..COMPO_NUM_TRUE_AA {
        for j in 0..COMPO_NUM_TRUE_AA {
            let freq = target_freq[i][j];
            entropy += freq * (freq / row_prob[i] / col_prob[j]).ln();
        }
    }
    entropy
}

/// NCBI: Blast_MatrixEntropy (composition_adjustment.c).
pub fn blast_matrix_entropy(
    matrix: &[[f64; COMPO_NUM_TRUE_AA]; COMPO_NUM_TRUE_AA],
    alphsize: usize,
    row_prob: &[f64; COMPO_NUM_TRUE_AA],
    col_prob: &[f64; COMPO_NUM_TRUE_AA],
    lambda: f64,
) -> f64 {
    let mut entropy = 0.0f64;
    for i in 0..alphsize {
        for j in 0..alphsize {
            let nat_score = lambda * matrix[i][j];
            entropy += nat_score * nat_score.exp() * row_prob[i] * col_prob[j];
        }
    }
    entropy
}

/// NCBI: Blast_EntropyOldFreqNewContext (composition_adjustment.c).
pub fn blast_entropy_old_freq_new_context(
    target_freq: &[[f64; COMPO_NUM_TRUE_AA]; COMPO_NUM_TRUE_AA],
    row_prob: &[f64; COMPO_NUM_TRUE_AA],
    col_prob: &[f64; COMPO_NUM_TRUE_AA],
) -> Option<(f64, f64)> {
    // Calculate old row/col probs from target freqs
    let mut old_row_prob = [0.0f64; COMPO_NUM_TRUE_AA];
    let mut old_col_prob = [0.0f64; COMPO_NUM_TRUE_AA];
    for i in 0..COMPO_NUM_TRUE_AA {
        for j in 0..COMPO_NUM_TRUE_AA {
            old_row_prob[i] += target_freq[i][j];
            old_col_prob[j] += target_freq[i][j];
        }
    }

    // Copy target_freq to scores, compute freq ratios, convert to scores
    let mut scores = *target_freq;
    blast_calc_freq_ratios(&mut scores, COMPO_NUM_TRUE_AA, &old_row_prob, &old_col_prob);
    blast_freq_ratio_to_score(&mut scores, COMPO_NUM_TRUE_AA, COMPO_NUM_TRUE_AA, 1.0);

    let (lambda, iter_count) =
        blast_calc_lambda_full_precision(&scores, row_prob, col_prob, COMPO_NUM_TRUE_AA);
    if iter_count >= LAMBDA_ITERATION_LIMIT {
        return None;
    }
    let entropy = blast_matrix_entropy(&scores, COMPO_NUM_TRUE_AA, row_prob, col_prob, lambda);
    Some((entropy, lambda))
}

/// NCBI: Blast_TrueAaToStdTargetFreqs (composition_adjustment.c).
/// Convert 20x20 target frequencies to NCBIstdaa alphabet matrix.
///
/// The `E_*CHAR < std_alphsize` bounds checks below are verbatim from
/// NCBI `Blast_TrueAaToStdTargetFreqs` where `Alphsize` is a runtime
/// parameter. In our port it's always `AA_SIZE = 28` and every
/// `E_*CHAR` is `< 28`, so the checks fold to true — `redundant_comparisons`
/// is suppressed locally rather than removing the NCBI-verbatim guards.
#[allow(clippy::redundant_comparisons)]
pub fn blast_true_aa_to_std_target_freqs(
    std_freq: &mut [[f64; AA_SIZE]; AA_SIZE],
    freq: &[[f64; COMPO_NUM_TRUE_AA]; COMPO_NUM_TRUE_AA],
) {
    let std_alphsize = AA_SIZE;
    // Compute normalization sum
    let mut sum = 0.0f64;
    for a in 0..COMPO_NUM_TRUE_AA {
        for b in 0..COMPO_NUM_TRUE_AA {
            sum += freq[a][b];
        }
    }

    for big_a in 0..std_alphsize {
        if big_a >= 28 || ALPHA_CONVERT[big_a] < 0 {
            for big_b in 0..std_alphsize {
                std_freq[big_a][big_b] = 0.0;
            }
        } else {
            let a = ALPHA_CONVERT[big_a] as usize;
            for big_b in 0..std_alphsize {
                if big_b >= 28 || ALPHA_CONVERT[big_b] < 0 {
                    std_freq[big_a][big_b] = 0.0;
                } else {
                    let b = ALPHA_CONVERT[big_b] as usize;
                    std_freq[big_a][big_b] = freq[a][b] / sum;
                }
            }
            // Set pairwise ambiguities
            if E_DCHAR < std_alphsize && E_NCHAR < std_alphsize && E_BCHAR < std_alphsize {
                std_freq[big_a][E_BCHAR] = std_freq[big_a][E_DCHAR] + std_freq[big_a][E_NCHAR];
            }
            if E_ECHAR < std_alphsize && E_QCHAR < std_alphsize && E_ZCHAR < std_alphsize {
                std_freq[big_a][E_ZCHAR] = std_freq[big_a][E_ECHAR] + std_freq[big_a][E_QCHAR];
            }
            if std_alphsize > E_JCHAR && E_ICHAR < std_alphsize && E_LCHAR < std_alphsize {
                std_freq[big_a][E_JCHAR] = std_freq[big_a][E_ICHAR] + std_freq[big_a][E_LCHAR];
            }
        }
    }
    // Add rows for ambiguity characters B, Z, J
    if E_BCHAR < std_alphsize && E_DCHAR < std_alphsize && E_NCHAR < std_alphsize {
        for j in 0..std_alphsize {
            std_freq[E_BCHAR][j] = std_freq[E_DCHAR][j] + std_freq[E_NCHAR][j];
        }
    }
    if E_ZCHAR < std_alphsize && E_ECHAR < std_alphsize && E_QCHAR < std_alphsize {
        for j in 0..std_alphsize {
            std_freq[E_ZCHAR][j] = std_freq[E_ECHAR][j] + std_freq[E_QCHAR][j];
        }
    }
    if std_alphsize > E_JCHAR && E_ICHAR < std_alphsize && E_LCHAR < std_alphsize {
        for j in 0..std_alphsize {
            std_freq[E_JCHAR][j] = std_freq[E_ICHAR][j] + std_freq[E_LCHAR][j];
        }
    }
}

/// NCBI: s_CalcAvgScore (composition_adjustment.c).
pub fn s_calc_avg_score(m: &[f64], alphsize: usize, inc_m: usize, probs: &[f64]) -> f64 {
    let mut score_ix = 0.0f64;
    for j in 0..alphsize {
        if j < 28 && ALPHA_CONVERT[j] >= 0 {
            score_ix += m[j * inc_m] * probs[j];
        }
    }
    score_ix
}

/// NCBI: s_SetXUOScores (composition_adjustment.c).
///
/// The `E_*CHAR < alphsize` guards are verbatim from NCBI; callers
/// always pass `alphsize = AA_SIZE = 28` so every compile-time-constant
/// E_*CHAR check folds to true. Suppress `redundant_comparisons`
/// locally rather than removing the defensive NCBI guards.
#[allow(clippy::redundant_comparisons)]
pub fn s_set_xuo_scores(
    m: &mut [[f64; AA_SIZE]; AA_SIZE],
    alphsize: usize,
    row_probs: &[f64],
    col_probs: &[f64],
) {
    let mut score_xx = 0.0f64;
    for i in 0..alphsize {
        if i < 28 && ALPHA_CONVERT[i] >= 0 {
            let avg_ix = s_calc_avg_score(&m[i], alphsize, 1, col_probs);
            if E_XCHAR < alphsize {
                m[i][E_XCHAR] = avg_ix.min(MAXIMUM_X_SCORE);
            }
            score_xx += avg_ix * row_probs[i];

            // Column X score: need to gather column i values
            let mut col_vals = vec![0.0f64; alphsize];
            for r in 0..alphsize {
                col_vals[r] = m[r][i];
            }
            if E_XCHAR < alphsize {
                m[E_XCHAR][i] =
                    s_calc_avg_score(&col_vals, alphsize, 1, row_probs).min(MAXIMUM_X_SCORE);
            }
        }
    }
    if E_XCHAR < alphsize {
        m[E_XCHAR][E_XCHAR] = score_xx.min(MAXIMUM_X_SCORE);
    }

    // Set X scores for pairwise ambiguities
    if E_BCHAR < alphsize && E_XCHAR < alphsize {
        m[E_BCHAR][E_XCHAR] =
            s_calc_avg_score(&m[E_BCHAR], alphsize, 1, col_probs).min(MAXIMUM_X_SCORE);
        let mut col_b = vec![0.0f64; alphsize];
        for r in 0..alphsize {
            col_b[r] = m[r][E_BCHAR];
        }
        m[E_XCHAR][E_BCHAR] = s_calc_avg_score(&col_b, alphsize, 1, row_probs).min(MAXIMUM_X_SCORE);
    }
    if E_ZCHAR < alphsize && E_XCHAR < alphsize {
        m[E_ZCHAR][E_XCHAR] =
            s_calc_avg_score(&m[E_ZCHAR], alphsize, 1, col_probs).min(MAXIMUM_X_SCORE);
        let mut col_z = vec![0.0f64; alphsize];
        for r in 0..alphsize {
            col_z[r] = m[r][E_ZCHAR];
        }
        m[E_XCHAR][E_ZCHAR] = s_calc_avg_score(&col_z, alphsize, 1, row_probs).min(MAXIMUM_X_SCORE);
    }
    if alphsize > E_JCHAR {
        m[E_JCHAR][E_XCHAR] =
            s_calc_avg_score(&m[E_JCHAR], alphsize, 1, col_probs).min(MAXIMUM_X_SCORE);
        let mut col_j = vec![0.0f64; alphsize];
        for r in 0..alphsize {
            col_j[r] = m[r][E_JCHAR];
        }
        m[E_XCHAR][E_JCHAR] = s_calc_avg_score(&col_j, alphsize, 1, row_probs).min(MAXIMUM_X_SCORE);
    }

    // Copy C scores to U (Selenocysteine)
    if E_SELENOCYSTEINE < alphsize && E_CCHAR < alphsize {
        m[E_SELENOCYSTEINE] = m[E_CCHAR];
        for i in 0..alphsize {
            m[i][E_SELENOCYSTEINE] = m[i][E_CCHAR];
        }
    }
    // Copy X scores to O (pyrrolysine)
    if alphsize > E_OCHAR && E_XCHAR < alphsize {
        m[E_OCHAR] = m[E_XCHAR];
        for i in 0..alphsize {
            m[i][E_OCHAR] = m[i][E_XCHAR];
        }
    }
}

/// NCBI: s_ScoresStdAlphabet (composition_adjustment.c).
/// Build integer scoring matrix in NCBIstdaa alphabet from optimized target freqs.
pub fn s_scores_std_alphabet(
    matrix: &mut [[i32; AA_SIZE]; AA_SIZE],
    alphsize: usize,
    target_freq: &[[f64; COMPO_NUM_TRUE_AA]; COMPO_NUM_TRUE_AA],
    start_matrix: &[[i32; AA_SIZE]; AA_SIZE],
    row_prob: &[f64; COMPO_NUM_TRUE_AA],
    col_prob: &[f64; COMPO_NUM_TRUE_AA],
    lambda: f64,
) {
    let mut row_prob_std = [0.0f64; AA_SIZE];
    let mut col_prob_std = [0.0f64; AA_SIZE];
    s_unpack_letter_probs(&mut row_prob_std, alphsize, row_prob);
    // Set pair ambiguity probs
    if E_BCHAR < alphsize {
        row_prob_std[E_BCHAR] = row_prob_std[E_DCHAR] + row_prob_std[E_NCHAR];
    }
    if E_ZCHAR < alphsize {
        row_prob_std[E_ZCHAR] = row_prob_std[E_ECHAR] + row_prob_std[E_QCHAR];
    }
    if alphsize > E_JCHAR {
        row_prob_std[E_JCHAR] = row_prob_std[E_ICHAR] + row_prob_std[E_LCHAR];
    }

    s_unpack_letter_probs(&mut col_prob_std, alphsize, col_prob);
    if E_BCHAR < alphsize {
        col_prob_std[E_BCHAR] = col_prob_std[E_DCHAR] + col_prob_std[E_NCHAR];
    }
    if E_ZCHAR < alphsize {
        col_prob_std[E_ZCHAR] = col_prob_std[E_ECHAR] + col_prob_std[E_QCHAR];
    }
    if alphsize > E_JCHAR {
        col_prob_std[E_JCHAR] = col_prob_std[E_ICHAR] + col_prob_std[E_LCHAR];
    }

    // Build double matrix: target freqs → freq ratios → scores
    let mut scores_f = [[0.0f64; AA_SIZE]; AA_SIZE];
    blast_true_aa_to_std_target_freqs(&mut scores_f, target_freq);

    // Freq ratios
    for i in 0..alphsize {
        if row_prob_std[i] > 0.0 {
            for j in 0..alphsize {
                if col_prob_std[j] > 0.0 {
                    scores_f[i][j] /= row_prob_std[i] * col_prob_std[j];
                }
            }
        }
    }
    // Freq ratio to score
    for i in 0..alphsize {
        for j in 0..alphsize {
            if scores_f[i][j] == 0.0 {
                scores_f[i][j] = COMPO_SCORE_MIN;
            } else {
                scores_f[i][j] = scores_f[i][j].ln() / lambda;
            }
        }
    }
    // Set X/U/O scores
    s_set_xuo_scores(&mut scores_f, alphsize, &row_prob_std, &col_prob_std);

    // Round to integer matrix
    for i in 0..alphsize {
        for j in 0..alphsize {
            if scores_f[i][j] < i32::MIN as f64 {
                matrix[i][j] = i32::MIN;
            } else {
                matrix[i][j] = crate::math::blast_nint(scores_f[i][j]) as i32;
            }
        }
    }
    // Preserve stop character scores from start matrix
    for i in 0..alphsize {
        matrix[i][E_STOP_CHAR] = start_matrix[i][E_STOP_CHAR];
        matrix[E_STOP_CHAR][i] = start_matrix[E_STOP_CHAR][i];
    }
}

/// NCBI: `Blast_CompositionMatrixAdj` (`composition_adjustment.c:1344`).
pub fn blast_composition_matrix_adj(
    matrix: &mut [[i32; AA_SIZE]; AA_SIZE],
    alphsize: usize,
    matrix_adjust_rule: crate::algo::blast::composition_adjustment::compo_mode_condition::MatrixAdjustRule,
    length1: i32,
    length2: i32,
    stdaa_row_probs: &[f64],
    stdaa_col_probs: &[f64],
    pseudocounts: i32,
    specified_re: f64,
    nr_record: &mut BlastCompositionWorkspace,
    matrix_info: &BlastMatrixInfo,
) -> i32 {
    use crate::algo::blast::composition_adjustment::compo_mode_condition::MatrixAdjustRule;
    use crate::algo::blast::composition_adjustment::optimize_target_freq::optimize_target_frequencies;

    let mut row_probs = [0.0f64; COMPO_NUM_TRUE_AA];
    let mut col_probs = [0.0f64; COMPO_NUM_TRUE_AA];
    s_gather_letter_probs(&mut row_probs, stdaa_row_probs, alphsize);
    s_gather_letter_probs(&mut col_probs, stdaa_col_probs, alphsize);

    let desired_re = match matrix_adjust_rule {
        MatrixAdjustRule::UnconstrainedRelEntropy => 0.0,
        MatrixAdjustRule::RelEntropyOldMatrixNewContext => {
            match blast_entropy_old_freq_new_context(&nr_record.mat_b, &row_probs, &col_probs) {
                Some((re, _lambda)) => re,
                None => 0.0,
            }
        }
        MatrixAdjustRule::RelEntropyOldMatrixOldContext => {
            blast_target_freq_entropy(&nr_record.mat_b)
        }
        MatrixAdjustRule::UserSpecifiedRelEntropy => specified_re,
        _ => return -1,
    };

    blast_apply_pseudocounts(
        &mut row_probs,
        length1.max(0) as usize,
        &nr_record.first_standard_freq,
        pseudocounts,
    );
    blast_apply_pseudocounts(
        &mut col_probs,
        length2.max(0) as usize,
        &nr_record.second_standard_freq,
        pseudocounts,
    );

    let n = COMPO_NUM_TRUE_AA * COMPO_NUM_TRUE_AA;
    let mut mat_b_flat = vec![0.0f64; n];
    for i in 0..COMPO_NUM_TRUE_AA {
        for j in 0..COMPO_NUM_TRUE_AA {
            mat_b_flat[i * COMPO_NUM_TRUE_AA + j] = nr_record.mat_b[i][j];
        }
    }
    let mut mat_final_flat = vec![0.0f64; n];

    let (status, _iterations) = optimize_target_frequencies(
        &mut mat_final_flat,
        COMPO_NUM_TRUE_AA,
        &mat_b_flat,
        &row_probs,
        &col_probs,
        desired_re > 0.0,
        desired_re,
        COMPO_ADJUST_ERR_TOLERANCE,
        COMPO_ADJUST_ITERATION_LIMIT,
    );

    if status != 0 {
        return status;
    }

    for i in 0..COMPO_NUM_TRUE_AA {
        for j in 0..COMPO_NUM_TRUE_AA {
            nr_record.mat_final[i][j] = mat_final_flat[i * COMPO_NUM_TRUE_AA + j];
        }
    }

    if matrix_info.start_matrix.len() < AA_SIZE {
        return -1;
    }
    let mut start_matrix = [[0i32; AA_SIZE]; AA_SIZE];
    for i in 0..AA_SIZE {
        if matrix_info.start_matrix[i].len() < AA_SIZE {
            return -1;
        }
        start_matrix[i].copy_from_slice(&matrix_info.start_matrix[i][..AA_SIZE]);
    }
    s_scores_std_alphabet(
        matrix,
        alphsize,
        &nr_record.mat_final,
        &start_matrix,
        &row_probs,
        &col_probs,
        matrix_info.ungapped_lambda,
    );

    0
}

/// NCBI: `Blast_CompositionBasedStats` (`composition_adjustment.c:1103`).
pub fn blast_composition_based_stats(
    matrix: &mut [Vec<i32>],
    lambda_ratio: &mut f64,
    ss: &BlastMatrixInfo,
    query_prob: &[f64],
    res_prob: &[f64],
    calc_lambda: fn(&[f64], i32, i32, f64) -> f64,
    p_value_adjustment: i32,
) -> i32 {
    let rows = ss.rows.max(0) as usize;
    let cols = ss.cols.max(0) as usize;
    if cols == 0 {
        return -1;
    }

    let (score_array, obs_min, obs_max) = if ss.position_based != 0 {
        if rows == 0 || ss.start_matrix.len() < rows {
            return -1;
        }
        let mut obs_min = 0i32;
        let mut obs_max = 0i32;
        for irow in 0..rows {
            if ss.start_matrix[irow].len() < cols {
                return -1;
            }
            for &jcol in &TRUE_CHAR_POSITIONS {
                if jcol >= cols {
                    continue;
                }
                let score = ss.start_matrix[irow][jcol];
                if score < obs_min && score > COMPO_SCORE_MIN as i32 {
                    obs_min = score;
                }
                if score > obs_max {
                    obs_max = score;
                }
            }
        }

        let range = obs_max - obs_min + 1;
        if range <= 0 {
            return -1;
        }
        let mut score_array = vec![0.0f64; range as usize];
        let one_pos_frac = 1.0 / rows as f64;
        for irow in 0..rows {
            for &jcol in &TRUE_CHAR_POSITIONS {
                if jcol >= cols {
                    continue;
                }
                let score = ss.start_matrix[irow][jcol];
                if score >= obs_min {
                    score_array[(score - obs_min) as usize] +=
                        one_pos_frac * res_prob.get(jcol).copied().unwrap_or(0.0);
                }
            }
        }
        (score_array, obs_min, obs_max)
    } else {
        if cols > AA_SIZE || ss.start_matrix.len() < cols {
            return -1;
        }
        let mut start_matrix = [[0i32; AA_SIZE]; AA_SIZE];
        for i in 0..cols {
            if ss.start_matrix[i].len() < cols {
                return -1;
            }
            start_matrix[i][..cols].copy_from_slice(&ss.start_matrix[i][..cols]);
        }
        s_get_matrix_score_probs(&start_matrix, cols, res_prob, query_prob)
    };
    let correct_ungapped_lambda = calc_lambda(&score_array, obs_min, obs_max, ss.ungapped_lambda);

    *lambda_ratio = correct_ungapped_lambda / ss.ungapped_lambda;
    if p_value_adjustment == 0 {
        *lambda_ratio = (*lambda_ratio).min(1.0);
    }
    *lambda_ratio = (*lambda_ratio).max(LAMBDA_RATIO_LOWER_BOUND);

    if *lambda_ratio > 0.0 {
        let scaled_lambda = ss.ungapped_lambda / *lambda_ratio;
        if ss.position_based != 0 {
            if matrix.len() < rows
                || ss.start_freq_ratios.len() < rows
                || ss.start_matrix.len() < rows
            {
                return -1;
            }
            for p in 0..rows {
                if matrix[p].len() < cols
                    || ss.start_freq_ratios[p].len() < cols
                    || ss.start_matrix[p].len() < cols
                {
                    return -1;
                }
                let mut row = [0.0f64; COMPO_LARGEST_ALPHABET];
                row[..cols].copy_from_slice(&ss.start_freq_ratios[p][..cols]);
                for value in row.iter_mut().take(cols) {
                    *value = if *value == 0.0 {
                        COMPO_SCORE_MIN
                    } else {
                        value.ln() / scaled_lambda
                    };
                }
                if E_XCHAR < cols {
                    let x_score = s_calc_avg_score(&row, cols, 1, res_prob).min(MAXIMUM_X_SCORE);
                    row[E_XCHAR] = x_score;
                    if E_SELENOCYSTEINE < cols && E_CCHAR < cols {
                        row[E_SELENOCYSTEINE] = row[E_CCHAR];
                    }
                    if cols > E_OCHAR {
                        row[E_OCHAR] = x_score;
                    }
                }
                for (j, value) in row.iter().enumerate().take(cols) {
                    matrix[p][j] = crate::math::blast_nint(*value) as i32;
                }
                if E_STOP_CHAR < cols {
                    matrix[p][E_STOP_CHAR] = ss.start_matrix[p][E_STOP_CHAR];
                }
            }
        } else {
            if cols > AA_SIZE || matrix.len() < cols || ss.start_freq_ratios.len() < cols {
                return -1;
            }
            let mut scores = [[0.0f64; AA_SIZE]; AA_SIZE];
            for i in 0..cols {
                if matrix[i].len() < cols
                    || ss.start_matrix[i].len() < cols
                    || ss.start_freq_ratios[i].len() < cols
                {
                    return -1;
                }
                scores[i][..cols].copy_from_slice(&ss.start_freq_ratios[i][..cols]);
            }

            for i in 0..cols {
                for j in 0..cols {
                    scores[i][j] = if scores[i][j] == 0.0 {
                        COMPO_SCORE_MIN
                    } else {
                        scores[i][j].ln() / scaled_lambda
                    };
                }
            }
            s_set_xuo_scores(&mut scores, cols, query_prob, res_prob);
            for i in 0..cols {
                for j in 0..cols {
                    matrix[i][j] = if scores[i][j] < i32::MIN as f64 {
                        i32::MIN
                    } else {
                        crate::math::blast_nint(scores[i][j]) as i32
                    };
                }
            }
            if E_STOP_CHAR < cols {
                for i in 0..cols {
                    matrix[i][E_STOP_CHAR] = ss.start_matrix[i][E_STOP_CHAR];
                    matrix[E_STOP_CHAR][i] = ss.start_matrix[E_STOP_CHAR][i];
                }
            }
        }
    }

    0
}

/// NCBI: `Blast_AdjustScores` (`composition_adjustment.c:1414`).
#[allow(clippy::too_many_arguments)]
pub fn blast_adjust_scores(
    matrix: &mut [Vec<i32>],
    query_composition: &BlastAminoAcidComposition,
    query_length: i32,
    subject_composition: &BlastAminoAcidComposition,
    subject_length: i32,
    matrix_info: &BlastMatrixInfo,
    composition_adjust_mode: ECompoAdjustModes,
    re_pseudocounts: i32,
    nr_record: &mut BlastCompositionWorkspace,
    matrix_adjust_rule: &mut crate::algo::blast::composition_adjustment::compo_mode_condition::MatrixAdjustRule,
    calc_lambda: fn(&[f64], i32, i32, f64) -> f64,
    pvalue_for_this_pair: &mut f64,
    composition_test_index: i32,
    ratio_to_pass_back: &mut f64,
) -> i32 {
    use crate::algo::blast::composition_adjustment::compo_mode_condition::{
        blast_choose_matrix_adjust_rule, MatrixAdjustRule,
    };

    let alphsize = matrix_info.cols.max(0) as usize;
    let mut lambda_for_pair;
    let iter_count;
    let mut permuted_query_probs = [0.0f64; COMPO_NUM_TRUE_AA];
    let mut permuted_match_probs = [0.0f64; COMPO_NUM_TRUE_AA];

    if query_composition.num_true_amino_acids == 0 || subject_composition.num_true_amino_acids == 0
    {
        return 1;
    }

    if composition_test_index > 0
        || (matrix_info.position_based == 0
            && composition_adjust_mode != ECompoAdjustModes::CompositionBasedStats)
    {
        s_gather_letter_probs(&mut permuted_query_probs, &query_composition.prob, alphsize);
        s_gather_letter_probs(
            &mut permuted_match_probs,
            &subject_composition.prob,
            alphsize,
        );
    }

    if composition_test_index > 0 {
        let (computed_lambda, computed_iter_count) = blast_calc_lambda_full_precision(
            &BLOS62,
            &permuted_query_probs,
            &permuted_match_probs,
            COMPO_NUM_TRUE_AA,
        );
        lambda_for_pair = computed_lambda;
        iter_count = computed_iter_count;
        if iter_count >= LAMBDA_ITERATION_LIMIT {
            lambda_for_pair = COMPO_MIN_LAMBDA;
        }
        *pvalue_for_this_pair = blast_composition_pvalue(lambda_for_pair);
    }

    if matrix_info.position_based != 0
        || composition_adjust_mode == ECompoAdjustModes::CompositionBasedStats
    {
        *matrix_adjust_rule = MatrixAdjustRule::ScaleOldMatrix;
    } else {
        *matrix_adjust_rule = blast_choose_matrix_adjust_rule(
            query_length.max(0) as usize,
            subject_length.max(0) as usize,
            &permuted_query_probs,
            &permuted_match_probs,
            composition_adjust_mode as u8,
        );
    }

    if *matrix_adjust_rule != MatrixAdjustRule::ScaleOldMatrix {
        if alphsize > AA_SIZE || matrix.len() < alphsize {
            return -1;
        }
        let mut square_matrix = [[0i32; AA_SIZE]; AA_SIZE];
        for i in 0..alphsize {
            if matrix[i].len() < alphsize {
                return -1;
            }
            square_matrix[i][..alphsize].copy_from_slice(&matrix[i][..alphsize]);
        }
        let status = blast_composition_matrix_adj(
            &mut square_matrix,
            alphsize,
            *matrix_adjust_rule,
            query_composition.num_true_amino_acids,
            subject_composition.num_true_amino_acids,
            &query_composition.prob,
            &subject_composition.prob,
            re_pseudocounts,
            FIXED_RE_BLOSUM62,
            nr_record,
            matrix_info,
        );
        *ratio_to_pass_back = 1.0;
        if status <= 0 {
            if status == 0 {
                for i in 0..alphsize {
                    matrix[i][..alphsize].copy_from_slice(&square_matrix[i][..alphsize]);
                }
            }
            return status;
        }
    }

    *matrix_adjust_rule = MatrixAdjustRule::ScaleOldMatrix;
    blast_composition_based_stats(
        matrix,
        ratio_to_pass_back,
        matrix_info,
        &query_composition.prob,
        &subject_composition.prob,
        calc_lambda,
        (composition_test_index > 0) as i32,
    )
}

/// NCBI: `Blast_CompositionWorkspaceFree` (`composition_adjustment.c:1253`).
pub fn blast_composition_workspace_free(pnr_record: &mut Option<Box<BlastCompositionWorkspace>>) {
    *pnr_record = None;
}

/// NCBI: `Blast_CompositionWorkspaceNew` (`composition_adjustment.c:1271`).
pub fn blast_composition_workspace_new() -> Option<Box<BlastCompositionWorkspace>> {
    Some(Box::new(BlastCompositionWorkspace {
        mat_b: [[0.0; COMPO_NUM_TRUE_AA]; COMPO_NUM_TRUE_AA],
        mat_final: [[0.0; COMPO_NUM_TRUE_AA]; COMPO_NUM_TRUE_AA],
        first_standard_freq: [0.0; COMPO_NUM_TRUE_AA],
        second_standard_freq: [0.0; COMPO_NUM_TRUE_AA],
    }))
}

/// NCBI: `Blast_CompositionWorkspaceInit` (`composition_adjustment.c:1317`).
pub fn blast_composition_workspace_init(
    nr_record: &mut BlastCompositionWorkspace,
    matrix_name: &str,
) -> i32 {
    let status = crate::algo::blast::composition_adjustment::matrix_frequency_data::blast_get_joint_probs_for_matrix(
        &mut nr_record.mat_b,
        &mut nr_record.first_standard_freq,
        &mut nr_record.second_standard_freq,
        matrix_name,
    );
    if status == 0 {
        0
    } else {
        eprintln!(
            "Matrix {} not currently supported for RE based adjustment",
            matrix_name
        );
        -1
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::encoding::{aminoacid_to_ncbistdaa_base, ncbistdaa_to_aminoacid_base};

    /// Pin the NCBIstdaa character indices to NCBI's enum values
    /// (`composition_adjustment.h:46-49`). Also cross-check them
    /// against the `encoding` tables so a rename there surfaces here.
    #[test]
    fn test_ncbistdaa_indices_match_ncbi_enum() {
        assert_eq!(E_GAP_CHAR, crate::encoding::NCBISTDAA_GAP as usize);
        assert_eq!(E_BCHAR, crate::encoding::NCBISTDAA_B as usize);
        assert_eq!(E_CCHAR, crate::encoding::NCBISTDAA_C as usize);
        assert_eq!(E_DCHAR, crate::encoding::NCBISTDAA_D as usize);
        assert_eq!(E_ECHAR, crate::encoding::NCBISTDAA_E as usize);
        assert_eq!(E_ICHAR, crate::encoding::NCBISTDAA_I as usize);
        assert_eq!(E_LCHAR, crate::encoding::NCBISTDAA_L as usize);
        assert_eq!(E_NCHAR, crate::encoding::NCBISTDAA_N as usize);
        assert_eq!(E_QCHAR, crate::encoding::NCBISTDAA_Q as usize);
        assert_eq!(E_XCHAR, crate::encoding::NCBISTDAA_X as usize);
        assert_eq!(E_ZCHAR, crate::encoding::NCBISTDAA_Z as usize);
        assert_eq!(E_SELENOCYSTEINE, crate::encoding::NCBISTDAA_U as usize);
        assert_eq!(E_STOP_CHAR, crate::encoding::NCBISTDAA_STOP as usize);
        assert_eq!(E_OCHAR, crate::encoding::NCBISTDAA_O as usize);
        assert_eq!(E_JCHAR, crate::encoding::NCBISTDAA_J as usize);

        // Cross-check with encoding helpers.
        assert_eq!(ncbistdaa_to_aminoacid_base(E_STOP_CHAR as u8), b'*');
        assert_eq!(ncbistdaa_to_aminoacid_base(E_OCHAR as u8), b'O');
        assert_eq!(ncbistdaa_to_aminoacid_base(E_GAP_CHAR as u8), b'-');
        assert_eq!(aminoacid_to_ncbistdaa_base(b'*') as usize, E_STOP_CHAR);
        assert_eq!(aminoacid_to_ncbistdaa_base(b'O') as usize, E_OCHAR);
        assert_eq!(aminoacid_to_ncbistdaa_base(b'J') as usize, E_JCHAR);
        assert_eq!(aminoacid_to_ncbistdaa_base(b'B') as usize, E_BCHAR);
        assert_eq!(aminoacid_to_ncbistdaa_base(b'Z') as usize, E_ZCHAR);
    }

    #[test]
    fn blast_karlin_lambda_nr_owns_nonnegative_expected_score_rejection() {
        let sprob = [0.0, 1.0];
        assert_eq!(blast_karlin_lambda_nr(&sprob, 0, 1, 0.5), -1.0);
    }

    #[test]
    fn blast_composition_pvalue_uses_ncbi_table_interpolation() {
        let lambda_at_low_cut = COMPO_MIN_LAMBDA + 35.0 * 0.001;
        assert_eq!(
            blast_composition_pvalue(lambda_at_low_cut - 0.000_5),
            P_LAMBDA_TABLE[35]
        );
        assert_eq!(
            blast_composition_pvalue(lambda_at_low_cut),
            P_LAMBDA_TABLE[35]
        );

        let interpolated = blast_composition_pvalue(lambda_at_low_cut + 0.000_25);
        let expected = 0.75 * P_LAMBDA_TABLE[35] + 0.25 * P_LAMBDA_TABLE[36];
        assert_eq!(interpolated, expected);
        assert_eq!(
            blast_composition_pvalue(COMPO_MIN_LAMBDA + 565.0 * 0.001),
            1.0
        );
    }

    #[test]
    fn blast_overall_pvalue_and_pe_transforms_match_ncbi_formulas() {
        let p_comp = 0.25;
        let p_alignment = 0.125;
        let product = p_comp * p_alignment;
        assert_eq!(
            blast_overall_p_value(p_comp, p_alignment),
            product * (1.0 - product.ln())
        );
        assert_eq!(blast_karlin_eto_p(0.0), 0.0);
        assert_eq!(blast_karlin_pto_e(-0.1), i32::MIN as f64);
        assert_eq!(blast_karlin_pto_e(1.1), i32::MIN as f64);
        assert_eq!(blast_karlin_pto_e(1.0), i32::MAX as f64);
        assert_eq!(blast_karlin_pto_e(0.25), -crate::math::blast_log1p(-0.25));
    }
}
