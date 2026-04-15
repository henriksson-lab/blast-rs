//! Verbatim port of NCBI compo_mode_condition.c.
//! Decides which matrix adjustment rule to use based on composition.

const COMPO_NUM_TRUE_AA: usize = 20;
const HIGH_PAIR_THRESHOLD: f64 = 0.4;
const LENGTH_LOWER_THRESHOLD: usize = 50;
const QUERY_MATCH_DISTANCE_THRESHOLD: f64 = 0.16;
const LENGTH_RATIO_THRESHOLD: f64 = 3.0;
const ANGLE_DEGREE_THRESHOLD: f64 = 70.0;
const PI: f64 = std::f64::consts::PI;

/// Port of NCBI EMatrixAdjustRule.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MatrixAdjustRule {
    DontAdjust = -1,
    ScaleOldMatrix = 0,
    UnconstrainedRelEntropy = 1,
    RelEntropyOldMatrixNewContext = 2,
    RelEntropyOldMatrixOldContext = 3,
    UserSpecifiedRelEntropy = 4,
}

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

/// Relative entropy (Jensen-Shannon divergence square root).
/// Port of NCBI Blast_GetRelativeEntropy.
pub fn relative_entropy(a: &[f64], b: &[f64]) -> f64 {
    let mut value = 0.0f64;
    for i in 0..COMPO_NUM_TRUE_AA {
        let temp = (a[i] + b[i]) / 2.0;
        if temp > 0.0 {
            if a[i] > 0.0 {
                value += a[i] * (a[i] / temp).ln() / 2.0;
            }
            if b[i] > 0.0 {
                value += b[i] * (b[i] / temp).ln() / 2.0;
            }
        }
    }
    if value < 0.0 {
        value = 0.0;
    }
    value.sqrt()
}

/// Port of NCBI s_HighPairFrequencies.
fn high_pair_frequencies(letter_probs: &[f64], length: usize) -> bool {
    if length <= LENGTH_LOWER_THRESHOLD {
        return false;
    }
    let mut max = 0.0f64;
    let mut second = 0.0f64;
    for i in 0..COMPO_NUM_TRUE_AA {
        if letter_probs[i] > second {
            second = letter_probs[i];
            if letter_probs[i] > max {
                second = max;
                max = letter_probs[i];
            }
        }
    }
    (max + second) > HIGH_PAIR_THRESHOLD
}

/// Port of NCBI s_HighPairEitherSeq.
fn high_pair_either_seq(p_query: &[f64], len1: usize, p_match: &[f64], len2: usize) -> bool {
    high_pair_frequencies(p_query, len1) || high_pair_frequencies(p_match, len2)
}

/// Port of NCBI s_TestToApplyREAdjustmentConditional.
/// Decides whether to use full matrix optimization or just scaling.
fn test_re_adjustment_conditional(
    len_query: usize,
    len_match: usize,
    p_query: &[f64; COMPO_NUM_TRUE_AA],
    p_match: &[f64; COMPO_NUM_TRUE_AA],
    p_matrix: &[f64; COMPO_NUM_TRUE_AA],
) -> MatrixAdjustRule {
    let mut corr_factor = 0.0f64;
    for i in 0..COMPO_NUM_TRUE_AA {
        corr_factor += (p_query[i] - p_matrix[i]) * (p_match[i] - p_matrix[i]);
    }

    let d_m_mat = relative_entropy(p_match, p_matrix);
    let d_q_mat = relative_entropy(p_query, p_matrix);
    let d_m_q = relative_entropy(p_match, p_query);

    let mut angle = ((d_m_mat * d_m_mat + d_q_mat * d_q_mat - d_m_q * d_m_q)
        / (2.0 * d_m_mat * d_q_mat))
        .acos();
    // Convert radians to degrees
    angle = angle * 180.0 / PI;

    let len_q = len_query as f64;
    let len_m = len_match as f64;
    let (len_large, len_small) = if len_q > len_m {
        (len_q, len_m)
    } else {
        (len_m, len_q)
    };

    if high_pair_either_seq(p_query, len_query, p_match, len_match) {
        MatrixAdjustRule::UserSpecifiedRelEntropy
    } else if d_m_q > QUERY_MATCH_DISTANCE_THRESHOLD
        && len_large / len_small > LENGTH_RATIO_THRESHOLD
        && angle > ANGLE_DEGREE_THRESHOLD
    {
        MatrixAdjustRule::ScaleOldMatrix
    } else {
        MatrixAdjustRule::UserSpecifiedRelEntropy
    }
}

/// Port of NCBI Blast_ChooseMatrixAdjustRule.
/// Returns the adjustment rule for a query-subject pair.
pub fn choose_matrix_adjust_rule(
    len_query: usize,
    len_match: usize,
    p_query: &[f64; COMPO_NUM_TRUE_AA],
    p_match: &[f64; COMPO_NUM_TRUE_AA],
    mode: u8, // ECompoAdjustModes
) -> MatrixAdjustRule {
    match mode {
        0 => MatrixAdjustRule::DontAdjust,     // eNoCompositionBasedStats
        1 => MatrixAdjustRule::ScaleOldMatrix, // eCompositionBasedStats
        2 => test_re_adjustment_conditional(
            // eCompositionMatrixAdjust
            len_query,
            len_match,
            p_query,
            p_match,
            &BLOSUM62_BG,
        ),
        3 => MatrixAdjustRule::UserSpecifiedRelEntropy, // eCompoForceFullMatrixAdjust
        _ => MatrixAdjustRule::DontAdjust,
    }
}

/// Map 28-element NCBIstdaa probability array to 20-element true AA array.
/// Port of NCBI s_GatherLetterProbs.
pub fn gather_letter_probs(prob28: &[f64], out20: &mut [f64; COMPO_NUM_TRUE_AA]) {
    use crate::composition::TRUE_CHAR_POSITIONS;
    for (k, &idx) in TRUE_CHAR_POSITIONS.iter().enumerate() {
        out20[k] = if idx < prob28.len() { prob28[idx] } else { 0.0 };
    }
}
