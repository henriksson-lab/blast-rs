//! Position-Specific Scoring Matrix (PSSM) for PSI-BLAST.
//! A PSSM represents the amino acid preferences at each position
//! of a multiple sequence alignment, used for iterative search.
//!
//! This module implements the NCBI PSI-BLAST PSSM computation algorithm:
//!   1. Henikoff position-based sequence weighting
//!   2. Weighted residue frequency (match weight) computation
//!   3. Pseudocount blending using matrix frequency ratios
//!   4. Conversion to integer log-odds scores
//!
//! References:
//!   - Altschul et al. (1997) Gapped BLAST and PSI-BLAST
//!   - Henikoff & Henikoff (1994) Position-based sequence weights
//!   - NCBI C source: blast_psi_priv.c

use crate::matrix::AA_SIZE;

/// Number of standard amino acid residue types used in pseudocount calculations.
const EFFECTIVE_ALPHABET: usize = 20;

/// Maximum number of independent observations for pseudocount calculation.
const MAX_IND_OBSERVATIONS: usize = 400;

/// Effective infinity for pseudocount weight.
const PSEUDO_MAX: f64 = 1_000_000.0;

/// Minimum value for relative entropy and other small quantities.
const POS_EPSILON: f64 = 0.0001;

/// Multiplier for entropy-based pseudocounts (PSEUDO_MULTIPLIER in posit.c).
const PSEUDO_MULT: f64 = 500.0;

/// Numerator of entropy-based pseudocount method.
const PSEUDO_NUMERATOR: f64 = 0.0457;

/// Exponent of denominator in entropy-based pseudocount method.
const PSEUDO_EXPONENT: f64 = 0.8;

/// Small initial pseudocount to avoid zero probabilities.
const PSEUDO_SMALL_INITIAL: f64 = 5.5;

/// Arbitrary constant for columns with zero observations.
const ZERO_OBS_PSEUDO: f64 = 30.0;

/// Scale factor for the scaled PSSM (used in _PSIConvertFreqRatiosToPSSM).
#[allow(dead_code)]
const PSI_SCALE_FACTOR: i32 = 200;

/// NCBIstdaa code for gap.
const GAP_RESIDUE: u8 = 0;

/// NCBIstdaa code for X (unknown).
const X_RESIDUE: u8 = 21;

/// Mapping from the 20-element effective alphabet index to NCBIstdaa residue
/// codes. This is the `charOrder` array from the NCBI C code
/// (s_fillColumnProbabilities in blast_psi_priv.c).
const CHAR_ORDER: [u8; EFFECTIVE_ALPHABET] = [
    1,  // A
    16, // R
    13, // N
    4,  // D
    3,  // C
    15, // Q
    5,  // E
    7,  // G
    8,  // H
    9,  // I
    11, // L
    10, // K
    12, // M
    6,  // F
    14, // P
    17, // S
    18, // T
    20, // W
    22, // Y
    19, // V
];

/// The 20 standard amino acid residue codes in NCBIstdaa encoding.
/// These are the residues that participate in scoring (indices 1..=20 and 22,
/// but excluding 21=X, so we list the actual standard 20).
#[cfg(test)]
const STD_RESIDUES: [u8; EFFECTIVE_ALPHABET] = [
    1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22,
];

/// BLOSUM62 underlying frequency ratios (28x28, NCBIstdaa order).
/// Copied from matrix_freq_ratios.c in the NCBI BLAST C core.
#[rustfmt::skip]
static BLOSUM62_FREQ_RATIOS: [[f64; AA_SIZE]; AA_SIZE] = [
 // row 0: gap
 [0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00],
 // row 1: A
 [0.00000000e+00, 3.90294070e+00, 5.64459671e-01, 8.67987664e-01,
  5.44605275e-01, 7.41264113e-01, 4.64893827e-01, 1.05686961e+00,
  5.69364849e-01, 6.32481035e-01, 7.75390239e-01, 6.01945975e-01,
  7.23150342e-01, 5.88307640e-01, 7.54121369e-01, 7.56803943e-01,
  6.12698600e-01, 1.47210399e+00, 9.84401956e-01, 9.36458396e-01,
  4.16548781e-01, 7.50000000e-01, 5.42611869e-01, 7.47274948e-01,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 6.14377313e-01],
 // row 2: B
 [0.00000000e+00, 5.64459671e-01, 4.43758048e+00, 3.45226274e-01,
  4.74290926e+00, 1.33503378e+00, 3.24101420e-01, 7.38524318e-01,
  9.25449581e-01, 3.33981361e-01, 8.54849426e-01, 2.97257620e-01,
  4.04640322e-01, 4.07083696e+00, 5.53838329e-01, 9.44103648e-01,
  7.02873767e-01, 1.05798620e+00, 8.26250098e-01, 3.51280513e-01,
  2.52855433e-01, 7.50000000e-01, 4.09444638e-01, 1.18382127e+00,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 3.12208474e-01],
 // row 3: C
 [0.00000000e+00, 8.67987664e-01, 3.45226274e-01, 1.95765857e+01,
  3.01454345e-01, 2.85934574e-01, 4.38990118e-01, 4.20387870e-01,
  3.55049505e-01, 6.53458801e-01, 3.49128465e-01, 6.42275633e-01,
  6.11354340e-01, 3.97802620e-01, 3.79562691e-01, 3.65781531e-01,
  3.08939296e-01, 7.38415701e-01, 7.40551692e-01, 7.55844055e-01,
  4.49983903e-01, 7.50000000e-01, 4.34203398e-01, 3.16819526e-01,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 6.46828489e-01],
 // row 4: D
 [0.00000000e+00, 5.44605275e-01, 4.74290926e+00, 3.01454345e-01,
  7.39792738e+00, 1.68781075e+00, 2.98969081e-01, 6.34301019e-01,
  6.78558839e-01, 3.39015407e-01, 7.84090406e-01, 2.86613046e-01,
  3.46454634e-01, 1.55385281e+00, 5.98716826e-01, 8.97081129e-01,
  5.73200024e-01, 9.13504624e-01, 6.94789868e-01, 3.36500142e-01,
  2.32102315e-01, 7.50000000e-01, 3.45683565e-01, 1.38195506e+00,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 3.07946931e-01],
 // row 5: E
 [0.00000000e+00, 7.41264113e-01, 1.33503378e+00, 2.85934574e-01,
  1.68781075e+00, 5.46952608e+00, 3.30743991e-01, 4.81267655e-01,
  9.60040718e-01, 3.30522558e-01, 1.30827885e+00, 3.72873704e-01,
  5.00342289e-01, 9.11298183e-01, 6.79202587e-01, 1.90173784e+00,
  9.60797602e-01, 9.50357185e-01, 7.41425610e-01, 4.28943130e-01,
  3.74300212e-01, 7.50000000e-01, 4.96467354e-01, 4.08949895e+00,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 3.55631838e-01],
 // row 6: F
 [0.00000000e+00, 4.64893827e-01, 3.24101420e-01, 4.38990118e-01,
  2.98969081e-01, 3.30743991e-01, 8.12879702e+00, 3.40640908e-01,
  6.51990521e-01, 9.45769883e-01, 3.44043119e-01, 1.15459749e+00,
  1.00437163e+00, 3.54288952e-01, 2.87444758e-01, 3.33972402e-01,
  3.80726330e-01, 4.39973597e-01, 4.81693683e-01, 7.45089738e-01,
  1.37437942e+00, 7.50000000e-01, 2.76938063e+00, 3.31992746e-01,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 1.06958025e+00],
 // row 7: G
 [0.00000000e+00, 1.05686961e+00, 7.38524318e-01, 4.20387870e-01,
  6.34301019e-01, 4.81267655e-01, 3.40640908e-01, 6.87630691e+00,
  4.92966576e-01, 2.75009722e-01, 5.88871736e-01, 2.84504012e-01,
  3.95486600e-01, 8.63711406e-01, 4.77385507e-01, 5.38649627e-01,
  4.49983999e-01, 9.03596525e-01, 5.79271582e-01, 3.36954912e-01,
  4.21690355e-01, 7.50000000e-01, 3.48714366e-01, 5.03463109e-01,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 2.80638726e-01],
 // row 8: H
 [0.00000000e+00, 5.69364849e-01, 9.25449581e-01, 3.55049505e-01,
  6.78558839e-01, 9.60040718e-01, 6.51990521e-01, 4.92966576e-01,
  1.35059997e+01, 3.26288125e-01, 7.78887490e-01, 3.80675486e-01,
  5.84132623e-01, 1.22200067e+00, 4.72879831e-01, 1.16798104e+00,
  9.17048021e-01, 7.36731740e-01, 5.57503254e-01, 3.39447442e-01,
  4.44088955e-01, 7.50000000e-01, 1.79790413e+00, 1.04047242e+00,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 3.58533474e-01],
 // row 9: I
 [0.00000000e+00, 6.32481035e-01, 3.33981361e-01, 6.53458801e-01,
  3.39015407e-01, 3.30522558e-01, 9.45769883e-01, 2.75009722e-01,
  3.26288125e-01, 3.99792994e+00, 3.96372934e-01, 1.69443475e+00,
  1.47774450e+00, 3.27934752e-01, 3.84662860e-01, 3.82937802e-01,
  3.54751311e-01, 4.43163582e-01, 7.79816110e-01, 2.41751209e+00,
  4.08874390e-01, 7.50000000e-01, 6.30388931e-01, 3.50796872e-01,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 2.63222650e+00],
 // row 10: K
 [0.00000000e+00, 7.75390239e-01, 8.54849426e-01, 3.49128465e-01,
  7.84090406e-01, 1.30827885e+00, 3.44043119e-01, 5.88871736e-01,
  7.78887490e-01, 3.96372934e-01, 4.76433717e+00, 4.28270363e-01,
  6.25302816e-01, 9.39841129e-01, 7.03774479e-01, 1.55432308e+00,
  2.07680867e+00, 9.31919141e-01, 7.92905803e-01, 4.56542720e-01,
  3.58930071e-01, 7.50000000e-01, 5.32179333e-01, 1.40344922e+00,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 4.15284382e-01],
 // row 11: L
 [0.00000000e+00, 6.01945975e-01, 2.97257620e-01, 6.42275633e-01,
  2.86613046e-01, 3.72873704e-01, 1.15459749e+00, 2.84504012e-01,
  3.80675486e-01, 1.69443475e+00, 4.28270363e-01, 3.79662137e+00,
  1.99429557e+00, 3.10043276e-01, 3.71121724e-01, 4.77325586e-01,
  4.73919278e-01, 4.28893743e-01, 6.60328975e-01, 1.31423573e+00,
  5.68037074e-01, 7.50000000e-01, 6.92059423e-01, 4.13275887e-01,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 2.94078574e+00],
 // row 12: M
 [0.00000000e+00, 7.23150342e-01, 4.04640322e-01, 6.11354340e-01,
  3.46454634e-01, 5.00342289e-01, 1.00437163e+00, 3.95486600e-01,
  5.84132623e-01, 1.47774450e+00, 6.25302816e-01, 1.99429557e+00,
  6.48145121e+00, 4.74529655e-01, 4.23898024e-01, 8.64250293e-01,
  6.22623369e-01, 5.98558924e-01, 7.93801616e-01, 1.26893679e+00,
  6.10296214e-01, 7.50000000e-01, 7.08364628e-01, 6.41102583e-01,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 1.78399892e+00],
 // row 13: N
 [0.00000000e+00, 5.88307640e-01, 4.07083696e+00, 3.97802620e-01,
  1.55385281e+00, 9.11298183e-01, 3.54288952e-01, 8.63711406e-01,
  1.22200067e+00, 3.27934752e-01, 9.39841129e-01, 3.10043276e-01,
  4.74529655e-01, 7.09409488e+00, 4.99932836e-01, 1.00058442e+00,
  8.58630478e-01, 1.23152924e+00, 9.84152635e-01, 3.69033853e-01,
  2.77782896e-01, 7.50000000e-01, 4.86030806e-01, 9.45834265e-01,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 3.17327197e-01],
 // row 14: P
 [0.00000000e+00, 7.54121369e-01, 5.53838329e-01, 3.79562691e-01,
  5.98716826e-01, 6.79202587e-01, 2.87444758e-01, 4.77385507e-01,
  4.72879831e-01, 3.84662860e-01, 7.03774479e-01, 3.71121724e-01,
  4.23898024e-01, 4.99932836e-01, 1.28375437e+01, 6.41280589e-01,
  4.81534905e-01, 7.55503259e-01, 6.88897122e-01, 4.43082984e-01,
  2.81833164e-01, 7.50000000e-01, 3.63521119e-01, 6.64534287e-01,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 3.76634549e-01],
 // row 15: Q
 [0.00000000e+00, 7.56803943e-01, 9.44103648e-01, 3.65781531e-01,
  8.97081129e-01, 1.90173784e+00, 3.33972402e-01, 5.38649627e-01,
  1.16798104e+00, 3.82937802e-01, 1.55432308e+00, 4.77325586e-01,
  8.64250293e-01, 1.00058442e+00, 6.41280589e-01, 6.24442175e+00,
  1.40579606e+00, 9.65555228e-01, 7.91320741e-01, 4.66777931e-01,
  5.09360272e-01, 7.50000000e-01, 6.11094097e-01, 3.58149606e+00,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 4.38898727e-01],
 // row 16: R
 [0.00000000e+00, 6.12698600e-01, 7.02873767e-01, 3.08939296e-01,
  5.73200024e-01, 9.60797602e-01, 3.80726330e-01, 4.49983999e-01,
  9.17048021e-01, 3.54751311e-01, 2.07680867e+00, 4.73919278e-01,
  6.22623369e-01, 8.58630478e-01, 4.81534905e-01, 1.40579606e+00,
  6.66557707e+00, 7.67165633e-01, 6.77754679e-01, 4.20072316e-01,
  3.95102106e-01, 7.50000000e-01, 5.55965425e-01, 1.13292384e+00,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 4.25403989e-01],
 // row 17: S
 [0.00000000e+00, 1.47210399e+00, 1.05798620e+00, 7.38415701e-01,
  9.13504624e-01, 9.50357185e-01, 4.39973597e-01, 9.03596525e-01,
  7.36731740e-01, 4.43163582e-01, 9.31919141e-01, 4.28893743e-01,
  5.98558924e-01, 1.23152924e+00, 7.55503259e-01, 9.65555228e-01,
  7.67165633e-01, 3.84284741e+00, 1.61392097e+00, 5.65223766e-01,
  3.85303035e-01, 7.50000000e-01, 5.57520051e-01, 9.56235816e-01,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 4.34703235e-01],
 // row 18: T
 [0.00000000e+00, 9.84401956e-01, 8.26250098e-01, 7.40551692e-01,
  6.94789868e-01, 7.41425610e-01, 4.81693683e-01, 5.79271582e-01,
  5.57503254e-01, 7.79816110e-01, 7.92905803e-01, 6.60328975e-01,
  7.93801616e-01, 9.84152635e-01, 6.88897122e-01, 7.91320741e-01,
  6.77754679e-01, 1.61392097e+00, 4.83210516e+00, 9.80943005e-01,
  4.30934144e-01, 7.50000000e-01, 5.73156574e-01, 7.60725140e-01,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 7.08974203e-01],
 // row 19: V
 [0.00000000e+00, 9.36458396e-01, 3.51280513e-01, 7.55844055e-01,
  3.36500142e-01, 4.28943130e-01, 7.45089738e-01, 3.36954912e-01,
  3.39447442e-01, 2.41751209e+00, 4.56542720e-01, 1.31423573e+00,
  1.26893679e+00, 3.69033853e-01, 4.43082984e-01, 4.66777931e-01,
  4.20072316e-01, 5.65223766e-01, 9.80943005e-01, 3.69215640e+00,
  3.74456332e-01, 7.50000000e-01, 6.58038693e-01, 4.43577702e-01,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 1.76339815e+00],
 // row 20: W
 [0.00000000e+00, 4.16548781e-01, 2.52855433e-01, 4.49983903e-01,
  2.32102315e-01, 3.74300212e-01, 1.37437942e+00, 4.21690355e-01,
  4.44088955e-01, 4.08874390e-01, 3.58930071e-01, 5.68037074e-01,
  6.10296214e-01, 2.77782896e-01, 2.81833164e-01, 5.09360272e-01,
  3.95102106e-01, 3.85303035e-01, 4.30934144e-01, 3.74456332e-01,
  3.81077833e+01, 7.50000000e-01, 2.10980812e+00, 4.26541694e-01,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 5.03239261e-01],
 // row 21: X
 [0.00000000e+00, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 7.50000000e-01],
 // row 22: Y
 [0.00000000e+00, 5.42611869e-01, 4.09444638e-01, 4.34203398e-01,
  3.45683565e-01, 4.96467354e-01, 2.76938063e+00, 3.48714366e-01,
  1.79790413e+00, 6.30388931e-01, 5.32179333e-01, 6.92059423e-01,
  7.08364628e-01, 4.86030806e-01, 3.63521119e-01, 6.11094097e-01,
  5.55965425e-01, 5.57520051e-01, 5.73156574e-01, 6.58038693e-01,
  2.10980812e+00, 7.50000000e-01, 9.83220341e+00, 5.40805192e-01,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 6.66952325e-01],
 // row 23: Z
 [0.00000000e+00, 7.47274948e-01, 1.18382127e+00, 3.16819526e-01,
  1.38195506e+00, 4.08949895e+00, 3.31992746e-01, 5.03463109e-01,
  1.04047242e+00, 3.50796872e-01, 1.40344922e+00, 4.13275887e-01,
  6.41102583e-01, 9.45834265e-01, 6.64534287e-01, 3.58149606e+00,
  1.13292384e+00, 9.56235816e-01, 7.60725140e-01, 4.43577702e-01,
  4.26541694e-01, 7.50000000e-01, 5.40805192e-01, 3.89300249e+00,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 3.87839626e-01],
 // row 24: U
 [0.00000000e+00, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 7.50000000e-01],
 // row 25: *
 [0.00000000e+00, 2.50000000e-01, 2.50000000e-01, 2.50000000e-01,
  2.50000000e-01, 2.50000000e-01, 2.50000000e-01, 2.50000000e-01,
  2.50000000e-01, 2.50000000e-01, 2.50000000e-01, 2.50000000e-01,
  2.50000000e-01, 2.50000000e-01, 2.50000000e-01, 2.50000000e-01,
  2.50000000e-01, 2.50000000e-01, 2.50000000e-01, 2.50000000e-01,
  2.50000000e-01, 2.50000000e-01, 2.50000000e-01, 2.50000000e-01,
  2.50000000e-01, 1.33300000e+00, 2.50000000e-01, 2.50000000e-01],
 // row 26: O
 [0.00000000e+00, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 7.50000000e-01],
 // row 27: J
 [0.00000000e+00, 6.14377313e-01, 3.12208474e-01, 6.46828489e-01,
  3.07946931e-01, 3.55631838e-01, 1.06958025e+00, 2.80638726e-01,
  3.58533474e-01, 2.63222650e+00, 4.15284382e-01, 2.94078574e+00,
  1.78399892e+00, 3.17327197e-01, 3.76634549e-01, 4.38898727e-01,
  4.25403989e-01, 4.34703235e-01, 7.08974203e-01, 1.76339815e+00,
  5.03239261e-01, 7.50000000e-01, 6.66952325e-01, 3.87839626e-01,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 2.81516607e+00],
];

/// Background frequencies for BLOSUM62, indexed by NCBIstdaa code (28 elements).
/// The standard 20 amino acid background probabilities from Robinson & Robinson
/// are placed at their NCBIstdaa positions; all other positions are 0.
///
/// The AA_FREQUENCIES in matrix.rs are in alphabetical order (A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y),
/// corresponding to NCBIstdaa codes [1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22].
fn std_prob_ncbistdaa() -> [f64; AA_SIZE] {
    use crate::matrix::AA_FREQUENCIES;
    let mut prob = [0.0f64; AA_SIZE];
    // Map from AA_FREQUENCIES (alphabetical) to NCBIstdaa codes
    let aa_to_stdaa: [usize; 20] = [1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22];
    for (i, &code) in aa_to_stdaa.iter().enumerate() {
        prob[code] = AA_FREQUENCIES[i];
    }
    prob
}

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
    ///
    /// This implements the full NCBI PSI-BLAST PSSM computation:
    ///   1. Henikoff position-based sequence weighting
    ///   2. Weighted residue frequency computation
    ///   3. Pseudocount blending via BLOSUM62 frequency ratios
    ///   4. Conversion to integer log-odds scores
    ///
    /// # Arguments
    /// * `aligned_seqs` - Aligned subject sequences in NCBIstdaa encoding.
    ///   Each sequence should be the same length as the query (self.length).
    ///   Gaps are represented as 0, unknown as 21.
    /// * `bg_freqs` - Background amino acid frequencies (20 values in alphabetical
    ///   order: A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y). Used for reference
    ///   but the actual NCBIstdaa-indexed probabilities are computed internally.
    /// * `_pseudocount_weight` - Legacy parameter (ignored; column-specific
    ///   pseudocounts are computed automatically per NCBI algorithm).
    pub fn update_from_alignment(
        &mut self,
        aligned_seqs: &[Vec<u8>],
        _bg_freqs: &[f64; 20],
        _pseudocount_weight: f64,
    ) {
        if aligned_seqs.is_empty() {
            return;
        }

        let std_prob = std_prob_ncbistdaa();
        let lambda = 0.3176; // BLOSUM62 ungapped lambda

        // Step 1: Compute Henikoff position-based sequence weights
        let match_weights = compute_sequence_weights_and_match_weights(
            aligned_seqs,
            self.length,
        );

        // Step 2: Compute effective observations and column-specific pseudocounts
        // Step 3: Blend observed frequencies with pseudocounts using frequency ratios
        // Step 4: Convert to integer PSSM scores

        for pos in 0..self.length {
            // Check if this position has any aligned residues
            let total_weight: f64 = (0..AA_SIZE)
                .map(|r| match_weights[pos][r])
                .sum();

            if total_weight < POS_EPSILON {
                // No alignment data at this position, keep original scores
                continue;
            }

            // Compute effective number of observations
            let observations = compute_effective_observations(
                aligned_seqs,
                pos,
                self.length,
                &std_prob,
            );

            // Compute column-specific pseudocounts
            let pseudo_weight = compute_column_pseudocounts(
                &match_weights[pos],
                &std_prob,
                observations,
            );

            // Compute frequency ratios with pseudocount blending
            for r in 0..AA_SIZE {
                if std_prob[r] <= POS_EPSILON {
                    // Score for non-standard residues: keep as-is or set to minimum
                    continue;
                }

                // pseudo = beta * sum_i(match_weights[pos][i] * freq_ratio[r][i])
                let mut pseudo = 0.0;
                for i in 0..AA_SIZE {
                    pseudo += match_weights[pos][i] * BLOSUM62_FREQ_RATIOS[r][i];
                }
                pseudo *= pseudo_weight;

                let numerator =
                    observations * match_weights[pos][r] / std_prob[r] + pseudo;
                let denominator = observations + pseudo_weight;

                if denominator <= 0.0 {
                    continue;
                }

                let q_over_p = numerator / denominator;

                // freq_ratio = q_over_p * std_prob[r] (the C code multiplies then
                // divides by std_prob in the next stage -- we go directly to score)
                // Score = round((1/lambda) * ln(q_over_p))
                if q_over_p > POS_EPSILON {
                    let score = (q_over_p.ln() / lambda).round() as i32;
                    self.scores[pos][r] = score;
                } else {
                    self.scores[pos][r] = i32::MIN / 2; // BLAST_SCORE_MIN equivalent
                }
            }
        }

        // Compute information content
        self.info_content = compute_information_content(&self.scores, &std_prob, self.length, lambda);
    }
}

/// Compute Henikoff position-based sequence weights and match weights.
///
/// For each column, count distinct residues (ignoring gaps and X),
/// then weight each sequence's contribution inversely proportional to
/// the number of distinct residues and count of that residue type.
/// Finally, normalize weights per column so they sum to 1, and
/// accumulate into match_weights[pos][residue].
fn compute_sequence_weights_and_match_weights(
    aligned_seqs: &[Vec<u8>],
    query_length: usize,
) -> Vec<[f64; AA_SIZE]> {
    let num_seqs = aligned_seqs.len();
    let mut match_weights = vec![[0.0f64; AA_SIZE]; query_length];

    for pos in 0..query_length {
        // Count occurrences of each residue at this position
        let mut residue_counts = [0u32; AA_SIZE];
        let mut num_aligned = 0u32;

        for seq in aligned_seqs.iter() {
            if pos < seq.len() {
                let r = seq[pos] as usize;
                if r < AA_SIZE {
                    residue_counts[r] += 1;
                    num_aligned += 1;
                }
            }
        }

        if num_aligned == 0 {
            continue;
        }

        // Count distinct residues (excluding gap=0 and X=21)
        let mut num_distinct = 0u32;
        for r in 0..AA_SIZE {
            if r == GAP_RESIDUE as usize || r == X_RESIDUE as usize {
                continue;
            }
            if residue_counts[r] > 0 {
                num_distinct += 1;
            }
        }

        if num_distinct == 0 {
            // All gaps or X -- distribute weight uniformly
            continue;
        }

        // Compute per-sequence weight contributions for this column
        // using simplified Henikoff weighting (without alignment extents)
        let mut seq_weights = vec![0.0f64; num_seqs];
        let mut weight_sum = 0.0;

        for (si, seq) in aligned_seqs.iter().enumerate() {
            if pos < seq.len() {
                let r = seq[pos] as usize;
                if r < AA_SIZE && r != GAP_RESIDUE as usize && r != X_RESIDUE as usize {
                    let w = 1.0 / (num_distinct as f64 * residue_counts[r] as f64);
                    seq_weights[si] = w;
                    weight_sum += w;
                }
            }
        }

        if weight_sum <= 0.0 {
            continue;
        }

        // Normalize and accumulate match weights
        for (si, seq) in aligned_seqs.iter().enumerate() {
            if pos < seq.len() && seq_weights[si] > 0.0 {
                let r = seq[pos] as usize;
                if r < AA_SIZE {
                    let norm_weight = seq_weights[si] / weight_sum;
                    match_weights[pos][r] += norm_weight;
                }
            }
        }
    }

    match_weights
}

/// Compute the effective number of independent observations at a column.
///
/// This is a simplified version of the NCBI s_effectiveObservations function.
/// It estimates how many independent sequences contribute to a column based
/// on the average number of distinct amino acid types seen.
fn compute_effective_observations(
    aligned_seqs: &[Vec<u8>],
    pos: usize,
    _query_length: usize,
    bg_prob: &[f64; AA_SIZE],
) -> f64 {
    // Count distinct residues at this position (excluding gaps and X)
    let mut residue_seen = [false; AA_SIZE];
    let mut num_participating = 0u32;

    for seq in aligned_seqs.iter() {
        if pos < seq.len() {
            let r = seq[pos] as usize;
            if r < AA_SIZE && r != GAP_RESIDUE as usize && r != X_RESIDUE as usize {
                residue_seen[r] = true;
                num_participating += 1;
            }
        }
    }

    let num_distinct: usize = residue_seen.iter().filter(|&&v| v).count();

    if num_distinct == 0 || num_participating == 0 {
        return 0.0;
    }

    // Build the expno table: expected number of distinct amino acids
    // as a function of number of independent trials
    // expno[j] = EFFECTIVE_ALPHABET - sum_k(exp(j * ln(1 - bg_prob[k])))
    // where bg_prob is over the 20 standard amino acids
    let bg20 = effective_bg_probs(bg_prob);

    // Find the number of independent observations that would give
    // the observed number of distinct residues
    let ave_distinct = num_distinct as f64;
    let mut indep = 0.0f64;

    // Binary/linear search through expected values
    for j in 1..=MAX_IND_OBSERVATIONS {
        let mut weighted_sum = 0.0;
        for k in 0..EFFECTIVE_ALPHABET {
            weighted_sum += (j as f64 * (1.0 - bg20[k]).ln()).exp();
        }
        let expected = EFFECTIVE_ALPHABET as f64 - weighted_sum;

        if expected >= ave_distinct {
            // Interpolate
            if j == 1 {
                indep = j as f64;
            } else {
                let mut prev_weighted_sum = 0.0;
                for k in 0..EFFECTIVE_ALPHABET {
                    prev_weighted_sum += ((j - 1) as f64 * (1.0 - bg20[k]).ln()).exp();
                }
                let prev_expected = EFFECTIVE_ALPHABET as f64 - prev_weighted_sum;
                indep = j as f64
                    - (expected - ave_distinct) / (expected - prev_expected);
            }
            break;
        }
        if j == MAX_IND_OBSERVATIONS {
            indep = MAX_IND_OBSERVATIONS as f64;
        }
    }

    // Cap at number of participating sequences, then subtract 1
    indep = indep.min(num_participating as f64);
    indep = (indep - 1.0).max(0.0);
    indep
}

/// Extract the 20 standard background probabilities in charOrder for
/// the effective alphabet computations.
fn effective_bg_probs(std_prob: &[f64; AA_SIZE]) -> [f64; EFFECTIVE_ALPHABET] {
    let mut bg = [0.0f64; EFFECTIVE_ALPHABET];
    for (i, &code) in CHAR_ORDER.iter().enumerate() {
        bg[i] = std_prob[code as usize];
    }
    // Normalize to sum to 1.0
    let sum: f64 = bg.iter().sum();
    if sum > 0.0 {
        for v in bg.iter_mut() {
            *v /= sum;
        }
    }
    bg
}

/// Compute column-specific pseudocounts following the NCBI entropy-based method.
///
/// Matches s_columnSpecificPseudocounts in blast_psi_priv.c:
///   1. Fill column probabilities from match weights
///   2. Adjust with small initial pseudocounts to avoid zeros
///   3. Compute relative entropy
///   4. pseudoWeight = PSEUDO_MULT * alpha / (1 - alpha)
///      where alpha = PSEUDO_NUMERATOR / (relativeEntropy ^ PSEUDO_EXPONENT)
fn compute_column_pseudocounts(
    match_weights_at_pos: &[f64; AA_SIZE],
    std_prob: &[f64; AA_SIZE],
    observations: f64,
) -> f64 {
    if observations < POS_EPSILON {
        return ZERO_OBS_PSEUDO;
    }

    // Step 1: Fill column probabilities (in charOrder)
    let mut col_probs = [0.0f64; EFFECTIVE_ALPHABET];
    for (c, &code) in CHAR_ORDER.iter().enumerate() {
        col_probs[c] = match_weights_at_pos[code as usize];
    }

    // Step 2: Get background probabilities in charOrder
    let bg = effective_bg_probs(std_prob);

    // Step 3: Adjust column probabilities with small initial pseudocounts
    // intermediateSums[c] = col_probs[c] * observations + bg[c] * PSEUDO_SMALL_INITIAL
    // then normalize
    let mut adjusted = [0.0f64; EFFECTIVE_ALPHABET];
    let mut overall_sum = 0.0;
    for c in 0..EFFECTIVE_ALPHABET {
        let v = col_probs[c] * observations + bg[c] * PSEUDO_SMALL_INITIAL;
        adjusted[c] = v;
        overall_sum += v;
    }
    if overall_sum > 0.0 {
        for c in 0..EFFECTIVE_ALPHABET {
            adjusted[c] /= overall_sum;
        }
    }

    // Step 4: Compute relative entropy
    let mut rel_entropy = 0.0;
    for c in 0..EFFECTIVE_ALPHABET {
        if adjusted[c] > POS_EPSILON {
            rel_entropy += adjusted[c] * (adjusted[c] / bg[c]).ln();
        }
    }
    if rel_entropy < POS_EPSILON {
        rel_entropy = POS_EPSILON;
    }

    // Step 5: Compute pseudocount weight
    let pseudo_denom = rel_entropy.powf(PSEUDO_EXPONENT);
    let alpha = PSEUDO_NUMERATOR / pseudo_denom;

    if alpha < 1.0 - POS_EPSILON {
        PSEUDO_MULT * alpha / (1.0 - alpha)
    } else {
        PSEUDO_MAX
    }
}

/// Compute information content per position from PSSM scores.
fn compute_information_content(
    scores: &[[i32; AA_SIZE]],
    std_prob: &[f64; AA_SIZE],
    query_length: usize,
    lambda: f64,
) -> Vec<f64> {
    let mut info = vec![0.0f64; query_length];

    for pos in 0..query_length {
        let mut sum = 0.0;
        for r in 0..AA_SIZE {
            if std_prob[r] > POS_EPSILON {
                let score = scores[pos][r] as f64;
                let freq_ratio = (lambda * score).exp();
                let weighted_freq = freq_ratio * std_prob[r];
                if weighted_freq > POS_EPSILON {
                    sum += weighted_freq * weighted_freq.ln();
                }
            }
        }
        info[pos] = sum / std::f64::consts::LN_2;
    }

    info
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

    #[test]
    fn test_henikoff_weights_identical_seqs() {
        // All identical sequences should produce uniform weights
        let seqs = vec![
            vec![1u8, 4, 7, 10],  // A, D, G, K
            vec![1u8, 4, 7, 10],
            vec![1u8, 4, 7, 10],
        ];
        let mw = compute_sequence_weights_and_match_weights(&seqs, 4);
        // Each position should have weight 1.0 for the single residue present
        for pos in 0..4 {
            let total: f64 = mw[pos].iter().sum();
            assert!((total - 1.0).abs() < 0.01,
                    "Match weights should sum to ~1 at pos {}, got {}", pos, total);
        }
    }

    #[test]
    fn test_henikoff_weights_diverse_seqs() {
        // Diverse sequences: rare residues should get more weight
        let seqs = vec![
            vec![1u8, 4],  // A, D
            vec![1u8, 4],  // A, D
            vec![1u8, 5],  // A, E
        ];
        let mw = compute_sequence_weights_and_match_weights(&seqs, 2);
        // At position 0: all have A, so weight should go to A=1
        assert!(mw[0][1] > 0.9, "All-A column should have weight ~1 on A");
        // At position 1: D appears 2x, E appears 1x
        // Henikoff: D seqs get 1/(2*2) = 0.25 each, E seq gets 1/(2*1) = 0.5
        // normalized: D total = 0.5, E total = 0.5
        assert!((mw[1][4] - 0.5).abs() < 0.1, "D weight should be ~0.5");
        assert!((mw[1][5] - 0.5).abs() < 0.1, "E weight should be ~0.5");
    }

    #[test]
    fn test_column_pseudocounts_conserved() {
        let std_prob = std_prob_ncbistdaa();
        // A highly conserved column (all weight on one residue) should have
        // low pseudocount weight (high information -> low pseudo)
        let mut mw = [0.0f64; AA_SIZE];
        mw[1] = 1.0; // All weight on Alanine
        let pseudo = compute_column_pseudocounts(&mw, &std_prob, 50.0);
        // Should be relatively small
        assert!(pseudo < 100.0, "Conserved column pseudo should be small, got {}", pseudo);
    }

    #[test]
    fn test_column_pseudocounts_uniform() {
        let std_prob = std_prob_ncbistdaa();
        // A uniform column should have high pseudocount weight (low information)
        let mut mw = [0.0f64; AA_SIZE];
        for &code in STD_RESIDUES.iter() {
            mw[code as usize] = 1.0 / EFFECTIVE_ALPHABET as f64;
        }
        let pseudo = compute_column_pseudocounts(&mw, &std_prob, 50.0);
        // Should be larger than for conserved column
        assert!(pseudo > 10.0, "Uniform column pseudo should be substantial, got {}", pseudo);
    }

    #[test]
    fn test_update_from_alignment_changes_scores() {
        let query = vec![1u8, 4, 7, 10, 13]; // A, D, G, K, N
        let matrix = simple_matrix();
        let mut pssm = Pssm::from_sequence(&query, &matrix);
        let original_scores = pssm.scores.clone();

        // Create aligned sequences that mostly match query but with some variation
        let aligned = vec![
            vec![1u8, 4, 7, 10, 13],  // exact match
            vec![1u8, 4, 7, 10, 13],  // exact match
            vec![1u8, 5, 7, 10, 4],   // E at pos 1, D at pos 4
            vec![1u8, 4, 8, 10, 13],  // H at pos 2
        ];

        let bg = crate::matrix::AA_FREQUENCIES;
        pssm.update_from_alignment(&aligned, &bg, 10.0);

        // Scores should have changed
        assert_ne!(pssm.scores, original_scores,
                   "PSSM scores should change after alignment update");
    }

    #[test]
    fn test_effective_observations() {
        let std_prob = std_prob_ncbistdaa();

        // With many diverse sequences, effective observations should be substantial
        let seqs: Vec<Vec<u8>> = (0..20).map(|i| {
            vec![STD_RESIDUES[i % 20]]
        }).collect();

        let obs = compute_effective_observations(&seqs, 0, 1, &std_prob);
        assert!(obs > 0.0, "Should have positive effective observations");
    }
}
