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

pub const PSI_SUCCESS: i32 = 0;
pub const PSIERR_BADPARAM: i32 = -1;
pub const PSIERR_OUTOFMEM: i32 = -2;
pub const PSIERR_BADSEQWEIGHTS: i32 = -3;
pub const PSIERR_POSITIVEAVGSCORE: i32 = -5;
pub const PSIERR_NOALIGNEDSEQS: i32 = -6;
pub const PSIERR_GAPINQUERY: i32 = -7;
pub const PSIERR_UNALIGNEDCOLUMN: i32 = -8;
pub const PSIERR_COLUMNOFGAPS: i32 = -9;
pub const PSIERR_STARTINGGAP: i32 = -10;
pub const PSIERR_ENDINGGAP: i32 = -11;
pub const PSIERR_BADPROFILE: i32 = -12;
pub const PSIERR_UNKNOWN: i32 = -255;

pub const K_PSI_NEAR_IDENTICAL: f64 = 0.94;
pub const K_PSI_IDENTICAL: f64 = 1.0;
pub const K_QUERY_INDEX: usize = 0;
pub const K_PSI_SCALE_FACTOR: f64 = 200.0;
pub const K_POSIT_SCALING_PERCENT: f64 = 0.05;
pub const K_POSIT_SCALING_NUM_ITERATIONS: usize = 10;

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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PSIMsaDimensions {
    pub query_length: u32,
    pub num_seqs: u32,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct PSIMsaCell {
    pub letter: u8,
    pub is_aligned: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PSIMsa {
    pub dimensions: PSIMsaDimensions,
    pub data: Vec<Vec<PSIMsaCell>>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct PSIMatrix {
    pub ncols: u32,
    pub nrows: u32,
    pub pssm: Vec<Vec<i32>>,
    pub lambda: f64,
    pub kappa: f64,
    pub h: f64,
    pub ung_lambda: f64,
    pub ung_kappa: f64,
    pub ung_h: f64,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct PSIDiagnosticsRequest {
    pub information_content: bool,
    pub residue_frequencies: bool,
    pub weighted_residue_frequencies: bool,
    pub frequency_ratios: bool,
    pub gapless_column_weights: bool,
    pub sigma: bool,
    pub interval_sizes: bool,
    pub num_matching_seqs: bool,
    pub independent_observations: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub struct PSIDiagnosticsResponse {
    pub information_content: Option<Vec<f64>>,
    pub residue_freqs: Option<Vec<Vec<u32>>>,
    pub weighted_residue_freqs: Option<Vec<Vec<f64>>>,
    pub frequency_ratios: Option<Vec<Vec<f64>>>,
    pub gapless_column_weights: Option<Vec<f64>>,
    pub sigma: Option<Vec<f64>>,
    pub interval_sizes: Option<Vec<u32>>,
    pub num_matching_seqs: Option<Vec<u32>>,
    pub query_length: u32,
    pub alphabet_size: u32,
    pub independent_observations: Option<Vec<f64>>,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct PsiPackedMsaCell {
    pub letter: u8,
    pub is_aligned: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PsiPackedMsa {
    pub dimensions: PSIMsaDimensions,
    pub data: Vec<Vec<PsiPackedMsaCell>>,
    pub use_sequence: Vec<bool>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PsiPurgeFsmState {
    Counting,
    Resting,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct PsiAlignmentTraits {
    pub start: u32,
    pub effective_length: u32,
    pub n_x_residues: u32,
    pub n_identical: u32,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SSeqRange {
    pub left: i32,
    pub right: i32,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PsiMsaCell {
    pub letter: u8,
    pub is_aligned: bool,
    pub extents: SSeqRange,
}

impl Default for PsiMsaCell {
    fn default() -> Self {
        Self {
            letter: 0,
            is_aligned: false,
            extents: SSeqRange { left: -1, right: 0 },
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PsiMsa {
    pub dimensions: PSIMsaDimensions,
    pub cell: Vec<Vec<PsiMsaCell>>,
    pub query: Vec<u8>,
    pub residue_counts: Vec<Vec<u32>>,
    pub alphabet_size: u32,
    pub num_matching_seqs: Vec<u32>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct PsiInternalPssmData {
    pub ncols: u32,
    pub nrows: u32,
    pub pssm: Vec<Vec<i32>>,
    pub scaled_pssm: Vec<Vec<i32>>,
    pub freq_ratios: Vec<Vec<f64>>,
    pub pseudocounts: Vec<f64>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct PSICdMsaCellData {
    pub wfreqs: Vec<f64>,
    pub iobsr: f64,
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct PSICdMsaCell {
    pub is_aligned: bool,
    pub data: Option<PSICdMsaCellData>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct PSICdMsa {
    pub dimensions: Option<PSIMsaDimensions>,
    pub query: Vec<u8>,
    pub msa: Vec<Vec<PSICdMsaCell>>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PsiAlignedBlock {
    pub pos_extnt: Vec<SSeqRange>,
    pub size: Vec<u32>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct PsiSequenceWeights {
    pub match_weights: Vec<Vec<f64>>,
    pub match_weights_size: u32,
    pub norm_seq_weights: Vec<f64>,
    pub row_sigma: Vec<f64>,
    pub sigma: Vec<f64>,
    pub std_prob: Vec<f64>,
    pub gapless_column_weights: Vec<f64>,
    pub pos_distinct_distrib: Vec<Vec<i32>>,
    pub pos_distinct_distrib_size: u32,
    pub pos_num_participating: Vec<i32>,
    pub independent_observations: Vec<f64>,
}

/// blast-rs public allocator for `PSIMsa`.
pub fn psi_public_msa_new(dimensions: Option<&PSIMsaDimensions>) -> Option<PSIMsa> {
    let dimensions = *dimensions?;
    let rows = dimensions.num_seqs.saturating_add(1) as usize;
    let cols = dimensions.query_length as usize;
    Some(PSIMsa {
        dimensions,
        data: vec![vec![PSIMsaCell::default(); cols]; rows],
    })
}

/// blast-rs public deallocator for `PSIMsa`.
pub fn psi_public_msa_free(mut msa: Option<PSIMsa>) -> Option<PSIMsa> {
    if let Some(msa) = msa.as_mut() {
        for row in &mut msa.data {
            row.clear();
        }
        msa.data.clear();
        msa.dimensions.query_length = 0;
        msa.dimensions.num_seqs = 0;
    }
    None
}

fn psi_public_msa_has_shape(msa: &PSIMsa) -> bool {
    let rows = msa.dimensions.num_seqs.saturating_add(1) as usize;
    let cols = msa.dimensions.query_length as usize;
    msa.data.len() >= rows
        && msa.data.iter().take(rows).all(|row| {
            row.len() >= cols
                && row
                    .iter()
                    .take(cols)
                    .all(|cell| cell.letter as usize <= crate::encoding::BLASTAA_SIZE)
        })
}

/// Port of NCBI `PSIMatrixNew` (`blast_psi.c:541`).
pub fn psi_matrix_new(query_length: u32, alphabet_size: u32) -> Option<PSIMatrix> {
    Some(PSIMatrix {
        ncols: query_length,
        nrows: alphabet_size,
        pssm: vec![vec![0; alphabet_size as usize]; query_length as usize],
        lambda: 0.0,
        kappa: 0.0,
        h: 0.0,
        ung_lambda: 0.0,
        ung_kappa: 0.0,
        ung_h: 0.0,
    })
}

/// Rust ownership equivalent of NCBI `PSIMatrixFree`.
pub fn psi_matrix_free(_: Option<PSIMatrix>) -> Option<PSIMatrix> {
    None
}

/// Port of NCBI static `s_PSISavePssm` (`blast_psi.c:440`).
pub fn s_psi_save_pssm(
    internal_pssm: Option<&PsiInternalPssmData>,
    sbp: Option<&crate::stat::BlastScoreBlk>,
    pssm: Option<&mut PSIMatrix>,
) -> i32 {
    let (Some(internal_pssm), Some(sbp), Some(pssm)) = (internal_pssm, sbp, pssm) else {
        return PSIERR_BADPARAM;
    };

    if internal_pssm.ncols != pssm.ncols || internal_pssm.nrows != pssm.nrows {
        return PSIERR_BADPARAM;
    }

    let rows = pssm.ncols as usize;
    let cols = pssm.nrows as usize;
    if internal_pssm.pssm.len() < rows || pssm.pssm.len() < rows {
        return PSIERR_BADPARAM;
    }
    for i in 0..rows {
        if internal_pssm.pssm[i].len() < cols || pssm.pssm[i].len() < cols {
            return PSIERR_BADPARAM;
        }
        pssm.pssm[i][..cols].copy_from_slice(&internal_pssm.pssm[i][..cols]);
    }

    let (Some(gap), Some(ungapped)) = (sbp.kbp_gap_psi.first(), sbp.kbp_psi.first()) else {
        return PSIERR_BADPARAM;
    };
    pssm.lambda = gap.lambda;
    pssm.kappa = gap.k;
    pssm.h = gap.h;
    pssm.ung_lambda = ungapped.lambda;
    pssm.ung_kappa = ungapped.k;
    pssm.ung_h = ungapped.h;

    PSI_SUCCESS
}

fn psi_create_pssm_cleanup(
    pssm: Option<&mut Option<PSIMatrix>>,
    diagnostics: Option<&mut Option<PSIDiagnosticsResponse>>,
) {
    if let Some(pssm) = pssm {
        *pssm = psi_matrix_free(pssm.take());
    }
    if let Some(diagnostics) = diagnostics {
        *diagnostics = psi_diagnostics_response_free(diagnostics.take());
    }
}

/// Port of NCBI `PSICreatePssmWithDiagnostics` (`blast_psi.c:105`).
///
/// This is the public PSI-BLAST PSSM constructor. It deliberately calls the
/// translated private stages in the same order as C instead of delegating to the
/// standalone [`Pssm`] builder.
pub fn psi_create_pssm_with_diagnostics(
    msap: Option<&PSIMsa>,
    options: Option<&crate::options::PSIBlastOptions>,
    sbp: Option<&mut crate::stat::BlastScoreBlk>,
    request: Option<&PSIDiagnosticsRequest>,
    pssm: Option<&mut Option<PSIMatrix>>,
    diagnostics: Option<&mut Option<PSIDiagnosticsResponse>>,
) -> i32 {
    let (Some(msap), Some(options), Some(sbp), Some(pssm_out)) = (msap, options, sbp, pssm) else {
        return PSIERR_BADPARAM;
    };
    *pssm_out = None;
    let mut diagnostics = diagnostics;
    if let Some(diags) = diagnostics.as_deref_mut() {
        *diags = None;
    }
    if !psi_public_msa_has_shape(msap) {
        return PSIERR_BADPARAM;
    }

    let Some(mut packed_msa) = psi_packed_msa_new(Some(msap)) else {
        return PSIERR_OUTOFMEM;
    };

    let mut status = psi_purge_biased_segments(Some(&mut packed_msa));
    if status != PSI_SUCCESS {
        return status;
    }

    let Some(mut msa) = psi_msa_new(Some(&packed_msa), sbp.alphabet_size as u32) else {
        return PSIERR_OUTOFMEM;
    };
    let Some(mut aligned_block) = psi_aligned_block_new(msa.dimensions.query_length) else {
        return PSIERR_OUTOFMEM;
    };
    let Some(mut seq_weights) = psi_sequence_weights_new(Some(&msa.dimensions), sbp.alphabet_size)
    else {
        return PSIERR_OUTOFMEM;
    };
    let Some(mut internal_pssm) =
        psi_internal_pssm_data_new(msa.dimensions.query_length, sbp.alphabet_size as u32)
    else {
        return PSIERR_OUTOFMEM;
    };
    let Some(out_matrix) = psi_matrix_new(msa.dimensions.query_length, sbp.alphabet_size as u32)
    else {
        return PSIERR_OUTOFMEM;
    };
    *pssm_out = Some(out_matrix);
    let _ = psi_packed_msa_free(Some(packed_msa));

    if options.nsg_compatibility_mode {
        psi_structure_group_customization(&mut msa);
        status = psi_validate_msa_structure_group(Some(&msa));
    } else {
        status = psi_validate_msa(Some(&msa), options.ignore_unaligned_positions);
    }
    if status != PSI_SUCCESS {
        psi_create_pssm_cleanup(Some(pssm_out), diagnostics.as_deref_mut());
        return status;
    }

    status = psi_compute_alignment_blocks(Some(&mut msa), Some(&mut aligned_block));
    if status != PSI_SUCCESS {
        psi_create_pssm_cleanup(Some(pssm_out), diagnostics.as_deref_mut());
        return status;
    }

    status = psi_compute_sequence_weights(
        Some(&msa),
        Some(&aligned_block),
        options.nsg_compatibility_mode,
        Some(&mut seq_weights),
    );
    if status != PSI_SUCCESS {
        psi_create_pssm_cleanup(Some(pssm_out), diagnostics.as_deref_mut());
        return status;
    }

    status = psi_compute_freq_ratios(
        Some(&msa),
        Some(&seq_weights),
        Some(&*sbp),
        Some(&aligned_block),
        options.pseudo_count,
        options.nsg_compatibility_mode,
        Some(&mut internal_pssm),
    );
    if status != PSI_SUCCESS {
        psi_create_pssm_cleanup(Some(pssm_out), diagnostics.as_deref_mut());
        return status;
    }

    status = psi_create_and_scale_pssm_from_frequency_ratios(
        Some(&mut internal_pssm),
        Some(&msa.query),
        Some(&seq_weights.std_prob),
        Some(&mut *sbp),
        options.impala_scaling_factor,
    );
    if status != PSI_SUCCESS {
        psi_create_pssm_cleanup(Some(pssm_out), diagnostics.as_deref_mut());
        return status;
    }

    status = s_psi_save_pssm(Some(&internal_pssm), Some(&*sbp), pssm_out.as_mut());
    if status != PSI_SUCCESS {
        psi_create_pssm_cleanup(Some(pssm_out), diagnostics.as_deref_mut());
        return status;
    }

    if let (Some(request), Some(diags_out)) = (request, diagnostics.as_deref_mut()) {
        let Some(mut response) = psi_diagnostics_response_new(
            msa.dimensions.query_length,
            sbp.alphabet_size as u32,
            Some(request),
        ) else {
            psi_create_pssm_cleanup(Some(pssm_out), Some(diags_out));
            return PSIERR_OUTOFMEM;
        };
        status = psi_save_diagnostics(
            Some(&msa),
            Some(&aligned_block),
            Some(&seq_weights),
            Some(&internal_pssm),
            Some(&mut response),
        );
        if status != PSI_SUCCESS {
            psi_create_pssm_cleanup(Some(pssm_out), Some(diags_out));
            return status;
        }
        *diags_out = Some(response);
    }

    PSI_SUCCESS
}

/// Port of NCBI `PSIDiagnosticsRequestNew`.
pub fn psi_diagnostics_request_new() -> PSIDiagnosticsRequest {
    PSIDiagnosticsRequest::default()
}

/// Port of NCBI `PSIDiagnosticsRequestNewEx`.
pub fn psi_diagnostics_request_new_ex(save_ascii_pssm: bool) -> PSIDiagnosticsRequest {
    let mut retval = psi_diagnostics_request_new();
    retval.frequency_ratios = true;
    if save_ascii_pssm {
        retval.information_content = true;
        retval.weighted_residue_frequencies = true;
        retval.gapless_column_weights = true;
        retval.sigma = true;
        retval.interval_sizes = true;
        retval.num_matching_seqs = true;
    }
    retval
}

/// Rust ownership equivalent of NCBI `PSIDiagnosticsRequestFree`.
pub fn psi_diagnostics_request_free(
    _: Option<PSIDiagnosticsRequest>,
) -> Option<PSIDiagnosticsRequest> {
    None
}

/// Port of NCBI `PSIDiagnosticsResponseNew`.
pub fn psi_diagnostics_response_new(
    query_length: u32,
    alphabet_size: u32,
    wants: Option<&PSIDiagnosticsRequest>,
) -> Option<PSIDiagnosticsResponse> {
    let wants = wants?;
    let q = query_length as usize;
    let a = alphabet_size as usize;
    Some(PSIDiagnosticsResponse {
        information_content: wants.information_content.then(|| vec![0.0; q]),
        residue_freqs: wants.residue_frequencies.then(|| vec![vec![0; a]; q]),
        weighted_residue_freqs: wants
            .weighted_residue_frequencies
            .then(|| vec![vec![0.0; a]; q]),
        frequency_ratios: wants.frequency_ratios.then(|| vec![vec![0.0; a]; q]),
        gapless_column_weights: wants.gapless_column_weights.then(|| vec![0.0; q]),
        sigma: wants.sigma.then(|| vec![0.0; q]),
        interval_sizes: wants.interval_sizes.then(|| vec![0; q]),
        num_matching_seqs: wants.num_matching_seqs.then(|| vec![0; q]),
        query_length,
        alphabet_size,
        independent_observations: wants.independent_observations.then(|| vec![0.0; q]),
    })
}

/// Rust ownership equivalent of NCBI `PSIDiagnosticsResponseFree`.
pub fn psi_diagnostics_response_free(
    mut diags: Option<PSIDiagnosticsResponse>,
) -> Option<PSIDiagnosticsResponse> {
    if let Some(diags) = diags.as_mut() {
        diags.information_content = None;
        diags.residue_freqs = None;
        diags.weighted_residue_freqs = None;
        diags.frequency_ratios = None;
        diags.gapless_column_weights = None;
        diags.sigma = None;
        diags.interval_sizes = None;
        diags.num_matching_seqs = None;
        diags.independent_observations = None;
        diags.query_length = 0;
        diags.alphabet_size = 0;
    }
    None
}

fn psi_diag_vec_has_len<T>(values: &Option<Vec<T>>, query_length: usize) -> bool {
    values
        .as_ref()
        .is_none_or(|values| values.len() >= query_length)
}

fn psi_diag_matrix_has_shape<T>(
    values: &Option<Vec<Vec<T>>>,
    query_length: usize,
    alphabet_size: usize,
) -> bool {
    values.as_ref().is_none_or(|rows| {
        rows.len() >= query_length
            && rows
                .iter()
                .take(query_length)
                .all(|row| row.len() >= alphabet_size)
    })
}

fn psi_diagnostics_response_buffers_are_large_enough(
    diagnostics: &PSIDiagnosticsResponse,
    query_length: usize,
    alphabet_size: usize,
) -> bool {
    psi_diag_vec_has_len(&diagnostics.information_content, query_length)
        && psi_diag_matrix_has_shape(&diagnostics.residue_freqs, query_length, alphabet_size)
        && psi_diag_matrix_has_shape(
            &diagnostics.weighted_residue_freqs,
            query_length,
            alphabet_size,
        )
        && psi_diag_matrix_has_shape(&diagnostics.frequency_ratios, query_length, alphabet_size)
        && psi_diag_vec_has_len(&diagnostics.gapless_column_weights, query_length)
        && psi_diag_vec_has_len(&diagnostics.sigma, query_length)
        && psi_diag_vec_has_len(&diagnostics.interval_sizes, query_length)
        && psi_diag_vec_has_len(&diagnostics.num_matching_seqs, query_length)
        && psi_diag_vec_has_len(&diagnostics.independent_observations, query_length)
}

fn psi_matrix_has_shape<T>(rows: &[Vec<T>], query_length: usize, alphabet_size: usize) -> bool {
    rows.len() >= query_length
        && rows
            .iter()
            .take(query_length)
            .all(|row| row.len() >= alphabet_size)
}

fn psi_save_diagnostics_sources_are_large_enough(
    msa: &PsiMsa,
    aligned_block: &PsiAlignedBlock,
    seq_weights: &PsiSequenceWeights,
    internal_pssm: &PsiInternalPssmData,
    diagnostics: &PSIDiagnosticsResponse,
    query_length: usize,
    alphabet_size: usize,
) -> bool {
    if diagnostics.information_content.is_some()
        && (seq_weights.std_prob.len() < alphabet_size
            || !psi_matrix_has_shape(&internal_pssm.freq_ratios, query_length, alphabet_size))
    {
        return false;
    }
    if diagnostics.residue_freqs.is_some()
        && !psi_matrix_has_shape(&msa.residue_counts, query_length, alphabet_size)
    {
        return false;
    }
    if diagnostics.weighted_residue_freqs.is_some()
        && !psi_matrix_has_shape(&seq_weights.match_weights, query_length, alphabet_size)
    {
        return false;
    }
    if diagnostics.frequency_ratios.is_some()
        && !psi_matrix_has_shape(&internal_pssm.freq_ratios, query_length, alphabet_size)
    {
        return false;
    }
    if diagnostics.gapless_column_weights.is_some()
        && (msa.num_matching_seqs.len() < query_length
            || msa.cell.is_empty()
            || msa.cell[0].len() < query_length
            || seq_weights.gapless_column_weights.len() < query_length
            || seq_weights.sigma.len() < query_length
            || aligned_block.size.len() < query_length
            || internal_pssm.pseudocounts.len() < query_length)
    {
        return false;
    }
    if diagnostics.sigma.is_some() && seq_weights.sigma.len() < query_length {
        return false;
    }
    if diagnostics.interval_sizes.is_some() && aligned_block.size.len() < query_length {
        return false;
    }
    if diagnostics.num_matching_seqs.is_some() && msa.num_matching_seqs.len() < query_length {
        return false;
    }
    if diagnostics.independent_observations.is_some()
        && seq_weights.independent_observations.len() < query_length
    {
        return false;
    }
    true
}

fn psi_save_cd_diagnostics_sources_are_large_enough(
    seq_weights: &PsiSequenceWeights,
    internal_pssm: &PsiInternalPssmData,
    diagnostics: &PSIDiagnosticsResponse,
    query_length: usize,
    alphabet_size: usize,
) -> bool {
    if diagnostics.information_content.is_some()
        && (seq_weights.std_prob.len() < alphabet_size
            || !psi_matrix_has_shape(&internal_pssm.freq_ratios, query_length, alphabet_size))
    {
        return false;
    }
    if diagnostics.weighted_residue_freqs.is_some()
        && !psi_matrix_has_shape(&seq_weights.match_weights, query_length, alphabet_size)
    {
        return false;
    }
    if diagnostics.frequency_ratios.is_some()
        && !psi_matrix_has_shape(&internal_pssm.freq_ratios, query_length, alphabet_size)
    {
        return false;
    }
    if diagnostics.independent_observations.is_some()
        && seq_weights.independent_observations.len() < query_length
    {
        return false;
    }
    true
}

/// Port of NCBI `_PSIPackedMsaNew` (`blast_psi_priv.c:129`).
pub fn psi_packed_msa_new(msa: Option<&PSIMsa>) -> Option<PsiPackedMsa> {
    let msa = msa?;
    Some(PsiPackedMsa {
        dimensions: msa.dimensions,
        data: msa
            .data
            .iter()
            .map(|row| {
                row.iter()
                    .map(|cell| PsiPackedMsaCell {
                        letter: cell.letter,
                        is_aligned: cell.is_aligned,
                    })
                    .collect()
            })
            .collect(),
        use_sequence: vec![true; msa.dimensions.num_seqs.saturating_add(1) as usize],
    })
}

/// Rust ownership equivalent of NCBI `_PSIPackedMsaFree`.
pub fn psi_packed_msa_free(mut msa: Option<PsiPackedMsa>) -> Option<PsiPackedMsa> {
    if let Some(msa) = msa.as_mut() {
        for row in &mut msa.data {
            row.clear();
        }
        msa.data.clear();
        msa.use_sequence.clear();
        msa.dimensions.query_length = 0;
        msa.dimensions.num_seqs = 0;
    }
    None
}

/// Port of NCBI `_PSIPackedMsaGetNumberOfAlignedSeqs`.
pub fn psi_packed_msa_get_number_of_aligned_seqs(msa: Option<&PsiPackedMsa>) -> u32 {
    msa.map(|msa| msa.use_sequence.iter().filter(|&&used| used).count() as u32)
        .unwrap_or(0)
        .saturating_sub(1)
}

/// Port of NCBI static `s_PSIDiscardIfUnused`.
pub fn s_psi_discard_if_unused(msa: &mut PsiPackedMsa, seq_index: usize) {
    let Some(row) = msa.data.get(seq_index) else {
        return;
    };
    let query_length = msa.dimensions.query_length as usize;
    let contains_aligned_regions = row.iter().take(query_length).any(|cell| cell.is_aligned);

    if !contains_aligned_regions {
        if let Some(use_sequence) = msa.use_sequence.get_mut(seq_index) {
            *use_sequence = false;
        }
    }
}

/// Port of NCBI `_PSIPurgeAlignedRegion`.
pub fn psi_purge_aligned_region(
    msa: Option<&mut PsiPackedMsa>,
    seq_index: u32,
    start: u32,
    stop: u32,
) -> i32 {
    let Some(msa) = msa else {
        return PSIERR_BADPARAM;
    };
    if seq_index == 0
        || seq_index > msa.dimensions.num_seqs.saturating_add(1)
        || stop > msa.dimensions.query_length
    {
        return PSIERR_BADPARAM;
    }

    let seq_index = seq_index as usize;
    let Some(row) = msa.data.get_mut(seq_index) else {
        return PSIERR_BADPARAM;
    };
    for cell in row
        .iter_mut()
        .skip(start as usize)
        .take(stop.saturating_sub(start) as usize)
    {
        cell.letter = 0;
        cell.is_aligned = false;
    }

    s_psi_discard_if_unused(msa, seq_index);
    PSI_SUCCESS
}

/// Debug string equivalent of C `DEBUG_printTraits`.
pub fn debug_print_traits(
    traits: &PsiAlignmentTraits,
    state: PsiPurgeFsmState,
    position: u32,
) -> String {
    let state = match state {
        PsiPurgeFsmState::Counting => "eCounting",
        PsiPurgeFsmState::Resting => "eResting",
    };
    format!(
        "Position: {position} - State: {state}\n\tstart: {}\n\teffective_length: {}\n\tn_x_residues: {}\n\tn_identical: {}\n",
        traits.start, traits.effective_length, traits.n_x_residues, traits.n_identical
    )
}

pub fn psi_reset_alignment_traits(traits: &mut PsiAlignmentTraits, position: u32) {
    *traits = PsiAlignmentTraits {
        start: position,
        ..PsiAlignmentTraits::default()
    };
}

/// Port of NCBI inline `_handleNeitherAligned`.
pub fn handle_neither_aligned(
    traits: &PsiAlignmentTraits,
    state: &mut PsiPurgeFsmState,
    msa: &mut PsiPackedMsa,
    seq_index: u32,
    max_percent_identity: f64,
) {
    if *state == PsiPurgeFsmState::Counting {
        if traits.effective_length > 0 {
            let percent_identity = traits.n_identical as f64 / traits.effective_length as f64;
            if percent_identity >= max_percent_identity {
                let align_stop = traits.start + traits.effective_length + traits.n_x_residues;
                let _ = psi_purge_aligned_region(Some(msa), seq_index, traits.start, align_stop);
            }
        }
        *state = PsiPurgeFsmState::Resting;
    }
}

pub fn handle_both_aligned_same_residue_no_x(
    traits: &mut PsiAlignmentTraits,
    state: PsiPurgeFsmState,
) {
    if state == PsiPurgeFsmState::Counting {
        traits.n_identical += 1;
    }
}

/// Port of NCBI inline `_handleEitherAlignedEitherX`.
pub fn handle_either_aligned_either_x(traits: &mut PsiAlignmentTraits, state: PsiPurgeFsmState) {
    if state == PsiPurgeFsmState::Counting {
        traits.n_x_residues += 1;
    }
}

pub fn handle_either_aligned_neither_x(
    traits: &mut PsiAlignmentTraits,
    state: &mut PsiPurgeFsmState,
    position: u32,
) {
    match *state {
        PsiPurgeFsmState::Counting => traits.effective_length += 1,
        PsiPurgeFsmState::Resting => {
            psi_reset_alignment_traits(traits, position);
            traits.effective_length = 1;
            *state = PsiPurgeFsmState::Counting;
        }
    }
}

/// Port of NCBI static `s_PSIPurgeSimilarAlignments`.
pub fn s_psi_purge_similar_alignments(
    msa: &mut PsiPackedMsa,
    seq_index1: u32,
    seq_index2: u32,
    max_percent_identity: f64,
) {
    let seq_index1 = seq_index1 as usize;
    let seq_index2 = seq_index2 as usize;
    if seq_index1 == seq_index2
        || !msa.use_sequence.get(seq_index1).copied().unwrap_or(false)
        || !msa.use_sequence.get(seq_index2).copied().unwrap_or(false)
    {
        return;
    }

    let Some(seq1) = msa.data.get(seq_index1).cloned() else {
        return;
    };
    let Some(seq2) = msa.data.get(seq_index2).cloned() else {
        return;
    };

    let x_residue = crate::encoding::AMINOACID_TO_NCBISTDAA[b'X' as usize];
    let query_length = msa.dimensions.query_length as usize;
    let mut state = PsiPurgeFsmState::Counting;
    let mut traits = PsiAlignmentTraits::default();
    psi_reset_alignment_traits(&mut traits, 0);

    for p in 0..query_length {
        let Some(cell1) = seq1.get(p).copied() else {
            break;
        };
        let Some(cell2) = seq2.get(p).copied() else {
            break;
        };
        let pos1_aligned = if seq_index1 == K_QUERY_INDEX {
            false
        } else {
            cell1.is_aligned
        };
        let pos2_aligned = cell2.is_aligned;

        if !pos1_aligned && !pos2_aligned {
            handle_neither_aligned(
                &traits,
                &mut state,
                msa,
                seq_index2 as u32,
                max_percent_identity,
            );
        } else {
            let neither_is_x = cell1.letter != x_residue && cell2.letter != x_residue;
            if neither_is_x {
                handle_either_aligned_neither_x(&mut traits, &mut state, p as u32);
            } else {
                handle_either_aligned_either_x(&mut traits, state);
            }

            if neither_is_x && pos2_aligned && cell1.is_aligned && cell1.letter == cell2.letter {
                handle_both_aligned_same_residue_no_x(&mut traits, state);
            }
        }
    }

    handle_neither_aligned(
        &traits,
        &mut state,
        msa,
        seq_index2 as u32,
        max_percent_identity,
    );
}

/// Port of NCBI static `s_PSIPurgeSelfHits`.
pub fn s_psi_purge_self_hits(msa: &mut PsiPackedMsa) {
    for s in (K_QUERY_INDEX + 1)..msa.dimensions.num_seqs.saturating_add(1) as usize {
        s_psi_purge_similar_alignments(msa, K_QUERY_INDEX as u32, s as u32, K_PSI_IDENTICAL);
    }
}

pub fn s_psi_purge_near_identical_alignments(msa: &mut PsiPackedMsa) {
    let limit = msa.dimensions.num_seqs.saturating_add(1) as usize;
    for i in 1..limit {
        for j in 1..limit.saturating_sub(i) {
            s_psi_purge_similar_alignments(msa, j as u32, (i + j) as u32, K_PSI_NEAR_IDENTICAL);
        }
    }
}

/// Port of NCBI `_PSIPurgeBiasedSegments`.
pub fn psi_purge_biased_segments(msa: Option<&mut PsiPackedMsa>) -> i32 {
    let Some(msa) = msa else {
        return PSIERR_BADPARAM;
    };
    s_psi_purge_self_hits(msa);
    s_psi_purge_near_identical_alignments(msa);
    PSI_SUCCESS
}

/// Port of NCBI `_PSIMsaNew` (`blast_psi_priv.c:308`).
pub fn psi_msa_new(msa: Option<&PsiPackedMsa>, alphabet_size: u32) -> Option<PsiMsa> {
    let msa = msa?;
    let num_seqs = psi_packed_msa_get_number_of_aligned_seqs(Some(msa));
    let dimensions = PSIMsaDimensions {
        query_length: msa.dimensions.query_length,
        num_seqs,
    };
    let rows = dimensions.num_seqs.saturating_add(1) as usize;
    let cols = dimensions.query_length as usize;
    let mut cell = vec![
        vec![
            PsiMsaCell {
                extents: SSeqRange {
                    left: -1,
                    right: dimensions.query_length as i32,
                },
                ..PsiMsaCell::default()
            };
            cols
        ];
        rows
    ];
    let mut out_row = 0usize;
    for (row_index, row) in msa.data.iter().enumerate() {
        if !msa.use_sequence.get(row_index).copied().unwrap_or(false) {
            continue;
        }
        for (col, packed) in row.iter().take(cols).enumerate() {
            cell[out_row][col].letter = packed.letter;
            cell[out_row][col].is_aligned = packed.is_aligned;
        }
        out_row += 1;
        if out_row >= rows {
            break;
        }
    }
    let query = msa
        .data
        .first()
        .map(|row| row.iter().take(cols).map(|cell| cell.letter).collect())
        .unwrap_or_else(|| vec![0; cols]);
    let mut retval = PsiMsa {
        dimensions,
        cell,
        query,
        residue_counts: vec![vec![0; alphabet_size as usize]; cols],
        alphabet_size,
        num_matching_seqs: vec![0; cols],
    };
    psi_update_position_counts(&mut retval);
    Some(retval)
}

/// Rust ownership equivalent of NCBI `_PSIMsaFree`.
pub fn psi_msa_free(mut msa: Option<PsiMsa>) -> Option<PsiMsa> {
    if let Some(msa) = msa.as_mut() {
        for row in &mut msa.cell {
            row.clear();
        }
        msa.cell.clear();
        msa.query.clear();
        for row in &mut msa.residue_counts {
            row.clear();
        }
        msa.residue_counts.clear();
        msa.num_matching_seqs.clear();
        msa.dimensions.query_length = 0;
        msa.dimensions.num_seqs = 0;
        msa.alphabet_size = 0;
    }
    None
}

fn psi_update_position_counts(msa: &mut PsiMsa) {
    for counts in &mut msa.residue_counts {
        counts.fill(0);
    }
    msa.num_matching_seqs.fill(0);

    for position in 0..msa.dimensions.query_length as usize {
        let mut matching = 0u32;
        for row in &msa.cell {
            if row.get(position).is_some_and(|cell| cell.is_aligned) {
                matching += 1;
                let residue = row[position].letter as usize;
                if residue < msa.alphabet_size as usize {
                    msa.residue_counts[position][residue] += 1;
                }
            }
        }
        msa.num_matching_seqs[position] = matching;
    }
}

/// Port of NCBI `_PSIStructureGroupCustomization`.
pub fn psi_structure_group_customization(msa: &mut PsiMsa) {
    let query_length = msa.dimensions.query_length as usize;
    if let Some(query_row) = msa.cell.get_mut(0) {
        for cell in query_row.iter_mut().take(query_length) {
            cell.letter = 0;
            cell.is_aligned = false;
        }
    }
    psi_update_position_counts(msa);
}

/// Port of NCBI static `_PSIGetLeftExtents`.
pub fn psi_get_left_extents(msa: &mut PsiMsa, seq_index: u32) {
    let seq_index = seq_index as usize;
    let query_length = msa.dimensions.query_length as usize;
    if query_length == 0 {
        return;
    }
    let gap = crate::encoding::AMINOACID_TO_NCBISTDAA[b'-' as usize];
    let Some(row) = msa.cell.get_mut(seq_index) else {
        return;
    };
    if row[0].is_aligned && row[0].letter != gap {
        row[0].extents.left = 0;
    }
    for curr in 1..query_length {
        let prev = curr - 1;
        if !row[curr].is_aligned {
            continue;
        }
        if row[prev].is_aligned {
            row[curr].extents.left = row[prev].extents.left;
        } else {
            row[curr].extents.left = curr as i32;
        }
    }
}

/// Port of NCBI static `_PSIGetRightExtents`.
pub fn psi_get_right_extents(msa: &mut PsiMsa, seq_index: u32) {
    let seq_index = seq_index as usize;
    let query_length = msa.dimensions.query_length as usize;
    if query_length == 0 {
        return;
    }
    let gap = crate::encoding::AMINOACID_TO_NCBISTDAA[b'-' as usize];
    let Some(row) = msa.cell.get_mut(seq_index) else {
        return;
    };
    let last = query_length - 1;
    if row[last].is_aligned && row[last].letter != gap {
        row[last].extents.right = last as i32;
    }
    for curr in (0..last).rev() {
        let next = curr + 1;
        if !row[curr].is_aligned {
            continue;
        }
        if row[next].is_aligned {
            row[curr].extents.right = row[next].extents.right;
        } else {
            row[curr].extents.right = curr as i32;
        }
    }
}

pub fn psi_compute_position_extents(
    msa: &PsiMsa,
    seq_index: u32,
    aligned_blocks: &mut PsiAlignedBlock,
) {
    let seq_index = seq_index as usize;
    let Some(row) = msa.cell.get(seq_index) else {
        return;
    };
    let query_length = msa.dimensions.query_length as usize;
    for i in 0..query_length {
        if row.get(i).is_some_and(|cell| cell.is_aligned) {
            aligned_blocks.pos_extnt[i].left =
                aligned_blocks.pos_extnt[i].left.max(row[i].extents.left);
            aligned_blocks.pos_extnt[i].right =
                aligned_blocks.pos_extnt[i].right.min(row[i].extents.right);
        }
    }
}

/// Port of NCBI static `_PSIComputeAlignedRegionLengths`.
pub fn psi_compute_aligned_region_lengths(msa: &PsiMsa, aligned_blocks: &mut PsiAlignedBlock) {
    let query_length = msa.dimensions.query_length as usize;
    for i in 0..query_length {
        aligned_blocks.size[i] =
            (aligned_blocks.pos_extnt[i].right - aligned_blocks.pos_extnt[i].left + 1) as u32;
    }

    let x_residue = crate::encoding::AMINOACID_TO_NCBISTDAA[b'X' as usize];
    for i in 0..query_length {
        if msa.query.get(i).copied() == Some(x_residue) {
            for idx in 0..i {
                if aligned_blocks.pos_extnt[idx].right >= i as i32
                    && msa.query.get(idx).copied() != Some(x_residue)
                {
                    aligned_blocks.size[idx] = aligned_blocks.size[idx].saturating_sub(1);
                }
            }
            for idx in ((i + 1)..query_length).rev() {
                if aligned_blocks.pos_extnt[idx].left <= i as i32
                    && msa.query.get(idx).copied() != Some(x_residue)
                {
                    aligned_blocks.size[idx] = aligned_blocks.size[idx].saturating_sub(1);
                }
            }
        }
    }
}

/// Port of NCBI `_PSIComputeAlignmentBlocks`.
pub fn psi_compute_alignment_blocks(
    msa: Option<&mut PsiMsa>,
    aligned_blocks: Option<&mut PsiAlignedBlock>,
) -> i32 {
    let (Some(msa), Some(aligned_blocks)) = (msa, aligned_blocks) else {
        return PSIERR_BADPARAM;
    };
    let q = msa.dimensions.query_length as usize;
    let rows = msa.dimensions.num_seqs.saturating_add(1) as usize;
    if msa.cell.len() < rows
        || !msa.cell.iter().take(rows).all(|row| row.len() >= q)
        || msa.query.len() < q
        || aligned_blocks.pos_extnt.len() < q
        || aligned_blocks.size.len() < q
    {
        return PSIERR_BADPARAM;
    }
    for s in (K_QUERY_INDEX + 1)..msa.dimensions.num_seqs.saturating_add(1) as usize {
        psi_get_left_extents(msa, s as u32);
        psi_get_right_extents(msa, s as u32);
        psi_compute_position_extents(msa, s as u32, aligned_blocks);
    }
    psi_compute_aligned_region_lengths(msa, aligned_blocks);
    PSI_SUCCESS
}

pub fn psi_get_aligned_sequences_for_position(msa: &PsiMsa, position: u32) -> Vec<u32> {
    let position = position as usize;
    let rows = msa.dimensions.num_seqs.saturating_add(1) as usize;
    let mut aligned_sequences = Vec::new();
    for i in 0..rows {
        if msa
            .cell
            .get(i)
            .and_then(|row| row.get(position))
            .is_some_and(|cell| cell.is_aligned)
        {
            aligned_sequences.push(i as u32);
        }
    }
    aligned_sequences
}

/// Port of NCBI static `_PSICalculateNormalizedSequenceWeights`.
pub fn psi_calculate_normalized_sequence_weights(
    msa: &PsiMsa,
    aligned_blocks: &PsiAlignedBlock,
    position: u32,
    aligned_seqs: &[u32],
    seq_weights: &mut PsiSequenceWeights,
) {
    let position = position as usize;
    if aligned_seqs.is_empty() {
        return;
    }
    let gap_residue = crate::encoding::AMINOACID_TO_NCBISTDAA[b'-' as usize] as usize;
    let x_residue = crate::encoding::AMINOACID_TO_NCBISTDAA[b'X' as usize] as usize;
    let left = aligned_blocks.pos_extnt[position].left.max(0) as usize;
    let right = aligned_blocks.pos_extnt[position]
        .right
        .min(msa.dimensions.query_length.saturating_sub(1) as i32) as usize;
    if left > right {
        return;
    }

    let mut distinct_residues_found = false;
    let mut sigma = 0u32;
    for i in left..=right {
        let mut residue_counts_for_column = vec![0u32; msa.alphabet_size as usize];
        let mut num_distinct_residues_for_column = 0u32;
        let mut num_local_std_letters = 0usize;

        for &seq_idx in aligned_seqs {
            let residue = msa.cell[seq_idx as usize][i].letter as usize;
            if residue >= residue_counts_for_column.len() {
                continue;
            }
            if residue_counts_for_column[residue] == 0 {
                num_distinct_residues_for_column += 1;
                if residue != gap_residue && residue != x_residue {
                    num_local_std_letters += 1;
                }
            }
            residue_counts_for_column[residue] += 1;
        }

        sigma += num_distinct_residues_for_column;
        let distrib_index = num_local_std_letters.min(EFFECTIVE_ALPHABET);
        seq_weights.pos_distinct_distrib[position][distrib_index] += 1;
        if num_distinct_residues_for_column > 1 {
            distinct_residues_found = true;
        }

        for &seq_idx in aligned_seqs {
            let seq_idx = seq_idx as usize;
            let residue = msa.cell[seq_idx][i].letter as usize;
            if residue < residue_counts_for_column.len() && num_distinct_residues_for_column > 0 {
                seq_weights.row_sigma[seq_idx] += 1.0
                    / (residue_counts_for_column[residue] as f64
                        * num_distinct_residues_for_column as f64);
            }
        }
    }

    seq_weights.sigma[position] = sigma as f64;
    if distinct_residues_found {
        let span = (right - left + 1) as f64;
        let mut weight_sum = 0.0;
        for &seq_idx in aligned_seqs {
            let seq_idx = seq_idx as usize;
            seq_weights.norm_seq_weights[seq_idx] = seq_weights.row_sigma[seq_idx] / span;
            weight_sum += seq_weights.norm_seq_weights[seq_idx];
        }
        if weight_sum > 0.0 {
            for &seq_idx in aligned_seqs {
                let seq_idx = seq_idx as usize;
                seq_weights.norm_seq_weights[seq_idx] /= weight_sum;
            }
        }
    } else {
        let weight = 1.0 / aligned_seqs.len() as f64;
        for &seq_idx in aligned_seqs {
            seq_weights.norm_seq_weights[seq_idx as usize] = weight;
        }
    }
}

/// Port of NCBI static `_PSICalculateMatchWeights`.
pub fn psi_calculate_match_weights(
    msa: &PsiMsa,
    position: u32,
    aligned_seqs: &[u32],
    seq_weights: &mut PsiSequenceWeights,
) {
    let position = position as usize;
    let gap_residue = crate::encoding::AMINOACID_TO_NCBISTDAA[b'-' as usize];
    for &seq_idx in aligned_seqs {
        let seq_idx = seq_idx as usize;
        let residue = msa.cell[seq_idx][position].letter as usize;
        if residue < seq_weights.match_weights[position].len() {
            seq_weights.match_weights[position][residue] += seq_weights.norm_seq_weights[seq_idx];
        }
        if residue != gap_residue as usize {
            seq_weights.gapless_column_weights[position] += seq_weights.norm_seq_weights[seq_idx];
        }
    }
}

/// Port of NCBI static `_PSISpreadGapWeights`.
pub fn psi_spread_gap_weights(
    msa: &PsiMsa,
    seq_weights: &mut PsiSequenceWeights,
    nsg_compatibility_mode: bool,
) {
    let gap_residue = crate::encoding::AMINOACID_TO_NCBISTDAA[b'-' as usize] as usize;
    let x_residue = crate::encoding::AMINOACID_TO_NCBISTDAA[b'X' as usize];
    let expected_num_matching_seqs = if nsg_compatibility_mode { 0 } else { 1 };
    for pos in 0..msa.dimensions.query_length as usize {
        if msa.num_matching_seqs[pos] <= expected_num_matching_seqs
            || msa.cell[K_QUERY_INDEX][pos].letter == x_residue
        {
            continue;
        }
        let gap_weight = seq_weights.match_weights[pos]
            .get(gap_residue)
            .copied()
            .unwrap_or(0.0);
        for res in 0..msa.alphabet_size as usize {
            if seq_weights.std_prob.get(res).copied().unwrap_or(0.0) > POS_EPSILON {
                seq_weights.match_weights[pos][res] += gap_weight * seq_weights.std_prob[res];
            }
        }
        if gap_residue < seq_weights.match_weights[pos].len() {
            seq_weights.match_weights[pos][gap_residue] = 0.0;
        }
    }
}

/// Port of NCBI static `_PSICheckSequenceWeights`
/// (`blast_psi_priv.c:1977`).
pub fn psi_check_sequence_weights(
    msa: Option<&PsiMsa>,
    seq_weights: Option<&PsiSequenceWeights>,
    nsg_compatibility_mode: bool,
) -> i32 {
    let (Some(msa), Some(seq_weights)) = (msa, seq_weights) else {
        return PSIERR_BADPARAM;
    };
    let x_residue = crate::encoding::AMINOACID_TO_NCBISTDAA[b'X' as usize];
    let expected_num_matching_seqs = if nsg_compatibility_mode { 0 } else { 1 };

    for pos in 0..msa.dimensions.query_length as usize {
        let num_matching = msa.num_matching_seqs.get(pos).copied().unwrap_or(0);
        let query_letter = msa
            .cell
            .get(K_QUERY_INDEX)
            .and_then(|row| row.get(pos))
            .map(|cell| cell.letter)
            .or_else(|| msa.query.get(pos).copied())
            .unwrap_or(0);
        if num_matching <= expected_num_matching_seqs || query_letter == x_residue {
            continue;
        }

        let Some(weights) = seq_weights.match_weights.get(pos) else {
            return PSIERR_BADSEQWEIGHTS;
        };
        let mut running_total = 0.0;
        for residue in 0..msa.alphabet_size as usize {
            running_total += weights.get(residue).copied().unwrap_or(0.0);
        }
        if !(0.99..=1.01).contains(&running_total) {
            return PSIERR_BADSEQWEIGHTS;
        }
    }
    0
}

/// Rust cleanup equivalent of C `s_PSIComputeFrequenciesFromCDsCleanup`.
pub fn s_psi_compute_frequencies_from_cds_cleanup(sum_weights: &mut Vec<f64>) {
    sum_weights.clear();
}

/// blast-rs spelling of the PSI CDs frequency routine; `CDs` is kept as one acronym.
pub fn psi_compute_frequencies_from_cds(
    cd_msa: Option<&PSICdMsa>,
    sbp: Option<&crate::stat::BlastScoreBlk>,
    options: Option<&crate::options::PSIBlastOptions>,
    seq_weights: Option<&mut PsiSequenceWeights>,
) -> i32 {
    let (Some(cd_msa), Some(sbp), Some(_options), Some(seq_weights)) =
        (cd_msa, sbp, options, seq_weights)
    else {
        return PSIERR_BADPARAM;
    };
    let Some(dimensions) = cd_msa.dimensions else {
        return PSIERR_BADPARAM;
    };
    if dimensions.num_seqs == 0 {
        return PSI_SUCCESS;
    }

    let alphabet_size = sbp.alphabet_size;
    let q = dimensions.query_length as usize;
    let rows = dimensions.num_seqs as usize;
    if alphabet_size > AA_SIZE
        || cd_msa.query.len() < q
        || cd_msa.msa.len() < rows
        || !cd_msa.msa.iter().take(rows).all(|row| row.len() >= q)
        || seq_weights.independent_observations.len() < q
        || !psi_matrix_has_shape(&seq_weights.match_weights, q, alphabet_size)
    {
        return PSIERR_BADPARAM;
    }
    let x_residue = crate::encoding::AMINOACID_TO_NCBISTDAA[b'X' as usize] as usize;
    let mut sum_weights = vec![0.0; alphabet_size];
    for pos in 0..q {
        let mut total_observations = 0.0;
        sum_weights.fill(0.0);

        for msa_index in 0..rows {
            let cell = &cd_msa.msa[msa_index][pos];
            if !cell.is_aligned {
                continue;
            }
            let Some(data) = cell.data.as_ref() else {
                return PSIERR_BADPROFILE;
            };
            if data.wfreqs.len() < alphabet_size {
                return PSIERR_BADPROFILE;
            }
            total_observations += data.iobsr;
            for (residue, value) in sum_weights.iter_mut().enumerate().take(alphabet_size) {
                *value += data.wfreqs[residue] * data.iobsr;
            }
        }

        let query_residue = cd_msa.query[pos] as usize;
        if total_observations > 0.0
            && query_residue != x_residue
            && query_residue < alphabet_size
            && sum_weights[query_residue] == 0.0
        {
            sum_weights[query_residue] = 1.0;
            total_observations += 1.0;
        }

        if total_observations > 0.0 {
            for (residue, weight) in sum_weights.iter().enumerate().take(alphabet_size) {
                seq_weights.match_weights[pos][residue] = *weight / total_observations;
            }
        }
        seq_weights.independent_observations[pos] =
            (MAX_IND_OBSERVATIONS as f64).min(total_observations);
    }

    s_psi_compute_frequencies_from_cds_cleanup(&mut sum_weights);
    PSI_SUCCESS
}

fn psi_scale_pssm_by_factor(internal_pssm: &mut PsiInternalPssmData, factor: f64) {
    let rows = internal_pssm.ncols as usize;
    let cols = internal_pssm.nrows as usize;
    for i in 0..rows {
        for j in 0..cols {
            let scaled = internal_pssm.scaled_pssm[i][j];
            internal_pssm.pssm[i][j] = if scaled != crate::stat::BLAST_SCORE_MIN {
                crate::math::blast_nint(factor * scaled as f64 / K_PSI_SCALE_FACTOR) as i32
            } else {
                crate::stat::BLAST_SCORE_MIN
            };
        }
    }
}

/// Port of NCBI `_PSIScaleMatrix`.
pub fn psi_scale_matrix(
    query: Option<&[u8]>,
    std_probs: Option<&[f64]>,
    internal_pssm: Option<&mut PsiInternalPssmData>,
    sbp: Option<&mut crate::stat::BlastScoreBlk>,
) -> i32 {
    let (Some(query), Some(std_probs), Some(internal_pssm), Some(sbp)) =
        (query, std_probs, internal_pssm, sbp)
    else {
        return PSIERR_BADPARAM;
    };
    let Some(ideal) = sbp.kbp_ideal.as_ref() else {
        return PSIERR_BADPARAM;
    };
    let query_length = internal_pssm.ncols as usize;
    let alphabet_size = internal_pssm.nrows as usize;
    if alphabet_size > AA_SIZE
        || query.len() < query_length
        || std_probs.len() < alphabet_size
        || sbp.kbp_psi.is_empty()
        || !psi_matrix_has_shape(&internal_pssm.pssm, query_length, alphabet_size)
        || !psi_matrix_has_shape(&internal_pssm.scaled_pssm, query_length, alphabet_size)
    {
        return PSIERR_BADPARAM;
    }

    let ideal_lambda = ideal.lambda;
    let mut first_time = true;
    let mut factor_low = 1.0;
    let mut factor_high = 1.0;
    let mut factor = 1.0;
    let mut too_high = true;

    loop {
        psi_scale_pssm_by_factor(internal_pssm, factor);
        if crate::blast_kappa::psi_update_lambda_k(
            &internal_pssm.pssm,
            query,
            query_length,
            std_probs,
            sbp,
        ) != 0
        {
            return PSIERR_POSITIVEAVGSCORE;
        }
        let new_lambda = sbp.kbp_psi[0].lambda;
        if new_lambda > ideal_lambda {
            if first_time {
                factor_high = 1.0 + K_POSIT_SCALING_PERCENT;
                factor = factor_high;
                factor_low = 1.0;
                too_high = true;
                first_time = false;
            } else {
                if !too_high {
                    break;
                }
                factor_high += factor_high - 1.0;
                factor = factor_high;
            }
        } else if new_lambda > 0.0 {
            if first_time {
                factor_high = 1.0;
                factor_low = 1.0 - K_POSIT_SCALING_PERCENT;
                factor = factor_low;
                too_high = false;
                first_time = false;
            } else {
                if too_high {
                    break;
                }
                factor_low += factor_low - 1.0;
                factor = factor_low;
            }
        } else {
            return PSIERR_POSITIVEAVGSCORE;
        }
    }

    for _ in 0..K_POSIT_SCALING_NUM_ITERATIONS {
        factor = (factor_high + factor_low) / 2.0;
        psi_scale_pssm_by_factor(internal_pssm, factor);
        if crate::blast_kappa::psi_update_lambda_k(
            &internal_pssm.pssm,
            query,
            query_length,
            std_probs,
            sbp,
        ) != 0
        {
            return PSIERR_POSITIVEAVGSCORE;
        }
        if sbp.kbp_psi[0].lambda > ideal_lambda {
            factor_low = factor;
        } else {
            factor_high = factor;
        }
    }

    PSI_SUCCESS
}

/// blast-rs spelling of the PSI CDs frequency-ratio routine; `CDs` is kept as one acronym.
pub fn psi_compute_freq_ratios_from_cds(
    cd_msa: Option<&PSICdMsa>,
    seq_weights: Option<&PsiSequenceWeights>,
    sbp: Option<&crate::stat::BlastScoreBlk>,
    pseudo_count: i32,
    internal_pssm: Option<&mut PsiInternalPssmData>,
) -> i32 {
    let (Some(cd_msa), Some(seq_weights), Some(sbp), Some(internal_pssm)) =
        (cd_msa, seq_weights, sbp, internal_pssm)
    else {
        return PSIERR_BADPARAM;
    };
    let Some(dimensions) = cd_msa.dimensions else {
        return PSIERR_BADPARAM;
    };
    if pseudo_count < 0 {
        return PSIERR_BADPARAM;
    }
    let q = dimensions.query_length as usize;
    let a = sbp.alphabet_size;
    if a > AA_SIZE
        || cd_msa.query.len() < q
        || seq_weights.independent_observations.len() < q
        || seq_weights.std_prob.len() < a
        || !psi_matrix_has_shape(&seq_weights.match_weights, q, a)
        || !psi_matrix_has_shape(&internal_pssm.freq_ratios, q, a)
    {
        return PSIERR_BADPARAM;
    }
    let Some(freq_ratios) =
        crate::matrix::get_matrix_freq_ratios(sbp.name.as_deref().unwrap_or("BLOSUM62"))
    else {
        return PSIERR_OUTOFMEM;
    };

    let x_residue = crate::encoding::AMINOACID_TO_NCBISTDAA[b'X' as usize];
    let mut std_prob_array = [0.0; AA_SIZE];
    for (i, value) in std_prob_array.iter_mut().enumerate() {
        *value = seq_weights.std_prob.get(i).copied().unwrap_or(0.0);
    }

    for p in 0..q {
        let mut observations = 0.0;
        let mut column_counts = 0.0;
        if cd_msa.query.get(p).copied() != Some(x_residue) {
            observations = (seq_weights.independent_observations[p] - 1.0).max(0.0);
            if pseudo_count == 0 {
                let mut match_weights = [0.0; AA_SIZE];
                for (r, value) in match_weights.iter_mut().enumerate() {
                    *value = seq_weights.match_weights[p].get(r).copied().unwrap_or(0.0);
                }
                column_counts =
                    compute_column_pseudocounts(&match_weights, &std_prob_array, observations);
            } else {
                column_counts = pseudo_count as f64;
            }
        }
        let pseudo_weight = if column_counts >= PSEUDO_MAX {
            observations = 0.0;
            ZERO_OBS_PSEUDO
        } else {
            column_counts
        };

        for r in 0..a {
            if cd_msa.query.get(p).copied() == Some(x_residue)
                || seq_weights.std_prob.get(r).copied().unwrap_or(0.0) <= POS_EPSILON
            {
                internal_pssm.freq_ratios[p][r] = 0.0;
            } else {
                let mut pseudo = 0.0;
                for i in 0..a {
                    if sbp
                        .matrix
                        .data
                        .get(r)
                        .and_then(|row| row.get(i))
                        .copied()
                        .unwrap_or(crate::stat::BLAST_SCORE_MIN)
                        != crate::stat::BLAST_SCORE_MIN
                    {
                        pseudo += seq_weights.match_weights[p][i] * freq_ratios[r][i];
                    }
                }
                pseudo *= pseudo_weight;
                let numerator = observations * seq_weights.match_weights[p][r]
                    / seq_weights.std_prob[r]
                    + pseudo;
                let denominator = observations + pseudo_weight;
                let q_over_p_estimate = if denominator == 0.0 {
                    0.0
                } else {
                    numerator / denominator
                };
                internal_pssm.freq_ratios[p][r] = q_over_p_estimate * seq_weights.std_prob[r];
            }
        }
    }

    PSI_SUCCESS
}

/// Rust ownership equivalent of NCBI `_PSIInternalPssmDataFree`.
pub fn psi_internal_pssm_data_free(
    mut pssm: Option<PsiInternalPssmData>,
) -> Option<PsiInternalPssmData> {
    if let Some(pssm) = pssm.as_mut() {
        for row in &mut pssm.pssm {
            row.clear();
        }
        pssm.pssm.clear();
        for row in &mut pssm.scaled_pssm {
            row.clear();
        }
        pssm.scaled_pssm.clear();
        for row in &mut pssm.freq_ratios {
            row.clear();
        }
        pssm.freq_ratios.clear();
        pssm.pseudocounts.clear();
        pssm.ncols = 0;
        pssm.nrows = 0;
    }
    None
}

/// Port of NCBI `_PSIInternalPssmDataNew`.
pub fn psi_internal_pssm_data_new(
    query_length: u32,
    alphabet_size: u32,
) -> Option<PsiInternalPssmData> {
    Some(PsiInternalPssmData {
        ncols: query_length,
        nrows: alphabet_size,
        pssm: vec![vec![0; alphabet_size as usize]; query_length as usize],
        scaled_pssm: vec![vec![0; alphabet_size as usize]; query_length as usize],
        freq_ratios: vec![vec![0.0; alphabet_size as usize]; query_length as usize],
        pseudocounts: vec![0.0; query_length as usize],
    })
}

/// Port of NCBI `_PSIAlignedBlockNew`.
pub fn psi_aligned_block_new(query_length: u32) -> Option<PsiAlignedBlock> {
    Some(PsiAlignedBlock {
        pos_extnt: vec![
            SSeqRange {
                left: -1,
                right: query_length as i32
            };
            query_length as usize
        ],
        size: vec![0; query_length as usize],
    })
}

/// Port of NCBI `_PSISequenceWeightsNew`.
pub fn psi_sequence_weights_new(
    dimensions: Option<&PSIMsaDimensions>,
    alphabet_size: usize,
) -> Option<PsiSequenceWeights> {
    let dimensions = *dimensions?;
    let rows = dimensions.num_seqs.saturating_add(1) as usize;
    let cols = dimensions.query_length as usize;
    Some(PsiSequenceWeights {
        match_weights: vec![vec![0.0; alphabet_size]; cols],
        match_weights_size: dimensions.query_length,
        norm_seq_weights: vec![0.0; rows],
        row_sigma: vec![0.0; rows],
        sigma: vec![0.0; cols],
        std_prob: crate::stat::protein_std_freq_ncbistdaa().to_vec(),
        gapless_column_weights: vec![0.0; cols],
        pos_distinct_distrib: vec![vec![0; EFFECTIVE_ALPHABET + 1]; cols + 1],
        pos_distinct_distrib_size: dimensions.query_length + 1,
        pos_num_participating: vec![0; cols + 1],
        independent_observations: vec![0.0; cols + 1],
    })
}

/// Rust ownership equivalent of NCBI `_PSISequenceWeightsFree`.
pub fn psi_sequence_weights_free(
    mut seq_weights: Option<PsiSequenceWeights>,
) -> Option<PsiSequenceWeights> {
    if let Some(seq_weights) = seq_weights.as_mut() {
        for row in &mut seq_weights.match_weights {
            row.clear();
        }
        seq_weights.match_weights.clear();
        seq_weights.norm_seq_weights.clear();
        seq_weights.row_sigma.clear();
        seq_weights.sigma.clear();
        seq_weights.std_prob.clear();
        seq_weights.gapless_column_weights.clear();
        for row in &mut seq_weights.pos_distinct_distrib {
            row.clear();
        }
        seq_weights.pos_distinct_distrib.clear();
        seq_weights.pos_num_participating.clear();
        seq_weights.independent_observations.clear();
        seq_weights.match_weights_size = 0;
        seq_weights.pos_distinct_distrib_size = 0;
    }
    None
}

/// Port of NCBI static `s_fillColumnProbabilities`.
pub fn s_fill_column_probabilities(
    pos_search: &PsiSequenceWeights,
    column_number: usize,
) -> Option<[f64; EFFECTIVE_ALPHABET]> {
    let row = pos_search.match_weights.get(column_number)?;
    let mut probabilities = [0.0; EFFECTIVE_ALPHABET];
    for (c, &residue) in CHAR_ORDER.iter().enumerate() {
        probabilities[c] = row.get(residue as usize).copied().unwrap_or(0.0);
    }
    Some(probabilities)
}

/// Port of NCBI `_PSIComputeSequenceWeights`.
pub fn psi_compute_sequence_weights(
    msa: Option<&PsiMsa>,
    aligned_blocks: Option<&PsiAlignedBlock>,
    nsg_compatibility_mode: bool,
    seq_weights: Option<&mut PsiSequenceWeights>,
) -> i32 {
    let (Some(msa), Some(aligned_blocks), Some(seq_weights)) = (msa, aligned_blocks, seq_weights)
    else {
        return PSIERR_BADPARAM;
    };
    let q = msa.dimensions.query_length as usize;
    let rows = msa.dimensions.num_seqs.saturating_add(1) as usize;
    let a = msa.alphabet_size as usize;
    if a > AA_SIZE
        || msa.cell.len() < rows
        || !msa.cell.iter().take(rows).all(|row| row.len() >= q)
        || msa.num_matching_seqs.len() < q
        || aligned_blocks.pos_extnt.len() < q
        || aligned_blocks.size.len() < q
        || seq_weights.norm_seq_weights.len() < rows
        || seq_weights.row_sigma.len() < rows
        || seq_weights.sigma.len() < q
        || seq_weights.gapless_column_weights.len() < q
        || seq_weights.std_prob.len() < a
        || seq_weights.pos_num_participating.len() < q
        || seq_weights.pos_distinct_distrib.len() < q
        || !seq_weights
            .pos_distinct_distrib
            .iter()
            .take(q)
            .all(|row| row.len() > EFFECTIVE_ALPHABET)
        || !psi_matrix_has_shape(&seq_weights.match_weights, q, a)
    {
        return PSIERR_BADPARAM;
    }
    let expected_num_matching = if nsg_compatibility_mode { 0 } else { 1 };

    for pos in 0..q {
        if aligned_blocks.size.get(pos).copied().unwrap_or(0) == 0
            || msa.num_matching_seqs.get(pos).copied().unwrap_or(0) <= expected_num_matching
        {
            continue;
        }
        let aligned_seqs = psi_get_aligned_sequences_for_position(msa, pos as u32);
        if aligned_seqs.len() as u32 <= expected_num_matching {
            continue;
        }

        seq_weights.norm_seq_weights.fill(0.0);
        seq_weights.row_sigma.fill(0.0);
        psi_calculate_normalized_sequence_weights(
            msa,
            aligned_blocks,
            pos as u32,
            &aligned_seqs,
            seq_weights,
        );
        if let Some(slot) = seq_weights.pos_num_participating.get_mut(pos) {
            *slot = aligned_seqs.len() as i32;
        }
        psi_calculate_match_weights(msa, pos as u32, &aligned_seqs, seq_weights);
    }

    let mut status =
        psi_check_sequence_weights(Some(msa), Some(seq_weights), nsg_compatibility_mode);
    if status != PSI_SUCCESS {
        return status;
    }
    psi_spread_gap_weights(msa, seq_weights, nsg_compatibility_mode);
    status = psi_check_sequence_weights(Some(msa), Some(seq_weights), nsg_compatibility_mode);
    status
}

fn psi_effective_observations_from_blocks(
    aligned_blocks: &PsiAlignedBlock,
    seq_weights: &PsiSequenceWeights,
    column_number: usize,
    query_length: usize,
    expno: &[f64; MAX_IND_OBSERVATIONS + 1],
) -> f64 {
    let Some(extent) = aligned_blocks.pos_extnt.get(column_number).copied() else {
        return 0.0;
    };
    if extent.left < 0 || extent.right >= query_length as i32 {
        return 0.0;
    }
    let half_num_columns = ((extent.right - extent.left + 2) / 2).max(1) as i32;
    let Some(distrib) = seq_weights.pos_distinct_distrib.get(column_number) else {
        return 0.0;
    };
    let mut columns_accounted_for = 0i32;
    let mut total_distinct_counts = 0i32;
    for k in (0..=EFFECTIVE_ALPHABET).rev() {
        let count = distrib.get(k).copied().unwrap_or(0);
        total_distinct_counts += count * k as i32;
        columns_accounted_for += count;
        if columns_accounted_for > half_num_columns {
            total_distinct_counts -= (columns_accounted_for - half_num_columns) * k as i32;
            columns_accounted_for = half_num_columns;
        }
        if columns_accounted_for >= half_num_columns {
            break;
        }
    }
    if columns_accounted_for == 0 {
        return 0.0;
    }
    let ave_distinct = total_distinct_counts as f64 / columns_accounted_for as f64;
    let mut i = 1usize;
    while i < MAX_IND_OBSERVATIONS && expno[i] <= ave_distinct {
        i += 1;
    }
    let indep = if i == MAX_IND_OBSERVATIONS {
        i as f64
    } else {
        i as f64 - (expno[i] - ave_distinct) / (expno[i] - expno[i - 1])
    };
    let participating = seq_weights
        .pos_num_participating
        .get(column_number)
        .copied()
        .unwrap_or(0) as f64;
    (indep.min(participating) - 1.0).max(0.0)
}

/// Port of NCBI `_PSIComputeFreqRatios`.
pub fn psi_compute_freq_ratios(
    msa: Option<&PsiMsa>,
    seq_weights: Option<&PsiSequenceWeights>,
    sbp: Option<&crate::stat::BlastScoreBlk>,
    aligned_blocks: Option<&PsiAlignedBlock>,
    pseudo_count: i32,
    nsg_compatibility_mode: bool,
    internal_pssm: Option<&mut PsiInternalPssmData>,
) -> i32 {
    let (Some(msa), Some(seq_weights), Some(sbp), Some(aligned_blocks), Some(internal_pssm)) =
        (msa, seq_weights, sbp, aligned_blocks, internal_pssm)
    else {
        return PSIERR_BADPARAM;
    };
    if pseudo_count < 0 || sbp.alphabet_size != msa.alphabet_size as usize {
        return PSIERR_BADPARAM;
    }
    let q = msa.dimensions.query_length as usize;
    let a = msa.alphabet_size as usize;
    if a > AA_SIZE
        || msa.cell.len() <= K_QUERY_INDEX
        || msa.cell[K_QUERY_INDEX].len() < q
        || msa.num_matching_seqs.len() < q
        || aligned_blocks.pos_extnt.len() < q
        || aligned_blocks.size.len() < q
        || seq_weights.std_prob.len() < a
        || seq_weights.pos_num_participating.len() < q
        || seq_weights.pos_distinct_distrib.len() < q
        || !seq_weights
            .pos_distinct_distrib
            .iter()
            .take(q)
            .all(|row| row.len() > EFFECTIVE_ALPHABET)
        || !psi_matrix_has_shape(&seq_weights.match_weights, q, a)
        || !psi_matrix_has_shape(&internal_pssm.freq_ratios, q, a)
        || internal_pssm.pseudocounts.len() < q
    {
        return PSIERR_BADPARAM;
    }
    let Some(freq_ratios) =
        crate::matrix::get_matrix_freq_ratios(sbp.name.as_deref().unwrap_or("BLOSUM62"))
    else {
        return PSIERR_OUTOFMEM;
    };
    let std_prob_array = crate::stat::protein_std_freq_ncbistdaa();
    let expno = s_initialize_exp_num_observations(&effective_bg_probs(&std_prob_array));
    let x_residue = crate::encoding::NCBISTDAA_X;

    for p in 0..q {
        let mut column_counts = 0.0;
        let mut observations = 0.0;
        if msa.cell[K_QUERY_INDEX][p].letter != x_residue {
            observations =
                psi_effective_observations_from_blocks(aligned_blocks, seq_weights, p, q, &expno);
            if pseudo_count == 0 {
                let mut match_weights = [0.0; AA_SIZE];
                for (r, value) in match_weights.iter_mut().enumerate() {
                    *value = seq_weights.match_weights[p].get(r).copied().unwrap_or(0.0);
                }
                column_counts =
                    compute_column_pseudocounts(&match_weights, &std_prob_array, observations);
            } else {
                column_counts = pseudo_count as f64;
            }
        }
        let pseudo_weight = if column_counts >= PSEUDO_MAX {
            observations = 0.0;
            ZERO_OBS_PSEUDO
        } else {
            column_counts
        };

        for r in 0..msa.alphabet_size as usize {
            if msa.cell[K_QUERY_INDEX][p].letter == x_residue
                || seq_weights.std_prob.get(r).copied().unwrap_or(0.0) <= POS_EPSILON
            {
                internal_pssm.freq_ratios[p][r] = 0.0;
            } else {
                let mut pseudo = 0.0;
                for i in 0..msa.alphabet_size as usize {
                    if sbp
                        .matrix
                        .data
                        .get(r)
                        .and_then(|row| row.get(i))
                        .copied()
                        .unwrap_or(crate::stat::BLAST_SCORE_MIN)
                        != crate::stat::BLAST_SCORE_MIN
                    {
                        pseudo += seq_weights.match_weights[p][i] * freq_ratios[r][i];
                    }
                }
                pseudo *= pseudo_weight;
                let denominator = observations + pseudo_weight;
                if nsg_compatibility_mode && denominator == 0.0 {
                    return PSIERR_UNKNOWN;
                }
                if denominator == 0.0 {
                    return PSIERR_BADPARAM;
                }
                let numerator = observations * seq_weights.match_weights[p][r]
                    / seq_weights.std_prob[r]
                    + pseudo;
                internal_pssm.freq_ratios[p][r] =
                    (numerator / denominator) * seq_weights.std_prob[r];
                if let Some(slot) = internal_pssm.pseudocounts.get_mut(p) {
                    *slot = pseudo_weight;
                }
            }
        }
    }

    PSI_SUCCESS
}

/// Port of NCBI `_PSIConvertFreqRatiosToPSSM`.
pub fn psi_convert_freq_ratios_to_pssm(
    internal_pssm: Option<&mut PsiInternalPssmData>,
    query: Option<&[u8]>,
    sbp: Option<&crate::stat::BlastScoreBlk>,
    std_probs: Option<&[f64]>,
) -> i32 {
    let (Some(internal_pssm), Some(query), Some(sbp), Some(std_probs)) =
        (internal_pssm, query, sbp, std_probs)
    else {
        return PSIERR_BADPARAM;
    };
    let Some(ideal) = sbp.kbp_ideal.as_ref() else {
        return PSIERR_BADPARAM;
    };
    let Some(freq_ratios) =
        crate::matrix::get_matrix_freq_ratios_with_scale(sbp.name.as_deref().unwrap_or("BLOSUM62"))
    else {
        return PSIERR_OUTOFMEM;
    };
    let x_residue = crate::encoding::NCBISTDAA_X as usize;
    let star_residue = crate::encoding::NCBISTDAA_STOP as usize;
    let ncols = internal_pssm.ncols as usize;
    let nrows = internal_pssm.nrows as usize;
    if nrows > AA_SIZE
        || query.len() < ncols
        || std_probs.len() < nrows
        || !psi_matrix_has_shape(&internal_pssm.freq_ratios, ncols, nrows)
        || !psi_matrix_has_shape(&internal_pssm.pssm, ncols, nrows)
        || !psi_matrix_has_shape(&internal_pssm.scaled_pssm, ncols, nrows)
    {
        return PSIERR_BADPARAM;
    }

    for i in 0..ncols {
        let mut is_unaligned_column = true;
        let residue = query[i] as usize;
        let matrix_residue = if residue < AA_SIZE {
            residue
        } else {
            x_residue
        };
        for j in 0..nrows {
            let q_over_p = if std_probs.get(j).copied().unwrap_or(0.0) > POS_EPSILON {
                internal_pssm.freq_ratios[i][j] / std_probs[j]
            } else {
                0.0
            };
            if q_over_p != 0.0 {
                is_unaligned_column = false;
            }
            internal_pssm.scaled_pssm[i][j] = if q_over_p == 0.0
                || std_probs.get(j).copied().unwrap_or(0.0) < POS_EPSILON
            {
                crate::stat::BLAST_SCORE_MIN
            } else {
                crate::math::blast_nint(K_PSI_SCALE_FACTOR * (q_over_p.ln() / ideal.lambda)) as i32
            };
            if (j == x_residue || j == star_residue)
                && sbp
                    .matrix
                    .data
                    .get(matrix_residue)
                    .and_then(|row| row.get(x_residue))
                    .copied()
                    .unwrap_or(crate::stat::BLAST_SCORE_MIN)
                    != crate::stat::BLAST_SCORE_MIN
            {
                internal_pssm.scaled_pssm[i][j] = sbp
                    .matrix
                    .data
                    .get(matrix_residue)
                    .and_then(|row| row.get(j))
                    .copied()
                    .unwrap_or(crate::stat::BLAST_SCORE_MIN)
                    * K_PSI_SCALE_FACTOR as i32;
            }
        }

        if is_unaligned_column {
            for j in 0..nrows {
                internal_pssm.pssm[i][j] = sbp
                    .matrix
                    .data
                    .get(matrix_residue)
                    .and_then(|row| row.get(j))
                    .copied()
                    .unwrap_or(crate::stat::BLAST_SCORE_MIN);
                internal_pssm.scaled_pssm[i][j] = if freq_ratios.data[matrix_residue][j] != 0.0 {
                    crate::math::blast_nint(
                        K_PSI_SCALE_FACTOR
                            * freq_ratios.bit_scale_factor as f64
                            * freq_ratios.data[matrix_residue][j].ln()
                            / crate::math::NCBIMATH_LN2,
                    ) as i32
                } else {
                    crate::stat::BLAST_SCORE_MIN
                };
            }
        }
    }

    PSI_SUCCESS
}

/// Port of NCBI `_PSICreateAndScalePssmFromFrequencyRatios`.
pub fn psi_create_and_scale_pssm_from_frequency_ratios(
    internal_pssm: Option<&mut PsiInternalPssmData>,
    query: Option<&[u8]>,
    std_prob: Option<&[f64]>,
    sbp: Option<&mut crate::stat::BlastScoreBlk>,
    impala_scaling_factor: f64,
) -> i32 {
    let (Some(internal_pssm), Some(query), Some(std_prob), Some(sbp)) =
        (internal_pssm, query, std_prob, sbp)
    else {
        return PSIERR_BADPARAM;
    };
    let status = psi_convert_freq_ratios_to_pssm(
        Some(internal_pssm),
        Some(query),
        Some(&*sbp),
        Some(std_prob),
    );
    if status != PSI_SUCCESS {
        return status;
    }
    if impala_scaling_factor == crate::options::K_PSSM_NO_IMPALA_SCALING {
        psi_scale_matrix(Some(query), Some(std_prob), Some(internal_pssm), Some(sbp))
    } else {
        let compact = crate::blast_kappa::KappaCompactSearchItems {
            standard_prob: std_prob.to_vec(),
            query: query.to_vec(),
            qlength: internal_pssm.ncols as usize,
            alphabet_size: internal_pssm.nrows as usize,
            matrix: sbp.matrix.data.clone(),
            kbp_std: sbp.kbp_std.clone(),
            kbp_psi: sbp.kbp_psi.clone(),
            kbp_gap_std: sbp.kbp_gap_std.clone(),
            kbp_gap_psi: sbp.kbp_gap_psi.clone(),
            lambda_ideal: sbp.kbp_ideal.as_ref().map(|kbp| kbp.lambda).unwrap_or(0.0),
            k_ideal: sbp.kbp_ideal.as_ref().map(|kbp| kbp.k).unwrap_or(0.0),
        };
        if crate::blast_kappa::impala_scale_matrix(
            &compact,
            &mut internal_pssm.pssm,
            &mut internal_pssm.scaled_pssm,
            impala_scaling_factor,
            true,
            sbp,
        ) {
            PSI_SUCCESS
        } else {
            PSIERR_POSITIVEAVGSCORE
        }
    }
}

/// Port of NCBI `_PSICalculateInformationContentFromFreqRatios`.
pub fn psi_calculate_information_content_from_freq_ratios(
    freq_ratios: &[Vec<f64>],
    std_prob: &[f64],
    query_length: u32,
    alphabet_size: u32,
) -> Option<Vec<f64>> {
    if std_prob.is_empty() {
        return None;
    }

    let q = query_length as usize;
    let a = alphabet_size as usize;
    let mut retval = vec![0.0; q];
    for p in 0..q {
        let row = freq_ratios.get(p)?;
        let mut info_sum = 0.0;
        for r in 0..a {
            let std = *std_prob.get(r)?;
            if std > POS_EPSILON {
                let ratio = *row.get(r)?;
                let q_over_p_estimate = ratio / std;
                if q_over_p_estimate > POS_EPSILON {
                    info_sum += ratio * q_over_p_estimate.ln() / crate::math::NCBIMATH_LN2;
                }
            }
        }
        retval[p] = info_sum;
    }
    Some(retval)
}

/// Port of NCBI `_PSISaveDiagnostics`.
pub fn psi_save_diagnostics(
    msa: Option<&PsiMsa>,
    aligned_block: Option<&PsiAlignedBlock>,
    seq_weights: Option<&PsiSequenceWeights>,
    internal_pssm: Option<&PsiInternalPssmData>,
    diagnostics: Option<&mut PSIDiagnosticsResponse>,
) -> i32 {
    let (Some(msa), Some(aligned_block), Some(seq_weights), Some(internal_pssm), Some(diagnostics)) =
        (msa, aligned_block, seq_weights, internal_pssm, diagnostics)
    else {
        return PSIERR_BADPARAM;
    };
    if internal_pssm.freq_ratios.is_empty() {
        return PSIERR_BADPARAM;
    }

    let q = diagnostics.query_length as usize;
    let a = diagnostics.alphabet_size as usize;
    if msa.dimensions.query_length as usize != q {
        return PSIERR_BADPARAM;
    }
    if !psi_diagnostics_response_buffers_are_large_enough(diagnostics, q, a) {
        return PSIERR_BADPARAM;
    }
    if !psi_save_diagnostics_sources_are_large_enough(
        msa,
        aligned_block,
        seq_weights,
        internal_pssm,
        diagnostics,
        q,
        a,
    ) {
        return PSIERR_BADPARAM;
    }

    if let Some(out) = diagnostics.information_content.as_mut() {
        let Some(info) = psi_calculate_information_content_from_freq_ratios(
            &internal_pssm.freq_ratios,
            &seq_weights.std_prob,
            diagnostics.query_length,
            diagnostics.alphabet_size,
        ) else {
            return PSIERR_BADPARAM;
        };
        out[..q].copy_from_slice(&info[..q]);
    }

    if let Some(out) = diagnostics.residue_freqs.as_mut() {
        for p in 0..q {
            for r in 0..a {
                out[p][r] = msa.residue_counts[p][r];
            }
        }
    }

    if let Some(out) = diagnostics.weighted_residue_freqs.as_mut() {
        for p in 0..q {
            for r in 0..a {
                out[p][r] = seq_weights.match_weights[p][r];
            }
        }
    }

    if let Some(out) = diagnostics.frequency_ratios.as_mut() {
        for p in 0..q {
            for r in 0..a {
                out[p][r] = internal_pssm.freq_ratios[p][r];
            }
        }
    }

    if let Some(out) = diagnostics.gapless_column_weights.as_mut() {
        let x_residue = crate::encoding::AMINOACID_TO_NCBISTDAA[b'X' as usize];
        for p in 0..q {
            if msa.num_matching_seqs[p] > 1 && msa.cell[0][p].letter != x_residue {
                out[p] = seq_weights.gapless_column_weights[p] / internal_pssm.pseudocounts[p];
                out[p] *= seq_weights.sigma[p] / aligned_block.size[p] as f64 - 1.0;
            } else {
                out[p] = 0.0;
            }
        }
    }

    if let Some(out) = diagnostics.sigma.as_mut() {
        out[..q].copy_from_slice(&seq_weights.sigma[..q]);
    }
    if let Some(out) = diagnostics.interval_sizes.as_mut() {
        out[..q].copy_from_slice(&aligned_block.size[..q]);
    }
    if let Some(out) = diagnostics.num_matching_seqs.as_mut() {
        out[..q].copy_from_slice(&msa.num_matching_seqs[..q]);
    }
    if let Some(out) = diagnostics.independent_observations.as_mut() {
        out[..q].copy_from_slice(&seq_weights.independent_observations[..q]);
    }

    PSI_SUCCESS
}

/// Port of NCBI `_PSISaveCDDiagnostics`.
pub fn psi_save_cd_diagnostics(
    cd_msa: Option<&PSICdMsa>,
    seq_weights: Option<&PsiSequenceWeights>,
    internal_pssm: Option<&PsiInternalPssmData>,
    diagnostics: Option<&mut PSIDiagnosticsResponse>,
) -> i32 {
    let (Some(cd_msa), Some(seq_weights), Some(internal_pssm), Some(diagnostics)) =
        (cd_msa, seq_weights, internal_pssm, diagnostics)
    else {
        return PSIERR_BADPARAM;
    };
    if internal_pssm.freq_ratios.is_empty() {
        return PSIERR_BADPARAM;
    }
    let Some(dimensions) = cd_msa.dimensions else {
        return PSIERR_BADPARAM;
    };
    if dimensions.query_length != diagnostics.query_length {
        return PSIERR_BADPARAM;
    }

    let q = diagnostics.query_length as usize;
    let a = diagnostics.alphabet_size as usize;
    let rows = dimensions.num_seqs as usize;
    if cd_msa.query.len() < q
        || cd_msa.msa.len() < rows
        || !cd_msa.msa.iter().take(rows).all(|row| row.len() >= q)
    {
        return PSIERR_BADPARAM;
    }
    if !psi_diagnostics_response_buffers_are_large_enough(diagnostics, q, a) {
        return PSIERR_BADPARAM;
    }
    if !psi_save_cd_diagnostics_sources_are_large_enough(
        seq_weights,
        internal_pssm,
        diagnostics,
        q,
        a,
    ) {
        return PSIERR_BADPARAM;
    }
    if let Some(out) = diagnostics.information_content.as_mut() {
        let Some(info) = psi_calculate_information_content_from_freq_ratios(
            &internal_pssm.freq_ratios,
            &seq_weights.std_prob,
            diagnostics.query_length,
            diagnostics.alphabet_size,
        ) else {
            return PSIERR_BADPARAM;
        };
        out[..q].copy_from_slice(&info[..q]);
    }
    if let Some(out) = diagnostics.weighted_residue_freqs.as_mut() {
        for p in 0..q {
            for r in 0..a {
                out[p][r] = seq_weights.match_weights[p][r];
            }
        }
    }
    if let Some(out) = diagnostics.frequency_ratios.as_mut() {
        for p in 0..q {
            for r in 0..a {
                out[p][r] = internal_pssm.freq_ratios[p][r];
            }
        }
    }
    if let Some(out) = diagnostics.independent_observations.as_mut() {
        out[..q].copy_from_slice(&seq_weights.independent_observations[..q]);
    }

    PSI_SUCCESS
}

/// Rust string-returning equivalent of debug `PrintMsaFP`.
pub fn print_msa_fp(msa: Option<&PSIMsa>) -> String {
    let Some(msa) = msa else {
        return String::new();
    };
    let mut out = String::new();
    for (row_index, row) in msa
        .data
        .iter()
        .take(msa.dimensions.num_seqs.saturating_add(1) as usize)
        .enumerate()
    {
        out.push_str(&format!("{row_index:3}\t"));
        for cell in row.iter().take(msa.dimensions.query_length as usize) {
            if cell.is_aligned {
                let letter = crate::encoding::NCBISTDAA_TO_AMINOACID
                    .get(cell.letter as usize)
                    .copied()
                    .unwrap_or(b'?' as i8) as u8;
                out.push(letter as char);
            } else {
                out.push('.');
            }
        }
        out.push('\n');
    }
    out
}

/// Rust string-returning equivalent of debug `PrintMsa`.
pub fn print_msa(_filename: Option<&str>, msa: Option<&PSIMsa>) -> String {
    print_msa_fp(msa)
}

/// Port of NCBI static `s_PSIValidateNoFlankingGaps`.
pub fn s_psi_validate_no_flanking_gaps(msa: Option<&PsiMsa>) -> i32 {
    let Some(msa) = msa else {
        return PSIERR_BADPARAM;
    };
    let gap = crate::encoding::NCBISTDAA_GAP;
    let cols = msa.dimensions.query_length as usize;
    for row in msa.cell.iter().take(msa.dimensions.num_seqs as usize + 1) {
        if let Some(cell) = row.iter().take(cols).find(|cell| cell.is_aligned) {
            if cell.letter == gap {
                return PSIERR_STARTINGGAP;
            }
        }
        if let Some(cell) = row.iter().take(cols).rev().find(|cell| cell.is_aligned) {
            if cell.letter == gap {
                return PSIERR_ENDINGGAP;
            }
        }
    }
    PSI_SUCCESS
}

/// Port of NCBI static `s_PSIValidateAlignedColumns`.
pub fn s_psi_validate_aligned_columns(msa: Option<&PsiMsa>) -> i32 {
    let Some(msa) = msa else {
        return PSIERR_BADPARAM;
    };
    let gap = crate::encoding::NCBISTDAA_GAP;
    let rows = msa.dimensions.num_seqs as usize + 1;
    for position in 0..msa.dimensions.query_length as usize {
        let mut found_aligned_sequence = false;
        let mut found_non_gap_residue = false;
        for row in msa.cell.iter().take(rows) {
            if let Some(cell) = row.get(position) {
                if cell.is_aligned {
                    found_aligned_sequence = true;
                    if cell.letter != gap {
                        found_non_gap_residue = true;
                        break;
                    }
                }
            }
        }
        if !found_aligned_sequence {
            return PSIERR_UNALIGNEDCOLUMN;
        }
        if !found_non_gap_residue {
            return PSIERR_COLUMNOFGAPS;
        }
    }
    PSI_SUCCESS
}

/// Port of NCBI static `s_PSIValidateParticipatingSequences`.
pub fn s_psi_validate_participating_sequences(msa: Option<&PsiMsa>) -> i32 {
    let Some(msa) = msa else {
        return PSIERR_BADPARAM;
    };
    if msa.dimensions.num_seqs > 0 {
        PSI_SUCCESS
    } else {
        PSIERR_NOALIGNEDSEQS
    }
}

fn psi_validate_no_gaps_in_query(msa: &PsiMsa) -> i32 {
    let gap = crate::encoding::NCBISTDAA_GAP;
    let Some(query_row) = msa.cell.first() else {
        return PSIERR_BADPARAM;
    };
    if msa
        .query
        .iter()
        .take(msa.dimensions.query_length as usize)
        .any(|&letter| letter == gap)
        || query_row
            .iter()
            .take(msa.dimensions.query_length as usize)
            .any(|cell| cell.is_aligned && cell.letter == gap)
    {
        PSIERR_GAPINQUERY
    } else {
        PSI_SUCCESS
    }
}

fn psi_msa_has_validation_shape(msa: &PsiMsa) -> bool {
    let q = msa.dimensions.query_length as usize;
    let rows = msa.dimensions.num_seqs.saturating_add(1) as usize;
    msa.cell.len() >= rows
        && msa.cell.iter().take(rows).all(|row| row.len() >= q)
        && msa.query.len() >= q
}

/// Port of NCBI `_PSIValidateMSA`.
pub fn psi_validate_msa(msa: Option<&PsiMsa>, ignore_unaligned_positions: bool) -> i32 {
    let Some(msa) = msa else {
        return PSIERR_BADPARAM;
    };
    if !psi_msa_has_validation_shape(msa) {
        return PSIERR_BADPARAM;
    }
    let mut retval = s_psi_validate_no_flanking_gaps(Some(msa));
    if retval != PSI_SUCCESS {
        return retval;
    }
    if !ignore_unaligned_positions {
        retval = s_psi_validate_aligned_columns(Some(msa));
        if retval != PSI_SUCCESS {
            return retval;
        }
    }
    retval = psi_validate_no_gaps_in_query(msa);
    if retval != PSI_SUCCESS {
        return retval;
    }
    s_psi_validate_participating_sequences(Some(msa))
}

/// Port of NCBI `_PSIValidateMSA_StructureGroup`.
pub fn psi_validate_msa_structure_group(msa: Option<&PsiMsa>) -> i32 {
    let Some(msa) = msa else {
        return PSIERR_BADPARAM;
    };
    if !psi_msa_has_validation_shape(msa) {
        return PSIERR_BADPARAM;
    }
    s_psi_validate_participating_sequences(Some(msa))
}

/// Port of NCBI `_PSIValidateCdMSA`.
pub fn psi_validate_cd_msa(cd_msa: Option<&PSICdMsa>, alphabet_size: u32) -> i32 {
    let Some(cd_msa) = cd_msa else {
        return PSIERR_BADPARAM;
    };
    let Some(dimensions) = cd_msa.dimensions else {
        return PSIERR_BADPARAM;
    };
    let q = dimensions.query_length as usize;
    let rows = dimensions.num_seqs as usize;
    if cd_msa.query.len() < q
        || cd_msa.msa.len() < rows
        || !cd_msa.msa.iter().take(rows).all(|row| row.len() >= q)
    {
        return PSIERR_BADPARAM;
    }
    let gap = crate::encoding::NCBISTDAA_GAP;
    for &letter in cd_msa.query.iter().take(q) {
        if letter == gap {
            return PSIERR_GAPINQUERY;
        }
    }
    for row in cd_msa.msa.iter().take(rows) {
        for cell in row.iter().take(q) {
            if !cell.is_aligned {
                continue;
            }
            let Some(data) = cell.data.as_ref() else {
                return PSIERR_BADPROFILE;
            };
            if data.wfreqs.len() < alphabet_size as usize || data.iobsr < POS_EPSILON {
                return PSIERR_BADPROFILE;
            }
            let mut sum = 0.0;
            for &weight in data.wfreqs.iter().take(alphabet_size as usize) {
                if weight < 0.0 {
                    return PSIERR_BADPROFILE;
                }
                sum += weight;
            }
            if (sum - 1.0).abs() > POS_EPSILON {
                return PSIERR_BADPROFILE;
            }
        }
    }
    PSI_SUCCESS
}

/// Port of NCBI `_PSISequenceLengthWithoutX`.
pub fn psi_sequence_length_without_x(seq: &[u8]) -> u32 {
    seq.iter()
        .filter(|&&residue| residue != crate::encoding::NCBISTDAA_X)
        .count() as u32
}

/// The 20 standard amino acid residue codes in NCBIstdaa encoding.
/// These are the residues that participate in scoring (indices 1..=20 and 22,
/// but excluding X, so we list the actual standard 20).
#[cfg(test)]
const STD_RESIDUES: [u8; EFFECTIVE_ALPHABET] = crate::encoding::NCBISTDAA_STANDARD_RESIDUES;

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
/// The AA_FREQUENCIES in matrix.rs are in alphabetical order
/// (A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y), corresponding to
/// [`crate::encoding::NCBISTDAA_STANDARD_RESIDUES`].
fn std_prob_ncbistdaa() -> [f64; AA_SIZE] {
    use crate::matrix::AA_FREQUENCIES;
    let mut prob = [0.0f64; AA_SIZE];
    for (freq, &code) in AA_FREQUENCIES
        .iter()
        .zip(crate::encoding::NCBISTDAA_STANDARD_RESIDUES.iter())
    {
        prob[code as usize] = *freq;
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

        Pssm {
            scores,
            length,
            info_content,
        }
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

    /// Return the consensus residue (highest-scoring amino acid) at a given position.
    pub fn consensus_at(&self, position: usize) -> Option<u8> {
        if position >= self.length {
            return None;
        }
        let row = &self.scores[position];
        // Only consider the 20 standard amino acids in NCBIstdaa.
        let (best_aa, _) = crate::encoding::NCBISTDAA_STANDARD_RESIDUES
            .iter()
            .copied()
            .map(|aa| (aa, row[aa as usize]))
            .max_by_key(|&(_, s)| s)?;
        Some(best_aa)
    }

    /// Update the PSSM from a set of aligned sequences (PSI-BLAST iteration).
    ///
    /// This implements the NCBI PSI-BLAST PSSM computation for BLOSUM62:
    ///   1. Henikoff position-based sequence weighting
    ///   2. Weighted residue frequency computation
    ///   3. Pseudocount blending via BLOSUM62 frequency ratios
    ///   4. Conversion to integer log-odds scores
    ///
    /// # Arguments
    /// * `aligned_seqs` - Aligned subject sequences in NCBIstdaa encoding.
    ///   Each sequence should be the same length as the query (self.length).
    ///   Gaps are represented as [`crate::encoding::NCBISTDAA_GAP`], unknown
    ///   as [`crate::encoding::NCBISTDAA_X`].
    /// * `bg_freqs` - Background amino acid frequencies (20 values in alphabetical
    ///   order: A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y). Used for reference
    ///   but the actual NCBIstdaa-indexed probabilities are computed internally.
    /// * `pseudocount_weight` - Fixed pseudocount weight when positive;
    ///   non-positive values use NCBI column-specific pseudocounts.
    pub fn update_from_alignment(
        &mut self,
        aligned_seqs: &[Vec<u8>],
        _bg_freqs: &[f64; 20],
        pseudocount_weight: f64,
    ) {
        self.update_from_alignment_with_matrix_and_pseudocount(
            aligned_seqs,
            "BLOSUM62",
            (pseudocount_weight > 0.0).then_some(pseudocount_weight),
        );
    }

    /// Update the PSSM using frequency-ratio and ungapped lambda tables for the
    /// named NCBI protein matrix.
    pub fn update_from_alignment_with_matrix(
        &mut self,
        aligned_seqs: &[Vec<u8>],
        matrix_name: &str,
    ) {
        self.update_from_alignment_with_matrix_and_pseudocount(aligned_seqs, matrix_name, None);
    }

    /// Update the PSSM with a named matrix and optional fixed pseudocount.
    ///
    /// When `fixed_pseudocount` is `None` or non-positive, the NCBI
    /// column-specific pseudocount method is used.
    pub fn update_from_alignment_with_matrix_and_pseudocount(
        &mut self,
        aligned_seqs: &[Vec<u8>],
        matrix_name: &str,
        fixed_pseudocount: Option<f64>,
    ) {
        if aligned_seqs.is_empty() {
            return;
        }

        let std_prob = std_prob_ncbistdaa();
        let freq_ratios = crate::matrix::get_matrix_freq_ratios_with_scale(matrix_name)
            .unwrap_or_else(|| crate::matrix::MatrixFreqRatios {
                data: BLOSUM62_FREQ_RATIOS,
                bit_scale_factor: 2,
            });
        let lambda = crate::stat::protein_ungapped_kbp_for_matrix(matrix_name).lambda;

        // Step 1: Compute Henikoff position-based sequence weights
        let match_weights = compute_sequence_weights_and_match_weights(aligned_seqs, self.length);

        // Step 2: Compute effective observations and column-specific pseudocounts
        // Step 3: Blend observed frequencies with pseudocounts using frequency ratios
        // Step 4: Convert to integer PSSM scores

        for pos in 0..self.length {
            // Check if this position has any aligned residues
            let total_weight: f64 = (0..AA_SIZE).map(|r| match_weights[pos][r]).sum();

            if total_weight < POS_EPSILON {
                // No alignment data at this position, keep original scores
                continue;
            }

            // Compute effective number of observations
            let observations = s_effective_observations(aligned_seqs, pos, self.length, &std_prob);

            let pseudo_weight = fixed_pseudocount
                .filter(|weight| *weight > 0.0)
                .unwrap_or_else(|| {
                    compute_column_pseudocounts(&match_weights[pos], &std_prob, observations)
                });

            // Compute frequency ratios with pseudocount blending
            for r in 0..AA_SIZE {
                if std_prob[r] <= POS_EPSILON {
                    // Score for non-standard residues: keep as-is or set to minimum
                    continue;
                }

                // pseudo = beta * sum_i(match_weights[pos][i] * freq_ratio[r][i])
                let mut pseudo = 0.0;
                for i in 0..AA_SIZE {
                    pseudo += match_weights[pos][i] * freq_ratios.data[r][i];
                }
                pseudo *= pseudo_weight;

                let numerator = observations * match_weights[pos][r] / std_prob[r] + pseudo;
                let denominator = observations + pseudo_weight;

                if denominator <= 0.0 {
                    continue;
                }

                let q_over_p = numerator / denominator;

                // freq_ratio = q_over_p * std_prob[r] (the C code multiplies then
                // divides by std_prob in the next stage -- we go directly to score)
                // Score = BLAST_Nint((1/lambda) * ln(q_over_p)).
                if q_over_p > POS_EPSILON {
                    let score = crate::math::blast_nint(q_over_p.ln() / lambda) as i32;
                    self.scores[pos][r] = score;
                } else {
                    self.scores[pos][r] = crate::stat::BLAST_SCORE_MIN;
                }
            }
        }

        // Compute information content
        self.info_content =
            compute_information_content(&self.scores, &std_prob, self.length, lambda);
    }
}

/// Compute Henikoff position-based sequence weights and match weights.
///
/// Each aligned row first receives one Henikoff sequence weight by summing
/// `1 / (distinct residues in column * count of this residue in column)` over
/// the row's aligned extent. The final per-column match weights are then the
/// normalized sequence weights for rows that contain a standard residue in that
/// column.
fn compute_sequence_weights_and_match_weights(
    aligned_seqs: &[Vec<u8>],
    query_length: usize,
) -> Vec<[f64; AA_SIZE]> {
    let num_seqs = aligned_seqs.len();
    let mut match_weights = vec![[0.0f64; AA_SIZE]; query_length];
    let extents: Vec<Option<(usize, usize)>> = aligned_seqs
        .iter()
        .map(|seq| aligned_sequence_extent(seq, query_length))
        .collect();
    let mut seq_weights = vec![0.0f64; num_seqs];

    for pos in 0..query_length {
        // Count occurrences of each residue at this position
        let mut residue_counts = [0u32; AA_SIZE];
        let mut num_aligned = 0u32;

        for (seq, extent) in aligned_seqs.iter().zip(extents.iter()) {
            if extent_contains(*extent, pos) && pos < seq.len() {
                let r = seq[pos] as usize;
                if is_effective_residue(seq[pos]) {
                    residue_counts[r] += 1;
                    num_aligned += 1;
                }
            }
        }

        if num_aligned == 0 {
            continue;
        }

        // Count distinct standard residues.
        let mut num_distinct = 0u32;
        for r in 0..AA_SIZE {
            if residue_counts[r] > 0 && is_effective_residue(r as u8) {
                num_distinct += 1;
            }
        }

        if num_distinct == 0 {
            // All gaps or X -- distribute weight uniformly
            continue;
        }

        for (si, seq) in aligned_seqs.iter().enumerate() {
            if extent_contains(extents[si], pos) && pos < seq.len() {
                let r = seq[pos] as usize;
                if is_effective_residue(seq[pos]) {
                    let w = 1.0 / (num_distinct as f64 * residue_counts[r] as f64);
                    seq_weights[si] += w;
                }
            }
        }
    }

    for pos in 0..query_length {
        let mut weight_sum = 0.0;
        for (si, seq) in aligned_seqs.iter().enumerate() {
            if extent_contains(extents[si], pos)
                && pos < seq.len()
                && is_effective_residue(seq[pos])
            {
                weight_sum += seq_weights[si];
            }
        }

        if weight_sum <= 0.0 {
            continue;
        }

        // Normalize and accumulate match weights
        for (si, seq) in aligned_seqs.iter().enumerate() {
            if extent_contains(extents[si], pos) && pos < seq.len() && seq_weights[si] > 0.0 {
                let r = seq[pos] as usize;
                if is_effective_residue(seq[pos]) {
                    let norm_weight = seq_weights[si] / weight_sum;
                    match_weights[pos][r] += norm_weight;
                }
            }
        }
    }

    match_weights
}

fn aligned_sequence_extent(seq: &[u8], query_length: usize) -> Option<(usize, usize)> {
    let end = seq.len().min(query_length);
    let left = seq
        .iter()
        .take(end)
        .position(|&r| r != crate::encoding::NCBISTDAA_X)?;
    let right = seq
        .iter()
        .take(end)
        .rposition(|&r| r != crate::encoding::NCBISTDAA_X)?;
    Some((left, right))
}

fn extent_contains(extent: Option<(usize, usize)>, pos: usize) -> bool {
    extent.is_some_and(|(left, right)| left <= pos && pos <= right)
}

/// Compute the effective number of independent observations at a column.
///
/// Port of NCBI `s_effectiveObservations` (`blast_psi_priv.c`). The C code
/// uses the aligned block spanning `columnNumber`; the compact Rust MSA rows
/// infer that block from non-X row extents.
fn s_effective_observations(
    aligned_seqs: &[Vec<u8>],
    pos: usize,
    query_length: usize,
    bg_prob: &[f64; AA_SIZE],
) -> f64 {
    if pos >= query_length {
        return 0.0;
    }

    let bg20 = effective_bg_probs(bg_prob);
    let expno = s_initialize_exp_num_observations(&bg20);
    let Some((left, right)) = aligned_block_extents(aligned_seqs, pos, query_length) else {
        return 0.0;
    };
    let pos_distinct_distrib = position_distinct_distribution(aligned_seqs, left, right);
    let num_participating = num_participating_at_position(aligned_seqs, pos) as f64;
    if num_participating <= 0.0 {
        return 0.0;
    }

    let ave_distinct = average_top_half_distinct_count(&pos_distinct_distrib, left, right);
    if ave_distinct <= 0.0 {
        return 0.0;
    }

    let indep = interpolate_independent_observations(ave_distinct, &expno);
    (indep.min(num_participating) - 1.0).max(0.0)
}

fn aligned_block_extents(
    aligned_seqs: &[Vec<u8>],
    pos: usize,
    query_length: usize,
) -> Option<(usize, usize)> {
    let mut left = usize::MAX;
    let mut right = 0usize;
    let mut found = false;

    for seq in aligned_seqs {
        let Some((seq_left, seq_right)) = aligned_sequence_extent(seq, query_length) else {
            continue;
        };
        if seq_left <= pos && pos <= seq_right {
            left = left.min(seq_left);
            right = right.max(seq_right);
            found = true;
        }
    }

    found.then_some((left, right))
}

/// Port of `s_initializeExpNumObservations`.
fn s_initialize_exp_num_observations(
    bg20: &[f64; EFFECTIVE_ALPHABET],
) -> [f64; MAX_IND_OBSERVATIONS + 1] {
    let mut expno = [0.0f64; MAX_IND_OBSERVATIONS + 1];
    for (j, slot) in expno
        .iter_mut()
        .enumerate()
        .take(MAX_IND_OBSERVATIONS)
        .skip(1)
    {
        let mut weighted_sum = 0.0;
        for &p in bg20 {
            weighted_sum += (j as f64 * (1.0 - p).ln()).exp();
        }
        *slot = EFFECTIVE_ALPHABET as f64 - weighted_sum;
    }
    expno
}

fn position_distinct_distribution(
    aligned_seqs: &[Vec<u8>],
    left: usize,
    right: usize,
) -> [i32; EFFECTIVE_ALPHABET + 1] {
    let mut distrib = [0i32; EFFECTIVE_ALPHABET + 1];
    for col in left..=right {
        let distinct = count_distinct_standard_residues_at(aligned_seqs, col);
        distrib[distinct] += 1;
    }
    distrib
}

fn count_distinct_standard_residues_at(aligned_seqs: &[Vec<u8>], pos: usize) -> usize {
    let mut residue_seen = [false; AA_SIZE];
    for seq in aligned_seqs {
        if let Some(&r) = seq.get(pos) {
            if is_effective_residue(r) {
                residue_seen[r as usize] = true;
            }
        }
    }
    residue_seen.iter().filter(|&&v| v).count()
}

fn num_participating_at_position(aligned_seqs: &[Vec<u8>], pos: usize) -> usize {
    aligned_seqs
        .iter()
        .filter(|seq| seq.get(pos).copied().is_some_and(is_effective_residue))
        .count()
}

fn is_effective_residue(r: u8) -> bool {
    crate::encoding::is_ncbistdaa_standard_residue(r)
}

fn average_top_half_distinct_count(
    pos_distinct_distrib: &[i32; EFFECTIVE_ALPHABET + 1],
    left: usize,
    right: usize,
) -> f64 {
    let half_num_columns = 1.max(((right as i32 - left as i32) + 2) / 2);
    let mut columns_accounted_for = 0i32;
    let mut total_distinct_counts = 0i32;
    let mut k = EFFECTIVE_ALPHABET as i32;

    while columns_accounted_for < half_num_columns && k >= 0 {
        let count = pos_distinct_distrib[k as usize];
        total_distinct_counts += count * k;
        columns_accounted_for += count;
        if columns_accounted_for > half_num_columns {
            total_distinct_counts -= (columns_accounted_for - half_num_columns) * k;
            columns_accounted_for = half_num_columns;
        }
        k -= 1;
    }

    if columns_accounted_for == 0 {
        0.0
    } else {
        total_distinct_counts as f64 / columns_accounted_for as f64
    }
}

fn interpolate_independent_observations(
    ave_distinct_aa: f64,
    expno: &[f64; MAX_IND_OBSERVATIONS + 1],
) -> f64 {
    let mut i = 1usize;
    while i < MAX_IND_OBSERVATIONS && expno[i] <= ave_distinct_aa {
        i += 1;
    }
    if i == MAX_IND_OBSERVATIONS {
        i as f64
    } else {
        let denom = expno[i] - expno[i - 1];
        if denom.abs() < f64::EPSILON {
            i as f64
        } else {
            i as f64 - (expno[i] - ave_distinct_aa) / denom
        }
    }
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
        info[pos] = sum / crate::math::NCBIMATH_LN2;
    }

    info
}

/// Run one iteration of PSI-BLAST.
/// Takes a query, a set of subject sequences, and the current PSSM.
/// Returns hits that pass the search e-value threshold.
/// Result of a PSI-BLAST hit.
#[derive(Debug, Clone)]
pub struct PsiBlastHit {
    pub subject_id: String,
    pub score: i32,
    pub evalue: f64,
    /// Start offset in query/PSSM where the best alignment begins.
    pub query_start: usize,
    /// End offset in query/PSSM, exclusive.
    pub query_end: usize,
    /// Start offset in subject where the best alignment begins.
    pub subject_start: usize,
    /// End offset in subject, exclusive.
    pub subject_end: usize,
    /// Length of the aligned query region, excluding query gaps.
    pub align_len: usize,
    /// Subject sequence length.
    pub subject_len: usize,
    /// Gapped query alignment, in ASCII amino-acid letters.
    pub query_aln: Vec<u8>,
    /// Gapped subject alignment, in ASCII amino-acid letters.
    pub subject_aln: Vec<u8>,
}

pub fn psi_blast_iteration(
    pssm: &Pssm,
    query: &[u8],
    subjects: &[(String, Vec<u8>)], // (id, NCBIstdaa sequence)
    evalue_threshold: f64,
    search_space: f64,
    kbp_lambda: f64,
    kbp_k: f64,
    gap_open: i32,
    gap_extend: i32,
    x_dropoff: i32,
    seed_cutoff: i32,
) -> Vec<PsiBlastHit> {
    let mut results = Vec::new();
    let pssm_rows = pssm_rows_for_alignment(pssm);
    let seed_cutoff = seed_cutoff.max(1);

    for (subj_id, subj_seq) in subjects {
        if subj_seq.is_empty() || pssm.length == 0 || query.len() < pssm.length {
            continue;
        }

        let candidates = pssm_ungapped_diagonal_candidates(pssm, subj_seq);
        let mut subject_hits = Vec::new();

        for candidate in candidates {
            if candidate.score < seed_cutoff {
                continue;
            }
            let gapped = crate::protein::protein_gapped_align_pssm(
                &query[..pssm.length],
                subj_seq,
                candidate.seed_query,
                candidate.seed_subject,
                0,
                &pssm_rows,
                gap_open,
                gap_extend,
                x_dropoff,
            );
            let (score, query_start, query_end, subject_start, subject_end, query_aln, subject_aln) =
                if let Some(gapped) = gapped {
                    let (query_aln, subject_aln) = gapped.edit_script.render_alignment(
                        &query[gapped.query_start..gapped.query_end],
                        &subj_seq[gapped.subject_start..gapped.subject_end],
                        crate::encoding::ncbistdaa_to_aminoacid_char,
                    );
                    (
                        gapped.score,
                        gapped.query_start,
                        gapped.query_end,
                        gapped.subject_start,
                        gapped.subject_end,
                        query_aln.into_bytes(),
                        subject_aln.into_bytes(),
                    )
                } else {
                    let query_end = candidate.query_start + candidate.align_len;
                    let subject_end = candidate.subject_start + candidate.align_len;
                    (
                        candidate.score,
                        candidate.query_start,
                        query_end,
                        candidate.subject_start,
                        subject_end,
                        crate::encoding::ncbistdaa_to_aminoacid_sequence(
                            &query[candidate.query_start..query_end],
                        ),
                        crate::encoding::ncbistdaa_to_aminoacid_sequence(
                            &subj_seq[candidate.subject_start..subject_end],
                        ),
                    )
                };

            let evalue = kbp_k * search_space * (-kbp_lambda * score as f64).exp();
            if evalue <= evalue_threshold {
                let hit = PsiBlastHit {
                    subject_id: subj_id.clone(),
                    score,
                    evalue,
                    query_start,
                    query_end,
                    subject_start,
                    subject_end,
                    align_len: query_end.saturating_sub(query_start),
                    subject_len: subj_seq.len(),
                    query_aln,
                    subject_aln,
                };
                if !subject_hits
                    .iter()
                    .any(|existing| same_pssm_hit(existing, &hit))
                {
                    subject_hits.push(hit);
                }
            }
        }

        subject_hits.sort_by(compare_pssm_hits);
        results.extend(subject_hits);
    }

    // Sort by evalue with NCBI `s_EvalueComp` denormal-region equivalence
    // (values < 1e-180 compare equal) — same helper used by HSP-list sorts.
    results.sort_by(compare_pssm_hits);
    results
}

fn same_pssm_hit(a: &PsiBlastHit, b: &PsiBlastHit) -> bool {
    a.score == b.score
        && a.query_start == b.query_start
        && a.query_end == b.query_end
        && a.subject_start == b.subject_start
        && a.subject_end == b.subject_end
        && a.query_aln == b.query_aln
        && a.subject_aln == b.subject_aln
}

fn compare_pssm_hits(a: &PsiBlastHit, b: &PsiBlastHit) -> std::cmp::Ordering {
    crate::hspstream::evalue_comp(a.evalue, b.evalue)
        .then_with(|| b.score.cmp(&a.score))
        .then_with(|| a.query_start.cmp(&b.query_start))
        .then_with(|| a.subject_start.cmp(&b.subject_start))
        .then_with(|| a.query_end.cmp(&b.query_end))
        .then_with(|| a.subject_end.cmp(&b.subject_end))
}

#[derive(Debug, Clone)]
struct PssmUngappedCandidate {
    score: i32,
    query_start: usize,
    subject_start: usize,
    align_len: usize,
    seed_query: usize,
    seed_subject: usize,
}

fn pssm_ungapped_diagonal_candidates(pssm: &Pssm, subj_seq: &[u8]) -> Vec<PssmUngappedCandidate> {
    let mut candidates = Vec::new();
    for diag_start_q in 0..pssm.length {
        collect_pssm_diagonal_candidate(pssm, subj_seq, diag_start_q, 0, &mut candidates);
    }
    for diag_start_s in 1..subj_seq.len() {
        collect_pssm_diagonal_candidate(pssm, subj_seq, 0, diag_start_s, &mut candidates);
    }
    candidates.sort_by(|a, b| {
        b.score
            .cmp(&a.score)
            .then_with(|| a.query_start.cmp(&b.query_start))
            .then_with(|| a.subject_start.cmp(&b.subject_start))
    });
    candidates
}

fn collect_pssm_diagonal_candidate(
    pssm: &Pssm,
    subj_seq: &[u8],
    diag_start_q: usize,
    diag_start_s: usize,
    candidates: &mut Vec<PssmUngappedCandidate>,
) {
    let diag_len = (pssm.length - diag_start_q).min(subj_seq.len() - diag_start_s);
    let mut score = 0i32;
    let mut run_query_start = diag_start_q;
    let mut run_subject_start = diag_start_s;
    let mut best: Option<PssmUngappedCandidate> = None;

    for offset in 0..diag_len {
        let q = diag_start_q + offset;
        let s = diag_start_s + offset;
        score += pssm.score_at(q, subj_seq[s]);
        if score <= 0 {
            if let Some(candidate) = best.take() {
                candidates.push(candidate);
            }
            score = 0;
            run_query_start = q + 1;
            run_subject_start = s + 1;
        } else if best
            .as_ref()
            .is_none_or(|candidate| score > candidate.score)
        {
            best = Some(PssmUngappedCandidate {
                score,
                query_start: run_query_start,
                subject_start: run_subject_start,
                align_len: q - run_query_start + 1,
                seed_query: q,
                seed_subject: s,
            });
        }
    }
    if let Some(candidate) = best {
        candidates.push(candidate);
    }
}

fn pssm_rows_for_alignment(pssm: &Pssm) -> Vec<Vec<i32>> {
    pssm.scores.iter().map(|row| row.to_vec()).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn simple_matrix() -> [[i32; AA_SIZE]; AA_SIZE] {
        let mut m = [[0i32; AA_SIZE]; AA_SIZE];
        for i in 1..21 {
            m[i][i] = 5;
        }
        for i in 1..21 {
            for j in 1..21 {
                if i != j {
                    m[i][j] = -2;
                }
            }
        }
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

        let results = psi_blast_iteration(
            &pssm, &query, &subjects, 1.0, 1000.0, 0.3176, 0.134, 11, 1, 50, 1,
        );
        assert!(!results.is_empty(), "Should find matching subject");
        assert_eq!(results[0].subject_id, "match");
        assert_eq!(results[0].query_start, 0);
        assert_eq!(results[0].subject_start, 0);
        assert_eq!(results[0].align_len, 5);
    }

    #[test]
    fn test_psi_blast_iteration_reports_query_offset_for_local_match() {
        let query = vec![1u8, 2, 3, 4, 5, 6];
        let matrix = simple_matrix();
        let pssm = Pssm::from_sequence(&query, &matrix);

        let subjects = vec![("suffix".to_string(), vec![3u8, 4, 5, 6])];
        let results = psi_blast_iteration(
            &pssm, &query, &subjects, 1.0, 1000.0, 0.3176, 0.134, 11, 1, 50, 1,
        );

        assert_eq!(results.len(), 1);
        assert_eq!(results[0].subject_id, "suffix");
        assert_eq!(results[0].query_start, 2);
        assert_eq!(results[0].subject_start, 0);
        assert_eq!(results[0].align_len, 4);
    }

    #[test]
    fn test_psi_blast_iteration_allows_single_residue_query() {
        let query = vec![crate::encoding::NCBISTDAA_A];
        let matrix = simple_matrix();
        let pssm = Pssm::from_sequence(&query, &matrix);
        let subjects = vec![("single".to_string(), query.clone())];

        let results = psi_blast_iteration(
            &pssm, &query, &subjects, 1.0e20, 1.0, 0.3176, 0.134, 11, 1, 50, 1,
        );

        assert_eq!(results.len(), 1);
        assert_eq!(results[0].subject_id, "single");
        assert_eq!(results[0].query_start, 0);
        assert_eq!(results[0].query_end, 1);
        assert_eq!(results[0].subject_start, 0);
        assert_eq!(results[0].subject_end, 1);
        assert_eq!(results[0].align_len, 1);
        assert_eq!(results[0].query_aln, b"A");
        assert_eq!(results[0].subject_aln, b"A");
    }

    #[test]
    fn test_psi_blast_iteration_applies_seed_cutoff_before_gapping() {
        let query = vec![1u8, 2, 3, 4, 5];
        let matrix = simple_matrix();
        let pssm = Pssm::from_sequence(&query, &matrix);
        let subjects = vec![("match".to_string(), query.clone())];

        let passing = psi_blast_iteration(
            &pssm, &query, &subjects, 1.0e20, 1.0, 0.3176, 0.134, 11, 1, 50, 1,
        );
        let filtered = psi_blast_iteration(
            &pssm, &query, &subjects, 1.0e20, 1.0, 0.3176, 0.134, 11, 1, 50, 10_000,
        );

        assert_eq!(passing.len(), 1);
        assert!(
            filtered.is_empty(),
            "gap-trigger seed cutoff should filter low ungapped-scoring candidates"
        );
    }

    #[test]
    fn test_psi_blast_iteration_keeps_multiple_diagonal_candidates() {
        let query = vec![
            crate::encoding::NCBISTDAA_A,
            crate::encoding::NCBISTDAA_C,
            crate::encoding::NCBISTDAA_D,
            crate::encoding::NCBISTDAA_E,
        ];
        let matrix = simple_matrix();
        let pssm = Pssm::from_sequence(&query, &matrix);
        let subject = vec![
            crate::encoding::NCBISTDAA_A,
            crate::encoding::NCBISTDAA_C,
            19, // V
            crate::encoding::NCBISTDAA_D,
            crate::encoding::NCBISTDAA_E,
        ];

        let candidates = pssm_ungapped_diagonal_candidates(&pssm, &subject);

        assert!(
            candidates
                .iter()
                .any(|c| c.query_start == 0 && c.subject_start == 0 && c.align_len == 2),
            "first positive diagonal run should be retained: {:?}",
            candidates
        );
        assert!(
            candidates
                .iter()
                .any(|c| c.query_start == 2 && c.subject_start == 3 && c.align_len == 2),
            "later positive diagonal run should also be retained: {:?}",
            candidates
        );
    }

    #[test]
    fn test_psi_blast_iteration_reports_multiple_hsps_per_subject() {
        let query = vec![
            crate::encoding::NCBISTDAA_A,
            crate::encoding::NCBISTDAA_C,
            crate::encoding::NCBISTDAA_D,
            crate::encoding::NCBISTDAA_A,
            crate::encoding::NCBISTDAA_C,
            crate::encoding::NCBISTDAA_D,
        ];
        let matrix = simple_matrix();
        let pssm = Pssm::from_sequence(&query, &matrix);
        let y = crate::encoding::aminoacid_to_ncbistdaa_base(b'Y');
        let subject = vec![
            crate::encoding::NCBISTDAA_A,
            crate::encoding::NCBISTDAA_C,
            crate::encoding::NCBISTDAA_D,
            y,
            y,
            crate::encoding::NCBISTDAA_A,
            crate::encoding::NCBISTDAA_C,
            crate::encoding::NCBISTDAA_D,
        ];
        let subjects = vec![("two_blocks".to_string(), subject)];

        let results = psi_blast_iteration(
            &pssm, &query, &subjects, 1.0e20, 1.0, 0.3176, 0.134, 100, 10, 10, 1,
        );
        let coords: Vec<(usize, usize, usize, usize)> = results
            .iter()
            .map(|h| (h.query_start, h.query_end, h.subject_start, h.subject_end))
            .collect();

        assert!(
            coords.contains(&(0, 3, 0, 3)),
            "first local block should be reported: {:?}",
            coords
        );
        assert!(
            coords.contains(&(3, 6, 5, 8)),
            "second local block should be reported: {:?}",
            coords
        );
    }

    #[test]
    fn test_henikoff_weights_identical_seqs() {
        // All identical sequences should produce uniform weights
        let seqs = vec![
            vec![1u8, 4, 7, 10], // A, D, G, K
            vec![1u8, 4, 7, 10],
            vec![1u8, 4, 7, 10],
        ];
        let mw = compute_sequence_weights_and_match_weights(&seqs, 4);
        // Each position should have weight 1.0 for the single residue present
        for pos in 0..4 {
            let total: f64 = mw[pos].iter().sum();
            assert!(
                (total - 1.0).abs() < 0.01,
                "Match weights should sum to ~1 at pos {}, got {}",
                pos,
                total
            );
        }
    }

    #[test]
    fn test_henikoff_weights_diverse_seqs() {
        // Diverse sequences: rare residues should get more weight
        let seqs = vec![
            vec![1u8, 4], // A, D
            vec![1u8, 4], // A, D
            vec![1u8, 5], // A, E
        ];
        let mw = compute_sequence_weights_and_match_weights(&seqs, 2);
        // At position 0: all have A, so weight should go to A=1
        assert!(mw[0][1] > 0.9, "All-A column should have weight ~1 on A");
        // Sequence-level Henikoff weights include the conserved A column, so
        // the two D rows remain slightly heavier than the one E row.
        assert!(
            (mw[1][4] - 7.0 / 12.0).abs() < 0.01,
            "D weight should be ~0.583"
        );
        assert!(
            (mw[1][5] - 5.0 / 12.0).abs() < 0.01,
            "E weight should be ~0.417"
        );
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
        assert!(
            pseudo < 100.0,
            "Conserved column pseudo should be small, got {}",
            pseudo
        );
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
        assert!(
            pseudo > 10.0,
            "Uniform column pseudo should be substantial, got {}",
            pseudo
        );
    }

    #[test]
    fn test_update_from_alignment_changes_scores() {
        let query = vec![1u8, 4, 7, 10, 13]; // A, D, G, K, N
        let matrix = simple_matrix();
        let mut pssm = Pssm::from_sequence(&query, &matrix);
        let original_scores = pssm.scores.clone();

        // Create aligned sequences that mostly match query but with some variation
        let aligned = vec![
            vec![1u8, 4, 7, 10, 13], // exact match
            vec![1u8, 4, 7, 10, 13], // exact match
            vec![1u8, 5, 7, 10, 4],  // E at pos 1, D at pos 4
            vec![1u8, 4, 8, 10, 13], // H at pos 2
        ];

        let bg = crate::matrix::AA_FREQUENCIES;
        pssm.update_from_alignment(&aligned, &bg, 10.0);

        // Scores should have changed
        assert_ne!(
            pssm.scores, original_scores,
            "PSSM scores should change after alignment update"
        );
    }

    #[test]
    fn test_update_from_alignment_uses_named_matrix_tables() {
        let query = vec![1u8, 4, 7, 10, 13]; // A, D, G, K, N
        let aligned = vec![
            vec![1u8, 4, 7, 10, 13],
            vec![1u8, 5, 7, 10, 13],
            vec![3u8, 4, 8, 11, 13],
        ];
        let mut blosum62_pssm = Pssm::from_sequence(&query, &crate::matrix::BLOSUM62);
        let mut blosum45_pssm = Pssm::from_sequence(&query, &crate::matrix::BLOSUM62);

        blosum62_pssm.update_from_alignment_with_matrix(&aligned, "BLOSUM62");
        blosum45_pssm.update_from_alignment_with_matrix(&aligned, "BLOSUM45");

        assert_ne!(
            blosum62_pssm.scores, blosum45_pssm.scores,
            "PSSM update should use the selected matrix frequency ratios and lambda"
        );
    }

    #[test]
    fn test_update_from_alignment_fixed_pseudocount_changes_scores() {
        let query = vec![1u8, 4, 7, 10, 13]; // A, D, G, K, N
        let aligned = vec![
            vec![1u8, 4, 7, 10, 13],
            vec![1u8, 5, 7, 10, 13],
            vec![3u8, 4, 8, 11, 13],
            vec![1u8, 4, 7, 10, 4],
        ];
        let mut dynamic_pssm = Pssm::from_sequence(&query, &crate::matrix::BLOSUM62);
        let mut fixed_pssm = Pssm::from_sequence(&query, &crate::matrix::BLOSUM62);
        let mut fixed_via_background_pssm = Pssm::from_sequence(&query, &crate::matrix::BLOSUM62);

        dynamic_pssm.update_from_alignment_with_matrix(&aligned, "BLOSUM62");
        fixed_pssm.update_from_alignment_with_matrix_and_pseudocount(
            &aligned,
            "BLOSUM62",
            Some(50.0),
        );
        fixed_via_background_pssm.update_from_alignment(
            &aligned,
            &crate::matrix::AA_FREQUENCIES,
            50.0,
        );

        assert_ne!(
            dynamic_pssm.scores, fixed_pssm.scores,
            "fixed pseudocount should use a distinct scoring path"
        );
        assert_eq!(
            fixed_pssm.scores, fixed_via_background_pssm.scores,
            "background-frequency PSSM update should select the fixed pseudocount path"
        );
    }

    #[test]
    fn test_effective_observations() {
        let std_prob = std_prob_ncbistdaa();

        // With many diverse sequences, effective observations should be substantial
        let seqs: Vec<Vec<u8>> = (0..20).map(|i| vec![STD_RESIDUES[i % 20]]).collect();

        let obs = s_effective_observations(&seqs, 0, 1, &std_prob);
        assert!(obs > 0.0, "Should have positive effective observations");
    }

    #[test]
    fn test_henikoff_weights_ignore_non_standard_residues() {
        let seqs = vec![
            vec![crate::encoding::NCBISTDAA_A],
            vec![crate::encoding::NCBISTDAA_B],
            vec![crate::encoding::NCBISTDAA_X],
            vec![crate::encoding::NCBISTDAA_U],
            vec![crate::encoding::NCBISTDAA_STOP],
        ];

        let mw = compute_sequence_weights_and_match_weights(&seqs, 1);

        assert!((mw[0][crate::encoding::NCBISTDAA_A as usize] - 1.0).abs() < 1e-12);
        assert_eq!(mw[0][crate::encoding::NCBISTDAA_B as usize], 0.0);
        assert_eq!(mw[0][crate::encoding::NCBISTDAA_X as usize], 0.0);
        assert_eq!(mw[0][crate::encoding::NCBISTDAA_U as usize], 0.0);
        assert_eq!(mw[0][crate::encoding::NCBISTDAA_STOP as usize], 0.0);
    }

    #[test]
    fn test_henikoff_weights_use_aligned_row_extents() {
        let seqs = vec![
            vec![
                crate::encoding::NCBISTDAA_A,
                crate::encoding::NCBISTDAA_C,
                crate::encoding::NCBISTDAA_X,
            ],
            vec![
                crate::encoding::NCBISTDAA_X,
                crate::encoding::NCBISTDAA_C,
                crate::encoding::NCBISTDAA_D,
            ],
        ];

        let mw = compute_sequence_weights_and_match_weights(&seqs, 3);

        assert_eq!(mw[0][crate::encoding::NCBISTDAA_A as usize], 1.0);
        assert_eq!(mw[2][crate::encoding::NCBISTDAA_D as usize], 1.0);
        assert!((mw[1][crate::encoding::NCBISTDAA_C as usize] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_effective_observation_counts_ignore_non_standard_residues() {
        let seqs = vec![
            vec![crate::encoding::NCBISTDAA_A],
            vec![crate::encoding::NCBISTDAA_B],
            vec![crate::encoding::NCBISTDAA_X],
            vec![crate::encoding::NCBISTDAA_U],
            vec![crate::encoding::NCBISTDAA_STOP],
        ];

        assert_eq!(count_distinct_standard_residues_at(&seqs, 0), 1);
        assert_eq!(num_participating_at_position(&seqs, 0), 1);
    }

    #[test]
    fn test_effective_observations_use_spanning_alignment_block() {
        let seqs = vec![
            vec![
                crate::encoding::NCBISTDAA_A,
                crate::encoding::NCBISTDAA_C,
                crate::encoding::NCBISTDAA_X,
                crate::encoding::NCBISTDAA_X,
            ],
            vec![
                crate::encoding::NCBISTDAA_X,
                crate::encoding::NCBISTDAA_C,
                crate::encoding::NCBISTDAA_D,
                crate::encoding::NCBISTDAA_X,
            ],
            vec![
                crate::encoding::NCBISTDAA_X,
                crate::encoding::NCBISTDAA_X,
                crate::encoding::NCBISTDAA_D,
                crate::encoding::NCBISTDAA_E,
            ],
        ];

        assert_eq!(aligned_block_extents(&seqs, 0, 4), Some((0, 1)));
        assert_eq!(aligned_block_extents(&seqs, 1, 4), Some((0, 2)));
        assert_eq!(aligned_block_extents(&seqs, 2, 4), Some((1, 3)));
        assert_eq!(aligned_block_extents(&seqs, 3, 4), Some((2, 3)));
    }

    // ---- Tests ported from pssmcreate_unit_test.cpp / pssmenginefreqratios_unit_test.cpp ----

    #[test]
    fn test_pssm_henikoff_weights_three_diverse_seqs() {
        // Three sequences with completely different residues at each position.
        // Henikoff weighting should assign more uniform weights (~1/3 each).
        let seqs = vec![
            vec![1u8, 4, 7, 10],  // A, D, G, K
            vec![3u8, 5, 8, 11],  // C, E, H, L
            vec![6u8, 9, 12, 13], // F, I, M, N
        ];
        let mw = compute_sequence_weights_and_match_weights(&seqs, 4);

        for pos in 0..4 {
            let total: f64 = mw[pos].iter().sum();
            assert!(
                (total - 1.0).abs() < 0.01,
                "Match weights should sum to ~1 at pos {}, got {}",
                pos,
                total
            );

            // Each of the 3 residues should get weight ~1/3
            let nonzero: Vec<f64> = mw[pos].iter().copied().filter(|&v| v > 0.0).collect();
            assert_eq!(
                nonzero.len(),
                3,
                "Should have 3 non-zero weights at pos {}",
                pos
            );
            for &w in &nonzero {
                assert!(
                    (w - 1.0 / 3.0).abs() < 0.05,
                    "Each weight should be ~0.333, got {} at pos {}",
                    w,
                    pos
                );
            }
        }
    }

    #[test]
    fn test_pssm_henikoff_weights_two_identical_one_different() {
        // Two identical seqs and one different. Henikoff weighting assigns
        // each distinct residue type equal total weight (1/num_distinct each),
        // then normalizes. So with 2 distinct types, each gets 0.5.
        // The key insight: the *per-sequence* weight for the rare residue's
        // sequence (0.5) is higher than each identical sequence's weight (0.25).
        let seqs = vec![
            vec![1u8, 4], // A, D
            vec![1u8, 4], // A, D  (duplicate)
            vec![3u8, 5], // C, E  (different)
        ];
        let mw = compute_sequence_weights_and_match_weights(&seqs, 2);

        // At each position, 2 distinct residues -> each residue type gets equal weight
        // Position 0: A (2 seqs) and C (1 seq) each get normalized weight 0.5
        let weight_a = mw[0][1]; // A = NCBIstdaa 1
        let weight_c = mw[0][3]; // C = NCBIstdaa 3
        assert!(
            (weight_a - 0.5).abs() < 0.01,
            "A weight should be ~0.5, got {}",
            weight_a
        );
        assert!(
            (weight_c - 0.5).abs() < 0.01,
            "C weight should be ~0.5, got {}",
            weight_c
        );

        // The per-sequence weights show the difference: the unique sequence
        // contributes its full 0.5 alone, while the two identical sequences
        // split their 0.5. Verify by checking a 4-sequence case where the
        // imbalance is more visible in a different way.
        let seqs2 = vec![
            vec![1u8], // A
            vec![1u8], // A
            vec![1u8], // A
            vec![3u8], // C (rare)
        ];
        let mw2 = compute_sequence_weights_and_match_weights(&seqs2, 1);
        // 2 distinct, A=3, C=1. Henikoff: A seqs each get 1/(2*3), C gets 1/(2*1)
        // Sum = 3/(2*3) + 1/2 = 0.5 + 0.5 = 1.0.  After normalization: A=0.5, C=0.5
        // Equal residue-type weights, but per-seq the rare one is 3x heavier.
        assert!((mw2[0][1] - 0.5).abs() < 0.01, "A total weight ~0.5");
        assert!((mw2[0][3] - 0.5).abs() < 0.01, "C total weight ~0.5");
    }

    #[test]
    fn test_pssm_frequency_computation() {
        // Given known sequences, verify the match weight (frequency) matrix values.
        // 4 seqs, length 2. Position 0: all Ala. Position 1: 2x Asp, 1x Glu, 1x Ala.
        let seqs = vec![
            vec![1u8, 4], // A, D
            vec![1u8, 4], // A, D
            vec![1u8, 5], // A, E
            vec![1u8, 1], // A, A
        ];
        let mw = compute_sequence_weights_and_match_weights(&seqs, 2);

        // Position 0: only A present, weight should be 1.0
        assert!(
            (mw[0][1] - 1.0).abs() < 0.01,
            "All-A column should have weight 1.0 on A"
        );

        // Position 1: 3 distinct residues (D=2, E=1, A=1)
        // Henikoff: D seqs get 1/(3*2) each, E seq gets 1/(3*1), A seq gets 1/(3*1)
        // sum = 2/(3*2) + 1/3 + 1/3 = 1/3 + 1/3 + 1/3 = 1.0
        // After normalization: D = 1/3, E = 1/3, A = 1/3
        let total_pos1: f64 = mw[1].iter().sum();
        assert!(
            (total_pos1 - 1.0).abs() < 0.01,
            "Position 1 weights should sum to 1.0, got {}",
            total_pos1
        );
        assert!(mw[1][4] > 0.0, "D should have nonzero weight");
        assert!(mw[1][5] > 0.0, "E should have nonzero weight");
        assert!(mw[1][1] > 0.0, "A should have nonzero weight");
    }

    #[test]
    fn test_pssm_pseudocount_blending() {
        // Verify pseudocounts blend PSSM frequencies with BLOSUM62 background.
        // A conserved column should have low pseudocount weight,
        // a diverse column should have higher pseudocount weight.
        let std_prob = std_prob_ncbistdaa();

        // Highly conserved: all weight on Ala
        let mut mw_conserved = [0.0f64; AA_SIZE];
        mw_conserved[1] = 1.0;
        let pseudo_conserved = compute_column_pseudocounts(&mw_conserved, &std_prob, 50.0);

        // More diverse: weight spread across 5 residues
        let mut mw_diverse = [0.0f64; AA_SIZE];
        for &r in &[1u8, 4, 7, 10, 13] {
            mw_diverse[r as usize] = 0.2;
        }
        let pseudo_diverse = compute_column_pseudocounts(&mw_diverse, &std_prob, 50.0);

        // Diverse column should get more pseudocounts (lower information content
        // leads to higher alpha and thus higher pseudocount weight)
        assert!(
            pseudo_diverse > pseudo_conserved,
            "Diverse column pseudo ({}) should exceed conserved ({})",
            pseudo_diverse,
            pseudo_conserved
        );
    }

    #[test]
    fn test_pssm_score_matrix_signs() {
        // In the final PSSM score matrix after alignment update,
        // self-match scores should generally be positive and mismatches negative.
        let query = vec![1u8, 4, 7, 10, 13]; // A, D, G, K, N
        let matrix = simple_matrix();
        let mut pssm = Pssm::from_sequence(&query, &matrix);

        // Create aligned sequences that strongly conserve the query residues
        let aligned = vec![
            vec![1u8, 4, 7, 10, 13],
            vec![1u8, 4, 7, 10, 13],
            vec![1u8, 4, 7, 10, 13],
            vec![1u8, 4, 7, 10, 13],
            vec![1u8, 4, 7, 10, 13],
            vec![1u8, 5, 7, 10, 13], // slight variation at pos 1
        ];

        let bg = crate::matrix::AA_FREQUENCIES;
        pssm.update_from_alignment(&aligned, &bg, 10.0);

        // For each position, the query residue score should be positive
        for (pos, &qr) in query.iter().enumerate() {
            let self_score = pssm.score_at(pos, qr);
            assert!(
                self_score > 0,
                "Self-match score at pos {} for residue {} should be positive, got {}",
                pos,
                qr,
                self_score
            );
        }

        // At least some mismatches should be negative
        let mut negative_count = 0;
        for pos in 0..query.len() {
            for &r in STD_RESIDUES.iter() {
                if r != query[pos] && pssm.score_at(pos, r) < 0 {
                    negative_count += 1;
                }
            }
        }
        assert!(
            negative_count > 0,
            "Should have at least some negative mismatch scores"
        );
    }

    #[test]
    fn test_pssm_from_single_sequence() {
        // Building a PSSM from a single sequence should approximate the scoring matrix.
        // After update_from_alignment with just 1 sequence, pseudocounts dominate
        // and the result should be driven by BLOSUM62 frequency ratios.
        let query = vec![1u8, 7, 11, 17]; // A, G, L, S
        let matrix = simple_matrix();
        let mut pssm = Pssm::from_sequence(&query, &matrix);

        let aligned = vec![
            vec![1u8, 7, 11, 17], // single sequence = query itself
        ];
        let bg = crate::matrix::AA_FREQUENCIES;
        pssm.update_from_alignment(&aligned, &bg, 10.0);

        // With a single sequence, the PSSM should still assign positive scores
        // to the query residues (BLOSUM62 self-scores are positive).
        for (pos, &qr) in query.iter().enumerate() {
            let score = pssm.score_at(pos, qr);
            assert!(
                score > 0,
                "Single-seq PSSM should have positive self-score at pos {} (residue {}), got {}",
                pos,
                qr,
                score
            );
        }

        // The overall pattern should resemble BLOSUM62: self > mismatch on average
        let mut self_total = 0i64;
        let mut mismatch_total = 0i64;
        let mut mismatch_count = 0i64;
        for (pos, &qr) in query.iter().enumerate() {
            self_total += pssm.score_at(pos, qr) as i64;
            for &r in STD_RESIDUES.iter() {
                if r != qr {
                    mismatch_total += pssm.score_at(pos, r) as i64;
                    mismatch_count += 1;
                }
            }
        }
        let avg_self = self_total as f64 / query.len() as f64;
        let avg_mismatch = mismatch_total as f64 / mismatch_count as f64;
        assert!(
            avg_self > avg_mismatch,
            "Average self-score ({}) should exceed average mismatch ({})",
            avg_self,
            avg_mismatch
        );
    }

    #[test]
    fn test_pssm_effective_observations() {
        // Test effective observations with known inputs.
        let std_prob = std_prob_ncbistdaa();

        // All identical residues: 1 distinct -> low effective observations
        let seqs_identical: Vec<Vec<u8>> = (0..10).map(|_| vec![1u8]).collect();
        let obs_identical = s_effective_observations(&seqs_identical, 0, 1, &std_prob);

        // All different residues: 10 distinct -> higher effective observations
        let seqs_diverse: Vec<Vec<u8>> = (0..10).map(|i| vec![STD_RESIDUES[i % 20]]).collect();
        let obs_diverse = s_effective_observations(&seqs_diverse, 0, 1, &std_prob);

        assert!(
            obs_diverse > obs_identical,
            "Diverse column ({}) should have more effective observations than identical ({})",
            obs_diverse,
            obs_identical
        );

        // Identical seqs: 1 distinct residue among 10 seqs
        // The formula caps at num_participating and subtracts 1, minimum 0.
        // With 1 distinct residue, the expected value lookup finds j=1 quickly,
        // so effective obs should be small.
        assert!(
            obs_identical < 5.0,
            "Identical seqs should have low effective observations, got {}",
            obs_identical
        );

        // 10 diverse residues should yield substantially more
        assert!(
            obs_diverse > 5.0,
            "10 diverse residues should yield substantial observations, got {}",
            obs_diverse
        );
    }

    #[test]
    fn test_effective_observations_uses_ncbi_top_half_distinct_columns() {
        let seqs = vec![
            vec![1u8, 1, 1, 1],
            vec![1u8, 3, 3, 3],
            vec![1u8, 1, 4, 4],
            vec![1u8, 3, 4, 5],
        ];
        let distrib = position_distinct_distribution(&seqs, 0, 3);

        assert_eq!(distrib[1], 1);
        assert_eq!(distrib[2], 1);
        assert_eq!(distrib[3], 1);
        assert_eq!(distrib[4], 1);
        assert_eq!(average_top_half_distinct_count(&distrib, 0, 3), 3.5);
    }

    #[test]
    fn test_pssm_purging_identical_sequences() {
        // If all sequences are identical, the effective observations should be
        // much lower than when sequences are diverse. This tests the conceptual
        // equivalent of NCBI's sequence purging: identical sequences contribute
        // less independent information.
        let std_prob = std_prob_ncbistdaa();

        // 20 identical sequences
        let seqs_identical: Vec<Vec<u8>> = (0..20)
            .map(|_| {
                vec![1u8, 4, 7, 10] // A, D, G, K
            })
            .collect();

        // 20 sequences with variation at each position
        let seqs_varied: Vec<Vec<u8>> = (0..20)
            .map(|i| {
                vec![
                    STD_RESIDUES[i % 20],
                    STD_RESIDUES[(i + 3) % 20],
                    STD_RESIDUES[(i + 7) % 20],
                    STD_RESIDUES[(i + 11) % 20],
                ]
            })
            .collect();

        // Check effective observations at each position
        for pos in 0..4 {
            let obs_ident = s_effective_observations(&seqs_identical, pos, 4, &std_prob);
            let obs_varied = s_effective_observations(&seqs_varied, pos, 4, &std_prob);

            assert!(
                obs_varied > obs_ident,
                "At pos {}: varied obs ({}) should exceed identical obs ({})",
                pos,
                obs_varied,
                obs_ident
            );
        }

        // Also verify via match weights: identical sequences produce less diversity
        let mw_ident = compute_sequence_weights_and_match_weights(&seqs_identical, 4);
        let mw_varied = compute_sequence_weights_and_match_weights(&seqs_varied, 4);

        // Identical seqs: each position has exactly 1 nonzero weight entry
        for pos in 0..4 {
            let nonzero_ident = mw_ident[pos].iter().filter(|&&v| v > 0.0).count();
            let nonzero_varied = mw_varied[pos].iter().filter(|&&v| v > 0.0).count();
            assert!(
                nonzero_varied > nonzero_ident,
                "Varied alignment should have more nonzero weights at pos {} ({} vs {})",
                pos,
                nonzero_varied,
                nonzero_ident
            );
        }
    }

    #[test]
    fn translated_psi_public_allocators_match_c_shapes() {
        let dims = PSIMsaDimensions {
            query_length: 3,
            num_seqs: 2,
        };
        let mut msa = psi_public_msa_new(Some(&dims)).expect("msa");
        assert_eq!(msa.data.len(), 3);
        assert_eq!(msa.data[0].len(), 3);
        assert_eq!(msa.data[0][0], PSIMsaCell::default());
        msa.data[0][0] = PSIMsaCell {
            letter: crate::encoding::NCBISTDAA_A,
            is_aligned: true,
        };
        assert!(psi_public_msa_free(Some(msa.clone())).is_none());
        assert!(psi_public_msa_new(None).is_none());

        let matrix = psi_matrix_new(3, AA_SIZE as u32).expect("matrix");
        assert_eq!(matrix.ncols, 3);
        assert_eq!(matrix.nrows, AA_SIZE as u32);
        assert_eq!(matrix.pssm.len(), 3);
        assert_eq!(matrix.pssm[0].len(), AA_SIZE);
        assert_eq!(matrix.lambda, 0.0);
        assert!(psi_matrix_free(Some(matrix)).is_none());
    }

    #[test]
    fn translated_psi_save_pssm_copies_matrix_and_karlin_fields() {
        let internal = PsiInternalPssmData {
            ncols: 2,
            nrows: 3,
            pssm: vec![vec![1, 2, 3], vec![4, 5, 6]],
            scaled_pssm: Vec::new(),
            freq_ratios: Vec::new(),
            pseudocounts: Vec::new(),
        };
        let mut sbp = crate::stat::blast_score_blk_new(crate::encoding::BLASTAA_SEQ_CODE, 1)
            .expect("score block");
        sbp.kbp_gap_psi[0] = crate::stat::KarlinBlk {
            lambda: 0.267,
            k: 0.041,
            log_k: 0.041_f64.ln(),
            h: 0.14,
            round_down: false,
        };
        sbp.kbp_psi[0] = crate::stat::KarlinBlk {
            lambda: 0.318,
            k: 0.134,
            log_k: 0.134_f64.ln(),
            h: 0.401,
            round_down: false,
        };
        let mut out = psi_matrix_new(2, 3).expect("matrix");

        assert_eq!(
            s_psi_save_pssm(Some(&internal), Some(&sbp), Some(&mut out)),
            PSI_SUCCESS
        );

        assert_eq!(out.pssm, internal.pssm);
        assert_eq!((out.lambda, out.kappa, out.h), (0.267, 0.041, 0.14));
        assert_eq!(
            (out.ung_lambda, out.ung_kappa, out.ung_h),
            (0.318, 0.134, 0.401)
        );
        assert_eq!(
            s_psi_save_pssm(None, Some(&sbp), Some(&mut out)),
            PSIERR_BADPARAM
        );
        assert_eq!(
            s_psi_save_pssm(Some(&internal), Some(&sbp), None),
            PSIERR_BADPARAM
        );

        let mut wrong_dims = psi_matrix_new(1, 3).expect("matrix");
        assert_eq!(
            s_psi_save_pssm(Some(&internal), Some(&sbp), Some(&mut wrong_dims)),
            PSIERR_BADPARAM
        );

        let mut truncated = internal.clone();
        truncated.pssm[1].truncate(2);
        assert_eq!(
            s_psi_save_pssm(Some(&truncated), Some(&sbp), Some(&mut out)),
            PSIERR_BADPARAM
        );

        let mut missing_gap_stats = sbp.clone();
        missing_gap_stats.kbp_gap_psi.clear();
        assert_eq!(
            s_psi_save_pssm(Some(&internal), Some(&missing_gap_stats), Some(&mut out)),
            PSIERR_BADPARAM
        );

        let mut missing_ungapped_stats = sbp.clone();
        missing_ungapped_stats.kbp_psi.clear();
        assert_eq!(
            s_psi_save_pssm(
                Some(&internal),
                Some(&missing_ungapped_stats),
                Some(&mut out)
            ),
            PSIERR_BADPARAM
        );
    }

    #[test]
    fn translated_psi_diagnostics_allocators_follow_request_flags() {
        let blank = psi_diagnostics_request_new();
        assert!(!blank.frequency_ratios);
        let ascii = psi_diagnostics_request_new_ex(true);
        assert!(ascii.frequency_ratios);
        assert!(ascii.information_content);
        assert!(ascii.weighted_residue_frequencies);
        assert!(ascii.num_matching_seqs);
        let mut with_observations = ascii.clone();
        with_observations.independent_observations = true;

        let response = psi_diagnostics_response_new(4, AA_SIZE as u32, Some(&with_observations))
            .expect("response");
        assert_eq!(response.query_length, 4);
        assert_eq!(response.alphabet_size, AA_SIZE as u32);
        assert_eq!(response.information_content.as_ref().unwrap().len(), 4);
        assert!(response.residue_freqs.is_none());
        assert_eq!(
            response.frequency_ratios.as_ref().unwrap()[0].len(),
            AA_SIZE
        );
        assert_eq!(response.independent_observations.as_ref().unwrap().len(), 4);
        let mut freed = Some(response);
        assert!(psi_diagnostics_response_free(freed.take()).is_none());
        assert!(psi_diagnostics_response_new(4, AA_SIZE as u32, None).is_none());
        assert!(psi_diagnostics_request_free(Some(ascii)).is_none());
    }

    #[test]
    fn translated_psi_internal_msa_and_weights_allocators_match_c_shapes() {
        let dims = PSIMsaDimensions {
            query_length: 2,
            num_seqs: 1,
        };
        let mut msa = psi_public_msa_new(Some(&dims)).expect("msa");
        msa.data[0][0] = PSIMsaCell {
            letter: crate::encoding::NCBISTDAA_A,
            is_aligned: true,
        };
        msa.data[0][1] = PSIMsaCell {
            letter: crate::encoding::NCBISTDAA_C,
            is_aligned: true,
        };
        msa.data[1][0] = PSIMsaCell {
            letter: crate::encoding::NCBISTDAA_A,
            is_aligned: true,
        };
        msa.data[1][1] = PSIMsaCell {
            letter: crate::encoding::NCBISTDAA_X,
            is_aligned: false,
        };

        let packed = psi_packed_msa_new(Some(&msa)).expect("packed");
        assert_eq!(packed.use_sequence, vec![true, true]);
        assert_eq!(psi_packed_msa_get_number_of_aligned_seqs(Some(&packed)), 1);
        let internal = psi_msa_new(Some(&packed), AA_SIZE as u32).expect("internal");
        assert_eq!(internal.dimensions.num_seqs, 1);
        assert_eq!(
            internal.query,
            vec![crate::encoding::NCBISTDAA_A, crate::encoding::NCBISTDAA_C]
        );
        assert_eq!(
            internal.residue_counts[0][crate::encoding::NCBISTDAA_A as usize],
            2
        );
        assert_eq!(internal.num_matching_seqs[1], 1);
        assert!(psi_msa_free(Some(internal)).is_none());
        assert!(psi_packed_msa_free(Some(packed)).is_none());

        let block = psi_aligned_block_new(2).expect("block");
        assert_eq!(block.pos_extnt[0], SSeqRange { left: -1, right: 2 });
        assert_eq!(block.size, vec![0, 0]);

        let weights = psi_sequence_weights_new(Some(&dims), AA_SIZE).expect("weights");
        assert_eq!(weights.match_weights.len(), 2);
        assert_eq!(weights.norm_seq_weights.len(), 2);
        assert_eq!(weights.pos_distinct_distrib.len(), 3);
        assert_eq!(weights.std_prob.len(), AA_SIZE);
        assert!(psi_sequence_weights_free(Some(weights)).is_none());
        assert!(psi_internal_pssm_data_free(None).is_none());
    }

    #[test]
    fn translated_psi_purge_aligned_region_matches_c_sequence_discard() {
        let dims = PSIMsaDimensions {
            query_length: 3,
            num_seqs: 2,
        };
        let mut msa = psi_public_msa_new(Some(&dims)).expect("msa");
        for row in &mut msa.data {
            for cell in row {
                *cell = PSIMsaCell {
                    letter: crate::encoding::NCBISTDAA_A,
                    is_aligned: true,
                };
            }
        }

        let mut packed = psi_packed_msa_new(Some(&msa)).expect("packed");
        assert_eq!(
            psi_purge_aligned_region(Some(&mut packed), 1, 1, 2),
            PSI_SUCCESS
        );
        assert_eq!(packed.data[1][1].letter, 0);
        assert!(!packed.data[1][1].is_aligned);
        assert!(packed.use_sequence[1]);

        assert_eq!(
            psi_purge_aligned_region(Some(&mut packed), 1, 0, 3),
            PSI_SUCCESS
        );
        assert!(!packed.use_sequence[1]);
        assert!(packed.use_sequence[2]);

        assert_eq!(
            psi_purge_aligned_region(Some(&mut packed), 0, 0, 1),
            PSIERR_BADPARAM
        );
        assert_eq!(
            psi_purge_aligned_region(Some(&mut packed), 1, 0, 4),
            PSIERR_BADPARAM
        );
        assert_eq!(psi_purge_aligned_region(None, 1, 0, 1), PSIERR_BADPARAM);
    }

    #[test]
    fn translated_psi_purge_similar_alignments_matches_fsm_regions() {
        let dims = PSIMsaDimensions {
            query_length: 5,
            num_seqs: 2,
        };
        let mut msa = psi_public_msa_new(Some(&dims)).expect("msa");
        let a = crate::encoding::NCBISTDAA_A;
        let c = crate::encoding::AMINOACID_TO_NCBISTDAA[b'C' as usize];
        let d = crate::encoding::AMINOACID_TO_NCBISTDAA[b'D' as usize];
        for (row, letters) in [
            vec![a, c, d, a, c],
            vec![a, c, d, a, c],
            vec![a, c, d, a, d],
        ]
        .into_iter()
        .enumerate()
        {
            for (col, letter) in letters.into_iter().enumerate() {
                msa.data[row][col] = PSIMsaCell {
                    letter,
                    is_aligned: true,
                };
            }
        }
        msa.data[1][3].is_aligned = false;
        msa.data[1][4].is_aligned = false;

        let mut packed = psi_packed_msa_new(Some(&msa)).expect("packed");
        s_psi_purge_similar_alignments(&mut packed, 0, 1, K_PSI_IDENTICAL);
        assert!(!packed.data[1][0].is_aligned);
        assert!(!packed.data[1][1].is_aligned);
        assert!(!packed.data[1][2].is_aligned);
        assert!(!packed.use_sequence[1]);
        assert!(packed.use_sequence[2]);

        s_psi_purge_similar_alignments(&mut packed, 0, 2, K_PSI_IDENTICAL);
        assert!(packed.data[2][0].is_aligned);
        assert!(packed.use_sequence[2]);
    }

    #[test]
    fn translated_psi_purge_biased_segments_runs_self_and_near_identical_passes() {
        let dims = PSIMsaDimensions {
            query_length: 4,
            num_seqs: 3,
        };
        let mut msa = psi_public_msa_new(Some(&dims)).expect("msa");
        let a = crate::encoding::NCBISTDAA_A;
        let c = crate::encoding::AMINOACID_TO_NCBISTDAA[b'C' as usize];
        let d = crate::encoding::AMINOACID_TO_NCBISTDAA[b'D' as usize];
        for (row, letters) in [
            vec![a, c, d, a],
            vec![a, c, d, a],
            vec![a, c, d, a],
            vec![a, c, d, c],
        ]
        .into_iter()
        .enumerate()
        {
            for (col, letter) in letters.into_iter().enumerate() {
                msa.data[row][col] = PSIMsaCell {
                    letter,
                    is_aligned: true,
                };
            }
        }

        let mut packed = psi_packed_msa_new(Some(&msa)).expect("packed");
        assert_eq!(psi_purge_biased_segments(Some(&mut packed)), PSI_SUCCESS);
        assert!(!packed.use_sequence[1]);
        assert!(!packed.use_sequence[2]);
        assert!(packed.use_sequence[3]);
        assert_eq!(psi_purge_biased_segments(None), PSIERR_BADPARAM);
    }

    #[test]
    fn translated_debug_print_traits_matches_c_debug_shape() {
        let traits = PsiAlignmentTraits {
            start: 2,
            effective_length: 3,
            n_x_residues: 1,
            n_identical: 2,
        };
        let text = debug_print_traits(&traits, PsiPurgeFsmState::Counting, 5);
        assert!(text.contains("Position: 5 - State: eCounting"));
        assert!(text.contains("effective_length: 3"));
    }

    #[test]
    fn translated_fill_column_probabilities_uses_ncbi_char_order() {
        let dims = PSIMsaDimensions {
            query_length: 1,
            num_seqs: 1,
        };
        let mut weights = psi_sequence_weights_new(Some(&dims), AA_SIZE).expect("weights");
        for (residue, value) in weights.match_weights[0].iter_mut().enumerate() {
            *value = residue as f64 + 0.25;
        }

        let probabilities = s_fill_column_probabilities(&weights, 0).expect("probabilities");
        assert_eq!(probabilities[0], crate::encoding::NCBISTDAA_A as f64 + 0.25);
        assert_eq!(
            probabilities[1],
            crate::encoding::AMINOACID_TO_NCBISTDAA[b'R' as usize] as f64 + 0.25
        );
        assert_eq!(
            probabilities[19],
            crate::encoding::AMINOACID_TO_NCBISTDAA[b'V' as usize] as f64 + 0.25
        );
        assert!(s_fill_column_probabilities(&weights, 1).is_none());
    }

    #[test]
    fn translated_psi_save_diagnostics_copies_requested_payloads() {
        let dims = PSIMsaDimensions {
            query_length: 2,
            num_seqs: 1,
        };
        let aligned = |letter| PsiMsaCell {
            letter,
            is_aligned: true,
            extents: SSeqRange { left: 0, right: 1 },
        };
        let msa = PsiMsa {
            dimensions: dims,
            cell: vec![
                vec![
                    aligned(crate::encoding::NCBISTDAA_A),
                    aligned(crate::encoding::AMINOACID_TO_NCBISTDAA[b'C' as usize]),
                ],
                vec![
                    aligned(crate::encoding::NCBISTDAA_A),
                    aligned(crate::encoding::AMINOACID_TO_NCBISTDAA[b'D' as usize]),
                ],
            ],
            query: vec![
                crate::encoding::NCBISTDAA_A,
                crate::encoding::AMINOACID_TO_NCBISTDAA[b'C' as usize],
            ],
            residue_counts: vec![vec![0, 2, 0], vec![1, 0, 1]],
            alphabet_size: 3,
            num_matching_seqs: vec![2, 1],
        };
        let aligned_block = PsiAlignedBlock {
            pos_extnt: vec![SSeqRange { left: 0, right: 1 }; 2],
            size: vec![2, 1],
        };
        let mut seq_weights = psi_sequence_weights_new(Some(&dims), 3).expect("weights");
        seq_weights.std_prob = vec![0.2, 0.3, 0.5];
        seq_weights.match_weights = vec![vec![0.1, 0.7, 0.2], vec![0.4, 0.2, 0.4]];
        seq_weights.gapless_column_weights = vec![6.0, 7.0];
        seq_weights.sigma = vec![4.0, 3.0];
        seq_weights.independent_observations = vec![8.0, 9.0, 0.0];
        let internal = PsiInternalPssmData {
            ncols: 2,
            nrows: 3,
            pssm: vec![vec![0; 3]; 2],
            scaled_pssm: vec![vec![0; 3]; 2],
            freq_ratios: vec![vec![0.2, 0.6, 0.2], vec![0.1, 0.4, 0.5]],
            pseudocounts: vec![3.0, 4.0],
        };
        let request = PSIDiagnosticsRequest {
            information_content: true,
            residue_frequencies: true,
            weighted_residue_frequencies: true,
            frequency_ratios: true,
            gapless_column_weights: true,
            sigma: true,
            interval_sizes: true,
            num_matching_seqs: true,
            independent_observations: true,
        };
        let mut diagnostics =
            psi_diagnostics_response_new(2, 3, Some(&request)).expect("diagnostics");

        assert_eq!(
            psi_save_diagnostics(
                Some(&msa),
                Some(&aligned_block),
                Some(&seq_weights),
                Some(&internal),
                Some(&mut diagnostics)
            ),
            PSI_SUCCESS
        );

        assert_eq!(
            diagnostics.residue_freqs.as_ref().unwrap()[0],
            vec![0, 2, 0]
        );
        assert_eq!(
            diagnostics.weighted_residue_freqs.as_ref().unwrap()[1],
            vec![0.4, 0.2, 0.4]
        );
        assert_eq!(
            diagnostics.frequency_ratios.as_ref().unwrap()[0],
            vec![0.2, 0.6, 0.2]
        );
        assert_eq!(diagnostics.sigma.as_ref().unwrap(), &vec![4.0, 3.0]);
        assert_eq!(diagnostics.interval_sizes.as_ref().unwrap(), &vec![2, 1]);
        assert_eq!(diagnostics.num_matching_seqs.as_ref().unwrap(), &vec![2, 1]);
        assert_eq!(
            diagnostics.independent_observations.as_ref().unwrap(),
            &vec![8.0, 9.0]
        );
        assert_eq!(diagnostics.gapless_column_weights.as_ref().unwrap()[0], 2.0);
        assert_eq!(diagnostics.gapless_column_weights.as_ref().unwrap()[1], 0.0);
        assert!(diagnostics.information_content.as_ref().unwrap()[0] > 0.0);
    }

    #[test]
    fn translated_psi_save_diagnostics_rejects_bad_inputs_and_dimensions() {
        let dims = PSIMsaDimensions {
            query_length: 1,
            num_seqs: 1,
        };
        let aligned = |letter| PsiMsaCell {
            letter,
            is_aligned: true,
            extents: SSeqRange { left: 0, right: 0 },
        };
        let msa = PsiMsa {
            dimensions: dims,
            cell: vec![vec![aligned(crate::encoding::NCBISTDAA_A)]],
            query: vec![crate::encoding::NCBISTDAA_A],
            residue_counts: vec![vec![1, 0]],
            alphabet_size: 2,
            num_matching_seqs: vec![1],
        };
        let aligned_block = PsiAlignedBlock {
            pos_extnt: vec![SSeqRange { left: 0, right: 0 }],
            size: vec![1],
        };
        let mut seq_weights = psi_sequence_weights_new(Some(&dims), 2).expect("weights");
        seq_weights.std_prob = vec![0.25, 0.75];
        seq_weights.match_weights = vec![vec![1.0, 0.0]];
        let internal = PsiInternalPssmData {
            ncols: 1,
            nrows: 2,
            pssm: vec![vec![0; 2]],
            scaled_pssm: vec![vec![0; 2]],
            freq_ratios: vec![vec![0.5, 0.5]],
            pseudocounts: vec![1.0],
        };
        let request = PSIDiagnosticsRequest {
            information_content: true,
            frequency_ratios: true,
            ..PSIDiagnosticsRequest::default()
        };
        let mut diagnostics =
            psi_diagnostics_response_new(1, 2, Some(&request)).expect("diagnostics");

        assert_eq!(
            psi_save_diagnostics(
                None,
                Some(&aligned_block),
                Some(&seq_weights),
                Some(&internal),
                Some(&mut diagnostics),
            ),
            PSIERR_BADPARAM
        );

        let empty_internal = PsiInternalPssmData {
            freq_ratios: Vec::new(),
            ..internal.clone()
        };
        assert_eq!(
            psi_save_diagnostics(
                Some(&msa),
                Some(&aligned_block),
                Some(&seq_weights),
                Some(&empty_internal),
                Some(&mut diagnostics),
            ),
            PSIERR_BADPARAM
        );

        let mut wrong_dims =
            psi_diagnostics_response_new(2, 2, Some(&request)).expect("diagnostics");
        assert_eq!(
            psi_save_diagnostics(
                Some(&msa),
                Some(&aligned_block),
                Some(&seq_weights),
                Some(&internal),
                Some(&mut wrong_dims),
            ),
            PSIERR_BADPARAM
        );

        let mut short_sources = seq_weights.clone();
        short_sources.sigma.clear();
        let sigma_request = PSIDiagnosticsRequest {
            sigma: true,
            ..PSIDiagnosticsRequest::default()
        };
        let mut sigma_diagnostics =
            psi_diagnostics_response_new(1, 2, Some(&sigma_request)).expect("diagnostics");
        assert_eq!(
            psi_save_diagnostics(
                Some(&msa),
                Some(&aligned_block),
                Some(&short_sources),
                Some(&internal),
                Some(&mut sigma_diagnostics),
            ),
            PSIERR_BADPARAM
        );

        let ratio_only_request = PSIDiagnosticsRequest {
            frequency_ratios: true,
            ..PSIDiagnosticsRequest::default()
        };
        let mut ratio_only_diagnostics =
            psi_diagnostics_response_new(1, 2, Some(&ratio_only_request)).expect("diagnostics");
        assert_eq!(
            psi_save_diagnostics(
                Some(&msa),
                Some(&aligned_block),
                Some(&short_sources),
                Some(&internal),
                Some(&mut ratio_only_diagnostics),
            ),
            PSI_SUCCESS
        );

        let empty_request = PSIDiagnosticsRequest::default();
        let mut no_payload_diagnostics =
            psi_diagnostics_response_new(1, 2, Some(&empty_request)).expect("diagnostics");
        let mut short_unrequested_sources = seq_weights.clone();
        short_unrequested_sources.match_weights.clear();
        short_unrequested_sources.sigma.clear();
        short_unrequested_sources.independent_observations.clear();
        assert_eq!(
            psi_save_diagnostics(
                Some(&msa),
                Some(&aligned_block),
                Some(&short_unrequested_sources),
                Some(&internal),
                Some(&mut no_payload_diagnostics),
            ),
            PSI_SUCCESS
        );

        let mut short_output =
            psi_diagnostics_response_new(1, 2, Some(&request)).expect("diagnostics");
        short_output.frequency_ratios = Some(vec![vec![0.0]]);
        assert_eq!(
            psi_save_diagnostics(
                Some(&msa),
                Some(&aligned_block),
                Some(&seq_weights),
                Some(&internal),
                Some(&mut short_output),
            ),
            PSIERR_BADPARAM
        );
    }

    #[test]
    fn translated_psi_save_cd_diagnostics_copies_cdd_subset() {
        let dims = PSIMsaDimensions {
            query_length: 1,
            num_seqs: 1,
        };
        let cd_msa = PSICdMsa {
            dimensions: Some(dims),
            query: vec![crate::encoding::NCBISTDAA_A],
            msa: vec![vec![PSICdMsaCell {
                is_aligned: true,
                data: Some(PSICdMsaCellData {
                    wfreqs: vec![1.0, 0.0],
                    iobsr: 1.0,
                }),
            }]],
        };
        let mut seq_weights = psi_sequence_weights_new(Some(&dims), 2).expect("weights");
        seq_weights.std_prob = vec![0.25, 0.75];
        seq_weights.match_weights = vec![vec![0.8, 0.2]];
        seq_weights.independent_observations = vec![5.0, 0.0];
        let internal = PsiInternalPssmData {
            ncols: 1,
            nrows: 2,
            pssm: vec![vec![0; 2]],
            scaled_pssm: vec![vec![0; 2]],
            freq_ratios: vec![vec![0.5, 0.5]],
            pseudocounts: vec![1.0],
        };
        let request = PSIDiagnosticsRequest {
            information_content: true,
            weighted_residue_frequencies: true,
            frequency_ratios: true,
            independent_observations: true,
            ..PSIDiagnosticsRequest::default()
        };
        let mut diagnostics =
            psi_diagnostics_response_new(1, 2, Some(&request)).expect("diagnostics");

        assert_eq!(
            psi_save_cd_diagnostics(
                Some(&cd_msa),
                Some(&seq_weights),
                Some(&internal),
                Some(&mut diagnostics)
            ),
            PSI_SUCCESS
        );
        assert_eq!(
            diagnostics.weighted_residue_freqs.as_ref().unwrap()[0],
            vec![0.8, 0.2]
        );
        assert_eq!(
            diagnostics.frequency_ratios.as_ref().unwrap()[0],
            vec![0.5, 0.5]
        );
        assert_eq!(
            diagnostics.independent_observations.as_ref().unwrap(),
            &vec![5.0]
        );
        assert!(diagnostics.information_content.as_ref().unwrap()[0] > 0.0);
        assert_eq!(
            psi_save_cd_diagnostics(
                None,
                Some(&seq_weights),
                Some(&internal),
                Some(&mut diagnostics)
            ),
            PSIERR_BADPARAM
        );

        let empty_internal = PsiInternalPssmData {
            freq_ratios: Vec::new(),
            ..internal.clone()
        };
        assert_eq!(
            psi_save_cd_diagnostics(
                Some(&cd_msa),
                Some(&seq_weights),
                Some(&empty_internal),
                Some(&mut diagnostics),
            ),
            PSIERR_BADPARAM
        );

        let wrong_dims_cd_msa = PSICdMsa {
            dimensions: Some(PSIMsaDimensions {
                query_length: 2,
                num_seqs: 1,
            }),
            ..cd_msa.clone()
        };
        assert_eq!(
            psi_save_cd_diagnostics(
                Some(&wrong_dims_cd_msa),
                Some(&seq_weights),
                Some(&internal),
                Some(&mut diagnostics),
            ),
            PSIERR_BADPARAM
        );

        let mut short_query_cd_msa = cd_msa.clone();
        short_query_cd_msa.query.clear();
        assert_eq!(
            psi_save_cd_diagnostics(
                Some(&short_query_cd_msa),
                Some(&seq_weights),
                Some(&internal),
                Some(&mut diagnostics),
            ),
            PSIERR_BADPARAM
        );

        let mut missing_rows_cd_msa = cd_msa.clone();
        missing_rows_cd_msa.msa.clear();
        assert_eq!(
            psi_save_cd_diagnostics(
                Some(&missing_rows_cd_msa),
                Some(&seq_weights),
                Some(&internal),
                Some(&mut diagnostics),
            ),
            PSIERR_BADPARAM
        );

        let mut short_row_cd_msa = cd_msa.clone();
        short_row_cd_msa.msa[0].clear();
        assert_eq!(
            psi_save_cd_diagnostics(
                Some(&short_row_cd_msa),
                Some(&seq_weights),
                Some(&internal),
                Some(&mut diagnostics),
            ),
            PSIERR_BADPARAM
        );

        let mut short_cd_sources = seq_weights.clone();
        short_cd_sources.match_weights[0].truncate(1);
        let mut cd_source_diagnostics =
            psi_diagnostics_response_new(1, 2, Some(&request)).expect("diagnostics");
        assert_eq!(
            psi_save_cd_diagnostics(
                Some(&cd_msa),
                Some(&short_cd_sources),
                Some(&internal),
                Some(&mut cd_source_diagnostics),
            ),
            PSIERR_BADPARAM
        );

        let cd_ratio_only_request = PSIDiagnosticsRequest {
            frequency_ratios: true,
            ..PSIDiagnosticsRequest::default()
        };
        let mut cd_ratio_only_diagnostics =
            psi_diagnostics_response_new(1, 2, Some(&cd_ratio_only_request)).expect("diagnostics");
        assert_eq!(
            psi_save_cd_diagnostics(
                Some(&cd_msa),
                Some(&short_cd_sources),
                Some(&internal),
                Some(&mut cd_ratio_only_diagnostics),
            ),
            PSI_SUCCESS
        );

        let no_payload_request = PSIDiagnosticsRequest::default();
        let mut no_payload_diagnostics =
            psi_diagnostics_response_new(1, 2, Some(&no_payload_request)).expect("diagnostics");
        let mut short_unrequested_sources = seq_weights.clone();
        short_unrequested_sources.match_weights.clear();
        short_unrequested_sources.independent_observations.clear();
        assert_eq!(
            psi_save_cd_diagnostics(
                Some(&cd_msa),
                Some(&short_unrequested_sources),
                Some(&internal),
                Some(&mut no_payload_diagnostics),
            ),
            PSI_SUCCESS
        );

        let mut short_cd_output =
            psi_diagnostics_response_new(1, 2, Some(&request)).expect("diagnostics");
        short_cd_output.weighted_residue_freqs = Some(vec![vec![0.0]]);
        assert_eq!(
            psi_save_cd_diagnostics(
                Some(&cd_msa),
                Some(&seq_weights),
                Some(&internal),
                Some(&mut short_cd_output),
            ),
            PSIERR_BADPARAM
        );
    }

    #[test]
    fn translated_psi_msa_validation_matches_c_error_order() {
        let dims = PSIMsaDimensions {
            query_length: 3,
            num_seqs: 1,
        };
        let aligned = |letter| PsiMsaCell {
            letter,
            is_aligned: true,
            extents: SSeqRange { left: 0, right: 1 },
        };
        let unaligned = || PsiMsaCell {
            letter: crate::encoding::NCBISTDAA_X,
            is_aligned: false,
            extents: SSeqRange { left: -1, right: 0 },
        };

        let valid = PsiMsa {
            dimensions: dims,
            cell: vec![
                vec![
                    aligned(crate::encoding::NCBISTDAA_A),
                    aligned(crate::encoding::NCBISTDAA_C),
                    aligned(crate::encoding::NCBISTDAA_D),
                ],
                vec![
                    aligned(crate::encoding::NCBISTDAA_A),
                    unaligned(),
                    aligned(crate::encoding::NCBISTDAA_D),
                ],
            ],
            query: vec![
                crate::encoding::NCBISTDAA_A,
                crate::encoding::NCBISTDAA_C,
                crate::encoding::NCBISTDAA_D,
            ],
            residue_counts: vec![vec![0; AA_SIZE]; 3],
            alphabet_size: AA_SIZE as u32,
            num_matching_seqs: vec![1, 1, 1],
        };
        assert_eq!(psi_validate_msa(Some(&valid), false), PSI_SUCCESS);
        assert_eq!(psi_validate_msa_structure_group(Some(&valid)), PSI_SUCCESS);
        assert_eq!(psi_validate_msa(Some(&valid), true), PSI_SUCCESS);

        let mut starting_gap = valid.clone();
        starting_gap.cell[1][0].letter = crate::encoding::NCBISTDAA_GAP;
        assert_eq!(
            psi_validate_msa(Some(&starting_gap), false),
            PSIERR_STARTINGGAP
        );

        let mut ending_gap = valid.clone();
        ending_gap.cell[1][2].letter = crate::encoding::NCBISTDAA_GAP;
        assert_eq!(psi_validate_msa(Some(&ending_gap), false), PSIERR_ENDINGGAP);

        let mut unaligned_column = valid.clone();
        unaligned_column.cell[0][1].is_aligned = false;
        assert_eq!(
            psi_validate_msa(Some(&unaligned_column), false),
            PSIERR_UNALIGNEDCOLUMN
        );
        assert_eq!(psi_validate_msa(Some(&unaligned_column), true), PSI_SUCCESS);

        let mut column_of_gaps = valid.clone();
        column_of_gaps.cell[0][1].letter = crate::encoding::NCBISTDAA_GAP;
        column_of_gaps.cell[1][1].letter = crate::encoding::NCBISTDAA_GAP;
        column_of_gaps.cell[1][1].is_aligned = true;
        assert_eq!(
            psi_validate_msa(Some(&column_of_gaps), false),
            PSIERR_COLUMNOFGAPS
        );

        let mut query_gap = valid.clone();
        query_gap.query[1] = crate::encoding::NCBISTDAA_GAP;
        assert_eq!(psi_validate_msa(Some(&query_gap), true), PSIERR_GAPINQUERY);

        let mut no_aligned_sequences = valid;
        no_aligned_sequences.dimensions.num_seqs = 0;
        assert_eq!(
            psi_validate_msa(Some(&no_aligned_sequences), false),
            PSIERR_NOALIGNEDSEQS
        );
    }

    #[test]
    fn translated_psi_msa_validation_rejects_malformed_shapes() {
        let dims = PSIMsaDimensions {
            query_length: 2,
            num_seqs: 1,
        };
        let aligned = |letter| PsiMsaCell {
            letter,
            is_aligned: true,
            extents: SSeqRange { left: 0, right: 1 },
        };
        let valid = PsiMsa {
            dimensions: dims,
            cell: vec![
                vec![
                    aligned(crate::encoding::NCBISTDAA_A),
                    aligned(crate::encoding::AMINOACID_TO_NCBISTDAA[b'C' as usize]),
                ],
                vec![
                    aligned(crate::encoding::NCBISTDAA_A),
                    aligned(crate::encoding::AMINOACID_TO_NCBISTDAA[b'D' as usize]),
                ],
            ],
            query: vec![
                crate::encoding::NCBISTDAA_A,
                crate::encoding::AMINOACID_TO_NCBISTDAA[b'C' as usize],
            ],
            residue_counts: vec![vec![0; AA_SIZE]; 2],
            alphabet_size: AA_SIZE as u32,
            num_matching_seqs: vec![1, 1],
        };

        let mut short_row = valid.clone();
        short_row.cell[1].truncate(1);
        assert_eq!(psi_validate_msa(Some(&short_row), false), PSIERR_BADPARAM);
        assert_eq!(psi_validate_msa(Some(&short_row), true), PSIERR_BADPARAM);
        assert_eq!(
            psi_validate_msa_structure_group(Some(&short_row)),
            PSIERR_BADPARAM
        );

        let mut short_query = valid;
        short_query.query.truncate(1);
        assert_eq!(psi_validate_msa(Some(&short_query), false), PSIERR_BADPARAM);
        assert_eq!(
            psi_validate_msa_structure_group(Some(&short_query)),
            PSIERR_BADPARAM
        );
    }

    #[test]
    fn translated_psi_structure_group_customization_discards_query_row() {
        let dims = PSIMsaDimensions {
            query_length: 2,
            num_seqs: 1,
        };
        let aligned = |letter| PsiMsaCell {
            letter,
            is_aligned: true,
            extents: SSeqRange { left: 0, right: 1 },
        };
        let mut msa = PsiMsa {
            dimensions: dims,
            cell: vec![
                vec![
                    aligned(crate::encoding::NCBISTDAA_A),
                    aligned(crate::encoding::AMINOACID_TO_NCBISTDAA[b'C' as usize]),
                ],
                vec![
                    aligned(crate::encoding::NCBISTDAA_A),
                    aligned(crate::encoding::AMINOACID_TO_NCBISTDAA[b'D' as usize]),
                ],
            ],
            query: vec![
                crate::encoding::NCBISTDAA_A,
                crate::encoding::AMINOACID_TO_NCBISTDAA[b'C' as usize],
            ],
            residue_counts: vec![vec![0; AA_SIZE]; 2],
            alphabet_size: AA_SIZE as u32,
            num_matching_seqs: vec![0; 2],
        };
        psi_update_position_counts(&mut msa);
        assert_eq!(msa.num_matching_seqs, vec![2, 2]);

        psi_structure_group_customization(&mut msa);
        assert_eq!(msa.cell[0][0].letter, 0);
        assert!(!msa.cell[0][0].is_aligned);
        assert_eq!(msa.num_matching_seqs, vec![1, 1]);
        assert_eq!(
            msa.residue_counts[0][crate::encoding::NCBISTDAA_A as usize],
            1
        );
    }

    #[test]
    fn translated_psi_alignment_blocks_compute_extents_and_x_adjusted_sizes() {
        let dims = PSIMsaDimensions {
            query_length: 4,
            num_seqs: 1,
        };
        let aligned = |letter| PsiMsaCell {
            letter,
            is_aligned: true,
            extents: SSeqRange { left: -1, right: 4 },
        };
        let unaligned = || PsiMsaCell {
            letter: 0,
            is_aligned: false,
            extents: SSeqRange { left: -1, right: 4 },
        };
        let mut msa = PsiMsa {
            dimensions: dims,
            cell: vec![
                vec![
                    aligned(crate::encoding::NCBISTDAA_A),
                    aligned(crate::encoding::AMINOACID_TO_NCBISTDAA[b'X' as usize]),
                    aligned(crate::encoding::AMINOACID_TO_NCBISTDAA[b'D' as usize]),
                    aligned(crate::encoding::AMINOACID_TO_NCBISTDAA[b'E' as usize]),
                ],
                vec![
                    unaligned(),
                    aligned(crate::encoding::AMINOACID_TO_NCBISTDAA[b'C' as usize]),
                    aligned(crate::encoding::AMINOACID_TO_NCBISTDAA[b'D' as usize]),
                    unaligned(),
                ],
            ],
            query: vec![
                crate::encoding::NCBISTDAA_A,
                crate::encoding::AMINOACID_TO_NCBISTDAA[b'X' as usize],
                crate::encoding::AMINOACID_TO_NCBISTDAA[b'D' as usize],
                crate::encoding::AMINOACID_TO_NCBISTDAA[b'E' as usize],
            ],
            residue_counts: vec![vec![0; AA_SIZE]; 4],
            alphabet_size: AA_SIZE as u32,
            num_matching_seqs: vec![0; 4],
        };
        let mut blocks = psi_aligned_block_new(4).expect("blocks");
        assert_eq!(
            psi_compute_alignment_blocks(Some(&mut msa), Some(&mut blocks)),
            PSI_SUCCESS
        );

        assert_eq!(msa.cell[1][1].extents.left, 1);
        assert_eq!(msa.cell[1][1].extents.right, 2);
        assert_eq!(blocks.pos_extnt[1], SSeqRange { left: 1, right: 2 });
        assert_eq!(blocks.pos_extnt[2], SSeqRange { left: 1, right: 2 });
        assert_eq!(blocks.size[2], 1);
        assert_eq!(
            psi_compute_alignment_blocks(None, Some(&mut blocks)),
            PSIERR_BADPARAM
        );
    }

    #[test]
    fn translated_psi_alignment_blocks_reject_malformed_shapes() {
        let dims = PSIMsaDimensions {
            query_length: 2,
            num_seqs: 1,
        };
        let aligned = |letter| PsiMsaCell {
            letter,
            is_aligned: true,
            extents: SSeqRange { left: 0, right: 1 },
        };
        let msa = PsiMsa {
            dimensions: dims,
            cell: vec![
                vec![
                    aligned(crate::encoding::NCBISTDAA_A),
                    aligned(crate::encoding::AMINOACID_TO_NCBISTDAA[b'C' as usize]),
                ],
                vec![
                    aligned(crate::encoding::NCBISTDAA_A),
                    aligned(crate::encoding::AMINOACID_TO_NCBISTDAA[b'D' as usize]),
                ],
            ],
            query: vec![
                crate::encoding::NCBISTDAA_A,
                crate::encoding::AMINOACID_TO_NCBISTDAA[b'C' as usize],
            ],
            residue_counts: vec![vec![0; AA_SIZE]; 2],
            alphabet_size: AA_SIZE as u32,
            num_matching_seqs: vec![0; 2],
        };
        let blocks = psi_aligned_block_new(2).expect("blocks");

        let mut short_row = msa.clone();
        short_row.cell[1].truncate(1);
        let mut blocks_for_short_row = blocks.clone();
        assert_eq!(
            psi_compute_alignment_blocks(Some(&mut short_row), Some(&mut blocks_for_short_row)),
            PSIERR_BADPARAM
        );

        let mut short_query = msa.clone();
        short_query.query.truncate(1);
        let mut blocks_for_short_query = blocks.clone();
        assert_eq!(
            psi_compute_alignment_blocks(Some(&mut short_query), Some(&mut blocks_for_short_query)),
            PSIERR_BADPARAM
        );

        let mut short_blocks = blocks;
        short_blocks.size.truncate(1);
        let mut valid_msa = msa;
        assert_eq!(
            psi_compute_alignment_blocks(Some(&mut valid_msa), Some(&mut short_blocks)),
            PSIERR_BADPARAM
        );
    }

    #[test]
    fn translated_psi_sequence_weight_helpers_match_henikoff_shape() {
        let dims = PSIMsaDimensions {
            query_length: 2,
            num_seqs: 2,
        };
        let aligned = |letter| PsiMsaCell {
            letter,
            is_aligned: true,
            extents: SSeqRange { left: 0, right: 1 },
        };
        let gap = crate::encoding::AMINOACID_TO_NCBISTDAA[b'-' as usize];
        let c = crate::encoding::AMINOACID_TO_NCBISTDAA[b'C' as usize];
        let mut msa = PsiMsa {
            dimensions: dims,
            cell: vec![
                vec![aligned(crate::encoding::NCBISTDAA_A), aligned(c)],
                vec![aligned(crate::encoding::NCBISTDAA_A), aligned(c)],
                vec![aligned(gap), aligned(c)],
            ],
            query: vec![crate::encoding::NCBISTDAA_A, c],
            residue_counts: vec![vec![0; AA_SIZE]; 2],
            alphabet_size: AA_SIZE as u32,
            num_matching_seqs: vec![3, 3],
        };
        psi_update_position_counts(&mut msa);
        let blocks = PsiAlignedBlock {
            pos_extnt: vec![SSeqRange { left: 0, right: 1 }; 2],
            size: vec![2, 2],
        };
        let mut weights = psi_sequence_weights_new(Some(&dims), AA_SIZE).expect("weights");
        let aligned_seqs = psi_get_aligned_sequences_for_position(&msa, 0);
        assert_eq!(aligned_seqs, vec![0, 1, 2]);
        psi_calculate_normalized_sequence_weights(&msa, &blocks, 0, &aligned_seqs, &mut weights);
        let norm_sum: f64 = aligned_seqs
            .iter()
            .map(|&idx| weights.norm_seq_weights[idx as usize])
            .sum();
        assert!((norm_sum - 1.0).abs() < 1.0e-12);

        psi_calculate_match_weights(&msa, 0, &aligned_seqs, &mut weights);
        assert!(weights.match_weights[0][crate::encoding::NCBISTDAA_A as usize] > 0.0);
        assert!(weights.match_weights[0][gap as usize] > 0.0);
        psi_spread_gap_weights(&msa, &mut weights, false);
        assert_eq!(weights.match_weights[0][gap as usize], 0.0);
        let row_sum: f64 = weights.match_weights[0].iter().sum();
        assert!((row_sum - 1.0).abs() < 1.0e-8);
    }

    #[test]
    fn translated_psi_compute_sequence_weights_rejects_malformed_shapes() {
        let dims = PSIMsaDimensions {
            query_length: 2,
            num_seqs: 1,
        };
        let aligned = |letter| PsiMsaCell {
            letter,
            is_aligned: true,
            extents: SSeqRange { left: 0, right: 1 },
        };
        let msa = PsiMsa {
            dimensions: dims,
            cell: vec![
                vec![
                    aligned(crate::encoding::NCBISTDAA_A),
                    aligned(crate::encoding::AMINOACID_TO_NCBISTDAA[b'C' as usize]),
                ],
                vec![
                    aligned(crate::encoding::NCBISTDAA_A),
                    aligned(crate::encoding::AMINOACID_TO_NCBISTDAA[b'D' as usize]),
                ],
            ],
            query: vec![
                crate::encoding::NCBISTDAA_A,
                crate::encoding::AMINOACID_TO_NCBISTDAA[b'C' as usize],
            ],
            residue_counts: vec![vec![0; AA_SIZE]; 2],
            alphabet_size: AA_SIZE as u32,
            num_matching_seqs: vec![2, 2],
        };
        let blocks = PsiAlignedBlock {
            pos_extnt: vec![SSeqRange { left: 0, right: 1 }; 2],
            size: vec![2, 2],
        };
        let mut weights = psi_sequence_weights_new(Some(&dims), AA_SIZE).expect("weights");

        assert_eq!(
            psi_compute_sequence_weights(Some(&msa), Some(&blocks), false, Some(&mut weights)),
            PSI_SUCCESS
        );

        let mut short_msa = msa.clone();
        short_msa.cell[1].truncate(1);
        assert_eq!(
            psi_compute_sequence_weights(
                Some(&short_msa),
                Some(&blocks),
                false,
                Some(&mut weights)
            ),
            PSIERR_BADPARAM
        );

        let mut short_blocks = blocks.clone();
        short_blocks.pos_extnt.truncate(1);
        assert_eq!(
            psi_compute_sequence_weights(
                Some(&msa),
                Some(&short_blocks),
                false,
                Some(&mut weights)
            ),
            PSIERR_BADPARAM
        );

        let mut short_norm = weights.clone();
        short_norm.norm_seq_weights.truncate(1);
        assert_eq!(
            psi_compute_sequence_weights(Some(&msa), Some(&blocks), false, Some(&mut short_norm)),
            PSIERR_BADPARAM
        );

        let mut short_distrib = weights.clone();
        short_distrib.pos_distinct_distrib[0].truncate(EFFECTIVE_ALPHABET);
        assert_eq!(
            psi_compute_sequence_weights(
                Some(&msa),
                Some(&blocks),
                false,
                Some(&mut short_distrib),
            ),
            PSIERR_BADPARAM
        );

        let mut short_match = weights;
        short_match.match_weights[1].truncate(AA_SIZE - 1);
        assert_eq!(
            psi_compute_sequence_weights(Some(&msa), Some(&blocks), false, Some(&mut short_match)),
            PSIERR_BADPARAM
        );
    }

    #[test]
    fn translated_psi_check_sequence_weights_matches_column_sum_rules() {
        let dims = PSIMsaDimensions {
            query_length: 2,
            num_seqs: 2,
        };
        let aligned = |letter| PsiMsaCell {
            letter,
            is_aligned: true,
            extents: SSeqRange { left: 0, right: 1 },
        };
        let x = crate::encoding::AMINOACID_TO_NCBISTDAA[b'X' as usize];
        let c = crate::encoding::AMINOACID_TO_NCBISTDAA[b'C' as usize];
        let mut msa = PsiMsa {
            dimensions: dims,
            cell: vec![
                vec![aligned(crate::encoding::NCBISTDAA_A), aligned(x)],
                vec![aligned(crate::encoding::NCBISTDAA_A), aligned(c)],
                vec![aligned(c), aligned(c)],
            ],
            query: vec![crate::encoding::NCBISTDAA_A, x],
            residue_counts: vec![vec![0; AA_SIZE]; 2],
            alphabet_size: AA_SIZE as u32,
            num_matching_seqs: vec![0; 2],
        };
        psi_update_position_counts(&mut msa);
        let mut weights = psi_sequence_weights_new(Some(&dims), AA_SIZE).expect("weights");
        weights.match_weights[0][crate::encoding::NCBISTDAA_A as usize] = 0.5;
        weights.match_weights[0][c as usize] = 0.5;
        weights.match_weights[1][c as usize] = 7.0;

        assert_eq!(
            psi_check_sequence_weights(Some(&msa), Some(&weights), false),
            0
        );
        weights.match_weights[0][c as usize] = 0.25;
        assert_eq!(
            psi_check_sequence_weights(Some(&msa), Some(&weights), false),
            PSIERR_BADSEQWEIGHTS
        );
        assert_eq!(
            psi_check_sequence_weights(None, Some(&weights), false),
            PSIERR_BADPARAM
        );
    }

    #[test]
    fn translated_psi_compute_frequencies_from_cds_matches_weighted_counts() {
        let dims = PSIMsaDimensions {
            query_length: 2,
            num_seqs: 2,
        };
        let a = crate::encoding::NCBISTDAA_A as usize;
        let c = crate::encoding::AMINOACID_TO_NCBISTDAA[b'C' as usize] as usize;
        let d = crate::encoding::AMINOACID_TO_NCBISTDAA[b'D' as usize] as usize;
        let mut row0_pos0 = vec![0.0; AA_SIZE];
        row0_pos0[a] = 0.5;
        row0_pos0[c] = 0.5;
        let mut row1_pos0 = vec![0.0; AA_SIZE];
        row1_pos0[c] = 1.0;
        let mut row0_pos1 = vec![0.0; AA_SIZE];
        row0_pos1[d] = 1.0;
        let cd_msa = PSICdMsa {
            dimensions: Some(dims),
            query: vec![crate::encoding::NCBISTDAA_A, c as u8],
            msa: vec![
                vec![
                    PSICdMsaCell {
                        is_aligned: true,
                        data: Some(PSICdMsaCellData {
                            wfreqs: row0_pos0,
                            iobsr: 2.0,
                        }),
                    },
                    PSICdMsaCell {
                        is_aligned: true,
                        data: Some(PSICdMsaCellData {
                            wfreqs: row0_pos1,
                            iobsr: 1.0,
                        }),
                    },
                ],
                vec![
                    PSICdMsaCell {
                        is_aligned: true,
                        data: Some(PSICdMsaCellData {
                            wfreqs: row1_pos0,
                            iobsr: 1.0,
                        }),
                    },
                    PSICdMsaCell {
                        is_aligned: false,
                        data: None,
                    },
                ],
            ],
        };
        let sbp =
            crate::stat::blast_score_blk_new(crate::encoding::BLASTAA_SEQ_CODE, 1).expect("sbp");
        let options = crate::options::PSIBlastOptions::default();
        let mut weights = psi_sequence_weights_new(Some(&dims), AA_SIZE).expect("weights");

        assert_eq!(
            psi_compute_frequencies_from_cds(
                Some(&cd_msa),
                Some(&sbp),
                Some(&options),
                Some(&mut weights),
            ),
            PSI_SUCCESS
        );
        assert!((weights.match_weights[0][a] - (1.0 / 3.0)).abs() < 1.0e-12);
        assert!((weights.match_weights[0][c] - (2.0 / 3.0)).abs() < 1.0e-12);
        assert_eq!(weights.match_weights[1][c], 0.5);
        assert_eq!(weights.match_weights[1][d], 0.5);
        assert_eq!(weights.independent_observations[0], 3.0);
        assert_eq!(
            psi_compute_frequencies_from_cds(None, Some(&sbp), Some(&options), Some(&mut weights)),
            PSIERR_BADPARAM
        );
    }

    #[test]
    fn translated_psi_compute_frequencies_from_cds_handles_edges_and_observation_cap() {
        let dims = PSIMsaDimensions {
            query_length: 1,
            num_seqs: 1,
        };
        let a = crate::encoding::NCBISTDAA_A as usize;
        let c = crate::encoding::AMINOACID_TO_NCBISTDAA[b'C' as usize] as usize;
        let sbp =
            crate::stat::blast_score_blk_new(crate::encoding::BLASTAA_SEQ_CODE, 1).expect("sbp");
        let options = crate::options::PSIBlastOptions::default();
        let mut weights = psi_sequence_weights_new(Some(&dims), AA_SIZE).expect("weights");

        let missing_dimensions = PSICdMsa {
            dimensions: None,
            query: vec![crate::encoding::NCBISTDAA_A],
            msa: Vec::new(),
        };
        assert_eq!(
            psi_compute_frequencies_from_cds(
                Some(&missing_dimensions),
                Some(&sbp),
                Some(&options),
                Some(&mut weights),
            ),
            PSIERR_BADPARAM
        );

        let no_sequences = PSICdMsa {
            dimensions: Some(PSIMsaDimensions {
                query_length: 1,
                num_seqs: 0,
            }),
            query: vec![crate::encoding::NCBISTDAA_A],
            msa: Vec::new(),
        };
        weights.match_weights[0][a] = 0.25;
        assert_eq!(
            psi_compute_frequencies_from_cds(
                Some(&no_sequences),
                Some(&sbp),
                Some(&options),
                Some(&mut weights),
            ),
            PSI_SUCCESS
        );
        assert_eq!(weights.match_weights[0][a], 0.25);

        let missing_profile = PSICdMsa {
            dimensions: Some(dims),
            query: vec![crate::encoding::NCBISTDAA_A],
            msa: vec![vec![PSICdMsaCell {
                is_aligned: true,
                data: None,
            }]],
        };
        assert_eq!(
            psi_compute_frequencies_from_cds(
                Some(&missing_profile),
                Some(&sbp),
                Some(&options),
                Some(&mut weights),
            ),
            PSIERR_BADPROFILE
        );

        let mut short_query = missing_profile.clone();
        short_query.query.clear();
        assert_eq!(
            psi_compute_frequencies_from_cds(
                Some(&short_query),
                Some(&sbp),
                Some(&options),
                Some(&mut weights),
            ),
            PSIERR_BADPARAM
        );

        let mut missing_rows = missing_profile.clone();
        missing_rows.msa.clear();
        assert_eq!(
            psi_compute_frequencies_from_cds(
                Some(&missing_rows),
                Some(&sbp),
                Some(&options),
                Some(&mut weights),
            ),
            PSIERR_BADPARAM
        );

        let mut short_row = missing_profile.clone();
        short_row.msa[0].clear();
        assert_eq!(
            psi_compute_frequencies_from_cds(
                Some(&short_row),
                Some(&sbp),
                Some(&options),
                Some(&mut weights),
            ),
            PSIERR_BADPARAM
        );

        let short_profile = PSICdMsa {
            dimensions: Some(dims),
            query: vec![crate::encoding::NCBISTDAA_A],
            msa: vec![vec![PSICdMsaCell {
                is_aligned: true,
                data: Some(PSICdMsaCellData {
                    wfreqs: vec![1.0],
                    iobsr: 1.0,
                }),
            }]],
        };
        assert_eq!(
            psi_compute_frequencies_from_cds(
                Some(&short_profile),
                Some(&sbp),
                Some(&options),
                Some(&mut weights),
            ),
            PSIERR_BADPROFILE
        );

        let mut short_observations = weights.clone();
        short_observations.independent_observations.clear();
        assert_eq!(
            psi_compute_frequencies_from_cds(
                Some(&missing_profile),
                Some(&sbp),
                Some(&options),
                Some(&mut short_observations),
            ),
            PSIERR_BADPARAM
        );

        let mut short_match_weights = weights.clone();
        short_match_weights.match_weights[0].truncate(AA_SIZE - 1);
        assert_eq!(
            psi_compute_frequencies_from_cds(
                Some(&missing_profile),
                Some(&sbp),
                Some(&options),
                Some(&mut short_match_weights),
            ),
            PSIERR_BADPARAM
        );

        let capped = PSICdMsa {
            dimensions: Some(dims),
            query: vec![crate::encoding::NCBISTDAA_A],
            msa: vec![vec![PSICdMsaCell {
                is_aligned: true,
                data: Some(PSICdMsaCellData {
                    wfreqs: {
                        let mut row = vec![0.0; AA_SIZE];
                        row[c] = 1.0;
                        row
                    },
                    iobsr: MAX_IND_OBSERVATIONS as f64 + 50.0,
                }),
            }]],
        };
        assert_eq!(
            psi_compute_frequencies_from_cds(
                Some(&capped),
                Some(&sbp),
                Some(&options),
                Some(&mut weights),
            ),
            PSI_SUCCESS
        );
        assert_eq!(
            weights.independent_observations[0],
            MAX_IND_OBSERVATIONS as f64
        );
        assert!(weights.match_weights[0][c] > weights.match_weights[0][a]);
    }

    #[test]
    fn translated_psi_scale_matrix_brackets_and_updates_lambda() {
        let query = vec![
            crate::encoding::NCBISTDAA_A,
            crate::encoding::AMINOACID_TO_NCBISTDAA[b'C' as usize],
        ];
        let mut sbp =
            crate::stat::blast_score_blk_new(crate::encoding::BLASTAA_SEQ_CODE, 1).expect("sbp");
        sbp.kbp_ideal = Some(crate::stat::KarlinBlk {
            lambda: 0.3176,
            k: 0.134,
            log_k: 0.134_f64.ln(),
            h: 0.401,
            round_down: false,
        });
        sbp.kbp_gap_std[0] = crate::stat::KarlinBlk {
            lambda: 0.267,
            k: 0.041,
            log_k: 0.041_f64.ln(),
            h: 0.14,
            round_down: false,
        };
        let mut scaled = vec![vec![-600; AA_SIZE]; 2];
        scaled[0][crate::encoding::NCBISTDAA_A as usize] = 1000;
        scaled[1][crate::encoding::AMINOACID_TO_NCBISTDAA[b'C' as usize] as usize] = 1000;
        let mut internal = PsiInternalPssmData {
            ncols: 2,
            nrows: AA_SIZE as u32,
            pssm: vec![vec![0; AA_SIZE]; 2],
            scaled_pssm: scaled,
            freq_ratios: vec![vec![0.0; AA_SIZE]; 2],
            pseudocounts: vec![0.0; 2],
        };
        let probs = crate::stat::protein_std_freq_ncbistdaa();

        assert_eq!(
            psi_scale_matrix(
                Some(&query),
                Some(&probs),
                Some(&mut internal),
                Some(&mut sbp)
            ),
            PSI_SUCCESS
        );
        assert!(sbp.kbp_psi[0].lambda > 0.0);
        assert_ne!(internal.pssm[0][crate::encoding::NCBISTDAA_A as usize], 0);
        assert_eq!(
            psi_scale_matrix(None, Some(&probs), Some(&mut internal), Some(&mut sbp)),
            PSIERR_BADPARAM
        );
    }

    #[test]
    fn translated_psi_scale_matrix_rejects_malformed_shapes() {
        let query = vec![
            crate::encoding::NCBISTDAA_A,
            crate::encoding::AMINOACID_TO_NCBISTDAA[b'C' as usize],
        ];
        let mut sbp =
            crate::stat::blast_score_blk_new(crate::encoding::BLASTAA_SEQ_CODE, 1).expect("sbp");
        sbp.kbp_ideal = Some(crate::stat::KarlinBlk {
            lambda: 0.3176,
            k: 0.134,
            log_k: 0.134_f64.ln(),
            h: 0.401,
            round_down: false,
        });
        let probs = crate::stat::protein_std_freq_ncbistdaa();
        let internal = PsiInternalPssmData {
            ncols: 2,
            nrows: AA_SIZE as u32,
            pssm: vec![vec![0; AA_SIZE]; 2],
            scaled_pssm: vec![vec![0; AA_SIZE]; 2],
            freq_ratios: vec![vec![0.0; AA_SIZE]; 2],
            pseudocounts: vec![0.0; 2],
        };

        let mut short_query = internal.clone();
        assert_eq!(
            psi_scale_matrix(
                Some(&query[..1]),
                Some(&probs),
                Some(&mut short_query),
                Some(&mut sbp),
            ),
            PSIERR_BADPARAM
        );

        let mut short_probs = internal.clone();
        assert_eq!(
            psi_scale_matrix(
                Some(&query),
                Some(&probs[..AA_SIZE - 1]),
                Some(&mut short_probs),
                Some(&mut sbp),
            ),
            PSIERR_BADPARAM
        );

        let mut short_pssm = PsiInternalPssmData {
            pssm: vec![vec![0; AA_SIZE], vec![0; AA_SIZE - 1]],
            ..internal.clone()
        };
        assert_eq!(
            psi_scale_matrix(
                Some(&query),
                Some(&probs),
                Some(&mut short_pssm),
                Some(&mut sbp),
            ),
            PSIERR_BADPARAM
        );

        let mut short_scaled = PsiInternalPssmData {
            scaled_pssm: vec![vec![0; AA_SIZE]],
            ..internal
        };
        assert_eq!(
            psi_scale_matrix(
                Some(&query),
                Some(&probs),
                Some(&mut short_scaled),
                Some(&mut sbp),
            ),
            PSIERR_BADPARAM
        );
    }

    #[test]
    fn translated_psi_convert_freq_ratios_to_pssm_rejects_malformed_shapes() {
        let query = vec![crate::encoding::NCBISTDAA_A, 250];
        let mut sbp =
            crate::stat::blast_score_blk_new(crate::encoding::BLASTAA_SEQ_CODE, 1).expect("sbp");
        sbp.name = Some("BLOSUM62".to_string());
        sbp.matrix.data = crate::matrix::BLOSUM62
            .iter()
            .map(|row| row.to_vec())
            .collect();
        sbp.kbp_ideal = Some(crate::stat::KarlinBlk {
            lambda: 0.3176,
            k: 0.134,
            log_k: 0.134_f64.ln(),
            h: 0.401,
            round_down: false,
        });
        let probs = crate::stat::protein_std_freq_ncbistdaa();
        let mut internal = PsiInternalPssmData {
            ncols: 2,
            nrows: AA_SIZE as u32,
            pssm: vec![vec![0; AA_SIZE]; 2],
            scaled_pssm: vec![vec![0; AA_SIZE]; 2],
            freq_ratios: vec![vec![0.0; AA_SIZE]; 2],
            pseudocounts: vec![0.0; 2],
        };

        assert_eq!(
            psi_convert_freq_ratios_to_pssm(
                Some(&mut internal),
                Some(&query),
                Some(&sbp),
                Some(&probs),
            ),
            PSI_SUCCESS
        );
        assert_eq!(
            internal.pssm[1][crate::encoding::NCBISTDAA_A as usize],
            sbp.matrix.data[crate::encoding::NCBISTDAA_X as usize]
                [crate::encoding::NCBISTDAA_A as usize]
        );

        let mut short_query_internal = internal.clone();
        assert_eq!(
            psi_convert_freq_ratios_to_pssm(
                Some(&mut short_query_internal),
                Some(&query[..1]),
                Some(&sbp),
                Some(&probs),
            ),
            PSIERR_BADPARAM
        );

        let mut short_probs_internal = internal.clone();
        assert_eq!(
            psi_convert_freq_ratios_to_pssm(
                Some(&mut short_probs_internal),
                Some(&query),
                Some(&sbp),
                Some(&probs[..AA_SIZE - 1]),
            ),
            PSIERR_BADPARAM
        );

        let mut short_freq_ratios = PsiInternalPssmData {
            freq_ratios: vec![vec![0.0; AA_SIZE]],
            ..internal.clone()
        };
        assert_eq!(
            psi_convert_freq_ratios_to_pssm(
                Some(&mut short_freq_ratios),
                Some(&query),
                Some(&sbp),
                Some(&probs),
            ),
            PSIERR_BADPARAM
        );

        let mut short_scaled = PsiInternalPssmData {
            scaled_pssm: vec![vec![0; AA_SIZE], vec![0; AA_SIZE - 1]],
            ..internal
        };
        assert_eq!(
            psi_convert_freq_ratios_to_pssm(
                Some(&mut short_scaled),
                Some(&query),
                Some(&sbp),
                Some(&probs),
            ),
            PSIERR_BADPARAM
        );
    }

    #[test]
    fn translated_psi_create_and_scale_pssm_from_frequency_ratios_reports_hard_errors() {
        let query = vec![crate::encoding::NCBISTDAA_A];
        let probs = crate::stat::protein_std_freq_ncbistdaa();
        let mut sbp =
            crate::stat::blast_score_blk_new(crate::encoding::BLASTAA_SEQ_CODE, 1).expect("sbp");
        sbp.name = Some("BLOSUM62".to_string());
        sbp.matrix.data = crate::matrix::BLOSUM62
            .iter()
            .map(|row| row.to_vec())
            .collect();
        sbp.kbp_ideal = Some(crate::stat::KarlinBlk {
            lambda: 0.3176,
            k: 0.134,
            log_k: 0.134_f64.ln(),
            h: 0.401,
            round_down: false,
        });
        let internal = PsiInternalPssmData {
            ncols: 1,
            nrows: AA_SIZE as u32,
            pssm: vec![vec![0; AA_SIZE]],
            scaled_pssm: vec![vec![0; AA_SIZE]],
            freq_ratios: vec![vec![0.0; AA_SIZE]],
            pseudocounts: vec![0.0],
        };

        let mut short_query = internal.clone();
        assert_eq!(
            psi_create_and_scale_pssm_from_frequency_ratios(
                Some(&mut short_query),
                Some(&[]),
                Some(&probs),
                Some(&mut sbp),
                crate::options::K_PSSM_NO_IMPALA_SCALING,
            ),
            PSIERR_BADPARAM
        );

        let mut missing_ideal = internal.clone();
        let mut sbp_without_ideal = sbp.clone();
        sbp_without_ideal.kbp_ideal = None;
        assert_eq!(
            psi_create_and_scale_pssm_from_frequency_ratios(
                Some(&mut missing_ideal),
                Some(&query),
                Some(&probs),
                Some(&mut sbp_without_ideal),
                crate::options::K_PSSM_NO_IMPALA_SCALING,
            ),
            PSIERR_BADPARAM
        );

        let mut zero_impala_scale = internal;
        assert_eq!(
            psi_create_and_scale_pssm_from_frequency_ratios(
                Some(&mut zero_impala_scale),
                Some(&query),
                Some(&probs),
                Some(&mut sbp),
                0.0,
            ),
            PSIERR_POSITIVEAVGSCORE
        );
    }

    #[test]
    fn translated_psi_compute_freq_ratios_from_cds_applies_pseudocount_formula() {
        let dims = PSIMsaDimensions {
            query_length: 2,
            num_seqs: 1,
        };
        let a = crate::encoding::NCBISTDAA_A as usize;
        let c = crate::encoding::AMINOACID_TO_NCBISTDAA[b'C' as usize] as usize;
        let x = crate::encoding::AMINOACID_TO_NCBISTDAA[b'X' as usize];
        let cd_msa = PSICdMsa {
            dimensions: Some(dims),
            query: vec![crate::encoding::NCBISTDAA_A, x],
            msa: vec![vec![
                PSICdMsaCell {
                    is_aligned: true,
                    data: Some(PSICdMsaCellData {
                        wfreqs: vec![0.0; AA_SIZE],
                        iobsr: 1.0,
                    }),
                },
                PSICdMsaCell {
                    is_aligned: true,
                    data: Some(PSICdMsaCellData {
                        wfreqs: vec![0.0; AA_SIZE],
                        iobsr: 1.0,
                    }),
                },
            ]],
        };
        let mut seq_weights = psi_sequence_weights_new(Some(&dims), AA_SIZE).expect("weights");
        seq_weights.match_weights[0][a] = 0.75;
        seq_weights.match_weights[0][c] = 0.25;
        seq_weights.independent_observations[0] = 3.0;
        seq_weights.match_weights[1][a] = 1.0;
        seq_weights.independent_observations[1] = 3.0;
        let mut sbp =
            crate::stat::blast_score_blk_new(crate::encoding::BLASTAA_SEQ_CODE, 1).expect("sbp");
        sbp.name = Some("BLOSUM62".to_string());
        sbp.matrix.data = crate::matrix::BLOSUM62
            .iter()
            .map(|row| row.to_vec())
            .collect();
        let mut internal = PsiInternalPssmData {
            ncols: 2,
            nrows: AA_SIZE as u32,
            pssm: vec![vec![0; AA_SIZE]; 2],
            scaled_pssm: vec![vec![0; AA_SIZE]; 2],
            freq_ratios: vec![vec![0.0; AA_SIZE]; 2],
            pseudocounts: vec![0.0; 2],
        };

        assert_eq!(
            psi_compute_freq_ratios_from_cds(
                Some(&cd_msa),
                Some(&seq_weights),
                Some(&sbp),
                5,
                Some(&mut internal),
            ),
            PSI_SUCCESS
        );
        assert!(internal.freq_ratios[0][a] > 0.0);
        assert!(internal.freq_ratios[0][c] > 0.0);
        assert_eq!(internal.freq_ratios[1][a], 0.0);
        assert_eq!(
            psi_compute_freq_ratios_from_cds(
                Some(&cd_msa),
                Some(&seq_weights),
                Some(&sbp),
                -1,
                Some(&mut internal),
            ),
            PSIERR_BADPARAM
        );
    }

    #[test]
    fn translated_psi_compute_freq_ratios_from_cds_rejects_malformed_shapes() {
        let dims = PSIMsaDimensions {
            query_length: 2,
            num_seqs: 1,
        };
        let cd_msa = PSICdMsa {
            dimensions: Some(dims),
            query: vec![
                crate::encoding::NCBISTDAA_A,
                crate::encoding::AMINOACID_TO_NCBISTDAA[b'C' as usize],
            ],
            msa: vec![vec![
                PSICdMsaCell {
                    is_aligned: true,
                    data: Some(PSICdMsaCellData {
                        wfreqs: vec![0.0; AA_SIZE],
                        iobsr: 1.0,
                    }),
                },
                PSICdMsaCell {
                    is_aligned: true,
                    data: Some(PSICdMsaCellData {
                        wfreqs: vec![0.0; AA_SIZE],
                        iobsr: 1.0,
                    }),
                },
            ]],
        };
        let mut seq_weights = psi_sequence_weights_new(Some(&dims), AA_SIZE).expect("weights");
        seq_weights.independent_observations[0] = 3.0;
        seq_weights.independent_observations[1] = 3.0;
        let mut sbp =
            crate::stat::blast_score_blk_new(crate::encoding::BLASTAA_SEQ_CODE, 1).expect("sbp");
        sbp.name = Some("BLOSUM62".to_string());
        let mut internal = PsiInternalPssmData {
            ncols: 2,
            nrows: AA_SIZE as u32,
            pssm: vec![vec![0; AA_SIZE]; 2],
            scaled_pssm: vec![vec![0; AA_SIZE]; 2],
            freq_ratios: vec![vec![0.0; AA_SIZE]; 2],
            pseudocounts: vec![0.0; 2],
        };

        let mut short_query = cd_msa.clone();
        short_query.query.pop();
        assert_eq!(
            psi_compute_freq_ratios_from_cds(
                Some(&short_query),
                Some(&seq_weights),
                Some(&sbp),
                5,
                Some(&mut internal),
            ),
            PSIERR_BADPARAM
        );

        let mut short_observations = seq_weights.clone();
        short_observations.independent_observations.truncate(1);
        assert_eq!(
            psi_compute_freq_ratios_from_cds(
                Some(&cd_msa),
                Some(&short_observations),
                Some(&sbp),
                5,
                Some(&mut internal),
            ),
            PSIERR_BADPARAM
        );

        let mut short_std_prob = seq_weights.clone();
        short_std_prob.std_prob.truncate(AA_SIZE - 1);
        assert_eq!(
            psi_compute_freq_ratios_from_cds(
                Some(&cd_msa),
                Some(&short_std_prob),
                Some(&sbp),
                5,
                Some(&mut internal),
            ),
            PSIERR_BADPARAM
        );

        let mut short_match_weights = seq_weights.clone();
        short_match_weights.match_weights[1].truncate(AA_SIZE - 1);
        assert_eq!(
            psi_compute_freq_ratios_from_cds(
                Some(&cd_msa),
                Some(&short_match_weights),
                Some(&sbp),
                5,
                Some(&mut internal),
            ),
            PSIERR_BADPARAM
        );

        let mut short_internal = PsiInternalPssmData {
            freq_ratios: vec![vec![0.0; AA_SIZE]],
            ..internal.clone()
        };
        assert_eq!(
            psi_compute_freq_ratios_from_cds(
                Some(&cd_msa),
                Some(&seq_weights),
                Some(&sbp),
                5,
                Some(&mut short_internal),
            ),
            PSIERR_BADPARAM
        );
    }

    #[test]
    fn translated_psi_compute_freq_ratios_rejects_malformed_shapes() {
        let dims = PSIMsaDimensions {
            query_length: 2,
            num_seqs: 1,
        };
        let aligned = |letter| PsiMsaCell {
            letter,
            is_aligned: true,
            extents: SSeqRange { left: 0, right: 1 },
        };
        let msa = PsiMsa {
            dimensions: dims,
            cell: vec![
                vec![
                    aligned(crate::encoding::NCBISTDAA_A),
                    aligned(crate::encoding::AMINOACID_TO_NCBISTDAA[b'C' as usize]),
                ],
                vec![
                    aligned(crate::encoding::NCBISTDAA_A),
                    aligned(crate::encoding::AMINOACID_TO_NCBISTDAA[b'D' as usize]),
                ],
            ],
            query: vec![
                crate::encoding::NCBISTDAA_A,
                crate::encoding::AMINOACID_TO_NCBISTDAA[b'C' as usize],
            ],
            residue_counts: vec![vec![0; AA_SIZE]; 2],
            alphabet_size: AA_SIZE as u32,
            num_matching_seqs: vec![2, 2],
        };
        let blocks = PsiAlignedBlock {
            pos_extnt: vec![SSeqRange { left: 0, right: 1 }; 2],
            size: vec![2, 2],
        };
        let mut seq_weights = psi_sequence_weights_new(Some(&dims), AA_SIZE).expect("weights");
        seq_weights.match_weights[0][crate::encoding::NCBISTDAA_A as usize] = 1.0;
        seq_weights.match_weights[1]
            [crate::encoding::AMINOACID_TO_NCBISTDAA[b'C' as usize] as usize] = 1.0;
        seq_weights.pos_num_participating = vec![2, 2];
        let mut sbp =
            crate::stat::blast_score_blk_new(crate::encoding::BLASTAA_SEQ_CODE, 1).expect("sbp");
        sbp.name = Some("BLOSUM62".to_string());
        let mut internal = PsiInternalPssmData {
            ncols: 2,
            nrows: AA_SIZE as u32,
            pssm: vec![vec![0; AA_SIZE]; 2],
            scaled_pssm: vec![vec![0; AA_SIZE]; 2],
            freq_ratios: vec![vec![0.0; AA_SIZE]; 2],
            pseudocounts: vec![0.0; 2],
        };

        let mut short_query_row = msa.clone();
        short_query_row.cell[0].truncate(1);
        assert_eq!(
            psi_compute_freq_ratios(
                Some(&short_query_row),
                Some(&seq_weights),
                Some(&sbp),
                Some(&blocks),
                5,
                false,
                Some(&mut internal),
            ),
            PSIERR_BADPARAM
        );

        let mut short_blocks = blocks.clone();
        short_blocks.size.truncate(1);
        assert_eq!(
            psi_compute_freq_ratios(
                Some(&msa),
                Some(&seq_weights),
                Some(&sbp),
                Some(&short_blocks),
                5,
                false,
                Some(&mut internal),
            ),
            PSIERR_BADPARAM
        );

        let mut short_std_prob = seq_weights.clone();
        short_std_prob.std_prob.truncate(AA_SIZE - 1);
        assert_eq!(
            psi_compute_freq_ratios(
                Some(&msa),
                Some(&short_std_prob),
                Some(&sbp),
                Some(&blocks),
                5,
                false,
                Some(&mut internal),
            ),
            PSIERR_BADPARAM
        );

        let mut short_distrib = seq_weights.clone();
        short_distrib.pos_distinct_distrib[1].truncate(EFFECTIVE_ALPHABET);
        assert_eq!(
            psi_compute_freq_ratios(
                Some(&msa),
                Some(&short_distrib),
                Some(&sbp),
                Some(&blocks),
                5,
                false,
                Some(&mut internal),
            ),
            PSIERR_BADPARAM
        );

        let mut short_match_weights = seq_weights.clone();
        short_match_weights.match_weights[1].truncate(AA_SIZE - 1);
        assert_eq!(
            psi_compute_freq_ratios(
                Some(&msa),
                Some(&short_match_weights),
                Some(&sbp),
                Some(&blocks),
                5,
                false,
                Some(&mut internal),
            ),
            PSIERR_BADPARAM
        );

        let mut short_internal = PsiInternalPssmData {
            pseudocounts: vec![0.0],
            ..internal.clone()
        };
        assert_eq!(
            psi_compute_freq_ratios(
                Some(&msa),
                Some(&seq_weights),
                Some(&sbp),
                Some(&blocks),
                5,
                false,
                Some(&mut short_internal),
            ),
            PSIERR_BADPARAM
        );
    }

    #[test]
    fn translated_psi_cd_msa_validation_checks_profile_data() {
        let dims = PSIMsaDimensions {
            query_length: 2,
            num_seqs: 1,
        };
        let good_cell = PSICdMsaCell {
            is_aligned: true,
            data: Some(PSICdMsaCellData {
                wfreqs: vec![0.25, 0.75],
                iobsr: 1.0,
            }),
        };
        let cd_msa = PSICdMsa {
            dimensions: Some(dims),
            query: vec![crate::encoding::NCBISTDAA_A, crate::encoding::NCBISTDAA_C],
            msa: vec![vec![good_cell.clone(), good_cell]],
        };
        assert_eq!(psi_validate_cd_msa(Some(&cd_msa), 2), PSI_SUCCESS);
        assert_eq!(psi_validate_cd_msa(None, 2), PSIERR_BADPARAM);

        let mut query_gap = cd_msa.clone();
        query_gap.query[0] = crate::encoding::NCBISTDAA_GAP;
        assert_eq!(psi_validate_cd_msa(Some(&query_gap), 2), PSIERR_GAPINQUERY);

        let mut missing_profile = cd_msa.clone();
        missing_profile.msa[0][0].data = None;
        assert_eq!(
            psi_validate_cd_msa(Some(&missing_profile), 2),
            PSIERR_BADPROFILE
        );

        let mut bad_sum = cd_msa;
        bad_sum.msa[0][1].data.as_mut().unwrap().wfreqs = vec![0.1, 0.1];
        assert_eq!(psi_validate_cd_msa(Some(&bad_sum), 2), PSIERR_BADPROFILE);
    }

    #[test]
    fn translated_psi_cd_msa_validation_rejects_malformed_profile_shapes() {
        let dims = PSIMsaDimensions {
            query_length: 1,
            num_seqs: 1,
        };
        let cd_msa = PSICdMsa {
            dimensions: Some(dims),
            query: vec![crate::encoding::NCBISTDAA_A],
            msa: vec![vec![PSICdMsaCell {
                is_aligned: true,
                data: Some(PSICdMsaCellData {
                    wfreqs: vec![0.25, 0.75],
                    iobsr: 1.0,
                }),
            }]],
        };

        let mut missing_dimensions = cd_msa.clone();
        missing_dimensions.dimensions = None;
        assert_eq!(
            psi_validate_cd_msa(Some(&missing_dimensions), 2),
            PSIERR_BADPARAM
        );

        let mut short_query = cd_msa.clone();
        short_query.query.clear();
        assert_eq!(psi_validate_cd_msa(Some(&short_query), 2), PSIERR_BADPARAM);

        let mut missing_rows = cd_msa.clone();
        missing_rows.msa.clear();
        assert_eq!(psi_validate_cd_msa(Some(&missing_rows), 2), PSIERR_BADPARAM);

        let mut short_row = cd_msa.clone();
        short_row.msa[0].clear();
        assert_eq!(psi_validate_cd_msa(Some(&short_row), 2), PSIERR_BADPARAM);

        let mut short_weights = cd_msa.clone();
        short_weights.msa[0][0].data.as_mut().unwrap().wfreqs = vec![1.0];
        assert_eq!(
            psi_validate_cd_msa(Some(&short_weights), 2),
            PSIERR_BADPROFILE
        );

        let mut negative_weight = cd_msa.clone();
        negative_weight.msa[0][0].data.as_mut().unwrap().wfreqs = vec![1.25, -0.25];
        assert_eq!(
            psi_validate_cd_msa(Some(&negative_weight), 2),
            PSIERR_BADPROFILE
        );

        let mut tiny_observation = cd_msa;
        tiny_observation.msa[0][0].data.as_mut().unwrap().iobsr = 0.0;
        assert_eq!(
            psi_validate_cd_msa(Some(&tiny_observation), 2),
            PSIERR_BADPROFILE
        );
    }

    #[test]
    fn translated_psi_debug_print_and_length_without_x_match_c_shape() {
        let dims = PSIMsaDimensions {
            query_length: 3,
            num_seqs: 1,
        };
        let msa = PSIMsa {
            dimensions: dims,
            data: vec![
                vec![
                    PSIMsaCell {
                        letter: crate::encoding::NCBISTDAA_A,
                        is_aligned: true,
                    },
                    PSIMsaCell {
                        letter: crate::encoding::NCBISTDAA_X,
                        is_aligned: false,
                    },
                    PSIMsaCell {
                        letter: crate::encoding::NCBISTDAA_C,
                        is_aligned: true,
                    },
                ],
                vec![
                    PSIMsaCell {
                        letter: crate::encoding::NCBISTDAA_D,
                        is_aligned: true,
                    },
                    PSIMsaCell {
                        letter: crate::encoding::NCBISTDAA_E,
                        is_aligned: true,
                    },
                    PSIMsaCell {
                        letter: crate::encoding::NCBISTDAA_L,
                        is_aligned: true,
                    },
                ],
            ],
        };
        assert_eq!(print_msa_fp(Some(&msa)), "  0\tA.C\n  1\tDEL\n");
        assert_eq!(
            print_msa(Some("ignored"), Some(&msa)),
            print_msa_fp(Some(&msa))
        );
        assert_eq!(
            psi_sequence_length_without_x(&[
                crate::encoding::NCBISTDAA_A,
                crate::encoding::NCBISTDAA_X,
                crate::encoding::NCBISTDAA_C,
            ]),
            2
        );
    }

    #[test]
    fn translated_matrix_frequency_ratios_with_scale_match_psi_name() {
        let ratios = crate::matrix::get_matrix_freq_ratios_with_scale("BLOSUM62").expect("ratios");
        assert_eq!(ratios.bit_scale_factor, 2);
        assert!(
            ratios.data[crate::encoding::NCBISTDAA_A as usize]
                [crate::encoding::NCBISTDAA_A as usize]
                > 1.0
        );
        assert!(crate::matrix::get_matrix_freq_ratios_with_scale("NOT_A_MATRIX").is_none());
    }

    fn psi_public_test_score_block() -> crate::stat::BlastScoreBlk {
        let mut sbp =
            crate::stat::blast_score_blk_new(crate::encoding::BLASTAA_SEQ_CODE, 1).expect("sbp");
        sbp.name = Some("BLOSUM62".to_string());
        sbp.matrix.data = crate::matrix::BLOSUM62
            .iter()
            .map(|row| row.to_vec())
            .collect();
        sbp.kbp_ideal = Some(crate::stat::KarlinBlk {
            lambda: 0.3176,
            k: 0.134,
            log_k: 0.134_f64.ln(),
            h: 0.401,
            round_down: false,
        });
        sbp.kbp_gap_std[0] = crate::stat::KarlinBlk {
            lambda: 0.267,
            k: 0.041,
            log_k: 0.041_f64.ln(),
            h: 0.14,
            round_down: false,
        };
        sbp
    }

    fn psi_public_test_msa() -> PSIMsa {
        let dims = PSIMsaDimensions {
            query_length: 3,
            num_seqs: 2,
        };
        let mut msa = psi_public_msa_new(Some(&dims)).expect("msa");
        let rows = [
            [
                crate::encoding::NCBISTDAA_A,
                crate::encoding::NCBISTDAA_C,
                crate::encoding::NCBISTDAA_D,
            ],
            [
                crate::encoding::NCBISTDAA_A,
                crate::encoding::NCBISTDAA_D,
                crate::encoding::NCBISTDAA_E,
            ],
            [
                crate::encoding::NCBISTDAA_C,
                crate::encoding::NCBISTDAA_D,
                crate::encoding::NCBISTDAA_E,
            ],
        ];
        for (row_idx, row) in rows.iter().enumerate() {
            for (col, &letter) in row.iter().enumerate() {
                msa.data[row_idx][col] = PSIMsaCell {
                    letter,
                    is_aligned: true,
                };
            }
        }
        msa
    }

    #[test]
    fn psi_create_pssm_with_diagnostics_runs_public_c_stage_order() {
        let msa = psi_public_test_msa();
        let options = crate::options::PSIBlastOptions {
            pseudo_count: 5,
            ..crate::options::PSIBlastOptions::default()
        };
        let mut sbp = psi_public_test_score_block();
        let request = psi_diagnostics_request_new_ex(true);
        let mut pssm = None;
        let mut diagnostics = None;

        assert_eq!(
            psi_create_pssm_with_diagnostics(
                Some(&msa),
                Some(&options),
                Some(&mut sbp),
                Some(&request),
                Some(&mut pssm),
                Some(&mut diagnostics),
            ),
            PSI_SUCCESS
        );

        let pssm = pssm.expect("pssm");
        assert_eq!(pssm.ncols, 3);
        assert_eq!(pssm.nrows, AA_SIZE as u32);
        assert!(pssm.pssm[0][crate::encoding::NCBISTDAA_A as usize] > crate::stat::BLAST_SCORE_MIN);
        assert!(sbp.kbp_psi[0].lambda > 0.0);

        let diagnostics = diagnostics.expect("diagnostics");
        assert_eq!(diagnostics.query_length, 3);
        assert!(
            diagnostics.frequency_ratios.as_ref().unwrap()[0]
                [crate::encoding::NCBISTDAA_A as usize]
                > 0.0
        );
        assert!(diagnostics.information_content.as_ref().unwrap()[0].is_finite());
        assert!(diagnostics.gapless_column_weights.as_ref().unwrap()[0] > 0.0);
    }

    #[test]
    fn psi_create_pssm_with_diagnostics_populates_requested_payloads() {
        let msa = psi_public_test_msa();
        let options = crate::options::PSIBlastOptions {
            pseudo_count: 5,
            ..crate::options::PSIBlastOptions::default()
        };
        let mut sbp = psi_public_test_score_block();
        let request = PSIDiagnosticsRequest {
            information_content: true,
            residue_frequencies: true,
            weighted_residue_frequencies: true,
            frequency_ratios: true,
            gapless_column_weights: true,
            sigma: true,
            interval_sizes: true,
            num_matching_seqs: true,
            independent_observations: true,
        };
        let mut pssm = None;
        let mut diagnostics = None;

        assert_eq!(
            psi_create_pssm_with_diagnostics(
                Some(&msa),
                Some(&options),
                Some(&mut sbp),
                Some(&request),
                Some(&mut pssm),
                Some(&mut diagnostics),
            ),
            PSI_SUCCESS
        );

        let diagnostics = diagnostics.expect("diagnostics");
        assert_eq!(
            diagnostics.residue_freqs.as_ref().unwrap()[0][crate::encoding::NCBISTDAA_A as usize],
            2
        );
        assert!(
            diagnostics.weighted_residue_freqs.as_ref().unwrap()[0]
                [crate::encoding::NCBISTDAA_A as usize]
                > 0.0
        );
        assert!(
            diagnostics.frequency_ratios.as_ref().unwrap()[0]
                [crate::encoding::NCBISTDAA_A as usize]
                > 0.0
        );
        assert!(diagnostics.information_content.as_ref().unwrap()[0].is_finite());
        assert!(diagnostics.gapless_column_weights.as_ref().unwrap()[0] > 0.0);
        assert!(diagnostics.sigma.as_ref().unwrap()[0] > 0.0);
        assert_eq!(diagnostics.interval_sizes.as_ref().unwrap()[0], 3);
        assert_eq!(diagnostics.num_matching_seqs.as_ref().unwrap()[0], 3);
        let independent_observations = diagnostics.independent_observations.as_ref().unwrap();
        assert_eq!(independent_observations.len(), 3);
        assert!(independent_observations
            .iter()
            .all(|value| value.is_finite()));
        assert!(pssm.expect("pssm").pssm[0][crate::encoding::NCBISTDAA_A as usize] > 0);
    }

    #[test]
    fn psi_create_pssm_with_diagnostics_ignores_request_without_diagnostics_output() {
        let msa = psi_public_test_msa();
        let options = crate::options::PSIBlastOptions {
            pseudo_count: 5,
            ..crate::options::PSIBlastOptions::default()
        };
        let mut sbp = psi_public_test_score_block();
        let request = psi_diagnostics_request_new_ex(true);
        let mut pssm = None;

        assert_eq!(
            psi_create_pssm_with_diagnostics(
                Some(&msa),
                Some(&options),
                Some(&mut sbp),
                Some(&request),
                Some(&mut pssm),
                None,
            ),
            PSI_SUCCESS
        );

        let pssm = pssm.expect("pssm");
        assert_eq!(pssm.ncols, 3);
        assert_eq!(pssm.nrows, AA_SIZE as u32);
    }

    #[test]
    fn psi_create_pssm_with_diagnostics_clears_diagnostics_output_without_request() {
        let msa = psi_public_test_msa();
        let options = crate::options::PSIBlastOptions {
            pseudo_count: 5,
            ..crate::options::PSIBlastOptions::default()
        };
        let mut sbp = psi_public_test_score_block();
        let request = psi_diagnostics_request_new_ex(true);
        let mut pssm = None;
        let mut diagnostics =
            Some(psi_diagnostics_response_new(1, AA_SIZE as u32, Some(&request)).unwrap());

        assert_eq!(
            psi_create_pssm_with_diagnostics(
                Some(&msa),
                Some(&options),
                Some(&mut sbp),
                None,
                Some(&mut pssm),
                Some(&mut diagnostics),
            ),
            PSI_SUCCESS
        );

        let pssm = pssm.expect("pssm");
        assert_eq!(pssm.ncols, 3);
        assert!(diagnostics.is_none());
    }

    #[test]
    fn psi_create_pssm_with_diagnostics_allows_no_payload_diagnostics_request() {
        let msa = psi_public_test_msa();
        let options = crate::options::PSIBlastOptions {
            pseudo_count: 5,
            ..crate::options::PSIBlastOptions::default()
        };
        let mut sbp = psi_public_test_score_block();
        let request = psi_diagnostics_request_new();
        let mut pssm = None;
        let mut diagnostics = None;

        assert_eq!(
            psi_create_pssm_with_diagnostics(
                Some(&msa),
                Some(&options),
                Some(&mut sbp),
                Some(&request),
                Some(&mut pssm),
                Some(&mut diagnostics),
            ),
            PSI_SUCCESS
        );

        let pssm = pssm.expect("pssm");
        assert_eq!(pssm.ncols, 3);
        let diagnostics = diagnostics.expect("diagnostics");
        assert_eq!(diagnostics.query_length, 3);
        assert_eq!(diagnostics.alphabet_size, AA_SIZE as u32);
        assert!(diagnostics.information_content.is_none());
        assert!(diagnostics.residue_freqs.is_none());
        assert!(diagnostics.weighted_residue_freqs.is_none());
        assert!(diagnostics.frequency_ratios.is_none());
        assert!(diagnostics.gapless_column_weights.is_none());
        assert!(diagnostics.sigma.is_none());
        assert!(diagnostics.interval_sizes.is_none());
        assert!(diagnostics.num_matching_seqs.is_none());
        assert!(diagnostics.independent_observations.is_none());
    }

    #[test]
    fn psi_create_pssm_with_diagnostics_rejects_missing_pssm_output_pointer() {
        let msa = psi_public_test_msa();
        let options = crate::options::PSIBlastOptions {
            pseudo_count: 5,
            ..crate::options::PSIBlastOptions::default()
        };
        let mut sbp = psi_public_test_score_block();
        let request = psi_diagnostics_request_new_ex(true);
        let mut diagnostics =
            Some(psi_diagnostics_response_new(1, AA_SIZE as u32, Some(&request)).unwrap());

        assert_eq!(
            psi_create_pssm_with_diagnostics(
                Some(&msa),
                Some(&options),
                Some(&mut sbp),
                Some(&request),
                None,
                Some(&mut diagnostics),
            ),
            PSIERR_BADPARAM
        );
        assert!(diagnostics.is_some());
    }

    #[test]
    fn psi_create_pssm_with_diagnostics_cleans_outputs_on_validation_error() {
        let mut msa = psi_public_test_msa();
        msa.data[0][1].letter = crate::encoding::NCBISTDAA_GAP;
        let options = crate::options::PSIBlastOptions::default();
        let mut sbp = psi_public_test_score_block();
        let request = psi_diagnostics_request_new_ex(true);
        let mut pssm = Some(psi_matrix_new(1, AA_SIZE as u32).expect("preexisting pssm"));
        let mut diagnostics =
            Some(psi_diagnostics_response_new(1, AA_SIZE as u32, Some(&request)).unwrap());

        assert_eq!(
            psi_create_pssm_with_diagnostics(
                Some(&msa),
                Some(&options),
                Some(&mut sbp),
                Some(&request),
                Some(&mut pssm),
                Some(&mut diagnostics),
            ),
            PSIERR_GAPINQUERY
        );
        assert!(pssm.is_none());
        assert!(diagnostics.is_none());
        assert_eq!(
            psi_create_pssm_with_diagnostics(
                None,
                Some(&options),
                Some(&mut sbp),
                None,
                Some(&mut pssm),
                None,
            ),
            PSIERR_BADPARAM
        );
    }

    #[test]
    fn psi_create_pssm_with_diagnostics_cleans_outputs_on_scaling_hard_error() {
        let msa = psi_public_test_msa();
        let options = crate::options::PSIBlastOptions {
            pseudo_count: 5,
            ..crate::options::PSIBlastOptions::default()
        };
        let mut sbp = psi_public_test_score_block();
        sbp.kbp_ideal = None;
        let request = psi_diagnostics_request_new_ex(true);
        let mut pssm = Some(psi_matrix_new(1, AA_SIZE as u32).expect("preexisting pssm"));
        let mut diagnostics =
            Some(psi_diagnostics_response_new(1, AA_SIZE as u32, Some(&request)).unwrap());

        assert_eq!(
            psi_create_pssm_with_diagnostics(
                Some(&msa),
                Some(&options),
                Some(&mut sbp),
                Some(&request),
                Some(&mut pssm),
                Some(&mut diagnostics),
            ),
            PSIERR_BADPARAM
        );
        assert!(pssm.is_none());
        assert!(diagnostics.is_none());
    }

    #[test]
    fn psi_create_pssm_with_diagnostics_cleans_outputs_on_frequency_ratio_error() {
        let msa = psi_public_test_msa();
        let options = crate::options::PSIBlastOptions {
            pseudo_count: -1,
            ..crate::options::PSIBlastOptions::default()
        };
        let mut sbp = psi_public_test_score_block();
        let request = psi_diagnostics_request_new_ex(true);
        let mut pssm = Some(psi_matrix_new(1, AA_SIZE as u32).expect("preexisting pssm"));
        let mut diagnostics =
            Some(psi_diagnostics_response_new(1, AA_SIZE as u32, Some(&request)).unwrap());

        assert_eq!(
            psi_create_pssm_with_diagnostics(
                Some(&msa),
                Some(&options),
                Some(&mut sbp),
                Some(&request),
                Some(&mut pssm),
                Some(&mut diagnostics),
            ),
            PSIERR_BADPARAM
        );
        assert!(pssm.is_none());
        assert!(diagnostics.is_none());
    }

    #[test]
    fn psi_create_pssm_with_diagnostics_cleans_outputs_on_bad_alphabet_size() {
        let msa = psi_public_test_msa();
        let options = crate::options::PSIBlastOptions {
            pseudo_count: 5,
            ..crate::options::PSIBlastOptions::default()
        };
        let mut sbp = psi_public_test_score_block();
        sbp.alphabet_size = crate::encoding::BLASTAA_SIZE + 1;
        let request = psi_diagnostics_request_new_ex(true);
        let mut pssm = Some(psi_matrix_new(1, AA_SIZE as u32).expect("preexisting pssm"));
        let mut diagnostics =
            Some(psi_diagnostics_response_new(1, AA_SIZE as u32, Some(&request)).unwrap());

        assert_eq!(
            psi_create_pssm_with_diagnostics(
                Some(&msa),
                Some(&options),
                Some(&mut sbp),
                Some(&request),
                Some(&mut pssm),
                Some(&mut diagnostics),
            ),
            PSIERR_BADPARAM
        );
        assert!(pssm.is_none());
        assert!(diagnostics.is_none());
    }

    #[test]
    fn psi_create_pssm_with_diagnostics_preserves_outputs_on_required_input_error() {
        let msa = psi_public_test_msa();
        let options = crate::options::PSIBlastOptions {
            pseudo_count: 5,
            ..crate::options::PSIBlastOptions::default()
        };
        let mut sbp = psi_public_test_score_block();
        let request = psi_diagnostics_request_new_ex(true);
        let mut pssm = Some(psi_matrix_new(1, AA_SIZE as u32).expect("preexisting pssm"));
        let mut diagnostics =
            Some(psi_diagnostics_response_new(1, AA_SIZE as u32, Some(&request)).unwrap());

        assert_eq!(
            psi_create_pssm_with_diagnostics(
                None,
                Some(&options),
                Some(&mut sbp),
                Some(&request),
                Some(&mut pssm),
                Some(&mut diagnostics),
            ),
            PSIERR_BADPARAM
        );
        assert!(pssm.is_some());
        assert!(diagnostics.is_some());

        assert_eq!(
            psi_create_pssm_with_diagnostics(
                Some(&msa),
                None,
                Some(&mut sbp),
                Some(&request),
                Some(&mut pssm),
                Some(&mut diagnostics),
            ),
            PSIERR_BADPARAM
        );
        assert!(pssm.is_some());
        assert!(diagnostics.is_some());

        assert_eq!(
            psi_create_pssm_with_diagnostics(
                Some(&msa),
                Some(&options),
                None,
                Some(&request),
                Some(&mut pssm),
                Some(&mut diagnostics),
            ),
            PSIERR_BADPARAM
        );
        assert!(pssm.is_some());
        assert!(diagnostics.is_some());
    }

    #[test]
    fn psi_create_pssm_with_diagnostics_reports_public_msa_validation_errors() {
        let options = crate::options::PSIBlastOptions::default();
        let request = psi_diagnostics_request_new_ex(true);
        let assert_public_validation_error = |msa: &PSIMsa, expected_status: i32| {
            let mut sbp = psi_public_test_score_block();
            let mut pssm = Some(psi_matrix_new(1, AA_SIZE as u32).expect("preexisting pssm"));
            let mut diagnostics =
                Some(psi_diagnostics_response_new(1, AA_SIZE as u32, Some(&request)).unwrap());
            assert_eq!(
                psi_create_pssm_with_diagnostics(
                    Some(msa),
                    Some(&options),
                    Some(&mut sbp),
                    Some(&request),
                    Some(&mut pssm),
                    Some(&mut diagnostics),
                ),
                expected_status
            );
            assert!(pssm.is_none());
            assert!(diagnostics.is_none());
        };

        let mut unaligned_column = psi_public_test_msa();
        for row in &mut unaligned_column.data {
            row[1].is_aligned = false;
        }
        assert_public_validation_error(&unaligned_column, PSIERR_UNALIGNEDCOLUMN);

        let mut starting_gap = psi_public_test_msa();
        starting_gap.data[1][0].letter = crate::encoding::NCBISTDAA_GAP;
        assert_public_validation_error(&starting_gap, PSIERR_STARTINGGAP);

        let mut ending_gap = psi_public_test_msa();
        ending_gap.data[1][2].letter = crate::encoding::NCBISTDAA_GAP;
        assert_public_validation_error(&ending_gap, PSIERR_ENDINGGAP);

        let mut column_of_gaps = psi_public_test_msa();
        for row in &mut column_of_gaps.data {
            row[1].letter = crate::encoding::NCBISTDAA_GAP;
            row[1].is_aligned = true;
        }
        assert_public_validation_error(&column_of_gaps, PSIERR_COLUMNOFGAPS);

        let mut no_aligned = psi_public_test_msa();
        no_aligned.dimensions.num_seqs = 0;
        no_aligned.data.truncate(1);
        let mut sbp = psi_public_test_score_block();
        let mut pssm = None;
        assert_eq!(
            psi_create_pssm_with_diagnostics(
                Some(&no_aligned),
                Some(&options),
                Some(&mut sbp),
                None,
                Some(&mut pssm),
                None,
            ),
            PSIERR_NOALIGNEDSEQS
        );
        assert!(pssm.is_none());
    }

    #[test]
    fn psi_create_pssm_with_diagnostics_rejects_malformed_public_msa_shapes() {
        let options = crate::options::PSIBlastOptions::default();
        let request = psi_diagnostics_request_new_ex(true);
        let assert_bad_public_msa = |msa: &PSIMsa| {
            let mut sbp = psi_public_test_score_block();
            let mut pssm = Some(psi_matrix_new(1, AA_SIZE as u32).expect("preexisting pssm"));
            let mut diagnostics =
                Some(psi_diagnostics_response_new(1, AA_SIZE as u32, Some(&request)).unwrap());

            assert_eq!(
                psi_create_pssm_with_diagnostics(
                    Some(msa),
                    Some(&options),
                    Some(&mut sbp),
                    Some(&request),
                    Some(&mut pssm),
                    Some(&mut diagnostics),
                ),
                PSIERR_BADPARAM
            );
            assert!(pssm.is_none());
            assert!(diagnostics.is_none());
        };

        let mut missing_rows = psi_public_test_msa();
        missing_rows.data.truncate(1);
        assert_bad_public_msa(&missing_rows);

        let mut short_query = psi_public_test_msa();
        short_query.data[0].truncate(2);
        assert_bad_public_msa(&short_query);

        let mut short_aligned_row = psi_public_test_msa();
        short_aligned_row.data[1].truncate(2);
        assert_bad_public_msa(&short_aligned_row);

        let invalid_letter = crate::encoding::BLASTAA_SIZE as u8 + 1;
        let mut invalid_query_letter = psi_public_test_msa();
        invalid_query_letter.data[0][0].letter = invalid_letter;
        assert_bad_public_msa(&invalid_query_letter);

        let mut invalid_aligned_letter = psi_public_test_msa();
        invalid_aligned_letter.data[1][1].letter = invalid_letter;
        assert_bad_public_msa(&invalid_aligned_letter);
    }

    #[test]
    fn psi_create_pssm_with_diagnostics_structure_group_accepts_query_unaligned_column() {
        let mut msa = psi_public_test_msa();
        msa.data[0][1].is_aligned = false;
        let options = crate::options::PSIBlastOptions {
            nsg_compatibility_mode: true,
            pseudo_count: 5,
            ..crate::options::PSIBlastOptions::default()
        };
        let mut sbp = psi_public_test_score_block();
        let mut pssm = None;

        assert_eq!(
            psi_create_pssm_with_diagnostics(
                Some(&msa),
                Some(&options),
                Some(&mut sbp),
                None,
                Some(&mut pssm),
                None,
            ),
            PSI_SUCCESS
        );
        assert_eq!(pssm.expect("pssm").ncols, 3);
    }

    #[test]
    fn psi_create_pssm_with_diagnostics_ignore_unaligned_positions_accepts_unaligned_column() {
        let mut msa = psi_public_test_msa();
        for row in &mut msa.data {
            row[1].is_aligned = false;
        }
        let mut sbp = psi_public_test_score_block();
        let mut pssm = Some(psi_matrix_new(1, AA_SIZE as u32).expect("preexisting pssm"));
        let mut diagnostics = None;
        let request = PSIDiagnosticsRequest {
            frequency_ratios: true,
            interval_sizes: true,
            num_matching_seqs: true,
            ..psi_diagnostics_request_new()
        };
        let strict_options = crate::options::PSIBlastOptions {
            pseudo_count: 5,
            ..crate::options::PSIBlastOptions::default()
        };

        assert_eq!(
            psi_create_pssm_with_diagnostics(
                Some(&msa),
                Some(&strict_options),
                Some(&mut sbp),
                Some(&request),
                Some(&mut pssm),
                Some(&mut diagnostics),
            ),
            PSIERR_UNALIGNEDCOLUMN
        );
        assert!(pssm.is_none());
        assert!(diagnostics.is_none());

        let ignore_options = crate::options::PSIBlastOptions {
            ignore_unaligned_positions: true,
            pseudo_count: 5,
            ..crate::options::PSIBlastOptions::default()
        };
        let mut sbp = psi_public_test_score_block();
        let mut pssm = None;
        let mut diagnostics = None;
        assert_eq!(
            psi_create_pssm_with_diagnostics(
                Some(&msa),
                Some(&ignore_options),
                Some(&mut sbp),
                Some(&request),
                Some(&mut pssm),
                Some(&mut diagnostics),
            ),
            PSI_SUCCESS
        );

        let pssm = pssm.expect("pssm");
        assert_eq!(pssm.ncols, 3);
        assert_eq!(pssm.nrows, AA_SIZE as u32);
        let diagnostics = diagnostics.expect("diagnostics");
        assert_eq!(diagnostics.num_matching_seqs.as_ref().unwrap()[1], 0);
        assert!(diagnostics.interval_sizes.as_ref().unwrap()[1] >= 3);
        assert!(diagnostics.frequency_ratios.as_ref().unwrap()[1]
            .iter()
            .all(|ratio| ratio.is_finite()));

        let mut query_gap = msa;
        query_gap.data[0][1].letter = crate::encoding::NCBISTDAA_GAP;
        let mut sbp = psi_public_test_score_block();
        let mut pssm = Some(psi_matrix_new(1, AA_SIZE as u32).expect("preexisting pssm"));
        let mut diagnostics =
            Some(psi_diagnostics_response_new(1, AA_SIZE as u32, Some(&request)).unwrap());
        assert_eq!(
            psi_create_pssm_with_diagnostics(
                Some(&query_gap),
                Some(&ignore_options),
                Some(&mut sbp),
                Some(&request),
                Some(&mut pssm),
                Some(&mut diagnostics),
            ),
            PSIERR_GAPINQUERY
        );
        assert!(pssm.is_none());
        assert!(diagnostics.is_none());
    }
}
