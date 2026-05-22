//! Rust equivalent of blast_parameters.c — computed search parameters.
//! These are derived from options + scoring blocks at search time.

use crate::link_hsps::LinkHSPParameters;
use crate::options::*;
use crate::pattern::SphiQueryInfo;
use crate::program::{self, ProgramType};
use crate::queryinfo::QueryInfo;
use crate::stat::{self, BlastScoreBlk, GumbelBlk, KarlinBlk};
use std::fmt::Write as _;

#[derive(Debug, Clone, Copy)]
struct BlastParameterScoreBlock<'a> {
    kbp: &'a [KarlinBlk],
    kbp_gap: Option<&'a [KarlinBlk]>,
    kbp_std: Option<&'a [KarlinBlk]>,
    gbp: Option<&'a GumbelBlk>,
    scale_factor: f64,
    reward: i32,
    penalty: i32,
    matrix_only_scoring: bool,
}

impl<'a> BlastParameterScoreBlock<'a> {
    /// blast-rs: Score-block view helper for Rust parameter calculation; not a
    /// direct NCBI C port.
    fn valid_kbp(&self, context: usize) -> Option<&KarlinBlk> {
        self.kbp
            .get(context)
            .filter(|kbp| stat::s_blast_karlin_blk_is_valid(Some(kbp)))
    }

    /// blast-rs: Selects the Karlin block array for update calculations; not a
    /// direct NCBI C port.
    fn update_kbp_array(&self) -> Option<(&[KarlinBlk], bool)> {
        if let Some(kbp_gap) = self.kbp_gap {
            Some((kbp_gap, true))
        } else {
            self.kbp_std.map(|kbp_std| (kbp_std, false))
        }
    }

    /// blast-rs: Selects the Karlin block array for hit-saving calculations;
    /// not a direct NCBI C port.
    fn hit_saving_kbp_array(&self) -> Option<(&[KarlinBlk], bool)> {
        if let Some(kbp_gap) = self.kbp_gap {
            Some((kbp_gap, true))
        } else if !self.kbp.is_empty() {
            Some((self.kbp, false))
        } else {
            None
        }
    }
}

impl<'a> From<&'a BlastScoreBlk> for BlastParameterScoreBlock<'a> {
    /// blast-rs: Builds a Rust view over `BlastScoreBlk`; not a direct NCBI C port.
    fn from(sbp: &'a BlastScoreBlk) -> Self {
        Self {
            kbp: &sbp.kbp,
            kbp_gap: (!sbp.kbp_gap.is_empty()).then_some(sbp.kbp_gap.as_slice()),
            kbp_std: (!sbp.kbp_std.is_empty()).then_some(sbp.kbp_std.as_slice()),
            gbp: sbp.gbp.as_ref(),
            scale_factor: sbp.scale_factor,
            reward: sbp.reward,
            penalty: sbp.penalty,
            matrix_only_scoring: sbp.matrix_only_scoring,
        }
    }
}

#[derive(Debug, Clone, Copy, Default)]
pub struct BlastGappedCutoffs {
    pub cutoff_score: i32,
    pub cutoff_score_max: i32,
}

#[derive(Debug, Clone, Copy, Default)]
pub struct BlastUngappedCutoffs {
    pub x_dropoff_init: i32,
    pub cutoff_score: i32,
}

/// Scoring parameters computed from options and score block.
#[derive(Debug, Clone)]
pub struct ScoringParameters {
    pub options: ScoringOptions,
    pub reward: i32,
    pub penalty: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub scale_factor: f64,
}

impl ScoringParameters {
    /// blast-rs: Native constructor from public scoring options; not a direct
    /// NCBI C port.
    pub fn from_options(opts: &ScoringOptions, scale_factor: f64) -> Self {
        ScoringParameters {
            options: opts.clone(),
            reward: opts.reward,
            penalty: opts.penalty,
            gap_open: opts.gap_open,
            gap_extend: opts.gap_extend,
            scale_factor,
        }
    }
}

/// Extension parameters computed from options and score block.
#[derive(Debug, Clone)]
pub struct ExtensionParameters {
    pub options: ExtensionOptions,
    pub gap_x_dropoff: i32,       // raw score x-dropoff for preliminary
    pub gap_x_dropoff_final: i32, // raw score x-dropoff for traceback
    pub gap_trigger: i32,         // raw score trigger for gapped extension
}

/// Hit saving parameters computed from options.
#[derive(Debug, Clone)]
pub struct HitSavingParameters {
    pub options: HitSavingOptions,
    pub cutoff_score_min: i32,
    pub low_score: Vec<i32>,
    pub cutoffs: Vec<BlastGappedCutoffs>,
    pub link_hsp_params: Option<LinkHSPParameters>,
    pub prelim_evalue: f64,
}

/// Initial word parameters computed from options and score block.
#[derive(Debug, Clone)]
pub struct InitialWordParameters {
    pub options: InitialWordOptions,
    pub x_dropoff_max: i32,
    pub cutoff_score_min: i32,
    pub cutoffs: Vec<BlastUngappedCutoffs>,
    pub ungapped_extension: bool,
    /// Scoring table for 4-base nucleotide words (256 entries).
    pub nucl_score_table: [i32; 256],
}

impl InitialWordParameters {
    /// blast-rs: Nucleotide packed-byte score-table builder for Rust lookup
    /// parameters; not a direct NCBI C port.
    ///
    /// Build the nucleotide score table for 4-base packed words.
    /// Indexed by XOR of query and subject packed bytes.
    /// Each bit pair in the XOR result: 00 = match (reward), != 00 = mismatch (penalty).
    pub fn build_nucl_score_table(reward: i32, penalty: i32) -> [i32; 256] {
        let mut table = [0i32; 256];
        for xor_val in 0..256u32 {
            let mut score = 0i32;
            for pos in 0..4 {
                let bits = (xor_val >> (6 - 2 * pos)) & 3;
                score += if bits == 0 { reward } else { penalty };
            }
            table[xor_val as usize] = score;
        }
        table
    }
}

/// blast-rs: Local scale-factor helper matching the integer conversion used by
/// BLAST parameter ports; not a standalone NCBI C function.
fn scaled_i32(value: i32, scale_factor: f64) -> i32 {
    value.saturating_mul(scale_factor as i32)
}

/// NCBI: s_GetCutoffEvalue (blast_parameters.c:134).
pub fn s_get_cutoff_evalue(program_number: ProgramType) -> f64 {
    match program_number {
        program::MAPPING | program::BLASTN | program::PHI_BLASTN => stat::CUTOFF_E_BLASTN,
        program::BLASTP | program::RPS_BLAST | program::PHI_BLASTP => stat::CUTOFF_E_BLASTP,
        program::BLASTX => stat::CUTOFF_E_BLASTX,
        program::TBLASTN | program::PSI_TBLASTN | program::RPS_TBLASTN => stat::CUTOFF_E_TBLASTN,
        program::TBLASTX => stat::CUTOFF_E_TBLASTX,
        _ => 0.0,
    }
}

/// NCBI: s_GetEstimatedPhiExpect (blast_parameters.c:659).
///
/// C reads `paramC`/`Lambda` from `sbp->kbp[0]` and pattern probability
/// from `query_info->pattern_info`. Rust keeps the PHI pattern block
/// separate from [`QueryInfo`], so those score-block values are explicit.
pub fn s_get_estimated_phi_expect(
    score: i32,
    query_info: Option<&QueryInfo>,
    pattern_info: Option<&SphiQueryInfo>,
    lambda: f64,
    param_c: f64,
    effective_num_patterns: i32,
) -> f64 {
    let (Some(query_info), Some(pattern_info)) = (query_info, pattern_info) else {
        return 0.0;
    };
    let Some(context) = query_info.contexts.first() else {
        return 0.0;
    };
    let score = score as f64;
    context.eff_searchsp as f64
        * param_c
        * (1.0 + lambda * score)
        * effective_num_patterns as f64
        * pattern_info.probability
        * (-lambda * score).exp()
}

/// NCBI: s_PhiBlastCutoffScore (blast_parameters.c:689).
pub fn s_phi_blast_cutoff_score(
    ethresh: f64,
    query_info: Option<&QueryInfo>,
    pattern_info: Option<&SphiQueryInfo>,
    lambda: f64,
    param_c: f64,
) -> i32 {
    let (Some(query_info), Some(pattern_info)) = (query_info, pattern_info) else {
        return 1;
    };
    let Some(context) = query_info.contexts.first() else {
        return 1;
    };
    let occurrence_offsets: Vec<i32> = pattern_info
        .occurrences
        .iter()
        .map(|occurrence| occurrence.offset)
        .collect();
    let effective_num_patterns = crate::hspstream::phi_blast_get_effective_number_of_patterns(
        &occurrence_offsets,
        context.length_adjustment,
    );

    let mut low_score = 1;
    let mut high_score = 100;
    for _ in 0..20 {
        let target_score = (low_score + high_score) / 2;
        let expect = s_get_estimated_phi_expect(
            target_score,
            Some(query_info),
            Some(pattern_info),
            lambda,
            param_c,
            effective_num_patterns,
        );
        if expect > ethresh {
            low_score = target_score;
        } else {
            high_score = target_score;
        }
        if high_score - low_score <= 1 {
            break;
        }
    }
    low_score
}

/// NCBI: BlastScoringParametersNew (blast_parameters.c:542).
pub fn blast_scoring_parameters_new(
    score_options: Option<&ScoringOptions>,
    sbp: Option<&BlastScoreBlk>,
    parameters: &mut Option<ScoringParameters>,
) -> i16 {
    let Some(score_options) = score_options else {
        return 1;
    };
    let scale_factor = sbp.map_or(1.0, |sbp| sbp.scale_factor);
    *parameters = Some(ScoringParameters {
        options: score_options.clone(),
        reward: score_options.reward,
        penalty: score_options.penalty,
        gap_open: scaled_i32(score_options.gap_open, scale_factor),
        gap_extend: scaled_i32(score_options.gap_extend, scale_factor),
        scale_factor,
    });
    0
}

/// blast-rs: Nullable pointer-shaped adapter for
/// [`blast_scoring_parameters_new`]; not a direct NCBI C port.
pub fn blast_scoring_parameters_new_c(
    score_options: Option<&ScoringOptions>,
    sbp: Option<&BlastScoreBlk>,
    parameters: Option<&mut Option<ScoringParameters>>,
) -> i16 {
    let Some(parameters) = parameters else {
        return 1;
    };
    blast_scoring_parameters_new(score_options, sbp, parameters)
}

/// blast-rs: Rust spelling of the effective-length parameters constructor; not
/// a direct NCBI C port.
pub fn blast_effective_lengths_parameters_new(
    options: &EffectiveLengthsOptions,
    db_length: i64,
    num_seqs: i32,
    parameters: &mut Option<EffectiveLengthsParameters>,
) -> i16 {
    *parameters = Some(EffectiveLengthsParameters {
        options: options.clone(),
        real_db_length: db_length,
        real_num_seqs: num_seqs,
    });
    0
}

/// blast-rs: Nullable pointer-shaped adapter for
/// [`blast_effective_lengths_parameters_new`]; not a direct NCBI C port.
pub fn blast_effective_lengths_parameters_new_c(
    options: Option<&EffectiveLengthsOptions>,
    db_length: i64,
    num_seqs: i32,
    parameters: Option<&mut Option<EffectiveLengthsParameters>>,
) -> i16 {
    let Some(parameters) = parameters else {
        return 1;
    };
    let Some(options) = options else {
        return 1;
    };
    blast_effective_lengths_parameters_new(options, db_length, num_seqs, parameters)
}

/// NCBI: BlastLinkHSPParametersNew (blast_parameters.c:598).
pub fn blast_link_hsp_parameters_new(
    program_number: ProgramType,
    gapped_calculation: bool,
    link_hsp_params: &mut Option<LinkHSPParameters>,
) -> i16 {
    let mut params = LinkHSPParameters::default();
    if program_number == program::BLASTN || !gapped_calculation {
        params.gap_prob = stat::BLAST_GAP_PROB;
        params.gap_decay_rate = stat::BLAST_GAP_DECAY_RATE;
    } else {
        params.gap_prob = stat::BLAST_GAP_PROB_GAPPED;
        params.gap_decay_rate = stat::BLAST_GAP_DECAY_RATE_GAPPED;
    }
    params.gap_size = stat::BLAST_GAP_SIZE;
    params.overlap_size = stat::BLAST_OVERLAP_SIZE;
    *link_hsp_params = Some(params);
    0
}

/// blast-rs: Nullable pointer-shaped adapter for
/// [`blast_link_hsp_parameters_new`]; not a direct NCBI C port.
pub fn blast_link_hsp_parameters_new_c(
    program_number: ProgramType,
    gapped_calculation: bool,
    link_hsp_params: Option<&mut Option<LinkHSPParameters>>,
) -> i16 {
    let Some(link_hsp_params) = link_hsp_params else {
        return -1;
    };
    blast_link_hsp_parameters_new(program_number, gapped_calculation, link_hsp_params)
}

/// NCBI: BlastLinkHSPParametersFree (blast_parameters.c:591).
pub fn blast_link_hsp_parameters_free(
    parameters: &mut Option<LinkHSPParameters>,
) -> Option<LinkHSPParameters> {
    *parameters = None;
    None
}

/// NCBI: BlastLinkHSPParametersUpdate (blast_parameters.c:625).
pub fn blast_link_hsp_parameters_update(
    word_params: Option<&InitialWordParameters>,
    hit_params: Option<&mut HitSavingParameters>,
    _gapped_calculation: bool,
) -> i16 {
    let (Some(word_params), Some(hit_params)) = (word_params, hit_params) else {
        return -1;
    };
    if let Some(link_hsp_params) = hit_params.link_hsp_params.as_mut() {
        link_hsp_params.cutoff_small_gap = word_params.cutoff_score_min;
    }
    0
}

/// blast-rs: Nullable pointer-shaped adapter for
/// [`blast_link_hsp_parameters_update`]; not a direct NCBI C port.
pub fn blast_link_hsp_parameters_update_c(
    word_params: Option<&InitialWordParameters>,
    hit_params: Option<&mut HitSavingParameters>,
    gapped_calculation: bool,
) -> i16 {
    blast_link_hsp_parameters_update(word_params, hit_params, gapped_calculation)
}

/// NCBI: CalculateLinkHSPCutoffs (blast_parameters.c:998).
pub fn calculate_link_hsp_cutoffs(
    program: ProgramType,
    query_info: &QueryInfo,
    sbp: &BlastScoreBlk,
    link_hsp_params: Option<&mut LinkHSPParameters>,
    word_params: &InitialWordParameters,
    mut db_length: i64,
    mut subject_length: i32,
) {
    let Some(link_hsp_params) = link_hsp_params else {
        return;
    };
    let Some((_, kbp)) = stat::s_blast_find_smallest_lambda(&sbp.kbp, query_info) else {
        return;
    };
    if query_info.contexts.is_empty() {
        return;
    }

    let window_size = link_hsp_params.gap_size + link_hsp_params.overlap_size + 1;
    let gap_prob = stat::BLAST_GAP_PROB;
    link_hsp_params.gap_prob = gap_prob;
    let gap_decay_rate = link_hsp_params.gap_decay_rate;

    let Some(last_context) = query_info.contexts.iter().rposition(|ctx| ctx.is_valid) else {
        return;
    };
    let last = &query_info.contexts[last_context];
    let mut query_length = (last.query_offset + last.query_length - 1) / (last_context as i32 + 1);

    if program::blast_subject_is_translated(program) || program == program::RPS_TBLASTN {
        subject_length /= crate::util::CODON_LENGTH as i32;
        db_length /= crate::util::CODON_LENGTH as i64;
    }

    let expected_length =
        crate::math::blast_nint((kbp.k * query_length as f64 * subject_length as f64).ln() / kbp.h)
            as i32;
    query_length = (query_length - expected_length).max(1);
    subject_length = (subject_length - expected_length).max(1);

    let y_variable = if db_length > subject_length as i64 {
        (db_length as f64 / subject_length as f64).ln() * kbp.k / gap_decay_rate
    } else {
        ((subject_length + expected_length) as f64 / subject_length as f64).ln() * kbp.k
            / gap_decay_rate
    };
    let search_sp = query_length as i64 * subject_length as i64;
    let mut x_variable = 0.25 * y_variable * search_sp as f64;
    const EPSILON: f64 = 1.0e-9;

    if search_sp > 8 * window_size as i64 * window_size as i64 {
        x_variable /= 1.0 - gap_prob + EPSILON;
        link_hsp_params.cutoff_big_gap = (x_variable.ln() / kbp.lambda).floor() as i32 + 1;
        x_variable = y_variable * (window_size * window_size) as f64;
        x_variable /= gap_prob + EPSILON;
        link_hsp_params.cutoff_small_gap = word_params
            .cutoff_score_min
            .max((x_variable.ln() / kbp.lambda).floor() as i32 + 1);
    } else {
        link_hsp_params.cutoff_big_gap = (x_variable.ln() / kbp.lambda).floor() as i32 + 1;
        link_hsp_params.gap_prob = 0.0;
        link_hsp_params.cutoff_small_gap = 0;
    }

    let scale_factor = sbp.scale_factor as i32;
    link_hsp_params.cutoff_big_gap *= scale_factor;
    link_hsp_params.cutoff_small_gap *= scale_factor;
}

/// blast-rs: Rust output-builder equivalent of NCBI debug
/// `printBlastScoringParameters`; not a direct C port.
pub fn print_blast_scoring_parameters(params: Option<&ScoringParameters>) -> String {
    let Some(params) = params else {
        return "parameters{ null }\n".to_string();
    };
    let mut out = String::new();
    let _ = writeln!(out, "BlastScoringParameters:");
    let _ = writeln!(out, "  options:");
    let _ = writeln!(
        out,
        "    matrix = {}",
        params.options.matrix_name.as_deref().unwrap_or("")
    );
    let _ = writeln!(out, "    reward = {}", params.options.reward);
    let _ = writeln!(out, "    penalty = {}", params.options.penalty);
    let _ = writeln!(
        out,
        "    gapped_calculation = {}",
        params.options.gapped_calculation as i32
    );
    let _ = writeln!(out, "    gap_open = {}", params.options.gap_open);
    let _ = writeln!(out, "    gap_extend = {}", params.options.gap_extend);
    let _ = writeln!(out, "  reward = {}", params.reward);
    let _ = writeln!(out, "  penalty = {}", params.penalty);
    let _ = writeln!(out, "  gap_open = {}", params.gap_open);
    let _ = writeln!(out, "  gap_extend = {}", params.gap_extend);
    let _ = writeln!(out, "  scale_factor = {}\n", params.scale_factor);
    out
}

/// blast-rs: Rust output-builder equivalent of NCBI debug
/// `printBlastInitialWordParamters`; not a direct C port.
pub fn print_blast_initial_word_paramters(
    word_params: &InitialWordParameters,
    query_info: &QueryInfo,
) -> String {
    let mut out = String::new();
    let _ = writeln!(out, "BlastInitialWordParamters:");
    let _ = writeln!(out, "  x_dropoff_max = {}", word_params.x_dropoff_max);
    let _ = writeln!(out, "  cutoff_score_min = {}", word_params.cutoff_score_min);
    let _ = writeln!(
        out,
        "  ungapped_extension = {}",
        word_params.ungapped_extension as i32
    );
    let _ = writeln!(out, "  cutoffs:");
    for (context, info) in query_info.contexts.iter().enumerate() {
        if !info.is_valid {
            continue;
        }
        if let Some(cutoffs) = word_params.cutoffs.get(context) {
            let _ = writeln!(
                out,
                "    {} x_dropoff_init = {}",
                context, cutoffs.x_dropoff_init
            );
            let _ = writeln!(
                out,
                "    {} cutoff_score = {}",
                context, cutoffs.cutoff_score
            );
        }
    }
    out
}

/// blast-rs: Rust output-builder equivalent of NCBI debug
/// `printBlastExtensionParameters`; not a direct C port.
pub fn print_blast_extension_parameters(ext_params: &ExtensionParameters) -> String {
    let mut out = String::new();
    let _ = writeln!(out, "BlastExtensionParameters:");
    let _ = writeln!(out, "  gap_x_dropoff = {}", ext_params.gap_x_dropoff);
    let _ = writeln!(
        out,
        "  gap_x_dropoff_final = {}",
        ext_params.gap_x_dropoff_final
    );
    out
}

/// blast-rs: Rust output-builder equivalent of NCBI debug
/// `printBlastHitSavingParameters`; not a direct C port.
pub fn print_blast_hit_saving_parameters(
    hit_params: &HitSavingParameters,
    query_info: &QueryInfo,
) -> String {
    let mut out = String::new();
    let _ = writeln!(out, "BlastHitSavingParameters:");
    let _ = writeln!(out, "  cutoff_score_min = {}", hit_params.cutoff_score_min);
    for (context, info) in query_info.contexts.iter().enumerate() {
        if !info.is_valid {
            continue;
        }
        if let Some(cutoffs) = hit_params.cutoffs.get(context) {
            let _ = writeln!(
                out,
                "    {} cutoff_score = {}",
                context, cutoffs.cutoff_score
            );
            let _ = writeln!(
                out,
                "    {} cutoff_score_max = {}",
                context, cutoffs.cutoff_score_max
            );
        }
    }
    out
}

/// blast-rs: Rust output-builder equivalent of NCBI debug `printAllParameters`;
/// not a direct C port.
pub fn print_all_parameters(
    hit_params: &HitSavingParameters,
    ext_params: &ExtensionParameters,
    word_params: &InitialWordParameters,
    query_info: &QueryInfo,
) -> String {
    let mut out = print_blast_initial_word_paramters(word_params, query_info);
    out.push_str(&print_blast_extension_parameters(ext_params));
    out.push_str(&print_blast_hit_saving_parameters(hit_params, query_info));
    out
}

/// NCBI: s_FillReturnCutoffsInfo (blast_parameters.c).
pub fn s_fill_return_cutoffs_info(
    return_cutoffs: Option<&mut crate::diagnostics::RawCutoffs>,
    _score_params: Option<&ScoringParameters>,
    word_params: Option<&InitialWordParameters>,
    ext_params: Option<&ExtensionParameters>,
    hit_params: Option<&HitSavingParameters>,
) -> i16 {
    let Some(return_cutoffs) = return_cutoffs else {
        return crate::util::BLASTERR_INVALIDPARAM;
    };
    if let Some(word_params) = word_params {
        return_cutoffs.x_drop_ungapped = word_params.x_dropoff_max;
        return_cutoffs.ungapped_cutoff = word_params.cutoff_score_min;
    }
    if let Some(ext_params) = ext_params {
        return_cutoffs.x_drop_gap = ext_params.gap_x_dropoff;
        return_cutoffs.x_drop_gap_final = ext_params.gap_x_dropoff_final;
    }
    if let Some(hit_params) = hit_params {
        return_cutoffs.cutoff_score = hit_params.cutoff_score_min;
    }
    0
}

/// NCBI: BlastExtensionParametersNew (blast_parameters.c:422).
pub fn blast_extension_parameters_new(
    program_number: ProgramType,
    options: &ExtensionOptions,
    sbp: &BlastScoreBlk,
    query_info: &QueryInfo,
    parameters: &mut Option<ExtensionParameters>,
) -> i16 {
    if stat::s_blast_find_valid_karlin_blk(&sbp.kbp, query_info).is_err() {
        *parameters = None;
        return -1;
    }
    let sbp_view = BlastParameterScoreBlock::from(sbp);

    let mut params = ExtensionParameters {
        options: options.clone(),
        gap_x_dropoff: 0,
        gap_x_dropoff_final: 0,
        gap_trigger: options.gap_trigger as i32,
    };

    if let Some(kbp_gap) = sbp_view.kbp_gap {
        let Some((min_lambda, _)) = stat::s_blast_find_smallest_lambda(kbp_gap, query_info) else {
            *parameters = None;
            return -1;
        };
        params.gap_x_dropoff =
            (options.gap_x_dropoff * crate::math::NCBIMATH_LN2 / min_lambda) as i32;
        params.gap_x_dropoff_final = ((options.gap_x_dropoff_final * crate::math::NCBIMATH_LN2
            / min_lambda)
            .max(params.gap_x_dropoff as f64)) as i32;
    }

    if sbp_view.scale_factor > 1.0 {
        params.gap_x_dropoff = scaled_i32(params.gap_x_dropoff, sbp_view.scale_factor);
        params.gap_x_dropoff_final = scaled_i32(params.gap_x_dropoff_final, sbp_view.scale_factor);
    }

    if (program_number == program::BLASTN && sbp_view.matrix_only_scoring)
        || program_number == program::MAPPING
    {
        params.gap_x_dropoff = options.gap_x_dropoff as i32;
        if program_number == program::BLASTN && sbp_view.matrix_only_scoring {
            params.gap_x_dropoff_final = options.gap_x_dropoff_final as i32;
        }
    }

    *parameters = Some(params);
    0
}

/// blast-rs: Nullable pointer-shaped adapter for
/// [`blast_extension_parameters_new`]; not a direct NCBI C port.
pub fn blast_extension_parameters_new_c(
    program_number: ProgramType,
    options: Option<&ExtensionOptions>,
    sbp: Option<&BlastScoreBlk>,
    query_info: Option<&QueryInfo>,
    parameters: Option<&mut Option<ExtensionParameters>>,
) -> i16 {
    let Some(parameters) = parameters else {
        return 0;
    };
    let (Some(options), Some(sbp), Some(query_info)) = (options, sbp, query_info) else {
        *parameters = None;
        return -1;
    };
    blast_extension_parameters_new(program_number, options, sbp, query_info, parameters)
}

/// NCBI: BlastHitSavingParametersUpdate (blast_parameters.c:831).
pub fn blast_hit_saving_parameters_update(
    program_number: ProgramType,
    sbp: &BlastScoreBlk,
    query_info: &QueryInfo,
    avg_subject_length: i32,
    composition_based_stats: i32,
    params: &mut HitSavingParameters,
) -> i16 {
    blast_hit_saving_parameters_update_c(
        program_number,
        Some(sbp),
        Some(query_info),
        avg_subject_length,
        composition_based_stats,
        Some(params),
    )
}

/// blast-rs: Nullable pointer-shaped adapter for
/// [`blast_hit_saving_parameters_update`]; not a direct NCBI C port.
pub fn blast_hit_saving_parameters_update_c(
    program_number: ProgramType,
    sbp: Option<&BlastScoreBlk>,
    query_info: Option<&QueryInfo>,
    avg_subject_length: i32,
    composition_based_stats: i32,
    params: Option<&mut HitSavingParameters>,
) -> i16 {
    let (Some(sbp), Some(query_info), Some(params)) = (sbp, query_info, params) else {
        return -1;
    };
    let sbp_view = BlastParameterScoreBlock::from(sbp);
    let Some((kbp_array, gapped_calculation)) = sbp_view.hit_saving_kbp_array() else {
        return -1;
    };
    let scale_factor = sbp_view.scale_factor;
    params.prelim_evalue = params.options.expect_value;
    if params.cutoffs.len() < query_info.num_contexts() {
        params
            .cutoffs
            .resize(query_info.num_contexts(), BlastGappedCutoffs::default());
    }

    if params.options.cutoff_score > 0 {
        let new_cutoff = scaled_i32(params.options.cutoff_score, scale_factor);
        for cutoff in &mut params.cutoffs {
            if program_number == program::BLASTN && sbp_view.matrix_only_scoring {
                cutoff.cutoff_score = params.options.cutoff_score;
                cutoff.cutoff_score_max = params.options.cutoff_score / 2;
            } else {
                cutoff.cutoff_score = new_cutoff;
                cutoff.cutoff_score_max = new_cutoff;
            }
        }
        params.cutoff_score_min = new_cutoff;
        return 0;
    }

    let mut cutoff_min = i32::MAX;
    for (context, context_info) in query_info.contexts.iter().enumerate() {
        if !context_info.is_valid {
            params.cutoffs[context].cutoff_score = i32::MAX;
            continue;
        }
        let Some(kbp) = kbp_array
            .get(context)
            .filter(|kbp| stat::s_blast_karlin_blk_is_valid(Some(kbp)))
        else {
            return -1;
        };
        let mut evalue = params.options.expect_value;
        let mut searchsp = context_info.eff_searchsp;
        if searchsp <= 0 {
            searchsp = (context_info.query_length.max(1) as i64)
                .saturating_mul(avg_subject_length.max(1) as i64);
        }
        if program_number == program::RPS_TBLASTN {
            searchsp /= crate::util::NUM_FRAMES as i64;
        }

        let new_cutoff =
            if let Some(gbp) = sbp_view.gbp.filter(|gbp| stat::gumbel_blk_is_filled(gbp)) {
                let cbs_stretch = if composition_based_stats > 1 {
                    5.0
                } else {
                    1.0
                };
                params.prelim_evalue = cbs_stretch * evalue;
                stat::spouge_etos(
                    cbs_stretch * evalue,
                    kbp,
                    gbp,
                    context_info.query_length,
                    avg_subject_length,
                )
            } else {
                let (cutoff, adjusted_evalue) =
                    stat::blast_cutoffs(1, evalue, kbp, searchsp as f64, false, 0.0);
                evalue = adjusted_evalue;
                params.prelim_evalue = evalue;
                cutoff
            };
        let _ = evalue; // suppress unused-modify; sum-stats path below uses its own.

        let new_cutoff = scaled_i32(new_cutoff, scale_factor);
        params.cutoffs[context].cutoff_score = new_cutoff;
        params.cutoffs[context].cutoff_score_max = new_cutoff;
        cutoff_min = cutoff_min.min(new_cutoff);
    }

    // NCBI `blast_parameters.c:950-976`: sum-statistics cutoff is computed
    // AFTER the per-context loop, using an average query length, the
    // sum-stats-specific `searchsp = MIN(avg_qlen, avg_subj) * avg_subj`,
    // and evalue = 1.0 (constant, NOT options->expect_value). Then
    // `cutoff[context] = MIN(new_cutoff, cutoff[context])`. Previously
    // our port inlined a buggy version inside the per-context loop with
    // wrong evalue and wrong searchsp.
    if params.link_hsp_params.is_some() && gapped_calculation && !query_info.contexts.is_empty() {
        let decay = params
            .link_hsp_params
            .as_ref()
            .map_or(0.0, |link| link.gap_decay_rate);
        let evalue_hsp = 1.0;
        // NCBI uses `query_info->last_context` — the index of the last
        // valid context.
        let Some(last_idx) = query_info.contexts.iter().rposition(|ctx| ctx.is_valid) else {
            return -1;
        };
        let last_ctx = &query_info.contexts[last_idx];
        let concat_qlen = last_ctx.query_offset + last_ctx.query_length;
        let avg_qlen = concat_qlen / (last_idx as i32 + 1);
        // NCBI's `MIN(avg_qlen, avg_subj) * avg_subj` — intentional
        // clamp for sum-stats search space (different from the simple
        // `avg_qlen * avg_subj` product used in `BLAST_Cutoffs`).
        let sum_searchsp = (avg_qlen.min(avg_subject_length).max(1) as i64)
            .saturating_mul(avg_subject_length.max(1) as i64);

        for (context, context_info) in query_info.contexts.iter().enumerate() {
            if !context_info.is_valid {
                continue;
            }
            let Some(kbp) = kbp_array
                .get(context)
                .filter(|kbp| stat::s_blast_karlin_blk_is_valid(Some(kbp)))
            else {
                return -1;
            };
            let (sum_cutoff, _) =
                stat::blast_cutoffs(1, evalue_hsp, kbp, sum_searchsp as f64, true, decay);
            let sum_cutoff = scaled_i32(sum_cutoff, scale_factor);
            let current = params.cutoffs[context].cutoff_score;
            let merged = sum_cutoff.min(current);
            params.cutoffs[context].cutoff_score = merged;
            cutoff_min = cutoff_min.min(merged);
        }
    }
    params.cutoff_score_min = cutoff_min;
    0
}

/// NCBI: BlastHitSavingParametersNew (blast_parameters.c:734).
pub fn blast_hit_saving_parameters_new(
    program_number: ProgramType,
    options: &HitSavingOptions,
    sbp: &BlastScoreBlk,
    query_info: &QueryInfo,
    avg_subject_length: i32,
    composition_based_stats: i32,
    parameters: &mut Option<HitSavingParameters>,
) -> i16 {
    let sbp_view = BlastParameterScoreBlock::from(sbp);
    let gapped_calculation = sbp_view.kbp_gap.is_some();
    if options.do_sum_stats && gapped_calculation && avg_subject_length <= 0 {
        *parameters = None;
        return 1;
    }

    let mut params_options = options.clone();
    let mut link_hsp_params = None;
    if options.do_sum_stats {
        blast_link_hsp_parameters_new(program_number, gapped_calculation, &mut link_hsp_params);
        if (program::blast_query_is_translated(program_number)
            || program::blast_subject_is_translated(program_number))
            && program_number != program::TBLASTX
        {
            let max_protein_gap = (options.longest_intron - 2) / crate::util::CODON_LENGTH as i32;
            if let Some(link) = link_hsp_params.as_mut() {
                if gapped_calculation {
                    if options.longest_intron == 0 {
                        link.longest_intron = (stat::DEFAULT_LONGEST_INTRON as i32 - 2)
                            / crate::util::CODON_LENGTH as i32;
                    } else if max_protein_gap <= 0 {
                        link_hsp_params = None;
                        params_options.do_sum_stats = false;
                    } else {
                        link.longest_intron = max_protein_gap;
                    }
                } else {
                    link.longest_intron = max_protein_gap.max(0);
                }
            }
        }
    }

    let mut params = HitSavingParameters {
        options: params_options,
        cutoff_score_min: 0,
        low_score: if options.low_score_perc > 0.00001 {
            vec![0; query_info.num_queries.max(0) as usize]
        } else {
            Vec::new()
        },
        cutoffs: vec![BlastGappedCutoffs::default(); query_info.num_contexts()],
        link_hsp_params,
        prelim_evalue: options.expect_value,
    };

    let status = blast_hit_saving_parameters_update(
        program_number,
        sbp,
        query_info,
        avg_subject_length,
        composition_based_stats,
        &mut params,
    );
    *parameters = Some(params);
    status
}

/// blast-rs: Nullable pointer-shaped adapter for
/// [`blast_hit_saving_parameters_new`]; not a direct NCBI C port.
pub fn blast_hit_saving_parameters_new_c(
    program_number: ProgramType,
    options: Option<&HitSavingOptions>,
    sbp: Option<&BlastScoreBlk>,
    query_info: Option<&QueryInfo>,
    avg_subject_length: i32,
    composition_based_stats: i32,
    parameters: Option<&mut Option<HitSavingParameters>>,
) -> i16 {
    let Some(parameters) = parameters else {
        return 0;
    };
    let (Some(options), Some(sbp), Some(query_info)) = (options, sbp, query_info) else {
        *parameters = None;
        return -1;
    };
    blast_hit_saving_parameters_new(
        program_number,
        options,
        sbp,
        query_info,
        avg_subject_length,
        composition_based_stats,
        parameters,
    )
}

/// NCBI: BlastInitialWordParametersUpdate (blast_parameters.c:280).
pub fn blast_initial_word_parameters_update(
    program_number: ProgramType,
    hit_params: &HitSavingParameters,
    sbp: &BlastScoreBlk,
    query_info: &QueryInfo,
    subject_length: u32,
    parameters: &mut InitialWordParameters,
) -> i16 {
    let sbp_view = BlastParameterScoreBlock::from(sbp);
    let Some((kbp_array, gapped_calculation)) = sbp_view.update_kbp_array() else {
        return -1;
    };
    if parameters.cutoffs.len() < query_info.num_contexts() {
        parameters
            .cutoffs
            .resize(query_info.num_contexts(), BlastUngappedCutoffs::default());
    }

    let gap_decay_rate = hit_params
        .link_hsp_params
        .as_ref()
        .map_or(0.0, |link| link.gap_decay_rate);
    let mut cutoff_min = i32::MAX;
    let mut xdrop_max = 0;
    for (context, context_info) in query_info.contexts.iter().enumerate() {
        let curr_cutoffs = &mut parameters.cutoffs[context];
        if !context_info.is_valid {
            curr_cutoffs.cutoff_score = i32::MAX;
            continue;
        }

        let mut gap_trigger = i32::MAX;
        if let Some(kbp_std) = sbp_view.kbp_std {
            if let Some(kbp) = kbp_std
                .get(context)
                .filter(|kbp| stat::s_blast_karlin_blk_is_valid(Some(kbp)))
            {
                gap_trigger = ((parameters.options.gap_trigger * crate::math::NCBIMATH_LN2
                    + kbp.log_k)
                    / kbp.lambda) as i32;
            }
        }

        let mut new_cutoff = if !gapped_calculation || sbp_view.matrix_only_scoring {
            let mut query_length = context_info.query_length;
            if program_number == program::BLASTN || program_number == program::MAPPING {
                query_length *= 2;
            }
            let Some(kbp) = kbp_array
                .get(context)
                .filter(|kbp| stat::s_blast_karlin_blk_is_valid(Some(kbp)))
            else {
                return -1;
            };
            // NCBI `BlastInitialWordParametersUpdate` uses
            // `subject_length * MIN(subject_length, query_length)` for
            // ungapped cutoff search space.
            let searchsp = (subject_length as u64)
                .saturating_mul((subject_length as u64).min(query_length.max(1) as u64))
                as f64;
            let (cutoff, _) = stat::blast_cutoffs(
                1,
                s_get_cutoff_evalue(program_number),
                kbp,
                searchsp,
                true,
                gap_decay_rate,
            );
            if program_number != program::BLASTN {
                cutoff.min(gap_trigger)
            } else {
                cutoff
            }
        } else {
            gap_trigger
        };
        new_cutoff = scaled_i32(new_cutoff, sbp_view.scale_factor);
        if let Some(hit_cutoff) = hit_params.cutoffs.get(context) {
            new_cutoff = new_cutoff.min(hit_cutoff.cutoff_score_max);
        }
        curr_cutoffs.cutoff_score = new_cutoff;
        cutoff_min = cutoff_min.min(new_cutoff);
        let x_dropoff = if curr_cutoffs.x_dropoff_init == 0 {
            new_cutoff
        } else {
            curr_cutoffs.x_dropoff_init
        };
        xdrop_max = xdrop_max.max(x_dropoff);
    }
    parameters.cutoff_score_min = cutoff_min;
    parameters.x_dropoff_max = xdrop_max;
    0
}

/// blast-rs: Nullable pointer-shaped adapter for
/// [`blast_initial_word_parameters_update`]; not a direct NCBI C port.
pub fn blast_initial_word_parameters_update_c(
    program_number: ProgramType,
    hit_params: Option<&HitSavingParameters>,
    sbp: Option<&BlastScoreBlk>,
    query_info: Option<&QueryInfo>,
    subject_length: u32,
    parameters: Option<&mut InitialWordParameters>,
) -> i16 {
    let (Some(hit_params), Some(sbp), Some(query_info), Some(parameters)) =
        (hit_params, sbp, query_info, parameters)
    else {
        return -1;
    };
    blast_initial_word_parameters_update(
        program_number,
        hit_params,
        sbp,
        query_info,
        subject_length,
        parameters,
    )
}

/// NCBI: BlastInitialWordParametersNew (blast_parameters.c:159).
pub fn blast_initial_word_parameters_new(
    program_number: ProgramType,
    word_options: &InitialWordOptions,
    hit_params: &HitSavingParameters,
    sbp: &BlastScoreBlk,
    query_info: &QueryInfo,
    subject_length: u32,
    parameters: &mut Option<InitialWordParameters>,
) -> i16 {
    if stat::s_blast_find_valid_karlin_blk(&sbp.kbp, query_info).is_err() {
        *parameters = None;
        return -1;
    }
    let sbp_view = BlastParameterScoreBlock::from(sbp);

    let mut params = InitialWordParameters {
        options: word_options.clone(),
        x_dropoff_max: 0,
        cutoff_score_min: 0,
        cutoffs: vec![BlastUngappedCutoffs::default(); query_info.num_contexts()],
        ungapped_extension: !program::blast_program_is_phi_blast(program_number),
        nucl_score_table: InitialWordParameters::build_nucl_score_table(
            sbp_view.reward,
            sbp_view.penalty,
        ),
    };

    for (context, context_info) in query_info.contexts.iter().enumerate() {
        if !context_info.is_valid {
            continue;
        }
        let Some(kbp) = sbp_view.valid_kbp(context) else {
            return -1;
        };
        params.cutoffs[context].x_dropoff_init =
            if program_number == program::BLASTN && sbp_view.matrix_only_scoring {
                word_options.x_dropoff as i32
            } else {
                (sbp_view.scale_factor
                    * (word_options.x_dropoff * crate::math::NCBIMATH_LN2 / kbp.lambda).ceil())
                    as i32
            };
    }

    let status = blast_initial_word_parameters_update(
        program_number,
        hit_params,
        sbp,
        query_info,
        subject_length,
        &mut params,
    );
    *parameters = Some(params);
    status
}

/// blast-rs: Nullable pointer-shaped adapter for
/// [`blast_initial_word_parameters_new`]; not a direct NCBI C port.
pub fn blast_initial_word_parameters_new_c(
    program_number: ProgramType,
    word_options: Option<&InitialWordOptions>,
    hit_params: Option<&HitSavingParameters>,
    sbp: Option<&BlastScoreBlk>,
    query_info: Option<&QueryInfo>,
    subject_length: u32,
    parameters: Option<&mut Option<InitialWordParameters>>,
) -> i16 {
    let Some(parameters) = parameters else {
        return 0;
    };
    let (Some(word_options), Some(hit_params), Some(sbp), Some(query_info)) =
        (word_options, hit_params, sbp, query_info)
    else {
        *parameters = None;
        return -1;
    };
    blast_initial_word_parameters_new(
        program_number,
        word_options,
        hit_params,
        sbp,
        query_info,
        subject_length,
        parameters,
    )
}

/// NCBI: BlastInitialWordParametersFree (blast_parameters.c:119).
pub fn blast_initial_word_parameters_free(
    parameters: &mut Option<InitialWordParameters>,
) -> Option<InitialWordParameters> {
    if let Some(parameters) = parameters.as_mut() {
        parameters.cutoffs.clear();
    }
    *parameters = None;
    None
}

/// NCBI: BlastExtensionParametersFree (blast_parameters.c:493).
pub fn blast_extension_parameters_free(
    parameters: &mut Option<ExtensionParameters>,
) -> Option<ExtensionParameters> {
    *parameters = None;
    None
}

/// NCBI: BlastScoringParametersFree (blast_parameters.c:500).
pub fn blast_scoring_parameters_free(
    parameters: &mut Option<ScoringParameters>,
) -> Option<ScoringParameters> {
    *parameters = None;
    None
}

/// NCBI: BlastHitSavingParametersFree (blast_parameters.c:718).
pub fn blast_hit_saving_parameters_free(
    parameters: &mut Option<HitSavingParameters>,
) -> Option<HitSavingParameters> {
    if let Some(parameters) = parameters.as_mut() {
        parameters.low_score.clear();
        parameters.cutoffs.clear();
        let _ = blast_link_hsp_parameters_free(&mut parameters.link_hsp_params);
    }
    *parameters = None;
    None
}

/// Effective length parameters.
#[derive(Debug, Clone)]
pub struct EffectiveLengthsParameters {
    pub options: EffectiveLengthsOptions,
    pub real_db_length: i64,
    pub real_num_seqs: i32,
}

/// blast-rs: Rust spelling of the effective-length parameters cleanup routine;
/// not a direct NCBI C port.
pub fn blast_effective_lengths_parameters_free(
    parameters: &mut Option<EffectiveLengthsParameters>,
) -> Option<EffectiveLengthsParameters> {
    *parameters = None;
    None
}

/// NCBI: s_BlastRunFullSearchCleanUp (blast_parameters.c).
pub fn s_blast_run_full_search_clean_up(
    gap_align: &mut Option<crate::blast_kappa::BlastGapAlignWorkspace>,
    score_params: &mut Option<ScoringParameters>,
    ext_params: &mut Option<ExtensionParameters>,
    hit_params: &mut Option<HitSavingParameters>,
    eff_len_params: &mut Option<EffectiveLengthsParameters>,
) {
    crate::blast_kappa::s_blast_gap_align_struct_free(gap_align);
    let _ = blast_scoring_parameters_free(score_params);
    let _ = blast_hit_saving_parameters_free(hit_params);
    let _ = blast_extension_parameters_free(ext_params);
    let _ = blast_effective_lengths_parameters_free(eff_len_params);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_scoring_params_from_options() {
        let opts = ScoringOptions::new_blastn();
        let params = ScoringParameters::from_options(&opts, 1.0);
        assert_eq!(params.reward, 1);
        assert_eq!(params.penalty, -3);
        assert_eq!(params.gap_open, 5);
        assert_eq!(params.gap_extend, 2);
        assert_eq!(params.scale_factor, 1.0);
    }

    /// Port of NCBI blastsetup_unit_test: derive parameters from blastn options.
    #[test]
    fn test_parameters_from_blastn_options() {
        let scoring_opts = ScoringOptions::new_blastn();
        let params = ScoringParameters::from_options(&scoring_opts, 1.0);

        // Parameters should faithfully carry over scoring options
        assert_eq!(params.reward, REWARD);
        assert_eq!(params.penalty, PENALTY);
        assert_eq!(params.gap_open, GAP_OPEN_NUCL);
        assert_eq!(params.gap_extend, GAP_EXTN_NUCL);
        assert_eq!(params.scale_factor, 1.0);

        // The stored options should match
        assert_eq!(params.options.reward, params.reward);
        assert_eq!(params.options.penalty, params.penalty);
        assert_eq!(params.options.gap_open, params.gap_open);
        assert_eq!(params.options.gap_extend, params.gap_extend);

        // Verify nucl score table for reward=1, penalty=-3
        let table = InitialWordParameters::build_nucl_score_table(
            scoring_opts.reward,
            scoring_opts.penalty,
        );
        // XOR 0x00 = all 4 bases match -> 4 * reward = 4
        assert_eq!(table[0x00], 4 * scoring_opts.reward);
        // XOR 0xFF = all 4 bases mismatch -> 4 * penalty = -12
        assert_eq!(table[0xFF], 4 * scoring_opts.penalty);
        // XOR 0x01 = 3 match + 1 mismatch -> 3*1 + 1*(-3) = 0
        assert_eq!(table[0x01], 3 * scoring_opts.reward + scoring_opts.penalty);
    }

    /// Port of NCBI blastsetup_unit_test: derive parameters from blastp options.
    #[test]
    fn test_parameters_from_blastp_options() {
        let scoring_opts = ScoringOptions::new_blastp();
        let params = ScoringParameters::from_options(&scoring_opts, 1.0);

        assert_eq!(params.gap_open, GAP_OPEN_PROT);
        assert_eq!(params.gap_extend, GAP_EXTN_PROT);
        assert_eq!(params.options.matrix_name.as_deref(), Some("BLOSUM62"));
        assert!(params.options.gapped_calculation);

        // With a non-default scale factor
        let scaled = ScoringParameters::from_options(&scoring_opts, 2.5);
        assert_eq!(scaled.scale_factor, 2.5);
        // Gap costs remain unscaled in the parameters struct
        assert_eq!(scaled.gap_open, GAP_OPEN_PROT);
        assert_eq!(scaled.gap_extend, GAP_EXTN_PROT);
    }

    #[test]
    fn translated_cutoff_evalue_matches_program_table() {
        assert_eq!(s_get_cutoff_evalue(program::BLASTN), stat::CUTOFF_E_BLASTN);
        assert_eq!(s_get_cutoff_evalue(program::MAPPING), stat::CUTOFF_E_BLASTN);
        assert_eq!(s_get_cutoff_evalue(program::BLASTP), stat::CUTOFF_E_BLASTP);
        assert_eq!(
            s_get_cutoff_evalue(program::RPS_BLAST),
            stat::CUTOFF_E_BLASTP
        );
        assert_eq!(s_get_cutoff_evalue(program::BLASTX), stat::CUTOFF_E_BLASTX);
        assert_eq!(
            s_get_cutoff_evalue(program::TBLASTN),
            stat::CUTOFF_E_TBLASTN
        );
        assert_eq!(
            s_get_cutoff_evalue(program::RPS_TBLASTN),
            stat::CUTOFF_E_TBLASTN
        );
        assert_eq!(
            s_get_cutoff_evalue(program::TBLASTX),
            stat::CUTOFF_E_TBLASTX
        );
    }

    #[test]
    fn translated_phi_cutoff_helpers_match_c_formula() {
        let mut query_info = QueryInfo::new_blastp(&[100]);
        query_info.contexts[0].eff_searchsp = 10_000;
        query_info.contexts[0].length_adjustment = 8;
        let pattern_info = SphiQueryInfo {
            num_patterns: 4,
            occurrences: vec![
                crate::pattern::SphiPatternInfo {
                    offset: 0,
                    length: 8,
                },
                crate::pattern::SphiPatternInfo {
                    offset: 3,
                    length: 8,
                },
                crate::pattern::SphiPatternInfo {
                    offset: 9,
                    length: 8,
                },
                crate::pattern::SphiPatternInfo {
                    offset: 20,
                    length: 8,
                },
            ],
            allocated_size: 8,
            probability: 0.025,
            pattern: Some("A-C".to_string()),
        };
        let lambda = 0.267;
        let param_c = 0.041;
        let effective = crate::hspstream::phi_blast_get_effective_number_of_patterns(
            &[0, 3, 9, 20],
            query_info.contexts[0].length_adjustment,
        );
        assert_eq!(effective, 3);

        let expect = s_get_estimated_phi_expect(
            30,
            Some(&query_info),
            Some(&pattern_info),
            lambda,
            param_c,
            effective,
        );
        let expected = 10_000.0
            * param_c
            * (1.0 + lambda * 30.0)
            * effective as f64
            * pattern_info.probability
            * (-lambda * 30.0).exp();
        assert!((expect - expected).abs() < 1e-12);

        let cutoff =
            s_phi_blast_cutoff_score(5.0, Some(&query_info), Some(&pattern_info), lambda, param_c);
        let high_expect = s_get_estimated_phi_expect(
            cutoff + 1,
            Some(&query_info),
            Some(&pattern_info),
            lambda,
            param_c,
            effective,
        );
        assert!(cutoff >= 1);
        assert!(high_expect <= 5.0);
        assert_eq!(
            s_phi_blast_cutoff_score(5.0, None, Some(&pattern_info), lambda, param_c),
            1
        );
    }

    #[test]
    fn translated_parameter_free_functions_clear_owned_parameters() {
        let scoring_opts = ScoringOptions::new_blastn();
        let extension_opts = ExtensionOptions::new_blastn();
        let hit_opts = HitSavingOptions::default();
        let initial_opts = InitialWordOptions::new_blastn();
        let effective_opts = EffectiveLengthsOptions::default();

        let mut scoring = Some(ScoringParameters::from_options(&scoring_opts, 1.0));
        let mut extension = Some(ExtensionParameters {
            options: extension_opts,
            gap_x_dropoff: 10,
            gap_x_dropoff_final: 20,
            gap_trigger: 30,
        });
        let mut hit_saving = Some(HitSavingParameters {
            options: hit_opts,
            cutoff_score_min: 12,
            low_score: vec![1, 2, 3],
            cutoffs: Vec::new(),
            link_hsp_params: None,
            prelim_evalue: 10.0,
        });
        let mut initial = Some(InitialWordParameters {
            options: initial_opts,
            x_dropoff_max: 7,
            cutoff_score_min: 4,
            cutoffs: Vec::new(),
            ungapped_extension: true,
            nucl_score_table: InitialWordParameters::build_nucl_score_table(1, -3),
        });
        let mut effective = Some(EffectiveLengthsParameters {
            options: effective_opts,
            real_db_length: 100,
            real_num_seqs: 2,
        });
        let mut link = Some(LinkHSPParameters::default());

        assert!(blast_scoring_parameters_free(&mut scoring).is_none());
        assert!(blast_extension_parameters_free(&mut extension).is_none());
        assert!(blast_hit_saving_parameters_free(&mut hit_saving).is_none());
        assert!(blast_initial_word_parameters_free(&mut initial).is_none());
        assert!(blast_effective_lengths_parameters_free(&mut effective).is_none());
        assert!(blast_link_hsp_parameters_free(&mut link).is_none());

        assert!(scoring.is_none());
        assert!(extension.is_none());
        assert!(hit_saving.is_none());
        assert!(initial.is_none());
        assert!(effective.is_none());
        assert!(link.is_none());
    }

    #[test]
    fn fill_return_cutoffs_info_copies_parameter_cutoffs() {
        let scoring_opts = ScoringOptions::new_blastn();
        let scoring = ScoringParameters::from_options(&scoring_opts, 1.0);
        let extension = ExtensionParameters {
            options: ExtensionOptions::new_blastn(),
            gap_x_dropoff: 17,
            gap_x_dropoff_final: 41,
            gap_trigger: 9,
        };
        let hit = HitSavingParameters {
            options: HitSavingOptions::default(),
            cutoff_score_min: 23,
            low_score: Vec::new(),
            cutoffs: Vec::new(),
            link_hsp_params: None,
            prelim_evalue: 10.0,
        };
        let word = InitialWordParameters {
            options: InitialWordOptions::new_blastn(),
            x_dropoff_max: 7,
            cutoff_score_min: 11,
            cutoffs: Vec::new(),
            ungapped_extension: true,
            nucl_score_table: InitialWordParameters::build_nucl_score_table(1, -3),
        };
        let mut raw = crate::diagnostics::RawCutoffs::default();

        assert_eq!(
            s_fill_return_cutoffs_info(
                Some(&mut raw),
                Some(&scoring),
                Some(&word),
                Some(&extension),
                Some(&hit),
            ),
            0
        );
        assert_eq!(raw.x_drop_ungapped, 7);
        assert_eq!(raw.ungapped_cutoff, 11);
        assert_eq!(raw.x_drop_gap, 17);
        assert_eq!(raw.x_drop_gap_final, 41);
        assert_eq!(raw.cutoff_score, 23);
        assert_eq!(
            s_fill_return_cutoffs_info(None, None, None, None, None),
            crate::util::BLASTERR_INVALIDPARAM
        );
    }

    #[test]
    fn blast_run_full_search_cleanup_clears_owned_slots() {
        let scoring_opts = ScoringOptions::new_blastn();
        let mut gap_align = crate::blast_kappa::blast_gap_align_struct_new(30);
        let mut scoring = Some(ScoringParameters::from_options(&scoring_opts, 1.0));
        let mut extension = Some(ExtensionParameters {
            options: ExtensionOptions::new_blastn(),
            gap_x_dropoff: 10,
            gap_x_dropoff_final: 20,
            gap_trigger: 30,
        });
        let mut hit = Some(HitSavingParameters {
            options: HitSavingOptions::default(),
            cutoff_score_min: 12,
            low_score: Vec::new(),
            cutoffs: Vec::new(),
            link_hsp_params: None,
            prelim_evalue: 10.0,
        });
        let mut eff = Some(EffectiveLengthsParameters {
            options: EffectiveLengthsOptions::default(),
            real_db_length: 100,
            real_num_seqs: 2,
        });

        s_blast_run_full_search_clean_up(
            &mut gap_align,
            &mut scoring,
            &mut extension,
            &mut hit,
            &mut eff,
        );
        assert!(gap_align.is_none());
        assert!(scoring.is_none());
        assert!(extension.is_none());
        assert!(hit.is_none());
        assert!(eff.is_none());
    }

    fn test_kbp() -> KarlinBlk {
        KarlinBlk {
            lambda: 0.267,
            k: 0.041,
            h: 0.14,
            log_k: 0.041_f64.ln(),
            round_down: false,
        }
    }

    fn test_score_blk(kbp: KarlinBlk, scale_factor: f64) -> BlastScoreBlk {
        let mut sbp =
            stat::blast_score_blk_new(crate::encoding::BLASTAA_SEQ_CODE, 1).expect("score block");
        sbp.kbp[0] = kbp.clone();
        sbp.kbp_std[0] = kbp.clone();
        sbp.kbp_gap = vec![kbp];
        sbp.scale_factor = scale_factor;
        sbp
    }

    #[test]
    fn blast_parameter_constructors_use_score_block_inputs() {
        let query_info = QueryInfo::new_blastp(&[100]);
        let kbp = vec![test_kbp()];
        let score_block = test_score_blk(kbp[0].clone(), 2.0);

        let scoring_opts = ScoringOptions::new_blastp();
        let mut scoring = None;
        assert_eq!(
            blast_scoring_parameters_new(Some(&scoring_opts), Some(&score_block), &mut scoring),
            0
        );
        let scoring = scoring.expect("scoring parameters");
        assert_eq!(scoring.gap_open, GAP_OPEN_PROT * 2);
        assert_eq!(scoring.gap_extend, GAP_EXTN_PROT * 2);
        assert_eq!(scoring.scale_factor, 2.0);
        assert_eq!(
            blast_scoring_parameters_new_c(Some(&scoring_opts), Some(&score_block), None),
            1
        );
        let mut scoring_c = None;
        assert_eq!(
            blast_scoring_parameters_new_c(
                Some(&scoring_opts),
                Some(&score_block),
                Some(&mut scoring_c),
            ),
            0
        );
        assert_eq!(
            scoring_c.expect("c scoring parameters").gap_open,
            GAP_OPEN_PROT * 2
        );

        let effective_opts = EffectiveLengthsOptions::default();
        let mut effective = None;
        assert_eq!(
            blast_effective_lengths_parameters_new_c(
                Some(&effective_opts),
                10_000,
                42,
                Some(&mut effective),
            ),
            0
        );
        let effective = effective.expect("effective lengths parameters");
        assert_eq!(effective.real_db_length, 10_000);
        assert_eq!(effective.real_num_seqs, 42);
        assert_eq!(
            blast_effective_lengths_parameters_new_c(Some(&effective_opts), 1, 1, None),
            1
        );

        let extension_opts = ExtensionOptions::new_blastp();
        let mut extension = None;
        assert_eq!(
            blast_extension_parameters_new(
                program::BLASTP,
                &extension_opts,
                &score_block,
                &query_info,
                &mut extension,
            ),
            0
        );
        let extension = extension.expect("extension parameters");
        let expected_xdrop =
            (extension_opts.gap_x_dropoff * crate::math::NCBIMATH_LN2 / kbp[0].lambda) as i32 * 2;
        assert_eq!(extension.gap_x_dropoff, expected_xdrop);
        assert!(extension.gap_x_dropoff_final >= extension.gap_x_dropoff);
        let mut extension_c = None;
        assert_eq!(
            blast_extension_parameters_new_c(
                program::BLASTP,
                Some(&extension_opts),
                Some(&score_block),
                Some(&query_info),
                Some(&mut extension_c),
            ),
            0
        );
        assert_eq!(
            extension_c.expect("c extension parameters").gap_x_dropoff,
            expected_xdrop
        );
        assert_eq!(
            blast_extension_parameters_new_c(
                program::BLASTP,
                Some(&extension_opts),
                Some(&score_block),
                Some(&query_info),
                None,
            ),
            0
        );

        let mut link = None;
        assert_eq!(
            blast_link_hsp_parameters_new_c(program::BLASTP, true, Some(&mut link)),
            0
        );
        assert_eq!(
            link.expect("link hsp parameters").gap_decay_rate,
            stat::BLAST_GAP_DECAY_RATE_GAPPED
        );
        assert_eq!(
            blast_link_hsp_parameters_new_c(program::BLASTP, true, None),
            -1
        );
    }

    #[test]
    fn blast_hit_and_word_parameter_updates_use_karlin_cutoffs() {
        let mut query_info = QueryInfo::new_blastp(&[80]);
        query_info.contexts[0].eff_searchsp = 80_000;
        let kbp = vec![test_kbp()];
        let score_block = test_score_blk(kbp[0].clone(), 1.0);

        let mut hit_opts = HitSavingOptions::default();
        hit_opts.do_sum_stats = true;
        let mut hit_params = None;
        assert_eq!(
            blast_hit_saving_parameters_new(
                program::BLASTP,
                &hit_opts,
                &score_block,
                &query_info,
                1000,
                0,
                &mut hit_params,
            ),
            0
        );
        let mut hit_params = hit_params.expect("hit saving parameters");
        assert_eq!(hit_params.cutoffs.len(), 1);
        assert!(hit_params.cutoffs[0].cutoff_score > 0);
        assert_eq!(
            hit_params.cutoff_score_min,
            hit_params.cutoffs[0].cutoff_score
        );
        assert!(hit_params.link_hsp_params.is_some());

        let mut no_sum_opts = hit_opts.clone();
        no_sum_opts.do_sum_stats = false;
        let mut no_sum_hit_params = None;
        assert_eq!(
            blast_hit_saving_parameters_new(
                program::BLASTP,
                &no_sum_opts,
                &score_block,
                &query_info,
                1000,
                0,
                &mut no_sum_hit_params,
            ),
            0
        );
        let no_sum_hit_params = no_sum_hit_params.expect("non-sum hit saving parameters");
        assert_eq!(
            hit_params.cutoffs[0].cutoff_score_max,
            no_sum_hit_params.cutoffs[0].cutoff_score_max
        );
        assert!(hit_params.cutoffs[0].cutoff_score <= no_sum_hit_params.cutoffs[0].cutoff_score);

        let mut too_short_gap_opts = hit_opts.clone();
        too_short_gap_opts.do_sum_stats = true;
        too_short_gap_opts.longest_intron = 2;
        let mut too_short_gap_params = None;
        assert_eq!(
            blast_hit_saving_parameters_new(
                program::BLASTX,
                &too_short_gap_opts,
                &score_block,
                &query_info,
                1000,
                0,
                &mut too_short_gap_params,
            ),
            0
        );
        let too_short_gap_params =
            too_short_gap_params.expect("too-short intron hit saving parameters");
        assert!(too_short_gap_params.link_hsp_params.is_none());
        assert!(!too_short_gap_params.options.do_sum_stats);

        assert_eq!(
            blast_hit_saving_parameters_new_c(
                program::BLASTP,
                Some(&hit_opts),
                Some(&score_block),
                Some(&query_info),
                1000,
                0,
                None,
            ),
            0
        );
        let mut hit_params_c = None;
        assert_eq!(
            blast_hit_saving_parameters_new_c(
                program::BLASTP,
                Some(&hit_opts),
                Some(&score_block),
                Some(&query_info),
                1000,
                0,
                Some(&mut hit_params_c),
            ),
            0
        );
        assert_eq!(
            hit_params_c
                .as_ref()
                .expect("c hit saving parameters")
                .cutoff_score_min,
            hit_params.cutoff_score_min
        );
        assert_eq!(
            blast_hit_saving_parameters_update_c(
                program::BLASTP,
                Some(&score_block),
                Some(&query_info),
                1000,
                0,
                hit_params_c.as_mut(),
            ),
            0
        );

        let word_opts = InitialWordOptions::new_blastp();
        let mut word_params = None;
        assert_eq!(
            blast_initial_word_parameters_new(
                program::BLASTP,
                &word_opts,
                &hit_params,
                &score_block,
                &query_info,
                1000,
                &mut word_params,
            ),
            0
        );
        let word_params = word_params.expect("initial word parameters");
        let gap_trigger = ((word_opts.gap_trigger * crate::math::NCBIMATH_LN2 + kbp[0].log_k)
            / kbp[0].lambda) as i32;
        assert_eq!(
            word_params.cutoff_score_min,
            gap_trigger.min(hit_params.cutoffs[0].cutoff_score_max)
        );
        assert_eq!(
            blast_initial_word_parameters_new_c(
                program::BLASTP,
                Some(&word_opts),
                Some(&hit_params),
                Some(&score_block),
                Some(&query_info),
                1000,
                None,
            ),
            0
        );
        let mut word_params_c = None;
        assert_eq!(
            blast_initial_word_parameters_new_c(
                program::BLASTP,
                Some(&word_opts),
                Some(&hit_params),
                Some(&score_block),
                Some(&query_info),
                1000,
                Some(&mut word_params_c),
            ),
            0
        );
        assert_eq!(
            word_params_c
                .as_ref()
                .expect("c initial word parameters")
                .cutoff_score_min,
            word_params.cutoff_score_min
        );
        assert_eq!(
            blast_initial_word_parameters_update_c(
                program::BLASTP,
                Some(&hit_params),
                Some(&score_block),
                Some(&query_info),
                1000,
                word_params_c.as_mut(),
            ),
            0
        );

        let mut zero_xdrop_opts = InitialWordOptions::new_blastp();
        zero_xdrop_opts.x_dropoff = 0.0;
        let mut zero_xdrop_params = None;
        assert_eq!(
            blast_initial_word_parameters_new(
                program::BLASTP,
                &zero_xdrop_opts,
                &hit_params,
                &score_block,
                &query_info,
                1000,
                &mut zero_xdrop_params,
            ),
            0
        );
        let zero_xdrop_params = zero_xdrop_params.expect("zero x-drop initial word parameters");
        assert_eq!(zero_xdrop_params.cutoffs[0].x_dropoff_init, 0);
        assert_eq!(
            zero_xdrop_params.x_dropoff_max,
            zero_xdrop_params.cutoffs[0].cutoff_score
        );

        assert_eq!(
            blast_link_hsp_parameters_update(Some(&word_params), Some(&mut hit_params), true),
            0
        );
        assert_eq!(
            blast_link_hsp_parameters_update_c(Some(&word_params), Some(&mut hit_params), true),
            0
        );
        assert_eq!(
            hit_params
                .link_hsp_params
                .as_ref()
                .expect("link params")
                .cutoff_small_gap,
            word_params.cutoff_score_min
        );
    }

    #[test]
    fn blastn_matrix_only_explicit_cutoff_uses_c_max_half_cutoff() {
        let query_info = QueryInfo::new_blastn(&[40]);
        let kbp = test_kbp();
        let mut score_block = test_score_blk(kbp, 2.0);
        score_block.matrix_only_scoring = true;

        let hit_opts = HitSavingOptions {
            cutoff_score: 18,
            ..HitSavingOptions::default()
        };
        let mut hit_params = None;
        assert_eq!(
            blast_hit_saving_parameters_new(
                program::BLASTN,
                &hit_opts,
                &score_block,
                &query_info,
                1000,
                0,
                &mut hit_params,
            ),
            0
        );
        let hit_params = hit_params.expect("blastn matrix-only hit parameters");
        assert_eq!(hit_params.cutoffs[0].cutoff_score, 18);
        assert_eq!(hit_params.cutoffs[0].cutoff_score_max, 9);
        assert_eq!(hit_params.cutoff_score_min, 36);
    }

    #[test]
    fn translated_calculate_link_hsp_cutoffs_matches_c_formula() {
        let mut query_info = QueryInfo::new_blastp(&[100]);
        query_info.contexts[0].query_offset = 0;
        let kbp = test_kbp();
        let score_block = test_score_blk(kbp, 1.0);
        let mut link_params = LinkHSPParameters {
            gap_prob: 1.0,
            gap_size: 40,
            overlap_size: 9,
            gap_decay_rate: stat::BLAST_GAP_DECAY_RATE,
            cutoff_small_gap: 0,
            cutoff_big_gap: 0,
            longest_intron: 0,
        };
        let word_params = InitialWordParameters {
            options: InitialWordOptions::new_blastp(),
            x_dropoff_max: 0,
            cutoff_score_min: 12,
            cutoffs: vec![BlastUngappedCutoffs {
                x_dropoff_init: 4,
                cutoff_score: 12,
            }],
            ungapped_extension: true,
            nucl_score_table: InitialWordParameters::build_nucl_score_table(1, -3),
        };

        calculate_link_hsp_cutoffs(
            program::BLASTP,
            &query_info,
            &score_block,
            Some(&mut link_params),
            &word_params,
            20_000,
            1_000,
        );

        assert_eq!(link_params.gap_prob, stat::BLAST_GAP_PROB);
        assert!(link_params.cutoff_big_gap > 0);
        assert!(link_params.cutoff_small_gap >= word_params.cutoff_score_min);
    }

    #[test]
    fn translated_parameter_debug_printers_emit_c_labels() {
        let scoring = ScoringParameters::from_options(&ScoringOptions::new_blastp(), 1.0);
        assert!(print_blast_scoring_parameters(Some(&scoring)).contains("BlastScoringParameters:"));
        assert_eq!(print_blast_scoring_parameters(None), "parameters{ null }\n");

        let query_info = QueryInfo::new_blastp(&[20]);
        let extension = ExtensionParameters {
            options: ExtensionOptions::new_blastp(),
            gap_x_dropoff: 3,
            gap_x_dropoff_final: 9,
            gap_trigger: 0,
        };
        let hit = HitSavingParameters {
            options: HitSavingOptions::default(),
            cutoff_score_min: 11,
            low_score: Vec::new(),
            cutoffs: vec![BlastGappedCutoffs {
                cutoff_score: 11,
                cutoff_score_max: 17,
            }],
            link_hsp_params: None,
            prelim_evalue: 10.0,
        };
        let word = InitialWordParameters {
            options: InitialWordOptions::new_blastp(),
            x_dropoff_max: 5,
            cutoff_score_min: 7,
            cutoffs: vec![BlastUngappedCutoffs {
                x_dropoff_init: 5,
                cutoff_score: 7,
            }],
            ungapped_extension: true,
            nucl_score_table: InitialWordParameters::build_nucl_score_table(1, -3),
        };
        let out = print_all_parameters(&hit, &extension, &word, &query_info);
        assert!(out.contains("BlastInitialWordParamters:"));
        assert!(out.contains("BlastExtensionParameters:"));
        assert!(out.contains("BlastHitSavingParameters:"));
    }
}
