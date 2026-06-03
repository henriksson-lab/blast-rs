//! Rust equivalent of `blast_setup.c` — search-space and effective-length
//! calculation, matrix-info initialisation, etc. Currently ports only
//! `BLAST_CalcEffLengths` and its private helpers; the rest of the file
//! will be filled in incrementally.

use crate::filter::{
    blast_complement_mask_locations, blast_filtering_options_from_string, blast_mask_loc_free,
    blast_mask_loc_protein_to_dna, blast_setup_get_filtering_locations, blast_setup_mask_query,
    sblast_filter_options_mask_at_hash, BlastMaskLoc,
};
use crate::options::{EffectiveLengthsOptions, QuerySetUpOptions, ScoringOptions};
use crate::options::{ExtensionOptions, HitSavingOptions, PrelimGapExt};
use crate::parameters::{
    blast_hit_saving_parameters_update, blast_initial_word_parameters_update,
    blast_link_hsp_parameters_update, EffectiveLengthsParameters, HitSavingParameters,
    InitialWordParameters,
};
use crate::pattern::{
    phi_get_pattern_occurrences, sphi_query_info_new, PhiPatternSearchBlk, SphiQueryInfo,
};
use crate::program::{
    blast_program_is_mapping, blast_program_is_phi_blast, blast_query_is_pssm,
    blast_query_is_translated, blast_subject_is_translated, ProgramType, BLASTN,
};
use crate::queryinfo::QueryInfo;
use crate::seqsrc::{
    blast_seq_src_get_max_seq_len, blast_seq_src_get_min_seq_len, blast_seq_src_get_num_seqs,
    blast_seq_src_get_num_seqs_stats, blast_seq_src_get_seq_len, blast_seq_src_get_tot_len,
    blast_seq_src_get_tot_len_stats, BlastSeqSource,
};
use crate::stat::{
    blast_compute_length_adjustment, blast_get_nucl_alpha_beta, blast_gumbel_blk_calc,
    blast_karlin_blk_gapped_calc, blast_karlin_blk_nucl_gapped_calc,
    blast_score_blk_kbp_ideal_calc, blast_score_blk_kbp_ungapped_calc, blast_score_blk_matrix_fill,
    blast_score_blk_new, blast_score_freq_new, blast_score_set_ambig_res, BlastScoreBlk, KarlinBlk,
    BLAST_GAP_EXTN_MEGABLAST, BLAST_GAP_OPEN_MEGABLAST,
};
use crate::util::{BlastSequenceBlk, SSeqRange};

pub use crate::extend::BlastGetGappedScoreType as BlastAuxGappedPath;
pub use crate::extend::BlastWordFinderType as BlastAuxWordFinder;
pub use crate::extend::JumperGappedType;

#[derive(Debug)]
pub struct BlastSetupAuxStructures {
    pub core: crate::extend::BlastCoreAuxStruct,
    pub na_extend_choice: Option<crate::lookup::NaExtendChoice>,
}

/// `BLAST_REWARD` (`blast_options.h:152`).
pub const BLAST_REWARD: i32 = 1;
/// `BLAST_PENALTY` (`blast_options.h:151`).
pub const BLAST_PENALTY: i32 = -3;

/// NCBI: BLAST_GapAlignSetUp (`blast_setup.c:888`).
#[allow(clippy::too_many_arguments)]
pub fn blast_gap_align_set_up(
    program_number: ProgramType,
    seq_src: Option<&dyn BlastSeqSource>,
    scoring_options: &ScoringOptions,
    eff_len_options: &EffectiveLengthsOptions,
    ext_options: &ExtensionOptions,
    hit_options: &HitSavingOptions,
    query_info: &mut QueryInfo,
    sbp: &mut BlastScoreBlk,
    score_params: &mut Option<crate::parameters::ScoringParameters>,
    ext_params: &mut Option<crate::parameters::ExtensionParameters>,
    hit_params: &mut Option<HitSavingParameters>,
    eff_len_params: &mut Option<EffectiveLengthsParameters>,
    gap_align: &mut Option<crate::blast_kappa::BlastGapAlignWorkspace>,
) -> i16 {
    let mut total_length = -1_i64;
    let mut num_seqs = -1_i32;

    if let Some(seq_src) = seq_src {
        total_length = blast_seq_src_get_tot_len_stats(seq_src);
        if total_length <= 0 {
            total_length = blast_seq_src_get_tot_len(seq_src);
        }

        if let Some(gbp) = sbp.gbp.as_mut() {
            let dbl = if eff_len_options.db_length > 0 {
                eff_len_options.db_length
            } else {
                total_length
            };
            gbp.db_length = if blast_subject_is_translated(program_number) {
                dbl / crate::util::CODON_LENGTH as i64
            } else {
                dbl
            };
        }

        if total_length > 0 {
            num_seqs = blast_seq_src_get_num_seqs_stats(seq_src);
            if num_seqs <= 0 {
                num_seqs = blast_seq_src_get_num_seqs(seq_src);
            }
        } else {
            let oid = 0;
            total_length = blast_seq_src_get_seq_len(seq_src, oid) as i64;
            if total_length < 0 {
                total_length = -1;
            }
            num_seqs = 1;
        }
    }

    let status = crate::parameters::blast_effective_lengths_parameters_new(
        eff_len_options,
        total_length,
        num_seqs,
        eff_len_params,
    );
    if status != 0 {
        return status;
    }
    let Some(eff_params) = eff_len_params.as_ref() else {
        return -1;
    };
    let kbp_array = if scoring_options.gapped_calculation && !sbp.kbp_gap_std.is_empty() {
        sbp.kbp_gap_std.as_slice()
    } else {
        sbp.kbp.as_slice()
    };
    let status = blast_calc_eff_lengths(
        program_number,
        scoring_options,
        eff_params,
        kbp_array,
        sbp.kbp_std.as_slice(),
        sbp.name.as_deref().unwrap_or(""),
        query_info,
    );
    if status != 0 {
        let _ = crate::parameters::blast_effective_lengths_parameters_free(eff_len_params);
        return status as i16;
    }

    let status = crate::parameters::blast_scoring_parameters_new(
        Some(scoring_options),
        Some(sbp),
        score_params,
    );
    if status != 0 {
        let _ = crate::parameters::blast_effective_lengths_parameters_free(eff_len_params);
        let _ = crate::parameters::blast_scoring_parameters_free(score_params);
        return status;
    }

    let status = crate::parameters::blast_extension_parameters_new(
        program_number,
        ext_options,
        sbp,
        query_info,
        ext_params,
    );
    if status != 0 {
        let _ = crate::parameters::blast_effective_lengths_parameters_free(eff_len_params);
        *score_params = None;
        let _ = crate::parameters::blast_scoring_parameters_free(score_params);
        let _ = crate::parameters::blast_extension_parameters_free(ext_params);
        return status;
    }

    let min_subject_length = if sbp.gbp.is_some() {
        let Some(seq_src) = seq_src else {
            return crate::diagnostics::BLASTERR_SUBJECT_LENGTH_INVALID;
        };
        let mut len = blast_seq_src_get_min_seq_len(seq_src);
        if blast_subject_is_translated(program_number) {
            len /= crate::util::CODON_LENGTH as i32;
        }
        len
    } else if num_seqs != 0 {
        (total_length / num_seqs as i64) as i32
    } else {
        -1
    };
    if min_subject_length <= 0 {
        return crate::diagnostics::BLASTERR_SUBJECT_LENGTH_INVALID;
    }

    let composition_based_stats = ext_params
        .as_ref()
        .map_or(ext_options.composition_based_stats, |params| {
            params.options.composition_based_stats
        });
    let status = crate::parameters::blast_hit_saving_parameters_new(
        program_number,
        hit_options,
        sbp,
        query_info,
        min_subject_length,
        composition_based_stats,
        hit_params,
    );
    if status != 0 {
        return status;
    }

    let max_subject_length = seq_src.map_or(0, blast_seq_src_get_max_seq_len).max(0) as u32;
    let Some(scoring_params) = score_params.as_ref() else {
        return -1;
    };
    let Some(extension_params) = ext_params.as_ref() else {
        return -1;
    };
    *gap_align = crate::blast_kappa::blast_gap_align_struct_new((
        scoring_params,
        extension_params,
        max_subject_length,
        blast_query_is_pssm(program_number),
    ));
    if gap_align.is_none() {
        return -1;
    }
    0
}

/// NCBI: BlastSetup_Validate (blast_setup.c:535).
///
/// C asserts that per-context statistics pointers are NULL for invalid
/// contexts. Rust stores most Karlin blocks by value, so only `sfp` keeps that
/// nullable shape; the observable return remains 0 if any context is valid,
/// otherwise 1.
pub fn blast_setup_validate(query_info: &QueryInfo, score_blk: Option<&BlastScoreBlk>) -> i16 {
    let mut valid_context_found = false;

    for (index, context) in query_info.contexts.iter().enumerate() {
        if context.is_valid {
            valid_context_found = true;
        } else if let Some(score_blk) = score_blk {
            if index < score_blk.sfp.len() {
                debug_assert!(score_blk.sfp[index].is_none());
            }
        }
    }

    if valid_context_found {
        0
    } else {
        1
    }
}

/// NCBI: Blast_ScoreBlkMatrixInit (blast_setup.c:330).
pub fn blast_score_blk_matrix_init(
    program_number: ProgramType,
    scoring_options: Option<&ScoringOptions>,
    sbp: Option<&mut BlastScoreBlk>,
) -> i16 {
    let (Some(scoring_options), Some(sbp)) = (scoring_options, sbp) else {
        return 1;
    };

    sbp.matrix_only_scoring = false;

    if program_number == crate::program::BLASTN || program_number == crate::program::MAPPING {
        let _ = blast_score_set_ambig_res(Some(&mut *sbp), b'N');
        let _ = blast_score_set_ambig_res(Some(&mut *sbp), b'-');

        if scoring_options.penalty == 0 && scoring_options.reward == 0 {
            sbp.matrix_only_scoring = true;
            sbp.penalty = BLAST_PENALTY;
            sbp.reward = BLAST_REWARD;
        } else {
            sbp.penalty = scoring_options.penalty;
            sbp.reward = scoring_options.reward;
        }

        if let Some(matrix) = scoring_options
            .matrix_name
            .as_deref()
            .filter(|matrix| !matrix.is_empty())
        {
            sbp.read_in_matrix = true;
            sbp.name = Some(matrix.to_string());
        } else {
            sbp.read_in_matrix = false;
            sbp.name = Some(format!("blastn matrix:{} {}", sbp.reward, sbp.penalty));
        }
    } else {
        sbp.read_in_matrix = true;
        let _ = blast_score_set_ambig_res(Some(&mut *sbp), b'X');
        sbp.name = Some(
            scoring_options
                .matrix_name
                .as_deref()
                .unwrap_or("BLOSUM62")
                .to_ascii_uppercase(),
        );
    }

    blast_score_blk_matrix_fill(sbp)
}

/// NCBI: Blast_ScoreBlkKbpGappedCalc (blast_setup.c:41).
pub fn blast_score_blk_kbp_gapped_calc(
    sbp: Option<&mut BlastScoreBlk>,
    scoring_options: Option<&ScoringOptions>,
    program: ProgramType,
    query_info: Option<&QueryInfo>,
) -> i16 {
    let (Some(sbp), Some(scoring_options), Some(query_info)) = (sbp, scoring_options, query_info)
    else {
        return 1;
    };

    if program == crate::program::BLASTN {
        // C leaves the preliminary nucleotide Gumbel path disabled.
    } else if let Some(gbp) = sbp.gbp.as_mut() {
        let retval = blast_gumbel_blk_calc(
            Some(gbp),
            scoring_options.gap_open,
            scoring_options.gap_extend,
            sbp.name.as_deref(),
            None,
        );
        if retval != 0 {
            return retval;
        }
    }

    let context_count = query_info.contexts.len();
    if sbp.kbp_gap_std.len() < context_count {
        sbp.kbp_gap_std.resize(context_count, KarlinBlk::default());
    }
    if sbp.kbp_gap_psi.len() < context_count {
        sbp.kbp_gap_psi.resize(context_count, KarlinBlk::default());
    }

    for (index, context) in query_info.contexts.iter().enumerate() {
        if !context.is_valid {
            continue;
        }

        let retval = if program == crate::program::BLASTN {
            let (reward, penalty) = if scoring_options.reward == 0 && scoring_options.penalty == 0 {
                (BLAST_REWARD, BLAST_PENALTY)
            } else {
                (scoring_options.reward, scoring_options.penalty)
            };
            let mut error_return = None;
            blast_karlin_blk_nucl_gapped_calc(
                Some(&mut sbp.kbp_gap_std[index]),
                scoring_options.gap_open,
                scoring_options.gap_extend,
                reward,
                penalty,
                sbp.kbp_std.get(index),
                Some(&mut sbp.round_down),
                Some(&mut error_return),
            )
        } else {
            blast_karlin_blk_gapped_calc(
                Some(&mut sbp.kbp_gap_std[index]),
                scoring_options.gap_open,
                scoring_options.gap_extend,
                sbp.name.as_deref(),
                None,
            )
        };
        if retval != 0 {
            return retval;
        }

        if program != crate::program::BLASTN && program != crate::program::MAPPING {
            sbp.kbp_gap_psi[index] = sbp.kbp_gap_std[index].clone();
        }
    }

    sbp.kbp_gap = if blast_query_is_pssm(program) {
        sbp.kbp_gap_psi.clone()
    } else {
        sbp.kbp_gap_std.clone()
    };

    0
}

type PhiKbpResult = Result<(f64, f64), i16>;

/// blast-rs: PHI-BLAST matrix/gap lookup table extracted from
/// `s_PHIScoreBlkFill`; not a standalone NCBI C function.
fn phi_blast_kbp_values(matrix: &str, gap_open: i32, gap_extend: i32) -> PhiKbpResult {
    let matrix = matrix.to_ascii_uppercase();
    let values = match matrix.as_str() {
        "BLOSUM62" => match (gap_open, gap_extend) {
            (11, 1) => Some((0.270, 0.047)),
            (9, 2) => Some((0.285, 0.075)),
            (8, 2) => Some((0.265, 0.046)),
            (7, 2) => Some((0.243, 0.032)),
            (12, 1) => Some((0.281, 0.057)),
            (10, 1) => Some((0.250, 0.033)),
            _ => None,
        },
        "PAM30" => match (gap_open, gap_extend) {
            (9, 1) => Some((0.295, 0.13)),
            (7, 2) => Some((0.306, 0.15)),
            (6, 2) => Some((0.292, 0.13)),
            (5, 2) => Some((0.263, 0.077)),
            (10, 1) => Some((0.309, 0.15)),
            (8, 1) => Some((0.270, 0.070)),
            _ => None,
        },
        "PAM70" => match (gap_open, gap_extend) {
            (10, 1) => Some((0.291, 0.089)),
            (8, 2) => Some((0.303, 0.13)),
            (7, 2) => Some((0.287, 0.095)),
            (6, 2) => Some((0.269, 0.079)),
            (11, 1) => Some((0.307, 0.13)),
            (9, 1) => Some((0.269, 0.058)),
            _ => None,
        },
        "BLOSUM80" => match (gap_open, gap_extend) {
            (10, 1) => Some((0.300, 0.072)),
            (8, 2) => Some((0.308, 0.089)),
            (7, 2) => Some((0.295, 0.077)),
            (6, 2) => Some((0.271, 0.051)),
            (11, 1) => Some((0.314, 0.096)),
            (9, 1) => Some((0.277, 0.046)),
            _ => None,
        },
        "BLOSUM45" => match (gap_open, gap_extend) {
            (14, 2) => Some((0.199, 0.040)),
            (13, 3) => Some((0.209, 0.057)),
            (12, 3) => Some((0.203, 0.049)),
            (11, 3) => Some((0.193, 0.037)),
            (10, 3) => Some((0.182, 0.029)),
            (15, 2) => Some((0.206, 0.049)),
            (13, 2) => Some((0.190, 0.032)),
            (12, 2) => Some((0.177, 0.023)),
            (19, 1) => Some((0.209, 0.049)),
            (18, 1) => Some((0.202, 0.041)),
            (17, 1) => Some((0.195, 0.034)),
            (16, 1) => Some((0.183, 0.024)),
            _ => None,
        },
        _ => return Err(-2),
    };
    values.ok_or(-1)
}

/// NCBI: s_PHIScoreBlkFill (blast_setup.c:129).
pub fn s_phi_score_blk_fill(
    sbp: Option<&mut BlastScoreBlk>,
    options: Option<&ScoringOptions>,
) -> i16 {
    let (Some(sbp), Some(options)) = (sbp, options) else {
        return 1;
    };
    sbp.read_in_matrix = true;
    let status = blast_score_blk_matrix_fill(sbp);
    if status != 0 {
        return status;
    }

    if sbp.kbp_gap_std.is_empty() {
        sbp.kbp_gap_std
            .resize(sbp.number_of_contexts.max(1) as usize, KarlinBlk::default());
    }
    sbp.kbp_gap_std[0] = KarlinBlk {
        h: 1.0,
        ..KarlinBlk::default()
    };
    sbp.kbp_gap = sbp.kbp_gap_std.clone();
    if !sbp.sfp.is_empty() {
        sbp.sfp[0] = blast_score_freq_new(sbp.loscore, sbp.hiscore);
    }

    let status = blast_score_blk_kbp_ideal_calc(Some(&mut *sbp));
    if status != 0 {
        return status;
    }

    let matrix = options.matrix_name.as_deref().unwrap_or("BLOSUM62");
    let (lambda, k) = match phi_blast_kbp_values(matrix, options.gap_open, options.gap_extend) {
        Ok(values) => values,
        Err(status) => return status,
    };
    sbp.kbp_gap_std[0].lambda = lambda;
    sbp.kbp_gap_std[0].k = k;
    sbp.kbp_gap_std[0].log_k = k.ln();
    sbp.kbp_gap_std[0].h = 1.0;

    let context_count = sbp.number_of_contexts.max(0) as usize;
    if sbp.kbp_gap_std.len() < context_count {
        sbp.kbp_gap_std.resize(context_count, KarlinBlk::default());
    }
    if sbp.kbp_std.len() < context_count {
        sbp.kbp_std.resize(context_count, KarlinBlk::default());
    }
    for index in 1..context_count {
        sbp.kbp_gap_std[index] = sbp.kbp_gap_std[0].clone();
    }
    for index in 0..context_count {
        sbp.kbp_std[index] = sbp.kbp_gap_std[0].clone();
    }
    sbp.kbp = sbp.kbp_std.clone();
    sbp.kbp_gap = sbp.kbp_gap_std.clone();

    0
}

/// NCBI: s_JumperScoreBlkFill (blast_setup.c:386).
pub fn s_jumper_score_blk_fill(
    sbp: Option<&mut BlastScoreBlk>,
    query_info: Option<&QueryInfo>,
) -> i16 {
    let (Some(sbp), Some(query_info)) = (sbp, query_info) else {
        return 1;
    };
    let status = blast_score_blk_kbp_ideal_calc(Some(&mut *sbp));
    if status != 0 {
        return status;
    }
    let Some(ideal) = sbp.kbp_ideal.clone() else {
        return 1;
    };

    let context_count = query_info.contexts.len();
    if sbp.kbp_std.len() < context_count {
        sbp.kbp_std.resize(context_count, KarlinBlk::default());
    }
    if sbp.kbp.len() < context_count {
        sbp.kbp.resize(context_count, KarlinBlk::default());
    }
    if sbp.kbp_gap_std.len() < context_count {
        sbp.kbp_gap_std.resize(context_count, KarlinBlk::default());
    }

    for (context, info) in query_info.contexts.iter().enumerate() {
        if info.is_valid {
            sbp.sfp[context] = None;
            sbp.kbp_std[context] = ideal.clone();
        }
    }
    sbp.kbp = sbp.kbp_std.clone();

    let Some(first_valid) = query_info
        .contexts
        .iter()
        .position(|context| context.is_valid)
    else {
        return 0;
    };
    let mut round_down = sbp.round_down;
    let status = blast_karlin_blk_nucl_gapped_calc(
        Some(&mut sbp.kbp_gap_std[first_valid]),
        BLAST_GAP_OPEN_MEGABLAST,
        BLAST_GAP_EXTN_MEGABLAST,
        BLAST_REWARD,
        BLAST_PENALTY,
        Some(&sbp.kbp_std[first_valid]),
        Some(&mut round_down),
        None,
    );
    sbp.round_down = round_down;
    if status != 0 {
        return status;
    }

    let gapped = sbp.kbp_gap_std[first_valid].clone();
    for context in first_valid + 1..context_count {
        if query_info.contexts[context].is_valid {
            sbp.kbp_gap_std[context] = gapped.clone();
        }
    }
    sbp.kbp_gap = sbp.kbp_gap_std.clone();
    0
}

/// NCBI: BlastSetup_ScoreBlkInit (blast_setup.c:456).
pub fn blast_setup_score_blk_init(
    query_blk: Option<&BlastSequenceBlk>,
    query_info: Option<&mut QueryInfo>,
    scoring_options: Option<&ScoringOptions>,
    program_number: ProgramType,
    sbpp: Option<&mut Option<BlastScoreBlk>>,
    scale_factor: f64,
) -> i16 {
    let Some(sbpp) = sbpp else {
        return 1;
    };
    let (Some(query_info), Some(scoring_options)) = (query_info, scoring_options) else {
        *sbpp = None;
        return 1;
    };

    let alphabet =
        if program_number == crate::program::BLASTN || program_number == crate::program::MAPPING {
            crate::encoding::BLASTNA_SEQ_CODE
        } else {
            crate::encoding::BLASTAA_SEQ_CODE
        };
    let Some(mut sbp) = blast_score_blk_new(alphabet, query_info.contexts.len() as i32) else {
        *sbpp = None;
        return 1;
    };
    if program_number == crate::program::BLASTN || program_number == crate::program::MAPPING {
        sbp.gbp = None;
    }
    sbp.scale_factor = scale_factor;

    let status = blast_score_blk_matrix_init(program_number, Some(scoring_options), Some(&mut sbp));
    if status != 0 {
        *sbpp = Some(sbp);
        return status;
    }

    let status = if blast_program_is_phi_blast(program_number) {
        s_phi_score_blk_fill(Some(&mut sbp), Some(scoring_options))
    } else if blast_program_is_mapping(program_number) {
        s_jumper_score_blk_fill(Some(&mut sbp), Some(query_info))
    } else {
        let Some(query_blk) = query_blk else {
            *sbpp = Some(sbp);
            return 1;
        };
        let status = blast_score_blk_kbp_ungapped_calc(
            program_number,
            Some(&mut sbp),
            query_blk.sequence.as_deref(),
            Some(query_info),
            None,
        );
        if status != 0 {
            status
        } else if scoring_options.gapped_calculation {
            blast_score_blk_kbp_gapped_calc(
                Some(&mut sbp),
                Some(scoring_options),
                program_number,
                Some(query_info),
            )
        } else {
            sbp.kbp_gap.clear();
            sbp.gbp = None;
            0
        }
    };

    *sbpp = Some(sbp);
    status
}

/// NCBI: s_BlastSetUpAuxStructures (blast_engine.c:903).
///
/// C stores word-finder and gapped-extension choices as function pointers on
/// `BlastCoreAuxStruct`. Rust keeps the owned scratch bundle in
/// [`crate::extend::BlastCoreAuxStruct`] and exposes the selected callbacks as
/// typed enums so call sites can dispatch explicitly.
pub fn s_blast_set_up_aux_structures(
    seq_src: Option<&dyn BlastSeqSource>,
    lookup_wrap: Option<&crate::lookup::LookupTableWrap>,
    word_params: Option<&InitialWordParameters>,
    ext_options: Option<&ExtensionOptions>,
    _hit_options: Option<&HitSavingOptions>,
    query: Option<&BlastSequenceBlk>,
    query_info: Option<&QueryInfo>,
    aux_struct_ptr: Option<&mut Option<BlastSetupAuxStructures>>,
    read_indexed_db: bool,
) -> i16 {
    let (
        Some(_seq_src),
        Some(lookup_wrap),
        Some(word_params),
        Some(ext_options),
        Some(query),
        Some(query_info),
        Some(aux_struct_ptr),
    ) = (
        seq_src,
        lookup_wrap,
        word_params,
        ext_options,
        query,
        query_info,
        aux_struct_ptr,
    )
    else {
        return -1;
    };

    *aux_struct_ptr = None;

    let table_type = lookup_wrap.table_type();
    let blastp = matches!(
        table_type,
        crate::lookup::LookupTableType::AaLookup
            | crate::lookup::LookupTableType::CompressedAaLookup
    );
    let rpsblast = table_type == crate::lookup::LookupTableType::RpsLookup;
    let phi_lookup = matches!(
        table_type,
        crate::lookup::LookupTableType::PhiLookup | crate::lookup::LookupTableType::PhiNaLookup
    );
    let smith_waterman = ext_options.prelim_gap_ext == PrelimGapExt::SmithWatermanScoreOnly;
    let jumper = ext_options.prelim_gap_ext == PrelimGapExt::JumperWithTraceback;

    let Some(ewp) = crate::extend::blast_extend_word_new(
        query.length.max(0) as u32,
        word_params,
        crate::extend::DiagTableType::DiagArray,
    ) else {
        return 1;
    };

    let mut na_extend_choice = None;
    let word_finder = if smith_waterman {
        BlastAuxWordFinder::None
    } else if phi_lookup {
        BlastAuxWordFinder::PhiBlast
    } else if blastp {
        BlastAuxWordFinder::Protein
    } else if rpsblast {
        BlastAuxWordFinder::Rps
    } else {
        if table_type != crate::lookup::LookupTableType::IndexedMbLookup {
            na_extend_choice = crate::lookup::blast_choose_na_extend(Some(lookup_wrap));
        }
        let selected = if read_indexed_db {
            BlastAuxWordFinder::IndexedMegablast
        } else {
            BlastAuxWordFinder::Nucleotide
        };
        if jumper {
            BlastAuxWordFinder::None
        } else {
            selected
        }
    };

    let offset_array_size = if phi_lookup {
        crate::pattern::PHI_MAX_HIT
    } else {
        crate::lookup::get_offset_array_size(lookup_wrap).max(0) as usize
    };

    let mapper_wordhits = if jumper && query_info.num_queries > 1000 {
        crate::lookup::mapper_word_hits_new(query, query_info)
    } else {
        None
    };

    let get_gapped_score = if phi_lookup {
        BlastAuxGappedPath::PhiBlast
    } else if smith_waterman {
        BlastAuxGappedPath::SmithWaterman
    } else {
        BlastAuxGappedPath::Standard
    };

    let jumper_gapped = if jumper && read_indexed_db {
        JumperGappedType::ShortReadIndexed
    } else if jumper {
        JumperGappedType::Jumper
    } else {
        JumperGappedType::None
    };

    *aux_struct_ptr = Some(BlastSetupAuxStructures {
        core: crate::extend::BlastCoreAuxStruct {
            WordFinder: word_finder,
            GetGappedScore: get_gapped_score,
            JumperGapped: jumper_gapped,
            ewp: Some(ewp),
            init_hitlist: crate::extend::blast_init_hit_list_new(),
            offset_pairs: vec![
                crate::lookup::OffsetPair {
                    query_offset: 0,
                    subject_offset: 0,
                };
                offset_array_size
            ],
            mapper_wordhits,
            translation_buffer: None,
            translation_table: None,
            translation_table_rc: None,
        },
        na_extend_choice,
    });

    0
}

/// NCBI: BLAST_MainSetUp (blast_setup.c:563).
///
/// Rust passes `translated_dna_lengths` explicitly for the optional
/// protein-to-DNA mask conversion because `QueryInfo` stores translated context
/// lengths, not the original nucleotide lengths.
pub fn blast_main_set_up(
    program_number: ProgramType,
    qsup_options: Option<&QuerySetUpOptions>,
    scoring_options: Option<&ScoringOptions>,
    query_blk: Option<&mut BlastSequenceBlk>,
    query_info: Option<&mut QueryInfo>,
    lookup_segments: Option<&mut Option<Box<crate::filter::BlastSeqLoc>>>,
    mask: Option<&mut Option<BlastMaskLoc>>,
    sbpp: Option<&mut Option<BlastScoreBlk>>,
    scale_factor: f64,
    translated_dna_lengths: Option<&[i32]>,
) -> i16 {
    let (Some(qsup_options), Some(scoring_options), Some(query_blk), Some(query_info), Some(sbpp)) =
        (qsup_options, scoring_options, query_blk, query_info, sbpp)
    else {
        return 1;
    };

    let mut parsed_filter_options = None;
    let filter_options = if let Some(options) = qsup_options.filtering_options.as_ref() {
        Some(options)
    } else if let Some(filter_string) = qsup_options.filter_string.as_deref() {
        let status = blast_filtering_options_from_string(
            program_number,
            Some(filter_string),
            &mut parsed_filter_options,
        );
        if status != 0 {
            return status;
        }
        parsed_filter_options.as_ref()
    } else {
        return 1;
    };

    let mut filter_maskloc = None;
    let status = blast_setup_get_filtering_locations(
        query_blk,
        query_info,
        program_number,
        filter_options,
        &mut filter_maskloc,
        None,
    );
    if status != 0 {
        return status;
    }

    let mask_at_hash = sblast_filter_options_mask_at_hash(filter_options);
    let filter_maskloc_ref = filter_maskloc.as_ref();
    if !mask_at_hash {
        if let Some(filter_maskloc_ref) = filter_maskloc_ref {
            blast_setup_mask_query(query_blk, query_info, filter_maskloc_ref, program_number);
        }
    }

    if let Some(lookup_segments) = lookup_segments {
        let _ = blast_complement_mask_locations(
            program_number,
            query_info,
            filter_maskloc.as_ref(),
            Some(lookup_segments),
        );
    }

    if let Some(mask_out) = mask {
        if blast_query_is_translated(program_number) {
            let Some(dna_lengths) = translated_dna_lengths else {
                return crate::util::BLASTERR_INVALIDPARAM;
            };
            if let Some(filter_maskloc) = filter_maskloc.as_mut() {
                let status =
                    blast_mask_loc_protein_to_dna(Some(filter_maskloc), query_info, dna_lengths);
                if status != 0 {
                    return status;
                }
            }
        }
        *mask_out = filter_maskloc.take();
    } else {
        let _ = blast_mask_loc_free(&mut filter_maskloc);
    }

    let status = blast_setup_score_blk_init(
        Some(&*query_blk),
        Some(&mut *query_info),
        Some(scoring_options),
        program_number,
        Some(&mut *sbpp),
        scale_factor,
    );
    if status != 0 {
        return status;
    }

    if let Some(score_blk) = sbpp.as_ref() {
        if blast_setup_validate(query_info, Some(score_blk)) != 0 {
            return 1;
        }
    }

    0
}

/// NCBI: s_GetEffectiveSearchSpaceForContext (blast_setup.c:676).
pub fn s_get_effective_search_space_for_context(
    eff_len_options: &EffectiveLengthsOptions,
    context_index: usize,
) -> i64 {
    if eff_len_options.num_searchspaces == 0 {
        0
    } else if eff_len_options.num_searchspaces == 1 {
        // C emits an `eBlastSevWarning` here when context_index != 0; we
        // omit the warning channel for now (no callers consume it).
        eff_len_options.searchsp_eff[0]
    } else {
        debug_assert!(context_index < eff_len_options.num_searchspaces as usize);
        eff_len_options.searchsp_eff[context_index]
    }
}

/// NCBI: BLAST_GetSubjectTotals (blast_setup.c:853).
pub fn blast_get_subject_totals(
    seqsrc: Option<&dyn BlastSeqSource>,
    total_length: Option<&mut i64>,
    num_seqs: Option<&mut i32>,
) {
    let (Some(total_length), Some(num_seqs)) = (total_length, num_seqs) else {
        return;
    };

    *total_length = -1;
    *num_seqs = -1;

    let Some(seqsrc) = seqsrc else {
        return;
    };

    *total_length = blast_seq_src_get_tot_len_stats(seqsrc);
    if *total_length <= 0 {
        *total_length = blast_seq_src_get_tot_len(seqsrc);
    }

    if *total_length > 0 {
        *num_seqs = blast_seq_src_get_num_seqs_stats(seqsrc);
        if *num_seqs <= 0 {
            *num_seqs = blast_seq_src_get_num_seqs(seqsrc);
        }
    } else {
        let oid = 0;
        *total_length = blast_seq_src_get_seq_len(seqsrc, oid) as i64;
        if *total_length < 0 {
            *total_length = -1;
            *num_seqs = -1;
            return;
        }
        *num_seqs = 1;
    }
}

/// NCBI: BlastEffectiveLengthsOptions_IsSearchSpaceSet (blast_options.c:1019).
///
/// Returns `true` iff any `searchsp_eff[0..num_searchspaces]` entry is non-zero.
pub fn blast_effective_lengths_options_is_search_space_set(
    options: &EffectiveLengthsOptions,
) -> bool {
    if options.searchsp_eff.is_empty() {
        return false;
    }
    options
        .searchsp_eff
        .iter()
        .take(options.num_searchspaces as usize)
        .any(|&s| s != 0)
}

/// NCBI: BLAST_GetAlphaBeta (blast_stat.c:3094).
///
/// Rust currently uses the supported protein-matrix subset needed by
/// `blast_calc_eff_lengths`.
fn blast_get_alpha_beta(
    matrix_name: &str,
    gap_open: i32,
    gap_extend: i32,
    gapped: bool,
    kbp_ungapped: &KarlinBlk,
) -> (f64, f64) {
    let mut alpha = kbp_ungapped.lambda / kbp_ungapped.h;
    let mut beta = 0.0;
    crate::stat::blast_get_alpha_beta(
        Some(matrix_name),
        &mut alpha,
        &mut beta,
        gapped,
        gap_open,
        gap_extend,
        Some(kbp_ungapped),
    );
    (alpha, beta)
}

/// NCBI: BLAST_CalcEffLengths (blast_setup.c:699).
///
/// Updates
/// `query_info` per-context `eff_searchsp` and `length_adjustment` in
/// place. Returns 0 on success, -1 if mandatory inputs are missing — the
/// same error code NCBI uses.
///
/// Inputs mirroring the C signature:
/// - `kbp_array`: per-context Karlin block array. The caller passes
///   `sbp->kbp_gap_std` when `scoring_options.gapped_calculation`,
///   `sbp->kbp` otherwise (matching the C `kbp_ptr` selection at
///   `blast_setup.c:768`).
/// - `kbp_std_array`: per-context **ungapped** Karlin block array
///   (`sbp->kbp_std`). Used by the nucleotide `Blast_GetNuclAlphaBeta`
///   path and protein `BLAST_GetAlphaBeta` path for alpha/beta lookup.
/// - `matrix_name`: protein matrix name. Ignored for nucleotide / mapping
///   / phi programs.
pub fn blast_calc_eff_lengths(
    program_number: ProgramType,
    scoring_options: &ScoringOptions,
    eff_len_params: &EffectiveLengthsParameters,
    kbp_array: &[KarlinBlk],
    kbp_std_array: &[KarlinBlk],
    matrix_name: &str,
    query_info: &mut QueryInfo,
) -> i32 {
    let eff_len_options = &eff_len_params.options;

    let mut db_length: i64 = if eff_len_options.db_length > 0 {
        eff_len_options.db_length
    } else {
        eff_len_params.real_db_length
    };

    if db_length == 0 && !blast_effective_lengths_options_is_search_space_set(eff_len_options) {
        return 0;
    }

    if blast_subject_is_translated(program_number) {
        db_length /= 3;
    }

    let db_num_seqs = if eff_len_options.dbseq_num > 0 {
        eff_len_options.dbseq_num
    } else {
        eff_len_params.real_num_seqs
    };

    if blast_program_is_mapping(program_number) {
        for ctx in &mut query_info.contexts {
            ctx.eff_searchsp = db_length;
        }
        return 0;
    }

    if blast_program_is_phi_blast(program_number) {
        for ctx in &mut query_info.contexts {
            let eff = db_length - (db_num_seqs as i64) * (ctx.length_adjustment as i64);
            ctx.eff_searchsp = eff;
        }
        return 0;
    }

    for (index, ctx) in query_info.contexts.iter_mut().enumerate() {
        let mut effective_search_space =
            s_get_effective_search_space_for_context(eff_len_options, index);

        let query_length = ctx.query_length;
        if !ctx.is_valid || query_length <= 0 {
            ctx.eff_searchsp = effective_search_space;
            ctx.length_adjustment = 0;
            continue;
        }

        let Some(kbp) = kbp_array.get(index) else {
            return -1;
        };
        let Some(kbp_std) = kbp_std_array.get(index) else {
            return -1;
        };

        let (alpha, beta) = if program_number == BLASTN {
            // C: when reward and penalty are zero, the matrix-only scoring
            // path substitutes BLAST_REWARD/BLAST_PENALTY before the
            // alpha/beta lookup so the KA calculation stays consistent
            // (`blast_setup.c:796`).
            let (rew, pen) = if scoring_options.reward == 0 && scoring_options.penalty == 0 {
                (BLAST_REWARD, BLAST_PENALTY)
            } else {
                (scoring_options.reward, scoring_options.penalty)
            };
            let (mut alpha, mut beta) = (0.0, 0.0);
            let _ = blast_get_nucl_alpha_beta(
                rew,
                pen,
                scoring_options.gap_open,
                scoring_options.gap_extend,
                kbp_std.lambda,
                kbp_std.h,
                scoring_options.gapped_calculation,
                &mut alpha,
                &mut beta,
            );
            (alpha, beta)
        } else {
            blast_get_alpha_beta(
                matrix_name,
                scoring_options.gap_open,
                scoring_options.gap_extend,
                scoring_options.gapped_calculation,
                kbp_std,
            )
        };

        let mut length_adjustment = 0;
        let _ = blast_compute_length_adjustment(
            kbp.k,
            kbp.log_k,
            alpha / kbp.lambda,
            beta,
            query_length,
            db_length,
            db_num_seqs,
            Some(&mut length_adjustment),
        );

        if effective_search_space == 0 {
            // C: `effective_db_length = db_length - db_num_seqs *
            // length_adjustment`, clamped to >= 1. Then
            // `effective_search_space = effective_db_length *
            // (query_length - length_adjustment)`.
            let mut effective_db_length =
                db_length - (db_num_seqs as i64) * (length_adjustment as i64);
            if effective_db_length <= 0 {
                effective_db_length = 1;
            }
            effective_search_space =
                effective_db_length * (query_length as i64 - length_adjustment as i64);
        }

        ctx.eff_searchsp = effective_search_space;
        ctx.length_adjustment = length_adjustment;
    }

    0
}

/// blast-rs: Nullable pointer-shaped adapter for effective-length calculation;
/// not a direct NCBI C port.
pub fn blast_calc_eff_lengths_c(
    program_number: ProgramType,
    scoring_options: Option<&ScoringOptions>,
    eff_len_params: Option<&EffectiveLengthsParameters>,
    kbp_array: Option<&[KarlinBlk]>,
    kbp_std_array: Option<&[KarlinBlk]>,
    matrix_name: Option<&str>,
    query_info: Option<&mut QueryInfo>,
) -> i32 {
    let (Some(scoring_options), Some(eff_len_params), Some(kbp_array), Some(kbp_std_array)) =
        (scoring_options, eff_len_params, kbp_array, kbp_std_array)
    else {
        return -1;
    };
    let Some(query_info) = query_info else {
        return -1;
    };
    blast_calc_eff_lengths(
        program_number,
        scoring_options,
        eff_len_params,
        kbp_array,
        kbp_std_array,
        matrix_name.unwrap_or(""),
        query_info,
    )
}

/// NCBI: BLAST_OneSubjectUpdateParameters (blast_setup.c:1001).
pub fn blast_one_subject_update_parameters(
    program_number: ProgramType,
    subject_length: u32,
    scoring_options: Option<&ScoringOptions>,
    query_info: Option<&mut QueryInfo>,
    sbp: Option<&crate::stat::BlastScoreBlk>,
    hit_params: Option<&mut HitSavingParameters>,
    word_params: Option<&mut InitialWordParameters>,
    eff_len_params: Option<&mut EffectiveLengthsParameters>,
) -> i16 {
    let (
        Some(scoring_options),
        Some(query_info),
        Some(sbp),
        Some(hit_params),
        Some(eff_len_params),
    ) = (scoring_options, query_info, sbp, hit_params, eff_len_params)
    else {
        return -1;
    };

    eff_len_params.real_db_length = subject_length as i64;
    let kbp_array = if scoring_options.gapped_calculation {
        &sbp.kbp_gap_std
    } else {
        &sbp.kbp
    };
    let status = blast_calc_eff_lengths(
        program_number,
        scoring_options,
        eff_len_params,
        kbp_array,
        &sbp.kbp_std,
        sbp.name.as_deref().unwrap_or(""),
        query_info,
    );
    if status != 0 {
        return status as i16;
    }

    let _ = blast_hit_saving_parameters_update(
        program_number,
        sbp,
        query_info,
        subject_length as i32,
        0,
        hit_params,
    );

    if let Some(word_params) = word_params {
        let _ = blast_initial_word_parameters_update(
            program_number,
            hit_params,
            sbp,
            query_info,
            subject_length,
            word_params,
        );
        let _ = blast_link_hsp_parameters_update(
            Some(word_params),
            Some(hit_params),
            scoring_options.gapped_calculation,
        );
    }

    status as i16
}

/// NCBI: Blast_SetPHIPatternInfo (blast_setup.c:1065).
///
/// Rust keeps PHI pattern metadata separate from [`QueryInfo`], so this
/// returns the allocated [`SphiQueryInfo`] through `pattern_info_out` while
/// preserving C's `query_info->contexts[0].length_adjustment` side effect.
pub fn blast_set_phi_pattern_info(
    program: ProgramType,
    pattern_blk: Option<&PhiPatternSearchBlk>,
    query_sequence: Option<&[u8]>,
    lookup_segments: Option<&[SSeqRange]>,
    query_info: Option<&mut QueryInfo>,
    pattern_info_out: Option<&mut Option<SphiQueryInfo>>,
) -> i16 {
    let (
        Some(pattern_blk),
        Some(query_sequence),
        Some(lookup_segments),
        Some(query_info),
        Some(pattern_info_out),
    ) = (
        pattern_blk,
        query_sequence,
        lookup_segments,
        query_info,
        pattern_info_out,
    )
    else {
        return -1;
    };

    let mut pattern_info = sphi_query_info_new().unwrap_or_default();
    let is_na = (program == crate::program::PHI_BLASTN) as u8;
    let num_patterns = phi_get_pattern_occurrences(
        pattern_blk,
        query_sequence,
        lookup_segments,
        is_na,
        &mut pattern_info,
        query_sequence.len() as i32,
    );

    if num_patterns == 0 || num_patterns == i32::MAX || num_patterns < 0 {
        *pattern_info_out = Some(pattern_info);
        return -1;
    }

    pattern_info.probability = pattern_blk.pattern_probability;
    pattern_info.pattern = pattern_blk.pattern.clone();
    if let Some(context) = query_info.contexts.get_mut(0) {
        context.length_adjustment = pattern_blk.min_pattern_match_length;
    }
    *pattern_info_out = Some(pattern_info);
    0
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::link_hsps::LinkHSPParameters;
    use crate::options::{HitSavingOptions, InitialWordOptions};
    use crate::parameters::{BlastGappedCutoffs, BlastUngappedCutoffs};
    use crate::pattern::sphi_pattern_search_blk_new;
    use crate::queryinfo::ContextInfo;
    use crate::seqsrc::{GetSeqArg, SeqData, SeqEncoding};

    struct SetupSeqSrc {
        seqs: Vec<Vec<u8>>,
        total_stats: i64,
        num_stats: i32,
    }

    impl BlastSeqSource for SetupSeqSrc {
        fn num_seqs(&self) -> i32 {
            self.seqs.len() as i32
        }

        fn num_seqs_stats(&self) -> i32 {
            self.num_stats
        }

        fn total_length(&self) -> i64 {
            self.seqs.iter().map(|seq| seq.len() as i64).sum()
        }

        fn total_length_stats(&self) -> i64 {
            self.total_stats
        }

        fn max_seq_len(&self) -> i32 {
            self.seqs
                .iter()
                .map(|seq| seq.len() as i32)
                .max()
                .unwrap_or(0)
        }

        fn avg_seq_len(&self) -> i32 {
            if self.seqs.is_empty() {
                0
            } else {
                (self.total_length() / self.seqs.len() as i64) as i32
            }
        }

        fn name(&self) -> &str {
            "setup-test"
        }

        fn is_protein(&self) -> bool {
            true
        }

        fn seq_len(&self, oid: i32) -> i32 {
            self.seqs
                .get(oid.max(0) as usize)
                .map_or(-1, |seq| seq.len() as i32)
        }

        fn get_sequence(&self, arg: &GetSeqArg) -> Option<SeqData> {
            let sequence = self.seqs.get(arg.oid.max(0) as usize)?.clone();
            Some(SeqData {
                length: sequence.len() as i32,
                sequence,
            })
        }

        fn iter_oids(&self) -> Box<dyn Iterator<Item = i32> + '_> {
            Box::new(0..self.seqs.len() as i32)
        }
    }

    #[test]
    fn local_blast_get_alpha_beta_uses_best_gapped_matrix_row() {
        let kbp = KarlinBlk {
            lambda: 0.318,
            k: 0.13,
            log_k: 0.13_f64.ln(),
            h: 0.4,
            round_down: false,
        };
        let fallback_alpha = kbp.lambda / kbp.h;

        let (alpha, beta) = blast_get_alpha_beta("BLOSUM62", 0, 0, true, &kbp);

        assert_ne!(alpha, fallback_alpha);
        assert_ne!(beta, 0.0);
    }

    #[test]
    fn gap_align_set_up_updates_gumbel_db_length_and_allocates_outputs() {
        let src = SetupSeqSrc {
            seqs: vec![vec![b'A'; 300], vec![b'C'; 150]],
            total_stats: 900,
            num_stats: 2,
        };
        let mut query_info = QueryInfo {
            num_queries: 1,
            contexts: vec![ContextInfo {
                query_offset: 0,
                query_length: 100,
                eff_searchsp: 0,
                length_adjustment: 0,
                query_index: 0,
                frame: 1,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: 100,
            min_length: 100,
        };
        let scoring = ScoringOptions {
            matrix_path: None,
            reward: 0,
            penalty: 0,
            gap_open: 11,
            gap_extend: 1,
            shift_pen: i16::MAX as i32,
            gapped_calculation: true,
            complexity_adjusted_scoring: false,
            matrix_name: Some("BLOSUM62".to_string()),
            is_ooframe: false,
            program_number: crate::program::TBLASTN,
        };
        let eff_options = EffectiveLengthsOptions {
            db_length: 1_200,
            ..EffectiveLengthsOptions::default()
        };
        let ext_options = ExtensionOptions::new_blastp();
        let hit_options = HitSavingOptions::default();
        let kbp = protein_kbp();
        let mut sbp =
            crate::stat::blast_score_blk_new(crate::encoding::BLASTAA_SEQ_CODE, 1).expect("sbp");
        sbp.name = Some("BLOSUM62".to_string());
        sbp.kbp = vec![kbp.clone()];
        sbp.kbp_std = vec![kbp.clone()];
        sbp.kbp_gap = vec![kbp.clone()];
        sbp.kbp_gap_std = vec![kbp.clone()];

        let mut score_params = None;
        let mut ext_params = None;
        let mut hit_params = None;
        let mut eff_len_params = None;
        let mut gap_align = None;

        assert_eq!(
            blast_gap_align_set_up(
                crate::program::TBLASTN,
                Some(&src),
                &scoring,
                &eff_options,
                &ext_options,
                &hit_options,
                &mut query_info,
                &mut sbp,
                &mut score_params,
                &mut ext_params,
                &mut hit_params,
                &mut eff_len_params,
                &mut gap_align,
            ),
            0
        );

        assert_eq!(
            sbp.gbp.as_ref().map(|gbp| gbp.db_length),
            Some(1_200 / crate::util::CODON_LENGTH as i64)
        );
        assert!(score_params.is_some());
        assert!(ext_params.is_some());
        assert!(hit_params.is_some());
        assert!(eff_len_params.is_some());
        let gap_align = gap_align.expect("gap align");
        let ext_params_ref = ext_params.as_ref().expect("extension parameters");
        assert_eq!(gap_align.gap_x_dropoff, ext_params_ref.gap_x_dropoff);
        assert_eq!(
            gap_align.max_mismatches,
            ext_params_ref.options.max_mismatches
        );
        assert_eq!(
            gap_align.mismatch_window,
            ext_params_ref.options.mismatch_window
        );
        assert!(gap_align.fwd_prelim_tback.is_some());
        assert!(gap_align.rev_prelim_tback.is_some());
        assert_eq!(
            gap_align.greedy_align_mem.as_ref().map(|mem| mem.max_dist),
            Some(151.min(1000))
        );
        assert!(query_info.contexts[0].eff_searchsp > 0);
    }

    #[test]
    fn setup_validate_reports_any_valid_context() {
        let mut query_info = QueryInfo {
            num_queries: 2,
            contexts: vec![
                ContextInfo {
                    query_offset: 0,
                    query_length: 10,
                    eff_searchsp: 0,
                    length_adjustment: 0,
                    query_index: 0,
                    frame: 0,
                    is_valid: false,
                    segment_flags: crate::queryinfo::E_NO_SEGMENTS,
                },
                ContextInfo {
                    query_offset: 11,
                    query_length: 12,
                    eff_searchsp: 0,
                    length_adjustment: 0,
                    query_index: 1,
                    frame: 0,
                    is_valid: true,
                    segment_flags: crate::queryinfo::E_NO_SEGMENTS,
                },
            ],
            max_length: 12,
            min_length: 0,
        };
        let mut sbp =
            crate::stat::blast_score_blk_new(crate::encoding::BLASTAA_SEQ_CODE, 2).expect("sbp");
        sbp.sfp[0] = None;

        assert_eq!(blast_setup_validate(&query_info, Some(&sbp)), 0);

        query_info.contexts[1].is_valid = false;
        assert_eq!(blast_setup_validate(&query_info, Some(&sbp)), 1);
        assert_eq!(blast_setup_validate(&query_info, None), 1);
    }

    #[test]
    fn score_blk_matrix_init_matches_nucleotide_and_protein_setup() {
        let nucleotide_options = ScoringOptions {
            matrix_path: None,
            reward: 0,
            penalty: 0,
            gap_open: 0,
            gap_extend: 0,
            shift_pen: 0,
            gapped_calculation: false,
            complexity_adjusted_scoring: false,
            matrix_name: None,
            is_ooframe: false,
            program_number: crate::program::UNDEFINED,
        };
        let mut nuc_sbp =
            crate::stat::blast_score_blk_new(crate::encoding::BLASTNA_SEQ_CODE, 2).expect("sbp");

        assert_eq!(
            blast_score_blk_matrix_init(
                crate::program::BLASTN,
                Some(&nucleotide_options),
                Some(&mut nuc_sbp),
            ),
            0
        );
        assert!(nuc_sbp.matrix_only_scoring);
        assert_eq!(nuc_sbp.reward, BLAST_REWARD);
        assert_eq!(nuc_sbp.penalty, BLAST_PENALTY);
        assert_eq!(nuc_sbp.name.as_deref(), Some("blastn matrix:1 -3"));
        assert!(!nuc_sbp.read_in_matrix);
        assert!(nuc_sbp
            .ambiguous_res
            .contains(&crate::encoding::IUPACNA_TO_BLASTNA[b'N' as usize]));

        let protein_options = ScoringOptions {
            matrix_path: None,
            reward: 0,
            penalty: 0,
            gap_open: 11,
            gap_extend: 1,
            shift_pen: i16::MAX as i32,
            gapped_calculation: true,
            complexity_adjusted_scoring: false,
            matrix_name: Some("blosum62".to_string()),
            is_ooframe: false,
            program_number: crate::program::UNDEFINED,
        };
        let mut protein_sbp =
            crate::stat::blast_score_blk_new(crate::encoding::BLASTAA_SEQ_CODE, 1).expect("sbp");

        assert_eq!(
            blast_score_blk_matrix_init(
                crate::program::BLASTP,
                Some(&protein_options),
                Some(&mut protein_sbp),
            ),
            0
        );
        assert_eq!(protein_sbp.name.as_deref(), Some("BLOSUM62"));
        assert!(protein_sbp.read_in_matrix);
        assert!(!protein_sbp.matrix_only_scoring);
        assert!(protein_sbp
            .ambiguous_res
            .contains(&crate::encoding::NCBISTDAA_X));
        assert_eq!(
            blast_score_blk_matrix_init(crate::program::BLASTP, None, Some(&mut protein_sbp)),
            1
        );
    }

    #[test]
    fn score_blk_kbp_gapped_calc_fills_valid_contexts_and_kbp_gap_vec() {
        let mut query_info = QueryInfo {
            num_queries: 2,
            contexts: vec![
                ContextInfo {
                    query_offset: 0,
                    query_length: 10,
                    eff_searchsp: 0,
                    length_adjustment: 0,
                    query_index: 0,
                    frame: 0,
                    is_valid: true,
                    segment_flags: crate::queryinfo::E_NO_SEGMENTS,
                },
                ContextInfo {
                    query_offset: 11,
                    query_length: 10,
                    eff_searchsp: 0,
                    length_adjustment: 0,
                    query_index: 1,
                    frame: 0,
                    is_valid: false,
                    segment_flags: crate::queryinfo::E_NO_SEGMENTS,
                },
            ],
            max_length: 10,
            min_length: 0,
        };
        let scoring = ScoringOptions {
            matrix_path: None,
            reward: 0,
            penalty: 0,
            gap_open: 11,
            gap_extend: 1,
            shift_pen: i16::MAX as i32,
            gapped_calculation: true,
            complexity_adjusted_scoring: false,
            matrix_name: Some("BLOSUM62".to_string()),
            is_ooframe: false,
            program_number: crate::program::UNDEFINED,
        };
        let mut protein_sbp =
            crate::stat::blast_score_blk_new(crate::encoding::BLASTAA_SEQ_CODE, 2).expect("sbp");
        protein_sbp.name = Some("BLOSUM62".to_string());

        assert_eq!(
            blast_score_blk_kbp_gapped_calc(
                Some(&mut protein_sbp),
                Some(&scoring),
                crate::program::BLASTP,
                Some(&query_info),
            ),
            0
        );
        assert!(protein_sbp.kbp_gap_std[0].lambda > 0.0);
        assert_eq!(protein_sbp.kbp_gap_std[1].lambda, 0.0);
        assert_eq!(
            protein_sbp.kbp_gap_psi[0].lambda,
            protein_sbp.kbp_gap_std[0].lambda
        );
        assert_eq!(
            protein_sbp.kbp_gap[0].lambda,
            protein_sbp.kbp_gap_std[0].lambda
        );

        query_info.contexts[1].is_valid = true;
        let mut psi_sbp =
            crate::stat::blast_score_blk_new(crate::encoding::BLASTAA_SEQ_CODE, 2).expect("sbp");
        psi_sbp.name = Some("BLOSUM62".to_string());
        assert_eq!(
            blast_score_blk_kbp_gapped_calc(
                Some(&mut psi_sbp),
                Some(&scoring),
                crate::program::PSI_BLAST,
                Some(&query_info),
            ),
            0
        );
        assert_eq!(psi_sbp.kbp_gap[0].lambda, psi_sbp.kbp_gap_psi[0].lambda);
        assert!(psi_sbp.kbp_gap[1].lambda > 0.0);

        assert_eq!(
            blast_score_blk_kbp_gapped_calc(
                None,
                Some(&scoring),
                crate::program::BLASTP,
                Some(&query_info)
            ),
            1
        );
    }

    #[test]
    fn score_blk_kbp_gapped_calc_uses_blastn_matrix_only_defaults() {
        let query_info = QueryInfo {
            num_queries: 1,
            contexts: vec![ContextInfo {
                query_offset: 0,
                query_length: 20,
                eff_searchsp: 0,
                length_adjustment: 0,
                query_index: 0,
                frame: 1,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: 20,
            min_length: 0,
        };
        let scoring = ScoringOptions {
            matrix_path: None,
            reward: 0,
            penalty: 0,
            gap_open: 5,
            gap_extend: 2,
            shift_pen: 0,
            gapped_calculation: true,
            complexity_adjusted_scoring: false,
            matrix_name: None,
            is_ooframe: false,
            program_number: crate::program::UNDEFINED,
        };
        let mut sbp =
            crate::stat::blast_score_blk_new(crate::encoding::BLASTNA_SEQ_CODE, 1).expect("sbp");
        sbp.kbp_std[0] = KarlinBlk {
            lambda: 1.37,
            k: 0.711,
            log_k: 0.711f64.ln(),
            h: 1.31,
            ..KarlinBlk::default()
        };

        assert_eq!(
            blast_score_blk_kbp_gapped_calc(
                Some(&mut sbp),
                Some(&scoring),
                crate::program::BLASTN,
                Some(&query_info),
            ),
            0
        );
        assert!(sbp.kbp_gap_std[0].lambda > 0.0);
        assert_eq!(sbp.kbp_gap[0].lambda, sbp.kbp_gap_std[0].lambda);
    }

    #[test]
    fn score_blk_kbp_gapped_calc_rejects_unsupported_blastn_gap_costs() {
        let query_info = QueryInfo {
            num_queries: 1,
            contexts: vec![ContextInfo {
                query_offset: 0,
                query_length: 20,
                eff_searchsp: 0,
                length_adjustment: 0,
                query_index: 0,
                frame: 1,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: 20,
            min_length: 0,
        };
        let scoring = ScoringOptions {
            matrix_path: None,
            reward: 1,
            penalty: -2,
            gap_open: 1,
            gap_extend: 3,
            shift_pen: 0,
            gapped_calculation: true,
            complexity_adjusted_scoring: false,
            matrix_name: None,
            is_ooframe: false,
            program_number: crate::program::UNDEFINED,
        };
        let mut sbp =
            crate::stat::blast_score_blk_new(crate::encoding::BLASTNA_SEQ_CODE, 1).expect("sbp");
        sbp.kbp_std[0] = KarlinBlk {
            lambda: 1.28,
            k: 0.46,
            log_k: 0.46f64.ln(),
            h: 0.85,
            ..KarlinBlk::default()
        };

        assert_eq!(
            blast_score_blk_kbp_gapped_calc(
                Some(&mut sbp),
                Some(&scoring),
                crate::program::BLASTN,
                Some(&query_info),
            ),
            1
        );
        assert_eq!(sbp.kbp_gap_std[0].lambda, 0.0);
        assert_eq!(sbp.kbp_gap_std[0].k, 0.0);
        assert_eq!(sbp.kbp_gap_std[0].h, 0.0);
    }

    #[test]
    fn phi_score_blk_fill_uses_phi_specific_gap_table() {
        let options = ScoringOptions {
            matrix_path: None,
            reward: 0,
            penalty: 0,
            gap_open: 11,
            gap_extend: 1,
            shift_pen: i16::MAX as i32,
            gapped_calculation: true,
            complexity_adjusted_scoring: false,
            matrix_name: Some("BLOSUM62".to_string()),
            is_ooframe: false,
            program_number: crate::program::UNDEFINED,
        };
        let mut sbp =
            crate::stat::blast_score_blk_new(crate::encoding::BLASTAA_SEQ_CODE, 2).expect("sbp");
        sbp.name = Some("BLOSUM62".to_string());

        assert_eq!(s_phi_score_blk_fill(Some(&mut sbp), Some(&options)), 0);
        assert!((sbp.kbp_gap_std[0].lambda - 0.270).abs() < 1e-12);
        assert!((sbp.kbp_gap_std[0].k - 0.047).abs() < 1e-12);
        assert_eq!(sbp.kbp_gap_std[0].h, 1.0);
        assert!(sbp.sfp[0].is_some());
        assert_eq!(sbp.kbp_gap[1].lambda, sbp.kbp_gap_std[0].lambda);
        assert_eq!(sbp.kbp[1].lambda, sbp.kbp_gap_std[0].lambda);

        let unsupported_gap = ScoringOptions {
            gap_open: 99,
            ..options.clone()
        };
        assert_eq!(
            s_phi_score_blk_fill(Some(&mut sbp), Some(&unsupported_gap)),
            -1
        );

        let unsupported_matrix = ScoringOptions {
            matrix_name: Some("BLOSUM50".to_string()),
            ..options
        };
        sbp.name = Some("BLOSUM50".to_string());
        assert_eq!(
            s_phi_score_blk_fill(Some(&mut sbp), Some(&unsupported_matrix)),
            -2
        );
    }

    #[test]
    fn jumper_score_blk_fill_creates_mapping_fake_blocks() {
        let query_info = QueryInfo {
            num_queries: 2,
            contexts: vec![
                ContextInfo {
                    query_offset: 0,
                    query_length: 20,
                    eff_searchsp: 0,
                    length_adjustment: 0,
                    query_index: 0,
                    frame: 1,
                    is_valid: true,
                    segment_flags: crate::queryinfo::E_NO_SEGMENTS,
                },
                ContextInfo {
                    query_offset: 21,
                    query_length: 20,
                    eff_searchsp: 0,
                    length_adjustment: 0,
                    query_index: 1,
                    frame: -1,
                    is_valid: false,
                    segment_flags: crate::queryinfo::E_NO_SEGMENTS,
                },
            ],
            max_length: 20,
            min_length: 0,
        };
        let mut sbp =
            crate::stat::blast_score_blk_new(crate::encoding::BLASTNA_SEQ_CODE, 2).expect("sbp");
        sbp.reward = BLAST_REWARD;
        sbp.penalty = BLAST_PENALTY;

        assert_eq!(
            s_jumper_score_blk_fill(Some(&mut sbp), Some(&query_info)),
            0
        );
        assert!(sbp.kbp_std[0].lambda > 0.0);
        assert_eq!(sbp.kbp_std[1].lambda, 0.0);
        assert!(sbp.kbp_gap_std[0].lambda > 0.0);
        assert_eq!(sbp.kbp_gap[0].lambda, sbp.kbp_gap_std[0].lambda);
        assert!(sbp.sfp[0].is_none());
    }

    #[test]
    fn setup_score_blk_init_dispatches_program_specific_paths() {
        let protein_query = crate::encoding::encode_ncbistdaa_sequence(b"ARNDCEQGHILKMFPSTWYV");
        let query_blk = BlastSequenceBlk {
            sequence: Some(protein_query.clone()),
            sequence_start: Some(protein_query),
            ..BlastSequenceBlk::default()
        };
        let mut query_info = QueryInfo::new_blastp(&[20]);
        let scoring = ScoringOptions {
            matrix_path: None,
            reward: 0,
            penalty: 0,
            gap_open: 11,
            gap_extend: 1,
            shift_pen: i16::MAX as i32,
            gapped_calculation: true,
            complexity_adjusted_scoring: false,
            matrix_name: Some("BLOSUM62".to_string()),
            is_ooframe: false,
            program_number: crate::program::UNDEFINED,
        };
        let mut sbp = None;
        assert_eq!(
            blast_setup_score_blk_init(
                Some(&query_blk),
                Some(&mut query_info),
                Some(&scoring),
                crate::program::BLASTP,
                Some(&mut sbp),
                1.25,
            ),
            0
        );
        let sbp_ref = sbp.as_ref().expect("score block");
        assert_eq!(sbp_ref.alphabet_code, crate::encoding::BLASTAA_SEQ_CODE);
        assert_eq!(sbp_ref.scale_factor, 1.25);
        assert!(sbp_ref.kbp[0].lambda > 0.0);
        assert!(sbp_ref.kbp_gap[0].lambda > 0.0);

        let mut mapping_query_info = QueryInfo::new_blastn(&[20]);
        let mapping_scoring = ScoringOptions {
            matrix_path: None,
            reward: BLAST_REWARD,
            penalty: BLAST_PENALTY,
            gap_open: BLAST_GAP_OPEN_MEGABLAST,
            gap_extend: BLAST_GAP_EXTN_MEGABLAST,
            shift_pen: 0,
            gapped_calculation: true,
            complexity_adjusted_scoring: false,
            matrix_name: None,
            is_ooframe: false,
            program_number: crate::program::UNDEFINED,
        };
        let mut mapping_sbp = None;
        assert_eq!(
            blast_setup_score_blk_init(
                None,
                Some(&mut mapping_query_info),
                Some(&mapping_scoring),
                crate::program::MAPPING,
                Some(&mut mapping_sbp),
                1.0,
            ),
            0
        );
        let mapping_ref = mapping_sbp.as_ref().expect("mapping score block");
        assert_eq!(mapping_ref.alphabet_code, crate::encoding::BLASTNA_SEQ_CODE);
        assert!(mapping_ref.kbp_gap[0].lambda > 0.0);

        let mut missing_out = None;
        assert_eq!(
            blast_setup_score_blk_init(
                Some(&query_blk),
                Some(&mut query_info),
                None,
                crate::program::BLASTP,
                Some(&mut missing_out),
                1.0,
            ),
            1
        );
        assert!(missing_out.is_none());
    }

    #[test]
    fn main_setup_runs_filter_lookup_mask_and_score_init() {
        let protein_query = crate::encoding::encode_ncbistdaa_sequence(b"ARNDCEQGHILKMFPSTWYV");
        let original = protein_query.clone();
        let mut query_blk = BlastSequenceBlk {
            sequence: Some(protein_query.clone()),
            sequence_start: Some(protein_query),
            ..BlastSequenceBlk::default()
        };
        let mut query_info = QueryInfo::new_blastp(&[20]);
        let qsup_options = QuerySetUpOptions {
            filtering_options: Some(crate::filter::SBlastFilterOptions {
                mask_at_hash: false,
                dust_options: None,
                seg_options: None,
                repeat_filter_options: None,
                window_masker_options: None,
                read_quality_options: None,
            }),
            filter_string: None,
            strand_option: 0,
            genetic_code: 1,
        };
        let scoring = ScoringOptions {
            matrix_path: None,
            reward: 0,
            penalty: 0,
            gap_open: 11,
            gap_extend: 1,
            shift_pen: i16::MAX as i32,
            gapped_calculation: true,
            complexity_adjusted_scoring: false,
            matrix_name: Some("BLOSUM62".to_string()),
            is_ooframe: false,
            program_number: crate::program::UNDEFINED,
        };
        let mut lookup_segments = None;
        let mut mask = None;
        let mut sbp = None;

        assert_eq!(
            blast_main_set_up(
                crate::program::BLASTP,
                Some(&qsup_options),
                Some(&scoring),
                Some(&mut query_blk),
                Some(&mut query_info),
                Some(&mut lookup_segments),
                Some(&mut mask),
                Some(&mut sbp),
                1.0,
                None,
            ),
            0
        );
        assert_eq!(query_blk.sequence.as_deref(), Some(original.as_slice()));
        assert_eq!(
            lookup_segments.as_ref().map(|node| node.ssl),
            Some(SSeqRange { left: 0, right: 19 })
        );
        assert_eq!(mask.as_ref().map(|mask| mask.masks.len()), Some(1));
        assert!(mask.as_ref().is_some_and(|mask| mask.masks[0].is_none()));
        assert!(sbp.as_ref().is_some_and(|sbp| sbp.kbp_gap[0].lambda > 0.0));
    }

    fn blastn_word_params() -> InitialWordParameters {
        InitialWordParameters {
            options: InitialWordOptions::new_blastn(),
            x_dropoff_max: 0,
            cutoff_score_min: 0,
            cutoffs: vec![BlastUngappedCutoffs::default()],
            ungapped_extension: true,
            nucl_score_table: InitialWordParameters::build_nucl_score_table(1, -3),
        }
    }

    #[test]
    fn setup_aux_structures_selects_nucleotide_callbacks_and_allocates_scratch() {
        let seq_src = SetupSeqSrc {
            seqs: vec![b"ACGTACGT".to_vec()],
            total_stats: 0,
            num_stats: 0,
        };
        let lookup = crate::lookup::LookupTableWrap::SmallNa(crate::lookup::SmallNaLookupTable {
            word_length: 7,
            backbone: vec![0; 16],
            overflow: Vec::new(),
            pv_array: Vec::new(),
            longest_chain: 7,
            scan_step: 1,
        });
        let word_params = blastn_word_params();
        let ext_options = ExtensionOptions::new_blastn();
        let hit_options = HitSavingOptions::default();
        let query = BlastSequenceBlk {
            length: 32,
            sequence: Some(vec![0; 32]),
            ..BlastSequenceBlk::default()
        };
        let query_info = QueryInfo::new_blastn(&[32]);
        let mut aux = None;

        assert_eq!(
            s_blast_set_up_aux_structures(
                Some(&seq_src),
                Some(&lookup),
                Some(&word_params),
                Some(&ext_options),
                Some(&hit_options),
                Some(&query),
                Some(&query_info),
                Some(&mut aux),
                false,
            ),
            0
        );

        let aux = aux.expect("aux structures");
        assert_eq!(aux.core.WordFinder, BlastAuxWordFinder::Nucleotide);
        assert_eq!(aux.core.GetGappedScore, BlastAuxGappedPath::Standard);
        assert_eq!(aux.core.JumperGapped, JumperGappedType::None);
        assert_eq!(
            aux.core.offset_pairs.len(),
            crate::lookup::get_offset_array_size(&lookup) as usize
        );
        assert!(aux
            .core
            .ewp
            .as_ref()
            .is_some_and(|ewp| ewp.diag_table.is_some()));
        assert!(aux.core.init_hitlist.is_some());
        assert!(aux.na_extend_choice.is_some_and(
            |choice| choice.lookup_callback == Some(crate::lookup::NaLookupCallback::SmallNa)
        ));
        assert!(aux.core.mapper_wordhits.is_none());
    }

    #[test]
    fn setup_aux_structures_jumper_uses_mapper_hits_for_many_queries() {
        let seq_src = SetupSeqSrc {
            seqs: vec![b"ACGTACGT".to_vec()],
            total_stats: 0,
            num_stats: 0,
        };
        let lookup = crate::lookup::LookupTableWrap::Megablast(crate::lookup::MbLookupTable {
            word_length: 18,
            lut_word_length: 11,
            discontiguous: false,
            template_length: 18,
            template_type: crate::lookup::DiscTemplateType::Contiguous,
            two_templates: false,
            second_template_type: crate::lookup::DiscTemplateType::Contiguous,
            hashtable: Vec::new(),
            hashtable2: Vec::new(),
            next_pos: Vec::new(),
            next_pos2: Vec::new(),
            pv_array: Vec::new(),
            pv_array_bts: 0,
            longest_chain: 3,
            scan_step: 18,
        });
        let word_params = blastn_word_params();
        let mut ext_options = ExtensionOptions::new_blastn();
        ext_options.prelim_gap_ext = PrelimGapExt::JumperWithTraceback;
        let hit_options = HitSavingOptions::default();
        let query = BlastSequenceBlk {
            length: 5000,
            sequence: Some(vec![0; 5000]),
            ..BlastSequenceBlk::default()
        };
        let query_lengths = vec![5; 1001];
        let query_info = QueryInfo::new_blastn(&query_lengths);
        let mut aux = None;

        assert_eq!(
            s_blast_set_up_aux_structures(
                Some(&seq_src),
                Some(&lookup),
                Some(&word_params),
                Some(&ext_options),
                Some(&hit_options),
                Some(&query),
                Some(&query_info),
                Some(&mut aux),
                true,
            ),
            0
        );

        let aux = aux.expect("aux structures");
        assert_eq!(aux.core.WordFinder, BlastAuxWordFinder::None);
        assert_eq!(aux.core.GetGappedScore, BlastAuxGappedPath::Standard);
        assert_eq!(aux.core.JumperGapped, JumperGappedType::ShortReadIndexed);
        let mapper = aux.core.mapper_wordhits.expect("mapper word hits");
        assert_eq!(mapper.num_arrays, 10);
        assert_eq!(mapper.divisor, 501);
    }

    fn protein_kbp() -> KarlinBlk {
        KarlinBlk {
            lambda: 0.267,
            k: 0.041,
            log_k: 0.041_f64.ln(),
            h: 0.140,
            round_down: false,
        }
    }

    #[test]
    fn calc_eff_lengths_blastp_basic() {
        // Single 100-aa query against a 1 Mb / 1000-seq DB with BLOSUM62 11/1.
        let mut qi = QueryInfo {
            num_queries: 1,
            contexts: vec![ContextInfo {
                query_offset: 0,
                query_length: 100,
                eff_searchsp: 0,
                length_adjustment: 0,
                query_index: 0,
                frame: 0,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: 100,
            min_length: 0,
        };
        let scoring = ScoringOptions {
            matrix_path: None,
            reward: 0,
            penalty: 0,
            gap_open: 11,
            gap_extend: 1,
            shift_pen: i16::MAX as i32,
            gapped_calculation: true,
            complexity_adjusted_scoring: false,
            matrix_name: Some("BLOSUM62".to_string()),
            is_ooframe: false,
            program_number: crate::program::UNDEFINED,
        };
        let eff = EffectiveLengthsParameters {
            options: EffectiveLengthsOptions::default(),
            real_db_length: 1_000_000,
            real_num_seqs: 1000,
        };
        let kbp = vec![protein_kbp()];
        let rc = blast_calc_eff_lengths(
            crate::program::BLASTP,
            &scoring,
            &eff,
            &kbp,
            &kbp,
            "BLOSUM62",
            &mut qi,
        );
        assert_eq!(rc, 0);
        assert!(qi.contexts[0].length_adjustment > 0);
        assert!(qi.contexts[0].eff_searchsp > 0);
    }

    #[test]
    fn calc_eff_lengths_returns_error_for_valid_context_missing_kbp() {
        let mut qi = QueryInfo {
            num_queries: 1,
            contexts: vec![ContextInfo {
                query_offset: 0,
                query_length: 100,
                eff_searchsp: 0,
                length_adjustment: 0,
                query_index: 0,
                frame: 0,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: 100,
            min_length: 0,
        };
        let scoring = ScoringOptions {
            matrix_path: None,
            reward: 0,
            penalty: 0,
            gap_open: 11,
            gap_extend: 1,
            shift_pen: i16::MAX as i32,
            gapped_calculation: true,
            complexity_adjusted_scoring: false,
            matrix_name: Some("BLOSUM62".to_string()),
            is_ooframe: false,
            program_number: crate::program::UNDEFINED,
        };
        let eff = EffectiveLengthsParameters {
            options: EffectiveLengthsOptions::default(),
            real_db_length: 1_000_000,
            real_num_seqs: 1000,
        };

        let rc = blast_calc_eff_lengths(
            crate::program::BLASTP,
            &scoring,
            &eff,
            &[],
            &[],
            "BLOSUM62",
            &mut qi,
        );

        assert_eq!(rc, -1);
    }

    #[test]
    fn calc_eff_lengths_translated_subject_div3() {
        // tblastn: db_length should be divided by 3 internally.
        let mut qi = QueryInfo {
            num_queries: 1,
            contexts: vec![ContextInfo {
                query_offset: 0,
                query_length: 100,
                eff_searchsp: 0,
                length_adjustment: 0,
                query_index: 0,
                frame: 0,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: 100,
            min_length: 0,
        };
        let scoring = ScoringOptions {
            matrix_path: None,
            reward: 0,
            penalty: 0,
            gap_open: 11,
            gap_extend: 1,
            shift_pen: i16::MAX as i32,
            gapped_calculation: true,
            complexity_adjusted_scoring: false,
            matrix_name: Some("BLOSUM62".to_string()),
            is_ooframe: false,
            program_number: crate::program::UNDEFINED,
        };
        let eff = EffectiveLengthsParameters {
            options: EffectiveLengthsOptions::default(),
            real_db_length: 3_000_000,
            real_num_seqs: 1,
        };
        let kbp = vec![protein_kbp()];
        let rc = blast_calc_eff_lengths(
            crate::program::TBLASTN,
            &scoring,
            &eff,
            &kbp,
            &kbp,
            "BLOSUM62",
            &mut qi,
        );
        assert_eq!(rc, 0);
        let eff_searchsp_translated = qi.contexts[0].eff_searchsp;
        assert!(eff_searchsp_translated > 0);

        // Same DB but as blastp (no translation): db_length is 3x larger,
        // so search space should be ~3x larger too.
        qi.contexts[0].eff_searchsp = 0;
        qi.contexts[0].length_adjustment = 0;
        let rc2 = blast_calc_eff_lengths(
            crate::program::BLASTP,
            &scoring,
            &eff,
            &kbp,
            &kbp,
            "BLOSUM62",
            &mut qi,
        );
        assert_eq!(rc2, 0);
        let eff_searchsp_full = qi.contexts[0].eff_searchsp;
        // Roughly 3x; exact ratio is between 2 and 3.5 depending on length
        // adjustment (the larger db_length implies a slightly larger
        // length_adjustment too, so the ratio is below 3).
        let ratio = eff_searchsp_full as f64 / eff_searchsp_translated as f64;
        assert!(ratio > 2.0 && ratio < 3.5, "ratio = {ratio}");
    }

    #[test]
    fn calc_eff_lengths_searchsp_override_takes_precedence() {
        let mut qi = QueryInfo {
            num_queries: 1,
            contexts: vec![ContextInfo {
                query_offset: 0,
                query_length: 100,
                eff_searchsp: 0,
                length_adjustment: 0,
                query_index: 0,
                frame: 0,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: 100,
            min_length: 0,
        };
        let scoring = ScoringOptions {
            matrix_path: None,
            reward: 0,
            penalty: 0,
            gap_open: 11,
            gap_extend: 1,
            shift_pen: i16::MAX as i32,
            gapped_calculation: true,
            complexity_adjusted_scoring: false,
            matrix_name: Some("BLOSUM62".to_string()),
            is_ooframe: false,
            program_number: crate::program::UNDEFINED,
        };
        let eff = EffectiveLengthsParameters {
            options: EffectiveLengthsOptions {
                db_length: 0,
                dbseq_num: 0,
                num_searchspaces: 1,
                searchsp_eff: vec![987_654_321],
            },
            real_db_length: 1_000_000,
            real_num_seqs: 1000,
        };
        let kbp = vec![protein_kbp()];
        let rc = blast_calc_eff_lengths(
            crate::program::BLASTP,
            &scoring,
            &eff,
            &kbp,
            &kbp,
            "BLOSUM62",
            &mut qi,
        );
        assert_eq!(rc, 0);
        assert_eq!(qi.contexts[0].eff_searchsp, 987_654_321);
    }

    #[test]
    fn calc_eff_lengths_clears_stale_adjustment_for_invalid_context() {
        let mut qi = QueryInfo {
            num_queries: 1,
            contexts: vec![
                ContextInfo {
                    query_offset: 0,
                    query_length: 100,
                    eff_searchsp: 0,
                    length_adjustment: 0,
                    query_index: 0,
                    frame: 0,
                    is_valid: true,
                    segment_flags: crate::queryinfo::E_NO_SEGMENTS,
                },
                ContextInfo {
                    query_offset: 100,
                    query_length: 100,
                    eff_searchsp: 123,
                    length_adjustment: 77,
                    query_index: 0,
                    frame: 0,
                    is_valid: false,
                    segment_flags: crate::queryinfo::E_NO_SEGMENTS,
                },
            ],
            max_length: 100,
            min_length: 0,
        };
        let scoring = ScoringOptions {
            matrix_path: None,
            reward: 0,
            penalty: 0,
            gap_open: 11,
            gap_extend: 1,
            shift_pen: i16::MAX as i32,
            gapped_calculation: true,
            complexity_adjusted_scoring: false,
            matrix_name: Some("BLOSUM62".to_string()),
            is_ooframe: false,
            program_number: crate::program::UNDEFINED,
        };
        let eff = EffectiveLengthsParameters {
            options: EffectiveLengthsOptions::default(),
            real_db_length: 1_000_000,
            real_num_seqs: 1000,
        };
        let kbp = vec![protein_kbp(), protein_kbp()];

        assert_eq!(
            blast_calc_eff_lengths(
                crate::program::BLASTP,
                &scoring,
                &eff,
                &kbp,
                &kbp,
                "BLOSUM62",
                &mut qi,
            ),
            0
        );
        assert!(qi.contexts[0].length_adjustment > 0);
        assert_eq!(qi.contexts[1].eff_searchsp, 0);
        assert_eq!(qi.contexts[1].length_adjustment, 0);
    }

    #[test]
    fn subject_totals_match_c_fallback_order() {
        let src = SetupSeqSrc {
            seqs: vec![b"AAAA".to_vec(), b"CC".to_vec()],
            total_stats: 100,
            num_stats: 5,
        };
        let mut total = 0;
        let mut num = 0;
        blast_get_subject_totals(Some(&src), Some(&mut total), Some(&mut num));
        assert_eq!((total, num), (100, 5));

        let src = SetupSeqSrc {
            seqs: vec![b"AAAA".to_vec(), b"CC".to_vec()],
            total_stats: 0,
            num_stats: 0,
        };
        blast_get_subject_totals(Some(&src), Some(&mut total), Some(&mut num));
        assert_eq!((total, num), (6, 2));

        let empty = SetupSeqSrc {
            seqs: Vec::new(),
            total_stats: 0,
            num_stats: 0,
        };
        blast_get_subject_totals(Some(&empty), Some(&mut total), Some(&mut num));
        assert_eq!((total, num), (-1, -1));

        total = 7;
        num = 3;
        blast_get_subject_totals(None, Some(&mut total), Some(&mut num));
        assert_eq!((total, num), (-1, -1));

        let mut arg = crate::seqsrc::BlastSeqSrcGetSeqArg {
            oid: 0,
            encoding: SeqEncoding::Protein.into(),
            ..crate::seqsrc::BlastSeqSrcGetSeqArg::default()
        };
        let data = crate::seqsrc::blast_seq_src_get_sequence(Some(&src), Some(&mut arg))
            .expect("sequence");
        assert_eq!(data.length, 4);
        assert!(arg.seq.is_some());
        assert!(crate::seqsrc::blast_seq_src_get_sequence(None, Some(&mut arg)).is_none());
    }

    #[test]
    fn calc_eff_lengths_c_wrapper_preserves_null_status() {
        let mut qi = QueryInfo {
            num_queries: 1,
            contexts: vec![ContextInfo {
                query_offset: 0,
                query_length: 100,
                eff_searchsp: 0,
                length_adjustment: 0,
                query_index: 0,
                frame: 0,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: 100,
            min_length: 0,
        };
        let scoring = ScoringOptions {
            matrix_path: None,
            reward: 0,
            penalty: 0,
            gap_open: 11,
            gap_extend: 1,
            shift_pen: i16::MAX as i32,
            gapped_calculation: true,
            complexity_adjusted_scoring: false,
            matrix_name: Some("BLOSUM62".to_string()),
            is_ooframe: false,
            program_number: crate::program::UNDEFINED,
        };
        let eff = EffectiveLengthsParameters {
            options: EffectiveLengthsOptions::default(),
            real_db_length: 1_000_000,
            real_num_seqs: 1000,
        };
        let kbp = vec![protein_kbp()];

        assert_eq!(
            blast_calc_eff_lengths_c(
                crate::program::BLASTP,
                Some(&scoring),
                Some(&eff),
                Some(&kbp),
                Some(&kbp),
                Some("BLOSUM62"),
                Some(&mut qi),
            ),
            0
        );
        assert!(qi.contexts[0].eff_searchsp > 0);
        assert_eq!(
            blast_calc_eff_lengths_c(
                crate::program::BLASTP,
                Some(&scoring),
                Some(&eff),
                Some(&kbp),
                Some(&kbp),
                Some("BLOSUM62"),
                None,
            ),
            -1
        );
        assert_eq!(
            blast_calc_eff_lengths_c(
                crate::program::BLASTP,
                None,
                Some(&eff),
                Some(&kbp),
                Some(&kbp),
                Some("BLOSUM62"),
                Some(&mut qi),
            ),
            -1
        );
    }

    #[test]
    fn one_subject_update_parameters_recomputes_dependent_cutoffs() {
        let mut qi = QueryInfo {
            num_queries: 1,
            contexts: vec![ContextInfo {
                query_offset: 0,
                query_length: 100,
                eff_searchsp: 0,
                length_adjustment: 0,
                query_index: 0,
                frame: 0,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: 100,
            min_length: 0,
        };
        let scoring = ScoringOptions {
            matrix_path: None,
            reward: 0,
            penalty: 0,
            gap_open: 11,
            gap_extend: 1,
            shift_pen: i16::MAX as i32,
            gapped_calculation: true,
            complexity_adjusted_scoring: false,
            matrix_name: Some("BLOSUM62".to_string()),
            is_ooframe: false,
            program_number: crate::program::UNDEFINED,
        };
        let kbp = protein_kbp();
        let mut sbp =
            crate::stat::blast_score_blk_new(crate::encoding::BLASTAA_SEQ_CODE, 1).expect("sbp");
        sbp.name = Some("BLOSUM62".to_string());
        sbp.kbp = vec![kbp.clone()];
        sbp.kbp_std = vec![kbp.clone()];
        sbp.kbp_gap = vec![kbp.clone()];
        sbp.kbp_gap_std = vec![kbp.clone()];

        let mut hit = HitSavingParameters {
            options: HitSavingOptions::default(),
            cutoff_score_min: 0,
            low_score: vec![0],
            cutoffs: vec![BlastGappedCutoffs::default()],
            link_hsp_params: Some(LinkHSPParameters::default()),
            prelim_evalue: 0.0,

            ..Default::default()
        };
        let mut word = InitialWordParameters {
            options: InitialWordOptions::new_blastp(),
            x_dropoff_max: 0,
            cutoff_score_min: 0,
            cutoffs: vec![BlastUngappedCutoffs::default()],
            ungapped_extension: true,
            nucl_score_table: InitialWordParameters::build_nucl_score_table(1, -3),
        };
        let mut eff = EffectiveLengthsParameters {
            options: EffectiveLengthsOptions::default(),
            real_db_length: 0,
            real_num_seqs: 1,
        };

        assert_eq!(
            blast_one_subject_update_parameters(
                crate::program::BLASTP,
                20_000,
                Some(&scoring),
                Some(&mut qi),
                Some(&sbp),
                Some(&mut hit),
                Some(&mut word),
                Some(&mut eff),
            ),
            0
        );
        assert_eq!(eff.real_db_length, 20_000);
        assert!(qi.contexts[0].eff_searchsp > 0);
        assert!(hit.cutoff_score_min > 0);
        assert!(word.cutoff_score_min > 0);
        assert_eq!(
            hit.link_hsp_params
                .as_ref()
                .expect("link params")
                .cutoff_small_gap,
            word.cutoff_score_min
        );

        assert_eq!(
            blast_one_subject_update_parameters(
                crate::program::BLASTP,
                20_000,
                None,
                Some(&mut qi),
                Some(&sbp),
                Some(&mut hit),
                None,
                Some(&mut eff),
            ),
            -1
        );
    }

    #[test]
    fn set_phi_pattern_info_populates_pattern_metadata_and_query_adjustment() {
        let pattern_blk = sphi_pattern_search_blk_new("A-B", false, None).expect("pattern");
        let mut query_info = QueryInfo::new_blastp(&[6]);
        let lookup = [SSeqRange { left: 0, right: 5 }];
        let query = crate::encoding::encode_ncbistdaa_sequence(b"XXABAB");
        let mut pattern_info = None;

        assert_eq!(
            blast_set_phi_pattern_info(
                crate::program::PHI_BLASTP,
                Some(&pattern_blk),
                Some(&query),
                Some(&lookup),
                Some(&mut query_info),
                Some(&mut pattern_info),
            ),
            0
        );
        let pattern_info = pattern_info.expect("pattern info");
        assert_eq!(pattern_info.num_patterns, 2);
        assert_eq!(pattern_info.occurrences[0].offset, 2);
        assert_eq!(pattern_info.occurrences[1].offset, 4);
        assert_eq!(pattern_info.pattern.as_deref(), Some("A-B"));
        assert!((pattern_info.probability - pattern_blk.pattern_probability).abs() < 1e-12);
        assert_eq!(
            query_info.contexts[0].length_adjustment,
            pattern_blk.min_pattern_match_length
        );

        let mut missing_info = None;
        assert_eq!(
            blast_set_phi_pattern_info(
                crate::program::PHI_BLASTP,
                Some(&pattern_blk),
                Some(&crate::encoding::encode_ncbistdaa_sequence(b"XXXXXX")),
                Some(&lookup),
                Some(&mut query_info),
                Some(&mut missing_info),
            ),
            -1
        );
        assert!(missing_info
            .as_ref()
            .is_some_and(|info| info.num_patterns == 0));
        assert_eq!(
            blast_set_phi_pattern_info(
                crate::program::PHI_BLASTP,
                None,
                Some(&query),
                Some(&lookup),
                Some(&mut query_info),
                Some(&mut missing_info),
            ),
            -1
        );
    }
}
