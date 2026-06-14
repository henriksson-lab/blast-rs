//! Rust equivalent of `blast_setup.c` — search-space and effective-length
//! calculation, matrix-info initialisation, etc. Currently ports only
//! `BLAST_CalcEffLengths` and its private helpers; the rest of the file
//! will be filled in incrementally.

use crate::algo::blast::core::blast_seqsrc::{
    blast_seq_src_get_max_seq_len, blast_seq_src_get_min_seq_len, blast_seq_src_get_num_seqs,
    blast_seq_src_get_num_seqs_stats, blast_seq_src_get_seq_len, blast_seq_src_get_tot_len,
    blast_seq_src_get_tot_len_stats, BlastSeqSrc,
};
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
use crate::stat::{
    blast_compute_length_adjustment, blast_get_nucl_alpha_beta, blast_gumbel_blk_calc,
    blast_karlin_blk_gapped_calc, blast_karlin_blk_nucl_gapped_calc,
    blast_score_blk_kbp_ideal_calc, blast_score_blk_kbp_ungapped_calc, blast_score_blk_matrix_fill,
    blast_score_blk_new, blast_score_freq_new, blast_score_set_ambig_res, BlastScoreBlk, KarlinBlk,
    BLAST_GAP_EXTN_MEGABLAST, BLAST_GAP_OPEN_MEGABLAST,
};
use crate::util::{BlastSequenceBlk, SSeqRange};

/// `BLAST_REWARD` (`blast_options.h:152`).
pub const BLAST_REWARD: i32 = 1;
/// `BLAST_PENALTY` (`blast_options.h:151`).
pub const BLAST_PENALTY: i32 = -3;

/// NCBI: BLAST_GapAlignSetUp (`blast_setup.c:888`).
#[allow(clippy::too_many_arguments)]
pub fn blast_gap_align_set_up(
    program_number: ProgramType,
    seq_src: Option<&BlastSeqSrc>,
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
            let mut oid = 0_i32;
            total_length = blast_seq_src_get_seq_len(seq_src, &mut oid) as i64;
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

    let matrix = options
        .matrix_name
        .as_deref()
        .unwrap_or("BLOSUM62")
        .to_ascii_uppercase();
    let values = match matrix.as_str() {
        "BLOSUM62" => match (options.gap_open, options.gap_extend) {
            (11, 1) => Some((0.270, 0.047)),
            (9, 2) => Some((0.285, 0.075)),
            (8, 2) => Some((0.265, 0.046)),
            (7, 2) => Some((0.243, 0.032)),
            (12, 1) => Some((0.281, 0.057)),
            (10, 1) => Some((0.250, 0.033)),
            _ => None,
        },
        "PAM30" => match (options.gap_open, options.gap_extend) {
            (9, 1) => Some((0.295, 0.13)),
            (7, 2) => Some((0.306, 0.15)),
            (6, 2) => Some((0.292, 0.13)),
            (5, 2) => Some((0.263, 0.077)),
            (10, 1) => Some((0.309, 0.15)),
            (8, 1) => Some((0.270, 0.070)),
            _ => None,
        },
        "PAM70" => match (options.gap_open, options.gap_extend) {
            (10, 1) => Some((0.291, 0.089)),
            (8, 2) => Some((0.303, 0.13)),
            (7, 2) => Some((0.287, 0.095)),
            (6, 2) => Some((0.269, 0.079)),
            (11, 1) => Some((0.307, 0.13)),
            (9, 1) => Some((0.269, 0.058)),
            _ => None,
        },
        "BLOSUM80" => match (options.gap_open, options.gap_extend) {
            (10, 1) => Some((0.300, 0.072)),
            (8, 2) => Some((0.308, 0.089)),
            (7, 2) => Some((0.295, 0.077)),
            (6, 2) => Some((0.271, 0.051)),
            (11, 1) => Some((0.314, 0.096)),
            (9, 1) => Some((0.277, 0.046)),
            _ => None,
        },
        "BLOSUM45" => match (options.gap_open, options.gap_extend) {
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
        _ => return -2,
    };
    let Some((lambda, k)) = values else {
        return -1;
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
    seqsrc: Option<&BlastSeqSrc>,
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
        let mut oid = 0_i32;
        *total_length = blast_seq_src_get_seq_len(seqsrc, &mut oid) as i64;
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
pub fn blast_get_alpha_beta(
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
    use std::any::Any;
    use std::sync::Arc;

    use crate::algo::blast::core::blast_encoding::EBlastEncoding;
    use crate::algo::blast::core::blast_seqsrc::{
        blast_seq_src_get_sequence, blast_seq_src_release_sequence, BlastSeqSrcGetSeqArg,
        BLAST_SEQSRC_ERROR, BLAST_SEQSRC_SUCCESS,
    };
    use crate::algo::blast::core::blast_util::BLAST_SequenceBlk;
    use crate::link_hsps::LinkHSPParameters;
    use crate::options::{HitSavingOptions, InitialWordOptions};
    use crate::parameters::{BlastGappedCutoffs, BlastUngappedCutoffs};
    use crate::pattern::sphi_pattern_search_blk_new;
    use crate::queryinfo::ContextInfo;

    struct SetupSeqSrc {
        seqs: Vec<Vec<u8>>,
        total_stats: i64,
        num_stats: i32,
    }

    fn setup_seq_src(seqs: Vec<Vec<u8>>, total_stats: i64, num_stats: i32) -> BlastSeqSrc {
        BlastSeqSrc {
            new_fn_ptr: None,
            delete_fn_ptr: None,
            copy_fn_ptr: None,
            set_number_of_threads: None,
            get_num_seqs: Some(setup_num_seqs),
            get_num_seqs_stats: Some(setup_num_seqs_stats),
            get_max_seq_len: Some(setup_max_seq_len),
            get_min_seq_len: Some(setup_min_seq_len),
            get_avg_seq_len: Some(setup_avg_seq_len),
            get_tot_len: Some(setup_tot_len),
            get_tot_len_stats: Some(setup_tot_len_stats),
            get_name: Some(setup_name),
            get_is_prot: Some(setup_is_prot),
            get_supports_partial_fetching: None,
            set_seq_range: None,
            get_sequence: Some(setup_get_sequence),
            get_seq_len: Some(setup_seq_len),
            release_sequence: Some(setup_release_sequence),
            iter_next: None,
            reset_chunk_iterator: None,
            data_structure: Some(Arc::new(SetupSeqSrc {
                seqs,
                total_stats,
                num_stats,
            })),
            init_error_str: None,
            get_gis: None,
        }
    }

    fn setup_data(data: &dyn Any) -> &SetupSeqSrc {
        data.downcast_ref::<SetupSeqSrc>()
            .expect("SetupSeqSrc callback data")
    }

    fn setup_num_seqs(data: &dyn Any, _arg: Option<&mut dyn Any>) -> i32 {
        setup_data(data).seqs.len() as i32
    }

    fn setup_num_seqs_stats(data: &dyn Any, _arg: Option<&mut dyn Any>) -> i32 {
        setup_data(data).num_stats
    }

    fn setup_tot_len(data: &dyn Any, _arg: Option<&mut dyn Any>) -> i64 {
        setup_data(data)
            .seqs
            .iter()
            .map(|seq| seq.len() as i64)
            .sum()
    }

    fn setup_tot_len_stats(data: &dyn Any, _arg: Option<&mut dyn Any>) -> i64 {
        setup_data(data).total_stats
    }

    fn setup_max_seq_len(data: &dyn Any, _arg: Option<&mut dyn Any>) -> i32 {
        setup_data(data)
            .seqs
            .iter()
            .map(|seq| seq.len() as i32)
            .max()
            .unwrap_or(0)
    }

    fn setup_min_seq_len(data: &dyn Any, _arg: Option<&mut dyn Any>) -> i32 {
        setup_data(data)
            .seqs
            .iter()
            .map(|seq| seq.len() as i32)
            .min()
            .unwrap_or(0)
    }

    fn setup_avg_seq_len(data: &dyn Any, _arg: Option<&mut dyn Any>) -> i32 {
        let data = setup_data(data);
        if data.seqs.is_empty() {
            0
        } else {
            (data.seqs.iter().map(|seq| seq.len() as i64).sum::<i64>() / data.seqs.len() as i64)
                as i32
        }
    }

    fn setup_name(_data: &dyn Any, _arg: Option<&mut dyn Any>) -> Option<String> {
        Some("setup-test".to_string())
    }

    fn setup_is_prot(_data: &dyn Any, _arg: Option<&mut dyn Any>) -> bool {
        true
    }

    fn setup_seq_len(data: &dyn Any, oid: Option<&mut dyn Any>) -> i32 {
        let Some(oid) = oid.and_then(|oid| oid.downcast_ref::<i32>()) else {
            return -1;
        };
        setup_data(data)
            .seqs
            .get((*oid).max(0) as usize)
            .map_or(-1, |seq| seq.len() as i32)
    }

    fn setup_get_sequence(data: &dyn Any, arg: &mut BlastSeqSrcGetSeqArg) -> i16 {
        let Some(sequence) = setup_data(data).seqs.get(arg.oid.max(0) as usize).cloned() else {
            arg.seq = std::ptr::null_mut();
            return BLAST_SEQSRC_ERROR;
        };
        arg.seq = Box::into_raw(Box::new(BLAST_SequenceBlk {
            sequence: Some(sequence),
            sequence_start: None,
            length: sequence.len() as i32,
            frame: 0,
            subject_strand: 0,
            oid: arg.oid,
            sequence_allocated: true,
            sequence_start_allocated: false,
            sequence_start_nomask: None,
            sequence_nomask: None,
            nomask_allocated: false,
            oof_sequence: None,
            oof_sequence_allocated: false,
            compressed_nuc_seq: None,
            compressed_nuc_seq_start: None,
            lcase_mask: None,
            lcase_mask_allocated: false,
            chunk: 0,
            gen_code_string: None,
            seq_ranges: None,
            num_seq_ranges: 0,
            seq_ranges_allocated: false,
            mask_type: crate::algo::blast::core::blast_util::SubjectMaskingType::NoSubjMasking,
            bases_offset: 0,
        }));
        BLAST_SEQSRC_SUCCESS
    }

    fn setup_release_sequence(_data: &dyn Any, arg: &mut BlastSeqSrcGetSeqArg) {
        if !arg.seq.is_null() {
            unsafe {
                drop(Box::from_raw(arg.seq));
            }
            arg.seq = std::ptr::null_mut();
        }
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
}
