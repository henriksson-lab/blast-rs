//! Rust equivalent of `blast_setup.c` — search-space and effective-length
//! calculation, matrix-info initialisation, etc. Currently ports only
//! `BLAST_CalcEffLengths` and its private helpers; the rest of the file
//! will be filled in incrementally.

use crate::options::{EffectiveLengthsOptions, ScoringOptions};
use crate::parameters::EffectiveLengthsParameters;
use crate::program::{
    is_mapping, is_phi_blast, subject_is_translated, ProgramType, BLASTN,
};
use crate::queryinfo::QueryInfo;
use crate::stat::{
    compute_length_adjustment_exact, lookup_protein_params, nucl_alpha_beta, KarlinBlk,
};

/// `BLAST_REWARD` (`blast_options.h:152`).
pub const BLAST_REWARD: i32 = 1;
/// `BLAST_PENALTY` (`blast_options.h:151`).
pub const BLAST_PENALTY: i32 = -3;

/// Port of `s_GetEffectiveSearchSpaceForContext` (`blast_setup.c:676`).
fn get_effective_search_space_for_context(
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

/// Port of `BlastEffectiveLengthsOptions_IsSearchSpaceSet`
/// (`blast_options.c:1019`). Returns `true` iff any
/// `searchsp_eff[0..num_searchspaces]` entry is non-zero.
pub fn is_search_space_set(options: &EffectiveLengthsOptions) -> bool {
    if options.searchsp_eff.is_empty() {
        return false;
    }
    options
        .searchsp_eff
        .iter()
        .take(options.num_searchspaces as usize)
        .any(|&s| s != 0)
}

/// Port of `BLAST_GetAlphaBeta` (`blast_stat.c:3094`) for the matrices we
/// currently support (BLOSUM62 — extend `lookup_protein_params` for more).
/// Used by `blast_calc_eff_lengths` for protein-on-protein scoring.
fn blast_get_alpha_beta(
    matrix_name: &str,
    gap_open: i32,
    gap_extend: i32,
    gapped: bool,
    kbp_ungapped: &KarlinBlk,
) -> (f64, f64) {
    if !gapped {
        // Ungapped path in NCBI uses `alpha_arr[0]` / `beta_arr[0]` for
        // the matrix; if no entry exists fall back to Lambda/H.
        return (kbp_ungapped.lambda / kbp_ungapped.h, 0.0);
    }
    if matrix_name.eq_ignore_ascii_case("BLOSUM62") || matrix_name.is_empty() {
        if let Some(p) = lookup_protein_params(gap_open, gap_extend) {
            return (p.alpha, p.beta);
        }
    }
    (kbp_ungapped.lambda / kbp_ungapped.h, 0.0)
}

/// 1-1 port of NCBI `BLAST_CalcEffLengths` (`blast_setup.c:699`). Updates
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
///   (`sbp->kbp_std`). Used only by `nucl_alpha_beta` /
///   `blast_get_alpha_beta` for the alpha/beta lookup.
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

    if db_length == 0 && !is_search_space_set(eff_len_options) {
        return 0;
    }

    if subject_is_translated(program_number) {
        db_length /= 3;
    }

    let db_num_seqs = if eff_len_options.dbseq_num > 0 {
        eff_len_options.dbseq_num
    } else {
        eff_len_params.real_num_seqs
    };

    if is_mapping(program_number) {
        for ctx in &mut query_info.contexts {
            ctx.eff_searchsp = db_length;
        }
        return 0;
    }

    if is_phi_blast(program_number) {
        for ctx in &mut query_info.contexts {
            let eff = db_length - (db_num_seqs as i64) * (ctx.length_adjustment as i64);
            ctx.eff_searchsp = eff;
        }
        return 0;
    }

    for (index, ctx) in query_info.contexts.iter_mut().enumerate() {
        let mut effective_search_space =
            get_effective_search_space_for_context(eff_len_options, index);

        let kbp = match kbp_array.get(index) {
            Some(k) => k,
            None => {
                ctx.eff_searchsp = effective_search_space;
                continue;
            }
        };
        let kbp_std = match kbp_std_array.get(index) {
            Some(k) => k,
            None => {
                ctx.eff_searchsp = effective_search_space;
                continue;
            }
        };

        let query_length = ctx.query_length;
        if !ctx.is_valid || query_length <= 0 {
            ctx.eff_searchsp = effective_search_space;
            continue;
        }

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
            nucl_alpha_beta(
                rew,
                pen,
                scoring_options.gap_open,
                scoring_options.gap_extend,
                kbp_std.lambda,
                kbp_std.h,
                scoring_options.gapped_calculation,
            )
        } else {
            blast_get_alpha_beta(
                matrix_name,
                scoring_options.gap_open,
                scoring_options.gap_extend,
                scoring_options.gapped_calculation,
                kbp_std,
            )
        };

        let (length_adjustment, _converged) = compute_length_adjustment_exact(
            kbp.k,
            kbp.log_k,
            alpha / kbp.lambda,
            beta,
            query_length,
            db_length,
            db_num_seqs,
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::queryinfo::ContextInfo;

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
            }],
            max_length: 100,
        };
        let scoring = ScoringOptions {
            reward: 0,
            penalty: 0,
            gap_open: 11,
            gap_extend: 1,
            gapped_calculation: true,
            matrix_name: Some("BLOSUM62".to_string()),
            is_ooframe: false,
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
            }],
            max_length: 100,
        };
        let scoring = ScoringOptions {
            reward: 0,
            penalty: 0,
            gap_open: 11,
            gap_extend: 1,
            gapped_calculation: true,
            matrix_name: Some("BLOSUM62".to_string()),
            is_ooframe: false,
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
            }],
            max_length: 100,
        };
        let scoring = ScoringOptions {
            reward: 0,
            penalty: 0,
            gap_open: 11,
            gap_extend: 1,
            gapped_calculation: true,
            matrix_name: Some("BLOSUM62".to_string()),
            is_ooframe: false,
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
}
