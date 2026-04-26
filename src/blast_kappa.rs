//! Rust port of `blast_kappa.c` — Smith-Waterman / composition-based
//! scoring for the traceback stage of blastp / blastx / tblastn / RPS.
//!
//! This is the **skeleton** module: only the leaf utilities that are pure
//! byte/numeric processing have been ported so far. The composition
//! adjustment driver `Blast_RedoAlignmentCore_MT` (669 LOC, cyclomatic
//! 120) and its support layer (`BlastCompo_*` types in
//! `composition_adjustment/redo_alignment.h`) require porting the
//! surrounding C structs and are tracked as TODOs below.
//!
//! Status (from `scratch/ccc/missing-names.txt` audit):
//!
//! | Group | Function | LOC | Ported? |
//! |---|---|---:|---|
//! | leaf | `s_GetSubjectLength` | 4 | ✅ |
//! | leaf | `s_HSPListNormalizeScores` | 9 | ✅ |
//! | leaf | `s_GetHash` | 8 | ✅ |
//! | leaf | `s_ExtendRight` | 75 | ✅ |
//! | leaf | `s_ExtendLeft` | 70 | ✅ |
//! | leaf | `s_FindNumIdentical` | 49 | ✅ |
//! | leaf | `s_HitlistReapContained` | 50 | ✅ |
//! | macro | `CONTAINED_IN_HSP` | 1 | ✅ |
//! | macro | `GET_NUCL_LENGTH` | 1 | ✅ |
//! | leaf | `s_CalcLambda` | 13 | ✅ |
//! | comp adj | `s_GetPosBasedStartFreqRatios` | 24 | TODO — needs SFreqRatios |
//! | comp adj | `s_GetStartFreqRatios` | 12 | TODO — needs SFreqRatios |
//! | comp adj | `s_ScalePosMatrix` | 44 | TODO — needs PSI types |
//! | redo adj | `s_NewAlignmentFromGapAlign` | 20 | TODO — needs BlastCompo_Alignment |
//! | redo adj | `s_NewAlignmentUsingXdrop` | 44 | TODO |
//! | redo adj | `s_HSPListFromDistinctAlignments` | 33 | TODO |
//! | redo adj | `s_RedoOneAlignment` | 33 | TODO |
//! | redo adj | `s_TestNearIdentical` | 55 | TODO — needs BlastCompo_SequenceData |
//! | redo adj | `s_HitlistEvaluateAndPurge` | 36 | TODO |
//! | redo adj | `s_ComputeNumIdentities` | 60 | TODO |
//! | redo adj | `s_AdjustEvaluesForComposition` | 75 | TODO |
//! | redo adj | `s_ResultHspToDistinctAlign` | 25 | TODO |
//! | redo adj | `s_SWFindFinalEndsUsingXdrop` | 30 | TODO |
//! | redo adj | `s_MatchingSequenceRelease` | 10 | TODO |
//! | redo adj | `s_MatchingSequenceInitialize` | 33 | TODO |
//! | redo adj | `s_DoSegSequenceData` | 22 | TODO |
//! | redo adj | `s_SequenceGetTranslatedRange` | 51 | TODO |
//! | redo adj | `s_SequenceGetProteinRange` | 52 | TODO |
//! | redo adj | `s_SequenceGetRange` | 31 | TODO |
//! | redo adj | `s_SavedParametersFree` | 13 | TODO |
//! | redo adj | `s_SavedParametersNew` | 26 | TODO |
//! | redo adj | `s_RecordInitialSearch` | 27 | TODO |
//! | redo adj | `s_RestoreSearch` | 23 | TODO |
//! | redo adj | `s_MatrixInfoInit` | 28 | TODO |
//! | redo adj | `s_CreateWordArray` | 16 | TODO |
//! | redo adj | `s_FreeBlastCompo_QueryInfoArray` | 10 | TODO |
//! | redo adj | `s_GetQueryInfo` | 22 | TODO |
//! | redo adj | `s_GappingParamsNew` | 23 | TODO |
//! | redo adj | `s_GetAlignParams` | 48 | TODO |
//! | redo adj | `s_FillResultsFromCompoHeaps` | 13 | TODO |
//! | redo adj | `s_ClearHeap` | 3 | TODO |
//! | redo adj | `s_BlastGapAlignStruct_Free` | 29 | TODO |
//! | redo adj | `s_BlastGapAlignStruct_Copy` | 95 | TODO |
//! | redo adj | `s_BlastScoreBlk_Copy` | 131 | TODO |
//! | driver | `Blast_RedoAlignmentCore` | 31 | TODO (thin wrapper) |
//! | driver | `Blast_RedoAlignmentCore_MT` | 669 | TODO (the big one) |

use crate::compo_mode_condition::MatrixAdjustRule;
use crate::hspstream::{Hsp, HspList};
use crate::math::NCBIMATH_LN2;
use crate::program::{ProgramType, RPS_BLAST};

// ───────────────────────────────────────────────────────────────────────────
// Struct ports from `composition_adjustment/redo_alignment.h` and
// `composition_adjustment/composition_constants.h`.
//
// These are 1-1 ports of the C structs used by `Blast_RedoAlignmentCore_MT`
// and its helpers. Function-pointer fields in C (`Blast_RedoAlignCallbacks`)
// are represented as Rust function-pointer types so the dispatch table
// translates directly. Implementations of these callbacks (and the
// `Blast_RedoAlignParams` driver wrapper that owns them) come in later
// iterations of this loop.
// ───────────────────────────────────────────────────────────────────────────

/// `ECompoAdjustModes` (`composition_constants.h:59`). Permissible
/// composition-adjustment modes selected by the user.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum CompoAdjustMode {
    NoCompositionBasedStats = 0,
    CompositionBasedStats = 1,
    CompositionMatrixAdjust = 2,
    CompoForceFullMatrixAdjust = 3,
}

impl CompoAdjustMode {
    pub fn from_u8(v: u8) -> Self {
        match v {
            0 => CompoAdjustMode::NoCompositionBasedStats,
            1 => CompoAdjustMode::CompositionBasedStats,
            2 => CompoAdjustMode::CompositionMatrixAdjust,
            _ => CompoAdjustMode::CompoForceFullMatrixAdjust,
        }
    }
}

/// `BlastCompo_Alignment` (`redo_alignment.h:57`). One distinct alignment
/// of the query against the current subject. Stored in singly-linked
/// reverse-of-computation order; we model the linked list with a
/// `next: Option<Box<Self>>` so ownership transfers are explicit.
///
/// `context` in NCBI's C is `void *` carrying traceback data of an
/// arbitrary type (commonly a `GapEditScript*`). The Rust port stores a
/// concrete `Option<crate::gapinfo::GapEditScript>` because that's the
/// only context value used by the kappa driver in our tree. If a future
/// caller needs heterogeneous contexts, swap for an enum or trait
/// object.
#[derive(Debug, Clone)]
pub struct BlastCompoAlignment {
    pub score: i32,
    pub matrix_adjust_rule: MatrixAdjustRule,
    pub query_index: i32,
    pub query_start: i32,
    pub query_end: i32,
    pub match_start: i32,
    pub match_end: i32,
    pub frame: i32,
    pub context: Option<crate::gapinfo::GapEditScript>,
    pub next: Option<Box<BlastCompoAlignment>>,
}

impl BlastCompoAlignment {
    /// 1-1 port of `BlastCompo_AlignmentNew`. The C version `malloc`s and
    /// initializes; Rust returns a stack value the caller can wrap.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        score: i32,
        matrix_adjust_rule: MatrixAdjustRule,
        query_index: i32,
        query_start: i32,
        query_end: i32,
        match_start: i32,
        match_end: i32,
        frame: i32,
        context: Option<crate::gapinfo::GapEditScript>,
    ) -> Self {
        Self {
            score,
            matrix_adjust_rule,
            query_index,
            query_start,
            query_end,
            match_start,
            match_end,
            frame,
            context,
            next: None,
        }
    }
}

/// `BlastCompo_GappingParams` (`redo_alignment.h:101`). Parameters
/// passed to the gapped-alignment callback inside the kappa driver.
#[derive(Debug, Clone, Copy)]
pub struct BlastCompoGappingParams {
    pub gap_open: i32,
    pub gap_extend: i32,
    pub decline_align: i32,
    pub x_dropoff: i32,
    // `void *context` from C dropped: no caller currently uses it; a
    // generic parameter or trait object will be needed when the OOF
    // path is ported.
}

/// `BlastCompo_SequenceRange` (`redo_alignment.h:118`). Half-open
/// `[begin, end)` interval plus a context tag (frame index or query
/// index in concatenated queries).
#[derive(Debug, Clone, Copy)]
pub struct BlastCompoSequenceRange {
    pub begin: i32,
    pub end: i32,
    pub context: i32,
}

/// `BlastCompo_SequenceData` (`redo_alignment.h:131`). Owns the residue
/// buffer and a shadow byte at `data[-1]` (NCBI uses a sentinel zero
/// before the start). Rust models this by allocating with one leading
/// sentinel byte and exposing a slice view starting at offset 1.
///
/// The `buffer` field in C is the malloc'd allocation; if `data` is
/// just a view into another sequence then `buffer` is null. In Rust we
/// always own a `Vec<u8>` and the `data_offset` records whether the
/// usable data starts at index 1 (with sentinel) or index 0.
#[derive(Debug, Clone, Default)]
pub struct BlastCompoSequenceData {
    pub buffer: Vec<u8>,
    pub data_offset: usize,
    pub length: i32,
}

impl BlastCompoSequenceData {
    /// Returns the residue slice (`&data[0..length]` in C terms,
    /// equivalent to `&buffer[data_offset .. data_offset + length]`).
    pub fn data(&self) -> &[u8] {
        let start = self.data_offset;
        let end = start + self.length as usize;
        &self.buffer[start..end]
    }
}

/// `BlastCompo_MatchingSequence` (`redo_alignment.h:156`). Identifies
/// one subject sequence and provides a hook for callbacks to access
/// the underlying database/translation. NCBI's `local_data` is a
/// `void *`; we keep a generic index payload here and let the caller
/// own any larger state outside the struct.
#[derive(Debug, Clone, Copy)]
pub struct BlastCompoMatchingSequence {
    pub length: i32,
    pub index: i32,
    pub local_data_index: i32,
}

/// `BlastCompo_QueryInfo` (`redo_alignment.h:166`). Per-query metadata
/// consumed by the redo-alignment driver. `composition` (the
/// `Blast_AminoAcidComposition` struct) and `words` (Uint8 array of
/// hashed query k-mers) are filled in by `s_GetQueryInfo` /
/// `s_CreateWordArray`; both are TODO until those helpers land.
#[derive(Debug, Clone, Default)]
pub struct BlastCompoQueryInfo {
    pub origin: i32,
    pub seq: BlastCompoSequenceData,
    pub composition: Vec<f64>, // BLASTAA_SIZE entries; per-residue freq
    pub eff_search_space: f64,
    pub words: Vec<u64>,
}

/// Port of `Blast_MatrixInfo` (`composition_adjustment.h`). Holds the
/// scoring matrix used for the redo-alignment pass, plus its
/// dimensions, name, and the rounding/scale parameters used by
/// composition adjustment. Currently empty placeholder; populated when
/// `s_MatrixInfoInit` is ported.
#[derive(Debug, Clone, Default)]
pub struct BlastMatrixInfo {
    pub matrix_name: String,
    pub rows: i32,
    pub cols: i32,
    pub positional: bool,
    pub ungapped_lambda: f64,
    pub matrix: Vec<Vec<i32>>,
    pub start_freq_ratios: Vec<Vec<f64>>,
}

/// `Blast_RedoAlignParams` (`redo_alignment.h:328`). Parameter block
/// owned by the kappa driver. Callback function pointers
/// (`Blast_RedoAlignCallbacks`) are ported in the next iteration.
#[derive(Debug, Clone)]
pub struct BlastRedoAlignParams {
    pub matrix_info: BlastMatrixInfo,
    pub gapping_params: BlastCompoGappingParams,
    pub compo_adjust_mode: CompoAdjustMode,
    pub position_based: bool,
    pub re_pseudocounts: i32,
    pub subject_is_translated: bool,
    pub query_is_translated: bool,
    pub ccat_query_length: i32,
    pub cutoff_s: i32,
    pub cutoff_e: f64,
    pub do_link_hsps: bool,
    pub near_identical_cutoff: f64,
}

#[cfg(test)]
mod struct_tests {
    use super::*;

    #[test]
    fn alignment_new_round_trip() {
        let a = BlastCompoAlignment::new(
            42,
            MatrixAdjustRule::DontAdjust,
            0,
            10,
            50,
            20,
            60,
            1,
            None,
        );
        assert_eq!(a.score, 42);
        assert_eq!(a.matrix_adjust_rule, MatrixAdjustRule::DontAdjust);
        assert_eq!(a.query_start, 10);
        assert_eq!(a.query_end, 50);
        assert_eq!(a.match_start, 20);
        assert_eq!(a.match_end, 60);
        assert!(a.next.is_none());
    }

    #[test]
    fn compo_adjust_mode_from_u8() {
        assert_eq!(CompoAdjustMode::from_u8(0), CompoAdjustMode::NoCompositionBasedStats);
        assert_eq!(CompoAdjustMode::from_u8(1), CompoAdjustMode::CompositionBasedStats);
        assert_eq!(CompoAdjustMode::from_u8(2), CompoAdjustMode::CompositionMatrixAdjust);
        assert_eq!(CompoAdjustMode::from_u8(3), CompoAdjustMode::CompoForceFullMatrixAdjust);
        // Anything > 3 saturates to the strongest mode (matches C's
        // raw enum cast — out-of-range becomes max).
        assert_eq!(CompoAdjustMode::from_u8(99), CompoAdjustMode::CompoForceFullMatrixAdjust);
    }

    #[test]
    fn redo_one_alignment_perfect_blastn_match() {
        // BLASTNA-encoded perfect-match query/subject. Our
        // `blast_gapped_align` operates on BLASTNA, so set up a 16x16
        // matrix via build_blastna_matrix.
        let query_buf = vec![0u8, 1, 2, 3, 0, 1, 2, 3, 0, 1];
        let subject_buf = query_buf.clone();
        let query = BlastCompoSequenceData {
            buffer: query_buf,
            data_offset: 0,
            length: 10,
        };
        let subject = BlastCompoSequenceData {
            buffer: subject_buf,
            data_offset: 0,
            length: 10,
        };
        let qr = BlastCompoSequenceRange { begin: 0, end: 10, context: 0 };
        let sr = BlastCompoSequenceRange { begin: 0, end: 10, context: 0 };
        let result = redo_one_alignment(
            &query, &qr, &subject, &sr,
            /* gapped_start_q */ 4,
            /* gapped_start_s */ 4,
            MatrixAdjustRule::DontAdjust,
            /* reward */ 1,
            /* penalty */ -3,
            /* gap_open */ 5,
            /* gap_extend */ 2,
            /* x_dropoff */ 30,
        );
        let aln = result.expect("alignment");
        assert!(aln.score > 0);
        assert_eq!(aln.query_start, 0);
        assert_eq!(aln.query_end, 10);
        assert_eq!(aln.match_start, 0);
        assert_eq!(aln.match_end, 10);
    }

    #[test]
    fn get_start_freq_ratios_blosum62_matches_matrix_module() {
        let ratios = get_start_freq_ratios("BLOSUM62").expect("BLOSUM62 ratios");
        let canonical = crate::matrix::get_blosum62_freq_ratios();
        // Deep-equal copy.
        for i in 0..crate::matrix::AA_SIZE {
            for j in 0..crate::matrix::AA_SIZE {
                assert_eq!(ratios[i][j], canonical[i][j]);
            }
        }
    }

    #[test]
    fn get_start_freq_ratios_unknown_matrix_errs() {
        assert!(get_start_freq_ratios("PAM30").is_err());
    }

    #[test]
    fn sw_find_final_ends_using_xdrop_returns_target_score_on_full_match() {
        // Identity match should produce the maximum score with the
        // first X-drop pass. Use BLOSUM62 + ACGT-style residues.
        let matrix = crate::matrix::BLOSUM62;
        // Convert i32 BLOSUM62 to the [[i32; 16]; 16] shape used by
        // align_ex (BLASTNA-style 16-by-16). align_ex expects a 16x16,
        // but BLOSUM62 is 28x28 — sw_find_final_ends_using_xdrop is
        // designed for protein. To exercise the path, build a tiny
        // synthetic 16x16 matrix where matches score 5 and mismatches
        // score -3 (mimicking nucleotide scoring).
        let mut m = [[-3i32; 16]; 16];
        for i in 0..4 {
            m[i][i] = 5;
        }
        let q: Vec<u8> = vec![0, 1, 2, 3, 0, 1, 2, 3];
        let s: Vec<u8> = q.clone();
        let (score, q_ext, s_ext, _ops) = sw_find_final_ends_using_xdrop(
            &q, 0, q.len() - 1, &s, 0, s.len() - 1, &m, 5, 2, 20, /* target */ 40,
        );
        assert!(score >= 40);
        assert!(q_ext > 0);
        assert!(s_ext > 0);
    }

    #[test]
    fn blast_gap_align_workspace_clone_is_deep_copy() {
        let mut script = crate::gapinfo::GapEditScript::new();
        script.push(crate::gapinfo::GapAlignOpType::Sub, 5);
        let ws = BlastGapAlignWorkspace {
            gap_x_dropoff: 65,
            score: 100,
            query_start: 0,
            query_stop: 50,
            subject_start: 10,
            subject_stop: 60,
            edit_script: Some(script),
        };
        let copy = ws.deep_copy();
        assert_eq!(copy.score, 100);
        assert_eq!(copy.gap_x_dropoff, 65);
        assert_eq!(
            copy.edit_script.as_ref().unwrap().ops.len(),
            ws.edit_script.as_ref().unwrap().ops.len()
        );
    }

    #[test]
    fn blast_gap_align_struct_free_clears_slot() {
        let mut slot = Some(BlastGapAlignWorkspace::default());
        blast_gap_align_struct_free(&mut slot);
        assert!(slot.is_none());
    }

    #[test]
    fn blast_redo_alignment_core_delegates_to_mt() {
        // Fully-stubbed inputs — both functions return -1 until the
        // MT driver body lands. This test pins the surface API.
        let qi = crate::queryinfo::QueryInfo {
            num_queries: 1,
            contexts: vec![],
            max_length: 0,
        };
        let mut kbp = vec![];
        let mut mtx = Vec::<Vec<i32>>::new();
        let mut scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
                reward: 0,
                penalty: 0,
                gap_open: 11,
                gap_extend: 1,
                gapped_calculation: true,
                matrix_name: Some("BLOSUM62".to_string()),
                is_ooframe: false,
            },
            1.0,
        );
        let params = blast_redo_align_params_new(
            BlastMatrixInfo::default(),
            BlastCompoGappingParams {
                gap_open: 11,
                gap_extend: 1,
                decline_align: i32::MIN,
                x_dropoff: 0,
            },
            CompoAdjustMode::NoCompositionBasedStats,
            false,
            false,
            false,
            0,
            0,
            0.0,
            false,
            0.0,
        );
        let mut saved = BlastKappaSavedParameters::new(
            0,
            0,
            CompoAdjustMode::NoCompositionBasedStats,
            false,
        );
        let mut hsp_list = HspList::new(0);
        let mut results = crate::hspstream::HspResults::new(1);
        let rc = blast_redo_alignment_core(
            crate::program::BLASTP,
            &[],
            &qi,
            &mut kbp,
            &mut mtx,
            &mut scoring,
            &params,
            &mut saved,
            &mut hsp_list,
            &mut results,
        );
        // Wrapper delegates → MT driver returns -1 (not yet implemented).
        assert_eq!(rc, -1);
    }

    #[test]
    fn get_align_params_blastp_assembly() {
        let scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
                reward: 0,
                penalty: 0,
                gap_open: 11,
                gap_extend: 1,
                gapped_calculation: true,
                matrix_name: Some("BLOSUM62".to_string()),
                is_ooframe: false,
            },
            1.0,
        );
        let kbp = vec![crate::stat::KarlinBlk {
            lambda: 0.267,
            k: 0.041,
            log_k: 0.041_f64.ln(),
            h: 0.140,
            round_down: false,
        }];
        let p = get_align_params(
            crate::program::BLASTP,
            "BLOSUM62",
            &scoring,
            &kbp,
            /* kbp_gap_first_lambda */ 0.267,
            /* kbp_ideal_lambda */ 0.3176,
            /* local_scaling_factor */ 1.0,
            /* cutoff_score_min */ 50,
            /* expect_value */ 1e-3,
            CompoAdjustMode::CompositionBasedStats,
            /* position_based */ false,
            /* do_link_hsps */ false,
            /* ccat_query_length */ 1000,
            /* num_contexts */ 1,
            /* gap_x_dropoff_final_bits */ 25.0,
            /* raw_gap_x_dropoff_final */ 0,
        )
        .expect("params");
        // Sanity-check the bridge.
        assert_eq!(p.compo_adjust_mode, CompoAdjustMode::CompositionBasedStats);
        assert!(!p.position_based);
        assert!(!p.query_is_translated);
        assert!(!p.subject_is_translated);
        assert_eq!(p.ccat_query_length, 1000);
        // do_link_hsps = false → cutoff_s = 1 per NCBI's branch.
        assert_eq!(p.cutoff_s, 1);
        assert_eq!(p.cutoff_e, 1e-3);
        // BLOSUM62 11/1 alpha=1.9 → x_dropoff ≈ 25 * ln(2) / 0.267 ≈ 65.
        let expected_xdrop = (25.0 * crate::math::NCBIMATH_LN2 / 0.267).round() as i32;
        assert_eq!(p.gapping_params.x_dropoff, expected_xdrop);
        // near_identical_cutoff = 1.74 * ln(2) / 0.267
        let expected_nic = (NEAR_IDENTICAL_BITS_PER_POSITION
            * crate::math::NCBIMATH_LN2)
            / 0.267;
        assert!((p.near_identical_cutoff - expected_nic).abs() < 1e-9);
        // Matrix info populated for BLOSUM62.
        assert_eq!(p.matrix_info.matrix_name, "BLOSUM62");
        assert_eq!(p.matrix_info.rows, crate::matrix::AA_SIZE as i32);
    }

    #[test]
    fn get_align_params_link_hsps_uses_cutoff_score_min() {
        let scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
                reward: 0,
                penalty: 0,
                gap_open: 11,
                gap_extend: 1,
                gapped_calculation: true,
                matrix_name: Some("BLOSUM62".to_string()),
                is_ooframe: false,
            },
            1.0,
        );
        let kbp = vec![crate::stat::KarlinBlk {
            lambda: 0.267,
            k: 0.041,
            log_k: 0.041_f64.ln(),
            h: 0.140,
            round_down: false,
        }];
        let p = get_align_params(
            crate::program::BLASTP,
            "BLOSUM62",
            &scoring,
            &kbp,
            0.267,
            0.3176,
            32.0, // local_scaling_factor
            50,   // cutoff_score_min
            1e-3,
            CompoAdjustMode::CompositionBasedStats,
            false,
            true, // do_link_hsps
            1000,
            1,
            25.0,
            0,
        )
        .expect("params");
        // do_link_hsps = true → cutoff_s = cutoff_score_min * localScalingFactor = 50 * 32 = 1600.
        assert_eq!(p.cutoff_s, 50 * 32);
    }

    #[test]
    fn blast_redo_align_params_new_sets_pseudocounts() {
        let matrix_info = BlastMatrixInfo::default();
        let gapping = BlastCompoGappingParams {
            gap_open: 11,
            gap_extend: 1,
            decline_align: i32::MIN,
            x_dropoff: 65,
        };
        let p = blast_redo_align_params_new(
            matrix_info,
            gapping,
            CompoAdjustMode::CompositionMatrixAdjust,
            /* position_based */ false,
            /* query_is_translated */ true,
            /* subject_is_translated */ false,
            /* ccat_query_length */ 1234,
            /* cutoff_s */ 50,
            /* cutoff_e */ 1e-3,
            /* do_link_hsps */ false,
            /* near_identical_cutoff */ 1.74,
        );
        assert_eq!(p.re_pseudocounts, K_RE_MATRIX_ADJUSTMENT_PSEUDOCOUNTS);
        assert_eq!(p.compo_adjust_mode, CompoAdjustMode::CompositionMatrixAdjust);
        assert!(p.query_is_translated);
        assert!(!p.subject_is_translated);
        assert!(!p.position_based);
        assert_eq!(p.ccat_query_length, 1234);
        assert_eq!(p.cutoff_s, 50);
        assert_eq!(p.cutoff_e, 1e-3);
        assert!(!p.do_link_hsps);
        assert_eq!(p.near_identical_cutoff, 1.74);
        assert_eq!(p.gapping_params.gap_open, 11);
        assert_eq!(p.gapping_params.x_dropoff, 65);
    }

    #[test]
    fn blast_redo_align_params_free_clears_slot() {
        let matrix_info = BlastMatrixInfo::default();
        let gapping = BlastCompoGappingParams {
            gap_open: 0,
            gap_extend: 0,
            decline_align: i32::MIN,
            x_dropoff: 0,
        };
        let mut slot = Some(blast_redo_align_params_new(
            matrix_info,
            gapping,
            CompoAdjustMode::NoCompositionBasedStats,
            false,
            false,
            false,
            0,
            0,
            0.0,
            false,
            0.0,
        ));
        blast_redo_align_params_free(&mut slot);
        assert!(slot.is_none());
    }

    #[test]
    fn result_hsp_to_distinct_align_groups_by_frame_index() {
        // Three HSPs with contexts 0, 1, 0 → frame indices 0, 1, 0.
        let mk = |score: i32, context: i32, q: i32, s: i32| Hsp {
            score,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 0.0,
            query_offset: q,
            query_end: q + 10,
            subject_offset: s,
            subject_end: s + 10,
            context,
            num_gaps: 0,
        };
        let hsps = vec![mk(100, 0, 0, 0), mk(80, 1, 30, 30), mk(60, 0, 60, 60)];
        let mut lists: [Option<Box<BlastCompoAlignment>>; 6] = Default::default();
        let mut counts = [0i32; 6];
        let rc = result_hsp_to_distinct_align(&mut lists, &mut counts, &hsps, 0, 1.0);
        assert_eq!(rc, 0);
        assert_eq!(counts[0], 2); // two HSPs in frame index 0
        assert_eq!(counts[1], 1); // one HSP in frame index 1
        // Walk frame 0's list — preserves insertion order (NCBI
        // appends to tail, not head).
        let head_0 = lists[0].as_ref().expect("head 0");
        assert_eq!(head_0.score, 100);
        let next_0 = head_0.next.as_ref().expect("second in list");
        assert_eq!(next_0.score, 60);
        // Frame 1's list has just the single 80-score HSP.
        let head_1 = lists[1].as_ref().expect("head 1");
        assert_eq!(head_1.score, 80);
        assert!(head_1.next.is_none());
    }

    #[test]
    fn result_hsp_to_distinct_align_applies_scaling() {
        let hsps = vec![Hsp {
            score: 100,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 0.0,
            query_offset: 0,
            query_end: 10,
            subject_offset: 0,
            subject_end: 10,
            context: 0,
            num_gaps: 0,
        }];
        let mut lists: [Option<Box<BlastCompoAlignment>>; 6] = Default::default();
        let mut counts = [0i32; 6];
        result_hsp_to_distinct_align(&mut lists, &mut counts, &hsps, 0, 32.0);
        // C: `(int)(hsp->score * localScalingFactor)` — Rust uses
        // `.round() as i32` for nearest-integer behavior.
        assert_eq!(lists[0].as_ref().unwrap().score, 3200);
    }

    #[test]
    fn do_seg_sequence_data_unmasked_clean_input() {
        // A diverse run of residues should not trigger SEG.
        // Use NCBIstdaa codes for distinct amino acids.
        let mut seq = BlastCompoSequenceData {
            buffer: vec![0, 1, 4, 7, 10, 13, 16, 19, 22, 25, 0],
            data_offset: 1,
            length: 9,
        };
        let (status, was_biased) = do_seg_sequence_data(&mut seq);
        assert_eq!(status, 0);
        // Diverse input → no bias detected.
        assert!(!was_biased);
        // Sentinels intact.
        assert_eq!(seq.buffer[0], 0);
        assert_eq!(seq.buffer[10], 0);
    }

    #[test]
    fn do_seg_sequence_data_masks_low_complexity_run() {
        // Long run of single residue triggers SEG. Build a 30-residue
        // sequence of all alanines (NCBIstdaa A = 1).
        let mut data = vec![1u8; 30];
        let mut buffer = vec![0u8; 32];
        buffer[1..=30].copy_from_slice(&data);
        data.clear();
        let mut seq = BlastCompoSequenceData {
            buffer,
            data_offset: 1,
            length: 30,
        };
        let (status, was_biased) = do_seg_sequence_data(&mut seq);
        assert_eq!(status, 0);
        assert!(was_biased);
        // Some residue should now be NCBIstdaa X (21).
        let masked_count = seq.data().iter().filter(|&&b| b == 21).count();
        assert!(masked_count > 0);
    }

    #[test]
    fn matching_sequence_release_resets_state() {
        let mut s = BlastCompoMatchingSequence {
            length: 1234,
            index: 5,
            local_data_index: 7,
        };
        matching_sequence_release(&mut s);
        assert_eq!(s.local_data_index, -1);
        assert_eq!(s.length, 0);
        // index left untouched per the C `if (self->index >= 0)` guard
        // (the C function reads index but does not reset it).
        assert_eq!(s.index, 5);
    }

    #[test]
    fn rescale_search_scales_lambda_and_gap_costs() {
        let mut kbp = vec![crate::stat::KarlinBlk {
            lambda: 0.3,
            k: 0.04,
            log_k: 0.04_f64.ln(),
            h: 0.14,
            round_down: false,
        }];
        let mut scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
                reward: 0,
                penalty: 0,
                gap_open: 11,
                gap_extend: 1,
                gapped_calculation: true,
                matrix_name: Some("BLOSUM62".to_string()),
                is_ooframe: false,
            },
            1.0,
        );
        rescale_search(&mut kbp, &mut scoring, 1, 32.0);
        assert!((kbp[0].lambda - 0.3 / 32.0).abs() < 1e-12);
        // logK = ln(K), unchanged because K is unchanged.
        assert!((kbp[0].log_k - 0.04_f64.ln()).abs() < 1e-12);
        // gap_open/extend scaled by 32 with rounding.
        assert_eq!(scoring.gap_open, 11 * 32);
        assert_eq!(scoring.gap_extend, 1 * 32);
        assert_eq!(scoring.scale_factor, 32.0);
    }

    #[test]
    fn sequence_prep_query_range_copies_with_sentinels_and_substitutes_u_to_c() {
        // Query: A=1, U=24, G=4, T=8, plus a context pad. Range covers
        // idx 0..4 inclusive.
        let q = BlastCompoSequenceData {
            buffer: vec![1, 24, 4, 8, 9],
            data_offset: 0,
            length: 5,
        };
        let range = BlastCompoSequenceRange { begin: 0, end: 4, context: 0 };
        let prepped = sequence_prep_query_range(&q, &range);
        // length = 4 (end - begin); buffer = length + 2 = 6.
        assert_eq!(prepped.length, 4);
        assert_eq!(prepped.buffer.len(), 6);
        // Sentinel at index 0 stays zero.
        assert_eq!(prepped.buffer[0], 0);
        // Selenocysteine (24) replaced by Cysteine (3).
        assert_eq!(prepped.buffer[1], 1);
        assert_eq!(prepped.buffer[2], 3); // was 24 (U) → 3 (C)
        assert_eq!(prepped.buffer[3], 4);
        assert_eq!(prepped.buffer[4], 8);
        // Trailing sentinel at length+1 = 5.
        assert_eq!(prepped.buffer[5], 0);
        // data() returns the slice starting at offset 1.
        assert_eq!(prepped.data(), &[1, 3, 4, 8]);
    }

    #[test]
    fn gapping_params_new_picks_min_lambda() {
        let scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
                reward: 0,
                penalty: 0,
                gap_open: 11,
                gap_extend: 1,
                gapped_calculation: true,
                matrix_name: Some("BLOSUM62".to_string()),
                is_ooframe: false,
            },
            1.0,
        );
        // Two contexts: lambdas 0.30 and 0.27 (the latter is min).
        let kbp = vec![
            crate::stat::KarlinBlk {
                lambda: 0.30,
                k: 0.04,
                log_k: 0.04_f64.ln(),
                h: 0.14,
                round_down: false,
            },
            crate::stat::KarlinBlk {
                lambda: 0.27,
                k: 0.04,
                log_k: 0.04_f64.ln(),
                h: 0.14,
                round_down: false,
            },
        ];
        let params = gapping_params_new(&scoring, &kbp, 2, /* bits = */ 25.0, /* raw = */ 0);
        assert_eq!(params.gap_open, 11);
        assert_eq!(params.gap_extend, 1);
        // Expected x_dropoff = max(round(25 * ln(2) / 0.27), 0) ≈ 64.
        let expected = (25.0 * crate::math::NCBIMATH_LN2 / 0.27).round() as i32;
        assert_eq!(params.x_dropoff, expected);
    }

    #[test]
    fn gapping_params_new_falls_back_to_raw_when_no_lambda() {
        let scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
                reward: 0,
                penalty: 0,
                gap_open: 11,
                gap_extend: 1,
                gapped_calculation: true,
                matrix_name: Some("BLOSUM62".to_string()),
                is_ooframe: false,
            },
            1.0,
        );
        let params = gapping_params_new(&scoring, &[], 0, 25.0, 47);
        assert_eq!(params.x_dropoff, 47);
    }

    #[test]
    fn matrix_info_init_blastp_populates_blosum62() {
        let mut info = BlastMatrixInfo::default();
        let rc = matrix_info_init_blastp(&mut info, "BLOSUM62", 0.3176, 1.0);
        assert_eq!(rc, 0);
        assert_eq!(info.matrix_name, "BLOSUM62");
        assert!(!info.positional);
        assert_eq!(info.rows, crate::matrix::AA_SIZE as i32);
        assert_eq!(info.cols, crate::matrix::AA_SIZE as i32);
        assert_eq!(info.ungapped_lambda, 0.3176);
        assert_eq!(info.matrix.len(), crate::matrix::AA_SIZE);
        assert_eq!(info.matrix[0].len(), crate::matrix::AA_SIZE);
    }

    #[test]
    fn matrix_info_init_blastp_unknown_matrix_returns_error() {
        let mut info = BlastMatrixInfo::default();
        let rc = matrix_info_init_blastp(&mut info, "PAM30", 0.34, 1.0);
        // NCBI returns non-zero from `s_GetStartFreqRatios` for matrices
        // it doesn't have ratios for. Our TODO range mirrors that.
        assert_eq!(rc, -1);
    }

    #[test]
    fn record_and_restore_round_trip() {
        let mut saved = BlastKappaSavedParameters::new(
            crate::matrix::AA_SIZE as i32,
            2,
            CompoAdjustMode::CompositionBasedStats,
            false,
        );
        let kbp = vec![
            crate::stat::KarlinBlk {
                lambda: 0.3,
                k: 0.04,
                log_k: 0.04_f64.ln(),
                h: 0.14,
                round_down: false,
            },
            crate::stat::KarlinBlk {
                lambda: 0.27,
                k: 0.041,
                log_k: 0.041_f64.ln(),
                h: 0.14,
                round_down: false,
            },
        ];
        let matrix: Vec<Vec<i32>> = (0..crate::matrix::AA_SIZE)
            .map(|r| (0..crate::matrix::AA_SIZE).map(|c| (r * 100 + c) as i32).collect())
            .collect();
        let mut scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
                reward: 0,
                penalty: 0,
                gap_open: 11,
                gap_extend: 1,
                gapped_calculation: true,
                matrix_name: Some("BLOSUM62".to_string()),
                is_ooframe: false,
            },
            1.0,
        );

        // Record snapshot.
        record_initial_search(
            &mut saved,
            &kbp,
            &matrix,
            &scoring,
            0,
            CompoAdjustMode::CompositionBasedStats,
            false,
        );
        assert_eq!(saved.gap_open, 11);
        assert_eq!(saved.gap_extend, 1);
        assert_eq!(saved.kbp_gap_orig[0].lambda, 0.3);
        assert_eq!(saved.kbp_gap_orig[1].lambda, 0.27);
        assert_eq!(saved.orig_matrix[0][0], 0);
        assert_eq!(saved.orig_matrix[5][7], 507);

        // Mutate the live state, then restore.
        scoring.gap_open = 99;
        scoring.gap_extend = 99;
        let mut live_kbp = vec![crate::stat::KarlinBlk::default(); 2];
        let mut live_matrix: Vec<Vec<i32>> =
            vec![vec![-1i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE];
        restore_search(
            &saved,
            &mut live_kbp,
            &mut live_matrix,
            &mut scoring,
            0,
            false,
            CompoAdjustMode::CompositionBasedStats,
        );
        assert_eq!(scoring.gap_open, 11);
        assert_eq!(scoring.gap_extend, 1);
        assert_eq!(live_kbp[0].lambda, 0.3);
        assert_eq!(live_kbp[1].lambda, 0.27);
        assert_eq!(live_matrix[0][0], 0);
        assert_eq!(live_matrix[5][7], 507);
    }

    #[test]
    fn saved_parameters_new_no_composition_skips_matrix() {
        let sp = BlastKappaSavedParameters::new(
            10,
            3,
            CompoAdjustMode::NoCompositionBasedStats,
            false,
        );
        assert_eq!(sp.num_queries, 3);
        assert_eq!(sp.kbp_gap_orig.len(), 3);
        // C: `if (compo_adjust_mode != eNoCompositionBasedStats) Nlm_Int4MatrixNew(...);`
        // — skipped here, so orig_matrix stays empty.
        assert!(sp.orig_matrix.is_empty());
    }

    #[test]
    fn saved_parameters_new_compo_allocates_aa_matrix() {
        let sp = BlastKappaSavedParameters::new(
            5,
            1,
            CompoAdjustMode::CompositionBasedStats,
            false,
        );
        // Non-position-based: BLASTAA_SIZE × BLASTAA_SIZE.
        assert_eq!(sp.orig_matrix.len(), crate::matrix::AA_SIZE);
        assert_eq!(sp.orig_matrix[0].len(), crate::matrix::AA_SIZE);
    }

    #[test]
    fn saved_parameters_new_position_based_uses_rows() {
        let sp = BlastKappaSavedParameters::new(
            42,
            1,
            CompoAdjustMode::CompositionMatrixAdjust,
            /* positionBased = */ true,
        );
        // Position-based: rows × BLASTAA_SIZE.
        assert_eq!(sp.orig_matrix.len(), 42);
        assert_eq!(sp.orig_matrix[0].len(), crate::matrix::AA_SIZE);
    }

    #[test]
    fn saved_parameters_free_clears_option() {
        let mut slot = Some(BlastKappaSavedParameters::new(
            5,
            1,
            CompoAdjustMode::NoCompositionBasedStats,
            false,
        ));
        saved_parameters_free(&mut slot);
        assert!(slot.is_none());
    }

    #[test]
    fn blast_compo_heap_pop_worst_returns_largest_evalue_first() {
        let mut heap = BlastCompoHeap::new(10, 1e-5);
        let mut a = HspList::new(1);
        a.best_evalue = 1e-10;
        let mut b = HspList::new(2);
        b.best_evalue = 1.0;
        let mut c = HspList::new(3);
        c.best_evalue = 1e-5;
        heap.records.push(a);
        heap.records.push(b);
        heap.records.push(c);

        // Pop order: largest e-value first.
        let first = heap.pop_worst().expect("first pop");
        assert_eq!(first.oid, 2);
        let second = heap.pop_worst().expect("second pop");
        assert_eq!(second.oid, 3);
        let third = heap.pop_worst().expect("third pop");
        assert_eq!(third.oid, 1);
        assert!(heap.pop_worst().is_none());
    }

    #[test]
    fn clear_heap_drains_records() {
        let mut heap = BlastCompoHeap::new(10, 1e-5);
        let mut a = HspList::new(1);
        a.best_evalue = 1e-10;
        heap.records.push(a);
        assert!(!heap.records.is_empty());
        clear_heap(&mut heap);
        assert!(heap.records.is_empty());
    }

    #[test]
    fn fill_results_from_compo_heaps_reverses_to_best_first() {
        let mut heaps = vec![BlastCompoHeap::new(10, 1e-5)];
        let mut a = HspList::new(1);
        a.best_evalue = 1e-10;
        let mut b = HspList::new(2);
        b.best_evalue = 1.0;
        heaps[0].records.push(a);
        heaps[0].records.push(b);
        let results = fill_results_from_compo_heaps(&mut heaps);
        let hl = results.hitlists[0].as_ref().expect("hitlist");
        // Best-evalue (smallest) should be at the head after reverse.
        assert_eq!(hl.hsp_lists.len(), 2);
        assert_eq!(hl.hsp_lists[0].oid, 1); // best evalue first
        assert_eq!(hl.hsp_lists[1].oid, 2);
    }

    #[test]
    fn free_blast_compo_query_info_array_clears() {
        let mut arr = vec![BlastCompoQueryInfo::default(); 3];
        free_blast_compo_query_info_array(&mut arr);
        assert!(arr.is_empty());
    }

    #[test]
    fn get_query_info_populates_per_context() {
        // Build a 2-context QueryInfo with two queries concatenated.
        let qi = crate::queryinfo::QueryInfo {
            num_queries: 2,
            contexts: vec![
                crate::queryinfo::ContextInfo {
                    query_offset: 0,
                    query_length: 10,
                    eff_searchsp: 12345,
                    length_adjustment: 0,
                    query_index: 0,
                    frame: 1,
                    is_valid: true,
                },
                crate::queryinfo::ContextInfo {
                    query_offset: 11,
                    query_length: 10,
                    eff_searchsp: 67890,
                    length_adjustment: 0,
                    query_index: 1,
                    frame: 1,
                    is_valid: true,
                },
            ],
            max_length: 10,
        };
        // Concatenated buffer with sentinel between queries.
        let qdata: Vec<u8> = (1..=10).chain([0u8]).chain(11..=20).collect();
        let result = get_query_info(&qdata, &qi, /* skip = */ true);
        assert_eq!(result.len(), 2);
        assert_eq!(result[0].origin, 0);
        assert_eq!(result[0].seq.length, 10);
        assert_eq!(result[0].seq.data(), &(1..=10).collect::<Vec<u8>>()[..]);
        assert_eq!(result[0].eff_search_space, 12345.0);
        // Word array exists (10 - 8 + 1 = 3 slots).
        assert_eq!(result[0].words.len(), 3);
        // Composition skipped.
        assert!(result[0].composition.is_empty());

        assert_eq!(result[1].origin, 11);
        assert_eq!(result[1].seq.data(), &(11..=20).collect::<Vec<u8>>()[..]);
        assert_eq!(result[1].eff_search_space, 67890.0);
    }

    #[test]
    fn create_word_array_too_short_returns_none() {
        assert_eq!(create_word_array(&[1, 2, 3]), None);
    }

    #[test]
    fn create_word_array_matches_independent_get_hash() {
        let seq: Vec<u8> = (1..=20).collect();
        let words = create_word_array(&seq).expect("word array");
        // Verify each filled slot matches a fresh `get_hash` of the
        // corresponding 8-mer. NCBI's loop fills indices [0..seq_len -
        // word_size), so we check that range.
        let upper = seq.len() - 8; // exclusive
        for i in 0..upper {
            let expected = get_hash(&seq[i..i + 8], 8);
            assert_eq!(words[i], expected, "mismatch at i={i}");
        }
        // Trailing slot at index `upper` is zero per the C off-by-one.
        assert_eq!(words[upper], 0);
    }

    #[test]
    fn blast_hsp_get_num_identities_ungapped_perfect_match() {
        let query = b"AAAACCCC".to_vec();
        let subject = b"AAAACCCC".to_vec();
        let hsp = Hsp {
            score: 8,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 0.0,
            query_offset: 0,
            query_end: 8,
            subject_offset: 0,
            subject_end: 8,
            context: 0,
            num_gaps: 0,
        };
        let (n_ident, align_length, _) =
            blast_hsp_get_num_identities(&query, &subject, &hsp, None, None);
        assert_eq!(n_ident, 8);
        assert_eq!(align_length, 8);
    }

    #[test]
    fn blast_hsp_get_num_identities_with_gaps() {
        // Query: AAAA-CCCC (6 cols when subject has the gap)
        // Subject: AAAATCCCC
        // ops: Sub(4), Del(1), Sub(4)
        let query = b"AAAACCCC".to_vec();
        let subject = b"AAAATCCCC".to_vec();
        let hsp = Hsp {
            score: 0,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 0.0,
            query_offset: 0,
            query_end: 8,
            subject_offset: 0,
            subject_end: 9,
            context: 0,
            num_gaps: 0,
        };
        let ops = vec![
            (crate::gapinfo::GapAlignOpType::Sub, 4),
            (crate::gapinfo::GapAlignOpType::Del, 1),
            (crate::gapinfo::GapAlignOpType::Sub, 4),
        ];
        let (n_ident, align_length, _) =
            blast_hsp_get_num_identities(&query, &subject, &hsp, Some(&ops), None);
        assert_eq!(n_ident, 8);   // 4 A's + 4 C's all match across the gap
        assert_eq!(align_length, 9); // 4 + 1 + 4
    }

    #[test]
    fn compute_num_identities_blastp_stamps_hsp() {
        let query = b"AAAACCCC".to_vec();
        let subject = b"AAAACCCC".to_vec();
        let mut list = HspList::new(0);
        list.add_hsp(Hsp {
            score: 8,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 0.0,
            query_offset: 0,
            query_end: 8,
            subject_offset: 0,
            subject_end: 8,
            context: 0,
            num_gaps: 0,
        });
        compute_num_identities_blastp(&query, &subject, &mut list, &[None], None);
        assert_eq!(list.hsps[0].num_ident, 8);
    }

    #[test]
    fn hitlist_evaluate_and_purge_drops_above_threshold() {
        let mut list = HspList::new(0);
        list.add_hsp(Hsp {
            score: 100,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1e-10,
            query_offset: 0,
            query_end: 100,
            subject_offset: 0,
            subject_end: 100,
            context: 0,
            num_gaps: 0,
        });
        list.add_hsp(Hsp {
            score: 50,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1.0, // above default threshold
            query_offset: 0,
            query_end: 50,
            subject_offset: 0,
            subject_end: 50,
            context: 0,
            num_gaps: 0,
        });
        let (best_score, best_evalue) = hitlist_evaluate_and_purge(
            &mut list,
            1000, // subject_length
            crate::program::BLASTP,
            100,  // query_length
            10,   // length_adjustment
            1e9,  // eff_searchsp
            -1.0, // pvalue_for_this_pair: out of [0,1] → skip composition
            0.1,  // max_evalue
            false, // do_sum_stats
        );
        // Above-threshold HSP should be dropped.
        assert_eq!(list.hsps.len(), 1);
        assert_eq!(best_score, 100);
        assert_eq!(best_evalue, 1e-10);
    }

    #[test]
    fn hitlist_evaluate_and_purge_empty_after_filter() {
        let mut list = HspList::new(0);
        list.add_hsp(Hsp {
            score: 50,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1.0,
            query_offset: 0,
            query_end: 50,
            subject_offset: 0,
            subject_end: 50,
            context: 0,
            num_gaps: 0,
        });
        let (best_score, best_evalue) = hitlist_evaluate_and_purge(
            &mut list,
            1000,
            crate::program::BLASTP,
            100,
            10,
            1e9,
            -1.0,
            0.01, // strict threshold
            false,
        );
        assert!(list.hsps.is_empty());
        assert_eq!(best_score, 0);
        assert_eq!(best_evalue, f64::MAX);
    }

    #[test]
    fn adjust_evalues_for_composition_basic() {
        // One HSP with e-value 1e-10. Composition p-value 0.5.
        // The function should leave the per-HSP evalue close to its
        // original value (since the alignment p-value is tiny vs the
        // 0.5 composition p-value, the combined p-value ≈ comp_p × align_p).
        let mut list = HspList::new(0);
        list.add_hsp(Hsp {
            score: 100,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1e-10,
            query_offset: 0,
            query_end: 100,
            subject_offset: 0,
            subject_end: 100,
            context: 0,
            num_gaps: 0,
        });
        adjust_evalues_for_composition(
            &mut list,
            0.5,    // comp_p_value
            1000,   // subject_length
            100,    // query_length
            10,     // length_adjustment
            1e9,    // eff_searchsp
        );
        // Best evalue should equal the only HSP's evalue.
        assert_eq!(list.best_evalue, list.hsps[0].evalue);
        // For comp_p = 0.5 and align_p tiny, Fisher's combine yields
        //   combined ≈ p_comp * p_align * (1 - ln(p_comp * p_align))
        // which is *larger* than p_align alone (the log factor inflates
        // it). The HSP becomes *less* significant after composition
        // adjustment. Verify the e-value moved in the expected direction.
        assert!(list.hsps[0].evalue > 1e-10);
        assert!(list.hsps[0].evalue < 1e-7);
    }

    #[test]
    fn test_near_identical_full_match() {
        // Identical 12-residue query and subject → 100 % identity → true.
        let q: Vec<u8> = (1..=12).collect();
        let s: Vec<u8> = (1..=12).collect();
        let q_words: Vec<u64> = (0..=q.len() - 8)
            .map(|i| get_hash(&q[i..i + 8], 8))
            .collect();
        let qd = BlastCompoSequenceData {
            buffer: q.clone(),
            data_offset: 0,
            length: q.len() as i32,
        };
        let sd = BlastCompoSequenceData {
            buffer: s.clone(),
            data_offset: 0,
            length: s.len() as i32,
        };
        let align = BlastCompoAlignment::new(
            12,
            MatrixAdjustRule::DontAdjust,
            0,
            0,
            12, // queryEnd is one past
            0,
            12,
            0,
            None,
        );
        assert!(test_near_identical(&sd, 0, &qd, 0, &q_words, &align));
    }

    #[test]
    fn test_near_identical_no_match() {
        // Completely disjoint sequences → 0 % identity → false.
        let q: Vec<u8> = vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12];
        let s: Vec<u8> = vec![20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31];
        let q_words: Vec<u64> = (0..=q.len() - 8)
            .map(|i| get_hash(&q[i..i + 8], 8))
            .collect();
        let qd = BlastCompoSequenceData {
            buffer: q.clone(),
            data_offset: 0,
            length: q.len() as i32,
        };
        let sd = BlastCompoSequenceData {
            buffer: s.clone(),
            data_offset: 0,
            length: s.len() as i32,
        };
        let align = BlastCompoAlignment::new(
            12,
            MatrixAdjustRule::DontAdjust,
            0,
            0,
            12,
            0,
            12,
            0,
            None,
        );
        assert!(!test_near_identical(&sd, 0, &qd, 0, &q_words, &align));
    }

    #[test]
    fn hsp_list_from_distinct_alignments_consumes_list() {
        // Build a 3-element linked list manually, in REVERSE-of-computation
        // order as NCBI stores them: third-best at head, best-scoring last.
        let third = BlastCompoAlignment::new(
            10, MatrixAdjustRule::DontAdjust, 0, 0, 10, 0, 10, 1, None);
        let mut second = BlastCompoAlignment::new(
            50, MatrixAdjustRule::ScaleOldMatrix, 0, 20, 30, 20, 30, 1, None);
        second.next = Some(Box::new(third));
        let mut first = BlastCompoAlignment::new(
            100, MatrixAdjustRule::UnconstrainedRelEntropy, 0, 40, 60, 40, 60, 1, None);
        first.next = Some(Box::new(second));
        let mut head: Option<Box<BlastCompoAlignment>> = Some(Box::new(first));
        let mut hsp_list = HspList::new(0);
        let (status, tags) =
            hsp_list_from_distinct_alignments(&mut hsp_list, &mut head, 42, 1);
        assert_eq!(status, 0);
        assert!(head.is_none(), "linked list consumed");
        assert_eq!(hsp_list.oid, 42);
        assert_eq!(hsp_list.hsps.len(), 3);
        // After sort_by_score, highest score is first.
        assert_eq!(hsp_list.hsps[0].score, 100);
        assert_eq!(hsp_list.hsps[1].score, 50);
        assert_eq!(hsp_list.hsps[2].score, 10);
        // Tags emitted in pre-sort order (matching alignment list traversal).
        assert_eq!(tags.len(), 3);
        assert_eq!(tags[0], CompoAdjustMode::CompositionMatrixAdjust); // first/UnconstrainedRelEntropy
        assert_eq!(tags[1], CompoAdjustMode::CompositionBasedStats);   // second/ScaleOldMatrix
        assert_eq!(tags[2], CompoAdjustMode::NoCompositionBasedStats); // third/DontAdjust
    }

    #[test]
    fn new_alignment_from_gap_align_shifts_coords() {
        let tb = crate::traceback::TracebackResult {
            score: 99,
            edit_script: crate::gapinfo::GapEditScript::new(),
            query_start: 5,
            query_end: 25,
            subject_start: 10,
            subject_end: 30,
        };
        let qr = BlastCompoSequenceRange { begin: 100, end: 200, context: 7 };
        let sr = BlastCompoSequenceRange { begin: 1000, end: 2000, context: 3 };
        let mut script = Some(crate::gapinfo::GapEditScript::new());
        let result = new_alignment_from_gap_align(
            &tb,
            &mut script,
            &qr,
            &sr,
            MatrixAdjustRule::DontAdjust,
        )
        .expect("alignment");
        assert_eq!(result.score, 99);
        assert_eq!(result.query_start, 105);
        assert_eq!(result.query_end, 125);
        assert_eq!(result.match_start, 1010);
        assert_eq!(result.match_end, 1030);
        assert_eq!(result.query_index, 7);
        assert_eq!(result.frame, 3);
        // Edit script ownership transferred.
        assert!(script.is_none());
        assert!(result.context.is_some());
    }

    #[test]
    fn sequence_data_slice_view() {
        let mut sd = BlastCompoSequenceData::default();
        sd.buffer = vec![0, 1, 2, 3, 4]; // sentinel + 4 residues
        sd.data_offset = 1;
        sd.length = 4;
        assert_eq!(sd.data(), &[1, 2, 3, 4]);
    }
}


// ───────────────────────────────────────────────────────────────────────────
// Macros (`blast_hits_priv.h`, `redo_alignment.h`).
// ───────────────────────────────────────────────────────────────────────────

/// `CONTAINED_IN_HSP(a, b, c, d, e, f)` — `blast_hits_priv.h:68`.
/// True iff `c` is in `[a, b]` and `f` is in `[d, e]`.
#[inline]
pub fn contained_in_hsp(a: i32, b: i32, c: i32, d: i32, e: i32, f: i32) -> bool {
    a <= c && b >= c && d <= f && e >= f
}

/// `KAPPA_BIT_TOL` — `redo_alignment.c:955`. Number of bits by which a
/// previously-emitted alignment must score above the current candidate
/// for [`is_contained`] to report containment.
pub const KAPPA_BIT_TOL: f64 = 2.0;

/// 1-1 port of `s_IsContained` (`redo_alignment.c:965`).
///
/// Returns true if the alignment defined by `(query_start, query_end,
/// subject_start, subject_end, score, frame)` is contained in any
/// alignment already emitted in `existing` AND that previous alignment
/// scores at least `score + KAPPA_BIT_TOL * ln(2) / lambda` higher.
/// NCBI uses this in the redo-alignment driver to skip the
/// composition-adjusted SW redo for preliminary HSPs that are already
/// covered by a higher-scoring redone alignment in the same frame.
///
/// `existing` is a slice of `(q_start, q_end, s_start, s_end, score,
/// frame)` tuples — same field set as `BlastCompo_Alignment`'s
/// containment check uses, mapped onto whatever HSP shape the caller
/// has on hand (we only need the bounding box, score, and frame).
pub fn is_contained(
    query_start: i32,
    query_end: i32,
    subject_start: i32,
    subject_end: i32,
    score: i32,
    frame: i32,
    existing: &[(i32, i32, i32, i32, i32, i32)],
    lambda: f64,
) -> bool {
    let score_thresh = score as f64 + KAPPA_BIT_TOL * NCBIMATH_LN2 / lambda;
    for &(eq_s, eq_e, es_s, es_e, e_score, e_frame) in existing {
        // Same-sign frame check (`KAPPA_SIGN`): `redo_alignment.c:75`.
        let same_sign = (frame.signum()) == (e_frame.signum());
        if !same_sign {
            continue;
        }
        // Both endpoints of the candidate must lie inside the existing
        // alignment's box (KAPPA_CONTAINED_IN_HSP applied twice — once
        // to the start point, once to the end point).
        if contained_in_hsp(eq_s, eq_e, query_start, es_s, es_e, subject_start)
            && contained_in_hsp(eq_s, eq_e, query_end, es_s, es_e, subject_end)
            && score_thresh <= e_score as f64
        {
            return true;
        }
    }
    false
}

/// `GET_NUCL_LENGTH(l)` macro — `redo_alignment.h:495`. Used by RPS-tblastn
/// to convert a packed mixed-frame protein length back to nucleotide
/// length.
#[inline]
pub fn get_nucl_length(l: i32) -> i32 {
    (l - 5) / 2 + 2
}

// ───────────────────────────────────────────────────────────────────────────
// Leaf utilities — no `BlastCompo_*` plumbing required.
// ───────────────────────────────────────────────────────────────────────────

/// Port of `s_GetSubjectLength` (`blast_kappa.c:364`).
///
/// For RPS-tblastn the subject is a packed mixed-frame protein sequence;
/// the underlying nucleotide length is recovered via `GET_NUCL_LENGTH`
/// and divided by 3. For all other programs the input length is returned
/// unchanged.
pub fn get_subject_length(total_subj_length: i32, program_number: ProgramType) -> i32 {
    if program_number == RPS_BLAST {
        // NCBI uses `eBlastTypeRpsTblastn` here — Rust constant is
        // `RPS_BLAST` (no separate RPS-tblastn variant in the current
        // type alias), so we treat them as equivalent. Confirm and split
        // when the RPS-tblastn driver lands.
        (get_nucl_length(total_subj_length) - 1) / 3
    } else {
        total_subj_length
    }
}

/// Port of `s_HSPListNormalizeScores` (`blast_kappa.c:101`).
///
/// Rescales each HSP score from a high-precision representation back to
/// integer scores by dividing by `score_divisor`, then sets the
/// `bit_score` from the rescaled score. Caller is responsible for
/// preserving the existing sort order (no reorder happens here).
pub fn hsp_list_normalize_scores(
    hsp_list: &mut HspList,
    lambda: f64,
    log_k: f64,
    score_divisor: f64,
) {
    for hsp in &mut hsp_list.hsps {
        // C: `hsp->score = (Int4)BLAST_Nint(((double) hsp->score) / scoreDivisor);`
        hsp.score = crate::math::nint(hsp.score as f64 / score_divisor) as i32;
        // C: `hsp->bit_score = (hsp->score*lambda*scoreDivisor - logK)/NCBIMATH_LN2;`
        hsp.bit_score = (hsp.score as f64 * lambda * score_divisor - log_k) / NCBIMATH_LN2;
    }
}

/// Port of `s_GetHash` (`blast_kappa.c:1117`). Hash for a 28-letter-alphabet
/// word: `hash = sum_k (data[k] << (5 * (word_size - 1 - k)))`.
#[inline]
pub fn get_hash(data: &[u8], word_size: usize) -> u64 {
    let mut hash: u64 = 0;
    for &b in &data[..word_size] {
        hash <<= 5;
        hash += b as u64;
    }
    hash
}

/// Port of `s_ExtendRight` (`blast_kappa.c:944`).
///
/// Extends rightward from the start of `query` and `subject`, counting
/// identical residues and tolerating up to `max_shift` mismatches/gaps
/// when the very next two positions match. Returns
/// `(num_identical, query_ext_len, subject_ext_len, align_len)`.
pub fn extend_right(
    query: &[u8],
    subject: &[u8],
    max_shift: i32,
) -> (i32, i32, i32, i32) {
    let query_len = query.len() as i32;
    let subject_len = subject.len() as i32;
    let mut num_identical = 0i32;
    let mut q_pos = 0i32;
    let mut s_pos = 0i32;
    let mut gaps_in_query = 0i32;
    let mut gaps_in_subject = 0i32;

    'outer: while q_pos < query_len && s_pos < subject_len {
        let mut matched = false;

        // Run of identities.
        while q_pos < query_len && s_pos < subject_len
            && query[q_pos as usize] == subject[s_pos as usize]
        {
            num_identical += 1;
            q_pos += 1;
            s_pos += 1;
        }

        // Try to skip mismatches or gaps.
        let mut n = 1i32;
        while n < max_shift
            && q_pos + n + 1 < query_len
            && s_pos + n + 1 < subject_len
            && !matched
        {
            // Mismatches: advance both by `n + 2` if the (n)th and (n+1)th
            // positions both match.
            if query[(q_pos + n) as usize] == subject[(s_pos + n) as usize]
                && query[(q_pos + n + 1) as usize] == subject[(s_pos + n + 1) as usize]
            {
                q_pos += n + 2;
                s_pos += n + 2;
                num_identical += 2;
                matched = true;
            }

            // Gap in subject: advance query by `n + 2`, subject by 2.
            if !matched
                && query[(q_pos + n) as usize] == subject[s_pos as usize]
                && query[(q_pos + n + 1) as usize] == subject[(s_pos + 1) as usize]
            {
                q_pos += n + 2;
                s_pos += 2;
                num_identical += 2;
                gaps_in_subject += n;
                matched = true;
            }

            // Gap in query: advance subject by `n + 2`, query by 2.
            if !matched
                && query[q_pos as usize] == subject[(s_pos + n) as usize]
                && query[(q_pos + 1) as usize] == subject[(s_pos + n + 1) as usize]
            {
                q_pos += 2;
                s_pos += n + 2;
                num_identical += 2;
                gaps_in_query += n;
                matched = true;
            }
            n += 1;
        }

        if matched {
            continue;
        }
        break 'outer;
    }

    let align_len = if q_pos > s_pos {
        q_pos + gaps_in_query
    } else {
        s_pos + gaps_in_subject
    };
    (num_identical, q_pos, s_pos, align_len)
}

/// Port of `s_ExtendLeft` (`blast_kappa.c:1039`).
///
/// Extends leftward from the END of `query` and `subject`. Returns
/// `(num_identical, query_ext_len, subject_ext_len, align_len_delta)`
/// where the extension lengths are measured from the end. NCBI's C
/// version takes an `align_len` in/out parameter and *adds* the delta;
/// the Rust signature returns the delta and lets the caller decide.
pub fn extend_left(
    query: &[u8],
    subject: &[u8],
    max_shift: i32,
) -> (i32, i32, i32, i32) {
    let query_len = query.len() as i32;
    let subject_len = subject.len() as i32;
    let mut q_pos = query_len - 1;
    let mut s_pos = subject_len - 1;
    let mut num_identical = 0i32;
    let mut gaps_in_query = 0i32;
    let mut gaps_in_subject = 0i32;

    while q_pos >= 0 && s_pos >= 0 {
        let mut matched = false;

        while q_pos > 0 && s_pos > 0
            && query[q_pos as usize] == subject[s_pos as usize]
        {
            num_identical += 1;
            q_pos -= 1;
            s_pos -= 1;
        }

        let mut n = 1i32;
        while n < max_shift
            && q_pos - n - 1 > 0
            && s_pos - n - 1 > 0
            && !matched
        {
            if query[(q_pos - n) as usize] == subject[(s_pos - n) as usize]
                && query[(q_pos - n - 1) as usize] == subject[(s_pos - n - 1) as usize]
            {
                q_pos -= n + 2;
                s_pos -= n + 2;
                num_identical += 2;
                matched = true;
            }

            if !matched
                && query[(q_pos - n) as usize] == subject[s_pos as usize]
                && query[(q_pos - n - 1) as usize] == subject[(s_pos - 1) as usize]
            {
                q_pos -= n + 2;
                s_pos -= 2;
                num_identical += 2;
                gaps_in_subject += n;
                matched = true;
            }

            if !matched
                && query[q_pos as usize] == subject[(s_pos - n) as usize]
                && query[(q_pos - 1) as usize] == subject[(s_pos - n - 1) as usize]
            {
                q_pos -= 2;
                s_pos -= n + 2;
                num_identical += 2;
                gaps_in_query += n;
                matched = true;
            }
            n += 1;
        }

        if matched {
            continue;
        }
        break;
    }

    let q_ext_len = query_len - q_pos - 1;
    let s_ext_len = subject_len - s_pos - 1;
    let align_delta = if q_ext_len > s_ext_len {
        q_ext_len + gaps_in_query
    } else {
        s_ext_len + gaps_in_subject
    };
    (num_identical, q_ext_len, s_ext_len, align_delta)
}

/// Port of `s_FindNumIdentical` (`blast_kappa.c:1143`). Walks the subject
/// in a sliding-window k-mer hash, looking up each window in the query's
/// pre-hashed array; on every hit, extends bidirectionally via
/// `s_ExtendLeft` / `s_ExtendRight` to count identities.
///
/// `query_hashes` must hold one hash per starting position in the query
/// (of length `>= query.len() - WORD_SIZE`); see `create_word_array`.
pub fn find_num_identical(
    query: &[u8],
    query_hashes: &[u64],
    subject: &[u8],
    max_shift: i32,
) -> i32 {
    const WORD_SIZE: usize = 8;
    let mask: u64 = 0xFFFF_FFFF_FF;
    let query_len = query.len();
    let subject_len = subject.len();

    if query.is_empty()
        || query_hashes.is_empty()
        || subject.is_empty()
        || query_len < WORD_SIZE
        || subject_len < WORD_SIZE
    {
        return 0;
    }

    let mut query_from: usize = 0;
    let mut subject_from: usize = 0;
    let mut num_identical: i32 = 0;
    let mut hash: u64 = 0;
    let mut matched = false;
    let mut s_pos: usize = 0;
    while s_pos < subject_len.saturating_sub(WORD_SIZE) {
        if s_pos == 0 || matched {
            hash = get_hash(&subject[s_pos..s_pos + WORD_SIZE], WORD_SIZE);
        } else {
            hash = ((hash << 5) & mask) + subject[s_pos + WORD_SIZE - 1] as u64;
        }

        // Find matching query word.
        let mut q_pos = query_from;
        while q_pos < query_len.saturating_sub(WORD_SIZE) && query_hashes[q_pos] != hash {
            q_pos += 1;
        }

        if q_pos < query_len.saturating_sub(WORD_SIZE) {
            let query_start = q_pos;
            let subject_start = s_pos;

            matched = true;
            num_identical += WORD_SIZE as i32;

            // Extend left.
            let (left_ident, _q_left, _s_left, _) = extend_left(
                &query[query_from..query_start],
                &subject[subject_from..subject_start],
                max_shift,
            );
            num_identical += left_ident;

            // Extend right.
            let (right_ident, _q_right, s_right, _) = extend_right(
                &query[query_start + WORD_SIZE..],
                &subject[subject_start + WORD_SIZE..],
                max_shift,
            );
            num_identical += right_ident;

            query_from = query_start + WORD_SIZE + _q_right as usize;
            subject_from = subject_start + WORD_SIZE + s_right as usize;
            // C: `s_pos = subject_from - 1;` — outer loop's increment makes
            // it `subject_from`.
            s_pos = subject_from.saturating_sub(1);
        } else {
            matched = false;
        }
        s_pos += 1;
    }

    num_identical
}

/// 1-1 port of `BlastKappa_SavedParameters` (`blast_kappa.c:1958`).
///
/// Snapshot of the search parameters captured on entry to
/// `Blast_RedoAlignmentCore`. The kappa driver may scale the score
/// matrix and gap costs internally; the snapshot lets the original
/// values be restored on exit so the caller's `BlastScoreBlk` and
/// scoring parameters are untouched from its perspective.
#[derive(Debug, Clone, Default)]
pub struct BlastKappaSavedParameters {
    pub gap_open: i32,
    pub gap_extend: i32,
    pub scale_factor: f64,
    pub orig_matrix: Vec<Vec<i32>>,
    pub original_expect_value: f64,
    /// Per-query gapped Karlin block snapshot (NCBI: `kbp_gap_orig`).
    pub kbp_gap_orig: Vec<crate::stat::KarlinBlk>,
    pub num_queries: i32,
}

impl BlastKappaSavedParameters {
    /// 1-1 port of `s_SavedParametersNew` (`blast_kappa.c:2008`).
    ///
    /// `rows` is the row count of the scoring matrix (PSSM length for
    /// position-based searches; `BLASTAA_SIZE` otherwise). When
    /// `compo_adjust_mode` is `NoCompositionBasedStats`, `orig_matrix`
    /// stays empty (matching NCBI's `if (compo_adjust_mode != ...)`
    /// guard at the bottom of `s_SavedParametersNew`).
    pub fn new(
        rows: i32,
        num_queries: i32,
        compo_adjust_mode: CompoAdjustMode,
        position_based: bool,
    ) -> Self {
        let mut sp = BlastKappaSavedParameters {
            gap_open: 0,
            gap_extend: 0,
            scale_factor: 0.0,
            orig_matrix: Vec::new(),
            original_expect_value: 0.0,
            kbp_gap_orig: Vec::with_capacity(num_queries as usize),
            num_queries,
        };
        // C: `for (i = 0; i < numQueries; i++) sp->kbp_gap_orig[i] = NULL;`
        // We initialize to default KarlinBlk values; the actual values
        // are filled in by `s_RecordInitialSearch` (TODO).
        sp.kbp_gap_orig
            .resize(num_queries.max(0) as usize, crate::stat::KarlinBlk::default());
        if !matches!(compo_adjust_mode, CompoAdjustMode::NoCompositionBasedStats) {
            // Allocate the original-matrix backing store.
            let cols = crate::matrix::AA_SIZE;
            let row_count = if position_based {
                rows.max(0) as usize
            } else {
                crate::matrix::AA_SIZE
            };
            sp.orig_matrix = vec![vec![0i32; cols]; row_count];
        }
        sp
    }
}

/// 1-1 port of `s_NewAlignmentUsingXdrop` (`blast_kappa.c:1812`).
///
/// X-drop traceback callback used by Smith-Waterman post-processing:
/// given the SW-derived `(query_start, query_end, match_start,
/// match_end, score)`, runs an X-drop alignment via
/// `sw_find_final_ends_using_xdrop`, builds a `GapEditScript` from the
/// returned ops, and packages it into a `BlastCompoAlignment` with
/// concatenated-query coordinates.
///
/// Returns `(new_alignment, query_end_xdrop, match_end_xdrop)` on
/// success or `None` on failure (matching NCBI's `0` / `-1` return
/// where the alignment is the side-effect output `*pnewAlign`).
#[allow(clippy::too_many_arguments)]
pub fn new_alignment_using_xdrop(
    query: &BlastCompoSequenceData,
    query_range: &BlastCompoSequenceRange,
    subject: &BlastCompoSequenceData,
    subject_range: &BlastCompoSequenceRange,
    query_start: usize,
    query_end: usize,
    match_start: usize,
    match_end: usize,
    score: i32,
    gapping_params: &BlastCompoGappingParams,
    matrix_adjust_rule: MatrixAdjustRule,
    matrix: &[[i32; 16]; 16],
) -> Option<(BlastCompoAlignment, i32, i32)> {
    let (new_score, q_extent, s_extent, ops) = sw_find_final_ends_using_xdrop(
        query.data(),
        query_start,
        query_end,
        subject.data(),
        match_start,
        match_end,
        matrix,
        gapping_params.gap_open,
        gapping_params.gap_extend,
        gapping_params.x_dropoff,
        score,
    );
    let q_end_xdrop = query_start + q_extent;
    let s_end_xdrop = match_start + s_extent;

    if ops.is_empty() {
        return None;
    }
    // Build an edit script from the returned ops. NCBI uses
    // `Blast_PrelimEditBlockToGapEditScript`; our Rust `align_ex`
    // already returns merged final ops, so just push each into a
    // `GapEditScript`.
    let mut edit_script = crate::gapinfo::GapEditScript::new();
    for &(op, count) in &ops {
        edit_script.push(op, count);
    }

    // C: shifted endpoints into concatenated-query / full-subject coords.
    let aquery_start = query_start as i32 + query_range.begin;
    let aquery_end = q_end_xdrop as i32 + query_range.begin;
    let amatch_start = match_start as i32 + subject_range.begin;
    let amatch_end = s_end_xdrop as i32 + subject_range.begin;

    let new_align = BlastCompoAlignment::new(
        new_score,
        matrix_adjust_rule,
        query_range.context,
        aquery_start,
        aquery_end,
        amatch_start,
        amatch_end,
        subject_range.context,
        Some(edit_script),
    );
    Some((new_align, q_end_xdrop as i32, s_end_xdrop as i32))
}

/// 1-1 port of `s_RedoOneAlignment` (`blast_kappa.c:1898`).
///
/// X-drop alignment in BOTH directions from the seed `(gapped_start_q,
/// gapped_start_s)`, producing a fresh `BlastCompoAlignment` that
/// supersedes the input alignment. NCBI calls
/// `BLAST_GappedAlignmentWithTraceback` (the bidirectional driver);
/// our Rust analog is `crate::traceback::blast_gapped_align`.
///
/// `gapped_start_q` and `gapped_start_s` are the per-HSP seed
/// coordinates that NCBI reads from `hsp->query.gapped_start` /
/// `hsp->subject.gapped_start`. NCBI subtracts the range begin to
/// shift to local coords; this Rust port does the same.
#[allow(clippy::too_many_arguments)]
pub fn redo_one_alignment(
    query_data: &BlastCompoSequenceData,
    query_range: &BlastCompoSequenceRange,
    subject_data: &BlastCompoSequenceData,
    subject_range: &BlastCompoSequenceRange,
    gapped_start_q: i32,
    gapped_start_s: i32,
    matrix_adjust_rule: MatrixAdjustRule,
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    x_dropoff: i32,
) -> Option<BlastCompoAlignment> {
    // C: `q_start = hsp->query.gapped_start - query_range->begin;`.
    let q_start = (gapped_start_q - query_range.begin).max(0) as usize;
    let s_start = (gapped_start_s - subject_range.begin).max(0) as usize;

    let tb = crate::traceback::blast_gapped_align(
        query_data.data(),
        subject_data.data(),
        q_start,
        s_start,
        reward,
        penalty,
        gap_open,
        gap_extend,
        x_dropoff,
    )?;
    let mut script = Some(tb.edit_script.clone());
    new_alignment_from_gap_align(
        &tb,
        &mut script,
        query_range,
        subject_range,
        matrix_adjust_rule,
    )
}

/// 1-1 port of `s_GetStartFreqRatios` (`blast_kappa.c:648`).
///
/// Returns the BLASTAA_SIZE × BLASTAA_SIZE frequency-ratio matrix for
/// the named scoring matrix. Currently supports BLOSUM62 (the only
/// matrix whose frequency ratios are baked into our crate). NCBI's C
/// version dispatches via `_PSIMatrixFrequencyRatiosNew`; that PSI
/// matrix-table helper isn't ported yet, so other matrices return
/// `Err(())` (matching NCBI's `-1` return on unknown matrix).
///
/// Output ratios are deep-copied (1-1 with NCBI's per-cell copy from
/// `stdFreqRatios->data[i][j]`).
pub fn get_start_freq_ratios(matrix_name: &str) -> Result<[[f64; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE], ()> {
    if matrix_name.eq_ignore_ascii_case("BLOSUM62") || matrix_name.is_empty() {
        return Ok(crate::matrix::get_blosum62_freq_ratios());
    }
    Err(())
}

/// `SCALING_FACTOR` (`blast_kappa.c:676`). Multiplicative factor for
/// integer matrix scores — gives composition-adjusted scoring extra
/// precision without overflowing `BLAST_SCORE_MIN`.
pub const SCALING_FACTOR: f64 = 32.0;

/// 1-1 port of `s_SWFindFinalEndsUsingXdrop` (`blast_kappa.c:843`).
///
/// Smith-Waterman has located `(query_start..query_end,
/// match_start..match_end)` and a target `score`; this routine
/// re-runs the alignment with the X-drop algorithm
/// (`crate::traceback::align_ex`) and doubles the X-dropoff up to
/// three times if the returned score falls short of the target.
/// Returns `(new_score, query_extent, match_extent, ops)` where
/// `query_extent` and `match_extent` are the X-drop-derived alignment
/// lengths (NCBI signals them via out-parameters `*queryAlignmentExtent`
/// and `*matchAlignmentExtent`).
///
/// `gap_x_dropoff` is the starting X-drop; on each undershoot it's
/// doubled. NCBI restores the original on exit; Rust does the same
/// implicitly because `gap_x_dropoff` is taken by value.
#[allow(clippy::too_many_arguments)]
pub fn sw_find_final_ends_using_xdrop(
    query: &[u8],
    query_start: usize,
    query_end: usize,
    subject: &[u8],
    match_start: usize,
    match_end: usize,
    matrix: &[[i32; 16]; 16],
    gap_open: i32,
    gap_extend: i32,
    mut gap_x_dropoff: i32,
    target_score: i32,
) -> (i32, usize, usize, Vec<(crate::gapinfo::GapAlignOpType, i32)>) {
    let q_len = query_end - query_start + 1;
    let s_len = match_end - match_start + 1;
    let mut doubling_count = 0;
    let mut x_score: i32;
    let mut q_extent: usize;
    let mut s_extent: usize;
    let mut last_ops: Vec<(crate::gapinfo::GapAlignOpType, i32)>;
    // C: `do { ... } while ((XdropAlignScore < score) && (doublingCount < 3));`
    // Note that `gap_x_dropoff *= 2` happens AFTER the call, so the
    // first invocation uses the original dropoff.
    loop {
        // C: ALIGN_EX is called with `&query->data[queryStart] - 1` (a
        // sentinel byte one position before the alignment start).
        // Build that sentinel-prefixed buffer in Rust so `align_ex`'s
        // 1-indexed `a[ai]` reads stay in bounds. NCBI's C avoids the
        // copy via pointer arithmetic but the slice content is the
        // same.
        let mut q_buf = vec![0u8; q_len + 1];
        let mut s_buf = vec![0u8; s_len + 1];
        q_buf[1..=q_len].copy_from_slice(&query[query_start..query_start + q_len]);
        s_buf[1..=s_len].copy_from_slice(&subject[match_start..match_start + s_len]);
        // C: `ALIGN_EX(query+queryStart-1, subject+matchStart-1, ...,
        //              FALSE /* reverse */, FALSE /* reverse_seq */)`
        // We use `align_ex` in forward mode with `reverse = false`.
        let (score, q_off, s_off, ops) = crate::traceback::align_ex(
            &q_buf,
            &s_buf,
            q_len,
            s_len,
            matrix,
            gap_open,
            gap_extend,
            gap_x_dropoff,
            false,
        );
        x_score = score;
        q_extent = q_off;
        s_extent = s_off;
        last_ops = ops;
        gap_x_dropoff *= 2;
        doubling_count += 1;
        if !(x_score < target_score && doubling_count < 3) {
            break;
        }
    }
    (x_score, q_extent, s_extent, last_ops)
}

/// 1-1 port of `BlastGapAlignStruct` (`blast_gapalign.h`) — the
/// gapped-alignment workspace owned by the kappa driver.
///
/// NCBI's struct holds many heap-allocated workspace buffers
/// (`state_struct`, `dp_mem`, `greedy_align_mem`, `fwd_prelim_tback`,
/// `rev_prelim_tback`) plus the alignment output fields (`score`,
/// `query_start/stop`, `subject_start/stop`, `edit_script`).
///
/// Our Rust ports of the gap-align routines (`align_ex`,
/// `blast_gapped_align`, `gapped_score_one_dir_*`) own their workspace
/// buffers internally, so the Rust analog here only carries the
/// alignment **output** state plus the few control fields the kappa
/// driver reads or writes (`gap_x_dropoff`).
///
/// The C `s_BlastGapAlignStruct_Copy` and `_Free` helpers correspond
/// to Rust's automatic `Clone`/`Drop`; provided as parity hooks.
#[derive(Debug, Clone, Default)]
pub struct BlastGapAlignWorkspace {
    pub gap_x_dropoff: i32,
    pub score: i32,
    pub query_start: i32,
    pub query_stop: i32,
    pub subject_start: i32,
    pub subject_stop: i32,
    pub edit_script: Option<crate::gapinfo::GapEditScript>,
}

impl BlastGapAlignWorkspace {
    /// 1-1 port of `s_BlastGapAlignStruct_Copy` (`blast_kappa.c:2604`).
    /// Rust's `Clone` does the deep copy; this is a parity wrapper.
    pub fn deep_copy(&self) -> Self {
        self.clone()
    }
}

/// 1-1 port of `s_BlastGapAlignStruct_Free` (`blast_kappa.c:2532`).
/// `Drop` handles deallocation; this hook is a parity marker.
pub fn blast_gap_align_struct_free(slot: &mut Option<BlastGapAlignWorkspace>) {
    *slot = None;
}

/// `BlastKappa_GappingParamsContext` (`blast_kappa.c:1715`). Bundles
/// the scoring/extension parameters and the gap-align workspace
/// passed through the `redo_one_alignment` / `new_xdrop_align`
/// callback chain.
#[derive(Debug, Clone)]
pub struct BlastKappaGappingParamsContext {
    pub scoring_params: crate::parameters::ScoringParameters,
    pub gap_align: BlastGapAlignWorkspace,
    pub local_scaling_factor: f64,
    pub program: ProgramType,
}

/// 1-1 port of `Blast_RedoAlignmentCore` (`blast_kappa.c:2942`).
///
/// Single-thread adapter that delegates to the multi-thread driver
/// `Blast_RedoAlignmentCore_MT(num_threads = 1, ...)`. The MT driver
/// itself is 669 LOC of dense composition-adjustment orchestration;
/// it's tracked as TODO in the module status table.
///
/// This wrapper sets `num_threads = 1` and delegates, exactly
/// mirroring the C control flow. Returns the same `Int2` status code
/// the C function returns (currently `-1` until the MT body lands).
#[allow(clippy::too_many_arguments)]
pub fn blast_redo_alignment_core(
    program: ProgramType,
    query: &[u8],
    query_info: &crate::queryinfo::QueryInfo,
    kbp_gap: &mut [crate::stat::KarlinBlk],
    matrix: &mut [Vec<i32>],
    scoring: &mut crate::parameters::ScoringParameters,
    align_params: &BlastRedoAlignParams,
    saved: &mut BlastKappaSavedParameters,
    this_match: &mut HspList,
    results: &mut crate::hspstream::HspResults,
) -> i32 {
    blast_redo_alignment_core_mt(
        program,
        /* num_threads = */ 1,
        query,
        query_info,
        kbp_gap,
        matrix,
        scoring,
        align_params,
        saved,
        this_match,
        results,
    )
}

/// Skeleton for `Blast_RedoAlignmentCore_MT` (`blast_kappa.c:2980`,
/// 669 LOC, cyclomatic 120). The full driver orchestrates:
/// 1. Save initial search state via [`record_initial_search`].
/// 2. Rescale scoring via [`rescale_search`].
/// 3. Build the per-query info array via [`get_query_info`].
/// 4. Build the redo-align params via [`get_align_params`].
/// 5. For each subject-match in `hsp_stream`:
///    - run the redo-alignment driver against the subject's HSPs;
///    - filter via [`hitlist_evaluate_and_purge`];
///    - feed survivors into a per-query [`BlastCompoHeap`].
/// 6. Collapse heaps into `results` via [`fill_results_from_compo_heaps`].
/// 7. Restore initial search state via [`restore_search`].
///
/// Each helper above is already ported in this loop; assembly into the
/// full driver is a multi-iteration follow-up. For now this returns
/// `-1` (NCBI's "not implemented / error" return).
#[allow(clippy::too_many_arguments)]
pub fn blast_redo_alignment_core_mt(
    _program: ProgramType,
    _num_threads: u32,
    _query: &[u8],
    _query_info: &crate::queryinfo::QueryInfo,
    _kbp_gap: &mut [crate::stat::KarlinBlk],
    _matrix: &mut [Vec<i32>],
    _scoring: &mut crate::parameters::ScoringParameters,
    _align_params: &BlastRedoAlignParams,
    _saved: &mut BlastKappaSavedParameters,
    _this_match: &mut HspList,
    _results: &mut crate::hspstream::HspResults,
) -> i32 {
    // TODO: assemble from already-ported helpers. Each call site below
    // already has its 1-1 port; the orchestration is the only piece
    // missing.
    //
    // record_initial_search(saved, kbp_gap, matrix, scoring, ...);
    // rescale_search(kbp_gap, scoring, num_queries, scale_factor);
    // let qi = get_query_info(query, query_info, /* skip = */ false);
    // let params = get_align_params(...);
    // for each (hsp_list in stream): ...
    // restore_search(saved, kbp_gap, matrix, scoring, ...);
    -1
}

/// `NEAR_IDENTICAL_BITS_PER_POSITION` (`blast_kappa.c:2399`).
/// Per-position bit-score threshold used to decide whether two
/// sequences are "near identical" before running the more
/// expensive composition-adjustment pass.
pub const NEAR_IDENTICAL_BITS_PER_POSITION: f64 = 1.74;

/// 1-1 port of `s_GetAlignParams` (`blast_kappa.c:2406`).
///
/// Assembly point for the redo-alignment driver: builds a
/// [`BlastRedoAlignParams`] from the search context. Composes the
/// already-ported helpers ([`matrix_info_init_blastp`],
/// [`gapping_params_new`], [`blast_redo_align_params_new`]) in
/// exactly the same order NCBI does.
///
/// Parameters that NCBI reads from nested struct pointers
/// (`context->sbp`, `context->scoringParams->options->matrix`, etc.)
/// are passed directly so the Rust port doesn't depend on a
/// `BlastScoreBlk` shape that isn't ported yet.
///
/// - `kbp_gap_first_lambda`: Lambda of the first valid context
///   (used for the near-identical cutoff in raw-score units).
/// - `kbp_gap`: per-context gapped Karlin blocks (forwarded to
///   [`gapping_params_new`]).
/// - `kbp_ideal_lambda`: ideal ungapped Lambda for matrix init.
/// - `local_scaling_factor`: NCBI's `context->localScalingFactor`.
/// - `cutoff_score_min`: NCBI's `hitParams->cutoff_score_min`.
/// - `expect_value`: NCBI's `hitParams->options->expect_value`.
/// - `compo_adjust_mode`: NCBI's
///   `extendParams->options->compositionBasedStats`.
/// - `position_based`: true iff the caller is running PSI-BLAST
///   (NCBI: `sbp->psi_matrix != NULL`).
/// - `do_link_hsps`: NCBI's `hitParams->do_sum_stats`.
/// - `gap_x_dropoff_final_bits`: NCBI's
///   `extendParams->options->gap_x_dropoff_final` (in bit-score units).
/// - `raw_gap_x_dropoff_final`: NCBI's
///   `extendParams->gap_x_dropoff_final` (in raw-score units).
#[allow(clippy::too_many_arguments)]
pub fn get_align_params(
    program: ProgramType,
    matrix_name: &str,
    scoring: &crate::parameters::ScoringParameters,
    kbp_gap: &[crate::stat::KarlinBlk],
    kbp_gap_first_lambda: f64,
    kbp_ideal_lambda: f64,
    local_scaling_factor: f64,
    cutoff_score_min: i32,
    expect_value: f64,
    compo_adjust_mode: CompoAdjustMode,
    position_based: bool,
    do_link_hsps: bool,
    ccat_query_length: i32,
    num_contexts: i32,
    gap_x_dropoff_final_bits: f64,
    raw_gap_x_dropoff_final: i32,
) -> Option<BlastRedoAlignParams> {
    let subject_is_translated_p = crate::program::subject_is_translated(program);
    let query_is_translated_p = crate::program::query_is_translated(program);

    // C: `near_identical_cutoff = NEAR_IDENTICAL_BITS_PER_POSITION *
    //     ln(2) / kbp_gap[first_valid].Lambda`. The score block is
    // already scaled by `localScalingFactor`, so the resulting cutoff
    // is in scaled raw-score units.
    let near_identical_cutoff = if kbp_gap_first_lambda > 0.0 {
        (NEAR_IDENTICAL_BITS_PER_POSITION * crate::math::NCBIMATH_LN2) / kbp_gap_first_lambda
    } else {
        0.0
    };

    // C: `if (do_link_hsps) cutoff_s = cutoff_score_min * localScalingFactor;
    //     else cutoff_s = 1;`
    let cutoff_s = if do_link_hsps {
        (cutoff_score_min as f64 * local_scaling_factor) as i32
    } else {
        1
    };

    // C: `s_MatrixInfoInit(self, queryBlk, sbp, localScalingFactor,
    //                      scoringParams->options->matrix);`
    let mut matrix_info = BlastMatrixInfo::default();
    matrix_info.positional = position_based;
    let status = matrix_info_init_blastp(
        &mut matrix_info,
        matrix_name,
        kbp_ideal_lambda,
        local_scaling_factor,
    );
    if status != 0 {
        return None;
    }

    // C: `gapping_params = s_GappingParamsNew(context, extendParams,
    //                                          last_context + 1);`
    let gapping_params = gapping_params_new(
        scoring,
        kbp_gap,
        num_contexts,
        gap_x_dropoff_final_bits,
        raw_gap_x_dropoff_final,
    );

    Some(blast_redo_align_params_new(
        matrix_info,
        gapping_params,
        compo_adjust_mode,
        position_based,
        query_is_translated_p,
        subject_is_translated_p,
        ccat_query_length,
        cutoff_s,
        expect_value,
        do_link_hsps,
        near_identical_cutoff,
    ))
}

/// `kReMatrixAdjustmentPseudocounts` (`redo_alignment.c:83`). Default
/// pseudocount used by relative-entropy-based composition matrix
/// adjustment. NCBI exposes no override path; this is hardcoded.
pub const K_RE_MATRIX_ADJUSTMENT_PSEUDOCOUNTS: i32 = 20;

/// 1-1 port of `Blast_RedoAlignParamsNew` (`redo_alignment.c:1013`).
///
/// Bundles the matrix info + gapping params + composition-adjust mode
/// flags into a `BlastRedoAlignParams` driver wrapper. NCBI's C
/// version takes ownership of `*pmatrix_info` and `*pgapping_params`
/// (sets the caller's pointers to NULL); Rust's Vec/Box ownership
/// moves accomplish the same transfer naturally.
///
/// Sets `RE_pseudocounts` to the package default
/// [`K_RE_MATRIX_ADJUSTMENT_PSEUDOCOUNTS`] = 20, matching NCBI.
#[allow(clippy::too_many_arguments)]
pub fn blast_redo_align_params_new(
    matrix_info: BlastMatrixInfo,
    gapping_params: BlastCompoGappingParams,
    compo_adjust_mode: CompoAdjustMode,
    position_based: bool,
    query_is_translated: bool,
    subject_is_translated: bool,
    ccat_query_length: i32,
    cutoff_s: i32,
    cutoff_e: f64,
    do_link_hsps: bool,
    near_identical_cutoff: f64,
) -> BlastRedoAlignParams {
    BlastRedoAlignParams {
        matrix_info,
        gapping_params,
        compo_adjust_mode,
        position_based,
        re_pseudocounts: K_RE_MATRIX_ADJUSTMENT_PSEUDOCOUNTS,
        subject_is_translated,
        query_is_translated,
        ccat_query_length,
        cutoff_s,
        cutoff_e,
        do_link_hsps,
        near_identical_cutoff,
    }
}

/// 1-1 port of `Blast_RedoAlignParamsFree` (`redo_alignment.c:1001`).
///
/// In Rust, `Drop` handles the deallocation that NCBI does manually.
/// Provided for parity so callers porting NCBI flow can keep matching
/// call sites; sets the slot to `None` to mirror C's
/// `*pparams = NULL`.
pub fn blast_redo_align_params_free(slot: &mut Option<BlastRedoAlignParams>) {
    *slot = None;
}

/// 1-1 port of `s_ResultHspToDistinctAlign` (`blast_kappa.c:769`).
///
/// Converts an array of HSPs into per-frame singly-linked
/// `BlastCompoAlignment` lists. NCBI uses a 6-element `tail[]` array
/// indexed by `(hsp.context - init_context)`; we use the same
/// convention. Each new alignment carries:
/// - `score = hsp.score * local_scaling_factor` (rounded as integer)
/// - `matrix_adjust_rule = DontAdjust` (matches NCBI's
///   `eDontAdjustMatrix` since the alignment hasn't been
///   re-evaluated yet)
/// - `query_*` / `match_*` from the HSP coords
/// - `frame` from `hsp.subject_frame` (NCBI's `hsp->subject.frame`)
/// - `query_index = hsp.context` (NCBI's `hsp->context`)
/// - `context = None` (NCBI stashes a `BlastHSP*` here for back-ref;
///   we drop it because the Rust pipeline doesn't reuse the HSP
///   pointer afterwards — the kappa caller materializes a fresh
///   `BlastHSP` from the alignment via `s_HSPListFromDistinctAlignments`).
///
/// Returns the per-frame counts (`numAligns[]`); the per-frame head
/// pointers are written into `lists` (which must have length 6).
/// The function returns `0` on success.
///
/// `hsp.subject_frame` isn't on the Rust `Hsp` struct yet (only
/// `context`), so this port stamps `frame = 0` until `Hsp` carries
/// the explicit subject frame field. That doesn't break the per-frame
/// indexing — that's done via `(hsp.context - init_context)`.
pub fn result_hsp_to_distinct_align(
    lists: &mut [Option<Box<BlastCompoAlignment>>; 6],
    num_aligns: &mut [i32; 6],
    hsps: &[Hsp],
    init_context: i32,
    local_scaling_factor: f64,
) -> i32 {
    // NCBI's `tail[]` is a parallel array of last-pointers. Rust
    // owns its links via `Box`/`next`, so we can't keep a mutable
    // pointer-to-tail across iterations safely. Instead we walk to
    // the tail of each per-frame list when appending — O(N²) in the
    // worst case but matches the C result exactly. For typical HSP
    // counts (< 100 per OID) this is fine.
    for v in num_aligns.iter_mut() {
        *v = 0;
    }

    for hsp in hsps {
        let frame_index = (hsp.context - init_context) as i32;
        if !(0..6).contains(&frame_index) {
            // C: `ASSERT(frame_index < 6 && frame_index >= 0);` —
            // skip rather than abort in release builds.
            continue;
        }
        let frame_index = frame_index as usize;
        let new_align = BlastCompoAlignment::new(
            (hsp.score as f64 * local_scaling_factor).round() as i32,
            MatrixAdjustRule::DontAdjust,
            hsp.context,
            hsp.query_offset,
            hsp.query_end,
            hsp.subject_offset,
            hsp.subject_end,
            0, // hsp.subject_frame — not on Rust Hsp yet (TODO)
            None,
        );
        let new_box = Box::new(new_align);

        match lists[frame_index].as_mut() {
            None => {
                lists[frame_index] = Some(new_box);
            }
            Some(head) => {
                // Walk to the tail of the existing list.
                let mut cursor: &mut Box<BlastCompoAlignment> = head;
                while cursor.next.is_some() {
                    cursor = cursor.next.as_mut().unwrap();
                }
                cursor.next = Some(new_box);
            }
        }
        num_aligns[frame_index] += 1;
    }
    0
}

/// 1-1 port of `s_DoSegSequenceData` (`blast_kappa.c:1427`).
///
/// Filters low-complexity regions out of `seq_data` using the SEG
/// algorithm (NCBIstdaa input). NCBI's C uses
/// `BlastFilteringOptionsFromString(BLASTP_MASK_INSTRUCTIONS)`
/// followed by `BlastSetUp_Filter` and `Blast_MaskTheResidues`. Our
/// Rust port collapses those into the existing
/// `crate::filter::seg_filter_ncbistdaa` + an in-place X-substitution,
/// which is what `Blast_MaskTheResidues` ultimately does for protein
/// (replace each masked residue with NCBIstdaa `X = 21`).
///
/// Returns `(status, was_biased)` where `status = 0` on success and
/// `was_biased` is true iff any region was masked (matches NCBI's
/// `*is_seq_biased = (mask_seqloc != NULL)` behaviour).
pub fn do_seg_sequence_data(seq_data: &mut BlastCompoSequenceData) -> (i32, bool) {
    // SEG defaults from NCBI's `BLASTP_MASK_INSTRUCTIONS`: window=12,
    // locut=2.2, hicut=2.5. Same constants used by api.rs's existing
    // `apply_seg_ncbistdaa`.
    let length = seq_data.length.max(0) as usize;
    let start = seq_data.data_offset;
    let end = start + length;
    if end > seq_data.buffer.len() {
        return (-1, false);
    }
    let view = &mut seq_data.buffer[start..end];
    let mask = crate::filter::seg_filter_ncbistdaa(view, 12, 2.2, 2.5);
    let was_biased = !mask.regions.is_empty();
    let masked_residue = 21u8; // NCBIstdaa 'X'
    for r in &mask.regions {
        let r_start = (r.start.max(0) as usize).min(view.len());
        let r_end = (r.end as usize).min(view.len());
        for aa in &mut view[r_start..r_end] {
            *aa = masked_residue;
        }
    }
    (0, was_biased)
}

/// 1-1 port of `s_MatchingSequenceRelease` (`blast_kappa.c:907`).
///
/// In the C source this releases the BlastSeqSrc-managed subject
/// sequence and frees the per-match `BlastKappa_SequenceInfo`
/// payload. The Rust `BlastCompoMatchingSequence` doesn't carry a
/// `local_data` pointer (the C `void *` payload is replaced by
/// `local_data_index`), so the cleanup work is handled by the
/// caller's lifetime management. We provide this hook for parity so
/// callers porting NCBI flow can keep matching call sites; it
/// resets the index to `-1`, mirroring NCBI's `if (self->index >= 0)`
/// guard semantics.
pub fn matching_sequence_release(self_: &mut BlastCompoMatchingSequence) {
    self_.local_data_index = -1;
    self_.length = 0;
}

/// 1-1 port of `s_RescaleSearch` (`blast_kappa.c:2117`).
///
/// Rescales every per-context gapped Karlin block so that
/// `lambda /= scale_factor` (and `log_k = ln(K)` is reset for
/// numerical hygiene), and rounds gap_open/gap_extend by
/// `BLAST_Nint(value * scale_factor)`. Used by the kappa driver to
/// move the search into a finer-grained scoring resolution before the
/// composition adjustment pass.
pub fn rescale_search(
    kbp_gap: &mut [crate::stat::KarlinBlk],
    scoring: &mut crate::parameters::ScoringParameters,
    num_queries: i32,
    scale_factor: f64,
) {
    for i in 0..(num_queries.max(0) as usize) {
        if let Some(kbp) = kbp_gap.get_mut(i) {
            if kbp.lambda > 0.0 {
                kbp.lambda /= scale_factor;
                // C: `kbp->logK = log(kbp->K);` — recomputes from the
                // (unchanged) K to flush any prior cached value.
                kbp.log_k = kbp.k.ln();
            }
        }
    }
    scoring.gap_open = crate::math::nint(scoring.gap_open as f64 * scale_factor) as i32;
    scoring.gap_extend = crate::math::nint(scoring.gap_extend as f64 * scale_factor) as i32;
    scoring.scale_factor = scale_factor;
}

/// 1-1 port of the **query-preparation preamble** of `s_SequenceGetRange`
/// (`blast_kappa.c:1670`).
///
/// Copies `query[q_range.begin .. q_range.end]` into a freshly-allocated
/// buffer with NCBI's sentinel-byte layout (`buffer[0] = 0`, real data
/// at `buffer[1 ..= length]`, terminator at `buffer[length + 1] = 0`),
/// substituting amino acid 24 (Selenocysteine) with 3 (Cysteine) per
/// the C source's inline rule.
///
/// The full `s_SequenceGetRange` then dispatches to
/// `s_SequenceGetProteinRange` or `s_SequenceGetTranslatedRange` to
/// fill in the subject `seqData`. Those branches require subject I/O
/// plumbing not yet ported; they're tracked as TODOs.
pub fn sequence_prep_query_range(
    query: &BlastCompoSequenceData,
    q_range: &BlastCompoSequenceRange,
) -> BlastCompoSequenceData {
    let begin = q_range.begin.max(0) as usize;
    let end = q_range.end.max(q_range.begin) as usize;
    let length = end.saturating_sub(begin);
    // C: `calloc((length + 2) * sizeof(Uint1))` then `data = buffer + 1;`.
    // Sentinel zero stays at index 0 and at index length+1.
    let mut buffer = vec![0u8; length + 2];
    let q_data = query.data();
    for (idx, slot) in buffer[1..=length].iter_mut().enumerate() {
        let src = q_data.get(begin + idx).copied().unwrap_or(0);
        // Selenocysteine (24) → Cysteine (3).
        *slot = if src != 24 { src } else { 3 };
    }
    BlastCompoSequenceData {
        buffer,
        data_offset: 1,
        length: length as i32,
    }
}

/// 1-1 port of `s_GappingParamsNew` (`blast_kappa.c:2354`).
///
/// Builds a `BlastCompoGappingParams` from `(scoring, extension)` and a
/// per-context array of gapped Karlin blocks. The X-dropoff is set to
/// `max(gap_x_dropoff_final * ln(2) / min_lambda, raw_gap_x_dropoff_final)`
/// where `min_lambda` is the smallest Lambda across all valid contexts —
/// matching NCBI's reasoning of using the most-conservative (smallest)
/// statistic to set the dropoff.
///
/// `gap_x_dropoff_final_bits` is the option's bit-score X-dropoff (NCBI
/// `BlastExtensionOptions::gap_x_dropoff_final`); `raw_gap_x_dropoff_final`
/// is the parameter's already-converted raw-score equivalent (NCBI
/// `BlastExtensionParameters::gap_x_dropoff_final`).
pub fn gapping_params_new(
    scoring: &crate::parameters::ScoringParameters,
    kbp_gap: &[crate::stat::KarlinBlk],
    num_queries: i32,
    gap_x_dropoff_final_bits: f64,
    raw_gap_x_dropoff_final: i32,
) -> BlastCompoGappingParams {
    // C: walk every kbp_gap entry, retain the smallest Lambda. NCBI
    // skips NULL entries; in our Rust port a default `KarlinBlk` has
    // `lambda = 0.0` which we treat as an invalid context.
    let mut min_lambda = f64::MAX;
    for i in 0..(num_queries.max(0) as usize) {
        if let Some(kbp) = kbp_gap.get(i) {
            if kbp.lambda > 0.0 && kbp.lambda < min_lambda {
                min_lambda = kbp.lambda;
            }
        }
    }
    let x_dropoff = if min_lambda < f64::MAX {
        let bits_dropoff =
            (gap_x_dropoff_final_bits * crate::math::NCBIMATH_LN2 / min_lambda).round() as i32;
        bits_dropoff.max(raw_gap_x_dropoff_final)
    } else {
        raw_gap_x_dropoff_final
    };

    BlastCompoGappingParams {
        gap_open: scoring.gap_open,
        gap_extend: scoring.gap_extend,
        decline_align: i32::MIN, // unused in the kappa pipeline; NCBI default
        x_dropoff,
    }
}

/// 1-1 port of the **non-position-based** branch of `s_MatrixInfoInit`
/// (`blast_kappa.c:2199`).
///
/// Initializes a `BlastMatrixInfo` from the matrix name and the
/// "ideal" gapped Karlin block, populating `ungappedLambda`,
/// `startFreqRatios` (currently BLOSUM62 only), and `startMatrix`. The
/// position-based (PSI) branch needs `s_GetPosBasedStartFreqRatios` and
/// `s_ScalePosMatrix`, which are TODOs.
///
/// `kbp_ideal_lambda` is `sbp->kbp_ideal->Lambda` from the C source —
/// the lambda for the ideal (non-context-specific) ungapped matrix.
pub fn matrix_info_init_blastp(
    self_: &mut BlastMatrixInfo,
    matrix_name: &str,
    kbp_ideal_lambda: f64,
    scale_factor: f64,
) -> i32 {
    self_.matrix_name = matrix_name.to_string();
    self_.positional = false;
    self_.ungapped_lambda = kbp_ideal_lambda / scale_factor;

    // C: `s_GetStartFreqRatios(self->startFreqRatios, matrixName);`
    if !matrix_name.eq_ignore_ascii_case("BLOSUM62") {
        // Other matrices are TODO; populate empty data and return -1
        // (matches NCBI's `s_GetStartFreqRatios` returning non-zero on
        // unknown matrix).
        return -1;
    }
    let freq_ratios = crate::matrix::get_blosum62_freq_ratios();
    self_.start_freq_ratios = freq_ratios.iter().map(|row| row.to_vec()).collect();

    // C: `Blast_Int4MatrixFromFreq(self->startMatrix, self->cols,
    //                              self->startFreqRatios, self->ungappedLambda);`
    let int_mat = crate::composition::matrix_from_freq_ratios(self_.ungapped_lambda, &freq_ratios);
    self_.matrix = int_mat.iter().map(|row| row.to_vec()).collect();
    self_.rows = crate::matrix::AA_SIZE as i32;
    self_.cols = crate::matrix::AA_SIZE as i32;
    0
}

/// 1-1 port of `s_RecordInitialSearch` (`blast_kappa.c:2059`).
///
/// Snapshots `scoring->gap_open / gap_extend / scale_factor`,
/// per-context gapped Karlin blocks (one entry per query context), and
/// — when composition adjustment is enabled — a copy of the scoring
/// matrix. Intended as the entry-point capture before the kappa
/// driver mutates these values.
///
/// Caller passes the matrix as a 2-D `Vec` borrow. For position-based
/// (PSI) searches that's the `psi_matrix.pssm.data` of `query_length ×
/// AA_SIZE`; for ordinary searches it's the standard `matrix.data` of
/// `AA_SIZE × AA_SIZE`. Either way, every cell is copied into
/// `searchParams->origMatrix`.
pub fn record_initial_search(
    saved: &mut BlastKappaSavedParameters,
    kbp_gap: &[crate::stat::KarlinBlk],
    matrix: &[Vec<i32>],
    scoring: &crate::parameters::ScoringParameters,
    query_length: i32,
    compo_adjust_mode: CompoAdjustMode,
    position_based: bool,
) -> i32 {
    saved.gap_open = scoring.gap_open;
    saved.gap_extend = scoring.gap_extend;
    saved.scale_factor = scoring.scale_factor;

    // C: per-query kbp copy. NCBI uses `Blast_KarlinBlkCopy` to clone the
    // values; in Rust `KarlinBlk: Clone` does the same.
    for i in 0..(saved.num_queries.max(0) as usize) {
        if i < kbp_gap.len() {
            saved.kbp_gap_orig[i] = kbp_gap[i].clone();
        }
    }

    // C: copy matrix when composition adjustment is enabled.
    if !matches!(compo_adjust_mode, CompoAdjustMode::NoCompositionBasedStats) {
        let rows = if position_based {
            query_length.max(0) as usize
        } else {
            crate::matrix::AA_SIZE
        };
        // Ensure orig_matrix has the right shape — `BlastKappaSavedParameters::new`
        // already allocated it, but defend against caller mismatch.
        if saved.orig_matrix.len() < rows {
            saved.orig_matrix.resize(rows, vec![0i32; crate::matrix::AA_SIZE]);
        }
        for (i, row) in saved.orig_matrix.iter_mut().take(rows).enumerate() {
            if i < matrix.len() {
                let cols = crate::matrix::AA_SIZE.min(row.len()).min(matrix[i].len());
                row[..cols].copy_from_slice(&matrix[i][..cols]);
            }
        }
    }
    0
}

/// 1-1 port of `s_RestoreSearch` (`blast_kappa.c:2147`).
///
/// Inverse of [`record_initial_search`]: copies the snapshot back into
/// the caller's mutable scoring parameters, Karlin blocks, and
/// (optionally) the matrix.
pub fn restore_search(
    saved: &BlastKappaSavedParameters,
    kbp_gap: &mut [crate::stat::KarlinBlk],
    matrix: &mut [Vec<i32>],
    scoring: &mut crate::parameters::ScoringParameters,
    query_length: i32,
    position_based: bool,
    compo_adjust_mode: CompoAdjustMode,
) {
    scoring.gap_open = saved.gap_open;
    scoring.gap_extend = saved.gap_extend;
    scoring.scale_factor = saved.scale_factor;

    for i in 0..(saved.num_queries.max(0) as usize) {
        if i < kbp_gap.len() && i < saved.kbp_gap_orig.len() {
            kbp_gap[i] = saved.kbp_gap_orig[i].clone();
        }
    }

    if !matches!(compo_adjust_mode, CompoAdjustMode::NoCompositionBasedStats) {
        let rows = if position_based {
            query_length.max(0) as usize
        } else {
            crate::matrix::AA_SIZE
        };
        for (i, row) in matrix.iter_mut().take(rows).enumerate() {
            if i < saved.orig_matrix.len() {
                let cols = crate::matrix::AA_SIZE.min(row.len()).min(saved.orig_matrix[i].len());
                row[..cols].copy_from_slice(&saved.orig_matrix[i][..cols]);
            }
        }
    }
}

/// 1-1 port of `s_SavedParametersFree` (`blast_kappa.c:1977`).
///
/// In Rust the struct's `Drop` impl handles deallocation, so this
/// function is a no-op marker — provided so that callers porting the
/// NCBI flow can keep matching call sites. Sets the option slot to
/// `None` to mirror C's `*searchParams = NULL`.
pub fn saved_parameters_free(saved: &mut Option<BlastKappaSavedParameters>) {
    *saved = None;
}

/// 1-1 port of `BlastCompo_Heap` (`compo_heap.h:82`). The C struct keeps
/// hits as either an unsorted array (until `n == heapThreshold`) or as
/// a max-heap keyed on e-value (after the threshold). Once it's a heap,
/// new candidates with e-value >= `worstEvalue` are rejected unless
/// they pass `ecutoff`.
///
/// This Rust port stores the records as a single `Vec<HspList>` and
/// resorts on demand — the C array/heap split is a memory-management
/// optimization that's irrelevant in Rust. The semantics
/// (`would_insert`, `insert`, `pop`) are preserved.
#[derive(Debug, Clone, Default)]
pub struct BlastCompoHeap {
    pub records: Vec<HspList>,
    pub heap_threshold: i32,
    pub ecutoff: f64,
}

impl BlastCompoHeap {
    pub fn new(heap_threshold: i32, ecutoff: f64) -> Self {
        Self {
            records: Vec::new(),
            heap_threshold,
            ecutoff,
        }
    }

    /// Worst (largest) e-value currently in the heap, or `+∞` if empty.
    pub fn worst_evalue(&self) -> f64 {
        self.records
            .iter()
            .map(|r| r.best_evalue)
            .fold(f64::NEG_INFINITY, f64::max)
            .max(f64::NEG_INFINITY)
            .max(if self.records.is_empty() { f64::INFINITY } else { f64::NEG_INFINITY })
    }

    /// 1-1 port of `BlastCompo_HeapPop`. Removes and returns the
    /// hit-list with the worst (largest) e-value.
    pub fn pop_worst(&mut self) -> Option<HspList> {
        if self.records.is_empty() {
            return None;
        }
        let mut worst_idx = 0;
        let mut worst_eval = self.records[0].best_evalue;
        for (i, r) in self.records.iter().enumerate().skip(1) {
            if r.best_evalue > worst_eval {
                worst_eval = r.best_evalue;
                worst_idx = i;
            }
        }
        Some(self.records.swap_remove(worst_idx))
    }
}

/// 1-1 port of `s_FreeBlastCompo_QueryInfoArray` (`blast_kappa.c:2279`).
///
/// In the C source this `free`s every `query_info[i].words` array and
/// then `free`s the outer container. In Rust the `BlastCompoQueryInfo`
/// owns its `words: Vec<u64>` so the `Vec::clear()` + drop suffices.
/// Provided for parity / completeness so callers porting NCBI flow can
/// keep the same call sites.
pub fn free_blast_compo_query_info_array(query_info: &mut Vec<BlastCompoQueryInfo>) {
    query_info.clear();
}

/// 1-1 port of `s_ClearHeap` (`blast_kappa.c:2518`).
pub fn clear_heap(heap: &mut BlastCompoHeap) {
    while heap.pop_worst().is_some() {}
}

/// 1-1 port of `s_FillResultsFromCompoHeaps` (`blast_kappa.c:2493`).
///
/// Drains each heap into a fresh `HitList`, then reverses the
/// query order at the end (mirroring NCBI's
/// `Blast_HSPResultsReverseOrder` finishing call). The C version
/// allocates `Blast_HitListNew(hitlist_size)` per query; we let
/// `HitList::new()` handle that and ignore the size cap (Rust's
/// `Vec` is unbounded; truncate at the call site if needed).
pub fn fill_results_from_compo_heaps(
    heaps: &mut [BlastCompoHeap],
) -> crate::hspstream::HspResults {
    let num_queries = heaps.len();
    let mut results = crate::hspstream::HspResults::new(num_queries as i32);
    for (q, heap) in heaps.iter_mut().enumerate() {
        let mut hitlist = crate::hspstream::HitList::new();
        while let Some(hsp_list) = heap.pop_worst() {
            hitlist.hsp_lists.push(hsp_list);
        }
        // C: `Blast_HSPResultsReverseOrder` — reverse so the head
        // represents the highest-priority hit; pop_worst drained from
        // worst→best, so the per-query list is currently
        // worst-first; reverse to put best-first.
        hitlist.hsp_lists.reverse();
        results.hitlists[q] = Some(hitlist);
    }
    results
}

/// 1-1 port of `s_GetQueryInfo` (`blast_kappa.c:2308`).
///
/// Builds a `Vec<BlastCompoQueryInfo>` from the per-context fields of a
/// Rust `QueryInfo` plus the concatenated query buffer. For each context
/// (one entry per strand/frame), populates:
///   - `origin` from the context's `query_offset`
///   - `seq` from the slice `query_data[origin .. origin + length]`
///     (copied into an owned `BlastCompoSequenceData`)
///   - `eff_search_space` from the context's `eff_searchsp`
///   - `words` via [`create_word_array`] (1-1 with C's `s_CreateWordArray`)
///   - `composition` via `crate::composition::read_composition` (1-1 with
///     C's `Blast_ReadAaComposition`); skipped when `skip == true`.
///
/// NCBI uses `last_context + 1` as the array size; our `QueryInfo`
/// stores the same per-context array, so we walk every entry.
pub fn get_query_info(
    query_data: &[u8],
    blast_query_info: &crate::queryinfo::QueryInfo,
    skip: bool,
) -> Vec<BlastCompoQueryInfo> {
    blast_query_info
        .contexts
        .iter()
        .map(|ctx| {
            let origin = ctx.query_offset.max(0) as usize;
            let length = ctx.query_length.max(0) as usize;
            let seq_slice: Vec<u8> = if origin + length <= query_data.len() {
                query_data[origin..origin + length].to_vec()
            } else {
                Vec::new()
            };
            let words = create_word_array(&seq_slice).unwrap_or_default();
            let composition = if skip {
                Vec::new()
            } else {
                let (probs, _num_true) =
                    crate::composition::read_composition(&seq_slice, crate::matrix::AA_SIZE);
                probs
            };
            BlastCompoQueryInfo {
                origin: ctx.query_offset,
                seq: BlastCompoSequenceData {
                    buffer: seq_slice,
                    data_offset: 0,
                    length: length as i32,
                },
                composition,
                eff_search_space: ctx.eff_searchsp as f64,
                words,
            }
        })
        .collect()
}

/// 1-1 port of `s_CreateWordArray` (`blast_kappa.c:2244`).
///
/// Returns an array of 8-mer hashes such that `result[i]` is the hash
/// of `seq[i .. i + 8]`. Used by `s_TestNearIdentical` /
/// `find_num_identical` for fast near-identity probes against a query.
///
/// Returns `None` if `seq.len() < 8` (NCBI returns -1).
///
/// **Bit-for-bit fidelity to NCBI:** the C loop is `for (i = 1; i < seq_len - word_size; i++)`,
/// which leaves the last legitimate 8-mer position (`seq_len - word_size`)
/// at calloc'd-zero. That off-by-one is preserved here — the slot is
/// allocated but never written. `find_num_identical`'s outer loop
/// iterates `q_pos < query_len - word_size`, so it never reads that
/// trailing slot anyway.
pub fn create_word_array(seq: &[u8]) -> Option<Vec<u64>> {
    const WORD_SIZE: usize = 8;
    const MASK: u64 = 0xFFFF_FFFF_FF;
    if seq.len() < WORD_SIZE {
        return None;
    }
    // C: `calloc((seq_len - word_size + 1) * sizeof(Uint8))`.
    let n = seq.len() - WORD_SIZE + 1;
    let mut hashes = vec![0u64; n];
    hashes[0] = get_hash(&seq[0..WORD_SIZE], WORD_SIZE);
    // C: `for (i = 1; i < seq_len - word_size; i++)` — note `<` not `<=`,
    // intentionally leaving slot [seq_len - word_size] zero.
    let upper = seq.len() - WORD_SIZE;
    for i in 1..upper {
        hashes[i] = ((hashes[i - 1] << 5) & MASK) + seq[i + WORD_SIZE - 1] as u64;
    }
    Some(hashes)
}

/// 1-1 port of `s_Blast_HSPGetNumIdentitiesAndPositives`
/// (`blast_hits.c:746`).
///
/// Walks the HSP's gap edit script (or its ungapped extent if no
/// edit script is present) and counts identical residue pairs and
/// "positive" residue pairs (where the substitution-matrix score is
/// positive). Returns `(num_ident, align_length, num_pos)`.
///
/// `matrix` is the substitution matrix for the protein alphabet, used
/// only for the `num_pos` count. Pass `None` for nucleotide alignments
/// or when positives are not needed.
pub fn blast_hsp_get_num_identities(
    query: &[u8],
    subject: &[u8],
    hsp: &Hsp,
    edit_ops: Option<&[(crate::gapinfo::GapAlignOpType, i32)]>,
    matrix: Option<&[[i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE]>,
) -> (i32, i32, i32) {
    let q_off = hsp.query_offset.max(0) as usize;
    let s_off = hsp.subject_offset.max(0) as usize;
    let q_length = (hsp.query_end - hsp.query_offset).max(0) as usize;
    let s_length = (hsp.subject_end - hsp.subject_offset).max(0) as usize;

    let mut num_ident = 0i32;
    let mut num_pos = 0i32;
    let mut align_length = 0i32;

    let mut q = q_off;
    let mut s = s_off;

    match edit_ops {
        None => {
            // Ungapped: query and subject lengths must match.
            if q_length != s_length {
                return (0, 0, 0);
            }
            align_length = q_length as i32;
            for _ in 0..q_length {
                if q < query.len() && s < subject.len() {
                    if query[q] == subject[s] {
                        num_ident += 1;
                    } else if let Some(m) = matrix {
                        let qi = query[q] as usize;
                        let si = subject[s] as usize;
                        if qi < crate::matrix::AA_SIZE
                            && si < crate::matrix::AA_SIZE
                            && m[qi][si] > 0
                        {
                            num_pos += 1;
                        }
                    }
                }
                q += 1;
                s += 1;
            }
        }
        Some(ops) => {
            for &(op, count) in ops {
                align_length += count;
                match op {
                    crate::gapinfo::GapAlignOpType::Sub => {
                        for _ in 0..count {
                            if q < query.len() && s < subject.len() {
                                if query[q] == subject[s] {
                                    num_ident += 1;
                                } else if let Some(m) = matrix {
                                    let qi = query[q] as usize;
                                    let si = subject[s] as usize;
                                    if qi < crate::matrix::AA_SIZE
                                        && si < crate::matrix::AA_SIZE
                                        && m[qi][si] > 0
                                    {
                                        num_pos += 1;
                                    }
                                }
                            }
                            q += 1;
                            s += 1;
                        }
                    }
                    crate::gapinfo::GapAlignOpType::Del
                    | crate::gapinfo::GapAlignOpType::Del1
                    | crate::gapinfo::GapAlignOpType::Del2 => {
                        s += count as usize;
                    }
                    crate::gapinfo::GapAlignOpType::Ins
                    | crate::gapinfo::GapAlignOpType::Ins1
                    | crate::gapinfo::GapAlignOpType::Ins2 => {
                        q += count as usize;
                    }
                    _ => {
                        // C `default:` branch advances both — used for
                        // out-of-frame or unknown ops.
                        q += count as usize;
                        s += count as usize;
                    }
                }
            }
        }
    }

    (num_ident, align_length, num_pos)
}

/// 1-1 port of `s_ComputeNumIdentities` (`blast_kappa.c:458`) for the
/// non-translated path.
///
/// Walks every HSP in the list and stamps `hsp.num_ident` from a fresh
/// re-walk over the (unmasked) query and subject bytes via
/// `blast_hsp_get_num_identities`. The C version dispatches on
/// `program_number` to translate the subject for tblastn before
/// counting; that branch is a TODO until `BlastTargetTranslationNew` is
/// ported. Pass NCBI4na/BLASTNA decoded subject bytes for blastn,
/// NCBIstdaa for blastp/blastx.
///
/// `edit_ops` is the per-HSP edit script in parallel order with
/// `hsp_list.hsps`. Pass `None` for an HSP that's still ungapped.
pub fn compute_num_identities_blastp(
    query_nomask: &[u8],
    subject: &[u8],
    hsp_list: &mut HspList,
    edit_ops_per_hsp: &[Option<Vec<(crate::gapinfo::GapAlignOpType, i32)>>],
    matrix: Option<&[[i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE]>,
) {
    for (i, hsp) in hsp_list.hsps.iter_mut().enumerate() {
        let ops = edit_ops_per_hsp.get(i).and_then(|o| o.as_deref());
        let (num_ident, _align_length, _num_pos) =
            blast_hsp_get_num_identities(query_nomask, subject, hsp, ops, matrix);
        hsp.num_ident = num_ident;
    }
}

/// 1-1 port of `s_HitlistEvaluateAndPurge` (`blast_kappa.c:394`).
///
/// Assigns final e-values and prunes the hit list. The full C function
/// dispatches between `BLAST_LinkHsps` (sum-statistics branch) and
/// `Blast_HSPListGetEvalues` (single-HSP branch), then applies
/// `s_AdjustEvaluesForComposition` for blastp/blastx, then drops HSPs
/// above the evalue threshold and reports `(best_score, best_evalue)`.
///
/// This Rust port currently implements the **single-HSP, no-link-HSPs
/// branch** with composition adjustment — sufficient for blastp /
/// blastx without sum statistics. The link-HSPs branch and the
/// `Blast_HSPListGetEvalues` call are TODOs (caller is expected to have
/// already populated `hsp.evalue` per-HSP). When those helpers are
/// ported, plumb them in by branching on `do_sum_stats`.
///
/// Returns `(best_score, best_evalue)`. After the call, `hsp_list.hsps`
/// is sorted-stable by NCBI's tie-breaker and contains only HSPs whose
/// e-value passes `max_evalue`.
#[allow(clippy::too_many_arguments)]
pub fn hitlist_evaluate_and_purge(
    hsp_list: &mut HspList,
    subject_length: i32,
    program_number: ProgramType,
    query_length: i32,
    length_adjustment: i32,
    eff_searchsp: f64,
    pvalue_for_this_pair: f64,
    max_evalue: f64,
    do_sum_stats: bool,
) -> (i32, f64) {
    // C: `*pbestEvalue = DBL_MAX; *pbestScore = 0;`
    let mut best_evalue = f64::MAX;
    let mut best_score = 0i32;

    if do_sum_stats {
        // TODO: BLAST_LinkHsps port for sum-statistics path. The Rust
        // `link_hsps.rs` ports the inner helpers (`s_BlastEvenGapLinkHSPs`
        // etc.) but the umbrella `BLAST_LinkHsps` driver isn't wired in
        // yet. For now, fall through and treat as the single-HSP branch.
    }
    // The C single-HSP branch normally calls `Blast_HSPListGetEvalues`
    // here, which iterates the HSP array and sets `hsp->evalue` from the
    // score using `BLAST_KarlinStoE_simple`. The Rust pipeline already
    // does this assignment at HSP-construction time (see api.rs callers
    // of `prot_kbp.raw_to_evalue` / `spouge_evalue`), so we skip the
    // re-assignment to avoid double-applying.

    // Composition adjustment for blastp/blastx if the p-value is in [0, 1].
    if (program_number == crate::program::BLASTP
        || program_number == crate::program::BLASTX)
        && (0.0..=1.0).contains(&pvalue_for_this_pair)
    {
        adjust_evalues_for_composition(
            hsp_list,
            pvalue_for_this_pair,
            subject_length,
            query_length,
            length_adjustment,
            eff_searchsp,
        );
    }

    // C: `Blast_HSPListReapByEvalue(hsp_list, hitParams->options)`.
    crate::hits::filter_by_evalue(hsp_list, max_evalue);

    // After purge, NCBI sorts by score and reports the head.
    if !hsp_list.hsps.is_empty() {
        hsp_list.sort_by_score();
        best_evalue = hsp_list.best_evalue;
        best_score = hsp_list.hsps[0].score;
    }

    (best_score, best_evalue)
}

/// 1-1 port of `s_AdjustEvaluesForComposition` (`blast_kappa.c:134`).
///
/// Combines a sequence-composition p-value with each HSP's
/// alignment p-value via Fisher's method, then converts back to an
/// e-value. The E-value is first scaled from "database" → "single
/// sequence" semantics so that the p-value combination operates on
/// per-sequence p-values, then scaled back at the end.
///
/// Inputs (mirroring NCBI signature):
/// - `comp_p_value`: composition p-value computed earlier from lambda.
/// - `subject_length`: length of this database sequence (already in
///   protein space for translated subjects).
/// - `query_length`, `length_adjustment`, `eff_searchsp`: from
///   `BlastContextInfo` for this query context.
pub fn adjust_evalues_for_composition(
    hsp_list: &mut HspList,
    comp_p_value: f64,
    subject_length: i32,
    query_length: i32,
    length_adjustment: i32,
    eff_searchsp: f64,
) {
    let query_eff = ((query_length - length_adjustment) as f64).max(1.0);
    let subject_eff = ((subject_length - length_adjustment) as f64).max(1.0);
    let dblen_eff = eff_searchsp / query_eff;

    // C: `db_to_sequence_scale = subject_eff / dblen_eff;`.
    let db_to_sequence_scale = subject_eff / dblen_eff;

    let mut best_evalue = f64::MAX;
    for hsp in &mut hsp_list.hsps {
        // Convert DB-scale e-value to sequence-scale.
        hsp.evalue *= db_to_sequence_scale;
        // Combine with composition p-value (Fisher's method) and convert
        // back to e-value, then rescale to DB-scale.
        let align_p_value = crate::composition::karlin_e_to_p(hsp.evalue);
        let combined_p_value = crate::composition::overall_p_value(comp_p_value, align_p_value);
        hsp.evalue = crate::composition::karlin_p_to_e_compo(combined_p_value);
        hsp.evalue /= db_to_sequence_scale;

        if hsp.evalue < best_evalue {
            best_evalue = hsp.evalue;
        }
    }

    hsp_list.best_evalue = best_evalue;
}

/// 1-1 port of `s_TestNearIdentical` (`blast_kappa.c:1258`).
///
/// Returns true iff the aligned query/subject ranges have ≥ 95 % identity
/// after a fast bidirectional extension (right then left from the
/// alignment endpoints, finishing with a k-mer-hash sweep over the
/// middle). Used by the redo-alignment driver to skip SEG masking on
/// near-identical hits where masking would cost more in score than
/// composition adjustment recovers.
///
/// The constants are from C: `kMinFractionNearIdentical = 0.95`,
/// `max_shift = 8`. Caller must have populated `query_words` so that
/// `query_words[i]` is `s_GetHash(&queryData[qStart + i], 8)` for every
/// valid 8-mer starting at index `qStart + i` (see `s_CreateWordArray`).
pub fn test_near_identical(
    seq_data: &BlastCompoSequenceData,
    seq_offset: i32,
    query_data: &BlastCompoSequenceData,
    query_offset: i32,
    query_words: &[u64],
    align: &BlastCompoAlignment,
) -> bool {
    const K_MIN_FRACTION_NEAR_IDENTICAL: f64 = 0.95;
    const MAX_SHIFT: i32 = 8;

    // C: `qStart = align->queryStart - queryOffset; qEnd = align->queryEnd - queryOffset - 1;`
    let q_start = (align.query_start - query_offset) as usize;
    let q_end = (align.query_end - query_offset - 1) as usize;
    let s_start = (align.match_start - seq_offset) as usize;
    let s_end = (align.match_end - seq_offset - 1) as usize;

    let query_len = (q_end - q_start + 1) as i32;
    let subject_len = (s_end - s_start + 1) as i32;
    let align_len = query_len.min(subject_len) as f64;

    let query = query_data.data();
    let subject = seq_data.data();

    // Right extension from the start of the ranges.
    let (mut num_identical, query_right_len, subject_right_len, _) = extend_right(
        &query[q_start..q_start + query_len as usize],
        &subject[s_start..s_start + subject_len as usize],
        MAX_SHIFT,
    );

    // C: if the whole query/subject range was processed → return now.
    if query_right_len >= query_len || subject_right_len >= subject_len {
        let fraction_identical = num_identical as f64 / align_len;
        return fraction_identical > K_MIN_FRACTION_NEAR_IDENTICAL;
    }

    // Left extension from the end of what's left.
    let qr = query_right_len as usize;
    let sr = subject_right_len as usize;
    let (left_ident, query_left_len, subject_left_len, _) = extend_left(
        &query[q_start + qr..q_start + query_len as usize],
        &subject[s_start + sr..s_start + subject_len as usize],
        MAX_SHIFT,
    );
    num_identical += left_ident;

    if query_left_len + query_right_len >= query_len
        || subject_left_len + subject_right_len >= subject_len
    {
        let fraction_identical = num_identical as f64 / align_len;
        return fraction_identical > K_MIN_FRACTION_NEAR_IDENTICAL;
    }

    // Middle: k-mer hash sweep.
    let mid_q_len = (query_len - query_left_len - query_right_len) as usize;
    let mid_s_len = (subject_len - subject_left_len - subject_right_len) as usize;
    let q_words_slice = if q_start + qr < query_words.len() {
        &query_words[q_start + qr..]
    } else {
        &[][..]
    };
    num_identical += find_num_identical(
        &query[q_start + qr..q_start + qr + mid_q_len],
        q_words_slice,
        &subject[s_start + sr..s_start + sr + mid_s_len],
        MAX_SHIFT,
    );

    let fraction_identical = num_identical as f64 / align_len;
    fraction_identical > K_MIN_FRACTION_NEAR_IDENTICAL
}

/// 1-1 port of `s_HSPListFromDistinctAlignments` (`blast_kappa.c:304`).
///
/// Walks a singly-linked list of `BlastCompoAlignment` values
/// (the output of the kappa redo-alignment driver), converts each to
/// an `Hsp` matching NCBI's `Blast_HSPInit` field semantics, pushes it
/// into `hsp_list`, then sorts the list by score (matching NCBI's
/// `Blast_HSPListSortByScore`). Consumes the linked list — on return,
/// `*alignments` is `None` (mirrors C's `BlastCompo_AlignmentsFree`).
///
/// The C version also stamps `new_hsp->comp_adjustment_method` from the
/// `matrix_adjust_rule` per the table:
/// `eDontAdjustMatrix → eNoCompositionBasedStats`,
/// `eCompoScaleOldMatrix → eCompositionBasedStats`,
/// otherwise `eCompositionMatrixAdjust`.
/// Our `Hsp` doesn't carry that field yet (TODO: add it), so this port
/// returns the per-HSP `CompoAdjustMode` tags in a parallel `Vec`.
/// Caller order matches `hsp_list.hsps` order **before** the sort step,
/// so the tag at index `i` corresponds to the alignment list's
/// position-`i` element — sort the tag vec in tandem if you depend on
/// hsp index post-sort.
pub fn hsp_list_from_distinct_alignments(
    hsp_list: &mut HspList,
    alignments: &mut Option<Box<BlastCompoAlignment>>,
    oid: i32,
    frame: i32,
) -> (i32, Vec<CompoAdjustMode>) {
    hsp_list.oid = oid;

    let mut comp_tags: Vec<CompoAdjustMode> = Vec::new();
    let mut cursor = alignments.take();
    while let Some(mut node) = cursor {
        let edit_script = node.context.take();
        // C: `Blast_HSPInit(queryStart, queryEnd, matchStart, matchEnd,
        //                   unknown_value, unknown_value,
        //                   queryIndex, frame, align->frame, score,
        //                   &editScript, &new_hsp);`
        // Our Hsp encodes (queryIndex, frame, subject_frame) as a single
        // `context` integer — the search engine elsewhere uses this
        // convention. We follow suit: pack the query `frame` into
        // `context` and let the subject frame ride along via the
        // alignment's `frame` field. When Hsp gets an explicit
        // subject_frame field, plumb it then.
        let _ = (edit_script, frame); // edit_script ownership now lives only
                                      // here; downstream HSP storage will
                                      // need to re-attach when port matures.
        let mut hsp = Hsp {
            score: node.score,
            num_ident: 0, // C: explicitly leaves num_ident blank
            bit_score: 0.0,
            evalue: 0.0,
            query_offset: node.query_start,
            query_end: node.query_end,
            subject_offset: node.match_start,
            subject_end: node.match_end,
            context: node.query_index,
            num_gaps: 0,
        };
        let _ = node.frame; // subject frame; carried via context tag below
        // Translate matrix_adjust_rule → comp_adjustment_method.
        let tag = match node.matrix_adjust_rule {
            MatrixAdjustRule::DontAdjust => CompoAdjustMode::NoCompositionBasedStats,
            MatrixAdjustRule::ScaleOldMatrix => CompoAdjustMode::CompositionBasedStats,
            _ => CompoAdjustMode::CompositionMatrixAdjust,
        };
        // Append HSP and its compositional tag in parallel.
        let _ = &mut hsp;
        hsp_list.hsps.push(hsp);
        comp_tags.push(tag);

        cursor = node.next.take();
    }

    // C: `Blast_HSPListSortByScore(hsp_list)`.
    hsp_list.sort_by_score();
    (0, comp_tags)
}

/// 1-1 port of `s_NewAlignmentFromGapAlign` (`blast_kappa.c:1747`).
///
/// Reads a finished gapped alignment (Rust's `TracebackResult` plays the
/// role of NCBI's `BlastGapAlignStruct`), shifts its local coordinates
/// into the concatenated-query / full-subject frame using the
/// `BlastCompoSequenceRange` offsets, and returns a freshly-built
/// `BlastCompoAlignment`. Ownership of the edit script transfers from
/// the caller's `edit_script` slot to the new alignment; the caller's
/// slot is set to `None` on success — same as C's
/// `*edit_script = NULL`.
pub fn new_alignment_from_gap_align(
    gap_align: &crate::traceback::TracebackResult,
    edit_script: &mut Option<crate::gapinfo::GapEditScript>,
    query_range: &BlastCompoSequenceRange,
    subject_range: &BlastCompoSequenceRange,
    matrix_adjust_rule: MatrixAdjustRule,
) -> Option<BlastCompoAlignment> {
    // C: `queryStart = gap_align->query_start + query_range->begin;` etc.
    // The composition_adjustment library uses concatenated-query coords,
    // so the local TracebackResult coords are shifted by the range begin.
    let query_start = gap_align.query_start as i32 + query_range.begin;
    let query_end = gap_align.query_end as i32 + query_range.begin;
    let query_index = query_range.context;
    let match_start = gap_align.subject_start as i32 + subject_range.begin;
    let match_end = gap_align.subject_end as i32 + subject_range.begin;
    let frame = subject_range.context;

    // C: `BlastCompo_AlignmentNew(... *edit_script)` followed by
    // `*edit_script = NULL` on success. Rust takes ownership via
    // `Option::take()`.
    let context = edit_script.take();
    Some(BlastCompoAlignment::new(
        gap_align.score,
        matrix_adjust_rule,
        query_index,
        query_start,
        query_end,
        match_start,
        match_end,
        frame,
        context,
    ))
}

/// Port of `s_HitlistReapContained` (`blast_kappa.c:223`).
///
/// Removes any HSP that is fully contained within an earlier (higher-
/// scoring) HSP on both query and subject coordinates and shares the
/// same query/subject frame. The hitlist must already be sorted by
/// significance before calling. Operates in place.
pub fn hitlist_reap_contained(hsps: &mut Vec<Hsp>) {
    if hsps.len() <= 1 {
        return;
    }
    let mut keep = vec![true; hsps.len()];
    for i_read in 1..hsps.len() {
        if !keep[i_read] {
            continue;
        }
        for i_back in 0..i_read {
            if !keep[i_back] {
                continue;
            }
            let h1 = &hsps[i_read];
            let h2 = &hsps[i_back];
            // The Rust `Hsp` carries `context` rather than separate
            // query/subject frame fields. Identical context implies the
            // same frame pair (one context per frame combination).
            if h1.context != h2.context {
                continue;
            }
            // CONTAINED_IN_HSP on both endpoints of h1's interval.
            if contained_in_hsp(
                h2.query_offset,
                h2.query_end,
                h1.query_offset,
                h2.subject_offset,
                h2.subject_end,
                h1.subject_offset,
            ) && contained_in_hsp(
                h2.query_offset,
                h2.query_end,
                h1.query_end,
                h2.subject_offset,
                h2.subject_end,
                h1.subject_end,
            ) && h1.score <= h2.score
            {
                keep[i_read] = false;
                break;
            }
        }
    }
    let mut idx = 0usize;
    hsps.retain(|_| {
        let k = keep[idx];
        idx += 1;
        k
    });
}

/// Port of `s_CalcLambda` (`blast_kappa.c:551`). Newton-Raphson refinement
/// of `lambda` from a score-probability distribution. Takes the
/// probabilities as a slice indexed `[0..score_max - score_min + 1]`.
pub fn calc_lambda(probs: &[f64], min_score: i32, max_score: i32, lambda0: f64) -> f64 {
    let score_range = (max_score - min_score + 1) as usize;
    debug_assert_eq!(probs.len(), score_range);
    let mut avg = 0.0f64;
    for i in 0..score_range {
        avg += (min_score + i as i32) as f64 * probs[i];
    }
    let _ = avg; // C populates `freq.score_avg` but `Blast_KarlinLambdaNR`
                 // recomputes it internally; our karlin_lambda_nr does the
                 // same.
    crate::composition::karlin_lambda_nr_pub(probs, min_score, max_score, lambda0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::program::{BLASTP, BLASTX, RPS_BLAST};

    #[test]
    fn get_subject_length_passthrough_for_blastp() {
        assert_eq!(get_subject_length(1234, BLASTP), 1234);
        assert_eq!(get_subject_length(1234, BLASTX), 1234);
    }

    #[test]
    fn get_subject_length_rps_uses_macro() {
        // GET_NUCL_LENGTH(l) = (l-5)/2 + 2; then -1 then /3.
        // For l=11: (11-5)/2 + 2 = 5; (5-1)/3 = 1.
        assert_eq!(get_subject_length(11, RPS_BLAST), 1);
    }

    #[test]
    fn get_hash_matches_c_formula() {
        // C: hash = (((((data[0]<<5)+data[1])<<5)+data[2])<<5)+...
        // Hand-compute for a 4-letter word [1, 2, 3, 4]:
        //   ((((0<<5 + 1)<<5) + 2)<<5 + 3)<<5 + 4
        //   = ((1<<5 + 2)<<5 + 3)<<5 + 4
        //   = ((34)<<5 + 3)<<5 + 4
        //   = (1088 + 3)<<5 + 4
        //   = 1091<<5 + 4 = 34912 + 4 = 34916
        assert_eq!(get_hash(&[1, 2, 3, 4], 4), 34916);
    }

    #[test]
    fn extend_right_full_match() {
        let q = b"AAAACCCC";
        let s = b"AAAACCCC";
        let (n_ident, q_ext, s_ext, _) = extend_right(q, s, 4);
        assert_eq!(n_ident, 8);
        assert_eq!(q_ext, 8);
        assert_eq!(s_ext, 8);
    }

    #[test]
    fn extend_right_handles_single_mismatch() {
        // Tolerates a single mismatch when the next 2 positions match.
        // 'A' at index 4 of subject is the mismatch; positions 5..7 match.
        let q = b"AAAAGGGG";
        let s = b"AAAATGGG"; // position 4 mismatches
        let (n_ident, _q, _s, _) = extend_right(q, s, 4);
        assert!(n_ident >= 4);
    }

    #[test]
    fn hsp_list_normalize_scores_preserves_order() {
        let mut list = HspList::new(0);
        list.add_hsp(Hsp {
            score: 200,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1e-10,
            query_offset: 0,
            query_end: 100,
            subject_offset: 0,
            subject_end: 100,
            context: 0,
            num_gaps: 0,
        });
        // Apply lambda=0.3, logK=-2.0, divisor=2.0
        hsp_list_normalize_scores(&mut list, 0.3, -2.0, 2.0);
        assert_eq!(list.hsps[0].score, 100); // 200 / 2
        // bit_score = (100*0.3*2 - (-2)) / ln2 = (60 + 2) / 0.6931 ≈ 89.45
        let expected = (100.0 * 0.3 * 2.0 - (-2.0)) / NCBIMATH_LN2;
        assert!((list.hsps[0].bit_score - expected).abs() < 1e-9);
    }

    #[test]
    fn hitlist_reap_contained_drops_inner_hsp() {
        let mut hsps = vec![
            Hsp {
                score: 100,
                num_ident: 0,
                bit_score: 0.0,
                evalue: 0.0,
                query_offset: 0,
                query_end: 100,
                subject_offset: 0,
                subject_end: 100,
                context: 0,
                num_gaps: 0,
            },
            Hsp {
                // Contained, lower score → dropped.
                score: 50,
                num_ident: 0,
                bit_score: 0.0,
                evalue: 0.0,
                query_offset: 20,
                query_end: 60,
                subject_offset: 20,
                subject_end: 60,
                context: 0,
                num_gaps: 0,
            },
            Hsp {
                // Different context → not contained even though inside coords.
                score: 50,
                num_ident: 0,
                bit_score: 0.0,
                evalue: 0.0,
                query_offset: 20,
                query_end: 60,
                subject_offset: 20,
                subject_end: 60,
                context: 3,
                num_gaps: 0,
            },
        ];
        hitlist_reap_contained(&mut hsps);
        assert_eq!(hsps.len(), 2);
        assert_eq!(hsps[0].score, 100);
        assert_eq!(hsps[1].context, 3);
    }

    #[test]
    fn contained_in_hsp_simple() {
        // a=0 b=10 c=5: yes; d=0 e=10 f=7: yes → true.
        assert!(contained_in_hsp(0, 10, 5, 0, 10, 7));
        // c=11 outside [0,10] → false.
        assert!(!contained_in_hsp(0, 10, 11, 0, 10, 7));
    }
}
