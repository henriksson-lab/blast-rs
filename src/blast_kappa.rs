//! Rust port of `blast_kappa.c` — Smith-Waterman / composition-based
//! scoring for the traceback stage of blastp / blastx / tblastn / RPS.
//!
//! Most leaf utilities and redo-alignment support helpers are ported here.
//! The remaining gaps are narrower integration points around full program
//! orchestration and fixture-level parity rather than the core redo helpers.
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
//! | comp adj | `s_GetPosBasedStartFreqRatios` | 24 | ✅ (BLOSUM62 family) |
//! | comp adj | `s_GetStartFreqRatios` | 12 | ✅ |
//! | comp adj | `s_ScalePosMatrix` | 44 | ✅ (position rows + scaled PSI scores + Lambda-ratio stat hook) |
//! | redo adj | `s_NewAlignmentFromGapAlign` | 20 | ✅ |
//! | redo adj | `s_NewAlignmentUsingXdrop` | 44 | ✅ |
//! | redo adj | `s_HSPListFromDistinctAlignments` | 33 | ✅ |
//! | redo adj | `s_RedoOneAlignment` | 33 | ✅ |
//! | redo adj | `Blast_RedoOneMatch` | 200 | ✅ callback + in-memory ordinary composition paths |
//! | redo adj | SW significance gate in `Blast_RedoOneMatchSmithWaterman` | 16 | ✅ |
//! | redo adj | `Blast_RedoOneMatchSmithWaterman` | 241 | ✅ callback + materialized ordinary/SW composition paths; BLASTN composition SW rejected |
//! | redo adj | `s_TestNearIdentical` | 55 | ✅ |
//! | redo adj | `s_HitlistEvaluateAndPurge` | 36 | ✅ (link-HSPs when context supplied) |
//! | redo adj | `s_ComputeNumIdentities` | 60 | ✅ (protein / in-memory translated subject) |
//! | redo adj | `s_AdjustEvaluesForComposition` | 75 | ✅ |
//! | redo adj | `s_ResultHspToDistinctAlign` | 25 | ✅ |
//! | redo adj | `s_WithDistinctEnds` | 86 | ✅ |
//! | redo adj | `s_WindowsFromProteinAligns` | 63 | ✅ |
//! | redo adj | `s_WindowsFromTranslatedAligns` | 122 | ✅ |
//! | redo adj | `s_WindowsFromAligns` | 31 | ✅ |
//! | redo adj | `s_DistinctAlignmentsSort` | 58 | ✅ |
//! | redo adj | `s_GetComposition` | 39 | ✅ |
//! | redo adj | `s_EvalueFromScore` | 4 | ✅ |
//! | redo adj | `s_preliminaryTestNearIdentical` | 42 | ✅ |
//! | redo adj | `BlastCompo_EarlyTermination` | 21 | ✅ |
//! | redo adj | `BlastCompo_AlignmentsFree` | 16 | ✅ |
//! | redo adj | `s_SequenceDataRelease` | 6 | ✅ |
//! | redo adj | `s_WindowInfoFree` | 8 | ✅ |
//! | redo adj | `s_SWFindFinalEndsUsingXdrop` | 30 | ✅ |
//! | redo adj | `s_MatchingSequenceRelease` | 10 | ✅ |
//! | redo adj | `s_MatchingSequenceInitialize` | 33 | ✅ |
//! | redo adj | `s_DoSegSequenceData` | 22 | ✅ |
//! | redo adj | `s_SequenceGetTranslatedRange` | 51 | ✅ (in-memory NCBI4na source) |
//! | redo adj | `s_SequenceGetProteinRange` | 52 | ✅ (byte-slice source) |
//! | redo adj | `s_SequenceGetRange` | 31 | ✅ (in-memory source) |
//! | redo adj | `s_SavedParametersFree` | 13 | ✅ |
//! | redo adj | `s_SavedParametersNew` | 26 | ✅ |
//! | redo adj | `s_RecordInitialSearch` | 27 | ✅ |
//! | redo adj | `s_RestoreSearch` | 23 | ✅ |
//! | redo adj | `s_MatrixInfoInit` | 28 | ✅ non-position + PSI setup, including private scaled rows and Lambda-ratio stat update hook |
//! | redo adj | `s_CreateWordArray` | 16 | ✅ |
//! | redo adj | `s_FreeBlastCompo_QueryInfoArray` | 10 | ✅ |
//! | redo adj | `s_GetQueryInfo` | 22 | ✅ |
//! | redo adj | `s_GappingParamsNew` | 23 | ✅ |
//! | redo adj | `s_GetAlignParams` | 48 | ✅ |
//! | redo adj | `s_FillResultsFromCompoHeaps` | 13 | ✅ |
//! | redo adj | `s_ClearHeap` | 3 | ✅ |
//! | redo adj | `s_BlastGapAlignStruct_Free` | 29 | ✅ |
//! | redo adj | `s_BlastGapAlignStruct_Copy` | 95 | ✅ |
//! | redo adj | `s_BlastScoreBlk_Copy` | 131 | ✅ (modeled fields) |
//! | redo adj | `Blast_RedoAlignCallbacks` | header struct | ✅ |
//! | driver | `Blast_RedoAlignmentCore` | 31 | ✅ single-thread translation |
//! | driver | `Blast_RedoAlignmentCore_MT` | 669 | ✅ materialized/stream/SeqSrc/callback redo, heap merge, sum-stat, translated, PSSM, and SW branches represented |

use crate::compo_mode_condition::MatrixAdjustRule;
use crate::hspstream::{Hsp, HspList};
use crate::math::NCBIMATH_LN2;
use crate::program::{ProgramType, RPS_TBLASTN};

/// `kWindowBorder` (`redo_alignment.c:112`).
pub const K_WINDOW_BORDER: i32 = 200;

// ───────────────────────────────────────────────────────────────────────────
// Struct ports from `composition_adjustment/redo_alignment.h` and
// `composition_adjustment/composition_constants.h`.
//
// These are 1-1 ports of the C structs used by `Blast_RedoAlignmentCore_MT`
// and its helpers. Function-pointer fields in C (`Blast_RedoAlignCallbacks`)
// are represented as Rust function-pointer types so the dispatch table
// translates directly. Executable callback implementations still live in the
// concrete in-memory helpers until the full external callback context is
// available to the driver.
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
    /// blast-rs: Numeric enum conversion helper; not a direct NCBI C port.
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
/// `context` in NCBI's C is `void *` carrying traceback data. The faithful
/// Rust slot is `BlastCompoAlignmentContext`: it can carry the copied
/// traceback script, the original HSP context used by `s_RedoOneAlignment`,
/// or both while the surrounding redo code is being made C-shaped.
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
    pub context: Option<BlastCompoAlignmentContext>,
    pub next: Option<Box<BlastCompoAlignment>>,
}

/// Rust owner for C `BlastCompo_Alignment.context`.
///
/// C stores an untyped pointer here. In the redo path it is either traceback
/// data that must be freed by the caller-provided deleter, or the original
/// `BlastHSP*` whose gapped starts seed `s_RedoOneAlignment`. Keeping those
/// members here lets the alignment struct keep the upstream field set instead
/// of carrying Rust-only top-level gapped-start fields.
#[derive(Debug, Clone, Default)]
pub struct BlastCompoAlignmentContext {
    pub edit_script: Option<crate::gapinfo::GapEditScript>,
    pub hsp: Option<BlastCompoAlignmentHspContext>,
}

/// Minimal owned mirror of the `BlastHSP*` data read through
/// `BlastCompo_Alignment.context` by NCBI `s_RedoOneAlignment`.
#[derive(Debug, Clone, Copy, Default)]
pub struct BlastCompoAlignmentHspContext {
    pub query_gapped_start: i32,
    pub subject_gapped_start: i32,
}

impl BlastCompoAlignment {
    /// NCBI: BlastCompo_AlignmentNew (redo_alignment.c:140).
    ///
    /// The C version `malloc`s and
    /// initializes; Rust returns a stack value the caller can wrap.
    /// naming: Rust exposes this as the associated constructor for
    /// `BlastCompoAlignment`.
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
        context: Option<BlastCompoAlignmentContext>,
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

    pub fn with_hsp_context(mut self, query_gapped_start: i32, subject_gapped_start: i32) -> Self {
        let ctx = self.context.get_or_insert_with(Default::default);
        ctx.hsp = Some(BlastCompoAlignmentHspContext {
            query_gapped_start,
            subject_gapped_start,
        });
        self
    }

    pub fn query_gapped_start(&self) -> i32 {
        self.context
            .as_ref()
            .and_then(|ctx| ctx.hsp)
            .map(|hsp| hsp.query_gapped_start)
            .unwrap_or(self.query_start)
    }

    pub fn match_gapped_start(&self) -> i32 {
        self.context
            .as_ref()
            .and_then(|ctx| ctx.hsp)
            .map(|hsp| hsp.subject_gapped_start)
            .unwrap_or(self.match_start)
    }
}

/// `BlastCompo_GappingParams` (`redo_alignment.h:101`). Parameters
/// passed to the gapped-alignment callback inside the kappa driver.
#[derive(Debug, Clone)]
pub struct BlastCompoGappingParams {
    pub gap_open: i32,
    pub gap_extend: i32,
    pub decline_align: i32,
    pub x_dropoff: i32,
    pub context: Option<BlastCompoGappingContext>,
}

impl Default for BlastCompoGappingParams {
    fn default() -> Self {
        Self {
            gap_open: 0,
            gap_extend: 0,
            decline_align: 0,
            x_dropoff: 0,
            context: None,
        }
    }
}

/// Rust owner for C `BlastCompo_GappingParams.context`.
#[derive(Debug, Clone, Default)]
pub struct BlastCompoGappingContext {
    pub oof_context_index: Option<i32>,
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

/// `s_WindowInfo` (`redo_alignment.c:118`). A subject/query range pair plus
/// the list of alignments assigned to that window.
#[derive(Debug, Clone)]
pub struct BlastCompoWindowInfo {
    pub query_range: BlastCompoSequenceRange,
    pub subject_range: BlastCompoSequenceRange,
    pub align: Option<Box<BlastCompoAlignment>>,
    pub hspcnt: i32,
}

impl BlastCompoWindowInfo {
    /// NCBI: s_WindowInfoNew (redo_alignment.c:527).
    /// naming: Associated constructor on `BlastCompoWindowInfo`.
    pub fn new(
        begin: i32,
        end: i32,
        context: i32,
        query_origin: i32,
        query_length: i32,
        query_index: i32,
        align: Option<Box<BlastCompoAlignment>>,
    ) -> Self {
        let hspcnt = alignment_list_len(align.as_deref()) as i32;
        Self {
            query_range: BlastCompoSequenceRange {
                begin: query_origin,
                end: query_origin + query_length,
                context: query_index,
            },
            subject_range: BlastCompoSequenceRange {
                begin,
                end,
                context,
            },
            align,
            hspcnt,
        }
    }

    /// NCBI: s_WindowSwapRange (redo_alignment.c:566).
    /// naming: Associated method on `BlastCompoWindowInfo`; type supplies `window_info`.
    pub fn swap_range(&mut self) {
        std::mem::swap(&mut self.query_range, &mut self.subject_range);
    }

    /// NCBI: s_WindowInfoJoin (redo_alignment.c:591).
    /// naming: Associated method on `BlastCompoWindowInfo`; type supplies `window_info`.
    ///
    /// Appends `other.align` to the tail of `self.align`, expands the subject
    /// range, and adds HSP counts. The caller is responsible for ensuring the
    /// subject/query contexts match, as NCBI asserts.
    pub fn join(&mut self, mut other: BlastCompoWindowInfo) {
        assert_eq!(self.subject_range.context, other.subject_range.context);
        assert_eq!(self.query_range.context, other.query_range.context);

        self.subject_range.begin = self.subject_range.begin.min(other.subject_range.begin);
        self.subject_range.end = self.subject_range.end.max(other.subject_range.end);
        self.hspcnt += other.hspcnt;

        if self.align.is_none() {
            self.align = other.align.take();
            return;
        }
        let mut cursor = self.align.as_mut().unwrap();
        while cursor.next.is_some() {
            cursor = cursor.next.as_mut().unwrap();
        }
        cursor.next = other.align.take();
    }
}

/// NCBI: s_WindowInfoFree (redo_alignment.c:552).
/// naming: Rust drops owned values, so the public parity hook omits the `s_` prefix.
///
/// Rust frees the window and its alignment list automatically. This parity
/// hook clears the caller's slot to mirror C's `*window = NULL`.
pub fn window_info_free(window: &mut Option<BlastCompoWindowInfo>) {
    if let Some(window) = window.as_mut() {
        alignments_free(&mut window.align);
    }
    *window = None;
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
    /// blast-rs: Slice accessor for Rust's sentinel-backed sequence storage; not a direct NCBI C port.
    ///
    /// Returns the residue slice (`&data[0..length]` in C terms,
    /// equivalent to `&buffer[data_offset .. data_offset + length]`).
    pub fn data(&self) -> &[u8] {
        let start = self.data_offset;
        let end = start + self.length as usize;
        &self.buffer[start..end]
    }
}

/// NCBI: s_SequenceDataRelease (redo_alignment.c:486).
/// naming: Public Rust parity hook omits the private `s_` prefix.
///
/// In C this frees `buffer`, then nulls both `data` and `buffer`. In Rust the
/// buffer owns the data; clearing it and resetting the offset/length gives the
/// same observable postcondition for later parity call sites.
pub fn sequence_data_release(self_: &mut BlastCompoSequenceData) {
    self_.buffer.clear();
    self_.data_offset = 0;
    self_.length = 0;
}

/// `BlastCompo_MatchingSequence` (`redo_alignment.h:156`). Identifies
/// one subject sequence and provides a hook for callbacks to access
/// the underlying database/translation.
#[derive(Debug, Clone)]
pub struct BlastCompoMatchingSequence {
    pub length: i32,
    pub index: i32,
    pub local_data: Option<BlastKappaSequenceInfo>,
}

/// `BlastKappa_SequenceInfo` (`blast_kappa.c:893`), stored through
/// `BlastCompo_MatchingSequence.local_data`.
#[derive(Clone)]
pub struct BlastKappaSequenceInfo {
    pub prog_number: ProgramType,
    pub seq_src: Option<*const dyn crate::seqsrc::BlastSeqSource>,
    pub seq_arg: crate::seqsrc::BlastSeqSrcGetSeqArg,
}

impl std::fmt::Debug for BlastKappaSequenceInfo {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("BlastKappaSequenceInfo")
            .field("prog_number", &self.prog_number)
            .field("seq_src", &self.seq_src.map(|_| "BlastSeqSource"))
            .field("seq_arg", &self.seq_arg)
            .finish()
    }
}

impl Default for BlastCompoMatchingSequence {
    fn default() -> Self {
        Self {
            length: 0,
            index: -1,
            local_data: None,
        }
    }
}

/// blast-rs: Materialized-subject convenience initializer; not a direct NCBI C port.
pub fn matching_sequence_initialize(
    length: i32,
    index: i32,
    local_data_index: i32,
) -> BlastCompoMatchingSequence {
    let mut matching_sequence = BlastCompoMatchingSequence {
        length: 0,
        index: 0,
        local_data: None,
    };
    matching_sequence.length = length;
    matching_sequence.index = index;
    matching_sequence.local_data = Some(BlastKappaSequenceInfo {
        prog_number: crate::program::UNDEFINED,
        seq_src: None,
        seq_arg: crate::seqsrc::BlastSeqSrcGetSeqArg {
            oid: local_data_index,
            check_oid_exclusion: true,
            ..crate::seqsrc::BlastSeqSrcGetSeqArg::default()
        },
    });
    matching_sequence
}

/// NCBI: `s_MatchingSequenceInitialize` (`blast_kappa.c:1357`).
#[allow(clippy::too_many_arguments)]
pub fn s_matching_sequence_initialize(
    self_: &mut BlastCompoMatchingSequence,
    program_number: ProgramType,
    seq_src: &dyn crate::seqsrc::BlastSeqSource,
    default_db_genetic_code: i32,
    subject_index: i32,
    ranges: Option<crate::seqsrc::BlastSeqSrcSetRangesArg>,
) -> i32 {
    self_.length = 0;
    self_.local_data = None;

    let mut seq_arg = crate::seqsrc::BlastSeqSrcGetSeqArg {
        oid: subject_index,
        check_oid_exclusion: true,
        ranges,
        encoding: if program_number == crate::program::TBLASTN {
            crate::seqsrc::EBlastEncoding::Ncbi4na
        } else {
            crate::seqsrc::EBlastEncoding::Protein
        },
        ..crate::seqsrc::BlastSeqSrcGetSeqArg::default()
    };
    self_.index = subject_index;

    if crate::seqsrc::blast_seq_src_get_sequence(Some(seq_src), Some(&mut seq_arg)).is_some() {
        self_.length = seq_arg.seq.as_ref().map(|seq| seq.length).unwrap_or(0);
        if crate::program::blast_subject_is_translated(program_number) {
            if let Some(seq) = seq_arg.seq.as_mut() {
                if seq.gen_code_string.is_none() {
                    let gen_code =
                        crate::util::gen_code_singleton_find(default_db_genetic_code as u32)
                            .and_then(|code| code.try_into().ok())
                            .unwrap_or_else(|| {
                                *crate::util::lookup_genetic_code(default_db_genetic_code as u8)
                            });
                    seq.gen_code_string = Some(gen_code);
                }
            }
        }
    }

    let seq_src_ptr = seq_src as *const dyn crate::seqsrc::BlastSeqSource;
    // C stores this as a borrowed BlastSeqSrc*. Rust's trait-object raw pointer
    // carries a lifetime, so erase it at the same boundary; callers must keep
    // the source alive until `s_MatchingSequenceRelease`.
    let seq_src_ptr: *const (dyn crate::seqsrc::BlastSeqSource + 'static) = unsafe {
        std::mem::transmute::<
            *const dyn crate::seqsrc::BlastSeqSource,
            *const (dyn crate::seqsrc::BlastSeqSource + 'static),
        >(seq_src_ptr)
    };

    self_.local_data = Some(BlastKappaSequenceInfo {
        prog_number: program_number,
        seq_src: Some(seq_src_ptr),
        seq_arg,
    });

    if self_.length == 0 {
        matching_sequence_release(self_);
        return -1;
    }

    0
}

/// `Blast_AminoAcidComposition` (`composition_adjustment.h:54`).
#[derive(Debug, Clone)]
pub struct BlastAminoAcidComposition {
    pub prob: Vec<f64>,
    pub num_true_amino_acids: i32,
}

impl Default for BlastAminoAcidComposition {
    fn default() -> Self {
        Self {
            prob: vec![0.0; crate::matrix::AA_SIZE],
            num_true_amino_acids: 0,
        }
    }
}

impl BlastAminoAcidComposition {
    pub fn from_prob(prob: Vec<f64>, num_true_amino_acids: i32) -> Self {
        Self {
            prob,
            num_true_amino_acids,
        }
    }
}

/// `BlastCompo_QueryInfo` (`redo_alignment.h:166`). Per-query metadata
/// consumed by the redo-alignment driver. `composition` (the
/// `Blast_AminoAcidComposition` struct) and `words` (Uint8 array of
/// hashed query k-mers) are filled in by `s_GetQueryInfo` /
/// `s_CreateWordArray`.
#[derive(Debug, Clone, Default)]
pub struct BlastCompoQueryInfo {
    pub origin: i32,
    pub seq: BlastCompoSequenceData,
    pub composition: BlastAminoAcidComposition,
    pub eff_search_space: f64,
    pub words: Vec<u64>,
}

/// NCBI: Blast_MatrixInfo (composition_adjustment.h).
///
/// Holds the
/// scoring matrix used for the redo-alignment pass, plus its
/// dimensions, name, and the rounding/scale parameters used by
/// composition adjustment.
#[derive(Debug, Clone, Default)]
pub struct BlastMatrixInfo {
    pub matrix_name: String,
    pub rows: i32,
    pub cols: i32,
    pub positional: bool,
    pub bit_scale_factor: i32,
    pub ungapped_lambda: f64,
    pub matrix: Vec<Vec<i32>>,
    pub scaled_matrix: Vec<Vec<i32>>,
    pub start_freq_ratios: Vec<Vec<f64>>,
}

/// `kPSIScaleFactor` (`blast_psi_priv.c:59`).
pub const K_PSI_SCALE_FACTOR: f64 = 200.0;
/// `kPositScalingPercent` (`blast_psi_priv.c:60`).
pub const K_POSIT_SCALING_PERCENT: f64 = 0.05;
/// `kPositScalingNumIterations` (`blast_psi_priv.c:61`).
pub const K_POSIT_SCALING_NUM_ITERATIONS: usize = 10;
/// `kScoreMatrixScoreRange` (`blast_posit.h:49`).
pub const K_SCORE_MATRIX_SCORE_RANGE: i32 = 10_000;

/// Rust-owned counterpart of NCBI `Kappa_posSearchItems`.
#[derive(Debug, Clone)]
pub struct KappaPosSearchItems {
    pub pos_matrix: Vec<Vec<i32>>,
    pub std_freq_ratios: Vec<Vec<f64>>,
    pub query_length: usize,
    pub pos_private_matrix: Vec<Vec<i32>>,
    pub pos_freqs: Vec<Vec<f64>>,
}

/// Rust-owned counterpart of NCBI `Kappa_compactSearchItems`.
#[derive(Debug, Clone)]
pub struct KappaCompactSearchItems {
    pub standard_prob: Vec<f64>,
    pub query: Vec<u8>,
    pub qlength: usize,
    pub alphabet_size: usize,
    pub matrix: Vec<Vec<i32>>,
    pub kbp_std: Vec<crate::stat::KarlinBlk>,
    pub kbp_psi: Vec<crate::stat::KarlinBlk>,
    pub kbp_gap_std: Vec<crate::stat::KarlinBlk>,
    pub kbp_gap_psi: Vec<crate::stat::KarlinBlk>,
    pub lambda_ideal: f64,
    pub k_ideal: f64,
}

/// Rust-owned port boundary for `Blast_CompositionWorkspace`
/// (`composition_adjustment.c:1285`).
///
/// NCBI stores `first_standard_freq`, `second_standard_freq`, and `mat_b`
/// inside this workspace after `Blast_CompositionWorkspaceInit`. The local
/// composition module currently exposes the BLOSUM62 initialization data, so
/// this struct carries those same fields for relative-entropy adjustment.
#[derive(Debug, Clone)]
pub struct BlastCompositionWorkspace {
    pub mat_b:
        [[f64; crate::composition::COMPO_NUM_TRUE_AA]; crate::composition::COMPO_NUM_TRUE_AA],
    pub mat_final:
        [[f64; crate::composition::COMPO_NUM_TRUE_AA]; crate::composition::COMPO_NUM_TRUE_AA],
    pub first_standard_freq: [f64; crate::composition::COMPO_NUM_TRUE_AA],
    pub second_standard_freq: [f64; crate::composition::COMPO_NUM_TRUE_AA],
}

impl BlastCompositionWorkspace {
    /// blast-rs: Constructs the local BLOSUM62 composition workspace; not a direct NCBI C port.
    pub fn blosum62() -> Self {
        let (mat_b, first_standard_freq, second_standard_freq) =
            crate::composition::blosum62_workspace();
        Self {
            mat_b,
            mat_final: [[0.0; crate::composition::COMPO_NUM_TRUE_AA];
                crate::composition::COMPO_NUM_TRUE_AA],
            first_standard_freq,
            second_standard_freq,
        }
    }
}

/// Rust-owned subset of NCBI `BlastScoreBlk` used by the kappa/composition
/// pipeline.
///
/// The C `s_BlastScoreBlk_Copy` deep-copies many optional fields from
/// `BlastScoreBlk` before composition adjustment mutates the live score
/// block. This crate does not model every C-side score-block member, so this
/// struct carries the fields that are represented locally and can be restored:
/// matrix data, PSI matrix data when present, per-context Karlin blocks, the
/// active matrix name, and scale/rounding metadata.
#[derive(Debug, Clone, Default)]
pub struct BlastScoreBlkSnapshot {
    pub matrix_name: String,
    pub matrix: Vec<Vec<i32>>,
    pub psi_matrix: Option<Vec<Vec<i32>>>,
    pub kbp: Vec<crate::stat::KarlinBlk>,
    pub kbp_gap: Vec<crate::stat::KarlinBlk>,
    pub scale_factor: f64,
    pub round_down: bool,
}

/// NCBI: s_BlastScoreBlk_Copy (blast_kappa.c:2609).
/// naming: Rust snapshot helper omits the private `s_` prefix and models only local fields.
///
/// Rust's `Clone` performs the deep copies that NCBI does field-by-field with
/// `Blast_KarlinBlkCopy`, matrix allocation/copy, and PSI-matrix row copies.
/// Keeping this translation makes call sites mirror the C flow and protects the
/// important semantic: mutating the returned snapshot must not alias the source.
pub fn blast_score_blk_copy(src: &BlastScoreBlkSnapshot) -> BlastScoreBlkSnapshot {
    src.clone()
}

/// Rust representation of `calc_lambda_type` (`redo_alignment.h:220`).
pub type BlastCalcLambdaFn = fn(probs: &[f64], min_score: i32, max_score: i32, lambda0: f64) -> f64;

/// Rust representation of `get_range_type` (`redo_alignment.h:247`).
///
/// NCBI's callback receives C pointers to sequence buffers plus the caller's
/// external context through `BlastCompo_MatchingSequence.local_data`. Rust
/// keeps the context outside this function pointer and passes owned structs and
/// slices directly.
pub type BlastGetRangeFn = fn(
    sequence: &BlastCompoMatchingSequence,
    range: &BlastCompoSequenceRange,
    data: &mut BlastCompoSequenceData,
    orig_query: &BlastCompoSequenceData,
    q_range: &BlastCompoSequenceRange,
    q_data: &mut BlastCompoSequenceData,
    query_words: &[u64],
    align: &BlastCompoAlignment,
    should_test_identical: bool,
    compo_adjust_mode: CompoAdjustMode,
    is_smith_waterman: bool,
    subject_maybe_biased: &mut bool,
) -> i32;

/// Rust representation of `redo_one_alignment_type`
/// (`redo_alignment.h:286`).
pub type BlastRedoOneAlignmentFn = fn(
    in_align: &BlastCompoAlignment,
    matrix_adjust_rule: MatrixAdjustRule,
    query_data: &BlastCompoSequenceData,
    query_range: &BlastCompoSequenceRange,
    ccat_query_length: i32,
    subject_data: &BlastCompoSequenceData,
    subject_range: &BlastCompoSequenceRange,
    full_subject_length: i32,
    gapping_params: &BlastCompoGappingParams,
) -> Option<BlastCompoAlignment>;

/// Rust representation of `new_xdrop_align_type`
/// (`redo_alignment.h:321`).
pub type BlastNewXdropAlignFn = fn(
    align: &mut Option<BlastCompoAlignment>,
    query_end: &mut i32,
    match_end: &mut i32,
    query_start: i32,
    match_start: i32,
    score: i32,
    query: &BlastCompoSequenceData,
    query_range: &BlastCompoSequenceRange,
    ccat_query_length: i32,
    subject: &BlastCompoSequenceData,
    subject_range: &BlastCompoSequenceRange,
    full_subject_length: i32,
    gapping_params: &BlastCompoGappingParams,
    matrix_adjust_rule: MatrixAdjustRule,
) -> i32;

/// Rust representation of `free_align_traceback_type`
/// (`redo_alignment.h:339`).
pub type BlastFreeAlignTracebackFn = fn(traceback_data: &mut Option<crate::gapinfo::GapEditScript>);

/// `Blast_RedoAlignCallbacks` (`redo_alignment.h:344`).
#[derive(Debug, Clone, Copy, Default)]
pub struct BlastRedoAlignCallbacks {
    pub calc_lambda: Option<BlastCalcLambdaFn>,
    pub get_range: Option<BlastGetRangeFn>,
    pub redo_one_alignment: Option<BlastRedoOneAlignmentFn>,
    pub new_xdrop_align: Option<BlastNewXdropAlignFn>,
    pub free_align_traceback: Option<BlastFreeAlignTracebackFn>,
}

/// `Blast_RedoAlignParams` (`redo_alignment.h:360`). Parameter block
/// owned by the kappa driver.
#[derive(Debug, Clone)]
pub struct BlastRedoAlignParams {
    pub matrix_info: BlastMatrixInfo,
    pub gapping_params: BlastCompoGappingParams,
    pub compo_adjust_mode: CompoAdjustMode,
    pub local_scaling_factor: f64,
    pub position_based: bool,
    pub re_pseudocounts: i32,
    pub subject_is_translated: bool,
    pub query_is_translated: bool,
    pub ccat_query_length: i32,
    pub cutoff_s: i32,
    pub cutoff_e: f64,
    pub do_link_hsps: bool,
    pub callbacks: Option<BlastRedoAlignCallbacks>,
    pub near_identical_cutoff: f64,
}

#[cfg(test)]
mod struct_tests {
    use super::*;

    #[test]
    fn alignment_new_round_trip() {
        let a =
            BlastCompoAlignment::new(42, MatrixAdjustRule::DontAdjust, 0, 10, 50, 20, 60, 1, None);
        assert_eq!(a.score, 42);
        assert_eq!(a.matrix_adjust_rule, MatrixAdjustRule::DontAdjust);
        assert_eq!(a.query_start, 10);
        assert_eq!(a.query_end, 50);
        assert_eq!(a.match_start, 20);
        assert_eq!(a.match_end, 60);
        assert!(a.next.is_none());
    }

    #[test]
    fn alignments_free_clears_entire_list() {
        let mut head = Some(Box::new(BlastCompoAlignment::new(
            42,
            MatrixAdjustRule::DontAdjust,
            0,
            10,
            50,
            20,
            60,
            1,
            None,
        )));
        head.as_mut().unwrap().next = Some(Box::new(BlastCompoAlignment::new(
            41,
            MatrixAdjustRule::ScaleOldMatrix,
            0,
            11,
            49,
            21,
            59,
            1,
            None,
        )));

        alignments_free(&mut head);
        assert!(head.is_none());
    }

    #[test]
    fn alignments_free_array_clears_prefix() {
        let mut alignments = vec![
            Some(Box::new(BlastCompoAlignment::new(
                42,
                MatrixAdjustRule::DontAdjust,
                0,
                10,
                50,
                20,
                60,
                1,
                None,
            ))),
            Some(Box::new(BlastCompoAlignment::new(
                41,
                MatrixAdjustRule::ScaleOldMatrix,
                1,
                11,
                49,
                21,
                59,
                1,
                None,
            ))),
        ];

        alignments_free_array(&mut alignments, 1);

        assert!(alignments[0].is_none());
        assert!(alignments[1].is_some());
    }

    #[test]
    fn free_edit_script_name_matched_wrapper_clears_owned_script() {
        let mut script = crate::gapinfo::GapEditScript::new();
        script.push(crate::gapinfo::GapAlignOpType::Sub, 3);

        assert!(s_free_edit_script(Some(script)).is_none());
        assert!(s_free_edit_script(None).is_none());
    }

    fn test_calc_lambda(_: &[f64], _: i32, _: i32, lambda0: f64) -> f64 {
        lambda0
    }

    #[test]
    fn redo_align_params_new_with_callbacks_preserves_callback_table() {
        let callbacks = BlastRedoAlignCallbacks {
            calc_lambda: Some(test_calc_lambda),
            ..Default::default()
        };

        let params = blast_redo_align_params_new_with_callbacks(
            BlastMatrixInfo::default(),
            BlastCompoGappingParams {
                gap_open: 5,
                gap_extend: 2,
                decline_align: 0,
                x_dropoff: 20,
                context: None,
            },
            CompoAdjustMode::NoCompositionBasedStats,
            1.0,
            false,
            false,
            false,
            10,
            1,
            10.0,
            false,
            Some(callbacks),
            0.0,
        );

        let stored = params.callbacks.expect("callbacks stored");
        assert!(stored.calc_lambda.is_some());
        assert!(stored.get_range.is_none());
        assert_eq!(params.re_pseudocounts, K_RE_MATRIX_ADJUSTMENT_PSEUDOCOUNTS);
    }

    #[test]
    fn window_info_free_clears_slot_and_alignment() {
        let align = Some(Box::new(BlastCompoAlignment::new(
            42,
            MatrixAdjustRule::DontAdjust,
            0,
            10,
            50,
            20,
            60,
            1,
            None,
        )));
        let mut window = Some(BlastCompoWindowInfo::new(20, 60, 1, 0, 100, 0, align));

        window_info_free(&mut window);
        assert!(window.is_none());
    }

    #[test]
    fn compo_adjust_mode_from_u8() {
        assert_eq!(
            CompoAdjustMode::from_u8(0),
            CompoAdjustMode::NoCompositionBasedStats
        );
        assert_eq!(
            CompoAdjustMode::from_u8(1),
            CompoAdjustMode::CompositionBasedStats
        );
        assert_eq!(
            CompoAdjustMode::from_u8(2),
            CompoAdjustMode::CompositionMatrixAdjust
        );
        assert_eq!(
            CompoAdjustMode::from_u8(3),
            CompoAdjustMode::CompoForceFullMatrixAdjust
        );
        // Anything > 3 saturates to the strongest mode (matches C's
        // raw enum cast — out-of-range becomes max).
        assert_eq!(
            CompoAdjustMode::from_u8(99),
            CompoAdjustMode::CompoForceFullMatrixAdjust
        );
    }

    #[test]
    fn redo_one_alignment_perfect_blastn_match() {
        // BLASTNA-encoded perfect-match query/subject. Our
        // `blast_gapped_alignment_with_traceback` operates on BLASTNA, so set up a 16x16
        // matrix via blast_score_blk_nucl_matrix_create.
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
        let qr = BlastCompoSequenceRange {
            begin: 0,
            end: 10,
            context: 0,
        };
        let sr = BlastCompoSequenceRange {
            begin: 0,
            end: 10,
            context: 0,
        };
        let result = redo_one_alignment(
            &query,
            &qr,
            &subject,
            &sr,
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

    fn redo_one_match_test_params(compo_adjust_mode: CompoAdjustMode) -> BlastRedoAlignParams {
        BlastRedoAlignParams {
            matrix_info: BlastMatrixInfo::default(),
            gapping_params: BlastCompoGappingParams {
                gap_open: 5,
                gap_extend: 2,
                decline_align: 0,
                x_dropoff: 30,
                context: None,
            },
            compo_adjust_mode,
            local_scaling_factor: 1.0,
            position_based: false,
            re_pseudocounts: K_RE_MATRIX_ADJUSTMENT_PSEUDOCOUNTS,
            subject_is_translated: false,
            query_is_translated: false,
            ccat_query_length: 10,
            cutoff_s: 1,
            cutoff_e: 10.0,
            do_link_hsps: false,
            callbacks: None,
            near_identical_cutoff: 0.0,
        }
    }

    fn callback_test_get_range(
        _sequence: &BlastCompoMatchingSequence,
        range: &BlastCompoSequenceRange,
        data: &mut BlastCompoSequenceData,
        orig_query: &BlastCompoSequenceData,
        _q_range: &BlastCompoSequenceRange,
        q_data: &mut BlastCompoSequenceData,
        _query_words: &[u64],
        _align: &BlastCompoAlignment,
        _should_test_identical: bool,
        _compo_adjust_mode: CompoAdjustMode,
        _is_smith_waterman: bool,
        subject_maybe_biased: &mut bool,
    ) -> i32 {
        let length = (range.end - range.begin).max(0) as usize;
        data.buffer = vec![0; length];
        data.data_offset = 0;
        data.length = length as i32;
        *q_data = orig_query.clone();
        *subject_maybe_biased = false;
        0
    }

    fn callback_test_get_range_alanine(
        _sequence: &BlastCompoMatchingSequence,
        range: &BlastCompoSequenceRange,
        data: &mut BlastCompoSequenceData,
        orig_query: &BlastCompoSequenceData,
        _q_range: &BlastCompoSequenceRange,
        q_data: &mut BlastCompoSequenceData,
        _query_words: &[u64],
        _align: &BlastCompoAlignment,
        _should_test_identical: bool,
        _compo_adjust_mode: CompoAdjustMode,
        _is_smith_waterman: bool,
        subject_maybe_biased: &mut bool,
    ) -> i32 {
        let length = (range.end - range.begin).max(0) as usize;
        data.buffer = vec![1; length];
        data.data_offset = 0;
        data.length = length as i32;
        *q_data = orig_query.clone();
        *subject_maybe_biased = false;
        0
    }

    fn callback_test_redo_one_alignment(
        in_align: &BlastCompoAlignment,
        matrix_adjust_rule: MatrixAdjustRule,
        _query_data: &BlastCompoSequenceData,
        _query_range: &BlastCompoSequenceRange,
        _ccat_query_length: i32,
        _subject_data: &BlastCompoSequenceData,
        _subject_range: &BlastCompoSequenceRange,
        _full_subject_length: i32,
        _gapping_params: &BlastCompoGappingParams,
    ) -> Option<BlastCompoAlignment> {
        Some(BlastCompoAlignment::new(
            in_align.score + 1,
            matrix_adjust_rule,
            in_align.query_index,
            in_align.query_start,
            in_align.query_end,
            in_align.match_start,
            in_align.match_end,
            in_align.frame,
            None,
        ))
    }

    fn callback_test_new_xdrop_align(
        align: &mut Option<BlastCompoAlignment>,
        query_end: &mut i32,
        match_end: &mut i32,
        query_start: i32,
        match_start: i32,
        score: i32,
        _query: &BlastCompoSequenceData,
        query_range: &BlastCompoSequenceRange,
        _ccat_query_length: i32,
        _subject: &BlastCompoSequenceData,
        subject_range: &BlastCompoSequenceRange,
        _full_subject_length: i32,
        _gapping_params: &BlastCompoGappingParams,
        matrix_adjust_rule: MatrixAdjustRule,
    ) -> i32 {
        *query_end += 1;
        *match_end += 1;
        *align = Some(BlastCompoAlignment::new(
            score,
            matrix_adjust_rule,
            query_range.context,
            query_range.begin + query_start,
            query_range.begin + *query_end,
            subject_range.begin + match_start,
            subject_range.begin + *match_end,
            subject_range.context,
            None,
        ));
        0
    }

    #[test]
    fn blast_redo_one_match_with_callbacks_redoes_no_composition() {
        let query_source = vec![0u8; 8];
        let query_info = vec![BlastCompoQueryInfo {
            origin: 0,
            seq: BlastCompoSequenceData {
                buffer: query_source,
                data_offset: 0,
                length: 8,
            },
            composition: BlastAminoAcidComposition::default(),
            eff_search_space: 100.0,
            words: Vec::new(),
        }];
        let incoming = Some(Box::new(BlastCompoAlignment::new(
            10,
            MatrixAdjustRule::DontAdjust,
            0,
            0,
            8,
            0,
            8,
            0,
            None,
        )));
        let mut params = redo_one_match_test_params(CompoAdjustMode::NoCompositionBasedStats);
        params.callbacks = Some(BlastRedoAlignCallbacks {
            get_range: Some(callback_test_get_range),
            redo_one_alignment: Some(callback_test_redo_one_alignment),
            ..Default::default()
        });
        let matching_seq = matching_sequence_initialize(8, 7, 0);
        let mut alignments = vec![None];

        let status = blast_redo_one_match_with_callbacks(
            &mut alignments,
            &params,
            &incoming,
            1,
            0.267,
            &matching_seq,
            &query_info,
        );

        assert_eq!(status, 0);
        let aln = alignments[0].as_ref().expect("callback alignment");
        assert_eq!(aln.score, 11);
        assert_eq!(aln.query_start, 0);
        assert_eq!(aln.match_end, 8);
    }

    #[test]
    fn blast_redo_one_match_with_callbacks_rejects_missing_callbacks() {
        let query_info = vec![BlastCompoQueryInfo {
            seq: BlastCompoSequenceData {
                buffer: vec![0u8; 4],
                data_offset: 0,
                length: 4,
            },
            ..Default::default()
        }];
        let incoming = Some(Box::new(BlastCompoAlignment::new(
            10,
            MatrixAdjustRule::DontAdjust,
            0,
            0,
            4,
            0,
            4,
            0,
            None,
        )));
        let params = redo_one_match_test_params(CompoAdjustMode::NoCompositionBasedStats);
        let matching_seq = matching_sequence_initialize(4, 7, 0);
        let mut alignments = vec![Some(Box::new(BlastCompoAlignment::new(
            9,
            MatrixAdjustRule::DontAdjust,
            0,
            0,
            4,
            0,
            4,
            0,
            None,
        )))];

        let status = blast_redo_one_match_with_callbacks(
            &mut alignments,
            &params,
            &incoming,
            1,
            0.267,
            &matching_seq,
            &query_info,
        );

        assert_eq!(status, -1);
        assert!(alignments[0].is_none());
    }

    #[test]
    fn blast_redo_one_match_with_callbacks_redoes_composition_adjusted() {
        let query_source = vec![1u8; 8];
        let (query_composition, query_num_true) =
            crate::composition::blast_read_aa_composition(&query_source, crate::matrix::AA_SIZE);
        let query_info = vec![BlastCompoQueryInfo {
            origin: 0,
            seq: BlastCompoSequenceData {
                buffer: query_source,
                data_offset: 0,
                length: 8,
            },
            composition: BlastAminoAcidComposition::from_prob(
                query_composition,
                query_num_true as i32,
            ),
            eff_search_space: 100.0,
            words: Vec::new(),
        }];
        let incoming = Some(Box::new(BlastCompoAlignment::new(
            10,
            MatrixAdjustRule::DontAdjust,
            0,
            0,
            8,
            0,
            8,
            0,
            None,
        )));
        let mut params = redo_one_match_test_params(CompoAdjustMode::CompositionBasedStats);
        assert_eq!(
            s_matrix_info_init(&mut params.matrix_info, "BLOSUM62", 0.3176, 1.0),
            0
        );
        params.callbacks = Some(BlastRedoAlignCallbacks {
            get_range: Some(callback_test_get_range_alanine),
            redo_one_alignment: Some(callback_test_redo_one_alignment),
            ..Default::default()
        });
        let matching_seq = matching_sequence_initialize(8, 7, 0);
        let mut alignments = vec![None];
        let mut matrix = square_matrix_from_vec(&params.matrix_info.matrix).expect("matrix");
        let mut pvalue = -1.0;
        let mut lambda_ratio = 0.0;

        let status = blast_redo_one_match_with_callbacks_and_adjustment(
            &mut alignments,
            &params,
            &incoming,
            1,
            0.267,
            &matching_seq,
            &query_info,
            &mut matrix,
            None,
            &mut pvalue,
            0,
            &mut lambda_ratio,
        );

        assert_eq!(status, 0);
        let aln = alignments[0].as_ref().expect("callback alignment");
        assert_eq!(aln.score, 11);
        assert_eq!(aln.matrix_adjust_rule, MatrixAdjustRule::ScaleOldMatrix);
        assert!((0.5..=1.0).contains(&lambda_ratio));
        assert_eq!(pvalue, -1.0);
    }

    #[test]
    fn blast_redo_one_match_with_callbacks_redoes_position_based_composition_adjusted() {
        let query_source = vec![1u8; 8];
        let (query_composition, query_num_true) =
            crate::composition::blast_read_aa_composition(&query_source, crate::matrix::AA_SIZE);
        let query_info = vec![BlastCompoQueryInfo {
            origin: 0,
            seq: BlastCompoSequenceData {
                buffer: query_source.clone(),
                data_offset: 0,
                length: query_source.len() as i32,
            },
            composition: BlastAminoAcidComposition::from_prob(
                query_composition,
                query_num_true as i32,
            ),
            eff_search_space: 100.0,
            words: Vec::new(),
        }];
        let incoming = Some(Box::new(BlastCompoAlignment::new(
            10,
            MatrixAdjustRule::DontAdjust,
            0,
            0,
            query_source.len() as i32,
            0,
            query_source.len() as i32,
            0,
            None,
        )));
        let mut params = redo_one_match_test_params(CompoAdjustMode::CompositionBasedStats);
        let start_numerator = vec![vec![0.0; crate::matrix::AA_SIZE]; query_source.len()];
        assert_eq!(
            matrix_info_init_psiblast_from_start_numerator(
                &mut params.matrix_info,
                &query_source,
                "BLOSUM62",
                &start_numerator,
                0.3176,
                1.0,
            ),
            0
        );
        params.position_based = true;
        params.gapping_params.gap_open = 11;
        params.gapping_params.gap_extend = 1;
        params.cutoff_s = 1;
        params.callbacks = Some(BlastRedoAlignCallbacks {
            get_range: Some(callback_test_get_range_alanine),
            ..Default::default()
        });
        let matching_seq = matching_sequence_initialize(query_source.len() as i32, 7, 0);
        let mut alignments = vec![None];
        let mut matrix = [[0i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE];
        let mut pvalue = -1.0;
        let mut lambda_ratio = 0.0;

        let status = blast_redo_one_match_with_callbacks_and_adjustment(
            &mut alignments,
            &params,
            &incoming,
            1,
            0.267,
            &matching_seq,
            &query_info,
            &mut matrix,
            None,
            &mut pvalue,
            0,
            &mut lambda_ratio,
        );

        assert_eq!(status, 0);
        let aln = alignments[0].as_ref().expect("callback PSSM alignment");
        assert_eq!(aln.query_start, 0);
        assert_eq!(aln.query_end, query_source.len() as i32);
        assert_eq!(aln.match_start, 0);
        assert_eq!(aln.match_end, query_source.len() as i32);
        assert_eq!(aln.matrix_adjust_rule, MatrixAdjustRule::ScaleOldMatrix);
        assert!((0.5..=1.0).contains(&lambda_ratio));
        assert_eq!(pvalue, -1.0);
    }

    #[test]
    fn blast_redo_one_match_smith_waterman_with_callbacks_redoes_no_composition() {
        let query_source = vec![0u8; 8];
        let query_info = vec![BlastCompoQueryInfo {
            origin: 0,
            seq: BlastCompoSequenceData {
                buffer: query_source,
                data_offset: 0,
                length: 8,
            },
            composition: BlastAminoAcidComposition::default(),
            eff_search_space: 100.0,
            words: Vec::new(),
        }];
        let incoming = Some(Box::new(BlastCompoAlignment::new(
            10,
            MatrixAdjustRule::DontAdjust,
            0,
            0,
            8,
            0,
            8,
            0,
            None,
        )));
        let mut params = redo_one_match_test_params(CompoAdjustMode::NoCompositionBasedStats);
        params.callbacks = Some(BlastRedoAlignCallbacks {
            get_range: Some(callback_test_get_range),
            new_xdrop_align: Some(callback_test_new_xdrop_align),
            ..Default::default()
        });
        let matching_seq = matching_sequence_initialize(8, 7, 0);
        let mut alignments = vec![None];
        let heaps = vec![BlastCompoHeap::new(10, 10.0)];
        let mut sw_matrix = [[-3i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE];
        sw_matrix[0][0] = 2;

        let status = blast_redo_one_match_smith_waterman_with_callbacks(
            &mut alignments,
            &params,
            &incoming,
            1,
            0.267,
            0.041_f64.ln(),
            &matching_seq,
            &query_info,
            &heaps,
            &sw_matrix,
        );

        assert_eq!(status, 0);
        let aln = alignments[0].as_ref().expect("callback SW alignment");
        assert_eq!(aln.score, 16);
        assert_eq!(aln.query_start, 0);
        assert_eq!(aln.query_end, 8);
        assert_eq!(aln.match_start, 0);
        assert_eq!(aln.match_end, 8);
    }

    #[test]
    fn blast_redo_one_match_smith_waterman_with_callbacks_redoes_composition_adjusted() {
        let query_source = vec![1u8; 8];
        let (query_composition, query_num_true) =
            crate::composition::blast_read_aa_composition(&query_source, crate::matrix::AA_SIZE);
        let query_info = vec![BlastCompoQueryInfo {
            origin: 0,
            seq: BlastCompoSequenceData {
                buffer: query_source,
                data_offset: 0,
                length: 8,
            },
            composition: BlastAminoAcidComposition::from_prob(
                query_composition,
                query_num_true as i32,
            ),
            eff_search_space: 100.0,
            words: Vec::new(),
        }];
        let incoming = Some(Box::new(BlastCompoAlignment::new(
            10,
            MatrixAdjustRule::DontAdjust,
            0,
            0,
            8,
            0,
            8,
            0,
            None,
        )));
        let mut params = redo_one_match_test_params(CompoAdjustMode::CompositionBasedStats);
        assert_eq!(
            s_matrix_info_init(&mut params.matrix_info, "BLOSUM62", 0.3176, 1.0),
            0
        );
        params.callbacks = Some(BlastRedoAlignCallbacks {
            get_range: Some(callback_test_get_range_alanine),
            new_xdrop_align: Some(callback_test_new_xdrop_align),
            ..Default::default()
        });
        let matching_seq = matching_sequence_initialize(8, 7, 0);
        let mut alignments = vec![None];
        let heaps = vec![BlastCompoHeap::new(10, 10.0)];
        let mut matrix = square_matrix_from_vec(&params.matrix_info.matrix).expect("matrix");
        let mut pvalue = -1.0;
        let mut lambda_ratio = 0.0;

        let status = blast_redo_one_match_smith_waterman_with_callbacks_and_adjustment(
            &mut alignments,
            &params,
            &incoming,
            1,
            0.267,
            0.041_f64.ln(),
            &matching_seq,
            &query_info,
            &heaps,
            &mut matrix,
            None,
            &mut pvalue,
            0,
            &mut lambda_ratio,
        );

        assert_eq!(status, 0);
        let aln = alignments[0].as_ref().expect("callback SW alignment");
        assert!(aln.score > 0);
        assert_eq!(aln.matrix_adjust_rule, MatrixAdjustRule::ScaleOldMatrix);
        assert!((0.5..=1.0).contains(&lambda_ratio));
        assert_eq!(pvalue, -1.0);
    }

    #[test]
    fn blast_redo_one_match_smith_waterman_with_callbacks_redoes_position_based_composition_adjusted(
    ) {
        let query_source = vec![1u8; 8];
        let (query_composition, query_num_true) =
            crate::composition::blast_read_aa_composition(&query_source, crate::matrix::AA_SIZE);
        let query_info = vec![BlastCompoQueryInfo {
            origin: 0,
            seq: BlastCompoSequenceData {
                buffer: query_source.clone(),
                data_offset: 0,
                length: query_source.len() as i32,
            },
            composition: BlastAminoAcidComposition::from_prob(
                query_composition,
                query_num_true as i32,
            ),
            eff_search_space: 100.0,
            words: Vec::new(),
        }];
        let incoming = Some(Box::new(BlastCompoAlignment::new(
            10,
            MatrixAdjustRule::DontAdjust,
            0,
            0,
            query_source.len() as i32,
            0,
            query_source.len() as i32,
            0,
            None,
        )));
        let mut params = redo_one_match_test_params(CompoAdjustMode::CompositionBasedStats);
        let start_numerator = vec![vec![0.0; crate::matrix::AA_SIZE]; query_source.len()];
        assert_eq!(
            matrix_info_init_psiblast_from_start_numerator(
                &mut params.matrix_info,
                &query_source,
                "BLOSUM62",
                &start_numerator,
                0.3176,
                1.0,
            ),
            0
        );
        params.position_based = true;
        params.gapping_params.gap_open = 11;
        params.gapping_params.gap_extend = 1;
        params.cutoff_s = 1;
        params.callbacks = Some(BlastRedoAlignCallbacks {
            get_range: Some(callback_test_get_range_alanine),
            ..Default::default()
        });
        let matching_seq = matching_sequence_initialize(query_source.len() as i32, 7, 0);
        let mut alignments = vec![None];
        let heaps = vec![BlastCompoHeap::new(10, 10.0)];
        let mut matrix = [[0i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE];
        let mut pvalue = -1.0;
        let mut lambda_ratio = 0.0;

        let status = blast_redo_one_match_smith_waterman_with_callbacks_and_adjustment(
            &mut alignments,
            &params,
            &incoming,
            1,
            0.267,
            0.041_f64.ln(),
            &matching_seq,
            &query_info,
            &heaps,
            &mut matrix,
            None,
            &mut pvalue,
            0,
            &mut lambda_ratio,
        );

        assert_eq!(status, 0);
        let aln = alignments[0].as_ref().expect("callback PSSM SW alignment");
        assert_eq!(aln.query_start, 0);
        assert_eq!(aln.query_end, query_source.len() as i32);
        assert_eq!(aln.match_start, 0);
        assert_eq!(aln.match_end, query_source.len() as i32);
        assert_eq!(aln.matrix_adjust_rule, MatrixAdjustRule::ScaleOldMatrix);
        assert!(aln.context.is_some());
        assert!((0.5..=1.0).contains(&lambda_ratio));
        assert_eq!(pvalue, -1.0);
    }

    #[test]
    fn blast_adjust_scores_scale_old_matrix_updates_matrix_and_ratio() {
        let mut matrix_info = BlastMatrixInfo::default();
        assert_eq!(
            s_matrix_info_init(&mut matrix_info, "BLOSUM62", 0.3176, 1.0),
            0
        );
        let mut matrix = square_matrix_from_vec(&matrix_info.matrix).expect("start matrix");
        let original = matrix;
        let mut query_prob = vec![0.0f64; crate::matrix::AA_SIZE];
        let mut subject_prob = vec![0.0f64; crate::matrix::AA_SIZE];
        for &idx in &crate::composition::TRUE_CHAR_POSITIONS {
            query_prob[idx] = 1.0 / crate::composition::COMPO_NUM_TRUE_AA as f64;
            subject_prob[idx] = 1.0 / crate::composition::COMPO_NUM_TRUE_AA as f64;
        }
        subject_prob[crate::composition::TRUE_CHAR_POSITIONS[0]] += 0.19;
        for &idx in &crate::composition::TRUE_CHAR_POSITIONS[1..] {
            subject_prob[idx] -= 0.01;
        }
        let mut pvalue = 0.0;
        let mut ratio = 0.0;

        let rule = blast_adjust_scores_scale_old_matrix(
            &mut matrix,
            &query_prob,
            crate::composition::COMPO_NUM_TRUE_AA,
            &subject_prob,
            crate::composition::COMPO_NUM_TRUE_AA,
            &matrix_info,
            CompoAdjustMode::CompositionBasedStats,
            0,
            &mut pvalue,
            &mut ratio,
        )
        .expect("scale old matrix");

        assert_eq!(rule, MatrixAdjustRule::ScaleOldMatrix);
        assert!((0.5..=1.0).contains(&ratio));
        assert_ne!(matrix, original);
        assert_eq!(pvalue, 0.0);
    }

    #[test]
    fn blast_adjust_scores_pvalue_test_computes_pvalue_and_scales_matrix() {
        let mut matrix_info = BlastMatrixInfo::default();
        assert_eq!(
            s_matrix_info_init(&mut matrix_info, "BLOSUM62", 0.3176, 1.0),
            0
        );
        let mut matrix = square_matrix_from_vec(&matrix_info.matrix).expect("start matrix");
        let original = matrix;
        let mut query_prob = vec![0.0f64; crate::matrix::AA_SIZE];
        let mut subject_prob = vec![0.0f64; crate::matrix::AA_SIZE];
        for &idx in &crate::composition::TRUE_CHAR_POSITIONS {
            query_prob[idx] = 1.0 / crate::composition::COMPO_NUM_TRUE_AA as f64;
            subject_prob[idx] = 1.0 / crate::composition::COMPO_NUM_TRUE_AA as f64;
        }
        query_prob[crate::composition::TRUE_CHAR_POSITIONS[4]] += 0.095;
        query_prob[crate::composition::TRUE_CHAR_POSITIONS[0]] -= 0.005;
        for &idx in &crate::composition::TRUE_CHAR_POSITIONS[1..4] {
            query_prob[idx] -= 0.03;
        }
        let mut pvalue = 0.0;
        let mut ratio = 0.0;

        let rule = blast_adjust_scores_with_workspace(
            &mut matrix,
            &query_prob,
            crate::composition::COMPO_NUM_TRUE_AA,
            &subject_prob,
            crate::composition::COMPO_NUM_TRUE_AA,
            &matrix_info,
            CompoAdjustMode::CompositionBasedStats,
            0,
            None,
            1,
            &mut pvalue,
            &mut ratio,
        )
        .expect("p-value test scale old matrix");

        assert_eq!(rule, MatrixAdjustRule::ScaleOldMatrix);
        assert!((0.0..=1.0).contains(&pvalue));
        assert!(pvalue > 0.0);
        assert!((0.5..=1.0).contains(&ratio));
        assert_ne!(matrix, original);
    }

    #[test]
    fn blast_adjust_position_based_scores_scales_pssm_rows() {
        let query = [1u8, 20u8];
        let start_numerator = vec![vec![0.0; crate::matrix::AA_SIZE]; query.len()];
        let mut matrix_info = BlastMatrixInfo::default();
        assert_eq!(
            matrix_info_init_psiblast_from_start_numerator(
                &mut matrix_info,
                &query,
                "BLOSUM62",
                &start_numerator,
                0.3176,
                1.0,
            ),
            0
        );
        let original = matrix_info.matrix.clone();
        let mut matrix = original.clone();
        let mut query_prob = vec![0.0f64; crate::matrix::AA_SIZE];
        let mut subject_prob = vec![0.0f64; crate::matrix::AA_SIZE];
        for &idx in &crate::composition::TRUE_CHAR_POSITIONS {
            query_prob[idx] = 1.0 / crate::composition::COMPO_NUM_TRUE_AA as f64;
            subject_prob[idx] = 1.0 / crate::composition::COMPO_NUM_TRUE_AA as f64;
        }
        subject_prob[crate::composition::TRUE_CHAR_POSITIONS[0]] += 0.19;
        for &idx in &crate::composition::TRUE_CHAR_POSITIONS[1..] {
            subject_prob[idx] -= 0.01;
        }
        let mut pvalue = 0.0;
        let mut ratio = 0.0;

        let rule = blast_adjust_position_based_scores(
            &mut matrix,
            &query_prob,
            crate::composition::COMPO_NUM_TRUE_AA,
            &subject_prob,
            crate::composition::COMPO_NUM_TRUE_AA,
            &matrix_info,
            1,
            &mut pvalue,
            &mut ratio,
        )
        .expect("position-based adjustment");

        assert_eq!(rule, MatrixAdjustRule::ScaleOldMatrix);
        assert_eq!(matrix.len(), query.len());
        assert_eq!(
            matrix[0][crate::encoding::NCBISTDAA_STOP as usize],
            original[0][crate::encoding::NCBISTDAA_STOP as usize]
        );
        assert_eq!(
            matrix[0][crate::encoding::NCBISTDAA_U as usize],
            matrix[0][crate::encoding::NCBISTDAA_C as usize]
        );
        assert!(matrix[0][crate::encoding::NCBISTDAA_X as usize] <= -1);
        assert!(pvalue > 0.0 && pvalue <= 1.0);
        assert!(ratio >= 0.5);
    }

    #[test]
    fn blast_adjust_scores_scale_old_matrix_returns_one_for_empty_composition() {
        let mut matrix_info = BlastMatrixInfo::default();
        assert_eq!(
            s_matrix_info_init(&mut matrix_info, "BLOSUM62", 0.3176, 1.0),
            0
        );
        let mut matrix = square_matrix_from_vec(&matrix_info.matrix).expect("start matrix");
        let probs = vec![0.0f64; crate::matrix::AA_SIZE];
        let mut pvalue = 0.0;
        let mut ratio = 0.0;

        let status = blast_adjust_scores_scale_old_matrix(
            &mut matrix,
            &probs,
            0,
            &probs,
            0,
            &matrix_info,
            CompoAdjustMode::CompositionBasedStats,
            0,
            &mut pvalue,
            &mut ratio,
        )
        .expect_err("empty composition skips adjustment");

        assert_eq!(status, 1);
    }

    #[test]
    fn blast_adjust_scores_with_workspace_runs_re_optimization() {
        let mut matrix_info = BlastMatrixInfo::default();
        assert_eq!(
            s_matrix_info_init(&mut matrix_info, "BLOSUM62", 0.3176, 1.0),
            0
        );
        let mut matrix = square_matrix_from_vec(&matrix_info.matrix).expect("start matrix");
        let original = matrix;
        let mut query_prob = vec![0.0f64; crate::matrix::AA_SIZE];
        let mut subject_prob = vec![0.0f64; crate::matrix::AA_SIZE];
        for &idx in &crate::composition::TRUE_CHAR_POSITIONS {
            query_prob[idx] = 1.0 / crate::composition::COMPO_NUM_TRUE_AA as f64;
            subject_prob[idx] = 1.0 / crate::composition::COMPO_NUM_TRUE_AA as f64;
        }
        let workspace = BlastCompositionWorkspace::blosum62();
        let mut pvalue = 0.0;
        let mut ratio = 0.0;

        let rule = blast_adjust_scores_with_workspace(
            &mut matrix,
            &query_prob,
            crate::composition::COMPO_NUM_TRUE_AA,
            &subject_prob,
            crate::composition::COMPO_NUM_TRUE_AA,
            &matrix_info,
            CompoAdjustMode::CompoForceFullMatrixAdjust,
            K_RE_MATRIX_ADJUSTMENT_PSEUDOCOUNTS,
            Some(&workspace),
            0,
            &mut pvalue,
            &mut ratio,
        )
        .expect("RE matrix adjustment");

        assert_eq!(rule, MatrixAdjustRule::UserSpecifiedRelEntropy);
        assert_eq!(ratio, 1.0);
        assert_ne!(matrix, original);
        assert_eq!(pvalue, 0.0);
    }

    #[test]
    fn blast_adjust_scores_with_workspace_requires_workspace_for_re_mode() {
        let mut matrix_info = BlastMatrixInfo::default();
        assert_eq!(
            s_matrix_info_init(&mut matrix_info, "BLOSUM62", 0.3176, 1.0),
            0
        );
        let mut matrix = square_matrix_from_vec(&matrix_info.matrix).expect("start matrix");
        let mut probs = vec![0.0f64; crate::matrix::AA_SIZE];
        for &idx in &crate::composition::TRUE_CHAR_POSITIONS {
            probs[idx] = 1.0 / crate::composition::COMPO_NUM_TRUE_AA as f64;
        }
        let mut pvalue = 0.0;
        let mut ratio = 0.0;

        let status = blast_adjust_scores_with_workspace(
            &mut matrix,
            &probs,
            crate::composition::COMPO_NUM_TRUE_AA,
            &probs,
            crate::composition::COMPO_NUM_TRUE_AA,
            &matrix_info,
            CompoAdjustMode::CompoForceFullMatrixAdjust,
            K_RE_MATRIX_ADJUSTMENT_PSEUDOCOUNTS,
            None,
            0,
            &mut pvalue,
            &mut ratio,
        )
        .expect_err("missing workspace");

        assert_eq!(status, -1);
    }

    #[test]
    fn blast_redo_one_match_in_memory_no_composition_redoes_alignment() {
        let query_source = vec![0u8, 1, 2, 3, 0, 1, 2, 3, 0, 1];
        let subject_source = query_source.clone();
        let query_info = vec![BlastCompoQueryInfo {
            origin: 0,
            seq: BlastCompoSequenceData {
                buffer: query_source.clone(),
                data_offset: 0,
                length: query_source.len() as i32,
            },
            composition: BlastAminoAcidComposition::default(),
            eff_search_space: 100.0,
            words: Vec::new(),
        }];
        let incoming = Some(Box::new(BlastCompoAlignment::new(
            10,
            MatrixAdjustRule::DontAdjust,
            0,
            0,
            10,
            0,
            10,
            0,
            None,
        )));
        let params = redo_one_match_test_params(CompoAdjustMode::NoCompositionBasedStats);
        let matching_seq = matching_sequence_initialize(subject_source.len() as i32, 7, 0);
        let mut alignments = vec![None];

        let status = blast_redo_one_match_in_memory(
            &mut alignments,
            &params,
            &incoming,
            1,
            0.267,
            &matching_seq,
            &query_info,
            crate::program::BLASTN,
            &subject_source,
            1,
            -3,
            &crate::util::STANDARD_GENETIC_CODE,
        );

        assert_eq!(status, 0);
        let aln = alignments[0].as_ref().expect("redone alignment");
        assert!(aln.score > 0);
        assert_eq!(aln.query_start, 0);
        assert_eq!(aln.query_end, 10);
        assert_eq!(aln.match_start, 0);
        assert_eq!(aln.match_end, 10);
        assert_eq!(aln.frame, 0);
    }

    #[test]
    fn blast_redo_one_match_in_memory_no_composition_redoes_protein_alignment() {
        let query_source = vec![12u8, 1, 12, 1, 12, 1, 12, 1];
        let subject_source = query_source.clone();
        let query_info = vec![BlastCompoQueryInfo {
            origin: 0,
            seq: BlastCompoSequenceData {
                buffer: query_source.clone(),
                data_offset: 0,
                length: query_source.len() as i32,
            },
            composition: BlastAminoAcidComposition::default(),
            eff_search_space: 100.0,
            words: Vec::new(),
        }];
        let incoming = Some(Box::new(BlastCompoAlignment::new(
            10,
            MatrixAdjustRule::DontAdjust,
            0,
            0,
            query_source.len() as i32,
            0,
            subject_source.len() as i32,
            0,
            None,
        )));
        let mut params = redo_one_match_test_params(CompoAdjustMode::NoCompositionBasedStats);
        assert_eq!(
            s_matrix_info_init(&mut params.matrix_info, "BLOSUM62", 0.3176, 1.0),
            0
        );
        params.gapping_params.gap_open = 11;
        params.gapping_params.gap_extend = 1;
        params.gapping_params.decline_align = i32::MIN;
        params.ccat_query_length = query_source.len() as i32;
        let matching_seq = matching_sequence_initialize(subject_source.len() as i32, 7, 0);
        let mut alignments = vec![None];

        let status = blast_redo_one_match_in_memory(
            &mut alignments,
            &params,
            &incoming,
            1,
            0.267,
            &matching_seq,
            &query_info,
            crate::program::BLASTP,
            &subject_source,
            0,
            0,
            &crate::util::STANDARD_GENETIC_CODE,
        );

        assert_eq!(status, 0);
        let aln = alignments[0].as_ref().expect("redone protein alignment");
        assert!(aln.score > 0);
        assert_eq!(aln.query_start, 0);
        assert_eq!(aln.query_end, query_source.len() as i32);
        assert_eq!(aln.match_start, 0);
        assert_eq!(aln.match_end, subject_source.len() as i32);
        assert_eq!(aln.frame, 0);
    }

    #[test]
    fn blast_redo_one_match_in_memory_no_composition_redoes_tblastn_alignment() {
        let query_source = vec![12u8, 1, 12, 1];
        let subject_source = vec![1u8, 8, 4, 4, 2, 8, 1, 8, 4, 4, 2, 8];
        let query_info = vec![BlastCompoQueryInfo {
            origin: 0,
            seq: BlastCompoSequenceData {
                buffer: query_source.clone(),
                data_offset: 0,
                length: query_source.len() as i32,
            },
            composition: BlastAminoAcidComposition::default(),
            eff_search_space: 100.0,
            words: Vec::new(),
        }];
        let incoming = Some(Box::new(BlastCompoAlignment::new(
            10,
            MatrixAdjustRule::DontAdjust,
            0,
            0,
            query_source.len() as i32,
            0,
            query_source.len() as i32,
            1,
            None,
        )));
        let mut params = redo_one_match_test_params(CompoAdjustMode::NoCompositionBasedStats);
        assert_eq!(
            s_matrix_info_init(&mut params.matrix_info, "BLOSUM62", 0.3176, 1.0),
            0
        );
        params.gapping_params.gap_open = 11;
        params.gapping_params.gap_extend = 1;
        params.gapping_params.decline_align = i32::MIN;
        params.subject_is_translated = true;
        params.ccat_query_length = query_source.len() as i32;
        let matching_seq = matching_sequence_initialize(subject_source.len() as i32, 7, 0);
        let mut alignments = vec![None];

        let status = blast_redo_one_match_in_memory(
            &mut alignments,
            &params,
            &incoming,
            1,
            0.267,
            &matching_seq,
            &query_info,
            crate::program::TBLASTN,
            &subject_source,
            0,
            0,
            &crate::util::STANDARD_GENETIC_CODE,
        );

        assert_eq!(status, 0);
        let aln = alignments[0]
            .as_ref()
            .expect("redone translated-subject alignment");
        assert!(aln.score > 0);
        assert_eq!(aln.query_start, 0);
        assert_eq!(aln.query_end, query_source.len() as i32);
        assert_eq!(aln.match_start, 0);
        assert_eq!(aln.match_end, query_source.len() as i32);
        assert_eq!(aln.frame, 1);
    }

    #[test]
    fn blast_redo_one_match_in_memory_no_composition_redoes_blastx_alignment() {
        let query_ncbi4na = vec![1u8, 8, 4, 4, 2, 8, 1, 8, 4, 4, 2, 8];
        let (query_source, offsets) = crate::util::blast_get_all_translations(
            &query_ncbi4na,
            query_ncbi4na.len(),
            &crate::util::STANDARD_GENETIC_CODE,
        );
        let qi = crate::queryinfo::QueryInfo::new_translated_query_from_offsets(&offsets);
        let query_info = get_query_info(&query_source, &qi, true);
        let subject_source = vec![12u8, 1, 12, 1];
        let incoming = Some(Box::new(BlastCompoAlignment::new(
            10,
            MatrixAdjustRule::DontAdjust,
            0,
            0,
            subject_source.len() as i32,
            0,
            subject_source.len() as i32,
            0,
            None,
        )));
        let mut params = redo_one_match_test_params(CompoAdjustMode::NoCompositionBasedStats);
        assert_eq!(
            s_matrix_info_init(&mut params.matrix_info, "BLOSUM62", 0.3176, 1.0),
            0
        );
        params.gapping_params.gap_open = 11;
        params.gapping_params.gap_extend = 1;
        params.gapping_params.decline_align = i32::MIN;
        params.query_is_translated = true;
        params.ccat_query_length = query_ncbi4na.len() as i32;
        let matching_seq = matching_sequence_initialize(subject_source.len() as i32, 7, 0);
        let mut alignments = vec![None; query_info.len()];

        let status = blast_redo_one_match_in_memory(
            &mut alignments,
            &params,
            &incoming,
            1,
            0.267,
            &matching_seq,
            &query_info,
            crate::program::BLASTX,
            &subject_source,
            0,
            0,
            &crate::util::STANDARD_GENETIC_CODE,
        );

        assert_eq!(status, 0);
        let aln = alignments[0].as_ref().expect("redone blastx alignment");
        assert!(aln.score > 0);
        assert_eq!(aln.query_start, 0);
        assert_eq!(aln.query_end, subject_source.len() as i32);
        assert_eq!(aln.match_start, 0);
        assert_eq!(aln.match_end, subject_source.len() as i32);
        assert_eq!(aln.frame, 0);
    }

    #[test]
    fn blast_redo_one_match_in_memory_no_composition_redoes_blastx_second_context() {
        let query_ncbi4na = vec![1u8, 8, 4, 4, 2, 8, 1, 8, 4, 4, 2, 8];
        let (query_source, offsets) = crate::util::blast_get_all_translations(
            &query_ncbi4na,
            query_ncbi4na.len(),
            &crate::util::STANDARD_GENETIC_CODE,
        );
        let qi = crate::queryinfo::QueryInfo::new_translated_query_from_offsets(&offsets);
        let query_info = get_query_info(&query_source, &qi, true);
        let subject_source = query_info[1].seq.data().to_vec();
        let incoming = Some(Box::new(BlastCompoAlignment::new(
            10,
            MatrixAdjustRule::DontAdjust,
            1,
            0,
            subject_source.len() as i32,
            0,
            subject_source.len() as i32,
            0,
            None,
        )));
        let mut params = redo_one_match_test_params(CompoAdjustMode::NoCompositionBasedStats);
        assert_eq!(
            s_matrix_info_init(&mut params.matrix_info, "BLOSUM62", 0.3176, 1.0),
            0
        );
        params.gapping_params.gap_open = 11;
        params.gapping_params.gap_extend = 1;
        params.gapping_params.decline_align = i32::MIN;
        params.query_is_translated = true;
        params.ccat_query_length = query_ncbi4na.len() as i32;
        let matching_seq = matching_sequence_initialize(subject_source.len() as i32, 7, 0);
        let mut alignments = vec![None; query_info.len()];

        let status = blast_redo_one_match_in_memory(
            &mut alignments,
            &params,
            &incoming,
            1,
            0.267,
            &matching_seq,
            &query_info,
            crate::program::BLASTX,
            &subject_source,
            0,
            0,
            &crate::util::STANDARD_GENETIC_CODE,
        );

        assert_eq!(status, 0);
        assert!(alignments[0].is_none());
        let aln = alignments[1]
            .as_ref()
            .expect("redone second-context blastx alignment");
        assert!(aln.score > 0);
        assert_eq!(aln.query_start, 0);
        assert_eq!(aln.query_end, subject_source.len() as i32);
        assert_eq!(aln.match_start, 0);
        assert_eq!(aln.match_end, subject_source.len() as i32);
        assert_eq!(aln.frame, 0);
    }

    #[test]
    fn blast_redo_one_match_in_memory_no_composition_redoes_tblastx_alignment() {
        let query_ncbi4na = vec![1u8, 8, 4, 4, 2, 8, 1, 8, 4, 4, 2, 8];
        let (query_source, offsets) = crate::util::blast_get_all_translations(
            &query_ncbi4na,
            query_ncbi4na.len(),
            &crate::util::STANDARD_GENETIC_CODE,
        );
        let qi = crate::queryinfo::QueryInfo::new_translated_query_from_offsets(&offsets);
        let query_info = get_query_info(&query_source, &qi, true);
        let subject_source = query_ncbi4na.clone();
        let incoming = Some(Box::new(BlastCompoAlignment::new(
            10,
            MatrixAdjustRule::DontAdjust,
            0,
            0,
            4,
            0,
            4,
            1,
            None,
        )));
        let mut params = redo_one_match_test_params(CompoAdjustMode::NoCompositionBasedStats);
        assert_eq!(
            s_matrix_info_init(&mut params.matrix_info, "BLOSUM62", 0.3176, 1.0),
            0
        );
        params.gapping_params.gap_open = 11;
        params.gapping_params.gap_extend = 1;
        params.gapping_params.decline_align = i32::MIN;
        params.query_is_translated = true;
        params.subject_is_translated = true;
        params.ccat_query_length = query_ncbi4na.len() as i32;
        let matching_seq = matching_sequence_initialize(subject_source.len() as i32, 7, 0);
        let mut alignments = vec![None; query_info.len()];

        let status = blast_redo_one_match_in_memory(
            &mut alignments,
            &params,
            &incoming,
            1,
            0.267,
            &matching_seq,
            &query_info,
            crate::program::TBLASTX,
            &subject_source,
            0,
            0,
            &crate::util::STANDARD_GENETIC_CODE,
        );

        assert_eq!(status, 0);
        let aln = alignments[0].as_ref().expect("redone tblastx alignment");
        assert!(aln.score > 0);
        assert_eq!(aln.query_start, 0);
        assert_eq!(aln.query_end, 4);
        assert_eq!(aln.match_start, 0);
        assert_eq!(aln.match_end, 4);
        assert_eq!(aln.frame, 1);
    }

    #[test]
    fn blast_redo_one_match_in_memory_no_composition_redoes_tblastx_second_context() {
        let query_ncbi4na = vec![1u8, 8, 4, 4, 2, 8, 1, 8, 4, 4, 2, 8];
        let (query_source, offsets) = crate::util::blast_get_all_translations(
            &query_ncbi4na,
            query_ncbi4na.len(),
            &crate::util::STANDARD_GENETIC_CODE,
        );
        let qi = crate::queryinfo::QueryInfo::new_translated_query_from_offsets(&offsets);
        let query_info = get_query_info(&query_source, &qi, true);
        let subject_source = query_ncbi4na.clone();
        let query_len = query_info[1].seq.length;
        let incoming = Some(Box::new(BlastCompoAlignment::new(
            10,
            MatrixAdjustRule::DontAdjust,
            1,
            0,
            query_len,
            0,
            query_len,
            2,
            None,
        )));
        let mut params = redo_one_match_test_params(CompoAdjustMode::NoCompositionBasedStats);
        assert_eq!(
            s_matrix_info_init(&mut params.matrix_info, "BLOSUM62", 0.3176, 1.0),
            0
        );
        params.gapping_params.gap_open = 11;
        params.gapping_params.gap_extend = 1;
        params.gapping_params.decline_align = i32::MIN;
        params.query_is_translated = true;
        params.subject_is_translated = true;
        params.ccat_query_length = query_ncbi4na.len() as i32;
        let matching_seq = matching_sequence_initialize(subject_source.len() as i32, 7, 0);
        let mut alignments = vec![None; query_info.len()];

        let status = blast_redo_one_match_in_memory(
            &mut alignments,
            &params,
            &incoming,
            1,
            0.267,
            &matching_seq,
            &query_info,
            crate::program::TBLASTX,
            &subject_source,
            0,
            0,
            &crate::util::STANDARD_GENETIC_CODE,
        );

        assert_eq!(status, 0);
        assert!(alignments[0].is_none());
        let aln = alignments[1]
            .as_ref()
            .expect("redone second-context tblastx alignment");
        assert!(aln.score > 0);
        assert_eq!(aln.query_start, 0);
        assert_eq!(aln.query_end, query_len);
        assert_eq!(aln.match_start, 0);
        assert_eq!(aln.match_end, query_len);
        assert_eq!(aln.frame, 2);
    }

    #[test]
    fn blast_redo_one_match_in_memory_rejects_composition_adjusted_path() {
        let query_source = vec![0u8, 1, 2, 3];
        let query_info = vec![BlastCompoQueryInfo {
            origin: 0,
            seq: BlastCompoSequenceData {
                buffer: query_source.clone(),
                data_offset: 0,
                length: query_source.len() as i32,
            },
            composition: BlastAminoAcidComposition::default(),
            eff_search_space: 100.0,
            words: Vec::new(),
        }];
        let incoming = Some(Box::new(BlastCompoAlignment::new(
            10,
            MatrixAdjustRule::DontAdjust,
            0,
            0,
            4,
            0,
            4,
            0,
            None,
        )));
        let params = redo_one_match_test_params(CompoAdjustMode::CompositionBasedStats);
        let matching_seq = matching_sequence_initialize(query_source.len() as i32, 7, 0);
        let mut alignments = vec![Some(Box::new(BlastCompoAlignment::new(
            9,
            MatrixAdjustRule::DontAdjust,
            0,
            0,
            4,
            0,
            4,
            0,
            None,
        )))];

        let status = blast_redo_one_match_in_memory(
            &mut alignments,
            &params,
            &incoming,
            1,
            0.267,
            &matching_seq,
            &query_info,
            crate::program::BLASTN,
            &query_source,
            1,
            -3,
            &crate::util::STANDARD_GENETIC_CODE,
        );

        assert_eq!(status, -1);
        assert!(alignments[0].is_none());
    }

    #[test]
    fn blast_redo_one_match_smith_waterman_in_memory_nucl_redoes_alignment() {
        let query_source = vec![0u8, 1, 2, 3, 0, 1, 2, 3, 0, 1];
        let subject_source = query_source.clone();
        let query_info = vec![BlastCompoQueryInfo {
            origin: 0,
            seq: BlastCompoSequenceData {
                buffer: query_source.clone(),
                data_offset: 0,
                length: query_source.len() as i32,
            },
            composition: BlastAminoAcidComposition::default(),
            eff_search_space: 100.0,
            words: Vec::new(),
        }];
        let incoming = Some(Box::new(BlastCompoAlignment::new(
            10,
            MatrixAdjustRule::DontAdjust,
            0,
            0,
            10,
            0,
            10,
            0,
            None,
        )));
        let params = redo_one_match_test_params(CompoAdjustMode::NoCompositionBasedStats);
        let matching_seq = matching_sequence_initialize(subject_source.len() as i32, 7, 0);
        let mut alignments = vec![None];
        let heaps = vec![BlastCompoHeap::new(10, 10.0)];
        let matrix = crate::matrix::nucleotide_matrix(1, -3);

        let status = blast_redo_one_match_smith_waterman_in_memory_nucl(
            &mut alignments,
            &params,
            &incoming,
            1,
            0.267,
            0.041_f64.ln(),
            &matching_seq,
            &query_info,
            &subject_source,
            &heaps,
            &matrix,
        );

        assert_eq!(status, 0);
        let aln = alignments[0].as_ref().expect("redone SW alignment");
        assert!(aln.score > 0);
        assert_eq!(aln.query_start, 0);
        assert_eq!(aln.query_end, 10);
        assert_eq!(aln.match_start, 0);
        assert_eq!(aln.match_end, 10);
        assert!(aln.context.is_some());
    }

    #[test]
    fn blast_redo_one_match_smith_waterman_in_memory_protein_redoes_alignment() {
        let query_source = vec![1u8; 8];
        let subject_source = query_source.clone();
        let query_info = vec![BlastCompoQueryInfo {
            origin: 0,
            seq: BlastCompoSequenceData {
                buffer: query_source.clone(),
                data_offset: 0,
                length: query_source.len() as i32,
            },
            composition: BlastAminoAcidComposition::default(),
            eff_search_space: 100.0,
            words: Vec::new(),
        }];
        let incoming = Some(Box::new(BlastCompoAlignment::new(
            10,
            MatrixAdjustRule::DontAdjust,
            0,
            0,
            query_source.len() as i32,
            0,
            subject_source.len() as i32,
            0,
            None,
        )));
        let mut params = redo_one_match_test_params(CompoAdjustMode::NoCompositionBasedStats);
        params.gapping_params.gap_open = 11;
        params.gapping_params.gap_extend = 1;
        assert_eq!(
            s_matrix_info_init(&mut params.matrix_info, "BLOSUM62", 0.3176, 1.0),
            0
        );
        let matching_seq = matching_sequence_initialize(subject_source.len() as i32, 7, 0);
        let mut alignments = vec![None];
        let heaps = vec![BlastCompoHeap::new(10, 10.0)];

        let status = blast_redo_one_match_smith_waterman_in_memory_protein(
            &mut alignments,
            &params,
            &incoming,
            1,
            0.267,
            0.041_f64.ln(),
            &matching_seq,
            &query_info,
            crate::program::BLASTP,
            &subject_source,
            &heaps,
            &crate::util::STANDARD_GENETIC_CODE,
        );

        assert_eq!(status, 0);
        let aln = alignments[0].as_ref().expect("redone protein SW alignment");
        assert!(aln.score > 0);
        assert_eq!(aln.query_start, 0);
        assert_eq!(aln.query_end, query_source.len() as i32);
        assert_eq!(aln.match_start, 0);
        assert_eq!(aln.match_end, subject_source.len() as i32);
        assert!(aln.context.is_some());
    }

    #[test]
    fn blast_redo_one_match_smith_waterman_keeps_going_while_significant() {
        let query_source = vec![1u8; 4];
        let subject_source = vec![1u8; 16];
        let query_info = vec![BlastCompoQueryInfo {
            origin: 0,
            seq: BlastCompoSequenceData {
                buffer: query_source.clone(),
                data_offset: 0,
                length: query_source.len() as i32,
            },
            composition: BlastAminoAcidComposition::default(),
            eff_search_space: 100.0,
            words: Vec::new(),
        }];
        let mut second = BlastCompoAlignment::new(
            10,
            MatrixAdjustRule::DontAdjust,
            0,
            0,
            query_source.len() as i32,
            4,
            8,
            0,
            None,
        );
        second.next = Some(Box::new(BlastCompoAlignment::new(
            10,
            MatrixAdjustRule::DontAdjust,
            0,
            0,
            query_source.len() as i32,
            0,
            4,
            0,
            None,
        )));
        let incoming = Some(Box::new(second));
        let mut params = redo_one_match_test_params(CompoAdjustMode::NoCompositionBasedStats);
        params.gapping_params.gap_open = 11;
        params.gapping_params.gap_extend = 1;
        params.cutoff_s = 1;
        params.do_link_hsps = true;
        assert_eq!(
            s_matrix_info_init(&mut params.matrix_info, "BLOSUM62", 0.3176, 1.0),
            0
        );
        let matching_seq = matching_sequence_initialize(subject_source.len() as i32, 7, 0);
        let mut alignments = vec![None];
        let heaps = vec![BlastCompoHeap::new(10, 10.0)];

        let status = blast_redo_one_match_smith_waterman_in_memory_protein(
            &mut alignments,
            &params,
            &incoming,
            2,
            0.267,
            0.041_f64.ln(),
            &matching_seq,
            &query_info,
            crate::program::BLASTP,
            &subject_source,
            &heaps,
            &crate::util::STANDARD_GENETIC_CODE,
        );

        assert_eq!(status, 0);
        assert!(alignment_list_len(alignments[0].as_deref()) > 2);
    }

    #[test]
    fn blast_redo_one_match_smith_waterman_in_memory_nucl_rejects_adjusted_path() {
        let query_source = vec![0u8, 1, 2, 3];
        let query_info = vec![BlastCompoQueryInfo {
            origin: 0,
            seq: BlastCompoSequenceData {
                buffer: query_source.clone(),
                data_offset: 0,
                length: query_source.len() as i32,
            },
            composition: BlastAminoAcidComposition::default(),
            eff_search_space: 100.0,
            words: Vec::new(),
        }];
        let incoming = Some(Box::new(BlastCompoAlignment::new(
            10,
            MatrixAdjustRule::DontAdjust,
            0,
            0,
            4,
            0,
            4,
            0,
            None,
        )));
        let params = redo_one_match_test_params(CompoAdjustMode::CompositionBasedStats);
        let matching_seq = matching_sequence_initialize(query_source.len() as i32, 7, 0);
        let mut alignments = vec![Some(Box::new(BlastCompoAlignment::new(
            9,
            MatrixAdjustRule::DontAdjust,
            0,
            0,
            4,
            0,
            4,
            0,
            None,
        )))];
        let heaps = vec![BlastCompoHeap::new(10, 10.0)];
        let matrix = crate::matrix::nucleotide_matrix(1, -3);

        let status = blast_redo_one_match_smith_waterman_in_memory_nucl(
            &mut alignments,
            &params,
            &incoming,
            1,
            0.267,
            0.041_f64.ln(),
            &matching_seq,
            &query_info,
            &query_source,
            &heaps,
            &matrix,
        );

        assert_eq!(status, -1);
        assert!(alignments[0].is_none());
    }

    #[test]
    fn smith_waterman_significance_link_hsps_uses_raw_score_cutoff() {
        let mut params = redo_one_match_test_params(CompoAdjustMode::NoCompositionBasedStats);
        params.do_link_hsps = true;
        params.cutoff_s = 50;
        let heap = BlastCompoHeap::new(1, 1e-10);

        assert!(smith_waterman_alignment_is_significant(
            50, 0.267, -3.0, 1000.0, &params, None, &heap, 7
        ));
        assert!(!smith_waterman_alignment_is_significant(
            49, 0.267, -3.0, 1000.0, &params, None, &heap, 7
        ));
    }

    #[test]
    fn smith_waterman_significance_first_alignment_checks_heap() {
        let mut params = redo_one_match_test_params(CompoAdjustMode::NoCompositionBasedStats);
        params.cutoff_e = 1.0;
        let mut heap = BlastCompoHeap::new(1, 1e-10);
        assert!(heap.insert(heap_hsp_list(1, 1e-10, 100)).is_none());

        // E-value is below cutoff_e, but the full heap would reject it, so the
        // first redone alignment for this query is not significant.
        assert!(!smith_waterman_alignment_is_significant(
            10,
            0.267,
            0.041_f64.ln(),
            1.0,
            &params,
            None,
            &heap,
            7,
        ));

        // Once a query already has an alignment, NCBI skips the heap gate.
        let existing =
            BlastCompoAlignment::new(100, MatrixAdjustRule::DontAdjust, 0, 0, 10, 0, 10, 0, None);
        assert!(smith_waterman_alignment_is_significant(
            10,
            0.267,
            0.041_f64.ln(),
            1.0,
            &params,
            Some(&existing),
            &heap,
            7,
        ));
    }

    #[test]
    fn smith_waterman_significance_non_link_uses_strict_evalue_cutoff() {
        let mut params = redo_one_match_test_params(CompoAdjustMode::NoCompositionBasedStats);
        let heap = BlastCompoHeap::new(10, 1.0);
        params.cutoff_e = evalue_from_score(20, 0.267, 0.041_f64.ln(), 1.0);

        assert!(!smith_waterman_alignment_is_significant(
            20,
            0.267,
            0.041_f64.ln(),
            1.0,
            &params,
            None,
            &heap,
            7,
        ));
        assert!(smith_waterman_alignment_is_significant(
            21,
            0.267,
            0.041_f64.ln(),
            1.0,
            &params,
            None,
            &heap,
            7,
        ));
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
    fn get_start_freq_ratios_blosum62_aliases_match_ncbi_scaling() {
        let canonical = crate::matrix::get_blosum62_freq_ratios();
        let b20 = get_start_freq_ratios("BLOSUM62_20").expect("BLOSUM62_20 ratios");
        let b20a = get_start_freq_ratios("BLOSUM62_20A").expect("BLOSUM62_20A ratios");
        let b20b = get_start_freq_ratios("BLOSUM62_20B").expect("BLOSUM62_20B ratios");

        assert_eq!(b20[1][1], canonical[1][1]);
        assert_eq!(b20a[1][1], canonical[1][1] * 0.9666);
        assert_eq!(b20b[1][1], canonical[1][1] * 0.9344);
    }

    #[test]
    fn get_start_freq_ratios_standard_matrices_are_ported() {
        let cases = [
            ("BLOSUM45", 2.95043377),
            ("BLOSUM50", 3.27354473),
            ("BLOSUM80", 4.77313697),
            ("BLOSUM90", 5.49812903),
            ("PAM30", 7.78912912),
            ("PAM70", 4.89972946),
            ("PAM250", 1.51578006),
        ];
        for (name, aa_ratio) in cases {
            let ratios = get_start_freq_ratios(name).expect(name);
            assert_eq!(ratios[1][1], aa_ratio);
        }
    }

    #[test]
    fn get_pos_based_start_freq_ratios_overlays_start_numerator() {
        let query = [1u8, 20u8]; // A, W in NCBIstdaa.
        let mut start_numerator = vec![vec![0.0; crate::matrix::AA_SIZE]; query.len()];
        start_numerator[0][3] = 0.25; // C column; standard C prob > epsilon.
        start_numerator[1][crate::encoding::NCBISTDAA_STOP as usize] = 99.0; // Stop column must not be overlaid.
        let ratios = get_pos_based_start_freq_ratios(&query, "BLOSUM62", &start_numerator)
            .expect("position ratios");
        let canonical = crate::matrix::get_blosum62_freq_ratios();
        let std_prob = crate::stat::protein_std_freq_ncbistdaa();

        assert_eq!(ratios[0][1], canonical[1][1]);
        assert_eq!(ratios[0][3], 0.25 / std_prob[3]);
        assert_eq!(
            ratios[1][crate::encoding::NCBISTDAA_STOP as usize],
            canonical[20][crate::encoding::NCBISTDAA_STOP as usize]
        );
    }

    #[test]
    fn get_pos_based_start_freq_ratios_rejects_short_numerator() {
        let query = [1u8, 20u8];
        let start_numerator = vec![vec![0.0; crate::matrix::AA_SIZE]];
        assert!(get_pos_based_start_freq_ratios(&query, "BLOSUM62", &start_numerator).is_err());
    }

    #[test]
    fn pos_freq_ratios_to_pssm_converts_log_odds_rows() {
        let ratios = vec![crate::matrix::get_blosum62_freq_ratios()[1].to_vec()];
        let pssm = pos_freq_ratios_to_pssm(&ratios, 0.3176).expect("pssm rows");
        assert_eq!(pssm.len(), 1);
        assert_eq!(pssm[0][1], 4);
        assert_eq!(pssm[0][3], 0);
        assert_eq!(pssm[0][4], -2);
    }

    #[test]
    fn pos_freq_ratios_to_scaled_pssm_uses_psi_private_scale() {
        let ratios = vec![crate::matrix::get_blosum62_freq_ratios()[1].to_vec()];
        let scaled = pos_freq_ratios_to_scaled_pssm(&ratios).expect("scaled pssm rows");
        assert_eq!(scaled.len(), 1);
        assert_eq!(scaled[0][1], 272);
        assert_eq!(scaled[0][3], -28);
        assert_eq!(scaled[0][4], -122);
    }

    #[test]
    fn scale_pos_matrix_populates_public_and_private_rows() {
        let ratios = vec![crate::matrix::get_blosum62_freq_ratios()[1].to_vec()];
        let mut info = BlastMatrixInfo {
            ungapped_lambda: 0.3176,
            ..Default::default()
        };
        let rc = scale_pos_matrix(&mut info, &ratios);
        assert_eq!(rc, 0);
        assert_eq!(info.matrix.len(), 1);
        assert_eq!(info.scaled_matrix.len(), 1);
        assert_eq!(info.matrix[0][1], 4);
        assert_eq!(info.scaled_matrix[0][1], 272);
    }

    #[test]
    fn scale_pos_matrix_rejects_invalid_lambda_or_short_rows() {
        let ratios = vec![crate::matrix::get_blosum62_freq_ratios()[1].to_vec()];
        let mut bad_lambda = BlastMatrixInfo::default();
        assert_eq!(scale_pos_matrix(&mut bad_lambda, &ratios), -1);

        let mut info = BlastMatrixInfo {
            ungapped_lambda: 0.3176,
            ..Default::default()
        };
        assert_eq!(scale_pos_matrix(&mut info, &[vec![1.0; 2]]), -1);
    }

    #[test]
    fn psi_private_update_lambda_statistics_applies_ratio_to_valid_blocks() {
        let mut kbp = vec![
            crate::stat::KarlinBlk {
                lambda: 0.267,
                k: 0.041,
                log_k: 0.041_f64.ln(),
                h: 0.14,
                round_down: true,
            },
            crate::stat::KarlinBlk {
                lambda: -1.0,
                k: -1.0,
                log_k: 0.0,
                h: -1.0,
                round_down: false,
            },
        ];

        assert_eq!(psi_private_update_lambda_statistics(&mut kbp, 0.75), 0);
        assert!((kbp[0].lambda - 0.267 * 0.75).abs() < 1e-12);
        assert_eq!(kbp[0].k, 0.041);
        assert_eq!(kbp[0].log_k, 0.041_f64.ln());
        assert_eq!(kbp[1].lambda, -1.0);
        assert_eq!(psi_private_update_lambda_statistics(&mut kbp, 0.0), -1);
        assert_eq!(psi_private_update_lambda_statistics(&mut kbp, f64::NAN), -1);
    }

    #[test]
    fn position_based_scaling_ratio_updates_private_lambda_statistics() {
        let query = vec![1u8, 20u8];
        let mut info = BlastMatrixInfo::default();
        let start_numerator = vec![vec![0.0; crate::matrix::AA_SIZE]; query.len()];
        assert_eq!(
            matrix_info_init_psiblast_from_start_numerator(
                &mut info,
                &query,
                "BLOSUM62",
                &start_numerator,
                0.3176,
                1.0,
            ),
            0
        );
        let mut query_composition = vec![0.0f64; crate::matrix::AA_SIZE];
        let mut subject_composition = vec![0.0f64; crate::matrix::AA_SIZE];
        for &idx in &crate::composition::TRUE_CHAR_POSITIONS {
            query_composition[idx] = 1.0 / crate::composition::COMPO_NUM_TRUE_AA as f64;
            subject_composition[idx] = 1.0 / crate::composition::COMPO_NUM_TRUE_AA as f64;
        }
        subject_composition[crate::composition::TRUE_CHAR_POSITIONS[0]] += 0.19;
        for &idx in &crate::composition::TRUE_CHAR_POSITIONS[1..] {
            subject_composition[idx] -= 0.01;
        }
        let mut adjusted = Vec::new();
        let mut pvalue = -1.0;
        let mut ratio = 0.0;
        assert_eq!(
            blast_adjust_position_based_scores(
                &mut adjusted,
                &query_composition,
                crate::composition::COMPO_NUM_TRUE_AA,
                &subject_composition,
                crate::composition::COMPO_NUM_TRUE_AA,
                &info,
                1,
                &mut pvalue,
                &mut ratio,
            ),
            Ok(MatrixAdjustRule::ScaleOldMatrix)
        );

        let mut kbp = vec![crate::stat::KarlinBlk {
            lambda: 0.267,
            k: 0.041,
            log_k: 0.041_f64.ln(),
            h: 0.14,
            round_down: true,
        }];
        assert_eq!(psi_private_update_lambda_statistics(&mut kbp, ratio), 0);
        assert!((0.5..=1.0).contains(&ratio));
        assert!((kbp[0].lambda - 0.267 * ratio).abs() < 1e-12);
    }

    fn kappa_test_score_block() -> crate::stat::BlastScoreBlk {
        let mut sbp = crate::stat::blast_score_blk_new(crate::encoding::BLASTAA_SEQ_CODE, 1)
            .expect("protein score block");
        sbp.matrix.data = crate::matrix::BLOSUM62
            .iter()
            .map(|row| row.to_vec())
            .collect();
        sbp.matrix.nrows = crate::encoding::BLASTAA_SIZE;
        sbp.matrix.ncols = crate::encoding::BLASTAA_SIZE;
        sbp.name = Some("BLOSUM62".to_string());
        sbp.loscore = -4;
        sbp.hiscore = 11;
        sbp.kbp_ideal = Some(crate::stat::KarlinBlk {
            lambda: 0.3176,
            k: 0.134,
            log_k: 0.134_f64.ln(),
            h: 0.401,
            round_down: false,
        });
        sbp.kbp_psi[0] = crate::stat::KarlinBlk {
            lambda: 0.267,
            k: 0.041,
            log_k: 0.041_f64.ln(),
            h: 0.14,
            round_down: false,
        };
        sbp
    }

    #[test]
    fn kappa_search_item_lifecycle_matches_c_shape() {
        let private = vec![vec![K_PSI_SCALE_FACTOR as i32; crate::encoding::BLASTAA_SIZE]; 2];
        let freqs = vec![vec![0.0; crate::encoding::BLASTAA_SIZE]; 2];
        let mut pos = Some(
            kappa_pos_search_items_new(2, "BLOSUM62", private, freqs).expect("pos search items"),
        );
        assert_eq!(pos.as_ref().unwrap().pos_matrix.len(), 2);
        assert_eq!(
            pos.as_ref().unwrap().std_freq_ratios.len(),
            crate::matrix::AA_SIZE
        );
        assert!(kappa_pos_search_items_free(&mut pos).is_none());
        assert!(pos.is_none());

        let sbp = kappa_test_score_block();
        let mut compact =
            Some(kappa_compact_search_items_new(&[1, 3], 2, &sbp).expect("compact search items"));
        assert_eq!(
            compact.as_ref().unwrap().standard_prob.len(),
            crate::encoding::BLASTAA_SIZE
        );
        assert_eq!(compact.as_ref().unwrap().lambda_ideal, 0.3176);
        assert!(kappa_compact_search_items_free(&mut compact).is_none());
        assert!(compact.is_none());
    }

    #[test]
    fn fill_sfp_matches_blast_posit_true_residue_distribution() {
        let mut row1 = vec![crate::stat::BLAST_SCORE_MIN; crate::encoding::BLASTAA_SIZE];
        let mut row2 = row1.clone();
        row1[crate::encoding::NCBISTDAA_A as usize] = 2;
        row1[crate::encoding::NCBISTDAA_C as usize] = -1;
        row2[crate::encoding::NCBISTDAA_A as usize] = 4;
        row2[crate::encoding::NCBISTDAA_C as usize] = -1;
        row2[crate::encoding::NCBISTDAA_X as usize] = 99;
        let mut probs = vec![0.0; crate::encoding::BLASTAA_SIZE];
        probs[crate::encoding::NCBISTDAA_A as usize] = 0.25;
        probs[crate::encoding::NCBISTDAA_C as usize] = 0.75;

        let sfp = fill_sfp(&[row1, row2], 2, &probs).expect("score frequencies");
        assert_eq!(sfp.obs_min, -1);
        assert_eq!(sfp.obs_max, 4);
        let at = |score: i32| sfp.sprob[(score - sfp.score_min) as usize];
        assert!((at(-1) - 0.75).abs() < 1e-12);
        assert!((at(2) - 0.125).abs() < 1e-12);
        assert!((at(4) - 0.125).abs() < 1e-12);
    }

    #[test]
    fn psi_update_lambda_k_recomputes_psi_and_gap_k() {
        let mut sbp = kappa_test_score_block();
        sbp.kbp_gap_std[0] = crate::stat::KarlinBlk {
            lambda: 0.267,
            k: 0.041,
            log_k: 0.041_f64.ln(),
            h: 0.14,
            round_down: false,
        };
        let query = vec![crate::encoding::NCBISTDAA_A, crate::encoding::NCBISTDAA_C];
        let mut pssm = vec![vec![-3; crate::encoding::BLASTAA_SIZE]; query.len()];
        pssm[0][crate::encoding::NCBISTDAA_A as usize] = 5;
        pssm[1][crate::encoding::NCBISTDAA_C as usize] = 5;
        let probs = crate::stat::protein_std_freq_ncbistdaa();

        assert_eq!(
            psi_update_lambda_k(&pssm, &query, query.len(), &probs, &mut sbp),
            0
        );
        assert!(sbp.kbp_psi[0].lambda > 0.0);
        assert!(sbp.kbp_psi[0].k > 0.0);
        let expected_gap_k =
            sbp.kbp_psi[0].k * sbp.kbp_gap_std[0].k / sbp.kbp_ideal.as_ref().unwrap().k;
        assert!((sbp.kbp_gap_psi[0].k - expected_gap_k).abs() < 1e-12);
        assert!((sbp.kbp_gap_psi[0].log_k - expected_gap_k.ln()).abs() < 1e-12);
    }

    #[test]
    fn kappa_impala_scaling_scales_public_and_private_pssm_rows() {
        let mut sbp = kappa_test_score_block();
        let private = vec![
            vec![K_PSI_SCALE_FACTOR as i32; crate::encoding::BLASTAA_SIZE],
            vec![-K_PSI_SCALE_FACTOR as i32; crate::encoding::BLASTAA_SIZE],
        ];
        let freqs = vec![vec![0.0; crate::encoding::BLASTAA_SIZE]; 2];
        let mut pos =
            kappa_pos_search_items_new(2, "BLOSUM62", private, freqs).expect("pos search items");
        let compact = kappa_compact_search_items_new(
            &[crate::encoding::NCBISTDAA_A, crate::encoding::NCBISTDAA_C],
            2,
            &sbp,
        )
        .expect("compact search items");

        assert_eq!(
            kappa_impala_scaling(&mut pos, &compact, 32.0, false, &mut sbp),
            0
        );
        assert_eq!(pos.pos_matrix[0][crate::encoding::NCBISTDAA_A as usize], 1);
        assert_eq!(pos.pos_matrix[1][crate::encoding::NCBISTDAA_A as usize], -1);
        assert_eq!(
            pos.pos_private_matrix[0][crate::encoding::NCBISTDAA_A as usize],
            32
        );
        assert_eq!(
            pos.pos_private_matrix[1][crate::encoding::NCBISTDAA_A as usize],
            -32
        );
    }

    #[test]
    fn matrix_info_init_psiblast_from_start_numerator_builds_position_rows() {
        let query = [1u8, 20u8];
        let mut start_numerator = vec![vec![0.0; crate::matrix::AA_SIZE]; query.len()];
        start_numerator[0][3] = 0.25;
        let mut info = BlastMatrixInfo::default();
        let rc = matrix_info_init_psiblast_from_start_numerator(
            &mut info,
            &query,
            "BLOSUM62",
            &start_numerator,
            0.3176,
            1.0,
        );
        assert_eq!(rc, 0);
        assert!(info.positional);
        assert_eq!(info.rows, 2);
        assert_eq!(info.cols, crate::matrix::AA_SIZE as i32);
        assert_eq!(info.bit_scale_factor, 2);
        assert_eq!(
            info.start_freq_ratios[0][3],
            0.25 / crate::stat::protein_std_freq_ncbistdaa()[3]
        );
        assert_eq!(info.matrix[1][20], 11);
        assert_eq!(info.scaled_matrix[1][20], 728);
    }

    #[test]
    fn get_start_freq_ratios_unknown_matrix_errs() {
        assert!(get_start_freq_ratios("UNKNOWN").is_err());
    }

    #[test]
    fn sw_find_final_ends_using_xdrop_returns_target_score_on_full_match() {
        // Identity match should produce the maximum score with the
        // first X-drop pass. Use BLOSUM62 + ACGT-style residues.
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
            &q,
            0,
            q.len() - 1,
            &s,
            0,
            s.len() - 1,
            &m,
            5,
            2,
            20,
            /* target */ 40,
        );
        assert!(score >= 40);
        assert!(q_ext > 0);
        assert!(s_ext > 0);
    }

    #[test]
    fn blast_gap_align_workspace_clone_is_deep_copy() {
        let mut script = crate::gapinfo::GapEditScript::new();
        script.push(crate::gapinfo::GapAlignOpType::Sub, 5);
        let mut ws = BlastGapAlignWorkspace {
            gap_x_dropoff: 65,
            score: 100,
            query_start: 0,
            query_stop: 50,
            subject_start: 10,
            subject_stop: 60,
            edit_script: Some(script),
            ..BlastGapAlignWorkspace::default()
        };
        let copy = ws.s_blast_gap_align_struct_copy();
        assert_eq!(copy.score, 100);
        assert_eq!(copy.gap_x_dropoff, 65);
        assert_eq!(
            copy.edit_script.as_ref().unwrap().len(),
            ws.edit_script.as_ref().unwrap().len()
        );
        ws.edit_script
            .as_mut()
            .unwrap()
            .push(crate::gapinfo::GapAlignOpType::Ins, 2);
        assert_eq!(copy.edit_script.as_ref().unwrap().len(), 1);
        assert_eq!(ws.edit_script.as_ref().unwrap().len(), 2);
    }

    #[test]
    fn s_blast_gap_align_struct_free_clears_slot() {
        let mut workspace = blast_gap_align_struct_new(42).expect("workspace");
        workspace.score = 100;
        workspace.query_start = 3;
        workspace.query_stop = 12;
        workspace.subject_start = 8;
        workspace.subject_stop = 17;
        workspace.edit_script = crate::gapinfo::gap_edit_script_new(2);
        workspace
            .edit_script
            .as_mut()
            .expect("edit script")
            .push(crate::gapinfo::GapAlignOpType::Sub, 9);
        let mut slot = Some(workspace);
        s_blast_gap_align_struct_free(&mut slot);
        assert!(slot.is_none());

        let workspace = blast_gap_align_struct_new(-5).expect("workspace");
        assert_eq!(workspace.gap_x_dropoff, 0);
        assert_eq!(workspace.score, 0);
        assert!(workspace.edit_script.is_none());

        let mut internal_slot = blast_gap_align_struct_new(7);
        s_blast_gap_align_struct_free(&mut internal_slot);
        assert!(internal_slot.is_none());
    }

    #[test]
    fn blast_gap_align_struct_new_initializes_dp_workspace_from_full_params() {
        let scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
                matrix_path: None,
                reward: 1,
                penalty: -3,
                gap_open: 5,
                gap_extend: 2,
                shift_pen: i16::MAX as i32,
                gapped_calculation: true,
                complexity_adjusted_scoring: false,
                matrix_name: None,
                is_ooframe: false,
                program_number: crate::program::BLASTN,
            },
            1.0,
        );
        let mut options = crate::options::ExtensionOptions::new_blastn();
        options.prelim_gap_ext = crate::options::PrelimGapExt::DynProgScoreOnly;
        options.max_mismatches = 4;
        options.mismatch_window = 9;
        let extension = crate::parameters::ExtensionParameters {
            options,
            gap_x_dropoff: 33,
            gap_x_dropoff_final: 80,
            gap_trigger: 22,
        };

        let workspace =
            blast_gap_align_struct_new((&scoring, &extension, 1234, true)).expect("workspace");

        assert!(workspace.position_based);
        assert_eq!(workspace.gap_x_dropoff, 33);
        assert_eq!(workspace.max_mismatches, 4);
        assert_eq!(workspace.mismatch_window, 9);
        assert_eq!(workspace.dp_mem_alloc, 1000);
        assert!(workspace.fwd_prelim_tback.is_some());
        assert!(workspace.rev_prelim_tback.is_some());
        assert!(workspace.greedy_align_mem.is_none());
    }

    #[test]
    fn blast_gap_align_struct_new_initializes_greedy_and_chaining_workspace() {
        let scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
                matrix_path: None,
                reward: 1,
                penalty: -3,
                gap_open: 5,
                gap_extend: 2,
                shift_pen: i16::MAX as i32,
                gapped_calculation: true,
                complexity_adjusted_scoring: false,
                matrix_name: None,
                is_ooframe: false,
                program_number: crate::program::BLASTN,
            },
            1.0,
        );
        let mut options = crate::options::ExtensionOptions::new_blastn();
        options.prelim_gap_ext = crate::options::PrelimGapExt::GreedyScoreOnly;
        options.chaining = true;
        let extension = crate::parameters::ExtensionParameters {
            options,
            gap_x_dropoff: 25,
            gap_x_dropoff_final: 80,
            gap_trigger: 22,
        };

        let workspace =
            blast_gap_align_struct_new((&scoring, &extension, 5_000, false)).expect("workspace");

        assert_eq!(workspace.gap_x_dropoff, 25);
        assert!(workspace.greedy_align_mem.is_some());
        assert!(workspace.chaining.is_some());
        assert_eq!(workspace.dp_mem_alloc, 0);
        assert!(workspace.jumper.is_none());
    }

    #[test]
    fn blast_gap_align_struct_new_initializes_jumper_zero_xdrop_like_c() {
        let scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
                matrix_path: None,
                reward: 1,
                penalty: -3,
                gap_open: 5,
                gap_extend: 2,
                shift_pen: i16::MAX as i32,
                gapped_calculation: true,
                complexity_adjusted_scoring: false,
                matrix_name: None,
                is_ooframe: false,
                program_number: crate::program::BLASTN,
            },
            1.0,
        );
        let mut options = crate::options::ExtensionOptions::new_blastn();
        options.prelim_gap_ext = crate::options::PrelimGapExt::JumperWithTraceback;
        let extension = crate::parameters::ExtensionParameters {
            options,
            gap_x_dropoff: 0,
            gap_x_dropoff_final: 80,
            gap_trigger: 22,
        };

        let workspace =
            blast_gap_align_struct_new((&scoring, &extension, 0, false)).expect("workspace");

        assert_eq!(workspace.gap_x_dropoff, 21);
        assert!(workspace.jumper.is_some());
        assert!(workspace.greedy_align_mem.is_none());
        assert_eq!(workspace.dp_mem_alloc, 0);
    }

    #[test]
    fn blast_score_blk_copy_is_deep_copy_for_modeled_fields() {
        let src = BlastScoreBlkSnapshot {
            matrix_name: "BLOSUM62".to_string(),
            matrix: vec![vec![1, 2], vec![3, 4]],
            psi_matrix: Some(vec![vec![5, 6, 7]]),
            kbp: vec![crate::stat::KarlinBlk {
                lambda: 0.3,
                k: 0.04,
                log_k: 0.04_f64.ln(),
                h: 0.14,
                round_down: false,
            }],
            kbp_gap: vec![crate::stat::KarlinBlk {
                lambda: 0.267,
                k: 0.041,
                log_k: 0.041_f64.ln(),
                h: 0.14,
                round_down: true,
            }],
            scale_factor: 32.0,
            round_down: true,
        };

        let mut copy = blast_score_blk_copy(&src);
        copy.matrix[0][0] = 99;
        copy.psi_matrix.as_mut().unwrap()[0][0] = 88;
        copy.kbp[0].lambda = 1.0;
        copy.kbp_gap[0].round_down = false;

        assert_eq!(src.matrix[0][0], 1);
        assert_eq!(src.psi_matrix.as_ref().unwrap()[0][0], 5);
        assert_eq!(src.kbp[0].lambda, 0.3);
        assert!(src.kbp_gap[0].round_down);
        assert_eq!(copy.matrix_name, "BLOSUM62");
        assert_eq!(copy.scale_factor, 32.0);
        assert!(copy.round_down);
    }

    #[test]
    fn blast_redo_alignment_core_delegates_to_mt() {
        // No-composition path: the single-thread entry delegates to the MT driver,
        // which preserves the incoming HSP list in per-query results.
        let qi = crate::queryinfo::QueryInfo {
            num_queries: 1,
            contexts: vec![crate::queryinfo::ContextInfo {
                query_offset: 0,
                query_length: 10,
                eff_searchsp: 100,
                length_adjustment: 0,
                query_index: 0,
                frame: 0,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: 0,
            min_length: 0,
        };
        let mut kbp = vec![];
        let mut mtx = Vec::<Vec<i32>>::new();
        let mut scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
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
                context: None,
            },
            CompoAdjustMode::NoCompositionBasedStats,
            1.0,
            false,
            false,
            false,
            0,
            0,
            0.0,
            false,
            0.0,
        );
        let mut saved = BlastKappaSavedParameters::s_saved_parameters_new(
            0,
            0,
            CompoAdjustMode::NoCompositionBasedStats,
            false,
        );
        let mut hsp_list = HspList::new(7);
        hsp_list.add_hsp(Hsp {
            score: 42,
            num_ident: 10,
            bit_score: 12.0,
            evalue: 1.0e-6,
            query_offset: 0,
            query_end: 10,
            query_gapped_start: 0,
            subject_offset: 3,
            subject_end: 13,
            subject_gapped_start: 3,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
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
            BlastRedoAlignmentSource::ExistingMatchOnly,
            &mut results,
        );
        assert_eq!(rc, 0);
        let hitlist = results.hitlists[0].as_ref().expect("query hitlist");
        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(hitlist.hsp_lists[0].oid, 7);
        assert_eq!(hitlist.hsp_lists[0].hsps[0].score, 42);
    }

    #[test]
    fn blast_redo_alignment_core_mt_restores_after_scaled_composition_setup() {
        let qi = crate::queryinfo::QueryInfo {
            num_queries: 1,
            contexts: vec![crate::queryinfo::ContextInfo {
                query_offset: 0,
                query_length: 10,
                eff_searchsp: 100,
                length_adjustment: 0,
                query_index: 0,
                frame: 0,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: 10,
            min_length: 0,
        };
        let original_kbp = crate::stat::KarlinBlk {
            lambda: 0.267,
            k: 0.041,
            log_k: 0.041_f64.ln(),
            h: 0.14,
            round_down: true,
        };
        let mut kbp = vec![original_kbp.clone()];
        let mut mtx = vec![vec![0i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE];
        mtx[1][1] = 4;
        let original_matrix_value = mtx[1][1];
        let mut scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
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
            },
            1.0,
        );
        let original_gap_open = scoring.gap_open;
        let original_gap_extend = scoring.gap_extend;
        let original_scale_factor = scoring.scale_factor;
        let params = blast_redo_align_params_new(
            BlastMatrixInfo::default(),
            BlastCompoGappingParams {
                gap_open: 352,
                gap_extend: 32,
                decline_align: i32::MIN,
                x_dropoff: 2078,
                context: None,
            },
            CompoAdjustMode::CompositionMatrixAdjust,
            SCALING_FACTOR,
            false,
            false,
            false,
            0,
            0,
            10.0,
            false,
            0.0,
        );
        let mut saved = BlastKappaSavedParameters::s_saved_parameters_new(
            0,
            1,
            CompoAdjustMode::CompositionMatrixAdjust,
            false,
        );
        let mut hsp_list = HspList::new(7);
        let mut results = crate::hspstream::HspResults::new(1);

        let rc = blast_redo_alignment_core_mt(
            crate::program::BLASTP,
            1,
            &[],
            &qi,
            &mut kbp,
            &mut mtx,
            &mut scoring,
            &params,
            &mut saved,
            &mut hsp_list,
            BlastRedoAlignmentSource::ExistingMatchOnly,
            &mut results,
        );

        assert_eq!(rc, -1);
        assert_eq!(scoring.gap_open, original_gap_open);
        assert_eq!(scoring.gap_extend, original_gap_extend);
        assert_eq!(scoring.scale_factor, original_scale_factor);
        assert_eq!(kbp[0].lambda, original_kbp.lambda);
        assert_eq!(kbp[0].log_k, original_kbp.log_k);
        assert_eq!(mtx[1][1], original_matrix_value);
    }

    #[test]
    fn blast_redo_alignment_core_mt_rejects_blosum62_20_without_composition_before_rescale() {
        let qi = crate::queryinfo::QueryInfo {
            num_queries: 1,
            contexts: vec![crate::queryinfo::ContextInfo {
                query_offset: 0,
                query_length: 10,
                eff_searchsp: 100,
                length_adjustment: 0,
                query_index: 0,
                frame: 0,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: 10,
            min_length: 0,
        };
        let original_kbp = crate::stat::KarlinBlk {
            lambda: 0.267,
            k: 0.041,
            log_k: 0.041_f64.ln(),
            h: 0.14,
            round_down: true,
        };
        let mut kbp = vec![original_kbp.clone()];
        let mut mtx = vec![vec![0i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE];
        mtx[1][1] = 4;
        let mut scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
                matrix_path: None,
                reward: 0,
                penalty: 0,
                gap_open: 11,
                gap_extend: 1,
                shift_pen: i16::MAX as i32,
                gapped_calculation: true,
                complexity_adjusted_scoring: false,
                matrix_name: Some("BLOSUM62_20".to_string()),
                is_ooframe: false,
                program_number: crate::program::UNDEFINED,
            },
            1.0,
        );
        let original_gap_open = scoring.gap_open;
        let params = blast_redo_align_params_new(
            BlastMatrixInfo::default(),
            BlastCompoGappingParams {
                gap_open: 352,
                gap_extend: 32,
                decline_align: i32::MIN,
                x_dropoff: 2078,
                context: None,
            },
            CompoAdjustMode::NoCompositionBasedStats,
            SCALING_FACTOR,
            false,
            false,
            false,
            0,
            0,
            10.0,
            false,
            0.0,
        );
        let mut saved = BlastKappaSavedParameters::s_saved_parameters_new(
            0,
            1,
            CompoAdjustMode::NoCompositionBasedStats,
            false,
        );
        let mut hsp_list = HspList::new(7);
        let mut results = crate::hspstream::HspResults::new(1);

        let rc = blast_redo_alignment_core_mt(
            crate::program::BLASTP,
            1,
            &[],
            &qi,
            &mut kbp,
            &mut mtx,
            &mut scoring,
            &params,
            &mut saved,
            &mut hsp_list,
            BlastRedoAlignmentSource::ExistingMatchOnly,
            &mut results,
        );

        assert_eq!(rc, -1);
        assert_eq!(scoring.gap_open, original_gap_open);
        assert_eq!(kbp[0].lambda, original_kbp.lambda);
        assert_eq!(mtx[1][1], 4);
    }

    #[test]
    fn blast_redo_alignment_core_mt_preflight_clamps_position_based_adjustment_mode() {
        let qi = crate::queryinfo::QueryInfo {
            num_queries: 1,
            contexts: vec![crate::queryinfo::ContextInfo {
                query_offset: 0,
                query_length: 3,
                eff_searchsp: 9,
                length_adjustment: 0,
                query_index: 0,
                frame: 0,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: 3,
            min_length: 0,
        };
        let scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
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
                context: None,
            },
            CompoAdjustMode::CompositionMatrixAdjust,
            SCALING_FACTOR,
            true,
            false,
            false,
            0,
            0,
            10.0,
            false,
            0.0,
        );
        let matrix = vec![vec![0i32; crate::matrix::AA_SIZE]; 3];

        let effective = blast_redo_alignment_core_mt_effective_params(
            &[1, 2, 3],
            &qi,
            &scoring,
            &matrix,
            &params,
        )
        .expect("effective params");

        assert_eq!(
            effective.compo_adjust_mode,
            CompoAdjustMode::CompositionBasedStats
        );
    }

    #[test]
    fn blast_redo_alignment_core_mt_with_callbacks_redoes_position_based_matches() {
        let query = vec![1u8; 8];
        let qi = crate::queryinfo::QueryInfo {
            num_queries: 1,
            contexts: vec![crate::queryinfo::ContextInfo {
                query_offset: 0,
                query_length: query.len() as i32,
                eff_searchsp: 100,
                length_adjustment: 0,
                query_index: 0,
                frame: 0,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: query.len() as u32,
            min_length: 0,
        };
        let original_kbp = crate::stat::KarlinBlk {
            lambda: 0.267,
            k: 0.041,
            log_k: 0.041_f64.ln(),
            h: 0.14,
            round_down: true,
        };
        let mut matrix_info = BlastMatrixInfo::default();
        let start_numerator = vec![vec![0.0; crate::matrix::AA_SIZE]; query.len()];
        assert_eq!(
            matrix_info_init_psiblast_from_start_numerator(
                &mut matrix_info,
                &query,
                "BLOSUM62",
                &start_numerator,
                0.3176,
                1.0,
            ),
            0
        );
        let mut params = blast_redo_align_params_new(
            matrix_info.clone(),
            BlastCompoGappingParams {
                gap_open: 11,
                gap_extend: 1,
                decline_align: i32::MIN,
                x_dropoff: 30,
                context: None,
            },
            CompoAdjustMode::CompositionBasedStats,
            1.0,
            true,
            false,
            false,
            query.len() as i32,
            1,
            10.0,
            false,
            0.0,
        );
        params.callbacks = Some(BlastRedoAlignCallbacks {
            get_range: Some(callback_test_get_range_alanine),
            ..Default::default()
        });
        let make_hsp_list = || {
            let mut list = HspList::new(7);
            list.add_hsp(Hsp {
                score: 10,
                num_ident: 0,
                bit_score: 0.0,
                evalue: 1.0,
                query_offset: 0,
                query_end: query.len() as i32,
                query_gapped_start: 0,
                subject_offset: 0,
                subject_end: query.len() as i32,
                subject_gapped_start: 0,
                context: 0,
                query_frame: 0,
                subject_frame: 0,
                num_gaps: 0,
                comp_adjustment_method: 0,
                edit_script: None,
                pat_info: None,
                map_info: None,
            });
            list
        };
        let run = |smith_waterman: bool| {
            let mut kbp = vec![original_kbp.clone()];
            let mut mtx = matrix_info.matrix.clone();
            let original_mtx = mtx.clone();
            let mut scoring = crate::parameters::ScoringParameters::from_options(
                &crate::options::ScoringOptions {
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
                },
                1.0,
            );
            let mut saved = BlastKappaSavedParameters::s_saved_parameters_new(
                query.len() as i32,
                1,
                CompoAdjustMode::CompositionBasedStats,
                true,
            );
            let mut hsp_list = make_hsp_list();
            let mut results = crate::hspstream::HspResults::new(1);

            let rc = blast_redo_alignment_core_mt_with_callbacks(
                crate::program::PSI_BLAST,
                1,
                &query,
                &qi,
                &mut kbp,
                &mut mtx,
                &mut scoring,
                &params,
                &mut saved,
                &mut hsp_list,
                BlastRedoCallbackSubjectConfig {
                    subject_length: query.len() as i32,
                    smith_waterman,
                    expect_value: 10.0,
                    hitlist_size: 10,
                    inclusion_ethresh: 10.0,
                    link_context: None,
                },
                &mut results,
            );

            assert_eq!(rc, 0);
            assert_eq!(kbp[0].lambda, original_kbp.lambda);
            assert_eq!(mtx, original_mtx);
            assert_eq!(
                hsp_list.hsps[0].comp_adjustment_method,
                CompoAdjustMode::CompositionBasedStats as i32
            );
            assert_eq!(hsp_list.hsps[0].num_ident, query.len() as i32);
            assert!(hsp_list.hsps[0].edit_script.is_some());
            let hitlist = results.hitlists[0].as_ref().expect("query hitlist");
            assert_eq!(hitlist.hsp_lists.len(), 1);
            assert_eq!(hitlist.hsp_lists[0].oid, 7);
            assert_eq!(
                hitlist.hsp_lists[0].hsps[0].comp_adjustment_method,
                CompoAdjustMode::CompositionBasedStats as i32
            );
            assert_eq!(hitlist.hsp_lists[0].hsps[0].num_ident, query.len() as i32);
        };

        run(false);
        run(true);
    }

    #[test]
    fn blast_redo_alignment_core_mt_with_callbacks_uses_link_hsp_context() {
        let query = vec![0u8; 8];
        let qi = crate::queryinfo::QueryInfo {
            num_queries: 1,
            contexts: vec![crate::queryinfo::ContextInfo {
                query_offset: 0,
                query_length: query.len() as i32,
                eff_searchsp: 100,
                length_adjustment: 0,
                query_index: 0,
                frame: 0,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: query.len() as u32,
            min_length: 0,
        };
        let karlin = crate::stat::KarlinBlk {
            lambda: 0.267,
            k: 0.041,
            log_k: 0.041_f64.ln(),
            h: 0.14,
            round_down: true,
        };
        let mut kbp = vec![karlin.clone()];
        let score_block = crate::link_hsps::LinkScoreBlock {
            kbp: vec![karlin],
            kbp_gap: Vec::new(),
            ..Default::default()
        };
        let link_params = crate::link_hsps::LinkHSPParameters::default();
        let link_context = HitlistLinkContext {
            query_info: &qi,
            query_context: 0,
            score_block: &score_block,
            link_params: &link_params,
            gapped_calculation: false,
        };
        let mut mtx = vec![vec![0i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE];
        let mut scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
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
            },
            1.0,
        );
        let mut matrix_info = BlastMatrixInfo::default();
        assert_eq!(
            s_matrix_info_init(&mut matrix_info, "BLOSUM62", 0.3176, 1.0),
            0
        );
        let mut params = blast_redo_align_params_new(
            matrix_info,
            BlastCompoGappingParams {
                gap_open: 11,
                gap_extend: 1,
                decline_align: i32::MIN,
                x_dropoff: 30,
                context: None,
            },
            CompoAdjustMode::NoCompositionBasedStats,
            1.0,
            false,
            false,
            false,
            query.len() as i32,
            1,
            10.0,
            true,
            0.0,
        );
        params.callbacks = Some(BlastRedoAlignCallbacks {
            get_range: Some(callback_test_get_range),
            redo_one_alignment: Some(callback_test_redo_one_alignment),
            ..Default::default()
        });
        let mut saved = BlastKappaSavedParameters::s_saved_parameters_new(
            0,
            1,
            CompoAdjustMode::NoCompositionBasedStats,
            false,
        );
        let mut hsp_list = HspList::new(7);
        hsp_list.add_hsp(Hsp {
            score: 10,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1.0,
            query_offset: 0,
            query_end: query.len() as i32,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: query.len() as i32,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        let mut results = crate::hspstream::HspResults::new(1);

        let rc = blast_redo_alignment_core_mt_with_callbacks(
            crate::program::BLASTP,
            1,
            &query,
            &qi,
            &mut kbp,
            &mut mtx,
            &mut scoring,
            &params,
            &mut saved,
            &mut hsp_list,
            BlastRedoCallbackSubjectConfig {
                subject_length: query.len() as i32,
                smith_waterman: false,
                expect_value: 10.0,
                hitlist_size: 10,
                inclusion_ethresh: 10.0,
                link_context: Some(&link_context),
            },
            &mut results,
        );

        assert_eq!(rc, 0);
        assert_eq!(hsp_list.hsps[0].num_ident, query.len() as i32);
        let hitlist = results.hitlists[0].as_ref().expect("query hitlist");
        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(hitlist.hsp_lists[0].oid, 7);
        assert_eq!(hitlist.hsp_lists[0].hsps[0].num_ident, query.len() as i32);
    }

    #[test]
    fn blast_redo_alignment_core_mt_in_memory_subject_redoes_single_match() {
        let query = vec![0u8, 1, 2, 3, 0, 1, 2, 3, 0, 1];
        let subject = query.clone();
        let qi = crate::queryinfo::QueryInfo {
            num_queries: 1,
            contexts: vec![crate::queryinfo::ContextInfo {
                query_offset: 0,
                query_length: query.len() as i32,
                eff_searchsp: 100,
                length_adjustment: 0,
                query_index: 0,
                frame: 0,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: query.len() as u32,
            min_length: 0,
        };
        let original_kbp = crate::stat::KarlinBlk {
            lambda: 0.267,
            k: 0.041,
            log_k: 0.041_f64.ln(),
            h: 0.14,
            round_down: true,
        };
        let mut kbp = vec![original_kbp.clone()];
        let mut mtx = vec![vec![0i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE];
        let mut scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
                matrix_path: None,
                reward: 1,
                penalty: -3,
                gap_open: 5,
                gap_extend: 2,
                shift_pen: i16::MAX as i32,
                gapped_calculation: true,
                complexity_adjusted_scoring: false,
                matrix_name: None,
                is_ooframe: false,
                program_number: crate::program::UNDEFINED,
            },
            1.0,
        );
        let params = blast_redo_align_params_new(
            BlastMatrixInfo::default(),
            BlastCompoGappingParams {
                gap_open: 5,
                gap_extend: 2,
                decline_align: i32::MIN,
                x_dropoff: 30,
                context: None,
            },
            CompoAdjustMode::NoCompositionBasedStats,
            1.0,
            false,
            false,
            false,
            query.len() as i32,
            1,
            10.0,
            false,
            0.0,
        );
        let mut saved = BlastKappaSavedParameters::s_saved_parameters_new(
            0,
            1,
            CompoAdjustMode::NoCompositionBasedStats,
            false,
        );
        let mut hsp_list = HspList::new(7);
        hsp_list.add_hsp(Hsp {
            score: 10,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1.0,
            query_offset: 2,
            query_end: 6,
            query_gapped_start: 4,
            subject_offset: 2,
            subject_end: 6,
            subject_gapped_start: 4,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        let mut results = crate::hspstream::HspResults::new(1);

        let rc = blast_redo_alignment_core_mt_in_memory_subject(
            crate::program::BLASTN,
            1,
            &query,
            &qi,
            &mut kbp,
            &mut mtx,
            &mut scoring,
            &params,
            &mut saved,
            &mut hsp_list,
            &mut results,
            BlastRedoInMemorySubject {
                subject_source: &subject,
                reward: 1,
                penalty: -3,
                genetic_code: &crate::util::STANDARD_GENETIC_CODE,
                smith_waterman: false,
                expect_value: 10.0,
                hitlist_size: 10,
                inclusion_ethresh: 10.0,
                link_context: None,
            },
        );

        assert_eq!(rc, 0);
        assert_eq!(kbp[0].lambda, original_kbp.lambda);
        assert_eq!(hsp_list.oid, 7);
        assert_eq!(hsp_list.hsps.len(), 1);
        assert_eq!(hsp_list.hsps[0].query_offset, 0);
        assert_eq!(hsp_list.hsps[0].subject_offset, 0);
        assert_eq!(hsp_list.hsps[0].num_ident, query.len() as i32);
        let hitlist = results.hitlists[0].as_ref().expect("query hitlist");
        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(hitlist.hsp_lists[0].hsps[0].query_end, query.len() as i32);
    }

    #[test]
    fn blast_redo_alignment_core_mt_in_memory_subject_redoes_smith_waterman_match() {
        let query = vec![0u8, 1, 2, 3, 0, 1, 2, 3, 0, 1];
        let subject = query.clone();
        let qi = crate::queryinfo::QueryInfo {
            num_queries: 1,
            contexts: vec![crate::queryinfo::ContextInfo {
                query_offset: 0,
                query_length: query.len() as i32,
                eff_searchsp: 100,
                length_adjustment: 0,
                query_index: 0,
                frame: 0,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: query.len() as u32,
            min_length: 0,
        };
        let original_kbp = crate::stat::KarlinBlk {
            lambda: 0.267,
            k: 0.041,
            log_k: 0.041_f64.ln(),
            h: 0.14,
            round_down: true,
        };
        let mut kbp = vec![original_kbp.clone()];
        let mut mtx = vec![vec![0i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE];
        let mut scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
                matrix_path: None,
                reward: 1,
                penalty: -3,
                gap_open: 5,
                gap_extend: 2,
                shift_pen: i16::MAX as i32,
                gapped_calculation: true,
                complexity_adjusted_scoring: false,
                matrix_name: None,
                is_ooframe: false,
                program_number: crate::program::UNDEFINED,
            },
            1.0,
        );
        let params = blast_redo_align_params_new(
            BlastMatrixInfo::default(),
            BlastCompoGappingParams {
                gap_open: 5,
                gap_extend: 2,
                decline_align: i32::MIN,
                x_dropoff: 30,
                context: None,
            },
            CompoAdjustMode::NoCompositionBasedStats,
            1.0,
            false,
            false,
            false,
            query.len() as i32,
            1,
            10.0,
            false,
            0.0,
        );
        let mut saved = BlastKappaSavedParameters::s_saved_parameters_new(
            0,
            1,
            CompoAdjustMode::NoCompositionBasedStats,
            false,
        );
        let mut hsp_list = HspList::new(7);
        hsp_list.add_hsp(Hsp {
            score: 10,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1.0,
            query_offset: 2,
            query_end: 6,
            query_gapped_start: 4,
            subject_offset: 2,
            subject_end: 6,
            subject_gapped_start: 4,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        let mut results = crate::hspstream::HspResults::new(1);

        let rc = blast_redo_alignment_core_mt_in_memory_subject(
            crate::program::BLASTN,
            1,
            &query,
            &qi,
            &mut kbp,
            &mut mtx,
            &mut scoring,
            &params,
            &mut saved,
            &mut hsp_list,
            &mut results,
            BlastRedoInMemorySubject {
                subject_source: &subject,
                reward: 1,
                penalty: -3,
                genetic_code: &crate::util::STANDARD_GENETIC_CODE,
                smith_waterman: true,
                expect_value: 10.0,
                hitlist_size: 10,
                inclusion_ethresh: 10.0,
                link_context: None,
            },
        );

        assert_eq!(rc, 0);
        assert_eq!(kbp[0].lambda, original_kbp.lambda);
        assert_eq!(hsp_list.hsps.len(), 1);
        assert_eq!(hsp_list.hsps[0].query_offset, 0);
        assert_eq!(hsp_list.hsps[0].subject_offset, 0);
        assert_eq!(hsp_list.hsps[0].query_end, query.len() as i32);
        assert!(hsp_list.hsps[0].edit_script.is_some());
        let hitlist = results.hitlists[0].as_ref().expect("query hitlist");
        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(hitlist.hsp_lists[0].hsps[0].num_ident, query.len() as i32);
    }

    #[test]
    fn blast_redo_alignment_core_mt_rejects_blastn_composition_smith_waterman() {
        let query = vec![0u8, 1, 2, 3, 0, 1, 2, 3, 0, 1];
        let subject = query.clone();
        let qi = crate::queryinfo::QueryInfo {
            num_queries: 1,
            contexts: vec![crate::queryinfo::ContextInfo {
                query_offset: 0,
                query_length: query.len() as i32,
                eff_searchsp: 100,
                length_adjustment: 0,
                query_index: 0,
                frame: 0,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: query.len() as u32,
            min_length: 0,
        };
        let original_kbp = crate::stat::KarlinBlk {
            lambda: 0.267,
            k: 0.041,
            log_k: 0.041_f64.ln(),
            h: 0.14,
            round_down: true,
        };
        let make_hsp_list = || {
            let mut list = HspList::new(7);
            list.add_hsp(Hsp {
                score: 10,
                num_ident: 0,
                bit_score: 0.0,
                evalue: 1.0,
                query_offset: 2,
                query_end: 6,
                query_gapped_start: 4,
                subject_offset: 2,
                subject_end: 6,
                subject_gapped_start: 4,
                context: 0,
                query_frame: 0,
                subject_frame: 0,
                num_gaps: 0,
                comp_adjustment_method: 0,
                edit_script: None,
                pat_info: None,
                map_info: None,
            });
            list
        };
        let mut kbp = vec![original_kbp.clone()];
        let mut mtx = vec![vec![0i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE];
        let mut scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
                matrix_path: None,
                reward: 1,
                penalty: -3,
                gap_open: 5,
                gap_extend: 2,
                shift_pen: i16::MAX as i32,
                gapped_calculation: true,
                complexity_adjusted_scoring: false,
                matrix_name: None,
                is_ooframe: false,
                program_number: crate::program::UNDEFINED,
            },
            1.0,
        );
        let params = blast_redo_align_params_new(
            BlastMatrixInfo::default(),
            BlastCompoGappingParams {
                gap_open: 5,
                gap_extend: 2,
                decline_align: i32::MIN,
                x_dropoff: 30,
                context: None,
            },
            CompoAdjustMode::CompositionBasedStats,
            1.0,
            false,
            false,
            false,
            query.len() as i32,
            1,
            10.0,
            false,
            0.0,
        );
        let mut saved = BlastKappaSavedParameters::s_saved_parameters_new(
            crate::matrix::AA_SIZE as i32,
            1,
            CompoAdjustMode::CompositionBasedStats,
            false,
        );
        let mut hsp_list = make_hsp_list();
        let original_hsp_list = hsp_list.clone();
        let mut results = crate::hspstream::HspResults::new(1);

        let rc = blast_redo_alignment_core_mt_in_memory_subject(
            crate::program::BLASTN,
            1,
            &query,
            &qi,
            &mut kbp,
            &mut mtx,
            &mut scoring,
            &params,
            &mut saved,
            &mut hsp_list,
            &mut results,
            BlastRedoInMemorySubject {
                subject_source: &subject,
                reward: 1,
                penalty: -3,
                genetic_code: &crate::util::STANDARD_GENETIC_CODE,
                smith_waterman: true,
                expect_value: 10.0,
                hitlist_size: 10,
                inclusion_ethresh: 10.0,
                link_context: None,
            },
        );

        assert_eq!(rc, -1);
        assert_eq!(kbp[0].lambda, original_kbp.lambda);
        assert_eq!(
            hsp_list.hsps[0].query_offset,
            original_hsp_list.hsps[0].query_offset
        );
        assert!(results.hitlists[0].is_none());

        let mut stream_matches = vec![BlastRedoInMemorySubjectMatch {
            hsp_list: make_hsp_list(),
            subject: BlastRedoInMemorySubject {
                subject_source: &subject,
                reward: 1,
                penalty: -3,
                genetic_code: &crate::util::STANDARD_GENETIC_CODE,
                smith_waterman: true,
                expect_value: 10.0,
                hitlist_size: 10,
                inclusion_ethresh: 10.0,
                link_context: None,
            },
        }];
        let mut saved_stream = BlastKappaSavedParameters::s_saved_parameters_new(
            crate::matrix::AA_SIZE as i32,
            1,
            CompoAdjustMode::CompositionBasedStats,
            false,
        );
        let mut results_stream = crate::hspstream::HspResults::new(1);

        let rc_stream = blast_redo_alignment_core_mt_in_memory_subjects(
            crate::program::BLASTN,
            2,
            &query,
            &qi,
            &mut kbp,
            &mut mtx,
            &mut scoring,
            &params,
            &mut saved_stream,
            &mut stream_matches,
            &mut results_stream,
        );

        assert_eq!(rc_stream, -1);
        assert_eq!(kbp[0].lambda, original_kbp.lambda);
        assert_eq!(stream_matches[0].hsp_list.hsps[0].query_offset, 2);
        assert!(results_stream.hitlists[0].is_none());
    }

    #[test]
    fn blast_redo_alignment_core_mt_in_memory_subject_redoes_protein_match() {
        let query = vec![1u8; 8];
        let subject = query.clone();
        let qi = crate::queryinfo::QueryInfo {
            num_queries: 1,
            contexts: vec![crate::queryinfo::ContextInfo {
                query_offset: 0,
                query_length: query.len() as i32,
                eff_searchsp: 100,
                length_adjustment: 0,
                query_index: 0,
                frame: 0,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: query.len() as u32,
            min_length: 0,
        };
        let mut kbp = vec![crate::stat::KarlinBlk {
            lambda: 0.267,
            k: 0.041,
            log_k: 0.041_f64.ln(),
            h: 0.14,
            round_down: true,
        }];
        let mut mtx = vec![vec![0i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE];
        let mut scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
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
            },
            1.0,
        );
        let mut matrix_info = BlastMatrixInfo::default();
        assert_eq!(
            s_matrix_info_init(&mut matrix_info, "BLOSUM62", 0.3176, 1.0),
            0
        );
        let params = blast_redo_align_params_new(
            matrix_info,
            BlastCompoGappingParams {
                gap_open: 11,
                gap_extend: 1,
                decline_align: i32::MIN,
                x_dropoff: 30,
                context: None,
            },
            CompoAdjustMode::NoCompositionBasedStats,
            1.0,
            false,
            false,
            false,
            query.len() as i32,
            1,
            10.0,
            false,
            0.0,
        );
        let mut saved = BlastKappaSavedParameters::s_saved_parameters_new(
            0,
            1,
            CompoAdjustMode::NoCompositionBasedStats,
            false,
        );
        let mut hsp_list = HspList::new(7);
        hsp_list.add_hsp(Hsp {
            score: 10,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1.0,
            query_offset: 0,
            query_end: 8,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: 8,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        let mut results = crate::hspstream::HspResults::new(1);

        let rc = blast_redo_alignment_core_mt_in_memory_subject(
            crate::program::BLASTP,
            1,
            &query,
            &qi,
            &mut kbp,
            &mut mtx,
            &mut scoring,
            &params,
            &mut saved,
            &mut hsp_list,
            &mut results,
            BlastRedoInMemorySubject {
                subject_source: &subject,
                reward: 0,
                penalty: 0,
                genetic_code: &crate::util::STANDARD_GENETIC_CODE,
                smith_waterman: false,
                expect_value: 10.0,
                hitlist_size: 10,
                inclusion_ethresh: 10.0,
                link_context: None,
            },
        );

        assert_eq!(rc, 0);
        assert_eq!(hsp_list.oid, 7);
        assert_eq!(hsp_list.hsps.len(), 1);
        assert_eq!(hsp_list.hsps[0].query_offset, 0);
        assert_eq!(hsp_list.hsps[0].subject_offset, 0);
        assert_eq!(hsp_list.hsps[0].num_ident, query.len() as i32);
        assert!(hsp_list.hsps[0].edit_script.is_some());
        let hitlist = results.hitlists[0].as_ref().expect("query hitlist");
        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(hitlist.hsp_lists[0].hsps[0].num_ident, query.len() as i32);
    }

    #[test]
    fn blast_redo_alignment_core_mt_in_memory_subject_redoes_protein_composition_adjusted_match() {
        let query = vec![1u8, 2, 3, 4, 1, 2, 3, 4];
        let subject = query.clone();
        let qi = crate::queryinfo::QueryInfo {
            num_queries: 1,
            contexts: vec![crate::queryinfo::ContextInfo {
                query_offset: 0,
                query_length: query.len() as i32,
                eff_searchsp: 100,
                length_adjustment: 0,
                query_index: 0,
                frame: 0,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: query.len() as u32,
            min_length: 0,
        };
        let original_kbp = crate::stat::KarlinBlk {
            lambda: 0.267,
            k: 0.041,
            log_k: 0.041_f64.ln(),
            h: 0.14,
            round_down: true,
        };
        let mut kbp = vec![original_kbp.clone()];
        let mut matrix_info = BlastMatrixInfo::default();
        assert_eq!(
            s_matrix_info_init(&mut matrix_info, "BLOSUM62", 0.3176, 1.0),
            0
        );
        let mut mtx = matrix_info.matrix.clone();
        let original_mtx = mtx.clone();
        let mut scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
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
            },
            1.0,
        );
        let params = blast_redo_align_params_new(
            matrix_info,
            BlastCompoGappingParams {
                gap_open: 11,
                gap_extend: 1,
                decline_align: i32::MIN,
                x_dropoff: 30,
                context: None,
            },
            CompoAdjustMode::CompositionBasedStats,
            1.0,
            false,
            false,
            false,
            query.len() as i32,
            1,
            10.0,
            false,
            0.0,
        );
        let mut saved = BlastKappaSavedParameters::s_saved_parameters_new(
            crate::matrix::AA_SIZE as i32,
            1,
            CompoAdjustMode::CompositionBasedStats,
            false,
        );
        let mut hsp_list = HspList::new(7);
        hsp_list.add_hsp(Hsp {
            score: 10,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1.0,
            query_offset: 0,
            query_end: query.len() as i32,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: subject.len() as i32,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        let mut results = crate::hspstream::HspResults::new(1);

        let rc = blast_redo_alignment_core_mt_in_memory_subject(
            crate::program::BLASTP,
            1,
            &query,
            &qi,
            &mut kbp,
            &mut mtx,
            &mut scoring,
            &params,
            &mut saved,
            &mut hsp_list,
            &mut results,
            BlastRedoInMemorySubject {
                subject_source: &subject,
                reward: 0,
                penalty: 0,
                genetic_code: &crate::util::STANDARD_GENETIC_CODE,
                smith_waterman: false,
                expect_value: 10.0,
                hitlist_size: 10,
                inclusion_ethresh: 10.0,
                link_context: None,
            },
        );

        assert_eq!(rc, 0);
        assert_eq!(kbp[0].lambda, original_kbp.lambda);
        assert_eq!(mtx, original_mtx);
        assert_eq!(hsp_list.hsps.len(), 1);
        assert_eq!(hsp_list.hsps[0].num_ident, query.len() as i32);
        assert_eq!(
            hsp_list.hsps[0].comp_adjustment_method,
            CompoAdjustMode::CompositionBasedStats as i32
        );
        assert!(hsp_list.hsps[0].edit_script.is_some());
        let hitlist = results.hitlists[0].as_ref().expect("query hitlist");
        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(
            hitlist.hsp_lists[0].hsps[0].comp_adjustment_method,
            CompoAdjustMode::CompositionBasedStats as i32
        );

        let mut hsp_list_sw = HspList::new(7);
        hsp_list_sw.add_hsp(Hsp {
            score: 10,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1.0,
            query_offset: 0,
            query_end: query.len() as i32,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: subject.len() as i32,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        let mut matches_sw = vec![BlastRedoInMemorySubjectMatch {
            hsp_list: hsp_list_sw,
            subject: BlastRedoInMemorySubject {
                subject_source: &subject,
                reward: 0,
                penalty: 0,
                genetic_code: &crate::util::STANDARD_GENETIC_CODE,
                smith_waterman: true,
                expect_value: 10.0,
                hitlist_size: 10,
                inclusion_ethresh: 10.0,
                link_context: None,
            },
        }];
        let mut saved_sw = BlastKappaSavedParameters::s_saved_parameters_new(
            crate::matrix::AA_SIZE as i32,
            1,
            CompoAdjustMode::CompositionBasedStats,
            false,
        );
        let mut results_sw = crate::hspstream::HspResults::new(1);

        let rc_sw = blast_redo_alignment_core_mt_in_memory_subjects(
            crate::program::BLASTP,
            2,
            &query,
            &qi,
            &mut kbp,
            &mut mtx,
            &mut scoring,
            &params,
            &mut saved_sw,
            &mut matches_sw,
            &mut results_sw,
        );

        assert_eq!(rc_sw, 0);
        assert_eq!(kbp[0].lambda, original_kbp.lambda);
        assert_eq!(mtx, original_mtx);
        assert_eq!(matches_sw[0].hsp_list.hsps[0].num_ident, query.len() as i32);
        assert_eq!(
            matches_sw[0].hsp_list.hsps[0].comp_adjustment_method,
            CompoAdjustMode::CompositionBasedStats as i32
        );
        let hitlist_sw = results_sw.hitlists[0].as_ref().expect("query hitlist");
        assert_eq!(hitlist_sw.hsp_lists.len(), 1);
        assert_eq!(hitlist_sw.hsp_lists[0].oid, 7);
        assert_eq!(
            hitlist_sw.hsp_lists[0].hsps[0].num_ident,
            query.len() as i32
        );
        assert_eq!(
            hitlist_sw.hsp_lists[0].hsps[0].comp_adjustment_method,
            CompoAdjustMode::CompositionBasedStats as i32
        );
    }

    #[test]
    fn blast_redo_alignment_core_mt_in_memory_subject_redoes_position_based_composition_adjusted_match(
    ) {
        let query = vec![1u8, 2, 3, 4, 1, 2, 3, 4];
        let subject = query.clone();
        let qi = crate::queryinfo::QueryInfo {
            num_queries: 1,
            contexts: vec![crate::queryinfo::ContextInfo {
                query_offset: 0,
                query_length: query.len() as i32,
                eff_searchsp: 100,
                length_adjustment: 0,
                query_index: 0,
                frame: 0,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: query.len() as u32,
            min_length: 0,
        };
        let original_kbp = crate::stat::KarlinBlk {
            lambda: 0.267,
            k: 0.041,
            log_k: 0.041_f64.ln(),
            h: 0.14,
            round_down: true,
        };
        let mut kbp = vec![original_kbp.clone()];
        let mut matrix_info = BlastMatrixInfo::default();
        let start_numerator = vec![vec![0.0; crate::matrix::AA_SIZE]; query.len()];
        assert_eq!(
            matrix_info_init_psiblast_from_start_numerator(
                &mut matrix_info,
                &query,
                "BLOSUM62",
                &start_numerator,
                0.3176,
                1.0,
            ),
            0
        );
        let mut mtx = matrix_info.matrix.clone();
        let original_mtx = mtx.clone();
        let mut scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
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
            },
            1.0,
        );
        let params = blast_redo_align_params_new(
            matrix_info,
            BlastCompoGappingParams {
                gap_open: 11,
                gap_extend: 1,
                decline_align: i32::MIN,
                x_dropoff: 30,
                context: None,
            },
            CompoAdjustMode::CompositionBasedStats,
            1.0,
            true,
            false,
            false,
            query.len() as i32,
            1,
            10.0,
            false,
            0.0,
        );
        let mut saved = BlastKappaSavedParameters::s_saved_parameters_new(
            query.len() as i32,
            1,
            CompoAdjustMode::CompositionBasedStats,
            true,
        );
        let mut hsp_list = HspList::new(7);
        hsp_list.add_hsp(Hsp {
            score: 10,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1.0,
            query_offset: 0,
            query_end: query.len() as i32,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: subject.len() as i32,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        let mut results = crate::hspstream::HspResults::new(1);

        let rc = blast_redo_alignment_core_mt_in_memory_subject(
            crate::program::PSI_BLAST,
            1,
            &query,
            &qi,
            &mut kbp,
            &mut mtx,
            &mut scoring,
            &params,
            &mut saved,
            &mut hsp_list,
            &mut results,
            BlastRedoInMemorySubject {
                subject_source: &subject,
                reward: 0,
                penalty: 0,
                genetic_code: &crate::util::STANDARD_GENETIC_CODE,
                smith_waterman: false,
                expect_value: 10.0,
                hitlist_size: 10,
                inclusion_ethresh: 10.0,
                link_context: None,
            },
        );

        assert_eq!(rc, 0);
        assert_eq!(kbp[0].lambda, original_kbp.lambda);
        assert_eq!(mtx, original_mtx);
        assert_eq!(hsp_list.hsps[0].num_ident, query.len() as i32);
        assert_eq!(
            hsp_list.hsps[0].comp_adjustment_method,
            CompoAdjustMode::CompositionBasedStats as i32
        );
        assert!(hsp_list.hsps[0].edit_script.is_some());
        let hitlist = results.hitlists[0].as_ref().expect("query hitlist");
        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(hitlist.hsp_lists[0].oid, 7);
        assert_eq!(hitlist.hsp_lists[0].hsps[0].num_ident, query.len() as i32);
    }

    #[test]
    fn blast_redo_alignment_core_mt_in_memory_subject_redoes_position_based_smith_waterman_match() {
        let query = vec![1u8, 2, 3, 4, 1, 2, 3, 4];
        let subject = query.clone();
        let qi = crate::queryinfo::QueryInfo {
            num_queries: 1,
            contexts: vec![crate::queryinfo::ContextInfo {
                query_offset: 0,
                query_length: query.len() as i32,
                eff_searchsp: 100,
                length_adjustment: 0,
                query_index: 0,
                frame: 0,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: query.len() as u32,
            min_length: 0,
        };
        let original_kbp = crate::stat::KarlinBlk {
            lambda: 0.267,
            k: 0.041,
            log_k: 0.041_f64.ln(),
            h: 0.14,
            round_down: true,
        };
        let mut kbp = vec![original_kbp.clone()];
        let mut matrix_info = BlastMatrixInfo::default();
        let start_numerator = vec![vec![0.0; crate::matrix::AA_SIZE]; query.len()];
        assert_eq!(
            matrix_info_init_psiblast_from_start_numerator(
                &mut matrix_info,
                &query,
                "BLOSUM62",
                &start_numerator,
                0.3176,
                1.0,
            ),
            0
        );
        let mut mtx = matrix_info.matrix.clone();
        let original_mtx = mtx.clone();
        let mut scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
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
            },
            1.0,
        );
        let params = blast_redo_align_params_new(
            matrix_info,
            BlastCompoGappingParams {
                gap_open: 11,
                gap_extend: 1,
                decline_align: i32::MIN,
                x_dropoff: 30,
                context: None,
            },
            CompoAdjustMode::CompositionBasedStats,
            1.0,
            true,
            false,
            false,
            query.len() as i32,
            1,
            10.0,
            false,
            0.0,
        );
        let mut saved = BlastKappaSavedParameters::s_saved_parameters_new(
            query.len() as i32,
            1,
            CompoAdjustMode::CompositionBasedStats,
            true,
        );
        let mut hsp_list = HspList::new(7);
        hsp_list.add_hsp(Hsp {
            score: 10,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1.0,
            query_offset: 0,
            query_end: query.len() as i32,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: subject.len() as i32,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        let mut results = crate::hspstream::HspResults::new(1);

        let rc = blast_redo_alignment_core_mt_in_memory_subject(
            crate::program::PSI_BLAST,
            1,
            &query,
            &qi,
            &mut kbp,
            &mut mtx,
            &mut scoring,
            &params,
            &mut saved,
            &mut hsp_list,
            &mut results,
            BlastRedoInMemorySubject {
                subject_source: &subject,
                reward: 0,
                penalty: 0,
                genetic_code: &crate::util::STANDARD_GENETIC_CODE,
                smith_waterman: true,
                expect_value: 10.0,
                hitlist_size: 10,
                inclusion_ethresh: 10.0,
                link_context: None,
            },
        );

        assert_eq!(rc, 0);
        assert_eq!(kbp[0].lambda, original_kbp.lambda);
        assert_eq!(mtx, original_mtx);
        assert_eq!(hsp_list.hsps.len(), 1);
        assert_eq!(hsp_list.hsps[0].num_ident, query.len() as i32);
        assert_eq!(
            hsp_list.hsps[0].comp_adjustment_method,
            CompoAdjustMode::CompositionBasedStats as i32
        );
        assert!(hsp_list.hsps[0].edit_script.is_some());
        let hitlist = results.hitlists[0].as_ref().expect("query hitlist");
        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(hitlist.hsp_lists[0].hsps[0].num_ident, query.len() as i32);
    }

    #[test]
    fn blast_redo_alignment_core_mt_in_memory_subject_redoes_protein_composition_adjusted_smith_waterman_match(
    ) {
        let query = vec![1u8, 2, 3, 4, 1, 2, 3, 4];
        let subject = query.clone();
        let qi = crate::queryinfo::QueryInfo {
            num_queries: 1,
            contexts: vec![crate::queryinfo::ContextInfo {
                query_offset: 0,
                query_length: query.len() as i32,
                eff_searchsp: 100,
                length_adjustment: 0,
                query_index: 0,
                frame: 0,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: query.len() as u32,
            min_length: 0,
        };
        let original_kbp = crate::stat::KarlinBlk {
            lambda: 0.267,
            k: 0.041,
            log_k: 0.041_f64.ln(),
            h: 0.14,
            round_down: true,
        };
        let mut kbp = vec![original_kbp.clone()];
        let mut matrix_info = BlastMatrixInfo::default();
        assert_eq!(
            s_matrix_info_init(&mut matrix_info, "BLOSUM62", 0.3176, 1.0),
            0
        );
        let mut mtx = matrix_info.matrix.clone();
        let original_mtx = mtx.clone();
        let mut scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
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
            },
            1.0,
        );
        let params = blast_redo_align_params_new(
            matrix_info,
            BlastCompoGappingParams {
                gap_open: 11,
                gap_extend: 1,
                decline_align: i32::MIN,
                x_dropoff: 30,
                context: None,
            },
            CompoAdjustMode::CompositionBasedStats,
            1.0,
            false,
            false,
            false,
            query.len() as i32,
            1,
            10.0,
            false,
            0.0,
        );
        let mut saved = BlastKappaSavedParameters::s_saved_parameters_new(
            crate::matrix::AA_SIZE as i32,
            1,
            CompoAdjustMode::CompositionBasedStats,
            false,
        );
        let mut hsp_list = HspList::new(7);
        hsp_list.add_hsp(Hsp {
            score: 10,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1.0,
            query_offset: 0,
            query_end: query.len() as i32,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: subject.len() as i32,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        let mut results = crate::hspstream::HspResults::new(1);

        let rc = blast_redo_alignment_core_mt_in_memory_subject(
            crate::program::BLASTP,
            1,
            &query,
            &qi,
            &mut kbp,
            &mut mtx,
            &mut scoring,
            &params,
            &mut saved,
            &mut hsp_list,
            &mut results,
            BlastRedoInMemorySubject {
                subject_source: &subject,
                reward: 0,
                penalty: 0,
                genetic_code: &crate::util::STANDARD_GENETIC_CODE,
                smith_waterman: true,
                expect_value: 10.0,
                hitlist_size: 10,
                inclusion_ethresh: 10.0,
                link_context: None,
            },
        );

        assert_eq!(rc, 0);
        assert_eq!(kbp[0].lambda, original_kbp.lambda);
        assert_eq!(mtx, original_mtx);
        assert_eq!(hsp_list.hsps.len(), 1);
        assert_eq!(hsp_list.hsps[0].num_ident, query.len() as i32);
        assert_eq!(
            hsp_list.hsps[0].comp_adjustment_method,
            CompoAdjustMode::CompositionBasedStats as i32
        );
        assert!(hsp_list.hsps[0].edit_script.is_some());
        let hitlist = results.hitlists[0].as_ref().expect("query hitlist");
        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(hitlist.hsp_lists[0].hsps[0].num_ident, query.len() as i32);
    }

    #[test]
    fn blast_redo_alignment_core_mt_in_memory_subject_redoes_protein_smith_waterman_match() {
        let query = vec![1u8; 8];
        let subject = query.clone();
        let qi = crate::queryinfo::QueryInfo {
            num_queries: 1,
            contexts: vec![crate::queryinfo::ContextInfo {
                query_offset: 0,
                query_length: query.len() as i32,
                eff_searchsp: 100,
                length_adjustment: 0,
                query_index: 0,
                frame: 0,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: query.len() as u32,
            min_length: 0,
        };
        let mut kbp = vec![crate::stat::KarlinBlk {
            lambda: 0.267,
            k: 0.041,
            log_k: 0.041_f64.ln(),
            h: 0.14,
            round_down: true,
        }];
        let mut mtx = vec![vec![0i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE];
        let mut scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
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
            },
            1.0,
        );
        let mut matrix_info = BlastMatrixInfo::default();
        assert_eq!(
            s_matrix_info_init(&mut matrix_info, "BLOSUM62", 0.3176, 1.0),
            0
        );
        let params = blast_redo_align_params_new(
            matrix_info,
            BlastCompoGappingParams {
                gap_open: 11,
                gap_extend: 1,
                decline_align: i32::MIN,
                x_dropoff: 30,
                context: None,
            },
            CompoAdjustMode::NoCompositionBasedStats,
            1.0,
            false,
            false,
            false,
            query.len() as i32,
            1,
            10.0,
            false,
            0.0,
        );
        let mut saved = BlastKappaSavedParameters::s_saved_parameters_new(
            0,
            1,
            CompoAdjustMode::NoCompositionBasedStats,
            false,
        );
        let mut hsp_list = HspList::new(7);
        hsp_list.add_hsp(Hsp {
            score: 10,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1.0,
            query_offset: 0,
            query_end: 8,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: 8,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        let mut results = crate::hspstream::HspResults::new(1);

        let rc = blast_redo_alignment_core_mt_in_memory_subject(
            crate::program::BLASTP,
            1,
            &query,
            &qi,
            &mut kbp,
            &mut mtx,
            &mut scoring,
            &params,
            &mut saved,
            &mut hsp_list,
            &mut results,
            BlastRedoInMemorySubject {
                subject_source: &subject,
                reward: 0,
                penalty: 0,
                genetic_code: &crate::util::STANDARD_GENETIC_CODE,
                smith_waterman: true,
                expect_value: 10.0,
                hitlist_size: 10,
                inclusion_ethresh: 10.0,
                link_context: None,
            },
        );

        assert_eq!(rc, 0);
        assert_eq!(hsp_list.hsps.len(), 1);
        assert_eq!(hsp_list.hsps[0].num_ident, query.len() as i32);
        assert!(hsp_list.hsps[0].edit_script.is_some());
        let hitlist = results.hitlists[0].as_ref().expect("query hitlist");
        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(hitlist.hsp_lists[0].hsps[0].num_ident, query.len() as i32);
    }

    #[test]
    fn blast_redo_alignment_core_mt_in_memory_subject_redoes_tblastn_match() {
        let query = vec![12u8, 1, 12, 1];
        let subject = vec![1u8, 8, 4, 4, 2, 8, 1, 8, 4, 4, 2, 8];
        let qi = crate::queryinfo::QueryInfo {
            num_queries: 1,
            contexts: vec![crate::queryinfo::ContextInfo {
                query_offset: 0,
                query_length: query.len() as i32,
                eff_searchsp: 100,
                length_adjustment: 0,
                query_index: 0,
                frame: 0,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: query.len() as u32,
            min_length: 0,
        };
        let mut kbp = vec![crate::stat::KarlinBlk {
            lambda: 0.267,
            k: 0.041,
            log_k: 0.041_f64.ln(),
            h: 0.14,
            round_down: true,
        }];
        let mut mtx = vec![vec![0i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE];
        let mut scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
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
            },
            1.0,
        );
        let mut matrix_info = BlastMatrixInfo::default();
        assert_eq!(
            s_matrix_info_init(&mut matrix_info, "BLOSUM62", 0.3176, 1.0),
            0
        );
        let params = blast_redo_align_params_new(
            matrix_info,
            BlastCompoGappingParams {
                gap_open: 11,
                gap_extend: 1,
                decline_align: i32::MIN,
                x_dropoff: 30,
                context: None,
            },
            CompoAdjustMode::NoCompositionBasedStats,
            1.0,
            false,
            false,
            true,
            query.len() as i32,
            1,
            10.0,
            false,
            0.0,
        );
        let mut saved = BlastKappaSavedParameters::s_saved_parameters_new(
            0,
            1,
            CompoAdjustMode::NoCompositionBasedStats,
            false,
        );
        let mut hsp_list = HspList::new(7);
        hsp_list.add_hsp(Hsp {
            score: 10,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1.0,
            query_offset: 0,
            query_end: query.len() as i32,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: query.len() as i32,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 0,
            subject_frame: 1,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        let mut results = crate::hspstream::HspResults::new(1);

        let rc = blast_redo_alignment_core_mt_in_memory_subject(
            crate::program::TBLASTN,
            1,
            &query,
            &qi,
            &mut kbp,
            &mut mtx,
            &mut scoring,
            &params,
            &mut saved,
            &mut hsp_list,
            &mut results,
            BlastRedoInMemorySubject {
                subject_source: &subject,
                reward: 0,
                penalty: 0,
                genetic_code: &crate::util::STANDARD_GENETIC_CODE,
                smith_waterman: false,
                expect_value: 10.0,
                hitlist_size: 10,
                inclusion_ethresh: 10.0,
                link_context: None,
            },
        );

        assert_eq!(rc, 0);
        assert_eq!(hsp_list.hsps.len(), 1);
        assert_eq!(hsp_list.hsps[0].subject_frame, 1);
        assert_eq!(hsp_list.hsps[0].num_ident, query.len() as i32);
        assert!(hsp_list.hsps[0].edit_script.is_some());
        let hitlist = results.hitlists[0].as_ref().expect("query hitlist");
        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(hitlist.hsp_lists[0].hsps[0].num_ident, query.len() as i32);
    }

    #[test]
    fn blast_redo_alignment_core_mt_in_memory_subject_redoes_tblastn_smith_waterman_match() {
        let query = vec![12u8, 1, 12, 1];
        let subject = vec![1u8, 8, 4, 4, 2, 8, 1, 8, 4, 4, 2, 8];
        let qi = crate::queryinfo::QueryInfo {
            num_queries: 1,
            contexts: vec![crate::queryinfo::ContextInfo {
                query_offset: 0,
                query_length: query.len() as i32,
                eff_searchsp: 100,
                length_adjustment: 0,
                query_index: 0,
                frame: 0,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: query.len() as u32,
            min_length: 0,
        };
        let mut kbp = vec![crate::stat::KarlinBlk {
            lambda: 0.267,
            k: 0.041,
            log_k: 0.041_f64.ln(),
            h: 0.14,
            round_down: true,
        }];
        let mut mtx = vec![vec![0i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE];
        let mut scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
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
            },
            1.0,
        );
        let mut matrix_info = BlastMatrixInfo::default();
        assert_eq!(
            s_matrix_info_init(&mut matrix_info, "BLOSUM62", 0.3176, 1.0),
            0
        );
        let params = blast_redo_align_params_new(
            matrix_info,
            BlastCompoGappingParams {
                gap_open: 11,
                gap_extend: 1,
                decline_align: i32::MIN,
                x_dropoff: 30,
                context: None,
            },
            CompoAdjustMode::NoCompositionBasedStats,
            1.0,
            false,
            false,
            true,
            query.len() as i32,
            1,
            10.0,
            false,
            0.0,
        );
        let mut saved = BlastKappaSavedParameters::s_saved_parameters_new(
            0,
            1,
            CompoAdjustMode::NoCompositionBasedStats,
            false,
        );
        let mut hsp_list = HspList::new(7);
        hsp_list.add_hsp(Hsp {
            score: 10,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1.0,
            query_offset: 0,
            query_end: query.len() as i32,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: query.len() as i32,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 0,
            subject_frame: 1,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        let mut results = crate::hspstream::HspResults::new(1);

        let rc = blast_redo_alignment_core_mt_in_memory_subject(
            crate::program::TBLASTN,
            1,
            &query,
            &qi,
            &mut kbp,
            &mut mtx,
            &mut scoring,
            &params,
            &mut saved,
            &mut hsp_list,
            &mut results,
            BlastRedoInMemorySubject {
                subject_source: &subject,
                reward: 0,
                penalty: 0,
                genetic_code: &crate::util::STANDARD_GENETIC_CODE,
                smith_waterman: true,
                expect_value: 10.0,
                hitlist_size: 10,
                inclusion_ethresh: 10.0,
                link_context: None,
            },
        );

        assert_eq!(rc, 0);
        assert_eq!(hsp_list.hsps.len(), 1);
        assert_eq!(hsp_list.hsps[0].subject_frame, 1);
        assert_eq!(hsp_list.hsps[0].num_ident, query.len() as i32);
        assert!(hsp_list.hsps[0].edit_script.is_some());
        let hitlist = results.hitlists[0].as_ref().expect("query hitlist");
        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(hitlist.hsp_lists[0].hsps[0].num_ident, query.len() as i32);
    }

    #[test]
    fn blast_redo_alignment_core_mt_in_memory_subject_redoes_tblastn_composition_adjusted_smith_waterman_match(
    ) {
        let query = vec![12u8, 1, 12, 1];
        let subject = vec![1u8, 8, 4, 4, 2, 8, 1, 8, 4, 4, 2, 8];
        let qi = crate::queryinfo::QueryInfo {
            num_queries: 1,
            contexts: vec![crate::queryinfo::ContextInfo {
                query_offset: 0,
                query_length: query.len() as i32,
                eff_searchsp: 100,
                length_adjustment: 0,
                query_index: 0,
                frame: 0,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: query.len() as u32,
            min_length: 0,
        };
        let original_kbp = crate::stat::KarlinBlk {
            lambda: 0.267,
            k: 0.041,
            log_k: 0.041_f64.ln(),
            h: 0.14,
            round_down: true,
        };
        let mut kbp = vec![original_kbp.clone()];
        let mut matrix_info = BlastMatrixInfo::default();
        assert_eq!(
            s_matrix_info_init(&mut matrix_info, "BLOSUM62", 0.3176, 1.0),
            0
        );
        let mut mtx = matrix_info.matrix.clone();
        let original_mtx = mtx.clone();
        let mut scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
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
            },
            1.0,
        );
        let params = blast_redo_align_params_new(
            matrix_info,
            BlastCompoGappingParams {
                gap_open: 11,
                gap_extend: 1,
                decline_align: i32::MIN,
                x_dropoff: 30,
                context: None,
            },
            CompoAdjustMode::CompositionBasedStats,
            1.0,
            false,
            false,
            true,
            query.len() as i32,
            1,
            10.0,
            false,
            0.0,
        );
        let mut saved = BlastKappaSavedParameters::s_saved_parameters_new(
            crate::matrix::AA_SIZE as i32,
            1,
            CompoAdjustMode::CompositionBasedStats,
            false,
        );
        let mut hsp_list = HspList::new(7);
        hsp_list.add_hsp(Hsp {
            score: 10,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1.0,
            query_offset: 0,
            query_end: query.len() as i32,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: query.len() as i32,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 0,
            subject_frame: 1,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        let mut results = crate::hspstream::HspResults::new(1);

        let rc = blast_redo_alignment_core_mt_in_memory_subject(
            crate::program::TBLASTN,
            1,
            &query,
            &qi,
            &mut kbp,
            &mut mtx,
            &mut scoring,
            &params,
            &mut saved,
            &mut hsp_list,
            &mut results,
            BlastRedoInMemorySubject {
                subject_source: &subject,
                reward: 0,
                penalty: 0,
                genetic_code: &crate::util::STANDARD_GENETIC_CODE,
                smith_waterman: true,
                expect_value: 10.0,
                hitlist_size: 10,
                inclusion_ethresh: 10.0,
                link_context: None,
            },
        );

        assert_eq!(rc, 0);
        assert_eq!(kbp[0].lambda, original_kbp.lambda);
        assert_eq!(mtx, original_mtx);
        assert_eq!(hsp_list.hsps.len(), 1);
        assert_eq!(hsp_list.hsps[0].subject_frame, 1);
        assert_eq!(hsp_list.hsps[0].num_ident, query.len() as i32);
        assert_eq!(
            hsp_list.hsps[0].comp_adjustment_method,
            CompoAdjustMode::CompositionBasedStats as i32
        );
        assert!(hsp_list.hsps[0].edit_script.is_some());
        let hitlist = results.hitlists[0].as_ref().expect("query hitlist");
        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(hitlist.hsp_lists[0].hsps[0].num_ident, query.len() as i32);
        assert_eq!(
            hitlist.hsp_lists[0].hsps[0].comp_adjustment_method,
            CompoAdjustMode::CompositionBasedStats as i32
        );
    }

    #[test]
    fn blast_redo_alignment_core_mt_in_memory_subject_redoes_blastx_match() {
        let query_ncbi4na = vec![1u8, 8, 4, 4, 2, 8, 1, 8, 4, 4, 2, 8];
        let (query, offsets) = crate::util::blast_get_all_translations(
            &query_ncbi4na,
            query_ncbi4na.len(),
            &crate::util::STANDARD_GENETIC_CODE,
        );
        let mut qi = crate::queryinfo::QueryInfo::new_translated_query_from_offsets(&offsets);
        for ctx in &mut qi.contexts {
            ctx.eff_searchsp = 100;
        }
        let subject = vec![12u8, 1, 12, 1];
        let mut kbp = vec![
            crate::stat::KarlinBlk {
                lambda: 0.267,
                k: 0.041,
                log_k: 0.041_f64.ln(),
                h: 0.14,
                round_down: true,
            };
            qi.contexts.len()
        ];
        let mut mtx = vec![vec![0i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE];
        let mut scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
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
            },
            1.0,
        );
        let mut matrix_info = BlastMatrixInfo::default();
        assert_eq!(
            s_matrix_info_init(&mut matrix_info, "BLOSUM62", 0.3176, 1.0),
            0
        );
        let params = blast_redo_align_params_new(
            matrix_info,
            BlastCompoGappingParams {
                gap_open: 11,
                gap_extend: 1,
                decline_align: i32::MIN,
                x_dropoff: 30,
                context: None,
            },
            CompoAdjustMode::NoCompositionBasedStats,
            1.0,
            false,
            true,
            false,
            query_ncbi4na.len() as i32,
            1,
            10.0,
            false,
            0.0,
        );
        let mut saved = BlastKappaSavedParameters::s_saved_parameters_new(
            0,
            1,
            CompoAdjustMode::NoCompositionBasedStats,
            false,
        );
        let mut hsp_list = HspList::new(7);
        hsp_list.add_hsp(Hsp {
            score: 10,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1.0,
            query_offset: 0,
            query_end: subject.len() as i32,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: subject.len() as i32,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 1,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        let mut results = crate::hspstream::HspResults::new(1);

        let rc = blast_redo_alignment_core_mt_in_memory_subject(
            crate::program::BLASTX,
            1,
            &query,
            &qi,
            &mut kbp,
            &mut mtx,
            &mut scoring,
            &params,
            &mut saved,
            &mut hsp_list,
            &mut results,
            BlastRedoInMemorySubject {
                subject_source: &subject,
                reward: 0,
                penalty: 0,
                genetic_code: &crate::util::STANDARD_GENETIC_CODE,
                smith_waterman: false,
                expect_value: 10.0,
                hitlist_size: 10,
                inclusion_ethresh: 10.0,
                link_context: None,
            },
        );

        assert_eq!(rc, 0);
        assert_eq!(hsp_list.hsps.len(), 1);
        assert_eq!(hsp_list.hsps[0].query_frame, 1);
        assert_eq!(hsp_list.hsps[0].subject_frame, 0);
        assert_eq!(hsp_list.hsps[0].num_ident, subject.len() as i32);
        let hitlist = results.hitlists[0].as_ref().expect("query hitlist");
        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(hitlist.hsp_lists[0].hsps[0].num_ident, subject.len() as i32);

        let mut hsp_list = HspList::new(8);
        hsp_list.add_hsp(Hsp {
            score: 10,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1.0,
            query_offset: 0,
            query_end: subject.len() as i32,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: subject.len() as i32,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 1,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        let mut results = crate::hspstream::HspResults::new(1);
        let rc = blast_redo_alignment_core_mt_in_memory_subject(
            crate::program::BLASTX,
            1,
            &query,
            &qi,
            &mut kbp,
            &mut mtx,
            &mut scoring,
            &params,
            &mut saved,
            &mut hsp_list,
            &mut results,
            BlastRedoInMemorySubject {
                subject_source: &subject,
                reward: 0,
                penalty: 0,
                genetic_code: &crate::util::STANDARD_GENETIC_CODE,
                smith_waterman: true,
                expect_value: 10.0,
                hitlist_size: 10,
                inclusion_ethresh: 10.0,
                link_context: None,
            },
        );

        assert_eq!(rc, 0);
        assert_eq!(hsp_list.oid, 8);
        assert_eq!(hsp_list.hsps.len(), 1);
        assert_eq!(hsp_list.hsps[0].query_frame, 1);
        assert_eq!(hsp_list.hsps[0].num_ident, subject.len() as i32);
        assert!(hsp_list.hsps[0].edit_script.is_some());
    }

    #[test]
    fn blast_redo_alignment_core_mt_in_memory_subject_redoes_blastx_composition_adjusted_match() {
        let query_ncbi4na = vec![1u8, 8, 4, 4, 2, 8, 1, 8, 4, 4, 2, 8];
        let (query, offsets) = crate::util::blast_get_all_translations(
            &query_ncbi4na,
            query_ncbi4na.len(),
            &crate::util::STANDARD_GENETIC_CODE,
        );
        let mut qi = crate::queryinfo::QueryInfo::new_translated_query_from_offsets(&offsets);
        for ctx in &mut qi.contexts {
            ctx.eff_searchsp = 100;
        }
        let frame_len = qi.contexts[0].query_length;
        let frame_begin = qi.contexts[0].query_offset as usize;
        let subject = query[frame_begin..frame_begin + frame_len as usize].to_vec();
        let mut kbp = vec![
            crate::stat::KarlinBlk {
                lambda: 0.267,
                k: 0.041,
                log_k: 0.041_f64.ln(),
                h: 0.14,
                round_down: true,
            };
            qi.contexts.len()
        ];
        let mut matrix_info = BlastMatrixInfo::default();
        assert_eq!(
            s_matrix_info_init(&mut matrix_info, "BLOSUM62", 0.3176, 1.0),
            0
        );
        let mut mtx = matrix_info.matrix.clone();
        let mut scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
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
            },
            1.0,
        );
        let params = blast_redo_align_params_new(
            matrix_info,
            BlastCompoGappingParams {
                gap_open: 11,
                gap_extend: 1,
                decline_align: i32::MIN,
                x_dropoff: 30,
                context: None,
            },
            CompoAdjustMode::CompositionBasedStats,
            1.0,
            false,
            true,
            false,
            query_ncbi4na.len() as i32,
            1,
            10.0,
            false,
            0.0,
        );
        let mut saved = BlastKappaSavedParameters::s_saved_parameters_new(
            crate::matrix::AA_SIZE as i32,
            1,
            CompoAdjustMode::CompositionBasedStats,
            false,
        );
        let mut hsp_list = HspList::new(7);
        hsp_list.add_hsp(Hsp {
            score: 10,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1.0,
            query_offset: 0,
            query_end: frame_len,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: frame_len,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 1,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        let mut results = crate::hspstream::HspResults::new(1);

        let rc = blast_redo_alignment_core_mt_in_memory_subject(
            crate::program::BLASTX,
            1,
            &query,
            &qi,
            &mut kbp,
            &mut mtx,
            &mut scoring,
            &params,
            &mut saved,
            &mut hsp_list,
            &mut results,
            BlastRedoInMemorySubject {
                subject_source: &subject,
                reward: 0,
                penalty: 0,
                genetic_code: &crate::util::STANDARD_GENETIC_CODE,
                smith_waterman: false,
                expect_value: 10.0,
                hitlist_size: 10,
                inclusion_ethresh: 10.0,
                link_context: None,
            },
        );

        assert_eq!(rc, 0);
        assert_eq!(hsp_list.hsps.len(), 1);
        assert_eq!(hsp_list.hsps[0].query_frame, 1);
        assert_eq!(hsp_list.hsps[0].num_ident, frame_len);
        assert_eq!(
            hsp_list.hsps[0].comp_adjustment_method,
            CompoAdjustMode::CompositionBasedStats as i32
        );
        assert!(hsp_list.hsps[0].edit_script.is_some());
        let hitlist = results.hitlists[0].as_ref().expect("query hitlist");
        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(
            hitlist.hsp_lists[0].hsps[0].comp_adjustment_method,
            CompoAdjustMode::CompositionBasedStats as i32
        );

        let mut hsp_list_sw = HspList::new(7);
        hsp_list_sw.add_hsp(Hsp {
            score: 10,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1.0,
            query_offset: 0,
            query_end: frame_len,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: frame_len,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 1,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        let mut results_sw = crate::hspstream::HspResults::new(1);

        let rc_sw = blast_redo_alignment_core_mt_in_memory_subject(
            crate::program::BLASTX,
            1,
            &query,
            &qi,
            &mut kbp,
            &mut mtx,
            &mut scoring,
            &params,
            &mut saved,
            &mut hsp_list_sw,
            &mut results_sw,
            BlastRedoInMemorySubject {
                subject_source: &subject,
                reward: 0,
                penalty: 0,
                genetic_code: &crate::util::STANDARD_GENETIC_CODE,
                smith_waterman: true,
                expect_value: 10.0,
                hitlist_size: 10,
                inclusion_ethresh: 10.0,
                link_context: None,
            },
        );

        assert_eq!(rc_sw, 0);
        assert_eq!(hsp_list_sw.hsps.len(), 1);
        assert_eq!(hsp_list_sw.hsps[0].query_frame, 1);
        assert_eq!(hsp_list_sw.hsps[0].num_ident, frame_len);
        assert_eq!(
            hsp_list_sw.hsps[0].comp_adjustment_method,
            CompoAdjustMode::CompositionBasedStats as i32
        );
        assert!(hsp_list_sw.hsps[0].edit_script.is_some());
        let hitlist_sw = results_sw.hitlists[0].as_ref().expect("query hitlist");
        assert_eq!(hitlist_sw.hsp_lists.len(), 1);
        assert_eq!(
            hitlist_sw.hsp_lists[0].hsps[0].comp_adjustment_method,
            CompoAdjustMode::CompositionBasedStats as i32
        );
    }

    #[test]
    fn blast_redo_alignment_core_mt_in_memory_subject_redoes_second_blastx_context() {
        let query_ncbi4na = vec![1u8, 8, 4, 4, 2, 8, 1, 8, 4, 4, 2, 8];
        let (query, offsets) = crate::util::blast_get_all_translations(
            &query_ncbi4na,
            query_ncbi4na.len(),
            &crate::util::STANDARD_GENETIC_CODE,
        );
        let mut qi = crate::queryinfo::QueryInfo::new_translated_query_from_offsets(&offsets);
        for ctx in &mut qi.contexts {
            ctx.eff_searchsp = 100;
        }
        let frame_len = qi.contexts[1].query_length;
        let frame_begin = qi.contexts[1].query_offset as usize;
        let subject = query[frame_begin..frame_begin + frame_len as usize].to_vec();
        let mut kbp = vec![
            crate::stat::KarlinBlk {
                lambda: 0.267,
                k: 0.041,
                log_k: 0.041_f64.ln(),
                h: 0.14,
                round_down: true,
            };
            qi.contexts.len()
        ];
        let mut mtx = vec![vec![0i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE];
        let mut scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
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
            },
            1.0,
        );
        let mut matrix_info = BlastMatrixInfo::default();
        assert_eq!(
            s_matrix_info_init(&mut matrix_info, "BLOSUM62", 0.3176, 1.0),
            0
        );
        let params = blast_redo_align_params_new(
            matrix_info,
            BlastCompoGappingParams {
                gap_open: 11,
                gap_extend: 1,
                decline_align: i32::MIN,
                x_dropoff: 30,
                context: None,
            },
            CompoAdjustMode::NoCompositionBasedStats,
            1.0,
            false,
            true,
            false,
            query_ncbi4na.len() as i32,
            1,
            10.0,
            false,
            0.0,
        );
        let mut saved = BlastKappaSavedParameters::s_saved_parameters_new(
            0,
            1,
            CompoAdjustMode::NoCompositionBasedStats,
            false,
        );
        let mut hsp_list = HspList::new(7);
        hsp_list.add_hsp(Hsp {
            score: 10,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1.0,
            query_offset: 0,
            query_end: frame_len,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: frame_len,
            subject_gapped_start: 0,
            context: 1,
            query_frame: 2,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        let mut results = crate::hspstream::HspResults::new(1);

        let rc = blast_redo_alignment_core_mt_in_memory_subject(
            crate::program::BLASTX,
            1,
            &query,
            &qi,
            &mut kbp,
            &mut mtx,
            &mut scoring,
            &params,
            &mut saved,
            &mut hsp_list,
            &mut results,
            BlastRedoInMemorySubject {
                subject_source: &subject,
                reward: 0,
                penalty: 0,
                genetic_code: &crate::util::STANDARD_GENETIC_CODE,
                smith_waterman: false,
                expect_value: 10.0,
                hitlist_size: 10,
                inclusion_ethresh: 10.0,
                link_context: None,
            },
        );

        assert_eq!(rc, 0);
        assert_eq!(hsp_list.hsps.len(), 1);
        assert_eq!(hsp_list.hsps[0].context, 1);
        assert_eq!(hsp_list.hsps[0].query_frame, 2);
        assert_eq!(hsp_list.hsps[0].subject_frame, 0);
        assert_eq!(hsp_list.hsps[0].num_ident, frame_len);
        let hitlist = results.hitlists[0].as_ref().expect("query hitlist");
        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(hitlist.hsp_lists[0].hsps[0].context, 1);
        assert_eq!(hitlist.hsp_lists[0].hsps[0].num_ident, frame_len);
    }

    #[test]
    fn blast_redo_alignment_core_mt_in_memory_subject_redoes_tblastx_match() {
        let query_ncbi4na = vec![1u8, 8, 4, 4, 2, 8, 1, 8, 4, 4, 2, 8];
        let (query, offsets) = crate::util::blast_get_all_translations(
            &query_ncbi4na,
            query_ncbi4na.len(),
            &crate::util::STANDARD_GENETIC_CODE,
        );
        let mut qi = crate::queryinfo::QueryInfo::new_translated_query_from_offsets(&offsets);
        for ctx in &mut qi.contexts {
            ctx.eff_searchsp = 100;
        }
        let subject = query_ncbi4na.clone();
        let mut kbp = vec![
            crate::stat::KarlinBlk {
                lambda: 0.267,
                k: 0.041,
                log_k: 0.041_f64.ln(),
                h: 0.14,
                round_down: true,
            };
            qi.contexts.len()
        ];
        let mut mtx = vec![vec![0i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE];
        let mut scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
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
            },
            1.0,
        );
        let mut matrix_info = BlastMatrixInfo::default();
        assert_eq!(
            s_matrix_info_init(&mut matrix_info, "BLOSUM62", 0.3176, 1.0),
            0
        );
        let params = blast_redo_align_params_new(
            matrix_info,
            BlastCompoGappingParams {
                gap_open: 11,
                gap_extend: 1,
                decline_align: i32::MIN,
                x_dropoff: 30,
                context: None,
            },
            CompoAdjustMode::NoCompositionBasedStats,
            1.0,
            false,
            true,
            true,
            query_ncbi4na.len() as i32,
            1,
            10.0,
            false,
            0.0,
        );
        let mut saved = BlastKappaSavedParameters::s_saved_parameters_new(
            0,
            1,
            CompoAdjustMode::NoCompositionBasedStats,
            false,
        );
        let mut hsp_list = HspList::new(7);
        hsp_list.add_hsp(Hsp {
            score: 10,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1.0,
            query_offset: 0,
            query_end: 4,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: 4,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 1,
            subject_frame: 1,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        let mut results = crate::hspstream::HspResults::new(1);

        let rc = blast_redo_alignment_core_mt_in_memory_subject(
            crate::program::TBLASTX,
            1,
            &query,
            &qi,
            &mut kbp,
            &mut mtx,
            &mut scoring,
            &params,
            &mut saved,
            &mut hsp_list,
            &mut results,
            BlastRedoInMemorySubject {
                subject_source: &subject,
                reward: 0,
                penalty: 0,
                genetic_code: &crate::util::STANDARD_GENETIC_CODE,
                smith_waterman: false,
                expect_value: 10.0,
                hitlist_size: 10,
                inclusion_ethresh: 10.0,
                link_context: None,
            },
        );

        assert_eq!(rc, 0);
        assert_eq!(hsp_list.hsps.len(), 1);
        assert_eq!(hsp_list.hsps[0].query_frame, 1);
        assert_eq!(hsp_list.hsps[0].subject_frame, 1);
        assert_eq!(hsp_list.hsps[0].num_ident, 4);
        let hitlist = results.hitlists[0].as_ref().expect("query hitlist");
        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(hitlist.hsp_lists[0].hsps[0].num_ident, 4);

        let mut hsp_list = HspList::new(8);
        hsp_list.add_hsp(Hsp {
            score: 10,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1.0,
            query_offset: 0,
            query_end: 4,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: 4,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 1,
            subject_frame: 1,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        let mut results = crate::hspstream::HspResults::new(1);
        let rc = blast_redo_alignment_core_mt_in_memory_subject(
            crate::program::TBLASTX,
            1,
            &query,
            &qi,
            &mut kbp,
            &mut mtx,
            &mut scoring,
            &params,
            &mut saved,
            &mut hsp_list,
            &mut results,
            BlastRedoInMemorySubject {
                subject_source: &subject,
                reward: 0,
                penalty: 0,
                genetic_code: &crate::util::STANDARD_GENETIC_CODE,
                smith_waterman: true,
                expect_value: 10.0,
                hitlist_size: 10,
                inclusion_ethresh: 10.0,
                link_context: None,
            },
        );

        assert_eq!(rc, 0);
        assert_eq!(hsp_list.oid, 8);
        assert_eq!(hsp_list.hsps.len(), 1);
        assert_eq!(hsp_list.hsps[0].query_frame, 1);
        assert_eq!(hsp_list.hsps[0].subject_frame, 1);
        assert_eq!(hsp_list.hsps[0].num_ident, 4);
        assert!(hsp_list.hsps[0].edit_script.is_some());

        let mut matrix_info_comp = BlastMatrixInfo::default();
        assert_eq!(
            s_matrix_info_init(&mut matrix_info_comp, "BLOSUM62", 0.3176, 1.0),
            0
        );
        let mut mtx_comp = matrix_info_comp.matrix.clone();
        let mut scoring_comp = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
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
            },
            1.0,
        );
        let params_comp = blast_redo_align_params_new(
            matrix_info_comp,
            BlastCompoGappingParams {
                gap_open: 11,
                gap_extend: 1,
                decline_align: i32::MIN,
                x_dropoff: 30,
                context: None,
            },
            CompoAdjustMode::CompositionBasedStats,
            1.0,
            false,
            true,
            true,
            query_ncbi4na.len() as i32,
            1,
            10.0,
            false,
            0.0,
        );
        let mut kbp_comp = vec![
            crate::stat::KarlinBlk {
                lambda: 0.267,
                k: 0.041,
                log_k: 0.041_f64.ln(),
                h: 0.14,
                round_down: true,
            };
            qi.contexts.len()
        ];
        let mut saved_comp = BlastKappaSavedParameters::s_saved_parameters_new(
            crate::matrix::AA_SIZE as i32,
            1,
            CompoAdjustMode::CompositionBasedStats,
            false,
        );
        let mut hsp_list_comp = HspList::new(9);
        hsp_list_comp.add_hsp(Hsp {
            score: 10,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1.0,
            query_offset: 0,
            query_end: 4,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: 4,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 1,
            subject_frame: 1,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        let mut results_comp = crate::hspstream::HspResults::new(1);

        let rc_comp = blast_redo_alignment_core_mt_in_memory_subject(
            crate::program::TBLASTX,
            1,
            &query,
            &qi,
            &mut kbp_comp,
            &mut mtx_comp,
            &mut scoring_comp,
            &params_comp,
            &mut saved_comp,
            &mut hsp_list_comp,
            &mut results_comp,
            BlastRedoInMemorySubject {
                subject_source: &subject,
                reward: 0,
                penalty: 0,
                genetic_code: &crate::util::STANDARD_GENETIC_CODE,
                smith_waterman: false,
                expect_value: 10.0,
                hitlist_size: 10,
                inclusion_ethresh: 10.0,
                link_context: None,
            },
        );

        assert_eq!(rc_comp, -1);
        assert!(results_comp.hitlists[0].is_none());

        let mut hsp_list_comp_sw = HspList::new(10);
        hsp_list_comp_sw.add_hsp(Hsp {
            score: 10,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1.0,
            query_offset: 0,
            query_end: 4,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: 4,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 1,
            subject_frame: 1,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        let mut results_comp_sw = crate::hspstream::HspResults::new(1);

        let rc_comp_sw = blast_redo_alignment_core_mt_in_memory_subject(
            crate::program::TBLASTX,
            1,
            &query,
            &qi,
            &mut kbp_comp,
            &mut mtx_comp,
            &mut scoring_comp,
            &params_comp,
            &mut saved_comp,
            &mut hsp_list_comp_sw,
            &mut results_comp_sw,
            BlastRedoInMemorySubject {
                subject_source: &subject,
                reward: 0,
                penalty: 0,
                genetic_code: &crate::util::STANDARD_GENETIC_CODE,
                smith_waterman: true,
                expect_value: 10.0,
                hitlist_size: 10,
                inclusion_ethresh: 10.0,
                link_context: None,
            },
        );

        assert_eq!(rc_comp_sw, -1);
        assert!(results_comp_sw.hitlists[0].is_none());
    }

    #[test]
    fn blast_redo_alignment_core_mt_in_memory_subject_redoes_second_tblastx_context() {
        let query_ncbi4na = vec![1u8, 8, 4, 4, 2, 8, 1, 8, 4, 4, 2, 8];
        let (query, offsets) = crate::util::blast_get_all_translations(
            &query_ncbi4na,
            query_ncbi4na.len(),
            &crate::util::STANDARD_GENETIC_CODE,
        );
        let mut qi = crate::queryinfo::QueryInfo::new_translated_query_from_offsets(&offsets);
        for ctx in &mut qi.contexts {
            ctx.eff_searchsp = 100;
        }
        let frame_len = qi.contexts[1].query_length;
        let subject = query_ncbi4na.clone();
        let mut kbp = vec![
            crate::stat::KarlinBlk {
                lambda: 0.267,
                k: 0.041,
                log_k: 0.041_f64.ln(),
                h: 0.14,
                round_down: true,
            };
            qi.contexts.len()
        ];
        let mut mtx = vec![vec![0i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE];
        let mut scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
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
            },
            1.0,
        );
        let mut matrix_info = BlastMatrixInfo::default();
        assert_eq!(
            s_matrix_info_init(&mut matrix_info, "BLOSUM62", 0.3176, 1.0),
            0
        );
        let params = blast_redo_align_params_new(
            matrix_info,
            BlastCompoGappingParams {
                gap_open: 11,
                gap_extend: 1,
                decline_align: i32::MIN,
                x_dropoff: 30,
                context: None,
            },
            CompoAdjustMode::NoCompositionBasedStats,
            1.0,
            false,
            true,
            true,
            query_ncbi4na.len() as i32,
            1,
            10.0,
            false,
            0.0,
        );
        let mut saved = BlastKappaSavedParameters::s_saved_parameters_new(
            0,
            1,
            CompoAdjustMode::NoCompositionBasedStats,
            false,
        );
        let mut hsp_list = HspList::new(7);
        hsp_list.add_hsp(Hsp {
            score: 10,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1.0,
            query_offset: 0,
            query_end: frame_len,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: frame_len,
            subject_gapped_start: 0,
            context: 1,
            query_frame: 2,
            subject_frame: 2,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        let mut results = crate::hspstream::HspResults::new(1);

        let rc = blast_redo_alignment_core_mt_in_memory_subject(
            crate::program::TBLASTX,
            1,
            &query,
            &qi,
            &mut kbp,
            &mut mtx,
            &mut scoring,
            &params,
            &mut saved,
            &mut hsp_list,
            &mut results,
            BlastRedoInMemorySubject {
                subject_source: &subject,
                reward: 0,
                penalty: 0,
                genetic_code: &crate::util::STANDARD_GENETIC_CODE,
                smith_waterman: false,
                expect_value: 10.0,
                hitlist_size: 10,
                inclusion_ethresh: 10.0,
                link_context: None,
            },
        );

        assert_eq!(rc, 0);
        assert_eq!(hsp_list.hsps.len(), 1);
        assert_eq!(hsp_list.hsps[0].context, 1);
        assert_eq!(hsp_list.hsps[0].query_frame, 2);
        assert_eq!(hsp_list.hsps[0].subject_frame, 2);
        assert_eq!(hsp_list.hsps[0].num_ident, frame_len);
        let hitlist = results.hitlists[0].as_ref().expect("query hitlist");
        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(hitlist.hsp_lists[0].hsps[0].context, 1);
        assert_eq!(hitlist.hsp_lists[0].hsps[0].num_ident, frame_len);
    }

    #[test]
    fn blast_redo_alignment_core_mt_in_memory_subjects_keeps_heap_best() {
        let query = vec![0u8, 1, 2, 3, 0, 1, 2, 3, 0, 1];
        let subject1 = query.clone();
        let subject2 = query.clone();
        let qi = crate::queryinfo::QueryInfo {
            num_queries: 1,
            contexts: vec![crate::queryinfo::ContextInfo {
                query_offset: 0,
                query_length: query.len() as i32,
                eff_searchsp: 100,
                length_adjustment: 0,
                query_index: 0,
                frame: 0,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: query.len() as u32,
            min_length: 0,
        };
        let original_kbp = crate::stat::KarlinBlk {
            lambda: 0.267,
            k: 0.041,
            log_k: 0.041_f64.ln(),
            h: 0.14,
            round_down: true,
        };
        let mut kbp = vec![original_kbp.clone()];
        let mut mtx = vec![vec![0i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE];
        let mut scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
                matrix_path: None,
                reward: 1,
                penalty: -3,
                gap_open: 5,
                gap_extend: 2,
                shift_pen: i16::MAX as i32,
                gapped_calculation: true,
                complexity_adjusted_scoring: false,
                matrix_name: None,
                is_ooframe: false,
                program_number: crate::program::UNDEFINED,
            },
            1.0,
        );
        let params = blast_redo_align_params_new(
            BlastMatrixInfo::default(),
            BlastCompoGappingParams {
                gap_open: 5,
                gap_extend: 2,
                decline_align: i32::MIN,
                x_dropoff: 30,
                context: None,
            },
            CompoAdjustMode::NoCompositionBasedStats,
            1.0,
            false,
            false,
            false,
            query.len() as i32,
            1,
            10.0,
            false,
            0.0,
        );
        let make_hsp_list = |oid| {
            let mut list = HspList::new(oid);
            list.add_hsp(Hsp {
                score: 10,
                num_ident: 0,
                bit_score: 0.0,
                evalue: 1.0,
                query_offset: 2,
                query_end: 6,
                query_gapped_start: 4,
                subject_offset: 2,
                subject_end: 6,
                subject_gapped_start: 4,
                context: 0,
                query_frame: 0,
                subject_frame: 0,
                num_gaps: 0,
                comp_adjustment_method: 0,
                edit_script: None,
                pat_info: None,
                map_info: None,
            });
            list
        };
        let mut matches = vec![
            BlastRedoInMemorySubjectMatch {
                hsp_list: make_hsp_list(1),
                subject: BlastRedoInMemorySubject {
                    subject_source: &subject1,
                    reward: 1,
                    penalty: -3,
                    genetic_code: &crate::util::STANDARD_GENETIC_CODE,
                    smith_waterman: false,
                    expect_value: 10.0,
                    hitlist_size: 1,
                    inclusion_ethresh: 0.0,
                    link_context: None,
                },
            },
            BlastRedoInMemorySubjectMatch {
                hsp_list: make_hsp_list(2),
                subject: BlastRedoInMemorySubject {
                    subject_source: &subject2,
                    reward: 1,
                    penalty: -3,
                    genetic_code: &crate::util::STANDARD_GENETIC_CODE,
                    smith_waterman: false,
                    expect_value: 10.0,
                    hitlist_size: 1,
                    inclusion_ethresh: 0.0,
                    link_context: None,
                },
            },
        ];
        let mut saved = BlastKappaSavedParameters::s_saved_parameters_new(
            0,
            1,
            CompoAdjustMode::NoCompositionBasedStats,
            false,
        );
        let mut results = crate::hspstream::HspResults::new(1);

        let rc = blast_redo_alignment_core_mt_in_memory_subjects(
            crate::program::BLASTN,
            2,
            &query,
            &qi,
            &mut kbp,
            &mut mtx,
            &mut scoring,
            &params,
            &mut saved,
            &mut matches,
            &mut results,
        );

        assert_eq!(rc, 0);
        assert_eq!(kbp[0].lambda, original_kbp.lambda);
        let hitlist = results.hitlists[0].as_ref().expect("query hitlist");
        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(hitlist.hsp_lists[0].oid, 2);
        assert_eq!(hitlist.hsp_lists[0].hsps[0].num_ident, query.len() as i32);
    }

    #[test]
    fn blast_redo_alignment_core_mt_in_memory_subjects_early_terminates_late_match() {
        let query = vec![0u8, 1, 2, 3, 0, 1, 2, 3, 0, 1];
        let subject1 = query.clone();
        let subject2 = query.clone();
        let qi = crate::queryinfo::QueryInfo {
            num_queries: 1,
            contexts: vec![crate::queryinfo::ContextInfo {
                query_offset: 0,
                query_length: query.len() as i32,
                eff_searchsp: 100,
                length_adjustment: 0,
                query_index: 0,
                frame: 0,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: query.len() as u32,
            min_length: 0,
        };
        let original_kbp = crate::stat::KarlinBlk {
            lambda: 0.267,
            k: 0.041,
            log_k: 0.041_f64.ln(),
            h: 0.14,
            round_down: true,
        };
        let mut kbp = vec![original_kbp.clone()];
        let mut mtx = vec![vec![0i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE];
        let mut scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
                matrix_path: None,
                reward: 1,
                penalty: -3,
                gap_open: 5,
                gap_extend: 2,
                shift_pen: i16::MAX as i32,
                gapped_calculation: true,
                complexity_adjusted_scoring: false,
                matrix_name: None,
                is_ooframe: false,
                program_number: crate::program::UNDEFINED,
            },
            1.0,
        );
        let params = blast_redo_align_params_new(
            BlastMatrixInfo::default(),
            BlastCompoGappingParams {
                gap_open: 5,
                gap_extend: 2,
                decline_align: i32::MIN,
                x_dropoff: 30,
                context: None,
            },
            CompoAdjustMode::NoCompositionBasedStats,
            1.0,
            false,
            false,
            false,
            query.len() as i32,
            1,
            10.0,
            false,
            0.0,
        );
        let make_hsp_list = |oid, evalue| {
            let mut list = HspList::new(oid);
            list.add_hsp(Hsp {
                score: 10,
                num_ident: 0,
                bit_score: 0.0,
                evalue,
                query_offset: 2,
                query_end: 6,
                query_gapped_start: 4,
                subject_offset: 2,
                subject_end: 6,
                subject_gapped_start: 4,
                context: 0,
                query_frame: 0,
                subject_frame: 0,
                num_gaps: 0,
                comp_adjustment_method: 0,
                edit_script: None,
                pat_info: None,
                map_info: None,
            });
            list
        };
        let mut matches = vec![
            BlastRedoInMemorySubjectMatch {
                hsp_list: make_hsp_list(1, 1.0),
                subject: BlastRedoInMemorySubject {
                    subject_source: &subject1,
                    reward: 1,
                    penalty: -3,
                    genetic_code: &crate::util::STANDARD_GENETIC_CODE,
                    smith_waterman: false,
                    expect_value: 100.0,
                    hitlist_size: 1,
                    inclusion_ethresh: 100.0,
                    link_context: None,
                },
            },
            BlastRedoInMemorySubjectMatch {
                hsp_list: make_hsp_list(2, 1000.0),
                subject: BlastRedoInMemorySubject {
                    subject_source: &subject2,
                    reward: 1,
                    penalty: -3,
                    genetic_code: &crate::util::STANDARD_GENETIC_CODE,
                    smith_waterman: false,
                    expect_value: 100.0,
                    hitlist_size: 1,
                    inclusion_ethresh: 100.0,
                    link_context: None,
                },
            },
        ];
        let mut saved = BlastKappaSavedParameters::s_saved_parameters_new(
            0,
            1,
            CompoAdjustMode::NoCompositionBasedStats,
            false,
        );
        let mut results = crate::hspstream::HspResults::new(1);

        let rc = blast_redo_alignment_core_mt_in_memory_subjects(
            crate::program::BLASTN,
            1,
            &query,
            &qi,
            &mut kbp,
            &mut mtx,
            &mut scoring,
            &params,
            &mut saved,
            &mut matches,
            &mut results,
        );

        assert_eq!(rc, 0);
        assert_eq!(matches[1].hsp_list.hsps[0].query_offset, 2);
        assert_eq!(matches[1].hsp_list.hsps[0].num_ident, 0);
        let hitlist = results.hitlists[0].as_ref().expect("query hitlist");
        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(hitlist.hsp_lists[0].oid, 1);
        assert_eq!(hitlist.hsp_lists[0].hsps[0].num_ident, query.len() as i32);
    }

    #[test]
    fn blast_redo_alignment_core_mt_in_memory_subjects_redoes_smith_waterman_matches() {
        let query = vec![0u8, 1, 2, 3, 0, 1, 2, 3, 0, 1];
        let subject1 = query.clone();
        let subject2 = query.clone();
        let qi = crate::queryinfo::QueryInfo {
            num_queries: 1,
            contexts: vec![crate::queryinfo::ContextInfo {
                query_offset: 0,
                query_length: query.len() as i32,
                eff_searchsp: 100,
                length_adjustment: 0,
                query_index: 0,
                frame: 0,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: query.len() as u32,
            min_length: 0,
        };
        let original_kbp = crate::stat::KarlinBlk {
            lambda: 0.267,
            k: 0.041,
            log_k: 0.041_f64.ln(),
            h: 0.14,
            round_down: true,
        };
        let mut kbp = vec![original_kbp.clone()];
        let mut mtx = vec![vec![0i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE];
        let mut scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
                matrix_path: None,
                reward: 1,
                penalty: -3,
                gap_open: 5,
                gap_extend: 2,
                shift_pen: i16::MAX as i32,
                gapped_calculation: true,
                complexity_adjusted_scoring: false,
                matrix_name: None,
                is_ooframe: false,
                program_number: crate::program::UNDEFINED,
            },
            1.0,
        );
        let params = blast_redo_align_params_new(
            BlastMatrixInfo::default(),
            BlastCompoGappingParams {
                gap_open: 5,
                gap_extend: 2,
                decline_align: i32::MIN,
                x_dropoff: 30,
                context: None,
            },
            CompoAdjustMode::NoCompositionBasedStats,
            1.0,
            false,
            false,
            false,
            query.len() as i32,
            1,
            10.0,
            false,
            0.0,
        );
        let make_hsp_list = |oid| {
            let mut list = HspList::new(oid);
            list.add_hsp(Hsp {
                score: 10,
                num_ident: 0,
                bit_score: 0.0,
                evalue: 1.0,
                query_offset: 2,
                query_end: 6,
                query_gapped_start: 4,
                subject_offset: 2,
                subject_end: 6,
                subject_gapped_start: 4,
                context: 0,
                query_frame: 0,
                subject_frame: 0,
                num_gaps: 0,
                comp_adjustment_method: 0,
                edit_script: None,
                pat_info: None,
                map_info: None,
            });
            list
        };
        let mut matches = vec![
            BlastRedoInMemorySubjectMatch {
                hsp_list: make_hsp_list(1),
                subject: BlastRedoInMemorySubject {
                    subject_source: &subject1,
                    reward: 1,
                    penalty: -3,
                    genetic_code: &crate::util::STANDARD_GENETIC_CODE,
                    smith_waterman: true,
                    expect_value: 10.0,
                    hitlist_size: 1,
                    inclusion_ethresh: 0.0,
                    link_context: None,
                },
            },
            BlastRedoInMemorySubjectMatch {
                hsp_list: make_hsp_list(2),
                subject: BlastRedoInMemorySubject {
                    subject_source: &subject2,
                    reward: 1,
                    penalty: -3,
                    genetic_code: &crate::util::STANDARD_GENETIC_CODE,
                    smith_waterman: true,
                    expect_value: 10.0,
                    hitlist_size: 1,
                    inclusion_ethresh: 0.0,
                    link_context: None,
                },
            },
        ];
        let mut saved = BlastKappaSavedParameters::s_saved_parameters_new(
            0,
            1,
            CompoAdjustMode::NoCompositionBasedStats,
            false,
        );
        let mut results = crate::hspstream::HspResults::new(1);

        let rc = blast_redo_alignment_core_mt_in_memory_subjects(
            crate::program::BLASTN,
            1,
            &query,
            &qi,
            &mut kbp,
            &mut mtx,
            &mut scoring,
            &params,
            &mut saved,
            &mut matches,
            &mut results,
        );

        assert_eq!(rc, 0);
        assert_eq!(kbp[0].lambda, original_kbp.lambda);
        assert_eq!(matches[0].hsp_list.hsps[0].query_offset, 0);
        assert_eq!(matches[1].hsp_list.hsps[0].query_offset, 0);
        let hitlist = results.hitlists[0].as_ref().expect("query hitlist");
        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(hitlist.hsp_lists[0].oid, 2);
        assert_eq!(hitlist.hsp_lists[0].hsps[0].num_ident, query.len() as i32);
        assert!(hitlist.hsp_lists[0].hsps[0].edit_script.is_some());
    }

    #[test]
    fn blast_redo_alignment_core_mt_in_memory_subjects_redoes_protein_stream_match() {
        let query = vec![1u8; 8];
        let subject = query.clone();
        let qi = crate::queryinfo::QueryInfo {
            num_queries: 1,
            contexts: vec![crate::queryinfo::ContextInfo {
                query_offset: 0,
                query_length: query.len() as i32,
                eff_searchsp: 100,
                length_adjustment: 0,
                query_index: 0,
                frame: 0,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: query.len() as u32,
            min_length: 0,
        };
        let mut kbp = vec![crate::stat::KarlinBlk {
            lambda: 0.267,
            k: 0.041,
            log_k: 0.041_f64.ln(),
            h: 0.14,
            round_down: true,
        }];
        let mut mtx = vec![vec![0i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE];
        let mut scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
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
            },
            1.0,
        );
        let mut matrix_info = BlastMatrixInfo::default();
        assert_eq!(
            s_matrix_info_init(&mut matrix_info, "BLOSUM62", 0.3176, 1.0),
            0
        );
        let params = blast_redo_align_params_new(
            matrix_info,
            BlastCompoGappingParams {
                gap_open: 11,
                gap_extend: 1,
                decline_align: i32::MIN,
                x_dropoff: 30,
                context: None,
            },
            CompoAdjustMode::NoCompositionBasedStats,
            1.0,
            false,
            false,
            false,
            query.len() as i32,
            1,
            10.0,
            false,
            0.0,
        );
        let mut hsp_list = HspList::new(7);
        hsp_list.add_hsp(Hsp {
            score: 10,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1.0,
            query_offset: 0,
            query_end: query.len() as i32,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: subject.len() as i32,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        let mut matches = vec![BlastRedoInMemorySubjectMatch {
            hsp_list,
            subject: BlastRedoInMemorySubject {
                subject_source: &subject,
                reward: 0,
                penalty: 0,
                genetic_code: &crate::util::STANDARD_GENETIC_CODE,
                smith_waterman: false,
                expect_value: 10.0,
                hitlist_size: 10,
                inclusion_ethresh: 10.0,
                link_context: None,
            },
        }];
        let mut saved = BlastKappaSavedParameters::s_saved_parameters_new(
            0,
            1,
            CompoAdjustMode::NoCompositionBasedStats,
            false,
        );
        let mut results = crate::hspstream::HspResults::new(1);

        let rc = blast_redo_alignment_core_mt_in_memory_subjects(
            crate::program::BLASTP,
            1,
            &query,
            &qi,
            &mut kbp,
            &mut mtx,
            &mut scoring,
            &params,
            &mut saved,
            &mut matches,
            &mut results,
        );

        assert_eq!(rc, 0);
        assert_eq!(matches[0].hsp_list.hsps[0].num_ident, query.len() as i32);
        let hitlist = results.hitlists[0].as_ref().expect("query hitlist");
        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(hitlist.hsp_lists[0].oid, 7);
        assert_eq!(hitlist.hsp_lists[0].hsps[0].num_ident, query.len() as i32);
    }

    #[test]
    fn blast_redo_alignment_core_mt_in_memory_subjects_redoes_protein_composition_adjusted_stream_match(
    ) {
        let query = vec![1u8, 2, 3, 4, 1, 2, 3, 4];
        let subject = query.clone();
        let qi = crate::queryinfo::QueryInfo {
            num_queries: 1,
            contexts: vec![crate::queryinfo::ContextInfo {
                query_offset: 0,
                query_length: query.len() as i32,
                eff_searchsp: 100,
                length_adjustment: 0,
                query_index: 0,
                frame: 0,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: query.len() as u32,
            min_length: 0,
        };
        let original_kbp = crate::stat::KarlinBlk {
            lambda: 0.267,
            k: 0.041,
            log_k: 0.041_f64.ln(),
            h: 0.14,
            round_down: true,
        };
        let mut kbp = vec![original_kbp.clone()];
        let mut matrix_info = BlastMatrixInfo::default();
        assert_eq!(
            s_matrix_info_init(&mut matrix_info, "BLOSUM62", 0.3176, 1.0),
            0
        );
        let mut mtx = matrix_info.matrix.clone();
        let original_mtx = mtx.clone();
        let mut scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
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
            },
            1.0,
        );
        let params = blast_redo_align_params_new(
            matrix_info,
            BlastCompoGappingParams {
                gap_open: 11,
                gap_extend: 1,
                decline_align: i32::MIN,
                x_dropoff: 30,
                context: None,
            },
            CompoAdjustMode::CompositionBasedStats,
            1.0,
            false,
            false,
            false,
            query.len() as i32,
            1,
            10.0,
            false,
            0.0,
        );
        let mut hsp_list = HspList::new(7);
        hsp_list.add_hsp(Hsp {
            score: 10,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1.0,
            query_offset: 0,
            query_end: query.len() as i32,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: subject.len() as i32,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        let mut matches = vec![BlastRedoInMemorySubjectMatch {
            hsp_list,
            subject: BlastRedoInMemorySubject {
                subject_source: &subject,
                reward: 0,
                penalty: 0,
                genetic_code: &crate::util::STANDARD_GENETIC_CODE,
                smith_waterman: false,
                expect_value: 10.0,
                hitlist_size: 10,
                inclusion_ethresh: 10.0,
                link_context: None,
            },
        }];
        let mut saved = BlastKappaSavedParameters::s_saved_parameters_new(
            crate::matrix::AA_SIZE as i32,
            1,
            CompoAdjustMode::CompositionBasedStats,
            false,
        );
        let mut results = crate::hspstream::HspResults::new(1);

        let rc = blast_redo_alignment_core_mt_in_memory_subjects(
            crate::program::BLASTP,
            2,
            &query,
            &qi,
            &mut kbp,
            &mut mtx,
            &mut scoring,
            &params,
            &mut saved,
            &mut matches,
            &mut results,
        );

        assert_eq!(rc, 0);
        assert_eq!(kbp[0].lambda, original_kbp.lambda);
        assert_eq!(mtx, original_mtx);
        assert_eq!(matches[0].hsp_list.hsps[0].num_ident, query.len() as i32);
        assert_eq!(
            matches[0].hsp_list.hsps[0].comp_adjustment_method,
            CompoAdjustMode::CompositionBasedStats as i32
        );
        let hitlist = results.hitlists[0].as_ref().expect("query hitlist");
        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(hitlist.hsp_lists[0].oid, 7);
        assert_eq!(hitlist.hsp_lists[0].hsps[0].num_ident, query.len() as i32);
        assert_eq!(
            hitlist.hsp_lists[0].hsps[0].comp_adjustment_method,
            CompoAdjustMode::CompositionBasedStats as i32
        );

        let mut hsp_list_sw = HspList::new(7);
        hsp_list_sw.add_hsp(Hsp {
            score: 10,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1.0,
            query_offset: 0,
            query_end: query.len() as i32,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: subject.len() as i32,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        let mut matches_sw = vec![BlastRedoInMemorySubjectMatch {
            hsp_list: hsp_list_sw,
            subject: BlastRedoInMemorySubject {
                subject_source: &subject,
                reward: 0,
                penalty: 0,
                genetic_code: &crate::util::STANDARD_GENETIC_CODE,
                smith_waterman: true,
                expect_value: 10.0,
                hitlist_size: 10,
                inclusion_ethresh: 10.0,
                link_context: None,
            },
        }];
        let mut saved_sw = BlastKappaSavedParameters::s_saved_parameters_new(
            crate::matrix::AA_SIZE as i32,
            1,
            CompoAdjustMode::CompositionBasedStats,
            false,
        );
        let mut results_sw = crate::hspstream::HspResults::new(1);

        let rc_sw = blast_redo_alignment_core_mt_in_memory_subjects(
            crate::program::BLASTP,
            2,
            &query,
            &qi,
            &mut kbp,
            &mut mtx,
            &mut scoring,
            &params,
            &mut saved_sw,
            &mut matches_sw,
            &mut results_sw,
        );

        assert_eq!(rc_sw, 0);
        assert_eq!(kbp[0].lambda, original_kbp.lambda);
        assert_eq!(mtx, original_mtx);
        assert_eq!(matches_sw[0].hsp_list.hsps[0].num_ident, query.len() as i32);
        assert_eq!(
            matches_sw[0].hsp_list.hsps[0].comp_adjustment_method,
            CompoAdjustMode::CompositionBasedStats as i32
        );
        assert!(matches_sw[0].hsp_list.hsps[0].edit_script.is_some());
        let hitlist_sw = results_sw.hitlists[0].as_ref().expect("query hitlist");
        assert_eq!(hitlist_sw.hsp_lists.len(), 1);
        assert_eq!(hitlist_sw.hsp_lists[0].oid, 7);
        assert_eq!(
            hitlist_sw.hsp_lists[0].hsps[0].num_ident,
            query.len() as i32
        );
        assert_eq!(
            hitlist_sw.hsp_lists[0].hsps[0].comp_adjustment_method,
            CompoAdjustMode::CompositionBasedStats as i32
        );
    }

    #[test]
    fn blast_redo_alignment_core_mt_in_memory_subjects_redoes_position_based_composition_adjusted_stream_match(
    ) {
        let query = vec![1u8, 2, 3, 4, 1, 2, 3, 4];
        let subject = query.clone();
        let qi = crate::queryinfo::QueryInfo {
            num_queries: 1,
            contexts: vec![crate::queryinfo::ContextInfo {
                query_offset: 0,
                query_length: query.len() as i32,
                eff_searchsp: 100,
                length_adjustment: 0,
                query_index: 0,
                frame: 0,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: query.len() as u32,
            min_length: 0,
        };
        let original_kbp = crate::stat::KarlinBlk {
            lambda: 0.267,
            k: 0.041,
            log_k: 0.041_f64.ln(),
            h: 0.14,
            round_down: true,
        };
        let mut kbp = vec![original_kbp.clone()];
        let mut matrix_info = BlastMatrixInfo::default();
        let start_numerator = vec![vec![0.0; crate::matrix::AA_SIZE]; query.len()];
        assert_eq!(
            matrix_info_init_psiblast_from_start_numerator(
                &mut matrix_info,
                &query,
                "BLOSUM62",
                &start_numerator,
                0.3176,
                1.0,
            ),
            0
        );
        let mut mtx = matrix_info.matrix.clone();
        let original_mtx = mtx.clone();
        let mut scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
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
            },
            1.0,
        );
        let params = blast_redo_align_params_new(
            matrix_info,
            BlastCompoGappingParams {
                gap_open: 11,
                gap_extend: 1,
                decline_align: i32::MIN,
                x_dropoff: 30,
                context: None,
            },
            CompoAdjustMode::CompositionBasedStats,
            1.0,
            true,
            false,
            false,
            query.len() as i32,
            1,
            10.0,
            false,
            0.0,
        );
        let mut hsp_list = HspList::new(7);
        hsp_list.add_hsp(Hsp {
            score: 10,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1.0,
            query_offset: 0,
            query_end: query.len() as i32,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: subject.len() as i32,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        let mut matches = vec![BlastRedoInMemorySubjectMatch {
            hsp_list,
            subject: BlastRedoInMemorySubject {
                subject_source: &subject,
                reward: 0,
                penalty: 0,
                genetic_code: &crate::util::STANDARD_GENETIC_CODE,
                smith_waterman: false,
                expect_value: 10.0,
                hitlist_size: 10,
                inclusion_ethresh: 10.0,
                link_context: None,
            },
        }];
        let mut saved = BlastKappaSavedParameters::s_saved_parameters_new(
            query.len() as i32,
            1,
            CompoAdjustMode::CompositionBasedStats,
            true,
        );
        let mut results = crate::hspstream::HspResults::new(1);

        let rc = blast_redo_alignment_core_mt_in_memory_subjects(
            crate::program::PSI_BLAST,
            2,
            &query,
            &qi,
            &mut kbp,
            &mut mtx,
            &mut scoring,
            &params,
            &mut saved,
            &mut matches,
            &mut results,
        );

        assert_eq!(rc, 0);
        assert_eq!(kbp[0].lambda, original_kbp.lambda);
        assert_eq!(mtx, original_mtx);
        assert_eq!(matches[0].hsp_list.hsps[0].num_ident, query.len() as i32);
        assert_eq!(
            matches[0].hsp_list.hsps[0].comp_adjustment_method,
            CompoAdjustMode::CompositionBasedStats as i32
        );
        let hitlist = results.hitlists[0].as_ref().expect("query hitlist");
        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(hitlist.hsp_lists[0].oid, 7);
        assert_eq!(hitlist.hsp_lists[0].hsps[0].num_ident, query.len() as i32);
        assert_eq!(
            hitlist.hsp_lists[0].hsps[0].comp_adjustment_method,
            CompoAdjustMode::CompositionBasedStats as i32
        );

        let mut kbp_sw = vec![original_kbp.clone()];
        let mut mtx_sw = original_mtx.clone();
        let mut scoring_sw = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
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
            },
            1.0,
        );
        let mut hsp_list_sw = HspList::new(8);
        hsp_list_sw.add_hsp(Hsp {
            score: 10,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1.0,
            query_offset: 0,
            query_end: query.len() as i32,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: subject.len() as i32,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        let mut matches_sw = vec![BlastRedoInMemorySubjectMatch {
            hsp_list: hsp_list_sw,
            subject: BlastRedoInMemorySubject {
                subject_source: &subject,
                reward: 0,
                penalty: 0,
                genetic_code: &crate::util::STANDARD_GENETIC_CODE,
                smith_waterman: true,
                expect_value: 10.0,
                hitlist_size: 10,
                inclusion_ethresh: 10.0,
                link_context: None,
            },
        }];
        let mut saved_sw = BlastKappaSavedParameters::s_saved_parameters_new(
            query.len() as i32,
            1,
            CompoAdjustMode::CompositionBasedStats,
            true,
        );
        let mut results_sw = crate::hspstream::HspResults::new(1);

        let rc_sw = blast_redo_alignment_core_mt_in_memory_subjects(
            crate::program::PSI_BLAST,
            2,
            &query,
            &qi,
            &mut kbp_sw,
            &mut mtx_sw,
            &mut scoring_sw,
            &params,
            &mut saved_sw,
            &mut matches_sw,
            &mut results_sw,
        );

        assert_eq!(rc_sw, 0);
        assert_eq!(kbp_sw[0].lambda, original_kbp.lambda);
        assert_eq!(mtx_sw, original_mtx);
        assert_eq!(matches_sw[0].hsp_list.hsps[0].num_ident, query.len() as i32);
        assert_eq!(
            matches_sw[0].hsp_list.hsps[0].comp_adjustment_method,
            CompoAdjustMode::CompositionBasedStats as i32
        );
        assert!(matches_sw[0].hsp_list.hsps[0].edit_script.is_some());
        let hitlist_sw = results_sw.hitlists[0].as_ref().expect("query hitlist");
        assert_eq!(hitlist_sw.hsp_lists.len(), 1);
        assert_eq!(hitlist_sw.hsp_lists[0].oid, 8);
        assert_eq!(
            hitlist_sw.hsp_lists[0].hsps[0].num_ident,
            query.len() as i32
        );
    }

    #[test]
    fn blast_redo_alignment_core_mt_in_memory_subjects_redoes_tblastn_stream_match() {
        let query = vec![12u8, 1, 12, 1];
        let subject = vec![1u8, 8, 4, 4, 2, 8, 1, 8, 4, 4, 2, 8];
        let qi = crate::queryinfo::QueryInfo {
            num_queries: 1,
            contexts: vec![crate::queryinfo::ContextInfo {
                query_offset: 0,
                query_length: query.len() as i32,
                eff_searchsp: 100,
                length_adjustment: 0,
                query_index: 0,
                frame: 0,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: query.len() as u32,
            min_length: 0,
        };
        let mut kbp = vec![crate::stat::KarlinBlk {
            lambda: 0.267,
            k: 0.041,
            log_k: 0.041_f64.ln(),
            h: 0.14,
            round_down: true,
        }];
        let mut mtx = vec![vec![0i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE];
        let mut scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
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
            },
            1.0,
        );
        let mut matrix_info = BlastMatrixInfo::default();
        assert_eq!(
            s_matrix_info_init(&mut matrix_info, "BLOSUM62", 0.3176, 1.0),
            0
        );
        let params = blast_redo_align_params_new(
            matrix_info,
            BlastCompoGappingParams {
                gap_open: 11,
                gap_extend: 1,
                decline_align: i32::MIN,
                x_dropoff: 30,
                context: None,
            },
            CompoAdjustMode::NoCompositionBasedStats,
            1.0,
            false,
            false,
            true,
            query.len() as i32,
            1,
            10.0,
            false,
            0.0,
        );
        let mut hsp_list = HspList::new(7);
        hsp_list.add_hsp(Hsp {
            score: 10,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1.0,
            query_offset: 0,
            query_end: query.len() as i32,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: query.len() as i32,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 0,
            subject_frame: 1,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        let mut matches = vec![BlastRedoInMemorySubjectMatch {
            hsp_list,
            subject: BlastRedoInMemorySubject {
                subject_source: &subject,
                reward: 0,
                penalty: 0,
                genetic_code: &crate::util::STANDARD_GENETIC_CODE,
                smith_waterman: false,
                expect_value: 10.0,
                hitlist_size: 10,
                inclusion_ethresh: 10.0,
                link_context: None,
            },
        }];
        let mut saved = BlastKappaSavedParameters::s_saved_parameters_new(
            0,
            1,
            CompoAdjustMode::NoCompositionBasedStats,
            false,
        );
        let mut results = crate::hspstream::HspResults::new(1);

        let rc = blast_redo_alignment_core_mt_in_memory_subjects(
            crate::program::TBLASTN,
            1,
            &query,
            &qi,
            &mut kbp,
            &mut mtx,
            &mut scoring,
            &params,
            &mut saved,
            &mut matches,
            &mut results,
        );

        assert_eq!(rc, 0);
        assert_eq!(matches[0].hsp_list.hsps[0].subject_frame, 1);
        assert_eq!(matches[0].hsp_list.hsps[0].num_ident, query.len() as i32);
        let hitlist = results.hitlists[0].as_ref().expect("query hitlist");
        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(hitlist.hsp_lists[0].oid, 7);
        assert_eq!(hitlist.hsp_lists[0].hsps[0].num_ident, query.len() as i32);
    }

    #[test]
    fn blast_redo_alignment_core_mt_in_memory_subjects_use_link_hsp_context() {
        let query = vec![1u8; 8];
        let subject = query.clone();
        let qi = crate::queryinfo::QueryInfo {
            num_queries: 1,
            contexts: vec![crate::queryinfo::ContextInfo {
                query_offset: 0,
                query_length: query.len() as i32,
                eff_searchsp: 100,
                length_adjustment: 0,
                query_index: 0,
                frame: 0,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: query.len() as u32,
            min_length: 0,
        };
        let karlin = crate::stat::KarlinBlk {
            lambda: 0.267,
            k: 0.041,
            log_k: 0.041_f64.ln(),
            h: 0.14,
            round_down: true,
        };
        let mut kbp = vec![karlin.clone()];
        let score_block = crate::link_hsps::LinkScoreBlock {
            kbp: vec![karlin],
            kbp_gap: Vec::new(),
            ..Default::default()
        };
        let link_params = crate::link_hsps::LinkHSPParameters::default();
        let link_context = HitlistLinkContext {
            query_info: &qi,
            query_context: 0,
            score_block: &score_block,
            link_params: &link_params,
            gapped_calculation: false,
        };
        let mut mtx = vec![vec![0i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE];
        let mut scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
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
            },
            1.0,
        );
        let mut matrix_info = BlastMatrixInfo::default();
        assert_eq!(
            s_matrix_info_init(&mut matrix_info, "BLOSUM62", 0.3176, 1.0),
            0
        );
        let params = blast_redo_align_params_new(
            matrix_info,
            BlastCompoGappingParams {
                gap_open: 11,
                gap_extend: 1,
                decline_align: i32::MIN,
                x_dropoff: 30,
                context: None,
            },
            CompoAdjustMode::NoCompositionBasedStats,
            1.0,
            false,
            false,
            false,
            query.len() as i32,
            1,
            10.0,
            true,
            0.0,
        );
        let mut hsp_list = HspList::new(7);
        hsp_list.add_hsp(Hsp {
            score: 10,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1.0,
            query_offset: 0,
            query_end: query.len() as i32,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: subject.len() as i32,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        let mut matches = vec![BlastRedoInMemorySubjectMatch {
            hsp_list,
            subject: BlastRedoInMemorySubject {
                subject_source: &subject,
                reward: 0,
                penalty: 0,
                genetic_code: &crate::util::STANDARD_GENETIC_CODE,
                smith_waterman: false,
                expect_value: 10.0,
                hitlist_size: 10,
                inclusion_ethresh: 10.0,
                link_context: Some(&link_context),
            },
        }];
        let mut saved = BlastKappaSavedParameters::s_saved_parameters_new(
            0,
            1,
            CompoAdjustMode::NoCompositionBasedStats,
            false,
        );
        let mut results = crate::hspstream::HspResults::new(1);

        let rc = blast_redo_alignment_core_mt_in_memory_subjects(
            crate::program::BLASTP,
            1,
            &query,
            &qi,
            &mut kbp,
            &mut mtx,
            &mut scoring,
            &params,
            &mut saved,
            &mut matches,
            &mut results,
        );

        assert_eq!(rc, 0);
        assert_eq!(matches[0].hsp_list.hsps[0].num_ident, query.len() as i32);
        let hitlist = results.hitlists[0].as_ref().expect("query hitlist");
        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(hitlist.hsp_lists[0].oid, 7);
        assert_eq!(hitlist.hsp_lists[0].hsps[0].num_ident, query.len() as i32);
    }

    #[test]
    fn blast_redo_alignment_core_mt_in_memory_subjects_redoes_blastx_stream_match() {
        let query_ncbi4na = vec![1u8, 8, 4, 4, 2, 8, 1, 8, 4, 4, 2, 8];
        let (query, offsets) = crate::util::blast_get_all_translations(
            &query_ncbi4na,
            query_ncbi4na.len(),
            &crate::util::STANDARD_GENETIC_CODE,
        );
        let mut qi = crate::queryinfo::QueryInfo::new_translated_query_from_offsets(&offsets);
        for ctx in &mut qi.contexts {
            ctx.eff_searchsp = 100;
        }
        let subject = vec![12u8, 1, 12, 1];
        let mut kbp = vec![
            crate::stat::KarlinBlk {
                lambda: 0.267,
                k: 0.041,
                log_k: 0.041_f64.ln(),
                h: 0.14,
                round_down: true,
            };
            qi.contexts.len()
        ];
        let mut mtx = vec![vec![0i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE];
        let mut scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
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
            },
            1.0,
        );
        let mut matrix_info = BlastMatrixInfo::default();
        assert_eq!(
            s_matrix_info_init(&mut matrix_info, "BLOSUM62", 0.3176, 1.0),
            0
        );
        let params = blast_redo_align_params_new(
            matrix_info,
            BlastCompoGappingParams {
                gap_open: 11,
                gap_extend: 1,
                decline_align: i32::MIN,
                x_dropoff: 30,
                context: None,
            },
            CompoAdjustMode::NoCompositionBasedStats,
            1.0,
            false,
            true,
            false,
            query_ncbi4na.len() as i32,
            1,
            10.0,
            false,
            0.0,
        );
        let mut hsp_list = HspList::new(7);
        hsp_list.add_hsp(Hsp {
            score: 10,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1.0,
            query_offset: 0,
            query_end: subject.len() as i32,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: subject.len() as i32,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 1,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        let mut matches = vec![BlastRedoInMemorySubjectMatch {
            hsp_list,
            subject: BlastRedoInMemorySubject {
                subject_source: &subject,
                reward: 0,
                penalty: 0,
                genetic_code: &crate::util::STANDARD_GENETIC_CODE,
                smith_waterman: false,
                expect_value: 10.0,
                hitlist_size: 10,
                inclusion_ethresh: 10.0,
                link_context: None,
            },
        }];
        let mut saved = BlastKappaSavedParameters::s_saved_parameters_new(
            0,
            1,
            CompoAdjustMode::NoCompositionBasedStats,
            false,
        );
        let mut results = crate::hspstream::HspResults::new(1);

        let rc = blast_redo_alignment_core_mt_in_memory_subjects(
            crate::program::BLASTX,
            1,
            &query,
            &qi,
            &mut kbp,
            &mut mtx,
            &mut scoring,
            &params,
            &mut saved,
            &mut matches,
            &mut results,
        );

        assert_eq!(rc, 0);
        assert_eq!(matches[0].hsp_list.hsps[0].query_frame, 1);
        assert_eq!(matches[0].hsp_list.hsps[0].num_ident, subject.len() as i32);
        let hitlist = results.hitlists[0].as_ref().expect("query hitlist");
        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(hitlist.hsp_lists[0].oid, 7);
        assert_eq!(hitlist.hsp_lists[0].hsps[0].num_ident, subject.len() as i32);

        let mut matrix_info_comp = BlastMatrixInfo::default();
        assert_eq!(
            s_matrix_info_init(&mut matrix_info_comp, "BLOSUM62", 0.3176, 1.0),
            0
        );
        let mut mtx_comp = matrix_info_comp.matrix.clone();
        let mut scoring_comp = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
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
            },
            1.0,
        );
        let params_comp = blast_redo_align_params_new(
            matrix_info_comp,
            BlastCompoGappingParams {
                gap_open: 11,
                gap_extend: 1,
                decline_align: i32::MIN,
                x_dropoff: 30,
                context: None,
            },
            CompoAdjustMode::CompositionBasedStats,
            1.0,
            false,
            true,
            false,
            query_ncbi4na.len() as i32,
            1,
            10.0,
            false,
            0.0,
        );
        let mut hsp_list_sw = HspList::new(7);
        hsp_list_sw.add_hsp(Hsp {
            score: 10,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1.0,
            query_offset: 0,
            query_end: subject.len() as i32,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: subject.len() as i32,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 1,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        let mut matches_sw = vec![BlastRedoInMemorySubjectMatch {
            hsp_list: hsp_list_sw,
            subject: BlastRedoInMemorySubject {
                subject_source: &subject,
                reward: 0,
                penalty: 0,
                genetic_code: &crate::util::STANDARD_GENETIC_CODE,
                smith_waterman: true,
                expect_value: 10.0,
                hitlist_size: 10,
                inclusion_ethresh: 10.0,
                link_context: None,
            },
        }];
        let mut kbp_comp = vec![
            crate::stat::KarlinBlk {
                lambda: 0.267,
                k: 0.041,
                log_k: 0.041_f64.ln(),
                h: 0.14,
                round_down: true,
            };
            qi.contexts.len()
        ];
        let mut saved_comp = BlastKappaSavedParameters::s_saved_parameters_new(
            crate::matrix::AA_SIZE as i32,
            1,
            CompoAdjustMode::CompositionBasedStats,
            false,
        );
        let mut results_sw = crate::hspstream::HspResults::new(1);

        let rc_sw = blast_redo_alignment_core_mt_in_memory_subjects(
            crate::program::BLASTX,
            2,
            &query,
            &qi,
            &mut kbp_comp,
            &mut mtx_comp,
            &mut scoring_comp,
            &params_comp,
            &mut saved_comp,
            &mut matches_sw,
            &mut results_sw,
        );

        assert_eq!(rc_sw, 0);
        assert_eq!(matches_sw[0].hsp_list.hsps[0].query_frame, 1);
        assert_eq!(
            matches_sw[0].hsp_list.hsps[0].num_ident,
            subject.len() as i32
        );
        assert_eq!(
            matches_sw[0].hsp_list.hsps[0].comp_adjustment_method,
            CompoAdjustMode::CompositionBasedStats as i32
        );
        assert!(matches_sw[0].hsp_list.hsps[0].edit_script.is_some());
        let hitlist_sw = results_sw.hitlists[0].as_ref().expect("query hitlist");
        assert_eq!(hitlist_sw.hsp_lists.len(), 1);
        assert_eq!(hitlist_sw.hsp_lists[0].oid, 7);
        assert_eq!(
            hitlist_sw.hsp_lists[0].hsps[0].comp_adjustment_method,
            CompoAdjustMode::CompositionBasedStats as i32
        );
    }

    #[test]
    fn blast_redo_alignment_core_mt_in_memory_subjects_redoes_tblastx_stream_match() {
        let query_ncbi4na = vec![1u8, 8, 4, 4, 2, 8, 1, 8, 4, 4, 2, 8];
        let (query, offsets) = crate::util::blast_get_all_translations(
            &query_ncbi4na,
            query_ncbi4na.len(),
            &crate::util::STANDARD_GENETIC_CODE,
        );
        let mut qi = crate::queryinfo::QueryInfo::new_translated_query_from_offsets(&offsets);
        for ctx in &mut qi.contexts {
            ctx.eff_searchsp = 100;
        }
        let subject = query_ncbi4na.clone();
        let mut kbp = vec![
            crate::stat::KarlinBlk {
                lambda: 0.267,
                k: 0.041,
                log_k: 0.041_f64.ln(),
                h: 0.14,
                round_down: true,
            };
            qi.contexts.len()
        ];
        let mut mtx = vec![vec![0i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE];
        let mut scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
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
            },
            1.0,
        );
        let mut matrix_info = BlastMatrixInfo::default();
        assert_eq!(
            s_matrix_info_init(&mut matrix_info, "BLOSUM62", 0.3176, 1.0),
            0
        );
        let params = blast_redo_align_params_new(
            matrix_info,
            BlastCompoGappingParams {
                gap_open: 11,
                gap_extend: 1,
                decline_align: i32::MIN,
                x_dropoff: 30,
                context: None,
            },
            CompoAdjustMode::NoCompositionBasedStats,
            1.0,
            false,
            true,
            true,
            query_ncbi4na.len() as i32,
            1,
            10.0,
            false,
            0.0,
        );
        let mut hsp_list = HspList::new(7);
        hsp_list.add_hsp(Hsp {
            score: 10,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1.0,
            query_offset: 0,
            query_end: 4,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: 4,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 1,
            subject_frame: 1,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        let mut matches = vec![BlastRedoInMemorySubjectMatch {
            hsp_list,
            subject: BlastRedoInMemorySubject {
                subject_source: &subject,
                reward: 0,
                penalty: 0,
                genetic_code: &crate::util::STANDARD_GENETIC_CODE,
                smith_waterman: false,
                expect_value: 10.0,
                hitlist_size: 10,
                inclusion_ethresh: 10.0,
                link_context: None,
            },
        }];
        let mut saved = BlastKappaSavedParameters::s_saved_parameters_new(
            0,
            1,
            CompoAdjustMode::NoCompositionBasedStats,
            false,
        );
        let mut results = crate::hspstream::HspResults::new(1);

        let rc = blast_redo_alignment_core_mt_in_memory_subjects(
            crate::program::TBLASTX,
            1,
            &query,
            &qi,
            &mut kbp,
            &mut mtx,
            &mut scoring,
            &params,
            &mut saved,
            &mut matches,
            &mut results,
        );

        assert_eq!(rc, 0);
        assert_eq!(matches[0].hsp_list.hsps[0].query_frame, 1);
        assert_eq!(matches[0].hsp_list.hsps[0].subject_frame, 1);
        assert_eq!(matches[0].hsp_list.hsps[0].num_ident, 4);
        let hitlist = results.hitlists[0].as_ref().expect("query hitlist");
        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(hitlist.hsp_lists[0].oid, 7);
        assert_eq!(hitlist.hsp_lists[0].hsps[0].num_ident, 4);

        let mut matrix_info_comp = BlastMatrixInfo::default();
        assert_eq!(
            s_matrix_info_init(&mut matrix_info_comp, "BLOSUM62", 0.3176, 1.0),
            0
        );
        let mut mtx_comp = matrix_info_comp.matrix.clone();
        let mut scoring_comp = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
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
            },
            1.0,
        );
        let params_comp = blast_redo_align_params_new(
            matrix_info_comp,
            BlastCompoGappingParams {
                gap_open: 11,
                gap_extend: 1,
                decline_align: i32::MIN,
                x_dropoff: 30,
                context: None,
            },
            CompoAdjustMode::CompositionBasedStats,
            1.0,
            false,
            true,
            true,
            query_ncbi4na.len() as i32,
            1,
            10.0,
            false,
            0.0,
        );
        let mut hsp_list_comp = HspList::new(9);
        hsp_list_comp.add_hsp(Hsp {
            score: 10,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1.0,
            query_offset: 0,
            query_end: 4,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: 4,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 1,
            subject_frame: 1,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        let mut matches_comp = vec![BlastRedoInMemorySubjectMatch {
            hsp_list: hsp_list_comp,
            subject: BlastRedoInMemorySubject {
                subject_source: &subject,
                reward: 0,
                penalty: 0,
                genetic_code: &crate::util::STANDARD_GENETIC_CODE,
                smith_waterman: false,
                expect_value: 10.0,
                hitlist_size: 10,
                inclusion_ethresh: 10.0,
                link_context: None,
            },
        }];
        let mut kbp_comp = vec![
            crate::stat::KarlinBlk {
                lambda: 0.267,
                k: 0.041,
                log_k: 0.041_f64.ln(),
                h: 0.14,
                round_down: true,
            };
            qi.contexts.len()
        ];
        let mut saved_comp = BlastKappaSavedParameters::s_saved_parameters_new(
            crate::matrix::AA_SIZE as i32,
            1,
            CompoAdjustMode::CompositionBasedStats,
            false,
        );
        let mut results_comp = crate::hspstream::HspResults::new(1);

        let rc_comp = blast_redo_alignment_core_mt_in_memory_subjects(
            crate::program::TBLASTX,
            2,
            &query,
            &qi,
            &mut kbp_comp,
            &mut mtx_comp,
            &mut scoring_comp,
            &params_comp,
            &mut saved_comp,
            &mut matches_comp,
            &mut results_comp,
        );

        assert_eq!(rc_comp, -1);
        assert!(results_comp.hitlists[0].is_none());

        let mut hsp_list_comp_sw = HspList::new(10);
        hsp_list_comp_sw.add_hsp(Hsp {
            score: 10,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1.0,
            query_offset: 0,
            query_end: 4,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: 4,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 1,
            subject_frame: 1,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        let mut matches_comp_sw = vec![BlastRedoInMemorySubjectMatch {
            hsp_list: hsp_list_comp_sw,
            subject: BlastRedoInMemorySubject {
                subject_source: &subject,
                reward: 0,
                penalty: 0,
                genetic_code: &crate::util::STANDARD_GENETIC_CODE,
                smith_waterman: true,
                expect_value: 10.0,
                hitlist_size: 10,
                inclusion_ethresh: 10.0,
                link_context: None,
            },
        }];
        let mut results_comp_sw = crate::hspstream::HspResults::new(1);

        let rc_comp_sw = blast_redo_alignment_core_mt_in_memory_subjects(
            crate::program::TBLASTX,
            2,
            &query,
            &qi,
            &mut kbp_comp,
            &mut mtx_comp,
            &mut scoring_comp,
            &params_comp,
            &mut saved_comp,
            &mut matches_comp_sw,
            &mut results_comp_sw,
        );

        assert_eq!(rc_comp_sw, -1);
        assert!(results_comp_sw.hitlists[0].is_none());
    }

    #[test]
    fn blast_redo_alignment_core_mt_in_memory_subjects_routes_by_query_index() {
        let query0 = vec![3u8, 3, 2, 2, 1, 1, 0, 0];
        let query1 = vec![0u8, 1, 2, 3, 0, 1, 2, 3, 0, 1];
        let mut query = query0.clone();
        query.extend_from_slice(&query1);
        let subject = query1.clone();
        let qi = crate::queryinfo::QueryInfo {
            num_queries: 2,
            contexts: vec![
                crate::queryinfo::ContextInfo {
                    query_offset: 0,
                    query_length: query0.len() as i32,
                    eff_searchsp: 80,
                    length_adjustment: 0,
                    query_index: 0,
                    frame: 0,
                    is_valid: true,
                    segment_flags: crate::queryinfo::E_NO_SEGMENTS,
                },
                crate::queryinfo::ContextInfo {
                    query_offset: query0.len() as i32,
                    query_length: query1.len() as i32,
                    eff_searchsp: 100,
                    length_adjustment: 0,
                    query_index: 1,
                    frame: 0,
                    is_valid: true,
                    segment_flags: crate::queryinfo::E_NO_SEGMENTS,
                },
            ],
            max_length: query1.len() as u32,
            min_length: 0,
        };
        let original_kbp = crate::stat::KarlinBlk {
            lambda: 0.267,
            k: 0.041,
            log_k: 0.041_f64.ln(),
            h: 0.14,
            round_down: true,
        };
        let context1_kbp = crate::stat::KarlinBlk {
            lambda: 0.321,
            k: 0.057,
            log_k: 0.057_f64.ln(),
            h: 0.16,
            round_down: true,
        };
        let mut kbp = vec![original_kbp.clone(), context1_kbp.clone()];
        let mut mtx = vec![vec![0i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE];
        let mut scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
                matrix_path: None,
                reward: 1,
                penalty: -3,
                gap_open: 5,
                gap_extend: 2,
                shift_pen: i16::MAX as i32,
                gapped_calculation: true,
                complexity_adjusted_scoring: false,
                matrix_name: None,
                is_ooframe: false,
                program_number: crate::program::UNDEFINED,
            },
            1.0,
        );
        let params = blast_redo_align_params_new(
            BlastMatrixInfo::default(),
            BlastCompoGappingParams {
                gap_open: 5,
                gap_extend: 2,
                decline_align: i32::MIN,
                x_dropoff: 30,
                context: None,
            },
            CompoAdjustMode::NoCompositionBasedStats,
            1.0,
            false,
            false,
            false,
            query1.len() as i32,
            2,
            10.0,
            false,
            0.0,
        );
        let mut hsp_list = HspList::new(11);
        hsp_list.add_hsp(Hsp {
            score: 10,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1.0,
            query_offset: 0,
            query_end: query1.len() as i32,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: subject.len() as i32,
            subject_gapped_start: 0,
            context: 1,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        let mut matches = vec![BlastRedoInMemorySubjectMatch {
            hsp_list,
            subject: BlastRedoInMemorySubject {
                subject_source: &subject,
                reward: 1,
                penalty: -3,
                genetic_code: &crate::util::STANDARD_GENETIC_CODE,
                smith_waterman: false,
                expect_value: 10.0,
                hitlist_size: 10,
                inclusion_ethresh: 10.0,
                link_context: None,
            },
        }];
        let mut saved = BlastKappaSavedParameters::s_saved_parameters_new(
            0,
            2,
            CompoAdjustMode::NoCompositionBasedStats,
            false,
        );
        let mut results = crate::hspstream::HspResults::new(2);

        let rc = blast_redo_alignment_core_mt_in_memory_subjects(
            crate::program::BLASTN,
            1,
            &query,
            &qi,
            &mut kbp,
            &mut mtx,
            &mut scoring,
            &params,
            &mut saved,
            &mut matches,
            &mut results,
        );

        assert_eq!(rc, 0);
        assert!(results.hitlists[0]
            .as_ref()
            .expect("query 0 hitlist")
            .hsp_lists
            .is_empty());
        let hitlist = results.hitlists[1].as_ref().expect("query 1 hitlist");
        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(hitlist.hsp_lists[0].oid, 11);
        assert_eq!(hitlist.hsp_lists[0].hsps[0].query_offset, 0);
        assert_eq!(hitlist.hsp_lists[0].hsps[0].num_ident, query1.len() as i32);
    }

    #[test]
    fn blast_redo_alignment_core_mt_seqsrc_subjects_fetches_and_redoes_blastn_match() {
        struct TestSeqSrc {
            seqs: Vec<Vec<u8>>,
            encodings: std::sync::Mutex<Vec<crate::seqsrc::SeqEncoding>>,
        }
        impl crate::seqsrc::BlastSeqSource for TestSeqSrc {
            fn num_seqs(&self) -> i32 {
                self.seqs.len() as i32
            }
            fn total_length(&self) -> i64 {
                self.seqs.iter().map(|seq| seq.len() as i64).sum()
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
                    (self.total_length() / self.num_seqs() as i64) as i32
                }
            }
            fn name(&self) -> &str {
                "test"
            }
            fn is_protein(&self) -> bool {
                false
            }
            fn seq_len(&self, oid: i32) -> i32 {
                self.seqs[oid as usize].len() as i32
            }
            fn get_sequence(
                &self,
                arg: &crate::seqsrc::GetSeqArg,
            ) -> Option<crate::seqsrc::SeqData> {
                self.encodings.lock().unwrap().push(arg.encoding);
                let sequence = self.seqs.get(arg.oid as usize)?.clone();
                Some(crate::seqsrc::SeqData {
                    length: sequence.len() as i32,
                    sequence,
                })
            }
            fn iter_oids(&self) -> Box<dyn Iterator<Item = i32> + '_> {
                Box::new(0..self.num_seqs())
            }
        }

        let query = vec![0u8, 1, 2, 3, 0, 1, 2, 3, 0, 1];
        let seqsrc = TestSeqSrc {
            seqs: vec![query.clone()],
            encodings: std::sync::Mutex::new(Vec::new()),
        };
        let qi = crate::queryinfo::QueryInfo {
            num_queries: 1,
            contexts: vec![crate::queryinfo::ContextInfo {
                query_offset: 0,
                query_length: query.len() as i32,
                eff_searchsp: 100,
                length_adjustment: 0,
                query_index: 0,
                frame: 0,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: query.len() as u32,
            min_length: 0,
        };
        let mut kbp = vec![crate::stat::KarlinBlk {
            lambda: 0.267,
            k: 0.041,
            log_k: 0.041_f64.ln(),
            h: 0.14,
            round_down: true,
        }];
        let mut mtx = vec![vec![0i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE];
        let mut scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
                matrix_path: None,
                reward: 1,
                penalty: -3,
                gap_open: 5,
                gap_extend: 2,
                shift_pen: i16::MAX as i32,
                gapped_calculation: true,
                complexity_adjusted_scoring: false,
                matrix_name: None,
                is_ooframe: false,
                program_number: crate::program::UNDEFINED,
            },
            1.0,
        );
        let params = blast_redo_align_params_new(
            BlastMatrixInfo::default(),
            BlastCompoGappingParams {
                gap_open: 5,
                gap_extend: 2,
                decline_align: i32::MIN,
                x_dropoff: 30,
                context: None,
            },
            CompoAdjustMode::NoCompositionBasedStats,
            1.0,
            false,
            false,
            false,
            query.len() as i32,
            1,
            10.0,
            false,
            0.0,
        );
        let mut hsp_list = HspList::new(0);
        hsp_list.add_hsp(Hsp {
            score: 10,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1.0,
            query_offset: 2,
            query_end: 6,
            query_gapped_start: 4,
            subject_offset: 2,
            subject_end: 6,
            subject_gapped_start: 4,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        let mut matches = vec![hsp_list];
        let mut saved = BlastKappaSavedParameters::s_saved_parameters_new(
            0,
            1,
            CompoAdjustMode::NoCompositionBasedStats,
            false,
        );
        let mut results = crate::hspstream::HspResults::new(1);

        let rc = blast_redo_alignment_core_mt_seqsrc_subjects(
            crate::program::BLASTN,
            1,
            &query,
            &qi,
            &mut kbp,
            &mut mtx,
            &mut scoring,
            &params,
            &mut saved,
            &seqsrc,
            &mut matches,
            BlastRedoSeqSrcSubjectConfig {
                reward: 1,
                penalty: -3,
                genetic_code: &crate::util::STANDARD_GENETIC_CODE,
                smith_waterman: false,
                expect_value: 10.0,
                hitlist_size: 10,
                inclusion_ethresh: 10.0,
                link_context: None,
            },
            &mut results,
        );

        assert_eq!(rc, 0);
        assert_eq!(
            *seqsrc.encodings.lock().unwrap(),
            vec![crate::seqsrc::SeqEncoding::Nucleotide]
        );
        assert_eq!(matches[0].hsps[0].num_ident, query.len() as i32);
        let hitlist = results.hitlists[0].as_ref().expect("query hitlist");
        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(hitlist.hsp_lists[0].oid, 0);

        seqsrc.encodings.lock().unwrap().clear();
        let params_adjusted = blast_redo_align_params_new(
            BlastMatrixInfo::default(),
            BlastCompoGappingParams {
                gap_open: 5,
                gap_extend: 2,
                decline_align: i32::MIN,
                x_dropoff: 30,
                context: None,
            },
            CompoAdjustMode::CompositionBasedStats,
            1.0,
            false,
            false,
            false,
            query.len() as i32,
            1,
            10.0,
            false,
            0.0,
        );
        let mut hsp_list_adjusted = HspList::new(0);
        hsp_list_adjusted.add_hsp(Hsp {
            score: 10,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1.0,
            query_offset: 2,
            query_end: 6,
            query_gapped_start: 4,
            subject_offset: 2,
            subject_end: 6,
            subject_gapped_start: 4,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        let mut adjusted_matches = vec![hsp_list_adjusted];
        let mut saved_adjusted = BlastKappaSavedParameters::s_saved_parameters_new(
            crate::matrix::AA_SIZE as i32,
            1,
            CompoAdjustMode::CompositionBasedStats,
            false,
        );
        let mut adjusted_results = crate::hspstream::HspResults::new(1);

        let rc_adjusted = blast_redo_alignment_core_mt_seqsrc_subjects(
            crate::program::BLASTN,
            2,
            &query,
            &qi,
            &mut kbp,
            &mut mtx,
            &mut scoring,
            &params_adjusted,
            &mut saved_adjusted,
            &seqsrc,
            &mut adjusted_matches,
            BlastRedoSeqSrcSubjectConfig {
                reward: 1,
                penalty: -3,
                genetic_code: &crate::util::STANDARD_GENETIC_CODE,
                smith_waterman: true,
                expect_value: 10.0,
                hitlist_size: 10,
                inclusion_ethresh: 10.0,
                link_context: None,
            },
            &mut adjusted_results,
        );

        assert_eq!(rc_adjusted, -1);
        assert_eq!(
            *seqsrc.encodings.lock().unwrap(),
            vec![crate::seqsrc::SeqEncoding::Nucleotide]
        );
        assert_eq!(adjusted_matches[0].hsps[0].query_offset, 2);
        assert!(adjusted_results.hitlists[0].is_none());
    }

    #[test]
    fn blast_redo_alignment_core_mt_seqsrc_subjects_fetches_and_redoes_protein_composition_adjusted_match(
    ) {
        struct TestSeqSrc {
            seqs: Vec<Vec<u8>>,
            encodings: std::sync::Mutex<Vec<crate::seqsrc::SeqEncoding>>,
        }
        impl crate::seqsrc::BlastSeqSource for TestSeqSrc {
            fn num_seqs(&self) -> i32 {
                self.seqs.len() as i32
            }
            fn total_length(&self) -> i64 {
                self.seqs.iter().map(|seq| seq.len() as i64).sum()
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
                    (self.total_length() / self.num_seqs() as i64) as i32
                }
            }
            fn name(&self) -> &str {
                "test"
            }
            fn is_protein(&self) -> bool {
                true
            }
            fn seq_len(&self, oid: i32) -> i32 {
                self.seqs[oid as usize].len() as i32
            }
            fn get_sequence(
                &self,
                arg: &crate::seqsrc::GetSeqArg,
            ) -> Option<crate::seqsrc::SeqData> {
                self.encodings.lock().unwrap().push(arg.encoding);
                let sequence = self.seqs.get(arg.oid as usize)?.clone();
                Some(crate::seqsrc::SeqData {
                    length: sequence.len() as i32,
                    sequence,
                })
            }
            fn iter_oids(&self) -> Box<dyn Iterator<Item = i32> + '_> {
                Box::new(0..self.num_seqs())
            }
        }

        let query = vec![1u8, 2, 3, 4, 1, 2, 3, 4];
        let seqsrc = TestSeqSrc {
            seqs: vec![query.clone()],
            encodings: std::sync::Mutex::new(Vec::new()),
        };
        let qi = crate::queryinfo::QueryInfo {
            num_queries: 1,
            contexts: vec![crate::queryinfo::ContextInfo {
                query_offset: 0,
                query_length: query.len() as i32,
                eff_searchsp: 100,
                length_adjustment: 0,
                query_index: 0,
                frame: 0,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: query.len() as u32,
            min_length: 0,
        };
        let original_kbp = crate::stat::KarlinBlk {
            lambda: 0.267,
            k: 0.041,
            log_k: 0.041_f64.ln(),
            h: 0.14,
            round_down: true,
        };
        let mut kbp = vec![original_kbp.clone()];
        let mut matrix_info = BlastMatrixInfo::default();
        assert_eq!(
            s_matrix_info_init(&mut matrix_info, "BLOSUM62", 0.3176, 1.0),
            0
        );
        let mut mtx = matrix_info.matrix.clone();
        let original_mtx = mtx.clone();
        let mut scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
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
            },
            1.0,
        );
        let params = blast_redo_align_params_new(
            matrix_info,
            BlastCompoGappingParams {
                gap_open: 11,
                gap_extend: 1,
                decline_align: i32::MIN,
                x_dropoff: 30,
                context: None,
            },
            CompoAdjustMode::CompositionBasedStats,
            1.0,
            false,
            false,
            false,
            query.len() as i32,
            1,
            10.0,
            false,
            0.0,
        );
        let mut hsp_list = HspList::new(0);
        hsp_list.add_hsp(Hsp {
            score: 10,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1.0,
            query_offset: 0,
            query_end: query.len() as i32,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: query.len() as i32,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        let mut matches = vec![hsp_list];
        let mut saved = BlastKappaSavedParameters::s_saved_parameters_new(
            crate::matrix::AA_SIZE as i32,
            1,
            CompoAdjustMode::CompositionBasedStats,
            false,
        );
        let mut results = crate::hspstream::HspResults::new(1);

        let rc = blast_redo_alignment_core_mt_seqsrc_subjects(
            crate::program::BLASTP,
            2,
            &query,
            &qi,
            &mut kbp,
            &mut mtx,
            &mut scoring,
            &params,
            &mut saved,
            &seqsrc,
            &mut matches,
            BlastRedoSeqSrcSubjectConfig {
                reward: 0,
                penalty: 0,
                genetic_code: &crate::util::STANDARD_GENETIC_CODE,
                smith_waterman: false,
                expect_value: 10.0,
                hitlist_size: 10,
                inclusion_ethresh: 10.0,
                link_context: None,
            },
            &mut results,
        );

        assert_eq!(rc, 0);
        assert_eq!(
            *seqsrc.encodings.lock().unwrap(),
            vec![crate::seqsrc::SeqEncoding::Protein]
        );
        assert_eq!(kbp[0].lambda, original_kbp.lambda);
        assert_eq!(mtx, original_mtx);
        assert_eq!(matches[0].hsps[0].num_ident, query.len() as i32);
        assert_eq!(
            matches[0].hsps[0].comp_adjustment_method,
            CompoAdjustMode::CompositionBasedStats as i32
        );
        let hitlist = results.hitlists[0].as_ref().expect("query hitlist");
        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(hitlist.hsp_lists[0].oid, 0);
        assert_eq!(hitlist.hsp_lists[0].hsps[0].num_ident, query.len() as i32);
        assert_eq!(
            hitlist.hsp_lists[0].hsps[0].comp_adjustment_method,
            CompoAdjustMode::CompositionBasedStats as i32
        );

        seqsrc.encodings.lock().unwrap().clear();
        let mut hsp_list_sw = HspList::new(0);
        hsp_list_sw.add_hsp(Hsp {
            score: 10,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1.0,
            query_offset: 0,
            query_end: query.len() as i32,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: query.len() as i32,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        let mut matches_sw = vec![hsp_list_sw];
        let mut saved_sw = BlastKappaSavedParameters::s_saved_parameters_new(
            crate::matrix::AA_SIZE as i32,
            1,
            CompoAdjustMode::CompositionBasedStats,
            false,
        );
        let mut results_sw = crate::hspstream::HspResults::new(1);

        let rc_sw = blast_redo_alignment_core_mt_seqsrc_subjects(
            crate::program::BLASTP,
            2,
            &query,
            &qi,
            &mut kbp,
            &mut mtx,
            &mut scoring,
            &params,
            &mut saved_sw,
            &seqsrc,
            &mut matches_sw,
            BlastRedoSeqSrcSubjectConfig {
                reward: 0,
                penalty: 0,
                genetic_code: &crate::util::STANDARD_GENETIC_CODE,
                smith_waterman: true,
                expect_value: 10.0,
                hitlist_size: 10,
                inclusion_ethresh: 10.0,
                link_context: None,
            },
            &mut results_sw,
        );

        assert_eq!(rc_sw, 0);
        assert_eq!(
            *seqsrc.encodings.lock().unwrap(),
            vec![crate::seqsrc::SeqEncoding::Protein]
        );
        assert_eq!(kbp[0].lambda, original_kbp.lambda);
        assert_eq!(mtx, original_mtx);
        assert_eq!(matches_sw[0].hsps[0].num_ident, query.len() as i32);
        assert_eq!(
            matches_sw[0].hsps[0].comp_adjustment_method,
            CompoAdjustMode::CompositionBasedStats as i32
        );
        assert!(matches_sw[0].hsps[0].edit_script.is_some());
        let hitlist_sw = results_sw.hitlists[0].as_ref().expect("query hitlist");
        assert_eq!(hitlist_sw.hsp_lists.len(), 1);
        assert_eq!(hitlist_sw.hsp_lists[0].oid, 0);
        assert_eq!(
            hitlist_sw.hsp_lists[0].hsps[0].num_ident,
            query.len() as i32
        );
        assert_eq!(
            hitlist_sw.hsp_lists[0].hsps[0].comp_adjustment_method,
            CompoAdjustMode::CompositionBasedStats as i32
        );
    }

    #[test]
    fn blast_redo_alignment_core_mt_seqsrc_subjects_fetches_and_redoes_position_based_composition_adjusted_match(
    ) {
        struct TestSeqSrc {
            seqs: Vec<Vec<u8>>,
            encodings: std::sync::Mutex<Vec<crate::seqsrc::SeqEncoding>>,
        }
        impl crate::seqsrc::BlastSeqSource for TestSeqSrc {
            fn num_seqs(&self) -> i32 {
                self.seqs.len() as i32
            }
            fn total_length(&self) -> i64 {
                self.seqs.iter().map(|seq| seq.len() as i64).sum()
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
                    (self.total_length() / self.num_seqs() as i64) as i32
                }
            }
            fn name(&self) -> &str {
                "test"
            }
            fn is_protein(&self) -> bool {
                true
            }
            fn seq_len(&self, oid: i32) -> i32 {
                self.seqs[oid as usize].len() as i32
            }
            fn get_sequence(
                &self,
                arg: &crate::seqsrc::GetSeqArg,
            ) -> Option<crate::seqsrc::SeqData> {
                self.encodings.lock().unwrap().push(arg.encoding);
                let sequence = self.seqs.get(arg.oid as usize)?.clone();
                Some(crate::seqsrc::SeqData {
                    length: sequence.len() as i32,
                    sequence,
                })
            }
            fn iter_oids(&self) -> Box<dyn Iterator<Item = i32> + '_> {
                Box::new(0..self.num_seqs())
            }
        }

        let query = vec![1u8, 2, 3, 4, 1, 2, 3, 4];
        let seqsrc = TestSeqSrc {
            seqs: vec![query.clone()],
            encodings: std::sync::Mutex::new(Vec::new()),
        };
        let qi = crate::queryinfo::QueryInfo {
            num_queries: 1,
            contexts: vec![crate::queryinfo::ContextInfo {
                query_offset: 0,
                query_length: query.len() as i32,
                eff_searchsp: 100,
                length_adjustment: 0,
                query_index: 0,
                frame: 0,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: query.len() as u32,
            min_length: 0,
        };
        let original_kbp = crate::stat::KarlinBlk {
            lambda: 0.267,
            k: 0.041,
            log_k: 0.041_f64.ln(),
            h: 0.14,
            round_down: true,
        };
        let mut kbp = vec![original_kbp.clone()];
        let mut matrix_info = BlastMatrixInfo::default();
        let start_numerator = vec![vec![0.0; crate::matrix::AA_SIZE]; query.len()];
        assert_eq!(
            matrix_info_init_psiblast_from_start_numerator(
                &mut matrix_info,
                &query,
                "BLOSUM62",
                &start_numerator,
                0.3176,
                1.0,
            ),
            0
        );
        let mut mtx = matrix_info.matrix.clone();
        let original_mtx = mtx.clone();
        let mut scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
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
            },
            1.0,
        );
        let params = blast_redo_align_params_new(
            matrix_info,
            BlastCompoGappingParams {
                gap_open: 11,
                gap_extend: 1,
                decline_align: i32::MIN,
                x_dropoff: 30,
                context: None,
            },
            CompoAdjustMode::CompositionBasedStats,
            1.0,
            true,
            false,
            false,
            query.len() as i32,
            1,
            10.0,
            false,
            0.0,
        );
        let mut hsp_list = HspList::new(0);
        hsp_list.add_hsp(Hsp {
            score: 10,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1.0,
            query_offset: 0,
            query_end: query.len() as i32,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: query.len() as i32,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        let mut matches = vec![hsp_list];
        let mut saved = BlastKappaSavedParameters::s_saved_parameters_new(
            query.len() as i32,
            1,
            CompoAdjustMode::CompositionBasedStats,
            true,
        );
        let mut results = crate::hspstream::HspResults::new(1);

        let rc = blast_redo_alignment_core_mt_seqsrc_subjects(
            crate::program::PSI_BLAST,
            2,
            &query,
            &qi,
            &mut kbp,
            &mut mtx,
            &mut scoring,
            &params,
            &mut saved,
            &seqsrc,
            &mut matches,
            BlastRedoSeqSrcSubjectConfig {
                reward: 0,
                penalty: 0,
                genetic_code: &crate::util::STANDARD_GENETIC_CODE,
                smith_waterman: false,
                expect_value: 10.0,
                hitlist_size: 10,
                inclusion_ethresh: 10.0,
                link_context: None,
            },
            &mut results,
        );

        assert_eq!(rc, 0);
        assert_eq!(
            *seqsrc.encodings.lock().unwrap(),
            vec![crate::seqsrc::SeqEncoding::Protein]
        );
        assert_eq!(kbp[0].lambda, original_kbp.lambda);
        assert_eq!(mtx, original_mtx);
        assert_eq!(matches[0].hsps[0].num_ident, query.len() as i32);
        assert_eq!(
            matches[0].hsps[0].comp_adjustment_method,
            CompoAdjustMode::CompositionBasedStats as i32
        );
        let hitlist = results.hitlists[0].as_ref().expect("query hitlist");
        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(hitlist.hsp_lists[0].oid, 0);
        assert_eq!(hitlist.hsp_lists[0].hsps[0].num_ident, query.len() as i32);
        assert_eq!(
            hitlist.hsp_lists[0].hsps[0].comp_adjustment_method,
            CompoAdjustMode::CompositionBasedStats as i32
        );

        seqsrc.encodings.lock().unwrap().clear();
        let mut kbp_sw = vec![original_kbp.clone()];
        let mut mtx_sw = original_mtx.clone();
        let mut scoring_sw = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
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
            },
            1.0,
        );
        let mut hsp_list_sw = HspList::new(0);
        hsp_list_sw.add_hsp(Hsp {
            score: 10,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1.0,
            query_offset: 0,
            query_end: query.len() as i32,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: query.len() as i32,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        let mut matches_sw = vec![hsp_list_sw];
        let mut saved_sw = BlastKappaSavedParameters::s_saved_parameters_new(
            query.len() as i32,
            1,
            CompoAdjustMode::CompositionBasedStats,
            true,
        );
        let mut results_sw = crate::hspstream::HspResults::new(1);

        let rc_sw = blast_redo_alignment_core_mt_seqsrc_subjects(
            crate::program::PSI_BLAST,
            2,
            &query,
            &qi,
            &mut kbp_sw,
            &mut mtx_sw,
            &mut scoring_sw,
            &params,
            &mut saved_sw,
            &seqsrc,
            &mut matches_sw,
            BlastRedoSeqSrcSubjectConfig {
                reward: 0,
                penalty: 0,
                genetic_code: &crate::util::STANDARD_GENETIC_CODE,
                smith_waterman: true,
                expect_value: 10.0,
                hitlist_size: 10,
                inclusion_ethresh: 10.0,
                link_context: None,
            },
            &mut results_sw,
        );

        assert_eq!(rc_sw, 0);
        assert_eq!(
            *seqsrc.encodings.lock().unwrap(),
            vec![crate::seqsrc::SeqEncoding::Protein]
        );
        assert_eq!(kbp_sw[0].lambda, original_kbp.lambda);
        assert_eq!(mtx_sw, original_mtx);
        assert_eq!(matches_sw[0].hsps[0].num_ident, query.len() as i32);
        assert_eq!(
            matches_sw[0].hsps[0].comp_adjustment_method,
            CompoAdjustMode::CompositionBasedStats as i32
        );
        assert!(matches_sw[0].hsps[0].edit_script.is_some());
        let hitlist_sw = results_sw.hitlists[0].as_ref().expect("query hitlist");
        assert_eq!(hitlist_sw.hsp_lists.len(), 1);
        assert_eq!(
            hitlist_sw.hsp_lists[0].hsps[0].num_ident,
            query.len() as i32
        );
    }

    #[test]
    fn seqsrc_encoding_for_redo_program_matches_subject_space() {
        assert_eq!(
            seqsrc_encoding_for_redo_program(crate::program::BLASTN),
            crate::seqsrc::SeqEncoding::Nucleotide
        );
        assert_eq!(
            seqsrc_encoding_for_redo_program(crate::program::BLASTP),
            crate::seqsrc::SeqEncoding::Protein
        );
        assert_eq!(
            seqsrc_encoding_for_redo_program(crate::program::BLASTX),
            crate::seqsrc::SeqEncoding::Protein
        );
        assert_eq!(
            seqsrc_encoding_for_redo_program(crate::program::PSI_BLAST),
            crate::seqsrc::SeqEncoding::Protein
        );
        assert_eq!(
            seqsrc_encoding_for_redo_program(crate::program::TBLASTN),
            crate::seqsrc::SeqEncoding::Ncbi4na
        );
        assert_eq!(
            seqsrc_encoding_for_redo_program(crate::program::TBLASTX),
            crate::seqsrc::SeqEncoding::Ncbi4na
        );
    }

    #[test]
    fn blast_redo_alignment_core_mt_seqsrc_subjects_fetches_and_redoes_tblastx_match() {
        struct TestSeqSrc {
            seqs: Vec<Vec<u8>>,
            encodings: std::sync::Mutex<Vec<crate::seqsrc::SeqEncoding>>,
        }
        impl crate::seqsrc::BlastSeqSource for TestSeqSrc {
            fn num_seqs(&self) -> i32 {
                self.seqs.len() as i32
            }
            fn total_length(&self) -> i64 {
                self.seqs.iter().map(|seq| seq.len() as i64).sum()
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
                    (self.total_length() / self.num_seqs() as i64) as i32
                }
            }
            fn name(&self) -> &str {
                "test"
            }
            fn is_protein(&self) -> bool {
                false
            }
            fn seq_len(&self, oid: i32) -> i32 {
                self.seqs[oid as usize].len() as i32
            }
            fn get_sequence(
                &self,
                arg: &crate::seqsrc::GetSeqArg,
            ) -> Option<crate::seqsrc::SeqData> {
                self.encodings.lock().unwrap().push(arg.encoding);
                let sequence = self.seqs.get(arg.oid as usize)?.clone();
                Some(crate::seqsrc::SeqData {
                    length: sequence.len() as i32,
                    sequence,
                })
            }
            fn iter_oids(&self) -> Box<dyn Iterator<Item = i32> + '_> {
                Box::new(0..self.num_seqs())
            }
        }

        let query_ncbi4na = vec![1u8, 8, 4, 4, 2, 8, 1, 8, 4, 4, 2, 8];
        let (query, offsets) = crate::util::blast_get_all_translations(
            &query_ncbi4na,
            query_ncbi4na.len(),
            &crate::util::STANDARD_GENETIC_CODE,
        );
        let mut qi = crate::queryinfo::QueryInfo::new_translated_query_from_offsets(&offsets);
        for ctx in &mut qi.contexts {
            ctx.eff_searchsp = 100;
        }
        let seqsrc = TestSeqSrc {
            seqs: vec![query_ncbi4na.clone()],
            encodings: std::sync::Mutex::new(Vec::new()),
        };
        let mut kbp = vec![
            crate::stat::KarlinBlk {
                lambda: 0.267,
                k: 0.041,
                log_k: 0.041_f64.ln(),
                h: 0.14,
                round_down: true,
            };
            qi.contexts.len()
        ];
        let mut mtx = vec![vec![0i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE];
        let mut scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
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
            },
            1.0,
        );
        let mut matrix_info = BlastMatrixInfo::default();
        assert_eq!(
            s_matrix_info_init(&mut matrix_info, "BLOSUM62", 0.3176, 1.0),
            0
        );
        let params = blast_redo_align_params_new(
            matrix_info,
            BlastCompoGappingParams {
                gap_open: 11,
                gap_extend: 1,
                decline_align: i32::MIN,
                x_dropoff: 30,
                context: None,
            },
            CompoAdjustMode::NoCompositionBasedStats,
            1.0,
            false,
            true,
            true,
            query_ncbi4na.len() as i32,
            1,
            10.0,
            false,
            0.0,
        );
        let mut hsp_list = HspList::new(0);
        hsp_list.add_hsp(Hsp {
            score: 10,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1.0,
            query_offset: 0,
            query_end: 4,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: 4,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 1,
            subject_frame: 1,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        let mut matches = vec![hsp_list];
        let mut saved = BlastKappaSavedParameters::s_saved_parameters_new(
            0,
            1,
            CompoAdjustMode::NoCompositionBasedStats,
            false,
        );
        let mut results = crate::hspstream::HspResults::new(1);

        let rc = blast_redo_alignment_core_mt_seqsrc_subjects(
            crate::program::TBLASTX,
            1,
            &query,
            &qi,
            &mut kbp,
            &mut mtx,
            &mut scoring,
            &params,
            &mut saved,
            &seqsrc,
            &mut matches,
            BlastRedoSeqSrcSubjectConfig {
                reward: 0,
                penalty: 0,
                genetic_code: &crate::util::STANDARD_GENETIC_CODE,
                smith_waterman: false,
                expect_value: 10.0,
                hitlist_size: 10,
                inclusion_ethresh: 10.0,
                link_context: None,
            },
            &mut results,
        );

        assert_eq!(rc, 0);
        assert_eq!(
            *seqsrc.encodings.lock().unwrap(),
            vec![crate::seqsrc::SeqEncoding::Ncbi4na]
        );
        assert_eq!(matches[0].hsps[0].query_frame, 1);
        assert_eq!(matches[0].hsps[0].subject_frame, 1);
        assert_eq!(matches[0].hsps[0].num_ident, 4);
        let hitlist = results.hitlists[0].as_ref().expect("query hitlist");
        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(hitlist.hsp_lists[0].oid, 0);

        seqsrc.encodings.lock().unwrap().clear();
        let mut matrix_info_comp = BlastMatrixInfo::default();
        assert_eq!(
            s_matrix_info_init(&mut matrix_info_comp, "BLOSUM62", 0.3176, 1.0),
            0
        );
        let mut mtx_comp = matrix_info_comp.matrix.clone();
        let mut scoring_comp = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
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
            },
            1.0,
        );
        let params_comp = blast_redo_align_params_new(
            matrix_info_comp,
            BlastCompoGappingParams {
                gap_open: 11,
                gap_extend: 1,
                decline_align: i32::MIN,
                x_dropoff: 30,
                context: None,
            },
            CompoAdjustMode::CompositionBasedStats,
            1.0,
            false,
            true,
            true,
            query_ncbi4na.len() as i32,
            1,
            10.0,
            false,
            0.0,
        );
        let mut hsp_list_comp = HspList::new(0);
        hsp_list_comp.add_hsp(Hsp {
            score: 10,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1.0,
            query_offset: 0,
            query_end: 4,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: 4,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 1,
            subject_frame: 1,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        let mut matches_comp = vec![hsp_list_comp];
        let mut kbp_comp = vec![
            crate::stat::KarlinBlk {
                lambda: 0.267,
                k: 0.041,
                log_k: 0.041_f64.ln(),
                h: 0.14,
                round_down: true,
            };
            qi.contexts.len()
        ];
        let mut saved_comp = BlastKappaSavedParameters::s_saved_parameters_new(
            crate::matrix::AA_SIZE as i32,
            1,
            CompoAdjustMode::CompositionBasedStats,
            false,
        );
        let mut results_comp = crate::hspstream::HspResults::new(1);

        let rc_comp = blast_redo_alignment_core_mt_seqsrc_subjects(
            crate::program::TBLASTX,
            2,
            &query,
            &qi,
            &mut kbp_comp,
            &mut mtx_comp,
            &mut scoring_comp,
            &params_comp,
            &mut saved_comp,
            &seqsrc,
            &mut matches_comp,
            BlastRedoSeqSrcSubjectConfig {
                reward: 0,
                penalty: 0,
                genetic_code: &crate::util::STANDARD_GENETIC_CODE,
                smith_waterman: false,
                expect_value: 10.0,
                hitlist_size: 10,
                inclusion_ethresh: 10.0,
                link_context: None,
            },
            &mut results_comp,
        );

        assert_eq!(rc_comp, -1);
        assert_eq!(
            *seqsrc.encodings.lock().unwrap(),
            vec![crate::seqsrc::SeqEncoding::Ncbi4na]
        );
        assert!(results_comp.hitlists[0].is_none());

        seqsrc.encodings.lock().unwrap().clear();
        let mut hsp_list_comp_sw = HspList::new(0);
        hsp_list_comp_sw.add_hsp(Hsp {
            score: 10,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1.0,
            query_offset: 0,
            query_end: 4,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: 4,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 1,
            subject_frame: 1,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        let mut matches_comp_sw = vec![hsp_list_comp_sw];
        let mut saved_comp_sw = BlastKappaSavedParameters::s_saved_parameters_new(
            crate::matrix::AA_SIZE as i32,
            1,
            CompoAdjustMode::CompositionBasedStats,
            false,
        );
        let mut results_comp_sw = crate::hspstream::HspResults::new(1);

        let rc_comp_sw = blast_redo_alignment_core_mt_seqsrc_subjects(
            crate::program::TBLASTX,
            2,
            &query,
            &qi,
            &mut kbp_comp,
            &mut mtx_comp,
            &mut scoring_comp,
            &params_comp,
            &mut saved_comp_sw,
            &seqsrc,
            &mut matches_comp_sw,
            BlastRedoSeqSrcSubjectConfig {
                reward: 0,
                penalty: 0,
                genetic_code: &crate::util::STANDARD_GENETIC_CODE,
                smith_waterman: true,
                expect_value: 10.0,
                hitlist_size: 10,
                inclusion_ethresh: 10.0,
                link_context: None,
            },
            &mut results_comp_sw,
        );

        assert_eq!(rc_comp_sw, -1);
        assert_eq!(
            *seqsrc.encodings.lock().unwrap(),
            vec![crate::seqsrc::SeqEncoding::Ncbi4na]
        );
        assert!(results_comp_sw.hitlists[0].is_none());
    }

    #[test]
    fn blast_redo_alignment_core_mt_in_memory_subject_redoes_second_blastn_context() {
        let query = vec![0u8, 1, 2, 3, 0, 1, 2, 3];
        let subject = query.clone();
        let qi = crate::queryinfo::QueryInfo {
            num_queries: 1,
            contexts: vec![
                crate::queryinfo::ContextInfo {
                    query_offset: 0,
                    query_length: query.len() as i32,
                    eff_searchsp: 100,
                    length_adjustment: 0,
                    query_index: 0,
                    frame: 1,
                    is_valid: true,
                    segment_flags: crate::queryinfo::E_NO_SEGMENTS,
                },
                crate::queryinfo::ContextInfo {
                    query_offset: 0,
                    query_length: query.len() as i32,
                    eff_searchsp: 100,
                    length_adjustment: 0,
                    query_index: 0,
                    frame: -1,
                    is_valid: true,
                    segment_flags: crate::queryinfo::E_NO_SEGMENTS,
                },
            ],
            max_length: query.len() as u32,
            min_length: 0,
        };
        let original_kbp = crate::stat::KarlinBlk {
            lambda: 0.267,
            k: 0.041,
            log_k: 0.041_f64.ln(),
            h: 0.14,
            round_down: true,
        };
        let context1_kbp = crate::stat::KarlinBlk {
            lambda: 0.321,
            k: 0.057,
            log_k: 0.057_f64.ln(),
            h: 0.16,
            round_down: true,
        };
        let mut kbp = vec![original_kbp.clone(), context1_kbp.clone()];
        let mut mtx = vec![vec![0i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE];
        let mut scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
                matrix_path: None,
                reward: 1,
                penalty: -3,
                gap_open: 5,
                gap_extend: 2,
                shift_pen: i16::MAX as i32,
                gapped_calculation: true,
                complexity_adjusted_scoring: false,
                matrix_name: None,
                is_ooframe: false,
                program_number: crate::program::UNDEFINED,
            },
            1.0,
        );
        let params = blast_redo_align_params_new(
            BlastMatrixInfo::default(),
            BlastCompoGappingParams {
                gap_open: 5,
                gap_extend: 2,
                decline_align: i32::MIN,
                x_dropoff: 30,
                context: None,
            },
            CompoAdjustMode::NoCompositionBasedStats,
            1.0,
            false,
            false,
            false,
            query.len() as i32,
            1,
            10.0,
            false,
            0.0,
        );
        let mut saved = BlastKappaSavedParameters::s_saved_parameters_new(
            0,
            2,
            CompoAdjustMode::NoCompositionBasedStats,
            false,
        );
        let mut hsp_list = HspList::new(7);
        hsp_list.add_hsp(Hsp {
            score: 10,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1.0,
            query_offset: 2,
            query_end: 6,
            query_gapped_start: 4,
            subject_offset: 2,
            subject_end: 6,
            subject_gapped_start: 4,
            context: 1,
            query_frame: -1,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        let mut results = crate::hspstream::HspResults::new(1);

        let rc = blast_redo_alignment_core_mt_in_memory_subject(
            crate::program::BLASTN,
            1,
            &query,
            &qi,
            &mut kbp,
            &mut mtx,
            &mut scoring,
            &params,
            &mut saved,
            &mut hsp_list,
            &mut results,
            BlastRedoInMemorySubject {
                subject_source: &subject,
                reward: 1,
                penalty: -3,
                genetic_code: &crate::util::STANDARD_GENETIC_CODE,
                smith_waterman: false,
                expect_value: 10.0,
                hitlist_size: 10,
                inclusion_ethresh: 10.0,
                link_context: None,
            },
        );

        assert_eq!(rc, 0);
        assert_eq!(kbp[1].lambda, context1_kbp.lambda);
        assert_eq!(hsp_list.hsps.len(), 1);
        assert_eq!(hsp_list.hsps[0].context, 1);
        assert_eq!(hsp_list.hsps[0].query_frame, -1);
        assert_eq!(hsp_list.hsps[0].query_offset, 0);
        assert_eq!(hsp_list.hsps[0].num_ident, query.len() as i32);
        let expected_bits = (hsp_list.hsps[0].score as f64 * context1_kbp.lambda
            - context1_kbp.log_k)
            / crate::math::NCBIMATH_LN2;
        assert!((hsp_list.hsps[0].bit_score - expected_bits).abs() < 1e-9);
        let hitlist = results.hitlists[0].as_ref().expect("query hitlist");
        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(hitlist.hsp_lists[0].hsps[0].context, 1);
    }

    #[test]
    fn blast_redo_alignment_core_mt_in_memory_subject_rejects_link_hsps_without_context() {
        let query = vec![1u8; 8];
        let subject = query.clone();
        let qi = crate::queryinfo::QueryInfo {
            num_queries: 1,
            contexts: vec![crate::queryinfo::ContextInfo {
                query_offset: 0,
                query_length: query.len() as i32,
                eff_searchsp: 100,
                length_adjustment: 0,
                query_index: 0,
                frame: 0,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: query.len() as u32,
            min_length: 0,
        };
        let mut kbp = vec![crate::stat::KarlinBlk {
            lambda: 0.267,
            k: 0.041,
            log_k: 0.041_f64.ln(),
            h: 0.14,
            round_down: true,
        }];
        let mut mtx = vec![vec![0i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE];
        let mut scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
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
            },
            1.0,
        );
        let mut matrix_info = BlastMatrixInfo::default();
        assert_eq!(
            s_matrix_info_init(&mut matrix_info, "BLOSUM62", 0.3176, 1.0),
            0
        );
        let params = blast_redo_align_params_new(
            matrix_info,
            BlastCompoGappingParams {
                gap_open: 11,
                gap_extend: 1,
                decline_align: i32::MIN,
                x_dropoff: 30,
                context: None,
            },
            CompoAdjustMode::NoCompositionBasedStats,
            1.0,
            false,
            false,
            false,
            query.len() as i32,
            1,
            10.0,
            true,
            0.0,
        );
        let mut saved = BlastKappaSavedParameters::s_saved_parameters_new(
            0,
            1,
            CompoAdjustMode::NoCompositionBasedStats,
            false,
        );
        let mut hsp_list = HspList::new(7);
        hsp_list.add_hsp(Hsp {
            score: 10,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1.0,
            query_offset: 0,
            query_end: query.len() as i32,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: subject.len() as i32,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        let mut results = crate::hspstream::HspResults::new(1);

        let rc = blast_redo_alignment_core_mt_in_memory_subject(
            crate::program::BLASTP,
            1,
            &query,
            &qi,
            &mut kbp,
            &mut mtx,
            &mut scoring,
            &params,
            &mut saved,
            &mut hsp_list,
            &mut results,
            BlastRedoInMemorySubject {
                subject_source: &subject,
                reward: 0,
                penalty: 0,
                genetic_code: &crate::util::STANDARD_GENETIC_CODE,
                smith_waterman: false,
                expect_value: 10.0,
                hitlist_size: 10,
                inclusion_ethresh: 10.0,
                link_context: None,
            },
        );

        assert_eq!(rc, -1);
        assert!(results.hitlists[0].is_none());
    }

    #[test]
    fn blast_redo_alignment_core_mt_in_memory_subject_uses_link_hsp_context() {
        let query = vec![1u8; 8];
        let subject = query.clone();
        let qi = crate::queryinfo::QueryInfo {
            num_queries: 1,
            contexts: vec![crate::queryinfo::ContextInfo {
                query_offset: 0,
                query_length: query.len() as i32,
                eff_searchsp: 100,
                length_adjustment: 0,
                query_index: 0,
                frame: 0,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: query.len() as u32,
            min_length: 0,
        };
        let karlin = crate::stat::KarlinBlk {
            lambda: 0.267,
            k: 0.041,
            log_k: 0.041_f64.ln(),
            h: 0.14,
            round_down: true,
        };
        let mut kbp = vec![karlin.clone()];
        let score_block = crate::link_hsps::LinkScoreBlock {
            kbp: vec![karlin],
            kbp_gap: Vec::new(),
            ..Default::default()
        };
        let link_params = crate::link_hsps::LinkHSPParameters::default();
        let link_context = HitlistLinkContext {
            query_info: &qi,
            query_context: 0,
            score_block: &score_block,
            link_params: &link_params,
            gapped_calculation: false,
        };
        let mut mtx = vec![vec![0i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE];
        let mut scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
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
            },
            1.0,
        );
        let mut matrix_info = BlastMatrixInfo::default();
        assert_eq!(
            s_matrix_info_init(&mut matrix_info, "BLOSUM62", 0.3176, 1.0),
            0
        );
        let params = blast_redo_align_params_new(
            matrix_info,
            BlastCompoGappingParams {
                gap_open: 11,
                gap_extend: 1,
                decline_align: i32::MIN,
                x_dropoff: 30,
                context: None,
            },
            CompoAdjustMode::NoCompositionBasedStats,
            1.0,
            false,
            false,
            false,
            query.len() as i32,
            1,
            10.0,
            true,
            0.0,
        );
        let mut saved = BlastKappaSavedParameters::s_saved_parameters_new(
            0,
            1,
            CompoAdjustMode::NoCompositionBasedStats,
            false,
        );
        let mut hsp_list = HspList::new(7);
        hsp_list.add_hsp(Hsp {
            score: 10,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1.0,
            query_offset: 0,
            query_end: query.len() as i32,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: subject.len() as i32,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        let mut results = crate::hspstream::HspResults::new(1);

        let rc = blast_redo_alignment_core_mt_in_memory_subject(
            crate::program::BLASTP,
            1,
            &query,
            &qi,
            &mut kbp,
            &mut mtx,
            &mut scoring,
            &params,
            &mut saved,
            &mut hsp_list,
            &mut results,
            BlastRedoInMemorySubject {
                subject_source: &subject,
                reward: 0,
                penalty: 0,
                genetic_code: &crate::util::STANDARD_GENETIC_CODE,
                smith_waterman: false,
                expect_value: 10.0,
                hitlist_size: 10,
                inclusion_ethresh: 10.0,
                link_context: Some(&link_context),
            },
        );

        assert_eq!(rc, 0);
        assert_eq!(hsp_list.hsps.len(), 1);
        assert_eq!(hsp_list.hsps[0].num_ident, query.len() as i32);
        let hitlist = results.hitlists[0].as_ref().expect("query hitlist");
        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(hitlist.hsp_lists[0].oid, 7);
    }

    #[test]
    fn blast_redo_alignment_core_mt_in_memory_subject_rejects_protein_without_matrix() {
        let query = vec![1u8; 8];
        let subject = query.clone();
        let qi = crate::queryinfo::QueryInfo {
            num_queries: 1,
            contexts: vec![crate::queryinfo::ContextInfo {
                query_offset: 0,
                query_length: query.len() as i32,
                eff_searchsp: 100,
                length_adjustment: 0,
                query_index: 0,
                frame: 0,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: query.len() as u32,
            min_length: 0,
        };
        let mut kbp = vec![crate::stat::KarlinBlk {
            lambda: 0.267,
            k: 0.041,
            log_k: 0.041_f64.ln(),
            h: 0.14,
            round_down: true,
        }];
        let mut mtx = vec![vec![0i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE];
        let mut scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
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
            },
            1.0,
        );
        let params = blast_redo_align_params_new(
            BlastMatrixInfo::default(),
            BlastCompoGappingParams {
                gap_open: 11,
                gap_extend: 1,
                decline_align: i32::MIN,
                x_dropoff: 30,
                context: None,
            },
            CompoAdjustMode::NoCompositionBasedStats,
            1.0,
            false,
            false,
            false,
            query.len() as i32,
            1,
            10.0,
            false,
            0.0,
        );
        let mut saved = BlastKappaSavedParameters::s_saved_parameters_new(
            0,
            1,
            CompoAdjustMode::NoCompositionBasedStats,
            false,
        );
        let mut hsp_list = HspList::new(7);
        hsp_list.add_hsp(Hsp {
            score: 10,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1.0,
            query_offset: 0,
            query_end: query.len() as i32,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: subject.len() as i32,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        let mut results = crate::hspstream::HspResults::new(1);

        let rc = blast_redo_alignment_core_mt_in_memory_subject(
            crate::program::BLASTP,
            1,
            &query,
            &qi,
            &mut kbp,
            &mut mtx,
            &mut scoring,
            &params,
            &mut saved,
            &mut hsp_list,
            &mut results,
            BlastRedoInMemorySubject {
                subject_source: &subject,
                reward: 0,
                penalty: 0,
                genetic_code: &crate::util::STANDARD_GENETIC_CODE,
                smith_waterman: false,
                expect_value: 10.0,
                hitlist_size: 10,
                inclusion_ethresh: 10.0,
                link_context: None,
            },
        );

        assert_eq!(rc, -1);
        assert!(results.hitlists[0].is_none());
    }

    #[test]
    fn get_align_params_blastp_assembly() {
        let scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
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
        let p = s_get_align_params(
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
        // C casts the final MAX expression to Int4, truncating toward zero.
        let expected_xdrop = (25.0 * crate::math::NCBIMATH_LN2 / 0.267) as i32;
        assert_eq!(p.gapping_params.x_dropoff, expected_xdrop);
        // near_identical_cutoff = 1.74 * ln(2) / 0.267
        let expected_nic = (NEAR_IDENTICAL_BITS_PER_POSITION * crate::math::NCBIMATH_LN2) / 0.267;
        assert!((p.near_identical_cutoff - expected_nic).abs() < 1e-9);
        // Matrix info populated for BLOSUM62.
        assert_eq!(p.matrix_info.matrix_name, "BLOSUM62");
        assert_eq!(p.matrix_info.rows, crate::matrix::AA_SIZE as i32);
    }

    #[test]
    fn get_align_params_link_hsps_uses_cutoff_score_min() {
        let scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
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
        let p = s_get_align_params(
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
    fn get_align_params_preserves_position_based_matrix_flag() {
        let scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
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
        let p = s_get_align_params(
            crate::program::PSI_BLAST,
            "BLOSUM62",
            &scoring,
            &kbp,
            0.267,
            0.3176,
            1.0,
            50,
            1e-3,
            CompoAdjustMode::CompositionBasedStats,
            true,
            false,
            1000,
            1,
            25.0,
            0,
        )
        .expect("params");

        assert!(p.position_based);
        assert!(p.matrix_info.positional);
    }

    #[test]
    fn blast_redo_align_params_new_sets_pseudocounts() {
        let matrix_info = BlastMatrixInfo::default();
        let gapping = BlastCompoGappingParams {
            gap_open: 11,
            gap_extend: 1,
            decline_align: i32::MIN,
            x_dropoff: 65,
            context: None,
        };
        let p = blast_redo_align_params_new(
            matrix_info,
            gapping,
            CompoAdjustMode::CompositionMatrixAdjust,
            32.0,
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
        assert_eq!(
            p.compo_adjust_mode,
            CompoAdjustMode::CompositionMatrixAdjust
        );
        assert!(p.query_is_translated);
        assert!(!p.subject_is_translated);
        assert!(!p.position_based);
        assert_eq!(p.ccat_query_length, 1234);
        assert_eq!(p.cutoff_s, 50);
        assert_eq!(p.cutoff_e, 1e-3);
        assert_eq!(p.local_scaling_factor, 32.0);
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
            context: None,
        };
        let mut slot = Some(blast_redo_align_params_new(
            matrix_info,
            gapping,
            CompoAdjustMode::NoCompositionBasedStats,
            1.0,
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
            query_gapped_start: q + 4,
            subject_offset: s,
            subject_end: s + 10,
            subject_gapped_start: s + 5,
            context,
            query_frame: context,
            subject_frame: context + 1,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
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
        assert_eq!(head_0.query_gapped_start(), 4);
        assert_eq!(head_0.match_gapped_start(), 5);
        let next_0 = head_0.next.as_ref().expect("second in list");
        assert_eq!(next_0.score, 60);
        assert_eq!(next_0.query_gapped_start(), 64);
        assert_eq!(next_0.match_gapped_start(), 65);
        // Frame 1's list has just the single 80-score HSP.
        let head_1 = lists[1].as_ref().expect("head 1");
        assert_eq!(head_1.score, 80);
        assert_eq!(head_1.frame, 2);
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
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: 10,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 0,
            subject_frame: -2,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        }];
        let mut lists: [Option<Box<BlastCompoAlignment>>; 6] = Default::default();
        let mut counts = [0i32; 6];
        result_hsp_to_distinct_align(&mut lists, &mut counts, &hsps, 0, 32.125);
        // C: `(int)(hsp->score * localScalingFactor)` — Rust uses
        // truncating cast semantics, not `BLAST_Nint`.
        assert_eq!(lists[0].as_ref().unwrap().score, 3212);
        assert_eq!(lists[0].as_ref().unwrap().frame, -2);
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
        // Some residue should now be NCBIstdaa X.
        let masked_count = seq
            .data()
            .iter()
            .filter(|&&b| b == crate::encoding::NCBISTDAA_X)
            .count();
        assert!(masked_count > 0);
    }

    #[test]
    fn matching_sequence_release_clears_local_data_only() {
        let mut s = matching_sequence_initialize(1234, 5, 7);
        assert_eq!(s.length, 1234);
        assert_eq!(s.index, 5);
        assert_eq!(s.local_data.as_ref().map(|data| data.seq_arg.oid), Some(7));

        matching_sequence_release(&mut s);
        assert!(s.local_data.is_none());
        // C nulls local_data but does not reset length or index.
        assert_eq!(s.length, 1234);
        assert_eq!(s.index, 5);
    }

    #[test]
    fn matching_sequence_initialize_fetches_seqsrc_subject() {
        struct TestSeqSrc {
            sequence: Vec<u8>,
            releases: std::sync::Mutex<i32>,
        }

        impl crate::seqsrc::BlastSeqSource for TestSeqSrc {
            fn num_seqs(&self) -> i32 {
                1
            }
            fn total_length(&self) -> i64 {
                self.sequence.len() as i64
            }
            fn max_seq_len(&self) -> i32 {
                self.sequence.len() as i32
            }
            fn avg_seq_len(&self) -> i32 {
                self.sequence.len() as i32
            }
            fn name(&self) -> &str {
                "test"
            }
            fn is_protein(&self) -> bool {
                false
            }
            fn seq_len(&self, _oid: i32) -> i32 {
                self.sequence.len() as i32
            }
            fn get_sequence(
                &self,
                arg: &crate::seqsrc::GetSeqArg,
            ) -> Option<crate::seqsrc::SeqData> {
                assert_eq!(arg.oid, 7);
                assert_eq!(arg.encoding, crate::seqsrc::SeqEncoding::Ncbi4na);
                assert!(arg.check_oid_exclusion);
                Some(crate::seqsrc::SeqData {
                    sequence: self.sequence.clone(),
                    length: self.sequence.len() as i32,
                })
            }
            fn release_sequence(&self, arg: &mut crate::seqsrc::BlastSeqSrcGetSeqArg) {
                *self.releases.lock().unwrap() += 1;
                arg.seq = None;
            }
            fn iter_oids(&self) -> Box<dyn Iterator<Item = i32> + '_> {
                Box::new(std::iter::once(7))
            }
        }

        let seqsrc = TestSeqSrc {
            sequence: vec![1, 2, 4, 8, 1],
            releases: std::sync::Mutex::new(0),
        };
        let mut s = BlastCompoMatchingSequence::default();
        let rc =
            s_matching_sequence_initialize(&mut s, crate::program::TBLASTN, &seqsrc, 1, 7, None);
        assert_eq!(rc, 0);
        assert_eq!(s.length, 5);
        assert_eq!(s.index, 7);
        let local_data = s.local_data.as_ref().expect("local data");
        assert_eq!(local_data.prog_number, crate::program::TBLASTN);
        assert_eq!(local_data.seq_arg.oid, 7);
        assert_eq!(
            local_data.seq_arg.encoding,
            crate::seqsrc::EBlastEncoding::Ncbi4na
        );
        assert!(local_data
            .seq_arg
            .seq
            .as_ref()
            .unwrap()
            .gen_code_string
            .is_some());

        matching_sequence_release(&mut s);
        assert!(s.local_data.is_none());
        assert_eq!(s.length, 5);
        assert_eq!(*seqsrc.releases.lock().unwrap(), 1);
    }

    #[test]
    fn matching_sequence_initialize_failed_fetch_does_not_release_sequence() {
        struct TestSeqSrc {
            releases: std::sync::Mutex<i32>,
        }

        impl crate::seqsrc::BlastSeqSource for TestSeqSrc {
            fn num_seqs(&self) -> i32 {
                1
            }
            fn total_length(&self) -> i64 {
                0
            }
            fn max_seq_len(&self) -> i32 {
                0
            }
            fn avg_seq_len(&self) -> i32 {
                0
            }
            fn name(&self) -> &str {
                "test"
            }
            fn is_protein(&self) -> bool {
                false
            }
            fn seq_len(&self, _oid: i32) -> i32 {
                0
            }
            fn get_sequence(
                &self,
                arg: &crate::seqsrc::GetSeqArg,
            ) -> Option<crate::seqsrc::SeqData> {
                assert_eq!(arg.oid, 7);
                None
            }
            fn release_sequence(&self, _arg: &mut crate::seqsrc::BlastSeqSrcGetSeqArg) {
                *self.releases.lock().unwrap() += 1;
            }
            fn iter_oids(&self) -> Box<dyn Iterator<Item = i32> + '_> {
                Box::new(std::iter::once(7))
            }
        }

        let seqsrc = TestSeqSrc {
            releases: std::sync::Mutex::new(0),
        };
        let mut s = BlastCompoMatchingSequence::default();
        let rc =
            s_matching_sequence_initialize(&mut s, crate::program::TBLASTN, &seqsrc, 1, 7, None);
        assert_eq!(rc, -1);
        assert!(s.local_data.is_none());
        assert_eq!(s.length, 0);
        assert_eq!(*seqsrc.releases.lock().unwrap(), 0);
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
        // Query uses NCBIstdaa residues, including U so the prep step can
        // exercise NCBI's U-to-C substitution. Range covers idx 0..4 inclusive.
        let q = BlastCompoSequenceData {
            buffer: vec![1, crate::encoding::NCBISTDAA_U, 4, 8, 9],
            data_offset: 0,
            length: 5,
        };
        let range = BlastCompoSequenceRange {
            begin: 0,
            end: 4,
            context: 0,
        };
        let prepped = sequence_prep_query_range(&q, &range);
        // length = 4 (end - begin); buffer = length + 2 = 6.
        assert_eq!(prepped.length, 4);
        assert_eq!(prepped.buffer.len(), 6);
        // Sentinel at index 0 stays zero.
        assert_eq!(prepped.buffer[0], 0);
        // Selenocysteine replaced by cysteine.
        assert_eq!(prepped.buffer[1], 1);
        assert_eq!(prepped.buffer[2], crate::encoding::NCBISTDAA_C);
        assert_eq!(prepped.buffer[3], 4);
        assert_eq!(prepped.buffer[4], 8);
        // Trailing sentinel at length+1 = 5.
        assert_eq!(prepped.buffer[5], 0);
        // data() returns the slice starting at offset 1.
        assert_eq!(prepped.data(), &[1, 3, 4, 8]);
    }

    #[test]
    fn sequence_get_protein_range_copies_subject_with_sentinels() {
        let source = vec![9u8, 8, 7, 6, 5, 4, 3];
        let range = BlastCompoSequenceRange {
            begin: 2,
            end: 6,
            context: 0,
        };
        let seq = sequence_get_protein_range(&source, &range);

        assert_eq!(seq.length, 4);
        assert_eq!(seq.data_offset, 1);
        assert_eq!(seq.buffer, vec![0, 7, 6, 5, 4, 0]);
        assert_eq!(seq.data(), &[7, 6, 5, 4]);
    }

    #[test]
    fn sequence_get_translated_range_uses_subject_frame_and_protein_range() {
        // NCBI4na: ATG GCT translates to M A in frame +1.
        let source = vec![1u8, 8, 4, 4, 2, 8];
        let range = BlastCompoSequenceRange {
            begin: 0,
            end: 2,
            context: 1,
        };
        let seq =
            sequence_get_translated_range(&source, &range, &crate::util::STANDARD_GENETIC_CODE)
                .expect("frame +1");

        assert_eq!(seq.length, 2);
        assert_eq!(seq.data_offset, 1);
        assert_eq!(seq.buffer, vec![0, 12, 1, 0]);
        assert_eq!(seq.data(), &[12, 1]);
    }

    #[test]
    fn sequence_get_translated_range_uses_minus_frame_translation() {
        // NCBI4na: ATG GCT reverse-complements to AGC CAT, so frame -1
        // translates to S H.
        let source = vec![1u8, 8, 4, 4, 2, 8];
        let range = BlastCompoSequenceRange {
            begin: 0,
            end: 2,
            context: -1,
        };
        let seq =
            sequence_get_translated_range(&source, &range, &crate::util::STANDARD_GENETIC_CODE)
                .expect("frame -1");

        assert_eq!(seq.buffer, vec![0, 17, 8, 0]);
        assert_eq!(seq.data(), &[17, 8]);
    }

    #[test]
    fn sequence_get_translated_range_honors_explicit_genetic_code() {
        // NCBI4na: TGA is Stop(25) in code 1 and W(20) in code 4.
        let source = vec![8u8, 4, 1];
        let range = BlastCompoSequenceRange {
            begin: 0,
            end: 1,
            context: 1,
        };
        let standard =
            sequence_get_translated_range(&source, &range, &crate::util::STANDARD_GENETIC_CODE)
                .expect("standard code");
        let code4 =
            sequence_get_translated_range(&source, &range, crate::util::lookup_genetic_code(4))
                .expect("code 4");

        assert_eq!(standard.data(), &[crate::encoding::NCBISTDAA_STOP]);
        assert_eq!(code4.data(), &[20]);
    }

    #[test]
    fn sequence_get_translated_range_name_matched_wrapper_matches_helper() {
        let source = vec![1u8, 8, 4, 4, 2, 8];
        let range = BlastCompoSequenceRange {
            begin: 0,
            end: 2,
            context: 1,
        };

        let seq =
            s_sequence_get_translated_range(&source, &range, &crate::util::STANDARD_GENETIC_CODE)
                .expect("frame +1");

        assert_eq!(seq.data(), &[12, 1]);
    }

    #[test]
    fn sequence_get_range_in_memory_dispatches_non_translated_subject() {
        let query = BlastCompoSequenceData {
            buffer: vec![1, crate::encoding::NCBISTDAA_U, 4, 8, 9],
            data_offset: 0,
            length: 5,
        };
        let q_range = BlastCompoSequenceRange {
            begin: 1,
            end: 4,
            context: 0,
        };
        let subject = vec![9u8, 8, 7, 6, 5, 4];
        let s_range = BlastCompoSequenceRange {
            begin: 2,
            end: 5,
            context: 0,
        };

        let (q_seq, s_seq) = sequence_get_range_in_memory(
            crate::program::BLASTP,
            &query,
            &q_range,
            &subject,
            &s_range,
        )
        .expect("non-translated subject");

        assert_eq!(q_seq.data(), &[crate::encoding::NCBISTDAA_C, 4, 8]);
        assert_eq!(s_seq.data(), &[7, 6, 5]);
        assert_eq!(q_seq.buffer[0], 0);
        assert_eq!(s_seq.buffer[0], 0);
    }

    #[test]
    fn sequence_get_range_in_memory_dispatches_translated_subject() {
        let query = BlastCompoSequenceData {
            buffer: vec![1, crate::encoding::NCBISTDAA_U, 4],
            data_offset: 0,
            length: 3,
        };
        let q_range = BlastCompoSequenceRange {
            begin: 0,
            end: 3,
            context: 0,
        };
        let s_range = BlastCompoSequenceRange {
            begin: 0,
            end: 2,
            context: 1,
        };
        let (q_seq, s_seq) = sequence_get_range_in_memory(
            crate::program::TBLASTN,
            &query,
            &q_range,
            &[1, 8, 4, 4, 2, 8],
            &s_range,
        )
        .expect("translated subject");

        assert_eq!(q_seq.data(), &[1, 3, 4]);
        assert_eq!(s_seq.data(), &[12, 1]);
    }

    #[test]
    fn sequence_get_range_for_redo_masks_low_complexity_protein_subject() {
        let query = BlastCompoSequenceData {
            buffer: vec![1u8; 30],
            data_offset: 0,
            length: 30,
        };
        let q_range = BlastCompoSequenceRange {
            begin: 0,
            end: 30,
            context: 0,
        };
        let subject = vec![1u8; 30];
        let s_range = BlastCompoSequenceRange {
            begin: 0,
            end: 30,
            context: 0,
        };
        let align =
            BlastCompoAlignment::new(100, MatrixAdjustRule::DontAdjust, 0, 0, 30, 0, 30, 0, None);
        let mut subject_maybe_biased = true;

        let (_q_seq, s_seq) = sequence_get_range_in_memory_with_code_for_redo(
            crate::program::BLASTP,
            &query,
            &q_range,
            &subject,
            &s_range,
            &crate::util::STANDARD_GENETIC_CODE,
            &[],
            &align,
            false,
            CompoAdjustMode::CompositionMatrixAdjust,
            &mut subject_maybe_biased,
        )
        .expect("protein redo range");

        assert!(subject_maybe_biased);
        assert!(s_seq
            .data()
            .iter()
            .any(|&aa| aa == crate::encoding::NCBISTDAA_X));
    }

    #[test]
    fn sequence_get_range_for_redo_honors_tblastn_no_seg_gate() {
        let query = BlastCompoSequenceData {
            buffer: vec![1u8; 30],
            data_offset: 0,
            length: 30,
        };
        let q_range = BlastCompoSequenceRange {
            begin: 0,
            end: 30,
            context: 0,
        };
        let subject: Vec<u8> = [4u8, 2, 8].repeat(30);
        let s_range = BlastCompoSequenceRange {
            begin: 0,
            end: 30,
            context: 1,
        };
        let align =
            BlastCompoAlignment::new(100, MatrixAdjustRule::DontAdjust, 0, 0, 30, 0, 30, 1, None);
        let mut subject_maybe_biased = true;

        let (_q_seq, s_seq) = sequence_get_range_in_memory_with_code_for_redo(
            crate::program::TBLASTN,
            &query,
            &q_range,
            &subject,
            &s_range,
            &crate::util::STANDARD_GENETIC_CODE,
            &[],
            &align,
            false,
            CompoAdjustMode::CompositionMatrixAdjust,
            &mut subject_maybe_biased,
        )
        .expect("translated redo range");

        assert!(subject_maybe_biased);
        assert!(!s_seq
            .data()
            .iter()
            .any(|&aa| aa == crate::encoding::NCBISTDAA_X));
    }

    #[test]
    fn sequence_get_range_for_redo_skips_seg_for_near_identical_subject() {
        let query = BlastCompoSequenceData {
            buffer: vec![1u8; 30],
            data_offset: 0,
            length: 30,
        };
        let q_range = BlastCompoSequenceRange {
            begin: 0,
            end: 30,
            context: 0,
        };
        let subject = vec![1u8; 30];
        let s_range = BlastCompoSequenceRange {
            begin: 0,
            end: 30,
            context: 0,
        };
        let align =
            BlastCompoAlignment::new(100, MatrixAdjustRule::DontAdjust, 0, 0, 30, 0, 30, 0, None);
        let mut subject_maybe_biased = true;

        let (_q_seq, s_seq) = sequence_get_range_in_memory_with_code_for_redo(
            crate::program::BLASTP,
            &query,
            &q_range,
            &subject,
            &s_range,
            &crate::util::STANDARD_GENETIC_CODE,
            &[],
            &align,
            true,
            CompoAdjustMode::CompositionMatrixAdjust,
            &mut subject_maybe_biased,
        )
        .expect("near-identical redo range");

        assert!(subject_maybe_biased);
        assert!(!s_seq
            .data()
            .iter()
            .any(|&aa| aa == crate::encoding::NCBISTDAA_X));
    }

    #[test]
    fn sequence_get_range_for_sw_redo_masks_without_subject_biased_state() {
        let query = BlastCompoSequenceData {
            buffer: vec![1u8; 30],
            data_offset: 0,
            length: 30,
        };
        let q_range = BlastCompoSequenceRange {
            begin: 0,
            end: 30,
            context: 0,
        };
        let subject = vec![1u8; 30];
        let s_range = BlastCompoSequenceRange {
            begin: 0,
            end: 30,
            context: 0,
        };
        let align =
            BlastCompoAlignment::new(100, MatrixAdjustRule::DontAdjust, 0, 0, 30, 0, 30, 0, None);

        let (_q_seq, s_seq) = sequence_get_range_in_memory_with_code_for_sw_redo(
            crate::program::BLASTP,
            &query,
            &q_range,
            &subject,
            &s_range,
            &crate::util::STANDARD_GENETIC_CODE,
            &[],
            &align,
            false,
            CompoAdjustMode::CompositionMatrixAdjust,
        )
        .expect("SW redo range");

        assert!(s_seq
            .data()
            .iter()
            .any(|&aa| aa == crate::encoding::NCBISTDAA_X));
    }

    #[test]
    fn sequence_get_range_for_sw_redo_skips_seg_for_near_identical_subject() {
        let query = BlastCompoSequenceData {
            buffer: vec![1u8; 30],
            data_offset: 0,
            length: 30,
        };
        let q_range = BlastCompoSequenceRange {
            begin: 0,
            end: 30,
            context: 0,
        };
        let subject = vec![1u8; 30];
        let s_range = BlastCompoSequenceRange {
            begin: 0,
            end: 30,
            context: 0,
        };
        let align =
            BlastCompoAlignment::new(100, MatrixAdjustRule::DontAdjust, 0, 0, 30, 0, 30, 0, None);

        let (_q_seq, s_seq) = sequence_get_range_in_memory_with_code_for_sw_redo(
            crate::program::BLASTP,
            &query,
            &q_range,
            &subject,
            &s_range,
            &crate::util::STANDARD_GENETIC_CODE,
            &[],
            &align,
            true,
            CompoAdjustMode::CompositionMatrixAdjust,
        )
        .expect("near-identical SW redo range");

        assert!(!s_seq
            .data()
            .iter()
            .any(|&aa| aa == crate::encoding::NCBISTDAA_X));
    }

    #[test]
    fn sequence_get_translated_range_rejects_invalid_frame() {
        let range = BlastCompoSequenceRange {
            begin: 0,
            end: 1,
            context: 0,
        };
        let err =
            sequence_get_translated_range(&[1, 8, 4], &range, &crate::util::STANDARD_GENETIC_CODE)
                .expect_err("invalid frame");
        assert_eq!(err, "invalid translated subject frame");
    }

    #[test]
    fn gapping_params_new_picks_min_lambda() {
        let scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
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
        let params = s_gapping_params_new(
            &scoring, &kbp, 2, /* bits = */ 25.0, /* raw = */ 0,
        );
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
            },
            1.0,
        );
        let params = s_gapping_params_new(&scoring, &[], 0, 25.0, 47);
        assert_eq!(params.x_dropoff, 47);
    }

    #[test]
    fn matrix_info_init_blastp_populates_blosum62() {
        let mut info = BlastMatrixInfo::default();
        let rc = s_matrix_info_init(&mut info, "BLOSUM62", 0.3176, 1.0);
        assert_eq!(rc, 0);
        assert_eq!(info.matrix_name, "BLOSUM62");
        assert!(!info.positional);
        assert_eq!(info.rows, crate::matrix::AA_SIZE as i32);
        assert_eq!(info.cols, crate::matrix::AA_SIZE as i32);
        assert_eq!(info.bit_scale_factor, 2);
        assert_eq!(info.ungapped_lambda, 0.3176);
        assert_eq!(info.matrix.len(), crate::matrix::AA_SIZE);
        assert_eq!(info.matrix[0].len(), crate::matrix::AA_SIZE);
        assert!(info.scaled_matrix.is_empty());
    }

    #[test]
    fn matrix_info_init_blastp_unknown_matrix_returns_error() {
        let mut info = BlastMatrixInfo::default();
        let rc = s_matrix_info_init(&mut info, "UNKNOWN", 0.34, 1.0);
        // NCBI returns non-zero from `s_GetStartFreqRatios` for unknown
        // matrix names.
        assert_eq!(rc, -1);
    }

    #[test]
    fn matrix_info_init_blastp_supports_pam30_freq_ratios() {
        let mut info = BlastMatrixInfo::default();
        let rc = s_matrix_info_init(&mut info, "PAM30", 0.34, 1.0);
        assert_eq!(rc, 0);
        assert_eq!(info.start_freq_ratios[1][1], 7.78912912);
        assert_eq!(info.bit_scale_factor, 2);
        assert_eq!(info.rows, crate::matrix::AA_SIZE as i32);
        assert_eq!(info.cols, crate::matrix::AA_SIZE as i32);
    }

    #[test]
    fn matrix_info_init_blastp_carries_blosum45_bit_scale() {
        let mut info = BlastMatrixInfo::default();
        let rc = s_matrix_info_init(&mut info, "BLOSUM45", 0.34, 1.0);
        assert_eq!(rc, 0);
        assert_eq!(info.start_freq_ratios[1][1], 2.95043377);
        assert_eq!(info.bit_scale_factor, 3);
    }

    #[test]
    fn record_and_restore_round_trip() {
        let mut saved = BlastKappaSavedParameters::s_saved_parameters_new(
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
            .map(|r| {
                (0..crate::matrix::AA_SIZE)
                    .map(|c| (r * 100 + c) as i32)
                    .collect()
            })
            .collect();
        let mut scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
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
    fn record_and_restore_position_based_round_trips_pssm_rows() {
        let query_length = 4;
        let mut saved = BlastKappaSavedParameters::s_saved_parameters_new(
            query_length,
            1,
            CompoAdjustMode::CompositionBasedStats,
            true,
        );
        let kbp = vec![crate::stat::KarlinBlk {
            lambda: 0.21,
            k: 0.03,
            log_k: 0.03_f64.ln(),
            h: 0.11,
            round_down: false,
        }];
        let matrix: Vec<Vec<i32>> = (0..query_length as usize)
            .map(|r| {
                (0..crate::matrix::AA_SIZE)
                    .map(|c| (1000 + r * 100 + c) as i32)
                    .collect()
            })
            .collect();
        let mut scoring = crate::parameters::ScoringParameters::from_options(
            &crate::options::ScoringOptions {
                matrix_path: None,
                reward: 0,
                penalty: 0,
                gap_open: 9,
                gap_extend: 2,
                shift_pen: i16::MAX as i32,
                gapped_calculation: true,
                complexity_adjusted_scoring: false,
                matrix_name: Some("BLOSUM62".to_string()),
                is_ooframe: false,
                program_number: crate::program::UNDEFINED,
            },
            1.0,
        );

        record_initial_search(
            &mut saved,
            &kbp,
            &matrix,
            &scoring,
            query_length,
            CompoAdjustMode::CompositionBasedStats,
            true,
        );
        assert_eq!(saved.orig_matrix.len(), query_length as usize);
        assert_eq!(saved.orig_matrix[3][7], 1307);

        scoring.gap_open = 99;
        scoring.gap_extend = 99;
        let mut live_kbp = vec![crate::stat::KarlinBlk::default()];
        let mut live_matrix = vec![vec![-1i32; crate::matrix::AA_SIZE]; 6];
        restore_search(
            &saved,
            &mut live_kbp,
            &mut live_matrix,
            &mut scoring,
            query_length,
            true,
            CompoAdjustMode::CompositionBasedStats,
        );

        assert_eq!(scoring.gap_open, 9);
        assert_eq!(scoring.gap_extend, 2);
        assert_eq!(live_kbp[0].lambda, 0.21);
        assert_eq!(live_matrix[0][0], 1000);
        assert_eq!(live_matrix[3][7], 1307);
        assert_eq!(live_matrix[4][7], -1);
    }

    #[test]
    fn saved_parameters_new_no_composition_skips_matrix() {
        let sp = BlastKappaSavedParameters::s_saved_parameters_new(
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
        let sp = BlastKappaSavedParameters::s_saved_parameters_new(
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
        let sp = BlastKappaSavedParameters::s_saved_parameters_new(
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
        let mut slot = Some(BlastKappaSavedParameters::s_saved_parameters_new(
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
        assert!(heap.insert(a).is_none());
        assert!(heap.insert(b).is_none());
        assert!(heap.insert(c).is_none());

        // Pop order: largest e-value first.
        let first = heap.pop_worst().expect("first pop");
        assert_eq!(first.oid, 2);
        let second = heap.pop_worst().expect("second pop");
        assert_eq!(second.oid, 3);
        let third = heap.pop_worst().expect("third pop");
        assert_eq!(third.oid, 1);
        assert!(heap.pop_worst().is_none());
    }

    fn heap_hsp_list(oid: i32, evalue: f64, score: i32) -> HspList {
        let mut list = HspList::new(oid);
        list.best_evalue = evalue;
        list.add_hsp(Hsp {
            score,
            num_ident: 0,
            bit_score: 0.0,
            evalue,
            query_offset: 0,
            query_end: 1,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: 1,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        list.best_evalue = evalue;
        list
    }

    #[test]
    fn blast_compo_heap_insert_replaces_only_worse_full_record() {
        let mut heap = BlastCompoHeap::new(2, 1e-5);
        assert!(heap.insert(heap_hsp_list(1, 1e-3, 10)).is_none());
        assert!(heap.insert(heap_hsp_list(2, 1e-2, 20)).is_none());

        assert!(!heap.would_insert(1e-1, 100, 99));
        let discarded = heap
            .insert(heap_hsp_list(99, 1e-1, 100))
            .expect("candidate discarded");
        assert_eq!(discarded.oid, 99);

        assert!(heap.would_insert(5e-4, 5, 3));
        let discarded = heap
            .insert(heap_hsp_list(3, 5e-4, 5))
            .expect("old worst discarded");
        assert_eq!(discarded.oid, 2);
        assert!(heap
            .records()
            .iter()
            .any(|record| record.subject_index == 3));
    }

    #[test]
    fn blast_compo_heap_ties_use_score_then_subject_index() {
        let mut heap = BlastCompoHeap::new(1, 0.0);
        assert!(heap.insert(heap_hsp_list(1, 1e-3, 10)).is_none());

        assert!(heap.would_insert(1e-3, 11, 0));
        let discarded = heap
            .insert(heap_hsp_list(0, 1e-3, 11))
            .expect("lower score record discarded");
        assert_eq!(discarded.oid, 1);

        assert!(heap.would_insert(1e-3, 11, 5));
        let discarded = heap
            .insert(heap_hsp_list(5, 1e-3, 11))
            .expect("lower subject-index record discarded");
        assert_eq!(discarded.oid, 0);

        assert!(!heap.would_insert(1e-3, 11, 2));
    }

    #[test]
    fn blast_compo_heap_keeps_all_records_once_filled_to_cutoff() {
        let mut heap = BlastCompoHeap::new(1, 1e-2);
        assert!(heap.insert(heap_hsp_list(1, 1e-3, 10)).is_none());
        assert!(heap.filled_to_cutoff());
        assert!(heap.insert(heap_hsp_list(2, 5e-3, 5)).is_none());
        assert_eq!(heap.records().len(), 2);
    }

    #[test]
    fn blast_compo_early_termination_requires_all_heaps_filled_to_cutoff() {
        let mut filled = BlastCompoHeap::new(1, 1e-3);
        assert!(filled.insert(heap_hsp_list(1, 1e-4, 10)).is_none());
        let unfilled = BlastCompoHeap::new(1, 1e-3);
        assert!(!blast_compo_early_termination(1e-2, &[filled, unfilled], 2));
    }

    #[test]
    fn blast_compo_early_termination_keeps_near_cutoff_candidates() {
        let mut heap = BlastCompoHeap::new(1, 1e-3);
        assert!(heap.insert(heap_hsp_list(1, 1e-4, 10)).is_none());
        assert!(!blast_compo_early_termination(5e-3, &[heap], 1));
    }

    #[test]
    fn blast_compo_early_termination_returns_true_past_stretched_cutoff() {
        let mut heap = BlastCompoHeap::new(1, 1e-3);
        assert!(heap.insert(heap_hsp_list(1, 1e-4, 10)).is_none());
        assert!(blast_compo_early_termination(6e-3, &[heap], 1));
    }

    #[test]
    fn blast_compo_early_termination_respects_num_queries() {
        let mut filled = BlastCompoHeap::new(1, 1e-3);
        assert!(filled.insert(heap_hsp_list(1, 1e-4, 10)).is_none());
        let extra_unfilled = BlastCompoHeap::new(1, 1e-3);

        assert!(blast_compo_early_termination(
            6e-3,
            &[filled.clone(), extra_unfilled],
            1
        ));
        assert!(!blast_compo_early_termination(6e-3, &[filled], 2));
    }

    #[test]
    fn clear_heap_drains_records() {
        let mut heap = BlastCompoHeap::new(10, 1e-5);
        let mut a = HspList::new(1);
        a.best_evalue = 1e-10;
        assert!(heap.insert(a).is_none());
        assert!(!heap.records().is_empty());
        clear_heap(&mut heap);
        assert!(heap.records().is_empty());
    }

    #[test]
    fn fill_results_from_compo_heaps_reverses_to_best_first() {
        let mut heaps = vec![BlastCompoHeap::new(10, 1e-5)];
        let mut a = HspList::new(1);
        a.best_evalue = 1e-10;
        let mut b = HspList::new(2);
        b.best_evalue = 1.0;
        assert!(heaps[0].insert(a).is_none());
        assert!(heaps[0].insert(b).is_none());
        let results = fill_results_from_compo_heaps(&mut heaps);
        let hl = results.hitlists[0].as_ref().expect("hitlist");
        // Best-evalue (smallest) should be at the head after reverse.
        assert_eq!(hl.hsp_lists.len(), 2);
        assert_eq!(hl.hsp_lists[0].oid, 1); // best evalue first
        assert_eq!(hl.hsp_lists[1].oid, 2);
    }

    #[test]
    fn merge_compo_thread_heaps_uses_destination_heap_retention_and_order() {
        let mut global = vec![BlastCompoHeap::new(2, 0.0)];
        let mut worker_heaps = vec![
            {
                let mut heaps = vec![BlastCompoHeap::new(2, 0.0)];
                assert!(heaps[0].insert(heap_hsp_list(10, 1e-2, 10)).is_none());
                assert!(heaps[0].insert(heap_hsp_list(20, 1e-4, 8)).is_none());
                heaps
            },
            {
                let mut heaps = vec![BlastCompoHeap::new(2, 0.0)];
                assert!(heaps[0].insert(heap_hsp_list(30, 1e-3, 100)).is_none());
                assert!(heaps[0].insert(heap_hsp_list(40, 1e-5, 5)).is_none());
                heaps
            },
        ];

        merge_compo_thread_heaps(&mut global, &mut worker_heaps);

        assert!(worker_heaps
            .iter()
            .all(|worker| worker.iter().all(|heap| heap.records().is_empty())));
        let results = fill_results_from_compo_heaps(&mut global);
        let hitlist = results.hitlists[0].as_ref().expect("query hitlist");
        let oids: Vec<i32> = hitlist.hsp_lists.iter().map(|list| list.oid).collect();
        assert_eq!(oids, vec![40, 20]);
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
                    segment_flags: crate::queryinfo::E_NO_SEGMENTS,
                },
                crate::queryinfo::ContextInfo {
                    query_offset: 11,
                    query_length: 10,
                    eff_searchsp: 67890,
                    length_adjustment: 0,
                    query_index: 1,
                    frame: 1,
                    is_valid: true,
                    segment_flags: crate::queryinfo::E_NO_SEGMENTS,
                },
            ],
            max_length: 10,
            min_length: 0,
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
        assert!(result[0].composition.prob.iter().all(|&p| p == 0.0));
        assert_eq!(result[0].composition.num_true_amino_acids, 0);

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
        // Verify each filled slot matches a fresh `s_get_hash` of the
        // corresponding 8-mer. NCBI's loop fills indices [0..seq_len -
        // word_size), so we check that range.
        let upper = seq.len() - 8; // exclusive
        for i in 0..upper {
            let expected = s_get_hash(&seq[i..i + 8], 8);
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
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: 8,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
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
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: 9,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        };
        let ops = vec![
            (crate::gapinfo::GapAlignOpType::Sub, 4),
            (crate::gapinfo::GapAlignOpType::Del, 1),
            (crate::gapinfo::GapAlignOpType::Sub, 4),
        ];
        let (n_ident, align_length, _) =
            blast_hsp_get_num_identities(&query, &subject, &hsp, Some(&ops), None);
        assert_eq!(n_ident, 8); // 4 A's + 4 C's all match across the gap
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
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: 8,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        compute_num_identities_blastp(&query, &subject, &mut list, &[None], None);
        assert_eq!(list.hsps[0].num_ident, 8);
    }

    #[test]
    fn compute_num_identities_translated_subject_uses_subject_frame() {
        // NCBI4na: ATG GCT translates to M A in frame +1.
        let query = vec![12u8, 1];
        let subject = vec![1u8, 8, 4, 4, 2, 8];
        let mut list = HspList::new(0);
        list.add_hsp(Hsp {
            score: 20,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 0.0,
            query_offset: 0,
            query_end: 2,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: 2,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 0,
            subject_frame: 1,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });

        compute_num_identities_translated_subject(
            &query,
            &subject,
            &mut list,
            &[None],
            None,
            &crate::util::STANDARD_GENETIC_CODE,
        )
        .expect("translated subject");
        assert_eq!(list.hsps[0].num_ident, 2);
    }

    #[test]
    fn compute_num_identities_translated_subject_supports_minus_frame() {
        // NCBI4na: ATG GCT reverse-complements to AGC CAT, so frame -1
        // translates to S H.
        let query = vec![17u8, 8];
        let subject = vec![1u8, 8, 4, 4, 2, 8];
        let mut list = HspList::new(0);
        list.add_hsp(Hsp {
            score: 20,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 0.0,
            query_offset: 0,
            query_end: 2,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: 2,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 0,
            subject_frame: -1,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });

        compute_num_identities_translated_subject(
            &query,
            &subject,
            &mut list,
            &[None],
            None,
            &crate::util::STANDARD_GENETIC_CODE,
        )
        .expect("minus frame");
        assert_eq!(list.hsps[0].num_ident, 2);
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
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: 100,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        list.add_hsp(Hsp {
            score: 50,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1.0, // above default threshold
            query_offset: 0,
            query_end: 50,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: 50,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        let (best_score, best_evalue) = s_hitlist_evaluate_and_purge(
            &mut list,
            1000, // subject_length
            crate::program::BLASTP,
            100,   // query_length
            10,    // length_adjustment
            1e9,   // eff_searchsp
            None,  // caller already populated e-values
            None,  // no link-HSP context
            -1.0,  // pvalue_for_this_pair: out of [0,1] → skip composition
            0.1,   // max_evalue
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
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: 50,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        let (best_score, best_evalue) = s_hitlist_evaluate_and_purge(
            &mut list,
            1000,
            crate::program::BLASTP,
            100,
            10,
            1e9,
            None,
            None,
            -1.0,
            0.01, // strict threshold
            false,
        );
        assert!(list.hsps.is_empty());
        assert_eq!(best_score, 0);
        assert_eq!(best_evalue, f64::MAX);
    }

    #[test]
    fn hitlist_evaluate_and_purge_recomputes_single_hsp_evalues_when_kbp_given() {
        let mut list = HspList::new(0);
        list.add_hsp(Hsp {
            score: 70,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 999.0,
            query_offset: 0,
            query_end: 20,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: 20,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        let kbp = crate::stat::KarlinBlk {
            lambda: 0.267,
            k: 0.041,
            log_k: 0.041_f64.ln(),
            h: 0.14,
            round_down: false,
        };
        let search_space = 10_000.0;
        let expected_evalue = kbp.raw_to_evalue(70, search_space);

        let (best_score, best_evalue) = s_hitlist_evaluate_and_purge(
            &mut list,
            1000,
            crate::program::BLASTP,
            100,
            10,
            search_space,
            Some(&kbp),
            None,
            -1.0,
            10.0,
            false,
        );

        assert_eq!(best_score, 70);
        assert_eq!(best_evalue, expected_evalue);
        assert_eq!(list.best_evalue, expected_evalue);
        assert_eq!(list.hsps[0].evalue, expected_evalue);
        assert_eq!(list.hsps[0].bit_score, 0.0);
    }

    #[test]
    fn blast_hsp_list_get_evalues_recomputes_all_hsps_and_best() {
        let mut list = HspList::new(0);
        for (score, old_evalue, old_bits) in [(60, 999.0, 0.0), (80, 500.0, 1.0), (40, 100.0, 2.0)]
        {
            list.add_hsp(Hsp {
                score,
                num_ident: 0,
                bit_score: old_bits,
                evalue: old_evalue,
                query_offset: 0,
                query_end: 20,
                query_gapped_start: 0,
                subject_offset: 0,
                subject_end: 20,
                subject_gapped_start: 0,
                context: 0,
                query_frame: 0,
                subject_frame: 0,
                num_gaps: 0,
                comp_adjustment_method: 0,
                edit_script: None,
                pat_info: None,
                map_info: None,
            });
        }
        let kbp = crate::stat::KarlinBlk {
            lambda: 0.267,
            k: 0.041,
            log_k: 0.041_f64.ln(),
            h: 0.14,
            round_down: false,
        };
        let search_space = 25_000.0;

        blast_hsp_list_get_evalues_simple(&mut list, &kbp, search_space);

        let expected: Vec<(f64, f64)> = [60, 80, 40]
            .into_iter()
            .map(|score| {
                (
                    kbp.raw_to_evalue(score, search_space),
                    kbp.raw_to_bit(score),
                )
            })
            .collect();
        for ((hsp, (evalue, _bit_score)), old_bits) in
            list.hsps.iter().zip(expected.iter()).zip([0.0, 1.0, 2.0])
        {
            assert_eq!(hsp.evalue, *evalue);
            assert_eq!(hsp.bit_score, old_bits);
        }
        assert_eq!(
            list.best_evalue,
            expected
                .iter()
                .map(|(evalue, _)| *evalue)
                // NCBI `s_BlastGetBestEvalue` (`blast_hits.c:1742`) seeds with `(double)INT4_MAX`.
                .fold(i32::MAX as f64, f64::min)
        );
        blast_hsp_list_get_bit_scores_simple(&mut list, &kbp);
        for (hsp, (_, bit_score)) in list.hsps.iter().zip(expected.iter()) {
            assert_eq!(hsp.bit_score, *bit_score);
        }
    }

    #[test]
    #[should_panic(expected = "do_sum_stats requires Blast_LinkHsps context")]
    fn hitlist_evaluate_and_purge_requires_link_context_for_sum_stats() {
        let mut list = HspList::new(0);
        list.add_hsp(Hsp {
            score: 50,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1.0,
            query_offset: 0,
            query_end: 50,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: 50,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });

        let _ = s_hitlist_evaluate_and_purge(
            &mut list,
            1000,
            crate::program::BLASTP,
            100,
            10,
            1e9,
            None,
            None,
            -1.0,
            10.0,
            true,
        );
    }

    #[test]
    fn blast_hsp_list_get_evalues_matches_c_context_rounding_and_decay() {
        let mut list = HspList::new(0);
        list.add_hsp(Hsp {
            score: 9,
            num_ident: 0,
            bit_score: 123.0,
            evalue: 999.0,
            query_offset: 0,
            query_end: 20,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: 20,
            subject_gapped_start: 0,
            context: 1,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });

        let kbp0 = crate::stat::KarlinBlk {
            lambda: 0.30,
            k: 0.10,
            log_k: 0.10_f64.ln(),
            h: 0.20,
            round_down: false,
        };
        let kbp1 = crate::stat::KarlinBlk {
            lambda: 0.50,
            k: 0.20,
            log_k: 0.20_f64.ln(),
            h: 0.40,
            round_down: false,
        };
        let mut sbp = crate::stat::blast_score_blk_new(crate::encoding::BLASTAA_SEQ_CODE, 2)
            .expect("score block");
        sbp.kbp_gap = vec![kbp0, kbp1.clone()];
        sbp.round_down = true;

        let query_info = crate::queryinfo::QueryInfo {
            num_queries: 1,
            contexts: vec![
                crate::queryinfo::ContextInfo {
                    query_offset: 0,
                    query_length: 10,
                    eff_searchsp: 1_000,
                    length_adjustment: 0,
                    query_index: 0,
                    frame: 0,
                    is_valid: true,
                    segment_flags: crate::queryinfo::E_NO_SEGMENTS,
                },
                crate::queryinfo::ContextInfo {
                    query_offset: 11,
                    query_length: 20,
                    eff_searchsp: 2_000,
                    length_adjustment: 0,
                    query_index: 0,
                    frame: 0,
                    is_valid: true,
                    segment_flags: crate::queryinfo::E_NO_SEGMENTS,
                },
            ],
            max_length: 20,
            min_length: 0,
        };

        let status = blast_hsp_list_get_evalues(
            crate::program::BLASTP,
            &query_info,
            100,
            &mut list,
            true,
            false,
            &sbp,
            0.5,
            2.0,
        );

        let mut scaled_kbp = kbp1.clone();
        scaled_kbp.lambda /= 2.0;
        scaled_kbp.round_down = false;
        let expected =
            scaled_kbp.raw_to_evalue(8, 2_000.0) / crate::stat::blast_gap_decay_divisor(0.5, 1);
        assert_eq!(status, 0);
        assert!((list.hsps[0].evalue - expected).abs() < 1e-12);
        assert_eq!(list.hsps[0].bit_score, 123.0);
        assert_eq!(list.best_evalue, list.hsps[0].evalue);
    }

    #[test]
    fn blast_link_hsps_for_kappa_preserves_frames_and_updates_best_evalue() {
        let mut list = HspList::new(7);
        list.add_hsp(Hsp {
            score: 50,
            num_ident: 10,
            bit_score: 20.0,
            evalue: 1e-4,
            query_offset: 0,
            query_end: 10,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: 10,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 1,
            subject_frame: -2,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        let mut linked_script = crate::gapinfo::GapEditScript::new();
        linked_script.push(crate::gapinfo::GapAlignOpType::Sub, 12);
        list.add_hsp(Hsp {
            score: 70,
            num_ident: 12,
            bit_score: 30.0,
            evalue: 1e-8,
            query_offset: 20,
            query_end: 32,
            query_gapped_start: 24,
            subject_offset: 40,
            subject_end: 52,
            subject_gapped_start: 44,
            context: 0,
            query_frame: 1,
            subject_frame: -2,
            num_gaps: 0,
            comp_adjustment_method: CompoAdjustMode::CompositionBasedStats as i32,
            edit_script: Some(linked_script),
            pat_info: None,
            map_info: None,
        });
        let qi = crate::queryinfo::QueryInfo {
            num_queries: 1,
            contexts: vec![crate::queryinfo::ContextInfo {
                query_offset: 0,
                query_length: 100,
                eff_searchsp: 10_000,
                length_adjustment: 0,
                query_index: 0,
                frame: 1,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: 100,
            min_length: 0,
        };
        let score_block = crate::link_hsps::LinkScoreBlock {
            kbp: vec![crate::stat::KarlinBlk {
                lambda: 0.267,
                k: 0.041,
                log_k: 0.041_f64.ln(),
                h: 0.14,
                round_down: false,
            }],
            kbp_gap: vec![],
            ..Default::default()
        };
        let rc = blast_link_hsps_for_kappa(
            &mut list,
            crate::program::BLASTP,
            200,
            &qi,
            0,
            &score_block,
            &crate::link_hsps::LinkHSPParameters::default(),
            false,
        );

        assert_eq!(rc, 0);
        assert_eq!(list.oid, 7);
        assert!(list.best_evalue > 1e-8);
        assert_eq!(list.hsps[0].score, 70);
        assert_eq!(list.hsps[0].evalue, list.best_evalue);
        assert_eq!(list.hsps[0].query_frame, 1);
        assert_eq!(list.hsps[0].subject_frame, -2);
        assert_eq!(list.hsps[0].query_gapped_start, 24);
        assert_eq!(list.hsps[0].subject_gapped_start, 44);
        assert_eq!(
            list.hsps[0].comp_adjustment_method,
            CompoAdjustMode::CompositionBasedStats as i32
        );
        assert_eq!(list.hsps[0].edit_script.as_ref().unwrap().len(), 1);
        assert_eq!(list.hsps[1].subject_frame, -2);
    }

    #[test]
    fn hitlist_evaluate_and_purge_takes_link_hsp_branch_when_context_supplied() {
        let mut list = HspList::new(7);
        list.add_hsp(Hsp {
            score: 50,
            num_ident: 10,
            bit_score: 20.0,
            evalue: 1e-4,
            query_offset: 0,
            query_end: 10,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: 10,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        list.add_hsp(Hsp {
            score: 70,
            num_ident: 12,
            bit_score: 30.0,
            evalue: 1e-8,
            query_offset: 20,
            query_end: 32,
            query_gapped_start: 20,
            subject_offset: 40,
            subject_end: 52,
            subject_gapped_start: 40,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        let qi = crate::queryinfo::QueryInfo {
            num_queries: 1,
            contexts: vec![crate::queryinfo::ContextInfo {
                query_offset: 0,
                query_length: 100,
                eff_searchsp: 10_000,
                length_adjustment: 0,
                query_index: 0,
                frame: 0,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: 100,
            min_length: 0,
        };
        let score_block = crate::link_hsps::LinkScoreBlock {
            kbp: vec![crate::stat::KarlinBlk {
                lambda: 0.267,
                k: 0.041,
                log_k: 0.041_f64.ln(),
                h: 0.14,
                round_down: false,
            }],
            kbp_gap: vec![],
            ..Default::default()
        };
        let link_params = crate::link_hsps::LinkHSPParameters::default();
        let link_context = HitlistLinkContext {
            query_info: &qi,
            query_context: 0,
            score_block: &score_block,
            link_params: &link_params,
            gapped_calculation: false,
        };

        let (best_score, best_evalue) = s_hitlist_evaluate_and_purge(
            &mut list,
            200,
            crate::program::BLASTP,
            100,
            0,
            10_000.0,
            None,
            Some(&link_context),
            -1.0,
            10.0,
            true,
        );

        assert_eq!(best_score, 70);
        assert!(best_evalue > 1e-8);
        assert_eq!(best_evalue, list.best_evalue);
        assert_eq!(list.hsps[0].score, 70);
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
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: 100,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        adjust_evalues_for_composition(
            &mut list, 0.5,  // comp_p_value
            1000, // subject_length
            100,  // query_length
            10,   // length_adjustment
            1e9,  // eff_searchsp
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
            .map(|i| s_get_hash(&q[i..i + 8], 8))
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
            .map(|i| s_get_hash(&q[i..i + 8], 8))
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
        let align =
            BlastCompoAlignment::new(12, MatrixAdjustRule::DontAdjust, 0, 0, 12, 0, 12, 0, None);
        assert!(!test_near_identical(&sd, 0, &qd, 0, &q_words, &align));
    }

    #[test]
    fn hsp_list_from_distinct_alignments_consumes_list() {
        // Build a 3-element linked list manually, in REVERSE-of-computation
        // order as NCBI stores them: third-best at head, best-scoring last.
        let third =
            BlastCompoAlignment::new(10, MatrixAdjustRule::DontAdjust, 0, 0, 10, 0, 10, -1, None);
        let mut second = BlastCompoAlignment::new(
            50,
            MatrixAdjustRule::ScaleOldMatrix,
            0,
            20,
            30,
            20,
            30,
            -2,
            None,
        );
        second.next = Some(Box::new(third));
        let mut first_script = crate::gapinfo::GapEditScript::new();
        first_script.push(crate::gapinfo::GapAlignOpType::Sub, 10);
        first_script.push(crate::gapinfo::GapAlignOpType::Del, 2);
        first_script.push(crate::gapinfo::GapAlignOpType::Sub, 8);
        let mut first = BlastCompoAlignment::new(
            100,
            MatrixAdjustRule::UnconstrainedRelEntropy,
            0,
            40,
            60,
            40,
            60,
            3,
            Some(BlastCompoAlignmentContext {
                edit_script: Some(first_script),
                hsp: None,
            }),
        )
        .with_hsp_context(47, 49);
        first.next = Some(Box::new(second));
        let mut head: Option<Box<BlastCompoAlignment>> = Some(Box::new(first));
        let mut hsp_list = HspList::new(0);
        let (status, tags) = hsp_list_from_distinct_alignments(&mut hsp_list, &mut head, 42, 1);
        assert_eq!(status, 0);
        assert!(head.is_none(), "linked list consumed");
        assert_eq!(hsp_list.oid, 42);
        assert_eq!(hsp_list.hsps.len(), 3);
        // After sort_by_score, highest score is first.
        assert_eq!(hsp_list.hsps[0].score, 100);
        assert_eq!(hsp_list.hsps[1].score, 50);
        assert_eq!(hsp_list.hsps[2].score, 10);
        assert_eq!(hsp_list.hsps[0].query_frame, 1);
        assert_eq!(hsp_list.hsps[0].subject_frame, 3);
        // C `s_HSPListFromDistinctAlignments` calls `Blast_HSPInit` with
        // `unknown_value` for both gapped-start arguments; it does not carry
        // the redo seed starts into the new HSP.
        assert_eq!(hsp_list.hsps[0].query_gapped_start, 0);
        assert_eq!(hsp_list.hsps[0].subject_gapped_start, 0);
        assert_eq!(hsp_list.hsps[1].subject_frame, -2);
        assert_eq!(hsp_list.hsps[2].subject_frame, -1);
        assert_eq!(
            hsp_list.hsps[0].comp_adjustment_method,
            CompoAdjustMode::CompositionMatrixAdjust as i32
        );
        assert_eq!(
            hsp_list.hsps[1].comp_adjustment_method,
            CompoAdjustMode::CompositionBasedStats as i32
        );
        assert_eq!(
            hsp_list.hsps[2].comp_adjustment_method,
            CompoAdjustMode::NoCompositionBasedStats as i32
        );
        assert_eq!(hsp_list.hsps[0].edit_script.as_ref().unwrap().len(), 3);
        assert_eq!(hsp_list.hsps[0].num_gaps, 1);
        assert!(hsp_list.hsps[1].edit_script.is_none());
        // Tags emitted in pre-sort order (matching alignment list traversal).
        assert_eq!(tags.len(), 3);
        assert_eq!(tags[0], CompoAdjustMode::CompositionMatrixAdjust); // first/UnconstrainedRelEntropy
        assert_eq!(tags[1], CompoAdjustMode::CompositionBasedStats); // second/ScaleOldMatrix
        assert_eq!(tags[2], CompoAdjustMode::NoCompositionBasedStats); // third/DontAdjust
    }

    fn linked_scores(head: &Option<Box<BlastCompoAlignment>>) -> Vec<i32> {
        let mut out = Vec::new();
        let mut cursor = head.as_deref();
        while let Some(node) = cursor {
            out.push(node.score);
            cursor = node.next.as_deref();
        }
        out
    }

    #[test]
    fn with_distinct_ends_discards_equal_or_lower_duplicate() {
        let old =
            BlastCompoAlignment::new(80, MatrixAdjustRule::DontAdjust, 0, 10, 30, 50, 70, 1, None);
        let mut head = Some(Box::new(old));
        let new_equal_end =
            BlastCompoAlignment::new(80, MatrixAdjustRule::DontAdjust, 0, 10, 35, 50, 75, 1, None);

        assert!(!with_distinct_ends(new_equal_end, &mut head, true));
        assert_eq!(linked_scores(&head), vec![80]);
    }

    #[test]
    fn with_distinct_ends_replaces_lower_duplicate_and_preserves_tail_order() {
        let mut second = BlastCompoAlignment::new(
            70,
            MatrixAdjustRule::DontAdjust,
            0,
            40,
            60,
            80,
            100,
            1,
            None,
        );
        second.next = Some(Box::new(BlastCompoAlignment::new(
            60,
            MatrixAdjustRule::DontAdjust,
            0,
            90,
            110,
            130,
            150,
            1,
            None,
        )));
        let mut first =
            BlastCompoAlignment::new(50, MatrixAdjustRule::DontAdjust, 0, 10, 30, 50, 70, 1, None);
        first.next = Some(Box::new(second));
        let mut head = Some(Box::new(first));
        let new_higher = BlastCompoAlignment::new(
            90,
            MatrixAdjustRule::DontAdjust,
            0,
            10,
            60,
            50,
            100,
            1,
            None,
        );

        assert!(with_distinct_ends(new_higher, &mut head, true));
        assert_eq!(linked_scores(&head), vec![90, 60]);
        let kept = head.as_ref().unwrap().next.as_ref().unwrap();
        assert_eq!((kept.query_start, kept.match_start), (90, 130));
    }

    #[test]
    fn with_distinct_ends_keeps_same_endpoint_in_different_frame() {
        let old = BlastCompoAlignment::new(
            100,
            MatrixAdjustRule::DontAdjust,
            0,
            10,
            30,
            50,
            70,
            2,
            None,
        );
        let mut head = Some(Box::new(old));
        let new_same_coords =
            BlastCompoAlignment::new(90, MatrixAdjustRule::DontAdjust, 0, 10, 30, 50, 70, 1, None);

        assert!(with_distinct_ends(new_same_coords, &mut head, true));
        assert_eq!(linked_scores(&head), vec![90, 100]);
    }

    #[test]
    fn with_distinct_ends_uses_similar_endpoint_for_mixed_adjustment_gate() {
        let old = BlastCompoAlignment::new(
            100,
            MatrixAdjustRule::DontAdjust,
            0,
            10,
            50,
            30,
            70,
            1,
            None,
        );
        let mut head = Some(Box::new(old));
        // Starts inside the old alignment on the same diagonal
        // (q - s == -20), but does not share an exact endpoint.
        let new_similar =
            BlastCompoAlignment::new(95, MatrixAdjustRule::DontAdjust, 0, 20, 60, 40, 80, 1, None);

        assert!(!with_distinct_ends(new_similar, &mut head, false));
        assert_eq!(linked_scores(&head), vec![100]);
    }

    #[test]
    fn windows_from_protein_aligns_groups_by_query_and_preserves_order() {
        let mut first_q1 =
            BlastCompoAlignment::new(10, MatrixAdjustRule::DontAdjust, 1, 0, 10, 5, 15, 0, None);
        let mut second_q0 =
            BlastCompoAlignment::new(20, MatrixAdjustRule::DontAdjust, 0, 0, 10, 20, 30, 0, None);
        let third_q1 =
            BlastCompoAlignment::new(30, MatrixAdjustRule::DontAdjust, 1, 20, 30, 40, 50, 0, None);
        second_q0.next = Some(Box::new(third_q1));
        first_q1.next = Some(Box::new(second_q0));
        let alignments = Some(Box::new(first_q1));
        let query_info = vec![
            BlastCompoQueryInfo {
                seq: BlastCompoSequenceData {
                    buffer: vec![0; 11],
                    data_offset: 0,
                    length: 11,
                },
                ..Default::default()
            },
            BlastCompoQueryInfo {
                seq: BlastCompoSequenceData {
                    buffer: vec![0; 17],
                    data_offset: 0,
                    length: 17,
                },
                ..Default::default()
            },
        ];

        let windows =
            windows_from_protein_aligns(&alignments, &query_info, 123).expect("protein windows");
        assert_eq!(windows.len(), 2);
        assert_eq!(windows[0].query_range.context, 0);
        assert_eq!(windows[0].query_range.end, 11);
        assert_eq!(windows[0].subject_range.begin, 0);
        assert_eq!(windows[0].subject_range.end, 123);
        assert_eq!(windows[0].hspcnt, 1);
        assert_eq!(linked_scores(&windows[0].align), vec![20]);

        assert_eq!(windows[1].query_range.context, 1);
        assert_eq!(windows[1].query_range.end, 17);
        assert_eq!(windows[1].hspcnt, 2);
        // NCBI prepends copies then reverses per window.
        assert_eq!(linked_scores(&windows[1].align), vec![10, 30]);
        // The source list was copied, not consumed.
        assert_eq!(linked_scores(&alignments), vec![10, 20, 30]);
    }

    #[test]
    fn windows_from_protein_aligns_rejects_bad_query_index() {
        let align =
            BlastCompoAlignment::new(10, MatrixAdjustRule::DontAdjust, 2, 0, 10, 5, 15, 0, None);
        let query_info = vec![BlastCompoQueryInfo::default()];
        let err = windows_from_protein_aligns(&Some(Box::new(align)), &query_info, 50)
            .expect_err("bad query index");
        assert_eq!(err, "alignment query index out of range");
    }

    #[test]
    fn translated_length_matches_ncbi_frame_formula() {
        assert_eq!(get_translated_length(101, 1, false), 33);
        assert_eq!(get_translated_length(101, 2, false), 33);
        assert_eq!(get_translated_length(101, 3, false), 33);
        assert_eq!(get_translated_length(100, -3, false), 32);
        // Position-based lengths are packed by GET_NUCL_LENGTH first.
        let packed = 2 * (101 - 2) + 5;
        assert_eq!(get_translated_length(packed, 1, true), 33);
        assert_eq!(get_translated_length(packed, -3, true), 33);
    }

    #[test]
    fn distinct_alignments_sort_matches_ncbi_tie_order() {
        let mut a =
            BlastCompoAlignment::new(80, MatrixAdjustRule::DontAdjust, 0, 10, 20, 30, 45, 0, None);
        let mut b =
            BlastCompoAlignment::new(90, MatrixAdjustRule::DontAdjust, 0, 40, 50, 80, 90, 0, None);
        let c =
            BlastCompoAlignment::new(80, MatrixAdjustRule::DontAdjust, 0, 5, 30, 30, 50, 0, None);
        b.next = Some(Box::new(c));
        a.next = Some(Box::new(b));
        let mut head = Some(Box::new(a));

        distinct_alignments_sort(&mut head);

        let mut values = Vec::new();
        let mut cur = head.as_deref();
        while let Some(node) = cur {
            values.push((
                node.score,
                node.match_start,
                node.match_end,
                node.query_start,
                node.query_end,
            ));
            cur = node.next.as_deref();
        }
        assert_eq!(
            values,
            vec![
                (90, 80, 90, 40, 50),
                (80, 30, 50, 5, 30),
                (80, 30, 45, 10, 20)
            ]
        );
    }

    #[test]
    fn windows_from_translated_aligns_joins_overlapping_subject_windows() {
        let mut first =
            BlastCompoAlignment::new(50, MatrixAdjustRule::DontAdjust, 0, 0, 10, 20, 30, 1, None);
        let second =
            BlastCompoAlignment::new(70, MatrixAdjustRule::DontAdjust, 0, 30, 40, 34, 42, 1, None);
        first.next = Some(Box::new(second));
        let alignments = Some(Box::new(first));
        let query_info = vec![BlastCompoQueryInfo {
            seq: BlastCompoSequenceData {
                buffer: vec![0; 100],
                data_offset: 0,
                length: 100,
            },
            ..Default::default()
        }];

        let windows =
            windows_from_translated_aligns(&alignments, &query_info, 2, 5, 150, true, false)
                .expect("translated windows");

        assert_eq!(windows.len(), 1);
        assert_eq!(windows[0].subject_range.begin, 15);
        assert_eq!(windows[0].subject_range.end, 47);
        assert_eq!(windows[0].subject_range.context, 1);
        assert_eq!(windows[0].query_range.begin, 0);
        assert_eq!(windows[0].query_range.end, 100);
        assert_eq!(windows[0].query_range.context, 0);
        assert_eq!(windows[0].hspcnt, 2);
        // Join appends current to next, then s_DistinctAlignmentsSort sorts by score.
        assert_eq!(linked_scores(&windows[0].align), vec![70, 50]);
        assert_eq!(linked_scores(&alignments), vec![50, 70]);
    }

    #[test]
    fn windows_from_translated_aligns_swaps_ranges_for_query_translation() {
        let align = BlastCompoAlignment::new(
            60,
            MatrixAdjustRule::DontAdjust,
            0,
            12,
            24,
            40,
            55,
            -2,
            None,
        );
        let query_info = vec![BlastCompoQueryInfo {
            seq: BlastCompoSequenceData {
                buffer: vec![0; 90],
                data_offset: 0,
                length: 90,
            },
            ..Default::default()
        }];

        let windows = windows_from_aligns(
            &Some(Box::new(align)),
            &query_info,
            1,
            1,
            4,
            200,
            true,
            false,
            false,
        )
        .expect("blastx windows");

        assert_eq!(windows.len(), 1);
        assert_eq!(windows[0].query_range.begin, 8);
        assert_eq!(windows[0].query_range.end, 28);
        assert_eq!(windows[0].query_range.context, 0);
        assert_eq!(windows[0].subject_range.begin, 0);
        assert_eq!(windows[0].subject_range.end, 200);
        assert_eq!(windows[0].subject_range.context, 0);
    }

    #[test]
    fn redo_seed_offset_rejects_seed_outside_fetched_range() {
        let range = BlastCompoSequenceRange {
            begin: 40,
            end: 80,
            context: 1,
        };

        assert_eq!(local_seed_offset(40, &range, 40), Some(0));
        assert_eq!(local_seed_offset(79, &range, 40), Some(39));
        assert_eq!(local_seed_offset(39, &range, 40), None);
        assert_eq!(local_seed_offset(80, &range, 40), None);
    }

    #[test]
    fn get_composition_range_uses_stop_margin_like_ncbi() {
        let mut seq = vec![1u8; 120];
        seq[9] = NCBISTDAA_STOP_CHAR;
        seq[100] = NCBISTDAA_STOP_CHAR;

        assert_eq!(get_composition_range(&seq, 40, 70), (30, 80));
        // Stops too close to the HSP boundary do not move the boundary.
        assert_eq!(get_composition_range(&seq, 20, 90), (20, 90));
    }

    #[test]
    fn get_composition_uses_full_range_for_non_translated_sequence() {
        let seq = BlastCompoSequenceData {
            buffer: vec![1, 1, 3, crate::encoding::NCBISTDAA_STOP, 4],
            data_offset: 0,
            length: 5,
        };
        let range = BlastCompoSequenceRange {
            begin: 0,
            end: 5,
            context: 0,
        };
        let align =
            BlastCompoAlignment::new(10, MatrixAdjustRule::DontAdjust, 0, 1, 3, 1, 3, 0, None);

        let (composition, num_true) =
            get_composition(&seq, &range, &align, crate::matrix::AA_SIZE, false, false);
        assert_eq!(num_true, 4);
        assert_eq!(composition[1], 0.5);
        assert_eq!(composition[3], 0.25);
        assert_eq!(composition[4], 0.25);
        assert_eq!(composition[crate::encoding::NCBISTDAA_STOP as usize], 0.0);
    }

    #[test]
    fn get_composition_trims_translated_range_at_stop_margin() {
        let mut data = vec![1u8; 120];
        data[9] = NCBISTDAA_STOP_CHAR;
        data[100] = NCBISTDAA_STOP_CHAR;
        data[30] = 3;
        data[79] = 4;
        let seq = BlastCompoSequenceData {
            buffer: data,
            data_offset: 0,
            length: 120,
        };
        let range = BlastCompoSequenceRange {
            begin: 0,
            end: 120,
            context: 1,
        };
        let align =
            BlastCompoAlignment::new(10, MatrixAdjustRule::DontAdjust, 0, 0, 10, 40, 70, 1, None);

        let (composition, num_true) =
            get_composition(&seq, &range, &align, crate::matrix::AA_SIZE, false, true);
        assert_eq!(num_true, 50);
        assert_eq!(composition[3], 1.0 / 50.0);
        assert_eq!(composition[4], 1.0 / 50.0);
        assert_eq!(composition[1], 48.0 / 50.0);
    }

    fn query_info_with_length(length: i32) -> Vec<BlastCompoQueryInfo> {
        vec![BlastCompoQueryInfo {
            seq: BlastCompoSequenceData {
                buffer: vec![0; length.max(0) as usize],
                data_offset: 0,
                length,
            },
            ..Default::default()
        }]
    }

    #[test]
    fn preliminary_test_near_identical_positive_cutoff_uses_score_density() {
        let query_info = query_info_with_length(80);
        let window = BlastCompoWindowInfo::new(0, 120, 0, 0, 80, 0, None);
        let high = BlastCompoAlignment::new(
            120,
            MatrixAdjustRule::DontAdjust,
            0,
            10,
            70,
            20,
            80,
            0,
            None,
        );
        let low =
            BlastCompoAlignment::new(90, MatrixAdjustRule::DontAdjust, 0, 10, 70, 20, 80, 0, None);

        assert!(preliminary_test_near_identical(
            &query_info,
            &window,
            &high,
            1.9
        ));
        assert!(!preliminary_test_near_identical(
            &query_info,
            &window,
            &low,
            1.9
        ));
    }

    #[test]
    fn preliminary_test_near_identical_positive_cutoff_requires_min_length() {
        let query_info = query_info_with_length(80);
        let window = BlastCompoWindowInfo::new(0, 120, 0, 0, 80, 0, None);
        let short =
            BlastCompoAlignment::new(200, MatrixAdjustRule::DontAdjust, 0, 0, 40, 10, 58, 0, None);

        assert!(!preliminary_test_near_identical(
            &query_info,
            &window,
            &short,
            1.0
        ));
    }

    #[test]
    fn preliminary_test_near_identical_legacy_requires_single_ungapped_equal_span() {
        let query_info = query_info_with_length(80);
        let mut single = BlastCompoWindowInfo::new(0, 120, 0, 0, 80, 0, None);
        single.hspcnt = 1;
        let mut multi = BlastCompoWindowInfo::new(0, 120, 0, 0, 80, 0, None);
        multi.hspcnt = 2;
        let equal =
            BlastCompoAlignment::new(50, MatrixAdjustRule::DontAdjust, 0, 0, 60, 10, 70, 0, None);
        let unequal =
            BlastCompoAlignment::new(50, MatrixAdjustRule::DontAdjust, 0, 0, 60, 10, 80, 0, None);

        assert!(preliminary_test_near_identical(
            &query_info,
            &single,
            &equal,
            0.0
        ));
        assert!(!preliminary_test_near_identical(
            &query_info,
            &multi,
            &equal,
            0.0
        ));
        assert!(!preliminary_test_near_identical(
            &query_info,
            &single,
            &unequal,
            0.0
        ));
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
        let qr = BlastCompoSequenceRange {
            begin: 100,
            end: 200,
            context: 7,
        };
        let sr = BlastCompoSequenceRange {
            begin: 1000,
            end: 2000,
            context: 3,
        };
        let mut script = Some(crate::gapinfo::GapEditScript::new());
        let result =
            new_alignment_from_gap_align(&tb, &mut script, &qr, &sr, MatrixAdjustRule::DontAdjust)
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

    #[test]
    fn sequence_data_release_resets_buffer_data_and_length() {
        let mut sd = BlastCompoSequenceData {
            buffer: vec![0, 1, 2, 3, 4],
            data_offset: 1,
            length: 4,
        };

        sequence_data_release(&mut sd);
        assert!(sd.buffer.is_empty());
        assert_eq!(sd.data_offset, 0);
        assert_eq!(sd.length, 0);
        assert!(sd.data().is_empty());
    }

    #[test]
    fn evalue_from_score_matches_ncbi_formula() {
        let score = 50;
        let lambda = 0.267;
        let log_k = 0.041_f64.ln();
        let searchsp = 12345.0;
        let expected = searchsp * (-(lambda * score as f64) + log_k).exp();

        assert_eq!(evalue_from_score(score, lambda, log_k, searchsp), expected);
    }
}

// ───────────────────────────────────────────────────────────────────────────
// Macros (`blast_hits_priv.h`, `redo_alignment.h`).
// ───────────────────────────────────────────────────────────────────────────

/// NCBI: CONTAINED_IN_HSP (blast_hits_priv.h:68).
/// naming: Rust helper uses snake_case for the C macro.
/// True iff `c` is in `[a, b]` and `f` is in `[d, e]`.
#[inline]
pub fn contained_in_hsp(a: i32, b: i32, c: i32, d: i32, e: i32, f: i32) -> bool {
    a <= c && b >= c && d <= f && e >= f
}

/// `KAPPA_BIT_TOL` — `redo_alignment.c:955`. Number of bits by which a
/// previously-emitted alignment must score above the current candidate
/// for [`is_contained`] to report containment.
pub const KAPPA_BIT_TOL: f64 = 2.0;
pub const COMPOSITION_MARGIN: i32 = 20;
pub const NCBISTDAA_STOP_CHAR: u8 = crate::encoding::NCBISTDAA_STOP;
pub const MINIMUM_LENGTH_NEAR_IDENTICAL: i32 = 50;
/// NCBI build-time kappa gate around TBLASTN translated-subject SEG masking.
///
/// The C source guards `s_SequenceGetTranslatedRange` with
/// `if (!(KAPPA_TBLASTN_NO_SEG_SEQUENCE))`. The vendored 2.17 reference binary
/// behaves as if this macro is enabled: composition-adjusted TBLASTN redo uses
/// the translated subject window without applying the BLASTP SEG mask.
pub const KAPPA_TBLASTN_NO_SEG_SEQUENCE: bool = true;

/// NCBI: s_IsContained (redo_alignment.c:965).
/// naming: Rust helper omits the private `s_` prefix.
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

#[inline]
/// blast-rs: Local endpoint predicate factored from `s_WithDistinctEnds`; not a direct NCBI C port.
fn is_same_endpoint(new_align: &BlastCompoAlignment, align: &BlastCompoAlignment) -> bool {
    (align.query_start == new_align.query_start && align.match_start == new_align.match_start)
        || (align.query_end == new_align.query_end && align.match_end == new_align.match_end)
}

#[inline]
/// blast-rs: Local endpoint predicate factored from `s_WithDistinctEnds`; not a direct NCBI C port.
fn is_similar_endpoint(new_align: &BlastCompoAlignment, align: &BlastCompoAlignment) -> bool {
    let start_contained = contained_in_hsp(
        align.query_start,
        align.query_end,
        new_align.query_start,
        align.match_start,
        align.match_end,
        new_align.match_start,
    );
    let end_contained = contained_in_hsp(
        align.query_start,
        align.query_end,
        new_align.query_end,
        align.match_start,
        align.match_end,
        new_align.match_end,
    );

    (start_contained
        && new_align.query_start - new_align.match_start == align.query_start - align.match_start)
        || (end_contained
            && new_align.query_end - new_align.match_end == align.query_end - align.match_end)
}

/// blast-rs: Linked-list length helper for Rust-owned alignment lists; not a direct NCBI C port.
fn alignment_list_len(mut head: Option<&BlastCompoAlignment>) -> usize {
    let mut len = 0;
    while let Some(node) = head {
        len += 1;
        head = node.next.as_deref();
    }
    len
}

/// blast-rs: Clone helper that drops the Rust linked-list tail; not a direct NCBI C port.
fn alignment_copy_without_next(align: &BlastCompoAlignment) -> BlastCompoAlignment {
    let mut copied = align.clone();
    copied.next = None;
    copied
}

/// blast-rs: In-place reversal helper for Rust-owned alignment lists; not a direct NCBI C port.
fn alignments_rev(head: &mut Option<Box<BlastCompoAlignment>>) {
    let mut prev = None;
    let mut cur = head.take();
    while let Some(mut node) = cur {
        let next = node.next.take();
        node.next = prev;
        prev = Some(node);
        cur = next;
    }
    *head = prev;
}

/// NCBI: BlastCompo_AlignmentsFree (redo_alignment.c:160).
/// naming: Public Rust parity hook omits the `blast_compo_` prefix.
///
/// Rust drops each node and its owned traceback context automatically. The C
/// callback parameter is unnecessary because `GapEditScript` is owned, so this
/// parity hook only takes and clears the list head.
pub fn alignments_free(alignments: &mut Option<Box<BlastCompoAlignment>>) {
    *alignments = None;
}

/// NCBI: s_FreeEditScript (blast_kappa.c:283).
pub fn s_free_edit_script(
    edit_script: Option<crate::gapinfo::GapEditScript>,
) -> Option<crate::gapinfo::GapEditScript> {
    crate::gapinfo::gap_edit_script_delete(edit_script)
}

/// blast-rs: Array cleanup helper for repeated `BlastCompo_AlignmentsFree`
/// loops in `Blast_RedoOneMatch*`; not a direct NCBI C port.
pub fn alignments_free_array(
    alignments: &mut [Option<Box<BlastCompoAlignment>>],
    num_queries: usize,
) {
    for slot in alignments.iter_mut().take(num_queries) {
        alignments_free(slot);
    }
}

/// NCBI: s_AlignmentCmp (redo_alignment.c:763).
/// naming: Rust comparator returns `Ordering` and omits the private `s_` prefix.
fn alignment_cmp(a: &BlastCompoAlignment, b: &BlastCompoAlignment) -> std::cmp::Ordering {
    b.score
        .cmp(&a.score)
        .then(a.match_start.cmp(&b.match_start))
        .then(b.match_end.cmp(&a.match_end))
        .then(a.query_start.cmp(&b.query_start))
        .then(b.query_end.cmp(&a.query_end))
}

/// NCBI: s_DistinctAlignmentsSort (redo_alignment.c:774).
/// naming: Rust helper omits the private `s_` prefix.
fn distinct_alignments_sort(head: &mut Option<Box<BlastCompoAlignment>>) {
    let mut nodes = Vec::new();
    let mut cur = head.take();
    while let Some(mut node) = cur {
        cur = node.next.take();
        nodes.push(node);
    }
    nodes.sort_by(|a, b| alignment_cmp(a, b));

    let mut out = None;
    for mut node in nodes.into_iter().rev() {
        node.next = out;
        out = Some(node);
    }
    *head = out;
}

/// NCBI: s_WithDistinctEnds (redo_alignment.c:395).
/// naming: Public Rust helper omits the private `s_` prefix.
///
/// The C routine conditionally inserts `new_align` at the head of a linked
/// list. If an old same-frame alignment has a matching endpoint and an equal
/// or higher score, the new alignment is discarded. Otherwise the new
/// alignment becomes the head, and any old same-frame alignments sharing an
/// exact start or end endpoint are removed while the remaining tail order is
/// preserved.
///
/// When `is_same_adjustment` is false, NCBI uses the looser
/// `s_IsSimilarEndPoint` only for the score gate because exact subjects can
/// extend farther than segged subjects; the removal pass still uses exact
/// endpoint equality. Returns true when `new_align` was inserted.
pub fn with_distinct_ends(
    mut new_align: BlastCompoAlignment,
    old_alignments: &mut Option<Box<BlastCompoAlignment>>,
    is_same_adjustment: bool,
) -> bool {
    new_align.next = None;

    let mut cursor = old_alignments.as_deref();
    while let Some(align) = cursor {
        if align.frame == new_align.frame {
            let shares_endpoint = if is_same_adjustment {
                is_same_endpoint(&new_align, align)
            } else {
                is_similar_endpoint(&new_align, align)
            };
            if shares_endpoint && new_align.score <= align.score {
                return false;
            }
        }
        cursor = align.next.as_deref();
    }

    let mut kept: Vec<Box<BlastCompoAlignment>> = Vec::new();
    let mut old = old_alignments.take();
    while let Some(mut align) = old {
        old = align.next.take();
        if align.frame == new_align.frame && is_same_endpoint(&new_align, &align) {
            continue;
        }
        kept.push(align);
    }

    let mut tail = None;
    for mut align in kept.into_iter().rev() {
        align.next = tail;
        tail = Some(align);
    }
    new_align.next = tail;
    *old_alignments = Some(Box::new(new_align));
    true
}

/// NCBI: s_SubjectCompareWindows (redo_alignment.c:620).
/// naming: Rust comparator returns `Ordering` and omits the private `s_` prefix.
fn subject_compare_windows(
    a: &BlastCompoWindowInfo,
    b: &BlastCompoWindowInfo,
) -> std::cmp::Ordering {
    a.subject_range
        .begin
        .cmp(&b.subject_range.begin)
        .then(a.subject_range.end.cmp(&b.subject_range.end))
        .then(a.subject_range.context.cmp(&b.subject_range.context))
        .then(a.query_range.begin.cmp(&b.query_range.begin))
        .then(a.query_range.end.cmp(&b.query_range.end))
        .then(a.query_range.context.cmp(&b.query_range.context))
}

/// NCBI: s_LocationCompareWindows (redo_alignment.c:643).
/// naming: Rust comparator returns `Ordering` and omits the private `s_` prefix.
fn location_compare_windows(
    a: &BlastCompoWindowInfo,
    b: &BlastCompoWindowInfo,
) -> std::cmp::Ordering {
    a.query_range
        .context
        .cmp(&b.query_range.context)
        .then(a.subject_range.context.cmp(&b.subject_range.context))
        .then(a.subject_range.begin.cmp(&b.subject_range.begin))
        .then(a.subject_range.end.cmp(&b.subject_range.end))
        .then(a.query_range.begin.cmp(&b.query_range.begin))
        .then(a.query_range.end.cmp(&b.query_range.end))
}

/// NCBI: s_GetTranslatedLength (redo_alignment.c:657).
/// naming: Public Rust helper omits the private `s_` prefix.
pub fn get_translated_length(length: i32, frame: i32, is_pos_based: bool) -> i32 {
    if is_pos_based {
        let seq_frame = if frame < 0 {
            frame.abs() + 2
        } else {
            frame - 1
        };
        let nucl_length = get_nucl_length(length);
        (nucl_length - seq_frame % 3) / 3
    } else {
        (length - frame.abs() + 1) / 3
    }
}

/// NCBI: s_WindowsFromTranslatedAligns (redo_alignment.c:676).
/// naming: Public Rust helper omits the private `s_` prefix.
///
/// Creates one bordered window per translated HSP, sorts by
/// `s_LocationCompareWindows`, joins overlapping windows with identical query
/// and subject contexts, swaps ranges for query-translated searches (blastx),
/// sorts each window's alignment list by `s_AlignmentCmp`, then returns the
/// subject-sorted window array.
pub fn windows_from_translated_aligns(
    alignments: &Option<Box<BlastCompoAlignment>>,
    query_info: &[BlastCompoQueryInfo],
    hspcnt: i32,
    border: i32,
    sequence_length: i32,
    subject_is_translated: bool,
    is_pos_based: bool,
) -> Result<Vec<BlastCompoWindowInfo>, String> {
    if hspcnt <= 0 {
        return Err("no translated alignment windows".to_string());
    }
    let mut windows = Vec::new();
    let mut align = alignments.as_deref();
    for _ in 0..hspcnt {
        let Some(node) = align else {
            return Err("translated alignment count exceeds list length".to_string());
        };
        let query_index = node.query_index;
        if query_index < 0 || query_index as usize >= query_info.len() {
            return Err("alignment query index out of range".to_string());
        }
        let query_index = query_index as usize;
        let query_length = query_info[query_index].seq.length;
        let translated_length = get_translated_length(sequence_length, node.frame, is_pos_based);
        let align_copy = Some(Box::new(alignment_copy_without_next(node)));

        let window = if subject_is_translated {
            let begin = (node.match_start - border).max(0);
            let end = (node.match_end + border).min(translated_length);
            BlastCompoWindowInfo::new(
                begin,
                end,
                node.frame,
                0,
                query_length,
                query_index as i32,
                align_copy,
            )
        } else {
            let begin = (node.query_start - border).max(0);
            let end = (node.query_end + border).min(query_length);
            BlastCompoWindowInfo::new(
                begin,
                end,
                query_index as i32,
                0,
                sequence_length,
                0,
                align_copy,
            )
        };
        windows.push(window);
        align = node.next.as_deref();
    }

    if windows.is_empty() {
        return Err("no translated alignment windows".to_string());
    }

    windows.sort_by(location_compare_windows);
    let mut joined: Vec<BlastCompoWindowInfo> = Vec::new();
    let k = 0;
    while k < windows.len() {
        if k + 1 < windows.len()
            && windows[k].subject_range.context == windows[k + 1].subject_range.context
            && windows[k].query_range.context == windows[k + 1].query_range.context
            && windows[k].subject_range.end >= windows[k + 1].subject_range.begin
        {
            let current = windows.remove(k);
            windows[k].join(current);
        } else {
            joined.push(windows.remove(k));
        }
    }

    if !subject_is_translated {
        for window in &mut joined {
            window.swap_range();
        }
    }
    for window in &mut joined {
        distinct_alignments_sort(&mut window.align);
    }
    joined.sort_by(subject_compare_windows);
    Ok(joined)
}

/// NCBI: s_WindowsFromAligns (redo_alignment.c:884).
/// naming: Public Rust helper omits the private `s_` prefix.
pub fn windows_from_aligns(
    alignments: &Option<Box<BlastCompoAlignment>>,
    query_info: &[BlastCompoQueryInfo],
    hspcnt: i32,
    num_queries: usize,
    border: i32,
    sequence_length: i32,
    query_is_translated: bool,
    subject_is_translated: bool,
    is_pos_based: bool,
) -> Result<Vec<BlastCompoWindowInfo>, String> {
    if subject_is_translated || query_is_translated {
        windows_from_translated_aligns(
            alignments,
            query_info,
            hspcnt,
            border,
            sequence_length,
            subject_is_translated,
            is_pos_based,
        )
    } else {
        windows_from_protein_aligns(
            alignments,
            &query_info[..num_queries.min(query_info.len())],
            sequence_length,
        )
    }
}

/// NCBI: Blast_GetCompositionRange (composition_adjustment.c:1236).
/// naming: Public Rust helper omits the `blast_` prefix.
pub fn get_composition_range(subject_data: &[u8], start: i32, finish: i32) -> (i32, i32) {
    let length = subject_data.len() as i32;
    let mut left = start.clamp(0, length);
    let finish = finish.clamp(0, length);

    let mut i = left;
    while i > 0 {
        if subject_data[(i - 1) as usize] == NCBISTDAA_STOP_CHAR {
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
        if subject_data[i as usize] == NCBISTDAA_STOP_CHAR {
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

/// NCBI: s_GetComposition (redo_alignment.c:930).
/// naming: Public Rust helper omits the private `s_` prefix.
pub fn get_composition(
    seq: &BlastCompoSequenceData,
    range: &BlastCompoSequenceRange,
    align: &BlastCompoAlignment,
    alphsize: usize,
    query_is_translated: bool,
    subject_is_translated: bool,
) -> (Vec<f64>, usize) {
    let data = seq.data();
    let length = (range.end - range.begin).max(0) as usize;
    let data = &data[..length.min(data.len())];
    let (left, right) = if query_is_translated || subject_is_translated {
        let start = if query_is_translated {
            align.query_start
        } else {
            align.match_start
        } - range.begin;
        let end = if query_is_translated {
            align.query_end
        } else {
            align.match_end
        } - range.begin;
        get_composition_range(data, start, end)
    } else {
        (0, data.len() as i32)
    };
    let left = left.max(0) as usize;
    let right = right.max(left as i32) as usize;
    crate::composition::blast_read_aa_composition(&data[left..right.min(data.len())], alphsize)
}

/// NCBI: s_EvalueFromScore (redo_alignment.c:976).
/// naming: Public Rust helper omits the private `s_` prefix.
pub fn evalue_from_score(score: i32, lambda: f64, log_k: f64, searchsp: f64) -> f64 {
    searchsp * (-(lambda * score as f64) + log_k).exp()
}

/// blast-rs: Significance decision factored from `Blast_RedoOneMatchSmithWaterman`
/// (`redo_alignment.c:1488`); not a direct NCBI C port.
///
/// Mirrors the C branch after `Blast_SmithWatermanScoreOnly`: link-HSP mode
/// saves by raw score cutoff; non-link mode saves by e-value cutoff, and for
/// the first alignment in a query it must also pass the significant-match heap
/// retention test.
pub fn smith_waterman_alignment_is_significant(
    sw_score: i32,
    lambda: f64,
    log_k: f64,
    searchsp: f64,
    params: &BlastRedoAlignParams,
    current_query_alignments: Option<&BlastCompoAlignment>,
    significant_matches: &BlastCompoHeap,
    subject_index: i32,
) -> bool {
    if params.do_link_hsps {
        return sw_score >= params.cutoff_s;
    }

    let sw_evalue = evalue_from_score(sw_score, lambda, log_k, searchsp);
    let mut significant = sw_evalue < params.cutoff_e;
    if current_query_alignments.is_none() {
        significant =
            significant && significant_matches.would_insert(sw_evalue, sw_score, subject_index);
    }
    significant
}

/// NCBI: s_preliminaryTestNearIdentical (redo_alignment.c:1087).
/// naming: Public Rust helper normalizes the mixed-case C static name to snake_case.
pub fn preliminary_test_near_identical(
    query_info: &[BlastCompoQueryInfo],
    window: &BlastCompoWindowInfo,
    align: &BlastCompoAlignment,
    cutoff: f64,
) -> bool {
    let query_index = align.query_index;
    if query_index < 0 || query_index as usize >= query_info.len() {
        return false;
    }
    let query_length = query_info[query_index as usize].seq.length;

    if cutoff > 0.0 {
        if align.match_end - align.match_start + 1 < query_length.min(MINIMUM_LENGTH_NEAR_IDENTICAL)
        {
            return false;
        }
        let align_len =
            (align.query_end - align.query_start).min(align.match_end - align.match_start);
        // NCBI `redo_alignment.c:1074` doesn't guard against align_len <= 0:
        //   `if ((double)align->score / (double)align_len < cutoff) return FALSE;`
        // Dividing positive score by 0 yields +inf, and `+inf < cutoff` is
        // false, so NCBI returns TRUE for zero-length alignments. Match that
        // behavior — don't bail out on align_len <= 0.
        if align.score as f64 / (align_len as f64) < cutoff {
            return false;
        }
    } else {
        if window.hspcnt > 1 || window.hspcnt < 1 {
            return false;
        }
        if align.query_end - align.query_start != align.match_end - align.match_start {
            return false;
        }
        if align.match_end - align.match_start + 1 < query_length.min(MINIMUM_LENGTH_NEAR_IDENTICAL)
        {
            return false;
        }
    }
    true
}

/// NCBI: s_WindowsFromProteinAligns (redo_alignment.c:807).
/// naming: Public Rust helper omits the private `s_` prefix.
///
/// Protein searches create one full-subject window per query that has at
/// least one alignment. Each alignment is copied into the window for its
/// `query_index`; NCBI prepends then reverses each per-window list so the
/// final order matches the incoming alignment order. The returned windows are
/// sorted with `s_SubjectCompareWindows`.
pub fn windows_from_protein_aligns(
    alignments: &Option<Box<BlastCompoAlignment>>,
    query_info: &[BlastCompoQueryInfo],
    sequence_length: i32,
) -> Result<Vec<BlastCompoWindowInfo>, String> {
    let mut windows: Vec<Option<BlastCompoWindowInfo>> = vec![None; query_info.len()];
    let mut align = alignments.as_deref();
    while let Some(node) = align {
        let query_index = node.query_index;
        if query_index < 0 || query_index as usize >= query_info.len() {
            return Err("alignment query index out of range".to_string());
        }
        let query_index = query_index as usize;
        if windows[query_index].is_none() {
            windows[query_index] = Some(BlastCompoWindowInfo::new(
                0,
                sequence_length,
                0,
                0,
                query_info[query_index].seq.length,
                query_index as i32,
                None,
            ));
        }

        let mut copied = Box::new(alignment_copy_without_next(node));
        let window = windows[query_index].as_mut().unwrap();
        copied.next = window.align.take();
        window.align = Some(copied);
        window.hspcnt += 1;
        align = node.next.as_deref();
    }

    let mut compacted = Vec::new();
    for mut window in windows.into_iter().flatten() {
        alignments_rev(&mut window.align);
        compacted.push(window);
    }
    if compacted.is_empty() {
        return Err("no protein alignment windows".to_string());
    }
    compacted.sort_by(subject_compare_windows);
    Ok(compacted)
}

/// NCBI: GET_NUCL_LENGTH (redo_alignment.h:495).
/// naming: Rust helper uses snake_case for the C macro.
///
/// Used by RPS-tblastn
/// to convert a packed mixed-frame protein length back to nucleotide
/// length.
#[inline]
pub fn get_nucl_length(l: i32) -> i32 {
    (l - 5) / 2 + 2
}

// ───────────────────────────────────────────────────────────────────────────
// Leaf utilities — no `BlastCompo_*` plumbing required.
// ───────────────────────────────────────────────────────────────────────────

/// NCBI: s_GetSubjectLength (blast_kappa.c:364).
///
/// For RPS-tblastn the subject is a packed mixed-frame protein sequence;
/// the underlying nucleotide length is recovered via `GET_NUCL_LENGTH`
/// and divided by 3. For all other programs (including ordinary RPS-blastp)
/// the input length is returned unchanged.
pub fn s_get_subject_length(total_subj_length: i32, program_number: ProgramType) -> i32 {
    // NCBI checks `eBlastTypeRpsTblastn` only — NOT `eBlastTypeRpsBlast`.
    // We mirror that here using the distinct `RPS_TBLASTN` constant.
    if program_number == RPS_TBLASTN {
        (get_nucl_length(total_subj_length) - 1) / 3
    } else {
        total_subj_length
    }
}

/// NCBI: s_HSPListNormalizeScores (blast_kappa.c:101).
///
/// Rescales each HSP score from a high-precision representation back to
/// integer scores by dividing by `score_divisor`, then sets the
/// `bit_score` from the rescaled score. Caller is responsible for
/// preserving the existing sort order (no reorder happens here).
pub fn s_hsp_list_normalize_scores(
    hsp_list: &mut HspList,
    lambda: f64,
    log_k: f64,
    score_divisor: f64,
) {
    for hsp in &mut hsp_list.hsps {
        // C: `hsp->score = (Int4)BLAST_Nint(((double) hsp->score) / scoreDivisor);`
        hsp.score = crate::math::blast_nint(hsp.score as f64 / score_divisor) as i32;
        // C: `hsp->bit_score = (hsp->score*lambda*scoreDivisor - logK)/NCBIMATH_LN2;`
        hsp.bit_score = (hsp.score as f64 * lambda * score_divisor - log_k) / NCBIMATH_LN2;
    }
}

fn s_hsp_normalize_score(hsp: &mut Hsp, lambda: f64, log_k: f64, score_divisor: f64) {
    hsp.score = crate::math::blast_nint(hsp.score as f64 / score_divisor) as i32;
    hsp.bit_score = (hsp.score as f64 * lambda * score_divisor - log_k) / NCBIMATH_LN2;
}

/// NCBI: s_GetHash (blast_kappa.c:1117).
///
/// Hash for a 28-letter-alphabet
/// word: `hash = sum_k (data[k] << (5 * (word_size - 1 - k)))`.
#[inline]
pub fn s_get_hash(data: &[u8], word_size: usize) -> u64 {
    let mut hash: u64 = 0;
    for &b in &data[..word_size] {
        hash <<= 5;
        hash += b as u64;
    }
    hash
}

/// NCBI: s_ExtendRight (blast_kappa.c:944).
///
/// Extends rightward from the start of `query` and `subject`, counting
/// identical residues and tolerating up to `max_shift` mismatches/gaps
/// when the very next two positions match. Returns
/// `(num_identical, query_ext_len, subject_ext_len, align_len)`.
pub fn s_extend_right(query: &[u8], subject: &[u8], max_shift: i32) -> (i32, i32, i32, i32) {
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
        while q_pos < query_len
            && s_pos < subject_len
            && query[q_pos as usize] == subject[s_pos as usize]
        {
            num_identical += 1;
            q_pos += 1;
            s_pos += 1;
        }

        // Try to skip mismatches or gaps.
        let mut n = 1i32;
        while n < max_shift && q_pos + n + 1 < query_len && s_pos + n + 1 < subject_len && !matched
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

/// NCBI: s_ExtendLeft (blast_kappa.c:1039).
///
/// Extends leftward from the END of `query` and `subject`. Returns
/// `(num_identical, query_ext_len, subject_ext_len, align_len_delta)`
/// where the extension lengths are measured from the end. NCBI's C
/// version takes an `align_len` in/out parameter and *adds* the delta;
/// the Rust signature returns the delta and lets the caller decide.
pub fn s_extend_left(query: &[u8], subject: &[u8], max_shift: i32) -> (i32, i32, i32, i32) {
    let query_len = query.len() as i32;
    let subject_len = subject.len() as i32;
    let mut q_pos = query_len - 1;
    let mut s_pos = subject_len - 1;
    let mut num_identical = 0i32;
    let mut gaps_in_query = 0i32;
    let mut gaps_in_subject = 0i32;

    while q_pos >= 0 && s_pos >= 0 {
        let mut matched = false;

        while q_pos > 0 && s_pos > 0 && query[q_pos as usize] == subject[s_pos as usize] {
            num_identical += 1;
            q_pos -= 1;
            s_pos -= 1;
        }

        let mut n = 1i32;
        while n < max_shift && q_pos - n - 1 > 0 && s_pos - n - 1 > 0 && !matched {
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

/// NCBI: s_FindNumIdentical (blast_kappa.c:1143).
///
/// Walks the subject
/// in a sliding-window k-mer hash, looking up each window in the query's
/// pre-hashed array; on every hit, extends bidirectionally via
/// `s_ExtendLeft` / `s_ExtendRight` to count identities.
///
/// `query_hashes` must hold one hash per starting position in the query
/// (of length `>= query.len() - WORD_SIZE`); see `create_word_array`.
pub fn s_find_num_identical(
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
            hash = s_get_hash(&subject[s_pos..s_pos + WORD_SIZE], WORD_SIZE);
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
            let (left_ident, _q_left, _s_left, _) = s_extend_left(
                &query[query_from..query_start],
                &subject[subject_from..subject_start],
                max_shift,
            );
            num_identical += left_ident;

            // Extend right.
            let (right_ident, _q_right, s_right, _) = s_extend_right(
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

/// NCBI: BlastKappa_SavedParameters (blast_kappa.c:1958).
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
    /// NCBI: s_SavedParametersNew (blast_kappa.c:2008).
    ///
    /// `rows` is the row count of the scoring matrix (PSSM length for
    /// position-based searches; `BLASTAA_SIZE` otherwise). When
    /// `compo_adjust_mode` is `NoCompositionBasedStats`, `orig_matrix`
    /// stays empty (matching NCBI's `if (compo_adjust_mode != ...)`
    /// guard at the bottom of `s_SavedParametersNew`).
    pub fn s_saved_parameters_new(
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
        // are filled in by `record_initial_search`.
        sp.kbp_gap_orig.resize(
            num_queries.max(0) as usize,
            crate::stat::KarlinBlk::default(),
        );
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

/// NCBI: s_NewAlignmentUsingXdrop (blast_kappa.c:1812).
/// naming: Public Rust helper omits the private `s_` prefix.
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
        Some(BlastCompoAlignmentContext {
            edit_script: Some(edit_script),
            hsp: None,
        }),
    );
    Some((new_align, q_end_xdrop as i32, s_end_xdrop as i32))
}

/// blast-rs: Protein-space counterpart of [`new_alignment_using_xdrop`];
/// not a direct NCBI C port.
///
/// NCBI's Smith-Waterman redo path calls the caller-supplied X-drop callback
/// after the SW score/start/end pass. For materialized protein and translated
/// searches, the equivalent callback is the matrix-scored protein X-drop
/// traceback over the SW-bounded rectangle.
#[allow(clippy::too_many_arguments)]
pub fn new_alignment_using_xdrop_protein(
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
    matrix: &[[i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE],
) -> Option<(BlastCompoAlignment, i32, i32)> {
    let q_extent = query_end.checked_sub(query_start)?;
    let s_extent = match_end.checked_sub(match_start)?;
    let tb = crate::protein::s_sw_find_final_ends_using_xdrop(
        query.data(),
        subject.data(),
        query_start,
        match_start,
        q_extent,
        s_extent,
        score,
        matrix,
        gapping_params.gap_open,
        gapping_params.gap_extend,
        gapping_params.x_dropoff,
    )?;

    let q_end_xdrop = tb.query_end;
    let s_end_xdrop = tb.subject_end;
    let aquery_start = tb.query_start as i32 + query_range.begin;
    let aquery_end = q_end_xdrop as i32 + query_range.begin;
    let amatch_start = tb.subject_start as i32 + subject_range.begin;
    let amatch_end = s_end_xdrop as i32 + subject_range.begin;

    let new_align = BlastCompoAlignment::new(
        tb.score,
        matrix_adjust_rule,
        query_range.context,
        aquery_start,
        aquery_end,
        amatch_start,
        amatch_end,
        subject_range.context,
        Some(BlastCompoAlignmentContext {
            edit_script: Some(tb.edit_script),
            hsp: None,
        }),
    );
    Some((new_align, q_end_xdrop as i32, s_end_xdrop as i32))
}

/// blast-rs: Position-specific/PSSM counterpart of [`new_alignment_using_xdrop_protein`];
/// not a direct NCBI C port.
#[allow(clippy::too_many_arguments)]
pub fn new_alignment_using_xdrop_protein_pssm(
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
    pssm: &[Vec<i32>],
) -> Option<(BlastCompoAlignment, i32, i32)> {
    let q_extent = query_end.checked_sub(query_start)?;
    let s_extent = match_end.checked_sub(match_start)?;
    let tb = crate::protein::protein_sw_bounded_xdrop_align_pssm(
        query.data(),
        subject.data(),
        query_start,
        match_start,
        q_extent,
        s_extent,
        score,
        query_range.begin.max(0) as usize,
        pssm,
        gapping_params.gap_open,
        gapping_params.gap_extend,
        gapping_params.x_dropoff,
    )?;

    let q_end_xdrop = tb.query_end;
    let s_end_xdrop = tb.subject_end;
    let aquery_start = tb.query_start as i32 + query_range.begin;
    let aquery_end = q_end_xdrop as i32 + query_range.begin;
    let amatch_start = tb.subject_start as i32 + subject_range.begin;
    let amatch_end = s_end_xdrop as i32 + subject_range.begin;

    let new_align = BlastCompoAlignment::new(
        tb.score,
        matrix_adjust_rule,
        query_range.context,
        aquery_start,
        aquery_end,
        amatch_start,
        amatch_end,
        subject_range.context,
        Some(BlastCompoAlignmentContext {
            edit_script: Some(tb.edit_script),
            hsp: None,
        }),
    );
    Some((new_align, q_end_xdrop as i32, s_end_xdrop as i32))
}

/// NCBI: s_RedoOneAlignment (blast_kappa.c:1898).
/// naming: Public Rust helper omits the private `s_` prefix.
///
/// X-drop alignment in BOTH directions from the seed `(gapped_start_q,
/// gapped_start_s)`, producing a fresh `BlastCompoAlignment` that
/// supersedes the input alignment. NCBI calls
/// `BLAST_GappedAlignmentWithTraceback` (the bidirectional driver);
/// our Rust analog is `crate::traceback::blast_gapped_alignment_with_traceback`.
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
    let q_start = local_seed_offset(gapped_start_q, query_range, query_data.length)?;
    let s_start = local_seed_offset(gapped_start_s, subject_range, subject_data.length)?;

    let tb = crate::traceback::blast_gapped_alignment_with_traceback(
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

/// blast-rs: Converts local `BlastMatrixInfo` rows to amino-acid arrays; not a direct NCBI C port.
fn matrix_info_to_aa_array(
    matrix_info: &BlastMatrixInfo,
) -> Option<[[i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE]> {
    if matrix_info.matrix.len() < crate::matrix::AA_SIZE {
        return None;
    }
    let mut matrix = [[0i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE];
    for (i, row) in matrix.iter_mut().enumerate() {
        let src = matrix_info.matrix.get(i)?;
        if src.len() < crate::matrix::AA_SIZE {
            return None;
        }
        row.copy_from_slice(&src[..crate::matrix::AA_SIZE]);
    }
    Some(matrix)
}

fn local_seed_offset(
    gapped_start: i32,
    range: &BlastCompoSequenceRange,
    data_length: i32,
) -> Option<usize> {
    let local = gapped_start.checked_sub(range.begin)?;
    if local < 0 || local >= data_length {
        return None;
    }
    Some(local as usize)
}

/// blast-rs: Matrix-scored counterpart of [`redo_one_alignment`] for protein-space
/// searches; not a direct NCBI C port.
///
/// Matrix-scored counterpart of [`redo_one_alignment`] for protein-space
/// searches. This is the materialized-subject analog of NCBI's
/// `BLAST_GappedAlignmentWithTraceback` callback when the query and subject
/// ranges are already protein residues, including translated-subject ranges.
#[allow(clippy::too_many_arguments)]
pub fn redo_one_alignment_protein(
    query_data: &BlastCompoSequenceData,
    query_range: &BlastCompoSequenceRange,
    subject_data: &BlastCompoSequenceData,
    subject_range: &BlastCompoSequenceRange,
    gapped_start_q: i32,
    gapped_start_s: i32,
    matrix_adjust_rule: MatrixAdjustRule,
    matrix_info: &BlastMatrixInfo,
    gap_open: i32,
    gap_extend: i32,
    x_dropoff: i32,
) -> Option<BlastCompoAlignment> {
    let q_start = local_seed_offset(gapped_start_q, query_range, query_data.length)?;
    let s_start = local_seed_offset(gapped_start_s, subject_range, subject_data.length)?;
    let matrix = matrix_info_to_aa_array(matrix_info)?;
    let tb = crate::protein::protein_gapped_align(
        query_data.data(),
        subject_data.data(),
        q_start,
        s_start,
        &matrix,
        gap_open,
        gap_extend,
        x_dropoff,
    )?;
    let mut edit_script = Some(tb.edit_script.clone());
    let shifted = crate::traceback::TracebackResult {
        score: tb.score,
        query_start: tb.query_start,
        query_end: tb.query_end,
        subject_start: tb.subject_start,
        subject_end: tb.subject_end,
        edit_script: tb.edit_script,
    };
    new_alignment_from_gap_align(
        &shifted,
        &mut edit_script,
        query_range,
        subject_range,
        matrix_adjust_rule,
    )
}

/// blast-rs: Matrix-scored protein redo with caller-owned adjusted matrix; not a direct NCBI C port.
///
/// Matrix-scored protein redo when the caller already owns the adjusted
/// square amino-acid matrix for this query/subject pair.
#[allow(clippy::too_many_arguments)]
pub fn redo_one_alignment_protein_with_matrix(
    query_data: &BlastCompoSequenceData,
    query_range: &BlastCompoSequenceRange,
    subject_data: &BlastCompoSequenceData,
    subject_range: &BlastCompoSequenceRange,
    gapped_start_q: i32,
    gapped_start_s: i32,
    matrix_adjust_rule: MatrixAdjustRule,
    matrix: &[[i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE],
    gap_open: i32,
    gap_extend: i32,
    x_dropoff: i32,
) -> Option<BlastCompoAlignment> {
    let q_start = local_seed_offset(gapped_start_q, query_range, query_data.length)?;
    let s_start = local_seed_offset(gapped_start_s, subject_range, subject_data.length)?;
    let tb = crate::protein::protein_gapped_align(
        query_data.data(),
        subject_data.data(),
        q_start,
        s_start,
        matrix,
        gap_open,
        gap_extend,
        x_dropoff,
    )?;
    let mut edit_script = Some(tb.edit_script.clone());
    let shifted = crate::traceback::TracebackResult {
        score: tb.score,
        query_start: tb.query_start,
        query_end: tb.query_end,
        subject_start: tb.subject_start,
        subject_end: tb.subject_end,
        edit_script: tb.edit_script,
    };
    new_alignment_from_gap_align(
        &shifted,
        &mut edit_script,
        query_range,
        subject_range,
        matrix_adjust_rule,
    )
}

/// blast-rs: Position-specific protein redo helper; not a direct NCBI C port.
///
/// Position-specific counterpart of [`redo_one_alignment_protein_with_matrix`]
/// for PSI/PSSM searches. Scores are taken from
/// `pssm[absolute_query_position][subject_residue]`.
#[allow(clippy::too_many_arguments)]
pub fn redo_one_alignment_protein_with_pssm(
    query_data: &BlastCompoSequenceData,
    query_range: &BlastCompoSequenceRange,
    subject_data: &BlastCompoSequenceData,
    subject_range: &BlastCompoSequenceRange,
    gapped_start_q: i32,
    gapped_start_s: i32,
    matrix_adjust_rule: MatrixAdjustRule,
    pssm: &[Vec<i32>],
    gap_open: i32,
    gap_extend: i32,
    x_dropoff: i32,
) -> Option<BlastCompoAlignment> {
    let q_start = local_seed_offset(gapped_start_q, query_range, query_data.length)?;
    let s_start = local_seed_offset(gapped_start_s, subject_range, subject_data.length)?;
    let tb = crate::protein::protein_gapped_align_pssm(
        query_data.data(),
        subject_data.data(),
        q_start,
        s_start,
        query_range.begin.max(0) as usize,
        pssm,
        gap_open,
        gap_extend,
        x_dropoff,
    )?;
    let mut edit_script = Some(tb.edit_script.clone());
    let shifted = crate::traceback::TracebackResult {
        score: tb.score,
        query_start: tb.query_start,
        query_end: tb.query_end,
        subject_start: tb.subject_start,
        subject_end: tb.subject_end,
        edit_script: tb.edit_script,
    };
    new_alignment_from_gap_align(
        &shifted,
        &mut edit_script,
        query_range,
        subject_range,
        matrix_adjust_rule,
    )
}

/// NCBI: s_GetStartFreqRatios (blast_kappa.c:648).
/// naming: Public Rust helper omits the private `s_` prefix.
///
/// Returns the BLASTAA_SIZE × BLASTAA_SIZE frequency-ratio matrix for
/// the named scoring matrix. NCBI's C version dispatches via
/// `_PSIMatrixFrequencyRatiosNew`; the same standard matrix names are routed
/// through [`crate::matrix::get_matrix_freq_ratios`].
///
/// Output ratios are deep-copied (1-1 with NCBI's per-cell copy from
/// `stdFreqRatios->data[i][j]`).
pub fn get_start_freq_ratios(
    matrix_name: &str,
) -> Result<[[f64; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE], ()> {
    crate::matrix::get_matrix_freq_ratios(matrix_name).ok_or(())
}

/// NCBI: s_GetPosBasedStartFreqRatios (blast_kappa.c:591).
/// naming: Public Rust helper omits the private `s_` prefix.
///
/// Builds the position-specific frequency-ratio rows underlying a PSSM. NCBI
/// first copies the standard matrix row selected by each query residue, then
/// overlays PSI `startNumerator[i][j] / standardProb[j]` when both standard
/// probabilities and the numerator are above `kPosEpsilon`, excluding stop and
/// X columns.
pub fn get_pos_based_start_freq_ratios(
    query: &[u8],
    matrix_name: &str,
    start_numerator: &[Vec<f64>],
) -> Result<Vec<Vec<f64>>, ()> {
    const K_POS_EPSILON: f64 = 0.0001;
    const E_X_CHAR: usize = crate::encoding::NCBISTDAA_X as usize;
    const E_STOP_CHAR: usize = crate::encoding::NCBISTDAA_STOP as usize;

    if start_numerator.len() < query.len() {
        return Err(());
    }
    let std_ratios = get_start_freq_ratios(matrix_name)?;
    let standard_prob = crate::stat::protein_std_freq_ncbistdaa();
    let mut return_ratios = vec![vec![0.0f64; crate::matrix::AA_SIZE]; query.len()];

    for (i, &q) in query.iter().enumerate() {
        let qi = q as usize;
        if qi >= crate::matrix::AA_SIZE || start_numerator[i].len() < crate::matrix::AA_SIZE {
            return Err(());
        }
        for j in 0..crate::matrix::AA_SIZE {
            return_ratios[i][j] = std_ratios[qi][j];
        }
    }

    for (i, &q) in query.iter().enumerate() {
        let qi = q as usize;
        for j in 0..crate::matrix::AA_SIZE {
            if standard_prob[qi] > K_POS_EPSILON
                && standard_prob[j] > K_POS_EPSILON
                && j != E_STOP_CHAR
                && j != E_X_CHAR
                && start_numerator[i][j] > K_POS_EPSILON
            {
                return_ratios[i][j] = start_numerator[i][j] / standard_prob[j];
            }
        }
    }

    Ok(return_ratios)
}

/// blast-rs: Position-frequency-ratio to PSSM converter; not a direct NCBI C port.
///
/// Convert position-specific frequency-ratio rows to integer PSSM rows.
///
/// This is the unscaled log-odds conversion stage used before PSI's
/// `_PSIScaleMatrix` / IMPALA scaling pass. It mirrors the same
/// `ln(freq_ratio) / lambda` conversion used by `Blast_Int4MatrixFromFreq`,
/// but works on `query_length × BLASTAA_SIZE` position-specific rows.
pub fn pos_freq_ratios_to_pssm(freq_ratios: &[Vec<f64>], lambda: f64) -> Result<Vec<Vec<i32>>, ()> {
    const COMPO_SCORE_MIN: i32 = -100_000;

    if lambda <= 0.0 {
        return Err(());
    }

    let mut pssm = Vec::with_capacity(freq_ratios.len());
    for row in freq_ratios {
        if row.len() < crate::matrix::AA_SIZE {
            return Err(());
        }
        let mut out = vec![0i32; crate::matrix::AA_SIZE];
        for j in 0..crate::matrix::AA_SIZE {
            out[j] = if row[j] <= 0.0 {
                COMPO_SCORE_MIN
            } else {
                crate::math::blast_nint(row[j].ln() / lambda) as i32
            };
        }
        pssm.push(out);
    }
    Ok(pssm)
}

/// PSI private scaled PSSM factor used by NCBI's
/// `_PSIConvertFreqRatiosToPSSM`.
pub const PSI_PSSM_SCALE_FACTOR: f64 = 200.0;

/// blast-rs: Build PSI-private scaled PSSM rows; not a direct NCBI C port.
///
/// Build the private scaled PSSM rows from position-specific frequency ratios.
///
/// NCBI keeps both an integer score matrix and `scaled_pssm` in
/// `_PSIInternalPssmData`; the latter is the natural-log odds score multiplied
/// by 200 for later PSI scaling/statistics routines.
pub fn pos_freq_ratios_to_scaled_pssm(freq_ratios: &[Vec<f64>]) -> Result<Vec<Vec<i32>>, ()> {
    const COMPO_SCORE_MIN: i32 = -100_000;

    let mut scaled = Vec::with_capacity(freq_ratios.len());
    for row in freq_ratios {
        if row.len() < crate::matrix::AA_SIZE {
            return Err(());
        }
        let mut out = vec![0i32; crate::matrix::AA_SIZE];
        for j in 0..crate::matrix::AA_SIZE {
            out[j] = if row[j] <= 0.0 {
                COMPO_SCORE_MIN
            } else {
                crate::math::blast_nint(row[j].ln() * PSI_PSSM_SCALE_FACTOR) as i32
            };
        }
        scaled.push(out);
    }
    Ok(scaled)
}

/// NCBI: s_ScalePosMatrix (blast_kappa.c:676).
/// naming: Public Rust helper omits the private `s_` prefix.
///
/// NCBI allocates a temporary `SFreqRatios`, runs
/// `_PSIConvertFreqRatiosToPSSM`, then copies the resulting public PSSM into
/// `self->startMatrix`. The PSI-private `scaled_pssm` is kept here as
/// `scaled_matrix` because later PSI scaling/statistics code consumes it in
/// the C pipeline.
pub fn scale_pos_matrix(self_: &mut BlastMatrixInfo, freq_ratios: &[Vec<f64>]) -> i32 {
    let matrix = match pos_freq_ratios_to_pssm(freq_ratios, self_.ungapped_lambda) {
        Ok(matrix) => matrix,
        Err(()) => return -1,
    };
    let scaled_matrix = match pos_freq_ratios_to_scaled_pssm(freq_ratios) {
        Ok(matrix) => matrix,
        Err(()) => return -1,
    };

    self_.matrix = matrix;
    self_.scaled_matrix = scaled_matrix;
    0
}

/// blast-rs: PSI-private statistic update after PSSM scaling; not a direct NCBI C port.
///
/// PSI-private statistic update after PSSM scaling.
///
/// NCBI's PSI scaling/statistics path applies the Lambda ratio returned by the
/// composition/PSSM scaling step to the active Karlin blocks
/// (`kbp->Lambda *= lambdaRatio`). The matrix rows themselves are updated by
/// [`scale_pos_matrix`] / [`composition_scale_pssm_with_ratio`]; this helper
/// exposes the matching statistics side effect so callers can keep bit-score
/// and e-value calculations tied to the scaled PSSM.
pub fn psi_private_update_lambda_statistics(
    kbp_gap: &mut [crate::stat::KarlinBlk],
    lambda_ratio: f64,
) -> i32 {
    if !lambda_ratio.is_finite() || lambda_ratio <= 0.0 {
        return -1;
    }
    for kbp in kbp_gap {
        if kbp.lambda > 0.0 {
            kbp.lambda *= lambda_ratio;
        }
    }
    0
}

/// `SCALING_FACTOR` (`blast_kappa.c:676`). Multiplicative factor for
/// integer matrix scores — gives composition-adjusted scoring extra
/// precision without overflowing `BLAST_SCORE_MIN`.
pub const SCALING_FACTOR: f64 = 32.0;

/// NCBI: s_SWFindFinalEndsUsingXdrop (blast_kappa.c:843).
/// naming: Public Rust helper omits the private `s_` prefix.
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
) -> (
    i32,
    usize,
    usize,
    Vec<(crate::gapinfo::GapAlignOpType, i32)>,
) {
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

/// NCBI: BlastGapAlignStruct (blast_gapalign.h).
///
/// The
/// gapped-alignment workspace owned by the kappa driver.
///
/// NCBI's struct holds many heap-allocated workspace buffers
/// (`state_struct`, `dp_mem`, `greedy_align_mem`, `fwd_prelim_tback`,
/// `rev_prelim_tback`) plus the alignment output fields (`score`,
/// `query_start/stop`, `subject_start/stop`, `edit_script`).
///
/// Our Rust ports of the gap-align routines (`align_ex`,
/// `blast_gapped_alignment_with_traceback`, `gapped_score_one_dir_*`) own their workspace
/// buffers internally, so the Rust analog here only carries the
/// alignment **output** state plus the few control fields the kappa
/// driver reads or writes (`gap_x_dropoff`).
///
/// The C `s_BlastGapAlignStruct_Copy` and `_Free` helpers correspond
/// to Rust's automatic `Clone`/`Drop`; these functions preserve the C
/// ownership transitions.
#[derive(Debug, Clone, Default)]
pub struct BlastGapAlignWorkspace {
    pub position_based: bool,
    pub gap_x_dropoff: i32,
    pub max_mismatches: i32,
    pub mismatch_window: i32,
    pub score: i32,
    pub query_start: i32,
    pub query_stop: i32,
    pub subject_start: i32,
    pub subject_stop: i32,
    pub edit_script: Option<crate::gapinfo::GapEditScript>,
    pub fwd_prelim_tback: Option<crate::gapinfo::GapPrelimEditBlock>,
    pub rev_prelim_tback: Option<crate::gapinfo::GapPrelimEditBlock>,
    pub greedy_align_mem: Option<crate::greedy::SGreedyAlignMem>,
    pub dp_mem_alloc: i32,
    pub jumper: Option<crate::gapinfo::JumperGapAlign>,
    pub chaining: Option<crate::extend::ChainingStruct>,
}

impl BlastGapAlignWorkspace {
    /// NCBI: s_BlastGapAlignStruct_Copy (blast_kappa.c:2604).
    pub fn s_blast_gap_align_struct_copy(&self) -> Self {
        let copied_script = if let Some(src_script) = self.edit_script.as_ref() {
            let mut dst_script = crate::gapinfo::GapEditScript::default();
            let op_count = src_script.len();
            dst_script.reserve(op_count);
            for index in 0..op_count {
                let (op_type, count) = src_script
                    .get(index)
                    .unwrap_or((crate::gapinfo::GapAlignOpType::Sub, 0));
                dst_script.push_raw(op_type, count);
            }
            Some(dst_script)
        } else {
            None
        };

        let mut dst = BlastGapAlignWorkspace::default();
        dst.position_based = self.position_based;
        dst.gap_x_dropoff = self.gap_x_dropoff;
        dst.max_mismatches = self.max_mismatches;
        dst.mismatch_window = self.mismatch_window;
        dst.score = self.score;
        dst.query_start = self.query_start;
        dst.query_stop = self.query_stop;
        dst.subject_start = self.subject_start;
        dst.subject_stop = self.subject_stop;
        dst.edit_script = None;
        dst.fwd_prelim_tback = self.fwd_prelim_tback.clone();
        dst.rev_prelim_tback = self.rev_prelim_tback.clone();
        dst.greedy_align_mem = self.greedy_align_mem.clone();
        dst.dp_mem_alloc = self.dp_mem_alloc;
        dst.jumper = self.jumper.clone();
        dst.chaining = self.chaining.clone();
        if let Some(script) = copied_script {
            dst.edit_script = Some(script);
        }
        dst
    }
}

/// Inputs consumed by [`blast_gap_align_struct_new`].
///
/// This is the Rust ownership equivalent of the non-output arguments to NCBI
/// `BLAST_GapAlignStructNew`: scoring parameters, extension parameters, max
/// subject length, score-block position-specific state, and the output slot.
#[derive(Clone, Copy)]
pub struct BlastGapAlignNewParams<'a> {
    pub scoring: Option<&'a crate::parameters::ScoringParameters>,
    pub extension: Option<&'a crate::parameters::ExtensionParameters>,
    pub max_subject_length: u32,
    pub position_based: bool,
    pub legacy_gap_x_dropoff: i32,
}

impl From<i32> for BlastGapAlignNewParams<'_> {
    fn from(gap_x_dropoff: i32) -> Self {
        Self {
            scoring: None,
            extension: None,
            max_subject_length: 0,
            position_based: false,
            legacy_gap_x_dropoff: gap_x_dropoff,
        }
    }
}

/// NCBI: `BLAST_GapAlignStructNew` (`blast_gapalign.c:293`).
///
/// Allocates and initializes the modeled `BlastGapAlignStruct` fields used by
/// the kappa/traceback ports. When full scoring/extension parameters are
/// supplied, this follows the C branch structure exactly: DP score-only gets a
/// 1000-cell allocation, greedy preliminary extension owns a greedy workspace
/// capped by the C max-subject formula, jumper traceback owns a jumper
/// workspace and computes the C fallback x-drop when the requested x-drop is
/// zero, chaining owns a chaining workspace, and both preliminary traceback
/// blocks are created on success.
///
/// Passing a raw `i32` preserves the historical blast-rs call shape for code
/// that only had an x-drop value available.
pub fn blast_gap_align_struct_new<'a, P>(params: P) -> Option<BlastGapAlignWorkspace>
where
    P: Into<BlastGapAlignNewParams<'a>>,
{
    let params = params.into();
    let mut workspace = BlastGapAlignWorkspace::default();
    workspace.position_based = params.position_based;

    if let Some(extension) = params.extension {
        workspace.gap_x_dropoff = extension.gap_x_dropoff;
        workspace.max_mismatches = extension.options.max_mismatches;
        workspace.mismatch_window = extension.options.mismatch_window;

        if extension.options.prelim_gap_ext != crate::options::PrelimGapExt::JumperWithTraceback {
            if extension.options.prelim_gap_ext == crate::options::PrelimGapExt::DynProgScoreOnly {
                workspace.dp_mem_alloc = 1000;
            } else {
                let scoring = params.scoring?;
                let max_subject_length =
                    blast_gap_align_greedy_max_subject_length(params.max_subject_length);
                workspace.greedy_align_mem = Some(crate::greedy::s_blast_greedy_align_mem_alloc(
                    scoring.reward,
                    scoring.penalty,
                    scoring.gap_open,
                    scoring.gap_extend,
                    max_subject_length as i32,
                    extension.gap_x_dropoff,
                )?);
            }
        } else {
            let scoring = params.scoring?;
            workspace.jumper = Some(crate::gapinfo::jumper_gap_align_new(200)?);
            if extension.gap_x_dropoff == 0 {
                workspace.gap_x_dropoff =
                    3 * (-scoring.penalty).max(scoring.gap_open + scoring.gap_extend);
            }
        }

        if extension.options.chaining {
            workspace.chaining = crate::extend::chaining_struct_new();
            workspace.chaining.as_ref()?;
        }
    } else {
        workspace.gap_x_dropoff = params.legacy_gap_x_dropoff.max(0);
    }

    workspace.fwd_prelim_tback = Some(crate::gapinfo::gap_prelim_edit_block_new());
    workspace.rev_prelim_tback = Some(crate::gapinfo::gap_prelim_edit_block_new());
    Some(workspace)
}

impl<'a> From<&'a crate::parameters::ExtensionParameters> for BlastGapAlignNewParams<'a> {
    fn from(extension: &'a crate::parameters::ExtensionParameters) -> Self {
        Self {
            scoring: None,
            extension: Some(extension),
            max_subject_length: 0,
            position_based: false,
            legacy_gap_x_dropoff: 0,
        }
    }
}

impl<'a>
    From<(
        &'a crate::parameters::ScoringParameters,
        &'a crate::parameters::ExtensionParameters,
        u32,
        bool,
    )> for BlastGapAlignNewParams<'a>
{
    fn from(
        value: (
            &'a crate::parameters::ScoringParameters,
            &'a crate::parameters::ExtensionParameters,
            u32,
            bool,
        ),
    ) -> Self {
        Self {
            scoring: Some(value.0),
            extension: Some(value.1),
            max_subject_length: value.2,
            position_based: value.3,
            legacy_gap_x_dropoff: 0,
        }
    }
}

fn blast_gap_align_greedy_max_subject_length(max_subject_length: u32) -> u32 {
    const MAX_DBSEQ_LEN: u32 = 5_000_000;
    const GREEDY_MAX_COST_FRACTION: u32 = 2;
    const GREEDY_MAX_COST: u32 = 1000;

    GREEDY_MAX_COST.min(max_subject_length.min(MAX_DBSEQ_LEN) / GREEDY_MAX_COST_FRACTION + 1)
}

/// NCBI: s_BlastGapAlignStruct_Free (blast_kappa.c:2532).
///
/// over the modeled Rust workspace fields.
pub fn s_blast_gap_align_struct_free(slot: &mut Option<BlastGapAlignWorkspace>) {
    if let Some(workspace) = slot.as_mut() {
        workspace.position_based = false;
        workspace.gap_x_dropoff = 0;
        workspace.max_mismatches = 0;
        workspace.mismatch_window = 0;
        workspace.score = 0;
        workspace.query_start = 0;
        workspace.query_stop = 0;
        workspace.subject_start = 0;
        workspace.subject_stop = 0;
        workspace.edit_script =
            crate::gapinfo::gap_edit_script_delete(workspace.edit_script.take());
        workspace.fwd_prelim_tback = None;
        workspace.rev_prelim_tback = None;
        workspace.greedy_align_mem = None;
        workspace.dp_mem_alloc = 0;
        workspace.jumper = crate::gapinfo::jumper_gap_align_free(workspace.jumper.take());
        workspace.chaining = crate::extend::chaining_struct_free(&mut workspace.chaining);
    }
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
    pub sbp: Option<BlastScoreBlkSnapshot>,
    pub local_scaling_factor: f64,
    pub program: ProgramType,
}

/// Explicit subject-source selector for `Blast_RedoAlignmentCore_MT`.
///
/// C makes this choice through `subjectBlk`, `seqSrc`, `thisMatch`, and
/// `hsp_stream` parameters. Keeping that branch decision in the Rust signature
/// prevents the MT entry point from guessing subject length or silently taking
/// the wrong redo path.
pub enum BlastRedoAlignmentSource<'a> {
    /// Compatibility branch for callers that already own a redone HSP list and
    /// only need the no-composition result-list update.
    ExistingMatchOnly,
    /// Callback-backed branch corresponding to NCBI's `Blast_RedoAlignCallbacks`.
    Callback(BlastRedoCallbackSubjectConfig<'a>),
    /// One materialized subject block and one HSP list.
    InMemorySubject(BlastRedoInMemorySubject<'a>),
    /// Materialized equivalent of C's stream loop over multiple subject HSP lists.
    InMemorySubjects(&'a mut [BlastRedoInMemorySubjectMatch<'a>]),
    /// `BlastSeqSrc`-backed branch: fetch each subject by OID, then redo it.
    SeqSrcSubjects {
        seqsrc: &'a dyn crate::seqsrc::BlastSeqSource,
        matches: &'a mut [HspList],
        config: BlastRedoSeqSrcSubjectConfig<'a>,
    },
}

/// NCBI: Blast_RedoAlignmentCore (blast_kappa.c:2942).
///
/// Single-thread adapter that delegates to the multi-thread driver
/// `Blast_RedoAlignmentCore_MT(num_threads = 1, ...)`. The MT driver
/// itself is 669 LOC of dense composition-adjustment orchestration; the
/// no-composition single-match path is implemented in the MT function.
///
/// Sets `num_threads = 1` and delegates, exactly mirroring the C control
/// flow. Returns the same `Int2` status code the MT function returns.
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
    source: BlastRedoAlignmentSource<'_>,
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
        source,
        results,
    )
}

/// Port boundary for `Blast_RedoAlignmentCore_MT` (`blast_kappa.c:2980`,
/// 669 LOC, cyclomatic 120). The full driver orchestrates:
/// 1. Save initial search state via [`record_initial_search`].
/// 2. Rescale scoring via [`rescale_search`].
/// 3. Build the per-query info array via [`get_query_info`].
/// 4. Build the redo-align params via [`s_get_align_params`].
/// 5. For each subject-match in `hsp_stream`:
///    - run the redo-alignment driver against the subject's HSPs;
///    - filter via [`s_hitlist_evaluate_and_purge`];
///    - feed survivors into a per-query [`BlastCompoHeap`].
/// 6. Collapse heaps into `results` via [`fill_results_from_compo_heaps`].
/// 7. Restore initial search state via [`restore_search`].
///
/// The subject branch is explicit in [`BlastRedoAlignmentSource`], mirroring
/// C's `subjectBlk` / `seqSrc` / `thisMatch` / `hsp_stream` split. The MT
/// entry point must not infer a subject block from HSP coordinates: callers
/// either supply the subject source needed for real redo or choose
/// [`BlastRedoAlignmentSource::ExistingMatchOnly`] for the no-subject,
/// no-composition compatibility branch.
#[allow(clippy::too_many_arguments)]
pub fn blast_redo_alignment_core_mt(
    program: ProgramType,
    num_threads: u32,
    query: &[u8],
    query_info: &crate::queryinfo::QueryInfo,
    kbp_gap: &mut [crate::stat::KarlinBlk],
    matrix: &mut [Vec<i32>],
    scoring: &mut crate::parameters::ScoringParameters,
    align_params: &BlastRedoAlignParams,
    saved: &mut BlastKappaSavedParameters,
    this_match: &mut HspList,
    source: BlastRedoAlignmentSource<'_>,
    results: &mut crate::hspstream::HspResults,
) -> i32 {
    let align_params = match blast_redo_alignment_core_mt_effective_params(
        query,
        query_info,
        scoring,
        matrix,
        align_params,
    ) {
        Ok(params) => params,
        Err(status) => return status,
    };
    let align_params = &align_params;

    let query_length = query_info.max_length as i32;
    let record_status = record_initial_search(
        saved,
        kbp_gap,
        matrix,
        scoring,
        query_length,
        align_params.compo_adjust_mode,
        align_params.position_based,
    );
    if record_status != 0 {
        return record_status;
    }

    // NCBI always rescales the search immediately after recording the
    // original state (`Blast_RedoAlignmentCore_MT`, after
    // `s_RecordInitialSearch`). For no-composition searches this factor is
    // 1.0; for composition adjustment it is normally `SCALING_FACTOR`.
    rescale_search(
        kbp_gap,
        scoring,
        query_info.contexts.len() as i32,
        align_params.local_scaling_factor,
    );

    let status = match source {
        BlastRedoAlignmentSource::Callback(subject) => {
            blast_redo_alignment_core_one_match_with_callbacks_inner(
                program,
                query,
                query_info,
                kbp_gap,
                align_params,
                this_match,
                subject,
                results,
            )
        }
        BlastRedoAlignmentSource::InMemorySubject(subject) => {
            blast_redo_alignment_core_one_match_in_memory_inner(
                program,
                query,
                query_info,
                kbp_gap,
                align_params,
                this_match,
                results,
                None,
                None,
                subject,
            )
        }
        BlastRedoAlignmentSource::InMemorySubjects(matches) => {
            if matches.is_empty() {
                *results = crate::hspstream::HspResults::new(query_info.num_queries as i32);
                0
            } else if !materialized_redo_program_is_supported(program) {
                -1
            } else {
                let num_queries = query_info.num_queries.max(0) as usize;
                let heap_cfg = matches[0].subject;
                let num_workers = (num_threads.max(1) as usize).min(matches.len().max(1));
                let mut worker_heaps =
                    vec![
                        vec![
                            BlastCompoHeap::new(heap_cfg.hitlist_size, heap_cfg.inclusion_ethresh);
                            num_queries
                        ];
                        num_workers
                    ];
                let mut status = 0;

                for (item_index, item) in matches.iter_mut().enumerate() {
                    let worker_index = item_index % num_workers;
                    let heaps = &mut worker_heaps[worker_index];
                    if blast_compo_early_termination(item.hsp_list.best_evalue, heaps, num_queries)
                    {
                        continue;
                    }

                    let mut local_results = crate::hspstream::HspResults::new(num_queries as i32);
                    status = blast_redo_alignment_core_one_match_in_memory_inner(
                        program,
                        query,
                        query_info,
                        kbp_gap,
                        align_params,
                        &mut item.hsp_list,
                        &mut local_results,
                        Some(&heaps[..]),
                        None,
                        item.subject,
                    );
                    if status != 0 {
                        break;
                    }

                    for (query_index, hitlist) in local_results.hitlists.into_iter().enumerate() {
                        if query_index >= heaps.len() {
                            continue;
                        }
                        if let Some(hitlist) = hitlist {
                            for hsp_list in hitlist.hsp_lists {
                                let _discarded = heaps[query_index].insert(hsp_list);
                            }
                        }
                    }
                }

                if status == 0 {
                    let mut heaps =
                        vec![
                            BlastCompoHeap::new(heap_cfg.hitlist_size, heap_cfg.inclusion_ethresh);
                            num_queries
                        ];
                    merge_compo_thread_heaps(&mut heaps, &mut worker_heaps);
                    *results = fill_results_from_compo_heaps(&mut heaps);
                }
                status
            }
        }
        BlastRedoAlignmentSource::SeqSrcSubjects {
            seqsrc,
            matches,
            config,
        } => {
            if matches.is_empty() {
                *results = crate::hspstream::HspResults::new(query_info.num_queries as i32);
                0
            } else if !materialized_redo_program_is_supported(program) {
                -1
            } else {
                let num_queries = query_info.num_queries.max(0) as usize;
                let num_workers = (num_threads.max(1) as usize).min(matches.len().max(1));
                let mut worker_heaps =
                    vec![
                        vec![
                            BlastCompoHeap::new(config.hitlist_size, config.inclusion_ethresh);
                            num_queries
                        ];
                        num_workers
                    ];
                let encoding = seqsrc_encoding_for_redo_program(program);
                let mut status = 0;

                for (item_index, hsp_list) in matches.iter_mut().enumerate() {
                    let worker_index = item_index % num_workers;
                    let heaps = &mut worker_heaps[worker_index];
                    if blast_compo_early_termination(hsp_list.best_evalue, heaps, num_queries) {
                        continue;
                    }

                    if program == crate::program::TBLASTN {
                        let mut matching_seq = BlastCompoMatchingSequence::default();
                        if s_matching_sequence_initialize(
                            &mut matching_seq,
                            program,
                            seqsrc,
                            1,
                            hsp_list.oid,
                            None,
                        ) != 0
                        {
                            status = -1;
                            break;
                        }
                        let subject_sequence = {
                            let Some(local_data) = matching_seq.local_data.as_ref() else {
                                matching_sequence_release(&mut matching_seq);
                                status = -1;
                                break;
                            };
                            let Some(subject_blk) = local_data.seq_arg.seq.as_ref() else {
                                matching_sequence_release(&mut matching_seq);
                                status = -1;
                                break;
                            };
                            let Some(sequence) = subject_blk.sequence.as_deref() else {
                                matching_sequence_release(&mut matching_seq);
                                status = -1;
                                break;
                            };
                            if subject_blk.length < 0
                                || subject_blk.length as usize > sequence.len()
                            {
                                matching_sequence_release(&mut matching_seq);
                                status = -1;
                                break;
                            }
                            &sequence[..subject_blk.length as usize]
                        };
                        let mut local_results =
                            crate::hspstream::HspResults::new(num_queries as i32);
                        status = blast_redo_alignment_core_one_match_in_memory_inner(
                            program,
                            query,
                            query_info,
                            kbp_gap,
                            align_params,
                            hsp_list,
                            &mut local_results,
                            Some(&heaps[..]),
                            None,
                            BlastRedoInMemorySubject {
                                subject_source: subject_sequence,
                                reward: config.reward,
                                penalty: config.penalty,
                                genetic_code: config.genetic_code,
                                smith_waterman: config.smith_waterman,
                                expect_value: config.expect_value,
                                hitlist_size: config.hitlist_size,
                                inclusion_ethresh: config.inclusion_ethresh,
                                link_context: config.link_context,
                            },
                        );
                        matching_sequence_release(&mut matching_seq);
                        if status != 0 {
                            break;
                        }

                        for (query_index, hitlist) in local_results.hitlists.into_iter().enumerate()
                        {
                            if query_index >= heaps.len() {
                                continue;
                            }
                            if let Some(hitlist) = hitlist {
                                for hsp_list in hitlist.hsp_lists {
                                    let _discarded = heaps[query_index].insert(hsp_list);
                                }
                            }
                        }
                        continue;
                    }

                    let Some(seq_data) = seqsrc.get_sequence(&crate::seqsrc::GetSeqArg {
                        oid: hsp_list.oid,
                        encoding,
                        ..crate::seqsrc::GetSeqArg::default()
                    }) else {
                        status = -1;
                        break;
                    };
                    if seq_data.length < 0 || seq_data.length as usize > seq_data.sequence.len() {
                        status = -1;
                        break;
                    }
                    let subject_sequence = &seq_data.sequence[..seq_data.length as usize];
                    let mut local_results = crate::hspstream::HspResults::new(num_queries as i32);
                    status = blast_redo_alignment_core_one_match_in_memory_inner(
                        program,
                        query,
                        query_info,
                        kbp_gap,
                        align_params,
                        hsp_list,
                        &mut local_results,
                        Some(&heaps[..]),
                        None,
                        BlastRedoInMemorySubject {
                            subject_source: subject_sequence,
                            reward: config.reward,
                            penalty: config.penalty,
                            genetic_code: config.genetic_code,
                            smith_waterman: config.smith_waterman,
                            expect_value: config.expect_value,
                            hitlist_size: config.hitlist_size,
                            inclusion_ethresh: config.inclusion_ethresh,
                            link_context: config.link_context,
                        },
                    );
                    if status != 0 {
                        break;
                    }

                    for (query_index, hitlist) in local_results.hitlists.into_iter().enumerate() {
                        if query_index >= heaps.len() {
                            continue;
                        }
                        if let Some(hitlist) = hitlist {
                            for hsp_list in hitlist.hsp_lists {
                                let _discarded = heaps[query_index].insert(hsp_list);
                            }
                        }
                    }
                }

                if status == 0 {
                    let mut heaps =
                        vec![
                            BlastCompoHeap::new(config.hitlist_size, config.inclusion_ethresh);
                            num_queries
                        ];
                    merge_compo_thread_heaps(&mut heaps, &mut worker_heaps);
                    *results = fill_results_from_compo_heaps(&mut heaps);
                }
                status
            }
        }
        BlastRedoAlignmentSource::ExistingMatchOnly => {
            if align_params.callbacks.is_some()
                || !matches!(
                    align_params.compo_adjust_mode,
                    CompoAdjustMode::NoCompositionBasedStats
                )
            {
                -1
            } else if this_match.hsps.is_empty() {
                0
            } else {
                this_match.sort_by_score();
                this_match.best_evalue = this_match
                    .hsps
                    .iter()
                    .map(|hsp| hsp.evalue)
                    .fold(i32::MAX as f64, f64::min);
                let per_query =
                    s_blast_redo_alignment_group_hsps_by_query(program, query_info, this_match);
                let mut touched_queries = s_blast_redo_alignment_update_hitlists(
                    per_query,
                    align_params,
                    this_match,
                    results,
                );
                s_blast_redo_alignment_finalize_touched_hitlists(&mut touched_queries, results);
                0
            }
        }
    };

    restore_search(
        saved,
        kbp_gap,
        matrix,
        scoring,
        query_length,
        align_params.position_based,
        align_params.compo_adjust_mode,
    );

    status
}

/// NCBI: entry validation and local composition-mode adjustment at the start
/// of `Blast_RedoAlignmentCore_MT`.
fn blast_redo_alignment_core_mt_effective_params(
    query: &[u8],
    query_info: &crate::queryinfo::QueryInfo,
    scoring: &crate::parameters::ScoringParameters,
    matrix: &[Vec<i32>],
    align_params: &BlastRedoAlignParams,
) -> Result<BlastRedoAlignParams, i32> {
    let mut effective = align_params.clone();
    let matrix_name = scoring.options.matrix_name.as_deref().unwrap_or("");

    if matrix_name.eq_ignore_ascii_case("BLOSUM62_20")
        && matches!(
            effective.compo_adjust_mode,
            CompoAdjustMode::NoCompositionBasedStats
        )
    {
        return Err(-1);
    }

    if effective.position_based {
        if effective.compo_adjust_mode as u8 > CompoAdjustMode::CompositionBasedStats as u8 {
            effective.compo_adjust_mode = CompoAdjustMode::CompositionBasedStats;
        }
        if query_info.num_queries != 1 {
            return Err(-1);
        }
        if !query.is_empty() && query.len() != query_info.max_length as usize {
            return Err(-1);
        }
        if !matrix.is_empty() && matrix.len() < query_info.max_length as usize {
            return Err(-1);
        }
    }

    if effective.compo_adjust_mode as u8 > CompoAdjustMode::CompositionBasedStats as u8
        && crate::matrix::get_matrix_freq_ratios_with_scale(matrix_name).is_none()
    {
        return Err(-1);
    }

    Ok(effective)
}

/// NCBI: score ordering and best-evalue refresh for an HSP list.
fn s_hsp_list_normalize_current_best_evalue(hsp_list: &mut HspList) {
    hsp_list.sort_by_score();
    hsp_list.best_evalue = hsp_list
        .hsps
        .iter()
        .map(|hsp| hsp.evalue)
        .fold(i32::MAX as f64, f64::min);
}

/// NCBI: split subject HSPs into per-query HSP lists.
fn s_blast_redo_alignment_group_hsps_by_query(
    program: ProgramType,
    query_info: &crate::queryinfo::QueryInfo,
    this_match: &HspList,
) -> Vec<Option<crate::hspstream::HspList>> {
    let mut per_query: Vec<Option<crate::hspstream::HspList>> = Vec::new();
    for hsp in this_match.hsps.iter().cloned() {
        let context = hsp.context;
        let context_info = query_info.contexts.get(context.max(0) as usize);
        if context_info.is_some_and(|ctx| !ctx.is_valid) {
            continue;
        }
        let query_index = context_info.map(|ctx| ctx.query_index).unwrap_or_else(|| {
            if context < 0 {
                -1
            } else {
                crate::queryinfo::blast_get_query_index_from_context(context, program)
            }
        });
        if query_index < 0 {
            continue;
        }

        let query_index = query_index as usize;
        if per_query.len() <= query_index {
            per_query.resize_with(query_index + 1, || None);
        }
        let query_hsp_list = per_query[query_index].get_or_insert_with(|| {
            let mut list = crate::hspstream::HspList::new(this_match.oid);
            list.hsp_max = this_match.hsp_max;
            list
        });
        query_hsp_list.hsps.push(hsp);
    }
    per_query
}

/// NCBI: `Blast_HitListUpdate`/merge phase for redo output lists.
fn s_blast_redo_alignment_update_hitlists(
    per_query: Vec<Option<crate::hspstream::HspList>>,
    align_params: &BlastRedoAlignParams,
    this_match: &HspList,
    results: &mut crate::hspstream::HspResults,
) -> Vec<usize> {
    let mut touched_queries = Vec::new();
    for (query_index, maybe_list) in per_query.into_iter().enumerate() {
        let Some(mut query_hsp_list) = maybe_list else {
            continue;
        };
        if query_hsp_list.hsps.is_empty() {
            continue;
        }

        s_blast_redo_alignment_filter_query_hsp_list(&mut query_hsp_list, align_params);
        if query_hsp_list.hsps.is_empty() {
            continue;
        }
        if query_index >= results.hitlists.len() {
            results.hitlists.resize_with(query_index + 1, || None);
        }
        touched_queries.push(query_index);
        let hitlist =
            results.hitlists[query_index].get_or_insert_with(crate::hspstream::HitList::new);
        if let Some(existing) = hitlist
            .hsp_lists
            .iter_mut()
            .find(|list| list.oid == query_hsp_list.oid)
        {
            let hsp_num_max = existing.hsp_max;
            let mut incoming = Some(query_hsp_list);
            let mut combined = Some(std::mem::replace(
                existing,
                crate::hspstream::HspList::new(-1),
            ));
            let _ = crate::hspstream::blast_hsp_lists_merge_simple(
                &mut incoming,
                &mut combined,
                hsp_num_max,
            );
            if let Some(merged) = combined {
                *existing = merged;
            }
        } else {
            hitlist.blast_hit_list_update(query_hsp_list);
        }
        hitlist.sort_by_evalue();
    }
    if this_match.hsps.is_empty() {
        touched_queries.clear();
    }
    touched_queries
}

/// NCBI: per-query HSP-list cutoff and size filtering.
fn s_blast_redo_alignment_filter_query_hsp_list(
    query_hsp_list: &mut HspList,
    align_params: &BlastRedoAlignParams,
) {
    query_hsp_list.sort_by_score();
    if query_hsp_list.hsp_max >= 0 && query_hsp_list.hsps.len() > query_hsp_list.hsp_max as usize {
        query_hsp_list
            .hsps
            .truncate(query_hsp_list.hsp_max as usize);
    }
    query_hsp_list.hsps.retain(|hsp| {
        hsp.score >= align_params.cutoff_s
            || hsp.evalue <= align_params.cutoff_e
            || align_params.cutoff_e <= 0.0
    });
    query_hsp_list.best_evalue = query_hsp_list
        .hsps
        .iter()
        .map(|hsp| hsp.evalue)
        .fold(i32::MAX as f64, f64::min);
}

/// NCBI: final result-list cleanup and score/e-value ordering.
fn s_blast_redo_alignment_finalize_touched_hitlists(
    touched_queries: &mut Vec<usize>,
    results: &mut crate::hspstream::HspResults,
) {
    touched_queries.sort_unstable();
    touched_queries.dedup();
    for query_index in touched_queries.iter().copied() {
        let Some(hitlist_slot) = results.hitlists.get_mut(query_index) else {
            continue;
        };
        let Some(hitlist) = hitlist_slot.as_mut() else {
            continue;
        };
        hitlist.hsp_lists.retain(|list| !list.hsps.is_empty());
        for hsp_list in &mut hitlist.hsp_lists {
            s_hsp_list_normalize_current_best_evalue(hsp_list);
        }
        if hitlist.hsp_lists.is_empty() {
            *hitlist_slot = None;
        } else {
            hitlist.sort_by_evalue();
        }
    }
}

/// Inputs for the materialized-subject branch of
/// [`blast_redo_alignment_core_mt_in_memory_subject`].
///
/// This models the C `subjectBlk != NULL && hsp_stream == NULL` use case:
/// the caller already owns the subject sequence bytes, so no `BlastSeqSrc`
/// fetch/copy layer is needed.
#[derive(Clone, Copy)]
pub struct BlastRedoInMemorySubject<'a> {
    pub subject_source: &'a [u8],
    pub reward: i32,
    pub penalty: i32,
    pub genetic_code: &'a [u8; 64],
    pub smith_waterman: bool,
    pub expect_value: f64,
    pub hitlist_size: i32,
    pub inclusion_ethresh: f64,
    pub link_context: Option<&'a HitlistLinkContext<'a>>,
}

/// One already-materialized subject/HSP-list pair for
/// [`blast_redo_alignment_core_mt_in_memory_subjects`].
pub struct BlastRedoInMemorySubjectMatch<'a> {
    pub hsp_list: HspList,
    pub subject: BlastRedoInMemorySubject<'a>,
}

/// Shared per-subject settings for the `BlastSeqSource` fetching branch of
/// [`blast_redo_alignment_core_mt_seqsrc_subjects`].
#[derive(Clone, Copy)]
pub struct BlastRedoSeqSrcSubjectConfig<'a> {
    pub reward: i32,
    pub penalty: i32,
    pub genetic_code: &'a [u8; 64],
    pub smith_waterman: bool,
    pub expect_value: f64,
    pub hitlist_size: i32,
    pub inclusion_ethresh: f64,
    pub link_context: Option<&'a HitlistLinkContext<'a>>,
}

/// Inputs for the represented RPS-BLAST profile/query redo boundary.
///
/// RPS traceback enters kappa redo with the profile as the position-specific
/// query and the user's protein query as the materialized subject. The returned
/// HSP list deliberately remains in that pre-`s_BlastHSPListRPSUpdate`
/// orientation so the existing RPS traceback driver can restore coordinates in
/// the same place as NCBI.
#[derive(Clone)]
pub struct BlastRpsRedoProfileConfig<'a> {
    pub matrix_name: &'a str,
    pub kbp_gap: crate::stat::KarlinBlk,
    pub kbp_ideal_lambda: f64,
    pub scoring: crate::parameters::ScoringParameters,
    pub gapping_params: BlastCompoGappingParams,
    pub compo_adjust_mode: CompoAdjustMode,
    pub local_scaling_factor: f64,
    pub cutoff_s: i32,
    pub expect_value: f64,
    pub hitlist_size: i32,
    pub inclusion_ethresh: f64,
    pub smith_waterman: bool,
}

/// External callback-backed subject settings for
/// [`blast_redo_alignment_core_mt_with_callbacks`].
///
/// This models the `Blast_RedoAlignmentCore_MT` branch where the redo driver
/// does not own a materialized subject buffer and must retrieve query/subject
/// windows through `Blast_RedoAlignCallbacks::get_range`.
#[derive(Clone, Copy)]
pub struct BlastRedoCallbackSubjectConfig<'a> {
    pub subject_length: i32,
    pub smith_waterman: bool,
    pub expect_value: f64,
    pub hitlist_size: i32,
    pub inclusion_ethresh: f64,
    pub link_context: Option<&'a HitlistLinkContext<'a>>,
}

/// blast-rs: Program-mode gate for materialized redo helpers; not a direct NCBI C port.
fn materialized_redo_program_is_supported(program: ProgramType) -> bool {
    matches!(
        program,
        crate::program::BLASTN
            | crate::program::BLASTP
            | crate::program::BLASTX
            | crate::program::TBLASTN
            | crate::program::TBLASTX
            | crate::program::PSI_BLAST
    )
}

/// NCBI `BlastExtensionOptionsValidate` permits composition adjustment only
/// for blastp, blastx, tblastn, and RPS/PSI variants; tblastx is ungapped-only
/// and never enters the kappa CBS traceback path.
fn composition_redo_program_is_supported(program: ProgramType) -> bool {
    matches!(
        program,
        crate::program::BLASTP
            | crate::program::BLASTX
            | crate::program::TBLASTN
            | crate::program::PSI_BLAST
            | crate::program::RPS_BLAST
            | crate::program::RPS_TBLASTN
    )
}

/// blast-rs: Context grouping helper for materialized redo dispatch; not a direct NCBI C port.
fn materialized_num_contexts_per_query(
    program: ProgramType,
    query_info: &crate::queryinfo::QueryInfo,
) -> usize {
    if (program == crate::program::BLASTN || crate::program::blast_query_is_translated(program))
        && query_info.num_queries > 0
    {
        (query_info.contexts.len() / query_info.num_queries as usize).max(1)
    } else {
        1
    }
}

/// blast-rs: Maps redo program type to the Rust `BlastSeqSource` encoding; not a direct NCBI C port.
fn seqsrc_encoding_for_redo_program(program: ProgramType) -> crate::seqsrc::SeqEncoding {
    if crate::program::blast_subject_is_translated(program) {
        crate::seqsrc::SeqEncoding::Ncbi4na
    } else if program == crate::program::BLASTN {
        crate::seqsrc::SeqEncoding::Nucleotide
    } else {
        crate::seqsrc::SeqEncoding::Protein
    }
}

fn matrix_info_init_rps_profile(
    self_: &mut BlastMatrixInfo,
    matrix_name: &str,
    profile_pssm: &[Vec<i32>],
    profile_freq_ratios: Option<&[Vec<f64>]>,
    profile_length: usize,
    kbp_ideal_lambda: f64,
    scale_factor: f64,
) -> i32 {
    if profile_length == 0 || profile_pssm.len() < profile_length || scale_factor <= 0.0 {
        return -1;
    }

    let Some(freq_ratio_info) = crate::matrix::get_matrix_freq_ratios_with_scale(matrix_name)
    else {
        return -1;
    };

    let mut matrix = Vec::with_capacity(profile_length);
    for row in profile_pssm.iter().take(profile_length) {
        if row.len() < crate::matrix::AA_SIZE {
            return -1;
        }
        matrix.push(row[..crate::matrix::AA_SIZE].to_vec());
    }

    let start_freq_ratios = if let Some(freq_rows) = profile_freq_ratios {
        if freq_rows.len() < profile_length {
            return -1;
        }
        let mut ratios = Vec::with_capacity(profile_length);
        for row in freq_rows.iter().take(profile_length) {
            if row.len() < crate::matrix::AA_SIZE {
                return -1;
            }
            ratios.push(row[..crate::matrix::AA_SIZE].to_vec());
        }
        ratios
    } else {
        vec![vec![1.0; crate::matrix::AA_SIZE]; profile_length]
    };

    self_.matrix_name = matrix_name.to_string();
    self_.positional = true;
    self_.bit_scale_factor = freq_ratio_info.bit_scale_factor;
    self_.ungapped_lambda = kbp_ideal_lambda / scale_factor;
    self_.rows = profile_length as i32;
    self_.cols = crate::matrix::AA_SIZE as i32;
    self_.matrix = matrix.clone();
    self_.scaled_matrix = matrix;
    self_.start_freq_ratios = start_freq_ratios;
    0
}

/// blast-rs: RPS profile/query integration boundary for
/// `Blast_RedoAlignmentCore` / `Blast_TracebackFromHSPList`; not a complete
/// direct NCBI C port.
///
/// This wires an owned RPS profile row selection into the existing
/// position-specific composition redo core. The profile consensus and PSSM
/// rows become the kappa "query", the user's protein query is the materialized
/// subject, and the incoming HSP list is updated in place with the redone
/// gapped HSPs when they pass the configured significance gates.
pub fn blast_rps_redo_alignment_core_profile(
    profile_consensus: &[u8],
    user_query: &[u8],
    profile_pssm: &[Vec<i32>],
    profile_freq_ratios: Option<&[Vec<f64>]>,
    hsp_list: &mut HspList,
    cfg: BlastRpsRedoProfileConfig<'_>,
) -> i32 {
    let profile_length = profile_consensus.len();
    if profile_length == 0 || user_query.is_empty() || hsp_list.hsps.is_empty() {
        return if hsp_list.hsps.is_empty() { 0 } else { -1 };
    }

    let mut matrix_info = BlastMatrixInfo::default();
    if matrix_info_init_rps_profile(
        &mut matrix_info,
        cfg.matrix_name,
        profile_pssm,
        profile_freq_ratios,
        profile_length,
        cfg.kbp_ideal_lambda,
        cfg.local_scaling_factor,
    ) != 0
    {
        return -1;
    }

    let mut local_list = hsp_list.clone();
    for hsp in &mut local_list.hsps {
        if hsp.query_offset < 0
            || hsp.query_end > profile_length as i32
            || hsp.subject_offset < 0
            || hsp.subject_end > user_query.len() as i32
        {
            return -1;
        }
        hsp.context = 0;
        hsp.query_frame = 0;
        hsp.subject_frame = 0;
    }

    let query_info = crate::queryinfo::QueryInfo {
        num_queries: 1,
        contexts: vec![crate::queryinfo::ContextInfo {
            query_offset: 0,
            query_length: profile_length as i32,
            eff_searchsp: (profile_length.max(1) * user_query.len().max(1)) as i64,
            length_adjustment: 0,
            query_index: 0,
            frame: 0,
            is_valid: true,
            segment_flags: crate::queryinfo::E_NO_SEGMENTS,
        }],
        max_length: profile_length as u32,
        min_length: profile_length as u32,
    };

    let mut kbp_gap = vec![cfg.kbp_gap.clone()];
    let mut matrix = matrix_info.matrix.clone();
    let mut scoring = cfg.scoring;
    let align_params = blast_redo_align_params_new(
        matrix_info,
        cfg.gapping_params,
        cfg.compo_adjust_mode,
        cfg.local_scaling_factor,
        true,
        false,
        false,
        profile_length as i32,
        cfg.cutoff_s,
        cfg.expect_value,
        false,
        0.0,
    );
    let mut saved = BlastKappaSavedParameters::s_saved_parameters_new(
        profile_length as i32,
        1,
        cfg.compo_adjust_mode,
        true,
    );
    let mut results = crate::hspstream::HspResults::new(1);

    let status = blast_redo_alignment_core_mt_in_memory_subject(
        crate::program::PSI_BLAST,
        1,
        profile_consensus,
        &query_info,
        &mut kbp_gap,
        &mut matrix,
        &mut scoring,
        &align_params,
        &mut saved,
        &mut local_list,
        &mut results,
        BlastRedoInMemorySubject {
            subject_source: user_query,
            reward: 0,
            penalty: 0,
            genetic_code: &crate::util::STANDARD_GENETIC_CODE,
            smith_waterman: cfg.smith_waterman,
            expect_value: cfg.expect_value,
            hitlist_size: cfg.hitlist_size,
            inclusion_ethresh: cfg.inclusion_ethresh,
            link_context: None,
        },
    );
    if status != 0 {
        return status;
    }

    if let Some(hitlist) = results.hitlists.into_iter().next().flatten() {
        if let Some(mut redone) = hitlist
            .hsp_lists
            .into_iter()
            .find(|list| list.oid == hsp_list.oid)
        {
            redone.hsp_max = hsp_list.hsp_max;
            *hsp_list = redone;
        } else {
            hsp_list.hsps.clear();
            hsp_list.best_evalue = i32::MAX as f64;
        }
    } else {
        hsp_list.hsps.clear();
        hsp_list.best_evalue = i32::MAX as f64;
    }

    0
}

/// blast-rs: Callback-backed branch of `Blast_RedoAlignmentCore_MT`; not a complete direct NCBI C port.
///
/// Callback-backed branch of `Blast_RedoAlignmentCore_MT`.
///
/// The compatibility translation [`blast_redo_alignment_core_mt`] preserves the
/// no-subject/no-composition pass-through. This helper covers the faithful
/// external-driver branch once callers provide the C-style callback table in
/// [`BlastRedoAlignParams::callbacks`]: convert the incoming HSP list to
/// distinct alignments, run ordinary or Smith-Waterman callback redo, rebuild
/// the HSP list, evaluate/purge, normalize, and feed the per-query results.
/// Position-specific/PSSM ordinary and Smith-Waterman callback redo are routed
/// through the already-ported PSSM helpers.
#[allow(clippy::too_many_arguments)]
pub fn blast_redo_alignment_core_mt_with_callbacks(
    program: ProgramType,
    _num_threads: u32,
    query: &[u8],
    query_info: &crate::queryinfo::QueryInfo,
    kbp_gap: &mut [crate::stat::KarlinBlk],
    matrix: &mut [Vec<i32>],
    scoring: &mut crate::parameters::ScoringParameters,
    align_params: &BlastRedoAlignParams,
    saved: &mut BlastKappaSavedParameters,
    this_match: &mut HspList,
    subject: BlastRedoCallbackSubjectConfig<'_>,
    results: &mut crate::hspstream::HspResults,
) -> i32 {
    let query_length = query_info.max_length as i32;
    let record_status = record_initial_search(
        saved,
        kbp_gap,
        matrix,
        scoring,
        query_length,
        align_params.compo_adjust_mode,
        align_params.position_based,
    );
    if record_status != 0 {
        return record_status;
    }

    rescale_search(
        kbp_gap,
        scoring,
        query_info.contexts.len() as i32,
        align_params.local_scaling_factor,
    );

    let status = blast_redo_alignment_core_one_match_with_callbacks_inner(
        program,
        query,
        query_info,
        kbp_gap,
        align_params,
        this_match,
        subject,
        results,
    );

    restore_search(
        saved,
        kbp_gap,
        matrix,
        scoring,
        query_length,
        align_params.position_based,
        align_params.compo_adjust_mode,
    );

    status
}

/// blast-rs: Materialized-subject branch of `Blast_RedoAlignmentCore_MT`; not a complete direct NCBI C port.
///
/// Materialized-subject branch of `Blast_RedoAlignmentCore_MT`.
///
/// This covers the single-match materialized-subject path for nucleotide,
/// protein, and translated-subject protein-space searches. It
/// follows the C loop shape for that use case:
/// `s_ResultHspToDistinctAlign`, `Blast_RedoOneMatch`,
/// `s_HSPListFromDistinctAlignments`, contained-HSP pruning,
/// `s_HitlistEvaluateAndPurge`, `s_HSPListNormalizeScores`, and
/// `s_ComputeNumIdentities`.
///
/// When `subject.smith_waterman` is true, this uses the in-memory nucleotide
/// Smith-Waterman redo branch for BLASTN and the protein-space Smith-Waterman
/// redo branch for BLASTP/BLASTX/TBLASTN/TBLASTX; otherwise it uses ordinary
/// gapped traceback with nucleotide reward/penalty scoring for BLASTN and
/// matrix scoring for protein-space searches, including materialized
/// translated-query buffers for BLASTX/TBLASTX.
/// Sum-stat/link-HSP evaluation is supported when
/// [`BlastRedoInMemorySubject::link_context`] is supplied.
/// Non-position protein-space, translated-subject protein-space, and
/// position-specific/PSSM composition adjustment are supported for ordinary
/// redo. Smith-Waterman redo covers no-composition materialized subjects,
/// non-position protein-space composition adjustment, translated-subject
/// composition adjustment, and position-specific/PSSM protein-space
/// adjustment. BLASTN/nucleotide composition adjustment is rejected before
/// Smith-Waterman redo because the adjusted redo branches operate in
/// protein/PSSM score space.
#[allow(clippy::too_many_arguments)]
pub fn blast_redo_alignment_core_mt_in_memory_subject(
    program: ProgramType,
    _num_threads: u32,
    query: &[u8],
    query_info: &crate::queryinfo::QueryInfo,
    kbp_gap: &mut [crate::stat::KarlinBlk],
    matrix: &mut [Vec<i32>],
    scoring: &mut crate::parameters::ScoringParameters,
    align_params: &BlastRedoAlignParams,
    saved: &mut BlastKappaSavedParameters,
    this_match: &mut HspList,
    results: &mut crate::hspstream::HspResults,
    subject: BlastRedoInMemorySubject<'_>,
) -> i32 {
    let query_length = query_info.max_length as i32;
    let record_status = record_initial_search(
        saved,
        kbp_gap,
        matrix,
        scoring,
        query_length,
        align_params.compo_adjust_mode,
        align_params.position_based,
    );
    if record_status != 0 {
        return record_status;
    }

    rescale_search(
        kbp_gap,
        scoring,
        query_info.contexts.len() as i32,
        align_params.local_scaling_factor,
    );

    let status = blast_redo_alignment_core_one_match_in_memory_inner(
        program,
        query,
        query_info,
        kbp_gap,
        align_params,
        this_match,
        results,
        None,
        None,
        subject,
    );

    restore_search(
        saved,
        kbp_gap,
        matrix,
        scoring,
        query_length,
        align_params.position_based,
        align_params.compo_adjust_mode,
    );

    status
}

/// blast-rs: Materialized-subject list branch of `Blast_RedoAlignmentCore_MT`; not a complete direct NCBI C port.
///
/// Materialized-subject list branch of `Blast_RedoAlignmentCore_MT`.
///
/// This models the `hsp_stream != NULL` control shape once the caller has
/// already fetched every subject sequence into memory: redo each match through
/// the single-subject helper, keep survivors in per-query `BlastCompo_Heap`
/// structures, apply `BlastCompo_EarlyTermination`, then collapse heaps via
/// [`fill_results_from_compo_heaps`]. The materialized branch supports
/// BLASTN, BLASTP, BLASTX, TBLASTN, TBLASTX, and ordinary PSI/PSSM redo when
/// translated-query callers supply the materialized flat translation buffer and
/// matching translated-query [`crate::queryinfo::QueryInfo`].
/// Sum-stat/link-HSP evaluation is supported per match when the subject entry
/// carries a [`HitlistLinkContext`].
#[allow(clippy::too_many_arguments)]
pub fn blast_redo_alignment_core_mt_in_memory_subjects(
    program: ProgramType,
    num_threads: u32,
    query: &[u8],
    query_info: &crate::queryinfo::QueryInfo,
    kbp_gap: &mut [crate::stat::KarlinBlk],
    matrix: &mut [Vec<i32>],
    scoring: &mut crate::parameters::ScoringParameters,
    align_params: &BlastRedoAlignParams,
    saved: &mut BlastKappaSavedParameters,
    matches: &mut [BlastRedoInMemorySubjectMatch<'_>],
    results: &mut crate::hspstream::HspResults,
) -> i32 {
    let query_length = query_info.max_length as i32;
    let record_status = record_initial_search(
        saved,
        kbp_gap,
        matrix,
        scoring,
        query_length,
        align_params.compo_adjust_mode,
        align_params.position_based,
    );
    if record_status != 0 {
        return record_status;
    }

    rescale_search(
        kbp_gap,
        scoring,
        query_info.contexts.len() as i32,
        align_params.local_scaling_factor,
    );

    let status = if matches.is_empty() {
        *results = crate::hspstream::HspResults::new(query_info.num_queries as i32);
        0
    } else if !materialized_redo_program_is_supported(program) {
        -1
    } else {
        let num_queries = query_info.num_queries.max(0) as usize;
        let heap_cfg = matches[0].subject;
        let num_workers = (num_threads.max(1) as usize).min(matches.len().max(1));
        let mut worker_heaps =
            vec![
                vec![
                    BlastCompoHeap::new(heap_cfg.hitlist_size, heap_cfg.inclusion_ethresh);
                    num_queries
                ];
                num_workers
            ];
        let mut status = 0;

        for (item_index, item) in matches.iter_mut().enumerate() {
            let worker_index = item_index % num_workers;
            let heaps = &mut worker_heaps[worker_index];
            if blast_compo_early_termination(item.hsp_list.best_evalue, heaps, num_queries) {
                continue;
            }

            let mut local_results = crate::hspstream::HspResults::new(num_queries as i32);
            status = blast_redo_alignment_core_one_match_in_memory_inner(
                program,
                query,
                query_info,
                kbp_gap,
                align_params,
                &mut item.hsp_list,
                &mut local_results,
                Some(&heaps[..]),
                None,
                item.subject,
            );
            if status != 0 {
                break;
            }

            for (query_index, hitlist) in local_results.hitlists.into_iter().enumerate() {
                if query_index >= heaps.len() {
                    continue;
                }
                if let Some(hitlist) = hitlist {
                    for hsp_list in hitlist.hsp_lists {
                        let _discarded = heaps[query_index].insert(hsp_list);
                    }
                }
            }
        }

        if status == 0 {
            let mut heaps =
                vec![
                    BlastCompoHeap::new(heap_cfg.hitlist_size, heap_cfg.inclusion_ethresh);
                    num_queries
                ];
            merge_compo_thread_heaps(&mut heaps, &mut worker_heaps);
            *results = fill_results_from_compo_heaps(&mut heaps);
        }
        status
    };

    restore_search(
        saved,
        kbp_gap,
        matrix,
        scoring,
        query_length,
        align_params.position_based,
        align_params.compo_adjust_mode,
    );

    status
}

/// blast-rs: `BlastSeqSrc` fetching branch of `Blast_RedoAlignmentCore_MT`; not a complete direct NCBI C port.
///
/// `BlastSeqSrc` fetching branch of `Blast_RedoAlignmentCore_MT`.
///
/// This mirrors the C `hsp_stream != NULL && subjectBlk == NULL` shape for
/// redo: each incoming subject HSP list names its database OID, the subject
/// sequence is fetched from the sequence source with the program's native
/// subject encoding, then the already-ported materialized one-match
/// redo/evaluate path feeds the per-query composition heaps. Non-position
/// protein-space and position-based/PSSM composition adjustment are handled by
/// the shared materialized one-match redo path. Program/mode combinations that
/// NCBI does not route through protein-space adjustment are rejected at that
/// boundary. Sum-stat/link-HSP evaluation is passed through when the caller
/// supplies [`BlastRedoSeqSrcSubjectConfig::link_context`].
#[allow(clippy::too_many_arguments)]
pub fn blast_redo_alignment_core_mt_seqsrc_subjects(
    program: ProgramType,
    num_threads: u32,
    query: &[u8],
    query_info: &crate::queryinfo::QueryInfo,
    kbp_gap: &mut [crate::stat::KarlinBlk],
    matrix: &mut [Vec<i32>],
    scoring: &mut crate::parameters::ScoringParameters,
    align_params: &BlastRedoAlignParams,
    saved: &mut BlastKappaSavedParameters,
    seqsrc: &dyn crate::seqsrc::BlastSeqSource,
    matches: &mut [HspList],
    subject_cfg: BlastRedoSeqSrcSubjectConfig<'_>,
    results: &mut crate::hspstream::HspResults,
) -> i32 {
    let query_length = query_info.max_length as i32;
    let record_status = record_initial_search(
        saved,
        kbp_gap,
        matrix,
        scoring,
        query_length,
        align_params.compo_adjust_mode,
        align_params.position_based,
    );
    if record_status != 0 {
        return record_status;
    }

    rescale_search(
        kbp_gap,
        scoring,
        query_info.contexts.len() as i32,
        align_params.local_scaling_factor,
    );

    let status = if matches.is_empty() {
        *results = crate::hspstream::HspResults::new(query_info.num_queries as i32);
        0
    } else if !materialized_redo_program_is_supported(program) {
        -1
    } else {
        let num_queries = query_info.num_queries.max(0) as usize;
        let num_workers = (num_threads.max(1) as usize).min(matches.len().max(1));
        let mut worker_heaps =
            vec![
                vec![
                    BlastCompoHeap::new(subject_cfg.hitlist_size, subject_cfg.inclusion_ethresh);
                    num_queries
                ];
                num_workers
            ];
        let encoding = seqsrc_encoding_for_redo_program(program);
        let mut status = 0;

        for (item_index, hsp_list) in matches.iter_mut().enumerate() {
            let worker_index = item_index % num_workers;
            let heaps = &mut worker_heaps[worker_index];
            if blast_compo_early_termination(hsp_list.best_evalue, heaps, num_queries) {
                continue;
            }

            if program == crate::program::TBLASTN {
                let mut matching_seq = BlastCompoMatchingSequence::default();
                if s_matching_sequence_initialize(
                    &mut matching_seq,
                    program,
                    seqsrc,
                    1,
                    hsp_list.oid,
                    None,
                ) != 0
                {
                    status = -1;
                    break;
                }
                let subject_sequence = {
                    let Some(local_data) = matching_seq.local_data.as_ref() else {
                        matching_sequence_release(&mut matching_seq);
                        status = -1;
                        break;
                    };
                    let Some(subject_blk) = local_data.seq_arg.seq.as_ref() else {
                        matching_sequence_release(&mut matching_seq);
                        status = -1;
                        break;
                    };
                    let Some(sequence) = subject_blk.sequence.as_deref() else {
                        matching_sequence_release(&mut matching_seq);
                        status = -1;
                        break;
                    };
                    if subject_blk.length < 0 || subject_blk.length as usize > sequence.len() {
                        matching_sequence_release(&mut matching_seq);
                        status = -1;
                        break;
                    }
                    &sequence[..subject_blk.length as usize]
                };
                let mut local_results = crate::hspstream::HspResults::new(num_queries as i32);
                status = blast_redo_alignment_core_one_match_in_memory_inner(
                    program,
                    query,
                    query_info,
                    kbp_gap,
                    align_params,
                    hsp_list,
                    &mut local_results,
                    Some(&heaps[..]),
                    None,
                    BlastRedoInMemorySubject {
                        subject_source: subject_sequence,
                        reward: subject_cfg.reward,
                        penalty: subject_cfg.penalty,
                        genetic_code: subject_cfg.genetic_code,
                        smith_waterman: subject_cfg.smith_waterman,
                        expect_value: subject_cfg.expect_value,
                        hitlist_size: subject_cfg.hitlist_size,
                        inclusion_ethresh: subject_cfg.inclusion_ethresh,
                        link_context: subject_cfg.link_context,
                    },
                );
                matching_sequence_release(&mut matching_seq);
                if status != 0 {
                    break;
                }

                for (query_index, hitlist) in local_results.hitlists.into_iter().enumerate() {
                    if query_index >= heaps.len() {
                        continue;
                    }
                    if let Some(hitlist) = hitlist {
                        for hsp_list in hitlist.hsp_lists {
                            let _discarded = heaps[query_index].insert(hsp_list);
                        }
                    }
                }
                continue;
            }

            let Some(seq_data) = seqsrc.get_sequence(&crate::seqsrc::GetSeqArg {
                oid: hsp_list.oid,
                encoding,
                ..crate::seqsrc::GetSeqArg::default()
            }) else {
                status = -1;
                break;
            };
            if seq_data.length < 0 || seq_data.length as usize > seq_data.sequence.len() {
                status = -1;
                break;
            }
            let subject_sequence = &seq_data.sequence[..seq_data.length as usize];
            let mut local_results = crate::hspstream::HspResults::new(num_queries as i32);
            status = blast_redo_alignment_core_one_match_in_memory_inner(
                program,
                query,
                query_info,
                kbp_gap,
                align_params,
                hsp_list,
                &mut local_results,
                Some(&heaps[..]),
                None,
                BlastRedoInMemorySubject {
                    subject_source: subject_sequence,
                    reward: subject_cfg.reward,
                    penalty: subject_cfg.penalty,
                    genetic_code: subject_cfg.genetic_code,
                    smith_waterman: subject_cfg.smith_waterman,
                    expect_value: subject_cfg.expect_value,
                    hitlist_size: subject_cfg.hitlist_size,
                    inclusion_ethresh: subject_cfg.inclusion_ethresh,
                    link_context: subject_cfg.link_context,
                },
            );
            if status != 0 {
                break;
            }

            for (query_index, hitlist) in local_results.hitlists.into_iter().enumerate() {
                if query_index >= heaps.len() {
                    continue;
                }
                if let Some(hitlist) = hitlist {
                    for hsp_list in hitlist.hsp_lists {
                        let _discarded = heaps[query_index].insert(hsp_list);
                    }
                }
            }
        }

        if status == 0 {
            let mut heaps =
                vec![
                    BlastCompoHeap::new(subject_cfg.hitlist_size, subject_cfg.inclusion_ethresh);
                    num_queries
                ];
            merge_compo_thread_heaps(&mut heaps, &mut worker_heaps);
            *results = fill_results_from_compo_heaps(&mut heaps);
        }
        status
    };

    restore_search(
        saved,
        kbp_gap,
        matrix,
        scoring,
        query_length,
        align_params.position_based,
        align_params.compo_adjust_mode,
    );

    status
}

/// NCBI-shaped materialized-subject branch when the caller still has a
/// `BlastHSPList`.
///
/// C `Blast_RedoAlignmentCore_MT` carries `BlastHSPList.query_index` into the
/// redo/update phase. This preserves that owner query while the remaining redo
/// internals still operate on legacy `HspList` alignment nodes.
#[allow(clippy::too_many_arguments)]
pub fn blast_redo_alignment_core_mt_in_memory_subject_hsp_list(
    program: ProgramType,
    _num_threads: u32,
    query: &[u8],
    query_info: &crate::queryinfo::QueryInfo,
    kbp_gap: &mut [crate::stat::KarlinBlk],
    matrix: &mut [Vec<i32>],
    scoring: &mut crate::parameters::ScoringParameters,
    align_params: &BlastRedoAlignParams,
    saved: &mut BlastKappaSavedParameters,
    this_match: &mut crate::hspstream::BlastHSPList,
    results: &mut crate::hspstream::HspResults,
    subject: BlastRedoInMemorySubject<'_>,
) -> i32 {
    let source_query_index = this_match.query_index;
    let placeholder = crate::hspstream::BlastHSPList {
        oid: -1,
        query_index: source_query_index,
        hsp_array: Vec::new(),
        hspcnt: 0,
        allocated: 0,
        hsp_max: 0,
        do_not_reallocate: false,
        best_evalue: i32::MAX as f64,
    };
    let mut legacy_match = std::mem::replace(this_match, placeholder).into_legacy_hsp_list();
    let query_length = query_info.max_length as i32;
    let record_status = record_initial_search(
        saved,
        kbp_gap,
        matrix,
        scoring,
        query_length,
        align_params.compo_adjust_mode,
        align_params.position_based,
    );
    if record_status != 0 {
        *this_match = crate::hspstream::BlastHSPList::from_legacy_hsp_list(legacy_match);
        this_match.query_index = source_query_index;
        return record_status;
    }

    rescale_search(
        kbp_gap,
        scoring,
        query_info.contexts.len() as i32,
        align_params.local_scaling_factor,
    );

    let status = blast_redo_alignment_core_one_match_in_memory_inner(
        program,
        query,
        query_info,
        kbp_gap,
        align_params,
        &mut legacy_match,
        results,
        None,
        Some(source_query_index),
        subject,
    );

    restore_search(
        saved,
        kbp_gap,
        matrix,
        scoring,
        query_length,
        align_params.position_based,
        align_params.compo_adjust_mode,
    );

    *this_match = crate::hspstream::BlastHSPList::from_legacy_hsp_list(legacy_match);
    this_match.query_index = source_query_index;
    status
}

#[derive(Clone, Copy)]
struct BlastRedoOneMatchCallbackPlan {
    num_contexts: usize,
    num_frames: usize,
    query_index: i32,
    context_index: usize,
}

/// NCBI: argument and context validation phase of `Blast_RedoAlignmentCore_MT`.
fn blast_redo_alignment_core_one_match_with_callbacks_validate(
    program: ProgramType,
    query_info: &crate::queryinfo::QueryInfo,
    align_params: &BlastRedoAlignParams,
    this_match: &HspList,
) -> Result<BlastRedoOneMatchCallbackPlan, i32> {
    if align_params.callbacks.is_none() {
        return Err(-1);
    }
    if this_match.hsps.is_empty() {
        return Err(0);
    }
    if !materialized_redo_program_is_supported(program) {
        return Err(-1);
    }
    if !matches!(
        align_params.compo_adjust_mode,
        CompoAdjustMode::NoCompositionBasedStats
    ) && program == crate::program::BLASTN
    {
        return Err(-1);
    }
    if !matches!(
        align_params.compo_adjust_mode,
        CompoAdjustMode::NoCompositionBasedStats
    ) && !composition_redo_program_is_supported(program)
    {
        return Err(-1);
    }
    if align_params.position_based
        && (align_params.query_is_translated
            || align_params.subject_is_translated
            || program == crate::program::BLASTN)
    {
        return Err(-1);
    }

    let num_contexts = query_info.contexts.len();
    if num_contexts == 0 {
        return Err(-1);
    }
    let num_frames = materialized_num_contexts_per_query(program, query_info);
    let first_context = this_match.hsps[0].context.max(0) as usize;
    let query_index = query_info
        .contexts
        .get(first_context)
        .map(|ctx| ctx.query_index)
        .unwrap_or_else(|| {
            crate::queryinfo::blast_get_query_index_from_context(
                this_match.hsps[0].context,
                program,
            )
        });
    if query_index < 0 {
        return Err(-1);
    }
    let context_index = query_index as usize * num_frames;
    if context_index >= num_contexts {
        return Err(-1);
    }

    Ok(BlastRedoOneMatchCallbackPlan {
        num_contexts,
        num_frames,
        query_index,
        context_index,
    })
}

/// NCBI: `ResultHspToDistinctAlign` phase of `Blast_RedoAlignmentCore_MT`.
fn blast_redo_alignment_core_one_match_with_callbacks_distinct_alignments(
    this_match: &HspList,
    context_index: usize,
    local_scaling_factor: f64,
) -> Result<([Option<Box<BlastCompoAlignment>>; 6], [i32; 6]), i32> {
    let mut incoming_align_set: [Option<Box<BlastCompoAlignment>>; 6] = Default::default();
    let mut num_aligns = [0i32; 6];
    let rc = result_hsp_to_distinct_align(
        &mut incoming_align_set,
        &mut num_aligns,
        &this_match.hsps,
        context_index as i32,
        local_scaling_factor,
    );
    if rc != 0 {
        return Err(rc);
    }
    Ok((incoming_align_set, num_aligns))
}

/// NCBI: per-frame composition redo plus `HspListFromDistinctAlignments`.
#[allow(clippy::too_many_arguments)]
fn blast_redo_alignment_core_one_match_with_callbacks_redo_frame(
    program: ProgramType,
    query_info: &crate::queryinfo::QueryInfo,
    kbp_gap: &[crate::stat::KarlinBlk],
    align_params: &BlastRedoAlignParams,
    subject: BlastRedoCallbackSubjectConfig<'_>,
    matching_seq: &BlastCompoMatchingSequence,
    compo_query_info: &[BlastCompoQueryInfo],
    incoming_aligns: &Option<Box<BlastCompoAlignment>>,
    num_aligns: i32,
    frame_index: usize,
    context_index: usize,
    num_contexts: usize,
    alignments: &mut [Option<Box<BlastCompoAlignment>>],
    redone: &mut HspList,
    pvalue_for_this_pair: &mut f64,
    lambda_ratio: &mut f64,
) -> i32 {
    if incoming_aligns.is_none() {
        return 0;
    }
    let context = context_index + frame_index;
    if context >= num_contexts || context >= kbp_gap.len() {
        return -1;
    }

    let status = if subject.smith_waterman {
        let significant_matches =
            vec![
                BlastCompoHeap::new(subject.hitlist_size, subject.inclusion_ethresh);
                compo_query_info.len()
            ];
        *pvalue_for_this_pair = -1.0;
        *lambda_ratio = 1.0;
        if align_params.position_based {
            let mut adjusted_pssm = Vec::new();
            blast_redo_one_match_smith_waterman_with_callbacks_and_position_adjustment(
                alignments,
                align_params,
                incoming_aligns,
                num_aligns,
                kbp_gap[context].lambda,
                kbp_gap[context].log_k,
                matching_seq,
                compo_query_info,
                &significant_matches,
                &mut adjusted_pssm,
                pvalue_for_this_pair,
                composition_test_index_for_redo(align_params),
                lambda_ratio,
            )
        } else {
            let mut adjusted_matrix = match square_matrix_from_vec(&align_params.matrix_info.matrix)
            {
                Some(matrix) => matrix,
                None => [[0i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE],
            };
            let workspace = if matches!(
                align_params.compo_adjust_mode,
                CompoAdjustMode::CompositionMatrixAdjust
            ) {
                Some(BlastCompositionWorkspace::blosum62())
            } else {
                None
            };
            blast_redo_one_match_smith_waterman_with_callbacks_and_adjustment(
                alignments,
                align_params,
                incoming_aligns,
                num_aligns,
                kbp_gap[context].lambda,
                kbp_gap[context].log_k,
                matching_seq,
                compo_query_info,
                &significant_matches,
                &mut adjusted_matrix,
                workspace.as_ref(),
                pvalue_for_this_pair,
                composition_test_index_for_redo(align_params),
                lambda_ratio,
            )
        }
    } else {
        *pvalue_for_this_pair = -1.0;
        *lambda_ratio = 1.0;
        if align_params.position_based {
            let mut adjusted_pssm = Vec::new();
            blast_redo_one_match_with_callbacks_and_position_adjustment(
                alignments,
                align_params,
                incoming_aligns,
                num_aligns,
                kbp_gap[context].lambda,
                matching_seq,
                compo_query_info,
                &mut adjusted_pssm,
                pvalue_for_this_pair,
                composition_test_index_for_redo(align_params),
                lambda_ratio,
            )
        } else {
            let mut adjusted_matrix = match square_matrix_from_vec(&align_params.matrix_info.matrix)
            {
                Some(matrix) => matrix,
                None => [[0i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE],
            };
            let workspace = if matches!(
                align_params.compo_adjust_mode,
                CompoAdjustMode::CompositionMatrixAdjust
            ) {
                Some(BlastCompositionWorkspace::blosum62())
            } else {
                None
            };
            blast_redo_one_match_with_callbacks_and_adjustment(
                alignments,
                align_params,
                incoming_aligns,
                num_aligns,
                kbp_gap[context].lambda,
                matching_seq,
                compo_query_info,
                &mut adjusted_matrix,
                workspace.as_ref(),
                pvalue_for_this_pair,
                composition_test_index_for_redo(align_params),
                lambda_ratio,
            )
        }
    };
    if status != 0 {
        return status;
    }

    if alignments[context].is_some() {
        let qframe = if crate::program::blast_query_is_translated(program) {
            let f = frame_index as i32;
            if f < 3 {
                f + 1
            } else {
                2 - f
            }
        } else {
            query_info.contexts[context].frame
        };
        let (status, _tags) = hsp_list_from_distinct_alignments(
            redone,
            &mut alignments[context],
            matching_seq.index,
            qframe,
        );
        if status != 0 {
            return status;
        }
    }

    0
}

fn composition_test_index_for_redo(params: &BlastRedoAlignParams) -> i32 {
    let _ = params;
    // C passes `extendParams->options->unifiedP` into
    // `Blast_RedoOneMatch`. `BlastExtensionOptionsNew` zero-initializes that
    // field, and this Rust parameter struct does not yet carry an explicit
    // `unifiedP` option. Do not derive it from `compositionBasedStats`; CBS2
    // (`compositionBasedStats == 2`) and the unified-p-value switch are
    // independent fields.
    0
}

/// NCBI: `s_HitlistEvaluateAndPurge` phase for callback-backed redo.
#[allow(clippy::too_many_arguments)]
fn blast_redo_alignment_core_one_match_with_callbacks_evaluate(
    program: ProgramType,
    query_info: &crate::queryinfo::QueryInfo,
    kbp_gap: &[crate::stat::KarlinBlk],
    align_params: &BlastRedoAlignParams,
    subject: BlastRedoCallbackSubjectConfig<'_>,
    redone: &mut HspList,
    context_index: usize,
    pvalue_for_this_pair: f64,
) -> (usize, f64) {
    if redone.hsps.len() > 1 {
        s_hitlist_reap_contained(&mut redone.hsps);
    }
    let eval_context = redone
        .hsps
        .first()
        .and_then(|hsp| {
            let context = hsp.context.max(0) as usize;
            (context < query_info.contexts.len()).then_some(context)
        })
        .unwrap_or_else(|| context_index.min(query_info.contexts.len() - 1));
    let ctx = &query_info.contexts[eval_context];
    let (_best_score, best_evalue) = s_hitlist_evaluate_and_purge(
        redone,
        subject.subject_length,
        program,
        ctx.query_length,
        ctx.length_adjustment,
        ctx.eff_searchsp as f64,
        kbp_gap.get(eval_context),
        subject.link_context,
        pvalue_for_this_pair,
        subject.expect_value,
        align_params.do_link_hsps,
    );
    (eval_context, best_evalue)
}

/// NCBI: accepted/rejected hit-list update phase of `Blast_RedoAlignmentCore_MT`.
#[allow(clippy::too_many_arguments)]
fn blast_redo_alignment_core_one_match_with_callbacks_update_results(
    query_index: i32,
    eval_context: usize,
    best_evalue: f64,
    query_info: &[BlastCompoQueryInfo],
    kbp_gap: &[crate::stat::KarlinBlk],
    align_params: &BlastRedoAlignParams,
    subject: BlastRedoCallbackSubjectConfig<'_>,
    matching_seq: &BlastCompoMatchingSequence,
    redone: &mut HspList,
    this_match: &mut HspList,
    results: &mut crate::hspstream::HspResults,
) -> i32 {
    if best_evalue <= subject.expect_value {
        if let Some(kbp) = kbp_gap.get(eval_context) {
            s_hsp_list_normalize_scores(
                redone,
                kbp.lambda,
                kbp.log_k,
                align_params.local_scaling_factor,
            );
        }
        if compute_num_identities_with_callbacks(
            redone,
            align_params,
            matching_seq,
            query_info,
            subject.subject_length,
        ) != 0
        {
            return -1;
        }
        *this_match = redone.clone();
        let idx = query_index as usize;
        if idx >= results.hitlists.len() {
            results.hitlists.resize_with(idx + 1, || None);
        }
        let hitlist = results.hitlists[idx].get_or_insert_with(crate::hspstream::HitList::new);
        if hitlist.hsp_lists.len() < subject.hitlist_size.max(0) as usize
            || best_evalue <= subject.inclusion_ethresh
        {
            hitlist.blast_hit_list_update(redone.clone());
            hitlist.sort_by_evalue();
        }
    } else {
        this_match.hsps.clear();
        // NCBI `s_BlastGetBestEvalue` seeds with `(double)INT4_MAX`.
        this_match.best_evalue = i32::MAX as f64;
    }

    0
}

/// blast-rs: Shared callback-backed one-match redo/evaluate path; not a direct NCBI C port.
#[allow(clippy::too_many_arguments)]
fn blast_redo_alignment_core_one_match_with_callbacks_inner(
    program: ProgramType,
    query: &[u8],
    query_info: &crate::queryinfo::QueryInfo,
    kbp_gap: &[crate::stat::KarlinBlk],
    align_params: &BlastRedoAlignParams,
    this_match: &mut HspList,
    subject: BlastRedoCallbackSubjectConfig<'_>,
    results: &mut crate::hspstream::HspResults,
) -> i32 {
    let plan = match blast_redo_alignment_core_one_match_with_callbacks_validate(
        program,
        query_info,
        align_params,
        this_match,
    ) {
        Ok(plan) => plan,
        Err(status) => return status,
    };
    let (incoming_align_set, num_aligns) =
        match blast_redo_alignment_core_one_match_with_callbacks_distinct_alignments(
            this_match,
            plan.context_index,
            align_params.local_scaling_factor,
        ) {
            Ok(converted) => converted,
            Err(status) => return status,
        };

    let compo_query_info = get_query_info(query, query_info, program == crate::program::BLASTX);
    let matching_seq = matching_sequence_initialize(subject.subject_length, this_match.oid, 0);
    let mut alignments = vec![None; plan.num_contexts];
    let mut redone = HspList::new(this_match.oid);
    let mut pvalue_for_this_pair = -1.0;
    let mut lambda_ratio = 1.0;

    for frame_index in 0..plan.num_frames {
        let status = blast_redo_alignment_core_one_match_with_callbacks_redo_frame(
            program,
            query_info,
            kbp_gap,
            align_params,
            subject,
            &matching_seq,
            &compo_query_info,
            &incoming_align_set[frame_index],
            num_aligns[frame_index],
            frame_index,
            plan.context_index,
            plan.num_contexts,
            &mut alignments,
            &mut redone,
            &mut pvalue_for_this_pair,
            &mut lambda_ratio,
        );
        if status != 0 {
            return status;
        }
    }

    let (eval_context, best_evalue) = blast_redo_alignment_core_one_match_with_callbacks_evaluate(
        program,
        query_info,
        kbp_gap,
        align_params,
        subject,
        &mut redone,
        plan.context_index,
        pvalue_for_this_pair,
    );

    blast_redo_alignment_core_one_match_with_callbacks_update_results(
        plan.query_index,
        eval_context,
        best_evalue,
        &compo_query_info,
        kbp_gap,
        align_params,
        subject,
        &matching_seq,
        &mut redone,
        this_match,
        results,
    )
}

/// blast-rs: Callback-backed identity recomputation helper; not a direct NCBI C port.
fn compute_num_identities_with_callbacks(
    hsp_list: &mut HspList,
    params: &BlastRedoAlignParams,
    matching_seq: &BlastCompoMatchingSequence,
    query_info: &[BlastCompoQueryInfo],
    subject_length: i32,
) -> i32 {
    let Some(callbacks) = params.callbacks else {
        return -1;
    };
    let Some(get_range) = callbacks.get_range else {
        return -1;
    };

    for hsp in &mut hsp_list.hsps {
        let context = hsp.context.max(0) as usize;
        let Some(query_context) = query_info.get(context) else {
            return -1;
        };
        let query_range = BlastCompoSequenceRange {
            begin: query_context.origin,
            end: query_context.origin + query_context.seq.length,
            context: hsp.context,
        };
        let subject_range = BlastCompoSequenceRange {
            begin: 0,
            end: subject_length,
            context: hsp.subject_frame,
        };
        let align = BlastCompoAlignment::new(
            hsp.score,
            MatrixAdjustRule::DontAdjust,
            hsp.context,
            hsp.query_offset,
            hsp.query_end,
            hsp.subject_offset,
            hsp.subject_end,
            hsp.subject_frame,
            hsp.edit_script
                .clone()
                .map(|edit_script| BlastCompoAlignmentContext {
                    edit_script: Some(edit_script),
                    hsp: None,
                }),
        )
        .with_hsp_context(hsp.query_gapped_start, hsp.subject_gapped_start);
        let mut subject = BlastCompoSequenceData::default();
        let mut query = BlastCompoSequenceData::default();
        let mut subject_maybe_biased = false;
        let status = get_range(
            matching_seq,
            &subject_range,
            &mut subject,
            &query_context.seq,
            &query_range,
            &mut query,
            &query_context.words,
            &align,
            false,
            params.compo_adjust_mode,
            params.local_scaling_factor != 1.0,
            &mut subject_maybe_biased,
        );
        if status != 0 {
            sequence_data_release(&mut subject);
            sequence_data_release(&mut query);
            return status;
        }

        let mut local_hsp = hsp.clone();
        local_hsp.query_offset -= query_range.begin;
        local_hsp.query_end -= query_range.begin;
        local_hsp.query_gapped_start -= query_range.begin;
        local_hsp.subject_offset -= subject_range.begin;
        local_hsp.subject_end -= subject_range.begin;
        local_hsp.subject_gapped_start -= subject_range.begin;
        let edit_ops = hsp.edit_script.as_ref().map(|script| script.ops_vec());
        let (num_ident, _align_length, _num_pos) = blast_hsp_get_num_identities(
            query.data(),
            subject.data(),
            &local_hsp,
            edit_ops.as_deref(),
            None,
        );
        hsp.num_ident = num_ident;
        sequence_data_release(&mut subject);
        sequence_data_release(&mut query);
    }

    0
}

/// blast-rs: Shared materialized one-match redo/evaluate path; not a direct NCBI C port.
#[allow(clippy::too_many_arguments)]
fn blast_redo_alignment_core_one_match_in_memory_inner(
    program: ProgramType,
    query: &[u8],
    query_info: &crate::queryinfo::QueryInfo,
    kbp_gap: &[crate::stat::KarlinBlk],
    align_params: &BlastRedoAlignParams,
    this_match: &mut HspList,
    results: &mut crate::hspstream::HspResults,
    significant_matches_for_sw: Option<&[BlastCompoHeap]>,
    explicit_query_index: Option<i32>,
    subject: BlastRedoInMemorySubject<'_>,
) -> i32 {
    let no_composition = matches!(
        align_params.compo_adjust_mode,
        CompoAdjustMode::NoCompositionBasedStats
    );
    if !no_composition && program == crate::program::BLASTN {
        return -1;
    }
    if !no_composition && !composition_redo_program_is_supported(program) {
        return -1;
    }
    if align_params.position_based
        && (align_params.query_is_translated || program == crate::program::BLASTN)
    {
        return -1;
    }
    if align_params.do_link_hsps && subject.link_context.is_none() {
        return -1;
    }
    if this_match.hsps.is_empty() {
        return 0;
    }
    if !materialized_redo_program_is_supported(program) {
        return -1;
    }

    let num_contexts = query_info.contexts.len();
    if num_contexts == 0 {
        return -1;
    }
    let num_frames = materialized_num_contexts_per_query(program, query_info);
    let first_context = this_match.hsps[0].context.max(0) as usize;
    let query_index = explicit_query_index.unwrap_or_else(|| {
        query_info
            .contexts
            .get(first_context)
            .map(|ctx| ctx.query_index)
            .unwrap_or_else(|| {
                crate::queryinfo::blast_get_query_index_from_context(
                    this_match.hsps[0].context,
                    program,
                )
            })
    });
    if query_index < 0 {
        return -1;
    }
    let context_index = query_index as usize * num_frames;
    if context_index >= num_contexts {
        return -1;
    }

    let mut incoming_align_set: [Option<Box<BlastCompoAlignment>>; 6] = Default::default();
    let mut num_aligns = [0i32; 6];
    let rc = result_hsp_to_distinct_align(
        &mut incoming_align_set,
        &mut num_aligns,
        &this_match.hsps,
        context_index as i32,
        align_params.local_scaling_factor,
    );
    if rc != 0 {
        return rc;
    }

    let compo_query_info = get_query_info(query, query_info, program == crate::program::BLASTX);
    let matching_seq =
        matching_sequence_initialize(subject.subject_source.len() as i32, this_match.oid, 0);
    let mut alignments = vec![None; num_contexts];
    let mut redone = HspList::new(this_match.oid);
    let mut pvalue_for_this_pair = -1.0;
    let mut _lambda_ratio = 1.0;

    for frame_index in 0..num_frames {
        let incoming_aligns = &incoming_align_set[frame_index];
        if incoming_aligns.is_none() {
            continue;
        }
        let context = context_index + frame_index;
        if context >= num_contexts || context >= kbp_gap.len() {
            return -1;
        }
        let status = if subject.smith_waterman {
            let owned_significant_matches;
            let significant_matches = if let Some(matches) = significant_matches_for_sw {
                if matches.len() >= compo_query_info.len() {
                    matches
                } else {
                    owned_significant_matches =
                        vec![
                            BlastCompoHeap::new(subject.hitlist_size, subject.inclusion_ethresh);
                            compo_query_info.len()
                        ];
                    &owned_significant_matches
                }
            } else {
                owned_significant_matches =
                    vec![
                        BlastCompoHeap::new(subject.hitlist_size, subject.inclusion_ethresh);
                        compo_query_info.len()
                    ];
                &owned_significant_matches
            };
            if program == crate::program::BLASTN {
                let nucleotide_matrix =
                    crate::matrix::nucleotide_matrix(subject.reward, subject.penalty);
                blast_redo_one_match_smith_waterman_in_memory_nucl(
                    &mut alignments,
                    align_params,
                    incoming_aligns,
                    num_aligns[frame_index],
                    kbp_gap[context].lambda,
                    kbp_gap[context].log_k,
                    &matching_seq,
                    &compo_query_info,
                    subject.subject_source,
                    &significant_matches,
                    &nucleotide_matrix,
                )
            } else if align_params.position_based {
                let mut adjusted_pssm = Vec::new();
                pvalue_for_this_pair = -1.0;
                _lambda_ratio = 1.0;
                blast_redo_one_match_smith_waterman_in_memory_protein_position_adjustment(
                    &mut alignments,
                    align_params,
                    incoming_aligns,
                    num_aligns[frame_index],
                    kbp_gap[context].lambda,
                    kbp_gap[context].log_k,
                    &matching_seq,
                    &compo_query_info,
                    program,
                    subject.subject_source,
                    &significant_matches,
                    subject.genetic_code,
                    &mut adjusted_pssm,
                    &mut pvalue_for_this_pair,
                    composition_test_index_for_redo(align_params),
                    &mut _lambda_ratio,
                )
            } else if no_composition {
                blast_redo_one_match_smith_waterman_in_memory_protein(
                    &mut alignments,
                    align_params,
                    incoming_aligns,
                    num_aligns[frame_index],
                    kbp_gap[context].lambda,
                    kbp_gap[context].log_k,
                    &matching_seq,
                    &compo_query_info,
                    program,
                    subject.subject_source,
                    &significant_matches,
                    subject.genetic_code,
                )
            } else {
                let mut adjusted_matrix = square_matrix_from_vec(&align_params.matrix_info.matrix)
                    .unwrap_or([[0i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE]);
                let workspace = if matches!(
                    align_params.compo_adjust_mode,
                    CompoAdjustMode::CompositionMatrixAdjust
                ) {
                    Some(BlastCompositionWorkspace::blosum62())
                } else {
                    None
                };
                pvalue_for_this_pair = -1.0;
                _lambda_ratio = 1.0;
                blast_redo_one_match_smith_waterman_in_memory_protein_with_adjustment(
                    &mut alignments,
                    align_params,
                    incoming_aligns,
                    num_aligns[frame_index],
                    kbp_gap[context].lambda,
                    kbp_gap[context].log_k,
                    &matching_seq,
                    &compo_query_info,
                    program,
                    subject.subject_source,
                    &significant_matches,
                    subject.genetic_code,
                    &mut adjusted_matrix,
                    workspace.as_ref(),
                    &mut pvalue_for_this_pair,
                    composition_test_index_for_redo(align_params),
                    &mut _lambda_ratio,
                )
            }
        } else if no_composition {
            blast_redo_one_match_in_memory(
                &mut alignments,
                align_params,
                incoming_aligns,
                num_aligns[frame_index],
                kbp_gap[context].lambda,
                &matching_seq,
                &compo_query_info,
                program,
                subject.subject_source,
                subject.reward,
                subject.penalty,
                subject.genetic_code,
            )
        } else if align_params.position_based {
            let mut adjusted_pssm = Vec::new();
            pvalue_for_this_pair = -1.0;
            _lambda_ratio = 1.0;
            blast_redo_one_match_in_memory_with_position_adjustment(
                &mut alignments,
                align_params,
                incoming_aligns,
                num_aligns[frame_index],
                kbp_gap[context].lambda,
                &matching_seq,
                &compo_query_info,
                program,
                subject.subject_source,
                subject.genetic_code,
                &mut adjusted_pssm,
                &mut pvalue_for_this_pair,
                composition_test_index_for_redo(align_params),
                &mut _lambda_ratio,
            )
        } else {
            let mut adjusted_matrix = square_matrix_from_vec(&align_params.matrix_info.matrix)
                .unwrap_or([[0i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE]);
            let workspace = if matches!(
                align_params.compo_adjust_mode,
                CompoAdjustMode::CompositionMatrixAdjust
            ) {
                Some(BlastCompositionWorkspace::blosum62())
            } else {
                None
            };
            pvalue_for_this_pair = -1.0;
            _lambda_ratio = 1.0;
            blast_redo_one_match_in_memory_with_adjustment(
                &mut alignments,
                align_params,
                incoming_aligns,
                num_aligns[frame_index],
                kbp_gap[context].lambda,
                &matching_seq,
                &compo_query_info,
                program,
                subject.subject_source,
                subject.genetic_code,
                &mut adjusted_matrix,
                workspace.as_ref(),
                &mut pvalue_for_this_pair,
                composition_test_index_for_redo(align_params),
                &mut _lambda_ratio,
            )
        };
        if status != 0 {
            return status;
        }

        if alignments[context].is_some() {
            let qframe = if crate::program::blast_query_is_translated(program) {
                let f = frame_index as i32;
                if f < 3 {
                    f + 1
                } else {
                    2 - f
                }
            } else {
                query_info.contexts[context].frame
            };
            let (status, _tags) = hsp_list_from_distinct_alignments(
                &mut redone,
                &mut alignments[context],
                matching_seq.index,
                qframe,
            );
            if status != 0 {
                return status;
            }
        }
    }

    if redone.hsps.len() > 1 {
        s_hitlist_reap_contained(&mut redone.hsps);
    }
    let fallback_eval_context = redone
        .hsps
        .first()
        .and_then(|hsp| {
            let context = hsp.context.max(0) as usize;
            (context < query_info.contexts.len()).then_some(context)
        })
        .unwrap_or_else(|| context_index.min(query_info.contexts.len() - 1));
    let ctx = &query_info.contexts[fallback_eval_context];
    let kbp_for_eval = if align_params.do_link_hsps {
        kbp_gap.get(fallback_eval_context)
    } else {
        let mut best_evalue = i32::MAX as f64;
        for hsp in &mut redone.hsps {
            let hsp_context = hsp.context.max(0) as usize;
            let Some(hsp_ctx) = query_info
                .contexts
                .get(hsp_context)
                .or_else(|| query_info.contexts.get(fallback_eval_context))
            else {
                continue;
            };
            let Some(kbp) = kbp_gap
                .get(hsp_context)
                .or_else(|| kbp_gap.get(fallback_eval_context))
            else {
                continue;
            };
            hsp.evalue = kbp.raw_to_evalue(hsp.score, hsp_ctx.eff_searchsp as f64);
            if hsp.evalue < best_evalue {
                best_evalue = hsp.evalue;
            }
        }
        redone.best_evalue = best_evalue;
        None
    };
    let (best_score, best_evalue) = s_hitlist_evaluate_and_purge(
        &mut redone,
        subject.subject_source.len() as i32,
        program,
        ctx.query_length,
        ctx.length_adjustment,
        ctx.eff_searchsp as f64,
        kbp_for_eval,
        subject.link_context,
        pvalue_for_this_pair,
        subject.expect_value,
        align_params.do_link_hsps,
    );
    let _ = best_score;

    if best_evalue <= subject.expect_value {
        for hsp in &mut redone.hsps {
            let hsp_context = hsp.context.max(0) as usize;
            if let Some(kbp) = kbp_gap
                .get(hsp_context)
                .or_else(|| kbp_gap.get(fallback_eval_context))
            {
                s_hsp_normalize_score(
                    hsp,
                    kbp.lambda,
                    kbp.log_k,
                    align_params.local_scaling_factor,
                );
            }
        }
        let edit_ops: Vec<Option<Vec<(crate::gapinfo::GapAlignOpType, i32)>>> = redone
            .hsps
            .iter()
            .map(|hsp| hsp.edit_script.as_ref().map(|script| script.ops_vec()))
            .collect();
        if crate::program::blast_subject_is_translated(program) {
            if compute_num_identities_translated_subject_by_context(
                &compo_query_info,
                query,
                subject.subject_source,
                &mut redone,
                &edit_ops,
                None,
                subject.genetic_code,
            )
            .is_err()
            {
                return -1;
            }
        } else {
            compute_num_identities_protein_by_context(
                &compo_query_info,
                query,
                subject.subject_source,
                &mut redone,
                &edit_ops,
                None,
            );
        }
        *this_match = redone.clone();
        let idx = query_index as usize;
        if idx >= results.hitlists.len() {
            results.hitlists.resize_with(idx + 1, || None);
        }
        let hitlist = results.hitlists[idx].get_or_insert_with(crate::hspstream::HitList::new);
        if hitlist.hsp_lists.len() < subject.hitlist_size.max(0) as usize
            || best_evalue <= subject.inclusion_ethresh
        {
            hitlist.blast_hit_list_update(redone);
            hitlist.sort_by_evalue();
        }
    } else {
        this_match.hsps.clear();
        // NCBI `s_BlastGetBestEvalue` seeds with `(double)INT4_MAX`.
        this_match.best_evalue = i32::MAX as f64;
    }

    0
}

/// blast-rs: Converts Rust alignment lists to `s_IsContained` tuple inputs; not a direct NCBI C port.
fn alignment_list_to_containment_tuples(
    head: Option<&BlastCompoAlignment>,
) -> Vec<(i32, i32, i32, i32, i32, i32)> {
    let mut out = Vec::new();
    let mut cursor = head;
    while let Some(align) = cursor {
        out.push((
            align.query_start,
            align.query_end,
            align.match_start,
            align.match_end,
            align.score,
            align.frame,
        ));
        cursor = align.next.as_deref();
    }
    out
}

/// blast-rs: In-memory no-composition port boundary for `Blast_RedoOneMatch`
/// (`redo_alignment.c:1133`); not a complete direct NCBI C port.
///
/// This covers the branch where no composition adjustment is requested, so the
/// C callback calls reduce to already-ported Rust helpers:
/// [`sequence_get_range_in_memory_with_code`], [`redo_one_alignment`],
/// [`is_contained`], and [`with_distinct_ends`]. Use
/// [`blast_redo_one_match_in_memory_with_adjustment`] for the non-position
/// protein-space composition-adjusted materialized path.
#[allow(clippy::too_many_arguments)]
pub fn blast_redo_one_match_in_memory(
    alignments: &mut [Option<Box<BlastCompoAlignment>>],
    params: &BlastRedoAlignParams,
    incoming_aligns: &Option<Box<BlastCompoAlignment>>,
    _hspcnt: i32,
    lambda: f64,
    matching_seq: &BlastCompoMatchingSequence,
    query_info: &[BlastCompoQueryInfo],
    program: ProgramType,
    subject_source: &[u8],
    reward: i32,
    penalty: i32,
    genetic_code: &[u8; 64],
) -> i32 {
    let num_queries = query_info.len();
    if alignments.len() < num_queries {
        return -1;
    }
    alignments_free_array(alignments, num_queries);

    if !matches!(
        params.compo_adjust_mode,
        CompoAdjustMode::NoCompositionBasedStats
    ) {
        return -1;
    }
    if program != crate::program::BLASTN && matrix_info_to_aa_array(&params.matrix_info).is_none() {
        return -1;
    }

    let windows = match windows_from_aligns(
        incoming_aligns,
        query_info,
        _hspcnt,
        num_queries,
        K_WINDOW_BORDER,
        matching_seq.length,
        params.query_is_translated,
        params.subject_is_translated,
        params.position_based,
    ) {
        Ok(windows) => windows,
        Err(_) => return -1,
    };

    for window in &windows {
        let query_index = match window.align.as_ref() {
            Some(align) if align.query_index >= 0 => align.query_index as usize,
            _ => {
                alignments_free_array(alignments, num_queries);
                return -1;
            }
        };
        if query_index >= num_queries {
            alignments_free_array(alignments, num_queries);
            return -1;
        }

        let (query, subject) = match sequence_get_range_in_memory_with_code(
            program,
            &query_info[query_index].seq,
            &window.query_range,
            subject_source,
            &window.subject_range,
            genetic_code,
        ) {
            Ok(pair) => pair,
            Err(_) => {
                alignments_free_array(alignments, num_queries);
                return -1;
            }
        };

        let mut in_align = window.align.as_deref();
        while let Some(align) = in_align {
            let existing = alignment_list_to_containment_tuples(alignments[query_index].as_deref());
            if !is_contained(
                align.query_start,
                align.query_end,
                align.match_start,
                align.match_end,
                align.score,
                align.frame,
                &existing,
                lambda,
            ) {
                let new_align = if program == crate::program::BLASTN {
                    redo_one_alignment(
                        &query,
                        &window.query_range,
                        &subject,
                        &window.subject_range,
                        align.query_gapped_start(),
                        align.match_gapped_start(),
                        MatrixAdjustRule::DontAdjust,
                        reward,
                        penalty,
                        params.gapping_params.gap_open,
                        params.gapping_params.gap_extend,
                        params.gapping_params.x_dropoff,
                    )
                } else {
                    redo_one_alignment_protein(
                        &query,
                        &window.query_range,
                        &subject,
                        &window.subject_range,
                        align.query_gapped_start(),
                        align.match_gapped_start(),
                        MatrixAdjustRule::DontAdjust,
                        &params.matrix_info,
                        params.gapping_params.gap_open,
                        params.gapping_params.gap_extend,
                        params.gapping_params.x_dropoff,
                    )
                };
                if let Some(new_align) = new_align {
                    if new_align.score >= params.cutoff_s {
                        with_distinct_ends(new_align, &mut alignments[query_index], false);
                    }
                }
            }
            in_align = align.next.as_deref();
        }
    }

    0
}

/// blast-rs: In-memory protein-space composition-adjusted counterpart of
/// [`blast_redo_one_match_in_memory`]; not a direct NCBI C port.
///
/// This covers the materialized-subject ordinary redo branch once query and
/// subject ranges can be copied locally: compute query/subject composition,
/// run `Blast_AdjustScores`, redo the alignment with the adjusted matrix, and
/// insert through `s_WithDistinctEnds`. Nucleotide composition adjustment is
/// rejected by the caller; position-specific/PSSM redo uses
/// [`blast_redo_one_match_in_memory_with_position_adjustment`].
#[allow(clippy::too_many_arguments)]
pub fn blast_redo_one_match_in_memory_with_adjustment(
    alignments: &mut [Option<Box<BlastCompoAlignment>>],
    params: &BlastRedoAlignParams,
    incoming_aligns: &Option<Box<BlastCompoAlignment>>,
    _hspcnt: i32,
    lambda: f64,
    matching_seq: &BlastCompoMatchingSequence,
    query_info: &[BlastCompoQueryInfo],
    program: ProgramType,
    subject_source: &[u8],
    genetic_code: &[u8; 64],
    matrix: &mut [[i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE],
    workspace: Option<&BlastCompositionWorkspace>,
    pvalue_for_this_pair: &mut f64,
    composition_test_index: i32,
    lambda_ratio: &mut f64,
) -> i32 {
    let num_queries = query_info.len();
    if alignments.len() < num_queries {
        return -1;
    }
    alignments_free_array(alignments, num_queries);

    if matches!(
        params.compo_adjust_mode,
        CompoAdjustMode::NoCompositionBasedStats
    ) || params.position_based
        || program == crate::program::BLASTN
        || (params.query_is_translated && params.subject_is_translated)
    {
        return -1;
    }
    if matrix_info_to_aa_array(&params.matrix_info).is_none() {
        return -1;
    }

    let windows = match windows_from_aligns(
        incoming_aligns,
        query_info,
        _hspcnt,
        num_queries,
        K_WINDOW_BORDER,
        matching_seq.length,
        params.query_is_translated,
        params.subject_is_translated,
        params.position_based,
    ) {
        Ok(windows) => windows,
        Err(_) => return -1,
    };

    for window in &windows {
        let query_index = match window.align.as_ref() {
            Some(align) if align.query_index >= 0 => align.query_index as usize,
            _ => {
                alignments_free_array(alignments, num_queries);
                return -1;
            }
        };
        if query_index >= num_queries {
            alignments_free_array(alignments, num_queries);
            return -1;
        }

        let range_align = match window.align.as_deref() {
            Some(align) => align,
            None => {
                alignments_free_array(alignments, num_queries);
                return -1;
            }
        };
        let mut in_align = window.align.as_deref();
        let mut query = BlastCompoSequenceData::default();
        let mut subject = BlastCompoSequenceData::default();
        let mut near_identical_status = true;
        let mut old_near_identical_status = false;
        let mut subject_maybe_biased = true;
        let mut hsp_index = 0;
        let mut num_adjustments = 0;
        let mut query_composition = query_info[query_index].composition.prob.clone();
        let mut query_num_true = if query_composition.iter().any(|&p| p > 0.0) {
            crate::composition::blast_read_aa_composition(
                query_info[query_index].seq.data(),
                crate::matrix::AA_SIZE,
            )
            .1
        } else {
            0
        };
        while let Some(align) = in_align {
            if hsp_index == 0 || subject_maybe_biased {
                near_identical_status = preliminary_test_near_identical(
                    query_info,
                    window,
                    align,
                    params.near_identical_cutoff,
                );
            }

            if hsp_index == 0
                || (subject_maybe_biased && near_identical_status != old_near_identical_status)
            {
                sequence_data_release(&mut subject);
                sequence_data_release(&mut query);
                match sequence_get_range_in_memory_with_code_for_redo(
                    program,
                    &query_info[query_index].seq,
                    &window.query_range,
                    subject_source,
                    &window.subject_range,
                    genetic_code,
                    &query_info[query_index].words,
                    range_align,
                    near_identical_status,
                    params.compo_adjust_mode,
                    &mut subject_maybe_biased,
                ) {
                    Ok((query_seq, subject_seq)) => {
                        query = query_seq;
                        subject = subject_seq;
                    }
                    Err(_) => {
                        alignments_free_array(alignments, num_queries);
                        return -1;
                    }
                }
            }

            if params.query_is_translated {
                let (composition, num_true) = get_composition(
                    &query,
                    &window.query_range,
                    align,
                    crate::matrix::AA_SIZE,
                    true,
                    false,
                );
                query_composition = composition;
                query_num_true = num_true;
            }

            let existing = alignment_list_to_containment_tuples(alignments[query_index].as_deref());
            if !is_contained(
                align.query_start,
                align.query_end,
                align.match_start,
                align.match_end,
                align.score,
                align.frame,
                &existing,
                lambda,
            ) {
                let mut matrix_adjust_rule = MatrixAdjustRule::DontAdjust;
                let mut adjust_search_failed = 0;
                if !matches!(
                    params.compo_adjust_mode,
                    CompoAdjustMode::NoCompositionBasedStats
                ) && (params.subject_is_translated
                    || hsp_index == 0
                    || near_identical_status != old_near_identical_status)
                {
                    let (subject_composition, subject_num_true) = get_composition(
                        &subject,
                        &window.subject_range,
                        align,
                        crate::matrix::AA_SIZE,
                        false,
                        params.subject_is_translated,
                    );
                    // NCBI `redo_alignment.c:1242-1245` passes `query.length`
                    // and `subject.length` (window-data buffer lengths,
                    // including any X residues), NOT the `numTrueAminoAcids`
                    // counts. Match exactly.
                    let query_window_len = query.data().len();
                    let subject_window_len = subject.data().len();
                    match blast_adjust_scores_with_workspace_v2(
                        matrix,
                        &query_composition,
                        query_window_len,
                        query_num_true,
                        &subject_composition,
                        subject_window_len,
                        subject_num_true,
                        &params.matrix_info,
                        params.compo_adjust_mode,
                        params.re_pseudocounts,
                        workspace,
                        composition_test_index,
                        pvalue_for_this_pair,
                        lambda_ratio,
                    ) {
                        Ok(rule) => {
                            matrix_adjust_rule = rule;
                            num_adjustments += 1;
                        }
                        Err(status) if status < 0 => {
                            alignments_free_array(alignments, num_queries);
                            return status;
                        }
                        Err(_) => {
                            adjust_search_failed = 1;
                            num_adjustments += 1;
                        }
                    }
                }

                if adjust_search_failed == 0 {
                    if let Some(new_align) = redo_one_alignment_protein_with_matrix(
                        &query,
                        &window.query_range,
                        &subject,
                        &window.subject_range,
                        align.query_gapped_start(),
                        align.match_gapped_start(),
                        matrix_adjust_rule,
                        matrix,
                        params.gapping_params.gap_open,
                        params.gapping_params.gap_extend,
                        params.gapping_params.x_dropoff,
                    ) {
                        if new_align.score >= params.cutoff_s {
                            with_distinct_ends(
                                new_align,
                                &mut alignments[query_index],
                                num_adjustments == 1,
                            );
                        }
                    }
                }
            }
            old_near_identical_status = near_identical_status;
            hsp_index += 1;
            in_align = align.next.as_deref();
        }
        sequence_data_release(&mut subject);
        sequence_data_release(&mut query);
    }

    0
}

/// blast-rs: Position-specific/PSSM composition-adjusted counterpart of
/// [`blast_redo_one_match_in_memory_with_adjustment`]; not a direct NCBI C port.
///
/// This covers the ordinary PSI-style materialized redo branch: copy local
/// query/subject windows, compute composition, run the position-based
/// `Blast_AdjustScores` scaling branch, and redo traceback against the scaled
/// PSSM rows.
#[allow(clippy::too_many_arguments)]
pub fn blast_redo_one_match_in_memory_with_position_adjustment(
    alignments: &mut [Option<Box<BlastCompoAlignment>>],
    params: &BlastRedoAlignParams,
    incoming_aligns: &Option<Box<BlastCompoAlignment>>,
    _hspcnt: i32,
    lambda: f64,
    matching_seq: &BlastCompoMatchingSequence,
    query_info: &[BlastCompoQueryInfo],
    program: ProgramType,
    subject_source: &[u8],
    genetic_code: &[u8; 64],
    adjusted_pssm: &mut Vec<Vec<i32>>,
    pvalue_for_this_pair: &mut f64,
    composition_test_index: i32,
    lambda_ratio: &mut f64,
) -> i32 {
    let num_queries = query_info.len();
    if alignments.len() < num_queries {
        return -1;
    }
    alignments_free_array(alignments, num_queries);

    if matches!(
        params.compo_adjust_mode,
        CompoAdjustMode::NoCompositionBasedStats
    ) || !params.position_based
        || program == crate::program::BLASTN
        || params.query_is_translated
    {
        return -1;
    }
    if !params.matrix_info.positional {
        return -1;
    }

    let windows = match windows_from_aligns(
        incoming_aligns,
        query_info,
        _hspcnt,
        num_queries,
        K_WINDOW_BORDER,
        matching_seq.length,
        params.query_is_translated,
        params.subject_is_translated,
        params.position_based,
    ) {
        Ok(windows) => windows,
        Err(_) => return -1,
    };

    for window in &windows {
        let query_index = match window.align.as_ref() {
            Some(align) if align.query_index >= 0 => align.query_index as usize,
            _ => {
                alignments_free_array(alignments, num_queries);
                return -1;
            }
        };
        if query_index >= num_queries {
            alignments_free_array(alignments, num_queries);
            return -1;
        }

        let (query, subject) = match sequence_get_range_in_memory_with_code(
            program,
            &query_info[query_index].seq,
            &window.query_range,
            subject_source,
            &window.subject_range,
            genetic_code,
        ) {
            Ok(pair) => pair,
            Err(_) => {
                alignments_free_array(alignments, num_queries);
                return -1;
            }
        };

        let mut in_align = window.align.as_deref();
        let mut query_composition = query_info[query_index].composition.prob.clone();
        let mut query_num_true = if query_composition.iter().any(|&p| p > 0.0) {
            crate::composition::blast_read_aa_composition(
                query_info[query_index].seq.data(),
                crate::matrix::AA_SIZE,
            )
            .1
        } else {
            0
        };
        while let Some(align) = in_align {
            if query_num_true == 0 {
                let (composition, num_true) = get_composition(
                    &query,
                    &window.query_range,
                    align,
                    crate::matrix::AA_SIZE,
                    true,
                    false,
                );
                query_composition = composition;
                query_num_true = num_true;
            }

            let existing = alignment_list_to_containment_tuples(alignments[query_index].as_deref());
            if !is_contained(
                align.query_start,
                align.query_end,
                align.match_start,
                align.match_end,
                align.score,
                align.frame,
                &existing,
                lambda,
            ) {
                let (subject_composition, subject_num_true) = get_composition(
                    &subject,
                    &window.subject_range,
                    align,
                    crate::matrix::AA_SIZE,
                    false,
                    params.subject_is_translated,
                );
                let matrix_adjust_rule = match blast_adjust_position_based_scores(
                    adjusted_pssm,
                    &query_composition,
                    query_num_true,
                    &subject_composition,
                    subject_num_true,
                    &params.matrix_info,
                    composition_test_index,
                    pvalue_for_this_pair,
                    lambda_ratio,
                ) {
                    Ok(rule) => rule,
                    Err(status) if status < 0 => {
                        alignments_free_array(alignments, num_queries);
                        return status;
                    }
                    Err(_) => {
                        in_align = align.next.as_deref();
                        continue;
                    }
                };

                if let Some(new_align) = redo_one_alignment_protein_with_pssm(
                    &query,
                    &window.query_range,
                    &subject,
                    &window.subject_range,
                    align.query_gapped_start(),
                    align.match_gapped_start(),
                    matrix_adjust_rule,
                    adjusted_pssm,
                    params.gapping_params.gap_open,
                    params.gapping_params.gap_extend,
                    params.gapping_params.x_dropoff,
                ) {
                    if new_align.score >= params.cutoff_s {
                        with_distinct_ends(new_align, &mut alignments[query_index], false);
                    }
                }
            }
            in_align = align.next.as_deref();
        }
    }

    0
}

/// blast-rs: Callback-driven no-composition port boundary for `Blast_RedoOneMatch`
/// (`redo_alignment.c:1133`); not a complete direct NCBI C port.
///
/// This follows the generic C driver shape: clear output lists, build windows,
/// rerun `get_range` when the near-identical state changes, skip HSPs already
/// contained in better alignments, call `redo_one_alignment`, and insert via
/// `s_WithDistinctEnds`. Use
/// [`blast_redo_one_match_with_callbacks_and_adjustment`] for the
/// callback-driven composition-adjusted path.
pub fn blast_redo_one_match_with_callbacks(
    alignments: &mut [Option<Box<BlastCompoAlignment>>],
    params: &BlastRedoAlignParams,
    incoming_aligns: &Option<Box<BlastCompoAlignment>>,
    _hspcnt: i32,
    lambda: f64,
    matching_seq: &BlastCompoMatchingSequence,
    query_info: &[BlastCompoQueryInfo],
) -> i32 {
    let mut matrix = [[0i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE];
    let mut pvalue_for_this_pair = -1.0;
    let mut lambda_ratio = 1.0;
    blast_redo_one_match_with_callbacks_and_adjustment(
        alignments,
        params,
        incoming_aligns,
        _hspcnt,
        lambda,
        matching_seq,
        query_info,
        &mut matrix,
        None,
        &mut pvalue_for_this_pair,
        0,
        &mut lambda_ratio,
    )
}

/// blast-rs: Callback-driven port boundary for `Blast_RedoOneMatch`
/// (`redo_alignment.c:1133`) with the composition-adjusted score path wired; not a complete direct NCBI C port.
///
/// The existing callback table supplies `get_range` and `redo_one_alignment`,
/// while this function owns the translated `Blast_AdjustScores` call between
/// them. This covers non-position square matrices and delegates
/// position-based/PSSM scoring to
/// [`blast_redo_one_match_with_callbacks_and_position_adjustment`].
#[allow(clippy::too_many_arguments)]
pub fn blast_redo_one_match_with_callbacks_and_adjustment(
    alignments: &mut [Option<Box<BlastCompoAlignment>>],
    params: &BlastRedoAlignParams,
    incoming_aligns: &Option<Box<BlastCompoAlignment>>,
    _hspcnt: i32,
    lambda: f64,
    matching_seq: &BlastCompoMatchingSequence,
    query_info: &[BlastCompoQueryInfo],
    matrix: &mut [[i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE],
    workspace: Option<&BlastCompositionWorkspace>,
    pvalue_for_this_pair: &mut f64,
    composition_test_index: i32,
    lambda_ratio: &mut f64,
) -> i32 {
    let num_queries = query_info.len();
    if alignments.len() < num_queries {
        return -1;
    }
    alignments_free_array(alignments, num_queries);

    let Some(callbacks) = params.callbacks else {
        return -1;
    };

    if params.position_based {
        let mut adjusted_pssm = Vec::new();
        return blast_redo_one_match_with_callbacks_and_position_adjustment(
            alignments,
            params,
            incoming_aligns,
            _hspcnt,
            lambda,
            matching_seq,
            query_info,
            &mut adjusted_pssm,
            pvalue_for_this_pair,
            composition_test_index,
            lambda_ratio,
        );
    }
    let (Some(get_range), Some(redo_one_alignment_cb)) =
        (callbacks.get_range, callbacks.redo_one_alignment)
    else {
        return -1;
    };
    if params.query_is_translated && params.subject_is_translated {
        return -1;
    }

    let windows = match windows_from_aligns(
        incoming_aligns,
        query_info,
        _hspcnt,
        num_queries,
        K_WINDOW_BORDER,
        matching_seq.length,
        params.query_is_translated,
        params.subject_is_translated,
        params.position_based,
    ) {
        Ok(windows) => windows,
        Err(_) => return -1,
    };

    for window in &windows {
        let query_index = match window.align.as_ref() {
            Some(align) if align.query_index >= 0 => align.query_index as usize,
            _ => {
                alignments_free_array(alignments, num_queries);
                return -1;
            }
        };
        if query_index >= num_queries {
            alignments_free_array(alignments, num_queries);
            return -1;
        }

        let range_align = match window.align.as_deref() {
            Some(align) => align,
            None => {
                alignments_free_array(alignments, num_queries);
                return -1;
            }
        };
        let mut subject = BlastCompoSequenceData::default();
        let mut query = BlastCompoSequenceData::default();
        let mut near_identical_status = true;
        let mut old_near_identical_status = false;
        let mut subject_maybe_biased = true;
        let mut num_adjustments = 0;
        let mut hsp_index = 0;
        let mut in_align = window.align.as_deref();
        let mut query_composition = query_info[query_index].composition.prob.clone();
        let mut query_num_true = if query_composition.iter().any(|&p| p > 0.0) {
            crate::composition::blast_read_aa_composition(
                query_info[query_index].seq.data(),
                crate::matrix::AA_SIZE,
            )
            .1
        } else {
            0
        };
        while let Some(align) = in_align {
            if hsp_index == 0 || subject_maybe_biased {
                near_identical_status = preliminary_test_near_identical(
                    query_info,
                    window,
                    align,
                    params.near_identical_cutoff,
                );
            }

            if hsp_index == 0
                || (subject_maybe_biased && near_identical_status != old_near_identical_status)
            {
                sequence_data_release(&mut subject);
                sequence_data_release(&mut query);
                let status = get_range(
                    matching_seq,
                    &window.subject_range,
                    &mut subject,
                    &query_info[query_index].seq,
                    &window.query_range,
                    &mut query,
                    &query_info[query_index].words,
                    range_align,
                    near_identical_status,
                    params.compo_adjust_mode,
                    false,
                    &mut subject_maybe_biased,
                );
                if status != 0 {
                    alignments_free_array(alignments, num_queries);
                    return status;
                }
            }

            if params.query_is_translated {
                let (composition, num_true) = get_composition(
                    &query,
                    &window.query_range,
                    align,
                    crate::matrix::AA_SIZE,
                    true,
                    false,
                );
                query_composition = composition;
                query_num_true = num_true;
            }

            let existing = alignment_list_to_containment_tuples(alignments[query_index].as_deref());
            if !is_contained(
                align.query_start,
                align.query_end,
                align.match_start,
                align.match_end,
                align.score,
                align.frame,
                &existing,
                lambda,
            ) {
                let mut matrix_adjust_rule = MatrixAdjustRule::DontAdjust;
                let mut adjust_search_failed = 0;
                if !matches!(
                    params.compo_adjust_mode,
                    CompoAdjustMode::NoCompositionBasedStats
                ) && (params.subject_is_translated
                    || hsp_index == 0
                    || near_identical_status != old_near_identical_status)
                {
                    let (subject_composition, subject_num_true) = get_composition(
                        &subject,
                        &window.subject_range,
                        align,
                        crate::matrix::AA_SIZE,
                        false,
                        params.subject_is_translated,
                    );
                    // NCBI passes full window-data lengths (not num_true) to
                    // `Blast_AdjustScores` — see `redo_alignment.c:1242`.
                    let query_window_len = query.data().len();
                    let subject_window_len = subject.data().len();
                    match blast_adjust_scores_with_workspace_v2(
                        matrix,
                        &query_composition,
                        query_window_len,
                        query_num_true,
                        &subject_composition,
                        subject_window_len,
                        subject_num_true,
                        &params.matrix_info,
                        params.compo_adjust_mode,
                        params.re_pseudocounts,
                        workspace,
                        composition_test_index,
                        pvalue_for_this_pair,
                        lambda_ratio,
                    ) {
                        Ok(rule) => {
                            matrix_adjust_rule = rule;
                            num_adjustments += 1;
                        }
                        Err(status) if status < 0 => {
                            alignments_free_array(alignments, num_queries);
                            return status;
                        }
                        Err(_) => {
                            adjust_search_failed = 1;
                            num_adjustments += 1;
                        }
                    }
                }

                if adjust_search_failed == 0 {
                    if let Some(new_align) = redo_one_alignment_cb(
                        align,
                        matrix_adjust_rule,
                        &query,
                        &window.query_range,
                        params.ccat_query_length,
                        &subject,
                        &window.subject_range,
                        matching_seq.length,
                        &params.gapping_params,
                    ) {
                        if new_align.score >= params.cutoff_s {
                            with_distinct_ends(
                                new_align,
                                &mut alignments[query_index],
                                num_adjustments == 1,
                            );
                        }
                    }
                }
            }

            old_near_identical_status = near_identical_status;
            hsp_index += 1;
            in_align = align.next.as_deref();
        }
        sequence_data_release(&mut subject);
        sequence_data_release(&mut query);
    }

    0
}

/// blast-rs: Callback-driven position-specific/PSSM counterpart of
/// [`blast_redo_one_match_with_callbacks_and_adjustment`]; not a direct NCBI C port.
///
/// The callback table supplies the C-style range fetcher; after
/// position-based score adjustment this helper runs the local PSSM traceback
/// against the scaled rows.
#[allow(clippy::too_many_arguments)]
pub fn blast_redo_one_match_with_callbacks_and_position_adjustment(
    alignments: &mut [Option<Box<BlastCompoAlignment>>],
    params: &BlastRedoAlignParams,
    incoming_aligns: &Option<Box<BlastCompoAlignment>>,
    _hspcnt: i32,
    lambda: f64,
    matching_seq: &BlastCompoMatchingSequence,
    query_info: &[BlastCompoQueryInfo],
    adjusted_pssm: &mut Vec<Vec<i32>>,
    pvalue_for_this_pair: &mut f64,
    composition_test_index: i32,
    lambda_ratio: &mut f64,
) -> i32 {
    let num_queries = query_info.len();
    if alignments.len() < num_queries {
        return -1;
    }
    alignments_free_array(alignments, num_queries);

    let Some(callbacks) = params.callbacks else {
        return -1;
    };
    let Some(get_range) = callbacks.get_range else {
        return -1;
    };

    if matches!(
        params.compo_adjust_mode,
        CompoAdjustMode::NoCompositionBasedStats
    ) || !params.position_based
        || !params.matrix_info.positional
        || params.query_is_translated
        || params.subject_is_translated
    {
        return -1;
    }

    let windows = match windows_from_aligns(
        incoming_aligns,
        query_info,
        _hspcnt,
        num_queries,
        K_WINDOW_BORDER,
        matching_seq.length,
        params.query_is_translated,
        params.subject_is_translated,
        params.position_based,
    ) {
        Ok(windows) => windows,
        Err(_) => return -1,
    };

    for window in &windows {
        let query_index = match window.align.as_ref() {
            Some(align) if align.query_index >= 0 => align.query_index as usize,
            _ => {
                alignments_free_array(alignments, num_queries);
                return -1;
            }
        };
        if query_index >= num_queries {
            alignments_free_array(alignments, num_queries);
            return -1;
        }

        let mut subject = BlastCompoSequenceData::default();
        let mut query = BlastCompoSequenceData::default();
        let mut near_identical_status = true;
        let mut old_near_identical_status = false;
        let mut subject_maybe_biased = true;
        let mut num_adjustments = 0;
        let mut hsp_index = 0;
        let mut in_align = window.align.as_deref();
        let mut query_composition = query_info[query_index].composition.prob.clone();
        let mut query_num_true = if query_composition.iter().any(|&p| p > 0.0) {
            query_info[query_index].seq.length.max(0) as usize
        } else {
            0
        };

        while let Some(align) = in_align {
            if hsp_index == 0 || subject_maybe_biased {
                near_identical_status = preliminary_test_near_identical(
                    query_info,
                    window,
                    align,
                    params.near_identical_cutoff,
                );
            }

            if hsp_index == 0
                || (subject_maybe_biased && near_identical_status != old_near_identical_status)
            {
                sequence_data_release(&mut subject);
                sequence_data_release(&mut query);
                let status = get_range(
                    matching_seq,
                    &window.subject_range,
                    &mut subject,
                    &query_info[query_index].seq,
                    &window.query_range,
                    &mut query,
                    &query_info[query_index].words,
                    align,
                    near_identical_status,
                    params.compo_adjust_mode,
                    false,
                    &mut subject_maybe_biased,
                );
                if status != 0 {
                    alignments_free_array(alignments, num_queries);
                    return status;
                }
            }

            if query_num_true == 0 {
                let (composition, num_true) = get_composition(
                    &query,
                    &window.query_range,
                    align,
                    crate::matrix::AA_SIZE,
                    true,
                    false,
                );
                query_composition = composition;
                query_num_true = num_true;
            }

            let existing = alignment_list_to_containment_tuples(alignments[query_index].as_deref());
            if !is_contained(
                align.query_start,
                align.query_end,
                align.match_start,
                align.match_end,
                align.score,
                align.frame,
                &existing,
                lambda,
            ) {
                let mut matrix_adjust_rule = MatrixAdjustRule::DontAdjust;
                let mut adjust_search_failed = 0;
                if hsp_index == 0 || near_identical_status != old_near_identical_status {
                    let (subject_composition, subject_num_true) = get_composition(
                        &subject,
                        &window.subject_range,
                        align,
                        crate::matrix::AA_SIZE,
                        false,
                        false,
                    );
                    match blast_adjust_position_based_scores(
                        adjusted_pssm,
                        &query_composition,
                        query_num_true,
                        &subject_composition,
                        subject_num_true,
                        &params.matrix_info,
                        composition_test_index,
                        pvalue_for_this_pair,
                        lambda_ratio,
                    ) {
                        Ok(rule) => {
                            matrix_adjust_rule = rule;
                            num_adjustments += 1;
                        }
                        Err(status) if status < 0 => {
                            alignments_free_array(alignments, num_queries);
                            return status;
                        }
                        Err(_) => {
                            adjust_search_failed = 1;
                            num_adjustments += 1;
                        }
                    }
                }

                if adjust_search_failed == 0 {
                    if let Some(new_align) = redo_one_alignment_protein_with_pssm(
                        &query,
                        &window.query_range,
                        &subject,
                        &window.subject_range,
                        align.query_gapped_start(),
                        align.match_gapped_start(),
                        matrix_adjust_rule,
                        adjusted_pssm,
                        params.gapping_params.gap_open,
                        params.gapping_params.gap_extend,
                        params.gapping_params.x_dropoff,
                    ) {
                        if new_align.score >= params.cutoff_s {
                            with_distinct_ends(
                                new_align,
                                &mut alignments[query_index],
                                num_adjustments == 1,
                            );
                        }
                    }
                }
            }

            old_near_identical_status = near_identical_status;
            hsp_index += 1;
            in_align = align.next.as_deref();
        }
        sequence_data_release(&mut subject);
        sequence_data_release(&mut query);
    }

    0
}

/// blast-rs: Converts the local nucleotide X-drop matrix into SW matrix shape; not a direct NCBI C port.
fn nucleotide_sw_matrix_from_xdrop_matrix(
    matrix: &[[i32; 16]; 16],
) -> [[i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE] {
    let mut out = [[0i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE];
    for i in 0..16 {
        for j in 0..16 {
            out[i][j] = matrix[i][j];
        }
    }
    out
}

/// blast-rs: In-memory no-composition nucleotide port boundary for
/// `Blast_RedoOneMatchSmithWaterman` (`redo_alignment.c:1342`); not a complete direct NCBI C port.
///
/// This mirrors the executable SW loop once sequence ranges are already
/// materialized locally: build windows, retrieve query/subject ranges, clear
/// forbidden ranges, run score-only SW, apply the NCBI significance gate, find
/// the start with reverse SW, finalize with the X-drop callback, prepend the
/// redone alignment, and mark the accepted rectangle forbidden when the input
/// window has multiple HSPs. Composition adjustment and translated subjects
/// are intentionally routed through the protein-space SW helpers instead.
#[allow(clippy::too_many_arguments)]
pub fn blast_redo_one_match_smith_waterman_in_memory_nucl(
    alignments: &mut [Option<Box<BlastCompoAlignment>>],
    params: &BlastRedoAlignParams,
    incoming_aligns: &Option<Box<BlastCompoAlignment>>,
    hspcnt: i32,
    lambda: f64,
    log_k: f64,
    matching_seq: &BlastCompoMatchingSequence,
    query_info: &[BlastCompoQueryInfo],
    subject_source: &[u8],
    significant_matches: &[BlastCompoHeap],
    matrix: &[[i32; 16]; 16],
) -> i32 {
    let num_queries = query_info.len();
    if alignments.len() < num_queries || significant_matches.len() < num_queries {
        return -1;
    }
    alignments_free_array(alignments, num_queries);

    if !matches!(
        params.compo_adjust_mode,
        CompoAdjustMode::NoCompositionBasedStats
    ) || params.query_is_translated
        || params.subject_is_translated
        || params.position_based
    {
        return -1;
    }

    let windows = match windows_from_aligns(
        incoming_aligns,
        query_info,
        hspcnt,
        num_queries,
        K_WINDOW_BORDER,
        matching_seq.length,
        false,
        false,
        false,
    ) {
        Ok(windows) => windows,
        Err(_) => return -1,
    };

    let sw_matrix = nucleotide_sw_matrix_from_xdrop_matrix(matrix);
    for window in &windows {
        let query_index = window.query_range.context;
        if query_index < 0 || query_index as usize >= num_queries {
            alignments_free_array(alignments, num_queries);
            return -1;
        }
        let query_index = query_index as usize;

        let (query, subject) = match sequence_get_range_in_memory(
            crate::program::BLASTN,
            &query_info[query_index].seq,
            &window.query_range,
            subject_source,
            &window.subject_range,
        ) {
            Ok(pair) => pair,
            Err(_) => {
                alignments_free_array(alignments, num_queries);
                return -1;
            }
        };
        let searchsp = query_info[query_index].eff_search_space;
        let mut forbidden = crate::smith_waterman::BlastForbiddenRanges::new(query.length);
        if hspcnt <= 0 {
            continue;
        }
        loop {
            let (sw_score, match_end, query_end) =
                crate::smith_waterman::blast_smith_waterman_score_only_with_forbidden(
                    subject.data(),
                    query.data(),
                    &sw_matrix,
                    params.gapping_params.gap_open,
                    params.gapping_params.gap_extend,
                    &forbidden,
                );
            let significant = smith_waterman_alignment_is_significant(
                sw_score,
                lambda,
                log_k,
                searchsp,
                params,
                alignments[query_index].as_deref(),
                &significant_matches[query_index],
                matching_seq.index,
            );
            if !significant {
                break;
            }

            let (_updated_score, match_start, query_start) =
                crate::smith_waterman::blast_smith_waterman_find_start_with_forbidden(
                    subject.data(),
                    query.data(),
                    &sw_matrix,
                    params.gapping_params.gap_open,
                    params.gapping_params.gap_extend,
                    match_end,
                    query_end,
                    sw_score,
                    &forbidden,
                );

            let Some((mut new_align, query_end_xdrop, match_end_xdrop)) = new_alignment_using_xdrop(
                &query,
                &window.query_range,
                &subject,
                &window.subject_range,
                query_start,
                query_end,
                match_start,
                match_end,
                sw_score,
                &params.gapping_params,
                MatrixAdjustRule::DontAdjust,
                matrix,
            ) else {
                alignments_free_array(alignments, num_queries);
                return -1;
            };
            new_align.next = alignments[query_index].take();
            alignments[query_index] = Some(Box::new(new_align));

            if window.hspcnt > 1 {
                forbidden.push(
                    query_start as i32,
                    query_end_xdrop + 1,
                    match_start as i32,
                    match_end_xdrop,
                );
            }
            if !(significant && window.hspcnt > 1) {
                break;
            }
        }
    }

    0
}

/// blast-rs: In-memory no-composition protein-space port boundary for
/// `Blast_RedoOneMatchSmithWaterman` (`redo_alignment.c:1342`); not a complete direct NCBI C port.
///
/// This mirrors [`blast_redo_one_match_smith_waterman_in_memory_nucl`] after
/// sequence ranges have been materialized as protein residues. It covers
/// BLASTP plus translated-query/translated-subject protein-space searches.
/// Composition-adjusted and position-specific matrix redo use the dedicated
/// materialized adjustment helpers.
#[allow(clippy::too_many_arguments)]
pub fn blast_redo_one_match_smith_waterman_in_memory_protein(
    alignments: &mut [Option<Box<BlastCompoAlignment>>],
    params: &BlastRedoAlignParams,
    incoming_aligns: &Option<Box<BlastCompoAlignment>>,
    hspcnt: i32,
    lambda: f64,
    log_k: f64,
    matching_seq: &BlastCompoMatchingSequence,
    query_info: &[BlastCompoQueryInfo],
    program: ProgramType,
    subject_source: &[u8],
    significant_matches: &[BlastCompoHeap],
    genetic_code: &[u8; 64],
) -> i32 {
    let num_queries = query_info.len();
    if alignments.len() < num_queries || significant_matches.len() < num_queries {
        return -1;
    }
    alignments_free_array(alignments, num_queries);

    if !matches!(
        params.compo_adjust_mode,
        CompoAdjustMode::NoCompositionBasedStats
    ) || params.position_based
        || program == crate::program::BLASTN
    {
        return -1;
    }
    let Some(sw_matrix) = matrix_info_to_aa_array(&params.matrix_info) else {
        return -1;
    };

    let windows = match windows_from_aligns(
        incoming_aligns,
        query_info,
        hspcnt,
        num_queries,
        K_WINDOW_BORDER,
        matching_seq.length,
        params.query_is_translated,
        params.subject_is_translated,
        params.position_based,
    ) {
        Ok(windows) => windows,
        Err(_) => return -1,
    };

    for window in &windows {
        let query_index = window.query_range.context;
        if query_index < 0 || query_index as usize >= num_queries {
            alignments_free_array(alignments, num_queries);
            return -1;
        }
        let query_index = query_index as usize;

        let (query, subject) = match sequence_get_range_in_memory_with_code(
            program,
            &query_info[query_index].seq,
            &window.query_range,
            subject_source,
            &window.subject_range,
            genetic_code,
        ) {
            Ok(pair) => pair,
            Err(_) => {
                alignments_free_array(alignments, num_queries);
                return -1;
            }
        };
        let searchsp = query_info[query_index].eff_search_space;
        let mut forbidden = crate::smith_waterman::BlastForbiddenRanges::new(query.length);
        if hspcnt <= 0 {
            continue;
        }
        loop {
            let (sw_score, match_end, query_end) =
                crate::smith_waterman::blast_smith_waterman_score_only_with_forbidden(
                    subject.data(),
                    query.data(),
                    &sw_matrix,
                    params.gapping_params.gap_open,
                    params.gapping_params.gap_extend,
                    &forbidden,
                );
            let significant = smith_waterman_alignment_is_significant(
                sw_score,
                lambda,
                log_k,
                searchsp,
                params,
                alignments[query_index].as_deref(),
                &significant_matches[query_index],
                matching_seq.index,
            );
            if !significant {
                break;
            }

            let (_updated_score, match_start, query_start) =
                crate::smith_waterman::blast_smith_waterman_find_start_with_forbidden(
                    subject.data(),
                    query.data(),
                    &sw_matrix,
                    params.gapping_params.gap_open,
                    params.gapping_params.gap_extend,
                    match_end,
                    query_end,
                    sw_score,
                    &forbidden,
                );

            let Some((mut new_align, query_end_xdrop, match_end_xdrop)) =
                new_alignment_using_xdrop_protein(
                    &query,
                    &window.query_range,
                    &subject,
                    &window.subject_range,
                    query_start,
                    query_end,
                    match_start,
                    match_end,
                    sw_score,
                    &params.gapping_params,
                    MatrixAdjustRule::DontAdjust,
                    &sw_matrix,
                )
            else {
                alignments_free_array(alignments, num_queries);
                return -1;
            };
            new_align.next = alignments[query_index].take();
            alignments[query_index] = Some(Box::new(new_align));

            if window.hspcnt > 1 {
                forbidden.push(
                    query_start as i32,
                    query_end_xdrop + 1,
                    match_start as i32,
                    match_end_xdrop,
                );
            }
            if !(significant && window.hspcnt > 1) {
                break;
            }
        }
    }

    0
}

/// blast-rs: In-memory non-position composition-adjusted protein-space port boundary for
/// `Blast_RedoOneMatchSmithWaterman`; not a complete direct NCBI C port.
///
/// This is the materialized-subject counterpart of the callback adjusted SW
/// path: copy protein-space query/subject ranges, compute compositions, run
/// square-matrix `Blast_AdjustScores`, then execute score-only SW, reverse
/// start finding, and bounded X-drop finalization with the adjusted matrix.
#[allow(clippy::too_many_arguments)]
pub fn blast_redo_one_match_smith_waterman_in_memory_protein_with_adjustment(
    alignments: &mut [Option<Box<BlastCompoAlignment>>],
    params: &BlastRedoAlignParams,
    incoming_aligns: &Option<Box<BlastCompoAlignment>>,
    hspcnt: i32,
    lambda: f64,
    log_k: f64,
    matching_seq: &BlastCompoMatchingSequence,
    query_info: &[BlastCompoQueryInfo],
    program: ProgramType,
    subject_source: &[u8],
    significant_matches: &[BlastCompoHeap],
    genetic_code: &[u8; 64],
    matrix: &mut [[i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE],
    workspace: Option<&BlastCompositionWorkspace>,
    pvalue_for_this_pair: &mut f64,
    composition_test_index: i32,
    lambda_ratio: &mut f64,
) -> i32 {
    let num_queries = query_info.len();
    if alignments.len() < num_queries || significant_matches.len() < num_queries {
        return -1;
    }
    alignments_free_array(alignments, num_queries);

    if matches!(
        params.compo_adjust_mode,
        CompoAdjustMode::NoCompositionBasedStats
    ) || params.position_based
        || program == crate::program::BLASTN
        || (params.query_is_translated && params.subject_is_translated)
    {
        return -1;
    }
    if matrix_info_to_aa_array(&params.matrix_info).is_none() {
        return -1;
    }

    let windows = match windows_from_aligns(
        incoming_aligns,
        query_info,
        hspcnt,
        num_queries,
        K_WINDOW_BORDER,
        matching_seq.length,
        params.query_is_translated,
        params.subject_is_translated,
        params.position_based,
    ) {
        Ok(windows) => windows,
        Err(_) => return -1,
    };

    for window in &windows {
        let query_index = window.query_range.context;
        if query_index < 0 || query_index as usize >= num_queries {
            alignments_free_array(alignments, num_queries);
            return -1;
        }
        let query_index = query_index as usize;
        let Some(window_align) = window.align.as_deref() else {
            alignments_free_array(alignments, num_queries);
            return -1;
        };

        let near_identical_status = preliminary_test_near_identical(
            query_info,
            window,
            window_align,
            params.near_identical_cutoff,
        );
        let (query, subject) = match sequence_get_range_in_memory_with_code_for_sw_redo(
            program,
            &query_info[query_index].seq,
            &window.query_range,
            subject_source,
            &window.subject_range,
            genetic_code,
            &query_info[query_index].words,
            window_align,
            near_identical_status,
            params.compo_adjust_mode,
        ) {
            Ok(pair) => pair,
            Err(_) => {
                alignments_free_array(alignments, num_queries);
                return -1;
            }
        };

        let mut query_composition = query_info[query_index].composition.prob.clone();
        let mut query_num_true = if query_composition.iter().any(|&p| p > 0.0) {
            query_info[query_index].seq.length.max(0) as usize
        } else {
            0
        };
        if params.query_is_translated || query_num_true == 0 {
            let (composition, num_true) = get_composition(
                &query,
                &window.query_range,
                window_align,
                crate::matrix::AA_SIZE,
                true,
                false,
            );
            query_composition = composition;
            query_num_true = num_true;
        }

        let (subject_composition, subject_num_true) = get_composition(
            &subject,
            &window.subject_range,
            window_align,
            crate::matrix::AA_SIZE,
            false,
            params.subject_is_translated,
        );
        // NCBI passes full window-data lengths to `Blast_AdjustScores`.
        let query_window_len = query.data().len();
        let subject_window_len = subject.data().len();
        let matrix_adjust_rule = match blast_adjust_scores_with_workspace_v2(
            matrix,
            &query_composition,
            query_window_len,
            query_num_true,
            &subject_composition,
            subject_window_len,
            subject_num_true,
            &params.matrix_info,
            params.compo_adjust_mode,
            params.re_pseudocounts,
            workspace,
            composition_test_index,
            pvalue_for_this_pair,
            lambda_ratio,
        ) {
            Ok(rule) => rule,
            Err(status) if status < 0 => {
                alignments_free_array(alignments, num_queries);
                return status;
            }
            Err(_) => continue,
        };

        let searchsp = query_info[query_index].eff_search_space;
        let mut forbidden = crate::smith_waterman::BlastForbiddenRanges::new(query.length);
        if hspcnt <= 0 {
            continue;
        }
        loop {
            let (sw_score, match_end, query_end) =
                crate::smith_waterman::blast_smith_waterman_score_only_with_forbidden(
                    subject.data(),
                    query.data(),
                    matrix,
                    params.gapping_params.gap_open,
                    params.gapping_params.gap_extend,
                    &forbidden,
                );
            let significant = smith_waterman_alignment_is_significant(
                sw_score,
                lambda,
                log_k,
                searchsp,
                params,
                alignments[query_index].as_deref(),
                &significant_matches[query_index],
                matching_seq.index,
            );
            if !significant {
                break;
            }

            let (_updated_score, match_start, query_start) =
                crate::smith_waterman::blast_smith_waterman_find_start_with_forbidden(
                    subject.data(),
                    query.data(),
                    matrix,
                    params.gapping_params.gap_open,
                    params.gapping_params.gap_extend,
                    match_end,
                    query_end,
                    sw_score,
                    &forbidden,
                );

            let Some((mut new_align, query_end_xdrop, match_end_xdrop)) =
                new_alignment_using_xdrop_protein(
                    &query,
                    &window.query_range,
                    &subject,
                    &window.subject_range,
                    query_start,
                    query_end,
                    match_start,
                    match_end,
                    sw_score,
                    &params.gapping_params,
                    matrix_adjust_rule,
                    matrix,
                )
            else {
                alignments_free_array(alignments, num_queries);
                return -1;
            };
            new_align.next = alignments[query_index].take();
            alignments[query_index] = Some(Box::new(new_align));

            if window.hspcnt > 1 {
                forbidden.push(
                    query_start as i32,
                    query_end_xdrop + 1,
                    match_start as i32,
                    match_end_xdrop,
                );
            }
            if !(significant && window.hspcnt > 1) {
                break;
            }
        }
    }

    0
}

/// blast-rs: In-memory position-specific/PSSM protein-space port boundary for
/// `Blast_RedoOneMatchSmithWaterman`; not a complete direct NCBI C port.
///
/// This mirrors the protein Smith-Waterman materialized loop, but adjusts the
/// position-specific rows for the current query/subject window before running
/// score-only SW, reverse start finding, and bounded X-drop finalization.
#[allow(clippy::too_many_arguments)]
pub fn blast_redo_one_match_smith_waterman_in_memory_protein_position_adjustment(
    alignments: &mut [Option<Box<BlastCompoAlignment>>],
    params: &BlastRedoAlignParams,
    incoming_aligns: &Option<Box<BlastCompoAlignment>>,
    hspcnt: i32,
    lambda: f64,
    log_k: f64,
    matching_seq: &BlastCompoMatchingSequence,
    query_info: &[BlastCompoQueryInfo],
    program: ProgramType,
    subject_source: &[u8],
    significant_matches: &[BlastCompoHeap],
    genetic_code: &[u8; 64],
    adjusted_pssm: &mut Vec<Vec<i32>>,
    pvalue_for_this_pair: &mut f64,
    composition_test_index: i32,
    lambda_ratio: &mut f64,
) -> i32 {
    let num_queries = query_info.len();
    if alignments.len() < num_queries || significant_matches.len() < num_queries {
        return -1;
    }
    alignments_free_array(alignments, num_queries);

    if matches!(
        params.compo_adjust_mode,
        CompoAdjustMode::NoCompositionBasedStats
    ) || !params.position_based
        || !params.matrix_info.positional
        || params.query_is_translated
        || program == crate::program::BLASTN
    {
        return -1;
    }

    let windows = match windows_from_aligns(
        incoming_aligns,
        query_info,
        hspcnt,
        num_queries,
        K_WINDOW_BORDER,
        matching_seq.length,
        params.query_is_translated,
        params.subject_is_translated,
        params.position_based,
    ) {
        Ok(windows) => windows,
        Err(_) => return -1,
    };

    for window in &windows {
        let query_index = match window.align.as_ref() {
            Some(align) if align.query_index >= 0 => align.query_index as usize,
            _ => {
                alignments_free_array(alignments, num_queries);
                return -1;
            }
        };
        if query_index >= num_queries {
            alignments_free_array(alignments, num_queries);
            return -1;
        }
        let Some(window_align) = window.align.as_deref() else {
            alignments_free_array(alignments, num_queries);
            return -1;
        };
        let near_identical_status = preliminary_test_near_identical(
            query_info,
            window,
            window_align,
            params.near_identical_cutoff,
        );
        let (query, subject) = match sequence_get_range_in_memory_with_code_for_sw_redo(
            program,
            &query_info[query_index].seq,
            &window.query_range,
            subject_source,
            &window.subject_range,
            genetic_code,
            &query_info[query_index].words,
            window_align,
            near_identical_status,
            params.compo_adjust_mode,
        ) {
            Ok(pair) => pair,
            Err(_) => {
                alignments_free_array(alignments, num_queries);
                return -1;
            }
        };

        let mut query_composition = query_info[query_index].composition.prob.clone();
        let mut query_num_true = if query_composition.iter().any(|&p| p > 0.0) {
            query_info[query_index].seq.length.max(0) as usize
        } else {
            0
        };
        if query_num_true == 0 {
            let (composition, num_true) = get_composition(
                &query,
                &window.query_range,
                window_align,
                crate::matrix::AA_SIZE,
                true,
                false,
            );
            query_composition = composition;
            query_num_true = num_true;
        }

        let (subject_composition, subject_num_true) = get_composition(
            &subject,
            &window.subject_range,
            window_align,
            crate::matrix::AA_SIZE,
            false,
            params.subject_is_translated,
        );
        let matrix_adjust_rule = match blast_adjust_position_based_scores(
            adjusted_pssm,
            &query_composition,
            query_num_true,
            &subject_composition,
            subject_num_true,
            &params.matrix_info,
            composition_test_index,
            pvalue_for_this_pair,
            lambda_ratio,
        ) {
            Ok(rule) => rule,
            Err(status) if status < 0 => {
                alignments_free_array(alignments, num_queries);
                return status;
            }
            Err(_) => continue,
        };

        let searchsp = query_info[query_index].eff_search_space;
        let mut forbidden = crate::smith_waterman::BlastForbiddenRanges::new(query.length);
        let query_offset = window.query_range.begin.max(0) as usize;
        if hspcnt <= 0 {
            continue;
        }

        loop {
            let (sw_score, match_end, query_end) =
                crate::smith_waterman::blast_smith_waterman_score_only_pssm_with_forbidden(
                    subject.data(),
                    query.data(),
                    adjusted_pssm,
                    query_offset,
                    params.gapping_params.gap_open,
                    params.gapping_params.gap_extend,
                    &forbidden,
                );
            let significant = smith_waterman_alignment_is_significant(
                sw_score,
                lambda,
                log_k,
                searchsp,
                params,
                alignments[query_index].as_deref(),
                &significant_matches[query_index],
                matching_seq.index,
            );
            if !significant {
                break;
            }

            let (_updated_score, match_start, query_start) =
                crate::smith_waterman::blast_smith_waterman_find_start_pssm_with_forbidden(
                    subject.data(),
                    query.data(),
                    adjusted_pssm,
                    query_offset,
                    params.gapping_params.gap_open,
                    params.gapping_params.gap_extend,
                    match_end,
                    query_end,
                    sw_score,
                    &forbidden,
                );

            let Some((mut new_align, query_end_xdrop, match_end_xdrop)) =
                new_alignment_using_xdrop_protein_pssm(
                    &query,
                    &window.query_range,
                    &subject,
                    &window.subject_range,
                    query_start,
                    query_end,
                    match_start,
                    match_end,
                    sw_score,
                    &params.gapping_params,
                    matrix_adjust_rule,
                    adjusted_pssm,
                )
            else {
                alignments_free_array(alignments, num_queries);
                return -1;
            };
            new_align.next = alignments[query_index].take();
            alignments[query_index] = Some(Box::new(new_align));

            if window.hspcnt > 1 {
                forbidden.push(
                    query_start as i32,
                    query_end_xdrop + 1,
                    match_start as i32,
                    match_end_xdrop,
                );
            }
            if !(significant && window.hspcnt > 1) {
                break;
            }
        }
    }

    0
}

/// blast-rs: Callback-driven no-composition port boundary for
/// `Blast_RedoOneMatchSmithWaterman` (`redo_alignment.c:1342`); not a complete direct NCBI C port.
///
/// The C routine obtains sequence ranges through `callbacks->get_range`, runs
/// Smith-Waterman, applies the significance gate, finds the start, finalizes
/// with `callbacks->new_xdrop_align`, prepends the result, and pushes accepted
/// rectangles into `Blast_ForbiddenRanges` for multi-HSP windows. This Rust
/// boundary performs that flow when no composition score adjustment is needed.
#[allow(clippy::too_many_arguments)]
pub fn blast_redo_one_match_smith_waterman_with_callbacks(
    alignments: &mut [Option<Box<BlastCompoAlignment>>],
    params: &BlastRedoAlignParams,
    incoming_aligns: &Option<Box<BlastCompoAlignment>>,
    hspcnt: i32,
    lambda: f64,
    log_k: f64,
    matching_seq: &BlastCompoMatchingSequence,
    query_info: &[BlastCompoQueryInfo],
    significant_matches: &[BlastCompoHeap],
    sw_matrix: &[[i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE],
) -> i32 {
    let mut matrix = *sw_matrix;
    let mut pvalue_for_this_pair = -1.0;
    let mut lambda_ratio = 1.0;
    blast_redo_one_match_smith_waterman_with_callbacks_and_adjustment(
        alignments,
        params,
        incoming_aligns,
        hspcnt,
        lambda,
        log_k,
        matching_seq,
        query_info,
        significant_matches,
        &mut matrix,
        None,
        &mut pvalue_for_this_pair,
        0,
        &mut lambda_ratio,
    )
}

/// blast-rs: Callback-driven port boundary for `Blast_RedoOneMatchSmithWaterman`
/// (`redo_alignment.c:1342`) with the non-position composition-adjusted
/// score path wired before the SW loop; not a complete direct NCBI C port.
#[allow(clippy::too_many_arguments)]
pub fn blast_redo_one_match_smith_waterman_with_callbacks_and_adjustment(
    alignments: &mut [Option<Box<BlastCompoAlignment>>],
    params: &BlastRedoAlignParams,
    incoming_aligns: &Option<Box<BlastCompoAlignment>>,
    hspcnt: i32,
    lambda: f64,
    log_k: f64,
    matching_seq: &BlastCompoMatchingSequence,
    query_info: &[BlastCompoQueryInfo],
    significant_matches: &[BlastCompoHeap],
    matrix: &mut [[i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE],
    workspace: Option<&BlastCompositionWorkspace>,
    pvalue_for_this_pair: &mut f64,
    composition_test_index: i32,
    lambda_ratio: &mut f64,
) -> i32 {
    let num_queries = query_info.len();
    if alignments.len() < num_queries || significant_matches.len() < num_queries {
        return -1;
    }
    alignments_free_array(alignments, num_queries);

    let Some(callbacks) = params.callbacks else {
        return -1;
    };
    if params.position_based {
        let mut adjusted_pssm = Vec::new();
        return blast_redo_one_match_smith_waterman_with_callbacks_and_position_adjustment(
            alignments,
            params,
            incoming_aligns,
            hspcnt,
            lambda,
            log_k,
            matching_seq,
            query_info,
            significant_matches,
            &mut adjusted_pssm,
            pvalue_for_this_pair,
            composition_test_index,
            lambda_ratio,
        );
    }
    let (Some(get_range), Some(new_xdrop_align)) = (callbacks.get_range, callbacks.new_xdrop_align)
    else {
        return -1;
    };
    if params.query_is_translated && params.subject_is_translated {
        return -1;
    }

    let windows = match windows_from_aligns(
        incoming_aligns,
        query_info,
        hspcnt,
        num_queries,
        K_WINDOW_BORDER,
        matching_seq.length,
        params.query_is_translated,
        params.subject_is_translated,
        params.position_based,
    ) {
        Ok(windows) => windows,
        Err(_) => return -1,
    };

    for window in &windows {
        let query_index = window.query_range.context;
        if query_index < 0 || query_index as usize >= num_queries {
            alignments_free_array(alignments, num_queries);
            return -1;
        }
        let query_index = query_index as usize;
        let Some(window_align) = window.align.as_deref() else {
            alignments_free_array(alignments, num_queries);
            return -1;
        };

        let near_identical_status = preliminary_test_near_identical(
            query_info,
            window,
            window_align,
            params.near_identical_cutoff,
        );
        let mut subject = BlastCompoSequenceData::default();
        let mut query = BlastCompoSequenceData::default();
        let mut subject_maybe_biased = true;
        let status = get_range(
            matching_seq,
            &window.subject_range,
            &mut subject,
            &query_info[query_index].seq,
            &window.query_range,
            &mut query,
            &query_info[query_index].words,
            window_align,
            near_identical_status,
            params.compo_adjust_mode,
            true,
            &mut subject_maybe_biased,
        );
        if status != 0 {
            alignments_free_array(alignments, num_queries);
            return status;
        }

        let mut query_composition = query_info[query_index].composition.prob.clone();
        let mut query_num_true = if query_composition.iter().any(|&p| p > 0.0) {
            query_info[query_index].seq.length.max(0) as usize
        } else {
            0
        };
        if params.query_is_translated {
            let (composition, num_true) = get_composition(
                &query,
                &window.query_range,
                window_align,
                crate::matrix::AA_SIZE,
                true,
                false,
            );
            query_composition = composition;
            query_num_true = num_true;
        }

        let mut matrix_adjust_rule = MatrixAdjustRule::DontAdjust;
        let mut adjust_search_failed = 0;
        if !matches!(
            params.compo_adjust_mode,
            CompoAdjustMode::NoCompositionBasedStats
        ) {
            let (subject_composition, subject_num_true) = get_composition(
                &subject,
                &window.subject_range,
                window_align,
                crate::matrix::AA_SIZE,
                false,
                params.subject_is_translated,
            );
            // NCBI passes full window-data lengths to `Blast_AdjustScores`.
            let query_window_len = query.data().len();
            let subject_window_len = subject.data().len();
            match blast_adjust_scores_with_workspace_v2(
                matrix,
                &query_composition,
                query_window_len,
                query_num_true,
                &subject_composition,
                subject_window_len,
                subject_num_true,
                &params.matrix_info,
                params.compo_adjust_mode,
                params.re_pseudocounts,
                workspace,
                composition_test_index,
                pvalue_for_this_pair,
                lambda_ratio,
            ) {
                Ok(rule) => {
                    matrix_adjust_rule = rule;
                }
                Err(status) if status < 0 => {
                    alignments_free_array(alignments, num_queries);
                    return status;
                }
                Err(_) => {
                    adjust_search_failed = 1;
                }
            }
        }

        if adjust_search_failed != 0 {
            sequence_data_release(&mut subject);
            sequence_data_release(&mut query);
            continue;
        }

        let searchsp = query_info[query_index].eff_search_space;
        let mut forbidden = crate::smith_waterman::BlastForbiddenRanges::new(query.length);
        if hspcnt <= 0 {
            sequence_data_release(&mut subject);
            sequence_data_release(&mut query);
            continue;
        }
        loop {
            let (sw_score, match_end, query_end) =
                crate::smith_waterman::blast_smith_waterman_score_only_with_forbidden(
                    subject.data(),
                    query.data(),
                    matrix,
                    params.gapping_params.gap_open,
                    params.gapping_params.gap_extend,
                    &forbidden,
                );
            let significant = smith_waterman_alignment_is_significant(
                sw_score,
                lambda,
                log_k,
                searchsp,
                params,
                alignments[query_index].as_deref(),
                &significant_matches[query_index],
                matching_seq.index,
            );
            if !significant {
                break;
            }

            let (_updated_score, match_start, query_start) =
                crate::smith_waterman::blast_smith_waterman_find_start_with_forbidden(
                    subject.data(),
                    query.data(),
                    matrix,
                    params.gapping_params.gap_open,
                    params.gapping_params.gap_extend,
                    match_end,
                    query_end,
                    sw_score,
                    &forbidden,
                );

            let mut new_align = None;
            let mut query_end_forbidden = query_end as i32;
            let mut match_end_forbidden = match_end as i32;
            let status = new_xdrop_align(
                &mut new_align,
                &mut query_end_forbidden,
                &mut match_end_forbidden,
                query_start as i32,
                match_start as i32,
                sw_score,
                &query,
                &window.query_range,
                params.ccat_query_length,
                &subject,
                &window.subject_range,
                matching_seq.length,
                &params.gapping_params,
                matrix_adjust_rule,
            );
            if status != 0 {
                alignments_free_array(alignments, num_queries);
                return status;
            }
            let Some(mut new_align) = new_align else {
                alignments_free_array(alignments, num_queries);
                return -1;
            };
            new_align.next = alignments[query_index].take();
            alignments[query_index] = Some(Box::new(new_align));

            if window.hspcnt > 1 {
                forbidden.push(
                    query_start as i32,
                    query_end_forbidden,
                    match_start as i32,
                    match_end_forbidden,
                );
            }
            if !(significant && window.hspcnt > 1) {
                break;
            }
        }
        sequence_data_release(&mut subject);
        sequence_data_release(&mut query);
    }

    0
}

/// blast-rs: Callback-driven position-specific/PSSM counterpart of
/// [`blast_redo_one_match_smith_waterman_with_callbacks_and_adjustment`]; not a direct NCBI C port.
#[allow(clippy::too_many_arguments)]
pub fn blast_redo_one_match_smith_waterman_with_callbacks_and_position_adjustment(
    alignments: &mut [Option<Box<BlastCompoAlignment>>],
    params: &BlastRedoAlignParams,
    incoming_aligns: &Option<Box<BlastCompoAlignment>>,
    hspcnt: i32,
    lambda: f64,
    log_k: f64,
    matching_seq: &BlastCompoMatchingSequence,
    query_info: &[BlastCompoQueryInfo],
    significant_matches: &[BlastCompoHeap],
    adjusted_pssm: &mut Vec<Vec<i32>>,
    pvalue_for_this_pair: &mut f64,
    composition_test_index: i32,
    lambda_ratio: &mut f64,
) -> i32 {
    let num_queries = query_info.len();
    if alignments.len() < num_queries || significant_matches.len() < num_queries {
        return -1;
    }
    alignments_free_array(alignments, num_queries);

    let Some(callbacks) = params.callbacks else {
        return -1;
    };
    let Some(get_range) = callbacks.get_range else {
        return -1;
    };

    if matches!(
        params.compo_adjust_mode,
        CompoAdjustMode::NoCompositionBasedStats
    ) || !params.position_based
        || !params.matrix_info.positional
        || params.query_is_translated
        || params.subject_is_translated
    {
        return -1;
    }

    let windows = match windows_from_aligns(
        incoming_aligns,
        query_info,
        hspcnt,
        num_queries,
        K_WINDOW_BORDER,
        matching_seq.length,
        params.query_is_translated,
        params.subject_is_translated,
        params.position_based,
    ) {
        Ok(windows) => windows,
        Err(_) => return -1,
    };

    for window in &windows {
        let query_index = window.query_range.context;
        if query_index < 0 || query_index as usize >= num_queries {
            alignments_free_array(alignments, num_queries);
            return -1;
        }
        let query_index = query_index as usize;
        let Some(window_align) = window.align.as_deref() else {
            alignments_free_array(alignments, num_queries);
            return -1;
        };

        let near_identical_status = preliminary_test_near_identical(
            query_info,
            window,
            window_align,
            params.near_identical_cutoff,
        );
        let mut subject = BlastCompoSequenceData::default();
        let mut query = BlastCompoSequenceData::default();
        let mut subject_maybe_biased = true;
        let status = get_range(
            matching_seq,
            &window.subject_range,
            &mut subject,
            &query_info[query_index].seq,
            &window.query_range,
            &mut query,
            &query_info[query_index].words,
            window_align,
            near_identical_status,
            params.compo_adjust_mode,
            true,
            &mut subject_maybe_biased,
        );
        if status != 0 {
            alignments_free_array(alignments, num_queries);
            return status;
        }

        let mut query_composition = query_info[query_index].composition.prob.clone();
        let mut query_num_true = if query_composition.iter().any(|&p| p > 0.0) {
            query_info[query_index].seq.length.max(0) as usize
        } else {
            0
        };
        if query_num_true == 0 {
            let (composition, num_true) = get_composition(
                &query,
                &window.query_range,
                window_align,
                crate::matrix::AA_SIZE,
                true,
                false,
            );
            query_composition = composition;
            query_num_true = num_true;
        }

        let (subject_composition, subject_num_true) = get_composition(
            &subject,
            &window.subject_range,
            window_align,
            crate::matrix::AA_SIZE,
            false,
            false,
        );
        let matrix_adjust_rule = match blast_adjust_position_based_scores(
            adjusted_pssm,
            &query_composition,
            query_num_true,
            &subject_composition,
            subject_num_true,
            &params.matrix_info,
            composition_test_index,
            pvalue_for_this_pair,
            lambda_ratio,
        ) {
            Ok(rule) => rule,
            Err(status) if status < 0 => {
                alignments_free_array(alignments, num_queries);
                return status;
            }
            Err(_) => {
                sequence_data_release(&mut subject);
                sequence_data_release(&mut query);
                continue;
            }
        };

        let searchsp = query_info[query_index].eff_search_space;
        let mut forbidden = crate::smith_waterman::BlastForbiddenRanges::new(query.length);
        let query_offset = window.query_range.begin.max(0) as usize;
        if hspcnt <= 0 {
            sequence_data_release(&mut subject);
            sequence_data_release(&mut query);
            continue;
        }

        loop {
            let (sw_score, match_end, query_end) =
                crate::smith_waterman::blast_smith_waterman_score_only_pssm_with_forbidden(
                    subject.data(),
                    query.data(),
                    adjusted_pssm,
                    query_offset,
                    params.gapping_params.gap_open,
                    params.gapping_params.gap_extend,
                    &forbidden,
                );
            let significant = smith_waterman_alignment_is_significant(
                sw_score,
                lambda,
                log_k,
                searchsp,
                params,
                alignments[query_index].as_deref(),
                &significant_matches[query_index],
                matching_seq.index,
            );
            if !significant {
                break;
            }

            let (_updated_score, match_start, query_start) =
                crate::smith_waterman::blast_smith_waterman_find_start_pssm_with_forbidden(
                    subject.data(),
                    query.data(),
                    adjusted_pssm,
                    query_offset,
                    params.gapping_params.gap_open,
                    params.gapping_params.gap_extend,
                    match_end,
                    query_end,
                    sw_score,
                    &forbidden,
                );

            let Some((mut new_align, query_end_xdrop, match_end_xdrop)) =
                new_alignment_using_xdrop_protein_pssm(
                    &query,
                    &window.query_range,
                    &subject,
                    &window.subject_range,
                    query_start,
                    query_end,
                    match_start,
                    match_end,
                    sw_score,
                    &params.gapping_params,
                    matrix_adjust_rule,
                    adjusted_pssm,
                )
            else {
                alignments_free_array(alignments, num_queries);
                return -1;
            };
            new_align.next = alignments[query_index].take();
            alignments[query_index] = Some(Box::new(new_align));

            if window.hspcnt > 1 {
                forbidden.push(
                    query_start as i32,
                    query_end_xdrop + 1,
                    match_start as i32,
                    match_end_xdrop,
                );
            }
            if !(significant && window.hspcnt > 1) {
                break;
            }
        }
        sequence_data_release(&mut subject);
        sequence_data_release(&mut query);
    }

    0
}

/// blast-rs: Converts modeled square matrix vectors into fixed amino-acid arrays; not a direct NCBI C port.
fn square_matrix_from_vec(
    matrix: &[Vec<i32>],
) -> Option<[[i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE]> {
    if matrix.len() < crate::matrix::AA_SIZE {
        return None;
    }
    let mut out = [[0i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE];
    for i in 0..crate::matrix::AA_SIZE {
        if matrix[i].len() < crate::matrix::AA_SIZE {
            return None;
        }
        out[i].copy_from_slice(&matrix[i][..crate::matrix::AA_SIZE]);
    }
    Some(out)
}

/// blast-rs: Converts modeled frequency-ratio rows into fixed amino-acid arrays; not a direct NCBI C port.
fn freq_ratios_from_vec(
    matrix: &[Vec<f64>],
) -> Option<[[f64; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE]> {
    if matrix.len() < crate::matrix::AA_SIZE {
        return None;
    }
    let mut out = [[0.0f64; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE];
    for i in 0..crate::matrix::AA_SIZE {
        if matrix[i].len() < crate::matrix::AA_SIZE {
            return None;
        }
        out[i].copy_from_slice(&matrix[i][..crate::matrix::AA_SIZE]);
    }
    Some(out)
}

/// blast-rs: Composition p-value helper for local `Blast_AdjustScores` adapters; not a direct NCBI C port.
fn composition_test_pvalue(query_prob: &[f64], subject_prob: &[f64]) -> f64 {
    let mut permuted_query = [0.0f64; crate::composition::COMPO_NUM_TRUE_AA];
    let mut permuted_subject = [0.0f64; crate::composition::COMPO_NUM_TRUE_AA];
    crate::compo_mode_condition::s_gather_letter_probs(query_prob, &mut permuted_query);
    crate::compo_mode_condition::s_gather_letter_probs(subject_prob, &mut permuted_subject);

    let (mut lambda_for_pair, iter_count) = crate::composition::blast_calc_lambda_full_precision(
        &crate::composition::BLOS62,
        &permuted_query,
        &permuted_subject,
        crate::composition::COMPO_NUM_TRUE_AA,
    );
    if iter_count >= crate::composition::LAMBDA_ITERATION_LIMIT {
        lambda_for_pair = crate::composition::COMPO_MIN_LAMBDA;
    }
    crate::composition::blast_composition_pvalue(lambda_for_pair)
}

/// blast-rs: PSSM score-probability helper for local position-based scaling; not a direct NCBI C port.
fn pssm_score_probs(
    matrix: &[Vec<i32>],
    rows: usize,
    subject_prob: &[f64],
) -> Option<(Vec<f64>, i32, i32)> {
    if rows == 0 || matrix.len() < rows {
        return None;
    }

    let mut obs_min = 0i32;
    let mut obs_max = 0i32;
    for row in matrix.iter().take(rows) {
        if row.len() < crate::matrix::AA_SIZE {
            return None;
        }
        for &jcol in &crate::composition::TRUE_CHAR_POSITIONS {
            let score = row[jcol];
            if score < obs_min && score > -100_000 {
                obs_min = score;
            }
            if score > obs_max {
                obs_max = score;
            }
        }
    }

    let mut score_probs = vec![0.0f64; (obs_max - obs_min + 1) as usize];
    let one_pos_frac = 1.0 / rows as f64;
    for row in matrix.iter().take(rows) {
        for &jcol in &crate::composition::TRUE_CHAR_POSITIONS {
            if jcol >= subject_prob.len() {
                return None;
            }
            let score = row[jcol];
            if score >= obs_min {
                score_probs[(score - obs_min) as usize] += one_pos_frac * subject_prob[jcol];
            }
        }
    }
    Some((score_probs, obs_min, obs_max))
}

/// blast-rs: PSSM lambda-ratio helper for local position-based scaling; not a direct NCBI C port.
fn pssm_lambda_ratio(
    matrix: &[Vec<i32>],
    rows: usize,
    subject_prob: &[f64],
    ungapped_lambda: f64,
    p_value_adjustment: bool,
) -> Option<f64> {
    let (score_probs, obs_min, obs_max) = pssm_score_probs(matrix, rows, subject_prob)?;
    let range = (obs_max - obs_min + 1) as usize;
    let mut avg = 0.0f64;
    for i in 0..range {
        avg += (obs_min + i as i32) as f64 * score_probs[i];
    }
    let adjusted_lambda = if avg >= 0.0 {
        -1.0
    } else {
        crate::composition::blast_karlin_lambda_nr(&score_probs, obs_min, obs_max, ungapped_lambda)
    };

    let mut ratio = adjusted_lambda / ungapped_lambda;
    if !p_value_adjustment {
        ratio = ratio.min(1.0);
    }
    Some(ratio.max(0.5))
}

/// blast-rs: PSSM X-residue score helper for local position-based scaling; not a direct NCBI C port.
fn pssm_x_score(row: &[f64], cols: usize, probs: &[f64]) -> f64 {
    let mut score = 0.0f64;
    for j in 0..cols.min(crate::composition::ALPHA_CONVERT.len()) {
        if crate::composition::ALPHA_CONVERT[j] >= 0 && j < probs.len() {
            score += row[j] * probs[j];
        }
    }
    score.min(-1.0)
}

/// blast-rs: Position-specific composition scaling helper; not a direct NCBI C
/// port.
///
/// This is the `s_GetPssmScoreProbs` + `s_ScalePSSM` branch used when
/// `Blast_MatrixInfo::positionBased` is true.
pub fn composition_scale_pssm_with_ratio(
    start_matrix: &[Vec<i32>],
    start_freq_ratios: &[Vec<f64>],
    rows: usize,
    cols: usize,
    subject_prob: &[f64],
    ungapped_lambda: f64,
    p_value_adjustment: bool,
) -> Option<(Vec<Vec<i32>>, f64)> {
    const E_CCHAR: usize = crate::encoding::NCBISTDAA_C as usize;
    const E_XCHAR: usize = crate::encoding::NCBISTDAA_X as usize;
    const E_SELENOCYSTEINE: usize = crate::encoding::NCBISTDAA_U as usize;
    const E_STOP_CHAR: usize = crate::encoding::NCBISTDAA_STOP as usize;
    const E_OCHAR: usize = crate::encoding::NCBISTDAA_O as usize;
    const COMPO_SCORE_MIN: f64 = -100_000.0;

    if start_matrix.len() < rows || start_freq_ratios.len() < rows || cols == 0 {
        return None;
    }
    let ratio = pssm_lambda_ratio(
        start_matrix,
        rows,
        subject_prob,
        ungapped_lambda,
        p_value_adjustment,
    )?;
    let scaled_lambda = ungapped_lambda / ratio;
    let mut matrix = vec![vec![0i32; cols]; rows];

    for p in 0..rows {
        if start_matrix[p].len() < cols || start_freq_ratios[p].len() < cols {
            return None;
        }
        let mut row = vec![0.0f64; cols];
        for (j, value) in row.iter_mut().enumerate() {
            *value = if start_freq_ratios[p][j] <= 0.0 {
                COMPO_SCORE_MIN
            } else {
                start_freq_ratios[p][j].ln() / scaled_lambda
            };
        }
        if E_XCHAR < cols {
            let x_score = pssm_x_score(&row, cols, subject_prob);
            row[E_XCHAR] = x_score;
            if E_OCHAR < cols {
                row[E_OCHAR] = x_score;
            }
        }
        if E_SELENOCYSTEINE < cols && E_CCHAR < cols {
            row[E_SELENOCYSTEINE] = row[E_CCHAR];
        }

        for j in 0..cols {
            matrix[p][j] = if row[j] < i32::MIN as f64 {
                i32::MIN
            } else {
                crate::math::blast_nint(row[j]) as i32
            };
        }
        if E_STOP_CHAR < cols {
            matrix[p][E_STOP_CHAR] = start_matrix[p][E_STOP_CHAR];
        }
    }

    Some((matrix, ratio))
}

/// blast-rs: Old-style scaling branch of `Blast_AdjustScores`
/// (`composition_adjustment.c:1446`); not a direct NCBI C port.
///
/// This covers the non-position-based `eCompositionBasedStats` path that
/// chooses `eCompoScaleOldMatrix` and calls `Blast_CompositionBasedStats`.
#[allow(clippy::too_many_arguments)]
pub fn blast_adjust_scores_scale_old_matrix(
    matrix: &mut [[i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE],
    query_prob: &[f64],
    query_num_true: usize,
    subject_prob: &[f64],
    subject_num_true: usize,
    matrix_info: &BlastMatrixInfo,
    composition_adjust_mode: CompoAdjustMode,
    composition_test_index: i32,
    pvalue_for_this_pair: &mut f64,
    ratio_to_pass_back: &mut f64,
) -> Result<MatrixAdjustRule, i32> {
    blast_adjust_scores_with_workspace(
        matrix,
        query_prob,
        query_num_true,
        subject_prob,
        subject_num_true,
        matrix_info,
        composition_adjust_mode,
        K_RE_MATRIX_ADJUSTMENT_PSEUDOCOUNTS,
        None,
        composition_test_index,
        pvalue_for_this_pair,
        ratio_to_pass_back,
    )
}

/// blast-rs: Position-based/PSSM old-style scaling branch of `Blast_AdjustScores`;
/// not a direct NCBI C port.
///
/// NCBI routes every `matrixInfo->positionBased` call to
/// `eCompoScaleOldMatrix`, then `Blast_CompositionBasedStats`, whose PSSM
/// branch computes score probabilities across PSSM rows and rescales each row
/// with `s_ScalePSSM`.
#[allow(clippy::too_many_arguments)]
pub fn blast_adjust_position_based_scores(
    matrix: &mut Vec<Vec<i32>>,
    query_prob: &[f64],
    query_num_true: usize,
    subject_prob: &[f64],
    subject_num_true: usize,
    matrix_info: &BlastMatrixInfo,
    composition_test_index: i32,
    pvalue_for_this_pair: &mut f64,
    ratio_to_pass_back: &mut f64,
) -> Result<MatrixAdjustRule, i32> {
    if query_num_true == 0 || subject_num_true == 0 {
        return Err(1);
    }
    if !matrix_info.positional {
        return Err(-1);
    }
    if composition_test_index > 0 {
        *pvalue_for_this_pair = composition_test_pvalue(query_prob, subject_prob);
    }

    let rows = matrix_info.rows.max(0) as usize;
    let cols = matrix_info.cols.max(0) as usize;
    let (scaled, ratio) = composition_scale_pssm_with_ratio(
        &matrix_info.matrix,
        &matrix_info.start_freq_ratios,
        rows,
        cols,
        subject_prob,
        matrix_info.ungapped_lambda,
        composition_test_index > 0,
    )
    .ok_or(1)?;

    *matrix = scaled;
    *ratio_to_pass_back = ratio;
    Ok(MatrixAdjustRule::ScaleOldMatrix)
}

/// NCBI: Kappa_posSearchItemsNew (blast_posit.c:43).
pub fn kappa_pos_search_items_new(
    query_length: usize,
    matrix_name: &str,
    pos_private_matrix: Vec<Vec<i32>>,
    pos_freqs: Vec<Vec<f64>>,
) -> Option<KappaPosSearchItems> {
    let freq_ratios = crate::matrix::get_matrix_freq_ratios(matrix_name)?;
    Some(KappaPosSearchItems {
        pos_matrix: vec![vec![0; crate::encoding::BLASTAA_SIZE]; query_length],
        std_freq_ratios: freq_ratios.iter().map(|row| row.to_vec()).collect(),
        query_length,
        pos_private_matrix,
        pos_freqs,
    })
}

/// NCBI: Kappa_posSearchItemsFree (blast_posit.c:86).
pub fn kappa_pos_search_items_free(
    pos_search: &mut Option<KappaPosSearchItems>,
) -> Option<KappaPosSearchItems> {
    if let Some(search) = pos_search.as_mut() {
        search.pos_matrix.clear();
        search.std_freq_ratios.clear();
        search.pos_private_matrix.clear();
        search.pos_freqs.clear();
        search.query_length = 0;
    }
    *pos_search = None;
    None
}

/// NCBI: Kappa_compactSearchItemsNew (blast_posit.c:101).
pub fn kappa_compact_search_items_new(
    query: &[u8],
    query_length: usize,
    sbp: &crate::stat::BlastScoreBlk,
) -> Option<KappaCompactSearchItems> {
    if !sbp.protein_alphabet
        || sbp.alphabet_code != crate::encoding::BLASTAA_SEQ_CODE
        || sbp.alphabet_size != crate::encoding::BLASTAA_SIZE
        || query.len() < query_length
    {
        return None;
    }
    let ideal = sbp.kbp_ideal.clone().unwrap_or_default();
    Some(KappaCompactSearchItems {
        standard_prob: crate::stat::protein_std_freq_ncbistdaa().to_vec(),
        query: query[..query_length].to_vec(),
        qlength: query_length,
        alphabet_size: crate::encoding::BLASTAA_SIZE,
        matrix: sbp.matrix.data.clone(),
        kbp_std: sbp.kbp_std.clone(),
        kbp_psi: sbp.kbp_psi.clone(),
        kbp_gap_std: sbp.kbp_gap_std.clone(),
        kbp_gap_psi: sbp.kbp_gap_psi.clone(),
        lambda_ideal: ideal.lambda,
        k_ideal: ideal.k,
    })
}

/// NCBI: Kappa_compactSearchItemsFree (blast_posit.c:141).
pub fn kappa_compact_search_items_free(
    compact_search: &mut Option<KappaCompactSearchItems>,
) -> Option<KappaCompactSearchItems> {
    if let Some(search) = compact_search.as_mut() {
        search.standard_prob.clear();
        search.query.clear();
        search.qlength = 0;
        search.alphabet_size = 0;
        search.matrix.clear();
        search.kbp_std.clear();
        search.kbp_psi.clear();
        search.kbp_gap_std.clear();
        search.kbp_gap_psi.clear();
        search.lambda_ideal = 0.0;
        search.k_ideal = 0.0;
    }
    *compact_search = None;
    None
}

/// NCBI: fillSfp (blast_posit.c:170).
pub fn fill_sfp(
    matrix: &[Vec<i32>],
    matrix_length: usize,
    query_prob_array: &[f64],
) -> Option<crate::stat::ScoreFreq> {
    if matrix_length == 0 || matrix.len() < matrix_length {
        return None;
    }

    let mut min_score = crate::stat::BLAST_SCORE_MAX;
    let mut max_score = crate::stat::BLAST_SCORE_MIN;
    for row in matrix.iter().take(matrix_length) {
        for &k in &crate::composition::TRUE_CHAR_POSITIONS {
            let score = *row.get(k)?;
            if score != crate::stat::BLAST_SCORE_MIN && score < min_score {
                min_score = score;
            }
            if score > max_score {
                max_score = score;
            }
        }
    }
    if min_score == crate::stat::BLAST_SCORE_MAX
        || max_score == crate::stat::BLAST_SCORE_MIN
        || max_score - min_score >= K_SCORE_MATRIX_SCORE_RANGE
    {
        return None;
    }

    let score_range = (max_score - min_score + 1) as usize;
    let mut sprob = Vec::new();
    sprob.try_reserve_exact(score_range).ok()?;
    sprob.resize(score_range, 0.0);
    let mut sfp = crate::stat::ScoreFreq {
        score_min: min_score,
        score_max: max_score,
        obs_min: min_score,
        obs_max: max_score,
        score_avg: 0.0,
        sprob,
    };
    let one_pos_frac = 1.0 / matrix_length as f64;
    for row in matrix.iter().take(matrix_length) {
        for &k in &crate::composition::TRUE_CHAR_POSITIONS {
            let score = row[k];
            if score >= min_score {
                sfp.sprob[(score - min_score) as usize] +=
                    one_pos_frac * query_prob_array.get(k).copied().unwrap_or(0.0);
            }
        }
    }
    sfp.score_avg = (min_score..=max_score)
        .map(|score| score as f64 * sfp.sprob[(score - min_score) as usize])
        .sum();
    Some(sfp)
}

/// NCBI: _PSIComputeScoreProbabilities (blast_psi_priv.c:2647).
pub fn psi_compute_score_probabilities(
    pssm: &[Vec<i32>],
    query: &[u8],
    query_length: usize,
    std_probs: &[f64],
    sbp: &crate::stat::BlastScoreBlk,
) -> Option<crate::stat::ScoreFreq> {
    if sbp.alphabet_code != crate::encoding::BLASTAA_SEQ_CODE
        || pssm.len() < query_length
        || query.len() < query_length
    {
        return None;
    }
    let effective_length = query
        .iter()
        .take(query_length)
        .filter(|&&residue| residue != crate::encoding::NCBISTDAA_X)
        .count();
    if effective_length == 0 {
        return None;
    }

    let mut min_score = crate::stat::BLAST_SCORE_MAX;
    let mut max_score = crate::stat::BLAST_SCORE_MIN;
    for p in 0..query_length {
        if query[p] == crate::encoding::NCBISTDAA_X {
            continue;
        }
        for &aa in &crate::composition::TRUE_CHAR_POSITIONS {
            let score = *pssm[p].get(aa)?;
            if score <= crate::stat::BLAST_SCORE_MIN || score >= crate::stat::BLAST_SCORE_MAX {
                continue;
            }
            min_score = min_score.min(score);
            max_score = max_score.max(score);
        }
    }
    if min_score == crate::stat::BLAST_SCORE_MAX || max_score == crate::stat::BLAST_SCORE_MIN {
        return None;
    }

    let score_range = (max_score - min_score + 1) as usize;
    let mut sprob = Vec::new();
    sprob.try_reserve_exact(score_range).ok()?;
    sprob.resize(score_range, 0.0);
    let mut sfp = crate::stat::ScoreFreq {
        score_min: min_score,
        score_max: max_score,
        obs_min: min_score,
        obs_max: max_score,
        score_avg: 0.0,
        sprob,
    };
    for p in 0..query_length {
        if query[p] == crate::encoding::NCBISTDAA_X {
            continue;
        }
        for &aa in &crate::composition::TRUE_CHAR_POSITIONS {
            let score = pssm[p][aa];
            if score <= crate::stat::BLAST_SCORE_MIN || score >= crate::stat::BLAST_SCORE_MAX {
                continue;
            }
            sfp.sprob[(score - min_score) as usize] +=
                std_probs.get(aa).copied().unwrap_or(0.0) / effective_length as f64;
        }
    }
    sfp.score_avg = (min_score..=max_score)
        .map(|score| score as f64 * sfp.sprob[(score - min_score) as usize])
        .sum();
    Some(sfp)
}

/// NCBI: _PSIUpdateLambdaK (blast_psi_priv.c:2732).
pub fn psi_update_lambda_k(
    pssm: &[Vec<i32>],
    query: &[u8],
    query_length: usize,
    std_probs: &[f64],
    sbp: &mut crate::stat::BlastScoreBlk,
) -> i32 {
    let Some(score_freqs) =
        psi_compute_score_probabilities(pssm, query, query_length, std_probs, sbp)
    else {
        return 1;
    };
    if sbp.kbp_psi.is_empty() {
        sbp.kbp_psi.resize_with(1, crate::stat::KarlinBlk::default);
    }
    if crate::stat::blast_karlin_blk_ungapped_calc(Some(&mut sbp.kbp_psi[0]), Some(&score_freqs))
        != 0
    {
        return 1;
    }
    let Some(ideal) = sbp.kbp_ideal.as_ref() else {
        return 1;
    };
    if sbp.kbp_gap_std.is_empty() || ideal.k == 0.0 {
        return 1;
    }
    if sbp.kbp_gap_psi.is_empty() {
        sbp.kbp_gap_psi
            .resize_with(1, crate::stat::KarlinBlk::default);
    }
    sbp.kbp_gap_psi[0].k = sbp.kbp_psi[0].k * sbp.kbp_gap_std[0].k / ideal.k;
    sbp.kbp_gap_psi[0].log_k = sbp.kbp_gap_psi[0].k.ln();
    0
}

/// blast-rs: Local row-scaling helper factored from `impalaScaleMatrix`; not a direct NCBI C port.
fn impala_apply_factor(
    matrix: &mut [Vec<i32>],
    private_matrix: &[Vec<i32>],
    dim1: usize,
    dim2: usize,
    factor: f64,
    div_factor: f64,
) -> Option<()> {
    for c in 0..dim1 {
        for a in 0..dim2 {
            let private = *private_matrix.get(c)?.get(a)?;
            let target = matrix.get_mut(c)?.get_mut(a)?;
            *target = if private == crate::stat::BLAST_SCORE_MIN {
                crate::stat::BLAST_SCORE_MIN
            } else {
                ((factor * private as f64) / div_factor) as i32
            };
        }
    }
    Some(())
}

/// NCBI: impalaScaleMatrix (blast_posit.c:230).
pub fn impala_scale_matrix(
    compact_search: &KappaCompactSearchItems,
    pos_matrix: &mut [Vec<i32>],
    pos_private_matrix: &mut [Vec<i32>],
    scaling_factor: f64,
    do_binary_search: bool,
    sbp: &mut crate::stat::BlastScoreBlk,
) -> bool {
    if scaling_factor == 0.0 || compact_search.qlength == 0 {
        return false;
    }
    let dim1 = compact_search.qlength;
    let dim2 = compact_search.alphabet_size;
    if pos_matrix.len() < dim1 || pos_private_matrix.len() < dim1 {
        return false;
    }

    let lambda = compact_search.lambda_ideal / scaling_factor;
    let div_factor = K_PSI_SCALE_FACTOR / scaling_factor;
    let initial_lambda = compact_search
        .kbp_psi
        .first()
        .map(|kbp| kbp.lambda / scaling_factor)
        .filter(|lambda| lambda.is_finite() && *lambda > 0.0)
        .unwrap_or(lambda);
    let mut factor = 1.0;

    if do_binary_search {
        let mut factor_low = 1.0;
        let mut factor_high = 1.0;
        let mut too_high = true;
        let mut first_time = true;

        loop {
            if impala_apply_factor(
                pos_matrix,
                pos_private_matrix,
                dim1,
                dim2,
                factor,
                div_factor,
            )
            .is_none()
            {
                return false;
            }
            let Some(sfp) = fill_sfp(pos_matrix, dim1, &compact_search.standard_prob) else {
                return false;
            };
            let new_lambda = crate::stat::blast_karlin_lambda_nr(Some(&sfp), initial_lambda);

            if new_lambda > lambda {
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
            } else if first_time {
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
        }

        for _ in 0..K_POSIT_SCALING_NUM_ITERATIONS {
            factor = 0.5 * (factor_high + factor_low);
            if impala_apply_factor(
                pos_matrix,
                pos_private_matrix,
                dim1,
                dim2,
                factor,
                div_factor,
            )
            .is_none()
            {
                return false;
            }
            let Some(sfp) = fill_sfp(pos_matrix, dim1, &compact_search.standard_prob) else {
                return false;
            };
            let new_lambda = crate::stat::blast_karlin_lambda_nr(Some(&sfp), initial_lambda);
            if new_lambda > lambda {
                factor_low = factor;
            } else {
                factor_high = factor;
            }
        }
    }

    for c in 0..dim1 {
        for a in 0..dim2 {
            if pos_private_matrix[c][a] != crate::stat::BLAST_SCORE_MIN {
                pos_matrix[c][a] = crate::math::blast_nint(
                    pos_private_matrix[c][a] as f64 * factor / K_PSI_SCALE_FACTOR,
                ) as i32;
            }
        }
    }

    if let Some(sfp) = fill_sfp(pos_matrix, dim1, &compact_search.standard_prob) {
        let new_lambda = crate::stat::blast_karlin_lambda_nr(Some(&sfp), initial_lambda);
        if new_lambda.is_finite() && new_lambda > 0.0 {
            if sbp.kbp_psi.is_empty() {
                sbp.kbp_psi.resize_with(1, crate::stat::KarlinBlk::default);
            }
            sbp.kbp_psi[0].lambda = new_lambda;
        }
    }

    let scale_factor = scaling_factor / K_PSI_SCALE_FACTOR;
    for c in 0..dim1 {
        for a in 0..dim2 {
            if pos_private_matrix[c][a] != crate::stat::BLAST_SCORE_MIN {
                pos_private_matrix[c][a] = crate::math::blast_nint(
                    pos_private_matrix[c][a] as f64 * factor * scale_factor,
                ) as i32;
            }
        }
    }
    true
}

/// NCBI: Kappa_impalaScaling (blast_posit.c:393).
pub fn kappa_impala_scaling(
    pos_search: &mut KappaPosSearchItems,
    compact_search: &KappaCompactSearchItems,
    scaling_factor: f64,
    do_binary_search: bool,
    sbp: &mut crate::stat::BlastScoreBlk,
) -> i32 {
    if impala_scale_matrix(
        compact_search,
        &mut pos_search.pos_matrix,
        &mut pos_search.pos_private_matrix,
        scaling_factor,
        do_binary_search,
        sbp,
    ) {
        0
    } else {
        1
    }
}

/// blast-rs: Port boundary for `Blast_AdjustScores` (`composition_adjustment.c:1446`);
/// not a complete direct NCBI C port.
///
/// Implements the non-position-based paths whose required state is available:
/// composition-based scale-old-matrix adjustment and relative-entropy matrix
/// optimization with an explicit [`BlastCompositionWorkspace`]. If RE
/// optimization returns a positive status, NCBI falls back to old-style
/// scaling; this helper does the same. The non-position-based p-value-test
/// branch mirrors `Blast_CalcLambdaFullPrecision` + `Blast_CompositionPvalue`.
/// Position-based callers should use [`blast_adjust_position_based_scores`],
/// since PSSMs have query-length rows rather than a square amino-acid matrix.
#[allow(clippy::too_many_arguments)]
pub fn blast_adjust_scores_with_workspace(
    matrix: &mut [[i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE],
    query_prob: &[f64],
    query_num_true: usize,
    subject_prob: &[f64],
    subject_num_true: usize,
    matrix_info: &BlastMatrixInfo,
    composition_adjust_mode: CompoAdjustMode,
    re_pseudocounts: i32,
    workspace: Option<&BlastCompositionWorkspace>,
    composition_test_index: i32,
    pvalue_for_this_pair: &mut f64,
    ratio_to_pass_back: &mut f64,
) -> Result<MatrixAdjustRule, i32> {
    // NCBI `Blast_AdjustScores` (`composition_adjustment.c:1414`) takes
    // both `queryLength`/`subjectLength` (full sequence-block lengths
    // including ambiguities) AND the `numTrueAminoAcids` field of the
    // composition struct. We collapse them here by approximating
    // `query_length = max(query_num_true, num residues in prob)` since
    // the prob array already encodes the residue distribution. For most
    // callers `query_num_true` equals the window length; only matters
    // when X-heavy windows reach `blast_choose_matrix_adjust_rule`. Backward-
    // compatible default; see [`blast_adjust_scores_with_workspace_v2`]
    // for the explicit-length variant that callers should prefer.
    blast_adjust_scores_with_workspace_v2(
        matrix,
        query_prob,
        query_num_true,
        query_num_true,
        subject_prob,
        subject_num_true,
        subject_num_true,
        matrix_info,
        composition_adjust_mode,
        re_pseudocounts,
        workspace,
        composition_test_index,
        pvalue_for_this_pair,
        ratio_to_pass_back,
    )
}

/// blast-rs: `Blast_AdjustScores` adapter with explicit
/// `queryLength`/`subjectLength` separate from `numTrueAminoAcids`. The
/// naming: Version suffix distinguishes this Rust adapter from the compatibility entry point.
/// length parameters are used for [`blast_choose_matrix_adjust_rule`] (which
/// applies the high-pair/length-ratio thresholds against the full
/// window length); the `num_true` counts are used for the
/// pseudocount-weighted matrix optimization path.
#[allow(clippy::too_many_arguments)]
pub fn blast_adjust_scores_with_workspace_v2(
    matrix: &mut [[i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE],
    query_prob: &[f64],
    query_length: usize,
    query_num_true: usize,
    subject_prob: &[f64],
    subject_length: usize,
    subject_num_true: usize,
    matrix_info: &BlastMatrixInfo,
    composition_adjust_mode: CompoAdjustMode,
    re_pseudocounts: i32,
    workspace: Option<&BlastCompositionWorkspace>,
    composition_test_index: i32,
    pvalue_for_this_pair: &mut f64,
    ratio_to_pass_back: &mut f64,
) -> Result<MatrixAdjustRule, i32> {
    if query_num_true == 0 || subject_num_true == 0 {
        return Err(1);
    }

    if matrix_info.positional {
        return Err(-1);
    }

    let mut permuted_query = [0.0f64; crate::composition::COMPO_NUM_TRUE_AA];
    let mut permuted_subject = [0.0f64; crate::composition::COMPO_NUM_TRUE_AA];
    if composition_test_index > 0
        || !matches!(
            composition_adjust_mode,
            CompoAdjustMode::CompositionBasedStats
        )
    {
        crate::compo_mode_condition::s_gather_letter_probs(query_prob, &mut permuted_query);
        crate::compo_mode_condition::s_gather_letter_probs(subject_prob, &mut permuted_subject);
    }

    if composition_test_index > 0 {
        *pvalue_for_this_pair = crate::composition::blast_composition_pvalue({
            let (mut lambda_for_pair, iter_count) =
                crate::composition::blast_calc_lambda_full_precision(
                    &crate::composition::BLOS62,
                    &permuted_query,
                    &permuted_subject,
                    crate::composition::COMPO_NUM_TRUE_AA,
                );
            if iter_count >= crate::composition::LAMBDA_ITERATION_LIMIT {
                lambda_for_pair = crate::composition::COMPO_MIN_LAMBDA;
            }
            lambda_for_pair
        });
    }

    // NCBI `composition_adjustment.c:1494` passes `queryLength` and
    // `subjectLength` (the full sequence-block lengths) to
    // `Blast_ChooseMatrixAdjustRule`, NOT `numTrueAminoAcids`. The
    // length-ratio and high-pair-frequency tests inside
    // `s_TestToApplyREAdjustmentConditional` use the lengths as bounds
    // (`length <= LENGTH_LOWER_THRESHOLD = 50` short-circuit) and as
    // ratio terms. Using `num_true` here was a subtle divergence —
    // shortens the effective length for X-heavy windows and biases the
    // decision toward `UserSpecifiedRelEntropy`.
    let mut matrix_adjust_rule = if matches!(
        composition_adjust_mode,
        CompoAdjustMode::CompositionBasedStats
    ) {
        MatrixAdjustRule::ScaleOldMatrix
    } else {
        crate::compo_mode_condition::blast_choose_matrix_adjust_rule(
            query_length,
            subject_length,
            &permuted_query,
            &permuted_subject,
            composition_adjust_mode as u8,
        )
    };

    if matrix_adjust_rule != MatrixAdjustRule::ScaleOldMatrix {
        let workspace = workspace.ok_or(-1)?;
        let start_matrix = square_matrix_from_vec(&matrix_info.matrix).ok_or(-1)?;
        let status = crate::composition::blast_composition_matrix_adj(
            matrix,
            matrix_info.cols.max(0) as usize,
            matrix_adjust_rule,
            query_num_true,
            subject_num_true,
            query_prob,
            subject_prob,
            re_pseudocounts,
            K_FIXED_RE_BLOSUM62,
            &workspace.mat_b,
            &workspace.first_standard_freq,
            &workspace.second_standard_freq,
            matrix_info.ungapped_lambda,
            &start_matrix,
        );
        *ratio_to_pass_back = 1.0;
        if status == 0 {
            return Ok(matrix_adjust_rule);
        }
        if status < 0 {
            return Err(status);
        }
        matrix_adjust_rule = MatrixAdjustRule::ScaleOldMatrix;
    }

    let start_matrix = square_matrix_from_vec(&matrix_info.matrix).ok_or(-1)?;
    let freq_ratios = freq_ratios_from_vec(&matrix_info.start_freq_ratios).ok_or(-1)?;
    let (scaled, ratio) = if composition_test_index > 0 {
        crate::composition::composition_scale_matrix_with_ratio_and_adjustment(
            &start_matrix,
            query_prob,
            subject_prob,
            matrix_info.ungapped_lambda,
            &freq_ratios,
            true,
        )
    } else {
        crate::composition::composition_scale_matrix_with_ratio(
            &start_matrix,
            query_prob,
            subject_prob,
            matrix_info.ungapped_lambda,
            &freq_ratios,
        )
    }
    .ok_or(1)?;
    *matrix = scaled;
    *ratio_to_pass_back = ratio;
    Ok(matrix_adjust_rule)
}

/// `NEAR_IDENTICAL_BITS_PER_POSITION` (`blast_kappa.c:2399`).
/// Per-position bit-score threshold used to decide whether two
/// sequences are "near identical" before running the more
/// expensive composition-adjustment pass.
pub const NEAR_IDENTICAL_BITS_PER_POSITION: f64 = 1.74;

/// NCBI: s_GetAlignParams (blast_kappa.c:2406).
///
/// Assembly point for the redo-alignment driver: builds a
/// [`BlastRedoAlignParams`] from the search context. Composes the
/// already-ported helpers ([`s_matrix_info_init`],
/// [`s_gapping_params_new`], [`blast_redo_align_params_new`]) in
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
///   [`s_gapping_params_new`]).
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
pub fn s_get_align_params(
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
    let subject_is_translated_p = crate::program::blast_subject_is_translated(program);
    let query_is_translated_p = crate::program::blast_query_is_translated(program);

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
    let status = s_matrix_info_init(
        &mut matrix_info,
        matrix_name,
        kbp_ideal_lambda,
        local_scaling_factor,
    );
    if status != 0 {
        return None;
    }
    // `s_matrix_info_init` fills the shared non-position matrix fields
    // and defaults this flag to false; preserve the caller's PSI/PSSM branch
    // flag at the assembly point, matching `sbp->psi_matrix != NULL`.
    matrix_info.positional = position_based;

    // C: `gapping_params = s_GappingParamsNew(context, extendParams,
    //                                          last_context + 1);`
    let gapping_params = s_gapping_params_new(
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
        local_scaling_factor,
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

/// `kFixedReBlosum62` (`composition_adjustment.c:77`).
pub const K_FIXED_RE_BLOSUM62: f64 = 0.44;

/// NCBI: Blast_RedoAlignParamsNew (redo_alignment.c:1013).
///
/// Bundles the matrix info + gapping params + composition-adjust mode
/// flags into `BlastRedoAlignParams` driver state. NCBI's C
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
    local_scaling_factor: f64,
    position_based: bool,
    query_is_translated: bool,
    subject_is_translated: bool,
    ccat_query_length: i32,
    cutoff_s: i32,
    cutoff_e: f64,
    do_link_hsps: bool,
    near_identical_cutoff: f64,
) -> BlastRedoAlignParams {
    blast_redo_align_params_new_with_callbacks(
        matrix_info,
        gapping_params,
        compo_adjust_mode,
        local_scaling_factor,
        position_based,
        query_is_translated,
        subject_is_translated,
        ccat_query_length,
        cutoff_s,
        cutoff_e,
        do_link_hsps,
        None,
        near_identical_cutoff,
    )
}

/// blast-rs: Variant of [`blast_redo_align_params_new`] that carries the
/// `Blast_RedoAlignCallbacks *callbacks` argument from NCBI's constructor; not a direct NCBI C port.
#[allow(clippy::too_many_arguments)]
pub fn blast_redo_align_params_new_with_callbacks(
    matrix_info: BlastMatrixInfo,
    gapping_params: BlastCompoGappingParams,
    compo_adjust_mode: CompoAdjustMode,
    local_scaling_factor: f64,
    position_based: bool,
    query_is_translated: bool,
    subject_is_translated: bool,
    ccat_query_length: i32,
    cutoff_s: i32,
    cutoff_e: f64,
    do_link_hsps: bool,
    callbacks: Option<BlastRedoAlignCallbacks>,
    near_identical_cutoff: f64,
) -> BlastRedoAlignParams {
    BlastRedoAlignParams {
        matrix_info,
        gapping_params,
        compo_adjust_mode,
        local_scaling_factor,
        position_based,
        re_pseudocounts: K_RE_MATRIX_ADJUSTMENT_PSEUDOCOUNTS,
        subject_is_translated,
        query_is_translated,
        ccat_query_length,
        cutoff_s,
        cutoff_e,
        do_link_hsps,
        callbacks,
        near_identical_cutoff,
    }
}

/// NCBI: Blast_RedoAlignParamsFree (redo_alignment.c:1001).
///
/// In Rust, `Drop` handles the deallocation that NCBI does manually.
/// Provided for parity so callers porting NCBI flow can keep matching
/// call sites; sets the slot to `None` to mirror C's
/// `*pparams = NULL`.
pub fn blast_redo_align_params_free(slot: &mut Option<BlastRedoAlignParams>) {
    *slot = None;
}

/// NCBI: s_ResultHspToDistinctAlign (blast_kappa.c:769).
/// naming: Public Rust helper omits the private `s_` prefix.
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
            (hsp.score as f64 * local_scaling_factor) as i32,
            MatrixAdjustRule::DontAdjust,
            hsp.context,
            hsp.query_offset,
            hsp.query_end,
            hsp.subject_offset,
            hsp.subject_end,
            hsp.subject_frame,
            None,
        )
        .with_hsp_context(hsp.query_gapped_start, hsp.subject_gapped_start);
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

/// NCBI: s_DoSegSequenceData (blast_kappa.c:1427).
/// naming: Public Rust helper omits the private `s_` prefix.
///
/// Filters low-complexity regions out of `seq_data` using the SEG
/// algorithm (NCBIstdaa input). NCBI's C uses
/// `BlastFilteringOptionsFromString(BLASTP_MASK_INSTRUCTIONS)`
/// followed by `BlastSetUp_Filter` and `Blast_MaskTheResidues`. Our
/// Rust port collapses those into the existing
/// `crate::filter::seg_filter_ncbistdaa` + an in-place X-substitution,
/// which is what `Blast_MaskTheResidues` ultimately does for protein
/// (replace each masked residue with NCBIstdaa `X`).
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
    for r in &mask.regions {
        let r_start = (r.start.max(0) as usize).min(view.len());
        let r_end = (r.end as usize).min(view.len());
        for aa in &mut view[r_start..r_end] {
            *aa = crate::encoding::NCBISTDAA_X;
        }
    }
    (0, was_biased)
}

/// NCBI: s_MatchingSequenceRelease (blast_kappa.c:907).
/// naming: Public Rust helper omits the private `s_` prefix.
///
/// In the C source this releases the BlastSeqSrc-managed subject
/// sequence and frees the per-match `BlastKappa_SequenceInfo`
/// payload. Rust now carries an owned `local_data` slot corresponding to C's
/// `void *`; the pointed-to subject-source state is still managed by the
/// caller's lifetime. We provide this hook for parity so
/// callers porting NCBI flow can keep matching call sites. The C function
/// releases through `BlastSeqSrc` only when `index >= 0 && length > 0`, then
/// nulls `local_data`; it does not reset `index` or `length`.
pub fn matching_sequence_release(self_: &mut BlastCompoMatchingSequence) {
    if let Some(mut local_data) = self_.local_data.take() {
        if self_.index >= 0 && self_.length > 0 {
            if let Some(seq_src) = local_data.seq_src {
                // C stores a borrowed BlastSeqSrc* in local_data and releases through
                // that vtable before freeing the payload.
                unsafe {
                    crate::seqsrc::blast_seq_src_release_sequence(
                        Some(&*seq_src),
                        Some(&mut local_data.seq_arg),
                    );
                }
            }
        }
    }
}

/// NCBI: s_RescaleSearch (blast_kappa.c:2117).
/// naming: Public Rust helper omits the private `s_` prefix.
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
        // NCBI gates only on `sbp->kbp_gap[i] != NULL` (`blast_kappa.c:2124`)
        // — there is no `Lambda > 0` guard before the divide/log. The Rust
        // analogue of the NULL check is simply "slot exists in the array".
        if let Some(kbp) = kbp_gap.get_mut(i) {
            kbp.lambda /= scale_factor;
            // C: `kbp->logK = log(kbp->K);` — recomputes from the
            // (unchanged) K to flush any prior cached value.
            kbp.log_k = kbp.k.ln();
        }
    }
    scoring.gap_open = crate::math::blast_nint(scoring.gap_open as f64 * scale_factor) as i32;
    scoring.gap_extend = crate::math::blast_nint(scoring.gap_extend as f64 * scale_factor) as i32;
    scoring.scale_factor = scale_factor;
}

/// blast-rs: Query-preparation preamble of `s_SequenceGetRange`
/// (`blast_kappa.c:1670`); not a complete direct NCBI C port.
///
/// Copies `query[q_range.begin .. q_range.end]` into a freshly-allocated
/// buffer with NCBI's sentinel-byte layout (`buffer[0] = 0`, real data
/// at `buffer[1 ..= length]`, terminator at `buffer[length + 1] = 0`),
/// substituting selenocysteine with cysteine per the C source's inline rule.
///
/// The full `s_SequenceGetRange` then dispatches to
/// `s_SequenceGetProteinRange` or `s_SequenceGetTranslatedRange` to
/// fill in the subject `seqData`; those branches are modeled by the
/// in-memory helpers below.
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
        *slot = if src != crate::encoding::NCBISTDAA_U {
            src
        } else {
            crate::encoding::NCBISTDAA_C
        };
    }
    BlastCompoSequenceData {
        buffer,
        data_offset: 1,
        length: length as i32,
    }
}

/// blast-rs: In-memory protein branch of `s_SequenceGetProteinRange`
/// (`blast_kappa.c:1586`); not a complete direct NCBI C port.
///
/// NCBI obtains the source residues from `BlastSeqSrcGetSequence`, then
/// creates the same sentinel-byte layout used by the query path:
/// leading zero, copied residues, trailing zero. This function covers
/// that byte-level work once the caller has already materialized the
/// subject protein sequence as NCBIstdaa bytes; the SeqSrc branch calls this
/// after fetching the program-specific subject encoding.
pub fn sequence_get_protein_range(
    source: &[u8],
    range: &BlastCompoSequenceRange,
) -> BlastCompoSequenceData {
    let begin = range.begin.max(0) as usize;
    let end = range.end.max(range.begin) as usize;
    let length = end.saturating_sub(begin);
    let mut buffer = vec![0u8; length + 2];
    if begin < source.len() {
        let copy_len = length.min(source.len() - begin);
        buffer[1..1 + copy_len].copy_from_slice(&source[begin..begin + copy_len]);
    }
    BlastCompoSequenceData {
        buffer,
        data_offset: 1,
        length: length as i32,
    }
}

/// In-memory counterpart of NCBI `s_SequenceGetProteinRange` used by
/// `Blast_RedoOneMatch`: copy the full subject, optionally SEG-mask it using
/// the C near-identical gate, then fit the data pointer to `range`.
#[allow(clippy::too_many_arguments)]
pub fn sequence_get_protein_range_for_redo(
    source: &[u8],
    range: &BlastCompoSequenceRange,
    query_range: &BlastCompoSequenceRange,
    query_data: &BlastCompoSequenceData,
    query_words: &[u64],
    align: &BlastCompoAlignment,
    should_test_identical: bool,
    compo_adjust_mode: CompoAdjustMode,
    subject_maybe_biased: &mut bool,
) -> Result<BlastCompoSequenceData, &'static str> {
    let mut seq_data = BlastCompoSequenceData {
        buffer: vec![0u8; source.len() + 2],
        data_offset: 1,
        length: source.len() as i32,
    };
    seq_data.buffer[1..1 + source.len()].copy_from_slice(source);

    if !matches!(compo_adjust_mode, CompoAdjustMode::NoCompositionBasedStats)
        && *subject_maybe_biased
        && (!should_test_identical
            || !test_near_identical(
                &seq_data,
                0,
                query_data,
                query_range.begin,
                query_words,
                align,
            ))
    {
        let (status, was_biased) = do_seg_sequence_data(&mut seq_data);
        if status != 0 {
            return Err("SEG filtering failed");
        }
        *subject_maybe_biased = was_biased;
    }

    let begin = range.begin.max(0) as usize;
    let end = range.end.max(range.begin) as usize;
    if end > source.len() {
        return Err("protein range outside subject");
    }
    seq_data.buffer[begin] = 0;
    seq_data.data_offset = begin + 1;
    seq_data.length = end.saturating_sub(begin) as i32;
    Ok(seq_data)
}

/// blast-rs: In-memory port of `s_SequenceGetTranslatedRange`
/// (`blast_kappa.c:1637`); not a complete direct NCBI C port.
///
/// The C path obtains a `BlastTargetTranslation` for the subject/frame and
/// copies the requested protein-space interval into sentinel-backed
/// `BlastCompo_SequenceData`. This helper covers the same byte-level behavior
/// when the caller has already materialized the translated subject source as
/// NCBI4na bytes.
///
/// `range.context` is the translated subject frame (`+1..+3` or `-1..-3`),
/// matching the `BlastCompo_Alignment.frame` value later produced by
/// `s_NewAlignmentFromGapAlign`.
pub fn sequence_get_translated_range(
    source_ncbi4na: &[u8],
    range: &BlastCompoSequenceRange,
    genetic_code: &[u8; 64],
) -> Result<BlastCompoSequenceData, &'static str> {
    s_sequence_get_translated_range(source_ncbi4na, range, genetic_code)
}

/// blast-rs: Name-matched wrapper for NCBI static `s_SequenceGetTranslatedRange`
/// (`blast_kappa.c:1475`) over a materialized NCBI4na subject; not a direct NCBI C port.
pub fn s_sequence_get_translated_range(
    source_ncbi4na: &[u8],
    range: &BlastCompoSequenceRange,
    genetic_code: &[u8; 64],
) -> Result<BlastCompoSequenceData, &'static str> {
    let frame = range.context;
    if !matches!(frame, 1..=3 | -3..=-1) {
        return Err("invalid translated subject frame");
    }
    let begin = range.begin.max(0);
    let end = range.end.max(range.begin);
    let translated_len = end.saturating_sub(begin);
    let translation_start = if frame > 0 {
        3 * begin
    } else {
        source_ncbi4na.len() as i32 - 3 * end + frame + 1
    };
    let num_nucleotides = 3 * translated_len + frame.abs() - 1;
    if translation_start < 0 || num_nucleotides < 0 {
        return Err("translated range outside subject");
    }
    let start = translation_start as usize;
    let nuc_len = num_nucleotides as usize;
    if start.saturating_add(nuc_len) > source_ncbi4na.len() {
        return Err("translated range outside subject");
    }

    let source_window = &source_ncbi4na[start..start + nuc_len];
    let reverse = if frame < 0 {
        crate::util::get_reverse_nucl_sequence(source_window, nuc_len)
    } else {
        Vec::new()
    };
    let mut buffer = vec![0u8; nuc_len / 3 + 2];
    let length = crate::util::blast_get_translation(
        source_window,
        &reverse,
        nuc_len,
        frame,
        &mut buffer,
        genetic_code,
    );
    Ok(BlastCompoSequenceData {
        buffer,
        data_offset: 1,
        length: length as i32,
    })
}

/// In-memory counterpart of NCBI `s_SequenceGetTranslatedRange` used by
/// `Blast_RedoOneMatch`: translate the requested subject window, then apply
/// the same composition/near-identical SEG gate as the C helper.
#[allow(clippy::too_many_arguments)]
pub fn sequence_get_translated_range_for_redo(
    source_ncbi4na: &[u8],
    range: &BlastCompoSequenceRange,
    genetic_code: &[u8; 64],
    query_range: &BlastCompoSequenceRange,
    query_data: &BlastCompoSequenceData,
    query_words: &[u64],
    align: &BlastCompoAlignment,
    should_test_identical: bool,
    compo_adjust_mode: CompoAdjustMode,
    subject_maybe_biased: &mut bool,
) -> Result<BlastCompoSequenceData, &'static str> {
    let mut seq_data = s_sequence_get_translated_range(source_ncbi4na, range, genetic_code)?;
    if !KAPPA_TBLASTN_NO_SEG_SEQUENCE
        && !matches!(compo_adjust_mode, CompoAdjustMode::NoCompositionBasedStats)
        && *subject_maybe_biased
        && (!should_test_identical
            || !test_near_identical(
                &seq_data,
                range.begin,
                query_data,
                query_range.begin,
                query_words,
                align,
            ))
    {
        let (status, was_biased) = do_seg_sequence_data(&mut seq_data);
        if status != 0 {
            return Err("SEG filtering failed");
        }
        *subject_maybe_biased = was_biased;
    }
    Ok(seq_data)
}

/// blast-rs: In-memory port boundary for `s_SequenceGetRange` (`blast_kappa.c:1670`);
/// not a complete direct NCBI C port.
///
/// The C helper prepares the query range, then dispatches subject extraction to
/// either `s_SequenceGetProteinRange` or `s_SequenceGetTranslatedRange`
/// depending on the program. This Rust helper covers the same dispatch when the
/// caller has already materialized the subject sequence bytes. Protein subjects
/// are NCBIstdaa bytes; translated subjects are NCBI4na bytes.
pub fn sequence_get_range_in_memory(
    program: ProgramType,
    query: &BlastCompoSequenceData,
    query_range: &BlastCompoSequenceRange,
    subject_source: &[u8],
    subject_range: &BlastCompoSequenceRange,
) -> Result<(BlastCompoSequenceData, BlastCompoSequenceData), &'static str> {
    sequence_get_range_in_memory_with_code(
        program,
        query,
        query_range,
        subject_source,
        subject_range,
        &crate::util::STANDARD_GENETIC_CODE,
    )
}

/// blast-rs: Same as [`sequence_get_range_in_memory`], with an explicit genetic code for
/// translated-subject programs; not a direct NCBI C port.
pub fn sequence_get_range_in_memory_with_code(
    program: ProgramType,
    query: &BlastCompoSequenceData,
    query_range: &BlastCompoSequenceRange,
    subject_source: &[u8],
    subject_range: &BlastCompoSequenceRange,
    genetic_code: &[u8; 64],
) -> Result<(BlastCompoSequenceData, BlastCompoSequenceData), &'static str> {
    let query_seq = sequence_prep_query_range(query, query_range);
    let subject_seq = if crate::program::blast_subject_is_translated(program) {
        sequence_get_translated_range(subject_source, subject_range, genetic_code)?
    } else {
        sequence_get_protein_range(subject_source, subject_range)
    };
    Ok((query_seq, subject_seq))
}

/// In-memory counterpart of NCBI `s_SequenceGetRange` for composition redo.
/// The plain helper above intentionally omits SEG; this helper mirrors the
/// redo callback path where the subject fetcher owns near-identical testing,
/// SEG masking, and `subject_maybe_biased` updates.
#[allow(clippy::too_many_arguments)]
pub fn sequence_get_range_in_memory_with_code_for_redo(
    program: ProgramType,
    query: &BlastCompoSequenceData,
    query_range: &BlastCompoSequenceRange,
    subject_source: &[u8],
    subject_range: &BlastCompoSequenceRange,
    genetic_code: &[u8; 64],
    query_words: &[u64],
    align: &BlastCompoAlignment,
    should_test_identical: bool,
    compo_adjust_mode: CompoAdjustMode,
    subject_maybe_biased: &mut bool,
) -> Result<(BlastCompoSequenceData, BlastCompoSequenceData), &'static str> {
    let query_seq = sequence_prep_query_range(query, query_range);
    let subject_seq = if crate::program::blast_subject_is_translated(program) {
        sequence_get_translated_range_for_redo(
            subject_source,
            subject_range,
            genetic_code,
            query_range,
            &query_seq,
            query_words,
            align,
            should_test_identical,
            compo_adjust_mode,
            subject_maybe_biased,
        )?
    } else {
        sequence_get_protein_range_for_redo(
            subject_source,
            subject_range,
            query_range,
            &query_seq,
            query_words,
            align,
            should_test_identical,
            compo_adjust_mode,
            subject_maybe_biased,
        )?
    };
    Ok((query_seq, subject_seq))
}

/// In-memory counterpart of the Smith-Waterman redo `get_range` call.
///
/// NCBI `Blast_RedoOneMatchSmithWaterman` passes `subject_maybe_biased = NULL`
/// while still setting `always_get_subject = TRUE`. That means the subject SEG
/// gate runs for composition-adjusted modes, but no per-match biased-state is
/// cached or updated across windows. The ordinary redo helper above deliberately
/// keeps the mutable state; SW redo should use this helper instead.
#[allow(clippy::too_many_arguments)]
pub fn sequence_get_range_in_memory_with_code_for_sw_redo(
    program: ProgramType,
    query: &BlastCompoSequenceData,
    query_range: &BlastCompoSequenceRange,
    subject_source: &[u8],
    subject_range: &BlastCompoSequenceRange,
    genetic_code: &[u8; 64],
    query_words: &[u64],
    align: &BlastCompoAlignment,
    near_identical_status: bool,
    compo_adjust_mode: CompoAdjustMode,
) -> Result<(BlastCompoSequenceData, BlastCompoSequenceData), &'static str> {
    let mut subject_maybe_biased = true;
    sequence_get_range_in_memory_with_code_for_redo(
        program,
        query,
        query_range,
        subject_source,
        subject_range,
        genetic_code,
        query_words,
        align,
        near_identical_status,
        compo_adjust_mode,
        &mut subject_maybe_biased,
    )
}

/// NCBI: s_GappingParamsNew (blast_kappa.c:2354).
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
pub fn s_gapping_params_new(
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
        (gap_x_dropoff_final_bits * crate::math::NCBIMATH_LN2 / min_lambda)
            .max(raw_gap_x_dropoff_final as f64) as i32
    } else {
        raw_gap_x_dropoff_final
    };

    BlastCompoGappingParams {
        gap_open: scoring.gap_open,
        gap_extend: scoring.gap_extend,
        decline_align: i32::MIN, // unused in the kappa pipeline; NCBI default
        x_dropoff,
        context: None,
    }
}

/// NCBI: s_MatrixInfoInit (blast_kappa.c:2199).
///
/// Initializes a `BlastMatrixInfo` from the matrix name and the
/// "ideal" gapped Karlin block, populating `ungappedLambda`,
/// `startFreqRatios`, and `startMatrix`. The position-based (PSI) branch has
/// a dedicated initializer plus `s_ScalePosMatrix` and Lambda-ratio statistic
/// update helpers.
///
/// `kbp_ideal_lambda` is `sbp->kbp_ideal->Lambda` from the C source —
/// the lambda for the ideal (non-context-specific) ungapped matrix.
pub fn s_matrix_info_init(
    self_: &mut BlastMatrixInfo,
    matrix_name: &str,
    kbp_ideal_lambda: f64,
    scale_factor: f64,
) -> i32 {
    self_.matrix_name = matrix_name.to_string();
    self_.positional = false;
    self_.ungapped_lambda = kbp_ideal_lambda / scale_factor;

    let freq_ratio_info = match crate::matrix::get_matrix_freq_ratios_with_scale(matrix_name) {
        Some(ratios) => ratios,
        None => return -1,
    };
    self_.bit_scale_factor = freq_ratio_info.bit_scale_factor;
    let freq_ratios = freq_ratio_info.data;
    self_.start_freq_ratios = freq_ratios.iter().map(|row| row.to_vec()).collect();

    // C: `Blast_Int4MatrixFromFreq(self->startMatrix, self->cols,
    //                              self->startFreqRatios, self->ungappedLambda);`
    let int_mat =
        crate::composition::blast_int4_matrix_from_freq(self_.ungapped_lambda, &freq_ratios);
    self_.matrix = int_mat.iter().map(|row| row.to_vec()).collect();
    self_.rows = crate::matrix::AA_SIZE as i32;
    self_.cols = crate::matrix::AA_SIZE as i32;
    0
}

/// blast-rs: Position-based counterpart to the frequency-ratio setup in
/// `s_MatrixInfoInit`; not a direct NCBI C port.
///
/// NCBI's full PSI branch calls `s_GetPosBasedStartFreqRatios`, then
/// `s_ScalePosMatrix`, which internally runs `_PSIConvertFreqRatiosToPSSM`.
/// Call [`psi_private_update_lambda_statistics`] after position-based score
/// adjustment when the private Lambda-ratio statistic side effect is needed.
pub fn matrix_info_init_psiblast_from_start_numerator(
    self_: &mut BlastMatrixInfo,
    query: &[u8],
    matrix_name: &str,
    start_numerator: &[Vec<f64>],
    kbp_ideal_lambda: f64,
    scale_factor: f64,
) -> i32 {
    self_.matrix_name = matrix_name.to_string();
    self_.positional = true;
    self_.ungapped_lambda = kbp_ideal_lambda / scale_factor;

    let freq_ratio_info = match crate::matrix::get_matrix_freq_ratios_with_scale(matrix_name) {
        Some(ratios) => ratios,
        None => return -1,
    };
    self_.bit_scale_factor = freq_ratio_info.bit_scale_factor;

    let freq_ratios = match get_pos_based_start_freq_ratios(query, matrix_name, start_numerator) {
        Ok(ratios) => ratios,
        Err(()) => return -1,
    };

    self_.start_freq_ratios = freq_ratios;
    self_.rows = query.len() as i32;
    self_.cols = crate::matrix::AA_SIZE as i32;
    let freq_ratios = self_.start_freq_ratios.clone();
    scale_pos_matrix(self_, &freq_ratios)
}

/// NCBI: s_RecordInitialSearch (blast_kappa.c:2059).
/// naming: Public Rust helper omits the private `s_` prefix.
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
        // Ensure orig_matrix has the right shape — `BlastKappaSavedParameters::s_saved_parameters_new`
        // already allocated it, but defend against caller mismatch.
        if saved.orig_matrix.len() < rows {
            saved
                .orig_matrix
                .resize(rows, vec![0i32; crate::matrix::AA_SIZE]);
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

/// NCBI: s_RestoreSearch (blast_kappa.c:2147).
/// naming: Public Rust helper omits the private `s_` prefix.
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
                let cols = crate::matrix::AA_SIZE
                    .min(row.len())
                    .min(saved.orig_matrix[i].len());
                row[..cols].copy_from_slice(&saved.orig_matrix[i][..cols]);
            }
        }
    }
}

/// NCBI: s_SavedParametersFree (blast_kappa.c:1977).
/// naming: Public Rust helper omits the private `s_` prefix.
///
/// In Rust the struct's `Drop` impl handles deallocation, so this
/// function is a no-op marker — provided so that callers porting the
/// NCBI flow can keep matching call sites. Sets the option slot to
/// `None` to mirror C's `*searchParams = NULL`.
pub fn saved_parameters_free(saved: &mut Option<BlastKappaSavedParameters>) {
    let Some(mut params) = saved.take() else {
        return;
    };
    params.orig_matrix.clear();
    params.kbp_gap_orig.clear();
    params.num_queries = 0;
    params.gap_open = 0;
    params.gap_extend = 0;
    params.scale_factor = 0.0;
    params.original_expect_value = 0.0;
}

/// Internal `BlastCompo_HeapRecord` (`compo_heap.c:75`).
#[derive(Debug, Clone)]
pub struct BlastCompoHeapRecord {
    pub best_evalue: f64,
    pub best_score: i32,
    pub subject_index: i32,
    pub these_alignments: HspList,
}

impl BlastCompoHeapRecord {
    fn from_hsp_list(record: HspList) -> Self {
        let best_score = record.hsps.iter().map(|hsp| hsp.score).max().unwrap_or(0);
        Self {
            best_evalue: record.best_evalue,
            best_score,
            subject_index: record.oid,
            these_alignments: record,
        }
    }

    fn into_hsp_list(self) -> HspList {
        self.these_alignments
    }
}

/// NCBI: BlastCompo_Heap (compo_heap.h:82).
///
/// The member set mirrors C: `n`, `capacity`, `heapThreshold`, `ecutoff`,
/// `worstEvalue`, plus the split `array` / `heapArray` storage. The remaining
/// method bodies are transitional and should be ported next to the C
/// heapify-up/down routines now that the fields exist.
#[derive(Debug, Clone, Default)]
pub struct BlastCompoHeap {
    pub n: i32,
    pub capacity: i32,
    pub heap_threshold: i32,
    pub ecutoff: f64,
    pub worst_evalue: f64,
    pub array: Vec<BlastCompoHeapRecord>,
    pub heap_array: Vec<BlastCompoHeapRecord>,
}

/// `EVALUE_STRETCH` (`redo_alignment.c:1588`).
pub const EVALUE_STRETCH: f64 = 5.0;

impl BlastCompoHeap {
    /// NCBI: BlastCompo_HeapNew (compo_heap.c:50).
    /// naming: Associated constructor on `BlastCompoHeap`.
    pub fn new(heap_threshold: i32, ecutoff: f64) -> Self {
        let capacity = heap_threshold.max(0);
        Self {
            n: 0,
            capacity,
            heap_threshold,
            ecutoff,
            worst_evalue: f64::INFINITY,
            array: Vec::with_capacity(capacity as usize),
            heap_array: Vec::new(),
        }
    }

    /// Transitional accessor while call sites move to C `array`/`heapArray`.
    fn records(&self) -> &[BlastCompoHeapRecord] {
        if self.heap_array.is_empty() {
            &self.array
        } else {
            &self.heap_array
        }
    }

    /// Transitional accessor while call sites move to C `array`/`heapArray`.
    fn records_mut(&mut self) -> &mut Vec<BlastCompoHeapRecord> {
        if self.heap_array.is_empty() {
            &mut self.array
        } else {
            &mut self.heap_array
        }
    }

    /// blast-rs: Worst-evalue accessor for Rust heap state; not a direct NCBI C port.
    ///
    /// Worst (largest) e-value currently in the heap, or `+∞` if empty.
    pub fn worst_evalue(&self) -> f64 {
        let records = self.records();
        records
            .iter()
            .map(|r| r.best_evalue)
            .fold(f64::NEG_INFINITY, f64::max)
            .max(f64::NEG_INFINITY)
            .max(if records.is_empty() {
                f64::INFINITY
            } else {
                f64::NEG_INFINITY
            })
    }

    /// blast-rs: Local heap tie-break helper; not a direct NCBI C port.
    fn best_score(record: &BlastCompoHeapRecord) -> i32 {
        record.best_score
    }

    /// blast-rs: Local heap tie-break helper; not a direct NCBI C port.
    fn record_worse_than(
        record: &BlastCompoHeapRecord,
        evalue: f64,
        score: i32,
        subject_index: i32,
    ) -> bool {
        if record.best_evalue > evalue {
            return true;
        }
        if record.best_evalue < evalue {
            return false;
        }
        let record_score = Self::best_score(record);
        if record_score < score {
            return true;
        }
        if record_score > score {
            return false;
        }
        record.subject_index < subject_index
    }

    /// blast-rs: Local heap scan helper replacing C heap internals; not a direct NCBI C port.
    fn worst_record_index(&self) -> Option<usize> {
        let records = self.records();
        if records.is_empty() {
            return None;
        }
        let mut worst_idx = 0usize;
        for i in 1..records.len() {
            let current = &records[i];
            let worst = &records[worst_idx];
            if Self::record_worse_than(
                current,
                worst.best_evalue,
                Self::best_score(worst),
                worst.subject_index,
            ) {
                worst_idx = i;
            }
        }
        Some(worst_idx)
    }

    /// NCBI: BlastCompo_HeapWouldInsert (compo_heap.c:190).
    /// naming: Rust exposes this as an associated method on the composition
    /// heap type.
    ///
    /// Returns true if a candidate with `(evalue, score, subject_index)` would
    /// be retained: the heap is not full, the candidate passes the inclusion
    /// cutoff, the candidate e-value is better than the current worst, or it
    /// wins the NCBI tie-breaker against the current worst record.
    pub fn would_insert(&self, evalue: f64, score: i32, subject_index: i32) -> bool {
        if self.records().len() < self.heap_threshold.max(0) as usize
            || evalue <= self.ecutoff
            || evalue < self.worst_evalue()
        {
            return true;
        }
        if let Some(worst_idx) = self.worst_record_index() {
            Self::record_worse_than(&self.records()[worst_idx], evalue, score, subject_index)
        } else {
            true
        }
    }

    /// NCBI: BlastCompo_HeapInsert (compo_heap.c:210).
    /// naming: Associated method on `BlastCompoHeap`; type supplies `blast_compo_heap`.
    ///
    /// Inserts `record` if it would be retained. When the heap is full and the
    /// candidate does not pass `ecutoff`, the worse of the candidate and the
    /// current worst record is returned as discarded.
    pub fn insert(&mut self, record: HspList) -> Option<HspList> {
        let new_record = BlastCompoHeapRecord::from_hsp_list(record);
        let evalue = new_record.best_evalue;
        let score = new_record.best_score;
        let subject_index = new_record.subject_index;
        if !self.would_insert(evalue, score, subject_index) {
            return Some(new_record.into_hsp_list());
        }

        let threshold = self.heap_threshold.max(0) as usize;
        if self.records().is_empty()
            || self.records().len() < threshold
            || (evalue <= self.ecutoff && self.worst_evalue() <= self.ecutoff)
        {
            self.records_mut().push(new_record);
            self.n = self.records().len() as i32;
            self.worst_evalue = self.worst_evalue();
            return None;
        }

        if let Some(worst_idx) = self.worst_record_index() {
            if Self::record_worse_than(&self.records()[worst_idx], evalue, score, subject_index) {
                let discarded = std::mem::replace(&mut self.records_mut()[worst_idx], new_record);
                self.n = self.records().len() as i32;
                self.worst_evalue = self.worst_evalue();
                return Some(discarded.into_hsp_list());
            }
        }
        Some(new_record.into_hsp_list())
    }

    /// NCBI: BlastCompo_HeapFilledToCutoff (compo_heap.c:270).
    /// naming: Rust exposes this as an associated method on the composition
    /// heap type.
    pub fn filled_to_cutoff(&self) -> bool {
        self.records().len() >= self.heap_threshold.max(0) as usize
            && self.worst_evalue() <= self.ecutoff
    }

    /// NCBI: BlastCompo_HeapPop (compo_heap.c:251).
    ///
    /// Removes and returns the
    /// hit-list with the worst (largest) e-value.
    /// naming: Rust exposes this as an associated method on the composition
    /// heap type.
    pub fn pop_worst(&mut self) -> Option<HspList> {
        let worst_idx = self.worst_record_index()?;
        let record = self.records_mut().swap_remove(worst_idx);
        self.n = self.records().len() as i32;
        self.worst_evalue = self.worst_evalue();
        Some(record.into_hsp_list())
    }
}

/// NCBI: BlastCompo_EarlyTermination (redo_alignment.c:1592).
pub fn blast_compo_early_termination(
    evalue: f64,
    significant_matches: &[BlastCompoHeap],
    num_queries: usize,
) -> bool {
    if significant_matches.len() < num_queries {
        return false;
    }

    for heap in &significant_matches[..num_queries] {
        if heap.filled_to_cutoff() {
            if evalue <= EVALUE_STRETCH * heap.ecutoff {
                return false;
            }
        } else {
            return false;
        }
    }
    true
}

/// NCBI: s_FreeBlastCompo_QueryInfoArray (blast_kappa.c:2279).
/// naming: Public Rust helper omits the private `s_` prefix.
///
/// In the C source this `free`s every `query_info[i].words` array and
/// then `free`s the outer container. In Rust the `BlastCompoQueryInfo`
/// owns its `words: Vec<u64>` so the `Vec::clear()` + drop suffices.
/// Provided for parity / completeness so callers porting NCBI flow can
/// keep the same call sites.
pub fn free_blast_compo_query_info_array(query_info: &mut Vec<BlastCompoQueryInfo>) {
    query_info.clear();
}

/// NCBI: s_ClearHeap (blast_kappa.c:2518).
/// naming: Public Rust helper omits the private `s_` prefix.
pub fn clear_heap(heap: &mut BlastCompoHeap) {
    while heap.pop_worst().is_some() {}
}

/// blast-rs: Merge worker-local composition heaps into the global per-query heaps;
/// not a direct NCBI C port.
///
/// NCBI's threaded redo core lets worker threads retain candidate HSP lists in
/// local `BlastCompo_Heap` arrays, then merges them by reinserting each worker
/// record into the final query heap. This function mirrors that semantic merge:
/// every worker heap is drained through [`BlastCompoHeap::pop_worst`] and each
/// record is retained or discarded only by the destination heap's normal
/// [`BlastCompoHeap::insert`] rules.
pub fn merge_compo_thread_heaps(
    global_heaps: &mut [BlastCompoHeap],
    worker_heaps: &mut [Vec<BlastCompoHeap>],
) {
    for worker in worker_heaps {
        for (query_index, worker_heap) in worker.iter_mut().enumerate() {
            if let Some(global_heap) = global_heaps.get_mut(query_index) {
                while let Some(hsp_list) = worker_heap.pop_worst() {
                    let _discarded = global_heap.insert(hsp_list);
                }
            } else {
                clear_heap(worker_heap);
            }
        }
    }
}

/// NCBI: s_FillResultsFromCompoHeaps (blast_kappa.c:2493).
/// naming: Public Rust helper omits the private `s_` prefix.
///
/// Drains each heap into a fresh `HitList`, then reverses the
/// query order at the end (mirroring NCBI's
/// `Blast_HSPResultsReverseOrder` finishing call). The C version
/// allocates `Blast_HitListNew(hitlist_size)` per query; we let
/// `HitList::new()` handle that and ignore the size cap (Rust's
/// `Vec` is unbounded; truncate at the call site if needed).
pub fn fill_results_from_compo_heaps(heaps: &mut [BlastCompoHeap]) -> crate::hspstream::HspResults {
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

/// NCBI: s_GetQueryInfo (blast_kappa.c:2308).
/// naming: Public Rust helper omits the private `s_` prefix.
///
/// Builds a `Vec<BlastCompoQueryInfo>` from the per-context fields of a
/// Rust `QueryInfo` plus the concatenated query buffer. For each context
/// (one entry per strand/frame), populates:
///   - `origin` from the context's `query_offset`
///   - `seq` from the slice `query_data[origin .. origin + length]`
///     (copied into an owned `BlastCompoSequenceData`)
///   - `eff_search_space` from the context's `eff_searchsp`
///   - `words` via [`create_word_array`] (1-1 with C's `s_CreateWordArray`)
///   - `composition` via `crate::composition::blast_read_aa_composition` (1-1 with
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
                BlastAminoAcidComposition::default()
            } else {
                let (probs, num_true) = crate::composition::blast_read_aa_composition(
                    &seq_slice,
                    crate::matrix::AA_SIZE,
                );
                BlastAminoAcidComposition::from_prob(probs, num_true as i32)
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

/// NCBI: s_CreateWordArray (blast_kappa.c:2244).
/// naming: Public Rust helper omits the private `s_` prefix.
///
/// Returns an array of 8-mer hashes such that `result[i]` is the hash
/// of `seq[i .. i + 8]`. Used by `s_TestNearIdentical` /
/// `s_find_num_identical` for fast near-identity probes against a query.
///
/// Returns `None` if `seq.len() < 8` (NCBI returns -1).
///
/// **Bit-for-bit fidelity to NCBI:** the C loop is `for (i = 1; i < seq_len - word_size; i++)`,
/// which leaves the last legitimate 8-mer position (`seq_len - word_size`)
/// at calloc'd-zero. That off-by-one is preserved here — the slot is
/// allocated but never written. `s_find_num_identical`'s outer loop
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
    hashes[0] = s_get_hash(&seq[0..WORD_SIZE], WORD_SIZE);
    // C: `for (i = 1; i < seq_len - word_size; i++)` — note `<` not `<=`,
    // intentionally leaving slot [seq_len - word_size] zero.
    let upper = seq.len() - WORD_SIZE;
    for i in 1..upper {
        hashes[i] = ((hashes[i - 1] << 5) & MASK) + seq[i + WORD_SIZE - 1] as u64;
    }
    Some(hashes)
}

/// NCBI: s_Blast_HSPGetNumIdentitiesAndPositives (blast_hits.c:746).
/// naming: Public Rust helper shortens the C suffix to the local identity-only API.
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
                    crate::gapinfo::GapAlignOpType::Del => {
                        s += count as usize;
                    }
                    crate::gapinfo::GapAlignOpType::Ins => {
                        q += count as usize;
                    }
                    _ => {
                        // NCBI `s_Blast_HSPGetNumIdentitiesAndPositives`
                        // (`blast_hits.c:818`) default branch advances both
                        // pointers — frame-shift ops (`Del1/Del2/Ins1/Ins2`)
                        // and `Decline` fall here. OOF alignments are
                        // counted via `s_Blast_HSPGetOOFNumIdentitiesAndPositives`
                        // instead.
                        q += count as usize;
                        s += count as usize;
                    }
                }
            }
        }
    }

    // NCBI `blast_hits.c:832-833`:
    //   if (NULL != matrix) *num_pos_ptr = num_pos + num_ident;
    // The "positives" count returned to callers is the SUM of
    // identical-residue matches AND non-identical-but-positive-matrix
    // matches. Our inner loop only incremented `num_pos` on the
    // non-identical-but-positive branch (matching NCBI's loop body), so
    // we need the final fold-in here when a matrix was supplied.
    if matrix.is_some() {
        num_pos += num_ident;
    }

    (num_ident, align_length, num_pos)
}

/// blast-rs: Protein/materialized-subject identity stamping helper; not a
/// direct NCBI C port.
///
/// Walks every HSP in the list and stamps `hsp.num_ident` from a fresh
/// re-walk over the (unmasked) query and subject bytes via
/// `blast_hsp_get_num_identities`. Pass NCBIstdaa subject bytes for
/// blastp/blastx. For tblastn-style NCBI4na subjects, use
/// [`compute_num_identities_translated_subject`] so the HSP's
/// `subject_frame` selects the translated target frame.
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
    let hsp_count = hsp_list.hsps.len();
    for hsp_index in 0..hsp_count {
        let hsp = &mut hsp_list.hsps[hsp_index];
        let ops = edit_ops_per_hsp
            .get(hsp_index)
            .and_then(|ops_for_hsp| ops_for_hsp.as_deref());
        let (num_ident, _align_length, _num_pos) =
            blast_hsp_get_num_identities(query_nomask, subject, hsp, ops, matrix);
        hsp.num_ident = num_ident;
    }
}

/// blast-rs: Translated-subject branch of `s_ComputeNumIdentities`; not a direct NCBI C port.
///
/// NCBI builds a `BlastTargetTranslation` for each tblastn/tblastx subject
/// frame before calling `Blast_HSPGetNumIdentitiesAndPositives`. This helper
/// mirrors that flow for an in-memory NCBI4na subject. `hsp.subject_offset` and
/// `hsp.subject_end` are expected to be protein-space offsets in the translated
/// frame selected by `hsp.subject_frame`.
pub fn compute_num_identities_translated_subject(
    query_nomask: &[u8],
    subject_ncbi4na: &[u8],
    hsp_list: &mut HspList,
    edit_ops_per_hsp: &[Option<Vec<(crate::gapinfo::GapAlignOpType, i32)>>],
    matrix: Option<&[[i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE]>,
    genetic_code: &[u8; 64],
) -> Result<(), &'static str> {
    let (translation, offsets) = crate::util::blast_get_all_translations(
        subject_ncbi4na,
        subject_ncbi4na.len(),
        genetic_code,
    );
    for (i, hsp) in hsp_list.hsps.iter_mut().enumerate() {
        let context = match hsp.subject_frame {
            1 => 0,
            2 => 1,
            3 => 2,
            -1 => 3,
            -2 => 4,
            -3 => 5,
            _ => return Err("invalid translated subject frame"),
        };
        let frame_begin = (offsets[context] + 1) as usize;
        let frame_end = offsets[context + 1] as usize;
        let subject = &translation[frame_begin..frame_end];
        let ops = edit_ops_per_hsp.get(i).and_then(|o| o.as_deref());
        let (num_ident, _align_length, _num_pos) =
            blast_hsp_get_num_identities(query_nomask, subject, hsp, ops, matrix);
        hsp.num_ident = num_ident;
    }
    Ok(())
}

fn compute_num_identities_translated_subject_by_context(
    query_info: &[BlastCompoQueryInfo],
    fallback_query: &[u8],
    subject_ncbi4na: &[u8],
    hsp_list: &mut HspList,
    edit_ops_per_hsp: &[Option<Vec<(crate::gapinfo::GapAlignOpType, i32)>>],
    matrix: Option<&[[i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE]>,
    genetic_code: &[u8; 64],
) -> Result<(), &'static str> {
    let (translation, offsets) = crate::util::blast_get_all_translations(
        subject_ncbi4na,
        subject_ncbi4na.len(),
        genetic_code,
    );
    for (i, hsp) in hsp_list.hsps.iter_mut().enumerate() {
        let subject_context = match hsp.subject_frame {
            1 => 0,
            2 => 1,
            3 => 2,
            -1 => 3,
            -2 => 4,
            -3 => 5,
            _ => return Err("invalid translated subject frame"),
        };
        let frame_begin = (offsets[subject_context] + 1) as usize;
        let frame_end = offsets[subject_context + 1] as usize;
        let subject = &translation[frame_begin..frame_end];
        let query_context = hsp.context.max(0) as usize;
        let query_nomask = query_info
            .get(query_context)
            .map(|info| info.seq.data())
            .unwrap_or(fallback_query);
        let ops = edit_ops_per_hsp.get(i).and_then(|o| o.as_deref());
        let (num_ident, _align_length, _num_pos) =
            blast_hsp_get_num_identities(query_nomask, subject, hsp, ops, matrix);
        hsp.num_ident = num_ident;
    }
    Ok(())
}

fn compute_num_identities_protein_by_context(
    query_info: &[BlastCompoQueryInfo],
    fallback_query: &[u8],
    subject: &[u8],
    hsp_list: &mut HspList,
    edit_ops_per_hsp: &[Option<Vec<(crate::gapinfo::GapAlignOpType, i32)>>],
    matrix: Option<&[[i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE]>,
) {
    for (i, hsp) in hsp_list.hsps.iter_mut().enumerate() {
        let query_context = hsp.context.max(0) as usize;
        let query_nomask = query_info
            .get(query_context)
            .map(|info| info.seq.data())
            .unwrap_or(fallback_query);
        let ops = edit_ops_per_hsp.get(i).and_then(|o| o.as_deref());
        let (num_ident, _align_length, _num_pos) =
            blast_hsp_get_num_identities(query_nomask, subject, hsp, ops, matrix);
        hsp.num_ident = num_ident;
    }
}

/// Link-HSP context needed by the sum-statistics branch of
/// `s_HitlistEvaluateAndPurge`.
///
/// NCBI has these pieces available through `BlastScoreBlk`,
/// `BlastQueryInfo`, and hit-saving parameters. Keeping them explicit here
/// makes the Rust call site choose whether it is truly ready to take the
/// `blast_link_hsps` path.
pub struct HitlistLinkContext<'a> {
    pub query_info: &'a crate::queryinfo::QueryInfo,
    pub query_context: usize,
    pub score_block: &'a crate::link_hsps::LinkScoreBlock,
    pub link_params: &'a crate::link_hsps::LinkHSPParameters,
    pub gapped_calculation: bool,
}

/// NCBI: s_HitlistEvaluateAndPurge (blast_kappa.c:394).
///
/// Assigns final e-values and prunes the hit list. The full C function
/// dispatches between `blast_link_hsps` (sum-statistics branch) and
/// `Blast_HSPListGetEvalues` (single-HSP branch), then applies
/// `s_AdjustEvaluesForComposition` for blastp/blastx, then drops HSPs
/// above the evalue threshold and reports `(best_score, best_evalue)`.
///
/// This Rust port implements the single-HSP branch and the link-HSPs branch.
/// For the single-HSP branch, pass `Some(kbp)` to recompute per-HSP
/// e-values exactly like `Blast_HSPListGetEvalues`; pass
/// `None` when the caller has already populated those fields with a
/// Spouge or translated-search-specific statistic.
///
/// Returns `(best_score, best_evalue)`. After the call, `hsp_list.hsps`
/// is sorted-stable by NCBI's tie-breaker and contains only HSPs whose
/// e-value passes `max_evalue`.
#[allow(clippy::too_many_arguments)]
pub fn s_hitlist_evaluate_and_purge(
    hsp_list: &mut HspList,
    subject_length: i32,
    program_number: ProgramType,
    query_length: i32,
    length_adjustment: i32,
    eff_searchsp: f64,
    kbp: Option<&crate::stat::KarlinBlk>,
    link_context: Option<&HitlistLinkContext<'_>>,
    pvalue_for_this_pair: f64,
    max_evalue: f64,
    do_sum_stats: bool,
) -> (i32, f64) {
    // C: `*pbestEvalue = DBL_MAX; *pbestScore = 0;`
    let mut best_evalue = f64::MAX;
    let mut best_score = 0i32;

    if do_sum_stats {
        let link_context = link_context
            .expect("s_hitlist_evaluate_and_purge: do_sum_stats requires Blast_LinkHsps context");
        let status = blast_link_hsps_for_kappa(
            hsp_list,
            program_number,
            subject_length,
            link_context.query_info,
            link_context.query_context,
            link_context.score_block,
            link_context.link_params,
            link_context.gapped_calculation,
        );
        if status != 0 {
            hsp_list.hsps.clear();
            hsp_list.best_evalue = i32::MAX as f64;
            return (best_score, best_evalue);
        }
    }
    if !do_sum_stats {
        if let Some(kbp) = kbp {
            blast_hsp_list_get_evalues_simple(hsp_list, kbp, eff_searchsp);
        }
    }

    // Composition adjustment for blastp/blastx if the p-value is in [0, 1].
    if (program_number == crate::program::BLASTP || program_number == crate::program::BLASTX)
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

/// blast-rs: Bridge from the legacy Rust `HspList` representation to the
/// upstream-shaped `BlastHSPList` linker and back.
pub fn blast_link_hsps_for_kappa(
    hsp_list: &mut HspList,
    program_number: ProgramType,
    subject_length: i32,
    query_info: &crate::queryinfo::QueryInfo,
    query_context: usize,
    score_block: &crate::link_hsps::LinkScoreBlock,
    link_params: &crate::link_hsps::LinkHSPParameters,
    gapped_calculation: bool,
) -> i32 {
    let source_hsp_max = hsp_list.hsp_max;
    let source_list = std::mem::replace(hsp_list, HspList::new(0));
    let mut blast_hsp_list = crate::hspstream::BlastHSPList::from_legacy_hsp_list(source_list);
    blast_hsp_list.query_index = query_info
        .contexts
        .get(query_context)
        .map(|ctx| ctx.query_index)
        .unwrap_or_else(|| {
            crate::queryinfo::blast_get_query_index_from_context(
                query_context as i32,
                program_number,
            )
        });
    let status = crate::link_hsps::blast_link_hsp_list(
        program_number,
        &mut blast_hsp_list,
        query_info,
        subject_length,
        score_block,
        link_params,
        gapped_calculation,
    );
    *hsp_list = blast_hsp_list.into_legacy_hsp_list();
    hsp_list.hsp_max = source_hsp_max;
    status
}

/// blast-rs: Count gap-opening operations in a `GapEditScript`; not a direct NCBI C port.
///
/// NCBI's HSP initialization stores the number of gap runs, not the total
/// number of gapped residues. Consecutive same-type gap ops are already merged
/// by `GapEditScript::push`, so each insertion/deletion op contributes one.
pub fn gap_edit_script_num_gap_opens(script: &crate::gapinfo::GapEditScript) -> i32 {
    script
        .iter()
        .filter(|(op, _)| {
            matches!(
                op,
                crate::gapinfo::GapAlignOpType::Del
                    | crate::gapinfo::GapAlignOpType::Del1
                    | crate::gapinfo::GapAlignOpType::Del2
                    | crate::gapinfo::GapAlignOpType::Ins
                    | crate::gapinfo::GapAlignOpType::Ins1
                    | crate::gapinfo::GapAlignOpType::Ins2
            )
        })
        .count() as i32
}

/// blast-rs: Compact compatibility helper for call sites that already collapsed
/// `BlastScoreBlk`/`BlastQueryInfo` to a single Karlin block and search space.
///
pub fn blast_hsp_list_get_evalues_simple(
    hsp_list: &mut HspList,
    kbp: &crate::stat::KarlinBlk,
    search_space: f64,
) {
    // NCBI `Blast_HSPListGetEvalues` (`blast_hits.c:1902`) closes with
    // `s_BlastGetBestEvalue` which seeds with `(double)INT4_MAX`.
    let mut best_evalue = i32::MAX as f64;
    for hsp in &mut hsp_list.hsps {
        hsp.evalue = kbp.raw_to_evalue(hsp.score, search_space);
        if hsp.evalue < best_evalue {
            best_evalue = hsp.evalue;
        }
    }
    hsp_list.best_evalue = best_evalue;
}

/// NCBI: Blast_HSPListGetBitScores (blast_hits.c:1907), compact single-Karlin
/// adapter for call sites that do not yet carry a full `BlastScoreBlk`.
pub fn blast_hsp_list_get_bit_scores_simple(hsp_list: &mut HspList, kbp: &crate::stat::KarlinBlk) {
    for hsp in &mut hsp_list.hsps {
        hsp.bit_score = kbp.raw_to_bit(hsp.score);
    }
}

/// NCBI: Blast_HSPListGetEvalues (blast_hits.c:1811).
///
/// This keeps the upstream-shaped inputs: program kind, query contexts,
/// subject length, gapped-vs-ungapped Karlin arrays, RPS preliminary context
/// selection, gap-decay division, RPS score scaling, Spouge finite-size
/// correction when `sbp->gbp` is present, and score rounding through
/// `sbp->round_down`.
#[allow(clippy::too_many_arguments)]
pub fn blast_hsp_list_get_evalues(
    program_number: ProgramType,
    query_info: &crate::queryinfo::QueryInfo,
    subject_length: i32,
    hsp_list: &mut HspList,
    gapped_calculation: bool,
    rps_prelim: bool,
    sbp: &crate::stat::BlastScoreBlk,
    gap_decay_rate: f64,
    scaling_factor: f64,
) -> i16 {
    if hsp_list.hsps.is_empty() {
        return 0;
    }

    let kbp_array = if gapped_calculation {
        &sbp.kbp_gap
    } else {
        &sbp.kbp
    };
    if kbp_array.is_empty() || scaling_factor == 0.0 {
        return 0;
    }

    let is_rps = crate::program::blast_program_is_rps_blast(program_number);
    let blast_gap_decay_divisor = if gap_decay_rate != 0.0 {
        crate::stat::blast_gap_decay_divisor(gap_decay_rate, 1)
    } else {
        1.0
    };

    // NCBI `Blast_HSPListGetEvalues` (`blast_hits.c:1902`) calls
    // `s_BlastGetBestEvalue` which seeds with `(double)INT4_MAX`.
    let mut best_evalue = i32::MAX as f64;
    for hsp in &mut hsp_list.hsps {
        let hsp_context = hsp.context.max(0) as usize;
        let mut kbp_context = hsp_context;
        if rps_prelim {
            kbp_context = kbp_array
                .iter()
                .position(crate::stat::KarlinBlk::is_valid)
                .unwrap_or(kbp_context);
        }

        let Some(kbp) = kbp_array.get(kbp_context).or_else(|| kbp_array.first()) else {
            continue;
        };
        let Some(context_info) = query_info
            .contexts
            .get(hsp_context)
            .or_else(|| query_info.contexts.first())
        else {
            continue;
        };

        let mut scaled_kbp = kbp.clone();
        scaled_kbp.lambda /= scaling_factor;
        scaled_kbp.round_down = false;

        let score = if gapped_calculation && sbp.round_down {
            hsp.score & !1
        } else {
            hsp.score
        };

        hsp.evalue = if let Some(gbp) = sbp
            .gbp
            .as_ref()
            .filter(|gbp| crate::stat::gumbel_blk_is_filled(gbp))
        {
            let query_length = context_info.query_length;
            if is_rps {
                crate::stat::blast_spouge_sto_e(
                    score,
                    Some(&scaled_kbp),
                    Some(gbp),
                    subject_length,
                    query_length,
                )
            } else {
                crate::stat::blast_spouge_sto_e(
                    score,
                    Some(&scaled_kbp),
                    Some(gbp),
                    query_length,
                    subject_length,
                )
            }
        } else {
            let search_space = context_info.eff_searchsp as f64;
            scaled_kbp.raw_to_evalue(score, search_space)
        };
        hsp.evalue /= blast_gap_decay_divisor;

        if hsp.evalue < best_evalue {
            best_evalue = hsp.evalue;
        }
    }
    hsp_list.best_evalue = best_evalue;
    0
}

/// NCBI: s_AdjustEvaluesForComposition (blast_kappa.c:134).
/// naming: Public Rust helper omits the private `s_` prefix.
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
        let align_p_value = crate::composition::blast_karlin_eto_p(hsp.evalue);
        let combined_p_value =
            crate::composition::blast_overall_p_value(comp_p_value, align_p_value);
        hsp.evalue = crate::composition::blast_karlin_pto_e(combined_p_value);
        hsp.evalue /= db_to_sequence_scale;

        if hsp.evalue < best_evalue {
            best_evalue = hsp.evalue;
        }
    }

    hsp_list.best_evalue = best_evalue;
}

/// NCBI: s_TestNearIdentical (blast_kappa.c:1258).
/// naming: Public Rust helper omits the private `s_` prefix.
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
    let (mut num_identical, query_right_len, subject_right_len, _) = s_extend_right(
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
    let (left_ident, query_left_len, subject_left_len, _) = s_extend_left(
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
    num_identical += s_find_num_identical(
        &query[q_start + qr..q_start + qr + mid_q_len],
        q_words_slice,
        &subject[s_start + sr..s_start + sr + mid_s_len],
        MAX_SHIFT,
    );

    let fraction_identical = num_identical as f64 / align_len;
    fraction_identical > K_MIN_FRACTION_NEAR_IDENTICAL
}

/// NCBI: s_HSPListFromDistinctAlignments (blast_kappa.c:304).
/// naming: Public Rust helper omits the private `s_` prefix.
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
/// The Rust `Hsp` stores `comp_adjustment_method` as the raw
/// `ECompoAdjustModes` integer and owns the transferred edit script.
/// The return value keeps the historical parallel tag vector for existing
/// call sites, but the primary state now lives on each HSP.
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
        let edit_script = node
            .context
            .take()
            .and_then(|mut ctx| ctx.edit_script.take());
        let num_gaps = edit_script
            .as_ref()
            .map(gap_edit_script_num_gap_opens)
            .unwrap_or(0);
        // C: `Blast_HSPInit(queryStart, queryEnd, matchStart, matchEnd,
        //                   unknown_value, unknown_value,
        //                   queryIndex, frame, align->frame, score,
        //                   &editScript, &new_hsp);`
        let tag = match node.matrix_adjust_rule {
            MatrixAdjustRule::DontAdjust => CompoAdjustMode::NoCompositionBasedStats,
            MatrixAdjustRule::ScaleOldMatrix => CompoAdjustMode::CompositionBasedStats,
            _ => CompoAdjustMode::CompositionMatrixAdjust,
        };
        let hsp = Hsp {
            score: node.score,
            num_ident: 0, // C: explicitly leaves num_ident blank
            bit_score: 0.0,
            evalue: 0.0,
            query_offset: node.query_start,
            query_end: node.query_end,
            query_gapped_start: 0,
            subject_offset: node.match_start,
            subject_end: node.match_end,
            subject_gapped_start: 0,
            context: node.query_index,
            query_frame: frame,
            subject_frame: node.frame,
            num_gaps,
            comp_adjustment_method: tag as i32,
            edit_script,
            pat_info: None,
            map_info: None,
        };
        // Append HSP and its compositional tag in parallel.
        hsp_list.hsps.push(hsp);
        comp_tags.push(tag);

        cursor = node.next.take();
    }

    // C: `Blast_HSPListSortByScore(hsp_list)`.
    hsp_list.sort_by_score();
    (0, comp_tags)
}

/// NCBI: s_NewAlignmentFromGapAlign (blast_kappa.c:1747).
/// naming: Public Rust helper omits the private `s_` prefix.
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
    let context = edit_script
        .take()
        .map(|edit_script| BlastCompoAlignmentContext {
            edit_script: Some(edit_script),
            hsp: None,
        });
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

/// NCBI: s_HitlistReapContained (blast_kappa.c:223).
///
/// Removes any HSP that is fully contained within an earlier (higher-
/// scoring) HSP on both query and subject coordinates and shares the
/// same query/subject frame. The hitlist must already be sorted by
/// significance before calling. Operates in place.
pub fn s_hitlist_reap_contained(hsps: &mut Vec<Hsp>) {
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
            // NCBI requires the same query and subject frame before an
            // earlier HSP can suppress a contained later HSP.
            if h1.query_frame != h2.query_frame || h1.subject_frame != h2.subject_frame {
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

/// NCBI: s_CalcLambda (blast_kappa.c:551).
///
/// Newton-Raphson refinement
/// of `lambda` from a score-probability distribution. Takes the
/// probabilities as a slice indexed `[0..score_max - score_min + 1]`.
pub fn s_calc_lambda(probs: &[f64], min_score: i32, max_score: i32, lambda0: f64) -> f64 {
    let score_range = (max_score - min_score + 1) as usize;
    debug_assert_eq!(probs.len(), score_range);
    let mut avg = 0.0f64;
    for i in 0..score_range {
        avg += (min_score + i as i32) as f64 * probs[i];
    }
    let _ = avg; // C populates `freq.score_avg` but `Blast_KarlinLambdaNR`
                 // recomputes it internally; our blast_karlin_lambda_nr does the
                 // same.
    crate::composition::blast_karlin_lambda_nr(probs, min_score, max_score, lambda0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::program::{BLASTP, BLASTX, RPS_BLAST, RPS_TBLASTN};

    #[test]
    fn get_subject_length_passthrough_for_blastp() {
        assert_eq!(s_get_subject_length(1234, BLASTP), 1234);
        assert_eq!(s_get_subject_length(1234, BLASTX), 1234);
        // NCBI: ordinary RPS-blastp does NOT divide by 3; only RPS-tblastn does.
        assert_eq!(s_get_subject_length(1234, RPS_BLAST), 1234);
    }

    #[test]
    fn get_subject_length_rps_tblastn_uses_macro() {
        // GET_NUCL_LENGTH(l) = (l-5)/2 + 2; then -1 then /3.
        // For l=11: (11-5)/2 + 2 = 5; (5-1)/3 = 1.
        assert_eq!(s_get_subject_length(11, RPS_TBLASTN), 1);
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
        assert_eq!(s_get_hash(&[1, 2, 3, 4], 4), 34916);
    }

    #[test]
    fn extend_right_full_match() {
        let q = b"AAAACCCC";
        let s = b"AAAACCCC";
        let (n_ident, q_ext, s_ext, _) = s_extend_right(q, s, 4);
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
        let (n_ident, _q, _s, _) = s_extend_right(q, s, 4);
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
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: 100,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        // Apply lambda=0.3, logK=-2.0, divisor=2.0
        s_hsp_list_normalize_scores(&mut list, 0.3, -2.0, 2.0);
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
                query_gapped_start: 0,
                subject_offset: 0,
                subject_end: 100,
                subject_gapped_start: 0,
                context: 0,
                query_frame: 0,
                subject_frame: 0,
                num_gaps: 0,
                comp_adjustment_method: 0,
                edit_script: None,
                pat_info: None,
                map_info: None,
            },
            Hsp {
                // Contained, lower score → dropped.
                score: 50,
                num_ident: 0,
                bit_score: 0.0,
                evalue: 0.0,
                query_offset: 20,
                query_end: 60,
                query_gapped_start: 20,
                subject_offset: 20,
                subject_end: 60,
                subject_gapped_start: 20,
                context: 0,
                query_frame: 0,
                subject_frame: 0,
                num_gaps: 0,
                comp_adjustment_method: 0,
                edit_script: None,
                pat_info: None,
                map_info: None,
            },
            Hsp {
                // Different frame pair -> not contained even though inside coords.
                score: 50,
                num_ident: 0,
                bit_score: 0.0,
                evalue: 0.0,
                query_offset: 20,
                query_end: 60,
                query_gapped_start: 20,
                subject_offset: 20,
                subject_end: 60,
                subject_gapped_start: 20,
                context: 0,
                query_frame: 0,
                subject_frame: -1,
                num_gaps: 0,
                comp_adjustment_method: 0,
                edit_script: None,
                pat_info: None,
                map_info: None,
            },
            Hsp {
                // Different context alone is not enough to avoid pruning once
                // the explicit frame pair matches.
                score: 40,
                num_ident: 0,
                bit_score: 0.0,
                evalue: 0.0,
                query_offset: 30,
                query_end: 50,
                query_gapped_start: 30,
                subject_offset: 30,
                subject_end: 50,
                subject_gapped_start: 30,
                context: 3,
                query_frame: 0,
                subject_frame: 0,
                num_gaps: 0,
                comp_adjustment_method: 0,
                edit_script: None,
                pat_info: None,
                map_info: None,
            },
        ];
        s_hitlist_reap_contained(&mut hsps);
        assert_eq!(hsps.len(), 2);
        assert_eq!(hsps[0].score, 100);
        assert_eq!(hsps[1].subject_frame, -1);
    }

    #[test]
    fn contained_in_hsp_simple() {
        // a=0 b=10 c=5: yes; d=0 e=10 f=7: yes → true.
        assert!(contained_in_hsp(0, 10, 5, 0, 10, 7));
        // c=11 outside [0,10] → false.
        assert!(!contained_in_hsp(0, 10, 11, 0, 10, 7));
    }
}
