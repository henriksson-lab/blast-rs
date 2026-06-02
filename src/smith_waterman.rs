//! 1-1 ports of NCBI's Smith-Waterman primitives from
//! `composition_adjustment/smith_waterman.c`.
//!
//! These are used by the composition-adjusted alignment redo flow
//! (`Blast_RedoAlignmentCore_MT` in `blast_kappa.c`) to compute the
//! optimal SW alignment under the adjusted scoring matrix, which can
//! differ in extent from the original ungapped/X-drop alignment found
//! in the preliminary search. The two-pass forward-then-reverse design
//! mirrors NCBI's flow exactly:
//!
//!   1. `blast_smith_waterman_score_only` — forward SW DP, finds
//!       optimal end position `(queryEnd, matchSeqEnd)` and score,
//!       no traceback.
//!   2. `blast_smith_waterman_find_start` — reverse SW DP from the
//!       optimal end, finds the matching start `(queryStart,
//!       matchSeqStart)`. Stops as soon as the previous score is
//!       reached (`bestScore >= score_in`).
//!
//! Forbidden ranges (used when multiple HSPs are taken from the same
//! match sequence) are tracked by [`BlastForbiddenRanges`]. The wrappers
//! [`blast_smith_waterman_score_only`] and
//! [`blast_smith_waterman_find_start`] dispatch to the basic variants
//! when the range set is empty and to [`bl_special_smith_waterman_score_only`]
//! / [`bl_special_smith_waterman_find_start`] otherwise — matching
//! NCBI's `Blast_SmithWaterman*` dispatch in `smith_waterman.c`.
//!
//! Both routines operate on NCBIstdaa-encoded protein sequences (or
//! 4-bit BLASTNA for nucleotide), where the encoded byte indexes
//! directly into the matrix row.

use crate::gapinfo::{GapAlignOpType, GapEditScript};
use crate::hspstream::{blast_hsp_list_new, blast_hsp_list_save_hsp, Hsp, HspList};
use crate::program::ProgramType;
use crate::queryinfo::QueryInfo;

/// Per-position gap state for the SW DP — 1-1 port of NCBI's
/// `SwGapInfo` (`smith_waterman.c:38`).
#[derive(Debug, Clone, Copy)]
struct SwGapInfo {
    /// score if no gap is open at this position
    no_gap: i32,
    /// score if a gap is open at this position
    gap_exists: i32,
}

#[derive(Debug, Clone, Copy, Default)]
struct BlastSwGapDp {
    best: i32,
    best_gap: i32,
}

/// Values for the editing script operations in `blast_sw.c`.
const EDIT_SUB: u8 = GapAlignOpType::Sub as u8;
const EDIT_GAP_IN_A: u8 = GapAlignOpType::Del as u8;
const EDIT_GAP_IN_B: u8 = GapAlignOpType::Ins as u8;
const EDIT_OP_MASK: u8 = 0x07;
const EDIT_START_GAP_A: u8 = 0x10;
const EDIT_START_GAP_B: u8 = 0x20;

#[derive(Debug, Clone, Copy, Default)]
struct BlastGapSwTraceDp {
    best: i32,
    best_gap: i32,
    path_score: i32,
    path_stop_i: usize,
    path_stop_j: usize,
}

/// Recovered local alignment from `SmithWatermanScoreWithTraceback`.
///
/// Coordinates are half-open offsets in the original `(A, B)` order passed by
/// the caller. `b_start`/`b_end` include `start_shift`, matching NCBI's final
/// `Blast_HSPAdjustSubjectOffset` step.
#[derive(Debug, Clone)]
pub struct SmithWatermanTracebackHit {
    pub score: i32,
    pub a_start: usize,
    pub a_end: usize,
    pub b_start: usize,
    pub b_end: usize,
    pub edit_script: GapEditScript,
}

pub struct BlastSmithWatermanGappedScoreParams<'a> {
    pub program_number: ProgramType,
    pub query: &'a [u8],
    pub query_info: &'a QueryInfo,
    pub subject: &'a [u8],
    pub subject_length: usize,
    pub subject_oid: i32,
    pub subject_frame: i32,
    pub protein_matrix: &'a [[i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE],
    pub pssm: Option<&'a [Vec<i32>]>,
    pub nucleotide_matrix:
        &'a [[i32; crate::encoding::BLASTNA_SIZE]; crate::encoding::BLASTNA_SIZE],
    pub gap_open: i32,
    pub gap_extend: i32,
    pub cutoffs: &'a [i32],
    pub hsp_num_max: i32,
}

/// 1-1 port of `s_SmithWatermanScoreOnly` from `blast_sw.c`.
///
/// This is the BLAST traceback-stage SW score-only primitive, distinct from
/// the composition-adjustment `Blast_SmithWatermanScoreOnly` helpers below.
/// For ordinary square matrices it mirrors C's sequence swap when `A` is
/// shorter than `B`; for PSSM input it keeps the query-position rows fixed.
pub fn s_smith_waterman_score_only(
    a: &[u8],
    b: &[u8],
    matrix: &[[i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE],
    pssm: Option<&[Vec<i32>]>,
    gap_open: i32,
    gap_extend: i32,
) -> i32 {
    if a.is_empty() || b.is_empty() {
        return 0;
    }

    let (a, b) = if pssm.is_none() && a.len() < b.len() {
        (b, a)
    } else {
        (a, b)
    };

    let gap_open_extend = gap_open + gap_extend;
    let mut scores = vec![BlastSwGapDp::default(); b.len() + 1];
    let mut final_best_score = 0;

    for i in 1..=a.len() {
        let matrix_row: &[i32] = match pssm {
            Some(pssm) => &pssm[i - 1],
            None => &matrix[a[i - 1] as usize],
        };
        let mut insert_score = 0;
        let mut row_score = 0;

        for j in 1..=b.len() {
            let mut best_score = scores[j].best_gap - gap_extend;
            if scores[j].best - gap_open_extend > best_score {
                best_score = scores[j].best - gap_open_extend;
            }
            scores[j].best_gap = best_score;

            best_score = insert_score - gap_extend;
            if row_score - gap_open_extend > best_score {
                best_score = row_score - gap_open_extend;
            }
            insert_score = best_score;

            best_score = (scores[j - 1].best + matrix_row[b[j - 1] as usize]).max(0);
            if insert_score > best_score {
                best_score = insert_score;
            }
            if scores[j].best_gap > best_score {
                best_score = scores[j].best_gap;
            }

            if best_score > final_best_score {
                final_best_score = best_score;
            }

            scores[j - 1].best = row_score;
            row_score = best_score;
        }

        scores[b.len()].best = row_score;
    }

    final_best_score
}

fn blast_sw_frame_to_context_blastx(frame: i32) -> i32 {
    match frame {
        1 => 0,
        2 => 1,
        3 => 2,
        -1 => 3,
        -2 => 4,
        -3 => 5,
        _ => 0,
    }
}

fn blast_sw_cutoff_for_context(
    params: &BlastSmithWatermanGappedScoreParams<'_>,
    context: usize,
) -> i32 {
    if crate::program::blast_program_is_rps_blast(params.program_number) {
        let mut rps_context = params.subject_oid.max(0) as usize;
        if params.program_number == crate::program::RPS_TBLASTN {
            rps_context = rps_context * crate::util::NUM_FRAMES
                + blast_sw_frame_to_context_blastx(params.subject_frame).max(0) as usize;
        }
        params.cutoffs.get(rps_context).copied().unwrap_or(0)
    } else {
        params.cutoffs.get(context).copied().unwrap_or(0)
    }
}

fn blast_sw_pssm_rows_for_context<'a>(
    pssm: Option<&'a [Vec<i32>]>,
    query_offset: i32,
    query_length: i32,
) -> Option<&'a [Vec<i32>]> {
    let pssm = pssm?;
    let start = query_offset.max(0) as usize;
    let len = query_length.max(0) as usize;
    if start + len <= pssm.len() {
        Some(&pssm[start..start + len])
    } else if len <= pssm.len() {
        Some(&pssm[..len])
    } else {
        None
    }
}

/// Rust ownership equivalent of NCBI `BLAST_SmithWatermanGetGappedScore`
/// (`blast_sw.c:629`).
///
/// The C function ignores initial ungapped hits, scans each valid query
/// context against the subject with score-only Smith-Waterman, and saves one
/// placeholder HSP per context whose score reaches the context cutoff. The
/// placeholder carries the score, context, query frame, and subject frame into
/// the later traceback phase; this port preserves that boundary and returns
/// the `HspList` directly.
pub fn blast_smith_waterman_get_gapped_score(
    params: &BlastSmithWatermanGappedScoreParams<'_>,
    hsp_list: Option<HspList>,
) -> (i16, HspList) {
    let mut hsp_list = hsp_list.unwrap_or_else(|| {
        let mut list = blast_hsp_list_new(crate::hspstream::blast_hsp_num_max(
            true,
            params.hsp_num_max,
        ));
        list.oid = params.subject_oid;
        list
    });

    let is_prot = params.program_number != crate::program::BLASTN
        && params.program_number != crate::program::PHI_BLASTN
        && params.program_number != crate::program::MAPPING;

    for (context, curr_ctx) in params.query_info.contexts.iter().enumerate() {
        if !curr_ctx.is_valid {
            continue;
        }
        let query_start = curr_ctx.query_offset.max(0) as usize;
        let query_len = curr_ctx.query_length.max(0) as usize;
        let Some(query_end) = query_start.checked_add(query_len) else {
            continue;
        };
        let Some(query_segment) = params.query.get(query_start..query_end) else {
            continue;
        };

        let cutoff_score = blast_sw_cutoff_for_context(params, context);
        let score = if is_prot {
            let pssm = blast_sw_pssm_rows_for_context(
                params.pssm,
                curr_ctx.query_offset,
                curr_ctx.query_length,
            );
            s_smith_waterman_score_only(
                query_segment,
                params.subject,
                params.protein_matrix,
                pssm,
                params.gap_open,
                params.gap_extend,
            )
        } else {
            s_nucl_smith_waterman(
                params.subject,
                params.subject_length,
                query_segment,
                params.nucleotide_matrix,
                params.gap_open,
                params.gap_extend,
            )
        };

        if score >= cutoff_score {
            let hsp = Hsp {
                score,
                num_ident: 0,
                bit_score: 0.0,
                evalue: f64::MAX,
                query_offset: 0,
                query_end: curr_ctx.query_length,
                query_gapped_start: 0,
                subject_offset: 0,
                subject_end: if is_prot {
                    params.subject.len() as i32
                } else {
                    params.subject_length as i32
                },
                subject_gapped_start: 0,
                context: context as i32,
                query_frame: curr_ctx.frame,
                subject_frame: params.subject_frame,
                num_gaps: 0,
                comp_adjustment_method: 0,
                edit_script: None,
                pat_info: None,
                map_info: None,
            };
            let _ = blast_hsp_list_save_hsp(&mut hsp_list, hsp);
        }
    }

    (0, hsp_list)
}

fn blast_sw_traceback_op(script: u8) -> Option<GapAlignOpType> {
    match script {
        EDIT_SUB => Some(GapAlignOpType::Sub),
        EDIT_GAP_IN_A => Some(GapAlignOpType::Del),
        EDIT_GAP_IN_B => Some(GapAlignOpType::Ins),
        _ => None,
    }
}

fn blast_sw_score(
    a: &[u8],
    b: &[u8],
    matrix: &[[i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE],
    pssm: Option<&[Vec<i32>]>,
    i: usize,
    j: usize,
) -> Option<i32> {
    if i == 0 || j == 0 {
        return None;
    }
    match pssm {
        Some(pssm) => pssm
            .get(i - 1)
            .and_then(|row| row.get(b[j - 1] as usize))
            .copied(),
        None => matrix
            .get(a[i - 1] as usize)
            .and_then(|row| row.get(b[j - 1] as usize))
            .copied(),
    }
}

/// Rust ownership equivalent of NCBI static `s_GetTraceback`
/// (`blast_sw.c:278`).
///
/// The C routine reconstructs one local alignment from the traceback byte
/// array, creates a `BlastHSP`, applies identity/length filters, adjusts the
/// subject offset, then saves into a `BlastHSPList`. This port keeps the exact
/// traceback walk and coordinate/script normalization, and returns the
/// recovered alignment for the caller to filter/save.
pub fn s_get_traceback(
    trace: &[u8],
    a: &[u8],
    b: &[u8],
    matrix: &[[i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE],
    pssm: Option<&[Vec<i32>]>,
    gap_open: i32,
    gap_extend: i32,
    a_end: usize,
    b_end: usize,
    best_score: i32,
    swapped: bool,
    start_shift: usize,
) -> Option<SmithWatermanTracebackHit> {
    let stride = b.len() + 1;
    if best_score <= 0
        || a_end > a.len()
        || b_end > b.len()
        || trace.len() < (a.len() + 1).checked_mul(stride)?
    {
        return None;
    }

    let mut i = a_end;
    let mut j = b_end;
    let mut script = trace[i * stride + j] & EDIT_OP_MASK;
    let mut curr_score = -best_score;
    let mut prelim = Vec::new();

    while curr_score != 0 {
        let op = blast_sw_traceback_op(script)?;
        prelim.push(op);

        let next_action = trace[i * stride + j];
        match script {
            EDIT_SUB => {
                curr_score += blast_sw_score(a, b, matrix, pssm, i, j)?;
                i = i.checked_sub(1)?;
                j = j.checked_sub(1)?;
                script = trace[i * stride + j] & EDIT_OP_MASK;
            }
            EDIT_GAP_IN_A => {
                j = j.checked_sub(1)?;
                if next_action & EDIT_START_GAP_A != 0 {
                    script = trace[i * stride + j] & EDIT_OP_MASK;
                    curr_score -= gap_open;
                }
                curr_score -= gap_extend;
            }
            EDIT_GAP_IN_B => {
                i = i.checked_sub(1)?;
                if next_action & EDIT_START_GAP_B != 0 {
                    script = trace[i * stride + j] & EDIT_OP_MASK;
                    curr_score -= gap_open;
                }
                curr_score -= gap_extend;
            }
            _ => return None,
        }
    }

    let mut a_start = i;
    let mut b_start = j;
    let mut a_end = a_end;
    let mut b_end = b_end;

    let mut edit_script = GapEditScript::with_capacity(prelim.len());
    for op in prelim.into_iter().rev() {
        let op = if swapped {
            match op {
                GapAlignOpType::Ins => GapAlignOpType::Del,
                GapAlignOpType::Del => GapAlignOpType::Ins,
                _ => op,
            }
        } else {
            op
        };
        edit_script.push(op, 1);
    }

    if swapped {
        std::mem::swap(&mut a_start, &mut b_start);
        std::mem::swap(&mut a_end, &mut b_end);
    }

    Some(SmithWatermanTracebackHit {
        score: best_score,
        a_start,
        a_end,
        b_start: b_start + start_shift,
        b_end: b_end + start_shift,
        edit_script,
    })
}

/// 1-1 port of `SmithWatermanScoreWithTraceback` from `blast_sw.c`.
///
/// C saves accepted HSPs into a `BlastHSPList`; this Rust translation returns
/// each recovered alignment so the surrounding translated traceback layer can
/// apply its own hit-saving policy.
pub fn smith_waterman_score_with_traceback(
    a: &[u8],
    b: &[u8],
    matrix: &[[i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE],
    pssm: Option<&[Vec<i32>]>,
    gap_open: i32,
    gap_extend: i32,
    start_shift: usize,
    cutoff: i32,
) -> Vec<SmithWatermanTracebackHit> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }

    let mut swapped = false;
    let (a_work, b_work) = if pssm.is_none() && a.len() < b.len() {
        swapped = true;
        (b, a)
    } else {
        (a, b)
    };

    let gap_open_extend = gap_open + gap_extend;
    let mut scores = vec![BlastGapSwTraceDp::default(); b_work.len() + 1];
    let stride = b_work.len() + 1;
    let mut traceback_array = vec![0u8; (a_work.len() + 1) * stride];
    for cell in &mut traceback_array[..stride] {
        *cell = EDIT_GAP_IN_A;
    }

    let mut hits = Vec::new();

    for i in 1..=a_work.len() {
        let mut insert_score = 0;
        let mut row_score = 0;
        let mut row_path_stop_i = 0usize;
        let mut row_path_stop_j = 0usize;
        let mut row_path_score = 0;
        let row_offset = i * stride;
        traceback_array[row_offset] = EDIT_GAP_IN_B;

        for j in 1..=b_work.len() {
            let mut best_score = scores[j].best_gap - gap_extend;
            let mut script = 0u8;
            if scores[j].best - gap_open_extend > best_score {
                script |= EDIT_START_GAP_B;
                best_score = scores[j].best - gap_open_extend;
            }
            scores[j].best_gap = best_score;

            best_score = insert_score - gap_extend;
            if row_score - gap_open_extend > best_score {
                script |= EDIT_START_GAP_A;
                best_score = row_score - gap_open_extend;
            }
            insert_score = best_score;

            let sub_score = blast_sw_score(a_work, b_work, matrix, pssm, i, j)
                .unwrap_or(crate::stat::MININT / 2);
            best_score = (scores[j - 1].best + sub_score).max(0);
            traceback_array[row_offset + j] = script | EDIT_SUB;
            let mut new_path_score = scores[j - 1].path_score;
            let mut new_path_stop_i = scores[j - 1].path_stop_i;
            let mut new_path_stop_j = scores[j - 1].path_stop_j;

            if insert_score > best_score {
                best_score = insert_score;
                traceback_array[row_offset + j] = script | EDIT_GAP_IN_A;
                new_path_score = row_path_score;
                new_path_stop_i = row_path_stop_i;
                new_path_stop_j = row_path_stop_j;
            }
            if scores[j].best_gap >= best_score {
                best_score = scores[j].best_gap;
                traceback_array[row_offset + j] = script | EDIT_GAP_IN_B;
                new_path_score = scores[j].path_score;
                new_path_stop_i = scores[j].path_stop_i;
                new_path_stop_j = scores[j].path_stop_j;
            }

            if best_score == 0 {
                if new_path_score >= cutoff {
                    if let Some(hit) = s_get_traceback(
                        &traceback_array,
                        a_work,
                        b_work,
                        matrix,
                        pssm,
                        gap_open,
                        gap_extend,
                        new_path_stop_i,
                        new_path_stop_j,
                        new_path_score,
                        swapped,
                        start_shift,
                    ) {
                        hits.push(hit);
                    }
                }
                new_path_score = 0;
            }

            if best_score > new_path_score {
                new_path_score = best_score;
                new_path_stop_i = i;
                new_path_stop_j = j;
            }

            scores[j - 1].best = row_score;
            scores[j - 1].path_score = row_path_score;
            scores[j - 1].path_stop_i = row_path_stop_i;
            scores[j - 1].path_stop_j = row_path_stop_j;

            row_score = best_score;
            row_path_score = new_path_score;
            row_path_stop_i = new_path_stop_i;
            row_path_stop_j = new_path_stop_j;
        }

        let last = b_work.len();
        scores[last].best = row_score;
        scores[last].path_score = row_path_score;
        scores[last].path_stop_i = row_path_stop_i;
        scores[last].path_stop_j = row_path_stop_j;

        if scores[last].path_score >= cutoff {
            if let Some(hit) = s_get_traceback(
                &traceback_array,
                a_work,
                b_work,
                matrix,
                pssm,
                gap_open,
                gap_extend,
                scores[last].path_stop_i,
                scores[last].path_stop_j,
                scores[last].path_score,
                swapped,
                start_shift,
            ) {
                hits.push(hit);
            }
        }
    }

    for i in 0..b_work.len() {
        if scores[i].best != 0 && scores[i].path_score >= cutoff {
            if let Some(hit) = s_get_traceback(
                &traceback_array,
                a_work,
                b_work,
                matrix,
                pssm,
                gap_open,
                gap_extend,
                scores[i].path_stop_i,
                scores[i].path_stop_j,
                scores[i].path_score,
                swapped,
                start_shift,
            ) {
                hits.push(hit);
            }
        }
    }

    hits
}

/// 1-1 port of `s_NuclSmithWaterman` from `blast_sw.c`.
///
/// `b_packed` is the packed NCBI2na sequence used by BLAST databases, while
/// `a` is BLASTNA/2-bit query data. The outer loop unpacks `B` exactly where
/// C uses `NCBI2NA_UNPACK_BASE`.
pub fn s_nucl_smith_waterman(
    b_packed: &[u8],
    b_size: usize,
    a: &[u8],
    matrix: &[[i32; crate::encoding::BLASTNA_SIZE]; crate::encoding::BLASTNA_SIZE],
    gap_open: i32,
    gap_extend: i32,
) -> i32 {
    if b_size == 0 || a.is_empty() {
        return 0;
    }

    let gap_open_extend = gap_open + gap_extend;
    let mut scores = vec![BlastSwGapDp::default(); a.len() + 1];
    let mut final_best_score = 0;

    for i in 1..=b_size {
        let base_pair = crate::encoding::ncbi2na_base_at(b_packed, i - 1) as usize;
        let matrix_row = &matrix[base_pair];
        let mut insert_score = 0;
        let mut row_score = 0;

        for j in 1..=a.len() {
            let mut best_score = scores[j].best_gap - gap_extend;
            if scores[j].best - gap_open_extend > best_score {
                best_score = scores[j].best - gap_open_extend;
            }
            scores[j].best_gap = best_score;

            best_score = insert_score - gap_extend;
            if row_score - gap_open_extend > best_score {
                best_score = row_score - gap_open_extend;
            }
            insert_score = best_score;

            best_score = (scores[j - 1].best + matrix_row[a[j - 1] as usize]).max(0);
            if insert_score > best_score {
                best_score = insert_score;
            }
            if scores[j].best_gap > best_score {
                best_score = scores[j].best_gap;
            }

            if best_score > final_best_score {
                final_best_score = best_score;
            }

            scores[j - 1].best = row_score;
            row_score = best_score;
        }

        scores[a.len()].best = row_score;
    }

    final_best_score
}

/// blast-rs: Rust-spelled wrapper for the C BLbasic Smith-Waterman
/// score-only primitive; not tagged as a direct NCBI C port because the
/// provenance audit mechanically expands `BLbasic` as `b_lbasic`.
///
/// Returns `(score, match_seq_end, query_end)` of the locally optimal
/// alignment ending at the highest-scoring cell. No traceback.
///
/// `matrix[query_byte][subject_byte]` indexed by encoded sequence bytes
/// for protein (NCBIstdaa). `gap_open` and `gap_extend` are the *positive*
/// gap costs (NCBI's convention is `-gapOpen`, `-gapExtend`); pass the
/// magnitudes here.
pub fn bl_basic_smith_waterman_score_only(
    match_seq: &[u8],
    query: &[u8],
    matrix: &[[i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE],
    gap_open: i32,
    gap_extend: i32,
) -> (i32, usize, usize) {
    let match_seq_length = match_seq.len();
    let query_length = query.len();
    if match_seq_length == 0 || query_length == 0 {
        return (0, 0, 0);
    }

    let mut score_vector: Vec<SwGapInfo> = (0..match_seq_length)
        .map(|_| SwGapInfo {
            no_gap: 0,
            gap_exists: -gap_open,
        })
        .collect();

    let mut best_score = 0i32;
    let mut best_match_seq_pos = 0usize;
    let mut best_query_pos = 0usize;
    let new_gap_cost = gap_open + gap_extend;

    for query_pos in 0..query_length {
        let matrix_row = &matrix[query[query_pos] as usize];
        let mut new_score = 0i32;
        let mut prev_score_no_gap_match_seq = 0i32;
        let mut prev_score_gap_match_seq = -gap_open;

        for match_seq_pos in 0..match_seq_length {
            // Gap in match_seq: start new or extend existing.
            new_score -= new_gap_cost;
            prev_score_gap_match_seq -= gap_extend;
            if new_score > prev_score_gap_match_seq {
                prev_score_gap_match_seq = new_score;
            }

            // Gap in query: start new or extend existing.
            new_score = score_vector[match_seq_pos].no_gap - new_gap_cost;
            let mut continue_gap_score = score_vector[match_seq_pos].gap_exists - gap_extend;
            if new_score > continue_gap_score {
                continue_gap_score = new_score;
            }

            // Substitution: extend by one position in match_seq + query.
            new_score = prev_score_no_gap_match_seq + matrix_row[match_seq[match_seq_pos] as usize];
            if new_score < 0 {
                new_score = 0; // Smith-Waterman locality
            }
            if new_score < prev_score_gap_match_seq {
                new_score = prev_score_gap_match_seq;
            }
            if new_score < continue_gap_score {
                new_score = continue_gap_score;
            }

            prev_score_no_gap_match_seq = score_vector[match_seq_pos].no_gap;
            score_vector[match_seq_pos].no_gap = new_score;
            score_vector[match_seq_pos].gap_exists = continue_gap_score;

            if new_score > best_score {
                best_score = new_score;
                best_query_pos = query_pos;
                best_match_seq_pos = match_seq_pos;
            }
        }
    }

    if best_score < 0 {
        best_score = 0;
    }

    (best_score, best_match_seq_pos, best_query_pos)
}

/// Reverse SW from a given endpoint — 1-1 port of `BLSmithWatermanFindStart`
/// (`smith_waterman.c:144`).
///
/// Walks the DP backward from `(query_end, match_seq_end)` until the
/// running best score reaches `score_in` (the score returned by the
/// forward score-only pass). Returns `(score, match_seq_start, query_start)`.
///
/// Same matrix/gap convention as `bl_basic_smith_waterman_score_only`.
pub fn bl_smith_waterman_find_start(
    match_seq: &[u8],
    query: &[u8],
    matrix: &[[i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE],
    gap_open: i32,
    gap_extend: i32,
    match_seq_end: usize,
    query_end: usize,
    score_in: i32,
) -> (i32, usize, usize) {
    let match_seq_length = match_seq_end + 1;
    if match_seq_length == 0 || query_end == usize::MAX {
        return (0, 0, 0);
    }

    let mut score_vector: Vec<SwGapInfo> = (0..match_seq_length)
        .map(|_| SwGapInfo {
            no_gap: 0,
            gap_exists: -gap_open,
        })
        .collect();

    let mut best_score = 0i32;
    let mut best_match_seq_pos = 0usize;
    let mut best_query_pos = 0usize;
    let new_gap_cost = gap_open + gap_extend;

    'outer: for query_pos in (0..=query_end).rev() {
        let matrix_row = &matrix[query[query_pos] as usize];
        let mut new_score = 0i32;
        let mut prev_score_no_gap_match_seq = 0i32;
        let mut prev_score_gap_match_seq = -gap_open;

        for match_seq_pos in (0..=match_seq_end).rev() {
            new_score -= new_gap_cost;
            prev_score_gap_match_seq -= gap_extend;
            if new_score > prev_score_gap_match_seq {
                prev_score_gap_match_seq = new_score;
            }

            new_score = score_vector[match_seq_pos].no_gap - new_gap_cost;
            let mut continue_gap_score = score_vector[match_seq_pos].gap_exists - gap_extend;
            if new_score > continue_gap_score {
                continue_gap_score = new_score;
            }

            new_score = prev_score_no_gap_match_seq + matrix_row[match_seq[match_seq_pos] as usize];
            if new_score < 0 {
                new_score = 0;
            }
            if new_score < prev_score_gap_match_seq {
                new_score = prev_score_gap_match_seq;
            }
            if new_score < continue_gap_score {
                new_score = continue_gap_score;
            }

            prev_score_no_gap_match_seq = score_vector[match_seq_pos].no_gap;
            score_vector[match_seq_pos].no_gap = new_score;
            score_vector[match_seq_pos].gap_exists = continue_gap_score;

            if new_score > best_score {
                best_score = new_score;
                best_query_pos = query_pos;
                best_match_seq_pos = match_seq_pos;
            }
            if best_score >= score_in {
                break 'outer;
            }
        }
    }

    if best_score < 0 {
        best_score = 0;
    }

    (best_score, best_match_seq_pos, best_query_pos)
}

/// Sentinel matching NCBI's `COMPO_SCORE_MIN`
/// (`composition_constants.h`). When a forbidden cell is encountered,
/// the substitution score is replaced by this value so the SW DP can
/// never claim a positive contribution from that position.
const COMPO_SCORE_MIN: i32 = i32::MIN / 2;

fn pssm_cell_score(pssm: &[Vec<i32>], query_offset: usize, query_pos: usize, subject: u8) -> i32 {
    pssm[query_offset + query_pos][subject as usize]
}

/// PSSM-backed variant of [`bl_basic_smith_waterman_score_only`].
///
/// `pssm[query_offset + query_pos][subject_byte]` supplies the
/// substitution score. The DP state and gap/local SW behavior are kept
/// identical to the square-matrix primitive.
pub fn bl_basic_smith_waterman_score_only_pssm(
    match_seq: &[u8],
    query: &[u8],
    pssm: &[Vec<i32>],
    query_offset: usize,
    gap_open: i32,
    gap_extend: i32,
) -> (i32, usize, usize) {
    let match_seq_length = match_seq.len();
    let query_length = query.len();
    if match_seq_length == 0 || query_length == 0 {
        return (0, 0, 0);
    }

    let mut score_vector: Vec<SwGapInfo> = (0..match_seq_length)
        .map(|_| SwGapInfo {
            no_gap: 0,
            gap_exists: -gap_open,
        })
        .collect();

    let mut best_score = 0i32;
    let mut best_match_seq_pos = 0usize;
    let mut best_query_pos = 0usize;
    let new_gap_cost = gap_open + gap_extend;

    for query_pos in 0..query_length {
        let mut new_score = 0i32;
        let mut prev_score_no_gap_match_seq = 0i32;
        let mut prev_score_gap_match_seq = -gap_open;

        for match_seq_pos in 0..match_seq_length {
            new_score -= new_gap_cost;
            prev_score_gap_match_seq -= gap_extend;
            if new_score > prev_score_gap_match_seq {
                prev_score_gap_match_seq = new_score;
            }

            new_score = score_vector[match_seq_pos].no_gap - new_gap_cost;
            let mut continue_gap_score = score_vector[match_seq_pos].gap_exists - gap_extend;
            if new_score > continue_gap_score {
                continue_gap_score = new_score;
            }

            new_score = prev_score_no_gap_match_seq
                + pssm_cell_score(pssm, query_offset, query_pos, match_seq[match_seq_pos]);
            if new_score < 0 {
                new_score = 0;
            }
            if new_score < prev_score_gap_match_seq {
                new_score = prev_score_gap_match_seq;
            }
            if new_score < continue_gap_score {
                new_score = continue_gap_score;
            }

            prev_score_no_gap_match_seq = score_vector[match_seq_pos].no_gap;
            score_vector[match_seq_pos].no_gap = new_score;
            score_vector[match_seq_pos].gap_exists = continue_gap_score;

            if new_score > best_score {
                best_score = new_score;
                best_query_pos = query_pos;
                best_match_seq_pos = match_seq_pos;
            }
        }
    }

    if best_score < 0 {
        best_score = 0;
    }

    (best_score, best_match_seq_pos, best_query_pos)
}

/// PSSM-backed variant of [`bl_smith_waterman_find_start`].
pub fn bl_smith_waterman_find_start_pssm(
    match_seq: &[u8],
    _query: &[u8],
    pssm: &[Vec<i32>],
    query_offset: usize,
    gap_open: i32,
    gap_extend: i32,
    match_seq_end: usize,
    query_end: usize,
    score_in: i32,
) -> (i32, usize, usize) {
    let match_seq_length = match_seq_end + 1;
    if match_seq_length == 0 || query_end == usize::MAX {
        return (0, 0, 0);
    }

    let mut score_vector: Vec<SwGapInfo> = (0..match_seq_length)
        .map(|_| SwGapInfo {
            no_gap: 0,
            gap_exists: -gap_open,
        })
        .collect();

    let mut best_score = 0i32;
    let mut best_match_seq_pos = 0usize;
    let mut best_query_pos = 0usize;
    let new_gap_cost = gap_open + gap_extend;

    'outer: for query_pos in (0..=query_end).rev() {
        let mut new_score = 0i32;
        let mut prev_score_no_gap_match_seq = 0i32;
        let mut prev_score_gap_match_seq = -gap_open;

        for match_seq_pos in (0..=match_seq_end).rev() {
            new_score -= new_gap_cost;
            prev_score_gap_match_seq -= gap_extend;
            if new_score > prev_score_gap_match_seq {
                prev_score_gap_match_seq = new_score;
            }

            new_score = score_vector[match_seq_pos].no_gap - new_gap_cost;
            let mut continue_gap_score = score_vector[match_seq_pos].gap_exists - gap_extend;
            if new_score > continue_gap_score {
                continue_gap_score = new_score;
            }

            new_score = prev_score_no_gap_match_seq
                + pssm_cell_score(pssm, query_offset, query_pos, match_seq[match_seq_pos]);
            if new_score < 0 {
                new_score = 0;
            }
            if new_score < prev_score_gap_match_seq {
                new_score = prev_score_gap_match_seq;
            }
            if new_score < continue_gap_score {
                new_score = continue_gap_score;
            }

            prev_score_no_gap_match_seq = score_vector[match_seq_pos].no_gap;
            score_vector[match_seq_pos].no_gap = new_score;
            score_vector[match_seq_pos].gap_exists = continue_gap_score;

            if new_score > best_score {
                best_score = new_score;
                best_query_pos = query_pos;
                best_match_seq_pos = match_seq_pos;
            }
            if best_score >= score_in {
                break 'outer;
            }
        }
    }

    if best_score < 0 {
        best_score = 0;
    }

    (best_score, best_match_seq_pos, best_query_pos)
}

/// 1-1 port of `Blast_ForbiddenRanges` (`smith_waterman.h`).
///
/// Tracks per-query-position lists of forbidden subject ranges, used by
/// the redo-alignment driver to suppress re-finding the same SW
/// alignment when multiple HSPs come from the same (query, subject)
/// pair. Each query position has a flat `Vec<i32>` of `[start, end,
/// start, end, …]` pairs — same memory layout as NCBI's
/// `int** ranges`.
#[derive(Debug, Clone, Default)]
pub struct BlastForbiddenRanges {
    /// `numForbidden[query_pos]` — number of forbidden ranges at this
    /// query position (each range is a `(start, end)` pair).
    pub num_forbidden: Vec<i32>,
    /// `ranges[query_pos]` — flat array of `[start0, end0, start1,
    /// end1, …]` for each forbidden range at this query position.
    pub ranges: Vec<Vec<i32>>,
    /// `isEmpty` — fast-path flag mirrored from NCBI; TRUE means no
    /// ranges have been pushed yet.
    pub is_empty: bool,
}

impl BlastForbiddenRanges {
    /// blast-rs: Associated constructor wrapper around the canonical
    /// forbidden-ranges initializer.
    ///
    /// `capacity` is the concatenated query length; we allocate one
    /// entry per query position. NCBI initializes `ranges[f]` to a
    /// 2-int slot (`[0, 0]`); we use empty `Vec`s instead since
    /// `num_forbidden[f] == 0` already signals "nothing here".
    pub fn new(capacity: i32) -> Self {
        blast_forbidden_ranges_initialize(capacity)
    }

    /// blast-rs: Associated method wrapper around the canonical
    /// forbidden-ranges clear function.
    pub fn clear(&mut self) {
        blast_forbidden_ranges_clear(self);
    }

    /// blast-rs: Associated method wrapper around the canonical
    /// forbidden-ranges push function.
    pub fn push(&mut self, query_start: i32, query_end: i32, match_start: i32, match_end: i32) {
        blast_forbidden_ranges_push(self, query_start, query_end, match_start, match_end);
    }
}

/// NCBI: Blast_ForbiddenRangesInitialize (smith_waterman.c:473).
///
/// `capacity` is the concatenated query length; we allocate one
/// entry per query position. NCBI initializes `ranges[f]` to a
/// 2-int slot (`[0, 0]`); we use empty `Vec`s instead since
/// `num_forbidden[f] == 0` already signals "nothing here".
pub fn blast_forbidden_ranges_initialize(capacity: i32) -> BlastForbiddenRanges {
    let cap = capacity.max(0) as usize;
    BlastForbiddenRanges {
        num_forbidden: vec![0i32; cap],
        ranges: vec![Vec::new(); cap],
        is_empty: true,
    }
}

/// NCBI: Blast_ForbiddenRangesClear (smith_waterman.c:505).
///
/// Resets all per-position counts to zero and flips `is_empty` back to true.
/// Backing arrays are kept around (matching NCBI).
pub fn blast_forbidden_ranges_clear(forbidden: &mut BlastForbiddenRanges) {
    for n in forbidden.num_forbidden.iter_mut() {
        *n = 0;
    }
    // NCBI doesn't truncate `ranges[f]` either — the reused backing storage is
    // harmless because `numForbidden[f] == 0` makes the inner loop in
    // `BLspecial*` skip them entirely.
    forbidden.is_empty = true;
}

/// NCBI: Blast_ForbiddenRangesPush (smith_waterman.c:516).
///
/// Marks the rectangular region `[query_start, query_end) x
/// [match_start, match_end]` forbidden. For every query position in
/// `[query_start, query_end)`, append one `(matchStart, matchEnd)` pair.
pub fn blast_forbidden_ranges_push(
    forbidden: &mut BlastForbiddenRanges,
    query_start: i32,
    query_end: i32,
    match_start: i32,
    match_end: i32,
) {
    let qs = query_start.max(0) as usize;
    let qe = query_end.max(0) as usize;
    if qe > forbidden.num_forbidden.len() {
        forbidden.num_forbidden.resize(qe, 0);
        forbidden.ranges.resize(qe, Vec::new());
    }
    for f in qs..qe {
        forbidden.ranges[f].push(match_start);
        forbidden.ranges[f].push(match_end);
        forbidden.num_forbidden[f] += 1;
    }
    forbidden.is_empty = false;
}

/// blast-rs: Rust-spelled wrapper for the C BLspecial Smith-Waterman
/// score-only primitive; not tagged as a direct NCBI C port because the
/// provenance audit mechanically expands `BLspecial` as `b_lspecial`.
///
/// Forward score-only SW with forbidden-range exclusion. For each cell
/// `(query_pos, match_pos)` we check whether `match_pos` falls in any
/// of the forbidden ranges registered at that query position; if so,
/// the substitution contribution is replaced by `COMPO_SCORE_MIN` —
/// matching NCBI's `forbidden ? COMPO_SCORE_MIN : score`.
pub fn bl_special_smith_waterman_score_only(
    match_seq: &[u8],
    query: &[u8],
    matrix: &[[i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE],
    gap_open: i32,
    gap_extend: i32,
    forbidden: &BlastForbiddenRanges,
) -> (i32, usize, usize) {
    let match_seq_length = match_seq.len();
    let query_length = query.len();
    if match_seq_length == 0 || query_length == 0 {
        return (0, 0, 0);
    }

    let mut score_vector: Vec<SwGapInfo> = (0..match_seq_length)
        .map(|_| SwGapInfo {
            no_gap: 0,
            gap_exists: -gap_open,
        })
        .collect();

    let mut best_score = 0i32;
    let mut best_match_seq_pos = 0usize;
    let mut best_query_pos = 0usize;
    let new_gap_cost = gap_open + gap_extend;

    for query_pos in 0..query_length {
        let matrix_row = &matrix[query[query_pos] as usize];
        let mut new_score = 0i32;
        let mut prev_score_no_gap_match_seq = 0i32;
        let mut prev_score_gap_match_seq = -gap_open;

        let nf = *forbidden.num_forbidden.get(query_pos).unwrap_or(&0) as usize;
        let row_ranges: &[i32] = forbidden
            .ranges
            .get(query_pos)
            .map(|v| v.as_slice())
            .unwrap_or(&[]);

        for match_seq_pos in 0..match_seq_length {
            new_score -= new_gap_cost;
            prev_score_gap_match_seq -= gap_extend;
            if new_score > prev_score_gap_match_seq {
                prev_score_gap_match_seq = new_score;
            }
            new_score = score_vector[match_seq_pos].no_gap - new_gap_cost;
            let mut continue_gap_score = score_vector[match_seq_pos].gap_exists - gap_extend;
            if new_score > continue_gap_score {
                continue_gap_score = new_score;
            }

            let mut is_forbidden = false;
            for f in 0..nf {
                let r0 = row_ranges[2 * f] as i64;
                let r1 = row_ranges[2 * f + 1] as i64;
                if (match_seq_pos as i64) >= r0 && (match_seq_pos as i64) <= r1 {
                    is_forbidden = true;
                    break;
                }
            }
            new_score = if is_forbidden {
                COMPO_SCORE_MIN
            } else {
                prev_score_no_gap_match_seq + matrix_row[match_seq[match_seq_pos] as usize]
            };
            if new_score < 0 {
                new_score = 0;
            }
            if new_score < prev_score_gap_match_seq {
                new_score = prev_score_gap_match_seq;
            }
            if new_score < continue_gap_score {
                new_score = continue_gap_score;
            }

            prev_score_no_gap_match_seq = score_vector[match_seq_pos].no_gap;
            score_vector[match_seq_pos].no_gap = new_score;
            score_vector[match_seq_pos].gap_exists = continue_gap_score;

            if new_score > best_score {
                best_score = new_score;
                best_query_pos = query_pos;
                best_match_seq_pos = match_seq_pos;
            }
        }
    }

    if best_score < 0 {
        best_score = 0;
    }
    (best_score, best_match_seq_pos, best_query_pos)
}

/// blast-rs: Rust-spelled wrapper for the C BLspecial Smith-Waterman
/// find-start primitive; not tagged as a direct NCBI C port because the
/// provenance audit mechanically expands `BLspecial` as `b_lspecial`.
///
/// Reverse SW from `(match_seq_end, query_end)` with forbidden-range
/// exclusion and an early-exit when `bestScore >= score_in` is reached.
pub fn bl_special_smith_waterman_find_start(
    match_seq: &[u8],
    query: &[u8],
    matrix: &[[i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE],
    gap_open: i32,
    gap_extend: i32,
    match_seq_end: usize,
    query_end: usize,
    score_in: i32,
    forbidden: &BlastForbiddenRanges,
) -> (i32, usize, usize) {
    let match_seq_length = match_seq_end + 1;
    if match_seq_length == 0 {
        return (0, 0, 0);
    }

    let mut score_vector: Vec<SwGapInfo> = (0..match_seq_length)
        .map(|_| SwGapInfo {
            no_gap: 0,
            gap_exists: -gap_open,
        })
        .collect();

    let mut best_score = 0i32;
    let mut best_match_seq_pos = 0usize;
    let mut best_query_pos = 0usize;
    let new_gap_cost = gap_open + gap_extend;

    'outer: for query_pos in (0..=query_end).rev() {
        let matrix_row = &matrix[query[query_pos] as usize];
        let mut new_score = 0i32;
        let mut prev_score_no_gap_match_seq = 0i32;
        let mut prev_score_gap_match_seq = -gap_open;

        let nf = *forbidden.num_forbidden.get(query_pos).unwrap_or(&0) as usize;
        let row_ranges: &[i32] = forbidden
            .ranges
            .get(query_pos)
            .map(|v| v.as_slice())
            .unwrap_or(&[]);

        for match_seq_pos in (0..=match_seq_end).rev() {
            new_score -= new_gap_cost;
            prev_score_gap_match_seq -= gap_extend;
            if new_score > prev_score_gap_match_seq {
                prev_score_gap_match_seq = new_score;
            }
            new_score = score_vector[match_seq_pos].no_gap - new_gap_cost;
            let mut continue_gap_score = score_vector[match_seq_pos].gap_exists - gap_extend;
            if new_score > continue_gap_score {
                continue_gap_score = new_score;
            }

            let mut is_forbidden = false;
            for f in 0..nf {
                let r0 = row_ranges[2 * f] as i64;
                let r1 = row_ranges[2 * f + 1] as i64;
                if (match_seq_pos as i64) >= r0 && (match_seq_pos as i64) <= r1 {
                    is_forbidden = true;
                    break;
                }
            }
            new_score = if is_forbidden {
                COMPO_SCORE_MIN
            } else {
                prev_score_no_gap_match_seq + matrix_row[match_seq[match_seq_pos] as usize]
            };
            if new_score < 0 {
                new_score = 0;
            }
            if new_score < prev_score_gap_match_seq {
                new_score = prev_score_gap_match_seq;
            }
            if new_score < continue_gap_score {
                new_score = continue_gap_score;
            }

            prev_score_no_gap_match_seq = score_vector[match_seq_pos].no_gap;
            score_vector[match_seq_pos].no_gap = new_score;
            score_vector[match_seq_pos].gap_exists = continue_gap_score;

            if new_score > best_score {
                best_score = new_score;
                best_query_pos = query_pos;
                best_match_seq_pos = match_seq_pos;
            }
            if best_score >= score_in {
                break 'outer;
            }
        }
    }

    if best_score < 0 {
        best_score = 0;
    }
    (best_score, best_match_seq_pos, best_query_pos)
}

/// PSSM-backed variant of [`bl_special_smith_waterman_score_only`].
pub fn bl_special_smith_waterman_score_only_pssm(
    match_seq: &[u8],
    query: &[u8],
    pssm: &[Vec<i32>],
    query_offset: usize,
    gap_open: i32,
    gap_extend: i32,
    forbidden: &BlastForbiddenRanges,
) -> (i32, usize, usize) {
    let match_seq_length = match_seq.len();
    let query_length = query.len();
    if match_seq_length == 0 || query_length == 0 {
        return (0, 0, 0);
    }

    let mut score_vector: Vec<SwGapInfo> = (0..match_seq_length)
        .map(|_| SwGapInfo {
            no_gap: 0,
            gap_exists: -gap_open,
        })
        .collect();

    let mut best_score = 0i32;
    let mut best_match_seq_pos = 0usize;
    let mut best_query_pos = 0usize;
    let new_gap_cost = gap_open + gap_extend;

    for query_pos in 0..query_length {
        let mut new_score = 0i32;
        let mut prev_score_no_gap_match_seq = 0i32;
        let mut prev_score_gap_match_seq = -gap_open;

        let nf = *forbidden.num_forbidden.get(query_pos).unwrap_or(&0) as usize;
        let row_ranges: &[i32] = forbidden
            .ranges
            .get(query_pos)
            .map(|v| v.as_slice())
            .unwrap_or(&[]);

        for match_seq_pos in 0..match_seq_length {
            new_score -= new_gap_cost;
            prev_score_gap_match_seq -= gap_extend;
            if new_score > prev_score_gap_match_seq {
                prev_score_gap_match_seq = new_score;
            }
            new_score = score_vector[match_seq_pos].no_gap - new_gap_cost;
            let mut continue_gap_score = score_vector[match_seq_pos].gap_exists - gap_extend;
            if new_score > continue_gap_score {
                continue_gap_score = new_score;
            }

            let mut is_forbidden = false;
            for f in 0..nf {
                let r0 = row_ranges[2 * f] as i64;
                let r1 = row_ranges[2 * f + 1] as i64;
                if (match_seq_pos as i64) >= r0 && (match_seq_pos as i64) <= r1 {
                    is_forbidden = true;
                    break;
                }
            }
            new_score = if is_forbidden {
                COMPO_SCORE_MIN
            } else {
                prev_score_no_gap_match_seq
                    + pssm_cell_score(pssm, query_offset, query_pos, match_seq[match_seq_pos])
            };
            if new_score < 0 {
                new_score = 0;
            }
            if new_score < prev_score_gap_match_seq {
                new_score = prev_score_gap_match_seq;
            }
            if new_score < continue_gap_score {
                new_score = continue_gap_score;
            }

            prev_score_no_gap_match_seq = score_vector[match_seq_pos].no_gap;
            score_vector[match_seq_pos].no_gap = new_score;
            score_vector[match_seq_pos].gap_exists = continue_gap_score;

            if new_score > best_score {
                best_score = new_score;
                best_query_pos = query_pos;
                best_match_seq_pos = match_seq_pos;
            }
        }
    }

    if best_score < 0 {
        best_score = 0;
    }
    (best_score, best_match_seq_pos, best_query_pos)
}

/// PSSM-backed variant of [`bl_special_smith_waterman_find_start`].
pub fn bl_special_smith_waterman_find_start_pssm(
    match_seq: &[u8],
    _query: &[u8],
    pssm: &[Vec<i32>],
    query_offset: usize,
    gap_open: i32,
    gap_extend: i32,
    match_seq_end: usize,
    query_end: usize,
    score_in: i32,
    forbidden: &BlastForbiddenRanges,
) -> (i32, usize, usize) {
    let match_seq_length = match_seq_end + 1;
    if match_seq_length == 0 {
        return (0, 0, 0);
    }

    let mut score_vector: Vec<SwGapInfo> = (0..match_seq_length)
        .map(|_| SwGapInfo {
            no_gap: 0,
            gap_exists: -gap_open,
        })
        .collect();

    let mut best_score = 0i32;
    let mut best_match_seq_pos = 0usize;
    let mut best_query_pos = 0usize;
    let new_gap_cost = gap_open + gap_extend;

    'outer: for query_pos in (0..=query_end).rev() {
        let mut new_score = 0i32;
        let mut prev_score_no_gap_match_seq = 0i32;
        let mut prev_score_gap_match_seq = -gap_open;

        let nf = *forbidden.num_forbidden.get(query_pos).unwrap_or(&0) as usize;
        let row_ranges: &[i32] = forbidden
            .ranges
            .get(query_pos)
            .map(|v| v.as_slice())
            .unwrap_or(&[]);

        for match_seq_pos in (0..=match_seq_end).rev() {
            new_score -= new_gap_cost;
            prev_score_gap_match_seq -= gap_extend;
            if new_score > prev_score_gap_match_seq {
                prev_score_gap_match_seq = new_score;
            }
            new_score = score_vector[match_seq_pos].no_gap - new_gap_cost;
            let mut continue_gap_score = score_vector[match_seq_pos].gap_exists - gap_extend;
            if new_score > continue_gap_score {
                continue_gap_score = new_score;
            }

            let mut is_forbidden = false;
            for f in 0..nf {
                let r0 = row_ranges[2 * f] as i64;
                let r1 = row_ranges[2 * f + 1] as i64;
                if (match_seq_pos as i64) >= r0 && (match_seq_pos as i64) <= r1 {
                    is_forbidden = true;
                    break;
                }
            }
            new_score = if is_forbidden {
                COMPO_SCORE_MIN
            } else {
                prev_score_no_gap_match_seq
                    + pssm_cell_score(pssm, query_offset, query_pos, match_seq[match_seq_pos])
            };
            if new_score < 0 {
                new_score = 0;
            }
            if new_score < prev_score_gap_match_seq {
                new_score = prev_score_gap_match_seq;
            }
            if new_score < continue_gap_score {
                new_score = continue_gap_score;
            }

            prev_score_no_gap_match_seq = score_vector[match_seq_pos].no_gap;
            score_vector[match_seq_pos].no_gap = new_score;
            score_vector[match_seq_pos].gap_exists = continue_gap_score;

            if new_score > best_score {
                best_score = new_score;
                best_query_pos = query_pos;
                best_match_seq_pos = match_seq_pos;
            }
            if best_score >= score_in {
                break 'outer;
            }
        }
    }

    if best_score < 0 {
        best_score = 0;
    }
    (best_score, best_match_seq_pos, best_query_pos)
}

/// 1-1 port of `Blast_SmithWatermanScoreOnly` (`smith_waterman.c:545`).
///
/// Public dispatch — falls through to [`bl_basic_smith_waterman_score_only`]
/// when `forbidden.is_empty`, otherwise calls the special variant. Both
/// variants are mirrors of NCBI's two paths.
pub fn blast_smith_waterman_score_only(
    match_seq: &[u8],
    query: &[u8],
    matrix: &[[i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE],
    gap_open: i32,
    gap_extend: i32,
    forbidden: Option<&BlastForbiddenRanges>,
) -> (i32, usize, usize) {
    if forbidden.is_none_or(|ranges| ranges.is_empty) {
        bl_basic_smith_waterman_score_only(match_seq, query, matrix, gap_open, gap_extend)
    } else {
        bl_special_smith_waterman_score_only(
            match_seq,
            query,
            matrix,
            gap_open,
            gap_extend,
            forbidden.expect("checked non-empty forbidden ranges"),
        )
    }
}

/// blast-rs: compatibility adapter for callers that already have a required
/// forbidden-range reference.
pub fn blast_smith_waterman_score_only_with_forbidden(
    match_seq: &[u8],
    query: &[u8],
    matrix: &[[i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE],
    gap_open: i32,
    gap_extend: i32,
    forbidden: &BlastForbiddenRanges,
) -> (i32, usize, usize) {
    blast_smith_waterman_score_only(
        match_seq,
        query,
        matrix,
        gap_open,
        gap_extend,
        Some(forbidden),
    )
}

/// PSSM-backed public score-only dispatch.
pub fn blast_smith_waterman_score_only_pssm(
    match_seq: &[u8],
    query: &[u8],
    pssm: &[Vec<i32>],
    query_offset: usize,
    gap_open: i32,
    gap_extend: i32,
) -> (i32, usize, usize) {
    bl_basic_smith_waterman_score_only_pssm(
        match_seq,
        query,
        pssm,
        query_offset,
        gap_open,
        gap_extend,
    )
}

/// PSSM-backed score-only dispatch with forbidden ranges.
pub fn blast_smith_waterman_score_only_pssm_with_forbidden(
    match_seq: &[u8],
    query: &[u8],
    pssm: &[Vec<i32>],
    query_offset: usize,
    gap_open: i32,
    gap_extend: i32,
    forbidden: &BlastForbiddenRanges,
) -> (i32, usize, usize) {
    if forbidden.is_empty {
        bl_basic_smith_waterman_score_only_pssm(
            match_seq,
            query,
            pssm,
            query_offset,
            gap_open,
            gap_extend,
        )
    } else {
        bl_special_smith_waterman_score_only_pssm(
            match_seq,
            query,
            pssm,
            query_offset,
            gap_open,
            gap_extend,
            forbidden,
        )
    }
}

/// 1-1 port of `Blast_SmithWatermanFindStart` (`smith_waterman.c:579`).
///
/// Public dispatch — falls through to [`bl_smith_waterman_find_start`]
/// for the empty-forbidden-range case.
pub fn blast_smith_waterman_find_start(
    match_seq: &[u8],
    query: &[u8],
    matrix: &[[i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE],
    gap_open: i32,
    gap_extend: i32,
    match_seq_end: usize,
    query_end: usize,
    score_in: i32,
    forbidden: Option<&BlastForbiddenRanges>,
) -> (i32, usize, usize) {
    if forbidden.is_none_or(|ranges| ranges.is_empty) {
        bl_smith_waterman_find_start(
            match_seq,
            query,
            matrix,
            gap_open,
            gap_extend,
            match_seq_end,
            query_end,
            score_in,
        )
    } else {
        bl_special_smith_waterman_find_start(
            match_seq,
            query,
            matrix,
            gap_open,
            gap_extend,
            match_seq_end,
            query_end,
            score_in,
            forbidden.expect("checked non-empty forbidden ranges"),
        )
    }
}

/// blast-rs: compatibility adapter for callers that already have a required
/// forbidden-range reference.
pub fn blast_smith_waterman_find_start_with_forbidden(
    match_seq: &[u8],
    query: &[u8],
    matrix: &[[i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE],
    gap_open: i32,
    gap_extend: i32,
    match_seq_end: usize,
    query_end: usize,
    score_in: i32,
    forbidden: &BlastForbiddenRanges,
) -> (i32, usize, usize) {
    blast_smith_waterman_find_start(
        match_seq,
        query,
        matrix,
        gap_open,
        gap_extend,
        match_seq_end,
        query_end,
        score_in,
        Some(forbidden),
    )
}

/// PSSM-backed public find-start dispatch.
pub fn blast_smith_waterman_find_start_pssm(
    match_seq: &[u8],
    query: &[u8],
    pssm: &[Vec<i32>],
    query_offset: usize,
    gap_open: i32,
    gap_extend: i32,
    match_seq_end: usize,
    query_end: usize,
    score_in: i32,
) -> (i32, usize, usize) {
    bl_smith_waterman_find_start_pssm(
        match_seq,
        query,
        pssm,
        query_offset,
        gap_open,
        gap_extend,
        match_seq_end,
        query_end,
        score_in,
    )
}

/// PSSM-backed find-start dispatch with forbidden ranges.
pub fn blast_smith_waterman_find_start_pssm_with_forbidden(
    match_seq: &[u8],
    query: &[u8],
    pssm: &[Vec<i32>],
    query_offset: usize,
    gap_open: i32,
    gap_extend: i32,
    match_seq_end: usize,
    query_end: usize,
    score_in: i32,
    forbidden: &BlastForbiddenRanges,
) -> (i32, usize, usize) {
    if forbidden.is_empty {
        bl_smith_waterman_find_start_pssm(
            match_seq,
            query,
            pssm,
            query_offset,
            gap_open,
            gap_extend,
            match_seq_end,
            query_end,
            score_in,
        )
    } else {
        bl_special_smith_waterman_find_start_pssm(
            match_seq,
            query,
            pssm,
            query_offset,
            gap_open,
            gap_extend,
            match_seq_end,
            query_end,
            score_in,
            forbidden,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::encoding::encode_ncbistdaa_sequence;
    use crate::matrix::AA_SIZE;

    fn blosum62() -> [[i32; AA_SIZE]; AA_SIZE] {
        *crate::api::get_matrix(crate::api::MatrixType::Blosum62)
    }

    fn offset_identity_pssm(query: &[u8], offset: usize) -> Vec<Vec<i32>> {
        let mut pssm = vec![vec![-5; AA_SIZE]; offset + query.len() + 1];
        for (i, &aa) in query.iter().enumerate() {
            pssm[offset + i][aa as usize] = 5;
        }
        pssm
    }

    #[test]
    fn blast_sw_score_only_name_matched_port_matches_public_score() {
        let q = encode_ncbistdaa_sequence(b"MKFLILLF");
        let s = encode_ncbistdaa_sequence(b"AAAAAMKFLILLF");
        let m = blosum62();

        let (expected, _, _) = blast_smith_waterman_score_only(&s, &q, &m, 11, 1, None);
        assert_eq!(
            s_smith_waterman_score_only(&q, &s, &m, None, 11, 1),
            expected
        );
        assert_eq!(
            s_smith_waterman_score_only(&s, &q, &m, None, 11, 1),
            expected
        );
    }

    #[test]
    fn blast_sw_traceback_name_matched_port_recovers_exact_hit() {
        let q = encode_ncbistdaa_sequence(b"MKFLILLF");
        let s = encode_ncbistdaa_sequence(b"AAAAAMKFLILLF");
        let m = blosum62();

        let hits = smith_waterman_score_with_traceback(&q, &s, &m, None, 11, 1, 0, 38);
        let hit = hits
            .iter()
            .find(|hit| hit.score == 38 && hit.a_start == 0 && hit.a_end == q.len())
            .expect("expected exact query hit");

        assert_eq!((hit.b_start, hit.b_end), (5, 13));
        assert_eq!(hit.edit_script.ops_vec(), vec![(GapAlignOpType::Sub, 8)]);
    }

    #[test]
    fn blast_sw_traceback_name_matched_port_swaps_square_matrix_inputs_back() {
        let q = encode_ncbistdaa_sequence(b"MKFLILLF");
        let s = encode_ncbistdaa_sequence(b"AAAAAMKFLILLF");
        let m = blosum62();

        let hits = smith_waterman_score_with_traceback(&q, &s, &m, None, 11, 1, 100, 38);
        let hit = hits
            .iter()
            .find(|hit| hit.score == 38 && hit.b_start == 105 && hit.b_end == 113)
            .expect("expected subject-offset hit after internal swap");

        assert_eq!((hit.a_start, hit.a_end), (0, 8));
        assert_eq!(hit.edit_script.ops_vec(), vec![(GapAlignOpType::Sub, 8)]);
    }

    #[test]
    fn blast_sw_traceback_name_matched_port_uses_pssm_without_swap() {
        let q = vec![1u8, 2, 3, 4];
        let s = q.clone();
        let pssm = offset_identity_pssm(&q, 0);

        let hits = smith_waterman_score_with_traceback(
            &q,
            &s,
            &crate::matrix::BLOSUM62,
            Some(&pssm),
            11,
            1,
            7,
            20,
        );
        let hit = hits
            .iter()
            .find(|hit| hit.score == 20)
            .expect("expected pssm traceback hit");

        assert_eq!((hit.a_start, hit.a_end), (0, 4));
        assert_eq!((hit.b_start, hit.b_end), (7, 11));
        assert_eq!(hit.edit_script.ops_vec(), vec![(GapAlignOpType::Sub, 4)]);
    }

    #[test]
    fn blast_smith_waterman_get_gapped_score_saves_placeholder_hsp() {
        let q = encode_ncbistdaa_sequence(b"MKFLILLF");
        let s = encode_ncbistdaa_sequence(b"AAAAAMKFLILLF");
        let m = blosum62();
        let qi = QueryInfo::new_blastp(&[q.len()]);
        let nuc_matrix = crate::traceback::blast_score_blk_nucl_matrix_create(1, -3);
        let params = BlastSmithWatermanGappedScoreParams {
            program_number: crate::program::BLASTP,
            query: &q,
            query_info: &qi,
            subject: &s,
            subject_length: s.len(),
            subject_oid: 42,
            subject_frame: 0,
            protein_matrix: &m,
            pssm: None,
            nucleotide_matrix: &nuc_matrix,
            gap_open: 11,
            gap_extend: 1,
            cutoffs: &[38],
            hsp_num_max: 10,
        };

        let (status, list) = blast_smith_waterman_get_gapped_score(&params, None);

        assert_eq!(status, 0);
        assert_eq!(list.oid, 42);
        assert_eq!(list.hsps.len(), 1);
        let hsp = &list.hsps[0];
        assert_eq!(hsp.score, 38);
        assert_eq!((hsp.query_offset, hsp.query_end), (0, q.len() as i32));
        assert_eq!((hsp.subject_offset, hsp.subject_end), (0, s.len() as i32));
        assert_eq!(hsp.context, 0);
        assert!(hsp.edit_script.is_none());
    }

    #[test]
    fn blast_smith_waterman_get_gapped_score_uses_pssm_context_rows() {
        let q = vec![1u8, 2, 3, 4, 1, 2, 3, 4];
        let s = vec![1u8, 2, 3, 4];
        let pssm = offset_identity_pssm(&q, 0);
        let qi = QueryInfo {
            num_queries: 1,
            contexts: vec![crate::queryinfo::ContextInfo {
                query_offset: 4,
                query_length: 4,
                eff_searchsp: 0,
                length_adjustment: 0,
                query_index: 0,
                frame: 0,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: 4,
            min_length: 0,
        };
        let nuc_matrix = crate::traceback::blast_score_blk_nucl_matrix_create(1, -3);
        let params = BlastSmithWatermanGappedScoreParams {
            program_number: crate::program::PSI_BLAST,
            query: &q,
            query_info: &qi,
            subject: &s,
            subject_length: s.len(),
            subject_oid: 7,
            subject_frame: 0,
            protein_matrix: &crate::matrix::BLOSUM62,
            pssm: Some(&pssm),
            nucleotide_matrix: &nuc_matrix,
            gap_open: 11,
            gap_extend: 1,
            cutoffs: &[20],
            hsp_num_max: 10,
        };

        let (_, list) = blast_smith_waterman_get_gapped_score(&params, None);

        assert_eq!(list.hsps.len(), 1);
        assert_eq!(list.hsps[0].score, 20);
        assert_eq!(list.hsps[0].context, 0);
    }

    #[test]
    fn blast_smith_waterman_get_gapped_score_scores_packed_nucleotide_subject() {
        let query = vec![0u8, 1, 2, 3, 0, 1, 2, 3];
        let subject = crate::encoding::pack_ncbi2na_bases(&query);
        let qi = QueryInfo::new_blastn(&[query.len()]);
        let prot_matrix = blosum62();
        let nuc_matrix = crate::traceback::blast_score_blk_nucl_matrix_create(1, -3);
        let params = BlastSmithWatermanGappedScoreParams {
            program_number: crate::program::BLASTN,
            query: &query,
            query_info: &qi,
            subject: &subject,
            subject_length: query.len(),
            subject_oid: 3,
            subject_frame: 1,
            protein_matrix: &prot_matrix,
            pssm: None,
            nucleotide_matrix: &nuc_matrix,
            gap_open: 5,
            gap_extend: 2,
            cutoffs: &[8, 100],
            hsp_num_max: 10,
        };

        let (_, list) = blast_smith_waterman_get_gapped_score(&params, None);

        assert_eq!(list.hsps.len(), 1);
        assert_eq!(list.hsps[0].score, 8);
        assert_eq!(list.hsps[0].subject_end, query.len() as i32);
    }

    #[test]
    fn blast_sw_score_only_name_matched_port_uses_pssm_rows_without_swap() {
        let q = vec![1u8, 2, 3, 4];
        let s = q.clone();
        let pssm = offset_identity_pssm(&q, 0);

        assert_eq!(
            s_smith_waterman_score_only(&q, &s, &crate::matrix::BLOSUM62, Some(&pssm), 11, 1),
            20
        );
    }

    #[test]
    fn nucl_smith_waterman_name_matched_port_scores_packed_subject() {
        let query = vec![0u8, 1, 2, 3, 0, 1, 2, 3];
        let subject = crate::encoding::pack_ncbi2na_bases(&query);
        let matrix = crate::traceback::blast_score_blk_nucl_matrix_create(1, -3);

        assert_eq!(
            s_nucl_smith_waterman(&subject, query.len(), &query, &matrix, 5, 2),
            query.len() as i32
        );
    }

    #[test]
    fn test_score_only_exact_match() {
        let q = encode_ncbistdaa_sequence(b"MKFLILLF");
        let s = encode_ncbistdaa_sequence(b"MKFLILLF");
        let m = blosum62();
        let (score, m_end, q_end) = blast_smith_waterman_score_only(&s, &q, &m, 11, 1, None);
        // Diag scores: 5+5+6+4+4+4+4+6 = 38
        assert_eq!(score, 38);
        assert_eq!(m_end, 7); // last position of subject
        assert_eq!(q_end, 7); // last position of query
    }

    #[test]
    fn test_find_start_exact_match() {
        let q = encode_ncbistdaa_sequence(b"MKFLILLF");
        let s = encode_ncbistdaa_sequence(b"MKFLILLF");
        let m = blosum62();
        let (score, m_start, q_start) =
            blast_smith_waterman_find_start(&s, &q, &m, 11, 1, 7, 7, 38, None);
        assert_eq!(score, 38);
        assert_eq!(m_start, 0);
        assert_eq!(q_start, 0);
    }

    #[test]
    fn test_score_only_offset_match() {
        // Subject has 5 leading X's then the query
        let q = encode_ncbistdaa_sequence(b"MKFLILLF");
        let s = encode_ncbistdaa_sequence(b"AAAAAMKFLILLF");
        let m = blosum62();
        let (score, m_end, q_end) = blast_smith_waterman_score_only(&s, &q, &m, 11, 1, None);
        assert_eq!(score, 38);
        assert_eq!(q_end, 7);
        assert_eq!(m_end, 12); // 5 + 7 = position of last F in subject

        // Reverse SW should find start at position 5 in subject
        let (score2, m_start, q_start) =
            blast_smith_waterman_find_start(&s, &q, &m, 11, 1, m_end, q_end, score, None);
        assert_eq!(score2, 38);
        assert_eq!(q_start, 0);
        assert_eq!(m_start, 5);
    }

    #[test]
    fn test_score_only_no_match() {
        let q = encode_ncbistdaa_sequence(b"AAAA");
        let s = encode_ncbistdaa_sequence(b"WWWW");
        let m = blosum62();
        let (score, _, _) = blast_smith_waterman_score_only(&s, &q, &m, 11, 1, None);
        // BLOSUM62 A-W = -3, so all-negative; SW returns 0
        assert_eq!(score, 0);
    }

    #[test]
    fn test_score_only_with_gap() {
        // Aligning MKF + gap + LILLF vs MKFXLILLF (X is a substitution penalty)
        let q = encode_ncbistdaa_sequence(b"MKFLILLF");
        let s = encode_ncbistdaa_sequence(b"MKFGGLILLF");
        let m = blosum62();
        let (score, _, _) = blast_smith_waterman_score_only(&s, &q, &m, 11, 1, None);
        // SW will find best alignment which includes a gap or extension.
        // Just verify it returns some positive score >= 0.
        assert!(score > 0);
    }

    #[test]
    fn test_forbidden_ranges_basic() {
        // Subject has the query motif twice; without forbidden ranges we
        // find one of them (whichever the SW picks). Pushing the first
        // match's range as forbidden should make the second SW pass find
        // the *other* occurrence.
        let q = encode_ncbistdaa_sequence(b"MKFLILLF");
        let s = encode_ncbistdaa_sequence(b"MKFLILLFAAAAAMKFLILLF");
        let m = blosum62();
        // First pass — empty forbidden set, picks the first occurrence.
        let mut fr = BlastForbiddenRanges::new(q.len() as i32);
        assert!(fr.is_empty);
        let (s1, m_end1, q_end1) =
            blast_smith_waterman_score_only_with_forbidden(&s, &q, &m, 11, 1, &fr);
        assert_eq!(s1, 38);
        let (s1_back, m_start1, q_start1) = blast_smith_waterman_find_start_with_forbidden(
            &s, &q, &m, 11, 1, m_end1, q_end1, s1, &fr,
        );
        assert_eq!(s1_back, 38);
        assert_eq!(q_start1, 0);

        // Mark the first occurrence as forbidden, do the second pass.
        fr.push(
            q_start1 as i32,
            q_end1 as i32 + 1,
            m_start1 as i32,
            m_end1 as i32,
        );
        assert!(!fr.is_empty);
        let (s2, m_end2, _q_end2) =
            blast_smith_waterman_score_only_with_forbidden(&s, &q, &m, 11, 1, &fr);
        assert_eq!(s2, 38, "second occurrence should give same score");
        // Second-pass end position must be inside the second occurrence
        // (positions 13..=20). The first occurrence ends at index 7.
        assert!(
            m_end2 >= 13,
            "expected second-pass end inside second occurrence, got {}",
            m_end2
        );
    }

    #[test]
    fn test_forbidden_ranges_clear_resets_is_empty() {
        let mut fr = BlastForbiddenRanges::new(10);
        fr.push(0, 5, 0, 7);
        assert!(!fr.is_empty);
        fr.clear();
        assert!(fr.is_empty);
        assert!(fr.num_forbidden.iter().all(|&n| n == 0));
    }

    #[test]
    fn test_smith_waterman_pssm_respects_query_offset() {
        let q = vec![1u8, 2, 3, 4];
        let s = q.clone();
        let pssm = offset_identity_pssm(&q, 3);

        let (score, m_end, q_end) = blast_smith_waterman_score_only_pssm(&s, &q, &pssm, 3, 11, 1);
        assert_eq!(score, 20);
        assert_eq!(m_end, 3);
        assert_eq!(q_end, 3);

        let (score2, m_start, q_start) =
            blast_smith_waterman_find_start_pssm(&s, &q, &pssm, 3, 11, 1, m_end, q_end, score);
        assert_eq!(score2, 20);
        assert_eq!(m_start, 0);
        assert_eq!(q_start, 0);
    }

    #[test]
    fn test_smith_waterman_pssm_forbidden_ranges() {
        let q = vec![1u8, 2, 3, 4];
        let s = vec![1u8, 2, 3, 4, 1, 2, 3, 4];
        let pssm = offset_identity_pssm(&q, 2);
        let mut fr = BlastForbiddenRanges::new(q.len() as i32);

        let (score1, m_end1, q_end1) =
            blast_smith_waterman_score_only_pssm_with_forbidden(&s, &q, &pssm, 2, 11, 1, &fr);
        assert_eq!(score1, 20);
        let (score1_back, m_start1, q_start1) = blast_smith_waterman_find_start_pssm_with_forbidden(
            &s, &q, &pssm, 2, 11, 1, m_end1, q_end1, score1, &fr,
        );
        assert_eq!(score1_back, 20);
        assert_eq!(q_start1, 0);

        fr.push(
            q_start1 as i32,
            q_end1 as i32 + 1,
            m_start1 as i32,
            m_end1 as i32,
        );
        let (score2, m_end2, _q_end2) =
            blast_smith_waterman_score_only_pssm_with_forbidden(&s, &q, &pssm, 2, 11, 1, &fr);
        assert_eq!(score2, 20);
        assert!(m_end2 >= 4, "expected second occurrence, got {m_end2}");

        let mut all_forbidden = BlastForbiddenRanges::new(q.len() as i32);
        all_forbidden.push(0, q.len() as i32, 0, (s.len() - 1) as i32);
        let (score3, _, _) = blast_smith_waterman_score_only_pssm_with_forbidden(
            &s,
            &q,
            &pssm,
            2,
            11,
            1,
            &all_forbidden,
        );
        assert_eq!(score3, 0);
    }
}
