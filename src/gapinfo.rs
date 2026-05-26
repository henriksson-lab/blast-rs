//! Rust equivalent of gapinfo.c — gap edit script structures.

/// Type of gap alignment operation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u32)]
pub enum GapAlignOpType {
    Del = 0, // Deletion in query (gap in query)
    Del2 = 1,
    Del1 = 2,
    Sub = 3,  // Substitution (aligned pair)
    Ins1 = 4, // Insertion in query (gap in subject)
    Ins2 = 5,
    Ins = 6,
    Decline = 7,
}

/// A gap edit script representing a gapped alignment.
#[derive(Debug, Clone, Default)]
pub struct GapEditScript {
    pub ops: Vec<(GapAlignOpType, i32)>, // (operation, count) pairs
}

#[derive(Debug, Clone, Default)]
pub struct GapStateArrayStruct {
    pub state_array: Vec<i32>,
    pub length: i32,
    pub used: bool,
    pub next: Option<Box<GapStateArrayStruct>>,
}

const GAP_STATE_CHUNKSIZE: usize = 2_097_152;
const GAP_STATE_MIN_CELLS: usize = GAP_STATE_CHUNKSIZE / std::mem::size_of::<i32>();

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GapPrelimEditScript {
    pub op_type: GapAlignOpType,
    pub num: i32,
}

#[derive(Debug, Clone)]
pub struct GapPrelimEditBlock {
    pub edit_ops: Vec<GapPrelimEditScript>,
    pub last_op: Option<GapAlignOpType>,
}

pub type JumperOpType = i32;
pub const JUMPER_MISMATCH: JumperOpType = 0;
pub const JUMPER_INSERTION: JumperOpType = -1;
pub const JUMPER_DELETION: JumperOpType = -2;
pub const MAPPER_EXON: u8 = 0x40;
pub const MAPPER_SPLICE_SIGNAL: u8 = 0x80;
pub const MAPPER_POLY_A: u8 = 0x20;
pub const MAPPER_ADAPTER: u8 = 0x10;

/// NCBI: s_GapCost (phi_gapalign.c:281).
pub fn s_gap_cost(gap_open: i32, gap_extend: i32, length: i32) -> i32 {
    if length <= 0 {
        0
    } else {
        gap_open + gap_extend * length
    }
}

#[derive(Debug, Clone, Copy, Default)]
struct PhiGapDp {
    best: i32,
    best_gap: i32,
}

const PHI_DIAGONAL_INSERT: i32 = 1;
const PHI_DIAGONAL_DELETE: i32 = 2;
const PHI_INSERT_CODE: i32 = 10;
const PHI_DELETE_CODE: i32 = 20;

/// blast-rs: Safe matrix lookup helper for PHI DP; not a direct NCBI C port.
fn phi_matrix_score(matrix: &[Vec<i32>], a: u8, b: u8) -> i32 {
    matrix
        .get(a as usize)
        .and_then(|row| row.get(b as usize))
        .copied()
        .unwrap_or(0)
}

/// NCBI: s_Align (phi_gapalign.c:93).
fn s_align(
    seq1: &[u8],
    seq2: &[u8],
    end1: i32,
    end2: i32,
    low_diag: i32,
    high_diag: i32,
    matrix: &[Vec<i32>],
    gap_open: i32,
    gap_extend: i32,
    align_script: &mut GapPrelimEditBlock,
) -> i32 {
    let band = high_diag - low_diag + 1;
    if end1 < 0 || end2 < 0 || band <= 0 {
        return 0;
    }

    let band_cells = (band + 2).max(0) as usize;
    let rows = (end1 + 1).max(0) as usize;
    let mut score_array = vec![PhiGapDp::default(); band_cells];
    let mut state = vec![vec![0i32; band_cells]; rows];
    let k_min_score = i32::MIN / 2;
    let gap_cost = gap_open + gap_extend;

    let mut leftd = 1 - low_diag;
    let mut rightd = high_diag - low_diag + 1;
    score_array[leftd as usize].best = 0;
    state[0][leftd as usize] = -1;
    let mut initial_score = -gap_open;
    for diag_index in leftd + 1..=rightd {
        initial_score -= gap_extend;
        score_array[diag_index as usize].best = initial_score;
        score_array[(diag_index - 1) as usize].best_gap = initial_score - gap_cost;
        state[0][diag_index as usize] = PHI_DIAGONAL_INSERT;
    }
    score_array[(rightd + 1) as usize].best = k_min_score;
    score_array[rightd as usize].best_gap = k_min_score;
    score_array[(leftd - 1) as usize].best_gap = -gap_cost;
    score_array[(leftd - 1) as usize].best = k_min_score;

    for i in 1..=end1 {
        if i > end2 - high_diag {
            rightd -= 1;
        }
        if leftd > 1 {
            leftd -= 1;
        }

        let mut temp_indel_score = score_array[leftd as usize].best_gap;
        let mut next_state = 0;
        let index2 = leftd + low_diag - 1 + i;
        let mut temp_sub_score = if index2 > 0 {
            score_array[leftd as usize].best
                + phi_matrix_score(matrix, seq1[i as usize], seq2[index2 as usize])
        } else {
            0
        };
        if temp_indel_score > temp_sub_score || index2 <= 0 {
            temp_sub_score = temp_indel_score;
            next_state = PHI_DIAGONAL_DELETE;
        }
        let mut temp_hor_score = temp_sub_score - gap_cost;
        if leftd >= 1 {
            temp_indel_score -= gap_extend;
            if temp_indel_score >= temp_hor_score {
                score_array[(leftd - 1) as usize].best_gap = temp_indel_score;
                next_state += PHI_DELETE_CODE;
            } else {
                score_array[(leftd - 1) as usize].best_gap = temp_hor_score;
            }
        }
        state[i as usize][leftd as usize] = next_state;
        score_array[leftd as usize].best = temp_sub_score;

        for curd in leftd + 1..=rightd {
            let curd_usize = curd as usize;
            temp_sub_score = score_array[curd_usize].best
                + phi_matrix_score(
                    matrix,
                    seq1[i as usize],
                    seq2[(curd + low_diag - 1 + i) as usize],
                );
            temp_indel_score = score_array[curd_usize].best_gap;
            if temp_indel_score > temp_sub_score {
                if temp_indel_score > temp_hor_score {
                    score_array[curd_usize].best = temp_indel_score;
                    temp_hor_score -= gap_extend;
                    score_array[(curd - 1) as usize].best_gap = temp_indel_score - gap_extend;
                    state[i as usize][curd_usize] =
                        PHI_DELETE_CODE + PHI_INSERT_CODE + PHI_DIAGONAL_DELETE;
                } else {
                    score_array[curd_usize].best = temp_hor_score;
                    temp_hor_score -= gap_extend;
                    score_array[(curd - 1) as usize].best_gap = temp_indel_score - gap_extend;
                    state[i as usize][curd_usize] =
                        PHI_DELETE_CODE + PHI_INSERT_CODE + PHI_DIAGONAL_INSERT;
                }
            } else if temp_hor_score > temp_sub_score {
                score_array[curd_usize].best = temp_hor_score;
                temp_hor_score -= gap_extend;
                score_array[(curd - 1) as usize].best_gap = temp_indel_score - gap_extend;
                state[i as usize][curd_usize] =
                    PHI_DELETE_CODE + PHI_INSERT_CODE + PHI_DIAGONAL_INSERT;
            } else {
                score_array[curd_usize].best = temp_sub_score;
                temp_sub_score -= gap_cost;
                temp_hor_score -= gap_extend;
                if temp_sub_score > temp_hor_score {
                    temp_hor_score = temp_sub_score;
                    next_state = 0;
                } else {
                    next_state = PHI_INSERT_CODE;
                }
                temp_indel_score -= gap_extend;
                if temp_sub_score > temp_indel_score {
                    state[i as usize][curd_usize] = next_state;
                    score_array[(curd - 1) as usize].best_gap = temp_sub_score;
                } else {
                    state[i as usize][curd_usize] = next_state + PHI_DELETE_CODE;
                    score_array[(curd - 1) as usize].best_gap = temp_indel_score;
                }
            }
        }
    }

    let score = score_array[rightd as usize].best;
    let mut edit_instructions = Vec::with_capacity((end1 + end2).max(0) as usize);
    let mut index1 = end1;
    let mut diag_index = rightd;
    let mut edit_step = 0;
    while index1 >= 0 {
        let state_decoder = state[index1 as usize][diag_index as usize];
        let mut next_edit_step = state_decoder % PHI_INSERT_CODE;
        if state_decoder == -1 {
            break;
        }
        if edit_step == PHI_DIAGONAL_INSERT && ((state_decoder / PHI_INSERT_CODE) % 2) == 1 {
            next_edit_step = PHI_DIAGONAL_INSERT;
        }
        if edit_step == PHI_DIAGONAL_DELETE && (state_decoder / PHI_DELETE_CODE) == 1 {
            next_edit_step = PHI_DIAGONAL_DELETE;
        }
        if next_edit_step == PHI_DIAGONAL_INSERT {
            diag_index -= 1;
            index1 += 1;
        } else if next_edit_step == PHI_DIAGONAL_DELETE {
            diag_index += 1;
        }
        edit_step = next_edit_step;
        edit_instructions.push(edit_step);
        index1 -= 1;
    }

    for &instruction in edit_instructions.iter().rev() {
        match instruction {
            0 => gap_prelim_edit_block_add(align_script, GapAlignOpType::Sub, 1),
            PHI_DIAGONAL_INSERT => gap_prelim_edit_block_add(align_script, GapAlignOpType::Del, 1),
            PHI_DIAGONAL_DELETE => gap_prelim_edit_block_add(align_script, GapAlignOpType::Ins, 1),
            _ => {}
        }
    }

    score
}

/// NCBI: s_BandedAlign (phi_gapalign.c:280).
pub fn s_banded_align(
    seq1: &[u8],
    seq2: &[u8],
    start1: i32,
    start2: i32,
    mut low_diag: i32,
    mut high_diag: i32,
    matrix: &[Vec<i32>],
    gap_open: i32,
    gap_extend: i32,
    align_script: &mut GapPrelimEditBlock,
) -> i32 {
    low_diag = (-start1).max(low_diag).min((start2 - start1).min(0));
    high_diag = start2.min(high_diag).max((start2 - start1).max(0));
    if start2 <= 0 {
        if start1 > 0 {
            gap_prelim_edit_block_add(align_script, GapAlignOpType::Ins, start1);
        }
        return -s_gap_cost(gap_open, gap_extend, start1);
    }
    if start1 <= 0 {
        gap_prelim_edit_block_add(align_script, GapAlignOpType::Del, start2);
        return -s_gap_cost(gap_open, gap_extend, start2);
    }
    if high_diag - low_diag < 1 {
        let mut score = 0;
        for i in 1..=start1 {
            gap_prelim_edit_block_add(align_script, GapAlignOpType::Sub, 1);
            score += phi_matrix_score(matrix, seq1[i as usize], seq2[i as usize]);
        }
        return score;
    }

    s_align(
        seq1,
        seq2,
        start1,
        start2,
        low_diag,
        high_diag,
        matrix,
        gap_open,
        gap_extend,
        align_script,
    )
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Jump {
    pub dcp: i32,
    pub dcq: i32,
    pub lng: i32,
    pub ok: i32,
}

/// NCBI `jumper_default` table (`jumper.c:41`).
pub const JUMPER_DEFAULT: [Jump; 15] = [
    Jump {
        dcp: 1,
        dcq: 1,
        lng: 9,
        ok: 0,
    },
    Jump {
        dcp: 1,
        dcq: 0,
        lng: 10,
        ok: 0,
    },
    Jump {
        dcp: 0,
        dcq: 1,
        lng: 10,
        ok: 0,
    },
    Jump {
        dcp: 2,
        dcq: 0,
        lng: 10,
        ok: 0,
    },
    Jump {
        dcp: 0,
        dcq: 2,
        lng: 10,
        ok: 0,
    },
    Jump {
        dcp: 3,
        dcq: 0,
        lng: 13,
        ok: 0,
    },
    Jump {
        dcp: 0,
        dcq: 3,
        lng: 13,
        ok: 0,
    },
    Jump {
        dcp: 2,
        dcq: 2,
        lng: 12,
        ok: 0,
    },
    Jump {
        dcp: 1,
        dcq: 0,
        lng: 10,
        ok: 2,
    },
    Jump {
        dcp: 0,
        dcq: 1,
        lng: 10,
        ok: 2,
    },
    Jump {
        dcp: 2,
        dcq: 0,
        lng: 10,
        ok: 2,
    },
    Jump {
        dcp: 0,
        dcq: 2,
        lng: 10,
        ok: 2,
    },
    Jump {
        dcp: 3,
        dcq: 0,
        lng: 13,
        ok: 2,
    },
    Jump {
        dcp: 0,
        dcq: 3,
        lng: 13,
        ok: 2,
    },
    Jump {
        dcp: 1,
        dcq: 1,
        lng: 0,
        ok: 0,
    },
];

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct JumperPrelimEditBlock {
    pub edit_ops: Vec<JumperOpType>,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct JumperGapAlign {
    pub left_prelim_block: Option<JumperPrelimEditBlock>,
    pub right_prelim_block: Option<JumperPrelimEditBlock>,
    pub table: Vec<u32>,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct JumperEdit {
    pub query_pos: i32,
    pub query_base: u8,
    pub subject_base: u8,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct JumperEditsBlock {
    pub edits: Vec<JumperEdit>,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct SequenceOverhangs {
    pub left: Option<Vec<u8>>,
    pub right: Option<Vec<u8>>,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct SubjectIndex {
    pub num_lookups: i32,
    pub width: i32,
    pub word_size: i32,
    pub positions_by_word: std::collections::BTreeMap<u32, Vec<i32>>,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct SubjectIndexIterator {
    pub subject_index: Option<SubjectIndex>,
    pub to: i32,
    pub lookup_index: i32,
    pub positions: Vec<i32>,
    pub word_index: i32,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct JumperAlignParams {
    pub max_mismatches: i32,
    pub mismatch_window: i32,
    pub gap_x_dropoff: i32,
}

/// blast-rs: Rust ownership equivalent of NCBI `GapStateFree` (gapinfo.c:38); not a direct NCBI C port.
pub fn gap_state_free(
    state_struct: Option<Box<GapStateArrayStruct>>,
) -> Option<Box<GapStateArrayStruct>> {
    let mut state_struct = state_struct;
    while let Some(mut current) = state_struct {
        current.state_array.clear();
        state_struct = current.next.take();
    }
    None
}

/// NCBI: s_GapGetState (blast_gapalign.c:71).
///
/// The C helper reuses the first unused state-array node whose allocation is
/// long enough; otherwise it appends a new node whose backing storage is at
/// least `CHUNKSIZE` bytes. The returned node is marked `used`.
pub fn s_gap_get_state(
    head: &mut Option<Box<GapStateArrayStruct>>,
    length: i32,
) -> Option<&mut GapStateArrayStruct> {
    let requested_len = length.max(0) as usize;
    if head.is_none() {
        *head = Some(Box::new(gap_state_array_new(requested_len)));
    }
    let mut cursor = head.as_deref_mut();
    while let Some(node) = cursor {
        if !node.used && node.length as usize >= requested_len {
            node.used = true;
            return Some(node);
        }
        if node.next.is_none() {
            node.next = Some(Box::new(gap_state_array_new(requested_len)));
        }
        cursor = node.next.as_deref_mut();
    }
    None
}

/// NCBI: s_GapPurgeState (blast_gapalign.c:128).
///
/// The C routine clears `used` on the supplied state node and following nodes,
/// making the arrays available to later `s_GapGetState` calls.
pub fn s_gap_purge_state(state_struct: &mut GapStateArrayStruct) -> bool {
    let mut cursor = Some(state_struct);
    while let Some(node) = cursor {
        node.used = false;
        cursor = node.next.as_deref_mut();
    }
    true
}

/// blast-rs: Allocates one Rust state-array pool node; not a direct NCBI C port.
fn gap_state_array_new(requested_len: usize) -> GapStateArrayStruct {
    let alloc_len = requested_len.max(GAP_STATE_MIN_CELLS);
    GapStateArrayStruct {
        state_array: vec![0; alloc_len],
        length: alloc_len as i32,
        used: false,
        next: None,
    }
}

/// NCBI: GapEditScriptNew (gapinfo.c:55).
pub fn gap_edit_script_new(size: i32) -> Option<GapEditScript> {
    if size <= 0 {
        return None;
    }
    Some(GapEditScript::with_capacity(size as usize))
}

/// blast-rs: Rust ownership equivalent of NCBI `GapEditScriptDelete` (gapinfo.c:74); not a direct NCBI C port.
pub fn gap_edit_script_delete(_: Option<GapEditScript>) -> Option<GapEditScript> {
    None
}

/// NCBI: GapEditScriptDup (gapinfo.c:88).
pub fn gap_edit_script_dup(old: Option<&GapEditScript>) -> Option<GapEditScript> {
    old.cloned()
}

/// NCBI: GapEditScriptPartialCopy (gapinfo.c:104).
pub fn gap_edit_script_partial_copy(
    new_script: &mut GapEditScript,
    offset: i32,
    old_script: Option<&GapEditScript>,
    start: i32,
    stop: i32,
) -> i16 {
    let Some(old_script) = old_script else {
        return -1;
    };
    if offset < 0 || start < 0 || stop < start {
        return -1;
    }

    let offset = offset as usize;
    let start = start as usize;
    let stop = stop as usize;
    let size = stop - start + 1;
    let needed = offset + size;
    if stop >= old_script.ops.len() || needed > new_script.ops.capacity() {
        return -1;
    }

    if new_script.ops.len() < needed {
        new_script.ops.resize(needed, (GapAlignOpType::Sub, 0));
    }
    let mut old_index = start;
    for new_index in offset..needed {
        new_script.ops[new_index] = old_script.ops[old_index];
        old_index += 1;
    }
    0
}

/// NCBI: s_GapPrelimEditBlockRealloc (gapinfo.c:132).
pub fn s_gap_prelim_edit_block_realloc(edit_block: &mut GapPrelimEditBlock, total_ops: i32) -> i16 {
    if total_ops < 0 {
        return -1;
    }
    let total_ops = total_ops as usize;
    if edit_block.edit_ops.capacity() <= total_ops {
        let new_size = total_ops.saturating_mul(2).max(1);
        let additional = new_size.saturating_sub(edit_block.edit_ops.capacity());
        if edit_block.edit_ops.try_reserve_exact(additional).is_err() {
            return -1;
        }
    }
    0
}

/// NCBI: s_GapPrelimEditBlockAddNew (gapinfo.c:157).
pub fn s_gap_prelim_edit_block_add_new(
    edit_block: &mut GapPrelimEditBlock,
    op_type: GapAlignOpType,
    num_ops: i32,
) -> i16 {
    if op_type == GapAlignOpType::Decline {
        return -1;
    }
    if s_gap_prelim_edit_block_realloc(edit_block, edit_block.edit_ops.len() as i32 + 2) != 0 {
        return -1;
    }
    edit_block.last_op = Some(op_type);
    edit_block.edit_ops.push(GapPrelimEditScript {
        op_type,
        num: num_ops,
    });
    0
}

/// NCBI: GapPrelimEditBlockAdd (gapinfo.c:174).
pub fn gap_prelim_edit_block_add(
    edit_block: &mut GapPrelimEditBlock,
    op_type: GapAlignOpType,
    num_ops: i32,
) {
    if num_ops == 0 {
        return;
    }
    if edit_block.last_op == Some(op_type) {
        if let Some(last) = edit_block.edit_ops.last_mut() {
            last.num += num_ops;
            return;
        }
    }
    let _ = s_gap_prelim_edit_block_add_new(edit_block, op_type, num_ops);
}

/// NCBI: GapPrelimEditBlockNew (gapinfo.c:187).
pub fn gap_prelim_edit_block_new() -> GapPrelimEditBlock {
    let mut edit_block = GapPrelimEditBlock {
        edit_ops: Vec::new(),
        last_op: None,
    };
    let _ = s_gap_prelim_edit_block_realloc(&mut edit_block, 100);
    edit_block
}

/// blast-rs: Rust ownership equivalent of NCBI `GapPrelimEditBlockFree`; not a direct NCBI C port.
pub fn gap_prelim_edit_block_free(_: Option<GapPrelimEditBlock>) -> Option<GapPrelimEditBlock> {
    None
}

/// NCBI: GapPrelimEditBlockReset (gapinfo.c:212).
pub fn gap_prelim_edit_block_reset(edit_block: Option<&mut GapPrelimEditBlock>) {
    if let Some(edit_block) = edit_block {
        edit_block.edit_ops.clear();
        edit_block.last_op = None;
    }
}

/// NCBI: GapPrelimEditBlockAppend (gapinfo.c:221).
pub fn gap_prelim_edit_block_append(
    edit_block1: &mut GapPrelimEditBlock,
    edit_block2: &GapPrelimEditBlock,
) {
    for op in &edit_block2.edit_ops {
        gap_prelim_edit_block_add(edit_block1, op.op_type, op.num);
    }
}

/// NCBI: Blast_PrelimEditBlockToGapEditScript (blast_gapalign.c:2481).
pub fn blast_prelim_edit_block_to_gap_edit_script(
    rev_prelim_tback: Option<&GapPrelimEditBlock>,
    fwd_prelim_tback: Option<&GapPrelimEditBlock>,
) -> Option<GapEditScript> {
    let rev_prelim_tback = rev_prelim_tback?;
    let fwd_prelim_tback = fwd_prelim_tback?;

    let merge_ops = fwd_prelim_tback
        .edit_ops
        .last()
        .zip(rev_prelim_tback.edit_ops.last())
        .is_some_and(|(fwd, rev)| fwd.op_type == rev.op_type);
    let mut size = fwd_prelim_tback.edit_ops.len() + rev_prelim_tback.edit_ops.len();
    if merge_ops {
        size = size.saturating_sub(1);
    }
    if size == 0 {
        return None;
    }

    let mut script = GapEditScript::with_capacity(size);
    for op in &rev_prelim_tback.edit_ops {
        script.ops.push((op.op_type, op.num));
    }

    if fwd_prelim_tback.edit_ops.is_empty() {
        return Some(script);
    }

    let mut fwd_ops = fwd_prelim_tback.edit_ops.as_slice();
    if merge_ops {
        if let (Some(last_script), Some(last_fwd)) = (script.ops.last_mut(), fwd_ops.last()) {
            last_script.1 += last_fwd.num;
        }
        fwd_ops = &fwd_ops[..fwd_ops.len() - 1];
    }

    for op in fwd_ops.iter().rev() {
        script.ops.push((op.op_type, op.num));
    }

    Some(script)
}

/// blast-rs: OOF nucleotide-span helper extracted from traceback conversion; not a direct NCBI C port.
fn oof_op_nucleotide_span(op: GapAlignOpType) -> i32 {
    match op {
        GapAlignOpType::Ins => GapAlignOpType::Sub as i32,
        other => other as i32,
    }
}

/// NCBI: s_BlastOOFTracebackToGapEditScript (blast_gapalign.c:4451).
///
/// OOF traceback stores frame-shift operations in the same numeric space as
/// `EGapAlignOpType` (`Del2`/`Del1`/`Ins1`/`Ins2`). This helper performs the C
/// post-processing: prepend the shifted substitution, flow frame-shift ops
/// through in-frame gaps, merge the left/right traceback halves, truncate by
/// translated nucleotide alignment length, split multi-count frame shifts into
/// single-op elements, and lengthen substitutions that follow a frame shift.
pub fn s_blast_oof_traceback_to_gap_edit_script(
    rev_prelim_tback: Option<&GapPrelimEditBlock>,
    fwd_prelim_tback: Option<&GapPrelimEditBlock>,
    nucl_align_length: i32,
    edit_script_ptr: Option<&mut Option<GapEditScript>>,
) -> i16 {
    let (Some(rev_prelim_tback), Some(fwd_prelim_tback), Some(edit_script_ptr)) =
        (rev_prelim_tback, fwd_prelim_tback, edit_script_ptr)
    else {
        return -1;
    };

    let mut tmp_prelim_tback = gap_prelim_edit_block_new();
    let mut last_op = GapAlignOpType::Sub;
    let mut last_num = 1;

    for next in &rev_prelim_tback.edit_ops {
        if next.op_type == last_op {
            last_num += next.num;
        } else if matches!(next.op_type, GapAlignOpType::Ins | GapAlignOpType::Del) {
            if last_num > 1 {
                gap_prelim_edit_block_add(&mut tmp_prelim_tback, last_op, last_num - 1);
            }
            gap_prelim_edit_block_add(&mut tmp_prelim_tback, next.op_type, next.num);
            last_num = 1;
        } else {
            gap_prelim_edit_block_add(&mut tmp_prelim_tback, last_op, last_num);
            last_op = next.op_type;
            last_num = next.num;
        }
    }

    last_num -= 1;
    if last_num > 0 {
        gap_prelim_edit_block_add(&mut tmp_prelim_tback, last_op, last_num);
    }

    let mut fwd_ops = fwd_prelim_tback.edit_ops.clone();
    if last_op != GapAlignOpType::Sub {
        let mut keep_len = fwd_ops.len();
        let mut merged_frame_shift = false;
        for idx in (0..fwd_ops.len()).rev() {
            let new_script = &mut fwd_ops[idx];
            if matches!(
                new_script.op_type,
                GapAlignOpType::Ins | GapAlignOpType::Del
            ) {
                gap_prelim_edit_block_add(
                    &mut tmp_prelim_tback,
                    new_script.op_type,
                    new_script.num,
                );
            } else {
                let merged_value =
                    last_op as u32 + new_script.op_type as u32 - GapAlignOpType::Sub as u32;
                let merged_op = match merged_value {
                    0 => GapAlignOpType::Del,
                    1 => GapAlignOpType::Del2,
                    2 => GapAlignOpType::Del1,
                    3 => GapAlignOpType::Sub,
                    4 => GapAlignOpType::Ins1,
                    5 => GapAlignOpType::Ins2,
                    6 => GapAlignOpType::Ins,
                    _ => GapAlignOpType::Sub,
                };
                gap_prelim_edit_block_add(&mut tmp_prelim_tback, merged_op, 1);
                new_script.num -= 1;
                keep_len = if new_script.num == 0 { idx } else { idx + 1 };
                merged_frame_shift = true;
                break;
            }
        }
        if !merged_frame_shift {
            keep_len = 0;
        }
        fwd_ops.truncate(keep_len);
    }

    let fwd_prelim = GapPrelimEditBlock {
        last_op: fwd_ops.last().map(|op| op.op_type),
        edit_ops: fwd_ops,
    };
    let Some(mut script) =
        blast_prelim_edit_block_to_gap_edit_script(Some(&tmp_prelim_tback), Some(&fwd_prelim))
    else {
        *edit_script_ptr = None;
        return 0;
    };
    if nucl_align_length <= 0 {
        *edit_script_ptr = None;
        return 0;
    }

    let mut num_nuc = 0;
    let mut keep_size = script.ops.len();
    for (idx, (op, num)) in script.ops.iter_mut().enumerate() {
        let span = oof_op_nucleotide_span(*op);
        let total_actions = span * *num;
        if num_nuc + total_actions >= nucl_align_length {
            *num = (nucl_align_length - num_nuc + span - 1) / span;
            keep_size = idx + 1;
            break;
        }
        num_nuc += total_actions;
    }
    script.ops.truncate(keep_size);

    let mut split_script = GapEditScript::new();
    for (op, num) in script.ops {
        if (op as u32) % 3 != 0 && num > 1 {
            for _ in 0..num {
                split_script.ops.push((op, 1));
            }
        } else {
            split_script.ops.push((op, num));
        }
    }

    if let Some((first, rest)) = split_script.ops.split_first_mut() {
        let mut last_op = first.0;
        for (op, num) in rest {
            if *op == GapAlignOpType::Sub && (last_op as u32) % 3 != 0 {
                *num += 1;
            }
            last_op = *op;
        }
    }

    *edit_script_ptr = Some(split_script);
    0
}

/// NCBI: JumperPrelimEditBlockNew (jumper.c:108).
pub fn jumper_prelim_edit_block_new(size: i32) -> Option<JumperPrelimEditBlock> {
    if size <= 0 {
        return None;
    }
    Some(JumperPrelimEditBlock {
        edit_ops: Vec::with_capacity(size as usize),
    })
}

/// blast-rs: Rust ownership equivalent of NCBI `JumperPrelimEditBlockFree`; not a direct NCBI C port.
pub fn jumper_prelim_edit_block_free(
    _: Option<JumperPrelimEditBlock>,
) -> Option<JumperPrelimEditBlock> {
    None
}

/// NCBI: JumperPrelimEditBlockAdd (jumper.c:139).
pub fn jumper_prelim_edit_block_add(block: &mut JumperPrelimEditBlock, op: JumperOpType) -> i32 {
    if block.edit_ops.len() >= block.edit_ops.capacity() {
        let grow_by = block.edit_ops.capacity().max(1);
        if block.edit_ops.try_reserve_exact(grow_by).is_err() {
            return -1;
        }
    }
    if op > 0 {
        if let Some(last) = block.edit_ops.last_mut() {
            if *last > 0 {
                *last += op;
                return 0;
            }
        }
    }
    block.edit_ops.push(op);
    0
}

/// NCBI: s_CreateTable (jumper.c:164).
pub fn s_create_table(table: &mut [u32]) {
    for (byte, entry) in table.iter_mut().enumerate() {
        let mut packed = 0u32;
        for pos in 0..4 {
            let base = ((byte >> (2 * (3 - pos))) & 3) as u32;
            packed |= base << (pos * 8);
        }
        *entry = packed;
    }
}

/// blast-rs: Rust ownership equivalent of NCBI `JumperGapAlignFree`; not a direct NCBI C port.
pub fn jumper_gap_align_free(_: Option<JumperGapAlign>) -> Option<JumperGapAlign> {
    None
}

/// NCBI: JumperGapAlignNew (jumper.c:199).
pub fn jumper_gap_align_new(size: i32) -> Option<JumperGapAlign> {
    if size <= 0 {
        return None;
    }
    let left_prelim_block = jumper_prelim_edit_block_new(size)?;
    let right_prelim_block = jumper_prelim_edit_block_new(size)?;
    let mut table = vec![0u32; 256];
    s_create_table(&mut table);
    Some(JumperGapAlign {
        left_prelim_block: Some(left_prelim_block),
        right_prelim_block: Some(right_prelim_block),
        table,
    })
}

/// NCBI: s_ResetJumperPrelimEditBlocks (jumper.c:229).
pub fn s_reset_jumper_prelim_edit_blocks(
    left: Option<&mut JumperPrelimEditBlock>,
    right: Option<&mut JumperPrelimEditBlock>,
) {
    let (Some(left), Some(right)) = (left, right) else {
        return;
    };
    left.edit_ops.clear();
    right.edit_ops.clear();
}

/// NCBI: s_GetSeqPositions (jumper.c:277).
pub fn s_get_seq_positions(
    edit_script: Option<&JumperPrelimEditBlock>,
    edit_index: i32,
    query_pos: &mut i32,
    subject_pos: &mut i32,
) -> i32 {
    let Some(edit_script) = edit_script else {
        return -1;
    };
    if edit_index < 0 || edit_index as usize > edit_script.edit_ops.len() {
        return -1;
    }

    for &op in edit_script.edit_ops.iter().take(edit_index as usize) {
        if op == JUMPER_MISMATCH {
            *query_pos += 1;
            *subject_pos += 1;
        } else if op == JUMPER_INSERTION {
            *query_pos += 1;
        } else if op == JUMPER_DELETION {
            *subject_pos += 1;
        } else {
            *query_pos += op;
            *subject_pos += op;
        }
    }

    0
}

/// blast-rs: Bounds-checked packed-subject base accessor; not a direct NCBI C port.
fn jumper_packed_subject_base(subject: &[u8], subject_pos: i32, subject_length: i32) -> Option<u8> {
    if subject_pos < 0 || subject_pos >= subject_length {
        return None;
    }
    if subject_pos as usize >= subject.len().saturating_mul(4) {
        return None;
    }
    Some(crate::encoding::ncbi2na_base_at(
        subject,
        subject_pos as usize,
    ))
}

/// blast-rs: Jumper base comparison wrapper over packed subject storage; not a direct NCBI C port.
fn jumper_bases_match(
    query: &[u8],
    subject: &[u8],
    query_pos: i32,
    subject_pos: i32,
    query_length: i32,
    subject_length: i32,
) -> bool {
    if query_pos < 0 || query_pos >= query_length {
        return false;
    }
    let Some(&query_base) = query.get(query_pos as usize) else {
        return false;
    };
    let Some(subject_base) = jumper_packed_subject_base(subject, subject_pos, subject_length)
    else {
        return false;
    };
    query_base == subject_base
}

/// NCBI: s_ShiftGapsRight (jumper.c:286).
pub fn s_shift_gaps_right(
    edit_script: &mut JumperPrelimEditBlock,
    query: &[u8],
    subject: &[u8],
    query_offset: i32,
    subject_offset: i32,
    query_length: i32,
    subject_length: i32,
    score: &mut i32,
    err_score: i32,
    num_identical: &mut i32,
) -> i32 {
    let mut i = 0usize;
    while i < edit_script.edit_ops.len() {
        if edit_script.edit_ops[i] < 0 {
            let mut k = i + 1;
            while k < edit_script.edit_ops.len()
                && edit_script.edit_ops[k] == edit_script.edit_ops[i]
            {
                k += 1;
            }

            while k < edit_script.edit_ops.len() {
                if edit_script.edit_ops[k] < 0 && edit_script.edit_ops[k] != edit_script.edit_ops[i]
                {
                    edit_script.edit_ops.remove(k);
                    edit_script.edit_ops[i] = JUMPER_MISMATCH;
                    i += 1;
                    *score -= err_score;

                    let mut q_pos = query_offset;
                    let mut s_pos = subject_offset;
                    let new_index = i as i32 - 1;
                    s_get_seq_positions(Some(edit_script), new_index, &mut q_pos, &mut s_pos);
                    if jumper_bases_match(
                        query,
                        subject,
                        q_pos,
                        s_pos,
                        query_length,
                        subject_length,
                    ) {
                        edit_script.edit_ops[i - 1] = 1;
                        *score += 1;
                        *num_identical += 1;
                    }
                    break;
                }

                if edit_script.edit_ops[k] == JUMPER_MISMATCH {
                    edit_script.edit_ops.swap(i, k);
                    i += 1;

                    let mut q_pos = query_offset;
                    let mut s_pos = subject_offset;
                    let new_index = i as i32 - 1;
                    s_get_seq_positions(Some(edit_script), new_index, &mut q_pos, &mut s_pos);
                    if jumper_bases_match(
                        query,
                        subject,
                        q_pos,
                        s_pos,
                        query_length,
                        subject_length,
                    ) {
                        edit_script.edit_ops[i - 1] = 1;
                        *score -= err_score;
                        *score += 1;
                        *num_identical += 1;
                    }
                }

                if edit_script.edit_ops[k] > 0 {
                    let num_matches = edit_script.edit_ops[k];
                    let mut moved = 0;
                    let mut q_pos = query_offset;
                    let mut s_pos = subject_offset;
                    s_get_seq_positions(Some(edit_script), i as i32, &mut q_pos, &mut s_pos);

                    while moved < edit_script.edit_ops[k]
                        && jumper_bases_match(
                            query,
                            subject,
                            q_pos,
                            s_pos,
                            query_length,
                            subject_length,
                        )
                    {
                        moved += 1;
                        q_pos += 1;
                        s_pos += 1;
                    }

                    if moved == 0 {
                        break;
                    }

                    debug_assert!(i > 0);
                    if i == 0 {
                        break;
                    }

                    if i > 0 && edit_script.edit_ops[i - 1] <= 0 {
                        edit_script.edit_ops.push(JUMPER_MISMATCH);
                        let mut j = edit_script.edit_ops.len() - 1;
                        while j > i {
                            edit_script.edit_ops[j] = edit_script.edit_ops[j - 1];
                            j -= 1;
                        }
                        k += 1;
                        i += 1;
                        edit_script.edit_ops[i - 1] = 0;
                    }

                    if i > 0 {
                        edit_script.edit_ops[i - 1] += moved;
                    }
                    edit_script.edit_ops[k] -= moved;

                    if moved < num_matches {
                        break;
                    }
                    edit_script.edit_ops.remove(k);
                    k = k.saturating_sub(1);
                }

                k += 1;
            }
            i = k.saturating_sub(1);
        }
        i += 1;
    }

    0
}

/// NCBI: s_ShiftGaps (jumper.c:457).
pub fn s_shift_gaps(
    jumper: &mut JumperGapAlign,
    query: &[u8],
    subject: &[u8],
    query_start: i32,
    subject_start: i32,
    query_stop: &mut i32,
    subject_stop: &mut i32,
    score: &mut i32,
    query_length: i32,
    subject_length: i32,
    err_score: i32,
    num_identical: &mut i32,
) -> i32 {
    let Some(left) = jumper.left_prelim_block.as_mut() else {
        return -1;
    };
    let Some(right) = jumper.right_prelim_block.as_ref() else {
        return -1;
    };

    let mut combined = JumperPrelimEditBlock {
        edit_ops: Vec::with_capacity(right.edit_ops.capacity().max(right.edit_ops.len())),
    };
    for &op in left.edit_ops.iter().rev() {
        combined.edit_ops.push(op);
    }
    for &op in &right.edit_ops {
        combined.edit_ops.push(op);
    }

    let mut i = 1usize;
    while i < combined.edit_ops.len() {
        if combined.edit_ops[i - 1] > 0 && combined.edit_ops[i] > 0 {
            combined.edit_ops[i - 1] += combined.edit_ops[i];
            combined.edit_ops.remove(i);
        } else {
            i += 1;
        }
    }

    s_shift_gaps_right(
        &mut combined,
        query,
        subject,
        query_start,
        subject_start,
        query_length,
        subject_length,
        score,
        err_score,
        num_identical,
    );

    while combined.edit_ops.last().is_some_and(|&op| op < 0) {
        let op = combined.edit_ops.pop().unwrap();
        if op == JUMPER_DELETION {
            *subject_stop -= 1;
        } else {
            *query_stop -= 1;
        }
        *score -= err_score;
    }

    left.edit_ops.clear();
    jumper.right_prelim_block = Some(combined);
    0
}

/// NCBI: s_TrimExtension (jumper.c:519).
pub fn s_trim_extension(
    jops: &mut JumperPrelimEditBlock,
    margin: i32,
    cp: &mut i32,
    cq: &mut i32,
    num_identical: &mut i32,
    is_right_ext: bool,
) {
    if jops.edit_ops.is_empty() || margin == 0 {
        return;
    }

    let mut num_matches = 0;
    let mut index = jops.edit_ops.len() as i32 - 1;
    while index >= 1 && jops.edit_ops[index as usize] > 0 {
        num_matches += jops.edit_ops[index as usize];
        index -= 1;
    }

    while jops.edit_ops.len() > 1 && num_matches < margin {
        let op = *jops.edit_ops.last().unwrap();
        if op >= 0 {
            if op > 0 {
                let delta = if is_right_ext { -op } else { op };
                *cp += delta;
                *cq += delta;
                *num_identical -= op;
            } else if is_right_ext {
                *cp -= 1;
                *cq -= 1;
            } else {
                *cp += 1;
                *cq += 1;
            }
        } else if op == JUMPER_INSERTION {
            if is_right_ext {
                *cp -= 1;
            } else {
                *cp += 1;
            }
        } else if op == JUMPER_DELETION {
            if is_right_ext {
                *cq -= 1;
            } else {
                *cq += 1;
            }
        }

        jops.edit_ops.pop();
        if index >= jops.edit_ops.len() as i32 {
            num_matches = 0;
            index = jops.edit_ops.len() as i32 - 1;
            while index >= 1 && jops.edit_ops[index as usize] > 0 {
                num_matches += jops.edit_ops[index as usize];
                index -= 1;
            }
        }
    }

    if jops.edit_ops.len() == 1 && jops.edit_ops[0] <= 0 {
        jops.edit_ops.clear();
    }
}

/// blast-rs: Rust ownership equivalent of NCBI `JumperEditsBlockFree`; not a direct NCBI C port.
pub fn jumper_edits_block_free(_: Option<JumperEditsBlock>) -> Option<JumperEditsBlock> {
    None
}

/// NCBI: JumperEditsBlockNew (jumper.c:2718).
pub fn jumper_edits_block_new(num: i32) -> Option<JumperEditsBlock> {
    if num < 0 {
        return None;
    }
    Some(JumperEditsBlock {
        edits: Vec::with_capacity(num as usize),
    })
}

/// NCBI: JumperEditsBlockDup (jumper.c:2737).
pub fn jumper_edits_block_dup(block: Option<&JumperEditsBlock>) -> Option<JumperEditsBlock> {
    block.cloned()
}

/// blast-rs: Convert one Jumper preliminary operation to the corresponding BLAST gap op; not a direct NCBI C port.
fn jumper_op_to_gap_op(op: JumperOpType) -> GapAlignOpType {
    if op >= 0 {
        GapAlignOpType::Sub
    } else if op == JUMPER_INSERTION {
        GapAlignOpType::Ins
    } else {
        GapAlignOpType::Del
    }
}

/// blast-rs: Convert one Jumper preliminary operation to its run length; not a direct NCBI C port.
fn jumper_op_to_num(op: JumperOpType) -> i32 {
    if op > 0 {
        op
    } else {
        1
    }
}

/// NCBI: JumperPrelimEditBlockToGapEditScript (jumper.c:610).
pub fn jumper_prelim_edit_block_to_gap_edit_script(
    rev_prelim_block: &JumperPrelimEditBlock,
    fwd_prelim_block: &JumperPrelimEditBlock,
) -> Option<GapEditScript> {
    if rev_prelim_block.edit_ops.is_empty() && fwd_prelim_block.edit_ops.is_empty() {
        return None;
    }

    let mut script = GapEditScript::with_capacity(
        rev_prelim_block.edit_ops.len() + fwd_prelim_block.edit_ops.len(),
    );

    for &op in rev_prelim_block.edit_ops.iter().rev() {
        script.push(jumper_op_to_gap_op(op), jumper_op_to_num(op));
    }
    for &op in &fwd_prelim_block.edit_ops {
        script.push(jumper_op_to_gap_op(op), jumper_op_to_num(op));
    }

    Some(script)
}

/// NCBI: JumperEditsBlockCombine (jumper.c:2860).
pub fn jumper_edits_block_combine(
    block: &mut Option<JumperEditsBlock>,
    append: &mut Option<JumperEditsBlock>,
) -> Option<JumperEditsBlock> {
    let Some(base) = block.as_mut() else {
        return None;
    };

    let Some(mut appended) = append.take() else {
        return Some(base.clone());
    };
    if appended.edits.is_empty() {
        return Some(base.clone());
    }

    if base.edits.try_reserve(appended.edits.len()).is_err() {
        *append = Some(appended);
        return None;
    }
    base.edits.append(&mut appended.edits);
    Some(base.clone())
}

/// NCBI: s_ComputeExtensionScore (jumper.c:722).
pub fn s_compute_extension_score(
    edit_script: &JumperPrelimEditBlock,
    match_score: i32,
    mismatch_score: i32,
    gap_open_score: i32,
    gap_extend_score: i32,
) -> i32 {
    let mut score = 0;
    let mut i = 0;
    while i < edit_script.edit_ops.len() {
        let op = edit_script.edit_ops[i];
        if op > 0 {
            score += op * match_score;
        } else if op == JUMPER_MISMATCH {
            score += mismatch_score;
        } else {
            score += gap_open_score;
            while i < edit_script.edit_ops.len() && edit_script.edit_ops[i] == op {
                score += gap_extend_score;
                i += 1;
            }
            i = i.saturating_sub(1);
        }
        i += 1;
    }
    score
}

/// blast-rs: Rust ownership equivalent of NCBI `SequenceOverhangsFree`; not a direct NCBI C port.
pub fn sequence_overhangs_free(_: Option<SequenceOverhangs>) -> Option<SequenceOverhangs> {
    None
}

/// NCBI: JumperGoodAlign (jumper.c:2650).
pub fn jumper_good_align(
    score: i32,
    query_start: i32,
    query_stop: i32,
    subject_start: i32,
    subject_stop: i32,
    num_identical: i32,
    hit_options: &crate::options::HitSavingOptions,
    query_length: i32,
) -> bool {
    let align_len = (query_stop - query_start).max(subject_stop - subject_start);
    if align_len <= 0 {
        return false;
    }

    if 100.0 * num_identical as f64 / (align_len as f64) < hit_options.percent_identity {
        return false;
    }

    if hit_options.splice {
        return true;
    }

    let cutoff_score = if hit_options.cutoff_score_fun[1] != 0 {
        (hit_options.cutoff_score_fun[0] + query_length * hit_options.cutoff_score_fun[1]) / 100
    } else if hit_options.cutoff_score == 0 {
        crate::filter::get_cutoff_score(query_length)
    } else {
        hit_options.cutoff_score
    };
    if score < cutoff_score {
        return false;
    }

    let edit_dist = align_len - num_identical;
    edit_dist <= hit_options.max_edit_distance
}

/// NCBI: JumperFindEdits (jumper.c:2755).
pub fn jumper_find_edits(
    query: &[u8],
    subject: &[u8],
    query_start: i32,
    subject_start: i32,
    query_stop: i32,
    subject_stop: i32,
    left_ext: &JumperPrelimEditBlock,
    right_ext: &JumperPrelimEditBlock,
) -> Option<JumperEditsBlock> {
    const GAP: u8 = 15;
    if query_start < 0
        || subject_start < 0
        || query_stop < query_start
        || subject_stop < subject_start
        || query_stop as usize > query.len()
        || subject_stop as usize > subject.len().saturating_mul(4)
    {
        return None;
    }
    let mut q_pos = query_start;
    let mut s_pos = subject_start;
    let mut edits = JumperEditsBlock {
        edits: Vec::with_capacity(left_ext.edit_ops.len() + right_ext.edit_ops.len()),
    };

    for &op in left_ext.edit_ops.iter().rev() {
        match op {
            JUMPER_MISMATCH => {
                edits.edits.push(JumperEdit {
                    query_pos: q_pos,
                    query_base: *query.get(q_pos as usize)?,
                    subject_base: crate::encoding::ncbi2na_base_at(subject, s_pos as usize),
                });
                q_pos += 1;
                s_pos += 1;
            }
            JUMPER_DELETION => {
                edits.edits.push(JumperEdit {
                    query_pos: q_pos,
                    query_base: GAP,
                    subject_base: crate::encoding::ncbi2na_base_at(subject, s_pos as usize),
                });
                s_pos += 1;
            }
            JUMPER_INSERTION => {
                edits.edits.push(JumperEdit {
                    query_pos: q_pos,
                    query_base: *query.get(q_pos as usize)?,
                    subject_base: GAP,
                });
                q_pos += 1;
            }
            run if run > 0 => {
                q_pos += run;
                s_pos += run;
            }
            _ => return None,
        }
    }

    for &op in &right_ext.edit_ops {
        match op {
            JUMPER_MISMATCH => {
                edits.edits.push(JumperEdit {
                    query_pos: q_pos,
                    query_base: *query.get(q_pos as usize)?,
                    subject_base: crate::encoding::ncbi2na_base_at(subject, s_pos as usize),
                });
                q_pos += 1;
                s_pos += 1;
            }
            JUMPER_DELETION => {
                edits.edits.push(JumperEdit {
                    query_pos: q_pos,
                    query_base: GAP,
                    subject_base: crate::encoding::ncbi2na_base_at(subject, s_pos as usize),
                });
                s_pos += 1;
            }
            JUMPER_INSERTION => {
                edits.edits.push(JumperEdit {
                    query_pos: q_pos,
                    query_base: *query.get(q_pos as usize)?,
                    subject_base: GAP,
                });
                q_pos += 1;
            }
            run if run > 0 => {
                q_pos += run;
                s_pos += run;
            }
            _ => return None,
        }
    }

    if q_pos != query_stop || s_pos != subject_stop {
        return None;
    }
    Some(edits)
}

/// NCBI: JumperFindSpliceSignals (jumper.c:2945).
pub fn jumper_find_splice_signals(
    hsp: &crate::hspstream::Hsp,
    map_info: &mut crate::hspstream::BlastHSPMappingInfo,
    query_len: i32,
    subject: &[u8],
    subject_len: i32,
) -> i32 {
    if query_len < 0
        || subject_len < 0
        || subject_len as usize > subject.len().saturating_mul(4)
        || hsp.query_offset < 0
        || hsp.query_end < hsp.query_offset
        || hsp.query_end > query_len
        || hsp.subject_offset < 0
        || hsp.subject_end < hsp.subject_offset
        || hsp.subject_end > subject_len
    {
        return -1;
    }

    if hsp.query_offset == 0 || hsp.subject_offset < 2 {
        map_info.left_edge = MAPPER_EXON;
    } else {
        map_info.left_edge =
            (crate::encoding::ncbi2na_base_at(subject, hsp.subject_offset as usize - 2) << 2)
                | crate::encoding::ncbi2na_base_at(subject, hsp.subject_offset as usize - 1);
    }

    if hsp.query_end == query_len || hsp.subject_end == subject_len {
        map_info.right_edge = MAPPER_EXON;
    } else {
        map_info.right_edge = (crate::encoding::ncbi2na_base_at(subject, hsp.subject_end as usize)
            << 2)
            | crate::encoding::ncbi2na_base_at(subject, hsp.subject_end as usize + 1);
    }
    0
}

/// NCBI: s_SaveSubjectOverhangs (jumper.c:2995).
pub fn s_save_subject_overhangs(
    hsp: &crate::hspstream::Hsp,
    map_info: &mut crate::hspstream::BlastHSPMappingInfo,
    subject: &[u8],
    subject_len: i32,
    query_len: i32,
) -> i32 {
    if query_len < 0
        || subject_len < 0
        || subject_len as usize > subject.len().saturating_mul(4)
        || hsp.subject_offset < 0
        || hsp.subject_end < 0
        || hsp.subject_offset > subject_len
        || hsp.subject_end > subject_len
    {
        return -1;
    }

    let max_subject_overhang = if query_len < 400 { 30 } else { 60 };
    let mut overhangs = SequenceOverhangs::default();

    if hsp.query_offset >= 0 {
        let mut len = hsp.query_offset.max(2).min(max_subject_overhang);
        if hsp.subject_offset < len {
            len = hsp.subject_offset;
        }
        if len > 0 {
            let start = hsp.subject_offset.saturating_sub(len);
            let left = (0..len)
                .map(|i| {
                    crate::encoding::ncbi2na_base_at(subject, start.saturating_add(i) as usize)
                })
                .collect::<Vec<_>>();
            overhangs.left = Some(left);
        }
    }

    if hsp.query_end <= query_len {
        let query_right = query_len.saturating_sub(hsp.query_end).saturating_add(1);
        let len = if query_right < 6 {
            query_right.max(2).min(max_subject_overhang)
        } else {
            max_subject_overhang
        }
        .min(subject_len.saturating_sub(hsp.subject_end));

        if len > 0 {
            let right = (0..len)
                .map(|i| {
                    crate::encoding::ncbi2na_base_at(
                        subject,
                        hsp.subject_end.saturating_add(i) as usize,
                    )
                })
                .collect::<Vec<_>>();
            overhangs.right = Some(right);
        }
    }

    map_info.subject_overhangs = Some(overhangs);
    0
}

/// blast-rs: Port-shaped equivalent of NCBI `s_CreateHSPForWordHit` (jumper.c:3098); not a direct NCBI C port.
///
/// C stores mapper data inside `BlastHSP::map_info`; Rust keeps it as a
/// separate [`crate::hspstream::BlastHSPMappingInfo`], so this returns both.
pub fn s_create_hsp_for_word_hit(
    q_offset: i32,
    s_offset: i32,
    length: i32,
    context: i32,
    query: &[u8],
    query_frame: i32,
    subject: &[u8],
    subject_len: i32,
    subject_frame: i32,
    query_len: i32,
) -> Option<(crate::hspstream::Hsp, crate::hspstream::BlastHSPMappingInfo)> {
    if q_offset < 0 || s_offset < 0 || length < 0 {
        return None;
    }
    let q_end = q_offset.checked_add(length)?;
    let s_end = s_offset.checked_add(length)?;
    if q_end as usize > query.len() || s_end > subject_len {
        return None;
    }
    if subject_len < 0 || subject_len as usize > subject.len().saturating_mul(4) {
        return None;
    }

    let mut edit_script = GapEditScript::new();
    edit_script.ops.push((GapAlignOpType::Sub, length));
    let mut hsp = crate::hspstream::Hsp {
        score: length,
        num_ident: length,
        bit_score: 0.0,
        evalue: 0.0,
        query_offset: q_offset,
        query_end: q_end,
        query_gapped_start: q_offset,
        subject_offset: s_offset,
        subject_end: s_end,
        subject_gapped_start: s_offset,
        context,
        query_frame,
        subject_frame,
        num_gaps: 0,
        comp_adjustment_method: 0,
        edit_script: Some(edit_script),
        pat_info: None,
        map_info: None,
    };
    let mut map_info = crate::hspstream::blast_hsp_mapping_info_new();
    let mut edits = JumperEditsBlock::default();
    for i in 0..length {
        let query_base = query[(q_offset + i) as usize];
        if query_base & 0xfc != 0 {
            edits.edits.push(JumperEdit {
                query_pos: q_offset + i,
                query_base,
                subject_base: crate::encoding::ncbi2na_base_at(subject, (s_offset + i) as usize),
            });
        }
    }
    hsp.num_ident = length - edits.edits.len() as i32;
    map_info.edits = Some(edits);
    jumper_find_splice_signals(&hsp, &mut map_info, query_len, subject, subject_len);
    s_save_subject_overhangs(&hsp, &mut map_info, subject, subject_len, query_len);
    hsp.map_info = Some(map_info.clone());
    Some((hsp, map_info))
}

/// blast-rs: Port-shaped equivalent of NCBI `s_CreateHSP` (jumper.c:3188); not a direct NCBI C port.
///
/// Rust keeps C `BlastHSP::map_info` as a separate value, so this returns the
/// created HSP with its mapping info side by side.
pub fn s_create_hsp(
    query_seq: &[u8],
    query_len: i32,
    context: i32,
    query_frame: i32,
    subject: &[u8],
    subject_len: i32,
    subject_frame: i32,
    jumper: &mut JumperGapAlign,
    query_start: i32,
    mut query_stop: i32,
    subject_start: i32,
    mut subject_stop: i32,
    mut score: i32,
    scoring_penalty: i32,
    splice: bool,
) -> Option<(crate::hspstream::Hsp, crate::hspstream::BlastHSPMappingInfo)> {
    if query_len < 0
        || subject_len < 0
        || query_len as usize > query_seq.len()
        || subject_len as usize > subject.len().saturating_mul(4)
    {
        return None;
    }

    let mut num_identical = 0;
    if std::env::var_os("MAPPER_NO_GAP_SHIFT").is_none() {
        s_shift_gaps(
            jumper,
            query_seq,
            subject,
            query_start,
            subject_start,
            &mut query_stop,
            &mut subject_stop,
            &mut score,
            query_len,
            subject_len,
            scoring_penalty,
            &mut num_identical,
        );
    }

    let left = jumper.left_prelim_block.as_ref()?;
    let right = jumper.right_prelim_block.as_ref()?;
    if query_start < 0
        || subject_start < 0
        || query_stop < query_start
        || subject_stop < subject_start
        || query_stop > query_len
        || subject_stop > subject_len
    {
        return None;
    }
    let edit_script = jumper_prelim_edit_block_to_gap_edit_script(left, right);
    let mut hsp = crate::hspstream::Hsp {
        score,
        num_ident: num_identical,
        bit_score: 0.0,
        evalue: 0.0,
        query_offset: query_start,
        query_end: query_stop,
        query_gapped_start: query_start,
        subject_offset: subject_start,
        subject_end: subject_stop,
        subject_gapped_start: subject_start,
        context,
        query_frame,
        subject_frame,
        num_gaps: 0,
        comp_adjustment_method: 0,
        edit_script,
        pat_info: None,
        map_info: None,
    };
    let mut map_info = crate::hspstream::blast_hsp_mapping_info_new();
    map_info.edits = jumper_find_edits(
        query_seq,
        subject,
        query_start,
        subject_start,
        query_stop,
        subject_stop,
        left,
        right,
    );

    if splice {
        jumper_find_splice_signals(&hsp, &mut map_info, query_len, subject, subject_len);
        s_save_subject_overhangs(&hsp, &mut map_info, subject, subject_len, query_len);
    }

    hsp.evalue = 0.0;
    hsp.map_info = Some(map_info.clone());
    Some((hsp, map_info))
}

#[allow(clippy::too_many_arguments)]
/// blast-rs: Indexed short-read HSP constructor for Rust mapper state; not a direct NCBI C port.
fn s_create_short_read_indexed_hsp(
    query_seq: &[u8],
    query_len: i32,
    context: i32,
    query_frame: i32,
    subject: &[u8],
    subject_len: i32,
    subject_frame: i32,
    jumper: &JumperGapAlign,
    query_start: i32,
    query_stop: i32,
    subject_start: i32,
    subject_stop: i32,
    seed_query_offset: i32,
    seed_subject_offset: i32,
    score: i32,
    num_identical: i32,
    splice: bool,
) -> Option<crate::hspstream::Hsp> {
    if query_len < 0
        || subject_len < 0
        || query_len as usize > query_seq.len()
        || subject_len as usize > subject.len().saturating_mul(4)
        || query_start < 0
        || subject_start < 0
        || query_stop < query_start
        || subject_stop < subject_start
        || query_stop > query_len
        || subject_stop > subject_len
        || seed_query_offset < 0
        || seed_subject_offset < 0
    {
        return None;
    }

    let left = jumper.left_prelim_block.as_ref()?;
    let right = jumper.right_prelim_block.as_ref()?;
    let edit_script = jumper_prelim_edit_block_to_gap_edit_script(left, right);
    let mut hsp = crate::hspstream::Hsp {
        score,
        num_ident: num_identical,
        bit_score: 0.0,
        evalue: 0.0,
        query_offset: query_start,
        query_end: query_stop,
        query_gapped_start: seed_query_offset,
        subject_offset: subject_start,
        subject_end: subject_stop,
        subject_gapped_start: seed_subject_offset,
        context,
        query_frame,
        subject_frame,
        num_gaps: 0,
        comp_adjustment_method: 0,
        edit_script,
        pat_info: None,
        map_info: None,
    };
    let mut map_info = crate::hspstream::blast_hsp_mapping_info_new();
    map_info.edits = jumper_find_edits(
        query_seq,
        subject,
        query_start,
        subject_start,
        query_stop,
        subject_stop,
        left,
        right,
    );

    if splice {
        jumper_find_splice_signals(&hsp, &mut map_info, query_len, subject, subject_len);
    }

    hsp.map_info = Some(map_info);
    Some(hsp)
}

/// blast-rs: Mapper rescue base matcher with ambiguous-query handling; not a direct NCBI C port.
fn jumper_base_matches_or_ambiguous(query_base: u8, subject: &[u8], subject_pos: i32) -> bool {
    query_base & 0xfc != 0
        || (subject_pos >= 0
            && (subject_pos as usize) < subject.len().saturating_mul(4)
            && query_base == crate::encoding::ncbi2na_base_at(subject, subject_pos as usize))
}

#[allow(clippy::too_many_arguments)]
/// blast-rs: Optional small-word rescue around saved Jumper HSPs; not a direct NCBI C port.
fn blast_na_extend_jumper_small_word_rescue(
    saved_hsp: &crate::hspstream::Hsp,
    s_index: &SubjectIndex,
    context: i32,
    query_slice: &[u8],
    query_frame: i32,
    query_len: i32,
    subject: &[u8],
    subject_length: i32,
    subject_frame: i32,
    score_params: &crate::parameters::ScoringParameters,
    longest_intron: i32,
    hsp_list: &mut crate::hspstream::HspList,
) -> i32 {
    const SUBJECT_INDEX_WORD_LENGTH: i32 = 4;
    const K_MIN_SUBJECT_OVERHANG: i32 = 100;
    if query_len < 0
        || query_len as usize > query_slice.len()
        || subject_length < 0
        || subject_length as usize > subject.len().saturating_mul(4)
        || saved_hsp.query_offset < 0
        || saved_hsp.query_end < saved_hsp.query_offset
        || saved_hsp.query_end > query_len
        || saved_hsp.subject_offset < 0
        || saved_hsp.subject_end < saved_hsp.subject_offset
        || saved_hsp.subject_end > subject_length
    {
        return 0;
    }
    let mut saved = 0;

    let right_query_overhang = query_len.saturating_sub(saved_hsp.query_end);
    let right_subject_overhang = subject_length.saturating_sub(saved_hsp.subject_end);
    if right_query_overhang < 16
        && right_query_overhang >= SUBJECT_INDEX_WORD_LENGTH
        && right_subject_overhang >= K_MIN_SUBJECT_OVERHANG
    {
        if right_query_overhang < -score_params.penalty {
            let mut i = 1;
            while i < right_query_overhang {
                if query_slice[(saved_hsp.query_end.saturating_add(i)) as usize]
                    != crate::encoding::ncbi2na_base_at(
                        subject,
                        saved_hsp.subject_end.saturating_add(i) as usize,
                    )
                {
                    break;
                }
                i += 1;
            }
            if i > SUBJECT_INDEX_WORD_LENGTH || i == right_query_overhang {
                if let Some((word_hsp, _)) = s_create_hsp_for_word_hit(
                    saved_hsp.query_end.saturating_add(1),
                    saved_hsp.subject_end.saturating_add(1),
                    i.saturating_sub(1),
                    context,
                    query_slice,
                    query_frame,
                    subject,
                    subject_length,
                    subject_frame,
                    query_len,
                ) {
                    if crate::hspstream::blast_hsp_list_save_hsp(hsp_list, word_hsp) != 0 {
                        return -1;
                    }
                    saved += 1;
                }
            }
        }

        let mut q = saved_hsp.query_end;
        let mut round = 0;
        while q + SUBJECT_INDEX_WORD_LENGTH <= query_len && round < 2 {
            while q + SUBJECT_INDEX_WORD_LENGTH <= query_len {
                let mut i = 0;
                while i < SUBJECT_INDEX_WORD_LENGTH {
                    if query_slice[(q + i) as usize] & 0xfc != 0 {
                        q += i + 1;
                        break;
                    }
                    i += 1;
                }
                if i == SUBJECT_INDEX_WORD_LENGTH {
                    break;
                }
            }
            if q.saturating_add(SUBJECT_INDEX_WORD_LENGTH)
                .saturating_sub(1)
                >= query_len
            {
                break;
            }
            let Some(word) = unpacked_query_word(query_slice, q, SUBJECT_INDEX_WORD_LENGTH) else {
                break;
            };
            let from = saved_hsp.subject_end;
            let query_tail = query_len.saturating_sub(q).saturating_add(1);
            let to = from
                .saturating_add(longest_intron)
                .min(subject_length.saturating_sub(query_tail));
            if to >= from {
                if let Some(mut it) = subject_index_iterator_new(s_index, word, from, to) {
                    let mut s_pos = subject_index_iterator_next(&mut it);
                    while s_pos >= 0 {
                        let mut qt = q;
                        let mut st = s_pos;
                        while qt < query_len
                            && st < subject_length
                            && jumper_base_matches_or_ambiguous(
                                query_slice[qt as usize],
                                subject,
                                st,
                            )
                        {
                            qt += 1;
                            st += 1;
                        }
                        if qt == query_len {
                            let mut qf = q;
                            let mut sf = s_pos;
                            while qf >= 0
                                && sf >= 0
                                && jumper_base_matches_or_ambiguous(
                                    query_slice[qf as usize],
                                    subject,
                                    sf,
                                )
                            {
                                qf -= 1;
                                sf -= 1;
                            }
                            qf += 1;
                            sf += 1;
                            if qt - qf >= 5
                                && qf <= saved_hsp.query_end.saturating_add(1)
                                && qf > saved_hsp.query_offset
                            {
                                while qf < qt && query_slice[qf as usize] & 0xfc != 0 {
                                    qf += 1;
                                    sf += 1;
                                }
                                while qt > qf && query_slice[(qt - 1) as usize] & 0xfc != 0 {
                                    qt -= 1;
                                }
                                if qf < qt {
                                    if let Some((word_hsp, _)) = s_create_hsp_for_word_hit(
                                        qf,
                                        sf,
                                        qt - qf,
                                        context,
                                        query_slice,
                                        query_frame,
                                        subject,
                                        subject_length,
                                        subject_frame,
                                        query_len,
                                    ) {
                                        if crate::hspstream::blast_hsp_list_save_hsp(
                                            hsp_list, word_hsp,
                                        ) != 0
                                        {
                                            return -1;
                                        }
                                        saved += 1;
                                    }
                                }
                            }
                        }
                        s_pos = subject_index_iterator_next(&mut it);
                    }
                }
            }
            q += 1;
            round += 1;
        }
    }

    if saved_hsp.query_offset < 16
        && saved_hsp.query_offset >= SUBJECT_INDEX_WORD_LENGTH
        && saved_hsp.subject_offset >= K_MIN_SUBJECT_OVERHANG
    {
        if saved_hsp.query_offset < -score_params.penalty {
            let mut i = 2;
            let left_query_limit = saved_hsp.query_offset.saturating_sub(1);
            while i < left_query_limit {
                if query_slice[(saved_hsp.query_offset.saturating_sub(i)) as usize]
                    != crate::encoding::ncbi2na_base_at(
                        subject,
                        saved_hsp.subject_offset.saturating_sub(i) as usize,
                    )
                {
                    break;
                }
                i += 1;
            }
            if i > SUBJECT_INDEX_WORD_LENGTH || i == left_query_limit {
                if let Some((word_hsp, _)) = s_create_hsp_for_word_hit(
                    saved_hsp.query_offset.saturating_sub(1).saturating_sub(i),
                    saved_hsp.subject_offset.saturating_sub(1).saturating_sub(i),
                    i,
                    context,
                    query_slice,
                    query_frame,
                    subject,
                    subject_length,
                    subject_frame,
                    query_len,
                ) {
                    if crate::hspstream::blast_hsp_list_save_hsp(hsp_list, word_hsp) != 0 {
                        return -1;
                    }
                    saved += 1;
                }
            }
        }

        let mut q = saved_hsp
            .query_offset
            .saturating_sub(SUBJECT_INDEX_WORD_LENGTH)
            .max(0);
        let mut round = 0;
        while q >= 0 && round < 2 {
            while q >= 0 {
                let mut i = 0;
                while i < SUBJECT_INDEX_WORD_LENGTH {
                    if query_slice[(q + i) as usize] & 0xfc != 0 {
                        q = q
                            .saturating_add(i)
                            .saturating_sub(SUBJECT_INDEX_WORD_LENGTH);
                        break;
                    }
                    i += 1;
                }
                if i == SUBJECT_INDEX_WORD_LENGTH {
                    break;
                }
            }
            if q < 0 {
                break;
            }
            let Some(word) = unpacked_query_word(query_slice, q, SUBJECT_INDEX_WORD_LENGTH) else {
                break;
            };
            let from = saved_hsp
                .subject_offset
                .saturating_sub(SUBJECT_INDEX_WORD_LENGTH);
            let to = from.saturating_sub(longest_intron).max(q + 1);
            if let Some(mut it) = subject_index_iterator_new(s_index, word, from, to) {
                let mut s_pos = subject_index_iterator_prev(&mut it);
                while s_pos >= 0 {
                    let mut k = q;
                    let mut s = s_pos;
                    while k >= 0
                        && s >= 0
                        && jumper_base_matches_or_ambiguous(query_slice[k as usize], subject, s)
                    {
                        k -= 1;
                        s -= 1;
                    }
                    if k == -1 {
                        k += 1;
                        s += 1;
                        let mut qt = q;
                        let mut st = s_pos;
                        while qt < query_len
                            && st < subject_length
                            && jumper_base_matches_or_ambiguous(
                                query_slice[qt as usize],
                                subject,
                                st,
                            )
                        {
                            qt += 1;
                            st += 1;
                        }
                        if qt - k >= 5
                            && s < saved_hsp.subject_offset
                            && saved_hsp.query_offset <= qt + 1
                        {
                            while k < qt && query_slice[k as usize] & 0xfc != 0 {
                                k += 1;
                                s += 1;
                            }
                            while qt > k && query_slice[(qt - 1) as usize] & 0xfc != 0 {
                                qt -= 1;
                            }
                            if k < qt {
                                if let Some((word_hsp, _)) = s_create_hsp_for_word_hit(
                                    k,
                                    s,
                                    qt - k,
                                    context,
                                    query_slice,
                                    query_frame,
                                    subject,
                                    subject_length,
                                    subject_frame,
                                    query_len,
                                ) {
                                    if crate::hspstream::blast_hsp_list_save_hsp(hsp_list, word_hsp)
                                        != 0
                                    {
                                        return -1;
                                    }
                                    saved += 1;
                                }
                            }
                        }
                    }
                    s_pos = subject_index_iterator_prev(&mut it);
                }
            }
            q -= 1;
            round += 1;
        }
    }

    saved
}

/// blast-rs: Rust ownership equivalent of NCBI `SubjectIndexIteratorFree`; not a direct NCBI C port.
pub fn subject_index_iterator_free(
    _: Option<SubjectIndexIterator>,
) -> Option<SubjectIndexIterator> {
    None
}

/// blast-rs: Packed-subject word extraction for Rust `SubjectIndex`; not a direct NCBI C port.
fn subject_word_at(subject: &[u8], subject_len: i32, pos: i32, word_size: i32) -> Option<u32> {
    if pos < 0 || word_size <= 0 {
        return None;
    }
    let end = pos.checked_add(word_size)?;
    if end > subject_len {
        return None;
    }
    if end as usize > subject.len().saturating_mul(4) {
        return None;
    }
    let mut word = 0u32;
    for i in 0..word_size {
        word = (word << 2)
            | crate::encoding::ncbi2na_base_at(subject, pos.saturating_add(i) as usize) as u32;
    }
    Some(word)
}

/// blast-rs: Port-shaped equivalent of NCBI `SubjectIndexNew` (jumper.c:3874); not a direct NCBI C port.
pub fn subject_index_new(
    subject: &[u8],
    subject_len: i32,
    width: i32,
    word_size: i32,
) -> Option<SubjectIndex> {
    if subject_len < 0 || width <= 0 || word_size <= 0 {
        return None;
    }
    let mut positions_by_word: std::collections::BTreeMap<u32, Vec<i32>> =
        std::collections::BTreeMap::new();
    let max_pos = subject_len.saturating_sub(word_size);
    for pos in 0..=max_pos.max(-1) {
        let word = subject_word_at(subject, subject_len, pos, word_size)?;
        positions_by_word.entry(word).or_default().push(pos);
    }
    Some(SubjectIndex {
        num_lookups: (subject_len / width).saturating_add(1),
        width,
        word_size,
        positions_by_word,
    })
}

/// NCBI: SubjectIndexIteratorNew (jumper.c:3982).
pub fn subject_index_iterator_new(
    s_index: &SubjectIndex,
    word: u32,
    from: i32,
    to: i32,
) -> Option<SubjectIndexIterator> {
    let positions = s_index
        .positions_by_word
        .get(&word)
        .cloned()
        .unwrap_or_default();
    let start = positions.partition_point(|&pos| pos < from) as i32;
    Some(SubjectIndexIterator {
        subject_index: Some(s_index.clone()),
        to,
        lookup_index: if s_index.width > 0 {
            from / s_index.width
        } else {
            0
        },
        positions,
        word_index: start,
    })
}

/// NCBI: SubjectIndexIteratorNext (jumper.c:4035).
pub fn subject_index_iterator_next(it: &mut SubjectIndexIterator) -> i32 {
    if it.word_index < 0 || it.word_index as usize >= it.positions.len() {
        return -1;
    }
    let pos = it.positions[it.word_index as usize];
    if pos > it.to {
        return -1;
    }
    it.word_index += 1;
    pos
}

/// NCBI: SubjectIndexIteratorPrev (jumper.c:4087).
pub fn subject_index_iterator_prev(it: &mut SubjectIndexIterator) -> i32 {
    if it.positions.is_empty() {
        return -1;
    }
    if it.word_index < 0 {
        return -1;
    }
    if it.word_index as usize >= it.positions.len() {
        it.word_index = it.positions.len() as i32 - 1;
    }
    let pos = it.positions[it.word_index as usize];
    if pos < it.to {
        return -1;
    }
    it.word_index -= 1;
    pos
}

/// NCBI: GapEditScriptCombine (jumper.c:2895).
pub fn gap_edit_script_combine(
    edit_script: &mut Option<GapEditScript>,
    append: &mut Option<GapEditScript>,
) -> Option<GapEditScript> {
    let Some(mut append_script) = append.take() else {
        return edit_script.clone();
    };
    if append_script.ops.is_empty() {
        return edit_script.clone();
    }
    let Some(script) = edit_script.as_mut() else {
        *edit_script = Some(append_script);
        return edit_script.clone();
    };

    for (op, count) in append_script.ops.drain(..) {
        script.push(op, count);
    }
    edit_script.clone()
}

/// blast-rs: Bit-mask helper for Jumper mismatch tracing; not a direct NCBI C port.
fn jumper_trace_mask(max_mismatches: i32) -> u32 {
    if max_mismatches <= 0 {
        0
    } else if max_mismatches >= 32 {
        u32::MAX
    } else {
        (1u32 << max_mismatches) - 1
    }
}

/// blast-rs: Bit-window update helper for Jumper mismatch tracing; not a direct NCBI C port.
fn jumper_shift_trace(trace: &mut u32, shift: i32, window: i32) {
    if *trace == 0 {
        return;
    }
    if shift < window {
        *trace = trace.wrapping_shl(shift.max(0) as u32);
    } else {
        *trace = 0;
    }
}

/// blast-rs: Unpacked-subject rightward jump validation helper; not a direct NCBI C port.
fn jumper_match_run_right(
    query: &[u8],
    subject: &[u8],
    mut cp: i32,
    mut cq: i32,
    len: i32,
    ok: i32,
) -> bool {
    let mut n = 0;
    let mut i = len;
    while i > 0 {
        if cp < 0
            || cq < 0
            || cp as usize >= query.len()
            || cq as usize >= subject.len()
            || query[cp as usize] != subject[cq as usize]
        {
            n += 1;
            if n > ok {
                return false;
            }
        }
        cp += 1;
        cq += 1;
        i -= 1;
    }
    true
}

/// blast-rs: Unpacked-subject leftward jump validation helper; not a direct NCBI C port.
fn jumper_match_run_left(
    query: &[u8],
    subject: &[u8],
    mut cp: i32,
    mut cq: i32,
    len: i32,
    ok: i32,
) -> bool {
    let mut n = 0;
    let mut i = len;
    while i > 0 {
        if cp < 0
            || cq < 0
            || cp as usize >= query.len()
            || cq as usize >= subject.len()
            || query[cp as usize] != subject[cq as usize]
        {
            n += 1;
            if n > ok {
                return false;
            }
        }
        cp -= 1;
        cq -= 1;
        i -= 1;
    }
    true
}

/// blast-rs: Finds the next rightward Jumper table entry for unpacked subject data; not a direct NCBI C port.
fn jumper_find_right(query: &[u8], subject: &[u8], cp: i32, cq: i32, jumps: &[Jump]) -> Jump {
    let cpmax = query.len() as i32;
    let cqmax = subject.len() as i32;
    for &jump in jumps {
        if jump.lng == 0 {
            return jump;
        }

        let cp1 = cp + jump.dcp;
        let cq1 = cq + jump.dcq;
        if cp1 >= cpmax || cq1 >= cqmax {
            continue;
        }
        if !jumper_match_run_right(query, subject, cp1, cq1, jump.ok, 0) {
            continue;
        }
        if jump.lng + cp1 >= cpmax || jump.lng + cq1 >= cqmax {
            continue;
        }
        if !jumper_match_run_right(query, subject, cp1, cq1, jump.lng, jump.ok) {
            continue;
        }
        return jump;
    }
    jumps.last().copied().unwrap_or(Jump {
        dcp: 1,
        dcq: 1,
        lng: 0,
        ok: 0,
    })
}

/// blast-rs: Finds the next leftward Jumper table entry for unpacked subject data; not a direct NCBI C port.
fn jumper_find_left(query: &[u8], subject: &[u8], cp: i32, cq: i32, jumps: &[Jump]) -> Jump {
    for &jump in jumps {
        if jump.lng == 0 {
            return jump;
        }

        let cp1 = cp - jump.dcp;
        let cq1 = cq - jump.dcq;
        if cp1 < 0 || cq1 < 0 {
            continue;
        }
        if !jumper_match_run_left(query, subject, cp1, cq1, jump.ok, 0) {
            continue;
        }
        if cp1 - jump.lng < 0 || cq1 - jump.lng < 0 {
            continue;
        }
        if !jumper_match_run_left(query, subject, cp1, cq1, jump.lng, jump.ok) {
            continue;
        }
        return jump;
    }
    jumps.last().copied().unwrap_or(Jump {
        dcp: 1,
        dcq: 1,
        lng: 0,
        ok: 0,
    })
}

/// blast-rs: Packed-subject rightward jump validation helper; not a direct NCBI C port.
fn jumper_packed_match_run_right(
    query: &[u8],
    subject: &[u8],
    mut cp: i32,
    mut cq: i32,
    len: i32,
    ok: i32,
    subject_length: i32,
) -> bool {
    let mut n = 0;
    let mut i = len;
    while i > 0 {
        if cp < 0
            || cq < 0
            || cp as usize >= query.len()
            || cq >= subject_length
            || query[cp as usize] != crate::encoding::ncbi2na_base_at(subject, cq as usize)
        {
            n += 1;
            if n > ok {
                return false;
            }
        }
        cp += 1;
        cq += 1;
        i -= 1;
    }
    true
}

/// blast-rs: Packed-subject leftward jump validation helper; not a direct NCBI C port.
fn jumper_packed_match_run_left(
    query: &[u8],
    subject: &[u8],
    mut cp: i32,
    mut cq: i32,
    len: i32,
    ok: i32,
) -> bool {
    let mut n = 0;
    let mut i = len;
    while i > 0 {
        if cp < 0
            || cq < 0
            || cp as usize >= query.len()
            || query[cp as usize] != crate::encoding::ncbi2na_base_at(subject, cq as usize)
        {
            n += 1;
            if n > ok {
                return false;
            }
        }
        cp -= 1;
        cq -= 1;
        i -= 1;
    }
    true
}

/// blast-rs: Finds the next rightward Jumper table entry for packed subject data; not a direct NCBI C port.
fn jumper_find_right_compressed(
    query: &[u8],
    subject: &[u8],
    subject_length: i32,
    cp: i32,
    cq: i32,
    jumps: &[Jump],
) -> Jump {
    let cpmax = query.len() as i32;
    for &jump in jumps {
        if jump.lng == 0 {
            return jump;
        }

        let cp1 = cp + jump.dcp;
        let cq1 = cq + jump.dcq;
        if cp1 >= cpmax || cq1 >= subject_length {
            continue;
        }
        if !jumper_packed_match_run_right(query, subject, cp1, cq1, jump.ok, 0, subject_length) {
            continue;
        }
        if jump.lng + cp1 >= cpmax || jump.lng + cq1 >= subject_length {
            continue;
        }
        if !jumper_packed_match_run_right(
            query,
            subject,
            cp1,
            cq1,
            jump.lng,
            jump.ok,
            subject_length,
        ) {
            continue;
        }
        return jump;
    }
    jumps.last().copied().unwrap_or(Jump {
        dcp: 1,
        dcq: 1,
        lng: 0,
        ok: 0,
    })
}

/// blast-rs: Finds the next leftward Jumper table entry for packed subject data; not a direct NCBI C port.
fn jumper_find_left_compressed(
    query: &[u8],
    subject: &[u8],
    cp: i32,
    cq: i32,
    jumps: &[Jump],
) -> Jump {
    for &jump in jumps {
        if jump.lng == 0 {
            return jump;
        }

        let cp1 = cp - jump.dcp;
        let cq1 = cq - jump.dcq;
        if cp1 < 0 || cq1 < 0 {
            continue;
        }
        if !jumper_packed_match_run_left(query, subject, cp1, cq1, jump.ok, 0) {
            continue;
        }
        if cp1 - jump.lng < 0 || cq1 - jump.lng < 0 {
            continue;
        }
        if !jumper_packed_match_run_left(query, subject, cp1, cq1, jump.lng, jump.ok) {
            continue;
        }
        return jump;
    }
    jumps.last().copied().unwrap_or(Jump {
        dcp: 1,
        dcq: 1,
        lng: 0,
        ok: 0,
    })
}

/// NCBI: JumperExtendRight (jumper.c:1384).
pub fn jumper_extend_right(
    query: &[u8],
    subject: &[u8],
    match_score: i32,
    mismatch_score: i32,
    gap_open_score: i32,
    gap_extend_score: i32,
    max_mismatches: i32,
    window: i32,
    edit_script: &mut GapPrelimEditBlock,
    left_extension: bool,
) -> (i32, i32, i32) {
    if query.is_empty() || subject.is_empty() {
        return (-1, 0, 0);
    }

    let mut cp = 1;
    let mut cq = 1;
    let mut score = 0;
    let mut num_mismatches = 0;
    let mut new_matches = if left_extension { 0 } else { 1 };
    let mut trace = 0u32;
    let trace_mask = jumper_trace_mask(max_mismatches);

    while cp < query.len() as i32 && cq < subject.len() as i32 && num_mismatches < max_mismatches {
        if query[cp as usize] == subject[cq as usize] {
            score += match_score;
            cp += 1;
            cq += 1;
            new_matches += 1;
            continue;
        }

        let jump = jumper_find_right(query, subject, cp, cq, &JUMPER_DEFAULT);
        if new_matches != 0 {
            gap_prelim_edit_block_add(edit_script, GapAlignOpType::Sub, new_matches);
            jumper_shift_trace(&mut trace, new_matches, window);
            new_matches = 0;
        }

        if jump.dcp == jump.dcq {
            score += mismatch_score * jump.dcp;
            if trace & trace_mask != 0 {
                num_mismatches += jump.dcp;
                trace = trace.wrapping_shl(jump.dcp as u32);
                trace |= (1u32 << jump.dcp) - 1;
            } else {
                num_mismatches = jump.dcp;
                trace = (1u32 << jump.dcp) - 1;
            }
            gap_prelim_edit_block_add(edit_script, GapAlignOpType::Sub, jump.dcp);
        } else if jump.dcp > jump.dcq {
            let gap_len = jump.dcp - jump.dcq;
            score += gap_open_score + gap_extend_score * gap_len;
            gap_prelim_edit_block_add(edit_script, GapAlignOpType::Ins, gap_len);
        } else {
            let gap_len = jump.dcq - jump.dcp;
            score += gap_open_score + gap_extend_score * gap_len;
            gap_prelim_edit_block_add(edit_script, GapAlignOpType::Del, gap_len);
        }

        cp += jump.dcp;
        cq += jump.dcq;

        if jump.ok == 0 && jump.lng != 0 {
            score += match_score * jump.lng;
            cp += jump.lng;
            cq += jump.lng;
            gap_prelim_edit_block_add(edit_script, GapAlignOpType::Sub, jump.lng);
            trace = trace.wrapping_shl(jump.lng as u32);
        }
    }

    if new_matches != 0 {
        gap_prelim_edit_block_add(edit_script, GapAlignOpType::Sub, new_matches);
    }

    (score, cp, cq)
}

/// NCBI: JumperExtendRightWithTraceback (jumper.c:1552).
pub fn jumper_extend_right_with_traceback(
    query: &[u8],
    subject: &[u8],
    match_score: i32,
    mismatch_score: i32,
    gap_open_score: i32,
    gap_extend_score: i32,
    max_mismatches: i32,
    window: i32,
    edit_script: &mut JumperPrelimEditBlock,
    num_identical: &mut i32,
    left_extension: bool,
    ungapped_ext_len: &mut i32,
) -> (i32, i32, i32) {
    if query.is_empty() || subject.is_empty() {
        return (-1, 0, 0);
    }

    let mut cp = 0;
    let mut cq = 0;
    let mut num_mismatches = 0;
    let mut new_matches = 0;
    let mut trace = 0u32;
    let trace_mask = jumper_trace_mask(max_mismatches);
    let mut is_ungapped = true;

    if left_extension {
        cp += 1;
        cq += 1;
    }

    while cp < query.len() as i32 && cq < subject.len() as i32 && num_mismatches < max_mismatches {
        if query[cp as usize] == subject[cq as usize] {
            cp += 1;
            cq += 1;
            new_matches += 1;
            continue;
        }

        let jump = jumper_find_right(query, subject, cp, cq, &JUMPER_DEFAULT);
        if new_matches != 0 {
            edit_script.edit_ops.push(new_matches);
            jumper_shift_trace(&mut trace, new_matches, window);
            *num_identical += new_matches;
            new_matches = 0;
        }

        if jump.dcp == jump.dcq {
            if trace & trace_mask != 0 {
                num_mismatches += jump.dcp;
                trace = trace.wrapping_shl(jump.dcp as u32);
                trace |= (1u32 << jump.dcp) - 1;
            } else {
                num_mismatches = jump.dcp;
                trace = (1u32 << jump.dcp) - 1;
            }
            for _ in 0..jump.dcp {
                edit_script.edit_ops.push(JUMPER_MISMATCH);
            }
        } else if jump.dcp > jump.dcq {
            for _ in 0..(jump.dcp - jump.dcq) {
                edit_script.edit_ops.push(JUMPER_INSERTION);
            }
        } else {
            for _ in 0..(jump.dcq - jump.dcp) {
                edit_script.edit_ops.push(JUMPER_DELETION);
            }
        }

        if is_ungapped && jump.dcp != jump.dcq {
            *ungapped_ext_len = cp - 1;
            is_ungapped = false;
        }

        cp += jump.dcp;
        cq += jump.dcq;

        if jump.ok == 0 && jump.lng != 0 {
            cp += jump.lng;
            cq += jump.lng;
            edit_script.edit_ops.push(jump.lng);
            trace = trace.wrapping_shl(jump.lng as u32);
            *num_identical += jump.lng;
            num_mismatches = 0;
        }
    }

    if new_matches != 0 {
        edit_script.edit_ops.push(new_matches);
        *num_identical += new_matches;
    }

    s_trim_extension(
        edit_script,
        -mismatch_score,
        &mut cp,
        &mut cq,
        num_identical,
        true,
    );

    if is_ungapped {
        *ungapped_ext_len = cp;
    }

    (
        s_compute_extension_score(
            edit_script,
            match_score,
            mismatch_score,
            gap_open_score,
            gap_extend_score,
        ),
        cp,
        cq,
    )
}

/// NCBI: JumperExtendLeft (jumper.c:2354).
pub fn jumper_extend_left(
    query: &[u8],
    subject: &[u8],
    query_offset: i32,
    subject_offset: i32,
    match_score: i32,
    mismatch_score: i32,
    gap_open_score: i32,
    gap_extend_score: i32,
    max_mismatches: i32,
    window: i32,
    edit_script: &mut GapPrelimEditBlock,
) -> (i32, i32, i32) {
    if query.is_empty()
        || subject.is_empty()
        || query_offset < 0
        || subject_offset < 0
        || query_offset as usize >= query.len()
        || subject_offset as usize >= subject.len()
    {
        return (-1, 0, 0);
    }

    let mut cp = query_offset;
    let mut cq = subject_offset;
    let mut score = 0;
    let mut num_mismatches = 0;
    let mut new_matches = 0;
    let mut trace = 0u32;
    let trace_mask = jumper_trace_mask(max_mismatches);

    while cp >= 0 && cq >= 0 && num_mismatches < max_mismatches {
        if query[cp as usize] == subject[cq as usize] {
            score += match_score;
            cp -= 1;
            cq -= 1;
            new_matches += 1;
            continue;
        }

        let jump = jumper_find_left(query, subject, cp, cq, &JUMPER_DEFAULT);
        if new_matches != 0 {
            gap_prelim_edit_block_add(edit_script, GapAlignOpType::Sub, new_matches);
            jumper_shift_trace(&mut trace, new_matches, window);
            new_matches = 0;
        }

        if jump.dcp == jump.dcq {
            score += mismatch_score * jump.dcp;
            if trace & trace_mask != 0 {
                num_mismatches += jump.dcp;
                trace = trace.wrapping_shl(jump.dcp as u32);
                trace |= (1u32 << jump.dcp) - 1;
            } else {
                num_mismatches = jump.dcp;
                trace = (1u32 << jump.dcp) - 1;
            }
            gap_prelim_edit_block_add(edit_script, GapAlignOpType::Sub, jump.dcp);
        } else if jump.dcp > jump.dcq {
            let gap_len = jump.dcp - jump.dcq;
            score += gap_open_score + gap_extend_score * gap_len;
            gap_prelim_edit_block_add(edit_script, GapAlignOpType::Ins, gap_len);
        } else {
            let gap_len = jump.dcq - jump.dcp;
            score += gap_open_score + gap_extend_score * gap_len;
            gap_prelim_edit_block_add(edit_script, GapAlignOpType::Del, gap_len);
        }

        cp -= jump.dcp;
        cq -= jump.dcq;

        if jump.ok == 0 && jump.lng != 0 {
            score += match_score * jump.lng;
            cp -= jump.lng;
            cq -= jump.lng;
            gap_prelim_edit_block_add(edit_script, GapAlignOpType::Sub, jump.lng);
            trace = trace.wrapping_shl(jump.lng as u32);
        }
    }

    if new_matches != 0 {
        gap_prelim_edit_block_add(edit_script, GapAlignOpType::Sub, new_matches);
    }

    (
        score,
        query_offset.saturating_sub(cp),
        subject_offset.saturating_sub(cq),
    )
}

/// NCBI: JumperExtendRightCompressed (jumper.c:734).
pub fn jumper_extend_right_compressed(
    query: &[u8],
    subject: &[u8],
    subject_length: i32,
    match_score: i32,
    mismatch_score: i32,
    max_mismatches: i32,
    window: i32,
    num_identical: &mut i32,
    ungapped_ext_len: &mut i32,
) -> (i32, i32, i32) {
    if query.is_empty() || subject.is_empty() || subject_length <= 0 {
        return (-1, 0, 0);
    }

    let mut cp = 1;
    let mut cq = 1;
    let mut cpstop = 0;
    let mut cqstop = 0;
    let mut num_mismatches = 0;
    let mut new_matches = 0;
    let mut trace = 0u32;
    let trace_mask = jumper_trace_mask(max_mismatches);
    let mut score = 0;
    let mut best_score = 0;
    let mut is_ungapped = true;

    while cp < query.len() as i32 && cq < subject_length && num_mismatches < max_mismatches {
        if query[cp as usize] == crate::encoding::ncbi2na_base_at(subject, cq as usize) {
            cp += 1;
            cq += 1;
            new_matches += 1;
            continue;
        }

        let jump =
            jumper_find_right_compressed(query, subject, subject_length, cp, cq, &JUMPER_DEFAULT);
        if new_matches != 0 {
            jumper_shift_trace(&mut trace, new_matches, window);
            *num_identical += new_matches;
            score += new_matches * match_score;
            new_matches = 0;
        }

        if jump.dcp == jump.dcq {
            score += mismatch_score * jump.dcp;
            if trace & trace_mask != 0 {
                num_mismatches += jump.dcp;
                trace = trace.wrapping_shl(jump.dcp as u32);
                trace |= (1u32 << jump.dcp) - 1;
            } else {
                num_mismatches = jump.dcp;
                trace = (1u32 << jump.dcp) - 1;
            }
        }

        if is_ungapped && jump.dcp != jump.dcq {
            *ungapped_ext_len = cp - 1;
            is_ungapped = false;
        }

        cp += jump.dcp;
        cq += jump.dcq;

        if jump.ok == 0 && jump.lng != 0 {
            cp += jump.lng;
            cq += jump.lng;
            trace = trace.wrapping_shl(jump.lng as u32);
            *num_identical += jump.lng;
            score += jump.lng * match_score;
        }

        if score >= best_score {
            cpstop = cp;
            cqstop = cq;
            best_score = score;
        }
    }

    if new_matches != 0 {
        *num_identical += new_matches;
        score += new_matches * match_score;
        if score >= best_score {
            cpstop = cp;
            cqstop = cq;
            best_score = score;
        }
    }

    if is_ungapped {
        *ungapped_ext_len = cpstop;
    }

    (best_score, cpstop, cqstop)
}

/// NCBI: JumperExtendLeftCompressed (jumper.c:1749).
pub fn jumper_extend_left_compressed(
    query: &[u8],
    subject: &[u8],
    query_offset: i32,
    subject_offset: i32,
    match_score: i32,
    mismatch_score: i32,
    max_mismatches: i32,
    window: i32,
    num_identical: &mut i32,
) -> (i32, i32, i32) {
    if query.is_empty()
        || subject.is_empty()
        || query_offset < 0
        || subject_offset < 0
        || query_offset as usize >= query.len()
    {
        return (-1, 0, 0);
    }

    let mut cp = query_offset;
    let mut cq = subject_offset;
    let mut cpstop = query_offset;
    let mut cqstop = subject_offset;
    let mut num_mismatches = 0;
    let mut new_matches = 0;
    let mut trace = 0u32;
    let trace_mask = jumper_trace_mask(max_mismatches);
    let mut score = 0;
    let mut best_score = 0;

    while cp >= 0 && cq >= 0 && num_mismatches < max_mismatches {
        if query[cp as usize] == crate::encoding::ncbi2na_base_at(subject, cq as usize) {
            cp -= 1;
            cq -= 1;
            new_matches += 1;
            continue;
        }

        let jump = jumper_find_left_compressed(query, subject, cp, cq, &JUMPER_DEFAULT);
        if new_matches != 0 {
            jumper_shift_trace(&mut trace, new_matches, window);
            *num_identical += new_matches;
            score = new_matches * match_score;
            new_matches = 0;
        }

        if jump.dcp == jump.dcq {
            score += mismatch_score * jump.dcp;
            if trace & trace_mask != 0 {
                num_mismatches += jump.dcp;
                trace = trace.wrapping_shl(jump.dcp as u32);
                trace |= (1u32 << jump.dcp) - 1;
            } else {
                num_mismatches = jump.dcp;
                trace = (1u32 << jump.dcp) - 1;
            }
        }

        cp -= jump.dcp;
        cq -= jump.dcq;

        if jump.ok == 0 && jump.lng != 0 {
            cp -= jump.lng;
            cq -= jump.lng;
            trace = trace.wrapping_shl(jump.lng as u32);
            *num_identical += jump.lng;
            score += jump.lng * match_score;
        }

        if score >= best_score {
            cpstop = cp;
            cqstop = cq;
            best_score = score;
        }
    }

    if new_matches != 0 {
        *num_identical += new_matches;
        score += new_matches * match_score;
        if score >= best_score {
            cpstop = cp;
            cqstop = cq;
            best_score = score;
        }
    }

    (
        best_score,
        query_offset.saturating_sub(cpstop),
        subject_offset.saturating_sub(cqstop),
    )
}

/// NCBI: JumperExtendRightCompressedWithTraceback (jumper.c:915).
pub fn jumper_extend_right_compressed_with_traceback(
    query: &[u8],
    subject: &[u8],
    subject_length: i32,
    match_score: i32,
    mismatch_score: i32,
    gap_open_score: i32,
    gap_extend_score: i32,
    max_mismatches: i32,
    window: i32,
    edit_script: &mut JumperPrelimEditBlock,
    num_identical: &mut i32,
    left_extension: bool,
    ungapped_ext_len: &mut i32,
) -> (i32, i32, i32) {
    if query.is_empty() || subject.is_empty() || subject_length <= 0 {
        return (-1, 0, 0);
    }

    let mut cp = 1;
    let mut cq = 1;
    let mut num_mismatches = 0;
    let mut new_matches = if left_extension { 0 } else { 1 };
    let mut trace = 0u32;
    let trace_mask = jumper_trace_mask(max_mismatches);
    let mut is_ungapped = true;

    while cp < query.len() as i32 && cq < subject_length && num_mismatches < max_mismatches {
        if query[cp as usize] == crate::encoding::ncbi2na_base_at(subject, cq as usize) {
            cp += 1;
            cq += 1;
            new_matches += 1;
            continue;
        }

        let jump =
            jumper_find_right_compressed(query, subject, subject_length, cp, cq, &JUMPER_DEFAULT);
        if new_matches != 0 {
            edit_script.edit_ops.push(new_matches);
            jumper_shift_trace(&mut trace, new_matches, window);
            *num_identical += new_matches;
            new_matches = 0;
        }

        if jump.dcp == jump.dcq {
            if trace & trace_mask != 0 {
                num_mismatches += jump.dcp;
                trace = trace.wrapping_shl(jump.dcp as u32);
                trace |= (1u32 << jump.dcp) - 1;
            } else {
                num_mismatches = jump.dcp;
                trace = (1u32 << jump.dcp) - 1;
            }
            for _ in 0..jump.dcp {
                edit_script.edit_ops.push(JUMPER_MISMATCH);
            }
        } else if jump.dcp > jump.dcq {
            for _ in 0..(jump.dcp - jump.dcq) {
                edit_script.edit_ops.push(JUMPER_INSERTION);
            }
        } else {
            for _ in 0..(jump.dcq - jump.dcp) {
                edit_script.edit_ops.push(JUMPER_DELETION);
            }
        }

        if is_ungapped && jump.dcp != jump.dcq {
            *ungapped_ext_len = cp - 1;
            is_ungapped = false;
        }

        cp += jump.dcp;
        cq += jump.dcq;

        if jump.ok == 0 && jump.lng != 0 {
            cp += jump.lng;
            cq += jump.lng;
            edit_script.edit_ops.push(jump.lng);
            trace = trace.wrapping_shl(jump.lng as u32);
            *num_identical += jump.lng;
            num_mismatches = 0;
        }
    }

    if new_matches != 0 {
        edit_script.edit_ops.push(new_matches);
        *num_identical += new_matches;
    }

    s_trim_extension(
        edit_script,
        -mismatch_score,
        &mut cp,
        &mut cq,
        num_identical,
        true,
    );

    if is_ungapped {
        *ungapped_ext_len = cp;
    }

    (
        s_compute_extension_score(
            edit_script,
            match_score,
            mismatch_score,
            gap_open_score,
            gap_extend_score,
        ),
        cp,
        cq,
    )
}

/// NCBI: JumperExtendLeftCompressedWithTraceback (jumper.c:1917).
pub fn jumper_extend_left_compressed_with_traceback(
    query: &[u8],
    subject: &[u8],
    query_offset: i32,
    subject_offset: i32,
    match_score: i32,
    mismatch_score: i32,
    gap_open_score: i32,
    gap_extend_score: i32,
    max_mismatches: i32,
    window: i32,
    edit_script: &mut JumperPrelimEditBlock,
    num_identical: &mut i32,
) -> (i32, i32, i32) {
    if query.is_empty()
        || subject.is_empty()
        || query_offset < 0
        || subject_offset < 0
        || query_offset as usize >= query.len()
    {
        return (-1, 0, 0);
    }

    let mut cp = query_offset;
    let mut cq = subject_offset;
    let mut num_mismatches = 0;
    let mut new_matches = 0;
    let mut trace = 0u32;
    let trace_mask = jumper_trace_mask(max_mismatches);

    while cp >= 0 && cq >= 0 && num_mismatches < max_mismatches {
        if query[cp as usize] == crate::encoding::ncbi2na_base_at(subject, cq as usize) {
            cp -= 1;
            cq -= 1;
            new_matches += 1;
            continue;
        }

        let jump = jumper_find_left_compressed(query, subject, cp, cq, &JUMPER_DEFAULT);
        if new_matches != 0 {
            edit_script.edit_ops.push(new_matches);
            jumper_shift_trace(&mut trace, new_matches, window);
            *num_identical += new_matches;
            new_matches = 0;
        }

        if jump.dcp == jump.dcq {
            if trace & trace_mask != 0 {
                num_mismatches += jump.dcp;
                trace = trace.wrapping_shl(jump.dcp as u32);
                trace |= (1u32 << jump.dcp) - 1;
            } else {
                num_mismatches = jump.dcp;
                trace = (1u32 << jump.dcp) - 1;
            }
            for _ in 0..jump.dcp {
                edit_script.edit_ops.push(JUMPER_MISMATCH);
            }
        } else if jump.dcp > jump.dcq {
            for _ in 0..(jump.dcp - jump.dcq) {
                edit_script.edit_ops.push(JUMPER_INSERTION);
            }
        } else {
            for _ in 0..(jump.dcq - jump.dcp) {
                edit_script.edit_ops.push(JUMPER_DELETION);
            }
        }

        cp -= jump.dcp;
        cq -= jump.dcq;

        if jump.ok == 0 && jump.lng != 0 {
            cp -= jump.lng;
            cq -= jump.lng;
            edit_script.edit_ops.push(jump.lng);
            trace = trace.wrapping_shl(jump.lng as u32);
            *num_identical += jump.lng;
            num_mismatches = 0;
        }
    }

    if new_matches != 0 {
        edit_script.edit_ops.push(new_matches);
        *num_identical += new_matches;
    }

    s_trim_extension(
        edit_script,
        -mismatch_score,
        &mut cp,
        &mut cq,
        num_identical,
        false,
    );

    (
        s_compute_extension_score(
            edit_script,
            match_score,
            mismatch_score,
            gap_open_score,
            gap_extend_score,
        ),
        query_offset.saturating_sub(cp),
        subject_offset.saturating_sub(cq),
    )
}

/// NCBI: JumperExtendRightCompressedWithTracebackOptimal (jumper.c:1124).
pub fn jumper_extend_right_compressed_with_traceback_optimal(
    query: &[u8],
    subject: &[u8],
    subject_length: i32,
    match_score: i32,
    mismatch_score: i32,
    gap_open_score: i32,
    gap_extend_score: i32,
    max_mismatches: i32,
    window: i32,
    x_drop: i32,
    edit_script: &mut JumperPrelimEditBlock,
    best_num_identical: &mut i32,
    left_extension: bool,
    ungapped_ext_len: &mut i32,
) -> (i32, i32, i32) {
    if query.is_empty() || subject.is_empty() || subject_length <= 0 {
        return (-1, 0, 0);
    }

    let mut cp = 1;
    let mut cq = 1;
    let mut cpstop = 0;
    let mut cqstop = 0;
    let mut num_mismatches = 0;
    let mut new_matches = if left_extension { 0 } else { 1 };
    let mut trace = 0u32;
    let trace_mask = jumper_trace_mask(max_mismatches);
    let mut is_ungapped = true;
    let mut score = 0;
    let mut best_score = 0;
    let mut num_ops = 0usize;
    let mut num_identical = *best_num_identical;
    let mut last_gap_open = 0;

    while cp < query.len() as i32 && cq < subject_length && num_mismatches < max_mismatches {
        if query[cp as usize] == crate::encoding::ncbi2na_base_at(subject, cq as usize) {
            cp += 1;
            cq += 1;
            new_matches += 1;
            continue;
        }

        let jump =
            jumper_find_right_compressed(query, subject, subject_length, cp, cq, &JUMPER_DEFAULT);
        if new_matches != 0 {
            edit_script.edit_ops.push(new_matches);
            jumper_shift_trace(&mut trace, new_matches, window);
            num_identical += new_matches;
            score += new_matches * match_score;
            new_matches = 0;
            last_gap_open = 0;
        }

        if score >= best_score {
            cpstop = cp;
            cqstop = cq;
            num_ops = edit_script.edit_ops.len();
            best_score = score;
            *best_num_identical = num_identical;
        }
        if best_score - score > x_drop {
            break;
        }

        if jump.dcp == jump.dcq {
            score += jump.dcp * mismatch_score;
            if trace & trace_mask != 0 {
                num_mismatches += jump.dcp;
                trace = trace.wrapping_shl(jump.dcp as u32);
                trace |= (1u32 << jump.dcp) - 1;
            } else {
                num_mismatches = jump.dcp;
                trace = (1u32 << jump.dcp) - 1;
            }
            for _ in 0..jump.dcp {
                edit_script.edit_ops.push(JUMPER_MISMATCH);
            }
        } else if jump.dcp > jump.dcq {
            for _ in 0..(jump.dcp - jump.dcq) {
                edit_script.edit_ops.push(JUMPER_INSERTION);
                score += gap_extend_score;
            }
            if last_gap_open != JUMPER_INSERTION {
                score += gap_open_score;
            }
            last_gap_open = JUMPER_INSERTION;
        } else {
            for _ in 0..(jump.dcq - jump.dcp) {
                edit_script.edit_ops.push(JUMPER_DELETION);
                score += gap_extend_score;
            }
            if last_gap_open != JUMPER_DELETION {
                score += gap_open_score;
            }
            last_gap_open = JUMPER_DELETION;
        }

        if is_ungapped && jump.dcp != jump.dcq {
            *ungapped_ext_len = cp - 1;
            is_ungapped = false;
        }

        cp += jump.dcp;
        cq += jump.dcq;

        if jump.ok == 0 && jump.lng != 0 {
            cp += jump.lng;
            cq += jump.lng;
            edit_script.edit_ops.push(jump.lng);
            trace = trace.wrapping_shl(jump.lng as u32);
            num_identical += jump.lng;
            score += jump.lng * match_score;
            last_gap_open = 0;
        }

        if score >= best_score {
            cpstop = cp;
            cqstop = cq;
            num_ops = edit_script.edit_ops.len();
            best_score = score;
            *best_num_identical = num_identical;
        }
    }

    if new_matches != 0 {
        edit_script.edit_ops.push(new_matches);
        num_identical += new_matches;
        score += new_matches;
    }

    if score >= best_score {
        cpstop = cp;
        cqstop = cq;
        num_ops = edit_script.edit_ops.len();
        best_score = score;
        *best_num_identical = num_identical;
    }

    edit_script.edit_ops.truncate(num_ops);
    if is_ungapped {
        *ungapped_ext_len = cpstop;
    }

    (best_score, cpstop, cqstop)
}

/// NCBI: JumperExtendLeftCompressedWithTracebackOptimal (jumper.c:2110).
pub fn jumper_extend_left_compressed_with_traceback_optimal(
    query: &[u8],
    subject: &[u8],
    query_offset: i32,
    subject_offset: i32,
    match_score: i32,
    mismatch_score: i32,
    gap_open_score: i32,
    gap_extend_score: i32,
    max_mismatches: i32,
    window: i32,
    x_drop: i32,
    edit_script: &mut JumperPrelimEditBlock,
    best_num_identical: &mut i32,
) -> (i32, i32, i32) {
    if query.is_empty()
        || subject.is_empty()
        || query_offset < 0
        || subject_offset < 0
        || query_offset as usize >= query.len()
    {
        return (-1, 0, 0);
    }

    let mut cp = query_offset;
    let mut cq = subject_offset;
    let mut cpstop = query_offset;
    let mut cqstop = subject_offset;
    let mut num_mismatches = 0;
    let mut new_matches = 0;
    let mut trace = 0u32;
    let trace_mask = jumper_trace_mask(max_mismatches);
    let mut score = 0;
    let mut best_score = 0;
    let mut num_ops = 0usize;
    let mut num_identical = *best_num_identical;
    let mut last_gap_open = 0;

    while cp >= 0 && cq >= 0 && num_mismatches < max_mismatches {
        if query[cp as usize] == crate::encoding::ncbi2na_base_at(subject, cq as usize) {
            cp -= 1;
            cq -= 1;
            new_matches += 1;
            continue;
        }

        let jump = jumper_find_left_compressed(query, subject, cp, cq, &JUMPER_DEFAULT);
        if new_matches != 0 {
            edit_script.edit_ops.push(new_matches);
            jumper_shift_trace(&mut trace, new_matches, window);
            num_identical += new_matches;
            score += new_matches * match_score;
            new_matches = 0;
            last_gap_open = 0;
        }

        if score >= best_score {
            cpstop = cp;
            cqstop = cq;
            best_score = score;
            num_ops = edit_script.edit_ops.len();
            *best_num_identical = num_identical;
        }
        if best_score - score > x_drop {
            break;
        }

        if jump.dcp == jump.dcq {
            score += jump.dcp * mismatch_score;
            if trace & trace_mask != 0 {
                num_mismatches += jump.dcp;
                trace = trace.wrapping_shl(jump.dcp as u32);
                trace |= (1u32 << jump.dcp) - 1;
            } else {
                num_mismatches = jump.dcp;
                trace = (1u32 << jump.dcp) - 1;
            }
            for _ in 0..jump.dcp {
                edit_script.edit_ops.push(JUMPER_MISMATCH);
            }
        } else if jump.dcp > jump.dcq {
            for _ in 0..(jump.dcp - jump.dcq) {
                edit_script.edit_ops.push(JUMPER_INSERTION);
                score += gap_extend_score;
            }
            if last_gap_open != JUMPER_INSERTION {
                score += gap_open_score;
            }
            last_gap_open = JUMPER_INSERTION;
        } else {
            for _ in 0..(jump.dcq - jump.dcp) {
                edit_script.edit_ops.push(JUMPER_DELETION);
                score += gap_extend_score;
            }
            if last_gap_open != JUMPER_DELETION {
                score += gap_open_score;
            }
            last_gap_open = JUMPER_DELETION;
        }

        cp -= jump.dcp;
        cq -= jump.dcq;

        if jump.ok == 0 && jump.lng != 0 {
            cp -= jump.lng;
            cq -= jump.lng;
            edit_script.edit_ops.push(jump.lng);
            trace = trace.wrapping_shl(jump.lng as u32);
            num_identical += jump.lng;
            score += jump.lng * match_score;
            last_gap_open = 0;
        }

        if score >= best_score {
            cpstop = cp;
            cqstop = cq;
            best_score = score;
            num_ops = edit_script.edit_ops.len();
            *best_num_identical = num_identical;
        }
    }

    if new_matches != 0 {
        edit_script.edit_ops.push(new_matches);
        num_identical += new_matches;
        score += new_matches * match_score;
    }

    if score >= best_score {
        cpstop = cp;
        cqstop = cq;
        best_score = score;
        num_ops = edit_script.edit_ops.len();
        *best_num_identical = num_identical;
    }

    edit_script.edit_ops.truncate(num_ops);
    (
        best_score,
        query_offset.saturating_sub(cpstop),
        subject_offset.saturating_sub(cqstop),
    )
}

/// blast-rs: Port-shaped equivalent of NCBI `JumperGappedAlignmentCompressedWithTraceback`
/// (jumper.c:2512); not a direct NCBI C port.
#[allow(clippy::too_many_arguments)]
pub fn jumper_gapped_alignment_compressed_with_traceback(
    query: &[u8],
    subject: &[u8],
    query_length: i32,
    subject_length: i32,
    query_start: i32,
    subject_start: i32,
    jumper: &mut JumperGapAlign,
    score_params: &crate::parameters::ScoringParameters,
    align_params: JumperAlignParams,
    num_identical: &mut i32,
    right_ungapped_ext_len: &mut i32,
) -> i32 {
    if query_start < 0
        || subject_start < 0
        || query_length <= 0
        || subject_length <= 0
        || query_length as usize > query.len()
        || subject_length as usize > subject.len().saturating_mul(4)
        || query_start >= query_length
        || subject_start >= subject_length
    {
        return -1;
    }

    *num_identical = 0;
    *right_ungapped_ext_len = 0;
    jumper
        .left_prelim_block
        .get_or_insert_with(JumperPrelimEditBlock::default)
        .edit_ops
        .clear();
    jumper
        .right_prelim_block
        .get_or_insert_with(JumperPrelimEditBlock::default)
        .edit_ops
        .clear();

    let mut score_left = 0;
    let mut score_right = 0;
    let mut left_ext_done = false;
    let offset_adjustment = 4 - (subject_start % 4);
    let q_length = query_start + offset_adjustment;
    let s_length = subject_start + offset_adjustment;

    let mut query_align_start = query_start;
    let mut subject_align_start = subject_start;
    let mut query_align_stop = query_start;
    let mut subject_align_stop = subject_start;

    if query_start > 0 && subject_start > 0 && q_length < query_length && s_length < subject_length
    {
        let rev_prelim_block = jumper
            .left_prelim_block
            .as_mut()
            .expect("left block initialized above");
        let (score, q_ext_len, s_ext_len) = jumper_extend_left_compressed_with_traceback_optimal(
            query,
            subject,
            q_length,
            s_length,
            score_params.reward,
            score_params.penalty,
            -score_params.gap_open,
            -score_params.gap_extend,
            align_params.max_mismatches,
            align_params.mismatch_window,
            align_params.gap_x_dropoff,
            rev_prelim_block,
            num_identical,
        );
        score_left = score;
        query_align_start = q_length.saturating_sub(q_ext_len).saturating_add(1);
        subject_align_start = s_length.saturating_sub(s_ext_len).saturating_add(1);
        left_ext_done = true;
    }

    if query_start < query_length - 1 && subject_start < subject_length - 1 {
        let subject_byte = ((s_length + 3) / 4).max(0) as usize;
        let subject_right = subject.get(subject_byte..).unwrap_or(&[]);
        let query_right = query.get(q_length.max(0) as usize..).unwrap_or(&[]);
        let fwd_prelim_block = jumper
            .right_prelim_block
            .as_mut()
            .expect("right block initialized above");
        let (score, q_ext_len, s_ext_len) = jumper_extend_right_compressed_with_traceback_optimal(
            query_right,
            subject_right,
            subject_length - s_length,
            score_params.reward,
            score_params.penalty,
            -score_params.gap_open,
            -score_params.gap_extend,
            align_params.max_mismatches,
            align_params.mismatch_window,
            align_params.gap_x_dropoff,
            fwd_prelim_block,
            num_identical,
            left_ext_done,
            right_ungapped_ext_len,
        );
        score_right = score;
        query_align_stop = q_length + q_ext_len;
        subject_align_stop = s_length + s_ext_len;
    }

    let mut score = score_left + score_right;
    if offset_adjustment != 0 && !left_ext_done {
        if let Some(block) = jumper.left_prelim_block.as_mut() {
            block.edit_ops.push(offset_adjustment);
        }
        *num_identical = num_identical.saturating_add(offset_adjustment);
        score = score.saturating_add(offset_adjustment.saturating_mul(score_params.reward));
    }
    if offset_adjustment != 0 && *right_ungapped_ext_len != 0 {
        *right_ungapped_ext_len = right_ungapped_ext_len.saturating_add(offset_adjustment);
    }

    const K_BASE_N: u8 = 14;
    let start = query_align_start.max(0) as usize;
    let stop = query_align_stop.max(query_align_start).min(query_length) as usize;
    for &base in query.get(start..stop).unwrap_or(&[]) {
        if base == K_BASE_N {
            score -= score_params.penalty;
        }
    }

    jumper.table.clear();
    jumper.table.extend_from_slice(&[
        score as u32,
        query_align_start as u32,
        query_align_stop as u32,
        subject_align_start as u32,
        subject_align_stop as u32,
    ]);
    0
}

/// blast-rs: Selects lookup word lengths for Jumper scanning; not a direct NCBI C port.
fn jumper_lookup_lengths(lookup_wrap: &crate::lookup::LookupTableWrap) -> Option<(i32, i32)> {
    jumper_scan_lengths(lookup_wrap)
}

/// blast-rs: Lookup word-length adapter for represented scan paths; not a direct NCBI C port.
fn jumper_scan_lengths(lookup_wrap: &crate::lookup::LookupTableWrap) -> Option<(i32, i32)> {
    match lookup_wrap {
        crate::lookup::LookupTableWrap::Na(table) => {
            Some((table.word_length, table.lut_word_length))
        }
        crate::lookup::LookupTableWrap::NaHash(table) => {
            Some((table.word_length, table.lut_word_length))
        }
        crate::lookup::LookupTableWrap::SmallNa(table) => Some((
            table.word_length,
            crate::lookup::small_na_lut_word_length(table),
        )),
        crate::lookup::LookupTableWrap::Megablast(table) => {
            if table.discontiguous {
                let template_length = table.template_length.max(table.word_length);
                Some((template_length, template_length))
            } else {
                Some((table.word_length, table.lut_word_length))
            }
        }
        crate::lookup::LookupTableWrap::Aa(_) | crate::lookup::LookupTableWrap::Rps(_) => None,
    }
}

/// blast-rs: Lookup word-length adapter for Jumper extension paths; not a direct NCBI C port.
fn jumper_extend_lengths(lookup_wrap: &crate::lookup::LookupTableWrap) -> Option<(i32, i32)> {
    match lookup_wrap {
        crate::lookup::LookupTableWrap::Na(table) => {
            Some((table.word_length, table.lut_word_length))
        }
        crate::lookup::LookupTableWrap::NaHash(table) => {
            Some((table.word_length, table.lut_word_length))
        }
        crate::lookup::LookupTableWrap::SmallNa(table) => Some((
            table.word_length,
            crate::lookup::small_na_lut_word_length(table),
        )),
        crate::lookup::LookupTableWrap::Megablast(table) => {
            Some((table.word_length, table.lut_word_length))
        }
        crate::lookup::LookupTableWrap::Aa(_) | crate::lookup::LookupTableWrap::Rps(_) => None,
    }
}

/// blast-rs: Decodes Rust scratch-table alignment outputs; not a direct NCBI C port.
fn jumper_alignment_outputs(jumper: &JumperGapAlign) -> Option<(i32, i32, i32, i32, i32)> {
    let data = &jumper.table;
    if data.len() < 5 {
        return None;
    }
    Some((
        data[0] as i32,
        data[1] as i32,
        data[2] as i32,
        data[3] as i32,
        data[4] as i32,
    ))
}

/// blast-rs: Packed-subject word extraction for lookup scanning; not a direct NCBI C port.
fn packed_subject_word(
    subject: &[u8],
    subject_length: i32,
    s_off: i32,
    word_size: i32,
) -> Option<i64> {
    if s_off < 0 || word_size <= 0 {
        return None;
    }
    let end = s_off.checked_add(word_size)?;
    if end > subject_length {
        return None;
    }
    if end as usize > subject.len().saturating_mul(4) {
        return None;
    }
    let mut word = 0i64;
    for i in 0..word_size {
        word = (word << 2) | crate::encoding::ncbi2na_base_at(subject, (s_off + i) as usize) as i64;
    }
    Some(word)
}

/// NCBI: s_DetermineScanningOffsets (masksubj.inl:43).
///
/// `range` stores the current subject range index, scan start, and inclusive
/// scan end. The function advances past masked/unusable ranges until scanning
/// can proceed or all ranges are exhausted.
fn s_determine_scanning_offsets(
    subject_ranges: &[crate::util::SSeqRange],
    word_length: i32,
    lut_word_length: i32,
    range: &mut [i32; 3],
) -> bool {
    if subject_ranges.is_empty() {
        return false;
    }
    while range[1] > range[2] {
        range[0] = range[0].saturating_add(1);
        if range[0] >= subject_ranges.len() as i32 {
            return false;
        }
        let subject_range = subject_ranges[range[0] as usize];
        range[1] = subject_range
            .left
            .saturating_add(word_length)
            .saturating_sub(lut_word_length);
        range[2] = subject_range.right.saturating_sub(lut_word_length);
    }
    true
}

struct JumperWordHitBatch {
    offset_pairs: Vec<crate::lookup::OffsetPair>,
    s_range: u32,
}

/// blast-rs: Collects represented lookup hits for one subject scan range; not a direct NCBI C port.
fn jumper_collect_word_hits_for_range(
    lookup_wrap: &crate::lookup::LookupTableWrap,
    subject: &[u8],
    subject_length: i32,
    lut_word_length: i32,
    start: i32,
    end: i32,
) -> Vec<crate::lookup::OffsetPair> {
    let mut offset_pairs = Vec::new();
    if start > end {
        return offset_pairs;
    }

    match lookup_wrap {
        crate::lookup::LookupTableWrap::Na(lookup) => {
            let lut_word = lut_word_length.max(1);
            let scan_step = lookup.scan_step.max(1);
            let mut s_off = start.max(0);
            let capped_end = end.min(subject_length - lut_word);
            while s_off <= capped_end {
                if let Some(index) = packed_subject_word(subject, subject_length, s_off, lut_word) {
                    crate::lookup::s_blast_lookup_retrieve(
                        lookup,
                        index as i32,
                        &mut offset_pairs,
                        s_off,
                    );
                }
                s_off = s_off.saturating_add(scan_step);
            }
        }
        crate::lookup::LookupTableWrap::NaHash(lookup) => {
            if let Some(mut hits) = crate::lookup::s_blast_na_hash_scan_subject_any(
                lookup,
                subject,
                subject_length,
                start.max(0),
                end.min(subject_length - lut_word_length),
            ) {
                offset_pairs.append(&mut hits);
            }
        }
        crate::lookup::LookupTableWrap::Megablast(lookup) => {
            if lookup.discontiguous {
                offset_pairs.extend(crate::lookup::s_mb_disc_word_scan_subject(
                    lookup,
                    subject,
                    subject_length.max(0) as usize,
                    start.max(0) as usize,
                    end.max(0) as usize,
                ));
                return offset_pairs;
            }
            let lut_word = lut_word_length.max(1);
            let scan_step = lookup.scan_step.max(1);
            let mut s_off = start.max(0);
            let capped_end = end.min(subject_length - lut_word);
            while s_off <= capped_end {
                if let Some(index) = packed_subject_word(subject, subject_length, s_off, lut_word) {
                    crate::lookup::s_blast_mb_lookup_retrieve(
                        lookup,
                        index,
                        &mut offset_pairs,
                        s_off,
                    );
                }
                s_off = s_off.saturating_add(scan_step);
            }
        }
        crate::lookup::LookupTableWrap::SmallNa(lookup) => {
            let lut_word = lut_word_length.max(1);
            let scan_step = lookup.scan_step.max(1);
            let mut s_off = start.max(0);
            let capped_end = end.min(subject_length - lut_word);
            while s_off <= capped_end {
                if let Some(index) = packed_subject_word(subject, subject_length, s_off, lut_word) {
                    let Ok(index) = i32::try_from(index) else {
                        s_off = s_off.saturating_add(scan_step);
                        continue;
                    };
                    crate::lookup::s_blast_small_na_lookup_retrieve(
                        lookup,
                        index,
                        &mut offset_pairs,
                        s_off,
                    );
                }
                s_off = s_off.saturating_add(scan_step);
            }
        }
        crate::lookup::LookupTableWrap::Aa(_) | crate::lookup::LookupTableWrap::Rps(_) => {}
    }
    offset_pairs
}

/// blast-rs: Batches represented lookup hits by subject scan range; not a direct NCBI C port.
fn jumper_collect_word_hit_batches_in_ranges(
    lookup_wrap: &crate::lookup::LookupTableWrap,
    subject: &[u8],
    subject_length: i32,
    word_length: i32,
    lut_word_length: i32,
    subject_ranges: Option<&[crate::util::SSeqRange]>,
    scan_all_ranges: bool,
) -> Vec<JumperWordHitBatch> {
    if subject_length <= 0 || word_length <= 0 || lut_word_length <= 0 {
        return Vec::new();
    }
    if subject_length as usize > subject.len().saturating_mul(4) {
        return Vec::new();
    }
    let scan_ranges: Vec<(i32, i32)> = if let Some(ranges) = subject_ranges {
        if ranges.is_empty() {
            return Vec::new();
        }
        let mut scan_ranges = Vec::new();
        let mut scan_range = [
            0,
            ranges[0]
                .left
                .saturating_add(word_length)
                .saturating_sub(lut_word_length),
            ranges[0].right.saturating_sub(lut_word_length),
        ];
        while s_determine_scanning_offsets(ranges, word_length, lut_word_length, &mut scan_range) {
            scan_ranges.push((scan_range[1], scan_range[2]));
            if !scan_all_ranges {
                break;
            }
            scan_range[1] = scan_range[2].saturating_add(1);
        }
        scan_ranges
    } else {
        vec![(0, subject_length - lut_word_length)]
    };

    scan_ranges
        .into_iter()
        .map(|(start, end)| {
            let capped_end = end.min(subject_length.saturating_sub(lut_word_length));
            JumperWordHitBatch {
                offset_pairs: jumper_collect_word_hits_for_range(
                    lookup_wrap,
                    subject,
                    subject_length,
                    lut_word_length,
                    start,
                    end,
                ),
                s_range: capped_end.saturating_add(lut_word_length).max(0) as u32,
            }
        })
        .collect()
}

#[cfg(test)]
/// blast-rs: Test-only flattened view of Jumper word-hit batches; not a direct NCBI C port.
fn jumper_collect_word_hits_in_ranges(
    lookup_wrap: &crate::lookup::LookupTableWrap,
    subject: &[u8],
    subject_length: i32,
    word_length: i32,
    lut_word_length: i32,
    subject_ranges: Option<&[crate::util::SSeqRange]>,
    scan_all_ranges: bool,
) -> Vec<crate::lookup::OffsetPair> {
    jumper_collect_word_hit_batches_in_ranges(
        lookup_wrap,
        subject,
        subject_length,
        word_length,
        lut_word_length,
        subject_ranges,
        scan_all_ranges,
    )
    .into_iter()
    .flat_map(|batch| batch.offset_pairs)
    .collect()
}

/// blast-rs: Port-shaped equivalent of NCBI `BlastNaExtendJumper` (jumper.c:3253); not a direct NCBI C port.
///
/// This is the mapper/short-read driver around the already translated Jumper
/// extension primitives: it sorts word hits by diagonal/query/subject, performs
/// C's exact left/right pre-extension to the full word length, skips already
/// explored diagonal spans, runs compressed Jumper traceback, applies
/// `JumperGoodAlign`, and saves accepted HSPs. The optional `SubjectIndex`
/// enables the environment-controlled small-word rescue branch.
#[allow(clippy::too_many_arguments)]
pub fn blast_na_extend_jumper(
    offset_pairs: &mut [crate::lookup::OffsetPair],
    word_params: &crate::parameters::InitialWordParameters,
    score_params: &crate::parameters::ScoringParameters,
    hit_params: &crate::parameters::HitSavingParameters,
    lookup_wrap: &crate::lookup::LookupTableWrap,
    query: &[u8],
    subject: &[u8],
    subject_length: i32,
    subject_frame: i32,
    query_info: &crate::queryinfo::QueryInfo,
    align_params: JumperAlignParams,
    jumper: &mut JumperGapAlign,
    hsp_list: &mut crate::hspstream::HspList,
    s_range: u32,
    subject_index: Option<&SubjectIndex>,
) -> i32 {
    let _ = word_params;
    let Some((word_length, lut_word_length)) = jumper_extend_lengths(lookup_wrap) else {
        return -1;
    };
    if word_length <= 0 || lut_word_length <= 0 || word_length < lut_word_length {
        return -1;
    }
    if subject_length < 0 || subject_length as usize > subject.len().saturating_mul(4) {
        return 0;
    }
    let ext_to = word_length - lut_word_length;
    let mut hits_extended = 0;
    let mut skip_until = 0;
    let mut last_diag = 0i64;

    offset_pairs.sort_by(|a, b| {
        let a_diag = i64::from(a.subject_offset) - i64::from(a.query_offset);
        let b_diag = i64::from(b.subject_offset) - i64::from(b.query_offset);
        a_diag
            .cmp(&b_diag)
            .then_with(|| a.query_offset.cmp(&b.query_offset))
            .then_with(|| a.subject_offset.cmp(&b.subject_offset))
    });

    for pair in offset_pairs.iter().copied() {
        let mut s_offset = pair.subject_offset;
        let mut q_offset = pair.query_offset;
        if q_offset < 0
            || s_offset < 0
            || q_offset >= query.len() as i32
            || s_offset >= subject_length
        {
            continue;
        }
        let diag = i64::from(s_offset) - i64::from(q_offset);
        if diag == last_diag && q_offset < skip_until {
            continue;
        }

        let mut ext_left = 0;
        while ext_left < ext_to && ext_left < s_offset && q_offset - ext_left > 0 {
            let probe_s = s_offset - ext_left - 1;
            let probe_q = q_offset - ext_left - 1;
            if crate::encoding::ncbi2na_base_at(subject, probe_s as usize)
                != query.get(probe_q as usize).copied().unwrap_or(0xff)
            {
                break;
            }
            ext_left += 1;
        }

        if ext_left < ext_to {
            let mut ext_right = 0;
            let mut s_off = s_offset.saturating_add(lut_word_length);
            if s_off.saturating_add(ext_to).saturating_sub(ext_left) > s_range as i32 {
                continue;
            }
            let mut q_off = q_offset.saturating_add(lut_word_length);
            while ext_right < ext_to.saturating_sub(ext_left) {
                if crate::encoding::ncbi2na_base_at(subject, s_off as usize)
                    != query.get(q_off as usize).copied().unwrap_or(0xff)
                {
                    break;
                }
                s_off += 1;
                q_off += 1;
                ext_right += 1;
            }
            if ext_left + ext_right < ext_to {
                continue;
            }
        }

        q_offset -= ext_left;
        s_offset -= ext_left;

        let context = crate::queryinfo::bsearch_context_info(q_offset, query_info);
        let Some(context_info) = query_info.contexts.get(context.max(0) as usize) else {
            continue;
        };
        let query_start = context_info.query_offset;
        let local_q_offset = q_offset - query_start;
        if local_q_offset < 0 {
            continue;
        }
        let query_len = context_info.query_length;
        let query_slice = query
            .get(query_start.max(0) as usize..)
            .and_then(|slice| slice.get(..query_len.max(0) as usize))
            .unwrap_or(&[]);
        if query_slice.len() < query_len.max(0) as usize {
            continue;
        }
        let mut num_identical = 0;
        let mut right_ungapped_ext_len = 0;

        let status = jumper_gapped_alignment_compressed_with_traceback(
            query_slice,
            subject,
            query_len,
            subject_length,
            local_q_offset,
            s_offset,
            jumper,
            score_params,
            align_params,
            &mut num_identical,
            &mut right_ungapped_ext_len,
        );
        if status != 0 {
            continue;
        }
        hits_extended += 1;
        skip_until = local_q_offset
            .saturating_add(query_start)
            .saturating_add(right_ungapped_ext_len);
        last_diag = diag;

        let Some((
            score,
            query_align_start,
            query_align_stop,
            subject_align_start,
            subject_align_stop,
        )) = jumper_alignment_outputs(jumper)
        else {
            continue;
        };
        if !jumper_good_align(
            score,
            query_align_start,
            query_align_stop,
            subject_align_start,
            subject_align_stop,
            num_identical,
            &hit_params.options,
            query_len,
        ) {
            continue;
        }

        let Some((new_hsp, _map_info)) = s_create_hsp(
            query_slice,
            query_len,
            context,
            context_info.frame,
            subject,
            subject_length,
            subject_frame,
            jumper,
            query_align_start,
            query_align_stop,
            subject_align_start,
            subject_align_stop,
            score,
            score_params.penalty,
            hit_params.options.splice,
        ) else {
            return -1;
        };
        let saved_hsp = new_hsp.clone();
        let status = crate::hspstream::blast_hsp_list_save_hsp(hsp_list, new_hsp);
        if status != 0 {
            break;
        }

        if std::env::var_os("MAPPER_USE_SMALL_WORDS").is_some() {
            if let Some(s_index) = subject_index {
                let status = blast_na_extend_jumper_small_word_rescue(
                    &saved_hsp,
                    s_index,
                    context,
                    query_slice,
                    context_info.frame,
                    query_len,
                    subject,
                    subject_length,
                    subject_frame,
                    score_params,
                    hit_params.options.longest_intron,
                    hsp_list,
                );
                if status < 0 {
                    return -1;
                }
            }
        }
    }

    hits_extended
}

/// blast-rs: Port-shaped equivalent of NCBI `JumperNaWordFinder`
/// (na_ungapped.c:1995) for the represented contiguous lookup-table path; not a direct NCBI C port.
///
/// C obtains hits by dispatching through the selected `scansub_callback`; Rust
/// performs the same represented work directly for the current typed
/// `LookupTableWrap`, then delegates all hit extension and HSP saving to
/// [`blast_na_extend_jumper`].
#[allow(clippy::too_many_arguments)]
pub fn jumper_na_word_finder(
    subject: &[u8],
    subject_length: i32,
    subject_frame: i32,
    query: &[u8],
    query_info: &crate::queryinfo::QueryInfo,
    lookup_wrap: &crate::lookup::LookupTableWrap,
    word_params: &crate::parameters::InitialWordParameters,
    score_params: &crate::parameters::ScoringParameters,
    hit_params: &crate::parameters::HitSavingParameters,
    align_params: JumperAlignParams,
    jumper: &mut JumperGapAlign,
    hsp_list: &mut crate::hspstream::HspList,
) -> i16 {
    jumper_na_word_finder_with_word_hits(
        subject,
        subject_length,
        subject_frame,
        query,
        query_info,
        lookup_wrap,
        word_params,
        score_params,
        hit_params,
        align_params,
        jumper,
        hsp_list,
        None,
        None,
        None,
    )
}

/// Variant of [`jumper_na_word_finder_with_word_hits`] that takes subject
/// unmasked ranges, matching the `subject->mask_type != eNoSubjMasking` branch
/// of the C jumper word finder.
/// blast-rs: Masked-subject adapter; not a direct NCBI C port.
#[allow(clippy::too_many_arguments)]
pub fn jumper_na_word_finder_with_subject_ranges(
    subject: &[u8],
    subject_length: i32,
    subject_frame: i32,
    subject_ranges: &[crate::util::SSeqRange],
    query: &[u8],
    query_info: &crate::queryinfo::QueryInfo,
    lookup_wrap: &crate::lookup::LookupTableWrap,
    word_params: &crate::parameters::InitialWordParameters,
    score_params: &crate::parameters::ScoringParameters,
    hit_params: &crate::parameters::HitSavingParameters,
    align_params: JumperAlignParams,
    jumper: &mut JumperGapAlign,
    hsp_list: &mut crate::hspstream::HspList,
    word_hits: Option<&mut crate::lookup::MapperWordHits>,
    ungapped_stats: Option<&mut crate::diagnostics::UngappedStats>,
    gapped_stats: Option<&mut crate::diagnostics::GappedStats>,
) -> i16 {
    jumper_na_word_finder_impl(
        subject,
        subject_length,
        subject_frame,
        Some(subject_ranges),
        query,
        query_info,
        lookup_wrap,
        word_params,
        score_params,
        hit_params,
        align_params,
        jumper,
        hsp_list,
        word_hits,
        ungapped_stats,
        gapped_stats,
    )
}

/// blast-rs: Flushes one MapperWordHits bucket into Jumper extension; not a direct NCBI C port.
fn mapper_word_hits_flush(
    word_hits: &mut crate::lookup::MapperWordHits,
    index: usize,
    word_params: &crate::parameters::InitialWordParameters,
    score_params: &crate::parameters::ScoringParameters,
    hit_params: &crate::parameters::HitSavingParameters,
    lookup_wrap: &crate::lookup::LookupTableWrap,
    query: &[u8],
    subject: &[u8],
    subject_length: i32,
    subject_frame: i32,
    query_info: &crate::queryinfo::QueryInfo,
    align_params: JumperAlignParams,
    jumper: &mut JumperGapAlign,
    hsp_list: &mut crate::hspstream::HspList,
    s_range: u32,
    subject_index: Option<&SubjectIndex>,
) -> i32 {
    let num = word_hits.num.get(index).copied().unwrap_or(0).max(0) as usize;
    if num == 0 {
        return 0;
    }
    let Some(pair_array) = word_hits.pair_arrays.get_mut(index) else {
        return -1;
    };
    let extended = blast_na_extend_jumper(
        &mut pair_array[..num],
        word_params,
        score_params,
        hit_params,
        lookup_wrap,
        query,
        subject,
        subject_length,
        subject_frame,
        query_info,
        align_params,
        jumper,
        hsp_list,
        s_range,
        subject_index,
    );
    word_hits.num[index] = 0;
    extended
}

/// blast-rs: MapperWordHits batching variant of the jumper word finder; not a
/// direct NCBI C port.
///
/// When `word_hits` is present, this mirrors C's per-query-bucket batching:
/// duplicate adjacent hits on the same diagonal are suppressed before batching,
/// full buckets are extended immediately, and remaining buckets are flushed at
/// the end. When `word_hits` is absent, this behaves like
/// [`jumper_na_word_finder`].
#[allow(clippy::too_many_arguments)]
pub fn jumper_na_word_finder_with_word_hits(
    subject: &[u8],
    subject_length: i32,
    subject_frame: i32,
    query: &[u8],
    query_info: &crate::queryinfo::QueryInfo,
    lookup_wrap: &crate::lookup::LookupTableWrap,
    word_params: &crate::parameters::InitialWordParameters,
    score_params: &crate::parameters::ScoringParameters,
    hit_params: &crate::parameters::HitSavingParameters,
    align_params: JumperAlignParams,
    jumper: &mut JumperGapAlign,
    hsp_list: &mut crate::hspstream::HspList,
    word_hits: Option<&mut crate::lookup::MapperWordHits>,
    ungapped_stats: Option<&mut crate::diagnostics::UngappedStats>,
    gapped_stats: Option<&mut crate::diagnostics::GappedStats>,
) -> i16 {
    jumper_na_word_finder_impl(
        subject,
        subject_length,
        subject_frame,
        None,
        query,
        query_info,
        lookup_wrap,
        word_params,
        score_params,
        hit_params,
        align_params,
        jumper,
        hsp_list,
        word_hits,
        ungapped_stats,
        gapped_stats,
    )
}

#[allow(clippy::too_many_arguments)]
/// blast-rs: Shared implementation for Jumper word-finder adapters; not a direct NCBI C port.
fn jumper_na_word_finder_impl(
    subject: &[u8],
    subject_length: i32,
    subject_frame: i32,
    subject_ranges: Option<&[crate::util::SSeqRange]>,
    query: &[u8],
    query_info: &crate::queryinfo::QueryInfo,
    lookup_wrap: &crate::lookup::LookupTableWrap,
    word_params: &crate::parameters::InitialWordParameters,
    score_params: &crate::parameters::ScoringParameters,
    hit_params: &crate::parameters::HitSavingParameters,
    align_params: JumperAlignParams,
    jumper: &mut JumperGapAlign,
    hsp_list: &mut crate::hspstream::HspList,
    word_hits: Option<&mut crate::lookup::MapperWordHits>,
    mut ungapped_stats: Option<&mut crate::diagnostics::UngappedStats>,
    gapped_stats: Option<&mut crate::diagnostics::GappedStats>,
) -> i16 {
    if subject_length <= 0 {
        return 0;
    }
    let Some((word_length, lut_word_length)) = jumper_lookup_lengths(lookup_wrap) else {
        return -1;
    };
    if word_length <= 0 || lut_word_length <= 0 || word_length < lut_word_length {
        return -1;
    }
    if subject_length < lut_word_length {
        return 0;
    }
    let read_is_query = (query_info.max_length as i64) < (subject_length as i64);
    let mut word_hit_batches = jumper_collect_word_hit_batches_in_ranges(
        lookup_wrap,
        subject,
        subject_length,
        word_length,
        lut_word_length,
        subject_ranges,
        read_is_query,
    );
    if !word_hit_batches
        .iter()
        .any(|batch| !batch.offset_pairs.is_empty())
    {
        return 0;
    }

    let subject_index = if std::env::var_os("MAPPER_USE_SMALL_WORDS").is_some() {
        subject_index_new(subject, subject_length, 10000, 4)
    } else {
        None
    };
    let subject_index_ref = subject_index.as_ref();
    let mut total_hits = 0;
    let mut hits_extended = 0;

    if let Some(word_hits) = word_hits {
        let num_arrays = word_hits.num_arrays.max(0) as usize;
        if word_hits.array_size <= 0
            || word_hits.num_arrays <= 0
            || num_arrays > word_hits.pair_arrays.len()
            || num_arrays > word_hits.num.len()
            || word_hits.last_diag.len() < query_info.contexts.len()
            || word_hits.last_pos.len() < query_info.contexts.len()
        {
            return -1;
        }
        for value in &mut word_hits.num {
            *value = 0;
        }
        for value in &mut word_hits.last_pos {
            *value = 0;
        }

        let mut last_s_range = word_hit_batches
            .last()
            .map(|batch| batch.s_range)
            .unwrap_or(subject_length.max(0) as u32);
        for batch in word_hit_batches {
            last_s_range = batch.s_range;
            total_hits += batch.offset_pairs.len() as i32;
            for pair in batch.offset_pairs {
                let q_off = pair.query_offset;
                let s_off = pair.subject_offset;
                if q_off < 0
                    || s_off < 0
                    || q_off as usize >= query.len()
                    || s_off >= subject_length
                {
                    continue;
                }
                let context = crate::queryinfo::bsearch_context_info(q_off, query_info);
                let context_index = context.max(0) as usize;
                if context_index >= query_info.contexts.len() {
                    continue;
                }
                let diag = s_off - q_off;
                let last_d = word_hits.last_diag[context_index];
                let last_p = word_hits.last_pos[context_index];
                word_hits.last_diag[context_index] = diag;
                word_hits.last_pos[context_index] = s_off;
                if last_p != 0
                    && last_d == diag
                    && s_off.saturating_sub(last_p) < lut_word_length.saturating_add(1)
                {
                    continue;
                }

                if word_hits.divisor <= 0 {
                    return -1;
                }
                let index = (q_off / word_hits.divisor) as usize;
                if index >= num_arrays {
                    return -1;
                }
                if word_hits.num[index] >= word_hits.array_size {
                    let extended = mapper_word_hits_flush(
                        word_hits,
                        index,
                        word_params,
                        score_params,
                        hit_params,
                        lookup_wrap,
                        query,
                        subject,
                        subject_length,
                        subject_frame,
                        query_info,
                        align_params,
                        jumper,
                        hsp_list,
                        batch.s_range,
                        subject_index_ref,
                    );
                    if extended < 0 {
                        return -1;
                    }
                    hits_extended += extended;
                }

                let slot = word_hits.num[index].max(0) as usize;
                if slot < word_hits.pair_arrays[index].len() {
                    word_hits.pair_arrays[index][slot] = pair;
                } else {
                    word_hits.pair_arrays[index].push(pair);
                }
                word_hits.num[index] += 1;
            }
        }

        for index in 0..num_arrays {
            let extended = mapper_word_hits_flush(
                word_hits,
                index,
                word_params,
                score_params,
                hit_params,
                lookup_wrap,
                query,
                subject,
                subject_length,
                subject_frame,
                query_info,
                align_params,
                jumper,
                hsp_list,
                last_s_range,
                subject_index_ref,
            );
            if extended < 0 {
                return -1;
            }
            hits_extended += extended;
        }
    } else {
        for batch in &mut word_hit_batches {
            if batch.offset_pairs.is_empty() {
                continue;
            }
            total_hits += batch.offset_pairs.len() as i32;
            hits_extended += blast_na_extend_jumper(
                &mut batch.offset_pairs,
                word_params,
                score_params,
                hit_params,
                lookup_wrap,
                query,
                subject,
                subject_length,
                subject_frame,
                query_info,
                align_params,
                jumper,
                hsp_list,
                batch.s_range,
                subject_index_ref,
            );
        }
    }

    crate::diagnostics::blast_ungapped_stats_update(
        ungapped_stats.as_deref_mut(),
        total_hits,
        0,
        0,
    );
    if let Some(gapped_stats) = gapped_stats {
        gapped_stats.extensions = hits_extended;
        if let Some(ungapped_stats) = ungapped_stats.as_deref_mut() {
            ungapped_stats.good_init_extends = hits_extended;
        }
    }
    0
}

/// blast-rs: Extracts an unambiguous word from unpacked query bases; not a direct NCBI C port.
fn unpacked_query_word(query: &[u8], start: i32, word_size: i32) -> Option<u32> {
    if start < 0 || word_size <= 0 {
        return None;
    }
    let end = start.checked_add(word_size)?;
    if end > query.len() as i32 {
        return None;
    }
    let mut word = 0u32;
    for i in 0..word_size {
        let base = query[(start + i) as usize];
        if base & 0xfc != 0 {
            return None;
        }
        word = (word << 2) | base as u32;
    }
    Some(word)
}

/// blast-rs: Port-shaped equivalent of NCBI `DoAnchoredScan`
/// (jumper.c:4139); not a direct NCBI C port.
#[allow(clippy::too_many_arguments)]
pub fn do_anchored_scan(
    query_seq: &[u8],
    query_len: i32,
    query_from: i32,
    context: i32,
    subject: &[u8],
    subject_len: i32,
    subject_frame: i32,
    subject_from: i32,
    subject_to: i32,
    query_info: &crate::queryinfo::QueryInfo,
    jumper: &mut JumperGapAlign,
    score_params: &crate::parameters::ScoringParameters,
    hit_params: &crate::parameters::HitSavingParameters,
    align_params: JumperAlignParams,
    hsp_list: &mut crate::hspstream::HspList,
) -> i32 {
    const WORD_SIZE: i32 = 12;
    const MAX_NUM_MATCHES: usize = 100;
    if query_len <= 0
        || query_len as usize > query_seq.len()
        || query_from < 0
        || query_from > query_len
        || subject_len < 0
        || subject_len as usize > subject.len().saturating_mul(4)
    {
        return 0;
    }
    let is_right = subject_from < subject_to;
    let big_word_size = if is_right {
        query_len
            .saturating_sub(query_from)
            .saturating_sub(5)
            .max(WORD_SIZE)
            .min(24)
    } else {
        query_from.saturating_sub(5).max(WORD_SIZE).min(24)
    };
    let scan_step = if is_right {
        big_word_size.saturating_sub(WORD_SIZE).saturating_add(1)
    } else {
        big_word_size
            .saturating_sub(WORD_SIZE)
            .saturating_add(1)
            .saturating_neg()
    };
    let scan_to = if is_right {
        subject_to.min(subject_len.saturating_sub(1))
    } else {
        subject_to
    };

    if (is_right
        && (query_len.saturating_sub(query_from).saturating_add(1) < big_word_size
            || scan_to.saturating_sub(subject_from) < big_word_size))
        || (!is_right
            && (query_from < big_word_size || subject_from.saturating_sub(scan_to) < big_word_size))
    {
        return 0;
    }

    let mut words = Vec::<(u32, i32)>::new();
    if is_right {
        let mut q = query_from;
        while q + big_word_size < query_len && words.len() < MAX_NUM_MATCHES {
            while q + big_word_size <= query_len {
                let mut ok = true;
                for i in 0..big_word_size {
                    if query_seq[(q + i) as usize] & 0xfc != 0 {
                        q += i + 1;
                        ok = false;
                        break;
                    }
                }
                if ok {
                    break;
                }
                q += 1;
            }
            if q.saturating_add(big_word_size).saturating_sub(1) >= query_len {
                break;
            }
            if let Some(word) = unpacked_query_word(query_seq, q, WORD_SIZE) {
                if word != 0 && word != 0x00ff_ffff {
                    words.push((word, q));
                }
            }
            q += 1;
        }
    } else {
        let mut q = query_from - big_word_size;
        while q >= 0 && words.len() < MAX_NUM_MATCHES {
            while q >= 0 {
                let mut ok = true;
                for i in 0..big_word_size {
                    if query_seq[(q + i) as usize] & 0xfc != 0 {
                        q = q.saturating_sub(big_word_size).saturating_add(i);
                        ok = false;
                        break;
                    }
                }
                if ok {
                    break;
                }
                q -= 1;
            }
            if q < 0 {
                break;
            }
            if let Some(word) = unpacked_query_word(query_seq, q, WORD_SIZE) {
                if word != 0 && word != 0x00ff_ffff {
                    words.push((word, q));
                }
            }
            q -= 1;
        }
    }

    if words.is_empty() {
        return 0;
    }

    let mut best_score = 0;
    let mut best_hsp = None;
    let mut i = subject_from;
    while (subject_from < scan_to && i < scan_to) || (subject_from > scan_to && i > scan_to) {
        let Some(index) = packed_subject_word(subject, subject_len, i, WORD_SIZE) else {
            i += scan_step;
            continue;
        };
        let Some((_, q_offset)) = words
            .iter()
            .find(|(word, _)| *word as i64 == index)
            .copied()
        else {
            i += scan_step;
            continue;
        };
        let mut k = WORD_SIZE;
        while k < big_word_size {
            if query_seq[(q_offset + k) as usize]
                != crate::encoding::ncbi2na_base_at(subject, (i + k) as usize)
            {
                break;
            }
            k += 1;
        }
        if k < big_word_size {
            i += scan_step;
            continue;
        }

        let mut num_identical = 0;
        let mut right_ungapped_ext_len = 0;
        if jumper_gapped_alignment_compressed_with_traceback(
            query_seq,
            subject,
            query_len,
            subject_len,
            q_offset,
            i,
            jumper,
            score_params,
            align_params,
            &mut num_identical,
            &mut right_ungapped_ext_len,
        ) != 0
        {
            i += scan_step;
            continue;
        }
        let Some((
            score,
            query_align_start,
            query_align_stop,
            subject_align_start,
            subject_align_stop,
        )) = jumper_alignment_outputs(jumper)
        else {
            i += scan_step;
            continue;
        };
        if score > best_score {
            best_score = score;
            let query_frame = query_info
                .contexts
                .get(context.max(0) as usize)
                .map_or(0, |ctx| ctx.frame);
            if let Some((hsp, _)) = s_create_hsp(
                query_seq,
                query_len,
                context,
                query_frame,
                subject,
                subject_len,
                subject_frame,
                jumper,
                query_align_start,
                query_align_stop,
                subject_align_start,
                subject_align_stop,
                score,
                score_params.penalty,
                hit_params.options.splice,
            ) {
                if hsp.score >= query_len.saturating_sub(query_from) {
                    best_hsp = Some(hsp);
                    break;
                }
                best_hsp = Some(hsp);
            }
        }

        i += scan_step;
    }

    if let Some(hsp) = best_hsp {
        let _ = crate::hspstream::blast_hsp_list_save_hsp(hsp_list, hsp);
        1
    } else {
        0
    }
}

/// blast-rs: Port-shaped equivalent of NCBI `DoAnchoredSearch`
/// (jumper.c:4394); not a direct NCBI C port.
///
/// C pulls partially covered mapper chains from the stream writer state, runs
/// `DoAnchoredScan` on uncovered query flanks, and writes a subject HSP list
/// back to the stream. Rust takes the mapper data directly and returns the
/// assembled list so callers and tests can inspect the same side effect.
#[allow(clippy::too_many_arguments)]
pub fn do_anchored_search(
    query: &[u8],
    subject: &[u8],
    subject_len: i32,
    subject_oid: i32,
    subject_frame: i32,
    word_size: i32,
    longest_intron: i32,
    mapper_data: Option<&crate::spliced_hits::BlastHSPMapperData>,
    query_info: &crate::queryinfo::QueryInfo,
    jumper: &mut JumperGapAlign,
    score_params: &crate::parameters::ScoringParameters,
    hit_params: &crate::parameters::HitSavingParameters,
    align_params: JumperAlignParams,
    hsp_stream: Option<&crate::hspstream::HspStream>,
) -> Option<crate::hspstream::HspList> {
    let hsp_list = do_anchored_search_build_list(
        query,
        subject,
        subject_len,
        subject_oid,
        subject_frame,
        word_size,
        longest_intron,
        mapper_data,
        query_info,
        jumper,
        score_params,
        hit_params,
        align_params,
    )?;

    if let Some(stream) = hsp_stream {
        let query_index = query_index_for_hsp_list(&hsp_list, query_info);
        let _ = stream.blast_hspstream_write(query_index, hsp_list.clone());
    }

    if hsp_list.hsps.is_empty() {
        None
    } else {
        Some(hsp_list)
    }
}

/// blast-rs: C-shaped status wrapper around anchored search; not a direct NCBI
/// C port.
///
/// Unlike [`do_anchored_search`], this keeps the C-shaped status contract: a
/// missing stream is an error, and an empty anchored search is a successful
/// stream write attempt. The Rust stream drops empty lists when materializing
/// `HspResults`, so no readable list is retained in that case.
#[allow(clippy::too_many_arguments)]
pub fn do_anchored_search_to_stream(
    query: &[u8],
    subject: &[u8],
    subject_len: i32,
    subject_oid: i32,
    subject_frame: i32,
    word_size: i32,
    longest_intron: i32,
    mapper_data: Option<&crate::spliced_hits::BlastHSPMapperData>,
    query_info: &crate::queryinfo::QueryInfo,
    jumper: &mut JumperGapAlign,
    score_params: &crate::parameters::ScoringParameters,
    hit_params: &crate::parameters::HitSavingParameters,
    align_params: JumperAlignParams,
    hsp_stream: Option<&crate::hspstream::HspStream>,
) -> i16 {
    let Some(stream) = hsp_stream else {
        return -1;
    };
    let Some(hsp_list) = do_anchored_search_build_list(
        query,
        subject,
        subject_len,
        subject_oid,
        subject_frame,
        word_size,
        longest_intron,
        mapper_data,
        query_info,
        jumper,
        score_params,
        hit_params,
        align_params,
    ) else {
        return -1;
    };
    let query_index = query_index_for_hsp_list(&hsp_list, query_info);
    stream.blast_hspstream_write(query_index, hsp_list) as i16
}

/// blast-rs: Maps a recovered HSP list back to a query index; not a direct NCBI C port.
fn query_index_for_hsp_list(
    hsp_list: &crate::hspstream::HspList,
    query_info: &crate::queryinfo::QueryInfo,
) -> i32 {
    hsp_list
        .hsps
        .first()
        .and_then(|hsp| query_info.contexts.get(hsp.context.max(0) as usize))
        .map_or(0, |ctx| ctx.query_index)
}

#[allow(clippy::too_many_arguments)]
/// blast-rs: Builds anchored-search HSP lists from Rust mapper data; not a direct NCBI C port.
fn do_anchored_search_build_list(
    query: &[u8],
    subject: &[u8],
    subject_len: i32,
    subject_oid: i32,
    subject_frame: i32,
    word_size: i32,
    longest_intron: i32,
    mapper_data: Option<&crate::spliced_hits::BlastHSPMapperData>,
    query_info: &crate::queryinfo::QueryInfo,
    jumper: &mut JumperGapAlign,
    score_params: &crate::parameters::ScoringParameters,
    hit_params: &crate::parameters::HitSavingParameters,
    align_params: JumperAlignParams,
) -> Option<crate::hspstream::HspList> {
    let mut hsp_list = crate::hspstream::HspList::new(subject_oid);
    let chains =
        crate::spliced_hits::find_partialy_covered_queries(mapper_data, subject_oid, word_size);
    let Some(chains) = chains else {
        return mapper_data.map(|_| hsp_list);
    };
    let mut chain = Some(&*chains);
    let flank_word_size = word_size.max(0);

    while let Some(current) = chain {
        let Some(first_container) = current.hsps.as_deref() else {
            chain = current.next.as_deref();
            continue;
        };
        let context = first_container.hsp.context;
        let Some(ctx) = query_info.contexts.get(context.max(0) as usize) else {
            chain = current.next.as_deref();
            continue;
        };
        if ctx.query_offset < 0 || ctx.query_length <= 0 {
            chain = current.next.as_deref();
            continue;
        }
        let query_start = ctx.query_offset as usize;
        let query_stop = query_start.saturating_add(ctx.query_length.max(0) as usize);
        let Some(query_seq) = query.get(query_start..query_stop.min(query.len())) else {
            chain = current.next.as_deref();
            continue;
        };
        if query_seq.len() < ctx.query_length as usize {
            chain = current.next.as_deref();
            continue;
        }
        let before_count = hsp_list.hsps.len();
        if first_container.hsp.query_offset > flank_word_size {
            let subject_from = first_container.hsp.subject_offset.saturating_sub(1);
            let subject_to = first_container
                .hsp
                .subject_offset
                .saturating_sub(1)
                .saturating_sub(longest_intron);
            let _ = do_anchored_scan(
                query_seq,
                ctx.query_length,
                first_container.hsp.query_offset.saturating_sub(1),
                context,
                subject,
                subject_len,
                subject_frame,
                subject_from,
                subject_to,
                query_info,
                jumper,
                score_params,
                hit_params,
                align_params,
                &mut hsp_list,
            );
        }

        let mut last_container = first_container;
        while let Some(next) = last_container.next.as_deref() {
            last_container = next;
        }
        if ctx
            .query_length
            .saturating_sub(last_container.hsp.query_end)
            > flank_word_size
        {
            let subject_to = last_container
                .hsp
                .subject_end
                .saturating_add(longest_intron);
            let _ = do_anchored_scan(
                query_seq,
                ctx.query_length,
                last_container.hsp.query_end,
                context,
                subject,
                subject_len,
                subject_frame,
                last_container.hsp.subject_end,
                subject_to,
                query_info,
                jumper,
                score_params,
                hit_params,
                align_params,
                &mut hsp_list,
            );
        }

        if hsp_list.hsps.len() > before_count {
            let mut container = Some(first_container);
            while let Some(hsp_container) = container {
                let _ = crate::hspstream::blast_hsp_list_save_hsp(
                    &mut hsp_list,
                    hsp_container.hsp.clone(),
                );
                container = hsp_container.next.as_deref();
            }
        }

        chain = current.next.as_deref();
    }

    Some(hsp_list)
}

/// blast-rs: Port-shaped equivalent of NCBI `MB_IndexedWordFinder`
/// (na_ungapped.c:1762) for subjects already known to be indexed; not a direct NCBI C port.
///
/// The C function obtains initial hits from indexed database callbacks, then
/// runs the same diagonal hash suppression and ungapped extension filter as the
/// normal megablast path. Rust models the callbacks as closures that fill the
/// supplied [`crate::extend::InitHitList`].
#[allow(clippy::too_many_arguments)]
pub fn mb_indexed_word_finder<CheckOid, GetResults>(
    subject: &[u8],
    subject_len: i32,
    subject_oid: i32,
    subject_chunk: i32,
    query: &[u8],
    query_info: &crate::queryinfo::QueryInfo,
    word_params: &crate::parameters::InitialWordParameters,
    score_params: &crate::parameters::ScoringParameters,
    init_hitlist: &mut crate::extend::InitHitList,
    check_index_oid: CheckOid,
    get_results: GetResults,
    ungapped_stats: Option<&mut crate::diagnostics::UngappedStats>,
) -> i16
where
    CheckOid: FnMut(i32) -> bool,
    GetResults: FnMut(i32, i32, &mut crate::extend::InitHitList) -> u32,
{
    mb_indexed_word_finder_with_fallback(
        subject,
        subject_len,
        subject_oid,
        subject_chunk,
        query,
        query_info,
        word_params,
        score_params,
        init_hitlist,
        check_index_oid,
        get_results,
        || 0,
        ungapped_stats,
    )
}

/// blast-rs: C-shaped indexed megablast dispatcher with explicit fallback; not
/// a direct NCBI C port.
///
/// Upstream falls back to `BlastNaWordFinder` when the indexed database reports
/// that an OID is not indexed. The narrow [`mb_indexed_word_finder`] helper
/// preserves its historical return-zero behavior by passing an empty fallback;
/// callers that carry the normal nucleotide word-finder state should use this
/// variant and route `fallback` to their represented `BlastNaWordFinder` path.
#[allow(clippy::too_many_arguments)]
pub fn mb_indexed_word_finder_with_fallback<CheckOid, GetResults, Fallback>(
    subject: &[u8],
    subject_len: i32,
    subject_oid: i32,
    subject_chunk: i32,
    query: &[u8],
    query_info: &crate::queryinfo::QueryInfo,
    word_params: &crate::parameters::InitialWordParameters,
    score_params: &crate::parameters::ScoringParameters,
    init_hitlist: &mut crate::extend::InitHitList,
    mut check_index_oid: CheckOid,
    mut get_results: GetResults,
    mut fallback: Fallback,
    ungapped_stats: Option<&mut crate::diagnostics::UngappedStats>,
) -> i16
where
    CheckOid: FnMut(i32) -> bool,
    GetResults: FnMut(i32, i32, &mut crate::extend::InitHitList) -> u32,
    Fallback: FnMut() -> i16,
{
    if !check_index_oid(subject_oid) {
        return fallback();
    }

    init_hitlist.reset();
    let word_size = get_results(subject_oid, subject_chunk, init_hitlist);
    let total_hits = init_hitlist.total() as i32;
    let mut hits_extended = 0;
    let mut saved_hits = 0;

    if subject_len < 0 || subject_len as usize > subject.len().saturating_mul(4) {
        if word_params.ungapped_extension {
            init_hitlist.hits.clear();
            crate::diagnostics::blast_ungapped_stats_update(ungapped_stats, total_hits, 0, 0);
        }
        return 0;
    }

    if word_size > 0 && word_params.ungapped_extension {
        let mut hash = crate::index_ungapped::ir_hash_create().expect("indexed diagonal hash");
        let mut kept = Vec::with_capacity(init_hitlist.hits.len());

        for hsp in init_hitlist.hits.iter().cloned() {
            let q_off = hsp.query_offset.max(0) as u32;
            let s_off = hsp.subject_offset.max(0) as u32;
            let diag = crate::index_ungapped::ir_diag(q_off, s_off);
            let key = crate::index_ungapped::ir_key(diag);
            let Some(entry_index) = crate::index_ungapped::ir_locate_macro(&mut hash, diag, key)
            else {
                kept.push(hsp);
                continue;
            };
            let qend = crate::index_ungapped::ir_hash_entry_mut(&mut hash, entry_index)
                .map_or(0, |entry| entry.diag_data.qend);
            let hit_qend = q_off
                .checked_add(word_size.saturating_sub(1))
                .unwrap_or(u32::MAX);
            if hit_qend <= qend {
                continue;
            }

            let context = crate::queryinfo::bsearch_context_info(hsp.query_offset, query_info);
            let cutoffs = word_params.cutoffs.get(context.max(0) as usize);
            let x_dropoff = cutoffs
                .map(|cutoff| cutoff.x_dropoff_init)
                .unwrap_or(word_params.x_dropoff_max)
                .max(0);
            hits_extended += 1;

            if hsp.query_offset < 0
                || hsp.subject_offset < 0
                || hsp.query_offset >= query.len() as i32
                || hsp.subject_offset >= subject_len
                || hsp.subject_offset as usize >= subject.len().saturating_mul(4)
            {
                if let Some(entry) =
                    crate::index_ungapped::ir_hash_entry_mut(&mut hash, entry_index)
                {
                    entry.diag_data.diag = diag;
                    entry.diag_data.qend = hit_qend;
                }
                continue;
            }

            let Some(ungapped_data) = crate::extend::na_ungapped_extend_len(
                query,
                subject,
                subject_len.max(0) as usize,
                hsp.query_offset,
                hsp.subject_offset,
                score_params.reward,
                score_params.penalty,
                x_dropoff,
            ) else {
                if let Some(entry) =
                    crate::index_ungapped::ir_hash_entry_mut(&mut hash, entry_index)
                {
                    entry.diag_data.diag = diag;
                    entry.diag_data.qend = hit_qend;
                }
                continue;
            };

            let cutoff_score = cutoffs
                .map(|cutoff| cutoff.cutoff_score)
                .unwrap_or(word_params.cutoff_score_min);
            let qend_new = i64::from(ungapped_data.q_start)
                .saturating_add(i64::from(ungapped_data.length))
                .saturating_sub(1)
                .clamp(0, i64::from(u32::MAX)) as u32;
            if let Some(entry) = crate::index_ungapped::ir_hash_entry_mut(&mut hash, entry_index) {
                entry.diag_data.diag = diag;
                entry.diag_data.qend = qend_new;
            }

            if ungapped_data.score >= cutoff_score {
                let mut kept_hsp = hsp;
                kept_hsp.ungapped_data = Some(ungapped_data);
                kept.push(kept_hsp);
                saved_hits += 1;
            }
        }

        init_hitlist.hits = kept;
        let _ = crate::index_ungapped::ir_hash_destroy(Some(hash));
    }

    if word_params.ungapped_extension {
        crate::diagnostics::blast_ungapped_stats_update(
            ungapped_stats,
            total_hits,
            hits_extended,
            saved_hits,
        );
        crate::extend::blast_init_hit_list_sort_by_score(init_hitlist);
    }
    0
}

/// Query-side native MegaBLAST index lookup key.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MegablastIndexQueryKey {
    pub hash_key: u32,
    pub query_offset: i32,
}

/// Minimal source interface for native MegaBLAST indexed offset-list hits.
pub trait MegablastIndexedOffsetSource {
    fn hkey_width(&self) -> u32;
    fn start_oid(&self) -> u32;
    fn stop_oid(&self) -> u32;
    fn offset_list_words(&self, hash_key: u32) -> Option<&[u32]>;
    fn logical_chunk_id_for_subject_chunk(&self, local_oid: u32, chunk: u32) -> Option<u32>;
    fn chunk_offset_for_encoded_offset(&self, encoded_offset: u32) -> Option<(u32, u32)>;
    fn chunk_seed_offset_for_encoded_offset(&self, encoded_offset: u32) -> Option<(u32, u32)>;
    fn min_encoded_offset(&self) -> Option<u32>;
}

/// Return whether a database OID is covered by a parsed native MegaBLAST index
/// volume.
pub fn megablast_index_volume_has_oid(
    volume: &impl MegablastIndexedOffsetSource,
    subject_oid: i32,
) -> bool {
    let Ok(subject_oid) = u32::try_from(subject_oid) else {
        return false;
    };
    subject_oid >= volume.start_oid() && subject_oid < volume.stop_oid()
}

/// Build query hash keys for a parsed native MegaBLAST index volume.
///
/// The native `.idx` offset lists are keyed by `hkey_width` 2-bit nucleotide
/// windows. Ambiguous query bytes are skipped, matching the indexed search
/// precondition before seed extension.
pub fn megablast_index_query_keys(query: &[u8], hkey_width: u32) -> Vec<MegablastIndexQueryKey> {
    let right = query.len().saturating_sub(1);
    let all_query = [crate::util::SSeqRange {
        left: 0,
        right: right as i32,
    }];
    megablast_index_query_keys_for_locations(query, hkey_width, &all_query)
}

/// Build native MegaBLAST index query keys for explicit unmasked query ranges.
///
/// NCBI's indexed pre-search consumes the same query locations used to build
/// ordinary nucleotide lookup tables. Ranges are inclusive `SSeqRange`
/// intervals; seed windows are never allowed to cross a range boundary.
pub fn megablast_index_query_keys_for_locations(
    query: &[u8],
    hkey_width: u32,
    locations: &[crate::util::SSeqRange],
) -> Vec<MegablastIndexQueryKey> {
    let Ok(width) = usize::try_from(hkey_width) else {
        return Vec::new();
    };
    if width == 0 || width > 15 || query.len() < width {
        return Vec::new();
    }

    let mut keys = Vec::new();
    for loc in locations {
        if loc.left < 0 || loc.right < loc.left {
            continue;
        }
        let from = loc.left as usize;
        let right = loc.right as usize;
        let loc_len = right - from + 1;
        if loc_len < width || from >= query.len() {
            continue;
        }
        let Some(last_start) = right.checked_sub(width - 1) else {
            continue;
        };
        for start in from..=last_start {
            let Some(word) = query.get(start..start + width) else {
                break;
            };
            let mut hash_key = 0u32;
            let mut valid = true;
            for &base in word {
                if base > 3 {
                    valid = false;
                    break;
                }
                hash_key = (hash_key << 2) | u32::from(base);
            }
            if valid {
                keys.push(MegablastIndexQueryKey {
                    hash_key,
                    query_offset: start as i32,
                });
            }
        }
    }
    keys
}

/// Callback adapter over one or more parsed native MegaBLAST index volumes.
///
/// NCBI's indexed search keeps a result object across all index volumes and
/// answers per-subject `CheckOid`/`GetResults` calls from that holder. This
/// struct is the parsed Rust boundary for the same shape: query keys are cached
/// by hash-key width, OID lookup selects the owning volume, and result filling
/// delegates to the single-volume offset-list adapter.
pub struct MegablastIndexResultHolder<'a, V: MegablastIndexedOffsetSource> {
    volumes: Vec<&'a V>,
    query_keys_by_width: Vec<(u32, Vec<MegablastIndexQueryKey>)>,
}

impl<'a, V: MegablastIndexedOffsetSource> MegablastIndexResultHolder<'a, V> {
    pub fn new(query: &[u8], volumes: impl IntoIterator<Item = &'a V>) -> Self {
        let right = query.len().saturating_sub(1);
        let locations = [crate::util::SSeqRange {
            left: 0,
            right: right as i32,
        }];
        Self::with_locations(query, &locations, volumes)
    }

    pub fn with_locations(
        query: &[u8],
        locations: &[crate::util::SSeqRange],
        volumes: impl IntoIterator<Item = &'a V>,
    ) -> Self {
        let volumes: Vec<&'a V> = volumes.into_iter().collect();
        let mut query_keys_by_width = Vec::new();
        for volume in &volumes {
            let width = volume.hkey_width();
            if !query_keys_by_width
                .iter()
                .any(|(known_width, _)| *known_width == width)
            {
                query_keys_by_width.push((
                    width,
                    megablast_index_query_keys_for_locations(query, width, locations),
                ));
            }
        }

        Self {
            volumes,
            query_keys_by_width,
        }
    }

    pub fn has_oid(&self, subject_oid: i32) -> bool {
        self.volume_for_oid(subject_oid).is_some()
    }

    pub fn fill_init_hitlist(
        &self,
        subject_oid: i32,
        subject_chunk: i32,
        init_hitlist: &mut crate::extend::InitHitList,
    ) -> u32 {
        let Some(volume) = self.volume_for_oid(subject_oid) else {
            init_hitlist.reset();
            return 0;
        };
        let Some((_, query_keys)) = self
            .query_keys_by_width
            .iter()
            .find(|(width, _)| *width == volume.hkey_width())
        else {
            init_hitlist.reset();
            return 0;
        };
        megablast_index_fill_init_hitlist(
            volume,
            subject_oid,
            subject_chunk,
            query_keys,
            init_hitlist,
        )
    }

    fn volume_for_oid(&self, subject_oid: i32) -> Option<&'a V> {
        self.volumes
            .iter()
            .copied()
            .find(|volume| megablast_index_volume_has_oid(*volume, subject_oid))
    }
}

/// Fill an `InitHitList` from parsed native MegaBLAST offset lists for the
/// existing `MB_IndexedWordFinder` callback shape.
///
/// This is a deliberately narrow adapter over decoded offset-list words. It
/// maps raw index offsets through the parsed `CSubjectMap` data and saves only
/// hits belonging to the requested local subject/chunk. Boundary control words
/// below NCBI `GetMinOffset(stride)` are consumed with their following real
/// offset so the native list shape matches `CSeedRoots::Add2` at this callback
/// boundary; full pre-search tracked-seed extension remains above this adapter.
pub fn megablast_index_fill_init_hitlist(
    volume: &impl MegablastIndexedOffsetSource,
    subject_oid: i32,
    subject_chunk: i32,
    query_keys: &[MegablastIndexQueryKey],
    init_hitlist: &mut crate::extend::InitHitList,
) -> u32 {
    init_hitlist.reset();
    let Ok(subject_oid) = u32::try_from(subject_oid) else {
        return 0;
    };
    let Some(local_oid) = subject_oid.checked_sub(volume.start_oid()) else {
        return 0;
    };
    let Ok(subject_chunk) = u32::try_from(subject_chunk) else {
        return 0;
    };
    let Some(target_chunk_id) = volume.logical_chunk_id_for_subject_chunk(local_oid, subject_chunk)
    else {
        return 0;
    };
    let Some(min_encoded_offset) = volume.min_encoded_offset() else {
        return 0;
    };

    let mut saved = 0u32;
    for query_key in query_keys {
        let Some(offsets) = volume.offset_list_words(query_key.hash_key) else {
            continue;
        };
        let mut index = 0usize;
        while index < offsets.len() {
            let encoded_offset = offsets[index];
            let real_encoded_offset = if encoded_offset < min_encoded_offset {
                let Some(&next_offset) = offsets.get(index + 1) else {
                    break;
                };
                index += 2;
                if next_offset < min_encoded_offset {
                    continue;
                }
                next_offset
            } else {
                index += 1;
                encoded_offset
            };

            let Some((logical_chunk_id, chunk_offset)) =
                volume.chunk_seed_offset_for_encoded_offset(real_encoded_offset)
            else {
                continue;
            };
            if logical_chunk_id == target_chunk_id {
                if crate::extend::blast_save_initial_hit(
                    init_hitlist,
                    query_key.query_offset,
                    chunk_offset as i32,
                    None,
                ) {
                    saved = saved.saturating_add(1);
                }
            }
        }
    }

    if saved == 0 {
        0
    } else {
        volume.hkey_width()
    }
}

/// blast-rs: Port-shaped equivalent of NCBI `ShortRead_IndexedWordFinder`
/// (na_ungapped.c:2227) for subjects already known to be indexed; not a direct NCBI C port.
///
/// For indexed subjects, C consumes indexed initial hits and immediately runs
/// Jumper gapped traceback per non-redundant diagonal.
#[allow(clippy::too_many_arguments)]
pub fn short_read_indexed_word_finder<CheckOid, GetResults>(
    subject: &[u8],
    subject_len: i32,
    subject_frame: i32,
    subject_oid: i32,
    subject_chunk: i32,
    query: &[u8],
    query_info: &crate::queryinfo::QueryInfo,
    score_params: &crate::parameters::ScoringParameters,
    hit_params: &crate::parameters::HitSavingParameters,
    align_params: JumperAlignParams,
    init_hitlist: &mut crate::extend::InitHitList,
    hsp_list: &mut crate::hspstream::HspList,
    jumper: &mut JumperGapAlign,
    check_index_oid: CheckOid,
    get_results: GetResults,
    gapped_stats: Option<&mut crate::diagnostics::GappedStats>,
) -> i16
where
    CheckOid: FnMut(i32) -> bool,
    GetResults: FnMut(i32, i32, &mut crate::extend::InitHitList) -> u32,
{
    short_read_indexed_word_finder_with_fallback(
        subject,
        subject_len,
        subject_frame,
        subject_oid,
        subject_chunk,
        query,
        query_info,
        score_params,
        hit_params,
        align_params,
        init_hitlist,
        hsp_list,
        jumper,
        check_index_oid,
        get_results,
        || 0,
        gapped_stats,
    )
}

/// blast-rs: C-shaped indexed short-read dispatcher with explicit fallback; not
/// a direct NCBI C port.
///
/// Upstream falls back to `JumperNaWordFinder` when an OID is not indexed. This
/// wrapper exposes that branch explicitly without forcing the indexed callback
/// helper to know about lookup-table state.
#[allow(clippy::too_many_arguments)]
pub fn short_read_indexed_word_finder_with_fallback<CheckOid, GetResults, Fallback>(
    subject: &[u8],
    subject_len: i32,
    subject_frame: i32,
    subject_oid: i32,
    subject_chunk: i32,
    query: &[u8],
    query_info: &crate::queryinfo::QueryInfo,
    score_params: &crate::parameters::ScoringParameters,
    hit_params: &crate::parameters::HitSavingParameters,
    align_params: JumperAlignParams,
    init_hitlist: &mut crate::extend::InitHitList,
    hsp_list: &mut crate::hspstream::HspList,
    jumper: &mut JumperGapAlign,
    mut check_index_oid: CheckOid,
    mut get_results: GetResults,
    mut fallback: Fallback,
    gapped_stats: Option<&mut crate::diagnostics::GappedStats>,
) -> i16
where
    CheckOid: FnMut(i32) -> bool,
    GetResults: FnMut(i32, i32, &mut crate::extend::InitHitList) -> u32,
    Fallback: FnMut() -> i16,
{
    if !check_index_oid(subject_oid) {
        return fallback();
    }

    init_hitlist.reset();
    let word_size = get_results(subject_oid, subject_chunk, init_hitlist);
    if word_size == 0 {
        return 0;
    }
    if subject_len < 0 || subject_len as usize > subject.len().saturating_mul(4) {
        return 0;
    }

    let mut hash = crate::index_ungapped::ir_hash_create().expect("indexed diagonal hash");
    let mut extensions = 0;
    let mut good_extensions = 0;

    for hsp in init_hitlist.hits.iter() {
        let q_off = hsp.query_offset.max(0) as u32;
        let s_off = hsp.subject_offset.max(0) as u32;
        let diag = crate::index_ungapped::ir_diag(q_off, s_off);
        let key = crate::index_ungapped::ir_key(diag);
        let Some(entry_index) = crate::index_ungapped::ir_locate_macro(&mut hash, diag, key) else {
            continue;
        };
        let qend = crate::index_ungapped::ir_hash_entry_mut(&mut hash, entry_index)
            .map_or(0, |entry| entry.diag_data.qend);
        let hit_qend = q_off
            .checked_add(word_size.saturating_sub(1))
            .unwrap_or(u32::MAX);
        if hit_qend <= qend {
            continue;
        }

        let context = crate::queryinfo::bsearch_context_info(hsp.query_offset, query_info);
        let Some(ctx) = query_info.contexts.get(context.max(0) as usize) else {
            continue;
        };
        let query_start = ctx.query_offset.max(0) as usize;
        let query_stop = query_start.saturating_add(ctx.query_length.max(0) as usize);
        let Some(query_seq) = query.get(query_start..query_stop.min(query.len())) else {
            continue;
        };
        if query_seq.len() < ctx.query_length.max(0) as usize {
            continue;
        }

        let mut num_identical = 0;
        let mut right_ungapped_ext_len = 0;
        extensions += 1;
        if jumper_gapped_alignment_compressed_with_traceback(
            query_seq,
            subject,
            ctx.query_length,
            subject_len,
            hsp.query_offset.saturating_sub(ctx.query_offset),
            hsp.subject_offset,
            jumper,
            score_params,
            align_params,
            &mut num_identical,
            &mut right_ungapped_ext_len,
        ) == 0
        {
            if let Some((
                score,
                query_align_start,
                query_align_stop,
                subject_align_start,
                subject_align_stop,
            )) = jumper_alignment_outputs(jumper)
            {
                if jumper_good_align(
                    score,
                    query_align_start,
                    query_align_stop,
                    subject_align_start,
                    subject_align_stop,
                    num_identical,
                    &hit_params.options,
                    ctx.query_length,
                ) {
                    if let Some(new_hsp) = s_create_short_read_indexed_hsp(
                        query_seq,
                        ctx.query_length,
                        context,
                        ctx.frame,
                        subject,
                        subject_len,
                        subject_frame,
                        jumper,
                        query_align_start,
                        query_align_stop,
                        subject_align_start,
                        subject_align_stop,
                        hsp.query_offset.saturating_sub(ctx.query_offset),
                        hsp.subject_offset,
                        score,
                        num_identical,
                        hit_params.options.splice,
                    ) {
                        let _ = crate::hspstream::blast_hsp_list_save_hsp(hsp_list, new_hsp);
                        good_extensions += 1;
                    }
                }
            }
        }

        if let Some(entry) = crate::index_ungapped::ir_hash_entry_mut(&mut hash, entry_index) {
            entry.diag_data.diag = diag;
            entry.diag_data.qend = i64::from(hsp.query_offset)
                .saturating_add(i64::from(right_ungapped_ext_len))
                .saturating_sub(1)
                .clamp(0, i64::from(u32::MAX)) as u32;
        }
    }

    if let Some(stats) = gapped_stats {
        stats.extensions += extensions;
        stats.good_extensions += good_extensions;
    }
    let _ = crate::index_ungapped::ir_hash_destroy(Some(hash));
    0
}

impl GapEditScript {
    /// blast-rs: Rust constructor wrapper; not a direct NCBI C port.
    pub fn new() -> Self {
        Self::default()
    }

    /// blast-rs: Rust constructor wrapper with reserved storage; not a direct NCBI C port.
    pub fn with_capacity(size: usize) -> Self {
        GapEditScript {
            ops: Vec::with_capacity(size),
        }
    }

    /// blast-rs: Append helper that preserves C edit-script merge semantics; not a direct NCBI C port.
    ///
    /// Append an edit op, merging with the previous op if the type matches.
    /// Mirrors NCBI's `Blast_PrelimEditBlockToGapEditScript`
    /// (`blast_gapalign.c:2482`) which collapses consecutive same-type ops
    /// when concatenating the reverse + forward halves into a single
    /// edit script. Without this merge, our left-half-end and right-half-
    /// start SUB runs become two ops (`3:6 3:14`) instead of NCBI's
    /// single op (`3:20`), and downstream reevaluation walks them
    /// differently — affecting Kadane's-best-subarray decisions.
    pub fn push(&mut self, op: GapAlignOpType, count: i32) {
        if let Some(last) = self.ops.last_mut() {
            if last.0 == op {
                last.1 += count;
                return;
            }
        }
        self.ops.push((op, count));
    }

    /// blast-rs: Native edit-script summary helper; not a direct NCBI C port.
    ///
    /// Total alignment length (sum of all op counts).
    pub fn alignment_length(&self) -> i32 {
        self.ops.iter().map(|(_, n)| *n).sum()
    }

    /// blast-rs: Native alignment rendering helper; not a direct NCBI C port.
    ///
    /// Render aligned query and subject strings from edit script.
    /// `query` and `subject` are the byte slices covering the aligned region.
    /// `to_char` converts a single encoded byte to a display character.
    pub fn render_alignment(
        &self,
        query: &[u8],
        subject: &[u8],
        to_char: fn(u8) -> char,
    ) -> (String, String) {
        let mut q_str = String::new();
        let mut s_str = String::new();
        let mut q_pos = 0usize;
        let mut s_pos = 0usize;

        for &(op, count) in &self.ops {
            let count = count as usize;
            match op {
                GapAlignOpType::Sub => {
                    for _ in 0..count {
                        q_str.push(if q_pos < query.len() {
                            to_char(query[q_pos])
                        } else {
                            'N'
                        });
                        s_str.push(if s_pos < subject.len() {
                            to_char(subject[s_pos])
                        } else {
                            'N'
                        });
                        q_pos += 1;
                        s_pos += 1;
                    }
                }
                GapAlignOpType::Del | GapAlignOpType::Del1 | GapAlignOpType::Del2 => {
                    // Gap in query, bases in subject
                    for _ in 0..count {
                        q_str.push('-');
                        s_str.push(if s_pos < subject.len() {
                            to_char(subject[s_pos])
                        } else {
                            'N'
                        });
                        s_pos += 1;
                    }
                }
                GapAlignOpType::Ins | GapAlignOpType::Ins1 | GapAlignOpType::Ins2 => {
                    // Bases in query, gap in subject
                    for _ in 0..count {
                        q_str.push(if q_pos < query.len() {
                            to_char(query[q_pos])
                        } else {
                            'N'
                        });
                        s_str.push('-');
                        q_pos += 1;
                    }
                }
                _ => {
                    // Decline or unknown — treat as substitution
                    for _ in 0..count {
                        q_str.push(if q_pos < query.len() {
                            to_char(query[q_pos])
                        } else {
                            'N'
                        });
                        s_str.push(if s_pos < subject.len() {
                            to_char(subject[s_pos])
                        } else {
                            'N'
                        });
                        q_pos += 1;
                        s_pos += 1;
                    }
                }
            }
        }

        (q_str, s_str)
    }

    /// blast-rs: Native edit-script identity counter; not a direct NCBI C port.
    ///
    /// Count identities given query and subject byte slices.
    pub fn count_identities(&self, query: &[u8], subject: &[u8]) -> (i32, i32, i32) {
        let mut q_pos = 0usize;
        let mut s_pos = 0usize;
        let mut align_len = 0i32;
        let mut num_ident = 0i32;
        let mut gap_opens = 0i32;

        for &(op, count) in &self.ops {
            let count = count as usize;
            align_len += count as i32;

            match op {
                GapAlignOpType::Sub => {
                    for _ in 0..count {
                        if q_pos < query.len()
                            && s_pos < subject.len()
                            && query[q_pos] == subject[s_pos]
                        {
                            num_ident += 1;
                        }
                        q_pos += 1;
                        s_pos += 1;
                    }
                }
                GapAlignOpType::Del | GapAlignOpType::Del1 | GapAlignOpType::Del2 => {
                    s_pos += count;
                    gap_opens += 1;
                }
                GapAlignOpType::Ins | GapAlignOpType::Ins1 | GapAlignOpType::Ins2 => {
                    q_pos += count;
                    gap_opens += 1;
                }
                _ => {
                    q_pos += count;
                    s_pos += count;
                }
            }
        }
        (align_len, num_ident, gap_opens)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_perfect_match() {
        let mut esp = GapEditScript::new();
        esp.push(GapAlignOpType::Sub, 10);
        let query = &[0u8, 1, 2, 3, 0, 1, 2, 3, 0, 1];
        let subject = &[0u8, 1, 2, 3, 0, 1, 2, 3, 0, 1];
        let (len, ident, gaps) = esp.count_identities(query, subject);
        assert_eq!(len, 10);
        assert_eq!(ident, 10);
        assert_eq!(gaps, 0);
    }

    #[test]
    fn test_with_gap() {
        let mut esp = GapEditScript::new();
        esp.push(GapAlignOpType::Sub, 5);
        esp.push(GapAlignOpType::Del, 2);
        esp.push(GapAlignOpType::Sub, 3);
        let query = &[0u8, 1, 2, 3, 0, 1, 2, 3];
        let subject = &[0u8, 1, 2, 3, 0, 9, 9, 1, 2, 3];
        let (len, ident, gaps) = esp.count_identities(query, subject);
        assert_eq!(len, 10);
        assert_eq!(ident, 8);
        assert_eq!(gaps, 1);
    }

    #[test]
    fn gapinfo_c_lifecycle_helpers_match_allocation_and_copy_rules() {
        let state = Some(Box::new(GapStateArrayStruct {
            state_array: vec![1, 2, 3],
            length: 3,
            used: true,
            next: Some(Box::new(GapStateArrayStruct {
                state_array: vec![4],
                length: 1,
                used: false,
                next: None,
            })),
        }));
        assert!(gap_state_free(state).is_none());

        let mut state_pool = None;
        {
            let state = s_gap_get_state(&mut state_pool, 16).expect("allocated state");
            assert!(state.used);
            assert!(state.length as usize >= GAP_STATE_MIN_CELLS);
            assert_eq!(state.state_array.len(), state.length as usize);
        }
        {
            let state = s_gap_get_state(&mut state_pool, 32).expect("second state");
            assert!(state.used);
            assert!(state.length as usize >= GAP_STATE_MIN_CELLS);
        }
        let mut first = state_pool.as_mut().expect("head state");
        assert!(first.used);
        assert!(first.next.as_ref().expect("next state").used);
        assert!(s_gap_purge_state(&mut first));
        assert!(!first.used);
        assert!(!first.next.as_ref().expect("next state").used);
        let reused = s_gap_get_state(&mut state_pool, 8).expect("reused state");
        assert!(reused.used);

        assert!(gap_edit_script_new(0).is_none());
        let mut source = gap_edit_script_new(3).expect("source script");
        source.push(GapAlignOpType::Sub, 5);
        source.push(GapAlignOpType::Del, 2);
        source.push(GapAlignOpType::Ins, 1);

        let duplicate = gap_edit_script_dup(Some(&source)).expect("duplicate script");
        assert_eq!(duplicate.ops, source.ops);
        assert!(gap_edit_script_dup(None).is_none());

        let mut dest = gap_edit_script_new(3).expect("destination script");
        assert_eq!(
            gap_edit_script_partial_copy(&mut dest, 1, Some(&source), 0, 1),
            0
        );
        assert_eq!(dest.ops[1..3], source.ops[0..2]);
        assert_eq!(
            gap_edit_script_partial_copy(&mut dest, 2, Some(&source), 0, 2),
            -1
        );
        assert_eq!(
            gap_edit_script_partial_copy(&mut dest, 0, Some(&source), 2, 1),
            -1
        );

        assert!(gap_edit_script_delete(Some(source)).is_none());
        assert!(gap_edit_script_delete(None).is_none());
    }

    #[test]
    fn prelim_edit_block_helpers_merge_reset_append_and_reserve() {
        let mut block = gap_prelim_edit_block_new();
        assert!(block.edit_ops.capacity() >= 100);
        assert_eq!(s_gap_prelim_edit_block_realloc(&mut block, -1), -1);
        assert_eq!(s_gap_prelim_edit_block_realloc(&mut block, 150), 0);
        assert!(block.edit_ops.capacity() >= 150);

        gap_prelim_edit_block_add(&mut block, GapAlignOpType::Sub, 0);
        assert!(block.edit_ops.is_empty());
        gap_prelim_edit_block_add(&mut block, GapAlignOpType::Sub, 5);
        gap_prelim_edit_block_add(&mut block, GapAlignOpType::Sub, 3);
        gap_prelim_edit_block_add(&mut block, GapAlignOpType::Del, 2);
        assert_eq!(
            block.edit_ops,
            vec![
                GapPrelimEditScript {
                    op_type: GapAlignOpType::Sub,
                    num: 8,
                },
                GapPrelimEditScript {
                    op_type: GapAlignOpType::Del,
                    num: 2,
                },
            ]
        );

        let mut other = gap_prelim_edit_block_new();
        gap_prelim_edit_block_add(&mut other, GapAlignOpType::Del, 4);
        gap_prelim_edit_block_add(&mut other, GapAlignOpType::Ins, 1);
        gap_prelim_edit_block_append(&mut block, &other);
        assert_eq!(block.edit_ops[1].num, 6);
        assert_eq!(
            block.edit_ops.last().map(|op| (op.op_type, op.num)),
            Some((GapAlignOpType::Ins, 1))
        );

        gap_prelim_edit_block_reset(Some(&mut block));
        assert!(block.edit_ops.is_empty());
        assert_eq!(block.last_op, None);
        gap_prelim_edit_block_reset(None);
        assert!(gap_prelim_edit_block_free(Some(block)).is_none());
    }

    #[test]
    fn translated_phi_gap_cost_matches_c_formula() {
        assert_eq!(s_gap_cost(5, 2, -3), 0);
        assert_eq!(s_gap_cost(5, 2, 0), 0);
        assert_eq!(s_gap_cost(5, 2, 1), 7);
        assert_eq!(s_gap_cost(5, 2, 4), 13);
    }

    fn phi_test_matrix() -> Vec<Vec<i32>> {
        let mut matrix = vec![vec![-3; 8]; 8];
        for i in 0..8 {
            matrix[i][i] = 2;
        }
        matrix
    }

    #[test]
    fn translated_phi_banded_align_matches_c_boundary_cases() {
        let matrix = phi_test_matrix();

        let mut script = gap_prelim_edit_block_new();
        assert_eq!(
            s_banded_align(&[0, 1, 2, 3], &[0], 3, 0, -4, 4, &matrix, 5, 2, &mut script),
            -11
        );
        assert_eq!(
            script.edit_ops,
            vec![GapPrelimEditScript {
                op_type: GapAlignOpType::Ins,
                num: 3,
            }]
        );

        let mut script = gap_prelim_edit_block_new();
        assert_eq!(
            s_banded_align(&[0], &[0, 1, 2], 0, 2, -4, 4, &matrix, 5, 2, &mut script),
            -9
        );
        assert_eq!(
            script.edit_ops,
            vec![GapPrelimEditScript {
                op_type: GapAlignOpType::Del,
                num: 2,
            }]
        );
    }

    #[test]
    fn translated_phi_banded_align_scores_single_diagonal_and_dp_paths() {
        let matrix = phi_test_matrix();
        let seq = [0, 1, 2, 3, 4];
        let mut script = gap_prelim_edit_block_new();

        assert_eq!(
            s_banded_align(&seq, &seq, 4, 4, 0, 0, &matrix, 5, 2, &mut script),
            8
        );
        assert_eq!(
            script.edit_ops,
            vec![GapPrelimEditScript {
                op_type: GapAlignOpType::Sub,
                num: 4,
            }]
        );

        let seq1 = [0, 1, 2, 3];
        let seq2 = [0, 1, 7, 2, 3];
        let mut script = gap_prelim_edit_block_new();
        assert_eq!(
            s_banded_align(&seq1, &seq2, 3, 4, -1, 1, &matrix, 5, 2, &mut script),
            -1
        );
        assert_eq!(
            script.edit_ops,
            vec![
                GapPrelimEditScript {
                    op_type: GapAlignOpType::Sub,
                    num: 1,
                },
                GapPrelimEditScript {
                    op_type: GapAlignOpType::Del,
                    num: 1,
                },
                GapPrelimEditScript {
                    op_type: GapAlignOpType::Sub,
                    num: 2,
                },
            ]
        );
    }

    #[test]
    fn blast_prelim_edit_block_to_gap_edit_script_reverses_forward_half() {
        let rev = GapPrelimEditBlock {
            edit_ops: vec![
                GapPrelimEditScript {
                    op_type: GapAlignOpType::Sub,
                    num: 4,
                },
                GapPrelimEditScript {
                    op_type: GapAlignOpType::Del,
                    num: 2,
                },
            ],
            last_op: Some(GapAlignOpType::Del),
        };
        let fwd = GapPrelimEditBlock {
            edit_ops: vec![
                GapPrelimEditScript {
                    op_type: GapAlignOpType::Ins,
                    num: 1,
                },
                GapPrelimEditScript {
                    op_type: GapAlignOpType::Sub,
                    num: 3,
                },
            ],
            last_op: Some(GapAlignOpType::Sub),
        };

        let script = blast_prelim_edit_block_to_gap_edit_script(Some(&rev), Some(&fwd))
            .expect("edit script");

        assert_eq!(
            script.ops,
            vec![
                (GapAlignOpType::Sub, 4),
                (GapAlignOpType::Del, 2),
                (GapAlignOpType::Sub, 3),
                (GapAlignOpType::Ins, 1),
            ]
        );
        assert!(blast_prelim_edit_block_to_gap_edit_script(None, Some(&fwd)).is_none());
    }

    #[test]
    fn blast_prelim_edit_block_to_gap_edit_script_merges_boundary_op() {
        let rev = GapPrelimEditBlock {
            edit_ops: vec![GapPrelimEditScript {
                op_type: GapAlignOpType::Sub,
                num: 4,
            }],
            last_op: Some(GapAlignOpType::Sub),
        };
        let fwd = GapPrelimEditBlock {
            edit_ops: vec![
                GapPrelimEditScript {
                    op_type: GapAlignOpType::Del,
                    num: 2,
                },
                GapPrelimEditScript {
                    op_type: GapAlignOpType::Sub,
                    num: 3,
                },
            ],
            last_op: Some(GapAlignOpType::Sub),
        };

        let script = blast_prelim_edit_block_to_gap_edit_script(Some(&rev), Some(&fwd))
            .expect("edit script");

        assert_eq!(
            script.ops,
            vec![(GapAlignOpType::Sub, 7), (GapAlignOpType::Del, 2)]
        );
        let empty = GapPrelimEditBlock {
            edit_ops: Vec::new(),
            last_op: None,
        };
        assert!(blast_prelim_edit_block_to_gap_edit_script(Some(&empty), Some(&empty)).is_none());
    }

    #[test]
    fn oof_traceback_to_gap_edit_script_flows_frameshift_through_in_frame_gap() {
        let rev = GapPrelimEditBlock {
            edit_ops: vec![GapPrelimEditScript {
                op_type: GapAlignOpType::Ins1,
                num: 1,
            }],
            last_op: Some(GapAlignOpType::Ins1),
        };
        let fwd = GapPrelimEditBlock {
            edit_ops: vec![
                GapPrelimEditScript {
                    op_type: GapAlignOpType::Sub,
                    num: 2,
                },
                GapPrelimEditScript {
                    op_type: GapAlignOpType::Del,
                    num: 1,
                },
            ],
            last_op: Some(GapAlignOpType::Del),
        };
        let mut script = None;

        assert_eq!(
            s_blast_oof_traceback_to_gap_edit_script(
                Some(&rev),
                Some(&fwd),
                100,
                Some(&mut script),
            ),
            0
        );

        assert_eq!(
            script.expect("script").ops,
            vec![
                (GapAlignOpType::Sub, 1),
                (GapAlignOpType::Del, 1),
                (GapAlignOpType::Ins1, 1),
                (GapAlignOpType::Sub, 2),
            ]
        );
    }

    #[test]
    fn oof_traceback_to_gap_edit_script_does_not_duplicate_all_gap_right_tail() {
        let rev = GapPrelimEditBlock {
            edit_ops: vec![GapPrelimEditScript {
                op_type: GapAlignOpType::Ins1,
                num: 1,
            }],
            last_op: Some(GapAlignOpType::Ins1),
        };
        let fwd = GapPrelimEditBlock {
            edit_ops: vec![GapPrelimEditScript {
                op_type: GapAlignOpType::Del,
                num: 2,
            }],
            last_op: Some(GapAlignOpType::Del),
        };
        let mut script = None;

        assert_eq!(
            s_blast_oof_traceback_to_gap_edit_script(
                Some(&rev),
                Some(&fwd),
                100,
                Some(&mut script),
            ),
            0
        );

        assert_eq!(
            script.expect("script").ops,
            vec![(GapAlignOpType::Sub, 1), (GapAlignOpType::Del, 2)]
        );
    }

    #[test]
    fn oof_traceback_to_gap_edit_script_rejects_missing_required_inputs() {
        let rev = GapPrelimEditBlock {
            edit_ops: vec![GapPrelimEditScript {
                op_type: GapAlignOpType::Sub,
                num: 1,
            }],
            last_op: Some(GapAlignOpType::Sub),
        };
        let fwd = GapPrelimEditBlock {
            edit_ops: vec![GapPrelimEditScript {
                op_type: GapAlignOpType::Sub,
                num: 1,
            }],
            last_op: Some(GapAlignOpType::Sub),
        };
        let mut script = Some(GapEditScript {
            ops: vec![(GapAlignOpType::Sub, 1)],
        });

        assert_eq!(
            s_blast_oof_traceback_to_gap_edit_script(None, Some(&fwd), 3, Some(&mut script)),
            -1
        );
        assert_eq!(
            s_blast_oof_traceback_to_gap_edit_script(Some(&rev), None, 3, Some(&mut script)),
            -1
        );
        assert_eq!(
            s_blast_oof_traceback_to_gap_edit_script(Some(&rev), Some(&fwd), 3, None),
            -1
        );
        assert_eq!(
            script.expect("existing script").ops,
            vec![(GapAlignOpType::Sub, 1)]
        );
    }

    #[test]
    fn oof_traceback_to_gap_edit_script_empty_halves_clear_output_script() {
        let rev = GapPrelimEditBlock {
            edit_ops: Vec::new(),
            last_op: None,
        };
        let fwd = GapPrelimEditBlock {
            edit_ops: Vec::new(),
            last_op: None,
        };
        let mut script = Some(GapEditScript {
            ops: vec![(GapAlignOpType::Sub, 1)],
        });

        assert_eq!(
            s_blast_oof_traceback_to_gap_edit_script(Some(&rev), Some(&fwd), 0, Some(&mut script),),
            0
        );
        assert!(script.is_none());

        script = Some(GapEditScript {
            ops: vec![(GapAlignOpType::Sub, 1)],
        });
        assert_eq!(
            s_blast_oof_traceback_to_gap_edit_script(Some(&rev), Some(&fwd), -3, Some(&mut script),),
            0
        );
        assert!(script.is_none());
    }

    #[test]
    fn oof_traceback_to_gap_edit_script_zero_length_clears_nonempty_output_script() {
        let rev = GapPrelimEditBlock {
            edit_ops: vec![GapPrelimEditScript {
                op_type: GapAlignOpType::Sub,
                num: 2,
            }],
            last_op: Some(GapAlignOpType::Sub),
        };
        let fwd = GapPrelimEditBlock {
            edit_ops: vec![GapPrelimEditScript {
                op_type: GapAlignOpType::Del,
                num: 1,
            }],
            last_op: Some(GapAlignOpType::Del),
        };
        let mut script = Some(GapEditScript {
            ops: vec![(GapAlignOpType::Sub, 1)],
        });

        assert_eq!(
            s_blast_oof_traceback_to_gap_edit_script(Some(&rev), Some(&fwd), 0, Some(&mut script),),
            0
        );
        assert!(script.is_none());
    }

    #[test]
    fn oof_traceback_to_gap_edit_script_lengthens_substitution_after_split_frameshift() {
        let rev = GapPrelimEditBlock {
            edit_ops: vec![GapPrelimEditScript {
                op_type: GapAlignOpType::Ins1,
                num: 2,
            }],
            last_op: Some(GapAlignOpType::Ins1),
        };
        let fwd = GapPrelimEditBlock {
            edit_ops: vec![GapPrelimEditScript {
                op_type: GapAlignOpType::Sub,
                num: 3,
            }],
            last_op: Some(GapAlignOpType::Sub),
        };
        let mut script = None;

        assert_eq!(
            s_blast_oof_traceback_to_gap_edit_script(
                Some(&rev),
                Some(&fwd),
                100,
                Some(&mut script),
            ),
            0
        );

        assert_eq!(
            script.expect("script").ops,
            vec![
                (GapAlignOpType::Sub, 1),
                (GapAlignOpType::Ins1, 1),
                (GapAlignOpType::Ins1, 1),
                (GapAlignOpType::Sub, 3),
            ]
        );
    }

    #[test]
    fn oof_traceback_to_gap_edit_script_splits_deletion_frameshifts_before_substitution() {
        let rev = GapPrelimEditBlock {
            edit_ops: vec![GapPrelimEditScript {
                op_type: GapAlignOpType::Del1,
                num: 2,
            }],
            last_op: Some(GapAlignOpType::Del1),
        };
        let fwd = GapPrelimEditBlock {
            edit_ops: vec![GapPrelimEditScript {
                op_type: GapAlignOpType::Sub,
                num: 3,
            }],
            last_op: Some(GapAlignOpType::Sub),
        };
        let mut script = None;

        assert_eq!(
            s_blast_oof_traceback_to_gap_edit_script(
                Some(&rev),
                Some(&fwd),
                100,
                Some(&mut script),
            ),
            0
        );

        assert_eq!(
            script.expect("script").ops,
            vec![
                (GapAlignOpType::Sub, 1),
                (GapAlignOpType::Del1, 1),
                (GapAlignOpType::Del1, 1),
                (GapAlignOpType::Sub, 3),
            ]
        );
    }

    #[test]
    fn oof_traceback_to_gap_edit_script_cancels_opposing_frameshifts_at_join() {
        let rev = GapPrelimEditBlock {
            edit_ops: vec![GapPrelimEditScript {
                op_type: GapAlignOpType::Del2,
                num: 1,
            }],
            last_op: Some(GapAlignOpType::Del2),
        };
        let fwd = GapPrelimEditBlock {
            edit_ops: vec![
                GapPrelimEditScript {
                    op_type: GapAlignOpType::Sub,
                    num: 2,
                },
                GapPrelimEditScript {
                    op_type: GapAlignOpType::Ins2,
                    num: 1,
                },
            ],
            last_op: Some(GapAlignOpType::Ins2),
        };
        let mut script = None;

        assert_eq!(
            s_blast_oof_traceback_to_gap_edit_script(
                Some(&rev),
                Some(&fwd),
                100,
                Some(&mut script),
            ),
            0
        );

        assert_eq!(script.expect("script").ops, vec![(GapAlignOpType::Sub, 4)]);
    }

    #[test]
    fn oof_traceback_to_gap_edit_script_combines_same_direction_frameshifts_at_join() {
        let rev = GapPrelimEditBlock {
            edit_ops: vec![GapPrelimEditScript {
                op_type: GapAlignOpType::Del1,
                num: 1,
            }],
            last_op: Some(GapAlignOpType::Del1),
        };
        let fwd = GapPrelimEditBlock {
            edit_ops: vec![
                GapPrelimEditScript {
                    op_type: GapAlignOpType::Sub,
                    num: 2,
                },
                GapPrelimEditScript {
                    op_type: GapAlignOpType::Del2,
                    num: 1,
                },
            ],
            last_op: Some(GapAlignOpType::Del2),
        };
        let mut script = None;

        assert_eq!(
            s_blast_oof_traceback_to_gap_edit_script(
                Some(&rev),
                Some(&fwd),
                100,
                Some(&mut script),
            ),
            0
        );

        assert_eq!(
            script.expect("script").ops,
            vec![
                (GapAlignOpType::Sub, 1),
                (GapAlignOpType::Del, 1),
                (GapAlignOpType::Sub, 2)
            ]
        );
    }

    #[test]
    fn oof_traceback_to_gap_edit_script_combines_insertion_frameshifts_at_join() {
        let rev = GapPrelimEditBlock {
            edit_ops: vec![GapPrelimEditScript {
                op_type: GapAlignOpType::Ins1,
                num: 1,
            }],
            last_op: Some(GapAlignOpType::Ins1),
        };
        let fwd = GapPrelimEditBlock {
            edit_ops: vec![
                GapPrelimEditScript {
                    op_type: GapAlignOpType::Sub,
                    num: 2,
                },
                GapPrelimEditScript {
                    op_type: GapAlignOpType::Ins2,
                    num: 1,
                },
            ],
            last_op: Some(GapAlignOpType::Ins2),
        };
        let mut script = None;

        assert_eq!(
            s_blast_oof_traceback_to_gap_edit_script(
                Some(&rev),
                Some(&fwd),
                100,
                Some(&mut script),
            ),
            0
        );

        assert_eq!(
            script.expect("script").ops,
            vec![
                (GapAlignOpType::Sub, 1),
                (GapAlignOpType::Ins, 1),
                (GapAlignOpType::Sub, 2)
            ]
        );
    }

    #[test]
    fn oof_traceback_to_gap_edit_script_keeps_gap_between_joined_frameshift_and_substitution() {
        let rev = GapPrelimEditBlock {
            edit_ops: vec![GapPrelimEditScript {
                op_type: GapAlignOpType::Ins1,
                num: 1,
            }],
            last_op: Some(GapAlignOpType::Ins1),
        };
        let fwd = GapPrelimEditBlock {
            edit_ops: vec![
                GapPrelimEditScript {
                    op_type: GapAlignOpType::Sub,
                    num: 2,
                },
                GapPrelimEditScript {
                    op_type: GapAlignOpType::Del,
                    num: 1,
                },
                GapPrelimEditScript {
                    op_type: GapAlignOpType::Ins1,
                    num: 1,
                },
            ],
            last_op: Some(GapAlignOpType::Ins1),
        };
        let mut script = None;

        assert_eq!(
            s_blast_oof_traceback_to_gap_edit_script(
                Some(&rev),
                Some(&fwd),
                100,
                Some(&mut script),
            ),
            0
        );

        assert_eq!(
            script.expect("script").ops,
            vec![
                (GapAlignOpType::Sub, 1),
                (GapAlignOpType::Ins2, 1),
                (GapAlignOpType::Del, 1),
                (GapAlignOpType::Sub, 2)
            ]
        );
    }

    #[test]
    fn oof_traceback_to_gap_edit_script_preserves_gap_skipped_before_right_frameshift_merge() {
        let rev = GapPrelimEditBlock {
            edit_ops: vec![GapPrelimEditScript {
                op_type: GapAlignOpType::Del1,
                num: 1,
            }],
            last_op: Some(GapAlignOpType::Del1),
        };
        let fwd = GapPrelimEditBlock {
            edit_ops: vec![
                GapPrelimEditScript {
                    op_type: GapAlignOpType::Sub,
                    num: 2,
                },
                GapPrelimEditScript {
                    op_type: GapAlignOpType::Ins2,
                    num: 1,
                },
                GapPrelimEditScript {
                    op_type: GapAlignOpType::Del,
                    num: 1,
                },
            ],
            last_op: Some(GapAlignOpType::Del),
        };
        let mut script = None;

        assert_eq!(
            s_blast_oof_traceback_to_gap_edit_script(
                Some(&rev),
                Some(&fwd),
                100,
                Some(&mut script),
            ),
            0
        );

        assert_eq!(
            script.expect("script").ops,
            vec![
                (GapAlignOpType::Sub, 1),
                (GapAlignOpType::Del, 1),
                (GapAlignOpType::Ins1, 1),
                (GapAlignOpType::Sub, 3)
            ]
        );
    }

    #[test]
    fn oof_traceback_to_gap_edit_script_keeps_residual_right_frameshift_and_sub_tail() {
        let rev = GapPrelimEditBlock {
            edit_ops: vec![GapPrelimEditScript {
                op_type: GapAlignOpType::Del1,
                num: 1,
            }],
            last_op: Some(GapAlignOpType::Del1),
        };
        let fwd = GapPrelimEditBlock {
            edit_ops: vec![
                GapPrelimEditScript {
                    op_type: GapAlignOpType::Sub,
                    num: 2,
                },
                GapPrelimEditScript {
                    op_type: GapAlignOpType::Del2,
                    num: 2,
                },
            ],
            last_op: Some(GapAlignOpType::Del2),
        };
        let mut script = None;

        assert_eq!(
            s_blast_oof_traceback_to_gap_edit_script(
                Some(&rev),
                Some(&fwd),
                100,
                Some(&mut script),
            ),
            0
        );

        assert_eq!(
            script.expect("script").ops,
            vec![
                (GapAlignOpType::Sub, 1),
                (GapAlignOpType::Del, 1),
                (GapAlignOpType::Del2, 1),
                (GapAlignOpType::Sub, 3)
            ]
        );
    }

    #[test]
    fn oof_traceback_to_gap_edit_script_keeps_exact_nucleotide_boundary() {
        let rev = GapPrelimEditBlock {
            edit_ops: Vec::new(),
            last_op: None,
        };
        let fwd = GapPrelimEditBlock {
            edit_ops: vec![
                GapPrelimEditScript {
                    op_type: GapAlignOpType::Sub,
                    num: 2,
                },
                GapPrelimEditScript {
                    op_type: GapAlignOpType::Del,
                    num: 1,
                },
            ],
            last_op: Some(GapAlignOpType::Del),
        };
        let mut script = None;

        assert_eq!(
            s_blast_oof_traceback_to_gap_edit_script(Some(&rev), Some(&fwd), 6, Some(&mut script),),
            0
        );

        assert_eq!(
            script.expect("script").ops,
            vec![(GapAlignOpType::Del, 1), (GapAlignOpType::Sub, 2)]
        );
    }

    #[test]
    fn oof_traceback_to_gap_edit_script_truncates_partial_substitution_codon() {
        let rev = GapPrelimEditBlock {
            edit_ops: Vec::new(),
            last_op: None,
        };
        let fwd = GapPrelimEditBlock {
            edit_ops: vec![GapPrelimEditScript {
                op_type: GapAlignOpType::Sub,
                num: 3,
            }],
            last_op: Some(GapAlignOpType::Sub),
        };
        let mut script = None;

        assert_eq!(
            s_blast_oof_traceback_to_gap_edit_script(Some(&rev), Some(&fwd), 4, Some(&mut script),),
            0
        );

        assert_eq!(script.expect("script").ops, vec![(GapAlignOpType::Sub, 2)]);
    }

    #[test]
    fn oof_traceback_to_gap_edit_script_truncates_ordinary_insertion_by_codon_span() {
        let rev = GapPrelimEditBlock {
            edit_ops: Vec::new(),
            last_op: None,
        };
        let fwd = GapPrelimEditBlock {
            edit_ops: vec![
                GapPrelimEditScript {
                    op_type: GapAlignOpType::Sub,
                    num: 2,
                },
                GapPrelimEditScript {
                    op_type: GapAlignOpType::Ins,
                    num: 2,
                },
                GapPrelimEditScript {
                    op_type: GapAlignOpType::Sub,
                    num: 1,
                },
            ],
            last_op: Some(GapAlignOpType::Sub),
        };
        let mut script = None;

        assert_eq!(
            s_blast_oof_traceback_to_gap_edit_script(Some(&rev), Some(&fwd), 7, Some(&mut script),),
            0
        );

        assert_eq!(
            script.expect("script").ops,
            vec![(GapAlignOpType::Sub, 1), (GapAlignOpType::Ins, 2)]
        );
    }

    #[test]
    fn oof_traceback_to_gap_edit_script_preserves_mixed_frameshifts_and_ordinary_gaps() {
        let rev = GapPrelimEditBlock {
            edit_ops: Vec::new(),
            last_op: None,
        };
        let fwd = GapPrelimEditBlock {
            edit_ops: vec![
                GapPrelimEditScript {
                    op_type: GapAlignOpType::Sub,
                    num: 2,
                },
                GapPrelimEditScript {
                    op_type: GapAlignOpType::Ins1,
                    num: 1,
                },
                GapPrelimEditScript {
                    op_type: GapAlignOpType::Ins,
                    num: 1,
                },
                GapPrelimEditScript {
                    op_type: GapAlignOpType::Del,
                    num: 1,
                },
                GapPrelimEditScript {
                    op_type: GapAlignOpType::Ins1,
                    num: 1,
                },
                GapPrelimEditScript {
                    op_type: GapAlignOpType::Sub,
                    num: 2,
                },
            ],
            last_op: Some(GapAlignOpType::Sub),
        };
        let mut script = None;

        assert_eq!(
            s_blast_oof_traceback_to_gap_edit_script(
                Some(&rev),
                Some(&fwd),
                100,
                Some(&mut script),
            ),
            0
        );

        assert_eq!(
            script.expect("script").ops,
            vec![
                (GapAlignOpType::Sub, 2),
                (GapAlignOpType::Ins1, 1),
                (GapAlignOpType::Del, 1),
                (GapAlignOpType::Ins, 1),
                (GapAlignOpType::Ins1, 1),
                (GapAlignOpType::Sub, 3)
            ]
        );
    }

    #[test]
    fn oof_traceback_to_gap_edit_script_lengthens_post_truncation_substitution_after_frameshift() {
        let rev = GapPrelimEditBlock {
            edit_ops: Vec::new(),
            last_op: None,
        };
        let fwd = GapPrelimEditBlock {
            edit_ops: vec![
                GapPrelimEditScript {
                    op_type: GapAlignOpType::Sub,
                    num: 2,
                },
                GapPrelimEditScript {
                    op_type: GapAlignOpType::Ins1,
                    num: 1,
                },
                GapPrelimEditScript {
                    op_type: GapAlignOpType::Sub,
                    num: 1,
                },
            ],
            last_op: Some(GapAlignOpType::Sub),
        };
        let mut script = None;

        assert_eq!(
            s_blast_oof_traceback_to_gap_edit_script(Some(&rev), Some(&fwd), 8, Some(&mut script),),
            0
        );

        assert_eq!(
            script.expect("script").ops,
            vec![
                (GapAlignOpType::Sub, 1),
                (GapAlignOpType::Ins1, 1),
                (GapAlignOpType::Sub, 2)
            ]
        );
    }

    #[test]
    fn oof_traceback_to_gap_edit_script_truncates_inside_frameshift_span() {
        let rev = GapPrelimEditBlock {
            edit_ops: Vec::new(),
            last_op: None,
        };
        let fwd = GapPrelimEditBlock {
            edit_ops: vec![
                GapPrelimEditScript {
                    op_type: GapAlignOpType::Sub,
                    num: 2,
                },
                GapPrelimEditScript {
                    op_type: GapAlignOpType::Ins1,
                    num: 3,
                },
            ],
            last_op: Some(GapAlignOpType::Ins1),
        };
        let mut script = None;

        assert_eq!(
            s_blast_oof_traceback_to_gap_edit_script(Some(&rev), Some(&fwd), 5, Some(&mut script),),
            0
        );

        assert_eq!(
            script.expect("script").ops,
            vec![(GapAlignOpType::Ins1, 1), (GapAlignOpType::Ins1, 1)]
        );
    }

    #[test]
    fn jumper_lifecycle_helpers_allocate_copy_and_free_owned_buffers() {
        assert!(jumper_prelim_edit_block_new(0).is_none());
        let mut prelim = jumper_prelim_edit_block_new(1).expect("prelim block");
        assert_eq!(
            jumper_prelim_edit_block_add(&mut prelim, JUMPER_MISMATCH),
            0
        );
        assert_eq!(
            jumper_prelim_edit_block_add(&mut prelim, JUMPER_INSERTION),
            0
        );
        assert_eq!(jumper_prelim_edit_block_add(&mut prelim, 9), 0);
        assert_eq!(jumper_prelim_edit_block_add(&mut prelim, 4), 0);
        assert_eq!(prelim.edit_ops, vec![JUMPER_MISMATCH, JUMPER_INSERTION, 13]);
        assert!(jumper_prelim_edit_block_free(Some(prelim)).is_none());

        let mut table = vec![0u32; 256];
        s_create_table(&mut table);
        assert_eq!(table[0], 0);
        assert_eq!(table[0b01_10_11_00], 0x0003_0201);

        let gap_align = jumper_gap_align_new(4).expect("jumper gap align");
        assert!(gap_align.left_prelim_block.is_some());
        assert!(gap_align.right_prelim_block.is_some());
        assert_eq!(gap_align.table.len(), 256);
        assert!(jumper_gap_align_new(0).is_none());
        assert!(jumper_gap_align_free(Some(gap_align)).is_none());

        let mut edits = jumper_edits_block_new(2).expect("edits");
        edits.edits.push(JumperEdit {
            query_pos: 3,
            query_base: 1,
            subject_base: 2,
        });
        let dup = jumper_edits_block_dup(Some(&edits)).expect("dup");
        assert_eq!(dup, edits);
        assert!(jumper_edits_block_dup(None).is_none());
        assert!(jumper_edits_block_free(Some(edits)).is_none());

        let overhangs = SequenceOverhangs {
            left: Some(vec![1, 2]),
            right: Some(vec![3]),
        };
        assert!(sequence_overhangs_free(Some(overhangs)).is_none());
    }

    #[test]
    fn jumper_reset_convert_and_combine_match_c_edit_rules() {
        let mut left = jumper_prelim_edit_block_new(4).expect("left block");
        let mut right = jumper_prelim_edit_block_new(4).expect("right block");
        left.edit_ops.extend([JUMPER_INSERTION, 3, JUMPER_MISMATCH]);
        right.edit_ops.extend([2, JUMPER_DELETION]);

        let script =
            jumper_prelim_edit_block_to_gap_edit_script(&left, &right).expect("gap script");
        assert_eq!(
            script.ops,
            vec![
                (GapAlignOpType::Sub, 4),
                (GapAlignOpType::Ins, 1),
                (GapAlignOpType::Sub, 2),
                (GapAlignOpType::Del, 1),
            ]
        );

        s_reset_jumper_prelim_edit_blocks(Some(&mut left), Some(&mut right));
        assert!(left.edit_ops.is_empty());
        assert!(right.edit_ops.is_empty());
        s_reset_jumper_prelim_edit_blocks(Some(&mut left), None);

        right
            .edit_ops
            .extend([2, JUMPER_MISMATCH, JUMPER_INSERTION, JUMPER_DELETION, 3]);
        let mut query_pos = 10;
        let mut subject_pos = 20;
        assert_eq!(
            s_get_seq_positions(Some(&right), 4, &mut query_pos, &mut subject_pos),
            0
        );
        assert_eq!((query_pos, subject_pos), (14, 24));
        assert_eq!(
            s_get_seq_positions(Some(&right), 6, &mut query_pos, &mut subject_pos),
            -1
        );
        assert_eq!(
            s_get_seq_positions(Some(&right), -1, &mut query_pos, &mut subject_pos),
            -1
        );
        assert_eq!(
            s_get_seq_positions(None, 0, &mut query_pos, &mut subject_pos),
            -1
        );

        query_pos = 0;
        subject_pos = 0;
        assert_eq!(
            s_get_seq_positions(
                Some(&right),
                right.edit_ops.len() as i32,
                &mut query_pos,
                &mut subject_pos,
            ),
            0
        );
        assert_eq!((query_pos, subject_pos), (7, 7));

        assert_eq!(s_compute_extension_score(&right, 2, -3, -5, -1), -5);
        let scored = JumperPrelimEditBlock {
            edit_ops: vec![
                2,
                JUMPER_INSERTION,
                JUMPER_INSERTION,
                JUMPER_DELETION,
                JUMPER_MISMATCH,
            ],
        };
        assert_eq!(s_compute_extension_score(&scored, 2, -3, -5, -1), -12);

        let mut base = Some(JumperEditsBlock {
            edits: vec![JumperEdit {
                query_pos: 1,
                query_base: 2,
                subject_base: 3,
            }],
        });
        let mut append = Some(JumperEditsBlock {
            edits: vec![JumperEdit {
                query_pos: 4,
                query_base: 5,
                subject_base: 6,
            }],
        });
        let combined = jumper_edits_block_combine(&mut base, &mut append).expect("combined");
        assert_eq!(combined.edits.len(), 2);
        assert!(append.is_none());
        assert_eq!(base.as_ref().unwrap().edits, combined.edits);

        let mut missing_base = None;
        let mut append_only = Some(JumperEditsBlock::default());
        assert!(jumper_edits_block_combine(&mut missing_base, &mut append_only).is_none());

        let it = SubjectIndexIterator {
            subject_index: Some(SubjectIndex {
                num_lookups: 3,
                width: 100,
                ..Default::default()
            }),
            to: 250,
            lookup_index: 1,
            ..Default::default()
        };
        assert!(subject_index_iterator_free(Some(it)).is_none());
    }

    #[test]
    fn subject_index_new_and_iterators_find_packed_words() {
        let subject_bases = [0u8, 1, 2, 3, 0, 1, 2, 3, 0];
        let packed = crate::encoding::pack_ncbi2na_bases(&subject_bases);
        let index = subject_index_new(&packed, subject_bases.len() as i32, 4, 2).expect("index");
        assert_eq!(index.num_lookups, 3);
        let word_ac = (0u32 << 2) | 1;

        let mut it = subject_index_iterator_new(&index, word_ac, 1, 7).expect("iterator");
        assert_eq!(subject_index_iterator_next(&mut it), 4);
        assert_eq!(subject_index_iterator_next(&mut it), -1);

        let mut rev = subject_index_iterator_new(&index, word_ac, 99, 0).expect("iterator");
        assert_eq!(subject_index_iterator_prev(&mut rev), 4);
        assert_eq!(subject_index_iterator_prev(&mut rev), 0);
        assert_eq!(subject_index_iterator_prev(&mut rev), -1);
    }

    #[test]
    fn subject_index_rejects_invalid_shapes_and_empty_iterators_stop() {
        let subject_bases = [0u8, 1, 2, 3];
        let packed = crate::encoding::pack_ncbi2na_bases(&subject_bases);

        assert!(subject_index_new(&packed, -1, 4, 2).is_none());
        assert!(subject_index_new(&packed, subject_bases.len() as i32, 0, 2).is_none());
        assert!(subject_index_new(&packed, subject_bases.len() as i32, 4, 0).is_none());
        assert!(subject_index_new(&[], subject_bases.len() as i32, 4, 2).is_none());
        assert!(subject_word_at(&packed, i32::MAX, i32::MAX, 4).is_none());

        let short = subject_index_new(&packed, 1, 4, 2).expect("short index");
        assert!(short.positions_by_word.is_empty());
        let mut missing_word = subject_index_iterator_new(&short, 0b0001, 0, 10).unwrap();
        assert_eq!(subject_index_iterator_next(&mut missing_word), -1);
        assert_eq!(subject_index_iterator_prev(&mut missing_word), -1);

        let index = subject_index_new(&packed, subject_bases.len() as i32, 4, 2).expect("index");
        let word_ac = (0u32 << 2) | 1;
        let mut out_of_bounds_next = subject_index_iterator_new(&index, word_ac, 99, 100).unwrap();
        assert_eq!(subject_index_iterator_next(&mut out_of_bounds_next), -1);

        let mut bounded_prev = subject_index_iterator_new(&index, word_ac, 99, 3).unwrap();
        assert_eq!(subject_index_iterator_prev(&mut bounded_prev), -1);

        let lookup_window = subject_index_iterator_new(&index, word_ac, 7, 8).unwrap();
        assert_eq!(lookup_window.lookup_index, 1);

        let mut negative_next = SubjectIndexIterator {
            positions: vec![0],
            word_index: -1,
            to: 10,
            ..Default::default()
        };
        assert_eq!(subject_index_iterator_next(&mut negative_next), -1);

        let mut negative_prev = SubjectIndexIterator {
            positions: vec![0],
            word_index: -1,
            to: 0,
            ..Default::default()
        };
        assert_eq!(subject_index_iterator_prev(&mut negative_prev), -1);
    }

    fn mapper_saved_hsp(
        query_offset: i32,
        query_end: i32,
        subject_offset: i32,
        subject_end: i32,
    ) -> crate::hspstream::Hsp {
        let span = query_end.saturating_sub(query_offset);
        crate::hspstream::Hsp {
            score: span,
            num_ident: span,
            bit_score: 0.0,
            evalue: 0.0,
            query_offset,
            query_end,
            query_gapped_start: query_offset,
            subject_offset,
            subject_end,
            subject_gapped_start: subject_offset,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        }
    }

    #[test]
    fn blast_na_extend_jumper_small_word_rescue_finds_right_indexed_overhang() {
        let query: Vec<u8> = (0..24).map(|i| (i % 4) as u8).collect();
        let mut subject_bases = vec![3; 150];
        subject_bases[39] = 0;
        subject_bases[40..52].copy_from_slice(&query[12..24]);
        let subject = crate::encoding::pack_ncbi2na_bases(&subject_bases);
        let index =
            subject_index_new(&subject, subject_bases.len() as i32, 10000, 4).expect("index");
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let saved_hsp = mapper_saved_hsp(4, 12, 8, 20);
        let mut list = crate::hspstream::HspList::new(1);

        assert!(
            blast_na_extend_jumper_small_word_rescue(
                &saved_hsp,
                &index,
                0,
                &query,
                0,
                query.len() as i32,
                &subject,
                subject_bases.len() as i32,
                0,
                &score_params,
                60,
                &mut list,
            ) >= 1
        );
        assert!(list.hsps.iter().any(|hsp| {
            hsp.query_offset == 12
                && hsp.query_end == 24
                && hsp.subject_offset == 40
                && hsp.subject_end == 52
        }));
    }

    #[test]
    fn blast_na_extend_jumper_small_word_rescue_trims_ambiguous_right_overhang_edges() {
        let mut query: Vec<u8> = (0..24).map(|i| (i % 4) as u8).collect();
        query[12] = 14;
        query[23] = 14;
        let mut subject_bases = vec![3; 150];
        subject_bases[39] = 0;
        subject_bases[40..52].copy_from_slice(&(12..24).map(|i| (i % 4) as u8).collect::<Vec<_>>());
        let subject = crate::encoding::pack_ncbi2na_bases(&subject_bases);
        let index =
            subject_index_new(&subject, subject_bases.len() as i32, 10000, 4).expect("index");
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let saved_hsp = mapper_saved_hsp(4, 12, 8, 20);
        let mut list = crate::hspstream::HspList::new(1);

        assert!(
            blast_na_extend_jumper_small_word_rescue(
                &saved_hsp,
                &index,
                0,
                &query,
                0,
                query.len() as i32,
                &subject,
                subject_bases.len() as i32,
                0,
                &score_params,
                60,
                &mut list,
            ) >= 1
        );
        assert!(list.hsps.iter().any(|hsp| {
            hsp.query_offset == 13
                && hsp.query_end == 23
                && hsp.subject_offset == 41
                && hsp.subject_end == 51
        }));
    }

    #[test]
    fn blast_na_extend_jumper_small_word_rescue_finds_left_indexed_overhang() {
        let query: Vec<u8> = (0..24).map(|i| (i % 4) as u8).collect();
        let mut subject_bases = vec![3; 160];
        subject_bases[72..96].copy_from_slice(&query);
        let subject = crate::encoding::pack_ncbi2na_bases(&subject_bases);
        let index =
            subject_index_new(&subject, subject_bases.len() as i32, 10000, 4).expect("index");
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let saved_hsp = mapper_saved_hsp(12, 20, 120, 128);
        let mut list = crate::hspstream::HspList::new(1);

        assert!(
            blast_na_extend_jumper_small_word_rescue(
                &saved_hsp,
                &index,
                0,
                &query,
                0,
                query.len() as i32,
                &subject,
                subject_bases.len() as i32,
                0,
                &score_params,
                60,
                &mut list,
            ) >= 1
        );
        assert!(list.hsps.iter().any(|hsp| {
            hsp.query_offset == 0
                && hsp.query_end == 24
                && hsp.subject_offset == 72
                && hsp.subject_end == 96
        }));
    }

    #[test]
    fn blast_na_extend_jumper_small_word_rescue_trims_ambiguous_left_overhang_edges() {
        let mut query: Vec<u8> = (0..24).map(|i| (i % 4) as u8).collect();
        query[0] = 14;
        query[23] = 14;
        let mut subject_bases = vec![3; 160];
        subject_bases[72..96].copy_from_slice(&(0..24).map(|i| (i % 4) as u8).collect::<Vec<_>>());
        let subject = crate::encoding::pack_ncbi2na_bases(&subject_bases);
        let index =
            subject_index_new(&subject, subject_bases.len() as i32, 10000, 4).expect("index");
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let saved_hsp = mapper_saved_hsp(12, 20, 120, 128);
        let mut list = crate::hspstream::HspList::new(1);

        assert!(
            blast_na_extend_jumper_small_word_rescue(
                &saved_hsp,
                &index,
                0,
                &query,
                0,
                query.len() as i32,
                &subject,
                subject_bases.len() as i32,
                0,
                &score_params,
                60,
                &mut list,
            ) >= 1
        );
        assert!(list.hsps.iter().any(|hsp| {
            hsp.query_offset == 1
                && hsp.query_end == 23
                && hsp.subject_offset == 73
                && hsp.subject_end == 95
        }));
    }

    #[test]
    fn blast_na_extend_jumper_small_word_rescue_obeys_subject_overhang_thresholds() {
        let query: Vec<u8> = (0..24).map(|i| (i % 4) as u8).collect();
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);

        let mut short_right_subject_bases = vec![3; 115];
        short_right_subject_bases[40..52].copy_from_slice(&query[12..24]);
        let short_right_subject = crate::encoding::pack_ncbi2na_bases(&short_right_subject_bases);
        let short_right_index = subject_index_new(
            &short_right_subject,
            short_right_subject_bases.len() as i32,
            10000,
            4,
        )
        .expect("index");
        let saved_right_hsp = mapper_saved_hsp(4, 12, 8, 20);
        let mut right_list = crate::hspstream::HspList::new(1);

        assert_eq!(
            blast_na_extend_jumper_small_word_rescue(
                &saved_right_hsp,
                &short_right_index,
                0,
                &query,
                0,
                query.len() as i32,
                &short_right_subject,
                short_right_subject_bases.len() as i32,
                0,
                &score_params,
                60,
                &mut right_list,
            ),
            0
        );
        assert!(right_list.hsps.is_empty());

        let mut short_left_subject_bases = vec![3; 160];
        short_left_subject_bases[72..96].copy_from_slice(&query);
        let short_left_subject = crate::encoding::pack_ncbi2na_bases(&short_left_subject_bases);
        let short_left_index = subject_index_new(
            &short_left_subject,
            short_left_subject_bases.len() as i32,
            10000,
            4,
        )
        .expect("index");
        let saved_left_hsp = mapper_saved_hsp(12, 20, 99, 107);
        let mut left_list = crate::hspstream::HspList::new(1);

        assert_eq!(
            blast_na_extend_jumper_small_word_rescue(
                &saved_left_hsp,
                &short_left_index,
                0,
                &query,
                0,
                query.len() as i32,
                &short_left_subject,
                short_left_subject_bases.len() as i32,
                0,
                &score_params,
                60,
                &mut left_list,
            ),
            0
        );
        assert!(left_list.hsps.is_empty());
    }

    #[test]
    fn blast_na_extend_jumper_small_word_rescue_saves_direct_short_overhangs() {
        let query: Vec<u8> = (0..24).map(|i| (i % 4) as u8).collect();
        let mut subject_bases = vec![3; 180];
        subject_bases[21..26].copy_from_slice(&query[13..18]);
        subject_bases[112..119].copy_from_slice(&query[0..7]);
        let subject = crate::encoding::pack_ncbi2na_bases(&subject_bases);
        let index =
            subject_index_new(&subject, subject_bases.len() as i32, 10000, 4).expect("index");
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let mut score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        score_params.penalty = -10;

        let saved_right_hsp = mapper_saved_hsp(4, 12, 8, 20);
        let mut right_list = crate::hspstream::HspList::new(4);
        assert!(
            blast_na_extend_jumper_small_word_rescue(
                &saved_right_hsp,
                &index,
                0,
                &query,
                0,
                18,
                &subject,
                subject_bases.len() as i32,
                0,
                &score_params,
                60,
                &mut right_list,
            ) >= 1
        );
        assert!(right_list.hsps.iter().any(|hsp| {
            hsp.query_offset == 13
                && hsp.query_end == 18
                && hsp.subject_offset == 21
                && hsp.subject_end == 26
        }));

        let saved_left_hsp = mapper_saved_hsp(8, 16, 120, 128);
        let mut left_list = crate::hspstream::HspList::new(4);
        assert!(
            blast_na_extend_jumper_small_word_rescue(
                &saved_left_hsp,
                &index,
                0,
                &query,
                0,
                query.len() as i32,
                &subject,
                subject_bases.len() as i32,
                0,
                &score_params,
                60,
                &mut left_list,
            ) >= 1
        );
        assert!(left_list.hsps.iter().any(|hsp| {
            hsp.query_offset == 0
                && hsp.query_end == 7
                && hsp.subject_offset == 112
                && hsp.subject_end == 119
        }));

        let mut boundary_hsp = mapper_saved_hsp(4, 12, 8, 20);
        boundary_hsp.query_end = 18;
        let mut boundary_list = crate::hspstream::HspList::new(4);
        assert_eq!(
            blast_na_extend_jumper_small_word_rescue(
                &boundary_hsp,
                &index,
                0,
                &query,
                0,
                18,
                &subject,
                subject_bases.len() as i32,
                0,
                &score_params,
                60,
                &mut boundary_list,
            ),
            0
        );
    }

    #[test]
    fn blast_na_extend_jumper_small_word_rescue_noops_when_index_has_no_match() {
        let query: Vec<u8> = (0..24).map(|i| (i % 4) as u8).collect();
        let subject_bases = vec![3; 180];
        let subject = crate::encoding::pack_ncbi2na_bases(&subject_bases);
        let index =
            subject_index_new(&subject, subject_bases.len() as i32, 10000, 4).expect("index");
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let saved_hsp = mapper_saved_hsp(4, 12, 8, 20);
        let mut list = crate::hspstream::HspList::new(1);

        assert_eq!(
            blast_na_extend_jumper_small_word_rescue(
                &saved_hsp,
                &index,
                0,
                &query,
                0,
                query.len() as i32,
                &subject,
                subject_bases.len() as i32,
                0,
                &score_params,
                60,
                &mut list,
            ),
            0
        );
        assert!(list.hsps.is_empty());

        let mut huge_intron_list = crate::hspstream::HspList::new(1);
        assert_eq!(
            blast_na_extend_jumper_small_word_rescue(
                &saved_hsp,
                &index,
                0,
                &query,
                0,
                query.len() as i32,
                &subject,
                subject_bases.len() as i32,
                0,
                &score_params,
                i32::MAX,
                &mut huge_intron_list,
            ),
            0
        );
        assert_eq!(
            blast_na_extend_jumper_small_word_rescue(
                &saved_hsp,
                &index,
                0,
                &query,
                0,
                query.len() as i32,
                &subject,
                subject_bases.len() as i32,
                0,
                &score_params,
                i32::MIN,
                &mut huge_intron_list,
            ),
            0
        );
        assert!(huge_intron_list.hsps.is_empty());

        let mut malformed_list = crate::hspstream::HspList::new(1);
        for bad_hsp in [
            mapper_saved_hsp(-1, 12, 8, 20),
            mapper_saved_hsp(12, 4, 8, 20),
            mapper_saved_hsp(4, 99, 8, 20),
            mapper_saved_hsp(4, 12, -1, 20),
            mapper_saved_hsp(4, 12, 20, 8),
            mapper_saved_hsp(4, 12, 8, 999),
        ] {
            assert_eq!(
                blast_na_extend_jumper_small_word_rescue(
                    &bad_hsp,
                    &index,
                    0,
                    &query,
                    0,
                    query.len() as i32,
                    &subject,
                    subject_bases.len() as i32,
                    0,
                    &score_params,
                    60,
                    &mut malformed_list,
                ),
                0
            );
        }
        assert_eq!(
            blast_na_extend_jumper_small_word_rescue(
                &saved_hsp,
                &index,
                0,
                &query[..8],
                0,
                query.len() as i32,
                &subject,
                subject_bases.len() as i32,
                0,
                &score_params,
                60,
                &mut malformed_list,
            ),
            0
        );
        assert_eq!(
            blast_na_extend_jumper_small_word_rescue(
                &saved_hsp,
                &index,
                0,
                &query,
                0,
                query.len() as i32,
                &[],
                subject_bases.len() as i32,
                0,
                &score_params,
                60,
                &mut malformed_list,
            ),
            0
        );
        assert!(malformed_list.hsps.is_empty());
    }

    #[test]
    fn jumper_base_matches_or_ambiguous_handles_bounds_and_wildcards() {
        let subject = crate::encoding::pack_ncbi2na_bases(&[0u8, 1, 2, 3]);

        assert!(jumper_base_matches_or_ambiguous(1, &subject, 1));
        assert!(!jumper_base_matches_or_ambiguous(2, &subject, 1));
        assert!(!jumper_base_matches_or_ambiguous(1, &subject, -1));
        assert!(!jumper_base_matches_or_ambiguous(1, &[], 0));
        assert!(jumper_base_matches_or_ambiguous(0xff, &[], -99));
    }

    #[test]
    fn jumper_shift_gaps_right_collapses_opposite_indels_into_match() {
        let query = [0u8, 1, 2, 3];
        let subject = crate::encoding::pack_ncbi2na_bases(&query);
        let mut script = JumperPrelimEditBlock {
            edit_ops: vec![JUMPER_INSERTION, JUMPER_DELETION],
        };
        let mut score = 10;
        let mut num_identical = 0;

        assert_eq!(
            s_shift_gaps_right(
                &mut script,
                &query,
                &subject,
                0,
                0,
                query.len() as i32,
                query.len() as i32,
                &mut score,
                -3,
                &mut num_identical,
            ),
            0
        );

        assert_eq!(script.edit_ops, vec![1]);
        assert_eq!(score, 14);
        assert_eq!(num_identical, 1);
    }

    #[test]
    fn jumper_shift_gaps_merges_halves_and_trims_flanking_gap() {
        let query = [0u8, 1, 2, 3];
        let subject = crate::encoding::pack_ncbi2na_bases(&query);
        let mut jumper = JumperGapAlign {
            left_prelim_block: Some(JumperPrelimEditBlock { edit_ops: vec![2] }),
            right_prelim_block: Some(JumperPrelimEditBlock {
                edit_ops: vec![3, JUMPER_DELETION],
            }),
            table: Vec::new(),
        };
        let mut query_stop = 4;
        let mut subject_stop = 5;
        let mut score = 20;
        let mut num_identical = 5;

        assert_eq!(
            s_shift_gaps(
                &mut jumper,
                &query,
                &subject,
                0,
                0,
                &mut query_stop,
                &mut subject_stop,
                &mut score,
                query.len() as i32,
                query.len() as i32,
                -3,
                &mut num_identical,
            ),
            0
        );

        assert_eq!(
            jumper.left_prelim_block.as_ref().unwrap().edit_ops,
            Vec::<JumperOpType>::new()
        );
        assert_eq!(
            jumper.right_prelim_block.as_ref().unwrap().edit_ops,
            vec![5]
        );
        assert_eq!((query_stop, subject_stop), (4, 4));
        assert_eq!(score, 23);
        assert_eq!(num_identical, 5);
    }

    #[test]
    fn jumper_shift_gaps_rejects_missing_prelim_blocks() {
        let query = [0u8, 1, 2, 3];
        let subject = crate::encoding::pack_ncbi2na_bases(&query);
        let mut query_stop = 4;
        let mut subject_stop = 4;
        let mut score = 8;
        let mut num_identical = 4;

        let mut missing_left = JumperGapAlign {
            left_prelim_block: None,
            right_prelim_block: Some(JumperPrelimEditBlock { edit_ops: vec![4] }),
            table: Vec::new(),
        };
        assert_eq!(
            s_shift_gaps(
                &mut missing_left,
                &query,
                &subject,
                0,
                0,
                &mut query_stop,
                &mut subject_stop,
                &mut score,
                query.len() as i32,
                query.len() as i32,
                -3,
                &mut num_identical,
            ),
            -1
        );
        assert_eq!(
            (query_stop, subject_stop, score, num_identical),
            (4, 4, 8, 4)
        );

        let mut missing_right = JumperGapAlign {
            left_prelim_block: Some(JumperPrelimEditBlock { edit_ops: vec![4] }),
            right_prelim_block: None,
            table: Vec::new(),
        };
        assert_eq!(
            s_shift_gaps(
                &mut missing_right,
                &query,
                &subject,
                0,
                0,
                &mut query_stop,
                &mut subject_stop,
                &mut score,
                query.len() as i32,
                query.len() as i32,
                -3,
                &mut num_identical,
            ),
            -1
        );
        assert_eq!(
            (query_stop, subject_stop, score, num_identical),
            (4, 4, 8, 4)
        );
    }

    #[test]
    fn jumper_packed_subject_base_and_match_reject_bad_bounds() {
        let query = [0u8, 1, 2, 3];
        let subject = crate::encoding::pack_ncbi2na_bases(&query);

        assert_eq!(jumper_packed_subject_base(&subject, 2, 4), Some(2));
        assert_eq!(jumper_packed_subject_base(&subject, -1, 4), None);
        assert_eq!(jumper_packed_subject_base(&subject, 4, 4), None);
        assert_eq!(jumper_packed_subject_base(&[], 0, 4), None);

        assert!(jumper_bases_match(&query, &subject, 1, 1, 4, 4));
        assert!(!jumper_bases_match(&query, &subject, -1, 1, 4, 4));
        assert!(!jumper_bases_match(&query, &subject, 4, 1, 4, 4));
        assert!(!jumper_bases_match(&query, &[], 1, 1, 4, 4));
    }

    #[test]
    fn jumper_trim_extension_tracks_c_pointer_arithmetic() {
        let mut script = JumperPrelimEditBlock {
            edit_ops: vec![1, JUMPER_MISMATCH, 2],
        };
        let mut cp = 10;
        let mut cq = 20;
        let mut num_identical = 3;

        s_trim_extension(&mut script, 3, &mut cp, &mut cq, &mut num_identical, true);

        assert_eq!(script.edit_ops, vec![1]);
        assert_eq!((cp, cq), (7, 17));
        assert_eq!(num_identical, 1);

        let mut left_script = JumperPrelimEditBlock {
            edit_ops: vec![1, JUMPER_INSERTION],
        };
        s_trim_extension(
            &mut left_script,
            2,
            &mut cp,
            &mut cq,
            &mut num_identical,
            false,
        );
        assert_eq!(left_script.edit_ops, vec![1]);
        assert_eq!((cp, cq), (8, 17));
    }

    #[test]
    fn jumper_trim_extension_handles_noop_and_single_gap_edges() {
        let mut script = JumperPrelimEditBlock {
            edit_ops: vec![JUMPER_DELETION],
        };
        let mut cp = 3;
        let mut cq = 4;
        let mut num_identical = 0;

        s_trim_extension(&mut script, 0, &mut cp, &mut cq, &mut num_identical, true);
        assert_eq!(script.edit_ops, vec![JUMPER_DELETION]);
        assert_eq!((cp, cq, num_identical), (3, 4, 0));

        s_trim_extension(&mut script, 2, &mut cp, &mut cq, &mut num_identical, true);
        assert!(script.edit_ops.is_empty());
        assert_eq!((cp, cq, num_identical), (3, 4, 0));

        let mut empty = JumperPrelimEditBlock::default();
        s_trim_extension(&mut empty, 2, &mut cp, &mut cq, &mut num_identical, false);
        assert!(empty.edit_ops.is_empty());
        assert_eq!((cp, cq, num_identical), (3, 4, 0));
    }

    #[test]
    fn jumper_extend_right_matches_c_seed_and_trace_shape() {
        let query = [0u8, 9, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0];
        let subject = [0u8, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0];
        let mut script = gap_prelim_edit_block_new();

        let (score, q_len, s_len) =
            jumper_extend_right(&query, &subject, 2, -3, -5, -1, 3, 8, &mut script, false);

        assert_eq!(score, 18);
        assert_eq!((q_len, s_len), (14, 13));
        assert_eq!(
            script
                .edit_ops
                .iter()
                .map(|op| (op.op_type, op.num))
                .collect::<Vec<_>>(),
            vec![
                (GapAlignOpType::Sub, 1),
                (GapAlignOpType::Ins, 1),
                (GapAlignOpType::Sub, 12),
            ]
        );
    }

    #[test]
    fn jumper_extend_right_with_traceback_records_raw_c_ops() {
        let query = [0u8, 9, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0];
        let subject = [0u8, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0];
        let mut script = JumperPrelimEditBlock {
            edit_ops: Vec::new(),
        };
        let mut num_identical = 0;
        let mut ungapped_ext_len = -1;

        let (score, q_len, s_len) = jumper_extend_right_with_traceback(
            &query,
            &subject,
            2,
            -3,
            -5,
            -1,
            3,
            8,
            &mut script,
            &mut num_identical,
            false,
            &mut ungapped_ext_len,
        );

        assert_eq!(score, 20);
        assert_eq!((q_len, s_len), (14, 13));
        assert_eq!(script.edit_ops, vec![1, JUMPER_INSERTION, 10, 2]);
        assert_eq!(num_identical, 13);
        assert_eq!(ungapped_ext_len, 0);
    }

    #[test]
    fn jumper_extend_left_matches_c_reverse_walk() {
        let query = [0u8, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 9];
        let subject = [0u8, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
        let mut script = gap_prelim_edit_block_new();

        let (score, q_len, s_len) =
            jumper_extend_left(&query, &subject, 12, 11, 2, -3, -5, -1, 3, 8, &mut script);

        assert_eq!(score, 18);
        assert_eq!((q_len, s_len), (13, 12));
        assert_eq!(
            script
                .edit_ops
                .iter()
                .map(|op| (op.op_type, op.num))
                .collect::<Vec<_>>(),
            vec![(GapAlignOpType::Ins, 1), (GapAlignOpType::Sub, 12)]
        );
    }

    #[test]
    fn jumper_extend_right_compressed_tracks_best_prefix_and_ungapped_len() {
        let query = [0u8, 9, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0];
        let subject_bases = [0u8, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0];
        let subject = crate::encoding::pack_ncbi2na_bases(&subject_bases);
        let mut num_identical = 0;
        let mut ungapped_ext_len = -1;

        let (score, q_len, s_len) = jumper_extend_right_compressed(
            &query,
            &subject,
            subject_bases.len() as i32,
            2,
            -3,
            3,
            8,
            &mut num_identical,
            &mut ungapped_ext_len,
        );

        assert_eq!(score, 24);
        assert_eq!((q_len, s_len), (14, 13));
        assert_eq!(num_identical, 12);
        assert_eq!(ungapped_ext_len, 0);
    }

    #[test]
    fn jumper_extend_left_compressed_walks_packed_subject_backwards() {
        let query = [0u8, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 9];
        let subject_bases = [0u8, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
        let subject = crate::encoding::pack_ncbi2na_bases(&subject_bases);
        let mut num_identical = 0;

        let (score, q_len, s_len) = jumper_extend_left_compressed(
            &query,
            &subject,
            12,
            11,
            2,
            -3,
            3,
            8,
            &mut num_identical,
        );

        assert_eq!(score, 24);
        assert_eq!((q_len, s_len), (13, 12));
        assert_eq!(num_identical, 12);
    }

    #[test]
    fn jumper_extend_right_compressed_with_traceback_records_raw_ops() {
        let query = [0u8, 9, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0];
        let subject_bases = [0u8, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0];
        let subject = crate::encoding::pack_ncbi2na_bases(&subject_bases);
        let mut script = JumperPrelimEditBlock::default();
        let mut num_identical = 0;
        let mut ungapped_ext_len = -1;

        let (score, q_len, s_len) = jumper_extend_right_compressed_with_traceback(
            &query,
            &subject,
            subject_bases.len() as i32,
            2,
            -3,
            -5,
            -1,
            3,
            8,
            &mut script,
            &mut num_identical,
            false,
            &mut ungapped_ext_len,
        );

        assert_eq!(score, 20);
        assert_eq!((q_len, s_len), (14, 13));
        assert_eq!(script.edit_ops, vec![1, JUMPER_INSERTION, 10, 2]);
        assert_eq!(num_identical, 13);
        assert_eq!(ungapped_ext_len, 0);
    }

    #[test]
    fn jumper_extend_left_compressed_with_traceback_records_raw_ops() {
        let query = [0u8, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 9];
        let subject_bases = [0u8, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
        let subject = crate::encoding::pack_ncbi2na_bases(&subject_bases);
        let mut script = JumperPrelimEditBlock::default();
        let mut num_identical = 0;

        let (score, q_len, s_len) = jumper_extend_left_compressed_with_traceback(
            &query,
            &subject,
            12,
            11,
            2,
            -3,
            -5,
            -1,
            3,
            8,
            &mut script,
            &mut num_identical,
        );

        assert_eq!(score, 18);
        assert_eq!((q_len, s_len), (13, 12));
        assert_eq!(script.edit_ops, vec![JUMPER_INSERTION, 10, 2]);
        assert_eq!(num_identical, 12);
    }

    #[test]
    fn jumper_extend_right_compressed_optimal_keeps_best_prefix() {
        let query = [0u8, 9, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0];
        let subject_bases = [0u8, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0];
        let subject = crate::encoding::pack_ncbi2na_bases(&subject_bases);
        let mut script = JumperPrelimEditBlock::default();
        let mut best_num_identical = 0;
        let mut ungapped_ext_len = -1;

        let (score, q_len, s_len) = jumper_extend_right_compressed_with_traceback_optimal(
            &query,
            &subject,
            subject_bases.len() as i32,
            2,
            -3,
            -5,
            -1,
            3,
            8,
            100,
            &mut script,
            &mut best_num_identical,
            false,
            &mut ungapped_ext_len,
        );

        assert_eq!(score, 18);
        assert_eq!((q_len, s_len), (14, 13));
        assert_eq!(script.edit_ops, vec![1, JUMPER_INSERTION, 10, 2]);
        assert_eq!(best_num_identical, 13);
        assert_eq!(ungapped_ext_len, 0);
    }

    #[test]
    fn jumper_extend_left_compressed_optimal_keeps_best_prefix() {
        let query = [0u8, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 9];
        let subject_bases = [0u8, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
        let subject = crate::encoding::pack_ncbi2na_bases(&subject_bases);
        let mut script = JumperPrelimEditBlock::default();
        let mut best_num_identical = 0;

        let (score, q_len, s_len) = jumper_extend_left_compressed_with_traceback_optimal(
            &query,
            &subject,
            12,
            11,
            2,
            -3,
            -5,
            -1,
            3,
            8,
            100,
            &mut script,
            &mut best_num_identical,
        );

        assert_eq!(score, 18);
        assert_eq!((q_len, s_len), (13, 12));
        assert_eq!(script.edit_ops, vec![JUMPER_INSERTION, 10, 2]);
        assert_eq!(best_num_identical, 12);
    }

    #[test]
    fn jumper_gapped_alignment_compressed_rejects_malformed_sequence_shapes() {
        let query = [0u8, 1, 2, 3, 0, 1, 2, 3];
        let subject = crate::encoding::pack_ncbi2na_bases(&query);
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let align_params = JumperAlignParams {
            max_mismatches: 5,
            mismatch_window: 10,
            gap_x_dropoff: 30,
        };
        let mut jumper = JumperGapAlign {
            left_prelim_block: Some(JumperPrelimEditBlock::default()),
            right_prelim_block: Some(JumperPrelimEditBlock::default()),
            table: Vec::new(),
        };
        let mut num_identical = 17;
        let mut right_ungapped_ext_len = 19;

        assert_eq!(
            jumper_gapped_alignment_compressed_with_traceback(
                &query[..4],
                &subject,
                query.len() as i32,
                query.len() as i32,
                0,
                0,
                &mut jumper,
                &score_params,
                align_params,
                &mut num_identical,
                &mut right_ungapped_ext_len,
            ),
            -1
        );
        assert_eq!((num_identical, right_ungapped_ext_len), (17, 19));

        assert_eq!(
            jumper_gapped_alignment_compressed_with_traceback(
                &query,
                &subject[..1],
                query.len() as i32,
                query.len() as i32,
                0,
                0,
                &mut jumper,
                &score_params,
                align_params,
                &mut num_identical,
                &mut right_ungapped_ext_len,
            ),
            -1
        );
        assert_eq!((num_identical, right_ungapped_ext_len), (17, 19));
    }

    #[test]
    fn jumper_good_align_applies_identity_score_and_edit_distance() {
        let mut options = crate::options::HitSavingOptions {
            cutoff_score: 10,
            percent_identity: 80.0,
            max_edit_distance: 2,
            ..Default::default()
        };

        assert!(!jumper_good_align(12, 5, 5, 7, 7, 0, &options, 10));
        assert!(jumper_good_align(12, 0, 10, 5, 15, 8, &options, 10));
        assert!(!jumper_good_align(12, 0, 10, 5, 15, 7, &options, 10));
        assert!(!jumper_good_align(9, 0, 10, 5, 15, 8, &options, 10));
        assert!(!jumper_good_align(12, 0, 10, 5, 15, 7, &options, 10));

        options.cutoff_score = 0;
        options.cutoff_score_fun = [500, 100];
        assert!(jumper_good_align(15, 0, 10, 5, 15, 8, &options, 10));
        assert!(!jumper_good_align(14, 0, 10, 5, 15, 8, &options, 10));

        options.splice = true;
        options.cutoff_score = 1000;
        options.max_edit_distance = 0;
        assert!(jumper_good_align(12, 0, 10, 5, 15, 8, &options, 10));
    }

    #[test]
    fn jumper_find_edits_extracts_mismatch_insertions_and_deletions() {
        let query = [0u8, 1, 2, 3, 0];
        let subject_bases = [0u8, 2, 2, 3, 1];
        let subject = crate::encoding::pack_ncbi2na_bases(&subject_bases);
        let left = JumperPrelimEditBlock::default();
        let right = JumperPrelimEditBlock {
            edit_ops: vec![1, JUMPER_MISMATCH, JUMPER_INSERTION, JUMPER_DELETION, 1],
        };

        let edits = jumper_find_edits(&query, &subject, 0, 0, 4, 4, &left, &right).expect("edits");

        assert_eq!(
            edits.edits,
            vec![
                JumperEdit {
                    query_pos: 1,
                    query_base: 1,
                    subject_base: 2,
                },
                JumperEdit {
                    query_pos: 2,
                    query_base: 2,
                    subject_base: 15,
                },
                JumperEdit {
                    query_pos: 3,
                    query_base: 15,
                    subject_base: 2,
                },
            ]
        );
    }

    #[test]
    fn jumper_find_edits_rejects_invalid_shapes_ops_and_endpoints() {
        let query = [0u8, 1, 2, 3];
        let subject = crate::encoding::pack_ncbi2na_bases(&query);
        let empty = JumperPrelimEditBlock::default();
        let one_match = JumperPrelimEditBlock { edit_ops: vec![1] };

        assert!(jumper_find_edits(&query, &subject, -1, 0, 1, 1, &empty, &empty).is_none());
        assert!(jumper_find_edits(&query, &subject, 0, -1, 1, 1, &empty, &empty).is_none());
        assert!(jumper_find_edits(&query, &subject, 2, 0, 1, 1, &empty, &empty).is_none());
        assert!(jumper_find_edits(&query, &[], 0, 0, 1, 1, &empty, &empty).is_none());

        let invalid = JumperPrelimEditBlock {
            edit_ops: vec![-99],
        };
        assert!(jumper_find_edits(&query, &subject, 0, 0, 0, 0, &empty, &invalid).is_none());
        assert!(jumper_find_edits(&query, &subject, 0, 0, 2, 1, &empty, &one_match).is_none());
    }

    #[test]
    fn jumper_word_hit_hsp_records_mapper_info_edges_and_overhangs() {
        let query = [0u8, 1, 0, 0xff, 2, 3, 1, 1, 2, 3];
        let subject_bases = [0u8, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
        let subject = crate::encoding::pack_ncbi2na_bases(&subject_bases);

        let (hsp, map_info) =
            s_create_hsp_for_word_hit(2, 4, 4, 1, &query, 1, &subject, 16, 0, 10).expect("hsp");

        assert_eq!((hsp.query_offset, hsp.query_end), (2, 6));
        assert_eq!((hsp.subject_offset, hsp.subject_end), (4, 8));
        assert_eq!(hsp.num_ident, 3);
        assert_eq!(
            map_info.edits.as_ref().unwrap().edits,
            vec![JumperEdit {
                query_pos: 3,
                query_base: 0xff,
                subject_base: 1,
            }]
        );
        assert_eq!(map_info.left_edge, (2 << 2) | 3);
        assert_eq!(map_info.right_edge, 1);
        assert_eq!(
            map_info.subject_overhangs.as_ref().unwrap().left,
            Some(vec![2, 3])
        );
        assert_eq!(
            map_info.subject_overhangs.as_ref().unwrap().right,
            Some(vec![0, 1, 2, 3, 0])
        );
        assert_eq!(hsp.map_info.as_ref(), Some(&map_info));
    }

    #[test]
    fn jumper_word_hit_hsp_rejects_invalid_offsets_and_subject_shape() {
        let query = [0u8, 1, 2, 3];
        let subject = crate::encoding::pack_ncbi2na_bases(&query);

        assert!(s_create_hsp_for_word_hit(-1, 0, 4, 0, &query, 0, &subject, 4, 0, 4).is_none());
        assert!(s_create_hsp_for_word_hit(0, -1, 4, 0, &query, 0, &subject, 4, 0, 4).is_none());
        assert!(s_create_hsp_for_word_hit(0, 0, -1, 0, &query, 0, &subject, 4, 0, 4).is_none());
        assert!(s_create_hsp_for_word_hit(2, 0, 4, 0, &query, 0, &subject, 4, 0, 4).is_none());
        assert!(s_create_hsp_for_word_hit(0, 2, 4, 0, &query, 0, &subject, 4, 0, 4).is_none());
        assert!(s_create_hsp_for_word_hit(0, 0, 4, 0, &query, 0, &[], 4, 0, 4).is_none());
        assert!(s_create_hsp_for_word_hit(0, 0, 4, 0, &query, 0, &subject, -1, 0, 4).is_none());
    }

    #[test]
    fn save_subject_overhangs_keeps_terminal_right_bases_and_long_query_cap() {
        let subject_bases: Vec<u8> = (0..180).map(|idx| (idx % 4) as u8).collect();
        let subject = crate::encoding::pack_ncbi2na_bases(&subject_bases);
        let mut hsp = crate::hspstream::Hsp {
            score: 10,
            num_ident: 10,
            bit_score: 0.0,
            evalue: 0.0,
            query_offset: 0,
            query_end: 100,
            query_gapped_start: 0,
            subject_offset: 20,
            subject_end: 40,
            subject_gapped_start: 20,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        };
        let mut map_info = crate::hspstream::blast_hsp_mapping_info_new();

        assert_eq!(
            s_save_subject_overhangs(
                &hsp,
                &mut map_info,
                &subject,
                subject_bases.len() as i32,
                100,
            ),
            0
        );
        let overhangs = map_info.subject_overhangs.as_ref().expect("overhangs");
        assert_eq!(overhangs.left, Some(vec![2, 3]));
        assert_eq!(overhangs.right, Some(vec![0, 1]));

        hsp.query_offset = 90;
        hsp.query_end = 390;
        hsp.subject_offset = 80;
        hsp.subject_end = 100;
        let mut map_info = crate::hspstream::blast_hsp_mapping_info_new();
        assert_eq!(
            s_save_subject_overhangs(
                &hsp,
                &mut map_info,
                &subject,
                subject_bases.len() as i32,
                500,
            ),
            0
        );
        let overhangs = map_info.subject_overhangs.as_ref().expect("overhangs");
        assert_eq!(overhangs.left.as_ref().unwrap().len(), 60);
        assert_eq!(overhangs.right.as_ref().unwrap().len(), 60);
    }

    #[test]
    fn save_subject_overhangs_replaces_with_empty_container_when_no_bases_exist() {
        let subject_bases: Vec<u8> = (0..40).map(|idx| (idx % 4) as u8).collect();
        let subject = crate::encoding::pack_ncbi2na_bases(&subject_bases);
        let hsp = crate::hspstream::Hsp {
            score: 10,
            num_ident: 40,
            bit_score: 0.0,
            evalue: 0.0,
            query_offset: 0,
            query_end: 40,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: 40,
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
        let mut map_info = crate::hspstream::blast_hsp_mapping_info_new();
        map_info.subject_overhangs = Some(SequenceOverhangs {
            left: Some(vec![1, 2]),
            right: Some(vec![3]),
        });

        assert_eq!(
            s_save_subject_overhangs(
                &hsp,
                &mut map_info,
                &subject,
                subject_bases.len() as i32,
                40,
            ),
            0
        );

        let overhangs = map_info.subject_overhangs.as_ref().expect("overhangs");
        assert_eq!(overhangs.left, None);
        assert_eq!(overhangs.right, None);
    }

    #[test]
    fn save_subject_overhangs_rejects_invalid_subject_shape_and_offsets() {
        let subject_bases: Vec<u8> = (0..20).map(|idx| (idx % 4) as u8).collect();
        let subject = crate::encoding::pack_ncbi2na_bases(&subject_bases);
        let mut hsp = crate::hspstream::Hsp {
            score: 10,
            num_ident: 10,
            bit_score: 0.0,
            evalue: 0.0,
            query_offset: 2,
            query_end: 10,
            query_gapped_start: 2,
            subject_offset: 4,
            subject_end: 12,
            subject_gapped_start: 4,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        };
        let mut map_info = crate::hspstream::blast_hsp_mapping_info_new();

        assert_eq!(
            s_save_subject_overhangs(&hsp, &mut map_info, &[], 20, 20),
            -1
        );
        assert_eq!(
            s_save_subject_overhangs(&hsp, &mut map_info, &subject, -1, 20),
            -1
        );
        assert_eq!(
            s_save_subject_overhangs(&hsp, &mut map_info, &subject, 20, -1),
            -1
        );

        hsp.subject_offset = -1;
        assert_eq!(
            s_save_subject_overhangs(&hsp, &mut map_info, &subject, 20, 20),
            -1
        );
        hsp.subject_offset = 4;
        hsp.subject_end = 21;
        assert_eq!(
            s_save_subject_overhangs(&hsp, &mut map_info, &subject, 20, 20),
            -1
        );

        hsp.subject_end = 12;
        hsp.query_end = i32::MIN;
        assert_eq!(
            s_save_subject_overhangs(&hsp, &mut map_info, &subject, 20, 20),
            0
        );
        assert_eq!(
            map_info.subject_overhangs.as_ref().unwrap().right,
            Some(subject_bases[12..20].to_vec())
        );
    }

    #[test]
    fn jumper_find_splice_signals_marks_terminal_edges_as_exons() {
        let subject_bases: Vec<u8> = (0..20).map(|idx| (idx % 4) as u8).collect();
        let subject = crate::encoding::pack_ncbi2na_bases(&subject_bases);
        let hsp = crate::hspstream::Hsp {
            score: 10,
            num_ident: 20,
            bit_score: 0.0,
            evalue: 0.0,
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
        };
        let mut map_info = crate::hspstream::blast_hsp_mapping_info_new();

        assert_eq!(
            jumper_find_splice_signals(&hsp, &mut map_info, 20, &subject, 20),
            0
        );

        assert_eq!(map_info.left_edge, MAPPER_EXON);
        assert_eq!(map_info.right_edge, MAPPER_EXON);
    }

    #[test]
    fn jumper_find_splice_signals_records_internal_edge_bases() {
        let subject_bases = [0u8, 1, 2, 3, 3, 2, 1, 0, 0, 1, 2, 3];
        let subject = crate::encoding::pack_ncbi2na_bases(&subject_bases);
        let hsp = crate::hspstream::Hsp {
            score: 6,
            num_ident: 6,
            bit_score: 0.0,
            evalue: 0.0,
            query_offset: 2,
            query_end: 8,
            query_gapped_start: 2,
            subject_offset: 4,
            subject_end: 8,
            subject_gapped_start: 4,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        };
        let mut map_info = crate::hspstream::blast_hsp_mapping_info_new();

        assert_eq!(
            jumper_find_splice_signals(
                &hsp,
                &mut map_info,
                12,
                &subject,
                subject_bases.len() as i32,
            ),
            0
        );

        assert_eq!(map_info.left_edge, (2 << 2) | 3);
        assert_eq!(map_info.right_edge, (0 << 2) | 1);
    }

    #[test]
    fn jumper_find_splice_signals_rejects_invalid_shapes_and_ranges() {
        let subject_bases = [0u8, 1, 2, 3, 0, 1, 2, 3];
        let subject = crate::encoding::pack_ncbi2na_bases(&subject_bases);
        let mut hsp = crate::hspstream::Hsp {
            score: 4,
            num_ident: 4,
            bit_score: 0.0,
            evalue: 0.0,
            query_offset: 2,
            query_end: 6,
            query_gapped_start: 2,
            subject_offset: 2,
            subject_end: 6,
            subject_gapped_start: 2,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        };
        let mut map_info = crate::hspstream::blast_hsp_mapping_info_new();

        assert_eq!(
            jumper_find_splice_signals(&hsp, &mut map_info, 8, &[], 8),
            -1
        );
        assert_eq!(
            jumper_find_splice_signals(&hsp, &mut map_info, -1, &subject, 8),
            -1
        );
        assert_eq!(
            jumper_find_splice_signals(&hsp, &mut map_info, 8, &subject, -1),
            -1
        );
        hsp.query_offset = -1;
        assert_eq!(
            jumper_find_splice_signals(&hsp, &mut map_info, 8, &subject, 8),
            -1
        );
        hsp.query_offset = 2;
        hsp.subject_end = 9;
        assert_eq!(
            jumper_find_splice_signals(&hsp, &mut map_info, 8, &subject, 8),
            -1
        );
    }

    #[test]
    fn jumper_create_hsp_returns_hsp_and_separate_mapping_info() {
        let query = [0u8, 1, 2, 3, 0, 1, 2, 3];
        let subject_bases = [0u8, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
        let subject = crate::encoding::pack_ncbi2na_bases(&subject_bases);
        let mut jumper = JumperGapAlign {
            left_prelim_block: Some(JumperPrelimEditBlock::default()),
            right_prelim_block: Some(JumperPrelimEditBlock { edit_ops: vec![4] }),
            table: Vec::new(),
        };

        let (hsp, map_info) = s_create_hsp(
            &query,
            query.len() as i32,
            2,
            1,
            &subject,
            subject_bases.len() as i32,
            0,
            &mut jumper,
            2,
            6,
            3,
            7,
            8,
            -3,
            true,
        )
        .expect("hsp");

        assert_eq!(hsp.score, 8);
        assert_eq!((hsp.query_offset, hsp.query_end), (2, 6));
        assert_eq!((hsp.subject_offset, hsp.subject_end), (3, 7));
        assert_eq!(
            hsp.edit_script.as_ref().unwrap().ops,
            vec![(GapAlignOpType::Sub, 4)]
        );
        assert_eq!(
            map_info.edits.as_ref().unwrap().edits,
            Vec::<JumperEdit>::new()
        );
        assert_eq!(map_info.left_edge, (1 << 2) | 2);
        assert_eq!(map_info.right_edge, (3 << 2) | 0);
        assert!(map_info.subject_overhangs.is_some());
        assert_eq!(hsp.map_info.as_ref(), Some(&map_info));
    }

    #[test]
    fn jumper_create_hsp_without_splice_skips_edges_and_overhangs() {
        let query = [0u8, 1, 2, 3, 0, 1];
        let subject_bases = [0u8, 1, 2, 3, 0, 1, 2, 3];
        let subject = crate::encoding::pack_ncbi2na_bases(&subject_bases);
        let mut jumper = JumperGapAlign {
            left_prelim_block: Some(JumperPrelimEditBlock::default()),
            right_prelim_block: Some(JumperPrelimEditBlock { edit_ops: vec![6] }),
            table: Vec::new(),
        };

        let (hsp, map_info) = s_create_hsp(
            &query,
            query.len() as i32,
            0,
            0,
            &subject,
            subject_bases.len() as i32,
            0,
            &mut jumper,
            0,
            6,
            0,
            6,
            12,
            -3,
            false,
        )
        .expect("hsp");

        assert_eq!(hsp.score, 12);
        assert_eq!(hsp.map_info.as_ref(), Some(&map_info));
        assert_eq!(map_info.left_edge, 0);
        assert_eq!(map_info.right_edge, 0);
        assert!(map_info.subject_overhangs.is_none());
        assert_eq!(
            map_info.edits.as_ref().unwrap().edits,
            Vec::<JumperEdit>::new()
        );
    }

    #[test]
    fn jumper_create_hsp_rejects_missing_prelim_blocks() {
        let query = [0u8, 1, 2, 3];
        let subject = crate::encoding::pack_ncbi2na_bases(&query);

        let mut missing_left = JumperGapAlign {
            left_prelim_block: None,
            right_prelim_block: Some(JumperPrelimEditBlock { edit_ops: vec![4] }),
            table: Vec::new(),
        };
        assert!(s_create_hsp(
            &query,
            query.len() as i32,
            0,
            0,
            &subject,
            query.len() as i32,
            0,
            &mut missing_left,
            0,
            4,
            0,
            4,
            8,
            -3,
            true,
        )
        .is_none());

        let mut missing_right = JumperGapAlign {
            left_prelim_block: Some(JumperPrelimEditBlock::default()),
            right_prelim_block: None,
            table: Vec::new(),
        };
        assert!(s_create_hsp(
            &query,
            query.len() as i32,
            0,
            0,
            &subject,
            query.len() as i32,
            0,
            &mut missing_right,
            0,
            4,
            0,
            4,
            8,
            -3,
            true,
        )
        .is_none());
    }

    #[test]
    fn jumper_create_hsp_rejects_invalid_ranges_and_sequence_shapes() {
        let query = [0u8, 1, 2, 3];
        let subject = crate::encoding::pack_ncbi2na_bases(&query);
        let mut jumper = JumperGapAlign {
            left_prelim_block: Some(JumperPrelimEditBlock::default()),
            right_prelim_block: Some(JumperPrelimEditBlock { edit_ops: vec![4] }),
            table: Vec::new(),
        };

        assert!(s_create_hsp(
            &query,
            5,
            0,
            0,
            &subject,
            4,
            0,
            &mut jumper,
            0,
            4,
            0,
            4,
            8,
            -3,
            true,
        )
        .is_none());
        assert!(s_create_hsp(
            &query,
            4,
            0,
            0,
            &[],
            4,
            0,
            &mut jumper,
            0,
            4,
            0,
            4,
            8,
            -3,
            true,
        )
        .is_none());
        assert!(s_create_hsp(
            &query,
            4,
            0,
            0,
            &subject,
            4,
            0,
            &mut jumper,
            -1,
            4,
            0,
            4,
            8,
            -3,
            true,
        )
        .is_none());
        assert!(s_create_hsp(
            &query,
            4,
            0,
            0,
            &subject,
            4,
            0,
            &mut jumper,
            3,
            2,
            0,
            4,
            8,
            -3,
            true,
        )
        .is_none());
        assert!(s_create_hsp(
            &query,
            4,
            0,
            0,
            &subject,
            4,
            0,
            &mut jumper,
            0,
            4,
            0,
            5,
            8,
            -3,
            true,
        )
        .is_none());
    }

    #[test]
    fn blast_na_extend_jumper_sorts_extends_and_saves_hsp() {
        let query = [0u8, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
        let subject = crate::encoding::pack_ncbi2na_bases(&query);
        let mut pairs = vec![
            crate::lookup::OffsetPair {
                query_offset: 4,
                subject_offset: 4,
            },
            crate::lookup::OffsetPair {
                query_offset: 0,
                subject_offset: 0,
            },
        ];
        let lookup = crate::lookup::LookupTableWrap::Megablast(crate::lookup::MbLookupTable {
            word_length: 4,
            lut_word_length: 4,
            discontiguous: false,
            template_length: 0,
            template_type: crate::lookup::DiscTemplateType::Contiguous,
            two_templates: false,
            second_template_type: crate::lookup::DiscTemplateType::Contiguous,
            hashtable: Vec::new(),
            hashtable2: Vec::new(),
            next_pos: Vec::new(),
            next_pos2: Vec::new(),
            pv_array: Vec::new(),
            pv_array_bts: 0,
            longest_chain: 0,
            scan_step: 1,
        });
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let word_options = crate::options::InitialWordOptions::new_blastn();
        let word_params = crate::parameters::InitialWordParameters {
            options: word_options,
            x_dropoff_max: 0,
            cutoff_score_min: 0,
            cutoffs: Vec::new(),
            ungapped_extension: true,
            nucl_score_table: crate::parameters::InitialWordParameters::build_nucl_score_table(
                score_params.reward,
                score_params.penalty,
            ),
        };
        let hit_params = crate::parameters::HitSavingParameters {
            options: crate::options::HitSavingOptions {
                cutoff_score: 1,
                ..Default::default()
            },
            cutoff_score_min: 0,
            low_score: Vec::new(),
            cutoffs: Vec::new(),
            link_hsp_params: None,
            prelim_evalue: 10.0,
        };
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[query.len()]);
        let mut jumper = JumperGapAlign {
            left_prelim_block: Some(JumperPrelimEditBlock::default()),
            right_prelim_block: Some(JumperPrelimEditBlock::default()),
            table: Vec::new(),
        };
        let mut list = crate::hspstream::HspList::new(7);

        let extended = blast_na_extend_jumper(
            &mut pairs,
            &word_params,
            &score_params,
            &hit_params,
            &lookup,
            &query,
            &subject,
            query.len() as i32,
            0,
            &query_info,
            JumperAlignParams {
                max_mismatches: 5,
                mismatch_window: 10,
                gap_x_dropoff: 30,
            },
            &mut jumper,
            &mut list,
            query.len() as u32,
            None,
        );

        assert!(extended >= 1);
        assert!(!list.hsps.is_empty());
        assert_eq!(list.hsps[0].query_offset, 0);
        assert_eq!(list.hsps[0].subject_offset, 0);
        assert!(list.hsps[0].score > 0);
    }

    #[test]
    fn blast_na_extend_jumper_uses_mb_lookup_width_for_discontiguous_preextension() {
        let query = [0u8, 1, 2, 3, 0, 1, 2, 3];
        let subject_bases = [1u8, 1, 2, 3, 0, 1, 2, 3];
        let subject = crate::encoding::pack_ncbi2na_bases(&subject_bases);
        let mut pairs = vec![crate::lookup::OffsetPair {
            query_offset: 2,
            subject_offset: 2,
        }];
        let lookup = crate::lookup::LookupTableWrap::Megablast(crate::lookup::MbLookupTable {
            word_length: 6,
            lut_word_length: 4,
            discontiguous: true,
            template_length: 6,
            template_type: crate::lookup::DiscTemplateType::Template11_18Coding,
            two_templates: false,
            second_template_type: crate::lookup::DiscTemplateType::Contiguous,
            hashtable: Vec::new(),
            hashtable2: Vec::new(),
            next_pos: Vec::new(),
            next_pos2: Vec::new(),
            pv_array: Vec::new(),
            pv_array_bts: 0,
            longest_chain: 0,
            scan_step: 1,
        });
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let word_params = crate::parameters::InitialWordParameters {
            options: crate::options::InitialWordOptions::new_blastn(),
            x_dropoff_max: 0,
            cutoff_score_min: 0,
            cutoffs: Vec::new(),
            ungapped_extension: true,
            nucl_score_table: crate::parameters::InitialWordParameters::build_nucl_score_table(
                score_params.reward,
                score_params.penalty,
            ),
        };
        let hit_params = crate::parameters::HitSavingParameters {
            options: crate::options::HitSavingOptions {
                cutoff_score: 1,
                ..Default::default()
            },
            cutoff_score_min: 0,
            low_score: Vec::new(),
            cutoffs: Vec::new(),
            link_hsp_params: None,
            prelim_evalue: 10.0,
        };
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[query.len()]);
        let mut jumper = JumperGapAlign {
            left_prelim_block: Some(JumperPrelimEditBlock::default()),
            right_prelim_block: Some(JumperPrelimEditBlock::default()),
            table: Vec::new(),
        };
        let mut list = crate::hspstream::HspList::new(7);

        let extended = blast_na_extend_jumper(
            &mut pairs,
            &word_params,
            &score_params,
            &hit_params,
            &lookup,
            &query,
            &subject,
            query.len() as i32,
            0,
            &query_info,
            JumperAlignParams {
                max_mismatches: 5,
                mismatch_window: 10,
                gap_x_dropoff: 30,
            },
            &mut jumper,
            &mut list,
            6,
            None,
        );

        assert_eq!(extended, 0);
        assert!(list.hsps.is_empty());
    }

    #[test]
    fn blast_na_extend_jumper_skips_invalid_initial_word_offsets() {
        let query = [0u8, 1, 2, 3, 0, 1, 2, 3];
        let subject = crate::encoding::pack_ncbi2na_bases(&query);
        let mut pairs = vec![
            crate::lookup::OffsetPair {
                query_offset: -1,
                subject_offset: 0,
            },
            crate::lookup::OffsetPair {
                query_offset: 0,
                subject_offset: -1,
            },
            crate::lookup::OffsetPair {
                query_offset: 99,
                subject_offset: 0,
            },
            crate::lookup::OffsetPair {
                query_offset: 0,
                subject_offset: 99,
            },
            crate::lookup::OffsetPair {
                query_offset: i32::MIN,
                subject_offset: i32::MAX,
            },
            crate::lookup::OffsetPair {
                query_offset: i32::MAX,
                subject_offset: i32::MIN,
            },
        ];
        let lookup = crate::lookup::LookupTableWrap::Megablast(crate::lookup::MbLookupTable {
            word_length: 4,
            lut_word_length: 4,
            discontiguous: false,
            template_length: 0,
            template_type: crate::lookup::DiscTemplateType::Contiguous,
            two_templates: false,
            second_template_type: crate::lookup::DiscTemplateType::Contiguous,
            hashtable: Vec::new(),
            hashtable2: Vec::new(),
            next_pos: Vec::new(),
            next_pos2: Vec::new(),
            pv_array: Vec::new(),
            pv_array_bts: 0,
            longest_chain: 0,
            scan_step: 1,
        });
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let word_params = crate::parameters::InitialWordParameters {
            options: crate::options::InitialWordOptions::new_blastn(),
            x_dropoff_max: 0,
            cutoff_score_min: 0,
            cutoffs: Vec::new(),
            ungapped_extension: true,
            nucl_score_table: crate::parameters::InitialWordParameters::build_nucl_score_table(
                score_params.reward,
                score_params.penalty,
            ),
        };
        let hit_params = crate::parameters::HitSavingParameters {
            options: crate::options::HitSavingOptions {
                cutoff_score: 1,
                ..Default::default()
            },
            cutoff_score_min: 0,
            low_score: Vec::new(),
            cutoffs: Vec::new(),
            link_hsp_params: None,
            prelim_evalue: 10.0,
        };
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[query.len()]);
        let mut jumper = JumperGapAlign {
            left_prelim_block: Some(JumperPrelimEditBlock::default()),
            right_prelim_block: Some(JumperPrelimEditBlock::default()),
            table: Vec::new(),
        };
        let mut list = crate::hspstream::HspList::new(7);

        assert_eq!(
            blast_na_extend_jumper(
                &mut pairs,
                &word_params,
                &score_params,
                &hit_params,
                &lookup,
                &query,
                &subject,
                query.len() as i32,
                0,
                &query_info,
                JumperAlignParams {
                    max_mismatches: 5,
                    mismatch_window: 10,
                    gap_x_dropoff: 30,
                },
                &mut jumper,
                &mut list,
                query.len() as u32,
                None,
            ),
            0
        );
        assert!(list.hsps.is_empty());

        let mut malformed_subject_pairs = vec![crate::lookup::OffsetPair {
            query_offset: 4,
            subject_offset: 4,
        }];
        let split_lookup =
            crate::lookup::LookupTableWrap::Megablast(crate::lookup::MbLookupTable {
                word_length: 8,
                lut_word_length: 4,
                discontiguous: false,
                template_length: 0,
                template_type: crate::lookup::DiscTemplateType::Contiguous,
                two_templates: false,
                second_template_type: crate::lookup::DiscTemplateType::Contiguous,
                hashtable: Vec::new(),
                hashtable2: Vec::new(),
                next_pos: Vec::new(),
                next_pos2: Vec::new(),
                pv_array: Vec::new(),
                pv_array_bts: 0,
                longest_chain: 0,
                scan_step: 1,
            });
        let mut malformed_subject_list = crate::hspstream::HspList::new(7);

        assert_eq!(
            blast_na_extend_jumper(
                &mut malformed_subject_pairs,
                &word_params,
                &score_params,
                &hit_params,
                &split_lookup,
                &query,
                &[],
                query.len() as i32,
                0,
                &query_info,
                JumperAlignParams {
                    max_mismatches: 5,
                    mismatch_window: 10,
                    gap_x_dropoff: 30,
                },
                &mut jumper,
                &mut malformed_subject_list,
                query.len() as u32,
                None,
            ),
            0
        );
        assert!(malformed_subject_list.hsps.is_empty());

        let mut invalid_lookup_pairs = vec![crate::lookup::OffsetPair {
            query_offset: 0,
            subject_offset: 0,
        }];
        let invalid_lookup =
            crate::lookup::LookupTableWrap::Megablast(crate::lookup::MbLookupTable {
                word_length: 3,
                lut_word_length: 4,
                discontiguous: false,
                template_length: 0,
                template_type: crate::lookup::DiscTemplateType::Contiguous,
                two_templates: false,
                second_template_type: crate::lookup::DiscTemplateType::Contiguous,
                hashtable: Vec::new(),
                hashtable2: Vec::new(),
                next_pos: Vec::new(),
                next_pos2: Vec::new(),
                pv_array: Vec::new(),
                pv_array_bts: 0,
                longest_chain: 0,
                scan_step: 1,
            });
        let mut invalid_lookup_list = crate::hspstream::HspList::new(7);

        assert_eq!(
            blast_na_extend_jumper(
                &mut invalid_lookup_pairs,
                &word_params,
                &score_params,
                &hit_params,
                &invalid_lookup,
                &query,
                &subject,
                query.len() as i32,
                0,
                &query_info,
                JumperAlignParams {
                    max_mismatches: 5,
                    mismatch_window: 10,
                    gap_x_dropoff: 30,
                },
                &mut jumper,
                &mut invalid_lookup_list,
                query.len() as u32,
                None,
            ),
            -1
        );
        assert!(invalid_lookup_list.hsps.is_empty());

        let mut truncated_query_pairs = vec![crate::lookup::OffsetPair {
            query_offset: 0,
            subject_offset: 0,
        }];
        let truncated_query_info = crate::queryinfo::QueryInfo::new_blastp(&[16]);
        let mut truncated_query_list = crate::hspstream::HspList::new(7);

        assert_eq!(
            blast_na_extend_jumper(
                &mut truncated_query_pairs,
                &word_params,
                &score_params,
                &hit_params,
                &lookup,
                &query,
                &subject,
                query.len() as i32,
                0,
                &truncated_query_info,
                JumperAlignParams {
                    max_mismatches: 5,
                    mismatch_window: 10,
                    gap_x_dropoff: 30,
                },
                &mut jumper,
                &mut truncated_query_list,
                query.len() as u32,
                None,
            ),
            0
        );
        assert!(truncated_query_list.hsps.is_empty());
    }

    #[test]
    fn jumper_na_word_finder_scans_lookup_and_extends_hits() {
        let query = [0u8, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
        let subject = crate::encoding::pack_ncbi2na_bases(&query);
        let mut hashtable = vec![0; 1 << 8];
        hashtable[0b00011011] = 1;
        let lookup = crate::lookup::LookupTableWrap::Megablast(crate::lookup::MbLookupTable {
            word_length: 4,
            lut_word_length: 4,
            discontiguous: false,
            template_length: 0,
            template_type: crate::lookup::DiscTemplateType::Contiguous,
            two_templates: false,
            second_template_type: crate::lookup::DiscTemplateType::Contiguous,
            hashtable,
            hashtable2: Vec::new(),
            next_pos: vec![0, 0],
            next_pos2: Vec::new(),
            pv_array: Vec::new(),
            pv_array_bts: 0,
            longest_chain: 1,
            scan_step: 1,
        });
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let word_params = crate::parameters::InitialWordParameters {
            options: crate::options::InitialWordOptions::new_blastn(),
            x_dropoff_max: 0,
            cutoff_score_min: 0,
            cutoffs: Vec::new(),
            ungapped_extension: true,
            nucl_score_table: crate::parameters::InitialWordParameters::build_nucl_score_table(
                score_params.reward,
                score_params.penalty,
            ),
        };
        let hit_params = crate::parameters::HitSavingParameters {
            options: crate::options::HitSavingOptions {
                cutoff_score: 1,
                ..Default::default()
            },
            cutoff_score_min: 0,
            low_score: Vec::new(),
            cutoffs: Vec::new(),
            link_hsp_params: None,
            prelim_evalue: 10.0,
        };
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[query.len()]);
        let mut jumper = JumperGapAlign {
            left_prelim_block: Some(JumperPrelimEditBlock::default()),
            right_prelim_block: Some(JumperPrelimEditBlock::default()),
            table: Vec::new(),
        };
        let mut list = crate::hspstream::HspList::new(9);

        assert_eq!(
            jumper_na_word_finder(
                &subject,
                query.len() as i32,
                0,
                &query,
                &query_info,
                &lookup,
                &word_params,
                &score_params,
                &hit_params,
                JumperAlignParams {
                    max_mismatches: 5,
                    mismatch_window: 10,
                    gap_x_dropoff: 30,
                },
                &mut jumper,
                &mut list,
            ),
            0
        );
        assert!(!list.hsps.is_empty());
    }

    #[test]
    fn jumper_na_word_finder_scans_small_na_lookup_branch() {
        let query = [0u8, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
        let subject = crate::encoding::pack_ncbi2na_bases(&query);
        let word_at_one = 0b0110_1100usize;
        let mut backbone = vec![0; 1 << 8];
        backbone[word_at_one] = 1;
        let lookup = crate::lookup::LookupTableWrap::SmallNa(crate::lookup::SmallNaLookupTable {
            word_length: 4,
            backbone,
            overflow: Vec::new(),
            pv_array: Vec::new(),
            longest_chain: 1,
            scan_step: 1,
        });
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let word_params = crate::parameters::InitialWordParameters {
            options: crate::options::InitialWordOptions::new_blastn(),
            x_dropoff_max: 0,
            cutoff_score_min: 0,
            cutoffs: Vec::new(),
            ungapped_extension: true,
            nucl_score_table: crate::parameters::InitialWordParameters::build_nucl_score_table(
                score_params.reward,
                score_params.penalty,
            ),
        };
        let hit_params = crate::parameters::HitSavingParameters {
            options: crate::options::HitSavingOptions {
                cutoff_score: 1,
                ..Default::default()
            },
            cutoff_score_min: 0,
            low_score: Vec::new(),
            cutoffs: Vec::new(),
            link_hsp_params: None,
            prelim_evalue: 10.0,
        };
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[query.len()]);
        let mut jumper = JumperGapAlign {
            left_prelim_block: Some(JumperPrelimEditBlock::default()),
            right_prelim_block: Some(JumperPrelimEditBlock::default()),
            table: Vec::new(),
        };
        let mut list = crate::hspstream::HspList::new(9);

        assert_eq!(
            jumper_na_word_finder(
                &subject,
                query.len() as i32,
                0,
                &query,
                &query_info,
                &lookup,
                &word_params,
                &score_params,
                &hit_params,
                JumperAlignParams {
                    max_mismatches: 5,
                    mismatch_window: 10,
                    gap_x_dropoff: 30,
                },
                &mut jumper,
                &mut list,
            ),
            0
        );
        assert!(!list.hsps.is_empty());
        assert_eq!(list.hsps[0].query_offset, 0);
        assert_eq!(list.hsps[0].subject_offset, 0);
    }

    #[test]
    fn jumper_collect_word_hits_preserves_small_na_overflow_chains() {
        let query = [0u8, 1, 2, 3, 0, 1, 2, 3];
        let subject = crate::encoding::pack_ncbi2na_bases(&query);
        let word = 0b0001_1011usize;
        let mut backbone = vec![-1; 1 << 8];
        backbone[word] = -2;
        let lookup = crate::lookup::LookupTableWrap::SmallNa(crate::lookup::SmallNaLookupTable {
            word_length: 4,
            backbone,
            overflow: vec![0, 0, 0, 4, -1],
            pv_array: Vec::new(),
            longest_chain: 2,
            scan_step: 1,
        });

        let hits = jumper_collect_word_hits_in_ranges(
            &lookup,
            &subject,
            query.len() as i32,
            4,
            4,
            None,
            true,
        );
        assert!(hits.contains(&crate::lookup::OffsetPair {
            query_offset: 0,
            subject_offset: 0,
        }));
        assert!(hits.contains(&crate::lookup::OffsetPair {
            query_offset: 4,
            subject_offset: 0,
        }));
    }

    #[test]
    fn jumper_collect_word_hits_honors_standard_na_scan_step() {
        let subject_bases = [0u8; 8];
        let subject = crate::encoding::pack_ncbi2na_bases(&subject_bases);
        let mut thick_backbone = vec![crate::lookup::NaLookupBackboneCell::default(); 1 << 8];
        thick_backbone[0].num_used = 1;
        thick_backbone[0].entries[0] = 7;
        let lookup = crate::lookup::LookupTableWrap::Na(crate::lookup::BlastNaLookupTable {
            mask: 255,
            word_length: 5,
            lut_word_length: 4,
            scan_step: 2,
            backbone_size: 1 << 8,
            longest_chain: 1,
            thick_backbone,
            overflow: Vec::new(),
            overflow_size: 0,
            pv: Vec::new(),
        });

        let hits = jumper_collect_word_hits_in_ranges(
            &lookup,
            &subject,
            subject_bases.len() as i32,
            5,
            4,
            None,
            true,
        );
        let subject_offsets: Vec<i32> = hits.iter().map(|hit| hit.subject_offset).collect();
        assert_eq!(subject_offsets, vec![0, 2, 4]);
        assert!(hits.iter().all(|hit| hit.query_offset == 7));
    }

    #[test]
    fn jumper_collect_word_hits_honors_megablast_scan_step() {
        let subject_bases = [0u8; 8];
        let subject = crate::encoding::pack_ncbi2na_bases(&subject_bases);
        let mut hashtable = vec![0; 1 << 8];
        hashtable[0] = 1;
        let lookup = crate::lookup::LookupTableWrap::Megablast(crate::lookup::MbLookupTable {
            word_length: 4,
            lut_word_length: 4,
            discontiguous: false,
            template_length: 0,
            template_type: crate::lookup::DiscTemplateType::Contiguous,
            two_templates: false,
            second_template_type: crate::lookup::DiscTemplateType::Contiguous,
            hashtable,
            hashtable2: Vec::new(),
            next_pos: vec![0, 0],
            next_pos2: Vec::new(),
            pv_array: Vec::new(),
            pv_array_bts: 0,
            longest_chain: 1,
            scan_step: 4,
        });

        let hits = jumper_collect_word_hits_in_ranges(
            &lookup,
            &subject,
            subject_bases.len() as i32,
            4,
            4,
            None,
            true,
        );
        let subject_offsets: Vec<i32> = hits.iter().map(|hit| hit.subject_offset).collect();
        assert_eq!(subject_offsets, vec![0, 4]);
        assert!(hits.iter().all(|hit| hit.query_offset == 0));
    }

    #[test]
    fn jumper_collect_word_hits_can_stop_after_first_mask_range_like_c() {
        let subject_bases = [0u8; 12];
        let subject = crate::encoding::pack_ncbi2na_bases(&subject_bases);
        let mut hashtable = vec![0; 1 << 8];
        hashtable[0] = 1;
        let lookup = crate::lookup::LookupTableWrap::Megablast(crate::lookup::MbLookupTable {
            word_length: 4,
            lut_word_length: 4,
            discontiguous: false,
            template_length: 0,
            template_type: crate::lookup::DiscTemplateType::Contiguous,
            two_templates: false,
            second_template_type: crate::lookup::DiscTemplateType::Contiguous,
            hashtable,
            hashtable2: Vec::new(),
            next_pos: vec![0, 0],
            next_pos2: Vec::new(),
            pv_array: Vec::new(),
            pv_array_bts: 0,
            longest_chain: 1,
            scan_step: 1,
        });
        let ranges = [
            crate::util::SSeqRange { left: 0, right: 4 },
            crate::util::SSeqRange { left: 8, right: 12 },
        ];

        let first_range_hits = jumper_collect_word_hits_in_ranges(
            &lookup,
            &subject,
            subject_bases.len() as i32,
            4,
            4,
            Some(&ranges),
            false,
        );
        let all_range_hits = jumper_collect_word_hits_in_ranges(
            &lookup,
            &subject,
            subject_bases.len() as i32,
            4,
            4,
            Some(&ranges),
            true,
        );
        let all_range_batches = jumper_collect_word_hit_batches_in_ranges(
            &lookup,
            &subject,
            subject_bases.len() as i32,
            4,
            4,
            Some(&ranges),
            true,
        );

        assert_eq!(
            first_range_hits
                .iter()
                .map(|hit| hit.subject_offset)
                .collect::<Vec<_>>(),
            vec![0]
        );
        assert_eq!(
            all_range_hits
                .iter()
                .map(|hit| hit.subject_offset)
                .collect::<Vec<_>>(),
            vec![0, 8]
        );
        assert_eq!(
            all_range_batches
                .iter()
                .map(|batch| batch.s_range)
                .collect::<Vec<_>>(),
            vec![4, 12]
        );
    }

    #[test]
    fn jumper_collect_word_hits_honors_small_na_scan_step() {
        let query = [0u8, 1, 2, 3, 0, 1, 2, 3];
        let subject = crate::encoding::pack_ncbi2na_bases(&query);
        let word = 0b0001_1011usize;
        let mut backbone = vec![-1; 1 << 8];
        backbone[word] = -2;
        let lookup = crate::lookup::LookupTableWrap::SmallNa(crate::lookup::SmallNaLookupTable {
            word_length: 4,
            backbone,
            overflow: vec![0, 0, 0, 4, -1],
            pv_array: Vec::new(),
            longest_chain: 2,
            scan_step: 4,
        });

        let hits = jumper_collect_word_hits_in_ranges(
            &lookup,
            &subject,
            query.len() as i32,
            4,
            4,
            None,
            true,
        );
        let subject_offsets: Vec<i32> = hits.iter().map(|hit| hit.subject_offset).collect();
        assert_eq!(subject_offsets, vec![0, 0, 4, 4]);
        assert!(hits.contains(&crate::lookup::OffsetPair {
            query_offset: 0,
            subject_offset: 4,
        }));
        assert!(hits.contains(&crate::lookup::OffsetPair {
            query_offset: 4,
            subject_offset: 4,
        }));
    }

    #[test]
    fn jumper_na_word_finder_handles_empty_subject_and_unsupported_lookup() {
        let query = [0u8, 1, 2, 3, 0, 1, 2, 3];
        let subject = crate::encoding::pack_ncbi2na_bases(&query);
        let lookup = crate::lookup::LookupTableWrap::Megablast(crate::lookup::MbLookupTable {
            word_length: 4,
            lut_word_length: 4,
            discontiguous: false,
            template_length: 0,
            template_type: crate::lookup::DiscTemplateType::Contiguous,
            two_templates: false,
            second_template_type: crate::lookup::DiscTemplateType::Contiguous,
            hashtable: vec![0; 1 << 8],
            hashtable2: Vec::new(),
            next_pos: Vec::new(),
            next_pos2: Vec::new(),
            pv_array: Vec::new(),
            pv_array_bts: 0,
            longest_chain: 0,
            scan_step: 1,
        });
        let aa_lookup = crate::lookup::LookupTableWrap::Aa(crate::lookup::AaLookupTable {
            word_length: 3,
            threshold: 11.0,
            backbone: Vec::new(),
            pv_array: Vec::new(),
        });
        let invalid_metadata_lookup =
            crate::lookup::LookupTableWrap::Megablast(crate::lookup::MbLookupTable {
                word_length: 3,
                lut_word_length: 4,
                discontiguous: false,
                template_length: 0,
                template_type: crate::lookup::DiscTemplateType::Contiguous,
                two_templates: false,
                second_template_type: crate::lookup::DiscTemplateType::Contiguous,
                hashtable: vec![0; 1 << 8],
                hashtable2: Vec::new(),
                next_pos: Vec::new(),
                next_pos2: Vec::new(),
                pv_array: Vec::new(),
                pv_array_bts: 0,
                longest_chain: 0,
                scan_step: 1,
            });
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let word_params = crate::parameters::InitialWordParameters {
            options: crate::options::InitialWordOptions::new_blastn(),
            x_dropoff_max: 0,
            cutoff_score_min: 0,
            cutoffs: Vec::new(),
            ungapped_extension: true,
            nucl_score_table: crate::parameters::InitialWordParameters::build_nucl_score_table(
                score_params.reward,
                score_params.penalty,
            ),
        };
        let hit_params = crate::parameters::HitSavingParameters {
            options: crate::options::HitSavingOptions {
                cutoff_score: 1,
                ..Default::default()
            },
            cutoff_score_min: 0,
            low_score: Vec::new(),
            cutoffs: Vec::new(),
            link_hsp_params: None,
            prelim_evalue: 10.0,
        };
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[query.len()]);
        let mut jumper = JumperGapAlign {
            left_prelim_block: Some(JumperPrelimEditBlock::default()),
            right_prelim_block: Some(JumperPrelimEditBlock::default()),
            table: Vec::new(),
        };
        let align_params = JumperAlignParams {
            max_mismatches: 5,
            mismatch_window: 10,
            gap_x_dropoff: 30,
        };
        let mut list = crate::hspstream::HspList::new(9);

        assert_eq!(
            jumper_na_word_finder(
                &subject,
                0,
                0,
                &query,
                &query_info,
                &lookup,
                &word_params,
                &score_params,
                &hit_params,
                align_params,
                &mut jumper,
                &mut list,
            ),
            0
        );
        assert!(list.hsps.is_empty());

        assert_eq!(
            jumper_na_word_finder(
                &subject,
                query.len() as i32,
                0,
                &query,
                &query_info,
                &aa_lookup,
                &word_params,
                &score_params,
                &hit_params,
                align_params,
                &mut jumper,
                &mut list,
            ),
            -1
        );
        assert!(list.hsps.is_empty());

        assert_eq!(
            jumper_na_word_finder(
                &subject,
                query.len() as i32,
                0,
                &query,
                &query_info,
                &invalid_metadata_lookup,
                &word_params,
                &score_params,
                &hit_params,
                align_params,
                &mut jumper,
                &mut list,
            ),
            -1
        );
        assert!(list.hsps.is_empty());
    }

    #[test]
    fn packed_subject_word_rejects_invalid_shapes_and_bounds() {
        let subject_bases = [0u8, 1, 2, 3, 0];
        let subject = crate::encoding::pack_ncbi2na_bases(&subject_bases);

        assert_eq!(
            packed_subject_word(&subject, subject_bases.len() as i32, 0, 4),
            Some(27)
        );
        assert_eq!(
            packed_subject_word(&subject, subject_bases.len() as i32, -1, 4),
            None
        );
        assert_eq!(
            packed_subject_word(&subject, subject_bases.len() as i32, 0, 0),
            None
        );
        assert_eq!(
            packed_subject_word(&subject, subject_bases.len() as i32, 2, 4),
            None
        );
        assert_eq!(packed_subject_word(&subject, i32::MAX, i32::MAX, 4), None);
        assert_eq!(
            packed_subject_word(&[], subject_bases.len() as i32, 0, 4),
            None
        );
    }

    #[test]
    fn jumper_collect_word_hits_in_ranges_noops_on_invalid_scan_inputs() {
        let query = [0u8, 1, 2, 3];
        let subject = crate::encoding::pack_ncbi2na_bases(&query);
        let lookup = crate::lookup::LookupTableWrap::Megablast(crate::lookup::MbLookupTable {
            word_length: 4,
            lut_word_length: 4,
            discontiguous: false,
            template_length: 0,
            template_type: crate::lookup::DiscTemplateType::Contiguous,
            two_templates: false,
            second_template_type: crate::lookup::DiscTemplateType::Contiguous,
            hashtable: vec![0; 1 << 8],
            hashtable2: Vec::new(),
            next_pos: Vec::new(),
            next_pos2: Vec::new(),
            pv_array: Vec::new(),
            pv_array_bts: 0,
            longest_chain: 0,
            scan_step: 1,
        });

        assert!(
            jumper_collect_word_hits_in_ranges(&lookup, &subject, 0, 4, 4, None, true).is_empty()
        );
        assert!(jumper_collect_word_hits_in_ranges(
            &lookup,
            &subject,
            query.len() as i32,
            0,
            4,
            None,
            true,
        )
        .is_empty());
        assert!(jumper_collect_word_hits_in_ranges(
            &lookup,
            &subject,
            query.len() as i32,
            4,
            0,
            None,
            true,
        )
        .is_empty());
        assert!(jumper_collect_word_hits_in_ranges(
            &lookup,
            &subject,
            query.len() as i32,
            4,
            4,
            Some(&[]),
            true,
        )
        .is_empty());
        assert!(jumper_collect_word_hits_in_ranges(&lookup, &[], 4, 4, 4, None, true).is_empty());
        assert!(jumper_collect_word_hits_in_ranges(
            &lookup,
            &subject,
            query.len() as i32,
            4,
            4,
            Some(&[
                crate::util::SSeqRange {
                    left: i32::MAX,
                    right: i32::MAX,
                },
                crate::util::SSeqRange {
                    left: i32::MIN,
                    right: i32::MIN,
                },
            ]),
            true,
        )
        .is_empty());

        let discontig_lookup =
            crate::lookup::LookupTableWrap::Megablast(crate::lookup::MbLookupTable {
                word_length: 11,
                lut_word_length: 11,
                discontiguous: true,
                template_length: 18,
                template_type: crate::lookup::DiscTemplateType::Template11_18Coding,
                two_templates: false,
                second_template_type: crate::lookup::DiscTemplateType::Contiguous,
                hashtable: Vec::new(),
                hashtable2: Vec::new(),
                next_pos: Vec::new(),
                next_pos2: Vec::new(),
                pv_array: Vec::new(),
                pv_array_bts: 0,
                longest_chain: 0,
                scan_step: 1,
            });
        assert!(
            jumper_collect_word_hits_in_ranges(&discontig_lookup, &[], 32, 18, 18, None, true,)
                .is_empty()
        );
    }

    #[test]
    fn jumper_na_word_finder_batches_mapper_word_hits_like_c() {
        let query = [0u8, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
        let subject = crate::encoding::pack_ncbi2na_bases(&query);
        let mut hashtable = vec![0; 1 << 8];
        hashtable[0b00011011] = 1;
        let lookup = crate::lookup::LookupTableWrap::Megablast(crate::lookup::MbLookupTable {
            word_length: 4,
            lut_word_length: 4,
            discontiguous: false,
            template_length: 0,
            template_type: crate::lookup::DiscTemplateType::Contiguous,
            two_templates: false,
            second_template_type: crate::lookup::DiscTemplateType::Contiguous,
            hashtable,
            hashtable2: Vec::new(),
            next_pos: vec![0, 0],
            next_pos2: Vec::new(),
            pv_array: Vec::new(),
            pv_array_bts: 0,
            longest_chain: 1,
            scan_step: 1,
        });
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let word_params = crate::parameters::InitialWordParameters {
            options: crate::options::InitialWordOptions::new_blastn(),
            x_dropoff_max: 0,
            cutoff_score_min: 0,
            cutoffs: Vec::new(),
            ungapped_extension: true,
            nucl_score_table: crate::parameters::InitialWordParameters::build_nucl_score_table(
                score_params.reward,
                score_params.penalty,
            ),
        };
        let hit_params = crate::parameters::HitSavingParameters {
            options: crate::options::HitSavingOptions {
                cutoff_score: 1,
                ..Default::default()
            },
            cutoff_score_min: 0,
            low_score: Vec::new(),
            cutoffs: Vec::new(),
            link_hsp_params: None,
            prelim_evalue: 10.0,
        };
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[query.len()]);
        let mut jumper = JumperGapAlign {
            left_prelim_block: Some(JumperPrelimEditBlock::default()),
            right_prelim_block: Some(JumperPrelimEditBlock::default()),
            table: Vec::new(),
        };
        let mut list = crate::hspstream::HspList::new(10);
        let mut word_hits = crate::lookup::MapperWordHits {
            pair_arrays: vec![vec![crate::lookup::OffsetPair {
                query_offset: 0,
                subject_offset: 0,
            }]],
            num: vec![0],
            num_arrays: 1,
            array_size: 1,
            divisor: query.len() as i32 + 1,
            last_diag: vec![0; query_info.contexts.len()],
            last_pos: vec![0; query_info.contexts.len()],
        };
        let mut ungapped_stats = crate::diagnostics::UngappedStats::default();
        let mut gapped_stats = crate::diagnostics::GappedStats::default();

        assert_eq!(
            jumper_na_word_finder_with_word_hits(
                &subject,
                query.len() as i32,
                0,
                &query,
                &query_info,
                &lookup,
                &word_params,
                &score_params,
                &hit_params,
                JumperAlignParams {
                    max_mismatches: 5,
                    mismatch_window: 10,
                    gap_x_dropoff: 30,
                },
                &mut jumper,
                &mut list,
                Some(&mut word_hits),
                Some(&mut ungapped_stats),
                Some(&mut gapped_stats),
            ),
            0
        );
        assert!(gapped_stats.extensions > 0);
        assert!(ungapped_stats.lookup_hits > 0);
        assert_eq!(ungapped_stats.good_init_extends, gapped_stats.extensions);
        assert!(word_hits.num.iter().all(|&num| num == 0));
        assert!(!list.hsps.is_empty());
    }

    #[test]
    fn jumper_na_word_finder_rejects_malformed_mapper_word_hits() {
        let query = [0u8, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
        let subject = crate::encoding::pack_ncbi2na_bases(&query);
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let word_params = crate::parameters::InitialWordParameters {
            options: crate::options::InitialWordOptions::new_blastn(),
            x_dropoff_max: 0,
            cutoff_score_min: 0,
            cutoffs: Vec::new(),
            ungapped_extension: true,
            nucl_score_table: crate::parameters::InitialWordParameters::build_nucl_score_table(
                score_params.reward,
                score_params.penalty,
            ),
        };
        let hit_params = crate::parameters::HitSavingParameters {
            options: crate::options::HitSavingOptions {
                cutoff_score: 1,
                ..Default::default()
            },
            cutoff_score_min: 0,
            low_score: Vec::new(),
            cutoffs: Vec::new(),
            link_hsp_params: None,
            prelim_evalue: 10.0,
        };
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[query.len()]);
        let align_params = JumperAlignParams {
            max_mismatches: 5,
            mismatch_window: 10,
            gap_x_dropoff: 30,
        };

        let mut divisor_zero_hashtable = vec![0; 1 << 8];
        divisor_zero_hashtable[0b00011011] = 1;
        let divisor_zero_lookup =
            crate::lookup::LookupTableWrap::Megablast(crate::lookup::MbLookupTable {
                word_length: 4,
                lut_word_length: 4,
                discontiguous: false,
                template_length: 0,
                template_type: crate::lookup::DiscTemplateType::Contiguous,
                two_templates: false,
                second_template_type: crate::lookup::DiscTemplateType::Contiguous,
                hashtable: divisor_zero_hashtable,
                hashtable2: Vec::new(),
                next_pos: vec![0, 0],
                next_pos2: Vec::new(),
                pv_array: Vec::new(),
                pv_array_bts: 0,
                longest_chain: 1,
                scan_step: 1,
            });
        let mut divisor_zero_hits = crate::lookup::MapperWordHits {
            pair_arrays: vec![vec![crate::lookup::OffsetPair {
                query_offset: 0,
                subject_offset: 0,
            }]],
            num: vec![0],
            num_arrays: 1,
            array_size: 1,
            divisor: 0,
            last_diag: vec![0; query_info.contexts.len()],
            last_pos: vec![0; query_info.contexts.len()],
        };
        let mut jumper = JumperGapAlign {
            left_prelim_block: Some(JumperPrelimEditBlock::default()),
            right_prelim_block: Some(JumperPrelimEditBlock::default()),
            table: Vec::new(),
        };
        let mut list = crate::hspstream::HspList::new(10);

        assert_eq!(
            jumper_na_word_finder_with_word_hits(
                &subject,
                query.len() as i32,
                0,
                &query,
                &query_info,
                &divisor_zero_lookup,
                &word_params,
                &score_params,
                &hit_params,
                align_params,
                &mut jumper,
                &mut list,
                Some(&mut divisor_zero_hits),
                None,
                None,
            ),
            -1
        );
        assert!(list.hsps.is_empty());

        let mut zero_array_size_hashtable = vec![0; 1 << 8];
        zero_array_size_hashtable[0b00011011] = 1;
        let zero_array_size_lookup =
            crate::lookup::LookupTableWrap::Megablast(crate::lookup::MbLookupTable {
                word_length: 4,
                lut_word_length: 4,
                discontiguous: false,
                template_length: 0,
                template_type: crate::lookup::DiscTemplateType::Contiguous,
                two_templates: false,
                second_template_type: crate::lookup::DiscTemplateType::Contiguous,
                hashtable: zero_array_size_hashtable,
                hashtable2: Vec::new(),
                next_pos: vec![0, 0],
                next_pos2: Vec::new(),
                pv_array: Vec::new(),
                pv_array_bts: 0,
                longest_chain: 1,
                scan_step: 1,
            });
        let mut zero_array_size_hits = crate::lookup::MapperWordHits {
            pair_arrays: vec![Vec::new()],
            num: vec![0],
            num_arrays: 1,
            array_size: 0,
            divisor: 1,
            last_diag: vec![0; query_info.contexts.len()],
            last_pos: vec![0; query_info.contexts.len()],
        };

        assert_eq!(
            jumper_na_word_finder_with_word_hits(
                &subject,
                query.len() as i32,
                0,
                &query,
                &query_info,
                &zero_array_size_lookup,
                &word_params,
                &score_params,
                &hit_params,
                align_params,
                &mut jumper,
                &mut list,
                Some(&mut zero_array_size_hits),
                None,
                None,
            ),
            -1
        );
        assert!(list.hsps.is_empty());

        let mut short_num_hashtable = vec![0; 1 << 8];
        short_num_hashtable[0b00011011] = 1;
        let short_num_lookup =
            crate::lookup::LookupTableWrap::Megablast(crate::lookup::MbLookupTable {
                word_length: 4,
                lut_word_length: 4,
                discontiguous: false,
                template_length: 0,
                template_type: crate::lookup::DiscTemplateType::Contiguous,
                two_templates: false,
                second_template_type: crate::lookup::DiscTemplateType::Contiguous,
                hashtable: short_num_hashtable,
                hashtable2: Vec::new(),
                next_pos: vec![0, 0],
                next_pos2: Vec::new(),
                pv_array: Vec::new(),
                pv_array_bts: 0,
                longest_chain: 1,
                scan_step: 1,
            });
        let mut short_num_hits = crate::lookup::MapperWordHits {
            pair_arrays: vec![vec![crate::lookup::OffsetPair {
                query_offset: 0,
                subject_offset: 0,
            }]],
            num: Vec::new(),
            num_arrays: 1,
            array_size: 1,
            divisor: 1,
            last_diag: vec![0; query_info.contexts.len()],
            last_pos: vec![0; query_info.contexts.len()],
        };

        assert_eq!(
            jumper_na_word_finder_with_word_hits(
                &subject,
                query.len() as i32,
                0,
                &query,
                &query_info,
                &short_num_lookup,
                &word_params,
                &score_params,
                &hit_params,
                align_params,
                &mut jumper,
                &mut list,
                Some(&mut short_num_hits),
                None,
                None,
            ),
            -1
        );
        assert!(list.hsps.is_empty());

        let mut short_pair_arrays_hashtable = vec![0; 1 << 8];
        short_pair_arrays_hashtable[0b00011011] = 1;
        let short_pair_arrays_lookup =
            crate::lookup::LookupTableWrap::Megablast(crate::lookup::MbLookupTable {
                word_length: 4,
                lut_word_length: 4,
                discontiguous: false,
                template_length: 0,
                template_type: crate::lookup::DiscTemplateType::Contiguous,
                two_templates: false,
                second_template_type: crate::lookup::DiscTemplateType::Contiguous,
                hashtable: short_pair_arrays_hashtable,
                hashtable2: Vec::new(),
                next_pos: vec![0, 0],
                next_pos2: Vec::new(),
                pv_array: Vec::new(),
                pv_array_bts: 0,
                longest_chain: 1,
                scan_step: 1,
            });
        let mut short_pair_arrays_hits = crate::lookup::MapperWordHits {
            pair_arrays: Vec::new(),
            num: vec![0],
            num_arrays: 1,
            array_size: 1,
            divisor: 1,
            last_diag: vec![0; query_info.contexts.len()],
            last_pos: vec![0; query_info.contexts.len()],
        };

        assert_eq!(
            jumper_na_word_finder_with_word_hits(
                &subject,
                query.len() as i32,
                0,
                &query,
                &query_info,
                &short_pair_arrays_lookup,
                &word_params,
                &score_params,
                &hit_params,
                align_params,
                &mut jumper,
                &mut list,
                Some(&mut short_pair_arrays_hits),
                None,
                None,
            ),
            -1
        );
        assert!(list.hsps.is_empty());

        let mut short_last_diag_hashtable = vec![0; 1 << 8];
        short_last_diag_hashtable[0b00011011] = 1;
        let short_last_diag_lookup =
            crate::lookup::LookupTableWrap::Megablast(crate::lookup::MbLookupTable {
                word_length: 4,
                lut_word_length: 4,
                discontiguous: false,
                template_length: 0,
                template_type: crate::lookup::DiscTemplateType::Contiguous,
                two_templates: false,
                second_template_type: crate::lookup::DiscTemplateType::Contiguous,
                hashtable: short_last_diag_hashtable,
                hashtable2: Vec::new(),
                next_pos: vec![0, 0],
                next_pos2: Vec::new(),
                pv_array: Vec::new(),
                pv_array_bts: 0,
                longest_chain: 1,
                scan_step: 1,
            });
        let mut short_last_diag_hits = crate::lookup::MapperWordHits {
            pair_arrays: vec![vec![crate::lookup::OffsetPair {
                query_offset: 0,
                subject_offset: 0,
            }]],
            num: vec![0],
            num_arrays: 1,
            array_size: 1,
            divisor: 1,
            last_diag: Vec::new(),
            last_pos: vec![0; query_info.contexts.len()],
        };

        assert_eq!(
            jumper_na_word_finder_with_word_hits(
                &subject,
                query.len() as i32,
                0,
                &query,
                &query_info,
                &short_last_diag_lookup,
                &word_params,
                &score_params,
                &hit_params,
                align_params,
                &mut jumper,
                &mut list,
                Some(&mut short_last_diag_hits),
                None,
                None,
            ),
            -1
        );
        assert!(list.hsps.is_empty());

        let mut short_last_pos_hashtable = vec![0; 1 << 8];
        short_last_pos_hashtable[0b00011011] = 1;
        let short_last_pos_lookup =
            crate::lookup::LookupTableWrap::Megablast(crate::lookup::MbLookupTable {
                word_length: 4,
                lut_word_length: 4,
                discontiguous: false,
                template_length: 0,
                template_type: crate::lookup::DiscTemplateType::Contiguous,
                two_templates: false,
                second_template_type: crate::lookup::DiscTemplateType::Contiguous,
                hashtable: short_last_pos_hashtable,
                hashtable2: Vec::new(),
                next_pos: vec![0, 0],
                next_pos2: Vec::new(),
                pv_array: Vec::new(),
                pv_array_bts: 0,
                longest_chain: 1,
                scan_step: 1,
            });
        let mut short_last_pos_hits = crate::lookup::MapperWordHits {
            pair_arrays: vec![vec![crate::lookup::OffsetPair {
                query_offset: 0,
                subject_offset: 0,
            }]],
            num: vec![0],
            num_arrays: 1,
            array_size: 1,
            divisor: 1,
            last_diag: vec![0; query_info.contexts.len()],
            last_pos: Vec::new(),
        };

        assert_eq!(
            jumper_na_word_finder_with_word_hits(
                &subject,
                query.len() as i32,
                0,
                &query,
                &query_info,
                &short_last_pos_lookup,
                &word_params,
                &score_params,
                &hit_params,
                align_params,
                &mut jumper,
                &mut list,
                Some(&mut short_last_pos_hits),
                None,
                None,
            ),
            -1
        );
        assert!(list.hsps.is_empty());

        let mut extreme_offset_hashtable = vec![0; 1 << 8];
        extreme_offset_hashtable[0b00011011] = i32::MAX;
        let extreme_offset_lookup =
            crate::lookup::LookupTableWrap::Megablast(crate::lookup::MbLookupTable {
                word_length: 4,
                lut_word_length: 4,
                discontiguous: false,
                template_length: 0,
                template_type: crate::lookup::DiscTemplateType::Contiguous,
                two_templates: false,
                second_template_type: crate::lookup::DiscTemplateType::Contiguous,
                hashtable: extreme_offset_hashtable,
                hashtable2: Vec::new(),
                next_pos: Vec::new(),
                next_pos2: Vec::new(),
                pv_array: Vec::new(),
                pv_array_bts: 0,
                longest_chain: 1,
                scan_step: 1,
            });
        let mut extreme_offset_hits = crate::lookup::MapperWordHits {
            pair_arrays: vec![Vec::new()],
            num: vec![0],
            num_arrays: 1,
            array_size: 1,
            divisor: 1,
            last_diag: vec![0; query_info.contexts.len()],
            last_pos: vec![0; query_info.contexts.len()],
        };

        assert_eq!(
            jumper_na_word_finder_with_word_hits(
                &subject,
                query.len() as i32,
                0,
                &query,
                &query_info,
                &extreme_offset_lookup,
                &word_params,
                &score_params,
                &hit_params,
                align_params,
                &mut jumper,
                &mut list,
                Some(&mut extreme_offset_hits),
                None,
                None,
            ),
            0
        );
        assert!(extreme_offset_hits.num.iter().all(|&num| num == 0));
        assert!(list.hsps.is_empty());

        let mut out_of_range_hashtable = vec![0; 1 << 8];
        out_of_range_hashtable[0b00011011] = 5;
        let out_of_range_lookup =
            crate::lookup::LookupTableWrap::Megablast(crate::lookup::MbLookupTable {
                word_length: 4,
                lut_word_length: 4,
                discontiguous: false,
                template_length: 0,
                template_type: crate::lookup::DiscTemplateType::Contiguous,
                two_templates: false,
                second_template_type: crate::lookup::DiscTemplateType::Contiguous,
                hashtable: out_of_range_hashtable,
                hashtable2: Vec::new(),
                next_pos: vec![0; 6],
                next_pos2: Vec::new(),
                pv_array: Vec::new(),
                pv_array_bts: 0,
                longest_chain: 1,
                scan_step: 1,
            });
        let mut out_of_range_hits = crate::lookup::MapperWordHits {
            pair_arrays: vec![vec![crate::lookup::OffsetPair {
                query_offset: 0,
                subject_offset: 0,
            }]],
            num: vec![0],
            num_arrays: 1,
            array_size: 1,
            divisor: 1,
            last_diag: vec![0; query_info.contexts.len()],
            last_pos: vec![0; query_info.contexts.len()],
        };

        assert_eq!(
            jumper_na_word_finder_with_word_hits(
                &subject,
                query.len() as i32,
                0,
                &query,
                &query_info,
                &out_of_range_lookup,
                &word_params,
                &score_params,
                &hit_params,
                align_params,
                &mut jumper,
                &mut list,
                Some(&mut out_of_range_hits),
                None,
                None,
            ),
            -1
        );
        assert!(list.hsps.is_empty());
    }

    #[test]
    fn jumper_na_word_finder_honors_masked_subject_scan_ranges() {
        let query = [0u8, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
        let mut subject_bases = vec![0u8, 1, 2, 3];
        subject_bases.extend_from_slice(&query);
        let subject = crate::encoding::pack_ncbi2na_bases(&subject_bases);
        let mut hashtable = vec![0; 1 << 8];
        hashtable[0b00011011] = 1;
        let lookup = crate::lookup::LookupTableWrap::Megablast(crate::lookup::MbLookupTable {
            word_length: 4,
            lut_word_length: 4,
            discontiguous: false,
            template_length: 0,
            template_type: crate::lookup::DiscTemplateType::Contiguous,
            two_templates: false,
            second_template_type: crate::lookup::DiscTemplateType::Contiguous,
            hashtable,
            hashtable2: Vec::new(),
            next_pos: vec![0, 0],
            next_pos2: Vec::new(),
            pv_array: Vec::new(),
            pv_array_bts: 0,
            longest_chain: 1,
            scan_step: 1,
        });
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let word_params = crate::parameters::InitialWordParameters {
            options: crate::options::InitialWordOptions::new_blastn(),
            x_dropoff_max: 0,
            cutoff_score_min: 0,
            cutoffs: Vec::new(),
            ungapped_extension: true,
            nucl_score_table: crate::parameters::InitialWordParameters::build_nucl_score_table(
                score_params.reward,
                score_params.penalty,
            ),
        };
        let hit_params = crate::parameters::HitSavingParameters {
            options: crate::options::HitSavingOptions {
                cutoff_score: 1,
                ..Default::default()
            },
            cutoff_score_min: 0,
            low_score: Vec::new(),
            cutoffs: Vec::new(),
            link_hsp_params: None,
            prelim_evalue: 10.0,
        };
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[query.len()]);
        let mut jumper = JumperGapAlign {
            left_prelim_block: Some(JumperPrelimEditBlock::default()),
            right_prelim_block: Some(JumperPrelimEditBlock::default()),
            table: Vec::new(),
        };
        let mut list = crate::hspstream::HspList::new(11);
        let ranges = [crate::util::SSeqRange {
            left: 4,
            right: subject_bases.len() as i32,
        }];

        assert_eq!(
            jumper_na_word_finder_with_subject_ranges(
                &subject,
                subject_bases.len() as i32,
                0,
                &ranges,
                &query,
                &query_info,
                &lookup,
                &word_params,
                &score_params,
                &hit_params,
                JumperAlignParams {
                    max_mismatches: 5,
                    mismatch_window: 10,
                    gap_x_dropoff: 30,
                },
                &mut jumper,
                &mut list,
                None,
                None,
                None,
            ),
            0
        );
        assert!(!list.hsps.is_empty());
        assert!(list.hsps.iter().all(|hsp| hsp.subject_offset >= 4));
    }

    #[test]
    fn determine_scanning_offsets_skips_short_mask_ranges_like_c() {
        let ranges = [
            crate::util::SSeqRange { left: 0, right: 3 },
            crate::util::SSeqRange { left: 5, right: 6 },
            crate::util::SSeqRange { left: 8, right: 16 },
        ];
        let mut scan_range = [0, ranges[0].left, ranges[0].right - 4];

        assert!(s_determine_scanning_offsets(&ranges, 4, 4, &mut scan_range));
        assert_eq!(scan_range, [2, 8, 12]);

        scan_range[1] = scan_range[2] + 1;
        assert!(!s_determine_scanning_offsets(
            &ranges,
            4,
            4,
            &mut scan_range
        ));
    }

    #[test]
    fn jumper_na_word_finder_scans_discontiguous_mb_templates() {
        let query: Vec<u8> = (0..18).map(|i| (i % 4) as u8).collect();
        let subject = crate::encoding::pack_ncbi2na_bases(&query);
        let mut mb_table = crate::lookup::MbLookupTable {
            word_length: 18,
            lut_word_length: 11,
            discontiguous: false,
            template_length: 0,
            template_type: crate::lookup::DiscTemplateType::Contiguous,
            two_templates: false,
            second_template_type: crate::lookup::DiscTemplateType::Contiguous,
            hashtable: vec![0; 1 << 22],
            hashtable2: Vec::new(),
            next_pos: Vec::new(),
            next_pos2: Vec::new(),
            pv_array: vec![0; (1 << 22) / 32],
            pv_array_bts: crate::stat::PV_ARRAY_BTS as i32,
            longest_chain: 0,
            scan_step: 0,
        };
        let lookup_options = crate::options::LookupTableOptions {
            word_size: 11,
            mb_template_length: 18,
            mb_template_type: crate::lookup::DiscWordType::Coding as i32,
            ..crate::options::LookupTableOptions::default()
        };
        assert_eq!(
            crate::lookup::s_fill_disc_mb_table(
                &query,
                &[crate::util::SSeqRange { left: 0, right: 17 }],
                &mut mb_table,
                &lookup_options,
            ),
            0
        );
        let lookup = crate::lookup::LookupTableWrap::Megablast(mb_table);

        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let word_params = crate::parameters::InitialWordParameters {
            options: crate::options::InitialWordOptions::new_blastn(),
            x_dropoff_max: 0,
            cutoff_score_min: 0,
            cutoffs: Vec::new(),
            ungapped_extension: true,
            nucl_score_table: crate::parameters::InitialWordParameters::build_nucl_score_table(
                score_params.reward,
                score_params.penalty,
            ),
        };
        let hit_params = crate::parameters::HitSavingParameters {
            options: crate::options::HitSavingOptions {
                cutoff_score: 1,
                ..Default::default()
            },
            cutoff_score_min: 0,
            low_score: Vec::new(),
            cutoffs: Vec::new(),
            link_hsp_params: None,
            prelim_evalue: 10.0,
        };
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[query.len()]);
        let mut jumper = JumperGapAlign {
            left_prelim_block: Some(JumperPrelimEditBlock::default()),
            right_prelim_block: Some(JumperPrelimEditBlock::default()),
            table: Vec::new(),
        };
        let mut list = crate::hspstream::HspList::new(12);

        assert_eq!(
            jumper_na_word_finder(
                &subject,
                query.len() as i32,
                0,
                &query,
                &query_info,
                &lookup,
                &word_params,
                &score_params,
                &hit_params,
                JumperAlignParams {
                    max_mismatches: 5,
                    mismatch_window: 10,
                    gap_x_dropoff: 30,
                },
                &mut jumper,
                &mut list,
            ),
            0
        );
        assert!(!list.hsps.is_empty());
        assert_eq!(list.hsps[0].query_offset, 0);
        assert_eq!(list.hsps[0].subject_offset, 0);
    }

    #[test]
    fn jumper_na_word_finder_extends_two_template_secondary_diagonals() {
        let query = crate::encoding::encode_blastna_sequence(
            b"ACGTACGTGACTTACCGTACGTACGTGACTTACCGTACGTACGTGACTTACCGTACGT",
        );
        let subject_bases = crate::encoding::encode_blastna_sequence(
            b"GGGGACGTACGTGACTTACCGTACGTACGTGACTTACCGTACGTACGTGACTTACCGTACGTCCCC",
        );
        let subject = crate::encoding::pack_ncbi2na_bases(&subject_bases);
        let mut mb_table = crate::lookup::MbLookupTable {
            word_length: 16,
            lut_word_length: 11,
            discontiguous: false,
            template_length: 0,
            template_type: crate::lookup::DiscTemplateType::Contiguous,
            two_templates: false,
            second_template_type: crate::lookup::DiscTemplateType::Contiguous,
            hashtable: vec![0; 1 << 22],
            hashtable2: Vec::new(),
            next_pos: Vec::new(),
            next_pos2: Vec::new(),
            pv_array: vec![0; (1 << 22) / 32],
            pv_array_bts: crate::stat::PV_ARRAY_BTS as i32,
            longest_chain: 0,
            scan_step: 0,
        };
        let lookup_options = crate::options::LookupTableOptions {
            word_size: 11,
            mb_template_length: 16,
            mb_template_type: crate::lookup::DiscWordType::TwoTemplates as i32,
            ..crate::options::LookupTableOptions::default()
        };
        assert_eq!(
            crate::lookup::s_fill_disc_mb_table(
                &query,
                &[crate::util::SSeqRange {
                    left: 0,
                    right: query.len() as i32 - 1,
                }],
                &mut mb_table,
                &lookup_options,
            ),
            0
        );
        let lookup = crate::lookup::LookupTableWrap::Megablast(mb_table);

        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let word_params = crate::parameters::InitialWordParameters {
            options: crate::options::InitialWordOptions::new_blastn(),
            x_dropoff_max: 0,
            cutoff_score_min: 0,
            cutoffs: Vec::new(),
            ungapped_extension: true,
            nucl_score_table: crate::parameters::InitialWordParameters::build_nucl_score_table(
                score_params.reward,
                score_params.penalty,
            ),
        };
        let hit_params = crate::parameters::HitSavingParameters {
            options: crate::options::HitSavingOptions {
                cutoff_score: 1,
                ..Default::default()
            },
            cutoff_score_min: 0,
            low_score: Vec::new(),
            cutoffs: Vec::new(),
            link_hsp_params: None,
            prelim_evalue: 10.0,
        };
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[query.len()]);
        let mut jumper = JumperGapAlign {
            left_prelim_block: Some(JumperPrelimEditBlock::default()),
            right_prelim_block: Some(JumperPrelimEditBlock::default()),
            table: Vec::new(),
        };
        let mut list = crate::hspstream::HspList::new(12);

        assert_eq!(
            jumper_na_word_finder(
                &subject,
                subject_bases.len() as i32,
                0,
                &query,
                &query_info,
                &lookup,
                &word_params,
                &score_params,
                &hit_params,
                JumperAlignParams {
                    max_mismatches: 5,
                    mismatch_window: 10,
                    gap_x_dropoff: 30,
                },
                &mut jumper,
                &mut list,
            ),
            0
        );
        assert!(
            list.hsps.iter().any(|hsp| hsp.query_offset == 18
                && hsp.query_end == 58
                && hsp.subject_offset == 4
                && hsp.subject_end == 44),
            "left-overlap HSP missing after two-template Jumper extension: {:?}",
            list.hsps
        );
        assert!(
            list.hsps.iter().any(|hsp| hsp.query_offset == 0
                && hsp.query_end == 40
                && hsp.subject_offset == 22
                && hsp.subject_end == 62),
            "right-overlap HSP missing after two-template Jumper extension: {:?}",
            list.hsps
        );
    }

    #[test]
    fn jumper_alignment_outputs_rejects_short_scratch_table() {
        let mut jumper = JumperGapAlign {
            left_prelim_block: Some(JumperPrelimEditBlock::default()),
            right_prelim_block: Some(JumperPrelimEditBlock::default()),
            table: vec![10, 1, 2, 3],
        };
        assert!(jumper_alignment_outputs(&jumper).is_none());

        jumper.table.push(4);
        assert_eq!(jumper_alignment_outputs(&jumper), Some((10, 1, 2, 3, 4)));
    }

    #[test]
    fn do_anchored_scan_saves_best_matching_flank_hsp() {
        let query: Vec<u8> = (0..32).map(|i| (i % 4) as u8).collect();
        let subject = crate::encoding::pack_ncbi2na_bases(&query);
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let hit_params = crate::parameters::HitSavingParameters {
            options: crate::options::HitSavingOptions {
                cutoff_score: 1,
                ..Default::default()
            },
            cutoff_score_min: 0,
            low_score: Vec::new(),
            cutoffs: Vec::new(),
            link_hsp_params: None,
            prelim_evalue: 10.0,
        };
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[query.len()]);
        let mut jumper = JumperGapAlign {
            left_prelim_block: Some(JumperPrelimEditBlock::default()),
            right_prelim_block: Some(JumperPrelimEditBlock::default()),
            table: Vec::new(),
        };
        let mut list = crate::hspstream::HspList::new(11);

        assert_eq!(
            do_anchored_scan(
                &query,
                query.len() as i32,
                0,
                0,
                &subject,
                query.len() as i32,
                0,
                0,
                query.len() as i32 - 1,
                &query_info,
                &mut jumper,
                &score_params,
                &hit_params,
                JumperAlignParams {
                    max_mismatches: 5,
                    mismatch_window: 10,
                    gap_x_dropoff: 30,
                },
                &mut list,
            ),
            1
        );
        assert_eq!(list.hsps[0].query_offset, 0);
        assert_eq!(list.hsps[0].subject_offset, 0);
    }

    #[test]
    fn anchored_scan_query_word_rejects_invalid_bounds_and_ambiguity() {
        let query = [0u8, 1, 2, 3, 0, 1, 2, 3];

        assert_eq!(unpacked_query_word(&query, 0, 4), Some(27));
        assert_eq!(unpacked_query_word(&query, -1, 4), None);
        assert_eq!(unpacked_query_word(&query, 0, 0), None);
        assert_eq!(unpacked_query_word(&query, 5, 4), None);
        assert_eq!(unpacked_query_word(&query, i32::MAX, 4), None);

        let ambiguous = [0u8, 1, 0xff, 3];
        assert_eq!(unpacked_query_word(&ambiguous, 0, 4), None);
    }

    #[test]
    fn do_anchored_scan_handles_left_flank_and_too_short_windows() {
        let query: Vec<u8> = (0..40).map(|i| (i % 4) as u8).collect();
        let subject_bases: Vec<u8> = (0..48).map(|i| (i % 4) as u8).collect();
        let subject = crate::encoding::pack_ncbi2na_bases(&subject_bases);
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let hit_params = crate::parameters::HitSavingParameters {
            options: crate::options::HitSavingOptions {
                cutoff_score: 1,
                ..Default::default()
            },
            cutoff_score_min: 0,
            low_score: Vec::new(),
            cutoffs: Vec::new(),
            link_hsp_params: None,
            prelim_evalue: 10.0,
        };
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[query.len()]);
        let mut jumper = JumperGapAlign {
            left_prelim_block: Some(JumperPrelimEditBlock::default()),
            right_prelim_block: Some(JumperPrelimEditBlock::default()),
            table: Vec::new(),
        };
        let align_params = JumperAlignParams {
            max_mismatches: 5,
            mismatch_window: 10,
            gap_x_dropoff: 30,
        };

        let mut list = crate::hspstream::HspList::new(11);
        assert_eq!(
            do_anchored_scan(
                &query,
                query.len() as i32,
                24,
                0,
                &subject,
                query.len() as i32,
                0,
                23,
                0,
                &query_info,
                &mut jumper,
                &score_params,
                &hit_params,
                align_params,
                &mut list,
            ),
            1
        );
        assert!(list.hsps.iter().any(|hsp| hsp.query_offset <= 12));

        let mut too_short = crate::hspstream::HspList::new(11);
        assert_eq!(
            do_anchored_scan(
                &query,
                query.len() as i32,
                8,
                0,
                &subject,
                query.len() as i32,
                0,
                7,
                0,
                &query_info,
                &mut jumper,
                &score_params,
                &hit_params,
                align_params,
                &mut too_short,
            ),
            0
        );
        assert!(too_short.hsps.is_empty());
    }

    #[test]
    fn do_anchored_scan_skips_ambiguous_and_polya_query_words_like_c() {
        let subject_bases = vec![0u8; 32];
        let subject = crate::encoding::pack_ncbi2na_bases(&subject_bases);
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let hit_params = crate::parameters::HitSavingParameters {
            options: crate::options::HitSavingOptions {
                cutoff_score: 1,
                ..Default::default()
            },
            cutoff_score_min: 0,
            low_score: Vec::new(),
            cutoffs: Vec::new(),
            link_hsp_params: None,
            prelim_evalue: 10.0,
        };
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[subject_bases.len()]);
        let mut jumper = JumperGapAlign {
            left_prelim_block: Some(JumperPrelimEditBlock::default()),
            right_prelim_block: Some(JumperPrelimEditBlock::default()),
            table: Vec::new(),
        };

        for query in [vec![0u8; 32], vec![4u8; 32]] {
            let mut list = crate::hspstream::HspList::new(11);
            assert_eq!(
                do_anchored_scan(
                    &query,
                    query.len() as i32,
                    0,
                    0,
                    &subject,
                    subject_bases.len() as i32,
                    0,
                    0,
                    subject_bases.len() as i32 - 1,
                    &query_info,
                    &mut jumper,
                    &score_params,
                    &hit_params,
                    JumperAlignParams {
                        max_mismatches: 5,
                        mismatch_window: 10,
                        gap_x_dropoff: 30,
                    },
                    &mut list,
                ),
                0
            );
            assert!(list.hsps.is_empty());
        }
    }

    #[test]
    fn do_anchored_scan_returns_empty_when_subject_words_do_not_match() {
        let query: Vec<u8> = (0..36).map(|i| (i % 4) as u8).collect();
        let subject_bases = vec![3u8; 36];
        let subject = crate::encoding::pack_ncbi2na_bases(&subject_bases);
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let hit_params = crate::parameters::HitSavingParameters {
            options: crate::options::HitSavingOptions {
                cutoff_score: 1,
                ..Default::default()
            },
            cutoff_score_min: 0,
            low_score: Vec::new(),
            cutoffs: Vec::new(),
            link_hsp_params: None,
            prelim_evalue: 10.0,
        };
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[query.len()]);
        let mut jumper = JumperGapAlign {
            left_prelim_block: Some(JumperPrelimEditBlock::default()),
            right_prelim_block: Some(JumperPrelimEditBlock::default()),
            table: Vec::new(),
        };
        let mut list = crate::hspstream::HspList::new(11);

        assert_eq!(
            do_anchored_scan(
                &query,
                query.len() as i32,
                0,
                0,
                &subject,
                subject_bases.len() as i32,
                0,
                0,
                subject_bases.len() as i32 - 1,
                &query_info,
                &mut jumper,
                &score_params,
                &hit_params,
                JumperAlignParams {
                    max_mismatches: 5,
                    mismatch_window: 10,
                    gap_x_dropoff: 30,
                },
                &mut list,
            ),
            0
        );
        assert!(list.hsps.is_empty());

        let mut malformed_subject = crate::hspstream::HspList::new(11);
        assert_eq!(
            do_anchored_scan(
                &query,
                query.len() as i32,
                0,
                0,
                &[],
                subject_bases.len() as i32,
                0,
                0,
                subject_bases.len() as i32 - 1,
                &query_info,
                &mut jumper,
                &score_params,
                &hit_params,
                JumperAlignParams {
                    max_mismatches: 5,
                    mismatch_window: 10,
                    gap_x_dropoff: 30,
                },
                &mut malformed_subject,
            ),
            0
        );
        assert!(malformed_subject.hsps.is_empty());

        let mut malformed_query = crate::hspstream::HspList::new(11);
        assert_eq!(
            do_anchored_scan(
                &query[..8],
                query.len() as i32,
                0,
                0,
                &subject,
                subject_bases.len() as i32,
                0,
                0,
                subject_bases.len() as i32 - 1,
                &query_info,
                &mut jumper,
                &score_params,
                &hit_params,
                JumperAlignParams {
                    max_mismatches: 5,
                    mismatch_window: 10,
                    gap_x_dropoff: 30,
                },
                &mut malformed_query,
            ),
            0
        );
        assert!(malformed_query.hsps.is_empty());

        let mut bad_query_from = crate::hspstream::HspList::new(11);
        assert_eq!(
            do_anchored_scan(
                &query,
                query.len() as i32,
                -1,
                0,
                &subject,
                subject_bases.len() as i32,
                0,
                0,
                subject_bases.len() as i32 - 1,
                &query_info,
                &mut jumper,
                &score_params,
                &hit_params,
                JumperAlignParams {
                    max_mismatches: 5,
                    mismatch_window: 10,
                    gap_x_dropoff: 30,
                },
                &mut bad_query_from,
            ),
            0
        );
        assert!(bad_query_from.hsps.is_empty());
    }

    #[test]
    fn do_anchored_search_scans_partial_mapper_chain_and_copies_original_hsps() {
        let query: Vec<u8> = (0..40).map(|i| (i % 4) as u8).collect();
        let subject = crate::encoding::pack_ncbi2na_bases(&query);
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let hit_params = crate::parameters::HitSavingParameters {
            options: crate::options::HitSavingOptions {
                cutoff_score: 1,
                ..Default::default()
            },
            cutoff_score_min: 0,
            low_score: Vec::new(),
            cutoffs: Vec::new(),
            link_hsp_params: None,
            prelim_evalue: 10.0,
        };
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[query.len()]);
        let chain_hsp = crate::hspstream::Hsp {
            score: 35,
            num_ident: 10,
            bit_score: 0.0,
            evalue: 0.0,
            query_offset: 13,
            query_end: 23,
            query_gapped_start: 13,
            subject_offset: 13,
            subject_end: 23,
            subject_gapped_start: 13,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        };
        let chain = Box::new(crate::spliced_hits::HSPChain {
            context: 0,
            oid: 21,
            score: 35,
            hsps: Some(Box::new(crate::spliced_hits::HSPContainer {
                hsp: chain_hsp,
                next: None,
            })),
            count: 1,
            pair: None,
            pair_score: None,
            pair_conf: crate::spliced_hits::PAIR_NONE,
            adapter: -1,
            poly_a: 0,
            next: None,
        });
        let mapper_data = crate::spliced_hits::BlastHSPMapperData {
            query_info: Some(query_info.clone()),
            saved_chains: vec![Some(chain)],
            ..Default::default()
        };
        let mut jumper = JumperGapAlign {
            left_prelim_block: Some(JumperPrelimEditBlock::default()),
            right_prelim_block: Some(JumperPrelimEditBlock::default()),
            table: Vec::new(),
        };

        let list = do_anchored_search(
            &query,
            &subject,
            query.len() as i32,
            21,
            0,
            12,
            12,
            Some(&mapper_data),
            &query_info,
            &mut jumper,
            &score_params,
            &hit_params,
            JumperAlignParams {
                max_mismatches: 5,
                mismatch_window: 10,
                gap_x_dropoff: 30,
            },
            None,
        )
        .expect("anchored search should produce a list");

        assert_eq!(list.oid, 21);
        assert!(list.hsps.iter().any(|hsp| hsp.query_offset == 0));
        assert!(list.hsps.iter().any(|hsp| hsp.query_offset == 13));
    }

    #[test]
    fn do_anchored_search_copies_multi_hsp_chain_after_flank_recovery() {
        let query: Vec<u8> = (0..44).map(|i| (i % 4) as u8).collect();
        let subject = crate::encoding::pack_ncbi2na_bases(&query);
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let hit_params = crate::parameters::HitSavingParameters {
            options: crate::options::HitSavingOptions {
                cutoff_score: 1,
                longest_intron: 16,
                ..Default::default()
            },
            cutoff_score_min: 0,
            low_score: Vec::new(),
            cutoffs: Vec::new(),
            link_hsp_params: None,
            prelim_evalue: 10.0,
        };
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[query.len()]);
        let second = Box::new(crate::spliced_hits::HSPContainer {
            hsp: mapper_saved_hsp(24, 28, 24, 28),
            next: None,
        });
        let first = Box::new(crate::spliced_hits::HSPContainer {
            hsp: mapper_saved_hsp(13, 23, 13, 23),
            next: Some(second),
        });
        let chain = Box::new(crate::spliced_hits::HSPChain {
            context: 0,
            oid: 28,
            score: 35,
            hsps: Some(first),
            count: 2,
            pair: None,
            pair_score: None,
            pair_conf: crate::spliced_hits::PAIR_NONE,
            adapter: -1,
            poly_a: 0,
            next: None,
        });
        let mapper_data = crate::spliced_hits::BlastHSPMapperData {
            query_info: Some(query_info.clone()),
            saved_chains: vec![Some(chain)],
            ..Default::default()
        };
        let mut jumper = JumperGapAlign {
            left_prelim_block: Some(JumperPrelimEditBlock::default()),
            right_prelim_block: Some(JumperPrelimEditBlock::default()),
            table: Vec::new(),
        };

        let list = do_anchored_search(
            &query,
            &subject,
            query.len() as i32,
            28,
            0,
            12,
            hit_params.options.longest_intron,
            Some(&mapper_data),
            &query_info,
            &mut jumper,
            &score_params,
            &hit_params,
            JumperAlignParams {
                max_mismatches: 5,
                mismatch_window: 10,
                gap_x_dropoff: 30,
            },
            None,
        )
        .expect("anchored search should recover a flank and copy the chain");

        assert_eq!(list.oid, 28);
        assert!(list
            .hsps
            .iter()
            .any(|hsp| hsp.query_offset == 0 && hsp.query_end >= 40));
        assert!(list
            .hsps
            .iter()
            .any(|hsp| hsp.query_offset == 13 && hsp.query_end == 23));
        assert!(list
            .hsps
            .iter()
            .any(|hsp| hsp.query_offset == 24 && hsp.query_end == 28));
    }

    #[test]
    fn do_anchored_search_uses_caller_word_size_for_opposite_flank_scan() {
        let query: Vec<u8> = (0..52).map(|i| (i % 4) as u8).collect();
        let subject = crate::encoding::pack_ncbi2na_bases(&query);
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let hit_params = crate::parameters::HitSavingParameters {
            options: crate::options::HitSavingOptions {
                cutoff_score: 1,
                longest_intron: 24,
                ..Default::default()
            },
            cutoff_score_min: 0,
            low_score: Vec::new(),
            cutoffs: Vec::new(),
            link_hsp_params: None,
            prelim_evalue: 10.0,
        };
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[query.len()]);
        let chain = Box::new(crate::spliced_hits::HSPChain {
            context: 0,
            oid: 28,
            score: 35,
            hsps: Some(Box::new(crate::spliced_hits::HSPContainer {
                hsp: mapper_saved_hsp(24, 35, 24, 35),
                next: None,
            })),
            count: 1,
            pair: None,
            pair_score: None,
            pair_conf: crate::spliced_hits::PAIR_NONE,
            adapter: -1,
            poly_a: 0,
            next: None,
        });
        let mapper_data = crate::spliced_hits::BlastHSPMapperData {
            query_info: Some(query_info.clone()),
            saved_chains: vec![Some(chain)],
            ..Default::default()
        };
        let mut jumper = JumperGapAlign {
            left_prelim_block: Some(JumperPrelimEditBlock::default()),
            right_prelim_block: Some(JumperPrelimEditBlock::default()),
            table: Vec::new(),
        };

        let list = do_anchored_search(
            &query,
            &subject,
            query.len() as i32,
            28,
            0,
            20,
            hit_params.options.longest_intron,
            Some(&mapper_data),
            &query_info,
            &mut jumper,
            &score_params,
            &hit_params,
            JumperAlignParams {
                max_mismatches: 5,
                mismatch_window: 10,
                gap_x_dropoff: 30,
            },
            None,
        )
        .expect("leading flank recovery should keep the original chain");

        assert!(list
            .hsps
            .iter()
            .any(|hsp| hsp.query_offset == 24 && hsp.query_end == 35));
        assert!(list.hsps.iter().any(|hsp| hsp.query_offset == 0));
        assert!(!list.hsps.iter().any(|hsp| hsp.query_offset == 35));
    }

    #[test]
    fn do_anchored_search_recovers_trailing_flank_and_copies_saved_chain() {
        let query: Vec<u8> = (0..52).map(|i| (i % 4) as u8).collect();
        let subject = crate::encoding::pack_ncbi2na_bases(&query);
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let hit_params = crate::parameters::HitSavingParameters {
            options: crate::options::HitSavingOptions {
                cutoff_score: 1,
                longest_intron: 40,
                ..Default::default()
            },
            cutoff_score_min: 0,
            low_score: Vec::new(),
            cutoffs: Vec::new(),
            link_hsp_params: None,
            prelim_evalue: 10.0,
        };
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[query.len()]);
        let chain = Box::new(crate::spliced_hits::HSPChain {
            context: 0,
            oid: 34,
            score: 35,
            hsps: Some(Box::new(crate::spliced_hits::HSPContainer {
                hsp: mapper_saved_hsp(0, 20, 0, 20),
                next: None,
            })),
            count: 1,
            pair: None,
            pair_score: None,
            pair_conf: crate::spliced_hits::PAIR_NONE,
            adapter: -1,
            poly_a: 0,
            next: None,
        });
        let mapper_data = crate::spliced_hits::BlastHSPMapperData {
            query_info: Some(query_info.clone()),
            saved_chains: vec![Some(chain)],
            ..Default::default()
        };
        let mut jumper = JumperGapAlign {
            left_prelim_block: Some(JumperPrelimEditBlock::default()),
            right_prelim_block: Some(JumperPrelimEditBlock::default()),
            table: Vec::new(),
        };

        let list = do_anchored_search(
            &query,
            &subject,
            query.len() as i32,
            34,
            0,
            12,
            hit_params.options.longest_intron,
            Some(&mapper_data),
            &query_info,
            &mut jumper,
            &score_params,
            &hit_params,
            JumperAlignParams {
                max_mismatches: 5,
                mismatch_window: 10,
                gap_x_dropoff: 30,
            },
            None,
        )
        .expect("anchored search should recover the uncovered trailing flank");

        assert_eq!(list.oid, 34);
        assert!(list
            .hsps
            .iter()
            .any(|hsp| hsp.query_offset == 0 && hsp.query_end == 20));
        assert!(list
            .hsps
            .iter()
            .any(|hsp| hsp.query_offset <= 20 && hsp.query_end == query.len() as i32));
    }

    #[test]
    fn do_anchored_search_does_not_copy_chain_without_recovered_flank() {
        let query = vec![1; 44];
        let subject_bases = vec![0; query.len()];
        let subject = crate::encoding::pack_ncbi2na_bases(&subject_bases);
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let hit_params = crate::parameters::HitSavingParameters {
            options: crate::options::HitSavingOptions {
                cutoff_score: 1,
                longest_intron: 16,
                ..Default::default()
            },
            cutoff_score_min: 0,
            low_score: Vec::new(),
            cutoffs: Vec::new(),
            link_hsp_params: None,
            prelim_evalue: 10.0,
        };
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[query.len()]);
        let chain = Box::new(crate::spliced_hits::HSPChain {
            context: 0,
            oid: 35,
            score: 35,
            hsps: Some(Box::new(crate::spliced_hits::HSPContainer {
                hsp: mapper_saved_hsp(13, 23, 13, 23),
                next: None,
            })),
            count: 1,
            pair: None,
            pair_score: None,
            pair_conf: crate::spliced_hits::PAIR_NONE,
            adapter: -1,
            poly_a: 0,
            next: None,
        });
        let mapper_data = crate::spliced_hits::BlastHSPMapperData {
            query_info: Some(query_info.clone()),
            saved_chains: vec![Some(chain)],
            ..Default::default()
        };
        let mut jumper = JumperGapAlign {
            left_prelim_block: Some(JumperPrelimEditBlock::default()),
            right_prelim_block: Some(JumperPrelimEditBlock::default()),
            table: Vec::new(),
        };

        assert!(do_anchored_search(
            &query,
            &subject,
            subject_bases.len() as i32,
            35,
            0,
            8,
            hit_params.options.longest_intron,
            Some(&mapper_data),
            &query_info,
            &mut jumper,
            &score_params,
            &hit_params,
            JumperAlignParams {
                max_mismatches: 0,
                mismatch_window: 8,
                gap_x_dropoff: 30,
            },
            None,
        )
        .is_none());
    }

    #[test]
    fn do_anchored_search_to_stream_drops_unrecovered_partial_chain() {
        let query = vec![1; 44];
        let subject_bases = vec![0; query.len()];
        let subject = crate::encoding::pack_ncbi2na_bases(&subject_bases);
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let hit_params = crate::parameters::HitSavingParameters {
            options: crate::options::HitSavingOptions {
                cutoff_score: 1,
                longest_intron: 16,
                ..Default::default()
            },
            cutoff_score_min: 0,
            low_score: Vec::new(),
            cutoffs: Vec::new(),
            link_hsp_params: None,
            prelim_evalue: 10.0,
        };
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[query.len()]);
        let chain = Box::new(crate::spliced_hits::HSPChain {
            context: 0,
            oid: 36,
            score: 35,
            hsps: Some(Box::new(crate::spliced_hits::HSPContainer {
                hsp: mapper_saved_hsp(13, 23, 13, 23),
                next: None,
            })),
            count: 1,
            pair: None,
            pair_score: None,
            pair_conf: crate::spliced_hits::PAIR_NONE,
            adapter: -1,
            poly_a: 0,
            next: None,
        });
        let mapper_data = crate::spliced_hits::BlastHSPMapperData {
            query_info: Some(query_info.clone()),
            saved_chains: vec![Some(chain)],
            ..Default::default()
        };
        let mut jumper = JumperGapAlign {
            left_prelim_block: Some(JumperPrelimEditBlock::default()),
            right_prelim_block: Some(JumperPrelimEditBlock::default()),
            table: Vec::new(),
        };
        let stream = crate::hspstream::blast_hsp_stream_new(1);

        assert_eq!(
            do_anchored_search_to_stream(
                &query,
                &subject,
                subject_bases.len() as i32,
                36,
                0,
                8,
                hit_params.options.longest_intron,
                Some(&mapper_data),
                &query_info,
                &mut jumper,
                &score_params,
                &hit_params,
                JumperAlignParams {
                    max_mismatches: 0,
                    mismatch_window: 8,
                    gap_x_dropoff: 30,
                },
                Some(&stream),
            ),
            0
        );
        let (status, written) = crate::hspstream::blast_hsp_stream_read(Some(&stream));
        assert_eq!(status, crate::hspstream::K_BLAST_HSP_STREAM_EOF);
        assert!(written.is_none());
    }

    #[test]
    fn do_anchored_search_recovers_and_copies_multiple_saved_chains_for_oid() {
        let query: Vec<u8> = (0..52).map(|i| (i % 4) as u8).collect();
        let subject = crate::encoding::pack_ncbi2na_bases(&query);
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let hit_params = crate::parameters::HitSavingParameters {
            options: crate::options::HitSavingOptions {
                cutoff_score: 1,
                longest_intron: 16,
                ..Default::default()
            },
            cutoff_score_min: 0,
            low_score: Vec::new(),
            cutoffs: Vec::new(),
            link_hsp_params: None,
            prelim_evalue: 10.0,
        };
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[query.len()]);
        let second_chain = Box::new(crate::spliced_hits::HSPChain {
            context: 0,
            oid: 29,
            score: 35,
            hsps: Some(Box::new(crate::spliced_hits::HSPContainer {
                hsp: mapper_saved_hsp(24, 34, 24, 34),
                next: None,
            })),
            count: 1,
            pair: None,
            pair_score: None,
            pair_conf: crate::spliced_hits::PAIR_NONE,
            adapter: -1,
            poly_a: 0,
            next: None,
        });
        let first_chain = Box::new(crate::spliced_hits::HSPChain {
            context: 0,
            oid: 29,
            score: 35,
            hsps: Some(Box::new(crate::spliced_hits::HSPContainer {
                hsp: mapper_saved_hsp(13, 23, 13, 23),
                next: None,
            })),
            count: 1,
            pair: None,
            pair_score: None,
            pair_conf: crate::spliced_hits::PAIR_NONE,
            adapter: -1,
            poly_a: 0,
            next: Some(second_chain),
        });
        let mapper_data = crate::spliced_hits::BlastHSPMapperData {
            query_info: Some(query_info.clone()),
            saved_chains: vec![Some(first_chain)],
            ..Default::default()
        };
        let mut jumper = JumperGapAlign {
            left_prelim_block: Some(JumperPrelimEditBlock::default()),
            right_prelim_block: Some(JumperPrelimEditBlock::default()),
            table: Vec::new(),
        };

        let list = do_anchored_search(
            &query,
            &subject,
            query.len() as i32,
            29,
            0,
            12,
            hit_params.options.longest_intron,
            Some(&mapper_data),
            &query_info,
            &mut jumper,
            &score_params,
            &hit_params,
            JumperAlignParams {
                max_mismatches: 5,
                mismatch_window: 10,
                gap_x_dropoff: 30,
            },
            None,
        )
        .expect("anchored search should recover flanks for linked chains");

        assert_eq!(list.oid, 29);
        assert!(list
            .hsps
            .iter()
            .any(|hsp| hsp.query_offset == 13 && hsp.query_end == 23));
        assert!(list
            .hsps
            .iter()
            .any(|hsp| hsp.query_offset == 24 && hsp.query_end == 34));
        assert!(list.hsps.iter().filter(|hsp| hsp.query_offset == 0).count() >= 1);
    }

    #[test]
    fn do_anchored_search_to_stream_writes_recovered_subject_list() {
        let query: Vec<u8> = (0..40).map(|i| (i % 4) as u8).collect();
        let subject = crate::encoding::pack_ncbi2na_bases(&query);
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let hit_params = crate::parameters::HitSavingParameters {
            options: crate::options::HitSavingOptions {
                cutoff_score: 1,
                ..Default::default()
            },
            cutoff_score_min: 0,
            low_score: Vec::new(),
            cutoffs: Vec::new(),
            link_hsp_params: None,
            prelim_evalue: 10.0,
        };
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[query.len()]);
        let chain_hsp = crate::hspstream::Hsp {
            score: 35,
            num_ident: 10,
            bit_score: 0.0,
            evalue: 0.0,
            query_offset: 13,
            query_end: 23,
            query_gapped_start: 13,
            subject_offset: 13,
            subject_end: 23,
            subject_gapped_start: 13,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        };
        let chain = Box::new(crate::spliced_hits::HSPChain {
            context: 0,
            oid: 22,
            score: 35,
            hsps: Some(Box::new(crate::spliced_hits::HSPContainer {
                hsp: chain_hsp,
                next: None,
            })),
            count: 1,
            pair: None,
            pair_score: None,
            pair_conf: crate::spliced_hits::PAIR_NONE,
            adapter: -1,
            poly_a: 0,
            next: None,
        });
        let mapper_data = crate::spliced_hits::BlastHSPMapperData {
            query_info: Some(query_info.clone()),
            saved_chains: vec![Some(chain)],
            ..Default::default()
        };
        let mut jumper = JumperGapAlign {
            left_prelim_block: Some(JumperPrelimEditBlock::default()),
            right_prelim_block: Some(JumperPrelimEditBlock::default()),
            table: Vec::new(),
        };
        let stream = crate::hspstream::blast_hsp_stream_new(1);

        assert_eq!(
            do_anchored_search_to_stream(
                &query,
                &subject,
                query.len() as i32,
                22,
                0,
                12,
                12,
                Some(&mapper_data),
                &query_info,
                &mut jumper,
                &score_params,
                &hit_params,
                JumperAlignParams {
                    max_mismatches: 5,
                    mismatch_window: 10,
                    gap_x_dropoff: 30,
                },
                Some(&stream),
            ),
            0
        );

        let (status, written) = crate::hspstream::blast_hsp_stream_read(Some(&stream));
        assert_eq!(status, crate::hspstream::K_BLAST_HSP_STREAM_SUCCESS);
        let written = written.expect("anchored search should write a list");
        assert_eq!(written.oid, 22);
        assert!(written.hsps.iter().any(|hsp| hsp.query_offset == 0));
        assert!(written.hsps.iter().any(|hsp| hsp.query_offset == 13));
    }

    #[test]
    fn do_anchored_search_to_stream_writes_linked_saved_chains_for_oid() {
        let query: Vec<u8> = (0..52).map(|i| (i % 4) as u8).collect();
        let subject = crate::encoding::pack_ncbi2na_bases(&query);
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let hit_params = crate::parameters::HitSavingParameters {
            options: crate::options::HitSavingOptions {
                cutoff_score: 1,
                longest_intron: 16,
                ..Default::default()
            },
            cutoff_score_min: 0,
            low_score: Vec::new(),
            cutoffs: Vec::new(),
            link_hsp_params: None,
            prelim_evalue: 10.0,
        };
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[query.len()]);
        let second_chain = Box::new(crate::spliced_hits::HSPChain {
            context: 0,
            oid: 30,
            score: 35,
            hsps: Some(Box::new(crate::spliced_hits::HSPContainer {
                hsp: mapper_saved_hsp(24, 34, 24, 34),
                next: None,
            })),
            count: 1,
            pair: None,
            pair_score: None,
            pair_conf: crate::spliced_hits::PAIR_NONE,
            adapter: -1,
            poly_a: 0,
            next: None,
        });
        let first_chain = Box::new(crate::spliced_hits::HSPChain {
            context: 0,
            oid: 30,
            score: 35,
            hsps: Some(Box::new(crate::spliced_hits::HSPContainer {
                hsp: mapper_saved_hsp(13, 23, 13, 23),
                next: None,
            })),
            count: 1,
            pair: None,
            pair_score: None,
            pair_conf: crate::spliced_hits::PAIR_NONE,
            adapter: -1,
            poly_a: 0,
            next: Some(second_chain),
        });
        let mapper_data = crate::spliced_hits::BlastHSPMapperData {
            query_info: Some(query_info.clone()),
            saved_chains: vec![Some(first_chain)],
            ..Default::default()
        };
        let mut jumper = JumperGapAlign {
            left_prelim_block: Some(JumperPrelimEditBlock::default()),
            right_prelim_block: Some(JumperPrelimEditBlock::default()),
            table: Vec::new(),
        };
        let stream = crate::hspstream::blast_hsp_stream_new(1);

        assert_eq!(
            do_anchored_search_to_stream(
                &query,
                &subject,
                query.len() as i32,
                30,
                0,
                12,
                hit_params.options.longest_intron,
                Some(&mapper_data),
                &query_info,
                &mut jumper,
                &score_params,
                &hit_params,
                JumperAlignParams {
                    max_mismatches: 5,
                    mismatch_window: 10,
                    gap_x_dropoff: 30,
                },
                Some(&stream),
            ),
            0
        );

        let (status, written) = crate::hspstream::blast_hsp_stream_read(Some(&stream));
        assert_eq!(status, crate::hspstream::K_BLAST_HSP_STREAM_SUCCESS);
        let written = written.expect("anchored search should write linked-chain list");
        assert_eq!(written.oid, 30);
        assert!(written
            .hsps
            .iter()
            .any(|hsp| hsp.query_offset == 13 && hsp.query_end == 23));
        assert!(written
            .hsps
            .iter()
            .any(|hsp| hsp.query_offset == 24 && hsp.query_end == 34));
        assert!(
            written
                .hsps
                .iter()
                .filter(|hsp| hsp.query_offset == 0)
                .count()
                >= 1
        );
    }

    #[test]
    fn do_anchored_search_to_stream_writes_nonzero_query_index_linked_chains() {
        let query: Vec<u8> = (0..52).map(|i| (i % 4) as u8).collect();
        let subject = crate::encoding::pack_ncbi2na_bases(&query);
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let hit_params = crate::parameters::HitSavingParameters {
            options: crate::options::HitSavingOptions {
                cutoff_score: 1,
                longest_intron: 16,
                ..Default::default()
            },
            cutoff_score_min: 0,
            low_score: Vec::new(),
            cutoffs: Vec::new(),
            link_hsp_params: None,
            prelim_evalue: 10.0,
        };
        let query_info = crate::queryinfo::QueryInfo {
            num_queries: 2,
            contexts: vec![
                crate::queryinfo::ContextInfo {
                    query_offset: 0,
                    query_length: query.len() as i32,
                    eff_searchsp: 0,
                    length_adjustment: 0,
                    query_index: 0,
                    frame: 0,
                    is_valid: true,
                    segment_flags: crate::queryinfo::E_NO_SEGMENTS,
                },
                crate::queryinfo::ContextInfo {
                    query_offset: 0,
                    query_length: query.len() as i32,
                    eff_searchsp: 0,
                    length_adjustment: 0,
                    query_index: 1,
                    frame: 0,
                    is_valid: true,
                    segment_flags: crate::queryinfo::E_NO_SEGMENTS,
                },
            ],
            max_length: query.len() as u32,
            min_length: 0,
        };
        let mut second_hsp = mapper_saved_hsp(24, 34, 24, 34);
        second_hsp.context = 1;
        let second_chain = Box::new(crate::spliced_hits::HSPChain {
            context: 1,
            oid: 32,
            score: 35,
            hsps: Some(Box::new(crate::spliced_hits::HSPContainer {
                hsp: second_hsp,
                next: None,
            })),
            count: 1,
            pair: None,
            pair_score: None,
            pair_conf: crate::spliced_hits::PAIR_NONE,
            adapter: -1,
            poly_a: 0,
            next: None,
        });
        let mut first_hsp = mapper_saved_hsp(13, 23, 13, 23);
        first_hsp.context = 1;
        let first_chain = Box::new(crate::spliced_hits::HSPChain {
            context: 1,
            oid: 32,
            score: 35,
            hsps: Some(Box::new(crate::spliced_hits::HSPContainer {
                hsp: first_hsp,
                next: None,
            })),
            count: 1,
            pair: None,
            pair_score: None,
            pair_conf: crate::spliced_hits::PAIR_NONE,
            adapter: -1,
            poly_a: 0,
            next: Some(second_chain),
        });
        let mapper_data = crate::spliced_hits::BlastHSPMapperData {
            query_info: Some(query_info.clone()),
            saved_chains: vec![None, Some(first_chain)],
            ..Default::default()
        };
        let mut jumper = JumperGapAlign {
            left_prelim_block: Some(JumperPrelimEditBlock::default()),
            right_prelim_block: Some(JumperPrelimEditBlock::default()),
            table: Vec::new(),
        };
        let stream = crate::hspstream::blast_hsp_stream_new(2);

        assert_eq!(
            do_anchored_search_to_stream(
                &query,
                &subject,
                query.len() as i32,
                32,
                0,
                12,
                hit_params.options.longest_intron,
                Some(&mapper_data),
                &query_info,
                &mut jumper,
                &score_params,
                &hit_params,
                JumperAlignParams {
                    max_mismatches: 5,
                    mismatch_window: 10,
                    gap_x_dropoff: 30,
                },
                Some(&stream),
            ),
            0
        );

        let (status, written) = crate::hspstream::blast_hsp_stream_read(Some(&stream));
        assert_eq!(status, crate::hspstream::K_BLAST_HSP_STREAM_SUCCESS);
        let written = written.expect("anchored search should write nonzero-query list");
        assert_eq!(written.oid, 32);
        assert!(written
            .hsps
            .iter()
            .any(|hsp| hsp.context == 1 && hsp.query_offset == 13 && hsp.query_end == 23));
        assert!(written
            .hsps
            .iter()
            .any(|hsp| hsp.context == 1 && hsp.query_offset == 24 && hsp.query_end == 34));
        assert!(
            written
                .hsps
                .iter()
                .filter(|hsp| hsp.context == 1 && hsp.query_offset == 0)
                .count()
                >= 1
        );
    }

    #[test]
    fn do_anchored_search_to_stream_filters_low_score_linked_chains() {
        let query: Vec<u8> = (0..52).map(|i| (i % 4) as u8).collect();
        let subject = crate::encoding::pack_ncbi2na_bases(&query);
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let hit_params = crate::parameters::HitSavingParameters {
            options: crate::options::HitSavingOptions {
                cutoff_score: 1,
                longest_intron: 16,
                ..Default::default()
            },
            cutoff_score_min: 0,
            low_score: Vec::new(),
            cutoffs: Vec::new(),
            link_hsp_params: None,
            prelim_evalue: 10.0,
        };
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[query.len()]);
        let high_score_chain = Box::new(crate::spliced_hits::HSPChain {
            context: 0,
            oid: 33,
            score: 35,
            hsps: Some(Box::new(crate::spliced_hits::HSPContainer {
                hsp: mapper_saved_hsp(24, 34, 24, 34),
                next: None,
            })),
            count: 1,
            pair: None,
            pair_score: None,
            pair_conf: crate::spliced_hits::PAIR_NONE,
            adapter: -1,
            poly_a: 0,
            next: None,
        });
        let low_score_chain = Box::new(crate::spliced_hits::HSPChain {
            context: 0,
            oid: 33,
            score: 29,
            hsps: Some(Box::new(crate::spliced_hits::HSPContainer {
                hsp: mapper_saved_hsp(13, 23, 13, 23),
                next: None,
            })),
            count: 1,
            pair: None,
            pair_score: None,
            pair_conf: crate::spliced_hits::PAIR_NONE,
            adapter: -1,
            poly_a: 0,
            next: Some(high_score_chain),
        });
        let mapper_data = crate::spliced_hits::BlastHSPMapperData {
            query_info: Some(query_info.clone()),
            saved_chains: vec![Some(low_score_chain)],
            ..Default::default()
        };
        let mut jumper = JumperGapAlign {
            left_prelim_block: Some(JumperPrelimEditBlock::default()),
            right_prelim_block: Some(JumperPrelimEditBlock::default()),
            table: Vec::new(),
        };
        let stream = crate::hspstream::blast_hsp_stream_new(1);

        assert_eq!(
            do_anchored_search_to_stream(
                &query,
                &subject,
                query.len() as i32,
                33,
                0,
                12,
                hit_params.options.longest_intron,
                Some(&mapper_data),
                &query_info,
                &mut jumper,
                &score_params,
                &hit_params,
                JumperAlignParams {
                    max_mismatches: 5,
                    mismatch_window: 10,
                    gap_x_dropoff: 30,
                },
                Some(&stream),
            ),
            0
        );

        let (status, written) = crate::hspstream::blast_hsp_stream_read(Some(&stream));
        assert_eq!(status, crate::hspstream::K_BLAST_HSP_STREAM_SUCCESS);
        let written = written.expect("anchored search should write eligible linked chain");
        assert_eq!(written.oid, 33);
        assert!(!written
            .hsps
            .iter()
            .any(|hsp| hsp.query_offset == 13 && hsp.query_end == 23));
        assert!(written
            .hsps
            .iter()
            .any(|hsp| hsp.query_offset == 24 && hsp.query_end == 34));
    }

    #[test]
    fn do_anchored_search_to_stream_reports_linked_chain_query_index_insert_error() {
        let query: Vec<u8> = (0..52).map(|i| (i % 4) as u8).collect();
        let subject = crate::encoding::pack_ncbi2na_bases(&query);
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let hit_params = crate::parameters::HitSavingParameters {
            options: crate::options::HitSavingOptions {
                cutoff_score: 1,
                longest_intron: 16,
                ..Default::default()
            },
            cutoff_score_min: 0,
            low_score: Vec::new(),
            cutoffs: Vec::new(),
            link_hsp_params: None,
            prelim_evalue: 10.0,
        };
        let query_info = crate::queryinfo::QueryInfo {
            num_queries: 2,
            contexts: vec![
                crate::queryinfo::ContextInfo {
                    query_offset: 0,
                    query_length: query.len() as i32,
                    eff_searchsp: 0,
                    length_adjustment: 0,
                    query_index: 0,
                    frame: 0,
                    is_valid: true,
                    segment_flags: crate::queryinfo::E_NO_SEGMENTS,
                },
                crate::queryinfo::ContextInfo {
                    query_offset: 0,
                    query_length: query.len() as i32,
                    eff_searchsp: 0,
                    length_adjustment: 0,
                    query_index: 1,
                    frame: 0,
                    is_valid: true,
                    segment_flags: crate::queryinfo::E_NO_SEGMENTS,
                },
            ],
            max_length: query.len() as u32,
            min_length: 0,
        };
        let mut second_hsp = mapper_saved_hsp(24, 34, 24, 34);
        second_hsp.context = 1;
        let second_chain = Box::new(crate::spliced_hits::HSPChain {
            context: 1,
            oid: 31,
            score: 35,
            hsps: Some(Box::new(crate::spliced_hits::HSPContainer {
                hsp: second_hsp,
                next: None,
            })),
            count: 1,
            pair: None,
            pair_score: None,
            pair_conf: crate::spliced_hits::PAIR_NONE,
            adapter: -1,
            poly_a: 0,
            next: None,
        });
        let mut first_hsp = mapper_saved_hsp(13, 23, 13, 23);
        first_hsp.context = 1;
        let first_chain = Box::new(crate::spliced_hits::HSPChain {
            context: 1,
            oid: 31,
            score: 35,
            hsps: Some(Box::new(crate::spliced_hits::HSPContainer {
                hsp: first_hsp,
                next: None,
            })),
            count: 1,
            pair: None,
            pair_score: None,
            pair_conf: crate::spliced_hits::PAIR_NONE,
            adapter: -1,
            poly_a: 0,
            next: Some(second_chain),
        });
        let mapper_data = crate::spliced_hits::BlastHSPMapperData {
            query_info: Some(query_info.clone()),
            saved_chains: vec![None, Some(first_chain)],
            ..Default::default()
        };
        let mut jumper = JumperGapAlign {
            left_prelim_block: Some(JumperPrelimEditBlock::default()),
            right_prelim_block: Some(JumperPrelimEditBlock::default()),
            table: Vec::new(),
        };
        let stream = crate::hspstream::blast_hsp_stream_new(1);

        assert_eq!(
            do_anchored_search_to_stream(
                &query,
                &subject,
                query.len() as i32,
                31,
                0,
                12,
                hit_params.options.longest_intron,
                Some(&mapper_data),
                &query_info,
                &mut jumper,
                &score_params,
                &hit_params,
                JumperAlignParams {
                    max_mismatches: 5,
                    mismatch_window: 10,
                    gap_x_dropoff: 30,
                },
                Some(&stream),
            ),
            1
        );
        let (status, written) = crate::hspstream::blast_hsp_stream_read(Some(&stream));
        assert_eq!(status, crate::hspstream::K_BLAST_HSP_STREAM_EOF);
        assert!(written.is_none());
    }

    #[test]
    fn do_anchored_search_to_stream_matches_c_empty_success_and_missing_stream_error() {
        let query: Vec<u8> = (0..40).map(|i| (i % 4) as u8).collect();
        let subject = crate::encoding::pack_ncbi2na_bases(&query);
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let hit_params = crate::parameters::HitSavingParameters {
            options: crate::options::HitSavingOptions {
                cutoff_score: 1,
                ..Default::default()
            },
            cutoff_score_min: 0,
            low_score: Vec::new(),
            cutoffs: Vec::new(),
            link_hsp_params: None,
            prelim_evalue: 10.0,
        };
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[query.len()]);
        let mapper_data = crate::spliced_hits::BlastHSPMapperData {
            query_info: Some(query_info.clone()),
            saved_chains: Vec::new(),
            ..Default::default()
        };
        let mut jumper = JumperGapAlign {
            left_prelim_block: Some(JumperPrelimEditBlock::default()),
            right_prelim_block: Some(JumperPrelimEditBlock::default()),
            table: Vec::new(),
        };

        assert_eq!(
            do_anchored_search_to_stream(
                &query,
                &subject,
                query.len() as i32,
                23,
                0,
                12,
                12,
                Some(&mapper_data),
                &query_info,
                &mut jumper,
                &score_params,
                &hit_params,
                JumperAlignParams {
                    max_mismatches: 5,
                    mismatch_window: 10,
                    gap_x_dropoff: 30,
                },
                None,
            ),
            -1
        );

        let stream = crate::hspstream::blast_hsp_stream_new(1);
        assert_eq!(
            do_anchored_search_to_stream(
                &query,
                &subject,
                query.len() as i32,
                23,
                0,
                12,
                12,
                Some(&mapper_data),
                &query_info,
                &mut jumper,
                &score_params,
                &hit_params,
                JumperAlignParams {
                    max_mismatches: 5,
                    mismatch_window: 10,
                    gap_x_dropoff: 30,
                },
                Some(&stream),
            ),
            0
        );
        let (status, written) = crate::hspstream::blast_hsp_stream_read(Some(&stream));
        assert_eq!(status, crate::hspstream::K_BLAST_HSP_STREAM_EOF);
        assert!(written.is_none());

        let port_stream = crate::hspstream::blast_hsp_stream_new(1);
        assert!(do_anchored_search(
            &query,
            &subject,
            query.len() as i32,
            23,
            0,
            12,
            12,
            Some(&mapper_data),
            &query_info,
            &mut jumper,
            &score_params,
            &hit_params,
            JumperAlignParams {
                max_mismatches: 5,
                mismatch_window: 10,
                gap_x_dropoff: 30,
            },
            Some(&port_stream),
        )
        .is_none());
        let (status, written) = crate::hspstream::blast_hsp_stream_read(Some(&port_stream));
        assert_eq!(status, crate::hspstream::K_BLAST_HSP_STREAM_EOF);
        assert!(written.is_none());

        let align_params = JumperAlignParams {
            max_mismatches: 5,
            mismatch_window: 10,
            gap_x_dropoff: 30,
        };
        let out_of_range_query_info = crate::queryinfo::QueryInfo {
            num_queries: 1,
            contexts: vec![crate::queryinfo::ContextInfo {
                query_offset: 100,
                query_length: query.len() as i32,
                eff_searchsp: 0,
                length_adjustment: 0,
                query_index: 0,
                frame: 0,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: query.len() as u32,
            min_length: 0,
        };
        let out_of_range_chain = Box::new(crate::spliced_hits::HSPChain {
            context: 0,
            oid: 26,
            score: 35,
            hsps: Some(Box::new(crate::spliced_hits::HSPContainer {
                hsp: mapper_saved_hsp(13, 23, 13, 23),
                next: None,
            })),
            count: 1,
            pair: None,
            pair_score: None,
            pair_conf: crate::spliced_hits::PAIR_NONE,
            adapter: -1,
            poly_a: 0,
            next: None,
        });
        let out_of_range_mapper_data = crate::spliced_hits::BlastHSPMapperData {
            query_info: Some(out_of_range_query_info.clone()),
            saved_chains: vec![Some(out_of_range_chain)],
            ..Default::default()
        };
        let out_of_range_stream = crate::hspstream::blast_hsp_stream_new(1);

        assert!(do_anchored_search(
            &query,
            &subject,
            query.len() as i32,
            26,
            0,
            12,
            12,
            Some(&out_of_range_mapper_data),
            &out_of_range_query_info,
            &mut jumper,
            &score_params,
            &hit_params,
            align_params,
            None,
        )
        .is_none());
        assert_eq!(
            do_anchored_search_to_stream(
                &query,
                &subject,
                query.len() as i32,
                26,
                0,
                12,
                12,
                Some(&out_of_range_mapper_data),
                &out_of_range_query_info,
                &mut jumper,
                &score_params,
                &hit_params,
                align_params,
                Some(&out_of_range_stream),
            ),
            0
        );
        let (status, written) = crate::hspstream::blast_hsp_stream_read(Some(&out_of_range_stream));
        assert_eq!(status, crate::hspstream::K_BLAST_HSP_STREAM_EOF);
        assert!(written.is_none());
    }

    #[test]
    fn do_anchored_search_rejects_missing_mapper_data_like_c_state() {
        let query: Vec<u8> = (0..40).map(|i| (i % 4) as u8).collect();
        let subject = crate::encoding::pack_ncbi2na_bases(&query);
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let hit_params = crate::parameters::HitSavingParameters {
            options: crate::options::HitSavingOptions {
                cutoff_score: 1,
                ..Default::default()
            },
            cutoff_score_min: 0,
            low_score: Vec::new(),
            cutoffs: Vec::new(),
            link_hsp_params: None,
            prelim_evalue: 10.0,
        };
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[query.len()]);
        let mut jumper = JumperGapAlign {
            left_prelim_block: Some(JumperPrelimEditBlock::default()),
            right_prelim_block: Some(JumperPrelimEditBlock::default()),
            table: Vec::new(),
        };
        let align_params = JumperAlignParams {
            max_mismatches: 5,
            mismatch_window: 10,
            gap_x_dropoff: 30,
        };
        let stream = crate::hspstream::blast_hsp_stream_new(1);

        assert!(do_anchored_search(
            &query,
            &subject,
            query.len() as i32,
            24,
            0,
            12,
            12,
            None,
            &query_info,
            &mut jumper,
            &score_params,
            &hit_params,
            align_params,
            Some(&stream),
        )
        .is_none());
        assert_eq!(
            do_anchored_search_to_stream(
                &query,
                &subject,
                query.len() as i32,
                24,
                0,
                12,
                12,
                None,
                &query_info,
                &mut jumper,
                &score_params,
                &hit_params,
                align_params,
                Some(&stream),
            ),
            -1
        );
        let (status, written) = crate::hspstream::blast_hsp_stream_read(Some(&stream));
        assert_eq!(status, crate::hspstream::K_BLAST_HSP_STREAM_EOF);
        assert!(written.is_none());

        let mismatch_subject = crate::encoding::pack_ncbi2na_bases(&vec![3u8; query.len()]);
        let mut huge_left_hsp = mapper_saved_hsp(13, 23, 13, 23);
        huge_left_hsp.context = 0;
        let huge_left_chain = Box::new(crate::spliced_hits::HSPChain {
            context: 0,
            oid: 25,
            score: 35,
            hsps: Some(Box::new(crate::spliced_hits::HSPContainer {
                hsp: huge_left_hsp,
                next: None,
            })),
            count: 1,
            pair: None,
            pair_score: None,
            pair_conf: crate::spliced_hits::PAIR_NONE,
            adapter: -1,
            poly_a: 0,
            next: None,
        });
        let huge_left_mapper_data = crate::spliced_hits::BlastHSPMapperData {
            query_info: Some(query_info.clone()),
            saved_chains: vec![Some(huge_left_chain)],
            ..Default::default()
        };
        assert!(do_anchored_search(
            &query,
            &mismatch_subject,
            query.len() as i32,
            25,
            0,
            12,
            i32::MIN,
            Some(&huge_left_mapper_data),
            &query_info,
            &mut jumper,
            &score_params,
            &hit_params,
            align_params,
            None,
        )
        .is_none());

        let huge_right_hsp = mapper_saved_hsp(0, 4, i32::MAX - 4, i32::MAX);
        let huge_right_chain = Box::new(crate::spliced_hits::HSPChain {
            context: 0,
            oid: 25,
            score: 35,
            hsps: Some(Box::new(crate::spliced_hits::HSPContainer {
                hsp: huge_right_hsp,
                next: None,
            })),
            count: 1,
            pair: None,
            pair_score: None,
            pair_conf: crate::spliced_hits::PAIR_NONE,
            adapter: -1,
            poly_a: 0,
            next: None,
        });
        let huge_right_mapper_data = crate::spliced_hits::BlastHSPMapperData {
            query_info: Some(query_info.clone()),
            saved_chains: vec![Some(huge_right_chain)],
            ..Default::default()
        };
        assert!(do_anchored_search(
            &query,
            &mismatch_subject,
            query.len() as i32,
            25,
            0,
            12,
            1,
            Some(&huge_right_mapper_data),
            &query_info,
            &mut jumper,
            &score_params,
            &hit_params,
            align_params,
            None,
        )
        .is_none());

        let huge_negative_subject_hsp = mapper_saved_hsp(13, 23, i32::MIN, i32::MIN + 10);
        let huge_negative_subject_chain = Box::new(crate::spliced_hits::HSPChain {
            context: 0,
            oid: 25,
            score: 35,
            hsps: Some(Box::new(crate::spliced_hits::HSPContainer {
                hsp: huge_negative_subject_hsp,
                next: None,
            })),
            count: 1,
            pair: None,
            pair_score: None,
            pair_conf: crate::spliced_hits::PAIR_NONE,
            adapter: -1,
            poly_a: 0,
            next: None,
        });
        let huge_negative_subject_mapper_data = crate::spliced_hits::BlastHSPMapperData {
            query_info: Some(query_info.clone()),
            saved_chains: vec![Some(huge_negative_subject_chain)],
            ..Default::default()
        };
        assert!(do_anchored_search(
            &query,
            &mismatch_subject,
            query.len() as i32,
            25,
            0,
            12,
            i32::MAX,
            Some(&huge_negative_subject_mapper_data),
            &query_info,
            &mut jumper,
            &score_params,
            &hit_params,
            align_params,
            None,
        )
        .is_none());
    }

    #[test]
    fn do_anchored_search_skips_empty_saved_chains() {
        let query: Vec<u8> = (0..40).map(|i| (i % 4) as u8).collect();
        let subject = crate::encoding::pack_ncbi2na_bases(&query);
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let hit_params = crate::parameters::HitSavingParameters {
            options: crate::options::HitSavingOptions {
                cutoff_score: 1,
                ..Default::default()
            },
            cutoff_score_min: 0,
            low_score: Vec::new(),
            cutoffs: Vec::new(),
            link_hsp_params: None,
            prelim_evalue: 10.0,
        };
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[query.len()]);
        let empty_chain = Box::new(crate::spliced_hits::HSPChain {
            context: 0,
            oid: 27,
            score: 35,
            hsps: None,
            count: 0,
            pair: None,
            pair_score: None,
            pair_conf: crate::spliced_hits::PAIR_NONE,
            adapter: -1,
            poly_a: 0,
            next: None,
        });
        let mapper_data = crate::spliced_hits::BlastHSPMapperData {
            query_info: Some(query_info.clone()),
            saved_chains: vec![Some(empty_chain)],
            ..Default::default()
        };
        let mut jumper = JumperGapAlign {
            left_prelim_block: Some(JumperPrelimEditBlock::default()),
            right_prelim_block: Some(JumperPrelimEditBlock::default()),
            table: Vec::new(),
        };
        let align_params = JumperAlignParams {
            max_mismatches: 5,
            mismatch_window: 10,
            gap_x_dropoff: 30,
        };

        assert!(do_anchored_search(
            &query,
            &subject,
            query.len() as i32,
            27,
            0,
            12,
            12,
            Some(&mapper_data),
            &query_info,
            &mut jumper,
            &score_params,
            &hit_params,
            align_params,
            None,
        )
        .is_none());

        let stream = crate::hspstream::blast_hsp_stream_new(1);
        assert_eq!(
            do_anchored_search_to_stream(
                &query,
                &subject,
                query.len() as i32,
                27,
                0,
                12,
                12,
                Some(&mapper_data),
                &query_info,
                &mut jumper,
                &score_params,
                &hit_params,
                align_params,
                Some(&stream),
            ),
            0
        );
        let (status, written) = crate::hspstream::blast_hsp_stream_read(Some(&stream));
        assert_eq!(status, crate::hspstream::K_BLAST_HSP_STREAM_EOF);
        assert!(written.is_none());
    }

    #[test]
    fn anchored_search_query_index_selection_defaults_to_zero() {
        let query_info = crate::queryinfo::QueryInfo {
            num_queries: 2,
            contexts: vec![
                crate::queryinfo::ContextInfo {
                    query_offset: 0,
                    query_length: 10,
                    eff_searchsp: 0,
                    length_adjustment: 0,
                    query_index: 3,
                    frame: 0,
                    is_valid: true,
                    segment_flags: crate::queryinfo::E_NO_SEGMENTS,
                },
                crate::queryinfo::ContextInfo {
                    query_offset: 10,
                    query_length: 10,
                    eff_searchsp: 0,
                    length_adjustment: 0,
                    query_index: 4,
                    frame: 0,
                    is_valid: true,
                    segment_flags: crate::queryinfo::E_NO_SEGMENTS,
                },
            ],
            max_length: 10,
            min_length: 0,
        };
        let mut list = crate::hspstream::HspList::new(1);
        assert_eq!(query_index_for_hsp_list(&list, &query_info), 0);

        let mut hsp = mapper_saved_hsp(0, 4, 0, 4);
        hsp.context = 1;
        list.add_hsp(hsp);
        assert_eq!(query_index_for_hsp_list(&list, &query_info), 4);

        list.hsps[0].context = -3;
        assert_eq!(query_index_for_hsp_list(&list, &query_info), 3);

        list.hsps[0].context = 9;
        assert_eq!(query_index_for_hsp_list(&list, &query_info), 0);
    }

    #[test]
    fn do_anchored_search_skips_selected_chains_with_unusable_query_contexts() {
        let query: Vec<u8> = (0..40).map(|i| (i % 4) as u8).collect();
        let subject = crate::encoding::pack_ncbi2na_bases(&query);
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let hit_params = crate::parameters::HitSavingParameters {
            options: crate::options::HitSavingOptions {
                cutoff_score: 1,
                ..Default::default()
            },
            cutoff_score_min: 0,
            low_score: Vec::new(),
            cutoffs: Vec::new(),
            link_hsp_params: None,
            prelim_evalue: 10.0,
        };
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[query.len()]);
        let mut missing_context_hsp = mapper_saved_hsp(13, 23, 13, 23);
        missing_context_hsp.context = 4;
        let missing_context_chain = Box::new(crate::spliced_hits::HSPChain {
            context: 4,
            oid: 25,
            score: 35,
            hsps: Some(Box::new(crate::spliced_hits::HSPContainer {
                hsp: missing_context_hsp,
                next: None,
            })),
            count: 1,
            pair: None,
            pair_score: None,
            pair_conf: crate::spliced_hits::PAIR_NONE,
            adapter: -1,
            poly_a: 0,
            next: None,
        });
        let mapper_data = crate::spliced_hits::BlastHSPMapperData {
            query_info: Some(query_info.clone()),
            saved_chains: vec![Some(missing_context_chain)],
            ..Default::default()
        };
        let mut jumper = JumperGapAlign {
            left_prelim_block: Some(JumperPrelimEditBlock::default()),
            right_prelim_block: Some(JumperPrelimEditBlock::default()),
            table: Vec::new(),
        };
        let align_params = JumperAlignParams {
            max_mismatches: 5,
            mismatch_window: 10,
            gap_x_dropoff: 30,
        };

        assert!(do_anchored_search(
            &query,
            &subject,
            query.len() as i32,
            25,
            0,
            12,
            12,
            Some(&mapper_data),
            &query_info,
            &mut jumper,
            &score_params,
            &hit_params,
            align_params,
            None,
        )
        .is_none());

        let stream = crate::hspstream::blast_hsp_stream_new(1);
        assert_eq!(
            do_anchored_search_to_stream(
                &query,
                &subject,
                query.len() as i32,
                25,
                0,
                12,
                12,
                Some(&mapper_data),
                &query_info,
                &mut jumper,
                &score_params,
                &hit_params,
                align_params,
                Some(&stream),
            ),
            0
        );
        let (status, written) = crate::hspstream::blast_hsp_stream_read(Some(&stream));
        assert_eq!(status, crate::hspstream::K_BLAST_HSP_STREAM_EOF);
        assert!(written.is_none());
    }

    struct SyntheticMegablastIndexSource {
        start_oid: u32,
        stop_oid: u32,
        hash_table: Vec<u32>,
        offsets: Vec<u32>,
    }

    impl MegablastIndexedOffsetSource for SyntheticMegablastIndexSource {
        fn hkey_width(&self) -> u32 {
            2
        }

        fn start_oid(&self) -> u32 {
            self.start_oid
        }

        fn stop_oid(&self) -> u32 {
            self.stop_oid
        }

        fn offset_list_words(&self, hash_key: u32) -> Option<&[u32]> {
            let start = *self.hash_table.get(hash_key as usize)? as usize;
            let rel_end = self.offsets[start..].iter().position(|&word| word == 0)?;
            Some(&self.offsets[start..start + rel_end])
        }

        fn logical_chunk_id_for_subject_chunk(&self, local_oid: u32, chunk: u32) -> Option<u32> {
            [1u32, 2, 0]
                .get(local_oid as usize)
                .copied()
                .and_then(|first_chunk| first_chunk.checked_add(chunk))
                .filter(|&logical_chunk_id| (1..=2).contains(&logical_chunk_id))
        }

        fn chunk_offset_for_encoded_offset(&self, encoded_offset: u32) -> Option<(u32, u32)> {
            match encoded_offset {
                64 => Some((1, 64)),
                80 => Some((2, 64)),
                _ => None,
            }
        }

        fn chunk_seed_offset_for_encoded_offset(&self, encoded_offset: u32) -> Option<(u32, u32)> {
            match encoded_offset {
                64 => Some((1, 0)),
                80 => Some((2, 0)),
                _ => None,
            }
        }

        fn min_encoded_offset(&self) -> Option<u32> {
            Some(64)
        }
    }

    fn synthetic_megablast_index_volume_for_callback() -> SyntheticMegablastIndexSource {
        let mut hash_table = vec![0u32; (1usize << 4) + 1];
        hash_table[6] = 1;
        hash_table[16] = 4;
        SyntheticMegablastIndexSource {
            start_oid: 10,
            stop_oid: 12,
            hash_table,
            offsets: vec![0, 3, 64, 0],
        }
    }

    #[test]
    fn megablast_index_adapter_feeds_mb_indexed_callback_shape() {
        let volume = synthetic_megablast_index_volume_for_callback();
        let query = vec![1u8, 2, 3, 0];
        let query_keys = megablast_index_query_keys(&query, volume.hkey_width());
        let subject_bases: Vec<u8> = (0..120).map(|i| ((i + 1) % 4) as u8).collect();
        let subject = crate::encoding::pack_ncbi2na_bases(&subject_bases);
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let word_params = crate::parameters::InitialWordParameters {
            options: crate::options::InitialWordOptions::new_blastn(),
            x_dropoff_max: 8,
            cutoff_score_min: 1,
            cutoffs: Vec::new(),
            ungapped_extension: false,
            nucl_score_table: crate::parameters::InitialWordParameters::build_nucl_score_table(
                score_params.reward,
                score_params.penalty,
            ),
        };
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[query.len()]);
        let mut init_hitlist = crate::extend::InitHitList::new();

        assert_eq!(
            mb_indexed_word_finder(
                &subject,
                subject_bases.len() as i32,
                10,
                0,
                &query,
                &query_info,
                &word_params,
                &score_params,
                &mut init_hitlist,
                |oid| megablast_index_volume_has_oid(&volume, oid),
                |oid, chunk, list| {
                    megablast_index_fill_init_hitlist(&volume, oid, chunk, &query_keys, list)
                },
                None,
            ),
            0
        );

        assert_eq!(init_hitlist.total(), 1);
        assert_eq!(init_hitlist.hits[0].query_offset, 0);
        assert_eq!(init_hitlist.hits[0].subject_offset, 0);
        assert_eq!(
            megablast_index_fill_init_hitlist(&volume, 11, 0, &query_keys, &mut init_hitlist),
            0
        );
        assert_eq!(init_hitlist.total(), 0);
    }

    #[test]
    fn megablast_index_result_holder_routes_oids_across_volumes() {
        let mut first = synthetic_megablast_index_volume_for_callback();
        first.start_oid = 10;
        first.stop_oid = 12;
        let mut second = synthetic_megablast_index_volume_for_callback();
        second.start_oid = 12;
        second.stop_oid = 14;
        let query = vec![1u8, 2, 3, 0];
        let holder = MegablastIndexResultHolder::new(&query, [&first, &second]);
        let subject_bases: Vec<u8> = (0..120).map(|i| ((i + 1) % 4) as u8).collect();
        let subject = crate::encoding::pack_ncbi2na_bases(&subject_bases);
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let word_params = crate::parameters::InitialWordParameters {
            options: crate::options::InitialWordOptions::new_blastn(),
            x_dropoff_max: 8,
            cutoff_score_min: 1,
            cutoffs: Vec::new(),
            ungapped_extension: false,
            nucl_score_table: crate::parameters::InitialWordParameters::build_nucl_score_table(
                score_params.reward,
                score_params.penalty,
            ),
        };
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[query.len()]);
        let mut init_hitlist = crate::extend::InitHitList::new();

        assert!(holder.has_oid(10));
        assert!(holder.has_oid(12));
        assert!(!holder.has_oid(14));
        assert_eq!(
            mb_indexed_word_finder(
                &subject,
                subject_bases.len() as i32,
                12,
                0,
                &query,
                &query_info,
                &word_params,
                &score_params,
                &mut init_hitlist,
                |oid| holder.has_oid(oid),
                |oid, chunk, list| holder.fill_init_hitlist(oid, chunk, list),
                None,
            ),
            0
        );

        assert_eq!(init_hitlist.total(), 1);
        assert_eq!(init_hitlist.hits[0].query_offset, 0);
        assert_eq!(init_hitlist.hits[0].subject_offset, 0);
        assert_eq!(holder.fill_init_hitlist(14, 0, &mut init_hitlist), 0);
        assert_eq!(init_hitlist.total(), 0);
    }

    #[test]
    fn megablast_index_query_keys_honor_unmasked_locations() {
        let query = vec![0u8, 1, 2, 3, 0, 1, 4, 2, 3];
        let locations = [
            crate::util::SSeqRange { left: 1, right: 3 },
            crate::util::SSeqRange { left: 5, right: 8 },
        ];

        let keys = megablast_index_query_keys_for_locations(&query, 2, &locations);

        assert_eq!(
            keys,
            vec![
                MegablastIndexQueryKey {
                    hash_key: 0b0110,
                    query_offset: 1,
                },
                MegablastIndexQueryKey {
                    hash_key: 0b1011,
                    query_offset: 2,
                },
                MegablastIndexQueryKey {
                    hash_key: 0b1011,
                    query_offset: 7,
                },
            ]
        );
        assert!(!keys.iter().any(|key| key.query_offset == 3));
        assert!(!keys.iter().any(|key| key.query_offset == 5));
        assert!(!keys.iter().any(|key| key.query_offset == 6));
    }

    #[test]
    fn mb_indexed_word_finder_filters_redundant_diagonal_hits() {
        let query: Vec<u8> = (0..24).map(|i| (i % 4) as u8).collect();
        let subject = crate::encoding::pack_ncbi2na_bases(&query);
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let mut word_params = crate::parameters::InitialWordParameters {
            options: crate::options::InitialWordOptions::new_blastn(),
            x_dropoff_max: 8,
            cutoff_score_min: 1,
            cutoffs: vec![crate::parameters::BlastUngappedCutoffs {
                x_dropoff_init: 8,
                cutoff_score: 1,
            }],
            ungapped_extension: true,
            nucl_score_table: crate::parameters::InitialWordParameters::build_nucl_score_table(
                score_params.reward,
                score_params.penalty,
            ),
        };
        word_params.options.word_size = 4;
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[query.len()]);
        let mut init_hitlist = crate::extend::InitHitList::new();
        let mut stats = crate::diagnostics::UngappedStats::default();

        assert_eq!(
            mb_indexed_word_finder(
                &subject,
                query.len() as i32,
                5,
                0,
                &query,
                &query_info,
                &word_params,
                &score_params,
                &mut init_hitlist,
                |oid| oid == 5,
                |_oid, _chunk, list| {
                    let _ = crate::extend::blast_save_initial_hit(list, 0, 0, None);
                    let _ = crate::extend::blast_save_initial_hit(list, 1, 1, None);
                    4
                },
                Some(&mut stats),
            ),
            0
        );

        assert_eq!(stats.lookup_hits, 2);
        assert_eq!(stats.init_extends, 1);
        assert_eq!(init_hitlist.total(), 1);
        assert!(init_hitlist.hits[0].ungapped_data.is_some());
    }

    #[test]
    fn mb_indexed_word_finder_drops_low_scoring_extensions_and_updates_stats() {
        let query: Vec<u8> = (0..24).map(|i| (i % 4) as u8).collect();
        let subject = crate::encoding::pack_ncbi2na_bases(&query);
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let word_params = crate::parameters::InitialWordParameters {
            options: crate::options::InitialWordOptions::new_blastn(),
            x_dropoff_max: 8,
            cutoff_score_min: 1000,
            cutoffs: vec![crate::parameters::BlastUngappedCutoffs {
                x_dropoff_init: 8,
                cutoff_score: 1000,
            }],
            ungapped_extension: true,
            nucl_score_table: crate::parameters::InitialWordParameters::build_nucl_score_table(
                score_params.reward,
                score_params.penalty,
            ),
        };
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[query.len()]);
        let mut init_hitlist = crate::extend::InitHitList::new();
        let mut stats = crate::diagnostics::UngappedStats::default();

        assert_eq!(
            mb_indexed_word_finder(
                &subject,
                query.len() as i32,
                5,
                0,
                &query,
                &query_info,
                &word_params,
                &score_params,
                &mut init_hitlist,
                |oid| oid == 5,
                |_oid, _chunk, list| {
                    let _ = crate::extend::blast_save_initial_hit(list, 0, 0, None);
                    4
                },
                Some(&mut stats),
            ),
            0
        );

        assert_eq!(init_hitlist.total(), 0);
        assert_eq!(stats.lookup_hits, 1);
        assert_eq!(stats.init_extends, 1);
        assert_eq!(stats.good_init_extends, 0);
        assert_eq!(stats.num_seqs_passed, 0);
    }

    #[test]
    fn mb_indexed_word_finder_zero_index_hits_leaves_stats_clear() {
        let query: Vec<u8> = (0..16).map(|i| (i % 4) as u8).collect();
        let subject = crate::encoding::pack_ncbi2na_bases(&query);
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let word_params = crate::parameters::InitialWordParameters {
            options: crate::options::InitialWordOptions::new_blastn(),
            x_dropoff_max: 8,
            cutoff_score_min: 1,
            cutoffs: Vec::new(),
            ungapped_extension: true,
            nucl_score_table: crate::parameters::InitialWordParameters::build_nucl_score_table(
                score_params.reward,
                score_params.penalty,
            ),
        };
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[query.len()]);
        let mut init_hitlist = crate::extend::InitHitList::new();
        let mut stats = crate::diagnostics::UngappedStats::default();

        assert_eq!(
            mb_indexed_word_finder(
                &subject,
                query.len() as i32,
                5,
                0,
                &query,
                &query_info,
                &word_params,
                &score_params,
                &mut init_hitlist,
                |oid| oid == 5,
                |_oid, _chunk, _list| 4,
                Some(&mut stats),
            ),
            0
        );

        assert_eq!(init_hitlist.total(), 0);
        assert_eq!(stats.lookup_hits, 0);
        assert_eq!(stats.init_extends, 0);
        assert_eq!(stats.good_init_extends, 0);
    }

    #[test]
    fn mb_indexed_word_finder_preserves_hits_when_ungapped_extension_disabled() {
        let query: Vec<u8> = (0..24).map(|i| (i % 4) as u8).collect();
        let subject = crate::encoding::pack_ncbi2na_bases(&query);
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let mut word_params = crate::parameters::InitialWordParameters {
            options: crate::options::InitialWordOptions::new_blastn(),
            x_dropoff_max: 8,
            cutoff_score_min: 1,
            cutoffs: vec![crate::parameters::BlastUngappedCutoffs {
                x_dropoff_init: 8,
                cutoff_score: 1,
            }],
            ungapped_extension: false,
            nucl_score_table: crate::parameters::InitialWordParameters::build_nucl_score_table(
                score_params.reward,
                score_params.penalty,
            ),
        };
        word_params.options.word_size = 4;
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[query.len()]);
        let mut init_hitlist = crate::extend::InitHitList::new();
        let mut stats = crate::diagnostics::UngappedStats::default();

        assert_eq!(
            mb_indexed_word_finder(
                &subject,
                query.len() as i32,
                5,
                0,
                &query,
                &query_info,
                &word_params,
                &score_params,
                &mut init_hitlist,
                |oid| oid == 5,
                |_oid, _chunk, list| {
                    let _ = crate::extend::blast_save_initial_hit(list, 3, 3, None);
                    let _ = crate::extend::blast_save_initial_hit(list, 0, 0, None);
                    4
                },
                Some(&mut stats),
            ),
            0
        );

        assert_eq!(stats.lookup_hits, 0);
        assert_eq!(stats.init_extends, 0);
        assert_eq!(init_hitlist.total(), 2);
        assert_eq!(init_hitlist.hits[0].query_offset, 3);
        assert_eq!(init_hitlist.hits[1].query_offset, 0);
        assert!(init_hitlist
            .hits
            .iter()
            .all(|hsp| hsp.ungapped_data.is_none()));
    }

    #[test]
    fn mb_indexed_word_finder_falls_back_for_nonindexed_oid() {
        let query: Vec<u8> = (0..16).map(|i| (i % 4) as u8).collect();
        let subject = crate::encoding::pack_ncbi2na_bases(&query);
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let word_params = crate::parameters::InitialWordParameters {
            options: crate::options::InitialWordOptions::new_blastn(),
            x_dropoff_max: 8,
            cutoff_score_min: 1,
            cutoffs: Vec::new(),
            ungapped_extension: true,
            nucl_score_table: crate::parameters::InitialWordParameters::build_nucl_score_table(
                score_params.reward,
                score_params.penalty,
            ),
        };
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[query.len()]);
        let mut init_hitlist = crate::extend::InitHitList::new();
        let _ = crate::extend::blast_save_initial_hit(&mut init_hitlist, 2, 2, None);
        let mut fallback_called = false;
        let mut get_results_called = false;

        let status = mb_indexed_word_finder_with_fallback(
            &subject,
            query.len() as i32,
            5,
            0,
            &query,
            &query_info,
            &word_params,
            &score_params,
            &mut init_hitlist,
            |_oid| false,
            |_oid, _chunk, _list| {
                get_results_called = true;
                0
            },
            || {
                fallback_called = true;
                7
            },
            None,
        );

        assert_eq!(status, 7);
        assert!(fallback_called);
        assert!(!get_results_called);
        assert_eq!(init_hitlist.total(), 1);
        assert_eq!(init_hitlist.hits[0].query_offset, 2);
    }

    #[test]
    fn mb_indexed_word_finder_passes_oid_chunk_and_skips_zero_word_size() {
        let query: Vec<u8> = (0..16).map(|i| (i % 4) as u8).collect();
        let subject = crate::encoding::pack_ncbi2na_bases(&query);
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let word_params = crate::parameters::InitialWordParameters {
            options: crate::options::InitialWordOptions::new_blastn(),
            x_dropoff_max: 8,
            cutoff_score_min: 1,
            cutoffs: Vec::new(),
            ungapped_extension: true,
            nucl_score_table: crate::parameters::InitialWordParameters::build_nucl_score_table(
                score_params.reward,
                score_params.penalty,
            ),
        };
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[query.len()]);
        let mut init_hitlist = crate::extend::InitHitList::new();
        let _ = crate::extend::blast_save_initial_hit(&mut init_hitlist, 3, 3, None);
        let mut stats = crate::diagnostics::UngappedStats::default();
        let mut callback_seen = None;

        assert_eq!(
            mb_indexed_word_finder(
                &subject,
                query.len() as i32,
                12,
                34,
                &query,
                &query_info,
                &word_params,
                &score_params,
                &mut init_hitlist,
                |oid| oid == 12,
                |oid, chunk, list| {
                    callback_seen = Some((oid, chunk));
                    let _ = crate::extend::blast_save_initial_hit(list, 1, 1, None);
                    0
                },
                Some(&mut stats),
            ),
            0
        );

        assert_eq!(callback_seen, Some((12, 34)));
        assert_eq!(init_hitlist.total(), 1);
        assert_eq!(init_hitlist.hits[0].query_offset, 1);
        assert_eq!(init_hitlist.hits[0].subject_offset, 1);
        assert!(init_hitlist.hits[0].ungapped_data.is_none());
        assert_eq!(stats.lookup_hits, 1);
        assert_eq!(stats.num_seqs_lookup_hits, 1);
        assert_eq!(stats.init_extends, 0);
        assert_eq!(stats.good_init_extends, 0);
        assert_eq!(stats.num_seqs_passed, 0);
    }

    #[test]
    fn mb_indexed_word_finder_counts_extension_failure_from_short_query() {
        let query: Vec<u8> = Vec::new();
        let subject_bases: Vec<u8> = (0..16).map(|i| (i % 4) as u8).collect();
        let subject = crate::encoding::pack_ncbi2na_bases(&subject_bases);
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let word_params = crate::parameters::InitialWordParameters {
            options: crate::options::InitialWordOptions::new_blastn(),
            x_dropoff_max: 8,
            cutoff_score_min: 1,
            cutoffs: Vec::new(),
            ungapped_extension: true,
            nucl_score_table: crate::parameters::InitialWordParameters::build_nucl_score_table(
                score_params.reward,
                score_params.penalty,
            ),
        };
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[0]);
        let mut init_hitlist = crate::extend::InitHitList::new();
        let mut stats = crate::diagnostics::UngappedStats::default();

        assert_eq!(
            mb_indexed_word_finder(
                &subject,
                subject_bases.len() as i32,
                12,
                34,
                &query,
                &query_info,
                &word_params,
                &score_params,
                &mut init_hitlist,
                |oid| oid == 12,
                |_oid, _chunk, list| {
                    let _ = crate::extend::blast_save_initial_hit(list, 0, 0, None);
                    4
                },
                Some(&mut stats),
            ),
            0
        );

        assert_eq!(init_hitlist.total(), 0);
        assert_eq!(stats.lookup_hits, 1);
        assert_eq!(stats.init_extends, 1);
        assert_eq!(stats.good_init_extends, 0);
        assert_eq!(stats.num_seqs_passed, 0);
    }

    #[test]
    fn mb_indexed_word_finder_counts_invalid_index_offsets_as_extension_failures() {
        let query: Vec<u8> = (0..8).map(|i| (i % 4) as u8).collect();
        let subject_bases: Vec<u8> = (0..8).map(|i| (i % 4) as u8).collect();
        let subject = crate::encoding::pack_ncbi2na_bases(&subject_bases);
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let word_params = crate::parameters::InitialWordParameters {
            options: crate::options::InitialWordOptions::new_blastn(),
            x_dropoff_max: 8,
            cutoff_score_min: 1,
            cutoffs: Vec::new(),
            ungapped_extension: true,
            nucl_score_table: crate::parameters::InitialWordParameters::build_nucl_score_table(
                score_params.reward,
                score_params.penalty,
            ),
        };
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[query.len()]);
        let mut init_hitlist = crate::extend::InitHitList::new();
        let mut stats = crate::diagnostics::UngappedStats::default();

        assert_eq!(
            mb_indexed_word_finder(
                &subject,
                subject_bases.len() as i32,
                12,
                34,
                &query,
                &query_info,
                &word_params,
                &score_params,
                &mut init_hitlist,
                |oid| oid == 12,
                |_oid, _chunk, list| {
                    let _ = crate::extend::blast_save_initial_hit(list, 99, 0, None);
                    let _ = crate::extend::blast_save_initial_hit(list, 0, 99, None);
                    4
                },
                Some(&mut stats),
            ),
            0
        );

        assert_eq!(init_hitlist.total(), 0);
        assert_eq!(stats.lookup_hits, 2);
        assert_eq!(stats.init_extends, 2);
        assert_eq!(stats.good_init_extends, 0);
        assert_eq!(stats.num_seqs_passed, 0);
    }

    #[test]
    fn mb_indexed_word_finder_rejects_malformed_packed_subject_before_extension() {
        let query: Vec<u8> = (0..8).map(|i| (i % 4) as u8).collect();
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let word_params = crate::parameters::InitialWordParameters {
            options: crate::options::InitialWordOptions::new_blastn(),
            x_dropoff_max: 8,
            cutoff_score_min: 1,
            cutoffs: Vec::new(),
            ungapped_extension: true,
            nucl_score_table: crate::parameters::InitialWordParameters::build_nucl_score_table(
                score_params.reward,
                score_params.penalty,
            ),
        };
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[query.len()]);
        let mut init_hitlist = crate::extend::InitHitList::new();
        let mut stats = crate::diagnostics::UngappedStats::default();

        assert_eq!(
            mb_indexed_word_finder(
                &[],
                query.len() as i32,
                12,
                34,
                &query,
                &query_info,
                &word_params,
                &score_params,
                &mut init_hitlist,
                |oid| oid == 12,
                |_oid, _chunk, list| {
                    let _ = crate::extend::blast_save_initial_hit(list, 0, 0, None);
                    4
                },
                Some(&mut stats),
            ),
            0
        );

        assert_eq!(init_hitlist.total(), 0);
        assert_eq!(stats.lookup_hits, 1);
        assert_eq!(stats.init_extends, 0);
        assert_eq!(stats.good_init_extends, 0);
        assert_eq!(stats.num_seqs_passed, 0);
    }

    #[test]
    fn mb_indexed_word_finder_handles_huge_callback_word_size_without_overflow() {
        let query: Vec<u8> = (0..8).map(|i| (i % 4) as u8).collect();
        let subject = crate::encoding::pack_ncbi2na_bases(&query);
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let word_params = crate::parameters::InitialWordParameters {
            options: crate::options::InitialWordOptions::new_blastn(),
            x_dropoff_max: 8,
            cutoff_score_min: 1,
            cutoffs: Vec::new(),
            ungapped_extension: true,
            nucl_score_table: crate::parameters::InitialWordParameters::build_nucl_score_table(
                score_params.reward,
                score_params.penalty,
            ),
        };
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[query.len()]);
        let mut init_hitlist = crate::extend::InitHitList::new();
        let mut stats = crate::diagnostics::UngappedStats::default();

        assert_eq!(
            mb_indexed_word_finder(
                &subject,
                query.len() as i32,
                12,
                34,
                &query,
                &query_info,
                &word_params,
                &score_params,
                &mut init_hitlist,
                |oid| oid == 12,
                |_oid, _chunk, list| {
                    let _ = crate::extend::blast_save_initial_hit(list, 99, 0, None);
                    u32::MAX
                },
                Some(&mut stats),
            ),
            0
        );

        assert_eq!(init_hitlist.total(), 0);
        assert_eq!(stats.lookup_hits, 1);
        assert_eq!(stats.init_extends, 1);
        assert_eq!(stats.good_init_extends, 0);
    }

    #[test]
    fn short_read_indexed_word_finder_runs_jumper_traceback_from_index_hits() {
        let query: Vec<u8> = (0..24).map(|i| (i % 4) as u8).collect();
        let subject = crate::encoding::pack_ncbi2na_bases(&query);
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let hit_params = crate::parameters::HitSavingParameters {
            options: crate::options::HitSavingOptions {
                cutoff_score: 1,
                ..Default::default()
            },
            cutoff_score_min: 0,
            low_score: Vec::new(),
            cutoffs: Vec::new(),
            link_hsp_params: None,
            prelim_evalue: 10.0,
        };
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[query.len()]);
        let mut init_hitlist = crate::extend::InitHitList::new();
        let mut hsp_list = crate::hspstream::HspList::new(8);
        let mut jumper = JumperGapAlign {
            left_prelim_block: Some(JumperPrelimEditBlock::default()),
            right_prelim_block: Some(JumperPrelimEditBlock::default()),
            table: Vec::new(),
        };
        let mut stats = crate::diagnostics::GappedStats::default();

        assert_eq!(
            short_read_indexed_word_finder(
                &subject,
                query.len() as i32,
                0,
                8,
                0,
                &query,
                &query_info,
                &score_params,
                &hit_params,
                JumperAlignParams {
                    max_mismatches: 5,
                    mismatch_window: 10,
                    gap_x_dropoff: 30,
                },
                &mut init_hitlist,
                &mut hsp_list,
                &mut jumper,
                |oid| oid == 8,
                |_oid, _chunk, list| {
                    let _ = crate::extend::blast_save_initial_hit(list, 0, 0, None);
                    4
                },
                Some(&mut stats),
            ),
            0
        );

        assert_eq!(stats.extensions, 1);
        assert!(!hsp_list.hsps.is_empty());
        assert_eq!(hsp_list.hsps[0].query_offset, 0);
        assert_eq!(hsp_list.hsps[0].subject_offset, 0);
    }

    #[test]
    fn short_read_indexed_word_finder_splice_path_does_not_save_subject_overhangs() {
        let query: Vec<u8> = (0..24).map(|i| (i % 4) as u8).collect();
        let subject = crate::encoding::pack_ncbi2na_bases(&query);
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let hit_params = crate::parameters::HitSavingParameters {
            options: crate::options::HitSavingOptions {
                cutoff_score: 1,
                splice: true,
                ..Default::default()
            },
            cutoff_score_min: 0,
            low_score: Vec::new(),
            cutoffs: Vec::new(),
            link_hsp_params: None,
            prelim_evalue: 10.0,
        };
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[query.len()]);
        let mut init_hitlist = crate::extend::InitHitList::new();
        let mut hsp_list = crate::hspstream::HspList::new(8);
        let mut jumper = JumperGapAlign {
            left_prelim_block: Some(JumperPrelimEditBlock::default()),
            right_prelim_block: Some(JumperPrelimEditBlock::default()),
            table: Vec::new(),
        };

        assert_eq!(
            short_read_indexed_word_finder(
                &subject,
                query.len() as i32,
                0,
                8,
                0,
                &query,
                &query_info,
                &score_params,
                &hit_params,
                JumperAlignParams {
                    max_mismatches: 5,
                    mismatch_window: 10,
                    gap_x_dropoff: 30,
                },
                &mut init_hitlist,
                &mut hsp_list,
                &mut jumper,
                |oid| oid == 8,
                |_oid, _chunk, list| {
                    let _ = crate::extend::blast_save_initial_hit(list, 0, 0, None);
                    4
                },
                None,
            ),
            0
        );

        let map_info = hsp_list.hsps[0]
            .map_info
            .as_ref()
            .expect("indexed short-read HSP map_info");
        assert!(map_info.subject_overhangs.is_none());
    }

    #[test]
    fn short_read_indexed_word_finder_suppresses_redundant_diagonal_hits() {
        let query: Vec<u8> = (0..32).map(|i| (i % 4) as u8).collect();
        let subject = crate::encoding::pack_ncbi2na_bases(&query);
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let hit_params = crate::parameters::HitSavingParameters {
            options: crate::options::HitSavingOptions {
                cutoff_score: 1,
                ..Default::default()
            },
            cutoff_score_min: 0,
            low_score: Vec::new(),
            cutoffs: Vec::new(),
            link_hsp_params: None,
            prelim_evalue: 10.0,
        };
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[query.len()]);
        let mut init_hitlist = crate::extend::InitHitList::new();
        let mut hsp_list = crate::hspstream::HspList::new(8);
        let mut jumper = JumperGapAlign {
            left_prelim_block: Some(JumperPrelimEditBlock::default()),
            right_prelim_block: Some(JumperPrelimEditBlock::default()),
            table: Vec::new(),
        };
        let mut stats = crate::diagnostics::GappedStats::default();

        assert_eq!(
            short_read_indexed_word_finder(
                &subject,
                query.len() as i32,
                0,
                8,
                0,
                &query,
                &query_info,
                &score_params,
                &hit_params,
                JumperAlignParams {
                    max_mismatches: 5,
                    mismatch_window: 10,
                    gap_x_dropoff: 30,
                },
                &mut init_hitlist,
                &mut hsp_list,
                &mut jumper,
                |oid| oid == 8,
                |_oid, _chunk, list| {
                    let _ = crate::extend::blast_save_initial_hit(list, 0, 0, None);
                    let _ = crate::extend::blast_save_initial_hit(list, 1, 1, None);
                    4
                },
                Some(&mut stats),
            ),
            0
        );

        assert_eq!(stats.extensions, 1);
        assert_eq!(stats.good_extensions, 1);
        assert_eq!(hsp_list.hsps.len(), 1);
        assert_eq!(hsp_list.hsps[0].query_offset, 0);
        assert_eq!(hsp_list.hsps[0].subject_offset, 0);
    }

    #[test]
    fn short_read_indexed_word_finder_counts_extension_rejected_by_good_align() {
        let query: Vec<u8> = (0..24).map(|i| (i % 4) as u8).collect();
        let subject = crate::encoding::pack_ncbi2na_bases(&query);
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let hit_params = crate::parameters::HitSavingParameters {
            options: crate::options::HitSavingOptions {
                cutoff_score: 1000,
                ..Default::default()
            },
            cutoff_score_min: 1000,
            low_score: Vec::new(),
            cutoffs: Vec::new(),
            link_hsp_params: None,
            prelim_evalue: 10.0,
        };
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[query.len()]);
        let mut init_hitlist = crate::extend::InitHitList::new();
        let mut hsp_list = crate::hspstream::HspList::new(8);
        let mut jumper = JumperGapAlign {
            left_prelim_block: Some(JumperPrelimEditBlock::default()),
            right_prelim_block: Some(JumperPrelimEditBlock::default()),
            table: Vec::new(),
        };
        let mut stats = crate::diagnostics::GappedStats::default();

        assert_eq!(
            short_read_indexed_word_finder(
                &subject,
                query.len() as i32,
                0,
                8,
                0,
                &query,
                &query_info,
                &score_params,
                &hit_params,
                JumperAlignParams {
                    max_mismatches: 5,
                    mismatch_window: 10,
                    gap_x_dropoff: 30,
                },
                &mut init_hitlist,
                &mut hsp_list,
                &mut jumper,
                |oid| oid == 8,
                |_oid, _chunk, list| {
                    let _ = crate::extend::blast_save_initial_hit(list, 0, 0, None);
                    4
                },
                Some(&mut stats),
            ),
            0
        );

        assert_eq!(stats.extensions, 1);
        assert_eq!(stats.good_extensions, 0);
        assert!(hsp_list.hsps.is_empty());
    }

    #[test]
    fn short_read_indexed_word_finder_counts_invalid_offsets_as_failed_extensions() {
        let query: Vec<u8> = (0..16).map(|i| (i % 4) as u8).collect();
        let subject = crate::encoding::pack_ncbi2na_bases(&query);
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let hit_params = crate::parameters::HitSavingParameters {
            options: crate::options::HitSavingOptions {
                cutoff_score: 1,
                ..Default::default()
            },
            cutoff_score_min: 0,
            low_score: Vec::new(),
            cutoffs: Vec::new(),
            link_hsp_params: None,
            prelim_evalue: 10.0,
        };
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[query.len()]);
        let mut init_hitlist = crate::extend::InitHitList::new();
        let mut hsp_list = crate::hspstream::HspList::new(8);
        let mut jumper = JumperGapAlign {
            left_prelim_block: Some(JumperPrelimEditBlock::default()),
            right_prelim_block: Some(JumperPrelimEditBlock::default()),
            table: Vec::new(),
        };
        let mut stats = crate::diagnostics::GappedStats::default();

        assert_eq!(
            short_read_indexed_word_finder(
                &subject,
                query.len() as i32,
                0,
                8,
                0,
                &query,
                &query_info,
                &score_params,
                &hit_params,
                JumperAlignParams {
                    max_mismatches: 5,
                    mismatch_window: 10,
                    gap_x_dropoff: 30,
                },
                &mut init_hitlist,
                &mut hsp_list,
                &mut jumper,
                |oid| oid == 8,
                |_oid, _chunk, list| {
                    let _ = crate::extend::blast_save_initial_hit(list, 99, 0, None);
                    let _ = crate::extend::blast_save_initial_hit(list, 0, 99, None);
                    let _ = crate::extend::blast_save_initial_hit(list, i32::MIN, 0, None);
                    4
                },
                Some(&mut stats),
            ),
            0
        );

        assert_eq!(init_hitlist.total(), 3);
        assert!(hsp_list.hsps.is_empty());
        assert_eq!(stats.extensions, 3);
        assert_eq!(stats.good_extensions, 0);
    }

    #[test]
    fn short_read_indexed_word_finder_handles_negative_hit_before_shifted_context() {
        let query: Vec<u8> = (0..16).map(|i| (i % 4) as u8).collect();
        let subject = crate::encoding::pack_ncbi2na_bases(&query);
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let hit_params = crate::parameters::HitSavingParameters {
            options: crate::options::HitSavingOptions {
                cutoff_score: 1,
                ..Default::default()
            },
            cutoff_score_min: 0,
            low_score: Vec::new(),
            cutoffs: Vec::new(),
            link_hsp_params: None,
            prelim_evalue: 10.0,
        };
        let query_info = crate::queryinfo::QueryInfo {
            num_queries: 1,
            contexts: vec![crate::queryinfo::ContextInfo {
                query_offset: 1,
                query_length: 15,
                eff_searchsp: 0,
                length_adjustment: 0,
                query_index: 0,
                frame: 0,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: 15,
            min_length: 0,
        };
        let mut init_hitlist = crate::extend::InitHitList::new();
        let mut hsp_list = crate::hspstream::HspList::new(8);
        let mut jumper = JumperGapAlign {
            left_prelim_block: Some(JumperPrelimEditBlock::default()),
            right_prelim_block: Some(JumperPrelimEditBlock::default()),
            table: Vec::new(),
        };
        let mut stats = crate::diagnostics::GappedStats::default();

        assert_eq!(
            short_read_indexed_word_finder(
                &subject,
                query.len() as i32,
                0,
                8,
                0,
                &query,
                &query_info,
                &score_params,
                &hit_params,
                JumperAlignParams {
                    max_mismatches: 5,
                    mismatch_window: 10,
                    gap_x_dropoff: 30,
                },
                &mut init_hitlist,
                &mut hsp_list,
                &mut jumper,
                |oid| oid == 8,
                |_oid, _chunk, list| {
                    let _ = crate::extend::blast_save_initial_hit(list, i32::MIN, 0, None);
                    4
                },
                Some(&mut stats),
            ),
            0
        );

        assert_eq!(init_hitlist.total(), 1);
        assert!(hsp_list.hsps.is_empty());
        assert_eq!(stats.extensions, 1);
        assert_eq!(stats.good_extensions, 0);
    }

    #[test]
    fn short_read_indexed_word_finder_skips_malformed_packed_subject() {
        let query: Vec<u8> = (0..16).map(|i| (i % 4) as u8).collect();
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let hit_params = crate::parameters::HitSavingParameters {
            options: crate::options::HitSavingOptions {
                cutoff_score: 1,
                ..Default::default()
            },
            cutoff_score_min: 0,
            low_score: Vec::new(),
            cutoffs: Vec::new(),
            link_hsp_params: None,
            prelim_evalue: 10.0,
        };
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[query.len()]);
        let mut init_hitlist = crate::extend::InitHitList::new();
        let mut hsp_list = crate::hspstream::HspList::new(8);
        let mut jumper = JumperGapAlign {
            left_prelim_block: Some(JumperPrelimEditBlock::default()),
            right_prelim_block: Some(JumperPrelimEditBlock::default()),
            table: Vec::new(),
        };
        let mut stats = crate::diagnostics::GappedStats::default();

        assert_eq!(
            short_read_indexed_word_finder(
                &[],
                query.len() as i32,
                0,
                8,
                0,
                &query,
                &query_info,
                &score_params,
                &hit_params,
                JumperAlignParams {
                    max_mismatches: 5,
                    mismatch_window: 10,
                    gap_x_dropoff: 30,
                },
                &mut init_hitlist,
                &mut hsp_list,
                &mut jumper,
                |oid| oid == 8,
                |_oid, _chunk, list| {
                    let _ = crate::extend::blast_save_initial_hit(list, 0, 0, None);
                    4
                },
                Some(&mut stats),
            ),
            0
        );

        assert_eq!(init_hitlist.total(), 1);
        assert!(hsp_list.hsps.is_empty());
        assert_eq!(stats.extensions, 0);
        assert_eq!(stats.good_extensions, 0);
        assert_eq!(stats.num_seqs_passed, 0);
    }

    #[test]
    fn short_read_indexed_word_finder_handles_huge_callback_word_size_without_overflow() {
        let query: Vec<u8> = (0..16).map(|i| (i % 4) as u8).collect();
        let subject = crate::encoding::pack_ncbi2na_bases(&query);
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let hit_params = crate::parameters::HitSavingParameters {
            options: crate::options::HitSavingOptions {
                cutoff_score: 1,
                ..Default::default()
            },
            cutoff_score_min: 0,
            low_score: Vec::new(),
            cutoffs: Vec::new(),
            link_hsp_params: None,
            prelim_evalue: 10.0,
        };
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[query.len()]);
        let mut init_hitlist = crate::extend::InitHitList::new();
        let mut hsp_list = crate::hspstream::HspList::new(8);
        let mut jumper = JumperGapAlign {
            left_prelim_block: Some(JumperPrelimEditBlock::default()),
            right_prelim_block: Some(JumperPrelimEditBlock::default()),
            table: Vec::new(),
        };
        let mut stats = crate::diagnostics::GappedStats::default();

        assert_eq!(
            short_read_indexed_word_finder(
                &subject,
                query.len() as i32,
                0,
                8,
                0,
                &query,
                &query_info,
                &score_params,
                &hit_params,
                JumperAlignParams {
                    max_mismatches: 5,
                    mismatch_window: 10,
                    gap_x_dropoff: 30,
                },
                &mut init_hitlist,
                &mut hsp_list,
                &mut jumper,
                |oid| oid == 8,
                |_oid, _chunk, list| {
                    let _ = crate::extend::blast_save_initial_hit(list, 0, 0, None);
                    u32::MAX
                },
                Some(&mut stats),
            ),
            0
        );

        assert_eq!(init_hitlist.total(), 1);
        assert_eq!(stats.extensions, 1);
        assert!(!hsp_list.hsps.is_empty());
    }

    #[test]
    fn short_read_indexed_word_finder_zero_index_hits_leaves_stats_clear() {
        let query: Vec<u8> = (0..24).map(|i| (i % 4) as u8).collect();
        let subject = crate::encoding::pack_ncbi2na_bases(&query);
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let hit_params = crate::parameters::HitSavingParameters {
            options: crate::options::HitSavingOptions {
                cutoff_score: 1,
                ..Default::default()
            },
            cutoff_score_min: 0,
            low_score: Vec::new(),
            cutoffs: Vec::new(),
            link_hsp_params: None,
            prelim_evalue: 10.0,
        };
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[query.len()]);
        let mut init_hitlist = crate::extend::InitHitList::new();
        let mut hsp_list = crate::hspstream::HspList::new(8);
        let mut jumper = JumperGapAlign {
            left_prelim_block: Some(JumperPrelimEditBlock::default()),
            right_prelim_block: Some(JumperPrelimEditBlock::default()),
            table: Vec::new(),
        };
        let mut stats = crate::diagnostics::GappedStats::default();

        assert_eq!(
            short_read_indexed_word_finder(
                &subject,
                query.len() as i32,
                0,
                8,
                0,
                &query,
                &query_info,
                &score_params,
                &hit_params,
                JumperAlignParams {
                    max_mismatches: 5,
                    mismatch_window: 10,
                    gap_x_dropoff: 30,
                },
                &mut init_hitlist,
                &mut hsp_list,
                &mut jumper,
                |oid| oid == 8,
                |_oid, _chunk, _list| 4,
                Some(&mut stats),
            ),
            0
        );

        assert_eq!(init_hitlist.total(), 0);
        assert!(hsp_list.hsps.is_empty());
        assert_eq!(stats.extensions, 0);
        assert_eq!(stats.good_extensions, 0);
        assert_eq!(stats.num_seqs_passed, 0);
    }

    #[test]
    fn short_read_indexed_word_finder_skips_hits_without_query_context() {
        let query: Vec<u8> = (0..16).map(|i| (i % 4) as u8).collect();
        let subject = crate::encoding::pack_ncbi2na_bases(&query);
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let hit_params = crate::parameters::HitSavingParameters {
            options: crate::options::HitSavingOptions {
                cutoff_score: 1,
                ..Default::default()
            },
            cutoff_score_min: 0,
            low_score: Vec::new(),
            cutoffs: Vec::new(),
            link_hsp_params: None,
            prelim_evalue: 10.0,
        };
        let query_info = crate::queryinfo::QueryInfo {
            num_queries: 0,
            contexts: Vec::new(),
            max_length: 0,
            min_length: 0,
        };
        let mut init_hitlist = crate::extend::InitHitList::new();
        let mut hsp_list = crate::hspstream::HspList::new(8);
        let mut jumper = JumperGapAlign {
            left_prelim_block: Some(JumperPrelimEditBlock::default()),
            right_prelim_block: Some(JumperPrelimEditBlock::default()),
            table: Vec::new(),
        };
        let mut stats = crate::diagnostics::GappedStats::default();

        assert_eq!(
            short_read_indexed_word_finder(
                &subject,
                query.len() as i32,
                0,
                8,
                0,
                &query,
                &query_info,
                &score_params,
                &hit_params,
                JumperAlignParams {
                    max_mismatches: 5,
                    mismatch_window: 10,
                    gap_x_dropoff: 30,
                },
                &mut init_hitlist,
                &mut hsp_list,
                &mut jumper,
                |oid| oid == 8,
                |_oid, _chunk, list| {
                    let _ = crate::extend::blast_save_initial_hit(list, 0, 0, None);
                    4
                },
                Some(&mut stats),
            ),
            0
        );

        assert_eq!(init_hitlist.total(), 1);
        assert!(hsp_list.hsps.is_empty());
        assert_eq!(stats.extensions, 0);
        assert_eq!(stats.good_extensions, 0);
        assert_eq!(stats.num_seqs_passed, 0);
    }

    #[test]
    fn short_read_indexed_word_finder_skips_truncated_query_context() {
        let query: Vec<u8> = (0..8).map(|i| (i % 4) as u8).collect();
        let subject_bases: Vec<u8> = (0..16).map(|i| (i % 4) as u8).collect();
        let subject = crate::encoding::pack_ncbi2na_bases(&subject_bases);
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let hit_params = crate::parameters::HitSavingParameters {
            options: crate::options::HitSavingOptions {
                cutoff_score: 1,
                ..Default::default()
            },
            cutoff_score_min: 0,
            low_score: Vec::new(),
            cutoffs: Vec::new(),
            link_hsp_params: None,
            prelim_evalue: 10.0,
        };
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[16]);
        let mut init_hitlist = crate::extend::InitHitList::new();
        let mut hsp_list = crate::hspstream::HspList::new(8);
        let mut jumper = JumperGapAlign {
            left_prelim_block: Some(JumperPrelimEditBlock::default()),
            right_prelim_block: Some(JumperPrelimEditBlock::default()),
            table: Vec::new(),
        };
        let mut stats = crate::diagnostics::GappedStats::default();

        assert_eq!(
            short_read_indexed_word_finder(
                &subject,
                subject_bases.len() as i32,
                0,
                8,
                0,
                &query,
                &query_info,
                &score_params,
                &hit_params,
                JumperAlignParams {
                    max_mismatches: 5,
                    mismatch_window: 10,
                    gap_x_dropoff: 30,
                },
                &mut init_hitlist,
                &mut hsp_list,
                &mut jumper,
                |oid| oid == 8,
                |_oid, _chunk, list| {
                    let _ = crate::extend::blast_save_initial_hit(list, 0, 0, None);
                    4
                },
                Some(&mut stats),
            ),
            0
        );

        assert_eq!(init_hitlist.total(), 1);
        assert!(hsp_list.hsps.is_empty());
        assert_eq!(stats.extensions, 0);
        assert_eq!(stats.good_extensions, 0);
        assert_eq!(stats.num_seqs_passed, 0);
    }

    #[test]
    fn short_read_indexed_word_finder_falls_back_for_nonindexed_oid() {
        let query: Vec<u8> = (0..16).map(|i| (i % 4) as u8).collect();
        let subject = crate::encoding::pack_ncbi2na_bases(&query);
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let hit_params = crate::parameters::HitSavingParameters {
            options: crate::options::HitSavingOptions {
                cutoff_score: 1,
                ..Default::default()
            },
            cutoff_score_min: 0,
            low_score: Vec::new(),
            cutoffs: Vec::new(),
            link_hsp_params: None,
            prelim_evalue: 10.0,
        };
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[query.len()]);
        let mut init_hitlist = crate::extend::InitHitList::new();
        let mut hsp_list = crate::hspstream::HspList::new(8);
        let _ = crate::extend::blast_save_initial_hit(&mut init_hitlist, 2, 2, None);
        hsp_list.add_hsp(mapper_saved_hsp(0, 11, 0, 11));
        let mut jumper = JumperGapAlign {
            left_prelim_block: Some(JumperPrelimEditBlock::default()),
            right_prelim_block: Some(JumperPrelimEditBlock::default()),
            table: Vec::new(),
        };
        let mut fallback_called = false;
        let mut get_results_called = false;

        let status = short_read_indexed_word_finder_with_fallback(
            &subject,
            query.len() as i32,
            0,
            8,
            0,
            &query,
            &query_info,
            &score_params,
            &hit_params,
            JumperAlignParams {
                max_mismatches: 5,
                mismatch_window: 10,
                gap_x_dropoff: 30,
            },
            &mut init_hitlist,
            &mut hsp_list,
            &mut jumper,
            |_oid| false,
            |_oid, _chunk, _list| {
                get_results_called = true;
                0
            },
            || {
                fallback_called = true;
                9
            },
            None,
        );

        assert_eq!(status, 9);
        assert!(fallback_called);
        assert!(!get_results_called);
        assert_eq!(init_hitlist.total(), 1);
        assert_eq!(init_hitlist.hits[0].query_offset, 2);
        assert_eq!(hsp_list.hsps.len(), 1);
        assert_eq!(hsp_list.hsps[0].score, 11);
    }

    #[test]
    fn short_read_indexed_word_finder_passes_oid_chunk_and_skips_zero_word_size() {
        let query: Vec<u8> = (0..16).map(|i| (i % 4) as u8).collect();
        let subject = crate::encoding::pack_ncbi2na_bases(&query);
        let scoring_options = crate::options::ScoringOptions::new_blastn();
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let hit_params = crate::parameters::HitSavingParameters {
            options: crate::options::HitSavingOptions {
                cutoff_score: 1,
                ..Default::default()
            },
            cutoff_score_min: 0,
            low_score: Vec::new(),
            cutoffs: Vec::new(),
            link_hsp_params: None,
            prelim_evalue: 10.0,
        };
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[query.len()]);
        let mut init_hitlist = crate::extend::InitHitList::new();
        let _ = crate::extend::blast_save_initial_hit(&mut init_hitlist, 3, 3, None);
        let mut hsp_list = crate::hspstream::HspList::new(8);
        let mut jumper = JumperGapAlign {
            left_prelim_block: Some(JumperPrelimEditBlock::default()),
            right_prelim_block: Some(JumperPrelimEditBlock::default()),
            table: Vec::new(),
        };
        let mut stats = crate::diagnostics::GappedStats::default();
        let mut callback_seen = None;

        assert_eq!(
            short_read_indexed_word_finder(
                &subject,
                query.len() as i32,
                0,
                12,
                34,
                &query,
                &query_info,
                &score_params,
                &hit_params,
                JumperAlignParams {
                    max_mismatches: 5,
                    mismatch_window: 10,
                    gap_x_dropoff: 30,
                },
                &mut init_hitlist,
                &mut hsp_list,
                &mut jumper,
                |oid| oid == 12,
                |oid, chunk, list| {
                    callback_seen = Some((oid, chunk));
                    let _ = crate::extend::blast_save_initial_hit(list, 1, 1, None);
                    0
                },
                Some(&mut stats),
            ),
            0
        );

        assert_eq!(callback_seen, Some((12, 34)));
        assert_eq!(init_hitlist.total(), 1);
        assert_eq!(init_hitlist.hits[0].query_offset, 1);
        assert_eq!(init_hitlist.hits[0].subject_offset, 1);
        assert!(hsp_list.hsps.is_empty());
        assert_eq!(stats.extensions, 0);
        assert_eq!(stats.good_extensions, 0);
    }

    #[test]
    fn gap_edit_script_combine_moves_append_and_merges_same_ops() {
        let mut base = Some(GapEditScript {
            ops: vec![(GapAlignOpType::Sub, 3)],
        });
        let mut append = Some(GapEditScript {
            ops: vec![(GapAlignOpType::Sub, 2), (GapAlignOpType::Del, 1)],
        });
        let combined = gap_edit_script_combine(&mut base, &mut append).expect("combined");
        assert_eq!(
            combined.ops,
            vec![(GapAlignOpType::Sub, 5), (GapAlignOpType::Del, 1)]
        );
        assert!(append.is_none());

        let mut empty_base = None;
        let mut append_only = Some(GapEditScript {
            ops: vec![(GapAlignOpType::Ins, 4)],
        });
        assert_eq!(
            gap_edit_script_combine(&mut empty_base, &mut append_only)
                .unwrap()
                .ops,
            vec![(GapAlignOpType::Ins, 4)]
        );
    }

    /// Helper: BLASTNA byte to IUPAC char (A=0,C=1,G=2,T=3).
    fn to_iupac(b: u8) -> char {
        crate::encoding::blastna_to_iupacna_char(b)
    }

    #[test]
    fn test_render_alignment_perfect_match() {
        let mut esp = GapEditScript::new();
        esp.push(GapAlignOpType::Sub, 5);
        let query = &[0u8, 1, 2, 3, 0]; // ACGTA
        let subject = &[0u8, 1, 2, 3, 0]; // ACGTA
        let (qseq, sseq) = esp.render_alignment(query, subject, to_iupac);
        assert_eq!(qseq, "ACGTA");
        assert_eq!(sseq, "ACGTA");
    }

    #[test]
    fn test_render_alignment_with_mismatch() {
        let mut esp = GapEditScript::new();
        esp.push(GapAlignOpType::Sub, 5);
        let query = &[0u8, 1, 2, 3, 0]; // ACGTA
        let subject = &[0u8, 3, 2, 1, 0]; // ATGCA
        let (qseq, sseq) = esp.render_alignment(query, subject, to_iupac);
        assert_eq!(qseq, "ACGTA");
        assert_eq!(sseq, "ATGCA");
    }

    #[test]
    fn test_render_alignment_with_gap_in_query() {
        // Del = gap in query
        let mut esp = GapEditScript::new();
        esp.push(GapAlignOpType::Sub, 3);
        esp.push(GapAlignOpType::Del, 2);
        esp.push(GapAlignOpType::Sub, 2);
        let query = &[0u8, 1, 2, 3, 0]; // ACGTA (5 bases, 2 gapped)
        let subject = &[0u8, 1, 2, 1, 1, 3, 0]; // ACGCCTA (7 bases)
        let (qseq, sseq) = esp.render_alignment(query, subject, to_iupac);
        assert_eq!(qseq, "ACG--TA");
        assert_eq!(sseq, "ACGCCTA");
    }

    #[test]
    fn test_render_alignment_with_gap_in_subject() {
        // Ins = gap in subject
        let mut esp = GapEditScript::new();
        esp.push(GapAlignOpType::Sub, 2);
        esp.push(GapAlignOpType::Ins, 2);
        esp.push(GapAlignOpType::Sub, 2);
        let query = &[0u8, 1, 2, 3, 0, 1]; // ACGTAC (6 bases)
        let subject = &[0u8, 1, 0, 1]; // ACAC (4 bases)
        let (qseq, sseq) = esp.render_alignment(query, subject, to_iupac);
        assert_eq!(qseq, "ACGTAC");
        assert_eq!(sseq, "AC--AC");
    }
}
