//! NCBI-style greedy nucleotide alignment for megablast.
//!
//! This module ports the non-affine greedy path from
//! `greedy_align.c`/`BLAST_GreedyGappedAlignment` that NCBI uses when
//! `gap_open == 0 && gap_extend == 0`.

use crate::gapinfo::{GapAlignOpType, GapEditScript};

const INVALID_OFFSET: i32 = -2;
const MAX_DBSEQ_LEN: usize = 5_000_000;
const GREEDY_MAX_COST_FRACTION: usize = 2;
const GREEDY_MAX_COST: usize = 1000;
const MBSPACE_MIN_SPACE: usize = 1_000_000;

/// Rust equivalent of NCBI `SGreedyOffset`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct SGreedyOffset {
    pub insert_off: i32,
    pub match_off: i32,
    pub delete_off: i32,
}

#[derive(Clone, Debug)]
struct SMBSpaceChunk {
    space_array: Vec<SGreedyOffset>,
    space_used: usize,
}

impl SMBSpaceChunk {
    fn new(num_space_arrays: usize) -> Self {
        let space_allocated = MBSPACE_MIN_SPACE.max(num_space_arrays);
        SMBSpaceChunk {
            space_array: vec![SGreedyOffset::default(); space_allocated],
            space_used: 0,
        }
    }

    fn space_allocated(&self) -> usize {
        self.space_array.len()
    }
}

/// Rust equivalent of the C `SMBSpace` linked pool.
#[derive(Clone, Debug)]
pub struct SMBSpace {
    chunks: Vec<SMBSpaceChunk>,
}

impl SMBSpace {
    pub fn space_used_per_chunk(&self) -> Vec<usize> {
        self.chunks.iter().map(|chunk| chunk.space_used).collect()
    }

    pub fn space_allocated_per_chunk(&self) -> Vec<usize> {
        self.chunks
            .iter()
            .map(SMBSpaceChunk::space_allocated)
            .collect()
    }
}

/// Port of NCBI `MBSpaceNew` (`greedy_align.c:44`).
/// naming: Rust keeps MBSpace as the established `mbspace` token.
pub fn mbspace_new(num_space_arrays: i32) -> SMBSpace {
    let requested = num_space_arrays.max(0) as usize;
    let space_allocated = MBSPACE_MIN_SPACE.max(requested);
    let space_array = vec![SGreedyOffset::default(); space_allocated];
    let first_chunk = SMBSpaceChunk {
        space_array,
        space_used: 0,
    };

    SMBSpace {
        chunks: vec![first_chunk],
    }
}

/// NCBI: s_RefreshMBSpace (greedy_align.c).
pub fn s_refresh_mb_space(space: Option<&mut SMBSpace>) {
    let Some(space) = space else {
        return;
    };
    for chunk in &mut space.chunks {
        chunk.space_used = 0;
    }
}

/// Rust ownership equivalent of NCBI `MBSpaceFree`.
/// naming: Rust keeps MBSpace as the established `mbspace` token.
pub fn mbspace_free(space: &mut Option<SMBSpace>) {
    *space = None;
}

/// NCBI: s_GetMBSpace (greedy_align.c).
pub fn s_get_mb_space(pool: Option<&mut SMBSpace>, num_alloc: i32) -> Option<&mut [SGreedyOffset]> {
    let pool = pool?;
    if num_alloc < 0 {
        return None;
    }
    let num_alloc = num_alloc as usize;
    let mut chunk_index = 0;

    loop {
        if chunk_index >= pool.chunks.len() {
            return None;
        }

        let chunk = &pool.chunks[chunk_index];
        if chunk.space_used + num_alloc <= chunk.space_allocated() {
            let chunk = &mut pool.chunks[chunk_index];
            let start = chunk.space_used;
            chunk.space_used += num_alloc;
            return Some(&mut chunk.space_array[start..start + num_alloc]);
        }

        if chunk_index + 1 == pool.chunks.len() {
            pool.chunks.push(SMBSpaceChunk::new(num_alloc));
        }
        chunk_index += 1;
    }
}

/// Rust-owned equivalent of NCBI `SGreedyAlignMem`.
#[derive(Clone, Debug)]
pub struct SGreedyAlignMem {
    pub max_dist: i32,
    pub xdrop: i32,
    pub last_seq2_off: Option<Vec<Vec<i32>>>,
    pub last_seq2_off_affine: Option<Vec<Vec<SGreedyOffset>>>,
    pub diag_bounds: Vec<i32>,
    pub max_score: Vec<i32>,
    pub space: SMBSpace,
}

// NCBI: BLAST_Gdb3 (ncbi_math.c:422).
fn blast_gdb3(a: &mut i32, b: &mut i32, c: &mut i32) -> i32 {
    let g = if *b == 0 {
        crate::math::blast_gcd(*a, *c)
    } else {
        crate::math::blast_gcd(*a, crate::math::blast_gcd(*b, *c))
    };
    if g > 1 {
        *a /= g;
        *b /= g;
        *c /= g;
    }
    g
}

/// NCBI: s_BlastGreedyAlignMemAlloc (blast_gapalign.c:172).
///
/// The C function receives score/extension parameter structs. This Rust port
/// takes the scalar fields used by the allocation logic and materializes the
/// same non-affine or affine work buffers in owned vectors.
pub fn s_blast_greedy_align_mem_alloc(
    score_reward: i32,
    score_penalty: i32,
    score_gap_open: i32,
    score_gap_extend: i32,
    max_d: i32,
    xdrop: i32,
) -> Option<SGreedyAlignMem> {
    if xdrop == 0 || max_d < 0 {
        return None;
    }

    let (reward, penalty, gap_open, mut gap_extend) = if score_reward % 2 == 1 {
        (
            2 * score_reward,
            -2 * score_penalty,
            2 * score_gap_open,
            2 * score_gap_extend,
        )
    } else {
        (
            score_reward,
            -score_penalty,
            score_gap_open,
            score_gap_extend,
        )
    };

    if gap_open == 0 && gap_extend == 0 {
        gap_extend = reward / 2 + penalty;
    }

    let max_d_usize = max_d as usize;
    let mut mem = SGreedyAlignMem {
        max_dist: max_d,
        xdrop,
        last_seq2_off: None,
        last_seq2_off_affine: None,
        diag_bounds: Vec::new(),
        max_score: Vec::new(),
        space: mbspace_new(0),
    };

    let mut max_score_max_d = max_d;
    let d_diff = if score_gap_open == 0 && score_gap_extend == 0 {
        let denom = penalty + reward;
        if denom <= 0 {
            return None;
        }
        let d_diff = (xdrop + reward / 2) / denom + 1;
        let width = max_d_usize * 2 + 6;
        mem.last_seq2_off = Some(vec![vec![0; width]; max_d_usize + 2]);
        d_diff
    } else {
        let mut mis_cost = reward + penalty;
        let mut gap_open_cost = gap_open;
        let mut ge_cost = gap_extend + reward / 2;
        let max_d_1 = max_d_usize;
        max_score_max_d = max_d * ge_cost;
        let scaled_max_d = (max_d * ge_cost).max(0) as usize;
        let max_cost = mis_cost.max(gap_open_cost + ge_cost).max(0) as usize;
        let gd = blast_gdb3(&mut mis_cost, &mut gap_open_cost, &mut ge_cost);
        if gd <= 0 {
            return None;
        }
        let d_diff = (xdrop + reward / 2) / gd + 1;
        mem.diag_bounds = vec![0; 2 * (scaled_max_d + 1 + max_cost)];
        let affine_rows = scaled_max_d.max(max_cost) + 2;
        mem.last_seq2_off_affine = Some(vec![
            vec![SGreedyOffset::default(); 2 * max_d_1 + 6];
            affine_rows
        ]);
        d_diff
    };

    let max_score_len = (max_score_max_d + 1 + d_diff).max(0) as usize;
    mem.max_score = vec![0; max_score_len];
    Some(mem)
}

#[derive(Clone, Copy, Debug, Default)]
struct GreedySeed {
    start_q: usize,
    start_s: usize,
    match_length: usize,
}

#[derive(Clone, Debug)]
struct GreedySideAlignment {
    dist: i32,
    q_ext: usize,
    s_ext: usize,
    prelim: PrelimEditBlock,
    seed: GreedySeed,
}

#[derive(Clone, Copy, Debug)]
struct GreedyAlignmentExtents {
    score: i32,
    query_start: usize,
    query_end: usize,
    subject_start: usize,
    subject_end: usize,
}

#[derive(Clone, Copy, Debug)]
struct GreedyScaledScores {
    reward: i32,
    penalty: i32,
    x_dropoff: i32,
}

#[derive(Clone, Copy, Debug)]
struct PrelimEditOp {
    op_type: GapAlignOpType,
    num: i32,
}

#[derive(Clone, Debug, Default)]
struct PrelimEditBlock {
    edit_ops: Vec<PrelimEditOp>,
    last_op: Option<GapAlignOpType>,
}

impl PrelimEditBlock {
    fn add(&mut self, op_type: GapAlignOpType, num_ops: i32) {
        if num_ops == 0 {
            return;
        }
        if self.last_op == Some(op_type) {
            if let Some(last) = self.edit_ops.last_mut() {
                last.num += num_ops;
                return;
            }
        }
        self.last_op = Some(op_type);
        self.edit_ops.push(PrelimEditOp {
            op_type,
            num: num_ops,
        });
    }
}

fn ceil_div_i32(num: i32, den: i32) -> i32 {
    debug_assert!(den > 0);
    if num >= 0 {
        (num + den - 1) / den
    } else {
        num / den
    }
}

fn initial_greedy_max_dist(max_subject_length: usize) -> usize {
    let max_subject_length = max_subject_length.min(MAX_DBSEQ_LEN);
    GREEDY_MAX_COST.min(max_subject_length / GREEDY_MAX_COST_FRACTION + 1)
}

// NCBI: s_FindFirstMismatch (greedy_align.c:315).
// naming: This helper covers the uncompressed path; packed-subject matching is
// kept in a separate native helper.
fn s_find_first_mismatch(
    seq1: &[u8],
    seq2: &[u8],
    len1: usize,
    len2: usize,
    mut seq1_index: usize,
    mut seq2_index: usize,
    reverse: bool,
) -> usize {
    let start = seq1_index;
    if reverse {
        while seq1_index < len1
            && seq2_index < len2
            && seq1[len1 - 1 - seq1_index] < 4
            && seq1[len1 - 1 - seq1_index] == seq2[len2 - 1 - seq2_index]
        {
            seq1_index += 1;
            seq2_index += 1;
        }
    } else {
        while seq1_index < len1
            && seq2_index < len2
            && seq1[seq1_index] < 4
            && seq1[seq1_index] == seq2[seq2_index]
        {
            seq1_index += 1;
            seq2_index += 1;
        }
    }
    seq1_index - start
}

// NCBI: s_GetNextNonAffineTback (greedy_align.c:278).
fn s_get_next_non_affine_tback(last_seq2_off: &[Vec<i32>], d: usize, diag: usize) -> (usize, i32) {
    let left = last_seq2_off[d - 1][diag - 1];
    let same = last_seq2_off[d - 1][diag];
    let right = last_seq2_off[d - 1][diag + 1];
    if left > same.max(right) {
        (diag - 1, left)
    } else if same > right {
        (diag, same)
    } else {
        (diag + 1, right)
    }
}

// NCBI: Blast_PrelimEditBlockToGapEditScript (blast_gapalign.c:2482).
fn blast_prelim_edit_block_to_gap_edit_script(
    rev_prelim_tback: &PrelimEditBlock,
    fwd_prelim_tback: &PrelimEditBlock,
) -> GapEditScript {
    let merge_ops = !fwd_prelim_tback.edit_ops.is_empty()
        && !rev_prelim_tback.edit_ops.is_empty()
        && fwd_prelim_tback.edit_ops.last().unwrap().op_type
            == rev_prelim_tback.edit_ops.last().unwrap().op_type;

    let size =
        fwd_prelim_tback.edit_ops.len() + rev_prelim_tback.edit_ops.len() - usize::from(merge_ops);
    let mut esp = GapEditScript::with_capacity(size);

    for op in &rev_prelim_tback.edit_ops {
        esp.push(op.op_type, op.num);
    }

    if fwd_prelim_tback.edit_ops.is_empty() {
        return esp;
    }

    if merge_ops {
        if let Some(last_index) = esp.len().checked_sub(1) {
            let (_, count) = esp.get(last_index).unwrap();
            esp.set_num(
                last_index,
                count + fwd_prelim_tback.edit_ops.last().unwrap().num,
            );
        }
    }

    let start = if merge_ops {
        fwd_prelim_tback.edit_ops.len().saturating_sub(1)
    } else {
        fwd_prelim_tback.edit_ops.len()
    };
    for idx in (0..start).rev() {
        let op = fwd_prelim_tback.edit_ops[idx];
        esp.push(op.op_type, op.num);
    }
    esp
}

// blast-rs: grow-and-retry wrapper around the local `BLAST_GreedyAlign` port;
// NCBI allocates the maximum distance before entry instead.
fn greedy_align_one_side_with_growth(
    query: &[u8],
    subject: &[u8],
    reverse: bool,
    scaled_xdrop: i32,
    scaled_reward: i32,
    scaled_penalty: i32,
    initial_max_dist: usize,
    max_possible_dist: usize,
) -> Option<GreedySideAlignment> {
    let mut max_dist = initial_max_dist;
    loop {
        if let Some((dist, q_ext, s_ext, prelim, seed)) = blast_greedy_align_one_side(
            query,
            subject,
            reverse,
            scaled_xdrop,
            scaled_reward,
            scaled_penalty,
            max_dist,
        ) {
            return Some(GreedySideAlignment {
                dist,
                q_ext,
                s_ext,
                prelim,
                seed,
            });
        }
        if max_dist >= max_possible_dist {
            return None;
        }
        max_dist = (max_dist * 2).min(max_possible_dist);
    }
}

// blast-rs: small value-builder for the two one-sided greedy extensions.
fn greedy_alignment_extents_from_extensions(
    q_seed: usize,
    s_seed: usize,
    reward: i32,
    penalty: i32,
    left: &GreedySideAlignment,
    right: &GreedySideAlignment,
) -> GreedyAlignmentExtents {
    let total_dist = right.dist + left.dist;
    GreedyAlignmentExtents {
        score: ((right.q_ext + right.s_ext + left.q_ext + left.s_ext) as i32 * reward) / 2
            - total_dist * (reward - penalty),
        query_start: q_seed.saturating_sub(left.q_ext),
        query_end: q_seed + right.q_ext,
        subject_start: s_seed.saturating_sub(left.s_ext),
        subject_end: s_seed + right.s_ext,
    }
}

// blast-rs: owns the Rust edit-script merge and final C-shaped gap reduction.
fn greedy_edit_script_from_extensions(
    left: &GreedySideAlignment,
    right: &GreedySideAlignment,
    query: &[u8],
    subject: &[u8],
    extents: GreedyAlignmentExtents,
) -> GapEditScript {
    let mut merged = blast_prelim_edit_block_to_gap_edit_script(&left.prelim, &right.prelim);
    s_reduce_gaps(
        &mut merged,
        &query[extents.query_start..extents.query_end],
        &subject[extents.subject_start..extents.subject_end],
    );
    merged
}

// NCBI: s_RebuildEditScript (blast_gapalign.c:2635).
fn s_rebuild_edit_script(esp: &mut GapEditScript) {
    let mut j: isize = -1;
    let old_size = esp.size.max(0) as usize;
    for i in 0..old_size {
        if esp.num[i] == 0 {
            continue;
        }
        if j >= 0 && esp.op_type[i] == esp.op_type[j as usize] {
            esp.num[j as usize] += esp.num[i];
        } else if j == -1
            || esp.op_type[i] == GapAlignOpType::Sub
            || esp.op_type[j as usize] == GapAlignOpType::Sub
        {
            j += 1;
            esp.op_type[j as usize] = esp.op_type[i];
            esp.num[j as usize] = esp.num[i];
        } else {
            let d = esp.num[j as usize] - esp.num[i];
            if d > 0 {
                esp.num[j as usize - 1] += esp.num[i];
                esp.num[j as usize] = d;
            } else if d < 0 {
                if j == 0 && i > 0 {
                    esp.op_type[j as usize] = GapAlignOpType::Sub;
                    j += 1;
                } else {
                    esp.num[j as usize - 1] += esp.num[j as usize];
                }
                esp.num[j as usize] = -d;
                esp.op_type[j as usize] = esp.op_type[i];
            } else {
                esp.num[j as usize - 1] += esp.num[j as usize];
                j -= 1;
            }
        }
    }
    esp.size = (j + 1).max(0) as i32;
    let new_size = esp.size as usize;
    esp.op_type.truncate(new_size);
    esp.num.truncate(new_size);
}

// NCBI: s_UpdateEditScript (blast_gapalign.c:2573).
fn s_update_edit_script(esp: &mut GapEditScript, pos: usize, bf: i32, af: i32) {
    // 1-1 port of NCBI `s_UpdateEditScript` (`blast_gapalign.c:2573`). The
    // `op_type[op] = Sub` write must precede the `op_type[pos-1] = Del/Ins`
    // (or `pos+1` for the af branch) write, mirroring the C statement
    // order: when `op == pos-1` (or `op == pos+1` for af), the C code
    // first stamps Sub and then immediately overwrites it with Del/Ins on
    // the very next statement, so Del/Ins wins. Reordering these in Rust
    // (Sub last) silently loses one gap_open per HSP whenever the search
    // window collapses to its neighbour.
    if bf > 0 {
        let mut op = pos as i32;
        let (mut qd, mut sd) = (bf, bf);
        loop {
            op -= 1;
            if op < 0 {
                return;
            }
            match esp.op_type[op as usize] {
                GapAlignOpType::Sub => {
                    qd -= esp.num[op as usize];
                    sd -= esp.num[op as usize];
                }
                GapAlignOpType::Ins => qd -= esp.num[op as usize],
                GapAlignOpType::Del => sd -= esp.num[op as usize],
                _ => {}
            }
            if qd <= 0 && sd <= 0 {
                break;
            }
        }
        esp.num[op as usize] = -qd.max(sd);
        // C: `esp->op_type[op++] = eGapAlignSub;` — stamp Sub at the old
        // op, then advance past it for the zero-out loop.
        esp.op_type[op as usize] = GapAlignOpType::Sub;
        let mut next = op as usize + 1;
        while next < pos.saturating_sub(1) {
            esp.num[next] = 0;
            next += 1;
        }
        esp.num[pos] += bf;
        qd -= sd;
        esp.op_type[pos - 1] = if qd > 0 {
            GapAlignOpType::Del
        } else {
            GapAlignOpType::Ins
        };
        esp.num[pos - 1] = qd.unsigned_abs() as i32;
    }

    if af > 0 {
        let mut op = pos as i32;
        let (mut qd, mut sd) = (af, af);
        loop {
            op += 1;
            if op >= esp.size {
                return;
            }
            match esp.op_type[op as usize] {
                GapAlignOpType::Sub => {
                    qd -= esp.num[op as usize];
                    sd -= esp.num[op as usize];
                }
                GapAlignOpType::Ins => qd -= esp.num[op as usize],
                GapAlignOpType::Del => sd -= esp.num[op as usize],
                _ => {}
            }
            if qd <= 0 && sd <= 0 {
                break;
            }
        }
        esp.num[op as usize] = -qd.max(sd);
        // C: `esp->op_type[op--] = eGapAlignSub;` — stamp Sub at op,
        // then step back through the gap-zeroing loop.
        esp.op_type[op as usize] = GapAlignOpType::Sub;
        let mut prev = op as usize - 1;
        while prev > pos + 1 {
            esp.num[prev] = 0;
            prev -= 1;
        }
        esp.num[pos] += af;
        qd -= sd;
        esp.op_type[pos + 1] = if qd > 0 {
            GapAlignOpType::Del
        } else {
            GapAlignOpType::Ins
        };
        esp.num[pos + 1] = qd.unsigned_abs() as i32;
    }
}

// NCBI: s_ReduceGaps (blast_gapalign.c:2687).
fn s_reduce_gaps(esp: &mut GapEditScript, query: &[u8], subject: &[u8]) {
    let (mut q1, mut s1) = (0usize, 0usize);
    let qf = query.len();
    let sf = subject.len();
    let old_size = esp.size.max(0) as usize;
    for i in 0..old_size {
        if esp.num[i] == 0 {
            continue;
        }
        match esp.op_type[i] {
            GapAlignOpType::Sub => {
                if esp.num[i] >= 12 {
                    let mut nm1 = 1i32;
                    if i > 0 {
                        while nm1 as usize <= q1
                            && nm1 as usize <= s1
                            && query[q1 - nm1 as usize] == subject[s1 - nm1 as usize]
                        {
                            nm1 += 1;
                        }
                    }
                    q1 += esp.num[i] as usize;
                    s1 += esp.num[i] as usize;
                    let mut nm2 = 0i32;
                    if i < old_size - 1 {
                        while q1 + 1 < qf && s1 + 1 < sf {
                            let is_match = query[q1] == subject[s1];
                            q1 += 1;
                            s1 += 1;
                            if is_match {
                                nm2 += 1;
                            } else {
                                break;
                            }
                        }
                    }
                    if nm1 > 1 || nm2 > 0 {
                        s_update_edit_script(esp, i, nm1 - 1, nm2);
                    }
                    q1 = q1.saturating_sub(1);
                    s1 = s1.saturating_sub(1);
                } else {
                    q1 += esp.num[i] as usize;
                    s1 += esp.num[i] as usize;
                }
            }
            GapAlignOpType::Ins => q1 += esp.num[i] as usize,
            GapAlignOpType::Del => s1 += esp.num[i] as usize,
            _ => {}
        }
    }
    s_rebuild_edit_script(esp);

    let (mut q, mut s) = (0usize, 0usize);
    let old_size = esp.size.max(0) as usize;
    for i in 0..old_size {
        if esp.op_type[i] == GapAlignOpType::Sub {
            q += esp.num[i].max(0) as usize;
            s += esp.num[i].max(0) as usize;
            continue;
        }
        if i > 1 && esp.op_type[i] != esp.op_type[i - 2] && esp.num[i - 2] > 0 {
            let mut d = esp.num[i] + esp.num[i - 1] + esp.num[i - 2];
            if d == 3 {
                esp.num[i - 2] = 0;
                esp.num[i - 1] = 2;
                esp.num[i] = 0;
                if esp.op_type[i] == GapAlignOpType::Ins {
                    q += 1;
                } else {
                    s += 1;
                }
            } else if d < 12 {
                let slide = esp.num[i].min(esp.num[i - 2]);
                let middle = esp.num[i - 1].max(0) as usize;
                q = q.saturating_sub(middle);
                s = s.saturating_sub(middle);
                let mut q1 = q;
                let mut s1 = s;
                if esp.op_type[i] == GapAlignOpType::Ins {
                    s = s.saturating_sub(slide.max(0) as usize);
                } else {
                    q = q.saturating_sub(slide.max(0) as usize);
                }
                let mut nm1 = 0i32;
                let mut nm2 = 0i32;
                for _ in 0..middle {
                    if query[q1] == subject[s1] {
                        nm1 += 1;
                    }
                    if query[q] == subject[s] {
                        nm2 += 1;
                    }
                    q1 += 1;
                    s1 += 1;
                    q += 1;
                    s += 1;
                }
                for _ in 0..slide.max(0) as usize {
                    if query[q] == subject[s] {
                        nm2 += 1;
                    }
                    q += 1;
                    s += 1;
                }
                d = slide;
                if nm2 >= nm1 - d {
                    esp.num[i - 2] -= d;
                    esp.num[i - 1] += d;
                    esp.num[i] -= d;
                } else {
                    q = q1;
                    s = s1;
                }
            }
        }
        match esp.op_type[i] {
            GapAlignOpType::Ins => q += esp.num[i].max(0) as usize,
            GapAlignOpType::Del => s += esp.num[i].max(0) as usize,
            _ => {}
        }
    }
    s_rebuild_edit_script(esp);
}

// NCBI: BLAST_GreedyAlign (greedy_align.c:348), specialized to one side.
// naming: This private helper implements one directional half of the C routine;
// the public `greedy_align` wrapper supplies the full two-sided alignment.
fn blast_greedy_align_one_side(
    seq1: &[u8],
    seq2: &[u8],
    reverse: bool,
    xdrop_threshold: i32,
    match_cost: i32,
    mismatch_cost: i32,
    max_dist: usize,
) -> Option<(i32, usize, usize, PrelimEditBlock, GreedySeed)> {
    let len1 = seq1.len();
    let len2 = seq2.len();
    let diag_origin = max_dist + 2;
    let width = max_dist * 2 + 6;
    let xdrop_offset =
        ((xdrop_threshold + match_cost / 2) / (match_cost + mismatch_cost) + 1).max(1) as usize;

    let mut last_seq2_off = vec![vec![INVALID_OFFSET; width]; max_dist + 2];
    let mut max_score = vec![0i32; max_dist + xdrop_offset + 1];

    let index = s_find_first_mismatch(seq1, seq2, len1, len2, 0, 0, reverse);
    let mut seq1_align_len = index;
    let mut seq2_align_len = index;
    let mut seed = GreedySeed {
        start_q: 0,
        start_s: 0,
        match_length: index,
    };
    if index == len1 || index == len2 {
        let mut ops = PrelimEditBlock::default();
        ops.add(GapAlignOpType::Sub, index as i32);
        return Some((0, seq1_align_len, seq2_align_len, ops, seed));
    }

    last_seq2_off[0][diag_origin] = index as i32;
    max_score[xdrop_offset] = (index as i32) * match_cost;

    let mut best_dist = 0usize;
    let mut best_diag = 0usize;
    let mut diag_lower = diag_origin - 1;
    let mut diag_upper = diag_origin + 1;
    let mut end1_reached = false;
    let mut end2_reached = false;
    let mut converged = false;

    for d in 1..=max_dist {
        let mut curr_extent = 0i32;
        let mut curr_seq2_index = 0i32;
        let mut curr_diag = 0usize;
        let tmp_diag_lower = diag_lower;
        let tmp_diag_upper = diag_upper;

        last_seq2_off[d - 1][diag_lower - 1] = INVALID_OFFSET;
        last_seq2_off[d - 1][diag_lower] = INVALID_OFFSET;
        last_seq2_off[d - 1][diag_upper] = INVALID_OFFSET;
        last_seq2_off[d - 1][diag_upper + 1] = INVALID_OFFSET;

        let numerator = max_score[d] + (match_cost + mismatch_cost) * d as i32 - xdrop_threshold;
        let xdrop_score = ceil_div_i32(numerator, match_cost / 2);

        for k in tmp_diag_lower..=tmp_diag_upper {
            let mut seq2_index = last_seq2_off[d - 1][k + 1].max(last_seq2_off[d - 1][k]) + 1;
            seq2_index = seq2_index.max(last_seq2_off[d - 1][k - 1]);
            let mut seq1_index = seq2_index + k as i32 - diag_origin as i32;

            if seq2_index < 0 || seq1_index + seq2_index < xdrop_score {
                if k == diag_lower {
                    diag_lower += 1;
                } else {
                    last_seq2_off[d][k] = INVALID_OFFSET;
                }
                continue;
            }
            diag_upper = k;

            let matched = s_find_first_mismatch(
                seq1,
                seq2,
                len1,
                len2,
                seq1_index as usize,
                seq2_index as usize,
                reverse,
            );
            if matched > seed.match_length {
                seed.start_q = seq1_index as usize;
                seed.start_s = seq2_index as usize;
                seed.match_length = matched;
            }
            seq1_index += matched as i32;
            seq2_index += matched as i32;

            last_seq2_off[d][k] = seq2_index;
            if seq1_index + seq2_index > curr_extent {
                curr_extent = seq1_index + seq2_index;
                curr_seq2_index = seq2_index;
                curr_diag = k;
            }

            if seq2_index as usize == len2 {
                diag_lower = k + 1;
                end2_reached = true;
            }
            if seq1_index as usize == len1 {
                diag_upper = k - 1;
                end1_reached = true;
            }
        }

        let curr_score = curr_extent * (match_cost / 2) - d as i32 * (match_cost + mismatch_cost);
        if curr_score >= max_score[d - 1 + xdrop_offset] {
            max_score[d + xdrop_offset] = curr_score;
            best_dist = d;
            best_diag = curr_diag;
            seq2_align_len = curr_seq2_index as usize;
            seq1_align_len = (curr_seq2_index + curr_diag as i32 - diag_origin as i32) as usize;
        } else {
            max_score[d + xdrop_offset] = max_score[d - 1 + xdrop_offset];
        }

        if diag_lower > diag_upper {
            converged = true;
            break;
        }
        if !end2_reached {
            diag_lower -= 1;
        }
        if !end1_reached {
            diag_upper += 1;
        }
    }

    if !converged {
        return None;
    }

    let mut ops = PrelimEditBlock::default();
    let mut d = best_dist;
    let mut seq2_index = seq2_align_len as i32;
    while d > 0 {
        let (new_diag, new_seq2_index) = s_get_next_non_affine_tback(&last_seq2_off, d, best_diag);
        if new_diag == best_diag {
            ops.add(GapAlignOpType::Sub, seq2_index - new_seq2_index);
        } else if new_diag < best_diag {
            ops.add(GapAlignOpType::Sub, seq2_index - new_seq2_index);
            ops.add(GapAlignOpType::Ins, 1);
        } else {
            ops.add(GapAlignOpType::Sub, seq2_index - new_seq2_index - 1);
            ops.add(GapAlignOpType::Del, 1);
        }
        d -= 1;
        best_diag = new_diag;
        seq2_index = new_seq2_index;
    }
    ops.add(GapAlignOpType::Sub, last_seq2_off[0][diag_origin].max(0));

    Some((best_dist as i32, seq1_align_len, seq2_align_len, ops, seed))
}

/// blast-rs: Seeded non-affine greedy alignment helper; not a direct NCBI C
/// port.
///
/// This is the exact megablast path for `gap_open == 0 && gap_extend == 0`.
pub fn greedy_align_with_seed(
    query: &[u8],
    subject: &[u8],
    q_seed: usize,
    s_seed: usize,
    reward: i32,
    penalty: i32,
    x_dropoff: i32,
) -> Option<(i32, usize, usize, usize, usize, GapEditScript, usize, usize)> {
    let scaled = blast_greedy_align_scale_scores(reward, penalty, x_dropoff);
    let (left, right) =
        blast_greedy_align_extend_both_sides(query, subject, q_seed, s_seed, scaled)?;
    let extents =
        greedy_alignment_extents_from_extensions(q_seed, s_seed, reward, penalty, &left, &right);
    let (adjusted_seed_q, adjusted_seed_s) =
        blast_greedy_align_adjust_seed(q_seed, s_seed, extents, left.seed, right.seed);
    let merged = blast_greedy_align_build_traceback(&left, &right, query, subject, extents);

    Some((
        extents.score,
        extents.query_start,
        extents.query_end,
        extents.subject_start,
        extents.subject_end,
        merged,
        adjusted_seed_q,
        adjusted_seed_s,
    ))
}

/// NCBI: greedy nucleotide score-scaling setup.
fn blast_greedy_align_scale_scores(
    reward: i32,
    penalty: i32,
    x_dropoff: i32,
) -> GreedyScaledScores {
    GreedyScaledScores {
        reward: if reward % 2 == 1 { reward * 2 } else { reward },
        penalty: if reward % 2 == 1 {
            -penalty * 2
        } else {
            -penalty
        },
        x_dropoff: if reward % 2 == 1 {
            x_dropoff * 2
        } else {
            x_dropoff
        },
    }
}

/// NCBI: right and left `BLAST_GreedyAlign` extension passes.
fn blast_greedy_align_extend_both_sides(
    query: &[u8],
    subject: &[u8],
    q_seed: usize,
    s_seed: usize,
    scaled: GreedyScaledScores,
) -> Option<(GreedySideAlignment, GreedySideAlignment)> {
    let initial_max_dist = initial_greedy_max_dist(subject.len());
    let max_possible_dist = query
        .len()
        .saturating_add(subject.len())
        .max(initial_max_dist);

    let right = greedy_align_one_side_with_growth(
        &query[q_seed..],
        &subject[s_seed..],
        false,
        scaled.x_dropoff,
        scaled.reward,
        scaled.penalty,
        initial_max_dist,
        max_possible_dist,
    )?;
    let left = greedy_align_one_side_with_growth(
        &query[..q_seed],
        &subject[..s_seed],
        true,
        scaled.x_dropoff,
        scaled.reward,
        scaled.penalty,
        initial_max_dist,
        max_possible_dist,
    )?;
    Some((left, right))
}

/// NCBI: seed adjustment after the two greedy extensions are merged.
fn blast_greedy_align_adjust_seed(
    q_seed: usize,
    s_seed: usize,
    extents: GreedyAlignmentExtents,
    left_seed: GreedySeed,
    right_seed: GreedySeed,
) -> (usize, usize) {
    adjusted_greedy_seed(
        q_seed,
        s_seed,
        extents.query_start,
        extents.subject_start,
        extents.query_end,
        extents.subject_end,
        left_seed,
        right_seed,
    )
}

/// NCBI: `GapEditScript` construction from merged greedy extensions.
fn blast_greedy_align_build_traceback(
    left: &GreedySideAlignment,
    right: &GreedySideAlignment,
    query: &[u8],
    subject: &[u8],
    extents: GreedyAlignmentExtents,
) -> GapEditScript {
    greedy_edit_script_from_extensions(left, right, query, subject, extents)
}

// blast-rs: choose a stable seed inside the merged greedy alignment extents.
fn adjusted_greedy_seed(
    q_seed: usize,
    s_seed: usize,
    query_start: usize,
    subject_start: usize,
    query_end: usize,
    subject_end: usize,
    left_seed: GreedySeed,
    right_seed: GreedySeed,
) -> (usize, usize) {
    let mut right_q = q_seed + right_seed.start_q;
    let mut right_s = s_seed + right_seed.start_s;
    let mut right_len = 0usize;
    if right_q < query_end && right_s < subject_end {
        right_len = (query_end - right_q)
            .min(subject_end - right_s)
            .min(right_seed.match_length)
            / 2;
    } else {
        right_q = q_seed;
        right_s = s_seed;
    }

    let mut left_q = q_seed.saturating_sub(left_seed.start_q);
    let mut left_s = s_seed.saturating_sub(left_seed.start_s);
    let mut left_len = 0usize;
    if left_q > query_start && left_s > subject_start {
        left_len = (left_q - query_start)
            .min(left_s - subject_start)
            .min(left_seed.match_length)
            / 2;
    } else {
        left_q = q_seed;
        left_s = s_seed;
    }

    if right_len > left_len {
        (right_q + right_len, right_s + right_len)
    } else {
        (
            left_q.saturating_sub(left_len),
            left_s.saturating_sub(left_len),
        )
    }
}

/// blast-rs: Rust convenience wrapper that preserves the older tuple return
/// shape by hiding the adjusted seed returned by `greedy_align_with_seed`.
///
/// `subject` must be a DECODED BLASTNA sequence (one base per byte, values
/// `< 16`, no NCBI2na 2-bit packing and no sentinels). NCBI's `s_FindFirstMismatch`
/// (`greedy_align.c:315`) has a separate compressed-subject (`rem` / `fence_hit`)
/// branch for the packed path; the wired Rust callers always decode the subject
/// first, so that branch is intentionally not ported here.
pub fn greedy_align(
    query: &[u8],
    subject: &[u8],
    q_seed: usize,
    s_seed: usize,
    reward: i32,
    penalty: i32,
    x_dropoff: i32,
) -> Option<(i32, usize, usize, usize, usize, GapEditScript)> {
    // The greedy extender assumes a decoded BLASTNA subject (see doc above).
    // NCBI's compressed `rem`/`fence_hit` mismatch branch is deliberately not
    // ported; this guards the invariant in debug builds.
    debug_assert!(
        subject.iter().all(|&b| b < 16),
        "greedy expects a decoded BLASTNA subject (no NCBI2na packing / sentinels)"
    );
    let (score, query_start, query_end, subject_start, subject_end, edit_script, _, _) =
        greedy_align_with_seed(query, subject, q_seed, s_seed, reward, penalty, x_dropoff)?;
    Some((
        score,
        query_start,
        query_end,
        subject_start,
        subject_end,
        edit_script,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn translated_mbspace_pool_allocates_refreshes_and_frees_like_c() {
        let mut pool = mbspace_new(3);
        assert_eq!(pool.space_allocated_per_chunk(), vec![MBSPACE_MIN_SPACE]);

        {
            let slice = s_get_mb_space(Some(&mut pool), 2).expect("first allocation");
            assert_eq!(slice.len(), 2);
            slice[0].match_off = 17;
        }
        assert_eq!(pool.space_used_per_chunk(), vec![2]);

        {
            let slice = s_get_mb_space(Some(&mut pool), MBSPACE_MIN_SPACE as i32)
                .expect("spillover allocation");
            assert_eq!(slice.len(), MBSPACE_MIN_SPACE);
            slice[0].insert_off = 9;
        }
        assert_eq!(
            pool.space_allocated_per_chunk(),
            vec![MBSPACE_MIN_SPACE, MBSPACE_MIN_SPACE]
        );
        assert_eq!(pool.space_used_per_chunk(), vec![2, MBSPACE_MIN_SPACE]);

        s_refresh_mb_space(Some(&mut pool));
        assert_eq!(pool.space_used_per_chunk(), vec![0, 0]);
        assert!(s_get_mb_space(Some(&mut pool), -1).is_none());

        let mut slot = Some(pool);
        mbspace_free(&mut slot);
        assert!(slot.is_none());
    }

    #[test]
    fn blast_greedy_align_mem_alloc_name_matched_non_affine_layout() {
        let mem = s_blast_greedy_align_mem_alloc(1, -3, 0, 0, 5, 10).expect("greedy mem");
        assert_eq!(mem.max_dist, 5);
        assert_eq!(mem.xdrop, 10);
        assert_eq!(
            mem.last_seq2_off
                .as_ref()
                .map(|rows| (rows.len(), rows[0].len())),
            Some((7, 16))
        );
        assert!(mem.last_seq2_off_affine.is_none());
        assert!(mem.diag_bounds.is_empty());
        assert_eq!(mem.max_score.len(), 8);
    }

    #[test]
    fn blast_greedy_align_mem_alloc_name_matched_affine_layout() {
        let mem = s_blast_greedy_align_mem_alloc(2, -3, 5, 2, 4, 20).expect("affine greedy mem");
        assert!(mem.last_seq2_off.is_none());
        assert_eq!(
            mem.last_seq2_off_affine
                .as_ref()
                .map(|rows| (rows.len(), rows[0].len())),
            Some((14, 14))
        );
        assert_eq!(mem.diag_bounds.len(), 2 * (12 + 1 + 8));
        assert_eq!(mem.max_score.len(), 35);
    }

    #[test]
    fn test_greedy_perfect_match() {
        let q = vec![0u8, 1, 2, 3, 0, 1, 2, 3];
        let s = vec![0u8, 1, 2, 3, 0, 1, 2, 3];
        let result = greedy_align(&q, &s, 4, 4, 1, -2, 20);
        assert!(result.is_some());
        let (score, qs, qe, ss, se, esp) = result.unwrap();
        assert_eq!(score, 8);
        assert_eq!((qs, qe), (0, 8));
        assert_eq!((ss, se), (0, 8));
        assert_eq!(esp.ops_vec(), vec![(GapAlignOpType::Sub, 8)]);
    }

    #[test]
    fn test_greedy_gap_in_query() {
        let q: Vec<u8> = [0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3].to_vec();
        let s: Vec<u8> = [
            0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 0, 0, 0, 0, 1, 2, 3, 0, 1, 2, 3,
        ]
        .to_vec();
        let result = greedy_align(&q, &s, 10, 10, 1, -2, 20).expect("greedy gap");
        assert!(result.0 > 0);
        assert!(!result.5.is_empty());
        assert!(result.2 > result.1);
        assert!(result.4 > result.3);
    }
}
