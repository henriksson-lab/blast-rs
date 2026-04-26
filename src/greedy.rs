//! NCBI-style greedy nucleotide alignment for megablast.
//!
//! This module ports the non-affine greedy path from
//! `greedy_align.c`/`BLAST_GreedyGappedAlignment` that NCBI uses when
//! `gap_open == 0 && gap_extend == 0`.

use crate::gapinfo::{GapAlignOpType, GapEditScript};
use crate::sequence::blastna_to_iupac;

const INVALID_OFFSET: i32 = -2;
const MAX_DBSEQ_LEN: usize = 5_000_000;
const GREEDY_MAX_COST_FRACTION: usize = 2;
const GREEDY_MAX_COST: usize = 1000;

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

fn trace_reduce_target() -> bool {
    std::env::var_os("BLAST_RS_TRACE_REDUCE").is_some()
}

fn trace_reduce_script(esp: &GapEditScript, query: &[u8], subject: &[u8]) -> bool {
    let Some(raw) = std::env::var_os("BLAST_RS_TRACE_REDUCE") else {
        return false;
    };
    if raw == "all" {
        return true;
    }
    let raw = raw.to_string_lossy();
    if let Some((q_range, s_range)) = raw.split_once(',') {
        if let (Some((q_lo, q_hi)), Some((s_lo, s_hi))) =
            (parse_trace_range(q_range), parse_trace_range(s_range))
        {
            return query.len() >= q_lo
                && query.len() <= q_hi
                && subject.len() >= s_lo
                && subject.len() <= s_hi;
        }
    }
    let (align_len, num_ident, gap_opens) = esp.count_identities(query, subject);
    (query.len() >= 635 && query.len() <= 645)
        && (subject.len() >= 645 && subject.len() <= 655)
        && align_len >= 650
        && num_ident >= 580
        && gap_opens >= 6
}

fn parse_trace_range(raw: &str) -> Option<(usize, usize)> {
    let (lo, hi) = raw.split_once("..")?;
    Some((lo.parse().ok()?, hi.parse().ok()?))
}

fn trace_btop(esp: &GapEditScript, query: &[u8], subject: &[u8]) -> String {
    let (qseq, sseq) = esp.render_alignment(query, subject, blastna_to_iupac);
    let mut out = String::new();
    let mut matches = 0usize;
    for (q, s) in qseq.bytes().zip(sseq.bytes()) {
        if q == s && q != b'-' {
            matches += 1;
        } else {
            if matches > 0 {
                out.push_str(&matches.to_string());
                matches = 0;
            }
            out.push(q as char);
            out.push(s as char);
        }
    }
    if matches > 0 {
        out.push_str(&matches.to_string());
    }
    out
}

fn trace_encoded_sequence(seq: &[u8]) -> String {
    let mut out = String::with_capacity(seq.len() * 2);
    for base in seq {
        out.push_str(&format!("{base:02x}"));
    }
    out
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

fn find_first_mismatch(
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

#[inline]
fn unpack_packed_base(seq2: &[u8], base_index: usize) -> u8 {
    let byte = seq2[base_index / 4];
    let shift = 3 - (base_index % 4);
    (byte >> (shift * 2)) & 0x03
}

fn find_first_mismatch_packed(
    seq1: &[u8],
    seq2_packed: &[u8],
    len1: usize,
    len2: usize,
    mut seq1_index: usize,
    mut seq2_index: usize,
    reverse: bool,
    rem: usize,
) -> usize {
    let start = seq1_index;
    if reverse {
        while seq1_index < len1
            && seq2_index < len2
            && seq1[len1 - 1 - seq1_index]
                == unpack_packed_base(seq2_packed, len2 - 1 - seq2_index)
        {
            seq1_index += 1;
            seq2_index += 1;
        }
    } else {
        while seq1_index < len1
            && seq2_index < len2
            && seq1[seq1_index]
                == unpack_packed_base(seq2_packed, seq2_index + rem)
        {
            seq1_index += 1;
            seq2_index += 1;
        }
    }
    seq1_index - start
}

fn get_next_non_affine_tback(last_seq2_off: &[Vec<i32>], d: usize, diag: usize) -> (usize, i32) {
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

fn prelim_to_gap_edit_script(
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
        if let Some(last) = esp.ops.last_mut() {
            last.1 += fwd_prelim_tback.edit_ops.last().unwrap().num;
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

fn greedy_prelim_hsp_from_extensions(
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

fn trace_greedy_seed_choice(
    q_seed: usize,
    s_seed: usize,
    extents: GreedyAlignmentExtents,
    left_seed: GreedySeed,
    right_seed: GreedySeed,
    adjusted_seed_q: usize,
    adjusted_seed_s: usize,
) {
    if std::env::var_os("BLAST_RS_TRACE_GREEDY_SEED").is_none() {
        return;
    }
    let right_q = q_seed + right_seed.start_q;
    let right_s = s_seed + right_seed.start_s;
    let right_len = if right_q < extents.query_end && right_s < extents.subject_end {
        (extents.query_end - right_q)
            .min(extents.subject_end - right_s)
            .min(right_seed.match_length)
            / 2
    } else {
        0
    };
    let left_q = q_seed.saturating_sub(left_seed.start_q);
    let left_s = s_seed.saturating_sub(left_seed.start_s);
    let left_len = if left_q > extents.query_start && left_s > extents.subject_start {
        (left_q - extents.query_start)
            .min(left_s - extents.subject_start)
            .min(left_seed.match_length)
            / 2
    } else {
        0
    };
    eprintln!(
        "[rs-gseed] q_seed={} s_seed={} q_box={}..{} s_box={}..{} left=({},{},{}) right=({},{},{}) left_valid={} right_valid={} adjusted=({}, {}) score={}",
        q_seed,
        s_seed,
        extents.query_start,
        extents.query_end,
        extents.subject_start,
        extents.subject_end,
        left_seed.start_q,
        left_seed.start_s,
        left_seed.match_length,
        right_seed.start_q,
        right_seed.start_s,
        right_seed.match_length,
        left_len,
        right_len,
        adjusted_seed_q,
        adjusted_seed_s,
        extents.score
    );
}

fn merge_greedy_prelim_ops(
    left: &GreedySideAlignment,
    right: &GreedySideAlignment,
    query: &[u8],
    subject: &[u8],
    extents: GreedyAlignmentExtents,
    q_seed: usize,
    s_seed: usize,
) -> GapEditScript {
    let mut merged = prelim_to_gap_edit_script(&left.prelim, &right.prelim);
    if trace_reduce_target() {
        let (align_len, num_ident, gap_opens) = merged.count_identities(
            &query[extents.query_start..extents.query_end],
            &subject[extents.subject_start..extents.subject_end],
        );
        if std::env::var_os("BLAST_RS_TRACE_TBACK").is_some()
            && extents.query_start == 692
            && extents.query_end == 1321
            && extents.subject_start == 17655079
            && extents.subject_end == 17655712
        {
            eprintln!(
                "[rs-split] q_seed={q_seed} s_seed={s_seed} q_ext_l={} s_ext_l={} q_ext_r={} s_ext_r={} left={:?} right={:?}",
                left.q_ext,
                left.s_ext,
                right.q_ext,
                right.s_ext,
                left.prelim.edit_ops,
                right.prelim.edit_ops
            );
        }
        eprintln!(
            "[reduce-trace] candidate q={}..{} s={}..{} score={} qlen={} slen={} align_len={} ident={} gaps={} ops={:?}",
            extents.query_start,
            extents.query_end,
            extents.subject_start,
            extents.subject_end,
            extents.score,
            extents.query_end.saturating_sub(extents.query_start),
            extents.subject_end.saturating_sub(extents.subject_start),
            align_len,
            num_ident,
            gap_opens,
            merged.ops
        );
    }
    reduce_gaps(
        &mut merged,
        &query[extents.query_start..extents.query_end],
        &subject[extents.subject_start..extents.subject_end],
    );
    merged
}

fn rebuild_edit_script(esp: &mut GapEditScript) {
    let len = esp.ops.len();
    let mut j: isize = -1;
    for i in 0..len {
        let (op, count) = esp.ops[i];
        if count == 0 {
            continue;
        }
        if j >= 0 && op == esp.ops[j as usize].0 {
            esp.ops[j as usize].1 += count;
        } else if j == -1
            || op == GapAlignOpType::Sub
            || esp.ops[j as usize].0 == GapAlignOpType::Sub
        {
            j += 1;
            esp.ops[j as usize] = (op, count);
        } else {
            let d = esp.ops[j as usize].1 - count;
            if d > 0 {
                esp.ops[j as usize - 1].1 += count;
                esp.ops[j as usize].1 = d;
            } else if d < 0 {
                if j == 0 && i > 0 {
                    esp.ops[j as usize].0 = GapAlignOpType::Sub;
                    j += 1;
                } else {
                    esp.ops[j as usize - 1].1 += esp.ops[j as usize].1;
                }
                esp.ops[j as usize] = (op, -d);
            } else {
                esp.ops[j as usize - 1].1 += esp.ops[j as usize].1;
                j -= 1;
            }
        }
    }
    esp.ops.truncate((j + 1).max(0) as usize);
}

fn update_edit_script(esp: &mut GapEditScript, pos: usize, bf: i32, af: i32) {
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
            match esp.ops[op as usize].0 {
                GapAlignOpType::Sub => {
                    qd -= esp.ops[op as usize].1;
                    sd -= esp.ops[op as usize].1;
                }
                GapAlignOpType::Ins => qd -= esp.ops[op as usize].1,
                GapAlignOpType::Del => sd -= esp.ops[op as usize].1,
                _ => {}
            }
            if qd <= 0 && sd <= 0 {
                break;
            }
        }
        esp.ops[op as usize].1 = -qd.max(sd);
        // C: `esp->op_type[op++] = eGapAlignSub;` — stamp Sub at the old
        // op, then advance past it for the zero-out loop.
        esp.ops[op as usize].0 = GapAlignOpType::Sub;
        let mut next = op as usize + 1;
        while next < pos.saturating_sub(1) {
            esp.ops[next].1 = 0;
            next += 1;
        }
        esp.ops[pos].1 += bf;
        qd -= sd;
        if pos > 0 {
            esp.ops[pos - 1].0 = if qd > 0 {
                GapAlignOpType::Del
            } else {
                GapAlignOpType::Ins
            };
            esp.ops[pos - 1].1 = qd.unsigned_abs() as i32;
        }
    }

    if af > 0 {
        let mut op = pos as i32;
        let (mut qd, mut sd) = (af, af);
        loop {
            op += 1;
            if op as usize >= esp.ops.len() {
                return;
            }
            match esp.ops[op as usize].0 {
                GapAlignOpType::Sub => {
                    qd -= esp.ops[op as usize].1;
                    sd -= esp.ops[op as usize].1;
                }
                GapAlignOpType::Ins => qd -= esp.ops[op as usize].1,
                GapAlignOpType::Del => sd -= esp.ops[op as usize].1,
                _ => {}
            }
            if qd <= 0 && sd <= 0 {
                break;
            }
        }
        esp.ops[op as usize].1 = -qd.max(sd);
        // C: `esp->op_type[op--] = eGapAlignSub;` — stamp Sub at op,
        // then step back through the gap-zeroing loop.
        esp.ops[op as usize].0 = GapAlignOpType::Sub;
        let mut prev = op as usize - 1;
        while prev > pos + 1 {
            esp.ops[prev].1 = 0;
            prev -= 1;
        }
        esp.ops[pos].1 += af;
        qd -= sd;
        if pos + 1 < esp.ops.len() {
            esp.ops[pos + 1].0 = if qd > 0 {
                GapAlignOpType::Del
            } else {
                GapAlignOpType::Ins
            };
            esp.ops[pos + 1].1 = qd.unsigned_abs() as i32;
        }
    }
}

fn reduce_gaps(esp: &mut GapEditScript, query: &[u8], subject: &[u8]) {
    if std::env::var_os("BLAST_RS_SKIP_REDUCE_GAPS").is_some() {
        return;
    }
    let trace_this = trace_reduce_script(esp, query, subject);
    if trace_this {
        eprintln!("[reduce-trace] begin ops={:?}", esp.ops);
        eprintln!(
            "[reduce-trace] begin qhex={} shex={}",
            trace_encoded_sequence(query),
            trace_encoded_sequence(subject)
        );
    }
    let (mut q1, mut s1) = (0usize, 0usize);
    let qf = query.len();
    let sf = subject.len();
    for i in 0..esp.ops.len() {
        if esp.ops[i].1 == 0 {
            continue;
        }
        match esp.ops[i].0 {
            GapAlignOpType::Sub => {
                if esp.ops[i].1 >= 12 {
                    let mut nm1 = 1i32;
                    if i > 0 {
                        while nm1 as usize <= q1
                            && nm1 as usize <= s1
                            && query[q1 - nm1 as usize] == subject[s1 - nm1 as usize]
                        {
                            nm1 += 1;
                        }
                    }
                    q1 += esp.ops[i].1 as usize;
                    s1 += esp.ops[i].1 as usize;
                    let mut nm2 = 0i32;
                    if i + 1 < esp.ops.len() {
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
                        if trace_this {
                            eprintln!(
                                "[reduce-trace] first-pass i={} op={:?} nm1={} nm2={} before={:?}",
                                i, esp.ops[i], nm1, nm2, esp.ops
                            );
                        }
                        update_edit_script(esp, i, nm1 - 1, nm2);
                        if trace_this {
                            eprintln!("[reduce-trace] first-pass after={:?}", esp.ops);
                        }
                    }
                    q1 = q1.saturating_sub(1);
                    s1 = s1.saturating_sub(1);
                } else {
                    q1 += esp.ops[i].1 as usize;
                    s1 += esp.ops[i].1 as usize;
                }
            }
            GapAlignOpType::Ins => q1 += esp.ops[i].1 as usize,
            GapAlignOpType::Del => s1 += esp.ops[i].1 as usize,
            _ => {}
        }
    }
    rebuild_edit_script(esp);

    let (mut q, mut s) = (0usize, 0usize);
    for i in 0..esp.ops.len() {
        if esp.ops[i].0 == GapAlignOpType::Sub {
            q += esp.ops[i].1.max(0) as usize;
            s += esp.ops[i].1.max(0) as usize;
            continue;
        }
        if i > 1 && esp.ops[i].0 != esp.ops[i - 2].0 && esp.ops[i - 2].1 > 0 {
            let mut d = esp.ops[i].1 + esp.ops[i - 1].1 + esp.ops[i - 2].1;
            if d == 3 {
                if trace_this {
                    eprintln!(
                        "[reduce-trace] second-pass i={} d=3 before={:?}",
                        i, esp.ops
                    );
                }
                esp.ops[i - 2].1 = 0;
                esp.ops[i - 1].1 = 2;
                esp.ops[i].1 = 0;
                if esp.ops[i].0 == GapAlignOpType::Ins {
                    q += 1;
                } else {
                    s += 1;
                }
                if trace_this {
                    eprintln!("[reduce-trace] second-pass d=3 after={:?}", esp.ops);
                }
            } else if d < 12 {
                let slide = esp.ops[i].1.min(esp.ops[i - 2].1);
                let middle = esp.ops[i - 1].1.max(0) as usize;
                let q_before = q;
                let s_before = s;
                q = q.saturating_sub(middle);
                s = s.saturating_sub(middle);
                let mut q1 = q;
                let mut s1 = s;
                if esp.ops[i].0 == GapAlignOpType::Ins {
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
                if trace_this {
                    eprintln!(
                        "[reduce-trace] second-pass i={} tri={:?}/{:?}/{:?} q_before={} s_before={} middle={} slide={} nm1={} nm2={} decision={}",
                        i,
                        esp.ops[i - 2],
                        esp.ops[i - 1],
                        esp.ops[i],
                        q_before,
                        s_before,
                        middle,
                        d,
                        nm1,
                        nm2,
                        if nm2 >= nm1 - d { "accept" } else { "reject" }
                    );
                }
                if nm2 >= nm1 - d {
                    esp.ops[i - 2].1 -= d;
                    esp.ops[i - 1].1 += d;
                    esp.ops[i].1 -= d;
                } else {
                    q = q1;
                    s = s1;
                }
            }
        }
        match esp.ops[i].0 {
            GapAlignOpType::Ins => q += esp.ops[i].1.max(0) as usize,
            GapAlignOpType::Del => s += esp.ops[i].1.max(0) as usize,
            _ => {}
        }
    }
    rebuild_edit_script(esp);
    if trace_this {
        let (align_len, num_ident, gap_opens) = esp.count_identities(query, subject);
        eprintln!(
            "[reduce-trace] end ident={} len={} gaps={} ops={:?} btop={}",
            num_ident,
            align_len,
            gap_opens,
            esp.ops,
            trace_btop(esp, query, subject)
        );
    }
}

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

    let index = find_first_mismatch(seq1, seq2, len1, len2, 0, 0, reverse);
    let mut seq1_align_len = index;
    let mut seq2_align_len = index;
    let mut seed = GreedySeed {
        start_q: 0,
        start_s: 0,
        match_length: index,
    };
    let trace_internal = std::env::var_os("BLAST_RS_TRACE_GREEDY_INTERNAL").is_some();

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

            let matched = find_first_mismatch(
                seq1,
                seq2,
                len1,
                len2,
                seq1_index as usize,
                seq2_index as usize,
                reverse,
            );
            if matched > seed.match_length {
                if trace_internal {
                    eprintln!(
                        "[rs-ginternal] reverse={} d={} k={} seq1_index={} seq2_index={} matched={} old_seed=({},{},{}) new_seed=({},{},{})",
                        reverse,
                        d,
                        k,
                        seq1_index,
                        seq2_index,
                        matched,
                        seed.start_q,
                        seed.start_s,
                        seed.match_length,
                        seq1_index,
                        seq2_index,
                        matched
                    );
                }
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
    let trace_tback = std::env::var_os("BLAST_RS_TRACE_TBACK").is_some()
        && ((seq1_align_len == 629 && seq2_align_len == 633)
            || (reverse
                && (500..=650).contains(&seq1_align_len)
                && (500..=650).contains(&seq2_align_len)));
    if trace_tback {
        eprintln!(
            "[rs-tback] start len1={len1} len2={len2} best_dist={best_dist} best_diag={best_diag} seq1_align_len={seq1_align_len} seq2_align_len={seq2_align_len} diag_origin={diag_origin}"
        );
        for row in 1..=best_dist.min(12) {
            let lo = diag_origin.saturating_sub(4);
            let hi = (diag_origin + 4).min(width - 1);
            let vals: Vec<String> = (lo..=hi)
                .map(|diag| format!("{diag}:{}", last_seq2_off[row][diag]))
                .collect();
            eprintln!("[rs-row] d={row} {}", vals.join(" "));
        }
    }
    while d > 0 {
        let left = last_seq2_off[d - 1][best_diag - 1];
        let same = last_seq2_off[d - 1][best_diag];
        let right = last_seq2_off[d - 1][best_diag + 1];
        let (new_diag, new_seq2_index) = get_next_non_affine_tback(&last_seq2_off, d, best_diag);
        if trace_tback && d <= 12 {
            eprintln!(
                "[rs-tback] d={d} diag={best_diag} seq2={seq2_index} left={left} same={same} right={right} -> new_diag={new_diag} new_seq2={new_seq2_index}"
            );
        }
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

    if trace_internal {
        eprintln!(
            "[rs-ginternal] done reverse={} best_dist={} best_diag={} q_align_len={} s_align_len={} seed=({},{},{})",
            reverse,
            best_dist,
            best_diag,
            seq1_align_len,
            seq2_align_len,
            seed.start_q,
            seed.start_s,
            seed.match_length
        );
    }

    Some((best_dist as i32, seq1_align_len, seq2_align_len, ops, seed))
}

fn blast_greedy_align_one_side_packed_subject(
    seq1: &[u8],
    seq2_packed: &[u8],
    len2: usize,
    reverse: bool,
    rem: usize,
    xdrop_threshold: i32,
    match_cost: i32,
    mismatch_cost: i32,
    max_dist: usize,
) -> Option<(i32, usize, usize, PrelimEditBlock, GreedySeed)> {
    let len1 = seq1.len();
    let diag_origin = max_dist + 2;
    let width = max_dist * 2 + 6;
    let xdrop_offset =
        ((xdrop_threshold + match_cost / 2) / (match_cost + mismatch_cost) + 1).max(1) as usize;

    let mut last_seq2_off = vec![vec![INVALID_OFFSET; width]; max_dist + 2];
    let mut max_score = vec![0i32; max_dist + xdrop_offset + 1];

    let index = find_first_mismatch_packed(seq1, seq2_packed, len1, len2, 0, 0, reverse, rem);
    let mut seq1_align_len = index;
    let mut seq2_align_len = index;
    let mut seed = GreedySeed {
        start_q: 0,
        start_s: 0,
        match_length: index,
    };
    let trace_internal = std::env::var_os("BLAST_RS_TRACE_GREEDY_INTERNAL").is_some();

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

            let matched = find_first_mismatch_packed(
                seq1,
                seq2_packed,
                len1,
                len2,
                seq1_index as usize,
                seq2_index as usize,
                reverse,
                rem,
            );
            if matched > seed.match_length {
                if trace_internal {
                    eprintln!(
                        "[rs-ginternal] reverse={} packed=1 rem={} d={} k={} seq1_index={} seq2_index={} matched={} old_seed=({},{},{}) new_seed=({},{},{})",
                        reverse,
                        rem,
                        d,
                        k,
                        seq1_index,
                        seq2_index,
                        matched,
                        seed.start_q,
                        seed.start_s,
                        seed.match_length,
                        seq1_index,
                        seq2_index,
                        matched
                    );
                }
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
        let (new_diag, new_seq2_index) = get_next_non_affine_tback(&last_seq2_off, d, best_diag);
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

    if trace_internal {
        eprintln!(
            "[rs-ginternal] done reverse={} packed=1 rem={} best_dist={} best_diag={} q_align_len={} s_align_len={} seed=({},{},{})",
            reverse,
            rem,
            best_dist,
            best_diag,
            seq1_align_len,
            seq2_align_len,
            seed.start_q,
            seed.start_s,
            seed.match_length
        );
    }

    Some((best_dist as i32, seq1_align_len, seq2_align_len, ops, seed))
}

fn greedy_align_one_side_with_growth_packed_subject(
    query: &[u8],
    subject_packed: &[u8],
    subject_len: usize,
    reverse: bool,
    rem: usize,
    scaled_xdrop: i32,
    scaled_reward: i32,
    scaled_penalty: i32,
    initial_max_dist: usize,
    max_possible_dist: usize,
) -> Option<GreedySideAlignment> {
    let mut max_dist = initial_max_dist;
    loop {
        if let Some((dist, q_ext, s_ext, prelim, seed)) = blast_greedy_align_one_side_packed_subject(
            query,
            subject_packed,
            subject_len,
            reverse,
            rem,
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

/// Port of the non-affine branch of NCBI `BLAST_GreedyGappedAlignment`.
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
    let scaled_reward = if reward % 2 == 1 { reward * 2 } else { reward };
    let scaled_penalty = if reward % 2 == 1 {
        -penalty * 2
    } else {
        -penalty
    };
    let scaled_xdrop = if reward % 2 == 1 {
        x_dropoff * 2
    } else {
        x_dropoff
    };

    let initial_max_dist = initial_greedy_max_dist(subject.len());
    let max_possible_dist = query.len().saturating_add(subject.len()).max(initial_max_dist);

    let right = greedy_align_one_side_with_growth(
        &query[q_seed..],
        &subject[s_seed..],
        false,
        scaled_xdrop,
        scaled_reward,
        scaled_penalty,
        initial_max_dist,
        max_possible_dist,
    )?;
    let left = greedy_align_one_side_with_growth(
        &query[..q_seed],
        &subject[..s_seed],
        true,
        scaled_xdrop,
        scaled_reward,
        scaled_penalty,
        initial_max_dist,
        max_possible_dist,
    )?;
    let extents = greedy_prelim_hsp_from_extensions(q_seed, s_seed, reward, penalty, &left, &right);
    let (adjusted_seed_q, adjusted_seed_s) = adjusted_greedy_seed(
        q_seed,
        s_seed,
        extents.query_start,
        extents.subject_start,
        extents.query_end,
        extents.subject_end,
        left.seed,
        right.seed,
    );
    trace_greedy_seed_choice(
        q_seed,
        s_seed,
        extents,
        left.seed,
        right.seed,
        adjusted_seed_q,
        adjusted_seed_s,
    );
    let merged = merge_greedy_prelim_ops(&left, &right, query, subject, extents, q_seed, s_seed);

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

#[inline(never)]
pub fn greedy_align_with_seed_packed_subject(
    query: &[u8],
    subject_packed: &[u8],
    subject_len: usize,
    q_seed: usize,
    s_seed: usize,
    reward: i32,
    penalty: i32,
    x_dropoff: i32,
) -> Option<(i32, usize, usize, usize, usize, usize, usize)> {
    let scaled_reward = if reward % 2 == 1 { reward * 2 } else { reward };
    let scaled_penalty = if reward % 2 == 1 {
        -penalty * 2
    } else {
        -penalty
    };
    let scaled_xdrop = if reward % 2 == 1 {
        x_dropoff * 2
    } else {
        x_dropoff
    };

    let initial_max_dist = initial_greedy_max_dist(subject_len);
    let max_possible_dist = query.len().saturating_add(subject_len).max(initial_max_dist);

    let right = greedy_align_one_side_with_growth_packed_subject(
        &query[q_seed..],
        &subject_packed[s_seed / 4..],
        subject_len - s_seed,
        false,
        s_seed % 4,
        scaled_xdrop,
        scaled_reward,
        scaled_penalty,
        initial_max_dist,
        max_possible_dist,
    )?;
    let left = greedy_align_one_side_with_growth_packed_subject(
        &query[..q_seed],
        subject_packed,
        s_seed,
        true,
        0,
        scaled_xdrop,
        scaled_reward,
        scaled_penalty,
        initial_max_dist,
        max_possible_dist,
    )?;
    let extents = greedy_prelim_hsp_from_extensions(q_seed, s_seed, reward, penalty, &left, &right);
    let (adjusted_seed_q, adjusted_seed_s) = adjusted_greedy_seed(
        q_seed,
        s_seed,
        extents.query_start,
        extents.subject_start,
        extents.query_end,
        extents.subject_end,
        left.seed,
        right.seed,
    );
    trace_greedy_seed_choice(
        q_seed,
        s_seed,
        extents,
        left.seed,
        right.seed,
        adjusted_seed_q,
        adjusted_seed_s,
    );

    Some((
        extents.score,
        extents.query_start,
        extents.query_end,
        extents.subject_start,
        extents.subject_end,
        adjusted_seed_q,
        adjusted_seed_s,
    ))
}

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
        (left_q.saturating_sub(left_len), left_s.saturating_sub(left_len))
    }
}

pub fn greedy_align(
    query: &[u8],
    subject: &[u8],
    q_seed: usize,
    s_seed: usize,
    reward: i32,
    penalty: i32,
    x_dropoff: i32,
) -> Option<(i32, usize, usize, usize, usize, GapEditScript)> {
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
    fn test_greedy_perfect_match() {
        let q = vec![0u8, 1, 2, 3, 0, 1, 2, 3];
        let s = vec![0u8, 1, 2, 3, 0, 1, 2, 3];
        let result = greedy_align(&q, &s, 4, 4, 1, -2, 20);
        assert!(result.is_some());
        let (score, qs, qe, ss, se, esp) = result.unwrap();
        assert_eq!(score, 8);
        assert_eq!((qs, qe), (0, 8));
        assert_eq!((ss, se), (0, 8));
        assert_eq!(esp.ops, vec![(GapAlignOpType::Sub, 8)]);
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
        assert!(!result.5.ops.is_empty());
        assert!(result.2 > result.1);
        assert!(result.4 > result.3);
    }
}
