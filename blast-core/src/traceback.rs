//! Rust equivalent of blast_traceback.c — traceback (gapped alignment with full edit script).

use crate::gapinfo::{GapAlignOpType, GapEditScript};

/// Perform traceback alignment between query and subject sequences.
/// This is the full dynamic programming alignment that produces
/// a complete GapEditScript for output formatting.
///
/// Unlike the preliminary gapped alignment (which only computes a score),
/// traceback produces the actual alignment path.
pub fn traceback_align(
    query: &[u8],
    subject: &[u8],
    q_start: usize,
    q_end: usize,
    s_start: usize,
    s_end: usize,
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
) -> (i32, GapEditScript, usize, usize, usize, usize) {
    let q = &query[q_start..q_end];
    let s = &subject[s_start..s_end];
    let m = q.len();
    let n = s.len();

    if m == 0 || n == 0 {
        return (0, GapEditScript::new(), 0, 0, 0, 0);
    }

    let gap_oe = gap_open + gap_extend;

    // DP matrices
    let cols = n + 1;
    let mut h = vec![vec![0i32; cols]; m + 1];
    let mut e = vec![vec![i32::MIN / 2; cols]; m + 1]; // gap in query
    let mut f = vec![vec![i32::MIN / 2; cols]; m + 1]; // gap in subject

    // Traceback direction matrix
    // 0=diag, 1=left(ins), 2=up(del)
    let mut trace = vec![vec![0u8; cols]; m + 1];

    for i in 1..=m {
        for j in 1..=n {
            let match_score = if q[i - 1] == s[j - 1] { reward } else { penalty };
            let diag = h[i - 1][j - 1] + match_score;

            e[i][j] = (h[i - 1][j] - gap_oe).max(e[i - 1][j] - gap_extend);
            f[i][j] = (h[i][j - 1] - gap_oe).max(f[i][j - 1] - gap_extend);

            h[i][j] = diag.max(e[i][j]).max(f[i][j]).max(0);

            if h[i][j] == 0 {
                trace[i][j] = 255; // reset
            } else if h[i][j] == diag {
                trace[i][j] = 0; // diagonal
            } else if h[i][j] == f[i][j] {
                trace[i][j] = 1; // left (gap in subject)
            } else {
                trace[i][j] = 2; // up (gap in query)
            }
        }
    }

    // Find best score position
    let mut best_score = 0;
    let mut best_i = 0;
    let mut best_j = 0;
    for i in 1..=m {
        for j in 1..=n {
            if h[i][j] > best_score {
                best_score = h[i][j];
                best_i = i;
                best_j = j;
            }
        }
    }

    // Traceback from best position
    let mut ops: Vec<(GapAlignOpType, i32)> = Vec::new();
    let mut i = best_i;
    let mut j = best_j;

    while i > 0 && j > 0 && h[i][j] > 0 {
        match trace[i][j] {
            0 => {
                // Diagonal (substitution)
                if let Some(last) = ops.last_mut() {
                    if last.0 == GapAlignOpType::Sub {
                        last.1 += 1;
                    } else {
                        ops.push((GapAlignOpType::Sub, 1));
                    }
                } else {
                    ops.push((GapAlignOpType::Sub, 1));
                }
                i -= 1;
                j -= 1;
            }
            1 => {
                // Left (gap in subject = insertion)
                if let Some(last) = ops.last_mut() {
                    if last.0 == GapAlignOpType::Ins {
                        last.1 += 1;
                    } else {
                        ops.push((GapAlignOpType::Ins, 1));
                    }
                } else {
                    ops.push((GapAlignOpType::Ins, 1));
                }
                j -= 1;
            }
            2 => {
                // Up (gap in query = deletion)
                if let Some(last) = ops.last_mut() {
                    if last.0 == GapAlignOpType::Del {
                        last.1 += 1;
                    } else {
                        ops.push((GapAlignOpType::Del, 1));
                    }
                } else {
                    ops.push((GapAlignOpType::Del, 1));
                }
                i -= 1;
            }
            _ => break,
        }
    }

    // Reverse the ops (we traced backwards)
    ops.reverse();
    let mut esp = GapEditScript::new();
    for (op, count) in ops {
        esp.push(op, count);
    }

    // i and j now point to the start of the local alignment (in 1-based DP coords)
    // Convert to 0-based positions in the original subsequences
    let align_q_start = i; // 0-based in sub-query
    let align_s_start = j;
    let align_q_end = best_i;
    let align_s_end = best_j;

    (best_score, esp, align_q_start, align_q_end, align_s_start, align_s_end)
}

/// Result of a traceback alignment.
pub struct TracebackResult {
    pub score: i32,
    pub edit_script: GapEditScript,
    pub query_start: usize,  // 0-based in original sequence
    pub query_end: usize,
    pub subject_start: usize,
    pub subject_end: usize,
}

/// Higher-level traceback that takes absolute coordinates and returns absolute coordinates.
pub fn traceback_align_abs(
    query: &[u8],
    subject: &[u8],
    seed_q: usize,
    seed_s: usize,
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    margin: usize,
) -> Option<TracebackResult> {
    let q_start = seed_q.saturating_sub(margin);
    let s_start = seed_s.saturating_sub(margin);
    let q_end = (seed_q + margin + 1).min(query.len());
    let s_end = (seed_s + margin + 1).min(subject.len());

    if q_end <= q_start || s_end <= s_start {
        return None;
    }

    let (score, esp, aq_start, aq_end, as_start, as_end) = traceback_align(
        query, subject, q_start, q_end, s_start, s_end,
        reward, penalty, gap_open, gap_extend,
    );

    if score <= 0 {
        return None;
    }

    Some(TracebackResult {
        score,
        edit_script: esp,
        query_start: q_start + aq_start,
        query_end: q_start + aq_end,
        subject_start: s_start + as_start,
        subject_end: s_start + as_end,
    })
}

// ---------------------------------------------------------------------------
// BLAST-style gapped alignment with X-dropoff (port of ALIGN_EX)
// ---------------------------------------------------------------------------

const MININT: i32 = i32::MIN / 2;
const SCRIPT_SUB: u8 = 3;
const SCRIPT_GAP_IN_A: u8 = 0; // gap in query = deletion
const SCRIPT_GAP_IN_B: u8 = 6; // gap in subject = insertion
const SCRIPT_OP_MASK: u8 = 0x07;
const SCRIPT_EXTEND_GAP_A: u8 = 0x10;
const SCRIPT_EXTEND_GAP_B: u8 = 0x40;

#[derive(Clone, Copy)]
struct GapDP { best: i32, best_gap: i32 }

fn script_to_op(s: u8) -> GapAlignOpType {
    match s { SCRIPT_GAP_IN_A => GapAlignOpType::Del, SCRIPT_GAP_IN_B => GapAlignOpType::Ins, _ => GapAlignOpType::Sub }
}

/// One-directional gapped extension with X-dropoff and traceback (port of ALIGN_EX).
fn align_ex(
    a: &[u8], b: &[u8], m: usize, n: usize,
    matrix: &[[i32; 16]; 16],
    gap_open: i32, gap_extend: i32, mut x_dropoff: i32,
    reverse: bool,
) -> (i32, usize, usize, Vec<(GapAlignOpType, i32)>) {
    let gap_oe = gap_open + gap_extend;
    if x_dropoff < gap_oe { x_dropoff = gap_oe; }
    if m == 0 || n == 0 { return (0, 0, 0, Vec::new()); }

    let extra = if gap_extend > 0 { (x_dropoff / gap_extend) as usize + 3 } else { n + 3 };
    // Cap allocation: band never exceeds extra width from diagonal
    let alloc = (extra * 2 + 20).min(n + extra + 10);
    let mut sa = vec![GapDP { best: MININT, best_gap: MININT }; alloc];

    // Row 0 initialization
    sa[0] = GapDP { best: 0, best_gap: -gap_oe };
    let mut score = -gap_oe;
    let mut b_size = 1usize;
    while b_size <= n && score >= -x_dropoff {
        sa[b_size] = GapDP { best: score, best_gap: score - gap_oe };
        score -= gap_extend;
        b_size += 1;
    }

    // Traceback storage: script[a_idx][b_idx - row_start]
    let mut scripts: Vec<Vec<u8>> = Vec::with_capacity(m + 1);
    let mut row_starts: Vec<usize> = Vec::with_capacity(m + 1);
    scripts.push(vec![SCRIPT_GAP_IN_A; b_size]);
    row_starts.push(0);

    let mut best_score = 0i32;
    let mut first_b = 0usize;
    let mut a_off = 0usize;
    let mut b_off = 0usize;

    for ai in 1..=m {
        let a_letter = if reverse { a[m - ai] } else { a[ai] };
        let mrow = &matrix[a_letter as usize & 0x0F];

        let mut row_script = vec![0u8; b_size + extra + 10];
        let mut sc = MININT;
        let mut sgr = MININT; // score_gap_row
        let mut last_b = first_b;

        for bi in first_b..b_size {
            // C code: b_ptr starts at B[first_b_index], then b_ptr += b_increment BEFORE access.
            // Forward: B[b_index + 1]. Reverse: B[N - b_index - 1].
            let b_idx = if reverse {
                n.checked_sub(1 + bi).unwrap_or(usize::MAX)
            } else {
                bi + 1
            };
            if b_idx >= b.len() { break; }
            let b_letter = b[b_idx];
            let sgc = sa[bi].best_gap; // score_gap_col
            let next_sc = sa[bi].best + mrow[b_letter as usize & 0x0F];

            // Best predecessor
            let mut script = SCRIPT_SUB;
            if sc < sgc { script = SCRIPT_GAP_IN_B; sc = sgc; }
            if sc < sgr { script = SCRIPT_GAP_IN_A; sc = sgr; }

            if best_score - sc > x_dropoff {
                if first_b == bi { first_b += 1; }
                sa[bi].best = MININT;
                // Note: C code does NOT reset best_gap here
            } else {
                last_b = bi;
                if sc > best_score { best_score = sc; a_off = ai; b_off = bi; }

                // Update gap scores
                if sgc - gap_extend < sc - gap_oe {
                    sa[bi].best_gap = sc - gap_oe;
                } else {
                    sa[bi].best_gap = sgc - gap_extend;
                    script |= SCRIPT_EXTEND_GAP_B;
                }
                if sgr - gap_extend < sc - gap_oe {
                    sgr = sc - gap_oe;
                } else {
                    sgr -= gap_extend;
                    script |= SCRIPT_EXTEND_GAP_A;
                }
                sa[bi].best = sc;
            }
            sc = next_sc;
            if bi < row_script.len() { row_script[bi] = script; }
        }

        // Handle the last diagonal value that extends one column past the loop bounds.
        // In the C code, this cell would be processed in the next bi iteration.
        // 'sc' now holds next_sc from the last iteration = diagonal for (ai, last_bi+1).
        if sc != MININT && best_score - sc <= x_dropoff && sc > best_score {
            let bi_extra = last_b + 1;
            if bi_extra < alloc {
                best_score = sc;
                a_off = ai;
                b_off = bi_extra;
                sa[bi_extra].best = sc;
                if bi_extra < row_script.len() {
                    row_script[bi_extra] = SCRIPT_SUB;
                }
                if bi_extra >= b_size { b_size = bi_extra + 1; }
                last_b = bi_extra;
            }
        }

        scripts.push(row_script);
        row_starts.push(first_b);

        if first_b >= b_size { break; }

        if last_b < b_size - 1 {
            b_size = last_b + 1;
        } else {
            // Extend band rightward
            while sgr >= best_score - x_dropoff && b_size <= n {
                if b_size >= sa.len() { sa.resize(b_size + 10, GapDP { best: MININT, best_gap: MININT }); }
                sa[b_size] = GapDP { best: sgr, best_gap: sgr - gap_oe };
                sgr -= gap_extend;
                b_size += 1;
            }
        }
        if b_size <= n {
            if b_size >= sa.len() { sa.resize(b_size + 10, GapDP { best: MININT, best_gap: MININT }); }
            sa[b_size] = GapDP { best: MININT, best_gap: MININT };
            b_size += 1;
        }
    }

    // Traceback
    let mut ops: Vec<(GapAlignOpType, i32)> = Vec::new();
    let mut ai = a_off;
    let mut bi = b_off;
    let mut cur_script = SCRIPT_SUB;

    while ai > 0 || bi > 0 {
        if ai >= scripts.len() { break; }
        if bi >= scripts[ai].len() { break; }
        let s = scripts[ai][bi];

        cur_script = match cur_script & SCRIPT_OP_MASK {
            SCRIPT_GAP_IN_A => if s & SCRIPT_EXTEND_GAP_A != 0 { SCRIPT_GAP_IN_A } else { s & SCRIPT_OP_MASK },
            SCRIPT_GAP_IN_B => if s & SCRIPT_EXTEND_GAP_B != 0 { SCRIPT_GAP_IN_B } else { s & SCRIPT_OP_MASK },
            _ => s & SCRIPT_OP_MASK,
        };

        let op = script_to_op(cur_script);
        match cur_script & SCRIPT_OP_MASK {
            SCRIPT_GAP_IN_A => { if bi == 0 { break; } bi -= 1; },
            SCRIPT_GAP_IN_B => { if ai == 0 { break; } ai -= 1; },
            _ => { if ai == 0 || bi == 0 { break; } ai -= 1; bi -= 1; },
        }
        if let Some(last) = ops.last_mut() {
            if last.0 == op { last.1 += 1; } else { ops.push((op, 1)); }
        } else { ops.push((op, 1)); }
    }

    (best_score, a_off, b_off, ops)
}

/// Build full BLASTNA 16x16 scoring matrix matching C engine's BlastScoreBlkNuclMatrixCreate.
fn build_blastna_matrix(reward: i32, penalty: i32) -> [[i32; 16]; 16] {
    const BLASTNA_TO_NCBI4NA: [u8; 16] = [1,2,4,8,5,10,3,12,9,6,14,13,11,7,15,0];
    let mut matrix = [[0i32; 16]; 16];
    // Compute degeneracy
    let mut deg = [0i32; 16];
    for i in 0..4 { deg[i] = 1; }
    for i in 4..16 {
        for j in 0..4 { if BLASTNA_TO_NCBI4NA[i] & BLASTNA_TO_NCBI4NA[j] != 0 { deg[i] += 1; } }
    }
    for i in 0..16 {
        for j in i..16 {
            let s = if BLASTNA_TO_NCBI4NA[i] & BLASTNA_TO_NCBI4NA[j] != 0 {
                let d = deg[j] as f64;
                (((d - 1.0) * penalty as f64 + reward as f64) / d).round() as i32
            } else { penalty };
            matrix[i][j] = s;
            matrix[j][i] = s;
        }
    }
    // Sentinel row/col (gap = index 15)
    for i in 0..16 { matrix[15][i] = i32::MIN / 2; matrix[i][15] = i32::MIN / 2; }
    matrix
}

/// BLAST-style gapped alignment with traceback — extends bidirectionally from seed.
/// Fast score-only gapped extension (no traceback).
/// Port of C engine's Blast_SemiGappedAlign — computes only the score.
pub fn blast_gapped_score_only(
    query: &[u8], subject: &[u8],
    seed_q: usize, seed_s: usize,
    reward: i32, penalty: i32,
    gap_open: i32, gap_extend: i32,
    x_dropoff: i32,
) -> i32 {
    let gap_oe = gap_open + gap_extend;
    let matrix = build_blastna_matrix(reward, penalty);

    // Left extension (score only)
    let score_l = gapped_score_one_dir(
        &query[..seed_q + 1], &subject[..seed_s + 1],
        seed_q + 1, seed_s + 1, &matrix, gap_oe, gap_extend, x_dropoff, true);

    // Right extension (score only)
    let score_r = if seed_q < query.len() - 1 && seed_s < subject.len() - 1 {
        gapped_score_one_dir(
            &query[seed_q..], &subject[seed_s..],
            query.len() - seed_q - 1, subject.len() - seed_s - 1,
            &matrix, gap_oe, gap_extend, x_dropoff, false)
    } else { 0 };

    score_l + score_r
}

/// One-directional score-only gapped extension with X-dropoff.
/// Much faster than align_ex because no traceback storage/recording.
fn gapped_score_one_dir(
    a: &[u8], b: &[u8], m: usize, n: usize,
    matrix: &[[i32; 16]; 16],
    gap_oe: i32, gap_extend: i32, mut x_dropoff: i32,
    reverse: bool,
) -> i32 {
    if x_dropoff < gap_oe { x_dropoff = gap_oe; }
    if m == 0 || n == 0 { return 0; }

    // Cap allocation to x_dropoff band width (like C engine's dp_mem_alloc = 1000)
    // The band never grows wider than x_dropoff/gap_extend + some margin
    let max_band = ((x_dropoff / gap_extend.max(1)) as usize + 10).min(n + 1);
    let mut sa = vec![GapDP { best: MININT, best_gap: MININT }; max_band];

    sa[0] = GapDP { best: 0, best_gap: -gap_oe };
    let mut score = -gap_oe;
    let mut b_size = 1usize;
    while b_size <= n && b_size < sa.len() && score >= -x_dropoff {
        sa[b_size] = GapDP { best: score, best_gap: score - gap_oe };
        score -= gap_extend;
        b_size += 1;
    }

    let mut best_score = 0i32;
    let mut first_b = 0usize;

    for ai in 1..=m {
        let a_letter = if reverse { a[m - ai] } else { a[ai] };
        let mrow = &matrix[a_letter as usize & 0x0F];
        let mut sc = MININT;
        let mut sgr = MININT;
        let mut last_b = first_b;

        for bi in first_b..b_size {
            let b_idx = if reverse {
                n.checked_sub(1 + bi).unwrap_or(usize::MAX)
            } else { bi + 1 };
            if b_idx >= b.len() { break; }
            let b_letter = b[b_idx];

            let sgc = sa[bi].best_gap;
            let next_sc = sa[bi].best + mrow[b_letter as usize & 0x0F];

            // Three-way max
            if sc < sgc { sc = sgc; }
            if sc < sgr { sc = sgr; }

            if best_score - sc > x_dropoff {
                if first_b == bi { first_b += 1; }
                sa[bi].best = MININT;
            } else {
                last_b = bi;
                if sc > best_score { best_score = sc; }
                if sgc - gap_extend < sc - gap_oe { sa[bi].best_gap = sc - gap_oe; }
                else { sa[bi].best_gap = sgc - gap_extend; }
                if sgr - gap_extend < sc - gap_oe { sgr = sc - gap_oe; }
                else { sgr -= gap_extend; }
                sa[bi].best = sc;
            }
            sc = next_sc;
        }

        // Handle last diagonal value
        if sc != MININT && best_score - sc <= x_dropoff && sc > best_score {
            best_score = sc;
        }

        if first_b >= b_size { break; }
        if last_b < b_size - 1 { b_size = last_b + 1; }
        else {
            while sgr >= best_score - x_dropoff && b_size <= n && b_size < sa.len() {
                sa[b_size] = GapDP { best: sgr, best_gap: sgr - gap_oe };
                sgr -= gap_extend;
                b_size += 1;
            }
        }
        if b_size <= n && b_size < sa.len() {
            sa[b_size] = GapDP { best: MININT, best_gap: MININT };
            b_size += 1;
        }
    }
    best_score
}

pub fn blast_gapped_align(
    query: &[u8], subject: &[u8],
    seed_q: usize, seed_s: usize,
    reward: i32, penalty: i32,
    gap_open: i32, gap_extend: i32,
    x_dropoff: i32,
) -> Option<TracebackResult> {
    // Build full BLASTNA scoring matrix matching C engine (handles ambiguous bases)
    let matrix = build_blastna_matrix(reward, penalty);

    // Left extension (reverse)
    let (score_l, ql, sl, left_ops) = align_ex(
        &query[..seed_q + 1], &subject[..seed_s + 1],
        seed_q + 1, seed_s + 1,
        &matrix, gap_open, gap_extend, x_dropoff, true,
    );
    let q_start = seed_q + 1 - ql;
    let s_start = seed_s + 1 - sl;

    // Right extension (forward)
    let (score_r, qr, sr, right_ops) = if seed_q < query.len() - 1 && seed_s < subject.len() - 1 {
        align_ex(
            &query[seed_q..], &subject[seed_s..],
            query.len() - seed_q - 1, subject.len() - seed_s - 1,
            &matrix, gap_open, gap_extend, x_dropoff, false,
        )
    } else { (0, 0, 0, Vec::new()) };

    let total_score = score_l + score_r;
    if total_score <= 0 { return None; }

    // Build edit script: left ops (already in forward order) + right ops (reversed)
    let mut esp = GapEditScript::new();
    for &(op, cnt) in &left_ops { esp.push(op, cnt); }
    for &(op, cnt) in right_ops.iter().rev() { esp.push(op, cnt); }

    // Prune terminal gaps
    while !esp.ops.is_empty() && esp.ops[0].0 != GapAlignOpType::Sub { esp.ops.remove(0); }
    while !esp.ops.is_empty() && esp.ops.last().unwrap().0 != GapAlignOpType::Sub { esp.ops.pop(); }

    Some(TracebackResult {
        score: total_score,
        edit_script: esp,
        query_start: q_start,
        query_end: seed_q + qr + 1,
        subject_start: s_start,
        subject_end: seed_s + sr + 1,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_traceback_perfect() {
        let q = vec![0u8, 1, 2, 3, 0, 1, 2, 3];
        let s = vec![0u8, 1, 2, 3, 0, 1, 2, 3];
        let (score, esp, _, _, _, _) = traceback_align(&q, &s, 0, 8, 0, 8, 2, -3, 5, 2);
        assert_eq!(score, 16);
        assert_eq!(esp.ops.len(), 1);
        assert_eq!(esp.ops[0], (GapAlignOpType::Sub, 8));
    }

    #[test]
    fn test_traceback_with_mismatch() {
        let q = vec![0u8, 1, 2, 3, 0, 1, 2, 3];
        let s = vec![0u8, 1, 2, 3, 3, 1, 2, 3]; // pos 4: 0→3
        let (score, esp, _, _, _, _) = traceback_align(&q, &s, 0, 8, 0, 8, 2, -3, 5, 2);
        assert!(score > 0);
        let total: i32 = esp.ops.iter().map(|(_, n)| *n).sum();
        assert!(total >= 7, "alignment should be at least 7 bases");
    }

    #[test]
    fn test_traceback_abs_coordinates() {
        // Subject has ACGTACGT at position 100
        let mut subject = vec![3u8; 200]; // all T
        for (i, &b) in [0u8, 1, 2, 3, 0, 1, 2, 3].iter().enumerate() {
            subject[100 + i] = b;
        }
        let query = vec![0u8, 1, 2, 3, 0, 1, 2, 3];

        let result = traceback_align_abs(&query, &subject, 4, 104, 2, -3, 5, 2, 20);
        assert!(result.is_some());
        let r = result.unwrap();
        assert!(r.score > 0);
        assert!(r.subject_start >= 100 && r.subject_start <= 104,
            "Subject start should be near 100, got {}", r.subject_start);
        assert!(r.subject_end >= 104 && r.subject_end <= 108,
            "Subject end should be near 108, got {}", r.subject_end);
    }

    #[test]
    fn test_blast_gapped_align_basic() {
        let q = vec![0u8, 1, 2, 3, 0, 1, 2, 3]; // ACGTACGT
        let s = vec![0u8, 1, 2, 3, 0, 1, 2, 3];
        let r = blast_gapped_align(&q, &s, 4, 4, 1, -3, 5, 2, 30);
        assert!(r.is_some(), "Should find alignment");
        let r = r.unwrap();
        eprintln!("blast_gapped: score={} q={}..{} s={}..{} ops={:?}", r.score, r.query_start, r.query_end, r.subject_start, r.subject_end, r.edit_script.ops);
        assert!(r.score > 0, "score={}", r.score);
        // Check edit script has content
        let total_ops: i32 = r.edit_script.ops.iter().map(|(_,n)| *n).sum();
        assert!(total_ops > 0, "edit script should have operations, got {:?}", r.edit_script.ops);
    }

    #[test]
    fn test_align_ex_forward() {
        let mut matrix = [[-3i32; 16]; 16];
        for i in 0..4 { matrix[i][i] = 1; }
        let a = vec![255u8, 0, 1, 2, 3]; // pad + ACGT
        let b = vec![255u8, 0, 1, 2, 3, 255]; // pad + ACGT + pad
        let (score, ao, bo, ops) = align_ex(&a, &b, 4, 4, &matrix, 5, 2, 30, false);
        assert_eq!(score, 4, "4 matches * reward 1");
        assert_eq!(ao, 4);
        assert_eq!(bo, 4);
        assert_eq!(ops, vec![(GapAlignOpType::Sub, 4)]);
    }

    #[test]
    fn test_align_ex_reverse() {
        let mut matrix = [[-3i32; 16]; 16];
        for i in 0..4 { matrix[i][i] = 2; } // reward=2 for stronger diagonal signal
        let a = vec![0u8, 1, 2, 3, 0];
        let b = vec![0u8, 1, 2, 3, 0];
        let (score, ao, _bo, ops) = align_ex(&a, &b, 5, 5, &matrix, 5, 2, 30, true);
        eprintln!("rev2: score={} ao={} ops={:?}", score, ao, ops);
        assert!(score >= 8, "should get at least 4 bases * 2, got {}", score);
    }
}
