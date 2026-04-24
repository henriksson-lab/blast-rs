//! Gapped alignment with traceback.
//!
//! Mirrors NCBI's `blast_gapalign.c` (the ALIGN_EX / gapped-alignment half) and
//! `blast_traceback.c` (the HSP-list traceback half). Function-name mapping
//! maintained for side-by-side auditing:
//!
//! | Rust                             | NCBI source (blast_gapalign.c)           |
//! | -------------------------------- | ---------------------------------------- |
//! | `align_ex`                       | `ALIGN_EX` (line 374)                    |
//! | `blast_gapped_align`             | `BLAST_GappedAlignmentWithTraceback`     |
//! | `blast_gapped_score_only`        | `Blast_SemiGappedAlign`                  |
//! | `blast_gapped_score_only_packed_subject` | `s_BlastDynProgNtGappedAlignment` / `s_BlastAlignPackedNucl` |
//! | `gapped_score_one_dir`           | inner recurrence of `Blast_SemiGappedAlign` |
//! | `build_blastna_matrix`           | `BLAST_ScoreBlk` matrix fill (blast_stat.c) |
//! | `traceback_align` / `_abs`       | standalone NW helper, no C analog        |
//!
//! Script/op constants match NCBI exactly (see `SCRIPT_SUB` etc.) —
//! `blast_gapalign.c:363-371` and `gapinfo.h:45-51`.

use crate::gapinfo::{GapAlignOpType, GapEditScript};
use std::time::Instant;

// Tracehash emission for `align_ex` — dev-only, used for parity debugging
// against NCBI's `ALIGN_EX` in `blast_gapalign.c`. See
// `tests/parity/align_ex.md` for the C side's matching emission schema.
#[cfg(test)]
pub(crate) static ALIGN_EX_TRACEHASH_ENABLED: std::sync::atomic::AtomicBool =
    std::sync::atomic::AtomicBool::new(false);

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
            let match_score = if q[i - 1] == s[j - 1] {
                reward
            } else {
                penalty
            };
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

    (
        best_score,
        esp,
        align_q_start,
        align_q_end,
        align_s_start,
        align_s_end,
    )
}

/// Result of a traceback alignment.
pub struct TracebackResult {
    pub score: i32,
    pub edit_script: GapEditScript,
    pub query_start: usize, // 0-based in original sequence
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
        query, subject, q_start, q_end, s_start, s_end, reward, penalty, gap_open, gap_extend,
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

use crate::stat::MININT;
const SCRIPT_SUB: u8 = 3;
const SCRIPT_GAP_IN_A: u8 = 0; // gap in query = deletion
const SCRIPT_GAP_IN_B: u8 = 6; // gap in subject = insertion
const SCRIPT_OP_MASK: u8 = 0x07;
const SCRIPT_EXTEND_GAP_A: u8 = 0x10;
const SCRIPT_EXTEND_GAP_B: u8 = 0x40;

#[derive(Clone, Copy)]
struct GapDP {
    best: i32,
    best_gap: i32,
}

fn script_to_op(s: u8) -> GapAlignOpType {
    match s {
        SCRIPT_GAP_IN_A => GapAlignOpType::Del,
        SCRIPT_GAP_IN_B => GapAlignOpType::Ins,
        _ => GapAlignOpType::Sub,
    }
}

/// One-directional gapped extension with X-dropoff and traceback.
/// Port of `ALIGN_EX` in `ncbi-blast-2.17.0+-src/c++/src/algo/blast/core/blast_gapalign.c:374`.
fn align_ex(
    a: &[u8],
    b: &[u8],
    m: usize,
    n: usize,
    matrix: &[[i32; 16]; 16],
    gap_open: i32,
    gap_extend: i32,
    mut x_dropoff: i32,
    reverse: bool,
) -> (i32, usize, usize, Vec<(GapAlignOpType, i32)>) {
    let profile_enabled = std::env::var_os("BLAST_RS_PROFILE").is_some();
    let align_start = if profile_enabled {
        Some(Instant::now())
    } else {
        None
    };
    let gap_oe = gap_open + gap_extend;
    if x_dropoff < gap_oe {
        x_dropoff = gap_oe;
    }
    if m == 0 || n == 0 {
        return (0, 0, 0, Vec::new());
    }

    let extra = if gap_extend > 0 {
        (x_dropoff / gap_extend) as usize + 3
    } else {
        n + 3
    };
    // Cap allocation: band never exceeds the X-drop-driven width. The
    // score-only path uses the same principle (`gapped_score_one_dir`), and
    // letting traceback grow unbounded here makes long near-self hits
    // pathological even when the intended X-drop band is narrow.
    let max_band = (extra * 2 + 20).min(n + extra + 10);
    let mut sa = vec![
        GapDP {
            best: MININT,
            best_gap: MININT
        };
        max_band
    ];

    // Row 0 initialization
    sa[0] = GapDP {
        best: 0,
        best_gap: -gap_oe,
    };
    let mut score = -gap_oe;
    let mut b_size = 1usize;
    while b_size <= n && score >= -x_dropoff {
        sa[b_size] = GapDP {
            best: score,
            best_gap: score - gap_oe,
        };
        score -= gap_extend;
        b_size += 1;
    }

    // Traceback storage: compact per-row scripts keyed by `row_starts`, mirroring
    // NCBI's `edit_script[a_index][b_index - edit_start_offset[a_index]]`.
    let mut scripts: Vec<Vec<u8>> = Vec::with_capacity(m + 1);
    let mut row_starts: Vec<usize> = Vec::with_capacity(m + 1);
    scripts.push(vec![SCRIPT_GAP_IN_A; b_size]);
    row_starts.push(0);

    let mut best_score = 0i32;
    let mut first_b = 0usize;
    let mut a_off = 0usize;
    let mut b_off = 0usize;
    let mut row_count = 0usize;
    let mut max_b_size = b_size;
    let mut max_row_script_len = scripts[0].len();

    for ai in 1..=m {
        row_count += 1;
        let a_letter = if reverse { a[m - ai] } else { a[ai] };
        let mrow = &matrix[a_letter as usize & 0x0F];

        let row_start = first_b;
        let mut row_script = vec![0u8; (b_size - row_start) + extra + 10];
        max_b_size = max_b_size.max(b_size);
        max_row_script_len = max_row_script_len.max(row_script.len());
        let mut sc = MININT;
        let mut sgr = MININT; // score_gap_row
        let mut last_b = first_b;

        // Mutating `first_b` inside the loop advances the DP band for the
        // NEXT outer iteration (NCBI `blast_gapalign.c:601`); the current
        // iterator was snapshotted at loop entry, so the mutation is safe.
        #[allow(clippy::mut_range_bound)]
        for bi in first_b..b_size {
            // NCBI blast_gapalign.c:566 — `b_ptr += b_increment` happens before
            // the letter is read. Iter `bi` reads B[bi + 1] (forward) or
            // B[N - bi - 1] (reverse). That letter is used only to compute
            // `next_score` for the NEXT cell — it does NOT affect the current
            // cell's decision. So when `b_idx` would be out of bounds we must
            // still process the current cell (sc is already set from the
            // previous iteration's next_score) and just skip setting next_sc.
            let b_idx = if reverse {
                n.checked_sub(1 + bi).unwrap_or(usize::MAX)
            } else {
                bi + 1
            };
            let b_in_range = b_idx < b.len();
            let sgc = sa[bi].best_gap; // score_gap_col
            let next_sc = if b_in_range {
                let b_letter = b[b_idx];
                sa[bi].best + mrow[b_letter as usize & 0x0F]
            } else {
                MININT
            };

            // Best predecessor (NCBI blast_gapalign.c:588-599).
            let mut script = SCRIPT_SUB;
            if sc < sgc {
                script = SCRIPT_GAP_IN_B;
                sc = sgc;
            }
            if sc < sgr {
                script = SCRIPT_GAP_IN_A;
                sc = sgr;
            }

            if best_score - sc > x_dropoff {
                // Match NCBI ALIGN_EX (blast_gapalign.c:601): when pruning the
                // leading column (`first_b_index == b_index`), advance
                // `first_b_index` and DO NOT overwrite `score_array[b_index]`.
                // The stale prev-row value remains so that next row's SUB
                // predecessor (at column new_first_b - 1) reads a meaningful
                // score, not MININT. Only non-leading pruned cells are set to
                // MININT.
                if first_b == bi {
                    first_b += 1;
                } else {
                    sa[bi].best = MININT;
                }
            } else {
                last_b = bi;
                if sc > best_score {
                    best_score = sc;
                    a_off = ai;
                    b_off = bi;
                }

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
            let script_idx = bi.saturating_sub(row_start);
            if script_idx < row_script.len() {
                row_script[script_idx] = script;
            }
        }

        // Dev-only tracehash emission: one event per DP row. The C-side
        // `ALIGN_EX` should emit the same-named event with the same fields
        // in the same order for a row-by-row diff.
        #[cfg(test)]
        if ALIGN_EX_TRACEHASH_ENABLED.load(std::sync::atomic::Ordering::Relaxed) {
            let mut th = tracehash::th_call!("align_ex_row");
            th.input_u64(ai as u64);
            th.input_u64(first_b as u64);
            th.input_u64(last_b as u64);
            th.input_u64(b_size as u64);
            th.input_i64(best_score as i64);
            th.input_i64(a_off as i64);
            th.input_i64(b_off as i64);
            // score_array slice: pack best/best_gap as interleaved i32 little-endian.
            let mut buf = Vec::with_capacity(b_size * 8);
            for cell in &sa[..b_size.min(sa.len())] {
                buf.extend_from_slice(&cell.best.to_le_bytes());
                buf.extend_from_slice(&cell.best_gap.to_le_bytes());
            }
            th.input_bytes(&buf);
            // Only the [first_b..b_size) range of the row script is written
            // inside the inner loop; the prefix [0..first_b) is uninitialized
            // in NCBI (pointer arithmetic trick on `edit_script_row`). Hash
            // the same range on both sides so the two traces are comparable.
            let used = b_size.saturating_sub(row_start).min(row_script.len());
            th.input_bytes(&row_script[..used]);
            th.finish();
        }

        if first_b >= b_size {
            scripts.push(row_script);
            row_starts.push(row_start);
            break;
        }

        if last_b < b_size - 1 {
            b_size = last_b + 1;
        } else {
            // Extend band rightward
            while sgr >= best_score - x_dropoff && b_size <= n && b_size < sa.len() {
                sa[b_size] = GapDP {
                    best: sgr,
                    best_gap: sgr - gap_oe,
                };
                sgr -= gap_extend;
                let script_idx = b_size.saturating_sub(row_start);
                if script_idx < row_script.len() {
                    row_script[script_idx] = SCRIPT_GAP_IN_A;
                }
                b_size += 1;
            }
        }
        if b_size <= n && b_size < sa.len() {
            sa[b_size] = GapDP {
                best: MININT,
                best_gap: MININT,
            };
            b_size += 1;
        }
        scripts.push(row_script);
        row_starts.push(row_start);
    }

    // Traceback
    let mut ops: Vec<(GapAlignOpType, i32)> = Vec::new();
    let mut ai = a_off;
    let mut bi = b_off;
    let mut cur_script = SCRIPT_SUB;
    let mut traceback_steps = 0usize;
    #[cfg(test)]
    let mut tb_step: u64 = 0;

    while ai > 0 || bi > 0 {
        traceback_steps += 1;
        if ai >= scripts.len() {
            break;
        }
        let row_start = row_starts[ai];
        if bi < row_start {
            break;
        }
        let script_idx = bi - row_start;
        if script_idx >= scripts[ai].len() {
            break;
        }
        let s = scripts[ai][script_idx];

        cur_script = match cur_script & SCRIPT_OP_MASK {
            SCRIPT_GAP_IN_A => {
                if s & SCRIPT_EXTEND_GAP_A != 0 {
                    SCRIPT_GAP_IN_A
                } else {
                    s & SCRIPT_OP_MASK
                }
            }
            SCRIPT_GAP_IN_B => {
                if s & SCRIPT_EXTEND_GAP_B != 0 {
                    SCRIPT_GAP_IN_B
                } else {
                    s & SCRIPT_OP_MASK
                }
            }
            _ => s & SCRIPT_OP_MASK,
        };

        #[cfg(test)]
        if ALIGN_EX_TRACEHASH_ENABLED.load(std::sync::atomic::Ordering::Relaxed) {
            let mut th = tracehash::th_call!("align_ex_tb");
            th.input_u64(tb_step);
            th.input_u64(ai as u64);
            th.input_u64(bi as u64);
            th.input_u64(cur_script as u64);
            th.finish();
            tb_step += 1;
        }

        let op = script_to_op(cur_script);
        match cur_script & SCRIPT_OP_MASK {
            SCRIPT_GAP_IN_A => {
                if bi == 0 {
                    break;
                }
                bi -= 1;
            }
            SCRIPT_GAP_IN_B => {
                if ai == 0 {
                    break;
                }
                ai -= 1;
            }
            _ => {
                if ai == 0 || bi == 0 {
                    break;
                }
                ai -= 1;
                bi -= 1;
            }
        }
        if let Some(last) = ops.last_mut() {
            if last.0 == op {
                last.1 += 1;
            } else {
                ops.push((op, 1));
            }
        } else {
            ops.push((op, 1));
        }
    }
    if let Some(start) = align_start {
        eprintln!(
            "[blastn-profile] align_ex reverse={} m={} n={} xdrop={} rows={} max_b_size={} max_row_script_len={} traceback_steps={} best_score={} total_ms={}",
            reverse,
            m,
            n,
            x_dropoff,
            row_count,
            max_b_size,
            max_row_script_len,
            traceback_steps,
            best_score,
            start.elapsed().as_millis()
        );
    }

    (best_score, a_off, b_off, ops)
}

/// Build full BLASTNA scoring matrix — verbatim port of NCBI
/// `BlastScoreBlkNuclMatrixCreate` (`blast_stat.c:1060`).
pub fn build_blastna_matrix(reward: i32, penalty: i32) -> [[i32; 16]; 16] {
    use crate::encoding::{BLASTNA_SIZE, BLASTNA_TO_NCBI4NA};
    const K_NUMBER_NON_AMBIG_BP: usize = 4;

    let mut matrix = [[0i32; 16]; BLASTNA_SIZE];

    // `degeneracy[i]` counts how many unambiguous bases (A,C,G,T) the
    // residue `i` can represent. NCBI `blast_stat.c:1087-1099`.
    let mut degeneracy = [0i32; BLASTNA_SIZE];
    for index1 in 0..K_NUMBER_NON_AMBIG_BP {
        degeneracy[index1] = 1;
    }
    for index1 in K_NUMBER_NON_AMBIG_BP..BLASTNA_SIZE {
        let mut degen = 0;
        for index2 in 0..K_NUMBER_NON_AMBIG_BP {
            if BLASTNA_TO_NCBI4NA[index1] & BLASTNA_TO_NCBI4NA[index2] != 0 {
                degen += 1;
            }
        }
        degeneracy[index1] = degen;
    }

    for index1 in 0..BLASTNA_SIZE {
        for index2 in index1..BLASTNA_SIZE {
            let s = if BLASTNA_TO_NCBI4NA[index1] & BLASTNA_TO_NCBI4NA[index2] != 0 {
                // NCBI `blast_stat.c:1107`:
                //   (Int4)BLAST_Nint((double)((degeneracy[index2]-1)*penalty + reward)
                //                    / (double)degeneracy[index2])
                crate::math::nint(
                    ((degeneracy[index2] - 1) * penalty + reward) as f64
                        / degeneracy[index2] as f64,
                ) as i32
            } else {
                penalty
            };
            matrix[index1][index2] = s;
            if index1 != index2 {
                matrix[index2][index1] = s;
            }
        }
    }

    // Gap-sentinel row/col — NCBI `blast_stat.c:1124-1127`.
    for index1 in 0..BLASTNA_SIZE {
        matrix[BLASTNA_SIZE - 1][index1] = MININT;
        matrix[index1][BLASTNA_SIZE - 1] = MININT;
    }
    matrix
}

/// BLAST-style gapped alignment with traceback — extends bidirectionally from seed.
/// Fast score-only gapped extension (no traceback).
/// Port of C engine's Blast_SemiGappedAlign — computes only the score.
pub fn blast_gapped_score_only(
    query: &[u8],
    subject: &[u8],
    seed_q: usize,
    seed_s: usize,
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    x_dropoff: i32,
) -> i32 {
    let gap_oe = gap_open + gap_extend;
    let matrix = build_blastna_matrix(reward, penalty);

    // Left extension (score only)
    let score_l = gapped_score_one_dir(
        &query[..seed_q + 1],
        &subject[..seed_s + 1],
        seed_q + 1,
        seed_s + 1,
        &matrix,
        gap_oe,
        gap_extend,
        x_dropoff,
        true,
    );

    // Right extension (score only)
    let score_r = if seed_q < query.len() - 1 && seed_s < subject.len() - 1 {
        gapped_score_one_dir(
            &query[seed_q..],
            &subject[seed_s..],
            query.len() - seed_q - 1,
            subject.len() - seed_s - 1,
            &matrix,
            gap_oe,
            gap_extend,
            x_dropoff,
            false,
        )
    } else {
        0
    };

    score_l + score_r
}

#[inline]
fn unpack_ncbi2na_base(packed: &[u8], pos: usize) -> u8 {
    let byte = packed[pos >> 2];
    (byte >> (6 - 2 * (pos & 3))) & 0x03
}

/// Score-only port of BLASTN's packed-subject preliminary DP path
/// (`s_BlastDynProgNtGappedAlignment` + `s_BlastAlignPackedNucl`).
pub fn blast_gapped_score_only_packed_subject(
    query: &[u8],
    subject_packed: &[u8],
    subject_len: usize,
    seed_q: usize,
    seed_s: usize,
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    x_dropoff: i32,
    ungapped_score: i32,
) -> i32 {
    blast_gapped_score_extents_packed_subject(
        query,
        subject_packed,
        subject_len,
        seed_q,
        seed_s,
        reward,
        penalty,
        gap_open,
        gap_extend,
        x_dropoff,
        ungapped_score,
    )
    .0
}

/// Packed-subject preliminary DP that returns the same information NCBI's
/// `s_BlastDynProgNtGappedAlignment` derives before saving a prelim HSP:
/// score plus preliminary extents around the chosen start.
pub fn blast_gapped_score_extents_packed_subject(
    query: &[u8],
    subject_packed: &[u8],
    subject_len: usize,
    seed_q: usize,
    seed_s: usize,
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    x_dropoff: i32,
    ungapped_score: i32,
) -> (i32, usize, usize, usize, usize, usize, usize, usize, usize, usize, usize, i32, i32) {
    let mut effective_x_dropoff = x_dropoff;
    if ungapped_score < effective_x_dropoff {
        effective_x_dropoff = ungapped_score;
    }

    let offset_adjustment = 4usize - (seed_s % 4);
    let mut q_length = seed_q + offset_adjustment;
    let mut s_length = seed_s + offset_adjustment;
    if q_length > query.len() || s_length > subject_len {
        q_length = q_length.saturating_sub(4);
        s_length = s_length.saturating_sub(4);
    }

    let matrix = build_blastna_matrix(reward, penalty);
    let gap_oe = gap_open + gap_extend;
    let (score_left, private_q_start, private_s_start) = gapped_score_one_dir_packed_subject(
        query,
        subject_packed,
        q_length,
        s_length,
        0,
        &matrix,
        gap_oe,
        gap_extend,
        effective_x_dropoff,
        true,
    );
    let query_start = q_length.saturating_sub(private_q_start);
    let subject_start = s_length.saturating_sub(private_s_start);

    let (score_right, query_stop, subject_stop) = if q_length < query.len() && s_length < subject_len {
        gapped_score_one_dir_packed_subject(
            &query[q_length - 1..],
            subject_packed,
            query.len() - q_length,
            subject_len - s_length,
            s_length,
            &matrix,
            gap_oe,
            gap_extend,
            effective_x_dropoff,
            false,
        )
    } else {
        (0, 0, 0)
    };

    let query_stop = if q_length < query.len() && s_length < subject_len {
        query_stop + q_length
    } else {
        q_length
    };
    let subject_stop = if q_length < query.len() && s_length < subject_len {
        subject_stop + s_length
    } else {
        s_length
    };

    (
        score_left + score_right,
        query_start,
        query_stop,
        subject_start,
        subject_stop,
        seed_q,
        seed_s,
        q_length,
        s_length,
        private_q_start,
        private_s_start,
        score_left,
        score_right,
    )
}

/// One-directional score-only gapped extension with X-dropoff.
/// Much faster than align_ex because no traceback storage/recording.
fn gapped_score_one_dir(
    a: &[u8],
    b: &[u8],
    m: usize,
    n: usize,
    matrix: &[[i32; 16]; 16],
    gap_oe: i32,
    gap_extend: i32,
    mut x_dropoff: i32,
    reverse: bool,
) -> i32 {
    if x_dropoff < gap_oe {
        x_dropoff = gap_oe;
    }
    if m == 0 || n == 0 {
        return 0;
    }

    let num_extra_cells = if gap_extend > 0 {
        (x_dropoff / gap_extend) as usize + 3
    } else {
        n + 3
    };
    let mut sa = vec![
        GapDP {
            best: MININT,
            best_gap: MININT
        };
        n + num_extra_cells + 10
    ];

    sa[0] = GapDP {
        best: 0,
        best_gap: -gap_oe,
    };
    let mut score = -gap_oe;
    let mut b_size = 1usize;
    while b_size <= n && score >= -x_dropoff {
        sa[b_size] = GapDP {
            best: score,
            best_gap: score - gap_oe,
        };
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

        // See the parallel forward loop for why `first_b` mutation is safe.
        #[allow(clippy::mut_range_bound)]
        for bi in first_b..b_size {
            let b_idx = if reverse {
                n.checked_sub(1 + bi).unwrap_or(usize::MAX)
            } else {
                bi + 1
            };
            if b_idx >= b.len() {
                break;
            }
            let b_letter = b[b_idx];

            let sgc = sa[bi].best_gap;
            let next_sc = sa[bi].best + mrow[b_letter as usize & 0x0F];

            // Three-way max
            if sc < sgc {
                sc = sgc;
            }
            if sc < sgr {
                sc = sgr;
            }

            if best_score - sc > x_dropoff {
                if first_b == bi {
                    first_b += 1;
                }
                sa[bi].best = MININT;
            } else {
                last_b = bi;
                if sc > best_score {
                    best_score = sc;
                }
                if sgc - gap_extend < sc - gap_oe {
                    sa[bi].best_gap = sc - gap_oe;
                } else {
                    sa[bi].best_gap = sgc - gap_extend;
                }
                if sgr - gap_extend < sc - gap_oe {
                    sgr = sc - gap_oe;
                } else {
                    sgr -= gap_extend;
                }
                sa[bi].best = sc;
            }
            sc = next_sc;
        }

        if first_b >= b_size {
            break;
        }

        if last_b + num_extra_cells + 3 >= sa.len() {
            sa.resize(
                (last_b + num_extra_cells + 100).max(sa.len() * 2),
                GapDP {
                    best: MININT,
                    best_gap: MININT,
                },
            );
        }

        if last_b < b_size - 1 {
            b_size = last_b + 1;
        } else {
            while sgr >= best_score - x_dropoff && b_size <= n {
                sa[b_size] = GapDP {
                    best: sgr,
                    best_gap: sgr - gap_oe,
                };
                sgr -= gap_extend;
                b_size += 1;
            }
        }
        if b_size <= n && b_size < sa.len() {
            sa[b_size] = GapDP {
                best: MININT,
                best_gap: MININT,
            };
            b_size += 1;
        }
    }
    best_score
}

fn gapped_score_one_dir_packed_subject(
    query: &[u8],
    subject_packed: &[u8],
    n: usize,
    m: usize,
    subject_base_offset: usize,
    matrix: &[[i32; 16]; 16],
    gap_oe: i32,
    gap_extend: i32,
    mut x_dropoff: i32,
    reverse: bool,
) -> (i32, usize, usize) {
    if x_dropoff < gap_oe {
        x_dropoff = gap_oe;
    }
    if m == 0 || n == 0 {
        return (0, 0, 0);
    }

    let num_extra_cells = if gap_extend > 0 {
        (x_dropoff / gap_extend) as usize + 3
    } else {
        n + 3
    };
    let mut sa = vec![
        GapDP {
            best: MININT,
            best_gap: MININT
        };
        n + num_extra_cells + 10
    ];

    sa[0] = GapDP {
        best: 0,
        best_gap: -gap_oe,
    };
    let mut score = -gap_oe;
    let mut b_size = 1usize;
    while b_size <= n && score >= -x_dropoff {
        sa[b_size] = GapDP {
            best: score,
            best_gap: score - gap_oe,
        };
        score -= gap_extend;
        b_size += 1;
    }

    let mut best_score = 0i32;
    let mut best_a = 0usize;
    let mut best_b = 0usize;
    let mut first_b = 0usize;

    for ai in 1..=m {
        let subj_base = if reverse {
            let abs_pos = subject_base_offset + (m - ai);
            let byte = subject_packed[abs_pos >> 2];
            (byte >> (2 * (3 - ((ai - 1) & 3)))) & 0x03
        } else {
            unpack_ncbi2na_base(subject_packed, subject_base_offset + ai - 1)
        };
        let mrow = &matrix[subj_base as usize & 0x0F];
        let mut sc = MININT;
        let mut sgr = MININT;
        let mut last_b = first_b;

        for bi in first_b..b_size {
            let q_idx = if reverse {
                n.checked_sub(1 + bi).unwrap_or(usize::MAX)
            } else {
                bi + 1
            };
            if q_idx >= query.len() {
                break;
            }
            let q_base = query[q_idx];
            let sgc = sa[bi].best_gap;
            let next_sc = sa[bi].best + mrow[q_base as usize & 0x0F];

            if sc < sgc {
                sc = sgc;
            }
            if sc < sgr {
                sc = sgr;
            }

            if best_score - sc > x_dropoff {
                if first_b == bi {
                    first_b += 1;
                } else {
                    sa[bi].best = MININT;
                }
            } else {
                last_b = bi;
                if sc > best_score {
                    best_score = sc;
                    best_a = ai;
                    best_b = bi;
                }
                sa[bi].best_gap = (sc - gap_oe).max(sgc - gap_extend);
                sgr = (sc - gap_oe).max(sgr - gap_extend);
                sa[bi].best = sc;
            }
            sc = next_sc;
        }

        if first_b >= b_size {
            break;
        }

        if last_b + num_extra_cells + 3 >= sa.len() {
            sa.resize(
                (last_b + num_extra_cells + 100).max(sa.len() * 2),
                GapDP {
                    best: MININT,
                    best_gap: MININT,
                },
            );
        }

        if last_b < b_size - 1 {
            b_size = last_b + 1;
        } else {
            while sgr >= best_score - x_dropoff && b_size <= n {
                sa[b_size] = GapDP {
                    best: sgr,
                    best_gap: sgr - gap_oe,
                };
                sgr -= gap_extend;
                b_size += 1;
            }
        }
        if b_size <= n && b_size < sa.len() {
            sa[b_size] = GapDP {
                best: MININT,
                best_gap: MININT,
            };
            b_size += 1;
        }
    }
    (best_score, best_b, best_a)
}

/// Port of NCBI `Blast_HSPReevaluateWithAmbiguitiesGapped`
/// (`blast_hits.c:479`). Walks the edit script, tracks the running sum,
/// trims below-cutoff prefix/suffix segments by moving `best_*_start/end`,
/// then extends the best-scoring subrange outward via exact-match walks.
///
/// Returns `true` if the HSP should be deleted (best score < cutoff or
/// empty alignment). On a successful refine, the input `TracebackResult`
/// is modified in place: endpoints updated, edit script pruned, score
/// recomputed.
///
/// Note: NCBI multiplies scores by `factor=2` only when the scoring is
/// non-affine greedy (`gap_open == 0 && gap_extend == 0 && reward % 2 == 1`).
/// For the affine blastn / blastn-short path the factor is 1 and the
/// implementation below only keeps the affine branch — mirroring
/// `blast_hits.c:517-526`.
#[allow(clippy::too_many_arguments)]
pub fn reevaluate_with_ambiguities_gapped(
    tb: &mut TracebackResult,
    query: &[u8],
    subject: &[u8],
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    cutoff_score: i32,
) -> bool {
    // Build the BLASTNA scoring matrix exactly like `align_ex` so that
    // substitution scoring matches what the DP used. NCBI does the same via
    // `sbp->matrix->data[*query & 0x0f][*subject]` (`blast_hits.c:548`).
    let matrix = build_blastna_matrix(reward, penalty);
    if tb.edit_script.ops.is_empty() {
        return true;
    }
    let (factor, effective_gap_open, effective_gap_extend) = if gap_open == 0 && gap_extend == 0 {
        // NCBI `Blast_HSPReevaluateWithAmbiguitiesGapped`: non-affine greedy
        // stores gap_open/gap_extend as 0, then reevaluates with factor=2 and
        // a synthetic linear gap cost.
        (2, 0, (reward - 2 * penalty))
    } else {
        (1, gap_open, gap_extend)
    };

    // Walk pointers (absolute coords into original query/subject).
    let mut qp = tb.query_start;
    let mut sp = tb.subject_start;

    // Track best-scoring subpath by position and by ESP index.
    let mut sum: i32 = 0;
    let mut score: i32 = 0;
    let mut best_q_start = qp;
    let mut best_s_start = sp;
    let mut best_q_end = qp;
    let mut best_s_end = sp;
    let mut current_q_start = qp;
    let mut current_s_start = sp;

    // ESP indices; we modify `ops` in place for the "split current esp on
    // sum<0 in the middle of a run" case.
    let mut best_start_esp_index: usize = 0;
    let mut best_end_esp_index: usize = 0;
    let mut current_start_esp_index: usize = 0;
    let mut best_end_esp_num: i32 = 0;

    // Clone the ops into a mutable buffer so we can split a run if needed.
    let mut ops = tb.edit_script.ops.clone();

    let mut index = 0;
    while index < ops.len() {
        let op_type = ops[index].0;
        let num_at_entry = ops[index].1;
        let mut op_index: i32 = 0;
        while op_index < ops[index].1 {
            match op_type {
                GapAlignOpType::Sub => {
                    // Apply one substitution using the BLASTNA matrix.
                    let q = (query.get(qp).copied().unwrap_or(15) & 0x0f) as usize;
                    let s = (subject.get(sp).copied().unwrap_or(15) & 0x0f) as usize;
                    sum += factor * matrix[q][s];
                    qp += 1;
                    sp += 1;
                    op_index += 1;
                }
                GapAlignOpType::Del | GapAlignOpType::Del1 | GapAlignOpType::Del2 => {
                    // Gap in query; advance subject by full run.
                    sum -= effective_gap_open + effective_gap_extend * num_at_entry;
                    sp += num_at_entry as usize;
                    op_index += num_at_entry;
                }
                GapAlignOpType::Ins | GapAlignOpType::Ins1 | GapAlignOpType::Ins2 => {
                    // Gap in subject; advance query by full run.
                    sum -= effective_gap_open + effective_gap_extend * num_at_entry;
                    qp += num_at_entry as usize;
                    op_index += num_at_entry;
                }
                GapAlignOpType::Decline => {
                    op_index += num_at_entry;
                }
            }

            if sum < 0 {
                // If we're mid-run (only possible for a Sub run), split
                // this esp entry: reduce its count, fall back to op_index=0
                // so the next iteration re-enters at the split point.
                if op_index < ops[index].1 {
                    ops[index].1 -= op_index;
                    current_start_esp_index = index;
                    op_index = 0;
                } else {
                    current_start_esp_index = index + 1;
                }
                sum = 0;
                current_q_start = qp;
                current_s_start = sp;
                // If cumulative score hasn't yet cleared the cutoff, drop
                // everything seen so far.
                if score < cutoff_score {
                    best_q_start = qp;
                    best_s_start = sp;
                    score = 0;
                    best_start_esp_index = current_start_esp_index;
                    best_end_esp_index = current_start_esp_index;
                }
            } else if sum > score {
                score = sum;
                best_q_start = current_q_start;
                best_s_start = current_s_start;
                best_q_end = qp;
                best_s_end = sp;
                best_start_esp_index = current_start_esp_index;
                best_end_esp_index = index;
                best_end_esp_num = op_index;
            }
        }
        index += 1;
    }

    score /= factor;

    if score < cutoff_score {
        return true;
    }

    // Post-processing: extend left and right from best_q_start/best_q_end
    // while exact matches continue. Only if both endpoints fall on Sub
    // runs (NCBI asserts this — blast_hits.c:616-617).
    let within = best_start_esp_index < ops.len() && best_end_esp_index < ops.len();
    if within
        && ops[best_start_esp_index].0 == GapAlignOpType::Sub
        && ops[best_end_esp_index].0 == GapAlignOpType::Sub
    {
        // Extend left.
        let mut ext: i32 = 0;
        let mut qpp = best_q_start;
        let mut spp = best_s_start;
        while qpp > 0 && spp > 0 {
            let q = query.get(qpp - 1).copied().unwrap_or(15) & 0x0f;
            let s = subject.get(spp - 1).copied().unwrap_or(15);
            // NCBI also requires q[qp]<4 (unambiguous) — blast_hits.c:622.
            if q != s || q >= 4 {
                break;
            }
            qpp -= 1;
            spp -= 1;
            ext += 1;
        }
        best_q_start = qpp;
        best_s_start = spp;
        ops[best_start_esp_index].1 += ext;
        if best_end_esp_index == best_start_esp_index {
            best_end_esp_num += ext;
        }
        score += ext * reward;

        // Extend right.
        let mut ext: i32 = 0;
        let mut qpp = best_q_end;
        let mut spp = best_s_end;
        while qpp < query.len() && spp < subject.len() {
            let q = query[qpp] & 0x0f;
            let s = subject[spp];
            if q != s || q >= 4 {
                break;
            }
            qpp += 1;
            spp += 1;
            ext += 1;
        }
        best_q_end = qpp;
        best_s_end = spp;
        ops[best_end_esp_index].1 += ext;
        best_end_esp_num += ext;
        score += ext * reward;
    }

    // Trim ops to [best_start_esp_index ..= best_end_esp_index] with
    // best_end_esp_num controlling the last run's count.
    if best_end_esp_index >= ops.len() || best_start_esp_index > best_end_esp_index {
        return true;
    }
    let mut trimmed: Vec<(GapAlignOpType, i32)> = Vec::new();
    for i in best_start_esp_index..=best_end_esp_index {
        let (op, n) = ops[i];
        if i == best_end_esp_index {
            let effective = if best_end_esp_num > 0 {
                best_end_esp_num
            } else {
                n
            };
            if effective > 0 {
                trimmed.push((op, effective));
            }
        } else {
            trimmed.push((op, n));
        }
    }

    tb.score = score;
    tb.edit_script.ops = trimmed;
    tb.query_start = best_q_start;
    tb.subject_start = best_s_start;
    tb.query_end = best_q_end;
    tb.subject_end = best_s_end;
    false
}

pub fn blast_gapped_align(
    query: &[u8],
    subject: &[u8],
    seed_q: usize,
    seed_s: usize,
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    x_dropoff: i32,
) -> Option<TracebackResult> {
    // Build full BLASTNA scoring matrix matching C engine (handles ambiguous bases)
    let matrix = build_blastna_matrix(reward, penalty);

    // Left extension (reverse)
    let (score_l, ql, sl, left_ops) = align_ex(
        &query[..seed_q + 1],
        &subject[..seed_s + 1],
        seed_q + 1,
        seed_s + 1,
        &matrix,
        gap_open,
        gap_extend,
        x_dropoff,
        true,
    );
    let q_start = seed_q + 1 - ql;
    let s_start = seed_s + 1 - sl;

    // Right extension (forward)
    let (score_r, qr, sr, right_ops) = if seed_q < query.len() - 1 && seed_s < subject.len() - 1 {
        align_ex(
            &query[seed_q..],
            &subject[seed_s..],
            query.len() - seed_q - 1,
            subject.len() - seed_s - 1,
            &matrix,
            gap_open,
            gap_extend,
            x_dropoff,
            false,
        )
    } else {
        (0, 0, 0, Vec::new())
    };

    let total_score = score_l + score_r;
    if total_score <= 0 {
        return None;
    }

    // Build edit script: left ops (already in forward order) + right ops (reversed)
    let mut esp = GapEditScript::new();
    for &(op, cnt) in &left_ops {
        esp.push(op, cnt);
    }
    for &(op, cnt) in right_ops.iter().rev() {
        esp.push(op, cnt);
    }

    // Prune terminal gaps
    while !esp.ops.is_empty() && esp.ops[0].0 != GapAlignOpType::Sub {
        esp.ops.remove(0);
    }
    while !esp.ops.is_empty() && esp.ops.last().unwrap().0 != GapAlignOpType::Sub {
        esp.ops.pop();
    }

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
    #[ignore = "parameter-sweep probe: find (seed, x_dropoff) that reproduces NCBI's BTOP"]
    fn probe_align_ex_adjacent_del_ins_params() {
        fn encode(s: &[u8]) -> Vec<u8> {
            s.iter()
                .map(|&c| match c {
                    b'A' => 0,
                    b'C' => 1,
                    b'G' => 2,
                    b'T' => 3,
                    _ => 15,
                })
                .collect()
        }
        let q = encode(
            b"ACGTTGCAACGATCGTACGATTCGAGCTTAGGCTAGGGTAATCGGATCCTAGCTAGGCTAATCGATCGTAGCTAGCATCGAT",
        );
        let s = encode(
            b"ACGTTGCAACGATCGTACGATTCGAGCTTAGGCTAATCGGATCCTAGCTAGGCTAATCGATCGTAGCTAGCATCGAT",
        );
        let seeds = [0, 10, 17, 20, 27, 30, 34];
        let xdrops = [10, 12, 14, 16, 18, 20, 22, 25, 30, 35, 40, 50];
        // NCBI reference: [(Sub, 35), (Ins, 5), (Sub, 42)]
        for &sp in &seeds {
            for &x in &xdrops {
                if let Some(rr) = blast_gapped_align(&q, &s, sp, sp, 1, -3, 5, 2, x) {
                    let ops: Vec<_> = rr.edit_script.ops.clone();
                    // Collapse consecutive same-op runs to match BTOP shape.
                    let mut collapsed: Vec<(GapAlignOpType, i32)> = Vec::new();
                    for (op, n) in ops {
                        if let Some(last) = collapsed.last_mut() {
                            if last.0 == op {
                                last.1 += n;
                                continue;
                            }
                        }
                        collapsed.push((op, n));
                    }
                    let is_ncbi = collapsed.len() == 3
                        && collapsed[0] == (GapAlignOpType::Sub, 35)
                        && collapsed[1] == (GapAlignOpType::Ins, 5)
                        && collapsed[2] == (GapAlignOpType::Sub, 42);
                    eprintln!(
                        "seed={:2} xdrop={:3} score={} {} ops={:?}",
                        sp,
                        x,
                        rr.score,
                        if is_ncbi { "== NCBI" } else { "" },
                        collapsed,
                    );
                }
            }
        }
    }

    #[test]
    #[ignore = "parity probe: emits a tracehash-rs trace for comparison against NCBI ALIGN_EX"]
    fn trace_align_ex_adjacent_del_ins() {
        // This test dumps the DP state of each row from blast_gapped_align
        // into a tracehash-rs trace file so it can be diffed row-by-row
        // against an NCBI C build that emits the same `align_ex_row` events.
        // Run with:
        //   TRACEHASH_OUT=/tmp/rust.tsv TRACEHASH_SIDE=rust \
        //   TRACEHASH_RUN_ID=adj_del_ins \
        //   cargo test --release --lib trace_align_ex_adjacent_del_ins \
        //     -- --ignored --nocapture
        fn encode(s: &[u8]) -> Vec<u8> {
            s.iter()
                .map(|&c| match c {
                    b'A' => 0,
                    b'C' => 1,
                    b'G' => 2,
                    b'T' => 3,
                    _ => 15,
                })
                .collect()
        }
        let q = encode(
            b"ACGTTGCAACGATCGTACGATTCGAGCTTAGGCTAGGGTAATCGGATCCTAGCTAGGCTAATCGATCGTAGCTAGCATCGAT",
        );
        let s = encode(
            b"ACGTTGCAACGATCGTACGATTCGAGCTTAGGCTAATCGGATCCTAGCTAGGCTAATCGATCGTAGCTAGCATCGAT",
        );
        ALIGN_EX_TRACEHASH_ENABLED.store(true, std::sync::atomic::Ordering::Relaxed);
        let r =
            blast_gapped_align(&q, &s, 17, 17, 1, -3, 5, 2, 16).expect("alignment should succeed");
        ALIGN_EX_TRACEHASH_ENABLED.store(false, std::sync::atomic::Ordering::Relaxed);
        eprintln!(
            "alignment: score={} q=[{}..{}] s=[{}..{}] ops={:?}",
            r.score, r.query_start, r.query_end, r.subject_start, r.subject_end, r.edit_script.ops,
        );
    }

    #[test]
    #[ignore = "low-level reproducer: output depends on which ungapped HSP's seed is picked. \
                The end-to-end parity test `blastn_subject_ncbi_parity_gapped_traceback_edge_matrix` \
                covers this fixture via the CLI path and is the authoritative regression guard."]
    fn test_adjacent_del_ins_gap_position_matches_ncbi() {
        // Historic reproducer (adjacent_del_ins). NCBI places the 5-gap block
        // after position 35 of the query (BTOP "35G-G-G-T-A-42"). Fixed on
        // 2026-04-18 by routing `x_dropoff_final` through the `_with_xdrops`
        // wrappers; see TODO.md for the investigation.
        //
        // At this low level the output depends on which seed from which
        // ungapped HSP you pick: seed=10 (first HSP, diagonal 0) yields
        // 33+5+44; seed=(81,76) (second HSP, diagonal 5) yields 35+5+42.
        // The CLI picks the latter after dedup, so end-to-end matches NCBI.
        //
        // Sequences from tests/integration.rs:3341-3345.
        fn encode(s: &[u8]) -> Vec<u8> {
            s.iter()
                .map(|&c| match c {
                    b'A' => 0,
                    b'C' => 1,
                    b'G' => 2,
                    b'T' => 3,
                    _ => 15,
                })
                .collect()
        }
        let q = encode(
            b"ACGTTGCAACGATCGTACGATTCGAGCTTAGGCTAGGGTAATCGGATCCTAGCTAGGCTAATCGATCGTAGCTAGCATCGAT",
        );
        let s = encode(
            b"ACGTTGCAACGATCGTACGATTCGAGCTTAGGCTAATCGGATCCTAGCTAGGCTAATCGATCGTAGCTAGCATCGAT",
        );
        assert_eq!(q.len(), 82);
        assert_eq!(s.len(), 77);
        // Seed = 10, matching what `blast_get_offsets_for_gapped_alignment`
        // (port of NCBI `BlastGetOffsetsForGappedAlignment`) returns for
        // this uniform-match ungapped HSP: the first sliding-window maxes
        // at index 10 (right edge of window [0..11]), strict `>` keeps the
        // first max.
        // x_dropoff=50 models blastn-short `xdrop_gap_final=100` bits after
        // `BlastExtensionParametersNew` converts with lambda≈1.374:
        // `(int)(100 * ln2 / 1.374) = 50`.
        let seed_q = 10;
        let seed_s = 10;
        let r = blast_gapped_align(&q, &s, seed_q, seed_s, 1, -3, 5, 2, 50)
            .expect("alignment should succeed");
        assert_eq!(r.score, 62, "score should be 62");
        assert_eq!((r.query_start, r.query_end), (0, 82));
        assert_eq!((r.subject_start, r.subject_end), (0, 77));
        // NCBI path: 35 Sub, 5 Ins (gap in B), 42 Sub.
        let ops = &r.edit_script.ops;
        assert_eq!(ops.len(), 3, "expected 3 op-runs, got {:?}", ops);
        assert_eq!(
            ops[0],
            (GapAlignOpType::Sub, 35),
            "first run should be 35 matches, got {:?}",
            ops[0]
        );
        assert_eq!(
            ops[1],
            (GapAlignOpType::Ins, 5),
            "middle run should be 5 insertions (gap in subject), got {:?}",
            ops[1]
        );
        assert_eq!(
            ops[2],
            (GapAlignOpType::Sub, 42),
            "last run should be 42 matches, got {:?}",
            ops[2]
        );
    }

    #[test]
    fn test_reevaluate_no_op_on_perfect_match() {
        // A 10-base perfect alignment should pass through unchanged (no
        // prefix/suffix to trim, no further matches to extend since the
        // alignment already covers the full sequences).
        let q: Vec<u8> = vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1]; // ACGTACGTAC
        let s = q.clone();
        let mut tb = TracebackResult {
            score: 10,
            edit_script: GapEditScript {
                ops: vec![(GapAlignOpType::Sub, 10)],
            },
            query_start: 0,
            query_end: 10,
            subject_start: 0,
            subject_end: 10,
        };
        let delete = reevaluate_with_ambiguities_gapped(&mut tb, &q, &s, 1, -3, 5, 2, 1);
        assert!(!delete);
        assert_eq!(tb.score, 10);
        assert_eq!(tb.edit_script.ops, vec![(GapAlignOpType::Sub, 10)]);
        assert_eq!((tb.query_start, tb.query_end), (0, 10));
        assert_eq!((tb.subject_start, tb.subject_end), (0, 10));
    }

    #[test]
    fn test_reevaluate_empty_edit_script_returns_true() {
        // An empty gap_info should always be marked for deletion.
        let q: Vec<u8> = vec![0, 1, 2, 3];
        let s: Vec<u8> = vec![0, 1, 2, 3];
        let mut tb = TracebackResult {
            score: 4,
            edit_script: GapEditScript { ops: vec![] },
            query_start: 0,
            query_end: 0,
            subject_start: 0,
            subject_end: 0,
        };
        let delete = reevaluate_with_ambiguities_gapped(&mut tb, &q, &s, 1, -3, 5, 2, 1);
        assert!(delete);
    }

    #[test]
    fn test_reevaluate_greedy_is_noop() {
        // NCBI's greedy branch (gap_open == 0 && gap_extend == 0) uses a
        // different factor=2 pricing that the Rust port does not implement.
        // The port must leave the alignment unchanged in that case.
        let q: Vec<u8> = vec![0, 1, 2, 3, 0];
        let s: Vec<u8> = vec![0, 1, 2, 3, 0];
        let mut tb = TracebackResult {
            score: 5,
            edit_script: GapEditScript {
                ops: vec![(GapAlignOpType::Sub, 5)],
            },
            query_start: 0,
            query_end: 5,
            subject_start: 0,
            subject_end: 5,
        };
        let delete = reevaluate_with_ambiguities_gapped(&mut tb, &q, &s, 1, -3, 0, 0, 1);
        assert!(!delete);
        assert_eq!(tb.score, 5);
        assert_eq!(tb.edit_script.ops, vec![(GapAlignOpType::Sub, 5)]);
    }

    #[test]
    fn test_reevaluate_deletes_if_best_score_below_cutoff() {
        // All mismatches → running sum never clears cutoff_score. The HSP
        // should be marked for deletion.
        let q: Vec<u8> = vec![0, 0, 0, 0]; // AAAA
        let s: Vec<u8> = vec![3, 3, 3, 3]; // TTTT
        let mut tb = TracebackResult {
            score: -12,
            edit_script: GapEditScript {
                ops: vec![(GapAlignOpType::Sub, 4)],
            },
            query_start: 0,
            query_end: 4,
            subject_start: 0,
            subject_end: 4,
        };
        let delete = reevaluate_with_ambiguities_gapped(&mut tb, &q, &s, 1, -3, 5, 2, 5);
        assert!(delete);
    }

    #[test]
    fn test_reevaluate_extends_right_on_exact_match() {
        // Alignment covers Q[0..4] only; two more matches available to the
        // right. Reevaluate's post-pass should extend the tail Sub run from
        // length 4 to length 6.
        let q: Vec<u8> = vec![0, 1, 2, 3, 0, 1];
        let s: Vec<u8> = vec![0, 1, 2, 3, 0, 1];
        let mut tb = TracebackResult {
            score: 4,
            edit_script: GapEditScript {
                ops: vec![(GapAlignOpType::Sub, 4)],
            },
            query_start: 0,
            query_end: 4,
            subject_start: 0,
            subject_end: 4,
        };
        let delete = reevaluate_with_ambiguities_gapped(&mut tb, &q, &s, 1, -3, 5, 2, 1);
        assert!(!delete);
        assert_eq!(tb.score, 6);
        assert_eq!(tb.edit_script.ops, vec![(GapAlignOpType::Sub, 6)]);
        assert_eq!(tb.query_end, 6);
        assert_eq!(tb.subject_end, 6);
    }

    #[test]
    fn test_reevaluate_trims_trailing_below_cutoff() {
        // Alignment = 10 perfect matches + 5 mismatches. Running sum peaks
        // at 10 after the first 10 Subs, then decreases. Reevaluate should
        // trim the trailing 5 mismatches and return a clean 10-Sub alignment.
        let q: Vec<u8> = vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2];
        let s: Vec<u8> = vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 0, 2, 1, 3, 0];
        let mut tb = TracebackResult {
            score: -5, // 10 - 5*3 = -5 with reward=1 penalty=-3
            edit_script: GapEditScript {
                ops: vec![(GapAlignOpType::Sub, 15)],
            },
            query_start: 0,
            query_end: 15,
            subject_start: 0,
            subject_end: 15,
        };
        let delete = reevaluate_with_ambiguities_gapped(&mut tb, &q, &s, 1, -3, 5, 2, 1);
        assert!(!delete);
        assert_eq!(tb.score, 10);
        assert_eq!(tb.edit_script.ops, vec![(GapAlignOpType::Sub, 10)]);
        assert_eq!((tb.query_start, tb.query_end), (0, 10));
    }

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
        assert!(
            r.subject_start >= 100 && r.subject_start <= 104,
            "Subject start should be near 100, got {}",
            r.subject_start
        );
        assert!(
            r.subject_end >= 104 && r.subject_end <= 108,
            "Subject end should be near 108, got {}",
            r.subject_end
        );
    }

    #[test]
    fn test_blast_gapped_align_basic() {
        let q = vec![0u8, 1, 2, 3, 0, 1, 2, 3]; // ACGTACGT
        let s = vec![0u8, 1, 2, 3, 0, 1, 2, 3];
        let r = blast_gapped_align(&q, &s, 4, 4, 1, -3, 5, 2, 30);
        assert!(r.is_some(), "Should find alignment");
        let r = r.unwrap();
        eprintln!(
            "blast_gapped: score={} q={}..{} s={}..{} ops={:?}",
            r.score, r.query_start, r.query_end, r.subject_start, r.subject_end, r.edit_script.ops
        );
        assert!(r.score > 0, "score={}", r.score);
        // Check edit script has content
        let total_ops: i32 = r.edit_script.ops.iter().map(|(_, n)| *n).sum();
        assert!(
            total_ops > 0,
            "edit script should have operations, got {:?}",
            r.edit_script.ops
        );
    }

    #[test]
    fn test_blast_gapped_align_exact_match_extends_to_edges() {
        let q = b"GAATCCATGCTGTGGGCCAGCAAGAGTTAAGGTGCTCATGGTTTTGAGAAAACATCTGAGGACTCTGACAGCACTCTCCCATCCTTGGTCTCCACAGTCT"
            .iter()
            .map(|b| match b {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => 15,
            })
            .collect::<Vec<u8>>();
        let r =
            blast_gapped_align(&q, &q, 50, 50, 1, -3, 5, 2, 20).expect("exact match should align");
        assert_eq!(r.score, 100);
        assert_eq!((r.query_start, r.query_end), (0, 100));
        assert_eq!((r.subject_start, r.subject_end), (0, 100));
    }

    #[test]
    fn test_blast_gapped_align_single_internal_gap() {
        let q: Vec<u8> = [
            0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
        ]
        .to_vec();
        let s: Vec<u8> = [
            0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 0, 0, 0, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
        ]
        .to_vec();

        let r = blast_gapped_align(&q, &s, 12, 12, 2, -3, 3, 1, 30)
            .expect("single-gap alignment should succeed");
        eprintln!(
            "blast_gapped_gap: score={} q={}..{} s={}..{} ops={:?}",
            r.score, r.query_start, r.query_end, r.subject_start, r.subject_end, r.edit_script.ops
        );
        assert!(
            r.edit_script.ops.iter().any(|(op, _)| matches!(
                op,
                GapAlignOpType::Del | GapAlignOpType::Del1 | GapAlignOpType::Del2
            )),
            "expected one gap-in-query operation, got {:?}",
            r.edit_script.ops
        );
        assert!(
            r.edit_script.ops.len() <= 5,
            "single-gap case should not fragment heavily, got {:?}",
            r.edit_script.ops
        );
    }

    #[test]
    fn test_align_ex_forward() {
        let mut matrix = [[-3i32; 16]; 16];
        for i in 0..4 {
            matrix[i][i] = 1;
        }
        let a = vec![255u8, 0, 1, 2, 3]; // pad + ACGT
        let b = vec![255u8, 0, 1, 2, 3, 255]; // pad + ACGT + pad
        let (score, ao, bo, ops) = align_ex(&a, &b, 4, 4, &matrix, 5, 2, 30, false);
        assert_eq!(score, 4, "4 matches * reward 1");
        assert_eq!(ao, 4);
        assert_eq!(bo, 4);
        assert_eq!(ops, vec![(GapAlignOpType::Sub, 4)]);
    }

    /// Alignment where query has a deletion relative to subject (gap in query).
    /// Query:   ACGTACGTACGT----ACGTACGTACGT
    /// Subject: ACGTACGTACGTAAAAACGTACGTACGT
    /// Long flanking matches make the gap worthwhile.
    #[test]
    fn test_traceback_with_gap_in_query() {
        // 12 matching bases + 4 extra in subject + 12 matching bases
        let q: Vec<u8> = [
            0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
        ]
        .to_vec(); // 24 bp
        let s: Vec<u8> = [
            0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 0, 0, 0, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
        ]
        .to_vec(); // 28 bp
        let (score, esp, _, _, _, _) = traceback_align(&q, &s, 0, q.len(), 0, s.len(), 2, -3, 3, 1);
        assert!(score > 0, "should find alignment, score={}", score);
        // The edit script should contain an Ins operation (extra subject bases => gap in query)
        // In gapinfo convention: Ins = "insertion in query" = advancing subject without query
        let has_ins = esp.ops.iter().any(|(op, _)| *op == GapAlignOpType::Ins);
        assert!(has_ins, "expected Ins (gap in query), ops={:?}", esp.ops);
    }

    /// Alignment where subject has a deletion relative to query (gap in subject).
    /// Query:   ACGTACGTACGTAAAAACGTACGTACGT  (28 bp, extra A's in middle)
    /// Subject: ACGTACGTACGT----ACGTACGTACGT  (24 bp)
    #[test]
    fn test_traceback_with_gap_in_subject() {
        let q: Vec<u8> = [
            0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 0, 0, 0, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
        ]
        .to_vec(); // 28 bp
        let s: Vec<u8> = [
            0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
        ]
        .to_vec(); // 24 bp
        let (score, esp, _, _, _, _) = traceback_align(&q, &s, 0, q.len(), 0, s.len(), 2, -3, 3, 1);
        assert!(score > 0, "should find alignment, score={}", score);
        // The edit script should contain a Del operation (extra query bases => gap in subject)
        // In gapinfo convention: Del = "deletion in query" = advancing query without subject
        let has_del = esp.ops.iter().any(|(op, _)| *op == GapAlignOpType::Del);
        assert!(has_del, "expected Del (gap in subject), ops={:?}", esp.ops);
    }

    /// Alignment requiring multiple gaps.
    /// Query:   ACG--TTA--CGT
    /// Subject: ACGAATTACCGCGT
    #[test]
    fn test_traceback_multiple_gaps() {
        // Query is shorter; subject has extra bases in two places
        let q = vec![0u8, 1, 2, 3, 3, 0, 1, 2, 3]; // ACGTTACGT (9 bp)
        let s = vec![0u8, 1, 2, 0, 0, 3, 3, 0, 1, 1, 2, 3]; // ACGAATTACCGT (12 bp)
        let (score, esp, _, _, _, _) = traceback_align(&q, &s, 0, q.len(), 0, s.len(), 2, -3, 3, 1);
        assert!(score > 0, "should find alignment, score={}", score);
        // Count gap operations (Del = gap in query)
        let gap_ops: Vec<_> = esp
            .ops
            .iter()
            .filter(|(op, _)| *op == GapAlignOpType::Del || *op == GapAlignOpType::Ins)
            .collect();
        assert!(
            !gap_ops.is_empty(),
            "expected at least one gap operation, ops={:?}",
            esp.ops
        );
    }

    /// Verify start/end coordinates are correct for a known alignment.
    #[test]
    fn test_traceback_coordinates() {
        // Embed a match region at positions 5..13 in both sequences, surrounded by mismatches
        let mut q = vec![3u8; 20]; // all T
        let mut s = vec![2u8; 20]; // all G (mismatches with T)
                                   // Place ACGTACGT at q[5..13] and s[5..13]
        for (i, &b) in [0u8, 1, 2, 3, 0, 1, 2, 3].iter().enumerate() {
            q[5 + i] = b;
            s[5 + i] = b;
        }
        let (score, esp, aq_start, aq_end, as_start, as_end) =
            traceback_align(&q, &s, 0, q.len(), 0, s.len(), 2, -3, 5, 2);
        assert_eq!(score, 16, "8 matches * reward 2");
        // The alignment should span the matching region [5..13)
        assert_eq!(aq_start, 5, "query alignment start");
        assert_eq!(aq_end, 13, "query alignment end");
        assert_eq!(as_start, 5, "subject alignment start");
        assert_eq!(as_end, 13, "subject alignment end");
        // Edit script should be a single Sub run of length 8
        assert_eq!(esp.ops.len(), 1);
        assert_eq!(esp.ops[0], (GapAlignOpType::Sub, 8));
    }

    /// The traceback score should match what the DP matrix says.
    /// Recompute expected score from edit script + sequences.
    #[test]
    fn test_traceback_score_matches_dp() {
        let q = vec![0u8, 1, 2, 3, 3, 0, 1, 2, 3]; // ACGTTACGT
        let s = vec![0u8, 1, 2, 0, 3, 0, 1, 2, 3]; // ACGATACGT
        let reward = 2i32;
        let penalty = -3i32;
        let gap_open = 5i32;
        let gap_ext = 2i32;
        let (score, esp, aq_start, _, as_start, _) = traceback_align(
            &q,
            &s,
            0,
            q.len(),
            0,
            s.len(),
            reward,
            penalty,
            gap_open,
            gap_ext,
        );
        assert!(score > 0);
        // Recompute score from edit script
        let mut computed_score = 0i32;
        let mut qi = aq_start;
        let mut si = as_start;
        for (op, count) in &esp.ops {
            match op {
                GapAlignOpType::Sub => {
                    for _ in 0..*count {
                        computed_score += if q[qi] == s[si] { reward } else { penalty };
                        qi += 1;
                        si += 1;
                    }
                }
                GapAlignOpType::Del => {
                    // Gap in query: consume subject
                    computed_score -= gap_open + gap_ext * count;
                    qi += *count as usize;
                }
                GapAlignOpType::Ins => {
                    // Gap in subject: consume query
                    computed_score -= gap_open + gap_ext * count;
                    si += *count as usize;
                }
                _ => {}
            }
        }
        assert_eq!(
            score, computed_score,
            "DP score {} should match recomputed score {} from edit script {:?}",
            score, computed_score, esp.ops
        );
    }

    #[test]
    fn test_align_ex_reverse() {
        let mut matrix = [[-3i32; 16]; 16];
        for i in 0..4 {
            matrix[i][i] = 2;
        } // reward=2 for stronger diagonal signal
        let a = vec![0u8, 1, 2, 3, 0];
        let b = vec![0u8, 1, 2, 3, 0];
        let (score, ao, _bo, ops) = align_ex(&a, &b, 5, 5, &matrix, 5, 2, 30, true);
        eprintln!("rev2: score={} ao={} ops={:?}", score, ao, ops);
        assert!(score >= 8, "should get at least 4 bases * 2, got {}", score);
    }

    fn pack_ncbi2na(bases: &[u8]) -> Vec<u8> {
        let mut packed = vec![0u8; bases.len().div_ceil(4)];
        for (i, &base) in bases.iter().enumerate() {
            packed[i >> 2] |= (base & 0x03) << (6 - 2 * (i & 3));
        }
        packed
    }

    #[test]
    fn test_gapped_score_one_dir_packed_subject_matches_decoded_basic() {
        let query = vec![0u8, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
        let subject = vec![0u8, 1, 2, 3, 0, 1, 3, 3, 0, 1, 2, 3];
        let packed = pack_ncbi2na(&subject);
        let matrix = build_blastna_matrix(1, -3);
        let left_decoded =
            gapped_score_one_dir(&query[..6], &subject[..6], 6, 6, &matrix, 7, 2, 12, true);
        let left_packed = gapped_score_one_dir_packed_subject(
            &query[..6],
            &packed,
            6,
            6,
            0,
            &matrix,
            7,
            2,
            12,
            true,
        );
        assert_eq!(left_packed, left_decoded);

    }
}
