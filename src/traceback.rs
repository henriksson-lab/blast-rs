//! Gapped alignment with traceback.
//!
//! Mirrors NCBI's `blast_gapalign.c` (the ALIGN_EX / gapped-alignment half) and
//! `blast_traceback.c` (the HSP-list traceback half). Function-name mapping
//! maintained for side-by-side auditing:
//!
//! | Rust                             | NCBI source (blast_gapalign.c)           |
//! | -------------------------------- | ---------------------------------------- |
//! | `align_ex`                       | `ALIGN_EX` (line 374)                    |
//! | `blast_gapped_alignment_with_traceback` | `BLAST_GappedAlignmentWithTraceback`     |
//! | `blast_semi_gapped_align`        | `Blast_SemiGappedAlign`                  |
//! | `s_blast_dyn_prog_nt_gapped_alignment` | `s_BlastDynProgNtGappedAlignment`        |
//! | `s_blast_align_packed_nucl`      | `s_BlastAlignPackedNucl`                 |
//! | `gapped_score_one_dir`           | inner recurrence of `Blast_SemiGappedAlign` |
//! | `blast_score_blk_nucl_matrix_create` | `BlastScoreBlkNuclMatrixCreate` (blast_stat.c) |
//! | `traceback_align` / `_abs`       | standalone NW helper, no C analog        |
//!
//! Script/op constants match NCBI exactly (see `SCRIPT_SUB` etc.) —
//! `blast_gapalign.c:363-371` and `gapinfo.h:45-51`.

pub use crate::extend::PhiInitialHit;
use crate::gapinfo::{GapAlignOpType, GapEditScript, GapPrelimEditBlock};

#[derive(Debug, Clone)]
pub struct OutOfFrameScoring {
    pub matrix: Vec<Vec<i32>>,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub shift_penalty: i32,
    pub x_dropoff: i32,
    pub pssm: Option<Vec<Vec<i32>>>,
}

impl OutOfFrameScoring {
    fn matrix_score(&self, a: u8, b: u8) -> i32 {
        self.matrix
            .get(a as usize)
            .and_then(|row| row.get(b as usize))
            .copied()
            .unwrap_or(if a == b { 1 } else { -1 })
    }

    fn oof_row_score(
        &self,
        query_offset: i32,
        a: &[u8],
        m: usize,
        a_index: usize,
        b_letter: u8,
        reversed: bool,
    ) -> i32 {
        if let Some(pssm) = &self.pssm {
            let row = if reversed {
                m.saturating_sub(a_index)
            } else {
                (a_index as i32 + query_offset).max(0) as usize
            };
            return pssm
                .get(row)
                .and_then(|row| row.get(b_letter as usize))
                .copied()
                .unwrap_or_else(|| self.matrix_score(*a.get(a_index).unwrap_or(&0), b_letter));
        }
        let row = if reversed {
            m.saturating_sub(a_index)
        } else {
            a_index
        };
        self.matrix_score(*a.get(row).unwrap_or(&0), b_letter)
    }
}

impl Default for OutOfFrameScoring {
    fn default() -> Self {
        Self {
            matrix: crate::matrix::BLOSUM62
                .iter()
                .map(|row| row.to_vec())
                .collect(),
            gap_open: 11,
            gap_extend: 1,
            shift_penalty: 1,
            x_dropoff: 25,
            pssm: None,
        }
    }
}

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
    let mut e = vec![vec![i32::MIN / 2; cols]; m + 1]; // gap in subject
    let mut f = vec![vec![i32::MIN / 2; cols]; m + 1]; // gap in query

    // Traceback direction matrix
    // 0=diag, 1=left(del), 2=up(ins)
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
                trace[i][j] = 1; // left (gap in query)
            } else {
                trace[i][j] = 2; // up (gap in subject)
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
                // Left consumes subject only: gap in query.
                if let Some(last) = ops.last_mut() {
                    if last.0 == GapAlignOpType::Del {
                        last.1 += 1;
                    } else {
                        ops.push((GapAlignOpType::Del, 1));
                    }
                } else {
                    ops.push((GapAlignOpType::Del, 1));
                }
                j -= 1;
            }
            2 => {
                // Up consumes query only: gap in subject.
                if let Some(last) = ops.last_mut() {
                    if last.0 == GapAlignOpType::Ins {
                        last.1 += 1;
                    } else {
                        ops.push((GapAlignOpType::Ins, 1));
                    }
                } else {
                    ops.push((GapAlignOpType::Ins, 1));
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
#[derive(Clone)]
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
const SCRIPT_AHEAD_ONE_FRAME: u8 = 1;
const SCRIPT_AHEAD_TWO_FRAMES: u8 = 2;
const SCRIPT_NEXT_IN_FRAME: u8 = 3;
const SCRIPT_NEXT_PLUS_ONE_FRAME: u8 = 4;
const SCRIPT_NEXT_PLUS_TWO_FRAMES: u8 = 5;
const SCRIPT_OP_MASK: u8 = 0x07;
const SCRIPT_OOF_OPEN_GAP: u8 = 0x08;
const SCRIPT_EXTEND_GAP_A: u8 = 0x10;
const SCRIPT_EXTEND_GAP_B: u8 = 0x40;

#[derive(Clone, Copy)]
struct GapDP {
    best: i32,
    best_gap: i32,
}

#[derive(Default)]
struct AlignExScratch {
    sa: Vec<GapDP>,
    scripts: Vec<u8>,
    row_offsets: Vec<usize>,
    row_lens: Vec<usize>,
    row_starts: Vec<usize>,
}

thread_local! {
    static ALIGN_EX_SCRATCH: std::cell::RefCell<AlignExScratch> =
        std::cell::RefCell::new(AlignExScratch::default());
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
///
/// `pub(crate)` so kappa-side callers (`blast_kappa::sw_find_final_ends_using_xdrop`)
/// can drive the traceback without re-implementing the DP. The
/// signature is identical to the C function's effective output:
/// `(best_score, a_offset, b_offset, edit_ops)`.
pub(crate) fn align_ex(
    a: &[u8],
    b: &[u8],
    m: usize,
    n: usize,
    matrix: &[[i32; 16]; 16],
    gap_open: i32,
    gap_extend: i32,
    x_dropoff: i32,
    reverse: bool,
) -> (i32, usize, usize, Vec<(GapAlignOpType, i32)>) {
    ALIGN_EX_SCRATCH.with(|scratch| {
        let mut scratch = scratch.borrow_mut();
        align_ex_with_scratch(
            a,
            b,
            m,
            n,
            matrix,
            gap_open,
            gap_extend,
            x_dropoff,
            reverse,
            &mut scratch,
        )
    })
}

#[allow(clippy::too_many_arguments)]
fn align_ex_with_scratch(
    a: &[u8],
    b: &[u8],
    m: usize,
    n: usize,
    matrix: &[[i32; 16]; 16],
    gap_open: i32,
    gap_extend: i32,
    mut x_dropoff: i32,
    reverse: bool,
    scratch: &mut AlignExScratch,
) -> (i32, usize, usize, Vec<(GapAlignOpType, i32)>) {
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
    // NCBI `ALIGN_EX` (`blast_gapalign.c:463-469`) sizes `dp_mem` to cover
    // any position that can survive the X-drop check, starting at
    // `num_extra_cells` and growing via `realloc` (same file, 3276-3284).
    // The score-only paths match that; `align_ex` must too, otherwise any
    // alignment where `b_size` would exceed the initial allocation gets
    // truncated silently.
    let AlignExScratch {
        sa,
        scripts,
        row_offsets,
        row_lens,
        row_starts,
    } = scratch;
    sa.clear();
    sa.resize(
        n + extra + 10,
        GapDP {
            best: MININT,
            best_gap: MININT,
        },
    );

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
    //
    // Store rows in one flat buffer. The older Vec<Vec<u8>> shape allocated once
    // per DP row, which is much heavier than the C gap_align scratch buffers used
    // by ALIGN_EX.
    let row_capacity_hint = (extra + 16).min(n + extra + 10);
    scripts.clear();
    scripts.reserve((m + 1).saturating_mul(row_capacity_hint));
    row_offsets.clear();
    row_offsets.reserve(m + 1);
    row_lens.clear();
    row_lens.reserve(m + 1);
    row_starts.clear();
    row_starts.reserve(m + 1);
    row_offsets.push(0);
    row_lens.push(b_size);
    scripts.resize(b_size, SCRIPT_GAP_IN_A);
    row_starts.push(0);

    let mut best_score = 0i32;
    let mut first_b = 0usize;
    let mut a_off = 0usize;
    let mut b_off = 0usize;

    for ai in 1..=m {
        let a_letter = if reverse { a[m - ai] } else { a[ai] };
        let mrow = &matrix[a_letter as usize & 0x0F];

        let row_start = first_b;
        let row_offset = scripts.len();
        let row_len = (b_size - row_start) + extra + 10;
        scripts.resize(row_offset + row_len, 0);
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
                // NCBI `ALIGN_EX` (`blast_gapalign.c:611`) uses strict `>`:
                // ties on `best_score` keep the FIRST cell encountered.
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
            if script_idx < row_len {
                scripts[row_offset + script_idx] = script;
            }
        }

        if first_b >= b_size {
            row_offsets.push(row_offset);
            row_lens.push(row_len);
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
                if script_idx < row_len {
                    scripts[row_offset + script_idx] = SCRIPT_GAP_IN_A;
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
        row_offsets.push(row_offset);
        row_lens.push(row_len);
        row_starts.push(row_start);
    }

    // Traceback
    let mut ops: Vec<(GapAlignOpType, i32)> = Vec::new();
    let mut ai = a_off;
    let mut bi = b_off;
    let mut cur_script = SCRIPT_SUB;

    while ai > 0 || bi > 0 {
        if ai >= row_starts.len() {
            break;
        }
        let row_start = row_starts[ai];
        if bi < row_start {
            break;
        }
        let script_idx = bi - row_start;
        if script_idx >= row_lens[ai] {
            break;
        }
        let s = scripts[row_offsets[ai] + script_idx];

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
    (best_score, a_off, b_off, ops)
}

/// Build full BLASTNA scoring matrix — verbatim port of NCBI
/// `BlastScoreBlkNuclMatrixCreate` (`blast_stat.c:1060`).
pub fn blast_score_blk_nucl_matrix_create(reward: i32, penalty: i32) -> [[i32; 16]; 16] {
    use crate::encoding::{blastna_pair_score, BLASTNA_SIZE};

    let mut matrix = [[0i32; 16]; BLASTNA_SIZE];

    for index1 in 0..BLASTNA_SIZE {
        for index2 in index1..BLASTNA_SIZE {
            let s = blastna_pair_score(index1 as u8, index2 as u8, reward, penalty);
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
pub fn blast_semi_gapped_align(
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
    let matrix = blast_score_blk_nucl_matrix_create(reward, penalty);

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

const RESTRICT_SIZE: usize = 10;

fn restricted_gapped_score_one_dir(
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
    let mut b_gap = 0usize;

    for ai in 1..=m {
        let a_letter = if reverse { a[m - ai] } else { a[ai] };
        let mrow = &matrix[a_letter as usize & 0x0F];
        let mut sc = MININT;
        let mut sgr = MININT;
        let mut last_b = first_b;
        let allow_row_gap_start = ai % RESTRICT_SIZE == 0;

        #[allow(clippy::mut_range_bound)]
        for bi in first_b..b_size {
            let b_idx = if reverse {
                n.checked_sub(1 + bi).unwrap_or(usize::MAX)
            } else {
                bi + 1
            };
            let next_sc = if b_idx < b.len() {
                sa[bi].best + mrow[b[b_idx] as usize & 0x0F]
            } else {
                MININT
            };
            let allow_col_gap_start = bi == b_gap;

            if allow_col_gap_start {
                b_gap += RESTRICT_SIZE;
                let sgc = sa[bi].best_gap;
                if sc < sgc {
                    sc = sgc;
                }
            }
            if allow_row_gap_start && sc < sgr {
                sc = sgr;
            }

            if best_score - sc > x_dropoff {
                sa[bi].best = MININT;
                if bi == first_b {
                    first_b += 1;
                }
            } else {
                last_b = bi;
                if sc > best_score {
                    best_score = sc;
                }

                if allow_col_gap_start {
                    let sgc = sa[bi].best_gap - gap_extend;
                    sa[bi].best_gap = (sc - gap_oe).max(sgc);
                }
                if allow_row_gap_start {
                    sgr = (sc - gap_oe).max(sgr - gap_extend);
                }
                sa[bi].best = sc;
            }
            sc = next_sc;
        }

        if first_b >= b_size {
            break;
        }

        b_gap = first_b;
        let remainder = first_b % RESTRICT_SIZE;
        if remainder > 0 {
            b_gap += RESTRICT_SIZE - remainder;
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
    }

    best_score
}

/// NCBI: s_RestrictedGappedAlign (`blast_gapalign.c`).
///
/// Score-only nucleotide extension with the Cameron/Williams/Cannane restricted
/// gap rule used by NCBI BLAST: a gap may start only at offsets divisible by
/// `RESTRICT_SIZE` from the seed.
#[allow(clippy::too_many_arguments)]
pub fn s_restricted_gapped_align(
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
    let matrix = blast_score_blk_nucl_matrix_create(reward, penalty);

    let score_l = restricted_gapped_score_one_dir(
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

    let score_r = if seed_q < query.len() - 1 && seed_s < subject.len() - 1 {
        restricted_gapped_score_one_dir(
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

/// Name-aligned wrapper for NCBI `s_BlastDynProgNtGappedAlignment`
/// (`blast_gapalign.c:2949`) on the packed-subject preliminary score path.
///
/// Audit note: this is a name-alignment shim over the existing Rust packed
/// preliminary DP path, not a line-for-line port of the full C function. Treat
/// it as potentially incomplete/hard to audit until the complete C control
/// flow around prelim HSP construction has been compared.
pub fn s_blast_dyn_prog_nt_gapped_alignment(
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
) -> (
    i32,
    usize,
    usize,
    usize,
    usize,
    usize,
    usize,
    usize,
    usize,
    usize,
    usize,
    i32,
    i32,
) {
    s_blast_dyn_prog_nt_gapped_alignment_extents(
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
}

/// blast-rs: shared preliminary-extent helper for
/// `s_blast_dyn_prog_nt_gapped_alignment`; not a direct NCBI C port.
fn s_blast_dyn_prog_nt_gapped_alignment_extents(
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
) -> (
    i32,
    usize,
    usize,
    usize,
    usize,
    usize,
    usize,
    usize,
    usize,
    usize,
    usize,
    i32,
    i32,
) {
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

    let matrix = blast_score_blk_nucl_matrix_create(reward, penalty);
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

    let (score_right, query_stop, subject_stop) =
        if q_length < query.len() && s_length < subject_len {
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

/// Name-aligned wrapper for NCBI `s_BlastAlignPackedNucl`
/// (`blast_gapalign.c:3034`).
///
/// Audit note: this delegates to the existing one-direction packed-subject DP.
/// It exists to make the C function name visible to the audit, but should not
/// by itself be considered proof that every C caller/context has been ported.
pub fn s_blast_align_packed_nucl(
    query: &[u8],
    subject_packed: &[u8],
    n: usize,
    m: usize,
    subject_base_offset: usize,
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    x_dropoff: i32,
    reverse: bool,
) -> (i32, usize, usize) {
    let matrix = blast_score_blk_nucl_matrix_create(reward, penalty);
    gapped_score_one_dir_packed_subject(
        query,
        subject_packed,
        n,
        m,
        subject_base_offset,
        &matrix,
        gap_open + gap_extend,
        gap_extend,
        x_dropoff,
        reverse,
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

            let sgc = sa[bi].best_gap;
            let next_sc = if b_idx < b.len() {
                sa[bi].best + mrow[b[b_idx] as usize & 0x0F]
            } else {
                MININT
            };

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
                } else {
                    sa[bi].best = MININT;
                }
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

/// Extent-returning variant of `gapped_score_one_dir` for the 16x16 nucleotide
/// matrix. Identical DP to `gapped_score_one_dir` (the faithful `ALIGN_EX`
/// score-only core, blast_gapalign.c `s_BlastAlignPackedNucl`), but additionally
/// tracks the offsets `(a_off, b_off)` of the maximum-scoring cell so the caller
/// can recover the gapped alignment extent — mirroring how the packed-subject
/// path (`gapped_score_one_dir_packed_subject`) and the generic-matrix path
/// (`gapped_score_one_dir_generic`) return `b_offset`/`a_offset` to
/// `s_BlastDynProgNtGappedAlignment`.
#[allow(clippy::too_many_arguments)]
fn gapped_score_one_dir_ex(
    a: &[u8],
    b: &[u8],
    m: usize,
    n: usize,
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
    let mut first_b = 0usize;
    let mut a_off = 0usize;
    let mut b_off = 0usize;

    for ai in 1..=m {
        let a_letter = if reverse { a[m - ai] } else { a[ai] };
        let mrow = &matrix[a_letter as usize & 0x0F];
        let mut sc = MININT;
        let mut sgr = MININT;
        let mut last_b = first_b;

        #[allow(clippy::mut_range_bound)]
        for bi in first_b..b_size {
            let b_idx = if reverse {
                n.checked_sub(1 + bi).unwrap_or(usize::MAX)
            } else {
                bi + 1
            };

            let sgc = sa[bi].best_gap;
            let next_sc = if b_idx < b.len() {
                sa[bi].best + mrow[b[b_idx] as usize & 0x0F]
            } else {
                MININT
            };

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
                    a_off = ai;
                    b_off = bi;
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

    (best_score, a_off, b_off)
}

/// Decoded-subject analog of `s_blast_dyn_prog_nt_gapped_alignment`
/// (NCBI `s_BlastDynProgNtGappedAlignment`, blast_gapalign.c:3006), returning
/// the gapped alignment score together with its query/subject extent. Unlike the
/// packed variant it operates on a one-base-per-byte subject, so the packed
/// 4-byte `offset_adjustment` (which exists only because `s_BlastAlignPackedNucl`
/// requires byte-aligned subject offsets) is NOT applied; the seed coordinates
/// are used directly. The X-drop clamp to `ungapped_score` and the
/// left/right two-directional extension match the C exactly.
#[allow(clippy::too_many_arguments)]
pub fn s_blast_dyn_prog_nt_gapped_alignment_decoded_extents(
    query: &[u8],
    subject: &[u8],
    seed_q: usize,
    seed_s: usize,
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    x_dropoff: i32,
    ungapped_score: i32,
) -> (i32, usize, usize, usize, usize) {
    let mut effective_x_dropoff = x_dropoff;
    if ungapped_score < effective_x_dropoff {
        effective_x_dropoff = ungapped_score;
    }

    let matrix = blast_score_blk_nucl_matrix_create(reward, penalty);
    let gap_oe = gap_open + gap_extend;

    // Left extension (reverse). `gapped_score_one_dir_ex` returns
    // (score, a_off, b_off) where `a` is the query and `b` the subject, so
    // a_off is the query-direction offset (1-based `ai`) and b_off the
    // subject-direction offset (0-based `bi`) of the best cell — mirroring the
    // `private_q_start`/`private_s_start` outputs of `s_BlastAlignPackedNucl`.
    let q_len_left = seed_q + 1;
    let s_len_left = seed_s + 1;
    let (score_left, private_q_start, private_s_start) = gapped_score_one_dir_ex(
        &query[..q_len_left],
        &subject[..s_len_left],
        q_len_left,
        s_len_left,
        &matrix,
        gap_oe,
        gap_extend,
        effective_x_dropoff,
        true,
    );
    let query_start = q_len_left.saturating_sub(private_q_start);
    let subject_start = s_len_left.saturating_sub(private_s_start);

    // Right extension (forward), score + extent.
    let (score_right, query_stop, subject_stop) =
        if seed_q < query.len() - 1 && seed_s < subject.len() - 1 {
            let (sr, q_stop, s_stop) = gapped_score_one_dir_ex(
                &query[seed_q..],
                &subject[seed_s..],
                query.len() - seed_q - 1,
                subject.len() - seed_s - 1,
                &matrix,
                gap_oe,
                gap_extend,
                effective_x_dropoff,
                false,
            );
            (sr, seed_q + q_stop, seed_s + s_stop)
        } else {
            (0, seed_q, seed_s)
        };

    (
        score_left + score_right,
        query_start,
        query_stop,
        subject_start,
        subject_stop,
    )
}

fn generic_matrix_score(matrix: &[Vec<i32>], a: u8, b: u8) -> i32 {
    matrix
        .get(a as usize)
        .and_then(|row| row.get(b as usize))
        .copied()
        .unwrap_or(0)
}

fn prelim_block_set_ops(
    block: &mut crate::gapinfo::GapPrelimEditBlock,
    ops: Vec<(GapAlignOpType, i32)>,
) {
    crate::gapinfo::gap_prelim_edit_block_reset(Some(block));
    for (op, count) in ops {
        crate::gapinfo::gap_prelim_edit_block_add(block, op, count);
    }
}

fn align_ex_generic(
    a: &[u8],
    b: &[u8],
    m: usize,
    n: usize,
    matrix: &[Vec<i32>],
    gap_open: i32,
    gap_extend: i32,
    mut x_dropoff: i32,
    reverse: bool,
) -> (i32, usize, usize, Vec<(GapAlignOpType, i32)>) {
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
    let mut sa = vec![
        GapDP {
            best: MININT,
            best_gap: MININT
        };
        n + extra + 10
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

    let mut scripts: Vec<Vec<u8>> = Vec::with_capacity(m + 1);
    let mut row_starts: Vec<usize> = Vec::with_capacity(m + 1);
    scripts.push(vec![SCRIPT_GAP_IN_A; b_size]);
    row_starts.push(0);

    let mut best_score = 0i32;
    let mut first_b = 0usize;
    let mut a_off = 0usize;
    let mut b_off = 0usize;

    for ai in 1..=m {
        let a_idx = if reverse { m - ai } else { ai };
        let Some(&a_letter) = a.get(a_idx) else {
            break;
        };
        let row_start = first_b;
        let mut row_script = vec![0u8; (b_size - row_start) + extra + 10];
        let mut sc = MININT;
        let mut sgr = MININT;
        let mut last_b = first_b;

        #[allow(clippy::mut_range_bound)]
        for bi in first_b..b_size {
            let b_idx = if reverse {
                n.checked_sub(1 + bi).unwrap_or(usize::MAX)
            } else {
                bi + 1
            };
            let sgc = sa[bi].best_gap;
            let next_sc = if let Some(&b_letter) = b.get(b_idx) {
                sa[bi].best + generic_matrix_score(matrix, a_letter, b_letter)
            } else {
                MININT
            };

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

        if first_b >= b_size {
            scripts.push(row_script);
            row_starts.push(row_start);
            break;
        }

        if last_b < b_size - 1 {
            b_size = last_b + 1;
        } else {
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

    let mut ops: Vec<(GapAlignOpType, i32)> = Vec::new();
    let mut ai = a_off;
    let mut bi = b_off;
    let mut cur_script = SCRIPT_SUB;

    while ai > 0 || bi > 0 {
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

    (best_score, a_off, b_off, ops)
}

fn gapped_score_one_dir_generic(
    a: &[u8],
    b: &[u8],
    m: usize,
    n: usize,
    matrix: &[Vec<i32>],
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
    let mut first_b = 0usize;
    let mut a_off = 0usize;
    let mut b_off = 0usize;

    for ai in 1..=m {
        let a_idx = if reverse { m - ai } else { ai };
        let Some(&a_letter) = a.get(a_idx) else {
            break;
        };
        let mut sc = MININT;
        let mut sgr = MININT;
        let mut last_b = first_b;

        #[allow(clippy::mut_range_bound)]
        for bi in first_b..b_size {
            let b_idx = if reverse {
                n.checked_sub(1 + bi).unwrap_or(usize::MAX)
            } else {
                bi + 1
            };

            let sgc = sa[bi].best_gap;
            let next_sc = b
                .get(b_idx)
                .map(|&b_letter| sa[bi].best + generic_matrix_score(matrix, a_letter, b_letter))
                .unwrap_or(MININT);

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
                    a_off = ai;
                    b_off = bi;
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

    (best_score, a_off, b_off)
}

/// NCBI: s_PHIGappedAlignment (`phi_gapalign.c:670`).
#[allow(clippy::too_many_arguments)]
pub fn s_phi_gapped_alignment(
    query: &[u8],
    subject: &[u8],
    gap_align: Option<&mut crate::blast_kappa::BlastGapAlignWorkspace>,
    matrix: &[Vec<i32>],
    gap_open: i32,
    gap_extend: i32,
    x_dropoff: i32,
    query_offset: i32,
    subject_offset: i32,
    query_length: i32,
    subject_length: i32,
) -> i16 {
    let Some(gap_align) = gap_align else {
        return -1;
    };
    if query_offset < 0 || subject_offset < 0 || query_length < 0 || subject_length < 0 {
        return -1;
    }

    let mut q_length = query_offset as usize;
    let mut s_length = subject_offset as usize;
    let gap_oe = gap_open + gap_extend;

    let mut found_start = false;
    let mut found_end = false;
    let mut score_left = 0;

    if q_length != 0 && s_length != 0 {
        found_start = true;
        let (score, private_q_start, private_s_start) = gapped_score_one_dir_generic(
            query, subject, q_length, s_length, matrix, gap_oe, gap_extend, x_dropoff, true,
        );
        score_left = score;
        gap_align.query_start = (q_length.saturating_sub(private_q_start) + 1) as i32;
        gap_align.subject_start = (s_length.saturating_sub(private_s_start) + 1) as i32;
    }

    q_length += query_length.saturating_sub(1) as usize;
    s_length += subject_length.saturating_sub(1) as usize;

    let mut score_right = 0;
    if q_length < query.len() && s_length < subject.len() {
        found_end = true;
        let (score, query_stop, subject_stop) = gapped_score_one_dir_generic(
            query.get(q_length..).unwrap_or(&[]),
            subject.get(s_length..).unwrap_or(&[]),
            query.len().saturating_sub(q_length + 1),
            subject.len().saturating_sub(s_length + 1),
            matrix,
            gap_oe,
            gap_extend,
            x_dropoff,
            false,
        );
        score_right = score;
        gap_align.query_stop = (query_stop + q_length) as i32;
        gap_align.subject_stop = (subject_stop + s_length) as i32;
    }

    if !found_start {
        gap_align.query_start = query_offset;
        gap_align.subject_start = subject_offset;
    }
    if !found_end {
        gap_align.query_stop = query_offset + query_length;
        gap_align.subject_stop = subject_offset + subject_length;
    }
    gap_align.score = score_right + score_left;

    0
}

/// NCBI: PHIGetGappedScore (`phi_gapalign.c:739`) over Rust's
/// explicit PHI initial-hit view.
#[allow(clippy::too_many_arguments)]
pub fn phi_get_gapped_score(
    query: &[u8],
    pattern_info: Option<&crate::pattern::SphiQueryInfo>,
    subject: &[u8],
    gap_align: Option<&mut crate::blast_kappa::BlastGapAlignWorkspace>,
    matrix: &[Vec<i32>],
    gap_open: i32,
    gap_extend: i32,
    x_dropoff: i32,
    cutoff_score_min: i32,
    phi_hits: &[PhiInitialHit],
    hsp_num_max: i32,
    query_frame: i32,
    subject_frame: i32,
    hsp_list: &mut Option<crate::hspstream::HspList>,
    mut gapped_extensions: Option<&mut i32>,
) -> i16 {
    let Some(pattern_info) = pattern_info else {
        return -1;
    };
    let Some(gap_align) = gap_align else {
        return -1;
    };
    if hsp_list.is_none() {
        *hsp_list = Some(crate::hspstream::blast_hsp_list_new(hsp_num_max));
    }
    let Some(hsp_list) = hsp_list.as_mut() else {
        return -1;
    };
    if phi_hits.is_empty() {
        return 0;
    }

    for pattern_index in 0..pattern_info.num_patterns.max(0) as usize {
        let Some(query_pattern) = pattern_info.occurrences.get(pattern_index) else {
            break;
        };
        for init_hsp in phi_hits {
            let subject_length = init_hsp.subject_end - init_hsp.subject_start + 1;
            if subject_length < 0 {
                continue;
            }
            if let Some(extensions) = gapped_extensions.as_deref_mut() {
                *extensions += 1;
            }

            let status = s_phi_gapped_alignment(
                query,
                subject,
                Some(gap_align),
                matrix,
                gap_open,
                gap_extend,
                x_dropoff,
                query_pattern.offset,
                init_hsp.subject_start,
                query_pattern.length,
                subject_length,
            );
            if status != 0 {
                return status;
            }

            if gap_align.score >= cutoff_score_min {
                let hsp = crate::hspstream::Hsp {
                    score: gap_align.score,
                    num_ident: 0,
                    bit_score: 0.0,
                    evalue: f64::MAX,
                    query_offset: gap_align.query_start,
                    query_end: gap_align.query_stop,
                    query_gapped_start: query_pattern.offset,
                    subject_offset: gap_align.subject_start,
                    subject_end: gap_align.subject_stop,
                    subject_gapped_start: init_hsp.subject_start,
                    context: 0,
                    query_frame,
                    subject_frame,
                    num_gaps: 0,
                    comp_adjustment_method: 0,
                    edit_script: gap_align.edit_script.clone(),
                    pat_info: Some(crate::hspstream::PhiPatInfo {
                        index: pattern_index,
                        length: subject_length,
                    }),
                    map_info: None,
                };
                let _ = crate::hspstream::blast_hsp_list_save_hsp(hsp_list, hsp);
            }
        }
    }

    let _ = crate::hspstream::blast_hsp_list_sort_by_score(Some(hsp_list));
    0
}

/// blast-rs: bridge PHI pattern segments to the shared banded-align script API.
fn align_phi_pattern_segment_blast_rs(
    query_segment: &[u8],
    subject_segment: &[u8],
    align_script: &mut crate::gapinfo::GapPrelimEditBlock,
    matrix: &[Vec<i32>],
    gap_open: i32,
    gap_extend: i32,
) -> i32 {
    if query_segment.is_empty() && subject_segment.is_empty() {
        return 0;
    }
    if query_segment.len() == subject_segment.len() {
        crate::gapinfo::gap_prelim_edit_block_add(
            align_script,
            crate::gapinfo::GapAlignOpType::Sub,
            query_segment.len() as i32,
        );
        return 0;
    }

    let mut query_prefixed = Vec::with_capacity(query_segment.len() + 1);
    query_prefixed.push(0);
    query_prefixed.extend_from_slice(query_segment);
    let mut subject_prefixed = Vec::with_capacity(subject_segment.len() + 1);
    subject_prefixed.push(0);
    subject_prefixed.extend_from_slice(subject_segment);
    crate::gapinfo::s_banded_align(
        &query_prefixed,
        &subject_prefixed,
        query_segment.len() as i32,
        subject_segment.len() as i32,
        -5,
        5,
        matrix,
        gap_open,
        gap_extend,
        align_script,
    )
}

/// NCBI: s_PHIBlastAlignPatterns (`phi_gapalign.c:497`).
#[allow(clippy::too_many_arguments)]
pub fn s_phi_blast_align_patterns(
    query_seq: &[u8],
    db_seq: &[u8],
    len_query_seq: i32,
    len_db_seq: i32,
    align_script: &mut crate::gapinfo::GapPrelimEditBlock,
    gap_open: i32,
    gap_extend: i32,
    matrix: &[Vec<i32>],
    pattern_blk: &crate::pattern::PhiPatternSearchBlk,
) -> i32 {
    let mut start_query_match = 0;
    let mut end_query_match = len_query_seq - 1;
    let mut start_db_match = 0;
    let mut end_db_match = len_db_seq - 1;

    match pattern_blk.flag_pattern_length {
        crate::pattern::PatternType::OneWord => {
            let _ = crate::pattern::s_phi_get_short_pattern(
                query_seq,
                len_query_seq,
                &mut start_query_match,
                &mut end_query_match,
                pattern_blk,
            );
            let _ = crate::pattern::s_phi_get_short_pattern(
                db_seq,
                len_db_seq,
                &mut start_db_match,
                &mut end_db_match,
                pattern_blk,
            );
        }
        crate::pattern::PatternType::MultiWord => {
            crate::pattern::s_phi_get_long_pattern(
                query_seq,
                len_query_seq,
                &mut start_query_match,
                &mut end_query_match,
                pattern_blk,
            );
            crate::pattern::s_phi_get_long_pattern(
                db_seq,
                len_db_seq,
                &mut start_db_match,
                &mut end_db_match,
                pattern_blk,
            );
        }
        crate::pattern::PatternType::VeryLong => {
            let mut query_offsets = Vec::new();
            let mut db_offsets = Vec::new();
            let _ = crate::pattern::s_phi_get_extra_long_pattern(
                query_seq,
                len_query_seq,
                &mut query_offsets,
                pattern_blk,
            );
            let _ = crate::pattern::s_phi_get_extra_long_pattern(
                db_seq,
                len_db_seq,
                &mut db_offsets,
                pattern_blk,
            );
            let mut score = 0;
            let mut query_cursor = 0usize;
            let mut db_cursor = 0usize;
            for (&query_end, &db_end) in query_offsets.iter().zip(db_offsets.iter()) {
                let query_end = query_end.max(0) as usize;
                let db_end = db_end.max(0) as usize;
                let query_len = query_end.saturating_sub(query_cursor);
                let db_len = db_end.saturating_sub(db_cursor);
                let query_slice = query_seq
                    .get(query_cursor..query_cursor + query_len)
                    .unwrap_or(&[]);
                let db_slice = db_seq.get(db_cursor..db_cursor + db_len).unwrap_or(&[]);
                score += align_phi_pattern_segment_blast_rs(
                    query_slice,
                    db_slice,
                    align_script,
                    matrix,
                    gap_open,
                    gap_extend,
                );
                query_cursor = query_end;
                db_cursor = db_end;
            }
            return score;
        }
    }

    let prefix_query = start_query_match.max(0) as usize;
    let prefix_db = start_db_match.max(0) as usize;
    let mut local_score = align_phi_pattern_segment_blast_rs(
        query_seq.get(..prefix_query).unwrap_or(&[]),
        db_seq.get(..prefix_db).unwrap_or(&[]),
        align_script,
        matrix,
        gap_open,
        gap_extend,
    );

    let query_match_len = (end_query_match - start_query_match + 1).max(0) as usize;
    let db_match_len = (end_db_match - start_db_match + 1).max(0) as usize;
    local_score += align_phi_pattern_segment_blast_rs(
        query_seq
            .get(prefix_query..prefix_query + query_match_len)
            .unwrap_or(&[]),
        db_seq
            .get(prefix_db..prefix_db + db_match_len)
            .unwrap_or(&[]),
        align_script,
        matrix,
        gap_open,
        gap_extend,
    );

    let query_suffix_start = (end_query_match + 1).max(0) as usize;
    let db_suffix_start = (end_db_match + 1).max(0) as usize;
    local_score += align_phi_pattern_segment_blast_rs(
        query_seq.get(query_suffix_start..).unwrap_or(&[]),
        db_seq.get(db_suffix_start..).unwrap_or(&[]),
        align_script,
        matrix,
        gap_open,
        gap_extend,
    );

    local_score
}

/// NCBI: PHIGappedAlignmentWithTraceback (`phi_gapalign.c:837`).
#[allow(clippy::too_many_arguments)]
pub fn phi_gapped_alignment_with_traceback(
    query: &[u8],
    subject: &[u8],
    gap_align: Option<&mut crate::blast_kappa::BlastGapAlignWorkspace>,
    matrix: &[Vec<i32>],
    gap_open: i32,
    gap_extend: i32,
    x_dropoff: i32,
    q_start: i32,
    s_start: i32,
    query_length: i32,
    subject_length: i32,
    q_pat_length: i32,
    s_pat_length: i32,
    pattern_blk: Option<&crate::pattern::PhiPatternSearchBlk>,
) -> i16 {
    let Some(gap_align) = gap_align else {
        return -1;
    };
    let Some(pattern_blk) = pattern_blk else {
        return -1;
    };
    if q_start < 0 || s_start < 0 || q_pat_length < 0 || s_pat_length < 0 {
        return -1;
    }

    let mut rev_prelim_tback = crate::gapinfo::gap_prelim_edit_block_new();
    let mut fwd_prelim_tback = crate::gapinfo::gap_prelim_edit_block_new();
    let mut pat_prelim_tback = crate::gapinfo::gap_prelim_edit_block_new();

    let q_start_usize = q_start as usize;
    let s_start_usize = s_start as usize;
    let (score_left, private_q_length, private_s_length, rev_ops) = align_ex_generic(
        query,
        subject,
        q_start_usize,
        s_start_usize,
        matrix,
        gap_open,
        gap_extend,
        x_dropoff,
        true,
    );
    gap_align.query_start = q_start - private_q_length as i32;
    gap_align.subject_start = s_start - private_s_length as i32;
    prelim_block_set_ops(&mut rev_prelim_tback, rev_ops);

    let query_pattern = query
        .get(q_start_usize..q_start_usize.saturating_add(q_pat_length as usize))
        .unwrap_or(&[]);
    let subject_pattern = subject
        .get(s_start_usize..s_start_usize.saturating_add(s_pat_length as usize))
        .unwrap_or(&[]);
    let _ = s_phi_blast_align_patterns(
        query_pattern,
        subject_pattern,
        q_pat_length,
        s_pat_length,
        &mut pat_prelim_tback,
        gap_open,
        gap_extend,
        matrix,
        pattern_blk,
    );
    crate::gapinfo::gap_prelim_edit_block_append(&mut rev_prelim_tback, &pat_prelim_tback);

    let q_right_start = q_start + q_pat_length - 1;
    let s_right_start = s_start + s_pat_length - 1;
    let mut score_right = 0;
    if q_right_start < query_length && s_right_start < subject_length {
        let q_right_usize = q_right_start.max(0) as usize;
        let s_right_usize = s_right_start.max(0) as usize;
        let (score, private_q_right, private_s_right, fwd_ops) = align_ex_generic(
            query.get(q_right_usize..).unwrap_or(&[]),
            subject.get(s_right_usize..).unwrap_or(&[]),
            (query_length - q_right_start - 1).max(0) as usize,
            (subject_length - s_right_start - 1).max(0) as usize,
            matrix,
            gap_open,
            gap_extend,
            x_dropoff,
            false,
        );
        score_right = score;
        gap_align.query_stop = q_right_start + private_q_right as i32 + 1;
        gap_align.subject_stop = s_right_start + private_s_right as i32 + 1;
        prelim_block_set_ops(&mut fwd_prelim_tback, fwd_ops);
    } else {
        gap_align.query_stop = q_right_start;
        gap_align.subject_stop = s_right_start;
    }

    gap_align.edit_script = crate::gapinfo::blast_prelim_edit_block_to_gap_edit_script(
        Some(&rev_prelim_tback),
        Some(&fwd_prelim_tback),
    );
    gap_align.score = score_right + score_left;
    0
}

/// NCBI: s_PHITracebackFromHSPList (`blast_traceback.c:752`).
#[allow(clippy::too_many_arguments)]
pub fn s_phi_traceback_from_hsp_list(
    program_number: crate::program::ProgramType,
    hsp_list: &mut crate::hspstream::HspList,
    query: &[u8],
    subject: &[u8],
    gap_align: &mut crate::blast_kappa::BlastGapAlignWorkspace,
    matrix: &[Vec<i32>],
    scoring_params: &crate::parameters::ScoringParameters,
    hit_params: &crate::parameters::HitSavingParameters,
    pattern_info: Option<&crate::pattern::SphiQueryInfo>,
    pattern_blk: Option<&crate::pattern::PhiPatternSearchBlk>,
    lambda: f64,
    param_c: f64,
) -> i16 {
    if !crate::program::blast_program_is_phi_blast(program_number) {
        return -1;
    }
    let Some(pattern_info) = pattern_info else {
        return -1;
    };
    let Some(pattern_blk) = pattern_blk else {
        return -1;
    };
    if hsp_list.hsps.is_empty() {
        return 0;
    }

    let mut retained = Vec::with_capacity(hsp_list.hsps.len());
    for mut hsp in hsp_list.hsps.drain(..) {
        let Some(pat_info) = hsp.pat_info else {
            continue;
        };
        let Some(query_pattern) = pattern_info.occurrences.get(pat_info.index) else {
            continue;
        };
        let gap_x_dropoff = gap_align.gap_x_dropoff;
        let status = phi_gapped_alignment_with_traceback(
            query,
            subject,
            Some(gap_align),
            matrix,
            scoring_params.gap_open,
            scoring_params.gap_extend,
            gap_x_dropoff,
            hsp.query_gapped_start,
            hsp.subject_gapped_start,
            query.len() as i32,
            subject.len() as i32,
            query_pattern.length,
            pat_info.length,
            Some(pattern_blk),
        );
        if status != 0 {
            return status;
        }
        if gap_align.score >= hit_params.cutoff_score_min {
            let _ =
                crate::hspstream::blast_hsp_update_with_traceback(Some(gap_align), Some(&mut hsp));
            retained.push(hsp);
        } else {
            gap_align.edit_script =
                crate::gapinfo::gap_edit_script_delete(gap_align.edit_script.take());
        }
    }
    hsp_list.hsps = retained;

    let _ = crate::hspstream::blast_hsp_list_sort_by_score(Some(hsp_list));
    let _ = crate::hspstream::blast_hsp_list_purge_null_hsps(Some(hsp_list));
    let occurrence_offsets: Vec<i32> = pattern_info
        .occurrences
        .iter()
        .map(|occurrence| occurrence.offset)
        .collect();
    crate::hspstream::blast_hsp_list_phi_get_evalues(
        hsp_list,
        lambda,
        param_c,
        &occurrence_offsets,
        pattern_blk.min_pattern_match_length,
        pattern_blk.num_patterns_db,
    );
    let _ = crate::hspstream::blast_hsp_list_reap_by_evalue(Some(hsp_list), &hit_params.options);
    crate::hspstream::blast_hsp_list_phi_get_bit_scores(hsp_list, lambda, param_c);
    0
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
        let subj_base = crate::encoding::ncbi2na_base_at(
            subject_packed,
            subject_base_offset + if reverse { m - ai } else { ai - 1 },
        );
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
            let sgc = sa[bi].best_gap;
            let next_sc = if q_idx < query.len() {
                sa[bi].best + mrow[query[q_idx] as usize & 0x0F]
            } else {
                MININT
            };

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
pub fn blast_hsp_reevaluate_with_ambiguities_gapped(
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
    let matrix = blast_score_blk_nucl_matrix_create(reward, penalty);
    if tb.edit_script.ops.is_empty() {
        return true;
    }
    let (factor, effective_gap_open, effective_gap_extend) = if gap_open == 0 && gap_extend == 0 {
        // NCBI `Blast_HSPReevaluateWithAmbiguitiesGapped` (blast_hits.c:517): for
        // non-affine greedy, factor=2 ONLY when reward is odd (else factor=1), and
        // gap_extend = (reward - 2*penalty) * factor / 2.
        let factor = if reward % 2 == 1 { 2 } else { 1 };
        (factor, 0, (reward - 2 * penalty) * factor / 2)
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
    // NCBI seeds with -1 (`blast_hits.c:538`); only set to a non-negative value
    // inside the `sum > score` branch. The post-loop right-extension does
    // `best_end_esp_num += ext`, so an initial -1 yields `ext - 1` in the
    // edge case where the loop never enters the better-score branch.
    let mut best_end_esp_num: i32 = -1;

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

    tb.edit_script.ops = ops;
    s_update_reevaluated_hsp(
        tb,
        true,
        cutoff_score,
        score,
        best_q_start,
        best_q_end,
        best_s_start,
        best_s_end,
        best_start_esp_index,
        best_end_esp_index,
        best_end_esp_num,
    )
}

/// Port of NCBI internal `s_UpdateReevaluatedHSP` (`blast_hits.c:440`).
///
/// Rust stores absolute offsets in `TracebackResult`; the C function receives
/// raw sequence pointers and computes these same offsets by pointer
/// subtraction.
#[allow(clippy::too_many_arguments)]
pub fn s_update_reevaluated_hsp(
    tb: &mut TracebackResult,
    gapped: bool,
    cutoff_score: i32,
    score: i32,
    best_q_start: usize,
    best_q_end: usize,
    best_s_start: usize,
    best_s_end: usize,
    best_start_esp_index: usize,
    best_end_esp_index: usize,
    best_end_esp_num: i32,
) -> bool {
    tb.score = score;

    if tb.score < cutoff_score {
        return true;
    }

    tb.query_start = best_q_start;
    tb.query_end = best_q_start + best_q_end.saturating_sub(best_q_start);
    tb.subject_start = best_s_start;
    tb.subject_end = best_s_start + best_s_end.saturating_sub(best_s_start);

    if gapped {
        let last_num = match tb.edit_script.ops.len().checked_sub(1) {
            Some(idx) => idx,
            None => return false,
        };

        if best_end_esp_index >= tb.edit_script.ops.len()
            || best_start_esp_index > best_end_esp_index
        {
            return false;
        }

        if best_end_esp_index != last_num || best_start_esp_index > 0 {
            tb.edit_script.ops = tb.edit_script.ops[best_start_esp_index..=best_end_esp_index]
                .iter()
                .copied()
                .collect();
        }

        if let Some(last) = tb.edit_script.ops.last_mut() {
            last.1 = best_end_esp_num;
        }
    }

    false
}

/// Port of NCBI internal `s_UpdateReevaluatedHSPUngapped` (`blast_hits.c:664`).
#[allow(clippy::too_many_arguments)]
pub fn s_update_reevaluated_hsp_ungapped(
    tb: &mut TracebackResult,
    cutoff_score: i32,
    score: i32,
    best_q_start: usize,
    best_q_end: usize,
    best_s_start: usize,
    best_s_end: usize,
) -> bool {
    s_update_reevaluated_hsp(
        tb,
        false,
        cutoff_score,
        score,
        best_q_start,
        best_q_end,
        best_s_start,
        best_s_end,
        0,
        0,
        0,
    )
}

pub fn blast_gapped_alignment_with_traceback(
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
    let matrix = blast_score_blk_nucl_matrix_create(reward, penalty);

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

    if score_l + score_r <= 0 {
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

    let mut query_start = q_start;
    let mut subject_start = s_start;
    let mut query_end = seed_q + qr + 1;
    let mut subject_end = seed_s + sr + 1;
    let mut score_left = score_l;
    let mut score_right = score_r;

    // NCBI `BLAST_GappedAlignmentWithTraceback` (`blast_gapalign.c:4771-4801`):
    // rarely (typically when the scoring system changes between the score-only
    // and traceback stages, as happens with composition-based statistics) it is
    // possible to compute an optimal alignment with a leading or trailing gap.
    // Prune these unneeded gaps here and update the score and alignment
    // boundaries. The penalty add-back uses positive-magnitude gap_open/gap_extend
    // (the DP applied them as negatives), so removing the op cancels them out.
    // GapAlignOpType convention: Del == gap in query == consumes subject;
    // Ins == gap in subject == consumes query (matches C eGapAlignDel/eGapAlignIns).
    while !esp.ops.is_empty() && esp.ops[0].0 != GapAlignOpType::Sub {
        let (op, num) = esp.ops[0];
        score_left += gap_open + num * gap_extend;
        if op == GapAlignOpType::Del {
            subject_start += num as usize;
        } else {
            query_start += num as usize;
        }
        esp.ops.remove(0);
    }
    while !esp.ops.is_empty() && esp.ops[esp.ops.len() - 1].0 != GapAlignOpType::Sub {
        let (op, num) = esp.ops[esp.ops.len() - 1];
        score_right += gap_open + num * gap_extend;
        if op == GapAlignOpType::Del {
            subject_end -= num as usize;
        } else {
            query_end -= num as usize;
        }
        esp.ops.pop();
    }

    let total_score = score_left + score_right;

    Some(TracebackResult {
        score: total_score,
        edit_script: esp,
        query_start,
        query_end,
        subject_start,
        subject_end,
    })
}

fn oof_script_to_op(script: u8) -> GapAlignOpType {
    match script & SCRIPT_OP_MASK {
        SCRIPT_GAP_IN_A => GapAlignOpType::Del,
        SCRIPT_AHEAD_ONE_FRAME => GapAlignOpType::Del2,
        SCRIPT_AHEAD_TWO_FRAMES => GapAlignOpType::Del1,
        SCRIPT_SUB => GapAlignOpType::Sub,
        SCRIPT_NEXT_PLUS_ONE_FRAME => GapAlignOpType::Ins1,
        SCRIPT_NEXT_PLUS_TWO_FRAMES => GapAlignOpType::Ins2,
        SCRIPT_GAP_IN_B => GapAlignOpType::Ins,
        _ => GapAlignOpType::Sub,
    }
}

fn oof_initial_trace_script(
    frame: u8,
    score_other_frame1: i32,
    score_other_frame2: i32,
    score_col1: i32,
    score_col2: i32,
    score_col3: i32,
    shift_penalty: i32,
) -> u8 {
    let score_col = match frame {
        0 => score_col1,
        1 => score_col2,
        _ => score_col3,
    };
    let shifted = score_other_frame1.max(score_other_frame2) - shift_penalty;
    if score_col >= shifted {
        return SCRIPT_NEXT_IN_FRAME;
    }

    match frame {
        0 if shifted + shift_penalty == score_other_frame1 => {
            if score_other_frame1 == score_col2 {
                SCRIPT_AHEAD_TWO_FRAMES
            } else {
                SCRIPT_NEXT_PLUS_TWO_FRAMES
            }
        }
        0 => {
            if score_other_frame2 == score_col3 {
                SCRIPT_AHEAD_ONE_FRAME
            } else {
                SCRIPT_NEXT_PLUS_ONE_FRAME
            }
        }
        1 if shifted + shift_penalty == score_other_frame1 => {
            if score_other_frame1 == score_col1 {
                SCRIPT_AHEAD_ONE_FRAME
            } else {
                SCRIPT_NEXT_PLUS_ONE_FRAME
            }
        }
        1 => {
            if score_other_frame2 == score_col3 {
                SCRIPT_AHEAD_TWO_FRAMES
            } else {
                SCRIPT_NEXT_PLUS_TWO_FRAMES
            }
        }
        _ if shifted + shift_penalty == score_other_frame1 => {
            if score_other_frame1 == score_col1 {
                SCRIPT_AHEAD_TWO_FRAMES
            } else {
                SCRIPT_NEXT_PLUS_TWO_FRAMES
            }
        }
        _ => {
            if score_other_frame2 == score_col2 {
                SCRIPT_AHEAD_ONE_FRAME
            } else {
                SCRIPT_NEXT_PLUS_ONE_FRAME
            }
        }
    }
}

/// Port-shaped low-level OOF traceback entry for NCBI
/// `s_OutOfFrameAlignWithTraceback` (`blast_gapalign.c:1334`).
#[allow(clippy::too_many_arguments)]
pub fn s_out_of_frame_align_with_traceback(
    a: &[u8],
    b: &[u8],
    m: i32,
    n: i32,
    a_offset: &mut i32,
    b_offset: &mut i32,
    edit_block: Option<&mut GapPrelimEditBlock>,
    scoring: &OutOfFrameScoring,
    query_offset: i32,
    reversed: bool,
) -> i32 {
    *a_offset = 0;
    *b_offset = -2;
    if n <= 0 || m <= 0 {
        return 0;
    }

    let m = m.min(a.len() as i32).max(0) as usize;
    let n = n.min(b.len() as i32).max(0) as usize;
    if m == 0 || n == 0 {
        return 0;
    }

    let minint = i32::MIN / 4;
    let gap_open_extend = scoring.gap_open + scoring.gap_extend;
    let x_dropoff = scoring.x_dropoff.max(gap_open_extend);
    let num_extra_cells = if scoring.gap_extend > 0 {
        3 * (x_dropoff / scoring.gap_extend + 5)
    } else {
        n as i32 + 5
    }
    .max(0) as usize;

    let mut score_array = vec![
        OofScoreCell {
            best: minint,
            best_gap: minint,
        };
        n + num_extra_cells + 8
    ];
    let mut edit_script: Vec<Vec<u8>> = Vec::with_capacity(m + 1);
    let mut edit_start_offset: Vec<usize> = Vec::with_capacity(m + 1);
    edit_script.push(vec![SCRIPT_GAP_IN_B; score_array.len()]);
    edit_start_offset.push(0);

    let mut score = -gap_open_extend;
    score_array[0].best = 0;
    score_array[0].best_gap = -gap_open_extend;
    let mut i = 3usize;
    while i <= n + 2 {
        if i + 2 >= score_array.len() {
            score_array.resize(
                i + 6,
                OofScoreCell {
                    best: minint,
                    best_gap: minint,
                },
            );
            edit_script[0].resize(i + 6, SCRIPT_GAP_IN_B);
        }
        score_array[i].best = score;
        score_array[i].best_gap = score - gap_open_extend;
        edit_script[0][i] = SCRIPT_GAP_IN_B;
        score_array[i - 1].best = minint;
        score_array[i - 1].best_gap = minint;
        edit_script[0][i - 1] = SCRIPT_GAP_IN_B;
        score_array[i - 2].best = minint;
        score_array[i - 2].best_gap = minint;
        edit_script[0][i - 2] = SCRIPT_GAP_IN_B;
        if score < -x_dropoff {
            break;
        }
        score -= scoring.gap_extend;
        i += 3;
    }

    let mut b_size = i.saturating_sub(2);
    if i >= score_array.len() {
        score_array.resize(
            i + 1,
            OofScoreCell {
                best: minint,
                best_gap: minint,
            },
        );
        edit_script[0].resize(i + 1, SCRIPT_GAP_IN_B);
    }
    score_array[i].best = minint;
    score_array[i].best_gap = minint;

    let mut best_score = 0;
    let mut first_b_index = 0usize;
    for a_index in 1..=m {
        let mut row = vec![0u8; b_size.saturating_sub(first_b_index) + num_extra_cells + 8];
        let row_offset = first_b_index;

        let mut score_row1 = minint;
        let mut score_row2 = minint;
        let mut score_row3 = minint;
        let mut score_col1 = minint;
        let mut score_col2 = minint;
        let mut score_col3 = minint;
        let mut score_other_frame1 = minint;
        let mut score_other_frame2 = minint;
        let mut last_b_index = first_b_index;
        let mut b_index = first_b_index;

        while b_index < b_size {
            let base = b_index - row_offset;
            if base >= row.len() {
                row.resize(base + 1, 0);
            }

            let frame_score = scoring.oof_row_score(
                query_offset,
                a,
                m,
                a_index,
                oof_subject_letter(b, b_index, reversed),
                reversed,
            );
            let script = oof_initial_trace_script(
                0,
                score_other_frame1,
                score_other_frame2,
                score_col1,
                score_col2,
                score_col3,
                scoring.shift_penalty,
            );
            let script = oof_step_frame_with_script(
                &mut score_array,
                b_index,
                &mut score_row1,
                &mut score_col1,
                &mut score_other_frame1,
                &mut score_other_frame2,
                &mut best_score,
                a_offset,
                b_offset,
                a_index,
                frame_score,
                gap_open_extend,
                scoring.gap_extend,
                scoring.shift_penalty,
                x_dropoff,
                &mut first_b_index,
                &mut last_b_index,
                script,
                false,
            );
            row[base] = script;
            b_index += 1;
            if b_index >= b_size {
                let score = score_row1;
                score_row1 = score_row2;
                score_row2 = score_row3;
                score_row3 = score;
                let _ = (score_row1, score_row2, score_row3);
                break;
            }

            let base = b_index - row_offset;
            if base >= row.len() {
                row.resize(base + 1, 0);
            }
            let frame_score = scoring.oof_row_score(
                query_offset,
                a,
                m,
                a_index,
                oof_subject_letter(b, b_index, reversed),
                reversed,
            );
            let script = oof_initial_trace_script(
                1,
                score_other_frame1,
                score_other_frame2,
                score_col1,
                score_col2,
                score_col3,
                scoring.shift_penalty,
            );
            let script = oof_step_frame_with_script(
                &mut score_array,
                b_index,
                &mut score_row2,
                &mut score_col2,
                &mut score_other_frame2,
                &mut score_other_frame1,
                &mut best_score,
                a_offset,
                b_offset,
                a_index,
                frame_score,
                gap_open_extend,
                scoring.gap_extend,
                scoring.shift_penalty,
                x_dropoff,
                &mut first_b_index,
                &mut last_b_index,
                script,
                true,
            );
            row[base] = script;
            b_index += 1;
            if b_index >= b_size {
                let score = score_row2;
                score_row2 = score_row1;
                score_row1 = score_row3;
                score_row3 = score;
                let _ = (score_row1, score_row2, score_row3);
                break;
            }

            let base = b_index - row_offset;
            if base >= row.len() {
                row.resize(base + 1, 0);
            }
            let frame_score = scoring.oof_row_score(
                query_offset,
                a,
                m,
                a_index,
                oof_subject_letter(b, b_index, reversed),
                reversed,
            );
            let frame2_old_other2 = score_other_frame2;
            let script = oof_initial_trace_script(
                2,
                score_other_frame1,
                score_other_frame2,
                score_col1,
                score_col2,
                score_col3,
                scoring.shift_penalty,
            );
            let script = oof_step_frame_with_script(
                &mut score_array,
                b_index,
                &mut score_row3,
                &mut score_col3,
                &mut score_other_frame1,
                &mut score_other_frame2,
                &mut best_score,
                a_offset,
                b_offset,
                a_index,
                frame_score,
                gap_open_extend,
                scoring.gap_extend,
                scoring.shift_penalty,
                x_dropoff,
                &mut first_b_index,
                &mut last_b_index,
                script,
                true,
            );
            score_other_frame1 = frame2_old_other2;
            row[base] = script;
            b_index += 1;
        }

        if first_b_index == b_size {
            edit_script.push(row);
            edit_start_offset.push(row_offset);
            break;
        }

        if b_size + num_extra_cells + 5 >= score_array.len() {
            score_array.resize(
                b_size + num_extra_cells + 100,
                OofScoreCell {
                    best: minint,
                    best_gap: minint,
                },
            );
        }

        if last_b_index < b_size.saturating_sub(1) {
            b_size = last_b_index + 1;
        } else {
            let mut score = score_row1.max(score_row2).max(score_row3);
            while score >= best_score - x_dropoff && b_size < n + 1 {
                let base = b_size - row_offset;
                if base + 2 >= row.len() {
                    row.resize(base + 3, 0);
                }
                score_array[b_size].best = score_row1;
                score_array[b_size].best_gap = score_row1 - gap_open_extend;
                score_row1 -= scoring.gap_extend;
                row[base] = SCRIPT_OOF_OPEN_GAP | SCRIPT_GAP_IN_B;

                score_array[b_size + 1].best = score_row2;
                score_array[b_size + 1].best_gap = score_row2 - gap_open_extend;
                score_row2 -= scoring.gap_extend;
                row[base + 1] = SCRIPT_OOF_OPEN_GAP | SCRIPT_GAP_IN_B;

                score_array[b_size + 2].best = score_row3;
                score_array[b_size + 2].best_gap = score_row3 - gap_open_extend;
                score_row3 -= scoring.gap_extend;
                row[base + 2] = SCRIPT_OOF_OPEN_GAP | SCRIPT_GAP_IN_B;

                b_size += 3;
                score -= scoring.gap_extend;
            }
        }

        b_size = b_size.min(n + 1);
        let used = b_index.max(b_size).saturating_sub(row_offset) + 1;
        if used > row.len() {
            row.resize(used, 0);
        }
        edit_script.push(row);
        edit_start_offset.push(row_offset);

        let last = (b_size + 4).min(n + 3);
        while b_size < last {
            if b_size >= score_array.len() {
                score_array.resize(
                    b_size + 1,
                    OofScoreCell {
                        best: minint,
                        best_gap: minint,
                    },
                );
            }
            score_array[b_size].best = minint;
            score_array[b_size].best_gap = minint;
            b_size += 1;
        }
    }

    if let Some(edit_block) = edit_block {
        crate::gapinfo::gap_prelim_edit_block_reset(Some(edit_block));
        let mut a_index = (*a_offset).max(0) as usize;
        let mut b_index = (*b_offset).max(0) as usize;
        let mut script = 1u8;

        while a_index > 0 || b_index > 0 {
            let row_offset = edit_start_offset.get(a_index).copied().unwrap_or(0);
            let next_script = edit_script
                .get(a_index)
                .and_then(|row| row.get(b_index.saturating_sub(row_offset)))
                .copied()
                .unwrap_or(SCRIPT_SUB);

            script = match script {
                SCRIPT_GAP_IN_A => {
                    let next_op = next_script & SCRIPT_OP_MASK;
                    if next_script & (SCRIPT_OOF_OPEN_GAP | SCRIPT_EXTEND_GAP_A) != 0 {
                        SCRIPT_GAP_IN_A
                    } else {
                        next_op
                    }
                }
                SCRIPT_GAP_IN_B => {
                    let next_op = next_script & SCRIPT_OP_MASK;
                    if next_script & (SCRIPT_OOF_OPEN_GAP | SCRIPT_EXTEND_GAP_B) != 0 {
                        SCRIPT_GAP_IN_B
                    } else {
                        next_op
                    }
                }
                _ => next_script & SCRIPT_OP_MASK,
            };

            if script == SCRIPT_GAP_IN_B {
                b_index = b_index.saturating_sub(3);
            } else {
                b_index = b_index.saturating_sub(script as usize);
                a_index = a_index.saturating_sub(1);
            }
            crate::gapinfo::gap_prelim_edit_block_add(edit_block, oof_script_to_op(script), 1);
        }
    }

    if !reversed {
        *b_offset -= 2;
    }
    best_score
}

#[derive(Clone, Copy)]
struct OofScoreCell {
    best: i32,
    best_gap: i32,
}

fn oof_subject_letter(b: &[u8], b_index: usize, reversed: bool) -> u8 {
    if reversed {
        b.get(b.len().saturating_sub(1).saturating_sub(b_index))
            .copied()
            .unwrap_or(0)
    } else {
        b_index
            .checked_sub(2)
            .and_then(|idx| b.get(idx))
            .copied()
            .unwrap_or(0)
    }
}

#[allow(clippy::too_many_arguments)]
fn oof_step_frame_with_script(
    score_array: &mut [OofScoreCell],
    b_index: usize,
    score_row: &mut i32,
    score_col: &mut i32,
    score_other_a: &mut i32,
    score_other_b: &mut i32,
    best_score: &mut i32,
    a_offset: &mut i32,
    b_offset: &mut i32,
    a_index: usize,
    frame_score: i32,
    gap_open_extend: i32,
    gap_extend: i32,
    shift_penalty: i32,
    x_dropoff: i32,
    first_b_index: &mut usize,
    last_b_index: &mut usize,
    initial_script: u8,
    gap_tie_prefers_column: bool,
) -> u8 {
    let minint = i32::MIN / 4;
    let mut script = initial_script;
    let mut score = (*score_other_a).max(*score_other_b) - shift_penalty;
    score = score.max(*score_col) + frame_score;
    *score_other_a = (*score_col).max(score_array[b_index].best);
    *score_col = score_array[b_index].best;
    let mut score_gap_col = score_array[b_index].best_gap;

    if score < score_gap_col.max(*score_row) {
        if score_gap_col > *score_row || (gap_tie_prefers_column && score_gap_col == *score_row) {
            score = score_gap_col;
            script = SCRIPT_OOF_OPEN_GAP | SCRIPT_GAP_IN_A;
        } else {
            score = *score_row;
            script = SCRIPT_OOF_OPEN_GAP | SCRIPT_GAP_IN_B;
        }

        if *best_score - score > x_dropoff {
            if *first_b_index == b_index {
                *first_b_index = b_index + 1;
            } else {
                score_array[b_index].best = minint;
            }
        } else {
            *last_b_index = b_index;
            score_array[b_index].best = score;
            score_array[b_index].best_gap = score_gap_col - gap_extend;
            *score_row -= gap_extend;
        }
    } else if *best_score - score > x_dropoff {
        if *first_b_index == b_index {
            *first_b_index = b_index + 1;
        } else {
            score_array[b_index].best = minint;
        }
    } else {
        *last_b_index = b_index;
        score_array[b_index].best = score;
        if score > *best_score {
            *best_score = score;
            *a_offset = a_index as i32;
            *b_offset = b_index as i32;
        }

        score -= gap_open_extend;
        *score_row -= gap_extend;
        if score > *score_row {
            *score_row = score;
        } else {
            script |= SCRIPT_EXTEND_GAP_B;
        }

        score_gap_col -= gap_extend;
        if score < score_gap_col {
            score_array[b_index].best_gap = score_gap_col;
            script |= SCRIPT_EXTEND_GAP_A;
        } else {
            score_array[b_index].best_gap = score;
        }
    }

    script
}

/// Port-shaped low-level OOF score/traceback dispatcher for NCBI
/// `s_OutOfFrameGappedAlign` (`blast_gapalign.c:1950`).
#[allow(clippy::too_many_arguments)]
pub fn s_out_of_frame_gapped_align(
    a: &[u8],
    b: &[u8],
    m: i32,
    n: i32,
    a_offset: &mut i32,
    b_offset: &mut i32,
    score_only: bool,
    edit_block: Option<&mut GapPrelimEditBlock>,
    scoring: &OutOfFrameScoring,
    query_offset: i32,
    reversed: bool,
) -> i32 {
    if !score_only {
        return s_out_of_frame_align_with_traceback(
            a,
            b,
            m,
            n,
            a_offset,
            b_offset,
            edit_block,
            scoring,
            query_offset,
            reversed,
        );
    }

    *a_offset = 0;
    *b_offset = -2;
    if n <= 0 || m <= 0 {
        return 0;
    }

    let m = m.min(a.len() as i32).max(0) as usize;
    let n = n.min(b.len() as i32).max(0) as usize;
    if m == 0 || n == 0 {
        return 0;
    }

    let minint = i32::MIN / 4;
    let gap_open_extend = scoring.gap_open + scoring.gap_extend;
    let x_dropoff = scoring.x_dropoff.max(gap_open_extend);
    let num_extra_cells = if scoring.gap_extend > 0 {
        3 * (x_dropoff / scoring.gap_extend + 5)
    } else {
        n as i32 + 5
    }
    .max(0) as usize;
    let mut score_array = vec![
        OofScoreCell {
            best: minint,
            best_gap: minint,
        };
        n + num_extra_cells + 8
    ];

    let mut score = -gap_open_extend;
    score_array[0].best = 0;
    score_array[0].best_gap = -gap_open_extend;
    let mut i = 3usize;
    while i <= n + 2 {
        if i + 2 >= score_array.len() {
            score_array.resize(
                i + 6,
                OofScoreCell {
                    best: minint,
                    best_gap: minint,
                },
            );
        }
        score_array[i].best = score;
        score_array[i].best_gap = score - gap_open_extend;
        score_array[i - 1].best = minint;
        score_array[i - 1].best_gap = minint;
        score_array[i - 2].best = minint;
        score_array[i - 2].best_gap = minint;
        if score < -x_dropoff {
            break;
        }
        score -= scoring.gap_extend;
        i += 3;
    }

    let mut b_size = i.saturating_sub(2);
    if i >= score_array.len() {
        score_array.resize(
            i + 1,
            OofScoreCell {
                best: minint,
                best_gap: minint,
            },
        );
    }
    score_array[i].best = minint;
    score_array[i].best_gap = minint;

    let mut best_score = 0;
    let mut first_b_index = 0usize;
    for a_index in 1..=m {
        let mut score_row1 = minint;
        let mut score_row2 = minint;
        let mut score_row3 = minint;
        let mut score_col1 = minint;
        let mut score_col2 = minint;
        let mut score_col3 = minint;
        let mut score_other_frame1 = minint;
        let mut score_other_frame2 = minint;
        let mut last_b_index = first_b_index;
        let mut b_index = first_b_index;

        while b_index < b_size {
            let frame_score = scoring.oof_row_score(
                query_offset,
                a,
                m,
                a_index,
                oof_subject_letter(b, b_index, reversed),
                reversed,
            );
            let _ = oof_step_frame_with_script(
                &mut score_array,
                b_index,
                &mut score_row1,
                &mut score_col1,
                &mut score_other_frame1,
                &mut score_other_frame2,
                &mut best_score,
                a_offset,
                b_offset,
                a_index,
                frame_score,
                gap_open_extend,
                scoring.gap_extend,
                scoring.shift_penalty,
                x_dropoff,
                &mut first_b_index,
                &mut last_b_index,
                SCRIPT_NEXT_IN_FRAME,
                false,
            );
            b_index += 1;
            if b_index >= b_size {
                let score = score_row1;
                score_row1 = score_row2;
                score_row2 = score_row3;
                score_row3 = score;
                let _ = (score_row1, score_row2, score_row3);
                break;
            }

            let frame_score = scoring.oof_row_score(
                query_offset,
                a,
                m,
                a_index,
                oof_subject_letter(b, b_index, reversed),
                reversed,
            );
            let _ = oof_step_frame_with_script(
                &mut score_array,
                b_index,
                &mut score_row2,
                &mut score_col2,
                &mut score_other_frame2,
                &mut score_other_frame1,
                &mut best_score,
                a_offset,
                b_offset,
                a_index,
                frame_score,
                gap_open_extend,
                scoring.gap_extend,
                scoring.shift_penalty,
                x_dropoff,
                &mut first_b_index,
                &mut last_b_index,
                SCRIPT_NEXT_IN_FRAME,
                true,
            );
            b_index += 1;
            if b_index >= b_size {
                let score = score_row2;
                score_row2 = score_row1;
                score_row1 = score_row3;
                score_row3 = score;
                let _ = (score_row1, score_row2, score_row3);
                break;
            }

            let frame_score = scoring.oof_row_score(
                query_offset,
                a,
                m,
                a_index,
                oof_subject_letter(b, b_index, reversed),
                reversed,
            );
            let old_other2 = score_other_frame2;
            let _ = oof_step_frame_with_script(
                &mut score_array,
                b_index,
                &mut score_row3,
                &mut score_col3,
                &mut score_other_frame1,
                &mut score_other_frame2,
                &mut best_score,
                a_offset,
                b_offset,
                a_index,
                frame_score,
                gap_open_extend,
                scoring.gap_extend,
                scoring.shift_penalty,
                x_dropoff,
                &mut first_b_index,
                &mut last_b_index,
                SCRIPT_NEXT_IN_FRAME,
                true,
            );
            score_other_frame1 = old_other2;
            b_index += 1;
        }

        if first_b_index == b_size {
            break;
        }

        if b_size + num_extra_cells + 5 >= score_array.len() {
            score_array.resize(
                b_size + num_extra_cells + 100,
                OofScoreCell {
                    best: minint,
                    best_gap: minint,
                },
            );
        }

        if last_b_index < b_size.saturating_sub(1) {
            b_size = last_b_index + 1;
        } else {
            let mut score = score_row1.max(score_row2).max(score_row3);
            while score >= best_score - x_dropoff && b_size < n + 1 {
                if b_size + 2 >= score_array.len() {
                    score_array.resize(
                        b_size + 6,
                        OofScoreCell {
                            best: minint,
                            best_gap: minint,
                        },
                    );
                }
                score_array[b_size].best = score_row1;
                score_array[b_size].best_gap = score_row1 - gap_open_extend;
                score_row1 -= scoring.gap_extend;

                score_array[b_size + 1].best = score_row2;
                score_array[b_size + 1].best_gap = score_row2 - gap_open_extend;
                score_row2 -= scoring.gap_extend;

                score_array[b_size + 2].best = score_row3;
                score_array[b_size + 2].best_gap = score_row3 - gap_open_extend;
                score_row3 -= scoring.gap_extend;

                b_size += 3;
                score -= scoring.gap_extend;
            }
        }

        b_size = b_size.min(n + 1);
        let last = (b_size + 4).min(n + 3);
        while b_size < last {
            if b_size >= score_array.len() {
                score_array.resize(
                    b_size + 1,
                    OofScoreCell {
                        best: minint,
                        best_gap: minint,
                    },
                );
            }
            score_array[b_size].best = minint;
            score_array[b_size].best_gap = minint;
            b_size += 1;
        }
    }

    if !reversed {
        *b_offset -= 2;
    }
    best_score
}

/// Port of NCBI wrapper `s_OutOfFrameSemiGappedAlignWrap`
/// (`blast_gapalign.c:4223`): BLASTX swaps query/subject for the underlying
/// OOF DP, TBLASTN does not.
#[allow(clippy::too_many_arguments)]
pub fn s_out_of_frame_semi_gapped_align_wrap(
    query: &[u8],
    subject: &[u8],
    q_off: i32,
    s_off: i32,
    private_q_start: &mut i32,
    private_s_start: &mut i32,
    score_only: bool,
    edit_block: Option<&mut GapPrelimEditBlock>,
    scoring: &OutOfFrameScoring,
    psi_offset: i32,
    reversed: bool,
    switch_seq: bool,
) -> i32 {
    if switch_seq {
        s_out_of_frame_gapped_align(
            subject,
            query,
            s_off,
            q_off,
            private_s_start,
            private_q_start,
            score_only,
            edit_block,
            scoring,
            psi_offset,
            reversed,
        )
    } else {
        s_out_of_frame_gapped_align(
            query,
            subject,
            q_off,
            s_off,
            private_q_start,
            private_s_start,
            score_only,
            edit_block,
            scoring,
            psi_offset,
            reversed,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::encoding::pack_ncbi2na_bases;

    #[test]
    fn test_adjacent_del_ins_gap_position_matches_ncbi() {
        // Historic reproducer (adjacent_del_ins). NCBI places the 5-gap block
        // after position 35 of the query (BTOP "35G-G-G-T-A-42"). Fixed on
        // 2026-04-18 by routing `x_dropoff_final` through the `_with_xdrops`
        // wrappers.
        //
        // At this low level the output depends on which seed from which
        // ungapped HSP you pick: early diagonal-0 seeds can yield 33+5+44,
        // while the CLI's selected seed for the NCBI-shaped HSP yields
        // 35+5+42.
        //
        // Sequences from tests/integration.rs:3341-3345.
        let q = crate::encoding::encode_blastna_sequence(
            b"ACGTTGCAACGATCGTACGATTCGAGCTTAGGCTAGGGTAATCGGATCCTAGCTAGGCTAATCGATCGTAGCTAGCATCGAT",
        );
        let s = crate::encoding::encode_blastna_sequence(
            b"ACGTTGCAACGATCGTACGATTCGAGCTTAGGCTAATCGGATCCTAGCTAGGCTAATCGATCGTAGCTAGCATCGAT",
        );
        assert_eq!(q.len(), 82);
        assert_eq!(s.len(), 77);
        // Seed = 34, matching the ungapped HSP that survives CLI dedup for
        // this fixture and reproduces NCBI's BTOP.
        // x_dropoff=50 models blastn-short `xdrop_gap_final=100` bits after
        // `BlastExtensionParametersNew` converts with lambda≈1.374:
        // `(int)(100 * ln2 / 1.374) = 50`.
        let seed_q = 34;
        let seed_s = 34;
        let r = blast_gapped_alignment_with_traceback(&q, &s, seed_q, seed_s, 1, -3, 5, 2, 50)
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
        let delete = blast_hsp_reevaluate_with_ambiguities_gapped(&mut tb, &q, &s, 1, -3, 5, 2, 1);
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
        let delete = blast_hsp_reevaluate_with_ambiguities_gapped(&mut tb, &q, &s, 1, -3, 5, 2, 1);
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
        let delete = blast_hsp_reevaluate_with_ambiguities_gapped(&mut tb, &q, &s, 1, -3, 0, 0, 1);
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
        let delete = blast_hsp_reevaluate_with_ambiguities_gapped(&mut tb, &q, &s, 1, -3, 5, 2, 5);
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
        let delete = blast_hsp_reevaluate_with_ambiguities_gapped(&mut tb, &q, &s, 1, -3, 5, 2, 1);
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
        let delete = blast_hsp_reevaluate_with_ambiguities_gapped(&mut tb, &q, &s, 1, -3, 5, 2, 1);
        assert!(!delete);
        assert_eq!(tb.score, 10);
        assert_eq!(tb.edit_script.ops, vec![(GapAlignOpType::Sub, 10)]);
        assert_eq!((tb.query_start, tb.query_end), (0, 10));
    }

    #[test]
    fn test_update_reevaluated_hsp_trims_gap_script_like_c() {
        let mut tb = TracebackResult {
            score: 0,
            edit_script: GapEditScript {
                ops: vec![
                    (GapAlignOpType::Sub, 3),
                    (GapAlignOpType::Ins, 2),
                    (GapAlignOpType::Sub, 5),
                    (GapAlignOpType::Del, 1),
                ],
            },
            query_start: 2,
            query_end: 13,
            subject_start: 10,
            subject_end: 21,
        };

        let delete = s_update_reevaluated_hsp(&mut tb, true, 20, 30, 4, 11, 12, 19, 1, 2, 4);

        assert!(!delete);
        assert_eq!(tb.score, 30);
        assert_eq!((tb.query_start, tb.query_end), (4, 11));
        assert_eq!((tb.subject_start, tb.subject_end), (12, 19));
        assert_eq!(
            tb.edit_script.ops,
            vec![(GapAlignOpType::Ins, 2), (GapAlignOpType::Sub, 4)]
        );
    }

    #[test]
    fn test_update_reevaluated_hsp_ungapped_leaves_script_and_deletes_below_cutoff() {
        let mut tb = TracebackResult {
            score: 0,
            edit_script: GapEditScript {
                ops: vec![(GapAlignOpType::Sub, 6)],
            },
            query_start: 0,
            query_end: 6,
            subject_start: 0,
            subject_end: 6,
        };

        assert!(s_update_reevaluated_hsp_ungapped(
            &mut tb, 10, 7, 1, 5, 2, 6
        ));
        assert_eq!(tb.score, 7);
        assert_eq!((tb.query_start, tb.query_end), (0, 6));
        assert_eq!(tb.edit_script.ops, vec![(GapAlignOpType::Sub, 6)]);

        assert!(!s_update_reevaluated_hsp_ungapped(
            &mut tb, 10, 11, 1, 5, 2, 6
        ));
        assert_eq!(tb.score, 11);
        assert_eq!((tb.query_start, tb.query_end), (1, 5));
        assert_eq!((tb.subject_start, tb.subject_end), (2, 6));
        assert_eq!(tb.edit_script.ops, vec![(GapAlignOpType::Sub, 6)]);
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
    fn test_blast_gapped_alignment_with_traceback_basic() {
        let q = vec![0u8, 1, 2, 3, 0, 1, 2, 3]; // ACGTACGT
        let s = vec![0u8, 1, 2, 3, 0, 1, 2, 3];
        let r = blast_gapped_alignment_with_traceback(&q, &s, 4, 4, 1, -3, 5, 2, 30);
        assert!(r.is_some(), "Should find alignment");
        let r = r.unwrap();
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
    fn test_blast_gapped_alignment_with_traceback_strips_terminal_gaps() {
        // NCBI BLAST_GappedAlignmentWithTraceback (blast_gapalign.c:4771-4801)
        // prunes leading/trailing non-Sub ops; the returned script must begin
        // and end on a Sub op.
        let q = vec![0u8, 1, 2, 3, 0, 1, 2, 3];
        let s = vec![0u8, 1, 2, 3, 0, 1, 2, 3];
        let r = blast_gapped_alignment_with_traceback(&q, &s, 4, 4, 1, -3, 5, 2, 30).unwrap();
        assert!(!r.edit_script.ops.is_empty());
        assert_eq!(r.edit_script.ops.first().unwrap().0, GapAlignOpType::Sub);
        assert_eq!(r.edit_script.ops.last().unwrap().0, GapAlignOpType::Sub);
    }

    #[test]
    fn test_blast_gapped_alignment_with_traceback_exact_match_extends_to_edges() {
        let q = crate::encoding::encode_blastna_sequence(
            b"GAATCCATGCTGTGGGCCAGCAAGAGTTAAGGTGCTCATGGTTTTGAGAAAACATCTGAGGACTCTGACAGCACTCTCCCATCCTTGGTCTCCACAGTCT",
        );
        let r = blast_gapped_alignment_with_traceback(&q, &q, 50, 50, 1, -3, 5, 2, 20)
            .expect("exact match should align");
        assert_eq!(r.score, 100);
        assert_eq!((r.query_start, r.query_end), (0, 100));
        assert_eq!((r.subject_start, r.subject_end), (0, 100));
    }

    #[test]
    fn test_blast_gapped_alignment_with_traceback_single_internal_gap() {
        let q: Vec<u8> = [
            0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
        ]
        .to_vec();
        let s: Vec<u8> = [
            0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 0, 0, 0, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
        ]
        .to_vec();

        let r = blast_gapped_alignment_with_traceback(&q, &s, 12, 12, 2, -3, 3, 1, 30)
            .expect("single-gap alignment should succeed");
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
        // Extra subject bases are represented as Del in gapinfo/C script
        // convention: deletion in query, i.e. a gap in query.
        let has_del = esp.ops.iter().any(|(op, _)| *op == GapAlignOpType::Del);
        assert!(has_del, "expected Del (gap in query), ops={:?}", esp.ops);
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
        // Extra query bases are represented as Ins in gapinfo/C script
        // convention: insertion in query, i.e. a gap in subject.
        let has_ins = esp.ops.iter().any(|(op, _)| *op == GapAlignOpType::Ins);
        assert!(has_ins, "expected Ins (gap in subject), ops={:?}", esp.ops);
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
        // Count gap operations (Del = gap in query, Ins = gap in subject).
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
                    si += *count as usize;
                }
                GapAlignOpType::Ins => {
                    // Gap in subject: consume query
                    computed_score -= gap_open + gap_ext * count;
                    qi += *count as usize;
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
        let (score, _ao, _bo, _ops) = align_ex(&a, &b, 5, 5, &matrix, 5, 2, 30, true);
        assert!(score >= 8, "should get at least 4 bases * 2, got {}", score);
    }

    #[test]
    fn test_gapped_score_one_dir_packed_subject_matches_decoded_basic() {
        let query = vec![0u8, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
        let subject = vec![0u8, 1, 2, 3, 0, 1, 3, 3, 0, 1, 2, 3];
        let packed = pack_ncbi2na_bases(&subject);
        let matrix = blast_score_blk_nucl_matrix_create(1, -3);
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
        assert_eq!(left_packed.0, left_decoded);
    }

    #[test]
    fn test_gapped_score_one_dir_processes_boundary_cell_without_next_subject() {
        let mut matrix = vec![vec![-3; 16]; 16];
        for i in 0..4 {
            matrix[i][i] = 2;
        }
        let query = vec![15u8, 0];
        let subject = vec![15u8, 0];

        let (score, a_off, b_off) =
            gapped_score_one_dir_generic(&query, &subject, 1, 1, &matrix, 7, 2, 12, false);

        assert_eq!(score, 2);
        assert_eq!(a_off, 1);
        assert_eq!(b_off, 1);
    }

    #[test]
    fn test_c_named_packed_gapped_alignment_wrappers_match_existing_paths() {
        let query = vec![0u8, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
        let subject = vec![0u8, 1, 2, 3, 0, 1, 3, 3, 0, 1, 2, 3];
        let packed = pack_ncbi2na_bases(&subject);

        let extents = s_blast_dyn_prog_nt_gapped_alignment_extents(
            &query,
            &packed,
            subject.len(),
            5,
            5,
            1,
            -3,
            5,
            2,
            12,
            8,
        );
        let wrapped_extents = s_blast_dyn_prog_nt_gapped_alignment(
            &query,
            &packed,
            subject.len(),
            5,
            5,
            1,
            -3,
            5,
            2,
            12,
            8,
        );
        assert_eq!(wrapped_extents.0, extents.0);
        assert_eq!(wrapped_extents.1, extents.1);
        assert_eq!(wrapped_extents.2, extents.2);
        assert_eq!(wrapped_extents.3, extents.3);
        assert_eq!(wrapped_extents.4, extents.4);
        assert_eq!(wrapped_extents.5, extents.5);
        assert_eq!(wrapped_extents.6, extents.6);
        assert_eq!(wrapped_extents.7, extents.7);
        assert_eq!(wrapped_extents.8, extents.8);
        assert_eq!(wrapped_extents.9, extents.9);
        assert_eq!(wrapped_extents.10, extents.10);
        assert_eq!(wrapped_extents.11, extents.11);
        assert_eq!(wrapped_extents.12, extents.12);

        let matrix = blast_score_blk_nucl_matrix_create(1, -3);
        assert_eq!(
            s_blast_align_packed_nucl(&query[..6], &packed, 6, 6, 0, 1, -3, 5, 2, 12, true),
            gapped_score_one_dir_packed_subject(
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
            )
        );
    }

    #[test]
    fn restricted_gapped_align_enforces_restrict_size_gap_starts() {
        let query = vec![0u8, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
        let subject = vec![0u8, 1, 2, 3, 0, 1, 3, 3, 0, 1, 2, 3];
        let unrestricted = blast_semi_gapped_align(&query, &subject, 5, 5, 1, -3, 5, 2, 12);
        let restricted = s_restricted_gapped_align(&query, &subject, 5, 5, 1, -3, 5, 2, 12);
        assert_eq!(unrestricted, 8);
        assert_eq!(restricted, 2);
        assert!(
            restricted <= unrestricted,
            "restricted gap starts should not improve over semigapped score"
        );
    }

    fn phi_test_matrix() -> Vec<Vec<i32>> {
        let mut matrix = vec![vec![-3; 256]; 256];
        for i in 0..256 {
            matrix[i][i] = 2;
        }
        matrix
    }

    #[test]
    fn translated_phi_gapped_alignment_matches_c_boundary_shape() {
        let matrix = phi_test_matrix();
        let mut gap_align = crate::blast_kappa::blast_gap_align_struct_new(12).expect("gap align");

        assert_eq!(
            s_phi_gapped_alignment(
                b"AB",
                b"AB",
                Some(&mut gap_align),
                &matrix,
                5,
                2,
                12,
                0,
                0,
                2,
                2,
            ),
            0
        );
        assert_eq!(gap_align.score, 0);
        assert_eq!((gap_align.query_start, gap_align.subject_start), (0, 0));
        assert_eq!((gap_align.query_stop, gap_align.subject_stop), (1, 1));

        assert_eq!(
            s_phi_gapped_alignment(b"AB", b"AB", None, &matrix, 5, 2, 12, 0, 0, 2, 2,),
            -1
        );
    }

    #[test]
    fn translated_phi_gapped_alignment_extends_around_pattern() {
        let matrix = phi_test_matrix();
        let mut gap_align = crate::blast_kappa::blast_gap_align_struct_new(12).expect("gap align");

        assert_eq!(
            s_phi_gapped_alignment(
                b"XXABYY",
                b"XXABYY",
                Some(&mut gap_align),
                &matrix,
                5,
                2,
                12,
                2,
                2,
                2,
                2,
            ),
            0
        );

        assert_eq!(gap_align.score, 8);
        assert_eq!((gap_align.query_start, gap_align.subject_start), (1, 1));
        assert_eq!((gap_align.query_stop, gap_align.subject_stop), (5, 5));
    }

    #[test]
    fn translated_phi_get_gapped_score_saves_pattern_hsps() {
        let matrix = phi_test_matrix();
        let mut gap_align = crate::blast_kappa::blast_gap_align_struct_new(12).expect("gap align");
        let pattern_info = crate::pattern::SphiQueryInfo {
            num_patterns: 1,
            occurrences: vec![crate::pattern::SphiPatternInfo {
                offset: 2,
                length: 2,
            }],
            allocated_size: 1,
            probability: 0.0,
            pattern: None,
        };
        let phi_hits = [PhiInitialHit {
            subject_start: 2,
            subject_end: 3,
        }];
        let mut list = None;
        let mut extensions = 0;

        assert_eq!(
            phi_get_gapped_score(
                b"XXABYY",
                Some(&pattern_info),
                b"XXABYY",
                Some(&mut gap_align),
                &matrix,
                5,
                2,
                12,
                1,
                &phi_hits,
                10,
                1,
                0,
                &mut list,
                Some(&mut extensions),
            ),
            0
        );

        let list = list.expect("hsp list");
        assert_eq!(extensions, 1);
        assert_eq!(list.hsps.len(), 1);
        let hsp = &list.hsps[0];
        assert_eq!(hsp.score, 8);
        assert_eq!((hsp.query_offset, hsp.query_end), (1, 5));
        assert_eq!((hsp.subject_offset, hsp.subject_end), (1, 5));
        assert_eq!(
            hsp.pat_info,
            Some(crate::hspstream::PhiPatInfo {
                index: 0,
                length: 2,
            })
        );
    }

    #[test]
    fn translated_phi_gapped_alignment_with_traceback_builds_pattern_script() {
        let matrix = phi_test_matrix();
        let pattern_blk =
            crate::pattern::sphi_pattern_search_blk_new("A-B", false, None).expect("pattern");
        let query = crate::encoding::encode_ncbistdaa_sequence(b"XXABYY");
        let subject = crate::encoding::encode_ncbistdaa_sequence(b"XXABYY");
        let mut gap_align = crate::blast_kappa::blast_gap_align_struct_new(12).expect("gap align");

        assert_eq!(
            phi_gapped_alignment_with_traceback(
                &query,
                &subject,
                Some(&mut gap_align),
                &matrix,
                5,
                2,
                12,
                2,
                2,
                query.len() as i32,
                subject.len() as i32,
                2,
                2,
                Some(&pattern_blk),
            ),
            0
        );

        assert_eq!(gap_align.score, 8);
        assert_eq!((gap_align.query_start, gap_align.subject_start), (0, 0));
        assert_eq!((gap_align.query_stop, gap_align.subject_stop), (6, 6));
        assert_eq!(
            gap_align
                .edit_script
                .as_ref()
                .map(|script| script.ops.as_slice()),
            Some(&[(GapAlignOpType::Sub, 6)][..])
        );
    }

    #[test]
    fn translated_phi_traceback_from_hsp_list_updates_and_reaps() {
        let matrix = phi_test_matrix();
        let pattern_blk =
            crate::pattern::sphi_pattern_search_blk_new("A-B", false, None).expect("pattern");
        let pattern_info = crate::pattern::SphiQueryInfo {
            num_patterns: 1,
            occurrences: vec![crate::pattern::SphiPatternInfo {
                offset: 2,
                length: 2,
            }],
            allocated_size: 1,
            probability: 0.0,
            pattern: None,
        };
        let query = crate::encoding::encode_ncbistdaa_sequence(b"XXABYY");
        let subject = crate::encoding::encode_ncbistdaa_sequence(b"XXABYY");
        let mut gap_align = crate::blast_kappa::blast_gap_align_struct_new(12).expect("gap align");
        let mut hsp_list = crate::hspstream::blast_hsp_list_new(10);
        crate::hspstream::blast_hsp_list_save_hsp(
            &mut hsp_list,
            crate::hspstream::Hsp {
                score: 4,
                num_ident: 0,
                bit_score: 0.0,
                evalue: f64::MAX,
                query_offset: 2,
                query_end: 4,
                query_gapped_start: 2,
                subject_offset: 2,
                subject_end: 4,
                subject_gapped_start: 2,
                context: 0,
                query_frame: 1,
                subject_frame: 0,
                num_gaps: 0,
                comp_adjustment_method: 0,
                edit_script: None,
                pat_info: Some(crate::hspstream::PhiPatInfo {
                    index: 0,
                    length: 2,
                }),
                map_info: None,
            },
        );
        let scoring_options = crate::options::ScoringOptions::new_blastp();
        let scoring_params = crate::parameters::ScoringParameters {
            options: scoring_options,
            reward: 0,
            penalty: 0,
            gap_open: 5,
            gap_extend: 2,
            scale_factor: 1.0,
        };
        let mut hit_options = crate::options::HitSavingOptions::default();
        hit_options.program_number = crate::program::PHI_BLASTP;
        hit_options.expect_value = f64::MAX;
        let hit_params = crate::parameters::HitSavingParameters {
            options: hit_options,
            cutoff_score_min: 1,
            low_score: Vec::new(),
            cutoffs: Vec::new(),
            link_hsp_params: None,
            prelim_evalue: 0.0,
        };

        assert_eq!(
            s_phi_traceback_from_hsp_list(
                crate::program::PHI_BLASTP,
                &mut hsp_list,
                &query,
                &subject,
                &mut gap_align,
                &matrix,
                &scoring_params,
                &hit_params,
                Some(&pattern_info),
                Some(&pattern_blk),
                0.267,
                0.041,
            ),
            0
        );

        assert_eq!(hsp_list.hsps.len(), 1);
        let hsp = &hsp_list.hsps[0];
        assert_eq!(hsp.score, 8);
        assert!(hsp.edit_script.is_some());
        assert!(hsp.evalue.is_finite());
        assert!(hsp.bit_score.is_finite());
    }

    #[test]
    fn out_of_frame_wrapper_switches_sequences_like_c() {
        let query = crate::encoding::encode_ncbistdaa_sequence(b"ACDE");
        let subject = crate::encoding::encode_ncbistdaa_sequence(b"ACDF");
        let scoring = OutOfFrameScoring::default();
        let mut q_start = 0;
        let mut s_start = 0;
        let mut block = crate::gapinfo::gap_prelim_edit_block_new();

        let score = s_out_of_frame_semi_gapped_align_wrap(
            &query,
            &subject,
            4,
            4,
            &mut q_start,
            &mut s_start,
            false,
            Some(&mut block),
            &scoring,
            0,
            false,
            true,
        );

        assert!(score > 0);
        assert!(q_start > 0);
        assert!(s_start > 0);
        assert!(!block.edit_ops.is_empty());
        assert!(block
            .edit_ops
            .iter()
            .any(|op| op.op_type == GapAlignOpType::Sub));
    }

    #[test]
    fn out_of_frame_wrapper_non_switch_matches_direct_tblastn_branch() {
        let query = crate::encoding::encode_ncbistdaa_sequence(b"ACDEFG");
        let subject = crate::encoding::encode_ncbistdaa_sequence(b"XXACDEFG");
        let scoring = OutOfFrameScoring::default();
        let mut direct_q_start = 0;
        let mut direct_s_start = 0;
        let direct = s_out_of_frame_gapped_align(
            &query,
            &subject,
            query.len() as i32,
            subject.len() as i32,
            &mut direct_q_start,
            &mut direct_s_start,
            true,
            None,
            &scoring,
            0,
            false,
        );
        let mut wrapped_q_start = 0;
        let mut wrapped_s_start = 0;
        let wrapped = s_out_of_frame_semi_gapped_align_wrap(
            &query,
            &subject,
            query.len() as i32,
            subject.len() as i32,
            &mut wrapped_q_start,
            &mut wrapped_s_start,
            true,
            None,
            &scoring,
            0,
            false,
            false,
        );

        assert_eq!(wrapped, direct);
        assert_eq!(wrapped_q_start, direct_q_start);
        assert_eq!(wrapped_s_start, direct_s_start);
    }

    #[test]
    fn out_of_frame_wrapper_non_switch_traceback_matches_direct_tblastn_branch() {
        let query = crate::encoding::encode_ncbistdaa_sequence(b"ACDEFG");
        let subject = crate::encoding::encode_ncbistdaa_sequence(b"XXACDEFG");
        let scoring = OutOfFrameScoring::default();
        let mut direct_q_start = 0;
        let mut direct_s_start = 0;
        let mut direct_block = crate::gapinfo::gap_prelim_edit_block_new();
        let direct = s_out_of_frame_gapped_align(
            &query,
            &subject,
            query.len() as i32,
            subject.len() as i32,
            &mut direct_q_start,
            &mut direct_s_start,
            false,
            Some(&mut direct_block),
            &scoring,
            0,
            false,
        );

        let mut wrapped_q_start = 0;
        let mut wrapped_s_start = 0;
        let mut wrapped_block = crate::gapinfo::gap_prelim_edit_block_new();
        let wrapped = s_out_of_frame_semi_gapped_align_wrap(
            &query,
            &subject,
            query.len() as i32,
            subject.len() as i32,
            &mut wrapped_q_start,
            &mut wrapped_s_start,
            false,
            Some(&mut wrapped_block),
            &scoring,
            0,
            false,
            false,
        );

        assert_eq!(wrapped, direct);
        assert_eq!(wrapped_q_start, direct_q_start);
        assert_eq!(wrapped_s_start, direct_s_start);
        assert_eq!(wrapped_block.edit_ops, direct_block.edit_ops);
        assert!(!wrapped_block.edit_ops.is_empty());
    }

    #[test]
    fn out_of_frame_wrapper_switch_score_only_matches_swapped_direct_branch() {
        let query = crate::encoding::encode_ncbistdaa_sequence(b"ACDEFG");
        let subject = crate::encoding::encode_ncbistdaa_sequence(b"XXACDEFG");
        let scoring = OutOfFrameScoring::default();
        let mut direct_subject_start = 0;
        let mut direct_query_start = 0;
        let direct = s_out_of_frame_gapped_align(
            &subject,
            &query,
            subject.len() as i32,
            query.len() as i32,
            &mut direct_subject_start,
            &mut direct_query_start,
            true,
            None,
            &scoring,
            0,
            false,
        );
        let mut wrapped_query_start = 0;
        let mut wrapped_subject_start = 0;
        let mut block = crate::gapinfo::gap_prelim_edit_block_new();
        let wrapped = s_out_of_frame_semi_gapped_align_wrap(
            &query,
            &subject,
            query.len() as i32,
            subject.len() as i32,
            &mut wrapped_query_start,
            &mut wrapped_subject_start,
            true,
            Some(&mut block),
            &scoring,
            0,
            false,
            true,
        );

        assert_eq!(wrapped, direct);
        assert_eq!(wrapped_query_start, direct_query_start);
        assert_eq!(wrapped_subject_start, direct_subject_start);
        assert!(block.edit_ops.is_empty());
    }

    #[test]
    fn out_of_frame_wrapper_switch_traceback_matches_swapped_direct_branch() {
        let query = crate::encoding::encode_ncbistdaa_sequence(b"ACDEFG");
        let subject = crate::encoding::encode_ncbistdaa_sequence(b"XXACDEFG");
        let scoring = OutOfFrameScoring::default();
        let mut direct_subject_start = 0;
        let mut direct_query_start = 0;
        let mut direct_block = crate::gapinfo::gap_prelim_edit_block_new();
        let direct = s_out_of_frame_gapped_align(
            &subject,
            &query,
            subject.len() as i32,
            query.len() as i32,
            &mut direct_subject_start,
            &mut direct_query_start,
            false,
            Some(&mut direct_block),
            &scoring,
            0,
            false,
        );

        let mut wrapped_query_start = 0;
        let mut wrapped_subject_start = 0;
        let mut wrapped_block = crate::gapinfo::gap_prelim_edit_block_new();
        let wrapped = s_out_of_frame_semi_gapped_align_wrap(
            &query,
            &subject,
            query.len() as i32,
            subject.len() as i32,
            &mut wrapped_query_start,
            &mut wrapped_subject_start,
            false,
            Some(&mut wrapped_block),
            &scoring,
            0,
            false,
            true,
        );

        assert_eq!(wrapped, direct);
        assert_eq!(wrapped_query_start, direct_query_start);
        assert_eq!(wrapped_subject_start, direct_subject_start);
        assert_eq!(wrapped_block.edit_ops, direct_block.edit_ops);
        assert!(!wrapped_block.edit_ops.is_empty());
    }

    #[test]
    fn out_of_frame_wrapper_switch_empty_input_maps_c_sentinel_offsets() {
        let scoring = OutOfFrameScoring::default();
        let mut q_start = 99;
        let mut s_start = 99;
        let mut block = crate::gapinfo::gap_prelim_edit_block_new();

        let score = s_out_of_frame_semi_gapped_align_wrap(
            &[],
            &[],
            0,
            0,
            &mut q_start,
            &mut s_start,
            false,
            Some(&mut block),
            &scoring,
            0,
            false,
            true,
        );

        assert_eq!(score, 0);
        assert_eq!(q_start, -2);
        assert_eq!(s_start, 0);
        assert!(block.edit_ops.is_empty());
    }

    #[test]
    fn out_of_frame_wrapper_non_switch_empty_input_keeps_c_sentinel_offsets() {
        let scoring = OutOfFrameScoring::default();
        let mut q_start = 99;
        let mut s_start = 99;
        let mut block = crate::gapinfo::gap_prelim_edit_block_new();

        let score = s_out_of_frame_semi_gapped_align_wrap(
            &[],
            &[],
            0,
            0,
            &mut q_start,
            &mut s_start,
            false,
            Some(&mut block),
            &scoring,
            0,
            false,
            false,
        );

        assert_eq!(score, 0);
        assert_eq!(q_start, 0);
        assert_eq!(s_start, -2);
        assert!(block.edit_ops.is_empty());
    }

    #[test]
    fn out_of_frame_wrapper_switch_score_only_empty_input_maps_c_sentinel_offsets() {
        let scoring = OutOfFrameScoring::default();
        let mut q_start = 99;
        let mut s_start = 99;
        let mut block = crate::gapinfo::gap_prelim_edit_block_new();

        let score = s_out_of_frame_semi_gapped_align_wrap(
            &[],
            &[],
            0,
            0,
            &mut q_start,
            &mut s_start,
            true,
            Some(&mut block),
            &scoring,
            0,
            false,
            true,
        );

        assert_eq!(score, 0);
        assert_eq!(q_start, -2);
        assert_eq!(s_start, 0);
        assert!(block.edit_ops.is_empty());
    }

    #[test]
    fn out_of_frame_reversed_matrix_branch_reads_query_backwards_like_c() {
        let a = crate::encoding::AMINOACID_TO_NCBISTDAA[b'A' as usize];
        let c = crate::encoding::AMINOACID_TO_NCBISTDAA[b'C' as usize];
        let d = crate::encoding::AMINOACID_TO_NCBISTDAA[b'D' as usize];
        let query = vec![a, c, d];
        let scoring = OutOfFrameScoring::default();

        assert_eq!(scoring.oof_row_score(0, &query, query.len(), 1, d, true), 6);
        assert_eq!(
            scoring.oof_row_score(0, &query, query.len(), 1, d, false),
            -3
        );
    }

    #[test]
    fn out_of_frame_subject_letter_selection_matches_c_pointer_offset() {
        let subject = [10u8, 11, 12, 13, 14];

        assert_eq!(oof_subject_letter(&subject, 0, false), 0);
        assert_eq!(oof_subject_letter(&subject, 1, false), 0);
        assert_eq!(oof_subject_letter(&subject, 2, false), 10);
        assert_eq!(oof_subject_letter(&subject, 4, false), 12);

        assert_eq!(oof_subject_letter(&subject, 0, true), 14);
        assert_eq!(oof_subject_letter(&subject, 2, true), 12);
    }

    #[test]
    fn out_of_frame_initial_trace_scripts_follow_c_frame_ties() {
        assert_eq!(
            oof_initial_trace_script(0, 10, 5, 0, 10, 0, 1),
            SCRIPT_AHEAD_TWO_FRAMES
        );
        assert_eq!(
            oof_initial_trace_script(0, 10, 5, 0, 3, 0, 1),
            SCRIPT_NEXT_PLUS_TWO_FRAMES
        );
        assert_eq!(
            oof_initial_trace_script(0, 5, 10, 0, 0, 10, 1),
            SCRIPT_AHEAD_ONE_FRAME
        );
        assert_eq!(
            oof_initial_trace_script(0, 5, 10, 0, 0, 3, 1),
            SCRIPT_NEXT_PLUS_ONE_FRAME
        );

        assert_eq!(
            oof_initial_trace_script(1, 10, 5, 10, 0, 0, 1),
            SCRIPT_AHEAD_ONE_FRAME
        );
        assert_eq!(
            oof_initial_trace_script(1, 10, 5, 3, 0, 0, 1),
            SCRIPT_NEXT_PLUS_ONE_FRAME
        );
        assert_eq!(
            oof_initial_trace_script(1, 5, 10, 0, 0, 10, 1),
            SCRIPT_AHEAD_TWO_FRAMES
        );
        assert_eq!(
            oof_initial_trace_script(1, 5, 10, 0, 0, 3, 1),
            SCRIPT_NEXT_PLUS_TWO_FRAMES
        );

        assert_eq!(
            oof_initial_trace_script(2, 10, 5, 10, 0, 0, 1),
            SCRIPT_AHEAD_TWO_FRAMES
        );
        assert_eq!(
            oof_initial_trace_script(2, 10, 5, 3, 0, 0, 1),
            SCRIPT_NEXT_PLUS_TWO_FRAMES
        );
        assert_eq!(
            oof_initial_trace_script(2, 5, 10, 0, 10, 0, 1),
            SCRIPT_AHEAD_ONE_FRAME
        );
        assert_eq!(
            oof_initial_trace_script(2, 5, 10, 0, 3, 0, 1),
            SCRIPT_NEXT_PLUS_ONE_FRAME
        );

        assert_eq!(
            oof_initial_trace_script(0, 10, 5, 9, 3, 0, 1),
            SCRIPT_NEXT_IN_FRAME
        );
    }

    #[test]
    fn out_of_frame_script_op_codes_map_to_gap_ops_like_c() {
        assert_eq!(oof_script_to_op(SCRIPT_GAP_IN_A), GapAlignOpType::Del);
        assert_eq!(
            oof_script_to_op(SCRIPT_AHEAD_ONE_FRAME),
            GapAlignOpType::Del2
        );
        assert_eq!(
            oof_script_to_op(SCRIPT_AHEAD_TWO_FRAMES),
            GapAlignOpType::Del1
        );
        assert_eq!(oof_script_to_op(SCRIPT_NEXT_IN_FRAME), GapAlignOpType::Sub);
        assert_eq!(
            oof_script_to_op(SCRIPT_NEXT_PLUS_ONE_FRAME),
            GapAlignOpType::Ins1
        );
        assert_eq!(
            oof_script_to_op(SCRIPT_NEXT_PLUS_TWO_FRAMES),
            GapAlignOpType::Ins2
        );
        assert_eq!(oof_script_to_op(SCRIPT_GAP_IN_B), GapAlignOpType::Ins);
    }

    #[test]
    fn out_of_frame_gapped_align_score_only_uses_specialized_oof_dp() {
        let query = crate::encoding::encode_ncbistdaa_sequence(b"ACDEFG");
        let subject = crate::encoding::encode_ncbistdaa_sequence(b"XXACDEFG");
        let scoring = OutOfFrameScoring::default();
        let mut a_offset = 0;
        let mut b_offset = 0;
        let mut block = crate::gapinfo::gap_prelim_edit_block_new();

        let score = s_out_of_frame_gapped_align(
            &query,
            &subject,
            query.len() as i32,
            subject.len() as i32,
            &mut a_offset,
            &mut b_offset,
            true,
            Some(&mut block),
            &scoring,
            0,
            false,
        );

        assert!(score > 0);
        assert!(a_offset > 0);
        assert!(b_offset >= -2);
        assert!(
            block.edit_ops.is_empty(),
            "score-only path must not populate traceback"
        );
    }

    #[test]
    fn out_of_frame_gapped_align_score_only_empty_input_keeps_c_sentinel_offsets() {
        let scoring = OutOfFrameScoring::default();
        let mut a_offset = 99;
        let mut b_offset = 99;
        let mut block = crate::gapinfo::gap_prelim_edit_block_new();

        let score = s_out_of_frame_gapped_align(
            &[],
            &[],
            0,
            0,
            &mut a_offset,
            &mut b_offset,
            true,
            Some(&mut block),
            &scoring,
            0,
            false,
        );

        assert_eq!(score, 0);
        assert_eq!(a_offset, 0);
        assert_eq!(b_offset, -2);
        assert!(block.edit_ops.is_empty());
    }

    #[test]
    fn out_of_frame_gapped_align_traceback_empty_input_keeps_c_sentinel_offsets() {
        let scoring = OutOfFrameScoring::default();
        let mut a_offset = 99;
        let mut b_offset = 99;
        let mut block = crate::gapinfo::gap_prelim_edit_block_new();

        let score = s_out_of_frame_gapped_align(
            &[],
            &[],
            0,
            0,
            &mut a_offset,
            &mut b_offset,
            false,
            Some(&mut block),
            &scoring,
            0,
            false,
        );

        assert_eq!(score, 0);
        assert_eq!(a_offset, 0);
        assert_eq!(b_offset, -2);
        assert!(block.edit_ops.is_empty());
    }

    #[test]
    fn out_of_frame_traceback_uses_same_oof_dp_score() {
        let query = crate::encoding::encode_ncbistdaa_sequence(b"ACDEFG");
        let subject = crate::encoding::encode_ncbistdaa_sequence(b"XXACDEFG");
        let scoring = OutOfFrameScoring::default();
        let mut score_a_offset = 0;
        let mut score_b_offset = 0;
        let score_only = s_out_of_frame_gapped_align(
            &query,
            &subject,
            query.len() as i32,
            subject.len() as i32,
            &mut score_a_offset,
            &mut score_b_offset,
            true,
            None,
            &scoring,
            0,
            false,
        );

        let mut trace_a_offset = 0;
        let mut trace_b_offset = 0;
        let mut block = crate::gapinfo::gap_prelim_edit_block_new();
        let trace_score = s_out_of_frame_gapped_align(
            &query,
            &subject,
            query.len() as i32,
            subject.len() as i32,
            &mut trace_a_offset,
            &mut trace_b_offset,
            false,
            Some(&mut block),
            &scoring,
            0,
            false,
        );

        assert_eq!(trace_score, score_only);
        assert_eq!(trace_a_offset, score_a_offset);
        assert_eq!(trace_b_offset, score_b_offset);
        assert!(!block.edit_ops.is_empty());
    }

    #[test]
    fn out_of_frame_reversed_traceback_uses_same_oof_dp_score() {
        let query = crate::encoding::encode_ncbistdaa_sequence(b"ACDEFG");
        let subject = crate::encoding::encode_ncbistdaa_sequence(b"ACDEFGXX");
        let scoring = OutOfFrameScoring::default();
        let mut score_a_offset = 0;
        let mut score_b_offset = 0;
        let score_only = s_out_of_frame_gapped_align(
            &query,
            &subject,
            query.len() as i32,
            subject.len() as i32,
            &mut score_a_offset,
            &mut score_b_offset,
            true,
            None,
            &scoring,
            0,
            true,
        );

        let mut trace_a_offset = 0;
        let mut trace_b_offset = 0;
        let mut block = crate::gapinfo::gap_prelim_edit_block_new();
        let trace_score = s_out_of_frame_gapped_align(
            &query,
            &subject,
            query.len() as i32,
            subject.len() as i32,
            &mut trace_a_offset,
            &mut trace_b_offset,
            false,
            Some(&mut block),
            &scoring,
            0,
            true,
        );

        assert_eq!(trace_score, score_only);
        assert_eq!(trace_a_offset, score_a_offset);
        assert_eq!(trace_b_offset, score_b_offset);
        assert!(!block.edit_ops.is_empty());
    }

    #[test]
    fn out_of_frame_traceback_without_edit_block_keeps_score_and_offsets() {
        let query = crate::encoding::encode_ncbistdaa_sequence(b"ACDEFG");
        let subject = crate::encoding::encode_ncbistdaa_sequence(b"XXACDEFG");
        let scoring = OutOfFrameScoring::default();
        let mut with_block_q_offset = 0;
        let mut with_block_s_offset = 0;
        let mut block = crate::gapinfo::gap_prelim_edit_block_new();
        let with_block = s_out_of_frame_align_with_traceback(
            &query,
            &subject,
            query.len() as i32,
            subject.len() as i32,
            &mut with_block_q_offset,
            &mut with_block_s_offset,
            Some(&mut block),
            &scoring,
            0,
            false,
        );

        let mut no_block_q_offset = 0;
        let mut no_block_s_offset = 0;
        let no_block = s_out_of_frame_align_with_traceback(
            &query,
            &subject,
            query.len() as i32,
            subject.len() as i32,
            &mut no_block_q_offset,
            &mut no_block_s_offset,
            None,
            &scoring,
            0,
            false,
        );

        assert_eq!(no_block, with_block);
        assert_eq!(no_block_q_offset, with_block_q_offset);
        assert_eq!(no_block_s_offset, with_block_s_offset);
        assert!(!block.edit_ops.is_empty());
    }

    #[test]
    fn out_of_frame_traceback_step_records_gap_flags_like_c() {
        let mut score_array = vec![OofScoreCell {
            best: 0,
            best_gap: 10,
        }];
        let mut score_row = 5;
        let mut score_col = 0;
        let mut other1 = i32::MIN / 4;
        let mut other2 = i32::MIN / 4;
        let mut best_score = 20;
        let mut a_offset = 0;
        let mut b_offset = 0;
        let mut first_b_index = 0;
        let mut last_b_index = 0;

        let gap_open = oof_step_frame_with_script(
            &mut score_array,
            0,
            &mut score_row,
            &mut score_col,
            &mut other1,
            &mut other2,
            &mut best_score,
            &mut a_offset,
            &mut b_offset,
            1,
            0,
            2,
            1,
            1,
            100,
            &mut first_b_index,
            &mut last_b_index,
            SCRIPT_SUB,
            false,
        );
        assert_eq!(gap_open, SCRIPT_OOF_OPEN_GAP | SCRIPT_GAP_IN_A);

        let mut score_array = vec![OofScoreCell {
            best: 0,
            best_gap: 10,
        }];
        let mut score_row = 0;
        let mut score_col = 10;
        let mut other1 = i32::MIN / 4;
        let mut other2 = i32::MIN / 4;
        let mut best_score = 0;
        let mut first_b_index = 0;
        let mut last_b_index = 0;
        let extended = oof_step_frame_with_script(
            &mut score_array,
            0,
            &mut score_row,
            &mut score_col,
            &mut other1,
            &mut other2,
            &mut best_score,
            &mut a_offset,
            &mut b_offset,
            1,
            0,
            2,
            1,
            1,
            100,
            &mut first_b_index,
            &mut last_b_index,
            SCRIPT_SUB,
            false,
        );
        assert_eq!(extended & SCRIPT_EXTEND_GAP_A, SCRIPT_EXTEND_GAP_A);

        let mut score_array = vec![OofScoreCell {
            best: 0,
            best_gap: 0,
        }];
        let mut score_row = 9;
        let mut score_col = 10;
        let mut other1 = i32::MIN / 4;
        let mut other2 = i32::MIN / 4;
        let mut best_score = 0;
        let mut first_b_index = 0;
        let mut last_b_index = 0;
        let extended = oof_step_frame_with_script(
            &mut score_array,
            0,
            &mut score_row,
            &mut score_col,
            &mut other1,
            &mut other2,
            &mut best_score,
            &mut a_offset,
            &mut b_offset,
            1,
            0,
            2,
            1,
            1,
            100,
            &mut first_b_index,
            &mut last_b_index,
            SCRIPT_SUB,
            false,
        );
        assert_eq!(extended & SCRIPT_EXTEND_GAP_B, SCRIPT_EXTEND_GAP_B);
    }

    #[test]
    fn out_of_frame_traceback_step_column_tie_preference_matches_c() {
        let run_tie = |gap_tie_prefers_column| {
            let mut score_array = vec![OofScoreCell {
                best: 0,
                best_gap: 10,
            }];
            let mut score_row = 10;
            let mut score_col = 0;
            let mut other1 = i32::MIN / 4;
            let mut other2 = i32::MIN / 4;
            let mut best_score = 30;
            let mut a_offset = 0;
            let mut b_offset = 0;
            let mut first_b_index = 0;
            let mut last_b_index = 0;

            oof_step_frame_with_script(
                &mut score_array,
                0,
                &mut score_row,
                &mut score_col,
                &mut other1,
                &mut other2,
                &mut best_score,
                &mut a_offset,
                &mut b_offset,
                1,
                0,
                2,
                1,
                1,
                100,
                &mut first_b_index,
                &mut last_b_index,
                SCRIPT_SUB,
                gap_tie_prefers_column,
            )
        };

        assert_eq!(run_tie(false), SCRIPT_OOF_OPEN_GAP | SCRIPT_GAP_IN_B);
        assert_eq!(run_tie(true), SCRIPT_OOF_OPEN_GAP | SCRIPT_GAP_IN_A);
    }

    #[test]
    fn oof_traceback_to_gap_edit_script_splits_frameshifts_and_truncates() {
        let mut rev = crate::gapinfo::gap_prelim_edit_block_new();
        crate::gapinfo::gap_prelim_edit_block_add(&mut rev, GapAlignOpType::Ins1, 2);
        crate::gapinfo::gap_prelim_edit_block_add(&mut rev, GapAlignOpType::Sub, 2);
        let mut fwd = crate::gapinfo::gap_prelim_edit_block_new();
        crate::gapinfo::gap_prelim_edit_block_add(&mut fwd, GapAlignOpType::Sub, 3);
        crate::gapinfo::gap_prelim_edit_block_add(&mut fwd, GapAlignOpType::Del, 1);
        let mut script = None;

        assert_eq!(
            crate::gapinfo::s_blast_oof_traceback_to_gap_edit_script(
                Some(&rev),
                Some(&fwd),
                18,
                Some(&mut script),
            ),
            0
        );
        let script = script.expect("script");
        assert!(script
            .ops
            .iter()
            .any(|(op, count)| *op == GapAlignOpType::Ins1 && *count == 1));
        assert!(script
            .ops
            .iter()
            .all(|(op, count)| (*op as u32) % 3 == 0 || *count == 1));
        assert!(script.ops.iter().map(|(_, count)| *count).sum::<i32>() > 0);
    }
}
