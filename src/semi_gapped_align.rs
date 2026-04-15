//! Verbatim port of NCBI Blast_SemiGappedAlign (score-only, no traceback).
//! Source: ncbi-cxx-toolkit-public/src/algo/blast/core/blast_gapalign.c lines 736-960
//!
//! This replaces protein_gapped_score_one_dir for accurate gapped extension
//! matching NCBI BLAST+ behavior.

use crate::matrix::AA_SIZE;

const MININT: i32 = i32::MIN / 2;

pub struct GapDP {
    pub best: i32,
    pub best_gap: i32,
}

/// Score-only gapped extension in one direction.
/// Verbatim port of NCBI Blast_SemiGappedAlign (score_only=TRUE path).
///
/// Returns (best_score, a_offset, b_offset) where offsets are how far
/// the extension reached in A and B respectively.
pub fn semi_gapped_align(
    a: &[u8], // query (or reversed query)
    b: &[u8], // subject (or reversed subject)
    m: usize, // length of A to consider
    n: usize, // length of B to consider
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
    gap_open: i32,
    gap_extend: i32,
    x_dropoff: i32,
    reverse: bool,
) -> (i32, usize, usize) {
    if m == 0 || n == 0 {
        return (0, 0, 0);
    }

    let gap_open_extend = gap_open + gap_extend;
    let mut x_dropoff = x_dropoff;
    if x_dropoff < gap_open_extend {
        x_dropoff = gap_open_extend;
    }

    // Initial band size (NCBI: num_extra_cells)
    let num_extra_cells = if gap_extend > 0 {
        (x_dropoff / gap_extend) as usize + 3
    } else {
        n + 3
    };

    let mut score_array: Vec<GapDP> = Vec::with_capacity(num_extra_cells + 100);
    score_array.push(GapDP {
        best: 0,
        best_gap: -gap_open_extend,
    });

    let mut score = -gap_open_extend;
    let mut i = 1;
    while i <= n && score >= -x_dropoff {
        score_array.push(GapDP {
            best: score,
            best_gap: score - gap_open_extend,
        });
        score -= gap_extend;
        i += 1;
    }

    let mut b_size = score_array.len();
    let mut best_score = 0i32;
    let mut a_off = 0usize;
    let mut b_off = 0usize;
    let mut first_b_index = 0usize;

    for a_index in 1..=m {
        // Pick the row of the score matrix for A[a_index]
        let a_letter = if reverse {
            a[m - a_index] as usize
        } else {
            a[a_index] as usize
        };
        if a_letter >= AA_SIZE {
            continue;
        }
        let matrix_row = &matrix[a_letter];

        // Initialize running-score variables
        let mut sc = MININT; // score
        let mut score_gap_row = MININT; // best score with gap in A (horizontal)
        let mut last_b_index = first_b_index;

        for b_index in first_b_index..b_size {
            let b_letter = if reverse {
                match n.checked_sub(1 + b_index) {
                    Some(idx) if idx < b.len() => b[idx] as usize,
                    _ => break,
                }
            } else {
                let idx = b_index + 1;
                if idx >= b.len() {
                    break;
                }
                b[idx] as usize
            };
            if b_letter >= AA_SIZE {
                continue;
            }

            let score_gap_col = score_array[b_index].best_gap;
            let next_score = score_array[b_index].best + matrix_row[b_letter];

            // sc = max(sc, score_gap_col, score_gap_row)
            if sc < score_gap_col {
                sc = score_gap_col;
            }
            if sc < score_gap_row {
                sc = score_gap_row;
            }

            if best_score - sc > x_dropoff {
                // Failed x-dropoff
                if b_index == first_b_index {
                    first_b_index += 1;
                }
                score_array[b_index].best = MININT;
                score_array[b_index].best_gap = MININT;
            } else {
                last_b_index = b_index;
                if sc > best_score {
                    best_score = sc;
                    a_off = a_index;
                    b_off = b_index;
                }

                // Update gap scores
                score_gap_row -= gap_extend;
                let score_gap_col_ext = score_gap_col - gap_extend;
                score_array[b_index].best_gap = (sc - gap_open_extend).max(score_gap_col_ext);
                score_gap_row = (sc - gap_open_extend).max(score_gap_row);
                score_array[b_index].best = sc;
            }

            sc = next_score;
        }

        // Check if all positions failed x-dropoff
        if first_b_index >= b_size {
            break;
        }

        if last_b_index < b_size - 1 {
            // This row ended early — shrink band
            b_size = last_b_index + 1;
        } else {
            // Expand band: NCBI grows b_size while score_gap_row is viable
            while score_gap_row >= (best_score - x_dropoff) && b_size <= n {
                if b_size >= score_array.len() {
                    score_array.push(GapDP {
                        best: MININT,
                        best_gap: MININT,
                    });
                }
                score_array[b_size].best = score_gap_row;
                score_array[b_size].best_gap = score_gap_row - gap_open_extend;
                score_gap_row -= gap_extend;
                b_size += 1;
            }
        }
    }

    (best_score, a_off, b_off + 1) // +1 to match NCBI's 1-based offset
}
