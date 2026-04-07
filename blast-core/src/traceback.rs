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
) -> (i32, GapEditScript) {
    let q = &query[q_start..q_end];
    let s = &subject[s_start..s_end];
    let m = q.len();
    let n = s.len();

    if m == 0 || n == 0 {
        return (0, GapEditScript::new());
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

    (best_score, esp)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_traceback_perfect() {
        let q = vec![0u8, 1, 2, 3, 0, 1, 2, 3];
        let s = vec![0u8, 1, 2, 3, 0, 1, 2, 3];
        let (score, esp) = traceback_align(&q, &s, 0, 8, 0, 8, 2, -3, 5, 2);
        assert_eq!(score, 16);
        assert_eq!(esp.ops.len(), 1);
        assert_eq!(esp.ops[0], (GapAlignOpType::Sub, 8));
    }

    #[test]
    fn test_traceback_with_mismatch() {
        let q = vec![0u8, 1, 2, 3, 0, 1, 2, 3];
        let s = vec![0u8, 1, 2, 3, 3, 1, 2, 3]; // pos 4: 0→3
        let (score, esp) = traceback_align(&q, &s, 0, 8, 0, 8, 2, -3, 5, 2);
        assert!(score > 0);
        // Should have 8 substitutions (one mismatch)
        let total: i32 = esp.ops.iter().map(|(_, n)| *n).sum();
        assert!(total >= 7, "alignment should be at least 7 bases");
    }
}
