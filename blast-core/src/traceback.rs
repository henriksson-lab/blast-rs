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
}
