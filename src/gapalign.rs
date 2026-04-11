//! Rust equivalent of blast_gapalign.c — gapped alignment.
//! Implements Smith-Waterman style gapped alignment with affine gap penalties.

use crate::gapinfo::{GapAlignOpType, GapEditScript};

/// Perform gapped alignment between two sequences using dynamic programming.
///
/// Returns (score, edit_script) or None if no significant alignment found.
pub fn gapped_align(
    query: &[u8],
    subject: &[u8],
    q_start: i32,
    s_start: i32,
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    x_dropoff: i32,
) -> Option<(i32, GapEditScript)> {
    let q = &query[q_start as usize..];
    let s = &subject[s_start as usize..];
    let q_len = q.len();
    let s_len = s.len();

    if q_len == 0 || s_len == 0 {
        return None;
    }

    // DP matrices: H[i][j] = best score ending at (i,j)
    // E[i][j] = best score with gap in query (deletion)
    // F[i][j] = best score with gap in subject (insertion)
    let cols = s_len + 1;
    let mut h_prev = vec![0i32; cols];
    let mut e_prev = vec![i32::MIN / 2; cols];
    let mut h_curr = vec![0i32; cols];
    let mut e_curr = vec![i32::MIN / 2; cols];

    let mut best_score = 0i32;
    let mut best_i = 0usize;
    let mut best_j = 0usize;

    let gap_open_extend = gap_open + gap_extend;

    for i in 1..=q_len {
        let mut f = i32::MIN / 2;
        h_curr[0] = 0;
        e_curr[0] = i32::MIN / 2;

        for j in 1..=cols - 1 {
            // Match/mismatch
            let match_score = if q[i - 1] == s[j - 1] { reward } else { penalty };
            let diag = h_prev[j - 1] + match_score;

            // Gap in query (deletion: advance subject)
            e_curr[j] = (h_prev[j] - gap_open_extend).max(e_prev[j] - gap_extend);

            // Gap in subject (insertion: advance query)
            f = (h_curr[j - 1] - gap_open_extend).max(f - gap_extend);

            // Best of all options (must be >= 0 for local alignment)
            h_curr[j] = diag.max(e_curr[j]).max(f).max(0);

            if h_curr[j] > best_score {
                best_score = h_curr[j];
                best_i = i;
                best_j = j;
            }
        }

        // X-dropoff check: if all values in this row are below best - x_dropoff, stop
        let row_max = *h_curr.iter().max().unwrap_or(&0);
        if best_score - row_max > x_dropoff {
            break;
        }

        std::mem::swap(&mut h_prev, &mut h_curr);
        std::mem::swap(&mut e_prev, &mut e_curr);
    }

    if best_score <= 0 {
        return None;
    }

    // Traceback from best position through the DP matrix
    // We need to recompute or store the full matrices for traceback.
    // Recompute H, E, F with full storage for traceback.
    let mut h_full = vec![vec![0i32; cols]; q_len + 1];
    let mut e_full = vec![vec![i32::MIN / 2; cols]; q_len + 1];
    let mut f_full = vec![vec![i32::MIN / 2; cols]; q_len + 1];
    // 0=diag, 1=left(gap in subject=ins), 2=up(gap in query=del), 255=reset
    let mut trace = vec![vec![255u8; cols]; q_len + 1];

    for i in 1..=q_len {
        let mut fi = i32::MIN / 2;
        for j in 1..=cols - 1 {
            let match_score = if q[i - 1] == s[j - 1] { reward } else { penalty };
            let diag = h_full[i - 1][j - 1] + match_score;

            e_full[i][j] = (h_full[i - 1][j] - gap_open_extend).max(e_full[i - 1][j] - gap_extend);
            fi = (h_full[i][j - 1] - gap_open_extend).max(fi - gap_extend);
            f_full[i][j] = fi;

            let val = diag.max(e_full[i][j]).max(f_full[i][j]).max(0);
            h_full[i][j] = val;

            if val == 0 {
                trace[i][j] = 255;
            } else if val == diag {
                trace[i][j] = 0;
            } else if val == f_full[i][j] {
                trace[i][j] = 1; // gap in subject
            } else {
                trace[i][j] = 2; // gap in query
            }
        }
    }

    let mut ops: Vec<(GapAlignOpType, i32)> = Vec::new();
    let mut i = best_i;
    let mut j = best_j;

    while i > 0 && j > 0 && h_full[i][j] > 0 {
        match trace[i][j] {
            0 => {
                if let Some(last) = ops.last_mut() {
                    if last.0 == GapAlignOpType::Sub { last.1 += 1; } else { ops.push((GapAlignOpType::Sub, 1)); }
                } else { ops.push((GapAlignOpType::Sub, 1)); }
                i -= 1;
                j -= 1;
            }
            1 => {
                if let Some(last) = ops.last_mut() {
                    if last.0 == GapAlignOpType::Ins { last.1 += 1; } else { ops.push((GapAlignOpType::Ins, 1)); }
                } else { ops.push((GapAlignOpType::Ins, 1)); }
                j -= 1;
            }
            2 => {
                if let Some(last) = ops.last_mut() {
                    if last.0 == GapAlignOpType::Del { last.1 += 1; } else { ops.push((GapAlignOpType::Del, 1)); }
                } else { ops.push((GapAlignOpType::Del, 1)); }
                i -= 1;
            }
            _ => break,
        }
    }

    ops.reverse();
    let mut esp = GapEditScript::new();
    for (op, count) in ops {
        esp.push(op, count);
    }

    Some((best_score, esp))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_perfect_match() {
        let query = vec![0u8, 1, 2, 3, 0, 1, 2, 3]; // ACGTACGT
        let subject = vec![0u8, 1, 2, 3, 0, 1, 2, 3]; // ACGTACGT

        let result = gapped_align(&query, &subject, 0, 0, 2, -3, 5, 2, 100);
        assert!(result.is_some());
        let (score, _esp) = result.unwrap();
        assert_eq!(score, 16); // 8 matches * 2
    }

    #[test]
    fn test_with_mismatch() {
        let query = vec![0u8, 1, 2, 3]; // ACGT
        let subject = vec![0u8, 0, 2, 3]; // AAGT

        let result = gapped_align(&query, &subject, 0, 0, 2, -3, 5, 2, 100);
        assert!(result.is_some());
        let (score, _) = result.unwrap();
        // Local alignment: G-G(2) + T-T(2) = 4 is better than full 2-3+2+2=3
        assert!(score >= 3, "score={} should be >= 3", score);
    }

    #[test]
    fn test_no_alignment() {
        let query = vec![0u8, 0, 0, 0]; // AAAA
        let subject = vec![3u8, 3, 3, 3]; // TTTT

        let result = gapped_align(&query, &subject, 0, 0, 1, -3, 5, 2, 5);
        // Score would be -12, so no alignment
        assert!(result.is_none());
    }
}
