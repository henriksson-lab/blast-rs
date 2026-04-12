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

    /// Gapped DP alignment with mismatches but no gaps needed.
    #[test]
    fn test_gapped_align_with_mismatch() {
        // Query:   ACGTACGT (0,1,2,3,0,1,2,3)
        // Subject: ACGAACGA (0,1,2,0,0,1,2,0) — mismatches at pos 3 and 7
        let query   = vec![0u8, 1, 2, 3, 0, 1, 2, 3];
        let subject = vec![0u8, 1, 2, 0, 0, 1, 2, 0];
        let result = gapped_align(&query, &subject, 0, 0, 2, -3, 5, 2, 100);
        assert!(result.is_some());
        let (score, esp) = result.unwrap();
        // Local alignment finds best scoring sub-alignment
        assert!(score > 0, "score={}", score);
        // All operations should be Sub (no gaps)
        for (op, _) in &esp.ops {
            assert_eq!(*op, GapAlignOpType::Sub,
                "expected only Sub ops (no gaps), got {:?}", esp.ops);
        }
        let total: i32 = esp.ops.iter().map(|(_, n)| *n).sum();
        // Local alignment may trim terminal mismatches, so aligned length >= 6
        assert!(total >= 6, "should align at least 6 positions, got {}", total);
    }

    /// Subject has an extra base relative to query (insertion in subject).
    #[test]
    fn test_gapped_align_with_insertion() {
        // Query:   ACGT-ACGT (8 bp)
        // Subject: ACGTAACGT (9 bp) — extra A at position 4
        let query   = vec![0u8, 1, 2, 3, 0, 1, 2, 3];
        let subject = vec![0u8, 1, 2, 3, 0, 0, 1, 2, 3];
        let result = gapped_align(&query, &subject, 0, 0, 2, -3, 3, 1, 100);
        assert!(result.is_some());
        let (score, esp) = result.unwrap();
        // Best alignment: 8 matches (16) - gap_open+gap_extend (3+1=4) = 12
        assert!(score >= 12, "score={} should be >= 12", score);
        // Should contain a Del operation (gap in query = deletion)
        let has_gap = esp.ops.iter().any(|(op, _)|
            *op == GapAlignOpType::Del || *op == GapAlignOpType::Ins);
        assert!(has_gap, "expected a gap operation, ops={:?}", esp.ops);
    }

    /// Query has an extra base relative to subject (deletion in subject).
    #[test]
    fn test_gapped_align_with_deletion() {
        // Query:   ACGTAACGT (9 bp) — extra A at position 4
        // Subject: ACGT-ACGT (8 bp)
        let query   = vec![0u8, 1, 2, 3, 0, 0, 1, 2, 3];
        let subject = vec![0u8, 1, 2, 3, 0, 1, 2, 3];
        let result = gapped_align(&query, &subject, 0, 0, 2, -3, 3, 1, 100);
        assert!(result.is_some());
        let (score, esp) = result.unwrap();
        assert!(score >= 12, "score={} should be >= 12", score);
        let has_gap = esp.ops.iter().any(|(op, _)|
            *op == GapAlignOpType::Del || *op == GapAlignOpType::Ins);
        assert!(has_gap, "expected a gap operation, ops={:?}", esp.ops);
    }

    /// Verify X-dropoff terminates alignment correctly.
    /// With a very tight x_dropoff, a long mismatch region should stop the alignment.
    #[test]
    fn test_gapped_align_xdrop_cutoff() {
        // 4 matches then 10 mismatches then 4 matches
        // With tight x_dropoff the second matching region should not be reached
        let query   = vec![0u8, 1, 2, 3,  3, 3, 3, 3, 3, 3, 3, 3, 3, 3,  0, 1, 2, 3];
        let subject = vec![0u8, 1, 2, 3,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 1, 2, 3];
        // x_dropoff = 5: after 4 matches (score 8), mismatches bring score to 8-3=5, 5-3=2, 2-3=-1
        // When best-current > 5, we stop
        let result_tight = gapped_align(&query, &subject, 0, 0, 2, -3, 5, 2, 5);
        // With generous x_dropoff, might recover more
        let result_generous = gapped_align(&query, &subject, 0, 0, 2, -3, 5, 2, 100);
        // Both should find something
        assert!(result_tight.is_some());
        assert!(result_generous.is_some());
        let (score_tight, _) = result_tight.unwrap();
        let (score_generous, _) = result_generous.unwrap();
        // Tight x_dropoff score should be <= generous score
        assert!(score_tight <= score_generous,
            "tight score {} should be <= generous score {}",
            score_tight, score_generous);
    }

    /// Alignment requiring a gap of length > 1.
    #[test]
    fn test_gapped_align_long_gap() {
        // Query:   ACGTACGTACGT------ACGTACGTACGT  (24 bp)
        // Subject: ACGTACGTACGTAAAAAACGTACGTACGT   (29 bp) — 5 extra A's
        // Long flanks make the gap worthwhile even with affine penalties.
        let query: Vec<u8>   = [0,1,2,3,0,1,2,3,0,1,2,3, 0,1,2,3,0,1,2,3,0,1,2,3].to_vec();
        let subject: Vec<u8> = [0,1,2,3,0,1,2,3,0,1,2,3, 0,0,0,0,0, 0,1,2,3,0,1,2,3,0,1,2,3].to_vec();
        let result = gapped_align(&query, &subject, 0, 0, 2, -3, 3, 1, 100);
        assert!(result.is_some());
        let (score, esp) = result.unwrap();
        assert!(score > 0, "score={}", score);
        // Find the gap operation and verify its length > 1
        let gap_ops: Vec<_> = esp.ops.iter()
            .filter(|(op, _)| *op == GapAlignOpType::Del || *op == GapAlignOpType::Ins)
            .collect();
        assert!(!gap_ops.is_empty(), "expected gap operations, ops={:?}", esp.ops);
        let max_gap_len = gap_ops.iter().map(|(_, n)| *n).max().unwrap_or(0);
        assert!(max_gap_len > 1, "expected gap length > 1, got {}, ops={:?}", max_gap_len, esp.ops);
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
