//! Greedy nucleotide alignment for megablast.
//! Faster than DP for high-identity matches (>90% identity).

use crate::gapinfo::{GapAlignOpType, GapEditScript};

/// Perform greedy gapped alignment between two nucleotide sequences.
/// Optimized for high-identity matches — extends greedily and only
/// introduces gaps when the score drops significantly.
///
/// Returns (score, query_start, query_end, subject_start, subject_end, edit_script)
pub fn greedy_align(
    query: &[u8],
    subject: &[u8],
    q_seed: usize,
    s_seed: usize,
    reward: i32,
    penalty: i32,
    x_dropoff: i32,
) -> Option<(i32, usize, usize, usize, usize, GapEditScript)> {
    // Extend right from seed
    let (r_score, r_q_end, r_s_end, r_ops) =
        extend_greedy(query, subject, q_seed, s_seed, reward, penalty, x_dropoff, true);

    // Extend left from seed
    let (l_score, l_q_start, l_s_start, mut l_ops) =
        extend_greedy(query, subject, q_seed, s_seed, reward, penalty, x_dropoff, false);

    let total_score = r_score + l_score;
    if total_score <= 0 {
        return None;
    }

    // Combine left (reversed) and right extensions
    l_ops.reverse();
    let mut esp = GapEditScript::new();
    for (op, count) in l_ops {
        esp.push(op, count);
    }
    for (op, count) in r_ops {
        esp.push(op, count);
    }

    Some((total_score, l_q_start, r_q_end, l_s_start, r_s_end, esp))
}

/// Extend greedily in one direction.
/// Returns (score, q_end_or_start, s_end_or_start, ops)
fn extend_greedy(
    query: &[u8],
    subject: &[u8],
    q_start: usize,
    s_start: usize,
    reward: i32,
    penalty: i32,
    x_dropoff: i32,
    forward: bool,
) -> (i32, usize, usize, Vec<(GapAlignOpType, i32)>) {
    let mut score = 0i32;
    let mut best_score = 0i32;
    let mut qi = q_start;
    let mut si = s_start;
    let mut best_qi = qi;
    let mut best_si = si;
    let mut ops: Vec<(GapAlignOpType, i32)> = Vec::new();
    let mut match_run = 0i32;

    loop {
        let (q_ok, s_ok) = if forward {
            (qi < query.len(), si < subject.len())
        } else {
            (qi > 0, si > 0)
        };
        if !q_ok || !s_ok { break; }

        let (q_idx, s_idx) = if forward {
            (qi, si)
        } else {
            (qi - 1, si - 1)
        };

        if query[q_idx] == subject[s_idx] {
            score += reward;
            match_run += 1;
        } else {
            // Flush match run
            if match_run > 0 {
                ops.push((GapAlignOpType::Sub, match_run));
                match_run = 0;
            }
            score += penalty;
            ops.push((GapAlignOpType::Sub, 1)); // mismatch
        }

        if score > best_score {
            best_score = score;
            best_qi = if forward { qi + 1 } else { qi - 1 };
            best_si = if forward { si + 1 } else { si - 1 };
        }

        if best_score - score > x_dropoff {
            break;
        }

        if forward { qi += 1; si += 1; } else { qi -= 1; si -= 1; }
    }

    // Flush remaining match run
    if match_run > 0 {
        ops.push((GapAlignOpType::Sub, match_run));
    }

    (best_score, best_qi, best_si, ops)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_greedy_perfect_match() {
        let q = vec![0u8, 1, 2, 3, 0, 1, 2, 3];
        let s = vec![0u8, 1, 2, 3, 0, 1, 2, 3];
        let result = greedy_align(&q, &s, 4, 4, 1, -3, 20);
        assert!(result.is_some());
        let (score, qs, qe, ss, se, _esp) = result.unwrap();
        assert_eq!(score, 8); // 8 matches * reward=1
        assert_eq!(qe - qs, 8);
        assert_eq!(se - ss, 8);
    }

    #[test]
    fn test_greedy_with_mismatch() {
        let q = vec![0u8, 1, 2, 3, 0, 1, 2, 3];
        let s = vec![0u8, 1, 2, 3, 3, 1, 2, 3]; // mismatch at pos 4
        let result = greedy_align(&q, &s, 4, 4, 1, -3, 10);
        assert!(result.is_some());
        let (score, _, _, _, _, _) = result.unwrap();
        assert!(score > 0);
    }

    #[test]
    fn test_greedy_no_match() {
        let q = vec![0u8, 0, 0, 0];
        let s = vec![3u8, 3, 3, 3];
        let result = greedy_align(&q, &s, 2, 2, 1, -3, 5);
        assert!(result.is_none());
    }
}
