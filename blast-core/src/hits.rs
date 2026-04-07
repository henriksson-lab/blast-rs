//! Rust equivalent of blast_hits.c — HSP list management and filtering.

use crate::hspstream::{Hsp, HspList};


/// Compute identities and alignment length from query/subject sequences
/// and a gap edit script. Used after traceback to populate HSP stats.
pub fn compute_identities(
    query: &[u8],
    subject: &[u8],
    q_offset: usize,
    s_offset: usize,
    gap_ops: &[(u32, i32)], // (op_type, count) pairs from GapEditScript
) -> (i32, i32, i32) {
    // op_type: 0=Del, 3=Sub, 4+=Ins
    let mut q_pos = q_offset;
    let mut s_pos = s_offset;
    let mut align_len = 0i32;
    let mut num_ident = 0i32;
    let mut gap_opens = 0i32;

    for &(op, count) in gap_ops {
        let count = count as usize;
        align_len += count as i32;
        match op {
            3 => {
                // Substitution
                for _ in 0..count {
                    if q_pos < query.len() && s_pos < subject.len()
                        && query[q_pos] == subject[s_pos]
                    {
                        num_ident += 1;
                    }
                    q_pos += 1;
                    s_pos += 1;
                }
            }
            0 | 1 | 2 => {
                // Deletion (gap in query, advance subject)
                s_pos += count;
                gap_opens += 1;
            }
            4 | 5 | 6 => {
                // Insertion (gap in subject, advance query)
                q_pos += count;
                gap_opens += 1;
            }
            _ => {
                q_pos += count;
                s_pos += count;
            }
        }
    }
    (align_len, num_ident, gap_opens)
}

/// Filter an HSP list by e-value threshold.
pub fn filter_by_evalue(list: &mut HspList, max_evalue: f64) {
    list.hsps.retain(|hsp| hsp.evalue <= max_evalue);
}

/// Sort HSP list by score (descending).
pub fn sort_by_score(list: &mut HspList) {
    list.hsps.sort_by(|a, b| b.score.cmp(&a.score));
}

/// Remove HSPs that are contained within higher-scoring HSPs.
/// Uses a simple O(n^2) containment check.
pub fn remove_contained(list: &mut HspList) {
    if list.hsps.len() <= 1 {
        return;
    }
    // Sort by score first
    list.hsps.sort_by(|a, b| b.score.cmp(&a.score));

    let mut keep = vec![true; list.hsps.len()];
    for i in 0..list.hsps.len() {
        if !keep[i] { continue; }
        for j in (i + 1)..list.hsps.len() {
            if !keep[j] { continue; }
            // Check if j is contained within i on both query and subject
            if list.hsps[j].query_offset >= list.hsps[i].query_offset
                && list.hsps[j].query_end <= list.hsps[i].query_end
                && list.hsps[j].subject_offset >= list.hsps[i].subject_offset
                && list.hsps[j].subject_end <= list.hsps[i].subject_end
            {
                keep[j] = false;
            }
        }
    }

    let mut idx = 0;
    list.hsps.retain(|_| {
        let k = keep[idx];
        idx += 1;
        k
    });
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compute_identities_perfect() {
        let query = vec![0u8, 1, 2, 3, 0, 1];
        let subject = vec![0u8, 1, 2, 3, 0, 1];
        let ops = vec![(3, 6)]; // 6 substitutions
        let (len, ident, gaps) = compute_identities(&query, &subject, 0, 0, &ops);
        assert_eq!(len, 6);
        assert_eq!(ident, 6);
        assert_eq!(gaps, 0);
    }

    #[test]
    fn test_compute_identities_with_gap() {
        let query = vec![0u8, 1, 2, 3];
        let subject = vec![0u8, 1, 9, 9, 2, 3];
        let ops = vec![(3, 2), (0, 2), (3, 2)]; // 2 sub, 2 del, 2 sub
        let (len, ident, gaps) = compute_identities(&query, &subject, 0, 0, &ops);
        assert_eq!(len, 6);
        assert_eq!(ident, 4);
        assert_eq!(gaps, 1);
    }

    #[test]
    fn test_filter_by_evalue() {
        let mut list = HspList::new(0);
        list.add_hsp(Hsp {
            score: 100, num_ident: 50, bit_score: 91.5, evalue: 1e-20,
            query_offset: 0, query_end: 50, subject_offset: 0, subject_end: 50,
            context: 0, num_gaps: 0,
        });
        list.add_hsp(Hsp {
            score: 20, num_ident: 10, bit_score: 22.3, evalue: 5.0,
            query_offset: 0, query_end: 10, subject_offset: 0, subject_end: 10,
            context: 0, num_gaps: 0,
        });
        filter_by_evalue(&mut list, 0.001);
        assert_eq!(list.hsps.len(), 1);
        assert_eq!(list.hsps[0].score, 100);
    }
}
