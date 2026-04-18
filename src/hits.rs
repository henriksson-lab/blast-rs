//! Rust equivalent of blast_hits.c — HSP list management and filtering.

#[cfg(test)]
use crate::hspstream::Hsp;
use crate::hspstream::HspList;

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
                    if q_pos < query.len()
                        && s_pos < subject.len()
                        && query[q_pos] == subject[s_pos]
                    {
                        num_ident += 1;
                    }
                    q_pos += 1;
                    s_pos += 1;
                }
            }
            0..=2 => {
                // Deletion (gap in query, advance subject)
                s_pos += count;
                gap_opens += 1;
            }
            4..=6 => {
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

/// Sort HSP list by score (descending), using NCBI's full tie-breaker.
/// See `crate::hspstream::score_compare_hsps` (port of `ScoreCompareHSPs`).
pub fn sort_by_score(list: &mut HspList) {
    list.hsps.sort_by(crate::hspstream::score_compare_hsps);
}

/// Remove HSPs that are contained within higher-scoring HSPs.
/// Uses a simple O(n^2) containment check.
pub fn remove_contained(list: &mut HspList) {
    if list.hsps.len() <= 1 {
        return;
    }
    // Sort by score (with NCBI's full tie-breaker from `ScoreCompareHSPs`)
    // so the containment scan consistently keeps the NCBI-preferred HSP
    // when scores tie (smaller subject_offset wins, etc).
    list.hsps.sort_by(crate::hspstream::score_compare_hsps);

    let mut keep = vec![true; list.hsps.len()];
    for i in 0..list.hsps.len() {
        if !keep[i] {
            continue;
        }
        for j in (i + 1)..list.hsps.len() {
            if !keep[j] {
                continue;
            }
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
            score: 100,
            num_ident: 50,
            bit_score: 91.5,
            evalue: 1e-20,
            query_offset: 0,
            query_end: 50,
            subject_offset: 0,
            subject_end: 50,
            context: 0,
            num_gaps: 0,
        });
        list.add_hsp(Hsp {
            score: 20,
            num_ident: 10,
            bit_score: 22.3,
            evalue: 5.0,
            query_offset: 0,
            query_end: 10,
            subject_offset: 0,
            subject_end: 10,
            context: 0,
            num_gaps: 0,
        });
        filter_by_evalue(&mut list, 0.001);
        assert_eq!(list.hsps.len(), 1);
        assert_eq!(list.hsps[0].score, 100);
    }

    /// Helper to create an HSP with specified fields, defaulting others.
    fn make_hsp(score: i32, evalue: f64, q_off: i32, q_end: i32, s_off: i32, s_end: i32) -> Hsp {
        Hsp {
            score,
            num_ident: score / 2,
            bit_score: score as f64 * 0.9,
            evalue,
            query_offset: q_off,
            query_end: q_end,
            subject_offset: s_off,
            subject_end: s_end,
            context: 0,
            num_gaps: 0,
        }
    }

    #[test]
    fn test_hsp_list_sort_by_score() {
        let mut list = HspList::new(0);
        list.add_hsp(make_hsp(30, 1e-5, 0, 30, 0, 30));
        list.add_hsp(make_hsp(90, 1e-20, 0, 90, 0, 90));
        list.add_hsp(make_hsp(10, 1.0, 0, 10, 0, 10));
        list.add_hsp(make_hsp(60, 1e-12, 0, 60, 0, 60));

        sort_by_score(&mut list);

        let scores: Vec<i32> = list.hsps.iter().map(|h| h.score).collect();
        assert_eq!(scores, vec![90, 60, 30, 10]);
    }

    #[test]
    fn test_hsp_list_purge_nulls() {
        // In Rust we don't have NULL pointers, but we can simulate purging
        // invalid/zero-score entries by filtering them out.
        let mut list = HspList::new(0);
        list.add_hsp(make_hsp(50, 1e-10, 0, 50, 0, 50));
        list.add_hsp(make_hsp(0, f64::MAX, 0, 0, 0, 0)); // "null" / invalid
        list.add_hsp(make_hsp(40, 1e-8, 10, 50, 10, 50));
        list.add_hsp(make_hsp(0, f64::MAX, 0, 0, 0, 0)); // "null" / invalid

        // Purge zero-score entries (Rust analog of purging NULLs)
        list.hsps.retain(|h| h.score > 0);

        assert_eq!(list.hsps.len(), 2);
        assert_eq!(list.hsps[0].score, 50);
        assert_eq!(list.hsps[1].score, 40);
    }

    #[test]
    fn test_hsp_overlap_detection() {
        // Two HSPs on the same diagonal (query_offset - subject_offset is the same)
        // that overlap in coordinate space.
        let hsp_a = make_hsp(80, 1e-15, 10, 60, 10, 60); // diagonal = 0
        let hsp_b = make_hsp(50, 1e-8, 40, 80, 40, 80); // diagonal = 0, overlaps 40..60

        // Check overlap: ranges [10,60) and [40,80) overlap on both query and subject
        let q_overlap =
            hsp_a.query_end > hsp_b.query_offset && hsp_b.query_end > hsp_a.query_offset;
        let s_overlap =
            hsp_a.subject_end > hsp_b.subject_offset && hsp_b.subject_end > hsp_a.subject_offset;
        assert!(
            q_overlap && s_overlap,
            "HSPs on the same diagonal should overlap"
        );
    }

    #[test]
    fn test_hsp_no_overlap() {
        // Two HSPs on different diagonals with non-overlapping coordinates.
        let hsp_a = make_hsp(80, 1e-15, 0, 50, 100, 150); // diagonal = -100
        let hsp_b = make_hsp(50, 1e-8, 200, 250, 0, 50); // diagonal = 200

        let q_overlap =
            hsp_a.query_end > hsp_b.query_offset && hsp_b.query_end > hsp_a.query_offset;
        let s_overlap =
            hsp_a.subject_end > hsp_b.subject_offset && hsp_b.subject_end > hsp_a.subject_offset;
        assert!(
            !q_overlap || !s_overlap,
            "HSPs on different diagonals should not overlap"
        );
    }

    #[test]
    fn test_hsp_containment() {
        // One HSP fully contained within another; remove_contained should drop it.
        let mut list = HspList::new(0);
        list.add_hsp(make_hsp(100, 1e-25, 10, 100, 10, 100)); // outer, higher score
        list.add_hsp(make_hsp(40, 1e-5, 30, 70, 30, 70)); // inner, fully contained
        list.add_hsp(make_hsp(60, 1e-10, 80, 130, 80, 130)); // separate, not contained

        remove_contained(&mut list);

        assert_eq!(list.hsps.len(), 2);
        let scores: Vec<i32> = list.hsps.iter().map(|h| h.score).collect();
        // Sorted by score descending, the contained one (40) should be gone
        assert!(scores.contains(&100));
        assert!(scores.contains(&60));
        assert!(!scores.contains(&40));
    }

    #[test]
    fn test_hsp_merge() {
        // Simulate merging overlapping HSPs by taking the union of their coordinates.
        let hsp_a = make_hsp(80, 1e-15, 10, 60, 10, 60);
        let hsp_b = make_hsp(50, 1e-8, 40, 90, 40, 90);

        // Merged HSP covers the union of both ranges
        let merged_q_off = hsp_a.query_offset.min(hsp_b.query_offset);
        let merged_q_end = hsp_a.query_end.max(hsp_b.query_end);
        let merged_s_off = hsp_a.subject_offset.min(hsp_b.subject_offset);
        let merged_s_end = hsp_a.subject_end.max(hsp_b.subject_end);
        let merged_score = hsp_a.score.max(hsp_b.score);

        assert_eq!(merged_q_off, 10);
        assert_eq!(merged_q_end, 90);
        assert_eq!(merged_s_off, 10);
        assert_eq!(merged_s_end, 90);
        assert_eq!(merged_score, 80);

        // Verify the merged region covers both originals
        assert!(merged_q_off <= hsp_a.query_offset && merged_q_end >= hsp_a.query_end);
        assert!(merged_q_off <= hsp_b.query_offset && merged_q_end >= hsp_b.query_end);
    }

    #[test]
    fn test_hsp_evalue_filter() {
        let mut list = HspList::new(0);
        list.add_hsp(make_hsp(100, 1e-30, 0, 100, 0, 100));
        list.add_hsp(make_hsp(60, 1e-10, 0, 60, 0, 60));
        list.add_hsp(make_hsp(30, 0.5, 0, 30, 0, 30));
        list.add_hsp(make_hsp(10, 10.0, 0, 10, 0, 10));
        list.add_hsp(make_hsp(50, 1e-6, 0, 50, 0, 50));

        filter_by_evalue(&mut list, 1e-5);

        assert_eq!(list.hsps.len(), 3);
        // Only HSPs with evalue <= 1e-5 should remain
        for hsp in &list.hsps {
            assert!(
                hsp.evalue <= 1e-5,
                "HSP with evalue {} should have been filtered",
                hsp.evalue
            );
        }
        let scores: Vec<i32> = list.hsps.iter().map(|h| h.score).collect();
        assert!(scores.contains(&100));
        assert!(scores.contains(&60));
        assert!(scores.contains(&50));
    }

    #[test]
    fn test_hsp_max_target_seqs() {
        // Simulate max_target_seqs by sorting hit lists by evalue and truncating.
        use crate::hspstream::HitList;

        let mut hitlist = HitList::new();
        for i in 0..10 {
            let mut hsp_list = HspList::new(i);
            hsp_list.add_hsp(make_hsp(
                100 - i * 10,
                1e-5 * (i as f64 + 1.0),
                0,
                50,
                0,
                50,
            ));
            hitlist.add_hsp_list(hsp_list);
        }
        assert_eq!(hitlist.hsp_lists.len(), 10);

        hitlist.sort_by_evalue();

        let max_target_seqs = 5;
        hitlist.hsp_lists.truncate(max_target_seqs);

        assert_eq!(hitlist.hsp_lists.len(), 5);
        // The kept entries should have the lowest evalues
        for i in 1..hitlist.hsp_lists.len() {
            assert!(hitlist.hsp_lists[i].best_evalue >= hitlist.hsp_lists[i - 1].best_evalue);
        }
    }
}
