//! Pure Rust BLAST search engine — no FFI.
//! This module implements the complete blastn search pipeline in Rust.


use crate::itree::{IntervalTree, Interval};
use crate::stat::KarlinBlk;
use crate::traceback::traceback_align_abs;


/// Result of a single HSP (High-Scoring Pair).
#[derive(Debug, Clone)]
pub struct SearchHsp {
    pub query_start: i32,
    pub query_end: i32,
    pub subject_start: i32,
    pub subject_end: i32,
    pub score: i32,
    pub bit_score: f64,
    pub evalue: f64,
    pub num_ident: i32,
    pub align_length: i32,
    pub mismatches: i32,
    pub gap_opens: i32,
    pub context: i32, // 0=plus, 1=minus
}

/// Result of searching one query against one subject.
#[derive(Debug)]
pub struct SubjectResult {
    pub oid: i32,
    pub hsps: Vec<SearchHsp>,
}

/// Perform a simple ungapped nucleotide word search.
/// Finds exact word matches between query and subject, extends them,
/// and returns significant HSPs.
pub fn blastn_ungapped_search(
    query_plus: &[u8],     // BLASTNA encoded, plus strand
    query_minus: &[u8],    // BLASTNA encoded, minus strand (RC)
    subject: &[u8],        // BLASTNA decoded subject
    word_size: usize,
    reward: i32,
    penalty: i32,
    x_dropoff: i32,
    kbp: &KarlinBlk,
    search_space: f64,
    evalue_threshold: f64,
) -> Vec<SearchHsp> {
    let mut hsps = Vec::new();

    // Use a small word (8 bases = 4^8 = 65536 entries) for the lookup table
    let lut_word = word_size.min(8);
    let lut_size = 1usize << (2 * lut_word);

    // Search both strands
    for (context, query) in [(0i32, query_plus), (1i32, query_minus)] {
        if query.len() < word_size || subject.len() < word_size {
            continue;
        }

        // Build lookup table from query (once per strand)
        let mut lut: Vec<i32> = vec![-1; lut_size]; // -1 = empty
        let mut next: Vec<i32> = vec![-1; query.len()]; // chain for collisions

        for i in (0..=(query.len() - word_size)).rev() {
            let key = word_hash_n(&query[i..i + lut_word], lut_word) as usize;
            next[i] = lut[key];
            lut[key] = i as i32;
        }

        // Scan subject
        for s_pos in 0..=(subject.len() - word_size) {
            let key = word_hash_n(&subject[s_pos..s_pos + lut_word], lut_word) as usize;
            let mut q_pos = lut[key];
            while q_pos >= 0 {
                let qp = q_pos as usize;
                // Verify remaining bases if word_size > lut_word
                let mut matches = true;
                if word_size > lut_word {
                    for k in lut_word..word_size {
                        if query[qp + k] != subject[s_pos + k] {
                            matches = false;
                            break;
                        }
                    }
                }
                if matches {
                    if let Some(hsp) = extend_seed(
                        query, subject, qp, s_pos,
                        reward, penalty, x_dropoff,
                        kbp, search_space, evalue_threshold, context,
                    ) {
                        hsps.push(hsp);
                    }
                }
                q_pos = next[qp];
            }
        }
    }

    // Sort by score descending
    hsps.sort_by(|a, b| b.score.cmp(&a.score));

    // Remove overlapping HSPs (simple dedup)
    dedup_hsps(&mut hsps);

    hsps
}

/// Hash the first n bases of a word for the lookup table.
#[inline(always)]
fn word_hash_n(word: &[u8], n: usize) -> u32 {
    let mut h = 0u32;
    for i in 0..n {
        h = (h << 2) | (word[i] & 3) as u32;
    }
    h
}

/// Extend a seed hit in both directions using ungapped extension.
#[inline]
fn extend_seed(
    query: &[u8],
    subject: &[u8],
    q_seed: usize,
    s_seed: usize,
    reward: i32,
    penalty: i32,
    x_dropoff: i32,
    kbp: &KarlinBlk,
    search_space: f64,
    evalue_threshold: f64,
    context: i32,
) -> Option<SearchHsp> {
    // Extend right
    let mut score = 0i32;
    let mut best_score = 0i32;
    let mut best_right = 0usize;
    let mut qi = q_seed;
    let mut si = s_seed;
    while qi < query.len() && si < subject.len() {
        score += if query[qi] == subject[si] { reward } else { penalty };
        if score > best_score {
            best_score = score;
            best_right = qi - q_seed + 1;
        }
        if best_score - score > x_dropoff {
            break;
        }
        qi += 1;
        si += 1;
    }

    // Extend left
    let mut score_l = 0i32;
    let mut best_score_l = 0i32;
    let mut best_left = 0usize;
    if q_seed > 0 && s_seed > 0 {
        qi = q_seed - 1;
        si = s_seed - 1;
        loop {
            score_l += if query[qi] == subject[si] { reward } else { penalty };
            if score_l > best_score_l {
                best_score_l = score_l;
                best_left = q_seed - qi;
            }
            if best_score_l - score_l > x_dropoff {
                break;
            }
            if qi == 0 || si == 0 {
                break;
            }
            qi -= 1;
            si -= 1;
        }
    }

    let total_score = best_score + best_score_l;
    if total_score <= 0 {
        return None;
    }

    // Compute statistics
    let evalue = kbp.raw_to_evalue(total_score, search_space);
    if evalue > evalue_threshold {
        return None;
    }
    let bit_score = kbp.raw_to_bit(total_score);

    let q_start = (q_seed - best_left) as i32;
    let q_end = (q_seed + best_right) as i32;
    let s_start = (s_seed - best_left) as i32;
    let s_end = (s_seed + best_right) as i32;
    let align_len = q_end - q_start;

    // Count identities
    let mut num_ident = 0;
    for i in 0..align_len as usize {
        if query[q_start as usize + i] == subject[s_start as usize + i] {
            num_ident += 1;
        }
    }

    Some(SearchHsp {
        query_start: q_start,
        query_end: q_end,
        subject_start: s_start,
        subject_end: s_end,
        score: total_score,
        bit_score,
        evalue,
        num_ident,
        align_length: align_len,
        mismatches: align_len - num_ident,
        gap_opens: 0,
        context,
    })
}

/// Remove HSPs that are contained within or significantly overlap higher-scoring ones.
/// Uses an interval tree for efficient containment checking.
fn dedup_hsps(hsps: &mut Vec<SearchHsp>) {
    if hsps.len() <= 1 {
        return;
    }
    hsps.sort_by(|a, b| b.score.cmp(&a.score));

    // Separate trees per context (strand)
    let q_max = hsps.iter().map(|h| h.query_end).max().unwrap_or(0) + 1;
    let s_max = hsps.iter().map(|h| h.subject_end).max().unwrap_or(0) + 1;
    let mut trees: std::collections::HashMap<i32, IntervalTree> = std::collections::HashMap::new();

    let mut keep = vec![false; hsps.len()];
    for (i, hsp) in hsps.iter().enumerate() {
        let tree = trees.entry(hsp.context).or_insert_with(|| IntervalTree::new(q_max, s_max));
        let qi = Interval::new(hsp.query_start, hsp.query_end);
        let si = Interval::new(hsp.subject_start, hsp.subject_end);
        if !tree.is_contained(qi, si) && !tree.has_significant_overlap(qi, si, 0.5) {
            keep[i] = true;
            tree.insert(qi, si, hsp.score);
        }
    }
    let mut idx = 0;
    hsps.retain(|_| { let k = keep[idx]; idx += 1; k });
}

/// Perform gapped blastn search with traceback.
/// First finds seeds via ungapped scanning, then does full gapped alignment
/// on seed hits that pass the ungapped score threshold.
pub fn blastn_gapped_search(
    query_plus: &[u8],
    query_minus: &[u8],
    subject: &[u8],
    word_size: usize,
    reward: i32,
    penalty: i32,
    _gap_open: i32,
    _gap_extend: i32,
    x_dropoff: i32,
    kbp: &KarlinBlk,
    search_space: f64,
    evalue_threshold: f64,
) -> Vec<SearchHsp> {
    // First do ungapped search to find seeds
    let ungapped = blastn_ungapped_search(
        query_plus, query_minus, subject,
        word_size, reward, penalty, x_dropoff,
        kbp, search_space, evalue_threshold * 100.0, // permissive threshold for seeds
    );

    let mut hsps = Vec::new();

    for seed in &ungapped {
        let query = if seed.context == 0 { query_plus } else { query_minus };

        // Determine the region to align (extend seed boundaries)
        let q_margin = 50.min(seed.query_start as usize);
        let s_margin = 50.min(seed.subject_start as usize);
        let q_start = seed.query_start as usize - q_margin;
        let s_start = seed.subject_start as usize - s_margin;
        let q_end = (seed.query_end as usize + 50).min(query.len());
        let s_end = (seed.subject_end as usize + 50).min(subject.len());

        if q_end <= q_start || s_end <= s_start {
            // Use the ungapped hit as-is
            hsps.push(seed.clone());
            continue;
        }

        // Do gapped traceback around the seed
        let seed_q = ((seed.query_start + seed.query_end) / 2) as usize;
        let seed_s = ((seed.subject_start + seed.subject_end) / 2) as usize;
        let margin = (seed.align_length as usize).max(50);

        if let Some(tb) = traceback_align_abs(
            query, subject, seed_q, seed_s,
            reward, penalty, _gap_open, _gap_extend, margin,
        ) {
            let evalue = kbp.raw_to_evalue(tb.score, search_space);
            if evalue > evalue_threshold { continue; }

            let (align_len, num_ident, gap_opens) = tb.edit_script.count_identities(
                &query[tb.query_start..tb.query_end],
                &subject[tb.subject_start..tb.subject_end],
            );

            hsps.push(SearchHsp {
                query_start: tb.query_start as i32,
                query_end: tb.query_end as i32,
                subject_start: tb.subject_start as i32,
                subject_end: tb.subject_end as i32,
                score: tb.score,
                bit_score: kbp.raw_to_bit(tb.score),
                evalue,
                num_ident,
                align_length: align_len,
                mismatches: (align_len - num_ident - gap_opens).max(0),
                gap_opens,
                context: seed.context,
            });
        } else {
            // Traceback failed, use ungapped seed
            hsps.push(seed.clone());
        }
    }

    hsps.sort_by(|a, b| b.score.cmp(&a.score));
    dedup_hsps(&mut hsps);
    hsps
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_kbp() -> KarlinBlk {
        KarlinBlk {
            lambda: 0.208,
            k: 0.049,
            log_k: 0.049_f64.ln(),
            h: 0.14,
        }
    }

    #[test]
    fn test_word_hash() {
        // ACGT = 0,1,2,3
        assert_eq!(word_hash_n(&[0, 1, 2, 3], 4), 0b00_01_10_11);
    }

    #[test]
    fn test_simple_search() {
        // Query: ACGTACGTACGT (12 bases)
        let query = vec![0u8, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
        let rc: Vec<u8> = query.iter().rev().map(|&b| 3 - b).collect();
        // Subject contains exact match at position 5
        let mut subject = vec![3u8; 30]; // TTTTTT...
        for (i, &b) in query.iter().enumerate() {
            subject[5 + i] = b;
        }

        let kbp = test_kbp();
        let results = blastn_ungapped_search(
            &query, &rc, &subject, 7, 2, -3, 20,
            &kbp, 1e6, 1e10, // very permissive e-value for testing
        );

        assert!(!results.is_empty(), "Should find at least one hit");
        assert_eq!(results[0].subject_start, 5);
        assert_eq!(results[0].subject_end, 17);
        assert_eq!(results[0].num_ident, 12);
    }

    #[test]
    fn test_ungapped_mismatch() {
        // Query with 2 mismatches should still find the hit
        let mut query = vec![0u8, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
        query[4] = 3; // mismatch at position 4
        query[8] = 1; // mismatch at position 8
        let rc: Vec<u8> = query.iter().rev().map(|&b| 3 - b).collect();
        let mut subject = vec![3u8; 30];
        let original = vec![0u8, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
        for (i, &b) in original.iter().enumerate() {
            subject[5 + i] = b;
        }

        let kbp = test_kbp();
        let results = blastn_ungapped_search(
            &query, &rc, &subject, 4, 2, -3, 20,
            &kbp, 1e6, 1e10,
        );
        assert!(!results.is_empty(), "Should find hit despite mismatches (word_size=4)");
        // The hit may be partial (avoiding the mismatch region)
        assert!(results[0].score > 0);
    }

    #[test]
    fn test_ungapped_both_strands() {
        // Subject contains ACGTACGTACGT at position 10
        // Query is the reverse complement — should match via the RC (context 1)
        let subject_insert = vec![0u8, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
        let mut subject = vec![3u8; 30];
        for (i, &b) in subject_insert.iter().enumerate() {
            subject[10 + i] = b;
        }
        // Query is RC of the insert: ACGTACGTACGT RC = 3,2,1,0,3,2,1,0,3,2,1,0
        let rc_query: Vec<u8> = subject_insert.iter().rev().map(|&b| 3 - b).collect();
        // The plus strand query won't match the subject, but the minus strand (rc of rc = original) will
        let query_plus = rc_query.clone(); // this is the RC
        let query_minus: Vec<u8> = query_plus.iter().rev().map(|&b| 3 - b).collect(); // RC of RC = original

        let kbp = test_kbp();
        let results = blastn_ungapped_search(
            &query_plus, &query_minus, &subject, 7, 2, -3, 20,
            &kbp, 1e6, 1e10,
        );
        assert!(!results.is_empty(), "Should find hit on minus strand");
        // The hit should be found via the minus strand query (context=1)
        assert!(results.iter().any(|h| h.context == 1), "Should have a minus-strand hit");
    }

    #[test]
    fn test_gapped_search() {
        let query = vec![0u8, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
        let rc: Vec<u8> = query.iter().rev().map(|&b| 3 - b).collect();
        let mut subject = vec![3u8; 30];
        for (i, &b) in query.iter().enumerate() {
            subject[5 + i] = b;
        }

        let kbp = test_kbp();
        let results = blastn_gapped_search(
            &query, &rc, &subject, 7, 2, -3, 5, 2, 20,
            &kbp, 1e6, 1e10,
        );
        assert!(!results.is_empty(), "Gapped search should find hit");
        assert_eq!(results[0].gap_opens, 0, "Perfect match should have no gaps");
    }

    #[test]
    fn test_no_match() {
        let query = vec![0u8; 20]; // AAAAAAA...
        let rc = vec![3u8; 20];    // TTTTTTT...
        let subject = vec![1u8; 50]; // CCCCCCC...

        let kbp = test_kbp();
        let results = blastn_ungapped_search(
            &query, &rc, &subject, 11, 1, -3, 20,
            &kbp, 1e6, 10.0,
        );
        assert!(results.is_empty(), "Should find no hits");
    }
}
