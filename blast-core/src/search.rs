//! Pure Rust BLAST search engine — no FFI.
//! This module implements the complete blastn search pipeline in Rust.

use crate::encoding;
use crate::gapinfo::{GapAlignOpType, GapEditScript};
use crate::stat::KarlinBlk;

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

    // Search both strands
    for (context, query) in [(0i32, query_plus), (1i32, query_minus)] {
        if query.len() < word_size || subject.len() < word_size {
            continue;
        }

        // Build lookup table: for word sizes up to 13, use direct array (4^w entries)
        // For larger word sizes, use the first 13 bases as key
        let lut_word = word_size.min(13);
        let lut_size = 1usize << (2 * lut_word); // 4^lut_word
        let mut lut: Vec<Vec<u32>> = vec![Vec::new(); lut_size];

        for i in 0..=(query.len() - word_size) {
            let key = word_hash_n(&query[i..i + lut_word], lut_word) as usize;
            lut[key].push(i as u32);
        }

        // Scan subject using the lookup table
        for s_pos in 0..=(subject.len() - word_size) {
            let key = word_hash_n(&subject[s_pos..s_pos + lut_word], lut_word) as usize;
            if lut[key].is_empty() { continue; }

            for &q_pos in &lut[key] {
                let q_pos = q_pos as usize;
                // For word_size > lut_word, verify remaining bases
                if word_size > lut_word {
                    if query[q_pos + lut_word..q_pos + word_size]
                        != subject[s_pos + lut_word..s_pos + word_size] {
                        continue;
                    }
                }
                // Extend the seed
                if let Some(hsp) = extend_seed(
                    query, subject, q_pos, s_pos,
                    reward, penalty, x_dropoff,
                    kbp, search_space, evalue_threshold, context,
                ) {
                    hsps.push(hsp);
                }
            }
        }
    }

    // Sort by score descending
    hsps.sort_by(|a, b| b.score.cmp(&a.score));

    // Remove overlapping HSPs (simple dedup)
    dedup_hsps(&mut hsps);

    hsps
}

/// Hash a word (sequence of BLASTNA values) for the lookup table.
fn word_hash(word: &[u8]) -> u64 {
    let mut h = 0u64;
    for &b in word {
        h = (h << 2) | (b & 3) as u64;
    }
    h
}

/// Extend a seed hit in both directions using ungapped extension.
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

/// Remove HSPs that overlap significantly with higher-scoring ones.
fn dedup_hsps(hsps: &mut Vec<SearchHsp>) {
    if hsps.len() <= 1 {
        return;
    }
    let mut keep = vec![true; hsps.len()];
    for i in 0..hsps.len() {
        if !keep[i] { continue; }
        for j in (i + 1)..hsps.len() {
            if !keep[j] { continue; }
            if hsps[i].context != hsps[j].context { continue; }
            // Check if j's query range overlaps significantly with i
            let q_overlap = (hsps[i].query_end.min(hsps[j].query_end)
                - hsps[i].query_start.max(hsps[j].query_start)).max(0);
            let q_len = hsps[j].query_end - hsps[j].query_start;
            if q_len > 0 && q_overlap as f64 / q_len as f64 > 0.5 {
                keep[j] = false;
            }
        }
    }
    let mut idx = 0;
    hsps.retain(|_| { let k = keep[idx]; idx += 1; k });
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
        assert_eq!(word_hash(&[0, 1, 2, 3]), 0b00_01_10_11);
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
