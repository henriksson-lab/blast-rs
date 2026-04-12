//! Pure Rust BLAST search engine — no FFI.
//! This module implements the complete blastn search pipeline in Rust.


use crate::itree::{IntervalTree, Interval};
use crate::sequence::blastna_to_iupac;
use crate::stat::KarlinBlk;
use crate::traceback::{blast_gapped_align, blast_gapped_score_only};


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
    pub qseq: Option<String>,
    pub sseq: Option<String>,
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

    let lut_word = word_size.min(8);
    let lut_size = 1usize << (2 * lut_word);

    for (context, query) in [(0i32, query_plus), (1i32, query_minus)] {
        if query.len() < word_size || subject.len() < word_size {
            continue;
        }

        // Build lookup table from query
        let mut lut: Vec<i32> = vec![-1; lut_size];
        let mut next: Vec<i32> = vec![-1; query.len()];

        for i in (0..=(query.len() - word_size)).rev() {
            let key = word_hash_n(&query[i..i + lut_word], lut_word) as usize;
            next[i] = lut[key];
            lut[key] = i as i32;
        }

        // Scan subject
        let mut s_pos = 0;
        while s_pos + word_size <= subject.len() {
            let key = word_hash_n(&subject[s_pos..s_pos + lut_word], lut_word) as usize;
            let mut q_pos = lut[key];
            while q_pos >= 0 {
                let qp = q_pos as usize;
                let mut matches = true;
                if word_size > lut_word {
                    for k in lut_word..word_size {
                        if qp + k >= query.len() || s_pos + k >= subject.len() {
                            matches = false;
                            break;
                        }
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
            s_pos += 1;
        }
    }

    hsps.sort_by(|a, b| b.score.cmp(&a.score));
    dedup_hsps(&mut hsps);
    hsps
}

/// Process a seed hit from the scan loop (verify word, diagonal check, extend).
#[inline]
fn process_hit(
    query: &[u8], subject_packed: &[u8], subject_len: usize, word_size: usize,
    lut: &[i32], next: &[i32], last_hit: &mut [i32], diag_mask: usize,
    h: usize, base_pos: usize,
    reward: i32, penalty: i32, x_dropoff: i32,
    kbp: &KarlinBlk, search_space: f64, evalue_threshold: f64,
    context: i32, hsps: &mut Vec<SearchHsp>,
) {
    let sp = base_pos;
    let mut q_pos = lut[h];
    while q_pos >= 0 {
        let qp = q_pos as usize;
        let extra = if word_size > 8 { word_size - 8 } else { 0 };
        let mut ext_left = 0usize;
        let mut ext_right = 0usize;

        if extra > 0 {
            let max_left = extra.min(sp).min(qp);
            while ext_left < max_left {
                if query[qp - ext_left - 1] != packed_base_at(subject_packed, sp - ext_left - 1) {
                    break;
                }
                ext_left += 1;
            }
            let need_right = extra - ext_left;
            let max_right = if sp + 8 + need_right <= subject_len && qp + 8 + need_right <= query.len() { need_right } else { 0 };
            while ext_right < max_right {
                if query[qp + 8 + ext_right] != packed_base_at(subject_packed, sp + 8 + ext_right) {
                    break;
                }
                ext_right += 1;
            }
        }

        if ext_left + ext_right >= extra {
            let wq = qp - ext_left;
            let ws = sp - ext_left;
            let diag = (ws + query.len() - wq) & diag_mask;
            if !(last_hit[diag] >= 0 && (last_hit[diag] as usize) > ws) {
                if let Some(hsp) = extend_seed_packed(
                    query, subject_packed, subject_len, wq, ws,
                    reward, penalty, x_dropoff,
                    kbp, search_space, evalue_threshold, context,
                ) {
                    last_hit[diag] = hsp.subject_end;
                    hsps.push(hsp);
                } else {
                    last_hit[diag] = (ws + word_size) as i32;
                }
            }
        }
        q_pos = next[qp];
    }
}

/// Fast ungapped search on packed NCBI2na subject data.
/// Port of C engine's s_BlastSmallNaScanSubject_8_4 algorithm:
/// processes one packed byte (4 bases) per iteration using shift+OR hash update.
pub fn blastn_ungapped_search_packed(
    query_plus: &[u8],
    query_minus: &[u8],
    subject_packed: &[u8],
    subject_len: usize,
    word_size: usize,
    reward: i32,
    penalty: i32,
    x_dropoff: i32,
    kbp: &KarlinBlk,
    search_space: f64,
    evalue_threshold: f64,
) -> Vec<SearchHsp> {
    let mut hsps = Vec::new();

    let lut_word = word_size.min(8);
    let lut_size = 1usize << (2 * lut_word);
    let lut_mask = (lut_size - 1) as u32;

    for (context, query) in [(0i32, query_plus), (1i32, query_minus)] {
        if query.len() < word_size || subject_len < word_size {
            continue;
        }

        // Build lookup table + PV
        let mut lut: Vec<i32> = vec![-1; lut_size];
        let mut next: Vec<i32> = vec![-1; query.len()];
        let pv_size = (lut_size + 63) / 64;
        let mut pv: Vec<u64> = vec![0; pv_size]; // 8KB, fits in L1

        for i in (0..=(query.len() - word_size)).rev() {
            let key = word_hash_n(&query[i..i + lut_word], lut_word) as usize;
            next[i] = lut[key];
            lut[key] = i as i32;
            pv[key >> 6] |= 1u64 << (key & 63);
        }

        if subject_len < word_size { continue; }
        let end = subject_len - word_size + 1;

        // Diagonal tracking: prevent re-extending the same query-subject diagonal
        // Use power-of-2 sized table with hash (like C engine's diag_mask)
        let diag_array_len = (query.len() * 2).next_power_of_two().max(256);
        let diag_mask = diag_array_len - 1;
        let mut last_hit: Vec<i32> = vec![-1; diag_array_len];

        if lut_word == 8 {
            // Fast path: 8-mer = 16 bits. Use per-base scanning with packed_base_at
            // but with the efficient hash update (shift by 2 bits per base).
            // This is simpler and correct — checks every position.
            let mut hash: u32 = 0;
            for i in 0..7 {
                hash = (hash << 2) | packed_base_at(subject_packed, i) as u32;
            }
            let mut s_pos = 0;
            while s_pos < end {
                hash = ((hash << 2) | packed_base_at(subject_packed, s_pos + 7) as u32) & lut_mask;
                let h = hash as usize;
                if pv[h >> 6] & (1u64 << (h & 63)) != 0 {
                    process_hit(query, subject_packed, subject_len, word_size,
                        &lut, &next, &mut last_hit, diag_mask,
                        h, s_pos, reward, penalty, x_dropoff,
                        kbp, search_space, evalue_threshold, context, &mut hsps);
                }
                s_pos += 1;
            }
        } else {
            // Generic path: per-base scanning using packed_base_at
            let mut hash: u32 = 0;
            for i in 0..lut_word - 1 {
                hash = (hash << 2) | packed_base_at(subject_packed, i) as u32;
            }
            let mut s_pos = 0;
            while s_pos < end {
                hash = ((hash << 2) | packed_base_at(subject_packed, s_pos + lut_word - 1) as u32) & lut_mask;
                let h = hash as usize;
                if pv[h >> 6] & (1u64 << (h & 63)) != 0 {
                    let mut q_pos = lut[h];
                    while q_pos >= 0 {
                        let qp = q_pos as usize;
                        let diag = (s_pos + query.len() - qp) & diag_mask;
                        if last_hit[diag] >= 0 && (last_hit[diag] as usize) > s_pos {
                            q_pos = next[qp];
                            continue;
                        }
                        let mut ok = true;
                        if word_size > lut_word {
                            for k in lut_word..word_size {
                                if query[qp + k] != packed_base_at(subject_packed, s_pos + k) {
                                    ok = false; break;
                                }
                            }
                        }
                        if ok {
                            if let Some(hsp) = extend_seed_packed(
                                query, subject_packed, subject_len, qp, s_pos,
                                reward, penalty, x_dropoff,
                                kbp, search_space, evalue_threshold, context,
                            ) {
                                last_hit[diag] = hsp.subject_end;
                                hsps.push(hsp);
                            } else {
                                last_hit[diag] = (s_pos + word_size) as i32;
                            }
                        }
                        q_pos = next[qp];
                    }
                }
                s_pos += 1;
            }
        }
    }

    hsps.sort_by(|a, b| b.score.cmp(&a.score));
    dedup_hsps(&mut hsps);
    hsps
}

/// Extract a single base from packed NCBI2na data (4 bases per byte).
#[inline(always)]
fn packed_base_at(packed: &[u8], pos: usize) -> u8 {
    let byte = packed[pos >> 2];
    (byte >> (6 - 2 * (pos & 3))) & 3
}

/// Ungapped extension on packed subject data.
#[inline]
fn extend_seed_packed(
    query: &[u8],
    subject_packed: &[u8],
    subject_len: usize,
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
    while qi < query.len() && si < subject_len {
        let sb = packed_base_at(subject_packed, si);
        score += if query[qi] == sb { reward } else { penalty };
        if score > best_score {
            best_score = score;
            best_right = qi - q_seed + 1;
        }
        if best_score - score > x_dropoff { break; }
        qi += 1;
        si += 1;
    }

    // Extend left
    let mut left_score = 0i32;
    let mut best_left_score = 0i32;
    let mut best_left = 0usize;
    if q_seed > 0 && s_seed > 0 {
        let mut qi = q_seed - 1;
        let mut si = s_seed - 1;
        loop {
            let sb = packed_base_at(subject_packed, si);
            left_score += if query[qi] == sb { reward } else { penalty };
            if left_score > best_left_score {
                best_left_score = left_score;
                best_left = q_seed - qi;
            }
            if best_left_score - left_score > x_dropoff { break; }
            if qi == 0 || si == 0 { break; }
            qi -= 1;
            si -= 1;
        }
    }

    let total_score = best_score + best_left_score;
    if total_score <= 0 { return None; }

    let evalue = kbp.raw_to_evalue(total_score, search_space);
    if evalue > evalue_threshold { return None; }

    let q_start = (q_seed - best_left) as i32;
    let q_end = (q_seed + best_right) as i32;
    let s_start = (s_seed - best_left) as i32;
    let s_end = (s_seed + best_right) as i32;
    let align_len = q_end - q_start;

    // Count identities and build aligned sequences
    let mut num_ident = 0i32;
    let mut qseq_str = String::with_capacity(align_len as usize);
    let mut sseq_str = String::with_capacity(align_len as usize);
    for i in 0..align_len as usize {
        let qb = query[q_start as usize + i];
        let sb = packed_base_at(subject_packed, s_start as usize + i);
        if qb == sb { num_ident += 1; }
        qseq_str.push(blastna_to_iupac(qb));
        sseq_str.push(blastna_to_iupac(sb));
    }

    Some(SearchHsp {
        query_start: q_start,
        query_end: q_end,
        subject_start: s_start,
        subject_end: s_end,
        score: total_score,
        bit_score: kbp.raw_to_bit(total_score),
        evalue,
        num_ident,
        align_length: align_len,
        mismatches: align_len - num_ident,
        gap_opens: 0,
        context,
        qseq: Some(qseq_str),
        sseq: Some(sseq_str),
    })
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

    // Count identities and build aligned sequences
    let mut num_ident = 0;
    let mut qseq_str = String::with_capacity(align_len as usize);
    let mut sseq_str = String::with_capacity(align_len as usize);
    for i in 0..align_len as usize {
        let qb = query[q_start as usize + i];
        let sb = subject[s_start as usize + i];
        if qb == sb { num_ident += 1; }
        qseq_str.push(blastna_to_iupac(qb));
        sseq_str.push(blastna_to_iupac(sb));
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
        qseq: Some(qseq_str),
        sseq: Some(sseq_str),
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
    // Delegate to extended version with nomask = same as query (no separate masking)
    blastn_gapped_search_nomask(
        query_plus, query_minus, query_plus, query_minus,
        subject, word_size, reward, penalty, _gap_open, _gap_extend,
        x_dropoff, kbp, search_space, evalue_threshold,
    )
}

/// Gapped search with separate masked (for seeds) and unmasked (for gapped alignment) queries.
pub fn blastn_gapped_search_nomask(
    query_plus: &[u8],        // masked, for seed finding
    query_minus: &[u8],       // masked, for seed finding
    query_plus_nomask: &[u8], // unmasked, for gapped alignment
    query_minus_nomask: &[u8],// unmasked, for gapped alignment
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
    // First do ungapped search to find seeds (uses masked query)
    let ungapped = blastn_ungapped_search(
        query_plus, query_minus, subject,
        word_size, reward, penalty, x_dropoff,
        kbp, search_space, evalue_threshold * 100.0, // permissive threshold for seeds
    );

    let mut hsps = Vec::new();

    for seed in &ungapped {
        // Use UNMASKED query for gapped alignment (matches C engine's sequence_nomask)
        let query = if seed.context == 0 { query_plus_nomask } else { query_minus_nomask };

        // BLAST-style gapped alignment: bidirectional X-dropoff extension from seed
        let seed_q = ((seed.query_start + seed.query_end) / 2) as usize;
        let seed_s = ((seed.subject_start + seed.subject_end) / 2) as usize;

        if let Some(tb) = blast_gapped_align(
            query, subject, seed_q, seed_s,
            reward, penalty, _gap_open, _gap_extend, x_dropoff,
        ) {
            let evalue = kbp.raw_to_evalue(tb.score, search_space);
            if evalue > evalue_threshold { continue; }

            let q_slice = &query[tb.query_start..tb.query_end];
            let s_slice = &subject[tb.subject_start..tb.subject_end];
            let (align_len, num_ident, gap_opens) = tb.edit_script.count_identities(
                q_slice, s_slice,
            );
            let (qseq, sseq) = tb.edit_script.render_alignment(
                q_slice, s_slice, blastna_to_iupac,
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
                qseq: Some(qseq),
                sseq: Some(sseq),
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

/// Fast gapped search on packed NCBI2na subject.
/// Decodes subject once, then does seed finding + gapped alignment.
pub fn blastn_gapped_search_packed(
    query_plus: &[u8],
    query_minus: &[u8],
    query_plus_nomask: &[u8],
    query_minus_nomask: &[u8],
    subject_packed: &[u8],
    subject_len: usize,
    word_size: usize,
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    x_dropoff: i32,
    kbp: &KarlinBlk,
    search_space: f64,
    evalue_threshold: f64,
) -> Vec<SearchHsp> {
    // Two-phase search:
    // 1. Fast packed scan to find seeds (no decode needed for scanning)
    // 2. Decode subject only if seeds found (for verify + gapped alignment)
    let ungapped = blastn_ungapped_search_packed(
        query_plus, query_minus, subject_packed, subject_len,
        word_size, reward, penalty, x_dropoff,
        kbp, search_space, evalue_threshold * 100.0,
    );
    let mut hsps = Vec::new();
    if ungapped.is_empty() { return hsps; }

    let cutoff_score = {
        let e = evalue_threshold.max(1.0e-297);
        ((kbp.k * search_space / e).ln() / kbp.lambda).ceil() as i32
    };

    // Filter seeds by ungapped cutoff first (cheap, no decode needed)
    let passing_seeds: Vec<&SearchHsp> = ungapped.iter()
        .filter(|s| s.score >= cutoff_score)
        .collect();
    if passing_seeds.is_empty() { return hsps; }

    // Decode subject only once, only if needed
    let subject_decoded = decode_packed_ncbi2na(subject_packed, subject_len);

    for seed in &passing_seeds {
        let query = if seed.context == 0 { query_plus_nomask } else { query_minus_nomask };
        let seed_q = ((seed.query_start + seed.query_end) / 2) as usize;
        let seed_s = ((seed.subject_start + seed.subject_end) / 2) as usize;

        // Score-only preliminary gapped extension
        let prelim_score = blast_gapped_score_only(
            query, &subject_decoded, seed_q, seed_s,
            reward, penalty, gap_open, gap_extend, x_dropoff,
        );
        if prelim_score < cutoff_score { continue; }

        // Full traceback only for seeds that pass the score cutoff
        if let Some(tb) = blast_gapped_align(
            query, &subject_decoded, seed_q, seed_s,
            reward, penalty, gap_open, gap_extend, x_dropoff,
        ) {
            let evalue = kbp.raw_to_evalue(tb.score, search_space);
            if evalue > evalue_threshold { continue; }

            let q_slice = &query[tb.query_start..tb.query_end];
            let s_slice = &subject_decoded[tb.subject_start..tb.subject_end];
            let (align_len, num_ident, gap_opens) = tb.edit_script.count_identities(
                q_slice, s_slice,
            );
            let (qseq, sseq) = tb.edit_script.render_alignment(
                q_slice, s_slice, blastna_to_iupac,
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
                qseq: Some(qseq),
                sseq: Some(sseq),
            });
        } else {
            hsps.push((*seed).clone());
        }
    }

    hsps.sort_by(|a, b| b.score.cmp(&a.score));
    dedup_hsps(&mut hsps);
    hsps
}

/// Decode packed NCBI2na to per-base (0=A, 1=C, 2=G, 3=T).
pub fn decode_packed_ncbi2na(packed: &[u8], len: usize) -> Vec<u8> {
    let mut decoded = Vec::with_capacity(len);
    let full_bytes = len / 4;
    let remainder = len % 4;
    for i in 0..full_bytes {
        let b = packed[i];
        decoded.push((b >> 6) & 3);
        decoded.push((b >> 4) & 3);
        decoded.push((b >> 2) & 3);
        decoded.push(b & 3);
    }
    if remainder > 0 && full_bytes < packed.len() {
        let b = packed[full_bytes];
        for j in 0..remainder {
            decoded.push((b >> (6 - 2 * j)) & 3);
        }
    }
    decoded
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
