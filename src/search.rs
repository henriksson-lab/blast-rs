//! Pure Rust BLAST search engine — no FFI.
//! This module implements the complete blastn search pipeline in Rust.

use crate::encoding::NCBI4NA_TO_BLASTNA;
use crate::itree::{Interval, IntervalTree};
use crate::sequence::blastna_to_iupac;
use crate::stat::KarlinBlk;
use crate::traceback::{
    blast_gapped_align, blast_gapped_score_only, reevaluate_with_ambiguities_gapped,
};

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

/// Port of NCBI `ScoreCompareHSPs` (`blast_hits.c:1330`) for `SearchHsp`.
/// Same tie-breaker structure as `hspstream::score_compare_hsps` but over
/// `SearchHsp`'s `*_start`/`*_end` fields instead of `*_offset`/`*_end`.
fn score_compare_search_hsps(a: &SearchHsp, b: &SearchHsp) -> std::cmp::Ordering {
    b.score
        .cmp(&a.score)
        .then_with(|| a.subject_start.cmp(&b.subject_start))
        .then_with(|| b.subject_end.cmp(&a.subject_end))
        .then_with(|| a.query_start.cmp(&b.query_start))
        .then_with(|| b.query_end.cmp(&a.query_end))
}

/// Result of searching one query against one subject.
#[derive(Debug)]
pub struct SubjectResult {
    pub oid: i32,
    pub hsps: Vec<SearchHsp>,
}

struct NaLookup<'a> {
    context: i32,
    query: &'a [u8],
    lut_word: usize,
    lut_mask: u32,
    scan_start: usize,
    scan_step: usize,
    lut: Vec<i32>,
    next: Vec<i32>,
    pv: Vec<u64>,
    diag_array_len: usize,
    diag_mask: usize,
}

/// Query-side nucleotide lookup tables prepared once per BLASTN query.
///
/// NCBI BLAST builds this lookup table before scanning database subjects.
/// Rebuilding it for every OID is prohibitively expensive on large databases.
pub struct PreparedBlastnQuery<'a> {
    word_size: usize,
    lookups: Vec<NaLookup<'a>>,
    paired_pv: Option<Vec<u64>>,
}

impl<'a> PreparedBlastnQuery<'a> {
    pub fn new(query_plus: &'a [u8], query_minus: &'a [u8], word_size: usize) -> Self {
        Self::new_with_selector(query_plus, query_minus, word_size, choose_small_na_lut_word)
    }

    pub fn new_megablast(query_plus: &'a [u8], query_minus: &'a [u8], word_size: usize) -> Self {
        Self::new_with_selector(
            query_plus,
            query_minus,
            word_size,
            choose_contiguous_mb_lut_word,
        )
    }

    fn new_with_selector(
        query_plus: &'a [u8],
        query_minus: &'a [u8],
        word_size: usize,
        choose_lut_word: fn(usize, usize) -> usize,
    ) -> Self {
        let mut lookups = Vec::with_capacity(2);
        let approx_entries = query_plus.len().saturating_sub(word_size).saturating_add(1)
            + query_minus
                .len()
                .saturating_sub(word_size)
                .saturating_add(1);
        if let Some(lookup) =
            NaLookup::new(0, query_plus, word_size, approx_entries, choose_lut_word)
        {
            lookups.push(lookup);
        }
        if let Some(lookup) =
            NaLookup::new(1, query_minus, word_size, approx_entries, choose_lut_word)
        {
            lookups.push(lookup);
        }
        let paired_pv = if lookups.len() == 2
            && lookups[0].lut_word == lookups[1].lut_word
            && lookups[0].lut_mask == lookups[1].lut_mask
            && lookups[0].pv.len() == lookups[1].pv.len()
        {
            Some(
                lookups[0]
                    .pv
                    .iter()
                    .zip(&lookups[1].pv)
                    .map(|(a, b)| a | b)
                    .collect(),
            )
        } else {
            None
        };
        PreparedBlastnQuery {
            word_size,
            lookups,
            paired_pv,
        }
    }

    pub fn is_empty(&self) -> bool {
        self.lookups.is_empty()
    }

    pub fn last_hit_scratch(&self) -> Vec<Vec<i32>> {
        self.lookups
            .iter()
            .map(|lookup| vec![-1; lookup.diag_array_len])
            .collect()
    }
}

impl<'a> NaLookup<'a> {
    fn new(
        context: i32,
        query: &'a [u8],
        word_size: usize,
        approx_entries: usize,
        choose_lut_word: fn(usize, usize) -> usize,
    ) -> Option<Self> {
        if query.len() < word_size {
            return None;
        }

        let lut_word = choose_lut_word(word_size, approx_entries);
        let lut_size = 1usize << (2 * lut_word);
        let lut_mask = (lut_size - 1) as u32;
        let scan_start = word_size - lut_word;
        let scan_step = scan_start + 1;

        let mut lut: Vec<i32> = vec![-1; lut_size];
        let mut next: Vec<i32> = vec![-1; query.len()];
        let pv_size = lut_size.div_ceil(64);
        let mut pv: Vec<u64> = vec![0; pv_size];

        for i in (0..=(query.len() - lut_word)).rev() {
            if i + word_size <= query.len() && has_ambiguous_base(&query[i..i + word_size]) {
                continue;
            }
            let key = word_hash_n(&query[i..i + lut_word], lut_word) as usize;
            next[i] = lut[key];
            lut[key] = i as i32;
            pv[key >> 6] |= 1u64 << (key & 63);
        }

        let diag_array_len = (query.len() * 2).next_power_of_two().max(256);
        let diag_mask = diag_array_len - 1;

        Some(NaLookup {
            context,
            query,
            lut_word,
            lut_mask,
            scan_start,
            scan_step,
            lut,
            next,
            pv,
            diag_array_len,
            diag_mask,
        })
    }
}

fn choose_small_na_lut_word(word_size: usize, approx_entries: usize) -> usize {
    match word_size {
        4..=6 => word_size,
        7 => {
            if approx_entries < 250 {
                6
            } else {
                7
            }
        }
        8 => {
            if approx_entries < 8500 {
                7
            } else {
                8
            }
        }
        9 => {
            if approx_entries < 1250 {
                7
            } else {
                8
            }
        }
        10 => {
            if approx_entries < 1250 {
                7
            } else {
                8
            }
        }
        11 | 12 => 8,
        _ => 8,
    }
}

fn choose_contiguous_mb_lut_word(word_size: usize, approx_entries: usize) -> usize {
    match word_size {
        0..=8 => word_size,
        9 | 10 => {
            if approx_entries < 18_000 {
                9
            } else {
                10
            }
        }
        11 => {
            if approx_entries < 180_000 {
                10
            } else {
                11
            }
        }
        12 => {
            if approx_entries < 18_000 {
                9
            } else if approx_entries < 60_000 {
                10
            } else if approx_entries < 900_000 {
                11
            } else {
                12
            }
        }
        _ => {
            if approx_entries < 300_000 {
                11
            } else {
                12
            }
        }
    }
}

/// Perform a simple ungapped nucleotide word search.
/// Finds exact word matches between query and subject, extends them,
/// and returns significant HSPs.
pub fn blastn_ungapped_search(
    query_plus: &[u8],  // BLASTNA encoded, plus strand
    query_minus: &[u8], // BLASTNA encoded, minus strand (RC)
    subject: &[u8],     // BLASTNA decoded subject
    word_size: usize,
    reward: i32,
    penalty: i32,
    x_dropoff: i32,
    kbp: &KarlinBlk,
    search_space: f64,
    evalue_threshold: f64,
) -> Vec<SearchHsp> {
    blastn_ungapped_search_inner(
        query_plus,
        query_minus,
        query_plus,
        query_minus,
        subject,
        word_size,
        reward,
        penalty,
        x_dropoff,
        kbp,
        search_space,
        evalue_threshold,
        true,
    )
}

/// Ungapped nucleotide word search without final HSP containment filtering.
/// This is used for `blastn -ungapped`, where NCBI reports overlapping
/// ungapped HSPs that are later controlled by hit-list limits.
pub fn blastn_ungapped_search_no_dedup(
    query_plus: &[u8],
    query_minus: &[u8],
    subject: &[u8],
    word_size: usize,
    reward: i32,
    penalty: i32,
    x_dropoff: i32,
    kbp: &KarlinBlk,
    search_space: f64,
    evalue_threshold: f64,
) -> Vec<SearchHsp> {
    blastn_ungapped_search_no_dedup_nomask(
        query_plus,
        query_minus,
        query_plus,
        query_minus,
        subject,
        word_size,
        reward,
        penalty,
        x_dropoff,
        kbp,
        search_space,
        evalue_threshold,
    )
}

pub fn blastn_ungapped_search_no_dedup_nomask(
    query_plus: &[u8],
    query_minus: &[u8],
    query_plus_nomask: &[u8],
    query_minus_nomask: &[u8],
    subject: &[u8],
    word_size: usize,
    reward: i32,
    penalty: i32,
    x_dropoff: i32,
    kbp: &KarlinBlk,
    search_space: f64,
    evalue_threshold: f64,
) -> Vec<SearchHsp> {
    blastn_ungapped_search_inner(
        query_plus,
        query_minus,
        query_plus_nomask,
        query_minus_nomask,
        subject,
        word_size,
        reward,
        penalty,
        x_dropoff,
        kbp,
        search_space,
        evalue_threshold,
        false,
    )
}

#[allow(clippy::too_many_arguments)]
fn blastn_ungapped_search_inner(
    query_plus: &[u8],
    query_minus: &[u8],
    query_plus_nomask: &[u8],
    query_minus_nomask: &[u8],
    subject: &[u8],
    word_size: usize,
    reward: i32,
    penalty: i32,
    x_dropoff: i32,
    kbp: &KarlinBlk,
    search_space: f64,
    evalue_threshold: f64,
    dedup: bool,
) -> Vec<SearchHsp> {
    let mut hsps = Vec::new();

    let lut_word = word_size.min(8);
    let lut_size = 1usize << (2 * lut_word);

    for (context, query, query_nomask) in [
        (0i32, query_plus, query_plus_nomask),
        (1i32, query_minus, query_minus_nomask),
    ] {
        if query.len() < word_size || subject.len() < word_size {
            continue;
        }
        let mut last_hit = vec![-1i32; subject.len() + query.len() + 1];

        // Build lookup table from query
        let mut lut: Vec<i32> = vec![-1; lut_size];
        let mut next: Vec<i32> = vec![-1; query.len()];

        for i in (0..=(query.len() - word_size)).rev() {
            if has_ambiguous_base(&query[i..i + word_size]) {
                continue;
            }
            let key = word_hash_n(&query[i..i + lut_word], lut_word) as usize;
            next[i] = lut[key];
            lut[key] = i as i32;
        }

        // Scan subject
        let mut s_pos = 0;
        while s_pos + word_size <= subject.len() {
            if has_ambiguous_base(&subject[s_pos..s_pos + word_size]) {
                s_pos += 1;
                continue;
            }
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
                    let diag = s_pos + query.len() - qp;
                    if last_hit[diag] >= 0 && (last_hit[diag] as usize) > s_pos {
                        q_pos = next[qp];
                        continue;
                    }
                    if let Some(hsp) = extend_seed(
                        query_nomask,
                        subject,
                        qp,
                        s_pos,
                        reward,
                        penalty,
                        x_dropoff,
                        kbp,
                        search_space,
                        evalue_threshold,
                        context,
                    ) {
                        last_hit[diag] = hsp.subject_end;
                        hsps.push(hsp);
                    } else {
                        last_hit[diag] = (s_pos + word_size) as i32;
                    }
                }
                q_pos = next[qp];
            }
            s_pos += 1;
        }
    }

    if dedup {
        hsps.sort_by(score_compare_search_hsps);
        dedup_hsps_with_min_diag_separation(
            &mut hsps,
            min_diag_separation_for_ungapped(word_size, reward, penalty),
        );
    } else {
        sort_ungapped_hsps_ncbi_no_dedup(&mut hsps);
    }
    hsps
}

/// Process a seed hit from the scan loop (verify word, diagonal check, extend).
#[inline]
fn process_hit(
    query: &[u8],
    subject_packed: &[u8],
    subject_len: usize,
    word_size: usize,
    lut_word: usize,
    lut: &[i32],
    next: &[i32],
    last_hit: &mut [i32],
    diag_mask: usize,
    h: usize,
    base_pos: usize,
    reward: i32,
    penalty: i32,
    x_dropoff: i32,
    kbp: &KarlinBlk,
    search_space: f64,
    evalue_threshold: f64,
    context: i32,
    hsps: &mut Vec<SearchHsp>,
) {
    let sp = base_pos;
    let mut q_pos = lut[h];
    while q_pos >= 0 {
        let qp = q_pos as usize;
        process_offset_pair(
            query,
            subject_packed,
            subject_len,
            word_size,
            lut_word,
            qp,
            sp,
            last_hit,
            diag_mask,
            reward,
            penalty,
            x_dropoff,
            kbp,
            search_space,
            evalue_threshold,
            context,
            hsps,
        );
        q_pos = next[qp];
    }
}

#[inline(always)]
fn pv_has_hash(pv: &[u64], h: usize) -> bool {
    let bucket = h >> 6;
    debug_assert!(bucket < pv.len());
    let value = unsafe { *pv.get_unchecked(bucket) };
    value & (1u64 << (h & 63)) != 0
}

#[inline(always)]
fn pv_bucket_unchecked(pv: &[u64], bucket: usize) -> u64 {
    debug_assert!(bucket < pv.len());
    unsafe { *pv.get_unchecked(bucket) }
}

#[allow(clippy::too_many_arguments)]
fn process_offset_pair(
    query: &[u8],
    subject_packed: &[u8],
    subject_len: usize,
    word_size: usize,
    lut_word: usize,
    qp: usize,
    sp: usize,
    last_hit: &mut [i32],
    diag_mask: usize,
    reward: i32,
    penalty: i32,
    x_dropoff: i32,
    kbp: &KarlinBlk,
    search_space: f64,
    evalue_threshold: f64,
    context: i32,
    hsps: &mut Vec<SearchHsp>,
) {
    let extra = word_size.saturating_sub(lut_word);
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
        let max_right = if sp + lut_word + need_right <= subject_len
            && qp + lut_word + need_right <= query.len()
        {
            need_right
        } else {
            0
        };
        while ext_right < max_right {
            if query[qp + lut_word + ext_right]
                != packed_base_at(subject_packed, sp + lut_word + ext_right)
            {
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
                query,
                subject_packed,
                subject_len,
                wq,
                ws,
                reward,
                penalty,
                x_dropoff,
                kbp,
                search_space,
                evalue_threshold,
                context,
            ) {
                last_hit[diag] = hsp.subject_end;
                hsps.push(hsp);
            } else {
                last_hit[diag] = (ws + word_size) as i32;
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn process_offset_pair_decoded(
    query: &[u8],
    subject: &[u8],
    word_size: usize,
    lut_word: usize,
    qp: usize,
    sp: usize,
    last_hit: &mut [i32],
    diag_mask: usize,
    reward: i32,
    penalty: i32,
    x_dropoff: i32,
    kbp: &KarlinBlk,
    search_space: f64,
    evalue_threshold: f64,
    context: i32,
    hsps: &mut Vec<SearchHsp>,
) {
    let extra = word_size.saturating_sub(lut_word);
    let mut ext_left = 0usize;
    let mut ext_right = 0usize;

    if extra > 0 {
        let max_left = extra.min(sp).min(qp);
        while ext_left < max_left {
            if query[qp - ext_left - 1] != subject[sp - ext_left - 1] {
                break;
            }
            ext_left += 1;
        }
        let need_right = extra - ext_left;
        let max_right = if sp + lut_word + need_right <= subject.len()
            && qp + lut_word + need_right <= query.len()
        {
            need_right
        } else {
            0
        };
        while ext_right < max_right {
            if query[qp + lut_word + ext_right] != subject[sp + lut_word + ext_right] {
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
            if let Some(hsp) = extend_seed(
                query,
                subject,
                wq,
                ws,
                reward,
                penalty,
                x_dropoff,
                kbp,
                search_space,
                evalue_threshold,
                context,
            ) {
                last_hit[diag] = hsp.subject_end;
                hsps.push(hsp);
            } else {
                last_hit[diag] = (ws + word_size) as i32;
            }
        }
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
    let prepared = PreparedBlastnQuery::new(query_plus, query_minus, word_size);
    blastn_ungapped_search_packed_prepared(
        &prepared,
        subject_packed,
        subject_len,
        reward,
        penalty,
        x_dropoff,
        kbp,
        search_space,
        evalue_threshold,
    )
}

pub fn blastn_ungapped_search_packed_prepared(
    prepared: &PreparedBlastnQuery<'_>,
    subject_packed: &[u8],
    subject_len: usize,
    reward: i32,
    penalty: i32,
    x_dropoff: i32,
    kbp: &KarlinBlk,
    search_space: f64,
    evalue_threshold: f64,
) -> Vec<SearchHsp> {
    let mut last_hit_scratch = prepared.last_hit_scratch();
    blastn_ungapped_search_packed_prepared_with_scratch(
        prepared,
        subject_packed,
        subject_len,
        reward,
        penalty,
        x_dropoff,
        kbp,
        search_space,
        evalue_threshold,
        &mut last_hit_scratch,
    )
}

pub fn blastn_ungapped_search_packed_prepared_with_scratch(
    prepared: &PreparedBlastnQuery<'_>,
    subject_packed: &[u8],
    subject_len: usize,
    reward: i32,
    penalty: i32,
    x_dropoff: i32,
    kbp: &KarlinBlk,
    search_space: f64,
    evalue_threshold: f64,
    last_hit_scratch: &mut [Vec<i32>],
) -> Vec<SearchHsp> {
    blastn_ungapped_search_packed_prepared_with_scratch_inner(
        prepared,
        subject_packed,
        subject_len,
        reward,
        penalty,
        x_dropoff,
        kbp,
        search_space,
        evalue_threshold,
        last_hit_scratch,
        true,
    )
}

pub fn blastn_ungapped_search_packed_prepared_with_scratch_no_dedup(
    prepared: &PreparedBlastnQuery<'_>,
    subject_packed: &[u8],
    subject_len: usize,
    reward: i32,
    penalty: i32,
    x_dropoff: i32,
    kbp: &KarlinBlk,
    search_space: f64,
    evalue_threshold: f64,
    last_hit_scratch: &mut [Vec<i32>],
) -> Vec<SearchHsp> {
    blastn_ungapped_search_packed_prepared_with_scratch_inner(
        prepared,
        subject_packed,
        subject_len,
        reward,
        penalty,
        x_dropoff,
        kbp,
        search_space,
        evalue_threshold,
        last_hit_scratch,
        false,
    )
}

#[allow(clippy::too_many_arguments)]
fn blastn_ungapped_search_packed_prepared_with_scratch_inner(
    prepared: &PreparedBlastnQuery<'_>,
    subject_packed: &[u8],
    subject_len: usize,
    reward: i32,
    penalty: i32,
    x_dropoff: i32,
    kbp: &KarlinBlk,
    search_space: f64,
    evalue_threshold: f64,
    last_hit_scratch: &mut [Vec<i32>],
    dedup: bool,
) -> Vec<SearchHsp> {
    let mut hsps = Vec::new();

    if prepared.is_empty() || subject_len < prepared.word_size {
        return hsps;
    }

    let end = subject_len - prepared.word_size + 1;

    assert_eq!(
        last_hit_scratch.len(),
        prepared.lookups.len(),
        "last-hit scratch must match prepared lookup contexts"
    );

    if prepared.lookups.len() == 2
        && prepared.lookups[0].lut_word == 6
        && prepared.lookups[1].lut_word == 6
        && prepared.lookups[0].scan_start == 1
        && prepared.lookups[1].scan_start == 1
        && prepared.lookups[0].scan_step == 2
        && prepared.lookups[1].scan_step == 2
        && prepared.paired_pv.is_some()
    {
        let (left_scratch, right_scratch) = last_hit_scratch.split_at_mut(1);
        let lookup0 = &prepared.lookups[0];
        let lookup1 = &prepared.lookups[1];
        let last_hit0 = &mut left_scratch[0];
        let last_hit1 = &mut right_scratch[0];

        if last_hit0.len() != lookup0.diag_array_len {
            last_hit0.resize(lookup0.diag_array_len, -1);
        } else {
            last_hit0.fill(-1);
        }
        if last_hit1.len() != lookup1.diag_array_len {
            last_hit1.resize(lookup1.diag_array_len, -1);
        } else {
            last_hit1.fill(-1);
        }

        scan_byte_oriented_lut6_step2_pair(
            subject_packed,
            subject_len,
            end,
            lookup0,
            last_hit0,
            lookup1,
            last_hit1,
            prepared
                .paired_pv
                .as_deref()
                .expect("paired PV required for paired lookup scan"),
            prepared.word_size,
            reward,
            penalty,
            x_dropoff,
            kbp,
            search_space,
            evalue_threshold,
            &mut hsps,
        );

        if dedup {
            hsps.sort_by(score_compare_search_hsps);
            dedup_hsps_with_min_diag_separation(
                &mut hsps,
                min_diag_separation_for_ungapped(prepared.word_size, reward, penalty),
            );
        } else {
            sort_ungapped_hsps_ncbi_no_dedup(&mut hsps);
        }
        return hsps;
    }

    for (lookup, last_hit) in prepared.lookups.iter().zip(last_hit_scratch.iter_mut()) {
        if last_hit.len() != lookup.diag_array_len {
            last_hit.resize(lookup.diag_array_len, -1);
        } else {
            last_hit.fill(-1);
        }

        if lookup.lut_word == 8 && prepared.word_size >= 8 {
            // Fast path: step-4 scanning, port of C engine's
            // s_BlastSmallNaScanSubject_8_4 from blast_nascan.c.
            // Processes whole packed bytes (4 bases at a time), unrolled 8x.
            let scan_step = prepared.word_size - lookup.lut_word + 1;
            if scan_step == 4 {
                scan_step4_unrolled(
                    subject_packed,
                    subject_len,
                    prepared.word_size,
                    lookup.query,
                    &lookup.lut,
                    &lookup.next,
                    &lookup.pv,
                    last_hit,
                    lookup.diag_mask,
                    reward,
                    penalty,
                    x_dropoff,
                    kbp,
                    search_space,
                    evalue_threshold,
                    lookup.context,
                    &mut hsps,
                );
            } else {
                scan_step1(
                    subject_packed,
                    subject_len,
                    end,
                    lookup.lut_word,
                    lookup.lut_mask,
                    lookup.scan_start,
                    lookup.scan_step,
                    lookup.query,
                    prepared.word_size,
                    &lookup.lut,
                    &lookup.next,
                    &lookup.pv,
                    last_hit,
                    lookup.diag_mask,
                    reward,
                    penalty,
                    x_dropoff,
                    kbp,
                    search_space,
                    evalue_threshold,
                    lookup.context,
                    &mut hsps,
                );
            }
        } else {
            scan_byte_oriented_step1(
                subject_packed,
                subject_len,
                end,
                lookup.lut_word,
                lookup.lut_mask,
                lookup.scan_start,
                lookup.scan_step,
                lookup.query,
                prepared.word_size,
                &lookup.lut,
                &lookup.next,
                &lookup.pv,
                last_hit,
                lookup.diag_mask,
                reward,
                penalty,
                x_dropoff,
                kbp,
                search_space,
                evalue_threshold,
                lookup.context,
                &mut hsps,
            );
        }
    }

    if dedup {
        hsps.sort_by(score_compare_search_hsps);
        dedup_hsps_with_min_diag_separation(
            &mut hsps,
            min_diag_separation_for_ungapped(prepared.word_size, reward, penalty),
        );
    } else {
        sort_ungapped_hsps_ncbi_no_dedup(&mut hsps);
    }
    hsps
}

#[allow(clippy::too_many_arguments)]
pub fn blastn_ungapped_search_decoded_prepared_with_scratch_no_dedup(
    prepared: &PreparedBlastnQuery<'_>,
    subject: &[u8],
    reward: i32,
    penalty: i32,
    x_dropoff: i32,
    kbp: &KarlinBlk,
    search_space: f64,
    evalue_threshold: f64,
    last_hit_scratch: &mut [Vec<i32>],
) -> Vec<SearchHsp> {
    let mut hsps = Vec::new();

    if prepared.is_empty() || subject.len() < prepared.word_size {
        return hsps;
    }

    assert_eq!(
        last_hit_scratch.len(),
        prepared.lookups.len(),
        "last-hit scratch must match prepared lookup contexts"
    );

    let end = subject.len() - prepared.word_size + 1;
    for (lookup, last_hit) in prepared.lookups.iter().zip(last_hit_scratch.iter_mut()) {
        if last_hit.len() != lookup.diag_array_len {
            last_hit.resize(lookup.diag_array_len, -1);
        } else {
            last_hit.fill(-1);
        }

        let mut s_pos = lookup.scan_start;
        while s_pos < end {
            if !has_ambiguous_base(&subject[s_pos..s_pos + prepared.word_size]) {
                let h =
                    word_hash_n(&subject[s_pos..s_pos + lookup.lut_word], lookup.lut_word) as usize;
                if pv_has_hash(&lookup.pv, h) {
                    let mut q_pos = lookup.lut[h];
                    while q_pos >= 0 {
                        process_offset_pair_decoded(
                            lookup.query,
                            subject,
                            prepared.word_size,
                            lookup.lut_word,
                            q_pos as usize,
                            s_pos,
                            last_hit,
                            lookup.diag_mask,
                            reward,
                            penalty,
                            x_dropoff,
                            kbp,
                            search_space,
                            evalue_threshold,
                            lookup.context,
                            &mut hsps,
                        );
                        q_pos = lookup.next[q_pos as usize];
                    }
                }
            }
            s_pos += lookup.scan_step;
        }
    }

    sort_ungapped_hsps_ncbi_no_dedup(&mut hsps);
    hsps
}

fn sort_ungapped_hsps_ncbi_no_dedup(hsps: &mut [SearchHsp]) {
    hsps.sort_by(|a, b| {
        b.score
            .cmp(&a.score)
            .then_with(|| a.subject_start.cmp(&b.subject_start))
            .then_with(|| a.query_start.cmp(&b.query_start))
            .then_with(|| a.subject_end.cmp(&b.subject_end))
            .then_with(|| a.query_end.cmp(&b.query_end))
    });
}

/// Step-4 unrolled scanner — port of C engine's s_BlastSmallNaScanSubject_8_4.
/// Processes whole packed bytes (4 bases at a time), 8 bytes (32 bases) per loop iteration.
/// Only checks every 4th subject position; the process_hit function handles
/// extending to verify the remaining bases of the full word.
#[allow(clippy::too_many_arguments)]
fn scan_step4_unrolled(
    subject_packed: &[u8],
    subject_len: usize,
    word_size: usize,
    query: &[u8],
    lut: &[i32],
    next: &[i32],
    pv: &[u64],
    last_hit: &mut [i32],
    diag_mask: usize,
    reward: i32,
    penalty: i32,
    x_dropoff: i32,
    kbp: &KarlinBlk,
    search_space: f64,
    evalue_threshold: f64,
    context: i32,
    hsps: &mut Vec<SearchHsp>,
) {
    if subject_len < word_size {
        return;
    }
    let end = subject_len - word_size + 1;
    // scan_range in subject positions (base coordinates)
    let scan_start: usize = 0;
    let scan_end = end.saturating_sub(1); // inclusive end

    // We need at least 2 packed bytes to form the initial 8-mer hash
    if subject_packed.len() < 2 {
        return;
    }

    let lut_mask: u32 = 0xFFFF; // 2 * 8 = 16 bits

    // Macro-like inline: check PV and process hits for a scan position
    // s_pos is the subject position (in base coordinates) of the 8-mer start
    macro_rules! check_and_process {
        ($hash:expr, $s_pos:expr) => {
            let h = ($hash & lut_mask) as usize;
            if pv_has_hash(pv, h) {
                process_hit(
                    query,
                    subject_packed,
                    subject_len,
                    word_size,
                    8,
                    lut,
                    next,
                    last_hit,
                    diag_mask,
                    h,
                    $s_pos,
                    reward,
                    penalty,
                    x_dropoff,
                    kbp,
                    search_space,
                    evalue_threshold,
                    context,
                    hsps,
                );
            }
        };
    }

    // The hash covers 8 bases = 2 packed bytes. With step-4 (= 1 byte), each
    // iteration shifts in one new byte and checks the resulting 16-bit index.
    //
    // Port of C's Duff's device unrolled loop from s_BlastSmallNaScanSubject_8_4.
    // Each packed byte encodes 4 bases. Advancing by 1 byte = 4 bases = step 4.
    //
    // s_pos = base position = byte_index * 4
    // hash  = (prev_hash << 8) | new_byte, masked to 16 bits

    let s_ptr = subject_packed.as_ptr();
    let s_len_bytes = subject_packed.len();

    // Number of step-4 positions to scan
    let num_steps = if scan_end >= scan_start {
        (scan_end - scan_start) / 4 + 1
    } else {
        0
    };
    if num_steps == 0 {
        return;
    }

    // Starting byte index (subject position / 4)
    let start_byte = scan_start / 4;

    // Initialize hash from first byte
    let mut init_index: u32 = if start_byte < s_len_bytes {
        unsafe { *s_ptr.add(start_byte) as u32 }
    } else {
        return;
    };

    // Process in groups of 8 bytes (32 bases) with Duff's device unrolling
    let full_groups = num_steps / 8;
    let remainder = num_steps % 8;

    // Handle remainder first (like C's Duff's device switch entry)
    let mut byte_idx = start_byte;
    let mut s_pos = scan_start;

    // Process remainder positions (0..remainder)
    for _ in 0..remainder {
        byte_idx += 1;
        if byte_idx >= s_len_bytes {
            break;
        }
        init_index = (init_index << 8) | unsafe { *s_ptr.add(byte_idx) as u32 };
        if s_pos < end {
            check_and_process!(init_index, s_pos);
        }
        s_pos += 4;
    }

    // Process full groups of 8
    for _ in 0..full_groups {
        if byte_idx + 8 >= s_len_bytes || s_pos + 28 >= subject_len {
            break;
        }

        // Unrolled 8 iterations: each shifts in one packed byte (4 bases)
        unsafe {
            init_index = (init_index << 8) | *s_ptr.add(byte_idx + 1) as u32;
            check_and_process!(init_index, s_pos);

            init_index = (init_index << 8) | *s_ptr.add(byte_idx + 2) as u32;
            check_and_process!(init_index, s_pos + 4);

            init_index = (init_index << 8) | *s_ptr.add(byte_idx + 3) as u32;
            check_and_process!(init_index, s_pos + 8);

            init_index = (init_index << 8) | *s_ptr.add(byte_idx + 4) as u32;
            check_and_process!(init_index, s_pos + 12);

            init_index = (init_index << 8) | *s_ptr.add(byte_idx + 5) as u32;
            check_and_process!(init_index, s_pos + 16);

            init_index = (init_index << 8) | *s_ptr.add(byte_idx + 6) as u32;
            check_and_process!(init_index, s_pos + 20);

            init_index = (init_index << 8) | *s_ptr.add(byte_idx + 7) as u32;
            check_and_process!(init_index, s_pos + 24);

            init_index = (init_index << 8) | *s_ptr.add(byte_idx + 8) as u32;
            check_and_process!(init_index, s_pos + 28);
        }
        byte_idx += 8;
        s_pos += 32;
    }
}

/// Byte-oriented step-1 scanner for lut_word < 8.
/// Port of C engine's s_BlastSmallNaScanSubject_7_1.
/// Reads packed bytes and extracts multiple hash values by shifting,
/// processing 4 subject positions per byte advance.
#[allow(clippy::too_many_arguments)]
fn scan_byte_oriented_step1(
    subject_packed: &[u8],
    subject_len: usize,
    end: usize,
    lut_word: usize,
    lut_mask: u32,
    scan_start: usize,
    scan_step: usize,
    query: &[u8],
    word_size: usize,
    lut: &[i32],
    next: &[i32],
    pv: &[u64],
    last_hit: &mut [i32],
    diag_mask: usize,
    reward: i32,
    penalty: i32,
    x_dropoff: i32,
    kbp: &KarlinBlk,
    search_space: f64,
    evalue_threshold: f64,
    context: i32,
    hsps: &mut Vec<SearchHsp>,
) {
    if end == 0 || subject_packed.is_empty() {
        return;
    }

    if lut_word == 6 && scan_start == 1 && scan_step == 2 {
        scan_byte_oriented_lut6_step2(
            subject_packed,
            subject_len,
            end,
            lut_mask,
            query,
            word_size,
            lut,
            next,
            pv,
            last_hit,
            diag_mask,
            reward,
            penalty,
            x_dropoff,
            kbp,
            search_space,
            evalue_threshold,
            context,
            hsps,
        );
        return;
    }

    if scan_step != 1 {
        let mut s_pos = scan_start;
        while s_pos < end {
            let h = packed_hash_at(subject_packed, s_pos, lut_word, lut_mask);
            if pv[h >> 6] & (1u64 << (h & 63)) != 0 {
                process_hit(
                    query,
                    subject_packed,
                    subject_len,
                    word_size,
                    lut_word,
                    lut,
                    next,
                    last_hit,
                    diag_mask,
                    h,
                    s_pos,
                    reward,
                    penalty,
                    x_dropoff,
                    kbp,
                    search_space,
                    evalue_threshold,
                    context,
                    hsps,
                );
            }
            s_pos += scan_step;
        }
        return;
    }

    let s_ptr = subject_packed.as_ptr();
    let s_len = subject_packed.len();

    // Number of bits in lut_word (= 2 * lut_word)
    let lut_bits = 2 * lut_word;

    // We process 4 positions per byte. For each group of 4 positions
    // (aligned to packed byte boundary), we read 2-3 bytes and extract
    // 4 different hash values by shifting.
    //
    // For lut_word=7 (14 bits): need ceil(14/8)=2 bytes to cover the hash.
    // Position 0 in a byte-group: hash = (byte0 << 8 | byte1) >> 2  (bits 16..2)
    // Position 1: hash = (byte0 << 8 | byte1) & mask                (bits 14..0)
    // Position 2: hash = (byte0 << 16 | byte1 << 8 | byte2) >> 6    (bits 22..8, masked)
    // Position 3: hash = (above >> 4) & mask                        (bits 18..4, shifted)
    //
    // General: bytes_needed = ceil((lut_bits + 6) / 8) to cover all 4 shifts

    let mut s_pos: usize = 0;
    let mut byte_idx: usize = 0;

    // Process 4 positions per iteration (one packed byte advance)
    while s_pos + 3 < end && byte_idx + 2 < s_len {
        // Read 3 bytes (enough for any lut_word <= 8)
        let b0: u32;
        let b1: u32;
        let b2: u32;
        unsafe {
            b0 = *s_ptr.add(byte_idx) as u32;
            b1 = *s_ptr.add(byte_idx + 1) as u32;
            b2 = if byte_idx + 2 < s_len {
                *s_ptr.add(byte_idx + 2) as u32
            } else {
                0
            };
        }
        let combined = (b0 << 16) | (b1 << 8) | b2;

        // Extract 4 hash values for positions s_pos, s_pos+1, s_pos+2, s_pos+3
        // Position within byte: hash shifts by 2 bits per position
        // Base shift for position 0 = 24 - lut_bits (since combined is 24 bits)
        let base_shift = 24 - lut_bits;

        macro_rules! check_pos {
            ($shift:expr, $off:expr) => {
                let pos = s_pos + $off;
                if pos < end {
                    let h = ((combined >> $shift) & lut_mask) as usize;
                    if pv_has_hash(pv, h) {
                        process_hit(
                            query,
                            subject_packed,
                            subject_len,
                            word_size,
                            lut_word,
                            lut,
                            next,
                            last_hit,
                            diag_mask,
                            h,
                            pos,
                            reward,
                            penalty,
                            x_dropoff,
                            kbp,
                            search_space,
                            evalue_threshold,
                            context,
                            hsps,
                        );
                    }
                }
            };
        }

        check_pos!(base_shift, 0); // Position 0: shift by base_shift
        check_pos!(base_shift - 2, 1); // Position 1: shift by base_shift - 2
        check_pos!(base_shift - 4, 2); // Position 2: shift by base_shift - 4
        check_pos!(base_shift - 6, 3); // Position 3: shift by base_shift - 6

        s_pos += 4;
        byte_idx += 1;
    }

    // Handle remaining positions (less than 4)
    while s_pos < end {
        let h = packed_hash_at(subject_packed, s_pos, lut_word, lut_mask);
        if pv[h >> 6] & (1u64 << (h & 63)) != 0 {
            process_hit(
                query,
                subject_packed,
                subject_len,
                word_size,
                lut_word,
                lut,
                next,
                last_hit,
                diag_mask,
                h,
                s_pos,
                reward,
                penalty,
                x_dropoff,
                kbp,
                search_space,
                evalue_threshold,
                context,
                hsps,
            );
        }
        s_pos += 1;
    }
}

#[allow(clippy::too_many_arguments)]
fn scan_byte_oriented_lut6_step2(
    subject_packed: &[u8],
    subject_len: usize,
    end: usize,
    lut_mask: u32,
    query: &[u8],
    word_size: usize,
    lut: &[i32],
    next: &[i32],
    pv: &[u64],
    last_hit: &mut [i32],
    diag_mask: usize,
    reward: i32,
    penalty: i32,
    x_dropoff: i32,
    kbp: &KarlinBlk,
    search_space: f64,
    evalue_threshold: f64,
    context: i32,
    hsps: &mut Vec<SearchHsp>,
) {
    let s_ptr = subject_packed.as_ptr();
    let s_len = subject_packed.len();
    let mut byte_idx = 0usize;
    let mut s_pos = 1usize;

    macro_rules! check_hash {
        ($combined:expr, $shift:expr, $pos:expr) => {
            if $pos < end {
                let h = (($combined >> $shift) & lut_mask) as usize;
                if pv_has_hash(pv, h) {
                    process_hit(
                        query,
                        subject_packed,
                        subject_len,
                        word_size,
                        6,
                        lut,
                        next,
                        last_hit,
                        diag_mask,
                        h,
                        $pos,
                        reward,
                        penalty,
                        x_dropoff,
                        kbp,
                        search_space,
                        evalue_threshold,
                        context,
                        hsps,
                    );
                }
            }
        };
    }

    while s_pos + 2 < end && byte_idx + 2 < s_len {
        let combined = unsafe {
            ((*s_ptr.add(byte_idx) as u32) << 16)
                | ((*s_ptr.add(byte_idx + 1) as u32) << 8)
                | (*s_ptr.add(byte_idx + 2) as u32)
        };

        check_hash!(combined, 10, s_pos);
        check_hash!(combined, 6, s_pos + 2);

        byte_idx += 1;
        s_pos += 4;
    }

    while s_pos < end {
        let h = packed_hash_at(subject_packed, s_pos, 6, lut_mask);
        if pv[h >> 6] & (1u64 << (h & 63)) != 0 {
            process_hit(
                query,
                subject_packed,
                subject_len,
                word_size,
                6,
                lut,
                next,
                last_hit,
                diag_mask,
                h,
                s_pos,
                reward,
                penalty,
                x_dropoff,
                kbp,
                search_space,
                evalue_threshold,
                context,
                hsps,
            );
        }
        s_pos += 2;
    }
}

#[allow(clippy::too_many_arguments)]
fn scan_byte_oriented_lut6_step2_pair(
    subject_packed: &[u8],
    subject_len: usize,
    end: usize,
    lookup0: &NaLookup<'_>,
    last_hit0: &mut [i32],
    lookup1: &NaLookup<'_>,
    last_hit1: &mut [i32],
    paired_pv: &[u64],
    word_size: usize,
    reward: i32,
    penalty: i32,
    x_dropoff: i32,
    kbp: &KarlinBlk,
    search_space: f64,
    evalue_threshold: f64,
    hsps: &mut Vec<SearchHsp>,
) {
    let s_ptr = subject_packed.as_ptr();
    let s_len = subject_packed.len();
    let lut_mask = lookup0.lut_mask;
    debug_assert_eq!(lut_mask, lookup1.lut_mask);
    debug_assert_eq!(paired_pv.len(), lookup0.pv.len());
    debug_assert_eq!(word_size, 7);

    let mut byte_idx = 0usize;
    let mut s_pos = 1usize;

    macro_rules! check_hash_pair {
        ($combined:expr, $shift:expr, $pos:expr) => {
            if $pos < end {
                let h = (($combined >> $shift) & lut_mask) as usize;
                let pv_bucket = h >> 6;
                let pv_bit = 1u64 << (h & 63);
                if pv_bucket_unchecked(paired_pv, pv_bucket) & pv_bit != 0 {
                    if pv_bucket_unchecked(&lookup0.pv, pv_bucket) & pv_bit != 0 {
                        process_hit(
                            lookup0.query,
                            subject_packed,
                            subject_len,
                            word_size,
                            6,
                            &lookup0.lut,
                            &lookup0.next,
                            last_hit0,
                            lookup0.diag_mask,
                            h,
                            $pos,
                            reward,
                            penalty,
                            x_dropoff,
                            kbp,
                            search_space,
                            evalue_threshold,
                            lookup0.context,
                            hsps,
                        );
                    }
                    if pv_bucket_unchecked(&lookup1.pv, pv_bucket) & pv_bit != 0 {
                        process_hit(
                            lookup1.query,
                            subject_packed,
                            subject_len,
                            word_size,
                            6,
                            &lookup1.lut,
                            &lookup1.next,
                            last_hit1,
                            lookup1.diag_mask,
                            h,
                            $pos,
                            reward,
                            penalty,
                            x_dropoff,
                            kbp,
                            search_space,
                            evalue_threshold,
                            lookup1.context,
                            hsps,
                        );
                    }
                }
            }
        };
    }

    if s_pos + 2 < end && byte_idx + 2 < s_len {
        let mut combined = unsafe {
            ((*s_ptr.add(byte_idx) as u32) << 16)
                | ((*s_ptr.add(byte_idx + 1) as u32) << 8)
                | (*s_ptr.add(byte_idx + 2) as u32)
        };

        loop {
            check_hash_pair!(combined, 10, s_pos);
            check_hash_pair!(combined, 6, s_pos + 2);

            byte_idx += 1;
            s_pos += 4;

            if !(s_pos + 2 < end && byte_idx + 2 < s_len) {
                break;
            }

            combined = ((combined & 0xffff) << 8) | unsafe { *s_ptr.add(byte_idx + 2) as u32 };
        }
    }

    while s_pos < end {
        let h = packed_hash_at(subject_packed, s_pos, 6, lut_mask);
        let pv_bucket = h >> 6;
        let pv_bit = 1u64 << (h & 63);
        if paired_pv[pv_bucket] & pv_bit != 0 {
            if lookup0.pv[pv_bucket] & pv_bit != 0 {
                process_hit(
                    lookup0.query,
                    subject_packed,
                    subject_len,
                    word_size,
                    6,
                    &lookup0.lut,
                    &lookup0.next,
                    last_hit0,
                    lookup0.diag_mask,
                    h,
                    s_pos,
                    reward,
                    penalty,
                    x_dropoff,
                    kbp,
                    search_space,
                    evalue_threshold,
                    lookup0.context,
                    hsps,
                );
            }
            if lookup1.pv[pv_bucket] & pv_bit != 0 {
                process_hit(
                    lookup1.query,
                    subject_packed,
                    subject_len,
                    word_size,
                    6,
                    &lookup1.lut,
                    &lookup1.next,
                    last_hit1,
                    lookup1.diag_mask,
                    h,
                    s_pos,
                    reward,
                    penalty,
                    x_dropoff,
                    kbp,
                    search_space,
                    evalue_threshold,
                    lookup1.context,
                    hsps,
                );
            }
        }
        s_pos += 2;
    }
}

/// Compute hash for a word starting at position `pos` in packed data.
#[inline(always)]
fn packed_hash_at(packed: &[u8], pos: usize, lut_word: usize, lut_mask: u32) -> usize {
    let mut h: u32 = 0;
    for i in 0..lut_word {
        h = (h << 2) | packed_base_at(packed, pos + i) as u32;
    }
    (h & lut_mask) as usize
}

/// Fallback step-1 scanning (every position). Used when scan_step != 4.
#[allow(clippy::too_many_arguments)]
fn scan_step1(
    subject_packed: &[u8],
    _subject_len: usize,
    end: usize,
    lut_word: usize,
    lut_mask: u32,
    scan_start: usize,
    scan_step: usize,
    query: &[u8],
    word_size: usize,
    lut: &[i32],
    next: &[i32],
    pv: &[u64],
    last_hit: &mut [i32],
    diag_mask: usize,
    reward: i32,
    penalty: i32,
    x_dropoff: i32,
    kbp: &KarlinBlk,
    search_space: f64,
    evalue_threshold: f64,
    context: i32,
    hsps: &mut Vec<SearchHsp>,
) {
    let subject_len = _subject_len;
    if scan_step != 1 || scan_start != 0 {
        let mut s_pos = scan_start;
        while s_pos < end {
            let h = packed_hash_at(subject_packed, s_pos, lut_word, lut_mask);
            if pv[h >> 6] & (1u64 << (h & 63)) != 0 {
                process_hit(
                    query,
                    subject_packed,
                    subject_len,
                    word_size,
                    lut_word,
                    lut,
                    next,
                    last_hit,
                    diag_mask,
                    h,
                    s_pos,
                    reward,
                    penalty,
                    x_dropoff,
                    kbp,
                    search_space,
                    evalue_threshold,
                    context,
                    hsps,
                );
            }
            s_pos += scan_step;
        }
        return;
    }

    let mut hash: u32 = 0;
    for i in 0..(lut_word - 1).min(7) {
        hash = (hash << 2) | packed_base_at(subject_packed, i) as u32;
    }
    let mut s_pos = 0;
    while s_pos < end {
        hash =
            ((hash << 2) | packed_base_at(subject_packed, s_pos + lut_word - 1) as u32) & lut_mask;
        let h = hash as usize;
        if pv[h >> 6] & (1u64 << (h & 63)) != 0 {
            process_hit(
                query,
                subject_packed,
                subject_len,
                word_size,
                lut_word,
                lut,
                next,
                last_hit,
                diag_mask,
                h,
                s_pos,
                reward,
                penalty,
                x_dropoff,
                kbp,
                search_space,
                evalue_threshold,
                context,
                hsps,
            );
        }
        s_pos += 1;
    }
}

/// Extract a single base from packed NCBI2na data (4 bases per byte).
#[inline(always)]
fn packed_base_at(packed: &[u8], pos: usize) -> u8 {
    let byte = packed[pos >> 2];
    (byte >> (6 - 2 * (pos & 3))) & 3
}

#[inline(always)]
fn blastna_to_iupac_byte(b: u8) -> u8 {
    const IUPAC: &[u8; 16] = b"ACGTRYMKWSBDHVN-";
    if b < 16 {
        // SAFETY: guarded by the range check above.
        unsafe { *IUPAC.get_unchecked(b as usize) }
    } else {
        b'-'
    }
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
    // Extend right using NCBI's exact nucleotide ungapped algorithm:
    // accumulate a temporary run score, fold positive runs into the HSP score,
    // and stop once the current run drops below max(-xdrop, -score).
    let mut score = 0i32;
    let mut sum = 0i32;
    let mut best_right = 0usize;
    let mut qi = q_seed;
    let mut si = s_seed;
    let x_dropoff_neg = -x_dropoff;
    let mut x_current = x_dropoff_neg;
    while qi < query.len() && si < subject_len {
        let sb = packed_base_at(subject_packed, si);
        sum += blastna_score(query[qi], sb, reward, penalty);
        if sum > 0 {
            score += sum;
            best_right = qi - q_seed + 1;
            x_current = (-score).max(x_dropoff_neg);
            sum = 0;
        } else if sum < x_current {
            break;
        }
        qi += 1;
        si += 1;
    }

    // Extend left
    let mut left_score = 0i32;
    let mut left_sum = 0i32;
    let mut best_left = 0usize;
    if q_seed > 0 && s_seed > 0 {
        let mut qi = q_seed - 1;
        let mut si = s_seed - 1;
        loop {
            let sb = packed_base_at(subject_packed, si);
            left_sum += blastna_score(query[qi], sb, reward, penalty);
            if left_sum > 0 {
                left_score += left_sum;
                best_left = q_seed - qi;
                left_sum = 0;
            } else if left_sum < x_dropoff_neg {
                break;
            }
            if qi == 0 || si == 0 {
                break;
            }
            qi -= 1;
            si -= 1;
        }
    }

    let total_score = score + left_score;
    if total_score <= 0 {
        return None;
    }

    let evalue = kbp.raw_to_evalue(total_score, search_space);
    if evalue > evalue_threshold {
        return None;
    }

    let q_start = (q_seed - best_left) as i32;
    let q_end = (q_seed + best_right) as i32;
    let s_start = (s_seed - best_left) as i32;
    let s_end = (s_seed + best_right) as i32;
    let align_len = q_end - q_start;

    // Count identities and build aligned sequences
    let mut num_ident = 0i32;
    let mut qseq = Vec::with_capacity(align_len as usize);
    let mut sseq = Vec::with_capacity(align_len as usize);
    for i in 0..align_len as usize {
        let qb = query[q_start as usize + i];
        let sb = packed_base_at(subject_packed, s_start as usize + i);
        if qb == sb {
            num_ident += 1;
        }
        qseq.push(blastna_to_iupac_byte(qb));
        sseq.push(blastna_to_iupac_byte(sb));
    }
    // SAFETY: blastna_to_iupac_byte emits ASCII only.
    let qseq_str = unsafe { String::from_utf8_unchecked(qseq) };
    let sseq_str = unsafe { String::from_utf8_unchecked(sseq) };

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
    // Extend right using NCBI's exact nucleotide ungapped algorithm:
    // accumulate a temporary run score, fold positive runs into the HSP score,
    // and stop once the current run drops below max(-xdrop, -score).
    let mut score = 0i32;
    let mut sum = 0i32;
    let mut best_right = 0usize;
    let mut qi = q_seed;
    let mut si = s_seed;
    let x_dropoff_neg = -x_dropoff;
    let mut x_current = x_dropoff_neg;
    while qi < query.len() && si < subject.len() {
        let sb = subject[si];
        if subject_base_blocks_ungapped_extension(sb) {
            break;
        }
        sum += blastna_score(query[qi], sb, reward, penalty);
        if sum > 0 {
            score += sum;
            best_right = qi - q_seed + 1;
            x_current = (-score).max(x_dropoff_neg);
            sum = 0;
        } else if sum < x_current {
            break;
        }
        qi += 1;
        si += 1;
    }

    // Extend left
    let mut score_l = 0i32;
    let mut sum_l = 0i32;
    let mut best_left = 0usize;
    if q_seed > 0 && s_seed > 0 {
        qi = q_seed - 1;
        si = s_seed - 1;
        loop {
            let sb = subject[si];
            if subject_base_blocks_ungapped_extension(sb) {
                break;
            }
            sum_l += blastna_score(query[qi], sb, reward, penalty);
            if sum_l > 0 {
                score_l += sum_l;
                best_left = q_seed - qi;
                sum_l = 0;
            } else if sum_l < x_dropoff_neg {
                break;
            }
            if qi == 0 || si == 0 {
                break;
            }
            qi -= 1;
            si -= 1;
        }
    }

    let total_score = score + score_l;
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
    let mut qseq = Vec::with_capacity(align_len as usize);
    let mut sseq = Vec::with_capacity(align_len as usize);
    for i in 0..align_len as usize {
        let qb = query[q_start as usize + i];
        let sb = subject[s_start as usize + i];
        if qb == sb {
            num_ident += 1;
        }
        qseq.push(blastna_to_iupac_byte(qb));
        sseq.push(blastna_to_iupac_byte(sb));
    }
    // SAFETY: blastna_to_iupac_byte emits ASCII only.
    let qseq_str = unsafe { String::from_utf8_unchecked(qseq) };
    let sseq_str = unsafe { String::from_utf8_unchecked(sseq) };

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
fn dedup_hsps_with_min_diag_separation(hsps: &mut Vec<SearchHsp>, min_diag_separation: i32) {
    if hsps.len() <= 1 {
        return;
    }
    hsps.sort_by(score_compare_search_hsps);

    // Separate trees per context (strand)
    let q_max = hsps.iter().map(|h| h.query_end).max().unwrap_or(0) + 1;
    let s_max = hsps.iter().map(|h| h.subject_end).max().unwrap_or(0) + 1;
    let mut trees: std::collections::HashMap<i32, IntervalTree> = std::collections::HashMap::new();

    let mut keep = vec![false; hsps.len()];
    for (i, hsp) in hsps.iter().enumerate() {
        let tree = trees
            .entry(hsp.context)
            .or_insert_with(|| IntervalTree::new(q_max, s_max));
        let qi = Interval::new(hsp.query_start, hsp.query_end);
        let si = Interval::new(hsp.subject_start, hsp.subject_end);
        if !tree.is_contained_with_min_diag_separation(qi, si, min_diag_separation)
            && !tree.has_significant_overlap_with_min_diag_separation(
                qi,
                si,
                0.99,
                min_diag_separation,
            )
        {
            keep[i] = true;
            tree.insert(qi, si, hsp.score);
        }
    }
    let mut idx = 0;
    hsps.retain(|_| {
        let k = keep[idx];
        idx += 1;
        k
    });
}

fn min_diag_separation_for_ungapped(word_size: usize, reward: i32, penalty: i32) -> i32 {
    if word_size >= 28 && reward == 1 && penalty == -2 {
        6
    } else {
        0
    }
}

#[inline(always)]
fn has_ambiguous_base(word: &[u8]) -> bool {
    word.iter().any(|&b| b >= 4)
}

#[inline(always)]
fn subject_base_blocks_ungapped_extension(b: u8) -> bool {
    b >= 4
}

fn blastna_score(a: u8, b: u8, reward: i32, penalty: i32) -> i32 {
    const MASKS: [u8; 16] = [1, 2, 4, 8, 5, 10, 3, 12, 9, 6, 14, 13, 11, 7, 15, 0];
    const DEGEN: [i32; 16] = [1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 0];

    let a = a as usize;
    let b = b as usize;
    if a >= 15 || b >= 15 {
        return i32::MIN / 2;
    }
    if MASKS[a] & MASKS[b] == 0 {
        return penalty;
    }
    let degen = DEGEN[a.max(b)];
    // NCBI `blast_stat.c:1106`:
    //   (Int4)BLAST_Nint((double)((degeneracy[i2]-1)*penalty + reward) / degeneracy[i2])
    crate::math::nint(((degen - 1) * penalty + reward) as f64 / degen as f64) as i32
}

use crate::stat::HSP_MAX_WINDOW;

/// Pick the gapped-alignment seed for a nucleotide ungapped HSP.
///
/// Port of `BlastGetOffsetsForGappedAlignment`
/// (`blast_gapalign.c:3248`): slide an `HSP_MAX_WINDOW`-wide window of
/// match/mismatch scores across the HSP, return the **strict** max-score
/// window's right edge. For uniform-score HSPs (perfect matches) the first
/// window wins and the returned offset is `q_start + HSP_MAX_WINDOW - 1`;
/// that matches NCBI's tie-break (NCBI uses `if (score > max_score)` with
/// strict `>`, so the first maximum is kept).
///
/// Returns `(q_off, s_off)` with the same semantics as NCBI's
/// `*q_retval, *s_retval`.
///
/// NCBI also has a short-HSP fast path (`q_length <= HSP_MAX_WINDOW` →
/// midpoint) and a `score > 0` fallback that checks the end window. Both
/// are ported verbatim. Ambiguity handling matches NCBI's `matrix[q][s]`
/// lookup (treating BLASTNA ambiguous codes via their matrix scores).
fn blast_get_offsets_for_gapped_alignment(
    query: &[u8],
    subject: &[u8],
    q_start: usize,
    q_end: usize,
    s_start: usize,
    s_end: usize,
    reward: i32,
    penalty: i32,
) -> Option<(usize, usize)> {
    let q_length = q_end.saturating_sub(q_start);
    let s_length = s_end.saturating_sub(s_start);
    if q_length == 0 || s_length == 0 {
        return None;
    }
    // Short-HSP fast path (blast_gapalign.c:3259-3263).
    if q_length <= HSP_MAX_WINDOW {
        let mid = q_length / 2;
        return Some((q_start + mid, s_start + mid));
    }

    // Initial window [q_start .. q_start + HSP_MAX_WINDOW).
    let mut score: i32 = 0;
    for i in 0..HSP_MAX_WINDOW {
        let q = query[q_start + i];
        let s = subject[s_start + i];
        score += blastna_score(q, s, reward, penalty);
    }
    let mut max_score = score;
    let mut max_offset = q_start + HSP_MAX_WINDOW - 1;

    // Slide the window to the end of the HSP.
    let hsp_end = q_start + q_length.min(s_length);
    for index1 in (q_start + HSP_MAX_WINDOW)..hsp_end {
        let q_out = query[index1 - HSP_MAX_WINDOW];
        let s_out = subject[(index1 - HSP_MAX_WINDOW) - q_start + s_start];
        let q_in = query[index1];
        let s_in = subject[index1 - q_start + s_start];
        score -= blastna_score(q_out, s_out, reward, penalty);
        score += blastna_score(q_in, s_in, reward, penalty);
        // Strict > tie-break, matches NCBI blast_gapalign.c:3287.
        if score > max_score {
            max_score = score;
            max_offset = index1;
        }
    }

    if max_score > 0 {
        return Some((max_offset, (max_offset - q_start) + s_start));
    }

    // Fallback: test the window at the end of the HSP
    // (blast_gapalign.c:3300-3318).
    let end_score: i32 = (0..HSP_MAX_WINDOW)
        .map(|i| {
            let q = query[q_end - HSP_MAX_WINDOW + i];
            let s = subject[s_end - HSP_MAX_WINDOW + i];
            blastna_score(q, s, reward, penalty)
        })
        .sum();
    if end_score > 0 {
        Some((q_end - HSP_MAX_WINDOW / 2, s_end - HSP_MAX_WINDOW / 2))
    } else {
        None
    }
}

fn min_diag_separation_for_gapped(
    word_size: usize,
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
) -> i32 {
    // NCBI `CBlastNucleotideOptionsHandle::SetMBHitSavingOptionsDefaults`
    // (`blast_nucl_options.cpp:259`) sets min_diag_separation = 6 for
    // megablast. Regular blastn (`SetHitSavingOptionsDefaults:239`)
    // uses 50 but Rust currently only applies the megablast value on
    // the matching scoring-system fingerprint.
    if word_size >= 28 && reward == 1 && penalty == -2 && gap_open == 0 && gap_extend == 0 {
        6
    } else {
        0
    }
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
        query_plus,
        query_minus,
        query_plus,
        query_minus,
        subject,
        word_size,
        reward,
        penalty,
        _gap_open,
        _gap_extend,
        x_dropoff,
        kbp,
        search_space,
        evalue_threshold,
    )
}

/// Gapped search with separate xdrop values for preliminary and final traceback.
/// Passes `max(x_dropoff, x_dropoff_final)` to the traceback DP, matching NCBI
/// `blast_parameters.c:462-463` which clamps the final x-drop to at least the
/// preliminary value before handing it to `ALIGN_EX`. The CLI's caller still
/// hardcodes the raw final value until the linear-gap lambda lookup matches
/// NCBI for non-affine scoring systems (see `main.rs` lambda-table note).
#[allow(clippy::too_many_arguments)]
pub fn blastn_gapped_search_nomask_with_xdrops(
    query_plus: &[u8],
    query_minus: &[u8],
    query_plus_nomask: &[u8],
    query_minus_nomask: &[u8],
    subject: &[u8],
    word_size: usize,
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    x_dropoff: i32,
    x_dropoff_final: i32,
    kbp: &KarlinBlk,
    search_space: f64,
    evalue_threshold: f64,
) -> Vec<SearchHsp> {
    // NCBI's `BlastExtensionParametersNew` guarantees
    // `gap_x_dropoff_final >= gap_x_dropoff`; traceback (`ALIGN_EX` in
    // `blast_traceback.c`) uses the final value. Pass the final x_dropoff
    // through so the traceback DP has the same exploration budget as NCBI.
    // Preliminary score-only filtering (separate pass, inside
    // `passes_preliminary_gapped_score`) still uses the stricter
    // `x_dropoff`.
    let effective = x_dropoff_final.max(x_dropoff);
    blastn_gapped_search_nomask(
        query_plus,
        query_minus,
        query_plus_nomask,
        query_minus_nomask,
        subject,
        word_size,
        reward,
        penalty,
        gap_open,
        gap_extend,
        effective,
        kbp,
        search_space,
        evalue_threshold,
    )
}

/// Gapped search with separate masked (for seeds) and unmasked (for gapped alignment) queries.
pub fn blastn_gapped_search_nomask(
    query_plus: &[u8],         // masked, for seed finding
    query_minus: &[u8],        // masked, for seed finding
    query_plus_nomask: &[u8],  // unmasked, for gapped alignment
    query_minus_nomask: &[u8], // unmasked, for gapped alignment
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
        query_plus,
        query_minus,
        subject,
        word_size,
        reward,
        penalty,
        x_dropoff,
        kbp,
        search_space,
        evalue_threshold * 100.0, // permissive threshold for seeds
    );

    let mut hsps = Vec::new();
    let min_diag_separation =
        min_diag_separation_for_gapped(word_size, reward, penalty, _gap_open, _gap_extend);

    for seed in &ungapped {
        // Use UNMASKED query for gapped alignment (matches C engine's sequence_nomask)
        let query = if seed.context == 0 {
            query_plus_nomask
        } else {
            query_minus_nomask
        };

        if min_diag_separation != 0 && is_perfect_ungapped_hsp(seed, reward) {
            hsps.push(seed.clone());
            continue;
        }

        // Seed position for gapped extension. Mirrors NCBI's nucleotide path in
        // blast_traceback.c:436-461: when `gapped_start` is unset (0,0), call
        // `BlastGetOffsetsForGappedAlignment` to slide an 11-wide window and
        // pick the max-score offset. Falls back to the midpoint if the seed
        // picker returns None (extremely short HSPs or all-ambiguous regions).
        let (seed_q, seed_s) = blast_get_offsets_for_gapped_alignment(
            query,
            subject,
            seed.query_start as usize,
            seed.query_end as usize,
            seed.subject_start as usize,
            seed.subject_end as usize,
            reward,
            penalty,
        )
        .unwrap_or_else(|| {
            (
                ((seed.query_start + seed.query_end) / 2) as usize,
                ((seed.subject_start + seed.subject_end) / 2) as usize,
            )
        });
        if !passes_preliminary_gapped_score(
            query,
            subject,
            seed_q,
            seed_s,
            reward,
            penalty,
            _gap_open,
            _gap_extend,
            x_dropoff,
            kbp,
            search_space,
            evalue_threshold,
        ) {
            continue;
        }

        if let Some(mut tb) = blast_gapped_align(
            query,
            subject,
            seed_q,
            seed_s,
            reward,
            penalty,
            _gap_open,
            _gap_extend,
            x_dropoff,
        ) {
            // NCBI runs `Blast_HSPReevaluateWithAmbiguitiesGapped` after
            // `ALIGN_EX` (`blast_traceback.c:657`). The cutoff_score comes
            // from `hit_params->cutoffs[context].cutoff_score` which is the
            // e-value-based cutoff derived from KBP, search_space, and the
            // user's evalue threshold (no decay adjustment — NCBI passes
            // `dodecay=FALSE` at `blast_parameters.c:943`).
            let cutoff_score = kbp.evalue_to_raw(evalue_threshold, search_space);
            let _delete = reevaluate_with_ambiguities_gapped(
                &mut tb,
                query,
                subject,
                reward,
                penalty,
                _gap_open,
                _gap_extend,
                cutoff_score.max(1),
            );
            if _delete {
                continue;
            }

            let evalue = kbp.raw_to_evalue(tb.score, search_space);
            if evalue > evalue_threshold {
                continue;
            }

            let q_slice = &query[tb.query_start..tb.query_end];
            let s_slice = &subject[tb.subject_start..tb.subject_end];
            let (align_len, num_ident, gap_opens) =
                tb.edit_script.count_identities(q_slice, s_slice);
            let (qseq, sseq) = tb
                .edit_script
                .render_alignment(q_slice, s_slice, blastna_to_iupac);

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

    hsps.sort_by(score_compare_search_hsps);
    dedup_hsps_with_min_diag_separation(&mut hsps, min_diag_separation);
    hsps
}

#[allow(clippy::too_many_arguments)]
fn blastn_gapped_search_decoded_prepared_with_scratch(
    prepared: &PreparedBlastnQuery<'_>,
    query_plus_nomask: &[u8],
    query_minus_nomask: &[u8],
    subject: &[u8],
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    x_dropoff: i32,
    kbp: &KarlinBlk,
    search_space: f64,
    evalue_threshold: f64,
    last_hit_scratch: &mut [Vec<i32>],
) -> Vec<SearchHsp> {
    let ungapped = blastn_ungapped_search_decoded_prepared_with_scratch_no_dedup(
        prepared,
        subject,
        reward,
        penalty,
        x_dropoff,
        kbp,
        search_space,
        evalue_threshold * 100.0,
        last_hit_scratch,
    );
    let mut hsps = Vec::new();
    let min_diag_separation =
        min_diag_separation_for_gapped(prepared.word_size, reward, penalty, gap_open, gap_extend);

    for seed in &ungapped {
        let query = if seed.context == 0 {
            query_plus_nomask
        } else {
            query_minus_nomask
        };

        if min_diag_separation != 0 && is_perfect_ungapped_hsp(seed, reward) {
            hsps.push(seed.clone());
            continue;
        }

        // Same NCBI-style seed picker as the non-prepared variant above.
        let (seed_q, seed_s) = blast_get_offsets_for_gapped_alignment(
            query,
            subject,
            seed.query_start as usize,
            seed.query_end as usize,
            seed.subject_start as usize,
            seed.subject_end as usize,
            reward,
            penalty,
        )
        .unwrap_or_else(|| {
            (
                ((seed.query_start + seed.query_end) / 2) as usize,
                ((seed.subject_start + seed.subject_end) / 2) as usize,
            )
        });

        if !passes_preliminary_gapped_score(
            query,
            subject,
            seed_q,
            seed_s,
            reward,
            penalty,
            gap_open,
            gap_extend,
            x_dropoff,
            kbp,
            search_space,
            evalue_threshold,
        ) {
            continue;
        }

        if let Some(mut tb) = blast_gapped_align(
            query, subject, seed_q, seed_s, reward, penalty, gap_open, gap_extend, x_dropoff,
        ) {
            // Same reevaluate pass as the non-prepared variant —
            // see `blast_traceback.c:657`. Matches NCBI `BLAST_Cutoffs`
            // with `dodecay=FALSE` (`blast_parameters.c:943`).
            let cutoff_score = kbp.evalue_to_raw(evalue_threshold, search_space);
            let _delete = reevaluate_with_ambiguities_gapped(
                &mut tb,
                query,
                subject,
                reward,
                penalty,
                gap_open,
                gap_extend,
                cutoff_score.max(1),
            );
            if _delete {
                continue;
            }

            let evalue = kbp.raw_to_evalue(tb.score, search_space);
            if evalue > evalue_threshold {
                continue;
            }

            let q_slice = &query[tb.query_start..tb.query_end];
            let s_slice = &subject[tb.subject_start..tb.subject_end];
            let (align_len, num_ident, gap_opens) =
                tb.edit_script.count_identities(q_slice, s_slice);
            let (qseq, sseq) = tb
                .edit_script
                .render_alignment(q_slice, s_slice, blastna_to_iupac);

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

    hsps.sort_by(score_compare_search_hsps);
    dedup_hsps_with_min_diag_separation(&mut hsps, min_diag_separation);
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
    // Keep the packed scan as a no-hit probe, but use the decoded seed path
    // for traceback so gapped DB mode sees the same seed shapes as subject mode.
    let ungapped_probe = blastn_ungapped_search_packed(
        query_plus,
        query_minus,
        subject_packed,
        subject_len,
        word_size,
        reward,
        penalty,
        x_dropoff,
        kbp,
        search_space,
        evalue_threshold * 100.0,
    );
    if ungapped_probe.is_empty() {
        return Vec::new();
    }

    let subject_decoded = decode_packed_ncbi2na(subject_packed, subject_len);
    blastn_gapped_search_nomask(
        query_plus,
        query_minus,
        query_plus_nomask,
        query_minus_nomask,
        &subject_decoded,
        word_size,
        reward,
        penalty,
        gap_open,
        gap_extend,
        x_dropoff,
        kbp,
        search_space,
        evalue_threshold,
    )
}

/// Variant of `blastn_gapped_search_packed_prepared` that accepts both the
/// preliminary and final gapped x-drop values. Pass the final value through
/// (`max(x_dropoff, x_dropoff_final)`) to match NCBI's `ALIGN_EX` traceback
/// exploration budget.
#[allow(clippy::too_many_arguments)]
pub fn blastn_gapped_search_packed_prepared_with_xdrops(
    prepared: &PreparedBlastnQuery<'_>,
    query_plus_nomask: &[u8],
    query_minus_nomask: &[u8],
    subject_packed: &[u8],
    subject_len: usize,
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    x_dropoff: i32,
    x_dropoff_final: i32,
    kbp: &KarlinBlk,
    search_space: f64,
    evalue_threshold: f64,
    last_hit_scratch: &mut [Vec<i32>],
) -> Vec<SearchHsp> {
    let effective = x_dropoff_final.max(x_dropoff);
    blastn_gapped_search_packed_prepared(
        prepared,
        query_plus_nomask,
        query_minus_nomask,
        subject_packed,
        subject_len,
        reward,
        penalty,
        gap_open,
        gap_extend,
        effective,
        kbp,
        search_space,
        evalue_threshold,
        last_hit_scratch,
    )
}

pub fn blastn_gapped_search_packed_prepared(
    prepared: &PreparedBlastnQuery<'_>,
    query_plus_nomask: &[u8],
    query_minus_nomask: &[u8],
    subject_packed: &[u8],
    subject_len: usize,
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    x_dropoff: i32,
    kbp: &KarlinBlk,
    search_space: f64,
    evalue_threshold: f64,
    last_hit_scratch: &mut [Vec<i32>],
) -> Vec<SearchHsp> {
    let ungapped_probe = blastn_ungapped_search_packed_prepared_with_scratch(
        prepared,
        subject_packed,
        subject_len,
        reward,
        penalty,
        x_dropoff,
        kbp,
        search_space,
        evalue_threshold * 100.0,
        last_hit_scratch,
    );
    if ungapped_probe.is_empty() {
        return Vec::new();
    }

    let subject_decoded = decode_packed_ncbi2na(subject_packed, subject_len);
    blastn_gapped_search_decoded_prepared_with_scratch(
        prepared,
        query_plus_nomask,
        query_minus_nomask,
        &subject_decoded,
        reward,
        penalty,
        gap_open,
        gap_extend,
        x_dropoff,
        kbp,
        search_space,
        evalue_threshold,
        last_hit_scratch,
    )
}

#[inline]
fn passes_preliminary_gapped_score(
    query: &[u8],
    subject: &[u8],
    seed_q: usize,
    seed_s: usize,
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    x_dropoff: i32,
    kbp: &KarlinBlk,
    search_space: f64,
    evalue_threshold: f64,
) -> bool {
    let prelim_score = blast_gapped_score_only(
        query, subject, seed_q, seed_s, reward, penalty, gap_open, gap_extend, x_dropoff,
    );
    prelim_score > 0 && kbp.raw_to_evalue(prelim_score, search_space) <= evalue_threshold
}

#[inline]
fn is_perfect_ungapped_hsp(seed: &SearchHsp, reward: i32) -> bool {
    reward > 0
        && seed.align_length > 0
        && seed.num_ident == seed.align_length
        && seed.mismatches == 0
        && seed.gap_opens == 0
        && seed.score as i64 == reward as i64 * seed.align_length as i64
}

/// Decode packed NCBI2na to per-base (0=A, 1=C, 2=G, 3=T).
pub fn decode_packed_ncbi2na(packed: &[u8], len: usize) -> Vec<u8> {
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if len >= 64 && std::arch::is_x86_feature_detected!("ssse3") {
            // SAFETY: guarded by runtime SSSE3 detection; the SIMD routine
            // only reads complete 16-byte chunks and handles the tail scalar.
            return unsafe { decode_packed_ncbi2na_ssse3(packed, len) };
        }
    }

    decode_packed_ncbi2na_scalar(packed, len)
}

fn decode_packed_ncbi2na_scalar(packed: &[u8], len: usize) -> Vec<u8> {
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

/// Decode packed NCBI2na and overlay BLAST DB ambiguity descriptors.
pub fn decode_packed_ncbi2na_with_ambiguity(
    packed: &[u8],
    len: usize,
    ambiguity_data: &[u8],
) -> Vec<u8> {
    let mut decoded = decode_packed_ncbi2na(packed, len);
    overlay_ncbi4na_ambiguity(&mut decoded, ambiguity_data);
    decoded
}

fn overlay_ncbi4na_ambiguity(decoded: &mut [u8], ambiguity_data: &[u8]) {
    if ambiguity_data.len() < 4 {
        return;
    }

    let words: Vec<u32> = ambiguity_data
        .chunks_exact(4)
        .map(|chunk| u32::from_be_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]))
        .collect();
    if words.is_empty() {
        return;
    }

    let mut amb_num = words[0];
    let new_format = (amb_num & 0x8000_0000) != 0;
    if new_format {
        amb_num &= 0x7fff_ffff;
    }

    let mut i = 1usize;
    let mut seen = 0u32;
    while seen < amb_num && i < words.len() {
        let word = words[i];
        let ncbi4na = ((word >> 28) & 0x0f) as usize;
        let run_len_minus_one: usize;
        let position: usize;

        if new_format {
            if i + 1 >= words.len() {
                break;
            }
            run_len_minus_one = ((word >> 16) & 0x0fff) as usize;
            position = words[i + 1] as usize;
            i += 2;
        } else {
            run_len_minus_one = ((word >> 24) & 0x0f) as usize;
            position = (word & 0x00ff_ffff) as usize;
            i += 1;
        }

        if ncbi4na < NCBI4NA_TO_BLASTNA.len() {
            let blastna = NCBI4NA_TO_BLASTNA[ncbi4na];
            let run_len = run_len_minus_one + 1;
            let end = position.saturating_add(run_len).min(decoded.len());
            if position < end {
                decoded[position..end].fill(blastna);
            }
        }
        seen += 1;
    }
}

#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::{
    __m128i, _mm_and_si128, _mm_loadu_si128, _mm_set1_epi8, _mm_set_epi8, _mm_shuffle_epi8,
    _mm_srli_epi16, _mm_storeu_si128, _mm_unpackhi_epi16, _mm_unpackhi_epi8, _mm_unpacklo_epi16,
    _mm_unpacklo_epi8,
};

#[cfg(target_arch = "x86")]
use std::arch::x86::{
    __m128i, _mm_and_si128, _mm_loadu_si128, _mm_set1_epi8, _mm_set_epi8, _mm_shuffle_epi8,
    _mm_srli_epi16, _mm_storeu_si128, _mm_unpackhi_epi16, _mm_unpackhi_epi8, _mm_unpacklo_epi16,
    _mm_unpacklo_epi8,
};

/// SSSE3 decoder for packed NCBI2na.
///
/// Each input byte holds four 2-bit bases. The SIMD path decodes 16 input bytes
/// into 64 output bytes by nibble lookup and byte interleaving:
/// high nibble => bases 0,1; low nibble => bases 2,3.
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "ssse3")]
unsafe fn decode_packed_ncbi2na_ssse3(packed: &[u8], len: usize) -> Vec<u8> {
    let mut decoded = vec![0u8; len];
    let full_bytes = len / 4;
    let simd_bytes = full_bytes & !15;

    let first_base_lut = _mm_set_epi8(3, 3, 3, 3, 2, 2, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0);
    let second_base_lut = _mm_set_epi8(3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0);
    let low_nibble_mask = _mm_set1_epi8(0x0f);

    let mut i = 0usize;
    while i < simd_bytes {
        let bytes = _mm_loadu_si128(packed.as_ptr().add(i) as *const __m128i);
        let hi_nibbles = _mm_and_si128(_mm_srli_epi16(bytes, 4), low_nibble_mask);
        let lo_nibbles = _mm_and_si128(bytes, low_nibble_mask);

        let b0 = _mm_shuffle_epi8(first_base_lut, hi_nibbles);
        let b1 = _mm_shuffle_epi8(second_base_lut, hi_nibbles);
        let b2 = _mm_shuffle_epi8(first_base_lut, lo_nibbles);
        let b3 = _mm_shuffle_epi8(second_base_lut, lo_nibbles);

        let p01_lo = _mm_unpacklo_epi8(b0, b1);
        let p01_hi = _mm_unpackhi_epi8(b0, b1);
        let p23_lo = _mm_unpacklo_epi8(b2, b3);
        let p23_hi = _mm_unpackhi_epi8(b2, b3);

        let out0 = _mm_unpacklo_epi16(p01_lo, p23_lo);
        let out1 = _mm_unpackhi_epi16(p01_lo, p23_lo);
        let out2 = _mm_unpacklo_epi16(p01_hi, p23_hi);
        let out3 = _mm_unpackhi_epi16(p01_hi, p23_hi);

        let out_ptr = decoded.as_mut_ptr().add(i * 4);
        _mm_storeu_si128(out_ptr as *mut __m128i, out0);
        _mm_storeu_si128(out_ptr.add(16) as *mut __m128i, out1);
        _mm_storeu_si128(out_ptr.add(32) as *mut __m128i, out2);
        _mm_storeu_si128(out_ptr.add(48) as *mut __m128i, out3);

        i += 16;
    }

    let decoded_bases = simd_bytes * 4;
    if decoded_bases < len {
        let tail = decode_packed_ncbi2na_scalar(&packed[simd_bytes..], len - decoded_bases);
        decoded[decoded_bases..].copy_from_slice(&tail);
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
            round_down: false,
        }
    }

    #[test]
    fn test_decode_packed_ncbi2na_with_ambiguity_overlays_runs() {
        let packed = pack_ncbi2na(&[0, 1, 2, 3, 0, 1, 2, 3]);
        let old_record = (15u32 << 28) | (2u32 << 24) | 2u32;
        let mut old = Vec::new();
        old.extend_from_slice(&1u32.to_be_bytes());
        old.extend_from_slice(&old_record.to_be_bytes());
        let decoded = decode_packed_ncbi2na_with_ambiguity(&packed, 8, &old);
        assert_eq!(&decoded[2..5], &[14, 14, 14]);

        let new_record = (5u32 << 28) | (1u32 << 16);
        let mut new = Vec::new();
        new.extend_from_slice(&(0x8000_0001u32).to_be_bytes());
        new.extend_from_slice(&new_record.to_be_bytes());
        new.extend_from_slice(&6u32.to_be_bytes());
        let decoded = decode_packed_ncbi2na_with_ambiguity(&packed, 8, &new);
        assert_eq!(&decoded[6..8], &[4, 4]);
    }

    #[test]
    fn test_word_hash() {
        // ACGT = 0,1,2,3
        assert_eq!(word_hash_n(&[0, 1, 2, 3], 4), 0b00_01_10_11);
    }

    #[test]
    fn test_lookup_width_selection_matches_ncbi_key_cases() {
        // From NCBI blast_nalookup.c:
        // small blastn word 7 with a small query uses a 6-base lookup.
        assert_eq!(choose_small_na_lut_word(7, 34), 6);
        // Contiguous megablast word 28 uses 11-base lookup below 300k
        // approximate table entries, and 12-base lookup above that.
        assert_eq!(choose_contiguous_mb_lut_word(28, 20_000), 11);
        assert_eq!(choose_contiguous_mb_lut_word(28, 300_000), 12);
    }

    #[test]
    fn test_megablast_prepared_lookup_uses_mb_width_and_stride() {
        let query = vec![0u8; 10_100];
        let prepared = PreparedBlastnQuery::new_megablast(&query, &[], 28);
        assert_eq!(prepared.lookups.len(), 1);
        assert_eq!(prepared.lookups[0].lut_word, 11);
        assert_eq!(prepared.lookups[0].scan_start, 17);
        assert_eq!(prepared.lookups[0].scan_step, 18);
    }

    #[test]
    fn test_decode_packed_ncbi2na_matches_scalar_for_various_lengths() {
        let packed: Vec<u8> = (0..80)
            .map(|i| (i as u8).wrapping_mul(37).wrapping_add(11))
            .collect();

        for len in [
            0usize, 1, 2, 3, 4, 5, 15, 16, 31, 32, 63, 64, 65, 127, 255, 319,
        ] {
            let decoded = decode_packed_ncbi2na(&packed, len);
            let scalar = decode_packed_ncbi2na_scalar(&packed, len);
            assert_eq!(decoded, scalar, "length {}", len);
            assert_eq!(decoded.len(), len);
        }
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
            &query, &rc, &subject, 7, 2, -3, 20, &kbp, 1e6,
            1e10, // very permissive e-value for testing
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
        let results = blastn_ungapped_search(&query, &rc, &subject, 4, 2, -3, 20, &kbp, 1e6, 1e10);
        assert!(
            !results.is_empty(),
            "Should find hit despite mismatches (word_size=4)"
        );
        // The hit may be partial (avoiding the mismatch region)
        assert!(results[0].score > 0);
    }

    #[test]
    fn test_extend_seed_right_stops_when_total_score_goes_negative() {
        let query = vec![0u8, 0, 0, 0];
        let subject = vec![0u8, 1, 0, 0];
        let kbp = test_kbp();

        let hsp = extend_seed(&query, &subject, 0, 0, 2, -3, 20, &kbp, 1e6, 1e10, 0)
            .expect("first matching base should produce an HSP");

        assert_eq!(hsp.score, 2);
        assert_eq!(hsp.query_start, 0);
        assert_eq!(hsp.query_end, 1);
        assert_eq!(hsp.subject_start, 0);
        assert_eq!(hsp.subject_end, 1);
        assert_eq!(hsp.align_length, 1);
    }

    #[test]
    fn test_extend_seed_packed_right_stops_when_total_score_goes_negative() {
        let query = vec![0u8, 0, 0, 0];
        let subject = pack_ncbi2na(&[0, 1, 0, 0]);
        let kbp = test_kbp();

        let hsp = extend_seed_packed(&query, &subject, 4, 0, 0, 2, -3, 20, &kbp, 1e6, 1e10, 0)
            .expect("first matching base should produce an HSP");

        assert_eq!(hsp.score, 2);
        assert_eq!(hsp.query_start, 0);
        assert_eq!(hsp.query_end, 1);
        assert_eq!(hsp.subject_start, 0);
        assert_eq!(hsp.subject_end, 1);
        assert_eq!(hsp.align_length, 1);
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
            &query_plus,
            &query_minus,
            &subject,
            7,
            2,
            -3,
            20,
            &kbp,
            1e6,
            1e10,
        );
        assert!(!results.is_empty(), "Should find hit on minus strand");
        // The hit should be found via the minus strand query (context=1)
        assert!(
            results.iter().any(|h| h.context == 1),
            "Should have a minus-strand hit"
        );
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
        let results =
            blastn_gapped_search(&query, &rc, &subject, 7, 2, -3, 5, 2, 20, &kbp, 1e6, 1e10);
        assert!(!results.is_empty(), "Gapped search should find hit");
        assert_eq!(results[0].gap_opens, 0, "Perfect match should have no gaps");
    }

    #[test]
    fn test_no_match() {
        let query = vec![0u8; 20]; // AAAAAAA...
        let rc = vec![3u8; 20]; // TTTTTTT...
        let subject = vec![1u8; 50]; // CCCCCCC...

        let kbp = test_kbp();
        let results = blastn_ungapped_search(&query, &rc, &subject, 11, 1, -3, 20, &kbp, 1e6, 10.0);
        assert!(results.is_empty(), "Should find no hits");
    }

    // --- Helper: pack a BLASTNA-encoded sequence into NCBI2na packed format ---
    fn pack_ncbi2na(bases: &[u8]) -> Vec<u8> {
        let full_bytes = bases.len() / 4;
        let remainder = bases.len() % 4;
        let total = full_bytes + if remainder > 0 { 1 } else { 0 };
        let mut packed = vec![0u8; total];
        for (i, &b) in bases.iter().enumerate() {
            let byte_idx = i / 4;
            let shift = 6 - 2 * (i % 4);
            packed[byte_idx] |= (b & 3) << shift;
        }
        packed
    }

    // ---- Tests ported from NCBI ntscan_unit_test, nuclwordfinder_unit_test, blastdiag_unit_test ----

    #[test]
    fn test_word_hash_all_bases() {
        // Single-base hashes: A=0, C=1, G=2, T=3
        assert_eq!(word_hash_n(&[0], 1), 0, "A should hash to 0");
        assert_eq!(word_hash_n(&[1], 1), 1, "C should hash to 1");
        assert_eq!(word_hash_n(&[2], 1), 2, "G should hash to 2");
        assert_eq!(word_hash_n(&[3], 1), 3, "T should hash to 3");
    }

    #[test]
    fn test_word_hash_8mer() {
        // ACGTACGT = [0,1,2,3,0,1,2,3]
        // Hash built as: h = (h<<2) | base for each base
        // = 0b_00_01_10_11_00_01_10_11 = 0x1B1B = 6939
        let word: Vec<u8> = vec![0, 1, 2, 3, 0, 1, 2, 3];
        let h = word_hash_n(&word, 8);
        assert_eq!(h, 0b00_01_10_11_00_01_10_11);
        assert_eq!(h, 0x1B1B);
        assert_eq!(h, 6939);
    }

    #[test]
    fn test_scan_finds_all_positions() {
        // Create an 11-mer pattern and place it at subject positions 0, 20, 40
        let pattern: Vec<u8> = vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2]; // 11 bases
        let query = pattern.clone();
        let rc: Vec<u8> = query.iter().rev().map(|&b| 3 - b).collect();

        // Subject of all T's (3) with pattern inserted at 0, 20, 40
        let mut subject = vec![3u8; 60];
        for &offset in &[0usize, 20, 40] {
            for (i, &b) in pattern.iter().enumerate() {
                subject[offset + i] = b;
            }
        }

        let kbp = test_kbp();
        let results = blastn_ungapped_search(&query, &rc, &subject, 11, 2, -3, 20, &kbp, 1e6, 1e10);

        // Should find hits at all 3 positions
        let plus_hits: Vec<&SearchHsp> = results.iter().filter(|h| h.context == 0).collect();
        assert!(
            plus_hits.len() >= 3,
            "Should find hits at all 3 positions, got {}",
            plus_hits.len()
        );
        let starts: Vec<i32> = plus_hits.iter().map(|h| h.subject_start).collect();
        for &pos in &[0i32, 20, 40] {
            assert!(
                starts.contains(&pos),
                "Should find hit at position {}, got starts {:?}",
                pos,
                starts
            );
        }
    }

    #[test]
    fn test_scan_no_hits_random() {
        // Query of all A's (0), subject of all C's (1). No word match possible.
        let query = vec![0u8; 20];
        let rc: Vec<u8> = query.iter().rev().map(|&b| 3 - b).collect(); // all T's
        let subject = vec![1u8; 50]; // all C's

        let kbp = test_kbp();
        let results = blastn_ungapped_search(&query, &rc, &subject, 11, 2, -3, 20, &kbp, 1e6, 1e10);
        assert!(
            results.is_empty(),
            "All-A query vs all-C subject should produce no hits"
        );
    }

    #[test]
    fn test_diagonal_tracking() {
        // Place the same 11-mer at subject positions 5 and 6 (overlapping).
        // The diagonal tracker should prevent redundant extensions but we should
        // still get a hit that covers the overlapping region.
        let pattern: Vec<u8> = vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2]; // 11 bases
        let query = pattern.clone();
        let rc: Vec<u8> = query.iter().rev().map(|&b| 3 - b).collect();

        // Subject: T's with pattern at positions 5 and 6 (overlap => positions 5..17 merge)
        let mut subject = vec![3u8; 30];
        for (i, &b) in pattern.iter().enumerate() {
            subject[5 + i] = b;
        }
        for (i, &b) in pattern.iter().enumerate() {
            subject[6 + i] = b;
        }

        let kbp = test_kbp();
        let results = blastn_ungapped_search(&query, &rc, &subject, 11, 2, -3, 20, &kbp, 1e6, 1e10);

        let plus_hits: Vec<&SearchHsp> = results.iter().filter(|h| h.context == 0).collect();
        assert!(
            !plus_hits.is_empty(),
            "Should find at least one hit covering the overlapping region"
        );
        // The overlapping patterns at 5 and 6 create a merged region 5..17.
        // The diagonal tracker may deduplicate, so we should get exactly one plus-strand hit.
        assert_eq!(
            plus_hits.len(),
            1,
            "Diagonal tracker should merge overlapping seeds into one hit, got {}",
            plus_hits.len()
        );
        // The hit should cover the overlapping region (starts at 5 or 6)
        assert!(
            plus_hits[0].subject_start <= 6,
            "Hit should start at or near position 5"
        );
    }

    #[test]
    fn test_different_word_sizes() {
        // Same query/subject with word_size=7 vs word_size=11.
        // Smaller word size should find at least as many hits.
        let query = vec![0u8, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3]; // 12 bases
        let rc: Vec<u8> = query.iter().rev().map(|&b| 3 - b).collect();

        // Subject with partial match (only first 9 bases match at pos 5)
        let mut subject = vec![3u8; 40];
        subject[5..5 + 9].copy_from_slice(&query[..9]);
        // Also a full match at position 25
        for (i, &b) in query.iter().enumerate() {
            subject[25 + i] = b;
        }

        let kbp = test_kbp();
        let results_7 =
            blastn_ungapped_search(&query, &rc, &subject, 7, 2, -3, 20, &kbp, 1e6, 1e10);
        let results_11 =
            blastn_ungapped_search(&query, &rc, &subject, 11, 2, -3, 20, &kbp, 1e6, 1e10);

        let plus_7: Vec<&SearchHsp> = results_7.iter().filter(|h| h.context == 0).collect();
        let plus_11: Vec<&SearchHsp> = results_11.iter().filter(|h| h.context == 0).collect();
        assert!(
            plus_7.len() >= plus_11.len(),
            "word_size=7 should find >= as many plus-strand hits as word_size=11 ({} vs {})",
            plus_7.len(),
            plus_11.len()
        );
    }

    #[test]
    fn test_gapped_search_with_insertion() {
        // Longer query for robust seed finding
        // Query: ACGTACGTACGTACGTACGTACGT (24 bases, 6 repeats of ACGT)
        let query: Vec<u8> = (0..24).map(|i| [0u8, 1, 2, 3][i % 4]).collect();
        let rc: Vec<u8> = query.iter().rev().map(|&b| 3 - b).collect();
        // Subject: same as query but with 1-base insertion after position 12
        // ACGTACGTACGT A ACGTACGTACGT  (insert A in the middle)
        let mut sub_seq: Vec<u8> = Vec::new();
        sub_seq.extend_from_slice(&query[..12]);
        sub_seq.push(0); // insertion
        sub_seq.extend_from_slice(&query[12..]);
        let mut subject = vec![3u8; 50];
        for (i, &b) in sub_seq.iter().enumerate() {
            subject[5 + i] = b;
        }

        let kbp = test_kbp();
        let results =
            blastn_gapped_search(&query, &rc, &subject, 7, 2, -3, 5, 2, 20, &kbp, 1e6, 1e10);

        assert!(
            !results.is_empty(),
            "Gapped search should find hit with 1-base insertion"
        );
        let best = &results[0];
        // With a 24-base query and 1 insertion, we expect high identity
        assert!(
            best.num_ident >= 20,
            "Should have high identity (got {} ident out of {} align_len)",
            best.num_ident,
            best.align_length
        );
    }

    #[test]
    fn test_gapped_search_with_deletion() {
        // Longer query for robust seed finding
        // Query: ACGTACGTACGTACGTACGTACGT (24 bases)
        let query: Vec<u8> = (0..24).map(|i| [0u8, 1, 2, 3][i % 4]).collect();
        let rc: Vec<u8> = query.iter().rev().map(|&b| 3 - b).collect();
        // Subject: same as query but with 1-base deletion at position 12
        let mut sub_seq: Vec<u8> = Vec::new();
        sub_seq.extend_from_slice(&query[..12]);
        // skip query[12], creating a deletion
        sub_seq.extend_from_slice(&query[13..]);
        let mut subject = vec![3u8; 50];
        for (i, &b) in sub_seq.iter().enumerate() {
            subject[5 + i] = b;
        }

        let kbp = test_kbp();
        let results =
            blastn_gapped_search(&query, &rc, &subject, 7, 2, -3, 5, 2, 20, &kbp, 1e6, 1e10);

        assert!(
            !results.is_empty(),
            "Gapped search should handle 1-base deletion"
        );
        let best = &results[0];
        // With a 24-base query and 1 deletion, we expect high identity
        assert!(
            best.num_ident >= 19,
            "Should have high identity (got {} ident out of {} align_len)",
            best.num_ident,
            best.align_length
        );
    }

    #[test]
    fn test_packed_search_basic() {
        // Verify packed search produces same results as unpacked search for the same input.
        let query: Vec<u8> = vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
        let rc: Vec<u8> = query.iter().rev().map(|&b| 3 - b).collect();
        let mut subject = vec![3u8; 30];
        for (i, &b) in query.iter().enumerate() {
            subject[5 + i] = b;
        }

        let kbp = test_kbp();
        let subject_packed = pack_ncbi2na(&subject);

        let unpacked_results =
            blastn_ungapped_search(&query, &rc, &subject, 7, 2, -3, 20, &kbp, 1e6, 1e10);
        let packed_results = blastn_ungapped_search_packed(
            &query,
            &rc,
            &subject_packed,
            subject.len(),
            7,
            2,
            -3,
            20,
            &kbp,
            1e6,
            1e10,
        );

        assert!(
            !unpacked_results.is_empty(),
            "Unpacked search should find hits"
        );
        assert!(!packed_results.is_empty(), "Packed search should find hits");

        // Both should find the same best hit
        let u = &unpacked_results[0];
        let p = &packed_results[0];
        assert_eq!(
            u.subject_start, p.subject_start,
            "Packed and unpacked should agree on subject_start"
        );
        assert_eq!(
            u.subject_end, p.subject_end,
            "Packed and unpacked should agree on subject_end"
        );
        assert_eq!(
            u.score, p.score,
            "Packed and unpacked should agree on score"
        );
    }

    #[test]
    fn test_search_near_sequence_boundary() {
        // Place matching region at the very start (position 0) and very end of subject.
        let pattern: Vec<u8> = vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2]; // 11 bases
        let query = pattern.clone();
        let rc: Vec<u8> = query.iter().rev().map(|&b| 3 - b).collect();

        // Subject: pattern at start, filler, pattern at end
        let filler_len = 20;
        let total_len = pattern.len() + filler_len + pattern.len();
        let mut subject = vec![3u8; total_len];
        // Place at position 0
        for (i, &b) in pattern.iter().enumerate() {
            subject[i] = b;
        }
        // Place at the very end
        let end_start = total_len - pattern.len();
        for (i, &b) in pattern.iter().enumerate() {
            subject[end_start + i] = b;
        }

        let kbp = test_kbp();
        let results = blastn_ungapped_search(&query, &rc, &subject, 11, 2, -3, 20, &kbp, 1e6, 1e10);

        let plus_hits: Vec<&SearchHsp> = results.iter().filter(|h| h.context == 0).collect();
        assert!(
            plus_hits.len() >= 2,
            "Should find hits at both boundaries, got {}",
            plus_hits.len()
        );
        let starts: Vec<i32> = plus_hits.iter().map(|h| h.subject_start).collect();
        assert!(starts.contains(&0), "Should find hit at position 0");
        assert!(
            starts.contains(&(end_start as i32)),
            "Should find hit at end position {}",
            end_start
        );
    }

    #[test]
    fn test_large_subject_search() {
        // 10,000-base subject with a known query embedded at a specific position.
        let subject_len = 10_000;
        let embed_pos = 4567;
        let query: Vec<u8> = vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2]; // 15 bases
        let rc: Vec<u8> = query.iter().rev().map(|&b| 3 - b).collect();

        // Fill subject with a pseudo-random-ish pattern that avoids matching the query
        // Use alternating C,G (1,2) which won't match ACGTACGT...
        let mut subject: Vec<u8> = (0..subject_len)
            .map(|i| if i % 2 == 0 { 1u8 } else { 2u8 })
            .collect();
        // Embed the query
        for (i, &b) in query.iter().enumerate() {
            subject[embed_pos + i] = b;
        }

        let kbp = test_kbp();
        let results = blastn_ungapped_search(&query, &rc, &subject, 11, 2, -3, 20, &kbp, 1e6, 1e10);

        assert!(
            !results.is_empty(),
            "Should find the embedded query in large subject"
        );
        let plus_hits: Vec<&SearchHsp> = results.iter().filter(|h| h.context == 0).collect();
        assert!(!plus_hits.is_empty(), "Should find plus-strand hit");
        assert_eq!(
            plus_hits[0].subject_start, embed_pos as i32,
            "Hit should be at the embedded position {}",
            embed_pos
        );
        assert_eq!(
            plus_hits[0].num_ident,
            query.len() as i32,
            "Should be a perfect match of {} bases",
            query.len()
        );
    }

    #[test]
    fn test_perfect_ungapped_hsp_fast_path_guard() {
        let perfect = SearchHsp {
            query_start: 0,
            query_end: 22,
            subject_start: 10,
            subject_end: 32,
            score: 22,
            bit_score: 0.0,
            evalue: 0.0,
            num_ident: 22,
            align_length: 22,
            mismatches: 0,
            gap_opens: 0,
            context: 0,
            qseq: Some("ACGTACGTACGTACGTACGTAC".to_string()),
            sseq: Some("ACGTACGTACGTACGTACGTAC".to_string()),
        };
        assert!(is_perfect_ungapped_hsp(&perfect, 1));

        let mut partial = perfect.clone();
        partial.query_start = 1;
        partial.query_end = 21;
        partial.align_length = 20;
        partial.num_ident = 20;
        partial.score = 20;
        assert!(is_perfect_ungapped_hsp(&partial, 1));

        let mut mismatch = perfect.clone();
        mismatch.num_ident = 21;
        mismatch.mismatches = 1;
        assert!(!is_perfect_ungapped_hsp(&mismatch, 1));
    }
}
