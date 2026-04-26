//! Pure Rust BLAST search engine — no FFI.
//! This module implements the complete blastn search pipeline in Rust.

use crate::encoding::NCBI4NA_TO_BLASTNA;
use crate::gapinfo::{GapAlignOpType, GapEditScript};
use crate::greedy::{greedy_align, greedy_align_with_seed, greedy_align_with_seed_packed_subject};
use crate::itree::{Interval, IntervalTree};
use crate::parameters::InitialWordParameters;
use crate::sequence::blastna_to_iupac;
use crate::stat::KarlinBlk;
use crate::traceback::{
    blast_gapped_align, blast_gapped_score_extents_packed_subject, blast_gapped_score_only,
    reevaluate_with_ambiguities_gapped, TracebackResult,
};
use std::time::Instant;

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

enum GappedCandidate {
    Final(SearchHsp),
    Traceback { context: i32, tb: TracebackResult },
}

fn prune_perfect_ungapped_gapped_candidates(
    perfect_seeds: &mut Vec<SearchHsp>,
    min_diag_separation: i32,
) {
    if perfect_seeds.len() <= 1 {
        return;
    }
    purge_common_endpoint_hsps(perfect_seeds);
    dedup_hsps_with_min_diag_separation(perfect_seeds, min_diag_separation);
}

struct NaLookup<'a> {
    context: i32,
    query: &'a [u8],
    query_nomask: &'a [u8],
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

#[derive(Clone, Copy)]
struct PreliminarySeed {
    query: usize,
    subject: usize,
}

#[derive(Clone, Copy)]
struct PreliminaryGappedResult {
    score: i32,
    prelim_q_start: usize,
    prelim_q_end: usize,
    prelim_s_start: usize,
    prelim_s_end: usize,
    gapped_start_q: usize,
    gapped_start_s: usize,
    q_length: usize,
    s_length: usize,
    private_q_start: usize,
    private_s_start: usize,
    score_left: i32,
    score_right: i32,
}

#[derive(Default)]
struct PackedGappedStageStats {
    filtered_seed_count: usize,
    preliminary_pass_count: usize,
    preliminary_score_ms: u128,
    traceback_ms: u128,
    traceback_count: usize,
}

#[derive(Clone, Copy)]
struct PackedTracebackSeed {
    raw_query: usize,
    raw_subject: usize,
    traceback_query: usize,
    traceback_subject: usize,
}

pub struct PackedDiagScratch {
    pub last_hit: Vec<i32>,
    pub flag: Vec<i32>,
    pub hit_len: Vec<i32>,
}

impl PackedDiagScratch {
    fn new(len: usize) -> Self {
        Self {
            last_hit: vec![0; len],
            flag: vec![0; len],
            hit_len: vec![0; len],
        }
    }

    fn resize_and_clear(&mut self, len: usize) {
        if self.last_hit.len() != len {
            self.last_hit.resize(len, 0);
            self.flag.resize(len, 0);
            self.hit_len.resize(len, 0);
        } else {
            self.last_hit.fill(0);
            self.flag.fill(0);
            self.hit_len.fill(0);
        }
    }
}

#[derive(Clone, Copy)]
struct ExactWordHit {
    query_start: usize,
    subject_start: usize,
    subject_match_end: usize,
}

#[inline]
fn find_exact_word_hit_packed_aligned(
    query: &[u8],
    subject_packed: &[u8],
    subject_len: usize,
    word_size: usize,
    lut_word: usize,
    qp: usize,
    sp: usize,
) -> Option<ExactWordHit> {
    let ext_to = word_size.saturating_sub(lut_word);
    let mut ext_left = 0usize;
    let ext_max = ext_to.min(sp).min(qp);

    while ext_left < ext_max {
        let q_base = qp - ext_left;
        let s_base = sp - ext_left;
        let byte_index = (s_base / 4).checked_sub(1)?;
        let byte = subject_packed[byte_index];

        if (byte & 3) != query[q_base - 1] {
            break;
        }
        ext_left += 1;
        if ext_left == ext_max {
            break;
        }
        if ((byte >> 2) & 3) != query[q_base - 2] {
            break;
        }
        ext_left += 1;
        if ext_left == ext_max {
            break;
        }
        if ((byte >> 4) & 3) != query[q_base - 3] {
            break;
        }
        ext_left += 1;
        if ext_left == ext_max {
            break;
        }
        if (byte >> 6) != query[q_base - 4] {
            break;
        }
        ext_left += 1;
    }

    if ext_left < ext_to {
        let mut ext_right = 0usize;
        let ext_max = ext_to - ext_left;
        if sp + lut_word + ext_max > subject_len || qp + lut_word + ext_max > query.len() {
            return None;
        }

        while ext_right < ext_max {
            let q_base = qp + lut_word + ext_right;
            let s_base = sp + lut_word + ext_right;
            let byte = subject_packed[s_base / 4];

            if (byte >> 6) != query[q_base] {
                break;
            }
            ext_right += 1;
            if ext_right == ext_max {
                break;
            }
            if ((byte >> 4) & 3) != query[q_base + 1] {
                break;
            }
            ext_right += 1;
            if ext_right == ext_max {
                break;
            }
            if ((byte >> 2) & 3) != query[q_base + 2] {
                break;
            }
            ext_right += 1;
            if ext_right == ext_max {
                break;
            }
            if (byte & 3) != query[q_base + 3] {
                break;
            }
            ext_right += 1;
        }

        if ext_left + ext_right < ext_to {
            return None;
        }
    }

    let query_start = qp - ext_left;
    let subject_start = sp - ext_left;
    Some(ExactWordHit {
        query_start,
        subject_start,
        subject_match_end: subject_start + word_size,
    })
}

fn traceback_window_padding(x_dropoff: i32, gap_open: i32, gap_extend: i32) -> usize {
    let gap_oe = gap_open.saturating_add(gap_extend).max(1) as usize;
    let x = x_dropoff.max(gap_oe as i32) as usize;
    // Heuristic tighter seed window for traceback. NCBI does not hand ALIGN_EX
    // the full remaining subject tail on long near-self hits; keeping a few
    // hundred bases of slack around the ungapped HSP is enough to preserve
    // local refinement while preventing pathological 100k+ DP widths.
    ((x / gap_oe) + 1) * 32 + 128
}

fn greedy_traceback_alignment(
    query: &[u8],
    subject: &[u8],
    seed_q: usize,
    seed_s: usize,
    reward: i32,
    penalty: i32,
    x_dropoff: i32,
) -> Option<TracebackResult> {
    let mut adjusted_seed_s = seed_s;
    let mut adjusted_subject_len = subject.len();
    let start_shift = adjust_subject_range(
        &mut adjusted_seed_s,
        &mut adjusted_subject_len,
        seed_q + 1,
        query.len(),
    );
    let adjusted_subject = &subject[start_shift..start_shift + adjusted_subject_len];
    let (score, query_start, query_end, subject_start, subject_end, edit_script) = greedy_align(
        query,
        adjusted_subject,
        seed_q,
        adjusted_seed_s,
        reward,
        penalty,
        x_dropoff,
    )?;
    Some(TracebackResult {
        score,
        edit_script,
        query_start,
        query_end,
        subject_start: subject_start + start_shift,
        subject_end: subject_end + start_shift,
    })
}

fn affine_traceback_alignment_windowed(
    query: &[u8],
    subject: &[u8],
    seed: &SearchHsp,
    seed_q: usize,
    seed_s: usize,
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    x_dropoff: i32,
) -> Option<TracebackResult> {
    let padding = traceback_window_padding(x_dropoff, gap_open, gap_extend);
    let q_start = seed.query_start.max(0) as usize;
    let q_end = seed.query_end.max(seed.query_start).max(0) as usize;
    let s_start = seed.subject_start.max(0) as usize;
    let s_end = seed.subject_end.max(seed.subject_start).max(0) as usize;

    let (q_lo, q_hi, s_lo, s_hi) = if std::env::var_os("BLAST_RS_FULL_TB_WINDOW").is_some() {
        (0, query.len(), 0, subject.len())
    } else {
        (
            q_start.saturating_sub(padding),
            query.len().min(q_end.saturating_add(padding)),
            s_start.saturating_sub(padding),
            subject.len().min(s_end.saturating_add(padding)),
        )
    };

    if q_lo >= q_hi || s_lo >= s_hi {
        return None;
    }

    let local_seed_q = seed_q.saturating_sub(q_lo).min(q_hi - q_lo - 1);
    let local_seed_s = seed_s.saturating_sub(s_lo).min(s_hi - s_lo - 1);
    let mut tb = blast_gapped_align(
        &query[q_lo..q_hi],
        &subject[s_lo..s_hi],
        local_seed_q,
        local_seed_s,
        reward,
        penalty,
        gap_open,
        gap_extend,
        x_dropoff,
    )?;
    tb.query_start += q_lo;
    tb.query_end += q_lo;
    tb.subject_start += s_lo;
    tb.subject_end += s_lo;
    Some(tb)
}

fn blast_gapped_align_windowed(
    query: &[u8],
    subject: &[u8],
    seed: &SearchHsp,
    seed_q: usize,
    seed_s: usize,
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    x_dropoff: i32,
) -> Option<TracebackResult> {
    if gap_open == 0 && gap_extend == 0 {
        return greedy_traceback_alignment(query, subject, seed_q, seed_s, reward, penalty, x_dropoff);
    }

    affine_traceback_alignment_windowed(
        query, subject, seed, seed_q, seed_s, reward, penalty, gap_open, gap_extend, x_dropoff,
    )
}

fn elapsed_ms_since(start: Option<Instant>) -> u128 {
    start.map(|t| t.elapsed().as_millis()).unwrap_or(0)
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
        Self::new_with_nomask(
            query_plus,
            query_minus,
            query_plus,
            query_minus,
            word_size,
        )
    }

    pub fn new_megablast(query_plus: &'a [u8], query_minus: &'a [u8], word_size: usize) -> Self {
        Self::new_megablast_with_nomask(query_plus, query_minus, query_plus, query_minus, word_size)
    }

    pub fn new_megablast_with_nomask(
        query_plus: &'a [u8],
        query_minus: &'a [u8],
        query_plus_nomask: &'a [u8],
        query_minus_nomask: &'a [u8],
        word_size: usize,
    ) -> Self {
        let mut lookups = Vec::with_capacity(2);
        if let Some(lookup) = NaLookup::new_contiguous(0, query_plus, query_plus_nomask, word_size)
        {
            lookups.push(lookup);
        }
        if let Some(lookup) =
            NaLookup::new_contiguous(1, query_minus, query_minus_nomask, word_size)
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

    pub fn new_with_nomask(
        query_plus: &'a [u8],
        query_minus: &'a [u8],
        query_plus_nomask: &'a [u8],
        query_minus_nomask: &'a [u8],
        word_size: usize,
    ) -> Self {
        Self::new_with_selector(
            query_plus,
            query_minus,
            query_plus_nomask,
            query_minus_nomask,
            word_size,
            choose_small_na_lut_word,
        )
    }

    fn new_with_selector(
        query_plus: &'a [u8],
        query_minus: &'a [u8],
        query_plus_nomask: &'a [u8],
        query_minus_nomask: &'a [u8],
        word_size: usize,
        choose_lut_word: fn(usize, usize) -> usize,
    ) -> Self {
        let mut lookups = Vec::with_capacity(2);
        let approx_entries = query_plus.len().saturating_sub(word_size).saturating_add(1)
            + query_minus
                .len()
                .saturating_sub(word_size)
                .saturating_add(1);
        if let Some(lookup) = NaLookup::new(
            0,
            query_plus,
            query_plus_nomask,
            word_size,
            approx_entries,
            choose_lut_word,
        )
        {
            lookups.push(lookup);
        }
        if let Some(lookup) = NaLookup::new(
            1,
            query_minus,
            query_minus_nomask,
            word_size,
            approx_entries,
            choose_lut_word,
        )
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

    pub fn last_hit_scratch(&self) -> Vec<PackedDiagScratch> {
        self.lookups
            .iter()
            .map(|lookup| PackedDiagScratch::new(lookup.diag_array_len))
            .collect()
    }

    pub fn decoded_last_hit_scratch(&self) -> Vec<Vec<i32>> {
        self.lookups
            .iter()
            .map(|lookup| vec![0; lookup.diag_array_len])
            .collect()
    }
}

impl<'a> NaLookup<'a> {
    fn new_contiguous(
        context: i32,
        query: &'a [u8],
        query_nomask: &'a [u8],
        word_size: usize,
    ) -> Option<Self> {
        let (approx_entries, max_q_off) = approximate_lookup_segments(query);
        let choose_lut_word = if should_use_small_na_lookup(word_size, approx_entries, max_q_off) {
            choose_small_na_lut_word
        } else {
            choose_contiguous_mb_lut_word
        };
        Self::new(
            context,
            query,
            query_nomask,
            word_size,
            approx_entries,
            choose_lut_word,
        )
    }

    fn new(
        context: i32,
        query: &'a [u8],
        query_nomask: &'a [u8],
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
        let scan_start = 0usize;
        let scan_step = if std::env::var_os("BLAST_RS_FORCE_SCAN_STEP1").is_some() {
            1
        } else {
            word_size - lut_word + 1
        };

        let mut lut: Vec<i32> = vec![-1; lut_size];
        let mut next: Vec<i32> = vec![-1; query.len()];
        let pv_size = lut_size.div_ceil(64);
        let mut pv: Vec<u64> = vec![0; pv_size];

        // Match NCBI's lookup admission behavior: indexing depends on the
        // lookup word width, while the remaining full-word verification happens
        // later during exact-hit extension.
        let lookup_mask_span = lut_word;
        let unmasked_run_ends = eligible_lookup_run_ends(query, lookup_mask_span);

        for i in (0..=(query.len() - lut_word)).rev() {
            if unmasked_run_ends[i] == 0 || i + lookup_mask_span > unmasked_run_ends[i] {
                continue;
            }
            let key = word_hash_n(&query[i..i + lut_word], lut_word) as usize;
            next[i] = lut[key];
            lut[key] = i as i32;
            pv[key >> 6] |= 1u64 << (key & 63);
        }

        let diag_array_len = (query.len() * 2).next_power_of_two().max(256);
        let diag_mask = diag_array_len - 1;

        if std::env::var_os("BLAST_RS_TRACE_ROW3_LOOKUP").is_some() {
            for &probe in &[873usize, 1018, 1025, 1048, 1116, 1228, 1253, 1294] {
                if probe + lut_word > query.len() {
                    continue;
                }
                let unmasked_end = unmasked_run_ends[probe];
                let full_word_ok = unmasked_end != 0 && probe + lookup_mask_span <= unmasked_end;
                let key = word_hash_n(&query[probe..probe + lut_word], lut_word) as usize;
                let mut indexed = false;
                let mut pos = lut[key];
                while pos >= 0 {
                    if pos as usize == probe {
                        indexed = true;
                        break;
                    }
                    pos = next[pos as usize];
                }
                eprintln!(
                    "[row3-lookup] ctx={} probe={} indexed={} full_word_ok={} unmasked_end={} key={}",
                    context, probe, indexed, full_word_ok, unmasked_end, key
                );
                if matches!(probe, 976 | 1048 | 1163 | 1253) {
                    for &sp in &[3_897_222usize, 3_897_726, 3_897_411, 3_897_899] {
                        eprintln!(
                            "[row3-lookup-key] ctx={} probe={} key={} subject_probe={}",
                            context, probe, key, sp
                        );
                    }
                }
            }
        }
        if let Some(probes) = trace_lookup_probes() {
            for probe in probes {
                if probe + lut_word > query.len() {
                    continue;
                }
                let unmasked_end = unmasked_run_ends[probe];
                let full_word_ok = unmasked_end != 0 && probe + lookup_mask_span <= unmasked_end;
                let key = word_hash_n(&query[probe..probe + lut_word], lut_word) as usize;
                let mut indexed = false;
                let mut pos = lut[key];
                while pos >= 0 {
                    if pos as usize == probe {
                        indexed = true;
                        break;
                    }
                    pos = next[pos as usize];
                }
                eprintln!(
                    "[trace-lookup] ctx={} probe={} indexed={} full_word_ok={} unmasked_end={} key={}",
                    context, probe, indexed, full_word_ok, unmasked_end, key
                );
            }
        }

        Some(NaLookup {
            context,
            query,
            query_nomask,
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

fn approximate_lookup_segments(query: &[u8]) -> (usize, usize) {
    let mut approx_entries = 0usize;
    let mut max_q_off = 0usize;
    let mut saw_unmasked = false;
    for (idx, &base) in query.iter().enumerate() {
        if base < 4 {
            approx_entries += 1;
            max_q_off = idx;
            saw_unmasked = true;
        }
    }
    if !saw_unmasked {
        max_q_off = 0;
    }
    (approx_entries, max_q_off)
}

fn should_use_small_na_lookup(word_size: usize, approx_entries: usize, max_q_off: usize) -> bool {
    let small_candidate = match word_size {
        4..=8 => true,
        9 => approx_entries < 21_000,
        10 => approx_entries < 8_500,
        11 => approx_entries < 12_000,
        12 => approx_entries < 8_500,
        _ => approx_entries < 8_500,
    };
    small_candidate && approx_entries < 32_767 && max_q_off < 32_768
}

fn eligible_lookup_run_ends(query: &[u8], word_size: usize) -> Vec<usize> {
    let mut run_ends = vec![0usize; query.len()];
    let mut pos = 0usize;
    while pos < query.len() {
        if query[pos] >= 4 {
            pos += 1;
            continue;
        }
        let start = pos;
        while pos < query.len() && query[pos] < 4 {
            pos += 1;
        }
        let end = pos;
        if end - start >= word_size {
            for slot in &mut run_ends[start..end] {
                *slot = end;
            }
        }
    }
    run_ends
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

#[allow(clippy::too_many_arguments)]
fn small_na_extend_packed_hits(
    query: &[u8],
    subject_packed: &[u8],
    subject_len: usize,
    word_size: usize,
    lut_word: usize,
    lut: &[i32],
    next: &[i32],
    last_hit: &mut [i32],
    flag: &mut [i32],
    hit_len: &mut [i32],
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
    nucl_score_table: &[i32; 256],
    reduced_nucl_cutoff_score: i32,
    hsps: &mut Vec<SearchHsp>,
    exact_word_hit_fn: fn(&[u8], &[u8], usize, usize, usize, usize, usize) -> Option<ExactWordHit>,
) {
    let sp = base_pos;
    let mut q_pos = lut[h];
    while q_pos >= 0 {
        let qp = q_pos as usize;
        let trace_seed_pair = matches_trace_seed_pair(qp, sp);
        let trace_subject_pos = trace_scan_positions()
            .map(|positions| positions.contains(&sp))
            .unwrap_or(false);
        if let Some((_, _, ss, se)) = trace_target_bounds() {
            if (sp as i32) <= se && (sp + word_size) as i32 >= ss {
                eprintln!(
                    "[trace-target] raw-hit ctx={} qp={} sp={} word_size={} lut_word={} hash={}",
                    context, qp, sp, word_size, lut_word, h
                );
            }
        }
        if trace_seed_pair {
            eprintln!(
                "[trace-seed] raw-hit ctx={} qp={} sp={} word_size={} lut_word={} hash={}",
                context, qp, sp, word_size, lut_word, h
            );
        }
        if trace_subject_pos {
            eprintln!(
                "[trace-hitlist] raw-hit ctx={} qp={} sp={} word_size={} lut_word={} hash={}",
                context, qp, sp, word_size, lut_word, h
            );
        }
        let hit = exact_word_hit_fn(query, subject_packed, subject_len, word_size, lut_word, qp, sp);
        if let Some(hit) = hit {
            if trace_seed_pair {
                eprintln!(
                    "[trace-seed] exact-hit ctx={} raw=({}, {}) hit_q={}..{} hit_s={}..{}",
                    context,
                    qp,
                    sp,
                    hit.query_start,
                    hit.query_start + word_size,
                    hit.subject_start,
                    hit.subject_start + word_size
                );
            }
            if trace_subject_pos {
                eprintln!(
                    "[trace-hitlist] exact-hit ctx={} raw=({}, {}) hit_q={}..{} hit_s={}..{}",
                    context,
                    qp,
                    sp,
                    hit.query_start,
                    hit.query_start + word_size,
                    hit.subject_start,
                    hit.subject_start + word_size
                );
            }
            let use_hash_initial_hit = false;
            if use_hash_initial_hit {
                diag_hash_extend_initial_hit_packed(
                    query,
                    subject_packed,
                    subject_len,
                    word_size,
                    last_hit,
                    flag,
                    hit_len,
                    diag_mask,
                    reward,
                    penalty,
                    x_dropoff,
                    kbp,
                    search_space,
                    evalue_threshold,
                    context,
                    nucl_score_table,
                    reduced_nucl_cutoff_score,
                    hit,
                    hsps,
                );
            } else {
                blast_na_extend_direct_packed(
                    query,
                    subject_packed,
                    subject_len,
                    word_size,
                    last_hit,
                    flag,
                    hit_len,
                    diag_mask,
                    reward,
                    penalty,
                    x_dropoff,
                    kbp,
                    search_space,
                    evalue_threshold,
                    context,
                    nucl_score_table,
                    reduced_nucl_cutoff_score,
                    hit,
                    hsps,
                );
            }
        } else if let Some((qs, qe, ss, se)) = trace_target_bounds() {
            if (qp as i32) <= qe
                && (qp + lut_word) as i32 >= qs
                && (sp as i32) <= se
                && (sp + lut_word) as i32 >= ss
            {
                let (ext_left, ext_right, ext_to) = exact_word_hit_packed_naive_extents(
                    query,
                    subject_packed,
                    subject_len,
                    word_size,
                    lut_word,
                    qp,
                    sp,
                );
                eprintln!(
                    "[trace-target] exact-hit-miss ctx={} raw_q={} raw_s={} word_size={} lut_word={} ext_left={} ext_right={} ext_to={}",
                    context, qp, sp, word_size, lut_word, ext_left, ext_right, ext_to
                );
            }
        }
        if trace_seed_pair && hit.is_none() {
            let (ext_left, ext_right, ext_to) = exact_word_hit_packed_naive_extents(
                query,
                subject_packed,
                subject_len,
                word_size,
                lut_word,
                qp,
                sp,
            );
            eprintln!(
                "[trace-seed] exact-hit-miss ctx={} raw=({}, {}) ext_left={} ext_right={} ext_to={}",
                context, qp, sp, ext_left, ext_right, ext_to
            );
        }
        if trace_subject_pos && hit.is_none() {
            let (ext_left, ext_right, ext_to) = exact_word_hit_packed_naive_extents(
                query,
                subject_packed,
                subject_len,
                word_size,
                lut_word,
                qp,
                sp,
            );
            eprintln!(
                "[trace-hitlist] exact-hit-miss ctx={} raw=({}, {}) ext_left={} ext_right={} ext_to={}",
                context, qp, sp, ext_left, ext_right, ext_to
            );
        }
        q_pos = next[qp];
    }
}

#[allow(clippy::too_many_arguments)]
fn blast_na_extend_packed_generic(
    query: &[u8],
    subject_packed: &[u8],
    subject_len: usize,
    word_size: usize,
    lut_word: usize,
    lut: &[i32],
    next: &[i32],
    last_hit: &mut [i32],
    flag: &mut [i32],
    hit_len: &mut [i32],
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
    nucl_score_table: &[i32; 256],
    reduced_nucl_cutoff_score: i32,
    hsps: &mut Vec<SearchHsp>,
) {
    small_na_extend_packed_hits(
        query,
        subject_packed,
        subject_len,
        word_size,
        lut_word,
        lut,
        next,
        last_hit,
        flag,
        hit_len,
        diag_mask,
        h,
        base_pos,
        reward,
        penalty,
        x_dropoff,
        kbp,
        search_space,
        evalue_threshold,
        context,
        nucl_score_table,
        reduced_nucl_cutoff_score,
        hsps,
        find_exact_word_hit_packed,
    );
}

#[allow(clippy::too_many_arguments)]
fn blast_na_extend_packed_aligned(
    query: &[u8],
    subject_packed: &[u8],
    subject_len: usize,
    word_size: usize,
    lut_word: usize,
    lut: &[i32],
    next: &[i32],
    last_hit: &mut [i32],
    flag: &mut [i32],
    hit_len: &mut [i32],
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
    nucl_score_table: &[i32; 256],
    reduced_nucl_cutoff_score: i32,
    hsps: &mut Vec<SearchHsp>,
) {
    small_na_extend_packed_hits(
        query,
        subject_packed,
        subject_len,
        word_size,
        lut_word,
        lut,
        next,
        last_hit,
        flag,
        hit_len,
        diag_mask,
        h,
        base_pos,
        reward,
        penalty,
        x_dropoff,
        kbp,
        search_space,
        evalue_threshold,
        context,
        nucl_score_table,
        reduced_nucl_cutoff_score,
        hsps,
        find_exact_word_hit_packed_aligned,
    );
}

/// Port boundary for NCBI's aligned small-query nucleotide extension callback
/// (`s_BlastSmallNaExtendAlignedOneByte` / aligned cases).
#[inline]
fn small_na_extend_packed_aligned_one_byte(
    query: &[u8],
    subject_packed: &[u8],
    subject_len: usize,
    word_size: usize,
    lut_word: usize,
    lut: &[i32],
    next: &[i32],
    last_hit: &mut [i32],
    flag: &mut [i32],
    hit_len: &mut [i32],
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
    nucl_score_table: &[i32; 256],
    reduced_nucl_cutoff_score: i32,
    hsps: &mut Vec<SearchHsp>,
) {
    blast_na_extend_packed_aligned(
        query,
        subject_packed,
        subject_len,
        word_size,
        lut_word,
        lut,
        next,
        last_hit,
        flag,
        hit_len,
        diag_mask,
        h,
        base_pos,
        reward,
        penalty,
        x_dropoff,
        kbp,
        search_space,
        evalue_threshold,
        context,
        nucl_score_table,
        reduced_nucl_cutoff_score,
        hsps,
    );
}

/// Port boundary for NCBI's generic small-query nucleotide extension callback
/// (`s_BlastSmallNaExtend`).
#[inline]
fn small_na_extend_packed_generic(
    query: &[u8],
    subject_packed: &[u8],
    subject_len: usize,
    word_size: usize,
    lut_word: usize,
    lut: &[i32],
    next: &[i32],
    last_hit: &mut [i32],
    flag: &mut [i32],
    hit_len: &mut [i32],
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
    nucl_score_table: &[i32; 256],
    reduced_nucl_cutoff_score: i32,
    hsps: &mut Vec<SearchHsp>,
) {
    blast_na_extend_packed_generic(
        query,
        subject_packed,
        subject_len,
        word_size,
        lut_word,
        lut,
        next,
        last_hit,
        flag,
        hit_len,
        diag_mask,
        h,
        base_pos,
        reward,
        penalty,
        x_dropoff,
        kbp,
        search_space,
        evalue_threshold,
        context,
        nucl_score_table,
        reduced_nucl_cutoff_score,
        hsps,
    );
}
#[inline]
fn small_na_extend_packed(
    query: &[u8],
    subject_packed: &[u8],
    subject_len: usize,
    word_size: usize,
    lut_word: usize,
    lut: &[i32],
    next: &[i32],
    last_hit: &mut [i32],
    flag: &mut [i32],
    hit_len: &mut [i32],
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
    nucl_score_table: &[i32; 256],
    reduced_nucl_cutoff_score: i32,
    hsps: &mut Vec<SearchHsp>,
) {
    let use_aligned_callback = (base_pos % 4) == 0 && (word_size % 4) == 0 && (lut_word % 4) == 0;
    if use_aligned_callback {
        small_na_extend_packed_aligned_one_byte(
            query,
            subject_packed,
            subject_len,
            word_size,
            lut_word,
            lut,
            next,
            last_hit,
            flag,
            hit_len,
            diag_mask,
            h,
            base_pos,
            reward,
            penalty,
            x_dropoff,
            kbp,
            search_space,
            evalue_threshold,
            context,
            nucl_score_table,
            reduced_nucl_cutoff_score,
            hsps,
        );
    } else {
        small_na_extend_packed_generic(
            query,
            subject_packed,
            subject_len,
            word_size,
            lut_word,
            lut,
            next,
            last_hit,
            flag,
            hit_len,
            diag_mask,
            h,
            base_pos,
            reward,
            penalty,
            x_dropoff,
            kbp,
            search_space,
            evalue_threshold,
            context,
            nucl_score_table,
            reduced_nucl_cutoff_score,
            hsps,
        );
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

#[inline]
fn find_exact_word_hit_packed(
    query: &[u8],
    subject_packed: &[u8],
    subject_len: usize,
    word_size: usize,
    lut_word: usize,
    qp: usize,
    sp: usize,
) -> Option<ExactWordHit> {
    if std::env::var_os("BLAST_RS_EXACT_WORD_NAIVE").is_some() {
        return find_exact_word_hit_packed_naive(
            query,
            subject_packed,
            subject_len,
            word_size,
            lut_word,
            qp,
            sp,
        );
    }

    if (qp % 4) == 0 && (sp % 4) == 0 && (word_size % 4) == 0 && (lut_word % 4) == 0 {
        return find_exact_word_hit_packed_aligned(
            query,
            subject_packed,
            subject_len,
            word_size,
            lut_word,
            qp,
            sp,
        );
    }

    let mut q_offset = qp;
    let mut s_offset = sp;
    let mut ext_left = 0usize;
    let mut ext_right = 0usize;
    let mut ext_max = word_size.saturating_sub(lut_word).min(s_offset).min(q_offset);

    // Match C s_BlastSmallNaExtend: start from the first 4-base boundary to the
    // right of the hit and count exact matches in 4-base groups.
    let rsdl = 4usize - (s_offset % 4);
    q_offset += rsdl;
    s_offset += rsdl;
    ext_max += rsdl;

    let mut q_off = q_offset;
    let mut s_off = s_offset;
    while ext_left < ext_max {
        let remaining = (ext_max - ext_left).min(4);
        let mut bases = 0usize;
        for i in 1..=remaining {
            if packed_base_at(subject_packed, s_off - i) != query[q_off - i] {
                break;
            }
            bases += 1;
        }
        ext_left += bases;
        if bases < 4 {
            break;
        }
        q_off -= 4;
        s_off -= 4;
    }
    ext_left = ext_left.min(ext_max);

    q_off = q_offset;
    s_off = s_offset;
    ext_max = (word_size - ext_left)
        .min(subject_len.saturating_sub(s_off))
        .min(query.len().saturating_sub(q_off));
    while ext_right < ext_max {
        let remaining = (ext_max - ext_right).min(4);
        let mut bases = 0usize;
        for i in 0..remaining {
            if packed_base_at(subject_packed, s_off + i) != query[q_off + i] {
                break;
            }
            bases += 1;
        }
        ext_right += bases;
        if bases < 4 {
            break;
        }
        q_off += 4;
        s_off += 4;
    }
    ext_right = ext_right.min(ext_max);

    if ext_left + ext_right < word_size {
        return None;
    }

    let query_start = q_offset - ext_left;
    let subject_start = s_offset - ext_left;
    Some(ExactWordHit {
        query_start,
        subject_start,
        subject_match_end: subject_start + word_size,
    })
}

#[inline]
fn find_exact_word_hit_packed_naive(
    query: &[u8],
    subject_packed: &[u8],
    subject_len: usize,
    word_size: usize,
    lut_word: usize,
    qp: usize,
    sp: usize,
) -> Option<ExactWordHit> {
    let ext_to = word_size.saturating_sub(lut_word);
    let mut ext_left = 0usize;
    while ext_left < ext_to && ext_left < qp && ext_left < sp {
        let q_idx = qp - ext_left - 1;
        let s_idx = sp - ext_left - 1;
        if packed_base_at(subject_packed, s_idx) != query[q_idx] {
            break;
        }
        ext_left += 1;
    }

    let mut ext_right = 0usize;
    while ext_left + ext_right < ext_to {
        let q_idx = qp + lut_word + ext_right;
        let s_idx = sp + lut_word + ext_right;
        if q_idx >= query.len() || s_idx >= subject_len {
            return None;
        }
        if packed_base_at(subject_packed, s_idx) != query[q_idx] {
            break;
        }
        ext_right += 1;
    }

    if ext_left + ext_right < ext_to {
        return None;
    }

    let query_start = qp - ext_left;
    let subject_start = sp - ext_left;
    Some(ExactWordHit {
        query_start,
        subject_start,
        subject_match_end: subject_start + word_size,
    })
}

#[inline]
fn exact_word_hit_packed_naive_extents(
    query: &[u8],
    subject_packed: &[u8],
    subject_len: usize,
    word_size: usize,
    lut_word: usize,
    qp: usize,
    sp: usize,
) -> (usize, usize, usize) {
    let ext_to = word_size.saturating_sub(lut_word);
    let mut ext_left = 0usize;
    while ext_left < ext_to && ext_left < qp && ext_left < sp {
        let q_idx = qp - ext_left - 1;
        let s_idx = sp - ext_left - 1;
        if packed_base_at(subject_packed, s_idx) != query[q_idx] {
            break;
        }
        ext_left += 1;
    }

    let mut ext_right = 0usize;
    while ext_left + ext_right < ext_to {
        let q_idx = qp + lut_word + ext_right;
        let s_idx = sp + lut_word + ext_right;
        if q_idx >= query.len() || s_idx >= subject_len {
            break;
        }
        if packed_base_at(subject_packed, s_idx) != query[q_idx] {
            break;
        }
        ext_right += 1;
    }

    (ext_left, ext_right, ext_to)
}

#[allow(clippy::too_many_arguments)]
fn diag_initial_hit_core_packed(
    query: &[u8],
    subject_packed: &[u8],
    subject_len: usize,
    word_size: usize,
    last_hit: &mut [i32],
    flag: &mut [i32],
    hit_len: &mut [i32],
    diag_mask: usize,
    reward: i32,
    penalty: i32,
    x_dropoff: i32,
    kbp: &KarlinBlk,
    search_space: f64,
    evalue_threshold: f64,
    context: i32,
    nucl_score_table: &[i32; 256],
    reduced_nucl_cutoff_score: i32,
    hit: ExactWordHit,
    hsps: &mut Vec<SearchHsp>,
) {
    if matches_trace_seed_pair(hit.query_start, hit.subject_start) {
        eprintln!(
            "[trace-seed] diag-entry ctx={} q={}..{} s={}..{} last_hit={} diag={}",
            context,
            hit.query_start,
            hit.query_start + word_size,
            hit.subject_start,
            hit.subject_start + word_size,
            last_hit[(hit.subject_start + (diag_mask + 1) - hit.query_start) & diag_mask],
            (hit.subject_start + (diag_mask + 1) - hit.query_start) & diag_mask
        );
    }
    let window_size = 0i32;
    let two_hits = window_size > 0;
    let real_diag = (hit.subject_start + (diag_mask + 1) - hit.query_start) & diag_mask;
    let hit_saved = flag[real_diag];
    let mut hit_ready = true;
    let mut s_end_pos = hit.subject_match_end as i32;
    let word_type = 1i32;
    if let Some((qs, qe, ss, se)) = trace_target_bounds() {
        if (hit.query_start as i32) <= qe
            && (hit.query_start + word_size) as i32 >= qs
            && (hit.subject_start as i32) <= se
            && (hit.subject_start + word_size) as i32 >= ss
        {
            eprintln!(
                "[trace-target] exact-hit ctx={} q={}..{} s={}..{} last_hit={} diag={}",
                context,
                hit.query_start,
                hit.query_start + word_size,
                hit.subject_start,
                hit.subject_start + word_size,
                last_hit[real_diag],
                real_diag
            );
        }
    }
    if (last_hit[real_diag] as usize) > hit.subject_start {
        if let Some((qs, qe, ss, se)) = trace_target_bounds() {
            if (hit.query_start as i32) <= qe
                && (hit.query_start + word_size) as i32 >= qs
                && (hit.subject_start as i32) <= se
                && (hit.subject_start + word_size) as i32 >= ss
            {
                eprintln!(
                    "[trace-target] exact-hit-suppressed ctx={} q={} s={} last_hit={}",
                    context, hit.query_start, hit.subject_start, last_hit[real_diag]
                );
            }
        }
        return;
    }

    if !two_hits || hit_saved != 0 || word_type != 1 {
        let extended = extend_seed_packed(
            query,
            subject_packed,
            subject_len,
            hit.query_start,
            hit.subject_start,
            hit.subject_match_end,
            reward,
            penalty,
            x_dropoff,
            word_size,
            nucl_score_table,
            reduced_nucl_cutoff_score,
        );

        if let Some(data) = extended {
            let Some(hsp) = build_packed_hsp(
                query,
                subject_packed,
                data,
                kbp,
                search_space,
                evalue_threshold,
                context,
            ) else {
                last_hit[real_diag] = hit.subject_match_end as i32;
                flag[real_diag] = 0;
                hit_len[real_diag] = if two_hits {
                    hit.subject_match_end as i32 - hit.subject_start as i32
                } else {
                    0
                };
                return;
            };
            if let Some((_, _, ss, se)) = trace_target_bounds() {
                if hsp.subject_start <= se && hsp.subject_end >= ss {
                    eprintln!(
                        "[trace-target] exact-hit-saved ctx={} q={}..{} s={}..{} score={}",
                        context,
                        hsp.query_start,
                        hsp.query_end,
                        hsp.subject_start,
                        hsp.subject_end,
                        hsp.score
                    );
                }
            }
            if matches_trace_seed_pair(hit.query_start, hit.subject_start) {
                eprintln!(
                    "[trace-seed] exact-hit-saved ctx={} q={}..{} s={}..{} score={}",
                    context,
                    hsp.query_start,
                    hsp.query_end,
                    hsp.subject_start,
                    hsp.subject_end,
                    hsp.score
                );
            }
            s_end_pos = hsp.subject_end;
            hsps.push(hsp);
        } else {
            hit_ready = false;
            if let Some((qs, qe, ss, se)) = trace_target_bounds() {
                if (hit.query_start as i32) <= qe
                    && (hit.query_start + word_size) as i32 >= qs
                    && (hit.subject_start as i32) <= se
                    && (hit.subject_start + word_size) as i32 >= ss
                {
                    eprintln!(
                        "[trace-target] exact-hit-rejected ctx={} q={}..{} s={}..{}",
                        context,
                        hit.query_start,
                        hit.query_start + word_size,
                        hit.subject_start,
                        hit.subject_start + word_size
                    );
                }
            }
            if matches_trace_seed_pair(hit.query_start, hit.subject_start) {
                eprintln!(
                    "[trace-seed] exact-hit-rejected ctx={} q={}..{} s={}..{}",
                    context,
                    hit.query_start,
                    hit.query_start + word_size,
                    hit.subject_start,
                    hit.subject_start + word_size
                );
            }
        }
    }
    last_hit[real_diag] = s_end_pos;
    flag[real_diag] = i32::from(hit_ready);
    hit_len[real_diag] = if two_hits && !hit_ready {
        s_end_pos - hit.subject_start as i32
    } else {
        0
    };
}

#[allow(clippy::too_many_arguments)]
fn diag_table_extend_initial_hit_packed(
    query: &[u8],
    subject_packed: &[u8],
    subject_len: usize,
    word_size: usize,
    last_hit: &mut [i32],
    flag: &mut [i32],
    hit_len: &mut [i32],
    diag_mask: usize,
    reward: i32,
    penalty: i32,
    x_dropoff: i32,
    kbp: &KarlinBlk,
    search_space: f64,
    evalue_threshold: f64,
    context: i32,
    nucl_score_table: &[i32; 256],
    reduced_nucl_cutoff_score: i32,
    hit: ExactWordHit,
    hsps: &mut Vec<SearchHsp>,
) {
    diag_initial_hit_core_packed(
        query,
        subject_packed,
        subject_len,
        word_size,
        last_hit,
        flag,
        hit_len,
        diag_mask,
        reward,
        penalty,
        x_dropoff,
        kbp,
        search_space,
        evalue_threshold,
        context,
        nucl_score_table,
        reduced_nucl_cutoff_score,
        hit,
        hsps,
    )
}

/// Port boundary for NCBI's direct packed nucleotide extension helper
/// (`s_BlastNaExtendDirect`) before the diagonal-table initial-hit helper.
#[allow(clippy::too_many_arguments)]
fn blast_na_extend_direct_packed(
    query: &[u8],
    subject_packed: &[u8],
    subject_len: usize,
    word_size: usize,
    last_hit: &mut [i32],
    flag: &mut [i32],
    hit_len: &mut [i32],
    diag_mask: usize,
    reward: i32,
    penalty: i32,
    x_dropoff: i32,
    kbp: &KarlinBlk,
    search_space: f64,
    evalue_threshold: f64,
    context: i32,
    nucl_score_table: &[i32; 256],
    reduced_nucl_cutoff_score: i32,
    hit: ExactWordHit,
    hsps: &mut Vec<SearchHsp>,
) {
    diag_table_extend_initial_hit_packed(
        query,
        subject_packed,
        subject_len,
        word_size,
        last_hit,
        flag,
        hit_len,
        diag_mask,
        reward,
        penalty,
        x_dropoff,
        kbp,
        search_space,
        evalue_threshold,
        context,
        nucl_score_table,
        reduced_nucl_cutoff_score,
        hit,
        hsps,
    )
}

#[allow(clippy::too_many_arguments)]
fn diag_hash_extend_initial_hit_packed(
    query: &[u8],
    subject_packed: &[u8],
    subject_len: usize,
    word_size: usize,
    last_hit: &mut [i32],
    flag: &mut [i32],
    hit_len: &mut [i32],
    diag_mask: usize,
    reward: i32,
    penalty: i32,
    x_dropoff: i32,
    kbp: &KarlinBlk,
    search_space: f64,
    evalue_threshold: f64,
    context: i32,
    nucl_score_table: &[i32; 256],
    reduced_nucl_cutoff_score: i32,
    hit: ExactWordHit,
    hsps: &mut Vec<SearchHsp>,
) {
    diag_initial_hit_core_packed(
        query,
        subject_packed,
        subject_len,
        word_size,
        last_hit,
        flag,
        hit_len,
        diag_mask,
        reward,
        penalty,
        x_dropoff,
        kbp,
        search_space,
        evalue_threshold,
        context,
        nucl_score_table,
        reduced_nucl_cutoff_score,
        hit,
        hsps,
    )
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
    let q_start = 0usize;
    let q_range = query.len();
    let rsdl = 4usize - (sp % 4);
    let shifted_sp = sp + rsdl;
    let shifted_qp = qp + rsdl;

    let mut ext_left = 0usize;
    let mut ext_right = 0usize;
    let mut ext_max = word_size
        .saturating_sub(lut_word)
        .min(sp)
        .min(qp.saturating_sub(q_start))
        .saturating_add(rsdl);

    while ext_left < ext_max {
        let q_idx = shifted_qp.saturating_sub(ext_left + 1);
        let s_idx = shifted_sp.saturating_sub(ext_left + 1);
        if q_idx >= query.len() || s_idx >= subject.len() || query[q_idx] != subject[s_idx] {
            break;
        }
        ext_left += 1;
    }
    ext_left = ext_left.min(ext_max).min(qp).min(sp);

    ext_max = word_size
        .saturating_sub(ext_left)
        .min(subject.len().saturating_sub(shifted_sp))
        .min(q_range.saturating_sub(shifted_qp));
    while ext_right < ext_max {
        let q_idx = shifted_qp + ext_right;
        let s_idx = shifted_sp + ext_right;
        if q_idx >= query.len() || s_idx >= subject.len() || query[q_idx] != subject[s_idx] {
            break;
        }
        ext_right += 1;
    }
    ext_right = ext_right.min(ext_max);

    if ext_left + ext_right < word_size {
        return;
    }

    let wq = qp - ext_left;
    let ws = sp - ext_left;
    let real_diag = (ws + (diag_mask + 1) - wq) & diag_mask;
    if (last_hit[real_diag] as usize) <= ws {
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
            last_hit[real_diag] = hsp.subject_end;
            hsps.push(hsp);
        } else {
            last_hit[real_diag] = (ws + word_size) as i32;
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
    last_hit_scratch: &mut [PackedDiagScratch],
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
    last_hit_scratch: &mut [PackedDiagScratch],
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
    last_hit_scratch: &mut [PackedDiagScratch],
    dedup: bool,
) -> Vec<SearchHsp> {
    let mut hsps = Vec::new();
    let nucl_score_table = InitialWordParameters::build_nucl_score_table(reward, penalty);
    let reduced_nucl_cutoff_score =
        ((kbp.evalue_to_raw(evalue_threshold, search_space) as f64) * 0.8) as i32;

    if prepared.is_empty() || subject_len < prepared.word_size {
        return hsps;
    }

    let end = subject_len - prepared.word_size + 1;

    assert_eq!(
        last_hit_scratch.len(),
        prepared.lookups.len(),
        "last-hit scratch must match prepared lookup contexts"
    );

    blast_na_word_finder_scan_subject_packed(
        prepared,
        subject_packed,
        subject_len,
        end,
        reward,
        penalty,
        x_dropoff,
        kbp,
        search_space,
        evalue_threshold,
        &nucl_score_table,
        reduced_nucl_cutoff_score,
        last_hit_scratch,
        &mut hsps,
    );

    if dedup {
        trace_matching_hsps("ungapped-pre-dedup", &hsps);
        hsps.sort_by(score_compare_search_hsps);
        dedup_hsps_with_min_diag_separation(
            &mut hsps,
            min_diag_separation_for_ungapped(prepared.word_size, reward, penalty),
        );
        trace_matching_hsps("ungapped-post-dedup", &hsps);
    } else {
        sort_ungapped_hsps_ncbi_no_dedup(&mut hsps);
        trace_matching_hsps("ungapped-no-dedup", &hsps);
    }
    hsps
}

#[allow(clippy::too_many_arguments)]
fn blast_na_word_finder_scan_subject_packed(
    prepared: &PreparedBlastnQuery<'_>,
    subject_packed: &[u8],
    subject_len: usize,
    end: usize,
    reward: i32,
    penalty: i32,
    x_dropoff: i32,
    kbp: &KarlinBlk,
    search_space: f64,
    evalue_threshold: f64,
    nucl_score_table: &[i32; 256],
    reduced_nucl_cutoff_score: i32,
    last_hit_scratch: &mut [PackedDiagScratch],
    hsps: &mut Vec<SearchHsp>,
) {
    if blast_na_hash_lookup_scan_subject_packed(
        prepared,
        subject_packed,
        subject_len,
        end,
        reward,
        penalty,
        x_dropoff,
        kbp,
        search_space,
        evalue_threshold,
        nucl_score_table,
        reduced_nucl_cutoff_score,
        last_hit_scratch,
        hsps,
    ) {
        return;
    }

    for (lookup, scratch) in prepared.lookups.iter().zip(last_hit_scratch.iter_mut()) {
        scratch.resize_and_clear(lookup.diag_array_len);
        blast_na_word_finder_scan_lookup_packed(
            prepared.word_size,
            lookup,
            scratch,
            subject_packed,
            subject_len,
            end,
            reward,
            penalty,
            x_dropoff,
            kbp,
            search_space,
            evalue_threshold,
            nucl_score_table,
            reduced_nucl_cutoff_score,
            hsps,
        );
    }
}

#[allow(clippy::too_many_arguments)]
fn blast_na_hash_lookup_scan_subject_packed(
    prepared: &PreparedBlastnQuery<'_>,
    subject_packed: &[u8],
    subject_len: usize,
    end: usize,
    reward: i32,
    penalty: i32,
    x_dropoff: i32,
    kbp: &KarlinBlk,
    search_space: f64,
    evalue_threshold: f64,
    nucl_score_table: &[i32; 256],
    reduced_nucl_cutoff_score: i32,
    last_hit_scratch: &mut [PackedDiagScratch],
    hsps: &mut Vec<SearchHsp>,
) -> bool {
    if !(prepared.lookups.len() == 2
        && prepared.lookups[0].lut_word == 6
        && prepared.lookups[1].lut_word == 6
        && prepared.lookups[0].scan_start == 1
        && prepared.lookups[1].scan_start == 1
        && prepared.lookups[0].scan_step == 2
        && prepared.lookups[1].scan_step == 2
        && prepared.paired_pv.is_some())
    {
        return false;
    }

    let (left_scratch, right_scratch) = last_hit_scratch.split_at_mut(1);
    let lookup0 = &prepared.lookups[0];
    let lookup1 = &prepared.lookups[1];
    let scratch0 = &mut left_scratch[0];
    let scratch1 = &mut right_scratch[0];

    scratch0.resize_and_clear(lookup0.diag_array_len);
    scratch1.resize_and_clear(lookup1.diag_array_len);

    scan_byte_oriented_lut6_step2_pair(
        subject_packed,
        subject_len,
        end,
        lookup0,
        &mut scratch0.last_hit,
        &mut scratch0.flag,
        &mut scratch0.hit_len,
        lookup1,
        &mut scratch1.last_hit,
        &mut scratch1.flag,
        &mut scratch1.hit_len,
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
        nucl_score_table,
        reduced_nucl_cutoff_score,
        hsps,
    );
    true
}

#[allow(clippy::too_many_arguments)]
fn blast_na_word_finder_scan_lookup_packed(
    word_size: usize,
    lookup: &NaLookup<'_>,
    scratch: &mut PackedDiagScratch,
    subject_packed: &[u8],
    subject_len: usize,
    end: usize,
    reward: i32,
    penalty: i32,
    x_dropoff: i32,
    kbp: &KarlinBlk,
    search_space: f64,
    evalue_threshold: f64,
    nucl_score_table: &[i32; 256],
    reduced_nucl_cutoff_score: i32,
    hsps: &mut Vec<SearchHsp>,
) {
    if lookup.lut_word == 8 && word_size >= 8 {
        let scan_step = word_size - lookup.lut_word + 1;
        if scan_step == 4 && lookup.scan_start == 0 {
            scan_step4_unrolled(
                subject_packed,
                subject_len,
                word_size,
                lookup.query_nomask,
                &lookup.lut,
                &lookup.next,
                &lookup.pv,
                &mut scratch.last_hit,
                &mut scratch.flag,
                &mut scratch.hit_len,
                lookup.diag_mask,
                reward,
                penalty,
                x_dropoff,
                kbp,
                search_space,
                evalue_threshold,
                lookup.context,
                nucl_score_table,
                reduced_nucl_cutoff_score,
                hsps,
            );
        } else if lookup.scan_step % 4 == 1 {
            scan_step1(
                subject_packed,
                subject_len,
                end,
                lookup.lut_word,
                lookup.lut_mask,
                lookup.scan_start,
                lookup.scan_step,
                lookup.query_nomask,
                word_size,
                &lookup.lut,
                &lookup.next,
                &lookup.pv,
                &mut scratch.last_hit,
                &mut scratch.flag,
                &mut scratch.hit_len,
                lookup.diag_mask,
                reward,
                penalty,
                x_dropoff,
                kbp,
                search_space,
                evalue_threshold,
                lookup.context,
                nucl_score_table,
                reduced_nucl_cutoff_score,
                hsps,
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
                lookup.query_nomask,
                word_size,
                &lookup.lut,
                &lookup.next,
                &lookup.pv,
                &mut scratch.last_hit,
                &mut scratch.flag,
                &mut scratch.hit_len,
                lookup.diag_mask,
                reward,
                penalty,
                x_dropoff,
                kbp,
                search_space,
                evalue_threshold,
                lookup.context,
                nucl_score_table,
                reduced_nucl_cutoff_score,
                hsps,
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
            lookup.query_nomask,
            word_size,
            &lookup.lut,
            &lookup.next,
            &lookup.pv,
            &mut scratch.last_hit,
            &mut scratch.flag,
            &mut scratch.hit_len,
            lookup.diag_mask,
            reward,
            penalty,
            x_dropoff,
            kbp,
            search_space,
            evalue_threshold,
            lookup.context,
            nucl_score_table,
            reduced_nucl_cutoff_score,
            hsps,
        );
    }
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
        let needed = lookup.diag_array_len;
        if last_hit.len() != needed {
            last_hit.resize(needed, 0);
        } else {
            last_hit.fill(0);
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
                            lookup.query_nomask,
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

fn trace_matching_hsps(label: &str, hsps: &[SearchHsp]) {
    for hsp in hsps {
        trace_hsp(label, hsp);
    }
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
    flag: &mut [i32],
    hit_len: &mut [i32],
    diag_mask: usize,
    reward: i32,
    penalty: i32,
    x_dropoff: i32,
    kbp: &KarlinBlk,
    search_space: f64,
    evalue_threshold: f64,
    context: i32,
    nucl_score_table: &[i32; 256],
    reduced_nucl_cutoff_score: i32,
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
            if let Some(positions) = trace_scan_positions() {
                if positions.contains(&$s_pos) {
                    let mut chain = Vec::new();
                    let mut q_pos = lut[h];
                    while q_pos >= 0 && chain.len() < 8 {
                        chain.push(q_pos);
                        q_pos = next[q_pos as usize];
                    }
                    eprintln!(
                        "[trace-scan] fn=step4-unrolled ctx={} s_pos={} hash={} pv={} chain={:?}",
                        context,
                        $s_pos,
                        h,
                        pv_has_hash(pv, h),
                        chain
                    );
                }
            }
            if pv_has_hash(pv, h) {
                small_na_extend_packed(
                    query,
                    subject_packed,
                    subject_len,
                    word_size,
                    8,
                    lut,
                    next,
                    last_hit,
                    flag,
                    hit_len,
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
                    nucl_score_table,
                    reduced_nucl_cutoff_score,
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
    flag: &mut [i32],
    hit_len: &mut [i32],
    diag_mask: usize,
    reward: i32,
    penalty: i32,
    x_dropoff: i32,
    kbp: &KarlinBlk,
    search_space: f64,
    evalue_threshold: f64,
    context: i32,
    nucl_score_table: &[i32; 256],
    reduced_nucl_cutoff_score: i32,
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
            flag,
            hit_len,
            diag_mask,
            reward,
            penalty,
            x_dropoff,
            kbp,
            search_space,
            evalue_threshold,
            context,
            nucl_score_table,
            reduced_nucl_cutoff_score,
            hsps,
        );
        return;
    }

    if scan_step != 1 {
        let mut s_pos = scan_start;
        while s_pos < end {
            let h = packed_hash_at(subject_packed, s_pos, lut_word, lut_mask);
            if pv[h >> 6] & (1u64 << (h & 63)) != 0 {
                small_na_extend_packed(
                    query,
                    subject_packed,
                    subject_len,
                    word_size,
                    lut_word,
                    lut,
                    next,
                    last_hit,
                    flag,
                    hit_len,
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
                    nucl_score_table,
                    reduced_nucl_cutoff_score,
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
                    if let Some(positions) = trace_scan_positions() {
                        if positions.contains(&pos) {
                            let mut chain = Vec::new();
                            let mut q_pos = lut[h];
                            while q_pos >= 0 && chain.len() < 8 {
                                chain.push(q_pos);
                                q_pos = next[q_pos as usize];
                            }
                            eprintln!(
                                "[trace-scan] fn=step4-main ctx={} s_pos={} hash={} pv={} chain={:?}",
                                context,
                                pos,
                                h,
                                pv_has_hash(pv, h),
                                chain
                            );
                        }
                    }
                    if pv_has_hash(pv, h) {
                        small_na_extend_packed(
                            query,
                            subject_packed,
                            subject_len,
                            word_size,
                            lut_word,
                            lut,
                            next,
                            last_hit,
                            flag,
                            hit_len,
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
                            nucl_score_table,
                            reduced_nucl_cutoff_score,
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
        if let Some(positions) = trace_scan_positions() {
            if positions.contains(&s_pos) {
                let mut chain = Vec::new();
                let mut q_pos = lut[h];
                while q_pos >= 0 && chain.len() < 8 {
                    chain.push(q_pos);
                    q_pos = next[q_pos as usize];
                }
                eprintln!(
                    "[trace-scan] fn=step4-tail ctx={} s_pos={} hash={} pv={} chain={:?}",
                    context,
                    s_pos,
                    h,
                    pv[h >> 6] & (1u64 << (h & 63)) != 0,
                    chain
                );
            }
        }
        if pv[h >> 6] & (1u64 << (h & 63)) != 0 {
            small_na_extend_packed(
                query,
                subject_packed,
                subject_len,
                word_size,
                lut_word,
                lut,
                next,
                last_hit,
                flag,
                hit_len,
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
                nucl_score_table,
                reduced_nucl_cutoff_score,
                hsps,
            );
        }
        s_pos += 1;
    }
}

#[allow(clippy::too_many_arguments, dead_code)]
fn scan_byte_oriented_lut8_mod1(
    subject_packed: &[u8],
    subject_len: usize,
    end: usize,
    lut_mask: u32,
    scan_start: usize,
    scan_step: usize,
    query: &[u8],
    word_size: usize,
    lut: &[i32],
    next: &[i32],
    pv: &[u64],
    last_hit: &mut [i32],
    flag: &mut [i32],
    hit_len: &mut [i32],
    diag_mask: usize,
    reward: i32,
    penalty: i32,
    x_dropoff: i32,
    kbp: &KarlinBlk,
    search_space: f64,
    evalue_threshold: f64,
    context: i32,
    nucl_score_table: &[i32; 256],
    reduced_nucl_cutoff_score: i32,
    hsps: &mut Vec<SearchHsp>,
) {
    if end == 0 || subject_packed.is_empty() {
        return;
    }

    let scan_step_byte = scan_step / 4;
    let s_len = subject_packed.len();
    let mut s_pos = scan_start;
    let mut byte_idx = scan_start / 4;
    let mut state = scan_start % 4;

    loop {
        if s_pos >= end {
            break;
        }

        let (h, byte_advance) = match state {
            0 => {
                if byte_idx + 1 >= s_len {
                    break;
                }
                let index = ((subject_packed[byte_idx] as u32) << 8) | subject_packed[byte_idx + 1] as u32;
                (index as usize, scan_step_byte)
            }
            1 => {
                if byte_idx + 2 >= s_len {
                    break;
                }
                let index = ((subject_packed[byte_idx] as u32) << 16)
                    | ((subject_packed[byte_idx + 1] as u32) << 8)
                    | subject_packed[byte_idx + 2] as u32;
                (((index >> 6) & lut_mask) as usize, scan_step_byte)
            }
            2 => {
                if byte_idx + 2 >= s_len {
                    break;
                }
                let index = ((subject_packed[byte_idx] as u32) << 16)
                    | ((subject_packed[byte_idx + 1] as u32) << 8)
                    | subject_packed[byte_idx + 2] as u32;
                (((index >> 4) & lut_mask) as usize, scan_step_byte)
            }
            3 => {
                if byte_idx + 2 >= s_len {
                    break;
                }
                let index = ((subject_packed[byte_idx] as u32) << 16)
                    | ((subject_packed[byte_idx + 1] as u32) << 8)
                    | subject_packed[byte_idx + 2] as u32;
                (((index >> 2) & lut_mask) as usize, scan_step_byte + 1)
            }
            _ => unreachable!(),
        };

        if let Some(positions) = trace_scan_positions() {
            if positions.contains(&s_pos) {
                let mut chain = Vec::new();
                let mut q_pos = lut[h];
                while q_pos >= 0 && chain.len() < 8 {
                    chain.push(q_pos);
                    q_pos = next[q_pos as usize];
                }
                let scalar_h = packed_hash_at(subject_packed, s_pos, 8, lut_mask);
                eprintln!(
                    "[trace-scan] fn=lut8_mod1 ctx={} s_pos={} state={} hash={} scalar_hash={} pv={} chain={:?}",
                    context,
                    s_pos,
                    state,
                    h,
                    scalar_h,
                    pv_has_hash(pv, h),
                    chain
                );
            }
        }

        if std::env::var_os("BLAST_RS_TRACE_ROW3_SCAN").is_some()
            && matches!(s_pos, 3_897_726 | 3_897_768 | 3_897_894 | 3_897_999)
        {
            let mut chain = Vec::new();
            let mut q_pos = lut[h];
            while q_pos >= 0 && chain.len() < 8 {
                chain.push(q_pos);
                q_pos = next[q_pos as usize];
            }
            let scalar_h = packed_hash_at(subject_packed, s_pos, 8, lut_mask);
            eprintln!(
                "[row3-scan] ctx={} s_pos={} state={} hash={} scalar_hash={} pv={} chain={:?}",
                context,
                s_pos,
                state,
                h,
                scalar_h,
                pv_has_hash(pv, h),
                chain
            );
        }

        if pv_has_hash(pv, h) {
            small_na_extend_packed(
                query,
                subject_packed,
                subject_len,
                word_size,
                8,
                lut,
                next,
                last_hit,
                flag,
                hit_len,
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
                nucl_score_table,
                reduced_nucl_cutoff_score,
                hsps,
            );
        }

        byte_idx += byte_advance;
        s_pos += scan_step;
        state = (state + scan_step) & 3;
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
    flag: &mut [i32],
    hit_len: &mut [i32],
    diag_mask: usize,
    reward: i32,
    penalty: i32,
    x_dropoff: i32,
    kbp: &KarlinBlk,
    search_space: f64,
    evalue_threshold: f64,
    context: i32,
    nucl_score_table: &[i32; 256],
    reduced_nucl_cutoff_score: i32,
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
                    small_na_extend_packed(
                        query,
                        subject_packed,
                        subject_len,
                        word_size,
                        6,
                        lut,
                        next,
                        last_hit,
                        flag,
                        hit_len,
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
                        nucl_score_table,
                        reduced_nucl_cutoff_score,
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
            small_na_extend_packed(
                query,
                subject_packed,
                subject_len,
                word_size,
                6,
                lut,
                next,
                last_hit,
                flag,
                hit_len,
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
                nucl_score_table,
                reduced_nucl_cutoff_score,
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
    flag0: &mut [i32],
    hit_len0: &mut [i32],
    lookup1: &NaLookup<'_>,
    last_hit1: &mut [i32],
    flag1: &mut [i32],
    hit_len1: &mut [i32],
    paired_pv: &[u64],
    word_size: usize,
    reward: i32,
    penalty: i32,
    x_dropoff: i32,
    kbp: &KarlinBlk,
    search_space: f64,
    evalue_threshold: f64,
    nucl_score_table: &[i32; 256],
    reduced_nucl_cutoff_score: i32,
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
                        small_na_extend_packed(
                            lookup0.query,
                            subject_packed,
                            subject_len,
                            word_size,
                            6,
                            &lookup0.lut,
                            &lookup0.next,
                            last_hit0,
                            flag0,
                            hit_len0,
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
                            nucl_score_table,
                            reduced_nucl_cutoff_score,
                            hsps,
                        );
                    }
                    if pv_bucket_unchecked(&lookup1.pv, pv_bucket) & pv_bit != 0 {
                        small_na_extend_packed(
                            lookup1.query,
                            subject_packed,
                            subject_len,
                            word_size,
                            6,
                            &lookup1.lut,
                            &lookup1.next,
                            last_hit1,
                            flag1,
                            hit_len1,
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
                            nucl_score_table,
                            reduced_nucl_cutoff_score,
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
                small_na_extend_packed(
                    lookup0.query,
                    subject_packed,
                    subject_len,
                    word_size,
                    6,
                    &lookup0.lut,
                    &lookup0.next,
                    last_hit0,
                    flag0,
                    hit_len0,
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
                    nucl_score_table,
                    reduced_nucl_cutoff_score,
                    hsps,
                );
            }
            if lookup1.pv[pv_bucket] & pv_bit != 0 {
                small_na_extend_packed(
                    lookup1.query,
                    subject_packed,
                    subject_len,
                    word_size,
                    6,
                    &lookup1.lut,
                    &lookup1.next,
                    last_hit1,
                    flag1,
                    hit_len1,
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
                    nucl_score_table,
                    reduced_nucl_cutoff_score,
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
    flag: &mut [i32],
    hit_len: &mut [i32],
    diag_mask: usize,
    reward: i32,
    penalty: i32,
    x_dropoff: i32,
    kbp: &KarlinBlk,
    search_space: f64,
    evalue_threshold: f64,
    context: i32,
    nucl_score_table: &[i32; 256],
    reduced_nucl_cutoff_score: i32,
    hsps: &mut Vec<SearchHsp>,
) {
    let subject_len = _subject_len;
    if scan_step != 1 || scan_start != 0 {
        let mut s_pos = scan_start;
        while s_pos < end {
            let h = packed_hash_at(subject_packed, s_pos, lut_word, lut_mask);
            if let Some(positions) = trace_scan_positions() {
                if positions.contains(&s_pos) {
                    let mut chain = Vec::new();
                    let mut q_pos = lut[h];
                    while q_pos >= 0 && chain.len() < 8 {
                        chain.push(q_pos);
                        q_pos = next[q_pos as usize];
                    }
                    eprintln!(
                        "[trace-scan] fn=step1-generic ctx={} s_pos={} hash={} pv={} chain={:?}",
                        context,
                        s_pos,
                        h,
                        pv[h >> 6] & (1u64 << (h & 63)) != 0,
                        chain
                    );
                }
            }
            if pv[h >> 6] & (1u64 << (h & 63)) != 0 {
                small_na_extend_packed(
                    query,
                    subject_packed,
                    subject_len,
                    word_size,
                    lut_word,
                    lut,
                    next,
                    last_hit,
                    flag,
                    hit_len,
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
                    nucl_score_table,
                    reduced_nucl_cutoff_score,
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
        if let Some(positions) = trace_scan_positions() {
            if positions.contains(&s_pos) {
                let mut chain = Vec::new();
                let mut q_pos = lut[h];
                while q_pos >= 0 && chain.len() < 8 {
                    chain.push(q_pos);
                    q_pos = next[q_pos as usize];
                }
                eprintln!(
                    "[trace-scan] fn=step1 ctx={} s_pos={} hash={} pv={} chain={:?}",
                    context,
                    s_pos,
                    h,
                    pv[h >> 6] & (1u64 << (h & 63)) != 0,
                    chain
                );
            }
        }
        if std::env::var_os("BLAST_RS_TRACE_ROW3_SCAN").is_some()
            && matches!(s_pos, 3_897_726 | 3_897_768 | 3_897_894 | 3_897_999)
        {
            let mut chain = Vec::new();
            let mut q_pos = lut[h];
            while q_pos >= 0 && chain.len() < 8 {
                chain.push(q_pos);
                q_pos = next[q_pos as usize];
            }
            eprintln!(
                "[row3-scan] ctx={} s_pos={} hash={} pv={} chain={:?}",
                context,
                s_pos,
                h,
                pv[h >> 6] & (1u64 << (h & 63)) != 0,
                chain
            );
        }
        if pv[h >> 6] & (1u64 << (h & 63)) != 0 {
            small_na_extend_packed(
                query,
                subject_packed,
                subject_len,
                word_size,
                lut_word,
                lut,
                next,
                last_hit,
                flag,
                hit_len,
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
                nucl_score_table,
                reduced_nucl_cutoff_score,
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
fn pack_query_byte(query: &[u8], pos: usize) -> u8 {
    (query[pos] << 6) | (query[pos + 1] << 4) | (query[pos + 2] << 2) | query[pos + 3]
}

struct PackedUngappedData {
    q_start: usize,
    s_start: usize,
    length: usize,
    score: i32,
}

#[derive(Clone, Copy, Debug)]
struct DecodedUngappedData {
    q_start: usize,
    s_start: usize,
    length: usize,
    score: i32,
}

fn build_packed_hsp(
    query: &[u8],
    subject_packed: &[u8],
    data: PackedUngappedData,
    kbp: &KarlinBlk,
    search_space: f64,
    evalue_threshold: f64,
    context: i32,
) -> Option<SearchHsp> {
    if data.score <= 0 || data.q_start + data.length > query.len() {
        return None;
    }
    let evalue = kbp.raw_to_evalue(data.score, search_space);
    if evalue > evalue_threshold {
        return None;
    }

    let q_start = data.q_start as i32;
    let q_end = (data.q_start + data.length) as i32;
    let s_start = data.s_start as i32;
    let s_end = (data.s_start + data.length) as i32;
    let align_len = data.length as i32;

    let mut num_ident = 0i32;
    let mut qseq = Vec::with_capacity(data.length);
    let mut sseq = Vec::with_capacity(data.length);
    for i in 0..data.length {
        let qb = query[data.q_start + i];
        let sb = packed_base_at(subject_packed, data.s_start + i);
        if qb == sb {
            num_ident += 1;
        }
        qseq.push(blastna_to_iupac_byte(qb));
        sseq.push(blastna_to_iupac_byte(sb));
    }
    let qseq_str = unsafe { String::from_utf8_unchecked(qseq) };
    let sseq_str = unsafe { String::from_utf8_unchecked(sseq) };

    let hsp = SearchHsp {
        query_start: q_start,
        query_end: q_end,
        subject_start: s_start,
        subject_end: s_end,
        score: data.score,
        bit_score: kbp.raw_to_bit(data.score),
        evalue,
        num_ident,
        align_length: align_len,
        mismatches: align_len - num_ident,
        gap_opens: 0,
        context,
        qseq: Some(qseq_str),
        sseq: Some(sseq_str),
    };
    trace_hsp("ungapped-packed", &hsp);
    Some(hsp)
}

#[inline]
fn build_decoded_hsp(
    query: &[u8],
    subject: &[u8],
    ungapped: DecodedUngappedData,
    kbp: &KarlinBlk,
    search_space: f64,
    evalue_threshold: f64,
    context: i32,
) -> Option<SearchHsp> {
    if ungapped.score <= 0 {
        return None;
    }

    let evalue = kbp.raw_to_evalue(ungapped.score, search_space);
    if evalue > evalue_threshold {
        return None;
    }
    let bit_score = kbp.raw_to_bit(ungapped.score);

    let q_start = ungapped.q_start as i32;
    let q_end = (ungapped.q_start + ungapped.length) as i32;
    let s_start = ungapped.s_start as i32;
    let s_end = (ungapped.s_start + ungapped.length) as i32;
    let align_len = q_end - q_start;

    let mut num_ident = 0;
    let mut qseq = Vec::with_capacity(align_len as usize);
    let mut sseq = Vec::with_capacity(align_len as usize);
    for i in 0..align_len as usize {
        let qb = query[ungapped.q_start + i];
        let sb = subject[ungapped.s_start + i];
        if qb == sb {
            num_ident += 1;
        }
        qseq.push(blastna_to_iupac_byte(qb));
        sseq.push(blastna_to_iupac_byte(sb));
    }
    let qseq_str = unsafe { String::from_utf8_unchecked(qseq) };
    let sseq_str = unsafe { String::from_utf8_unchecked(sseq) };

    let hsp = SearchHsp {
        query_start: q_start,
        query_end: q_end,
        subject_start: s_start,
        subject_end: s_end,
        score: ungapped.score,
        bit_score,
        evalue,
        num_ident,
        align_length: align_len,
        mismatches: align_len - num_ident,
        gap_opens: 0,
        context,
        qseq: Some(qseq_str),
        sseq: Some(sseq_str),
    };
    trace_hsp("ungapped-decoded", &hsp);
    Some(hsp)
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
#[allow(dead_code)]
#[inline]
fn extend_seed_packed(
    query: &[u8],
    subject_packed: &[u8],
    subject_len: usize,
    q_seed: usize,
    s_seed: usize,
    s_match_end: usize,
    reward: i32,
    penalty: i32,
    x_dropoff: i32,
    word_size: usize,
    nucl_score_table: &[i32; 256],
    reduced_nucl_cutoff_score: i32,
) -> Option<PackedUngappedData> {
    if word_size < 11 {
        return extend_seed_packed_exact(
            query,
            subject_packed,
            subject_len,
            q_seed,
            s_seed,
            reward,
            penalty,
            x_dropoff,
        );
    }

    let approx = extend_seed_packed_approx(
        query,
        subject_packed,
        subject_len,
        q_seed,
        s_seed,
        s_match_end,
        x_dropoff,
        nucl_score_table,
    )?;
    if approx.score >= reduced_nucl_cutoff_score {
        extend_seed_packed_exact(
            query,
            subject_packed,
            subject_len,
            q_seed,
            s_seed,
            reward,
            penalty,
            x_dropoff,
        )
    } else {
        Some(approx)
    }
}

#[inline]
fn extend_seed_packed_approx(
    query: &[u8],
    subject_packed: &[u8],
    subject_len: usize,
    q_seed: usize,
    s_seed: usize,
    s_match_end: usize,
    x_dropoff: i32,
    nucl_score_table: &[i32; 256],
) -> Option<PackedUngappedData> {
    let x_dropoff_neg = -x_dropoff;
    let len = (4 - (s_seed % 4)) % 4;
    let q_ext = q_seed + len;
    let s_ext = s_seed + len;

    let mut score = 0i32;
    let mut sum = 0i32;
    let mut new_q = q_ext;
    let mut q = q_ext;
    let mut left_len = q_ext.min(s_ext) / 4;
    while left_len > 0 {
        let s_byte = subject_packed[s_ext / 4 - left_len];
        let q_byte = pack_query_byte(query, q - 4);
        sum += nucl_score_table[(q_byte ^ s_byte) as usize];
        if sum > 0 {
            new_q = q - 4;
            score += sum;
            sum = 0;
        }
        if sum < x_dropoff_neg {
            break;
        }
        q -= 4;
        left_len -= 1;
    }

    let q_start = new_q;
    let s_start = s_ext - (q_ext - q_start);

    let mut right_q = q_ext;
    let mut right_sum = 0i32;
    let mut right_new_q = q_ext;
    let mut right_len = (query.len() - q_ext).min(subject_len - s_ext) / 4;
    let mut s_byte_idx = s_ext / 4;
    while right_len > 0 {
        let s_byte = subject_packed[s_byte_idx];
        let q_byte = pack_query_byte(query, right_q);
        right_sum += nucl_score_table[(q_byte ^ s_byte) as usize];
        if right_sum > 0 {
            right_new_q = right_q + 3;
            score += right_sum;
            right_sum = 0;
        }
        if right_sum < x_dropoff_neg {
            break;
        }
        right_q += 4;
        s_byte_idx += 1;
        right_len -= 1;
    }

    let length = s_match_end
        .saturating_sub(s_start)
        .max(right_new_q.saturating_sub(q_start) + 1);
    Some(PackedUngappedData {
        q_start,
        s_start,
        length,
        score,
    })
}

#[inline]
fn extend_seed_packed_exact(
    query: &[u8],
    subject_packed: &[u8],
    subject_len: usize,
    q_seed: usize,
    s_seed: usize,
    reward: i32,
    penalty: i32,
    x_dropoff: i32,
) -> Option<PackedUngappedData> {
    if std::env::var_os("BLAST_RS_PACKED_EXACT_DECODED").is_some() {
        let mut subject = Vec::with_capacity(subject_len);
        for i in 0..subject_len {
            subject.push(packed_base_at(subject_packed, i));
        }
        return extend_seed_data(query, &subject, q_seed, s_seed, reward, penalty, x_dropoff)
            .map(|d| PackedUngappedData {
                q_start: d.q_start,
                s_start: d.s_start,
                length: d.length,
                score: d.score,
            });
    }

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

    Some(PackedUngappedData {
        q_start: q_seed - best_left,
        s_start: s_seed - best_left,
        length: best_left + best_right,
        score: score + left_score,
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
fn extend_seed_data(
    query: &[u8],
    subject: &[u8],
    q_seed: usize,
    s_seed: usize,
    reward: i32,
    penalty: i32,
    x_dropoff: i32,
)-> Option<DecodedUngappedData> {
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

    Some(DecodedUngappedData {
        q_start: q_seed - best_left,
        s_start: s_seed - best_left,
        length: best_left + best_right,
        score: total_score,
    })
}

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
    let ungapped = extend_seed_data(query, subject, q_seed, s_seed, reward, penalty, x_dropoff)?;
    build_decoded_hsp(
        query,
        subject,
        ungapped,
        kbp,
        search_space,
        evalue_threshold,
        context,
    )
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
        let contained = tree.is_contained_with_min_diag_separation(qi, si, min_diag_separation);
        let overlap = tree.has_significant_overlap_with_min_diag_separation(
            qi,
            si,
            0.99,
            min_diag_separation,
        );
        if !contained && !overlap {
            keep[i] = true;
            tree.insert(qi, si, hsp.score);
            trace_hsp("ungapped-dedup-keep", hsp);
        } else {
            trace_hsp("ungapped-dedup-drop", hsp);
        }
    }
    let mut idx = 0;
    hsps.retain(|_| {
        let k = keep[idx];
        idx += 1;
        k
    });
}

fn purge_common_endpoint_hsps(hsps: &mut Vec<SearchHsp>) {
    if hsps.len() <= 1 {
        return;
    }

    hsps.sort_by(|a, b| {
        a.context
            .cmp(&b.context)
            .then_with(|| a.query_start.cmp(&b.query_start))
            .then_with(|| a.subject_start.cmp(&b.subject_start))
            .then_with(|| b.score.cmp(&a.score))
            .then_with(|| b.query_end.cmp(&a.query_end))
            .then_with(|| b.subject_end.cmp(&a.subject_end))
    });
    hsps.dedup_by(|b, a| {
        a.context == b.context
            && a.query_start == b.query_start
            && a.subject_start == b.subject_start
    });

    hsps.sort_by(|a, b| {
        a.context
            .cmp(&b.context)
            .then_with(|| a.query_end.cmp(&b.query_end))
            .then_with(|| a.subject_end.cmp(&b.subject_end))
            .then_with(|| b.score.cmp(&a.score))
            .then_with(|| b.query_start.cmp(&a.query_start))
            .then_with(|| b.subject_start.cmp(&a.subject_start))
    });
    hsps.dedup_by(|b, a| {
        a.context == b.context && a.query_end == b.query_end && a.subject_end == b.subject_end
    });
}

fn candidate_context(candidate: &GappedCandidate) -> i32 {
    match candidate {
        GappedCandidate::Final(hsp) => hsp.context,
        GappedCandidate::Traceback { context, .. } => *context,
    }
}

fn candidate_query_start(candidate: &GappedCandidate) -> i32 {
    match candidate {
        GappedCandidate::Final(hsp) => hsp.query_start,
        GappedCandidate::Traceback { tb, .. } => tb.query_start as i32,
    }
}

fn candidate_query_end(candidate: &GappedCandidate) -> i32 {
    match candidate {
        GappedCandidate::Final(hsp) => hsp.query_end,
        GappedCandidate::Traceback { tb, .. } => tb.query_end as i32,
    }
}

fn candidate_subject_start(candidate: &GappedCandidate) -> i32 {
    match candidate {
        GappedCandidate::Final(hsp) => hsp.subject_start,
        GappedCandidate::Traceback { tb, .. } => tb.subject_start as i32,
    }
}

fn candidate_subject_end(candidate: &GappedCandidate) -> i32 {
    match candidate {
        GappedCandidate::Final(hsp) => hsp.subject_end,
        GappedCandidate::Traceback { tb, .. } => tb.subject_end as i32,
    }
}

fn candidate_score(candidate: &GappedCandidate) -> i32 {
    match candidate {
        GappedCandidate::Final(hsp) => hsp.score,
        GappedCandidate::Traceback { tb, .. } => tb.score,
    }
}

fn candidate_query_for_context<'a>(
    context: i32,
    query_plus_nomask: &'a [u8],
    query_minus_nomask: &'a [u8],
) -> &'a [u8] {
    if context == 0 {
        query_plus_nomask
    } else {
        query_minus_nomask
    }
}

fn candidate_traceback_quality(
    candidate: &GappedCandidate,
    query_plus_nomask: &[u8],
    query_minus_nomask: &[u8],
    subject: &[u8],
) -> Option<(i32, i32, i32)> {
    let GappedCandidate::Traceback { context, tb } = candidate else {
        return None;
    };
    let query = candidate_query_for_context(*context, query_plus_nomask, query_minus_nomask);
    let q_slice = &query[tb.query_start..tb.query_end];
    let s_slice = &subject[tb.subject_start..tb.subject_end];
    let (align_len, num_ident, gap_opens) = tb.edit_script.count_identities(q_slice, s_slice);
    Some((num_ident, -gap_opens, -align_len))
}

fn collapse_exact_traceback_duplicates(
    candidates: &mut Vec<GappedCandidate>,
    query_plus_nomask: &[u8],
    query_minus_nomask: &[u8],
    subject: &[u8],
) {
    let mut owned = std::mem::take(candidates);
    owned.sort_by(|a, b| {
        candidate_context(a)
            .cmp(&candidate_context(b))
            .then_with(|| candidate_query_start(a).cmp(&candidate_query_start(b)))
            .then_with(|| candidate_query_end(a).cmp(&candidate_query_end(b)))
            .then_with(|| candidate_subject_start(a).cmp(&candidate_subject_start(b)))
            .then_with(|| candidate_subject_end(a).cmp(&candidate_subject_end(b)))
            .then_with(|| candidate_score(a).cmp(&candidate_score(b)))
    });

    let mut deduped = Vec::with_capacity(owned.len());
    let mut group: Vec<GappedCandidate> = Vec::new();

    for candidate in owned.into_iter() {
        let same_group = group.first().is_some_and(|first| {
            candidate_context(first) == candidate_context(&candidate)
                && candidate_query_start(first) == candidate_query_start(&candidate)
                && candidate_query_end(first) == candidate_query_end(&candidate)
                && candidate_subject_start(first) == candidate_subject_start(&candidate)
                && candidate_subject_end(first) == candidate_subject_end(&candidate)
                && candidate_score(first) == candidate_score(&candidate)
        });
        if group.is_empty() || same_group {
            group.push(candidate);
        } else {
            if group.len() > 1 && candidate_matches_trace_target(&group[0]) {
                eprintln!(
                    "[trace-target] collapse group size={} q={}..{} s={}..{} score={}",
                    group.len(),
                    candidate_query_start(&group[0]),
                    candidate_query_end(&group[0]),
                    candidate_subject_start(&group[0]),
                    candidate_subject_end(&group[0]),
                    candidate_score(&group[0])
                );
                for (idx, member) in group.iter().enumerate() {
                    eprintln!(
                        "[trace-target] collapse member#{} quality={:?}",
                        idx,
                        candidate_traceback_quality(
                            member,
                            query_plus_nomask,
                            query_minus_nomask,
                            subject,
                        )
                    );
                    trace_candidate("collapse-member", member);
                }
            }
            let mut best = 0usize;
            for j in 1..group.len() {
                let best_quality = candidate_traceback_quality(
                    &group[best],
                    query_plus_nomask,
                    query_minus_nomask,
                    subject,
                );
                let cand_quality = candidate_traceback_quality(
                    &group[j],
                    query_plus_nomask,
                    query_minus_nomask,
                    subject,
                );
                if cand_quality > best_quality {
                    best = j;
                }
            }
            deduped.push(group.swap_remove(best));
            group.clear();
            group.push(candidate);
        }
    }

    if !group.is_empty() {
        if group.len() > 1 && candidate_matches_trace_target(&group[0]) {
            eprintln!(
                "[trace-target] collapse group size={} q={}..{} s={}..{} score={}",
                group.len(),
                candidate_query_start(&group[0]),
                candidate_query_end(&group[0]),
                candidate_subject_start(&group[0]),
                candidate_subject_end(&group[0]),
                candidate_score(&group[0])
            );
            for (idx, member) in group.iter().enumerate() {
                eprintln!(
                    "[trace-target] collapse member#{} quality={:?}",
                    idx,
                    candidate_traceback_quality(
                        member,
                        query_plus_nomask,
                        query_minus_nomask,
                        subject,
                    )
                );
                trace_candidate("collapse-member", member);
            }
        }
        let mut best = 0usize;
        for j in 1..group.len() {
            let best_quality = candidate_traceback_quality(
                &group[best],
                query_plus_nomask,
                query_minus_nomask,
                subject,
            );
            let cand_quality = candidate_traceback_quality(
                &group[j],
                query_plus_nomask,
                query_minus_nomask,
                subject,
            );
            if cand_quality > best_quality {
                best = j;
            }
        }
        deduped.push(group.swap_remove(best));
    }
    *candidates = deduped;
}

fn trace_target_bounds() -> Option<(i32, i32, i32, i32)> {
    let raw = std::env::var("BLAST_RS_TRACE_TARGET").ok()?;
    let mut parts = raw.split(',').map(str::trim);
    let qs = parts.next()?.parse().ok()?;
    let qe = parts.next()?.parse().ok()?;
    let ss = parts.next()?.parse().ok()?;
    let se = parts.next()?.parse().ok()?;
    Some((qs, qe, ss, se))
}

fn trace_target_score() -> Option<i32> {
    std::env::var("BLAST_RS_TRACE_SCORE").ok()?.parse().ok()
}

fn trace_lookup_probes() -> Option<Vec<usize>> {
    let raw = std::env::var("BLAST_RS_TRACE_LOOKUP_PROBES").ok()?;
    let probes = raw
        .split(',')
        .filter_map(|s| s.trim().parse::<usize>().ok())
        .collect::<Vec<_>>();
    if probes.is_empty() {
        None
    } else {
        Some(probes)
    }
}

fn trace_scan_positions() -> Option<Vec<usize>> {
    let raw = std::env::var("BLAST_RS_TRACE_SCAN_POSITIONS").ok()?;
    let positions = raw
        .split(',')
        .filter_map(|s| s.trim().parse::<usize>().ok())
        .collect::<Vec<_>>();
    if positions.is_empty() {
        None
    } else {
        Some(positions)
    }
}

fn trace_seed_pairs() -> Option<Vec<(usize, usize)>> {
    let raw = std::env::var("BLAST_RS_TRACE_SEED_PAIRS").ok()?;
    let pairs = raw
        .split(',')
        .filter_map(|entry| {
            let (q, s) = entry.trim().split_once(':')?;
            Some((q.parse::<usize>().ok()?, s.parse::<usize>().ok()?))
        })
        .collect::<Vec<_>>();
    if pairs.is_empty() {
        None
    } else {
        Some(pairs)
    }
}

fn matches_trace_seed_pair(q: usize, s: usize) -> bool {
    trace_seed_pairs()
        .map(|pairs| pairs.into_iter().any(|(qq, ss)| qq == q && ss == s))
        .unwrap_or(false)
}

fn hsp_matches_trace_target(hsp: &SearchHsp) -> bool {
    if let Some(score) = trace_target_score() {
        return hsp.score == score;
    }
    let Some((qs, qe, ss, se)) = trace_target_bounds() else {
        return false;
    };
    hsp.query_start <= qe && hsp.query_end >= qs && hsp.subject_start <= se && hsp.subject_end >= ss
}

fn trace_hsp(label: &str, hsp: &SearchHsp) {
    if !hsp_matches_trace_target(hsp) {
        return;
    }
    eprintln!(
        "[trace-target] {} hsp q={}..{} s={}..{} score={} gaps={} ident={} len={} ctx={}",
        label,
        hsp.query_start,
        hsp.query_end,
        hsp.subject_start,
        hsp.subject_end,
        hsp.score,
        hsp.gap_opens,
        hsp.num_ident,
        hsp.align_length,
        hsp.context
    );
}

fn candidate_matches_trace_target(candidate: &GappedCandidate) -> bool {
    if let Some(score) = trace_target_score() {
        if candidate_score(candidate) != score {
            return false;
        }
        return true;
    }
    let Some((qs, qe, ss, se)) = trace_target_bounds() else {
        return false;
    };
    let cq0 = candidate_query_start(candidate);
    let cq1 = candidate_query_end(candidate);
    let cs0 = candidate_subject_start(candidate);
    let cs1 = candidate_subject_end(candidate);
    cq0 <= qe && cq1 >= qs && cs0 <= se && cs1 >= ss
}

fn trace_candidate(label: &str, candidate: &GappedCandidate) {
    if !candidate_matches_trace_target(candidate) {
        return;
    }
    match candidate {
        GappedCandidate::Final(hsp) => {
            eprintln!(
                "[trace-target] {} final q={}..{} s={}..{} score={} gaps={} ident={} len={}",
                label,
                hsp.query_start,
                hsp.query_end,
                hsp.subject_start,
                hsp.subject_end,
                hsp.score,
                hsp.gap_opens,
                hsp.num_ident,
                hsp.align_length
            );
        }
        GappedCandidate::Traceback { tb, .. } => {
            let preview_len = tb.edit_script.ops.len().min(24);
            eprintln!(
                "[trace-target] {} tb q={}..{} s={}..{} score={} ops_len={} ops={:?}",
                label,
                tb.query_start,
                tb.query_end,
                tb.subject_start,
                tb.subject_end,
                tb.score,
                tb.edit_script.ops.len(),
                &tb.edit_script.ops[..preview_len]
            );
        }
    }
}

fn trace_traceback(label: &str, tb: &TracebackResult) {
    if let Some(score) = trace_target_score() {
        if tb.score != score {
            return;
        }
    } else {
        let Some((qs, qe, ss, se)) = trace_target_bounds() else {
            return;
        };
        let cq0 = tb.query_start as i32;
        let cq1 = tb.query_end as i32;
        let cs0 = tb.subject_start as i32;
        let cs1 = tb.subject_end as i32;
        if !(cq0 <= qe && cq1 >= qs && cs0 <= se && cs1 >= ss) {
            return;
        }
    }
    let preview_len = tb.edit_script.ops.len().min(24);
    eprintln!(
        "[trace-target] {} tb q={}..{} s={}..{} score={} ops_len={} ops={:?}",
        label,
        tb.query_start,
        tb.query_end,
        tb.subject_start,
        tb.subject_end,
        tb.score,
        tb.edit_script.ops.len(),
        &tb.edit_script.ops[..preview_len]
    );
}

fn cut_traceback_at(tb: &mut TracebackResult, q_cut: usize, s_cut: usize, cut_begin: bool) -> bool {
    let mut qid = 0usize;
    let mut sid = 0usize;
    let mut found = None;

    for (index, &(op_type, count_i32)) in tb.edit_script.ops.iter().enumerate() {
        let count = count_i32 as usize;
        let mut opid = 0usize;
        match op_type {
            GapAlignOpType::Sub => {
                while opid < count {
                    qid += 1;
                    sid += 1;
                    opid += 1;
                    if qid >= q_cut && sid >= s_cut {
                        found = Some((index, opid));
                        break;
                    }
                }
            }
            GapAlignOpType::Del | GapAlignOpType::Del1 | GapAlignOpType::Del2 => {
                sid += count;
                opid = count;
                if qid >= q_cut && sid >= s_cut {
                    found = Some((index, opid));
                }
            }
            GapAlignOpType::Ins | GapAlignOpType::Ins1 | GapAlignOpType::Ins2 => {
                qid += count;
                opid = count;
                if qid >= q_cut && sid >= s_cut {
                    found = Some((index, opid));
                }
            }
            GapAlignOpType::Decline => {
                qid += count;
                sid += count;
                opid = count;
                if qid >= q_cut && sid >= s_cut {
                    found = Some((index, opid));
                }
            }
        }
        if found.is_some() {
            break;
        }
    }

    let Some((index, opid)) = found else {
        return false;
    };

    let mut new_ops = Vec::new();
    if cut_begin {
        let count = tb.edit_script.ops[index].1 as usize;
        if opid < count {
            debug_assert!(matches!(tb.edit_script.ops[index].0, GapAlignOpType::Sub));
            new_ops.push((tb.edit_script.ops[index].0, (count - opid) as i32));
        }
        new_ops.extend(tb.edit_script.ops.iter().skip(index + 1).copied());
        tb.query_start += qid;
        tb.subject_start += sid;
    } else {
        new_ops.extend(tb.edit_script.ops.iter().take(index).copied());
        if opid < tb.edit_script.ops[index].1 as usize {
            debug_assert!(matches!(tb.edit_script.ops[index].0, GapAlignOpType::Sub));
            new_ops.push((tb.edit_script.ops[index].0, opid as i32));
        } else {
            new_ops.push(tb.edit_script.ops[index]);
        }
        tb.query_end = tb.query_start + qid;
        tb.subject_end = tb.subject_start + sid;
    }
    tb.edit_script = GapEditScript { ops: new_ops };
    true
}

fn purge_common_endpoint_tracebacks(candidates: &mut Vec<GappedCandidate>) {
    if candidates.len() <= 1 {
        return;
    }

    if trace_target_bounds().is_some() || trace_target_score().is_some() {
        eprintln!(
            "[trace-target] purge_common_endpoint_tracebacks candidates={}",
            candidates.len()
        );
        for candidate in candidates.iter() {
            trace_candidate("pre-purge", candidate);
        }
    }

    candidates.sort_by(|a, b| {
        candidate_context(a)
            .cmp(&candidate_context(b))
            .then_with(|| candidate_query_start(a).cmp(&candidate_query_start(b)))
            .then_with(|| candidate_subject_start(a).cmp(&candidate_subject_start(b)))
            .then_with(|| candidate_score(b).cmp(&candidate_score(a)))
            .then_with(|| candidate_query_end(b).cmp(&candidate_query_end(a)))
            .then_with(|| candidate_subject_end(b).cmp(&candidate_subject_end(a)))
    });
    let mut i = 0usize;
    while i < candidates.len() {
        let mut j = i + 1;
        while j < candidates.len()
            && candidate_context(&candidates[i]) == candidate_context(&candidates[j])
            && candidate_query_start(&candidates[i]) == candidate_query_start(&candidates[j])
            && candidate_subject_start(&candidates[i]) == candidate_subject_start(&candidates[j])
        {
            let leader_q_end = candidate_query_end(&candidates[i]);
            let leader_s_end = candidate_subject_end(&candidates[i]);
            let keep = match &mut candidates[j] {
                GappedCandidate::Traceback { tb, .. }
                    if tb.query_end as i32 > leader_q_end
                        && tb.subject_end as i32 > leader_s_end =>
                {
                    trace_traceback("cut-begin-before", tb);
                    cut_traceback_at(
                        tb,
                        (leader_q_end - tb.query_start as i32).max(0) as usize,
                        (leader_s_end - tb.subject_start as i32).max(0) as usize,
                        true,
                    ) && !tb.edit_script.ops.is_empty()
                }
                _ => false,
            };
            if keep {
                trace_candidate("cut-begin-after", &candidates[j]);
            }
            if keep {
                j += 1;
            } else {
                candidates.remove(j);
            }
        }
        i = j;
    }

    candidates.sort_by(|a, b| {
        candidate_context(a)
            .cmp(&candidate_context(b))
            .then_with(|| candidate_query_end(a).cmp(&candidate_query_end(b)))
            .then_with(|| candidate_subject_end(a).cmp(&candidate_subject_end(b)))
            .then_with(|| candidate_score(b).cmp(&candidate_score(a)))
            .then_with(|| candidate_query_start(b).cmp(&candidate_query_start(a)))
            .then_with(|| candidate_subject_start(b).cmp(&candidate_subject_start(a)))
    });
    let mut i = 0usize;
    while i < candidates.len() {
        let mut j = i + 1;
        while j < candidates.len()
            && candidate_context(&candidates[i]) == candidate_context(&candidates[j])
            && candidate_query_end(&candidates[i]) == candidate_query_end(&candidates[j])
            && candidate_subject_end(&candidates[i]) == candidate_subject_end(&candidates[j])
        {
            let leader_q_start = candidate_query_start(&candidates[i]);
            let leader_s_start = candidate_subject_start(&candidates[i]);
            let keep = match &mut candidates[j] {
                GappedCandidate::Traceback { tb, .. }
                    if (tb.query_start as i32) < leader_q_start
                        && (tb.subject_start as i32) < leader_s_start =>
                {
                    trace_traceback("cut-end-before", tb);
                    cut_traceback_at(
                        tb,
                        (leader_q_start - tb.query_start as i32).max(0) as usize,
                        (leader_s_start - tb.subject_start as i32).max(0) as usize,
                        false,
                    ) && !tb.edit_script.ops.is_empty()
                }
                _ => false,
            };
            if keep {
                trace_candidate("cut-end-after", &candidates[j]);
            }
            if keep {
                j += 1;
            } else {
                candidates.remove(j);
            }
        }
        i = j;
    }
}

fn finalize_gapped_candidates(
    mut candidates: Vec<GappedCandidate>,
    query_plus_nomask: &[u8],
    query_minus_nomask: &[u8],
    subject: &[u8],
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    kbp: &KarlinBlk,
    search_space: f64,
    evalue_threshold: f64,
    min_diag_separation: i32,
) -> Vec<SearchHsp> {
    collapse_exact_traceback_duplicates(
        &mut candidates,
        query_plus_nomask,
        query_minus_nomask,
        subject,
    );
    trace_pre_finalize_tracebacks(&candidates, query_plus_nomask, query_minus_nomask, subject);
    purge_common_endpoint_tracebacks(&mut candidates);

    let mut hsps = Vec::new();
    for candidate in candidates {
        match candidate {
            GappedCandidate::Final(hsp) => hsps.push(hsp),
            GappedCandidate::Traceback { context, tb } => {
                if let Some(hsp) = render_traceback_candidate(
                    context,
                    tb,
                    query_plus_nomask,
                    query_minus_nomask,
                    subject,
                    reward,
                    penalty,
                    gap_open,
                    gap_extend,
                    kbp,
                    search_space,
                    evalue_threshold,
                ) {
                    hsps.push(hsp);
                }
            }
        }
    }
    finalize_rendered_hsps(hsps, min_diag_separation)
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

const MAX_SUBJECT_OFFSET: usize = 90_000;
const MAX_TOTAL_GAPS: usize = 3_000;

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

fn adjust_subject_range(
    subject_offset: &mut usize,
    subject_length: &mut usize,
    query_offset: usize,
    query_length: usize,
) -> usize {
    if *subject_length < MAX_SUBJECT_OFFSET {
        return 0;
    }

    let s_offset = *subject_offset;
    let max_extension_left = query_offset.saturating_add(MAX_TOTAL_GAPS);
    let max_extension_right = query_length
        .saturating_sub(query_offset)
        .saturating_add(MAX_TOTAL_GAPS);

    let start_shift = if s_offset <= max_extension_left {
        0
    } else {
        *subject_offset = max_extension_left;
        s_offset - max_extension_left
    };

    *subject_length = (*subject_length).min(s_offset.saturating_add(max_extension_right)) - start_shift;
    start_shift
}

/// Port of NCBI `BlastGetStartForGappedAlignmentNucl`
/// (`blast_gapalign.c:3323`). Given an initial diagonal seed, tighten it onto
/// a longer exact-match run along the same diagonal when possible.
fn blast_get_start_for_gapped_alignment_nucl(
    query: &[u8],
    subject: &[u8],
    hsp_q_start: usize,
    hsp_q_end: usize,
    hsp_s_start: usize,
    hsp_s_end: usize,
    seed_q: usize,
    seed_s: usize,
) -> (usize, usize) {
    let trace_this = trace_target_bounds()
        .map(|(qs, qe, ss, se)| {
            hsp_q_start as i32 <= qe
                && hsp_q_end as i32 >= qs
                && hsp_s_start as i32 <= se
                && hsp_s_end as i32 >= ss
        })
        .unwrap_or(false);
    let mut hsp_max_ident_run = 10usize;
    let offset = (seed_s.saturating_sub(hsp_s_start)).min(seed_q.saturating_sub(hsp_q_start));

    // First check whether the existing seed already sits inside a sufficiently
    // long exact run.
    let mut score = -1i32;
    let mut q = seed_q;
    let mut s = seed_s;
    while q < hsp_q_end && s < subject.len() && query[q] == subject[s] {
        score += 1;
        if score > hsp_max_ident_run as i32 {
            if trace_this {
                eprintln!(
                    "[trace-target] tbseed early-keep hsp_q={}..{} hsp_s={}..{} seed=({}, {}) score={} threshold={}",
                    hsp_q_start,
                    hsp_q_end,
                    hsp_s_start,
                    hsp_s_end,
                    seed_q,
                    seed_s,
                    score,
                    hsp_max_ident_run
                );
            }
            return (seed_q, seed_s);
        }
        q += 1;
        s += 1;
    }
    let mut q_back = seed_q as isize;
    let mut s_back = seed_s as isize;
    while q_back >= 0 && s_back >= 0 && query[q_back as usize] == subject[s_back as usize] {
        score += 1;
        if score > hsp_max_ident_run as i32 {
            if trace_this {
                eprintln!(
                    "[trace-target] tbseed early-keep hsp_q={}..{} hsp_s={}..{} seed=({}, {}) back_score={} threshold={}",
                    hsp_q_start,
                    hsp_q_end,
                    hsp_s_start,
                    hsp_s_end,
                    seed_q,
                    seed_s,
                    score,
                    hsp_max_ident_run
                );
            }
            return (seed_q, seed_s);
        }
        q_back -= 1;
        s_back -= 1;
    }

    // If we move the start point, require a longer exact run.
    hsp_max_ident_run = ((hsp_max_ident_run as f32) * 1.5) as usize;
    let q_start = seed_q.saturating_sub(offset);
    let s_start = seed_s.saturating_sub(offset);
    let q_len = (hsp_s_end.saturating_sub(s_start)).min(hsp_q_end.saturating_sub(q_start));
    if q_len == 0 {
        return (seed_q, seed_s);
    }

    let mut max_score = 0usize;
    let mut max_offset = q_start;
    let mut run_score = 0usize;
    let mut is_match = false;
    let mut prev_match = false;
    let mut last_index = q_start;

    for index in q_start..(q_start + q_len) {
        last_index = index;
        is_match = query[index] == subject[s_start + (index - q_start)];
        if is_match != prev_match {
            prev_match = is_match;
            if is_match {
                run_score = 1;
            } else if run_score > max_score {
                max_score = run_score;
                max_offset = index - run_score / 2;
            }
        } else if is_match {
            run_score += 1;
            if run_score > hsp_max_ident_run {
                let q_seed = index - hsp_max_ident_run / 2;
                if trace_this {
                    eprintln!(
                        "[trace-target] tbseed move-early hsp_q={}..{} hsp_s={}..{} seed=({}, {}) q_start={} s_start={} q_len={} index={} run_score={} threshold={} result=({}, {})",
                        hsp_q_start,
                        hsp_q_end,
                        hsp_s_start,
                        hsp_s_end,
                        seed_q,
                        seed_s,
                        q_start,
                        s_start,
                        q_len,
                        index,
                        run_score,
                        hsp_max_ident_run,
                        q_seed,
                        q_seed + s_start - q_start
                    );
                }
                return (q_seed, q_seed + s_start - q_start);
            }
        }
    }

    if is_match && run_score > max_score {
        max_score = run_score;
        max_offset = last_index + 1 - run_score / 2;
    }

    if max_score > 0 {
        if trace_this {
            eprintln!(
                "[trace-target] tbseed move-final hsp_q={}..{} hsp_s={}..{} seed=({}, {}) q_start={} s_start={} q_len={} max_score={} max_offset={} result=({}, {})",
                hsp_q_start,
                hsp_q_end,
                hsp_s_start,
                hsp_s_end,
                seed_q,
                seed_s,
                q_start,
                s_start,
                q_len,
                max_score,
                max_offset,
                max_offset,
                max_offset + s_start - q_start
            );
        }
        (max_offset, max_offset + s_start - q_start)
    } else {
        if trace_this {
            eprintln!(
                "[trace-target] tbseed no-move hsp_q={}..{} hsp_s={}..{} seed=({}, {}) q_start={} s_start={} q_len={}",
                hsp_q_start,
                hsp_q_end,
                hsp_s_start,
                hsp_s_end,
                seed_q,
                seed_s,
                q_start,
                s_start,
                q_len
            );
        }
        (seed_q, seed_s)
    }
}

#[inline]
fn adjust_nt_prelim_seed_like_ncbi(
    hsp_q_end: usize,
    hsp_s_end: usize,
    seed_q: usize,
    seed_s: usize,
    query_len: usize,
    subject_len: usize,
) -> (usize, usize) {
    if hsp_q_end >= seed_q.saturating_add(8)
        && hsp_s_end >= seed_s.saturating_add(8)
        && seed_q.saturating_add(3) < query_len
        && seed_s.saturating_add(3) < subject_len
    {
        (seed_q + 3, seed_s + 3)
    } else {
        (seed_q, seed_s)
    }
}

#[inline]
fn adjust_nt_packed_prelim_start_like_ncbi(
    seed_q: usize,
    seed_s: usize,
    query_len: usize,
    subject_len: usize,
) -> (usize, usize) {
    let offset_adjustment = 4usize.saturating_sub(seed_s % 4);
    let mut adjusted_q = seed_q + offset_adjustment.saturating_sub(1);
    let mut adjusted_s = seed_s + offset_adjustment.saturating_sub(1);
    if adjusted_q >= query_len || adjusted_s >= subject_len {
        adjusted_q = adjusted_q.saturating_sub(4);
        adjusted_s = adjusted_s.saturating_sub(4);
    }
    (
        adjusted_q.min(query_len.saturating_sub(1)),
        adjusted_s.min(subject_len.saturating_sub(1)),
    )
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

fn score_edit_script_affine(
    ops: &[(GapAlignOpType, i32)],
    query: &[u8],
    subject: &[u8],
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
) -> i32 {
    let mut qp = 0usize;
    let mut sp = 0usize;
    let mut score = 0i32;
    for &(op, count_i32) in ops {
        let count = count_i32.max(0) as usize;
        match op {
            GapAlignOpType::Sub => {
                for _ in 0..count {
                    let q = query.get(qp).copied().unwrap_or(15);
                    let s = subject.get(sp).copied().unwrap_or(15);
                    score += if q == s { reward } else { penalty };
                    qp += 1;
                    sp += 1;
                }
            }
            GapAlignOpType::Del | GapAlignOpType::Del1 | GapAlignOpType::Del2 => {
                if count > 0 {
                    score -= gap_open + gap_extend * count as i32;
                    sp += count;
                }
            }
            GapAlignOpType::Ins | GapAlignOpType::Ins1 | GapAlignOpType::Ins2 => {
                if count > 0 {
                    score -= gap_open + gap_extend * count as i32;
                    qp += count;
                }
            }
            GapAlignOpType::Decline => {
                qp += count;
                sp += count;
            }
        }
    }
    score
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
    blastn_gapped_search_nomask_with_split_xdrop(
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
        x_dropoff,
        x_dropoff,
        x_dropoff_final.max(x_dropoff),
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
    blastn_gapped_search_nomask_with_split_xdrop(
        query_plus,
        query_minus,
        query_plus_nomask,
        query_minus_nomask,
        subject,
        word_size,
        reward,
        penalty,
        _gap_open,
        _gap_extend,
        x_dropoff,
        x_dropoff,
        x_dropoff,
        kbp,
        search_space,
        evalue_threshold,
    )
}

fn blastn_gapped_search_nomask_with_split_xdrop(
    query_plus: &[u8],
    query_minus: &[u8],
    query_plus_nomask: &[u8],
    query_minus_nomask: &[u8],
    subject: &[u8],
    word_size: usize,
    reward: i32,
    penalty: i32,
    _gap_open: i32,
    _gap_extend: i32,
    ungapped_x_dropoff: i32,
    prelim_x_dropoff: i32,
    traceback_x_dropoff: i32,
    kbp: &KarlinBlk,
    search_space: f64,
    evalue_threshold: f64,
) -> Vec<SearchHsp> {
    let prepared = PreparedBlastnQuery::new(query_plus, query_minus, word_size);
    // First do ungapped search to find seeds (uses masked query)
    let ungapped = blastn_ungapped_search(
        query_plus,
        query_minus,
        subject,
        word_size,
        reward,
        penalty,
        ungapped_x_dropoff,
        kbp,
        search_space,
        evalue_threshold * 100.0, // permissive threshold for seeds
    );

    let mut candidates = Vec::new();
    let mut perfect_seeds = Vec::new();
    let min_diag_separation =
        min_diag_separation_for_gapped(word_size, reward, penalty, _gap_open, _gap_extend);

    collect_decoded_gapped_candidates(
        &prepared,
        &ungapped,
        &mut candidates,
        &mut perfect_seeds,
        query_plus_nomask,
        query_minus_nomask,
        subject,
        reward,
        penalty,
        _gap_open,
        _gap_extend,
        prelim_x_dropoff,
        traceback_x_dropoff,
        evalue_threshold,
        search_space,
        kbp,
        min_diag_separation,
    );
    if min_diag_separation != 0 && !perfect_seeds.is_empty() {
        prune_perfect_ungapped_gapped_candidates(&mut perfect_seeds, min_diag_separation);
    }
    candidates.extend(perfect_seeds.into_iter().map(GappedCandidate::Final));

    finalize_gapped_candidates(
        candidates,
        query_plus_nomask,
        query_minus_nomask,
        subject,
        reward,
        penalty,
        _gap_open,
        _gap_extend,
        kbp,
        search_space,
        evalue_threshold,
        min_diag_separation,
    )
}

#[allow(clippy::too_many_arguments)]
fn collect_decoded_gapped_candidates(
    prepared: &PreparedBlastnQuery<'_>,
    ungapped: &[SearchHsp],
    candidates: &mut Vec<GappedCandidate>,
    perfect_seeds: &mut Vec<SearchHsp>,
    query_plus_nomask: &[u8],
    query_minus_nomask: &[u8],
    subject: &[u8],
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    prelim_x_dropoff: i32,
    traceback_x_dropoff: i32,
    _evalue_threshold: f64,
    _search_space: f64,
    _kbp: &KarlinBlk,
    min_diag_separation: i32,
) {
    for seed in ungapped {
        let traceback_query = if seed.context == 0 {
            query_plus_nomask
        } else {
            query_minus_nomask
        };
        let prelim_query = prepared
            .lookups
            .iter()
            .find(|lookup| lookup.context == seed.context)
            .map(|lookup| lookup.query)
            .unwrap_or(traceback_query);

        if min_diag_separation != 0 && is_perfect_ungapped_hsp(seed, reward) {
            perfect_seeds.push(seed.clone());
            continue;
        }

        let (prelim_seed, seed_q, seed_s) = select_decoded_candidate_seed(
            prelim_query,
            subject,
            seed,
            reward,
            penalty,
            gap_open == 0 && gap_extend == 0,
        );
        if candidate_matches_trace_target(&GappedCandidate::Final(seed.clone())) {
            eprintln!(
                "[trace-target] seed q={} s={} hsp_q={}..{} hsp_s={}..{} score={}",
                seed_q,
                seed_s,
                seed.query_start,
                seed.query_end,
                seed.subject_start,
                seed.subject_end,
                seed.score
            );
        }
        let Some(prelim) = preliminary_gapped_score_and_seed(
            prelim_query,
            subject,
            seed,
            prelim_seed.query,
            prelim_seed.subject,
            seed.score,
            reward,
            penalty,
            gap_open,
            gap_extend,
            prelim_x_dropoff,
        )
        else {
            continue;
        };
        let traceback_seed = traceback_seed_from_preliminary_decoded(
            traceback_query,
            subject,
            seed,
            prelim,
            seed_q,
            seed_s,
        );
        let prelim_seed_hsp = prelim_hsp_from_result(seed, prelim);
        if let Some(tb) = run_decoded_traceback_for_seed(
            traceback_query,
            subject,
            &prelim_seed_hsp,
            traceback_seed,
            reward,
            penalty,
            gap_open,
            gap_extend,
            traceback_x_dropoff,
        ) {
            candidates.push(GappedCandidate::Traceback {
                context: prelim_seed_hsp.context,
                tb,
            });
        } else {
            candidates.push(GappedCandidate::Final(seed.clone()));
        }
    }
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
    ungapped_x_dropoff: i32,
    x_dropoff: i32,
    x_dropoff_final: i32,
    kbp: &KarlinBlk,
    search_space: f64,
    evalue_threshold: f64,
    last_hit_scratch: &mut [PackedDiagScratch],
) -> Vec<SearchHsp> {
    blastn_gapped_search_packed_prepared_with_split_xdrop(
        prepared,
        query_plus_nomask,
        query_minus_nomask,
        subject_packed,
        subject_len,
        reward,
        penalty,
        gap_open,
        gap_extend,
        ungapped_x_dropoff,
        x_dropoff,
        x_dropoff_final.max(x_dropoff),
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
    ungapped_x_dropoff: i32,
    x_dropoff: i32,
    kbp: &KarlinBlk,
    search_space: f64,
    evalue_threshold: f64,
    last_hit_scratch: &mut [PackedDiagScratch],
) -> Vec<SearchHsp> {
    blastn_gapped_search_packed_prepared_with_split_xdrop(
        prepared,
        query_plus_nomask,
        query_minus_nomask,
        subject_packed,
        subject_len,
        reward,
        penalty,
        gap_open,
        gap_extend,
        ungapped_x_dropoff,
        x_dropoff,
        x_dropoff,
        kbp,
        search_space,
        evalue_threshold,
        last_hit_scratch,
    )
}

fn trace_packed_prepared_search_path(prepared: &PreparedBlastnQuery<'_>) {
    if trace_target_bounds().is_some() || trace_target_score().is_some() {
        let lookup_desc = prepared
            .lookups
            .iter()
            .map(|l| format!("ctx{}:lut{} start{} step{}", l.context, l.lut_word, l.scan_start, l.scan_step))
            .collect::<Vec<_>>()
            .join(" ");
        eprintln!(
            "[trace-target] path=blastn_gapped_search_packed_prepared word_size={} force_decoded={} force_naive={} {}",
            prepared.word_size,
            std::env::var_os("BLAST_RS_FORCE_DECODED_UNGAPPED").is_some(),
            std::env::var_os("BLAST_RS_FORCE_NAIVE_UNGAPPED").is_some(),
            lookup_desc
        );
    }
}

fn collect_packed_prepared_ungapped(
    prepared: &PreparedBlastnQuery<'_>,
    subject_packed: &[u8],
    subject_len: usize,
    reward: i32,
    penalty: i32,
    prelim_x_dropoff: i32,
    kbp: &KarlinBlk,
    search_space: f64,
    evalue_threshold: f64,
    last_hit_scratch: &mut [PackedDiagScratch],
) -> Vec<SearchHsp> {
    let ungapped_collect_evalue = if std::env::var_os("BLAST_RS_WIDE_UNGAPPED_EVALUE").is_some() {
        evalue_threshold * 10_000.0
    } else {
        evalue_threshold * 100.0
    };
    if std::env::var_os("BLAST_RS_FORCE_DECODED_UNGAPPED").is_some() {
        let subject_decoded = decode_packed_ncbi2na(subject_packed, subject_len);
        let mut decoded_last_hit_scratch: Vec<Vec<i32>> = prepared
            .lookups
            .iter()
            .map(|lookup| vec![0; lookup.diag_array_len])
            .collect();
        blastn_ungapped_search_decoded_prepared_with_scratch_no_dedup(
            prepared,
            &subject_decoded,
            reward,
            penalty,
            prelim_x_dropoff,
            kbp,
            search_space,
            ungapped_collect_evalue,
            &mut decoded_last_hit_scratch,
        )
    } else if std::env::var_os("BLAST_RS_SKIP_UNGAPPED_DEDUP").is_some() {
        blastn_ungapped_search_packed_prepared_with_scratch_no_dedup(
            prepared,
            subject_packed,
            subject_len,
            reward,
            penalty,
            prelim_x_dropoff,
            kbp,
            search_space,
            ungapped_collect_evalue,
            last_hit_scratch,
        )
    } else {
        blastn_ungapped_search_packed_prepared_with_scratch(
            prepared,
            subject_packed,
            subject_len,
            reward,
            penalty,
            prelim_x_dropoff,
            kbp,
            search_space,
            ungapped_collect_evalue,
            last_hit_scratch,
        )
    }
}

#[allow(clippy::too_many_arguments)]
fn collect_packed_gapped_candidates(
    prepared: &PreparedBlastnQuery<'_>,
    ungapped: &[SearchHsp],
    candidates: &mut Vec<GappedCandidate>,
    perfect_seeds: &mut Vec<SearchHsp>,
    query_plus_nomask: &[u8],
    query_minus_nomask: &[u8],
    subject_packed: &[u8],
    subject_len: usize,
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    prelim_x_dropoff: i32,
    traceback_x_dropoff: i32,
    cutoff_score: i32,
    evalue_threshold: f64,
    search_space: f64,
    kbp: &KarlinBlk,
    min_diag_separation: i32,
    subject_decoded: &mut Option<Vec<u8>>,
    profile_enabled: bool,
    stats: &mut PackedGappedStageStats,
) {
    let filter_by_cutoff = std::env::var_os("BLAST_RS_SKIP_GAPPED_SEED_CUTOFF").is_none();
    for seed in ungapped
        .iter()
        .filter(|s| !filter_by_cutoff || s.score >= cutoff_score)
    {
        stats.filtered_seed_count += 1;
        let traceback_query = if seed.context == 0 {
            query_plus_nomask
        } else {
            query_minus_nomask
        };
        let prelim_query = prepared
            .lookups
            .iter()
            .find(|lookup| lookup.context == seed.context)
            .map(|lookup| lookup.query)
            .unwrap_or(traceback_query);

        if min_diag_separation != 0 && is_perfect_ungapped_hsp(seed, reward) {
            perfect_seeds.push(seed.clone());
            continue;
        }

        let subject_decoded =
            subject_decoded.get_or_insert_with(|| decode_packed_ncbi2na(subject_packed, subject_len));
        let (prelim_seed, seed_q, seed_s) = select_packed_candidate_seed(
            prelim_query,
            subject_decoded,
            subject_len,
            seed,
            reward,
            penalty,
            gap_open == 0 && gap_extend == 0,
        );
        if candidate_matches_trace_target(&GappedCandidate::Final(seed.clone())) {
            eprintln!(
                "[trace-target] seed q={} s={} hsp_q={}..{} hsp_s={}..{} score={}",
                seed_q,
                seed_s,
                seed.query_start,
                seed.query_end,
                seed.subject_start,
                seed.subject_end,
                seed.score
            );
        }

        let prelim_start = if profile_enabled {
            Some(Instant::now())
        } else {
            None
        };
        let Some(prelim) = (if std::env::var_os("BLAST_RS_USE_DECODED_PRELIM").is_some() {
            preliminary_gapped_score_and_seed(
                prelim_query,
                subject_decoded,
                seed,
                prelim_seed.query,
                prelim_seed.subject,
                seed.score,
                reward,
                penalty,
                gap_open,
                gap_extend,
                prelim_x_dropoff,
            )
        } else {
            preliminary_gapped_score_and_seed_packed(
                prelim_query,
                subject_packed,
                subject_len,
                prelim_seed.query,
                prelim_seed.subject,
                seed.score,
                reward,
                penalty,
                gap_open,
                gap_extend,
                prelim_x_dropoff,
            )
        }) else {
            if let Some(start) = prelim_start {
                stats.preliminary_score_ms += start.elapsed().as_millis();
            }
            continue;
        };
        if candidate_matches_trace_target(&GappedCandidate::Final(seed.clone())) {
            eprintln!(
                "[trace-target] prelim-raw score={} prelim_hsp q={}..{} s={}..{} raw_seed=({}, {}) gapped_start=({}, {})",
                prelim.score,
                prelim.prelim_q_start,
                prelim.prelim_q_end,
                prelim.prelim_s_start,
                prelim.prelim_s_end,
                seed_q,
                seed_s,
                prelim.gapped_start_q,
                prelim.gapped_start_s
            );
        }
        if std::env::var_os("BLAST_RS_SKIP_PRELIM_EVALUE_CUTOFF").is_none()
            && std::env::var_os("BLAST_RS_ENABLE_PRELIM_EVALUE_CUTOFF").is_some()
            && kbp.raw_to_evalue(prelim.score, search_space) > evalue_threshold
        {
            if let Some(start) = prelim_start {
                stats.preliminary_score_ms += start.elapsed().as_millis();
            }
            continue;
        }
        let traceback_seed = traceback_seed_from_preliminary(
            traceback_query,
            subject_decoded,
            seed,
            prelim,
            seed_q,
            seed_s,
        );
        if let Some(start) = prelim_start {
            stats.preliminary_score_ms += start.elapsed().as_millis();
        }
        stats.preliminary_pass_count += 1;

        let prelim_seed_hsp = prelim_hsp_from_result(seed, prelim);

        let traceback_start = if profile_enabled {
            Some(Instant::now())
        } else {
            None
        };
        if let Some(tb) = run_packed_traceback_for_seed(
            traceback_query,
            subject_decoded,
            seed,
            &prelim_seed_hsp,
            traceback_seed,
            reward,
            penalty,
            gap_open,
            gap_extend,
            traceback_x_dropoff,
        ) {
            candidates.push(GappedCandidate::Traceback {
                context: prelim_seed_hsp.context,
                tb,
            });
            if let Some(start) = traceback_start {
                stats.traceback_ms += start.elapsed().as_millis();
                stats.traceback_count += 1;
            }
        } else {
            candidates.push(GappedCandidate::Final(seed.clone()));
            if let Some(start) = traceback_start {
                stats.traceback_ms += start.elapsed().as_millis();
            }
        }
    }
}

fn select_packed_candidate_seed(
    query: &[u8],
    subject_decoded: &[u8],
    subject_len: usize,
    seed: &SearchHsp,
    reward: i32,
    penalty: i32,
    greedy_prelim: bool,
) -> (PreliminarySeed, usize, usize) {
    let (seed_q, seed_s) = if greedy_prelim {
        (
            ((seed.query_start + seed.query_end) / 2) as usize,
            ((seed.subject_start + seed.subject_end) / 2) as usize,
        )
    } else {
        blast_get_offsets_for_gapped_alignment(
            query,
            subject_decoded,
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
        })
    };
    let (prelim_q, prelim_s) = if greedy_prelim {
        (seed_q, seed_s)
    } else {
        adjust_nt_prelim_seed_like_ncbi(
            seed.query_end as usize,
            seed.subject_end as usize,
            seed_q,
            seed_s,
            query.len(),
            subject_len,
        )
    };
    (
        PreliminarySeed {
            query: prelim_q,
            subject: prelim_s,
        },
        seed_q,
        seed_s,
    )
}

fn prelim_hsp_from_result(seed: &SearchHsp, prelim: PreliminaryGappedResult) -> SearchHsp {
    SearchHsp {
        query_start: prelim.prelim_q_start as i32,
        query_end: prelim.prelim_q_end as i32,
        subject_start: prelim.prelim_s_start as i32,
        subject_end: prelim.prelim_s_end as i32,
        score: prelim.score,
        bit_score: 0.0,
        evalue: 0.0,
        num_ident: 0,
        align_length: (prelim.prelim_q_end.saturating_sub(prelim.prelim_q_start)) as i32,
        mismatches: 0,
        gap_opens: 0,
        context: seed.context,
        qseq: None,
        sseq: None,
    }
}

fn select_decoded_candidate_seed(
    query: &[u8],
    subject: &[u8],
    seed: &SearchHsp,
    reward: i32,
    penalty: i32,
    greedy_prelim: bool,
) -> (PreliminarySeed, usize, usize) {
    let (seed_q, seed_s) = if greedy_prelim {
        (
            ((seed.query_start + seed.query_end) / 2) as usize,
            ((seed.subject_start + seed.subject_end) / 2) as usize,
        )
    } else {
        blast_get_offsets_for_gapped_alignment(
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
        })
    };
    let (prelim_q, prelim_s) = if greedy_prelim {
        (seed_q, seed_s)
    } else {
        adjust_nt_prelim_seed_like_ncbi(
            seed.query_end as usize,
            seed.subject_end as usize,
            seed_q,
            seed_s,
            query.len(),
            subject.len(),
        )
    };
    (
        PreliminarySeed {
            query: prelim_q,
            subject: prelim_s,
        },
        seed_q,
        seed_s,
    )
}

fn trace_preliminary_decoded_candidate(
    seed: &SearchHsp,
    prelim: PreliminaryGappedResult,
    trace_seed: PackedTracebackSeed,
) {
    if !candidate_matches_trace_target(&GappedCandidate::Final(seed.clone())) {
        return;
    }
    eprintln!(
        "[trace-target] prelim score={} prelim_hsp q={}..{} s={}..{} raw_seed=({}, {}) gapped_start=({}, {}) traceback_seed=({}, {})",
        prelim.score,
        prelim.prelim_q_start,
        prelim.prelim_q_end,
        prelim.prelim_s_start,
        prelim.prelim_s_end,
        trace_seed.raw_query,
        trace_seed.raw_subject,
        prelim.gapped_start_q,
        prelim.gapped_start_s,
        trace_seed.traceback_query,
        trace_seed.traceback_subject
    );
}

fn traceback_seed_from_preliminary_decoded(
    query: &[u8],
    subject: &[u8],
    seed: &SearchHsp,
    prelim: PreliminaryGappedResult,
    raw_seed_q: usize,
    raw_seed_s: usize,
) -> PackedTracebackSeed {
    let (traceback_query, traceback_subject) = blast_get_start_for_gapped_alignment_nucl(
        query,
        subject,
        prelim.prelim_q_start,
        prelim.prelim_q_end,
        prelim.prelim_s_start,
        prelim.prelim_s_end,
        prelim.gapped_start_q,
        prelim.gapped_start_s,
    );
    let trace_seed = PackedTracebackSeed {
        raw_query: raw_seed_q,
        raw_subject: raw_seed_s,
        traceback_query,
        traceback_subject,
    };
    trace_preliminary_decoded_candidate(seed, prelim, trace_seed);
    trace_seed
}

fn run_decoded_traceback_for_seed(
    query: &[u8],
    subject: &[u8],
    prelim_seed_hsp: &SearchHsp,
    traceback_seed: PackedTracebackSeed,
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    traceback_x_dropoff: i32,
) -> Option<TracebackResult> {
    blast_gapped_align_windowed(
        query,
        subject,
        prelim_seed_hsp,
        traceback_seed.traceback_query,
        traceback_seed.traceback_subject,
        reward,
        penalty,
        gap_open,
        gap_extend,
        traceback_x_dropoff,
    )
}

fn trace_preliminary_packed_candidate(
    seed: &SearchHsp,
    prelim: PreliminaryGappedResult,
    trace_seed: PackedTracebackSeed,
) {
    if !candidate_matches_trace_target(&GappedCandidate::Final(seed.clone())) {
        return;
    }
    eprintln!(
        "[trace-target] prelim score={} prelim_hsp q={}..{} s={}..{} raw_seed=({}, {}) gapped_start=({}, {}) traceback_seed=({}, {}) qlen={} slen={} private=({}, {}) lr_scores=({}, {})",
        prelim.score,
        prelim.prelim_q_start,
        prelim.prelim_q_end,
        prelim.prelim_s_start,
        prelim.prelim_s_end,
        trace_seed.raw_query,
        trace_seed.raw_subject,
        prelim.gapped_start_q,
        prelim.gapped_start_s,
        trace_seed.traceback_query,
        trace_seed.traceback_subject,
        prelim.q_length,
        prelim.s_length,
        prelim.private_q_start,
        prelim.private_s_start,
        prelim.score_left,
        prelim.score_right
    );
}

fn traceback_seed_from_preliminary(
    query: &[u8],
    subject: &[u8],
    seed: &SearchHsp,
    prelim: PreliminaryGappedResult,
    raw_seed_q: usize,
    raw_seed_s: usize,
) -> PackedTracebackSeed {
    let (traceback_query, traceback_subject) = blast_get_start_for_gapped_alignment_nucl(
        query,
        subject,
        prelim.prelim_q_start,
        prelim.prelim_q_end,
        prelim.prelim_s_start,
        prelim.prelim_s_end,
        prelim.gapped_start_q,
        prelim.gapped_start_s,
    );
    trace_preliminary_packed_candidate(
        seed,
        prelim,
        PackedTracebackSeed {
            raw_query: raw_seed_q,
            raw_subject: raw_seed_s,
            traceback_query,
            traceback_subject,
        },
    );
    PackedTracebackSeed {
        raw_query: raw_seed_q,
        raw_subject: raw_seed_s,
        traceback_query,
        traceback_subject,
    }
}

fn run_packed_traceback_for_seed(
    query: &[u8],
    subject: &[u8],
    seed: &SearchHsp,
    prelim_seed_hsp: &SearchHsp,
    traceback_seed: PackedTracebackSeed,
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    traceback_x_dropoff: i32,
) -> Option<TracebackResult> {
    blast_gapped_align_windowed(
        query,
        subject,
        prelim_seed_hsp,
        traceback_seed.traceback_query,
        traceback_seed.traceback_subject,
        reward,
        penalty,
        gap_open,
        gap_extend,
        traceback_x_dropoff,
    )
    .or_else(|| {
        if candidate_matches_trace_target(&GappedCandidate::Final(seed.clone())) {
            eprintln!(
                "[trace-target] traceback-fallback q={}..{} s={}..{} score={}",
                seed.query_start, seed.query_end, seed.subject_start, seed.subject_end, seed.score
            );
        }
        None
    })
}

fn trace_pre_finalize_tracebacks(
    candidates: &[GappedCandidate],
    query_plus_nomask: &[u8],
    query_minus_nomask: &[u8],
    subject: &[u8],
) {
    if trace_target_bounds().is_none() {
        return;
    }
    for candidate in candidates {
        if !candidate_matches_trace_target(candidate) {
            continue;
        }
        if let GappedCandidate::Traceback { context, tb } = candidate {
            let query = candidate_query_for_context(*context, query_plus_nomask, query_minus_nomask);
            let q_slice = &query[tb.query_start..tb.query_end];
            let s_slice = &subject[tb.subject_start..tb.subject_end];
            let (align_len, num_ident, gap_opens) = tb.edit_script.count_identities(q_slice, s_slice);
            eprintln!(
                "[trace-target] pre-finalize tb q={}..{} s={}..{} score={} ident={} len={} gaps={}",
                tb.query_start,
                tb.query_end,
                tb.subject_start,
                tb.subject_end,
                tb.score,
                num_ident,
                align_len,
                gap_opens
            );
        }
    }
}

fn render_traceback_candidate(
    context: i32,
    mut tb: TracebackResult,
    query_plus_nomask: &[u8],
    query_minus_nomask: &[u8],
    subject: &[u8],
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    kbp: &KarlinBlk,
    search_space: f64,
    evalue_threshold: f64,
) -> Option<SearchHsp> {
    let query = candidate_query_for_context(context, query_plus_nomask, query_minus_nomask);
    let cutoff_score = kbp.evalue_to_raw(evalue_threshold, search_space);
    if reevaluate_with_ambiguities_gapped(
        &mut tb,
        query,
        subject,
        reward,
        penalty,
        gap_open,
        gap_extend,
        cutoff_score.max(1),
    ) {
        return None;
    }
    let evalue = kbp.raw_to_evalue(tb.score, search_space);
    if evalue > evalue_threshold {
        return None;
    }
    let q_slice = &query[tb.query_start..tb.query_end];
    let s_slice = &subject[tb.subject_start..tb.subject_end];
    let (align_len, num_ident, gap_opens) = tb.edit_script.count_identities(q_slice, s_slice);
    if std::env::var_os("BLAST_RS_PROFILE").is_some() && gap_opens >= 100 {
        let script_score = score_edit_script_affine(
            &tb.edit_script.ops,
            q_slice,
            s_slice,
            reward,
            penalty,
            gap_open,
            gap_extend,
        );
        let preview_len = tb.edit_script.ops.len().min(24);
        eprintln!(
            "[blastn-profile] fragmented_tb score={} script_score={} q={}..{} s={}..{} align_len={} ident={} gap_opens={} ops_preview={:?}",
            tb.score,
            script_score,
            tb.query_start,
            tb.query_end,
            tb.subject_start,
            tb.subject_end,
            align_len,
            num_ident,
            gap_opens,
            &tb.edit_script.ops[..preview_len]
        );
    }
    let (qseq, sseq) = tb.edit_script.render_alignment(q_slice, s_slice, blastna_to_iupac);

    let hsp = SearchHsp {
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
        context,
        qseq: Some(qseq),
        sseq: Some(sseq),
    };
    if let Some((qs, qe, ss, se)) = trace_target_bounds() {
        if hsp.query_start <= qe
            && hsp.query_end >= qs
            && hsp.subject_start <= se
            && hsp.subject_end >= ss
        {
            eprintln!(
                "[trace-target] final-hsp q={}..{} s={}..{} score={} gaps={} ident={} len={}",
                hsp.query_start,
                hsp.query_end,
                hsp.subject_start,
                hsp.subject_end,
                hsp.score,
                hsp.gap_opens,
                hsp.num_ident,
                hsp.align_length
            );
        }
    }
    Some(hsp)
}

fn finalize_rendered_hsps(mut hsps: Vec<SearchHsp>, min_diag_separation: i32) -> Vec<SearchHsp> {
    if std::env::var_os("BLAST_RS_SKIP_FINAL_ENDPOINT_PURGE").is_none() {
        purge_common_endpoint_hsps(&mut hsps);
    }
    hsps.sort_by(score_compare_search_hsps);
    if std::env::var_os("BLAST_RS_SKIP_FINAL_DEDUP").is_none() {
        dedup_hsps_with_min_diag_separation(&mut hsps, min_diag_separation);
    }
    hsps
}

#[allow(clippy::too_many_arguments)]
fn blastn_preliminary_search_packed_prepared_with_split_xdrop(
    prepared: &PreparedBlastnQuery<'_>,
    query_plus_nomask: &[u8],
    query_minus_nomask: &[u8],
    subject_packed: &[u8],
    subject_len: usize,
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    ungapped_x_dropoff: i32,
    prelim_x_dropoff: i32,
    traceback_x_dropoff: i32,
    kbp: &KarlinBlk,
    search_space: f64,
    evalue_threshold: f64,
    last_hit_scratch: &mut [PackedDiagScratch],
) -> Vec<SearchHsp> {
    trace_packed_prepared_search_path(prepared);
    if std::env::var_os("BLAST_RS_FORCE_NAIVE_UNGAPPED").is_some() {
        let subject_decoded = decode_packed_ncbi2na(subject_packed, subject_len);
        return blastn_gapped_search_nomask_with_split_xdrop(
            query_plus_nomask,
            query_minus_nomask,
            query_plus_nomask,
            query_minus_nomask,
            &subject_decoded,
            prepared.word_size,
            reward,
            penalty,
            gap_open,
            gap_extend,
            ungapped_x_dropoff,
            prelim_x_dropoff,
            traceback_x_dropoff,
            kbp,
            search_space,
            evalue_threshold,
        );
    }

    collect_packed_prepared_ungapped(
        prepared,
        subject_packed,
        subject_len,
        reward,
        penalty,
        ungapped_x_dropoff,
        kbp,
        search_space,
        evalue_threshold,
        last_hit_scratch,
    )
}

#[allow(clippy::too_many_arguments)]
fn blast_preliminary_search_engine_packed_prepared_with_split_xdrop(
    prepared: &PreparedBlastnQuery<'_>,
    query_plus_nomask: &[u8],
    query_minus_nomask: &[u8],
    subject_packed: &[u8],
    subject_len: usize,
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    ungapped_x_dropoff: i32,
    prelim_x_dropoff: i32,
    traceback_x_dropoff: i32,
    kbp: &KarlinBlk,
    search_space: f64,
    evalue_threshold: f64,
    last_hit_scratch: &mut [PackedDiagScratch],
) -> Vec<SearchHsp> {
    blastn_preliminary_search_packed_prepared_with_split_xdrop(
        prepared,
        query_plus_nomask,
        query_minus_nomask,
        subject_packed,
        subject_len,
        reward,
        penalty,
        gap_open,
        gap_extend,
        ungapped_x_dropoff,
        prelim_x_dropoff,
        traceback_x_dropoff,
        kbp,
        search_space,
        evalue_threshold,
        last_hit_scratch,
    )
}

#[allow(clippy::too_many_arguments)]
fn blast_run_preliminary_search_with_interrupt_packed_prepared_with_split_xdrop(
    prepared: &PreparedBlastnQuery<'_>,
    query_plus_nomask: &[u8],
    query_minus_nomask: &[u8],
    subject_packed: &[u8],
    subject_len: usize,
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    ungapped_x_dropoff: i32,
    prelim_x_dropoff: i32,
    traceback_x_dropoff: i32,
    kbp: &KarlinBlk,
    search_space: f64,
    evalue_threshold: f64,
    last_hit_scratch: &mut [PackedDiagScratch],
) -> Vec<SearchHsp> {
    blast_preliminary_search_engine_packed_prepared_with_split_xdrop(
        prepared,
        query_plus_nomask,
        query_minus_nomask,
        subject_packed,
        subject_len,
        reward,
        penalty,
        gap_open,
        gap_extend,
        ungapped_x_dropoff,
        prelim_x_dropoff,
        traceback_x_dropoff,
        kbp,
        search_space,
        evalue_threshold,
        last_hit_scratch,
    )
}

#[allow(clippy::too_many_arguments)]
fn blast_run_preliminary_search_packed_prepared_with_split_xdrop(
    prepared: &PreparedBlastnQuery<'_>,
    query_plus_nomask: &[u8],
    query_minus_nomask: &[u8],
    subject_packed: &[u8],
    subject_len: usize,
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    ungapped_x_dropoff: i32,
    prelim_x_dropoff: i32,
    traceback_x_dropoff: i32,
    kbp: &KarlinBlk,
    search_space: f64,
    evalue_threshold: f64,
    last_hit_scratch: &mut [PackedDiagScratch],
) -> Vec<SearchHsp> {
    blast_run_preliminary_search_with_interrupt_packed_prepared_with_split_xdrop(
        prepared,
        query_plus_nomask,
        query_minus_nomask,
        subject_packed,
        subject_len,
        reward,
        penalty,
        gap_open,
        gap_extend,
        ungapped_x_dropoff,
        prelim_x_dropoff,
        traceback_x_dropoff,
        kbp,
        search_space,
        evalue_threshold,
        last_hit_scratch,
    )
}

#[allow(clippy::too_many_arguments)]
fn blastn_full_search_packed_prepared_with_split_xdrop(
    prepared: &PreparedBlastnQuery<'_>,
    query_plus_nomask: &[u8],
    query_minus_nomask: &[u8],
    subject_packed: &[u8],
    subject_len: usize,
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    prelim_x_dropoff: i32,
    traceback_x_dropoff: i32,
    kbp: &KarlinBlk,
    search_space: f64,
    evalue_threshold: f64,
    ungapped: Vec<SearchHsp>,
) -> Vec<SearchHsp> {
    let profile_enabled = std::env::var_os("BLAST_RS_PROFILE").is_some();
    let gapped_start = if profile_enabled {
        Some(Instant::now())
    } else {
        None
    };
    if let Some(start) = gapped_start {
        eprintln!(
            "[blastn-profile] packed_gapped subject_len={} ungapped_count={} ungapped_ms={}",
            subject_len,
            ungapped.len(),
            start.elapsed().as_millis()
        );
    }
    let mut candidates = Vec::new();
    let mut perfect_seeds = Vec::new();
    if ungapped.is_empty() {
        return Vec::new();
    }

    let cutoff_score = {
        let e = evalue_threshold.max(1.0e-297);
        ((kbp.k * search_space / e).ln() / kbp.lambda).ceil() as i32
    };

    let mut subject_decoded = None;
    let min_diag_separation =
        min_diag_separation_for_gapped(prepared.word_size, reward, penalty, gap_open, gap_extend);
    let mut stats = PackedGappedStageStats::default();
    collect_packed_gapped_candidates(
        prepared,
        &ungapped,
        &mut candidates,
        &mut perfect_seeds,
        query_plus_nomask,
        query_minus_nomask,
        subject_packed,
        subject_len,
        reward,
        penalty,
        gap_open,
        gap_extend,
        prelim_x_dropoff,
        traceback_x_dropoff,
        cutoff_score,
        evalue_threshold,
        search_space,
        kbp,
        min_diag_separation,
        &mut subject_decoded,
        profile_enabled,
        &mut stats,
    );
    if min_diag_separation != 0 && !perfect_seeds.is_empty() {
        prune_perfect_ungapped_gapped_candidates(&mut perfect_seeds, min_diag_separation);
    }
    candidates.extend(perfect_seeds.into_iter().map(GappedCandidate::Final));
    if profile_enabled {
        eprintln!(
            "[blastn-profile] packed_gapped filtered_seed_count={} preliminary_pass_count={} candidate_count={} preliminary_score_ms={} traceback_ms={} traceback_count={} total_ms={}",
            stats.filtered_seed_count,
            stats.preliminary_pass_count,
            candidates.len(),
            stats.preliminary_score_ms,
            stats.traceback_ms,
            stats.traceback_count,
            elapsed_ms_since(gapped_start)
        );
    }

    finalize_gapped_candidates(
        candidates,
        query_plus_nomask,
        query_minus_nomask,
        subject_decoded.get_or_insert_with(|| decode_packed_ncbi2na(subject_packed, subject_len)),
        reward,
        penalty,
        gap_open,
        gap_extend,
        kbp,
        search_space,
        evalue_threshold,
        min_diag_separation,
    )
}

#[allow(clippy::too_many_arguments)]
fn blast_run_full_search_packed_prepared_with_split_xdrop(
    prepared: &PreparedBlastnQuery<'_>,
    query_plus_nomask: &[u8],
    query_minus_nomask: &[u8],
    subject_packed: &[u8],
    subject_len: usize,
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    ungapped_x_dropoff: i32,
    prelim_x_dropoff: i32,
    traceback_x_dropoff: i32,
    kbp: &KarlinBlk,
    search_space: f64,
    evalue_threshold: f64,
    last_hit_scratch: &mut [PackedDiagScratch],
) -> Vec<SearchHsp> {
    let ungapped = blast_run_preliminary_search_packed_prepared_with_split_xdrop(
        prepared,
        query_plus_nomask,
        query_minus_nomask,
        subject_packed,
        subject_len,
        reward,
        penalty,
        gap_open,
        gap_extend,
        ungapped_x_dropoff,
        prelim_x_dropoff,
        traceback_x_dropoff,
        kbp,
        search_space,
        evalue_threshold,
        last_hit_scratch,
    );
    blastn_full_search_packed_prepared_with_split_xdrop(
        prepared,
        query_plus_nomask,
        query_minus_nomask,
        subject_packed,
        subject_len,
        reward,
        penalty,
        gap_open,
        gap_extend,
        prelim_x_dropoff,
        traceback_x_dropoff,
        kbp,
        search_space,
        evalue_threshold,
        ungapped,
    )
}

fn blastn_gapped_search_packed_prepared_with_split_xdrop(
    prepared: &PreparedBlastnQuery<'_>,
    query_plus_nomask: &[u8],
    query_minus_nomask: &[u8],
    subject_packed: &[u8],
    subject_len: usize,
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    ungapped_x_dropoff: i32,
    prelim_x_dropoff: i32,
    traceback_x_dropoff: i32,
    kbp: &KarlinBlk,
    search_space: f64,
    evalue_threshold: f64,
    last_hit_scratch: &mut [PackedDiagScratch],
) -> Vec<SearchHsp> {
    blast_run_full_search_packed_prepared_with_split_xdrop(
        prepared,
        query_plus_nomask,
        query_minus_nomask,
        subject_packed,
        subject_len,
        reward,
        penalty,
        gap_open,
        gap_extend,
        ungapped_x_dropoff,
        prelim_x_dropoff,
        traceback_x_dropoff,
        kbp,
        search_space,
        evalue_threshold,
        last_hit_scratch,
    )
}

#[inline]
fn preliminary_gapped_score_and_seed_greedy(
    query: &[u8],
    subject: &[u8],
    _seed: &SearchHsp,
    seed_q: usize,
    seed_s: usize,
    reward: i32,
    penalty: i32,
    x_dropoff: i32,
) -> Option<PreliminaryGappedResult> {
    let (score, query_start, query_end, subject_start, subject_end, _, adjusted_q, adjusted_s) = greedy_align_with_seed(
        query,
        subject,
        seed_q,
        seed_s,
        reward,
        penalty,
        x_dropoff,
    )?;
    Some(PreliminaryGappedResult {
        score,
        prelim_q_start: query_start,
        prelim_q_end: query_end,
        prelim_s_start: subject_start,
        prelim_s_end: subject_end,
        gapped_start_q: adjusted_q,
        gapped_start_s: adjusted_s,
        q_length: 0,
        s_length: 0,
        private_q_start: 0,
        private_s_start: 0,
        score_left: score,
        score_right: 0,
    })
}

#[inline]
fn preliminary_gapped_score_and_seed_affine(
    query: &[u8],
    subject: &[u8],
    seed: &SearchHsp,
    seed_q: usize,
    seed_s: usize,
    ungapped_score: i32,
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    x_dropoff: i32,
) -> Option<PreliminaryGappedResult> {
    let effective_x_dropoff = x_dropoff.min(ungapped_score.max(0));
    let (seed_q, seed_s) =
        adjust_nt_packed_prelim_start_like_ncbi(seed_q, seed_s, query.len(), subject.len());
    let mut adjusted_seed_s = seed_s;
    let mut adjusted_subject_len = subject.len();
    let start_shift = adjust_subject_range(
        &mut adjusted_seed_s,
        &mut adjusted_subject_len,
        seed_q + 1,
        query.len(),
    );
    let adjusted_subject = &subject[start_shift..start_shift + adjusted_subject_len];
    let score = blast_gapped_score_only(
        query,
        adjusted_subject,
        seed_q,
        adjusted_seed_s,
        reward,
        penalty,
        gap_open,
        gap_extend,
        effective_x_dropoff,
    );
    if score > 0 {
        Some(PreliminaryGappedResult {
            score,
            prelim_q_start: seed.query_start as usize,
            prelim_q_end: seed.query_end as usize,
            prelim_s_start: seed.subject_start as usize,
            prelim_s_end: seed.subject_end as usize,
            gapped_start_q: seed_q,
            gapped_start_s: seed_s,
            q_length: 0,
            s_length: 0,
            private_q_start: 0,
            private_s_start: 0,
            score_left: score,
            score_right: 0,
        })
    } else {
        None
    }
}

#[inline]
fn preliminary_gapped_score_and_seed(
    query: &[u8],
    subject: &[u8],
    seed: &SearchHsp,
    seed_q: usize,
    seed_s: usize,
    ungapped_score: i32,
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    x_dropoff: i32,
) -> Option<PreliminaryGappedResult> {
    if gap_open == 0 && gap_extend == 0 {
        preliminary_gapped_score_and_seed_greedy(
            query, subject, seed, seed_q, seed_s, reward, penalty, x_dropoff,
        )
    } else {
        preliminary_gapped_score_and_seed_affine(
            query,
            subject,
            seed,
            seed_q,
            seed_s,
            ungapped_score,
            reward,
            penalty,
            gap_open,
            gap_extend,
            x_dropoff,
        )
    }
}

#[inline]
fn preliminary_gapped_score_and_seed_packed_greedy(
    query: &[u8],
    subject_packed: &[u8],
    subject_len: usize,
    seed_q: usize,
    seed_s: usize,
    reward: i32,
    penalty: i32,
    x_dropoff: i32,
) -> Option<PreliminaryGappedResult> {
    let packed =
        greedy_align_with_seed_packed_subject(
            query,
            subject_packed,
            subject_len,
            seed_q,
            seed_s,
            reward,
            penalty,
            x_dropoff,
        )?;
    if std::env::var_os("BLAST_RS_COMPARE_GREEDY_PACKED").is_some()
        && trace_target_bounds().is_some()
    {
        let subject_decoded = decode_packed_ncbi2na(subject_packed, subject_len);
        let decoded = greedy_align_with_seed(
            query,
            &subject_decoded,
            seed_q,
            seed_s,
            reward,
            penalty,
            x_dropoff,
        );
        eprintln!(
            "[trace-target] greedy-compare packed={:?} decoded={:?}",
            packed,
            decoded.as_ref().map(
                |(score, q0, q1, s0, s1, _esp, gq, gs)| (*score, *q0, *q1, *s0, *s1, *gq, *gs)
            )
        );
    }
    let (score, query_start, query_end, subject_start, subject_end, adjusted_q, adjusted_s) =
        packed;
    Some(PreliminaryGappedResult {
        score,
        prelim_q_start: query_start,
        prelim_q_end: query_end,
        prelim_s_start: subject_start,
        prelim_s_end: subject_end,
        gapped_start_q: adjusted_q,
        gapped_start_s: adjusted_s,
        q_length: 0,
        s_length: 0,
        private_q_start: 0,
        private_s_start: 0,
        score_left: score,
        score_right: 0,
    })
}

#[inline]
fn preliminary_gapped_score_and_seed_packed_affine(
    query: &[u8],
    subject_packed: &[u8],
    subject_len: usize,
    seed_q: usize,
    seed_s: usize,
    ungapped_score: i32,
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    x_dropoff: i32,
) -> Option<PreliminaryGappedResult> {
    let (
        score,
        query_start,
        query_stop,
        subject_start,
        subject_stop,
        gapped_q,
        gapped_s,
        q_length,
        s_length,
        private_q_start,
        private_s_start,
        score_left,
        score_right,
    ) =
        blast_gapped_score_extents_packed_subject(
        query,
        subject_packed,
        subject_len,
        seed_q,
        seed_s,
        reward,
        penalty,
        gap_open,
        gap_extend,
        x_dropoff,
        ungapped_score,
    );
    if score > 0 {
        Some(PreliminaryGappedResult {
            score,
            prelim_q_start: query_start,
            prelim_q_end: query_stop,
            prelim_s_start: subject_start,
            prelim_s_end: subject_stop,
            gapped_start_q: gapped_q,
            gapped_start_s: gapped_s,
            q_length,
            s_length,
            private_q_start,
            private_s_start,
            score_left,
            score_right,
        })
    } else {
        None
    }
}

#[inline]
fn preliminary_gapped_score_and_seed_packed(
    query: &[u8],
    subject_packed: &[u8],
    subject_len: usize,
    seed_q: usize,
    seed_s: usize,
    ungapped_score: i32,
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    x_dropoff: i32,
) -> Option<PreliminaryGappedResult> {
    if gap_open == 0 && gap_extend == 0 {
        preliminary_gapped_score_and_seed_packed_greedy(
            query,
            subject_packed,
            subject_len,
            seed_q,
            seed_s,
            reward,
            penalty,
            x_dropoff,
        )
    } else {
        preliminary_gapped_score_and_seed_packed_affine(
            query,
            subject_packed,
            subject_len,
            seed_q,
            seed_s,
            ungapped_score,
            reward,
            penalty,
            gap_open,
            gap_extend,
            x_dropoff,
        )
    }
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

pub fn ambiguity_data_overlaps_hsps(ambiguity_data: &[u8], hsps: &[SearchHsp]) -> bool {
    if ambiguity_data.len() < 4 || hsps.is_empty() {
        return false;
    }

    let words: Vec<u32> = ambiguity_data
        .chunks_exact(4)
        .map(|chunk| u32::from_be_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]))
        .collect();
    if words.is_empty() {
        return false;
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

        let end = position.saturating_add(run_len_minus_one + 1);
        if hsps.iter().any(|hsp| {
            let hsp_start = hsp.subject_start.max(0) as usize;
            let hsp_end = hsp.subject_end.max(0) as usize;
            position < hsp_end && end > hsp_start
        }) {
            return true;
        }

        seen += 1;
    }

    false
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
    use crate::input::iupacna_to_blastna;
    use crate::sequence::complement_blastna;

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
    fn test_ambiguity_data_overlaps_hsps_detects_old_and_new_format_runs() {
        let hsps = vec![SearchHsp {
            query_start: 0,
            query_end: 4,
            subject_start: 2,
            subject_end: 5,
            score: 8,
            bit_score: 0.0,
            evalue: 0.0,
            num_ident: 4,
            align_length: 4,
            mismatches: 0,
            gap_opens: 0,
            context: 0,
            qseq: None,
            sseq: None,
        }];

        let old_record = (15u32 << 28) | (2u32 << 24) | 2u32;
        let mut old = Vec::new();
        old.extend_from_slice(&1u32.to_be_bytes());
        old.extend_from_slice(&old_record.to_be_bytes());
        assert!(ambiguity_data_overlaps_hsps(&old, &hsps));

        let new_record = (5u32 << 28) | (1u32 << 16);
        let mut new = Vec::new();
        new.extend_from_slice(&(0x8000_0001u32).to_be_bytes());
        new.extend_from_slice(&new_record.to_be_bytes());
        new.extend_from_slice(&6u32.to_be_bytes());
        assert!(!ambiguity_data_overlaps_hsps(&new, &hsps));
    }

    #[test]
    fn test_purge_common_endpoint_hsps_matches_ncbi_tie_policy() {
        let mut hsps = vec![
            SearchHsp {
                query_start: 5,
                query_end: 25,
                subject_start: 10,
                subject_end: 30,
                score: 40,
                bit_score: 0.0,
                evalue: 0.0,
                num_ident: 20,
                align_length: 20,
                mismatches: 0,
                gap_opens: 0,
                context: 0,
                qseq: None,
                sseq: None,
            },
            SearchHsp {
                query_start: 5,
                query_end: 22,
                subject_start: 10,
                subject_end: 27,
                score: 35,
                bit_score: 0.0,
                evalue: 0.0,
                num_ident: 17,
                align_length: 17,
                mismatches: 0,
                gap_opens: 0,
                context: 0,
                qseq: None,
                sseq: None,
            },
            SearchHsp {
                query_start: 12,
                query_end: 30,
                subject_start: 17,
                subject_end: 40,
                score: 45,
                bit_score: 0.0,
                evalue: 0.0,
                num_ident: 18,
                align_length: 23,
                mismatches: 0,
                gap_opens: 1,
                context: 0,
                qseq: None,
                sseq: None,
            },
            SearchHsp {
                query_start: 15,
                query_end: 30,
                subject_start: 20,
                subject_end: 40,
                score: 41,
                bit_score: 0.0,
                evalue: 0.0,
                num_ident: 15,
                align_length: 20,
                mismatches: 0,
                gap_opens: 1,
                context: 0,
                qseq: None,
                sseq: None,
            },
        ];

        purge_common_endpoint_hsps(&mut hsps);
        hsps.sort_by(score_compare_search_hsps);

        assert_eq!(hsps.len(), 2);
        assert_eq!((hsps[0].query_start, hsps[0].query_end), (12, 30));
        assert_eq!((hsps[1].query_start, hsps[1].query_end), (5, 25));
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
        assert_eq!(prepared.lookups[0].scan_start, 0);
        assert_eq!(prepared.lookups[0].scan_step, 18);
    }

    #[test]
    fn test_megablast_prepared_lookup_can_downgrade_to_small_table() {
        let query = vec![0u8; 8_000];
        let prepared = PreparedBlastnQuery::new_megablast(&query, &[], 28);
        assert_eq!(prepared.lookups.len(), 1);
        assert_eq!(prepared.lookups[0].lut_word, 8);
        assert_eq!(prepared.lookups[0].scan_start, 0);
        assert_eq!(prepared.lookups[0].scan_step, 21);
    }

    #[test]
    fn test_masked_lookup_allows_lut_word_inside_long_unmasked_run() {
        // `BlastChooseNaLookupTable` (`blast_nalookup.c:162-174`) only selects
        // the megablast table (lut_word=11) when approx_table_entries >= 8500;
        // smaller queries downgrade to the small-NA table (lut_word=8). Size
        // the query past that threshold so we actually exercise the MB path.
        let mut query: Vec<u8> = (0..9000).map(|i| (i % 4) as u8).collect();
        let mask_pos = 4045usize;
        query[mask_pos] = 14;

        let prepared = PreparedBlastnQuery::new_megablast(&query, &[], 28);
        let lookup = &prepared.lookups[0];
        assert_eq!(lookup.lut_word, 11);

        let key_allowed = word_hash_n(&query[20..31], lookup.lut_word) as usize;
        assert!(
            lookup.lut[key_allowed] >= 0,
            "NCBI indexes lookup words inside any unmasked run at least word_size long"
        );

        let masked_word_start = mask_pos - 10;
        let key_masked =
            word_hash_n(&query[masked_word_start..masked_word_start + 11], lookup.lut_word)
                as usize;
        let mut pos = lookup.lut[key_masked];
        while pos >= 0 {
            assert_ne!(pos as usize, masked_word_start);
            pos = lookup.next[pos as usize];
        }
    }

    #[test]
    fn test_masked_lookup_allows_offsets_when_lut_word_stays_unmasked() {
        let mut query = vec![0u8; 80];
        query[45] = 14;

        let prepared = PreparedBlastnQuery::new_megablast(&query, &[], 28);
        let lookup = &prepared.lookups[0];
        assert_eq!(lookup.lut_word, 8);

        let key = word_hash_n(&query[19..30], lookup.lut_word) as usize;
        let mut pos = lookup.lut[key];
        let mut found = false;
        while pos >= 0 {
            if pos == 19 {
                found = true;
                break;
            }
            pos = lookup.next[pos as usize];
        }
        assert!(
            found,
            "offset 19 should be indexed when its lut_word span is unmasked even if the full word crosses the masked base"
        );
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
        let nucl_score_table = InitialWordParameters::build_nucl_score_table(2, -3);

        let data = extend_seed_packed(
            &query,
            &subject,
            4,
            0,
            0,
            1,
            2,
            -3,
            20,
            4,
            &nucl_score_table,
            0,
        )
        .expect("first matching base should produce an extension");
        let hsp = build_packed_hsp(&query, &subject, data, &kbp, 1e6, 1e10, 0)
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
    fn test_large_db_row3_exact_word_context_probe() {
        let base = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/large_db/celegans4");
        if !base.with_extension("nin").exists() {
            eprintln!("Skipping: large_db fixture not present at {}", base.display());
            return;
        }

        let query_path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/large_db/query_2000.fa");
        let query_file = std::fs::File::open(&query_path).expect("open query fasta");
        let records = crate::input::parse_fasta(query_file);
        let raw_query = &records[0].sequence;
        let query_plus: Vec<u8> = raw_query.iter().map(|&b| iupacna_to_blastna(b)).collect();
        let query_minus: Vec<u8> = query_plus
            .iter()
            .rev()
            .map(|&b| complement_blastna(b))
            .collect();

        let db = crate::db::BlastDb::open(&base).expect("open large_db blast db");
        let oid = (0..db.num_oids)
            .find(|&oid| db.get_accession(oid).as_deref() == Some("NC_003279.8"))
            .expect("find NC_003279.8 oid");
        let subject_packed = db.get_sequence(oid);
        let subject_len = db.get_seq_len(oid) as usize;

        for &(label, qp, sp) in &[
            ("c-hit-a", 1025usize, 3_897_768usize),
            ("c-hit-b", 1048usize, 3_897_726usize),
            ("c-hit-c", 1253usize, 3_897_999usize),
        ] {
            let plus_hit =
                find_exact_word_hit_packed_naive(&query_plus, subject_packed, subject_len, 28, 8, qp, sp);
            let minus_hit = find_exact_word_hit_packed_naive(
                &query_minus,
                subject_packed,
                subject_len,
                28,
                8,
                qp,
                sp,
            );
            let plus_ext = exact_word_hit_packed_naive_extents(
                &query_plus,
                subject_packed,
                subject_len,
                28,
                8,
                qp,
                sp,
            );
            let minus_ext = exact_word_hit_packed_naive_extents(
                &query_minus,
                subject_packed,
                subject_len,
                28,
                8,
                qp,
                sp,
            );
            eprintln!(
                "[row3-probe] {} qp={} sp={} plus_hit={:?} plus_ext={:?} minus_hit={:?} minus_ext={:?}",
                label, qp, sp, plus_hit.map(|h| (h.query_start, h.subject_start, h.subject_match_end)),
                plus_ext,
                minus_hit.map(|h| (h.query_start, h.subject_start, h.subject_match_end)),
                minus_ext
            );
        }
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
