//! High-level public API for running BLAST searches.
//!
//! Provides `blastp`, `blastn`, `blastx`, `tblastn`, and other search functions
//! with a simple builder-pattern interface, matching the API of the previous
//! blast-rs implementation.

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use crate::db::{defline::encode_defline_asn1, index_writer::write_index_file, BlastDb, DbType};
use crate::encoding::{
    encode_blastna_sequence, encode_ncbi2na_ambiguity_data, encode_ncbi2na_sequence,
    encode_ncbi4na_sequence, encode_ncbistdaa_sequence, ncbistdaa_to_aminoacid_base,
    ncbistdaa_to_aminoacid_char, ncbistdaa_to_aminoacid_sequence,
    reverse_complement_blastna_sequence, reverse_complement_iupacna_sequence, NCBISTDAA_GAP,
    NCBISTDAA_X,
};
use crate::matrix::AA_SIZE;
use crate::search::{blastn_gapped_search_nomask, SearchHsp};
use crate::stat::{
    blast_compute_length_adjustment, blast_get_nucl_alpha_beta, blast_karlin_blk_nucl_gapped_calc,
    ungapped_kbp_calc, GappedParams, GumbelBlk, KarlinBlk, UngappedKbpContext,
};
use crate::traceback::blast_score_blk_nucl_matrix_create;

/// Re-export of `hspstream::evalue_comp` (port of NCBI `s_EvalueComp`):
/// compare two e-values, treating both below `1e-180` as equal.
pub use crate::hspstream::evalue_comp;

const COMPO_ADJUST_SCALE_FACTOR: f64 = 32.0;

/// NCBI composition redo compares scaled alignment scores to
/// `cutoff_score * localScalingFactor` before normalizing the saved HSP score.
/// `localScalingFactor` is 32 for the matrix-adjustment path used here.
fn composition_scaled_cutoff(cutoff_score: i32) -> i32 {
    (cutoff_score as f64 * COMPO_ADJUST_SCALE_FACTOR) as i32
}

fn link_kbp_for_result(base: &KarlinBlk, result: &SearchResult) -> KarlinBlk {
    let Some(lambda) = result.hsps.iter().find_map(|hsp| hsp.link_lambda) else {
        return base.clone();
    };
    let mut kbp = base.clone();
    kbp.lambda = lambda;
    kbp.round_down = false;
    kbp
}

fn db_component_path(base_path: &Path, ext: &str) -> PathBuf {
    let mut path = base_path.as_os_str().to_os_string();
    path.push(".");
    path.push(ext);
    path.into()
}

/// Per-thread scratch space reused across all OIDs processed by one worker.
///
/// Mirrors NCBI's pattern of allocating a single set of preliminary/traceback
/// `BlastIntervalTree`s, a diagonal-tracking array, and a subject masking
/// buffer once per thread (e.g. inside `BlastGapAlignStruct` /
/// `BlastInitialWordResults`) and resetting them between subjects rather
/// than reallocating per OID.
pub struct ProteinScratch {
    diag_buf: Vec<(i32, bool)>,
    diag_table: crate::protein_lookup::BlastExtendWord,
    prelim_tree: crate::itree::IntervalTree,
    tb_tree: crate::itree::IntervalTree,
    mask_buf: Vec<u8>,
}

impl ProteinScratch {
    /// blast-rs: Native scratch constructor; not a direct NCBI C port.
    pub fn new() -> Self {
        Self {
            diag_buf: Vec::new(),
            diag_table: crate::protein_lookup::BlastExtendWord::new(0),
            prelim_tree: crate::itree::IntervalTree::new(1, 1),
            tb_tree: crate::itree::IntervalTree::new(1, 1),
            mask_buf: Vec::new(),
        }
    }
}

impl Default for ProteinScratch {
    /// blast-rs: Native trait implementation; not a direct NCBI C port.
    fn default() -> Self {
        Self::new()
    }
}

// ── Result types ────────────────────────────────────────────────────────────

/// A single high-scoring segment pair from a BLAST search.
#[derive(Debug, Clone)]
pub struct Hsp {
    pub score: i32,
    pub bit_score: f64,
    pub evalue: f64,
    pub query_start: usize,
    pub query_end: usize,
    pub subject_start: usize,
    pub subject_end: usize,
    /// Internal seed/traceback gapped start carried from NCBI `BlastHSP`.
    /// Output formatters use `query_start`; translated linkers use this field
    /// to preserve the HSP metadata instead of reconstructing it from bounds.
    #[doc(hidden)]
    pub query_gapped_start: usize,
    /// Internal seed/traceback gapped start carried from NCBI `BlastHSP`.
    #[doc(hidden)]
    pub subject_gapped_start: usize,
    /// Internal query offset in the coordinate space consumed by `BLAST_LinkHsps`.
    #[doc(hidden)]
    pub query_link_start: usize,
    /// Internal query end in the coordinate space consumed by `BLAST_LinkHsps`.
    #[doc(hidden)]
    pub query_link_end: usize,
    /// Internal query gapped start in the coordinate space consumed by `BLAST_LinkHsps`.
    #[doc(hidden)]
    pub query_link_gapped_start: usize,
    /// Internal subject offset in the coordinate space consumed by `BLAST_LinkHsps`.
    #[doc(hidden)]
    pub subject_link_start: usize,
    /// Internal subject end in the coordinate space consumed by `BLAST_LinkHsps`.
    #[doc(hidden)]
    pub subject_link_end: usize,
    /// Internal subject gapped start in the coordinate space consumed by `BLAST_LinkHsps`.
    #[doc(hidden)]
    pub subject_link_gapped_start: usize,
    /// Internal score consumed by `BLAST_LinkHsps`; composition-adjusted HSPs
    /// keep the public score normalized but link in the scaled score space.
    #[doc(hidden)]
    pub link_score: Option<i32>,
    /// Per-HSP lambda matching `link_score` when it differs from the score block.
    #[doc(hidden)]
    pub link_lambda: Option<f64>,
    pub num_identities: usize,
    pub num_gaps: usize,
    pub alignment_length: usize,
    pub query_aln: Vec<u8>,
    pub midline: Vec<u8>,
    pub subject_aln: Vec<u8>,
    pub query_frame: i32,
    pub subject_frame: i32,
    /// Number of HSPs in the linked-set chain (NCBI `BlastHSP::num`).
    /// 1 for unlinked / singleton HSPs; >=2 for sum-stats chains.
    /// Drives the pairwise `Expect(N) = ...` annotation.
    pub num_links: i32,
    /// Per-HSP composition adjustment method actually applied
    /// (NCBI `BlastHSP::comp_adjustment_method`, `blast_kappa.c:332-342`).
    /// 0 = eNoCompositionBasedStats; 1 = eCompositionBasedStats
    /// (ScaleOldMatrix, prints "Composition-based stats.");
    /// 2 = eCompositionMatrixAdjust (UserSpecifiedRelEntropy / full optimisation,
    /// prints "Compositional matrix adjust."). For comp_adjust=2 search params,
    /// `Blast_ChooseMatrixAdjustRule` may fall back to 1 per HSP.
    pub comp_adjust_method: u8,
}

impl Hsp {
    /// Percent identity of this HSP.
    /// blast-rs: Public API helper; not a direct NCBI C port.
    pub fn percent_identity(&self) -> f64 {
        if self.alignment_length == 0 {
            return 0.0;
        }
        100.0 * self.num_identities as f64 / self.alignment_length as f64
    }
}

/// Result for one subject sequence (may contain multiple HSPs).
#[derive(Debug, Clone)]
pub struct SearchResult {
    pub subject_oid: u32,
    pub subject_title: String,
    pub subject_accession: String,
    pub subject_len: usize,
    pub hsps: Vec<Hsp>,
    pub taxids: Vec<i32>,
}

impl SearchResult {
    /// blast-rs: Public API helper modeled on NCBI best-evalue semantics; not a direct C port.
    pub fn best_evalue(&self) -> f64 {
        // NCBI `s_BlastGetBestEvalue` (`blast_hits.c:1742`) seeds with
        // `(double)INT4_MAX` so an empty list returns 2147483647.0 rather
        // than infinity.
        self.hsps
            .iter()
            .map(|h| h.evalue)
            .fold(i32::MAX as f64, f64::min)
    }

    /// Score of the first HSP after C HSP-list e-value ordering, or
    /// `i32::MIN` if the result has no HSPs. Used as the secondary key in
    /// `compare_search_results`.
    /// blast-rs: SearchResult adapter; not a direct NCBI C port.
    fn first_hsp_score_by_evalue(&self) -> i32 {
        self.hsps
            .iter()
            .min_by(|a, b| compare_hsps_by_evalue_then_score(a, b))
            .map(|h| h.score)
            .unwrap_or(i32::MIN)
    }
}

/// blast-rs: SearchResult comparator modeled on the C e-value HSP-list
/// comparator, but not a direct NCBI C port. Primary key is `best_evalue` (compared via
/// `crate::stat::evalue_comp` — both <1e-180 compare equal). Ties break on
/// top HSP score descending, then `subject_oid` descending.
pub fn compare_search_results(a: &SearchResult, b: &SearchResult) -> std::cmp::Ordering {
    use std::cmp::Ordering::*;
    // Empty-hsps results sort to the end.
    match (a.hsps.is_empty(), b.hsps.is_empty()) {
        (true, true) => return Equal,
        (true, false) => return Greater,
        (false, true) => return Less,
        (false, false) => {}
    }
    let by_evalue = crate::hspstream::evalue_comp(a.best_evalue(), b.best_evalue());
    if by_evalue != Equal {
        return by_evalue;
    }
    let by_score = b
        .first_hsp_score_by_evalue()
        .cmp(&a.first_hsp_score_by_evalue());
    if by_score != Equal {
        return by_score;
    }
    b.subject_oid.cmp(&a.subject_oid)
}

/// blast-rs: Public API post-filter; not a direct NCBI C port.
fn apply_api_culling_limit(
    results: &mut Vec<SearchResult>,
    culling_limit: usize,
    program: crate::program::ProgramType,
) {
    if culling_limit == 0 {
        return;
    }
    let hsp_count: usize = results.iter().map(|result| result.hsps.len()).sum();
    if hsp_count <= 1 {
        return;
    }

    let mut ordered = Vec::with_capacity(hsp_count);
    for (result_idx, result) in results.iter().enumerate() {
        for (hsp_idx, hsp) in result.hsps.iter().enumerate() {
            ordered.push((result_idx, hsp_idx, hsp));
        }
    }
    ordered.sort_by(|a, b| {
        compare_api_culling_hsps(a.2, b.2, results[a.0].subject_oid, results[b.0].subject_oid)
    });

    let mut keep: Vec<Vec<bool>> = results
        .iter()
        .map(|result| vec![true; result.hsps.len()])
        .collect();

    let mut accepted: Vec<(usize, usize, crate::hspfilter_culling::LinkedHsp)> =
        Vec::with_capacity(hsp_count);
    for (result_idx, hsp_idx, hsp) in ordered {
        let mut candidate = api_hsp_as_culling_node(hsp, results[result_idx].subject_oid);
        candidate.merit = culling_limit as i32;
        let enveloping = accepted
            .iter()
            .filter(|(accepted_result_idx, accepted_hsp_idx, dominator)| {
                let existing = &results[*accepted_result_idx].hsps[*accepted_hsp_idx];
                if api_culling_context_id(existing, program) != api_culling_context_id(hsp, program)
                {
                    return false;
                }
                crate::hspfilter_culling::s_dominate_test(dominator, &candidate)
            })
            .take(culling_limit)
            .count();
        if enveloping >= culling_limit {
            keep[result_idx][hsp_idx] = false;
        } else {
            let mut idx = 0usize;
            while idx < accepted.len() {
                let (accepted_result_idx, accepted_hsp_idx, accepted_node) = &mut accepted[idx];
                let existing = &results[*accepted_result_idx].hsps[*accepted_hsp_idx];
                if api_culling_context_id(existing, program) == api_culling_context_id(hsp, program)
                    && api_culling_candidate_can_displace_accepted(hsp, existing)
                    && crate::hspfilter_culling::s_dominate_test(&candidate, accepted_node)
                {
                    accepted_node.merit -= 1;
                    if accepted_node.merit <= 0 {
                        keep[*accepted_result_idx][*accepted_hsp_idx] = false;
                        accepted.remove(idx);
                        continue;
                    }
                }
                idx += 1;
            }
            accepted.push((result_idx, hsp_idx, candidate));
        }
    }

    for (result_idx, result) in results.iter_mut().enumerate() {
        let mut idx = 0usize;
        result.hsps.retain(|_| {
            let retained = keep[result_idx][idx];
            idx += 1;
            retained
        });
    }
    results.retain(|result| !result.hsps.is_empty());
}

/// blast-rs: Public API culling tie helper; not a direct NCBI C port.
fn api_culling_candidate_can_displace_accepted(candidate: &Hsp, accepted: &Hsp) -> bool {
    candidate.score > accepted.score
        || (candidate.score == accepted.score
            && crate::hspstream::evalue_comp(candidate.evalue, accepted.evalue)
                == std::cmp::Ordering::Less)
}

/// blast-rs: Public API culling helper; not a direct NCBI C port.
fn api_culling_context_id(hsp: &Hsp, program: crate::program::ProgramType) -> i32 {
    if program == crate::program::BLASTN {
        0
    } else {
        hsp.query_frame
    }
}

/// blast-rs: Public API post-filter; not a direct NCBI C port.
fn apply_api_max_hsps_limit(results: &mut Vec<SearchResult>, max_hsps: Option<usize>) {
    let Some(max) = max_hsps else {
        return;
    };
    for result in results.iter_mut() {
        if max == 0 {
            result.hsps.clear();
            continue;
        }
        if result.hsps.len() <= max {
            continue;
        }
        let mut indices: Vec<usize> = (0..result.hsps.len()).collect();
        indices.sort_by(|&a, &b| compare_hsps_by_score(&result.hsps[a], &result.hsps[b]));
        let keep: std::collections::HashSet<usize> = indices.into_iter().take(max).collect();
        let mut idx = 0usize;
        result.hsps.retain(|_| {
            let retain = keep.contains(&idx);
            idx += 1;
            retain
        });
    }
    results.retain(|result| !result.hsps.is_empty());
}

/// blast-rs: Public API post-filter; not a direct NCBI C port.
fn apply_api_min_score_filter(results: &mut Vec<SearchResult>, min_score: i32) {
    if min_score <= 0 {
        return;
    }
    for result in results.iter_mut() {
        result.hsps.retain(|hsp| hsp.score >= min_score);
    }
    results.retain(|result| !result.hsps.is_empty());
}

/// blast-rs: Public report-threshold adapter after the engine prelim-evalue reap.
fn apply_api_evalue_filter(results: &mut Vec<SearchResult>, evalue_threshold: f64) {
    for result in results.iter_mut() {
        result.hsps.retain(|hsp| hsp.evalue <= evalue_threshold);
    }
    results.retain(|result| !result.hsps.is_empty());
}

#[derive(Clone, Copy)]
enum ProteinMidlineStyle<'a> {
    AnyMismatchPlus,
    PositiveMatrix(&'a [[i32; AA_SIZE]; AA_SIZE]),
}

fn fill_protein_midlines(results: &mut [SearchResult], style: ProteinMidlineStyle<'_>) {
    for result in results {
        for hsp in &mut result.hsps {
            hsp.midline = protein_midline(&hsp.query_aln, &hsp.subject_aln, style);
        }
    }
}

fn protein_midline(
    query_aln: &[u8],
    subject_aln: &[u8],
    style: ProteinMidlineStyle<'_>,
) -> Vec<u8> {
    query_aln
        .iter()
        .zip(subject_aln.iter())
        .map(|(&q, &s)| {
            let q_upper = q.to_ascii_uppercase();
            let s_upper = s.to_ascii_uppercase();
            if q_upper == s_upper {
                q_upper
            } else if q == b'-' || s == b'-' {
                b' '
            } else {
                match style {
                    ProteinMidlineStyle::AnyMismatchPlus => b'+',
                    ProteinMidlineStyle::PositiveMatrix(matrix) => {
                        let q_idx = crate::encoding::aminoacid_to_ncbistdaa_base(q) as usize;
                        let s_idx = crate::encoding::aminoacid_to_ncbistdaa_base(s) as usize;
                        if q_idx < AA_SIZE && s_idx < AA_SIZE && matrix[q_idx][s_idx] > 0 {
                            b'+'
                        } else {
                            b' '
                        }
                    }
                }
            }
        })
        .collect()
}

/// NCBI: s_Blast_HSPListReapByPrelimEvalue (`blast_engine.c`).
fn reap_hsps_by_prelim_evalue(results: &mut Vec<SearchResult>, prelim_evalue: f64) {
    for result in results.iter_mut() {
        result.hsps.retain(|hsp| hsp.evalue <= prelim_evalue);
    }
    results.retain(|result| !result.hsps.is_empty());
}

/// blast-rs: Public API culling ordering helper; not a direct NCBI C port.
fn compare_api_culling_hsps(
    a: &Hsp,
    b: &Hsp,
    a_subject_oid: u32,
    b_subject_oid: u32,
) -> std::cmp::Ordering {
    let a_subject_lo = a.subject_start.min(a.subject_end);
    let b_subject_lo = b.subject_start.min(b.subject_end);
    let a_query_lo = a.query_start.min(a.query_end);
    let b_query_lo = b.query_start.min(b.query_end);

    crate::hspstream::evalue_comp(a.evalue, b.evalue)
        .then_with(|| b.score.cmp(&a.score))
        .then_with(|| a_subject_oid.cmp(&b_subject_oid))
        .then_with(|| a_subject_lo.cmp(&b_subject_lo))
        .then_with(|| a_query_lo.cmp(&b_query_lo))
        .then_with(|| b.subject_frame.cmp(&a.subject_frame))
}

/// blast-rs: Adapter from public HSPs into culling internals; not a direct NCBI C port.
fn api_hsp_as_culling_node(hsp: &Hsp, subject_oid: u32) -> crate::hspfilter_culling::LinkedHsp {
    let query_start = hsp.query_start.min(hsp.query_end).saturating_sub(1) as i32;
    let query_end = hsp.query_start.max(hsp.query_end) as i32;
    let subject_start = hsp.subject_start.min(hsp.subject_end).saturating_sub(1) as i32;
    let subject_end = hsp.subject_start.max(hsp.subject_end) as i32;
    crate::hspfilter_culling::LinkedHsp {
        hsp: crate::hspstream::Hsp {
            score: hsp.score,
            num_ident: hsp.num_identities as i32,
            bit_score: hsp.bit_score,
            evalue: hsp.evalue,
            query_offset: query_start,
            query_end,
            query_gapped_start: hsp.query_gapped_start as i32,
            subject_offset: subject_start,
            subject_end,
            subject_gapped_start: hsp.subject_gapped_start as i32,
            context: hsp.query_frame,
            query_frame: hsp.query_frame,
            subject_frame: hsp.subject_frame,
            num_gaps: hsp.num_gaps as i32,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        },
        cid: hsp.query_frame,
        sid: subject_oid as i32,
        begin: query_start,
        end: query_end,
        merit: 1,
        next: None,
    }
}

/// Map a closure over all OIDs, using each worker thread's reusable
/// scratch value produced by `init`. The mapping closure receives `&mut S`
/// and the OID. Mirrors NCBI's per-thread `BlastGapAlignStruct` allocation
/// pattern: one scratch buffer per worker, reset between subjects instead
/// of reallocated per OID.
/// blast-rs: Rayon database-OID dispatch adapter; not a direct NCBI C port.
fn map_database_oids_init<T, S, FInit, F>(
    db: &BlastDb,
    params: &SearchParams,
    init: FInit,
    f: F,
) -> Vec<T>
where
    T: Send,
    S: Send,
    FInit: Fn() -> S + Sync + Send,
    F: Fn(&mut S, u32) -> T + Sync + Send,
{
    if params.thread_pool.is_none() && params.num_threads == 1 {
        let mut state = init();
        (0..db.num_oids).map(|oid| f(&mut state, oid)).collect()
    } else if let Some(pool) = params.thread_pool.as_deref() {
        use rayon::prelude::*;
        pool.install(|| {
            (0..db.num_oids)
                .into_par_iter()
                .map_init(&init, &f)
                .collect()
        })
    } else {
        use rayon::prelude::*;
        let num_threads = if params.num_threads == 0 {
            rayon::current_num_threads()
        } else {
            params.num_threads
        };
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .stack_size(64 * 1024 * 1024)
            .build()
            .unwrap_or_else(|_| rayon::ThreadPoolBuilder::new().build().unwrap());
        pool.install(|| (0..db.num_oids).into_par_iter().map_init(init, f).collect())
    }
}

struct TranslatedContextStats {
    context_id: usize,
    query_offset: i32,
    frame: i32,
    query_length: i32,
    eff_searchsp: i64,
    length_adjustment: i32,
    /// Per-context ungapped Karlin params. NCBI computes one Karlin block per
    /// query frame from the frame's translation composition; e-value sites
    /// must match the bit-score path, which already uses these.
    kbp: KarlinBlk,
}

fn api_hsp_to_blast_hsp(
    hsp: &Hsp,
    context: i32,
    score: i32,
    query: crate::hspstream::BlastSeg,
    subject: crate::hspstream::BlastSeg,
) -> crate::hspstream::BlastHSP {
    crate::hspstream::BlastHSP {
        score,
        num_ident: hsp.num_identities as i32,
        bit_score: hsp.bit_score,
        evalue: hsp.evalue,
        query,
        subject,
        context,
        gap_info: None,
        num: hsp.num_links.max(1),
        comp_adjustment_method: hsp.comp_adjust_method as i16,
        pat_info: None,
        num_positives: 0,
        map_info: None,
    }
}

fn linked_api_hsp_key(
    context: i32,
    score: i32,
    query: &crate::hspstream::BlastSeg,
    subject: &crate::hspstream::BlastSeg,
) -> (i32, i32, i16, i32, i32, i32, i16, i32, i32, i32) {
    (
        score,
        context,
        query.frame,
        query.offset,
        query.end,
        query.gapped_start,
        subject.frame,
        subject.offset,
        subject.end,
        subject.gapped_start,
    )
}

fn linked_api_hsp_same_head(
    context: i32,
    score: i32,
    query: &crate::hspstream::BlastSeg,
    subject: &crate::hspstream::BlastSeg,
    linked: &crate::hspstream::BlastHSP,
) -> bool {
    score == linked.score
        && context == linked.context
        && query.frame == linked.query.frame
        && query.offset == linked.query.offset
        && query.gapped_start == linked.query.gapped_start
        && subject.frame == linked.subject.frame
        && subject.offset == linked.subject.offset
        && subject.gapped_start == linked.subject.gapped_start
}

fn apply_blast_hsp_list_to_api_hsps(
    hsps: &mut Vec<Hsp>,
    linked_list: crate::hspstream::BlastHSPList,
    key_for_hsp: impl Fn(
        &Hsp,
    ) -> (
        i32,
        i32,
        crate::hspstream::BlastSeg,
        crate::hspstream::BlastSeg,
    ),
) {
    let mut remaining: Vec<Option<Hsp>> = std::mem::take(hsps).into_iter().map(Some).collect();
    let mut linked_hsps = Vec::with_capacity(linked_list.hsp_array.len());

    for linked in linked_list.hsp_array.into_iter().flatten() {
        let mut selected = None;
        let mut selected_same_head = None;
        for (idx, candidate) in remaining.iter().enumerate() {
            let Some(candidate) = candidate.as_ref() else {
                continue;
            };
            let (score, context, query, subject) = key_for_hsp(candidate);
            if linked_api_hsp_key(context, score, &query, &subject)
                == (
                    linked.score,
                    linked.context,
                    linked.query.frame,
                    linked.query.offset,
                    linked.query.end,
                    linked.query.gapped_start,
                    linked.subject.frame,
                    linked.subject.offset,
                    linked.subject.end,
                    linked.subject.gapped_start,
                )
            {
                selected = Some(idx);
                break;
            }
            if selected_same_head.is_none()
                && linked_api_hsp_same_head(context, score, &query, &subject, &linked)
            {
                selected_same_head = Some(idx);
            }
        }

        if let Some(idx) = selected.or(selected_same_head) {
            if let Some(mut hsp) = remaining[idx].take() {
                // C `s_CombineLinkedHSPSets` mutates only the existing
                // `BlastHSP` link fields (`evalue`, `num`, `xsum`). Keep API
                // coordinates from the pre-link HSP even if the Rust link
                // bridge's reduced model carried a changed end coordinate.
                hsp.evalue = linked.evalue;
                hsp.num_links = linked.num.max(1);
                linked_hsps.push(hsp);
            }
        }
    }

    hsps.extend(linked_hsps);
}

/// blast-rs: tblastx adapter into `BLAST_LinkHsps` over canonical `BlastHSPList`.
fn apply_tblastx_linked_sum_stats(
    results: &mut [Option<SearchResult>],
    query_contexts: &[TranslatedContextStats],
    max_intron_length: i32,
    db_length_nt: i64,
    cutoff_score_min: i32,
) {
    use crate::link_hsps::{blast_link_hsp_list, LinkHSPParameters, LinkScoreBlock};
    use crate::program::TBLASTX;
    use crate::queryinfo::{ContextInfo, QueryInfo, E_NO_SEGMENTS};

    if query_contexts.is_empty() {
        return;
    }

    let num_contexts = crate::util::NUM_FRAMES.max(
        query_contexts
            .iter()
            .map(|ctx| ctx.context_id)
            .max()
            .unwrap_or(0)
            + 1,
    );
    let mut contexts = vec![
        ContextInfo {
            query_offset: 0,
            query_length: 0,
            eff_searchsp: 0,
            length_adjustment: 0,
            query_index: 0,
            frame: 0,
            is_valid: false,
            segment_flags: E_NO_SEGMENTS,
        };
        num_contexts
    ];
    for ctx in query_contexts {
        contexts[ctx.context_id] = ContextInfo {
            query_offset: ctx.query_offset,
            query_length: ctx.query_length,
            eff_searchsp: ctx.eff_searchsp,
            length_adjustment: ctx.length_adjustment,
            query_index: 0,
            frame: ctx.frame,
            is_valid: true,
            segment_flags: E_NO_SEGMENTS,
        };
    }
    let (max_length, min_length) = {
        let mut max_l = 0u32;
        let mut min_l = u32::MAX;
        for ctx in &contexts {
            let len = ctx.query_length.max(0) as u32;
            if len > max_l {
                max_l = len;
            }
            if len < min_l {
                min_l = len;
            }
        }
        if min_l == u32::MAX {
            min_l = 0;
        }
        (max_l, min_l)
    };
    let query_info = QueryInfo {
        num_queries: 1,
        max_length,
        min_length,
        contexts,
    };
    let default_kbp = query_contexts[0].kbp.clone();
    let mut kbps: Vec<KarlinBlk> = vec![default_kbp; num_contexts];
    for ctx in query_contexts {
        kbps[ctx.context_id] = ctx.kbp.clone();
    }
    let score_block = LinkScoreBlock {
        kbp: kbps.clone(),
        kbp_gap: kbps,
        ..LinkScoreBlock::default()
    };
    let mut cutoff_score_block = crate::stat::blast_score_blk_new(
        crate::encoding::BLASTAA_SEQ_CODE,
        query_info.contexts.len() as i32,
    )
    .expect("protein score block");
    cutoff_score_block.kbp = score_block.kbp.clone();
    let word_params = crate::parameters::InitialWordParameters {
        options: crate::options::InitialWordOptions::new_blastp(),
        x_dropoff_max: 0,
        cutoff_score_min,
        cutoffs: Vec::new(),
        ungapped_extension: true,
        nucl_score_table: [0; 256],
    };

    for result in results.iter_mut() {
        let Some(result) = result.as_mut() else {
            continue;
        };
        if result.hsps.is_empty() {
            continue;
        }

        let mut link_params = LinkHSPParameters {
            longest_intron: translated_link_longest_intron(max_intron_length),
            ..LinkHSPParameters::default()
        };
        crate::parameters::calculate_link_hsp_cutoffs(
            TBLASTX,
            &query_info,
            &cutoff_score_block,
            Some(&mut link_params),
            &word_params,
            db_length_nt,
            result.subject_len as i32,
        );

        let mut hsp_list = crate::hspstream::BlastHSPList {
            oid: result.subject_oid as i32,
            query_index: 0,
            hsp_array: result
                .hsps
                .iter()
                .map(|hsp| {
                    let context = query_contexts
                        .iter()
                        .position(|ctx| ctx.frame == hsp.query_frame)
                        .map(|idx| query_contexts[idx].context_id)
                        .unwrap_or(0) as i32;
                    Some(api_hsp_to_blast_hsp(
                        hsp,
                        context,
                        hsp.score,
                        crate::hspstream::BlastSeg {
                            frame: hsp.query_frame as i16,
                            offset: hsp.query_link_start as i32,
                            end: hsp.query_link_end as i32,
                            gapped_start: hsp.query_link_gapped_start as i32,
                        },
                        crate::hspstream::BlastSeg {
                            frame: hsp.subject_frame as i16,
                            offset: hsp.subject_link_start as i32,
                            end: hsp.subject_link_end as i32,
                            gapped_start: hsp.subject_link_gapped_start as i32,
                        },
                    ))
                })
                .collect(),
            hspcnt: result.hsps.len() as i32,
            allocated: result.hsps.len() as i32,
            hsp_max: i32::MAX,
            do_not_reallocate: false,
            best_evalue: result.best_evalue(),
        };

        // tblastx is ungapped (`blast_options.c:869`); pass
        // `gapped_calculation=false` so `blast_link_hsps` reads `kbp` (ungapped),
        // matching `link_hsps.c:457`. The `score_block` already populates both
        // arrays with ungapped Karlin params for this reason.
        blast_link_hsp_list(
            TBLASTX,
            &mut hsp_list,
            &query_info,
            result.subject_len as i32,
            &score_block,
            &link_params,
            false,
        );

        apply_blast_hsp_list_to_api_hsps(&mut result.hsps, hsp_list, |hsp| {
            let context = query_contexts
                .iter()
                .position(|ctx| ctx.frame == hsp.query_frame)
                .map(|idx| query_contexts[idx].context_id)
                .unwrap_or(0) as i32;
            (
                hsp.score,
                context,
                crate::hspstream::BlastSeg {
                    frame: hsp.query_frame as i16,
                    offset: hsp.query_link_start as i32,
                    end: hsp.query_link_end as i32,
                    gapped_start: hsp.query_link_gapped_start as i32,
                },
                crate::hspstream::BlastSeg {
                    frame: hsp.subject_frame as i16,
                    offset: hsp.subject_link_start as i32,
                    end: hsp.subject_link_end as i32,
                    gapped_start: hsp.subject_link_gapped_start as i32,
                },
            )
        });
    }
}

fn blastx_link_key(
    hsp: &Hsp,
    query_info: &crate::queryinfo::QueryInfo,
) -> (
    i32,
    i32,
    crate::hspstream::BlastSeg,
    crate::hspstream::BlastSeg,
) {
    let context = query_info
        .contexts
        .iter()
        .position(|ctx| ctx.frame == hsp.query_frame)
        .unwrap_or(0) as i32;
    (
        hsp.link_score.unwrap_or(hsp.score),
        context,
        crate::hspstream::BlastSeg {
            frame: hsp.query_frame as i16,
            offset: hsp.query_link_start as i32,
            end: hsp.query_link_end as i32,
            gapped_start: hsp.query_link_gapped_start as i32,
        },
        crate::hspstream::BlastSeg {
            frame: 0,
            offset: hsp.subject_link_start as i32,
            end: hsp.subject_link_end as i32,
            gapped_start: hsp.subject_link_gapped_start as i32,
        },
    )
}

fn tblastn_link_key(
    hsp: &Hsp,
) -> (
    i32,
    i32,
    crate::hspstream::BlastSeg,
    crate::hspstream::BlastSeg,
) {
    (
        hsp.link_score.unwrap_or(hsp.score),
        0,
        crate::hspstream::BlastSeg {
            frame: 0,
            offset: hsp.query_start as i32,
            end: hsp.query_end as i32,
            gapped_start: hsp.query_gapped_start as i32,
        },
        crate::hspstream::BlastSeg {
            frame: hsp.subject_frame as i16,
            offset: hsp.subject_link_start as i32,
            end: hsp.subject_link_end as i32,
            gapped_start: hsp.subject_link_gapped_start as i32,
        },
    )
}

/// blast-rs: blastx adapter into `BLAST_LinkHsps` over canonical `BlastHSPList`.
fn apply_blastx_linked_sum_stats(
    results: &mut [SearchResult],
    query_info: &crate::queryinfo::QueryInfo,
    prot_kbp: &KarlinBlk,
    gumbel_blk: Option<&crate::stat::GumbelBlk>,
    statistical_db_length: i64,
    recompute_evalues_before_linking: bool,
    max_intron_length: i32,
) {
    use crate::link_hsps::{blast_link_hsp_list, LinkHSPParameters, LinkScoreBlock};
    use crate::program::BLASTX;

    let longest_intron = translated_gapped_link_longest_intron(max_intron_length);
    if longest_intron <= 0 {
        for result in results.iter_mut() {
            let link_kbp = link_kbp_for_result(prot_kbp, result);
            for hsp in &mut result.hsps {
                let context = query_info
                    .contexts
                    .iter()
                    .position(|ctx| ctx.frame == hsp.query_frame)
                    .unwrap_or(0);
                let ctx = &query_info.contexts[context];
                hsp.evalue = if let Some(mut gbp) = gumbel_blk.cloned() {
                    gbp.db_length = statistical_db_length.max(1);
                    crate::stat::blast_spouge_sto_e(
                        hsp.link_score.unwrap_or(hsp.score),
                        Some(&link_kbp),
                        Some(&gbp),
                        ctx.query_length.max(1),
                        result.subject_len.max(1) as i32,
                    )
                } else {
                    let searchsp = if ctx.eff_searchsp > 0 {
                        ctx.eff_searchsp as f64
                    } else {
                        ctx.query_length.max(1) as f64 * result.subject_len.max(1) as f64
                    };
                    link_kbp.raw_to_evalue(hsp.link_score.unwrap_or(hsp.score), searchsp)
                };
                let decay = crate::stat::blast_gap_decay_divisor(
                    crate::stat::BLAST_GAP_DECAY_RATE_GAPPED,
                    1,
                );
                hsp.evalue /= decay;
                if !recompute_evalues_before_linking {
                    hsp.evalue /= decay;
                }
                hsp.num_links = 1;
            }
        }
        return;
    }

    let link_params = LinkHSPParameters {
        gap_prob: crate::stat::BLAST_GAP_PROB_GAPPED,
        gap_decay_rate: crate::stat::BLAST_GAP_DECAY_RATE_GAPPED,
        longest_intron,
        ..LinkHSPParameters::default()
    };

    for result in results.iter_mut() {
        if result.hsps.is_empty() {
            continue;
        }
        let link_kbp = link_kbp_for_result(prot_kbp, result);
        let score_block = LinkScoreBlock {
            kbp: vec![link_kbp.clone(); query_info.contexts.len()],
            kbp_gap: vec![link_kbp; query_info.contexts.len()],
            gbp: gumbel_blk.cloned(),
            // blastx's subject is an UNTRANSLATED protein, so NCBI sets
            // `gbp->db_length = dbl` (no /3) and the linker recompute uses
            // `n_ = subject_length` directly (`Blast_SubjectIsTranslated` is false).
            link_gbp_db_length: Some(statistical_db_length.max(1)),
            recompute_evalues_before_uneven_linking: recompute_evalues_before_linking,
        };
        let mut hsp_list = crate::hspstream::BlastHSPList {
            oid: result.subject_oid as i32,
            query_index: 0,
            hsp_array: result
                .hsps
                .iter()
                .map(|hsp| {
                    let (score, context, query, subject) = blastx_link_key(hsp, query_info);
                    Some(api_hsp_to_blast_hsp(hsp, context, score, query, subject))
                })
                .collect(),
            hspcnt: result.hsps.len() as i32,
            allocated: result.hsps.len() as i32,
            hsp_max: i32::MAX,
            do_not_reallocate: false,
            best_evalue: result.best_evalue(),
        };

        blast_link_hsp_list(
            BLASTX,
            &mut hsp_list,
            query_info,
            result.subject_len as i32,
            &score_block,
            &link_params,
            true,
        );

        apply_blast_hsp_list_to_api_hsps(&mut result.hsps, hsp_list, |hsp| {
            blastx_link_key(hsp, query_info)
        });
    }
}

/// blast-rs: tblastn adapter into `BLAST_LinkHsps` over canonical `BlastHSPList`.
fn apply_tblastn_linked_sum_stats(
    results: &mut [SearchResult],
    subject_stat_lengths: Option<&[i32]>,
    query_info: &crate::queryinfo::QueryInfo,
    prot_kbp: &KarlinBlk,
    gumbel_blk: Option<&crate::stat::GumbelBlk>,
    statistical_db_length: i64,
    recompute_evalues_before_linking: bool,
    word_cutoff_score_min: i32,
    _gap_open: i32,
    _gap_extend: i32,
    max_intron_length: i32,
) {
    use crate::link_hsps::{blast_link_hsp_list, LinkHSPParameters, LinkScoreBlock};
    use crate::program::TBLASTN;

    let link_params = LinkHSPParameters {
        gap_prob: crate::stat::BLAST_GAP_PROB_GAPPED,
        // TBLASTN's gapped translated-subject linker now feeds in per-HSP
        // evalues that have ALREADY been divided by 0.9 (sum-stats decay
        // applied at Spouge time). Use the standard gapped rate 0.1
        // (divisor 0.9) so the linker doesn't double-apply the decay; the
        // multi-HSP linker output then matches NCBI's empirical evalues.
        gap_decay_rate: crate::stat::BLAST_GAP_DECAY_RATE_GAPPED,
        longest_intron: translated_gapped_link_longest_intron(max_intron_length),
        cutoff_small_gap: word_cutoff_score_min,
        ..LinkHSPParameters::default()
    };

    for (result_index, result) in results.iter_mut().enumerate() {
        if result.hsps.is_empty() {
            continue;
        }
        let link_kbp = link_kbp_for_result(prot_kbp, result);
        let score_block = LinkScoreBlock {
            kbp: vec![link_kbp.clone()],
            kbp_gap: vec![link_kbp],
            gbp: gumbel_blk.cloned(),
            // NCBI setup stores `sbp->gbp->db_length = total_db_length / 3`
            // for translated-subject programs before effective lengths are
            // calculated (`blast_setup.c:913-919`).
            link_gbp_db_length: Some(
                (statistical_db_length / crate::stat::CODON_LENGTH as i64).max(1),
            ),
            recompute_evalues_before_uneven_linking: recompute_evalues_before_linking,
        };
        // The Rust `blast_link_hsps` port processes ALL HSPs including singletons —
        // for num=1 the uneven-gap formula `searchsp_eff * exp(-xsum)` is
        // still applied, overriding the pre-link Spouge evalue
        // (`link_hsps.c:1129` enters s_BlastUnevenGapLinkHSPs for tblastn).
        // The previous skip-singleton shortcut matched celegans by accident
        // (Spouge ≈ simple-Karlin for long alignments) but diverged for
        // short hits where FSC correction is non-trivial.
        let mut hsp_list = crate::hspstream::BlastHSPList {
            oid: result.subject_oid as i32,
            query_index: 0,
            hsp_array: result
                .hsps
                .iter()
                .map(|hsp| {
                    let (score, context, query, subject) = tblastn_link_key(hsp);
                    Some(api_hsp_to_blast_hsp(hsp, context, score, query, subject))
                })
                .collect(),
            hspcnt: result.hsps.len() as i32,
            allocated: result.hsps.len() as i32,
            hsp_max: i32::MAX,
            do_not_reallocate: false,
            best_evalue: result.best_evalue(),
        };

        blast_link_hsp_list(
            TBLASTN,
            &mut hsp_list,
            query_info,
            subject_stat_lengths
                .and_then(|lengths| lengths.get(result_index).copied())
                .unwrap_or_else(|| {
                    result
                        .hsps
                        .iter()
                        .map(|hsp| hsp.subject_link_end.max(hsp.subject_link_start) as i32)
                        .max()
                        .unwrap_or(result.subject_len as i32)
                }),
            &score_block,
            &link_params,
            true,
        );

        apply_blast_hsp_list_to_api_hsps(&mut result.hsps, hsp_list, tblastn_link_key);
    }
}

/// blast-rs: Public-parameter conversion for translated link-HSP setup; not a
/// direct NCBI C port.
fn translated_link_longest_intron(max_intron_length: i32) -> i32 {
    if max_intron_length <= 0 {
        0
    } else {
        max_intron_length.saturating_sub(2) / crate::stat::CODON_LENGTH as i32
    }
}

/// blast-rs: Gapped translated link-HSP intron default adapter; not a direct
/// NCBI C port.
fn translated_gapped_link_longest_intron(max_intron_length: i32) -> i32 {
    if max_intron_length == 0 {
        // NCBI treats translated gapped blastx/tblastn zero as "use the
        // default intron", not as even-gap linking (`blast_parameters.c:795`).
        (crate::stat::DEFAULT_LONGEST_INTRON as i32 - 2) / crate::stat::CODON_LENGTH as i32
    } else {
        translated_link_longest_intron(max_intron_length)
    }
}

/// blast-rs: SearchParams predicate for translated sum-statistics; not a
/// direct NCBI C port.
fn translated_sum_stats_enabled(params: &SearchParams) -> bool {
    params.sum_stats
        && (params.max_intron_length <= 0
            || translated_gapped_link_longest_intron(params.max_intron_length) > 0)
}

/// blast-rs: SearchParams adapter for preliminary composition e-value stretch;
/// not a direct NCBI C port.
fn composition_prelim_evalue(params: &SearchParams) -> f64 {
    if params.comp_adjust > 1 {
        params.evalue_threshold * 5.0
    } else {
        params.evalue_threshold
    }
}

/// blast-rs: SearchParams database-length override helper; not a direct NCBI C port.
fn statistical_db_length(params: &SearchParams, actual_db_length: i64) -> i64 {
    if params.db_length > 0 {
        params.db_length
    } else {
        actual_db_length
    }
}

/// blast-rs: Converts public SearchParams into setup options; not a direct
/// NCBI C port.
fn effective_lengths_options(params: &SearchParams) -> crate::options::EffectiveLengthsOptions {
    let mut options = crate::options::EffectiveLengthsOptions::default();
    if params.db_length > 0 {
        options.db_length = params.db_length;
    }
    if params.effective_search_space > 0 {
        options.num_searchspaces = 1;
        options.searchsp_eff = vec![params.effective_search_space];
    }
    options
}

fn blastn_search_hsp_link_key(
    hsp: &SearchHsp,
) -> (
    i32,
    i32,
    crate::hspstream::BlastSeg,
    crate::hspstream::BlastSeg,
) {
    (
        hsp.score,
        hsp.context,
        crate::hspstream::BlastSeg {
            frame: 1,
            offset: hsp.query_start,
            end: hsp.query_end,
            gapped_start: hsp.query_start,
        },
        crate::hspstream::BlastSeg {
            frame: if hsp.context == 1 { -1 } else { 1 },
            offset: hsp.subject_start,
            end: hsp.subject_end,
            gapped_start: hsp.subject_start,
        },
    )
}

/// blast-rs: blastn SearchHsp adapter into `BLAST_LinkHsps` over canonical `BlastHSPList`.
fn apply_blastn_linked_sum_stats_to_search_hsps(
    hsps: &mut Vec<SearchHsp>,
    query_len: i32,
    subject_len: i32,
    kbp_plus: &KarlinBlk,
    kbp_minus: &KarlinBlk,
    searchsp_plus: f64,
    searchsp_minus: f64,
    len_adj_plus: i32,
    len_adj_minus: i32,
) {
    use crate::link_hsps::{blast_link_hsp_list, LinkHSPParameters, LinkScoreBlock};
    use crate::program::BLASTN;
    use crate::queryinfo::QueryInfo;

    if hsps.len() <= 1 {
        return;
    }

    let mut query_info = QueryInfo::new_blastn(&[query_len.max(0) as usize]);
    if let Some(ctx) = query_info.contexts.get_mut(0) {
        ctx.eff_searchsp = searchsp_plus.round() as i64;
        ctx.length_adjustment = len_adj_plus;
    }
    if let Some(ctx) = query_info.contexts.get_mut(1) {
        ctx.eff_searchsp = searchsp_minus.round() as i64;
        ctx.length_adjustment = len_adj_minus;
    }

    let score_block = LinkScoreBlock {
        kbp: vec![kbp_plus.clone(), kbp_minus.clone()],
        kbp_gap: vec![kbp_plus.clone(), kbp_minus.clone()],
        ..LinkScoreBlock::default()
    };
    let link_params = LinkHSPParameters::default();
    let mut hsp_list = crate::hspstream::BlastHSPList {
        oid: 0,
        query_index: 0,
        hsp_array: hsps
            .iter()
            .map(|hsp| {
                let (score, context, query, subject) = blastn_search_hsp_link_key(hsp);
                Some(crate::hspstream::BlastHSP {
                    score,
                    num_ident: hsp.num_ident,
                    bit_score: hsp.bit_score,
                    evalue: hsp.evalue,
                    query,
                    subject,
                    context,
                    gap_info: None,
                    num: 1,
                    comp_adjustment_method: 0,
                    pat_info: None,
                    num_positives: 0,
                    map_info: None,
                })
            })
            .collect(),
        hspcnt: hsps.len() as i32,
        allocated: hsps.len() as i32,
        hsp_max: i32::MAX,
        do_not_reallocate: false,
        best_evalue: f64::INFINITY,
    };

    blast_link_hsp_list(
        BLASTN,
        &mut hsp_list,
        &query_info,
        subject_len,
        &score_block,
        &link_params,
        false,
    );

    let mut remaining: Vec<Option<SearchHsp>> =
        std::mem::take(hsps).into_iter().map(Some).collect();
    for linked in hsp_list.hsp_array.into_iter().flatten() {
        let mut selected = None;
        for (idx, candidate) in remaining.iter().enumerate() {
            let Some(candidate) = candidate.as_ref() else {
                continue;
            };
            let (score, context, query, subject) = blastn_search_hsp_link_key(candidate);
            if linked_api_hsp_key(context, score, &query, &subject)
                == (
                    linked.score,
                    linked.context,
                    linked.query.frame,
                    linked.query.offset,
                    linked.query.end,
                    linked.query.gapped_start,
                    linked.subject.frame,
                    linked.subject.offset,
                    linked.subject.end,
                    linked.subject.gapped_start,
                )
            {
                selected = Some(idx);
                break;
            }
        }
        if let Some(idx) = selected {
            if let Some(mut hsp) = remaining[idx].take() {
                hsp.evalue = linked.evalue;
                hsps.push(hsp);
            }
        }
    }
}

/// blast-rs: Aggregates API HSPs by subject OID; not a direct NCBI C port.
fn push_hsp_for_subject(
    results: &mut [Option<SearchResult>],
    oid: u32,
    title: &str,
    accession: &str,
    subject_len: usize,
    taxids: &[i32],
    hsp: Hsp,
) {
    let slot = &mut results[oid as usize];
    match slot {
        Some(existing) => existing.hsps.push(hsp),
        None => {
            *slot = Some(SearchResult {
                subject_oid: oid,
                subject_title: title.to_string(),
                subject_accession: accession.to_string(),
                subject_len,
                hsps: vec![hsp],
                taxids: taxids.to_vec(),
            });
        }
    }
}

/// blast-rs: Public-HSP e-value/score ordering adapter; not a direct NCBI C port.
fn compare_hsps_by_evalue_then_score(a: &Hsp, b: &Hsp) -> std::cmp::Ordering {
    crate::hspstream::evalue_comp(a.evalue, b.evalue).then_with(|| compare_hsps_by_score(a, b))
}

/// blast-rs: Public-HSP score ordering adapter; not a direct NCBI C port.
fn compare_hsps_by_score(a: &Hsp, b: &Hsp) -> std::cmp::Ordering {
    let a_subject_offset = a.subject_start.min(a.subject_end).saturating_sub(1);
    let b_subject_offset = b.subject_start.min(b.subject_end).saturating_sub(1);
    let a_subject_end = a.subject_start.max(a.subject_end);
    let b_subject_end = b.subject_start.max(b.subject_end);
    let a_query_offset = a.query_start.min(a.query_end).saturating_sub(1);
    let b_query_offset = b.query_start.min(b.query_end).saturating_sub(1);
    let a_query_end = a.query_start.max(a.query_end);
    let b_query_end = b.query_start.max(b.query_end);

    b.score
        .cmp(&a.score)
        .then_with(|| a_subject_offset.cmp(&b_subject_offset))
        .then_with(|| b_subject_end.cmp(&a_subject_end))
        .then_with(|| a_query_offset.cmp(&b_query_offset))
        .then_with(|| b_query_end.cmp(&a_query_end))
}

/// blast-rs: Protein cutoff adapter over Karlin/Spouge helpers; not a direct
/// NCBI C port.
fn protein_eval_cutoff(
    evalue_threshold: f64,
    prot_kbp: &crate::stat::KarlinBlk,
    gumbel_blk: Option<&crate::stat::GumbelBlk>,
    query_length: i32,
    min_subject_length: i32,
    search_space: f64,
) -> i32 {
    if let Some(gbp) = gumbel_blk {
        crate::stat::blast_spouge_eto_s(
            evalue_threshold,
            Some(prot_kbp),
            Some(gbp),
            query_length.max(1),
            min_subject_length.max(1),
        )
        .max(1)
    } else {
        prot_kbp
            .evalue_to_raw(evalue_threshold, search_space.max(1.0))
            .max(1)
    }
}

/// blast-rs: Preliminary protein seed cutoff adapter; not a direct NCBI C port.
fn protein_prelim_seed_cutoff(
    gap_trigger_raw: i32,
    evalue_threshold: f64,
    prot_kbp: &crate::stat::KarlinBlk,
    gumbel_blk: Option<&crate::stat::GumbelBlk>,
    query_length: i32,
    min_subject_length: i32,
    search_space: f64,
) -> i32 {
    let eval_cutoff = protein_eval_cutoff(
        evalue_threshold,
        prot_kbp,
        gumbel_blk,
        query_length,
        min_subject_length,
        search_space,
    );
    gap_trigger_raw.min(eval_cutoff).max(1)
}

/// NCBI: translated gapped sum-statistics relaxes `hit_params->cutoffs[context]`
/// after the normal preliminary cutoff is computed (`blast_parameters.c:951`).
/// `cutoff_score_max` is not relaxed, so callers should keep word-finder
/// cutoffs separate from this final hit-saving cutoff.
fn protein_sum_stats_hit_cutoff(
    normal_hit_cutoff: i32,
    prot_kbp: &crate::stat::KarlinBlk,
    avg_query_length: i32,
    avg_subject_length: i32,
) -> i32 {
    let searchsp = avg_query_length.max(1).min(avg_subject_length.max(1)) as f64
        * avg_subject_length.max(1) as f64;
    let (relaxed_cutoff, _) = crate::stat::blast_cutoffs(
        1,
        1.0,
        prot_kbp,
        searchsp.max(1.0),
        true,
        crate::stat::BLAST_GAP_DECAY_RATE_GAPPED,
    );
    normal_hit_cutoff.min(relaxed_cutoff).max(1)
}

fn query_info_average_context_length(query_info: &crate::queryinfo::QueryInfo) -> i32 {
    let Some(last) = query_info.contexts.last() else {
        return 1;
    };
    let denom = query_info.contexts.len().max(1) as i32;
    ((last.query_offset + last.query_length) / denom).max(1)
}

fn tblastx_initial_word_cutoff(
    gap_trigger_raw: i32,
    ctx_kbp: &crate::stat::KarlinBlk,
    query_length: i32,
    subject_length: i32,
    hit_cutoff: i32,
) -> i32 {
    let q = query_length.max(1) as u64;
    let s = subject_length.max(1) as u64;
    let searchsp = s.saturating_mul(s.min(q)) as f64;
    let (initial_cutoff, _) = crate::stat::blast_cutoffs(
        1,
        crate::stat::CUTOFF_E_TBLASTX,
        ctx_kbp,
        searchsp.max(1.0),
        true,
        crate::stat::BLAST_GAP_DECAY_RATE,
    );
    gap_trigger_raw.min(initial_cutoff).min(hit_cutoff).max(1)
}

/// blast-rs: Converts public bit gap-trigger to raw score; not a direct NCBI C port.
fn protein_gap_trigger_raw(gap_trigger_bits: f64, ungapped_kbp: &crate::stat::KarlinBlk) -> i32 {
    raw_gap_trigger_bits(gap_trigger_bits, ungapped_kbp)
}

/// NCBI `BlastInitialWordParametersNew`: ungapped x-drop is converted with
/// `ceil(bits * ln(2) / Lambda)` before the `Int4` cast.
fn raw_ungapped_xdrop_bits(bits: i32, kbp: &crate::stat::KarlinBlk) -> i32 {
    (bits as f64 * crate::math::NCBIMATH_LN2 / kbp.lambda).ceil() as i32
}

/// NCBI `BlastScoringParametersNew`: gapped x-drop values use a plain `Int4`
/// cast, i.e. truncation toward zero.
fn raw_gapped_xdrop_bits(bits: i32, kbp: &crate::stat::KarlinBlk) -> i32 {
    raw_gapped_xdrop_bits_f64(bits as f64, kbp)
}

fn raw_gapped_xdrop_bits_f64(bits: f64, kbp: &crate::stat::KarlinBlk) -> i32 {
    (bits * crate::math::NCBIMATH_LN2 / kbp.lambda) as i32
}

/// NCBI `BlastInitialWordParametersUpdate`: gap trigger uses the ungapped KBP
/// and a plain `Int4` cast.
fn raw_gap_trigger_bits(bits: f64, ungapped_kbp: &crate::stat::KarlinBlk) -> i32 {
    ((bits * crate::math::NCBIMATH_LN2 + ungapped_kbp.log_k) / ungapped_kbp.lambda) as i32
}

/// blast-rs: Kappa redo near-identical predicate adapter; not a direct NCBI C port.
fn kappa_redo_near_identical(
    ph: &crate::protein_lookup::ProteinHit,
    query_len: usize,
    subject_len: usize,
    gapped_lambda: f64,
) -> bool {
    const NEAR_IDENTICAL_BITS_PER_POSITION: f64 = 1.74;
    const MINIMUM_LENGTH_NEAR_IDENTICAL: usize = 50;

    if gapped_lambda <= 0.0 {
        return false;
    }
    let query_span = ph.query_end.saturating_sub(ph.query_start);
    let subject_span = ph.subject_end.saturating_sub(ph.subject_start);
    // C `s_preliminaryTestNearIdentical` (redo_alignment.c:1066) tests
    // `(matchEnd - matchStart + 1) < MIN(queryLength, MIN_LEN)`. `subject_span`
    // here is `matchEnd - matchStart` (no +1), so add 1 to match C exactly.
    if subject_span + 1 < query_len.min(MINIMUM_LENGTH_NEAR_IDENTICAL) {
        return false;
    }
    let align_len = query_span.min(subject_span);
    if align_len == 0 || subject_len == 0 {
        return false;
    }
    let cutoff = (NEAR_IDENTICAL_BITS_PER_POSITION * crate::math::NCBIMATH_LN2) / gapped_lambda;
    (ph.score as f64 / align_len as f64) >= cutoff
}

/// blast-rs: Kappa redo subject SEG-mask helper; not a direct NCBI C port.
fn kappa_seg_mask_subject_for_redo_into(subject: &[u8], buf: &mut Vec<u8>) {
    buf.clear();
    buf.extend_from_slice(subject);
    let mask = crate::filter::seg_filter_ncbistdaa(buf, 10, 1.8, 2.1);
    for region in &mask.regions {
        let start = region.start.max(0) as usize;
        let end = (region.end as usize).min(buf.len());
        for aa in &mut buf[start..end] {
            *aa = NCBISTDAA_X;
        }
    }
}

/// blast-rs: Kappa redo stage-2 near-identical predicate adapter; not a direct
/// NCBI C port.
///
/// NCBI's `s_SequenceGetProteinRange` (blast_kappa.c:1629-1637) skips SEG
/// masking of the subject only when BOTH the cheap stage-1 density test
/// (`s_preliminaryTestNearIdentical`) AND the expensive stage-2 fraction-
/// identical test (`s_TestNearIdentical`, true identity > 0.95) pass:
/// `if (shouldTestIdentical && !s_TestNearIdentical(...)) { SEG the subject; }`.
/// Stage-1 alone passes for ~94% paralogs, so omitting stage-2 wrongly skips
/// SEG and over-scores the composition redo by ~20 bits. This wraps the faithful
/// stage-2 port `blast_kappa::test_near_identical` for the production path,
/// building the `BlastCompo*` views from the HSP coordinates (query/subject
/// `*_end` are one-past, matching `BlastCompoAlignment::query_end/match_end`).
fn kappa_redo_stage2_near_identical(
    ph: &crate::protein_lookup::ProteinHit,
    query: &[u8],
    subject: &[u8],
) -> bool {
    // 8-mer rolling hashes over the full query (MAX_SHIFT/word_size = 8 inside
    // s_FindNumIdentical / s_TestNearIdentical).
    const WORD_SIZE: usize = 8;
    if query.len() < WORD_SIZE {
        // Too short for the k-mer sweep; treat as "not near-identical" so SEG
        // is applied (conservative, matches NCBI which still runs the extend
        // tests but cannot satisfy the >0.95 fraction on a tiny range).
        return false;
    }
    let query_words: Vec<u64> = (0..=query.len() - WORD_SIZE)
        .map(|i| crate::blast_kappa::s_get_hash(&query[i..i + WORD_SIZE], WORD_SIZE))
        .collect();
    let query_data = crate::blast_kappa::BlastCompoSequenceData {
        buffer: query.to_vec(),
        data_offset: 0,
        length: query.len() as i32,
    };
    let seq_data = crate::blast_kappa::BlastCompoSequenceData {
        buffer: subject.to_vec(),
        data_offset: 0,
        length: subject.len() as i32,
    };
    let align = crate::blast_kappa::BlastCompoAlignment::new(
        ph.score,
        crate::compo_mode_condition::MatrixAdjustRule::DontAdjust,
        0,
        ph.query_start as i32,
        ph.query_end as i32, // one-past
        ph.subject_start as i32,
        ph.subject_end as i32, // one-past
        0,
        None,
    );
    crate::blast_kappa::test_near_identical(&seq_data, 0, &query_data, 0, &query_words, &align)
}

/// blast-rs: Selects subject sequence for kappa redo alignment; not a direct
/// NCBI C port.
fn kappa_redo_subject_sequence<'a>(
    query: &[u8],
    subject: &'a [u8],
    phits: &[crate::protein_lookup::ProteinHit],
    gapped_lambda: f64,
    mask_buf: &'a mut Vec<u8>,
) -> &'a [u8] {
    // NCBI gates SEG-skipping on BOTH stages: the cheap stage-1 density test
    // AND the expensive stage-2 true-fraction-identical test (> 0.95). Only
    // skip SEG when BOTH pass; otherwise SEG the subject for the redo.
    let near_identical = phits
        .first()
        .map(|ph| {
            kappa_redo_near_identical(ph, query.len(), subject.len(), gapped_lambda)
                && kappa_redo_stage2_near_identical(ph, query, subject)
        })
        .unwrap_or(false);
    if near_identical {
        subject
    } else {
        kappa_seg_mask_subject_for_redo_into(subject, mask_buf);
        mask_buf.as_slice()
    }
}

/// blast-rs: ProteinHit score ordering adapter; not a direct NCBI C port.
fn compare_protein_hits_by_score(
    a: &crate::protein_lookup::ProteinHit,
    b: &crate::protein_lookup::ProteinHit,
) -> std::cmp::Ordering {
    b.score
        .cmp(&a.score)
        .then_with(|| a.subject_start.cmp(&b.subject_start))
        .then_with(|| b.subject_end.cmp(&a.subject_end))
        .then_with(|| a.query_start.cmp(&b.query_start))
        .then_with(|| b.query_end.cmp(&a.query_end))
}

/// NCBI: s_ChainingAlignment (blast_gapalign.c:3592).
/// naming: Public Rust wrapper keeps the C static symbol name in snake_case.
pub fn s_chaining_alignment(
    ungapped_hits: &[crate::protein_lookup::ProteinHit],
    gap_open: i32,
    gap_extend: i32,
    word_cutoff: i32,
    hit_cutoff: i32,
) -> Vec<bool> {
    #[derive(Clone, Copy)]
    struct ChainingNode<'a> {
        original_index: usize,
        hit: &'a crate::protein_lookup::ProteinHit,
    }

    let mut chained_hits: Vec<ChainingNode<'_>> = ungapped_hits
        .iter()
        .enumerate()
        .map(|(original_index, hit)| ChainingNode {
            original_index,
            hit,
        })
        .collect();
    chained_hits.sort_by(|a, b| {
        a.hit
            .query_start
            .cmp(&b.hit.query_start)
            .then_with(|| b.hit.score.cmp(&a.hit.score))
            .then_with(|| a.original_index.cmp(&b.original_index))
    });
    let gap_score = gap_open + gap_extend;
    let n = chained_hits.len();
    let mut best_score: Vec<i32> = chained_hits.iter().map(|node| node.hit.score).collect();
    let mut kept = vec![true; ungapped_hits.len()];

    for k in (0..n).rev() {
        let self_score = best_score[k];
        for j in (k + 1)..n {
            let q_diff = chained_hits[j].hit.query_start as i32
                - chained_hits[k].hit.query_start as i32
                + chained_hits[k]
                    .hit
                    .query_end
                    .saturating_sub(chained_hits[k].hit.query_start) as i32;
            let s_diff = chained_hits[j].hit.subject_start as i32
                - chained_hits[k].hit.subject_start as i32
                + chained_hits[k]
                    .hit
                    .subject_end
                    .saturating_sub(chained_hits[k].hit.subject_start) as i32;
            if s_diff < 0 {
                continue;
            }
            let bridge = (q_diff.min(s_diff) * 3).min(word_cutoff);
            let gap_penalty = q_diff.abs_diff(s_diff).max(1) as i32 + gap_open;
            let new_score = self_score + best_score[j] + bridge - gap_penalty;
            if new_score > best_score[k] {
                best_score[k] = new_score;
            }
        }
    }

    for (node, &score) in chained_hits.iter().zip(best_score.iter()) {
        if score - gap_score + word_cutoff - 1 < hit_cutoff {
            kept[node.original_index] = false;
        }
    }
    kept
}

/// NCBI: s_BlastProtGappedAlignment (`blast_gapalign.c:4298`).
///
/// Rust's protein path stores the preliminary result directly as a
/// `ProteinHit`, but the original callgraph has this wrapper between
/// `BLAST_GetGappedScore` and `Blast_SemiGappedAlign`.
fn s_blast_prot_gapped_alignment(
    query_aa: &[u8],
    subj_aa: &[u8],
    seed_q: usize,
    seed_s: usize,
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
    gap_open: i32,
    gap_extend: i32,
    x_drop_gapped: i32,
) -> Option<crate::protein::ProteinGappedResult> {
    crate::protein::protein_gapped_align(
        query_aa,
        subj_aa,
        seed_q,
        seed_s,
        matrix,
        gap_open,
        gap_extend,
        x_drop_gapped,
    )
}

/// NCBI: BLAST_GappedAlignmentWithTraceback (`blast_gapalign.c`).
///
/// The Rust aligner returns the edit script together with the score, so this is
/// a callgraph-preserving wrapper over the same DP used by traceback.
fn blast_gapped_alignment_with_traceback(
    query_aa: &[u8],
    subj_aa: &[u8],
    seed_q: usize,
    seed_s: usize,
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
    gap_open: i32,
    gap_extend: i32,
    x_drop_final: i32,
) -> Option<crate::protein::ProteinGappedResult> {
    crate::protein::protein_gapped_align(
        query_aa,
        subj_aa,
        seed_q,
        seed_s,
        matrix,
        gap_open,
        gap_extend,
        x_drop_final,
    )
}

fn protein_hit_query_interval(hit: &crate::protein_lookup::ProteinHit) -> crate::itree::Interval {
    crate::itree::Interval::new(hit.query_start as i32, hit.query_end as i32)
}

fn protein_hit_subject_interval(hit: &crate::protein_lookup::ProteinHit) -> crate::itree::Interval {
    crate::itree::Interval::new(hit.subject_start as i32, hit.subject_end as i32)
}

/// NCBI: BlastIntervalTreeContainsHSP (`blast_itree.c`), preserving the
/// query-index and subject-frame metadata used by translated searches.
fn blast_interval_tree_contains_protein_hit_with_metadata(
    tree: &crate::itree::IntervalTree,
    hit: &crate::protein_lookup::ProteinHit,
    query_index: i32,
    subject_frame: i32,
) -> bool {
    tree.is_contained_with_metadata(
        protein_hit_query_interval(hit),
        protein_hit_subject_interval(hit),
        hit.score,
        query_index,
        subject_frame,
    )
}

/// NCBI: BlastIntervalTreeAddHSP (`blast_itree.c:511`), preserving translated
/// search metadata.
fn blast_interval_tree_add_protein_hit_with_metadata(
    tree: &mut crate::itree::IntervalTree,
    hit: &crate::protein_lookup::ProteinHit,
    query_index: i32,
    subject_frame: i32,
) {
    tree.insert_with_metadata(
        protein_hit_query_interval(hit),
        protein_hit_subject_interval(hit),
        hit.score,
        query_index,
        subject_frame,
    );
}

/// NCBI: Blast_HSPInit (`blast_hits.c:151`), specialized to the local
/// `ProteinHit` representation.
fn blast_hsp_init_protein_hit(
    gapped: &crate::protein::ProteinGappedResult,
    gapped_start_q: usize,
    gapped_start_s: usize,
    _needs_traceback: bool,
) -> crate::protein_lookup::ProteinHit {
    crate::protein_lookup::ProteinHit {
        query_start: gapped.query_start,
        query_end: gapped.query_end,
        subject_start: gapped.subject_start,
        subject_end: gapped.subject_end,
        score: gapped.score,
        num_ident: gapped.num_ident,
        align_length: gapped.align_length,
        mismatches: gapped.mismatches,
        gap_opens: gapped.gap_opens,
        qseq: None,
        sseq: None,
        scaled_score: None,
        adjusted_evalue: None,
        gapped_start_q,
        gapped_start_s,
    }
}

/// NCBI: Blast_HSPGetNumIdentitiesAndPositives (`blast_hits.c`).
fn blast_hsp_get_num_identities_and_positives_protein_hit(
    gapped: &crate::protein::ProteinGappedResult,
) -> (i32, i32, i32, i32) {
    (
        gapped.num_ident,
        gapped.align_length,
        gapped.mismatches,
        gapped.gap_opens,
    )
}

/// NCBI: Blast_HSPUpdateWithTraceback (`blast_traceback.c:78`), specialized to
/// the local `ProteinHit` representation.
fn blast_hsp_update_with_traceback_protein_hit(
    hit: &mut crate::protein_lookup::ProteinHit,
    gapped: crate::protein::ProteinGappedResult,
    qseq: String,
    sseq: String,
) {
    let (num_ident, align_length, mismatches, gap_opens) =
        blast_hsp_get_num_identities_and_positives_protein_hit(&gapped);
    hit.query_start = gapped.query_start;
    hit.query_end = gapped.query_end;
    hit.subject_start = gapped.subject_start;
    hit.subject_end = gapped.subject_end;
    hit.score = gapped.score;
    hit.num_ident = num_ident;
    hit.align_length = align_length;
    hit.mismatches = mismatches;
    hit.gap_opens = gap_opens;
    hit.qseq = Some(qseq);
    hit.sseq = Some(sseq);
}

/// NCBI: Blast_HSPListPurgeHSPsWithCommonEndpoints (`blast_hits.c`).
fn blast_hsp_list_purge_hsps_with_common_endpoints_protein_hits(
    hits: &mut Vec<crate::protein_lookup::ProteinHit>,
) {
    purge_hsps_with_common_endpoints(hits);
}

/// blast-rs: API-local preliminary protein gapping over `ProteinHit` values.
///
/// The C-shaped `BLAST_GetGappedScore` translation lives in `crate::extend`;
/// this helper remains only while the public translated-search paths still
/// carry preliminary hits as `ProteinHit` rather than `BlastHSPList`.
fn protein_prelim_gapped_hits(
    query_aa: &[u8],
    subj_aa: &[u8],
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
    ungapped_hits: &[crate::protein_lookup::ProteinHit],
    gap_open: i32,
    gap_extend: i32,
    x_drop_gapped: i32,
    _x_drop_final: i32,
    word_cutoff: i32,
    hit_cutoff: i32,
    chaining_enabled: bool,
    scratch: &mut ProteinScratch,
    reset_tree: bool,
    query_index: i32,
    subject_frame: i32,
) -> Vec<crate::protein_lookup::ProteinHit> {
    let mut hits = Vec::new();
    // Mirror NCBI's `BLAST_GetGappedScore` flow: sort the init-hit list inside
    // this function (`Blast_InitHitListSortByScore`). NCBI only
    // runs `s_ChainingAlignment` when protein chaining is enabled and the
    // matrix is BLOSUM62 (`blast_gapalign.c:3793`); otherwise every saved init
    // HSP remains eligible for the preliminary gapped DP.
    let mut sorted_hits = ungapped_hits.to_vec();
    sorted_hits.sort_by(compare_protein_hits_by_score);
    let chained_kept = if chaining_enabled {
        s_chaining_alignment(
            &sorted_hits,
            gap_open,
            gap_extend,
            word_cutoff.max(1),
            hit_cutoff.max(1),
        )
    } else {
        vec![true; sorted_hits.len()]
    };
    if reset_tree {
        scratch
            .prelim_tree
            .reset_for_query(query_aa.len() as i32 + 1, subj_aa.len() as i32 + 1);
    }
    let tree = &mut scratch.prelim_tree;
    for (seed, keep) in sorted_hits.iter().zip(chained_kept.iter()) {
        if seed.score < word_cutoff.max(1) || !*keep {
            continue;
        }
        if blast_interval_tree_contains_protein_hit_with_metadata(
            tree,
            seed,
            query_index,
            subject_frame,
        ) {
            continue;
        }
        let (seed_q, seed_s) = crate::protein::blast_get_start_for_gapped_alignment(
            query_aa,
            subj_aa,
            seed.query_start,
            seed.query_end.saturating_sub(seed.query_start),
            seed.subject_start,
            seed.subject_end.saturating_sub(seed.subject_start),
            matrix,
        );
        // PRELIMINARY gapped DP only (matches NCBI engine flow:
        // `s_BlastProtGappedAlignment` calls `Blast_SemiGappedAlign` with
        // preliminary `gap_x_dropoff` only — the larger
        // `gap_x_dropoff_final` is used in `Blast_TracebackFromHSPList`).
        // Pre/post-gapped containment uses preliminary bounds (interval
        // tree carries preliminary HSPs only). Final-xdrop traceback
        // happens in the post-loop block below.
        let Some(prelim) = s_blast_prot_gapped_alignment(
            query_aa,
            subj_aa,
            seed_q,
            seed_s,
            matrix,
            gap_open,
            gap_extend,
            x_drop_gapped,
        ) else {
            continue;
        };
        if prelim.score < hit_cutoff.max(1) {
            continue;
        }
        let prelim_hit = blast_hsp_init_protein_hit(&prelim, seed_q, seed_s, true);
        if blast_interval_tree_contains_protein_hit_with_metadata(
            tree,
            &prelim_hit,
            query_index,
            subject_frame,
        ) {
            continue;
        }
        blast_interval_tree_add_protein_hit_with_metadata(
            tree,
            &prelim_hit,
            query_index,
            subject_frame,
        );
        hits.push(prelim_hit);
    }
    hits
}

/// NCBI: Blast_TracebackFromHSPList (blast_traceback.c:345).
pub fn blast_traceback_from_hsp_list(
    program_number: crate::program::ProgramType,
    hsp_list: &mut crate::hspstream::BlastHSPList,
    query_blk: &crate::util::BlastSequenceBlk,
    subject_blk: &mut crate::util::BlastSequenceBlk,
    query_info: &crate::queryinfo::QueryInfo,
    gap_align: &mut crate::blast_kappa::BlastGapAlignWorkspace,
    sbp: &crate::stat::BlastScoreBlk,
    score_params: &crate::parameters::ScoringParameters,
    ext_options: &crate::options::ExtensionOptions,
    hit_params: &crate::parameters::HitSavingParameters,
    gen_code_string: Option<&[u8; 64]>,
    fence_hit: Option<&mut bool>,
    scratch: &mut ProteinScratch,
    stat_length_out: Option<&mut i32>,
) -> i32 {
    if hsp_list.hspcnt == 0 {
        return 0;
    }

    let _ = (query_info, ext_options, gen_code_string, fence_hit);
    let stat_length = subject_blk.length;

    let Some(query) = query_blk.sequence.as_deref() else {
        return -1;
    };
    let Some(subject) = subject_blk.sequence.as_deref() else {
        return -1;
    };
    if query_blk.length < 0
        || subject_blk.length < 0
        || query_blk.length as usize > query.len()
        || subject_blk.length as usize > subject.len()
    {
        return -1;
    }
    if sbp.matrix.data.len() < AA_SIZE || sbp.matrix.data.iter().any(|row| row.len() < AA_SIZE) {
        return -1;
    }

    let mut matrix = [[0i32; AA_SIZE]; AA_SIZE];
    for (row_index, row) in matrix.iter_mut().enumerate() {
        row.copy_from_slice(&sbp.matrix.data[row_index][..AA_SIZE]);
    }

    if crate::program::blast_subject_is_translated(program_number) {
        let gen_code = gen_code_string
            .copied()
            .or(subject_blk.gen_code_string)
            .unwrap_or(crate::util::STANDARD_GENETIC_CODE);
        let mut target_t = None;
        if crate::util::blast_target_translation_new(
            subject_blk,
            &gen_code,
            program_number,
            false,
            &mut target_t,
        ) != 0
        {
            return -1;
        }
        let Some(mut target_t) = target_t else {
            return -1;
        };

        hsp_list
            .hsp_array
            .sort_by(|a, b| match (a.as_ref(), b.as_ref()) {
                (Some(a), Some(b)) => compare_blast_hsps_by_score(a, b),
                (Some(_), None) => std::cmp::Ordering::Less,
                (None, Some(_)) => std::cmp::Ordering::Greater,
                (None, None) => std::cmp::Ordering::Equal,
            });
        scratch
            .tb_tree
            .reset_for_query(query_blk.length + 1, (subject_blk.length / 3).max(1) + 1);
        let tb_tree = &mut scratch.tb_tree;
        let mut keep = vec![true; hsp_list.hsp_array.len()];
        for (idx, hsp_slot) in hsp_list.hsp_array.iter_mut().enumerate() {
            let Some(hsp) = hsp_slot.as_mut() else {
                keep[idx] = false;
                continue;
            };
            let mut hit = blast_hsp_to_protein_hit(hsp);
            let subject_frame = hsp.subject.frame as i32;
            if blast_interval_tree_contains_protein_hit_with_metadata(
                tb_tree,
                &hit,
                hsp_list.query_index,
                subject_frame,
            ) {
                keep[idx] = false;
                continue;
            }

            let legacy_hsp = hsp.clone().into_legacy_hsp();
            let Some(view) = crate::hspstream::blast_hsp_get_target_translation(
                &mut target_t,
                Some(&legacy_hsp),
            ) else {
                keep[idx] = false;
                continue;
            };
            let translation_start = (1isize - view.pointer_offset).max(0) as i32;
            let translation_stop = view.translated_length.max(translation_start);
            let mut translated_subject =
                Vec::with_capacity((translation_stop - translation_start) as usize);
            for offset in translation_start..translation_stop {
                let Some(aa) = view.get(offset) else {
                    translated_subject.clear();
                    break;
                };
                translated_subject.push(aa);
            }
            if translated_subject.is_empty() {
                keep[idx] = false;
                continue;
            }

            let local_seed_s = hit.gapped_start_s as i32 - translation_start;
            if local_seed_s < 0 {
                keep[idx] = false;
                continue;
            }
            let Some(gr) = blast_gapped_alignment_with_traceback(
                &query[..query_blk.length as usize],
                &translated_subject,
                hit.gapped_start_q,
                local_seed_s as usize,
                &matrix,
                score_params.gap_open,
                score_params.gap_extend,
                gap_align.gap_x_dropoff,
            ) else {
                keep[idx] = false;
                continue;
            };

            let q_slice = &query[gr.query_start..gr.query_end];
            let s_slice = &translated_subject[gr.subject_start..gr.subject_end];
            let edit_script = gr.edit_script.clone();
            let (qseq, sseq) =
                gr.edit_script
                    .render_alignment(q_slice, s_slice, ncbistdaa_to_aminoacid_char);
            blast_hsp_update_with_traceback_protein_hit(&mut hit, gr, qseq, sseq);
            blast_hsp_update_from_protein_hit(hsp, hit.clone());
            hsp.gap_info = Some(edit_script);
            hsp.subject.frame = subject_frame as i16;
            s_compute_num_identities_blast_hsp(
                hsp,
                &query[..query_blk.length as usize],
                &translated_subject,
                Some(&matrix),
            );
            hit.subject_start = hit
                .subject_start
                .saturating_add(translation_start.max(0) as usize);
            hit.subject_end = hit
                .subject_end
                .saturating_add(translation_start.max(0) as usize);
            hit.gapped_start_s = hit
                .gapped_start_s
                .saturating_add(translation_start.max(0) as usize);
            blast_hsp_update_from_protein_hit(hsp, hit);
            hsp.subject.frame = subject_frame as i16;
            if s_hsp_test_blast_hsp(
                hsp,
                &hit_params.options,
                hsp.gap_info
                    .as_ref()
                    .map(|script| script.alignment_length())
                    .unwrap_or_else(|| (hsp.query.end - hsp.query.offset).max(0)),
            ) {
                keep[idx] = false;
                continue;
            }
            let accepted_hit = blast_hsp_to_protein_hit(hsp);
            blast_interval_tree_add_protein_hit_with_metadata(
                tb_tree,
                &accepted_hit,
                hsp_list.query_index,
                subject_frame,
            );
        }

        let mut idx = 0usize;
        hsp_list.hsp_array.retain(|_| {
            let k = keep[idx];
            idx += 1;
            k
        });
        hsp_list.hspcnt = hsp_list.hsp_array.len() as i32;
        hsp_list.allocated = hsp_list.hsp_array.len() as i32;
        crate::hspstream::blast_hsp_list_purge_blast_hsps_with_common_endpoints(
            program_number,
            Some(hsp_list),
            true,
        );
        hsp_list.best_evalue = hsp_list
            .hsp_array
            .iter()
            .filter_map(|hsp| hsp.as_ref().map(|hsp| hsp.evalue))
            .fold(i32::MAX as f64, f64::min);
        if let Some(out) = stat_length_out {
            *out = stat_length;
        }
        return 0;
    }

    if crate::program::blast_query_is_translated(program_number)
        || program_number == crate::program::RPS_TBLASTN
    {
        return -1;
    }

    sort_blast_hsp_list_by_score(hsp_list);
    scratch
        .tb_tree
        .reset_for_query(query_blk.length + 1, subject_blk.length + 1);
    let tb_tree = &mut scratch.tb_tree;
    let mut keep = vec![true; hsp_list.hsp_array.len()];
    for (idx, hsp_slot) in hsp_list.hsp_array.iter_mut().enumerate() {
        let Some(hsp) = hsp_slot.as_mut() else {
            keep[idx] = false;
            continue;
        };
        let mut hit = blast_hsp_to_protein_hit(hsp);
        let subject_frame = hsp.subject.frame as i32;
        if blast_interval_tree_contains_protein_hit_with_metadata(
            tb_tree,
            &hit,
            hsp_list.query_index,
            subject_frame,
        ) {
            keep[idx] = false;
            continue;
        }
        let Some(gr) = blast_gapped_alignment_with_traceback(
            &query[..query_blk.length as usize],
            &subject[..subject_blk.length as usize],
            hit.gapped_start_q,
            hit.gapped_start_s,
            &matrix,
            score_params.gap_open,
            score_params.gap_extend,
            gap_align.gap_x_dropoff,
        ) else {
            keep[idx] = false;
            continue;
        };
        let q_slice = &query[gr.query_start..gr.query_end];
        let s_slice = &subject[gr.subject_start..gr.subject_end];
        let edit_script = gr.edit_script.clone();
        let (qseq, sseq) =
            gr.edit_script
                .render_alignment(q_slice, s_slice, ncbistdaa_to_aminoacid_char);
        blast_hsp_update_with_traceback_protein_hit(&mut hit, gr, qseq, sseq);
        blast_hsp_update_from_protein_hit(hsp, hit);
        hsp.gap_info = Some(edit_script);
        s_compute_num_identities_blast_hsp(
            hsp,
            &query[..query_blk.length as usize],
            &subject[..subject_blk.length as usize],
            Some(&matrix),
        );
        if s_hsp_test_blast_hsp(
            hsp,
            &hit_params.options,
            hsp.gap_info
                .as_ref()
                .map(|script| script.alignment_length())
                .unwrap_or_else(|| (hsp.query.end - hsp.query.offset).max(0)),
        ) {
            keep[idx] = false;
            continue;
        }
        let accepted_hit = blast_hsp_to_protein_hit(hsp);
        blast_interval_tree_add_protein_hit_with_metadata(
            tb_tree,
            &accepted_hit,
            hsp_list.query_index,
            subject_frame,
        );
    }
    let mut idx = 0usize;
    hsp_list.hsp_array.retain(|_| {
        let k = keep[idx];
        idx += 1;
        k
    });
    hsp_list.hspcnt = hsp_list.hsp_array.len() as i32;
    hsp_list.allocated = hsp_list.hspcnt;
    crate::hspstream::blast_hsp_list_purge_blast_hsps_with_common_endpoints(
        program_number,
        Some(hsp_list),
        true,
    );
    sort_blast_hsp_list_by_score(hsp_list);
    hsp_list.hspcnt = hsp_list.hsp_array.len() as i32;
    hsp_list.allocated = hsp_list.hsp_array.len() as i32;
    hsp_list.best_evalue = hsp_list
        .hsp_array
        .iter()
        .filter_map(|hsp| hsp.as_ref().map(|hsp| hsp.evalue))
        .fold(i32::MAX as f64, f64::min);
    if let Some(out) = stat_length_out {
        *out = stat_length;
    }

    0
}

/// NCBI: Blast_TracebackFromHSPList (blast_traceback.c:345), specialized to
/// the local `ProteinHit` representation used by ordinary protein searches.
fn blast_traceback_from_protein_hit_hsp_list(
    query_aa: &[u8],
    subj_aa: &[u8],
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
    hits: &mut Vec<crate::protein_lookup::ProteinHit>,
    gap_open: i32,
    gap_extend: i32,
    x_drop_final: i32,
    scratch: &mut ProteinScratch,
) {
    // Re-run gapped DP with final x-drop, matching the ordinary protein
    // `Blast_TracebackFromHSPList` flow: process preliminary HSPs in score
    // order, skip any already contained by a traceback-accepted HSP, then add
    // tracebacked bounds to a fresh interval tree.
    hits.sort_by(compare_protein_hits_by_score);
    scratch
        .tb_tree
        .reset_for_query(query_aa.len() as i32 + 1, subj_aa.len() as i32 + 1);
    let tb_tree = &mut scratch.tb_tree;
    let mut keep = vec![true; hits.len()];
    for (idx, ph) in hits.iter_mut().enumerate() {
        if blast_interval_tree_contains_protein_hit_with_metadata(tb_tree, ph, 0, 0) {
            keep[idx] = false;
            continue;
        }
        let seed_q = ph.gapped_start_q;
        let seed_s = ph.gapped_start_s;
        let Some(gr) = blast_gapped_alignment_with_traceback(
            query_aa,
            subj_aa,
            seed_q,
            seed_s,
            matrix,
            gap_open,
            gap_extend,
            x_drop_final,
        ) else {
            keep[idx] = false;
            continue;
        };
        let q_slice = &query_aa[gr.query_start..gr.query_end];
        let s_slice = &subj_aa[gr.subject_start..gr.subject_end];
        let (qseq, sseq) =
            gr.edit_script
                .render_alignment(q_slice, s_slice, ncbistdaa_to_aminoacid_char);
        blast_hsp_update_with_traceback_protein_hit(ph, gr, qseq, sseq);
        blast_interval_tree_add_protein_hit_with_metadata(tb_tree, ph, 0, 0);
    }
    let mut idx = 0usize;
    hits.retain(|_| {
        let k = keep[idx];
        idx += 1;
        k
    });
    // Mirror the common-endpoint purge after traceback.
    blast_hsp_list_purge_hsps_with_common_endpoints_protein_hits(hits);
    hits.sort_by(compare_protein_hits_by_score);
}

fn compare_blast_hsps_by_score(
    a: &crate::hspstream::BlastHSP,
    b: &crate::hspstream::BlastHSP,
) -> std::cmp::Ordering {
    b.score
        .cmp(&a.score)
        .then_with(|| a.subject.offset.cmp(&b.subject.offset))
        .then_with(|| b.subject.end.cmp(&a.subject.end))
        .then_with(|| a.query.offset.cmp(&b.query.offset))
        .then_with(|| b.query.end.cmp(&a.query.end))
}

#[allow(dead_code)]
fn blast_seg_get_translated_offsets_api(
    offset: i32,
    end: i32,
    frame: i32,
    seq_length: i32,
) -> (i32, i32) {
    if frame > 0 {
        (
            offset * crate::util::CODON_LENGTH as i32 + frame,
            end * crate::util::CODON_LENGTH as i32 + frame - 1,
        )
    } else if frame < 0 {
        (
            seq_length - offset * crate::util::CODON_LENGTH as i32 + frame + 1,
            seq_length - end * crate::util::CODON_LENGTH as i32 + frame + 2,
        )
    } else {
        (offset + 1, end)
    }
}

/// NCBI: Blast_HSPGetAdjustedOffsets (blast_hits.c:1109), specialized to
/// internal `BlastHSP` records used by this API traceback path.
#[allow(dead_code)]
fn blast_hsp_get_adjusted_offsets(
    program: crate::program::ProgramType,
    hsp: &crate::hspstream::BlastHSP,
    query_length: i32,
    subject_length: i32,
) -> (i32, i32, i32, i32) {
    if hsp.gap_info.is_none() {
        return (
            hsp.query.offset + 1,
            hsp.query.end,
            hsp.subject.offset + 1,
            hsp.subject.end,
        );
    }

    if !crate::program::blast_query_is_translated(program)
        && !crate::program::blast_subject_is_translated(program)
    {
        if hsp.query.frame != hsp.subject.frame {
            let q_end = query_length - hsp.query.offset;
            let q_start = q_end - hsp.query.end + hsp.query.offset + 1;
            return (q_start, q_end, hsp.subject.end, hsp.subject.offset + 1);
        }
        return (
            hsp.query.offset + 1,
            hsp.query.end,
            hsp.subject.offset + 1,
            hsp.subject.end,
        );
    }

    let (q_start, q_end) = if crate::program::blast_query_is_translated(program) {
        blast_seg_get_translated_offsets_api(
            hsp.query.offset,
            hsp.query.end,
            hsp.query.frame as i32,
            query_length,
        )
    } else {
        (hsp.query.offset + 1, hsp.query.end)
    };
    let (s_start, s_end) = if crate::program::blast_subject_is_translated(program) {
        blast_seg_get_translated_offsets_api(
            hsp.subject.offset,
            hsp.subject.end,
            hsp.subject.frame as i32,
            subject_length,
        )
    } else {
        (hsp.subject.offset + 1, hsp.subject.end)
    };

    (q_start, q_end, s_start, s_end)
}

/// NCBI: s_HSPTest (blast_hits.c:993), specialized to `BlastHSP`.
fn s_hsp_test_blast_hsp(
    hsp: &crate::hspstream::BlastHSP,
    hit_options: &crate::options::HitSavingOptions,
    align_length: i32,
) -> bool {
    (hsp.num_ident as f64 * 100.0 < align_length as f64 * hit_options.percent_identity)
        || align_length < hit_options.min_hit_length
}

/// NCBI: s_Blast_HSPGetNumIdentitiesAndPositives / s_ComputeNumIdentities.
fn s_compute_num_identities_blast_hsp(
    hsp: &mut crate::hspstream::BlastHSP,
    query: &[u8],
    subject: &[u8],
    matrix: Option<&[[i32; AA_SIZE]; AA_SIZE]>,
) {
    let q_off = hsp.query.offset.max(0) as usize;
    let s_off = hsp.subject.offset.max(0) as usize;
    let q_length = (hsp.query.end - hsp.query.offset).max(0) as usize;
    let s_length = (hsp.subject.end - hsp.subject.offset).max(0) as usize;
    let mut q = q_off;
    let mut s = s_off;
    let mut num_ident = 0i32;
    let mut num_pos = 0i32;
    let mut align_length = 0i32;

    if let Some(script) = hsp.gap_info.as_ref() {
        for (op, count) in script.iter() {
            align_length += count;
            match op {
                crate::gapinfo::GapAlignOpType::Sub => {
                    for _ in 0..count {
                        if q < query.len() && s < subject.len() {
                            if query[q] == subject[s] {
                                num_ident += 1;
                            } else if let Some(matrix) = matrix {
                                let qi = query[q] as usize;
                                let si = subject[s] as usize;
                                if qi < AA_SIZE && si < AA_SIZE && matrix[qi][si] > 0 {
                                    num_pos += 1;
                                }
                            }
                        }
                        q += 1;
                        s += 1;
                    }
                }
                crate::gapinfo::GapAlignOpType::Del
                | crate::gapinfo::GapAlignOpType::Del1
                | crate::gapinfo::GapAlignOpType::Del2 => {
                    s += count.max(0) as usize;
                }
                crate::gapinfo::GapAlignOpType::Ins
                | crate::gapinfo::GapAlignOpType::Ins1
                | crate::gapinfo::GapAlignOpType::Ins2 => {
                    q += count.max(0) as usize;
                }
                crate::gapinfo::GapAlignOpType::Decline => {}
            }
        }
    } else if q_length == s_length {
        align_length = q_length as i32;
        for _ in 0..q_length {
            if q < query.len() && s < subject.len() {
                if query[q] == subject[s] {
                    num_ident += 1;
                } else if let Some(matrix) = matrix {
                    let qi = query[q] as usize;
                    let si = subject[s] as usize;
                    if qi < AA_SIZE && si < AA_SIZE && matrix[qi][si] > 0 {
                        num_pos += 1;
                    }
                }
            }
            q += 1;
            s += 1;
        }
    }

    hsp.num_ident = num_ident;
    hsp.num_positives = num_ident + num_pos;
    let _ = align_length;
}

/// NCBI: s_UpdateReevaluatedHSPUngapped (blast_hits.c:664), specialized to
/// the local preliminary protein-hit record.
fn s_update_reevaluated_hsp_ungapped(
    hsp: &mut crate::protein_lookup::ProteinHit,
    cutoff_score: i32,
    score: i32,
    best_q_start: usize,
    best_q_end: usize,
    best_s_start: usize,
    best_s_end: usize,
) -> bool {
    hsp.score = score;
    if hsp.score < cutoff_score {
        return true;
    }
    hsp.query_start = best_q_start;
    hsp.query_end = best_q_end;
    hsp.subject_start = best_s_start;
    hsp.subject_end = best_s_end;
    hsp.gapped_start_q = best_q_start;
    hsp.gapped_start_s = best_s_start;
    hsp.align_length = best_q_end.saturating_sub(best_q_start) as i32;
    false
}

fn blast_hsp_list_new(query_index: i32) -> crate::hspstream::BlastHSPList {
    crate::hspstream::BlastHSPList {
        oid: -1,
        query_index,
        hsp_array: Vec::new(),
        hspcnt: 0,
        allocated: 0,
        hsp_max: i32::MAX,
        do_not_reallocate: false,
        best_evalue: i32::MAX as f64,
    }
}

fn blast_hsp_list_push(
    hsp_list: &mut crate::hspstream::BlastHSPList,
    hsp: crate::hspstream::BlastHSP,
) {
    hsp_list.hsp_array.push(Some(hsp));
    hsp_list.hspcnt = hsp_list.hsp_array.len() as i32;
    hsp_list.allocated = hsp_list.hspcnt;
}

fn translated_frame_to_context(frame: i32) -> Option<usize> {
    match frame {
        1 => Some(0),
        2 => Some(1),
        3 => Some(2),
        -1 => Some(3),
        -2 => Some(4),
        -3 => Some(5),
        _ => None,
    }
}

fn protein_hit_to_blast_hsp(
    hit: crate::protein_lookup::ProteinHit,
    context: i32,
    query_frame: i32,
    subject_frame: i32,
) -> crate::hspstream::BlastHSP {
    crate::hspstream::BlastHSP {
        score: hit.score,
        num_ident: hit.num_ident,
        bit_score: 0.0,
        evalue: i32::MAX as f64,
        query: crate::hspstream::BlastSeg {
            frame: query_frame as i16,
            offset: hit.query_start as i32,
            end: hit.query_end as i32,
            gapped_start: hit.gapped_start_q as i32,
        },
        subject: crate::hspstream::BlastSeg {
            frame: subject_frame as i16,
            offset: hit.subject_start as i32,
            end: hit.subject_end as i32,
            gapped_start: hit.gapped_start_s as i32,
        },
        context,
        gap_info: None,
        num: 0,
        comp_adjustment_method: 0,
        pat_info: None,
        num_positives: 0,
        map_info: None,
    }
}

fn blast_hsp_to_protein_hit(hsp: &crate::hspstream::BlastHSP) -> crate::protein_lookup::ProteinHit {
    crate::protein_lookup::ProteinHit {
        query_start: hsp.query.offset.max(0) as usize,
        query_end: hsp.query.end.max(hsp.query.offset).max(0) as usize,
        subject_start: hsp.subject.offset.max(0) as usize,
        subject_end: hsp.subject.end.max(hsp.subject.offset).max(0) as usize,
        score: hsp.score,
        num_ident: hsp.num_ident,
        align_length: (hsp.query.end - hsp.query.offset).max(0),
        mismatches: 0,
        gap_opens: 0,
        qseq: None,
        sseq: None,
        scaled_score: None,
        adjusted_evalue: None,
        gapped_start_q: hsp.query.gapped_start.max(0) as usize,
        gapped_start_s: hsp.subject.gapped_start.max(0) as usize,
    }
}

fn blast_hsp_update_from_protein_hit(
    hsp: &mut crate::hspstream::BlastHSP,
    hit: crate::protein_lookup::ProteinHit,
) {
    hsp.score = hit.score;
    hsp.num_ident = hit.num_ident;
    hsp.query.offset = hit.query_start as i32;
    hsp.query.end = hit.query_end as i32;
    hsp.query.gapped_start = hit.gapped_start_q as i32;
    hsp.subject.offset = hit.subject_start as i32;
    hsp.subject.end = hit.subject_end as i32;
    hsp.subject.gapped_start = hit.gapped_start_s as i32;
}

fn tblastn_link_internal_hsps_to_api_hsps(
    hsp_list: &mut crate::hspstream::BlastHSPList,
    query_info: &crate::queryinfo::QueryInfo,
    query_context: usize,
    subject_stat_length: i32,
    prot_kbp: &KarlinBlk,
    gumbel_blk: Option<&GumbelBlk>,
    link_gbp_db_length: i64,
    word_cutoff_score_min: i32,
    max_intron_length: i32,
    max_evalue: f64,
) -> Vec<Hsp> {
    hsp_list.query_index = query_context as i32;
    let original_hsps: Vec<crate::hspstream::BlastHSP> =
        hsp_list.hsp_array.iter().filter_map(Clone::clone).collect();
    let score_block = crate::link_hsps::LinkScoreBlock {
        kbp: vec![prot_kbp.clone()],
        kbp_gap: vec![prot_kbp.clone()],
        gbp: gumbel_blk.cloned(),
        link_gbp_db_length: Some(link_gbp_db_length.max(1)),
        recompute_evalues_before_uneven_linking: true,
    };
    let link_params = crate::link_hsps::LinkHSPParameters {
        gap_prob: crate::stat::BLAST_GAP_PROB_GAPPED,
        gap_decay_rate: crate::stat::BLAST_GAP_DECAY_RATE_GAPPED,
        longest_intron: translated_gapped_link_longest_intron(max_intron_length),
        cutoff_small_gap: word_cutoff_score_min,
        ..crate::link_hsps::LinkHSPParameters::default()
    };
    crate::link_hsps::blast_link_hsp_list(
        crate::program::TBLASTN,
        hsp_list,
        query_info,
        subject_stat_length,
        &score_block,
        &link_params,
        true,
    );
    hsp_list
        .hsp_array
        .retain(|hsp| hsp.as_ref().is_some_and(|hsp| hsp.evalue <= max_evalue));
    hsp_list.hspcnt = hsp_list.hsp_array.len() as i32;

    hsp_list
        .hsp_array
        .iter()
        .filter_map(|hsp| hsp.as_ref())
        .filter_map(|linked| {
            let original = original_hsps
                .iter()
                .find(|original| same_tblastn_link_head(original, &linked))?;
            let mut hsp = tblastn_blast_hsp_to_api_hsp(original, prot_kbp);
            hsp.evalue = linked.evalue;
            hsp.num_links = linked.num.max(1);
            Some(hsp)
        })
        .collect()
}

fn same_tblastn_link_head(
    original: &crate::hspstream::BlastHSP,
    linked: &crate::hspstream::BlastHSP,
) -> bool {
    original.score == linked.score
        && original.context == linked.context
        && original.query.frame == linked.query.frame
        && original.query.offset == linked.query.offset
        && original.query.gapped_start == linked.query.gapped_start
        && original.subject.frame == linked.subject.frame
        && original.subject.offset == linked.subject.offset
        && original.subject.gapped_start == linked.subject.gapped_start
}

fn tblastn_blast_hsp_to_api_hsp(hsp: &crate::hspstream::BlastHSP, prot_kbp: &KarlinBlk) -> Hsp {
    let (subject_start, subject_end) = crate::util::protein_to_oriented_nuc_coords(
        hsp.subject.offset.max(0) as usize,
        hsp.subject.end.max(hsp.subject.offset).max(0) as usize,
        hsp.subject.frame as i32,
    );
    let (subject_gapped_start, _) = crate::util::protein_to_oriented_nuc_coords(
        hsp.subject.gapped_start.max(0) as usize,
        hsp.subject.gapped_start.max(0) as usize + 1,
        hsp.subject.frame as i32,
    );
    let query_start = hsp.query.offset.max(0) as usize;
    let query_end = hsp.query.end.max(hsp.query.offset).max(0) as usize;
    let alignment_length = hsp
        .gap_info
        .as_ref()
        .map(|script| script.alignment_length().max(0) as usize)
        .unwrap_or(query_end.saturating_sub(query_start));
    Hsp {
        score: hsp.score,
        bit_score: prot_kbp.raw_to_bit(hsp.score),
        evalue: hsp.evalue,
        query_start,
        query_end,
        subject_start,
        subject_end,
        query_gapped_start: hsp.query.gapped_start.max(0) as usize,
        subject_gapped_start,
        query_link_start: query_start,
        query_link_end: query_end,
        query_link_gapped_start: hsp.query.gapped_start.max(0) as usize,
        subject_link_start: hsp.subject.offset.max(0) as usize,
        subject_link_end: hsp.subject.end.max(hsp.subject.offset).max(0) as usize,
        subject_link_gapped_start: hsp.subject.gapped_start.max(0) as usize,
        link_score: None,
        link_lambda: None,
        num_identities: hsp.num_ident.max(0) as usize,
        num_gaps: hsp
            .gap_info
            .as_ref()
            .map(crate::blast_kappa::gap_edit_script_num_gap_opens)
            .unwrap_or(0)
            .max(0) as usize,
        alignment_length,
        query_aln: Vec::new(),
        midline: Vec::new(),
        subject_aln: Vec::new(),
        query_frame: 0,
        subject_frame: hsp.subject.frame as i32,
        num_links: hsp.num.max(1),
        comp_adjust_method: hsp.comp_adjustment_method as u8,
    }
}

fn sort_blast_hsp_list_by_score(hsp_list: &mut crate::hspstream::BlastHSPList) {
    hsp_list
        .hsp_array
        .sort_by(|a, b| match (a.as_ref(), b.as_ref()) {
            (Some(a), Some(b)) => compare_blast_hsps_by_score(a, b),
            (Some(_), None) => std::cmp::Ordering::Less,
            (None, Some(_)) => std::cmp::Ordering::Greater,
            (None, None) => std::cmp::Ordering::Equal,
        });
    hsp_list.hspcnt = hsp_list
        .hsp_array
        .iter()
        .filter(|hsp| hsp.is_some())
        .count() as i32;
}

/// NCBI: `Blast_HSPListPurgeHSPsWithCommonEndpoints`, specialized to BLASTX
/// internal HSPs. The C comparators sort by context/offsets/score/range; the
/// duplicate test also includes subject frame.
fn purge_blastx_hsp_list_with_common_endpoints(hsp_list: &mut crate::hspstream::BlastHSPList) {
    let hsps = &mut hsp_list.hsp_array;
    if hsps.len() <= 1 {
        return;
    }
    hsps.sort_by(|a, b| {
        let Some(a) = a.as_ref() else {
            return std::cmp::Ordering::Greater;
        };
        let Some(b) = b.as_ref() else {
            return std::cmp::Ordering::Less;
        };
        a.context
            .cmp(&b.context)
            .then(a.query.offset.cmp(&b.query.offset))
            .then(a.subject.offset.cmp(&b.subject.offset))
            .then(b.score.cmp(&a.score))
            .then(b.query.end.cmp(&a.query.end))
            .then(b.subject.end.cmp(&a.subject.end))
    });
    let mut keep = vec![true; hsps.len()];
    let mut i = 0usize;
    while i < hsps.len() {
        if !keep[i] {
            i += 1;
            continue;
        }
        let mut j = i + 1;
        while j < hsps.len()
            && hsps[i]
                .as_ref()
                .zip(hsps[j].as_ref())
                .is_some_and(|(a, b)| {
                    a.context == b.context
                        && a.query.offset == b.query.offset
                        && a.subject.offset == b.subject.offset
                        && a.subject.frame == b.subject.frame
                })
        {
            keep[j] = false;
            j += 1;
        }
        i = j;
    }
    let mut idx = 0usize;
    hsps.retain(|_| {
        let k = keep[idx];
        idx += 1;
        k
    });

    if hsps.len() <= 1 {
        return;
    }
    hsps.sort_by(|a, b| {
        let Some(a) = a.as_ref() else {
            return std::cmp::Ordering::Greater;
        };
        let Some(b) = b.as_ref() else {
            return std::cmp::Ordering::Less;
        };
        a.context
            .cmp(&b.context)
            .then(a.query.end.cmp(&b.query.end))
            .then(a.subject.end.cmp(&b.subject.end))
            .then(b.score.cmp(&a.score))
            .then(b.query.offset.cmp(&a.query.offset))
            .then(b.subject.offset.cmp(&a.subject.offset))
    });
    let mut keep = vec![true; hsps.len()];
    let mut i = 0usize;
    while i < hsps.len() {
        if !keep[i] {
            i += 1;
            continue;
        }
        let mut j = i + 1;
        while j < hsps.len()
            && hsps[i]
                .as_ref()
                .zip(hsps[j].as_ref())
                .is_some_and(|(a, b)| {
                    a.context == b.context
                        && a.query.end == b.query.end
                        && a.subject.end == b.subject.end
                        && a.subject.frame == b.subject.frame
                })
        {
            keep[j] = false;
            j += 1;
        }
        i = j;
    }
    let mut idx = 0usize;
    hsps.retain(|_| {
        let k = keep[idx];
        idx += 1;
        k
    });
    hsp_list.hspcnt = hsp_list.hsp_array.len() as i32;
    hsp_list.allocated = hsp_list.hspcnt;
}

/// NCBI: `Blast_HSPListPurgeHSPsWithCommonEndpoints`, specialized to TBLASTN
/// internal HSPs. The C comparators sort by context/offsets/score/range; the
/// duplicate test also includes subject frame.
fn purge_tblastn_hsp_list_with_common_endpoints(hsp_list: &mut crate::hspstream::BlastHSPList) {
    let hsps = &mut hsp_list.hsp_array;
    if hsps.len() <= 1 {
        return;
    }
    hsps.sort_by(|a, b| {
        let Some(a) = a.as_ref() else {
            return std::cmp::Ordering::Greater;
        };
        let Some(b) = b.as_ref() else {
            return std::cmp::Ordering::Less;
        };
        a.context
            .cmp(&b.context)
            .then(a.query.offset.cmp(&b.query.offset))
            .then(a.subject.offset.cmp(&b.subject.offset))
            .then(b.score.cmp(&a.score))
            .then(b.query.end.cmp(&a.query.end))
            .then(b.subject.end.cmp(&a.subject.end))
    });
    let mut keep = vec![true; hsps.len()];
    let mut i = 0usize;
    while i < hsps.len() {
        if !keep[i] {
            i += 1;
            continue;
        }
        let mut j = i + 1;
        while j < hsps.len()
            && hsps[i]
                .as_ref()
                .zip(hsps[j].as_ref())
                .is_some_and(|(a, b)| {
                    a.context == b.context
                        && a.query.offset == b.query.offset
                        && a.subject.offset == b.subject.offset
                        && a.subject.frame == b.subject.frame
                })
        {
            keep[j] = false;
            j += 1;
        }
        i = j;
    }
    let mut idx = 0usize;
    hsps.retain(|_| {
        let k = keep[idx];
        idx += 1;
        k
    });

    if hsps.len() <= 1 {
        return;
    }
    hsps.sort_by(|a, b| {
        let Some(a) = a.as_ref() else {
            return std::cmp::Ordering::Greater;
        };
        let Some(b) = b.as_ref() else {
            return std::cmp::Ordering::Less;
        };
        a.context
            .cmp(&b.context)
            .then(a.query.end.cmp(&b.query.end))
            .then(a.subject.end.cmp(&b.subject.end))
            .then(b.score.cmp(&a.score))
            .then(b.query.offset.cmp(&a.query.offset))
            .then(b.subject.offset.cmp(&a.subject.offset))
    });
    let mut keep = vec![true; hsps.len()];
    let mut i = 0usize;
    while i < hsps.len() {
        if !keep[i] {
            i += 1;
            continue;
        }
        let mut j = i + 1;
        while j < hsps.len()
            && hsps[i]
                .as_ref()
                .zip(hsps[j].as_ref())
                .is_some_and(|(a, b)| {
                    a.context == b.context
                        && a.query.end == b.query.end
                        && a.subject.end == b.subject.end
                        && a.subject.frame == b.subject.frame
                })
        {
            keep[j] = false;
            j += 1;
        }
        i = j;
    }
    let mut idx = 0usize;
    hsps.retain(|_| {
        let k = keep[idx];
        idx += 1;
        k
    });
    hsp_list.hspcnt = hsp_list.hsp_array.len() as i32;
    hsp_list.allocated = hsp_list.hspcnt;
}

/// NCBI: `Blast_TracebackFromHSPList` over one BLASTX subject HSP list.
/// Query frame/context stays attached until after traceback containment and
/// common-endpoint purge; public HSP conversion is a final formatting step.
#[allow(clippy::too_many_arguments)]
fn blastx_traceback_subject_hsp_list(
    translation_buffer: &[u8],
    frame_offsets: &[u32],
    subj_aa: &[u8],
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
    hsp_list: &mut crate::hspstream::BlastHSPList,
    gap_open: i32,
    gap_extend: i32,
    x_drop_final: i32,
    scratch: &mut ProteinScratch,
    max_query_prot_len: usize,
) {
    sort_blast_hsp_list_by_score(hsp_list);
    scratch
        .tb_tree
        .reset_for_query(max_query_prot_len as i32 + 1, subj_aa.len() as i32 + 1);
    let tb_tree = &mut scratch.tb_tree;
    let mut keep = vec![true; hsp_list.hsp_array.len()];
    for (idx, hsp_slot) in hsp_list.hsp_array.iter_mut().enumerate() {
        let Some(hsp) = hsp_slot.as_mut() else {
            keep[idx] = false;
            continue;
        };
        let item_context = hsp.context.max(0) as usize;
        if item_context + 1 >= frame_offsets.len() {
            keep[idx] = false;
            continue;
        }
        let item_query_range = (frame_offsets[item_context] + 1).max(0) as usize
            ..frame_offsets[item_context + 1].max(0) as usize;
        let mut hit = blast_hsp_to_protein_hit(hsp);
        if blast_interval_tree_contains_protein_hit_with_metadata(
            tb_tree,
            &hit,
            item_context as i32,
            0,
        ) {
            keep[idx] = false;
            continue;
        }
        let query_prot = &translation_buffer[item_query_range];
        let seed_q = hit.gapped_start_q;
        let seed_s = hit.gapped_start_s;
        let Some(gr) = blast_gapped_alignment_with_traceback(
            query_prot,
            subj_aa,
            seed_q,
            seed_s,
            matrix,
            gap_open,
            gap_extend,
            x_drop_final,
        ) else {
            keep[idx] = false;
            continue;
        };
        let q_slice = &query_prot[gr.query_start..gr.query_end];
        let s_slice = &subj_aa[gr.subject_start..gr.subject_end];
        let (qseq, sseq) =
            gr.edit_script
                .render_alignment(q_slice, s_slice, ncbistdaa_to_aminoacid_char);
        blast_hsp_update_with_traceback_protein_hit(&mut hit, gr, qseq, sseq);
        blast_interval_tree_add_protein_hit_with_metadata(tb_tree, &hit, item_context as i32, 0);
        blast_hsp_update_from_protein_hit(hsp, hit.clone());
    }
    let mut idx = 0usize;
    hsp_list.hsp_array.retain(|_| {
        let k = keep[idx];
        idx += 1;
        k
    });

    purge_blastx_hsp_list_with_common_endpoints(hsp_list);
    sort_blast_hsp_list_by_score(hsp_list);
}

/// NCBI: Blast_TracebackFromHSPList over one translated subject HSP list.
/// Each element remembers the translated subject frame and protein slice, but
/// containment is owned by one subject-level traceback tree.
#[allow(clippy::too_many_arguments)]
fn tblastn_traceback_subject_hsp_list(
    query_aa: &[u8],
    translation_buffer: &[u8],
    frame_offsets: &[u32],
    subject_source_ncbi4na: &[u8],
    gen_code_string: &[u8; 64],
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
    hsp_list: &mut crate::hspstream::BlastHSPList,
    gap_open: i32,
    gap_extend: i32,
    x_drop_final: i32,
    scratch: &mut ProteinScratch,
    query_index: i32,
    _max_subj_prot_len: usize,
) -> i32 {
    hsp_list.query_index = query_index;
    let mut stat_length = subject_source_ncbi4na.len() as i32;
    let mut sbp = crate::stat::blast_score_blk_new(crate::encoding::BLASTAA_SEQ_CODE, 1)
        .expect("protein score block");
    for row in 0..AA_SIZE {
        sbp.matrix.data[row][..AA_SIZE].copy_from_slice(&matrix[row]);
    }
    let scoring_options = crate::options::ScoringOptions {
        matrix_path: None,
        reward: 0,
        penalty: 0,
        gap_open,
        gap_extend,
        shift_pen: i16::MAX as i32,
        gapped_calculation: true,
        complexity_adjusted_scoring: false,
        matrix_name: None,
        is_ooframe: false,
        program_number: crate::program::UNDEFINED,
    };
    let score_params = crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
    let ext_options = crate::options::ExtensionOptions::new_blastp();
    let hit_params = crate::parameters::HitSavingParameters {
        options: crate::options::HitSavingOptions::default(),
        cutoff_score_min: 0,
        low_score: Vec::new(),
        cutoffs: Vec::new(),
        link_hsp_params: None,
        restricted_align: false,
        do_sum_stats: false,
        mask_level: 0,
        prelim_evalue: 0.0,
    };
    let query_blk = crate::util::BlastSequenceBlk {
        sequence: Some(query_aa.to_vec()),
        sequence_start: Some(query_aa.to_vec()),
        length: query_aa.len() as i32,
        ..Default::default()
    };
    let mut subject_blk = crate::util::BlastSequenceBlk {
        sequence: Some(subject_source_ncbi4na.to_vec()),
        sequence_start: Some(subject_source_ncbi4na.to_vec()),
        length: subject_source_ncbi4na.len() as i32,
        gen_code_string: Some(*gen_code_string),
        ..Default::default()
    };
    let mut gap_align = crate::blast_kappa::BlastGapAlignWorkspace {
        gap_x_dropoff: x_drop_final,
        ..Default::default()
    };
    if blast_traceback_from_hsp_list(
        crate::program::TBLASTN,
        hsp_list,
        &query_blk,
        &mut subject_blk,
        &crate::queryinfo::QueryInfo::new_blastp(&[query_aa.len()]),
        &mut gap_align,
        &sbp,
        &score_params,
        &ext_options,
        &hit_params,
        Some(gen_code_string),
        None,
        scratch,
        Some(&mut stat_length),
    ) == 0
    {
        for hsp_slot in hsp_list.hsp_array.iter_mut() {
            let Some(hsp) = hsp_slot.as_mut() else {
                continue;
            };
            let Some(ctx) = translated_frame_to_context(hsp.subject.frame as i32) else {
                continue;
            };
            if ctx + 1 >= frame_offsets.len() {
                continue;
            }
            let prot_range =
                (frame_offsets[ctx] + 1).max(0) as usize..frame_offsets[ctx + 1].max(0) as usize;
            let mut hit = blast_hsp_to_protein_hit(&hsp);
            if let Some(edit_script) = hsp.gap_info.as_ref() {
                let q_start = hsp.query.offset.max(0) as usize;
                let q_end = hsp.query.end.max(hsp.query.offset).max(0) as usize;
                let prot = &translation_buffer[prot_range.clone()];
                let range_start = prot_range.start as i32;
                let mut s_start_i = hsp.subject.offset;
                let mut s_end_i = hsp.subject.end.max(hsp.subject.offset);
                if s_end_i as usize > prot.len() && s_start_i >= range_start {
                    s_start_i -= range_start;
                    s_end_i -= range_start;
                }
                let s_start = s_start_i.max(0) as usize;
                let s_end = s_end_i.max(s_start_i).max(0) as usize;
                if q_end <= query_aa.len() && s_end <= prot.len() {
                    let q_slice = &query_aa[q_start..q_end];
                    let s_slice = &prot[s_start..s_end];
                    let (qseq, sseq) =
                        edit_script.render_alignment(q_slice, s_slice, ncbistdaa_to_aminoacid_char);
                    let (align_length, num_ident, gap_opens) =
                        edit_script.count_identities(q_slice, s_slice);
                    hit.qseq = Some(qseq);
                    hit.sseq = Some(sseq);
                    hit.align_length = align_length;
                    hit.num_ident = num_ident;
                    hit.gap_opens = gap_opens;
                    hit.mismatches = (align_length - num_ident - gap_opens).max(0);
                }
            }
            hsp.num = 0;
        }
        purge_tblastn_hsp_list_with_common_endpoints(hsp_list);
        sort_blast_hsp_list_by_score(hsp_list);
        return stat_length;
    }

    stat_length
}

/// Ungapped-only variant: skip `blast_get_gapped_score` and return each
/// surviving ungapped HSP for callers that implement NCBI's `-ungapped`
/// mode (blastx/tblastn). Mirrors `protein_alignment_hits` minus the
/// gapped phase.
/// blast-rs: Protein search adapter for public ungapped mode; not a direct NCBI C port.
fn protein_alignment_hits_ungapped_only(
    query_aa: &[u8],
    subj_aa: &[u8],
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
    lookup_table: &crate::protein_lookup::ProteinLookupTable,
    x_drop_ungapped: i32,
    window: i32,
    seed_cutoff: i32,
    scratch: &mut ProteinScratch,
) -> Vec<crate::protein_lookup::ProteinHit> {
    // NCBI gates the ungapped init-HSP save on `score >= cutoffs->cutoff_score`
    // (`aa_ungapped.c:588`). In `-ungapped` mode that cutoff is `seed_cutoff`.
    let mut ungapped_hits = crate::protein_lookup::protein_scan_with_extend_word(
        query_aa,
        subj_aa,
        matrix,
        lookup_table,
        x_drop_ungapped,
        window,
        seed_cutoff.max(1),
        &mut scratch.diag_table,
    );
    scratch.diag_table.exit(subj_aa.len());
    ungapped_hits
        .iter_mut()
        .filter(|uh| uh.score >= seed_cutoff)
        .for_each(|uh| {
            let qseq = query_aa[uh.query_start..uh.query_end]
                .iter()
                .map(|&aa| ncbistdaa_to_aminoacid_char(aa))
                .collect();
            let sseq = subj_aa[uh.subject_start..uh.subject_end]
                .iter()
                .map(|&aa| ncbistdaa_to_aminoacid_char(aa))
                .collect();
            uh.qseq = Some(qseq);
            uh.sseq = Some(sseq);
        });
    ungapped_hits
        .into_iter()
        .filter(|uh| uh.score >= seed_cutoff)
        .collect()
}

// ── Search parameters ───────────────────────────────────────────────────────

/// Scoring matrix type for protein searches.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MatrixType {
    Blosum45,
    Blosum50,
    Blosum62,
    Blosum80,
    Blosum90,
    Pam30,
    Pam70,
    Pam250,
    Identity,
}

/// Scoring matrix wrapper (28x28 NCBIstdaa-indexed scores with metadata).
#[derive(Debug, Clone)]
pub struct ScoringMatrix {
    pub scores: [[i32; AA_SIZE]; AA_SIZE],
    pub min_score: i32,
    pub name: String,
}

impl ScoringMatrix {
    /// blast-rs: Public scoring-matrix constructor; not a direct NCBI C port.
    pub fn from_type(mt: MatrixType) -> Self {
        let scores = *get_matrix(mt);
        let min_score = scores
            .iter()
            .flat_map(|row| row.iter())
            .filter(|&&v| v > i32::MIN / 4)
            .copied()
            .min()
            .unwrap_or(-4);
        let name = format!("{:?}", mt).to_uppercase();
        ScoringMatrix {
            scores,
            min_score,
            name,
        }
    }
    /// blast-rs: Public BLOSUM62 constructor; not a direct NCBI C port.
    pub fn blosum62() -> Self {
        Self::from_type(MatrixType::Blosum62)
    }
    /// blast-rs: Public BLOSUM45 constructor; not a direct NCBI C port.
    pub fn blosum45() -> Self {
        Self::from_type(MatrixType::Blosum45)
    }
    /// blast-rs: Public BLOSUM50 constructor; not a direct NCBI C port.
    pub fn blosum50() -> Self {
        Self::from_type(MatrixType::Blosum50)
    }
    /// blast-rs: Public BLOSUM80 constructor; not a direct NCBI C port.
    pub fn blosum80() -> Self {
        Self::from_type(MatrixType::Blosum80)
    }
    /// blast-rs: Public BLOSUM90 constructor; not a direct NCBI C port.
    pub fn blosum90() -> Self {
        Self::from_type(MatrixType::Blosum90)
    }
    /// blast-rs: Public PAM30 constructor; not a direct NCBI C port.
    pub fn pam30() -> Self {
        Self::from_type(MatrixType::Pam30)
    }
    /// blast-rs: Public PAM70 constructor; not a direct NCBI C port.
    pub fn pam70() -> Self {
        Self::from_type(MatrixType::Pam70)
    }
    /// blast-rs: Public PAM250 constructor; not a direct NCBI C port.
    pub fn pam250() -> Self {
        Self::from_type(MatrixType::Pam250)
    }
    /// blast-rs: Public identity-matrix constructor; not a direct NCBI C port.
    pub fn identity() -> Self {
        Self::from_type(MatrixType::Identity)
    }
    /// blast-rs: Public scoring-matrix accessor; not a direct NCBI C port.
    pub fn score(&self, a: u8, b: u8) -> i32 {
        self.scores[a as usize & 0x1F][b as usize & 0x1F]
    }
}

/// Search configuration with builder pattern.
#[derive(Debug, Clone)]
pub struct SearchParams {
    pub word_size: usize,
    pub matrix: MatrixType,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub evalue_threshold: f64,
    pub max_target_seqs: usize,
    pub x_drop_ungapped: i32,
    pub x_drop_gapped: i32,
    pub ungapped_cutoff: i32,
    pub min_score: i32,
    pub match_score: i32,
    pub mismatch: i32,
    pub num_threads: usize,
    pub filter_low_complexity: bool,
    pub seg_window: usize,
    pub seg_locut: f64,
    pub seg_hicut: f64,
    pub word_threshold: Option<f64>,
    pub comp_adjust: u8,
    pub strand: String,
    pub query_gencode: u8,
    pub db_gencode: u8,
    pub db_length: i64,
    pub effective_search_space: i64,
    pub max_intron_length: i32,
    pub max_hsps: Option<usize>,
    pub culling_limit: Option<usize>,
    pub two_hit: bool,
    pub two_hit_window: usize,
    /// Nucleotide two-hit window size (NCBI `-window_size`). 0 (the default for
    /// blastn/megablast) disables the two-hit gate, leaving the scan one-hit.
    /// Only the nucleotide search path consumes this; protein paths use
    /// `two_hit_window`.
    pub window_size: i32,
    pub x_drop_final: i32,
    /// NCBI `BlastExtensionOptions::chaining`. Protein gapped extension only
    /// runs `s_ChainingAlignment` when this is enabled and the matrix is
    /// BLOSUM62.
    pub chaining: bool,
    pub soft_masking: bool,
    pub lcase_masking: bool,
    /// Enable NCBI linked sum-statistics e-value adjustment where the program
    /// supports it. BLASTN enables this by default, matching CLI BLAST+.
    pub sum_stats: bool,
    /// Run ungapped-only mode (NCBI `-ungapped`). Skips the gapped DP / redo
    /// alignment phase and emits each surviving ungapped HSP from the
    /// preliminary scan with `kbp_std` (ungapped Karlin) parameters.
    pub ungapped: bool,
    /// Optional reusable Rayon pool for API searches that parallelize over
    /// database subjects. When present, it is used instead of constructing a
    /// per-search pool from `num_threads`.
    pub thread_pool: Option<Arc<rayon::ThreadPool>>,
}

impl SearchParams {
    /// blast-rs: Public blastp SearchParams constructor; not a direct NCBI C port.
    pub fn blastp() -> Self {
        Self::blastp_defaults()
    }
    /// blast-rs: Public blastn SearchParams constructor; not a direct NCBI C port.
    pub fn blastn() -> Self {
        Self::blastn_defaults()
    }

    /// blast-rs: Public blastx SearchParams constructor; not a direct NCBI C port.
    pub fn blastx() -> Self {
        let mut params = Self::blastp_defaults();
        params.max_intron_length = 0;
        params
    }
    /// blast-rs: Public tblastn SearchParams constructor; not a direct NCBI C port.
    pub fn tblastn() -> Self {
        let mut params = Self::blastp_defaults();
        params.max_intron_length = 0;
        params
    }
    /// blast-rs: Public tblastx SearchParams constructor; not a direct NCBI C port.
    pub fn tblastx() -> Self {
        let mut params = Self::blastp_defaults();
        params.comp_adjust = 0;
        params.max_intron_length = 0;
        params.x_drop_gapped = crate::stat::BLAST_GAP_X_DROPOFF_TBLASTX;
        params.x_drop_final = crate::stat::BLAST_GAP_X_DROPOFF_FINAL_TBLASTX;
        params
    }

    /// blast-rs: Public blastp default-parameter bundle; not a direct NCBI C port.
    pub fn blastp_defaults() -> Self {
        SearchParams {
            word_size: crate::stat::BLAST_WORDSIZE_PROT as usize,
            matrix: MatrixType::Blosum62,
            gap_open: crate::stat::BLAST_GAP_OPEN_PROT,
            gap_extend: crate::stat::BLAST_GAP_EXTN_PROT,
            evalue_threshold: crate::stat::BLAST_EXPECT_VALUE,
            max_target_seqs: crate::stat::BLAST_HITLIST_SIZE,
            x_drop_ungapped: crate::stat::BLAST_UNGAPPED_X_DROPOFF_PROT,
            x_drop_gapped: crate::stat::BLAST_GAP_X_DROPOFF_PROT,
            ungapped_cutoff: 0,
            min_score: 0,
            match_score: 1,
            mismatch: -2,
            num_threads: 1,
            filter_low_complexity: true,
            seg_window: 12,
            seg_locut: 2.2,
            seg_hicut: 2.5,
            word_threshold: None,
            comp_adjust: 2, // NCBI default: conditional compositional matrix adjust
            strand: "both".to_string(),
            query_gencode: 1,
            db_gencode: 1,
            db_length: 0,
            effective_search_space: 0,
            max_intron_length: crate::stat::DEFAULT_LONGEST_INTRON as i32,
            max_hsps: None,
            culling_limit: None,
            two_hit: false,
            two_hit_window: crate::stat::BLAST_WINDOW_SIZE_PROT as usize,
            window_size: 0,
            x_drop_final: crate::stat::BLAST_GAP_X_DROPOFF_FINAL_PROT,
            chaining: true,
            soft_masking: false,
            lcase_masking: false,
            sum_stats: true,
            ungapped: false,
            thread_pool: None,
        }
    }

    /// blast-rs: Public blastn default-parameter bundle; not a direct NCBI C port.
    pub fn blastn_defaults() -> Self {
        SearchParams {
            word_size: crate::stat::BLAST_WORDSIZE_NUCL as usize,
            matrix: MatrixType::Blosum62,
            gap_open: crate::stat::BLAST_GAP_OPEN_NUCL,
            gap_extend: crate::stat::BLAST_GAP_EXTN_NUCL,
            evalue_threshold: crate::stat::BLAST_EXPECT_VALUE,
            max_target_seqs: crate::stat::BLAST_HITLIST_SIZE,
            x_drop_ungapped: crate::stat::BLAST_UNGAPPED_X_DROPOFF_NUCL,
            x_drop_gapped: crate::stat::BLAST_GAP_X_DROPOFF_NUCL,
            ungapped_cutoff: 0,
            min_score: 0,
            match_score: 2,
            mismatch: -3,
            num_threads: 1,
            filter_low_complexity: true,
            seg_window: 12,
            seg_locut: 2.2,
            seg_hicut: 2.5,
            word_threshold: None,
            comp_adjust: 0,
            strand: "both".to_string(),
            query_gencode: 1,
            db_gencode: 1,
            db_length: 0,
            effective_search_space: 0,
            max_intron_length: crate::stat::DEFAULT_LONGEST_INTRON as i32,
            max_hsps: None,
            culling_limit: None,
            two_hit: false,
            // NCBI `BLAST_WINDOW_SIZE_NUCL = 0` would disable two-hit
            // regardless of `two_hit`; Rust keeps the protein default
            // so enabling two-hit on the blastn path still has a usable
            // window (40). `two_hit=false` by default, so unused.
            two_hit_window: crate::stat::BLAST_WINDOW_SIZE_PROT as usize,
            // NCBI `BLAST_WINDOW_SIZE_NUCL = 0`: blastn/megablast default to a
            // one-hit scan (the two-hit gate in `search.rs` is inert at 0).
            window_size: 0,
            x_drop_final: crate::stat::BLAST_GAP_X_DROPOFF_FINAL_NUCL,
            chaining: true,
            soft_masking: true,
            lcase_masking: false,
            sum_stats: true,
            ungapped: false,
            thread_pool: None,
        }
    }

    // Builder methods
    /// blast-rs: SearchParams builder setter; not a direct NCBI C port.
    pub fn evalue(mut self, v: f64) -> Self {
        self.evalue_threshold = v;
        self
    }
    /// blast-rs: SearchParams builder setter; not a direct NCBI C port.
    pub fn max_target_seqs(mut self, v: usize) -> Self {
        self.max_target_seqs = v;
        self
    }
    /// blast-rs: SearchParams builder setter; not a direct NCBI C port.
    pub fn num_threads(mut self, v: usize) -> Self {
        self.num_threads = v;
        self
    }
    /// Use an existing Rayon thread pool for parallel API searches.
    ///
    /// Supplying a pool enables the parallel path even when `num_threads` is
    /// left at its default value; the pool's own size controls concurrency.
    /// blast-rs: SearchParams builder setter; not a direct NCBI C port.
    pub fn thread_pool(mut self, pool: Arc<rayon::ThreadPool>) -> Self {
        self.thread_pool = Some(pool);
        self
    }
    /// blast-rs: SearchParams builder setter; not a direct NCBI C port.
    pub fn filter_low_complexity(mut self, v: bool) -> Self {
        self.filter_low_complexity = v;
        self
    }
    /// blast-rs: SearchParams builder setter; not a direct NCBI C port.
    pub fn seg_options(mut self, window: usize, locut: f64, hicut: f64) -> Self {
        self.seg_window = window;
        self.seg_locut = locut;
        self.seg_hicut = hicut;
        self
    }
    /// blast-rs: SearchParams builder setter; not a direct NCBI C port.
    pub fn word_threshold(mut self, v: f64) -> Self {
        self.word_threshold = Some(v);
        self
    }
    /// blast-rs: SearchParams builder setter; not a direct NCBI C port.
    pub fn comp_adjust(mut self, v: u8) -> Self {
        self.comp_adjust = v;
        self
    }
    /// blast-rs: SearchParams builder setter; not a direct NCBI C port.
    pub fn strand(mut self, v: &str) -> Self {
        self.strand = v.to_string();
        self
    }
    /// blast-rs: SearchParams builder setter; not a direct NCBI C port.
    pub fn word_size(mut self, v: usize) -> Self {
        self.word_size = v;
        self
    }
    /// blast-rs: SearchParams builder setter; not a direct NCBI C port.
    pub fn matrix(mut self, v: MatrixType) -> Self {
        self.matrix = v;
        self
    }
    /// blast-rs: SearchParams builder setter; not a direct NCBI C port.
    pub fn gap_open(mut self, v: i32) -> Self {
        self.gap_open = v;
        self
    }
    /// blast-rs: SearchParams builder setter; not a direct NCBI C port.
    pub fn gap_extend(mut self, v: i32) -> Self {
        self.gap_extend = v;
        self
    }
    /// blast-rs: SearchParams builder setter; not a direct NCBI C port.
    pub fn match_score(mut self, v: i32) -> Self {
        self.match_score = v;
        self
    }
    /// blast-rs: SearchParams builder setter; not a direct NCBI C port.
    pub fn mismatch(mut self, v: i32) -> Self {
        self.mismatch = v;
        self
    }
    /// blast-rs: SearchParams builder setter; not a direct NCBI C port.
    pub fn min_score(mut self, v: i32) -> Self {
        self.min_score = v;
        self
    }
    /// blast-rs: SearchParams builder setter; not a direct NCBI C port.
    pub fn query_gencode(mut self, v: u8) -> Self {
        self.query_gencode = v;
        self
    }
    /// blast-rs: SearchParams builder setter; not a direct NCBI C port.
    pub fn db_gencode(mut self, v: u8) -> Self {
        self.db_gencode = v;
        self
    }
    /// blast-rs: SearchParams builder setter; not a direct NCBI C port.
    pub fn db_length(mut self, v: i64) -> Self {
        self.db_length = v;
        self
    }
    /// blast-rs: SearchParams builder setter; not a direct NCBI C port.
    pub fn effective_search_space(mut self, v: i64) -> Self {
        self.effective_search_space = v;
        self
    }
    /// blast-rs: SearchParams builder setter; not a direct NCBI C port.
    pub fn max_intron_length(mut self, v: i32) -> Self {
        self.max_intron_length = v;
        self
    }
    /// blast-rs: SearchParams builder setter; not a direct NCBI C port.
    pub fn max_hsps(mut self, v: Option<usize>) -> Self {
        self.max_hsps = v;
        self
    }
    /// blast-rs: SearchParams builder setter; not a direct NCBI C port.
    pub fn culling_limit(mut self, v: Option<usize>) -> Self {
        self.culling_limit = v;
        self
    }
    /// blast-rs: SearchParams builder setter; not a direct NCBI C port.
    pub fn two_hit(mut self, v: bool) -> Self {
        self.two_hit = v;
        self
    }
    /// blast-rs: SearchParams builder setter; not a direct NCBI C port.
    pub fn two_hit_window(mut self, v: usize) -> Self {
        self.two_hit_window = v;
        self
    }
    /// blast-rs: SearchParams builder setter; not a direct NCBI C port.
    pub fn x_drop_ungapped(mut self, v: i32) -> Self {
        self.x_drop_ungapped = v;
        self
    }
    /// blast-rs: SearchParams builder setter; not a direct NCBI C port.
    pub fn x_drop_gapped(mut self, v: i32) -> Self {
        self.x_drop_gapped = v;
        self
    }
    /// blast-rs: SearchParams builder setter; not a direct NCBI C port.
    pub fn x_drop_final(mut self, v: i32) -> Self {
        self.x_drop_final = v;
        self
    }
    /// blast-rs: SearchParams builder setter; not a direct NCBI C port.
    pub fn chaining(mut self, v: bool) -> Self {
        self.chaining = v;
        self
    }
    /// blast-rs: SearchParams builder setter; not a direct NCBI C port.
    pub fn soft_masking(mut self, v: bool) -> Self {
        self.soft_masking = v;
        self
    }
    /// blast-rs: SearchParams builder setter; not a direct NCBI C port.
    pub fn lcase_masking(mut self, v: bool) -> Self {
        self.lcase_masking = v;
        self
    }
    /// blast-rs: SearchParams builder setter; not a direct NCBI C port.
    pub fn sum_stats(mut self, v: bool) -> Self {
        self.sum_stats = v;
        self
    }
}

// ── Database builder ────────────────────────────────────────────────────────

/// Entry for building a database.
#[derive(Debug, Clone)]
pub struct SequenceEntry {
    pub title: String,
    pub accession: String,
    pub sequence: Vec<u8>,
    pub taxid: Option<u32>,
}

/// Builder for creating BLAST databases.
pub struct BlastDbBuilder {
    pub seq_type: DbType,
    pub db_title: String,
    pub entries: Vec<SequenceEntry>,
}

impl BlastDbBuilder {
    /// blast-rs: Public BLAST database builder constructor; not a direct NCBI C port.
    pub fn new(seq_type: DbType, db_title: impl Into<String>) -> Self {
        BlastDbBuilder {
            seq_type,
            db_title: db_title.into(),
            entries: Vec::new(),
        }
    }

    /// blast-rs: Public BLAST database builder mutator; not a direct NCBI C port.
    pub fn add(&mut self, entry: SequenceEntry) {
        self.entries.push(entry);
    }

    /// blast-rs: Public BLAST database writer dispatcher; not a direct NCBI C port.
    pub fn write(&self, base_path: &Path) -> io::Result<()> {
        match self.seq_type {
            DbType::Nucleotide => self.write_nucleotide(base_path),
            DbType::Protein => self.write_protein(base_path),
        }
    }

    /// blast-rs: Writes the Rust v4 nucleotide BLAST database files; not a
    /// direct NCBI C port.
    fn write_nucleotide(&self, base_path: &Path) -> io::Result<()> {
        // Write .nsq
        let mut nsq = BufWriter::new(File::create(db_component_path(base_path, "nsq"))?);
        nsq.write_all(&[0u8])?; // sentinel

        let mut seq_offsets = vec![1u32];
        let mut amb_offsets = Vec::new();
        let mut total_length = 0u64;
        let mut max_seq_len = 0u32;

        for entry in &self.entries {
            let seq = &entry.sequence;
            total_length = total_length.saturating_add(seq.len() as u64);
            max_seq_len = max_seq_len.max(seq.len() as u32);
            let packed = encode_ncbi2na_sequence(seq);

            nsq.write_all(&packed)?;
            let seq_start = *seq_offsets.last().unwrap();
            let amb_offset = seq_start + packed.len() as u32;
            let ambiguity_data = encode_ncbi2na_ambiguity_data(seq);
            nsq.write_all(&ambiguity_data)?;
            amb_offsets.push(amb_offset);
            seq_offsets.push(amb_offset + ambiguity_data.len() as u32);
        }
        amb_offsets.push(*seq_offsets.last().unwrap_or(&0));
        nsq.flush()?;

        // Write .nhr
        let mut nhr = BufWriter::new(File::create(db_component_path(base_path, "nhr"))?);
        let mut hdr_offsets = vec![0u32];
        for (oid, entry) in self.entries.iter().enumerate() {
            let hdr = format!("{} {}", entry.accession, entry.title);
            let asn1 = encode_defline_asn1(&hdr, oid as i32);
            nhr.write_all(&asn1)?;
            hdr_offsets.push(hdr_offsets.last().unwrap() + asn1.len() as u32);
        }
        nhr.flush()?;

        // Write .nin
        write_index_file(
            &db_component_path(base_path, "nin"),
            4, // format version
            DbType::Nucleotide,
            &self.db_title,
            self.entries.len() as u32,
            total_length,
            max_seq_len,
            &hdr_offsets,
            &seq_offsets,
            Some(&amb_offsets),
        )
    }

    /// blast-rs: Writes the Rust v4 protein BLAST database files; not a direct
    /// NCBI C port.
    fn write_protein(&self, base_path: &Path) -> io::Result<()> {
        // Write .psq (protein sequences in NCBIstdaa)
        let mut psq = BufWriter::new(File::create(db_component_path(base_path, "psq"))?);
        psq.write_all(&[0u8])?; // sentinel

        let mut seq_offsets = vec![1u32];
        let mut total_length = 0u64;
        let mut max_seq_len = 0u32;
        for entry in &self.entries {
            total_length = total_length.saturating_add(entry.sequence.len() as u64);
            max_seq_len = max_seq_len.max(entry.sequence.len() as u32);
            let encoded = encode_ncbistdaa_sequence(&entry.sequence);
            psq.write_all(&encoded)?;
            psq.write_all(&[0u8])?; // sentinel between sequences
            let prev = *seq_offsets.last().unwrap();
            seq_offsets.push(prev + encoded.len() as u32 + 1);
        }
        psq.flush()?;

        // Write .phr
        let mut phr = BufWriter::new(File::create(db_component_path(base_path, "phr"))?);
        let mut hdr_offsets = vec![0u32];
        for (oid, entry) in self.entries.iter().enumerate() {
            let hdr = format!("{} {}", entry.accession, entry.title);
            let asn1 = encode_defline_asn1(&hdr, oid as i32);
            phr.write_all(&asn1)?;
            hdr_offsets.push(hdr_offsets.last().unwrap() + asn1.len() as u32);
        }
        phr.flush()?;

        // Write .pin
        write_index_file(
            &db_component_path(base_path, "pin"),
            4,
            DbType::Protein,
            &self.db_title,
            self.entries.len() as u32,
            total_length,
            max_seq_len,
            &hdr_offsets,
            &seq_offsets,
            None,
        )
    }
}

// ── Search functions ────────────────────────────────────────────────────────

struct BlastPreliminarySearchResult {
    prelim_evalue: f64,
    saved_ungapped_hits: Vec<crate::protein_lookup::ProteinHit>,
    preliminary_hsps: Vec<crate::protein_lookup::ProteinHit>,
}

/// NCBI: Blast_RunPreliminarySearch (`blast_engine.c`), specialized to one
/// protein query/subject pair after the AA word finder has produced init HSPs.
#[allow(clippy::too_many_arguments)]
fn blast_run_preliminary_search(
    query_aa: &[u8],
    gap_trigger_raw: i32,
    search_space: f64,
    subj_aa: &[u8],
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
    prot_kbp: &crate::stat::KarlinBlk,
    gumbel_blk: &Option<crate::stat::GumbelBlk>,
    gap_open: i32,
    gap_extend: i32,
    x_drop_gapped: i32,
    x_drop_final: i32,
    min_subject_length: i32,
    params: &SearchParams,
    ungapped_hits: Vec<crate::protein_lookup::ProteinHit>,
    scratch: &mut ProteinScratch,
) -> BlastPreliminarySearchResult {
    if ungapped_hits.is_empty() {
        return BlastPreliminarySearchResult {
            prelim_evalue: composition_prelim_evalue(params),
            saved_ungapped_hits: Vec::new(),
            preliminary_hsps: Vec::new(),
        };
    }
    let prelim_evalue = composition_prelim_evalue(params);
    let eval_cutoff = protein_eval_cutoff(
        prelim_evalue,
        &prot_kbp,
        gumbel_blk.as_ref(),
        query_aa.len() as i32,
        min_subject_length,
        search_space,
    );
    let adjusted_cutoff = gap_trigger_raw.min(eval_cutoff).max(1);
    // word_cutoff = MIN(gap_trigger, eval_cutoff) per NCBI
    // `BlastInitialWordParametersUpdate` (`blast_parameters.c:367`).
    let word_cutoff = adjusted_cutoff;
    // hit_cutoff = SpougeEtoS(eval, ...) per NCBI
    // `BlastHitSavingParametersUpdate` (`blast_parameters.c:939`).
    let hit_cutoff = eval_cutoff;
    let chained_kept = s_chaining_alignment(
        &ungapped_hits,
        gap_open,
        gap_extend,
        word_cutoff,
        hit_cutoff,
    );
    let chained_ungapped_hits: Vec<crate::protein_lookup::ProteinHit> = ungapped_hits
        .iter()
        .zip(chained_kept.iter())
        .filter(|(uh, keep)| uh.score >= adjusted_cutoff && **keep)
        .map(|(uh, _)| uh)
        .cloned()
        .collect();
    let preliminary_hsps = if params.ungapped {
        Vec::new()
    } else {
        protein_prelim_gapped_hits(
            query_aa,
            subj_aa,
            matrix,
            &chained_ungapped_hits,
            gap_open,
            gap_extend,
            x_drop_gapped,
            x_drop_final,
            adjusted_cutoff,
            hit_cutoff,
            false,
            scratch,
            true,
            0,
            0,
        )
    };
    BlastPreliminarySearchResult {
        prelim_evalue,
        saved_ungapped_hits: chained_ungapped_hits,
        preliminary_hsps,
    }
}

/// NCBI: Blast_RunFullSearch (`blast_engine.c`), specialized to the ordinary
/// protein path between preliminary gapped scoring and composition/e-value
/// processing.
#[allow(clippy::too_many_arguments)]
fn blast_run_full_search(
    query_aa: &[u8],
    subj_aa: &[u8],
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
    preliminary_hsps: Vec<crate::protein_lookup::ProteinHit>,
    gap_open: i32,
    gap_extend: i32,
    x_drop_final: i32,
    scratch: &mut ProteinScratch,
) -> Vec<crate::protein_lookup::ProteinHit> {
    let mut phits = preliminary_hsps;
    blast_traceback_from_protein_hit_hsp_list(
        query_aa,
        subj_aa,
        matrix,
        &mut phits,
        gap_open,
        gap_extend,
        x_drop_final,
        scratch,
    );
    for ph in &mut phits {
        ph.score = rescore_protein_hit(ph, query_aa, subj_aa, matrix, gap_open, gap_extend);
    }
    phits
}

/// NCBI: s_ProcessHSPList (`blast_engine.c`), specialized to final protein HSP
/// e-value filtering and API-HSP materialization after traceback/redo.
#[allow(clippy::too_many_arguments)]
fn s_process_hsp_list(
    final_phits: Vec<crate::protein_lookup::ProteinHit>,
    query_aa: &[u8],
    subj_aa: &[u8],
    prot_kbp: &crate::stat::KarlinBlk,
    gumbel_blk: &Option<crate::stat::GumbelBlk>,
    search_space: f64,
    evalue_threshold: f64,
    comp_mode: u8,
    comp_scale: f64,
    use_adj_matrix: bool,
    lambda_ratio_opt: Option<f64>,
    comp_adjust_method_id: u8,
) -> Vec<Hsp> {
    let sl = subj_aa.len();
    final_phits
        .into_iter()
        .filter_map(|mut ph| {
            // Mirror NCBI's `Blast_HSPListGetEvalues` (`blast_hits.c:1873`):
            // when comp_adjust is on, NCBI computes the e-value with the
            // *scaled* score and Lambda/scale_factor, then later rounds the
            // score back via `s_HSPListNormalizeScores`.
            let (e_score_i32, e_kbp_lambda) = match ph.scaled_score {
                Some(s) if comp_mode > 0 => (s, prot_kbp.lambda / comp_scale),
                _ => (ph.score, prot_kbp.lambda),
            };
            let e_kbp = crate::stat::KarlinBlk {
                lambda: e_kbp_lambda,
                k: prot_kbp.k,
                log_k: prot_kbp.log_k,
                h: prot_kbp.h,
                round_down: prot_kbp.round_down,
            };
            let evalue = if let Some(ref gbp) = gumbel_blk {
                let base_ev = crate::stat::blast_spouge_sto_e(
                    e_score_i32,
                    Some(&e_kbp),
                    Some(gbp),
                    query_aa.len() as i32,
                    sl as i32,
                );
                if use_adj_matrix {
                    base_ev
                } else if let Some(lr) = lambda_ratio_opt {
                    let scaled_kbp = crate::stat::KarlinBlk {
                        lambda: e_kbp.lambda / lr,
                        k: e_kbp.k,
                        log_k: e_kbp.log_k,
                        h: e_kbp.h,
                        round_down: e_kbp.round_down,
                    };
                    crate::stat::blast_spouge_sto_e(
                        e_score_i32,
                        Some(&scaled_kbp),
                        Some(gbp),
                        query_aa.len() as i32,
                        sl as i32,
                    )
                } else {
                    base_ev
                }
            } else {
                let raw_evalue = e_kbp.raw_to_evalue(e_score_i32, search_space);
                if use_adj_matrix {
                    raw_evalue
                } else if let Some(lr) = lambda_ratio_opt {
                    let scaled_lambda = e_kbp.lambda / lr;
                    search_space * e_kbp.k * (-scaled_lambda * e_score_i32 as f64).exp()
                } else {
                    raw_evalue
                }
            };
            if evalue > evalue_threshold {
                return None;
            }
            let (q_aln, s_aln) = match (ph.qseq.take(), ph.sseq.take()) {
                (Some(qs), Some(ss)) => (qs.into_bytes(), ss.into_bytes()),
                _ => {
                    let q_aln: Vec<u8> = (0..ph.align_length as usize)
                        .map(|i| {
                            let idx = ph.query_start + i;
                            if idx < query_aa.len() {
                                ncbistdaa_to_aminoacid_base(query_aa[idx])
                            } else {
                                b'-'
                            }
                        })
                        .collect();
                    let s_aln: Vec<u8> = (0..ph.align_length as usize)
                        .map(|i| {
                            let idx = ph.subject_start + i;
                            if idx < sl {
                                ncbistdaa_to_aminoacid_base(subj_aa[idx])
                            } else {
                                b'-'
                            }
                        })
                        .collect();
                    (q_aln, s_aln)
                }
            };
            Some(Hsp {
                score: ph.score,
                bit_score: prot_kbp.raw_to_bit(ph.score),
                evalue,
                query_start: ph.query_start,
                query_end: ph.query_end,
                subject_start: ph.subject_start,
                subject_end: ph.subject_end,
                query_gapped_start: ph.gapped_start_q,
                subject_gapped_start: ph.gapped_start_s,
                query_link_start: ph.query_start,
                query_link_end: ph.query_end,
                query_link_gapped_start: ph.gapped_start_q,
                subject_link_start: ph.subject_start,
                subject_link_end: ph.subject_end,
                subject_link_gapped_start: ph.gapped_start_s,
                link_score: ph.scaled_score,
                link_lambda: ph
                    .scaled_score
                    .map(|_| prot_kbp.lambda / COMPO_ADJUST_SCALE_FACTOR),
                num_identities: ph.num_ident as usize,
                num_gaps: ph.gap_opens as usize,
                alignment_length: ph.align_length as usize,
                query_aln: q_aln,
                midline: Vec::new(),
                subject_aln: s_aln,
                query_frame: 0,
                subject_frame: 0,
                num_links: 1,
                comp_adjust_method: comp_adjust_method_id,
            })
        })
        .collect()
}

/// Run a blastp search (protein query vs protein database).
/// blast-rs: Public blastp API pipeline assembled from ported lower-level pieces;
/// not a direct NCBI C port.
#[allow(clippy::too_many_arguments)]
fn process_protein_oid(
    query_aa: &[u8],
    ungapped_kbp: &crate::stat::KarlinBlk,
    gap_trigger_raw: i32,
    search_space: f64,
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
    prot_kbp: &crate::stat::KarlinBlk,
    gumbel_blk: &Option<crate::stat::GumbelBlk>,
    gap_open: i32,
    gap_extend: i32,
    x_drop_gapped: i32,
    x_drop_final: i32,
    min_subject_length: i32,
    evalue_threshold: f64,
    params: &SearchParams,
    db: &BlastDb,
    oid: u32,
    subj_aa: &[u8],
    ungapped_hits: Vec<crate::protein_lookup::ProteinHit>,
    scratch: &mut ProteinScratch,
) -> Option<SearchResult> {
    let subj_len = subj_aa.len();
    if ungapped_hits.is_empty() {
        return None;
    }

    // NCBI uses `cbs_stretch * evalue` (= 5*evalue when comp_adjust>1) for
    // preliminary composition-based searches (`blast_engine.c:653`).
    let preliminary = blast_run_preliminary_search(
        query_aa,
        gap_trigger_raw,
        search_space,
        subj_aa,
        matrix,
        prot_kbp,
        gumbel_blk,
        gap_open,
        gap_extend,
        x_drop_gapped,
        x_drop_final,
        min_subject_length,
        params,
        ungapped_hits,
        scratch,
    );
    // NCBI `-ungapped` mode (blastp/blastx/tblastn/tblastx): skip the
    // gapped DP entirely and emit each surviving ungapped HSP using
    // ungapped Karlin params for bit/evalue. Mirrors NCBI's
    // `s_BlastGetUngappedAlignmentResults` behavior — one HSP per
    // qualifying seed, no gapped-DP dedup.
    if params.ungapped {
        let mut ungapped_hsps: Vec<Hsp> = Vec::new();
        for uh in &preliminary.saved_ungapped_hits {
            // NCBI's `Blast_ScoreBlkMatrixInit` frees `sbp->gbp` when
            // `gapped_calculation=FALSE` (`blast_setup.c:524`), so the
            // ungapped path uses `BLAST_KarlinStoE_simple` (simple
            // Karlin with the per-query effective search space), NOT
            // Spouge FSC.
            let evalue = ungapped_kbp.raw_to_evalue(uh.score, search_space);
            if evalue > evalue_threshold {
                continue;
            }
            let q_slice = &query_aa[uh.query_start..uh.query_end];
            let s_slice = &subj_aa[uh.subject_start..uh.subject_end];
            let q_aln: Vec<u8> = q_slice
                .iter()
                .map(|&b| ncbistdaa_to_aminoacid_base(b))
                .collect();
            let s_aln: Vec<u8> = s_slice
                .iter()
                .map(|&b| ncbistdaa_to_aminoacid_base(b))
                .collect();
            ungapped_hsps.push(Hsp {
                score: uh.score,
                bit_score: ungapped_kbp.raw_to_bit(uh.score),
                evalue,
                query_start: uh.query_start,
                query_end: uh.query_end,
                subject_start: uh.subject_start,
                subject_end: uh.subject_end,
                query_gapped_start: uh.gapped_start_q,
                subject_gapped_start: uh.gapped_start_s,
                query_link_start: uh.query_start,
                query_link_end: uh.query_end,
                query_link_gapped_start: uh.gapped_start_q,
                subject_link_start: uh.subject_start,
                subject_link_end: uh.subject_end,
                subject_link_gapped_start: uh.gapped_start_s,
                link_score: None,
                link_lambda: None,
                num_identities: uh.num_ident as usize,
                num_gaps: 0,
                alignment_length: uh.align_length as usize,
                query_aln: q_aln,
                midline: Vec::new(),
                subject_aln: s_aln,
                query_frame: 0,
                subject_frame: 0,
                num_links: 1,
                comp_adjust_method: 0,
            });
        }
        if ungapped_hsps.is_empty() {
            return None;
        }
        let accession = db
            .get_accession(oid)
            .unwrap_or_else(|| format!("oid_{}", oid));
        let title = String::from_utf8_lossy(db.get_header(oid)).to_string();
        return Some(SearchResult {
            subject_oid: oid,
            subject_title: title,
            subject_accession: accession,
            subject_len: subj_aa.len(),
            hsps: ungapped_hsps,
            taxids: vec![],
        });
    }
    let phits = blast_run_full_search(
        query_aa,
        subj_aa,
        matrix,
        preliminary.preliminary_hsps,
        gap_open,
        gap_extend,
        x_drop_final,
        scratch,
    );
    if phits.is_empty() {
        return None;
    }

    // Pre-filter: check e-value with Spouge FSC if available, else simple
    // Karlin. NCBI uses `hit_params->prelim_evalue` here, which already
    // includes the composition-based stretch factor above.
    let best_raw_ev = phits
        .iter()
        .map(|ph| {
            if let Some(ref gbp) = gumbel_blk {
                crate::stat::blast_spouge_sto_e(
                    ph.score,
                    Some(&prot_kbp),
                    Some(gbp),
                    query_aa.len() as i32,
                    subj_len as i32,
                )
            } else {
                prot_kbp.raw_to_evalue(ph.score, search_space)
            }
        })
        // NCBI `s_BlastGetBestEvalue` (`blast_hits.c:1742`) seeds with `(double)INT4_MAX`.
        .fold(i32::MAX as f64, f64::min);
    if best_raw_ev > preliminary.prelim_evalue {
        return None;
    }

    let accession = db
        .get_accession(oid)
        .unwrap_or_else(|| format!("oid_{}", oid));
    let title = String::from_utf8_lossy(db.get_header(oid)).to_string();
    let sl = subj_aa.len();
    let matrix_name = protein_matrix_name(params.matrix);

    // Composition-based e-value adjustment (NCBI comp_based_stats).
    // Verbatim port of Blast_AdjustScores + s_AdjustEvaluesForComposition.
    let comp_mode = params.comp_adjust;
    let redo_subj_aa: &[u8] = if comp_mode > 0 {
        kappa_redo_subject_sequence(
            query_aa,
            subj_aa,
            &phits,
            prot_kbp.lambda,
            &mut scratch.mask_buf,
        )
    } else {
        subj_aa
    };
    let comp_scale = if comp_mode > 0 {
        COMPO_ADJUST_SCALE_FACTOR
    } else {
        1.0
    };
    let scaled_gap_open = crate::math::blast_nint(gap_open as f64 * comp_scale) as i32;
    let scaled_gap_extend = crate::math::blast_nint(gap_extend as f64 * comp_scale) as i32;
    let scaled_x_drop_final = crate::math::blast_nint(x_drop_final as f64 * comp_scale) as i32;

    // Determine adjustment rule and optionally build adjusted matrix.
    // adj_result: None = no adjustment, Some((adjusted_matrix_opt, lambda_ratio_opt))
    type AdjResult = Option<(Option<[[i32; AA_SIZE]; AA_SIZE]>, Option<f64>)>;
    // Track which adjustment rule was actually applied per subject; mirrors
    // NCBI's `BlastHSP::comp_adjustment_method` (`blast_kappa.c:332-342`)
    // for pairwise `, Method: ...` annotation.
    let mut comp_adjust_method_id: u8 = 0;
    let adj_result: AdjResult = if comp_mode > 0 {
        let (qcomp28, qn) = crate::composition::blast_read_aa_composition(&query_aa, AA_SIZE);
        let (scomp28, sn) = crate::composition::blast_read_aa_composition(redo_subj_aa, AA_SIZE);
        if qn == 0 || sn == 0 {
            Some((None, None))
        } else {
            let mut qp20 = [0.0f64; 20];
            let mut sp20 = [0.0f64; 20];
            crate::compo_mode_condition::s_gather_letter_probs(&qcomp28, &mut qp20);
            crate::compo_mode_condition::s_gather_letter_probs(&scomp28, &mut sp20);

            let rule = if comp_mode == 1 {
                crate::compo_mode_condition::MatrixAdjustRule::ScaleOldMatrix
            } else {
                crate::compo_mode_condition::blast_choose_matrix_adjust_rule(
                    query_aa.len(),
                    subj_aa.len(),
                    &qp20,
                    &sp20,
                    comp_mode,
                )
            };
            use crate::compo_mode_condition::MatrixAdjustRule::*;
            comp_adjust_method_id = match rule {
                DontAdjust => 0,
                ScaleOldMatrix => 1,
                UserSpecifiedRelEntropy
                | UnconstrainedRelEntropy
                | RelEntropyOldMatrixNewContext
                | RelEntropyOldMatrixOldContext => 2,
            };

            use crate::compo_mode_condition::MatrixAdjustRule;
            match rule {
                MatrixAdjustRule::DontAdjust => None,
                MatrixAdjustRule::ScaleOldMatrix => {
                    // Port of NCBI Blast_CompositionBasedStats: rescale matrix
                    // using composition-specific lambda ratio, then re-align.
                    let ungapped_lambda =
                        crate::stat::protein_ideal_ungapped_kbp_for_matrix(matrix_name).lambda
                            / comp_scale;
                    let Some(freq_ratios) = crate::matrix::get_matrix_freq_ratios(matrix_name)
                    else {
                        return None;
                    };
                    let start_matrix = crate::composition::blast_int4_matrix_from_freq(
                        ungapped_lambda,
                        &freq_ratios,
                    );
                    crate::composition::composition_scale_matrix_with_ratio(
                        &start_matrix,
                        &qcomp28,
                        &scomp28,
                        ungapped_lambda,
                        &freq_ratios,
                    )
                    .map(|(adj_mat, lambda_ratio)| (Some(adj_mat), Some(lambda_ratio)))
                }
                MatrixAdjustRule::UserSpecifiedRelEntropy
                | MatrixAdjustRule::UnconstrainedRelEntropy
                | MatrixAdjustRule::RelEntropyOldMatrixNewContext
                | MatrixAdjustRule::RelEntropyOldMatrixOldContext => {
                    // Full matrix optimization (Blast_CompositionMatrixAdj)
                    let (joint_probs, first_std, second_std) =
                        crate::composition::blosum62_workspace();
                    let mut adj_matrix = *matrix;
                    // NCBI uses ungappedLambda (0.3176 for BLOSUM62) for matrix scaling,
                    // NOT the gapped lambda. See matrixInfo->ungappedLambda.
                    let ungapped_lambda =
                        crate::stat::protein_ideal_ungapped_kbp_for_matrix(matrix_name).lambda
                            / comp_scale;
                    let Some(freq_ratios) = crate::matrix::get_matrix_freq_ratios(matrix_name)
                    else {
                        return None;
                    };
                    let start_matrix = crate::composition::blast_int4_matrix_from_freq(
                        ungapped_lambda,
                        &freq_ratios,
                    );
                    let status = crate::composition::blast_composition_matrix_adj(
                        &mut adj_matrix,
                        AA_SIZE,
                        rule,
                        qn,
                        sn,
                        &qcomp28,
                        &scomp28,
                        20,   // RE_pseudocounts
                        0.44, // kFixedReBlosum62
                        &joint_probs,
                        &first_std,
                        &second_std,
                        ungapped_lambda,
                        &start_matrix,
                    );
                    if status == 0 {
                        // Optimization succeeded — use adjusted matrix
                        Some((Some(adj_matrix), None))
                    } else {
                        // 1-1 with `Blast_AdjustScores` (composition_adjustment.c:1501–1530):
                        // when matrix optimization returns a non-fatal failure,
                        // NCBI resets the rule to `eCompoScaleOldMatrix` and
                        // calls `Blast_CompositionBasedStats` to rescale the
                        // matrix using lambda_ratio. We mirror that by
                        // invoking `composition_scale_matrix` rather than
                        // returning bare `lambda_ratio`.
                        comp_adjust_method_id = 1;
                        crate::composition::composition_scale_matrix_with_ratio(
                            &start_matrix,
                            &qcomp28,
                            &scomp28,
                            ungapped_lambda,
                            &freq_ratios,
                        )
                        .map(|(adj_mat, lambda_ratio)| (Some(adj_mat), Some(lambda_ratio)))
                    }
                }
            }
        }
    } else {
        None
    };

    // If we have an adjusted matrix, re-align and recompute scores.
    // NCBI uses the SW-bounded forward traceback (`s_SWFindFinalEndsUsingXdrop`)
    // ONLY in Smith-Waterman traceback mode (`-use_sw_tback`). The DEFAULT
    // composition-adjusted path (`comp_based_stats` 1/2, no SW) re-aligns
    // BIDIRECTIONALLY from the seed with the adjusted matrix via
    // `protein_gapped_align`. Forcing the SW-bounded forward path here
    // produced a different (suboptimal) co-optimal gap placement than NCBI —
    // correct score but wrong alignment/identity (e.g. pident 32% vs NCBI's
    // 68% on HIV env hits), corrupting `pident`/`mismatch`/`nident` and
    // shifting recall on every default blastp/blastx/tblastn run. blast-rs
    // does not implement `-use_sw_tback`, so the bidirectional path is always
    // correct here.
    let sw_bounded_blastp: Option<(usize, usize, usize, usize, i32)> = None;
    let (final_phits, use_adj_matrix) = if let Some((None, _)) = adj_result {
        (Vec::new(), true)
    } else if let Some((Some(ref adj_mat), _)) = adj_result {
        // Re-do gapped alignment with adjusted matrix
        let mut new_phits = Vec::new();
        for ph in &phits {
            // In the single-HSP Smith-Waterman redo arm, NCBI takes the
            // SW-derived bounds through `s_SWFindFinalEndsUsingXdrop` and
            // converts that traceback directly. Only fall back to the
            // ordinary bidirectional redo when the bounded path is not
            // available for this subject/HSP set.
            let bounded_gr = if let Some((q_start, m_start, q_extent, s_extent, target_score)) =
                sw_bounded_blastp
            {
                crate::protein::s_sw_find_final_ends_using_xdrop(
                    &query_aa,
                    redo_subj_aa,
                    q_start,
                    m_start,
                    q_extent,
                    s_extent,
                    target_score,
                    adj_mat,
                    scaled_gap_open,
                    scaled_gap_extend,
                    scaled_x_drop_final,
                )
            } else {
                None
            };
            // NCBI's `s_RedoOneAlignment` (`blast_kappa.c:1924-1926`) uses
            // `hsp->query.gapped_start` (the seed midpoint) — NOT the
            // alignment's `query.offset` — when re-running the rescaled
            // gapped DP. Using `query_start` (alignment left edge) here
            // shifts the WHOLE alignment off-center: with X-drop=2078 in
            // scaled units, both directions extend further, and
            // (`a_offset`, `b_offset`) lock onto a maximum cell that
            // produces longer-but-suboptimal bounds vs NCBI's
            // seed-centered DP.
            let gr_opt = if let Some(b) = bounded_gr {
                Some(b)
            } else {
                crate::protein::protein_gapped_align(
                    &query_aa,
                    redo_subj_aa,
                    ph.gapped_start_q,
                    ph.gapped_start_s,
                    adj_mat,
                    scaled_gap_open,
                    scaled_gap_extend,
                    scaled_x_drop_final,
                )
            };
            if let Some(gr) = gr_opt {
                let q_slice = &query_aa[gr.query_start..gr.query_end];
                // NCBI masks the subject ONLY for scoring/edit-script (the
                // composition-adjusted `gr.score` is correct and must NOT
                // change). The DISPLAYED sequence and identities are rebuilt
                // from the UNMASKED subject (blast_kappa.c: sseq rebuilt from
                // the real DB residues during traceback, num_ident recomputed
                // from those). Render from `subj_aa` (unmasked), not the
                // SEG-masked `redo_subj_aa`, so `sseq` shows the real residues
                // and `pident`/`mismatch`/`nident` reflect them.
                let s_slice = &subj_aa[gr.subject_start..gr.subject_end];
                let (qs, ss) =
                    gr.edit_script
                        .render_alignment(q_slice, s_slice, ncbistdaa_to_aminoacid_char);
                let (align_length, num_ident, gap_opens) =
                    gr.edit_script.count_identities(q_slice, s_slice);
                let mismatches = (align_length - num_ident - gap_opens).max(0);
                new_phits.push(crate::protein_lookup::ProteinHit {
                    query_start: gr.query_start,
                    query_end: gr.query_end,
                    subject_start: gr.subject_start,
                    subject_end: gr.subject_end,
                    score: crate::math::blast_nint(gr.score as f64 / comp_scale) as i32,
                    num_ident,
                    align_length,
                    mismatches,
                    gap_opens,
                    qseq: Some(qs),
                    sseq: Some(ss),
                    scaled_score: Some(gr.score),
                    adjusted_evalue: None,
                    gapped_start_q: ph.gapped_start_q,
                    gapped_start_s: ph.gapped_start_s,
                });
            }
        }
        (new_phits, true)
    } else {
        (phits.clone(), false)
    };

    let lambda_ratio_opt = adj_result.as_ref().and_then(|(_, lr)| *lr);

    let hsps = s_process_hsp_list(
        final_phits,
        query_aa,
        subj_aa,
        prot_kbp,
        gumbel_blk,
        search_space,
        evalue_threshold,
        comp_mode,
        comp_scale,
        use_adj_matrix,
        lambda_ratio_opt,
        comp_adjust_method_id,
    );

    if hsps.is_empty() {
        return None;
    }
    Some(SearchResult {
        subject_oid: oid,
        subject_title: title,
        subject_accession: accession,
        subject_len: sl,
        hsps,
        taxids: vec![],
    })
}

pub fn blastp(db: &BlastDb, query: &[u8], params: &SearchParams) -> Vec<SearchResult> {
    if query.is_empty() {
        return Vec::new();
    }

    let query_aa = encode_ncbistdaa_sequence(query);
    let query_aa_masked = if params.filter_low_complexity {
        encode_protein_query(
            query,
            true,
            params.seg_window,
            params.seg_locut,
            params.seg_hicut,
        )
    } else {
        query_aa.clone()
    };
    if crate::composition::blast_read_aa_composition(&query_aa_masked, AA_SIZE).1 == 0 {
        return Vec::new();
    }

    let matrix = *get_matrix(params.matrix);
    let word_size = params.word_size.clamp(2, 7);
    let threshold = params.word_threshold.unwrap_or_else(|| {
        suggested_word_threshold(params.matrix, crate::program::BLASTP, params.word_size)
    });

    let prot_kbp = protein_kbp_for_matrix(params.matrix, params.gap_open, params.gap_extend);

    // Convert bit-score x_dropoff values to raw scores.
    // NCBI uses UNGAPPED KBP for ungapped x_drop, GAPPED KBP for gapped x_drop.
    // Ungapped path: `ceil()` then Int4 cast (`blast_parameters.c:221`).
    // Gapped path: plain `(Int4)` cast — truncation toward zero
    // (`blast_parameters.c:457-463`).
    // NCBI's `Blast_ScoreBlkKbpUngappedCalc` (`blast_stat.c:2737`) populates
    // `sbp->kbp_std[context]` from the **query's actual amino-acid composition**,
    // not the ideal Robinson background. For most queries the result is close
    // to the ideal, but the small drift in lambda/logK shifts `gap_trigger` by
    // 1 raw score unit — enough to swing boundary hits like seqp's DAA02208
    // (max ungapped score 40 vs ideal cutoff 41 vs query-specific cutoff 40).
    let matrix_name = protein_matrix_name(params.matrix);
    let ungapped_kbp = crate::stat::query_specific_protein_ungapped_kbp_for_matrix(
        &query_aa,
        matrix_name,
        &matrix,
    );
    let x_drop_ungapped = raw_ungapped_xdrop_bits(params.x_drop_ungapped, &ungapped_kbp);
    // NCBI's PRELIMINARY gapped extension uses `gap_x_dropoff` (15 bits for
    // protein default) — `blast_parameters.c:457-463` truncates with `(Int4)`.
    // The TRACEBACK phase then re-extends with the larger `gap_x_dropoff_final`
    // (25 bits). Using `x_drop_final` here makes our preliminary gapped DP
    // far more permissive than NCBI's, finding gap-rich alignments NCBI rejects
    // (e.g. NP_777001 score=50 vs NCBI's 46).
    let x_drop_gapped = raw_gapped_xdrop_bits(params.x_drop_gapped, &prot_kbp);
    let x_drop_final = raw_gapped_xdrop_bits(params.x_drop_final, &prot_kbp);

    // Gap trigger: NCBI `blast_parameters.c:343`:
    //   (Int4)((gap_trigger * NCBIMATH_LN2 + kbp->logK) / kbp->Lambda)
    // where `kbp` is the UNGAPPED KBP (per-context, derived from query composition).
    let gap_trigger_raw = raw_gap_trigger_bits(crate::stat::BLAST_GAP_TRIGGER_PROT, &ungapped_kbp);

    let total_subj_len: usize = (0..db.num_oids)
        .map(|oid| db.get_seq_len(oid) as usize)
        .sum();
    let stats_db_len = statistical_db_length(params, total_subj_len as i64);
    // NCBI uses MIN subject length (not average) for `cutoff_score_max`
    // calculation when gbp is filled (`blast_setup.c:970`). Compute it here
    // and reuse for the per-OID seed cutoff.
    // NCBI's `BlastSeqSrcGetMinSeqLen` returns the DB's `m_MinLen` from
    // metadata, defaulting to `BLAST_SEQSRC_MINLENGTH = 10` when the
    // metadata isn't stored (V4 DBs don't store min_seq_len, only
    // max_seq_len — see `seqdbimpl.cpp:131-136`). Our V4 DB matches:
    // we don't have a stored min_seq_len, so faithfully default to 10.
    // (Computing actual min from scan diverges from NCBI: e.g. seqp's
    // shortest seq is 7 aa but NCBI uses 10, giving different
    // SpougeEtoS cutoffs.)
    const BLAST_SEQSRC_MINLENGTH: i32 = 10;
    let min_subject_length: i32 = BLAST_SEQSRC_MINLENGTH;

    let (_len_adj, search_space) = protein_api_search_space(
        params,
        query_aa.len(),
        stats_db_len,
        db.num_oids as i32,
        &prot_kbp,
    );

    // Build Gumbel block for Spouge FSC e-value (per-subject length correction).
    // When the user overrides `-searchsp`, NCBI bypasses the per-subject Spouge
    // calculation and uses the simple Karlin formula `K * S * exp(-Lambda*y)`
    // so the user-supplied search space drives the e-value directly. We mirror
    // this by zeroing the Gumbel block, which routes every per-pair evalue
    // through the `prot_kbp.raw_to_evalue(score, search_space)` branch.
    let gumbel_blk = if params.effective_search_space > 0 {
        None
    } else {
        protein_gumbel_for_matrix(
            params.matrix,
            params.gap_open,
            params.gap_extend,
            stats_db_len,
        )
    };

    // Build lookup table once per query (not per subject).
    let lookup_table = crate::protein_lookup::ProteinLookupTable::build(
        &query_aa_masked,
        word_size,
        &matrix,
        threshold,
    );

    let max_hsps = params.max_hsps;
    let evalue_threshold = params.evalue_threshold;
    let gap_open = params.gap_open;
    let gap_extend = params.gap_extend;
    let two_hit_window = params.two_hit_window as i32;

    // NCBI's per-context ungapped `cutoffs->cutoff_score` is
    // `MIN(gap_trigger, cutoff_score_max)` (`blast_parameters.c:367`), which is
    // exactly `protein_prelim_seed_cutoff` here. `s_BlastAaWordFinder_TwoHit`
    // (`aa_ungapped.c:588`) only saves init HSPs scoring `>= cutoff_score`. All
    // inputs (gap_trigger, prot_kbp, gumbel, min_subject_length, search_space)
    // are per-query, so this is computed once and threaded into the scan. Using
    // the same `composition_prelim_evalue` that `process_protein_oid` uses keeps
    // the finder gate identical to the later gapped filter.
    let ungapped_seed_cutoff = protein_prelim_seed_cutoff(
        gap_trigger_raw,
        composition_prelim_evalue(params),
        &prot_kbp,
        gumbel_blk.as_ref(),
        query_aa.len() as i32,
        min_subject_length,
        search_space,
    );

    // Process a single subject OID — extracted so it can be called
    // either sequentially or in parallel without per-call pool overhead.
    // `scratch` is a per-thread reusable scratch (see [`ProteinScratch`]).
    let search_oid = |scratch: &mut ProteinScratch, oid: u32| -> Option<SearchResult> {
        let subj_raw = db.get_sequence(oid);
        let subj_len = db.get_seq_len(oid) as usize;
        if subj_len < word_size {
            return None;
        }
        let subj_aa = &subj_raw[..subj_len];
        scratch.diag_buf.clear();
        let ungapped_hits = crate::protein_lookup::protein_scan_with_table_reuse(
            &query_aa,
            subj_aa,
            &matrix,
            &lookup_table,
            x_drop_ungapped,
            two_hit_window,
            ungapped_seed_cutoff,
            &mut scratch.diag_buf,
        );
        process_protein_oid(
            &query_aa,
            &ungapped_kbp,
            gap_trigger_raw,
            search_space,
            &matrix,
            &prot_kbp,
            &gumbel_blk,
            gap_open,
            gap_extend,
            x_drop_gapped,
            x_drop_final,
            min_subject_length,
            evalue_threshold,
            params,
            db,
            oid,
            subj_aa,
            ungapped_hits,
            scratch,
        )
    };

    // Run sequentially or in parallel depending on num_threads/pool. Each
    // worker thread gets one [`ProteinScratch`] reused across all subjects
    // it processes (NCBI's per-thread `BlastGapAlignStruct` pattern).
    let mut results: Vec<SearchResult> =
        map_database_oids_init(db, params, ProteinScratch::new, |scratch, oid| {
            search_oid(scratch, oid)
        })
        .into_iter()
        .flatten()
        .collect();

    apply_api_min_score_filter(&mut results, params.min_score);
    if let Some(culling_limit) = params.culling_limit {
        apply_api_culling_limit(&mut results, culling_limit, crate::program::BLASTP);
    }
    apply_api_max_hsps_limit(&mut results, max_hsps);
    results.sort_by(compare_search_results);
    if results.len() > params.max_target_seqs {
        results.truncate(params.max_target_seqs);
    }
    let midline_style = if params.ungapped {
        ProteinMidlineStyle::PositiveMatrix(&matrix)
    } else {
        ProteinMidlineStyle::AnyMismatchPlus
    };
    fill_protein_midlines(&mut results, midline_style);
    results
}

/// Run a batch blastp search: multiple protein queries vs one protein database.
///
/// This is much more efficient than calling `blastp()` per query because
/// subjects are scanned once and checked against all query lookup tables.
/// Each subject is loaded into cache once, then all queries are checked.
///
/// Returns one `Vec<SearchResult>` per query, in the same order as `queries`.
/// blast-rs: Batched public blastp API pipeline; not a direct NCBI C port.
///
/// Follows NCBI's multi-query structure: one concatenated query (single
/// inter-query sentinel, `QueryInfo::new_blastp` offset layout) + one lookup
/// table, scanned once per subject, with hits attributed via
/// `bsearch_context_info` (NCBI `BSearchContextInfo`).
///
/// The shared scan uses a single `scan_x_drop` (max over queries). NCBI's
/// word-params x-dropoff is likewise one per-search value, so this matches NCBI
/// rather than blast-rs's own per-context x-drop refinement on the single-query
/// path.
pub fn blastp_batch(
    db: &BlastDb,
    queries: &[&[u8]],
    params: &SearchParams,
) -> Vec<Vec<SearchResult>> {
    if queries.is_empty() {
        return Vec::new();
    }
    // comp-based stats / SEG masking fall back to the per-query path.
    if params.filter_low_complexity || params.comp_adjust != 0 {
        return queries.iter().map(|q| blastp(db, q, params)).collect();
    }

    let matrix = *get_matrix(params.matrix);
    let matrix_name = protein_matrix_name(params.matrix);
    let word_size = params.word_size.clamp(2, 7);
    let threshold = params.word_threshold.unwrap_or_else(|| {
        suggested_word_threshold(params.matrix, crate::program::BLASTP, params.word_size)
    });
    let prot_kbp = protein_kbp_for_matrix(params.matrix, params.gap_open, params.gap_extend);
    let x_drop_gapped = raw_gapped_xdrop_bits(params.x_drop_gapped, &prot_kbp);
    let x_drop_final = raw_gapped_xdrop_bits(params.x_drop_final, &prot_kbp);
    let gap_open = params.gap_open;
    let gap_extend = params.gap_extend;
    let evalue_threshold = params.evalue_threshold;
    let max_hsps = params.max_hsps;
    let two_hit_window = params.two_hit_window as i32;

    let total_subj_len: usize = (0..db.num_oids)
        .map(|oid| db.get_seq_len(oid) as usize)
        .sum();
    let stats_db_len = statistical_db_length(params, total_subj_len as i64);
    const BLAST_SEQSRC_MINLENGTH: i32 = 10;
    let min_subject_length: i32 = BLAST_SEQSRC_MINLENGTH;
    let gumbel_blk = if params.effective_search_space > 0 {
        None
    } else {
        protein_gumbel_for_matrix(
            params.matrix,
            params.gap_open,
            params.gap_extend,
            stats_db_len,
        )
    };

    struct PreparedQuery {
        aa: Vec<u8>,
        ungapped_kbp: crate::stat::KarlinBlk,
        search_space: f64,
        gap_trigger_raw: i32,
        x_drop_ungapped: i32,
    }
    let prepared: Vec<PreparedQuery> = queries
        .iter()
        .map(|q| {
            let aa = encode_ncbistdaa_sequence(q);
            let ungapped_kbp = crate::stat::query_specific_protein_ungapped_kbp_for_matrix(
                &aa,
                matrix_name,
                &matrix,
            );
            let x_drop_ungapped = raw_ungapped_xdrop_bits(params.x_drop_ungapped, &ungapped_kbp);
            let gap_trigger_raw =
                raw_gap_trigger_bits(crate::stat::BLAST_GAP_TRIGGER_PROT, &ungapped_kbp);
            let search_space = if params.effective_search_space > 0 {
                params.effective_search_space as f64
            } else {
                protein_api_search_space(
                    params,
                    aa.len(),
                    stats_db_len,
                    db.num_oids as i32,
                    &prot_kbp,
                )
                .1
            };
            PreparedQuery {
                aa,
                ungapped_kbp,
                search_space,
                gap_trigger_raw,
                x_drop_ungapped,
            }
        })
        .collect();
    let num_queries = queries.len();

    // True concatenated scan, following NCBI's multi-query model: all queries
    // live in ONE concatenated buffer separated by a single inter-query sentinel,
    // indexed by ONE lookup table, scanned once per subject. The offset layout is
    // `QueryInfo::new_blastp` (offset += len + 1), so a hit's concatenated offset
    // maps back to its query via `bsearch_context_info` — the same machinery NCBI
    // uses (`BSearchContextInfo`). `PROTEIN_CONCAT_SENTINEL` (>= AA_SIZE) makes the
    // builder skip boundary-spanning words, and the ungapped extender stops at it,
    // so each HSP's bounds are identical to a standalone per-query scan.
    let scan_x_drop = prepared
        .iter()
        .map(|p| p.x_drop_ungapped)
        .max()
        .unwrap_or(0);
    let lengths: Vec<usize> = prepared.iter().map(|p| p.aa.len()).collect();
    let query_info = crate::queryinfo::QueryInfo::new_blastp(&lengths);
    let mut concat: Vec<u8> =
        vec![crate::protein_lookup::PROTEIN_CONCAT_SENTINEL; query_info.seq_buf_len()];
    for (qi, p) in prepared.iter().enumerate() {
        let off = query_info.contexts[qi].query_offset as usize;
        concat[off..off + p.aa.len()].copy_from_slice(&p.aa);
    }
    let context_starts: Vec<usize> = query_info
        .contexts
        .iter()
        .map(|ctx| ctx.query_offset as usize)
        .collect();
    let table =
        crate::protein_lookup::ProteinLookupTable::build(&concat, word_size, &matrix, threshold);

    // NCBI's two-hit word finder gates each saved init HSP on the per-context
    // `cutoffs->cutoff_score` (`aa_ungapped.c:588`). This concatenated scan uses
    // ONE cutoff for all contexts, so use the minimum per-query ungapped seed
    // cutoff (NCBI's `cutoff_score_min`); each query is re-filtered downstream by
    // its own cutoff inside `process_protein_oid`.
    let scan_cutoff_min = prepared
        .iter()
        .map(|p| {
            protein_prelim_seed_cutoff(
                p.gap_trigger_raw,
                composition_prelim_evalue(params),
                &prot_kbp,
                gumbel_blk.as_ref(),
                p.aa.len() as i32,
                min_subject_length,
                p.search_space,
            )
        })
        .min()
        .unwrap_or(1)
        .max(1);

    let process = |scratch: &mut ProteinScratch, oid: u32| -> Vec<(usize, SearchResult)> {
        let subj_raw = db.get_sequence(oid);
        let subj_len = db.get_seq_len(oid) as usize;
        if subj_len < word_size {
            return Vec::new();
        }
        let subj_aa = &subj_raw[..subj_len];
        scratch.diag_buf.clear();
        let hits = crate::protein_lookup::protein_scan_with_table_reuse_contexts(
            &concat,
            subj_aa,
            &matrix,
            &table,
            scan_x_drop,
            two_hit_window,
            scan_cutoff_min,
            &mut scratch.diag_buf,
            &context_starts,
        );
        if hits.is_empty() {
            return Vec::new();
        }
        let mut per_query: Vec<Vec<crate::protein_lookup::ProteinHit>> =
            (0..num_queries).map(|_| Vec::new()).collect();
        for mut h in hits {
            // Map the concatenated query offset back to its query (NCBI
            // BSearchContextInfo); blastp has one context per query.
            let qi = crate::queryinfo::bsearch_context_info(h.query_start as i32, &query_info);
            if qi < 0 || qi as usize >= num_queries {
                continue;
            }
            let qi = qi as usize;
            let off = query_info.contexts[qi].query_offset as usize;
            h.query_start -= off;
            h.query_end = h.query_end.saturating_sub(off);
            h.gapped_start_q = h.gapped_start_q.saturating_sub(off);
            per_query[qi].push(h);
        }
        let mut out = Vec::new();
        for qi in 0..num_queries {
            if per_query[qi].is_empty() {
                continue;
            }
            let p = &prepared[qi];
            if let Some(sr) = process_protein_oid(
                &p.aa,
                &p.ungapped_kbp,
                p.gap_trigger_raw,
                p.search_space,
                &matrix,
                &prot_kbp,
                &gumbel_blk,
                gap_open,
                gap_extend,
                x_drop_gapped,
                x_drop_final,
                min_subject_length,
                evalue_threshold,
                params,
                db,
                oid,
                subj_aa,
                std::mem::take(&mut per_query[qi]),
                scratch,
            ) {
                out.push((qi, sr));
            }
        }
        out
    };

    let all_hits: Vec<Vec<(usize, SearchResult)>> =
        map_database_oids_init(db, params, ProteinScratch::new, process);

    let mut results: Vec<Vec<SearchResult>> = (0..num_queries).map(|_| Vec::new()).collect();
    for oid_hits in all_hits {
        for (qi, sr) in oid_hits {
            results[qi].push(sr);
        }
    }
    for r in &mut results {
        apply_api_min_score_filter(r, params.min_score);
        if let Some(culling_limit) = params.culling_limit {
            apply_api_culling_limit(r, culling_limit, crate::program::BLASTP);
        }
        apply_api_max_hsps_limit(r, max_hsps);
        r.sort_by(compare_search_results);
        if r.len() > params.max_target_seqs {
            r.truncate(params.max_target_seqs);
        }
        let midline_style = if params.ungapped {
            ProteinMidlineStyle::PositiveMatrix(&matrix)
        } else {
            ProteinMidlineStyle::AnyMismatchPlus
        };
        fill_protein_midlines(r, midline_style);
    }
    results
}

/// Run a blastn search (nucleotide query vs nucleotide database).
/// blast-rs: Public blastn API pipeline assembled from ported lower-level pieces;
/// not a direct NCBI C port.
pub fn blastn_search(db: &BlastDb, query: &[u8], params: &SearchParams) -> Vec<SearchResult> {
    if query.is_empty() {
        return Vec::new();
    }

    let query_plus_nomask = encode_blastna_sequence(query);
    let mut query_plus_lookup = query_plus_nomask.clone();
    apply_blastn_query_lookup_masks(
        &mut query_plus_lookup,
        query,
        params.filter_low_complexity,
        params.lcase_masking,
    );
    let query_minus_lookup = reverse_complement_blastna_sequence(&query_plus_lookup);
    let query_minus_nomask = reverse_complement_blastna_sequence(&query_plus_nomask);

    let reward = params.match_score;
    let penalty = params.mismatch;

    let total_subj_len = db.total_length;
    let (kbp, search_space, len_adj) = blastn_api_stats(
        params,
        &query_plus_nomask,
        total_subj_len as i64,
        db.num_oids as i32,
    );

    let x_dropoff = params.x_drop_ungapped;

    let (q_plus_lookup, q_minus_lookup) = match params.strand.as_str() {
        "plus" => (query_plus_lookup.as_slice(), &[] as &[u8]),
        "minus" => (&[] as &[u8], query_minus_lookup.as_slice()),
        _ => (query_plus_lookup.as_slice(), query_minus_lookup.as_slice()),
    };
    let (q_plus_nomask, q_minus_nomask) = match params.strand.as_str() {
        "plus" => (query_plus_nomask.as_slice(), &[] as &[u8]),
        "minus" => (&[] as &[u8], query_minus_nomask.as_slice()),
        _ => (query_plus_nomask.as_slice(), query_minus_nomask.as_slice()),
    };
    let (q_plus_runtime, q_minus_runtime) = if params.soft_masking {
        (q_plus_nomask, q_minus_nomask)
    } else {
        (q_plus_lookup, q_minus_lookup)
    };
    let prepared_query = crate::search::PreparedBlastnQuery::new_megablast_with_nomask(
        q_plus_lookup,
        q_minus_lookup,
        q_plus_runtime,
        q_minus_runtime,
        params.word_size,
    );

    let search_oid = |last_hit_scratch: &mut Vec<crate::search::PackedDiagScratch>,
                      oid: u32|
     -> Option<SearchResult> {
        let subject_packed = db.get_sequence(oid);
        let subject_len = db.get_seq_len(oid) as usize;
        if subject_len < params.word_size {
            return None;
        }

        let mut hsps = crate::search::blastn_gapped_search_packed_prepared_with_xdrops(
            &prepared_query,
            q_plus_runtime,
            q_minus_runtime,
            subject_packed,
            subject_len,
            reward,
            penalty,
            params.gap_open,
            params.gap_extend,
            x_dropoff,
            params.x_drop_gapped,
            params.x_drop_final,
            &kbp,
            search_space,
            params.evalue_threshold,
            last_hit_scratch,
        );

        if hsps.is_empty() {
            return None;
        }
        if params.sum_stats && hsps.len() > 1 {
            apply_blastn_linked_sum_stats_to_search_hsps(
                &mut hsps,
                query.len() as i32,
                subject_len as i32,
                &kbp,
                &kbp,
                search_space,
                search_space,
                len_adj,
                len_adj,
            );
        }

        let accession = db
            .get_accession(oid)
            .unwrap_or_else(|| format!("oid_{}", oid));
        let title = String::from_utf8_lossy(db.get_header(oid)).to_string();

        let api_hsps: Vec<Hsp> = hsps
            .iter()
            .map(|h| Hsp {
                score: h.score,
                bit_score: h.bit_score,
                evalue: h.evalue,
                query_start: h.query_start as usize,
                query_end: h.query_end as usize,
                subject_start: h.subject_start as usize,
                subject_end: h.subject_end as usize,
                query_gapped_start: h.query_gapped_start.max(0) as usize,
                subject_gapped_start: h.subject_gapped_start.max(0) as usize,
                query_link_start: h.query_start.max(0) as usize,
                query_link_end: h.query_end.max(0) as usize,
                query_link_gapped_start: h.query_gapped_start.max(0) as usize,
                subject_link_start: h.subject_start.max(0) as usize,
                subject_link_end: h.subject_end.max(0) as usize,
                subject_link_gapped_start: h.subject_gapped_start.max(0) as usize,
                link_score: None,
                link_lambda: None,
                num_identities: h.num_ident as usize,
                num_gaps: h.gap_opens as usize,
                alignment_length: h.align_length as usize,
                query_aln: h
                    .qseq
                    .as_ref()
                    .map(|s| s.as_bytes().to_vec())
                    .unwrap_or_default(),
                subject_aln: h
                    .sseq
                    .as_ref()
                    .map(|s| s.as_bytes().to_vec())
                    .unwrap_or_default(),
                midline: h
                    .qseq
                    .as_deref()
                    .unwrap_or("")
                    .bytes()
                    .zip(h.sseq.as_deref().unwrap_or("").bytes())
                    .map(|(q, s)| if q == s { b'|' } else { b' ' })
                    .collect(),
                query_frame: if h.context == 1 { -1 } else { 1 },
                subject_frame: 0,
                num_links: 1,
                comp_adjust_method: 0,
            })
            .collect();

        Some(SearchResult {
            subject_oid: oid,
            subject_title: title,
            subject_accession: accession,
            subject_len,
            hsps: api_hsps,
            taxids: db.get_taxids(oid),
        })
    };

    // Wire NCBI `-window_size` into the contiguous-blastn two-hit gate. The gate
    // lives in `search.rs` (`diag_initial_hit_core_packed`) and reads a
    // per-thread window via `blastn_window_size()`; window 0 (the blastn/
    // megablast default) keeps the one-hit path byte-identical. Because the scan
    // is parallelized per-OID across rayon workers, the per-thread value must be
    // set on each worker — the per-thread scratch `init` closure runs once on
    // every worker, so set it there; the main thread is also covered for the
    // single-threaded path. Restore/clear after the search completes.
    let window_size = params.window_size.max(0);
    let prev_window = crate::search::set_blastn_window_size(window_size);

    // One reusable diagonal-scratch per worker thread (NCBI reuses `aux_struct->ewp`
    // across all subjects instead of reallocating per OID); the scan resizes/clears
    // it per call, so reuse is byte-identical.
    let mut results: Vec<SearchResult> = map_database_oids_init(
        db,
        params,
        || {
            // Each rayon worker thread sets its own copy of the window-size
            // thread-local before allocating scratch (thread-locals are not
            // inherited from the spawning thread).
            crate::search::set_blastn_window_size(window_size);
            prepared_query.last_hit_scratch(prepared_query.use_diag_hash())
        },
        |scratch, oid| search_oid(scratch, oid),
    )
    .into_iter()
    .flatten()
    .collect();

    // Restore the previous window on this thread and drop any two-hit companion
    // store left over from the scan.
    crate::search::set_blastn_window_size(prev_window);
    crate::search::clear_blastn_two_hit_store();

    apply_api_min_score_filter(&mut results, params.min_score);
    if let Some(culling_limit) = params.culling_limit {
        apply_api_culling_limit(&mut results, culling_limit, crate::program::BLASTN);
    }
    apply_api_max_hsps_limit(&mut results, params.max_hsps);
    results.sort_by(compare_search_results);
    if results.len() > params.max_target_seqs {
        results.truncate(params.max_target_seqs);
    }
    results
}

/// blast-rs: Public API blastn lookup masking helper; not a direct NCBI C port.
fn apply_blastn_query_lookup_masks(
    encoded_query: &mut [u8],
    raw_query: &[u8],
    filter_low_complexity: bool,
    lcase_masking: bool,
) {
    if filter_low_complexity {
        let mask = crate::filter::dust_filter(encoded_query, 20, 64, 1);
        mask.apply(encoded_query, crate::filter::K_NUCL_MASK);
    }
    if lcase_masking {
        for (encoded, raw) in encoded_query.iter_mut().zip(raw_query.iter()) {
            if raw.is_ascii_lowercase() {
                *encoded = crate::filter::K_NUCL_MASK;
            }
        }
    }
}

/// blast-rs: Computes blastn API Karlin blocks from public parameters; not a
/// direct NCBI C port.
fn blastn_api_kbps(
    query_plus: &[u8],
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
) -> (KarlinBlk, KarlinBlk) {
    const BLASTNA_SIZE: usize = 16;

    let matrix = blast_score_blk_nucl_matrix_create(reward, penalty);
    let matrix_fn = |i: usize, j: usize| -> i32 { matrix[i][j] };
    let mut lo = i32::MAX;
    let mut hi = i32::MIN;
    for row in matrix.iter().take(BLASTNA_SIZE) {
        for &score in row.iter().take(BLASTNA_SIZE) {
            if score <= -100000000 || score >= 100000000 {
                continue;
            }
            lo = lo.min(score);
            hi = hi.max(score);
        }
    }

    let ctx = UngappedKbpContext {
        query_offset: 0,
        query_length: query_plus.len() as i32,
        is_valid: true,
    };
    let ambiguous_residues: &[u8] = &[14, 15];
    let ungapped = ungapped_kbp_calc(
        query_plus,
        &[ctx],
        lo,
        hi,
        BLASTNA_SIZE,
        ambiguous_residues,
        &matrix_fn,
    )[0]
    .clone()
    .unwrap_or(KarlinBlk {
        lambda: 1.374,
        k: 0.621,
        log_k: 0.621_f64.ln(),
        h: 1.286,
        round_down: false,
    });
    let mut gapped = KarlinBlk::default();
    let mut round_down = false;
    let mut matrix_error = None;
    if blast_karlin_blk_nucl_gapped_calc(
        Some(&mut gapped),
        gap_open,
        gap_extend,
        reward,
        penalty,
        Some(&ungapped),
        Some(&mut round_down),
        Some(&mut matrix_error),
    ) != 0
    {
        gapped = ungapped.clone();
    }
    (ungapped, gapped)
}

/// blast-rs: Computes blastn API effective search statistics; not a direct
/// NCBI C port.
fn blastn_api_stats(
    params: &SearchParams,
    query_plus: &[u8],
    total_subject_len: i64,
    num_subjects: i32,
) -> (KarlinBlk, f64, i32) {
    let (ungapped, gapped) = blastn_api_kbps(
        query_plus,
        params.match_score,
        params.mismatch,
        params.gap_open,
        params.gap_extend,
    );
    if params.effective_search_space > 0 {
        return (gapped, params.effective_search_space as f64, 0);
    }

    let database_length = if params.db_length > 0 {
        params.db_length
    } else {
        total_subject_len
    };
    let (mut alpha, mut beta) = (0.0, 0.0);
    let _ = blast_get_nucl_alpha_beta(
        params.match_score,
        params.mismatch,
        params.gap_open,
        params.gap_extend,
        ungapped.lambda,
        ungapped.h,
        true,
        &mut alpha,
        &mut beta,
    );
    let qlen = query_plus.len() as i32;
    let mut len_adj = 0;
    let _ = blast_compute_length_adjustment(
        gapped.k,
        gapped.log_k,
        alpha / gapped.lambda,
        beta,
        qlen,
        database_length,
        num_subjects,
        Some(&mut len_adj),
    );
    let eff_db = (database_length - num_subjects as i64 * len_adj as i64).max(1);
    let search_space = eff_db as f64 * (qlen - len_adj).max(1) as f64;
    (gapped, search_space, len_adj)
}

/// Run a blastx search (translated nucleotide query vs protein database).
/// blast-rs: Public blastx API pipeline assembled from ported lower-level pieces;
/// not a direct NCBI C port.
/// blast-rs: per-(frame, subject) blastx result body, extracted from `blastx` so
/// the cross-query batch driver can reuse it. Returns the HSPs for this frame x
/// subject (frame-oriented nucleotide query coords); the caller pushes them.
#[allow(clippy::too_many_arguments)]
fn process_blastx_frame_hsps(
    phits: Vec<crate::protein_lookup::ProteinHit>,
    query_data: &[u8],
    prot: &[u8],
    frame: i32,
    subj_aa: &[u8],
    subj_len: usize,
    query_info: &crate::queryinfo::QueryInfo,
    query_context: usize,
    search_space: f64,
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
    prot_kbp: &crate::stat::KarlinBlk,
    ungapped_kbp: &crate::stat::KarlinBlk,
    gumbel_blk: &Option<crate::stat::GumbelBlk>,
    x_drop_final: i32,
    hit_cutoff: i32,
    redo_cutoff_score_min: i32,
    link_word_cutoff_score_min: i32,
    translated_sum_stats: bool,
    params: &SearchParams,
    _scratch: &mut ProteinScratch,
) -> Vec<Hsp> {
    let mut out_hsps: Vec<Hsp> = Vec::new();
    let matrix_name = protein_matrix_name(params.matrix);
    if phits.is_empty() {
        return out_hsps;
    }
    if params.comp_adjust > 0 {
        let query_index = query_info
            .contexts
            .get(query_context)
            .map(|ctx| ctx.query_index)
            .unwrap_or(0);
        let mut this_match = blast_hsp_list_new(query_index);
        for hit in phits {
            blast_hsp_list_push(
                &mut this_match,
                protein_hit_to_blast_hsp(hit, query_context as i32, frame, 0),
            );
        }

        let mut matrix_info = crate::blast_kappa::BlastMatrixInfo::default();
        let ideal_lambda = crate::stat::protein_ideal_ungapped_kbp_for_matrix(matrix_name).lambda;
        if crate::blast_kappa::s_matrix_info_init(
            &mut matrix_info,
            matrix_name,
            ideal_lambda,
            COMPO_ADJUST_SCALE_FACTOR,
        ) != 0
        {
            return out_hsps;
        }

        let scaled_kbp = crate::stat::KarlinBlk {
            lambda: prot_kbp.lambda / COMPO_ADJUST_SCALE_FACTOR,
            k: prot_kbp.log_k.exp(),
            log_k: prot_kbp.log_k,
            h: prot_kbp.h,
            round_down: prot_kbp.round_down,
        };
        let score_block = crate::link_hsps::LinkScoreBlock {
            kbp: vec![scaled_kbp.clone(); query_info.contexts.len().max(1)],
            kbp_gap: vec![scaled_kbp; query_info.contexts.len().max(1)],
            gbp: gumbel_blk.clone(),
            link_gbp_db_length: gumbel_blk.as_ref().map(|gbp| gbp.db_length.max(1)),
            recompute_evalues_before_uneven_linking: true,
        };
        let link_params = crate::link_hsps::LinkHSPParameters {
            gap_prob: crate::stat::BLAST_GAP_PROB_GAPPED,
            gap_decay_rate: crate::stat::BLAST_GAP_DECAY_RATE_GAPPED,
            cutoff_small_gap: link_word_cutoff_score_min,
            longest_intron: translated_gapped_link_longest_intron(params.max_intron_length),
            ..crate::link_hsps::LinkHSPParameters::default()
        };
        let link_context = crate::blast_kappa::HitlistLinkContext {
            query_info,
            query_context,
            score_block: &score_block,
            link_params: &link_params,
            gapped_calculation: true,
        };
        let scaled_gap_open =
            crate::math::blast_nint(params.gap_open as f64 * COMPO_ADJUST_SCALE_FACTOR) as i32;
        let scaled_gap_extend =
            crate::math::blast_nint(params.gap_extend as f64 * COMPO_ADJUST_SCALE_FACTOR) as i32;
        let scaled_x_drop_final =
            crate::math::blast_nint(x_drop_final as f64 * COMPO_ADJUST_SCALE_FACTOR) as i32;
        let do_link_hsps = translated_sum_stats;
        let align_params = crate::blast_kappa::blast_redo_align_params_new(
            matrix_info,
            crate::blast_kappa::BlastCompoGappingParams {
                gap_open: scaled_gap_open,
                gap_extend: scaled_gap_extend,
                decline_align: i32::MIN,
                x_dropoff: scaled_x_drop_final,
                context: None,
            },
            crate::blast_kappa::CompoAdjustMode::from_u8(params.comp_adjust),
            COMPO_ADJUST_SCALE_FACTOR,
            false,
            true,
            false,
            query_info.max_length as i32,
            if do_link_hsps {
                composition_scaled_cutoff(redo_cutoff_score_min)
            } else {
                1
            },
            params.evalue_threshold,
            do_link_hsps,
            0.0,
        );
        let scoring_options = crate::options::ScoringOptions {
            matrix_path: None,
            reward: 0,
            penalty: 0,
            gap_open: params.gap_open,
            gap_extend: params.gap_extend,
            shift_pen: i16::MAX as i32,
            gapped_calculation: true,
            complexity_adjusted_scoring: false,
            matrix_name: Some(matrix_name.to_string()),
            is_ooframe: false,
            program_number: crate::program::UNDEFINED,
        };
        let mut scoring = crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let mut matrix_vec: Vec<Vec<i32>> = matrix.iter().map(|row| row.to_vec()).collect();
        let mut kbp_gap = vec![prot_kbp.clone(); query_info.contexts.len().max(1)];
        let mut saved = crate::blast_kappa::BlastKappaSavedParameters::default();
        let mut results =
            crate::hspstream::HspResults::new(query_info.contexts.len().max(1) as i32);
        let subject = crate::blast_kappa::BlastRedoInMemorySubject {
            subject_source: subj_aa,
            reward: 0,
            penalty: 0,
            genetic_code: crate::util::lookup_genetic_code(params.db_gencode),
            smith_waterman: false,
            expect_value: params.evalue_threshold,
            hitlist_size: params
                .max_hsps
                .unwrap_or(crate::options::HITLIST_SIZE as usize) as i32,
            inclusion_ethresh: f64::INFINITY,
            link_context: do_link_hsps.then_some(&link_context),
        };
        let mut legacy_this_match = crate::hspstream::HspList::new(this_match.oid);
        let status = crate::blast_kappa::blast_redo_alignment_core_mt(
            crate::program::BLASTX,
            1,
            query_data,
            query_info,
            &mut kbp_gap,
            &mut matrix_vec,
            &mut scoring,
            &align_params,
            &mut saved,
            &mut legacy_this_match,
            crate::blast_kappa::BlastRedoAlignmentSource::InMemorySubjectHspList {
                hsp_list: &mut this_match,
                subject,
            },
            &mut results,
        );
        if status != 0 {
            return out_hsps;
        }

        for item in this_match.hsp_array.into_iter().flatten() {
            if item.evalue > params.evalue_threshold {
                continue;
            }
            let q_start_aa = item.query.offset.max(0) as usize;
            let q_end_aa = item.query.end.max(item.query.offset).max(0) as usize;
            let s_start = item.subject.offset.max(0) as usize;
            let s_end = item.subject.end.max(item.subject.offset).max(0) as usize;
            if q_end_aa > prot.len() || s_end > subj_aa.len() {
                continue;
            }
            let (query_start, query_end) =
                crate::util::protein_to_oriented_nuc_coords(q_start_aa, q_end_aa, frame);
            let (query_gapped_start, _) = crate::util::protein_to_oriented_nuc_coords(
                item.query.gapped_start.max(0) as usize,
                item.query.gapped_start.max(0) as usize + 1,
                frame,
            );
            let (q_aln, s_aln, align_length, num_gaps) =
                if let Some(script) = item.gap_info.as_ref() {
                    let q_slice = &prot[q_start_aa..q_end_aa];
                    let s_slice = &subj_aa[s_start..s_end];
                    let (qs, ss) =
                        script.render_alignment(q_slice, s_slice, ncbistdaa_to_aminoacid_char);
                    (
                        qs.into_bytes(),
                        ss.into_bytes(),
                        script.alignment_length().max(0) as usize,
                        crate::blast_kappa::gap_edit_script_num_gap_opens(script).max(0) as usize,
                    )
                } else {
                    (
                        Vec::new(),
                        Vec::new(),
                        q_end_aa.saturating_sub(q_start_aa),
                        0,
                    )
                };
            out_hsps.push(Hsp {
                score: item.score,
                bit_score: prot_kbp.raw_to_bit(item.score),
                evalue: item.evalue,
                query_start,
                query_end,
                subject_start: s_start,
                subject_end: s_end,
                query_gapped_start,
                subject_gapped_start: item.subject.gapped_start.max(0) as usize,
                query_link_start: q_start_aa,
                query_link_end: q_end_aa,
                query_link_gapped_start: item.query.gapped_start.max(0) as usize,
                subject_link_start: s_start,
                subject_link_end: s_end,
                subject_link_gapped_start: item.subject.gapped_start.max(0) as usize,
                link_score: None,
                link_lambda: None,
                num_identities: item.num_ident.max(0) as usize,
                num_gaps,
                alignment_length: align_length,
                query_aln: q_aln,
                midline: Vec::new(),
                subject_aln: s_aln,
                query_frame: frame,
                subject_frame: 0,
                num_links: item.num.max(1),
                comp_adjust_method: item.comp_adjustment_method as u8,
            });
        }
        return out_hsps;
    }
    let (final_phits, use_adj_matrix, lambda_ratio_opt, comp_adjust_method_id) =
        (phits, false, None::<f64>, 0);
    for mut ph in final_phits {
        if params.comp_adjust > 0 {
            if let Some(scaled_score) = ph.scaled_score {
                if scaled_score < composition_scaled_cutoff(hit_cutoff) {
                    continue;
                }
            } else if ph.score < hit_cutoff {
                continue;
            }
        }
        // H2: for `-ungapped`, NCBI frees `sbp->gbp` (gapped_calculation=FALSE)
        // and reports bit scores / e-values from the PER-FRAME ungapped Karlin
        // block via simple Karlin (not the gapped kbp, not Spouge). The gapped
        // path keeps the gapped `prot_kbp`.
        let bit_score = if params.ungapped {
            ungapped_kbp.raw_to_bit(ph.score)
        } else {
            prot_kbp.raw_to_bit(ph.score)
        };
        let (e_score_i32, e_lambda) = match ph.scaled_score {
            Some(s) if params.comp_adjust > 0 => (s, prot_kbp.lambda / COMPO_ADJUST_SCALE_FACTOR),
            _ => (ph.score, prot_kbp.lambda),
        };
        let e_kbp = crate::stat::KarlinBlk {
            lambda: e_lambda,
            k: prot_kbp.k,
            log_k: prot_kbp.log_k,
            h: prot_kbp.h,
            round_down: prot_kbp.round_down,
        };
        let evalue = if let Some(adjusted_evalue) = ph.adjusted_evalue {
            adjusted_evalue
        } else if params.ungapped {
            // Ungapped: simple Karlin with the per-frame ungapped kbp.
            ungapped_kbp.raw_to_evalue(ph.score, search_space)
        } else if let Some(ref gbp) = gumbel_blk {
            let mut base_ev = crate::stat::blast_spouge_sto_e(
                e_score_i32,
                Some(&e_kbp),
                Some(gbp),
                prot.len() as i32,
                subj_len as i32,
            );
            if translated_sum_stats {
                base_ev /= crate::stat::blast_gap_decay_divisor(
                    crate::stat::BLAST_GAP_DECAY_RATE_GAPPED,
                    1,
                );
            }
            if use_adj_matrix {
                base_ev
            } else if let Some(lr) = lambda_ratio_opt {
                let scaled_kbp = crate::stat::KarlinBlk {
                    lambda: e_kbp.lambda / lr,
                    k: e_kbp.k,
                    log_k: e_kbp.log_k,
                    h: e_kbp.h,
                    round_down: e_kbp.round_down,
                };
                let mut evalue = crate::stat::blast_spouge_sto_e(
                    e_score_i32,
                    Some(&scaled_kbp),
                    Some(gbp),
                    prot.len() as i32,
                    subj_len as i32,
                );
                if translated_sum_stats {
                    evalue /= crate::stat::blast_gap_decay_divisor(
                        crate::stat::BLAST_GAP_DECAY_RATE_GAPPED,
                        1,
                    );
                }
                evalue
            } else {
                base_ev
            }
        } else {
            let raw_evalue = e_kbp.raw_to_evalue(e_score_i32, search_space);
            if use_adj_matrix {
                raw_evalue
            } else if let Some(lr) = lambda_ratio_opt {
                let scaled_lambda = e_kbp.lambda / lr;
                let scaled_kbp = crate::stat::KarlinBlk {
                    lambda: scaled_lambda,
                    k: e_kbp.k,
                    log_k: e_kbp.log_k,
                    h: e_kbp.h,
                    round_down: e_kbp.round_down,
                };
                scaled_kbp.raw_to_evalue(e_score_i32, search_space)
            } else {
                raw_evalue
            }
        };
        let (query_start, query_end) =
            crate::util::protein_to_oriented_nuc_coords(ph.query_start, ph.query_end, frame);
        let (query_gapped_start, _) = crate::util::protein_to_oriented_nuc_coords(
            ph.gapped_start_q,
            ph.gapped_start_q.saturating_add(1),
            frame,
        );
        let (q_aln, s_aln) = match (ph.qseq.take(), ph.sseq.take()) {
            (Some(qs), Some(ss)) => (qs.into_bytes(), ss.into_bytes()),
            _ => (Vec::new(), Vec::new()),
        };
        out_hsps.push(Hsp {
            score: ph.score,
            bit_score,
            evalue,
            query_start,
            query_end,
            subject_start: ph.subject_start,
            subject_end: ph.subject_end,
            query_gapped_start,
            subject_gapped_start: ph.gapped_start_s,
            query_link_start: ph.query_start,
            query_link_end: ph.query_end,
            query_link_gapped_start: ph.gapped_start_q,
            subject_link_start: ph.subject_start,
            subject_link_end: ph.subject_end,
            subject_link_gapped_start: ph.gapped_start_s,
            link_score: ph.scaled_score,
            link_lambda: ph
                .scaled_score
                .map(|_| prot_kbp.lambda / COMPO_ADJUST_SCALE_FACTOR),
            num_identities: ph.num_ident as usize,
            num_gaps: ph.gap_opens as usize,
            alignment_length: ph.align_length as usize,
            query_aln: q_aln,
            midline: Vec::new(),
            subject_aln: s_aln,
            query_frame: frame,
            subject_frame: 0,
            num_links: 1,
            comp_adjust_method: comp_adjust_method_id,
        });
    }
    out_hsps
}

pub fn blastx(db: &BlastDb, query: &[u8], params: &SearchParams) -> Vec<SearchResult> {
    if query.len() < 3 {
        return Vec::new();
    }
    let query_ncbi4na = encode_ncbi4na_sequence(query);
    let code = crate::util::lookup_genetic_code(params.query_gencode);
    // 1-1 with blast_engine.c:s_BlastSearchEngineCore: call
    // BLAST_GetAllTranslations once, then iterate the 6 contexts via the
    // returned (translation_buffer, frame_offsets) pair.
    let (mut translation_buffer, frame_offsets) =
        crate::util::blast_get_all_translations(&query_ncbi4na, query_ncbi4na.len(), code);
    if params.filter_low_complexity {
        for ctx in 0..crate::util::NUM_FRAMES {
            let begin = (frame_offsets[ctx] + 1) as usize;
            let end = frame_offsets[ctx + 1] as usize;
            if begin < end {
                apply_seg_ncbistdaa_with_options(
                    &mut translation_buffer[begin..end],
                    params.seg_window,
                    params.seg_locut,
                    params.seg_hicut,
                );
            }
        }
    }
    let matrix = *get_matrix(params.matrix);
    let matrix_name = protein_matrix_name(params.matrix);
    let prot_kbp = protein_kbp_for_matrix(params.matrix, params.gap_open, params.gap_extend);

    let ungapped_kbp = crate::stat::protein_ungapped_kbp_for_matrix(matrix_name);
    let x_drop_gapped = raw_gapped_xdrop_bits(params.x_drop_gapped, &prot_kbp);
    let x_drop_final = raw_gapped_xdrop_bits(params.x_drop_final, &prot_kbp);
    // H2: gap_trigger is now computed per-frame from the per-context ungapped
    // Karlin block (see `ungapped_kbp_array` below); no single rounded value.

    let total_subj_len: usize = (0..db.num_oids)
        .map(|oid| db.get_seq_len(oid) as usize)
        .sum();
    let stats_db_len = statistical_db_length(params, total_subj_len as i64);
    let avg_subject_length = (total_subj_len / db.num_oids.max(1) as usize).max(1);
    // Zero out Gumbel when -searchsp overrides the search space — see the
    // blastp branch for the rationale (NCBI uses Karlin in that case).
    let gumbel_blk = if params.effective_search_space > 0 {
        None
    } else {
        protein_gumbel_for_matrix(
            params.matrix,
            params.gap_open,
            params.gap_extend,
            stats_db_len,
        )
    };
    let translated_sum_stats = translated_sum_stats_enabled(params);
    let prelim_evalue = composition_prelim_evalue(params);

    let word_size = params.word_size.clamp(2, 7);
    let threshold = params.word_threshold.unwrap_or_else(|| {
        suggested_word_threshold(params.matrix, crate::program::BLASTX, params.word_size)
    });
    let max_hsps = params.max_hsps;

    // Build a 6-context QueryInfo from the translation offsets and call
    // BLAST_CalcEffLengths once (mirrors NCBI's setup-time invocation in
    // BLAST_GapAlignSetUp). The per-context `eff_searchsp` /
    // `length_adjustment` are then read inside the loop below.
    let mut query_info =
        crate::queryinfo::QueryInfo::new_translated_query_from_offsets(&frame_offsets);
    let scoring_options = crate::options::ScoringOptions {
        matrix_path: None,
        reward: 0,
        penalty: 0,
        gap_open: params.gap_open,
        gap_extend: params.gap_extend,
        shift_pen: i16::MAX as i32,
        // NCBI propagates `-ungapped` as `gapped_calculation = FALSE` into
        // BLAST_CalcEffLengths, which uses ungapped alpha/beta for length
        // adjustment (`blast_stat.c:3128`).
        gapped_calculation: !params.ungapped,
        complexity_adjusted_scoring: false,
        matrix_name: Some(matrix_name.to_string()),
        is_ooframe: false,
        program_number: crate::program::UNDEFINED,
    };
    let eff_params = crate::parameters::EffectiveLengthsParameters {
        options: effective_lengths_options(params),
        real_db_length: stats_db_len,
        real_num_seqs: db.num_oids as i32,
    };
    // H2: NCBI `Blast_ScoreBlkKbpUngappedCalc` (blast_stat.c:2737) computes the
    // ungapped Karlin block PER CONTEXT from the translated query-frame's amino
    // acid composition, applying the ideal-Lambda cap for translated queries
    // (`check_ideal` is TRUE for blastx). These per-frame ungapped kbp feed the
    // effective-length std array, the per-frame gap_trigger, and (in -ungapped
    // mode) the ungapped-HSP bit scores and e-values. Build that array here.
    let ideal_ungapped_kbp = crate::stat::protein_ideal_ungapped_kbp_for_matrix(matrix_name);
    let mut ungapped_kbp_array = vec![ungapped_kbp.clone(); crate::util::NUM_FRAMES];
    for ctx in 0..crate::util::NUM_FRAMES {
        let begin = (frame_offsets[ctx] + 1) as usize;
        let end = frame_offsets[ctx + 1] as usize;
        if begin >= end {
            continue;
        }
        let prot: &[u8] = &translation_buffer[begin..end];
        if crate::composition::blast_read_aa_composition(prot, AA_SIZE).1 == 0 {
            continue;
        }
        let mut ctx_kbp =
            crate::stat::query_specific_protein_ungapped_kbp_for_matrix(prot, matrix_name, &matrix);
        // Translated-query ideal-Lambda cap (`blast_stat.c:2796`).
        if ctx_kbp.lambda >= ideal_ungapped_kbp.lambda {
            ctx_kbp = ideal_ungapped_kbp.clone();
        }
        ungapped_kbp_array[ctx] = ctx_kbp;
    }

    let kbp_array = vec![prot_kbp.clone(); crate::util::NUM_FRAMES];
    let kbp_std_array = ungapped_kbp_array.clone();
    crate::blast_setup::blast_calc_eff_lengths(
        crate::program::BLASTX,
        &scoring_options,
        &eff_params,
        &kbp_array,
        &kbp_std_array,
        matrix_name,
        &mut query_info,
    );
    let link_avg_query_length = query_info_average_context_length(&query_info);

    // Concatenate all six frames into ONE lookup table and scan each subject
    // once, replacing the per-frame DB rescan. Frames are separated by
    // PROTEIN_CONCAT_SENTINEL (>= AA_SIZE, so boundary-spanning words are
    // skipped and ungapped extension stops at the boundary); each hit is
    // attributed back to its frame via bsearch_context_info, then the unchanged
    // per-frame gapped + coordinate pipeline runs. blastx uses a single
    // x_drop_ungapped, so the shared scan is byte-identical to the per-frame one.
    let mut frame_plans: Vec<(usize, usize, i32, f64, usize)> = Vec::new();
    for ctx in 0..crate::util::NUM_FRAMES {
        let frame = crate::util::blast_context_to_frame(ctx as u32);
        let begin = (frame_offsets[ctx] + 1) as usize;
        let end = frame_offsets[ctx + 1] as usize;
        if begin >= end {
            continue;
        }
        let prot: &[u8] = &translation_buffer[begin..end];
        if prot.len() < word_size {
            continue;
        }
        if crate::composition::blast_read_aa_composition(prot, AA_SIZE).1 == 0 {
            continue;
        }
        let search_space = query_info.contexts[ctx].eff_searchsp.max(1) as f64;
        frame_plans.push((begin, end, frame, search_space, ctx));
    }
    let frame_lengths: Vec<usize> = frame_plans.iter().map(|&(b, e, _, _, _)| e - b).collect();
    let attrib_info = crate::queryinfo::QueryInfo::new_blastp(&frame_lengths);
    let mut concat: Vec<u8> =
        vec![crate::protein_lookup::PROTEIN_CONCAT_SENTINEL; attrib_info.seq_buf_len().max(1)];
    for (fi, &(b, e, _, _, _)) in frame_plans.iter().enumerate() {
        let off = attrib_info.contexts[fi].query_offset as usize;
        concat[off..off + (e - b)].copy_from_slice(&translation_buffer[b..e]);
    }
    let concat_lookup = if frame_plans.is_empty() {
        None
    } else {
        Some(crate::protein_lookup::ProteinLookupTable::build(
            &concat, word_size, &matrix, threshold,
        ))
    };

    // NCBI's two-hit finder reads x-drop and save cutoff from the active query
    // context (`aa_ungapped.c:560,588`), so the concatenated-frame lookup must
    // carry those per-context word parameters into the scanner.
    let context_scan_params: Vec<crate::protein_lookup::ProteinContextScanParams> = frame_plans
        .iter()
        .enumerate()
        .map(|(fi, &(begin, end, _frame, search_space, ctx))| {
            let frame_ungapped_kbp = &ungapped_kbp_array[ctx];
            let frame_gap_trigger_raw =
                raw_gap_trigger_bits(crate::stat::BLAST_GAP_TRIGGER_PROT, frame_ungapped_kbp);
            let cutoff_score = protein_prelim_seed_cutoff(
                frame_gap_trigger_raw,
                prelim_evalue,
                &prot_kbp,
                gumbel_blk.as_ref(),
                (end - begin) as i32,
                avg_subject_length as i32,
                search_space,
            );
            crate::protein_lookup::ProteinContextScanParams {
                query_start: attrib_info.contexts[fi].query_offset as usize,
                x_dropoff: raw_ungapped_xdrop_bits(params.x_drop_ungapped, frame_ungapped_kbp),
                cutoff_score: cutoff_score.max(1),
            }
        })
        .collect();

    let mut results = vec![None; db.num_oids as usize];
    let mut scratch = ProteinScratch::new();
    if let Some(concat_lookup) = concat_lookup.as_ref() {
        for oid in 0..db.num_oids {
            let subj_raw = db.get_sequence(oid);
            let subj_len = db.get_seq_len(oid) as usize;
            if subj_len == 0 {
                continue;
            }
            let subj_aa = &subj_raw[..subj_len];
            // Defer defline parsing until this subject produces a reported HSP.
            let mut accession: Option<String> = None;
            let mut title: Option<String> = None;
            let mut subject_internal_hsps = blast_hsp_list_new(0);
            let max_query_prot_len = frame_plans
                .iter()
                .map(|&(begin, end, _, _, _)| end.saturating_sub(begin))
                .max()
                .unwrap_or(0)
                .max(1);
            if !params.ungapped {
                scratch
                    .prelim_tree
                    .reset_for_query(max_query_prot_len as i32 + 1, subj_aa.len() as i32 + 1);
            }
            scratch.diag_buf.clear();
            let scan_hits = crate::protein_lookup::protein_scan_with_table_reuse_context_params(
                &concat,
                subj_aa,
                &matrix,
                concat_lookup,
                params.two_hit_window as i32,
                &mut scratch.diag_buf,
                &context_scan_params,
            );
            if scan_hits.is_empty() {
                continue;
            }
            let mut per_frame: Vec<Vec<crate::protein_lookup::ProteinHit>> =
                (0..frame_plans.len()).map(|_| Vec::new()).collect();
            for mut h in scan_hits {
                let fi = crate::queryinfo::bsearch_context_info(h.query_start as i32, &attrib_info);
                if fi < 0 || fi as usize >= frame_plans.len() {
                    continue;
                }
                let fi = fi as usize;
                let off = attrib_info.contexts[fi].query_offset as usize;
                h.query_start -= off;
                h.query_end = h.query_end.saturating_sub(off);
                h.gapped_start_q = h.gapped_start_q.saturating_sub(off);
                per_frame[fi].push(h);
            }
            for (fi, &(begin, end, frame, search_space, ctx)) in frame_plans.iter().enumerate() {
                let frame_ungapped = std::mem::take(&mut per_frame[fi]);
                if frame_ungapped.is_empty() {
                    continue;
                }
                let prot: &[u8] = &translation_buffer[begin..end];
                // H2: per-frame ungapped Karlin block for gap_trigger and
                // ungapped-HSP statistics.
                let frame_ungapped_kbp = &ungapped_kbp_array[ctx];
                let frame_gap_trigger_raw =
                    raw_gap_trigger_bits(crate::stat::BLAST_GAP_TRIGGER_PROT, frame_ungapped_kbp);
                let blastx_seed_cutoff = protein_prelim_seed_cutoff(
                    frame_gap_trigger_raw,
                    prelim_evalue,
                    &prot_kbp,
                    gumbel_blk.as_ref(),
                    prot.len() as i32,
                    avg_subject_length as i32,
                    search_space,
                );
                let blastx_hit_cutoff = protein_eval_cutoff(
                    prelim_evalue,
                    &prot_kbp,
                    gumbel_blk.as_ref(),
                    prot.len() as i32,
                    avg_subject_length as i32,
                    search_space,
                );
                let blastx_redo_cutoff = blastx_hit_cutoff;
                let blastx_hit_cutoff = if translated_sum_stats && !params.ungapped {
                    protein_sum_stats_hit_cutoff(
                        blastx_hit_cutoff,
                        &prot_kbp,
                        link_avg_query_length,
                        avg_subject_length as i32,
                    )
                } else {
                    blastx_hit_cutoff
                };
                if params.ungapped {
                    let mut ug = frame_ungapped;
                    ug.iter_mut().filter(|uh| uh.score >= 1).for_each(|uh| {
                        uh.qseq = Some(
                            prot[uh.query_start..uh.query_end]
                                .iter()
                                .map(|&aa| ncbistdaa_to_aminoacid_char(aa))
                                .collect(),
                        );
                        uh.sseq = Some(
                            subj_aa[uh.subject_start..uh.subject_end]
                                .iter()
                                .map(|&aa| ncbistdaa_to_aminoacid_char(aa))
                                .collect(),
                        );
                    });
                    for hsp in process_blastx_frame_hsps(
                        ug,
                        &translation_buffer,
                        prot,
                        frame,
                        subj_aa,
                        subj_len,
                        &query_info,
                        ctx,
                        search_space,
                        &matrix,
                        &prot_kbp,
                        frame_ungapped_kbp,
                        &gumbel_blk,
                        x_drop_final,
                        blastx_hit_cutoff,
                        blastx_redo_cutoff,
                        blastx_seed_cutoff,
                        translated_sum_stats,
                        params,
                        &mut scratch,
                    ) {
                        if accession.is_none() {
                            accession = Some(
                                db.get_accession(oid)
                                    .unwrap_or_else(|| format!("oid_{}", oid)),
                            );
                            title = Some(String::from_utf8_lossy(db.get_header(oid)).to_string());
                        }
                        push_hsp_for_subject(
                            &mut results,
                            oid,
                            title.as_deref().unwrap(),
                            accession.as_deref().unwrap(),
                            subj_len,
                            &[],
                            hsp,
                        );
                    }
                } else {
                    let phits = protein_prelim_gapped_hits(
                        prot,
                        subj_aa,
                        &matrix,
                        &frame_ungapped,
                        params.gap_open,
                        params.gap_extend,
                        x_drop_gapped,
                        x_drop_final,
                        blastx_seed_cutoff.max(1),
                        blastx_hit_cutoff.max(1),
                        params.chaining && params.matrix == MatrixType::Blosum62,
                        &mut scratch,
                        false,
                        ctx as i32,
                        0,
                    );
                    for hit in phits {
                        blast_hsp_list_push(
                            &mut subject_internal_hsps,
                            protein_hit_to_blast_hsp(hit, ctx as i32, frame, 0),
                        );
                    }
                }
            }
            if !params.ungapped && subject_internal_hsps.hspcnt > 0 {
                blastx_traceback_subject_hsp_list(
                    &translation_buffer,
                    &frame_offsets,
                    subj_aa,
                    &matrix,
                    &mut subject_internal_hsps,
                    params.gap_open,
                    params.gap_extend,
                    x_drop_final,
                    &mut scratch,
                    max_query_prot_len,
                );
                subject_internal_hsps.hsp_array.sort_by(|a, b| {
                    let Some(a) = a.as_ref() else {
                        return std::cmp::Ordering::Greater;
                    };
                    let Some(b) = b.as_ref() else {
                        return std::cmp::Ordering::Less;
                    };
                    a.context
                        .cmp(&b.context)
                        .then((a.query.frame as i32).cmp(&(b.query.frame as i32)))
                        .then_with(|| compare_blast_hsps_by_score(a, b))
                });
                let mut group_start = 0usize;
                while group_start < subject_internal_hsps.hsp_array.len() {
                    let Some(first) = subject_internal_hsps.hsp_array[group_start].as_ref() else {
                        group_start += 1;
                        continue;
                    };
                    let context = first.context.max(0) as usize;
                    if context + 1 >= frame_offsets.len() {
                        group_start += 1;
                        continue;
                    }
                    let frame = first.query.frame as i32;
                    let query_range = (frame_offsets[context] + 1).max(0) as usize
                        ..frame_offsets[context + 1].max(0) as usize;
                    let search_space = query_info.contexts[context].eff_searchsp.max(1) as f64;
                    let frame_ungapped_kbp = &ungapped_kbp_array[context];
                    let frame_gap_trigger_raw = raw_gap_trigger_bits(
                        crate::stat::BLAST_GAP_TRIGGER_PROT,
                        frame_ungapped_kbp,
                    );
                    let seed_cutoff = protein_prelim_seed_cutoff(
                        frame_gap_trigger_raw,
                        prelim_evalue,
                        &prot_kbp,
                        gumbel_blk.as_ref(),
                        query_range.len() as i32,
                        avg_subject_length as i32,
                        search_space,
                    );
                    let mut hit_cutoff = protein_eval_cutoff(
                        prelim_evalue,
                        &prot_kbp,
                        gumbel_blk.as_ref(),
                        query_range.len() as i32,
                        avg_subject_length as i32,
                        search_space,
                    );
                    if translated_sum_stats {
                        hit_cutoff = protein_sum_stats_hit_cutoff(
                            hit_cutoff,
                            &prot_kbp,
                            link_avg_query_length,
                            avg_subject_length as i32,
                        );
                    }
                    let mut group_end = group_start + 1;
                    while group_end < subject_internal_hsps.hsp_array.len()
                        && subject_internal_hsps.hsp_array[group_end]
                            .as_ref()
                            .is_some_and(|hsp| {
                                hsp.context.max(0) as usize == context
                                    && hsp.query.frame as i32 == frame
                            })
                    {
                        group_end += 1;
                    }
                    let prot = &translation_buffer[query_range.clone()];
                    let phits: Vec<_> = subject_internal_hsps.hsp_array[group_start..group_end]
                        .iter()
                        .filter_map(|item| item.as_ref().map(blast_hsp_to_protein_hit))
                        .collect();
                    for hsp in process_blastx_frame_hsps(
                        phits,
                        &translation_buffer,
                        prot,
                        frame,
                        subj_aa,
                        subj_len,
                        &query_info,
                        context,
                        search_space,
                        &matrix,
                        &prot_kbp,
                        frame_ungapped_kbp,
                        &gumbel_blk,
                        x_drop_final,
                        hit_cutoff,
                        hit_cutoff,
                        seed_cutoff,
                        translated_sum_stats,
                        params,
                        &mut scratch,
                    ) {
                        if accession.is_none() {
                            accession = Some(
                                db.get_accession(oid)
                                    .unwrap_or_else(|| format!("oid_{}", oid)),
                            );
                            title = Some(String::from_utf8_lossy(db.get_header(oid)).to_string());
                        }
                        push_hsp_for_subject(
                            &mut results,
                            oid,
                            title.as_deref().unwrap(),
                            accession.as_deref().unwrap(),
                            subj_len,
                            &[],
                            hsp,
                        );
                    }
                    group_start = group_end;
                }
            }
        }
    }

    let mut results: Vec<SearchResult> = results.into_iter().flatten().collect();
    if translated_sum_stats && !(params.comp_adjust == 0 && !params.ungapped) {
        apply_blastx_linked_sum_stats(
            &mut results,
            &query_info,
            &prot_kbp,
            gumbel_blk.as_ref(),
            stats_db_len,
            params.comp_adjust == 0,
            params.max_intron_length,
        );
    }
    reap_hsps_by_prelim_evalue(&mut results, prelim_evalue);
    apply_api_evalue_filter(&mut results, params.evalue_threshold);
    for result in &mut results {
        result.hsps.sort_by(compare_hsps_by_evalue_then_score);
    }
    apply_api_min_score_filter(&mut results, params.min_score);
    if let Some(culling_limit) = params.culling_limit {
        apply_api_culling_limit(&mut results, culling_limit, crate::program::BLASTX);
    }
    apply_api_max_hsps_limit(&mut results, max_hsps);
    results.sort_by(compare_search_results);
    if results.len() > params.max_target_seqs {
        results.truncate(params.max_target_seqs);
    }
    fill_protein_midlines(&mut results, ProteinMidlineStyle::AnyMismatchPlus);
    results
}

/// blast-rs: cross-query blastx. Concatenates the six frames of ALL queries into
/// ONE lookup table and scans the database once (NCBI's concatenated-query engine
/// model), instead of translating + re-scanning per query. Byte-identical to
/// looping `blastx` per query; reuses `process_blastx_frame_hsps`. Returns one
/// `Vec<SearchResult>` per query.
#[allow(clippy::too_many_arguments)]
pub fn blastx_batch(
    db: &BlastDb,
    queries: &[&[u8]],
    params: &SearchParams,
) -> Vec<Vec<SearchResult>> {
    if queries.is_empty() {
        return Vec::new();
    }
    let matrix = *get_matrix(params.matrix);
    let matrix_name = protein_matrix_name(params.matrix);
    let prot_kbp = protein_kbp_for_matrix(params.matrix, params.gap_open, params.gap_extend);
    let ungapped_kbp = crate::stat::protein_ungapped_kbp_for_matrix(matrix_name);
    let x_drop_gapped = raw_gapped_xdrop_bits(params.x_drop_gapped, &prot_kbp);
    let x_drop_final = raw_gapped_xdrop_bits(params.x_drop_final, &prot_kbp);
    let total_subj_len: usize = (0..db.num_oids)
        .map(|oid| db.get_seq_len(oid) as usize)
        .sum();
    let stats_db_len = statistical_db_length(params, total_subj_len as i64);
    let avg_subject_length = (total_subj_len / db.num_oids.max(1) as usize).max(1);
    let gumbel_blk = if params.effective_search_space > 0 {
        None
    } else {
        protein_gumbel_for_matrix(
            params.matrix,
            params.gap_open,
            params.gap_extend,
            stats_db_len,
        )
    };
    let translated_sum_stats = translated_sum_stats_enabled(params);
    let prelim_evalue = composition_prelim_evalue(params);
    let word_size = params.word_size.clamp(2, 7);
    let threshold = params.word_threshold.unwrap_or_else(|| {
        suggested_word_threshold(params.matrix, crate::program::BLASTX, params.word_size)
    });
    let max_hsps = params.max_hsps;
    let code = crate::util::lookup_genetic_code(params.query_gencode);

    let scoring_options = crate::options::ScoringOptions {
        matrix_path: None,
        reward: 0,
        penalty: 0,
        gap_open: params.gap_open,
        gap_extend: params.gap_extend,
        shift_pen: i16::MAX as i32,
        gapped_calculation: !params.ungapped,
        complexity_adjusted_scoring: false,
        matrix_name: Some(matrix_name.to_string()),
        is_ooframe: false,
        program_number: crate::program::UNDEFINED,
    };
    let eff_params = crate::parameters::EffectiveLengthsParameters {
        options: effective_lengths_options(params),
        real_db_length: stats_db_len,
        real_num_seqs: db.num_oids as i32,
    };
    let kbp_array = vec![prot_kbp.clone(); crate::util::NUM_FRAMES];

    struct QPlan {
        translation: Vec<u8>,
        frame_offsets: Vec<u32>,
        query_info: crate::queryinfo::QueryInfo,
        frames: Vec<(usize, usize, i32, f64, usize)>,
        ungapped_kbp_array: Vec<KarlinBlk>,
        link_avg_query_length: i32,
    }
    let qplans: Vec<QPlan> = queries
        .iter()
        .map(|q| {
            if q.len() < 3 {
                return QPlan {
                    translation: Vec::new(),
                    frame_offsets: vec![0; crate::util::NUM_FRAMES + 1],
                    query_info: crate::queryinfo::QueryInfo::new_blastp(&[]),
                    frames: Vec::new(),
                    ungapped_kbp_array: vec![ungapped_kbp.clone(); crate::util::NUM_FRAMES],
                    link_avg_query_length: 1,
                };
            }
            let query_ncbi4na = encode_ncbi4na_sequence(q);
            let (mut translation, frame_offsets) =
                crate::util::blast_get_all_translations(&query_ncbi4na, query_ncbi4na.len(), code);
            if params.filter_low_complexity {
                for ctx in 0..crate::util::NUM_FRAMES {
                    let begin = (frame_offsets[ctx] + 1) as usize;
                    let end = frame_offsets[ctx + 1] as usize;
                    if begin < end {
                        apply_seg_ncbistdaa_with_options(
                            &mut translation[begin..end],
                            params.seg_window,
                            params.seg_locut,
                            params.seg_hicut,
                        );
                    }
                }
            }
            let mut query_info =
                crate::queryinfo::QueryInfo::new_translated_query_from_offsets(&frame_offsets);
            let ideal_ungapped_kbp =
                crate::stat::protein_ideal_ungapped_kbp_for_matrix(matrix_name);
            let mut ungapped_kbp_array = vec![ungapped_kbp.clone(); crate::util::NUM_FRAMES];
            for ctx in 0..crate::util::NUM_FRAMES {
                let begin = (frame_offsets[ctx] + 1) as usize;
                let end = frame_offsets[ctx + 1] as usize;
                if begin >= end {
                    continue;
                }
                let prot: &[u8] = &translation[begin..end];
                if crate::composition::blast_read_aa_composition(prot, AA_SIZE).1 == 0 {
                    continue;
                }
                let mut ctx_kbp = crate::stat::query_specific_protein_ungapped_kbp_for_matrix(
                    prot,
                    matrix_name,
                    &matrix,
                );
                if ctx_kbp.lambda >= ideal_ungapped_kbp.lambda {
                    ctx_kbp = ideal_ungapped_kbp.clone();
                }
                ungapped_kbp_array[ctx] = ctx_kbp;
            }
            let kbp_std_array = ungapped_kbp_array.clone();
            crate::blast_setup::blast_calc_eff_lengths(
                crate::program::BLASTX,
                &scoring_options,
                &eff_params,
                &kbp_array,
                &kbp_std_array,
                matrix_name,
                &mut query_info,
            );
            let mut frames = Vec::new();
            for ctx in 0..crate::util::NUM_FRAMES {
                let frame = crate::util::blast_context_to_frame(ctx as u32);
                let begin = (frame_offsets[ctx] + 1) as usize;
                let end = frame_offsets[ctx + 1] as usize;
                if begin >= end {
                    continue;
                }
                let prot: &[u8] = &translation[begin..end];
                if prot.len() < word_size {
                    continue;
                }
                if crate::composition::blast_read_aa_composition(prot, AA_SIZE).1 == 0 {
                    continue;
                }
                let search_space = query_info.contexts[ctx].eff_searchsp.max(1) as f64;
                frames.push((begin, end, frame, search_space, ctx));
            }
            QPlan {
                translation,
                frame_offsets: frame_offsets.to_vec(),
                link_avg_query_length: query_info_average_context_length(&query_info),
                query_info,
                frames,
                ungapped_kbp_array,
            }
        })
        .collect();

    let mut flat: Vec<(usize, usize)> = Vec::new();
    for (qi, qp) in qplans.iter().enumerate() {
        for fj in 0..qp.frames.len() {
            flat.push((qi, fj));
        }
    }
    let flat_lengths: Vec<usize> = flat
        .iter()
        .map(|&(qi, fj)| {
            let (b, e, _, _, _) = qplans[qi].frames[fj];
            e - b
        })
        .collect();
    let attrib_info = crate::queryinfo::QueryInfo::new_blastp(&flat_lengths);
    let mut concat: Vec<u8> =
        vec![crate::protein_lookup::PROTEIN_CONCAT_SENTINEL; attrib_info.seq_buf_len().max(1)];
    for (gi, &(qi, fj)) in flat.iter().enumerate() {
        let (b, e, _, _, _) = qplans[qi].frames[fj];
        let off = attrib_info.contexts[gi].query_offset as usize;
        concat[off..off + (e - b)].copy_from_slice(&qplans[qi].translation[b..e]);
    }
    let lookup = if flat.is_empty() {
        None
    } else {
        Some(crate::protein_lookup::ProteinLookupTable::build(
            &concat, word_size, &matrix, threshold,
        ))
    };

    let context_scan_params: Vec<crate::protein_lookup::ProteinContextScanParams> = flat
        .iter()
        .enumerate()
        .map(|(gi, &(qi, fj))| {
            let (begin, end, _frame, search_space, ctx) = qplans[qi].frames[fj];
            let frame_ungapped_kbp = &qplans[qi].ungapped_kbp_array[ctx];
            let frame_gap_trigger_raw =
                raw_gap_trigger_bits(crate::stat::BLAST_GAP_TRIGGER_PROT, frame_ungapped_kbp);
            let cutoff_score = protein_prelim_seed_cutoff(
                frame_gap_trigger_raw,
                prelim_evalue,
                &prot_kbp,
                gumbel_blk.as_ref(),
                (end - begin) as i32,
                avg_subject_length as i32,
                search_space,
            );
            crate::protein_lookup::ProteinContextScanParams {
                query_start: attrib_info.contexts[gi].query_offset as usize,
                x_dropoff: raw_ungapped_xdrop_bits(params.x_drop_ungapped, frame_ungapped_kbp),
                cutoff_score: cutoff_score.max(1),
            }
        })
        .collect();

    let mut results: Vec<Vec<Option<SearchResult>>> = (0..queries.len())
        .map(|_| vec![None; db.num_oids as usize])
        .collect();
    let mut scratch = ProteinScratch::new();
    if let Some(lookup) = lookup.as_ref() {
        for oid in 0..db.num_oids {
            let subj_raw = db.get_sequence(oid);
            let subj_len = db.get_seq_len(oid) as usize;
            if subj_len == 0 {
                continue;
            }
            let subj_aa = &subj_raw[..subj_len];
            scratch.diag_buf.clear();
            let scan_hits = crate::protein_lookup::protein_scan_with_table_reuse_context_params(
                &concat,
                subj_aa,
                &matrix,
                lookup,
                params.two_hit_window as i32,
                &mut scratch.diag_buf,
                &context_scan_params,
            );
            if scan_hits.is_empty() {
                continue;
            }
            let mut per_flat: Vec<Vec<crate::protein_lookup::ProteinHit>> =
                (0..flat.len()).map(|_| Vec::new()).collect();
            for mut h in scan_hits {
                let gi = crate::queryinfo::bsearch_context_info(h.query_start as i32, &attrib_info);
                if gi < 0 || gi as usize >= flat.len() {
                    continue;
                }
                let gi = gi as usize;
                let off = attrib_info.contexts[gi].query_offset as usize;
                h.query_start -= off;
                h.query_end = h.query_end.saturating_sub(off);
                h.gapped_start_q = h.gapped_start_q.saturating_sub(off);
                per_flat[gi].push(h);
            }
            let mut accession: Option<String> = None;
            let mut title: Option<String> = None;
            let mut per_query_internal_hsps: Vec<crate::hspstream::BlastHSPList> = (0..queries
                .len())
                .map(|qi| blast_hsp_list_new(qi as i32))
                .collect();
            let max_query_prot_lens: Vec<usize> = qplans
                .iter()
                .map(|qp| {
                    qp.frames
                        .iter()
                        .map(|&(begin, end, _, _, _)| end.saturating_sub(begin))
                        .max()
                        .unwrap_or(0)
                        .max(1)
                })
                .collect();
            let mut prelim_tree_query: Option<usize> = None;
            for (gi, &(qi, fj)) in flat.iter().enumerate() {
                let frame_ungapped = std::mem::take(&mut per_flat[gi]);
                if frame_ungapped.is_empty() {
                    continue;
                }
                let (begin, end, frame, search_space, ctx) = qplans[qi].frames[fj];
                let prot: &[u8] = &qplans[qi].translation[begin..end];
                let frame_gap_trigger_raw = raw_gap_trigger_bits(
                    crate::stat::BLAST_GAP_TRIGGER_PROT,
                    &qplans[qi].ungapped_kbp_array[ctx],
                );
                let blastx_seed_cutoff = protein_prelim_seed_cutoff(
                    frame_gap_trigger_raw,
                    prelim_evalue,
                    &prot_kbp,
                    gumbel_blk.as_ref(),
                    prot.len() as i32,
                    avg_subject_length as i32,
                    search_space,
                );
                let blastx_hit_cutoff = protein_eval_cutoff(
                    prelim_evalue,
                    &prot_kbp,
                    gumbel_blk.as_ref(),
                    prot.len() as i32,
                    avg_subject_length as i32,
                    search_space,
                );
                let blastx_redo_cutoff = blastx_hit_cutoff;
                let blastx_hit_cutoff = if translated_sum_stats && !params.ungapped {
                    protein_sum_stats_hit_cutoff(
                        blastx_hit_cutoff,
                        &prot_kbp,
                        qplans[qi].link_avg_query_length,
                        avg_subject_length as i32,
                    )
                } else {
                    blastx_hit_cutoff
                };
                if params.ungapped {
                    let mut ug = frame_ungapped;
                    ug.iter_mut().filter(|uh| uh.score >= 1).for_each(|uh| {
                        uh.qseq = Some(
                            prot[uh.query_start..uh.query_end]
                                .iter()
                                .map(|&aa| ncbistdaa_to_aminoacid_char(aa))
                                .collect(),
                        );
                        uh.sseq = Some(
                            subj_aa[uh.subject_start..uh.subject_end]
                                .iter()
                                .map(|&aa| ncbistdaa_to_aminoacid_char(aa))
                                .collect(),
                        );
                    });
                    if accession.is_none() {
                        accession = Some(
                            db.get_accession(oid)
                                .unwrap_or_else(|| format!("oid_{}", oid)),
                        );
                        title = Some(String::from_utf8_lossy(db.get_header(oid)).to_string());
                    }
                    for hsp in process_blastx_frame_hsps(
                        ug,
                        &qplans[qi].translation,
                        prot,
                        frame,
                        subj_aa,
                        subj_len,
                        &qplans[qi].query_info,
                        ctx,
                        search_space,
                        &matrix,
                        &prot_kbp,
                        &qplans[qi].ungapped_kbp_array[ctx],
                        &gumbel_blk,
                        x_drop_final,
                        blastx_hit_cutoff,
                        blastx_redo_cutoff,
                        blastx_seed_cutoff,
                        translated_sum_stats,
                        params,
                        &mut scratch,
                    ) {
                        push_hsp_for_subject(
                            &mut results[qi],
                            oid,
                            title.as_deref().unwrap(),
                            accession.as_deref().unwrap(),
                            subj_len,
                            &[],
                            hsp,
                        );
                    }
                } else {
                    if prelim_tree_query != Some(qi) {
                        scratch.prelim_tree.reset_for_query(
                            max_query_prot_lens[qi] as i32 + 1,
                            subj_aa.len() as i32 + 1,
                        );
                        prelim_tree_query = Some(qi);
                    }
                    let phits = protein_prelim_gapped_hits(
                        prot,
                        subj_aa,
                        &matrix,
                        &frame_ungapped,
                        params.gap_open,
                        params.gap_extend,
                        x_drop_gapped,
                        x_drop_final,
                        blastx_seed_cutoff.max(1),
                        blastx_hit_cutoff.max(1),
                        params.chaining && params.matrix == MatrixType::Blosum62,
                        &mut scratch,
                        false,
                        ctx as i32,
                        0,
                    );
                    for hit in phits {
                        blast_hsp_list_push(
                            &mut per_query_internal_hsps[qi],
                            protein_hit_to_blast_hsp(hit, ctx as i32, frame, 0),
                        );
                    }
                }
            }
            if !params.ungapped {
                for qi in 0..queries.len() {
                    if per_query_internal_hsps[qi].hspcnt == 0 {
                        continue;
                    }
                    blastx_traceback_subject_hsp_list(
                        &qplans[qi].translation,
                        &qplans[qi].frame_offsets,
                        subj_aa,
                        &matrix,
                        &mut per_query_internal_hsps[qi],
                        params.gap_open,
                        params.gap_extend,
                        x_drop_final,
                        &mut scratch,
                        max_query_prot_lens[qi],
                    );
                    per_query_internal_hsps[qi].hsp_array.sort_by(|a, b| {
                        let Some(a) = a.as_ref() else {
                            return std::cmp::Ordering::Greater;
                        };
                        let Some(b) = b.as_ref() else {
                            return std::cmp::Ordering::Less;
                        };
                        a.context
                            .cmp(&b.context)
                            .then((a.query.frame as i32).cmp(&(b.query.frame as i32)))
                            .then_with(|| compare_blast_hsps_by_score(a, b))
                    });
                    let items = std::mem::take(&mut per_query_internal_hsps[qi].hsp_array);
                    let mut group_start = 0usize;
                    while group_start < items.len() {
                        let Some(first) = items[group_start].as_ref() else {
                            group_start += 1;
                            continue;
                        };
                        let context = first.context.max(0) as usize;
                        if context + 1 >= qplans[qi].frame_offsets.len() {
                            group_start += 1;
                            continue;
                        }
                        let frame = first.query.frame as i32;
                        let query_range = (qplans[qi].frame_offsets[context] + 1).max(0) as usize
                            ..qplans[qi].frame_offsets[context + 1].max(0) as usize;
                        let search_space =
                            qplans[qi].query_info.contexts[context].eff_searchsp.max(1) as f64;
                        let frame_ungapped_kbp = &qplans[qi].ungapped_kbp_array[context];
                        let frame_gap_trigger_raw = raw_gap_trigger_bits(
                            crate::stat::BLAST_GAP_TRIGGER_PROT,
                            frame_ungapped_kbp,
                        );
                        let seed_cutoff = protein_prelim_seed_cutoff(
                            frame_gap_trigger_raw,
                            prelim_evalue,
                            &prot_kbp,
                            gumbel_blk.as_ref(),
                            query_range.len() as i32,
                            avg_subject_length as i32,
                            search_space,
                        );
                        let mut hit_cutoff = protein_eval_cutoff(
                            prelim_evalue,
                            &prot_kbp,
                            gumbel_blk.as_ref(),
                            query_range.len() as i32,
                            avg_subject_length as i32,
                            search_space,
                        );
                        if translated_sum_stats {
                            hit_cutoff = protein_sum_stats_hit_cutoff(
                                hit_cutoff,
                                &prot_kbp,
                                qplans[qi].link_avg_query_length,
                                avg_subject_length as i32,
                            );
                        }
                        let mut group_end = group_start + 1;
                        while group_end < items.len()
                            && items[group_end].as_ref().is_some_and(|hsp| {
                                hsp.context.max(0) as usize == context
                                    && hsp.query.frame as i32 == frame
                            })
                        {
                            group_end += 1;
                        }
                        let prot = &qplans[qi].translation[query_range.clone()];
                        let phits: Vec<_> = items[group_start..group_end]
                            .iter()
                            .filter_map(|item| item.as_ref().map(blast_hsp_to_protein_hit))
                            .collect();
                        for hsp in process_blastx_frame_hsps(
                            phits,
                            &qplans[qi].translation,
                            prot,
                            frame,
                            subj_aa,
                            subj_len,
                            &qplans[qi].query_info,
                            context,
                            search_space,
                            &matrix,
                            &prot_kbp,
                            &qplans[qi].ungapped_kbp_array[context],
                            &gumbel_blk,
                            x_drop_final,
                            hit_cutoff,
                            hit_cutoff,
                            seed_cutoff,
                            translated_sum_stats,
                            params,
                            &mut scratch,
                        ) {
                            if accession.is_none() {
                                accession = Some(
                                    db.get_accession(oid)
                                        .unwrap_or_else(|| format!("oid_{}", oid)),
                                );
                                title =
                                    Some(String::from_utf8_lossy(db.get_header(oid)).to_string());
                            }
                            push_hsp_for_subject(
                                &mut results[qi],
                                oid,
                                title.as_deref().unwrap(),
                                accession.as_deref().unwrap(),
                                subj_len,
                                &[],
                                hsp,
                            );
                        }
                        group_start = group_end;
                    }
                }
            }
        }
    }

    (0..queries.len())
        .map(|qi| {
            let mut qres: Vec<SearchResult> = std::mem::take(&mut results[qi])
                .into_iter()
                .flatten()
                .collect();
            if translated_sum_stats && !(params.comp_adjust == 0 && !params.ungapped) {
                apply_blastx_linked_sum_stats(
                    &mut qres,
                    &qplans[qi].query_info,
                    &prot_kbp,
                    gumbel_blk.as_ref(),
                    stats_db_len,
                    params.comp_adjust == 0,
                    params.max_intron_length,
                );
            }
            reap_hsps_by_prelim_evalue(&mut qres, prelim_evalue);
            apply_api_evalue_filter(&mut qres, params.evalue_threshold);
            for result in &mut qres {
                result.hsps.sort_by(compare_hsps_by_evalue_then_score);
            }
            apply_api_min_score_filter(&mut qres, params.min_score);
            if let Some(culling_limit) = params.culling_limit {
                apply_api_culling_limit(&mut qres, culling_limit, crate::program::BLASTX);
            }
            apply_api_max_hsps_limit(&mut qres, max_hsps);
            qres.sort_by(compare_search_results);
            if qres.len() > params.max_target_seqs {
                qres.truncate(params.max_target_seqs);
            }
            fill_protein_midlines(&mut qres, ProteinMidlineStyle::AnyMismatchPlus);
            qres
        })
        .collect()
}

/// Run a tblastn search (protein query vs translated nucleotide database).
/// blast-rs: Public tblastn API pipeline assembled from ported lower-level pieces;
/// not a direct NCBI C port.
/// blast-rs: per-(subject frame) tblastn result body, extracted from `tblastn`
/// so the cross-query batch driver can reuse it. Returns this query x subject-
/// frame's HSPs (subject coords frame-oriented to nucleotide); caller accumulates.
#[allow(clippy::too_many_arguments)]
fn process_tblastn_frame_hsps(
    phits: Vec<crate::protein_lookup::ProteinHit>,
    query_aa: &[u8],
    _prot: &[u8],
    subject_source_ncbi4na: Option<&[u8]>,
    frame: i32,
    _subject_len: usize,
    subj_prot_len: usize,
    tblastn_spouge_subject_len: usize,
    query_info: &crate::queryinfo::QueryInfo,
    query_context: usize,
    search_space: f64,
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
    prot_kbp: &crate::stat::KarlinBlk,
    gumbel_blk: &Option<crate::stat::GumbelBlk>,
    x_drop_final: i32,
    hit_cutoff: i32,
    redo_cutoff_score_min: i32,
    link_word_cutoff_score_min: i32,
    translated_sum_stats: bool,
    params: &SearchParams,
    _scratch: &mut ProteinScratch,
) -> Vec<Hsp> {
    let mut out_hsps: Vec<Hsp> = Vec::new();
    let matrix_name = protein_matrix_name(params.matrix);
    if phits.is_empty() {
        return out_hsps;
    }
    if !params.ungapped && params.comp_adjust > 0 {
        let Some(subject_ncbi4na) = subject_source_ncbi4na else {
            return out_hsps;
        };
        let mut this_match = blast_hsp_list_new(query_context as i32);
        for hit in phits {
            blast_hsp_list_push(
                &mut this_match,
                protein_hit_to_blast_hsp(hit, query_context as i32, 0, frame),
            );
        }
        this_match.query_index = query_context as i32;
        let mut matrix_info = crate::blast_kappa::BlastMatrixInfo::default();
        let ideal_lambda = crate::stat::protein_ideal_ungapped_kbp_for_matrix(matrix_name).lambda;
        if crate::blast_kappa::s_matrix_info_init(
            &mut matrix_info,
            matrix_name,
            ideal_lambda,
            COMPO_ADJUST_SCALE_FACTOR,
        ) != 0
        {
            return out_hsps;
        }
        let scaled_gap_open =
            crate::math::blast_nint(params.gap_open as f64 * COMPO_ADJUST_SCALE_FACTOR);
        let scaled_gap_extend =
            crate::math::blast_nint(params.gap_extend as f64 * COMPO_ADJUST_SCALE_FACTOR);
        let scaled_x_drop_final =
            crate::math::blast_nint(x_drop_final as f64 * COMPO_ADJUST_SCALE_FACTOR);
        let align_params = crate::blast_kappa::blast_redo_align_params_new(
            matrix_info,
            crate::blast_kappa::BlastCompoGappingParams {
                gap_open: scaled_gap_open as i32,
                gap_extend: scaled_gap_extend as i32,
                decline_align: i32::MIN,
                x_dropoff: scaled_x_drop_final as i32,
                context: None,
            },
            crate::blast_kappa::CompoAdjustMode::from_u8(params.comp_adjust),
            COMPO_ADJUST_SCALE_FACTOR,
            false,
            false,
            true,
            query_aa.len() as i32,
            composition_scaled_cutoff(redo_cutoff_score_min),
            params.evalue_threshold,
            true,
            0.0,
        );
        let scaled_kbp = crate::stat::KarlinBlk {
            lambda: prot_kbp.lambda / COMPO_ADJUST_SCALE_FACTOR,
            k: prot_kbp.log_k.exp(),
            log_k: prot_kbp.log_k,
            h: prot_kbp.h,
            round_down: prot_kbp.round_down,
        };
        let link_gbp_db_length = gumbel_blk
            .as_ref()
            .map(|gbp| gbp.db_length.max(1))
            .unwrap_or_else(|| _subject_len.max(1) as i64);
        let score_block = crate::link_hsps::LinkScoreBlock {
            kbp: vec![scaled_kbp.clone(); query_info.contexts.len().max(1)],
            kbp_gap: vec![scaled_kbp; query_info.contexts.len().max(1)],
            gbp: gumbel_blk.clone(),
            link_gbp_db_length: Some(link_gbp_db_length.max(1)),
            recompute_evalues_before_uneven_linking: true,
        };
        let link_params = crate::link_hsps::LinkHSPParameters {
            gap_prob: crate::stat::BLAST_GAP_PROB_GAPPED,
            gap_decay_rate: crate::stat::BLAST_GAP_DECAY_RATE_GAPPED,
            longest_intron: translated_gapped_link_longest_intron(params.max_intron_length),
            cutoff_small_gap: link_word_cutoff_score_min,
            ..crate::link_hsps::LinkHSPParameters::default()
        };
        let link_context = crate::blast_kappa::HitlistLinkContext {
            query_info,
            query_context,
            score_block: &score_block,
            link_params: &link_params,
            gapped_calculation: true,
        };
        let scoring_options = crate::options::ScoringOptions {
            matrix_path: None,
            reward: 0,
            penalty: 0,
            gap_open: params.gap_open,
            gap_extend: params.gap_extend,
            shift_pen: i16::MAX as i32,
            gapped_calculation: true,
            complexity_adjusted_scoring: false,
            matrix_name: Some(matrix_name.to_string()),
            is_ooframe: false,
            program_number: crate::program::UNDEFINED,
        };
        let mut scoring = crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let mut matrix_vec: Vec<Vec<i32>> = matrix.iter().map(|row| row.to_vec()).collect();
        let mut kbp_gap = vec![prot_kbp.clone(); query_info.contexts.len().max(1)];
        let mut saved = crate::blast_kappa::BlastKappaSavedParameters::default();
        let mut results =
            crate::hspstream::HspResults::new(query_info.contexts.len().max(1) as i32);
        let subject = crate::blast_kappa::BlastRedoInMemorySubject {
            subject_source: subject_ncbi4na,
            reward: 0,
            penalty: 0,
            genetic_code: crate::util::lookup_genetic_code(params.db_gencode),
            smith_waterman: false,
            expect_value: params.evalue_threshold,
            hitlist_size: params
                .max_hsps
                .unwrap_or(crate::options::HITLIST_SIZE as usize) as i32,
            inclusion_ethresh: f64::INFINITY,
            link_context: Some(&link_context),
        };
        let mut legacy_this_match = crate::hspstream::HspList::new(this_match.oid);
        let status = crate::blast_kappa::blast_redo_alignment_core_mt(
            crate::program::TBLASTN,
            1,
            query_aa,
            query_info,
            &mut kbp_gap,
            &mut matrix_vec,
            &mut scoring,
            &align_params,
            &mut saved,
            &mut legacy_this_match,
            crate::blast_kappa::BlastRedoAlignmentSource::InMemorySubjectHspList {
                hsp_list: &mut this_match,
                subject,
            },
            &mut results,
        );
        if status != 0 {
            return out_hsps;
        }
        this_match.hsp_array.retain(|hsp| {
            hsp.as_ref()
                .is_some_and(|hsp| hsp.evalue <= params.evalue_threshold)
        });
        this_match.hspcnt = this_match.hsp_array.len() as i32;
        return this_match
            .hsp_array
            .into_iter()
            .flatten()
            .map(|hsp| tblastn_blast_hsp_to_api_hsp(&hsp, prot_kbp))
            .collect();
    }
    let final_phits = phits;
    let use_adj_matrix = false;
    let lambda_ratio_opt = None::<f64>;
    let comp_adjust_method_id = 0;
    for mut ph in final_phits {
        if !params.ungapped {
            if let Some(scaled_score) = ph.scaled_score {
                if scaled_score < composition_scaled_cutoff(hit_cutoff) {
                    continue;
                }
            } else if ph.score < hit_cutoff {
                continue;
            }
        }
        let bit_score = prot_kbp.raw_to_bit(ph.score);
        let (e_score_i32, e_lambda) = match ph.scaled_score {
            Some(s) if params.comp_adjust > 0 => (s, prot_kbp.lambda / COMPO_ADJUST_SCALE_FACTOR),
            _ => (ph.score, prot_kbp.lambda),
        };
        let e_kbp = crate::stat::KarlinBlk {
            lambda: e_lambda,
            k: prot_kbp.k,
            log_k: prot_kbp.log_k,
            h: prot_kbp.h,
            round_down: prot_kbp.round_down,
        };
        let evalue = if let Some(adjusted_evalue) = ph.adjusted_evalue {
            adjusted_evalue
        } else if let Some(ref gbp) = gumbel_blk {
            let base_ev = if params.comp_adjust == 0 {
                // NCBI applies sum-stats gap-decay divisor for
                // translated programs regardless of comp_adjust
                // (`BLAST_GapDecayDivisor` from
                // `blast_stat.c:5360`). Without this, single-HSP
                // tblastn evalues are ~10% smaller than NCBI's.
                let mut evalue = crate::stat::blast_spouge_sto_e(
                    e_score_i32,
                    Some(&e_kbp),
                    Some(gbp),
                    query_aa.len() as i32,
                    tblastn_spouge_subject_len as i32,
                );
                if translated_sum_stats {
                    evalue /= crate::stat::blast_gap_decay_divisor(
                        crate::stat::BLAST_GAP_DECAY_RATE_GAPPED,
                        1,
                    );
                }
                evalue
            } else {
                let mut evalue = crate::stat::blast_spouge_sto_e(
                    e_score_i32,
                    Some(&e_kbp),
                    Some(gbp),
                    query_aa.len() as i32,
                    subj_prot_len as i32,
                );
                if translated_sum_stats {
                    evalue /= crate::stat::blast_gap_decay_divisor(
                        crate::stat::BLAST_GAP_DECAY_RATE_GAPPED,
                        1,
                    );
                }
                evalue
            };
            if use_adj_matrix {
                base_ev
            } else if let Some(lr) = lambda_ratio_opt {
                let scaled_kbp = crate::stat::KarlinBlk {
                    lambda: e_kbp.lambda / lr,
                    k: e_kbp.k,
                    log_k: e_kbp.log_k,
                    h: e_kbp.h,
                    round_down: e_kbp.round_down,
                };
                if params.comp_adjust == 0 {
                    crate::stat::blast_spouge_sto_e(
                        e_score_i32,
                        Some(&scaled_kbp),
                        Some(gbp),
                        query_aa.len() as i32,
                        tblastn_spouge_subject_len as i32,
                    )
                } else {
                    let mut evalue = crate::stat::blast_spouge_sto_e(
                        e_score_i32,
                        Some(&scaled_kbp),
                        Some(gbp),
                        query_aa.len() as i32,
                        subj_prot_len as i32,
                    );
                    if translated_sum_stats {
                        evalue /= crate::stat::blast_gap_decay_divisor(
                            crate::stat::BLAST_GAP_DECAY_RATE_GAPPED,
                            1,
                        );
                    }
                    evalue
                }
            } else {
                base_ev
            }
        } else {
            let raw_evalue = e_kbp.raw_to_evalue(e_score_i32, search_space);
            if use_adj_matrix {
                raw_evalue
            } else if let Some(lr) = lambda_ratio_opt {
                let scaled_lambda = e_kbp.lambda / lr;
                let scaled_kbp = crate::stat::KarlinBlk {
                    lambda: scaled_lambda,
                    k: e_kbp.k,
                    log_k: e_kbp.log_k,
                    h: e_kbp.h,
                    round_down: e_kbp.round_down,
                };
                scaled_kbp.raw_to_evalue(e_score_i32, search_space)
            } else {
                raw_evalue
            }
        };
        let (subject_start, subject_end) =
            crate::util::protein_to_oriented_nuc_coords(ph.subject_start, ph.subject_end, frame);
        let (subject_gapped_start, _) = crate::util::protein_to_oriented_nuc_coords(
            ph.gapped_start_s,
            ph.gapped_start_s.saturating_add(1),
            frame,
        );
        let (q_aln, s_aln) = match (ph.qseq.take(), ph.sseq.take()) {
            (Some(qs), Some(ss)) => (qs.into_bytes(), ss.into_bytes()),
            _ => (Vec::new(), Vec::new()),
        };
        let hsp = Hsp {
            score: ph.score,
            bit_score,
            evalue,
            query_start: ph.query_start,
            query_end: ph.query_end,
            subject_start,
            subject_end,
            query_gapped_start: ph.gapped_start_q,
            subject_gapped_start,
            query_link_start: ph.query_start,
            query_link_end: ph.query_end,
            query_link_gapped_start: ph.gapped_start_q,
            subject_link_start: ph.subject_start,
            subject_link_end: ph.subject_end,
            subject_link_gapped_start: ph.gapped_start_s,
            link_score: ph.scaled_score,
            link_lambda: ph
                .scaled_score
                .map(|_| prot_kbp.lambda / COMPO_ADJUST_SCALE_FACTOR),
            num_identities: ph.num_ident as usize,
            num_gaps: ph.gap_opens as usize,
            alignment_length: ph.align_length as usize,
            query_aln: q_aln,
            midline: Vec::new(),
            subject_aln: s_aln,
            query_frame: 0,
            subject_frame: frame,
            num_links: 1,
            comp_adjust_method: comp_adjust_method_id,
        };
        out_hsps.push(hsp);
    }
    out_hsps
}

pub fn tblastn(db: &BlastDb, query: &[u8], params: &SearchParams) -> Vec<SearchResult> {
    if query.is_empty() {
        return Vec::new();
    }
    let query_aa = encode_protein_query(
        query,
        params.filter_low_complexity,
        params.seg_window,
        params.seg_locut,
        params.seg_hicut,
    );
    if crate::composition::blast_read_aa_composition(&query_aa, AA_SIZE).1 == 0 {
        return Vec::new();
    }
    let matrix = *get_matrix(params.matrix);
    let matrix_name = protein_matrix_name(params.matrix);
    let word_size = params.word_size.clamp(2, 7);
    let threshold = params.word_threshold.unwrap_or_else(|| {
        suggested_word_threshold(params.matrix, crate::program::TBLASTN, params.word_size)
    });

    let prot_kbp = protein_kbp_for_matrix(params.matrix, params.gap_open, params.gap_extend);
    // H2: NCBI `Blast_ScoreBlkKbpUngappedCalc` (blast_stat.c:2737) computes the
    // ungapped Karlin block from the (single, untranslated) protein query's
    // amino-acid composition. tblastn is NOT in the `check_ideal` set
    // (blastx/tblastx/rps_tblastn only), so the ideal-Lambda cap is NOT applied
    // here. This feeds x_drop_ungapped, gap_trigger, and the ungapped-HSP
    // statistics.
    let ungapped_kbp = crate::stat::query_specific_protein_ungapped_kbp_for_matrix(
        &query_aa,
        matrix_name,
        &matrix,
    );
    let x_drop_ungapped = raw_ungapped_xdrop_bits(params.x_drop_ungapped, &ungapped_kbp);
    let x_drop_gapped = raw_gapped_xdrop_bits(params.x_drop_gapped, &prot_kbp);
    let x_drop_final = raw_gapped_xdrop_bits(params.x_drop_final, &prot_kbp);
    let gap_trigger_raw = raw_gap_trigger_bits(crate::stat::BLAST_GAP_TRIGGER_PROT, &ungapped_kbp);
    let total_subj_len: usize = (0..db.num_oids)
        .map(|oid| db.get_seq_len(oid) as usize)
        .sum();
    let min_subject_length = (0..db.num_oids)
        .map(|oid| db.get_seq_len(oid) as usize)
        .min()
        .unwrap_or(1);
    let stats_db_len = statistical_db_length(params, total_subj_len as i64);
    // BLAST_CalcEffLengths handles `db_length /= 3` internally for tblastn,
    // so effective search space uses translated subject letters. The
    // no-composition Spouge path keeps the raw database length in the Gumbel
    // block for the covered subject-mode parity fixtures; the
    // composition-adjusted path preserves the translated DB length used by the
    // existing parity fixtures.
    let translated_total_subj_len = (stats_db_len as usize / 3).max(1);
    // NCBI passes BlastSeqSrcGetMinSeqLen (translated to protein units for
    // tblastn) into BlastHitSavingParametersNew. The downstream cutoff code
    // names that parameter avg_subject_length, but it is the minimum subject
    // length in the database.
    let cutoff_subject_length = (min_subject_length / crate::util::CODON_LENGTH).max(1);
    let mut query_info = crate::queryinfo::QueryInfo::new_blastp(&[query_aa.len()]);
    let scoring_options = crate::options::ScoringOptions {
        matrix_path: None,
        reward: 0,
        penalty: 0,
        gap_open: params.gap_open,
        gap_extend: params.gap_extend,
        shift_pen: i16::MAX as i32,
        // NCBI `-ungapped` propagates `gapped_calculation=FALSE` into
        // `BLAST_CalcEffLengths`, which uses ungapped alpha/beta for
        // length_adjustment (`blast_stat.c:3128`). Without this the
        // effective search space comes from gapped alpha/beta and our
        // ungapped evalues are systematically lower than NCBI's.
        gapped_calculation: !params.ungapped,
        complexity_adjusted_scoring: false,
        matrix_name: Some(matrix_name.to_string()),
        is_ooframe: false,
        program_number: crate::program::UNDEFINED,
    };
    let eff_params = crate::parameters::EffectiveLengthsParameters {
        options: effective_lengths_options(params),
        real_db_length: stats_db_len,
        real_num_seqs: db.num_oids as i32,
    };
    let kbp_array = vec![prot_kbp.clone()];
    let kbp_std_array = vec![ungapped_kbp.clone()];
    crate::blast_setup::blast_calc_eff_lengths(
        crate::program::TBLASTN,
        &scoring_options,
        &eff_params,
        &kbp_array,
        &kbp_std_array,
        matrix_name,
        &mut query_info,
    );
    let search_space = query_info.contexts[0].eff_searchsp.max(1) as f64;
    let translated_sum_stats = translated_sum_stats_enabled(params);
    let prelim_evalue = composition_prelim_evalue(params);
    let translated_gumbel_db_len = translated_total_subj_len as i64;
    // NCBI sets `sbp->gbp->db_length = total_db_length / 3` for TBLASTN
    // before `BLAST_CalcEffLengths` (`blast_setup.c:913-919`).
    // Zero out Gumbel when -searchsp overrides the search space (see blastp).
    let gumbel_blk = if params.effective_search_space > 0 {
        None
    } else {
        protein_gumbel_for_matrix(
            params.matrix,
            params.gap_open,
            params.gap_extend,
            translated_gumbel_db_len,
        )
    };
    let tblastn_link_word_cutoff = protein_prelim_seed_cutoff(
        gap_trigger_raw,
        prelim_evalue,
        &prot_kbp,
        gumbel_blk.as_ref(),
        query_aa.len() as i32,
        cutoff_subject_length as i32,
        search_space,
    );
    let _ = translated_sum_stats; // still used below in Spouge selection
    let max_hsps = params.max_hsps;

    // Build lookup table once per query.
    let lookup_table =
        crate::protein_lookup::ProteinLookupTable::build(&query_aa, word_size, &matrix, threshold);

    let results = map_database_oids_init(db, params, ProteinScratch::new, |scratch, oid| {
        let mut result: Option<SearchResult> = None;
        // Defer defline parsing until this subject produces a reported HSP.
        let mut accession: Option<String> = None;
        let mut title: Option<String> = None;
        let subject_packed = db.get_sequence(oid);
        let subject_len = db.get_seq_len(oid) as usize;
        let mut subject_stat_length = subject_len as i32;

        let (translation_buffer, frame_offsets, subject_ncbi4na) = match db.get_ambiguity_data(oid)
        {
            Some(amb) => {
                let subj_blastna = crate::search::decode_packed_ncbi2na_with_ambiguity(
                    subject_packed,
                    subject_len,
                    amb,
                );
                let subj_ncbi4na = crate::encoding::blastna_to_ncbi4na_sequence(&subj_blastna);
                let translations = crate::util::blast_get_all_translations(
                    &subj_ncbi4na,
                    subj_ncbi4na.len(),
                    crate::util::lookup_genetic_code(params.db_gencode),
                );
                (translations.0, translations.1, subj_ncbi4na)
            }
            None => {
                let subj_blastna =
                    crate::search::decode_packed_ncbi2na(subject_packed, subject_len);
                let subj_ncbi4na = crate::encoding::blastna_to_ncbi4na_sequence(&subj_blastna);
                let translations = crate::util::blast_get_all_translations_packed(
                    subject_packed,
                    subject_len,
                    crate::util::lookup_genetic_code(params.db_gencode),
                );
                (translations.0, translations.1, subj_ncbi4na)
            }
        };
        let max_subj_prot_len = (0..crate::util::NUM_FRAMES)
            .map(|ctx| {
                let begin = (frame_offsets[ctx] + 1) as usize;
                let end = frame_offsets[ctx + 1] as usize;
                end.saturating_sub(begin)
            })
            .max()
            .unwrap_or(0)
            .max(1);
        if !params.ungapped {
            scratch
                .prelim_tree
                .reset_for_query(query_aa.len() as i32 + 1, max_subj_prot_len as i32 + 1);
        }
        let mut subject_internal_hsps = blast_hsp_list_new(0);
        let tblastn_spouge_subject_len = (subject_len / crate::util::CODON_LENGTH).max(1);
        let tblastn_seed_cutoff = protein_prelim_seed_cutoff(
            gap_trigger_raw,
            prelim_evalue,
            &prot_kbp,
            gumbel_blk.as_ref(),
            query_aa.len() as i32,
            cutoff_subject_length as i32,
            search_space,
        );
        let tblastn_hit_cutoff = protein_eval_cutoff(
            prelim_evalue,
            &prot_kbp,
            gumbel_blk.as_ref(),
            query_aa.len() as i32,
            cutoff_subject_length as i32,
            search_space,
        );
        let tblastn_redo_cutoff = tblastn_hit_cutoff;
        let tblastn_hit_cutoff = if translated_sum_stats && !params.ungapped {
            protein_sum_stats_hit_cutoff(
                tblastn_hit_cutoff,
                &prot_kbp,
                query_aa.len() as i32,
                cutoff_subject_length as i32,
            )
        } else {
            tblastn_hit_cutoff
        };

        for ctx in 0..crate::util::NUM_FRAMES {
            let frame = crate::util::blast_context_to_frame(ctx as u32);
            let begin = (frame_offsets[ctx] + 1) as usize;
            let end = frame_offsets[ctx + 1] as usize;
            if begin >= end {
                continue;
            }
            let prot: &[u8] = &translation_buffer[begin..end];
            if prot.len() < word_size {
                continue;
            }
            if crate::composition::blast_read_aa_composition(prot, AA_SIZE).1 == 0 {
                continue;
            }
            let prot_eval: &[u8] = prot;
            let subj_prot_len = prot.len();
            if params.ungapped {
                // NCBI `-ungapped` stops after the protein word-finder's
                // ungapped extensions. Keep all positive HSPs here and let
                // the final e-value threshold decide reporting.
                let phits = protein_alignment_hits_ungapped_only(
                    &query_aa,
                    prot,
                    &matrix,
                    &lookup_table,
                    x_drop_ungapped,
                    params.two_hit_window as i32,
                    1,
                    scratch,
                );
                if phits.is_empty() {
                    continue;
                }
                for hsp in process_tblastn_frame_hsps(
                    phits,
                    &query_aa,
                    prot_eval,
                    Some(&subject_ncbi4na),
                    frame,
                    subject_len,
                    subj_prot_len,
                    tblastn_spouge_subject_len,
                    &query_info,
                    0,
                    search_space,
                    &matrix,
                    &prot_kbp,
                    &gumbel_blk,
                    x_drop_final,
                    tblastn_hit_cutoff,
                    tblastn_redo_cutoff,
                    tblastn_seed_cutoff,
                    translated_sum_stats,
                    params,
                    scratch,
                ) {
                    match &mut result {
                        Some(existing) => existing.hsps.push(hsp),
                        None => {
                            let acc = accession.get_or_insert_with(|| {
                                db.get_accession(oid)
                                    .unwrap_or_else(|| format!("oid_{}", oid))
                            });
                            let ttl = title.get_or_insert_with(|| {
                                String::from_utf8_lossy(db.get_header(oid)).to_string()
                            });
                            result = Some(SearchResult {
                                subject_oid: oid,
                                subject_title: ttl.clone(),
                                subject_accession: acc.clone(),
                                subject_len,
                                hsps: vec![hsp],
                                taxids: vec![],
                            });
                        }
                    }
                }
            } else {
                let ungapped_hits = crate::protein_lookup::protein_scan_with_extend_word(
                    &query_aa,
                    prot,
                    &matrix,
                    &lookup_table,
                    x_drop_ungapped,
                    params.two_hit_window as i32,
                    tblastn_seed_cutoff.max(1),
                    &mut scratch.diag_table,
                );
                scratch.diag_table.exit(prot.len());
                let phits = protein_prelim_gapped_hits(
                    &query_aa,
                    prot_eval,
                    &matrix,
                    &ungapped_hits,
                    params.gap_open,
                    params.gap_extend,
                    x_drop_gapped,
                    x_drop_final,
                    tblastn_seed_cutoff,
                    tblastn_hit_cutoff,
                    params.chaining && params.matrix == MatrixType::Blosum62,
                    scratch,
                    false,
                    0,
                    frame,
                );
                for hit in phits {
                    blast_hsp_list_push(
                        &mut subject_internal_hsps,
                        protein_hit_to_blast_hsp(hit, 0, 0, frame),
                    );
                }
            }
        }
        if !params.ungapped && subject_internal_hsps.hspcnt > 0 {
            if params.comp_adjust > 0 {
                let db_code = crate::util::lookup_genetic_code(params.db_gencode);
                let api_hsps = 'redo: {
                    subject_internal_hsps.query_index = 0;
                    let mut matrix_info = crate::blast_kappa::BlastMatrixInfo::default();
                    let ideal_lambda =
                        crate::stat::protein_ideal_ungapped_kbp_for_matrix(matrix_name).lambda;
                    if crate::blast_kappa::s_matrix_info_init(
                        &mut matrix_info,
                        matrix_name,
                        ideal_lambda,
                        COMPO_ADJUST_SCALE_FACTOR,
                    ) != 0
                    {
                        break 'redo Vec::new();
                    }
                    let scaled_gap_open =
                        crate::math::blast_nint(params.gap_open as f64 * COMPO_ADJUST_SCALE_FACTOR);
                    let scaled_gap_extend = crate::math::blast_nint(
                        params.gap_extend as f64 * COMPO_ADJUST_SCALE_FACTOR,
                    );
                    let scaled_x_drop_final =
                        crate::math::blast_nint(x_drop_final as f64 * COMPO_ADJUST_SCALE_FACTOR);
                    let align_params = crate::blast_kappa::blast_redo_align_params_new(
                        matrix_info,
                        crate::blast_kappa::BlastCompoGappingParams {
                            gap_open: scaled_gap_open as i32,
                            gap_extend: scaled_gap_extend as i32,
                            decline_align: i32::MIN,
                            x_dropoff: scaled_x_drop_final as i32,
                            context: None,
                        },
                        crate::blast_kappa::CompoAdjustMode::from_u8(params.comp_adjust),
                        COMPO_ADJUST_SCALE_FACTOR,
                        false,
                        false,
                        true,
                        query_aa.len() as i32,
                        composition_scaled_cutoff(tblastn_redo_cutoff),
                        params.evalue_threshold,
                        true,
                        0.0,
                    );
                    let scaled_kbp = crate::stat::KarlinBlk {
                        lambda: prot_kbp.lambda / COMPO_ADJUST_SCALE_FACTOR,
                        k: prot_kbp.log_k.exp(),
                        log_k: prot_kbp.log_k,
                        h: prot_kbp.h,
                        round_down: prot_kbp.round_down,
                    };
                    let score_block = crate::link_hsps::LinkScoreBlock {
                        kbp: vec![scaled_kbp.clone(); query_info.contexts.len().max(1)],
                        kbp_gap: vec![scaled_kbp; query_info.contexts.len().max(1)],
                        gbp: gumbel_blk.clone(),
                        link_gbp_db_length: Some(translated_gumbel_db_len.max(1)),
                        recompute_evalues_before_uneven_linking: true,
                    };
                    let link_params = crate::link_hsps::LinkHSPParameters {
                        gap_prob: crate::stat::BLAST_GAP_PROB_GAPPED,
                        gap_decay_rate: crate::stat::BLAST_GAP_DECAY_RATE_GAPPED,
                        longest_intron: translated_gapped_link_longest_intron(
                            params.max_intron_length,
                        ),
                        cutoff_small_gap: tblastn_seed_cutoff,
                        ..crate::link_hsps::LinkHSPParameters::default()
                    };
                    let link_context = crate::blast_kappa::HitlistLinkContext {
                        query_info: &query_info,
                        query_context: 0,
                        score_block: &score_block,
                        link_params: &link_params,
                        gapped_calculation: true,
                    };
                    let scoring_options = crate::options::ScoringOptions {
                        matrix_path: None,
                        reward: 0,
                        penalty: 0,
                        gap_open: params.gap_open,
                        gap_extend: params.gap_extend,
                        shift_pen: i16::MAX as i32,
                        gapped_calculation: true,
                        complexity_adjusted_scoring: false,
                        matrix_name: Some(matrix_name.to_string()),
                        is_ooframe: false,
                        program_number: crate::program::UNDEFINED,
                    };
                    let mut scoring =
                        crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
                    let mut matrix_vec: Vec<Vec<i32>> =
                        matrix.iter().map(|row| row.to_vec()).collect();
                    let mut kbp_gap = vec![prot_kbp.clone(); query_info.contexts.len().max(1)];
                    let mut saved = crate::blast_kappa::BlastKappaSavedParameters::default();
                    let mut results =
                        crate::hspstream::HspResults::new(query_info.contexts.len().max(1) as i32);
                    let subject = crate::blast_kappa::BlastRedoInMemorySubject {
                        subject_source: &subject_ncbi4na,
                        reward: 0,
                        penalty: 0,
                        genetic_code: db_code,
                        smith_waterman: false,
                        expect_value: params.evalue_threshold,
                        hitlist_size: max_hsps.unwrap_or(crate::options::HITLIST_SIZE as usize)
                            as i32,
                        inclusion_ethresh: f64::INFINITY,
                        link_context: Some(&link_context),
                    };
                    let mut legacy_this_match =
                        crate::hspstream::HspList::new(subject_internal_hsps.oid);
                    let status = crate::blast_kappa::blast_redo_alignment_core_mt(
                        crate::program::TBLASTN,
                        1,
                        &query_aa,
                        &query_info,
                        &mut kbp_gap,
                        &mut matrix_vec,
                        &mut scoring,
                        &align_params,
                        &mut saved,
                        &mut legacy_this_match,
                        crate::blast_kappa::BlastRedoAlignmentSource::InMemorySubjectHspList {
                            hsp_list: &mut subject_internal_hsps,
                            subject,
                        },
                        &mut results,
                    );
                    if status != 0 {
                        break 'redo Vec::new();
                    }
                    let mut hsp_list =
                        std::mem::replace(&mut subject_internal_hsps, blast_hsp_list_new(0));
                    hsp_list.hsp_array.retain(|hsp| {
                        hsp.as_ref()
                            .is_some_and(|hsp| hsp.evalue <= params.evalue_threshold)
                    });
                    hsp_list.hspcnt = hsp_list.hsp_array.len() as i32;
                    hsp_list
                        .hsp_array
                        .into_iter()
                        .flatten()
                        .map(|hsp| tblastn_blast_hsp_to_api_hsp(&hsp, &prot_kbp))
                        .collect()
                };
                for api_hsp in api_hsps {
                    match &mut result {
                        Some(existing) => existing.hsps.push(api_hsp),
                        None => {
                            let acc = accession.get_or_insert_with(|| {
                                db.get_accession(oid)
                                    .unwrap_or_else(|| format!("oid_{}", oid))
                            });
                            let ttl = title.get_or_insert_with(|| {
                                String::from_utf8_lossy(db.get_header(oid)).to_string()
                            });
                            result = Some(SearchResult {
                                subject_oid: oid,
                                subject_title: ttl.clone(),
                                subject_accession: acc.clone(),
                                subject_len,
                                hsps: vec![api_hsp],
                                taxids: vec![],
                            });
                        }
                    }
                }
                return result.map(|result| (result, subject_stat_length));
            }
            subject_stat_length = tblastn_traceback_subject_hsp_list(
                &query_aa,
                &translation_buffer,
                &frame_offsets,
                &subject_ncbi4na,
                crate::util::lookup_genetic_code(params.db_gencode),
                &matrix,
                &mut subject_internal_hsps,
                params.gap_open,
                params.gap_extend,
                x_drop_final,
                scratch,
                0,
                max_subj_prot_len,
            );
            if translated_sum_stats && params.comp_adjust == 0 {
                for api_hsp in tblastn_link_internal_hsps_to_api_hsps(
                    &mut subject_internal_hsps,
                    &query_info,
                    0,
                    subject_stat_length,
                    &prot_kbp,
                    gumbel_blk.as_ref(),
                    translated_gumbel_db_len,
                    tblastn_seed_cutoff,
                    params.max_intron_length,
                    params.evalue_threshold,
                ) {
                    match &mut result {
                        Some(existing) => existing.hsps.push(api_hsp),
                        None => {
                            let acc = accession.get_or_insert_with(|| {
                                db.get_accession(oid)
                                    .unwrap_or_else(|| format!("oid_{}", oid))
                            });
                            let ttl = title.get_or_insert_with(|| {
                                String::from_utf8_lossy(db.get_header(oid)).to_string()
                            });
                            result = Some(SearchResult {
                                subject_oid: oid,
                                subject_title: ttl.clone(),
                                subject_accession: acc.clone(),
                                subject_len,
                                hsps: vec![api_hsp],
                                taxids: vec![],
                            });
                        }
                    }
                }
                return result.map(|result| (result, subject_stat_length));
            }
            subject_internal_hsps.hsp_array.sort_by(|a, b| {
                let Some(a) = a.as_ref() else {
                    return std::cmp::Ordering::Greater;
                };
                let Some(b) = b.as_ref() else {
                    return std::cmp::Ordering::Less;
                };
                (a.subject.frame as i32)
                    .cmp(&(b.subject.frame as i32))
                    .then_with(|| compare_blast_hsps_by_score(a, b))
            });
            let mut group_start = 0usize;
            while group_start < subject_internal_hsps.hsp_array.len() {
                let Some(first) = subject_internal_hsps.hsp_array[group_start].as_ref() else {
                    group_start += 1;
                    continue;
                };
                let frame = first.subject.frame as i32;
                let Some(ctx) = translated_frame_to_context(frame) else {
                    group_start += 1;
                    continue;
                };
                let prot_range = (frame_offsets[ctx] + 1).max(0) as usize
                    ..frame_offsets[ctx + 1].max(0) as usize;
                let subj_prot_len = prot_range.len();
                let mut group_end = group_start + 1;
                while group_end < subject_internal_hsps.hsp_array.len()
                    && subject_internal_hsps.hsp_array[group_end]
                        .as_ref()
                        .is_some_and(|hsp| hsp.subject.frame as i32 == frame)
                {
                    group_end += 1;
                }
                let prot_eval = &translation_buffer[prot_range.clone()];
                let phits: Vec<_> = subject_internal_hsps.hsp_array[group_start..group_end]
                    .iter()
                    .filter_map(|item| item.as_ref().map(blast_hsp_to_protein_hit))
                    .collect();
                for hsp in process_tblastn_frame_hsps(
                    phits,
                    &query_aa,
                    prot_eval,
                    Some(&subject_ncbi4na),
                    frame,
                    subject_len,
                    subj_prot_len,
                    tblastn_spouge_subject_len,
                    &query_info,
                    0,
                    search_space,
                    &matrix,
                    &prot_kbp,
                    &gumbel_blk,
                    x_drop_final,
                    tblastn_hit_cutoff,
                    tblastn_redo_cutoff,
                    tblastn_seed_cutoff,
                    translated_sum_stats,
                    params,
                    scratch,
                ) {
                    match &mut result {
                        Some(existing) => existing.hsps.push(hsp),
                        None => {
                            let acc = accession.get_or_insert_with(|| {
                                db.get_accession(oid)
                                    .unwrap_or_else(|| format!("oid_{}", oid))
                            });
                            let ttl = title.get_or_insert_with(|| {
                                String::from_utf8_lossy(db.get_header(oid)).to_string()
                            });
                            result = Some(SearchResult {
                                subject_oid: oid,
                                subject_title: ttl.clone(),
                                subject_accession: acc.clone(),
                                subject_len,
                                hsps: vec![hsp],
                                taxids: vec![],
                            });
                        }
                    }
                }
                group_start = group_end;
            }
        }
        result.map(|result| (result, subject_stat_length))
    });

    let mut results_with_stat_lengths: Vec<(SearchResult, i32)> =
        results.into_iter().flatten().collect();
    let mut subject_stat_lengths: Vec<i32> = Vec::with_capacity(results_with_stat_lengths.len());
    let mut results: Vec<SearchResult> = Vec::with_capacity(results_with_stat_lengths.len());
    for (result, stat_length) in results_with_stat_lengths.drain(..) {
        results.push(result);
        subject_stat_lengths.push(stat_length);
    }
    if translated_sum_stats && !(params.comp_adjust == 0 && !params.ungapped) {
        apply_tblastn_linked_sum_stats(
            &mut results,
            Some(&subject_stat_lengths),
            &query_info,
            &prot_kbp,
            gumbel_blk.as_ref(),
            stats_db_len,
            params.comp_adjust == 0,
            tblastn_link_word_cutoff,
            params.gap_open,
            params.gap_extend,
            params.max_intron_length,
        );
    }
    reap_hsps_by_prelim_evalue(&mut results, prelim_evalue);
    apply_api_evalue_filter(&mut results, params.evalue_threshold);
    for result in &mut results {
        result.hsps.sort_by(compare_hsps_by_evalue_then_score);
    }
    apply_api_min_score_filter(&mut results, params.min_score);
    if let Some(culling_limit) = params.culling_limit {
        apply_api_culling_limit(&mut results, culling_limit, crate::program::TBLASTN);
    }
    apply_api_max_hsps_limit(&mut results, max_hsps);
    results.sort_by(compare_search_results);
    if results.len() > params.max_target_seqs {
        results.truncate(params.max_target_seqs);
    }
    fill_protein_midlines(&mut results, ProteinMidlineStyle::AnyMismatchPlus);
    results
}

/// blast-rs: cross-query tblastn. Concatenates ALL protein queries into ONE
/// lookup table and translates each database subject ONCE (NCBI's concatenated-
/// query engine), instead of re-building + re-translating per query. Byte-
/// identical to looping `tblastn` per query; reuses `process_tblastn_frame_hsps`.
#[allow(clippy::too_many_arguments)]
pub fn tblastn_batch(
    db: &BlastDb,
    queries: &[&[u8]],
    params: &SearchParams,
) -> Vec<Vec<SearchResult>> {
    if queries.is_empty() {
        return Vec::new();
    }
    let matrix = *get_matrix(params.matrix);
    let matrix_name = protein_matrix_name(params.matrix);
    let word_size = params.word_size.clamp(2, 7);
    let threshold = params.word_threshold.unwrap_or_else(|| {
        suggested_word_threshold(params.matrix, crate::program::TBLASTN, params.word_size)
    });
    let prot_kbp = protein_kbp_for_matrix(params.matrix, params.gap_open, params.gap_extend);
    let x_drop_gapped = raw_gapped_xdrop_bits(params.x_drop_gapped, &prot_kbp);
    let x_drop_final = raw_gapped_xdrop_bits(params.x_drop_final, &prot_kbp);
    // H2: gap_trigger, x-drop, and eff-length std array are computed per query
    // from each query's ungapped Karlin block (see `QPlan`).
    let total_subj_len: usize = (0..db.num_oids)
        .map(|oid| db.get_seq_len(oid) as usize)
        .sum();
    let min_subject_length = (0..db.num_oids)
        .map(|oid| db.get_seq_len(oid) as usize)
        .min()
        .unwrap_or(1);
    let stats_db_len = statistical_db_length(params, total_subj_len as i64);
    let translated_total_subj_len = (stats_db_len as usize / 3).max(1);
    let cutoff_subject_length = (min_subject_length / crate::util::CODON_LENGTH).max(1);
    let translated_sum_stats = translated_sum_stats_enabled(params);
    let prelim_evalue = composition_prelim_evalue(params);
    let max_hsps = params.max_hsps;
    let db_code = crate::util::lookup_genetic_code(params.db_gencode);

    let scoring_options = crate::options::ScoringOptions {
        matrix_path: None,
        reward: 0,
        penalty: 0,
        gap_open: params.gap_open,
        gap_extend: params.gap_extend,
        shift_pen: i16::MAX as i32,
        gapped_calculation: !params.ungapped,
        complexity_adjusted_scoring: false,
        matrix_name: Some(matrix_name.to_string()),
        is_ooframe: false,
        program_number: crate::program::UNDEFINED,
    };
    let eff_params = crate::parameters::EffectiveLengthsParameters {
        options: effective_lengths_options(params),
        real_db_length: stats_db_len,
        real_num_seqs: db.num_oids as i32,
    };
    let kbp_array = vec![prot_kbp.clone()];

    struct QPlan {
        aa: Vec<u8>,
        query_info: crate::queryinfo::QueryInfo,
        search_space: f64,
        seed_cutoff: i32,
        hit_cutoff: i32,
        redo_cutoff: i32,
        x_drop_ungapped: i32,
        gumbel_blk: Option<crate::stat::GumbelBlk>,
        link_gbp_db_length: i64,
    }
    let qplans: Vec<QPlan> = queries
        .iter()
        .map(|q| {
            let aa = encode_protein_query(
                q,
                params.filter_low_complexity,
                params.seg_window,
                params.seg_locut,
                params.seg_hicut,
            );
            let q_ungapped_kbp = crate::stat::query_specific_protein_ungapped_kbp_for_matrix(
                &aa,
                matrix_name,
                &matrix,
            );
            let kbp_std_array = vec![q_ungapped_kbp.clone()];
            let mut query_info = crate::queryinfo::QueryInfo::new_blastp(&[aa.len()]);
            crate::blast_setup::blast_calc_eff_lengths(
                crate::program::TBLASTN,
                &scoring_options,
                &eff_params,
                &kbp_array,
                &kbp_std_array,
                matrix_name,
                &mut query_info,
            );
            let search_space = query_info
                .contexts
                .first()
                .map(|c| c.eff_searchsp.max(1) as f64)
                .unwrap_or(1.0);
            let translated_gumbel_db_len = translated_total_subj_len as i64;
            let gumbel_blk = if params.effective_search_space > 0 {
                None
            } else {
                protein_gumbel_for_matrix(
                    params.matrix,
                    params.gap_open,
                    params.gap_extend,
                    translated_gumbel_db_len,
                )
            };
            let q_gap_trigger_raw =
                raw_gap_trigger_bits(crate::stat::BLAST_GAP_TRIGGER_PROT, &q_ungapped_kbp);
            let seed_cutoff = protein_prelim_seed_cutoff(
                q_gap_trigger_raw,
                prelim_evalue,
                &prot_kbp,
                gumbel_blk.as_ref(),
                aa.len() as i32,
                cutoff_subject_length as i32,
                search_space,
            );
            let hit_cutoff = protein_eval_cutoff(
                prelim_evalue,
                &prot_kbp,
                gumbel_blk.as_ref(),
                aa.len() as i32,
                cutoff_subject_length as i32,
                search_space,
            );
            let redo_cutoff = hit_cutoff;
            let hit_cutoff = if translated_sum_stats && !params.ungapped {
                protein_sum_stats_hit_cutoff(
                    hit_cutoff,
                    &prot_kbp,
                    aa.len() as i32,
                    cutoff_subject_length as i32,
                )
            } else {
                hit_cutoff
            };
            QPlan {
                aa,
                query_info,
                search_space,
                seed_cutoff,
                hit_cutoff,
                redo_cutoff,
                x_drop_ungapped: raw_ungapped_xdrop_bits(params.x_drop_ungapped, &q_ungapped_kbp),
                gumbel_blk,
                link_gbp_db_length: translated_gumbel_db_len,
            }
        })
        .collect();

    // Concatenate the protein queries (one inter-query sentinel) into one lookup.
    let lengths: Vec<usize> = qplans.iter().map(|p| p.aa.len()).collect();
    let attrib_info = crate::queryinfo::QueryInfo::new_blastp(&lengths);
    let mut concat: Vec<u8> =
        vec![crate::protein_lookup::PROTEIN_CONCAT_SENTINEL; attrib_info.seq_buf_len().max(1)];
    for (qi, p) in qplans.iter().enumerate() {
        let off = attrib_info.contexts[qi].query_offset as usize;
        concat[off..off + p.aa.len()].copy_from_slice(&p.aa);
    }
    let lookup =
        crate::protein_lookup::ProteinLookupTable::build(&concat, word_size, &matrix, threshold);

    // NCBI's AA word finder uses the active query context's x-drop and saved
    // init-HSP cutoff (`aa_ungapped.c:560,575,588`). Preserve those per-query
    // word parameters through the concatenated-query scanner.
    let context_scan_params: Vec<crate::protein_lookup::ProteinContextScanParams> = qplans
        .iter()
        .enumerate()
        .map(|(qi, p)| crate::protein_lookup::ProteinContextScanParams {
            query_start: attrib_info.contexts[qi].query_offset as usize,
            x_dropoff: p.x_drop_ungapped,
            cutoff_score: p.seed_cutoff.max(1),
        })
        .collect();

    // Per OID (parallel): translate subject once, scan 6 frames against the
    // concatenated-query lookup, attribute each hit to its query, gapped + body.
    let per_oid: Vec<Vec<(usize, SearchResult, i32)>> = map_database_oids_init(
        db,
        params,
        ProteinScratch::new,
        |scratch, oid| {
            let mut out: Vec<(usize, SearchResult, i32)> = Vec::new();
            let subject_packed = db.get_sequence(oid);
            let subject_len = db.get_seq_len(oid) as usize;
            let (translation_buffer, frame_offsets, subject_ncbi4na) = match db
                .get_ambiguity_data(oid)
            {
                Some(amb) => {
                    let subj_blastna = crate::search::decode_packed_ncbi2na_with_ambiguity(
                        subject_packed,
                        subject_len,
                        amb,
                    );
                    let subj_ncbi4na = crate::encoding::blastna_to_ncbi4na_sequence(&subj_blastna);
                    let translations = crate::util::blast_get_all_translations(
                        &subj_ncbi4na,
                        subj_ncbi4na.len(),
                        db_code,
                    );
                    (translations.0, translations.1, subj_ncbi4na)
                }
                None => {
                    let subj_blastna =
                        crate::search::decode_packed_ncbi2na(subject_packed, subject_len);
                    let subj_ncbi4na = crate::encoding::blastna_to_ncbi4na_sequence(&subj_blastna);
                    let translations = crate::util::blast_get_all_translations_packed(
                        subject_packed,
                        subject_len,
                        db_code,
                    );
                    (translations.0, translations.1, subj_ncbi4na)
                }
            };
            let mut per_query_result: Vec<Option<SearchResult>> = vec![None; queries.len()];
            let mut per_query_stat_length: Vec<i32> = vec![subject_len as i32; queries.len()];
            // Defer defline parsing until a subject actually produces a hit
            // (NCBI builds deflines only for reported OIDs, not every subject).
            let mut accession: Option<String> = None;
            let mut title: Option<String> = None;
            let max_query_len = qplans.iter().map(|p| p.aa.len()).max().unwrap_or(0).max(1);
            let max_subj_prot_len = (0..crate::util::NUM_FRAMES)
                .map(|ctx| {
                    let begin = (frame_offsets[ctx] + 1) as usize;
                    let end = frame_offsets[ctx + 1] as usize;
                    end.saturating_sub(begin)
                })
                .max()
                .unwrap_or(0)
                .max(1);
            if !params.ungapped {
                scratch
                    .prelim_tree
                    .reset_for_query(max_query_len as i32 + 1, max_subj_prot_len as i32 + 1);
            }
            let mut per_query_internal_hsps: Vec<crate::hspstream::BlastHSPList> = (0..queries
                .len())
                .map(|qi| blast_hsp_list_new(qi as i32))
                .collect();
            for ctx in 0..crate::util::NUM_FRAMES {
                let frame = crate::util::blast_context_to_frame(ctx as u32);
                let begin = (frame_offsets[ctx] + 1) as usize;
                let end = frame_offsets[ctx + 1] as usize;
                if begin >= end {
                    continue;
                }
                let prot: &[u8] = &translation_buffer[begin..end];
                if prot.len() < word_size {
                    continue;
                }
                if crate::composition::blast_read_aa_composition(prot, AA_SIZE).1 == 0 {
                    continue;
                }
                let subj_prot_len = prot.len();
                let scan_hits = crate::protein_lookup::protein_scan_with_extend_word_context_params(
                    &concat,
                    prot,
                    &matrix,
                    &lookup,
                    params.two_hit_window as i32,
                    &mut scratch.diag_table,
                    &context_scan_params,
                );
                scratch.diag_table.exit(prot.len());
                if scan_hits.is_empty() {
                    continue;
                }
                let mut per_query: Vec<Vec<crate::protein_lookup::ProteinHit>> =
                    (0..queries.len()).map(|_| Vec::new()).collect();
                for mut h in scan_hits {
                    let qi =
                        crate::queryinfo::bsearch_context_info(h.query_start as i32, &attrib_info);
                    if qi < 0 || qi as usize >= queries.len() {
                        continue;
                    }
                    let qi = qi as usize;
                    let off = attrib_info.contexts[qi].query_offset as usize;
                    h.query_start -= off;
                    h.query_end = h.query_end.saturating_sub(off);
                    h.gapped_start_q = h.gapped_start_q.saturating_sub(off);
                    per_query[qi].push(h);
                }
                for qi in 0..queries.len() {
                    let ungapped = std::mem::take(&mut per_query[qi]);
                    if ungapped.is_empty() {
                        continue;
                    }
                    let p = &qplans[qi];
                    if params.ungapped {
                        let mut ug = ungapped;
                        ug.iter_mut().filter(|uh| uh.score >= 1).for_each(|uh| {
                            uh.qseq = Some(
                                p.aa[uh.query_start..uh.query_end]
                                    .iter()
                                    .map(|&aa| ncbistdaa_to_aminoacid_char(aa))
                                    .collect(),
                            );
                            uh.sseq = Some(
                                prot[uh.subject_start..uh.subject_end]
                                    .iter()
                                    .map(|&aa| ncbistdaa_to_aminoacid_char(aa))
                                    .collect(),
                            );
                        });
                        if ug.is_empty() {
                            continue;
                        }
                        for hsp in process_tblastn_frame_hsps(
                            ug,
                            &p.aa,
                            prot,
                            Some(&subject_ncbi4na),
                            frame,
                            subject_len,
                            subj_prot_len,
                            (subject_len / crate::util::CODON_LENGTH).max(1),
                            &p.query_info,
                            0,
                            p.search_space,
                            &matrix,
                            &prot_kbp,
                            &p.gumbel_blk,
                            x_drop_final,
                            p.hit_cutoff,
                            p.redo_cutoff,
                            p.seed_cutoff,
                            translated_sum_stats,
                            params,
                            scratch,
                        ) {
                            match &mut per_query_result[qi] {
                                Some(existing) => existing.hsps.push(hsp),
                                None => {
                                    let acc = accession.get_or_insert_with(|| {
                                        db.get_accession(oid)
                                            .unwrap_or_else(|| format!("oid_{}", oid))
                                    });
                                    let ttl = title.get_or_insert_with(|| {
                                        String::from_utf8_lossy(db.get_header(oid)).to_string()
                                    });
                                    per_query_result[qi] = Some(SearchResult {
                                        subject_oid: oid,
                                        subject_title: ttl.clone(),
                                        subject_accession: acc.clone(),
                                        subject_len,
                                        hsps: vec![hsp],
                                        taxids: vec![],
                                    });
                                }
                            }
                        }
                    } else {
                        let phits = protein_prelim_gapped_hits(
                            &p.aa,
                            prot,
                            &matrix,
                            &ungapped,
                            params.gap_open,
                            params.gap_extend,
                            x_drop_gapped,
                            x_drop_final,
                            p.seed_cutoff,
                            p.hit_cutoff,
                            params.chaining && params.matrix == MatrixType::Blosum62,
                            scratch,
                            false,
                            qi as i32,
                            frame,
                        );
                        for hit in phits {
                            blast_hsp_list_push(
                                &mut per_query_internal_hsps[qi],
                                protein_hit_to_blast_hsp(hit, qi as i32, 0, frame),
                            );
                        }
                    }
                }
            }
            if !params.ungapped {
                for qi in 0..queries.len() {
                    if per_query_internal_hsps[qi].hspcnt == 0 {
                        continue;
                    }
                    for hsp in &mut per_query_internal_hsps[qi].hsp_array {
                        if let Some(hsp) = hsp.as_mut() {
                            hsp.context = 0;
                        }
                    }
                    let p = &qplans[qi];
                    if params.comp_adjust > 0 {
                        let api_hsps = 'redo: {
                            per_query_internal_hsps[qi].query_index = 0;
                            let mut matrix_info = crate::blast_kappa::BlastMatrixInfo::default();
                            let ideal_lambda =
                                crate::stat::protein_ideal_ungapped_kbp_for_matrix(matrix_name)
                                    .lambda;
                            if crate::blast_kappa::s_matrix_info_init(
                                &mut matrix_info,
                                matrix_name,
                                ideal_lambda,
                                COMPO_ADJUST_SCALE_FACTOR,
                            ) != 0
                            {
                                break 'redo Vec::new();
                            }
                            let scaled_gap_open = crate::math::blast_nint(
                                params.gap_open as f64 * COMPO_ADJUST_SCALE_FACTOR,
                            );
                            let scaled_gap_extend = crate::math::blast_nint(
                                params.gap_extend as f64 * COMPO_ADJUST_SCALE_FACTOR,
                            );
                            let scaled_x_drop_final = crate::math::blast_nint(
                                x_drop_final as f64 * COMPO_ADJUST_SCALE_FACTOR,
                            );
                            let align_params = crate::blast_kappa::blast_redo_align_params_new(
                                matrix_info,
                                crate::blast_kappa::BlastCompoGappingParams {
                                    gap_open: scaled_gap_open as i32,
                                    gap_extend: scaled_gap_extend as i32,
                                    decline_align: i32::MIN,
                                    x_dropoff: scaled_x_drop_final as i32,
                                    context: None,
                                },
                                crate::blast_kappa::CompoAdjustMode::from_u8(params.comp_adjust),
                                COMPO_ADJUST_SCALE_FACTOR,
                                false,
                                false,
                                true,
                                p.aa.len() as i32,
                                composition_scaled_cutoff(p.redo_cutoff),
                                params.evalue_threshold,
                                true,
                                0.0,
                            );
                            let scaled_kbp = crate::stat::KarlinBlk {
                                lambda: prot_kbp.lambda / COMPO_ADJUST_SCALE_FACTOR,
                                k: prot_kbp.log_k.exp(),
                                log_k: prot_kbp.log_k,
                                h: prot_kbp.h,
                                round_down: prot_kbp.round_down,
                            };
                            let score_block = crate::link_hsps::LinkScoreBlock {
                                kbp: vec![scaled_kbp.clone(); p.query_info.contexts.len().max(1)],
                                kbp_gap: vec![scaled_kbp; p.query_info.contexts.len().max(1)],
                                gbp: p.gumbel_blk.clone(),
                                link_gbp_db_length: Some(p.link_gbp_db_length.max(1)),
                                recompute_evalues_before_uneven_linking: true,
                            };
                            let link_params = crate::link_hsps::LinkHSPParameters {
                                gap_prob: crate::stat::BLAST_GAP_PROB_GAPPED,
                                gap_decay_rate: crate::stat::BLAST_GAP_DECAY_RATE_GAPPED,
                                longest_intron: translated_gapped_link_longest_intron(
                                    params.max_intron_length,
                                ),
                                cutoff_small_gap: p.seed_cutoff,
                                ..crate::link_hsps::LinkHSPParameters::default()
                            };
                            let link_context = crate::blast_kappa::HitlistLinkContext {
                                query_info: &p.query_info,
                                query_context: 0,
                                score_block: &score_block,
                                link_params: &link_params,
                                gapped_calculation: true,
                            };
                            let scoring_options = crate::options::ScoringOptions {
                                matrix_path: None,
                                reward: 0,
                                penalty: 0,
                                gap_open: params.gap_open,
                                gap_extend: params.gap_extend,
                                shift_pen: i16::MAX as i32,
                                gapped_calculation: true,
                                complexity_adjusted_scoring: false,
                                matrix_name: Some(matrix_name.to_string()),
                                is_ooframe: false,
                                program_number: crate::program::UNDEFINED,
                            };
                            let mut scoring = crate::parameters::ScoringParameters::from_options(
                                &scoring_options,
                                1.0,
                            );
                            let mut matrix_vec: Vec<Vec<i32>> =
                                matrix.iter().map(|row| row.to_vec()).collect();
                            let mut kbp_gap =
                                vec![prot_kbp.clone(); p.query_info.contexts.len().max(1)];
                            let mut saved =
                                crate::blast_kappa::BlastKappaSavedParameters::default();
                            let mut results = crate::hspstream::HspResults::new(
                                p.query_info.contexts.len().max(1) as i32,
                            );
                            let subject = crate::blast_kappa::BlastRedoInMemorySubject {
                                subject_source: &subject_ncbi4na,
                                reward: 0,
                                penalty: 0,
                                genetic_code: db_code,
                                smith_waterman: false,
                                expect_value: params.evalue_threshold,
                                hitlist_size: max_hsps
                                    .unwrap_or(crate::options::HITLIST_SIZE as usize)
                                    as i32,
                                inclusion_ethresh: f64::INFINITY,
                                link_context: Some(&link_context),
                            };
                            let mut legacy_this_match =
                                crate::hspstream::HspList::new(per_query_internal_hsps[qi].oid);
                            let status = crate::blast_kappa::blast_redo_alignment_core_mt(
                                crate::program::TBLASTN,
                                1,
                                &p.aa,
                                &p.query_info,
                                &mut kbp_gap,
                                &mut matrix_vec,
                                &mut scoring,
                                &align_params,
                                &mut saved,
                                &mut legacy_this_match,
                                crate::blast_kappa::BlastRedoAlignmentSource::InMemorySubjectHspList {
                                    hsp_list: &mut per_query_internal_hsps[qi],
                                    subject,
                                },
                                &mut results,
                            );
                            if status != 0 {
                                break 'redo Vec::new();
                            }
                            let mut hsp_list = std::mem::replace(
                                &mut per_query_internal_hsps[qi],
                                blast_hsp_list_new(0),
                            );
                            hsp_list.hsp_array.retain(|hsp| {
                                hsp.as_ref()
                                    .is_some_and(|hsp| hsp.evalue <= params.evalue_threshold)
                            });
                            hsp_list.hspcnt = hsp_list.hsp_array.len() as i32;
                            hsp_list
                                .hsp_array
                                .into_iter()
                                .flatten()
                                .map(|hsp| tblastn_blast_hsp_to_api_hsp(&hsp, &prot_kbp))
                                .collect()
                        };
                        for api_hsp in api_hsps {
                            match &mut per_query_result[qi] {
                                Some(existing) => existing.hsps.push(api_hsp),
                                None => {
                                    let acc = accession.get_or_insert_with(|| {
                                        db.get_accession(oid)
                                            .unwrap_or_else(|| format!("oid_{}", oid))
                                    });
                                    let ttl = title.get_or_insert_with(|| {
                                        String::from_utf8_lossy(db.get_header(oid)).to_string()
                                    });
                                    per_query_result[qi] = Some(SearchResult {
                                        subject_oid: oid,
                                        subject_title: ttl.clone(),
                                        subject_accession: acc.clone(),
                                        subject_len,
                                        hsps: vec![api_hsp],
                                        taxids: vec![],
                                    });
                                }
                            }
                        }
                        continue;
                    }
                    per_query_stat_length[qi] = tblastn_traceback_subject_hsp_list(
                        &p.aa,
                        &translation_buffer,
                        &frame_offsets,
                        &subject_ncbi4na,
                        db_code,
                        &matrix,
                        &mut per_query_internal_hsps[qi],
                        params.gap_open,
                        params.gap_extend,
                        x_drop_final,
                        scratch,
                        qi as i32,
                        max_subj_prot_len,
                    );
                    if translated_sum_stats && params.comp_adjust == 0 {
                        for api_hsp in tblastn_link_internal_hsps_to_api_hsps(
                            &mut per_query_internal_hsps[qi],
                            &p.query_info,
                            0,
                            per_query_stat_length[qi],
                            &prot_kbp,
                            p.gumbel_blk.as_ref(),
                            p.link_gbp_db_length,
                            p.seed_cutoff,
                            params.max_intron_length,
                            params.evalue_threshold,
                        ) {
                            match &mut per_query_result[qi] {
                                Some(existing) => existing.hsps.push(api_hsp),
                                None => {
                                    let acc = accession.get_or_insert_with(|| {
                                        db.get_accession(oid)
                                            .unwrap_or_else(|| format!("oid_{}", oid))
                                    });
                                    let ttl = title.get_or_insert_with(|| {
                                        String::from_utf8_lossy(db.get_header(oid)).to_string()
                                    });
                                    per_query_result[qi] = Some(SearchResult {
                                        subject_oid: oid,
                                        subject_title: ttl.clone(),
                                        subject_accession: acc.clone(),
                                        subject_len,
                                        hsps: vec![api_hsp],
                                        taxids: vec![],
                                    });
                                }
                            }
                        }
                        continue;
                    }
                    per_query_internal_hsps[qi].hsp_array.sort_by(|a, b| {
                        let Some(a) = a.as_ref() else {
                            return std::cmp::Ordering::Greater;
                        };
                        let Some(b) = b.as_ref() else {
                            return std::cmp::Ordering::Less;
                        };
                        (a.subject.frame as i32)
                            .cmp(&(b.subject.frame as i32))
                            .then_with(|| compare_blast_hsps_by_score(a, b))
                    });
                    let items = std::mem::take(&mut per_query_internal_hsps[qi].hsp_array);
                    let mut group_start = 0usize;
                    while group_start < items.len() {
                        let Some(first) = items[group_start].as_ref() else {
                            group_start += 1;
                            continue;
                        };
                        let frame = first.subject.frame as i32;
                        let Some(ctx) = translated_frame_to_context(frame) else {
                            group_start += 1;
                            continue;
                        };
                        let prot_range = (frame_offsets[ctx] + 1).max(0) as usize
                            ..frame_offsets[ctx + 1].max(0) as usize;
                        let subj_prot_len = prot_range.len();
                        let mut group_end = group_start + 1;
                        while group_end < items.len()
                            && items[group_end]
                                .as_ref()
                                .is_some_and(|hsp| hsp.subject.frame as i32 == frame)
                        {
                            group_end += 1;
                        }
                        let prot = &translation_buffer[prot_range.clone()];
                        let phits: Vec<_> = items[group_start..group_end]
                            .iter()
                            .filter_map(|item| item.as_ref().map(blast_hsp_to_protein_hit))
                            .collect();
                        for hsp in process_tblastn_frame_hsps(
                            phits,
                            &p.aa,
                            prot,
                            Some(&subject_ncbi4na),
                            frame,
                            subject_len,
                            subj_prot_len,
                            (subject_len / crate::util::CODON_LENGTH).max(1),
                            &p.query_info,
                            0,
                            p.search_space,
                            &matrix,
                            &prot_kbp,
                            &p.gumbel_blk,
                            x_drop_final,
                            p.hit_cutoff,
                            p.redo_cutoff,
                            p.seed_cutoff,
                            translated_sum_stats,
                            params,
                            scratch,
                        ) {
                            match &mut per_query_result[qi] {
                                Some(existing) => existing.hsps.push(hsp),
                                None => {
                                    let acc = accession.get_or_insert_with(|| {
                                        db.get_accession(oid)
                                            .unwrap_or_else(|| format!("oid_{}", oid))
                                    });
                                    let ttl = title.get_or_insert_with(|| {
                                        String::from_utf8_lossy(db.get_header(oid)).to_string()
                                    });
                                    per_query_result[qi] = Some(SearchResult {
                                        subject_oid: oid,
                                        subject_title: ttl.clone(),
                                        subject_accession: acc.clone(),
                                        subject_len,
                                        hsps: vec![hsp],
                                        taxids: vec![],
                                    });
                                }
                            }
                        }
                        group_start = group_end;
                    }
                }
            }
            for (qi, r) in per_query_result.into_iter().enumerate() {
                if let Some(sr) = r {
                    out.push((qi, sr, per_query_stat_length[qi]));
                }
            }
            out
        },
    );

    let mut results: Vec<Vec<SearchResult>> = (0..queries.len()).map(|_| Vec::new()).collect();
    let mut subject_stat_lengths: Vec<Vec<i32>> = (0..queries.len()).map(|_| Vec::new()).collect();
    for oid_hits in per_oid {
        for (qi, sr, stat_length) in oid_hits {
            results[qi].push(sr);
            subject_stat_lengths[qi].push(stat_length);
        }
    }
    (0..queries.len())
        .map(|qi| {
            let mut qres = std::mem::take(&mut results[qi]);
            if translated_sum_stats && !(params.comp_adjust == 0 && !params.ungapped) {
                apply_tblastn_linked_sum_stats(
                    &mut qres,
                    Some(&subject_stat_lengths[qi]),
                    &qplans[qi].query_info,
                    &prot_kbp,
                    qplans[qi].gumbel_blk.as_ref(),
                    stats_db_len,
                    params.comp_adjust == 0,
                    qplans[qi].seed_cutoff,
                    params.gap_open,
                    params.gap_extend,
                    params.max_intron_length,
                );
            }
            reap_hsps_by_prelim_evalue(&mut qres, prelim_evalue);
            apply_api_evalue_filter(&mut qres, params.evalue_threshold);
            for result in &mut qres {
                result.hsps.sort_by(compare_hsps_by_evalue_then_score);
            }
            apply_api_min_score_filter(&mut qres, params.min_score);
            if let Some(culling_limit) = params.culling_limit {
                apply_api_culling_limit(&mut qres, culling_limit, crate::program::TBLASTN);
            }
            apply_api_max_hsps_limit(&mut qres, max_hsps);
            qres.sort_by(compare_search_results);
            if qres.len() > params.max_target_seqs {
                qres.truncate(params.max_target_seqs);
            }
            fill_protein_midlines(&mut qres, ProteinMidlineStyle::AnyMismatchPlus);
            qres
        })
        .collect()
}

/// Convert one (query-frame, subject-frame) batch of ungapped tblastx hits into
/// internal `BlastHSP`-shaped records. Public `Hsp` conversion happens only
/// after the subject-level HSP list has been reevaluated/scored/reaped.
/// blast-rs: refactor-extracted helper around ported lower-level pieces.
#[allow(clippy::too_many_arguments)]
fn process_tblastx_pair_internal_hsps(
    ungapped_hits: Vec<crate::protein_lookup::ProteinHit>,
    q_prot: &[u8],
    q_context: usize,
    qframe: i32,
    sframe: i32,
    target_t: &mut Option<crate::util::SBlastTargetTranslation>,
    save_cutoff: i32,
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
) -> Vec<crate::hspstream::Hsp> {
    let mut out = Vec::new();
    for mut ph in ungapped_hits {
        if ph.score < save_cutoff {
            continue;
        }
        let mut target_hsp = crate::hspstream::Hsp {
            score: ph.score,
            num_ident: ph.num_ident,
            bit_score: 0.0,
            evalue: 0.0,
            query_offset: ph.query_start as i32,
            query_end: ph.query_end as i32,
            query_gapped_start: ph.gapped_start_q as i32,
            subject_offset: ph.subject_start as i32,
            subject_end: ph.subject_end as i32,
            subject_gapped_start: ph.gapped_start_s as i32,
            context: 0,
            query_frame: qframe,
            subject_frame: sframe,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        };
        let target_view = target_t.as_mut().and_then(|target| {
            crate::hspstream::blast_hsp_get_target_translation(target, Some(&target_hsp))
        });
        let Some(target_view) = target_view else {
            continue;
        };
        if !reevaluate_ungapped_translated_hsp_target(
            &mut ph,
            q_prot,
            &target_view,
            matrix,
            save_cutoff,
        ) {
            continue;
        }
        target_hsp.subject_offset = ph.subject_start as i32;
        target_hsp.subject_end = ph.subject_end as i32;
        target_hsp.subject_gapped_start = ph.gapped_start_s as i32;
        out.push(crate::hspstream::Hsp {
            score: ph.score,
            num_ident: ph.num_ident,
            bit_score: 0.0,
            evalue: 0.0,
            query_offset: ph.query_start as i32,
            query_end: ph.query_end as i32,
            query_gapped_start: ph.gapped_start_q as i32,
            subject_offset: ph.subject_start as i32,
            subject_end: ph.subject_end as i32,
            subject_gapped_start: ph.gapped_start_s as i32,
            context: q_context as i32,
            query_frame: qframe,
            subject_frame: sframe,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
    }
    out
}

fn finalize_tblastx_internal_hsp_list(
    hsp_list: &mut crate::hspstream::HspList,
    query_info: &crate::queryinfo::QueryInfo,
    kbp_array: &[KarlinBlk],
    subject_len: usize,
    evalue_threshold: f64,
    sum_stats: bool,
) {
    if hsp_list.hsps.is_empty() {
        return;
    }
    let mut sbp = crate::stat::blast_score_blk_new(
        crate::encoding::BLASTAA_SEQ_CODE,
        query_info.contexts.len() as i32,
    )
    .expect("protein score block");
    sbp.kbp = kbp_array.to_vec();
    sbp.gbp = None;
    let translated_subject_len = (subject_len / crate::util::CODON_LENGTH).max(1);
    let _ = crate::blast_kappa::blast_hsp_list_get_evalues(
        crate::program::TBLASTX,
        query_info,
        translated_subject_len as i32,
        hsp_list,
        false,
        false,
        &sbp,
        0.0,
        1.0,
    );
    for hsp in &mut hsp_list.hsps {
        let context = hsp.context.max(0) as usize;
        if let Some(kbp) = kbp_array.get(context).or_else(|| kbp_array.first()) {
            hsp.bit_score = kbp.raw_to_bit(hsp.score);
        }
    }
    if !sum_stats {
        crate::hits::filter_by_evalue(hsp_list, evalue_threshold);
    }
    hsp_list.sort_by_score();
}

#[allow(clippy::too_many_arguments)]
fn tblastx_internal_hsp_to_api_hsp(
    hsp: &crate::hspstream::Hsp,
    q_prot: &[u8],
    q_prot_unmasked: &[u8],
    target_t: &mut Option<crate::util::SBlastTargetTranslation>,
) -> Option<Hsp> {
    let target_view = target_t
        .as_mut()
        .and_then(|target| crate::hspstream::blast_hsp_get_target_translation(target, Some(hsp)))?;
    let q_start = hsp.query_offset.max(0) as usize;
    let q_end = hsp.query_end.max(hsp.query_offset) as usize;
    let s_start = hsp.subject_offset.max(0) as usize;
    let s_end = hsp.subject_end.max(hsp.subject_offset) as usize;
    if q_end > q_prot_unmasked.len() || q_start > q_end || s_start > s_end {
        return None;
    }

    let subject_residues: Option<Vec<u8>> = (s_start..s_end)
        .map(|offset| target_view.get(offset as i32))
        .collect();
    let subject_residues = subject_residues?;
    let recount_ident: i32 = q_prot_unmasked[q_start..q_end]
        .iter()
        .zip(subject_residues.iter())
        .filter(|(q, s)| q == s)
        .count() as i32;
    let mut q_aln = ncbistdaa_to_aminoacid_sequence(&q_prot_unmasked[q_start..q_end]);
    let s_aln = ncbistdaa_to_aminoacid_sequence(&subject_residues);
    for (i, idx) in (q_start..q_end).enumerate() {
        if i < q_aln.len()
            && q_prot[idx] == crate::encoding::NCBISTDAA_X
            && q_prot_unmasked[idx] != crate::encoding::NCBISTDAA_X
        {
            q_aln[i] = q_aln[i].to_ascii_lowercase();
        }
    }

    let (query_start, query_end) =
        crate::util::protein_to_oriented_nuc_coords(q_start, q_end, hsp.query_frame);
    let (subject_start, subject_end) =
        crate::util::protein_to_oriented_nuc_coords(s_start, s_end, hsp.subject_frame);
    let (query_gapped_start, _) = crate::util::protein_to_oriented_nuc_coords(
        hsp.query_gapped_start.max(0) as usize,
        hsp.query_gapped_start.max(0) as usize + 1,
        hsp.query_frame,
    );
    let (subject_gapped_start, _) = crate::util::protein_to_oriented_nuc_coords(
        hsp.subject_gapped_start.max(0) as usize,
        hsp.subject_gapped_start.max(0) as usize + 1,
        hsp.subject_frame,
    );

    Some(Hsp {
        score: hsp.score,
        bit_score: hsp.bit_score,
        evalue: hsp.evalue,
        query_start,
        query_end,
        subject_start,
        subject_end,
        query_gapped_start,
        subject_gapped_start,
        query_link_start: q_start,
        query_link_end: q_end,
        query_link_gapped_start: hsp.query_gapped_start.max(0) as usize,
        subject_link_start: s_start,
        subject_link_end: s_end,
        subject_link_gapped_start: hsp.subject_gapped_start.max(0) as usize,
        link_score: None,
        link_lambda: None,
        num_identities: recount_ident as usize,
        num_gaps: hsp.num_gaps.max(0) as usize,
        alignment_length: q_end.saturating_sub(q_start),
        query_aln: q_aln,
        midline: Vec::new(),
        subject_aln: s_aln,
        query_frame: hsp.query_frame,
        subject_frame: hsp.subject_frame,
        num_links: 1,
        comp_adjust_method: hsp.comp_adjustment_method as u8,
    })
}

/// Port of NCBI `Blast_HSPReevaluateWithAmbiguitiesUngapped` for translated
/// ungapped HSPs (`blast_hits.c:676`). NCBI runs this for ungapped searches
/// against nucleotide subjects after initial collection, including tblastx.
#[cfg(test)]
fn reevaluate_ungapped_translated_hsp(
    ph: &mut crate::protein_lookup::ProteinHit,
    query: &[u8],
    subject: &[u8],
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
    cutoff_score: i32,
) -> bool {
    let hsp_len = ph.query_end.saturating_sub(ph.query_start);
    if hsp_len == 0 || ph.subject_end.saturating_sub(ph.subject_start) != hsp_len {
        return false;
    }

    let mut q = ph.query_start;
    let mut s = ph.subject_start;
    let mut score = 0i32;
    let mut sum = 0i32;
    let mut current_q_start = q;
    let mut current_s_start = s;
    let mut best_q_start = q;
    let mut best_q_end = q;
    let mut best_s_start = s;
    let mut best_s_end = s;

    for _ in 0..hsp_len {
        let Some(&q_res) = query.get(q) else {
            return false;
        };
        let Some(&s_res) = subject.get(s) else {
            return false;
        };
        if q_res as usize >= AA_SIZE || s_res as usize >= AA_SIZE {
            return false;
        }

        sum += matrix[q_res as usize][s_res as usize];
        q += 1;
        s += 1;

        if sum < 0 {
            sum = 0;
            current_q_start = q;
            current_s_start = s;
            if score < cutoff_score {
                best_q_start = q;
                best_q_end = q;
                best_s_start = s;
                best_s_end = s;
                score = 0;
            }
        } else if sum > score {
            score = sum;
            best_q_start = current_q_start;
            best_q_end = q;
            best_s_start = current_s_start;
            best_s_end = s;
        }
    }

    if s_update_reevaluated_hsp_ungapped(
        ph,
        cutoff_score,
        score,
        best_q_start,
        best_q_end,
        best_s_start,
        best_s_end,
    ) {
        return false;
    }

    ph.num_ident = query[ph.query_start..ph.query_end]
        .iter()
        .zip(subject[ph.subject_start..ph.subject_end].iter())
        .filter(|(q, s)| q == s)
        .count() as i32;
    ph.mismatches = ph.align_length - ph.num_ident;
    true
}

fn reevaluate_ungapped_translated_hsp_target(
    ph: &mut crate::protein_lookup::ProteinHit,
    query: &[u8],
    subject: &crate::hspstream::BlastTargetTranslationView<'_>,
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
    cutoff_score: i32,
) -> bool {
    let hsp_len = ph.query_end.saturating_sub(ph.query_start);
    if hsp_len == 0 || ph.subject_end.saturating_sub(ph.subject_start) != hsp_len {
        return false;
    }

    let mut q = ph.query_start;
    let mut s = ph.subject_start;
    let mut score = 0i32;
    let mut sum = 0i32;
    let mut current_q_start = q;
    let mut current_s_start = s;
    let mut best_q_start = q;
    let mut best_q_end = q;
    let mut best_s_start = s;
    let mut best_s_end = s;

    for _ in 0..hsp_len {
        let Some(&q_res) = query.get(q) else {
            return false;
        };
        let Some(s_res) = subject.get(s as i32) else {
            return false;
        };
        if q_res as usize >= AA_SIZE || s_res as usize >= AA_SIZE {
            return false;
        }

        sum += matrix[q_res as usize][s_res as usize];
        q += 1;
        s += 1;

        if sum < 0 {
            sum = 0;
            current_q_start = q;
            current_s_start = s;
            if score < cutoff_score {
                best_q_start = q;
                best_q_end = q;
                best_s_start = s;
                best_s_end = s;
                score = 0;
            }
        } else if sum > score {
            score = sum;
            best_q_start = current_q_start;
            best_q_end = q;
            best_s_start = current_s_start;
            best_s_end = s;
        }
    }

    if s_update_reevaluated_hsp_ungapped(
        ph,
        cutoff_score,
        score,
        best_q_start,
        best_q_end,
        best_s_start,
        best_s_end,
    ) {
        return false;
    }

    ph.num_ident = (ph.query_start..ph.query_end)
        .zip(ph.subject_start..ph.subject_end)
        .filter(|&(qi, si)| query.get(qi).copied() == subject.get(si as i32))
        .count() as i32;
    ph.mismatches = ph.align_length - ph.num_ident;
    true
}

/// Run a tblastx search (translated nt query vs translated nt database).
/// blast-rs: Public tblastx API pipeline assembled from ported lower-level pieces;
/// not a direct NCBI C port.
pub fn tblastx(db: &BlastDb, query: &[u8], params: &SearchParams) -> Vec<SearchResult> {
    if query.len() < 3 {
        return Vec::new();
    }
    let query_ncbi4na = encode_ncbi4na_sequence(query);
    let q_code = crate::util::lookup_genetic_code(params.query_gencode);
    let (query_translation, query_offsets) =
        crate::util::blast_get_all_translations(&query_ncbi4na, query_ncbi4na.len(), q_code);
    // NCBI tblastx applies SEG to the translated query: extension uses the
    // masked bytes (X has slightly-negative matrix entries which contain
    // X-drop in low-complexity), but identity-counting uses the ORIGINAL
    // unmasked residues (NCBI's pairwise output shows `i*k*iik` lowercase,
    // and `s_Blast_HSPGetNumIdentitiesAndPositives` reads the query buffer
    // which was NOT modified by NCBI's soft-mask path). We mirror that
    // by keeping the unmasked copy for post-extension identity recount and
    // mutating the working buffer with X for the scan path.
    let query_translation_unmasked = query_translation.clone();
    let mut query_translation = query_translation;
    if params.filter_low_complexity {
        for ctx in 0..crate::util::NUM_FRAMES {
            let begin = (query_offsets[ctx] + 1) as usize;
            let end = query_offsets[ctx + 1] as usize;
            if begin < end {
                apply_seg_ncbistdaa_with_options(
                    &mut query_translation[begin..end],
                    params.seg_window,
                    params.seg_locut,
                    params.seg_hicut,
                );
            }
        }
    }
    let matrix = *get_matrix(params.matrix);
    let matrix_name = protein_matrix_name(params.matrix);
    let word_size = params.word_size.clamp(2, 7);
    let threshold = params.word_threshold.unwrap_or_else(|| {
        suggested_word_threshold(params.matrix, crate::program::TBLASTX, params.word_size)
    });

    // H2: x_drop_ungapped is computed per query frame inside the scan loop from
    // each context's ungapped Karlin block; `ungapped_kbp` here is only the
    // default fill for `kbp_array` before the per-context values are computed.
    let ungapped_kbp = crate::stat::protein_ungapped_kbp_for_matrix(matrix_name);
    let total_subj_len: usize = (0..db.num_oids)
        .map(|oid| db.get_seq_len(oid) as usize)
        .sum();
    let stats_db_len = statistical_db_length(params, total_subj_len as i64);
    let translated_total_subj_len = (stats_db_len as usize / 3).max(1);

    // BLAST_CalcEffLengths once for all 6 query contexts; produces
    // per-context eff_searchsp + length_adjustment that the search loop
    // reads directly. Mirrors NCBI's setup-time invocation in
    // blast_engine.c. db_length /= 3 happens inside CalcEffLengths because
    // tblastx's subject_is_translated is true.
    let mut query_info_calc =
        crate::queryinfo::QueryInfo::new_translated_query_from_offsets(&query_offsets);
    struct TblastxQueryPrecalc<'a> {
        ctx: usize,
        frame: i32,
        prot: &'a [u8],
        prot_unmasked: &'a [u8],
        kbp: KarlinBlk,
        gap_trigger_raw: i32,
    }

    let mut query_precalc = Vec::new();
    let mut kbp_array = vec![ungapped_kbp.clone(); crate::util::NUM_FRAMES];
    for q_ctx in 0..crate::util::NUM_FRAMES {
        let qframe = crate::util::blast_context_to_frame(q_ctx as u32);
        let q_begin = (query_offsets[q_ctx] + 1) as usize;
        let q_end = query_offsets[q_ctx + 1] as usize;
        if q_begin >= q_end {
            continue;
        }
        let q_prot: &[u8] = &query_translation[q_begin..q_end];
        let q_prot_unmasked: &[u8] = &query_translation_unmasked[q_begin..q_end];
        if q_prot.len() < word_size {
            continue;
        }
        if crate::composition::blast_read_aa_composition(q_prot, AA_SIZE).1 == 0 {
            continue;
        }

        // Per-context ungapped Karlin params from this query frame's
        // amino-acid composition. For translated queries (blastx/tblastx/
        // rps-tblastn) NCBI also applies an "ideal-Lambda cap"
        // (`blast_stat.c:2796`) that we mirror here.
        let ideal = crate::stat::protein_ideal_ungapped_kbp_for_matrix(matrix_name);
        let mut ctx_kbp = crate::stat::query_specific_protein_ungapped_kbp_for_matrix(
            q_prot,
            matrix_name,
            &matrix,
        );
        if ctx_kbp.lambda >= ideal.lambda {
            ctx_kbp = ideal;
        }
        let ctx_gap_trigger_raw =
            raw_gap_trigger_bits(crate::stat::BLAST_GAP_TRIGGER_PROT, &ctx_kbp);
        kbp_array[q_ctx] = ctx_kbp.clone();
        query_precalc.push(TblastxQueryPrecalc {
            ctx: q_ctx,
            frame: qframe,
            prot: q_prot,
            prot_unmasked: q_prot_unmasked,
            kbp: ctx_kbp,
            gap_trigger_raw: ctx_gap_trigger_raw,
        });
    }

    // tblastx is ungapped-only (`blast_options.c:869`). NCBI propagates
    // `gapped_calculation = FALSE` into `BLAST_CalcEffLengths`, which routes
    // alpha/beta through the ungapped lookup (`blast_setup.c:806`). Pass the
    // query-specific `sbp->kbp` equivalent so length adjustment and raw HSP
    // e-values use the same per-frame Karlin block.
    let scoring_options = crate::options::ScoringOptions {
        matrix_path: None,
        reward: 0,
        penalty: 0,
        gap_open: params.gap_open,
        gap_extend: params.gap_extend,
        shift_pen: i16::MAX as i32,
        gapped_calculation: false,
        complexity_adjusted_scoring: false,
        matrix_name: Some(matrix_name.to_string()),
        is_ooframe: false,
        program_number: crate::program::UNDEFINED,
    };
    let eff_params = crate::parameters::EffectiveLengthsParameters {
        options: effective_lengths_options(params),
        real_db_length: stats_db_len,
        real_num_seqs: db.num_oids as i32,
    };
    // H2: tblastx is ungapped-only; the eff-lengths std array uses the
    // per-context ungapped Karlin blocks (`Blast_ScoreBlkKbpUngappedCalc`),
    // matching `kbp_array` (also per-context here), not a single rounded row.
    let kbp_std_array = kbp_array.clone();
    crate::blast_setup::blast_calc_eff_lengths(
        crate::program::TBLASTX,
        &scoring_options,
        &eff_params,
        &kbp_array,
        &kbp_std_array,
        matrix_name,
        &mut query_info_calc,
    );

    // tblastx is ungapped-only; NCBI's `Blast_ScoreBlkMatrixInit` frees
    // `sbp->gbp` when `gapped_calculation` is FALSE (`blast_setup.c:524`),
    // so e-values fall back to `BLAST_KarlinStoE_simple`, NOT Spouge FSC.
    let gumbel_blk: Option<crate::stat::GumbelBlk> = None;
    let _ = translated_total_subj_len;
    let max_hsps = params.max_hsps;

    struct TblastxQueryFrame<'a> {
        ctx: usize,
        query_offset: usize,
        frame: i32,
        prot: &'a [u8],
        /// Original (unmasked) translation for this frame. Used for
        /// post-extension identity counting to mirror NCBI's
        /// `s_Blast_HSPGetNumIdentitiesAndPositives` which reads
        /// `query_blk->sequence` (untouched by soft-mask).
        prot_unmasked: &'a [u8],
        search_space: f64,
        kbp: KarlinBlk,
        gap_trigger_raw: i32,
    }

    let mut query_contexts = Vec::new();
    let mut query_frames = Vec::new();
    for q_plan in query_precalc {
        let len_adj = query_info_calc.contexts[q_plan.ctx].length_adjustment;
        let search_space = query_info_calc.contexts[q_plan.ctx].eff_searchsp.max(1) as f64;

        query_contexts.push(TranslatedContextStats {
            context_id: q_plan.ctx,
            query_offset: query_info_calc.contexts[q_plan.ctx].query_offset,
            frame: q_plan.frame,
            query_length: q_plan.prot.len() as i32,
            eff_searchsp: search_space.max(1.0) as i64,
            length_adjustment: len_adj,
            kbp: q_plan.kbp.clone(),
        });
        query_frames.push(TblastxQueryFrame {
            ctx: q_plan.ctx,
            query_offset: query_info_calc.contexts[q_plan.ctx].query_offset as usize,
            frame: q_plan.frame,
            prot: q_plan.prot,
            prot_unmasked: q_plan.prot_unmasked,
            search_space,
            kbp: q_plan.kbp,
            gap_trigger_raw: q_plan.gap_trigger_raw,
        });
    }

    if query_frames.is_empty() {
        return Vec::new();
    }

    let mut concat: Vec<u8> =
        vec![crate::protein_lookup::PROTEIN_CONCAT_SENTINEL; query_info_calc.seq_buf_len().max(1)];
    for q_plan in &query_frames {
        let off = q_plan.query_offset;
        concat[off..off + q_plan.prot.len()].copy_from_slice(q_plan.prot);
    }
    let lookup =
        crate::protein_lookup::ProteinLookupTable::build(&concat, word_size, &matrix, threshold);

    let mut results = vec![None; db.num_oids as usize];
    let avg_subj_len_nt = if db.num_oids > 0 {
        (total_subj_len as i64 / db.num_oids as i64).max(1) as i32
    } else {
        1
    };
    let cutoff_score_min = query_frames
        .iter()
        .map(|q_plan| {
            let hit_cutoff = protein_eval_cutoff(
                params.evalue_threshold,
                &q_plan.kbp,
                gumbel_blk.as_ref(),
                q_plan.prot.len() as i32,
                1,
                q_plan.search_space,
            );
            tblastx_initial_word_cutoff(
                q_plan.gap_trigger_raw,
                &q_plan.kbp,
                q_plan.prot.len() as i32,
                avg_subj_len_nt,
                hit_cutoff,
            )
        })
        .min()
        .unwrap_or(0);
    // Reuse a single diagonal tracking buffer across all
    // (oid, s_frame) calls of the concatenated-context word finder.
    // The scanner clears + resizes it on each call, so it's safe to share.
    let mut scan_diag_table =
        crate::protein_lookup::BlastExtendWord::new(params.two_hit_window as i32);
    for oid in 0..db.num_oids {
        let subject_packed = db.get_sequence(oid);
        let subject_len = db.get_seq_len(oid) as usize;
        // Defer defline parsing until this subject produces a reported HSP,
        // and build it at most once per OID (not once per frame-pair).
        let mut accession: Option<String> = None;
        let mut title: Option<String> = None;

        let (subj_translation, subj_offsets, subj_ncbi4na) = match db.get_ambiguity_data(oid) {
            Some(amb) => {
                let subj_blastna = crate::search::decode_packed_ncbi2na_with_ambiguity(
                    subject_packed,
                    subject_len,
                    amb,
                );
                let subj_ncbi4na = crate::encoding::blastna_to_ncbi4na_sequence(&subj_blastna);
                let translations = crate::util::blast_get_all_translations(
                    &subj_ncbi4na,
                    subj_ncbi4na.len(),
                    crate::util::lookup_genetic_code(params.db_gencode),
                );
                (translations.0, translations.1, subj_ncbi4na)
            }
            None => {
                let subj_blastna =
                    crate::search::decode_packed_ncbi2na(subject_packed, subject_len);
                let subj_ncbi4na = crate::encoding::blastna_to_ncbi4na_sequence(&subj_blastna);
                let translations = crate::util::blast_get_all_translations_packed(
                    subject_packed,
                    subject_len,
                    crate::util::lookup_genetic_code(params.db_gencode),
                );
                (translations.0, translations.1, subj_ncbi4na)
            }
        };
        let mut subject_blk = crate::util::BlastSequenceBlk {
            sequence: Some(subj_ncbi4na.clone()),
            sequence_start: Some(subj_ncbi4na),
            length: subject_len as i32,
            oid: oid as i32,
            ..Default::default()
        };
        let mut target_t = None;
        crate::util::blast_target_translation_new(
            &mut subject_blk,
            crate::util::lookup_genetic_code(params.db_gencode),
            crate::program::TBLASTX,
            false,
            &mut target_t,
        );
        let mut subject_hsp_list = crate::hspstream::HspList::new(oid as i32);

        for s_ctx in 0..crate::util::NUM_FRAMES {
            let sframe = crate::util::blast_context_to_frame(s_ctx as u32);
            let s_begin = (subj_offsets[s_ctx] + 1) as usize;
            let s_end = subj_offsets[s_ctx + 1] as usize;
            if s_begin >= s_end {
                continue;
            }
            let s_prot: &[u8] = &subj_translation[s_begin..s_end];
            if s_prot.len() < word_size {
                continue;
            }
            if crate::composition::blast_read_aa_composition(s_prot, AA_SIZE).1 == 0 {
                continue;
            }
            let mut context_scan_params = Vec::with_capacity(query_frames.len());
            let mut save_cutoffs = Vec::with_capacity(query_frames.len());
            for q_plan in &query_frames {
                let q_prot = q_plan.prot;
                let search_space = q_plan.search_space;
                let ctx_kbp = &q_plan.kbp;
                let hit_cutoff = protein_eval_cutoff(
                    params.evalue_threshold,
                    ctx_kbp,
                    gumbel_blk.as_ref(),
                    q_prot.len() as i32,
                    1,
                    search_space,
                );
                let save_cutoff = tblastx_initial_word_cutoff(
                    q_plan.gap_trigger_raw,
                    ctx_kbp,
                    q_prot.len() as i32,
                    avg_subj_len_nt,
                    hit_cutoff,
                );
                save_cutoffs.push(save_cutoff);
                context_scan_params.push(crate::protein_lookup::ProteinContextScanParams {
                    query_start: q_plan.query_offset,
                    x_dropoff: raw_ungapped_xdrop_bits(params.x_drop_ungapped, ctx_kbp),
                    cutoff_score: save_cutoff.max(1),
                });
            }
            let scan_hits = crate::protein_lookup::protein_scan_with_extend_word_context_params(
                &concat,
                s_prot,
                &matrix,
                &lookup,
                params.two_hit_window as i32,
                &mut scan_diag_table,
                &context_scan_params,
            );
            scan_diag_table.exit(s_prot.len());
            if scan_hits.is_empty() {
                continue;
            }
            let mut per_ctx: Vec<Vec<crate::protein_lookup::ProteinHit>> =
                (0..query_frames.len()).map(|_| Vec::new()).collect();
            for mut h in scan_hits {
                let ctx =
                    crate::queryinfo::bsearch_context_info(h.query_start as i32, &query_info_calc);
                if ctx < 0 {
                    continue;
                }
                let Some(ci) = query_frames.iter().position(|q| q.ctx == ctx as usize) else {
                    continue;
                };
                let off = query_frames[ci].query_offset;
                h.query_start -= off;
                h.query_end = h.query_end.saturating_sub(off);
                h.gapped_start_q = h.gapped_start_q.saturating_sub(off);
                per_ctx[ci].push(h);
            }
            for ci in 0..query_frames.len() {
                let ungapped_hits = std::mem::take(&mut per_ctx[ci]);
                if ungapped_hits.is_empty() {
                    continue;
                }
                let q_plan = &query_frames[ci];
                let qframe = q_plan.frame;
                let q_prot = q_plan.prot;
                let save_cutoff = save_cutoffs[ci];
                let pair_hsps = process_tblastx_pair_internal_hsps(
                    ungapped_hits,
                    q_prot,
                    q_plan.ctx,
                    qframe,
                    sframe,
                    &mut target_t,
                    save_cutoff,
                    &matrix,
                );
                subject_hsp_list.hsps.extend(pair_hsps);
            }
        }
        finalize_tblastx_internal_hsp_list(
            &mut subject_hsp_list,
            &query_info_calc,
            &kbp_array,
            subject_len,
            params.evalue_threshold,
            params.sum_stats,
        );
        let mut subject_hsps = Vec::new();
        for hsp in &subject_hsp_list.hsps {
            let context = hsp.context.max(0) as usize;
            let Some(q_plan) = query_frames.iter().find(|plan| plan.ctx == context) else {
                continue;
            };
            if let Some(api_hsp) = tblastx_internal_hsp_to_api_hsp(
                hsp,
                q_plan.prot,
                q_plan.prot_unmasked,
                &mut target_t,
            ) {
                subject_hsps.push(api_hsp);
            }
        }
        if !subject_hsps.is_empty() {
            let acc = accession.get_or_insert_with(|| {
                db.get_accession(oid)
                    .unwrap_or_else(|| format!("oid_{}", oid))
            });
            let ttl = title
                .get_or_insert_with(|| String::from_utf8_lossy(db.get_header(oid)).to_string());
            results[oid as usize] = Some(SearchResult {
                subject_oid: oid,
                subject_title: ttl.clone(),
                subject_accession: acc.clone(),
                subject_len,
                hsps: subject_hsps,
                taxids: vec![],
            });
        }
    }

    if params.sum_stats {
        apply_tblastx_linked_sum_stats(
            &mut results,
            &query_contexts,
            0,
            stats_db_len,
            cutoff_score_min,
        );
    }

    let mut results: Vec<SearchResult> = results.into_iter().flatten().collect();
    for result in &mut results {
        result.hsps.sort_by(compare_hsps_by_evalue_then_score);
    }
    // Re-apply the evalue threshold after sum-stats linking: the linker can
    // raise per-HSP evalues above the user threshold (when a multi-HSP
    // combination is required to reach significance, the surviving
    // single-HSP evalues are larger than what passed the pre-link cutoff).
    // NCBI emits only HSPs whose post-link evalue is still <= threshold.
    let evalue_cutoff = params.evalue_threshold;
    for result in &mut results {
        result.hsps.retain(|hsp| hsp.evalue <= evalue_cutoff);
    }
    results.retain(|r| !r.hsps.is_empty());
    apply_api_min_score_filter(&mut results, params.min_score);
    if let Some(culling_limit) = params.culling_limit {
        apply_api_culling_limit(&mut results, culling_limit, crate::program::TBLASTX);
    }
    apply_api_max_hsps_limit(&mut results, max_hsps);
    results.sort_by(compare_search_results);
    if results.len() > params.max_target_seqs {
        results.truncate(params.max_target_seqs);
    }
    fill_protein_midlines(&mut results, ProteinMidlineStyle::AnyMismatchPlus);
    results
}

/// Concatenated-query tblastx: translate every query into its 6 frames, SEG +
/// per-frame Karlin params exactly as the single-query `tblastx`, then merge
/// ALL query frames (across all queries) into ONE protein lookup table with
/// inter-context sentinels (NCBI's `BlastQueryInfo` concatenation). Each DB
/// subject is decoded + 6-frame translated ONCE and scanned once per subject
/// frame against the single lookup; hits are attributed back to their
/// (query, frame) via `bsearch_context_info` and the per-context offset. This
/// removes the per-query database rescan (the multi-query slowdown) while
/// staying byte-identical with `tblastx` (the per-hit body is shared via
/// `process_tblastx_pair_hsps`, and the post-search sum-stats + filter tail
/// is the same as the single-query path, run per query).
/// blast-rs: multi-query orchestration over ported pieces; not a direct NCBI C port.
pub fn tblastx_batch(
    db: &BlastDb,
    queries: &[&[u8]],
    params: &SearchParams,
) -> Vec<Vec<SearchResult>> {
    if queries.is_empty() {
        return Vec::new();
    }
    let empty = || (0..queries.len()).map(|_| Vec::new()).collect();

    let matrix = *get_matrix(params.matrix);
    let matrix_name = protein_matrix_name(params.matrix);
    let word_size = params.word_size.clamp(2, 7);
    let threshold = params.word_threshold.unwrap_or_else(|| {
        suggested_word_threshold(params.matrix, crate::program::TBLASTX, params.word_size)
    });
    let ungapped_kbp = crate::stat::protein_ungapped_kbp_for_matrix(matrix_name);
    let total_subj_len: usize = (0..db.num_oids)
        .map(|oid| db.get_seq_len(oid) as usize)
        .sum();
    let stats_db_len = statistical_db_length(params, total_subj_len as i64);
    // tblastx is ungapped-only; NCBI frees `sbp->gbp`, so e-values use the
    // simple Karlin path (no Spouge FSC).
    let gumbel_blk: Option<crate::stat::GumbelBlk> = None;
    let max_hsps = params.max_hsps;
    let q_code = crate::util::lookup_genetic_code(params.query_gencode);
    let db_code = crate::util::lookup_genetic_code(params.db_gencode);

    // tblastx is ungapped-only (`blast_options.c:869`): route alpha/beta through
    // the ungapped lookup in `BLAST_CalcEffLengths` (gapped_calculation=false).
    let scoring_options = crate::options::ScoringOptions {
        matrix_path: None,
        reward: 0,
        penalty: 0,
        gap_open: params.gap_open,
        gap_extend: params.gap_extend,
        shift_pen: i16::MAX as i32,
        gapped_calculation: false,
        complexity_adjusted_scoring: false,
        matrix_name: Some(matrix_name.to_string()),
        is_ooframe: false,
        program_number: crate::program::UNDEFINED,
    };
    let eff_params = crate::parameters::EffectiveLengthsParameters {
        options: effective_lengths_options(params),
        real_db_length: stats_db_len,
        real_num_seqs: db.num_oids as i32,
    };

    // One concatenated context per (query, valid frame). `masked` is the SEG'd
    // working buffer used for the scan/extension; `unmasked` is the original
    // translation used for identity recount + alignment text (NCBI reads the
    // untouched `query_blk->sequence`).
    struct PreCtx {
        query_idx: usize,
        ctx: usize,
        frame: i32,
        global_query_offset: usize,
        masked: Vec<u8>,
        unmasked: Vec<u8>,
        kbp: KarlinBlk,
        gap_trigger_raw: i32,
        search_space: f64,
    }
    let mut pre: Vec<PreCtx> = Vec::new();
    // Per-query frame stats for the post-search sum-stats linker.
    let mut per_query_contexts: Vec<Vec<TranslatedContextStats>> =
        (0..queries.len()).map(|_| Vec::new()).collect();
    let mut per_query_infos: Vec<Option<crate::queryinfo::QueryInfo>> =
        (0..queries.len()).map(|_| None).collect();
    let mut per_query_kbps: Vec<Vec<KarlinBlk>> = (0..queries.len())
        .map(|_| vec![ungapped_kbp.clone(); crate::util::NUM_FRAMES])
        .collect();
    let mut global_query_base = 0usize;

    for (qi, q) in queries.iter().enumerate() {
        if q.len() < 3 {
            continue;
        }
        let query_base = global_query_base;
        let query_ncbi4na = encode_ncbi4na_sequence(q);
        let (query_translation, query_offsets) =
            crate::util::blast_get_all_translations(&query_ncbi4na, query_ncbi4na.len(), q_code);
        let query_translation_unmasked = query_translation.clone();
        let mut query_translation = query_translation;
        if params.filter_low_complexity {
            for ctx in 0..crate::util::NUM_FRAMES {
                let begin = (query_offsets[ctx] + 1) as usize;
                let end = query_offsets[ctx + 1] as usize;
                if begin < end {
                    apply_seg_ncbistdaa_with_options(
                        &mut query_translation[begin..end],
                        params.seg_window,
                        params.seg_locut,
                        params.seg_hicut,
                    );
                }
            }
        }

        let mut query_info_calc =
            crate::queryinfo::QueryInfo::new_translated_query_from_offsets(&query_offsets);
        let mut kbp_array = vec![ungapped_kbp.clone(); crate::util::NUM_FRAMES];
        struct Tmp {
            ctx: usize,
            frame: i32,
            masked: Vec<u8>,
            unmasked: Vec<u8>,
            kbp: KarlinBlk,
            gap_trigger_raw: i32,
        }
        let mut tmps: Vec<Tmp> = Vec::new();
        for q_ctx in 0..crate::util::NUM_FRAMES {
            let qframe = crate::util::blast_context_to_frame(q_ctx as u32);
            let q_begin = (query_offsets[q_ctx] + 1) as usize;
            let q_end = query_offsets[q_ctx + 1] as usize;
            if q_begin >= q_end {
                continue;
            }
            let q_prot: &[u8] = &query_translation[q_begin..q_end];
            let q_prot_unmasked: &[u8] = &query_translation_unmasked[q_begin..q_end];
            if q_prot.len() < word_size {
                continue;
            }
            if crate::composition::blast_read_aa_composition(q_prot, AA_SIZE).1 == 0 {
                continue;
            }
            // Per-frame ungapped Karlin params from this frame's AA composition,
            // with the translated-search "ideal-Lambda cap" (`blast_stat.c:2796`).
            let ideal = crate::stat::protein_ideal_ungapped_kbp_for_matrix(matrix_name);
            let mut ctx_kbp = crate::stat::query_specific_protein_ungapped_kbp_for_matrix(
                q_prot,
                matrix_name,
                &matrix,
            );
            if ctx_kbp.lambda >= ideal.lambda {
                ctx_kbp = ideal;
            }
            let ctx_gap_trigger_raw =
                raw_gap_trigger_bits(crate::stat::BLAST_GAP_TRIGGER_PROT, &ctx_kbp);
            kbp_array[q_ctx] = ctx_kbp.clone();
            tmps.push(Tmp {
                ctx: q_ctx,
                frame: qframe,
                masked: q_prot.to_vec(),
                unmasked: q_prot_unmasked.to_vec(),
                kbp: ctx_kbp,
                gap_trigger_raw: ctx_gap_trigger_raw,
            });
        }
        let kbp_std_array = kbp_array.clone();
        crate::blast_setup::blast_calc_eff_lengths(
            crate::program::TBLASTX,
            &scoring_options,
            &eff_params,
            &kbp_array,
            &kbp_std_array,
            matrix_name,
            &mut query_info_calc,
        );
        per_query_infos[qi] = Some(query_info_calc.clone());
        per_query_kbps[qi] = kbp_array.clone();
        global_query_base += query_info_calc.seq_buf_len().max(1);
        for t in tmps {
            let len_adj = query_info_calc.contexts[t.ctx].length_adjustment;
            let search_space = query_info_calc.contexts[t.ctx].eff_searchsp.max(1) as f64;
            let query_offset = query_info_calc.contexts[t.ctx].query_offset as usize;
            per_query_contexts[qi].push(TranslatedContextStats {
                context_id: t.ctx,
                query_offset: query_info_calc.contexts[t.ctx].query_offset,
                frame: t.frame,
                query_length: t.masked.len() as i32,
                eff_searchsp: search_space.max(1.0) as i64,
                length_adjustment: len_adj,
                kbp: t.kbp.clone(),
            });
            pre.push(PreCtx {
                query_idx: qi,
                ctx: t.ctx,
                frame: t.frame,
                global_query_offset: query_base + query_offset,
                masked: t.masked,
                unmasked: t.unmasked,
                kbp: t.kbp,
                gap_trigger_raw: t.gap_trigger_raw,
                search_space,
            });
        }
    }

    if pre.is_empty() {
        return empty();
    }
    let avg_subj_len_nt = if db.num_oids > 0 {
        (total_subj_len as i64 / db.num_oids as i64).max(1) as i32
    } else {
        1
    };

    // Concatenate all query-frame contexts at their original translated-query
    // offsets, shifted by a per-query base. This preserves the query-offset
    // structure that NCBI's BlastQueryInfo exposes to BlastAaWordFinder.
    let concat_len = pre
        .iter()
        .map(|p| p.global_query_offset + p.masked.len() + 2)
        .max()
        .unwrap_or(1);
    let mut concat: Vec<u8> = vec![crate::protein_lookup::PROTEIN_CONCAT_SENTINEL; concat_len];
    for p in &pre {
        let off = p.global_query_offset;
        concat[off..off + p.masked.len()].copy_from_slice(&p.masked);
    }
    let contexts: Vec<crate::queryinfo::ContextInfo> = pre
        .iter()
        .map(|p| crate::queryinfo::ContextInfo {
            query_offset: p.global_query_offset as i32,
            query_length: p.masked.len() as i32,
            eff_searchsp: p.search_space.max(1.0) as i64,
            length_adjustment: 0,
            query_index: p.query_idx as i32,
            frame: p.frame,
            is_valid: true,
            segment_flags: crate::queryinfo::E_NO_SEGMENTS,
        })
        .collect();
    let max_length = contexts
        .iter()
        .map(|c| c.query_length.max(0) as u32)
        .max()
        .unwrap_or(0);
    let min_length = contexts
        .iter()
        .map(|c| c.query_length.max(0) as u32)
        .min()
        .unwrap_or(0);
    let attrib_info = crate::queryinfo::QueryInfo {
        num_queries: queries.len() as i32,
        contexts,
        max_length,
        min_length,
    };
    let lookup =
        crate::protein_lookup::ProteinLookupTable::build(&concat, word_size, &matrix, threshold);

    // Per OID (parallel): translate subject once, scan each subject frame once
    // against the concatenated lookup, attribute hits, run the shared per-pair
    // body, scatter to the owning query.
    let per_oid: Vec<Vec<(usize, SearchResult)>> = map_database_oids_init(
        db,
        params,
        || crate::protein_lookup::BlastExtendWord::new(params.two_hit_window as i32),
        |scan_diag_table, oid| {
            let mut out: Vec<(usize, SearchResult)> = Vec::new();
            let subject_packed = db.get_sequence(oid);
            let subject_len = db.get_seq_len(oid) as usize;
            let (subj_translation, subj_offsets, subj_ncbi4na) = match db.get_ambiguity_data(oid) {
                Some(amb) => {
                    let subj_blastna = crate::search::decode_packed_ncbi2na_with_ambiguity(
                        subject_packed,
                        subject_len,
                        amb,
                    );
                    let subj_ncbi4na = crate::encoding::blastna_to_ncbi4na_sequence(&subj_blastna);
                    let translations = crate::util::blast_get_all_translations(
                        &subj_ncbi4na,
                        subj_ncbi4na.len(),
                        db_code,
                    );
                    (translations.0, translations.1, subj_ncbi4na)
                }
                None => {
                    let subj_blastna =
                        crate::search::decode_packed_ncbi2na(subject_packed, subject_len);
                    let subj_ncbi4na = crate::encoding::blastna_to_ncbi4na_sequence(&subj_blastna);
                    let translations = crate::util::blast_get_all_translations_packed(
                        subject_packed,
                        subject_len,
                        db_code,
                    );
                    (translations.0, translations.1, subj_ncbi4na)
                }
            };
            let mut subject_blk = crate::util::BlastSequenceBlk {
                sequence: Some(subj_ncbi4na.clone()),
                sequence_start: Some(subj_ncbi4na),
                length: subject_len as i32,
                oid: oid as i32,
                ..Default::default()
            };
            let mut target_t = None;
            crate::util::blast_target_translation_new(
                &mut subject_blk,
                db_code,
                crate::program::TBLASTX,
                false,
                &mut target_t,
            );
            let mut per_query_hsp_lists: Vec<crate::hspstream::HspList> = (0..queries.len())
                .map(|_| crate::hspstream::HspList::new(oid as i32))
                .collect();
            // Defer defline parsing until a subject produces a hit.
            let mut accession: Option<String> = None;
            let mut title: Option<String> = None;

            for s_ctx in 0..crate::util::NUM_FRAMES {
                let sframe = crate::util::blast_context_to_frame(s_ctx as u32);
                let s_begin = (subj_offsets[s_ctx] + 1) as usize;
                let s_end = subj_offsets[s_ctx + 1] as usize;
                if s_begin >= s_end {
                    continue;
                }
                let s_prot: &[u8] = &subj_translation[s_begin..s_end];
                if s_prot.len() < word_size {
                    continue;
                }
                if crate::composition::blast_read_aa_composition(s_prot, AA_SIZE).1 == 0 {
                    continue;
                }
                let context_scan_params: Vec<crate::protein_lookup::ProteinContextScanParams> = pre
                    .iter()
                    .map(|p| {
                        let q_prot_len = p.masked.len() as i32;
                        let hit_cutoff = protein_eval_cutoff(
                            params.evalue_threshold,
                            &p.kbp,
                            gumbel_blk.as_ref(),
                            q_prot_len,
                            1,
                            p.search_space,
                        );
                        let save_cutoff = tblastx_initial_word_cutoff(
                            p.gap_trigger_raw,
                            &p.kbp,
                            q_prot_len,
                            avg_subj_len_nt,
                            hit_cutoff,
                        );
                        crate::protein_lookup::ProteinContextScanParams {
                            query_start: p.global_query_offset,
                            x_dropoff: raw_ungapped_xdrop_bits(params.x_drop_ungapped, &p.kbp),
                            cutoff_score: save_cutoff.max(1),
                        }
                    })
                    .collect();
                let scan_hits = crate::protein_lookup::protein_scan_with_extend_word_context_params(
                    &concat,
                    s_prot,
                    &matrix,
                    &lookup,
                    params.two_hit_window as i32,
                    scan_diag_table,
                    &context_scan_params,
                );
                scan_diag_table.exit(s_prot.len());
                if scan_hits.is_empty() {
                    continue;
                }
                let mut per_ctx: Vec<Vec<crate::protein_lookup::ProteinHit>> =
                    (0..pre.len()).map(|_| Vec::new()).collect();
                for mut h in scan_hits {
                    let ci =
                        crate::queryinfo::bsearch_context_info(h.query_start as i32, &attrib_info);
                    if ci < 0 || ci as usize >= pre.len() {
                        continue;
                    }
                    let ci = ci as usize;
                    let off = pre[ci].global_query_offset;
                    h.query_start -= off;
                    h.query_end = h.query_end.saturating_sub(off);
                    h.gapped_start_q = h.gapped_start_q.saturating_sub(off);
                    per_ctx[ci].push(h);
                }
                for ci in 0..pre.len() {
                    let ungapped = std::mem::take(&mut per_ctx[ci]);
                    if ungapped.is_empty() {
                        continue;
                    }
                    let p = &pre[ci];
                    let q_prot: &[u8] = &p.masked;
                    let ctx_kbp = &p.kbp;
                    // Same cutoffs as single-query tblastx, per (frame, subject-frame).
                    let hit_cutoff = protein_eval_cutoff(
                        params.evalue_threshold,
                        ctx_kbp,
                        gumbel_blk.as_ref(),
                        q_prot.len() as i32,
                        1,
                        p.search_space,
                    );
                    let save_cutoff = tblastx_initial_word_cutoff(
                        p.gap_trigger_raw,
                        ctx_kbp,
                        q_prot.len() as i32,
                        avg_subj_len_nt,
                        hit_cutoff,
                    );
                    let hsps = process_tblastx_pair_internal_hsps(
                        ungapped,
                        q_prot,
                        p.ctx,
                        p.frame,
                        sframe,
                        &mut target_t,
                        save_cutoff,
                        &matrix,
                    );
                    if hsps.is_empty() {
                        continue;
                    }
                    per_query_hsp_lists[p.query_idx].hsps.extend(hsps);
                }
            }
            for (qi, mut hsp_list) in per_query_hsp_lists.into_iter().enumerate() {
                let Some(query_info) = per_query_infos.get(qi).and_then(Option::as_ref) else {
                    continue;
                };
                let Some(query_kbps) = per_query_kbps.get(qi).map(Vec::as_slice) else {
                    continue;
                };
                finalize_tblastx_internal_hsp_list(
                    &mut hsp_list,
                    query_info,
                    query_kbps,
                    subject_len,
                    params.evalue_threshold,
                    params.sum_stats,
                );
                let mut hsps = Vec::new();
                for hsp in &hsp_list.hsps {
                    let context = hsp.context.max(0) as usize;
                    let Some(p) = pre.iter().find(|p| p.query_idx == qi && p.ctx == context) else {
                        continue;
                    };
                    if let Some(api_hsp) =
                        tblastx_internal_hsp_to_api_hsp(hsp, &p.masked, &p.unmasked, &mut target_t)
                    {
                        hsps.push(api_hsp);
                    }
                }
                if !hsps.is_empty() {
                    let acc = accession
                        .get_or_insert_with(|| {
                            db.get_accession(oid)
                                .unwrap_or_else(|| format!("oid_{}", oid))
                        })
                        .clone();
                    let ttl = title
                        .get_or_insert_with(|| {
                            String::from_utf8_lossy(db.get_header(oid)).to_string()
                        })
                        .clone();
                    out.push((
                        qi,
                        SearchResult {
                            subject_oid: oid,
                            subject_title: ttl,
                            subject_accession: acc,
                            subject_len,
                            hsps,
                            taxids: vec![],
                        },
                    ));
                }
            }
            out
        },
    );

    // Gather per-query subject lists, then run the single-query post-search tail
    // (sum-stats linking → re-sort → evalue re-filter → culling → truncate).
    let mut by_query: Vec<Vec<Option<SearchResult>>> =
        (0..queries.len()).map(|_| Vec::new()).collect();
    for oid_hits in per_oid {
        for (qi, sr) in oid_hits {
            by_query[qi].push(Some(sr));
        }
    }
    let avg_subj_len_nt = if db.num_oids > 0 {
        (total_subj_len as i64 / db.num_oids as i64).max(1) as i32
    } else {
        1
    };

    (0..queries.len())
        .map(|qi| {
            let mut opt_results = std::mem::take(&mut by_query[qi]);
            if params.sum_stats {
                // tblastx linker `cutoff_score_min` is the lowest initial-word
                // cutoff (NCBI word-params.cutoff_score_min).
                let cutoff_score_min = pre
                    .iter()
                    .filter(|p| p.query_idx == qi)
                    .map(|p| {
                        let hit_cutoff = protein_eval_cutoff(
                            params.evalue_threshold,
                            &p.kbp,
                            gumbel_blk.as_ref(),
                            p.masked.len() as i32,
                            1,
                            p.search_space,
                        );
                        tblastx_initial_word_cutoff(
                            p.gap_trigger_raw,
                            &p.kbp,
                            p.masked.len() as i32,
                            avg_subj_len_nt,
                            hit_cutoff,
                        )
                    })
                    .min()
                    .unwrap_or(0);
                apply_tblastx_linked_sum_stats(
                    &mut opt_results,
                    &per_query_contexts[qi],
                    0,
                    stats_db_len,
                    cutoff_score_min,
                );
            }
            let mut qres: Vec<SearchResult> = opt_results.into_iter().flatten().collect();
            for result in &mut qres {
                result.hsps.sort_by(compare_hsps_by_evalue_then_score);
            }
            let evalue_cutoff = params.evalue_threshold;
            for result in &mut qres {
                result.hsps.retain(|hsp| hsp.evalue <= evalue_cutoff);
            }
            qres.retain(|r| !r.hsps.is_empty());
            apply_api_min_score_filter(&mut qres, params.min_score);
            if let Some(culling_limit) = params.culling_limit {
                apply_api_culling_limit(&mut qres, culling_limit, crate::program::TBLASTX);
            }
            apply_api_max_hsps_limit(&mut qres, max_hsps);
            qres.sort_by(compare_search_results);
            if qres.len() > params.max_target_seqs {
                qres.truncate(params.max_target_seqs);
            }
            fill_protein_midlines(&mut qres, ProteinMidlineStyle::AnyMismatchPlus);
            qres
        })
        .collect()
}

// ── Utility functions ───────────────────────────────────────────────────────

#[cfg(test)]
fn translated_composition_window_bounds(
    phits: &[crate::protein_lookup::ProteinHit],
    seq_len: usize,
    translated_side_is_query: bool,
    side_is_translated: bool,
) -> (usize, usize) {
    if !side_is_translated || phits.is_empty() {
        return (0, seq_len);
    }

    let border = crate::blast_kappa::K_WINDOW_BORDER.max(0) as usize;
    let mut begin = seq_len;
    let mut end = 0usize;
    for ph in phits {
        let (start, stop) = if translated_side_is_query {
            (ph.query_start, ph.query_end)
        } else {
            (ph.subject_start, ph.subject_end)
        };
        begin = begin.min(start.saturating_sub(border));
        end = end.max(stop.saturating_add(border).min(seq_len));
    }

    if begin >= end {
        (0, seq_len)
    } else {
        (begin, end)
    }
}

/// blast-rs: ProteinHit endpoint dedup helper modeled on the C common-endpoint
/// HSP purge; not a direct NCBI C port. Two-pass dedup:
/// 1. Sort by (query.offset, subject.offset) ASC, score DESC, query.end DESC,
///    subject.end DESC. Within consecutive HSPs sharing same query.offset and
///    same subject.offset: keep the FIRST (highest score), drop the rest.
/// 2. Sort by (query.end, subject.end) ASC, score DESC, query.offset DESC,
///    subject.offset DESC. Same rule on (query.end, subject.end).
///
/// For protein BLAST the engine calls this purge after gapped scoring.
fn purge_hsps_with_common_endpoints(phits: &mut Vec<crate::protein_lookup::ProteinHit>) {
    if phits.len() <= 1 {
        return;
    }
    // Pass 1: sort by query.offset / subject.offset.
    phits.sort_by(|a, b| {
        a.query_start
            .cmp(&b.query_start)
            .then(a.subject_start.cmp(&b.subject_start))
            .then(b.score.cmp(&a.score)) // higher score first on ties
            .then(b.query_end.cmp(&a.query_end))
            .then(b.subject_end.cmp(&a.subject_end))
    });
    let mut keep = vec![true; phits.len()];
    let mut i = 0;
    while i < phits.len() {
        if !keep[i] {
            i += 1;
            continue;
        }
        let mut j = i + 1;
        while j < phits.len()
            && phits[i].query_start == phits[j].query_start
            && phits[i].subject_start == phits[j].subject_start
        {
            keep[j] = false;
            j += 1;
        }
        i = j;
    }
    let mut idx = 0;
    phits.retain(|_| {
        let k = keep[idx];
        idx += 1;
        k
    });
    if phits.len() <= 1 {
        return;
    }
    // Pass 2: sort by query.end / subject.end.
    phits.sort_by(|a, b| {
        a.query_end
            .cmp(&b.query_end)
            .then(a.subject_end.cmp(&b.subject_end))
            .then(b.score.cmp(&a.score))
            .then(b.query_start.cmp(&a.query_start))
            .then(b.subject_start.cmp(&a.subject_start))
    });
    let mut keep = vec![true; phits.len()];
    let mut i = 0;
    while i < phits.len() {
        if !keep[i] {
            i += 1;
            continue;
        }
        let mut j = i + 1;
        while j < phits.len()
            && phits[i].query_end == phits[j].query_end
            && phits[i].subject_end == phits[j].subject_end
        {
            keep[j] = false;
            j += 1;
        }
        i = j;
    }
    let mut idx = 0;
    phits.retain(|_| {
        let k = keep[idx];
        idx += 1;
        k
    });
}

/// blast-rs: Recomputes ProteinHit score from rendered alignment; not a direct
/// NCBI C port.
fn rescore_protein_hit(
    ph: &crate::protein_lookup::ProteinHit,
    query_aa: &[u8],
    subj_aa: &[u8],
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
    gap_open: i32,
    gap_extend: i32,
) -> i32 {
    if let (Some(qs), Some(ss)) = (&ph.qseq, &ph.sseq) {
        let mut score = 0;
        let mut q_idx = ph.query_start;
        let mut s_idx = ph.subject_start;
        let mut in_q_gap = false;
        let mut in_s_gap = false;
        for (qch, sch) in qs.as_bytes().iter().zip(ss.as_bytes().iter()) {
            match (*qch, *sch) {
                (b'-', b'-') => {}
                (b'-', _) => {
                    score -= if in_q_gap {
                        gap_extend
                    } else {
                        gap_open + gap_extend
                    };
                    in_q_gap = true;
                    in_s_gap = false;
                    s_idx += 1;
                }
                (_, b'-') => {
                    score -= if in_s_gap {
                        gap_extend
                    } else {
                        gap_open + gap_extend
                    };
                    in_s_gap = true;
                    in_q_gap = false;
                    q_idx += 1;
                }
                _ => {
                    if q_idx < query_aa.len() && s_idx < subj_aa.len() {
                        score += matrix[query_aa[q_idx] as usize][subj_aa[s_idx] as usize];
                    }
                    in_q_gap = false;
                    in_s_gap = false;
                    q_idx += 1;
                    s_idx += 1;
                }
            }
        }
        score
    } else {
        let q = &query_aa[ph.query_start..ph.query_end];
        let s = &subj_aa[ph.subject_start..ph.subject_end];
        q.iter()
            .zip(s.iter())
            .map(|(&qa, &sa)| matrix[qa as usize][sa as usize])
            .sum()
    }
}

/// Parse a multi-FASTA byte slice into (title, sequence) pairs.
/// blast-rs: Public FASTA parser helper; not a direct NCBI C port.
pub fn parse_fasta(input: &[u8]) -> Vec<(String, Vec<u8>)> {
    let mut sequences = Vec::new();
    let mut current_title = String::new();
    let mut current_seq: Vec<u8> = Vec::new();

    for line in input.split(|&b| b == b'\n') {
        let line = line.strip_suffix(b"\r").unwrap_or(line);
        if line.is_empty() {
            continue;
        }
        if line.starts_with(b">") {
            if !current_title.is_empty() {
                sequences.push((current_title.clone(), current_seq.clone()));
                current_seq.clear();
            }
            current_title = String::from_utf8_lossy(&line[1..]).trim().to_string();
        } else {
            current_seq.extend_from_slice(line);
        }
    }
    if !current_title.is_empty() || !current_seq.is_empty() {
        sequences.push((current_title, current_seq));
    }
    sequences
}

/// Reverse complement an ASCII nucleotide sequence.
/// blast-rs: Public reverse-complement wrapper; not a direct NCBI C port.
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    reverse_complement_iupacna_sequence(seq)
}

/// Six-frame translation of an ASCII nucleotide sequence (public API).
/// blast-rs: Public wrapper over the direct translation helper; not a direct
/// NCBI C port. Unpacks the resulting flat translation buffer
/// into per-frame `TranslatedFrame` structs for caller convenience.
pub fn six_frame_translate(nt_seq: &[u8]) -> [TranslatedFrame; 6] {
    six_frame_translate_with_table(nt_seq, &crate::util::STANDARD_GENETIC_CODE)
}

/// blast-rs: Six-frame translation adapter over `BLAST_GetAllTranslations`;
/// not a direct NCBI C port.
fn six_frame_translate_with_table(
    nt_seq: &[u8],
    genetic_code: &'static [u8; 64],
) -> [TranslatedFrame; 6] {
    let nt_ncbi4na = encode_ncbi4na_sequence(nt_seq);
    let (translation_buffer, frame_offsets) =
        crate::util::blast_get_all_translations(&nt_ncbi4na, nt_ncbi4na.len(), genetic_code);
    let mut frames: [TranslatedFrame; 6] = std::array::from_fn(|ctx| {
        let frame = crate::util::blast_context_to_frame(ctx as u32);
        let begin = (frame_offsets[ctx] + 1) as usize;
        let end = frame_offsets[ctx + 1] as usize;
        let ascii: Vec<u8> = if begin < end {
            ncbistdaa_to_aminoacid_sequence(&translation_buffer[begin..end])
        } else {
            Vec::new()
        };
        TranslatedFrame {
            frame,
            protein: ascii,
            nt_len: nt_seq.len(),
        }
    });
    // C ordering: frames 1, 2, 3, -1, -2, -3 — already produced in that order
    // by `blast_context_to_frame`. Keep the array as-is.
    let _ = &mut frames;
    frames
}

/// A translated reading frame.
pub struct TranslatedFrame {
    pub frame: i32,
    pub protein: Vec<u8>,
    pub nt_len: usize,
}

// ── Additional types ────────────────────────────────────────────────────────

/// Parsed BLAST database header/defline.
#[derive(Debug, Clone)]
pub struct BlastDefLine {
    pub title: String,
    pub accession: String,
    pub taxid: u32,
}

/// Apply SEG masking on NCBIstdaa-encoded sequence in place.
/// blast-rs: Public SEG masking wrapper using default options; not a direct
/// NCBI C port.
pub fn apply_seg_ncbistdaa(seq: &mut [u8]) {
    apply_seg_ncbistdaa_with_options(seq, 12, 2.2, 2.5)
}

/// blast-rs: Public SEG masking wrapper with explicit options; not a direct
/// NCBI C port.
pub fn apply_seg_ncbistdaa_with_options(seq: &mut [u8], window: usize, locut: f64, hicut: f64) {
    if is_single_residue_low_complexity(seq) {
        for aa in seq {
            if *aa != NCBISTDAA_X {
                *aa = NCBISTDAA_X;
            }
        }
        return;
    }
    let mask = crate::filter::seg_filter_ncbistdaa(seq, window, locut, hicut);
    for r in &mask.regions {
        let start = r.start.max(0) as usize;
        let end = (r.end as usize).min(seq.len());
        for aa in &mut seq[start..end] {
            if *aa != NCBISTDAA_X {
                *aa = NCBISTDAA_X;
            }
        }
    }
}

/// blast-rs: SEG edge-case helper for homopolymer protein queries; not a direct
/// NCBI C port.
fn is_single_residue_low_complexity(seq: &[u8]) -> bool {
    if seq.len() < 12 {
        return false;
    }
    let mut residue = None;
    let mut true_count = 0usize;
    for &aa in seq {
        if !crate::encoding::is_ncbistdaa_composition_residue(aa) {
            continue;
        }
        true_count += 1;
        match residue {
            Some(prev) if prev != aa => return false,
            None => residue = Some(aa),
            _ => {}
        }
    }
    true_count >= 12
}

/// Compute a boolean mask indicating which positions are lowercase.
/// blast-rs: Public lowercase-mask helper; not a direct NCBI C port.
pub fn lowercase_mask(seq: &[u8]) -> Vec<bool> {
    seq.iter().map(|b| b.is_ascii_lowercase()).collect()
}

/// Apply repeat masking based on n-mer frequency.
/// blast-rs: Public repeat-mask mutator; not a direct NCBI C port.
pub fn apply_repeat_mask(seq: &mut [u8]) {
    let mask = repeat_mask(seq, 11, 2.0);
    for (i, masked) in mask.iter().enumerate() {
        if *masked {
            seq[i] = b'N';
        }
    }
}

/// Compute repeat mask: positions where n-mer frequency exceeds threshold.
/// blast-rs: Public repeat-mask helper; not a direct NCBI C port.
pub fn repeat_mask(seq: &[u8], nmer_size: usize, threshold: f64) -> Vec<bool> {
    let mut mask = vec![false; seq.len()];
    if seq.len() < nmer_size {
        return mask;
    }

    // Count n-mer occurrences
    let mut counts = std::collections::HashMap::new();
    for i in 0..=seq.len() - nmer_size {
        *counts.entry(&seq[i..i + nmer_size]).or_insert(0u32) += 1;
    }

    // Mark positions where n-mer count exceeds threshold * average
    let avg = seq.len() as f64 / 4f64.powi(nmer_size as i32).max(1.0);
    let cutoff = (threshold * avg).max(2.0) as u32;
    for i in 0..=seq.len() - nmer_size {
        if counts[&seq[i..i + nmer_size]] >= cutoff {
            for j in i..i + nmer_size {
                mask[j] = true;
            }
        }
    }
    mask
}

// ── BlastnSearch builder (moved from blastn.rs) ─────────────────────────────

/// Which query strand(s) to search.
#[derive(Clone, Copy, PartialEq)]
pub enum Strand {
    Both,
    Plus,
    Minus,
}

/// Builder for configuring and running a blastn search.
///
/// # Example
///
/// ```no_run
/// use blast_rs::BlastnSearch;
///
/// let results = BlastnSearch::new()
///     .query(b"ACGTACGTACGTACGTACGTACGTACGT")
///     .subject(b"NNNNNACGTACGTACGTACGTACGTACGTNNNNN")
///     .run();
///
/// for hit in &results {
///     println!("score={} evalue={:.2e}", hit.score, hit.evalue);
/// }
/// ```
pub struct BlastnSearch {
    pub word_size: usize,
    pub reward: i32,
    pub penalty: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub evalue: f64,
    pub dust: bool,
    pub strand: Strand,
    pub xdrop_gap_final: f64,
    query_raw: Vec<u8>,
    subject_raw: Vec<u8>,
}

impl Default for BlastnSearch {
    /// blast-rs: Native BlastnSearch default implementation; not a direct NCBI C port.
    fn default() -> Self {
        Self::new()
    }
}

impl BlastnSearch {
    /// blast-rs: Public BlastnSearch builder constructor; not a direct NCBI C port.
    pub fn new() -> Self {
        BlastnSearch {
            word_size: crate::stat::BLAST_WORDSIZE_NUCL as usize,
            reward: crate::stat::BLAST_REWARD,
            penalty: crate::stat::BLAST_PENALTY,
            gap_open: crate::stat::BLAST_GAP_OPEN_NUCL,
            gap_extend: crate::stat::BLAST_GAP_EXTN_NUCL,
            evalue: crate::stat::BLAST_EXPECT_VALUE,
            dust: true,
            strand: Strand::Both,
            xdrop_gap_final: crate::stat::BLAST_GAP_X_DROPOFF_FINAL_NUCL as f64,
            query_raw: Vec::new(),
            subject_raw: Vec::new(),
        }
    }

    /// blast-rs: BlastnSearch builder setter; not a direct NCBI C port.
    pub fn query(mut self, seq: &[u8]) -> Self {
        self.query_raw = seq.to_vec();
        self
    }
    /// blast-rs: BlastnSearch builder setter; not a direct NCBI C port.
    pub fn subject(mut self, seq: &[u8]) -> Self {
        self.subject_raw = seq.to_vec();
        self
    }
    /// blast-rs: BlastnSearch builder setter; not a direct NCBI C port.
    pub fn word_size(mut self, ws: usize) -> Self {
        self.word_size = ws;
        self
    }
    /// blast-rs: BlastnSearch builder setter; not a direct NCBI C port.
    pub fn reward(mut self, r: i32) -> Self {
        self.reward = r;
        self
    }
    /// blast-rs: BlastnSearch builder setter; not a direct NCBI C port.
    pub fn penalty(mut self, p: i32) -> Self {
        self.penalty = p;
        self
    }
    /// blast-rs: BlastnSearch builder setter; not a direct NCBI C port.
    pub fn gap_open(mut self, g: i32) -> Self {
        self.gap_open = g;
        self
    }
    /// blast-rs: BlastnSearch builder setter; not a direct NCBI C port.
    pub fn gap_extend(mut self, g: i32) -> Self {
        self.gap_extend = g;
        self
    }
    /// blast-rs: BlastnSearch builder setter; not a direct NCBI C port.
    pub fn evalue(mut self, e: f64) -> Self {
        self.evalue = e;
        self
    }
    /// blast-rs: BlastnSearch builder setter; not a direct NCBI C port.
    pub fn dust(mut self, d: bool) -> Self {
        self.dust = d;
        self
    }
    /// blast-rs: BlastnSearch builder setter; not a direct NCBI C port.
    pub fn strand(mut self, s: Strand) -> Self {
        self.strand = s;
        self
    }

    /// blast-rs: Public BlastnSearch execution wrapper; not a direct NCBI C port.
    pub fn run(&self) -> Vec<SearchHsp> {
        if self.query_raw.is_empty() || self.subject_raw.is_empty() {
            return Vec::new();
        }

        let mut query_plus = encode_blastna_sequence(&self.query_raw);

        if self.dust {
            let mask = crate::filter::dust_filter(&query_plus, 20, 64, 1);
            mask.apply(&mut query_plus, 14);
        }

        let query_plus_nomask = encode_blastna_sequence(&self.query_raw);
        let query_minus = reverse_complement_blastna_sequence(&query_plus);
        let query_minus_nomask = reverse_complement_blastna_sequence(&query_plus_nomask);

        let qp = if self.strand != Strand::Minus {
            &query_plus[..]
        } else {
            &[]
        };
        let qm = if self.strand != Strand::Plus {
            &query_minus[..]
        } else {
            &[]
        };
        let qpn = if self.strand != Strand::Minus {
            &query_plus_nomask[..]
        } else {
            &[]
        };
        let qmn = if self.strand != Strand::Plus {
            &query_minus_nomask[..]
        } else {
            &[]
        };

        let subject = encode_blastna_sequence(&self.subject_raw);

        let (ungapped_kbp, gapped_kbp) = blastn_api_kbps(
            &query_plus_nomask,
            self.reward,
            self.penalty,
            self.gap_open,
            self.gap_extend,
        );

        let (mut alpha, mut beta) = (0.0, 0.0);
        let _ = blast_get_nucl_alpha_beta(
            self.reward,
            self.penalty,
            self.gap_open,
            self.gap_extend,
            ungapped_kbp.lambda,
            ungapped_kbp.h,
            true,
            &mut alpha,
            &mut beta,
        );
        let db_length = self.subject_raw.len() as i64;
        let mut len_adj = 0;
        let _ = blast_compute_length_adjustment(
            gapped_kbp.k,
            gapped_kbp.log_k,
            alpha / gapped_kbp.lambda,
            beta,
            self.query_raw.len() as i32,
            db_length,
            1,
            Some(&mut len_adj),
        );
        let eff_db = (db_length - len_adj as i64).max(1);
        let search_space = eff_db as f64 * (self.query_raw.len() as i32 - len_adj).max(1) as f64;
        let x_dropoff = raw_gapped_xdrop_bits_f64(self.xdrop_gap_final, &gapped_kbp);

        blastn_gapped_search_nomask(
            qp,
            qm,
            qpn,
            qmn,
            &subject,
            self.word_size,
            self.reward,
            self.penalty,
            self.gap_open,
            self.gap_extend,
            x_dropoff,
            &gapped_kbp,
            search_space,
            self.evalue,
        )
    }
}

// ── Internal helpers ────────────────────────────────────────────────────────

/// blast-rs: Private API protein-query encoding/masking helper; not a direct
/// NCBI C port.
fn encode_protein_query(
    sequence: &[u8],
    filter_low_complexity: bool,
    seg_window: usize,
    seg_locut: f64,
    seg_hicut: f64,
) -> Vec<u8> {
    let mut query_aa = encode_ncbistdaa_sequence(sequence);
    if filter_low_complexity {
        apply_seg_ncbistdaa_with_options(&mut query_aa, seg_window, seg_locut, seg_hicut);
    }
    query_aa
}

#[cfg(test)]
mod low_complexity_tests {
    use super::*;

    #[test]
    fn seg_masks_short_protein_homopolymer_completely() {
        let mut seq = encode_ncbistdaa_sequence(b"AAAAAAAAAAAAAAAAAAAA");
        apply_seg_ncbistdaa(&mut seq);
        assert_eq!(
            crate::composition::blast_read_aa_composition(&seq, AA_SIZE).1,
            0
        );
    }

    #[test]
    fn blastp_filtered_homopolymer_query_returns_no_hits() {
        let tmp = tempfile::tempdir().unwrap();
        let base = tmp.path().join("db");
        let mut builder = BlastDbBuilder::new(DbType::Protein, "db");
        builder.add(SequenceEntry {
            title: "s1".to_string(),
            accession: "s1".to_string(),
            sequence: b"AAAAAAAAAAAAAAAAAAAA".to_vec(),
            taxid: None,
        });
        builder.write(&base).unwrap();
        let db = BlastDb::open(&base).unwrap();
        let params = SearchParams::blastp();
        assert!(blastp(&db, b"AAAAAAAAAAAAAAAAAAAA", &params).is_empty());
    }

    #[test]
    fn blastdb_builder_writes_v4_length_metadata() {
        let tmp = tempfile::tempdir().unwrap();

        let nt_base = tmp.path().join("ntdb");
        let mut nt_builder = BlastDbBuilder::new(DbType::Nucleotide, "ntdb");
        nt_builder.add(SequenceEntry {
            title: "nt1".to_string(),
            accession: "nt1".to_string(),
            sequence: b"ACGTAC".to_vec(),
            taxid: None,
        });
        nt_builder.add(SequenceEntry {
            title: "nt2".to_string(),
            accession: "nt2".to_string(),
            sequence: b"ACGTACGTACGT".to_vec(),
            taxid: None,
        });
        nt_builder.write(&nt_base).unwrap();
        let nt_db = BlastDb::open(&nt_base).unwrap();
        assert_eq!(nt_db.total_length, 18);
        assert_eq!(nt_db.max_seq_len, 12);
        assert!(is_blastdb_date(&nt_db.date));

        let prot_base = tmp.path().join("protdb");
        let mut prot_builder = BlastDbBuilder::new(DbType::Protein, "protdb");
        prot_builder.add(SequenceEntry {
            title: "p1".to_string(),
            accession: "p1".to_string(),
            sequence: b"MTEYK".to_vec(),
            taxid: None,
        });
        prot_builder.add(SequenceEntry {
            title: "p2".to_string(),
            accession: "p2".to_string(),
            sequence: b"MTEYKLVVVG".to_vec(),
            taxid: None,
        });
        prot_builder.write(&prot_base).unwrap();
        let prot_db = BlastDb::open(&prot_base).unwrap();
        assert_eq!(prot_db.total_length, 15);
        assert_eq!(prot_db.max_seq_len, 10);
        assert!(is_blastdb_date(&prot_db.date));
    }

    #[test]
    fn blastdb_builder_writes_ncbi_defline_set_with_ord_id() {
        let tmp = tempfile::tempdir().unwrap();
        let base = tmp.path().join("protdb");
        let mut builder = BlastDbBuilder::new(DbType::Protein, "protdb");
        builder.add(SequenceEntry {
            title: "first protein".to_string(),
            accession: "p1".to_string(),
            sequence: b"MTEYK".to_vec(),
            taxid: None,
        });
        builder.write(&base).unwrap();

        let raw_header = std::fs::read(db_component_path(&base, "phr")).unwrap();
        assert!(raw_header.starts_with(&[0x30, 0x80, 0x30, 0x80]));
        assert!(raw_header
            .windows(b"BL_ORD_ID".len())
            .any(|w| w == b"BL_ORD_ID"));

        let db = BlastDb::open(&base).unwrap();
        assert_eq!(db.get_accession(0).as_deref(), Some("p1"));
        assert_eq!(db.get_defline(0).as_deref(), Some("p1 first protein"));
    }

    #[test]
    fn blastdb_builder_appends_component_extensions() {
        let tmp = tempfile::tempdir().unwrap();

        let nt_base = tmp.path().join("ntdb.00");
        let mut nt_builder = BlastDbBuilder::new(DbType::Nucleotide, "ntdb");
        nt_builder.add(SequenceEntry {
            title: "nt1".to_string(),
            accession: "nt1".to_string(),
            sequence: b"ACGT".to_vec(),
            taxid: None,
        });
        nt_builder.write(&nt_base).unwrap();
        assert!(db_component_path(&nt_base, "nin").exists());
        assert!(db_component_path(&nt_base, "nsq").exists());
        assert!(db_component_path(&nt_base, "nhr").exists());
        assert!(!nt_base.with_extension("nin").exists());

        let prot_base = tmp.path().join("protdb.00");
        let mut prot_builder = BlastDbBuilder::new(DbType::Protein, "protdb");
        prot_builder.add(SequenceEntry {
            title: "p1".to_string(),
            accession: "p1".to_string(),
            sequence: b"MTEYK".to_vec(),
            taxid: None,
        });
        prot_builder.write(&prot_base).unwrap();
        assert!(db_component_path(&prot_base, "pin").exists());
        assert!(db_component_path(&prot_base, "psq").exists());
        assert!(db_component_path(&prot_base, "phr").exists());
        assert!(!prot_base.with_extension("pin").exists());
    }

    #[test]
    fn blastdb_builder_writes_nucleotide_ambiguity_data() {
        let tmp = tempfile::tempdir().unwrap();
        let base = tmp.path().join("ntdb");
        let mut builder = BlastDbBuilder::new(DbType::Nucleotide, "ntdb");
        builder.add(SequenceEntry {
            title: "ambiguous sequence".to_string(),
            accession: "amb1".to_string(),
            sequence: b"ACGTRRNNY".to_vec(),
            taxid: None,
        });
        builder.write(&base).unwrap();

        let db = BlastDb::open(&base).unwrap();
        let ambiguity = db.get_ambiguity_data(0).expect("ambiguity data");
        let decoded = crate::search::decode_packed_ncbi2na_with_ambiguity(
            db.get_sequence(0),
            db.get_seq_len(0) as usize,
            ambiguity,
        );

        assert_eq!(decoded, vec![0, 1, 2, 3, 4, 4, 14, 14, 5]);
    }

    fn is_blastdb_date(date: &str) -> bool {
        date.len() == 22
            && matches!(
                &date[..3],
                "Jan"
                    | "Feb"
                    | "Mar"
                    | "Apr"
                    | "May"
                    | "Jun"
                    | "Jul"
                    | "Aug"
                    | "Sep"
                    | "Oct"
                    | "Nov"
                    | "Dec"
            )
            && date.as_bytes()[3] == b' '
            && &date[6..8] == ", "
            && &date[12..14] == "  "
            && date.as_bytes()[16] == b':'
            && date.as_bytes()[19] == b' '
            && matches!(&date[20..], "AM" | "PM")
            && date.bytes().enumerate().all(|(i, b)| {
                matches!(i, 0..=2 | 3 | 6 | 7 | 12 | 13 | 16 | 19 | 20 | 21) || b.is_ascii_digit()
            })
    }
}

// ── Aliases for API compatibility ───────────────────────────────────────────

/// Get the scoring matrix for a given MatrixType.
/// blast-rs: Public matrix lookup wrapper; not a direct NCBI C port.
pub fn get_matrix(mt: MatrixType) -> &'static [[i32; AA_SIZE]; AA_SIZE] {
    match mt {
        MatrixType::Blosum45 => &crate::matrix::BLOSUM45,
        MatrixType::Blosum50 => &crate::matrix::BLOSUM50,
        MatrixType::Blosum62 => &crate::matrix::BLOSUM62,
        MatrixType::Blosum80 => &crate::matrix::BLOSUM80,
        MatrixType::Blosum90 => &crate::matrix::BLOSUM90,
        MatrixType::Pam30 => &crate::matrix::PAM30,
        MatrixType::Pam70 => &crate::matrix::PAM70,
        MatrixType::Pam250 => &crate::matrix::PAM250,
        MatrixType::Identity => &crate::matrix::IDENTITY,
    }
}

/// blast-rs: Matrix/gap lookup helper for public API parameters; not a direct
/// NCBI C port.
fn protein_gapped_params_for_matrix(
    matrix: MatrixType,
    gap_open: i32,
    gap_extend: i32,
) -> Option<GappedParams> {
    crate::stat::lookup_matrix_params(protein_matrix_name(matrix), gap_open, gap_extend)
}

/// blast-rs: Matrix/gap Karlin-block lookup helper; not a direct NCBI C port.
fn protein_kbp_for_matrix(matrix: MatrixType, gap_open: i32, gap_extend: i32) -> KarlinBlk {
    protein_gapped_params_for_matrix(matrix, gap_open, gap_extend)
        .map(|p| KarlinBlk {
            lambda: p.lambda,
            k: p.k,
            log_k: p.k.ln(),
            h: p.h,
            round_down: false,
        })
        .unwrap_or_else(|| {
            crate::stat::protein_ungapped_kbp_for_matrix(protein_matrix_name(matrix))
        })
}

/// blast-rs: Protein API effective search-space adapter; not a direct NCBI C port.
fn protein_api_search_space(
    params: &SearchParams,
    query_len: usize,
    stats_db_len: i64,
    num_subjects: i32,
    prot_kbp: &KarlinBlk,
) -> (i32, f64) {
    if params.effective_search_space > 0 {
        return (0, params.effective_search_space as f64);
    }

    // NCBI's `BLAST_CalcEffLengths` uses ungapped alpha/beta when
    // `scoring_options->gapped_calculation = FALSE` (mirrors `-ungapped`).
    // Pulls (alpha, beta) from row 0 of the matrix table (`BLAST_GetAlphaBeta`
    // ungapped branch, `blast_stat.c:3128`). For BLOSUM62 row 0 has
    // alpha=0.7916 and beta=-3.2.
    if params.ungapped {
        let ungapped_kbp =
            crate::stat::protein_ungapped_kbp_for_matrix(protein_matrix_name(params.matrix));
        let (alpha, beta) =
            crate::stat::lookup_matrix_ungapped_alpha_beta(protein_matrix_name(params.matrix))
                .unwrap_or((ungapped_kbp.lambda / ungapped_kbp.h.max(1e-9), 0.0));
        let mut len_adj = 0;
        let _ = blast_compute_length_adjustment(
            ungapped_kbp.k,
            ungapped_kbp.log_k,
            alpha / ungapped_kbp.lambda,
            beta,
            query_len as i32,
            stats_db_len,
            num_subjects,
            Some(&mut len_adj),
        );
        let eff_query = (query_len as i64 - len_adj as i64).max(1);
        let eff_db = (stats_db_len - num_subjects as i64 * len_adj as i64).max(1);
        return (len_adj, eff_query as f64 * eff_db as f64);
    }

    if let Some(gapped) =
        protein_gapped_params_for_matrix(params.matrix, params.gap_open, params.gap_extend)
    {
        let mut len_adj = 0;
        let _ = blast_compute_length_adjustment(
            prot_kbp.k,
            prot_kbp.log_k,
            gapped.alpha / prot_kbp.lambda,
            gapped.beta,
            query_len as i32,
            stats_db_len,
            num_subjects,
            Some(&mut len_adj),
        );
        let eff_query = (query_len as i64 - len_adj as i64).max(1);
        let eff_db = (stats_db_len - num_subjects as i64 * len_adj as i64).max(1);
        (len_adj, eff_query as f64 * eff_db as f64)
    } else {
        let len_adj = crate::stat::compute_length_adjustment(
            query_len as i32,
            stats_db_len,
            num_subjects,
            prot_kbp,
        );
        let search_space = crate::stat::compute_search_space(
            query_len as i64,
            stats_db_len,
            num_subjects,
            len_adj,
        );
        (len_adj, search_space)
    }
}

/// blast-rs: Matrix/gap Gumbel-block lookup helper; not a direct NCBI C port.
fn protein_gumbel_for_matrix(
    matrix: MatrixType,
    gap_open: i32,
    gap_extend: i32,
    db_length: i64,
) -> Option<GumbelBlk> {
    crate::stat::matrix_gumbel_blk(protein_matrix_name(matrix), gap_open, gap_extend, db_length)
}

/// blast-rs: Private enum to BLAST matrix-name adapter; not a direct NCBI C port.
fn protein_matrix_name(matrix: MatrixType) -> &'static str {
    match matrix {
        MatrixType::Blosum45 => "BLOSUM45",
        MatrixType::Blosum50 => "BLOSUM50",
        MatrixType::Blosum62 => "BLOSUM62",
        MatrixType::Blosum80 => "BLOSUM80",
        MatrixType::Blosum90 => "BLOSUM90",
        MatrixType::Pam30 => "PAM30",
        MatrixType::Pam70 => "PAM70",
        MatrixType::Pam250 => "PAM250",
        MatrixType::Identity => "IDENTITY",
    }
}

/// blast-rs: Public API default word-threshold selector; not a direct NCBI C port.
///
/// Mirrors NCBI's CLI word-threshold override in `CBlastAppArgs::SetWordSize`
/// (blast_args.cpp:288-299): when the QUERY is protein (`m_QueryIsProtein`,
/// i.e. blastp / tblastn) AND word_size > 4, NCBI switches to the compressed
/// AA lookup table and OVERRIDES the threshold purely by word size
/// (matrix-independent), as a flat assignment WITHOUT the translated `+1/+2`
/// adjustment:
///   ws == 5 → 19.3, ws == 6 → 21.0, ws >= 7 → 20.25.
/// For nucleotide-query translated programs (blastx / tblastx, where
/// `m_QueryIsProtein` is false) NCBI applies NO such override, so the historical
/// matrix-based value plus the translated adjustment is kept. Word sizes <= 4
/// (including the default ws3) keep the historical matrix-based logic on every
/// program — that path is byte-for-byte unchanged.
fn suggested_word_threshold(
    matrix: MatrixType,
    program: crate::program::ProgramType,
    word_size: usize,
) -> f64 {
    // NCBI CLI override: protein query + word_size > 4 → matrix-independent,
    // flat threshold (no translated adjustment applied afterwards).
    if crate::program::blast_query_is_protein(program) && word_size > 4 {
        return match word_size {
            5 => 19.3,
            6 => 21.0,
            _ => 20.25, // word_size >= 7
        };
    }
    let mut threshold = match matrix {
        MatrixType::Blosum45 => 14.0,
        MatrixType::Blosum62 => 11.0,
        MatrixType::Blosum80 => 12.0,
        MatrixType::Pam30 => 16.0,
        MatrixType::Pam70 => 14.0,
        MatrixType::Identity => 27.0,
        _ => 11.0,
    };
    if crate::program::blast_subject_is_translated(program) {
        threshold += 2.0;
    } else if matches!(program, crate::program::BLASTX | crate::program::TBLASTX) {
        threshold += 1.0;
    }
    threshold
}

/// Get a genetic code translation table by NCBI code number.
/// Returns a 64-byte table mapping codons to NCBIstdaa amino acid codes.
/// Supports the NCBI genetic codes accepted by BLAST 2.12
/// (1-6, 9-16, 21-31, 33). Unknown codes fall back to the standard code.
/// blast-rs: Public genetic-code lookup wrapper; not a direct NCBI C port.
pub fn get_codon_table(code: u8) -> &'static [u8; 64] {
    crate::util::lookup_genetic_code(code)
}

// ── Masking wrappers ────────────────────────────────────────────────────────

/// Apply DUST low-complexity masking in place (replaces masked positions with N).
/// blast-rs: Public DUST masking wrapper; not a direct NCBI C port.
pub fn apply_dust(seq: &mut [u8]) {
    let mask = crate::filter::dust_filter(seq, 20, 64, 1);
    for r in &mask.regions {
        let start = r.start.max(0) as usize;
        let end = (r.end as usize).saturating_add(1).min(seq.len());
        for i in start..end {
            seq[i] = b'N';
        }
    }
}

/// Apply SEG low-complexity masking in place (replaces masked positions with X).
/// blast-rs: Public SEG masking wrapper; not a direct NCBI C port.
pub fn apply_seg(seq: &mut [u8]) {
    let mask = crate::filter::seg_filter(seq, 12, 2.2);
    for r in &mask.regions {
        let start = r.start.max(0) as usize;
        let end = (r.end as usize).min(seq.len());
        for i in start..end {
            seq[i] = b'X';
        }
    }
}

/// Replace lowercase characters with N (nucleotide masking).
/// blast-rs: Public lowercase nucleotide masking helper; not a direct NCBI C port.
pub fn apply_lowercase_mask_nucleotide(seq: &mut [u8]) {
    for b in seq.iter_mut() {
        if b.is_ascii_lowercase() {
            *b = b'N';
        }
    }
}

/// Replace lowercase characters with X (protein masking).
/// blast-rs: Public lowercase protein masking helper; not a direct NCBI C port.
pub fn apply_lowercase_mask_protein(seq: &mut [u8]) {
    for b in seq.iter_mut() {
        if b.is_ascii_lowercase() {
            *b = b'X';
        }
    }
}

// ── Composition statistics ──────────────────────────────────────────────────

/// Compute amino acid composition (frequencies) from NCBIstdaa-encoded sequence.
/// blast-rs: Public composition convenience helper; not a direct NCBI C port.
pub fn composition_ncbistdaa(seq: &[u8]) -> [f64; 28] {
    let mut counts = [0u64; 28];
    for &b in seq {
        if (b as usize) < 28 {
            counts[b as usize] += 1;
        }
    }
    let total = seq.len().max(1) as f64;
    let mut freqs = [0.0; 28];
    for i in 0..28 {
        freqs[i] = counts[i] as f64 / total;
    }
    freqs
}

/// Adjust an E-value using per-sequence composition correction (mode 1).
///
/// q and r are amino acid frequency vectors (28 elements, NCBIstdaa indexed).
/// Returns the adjusted E-value, or the original if correction is inapplicable.
/// blast-rs: Public composition-adjusted e-value helper; not a direct NCBI C port.
#[allow(clippy::too_many_arguments)]
pub fn adjust_evalue(
    raw_evalue: f64,
    score: i32,
    q: &[f64; 28],
    r: &[f64; 28],
    matrix: &ScoringMatrix,
    lambda_standard: f64,
    k: f64,
    eff_query_len: usize,
    eff_db_len: u64,
) -> f64 {
    match find_adjusted_lambda(q, r, matrix, lambda_standard) {
        None => raw_evalue,
        Some(lambda_prime) => {
            (eff_query_len as f64)
                * (eff_db_len as f64)
                * k
                * (-(lambda_prime * score as f64)).exp()
        }
    }
}

/// Adjust E-value with mode selection:
///   0 = no adjustment (returns raw_evalue)
///   1 = unconditional composition-based lambda adjustment
///   2 = conditional: only apply if composition diverges significantly
///   3 = unconditional (always applied, even if expected_score >= 0)
/// blast-rs: Public composition-adjusted e-value helper with mode selector;
/// not a direct NCBI C port.
#[allow(clippy::too_many_arguments)]
pub fn adjust_evalue_with_mode(
    raw_evalue: f64,
    score: i32,
    q: &[f64; 28],
    r: &[f64; 28],
    matrix: &ScoringMatrix,
    lambda_standard: f64,
    k: f64,
    eff_query_len: usize,
    eff_db_len: u64,
    mode: u8,
) -> f64 {
    match mode {
        0 => raw_evalue,
        1 => adjust_evalue(
            raw_evalue,
            score,
            q,
            r,
            matrix,
            lambda_standard,
            k,
            eff_query_len,
            eff_db_len,
        ),
        2 => {
            let mu_bg = expected_score_with_bg(matrix);
            let mu_actual = compo_expected_score(q, r, matrix);
            let threshold = -0.2 * lambda_standard;
            if (mu_actual - mu_bg).abs() > threshold.abs() {
                adjust_evalue(
                    raw_evalue,
                    score,
                    q,
                    r,
                    matrix,
                    lambda_standard,
                    k,
                    eff_query_len,
                    eff_db_len,
                )
            } else {
                raw_evalue
            }
        }
        3 => {
            let eval_sum = |lam: f64| -> f64 {
                let mut sum = 0.0f64;
                for &i in crate::encoding::NCBISTDAA_STANDARD_RESIDUES.iter() {
                    for &j in crate::encoding::NCBISTDAA_STANDARD_RESIDUES.iter() {
                        let i = i as usize;
                        let j = j as usize;
                        let s = matrix.score(i as u8, j as u8) as f64;
                        sum += q[i] * r[j] * (lam * s).exp();
                    }
                }
                sum
            };
            let mut lo = 0.0f64;
            let mut hi = lambda_standard * 4.0;
            if eval_sum(hi) >= 1.0 || eval_sum(lo) < 1.0 {
                return raw_evalue;
            }
            for _ in 0..60 {
                let mid = (lo + hi) / 2.0;
                if eval_sum(mid) > 1.0 {
                    lo = mid;
                } else {
                    hi = mid;
                }
                if hi - lo < 1e-10 {
                    break;
                }
            }
            let lambda_prime = (lo + hi) / 2.0;
            (eff_query_len as f64)
                * (eff_db_len as f64)
                * k
                * (-(lambda_prime * score as f64)).exp()
        }
        _ => raw_evalue,
    }
}

/// Find adjusted lambda via bisection for composition-based statistics.
/// blast-rs: Lambda root-finding helper for API composition adjustment; not a
/// direct NCBI C port.
fn find_adjusted_lambda(
    q: &[f64; 28],
    r: &[f64; 28],
    matrix: &ScoringMatrix,
    lambda_standard: f64,
) -> Option<f64> {
    if compo_expected_score(q, r, matrix) >= 0.0 {
        return None;
    }
    let eval_sum = |lam: f64| -> f64 {
        let mut sum = 0.0f64;
        for &i in crate::encoding::NCBISTDAA_STANDARD_RESIDUES.iter() {
            for &j in crate::encoding::NCBISTDAA_STANDARD_RESIDUES.iter() {
                let i = i as usize;
                let j = j as usize;
                sum += q[i] * r[j] * (lam * matrix.score(i as u8, j as u8) as f64).exp();
            }
        }
        sum
    };
    if eval_sum(0.0) < 1.0 {
        return None;
    }
    let mut lo = 0.0f64;
    let mut hi = lambda_standard * 4.0;
    if eval_sum(hi) >= 1.0 {
        return None;
    }
    for _ in 0..60 {
        let mid = (lo + hi) / 2.0;
        if eval_sum(mid) > 1.0 {
            lo = mid;
        } else {
            hi = mid;
        }
        if hi - lo < 1e-10 {
            break;
        }
    }
    Some((lo + hi) / 2.0)
}

/// blast-rs: Expected-score helper for API composition adjustment; not a direct
/// NCBI C port.
fn compo_expected_score(q: &[f64; 28], r: &[f64; 28], matrix: &ScoringMatrix) -> f64 {
    let mut mu = 0.0f64;
    for &i in crate::encoding::NCBISTDAA_STANDARD_RESIDUES.iter() {
        for &j in crate::encoding::NCBISTDAA_STANDARD_RESIDUES.iter() {
            let i = i as usize;
            let j = j as usize;
            mu += q[i] * r[j] * matrix.score(i as u8, j as u8) as f64;
        }
    }
    mu
}

/// blast-rs: Background expected-score helper; not a direct NCBI C port.
fn expected_score_with_bg(matrix: &ScoringMatrix) -> f64 {
    let bg = &crate::matrix::AA_FREQUENCIES;
    // AA_FREQUENCIES is [f64; 20] in ACDEFGHIKLMNPQRSTVWY order, need to map to NCBIstdaa
    let mut freq28 = [0.0f64; 28];
    for (k, &idx) in crate::encoding::NCBISTDAA_STANDARD_RESIDUES
        .iter()
        .enumerate()
    {
        if k < bg.len() {
            freq28[idx as usize] = bg[k];
        }
    }
    compo_expected_score(&freq28, &freq28, matrix)
}

// ── Low-level PSSM functions ────────────────────────────────────────────────

/// Build a PSSM from search results (for PSI-BLAST iteration).
/// blast-rs: Public PSSM-construction helper; not a direct NCBI C port.
pub fn build_pssm(
    query: &[u8],
    results: &[SearchResult],
    inclusion_evalue: f64,
    matrix_type: MatrixType,
    _lambda: f64,
) -> crate::pssm::Pssm {
    let query_aa = encode_ncbistdaa_sequence(query);
    let matrix = get_matrix(matrix_type);
    let mut pssm = crate::pssm::Pssm::from_sequence(&query_aa, matrix);

    // Collect aligned subject sequences from results that pass the inclusion threshold.
    let aligned: Vec<Vec<u8>> = results
        .iter()
        .flat_map(|r| r.hsps.iter())
        .filter(|hsp| hsp.evalue <= inclusion_evalue)
        .filter_map(|hsp| project_subject_alignment_to_query(hsp, query_aa.len()))
        .collect();

    if !aligned.is_empty() {
        pssm.update_from_alignment_with_matrix(&aligned, protein_matrix_name(matrix_type));
        let params = SearchParams {
            matrix: matrix_type,
            ..SearchParams::blastp()
        };
        attach_psi_ancillary_gap_kbp(&mut pssm, &query_aa, &params);
    }
    pssm
}

/// blast-rs: Projects public HSP alignment into PSSM query coordinates; not a
/// direct NCBI C port.
fn project_subject_alignment_to_query(hsp: &Hsp, query_len: usize) -> Option<Vec<u8>> {
    if hsp.subject_aln.is_empty() {
        return None;
    }

    let query_start = hsp.query_start.min(query_len);
    if query_start >= query_len {
        return None;
    }

    let mut aligned = vec![NCBISTDAA_X; query_len];
    if !hsp.query_aln.is_empty() && hsp.query_aln.len() == hsp.subject_aln.len() {
        let mut query_pos = query_start;
        let mut placed = false;
        for (&q, &s) in hsp.query_aln.iter().zip(hsp.subject_aln.iter()) {
            if q == b'-' {
                continue;
            }
            if query_pos >= query_len {
                break;
            }
            aligned[query_pos] = if s == b'-' {
                NCBISTDAA_GAP
            } else {
                crate::encoding::aminoacid_to_ncbistdaa_base(s)
            };
            placed = true;
            query_pos += 1;
        }
        return placed.then_some(aligned);
    }

    let subject_aln = encode_ncbistdaa_sequence(&hsp.subject_aln);
    let copy_len = subject_aln.len().min(query_len - query_start);
    if copy_len == 0 {
        return None;
    }
    aligned[query_start..query_start + copy_len].copy_from_slice(&subject_aln[..copy_len]);
    Some(aligned)
}

/// Search a database using a PSSM instead of a substitution matrix.
/// blast-rs: Public PSSM search wrapper; not a direct NCBI C port.
pub fn search_with_pssm(
    db: &BlastDb,
    query: &[u8],
    pssm: &crate::pssm::Pssm,
    params: &SearchParams,
) -> Vec<SearchResult> {
    let query_aa = encode_ncbistdaa_sequence(query);
    let prot_kbp = protein_kbp_for_matrix(params.matrix, params.gap_open, params.gap_extend);

    let total_subj_len: usize = (0..db.num_oids)
        .map(|oid| db.get_seq_len(oid) as usize)
        .sum();
    let stats_db_len = statistical_db_length(params, total_subj_len as i64);
    let search_space = protein_api_search_space(
        params,
        pssm.length,
        stats_db_len,
        db.num_oids as i32,
        &prot_kbp,
    )
    .1;

    let subj_pairs = pssm_subject_pairs(db);
    let gumbel_blk = if params.effective_search_space > 0 {
        None
    } else {
        protein_gumbel_for_matrix(
            params.matrix,
            params.gap_open,
            params.gap_extend,
            stats_db_len,
        )
    };

    let hits = crate::pssm::psi_blast_iteration(
        pssm,
        &query_aa,
        &subj_pairs,
        params.evalue_threshold,
        search_space,
        prot_kbp.lambda,
        prot_kbp.k,
        gumbel_blk.as_ref(),
        params.gap_open,
        params.gap_extend,
        raw_gapped_xdrop_bits(params.x_drop_gapped, &prot_kbp),
        1,
        protein_matrix_name(params.matrix),
        params.comp_adjust > 0,
    );

    pssm_hits_to_search_results(db, &hits, &query_aa, pssm, &prot_kbp, params)
}

/// blast-rs: Converts BlastDb subjects into PSSM iteration pairs; not a direct
/// NCBI C port.
fn pssm_subject_pairs(db: &BlastDb) -> Vec<(String, Vec<u8>)> {
    (0..db.num_oids)
        .map(|oid| {
            let raw = db.get_sequence(oid);
            let aa: Vec<u8> = raw.iter().filter(|&&b| b != 0).copied().collect();
            (oid.to_string(), aa)
        })
        .collect()
}

/// Convert blastp `SearchResult` HSPs into `PsiBlastHit` records so the
/// PSSM update step can build alignments uniformly regardless of whether
/// iter 1 came from blastp (no PSSM yet) or `psi_blast_iteration`.
/// blast-rs: Adapter from public blastp results into PSI-BLAST hits; not a
/// direct NCBI C port.
fn blastp_results_to_psi_hits(
    results: &[SearchResult],
    subj_pairs: &[(String, Vec<u8>)],
) -> Vec<crate::pssm::PsiBlastHit> {
    let mut out = Vec::new();
    for result in results {
        let oid = result.subject_oid;
        let subj_id = oid.to_string();
        let subj_len = subj_pairs
            .iter()
            .find(|(id, _)| id == &subj_id)
            .map(|(_, seq)| seq.len())
            .unwrap_or(result.subject_len);
        for hsp in &result.hsps {
            out.push(crate::pssm::PsiBlastHit {
                subject_id: subj_id.clone(),
                score: hsp.score,
                bit_score: hsp.bit_score,
                evalue: hsp.evalue,
                query_start: hsp.query_start,
                query_end: hsp.query_end,
                subject_start: hsp.subject_start,
                subject_end: hsp.subject_end,
                align_len: hsp.alignment_length,
                subject_len: subj_len,
                query_aln: hsp.query_aln.clone(),
                subject_aln: hsp.subject_aln.clone(),
            });
        }
    }
    out
}

/// blast-rs: Adapter from PSI-BLAST hits into public SearchResult values; not a
/// direct NCBI C port.
fn pssm_hits_to_search_results(
    db: &BlastDb,
    hits: &[crate::pssm::PsiBlastHit],
    query_aa: &[u8],
    pssm: &crate::pssm::Pssm,
    prot_kbp: &KarlinBlk,
    params: &SearchParams,
) -> Vec<SearchResult> {
    let mut grouped: HashMap<u32, SearchResult> = HashMap::new();
    for h in hits {
        let Some(oid) = h.subject_id.parse::<u32>().ok() else {
            continue;
        };
        let raw = db.get_sequence(oid);
        let subject_aa: Vec<u8> = raw.iter().filter(|&&b| b != 0).copied().collect();
        let Some(hsp) = pssm_hit_to_hsp(h, query_aa, &subject_aa, pssm, prot_kbp) else {
            continue;
        };
        grouped
            .entry(oid)
            .or_insert_with(|| SearchResult {
                subject_oid: oid,
                subject_title: String::from_utf8_lossy(db.get_header(oid)).to_string(),
                subject_accession: db
                    .get_accession(oid)
                    .unwrap_or_else(|| format!("oid_{}", oid)),
                subject_len: h.subject_len,
                hsps: Vec::new(),
                taxids: vec![],
            })
            .hsps
            .push(hsp);
    }
    let mut results: Vec<SearchResult> = grouped.into_values().collect();

    apply_api_min_score_filter(&mut results, params.min_score);
    if let Some(culling_limit) = params.culling_limit {
        apply_api_culling_limit(&mut results, culling_limit, crate::program::BLASTP);
    }
    apply_api_max_hsps_limit(&mut results, params.max_hsps);
    results.sort_by(compare_search_results);
    if results.len() > params.max_target_seqs {
        results.truncate(params.max_target_seqs);
    }
    results
}

/// blast-rs: Adapter from one PSI-BLAST hit into a public HSP; not a direct
/// NCBI C port.
fn pssm_hit_to_hsp(
    h: &crate::pssm::PsiBlastHit,
    query_aa: &[u8],
    subject_aa: &[u8],
    pssm: &crate::pssm::Pssm,
    _prot_kbp: &KarlinBlk,
) -> Option<Hsp> {
    let query_end = h.query_end.min(query_aa.len());
    let subject_end = h.subject_end.min(subject_aa.len());
    if query_end <= h.query_start || subject_end <= h.subject_start {
        return None;
    }
    let query_slice = &query_aa[h.query_start..query_end];
    let subject_slice = &subject_aa[h.subject_start..subject_end];
    let query_aln = if h.query_aln.is_empty() {
        ncbistdaa_to_aminoacid_sequence(query_slice)
    } else {
        h.query_aln.clone()
    };
    let subject_aln = if h.subject_aln.is_empty() {
        ncbistdaa_to_aminoacid_sequence(subject_slice)
    } else {
        h.subject_aln.clone()
    };
    let mut num_identities = 0usize;
    let mut num_gaps = 0usize;
    let mut query_pos = h.query_start;
    let midline: Vec<u8> = query_aln
        .iter()
        .zip(subject_aln.iter())
        .map(|(&q_ascii, &s_ascii)| {
            if q_ascii == b'-' || s_ascii == b'-' {
                num_gaps += 1;
                if q_ascii != b'-' {
                    query_pos += 1;
                }
                b' '
            } else {
                let q = crate::encoding::aminoacid_to_ncbistdaa_base(q_ascii);
                let s = crate::encoding::aminoacid_to_ncbistdaa_base(s_ascii);
                let c = if q == s {
                    num_identities += 1;
                    ncbistdaa_to_aminoacid_base(q)
                } else if query_pos < pssm.length && pssm.score_at(query_pos, s) > 0 {
                    b'+'
                } else {
                    b' '
                };
                query_pos += 1;
                c
            }
        })
        .collect();
    let alignment_length = query_aln.len().min(subject_aln.len());
    Some(Hsp {
        score: h.score,
        bit_score: h.bit_score,
        evalue: h.evalue,
        query_start: h.query_start,
        query_end,
        subject_start: h.subject_start,
        subject_end,
        query_gapped_start: h.query_start,
        subject_gapped_start: h.subject_start,
        query_link_start: h.query_start,
        query_link_end: query_end,
        query_link_gapped_start: h.query_start,
        subject_link_start: h.subject_start,
        subject_link_end: subject_end,
        subject_link_gapped_start: h.subject_start,
        link_score: None,
        link_lambda: None,
        num_identities,
        num_gaps,
        alignment_length,
        query_aln,
        midline,
        subject_aln,
        query_frame: 0,
        subject_frame: 0,
        num_links: 1,
        comp_adjust_method: 0,
    })
}

/// Alias for `blastn_search` matching the old API.
/// Note: at the crate root this is available as `blastn_search` since the `blastn` name
/// is used by the builder-pattern module. Use `blast_rs::api::blastn()` for the function form.
/// blast-rs: Public API compatibility alias; not a direct NCBI C port.
pub fn blastn(db: &BlastDb, query: &[u8], params: &SearchParams) -> Vec<SearchResult> {
    blastn_search(db, query, params)
}

/// Alias for `blastp` — low-level protein search backend.
/// blast-rs: Public API compatibility alias; not a direct NCBI C port.
pub fn blast_search(db: &BlastDb, query: &[u8], params: &SearchParams) -> Vec<SearchResult> {
    blastp(db, query, params)
}

/// Low-level blastx search (alias for `blastx`).
/// blast-rs: Public API compatibility alias; not a direct NCBI C port.
pub fn blastx_search(db: &BlastDb, query: &[u8], params: &SearchParams) -> Vec<SearchResult> {
    blastx(db, query, params)
}

/// Low-level tblastn search (alias for `tblastn`).
/// blast-rs: Public API compatibility alias; not a direct NCBI C port.
pub fn tblastn_search(db: &BlastDb, query: &[u8], params: &SearchParams) -> Vec<SearchResult> {
    tblastn(db, query, params)
}

/// Low-level tblastx search (alias for `tblastx`).
/// blast-rs: Public API compatibility alias; not a direct NCBI C port.
pub fn tblastx_search(db: &BlastDb, query: &[u8], params: &SearchParams) -> Vec<SearchResult> {
    tblastx(db, query, params)
}

/// Low-level PSI-BLAST iteration (alias for `psiblast`).
/// blast-rs: Public API compatibility alias; not a direct NCBI C port.
pub fn psiblast_search(
    db: &BlastDb,
    query: &[u8],
    params: &PsiblastParams,
) -> (Vec<SearchResult>, crate::pssm::Pssm) {
    psiblast(db, query, params)
}

/// Six-frame translation with a specific genetic code number
/// (1=standard, 2=mito, …). Public-API wrapper over the 1-1 port of NCBI
/// `BLAST_GetAllTranslations`.
/// blast-rs: Public genetic-code translation wrapper; not a direct NCBI C port.
pub fn six_frame_translate_with_code(nt_seq: &[u8], genetic_code: u8) -> [TranslatedFrame; 6] {
    let table = crate::util::lookup_genetic_code(genetic_code);
    six_frame_translate_with_table(nt_seq, table)
}

// ── PSI-BLAST ───────────────────────────────────────────────────────────────

/// PSI-BLAST specific parameters.
#[derive(Debug, Clone)]
pub struct PsiblastParams {
    pub search: SearchParams,
    /// Number of PSI-BLAST iterations (default: 3).
    pub num_iterations: u32,
    /// E-value threshold for including hits in PSSM construction (default: 0.002).
    pub inclusion_evalue: f64,
    /// Fixed pseudocount for PSSM construction. Zero uses column-specific pseudocounts.
    pub pseudocount: i32,
    /// Gap-trigger threshold in bits for starting gapped extension.
    pub gap_trigger: f64,
    /// Optional initial PSSM loaded from a checkpoint.
    pub initial_pssm: Option<crate::pssm::Pssm>,
    /// Optional restart alignment rows in NCBIstdaa encoding.
    pub restart_alignment: Vec<Vec<u8>>,
}

/// PSI-BLAST run output with optional per-round PSSM snapshots.
#[derive(Debug, Clone)]
pub struct PsiblastRun {
    pub results: Vec<SearchResult>,
    pub pssm: crate::pssm::Pssm,
    pub round_results: Vec<Vec<SearchResult>>,
    pub round_pssms: Vec<crate::pssm::Pssm>,
    pub converged: bool,
}

impl PsiblastParams {
    /// blast-rs: Public PSI-BLAST parameter constructor; not a direct NCBI C port.
    pub fn new(search: SearchParams) -> Self {
        PsiblastParams {
            search,
            num_iterations: 3,
            inclusion_evalue: crate::stat::PSI_INCLUSION_ETHRESH,
            pseudocount: crate::stat::PSI_PSEUDO_COUNT_CONST,
            gap_trigger: crate::stat::BLAST_GAP_TRIGGER_PROT,
            initial_pssm: None,
            restart_alignment: Vec::new(),
        }
    }
    /// blast-rs: PsiblastParams builder setter; not a direct NCBI C port.
    pub fn num_iterations(mut self, v: u32) -> Self {
        self.num_iterations = v;
        self
    }
    /// blast-rs: PsiblastParams builder setter; not a direct NCBI C port.
    pub fn inclusion_evalue(mut self, v: f64) -> Self {
        self.inclusion_evalue = v;
        self
    }
    /// blast-rs: PsiblastParams builder setter; not a direct NCBI C port.
    pub fn pseudocount(mut self, v: i32) -> Self {
        self.pseudocount = v;
        self
    }
    /// blast-rs: PsiblastParams builder setter; not a direct NCBI C port.
    pub fn gap_trigger(mut self, v: f64) -> Self {
        self.gap_trigger = v;
        self
    }
    /// blast-rs: PsiblastParams builder setter; not a direct NCBI C port.
    pub fn initial_pssm(mut self, pssm: crate::pssm::Pssm) -> Self {
        self.initial_pssm = Some(pssm);
        self
    }
    /// blast-rs: PsiblastParams builder setter; not a direct NCBI C port.
    pub fn restart_alignment(mut self, aligned_seqs: Vec<Vec<u8>>) -> Self {
        self.restart_alignment = aligned_seqs;
        self
    }
}

/// Run iterative PSI-BLAST search.
/// Returns final-round hits and the resulting PSSM.
/// blast-rs: Public PSI-BLAST convenience wrapper; not a direct NCBI C port.
pub fn psiblast(
    db: &BlastDb,
    query: &[u8],
    params: &PsiblastParams,
) -> (Vec<SearchResult>, crate::pssm::Pssm) {
    let run = psiblast_with_rounds(db, query, params);
    (run.results, run.pssm)
}

/// Run iterative PSI-BLAST search and keep PSSM snapshots after each update.
/// blast-rs: Public PSI-BLAST iterative API pipeline; not a direct NCBI C port.
pub fn psiblast_with_rounds(db: &BlastDb, query: &[u8], params: &PsiblastParams) -> PsiblastRun {
    let query_aa = encode_ncbistdaa_sequence(query);
    let matrix = *get_matrix(params.search.matrix);

    // Initial PSSM from query
    let mut pssm = params
        .initial_pssm
        .clone()
        .unwrap_or_else(|| crate::pssm::Pssm::from_sequence(&query_aa, &matrix));
    if !params.restart_alignment.is_empty() {
        pssm.update_from_alignment_with_matrix_and_pseudocount(
            &params.restart_alignment,
            protein_matrix_name(params.search.matrix),
            (params.pseudocount > 0).then_some(params.pseudocount as f64),
        );
    }
    attach_psi_ancillary_gap_kbp(&mut pssm, &query_aa, &params.search);

    let prot_kbp = protein_kbp_for_matrix(
        params.search.matrix,
        params.search.gap_open,
        params.search.gap_extend,
    );
    let matrix_name = protein_matrix_name(params.search.matrix);
    let ungapped_kbp = crate::stat::query_specific_protein_ungapped_kbp_for_matrix(
        &query_aa,
        matrix_name,
        &matrix,
    );
    let gap_trigger_raw = protein_gap_trigger_raw(params.gap_trigger, &ungapped_kbp).max(1);

    let total_subj_len: usize = (0..db.num_oids)
        .map(|oid| db.get_seq_len(oid) as usize)
        .sum();
    let stats_db_len = statistical_db_length(&params.search, total_subj_len as i64);
    let search_space = protein_api_search_space(
        &params.search,
        query_aa.len(),
        stats_db_len,
        db.num_oids as i32,
        &prot_kbp,
    )
    .1;

    let mut final_results = Vec::new();
    let mut round_results = Vec::new();
    let mut round_pssms = Vec::new();
    let mut converged = false;
    let gumbel_blk = if params.search.effective_search_space > 0 {
        None
    } else {
        protein_gumbel_for_matrix(
            params.search.matrix,
            params.search.gap_open,
            params.search.gap_extend,
            stats_db_len,
        )
    };

    // Build subject pairs for PSI-BLAST iteration
    let subj_pairs = pssm_subject_pairs(db);
    let until_convergence = params.num_iterations == 0;
    let iteration_limit = if until_convergence {
        100
    } else {
        params.num_iterations
    };
    // NCBI's PSI-BLAST iter 1 runs standard `blastp` when no user-supplied
    // PSSM/MSA is provided (`prelim_stage.cpp` falls into the blastp path
    // until the first PSSM is built). For iter 2+, the engine switches to
    // PSSM scoring. Detect "starting from query sequence" by whether the
    // caller supplied an initial PSSM or a restart alignment.
    let iter1_uses_blastp = params.initial_pssm.is_none() && params.restart_alignment.is_empty();

    for iter_idx in 0..iteration_limit {
        let hits_from_iter1_blastp = if iter_idx == 0 && iter1_uses_blastp {
            // Run standard blastp once. Convert to the PSI-BLAST hit shape so
            // downstream PSSM-update code sees uniform input.
            Some(blastp(db, query, &params.search))
        } else {
            None
        };

        let hits = if hits_from_iter1_blastp.is_some() {
            // Convert SearchResult HSPs into PsiBlastHit for PSSM building.
            blastp_results_to_psi_hits(hits_from_iter1_blastp.as_ref().unwrap(), &subj_pairs)
        } else {
            crate::pssm::psi_blast_iteration(
                &pssm,
                &query_aa,
                &subj_pairs,
                params.search.evalue_threshold,
                search_space,
                prot_kbp.lambda,
                prot_kbp.k,
                gumbel_blk.as_ref(),
                params.search.gap_open,
                params.search.gap_extend,
                raw_gapped_xdrop_bits(params.search.x_drop_gapped, &prot_kbp),
                gap_trigger_raw,
                matrix_name,
                params.search.comp_adjust > 0,
            )
        };

        if hits.is_empty() {
            break;
        }

        final_results = if let Some(blastp_results) = hits_from_iter1_blastp {
            blastp_results
        } else {
            pssm_hits_to_search_results(db, &hits, &query_aa, &pssm, &prot_kbp, &params.search)
        };
        round_results.push(final_results.clone());

        // Update PSSM from aligned sequences
        let aligned: Vec<Vec<u8>> = hits
            .iter()
            .filter(|h| h.evalue <= params.inclusion_evalue)
            .filter_map(|h| {
                if !h.query_aln.is_empty() && h.query_aln.len() == h.subject_aln.len() {
                    project_psiblast_hit_alignment_to_query(h, pssm.length)
                } else {
                    subj_pairs
                        .iter()
                        .find(|(id, _)| id == &h.subject_id)
                        .and_then(|(_, seq)| {
                            if h.subject_start >= seq.len() {
                                return None;
                            }
                            let end = h.subject_end.min(seq.len());
                            if h.query_start >= pssm.length {
                                return None;
                            }
                            let mut aligned = vec![NCBISTDAA_X; pssm.length];
                            let copy_len = (end - h.subject_start).min(pssm.length - h.query_start);
                            aligned[h.query_start..h.query_start + copy_len]
                                .copy_from_slice(&seq[h.subject_start..h.subject_start + copy_len]);
                            Some(aligned)
                        })
                }
            })
            .collect();

        if !aligned.is_empty() {
            let previous_scores = until_convergence.then(|| pssm.scores.clone());
            pssm.update_from_alignment_with_matrix_and_pseudocount(
                &aligned,
                protein_matrix_name(params.search.matrix),
                (params.pseudocount > 0).then_some(params.pseudocount as f64),
            );
            attach_psi_ancillary_gap_kbp(&mut pssm, &query_aa, &params.search);
            let converged_now = previous_scores
                .as_ref()
                .is_some_and(|scores| scores == &pssm.scores);
            round_pssms.push(pssm.clone());
            if converged_now {
                converged = true;
                break;
            }
        } else {
            break;
        }
    }

    PsiblastRun {
        results: final_results,
        pssm,
        round_results,
        round_pssms,
        converged,
    }
}

fn attach_psi_ancillary_gap_kbp(
    pssm: &mut crate::pssm::Pssm,
    query_aa: &[u8],
    params: &SearchParams,
) {
    let _ = crate::pssm::psi_blast_add_ancillary_pssm_data(
        Some(pssm),
        Some(query_aa),
        Some(protein_matrix_name(params.matrix)),
        params.gap_open,
        params.gap_extend,
    );
}

/// blast-rs: Projects PSI-BLAST hit alignment into query coordinates; not a
/// direct NCBI C port.
fn project_psiblast_hit_alignment_to_query(
    hit: &crate::pssm::PsiBlastHit,
    query_len: usize,
) -> Option<Vec<u8>> {
    let mut query_pos = hit.query_start;
    if query_pos >= query_len || hit.query_aln.len() != hit.subject_aln.len() {
        return None;
    }

    let mut aligned = vec![NCBISTDAA_X; query_len];
    let mut placed = false;
    for (&q, &s) in hit.query_aln.iter().zip(hit.subject_aln.iter()) {
        if q == b'-' {
            continue;
        }
        if query_pos >= query_len {
            break;
        }
        aligned[query_pos] = if s == b'-' {
            NCBISTDAA_GAP
        } else {
            crate::encoding::aminoacid_to_ncbistdaa_base(s)
        };
        placed = true;
        query_pos += 1;
    }
    placed.then_some(aligned)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_hsp(score: i32, evalue: f64, query_start: usize, query_end: usize) -> Hsp {
        Hsp {
            score,
            bit_score: score as f64,
            evalue,
            query_start,
            query_end,
            subject_start: query_start,
            subject_end: query_end,
            query_gapped_start: query_start,
            subject_gapped_start: query_start,
            query_link_start: query_start,
            query_link_end: query_end,
            query_link_gapped_start: query_start,
            subject_link_start: query_start,
            subject_link_end: query_end,
            subject_link_gapped_start: query_start,
            link_score: None,
            link_lambda: None,
            num_identities: query_end.saturating_sub(query_start),
            num_gaps: 0,
            alignment_length: query_end.saturating_sub(query_start),
            query_aln: Vec::new(),
            midline: Vec::new(),
            subject_aln: Vec::new(),
            query_frame: 1,
            subject_frame: 0,
            num_links: 1,
            comp_adjust_method: 0,
        }
    }

    fn protein_hit(
        score: i32,
        query_start: usize,
        query_end: usize,
        subject_start: usize,
        subject_end: usize,
    ) -> crate::protein_lookup::ProteinHit {
        crate::protein_lookup::ProteinHit {
            query_start,
            query_end,
            subject_start,
            subject_end,
            score,
            num_ident: 0,
            align_length: query_end.saturating_sub(query_start) as i32,
            mismatches: 0,
            gap_opens: 0,
            qseq: None,
            sseq: None,
            scaled_score: None,
            adjusted_evalue: None,
            gapped_start_q: query_start,
            gapped_start_s: subject_start,
        }
    }

    fn internal_blast_hsp(
        score: i32,
        query_offset: i32,
        query_end: i32,
        subject_offset: i32,
        subject_end: i32,
    ) -> crate::hspstream::BlastHSP {
        crate::hspstream::BlastHSP {
            score,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1.0,
            query: crate::hspstream::BlastSeg {
                frame: 0,
                offset: query_offset,
                end: query_end,
                gapped_start: query_offset,
            },
            subject: crate::hspstream::BlastSeg {
                frame: 0,
                offset: subject_offset,
                end: subject_end,
                gapped_start: subject_offset,
            },
            context: 0,
            gap_info: None,
            num: 0,
            comp_adjustment_method: 0,
            pat_info: None,
            num_positives: 0,
            map_info: None,
        }
    }

    #[test]
    fn c_leaf_adjusted_offsets_handles_translated_subject_frame() {
        let mut hsp = internal_blast_hsp(10, 2, 7, 1, 4);
        hsp.subject.frame = -2;
        hsp.gap_info = Some(crate::gapinfo::GapEditScript::from_ops(vec![(
            crate::gapinfo::GapAlignOpType::Sub,
            5,
        )]));

        assert_eq!(
            blast_hsp_get_adjusted_offsets(crate::program::TBLASTN, &hsp, 20, 30),
            (3, 7, 26, 18)
        );
    }

    #[test]
    fn c_leaf_compute_num_identities_walks_gap_script() {
        let mut hsp = internal_blast_hsp(10, 0, 4, 0, 3);
        hsp.gap_info = Some(crate::gapinfo::GapEditScript::from_ops(vec![
            (crate::gapinfo::GapAlignOpType::Sub, 2),
            (crate::gapinfo::GapAlignOpType::Ins, 1),
            (crate::gapinfo::GapAlignOpType::Sub, 1),
        ]));

        s_compute_num_identities_blast_hsp(&mut hsp, &[1, 2, 9, 3], &[1, 4, 3], None);

        assert_eq!(hsp.num_ident, 2);
        assert_eq!(hsp.num_positives, 2);
    }

    #[test]
    fn c_leaf_hsp_test_matches_c_thresholds() {
        let mut hsp = internal_blast_hsp(64, 0, 4, 0, 4);
        hsp.num_ident = 2;
        let options = crate::options::HitSavingOptions {
            percent_identity: 60.0,
            min_hit_length: 3,
            ..Default::default()
        };

        assert!(s_hsp_test_blast_hsp(&hsp, &options, 4));
    }

    #[test]
    fn c_leaf_update_reevaluated_hsp_ungapped_updates_or_rejects() {
        let mut hsp = protein_hit(10, 0, 5, 10, 15);

        assert!(!s_update_reevaluated_hsp_ungapped(
            &mut hsp, 5, 7, 1, 4, 11, 14
        ));
        assert_eq!((hsp.score, hsp.query_start, hsp.query_end), (7, 1, 4));
        assert_eq!((hsp.subject_start, hsp.subject_end), (11, 14));

        assert!(s_update_reevaluated_hsp_ungapped(
            &mut hsp, 10, 9, 2, 5, 12, 15
        ));
        assert_eq!(hsp.score, 9);
    }

    #[test]
    fn tblastn_subject_frame_maps_to_original_translation_context() {
        let mut hsp = protein_hit_to_blast_hsp(protein_hit(40, 0, 8, 0, 8), 0, 0, 1);
        hsp.num = 7;

        assert_eq!(
            translated_frame_to_context(hsp.subject.frame as i32),
            Some(0)
        );
        assert_eq!(hsp.num, 7);
    }

    #[test]
    fn tblastn_subject_frame_rejects_invalid_context() {
        let mut hsp = protein_hit_to_blast_hsp(protein_hit(40, 0, 8, 0, 8), 0, 0, 1);
        hsp.subject.frame = 0;

        assert_eq!(translated_frame_to_context(hsp.subject.frame as i32), None);
    }

    #[test]
    fn tblastn_hsp_list_preserves_query_index_and_hsp_fields() {
        let mut hsp_list = blast_hsp_list_new(3);
        blast_hsp_list_push(
            &mut hsp_list,
            protein_hit_to_blast_hsp(protein_hit(40, 0, 8, 0, 8), 3, 0, 1),
        );

        assert_eq!(hsp_list.query_index, 3);
        assert_eq!(hsp_list.hspcnt, 1);
        assert_eq!(hsp_list.hsp_array[0].as_ref().unwrap().num, 0);
        assert_eq!(hsp_list.hsp_array[0].as_ref().unwrap().context, 3);
        assert_eq!(hsp_list.hsp_array[0].as_ref().unwrap().subject.frame, 1);
    }

    #[test]
    fn blastx_common_endpoint_purge_keeps_distinct_subject_frames() {
        let mut hsp_list = blast_hsp_list_new(0);
        blast_hsp_list_push(
            &mut hsp_list,
            protein_hit_to_blast_hsp(protein_hit(50, 0, 8, 0, 8), 0, 1, 1),
        );
        blast_hsp_list_push(
            &mut hsp_list,
            protein_hit_to_blast_hsp(protein_hit(40, 0, 8, 0, 8), 0, 1, -1),
        );

        purge_blastx_hsp_list_with_common_endpoints(&mut hsp_list);

        assert_eq!(hsp_list.hspcnt, 2);
        let frames: Vec<i16> = hsp_list
            .hsp_array
            .iter()
            .flatten()
            .map(|hsp| hsp.subject.frame)
            .collect();
        assert!(frames.contains(&1));
        assert!(frames.contains(&-1));
    }

    #[test]
    fn tblastn_common_endpoint_purge_matches_c_comparator_frame_boundary() {
        let mut hsp_list = blast_hsp_list_new(0);
        blast_hsp_list_push(
            &mut hsp_list,
            protein_hit_to_blast_hsp(protein_hit(100, 0, 8, 0, 8), 0, 0, -1),
        );
        blast_hsp_list_push(
            &mut hsp_list,
            protein_hit_to_blast_hsp(protein_hit(90, 0, 8, 0, 8), 0, 0, 1),
        );
        blast_hsp_list_push(
            &mut hsp_list,
            protein_hit_to_blast_hsp(protein_hit(80, 0, 8, 0, 8), 0, 0, -1),
        );

        purge_tblastn_hsp_list_with_common_endpoints(&mut hsp_list);

        assert_eq!(hsp_list.hspcnt, 3);
    }

    #[test]
    fn tblastn_link_head_match_ignores_linked_endpoint_drift() {
        let original = protein_hit_to_blast_hsp(protein_hit(80, 112, 139, 44, 72), 0, 0, 1);
        let mut linked = original.clone();
        linked.query.end = 166;
        linked.subject.end = 101;
        linked.evalue = 1.0e-20;
        linked.num = 2;

        assert!(same_tblastn_link_head(&original, &linked));

        let mut kbp = KarlinBlk::default();
        kbp.lambda = 0.267;
        kbp.k = 0.041;
        let mut api_hsp = tblastn_blast_hsp_to_api_hsp(&original, &kbp);
        api_hsp.evalue = linked.evalue;
        api_hsp.num_links = linked.num.max(1);

        assert_eq!((api_hsp.query_start, api_hsp.query_end), (112, 139));
        assert_eq!((api_hsp.subject_start, api_hsp.subject_end), (132, 216));
        assert_eq!(api_hsp.evalue, linked.evalue);
        assert_eq!(api_hsp.num_links, 2);
    }

    #[test]
    fn api_link_handoff_preserves_original_bounds_when_link_bridge_changes_end() {
        let mut hsps = vec![Hsp {
            subject_frame: 1,
            subject_start: 132,
            subject_end: 216,
            subject_gapped_start: 132,
            subject_link_start: 44,
            subject_link_end: 72,
            subject_link_gapped_start: 44,
            query_frame: 0,
            query_start: 112,
            query_end: 139,
            query_gapped_start: 112,
            query_link_start: 112,
            query_link_end: 139,
            query_link_gapped_start: 112,
            score: 80,
            bit_score: 42.0,
            evalue: 1.0,
            link_score: None,
            link_lambda: None,
            num_identities: 0,
            num_gaps: 0,
            alignment_length: 27,
            query_aln: Vec::new(),
            midline: Vec::new(),
            subject_aln: Vec::new(),
            num_links: 1,
            comp_adjust_method: 0,
        }];
        let mut linked = api_hsp_to_blast_hsp(
            &hsps[0],
            0,
            80,
            crate::hspstream::BlastSeg {
                frame: 0,
                offset: 112,
                end: 166,
                gapped_start: 112,
            },
            crate::hspstream::BlastSeg {
                frame: 1,
                offset: 44,
                end: 101,
                gapped_start: 44,
            },
        );
        linked.evalue = 1.0e-20;
        linked.num = 2;
        let linked_list = crate::hspstream::BlastHSPList {
            oid: 0,
            query_index: 0,
            hsp_array: vec![Some(linked)],
            hspcnt: 1,
            allocated: 1,
            hsp_max: i32::MAX,
            do_not_reallocate: false,
            best_evalue: 1.0e-20,
        };

        apply_blast_hsp_list_to_api_hsps(&mut hsps, linked_list, tblastn_link_key);

        assert_eq!(hsps.len(), 1);
        assert_eq!((hsps[0].query_start, hsps[0].query_end), (112, 139));
        assert_eq!((hsps[0].subject_start, hsps[0].subject_end), (132, 216));
        assert_eq!(hsps[0].evalue, 1.0e-20);
        assert_eq!(hsps[0].num_links, 2);
    }

    #[test]
    fn translated_s_chaining_alignment_keeps_colinear_chain_members() {
        let hits = vec![
            protein_hit(20, 0, 10, 0, 10),
            protein_hit(20, 40, 50, 40, 50),
            protein_hit(5, 80, 90, 5, 15),
        ];

        let kept = s_chaining_alignment(&hits, 11, 1, 22, 45);

        assert!(kept[0]);
        assert!(!kept[1]);
        assert!(!kept[2]);
    }

    #[test]
    fn translated_s_chaining_alignment_preserves_same_start_identity() {
        let hits = vec![protein_hit(40, 0, 10, 0, 10), protein_hit(5, 0, 5, 0, 5)];

        let kept = s_chaining_alignment(&hits, 11, 1, 8, 35);

        assert!(kept[0]);
        assert!(!kept[1]);
    }

    #[test]
    fn test_api_culling_limit_removes_enveloped_hsps_across_subjects() {
        let mut results = vec![
            SearchResult {
                subject_oid: 0,
                subject_title: "s0".to_string(),
                subject_accession: "s0".to_string(),
                subject_len: 200,
                hsps: vec![test_hsp(80, 1.0e-20, 10, 100)],
                taxids: Vec::new(),
            },
            SearchResult {
                subject_oid: 1,
                subject_title: "s1".to_string(),
                subject_accession: "s1".to_string(),
                subject_len: 200,
                hsps: vec![test_hsp(60, 1.0e-10, 20, 90)],
                taxids: Vec::new(),
            },
            SearchResult {
                subject_oid: 2,
                subject_title: "s2".to_string(),
                subject_accession: "s2".to_string(),
                subject_len: 200,
                hsps: vec![test_hsp(60, 1.0e-10, 110, 150)],
                taxids: Vec::new(),
            },
        ];

        apply_api_culling_limit(&mut results, 1, crate::program::BLASTP);

        let oids: Vec<u32> = results.iter().map(|result| result.subject_oid).collect();
        assert_eq!(oids, vec![0, 2]);
    }

    #[test]
    fn test_api_culling_limit_accepted_hsp_removes_dominated_earlier_hsp() {
        let mut results = vec![
            SearchResult {
                subject_oid: 0,
                subject_title: "lower_evalue".to_string(),
                subject_accession: "lower_evalue".to_string(),
                subject_len: 200,
                hsps: vec![test_hsp(20, 1.0e-50, 10, 50)],
                taxids: Vec::new(),
            },
            SearchResult {
                subject_oid: 1,
                subject_title: "higher_score".to_string(),
                subject_accession: "higher_score".to_string(),
                subject_len: 200,
                hsps: vec![test_hsp(80, 1.0e-20, 10, 50)],
                taxids: Vec::new(),
            },
        ];

        apply_api_culling_limit(&mut results, 1, crate::program::BLASTP);

        let oids: Vec<u32> = results.iter().map(|result| result.subject_oid).collect();
        assert_eq!(oids, vec![1]);
    }

    #[test]
    fn test_api_culling_limit_two_dominators_exhaust_earlier_hsp_merit() {
        let mut results = vec![
            SearchResult {
                subject_oid: 0,
                subject_title: "lower_evalue".to_string(),
                subject_accession: "lower_evalue".to_string(),
                subject_len: 200,
                hsps: vec![test_hsp(20, 1.0e-50, 10, 50)],
                taxids: Vec::new(),
            },
            SearchResult {
                subject_oid: 1,
                subject_title: "higher_score_1".to_string(),
                subject_accession: "higher_score_1".to_string(),
                subject_len: 200,
                hsps: vec![test_hsp(80, 1.0e-20, 10, 50)],
                taxids: Vec::new(),
            },
            SearchResult {
                subject_oid: 2,
                subject_title: "higher_score_2".to_string(),
                subject_accession: "higher_score_2".to_string(),
                subject_len: 200,
                hsps: vec![test_hsp(70, 1.0e-18, 10, 50)],
                taxids: Vec::new(),
            },
        ];

        apply_api_culling_limit(&mut results, 2, crate::program::BLASTP);

        let oids: Vec<u32> = results.iter().map(|result| result.subject_oid).collect();
        assert_eq!(oids, vec![1, 2]);
    }

    #[test]
    fn test_api_min_score_filter_removes_low_scoring_hsps_and_empty_results() {
        let mut results = vec![
            SearchResult {
                subject_oid: 0,
                subject_title: "s0".to_string(),
                subject_accession: "s0".to_string(),
                subject_len: 200,
                hsps: vec![
                    test_hsp(80, 1.0e-20, 10, 100),
                    test_hsp(30, 1.0e-5, 110, 130),
                ],
                taxids: Vec::new(),
            },
            SearchResult {
                subject_oid: 1,
                subject_title: "s1".to_string(),
                subject_accession: "s1".to_string(),
                subject_len: 200,
                hsps: vec![test_hsp(20, 1.0e-3, 20, 40)],
                taxids: Vec::new(),
            },
        ];

        apply_api_min_score_filter(&mut results, 50);

        assert_eq!(results.len(), 1);
        assert_eq!(results[0].subject_oid, 0);
        assert_eq!(results[0].hsps.len(), 1);
        assert_eq!(results[0].hsps[0].score, 80);
    }

    #[test]
    fn test_late_protein_midline_preserves_translated_lowercase_match_text() {
        assert_eq!(
            protein_midline(b"i*k", b"IQK", ProteinMidlineStyle::AnyMismatchPlus),
            b"I+K"
        );
    }

    #[test]
    fn test_tblastx_search_params_use_ncbi_ungapped_program_defaults() {
        let params = SearchParams::tblastx();

        assert_eq!(params.comp_adjust, 0);
        assert_eq!(params.max_intron_length, 0);
        assert_eq!(
            params.x_drop_gapped,
            crate::stat::BLAST_GAP_X_DROPOFF_TBLASTX
        );
        assert_eq!(
            params.x_drop_final,
            crate::stat::BLAST_GAP_X_DROPOFF_FINAL_TBLASTX
        );
    }

    #[test]
    fn test_tblastx_ungapped_reevaluation_uses_masked_query_sequence() {
        let mut hit = protein_hit(15, 0, 3, 0, 3);
        let mut matrix = [[-1; AA_SIZE]; AA_SIZE];
        let a = crate::encoding::AMINOACID_TO_NCBISTDAA[b'A' as usize] as usize;
        let x = crate::encoding::NCBISTDAA_X as usize;
        matrix[a][a] = 5;
        matrix[x][a] = -5;

        // CTBlastxOptionsHandle inherits protein SEG defaults and does not set
        // mask_at_hash. In C, BLAST_MainSetUp therefore calls
        // BlastSetUp_MaskQuery: query_blk->sequence is masked and
        // sequence_nomask is saved separately. Blast_HSPListReevaluateUngapped
        // passes query_blk->sequence to Blast_HSPReevaluateWithAmbiguitiesUngapped.
        let query_masked = [crate::encoding::NCBISTDAA_X, a as u8, a as u8];
        let subject = [a as u8, a as u8, a as u8];

        assert!(reevaluate_ungapped_translated_hsp(
            &mut hit,
            &query_masked,
            &subject,
            &matrix,
            1,
        ));
        assert_eq!(hit.score, 10);
        assert_eq!((hit.query_start, hit.query_end), (1, 3));
        assert_eq!((hit.subject_start, hit.subject_end), (1, 3));
    }

    #[test]
    fn test_translated_search_params_use_ncbi_default_zero_max_intron() {
        assert_eq!(SearchParams::blastx().max_intron_length, 0);
        assert_eq!(SearchParams::tblastn().max_intron_length, 0);
        assert_eq!(SearchParams::tblastx().max_intron_length, 0);
    }

    #[test]
    fn translated_composition_window_bounds_use_bordered_query_window_for_blastx() {
        let hits = vec![
            protein_hit(40, 250, 270, 5, 25),
            protein_hit(35, 320, 340, 50, 70),
        ];

        assert_eq!(
            translated_composition_window_bounds(&hits, 1000, true, true),
            (50, 540)
        );
        assert_eq!(
            translated_composition_window_bounds(&hits, 500, false, false),
            (0, 500)
        );
    }

    #[test]
    fn translated_composition_window_bounds_use_bordered_subject_window_for_tblastn() {
        let hits = vec![
            protein_hit(40, 5, 25, 250, 270),
            protein_hit(35, 50, 70, 320, 340),
        ];

        assert_eq!(
            translated_composition_window_bounds(&hits, 500, true, false),
            (0, 500)
        );
        assert_eq!(
            translated_composition_window_bounds(&hits, 1000, false, true),
            (50, 540)
        );
    }

    #[test]
    fn test_blastn_defaults_use_soft_lookup_masking() {
        assert!(SearchParams::blastn().soft_masking);
    }

    #[test]
    fn test_blastn_lookup_masks_lowercase_without_mutating_nomask_copy() {
        let raw = b"ACGTaaaaACGT";
        let mut lookup = encode_blastna_sequence(raw);
        let nomask = lookup.clone();

        apply_blastn_query_lookup_masks(&mut lookup, raw, false, true);

        assert_eq!(nomask, encode_blastna_sequence(raw));
        assert_eq!(&lookup[..4], &nomask[..4]);
        assert!(lookup[4..8]
            .iter()
            .all(|&base| base == crate::filter::K_NUCL_MASK));
        assert_eq!(&lookup[8..], &nomask[8..]);
    }

    #[test]
    fn test_blastn_soft_masking_false_hard_masks_runtime_query() {
        let tmp = tempfile::tempdir().unwrap();
        let base = tmp.path().join("db");
        let query = b"ACGTACGTaaaaACGTACGT";
        let subject = b"ACGTACGTAAAACGTACGT";
        let mut builder = BlastDbBuilder::new(DbType::Nucleotide, "db");
        builder.add(SequenceEntry {
            title: "subject".to_string(),
            accession: "subject".to_string(),
            sequence: subject.to_vec(),
            taxid: None,
        });
        builder.write(&base).unwrap();
        let db = BlastDb::open(&base).unwrap();

        let base_params = SearchParams::blastn()
            .word_size(4)
            .evalue(1.0e20)
            .filter_low_complexity(false)
            .lcase_masking(true)
            .sum_stats(false);
        let soft = blastn_search(&db, query, &base_params.clone().soft_masking(true));
        let hard = blastn_search(&db, query, &base_params.soft_masking(false));

        let soft_best = soft[0].hsps.iter().max_by_key(|hsp| hsp.score).unwrap();
        let hard_best = hard[0].hsps.iter().max_by_key(|hsp| hsp.score).unwrap();
        assert!(
            soft_best.score > hard_best.score,
            "soft masking should score through lowercase bases; hard masking should not"
        );
        assert!(
            !soft_best.query_aln.contains(&b'N') && hard_best.query_aln.contains(&b'N'),
            "hard masking should expose masked lowercase bases as N in the runtime query"
        );
    }

    #[test]
    fn test_composition_prelim_evalue_uses_ncbi_cbs_stretch() {
        let base = SearchParams::blastp().evalue(2.0);
        let mode0 = base.clone().comp_adjust(0);
        let mode1 = base.clone().comp_adjust(1);
        let mode2 = base.clone().comp_adjust(2);
        let mode3 = base.comp_adjust(3);

        assert_eq!(composition_prelim_evalue(&mode0), 2.0);
        assert_eq!(composition_prelim_evalue(&mode1), 2.0);
        assert_eq!(composition_prelim_evalue(&mode2), 10.0);
        assert_eq!(composition_prelim_evalue(&mode3), 10.0);
    }

    #[test]
    fn compo_expected_score_ignores_ambiguous_ncbistdaa_slots() {
        let mut scores = [[0i32; AA_SIZE]; AA_SIZE];
        scores[1][1] = 7; // A/A, standard.
        scores[2][1] = 1000; // B/A, ambiguous and must not contribute.
        scores[crate::encoding::NCBISTDAA_X as usize][1] = 1000; // X/A, also excluded.
        let matrix = ScoringMatrix {
            scores,
            min_score: 0,
            name: "TEST".to_string(),
        };
        let mut q = [0.0f64; 28];
        let mut r = [0.0f64; 28];
        q[1] = 0.5;
        q[2] = 0.25;
        q[crate::encoding::NCBISTDAA_X as usize] = 0.25;
        r[1] = 1.0;

        assert_eq!(compo_expected_score(&q, &r, &matrix), 3.5);
    }

    #[test]
    fn test_blastn_search_applies_culling_limit() {
        let tmp = tempfile::tempdir().unwrap();
        let base = tmp.path().join("db");
        let query = b"ACGTCGATGCTAGCTAGGCTAACCGTATCGGATCCGTAAGCTTAGCTAGGATCCGATACGGTTAGCCTA";
        let mut builder = BlastDbBuilder::new(DbType::Nucleotide, "db");
        builder.add(SequenceEntry {
            title: "wide".to_string(),
            accession: "wide".to_string(),
            sequence: query.to_vec(),
            taxid: None,
        });
        builder.add(SequenceEntry {
            title: "contained".to_string(),
            accession: "contained".to_string(),
            sequence: query.to_vec(),
            taxid: None,
        });
        builder.write(&base).unwrap();
        let db = BlastDb::open(&base).unwrap();

        let base_params = SearchParams::blastn()
            .word_size(4)
            .evalue(1.0e20)
            .filter_low_complexity(false)
            .sum_stats(false);
        let uncull = blastn_search(&db, query, &base_params);
        let culled = blastn_search(&db, query, &base_params.clone().culling_limit(Some(1)));

        assert_eq!(uncull.len(), 2);
        assert_eq!(culled.len(), 1);
        assert!(!culled[0].hsps.is_empty());
    }

    #[test]
    fn test_builder_defaults() {
        let s = BlastnSearch::new();
        assert_eq!(s.word_size, 11);
        assert_eq!(s.reward, 1);
        assert_eq!(s.penalty, -3);
        assert_eq!(s.evalue, 10.0);
    }

    #[test]
    fn test_perfect_match() {
        let results = BlastnSearch::new()
            .word_size(7)
            .dust(false)
            .query(b"ACGTACGTACGTACGTACGTACGT")
            .subject(b"NNNNNACGTACGTACGTACGTACGTNNNNN")
            .run();
        assert!(!results.is_empty(), "Should find perfect match");
        assert!(results[0].score > 0);
    }

    #[test]
    fn test_no_match() {
        let results = BlastnSearch::new()
            .word_size(7)
            .dust(false)
            .strand(Strand::Plus)
            .query(b"AAAAAAAAAAAAAAAA")
            .subject(b"CCCCCCCCCCCCCCCC")
            .run();
        assert!(results.is_empty(), "Should find no match");
    }

    #[test]
    fn test_custom_scoring() {
        let results = BlastnSearch::new()
            .word_size(7)
            .reward(2)
            .penalty(-3)
            .gap_open(5)
            .gap_extend(2)
            .dust(false)
            .query(b"ACGTACGTACGTACGT")
            .subject(b"ACGTACGTACGTACGT")
            .run();
        assert!(!results.is_empty());
    }

    #[test]
    fn test_blastn_sum_stats_helper_updates_linked_evalues() {
        let kbp = KarlinBlk {
            lambda: 1.3,
            k: 0.71,
            log_k: 0.71_f64.ln(),
            h: 1.0,
            round_down: false,
        };
        let searchsp = 1000.0;
        let initial = kbp.raw_to_evalue(50, searchsp);
        let mk_hsp = |query_start, query_end, subject_start, subject_end| SearchHsp {
            query_start,
            query_end,
            subject_start,
            subject_end,
            query_gapped_start: query_start,
            subject_gapped_start: subject_start,
            score: 50,
            bit_score: kbp.raw_to_bit(50),
            evalue: initial,
            num_ident: 25,
            align_length: 25,
            mismatches: 0,
            gap_opens: 0,
            context: 0,
            qseq: None,
            sseq: None,
        };
        let mut hsps = vec![mk_hsp(10, 60, 100, 150), mk_hsp(70, 120, 170, 220)];

        apply_blastn_linked_sum_stats_to_search_hsps(
            &mut hsps, 500, 5000, &kbp, &kbp, searchsp, searchsp, 0, 0,
        );

        assert!(hsps.iter().any(|hsp| hsp.evalue != initial));
        assert!(hsps.iter().all(|hsp| hsp.evalue <= initial));
    }

    #[test]
    fn test_blastn_api_sum_stats_default_updates_multi_hsp_evalues() {
        let tmp = tempfile::tempdir().unwrap();
        let base = tmp.path().join("db");
        let query = b"ACGTCGATGCTAGCTAGGCTAACCGTATCGGATCCGTAAGCTTAGCTA";
        let mut subject = Vec::new();
        subject.extend_from_slice(query);
        subject.extend_from_slice(b"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT");
        subject.extend_from_slice(query);

        let mut builder = BlastDbBuilder::new(DbType::Nucleotide, "db");
        builder.add(SequenceEntry {
            title: "sum_stats_subject".to_string(),
            accession: "sum_stats_subject".to_string(),
            sequence: subject,
            taxid: None,
        });
        builder.write(&base).unwrap();
        let db = BlastDb::open(&base).unwrap();

        let raw_params = SearchParams::blastn()
            .word_size(7)
            .evalue(1.0e20)
            .filter_low_complexity(false)
            .sum_stats(false);
        let linked_params = raw_params.clone().sum_stats(true);

        let raw = blastn_search(&db, query, &raw_params);
        let linked = blastn_search(&db, query, &linked_params);
        assert_eq!(raw.len(), 1);
        assert_eq!(linked.len(), 1);
        assert!(raw[0].hsps.len() > 1);
        assert_eq!(raw[0].hsps.len(), linked[0].hsps.len());
        assert!(
            raw[0]
                .hsps
                .iter()
                .zip(&linked[0].hsps)
                .any(|(raw_hsp, linked_hsp)| linked_hsp.evalue < raw_hsp.evalue),
            "default linked sum-stat path should improve at least one HSP e-value"
        );
    }

    #[test]
    fn tblastn_singleton_sum_stats_repro_expected_ncbi_evalue() {
        let tmp = tempfile::tempdir().unwrap();
        let base = tmp.path().join("tblastn_db");
        let mut builder = BlastDbBuilder::new(DbType::Nucleotide, "tblastn_db");
        builder.add(SequenceEntry {
            title: "s1".to_string(),
            accession: "s1".to_string(),
            sequence: b"ATGAAATTTCTGATTCTGCTGTTTAAATTTATGAAATTTCTGATTCTGCTGTTT".to_vec(),
            taxid: None,
        });
        builder.write(&base).unwrap();
        let db = BlastDb::open(&base).unwrap();

        let params = SearchParams::tblastn()
            .filter_low_complexity(false)
            .comp_adjust(0)
            .sum_stats(true)
            .evalue(1000.0);
        let results = tblastn(&db, b"MKFLILLFMKFLILLF", &params);

        assert_eq!(results.len(), 1);
        assert_eq!(results[0].hsps.len(), 1);
        let hsp = &results[0].hsps[0];
        assert_eq!(hsp.score, 63);
        assert_eq!(hsp.alignment_length, 18);
        assert_eq!((hsp.query_start, hsp.query_end), (0, 16));
        assert_eq!((hsp.subject_start, hsp.subject_end), (0, 54));
        assert!(
            (hsp.evalue - 6.026535317927053e-09).abs() <= 1.0e-20,
            "public Blast_TracebackFromHSPList -> BLAST_LinkHsps handoff reports {}",
            hsp.evalue
        );
    }

    #[test]
    fn tblastn_batch_singleton_sum_stats_matches_single_query() {
        let tmp = tempfile::tempdir().unwrap();
        let base = tmp.path().join("tblastn_batch_db");
        let mut builder = BlastDbBuilder::new(DbType::Nucleotide, "tblastn_batch_db");
        builder.add(SequenceEntry {
            title: "s1".to_string(),
            accession: "s1".to_string(),
            sequence: b"ATGAAATTTCTGATTCTGCTGTTTAAATTTATGAAATTTCTGATTCTGCTGTTT".to_vec(),
            taxid: None,
        });
        builder.write(&base).unwrap();
        let db = BlastDb::open(&base).unwrap();

        let query = b"MKFLILLFMKFLILLF";
        let params = SearchParams::tblastn()
            .filter_low_complexity(false)
            .comp_adjust(0)
            .sum_stats(true)
            .evalue(1000.0);
        let single = tblastn(&db, query, &params);
        let batch = tblastn_batch(&db, &[query.as_slice()], &params);

        assert_eq!(batch.len(), 1);
        assert_eq!(single.len(), 1);
        assert_eq!(batch[0].len(), 1);
        let single_hsp = &single[0].hsps[0];
        let batch_hsp = &batch[0][0].hsps[0];
        assert_eq!(batch_hsp.score, single_hsp.score);
        assert_eq!(batch_hsp.alignment_length, single_hsp.alignment_length);
        assert_eq!(
            (batch_hsp.query_start, batch_hsp.query_end),
            (single_hsp.query_start, single_hsp.query_end)
        );
        assert_eq!(
            (batch_hsp.subject_start, batch_hsp.subject_end),
            (single_hsp.subject_start, single_hsp.subject_end)
        );
        assert!(
            (batch_hsp.evalue - single_hsp.evalue).abs()
                <= single_hsp.evalue.abs().max(1.0) * 1.0e-12,
            "batch evalue {} differed from single {}",
            batch_hsp.evalue,
            single_hsp.evalue
        );
    }

    #[test]
    fn test_blastp_batch_matches_single_query_statistics() {
        let tmp = tempfile::tempdir().unwrap();
        let base = tmp.path().join("db");
        let subject = b"MKKWLFGFLGVALIVVCSVFGYVSYQKHEGEVFKQNIEKKMPVDQINAHAKSY";
        let query_a = b"MKKWLFGFLGVALIVVCSVFGYVSYQKHEGEVFKQNIEKK";
        let query_b = b"GEVFKQNIEKKMPVDQINAHAKSY";

        let mut builder = BlastDbBuilder::new(DbType::Protein, "db");
        builder.add(SequenceEntry {
            title: "protein_subject".to_string(),
            accession: "protein_subject".to_string(),
            sequence: subject.to_vec(),
            taxid: None,
        });
        builder.write(&base).unwrap();
        let db = BlastDb::open(&base).unwrap();
        let params = SearchParams::blastp()
            .word_size(3)
            .word_threshold(9.0)
            .evalue(1.0e20)
            .filter_low_complexity(false)
            .comp_adjust(0)
            .sum_stats(false);

        let batch = blastp_batch(&db, &[query_a.as_slice(), query_b.as_slice()], &params);
        let singles = [blastp(&db, query_a, &params), blastp(&db, query_b, &params)];

        assert_eq!(batch.len(), singles.len());
        for (batch_results, single_results) in batch.iter().zip(singles.iter()) {
            assert_eq!(batch_results.len(), single_results.len());
            assert!(!batch_results.is_empty());
            let batch_hsp = &batch_results[0].hsps[0];
            let single_hsp = &single_results[0].hsps[0];
            assert_eq!(batch_hsp.score, single_hsp.score);
            assert_eq!(batch_hsp.alignment_length, single_hsp.alignment_length);
            assert!(
                (batch_hsp.evalue - single_hsp.evalue).abs()
                    <= single_hsp.evalue.abs().max(1.0) * 1.0e-12,
                "batch evalue {} differed from single {}",
                batch_hsp.evalue,
                single_hsp.evalue
            );
        }
    }

    #[test]
    fn test_blastp_batch_matches_single_query_default_filters() {
        let tmp = tempfile::tempdir().unwrap();
        let base = tmp.path().join("db");
        let subject = b"MKKWLFGFLGVALIVVCSVFGYVSYQKHEGEVFKQNIEKKMPVDQINAHAKSY";
        let query_a = b"MKKWLFGFLGVALIVVCSVFGYVSYQKHEGEVFKQNIEKK";
        let query_b = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";

        let mut builder = BlastDbBuilder::new(DbType::Protein, "db");
        builder.add(SequenceEntry {
            title: "protein_subject".to_string(),
            accession: "protein_subject".to_string(),
            sequence: subject.to_vec(),
            taxid: None,
        });
        builder.add(SequenceEntry {
            title: "low_complexity_subject".to_string(),
            accession: "low_complexity_subject".to_string(),
            sequence: query_b.to_vec(),
            taxid: None,
        });
        builder.write(&base).unwrap();
        let db = BlastDb::open(&base).unwrap();
        let params = SearchParams::blastp()
            .word_size(3)
            .word_threshold(9.0)
            .evalue(1.0e20)
            .sum_stats(false);

        let batch = blastp_batch(&db, &[query_a.as_slice(), query_b.as_slice()], &params);
        let singles = [blastp(&db, query_a, &params), blastp(&db, query_b, &params)];

        assert_eq!(batch.len(), singles.len());
        for (batch_results, single_results) in batch.iter().zip(singles.iter()) {
            assert_eq!(batch_results.len(), single_results.len());
            for (batch_hit, single_hit) in batch_results.iter().zip(single_results.iter()) {
                assert_eq!(batch_hit.subject_oid, single_hit.subject_oid);
                assert_eq!(batch_hit.hsps.len(), single_hit.hsps.len());
                for (batch_hsp, single_hsp) in batch_hit.hsps.iter().zip(single_hit.hsps.iter()) {
                    assert_eq!(batch_hsp.score, single_hsp.score);
                    assert_eq!(batch_hsp.alignment_length, single_hsp.alignment_length);
                    assert!(
                        (batch_hsp.evalue - single_hsp.evalue).abs()
                            <= single_hsp.evalue.abs().max(1.0) * 1.0e-12,
                        "batch evalue {} differed from single {}",
                        batch_hsp.evalue,
                        single_hsp.evalue
                    );
                }
            }
        }
    }

    #[test]
    fn test_search_with_pssm_honors_effective_search_space() {
        let tmp = tempfile::tempdir().unwrap();
        let base = tmp.path().join("db");
        let query = b"MKKWLFGFLGVALIVVCSVFGYVSYQKHEGEVFKQNIEKK";
        let mut builder = BlastDbBuilder::new(DbType::Protein, "db");
        builder.add(SequenceEntry {
            title: "protein_subject".to_string(),
            accession: "protein_subject".to_string(),
            sequence: query.to_vec(),
            taxid: None,
        });
        builder.write(&base).unwrap();
        let db = BlastDb::open(&base).unwrap();

        let pssm = crate::pssm::Pssm::from_sequence(
            &encode_ncbistdaa_sequence(query),
            get_matrix(MatrixType::Blosum62),
        );
        let small_space = SearchParams::blastp()
            .filter_low_complexity(false)
            .comp_adjust(0)
            .effective_search_space(1_000)
            .evalue(1.0e20);
        let large_space = small_space.clone().effective_search_space(1_000_000);

        let small = search_with_pssm(&db, query, &pssm, &small_space);
        let large = search_with_pssm(&db, query, &pssm, &large_space);
        assert_eq!(small.len(), 1);
        assert_eq!(large.len(), 1);
        assert_eq!(small[0].hsps[0].score, large[0].hsps[0].score);
        assert!(
            large[0].hsps[0].evalue > small[0].hsps[0].evalue * 100.0,
            "effective search space should scale PSSM e-values: small={} large={}",
            small[0].hsps[0].evalue,
            large[0].hsps[0].evalue
        );
    }

    #[test]
    fn test_search_with_pssm_preserves_duplicate_accession_oids() {
        let tmp = tempfile::tempdir().unwrap();
        let base = tmp.path().join("db");
        let query = b"MKKWLFGFLG";
        let mut builder = BlastDbBuilder::new(DbType::Protein, "db");
        for title in ["first", "second"] {
            builder.add(SequenceEntry {
                title: title.to_string(),
                accession: "same_accession".to_string(),
                sequence: query.to_vec(),
                taxid: None,
            });
        }
        builder.write(&base).unwrap();
        let db = BlastDb::open(&base).unwrap();

        let pssm = crate::pssm::Pssm::from_sequence(
            &encode_ncbistdaa_sequence(query),
            get_matrix(MatrixType::Blosum62),
        );
        let params = SearchParams::blastp()
            .filter_low_complexity(false)
            .comp_adjust(0)
            .evalue(1.0e20);

        let results = search_with_pssm(&db, query, &pssm, &params);
        let oids: Vec<u32> = results.iter().map(|result| result.subject_oid).collect();
        assert_eq!(oids, vec![1, 0]);
        assert!(results
            .iter()
            .all(|result| result.subject_accession == "same_accession"));
    }

    #[test]
    fn test_search_with_pssm_applies_max_target_seqs() {
        let tmp = tempfile::tempdir().unwrap();
        let base = tmp.path().join("db");
        let query = b"MKKWLFGFLG";
        let mut builder = BlastDbBuilder::new(DbType::Protein, "db");
        for idx in 0..3 {
            builder.add(SequenceEntry {
                title: format!("subject_{idx}"),
                accession: format!("subject_{idx}"),
                sequence: query.to_vec(),
                taxid: None,
            });
        }
        builder.write(&base).unwrap();
        let db = BlastDb::open(&base).unwrap();

        let pssm = crate::pssm::Pssm::from_sequence(
            &encode_ncbistdaa_sequence(query),
            get_matrix(MatrixType::Blosum62),
        );
        let params = SearchParams::blastp()
            .filter_low_complexity(false)
            .comp_adjust(0)
            .evalue(1.0e20)
            .max_target_seqs(2);

        let results = search_with_pssm(&db, query, &pssm, &params);
        assert_eq!(results.len(), 2);
    }

    #[test]
    fn test_search_with_pssm_populates_alignment_fields() {
        let tmp = tempfile::tempdir().unwrap();
        let base = tmp.path().join("db");
        let query = b"MKKWLFGFLG";
        let mut builder = BlastDbBuilder::new(DbType::Protein, "db");
        builder.add(SequenceEntry {
            title: "protein_subject".to_string(),
            accession: "protein_subject".to_string(),
            sequence: query.to_vec(),
            taxid: None,
        });
        builder.write(&base).unwrap();
        let db = BlastDb::open(&base).unwrap();

        let pssm = crate::pssm::Pssm::from_sequence(
            &encode_ncbistdaa_sequence(query),
            get_matrix(MatrixType::Blosum62),
        );
        let params = SearchParams::blastp()
            .filter_low_complexity(false)
            .comp_adjust(0)
            .evalue(1.0e20);

        let results = search_with_pssm(&db, query, &pssm, &params);
        assert_eq!(results.len(), 1);
        let hsp = &results[0].hsps[0];
        assert_eq!(hsp.num_identities, query.len());
        assert_eq!(hsp.alignment_length, query.len());
        assert_eq!(hsp.query_aln, query);
        assert_eq!(hsp.midline, query);
        assert_eq!(hsp.subject_aln, query);
    }

    #[test]
    fn test_pssm_hit_conversion_groups_multiple_hsps_by_oid() {
        let tmp = tempfile::tempdir().unwrap();
        let base = tmp.path().join("db");
        let query = b"ACDACD";
        let subject = b"ACDYYACD";
        let mut builder = BlastDbBuilder::new(DbType::Protein, "db");
        builder.add(SequenceEntry {
            title: "protein_subject".to_string(),
            accession: "protein_subject".to_string(),
            sequence: subject.to_vec(),
            taxid: None,
        });
        builder.write(&base).unwrap();
        let db = BlastDb::open(&base).unwrap();

        let query_aa = encode_ncbistdaa_sequence(query);
        let pssm = crate::pssm::Pssm::from_sequence(&query_aa, get_matrix(MatrixType::Blosum62));
        let prot_kbp = protein_kbp_for_matrix(MatrixType::Blosum62, 11, 1);
        let hits = vec![
            crate::pssm::PsiBlastHit {
                subject_id: "0".to_string(),
                score: 18,
                bit_score: 123.0,
                evalue: 1.0e-6,
                query_start: 0,
                query_end: 3,
                subject_start: 0,
                subject_end: 3,
                align_len: 3,
                subject_len: subject.len(),
                query_aln: b"ACD".to_vec(),
                subject_aln: b"ACD".to_vec(),
            },
            crate::pssm::PsiBlastHit {
                subject_id: "0".to_string(),
                score: 18,
                bit_score: 456.0,
                evalue: 1.0e-6,
                query_start: 3,
                query_end: 6,
                subject_start: 5,
                subject_end: 8,
                align_len: 3,
                subject_len: subject.len(),
                query_aln: b"ACD".to_vec(),
                subject_aln: b"ACD".to_vec(),
            },
        ];
        let params = SearchParams::blastp()
            .filter_low_complexity(false)
            .comp_adjust(0)
            .evalue(1.0e20);

        let results = pssm_hits_to_search_results(&db, &hits, &query_aa, &pssm, &prot_kbp, &params);

        assert_eq!(results.len(), 1);
        assert_eq!(results[0].subject_oid, 0);
        assert_eq!(results[0].hsps.len(), 2);
        assert!(results[0]
            .hsps
            .iter()
            .any(|hsp| (hsp.query_start, hsp.subject_start) == (0, 0)));
        assert!(results[0]
            .hsps
            .iter()
            .any(|hsp| (hsp.query_start, hsp.subject_start) == (3, 5)));
        assert!(results[0]
            .hsps
            .iter()
            .any(|hsp| (hsp.query_start, hsp.bit_score) == (0, 123.0)));
        assert!(results[0]
            .hsps
            .iter()
            .any(|hsp| (hsp.query_start, hsp.bit_score) == (3, 456.0)));
    }

    #[test]
    fn test_search_with_pssm_preserves_query_offset_for_local_hit() {
        let tmp = tempfile::tempdir().unwrap();
        let base = tmp.path().join("db");
        let query = b"ACDEFGHIK";
        let subject = b"DEFGHIK";
        let mut builder = BlastDbBuilder::new(DbType::Protein, "db");
        builder.add(SequenceEntry {
            title: "protein_subject".to_string(),
            accession: "protein_subject".to_string(),
            sequence: subject.to_vec(),
            taxid: None,
        });
        builder.write(&base).unwrap();
        let db = BlastDb::open(&base).unwrap();

        let pssm = crate::pssm::Pssm::from_sequence(
            &encode_ncbistdaa_sequence(query),
            get_matrix(MatrixType::Blosum62),
        );
        let params = SearchParams::blastp()
            .filter_low_complexity(false)
            .comp_adjust(0)
            .evalue(1.0e20);

        let results = search_with_pssm(&db, query, &pssm, &params);
        assert_eq!(results.len(), 1);
        let hsp = &results[0].hsps[0];
        assert_eq!(hsp.query_start, 2);
        assert_eq!(hsp.query_end, query.len());
        assert_eq!(hsp.subject_start, 0);
        assert_eq!(hsp.subject_end, subject.len());
        assert_eq!(hsp.query_aln, subject);
        assert_eq!(hsp.subject_aln, subject);
    }

    #[test]
    fn test_search_with_pssm_midline_uses_absolute_query_column() {
        let tmp = tempfile::tempdir().unwrap();
        let base = tmp.path().join("db");
        let query = b"AAA";
        let subject = b"EEE";
        let query_aa = encode_ncbistdaa_sequence(query);
        let subject_aa = encode_ncbistdaa_sequence(subject);
        let mut builder = BlastDbBuilder::new(DbType::Protein, "db");
        builder.add(SequenceEntry {
            title: "protein_subject".to_string(),
            accession: "protein_subject".to_string(),
            sequence: subject.to_vec(),
            taxid: None,
        });
        builder.write(&base).unwrap();
        let db = BlastDb::open(&base).unwrap();

        let mut pssm = crate::pssm::Pssm {
            scores: vec![[-10; AA_SIZE]; query_aa.len()],
            length: query_aa.len(),
            info_content: vec![0.0; query_aa.len()],
            start_numerator: None,
            ancillary_gap_kbp: None,
        };
        pssm.scores[2][subject_aa[0] as usize] = 12;
        let params = SearchParams::blastp()
            .filter_low_complexity(false)
            .comp_adjust(0)
            .evalue(1.0e20);

        let results = search_with_pssm(&db, query, &pssm, &params);
        assert_eq!(results.len(), 1);
        let hsp = &results[0].hsps[0];
        assert_eq!(hsp.query_start, 2);
        assert_eq!(hsp.subject_start, 0);
        assert_eq!(hsp.query_aln, b"A");
        assert_eq!(hsp.subject_aln, b"E");
        assert_eq!(hsp.midline, b"+");
    }

    #[test]
    fn test_search_with_pssm_reports_gapped_alignment() {
        let tmp = tempfile::tempdir().unwrap();
        let base = tmp.path().join("db");
        let query = b"ACD";
        let subject = b"AWCD";
        let query_aa = encode_ncbistdaa_sequence(query);
        let mut builder = BlastDbBuilder::new(DbType::Protein, "db");
        builder.add(SequenceEntry {
            title: "protein_subject".to_string(),
            accession: "protein_subject".to_string(),
            sequence: subject.to_vec(),
            taxid: None,
        });
        builder.write(&base).unwrap();
        let db = BlastDb::open(&base).unwrap();

        let mut pssm = crate::pssm::Pssm {
            scores: vec![[-20; AA_SIZE]; query_aa.len()],
            length: query_aa.len(),
            info_content: vec![0.0; query_aa.len()],
            start_numerator: None,
            ancillary_gap_kbp: None,
        };
        for (pos, &aa) in query_aa.iter().enumerate() {
            pssm.scores[pos][aa as usize] = 10;
        }
        let params = SearchParams::blastp()
            .filter_low_complexity(false)
            .comp_adjust(0)
            .gap_open(5)
            .gap_extend(1)
            .evalue(1.0e20);

        let results = search_with_pssm(&db, query, &pssm, &params);

        assert_eq!(results.len(), 1);
        let hsp = &results[0].hsps[0];
        assert_eq!(hsp.query_start, 0);
        assert_eq!(hsp.query_end, 3);
        assert_eq!(hsp.subject_start, 0);
        assert_eq!(hsp.subject_end, 4);
        assert_eq!(hsp.query_aln, b"A-CD");
        assert_eq!(hsp.subject_aln, b"AWCD");
        assert_eq!(hsp.num_gaps, 1);
        assert_eq!(hsp.alignment_length, 4);
    }

    #[test]
    fn test_search_with_pssm_allows_single_residue_query() {
        let tmp = tempfile::tempdir().unwrap();
        let base = tmp.path().join("db");
        let query = b"A";
        let query_aa = encode_ncbistdaa_sequence(query);
        let mut builder = BlastDbBuilder::new(DbType::Protein, "db");
        builder.add(SequenceEntry {
            title: "protein_subject".to_string(),
            accession: "protein_subject".to_string(),
            sequence: query.to_vec(),
            taxid: None,
        });
        builder.write(&base).unwrap();
        let db = BlastDb::open(&base).unwrap();

        let mut pssm = crate::pssm::Pssm {
            scores: vec![[-20; AA_SIZE]; 1],
            length: 1,
            info_content: vec![0.0],
            start_numerator: None,
            ancillary_gap_kbp: None,
        };
        pssm.scores[0][query_aa[0] as usize] = 8;
        let params = SearchParams::blastp()
            .filter_low_complexity(false)
            .comp_adjust(0)
            .evalue(1.0e20);

        let results = search_with_pssm(&db, query, &pssm, &params);

        assert_eq!(results.len(), 1);
        let hsp = &results[0].hsps[0];
        assert_eq!(hsp.query_start, 0);
        assert_eq!(hsp.query_end, 1);
        assert_eq!(hsp.subject_start, 0);
        assert_eq!(hsp.subject_end, 1);
        assert_eq!(hsp.query_aln, b"A");
        assert_eq!(hsp.subject_aln, b"A");
        assert_eq!(hsp.midline, b"A");
    }

    #[test]
    fn test_build_pssm_places_partial_hsp_at_query_offset() {
        let query = b"ACDEFGHIK";
        let subject_aln = b"DEFGHIK".to_vec();
        let result = SearchResult {
            subject_oid: 0,
            subject_title: "subject".to_string(),
            subject_accession: "subject".to_string(),
            subject_len: subject_aln.len(),
            hsps: vec![Hsp {
                score: 42,
                bit_score: 42.0,
                evalue: 1.0e-20,
                query_start: 2,
                query_end: query.len(),
                subject_start: 0,
                subject_end: subject_aln.len(),
                query_gapped_start: 2,
                subject_gapped_start: 0,
                query_link_start: 2,
                query_link_end: query.len(),
                query_link_gapped_start: 2,
                subject_link_start: 0,
                subject_link_end: subject_aln.len(),
                subject_link_gapped_start: 0,
                link_score: None,
                link_lambda: None,
                num_identities: subject_aln.len(),
                num_gaps: 0,
                alignment_length: subject_aln.len(),
                query_aln: subject_aln.clone(),
                midline: subject_aln.clone(),
                subject_aln,
                query_frame: 0,
                subject_frame: 0,
                num_links: 1,
                comp_adjust_method: 0,
            }],
            taxids: Vec::new(),
        };

        let pssm = build_pssm(query, &[result], 1.0e-3, MatrixType::Blosum62, 0.0);
        let encoded_query = encode_ncbistdaa_sequence(query);

        assert!(
            pssm.score_at(2, encoded_query[2]) > pssm.score_at(0, encoded_query[2]),
            "partial HSP residue should reinforce its true query column, not column zero"
        );
    }

    #[test]
    fn test_build_pssm_projects_gapped_alignment_to_query_columns() {
        let hsp = Hsp {
            score: 42,
            bit_score: 42.0,
            evalue: 1.0e-20,
            query_start: 0,
            query_end: 3,
            subject_start: 0,
            subject_end: 4,
            query_gapped_start: 0,
            subject_gapped_start: 0,
            query_link_start: 0,
            query_link_end: 3,
            query_link_gapped_start: 0,
            subject_link_start: 0,
            subject_link_end: 4,
            subject_link_gapped_start: 0,
            link_score: None,
            link_lambda: None,
            num_identities: 3,
            num_gaps: 1,
            alignment_length: 4,
            query_aln: b"A-CD".to_vec(),
            midline: b"A CD".to_vec(),
            subject_aln: b"AWCD".to_vec(),
            query_frame: 0,
            subject_frame: 0,
            num_links: 1,
            comp_adjust_method: 0,
        };

        let projected = project_subject_alignment_to_query(&hsp, 3).unwrap();

        assert_eq!(
            projected,
            encode_ncbistdaa_sequence(b"ACD"),
            "subject insertions relative to the query must not shift PSSM columns"
        );
    }

    #[test]
    fn test_build_pssm_preserves_subject_deletion_as_gap_column() {
        let hsp = Hsp {
            score: 42,
            bit_score: 42.0,
            evalue: 1.0e-20,
            query_start: 1,
            query_end: 4,
            subject_start: 0,
            subject_end: 2,
            query_gapped_start: 1,
            subject_gapped_start: 0,
            query_link_start: 1,
            query_link_end: 4,
            query_link_gapped_start: 1,
            subject_link_start: 0,
            subject_link_end: 2,
            subject_link_gapped_start: 0,
            link_score: None,
            link_lambda: None,
            num_identities: 2,
            num_gaps: 1,
            alignment_length: 3,
            query_aln: b"CDE".to_vec(),
            midline: b"C E".to_vec(),
            subject_aln: b"C-E".to_vec(),
            query_frame: 0,
            subject_frame: 0,
            num_links: 1,
            comp_adjust_method: 0,
        };

        let projected = project_subject_alignment_to_query(&hsp, 5).unwrap();

        assert_eq!(projected[0], NCBISTDAA_X);
        assert_eq!(projected[1], crate::encoding::NCBISTDAA_C);
        assert_eq!(projected[2], NCBISTDAA_GAP);
        assert_eq!(projected[3], crate::encoding::NCBISTDAA_E);
        assert_eq!(projected[4], NCBISTDAA_X);
    }

    #[test]
    fn test_psiblast_update_projects_gapped_hit_alignment_to_query_columns() {
        let hit = crate::pssm::PsiBlastHit {
            subject_id: "0".to_string(),
            score: 42,
            bit_score: 0.0,
            evalue: 1.0e-20,
            query_start: 0,
            query_end: 3,
            subject_start: 0,
            subject_end: 4,
            align_len: 3,
            subject_len: 4,
            query_aln: b"A-CD".to_vec(),
            subject_aln: b"AWCD".to_vec(),
        };

        let projected = project_psiblast_hit_alignment_to_query(&hit, 3).unwrap();

        assert_eq!(
            projected,
            encode_ncbistdaa_sequence(b"ACD"),
            "subject insertions from PSI-BLAST hits must not shift next-iteration PSSM columns"
        );
    }

    #[test]
    fn test_psiblast_update_preserves_hit_subject_deletion_as_gap() {
        let hit = crate::pssm::PsiBlastHit {
            subject_id: "0".to_string(),
            score: 42,
            bit_score: 0.0,
            evalue: 1.0e-20,
            query_start: 1,
            query_end: 4,
            subject_start: 0,
            subject_end: 2,
            align_len: 3,
            subject_len: 2,
            query_aln: b"CDE".to_vec(),
            subject_aln: b"C-E".to_vec(),
        };

        let projected = project_psiblast_hit_alignment_to_query(&hit, 5).unwrap();

        assert_eq!(projected[0], NCBISTDAA_X);
        assert_eq!(projected[1], crate::encoding::NCBISTDAA_C);
        assert_eq!(projected[2], NCBISTDAA_GAP);
        assert_eq!(projected[3], crate::encoding::NCBISTDAA_E);
        assert_eq!(projected[4], NCBISTDAA_X);
    }

    #[test]
    fn test_psiblast_reports_hits_above_inclusion_threshold() {
        let tmp = tempfile::tempdir().unwrap();
        let base = tmp.path().join("db");
        let query = b"MKKWLFGFLG";
        let mut builder = BlastDbBuilder::new(DbType::Protein, "db");
        builder.add(SequenceEntry {
            title: "protein_subject".to_string(),
            accession: "protein_subject".to_string(),
            sequence: query.to_vec(),
            taxid: None,
        });
        builder.write(&base).unwrap();
        let db = BlastDb::open(&base).unwrap();

        let search = SearchParams::blastp()
            .filter_low_complexity(false)
            .comp_adjust(0)
            .evalue(1.0e20);
        let params = PsiblastParams::new(search)
            .num_iterations(1)
            .inclusion_evalue(1.0e-100);

        let (results, _) = psiblast(&db, query, &params);
        assert_eq!(results.len(), 1);
        assert!(
            results[0].hsps[0].evalue > params.inclusion_evalue,
            "reported hit should be controlled by search e-value, not inclusion e-value"
        );
    }

    #[test]
    fn test_psiblast_params_default_inclusion_threshold_matches_ncbi() {
        let params = PsiblastParams::new(SearchParams::blastp());

        assert_eq!(params.inclusion_evalue, crate::stat::PSI_INCLUSION_ETHRESH);
    }

    #[test]
    fn test_psiblast_params_pseudocount_matches_ncbi_default_and_setter() {
        let params = PsiblastParams::new(SearchParams::blastp());

        assert_eq!(params.pseudocount, crate::stat::PSI_PSEUDO_COUNT_CONST);
        assert_eq!(params.pseudocount(7).pseudocount, 7);
    }

    #[test]
    fn test_psiblast_params_gap_trigger_matches_ncbi_default_and_setter() {
        let params = PsiblastParams::new(SearchParams::blastp());

        assert_eq!(params.gap_trigger, crate::stat::BLAST_GAP_TRIGGER_PROT);
        assert_eq!(params.gap_trigger(18.0).gap_trigger, 18.0);
    }

    #[test]
    fn test_psiblast_restart_alignment_updates_initial_pssm() {
        let tmp = tempfile::tempdir().unwrap();
        let base = tmp.path().join("db");
        let query = b"ACD";
        let mut builder = BlastDbBuilder::new(DbType::Protein, "db");
        builder.add(SequenceEntry {
            title: "protein_subject".to_string(),
            accession: "protein_subject".to_string(),
            sequence: query.to_vec(),
            taxid: None,
        });
        builder.write(&base).unwrap();
        let db = BlastDb::open(&base).unwrap();

        let query_aa = encode_ncbistdaa_sequence(query);
        let restart = vec![
            encode_ncbistdaa_sequence(b"ACD"),
            encode_ncbistdaa_sequence(b"ADD"),
            encode_ncbistdaa_sequence(b"AED"),
        ];
        let baseline =
            crate::pssm::Pssm::from_sequence(&query_aa, get_matrix(MatrixType::Blosum62));
        let params = PsiblastParams::new(SearchParams::blastp())
            .num_iterations(0)
            .restart_alignment(restart);

        let (_, pssm) = psiblast(&db, query, &params);

        assert_ne!(
            pssm.scores, baseline.scores,
            "restart alignment should update the initial PSI-BLAST PSSM before iteration"
        );
    }

    #[test]
    fn test_psiblast_with_rounds_collects_updated_pssms() {
        let tmp = tempfile::tempdir().unwrap();
        let base = tmp.path().join("db");
        let query = b"MKKWLFGFLG";
        let mut builder = BlastDbBuilder::new(DbType::Protein, "db");
        builder.add(SequenceEntry {
            title: "protein_subject".to_string(),
            accession: "protein_subject".to_string(),
            sequence: query.to_vec(),
            taxid: None,
        });
        builder.write(&base).unwrap();
        let db = BlastDb::open(&base).unwrap();

        let params = PsiblastParams::new(
            SearchParams::blastp()
                .filter_low_complexity(false)
                .comp_adjust(0)
                .evalue(1.0e20),
        )
        .num_iterations(1);

        let run = psiblast_with_rounds(&db, query, &params);

        assert_eq!(run.round_pssms.len(), 1);
        assert_eq!(run.round_results.len(), 1);
        assert_eq!(run.round_pssms[0].length, query.len());
        assert_eq!(run.pssm.scores, run.round_pssms[0].scores);
        assert!(!run.converged);
    }

    #[test]
    fn test_psiblast_num_iterations_zero_runs_until_convergence() {
        let tmp = tempfile::tempdir().unwrap();
        let base = tmp.path().join("db");
        let query = b"MKKWLFGFLG";
        let mut builder = BlastDbBuilder::new(DbType::Protein, "db");
        builder.add(SequenceEntry {
            title: "protein_subject".to_string(),
            accession: "protein_subject".to_string(),
            sequence: query.to_vec(),
            taxid: None,
        });
        builder.write(&base).unwrap();
        let db = BlastDb::open(&base).unwrap();

        let params = PsiblastParams::new(
            SearchParams::blastp()
                .filter_low_complexity(false)
                .comp_adjust(0)
                .evalue(1.0e20),
        )
        .num_iterations(0);

        let run = psiblast_with_rounds(&db, query, &params);

        assert!(
            !run.results.is_empty(),
            "num_iterations=0 should not degenerate into a zero-round search"
        );
        assert_eq!(
            run.round_pssms.len(),
            2,
            "exact repeat fixture should converge after the second update"
        );
        assert_eq!(run.round_results.len(), 2);
        assert!(run.converged);
        assert_eq!(run.round_pssms[0].length, query.len());
        assert_eq!(run.pssm.scores, run.round_pssms[1].scores);
    }
}
