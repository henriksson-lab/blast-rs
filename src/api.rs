//! High-level public API for running BLAST searches.
//!
//! Provides `blastp`, `blastn`, `blastx`, `tblastn`, and other search functions
//! with a simple builder-pattern interface, matching the API of the previous
//! blast-rs implementation.

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;
use std::sync::Arc;

use crate::db::{BlastDb, DbType};
use crate::encoding::{AMINOACID_TO_NCBISTDAA, IUPACNA_TO_BLASTNA, NCBISTDAA_TO_AMINOACID};
use crate::matrix::AA_SIZE;
use crate::search::{blastn_gapped_search_nomask, SearchHsp};
use crate::stat::{
    compute_length_adjustment_exact, nucl_alpha_beta, nucl_gapped_kbp_lookup, ungapped_kbp_calc,
    GappedParams, GumbelBlk, KarlinBlk, UngappedKbpContext,
};
use crate::traceback::build_blastna_matrix;

/// Re-export of `hspstream::evalue_comp` (port of NCBI `s_EvalueComp`):
/// compare two e-values, treating both below `1e-180` as equal.
pub use crate::hspstream::evalue_comp;

const COMPO_ADJUST_SCALE_FACTOR: f64 = 32.0;

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
    pub num_identities: usize,
    pub num_gaps: usize,
    pub alignment_length: usize,
    pub query_aln: Vec<u8>,
    pub midline: Vec<u8>,
    pub subject_aln: Vec<u8>,
    pub query_frame: i32,
    pub subject_frame: i32,
}

impl Hsp {
    /// Percent identity of this HSP.
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
    pub fn best_evalue(&self) -> f64 {
        self.hsps
            .iter()
            .map(|h| h.evalue)
            .fold(f64::INFINITY, f64::min)
    }

    /// Top HSP score, or `i32::MIN` if the result has no HSPs. Used as the
    /// secondary key in `compare_search_results` (NCBI `s_EvalueCompareHSPLists`).
    fn top_score(&self) -> i32 {
        self.hsps.iter().map(|h| h.score).max().unwrap_or(i32::MIN)
    }
}

/// Port of NCBI `s_EvalueCompareHSPLists` (`blast_hits.c:3078`) for
/// `SearchResult`. Primary key is `best_evalue` (compared via
/// `crate::stat::evalue_comp` — both <1e-180 compare equal). Ties break on
/// top HSP score descending, then `subject_oid` descending (NCBI's
/// `BLAST_CMP(h2->oid, h1->oid)`).
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
    let by_score = b.top_score().cmp(&a.top_score());
    if by_score != Equal {
        return by_score;
    }
    b.subject_oid.cmp(&a.subject_oid)
}

fn map_database_oids<T, F>(db: &BlastDb, params: &SearchParams, f: F) -> Vec<T>
where
    T: Send,
    F: Fn(u32) -> T + Sync + Send,
{
    if params.thread_pool.is_none() && params.num_threads == 1 {
        (0..db.num_oids).map(f).collect()
    } else if let Some(pool) = params.thread_pool.as_deref() {
        use rayon::prelude::*;
        pool.install(|| (0..db.num_oids).into_par_iter().map(&f).collect())
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
        pool.install(|| (0..db.num_oids).into_par_iter().map(f).collect())
    }
}

fn nearly_equal_evalues(a: f64, b: f64) -> bool {
    if a == b {
        return true;
    }
    let scale = a.abs().max(b.abs()).max(1.0e-300);
    ((a - b).abs() / scale) <= 0.02
}

fn compare_tblastx_hsps(a: &Hsp, b: &Hsp) -> std::cmp::Ordering {
    let rendered_equal =
        crate::format::format_evalue(a.evalue) == crate::format::format_evalue(b.evalue);
    let by_evalue =
        if a.score == b.score || rendered_equal || nearly_equal_evalues(a.evalue, b.evalue) {
            std::cmp::Ordering::Equal
        } else {
            crate::hspstream::evalue_comp(a.evalue, b.evalue)
        };
    by_evalue
        .then_with(|| b.score.cmp(&a.score))
        .then_with(|| a.query_frame.abs().cmp(&b.query_frame.abs()))
        .then_with(|| b.query_frame.cmp(&a.query_frame))
        .then_with(|| a.subject_frame.abs().cmp(&b.subject_frame.abs()))
        .then_with(|| b.subject_frame.cmp(&a.subject_frame))
}

struct TranslatedContextStats {
    frame: i32,
    query_length: i32,
    eff_searchsp: i64,
    length_adjustment: i32,
    /// Per-context ungapped Karlin params. NCBI computes one Karlin block per
    /// query frame from the frame's translation composition; e-value sites
    /// must match the bit-score path, which already uses these.
    kbp: KarlinBlk,
}

fn apply_tblastx_linked_sum_stats(
    results: &mut [Option<SearchResult>],
    query_contexts: &[TranslatedContextStats],
) {
    use crate::link_hsps::{
        BLAST_LinkHsps, LinkBlastHsp, LinkBlastHspList, LinkBlastSeg, LinkHSPParameters,
        LinkScoreBlock,
    };
    use crate::program::TBLASTX;
    use crate::queryinfo::{ContextInfo, QueryInfo};

    if query_contexts.is_empty() {
        return;
    }

    fn translated_coord_to_protein(coord: usize, frame: i32) -> i32 {
        let offset = frame.unsigned_abs() as usize - 1;
        ((coord.saturating_sub(offset)) / 3) as i32
    }

    let mut contexts = Vec::with_capacity(query_contexts.len());
    let mut offset = 0i32;
    for ctx in query_contexts {
        contexts.push(ContextInfo {
            query_offset: offset,
            query_length: ctx.query_length,
            eff_searchsp: ctx.eff_searchsp,
            length_adjustment: ctx.length_adjustment,
            query_index: 0,
            frame: ctx.frame,
            is_valid: true,
        });
        offset += ctx.query_length + 1;
    }
    let query_info = QueryInfo {
        num_queries: 1,
        max_length: contexts
            .iter()
            .map(|ctx| ctx.query_length.max(0) as u32)
            .max()
            .unwrap_or(0),
        contexts,
    };
    let kbps: Vec<KarlinBlk> = query_contexts.iter().map(|c| c.kbp.clone()).collect();
    let score_block = LinkScoreBlock {
        kbp: kbps.clone(),
        kbp_gap: kbps,
    };
    let link_params = LinkHSPParameters {
        gap_prob: 1.0,
        gap_decay_rate: crate::stat::BLAST_GAP_DECAY_RATE_GAPPED,
        cutoff_small_gap: 0,
        cutoff_big_gap: 0,
        longest_intron: (122 - 2) / 3,
        ..LinkHSPParameters::default()
    };

    for (oid, result) in results.iter_mut().enumerate() {
        let Some(result) = result.as_mut() else {
            continue;
        };
        if result.hsps.len() <= 1 {
            continue;
        }

        let mut hsp_list = LinkBlastHspList {
            oid: oid as i32,
            query_index: 0,
            hsp_array: result
                .hsps
                .iter()
                .map(|hsp| {
                    let context = query_contexts
                        .iter()
                        .position(|ctx| ctx.frame == hsp.query_frame)
                        .unwrap_or(0) as i32;
                    LinkBlastHsp {
                        score: hsp.score,
                        num_ident: hsp.num_identities as i32,
                        bit_score: hsp.bit_score,
                        evalue: hsp.evalue,
                        query: LinkBlastSeg {
                            frame: hsp.query_frame,
                            offset: translated_coord_to_protein(hsp.query_start, hsp.query_frame),
                            end: translated_coord_to_protein(hsp.query_end, hsp.query_frame),
                            gapped_start: translated_coord_to_protein(
                                hsp.query_start,
                                hsp.query_frame,
                            ),
                        },
                        subject: LinkBlastSeg {
                            frame: hsp.subject_frame,
                            offset: translated_coord_to_protein(
                                hsp.subject_start,
                                hsp.subject_frame,
                            ),
                            end: translated_coord_to_protein(hsp.subject_end, hsp.subject_frame),
                            gapped_start: translated_coord_to_protein(
                                hsp.subject_start,
                                hsp.subject_frame,
                            ),
                        },
                        context,
                        num: 1,
                    }
                })
                .collect(),
            best_evalue: result.best_evalue(),
        };
        let original_keys: Vec<(i32, i32, i32, i32, i32)> = result
            .hsps
            .iter()
            .map(|hsp| {
                (
                    hsp.score,
                    hsp.query_frame,
                    hsp.subject_frame,
                    translated_coord_to_protein(hsp.query_start, hsp.query_frame),
                    translated_coord_to_protein(hsp.subject_start, hsp.subject_frame),
                )
            })
            .collect();

        // tblastx is ungapped (`blast_options.c:869`); pass
        // `gapped_calculation=false` so `BLAST_LinkHsps` reads `kbp` (ungapped),
        // matching `link_hsps.c:457`. The `score_block` already populates both
        // arrays with ungapped Karlin params for this reason.
        BLAST_LinkHsps(
            TBLASTX,
            &mut hsp_list,
            &query_info,
            result.subject_len as i32,
            &score_block,
            &link_params,
            false,
        );

        let mut linked_evalues: std::collections::HashMap<(i32, i32, i32, i32, i32), Vec<f64>> =
            std::collections::HashMap::new();
        for linked in &hsp_list.hsp_array {
            linked_evalues
                .entry((
                    linked.score,
                    linked.query.frame,
                    linked.subject.frame,
                    linked.query.offset,
                    linked.subject.offset,
                ))
                .or_default()
                .push(linked.evalue);
        }

        for (hsp, key) in result.hsps.iter_mut().zip(original_keys) {
            if let Some(evalues) = linked_evalues.get_mut(&key) {
                if let Some(evalue) = evalues.pop() {
                    hsp.evalue = evalue;
                }
            }
        }
    }
}

fn apply_tblastn_linked_sum_stats(
    results: &mut [SearchResult],
    query_info: &crate::queryinfo::QueryInfo,
    prot_kbp: &KarlinBlk,
) {
    use crate::link_hsps::{
        BLAST_LinkHsps, LinkBlastHsp, LinkBlastHspList, LinkBlastSeg, LinkHSPParameters,
        LinkScoreBlock,
    };
    use crate::program::TBLASTN;

    fn nuc_coord_to_protein(coord: usize, frame: i32) -> i32 {
        let offset = frame.unsigned_abs() as usize - 1;
        ((coord.saturating_sub(offset)) / 3) as i32
    }

    let score_block = LinkScoreBlock {
        kbp: vec![prot_kbp.clone()],
        kbp_gap: vec![prot_kbp.clone()],
    };
    let link_params = LinkHSPParameters::default();

    for (oid, result) in results.iter_mut().enumerate() {
        if result.hsps.is_empty() {
            continue;
        }
        let mut hsp_list = LinkBlastHspList {
            oid: oid as i32,
            query_index: 0,
            hsp_array: result
                .hsps
                .iter()
                .map(|hsp| LinkBlastHsp {
                    score: hsp.score,
                    num_ident: hsp.num_identities as i32,
                    bit_score: hsp.bit_score,
                    evalue: hsp.evalue,
                    query: LinkBlastSeg {
                        frame: 0,
                        offset: hsp.query_start as i32,
                        end: hsp.query_end as i32,
                        gapped_start: hsp.query_start as i32,
                    },
                    subject: LinkBlastSeg {
                        frame: hsp.subject_frame,
                        offset: nuc_coord_to_protein(hsp.subject_start, hsp.subject_frame),
                        end: nuc_coord_to_protein(hsp.subject_end, hsp.subject_frame),
                        gapped_start: nuc_coord_to_protein(hsp.subject_start, hsp.subject_frame),
                    },
                    context: 0,
                    num: 1,
                })
                .collect(),
            best_evalue: result.best_evalue(),
        };
        let original_keys: Vec<(i32, i32, i32, i32)> = result
            .hsps
            .iter()
            .map(|hsp| {
                (
                    hsp.score,
                    hsp.subject_frame,
                    hsp.query_start as i32,
                    nuc_coord_to_protein(hsp.subject_start, hsp.subject_frame),
                )
            })
            .collect();

        BLAST_LinkHsps(
            TBLASTN,
            &mut hsp_list,
            query_info,
            result.subject_len as i32,
            &score_block,
            &link_params,
            true,
        );

        let mut linked_evalues: std::collections::HashMap<(i32, i32, i32, i32), Vec<f64>> =
            std::collections::HashMap::new();
        for linked in &hsp_list.hsp_array {
            linked_evalues
                .entry((
                    linked.score,
                    linked.subject.frame,
                    linked.query.offset,
                    linked.subject.offset,
                ))
                .or_default()
                .push(linked.evalue);
        }

        for (hsp, key) in result.hsps.iter_mut().zip(original_keys) {
            if let Some(evalues) = linked_evalues.get_mut(&key) {
                if let Some(evalue) = evalues.pop() {
                    hsp.evalue = evalue;
                }
            }
        }
    }
}

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
    use crate::link_hsps::{
        BLAST_LinkHsps, LinkBlastHsp, LinkBlastHspList, LinkBlastSeg, LinkHSPParameters,
        LinkScoreBlock,
    };
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
    };
    let link_params = LinkHSPParameters::default();
    let mut hsp_list = LinkBlastHspList {
        oid: 0,
        query_index: 0,
        hsp_array: hsps
            .iter()
            .map(|hsp| LinkBlastHsp {
                score: hsp.score,
                num_ident: hsp.num_ident,
                bit_score: hsp.bit_score,
                evalue: hsp.evalue,
                query: LinkBlastSeg {
                    frame: 1,
                    offset: hsp.query_start,
                    end: hsp.query_end,
                    gapped_start: hsp.query_start,
                },
                subject: LinkBlastSeg {
                    frame: if hsp.context == 1 { -1 } else { 1 },
                    offset: hsp.subject_start,
                    end: hsp.subject_end,
                    gapped_start: hsp.subject_start,
                },
                context: hsp.context,
                num: 1,
            })
            .collect(),
        best_evalue: f64::INFINITY,
    };
    let original_keys: Vec<(i32, i32, i32, i32, i32, i32)> = hsps
        .iter()
        .map(|hsp| {
            (
                hsp.score,
                hsp.context,
                hsp.query_start,
                hsp.query_end,
                hsp.subject_start,
                hsp.subject_end,
            )
        })
        .collect();

    BLAST_LinkHsps(
        BLASTN,
        &mut hsp_list,
        &query_info,
        subject_len,
        &score_block,
        &link_params,
        false,
    );

    let mut linked_evalues: HashMap<(i32, i32, i32, i32, i32, i32), Vec<f64>> = HashMap::new();
    for linked in &hsp_list.hsp_array {
        linked_evalues
            .entry((
                linked.score,
                linked.context,
                linked.query.offset,
                linked.query.end,
                linked.subject.offset,
                linked.subject.end,
            ))
            .or_default()
            .push(linked.evalue);
    }

    for (hsp, key) in hsps.iter_mut().zip(original_keys) {
        if let Some(evalues) = linked_evalues.get_mut(&key) {
            if let Some(evalue) = evalues.pop() {
                hsp.evalue = evalue;
            }
        }
    }
}

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

fn prune_translated_hsp_variants(hsps: &mut Vec<Hsp>) {
    let mut best_by_start: HashMap<(usize, usize, i32, i32), Hsp> = HashMap::new();
    for hsp in hsps.drain(..) {
        let key = (
            hsp.query_start,
            hsp.subject_start,
            hsp.query_frame,
            hsp.subject_frame,
        );
        match best_by_start.get_mut(&key) {
            Some(existing) => {
                let hsp_span = (
                    hsp.score,
                    usize::MAX
                        - (hsp.query_end.saturating_sub(hsp.query_start)
                            + hsp.subject_end.saturating_sub(hsp.subject_start)),
                    usize::MAX - hsp.query_end,
                    usize::MAX - hsp.subject_end,
                );
                let existing_span = (
                    existing.score,
                    usize::MAX
                        - (existing.query_end.saturating_sub(existing.query_start)
                            + existing.subject_end.saturating_sub(existing.subject_start)),
                    usize::MAX - existing.query_end,
                    usize::MAX - existing.subject_end,
                );
                if hsp_span > existing_span {
                    *existing = hsp;
                }
            }
            None => {
                best_by_start.insert(key, hsp);
            }
        }
    }
    hsps.extend(best_by_start.into_values());
}

fn compare_hsps_by_evalue_then_coords(a: &Hsp, b: &Hsp) -> std::cmp::Ordering {
    crate::hspstream::evalue_comp(a.evalue, b.evalue)
        .then_with(|| a.query_start.cmp(&b.query_start))
        .then_with(|| a.subject_start.cmp(&b.subject_start))
        .then_with(|| a.query_end.cmp(&b.query_end))
        .then_with(|| a.subject_end.cmp(&b.subject_end))
}

fn protein_eval_cutoff(
    evalue_threshold: f64,
    prot_kbp: &crate::stat::KarlinBlk,
    gumbel_blk: Option<&crate::stat::GumbelBlk>,
    query_length: i32,
    min_subject_length: i32,
    search_space: f64,
) -> i32 {
    if let Some(gbp) = gumbel_blk {
        crate::stat::spouge_etos(
            evalue_threshold,
            prot_kbp,
            gbp,
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
    if subject_span < query_len.min(MINIMUM_LENGTH_NEAR_IDENTICAL) {
        return false;
    }
    let align_len = query_span.min(subject_span);
    if align_len == 0 || subject_len == 0 {
        return false;
    }
    let cutoff = (NEAR_IDENTICAL_BITS_PER_POSITION * crate::math::NCBIMATH_LN2) / gapped_lambda;
    (ph.score as f64 / align_len as f64) >= cutoff
}

fn kappa_seg_mask_subject_for_redo(subject: &[u8]) -> Vec<u8> {
    let mut masked = subject.to_vec();
    let x = AMINOACID_TO_NCBISTDAA[b'X' as usize & 0x7F];
    let mask = crate::filter::seg_filter_ncbistdaa(&masked, 10, 1.8, 2.1);
    for region in &mask.regions {
        let start = region.start.max(0) as usize;
        let end = (region.end as usize).min(masked.len());
        for aa in &mut masked[start..end] {
            *aa = x;
        }
    }
    masked
}

fn kappa_redo_subject_sequence<'a>(
    query_len: usize,
    subject: &'a [u8],
    phits: &[crate::protein_lookup::ProteinHit],
    gapped_lambda: f64,
) -> std::borrow::Cow<'a, [u8]> {
    let near_identical = phits
        .first()
        .map(|ph| kappa_redo_near_identical(ph, query_len, subject.len(), gapped_lambda))
        .unwrap_or(false);
    if near_identical {
        std::borrow::Cow::Borrowed(subject)
    } else {
        std::borrow::Cow::Owned(kappa_seg_mask_subject_for_redo(subject))
    }
}

fn protein_alignment_hits(
    query_aa: &[u8],
    subj_aa: &[u8],
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
    lookup_table: &crate::protein_lookup::ProteinLookupTable,
    x_drop_ungapped: i32,
    gap_open: i32,
    gap_extend: i32,
    x_drop_gapped: i32,
    x_drop_final: i32,
    gap_trigger_raw: i32,
) -> Vec<crate::protein_lookup::ProteinHit> {
    use crate::itree::{Interval, IntervalTree};

    let ungapped_hits = crate::protein_lookup::protein_scan_with_table(
        query_aa,
        subj_aa,
        matrix,
        lookup_table,
        x_drop_ungapped,
    );
    let mut hits = Vec::new();
    // Caller passes the desired seed cutoff directly. Keep the program-specific
    // cutoff calculation at the call site so blastp, blastx, tblastn, and
    // tblastx can follow their own setup rules.
    let seed_cutoff = gap_trigger_raw.max(1);
    // Mirror NCBI's `BLAST_GetGappedScore` flow (iter 120/121 — seed iteration
    // by score desc, pre-gapped containment check, gapped DP, post-gapped
    // containment, then `Blast_HSPListPurgeHSPsWithCommonEndpoints` after the
    // loop).
    let mut passing_seeds: Vec<&crate::protein_lookup::ProteinHit> = ungapped_hits
        .iter()
        .filter(|uh| uh.score >= seed_cutoff)
        .collect();
    passing_seeds.sort_by(|a, b| b.score.cmp(&a.score));
    let mut tree = IntervalTree::new(query_aa.len() as i32 + 1, subj_aa.len() as i32 + 1);
    for seed in &passing_seeds {
        let pre_contained = tree.is_contained(
            Interval::new(seed.query_start as i32, seed.query_end as i32),
            Interval::new(seed.subject_start as i32, seed.subject_end as i32),
            seed.score,
        );
        if pre_contained {
            continue;
        }
        let (seed_q, seed_s) = crate::protein::get_start_for_gapped_alignment(
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
        let Some(prelim) = crate::protein::protein_gapped_align(
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
        if prelim.score < seed_cutoff {
            continue;
        }
        let post_contained = tree.is_contained(
            Interval::new(prelim.query_start as i32, prelim.query_end as i32),
            Interval::new(prelim.subject_start as i32, prelim.subject_end as i32),
            prelim.score,
        );
        if post_contained {
            continue;
        }
        hits.push(crate::protein_lookup::ProteinHit {
            query_start: prelim.query_start,
            query_end: prelim.query_end,
            subject_start: prelim.subject_start,
            subject_end: prelim.subject_end,
            score: prelim.score,
            num_ident: prelim.num_ident,
            align_length: prelim.align_length,
            mismatches: prelim.mismatches,
            gap_opens: prelim.gap_opens,
            qseq: None,
            sseq: None,
            // Pack seed for post-loop traceback re-run.
            scaled_score: Some(((seed_q as i32) << 16) | (seed_s as i32 & 0xffff)),
            gapped_start_q: seed_q,
            gapped_start_s: seed_s,
        });
        tree.insert(
            Interval::new(prelim.query_start as i32, prelim.query_end as i32),
            Interval::new(prelim.subject_start as i32, prelim.subject_end as i32),
            prelim.score,
        );
    }
    // TRACEBACK phase: re-run gapped DP with `x_drop_final` for each
    // accepted preliminary HSP. Mirrors NCBI's
    // `Blast_TracebackFromHSPList` loop (`blast_traceback.c:375-625`).
    for ph in hits.iter_mut() {
        let Some(packed) = ph.scaled_score.take() else {
            continue;
        };
        let seed_q = (packed >> 16) as usize;
        let seed_s = (packed & 0xffff) as usize;
        let Some(gr) = crate::protein::protein_gapped_align(
            query_aa,
            subj_aa,
            seed_q,
            seed_s,
            matrix,
            gap_open,
            gap_extend,
            x_drop_final,
        ) else {
            continue;
        };
        let q_slice = &query_aa[gr.query_start..gr.query_end];
        let s_slice = &subj_aa[gr.subject_start..gr.subject_end];
        let (qseq, sseq) =
            gr.edit_script
                .render_alignment(q_slice, s_slice, crate::protein::ncbistdaa_to_char);
        ph.query_start = gr.query_start;
        ph.query_end = gr.query_end;
        ph.subject_start = gr.subject_start;
        ph.subject_end = gr.subject_end;
        ph.score = gr.score;
        ph.num_ident = gr.num_ident;
        ph.align_length = gr.align_length;
        ph.mismatches = gr.mismatches;
        ph.gap_opens = gr.gap_opens;
        ph.qseq = Some(qseq);
        ph.sseq = Some(sseq);
    }
    // Mirror NCBI's `Blast_HSPListPurgeHSPsWithCommonEndpoints`
    // (`blast_hits.c:2455`), called after gapped extension at
    // `blast_engine.c:544`.
    purge_hsps_with_common_endpoints(&mut hits);
    hits.sort_by(|a, b| {
        b.score
            .cmp(&a.score)
            .then_with(|| a.subject_start.cmp(&b.subject_start))
            .then_with(|| b.subject_end.cmp(&a.subject_end))
            .then_with(|| a.query_start.cmp(&b.query_start))
            .then_with(|| b.query_end.cmp(&a.query_end))
    });
    hits
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
    pub fn blosum62() -> Self {
        Self::from_type(MatrixType::Blosum62)
    }
    pub fn blosum45() -> Self {
        Self::from_type(MatrixType::Blosum45)
    }
    pub fn blosum50() -> Self {
        Self::from_type(MatrixType::Blosum50)
    }
    pub fn blosum80() -> Self {
        Self::from_type(MatrixType::Blosum80)
    }
    pub fn blosum90() -> Self {
        Self::from_type(MatrixType::Blosum90)
    }
    pub fn pam30() -> Self {
        Self::from_type(MatrixType::Pam30)
    }
    pub fn pam70() -> Self {
        Self::from_type(MatrixType::Pam70)
    }
    pub fn pam250() -> Self {
        Self::from_type(MatrixType::Pam250)
    }
    pub fn identity() -> Self {
        Self::from_type(MatrixType::Identity)
    }
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
    pub max_hsps: Option<usize>,
    pub culling_limit: Option<usize>,
    pub two_hit: bool,
    pub two_hit_window: usize,
    pub x_drop_final: i32,
    pub soft_masking: bool,
    pub lcase_masking: bool,
    /// Enable NCBI linked sum-statistics e-value adjustment where the program
    /// supports it. BLASTN enables this by default, matching CLI BLAST+.
    pub sum_stats: bool,
    /// Optional reusable Rayon pool for API searches that parallelize over
    /// database subjects. When present, it is used instead of constructing a
    /// per-search pool from `num_threads`.
    pub thread_pool: Option<Arc<rayon::ThreadPool>>,
}

impl SearchParams {
    pub fn blastp() -> Self {
        Self::blastp_defaults()
    }
    pub fn blastn() -> Self {
        Self::blastn_defaults()
    }

    pub fn blastx() -> Self {
        Self::blastp_defaults()
    }
    pub fn tblastn() -> Self {
        Self::blastp_defaults()
    }
    pub fn tblastx() -> Self {
        let mut params = Self::blastp_defaults();
        params.comp_adjust = 0;
        params
    }

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
            max_hsps: None,
            culling_limit: None,
            two_hit: false,
            two_hit_window: crate::stat::BLAST_WINDOW_SIZE_PROT as usize,
            x_drop_final: crate::stat::BLAST_GAP_X_DROPOFF_FINAL_PROT,
            soft_masking: false,
            lcase_masking: false,
            sum_stats: true,
            thread_pool: None,
        }
    }

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
            max_hsps: None,
            culling_limit: None,
            two_hit: false,
            // NCBI `BLAST_WINDOW_SIZE_NUCL = 0` would disable two-hit
            // regardless of `two_hit`; Rust keeps the protein default
            // so enabling two-hit on the blastn path still has a usable
            // window (40). `two_hit=false` by default, so unused.
            two_hit_window: crate::stat::BLAST_WINDOW_SIZE_PROT as usize,
            x_drop_final: crate::stat::BLAST_GAP_X_DROPOFF_FINAL_NUCL,
            soft_masking: false,
            lcase_masking: false,
            sum_stats: true,
            thread_pool: None,
        }
    }

    // Builder methods
    pub fn evalue(mut self, v: f64) -> Self {
        self.evalue_threshold = v;
        self
    }
    pub fn max_target_seqs(mut self, v: usize) -> Self {
        self.max_target_seqs = v;
        self
    }
    pub fn num_threads(mut self, v: usize) -> Self {
        self.num_threads = v;
        self
    }
    /// Use an existing Rayon thread pool for parallel API searches.
    ///
    /// Supplying a pool enables the parallel path even when `num_threads` is
    /// left at its default value; the pool's own size controls concurrency.
    pub fn thread_pool(mut self, pool: Arc<rayon::ThreadPool>) -> Self {
        self.thread_pool = Some(pool);
        self
    }
    pub fn filter_low_complexity(mut self, v: bool) -> Self {
        self.filter_low_complexity = v;
        self
    }
    pub fn seg_options(mut self, window: usize, locut: f64, hicut: f64) -> Self {
        self.seg_window = window;
        self.seg_locut = locut;
        self.seg_hicut = hicut;
        self
    }
    pub fn word_threshold(mut self, v: f64) -> Self {
        self.word_threshold = Some(v);
        self
    }
    pub fn comp_adjust(mut self, v: u8) -> Self {
        self.comp_adjust = v;
        self
    }
    pub fn strand(mut self, v: &str) -> Self {
        self.strand = v.to_string();
        self
    }
    pub fn word_size(mut self, v: usize) -> Self {
        self.word_size = v;
        self
    }
    pub fn matrix(mut self, v: MatrixType) -> Self {
        self.matrix = v;
        self
    }
    pub fn gap_open(mut self, v: i32) -> Self {
        self.gap_open = v;
        self
    }
    pub fn gap_extend(mut self, v: i32) -> Self {
        self.gap_extend = v;
        self
    }
    pub fn match_score(mut self, v: i32) -> Self {
        self.match_score = v;
        self
    }
    pub fn mismatch(mut self, v: i32) -> Self {
        self.mismatch = v;
        self
    }
    pub fn query_gencode(mut self, v: u8) -> Self {
        self.query_gencode = v;
        self
    }
    pub fn db_gencode(mut self, v: u8) -> Self {
        self.db_gencode = v;
        self
    }
    pub fn max_hsps(mut self, v: Option<usize>) -> Self {
        self.max_hsps = v;
        self
    }
    pub fn culling_limit(mut self, v: Option<usize>) -> Self {
        self.culling_limit = v;
        self
    }
    pub fn two_hit(mut self, v: bool) -> Self {
        self.two_hit = v;
        self
    }
    pub fn two_hit_window(mut self, v: usize) -> Self {
        self.two_hit_window = v;
        self
    }
    pub fn x_drop_final(mut self, v: i32) -> Self {
        self.x_drop_final = v;
        self
    }
    pub fn soft_masking(mut self, v: bool) -> Self {
        self.soft_masking = v;
        self
    }
    pub fn lcase_masking(mut self, v: bool) -> Self {
        self.lcase_masking = v;
        self
    }
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
    pub fn new(seq_type: DbType, db_title: impl Into<String>) -> Self {
        BlastDbBuilder {
            seq_type,
            db_title: db_title.into(),
            entries: Vec::new(),
        }
    }

    pub fn add(&mut self, entry: SequenceEntry) {
        self.entries.push(entry);
    }

    pub fn write(&self, base_path: &Path) -> io::Result<()> {
        match self.seq_type {
            DbType::Nucleotide => self.write_nucleotide(base_path),
            DbType::Protein => self.write_protein(base_path),
        }
    }

    fn write_nucleotide(&self, base_path: &Path) -> io::Result<()> {
        // Write .nsq
        let mut nsq = BufWriter::new(File::create(base_path.with_extension("nsq"))?);
        nsq.write_all(&[0u8])?; // sentinel

        let mut seq_offsets = vec![1u32];
        let mut amb_offsets = Vec::new();

        for entry in &self.entries {
            let seq = &entry.sequence;
            let iupac_to_2na = |b: u8| -> u8 {
                match b {
                    b'A' | b'a' => 0,
                    b'C' | b'c' => 1,
                    b'G' | b'g' => 2,
                    b'T' | b't' => 3,
                    _ => 0,
                }
            };

            let mut packed = Vec::new();
            let full_bytes = seq.len() / 4;
            let remainder = seq.len() % 4;
            for i in 0..full_bytes {
                let b = (iupac_to_2na(seq[i * 4]) << 6)
                    | (iupac_to_2na(seq[i * 4 + 1]) << 4)
                    | (iupac_to_2na(seq[i * 4 + 2]) << 2)
                    | iupac_to_2na(seq[i * 4 + 3]);
                packed.push(b);
            }
            if remainder > 0 {
                let mut last = 0u8;
                for j in 0..remainder {
                    last |= iupac_to_2na(seq[full_bytes * 4 + j]) << (6 - 2 * j);
                }
                last |= remainder as u8;
                packed.push(last);
            } else {
                packed.push(0);
            }

            nsq.write_all(&packed)?;
            let seq_start = *seq_offsets.last().unwrap();
            let amb_offset = seq_start + packed.len() as u32;
            amb_offsets.push(amb_offset);
            seq_offsets.push(amb_offset);
        }
        amb_offsets.push(*seq_offsets.last().unwrap_or(&0));
        nsq.flush()?;

        // Write .nhr
        let mut nhr = BufWriter::new(File::create(base_path.with_extension("nhr"))?);
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
            &base_path.with_extension("nin"),
            4, // format version
            DbType::Nucleotide,
            &self.db_title,
            self.entries.len() as u32,
            &hdr_offsets,
            &seq_offsets,
            Some(&amb_offsets),
        )
    }

    fn write_protein(&self, base_path: &Path) -> io::Result<()> {
        // Write .psq (protein sequences in NCBIstdaa)
        let mut psq = BufWriter::new(File::create(base_path.with_extension("psq"))?);
        psq.write_all(&[0u8])?; // sentinel

        let mut seq_offsets = vec![1u32];
        for entry in &self.entries {
            let encoded: Vec<u8> = entry
                .sequence
                .iter()
                .map(|&b| AMINOACID_TO_NCBISTDAA[b as usize & 0x7F])
                .collect();
            psq.write_all(&encoded)?;
            psq.write_all(&[0u8])?; // sentinel between sequences
            let prev = *seq_offsets.last().unwrap();
            seq_offsets.push(prev + encoded.len() as u32 + 1);
        }
        psq.flush()?;

        // Write .phr
        let mut phr = BufWriter::new(File::create(base_path.with_extension("phr"))?);
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
            &base_path.with_extension("pin"),
            4,
            DbType::Protein,
            &self.db_title,
            self.entries.len() as u32,
            &hdr_offsets,
            &seq_offsets,
            None,
        )
    }
}

fn encode_defline_asn1(header: &str, _oid: i32) -> Vec<u8> {
    // Minimal ASN.1 BER encoding of Blast-def-line-set
    let title_bytes = header.as_bytes();
    let inner_len = 2 + title_bytes.len();
    let mut out = Vec::with_capacity(inner_len + 10);
    // SEQUENCE { VisibleString title }
    out.push(0x30); // SEQUENCE tag
    if inner_len < 128 {
        out.push(inner_len as u8);
    } else {
        out.push(0x81);
        out.push(inner_len as u8);
    }
    out.push(0x1A); // VisibleString tag
    out.push(title_bytes.len() as u8);
    out.extend_from_slice(title_bytes);
    out
}

fn write_index_file(
    path: &Path,
    format_version: u32,
    db_type: DbType,
    title: &str,
    num_oids: u32,
    hdr_offsets: &[u32],
    seq_offsets: &[u32],
    amb_offsets: Option<&[u32]>,
) -> io::Result<()> {
    use byteorder::{BigEndian, WriteBytesExt};
    let mut f = BufWriter::new(File::create(path)?);

    f.write_u32::<BigEndian>(format_version)?;
    let db_type_val: u32 = match db_type {
        DbType::Protein => 1,
        DbType::Nucleotide => 0,
    };
    f.write_u32::<BigEndian>(db_type_val)?;

    // Title (length-prefixed)
    f.write_u32::<BigEndian>(title.len() as u32)?;
    f.write_all(title.as_bytes())?;

    // Timestamp placeholder
    let ts = "2024-01-01T00:00:00";
    f.write_u32::<BigEndian>(ts.len() as u32)?;
    f.write_all(ts.as_bytes())?;

    f.write_u32::<BigEndian>(num_oids)?;

    // Total residue count (not critical for search)
    let total: u64 = 0;
    f.write_u64::<BigEndian>(total)?;
    // Max seq length
    f.write_u32::<BigEndian>(0)?;

    // Write offset arrays
    for &off in hdr_offsets {
        f.write_u32::<BigEndian>(off)?;
    }
    for &off in seq_offsets {
        f.write_u32::<BigEndian>(off)?;
    }
    if let Some(amb) = amb_offsets {
        for &off in amb {
            f.write_u32::<BigEndian>(off)?;
        }
    }

    f.flush()
}

// ── Search functions ────────────────────────────────────────────────────────

/// Run a blastp search (protein query vs protein database).
pub fn blastp(db: &BlastDb, query: &[u8], params: &SearchParams) -> Vec<SearchResult> {
    if query.is_empty() {
        return Vec::new();
    }

    let query_aa = encode_protein_query_nomask(query);
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
    if crate::composition::read_composition(&query_aa_masked, AA_SIZE).1 == 0 {
        return Vec::new();
    }

    let matrix = *get_matrix(params.matrix);
    let word_size = params.word_size.clamp(2, 6);
    let threshold = params
        .word_threshold
        .unwrap_or_else(|| suggested_word_threshold(params.matrix, crate::program::BLASTP));

    let prot_kbp = protein_kbp_for_matrix(params.matrix, params.gap_open, params.gap_extend);

    // Convert bit-score x_dropoff values to raw scores.
    // NCBI uses UNGAPPED KBP for ungapped x_drop, GAPPED KBP for gapped x_drop.
    // Ungapped path: `ceil()` then Int4 cast (`blast_parameters.c:221`).
    // Gapped path: plain `(Int4)` cast — truncation toward zero
    // (`blast_parameters.c:457-463`).
    let ln2 = crate::math::NCBIMATH_LN2;
    // NCBI's `Blast_ScoreBlkKbpUngappedCalc` (`blast_stat.c:2737`) populates
    // `sbp->kbp_std[context]` from the **query's actual amino-acid composition**,
    // not the ideal Robinson background. For most queries the result is close
    // to the ideal, but the small drift in lambda/logK shifts `gap_trigger` by
    // 1 raw score unit — enough to swing boundary hits like seqp's DAA02208
    // (max ungapped score 40 vs ideal cutoff 41 vs query-specific cutoff 40).
    let ungapped_kbp = crate::stat::query_specific_protein_ungapped_kbp(&query_aa, &matrix);
    let x_drop_ungapped = (params.x_drop_ungapped as f64 * ln2 / ungapped_kbp.lambda).ceil() as i32;
    // NCBI's PRELIMINARY gapped extension uses `gap_x_dropoff` (15 bits for
    // protein default) — `blast_parameters.c:457-463` truncates with `(Int4)`.
    // The TRACEBACK phase then re-extends with the larger `gap_x_dropoff_final`
    // (25 bits). Using `x_drop_final` here makes our preliminary gapped DP
    // far more permissive than NCBI's, finding gap-rich alignments NCBI rejects
    // (e.g. NP_777001 score=50 vs NCBI's 46).
    let x_drop_gapped = (params.x_drop_gapped as f64 * ln2 / prot_kbp.lambda) as i32;
    let x_drop_final = (params.x_drop_final as f64 * ln2 / prot_kbp.lambda) as i32;

    // Gap trigger: NCBI `blast_parameters.c:343`:
    //   (Int4)((gap_trigger * NCBIMATH_LN2 + kbp->logK) / kbp->Lambda)
    // where `kbp` is the UNGAPPED KBP (per-context, derived from query composition).
    let gap_trigger_raw = ((crate::stat::BLAST_GAP_TRIGGER_PROT * ln2 + ungapped_kbp.log_k)
        / ungapped_kbp.lambda) as i32;

    let total_subj_len: usize = (0..db.num_oids)
        .map(|oid| db.get_seq_len(oid) as usize)
        .sum();
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

    // Use exact length adjustment with alpha/beta from gapped params (matching NCBI C engine)
    let gapped_params =
        protein_gapped_params_for_matrix(params.matrix, params.gap_open, params.gap_extend);
    let (_len_adj, search_space) = if let Some(ref gp) = gapped_params {
        let alpha_d_lambda = gp.alpha / prot_kbp.lambda;
        let (adj, _) = crate::stat::compute_length_adjustment_exact(
            prot_kbp.k,
            prot_kbp.log_k,
            alpha_d_lambda,
            gp.beta,
            query_aa.len() as i32,
            total_subj_len as i64,
            db.num_oids as i32,
        );
        let eff_q = (query_aa.len() as i64 - adj as i64).max(1);
        let eff_db = (total_subj_len as i64 - db.num_oids as i64 * adj as i64).max(1);
        (adj, eff_q as f64 * eff_db as f64)
    } else {
        let adj = crate::stat::compute_length_adjustment(
            query_aa.len() as i32,
            total_subj_len as i64,
            db.num_oids as i32,
            &prot_kbp,
        );
        let ss = crate::stat::compute_search_space(
            query_aa.len() as i64,
            total_subj_len as i64,
            db.num_oids as i32,
            adj,
        );
        (adj, ss)
    };

    // Build Gumbel block for Spouge FSC e-value (per-subject length correction)
    let gumbel_blk = protein_gumbel_for_matrix(
        params.matrix,
        params.gap_open,
        params.gap_extend,
        total_subj_len as i64,
    );

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

    // Process a single subject OID — extracted so it can be called
    // either sequentially or in parallel without per-call pool overhead.
    let search_oid = |oid: u32| -> Option<SearchResult> {
        let subj_raw = db.get_sequence(oid);
        let subj_len = db.get_seq_len(oid) as usize;
        if subj_len < word_size {
            return None;
        }
        let trace_acc = std::env::var("NB_TRACE_ACC").ok();
        let subject_accession = db
            .get_accession(oid)
            .unwrap_or_else(|| format!("oid_{}", oid));
        let do_trace = trace_acc
            .as_deref()
            .map_or(false, |want| want == subject_accession);

        // Use length-based slice — no allocation (matches NCBI C approach).
        // get_sequence() includes trailing sentinel; subj_len excludes it.
        let subj_aa = &subj_raw[..subj_len];

        let ungapped_hits = crate::protein_lookup::protein_scan_with_table(
            &query_aa,
            subj_aa,
            &matrix,
            &lookup_table,
            x_drop_ungapped,
        );
        if ungapped_hits.is_empty() {
            return None;
        }

        // NCBI's per-context seed cutoff is `MIN(gap_trigger,
        // cutoff_score_max)` (`blast_parameters.c:367`). When gbp is filled,
        // `cutoff_score_max` is `BLAST_SpougeEtoS(..., query_length,
        // min_subject_length)` (`blast_parameters.c:935-941`), with the DB
        // minimum subject length from `BlastSeqSrcGetMinSeqLen`
        // (`blast_setup.c:970`), not the average.
        // NCBI also uses `cbs_stretch * evalue` (= 5*evalue when comp_adjust>1)
        // to make the prelim cutoff more permissive, on the theory that
        // composition adjustment will rescue the borderline seeds. Until our
        // `composition_matrix_adj` matches NCBI's exactly (iter-25 / iter-49
        // known divergence on certain compositions), enabling cbs_stretch
        // pushes through borderline subjects whose comp-adjusted scores
        // diverge from NCBI's, producing both FPs (when our adj DP scores
        // higher) and false-bound mismatches (when our adj DP scores lower).
        // Keep the strict cutoff for now.
        let _cbs_stretch: f64 = if params.comp_adjust > 1 { 5.0 } else { 1.0 };
        let eval_cutoff = protein_eval_cutoff(
            evalue_threshold,
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
        // Chaining: NCBI `s_ChainingAlignment` (`blast_gapalign.c:3592`) runs
        // before gapped DP for blastp+BLOSUM62 by default
        // (`ext_params->options->chaining`). It computes a chained-score
        // approximation for each ungapped seed and DROPS seeds whose
        // best chained score (minus a single gap penalty plus
        // word_cutoff - 1) cannot reach the hit_cutoff. Without this,
        // we send many more seeds to gapped DP than NCBI does — each
        // can produce slightly different bounds and scores via
        // X-drop-band differences in repeats, breaking parity even when
        // every other function is byte-faithful.
        //
        // ungapped_hits is sorted by score desc, but chaining wants
        // sort by query offset (asc) within each context. We have a
        // single context per call, so just sort by q_start.
        let mut chained_hits: Vec<crate::protein_lookup::ProteinHit> =
            ungapped_hits.iter().cloned().collect();
        chained_hits.sort_by_key(|h| h.query_start);
        let gap_score = gap_open + gap_extend;
        let n = chained_hits.len();
        // best_score[k] starts at chained_hits[k].score and gets updated
        // by chaining DP. Process k from last to first.
        let mut best_score: Vec<i32> = chained_hits.iter().map(|h| h.score).collect();
        for k in (0..n).rev() {
            let self_score = best_score[k];
            for j in (k + 1)..n {
                let q_diff = chained_hits[j].query_start as i32
                    - chained_hits[k].query_start as i32
                    + (chained_hits[k].query_end - chained_hits[k].query_start) as i32;
                let s_diff = chained_hits[j].subject_start as i32
                    - chained_hits[k].subject_start as i32
                    + (chained_hits[k].subject_end - chained_hits[k].subject_start) as i32;
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
        // Drop chained_hits[k] when its chained score can't reach hit_cutoff:
        //   best_score[k] - gap_score + word_cutoff - 1 < hit_cutoff
        // i.e. keep when best_score[k] >= hit_cutoff + gap_score - word_cutoff + 1.
        let mut keep = vec![true; n];
        for k in 0..n {
            if best_score[k] - gap_score + word_cutoff - 1 < hit_cutoff {
                keep[k] = false;
            }
        }
        let chained_kept: std::collections::HashSet<(usize, usize)> = chained_hits
            .iter()
            .zip(keep.iter())
            .filter_map(|(h, &k)| {
                if k {
                    Some((h.query_start, h.subject_start))
                } else {
                    None
                }
            })
            .collect();
        // Mirror NCBI's `BLAST_GetGappedScore` flow (`blast_gapalign.c:3925+`):
        //   1. Sort ungapped HSPs by score descending.
        //   2. For each ungapped HSP, check `BlastIntervalTreeContainsHSP`
        //      against already-accepted gapped HSPs (line 3993). If the
        //      UNGAPPED bounds are contained in any accepted GAPPED HSP and
        //      ungapped.score <= gapped.score, SKIP (no gapped DP run).
        //   3. Otherwise run preliminary gapped DP (`x_drop_gapped`) →
        //      cutoff filter → traceback re-extension (`x_drop_final`).
        //   4. After accepting, ALSO check the resulting gapped HSP isn't
        //      contained in another (`s_HSPIsContained`, `blast_itree.c:810`).
        //
        // The PRE-gapped containment check is what suppresses the iter-101
        // NP_777001 FP: seed[1]'s ungapped (q=59-68 score=43) is contained
        // in seed[0]'s gapped (q=57-75 score=46) on both axes with
        // 43 <= 46. NCBI skips. Our previous impl ran gapped on seed[1],
        // got score=50, and either kept it (iter-105 → FP) or rejected via
        // a stricter overlap rule (iter-104 → also rejects legit overlapping
        // HSPs). Mirroring NCBI's exact pre-check resolves both.
        let mut passing_seeds: Vec<&crate::protein_lookup::ProteinHit> = ungapped_hits
            .iter()
            .filter(|uh| {
                uh.score >= adjusted_cutoff
                    && chained_kept.contains(&(uh.query_start, uh.subject_start))
            })
            .collect();
        passing_seeds.sort_by(|a, b| b.score.cmp(&a.score));
        if do_trace {
            eprintln!(
                "[trace acc={}] ungapped={} passing={} cutoff={}",
                subject_accession,
                ungapped_hits.len(),
                passing_seeds.len(),
                adjusted_cutoff
            );
            for uh in &ungapped_hits {
                eprintln!(
                    "  ungapped q={}-{} s={}-{} score={}",
                    uh.query_start, uh.query_end, uh.subject_start, uh.subject_end, uh.score
                );
            }
        }
        let mut phits: Vec<crate::protein_lookup::ProteinHit> = Vec::new();
        let mut tree =
            crate::itree::IntervalTree::new(query_aa.len() as i32 + 1, subj_len as i32 + 1);
        for uh in &passing_seeds {
            // Ungapped span (NCBI's `tmp_hsp` constructed at line 3976-3990
            // from `init_hsp->ungapped_data`).
            let ungap_q_start = uh.query_start;
            let ungap_q_end = uh.query_end;
            let ungap_s_start = uh.subject_start;
            let ungap_s_end = uh.subject_end;
            let ungap_score = uh.score;
            // PRE-gapped containment check: skip if this ungapped HSP is
            // fully enveloped by an already-accepted gapped HSP with score
            // greater-or-equal (`s_HSPIsContained`, `blast_itree.c:810`).
            let pre_contained = tree.is_contained(
                crate::itree::Interval::new(ungap_q_start as i32, ungap_q_end as i32),
                crate::itree::Interval::new(ungap_s_start as i32, ungap_s_end as i32),
                ungap_score,
            );
            if pre_contained {
                continue;
            }
            let (seed_q, seed_s) = crate::protein::get_start_for_gapped_alignment(
                &query_aa,
                subj_aa,
                uh.query_start,
                uh.query_end.saturating_sub(uh.query_start),
                uh.subject_start,
                uh.subject_end.saturating_sub(uh.subject_start),
                &matrix,
            );
            // PRELIMINARY gapped DP only (matches NCBI engine flow:
            // `s_BlastProtGappedAlignment` calls `Blast_SemiGappedAlign`
            // with `gap_x_dropoff` only — the larger `gap_x_dropoff_final`
            // is used in `Blast_TracebackFromHSPList`, NOT in the engine).
            // Using x_drop_final here makes our containment tree contain
            // larger bounds than NCBI's, which envelops legitimate seeds for
            // weaker HSPs. The traceback (final-xdrop) re-runs after the
            // engine completes — see the post-loop block.
            let Some(prelim) = crate::protein::protein_gapped_align(
                &query_aa,
                subj_aa,
                seed_q,
                seed_s,
                &matrix,
                gap_open,
                gap_extend,
                x_drop_gapped,
            ) else {
                continue;
            };
            if do_trace {
                eprintln!(
                    "  prelim from ungapped q={}-{} s={}-{} score={} -> seed=({}, {}) prelim_score={} bounds q={}-{} s={}-{}",
                    uh.query_start, uh.query_end, uh.subject_start, uh.subject_end, uh.score,
                    seed_q, seed_s, prelim.score,
                    prelim.query_start, prelim.query_end, prelim.subject_start, prelim.subject_end
                );
            }
            if prelim.score < adjusted_cutoff {
                continue;
            }
            // Post-gapped containment check uses PRELIMINARY bounds
            // (`BlastIntervalTreeContainsHSP` in NCBI's engine sees
            // `Blast_SemiGappedAlign` bounds, not traceback bounds).
            let post_contained = tree.is_contained(
                crate::itree::Interval::new(prelim.query_start as i32, prelim.query_end as i32),
                crate::itree::Interval::new(prelim.subject_start as i32, prelim.subject_end as i32),
                prelim.score,
            );
            if post_contained {
                continue;
            }
            if do_trace {
                eprintln!(
                    "  prelim accepted q={}-{} s={}-{} score={}",
                    prelim.query_start,
                    prelim.query_end,
                    prelim.subject_start,
                    prelim.subject_end,
                    prelim.score
                );
            }
            // Stash the (seed_q, seed_s) so we can re-run final-xdrop DP after
            // the engine loop closes (NCBI's traceback phase).
            phits.push(crate::protein_lookup::ProteinHit {
                query_start: prelim.query_start,
                query_end: prelim.query_end,
                subject_start: prelim.subject_start,
                subject_end: prelim.subject_end,
                score: prelim.score,
                num_ident: prelim.num_ident,
                align_length: prelim.align_length,
                mismatches: prelim.mismatches,
                gap_opens: prelim.gap_opens,
                // Pack seed for later traceback re-run via qseq/sseq.
                // We don't have a real edit script yet — mark so the
                // post-loop traceback can detect and refresh.
                qseq: None,
                sseq: None,
                scaled_score: Some(((seed_q as i32) << 16) | (seed_s as i32 & 0xffff)),
                gapped_start_q: seed_q,
                gapped_start_s: seed_s,
            });
            tree.insert(
                crate::itree::Interval::new(prelim.query_start as i32, prelim.query_end as i32),
                crate::itree::Interval::new(prelim.subject_start as i32, prelim.subject_end as i32),
                prelim.score,
            );
        }
        // TRACEBACK phase: re-run gapped DP with `x_drop_final` for each
        // accepted preliminary HSP. Mirrors NCBI's `Blast_TracebackFromHSPList`
        // loop (`blast_traceback.c:375-625`):
        // 1. Sort HSPs by prelim score desc.
        // 2. For each HSP in order, check `BlastIntervalTreeContainsHSP`
        //    against tree of already-tracebacked HSPs (line 404).
        // 3. If contained: drop HSP without running ALIGN_EX.
        // 4. Else: run ALIGN_EX (final xdrop), update bounds, add to tree.
        //
        // Without the containment check, redundant traceback DPs from
        // different seeds can produce slightly different bounds and
        // scores in repetitive regions, leaving multiple HSPs that the
        // common-endpoints purge keeps the higher-scoring one of —
        // but NCBI never computed the higher-scoring one because its
        // traceback skipped the seed.
        phits.sort_by(|a, b| b.score.cmp(&a.score));
        let mut tb_tree =
            crate::itree::IntervalTree::new(query_aa.len() as i32 + 1, subj_len as i32 + 1);
        let mut keep = vec![true; phits.len()];
        for (idx, ph) in phits.iter_mut().enumerate() {
            let Some(packed) = ph.scaled_score.take() else {
                continue;
            };
            let pre_contained = tb_tree.is_contained(
                crate::itree::Interval::new(ph.query_start as i32, ph.query_end as i32),
                crate::itree::Interval::new(ph.subject_start as i32, ph.subject_end as i32),
                ph.score,
            );
            if pre_contained {
                keep[idx] = false;
                continue;
            }
            let seed_q = (packed >> 16) as usize;
            let seed_s = (packed & 0xffff) as usize;
            let Some(gr) = crate::protein::protein_gapped_align(
                &query_aa,
                subj_aa,
                seed_q,
                seed_s,
                &matrix,
                gap_open,
                gap_extend,
                x_drop_final,
            ) else {
                keep[idx] = false;
                continue;
            };
            let q_slice = &query_aa[gr.query_start..gr.query_end];
            let s_slice = &subj_aa[gr.subject_start..gr.subject_end];
            let (qs, ss) = gr.edit_script.render_alignment(
                q_slice,
                s_slice,
                crate::protein::ncbistdaa_to_char,
            );
            ph.query_start = gr.query_start;
            ph.query_end = gr.query_end;
            ph.subject_start = gr.subject_start;
            ph.subject_end = gr.subject_end;
            ph.score = gr.score;
            ph.num_ident = gr.num_ident;
            ph.align_length = gr.align_length;
            ph.mismatches = gr.mismatches;
            ph.gap_opens = gr.gap_opens;
            ph.qseq = Some(qs);
            ph.sseq = Some(ss);
            tb_tree.insert(
                crate::itree::Interval::new(ph.query_start as i32, ph.query_end as i32),
                crate::itree::Interval::new(ph.subject_start as i32, ph.subject_end as i32),
                ph.score,
            );
        }
        let mut idx = 0usize;
        phits.retain(|_| {
            let k = keep[idx];
            idx += 1;
            k
        });
        // Mirror NCBI's `Blast_HSPListPurgeHSPsWithCommonEndpoints`
        // (`blast_hits.c:2455`), which the engine calls at
        // `blast_engine.c:544` after `BLAST_GetGappedScore`. Removes HSPs
        // that share either (query.offset, subject.offset) OR
        // (query.end, subject.end). Within each group of duplicates, the
        // highest-scoring HSP wins (sort tiebreaker by score descending).
        // For protein this is the FULL purge variant (FREE duplicates,
        // not preserve via gap-edit cutoff).
        // NCBI forces `purge=TRUE` for protein at `blast_hits.c:2464`
        // (`purge |= (program != eBlastTypeBlastn)`). Pass `true` here so
        // common-endpoint duplicates are FREED, not preserved via
        // cut-off-edit-script (the cut-off path is the blastn-only branch).
        purge_hsps_with_common_endpoints(&mut phits);
        phits.sort_by(|a, b| {
            b.score
                .cmp(&a.score)
                .then_with(|| a.subject_start.cmp(&b.subject_start))
                .then_with(|| b.subject_end.cmp(&a.subject_end))
                .then_with(|| a.query_start.cmp(&b.query_start))
                .then_with(|| b.query_end.cmp(&a.query_end))
        });
        // Mirror NCBI's traceback-time containment pass
        // (`blast_traceback.c:371-405`). NCBI builds a fresh interval tree
        // and processes HSPs in **gapped score desc** order, dropping any
        // HSP enveloped by an already-accepted higher-scoring one.
        // The engine-loop containment (iter 120 / 122) processes seeds in
        // **ungapped score** order, so a low-ungapped-score seed that
        // produces a high-gapped-score HSP can arrive AFTER several
        // high-ungapped-score seeds whose gapped HSPs are smaller — those
        // smaller HSPs end up accepted, then the bigger one too.
        // Re-running containment on the gapped list in score order kills
        // those duplicates. Concrete case (iter 123): AAC46500 vs
        // XP_353871, where 4 seeds with ungapped 64/61 produced gapped
        // HSPs at scores 94/87/76/70 all contained in HSP_B (gapped=98
        // from a seed with ungapped=56).
        {
            let mut accepted: Vec<crate::protein_lookup::ProteinHit> =
                Vec::with_capacity(phits.len());
            let mut tree =
                crate::itree::IntervalTree::new(query_aa.len() as i32 + 1, subj_len as i32 + 1);
            for ph in phits.drain(..) {
                let contained = tree.is_contained(
                    crate::itree::Interval::new(ph.query_start as i32, ph.query_end as i32),
                    crate::itree::Interval::new(ph.subject_start as i32, ph.subject_end as i32),
                    ph.score,
                );
                if !contained {
                    tree.insert(
                        crate::itree::Interval::new(ph.query_start as i32, ph.query_end as i32),
                        crate::itree::Interval::new(ph.subject_start as i32, ph.subject_end as i32),
                        ph.score,
                    );
                    accepted.push(ph);
                }
            }
            phits = accepted;
        }
        if phits.is_empty() {
            return None;
        }
        for ph in &mut phits {
            ph.score = rescore_protein_hit(ph, &query_aa, subj_aa, &matrix, gap_open, gap_extend);
        }

        // Pre-filter: check e-value with Spouge FSC if available, else simple
        // Karlin. NCBI uses `hit_params->prelim_evalue` here
        // (`blast_engine.c:653`) which is `cbs_stretch * evalue`. Faithfully
        // matching NCBI requires our composition-adjusted DP to match NCBI's,
        // since the boosted score can either keep or drop the HSP at the final
        // evalue=10 gate. Our `composition_matrix_adj` (Newton optimization)
        // diverges from NCBI's on certain query/subject compositions
        // (iter-25 / iter-49 known issue), producing higher scores for some
        // borderline subjects (e.g. NP_982592→NP_777001 score=1358 ours vs
        // NCBI's 403 for similar bounds). So using `prelim_evalue` here lets
        // those FPs through. Until comp_adjust matrix matches NCBI's, gate
        // with `evalue_threshold` (= the strict user e-value). This keeps the
        // seed fix + cbs_stretch wins (more legitimate seeds) while preventing
        // FPs from the matrix divergence.
        let best_raw_ev = phits
            .iter()
            .map(|ph| {
                if let Some(ref gbp) = gumbel_blk {
                    crate::stat::spouge_evalue(
                        ph.score,
                        &prot_kbp,
                        gbp,
                        query_aa.len() as i32,
                        subj_len as i32,
                    )
                } else {
                    prot_kbp.raw_to_evalue(ph.score, search_space)
                }
            })
            .fold(f64::MAX, f64::min);
        if best_raw_ev > evalue_threshold {
            return None;
        }

        let accession = subject_accession;
        let title = String::from_utf8_lossy(db.get_header(oid)).to_string();
        let sl = subj_aa.len();

        // Composition-based e-value adjustment (NCBI comp_based_stats).
        // Verbatim port of Blast_AdjustScores + s_AdjustEvaluesForComposition.
        let comp_mode = params.comp_adjust;
        let redo_subj_storage = if comp_mode > 0 {
            kappa_redo_subject_sequence(query_aa.len(), subj_aa, &phits, prot_kbp.lambda)
        } else {
            std::borrow::Cow::Borrowed(subj_aa)
        };
        let redo_subj_aa: &[u8] = redo_subj_storage.as_ref();
        let comp_scale = if comp_mode > 0 {
            COMPO_ADJUST_SCALE_FACTOR
        } else {
            1.0
        };
        let scaled_gap_open = crate::math::nint(gap_open as f64 * comp_scale) as i32;
        let scaled_gap_extend = crate::math::nint(gap_extend as f64 * comp_scale) as i32;
        let scaled_x_drop_final = crate::math::nint(x_drop_final as f64 * comp_scale) as i32;

        // Determine adjustment rule and optionally build adjusted matrix.
        // adj_result: None = no adjustment, Some((adjusted_matrix_opt, lambda_ratio_opt))
        type AdjResult = Option<(Option<[[i32; AA_SIZE]; AA_SIZE]>, Option<f64>)>;
        let adj_result: AdjResult = if comp_mode > 0 {
            let (qcomp28, qn) = crate::composition::read_composition(&query_aa, AA_SIZE);
            let (scomp28, sn) = crate::composition::read_composition(redo_subj_aa, AA_SIZE);
            if qn == 0 || sn == 0 {
                Some((None, None))
            } else {
                let mut qp20 = [0.0f64; 20];
                let mut sp20 = [0.0f64; 20];
                crate::compo_mode_condition::gather_letter_probs(&qcomp28, &mut qp20);
                crate::compo_mode_condition::gather_letter_probs(&scomp28, &mut sp20);

                let rule = crate::compo_mode_condition::choose_matrix_adjust_rule(
                    query_aa.len(),
                    subj_aa.len(),
                    &qp20,
                    &sp20,
                    comp_mode,
                );

                use crate::compo_mode_condition::MatrixAdjustRule;
                match rule {
                    MatrixAdjustRule::DontAdjust => None,
                    MatrixAdjustRule::ScaleOldMatrix => {
                        // Port of NCBI Blast_CompositionBasedStats: rescale matrix
                        // using composition-specific lambda ratio, then re-align.
                        let ungapped_lambda =
                            crate::stat::protein_ungapped_kbp().lambda / comp_scale;
                        // Build the frequency ratio matrix from the standard BLOSUM62
                        // joint probs (used as freq_ratios in s_ScaleSquareMatrix).
                        // For non-position-based, startFreqRatios is initialized from
                        // the standard matrix frequency ratios.
                        let freq_ratios = crate::matrix::get_blosum62_freq_ratios();
                        let start_matrix = crate::composition::matrix_from_freq_ratios(
                            ungapped_lambda,
                            &freq_ratios,
                        );
                        crate::composition::composition_scale_matrix(
                            &start_matrix,
                            &qcomp28,
                            &scomp28,
                            ungapped_lambda,
                            &freq_ratios,
                        )
                        .map(|adj_mat| (Some(adj_mat), None))
                    }
                    MatrixAdjustRule::UserSpecifiedRelEntropy
                    | MatrixAdjustRule::UnconstrainedRelEntropy
                    | MatrixAdjustRule::RelEntropyOldMatrixNewContext
                    | MatrixAdjustRule::RelEntropyOldMatrixOldContext => {
                        // Full matrix optimization (Blast_CompositionMatrixAdj)
                        let (joint_probs, first_std, second_std) =
                            crate::composition::blosum62_workspace();
                        let mut adj_matrix = matrix;
                        // NCBI uses ungappedLambda (0.3176 for BLOSUM62) for matrix scaling,
                        // NOT the gapped lambda. See matrixInfo->ungappedLambda.
                        let ungapped_lambda =
                            crate::stat::protein_ungapped_kbp().lambda / comp_scale;
                        let freq_ratios = crate::matrix::get_blosum62_freq_ratios();
                        let start_matrix = crate::composition::matrix_from_freq_ratios(
                            ungapped_lambda,
                            &freq_ratios,
                        );
                        let status = crate::composition::composition_matrix_adj(
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
                            crate::composition::composition_scale_matrix(
                                &start_matrix,
                                &qcomp28,
                                &scomp28,
                                ungapped_lambda,
                                &freq_ratios,
                            )
                            .map(|adj_mat| (Some(adj_mat), None))
                        }
                    }
                }
            }
        } else {
            None
        };

        // If we have an adjusted matrix, re-align and recompute scores.
        // 1-1 with NCBI's `Blast_RedoAlignmentCore_MT` (single-HSP arm,
        // redo_alignment.c:1430-1530): full SW under the adjusted
        // matrix → reverse SW for the start → forward-only X-drop
        // bounded by the SW endpoints. The bounded forward-only X-drop
        // (vs bidirectional X-drop from the seed) prevents
        // alignment-end overshoot that produced longer-but-suboptimal
        // alignments and ~2-bit discrepancies vs NCBI on the
        // sortase / dinB fixtures.
        let sw_bounded_blastp = if let Some((Some(ref adj_mat), _)) = adj_result {
            if phits.len() == 1 {
                let (sw_score, m_end_sw, q_end_sw) =
                    crate::smith_waterman::blast_smith_waterman_score_only(
                        redo_subj_aa,
                        &query_aa,
                        adj_mat,
                        scaled_gap_open,
                        scaled_gap_extend,
                    );
                if sw_score > 0 {
                    let (_, m_start_sw, q_start_sw) =
                        crate::smith_waterman::blast_smith_waterman_find_start(
                            redo_subj_aa,
                            &query_aa,
                            adj_mat,
                            scaled_gap_open,
                            scaled_gap_extend,
                            m_end_sw,
                            q_end_sw,
                            sw_score,
                        );
                    let q_extent = q_end_sw - q_start_sw;
                    let s_extent = m_end_sw - m_start_sw;
                    Some((q_start_sw, m_start_sw, q_extent, s_extent, sw_score))
                } else {
                    None
                }
            } else {
                None
            }
        } else {
            None
        };
        let (final_phits, use_adj_matrix) = if let Some((None, _)) = adj_result {
            (Vec::new(), true)
        } else if let Some((Some(ref adj_mat), _)) = adj_result {
            // Re-do gapped alignment with adjusted matrix
            let mut new_phits = Vec::new();
            for ph in &phits {
                // Run BOTH bounded SW (NCBI's `s_SWFindFinalEndsUsingXdrop`-style)
                // and bidirectional X-drop. The bounded helper avoids
                // alignment-end overshoot but can pick a low-identity SW path
                // when multiple equal-score paths exist (iter 50 finding).
                // Picking the strictly-higher-scoring result keeps the
                // overshoot fix while preserving NCBI-matching identity
                // counts on cases where bidirectional already finds the
                // correct optimum.
                let bounded_gr = if let Some((q_start, m_start, q_extent, s_extent, target_score)) =
                    sw_bounded_blastp
                {
                    crate::protein::protein_sw_bounded_xdrop_align(
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
                let bidir_gr = crate::protein::protein_gapped_align(
                    &query_aa,
                    redo_subj_aa,
                    ph.gapped_start_q,
                    ph.gapped_start_s,
                    adj_mat,
                    scaled_gap_open,
                    scaled_gap_extend,
                    scaled_x_drop_final,
                );
                let gr_opt = match (bounded_gr, bidir_gr) {
                    (Some(b), Some(d)) => {
                        // Prefer bounded ONLY when it strictly outscores
                        // bidirectional — bidirectional's path tends to
                        // match NCBI's identity counts when scores tie.
                        if b.score > d.score {
                            Some(b)
                        } else {
                            Some(d)
                        }
                    }
                    (Some(b), None) => Some(b),
                    (None, d) => d,
                };
                if let Some(gr) = gr_opt {
                    let q_slice = &query_aa[gr.query_start..gr.query_end];
                    let s_slice = &redo_subj_aa[gr.subject_start..gr.subject_end];
                    let (qs, ss) = gr.edit_script.render_alignment(
                        q_slice,
                        s_slice,
                        crate::protein::ncbistdaa_to_char,
                    );
                    new_phits.push(crate::protein_lookup::ProteinHit {
                        query_start: gr.query_start,
                        query_end: gr.query_end,
                        subject_start: gr.subject_start,
                        subject_end: gr.subject_end,
                        score: crate::math::nint(gr.score as f64 / comp_scale) as i32,
                        num_ident: gr.num_ident,
                        align_length: gr.align_length,
                        mismatches: gr.mismatches,
                        gap_opens: gr.gap_opens,
                        qseq: Some(qs),
                        sseq: Some(ss),
                        scaled_score: Some(gr.score),
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

        let hsps: Vec<Hsp> = final_phits
            .iter()
            .filter_map(|ph| {
                // Mirror NCBI's `Blast_HSPListGetEvalues` (`blast_hits.c:1873`):
                // when comp_adjust is on, NCBI computes the e-value with the
                // *scaled* score and Lambda/scale_factor, then later rounds the
                // score back via `s_HSPListNormalizeScores`. Rounding before the
                // e-value loses up to 0.5 raw-score units of precision and
                // introduces a ~0.92× drift in `comp_adjust >= 1`. When
                // `scaled_score` is set, use the scaled values directly.
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
                // Use Spouge FSC for e-value when Gumbel params are available.
                // This gives per-subject-length corrected e-values matching NCBI.
                let evalue = if let Some(ref gbp) = gumbel_blk {
                    let base_ev = crate::stat::spouge_evalue(
                        e_score_i32,
                        &e_kbp,
                        gbp,
                        query_aa.len() as i32,
                        sl as i32,
                    );
                    if use_adj_matrix {
                        base_ev
                    } else if let Some(lr) = lambda_ratio_opt {
                        // Scale the e-value by the lambda ratio
                        let scaled_kbp = crate::stat::KarlinBlk {
                            lambda: e_kbp.lambda / lr,
                            k: e_kbp.k,
                            log_k: e_kbp.log_k,
                            h: e_kbp.h,
                            round_down: e_kbp.round_down,
                        };
                        crate::stat::spouge_evalue(
                            e_score_i32,
                            &scaled_kbp,
                            gbp,
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
                let (q_aln, s_aln) = if let (Some(ref qs), Some(ref ss)) = (&ph.qseq, &ph.sseq) {
                    (qs.as_bytes().to_vec(), ss.as_bytes().to_vec())
                } else {
                    let q_aln: Vec<u8> = (0..ph.align_length as usize)
                        .map(|i| {
                            let idx = ph.query_start + i;
                            if idx < query_aa.len() {
                                ncbistdaa_to_ascii(query_aa[idx])
                            } else {
                                b'-'
                            }
                        })
                        .collect();
                    let s_aln: Vec<u8> = (0..ph.align_length as usize)
                        .map(|i| {
                            let idx = ph.subject_start + i;
                            if idx < sl {
                                ncbistdaa_to_ascii(subj_aa[idx])
                            } else {
                                b'-'
                            }
                        })
                        .collect();
                    (q_aln, s_aln)
                };
                let midline: Vec<u8> = q_aln
                    .iter()
                    .zip(s_aln.iter())
                    .map(|(&q, &s)| {
                        if q == s {
                            q
                        } else if q != b'-' && s != b'-' {
                            b'+'
                        } else {
                            b' '
                        }
                    })
                    .collect();
                Some(Hsp {
                    score: ph.score,
                    bit_score: prot_kbp.raw_to_bit(ph.score),
                    evalue,
                    query_start: ph.query_start,
                    query_end: ph.query_end,
                    subject_start: ph.subject_start,
                    subject_end: ph.subject_end,
                    num_identities: ph.num_ident as usize,
                    num_gaps: ph.gap_opens as usize,
                    alignment_length: ph.align_length as usize,
                    query_aln: q_aln,
                    midline,
                    subject_aln: s_aln,
                    query_frame: 0,
                    subject_frame: 0,
                })
            })
            .collect();

        if hsps.is_empty() {
            return None;
        }
        let hsps = if let Some(max) = max_hsps {
            hsps.into_iter().take(max).collect()
        } else {
            hsps
        };
        Some(SearchResult {
            subject_oid: oid,
            subject_title: title,
            subject_accession: accession,
            subject_len: sl,
            hsps,
            taxids: vec![],
        })
    };

    // Run sequentially or in parallel depending on num_threads/pool.
    let mut results: Vec<SearchResult> = if params.thread_pool.is_none() && params.num_threads == 1
    {
        (0..db.num_oids).filter_map(search_oid).collect()
    } else if let Some(pool) = params.thread_pool.as_deref() {
        use rayon::prelude::*;
        pool.install(|| {
            (0..db.num_oids)
                .into_par_iter()
                .filter_map(search_oid)
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
        pool.install(|| {
            (0..db.num_oids)
                .into_par_iter()
                .filter_map(search_oid)
                .collect()
        })
    };

    results.sort_by(compare_search_results);
    if results.len() > params.max_target_seqs {
        results.truncate(params.max_target_seqs);
    }
    results
}

/// Run a batch blastp search: multiple protein queries vs one protein database.
///
/// This is much more efficient than calling `blastp()` per query because
/// subjects are scanned once and checked against all query lookup tables.
/// Each subject is loaded into cache once, then all queries are checked.
///
/// Returns one `Vec<SearchResult>` per query, in the same order as `queries`.
pub fn blastp_batch(
    db: &BlastDb,
    queries: &[&[u8]],
    params: &SearchParams,
) -> Vec<Vec<SearchResult>> {
    if queries.is_empty() {
        return Vec::new();
    }

    let matrix = *get_matrix(params.matrix);
    let word_size = params.word_size.clamp(2, 6);
    let threshold = crate::stat::BLAST_WORD_THRESHOLD_BLASTP;

    let prot_kbp = protein_kbp_for_matrix(params.matrix, params.gap_open, params.gap_extend);

    let total_subj_len: usize = (0..db.num_oids)
        .map(|oid| db.get_seq_len(oid) as usize)
        .sum();

    let max_hsps = params.max_hsps;
    let evalue_threshold = params.evalue_threshold;
    // Convert bit-score parameters to raw (same as blastp()):
    // ungapped uses `ceil()+Int4` with UNGAPPED KBP (blast_parameters.c:221),
    // gapped uses plain `(Int4)` truncation with GAPPED KBP
    // (blast_parameters.c:457-463). `gap_trigger` is
    // `(Int4)((bits*ln2 + ungapped_log_k) / ungapped_lambda)` per
    // blast_parameters.c:343-344. Each query gets its own ungapped kbp and
    // hence its own gap_trigger (NCBI computes kbp_std per-context — the small
    // composition-driven drift is what makes seed cutoffs match bit-for-bit;
    // see iter 99 fix in `blastp()` for the boundary case it resolves).
    let ln2_b = crate::math::NCBIMATH_LN2;
    let _x_drop_gapped = (params.x_drop_gapped as f64 * ln2_b / prot_kbp.lambda) as i32;
    let x_drop_final = (params.x_drop_final as f64 * ln2_b / prot_kbp.lambda) as i32;
    let gap_open = params.gap_open;
    let gap_extend = params.gap_extend;

    // Prepare all queries: encode, build lookup tables, compute search spaces
    struct PreparedQuery {
        aa: Vec<u8>,
        lookup: crate::protein_lookup::ProteinLookupTable,
        search_space: f64,
        gap_trigger_raw: i32,
        x_drop_ungapped: i32,
    }
    let prepared: Vec<PreparedQuery> = queries
        .iter()
        .map(|q| {
            let aa: Vec<u8> = q
                .iter()
                .map(|&b| AMINOACID_TO_NCBISTDAA[b as usize & 0x7F])
                .collect();
            let ungapped_kbp = crate::stat::query_specific_protein_ungapped_kbp(&aa, &matrix);
            let x_drop_ungapped =
                (params.x_drop_ungapped as f64 * ln2_b / ungapped_kbp.lambda).ceil() as i32;
            let gap_trigger_raw = ((crate::stat::BLAST_GAP_TRIGGER_PROT * ln2_b
                + ungapped_kbp.log_k)
                / ungapped_kbp.lambda) as i32;
            let len_adj = crate::stat::compute_length_adjustment(
                aa.len() as i32,
                total_subj_len as i64,
                db.num_oids as i32,
                &prot_kbp,
            );
            let search_space = crate::stat::compute_search_space(
                aa.len() as i64,
                total_subj_len as i64,
                db.num_oids as i32,
                len_adj,
            );
            let lookup = crate::protein_lookup::ProteinLookupTable::build(
                &aa, word_size, &matrix, threshold,
            );
            PreparedQuery {
                aa,
                lookup,
                search_space,
                gap_trigger_raw,
                x_drop_ungapped,
            }
        })
        .collect();

    let num_queries = queries.len();

    // Build merged PV — bitwise OR of all query PVs.
    // Subject positions that don't match the merged PV can't match ANY query.
    let table_refs: Vec<&crate::protein_lookup::ProteinLookupTable> =
        prepared.iter().map(|pq| &pq.lookup).collect();
    let merged_pv = crate::protein_lookup::merge_pv(&table_refs);

    let query_refs: Vec<&[u8]> = prepared.iter().map(|pq| pq.aa.as_slice()).collect();

    // Shared scan x_dropoff: take the MAX across queries so the batched scan
    // never rejects a hit that a per-query NCBI run would have accepted.
    // (Per-query gap_trigger still filters seeds in the gapped phase below.)
    let x_drop_ungapped = prepared
        .iter()
        .map(|pq| pq.x_drop_ungapped)
        .max()
        .unwrap_or(0);

    // Process one subject: use merged PV scan, then gapped alignment on hits.
    let process_oid = |oid: u32| -> Vec<(usize, SearchResult)> {
        let subj_raw = db.get_sequence(oid);
        let subj_len = db.get_seq_len(oid) as usize;
        if subj_len < word_size {
            return Vec::new();
        }

        let subj_aa: Vec<u8> = subj_raw.iter().filter(|&&b| b != 0).copied().collect();
        if subj_aa.is_empty() {
            return Vec::new();
        }

        // Batch scan: one pass over subject, checks merged PV first
        let batch_ungapped = crate::protein_lookup::batch_scan_subject(
            &query_refs,
            &table_refs,
            &merged_pv,
            &subj_aa,
            &matrix,
            x_drop_ungapped,
        );

        let mut hits_for_queries: Vec<(usize, SearchResult)> = Vec::new();

        for (qi, ungapped) in batch_ungapped {
            let pq = &prepared[qi];

            // Phase 2: gapped extension on best seed.
            // Per-query gap_trigger from query-specific kbp_std — see iter 99.
            let ungap_cutoff = pq.gap_trigger_raw;
            let best_seed = ungapped
                .iter()
                .filter(|uh| uh.score >= ungap_cutoff)
                .max_by_key(|uh| uh.score);

            let mut phits = Vec::new();
            if let Some(uh) =
                best_seed.filter(|uh| prot_kbp.raw_to_evalue(uh.score, pq.search_space) < 10.0)
            {
                let seed_q = (uh.query_start + uh.query_end) / 2;
                let seed_s = (uh.subject_start + uh.subject_end) / 2;
                if let Some(gr) = crate::protein::protein_gapped_align(
                    &pq.aa,
                    &subj_aa,
                    seed_q,
                    seed_s,
                    &matrix,
                    gap_open,
                    gap_extend,
                    x_drop_final,
                ) {
                    let q_slice = &pq.aa[gr.query_start..gr.query_end];
                    let s_slice = &subj_aa[gr.subject_start..gr.subject_end];
                    let (qs, ss) = gr.edit_script.render_alignment(
                        q_slice,
                        s_slice,
                        crate::protein::ncbistdaa_to_char,
                    );
                    phits.push(crate::protein_lookup::ProteinHit {
                        query_start: gr.query_start,
                        query_end: gr.query_end,
                        subject_start: gr.subject_start,
                        subject_end: gr.subject_end,
                        score: gr.score,
                        num_ident: gr.num_ident,
                        align_length: gr.align_length,
                        mismatches: gr.mismatches,
                        gap_opens: gr.gap_opens,
                        qseq: Some(qs),
                        sseq: Some(ss),
                        scaled_score: None,
                        gapped_start_q: seed_q,
                        gapped_start_s: seed_s,
                    });
                }
            }
            if phits.is_empty() {
                phits = ungapped;
            }
            if phits.is_empty() {
                continue;
            }

            let best_ev = phits
                .iter()
                .map(|ph| prot_kbp.raw_to_evalue(ph.score, pq.search_space))
                .fold(f64::MAX, f64::min);
            if best_ev > evalue_threshold {
                continue;
            }

            let accession = db
                .get_accession(oid)
                .unwrap_or_else(|| format!("oid_{}", oid));
            let title = String::from_utf8_lossy(db.get_header(oid)).to_string();
            let sl = subj_aa.len();

            let hsps: Vec<Hsp> = phits
                .iter()
                .filter_map(|ph| {
                    let evalue = prot_kbp.raw_to_evalue(ph.score, pq.search_space);
                    if evalue > evalue_threshold {
                        return None;
                    }
                    let (q_aln, s_aln) = if let (Some(ref qs), Some(ref ss)) = (&ph.qseq, &ph.sseq)
                    {
                        (qs.as_bytes().to_vec(), ss.as_bytes().to_vec())
                    } else {
                        let qa: Vec<u8> = (0..ph.align_length as usize)
                            .map(|i| {
                                let idx = ph.query_start + i;
                                if idx < pq.aa.len() {
                                    ncbistdaa_to_ascii(pq.aa[idx])
                                } else {
                                    b'-'
                                }
                            })
                            .collect();
                        let sa: Vec<u8> = (0..ph.align_length as usize)
                            .map(|i| {
                                let idx = ph.subject_start + i;
                                if idx < sl {
                                    ncbistdaa_to_ascii(subj_aa[idx])
                                } else {
                                    b'-'
                                }
                            })
                            .collect();
                        (qa, sa)
                    };
                    let midline: Vec<u8> = q_aln
                        .iter()
                        .zip(s_aln.iter())
                        .map(|(&q, &s)| {
                            if q == s {
                                q
                            } else if q != b'-' && s != b'-' {
                                b'+'
                            } else {
                                b' '
                            }
                        })
                        .collect();
                    Some(Hsp {
                        score: ph.score,
                        bit_score: prot_kbp.raw_to_bit(ph.score),
                        evalue,
                        query_start: ph.query_start,
                        query_end: ph.query_end,
                        subject_start: ph.subject_start,
                        subject_end: ph.subject_end,
                        num_identities: ph.num_ident as usize,
                        num_gaps: ph.gap_opens as usize,
                        alignment_length: ph.align_length as usize,
                        query_aln: q_aln,
                        midline,
                        subject_aln: s_aln,
                        query_frame: 0,
                        subject_frame: 0,
                    })
                })
                .collect();

            if hsps.is_empty() {
                continue;
            }
            let hsps = if let Some(max) = max_hsps {
                hsps.into_iter().take(max).collect()
            } else {
                hsps
            };
            hits_for_queries.push((
                qi,
                SearchResult {
                    subject_oid: oid,
                    subject_title: title,
                    subject_accession: accession,
                    subject_len: sl,
                    hsps,
                    taxids: vec![],
                },
            ));
        }

        hits_for_queries
    };

    // Dispatch: sequential or parallel over subjects.
    let all_hits: Vec<Vec<(usize, SearchResult)>> =
        if params.thread_pool.is_none() && params.num_threads == 1 {
            (0..db.num_oids).map(process_oid).collect()
        } else if let Some(pool) = params.thread_pool.as_deref() {
            use rayon::prelude::*;
            pool.install(|| (0..db.num_oids).into_par_iter().map(process_oid).collect())
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
            pool.install(|| (0..db.num_oids).into_par_iter().map(process_oid).collect())
        };

    // Scatter results into per-query buckets
    let mut results: Vec<Vec<SearchResult>> = vec![Vec::new(); num_queries];
    for oid_hits in all_hits {
        for (qi, sr) in oid_hits {
            results[qi].push(sr);
        }
    }

    // Sort each query's results by e-value and truncate
    for r in &mut results {
        r.sort_by(compare_search_results);
        if r.len() > params.max_target_seqs {
            r.truncate(params.max_target_seqs);
        }
    }

    results
}

/// Run a blastn search (nucleotide query vs nucleotide database).
pub fn blastn_search(db: &BlastDb, query: &[u8], params: &SearchParams) -> Vec<SearchResult> {
    if query.is_empty() {
        return Vec::new();
    }

    // Encode query to BLASTNA
    let query_plus: Vec<u8> = query
        .iter()
        .map(|&b| IUPACNA_TO_BLASTNA[b as usize & 0x7F])
        .collect();
    let query_minus: Vec<u8> = crate::sequence::reverse_complement(&query_plus);

    let reward = params.match_score;
    let penalty = params.mismatch;

    let ungapped_kbp = crate::stat::KarlinBlk {
        lambda: 1.28,
        k: 0.46,
        log_k: 0.46_f64.ln(),
        h: 0.85,
        round_down: false,
    };
    let kbp = crate::stat::nucl_gapped_kbp_lookup(
        params.gap_open,
        params.gap_extend,
        reward,
        penalty,
        &ungapped_kbp,
    )
    .map(|(k, _)| k)
    .unwrap_or(ungapped_kbp);

    let total_subj_len = db.total_length;
    let len_adj = crate::stat::compute_length_adjustment(
        query.len() as i32,
        total_subj_len as i64,
        db.num_oids as i32,
        &kbp,
    );
    let search_space = crate::stat::compute_search_space(
        query.len() as i64,
        total_subj_len as i64,
        db.num_oids as i32,
        len_adj,
    );

    let x_dropoff = params.x_drop_ungapped;

    let (q_plus, q_minus) = match params.strand.as_str() {
        "plus" => (query_plus.as_slice(), &[] as &[u8]),
        "minus" => (&[] as &[u8], query_minus.as_slice()),
        _ => (query_plus.as_slice(), query_minus.as_slice()),
    };
    let prepared_query =
        crate::search::PreparedBlastnQuery::new_megablast(q_plus, q_minus, params.word_size);

    let search_oid = |oid: u32| -> Option<SearchResult> {
        let subject_packed = db.get_sequence(oid);
        let subject_len = db.get_seq_len(oid) as usize;
        if subject_len < params.word_size {
            return None;
        }

        let mut last_hit_scratch = prepared_query.last_hit_scratch();
        let mut hsps = crate::search::blastn_gapped_search_packed_prepared_with_xdrops(
            &prepared_query,
            q_plus,
            q_minus,
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
            &mut last_hit_scratch,
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
                midline: build_midline(
                    h.qseq.as_deref().unwrap_or(""),
                    h.sseq.as_deref().unwrap_or(""),
                ),
                query_frame: if h.context == 1 { -1 } else { 1 },
                subject_frame: 0,
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

    let mut results: Vec<SearchResult> = if params.thread_pool.is_none() && params.num_threads == 1
    {
        (0..db.num_oids).filter_map(search_oid).collect()
    } else if let Some(pool) = params.thread_pool.as_deref() {
        use rayon::prelude::*;
        pool.install(|| {
            (0..db.num_oids)
                .into_par_iter()
                .filter_map(search_oid)
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
        pool.install(|| {
            (0..db.num_oids)
                .into_par_iter()
                .filter_map(search_oid)
                .collect()
        })
    };

    results.sort_by(compare_search_results);
    if results.len() > params.max_target_seqs {
        results.truncate(params.max_target_seqs);
    }
    results
}

/// Run a blastx search (translated nucleotide query vs protein database).
pub fn blastx(db: &BlastDb, query: &[u8], params: &SearchParams) -> Vec<SearchResult> {
    if query.len() < 3 {
        return Vec::new();
    }
    let query_ncbi4na = ascii_to_ncbi4na(query);
    let code = crate::util::lookup_genetic_code(params.query_gencode);
    // 1-1 with blast_engine.c:s_BlastSearchEngineCore: call
    // BLAST_GetAllTranslations once, then iterate the 6 contexts via the
    // returned (translation_buffer, frame_offsets) pair.
    let (mut translation_buffer, frame_offsets) =
        crate::util::blast_get_all_translations_ncbi4na(&query_ncbi4na, query_ncbi4na.len(), code);
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

    let ln2 = crate::math::NCBIMATH_LN2;
    let ungapped_kbp = crate::stat::protein_ungapped_kbp();
    let x_drop_ungapped = (params.x_drop_ungapped as f64 * ln2 / ungapped_kbp.lambda).ceil() as i32;
    let x_drop_gapped = (params.x_drop_gapped as f64 * ln2 / prot_kbp.lambda) as i32;
    let x_drop_final = (params.x_drop_final as f64 * ln2 / prot_kbp.lambda) as i32;
    let gap_trigger_raw = ((crate::stat::BLAST_GAP_TRIGGER_PROT * ln2 + ungapped_kbp.log_k)
        / ungapped_kbp.lambda) as i32;

    let total_subj_len: usize = (0..db.num_oids)
        .map(|oid| db.get_seq_len(oid) as usize)
        .sum();
    let min_subject_length = (0..db.num_oids)
        .map(|oid| db.get_seq_len(oid) as i32)
        .filter(|&len| len > 0)
        .min()
        .unwrap_or(1);
    let avg_subject_length = (total_subj_len / db.num_oids.max(1) as usize).max(1);
    let gumbel_blk = protein_gumbel_for_matrix(
        params.matrix,
        params.gap_open,
        params.gap_extend,
        total_subj_len as i64,
    );

    let word_size = params.word_size.clamp(2, 6);
    let threshold = params
        .word_threshold
        .unwrap_or_else(|| suggested_word_threshold(params.matrix, crate::program::BLASTX));
    let max_hsps = params.max_hsps;

    // Build a 6-context QueryInfo from the translation offsets and call
    // BLAST_CalcEffLengths once (mirrors NCBI's setup-time invocation in
    // BLAST_GapAlignSetUp). The per-context `eff_searchsp` /
    // `length_adjustment` are then read inside the loop below.
    let mut query_info = crate::queryinfo::QueryInfo {
        num_queries: 1,
        contexts: (0..crate::util::NUM_FRAMES)
            .map(|ctx| {
                let begin = (frame_offsets[ctx] + 1) as usize;
                let end = frame_offsets[ctx + 1] as usize;
                let qlen = end.saturating_sub(begin) as i32;
                crate::queryinfo::ContextInfo {
                    query_offset: begin as i32,
                    query_length: qlen,
                    eff_searchsp: 0,
                    length_adjustment: 0,
                    query_index: 0,
                    frame: crate::util::blast_context_to_frame_blastx(ctx as u32),
                    is_valid: qlen > 0,
                }
            })
            .collect(),
        max_length: query_ncbi4na.len() as u32,
    };
    let scoring_options = crate::options::ScoringOptions {
        reward: 0,
        penalty: 0,
        gap_open: params.gap_open,
        gap_extend: params.gap_extend,
        gapped_calculation: true,
        matrix_name: Some(matrix_name.to_string()),
        is_ooframe: false,
    };
    let eff_params = crate::parameters::EffectiveLengthsParameters {
        options: crate::options::EffectiveLengthsOptions::default(),
        real_db_length: total_subj_len as i64,
        real_num_seqs: db.num_oids as i32,
    };
    let kbp_array = vec![prot_kbp.clone(); crate::util::NUM_FRAMES];
    let kbp_std_array = vec![ungapped_kbp.clone(); crate::util::NUM_FRAMES];
    crate::blast_setup::blast_calc_eff_lengths(
        crate::program::BLASTX,
        &scoring_options,
        &eff_params,
        &kbp_array,
        &kbp_std_array,
        matrix_name,
        &mut query_info,
    );

    let mut results = vec![None; db.num_oids as usize];
    for ctx in 0..crate::util::NUM_FRAMES {
        let frame = crate::util::blast_context_to_frame_blastx(ctx as u32);
        let begin = (frame_offsets[ctx] + 1) as usize;
        let end = frame_offsets[ctx + 1] as usize;
        if begin >= end {
            continue;
        }
        let prot: &[u8] = &translation_buffer[begin..end];
        if prot.len() < word_size {
            continue;
        }
        if crate::composition::read_composition(prot, AA_SIZE).1 == 0 {
            continue;
        }
        let search_space = query_info.contexts[ctx].eff_searchsp.max(1) as f64;
        // Build lookup table once per frame (not per subject).
        let lookup_table =
            crate::protein_lookup::ProteinLookupTable::build(prot, word_size, &matrix, threshold);

        let blastx_seed_cutoff = protein_prelim_seed_cutoff(
            gap_trigger_raw,
            params.evalue_threshold,
            &prot_kbp,
            gumbel_blk.as_ref(),
            prot.len() as i32,
            min_subject_length,
            search_space,
        );
        for oid in 0..db.num_oids {
            let subj_raw = db.get_sequence(oid);
            let subj_len = db.get_seq_len(oid) as usize;
            if subj_len == 0 {
                continue;
            }
            let subj_aa = &subj_raw[..subj_len];

            let phits = protein_alignment_hits(
                prot,
                &subj_aa,
                &matrix,
                &lookup_table,
                x_drop_ungapped,
                params.gap_open,
                params.gap_extend,
                x_drop_gapped,
                x_drop_final,
                blastx_seed_cutoff,
            );
            if phits.is_empty() {
                continue;
            }
            let best_raw_ev = phits
                .iter()
                .map(|ph| {
                    if let Some(ref gbp) = gumbel_blk {
                        crate::stat::spouge_evalue_with_gap_decay(
                            ph.score,
                            &prot_kbp,
                            gbp,
                            prot.len() as i32,
                            subj_len as i32,
                        )
                    } else {
                        prot_kbp.raw_to_evalue(ph.score, search_space)
                    }
                })
                .fold(f64::INFINITY, f64::min);
            if best_raw_ev > params.evalue_threshold {
                continue;
            }
            let accession = db
                .get_accession(oid)
                .unwrap_or_else(|| format!("oid_{}", oid));
            let title = String::from_utf8_lossy(db.get_header(oid)).to_string();
            let cutoff_s_blastx = protein_eval_cutoff(
                params.evalue_threshold,
                &prot_kbp,
                gumbel_blk.as_ref(),
                prot.len() as i32,
                avg_subject_length as i32,
                search_space,
            );
            let (final_phits, use_adj_matrix, lambda_ratio_opt) =
                apply_compositional_adjustment_per_subject(
                    prot,
                    &subj_aa,
                    &matrix,
                    phits,
                    params.gap_open,
                    params.gap_extend,
                    x_drop_final,
                    params.comp_adjust,
                    cutoff_s_blastx,
                    prot_kbp.lambda,
                );
            for ph in final_phits {
                if ph.score < cutoff_s_blastx {
                    continue;
                }
                let bit_score = prot_kbp.raw_to_bit(ph.score);
                let (e_score_i32, e_lambda) = match ph.scaled_score {
                    Some(s) if params.comp_adjust > 0 => {
                        (s, prot_kbp.lambda / COMPO_ADJUST_SCALE_FACTOR)
                    }
                    _ => (ph.score, prot_kbp.lambda),
                };
                let e_kbp = crate::stat::KarlinBlk {
                    lambda: e_lambda,
                    k: prot_kbp.k,
                    log_k: prot_kbp.log_k,
                    h: prot_kbp.h,
                    round_down: prot_kbp.round_down,
                };
                let evalue = if let Some(ref gbp) = gumbel_blk {
                    let base_ev = crate::stat::spouge_evalue_with_gap_decay(
                        e_score_i32,
                        &e_kbp,
                        gbp,
                        prot.len() as i32,
                        subj_len as i32,
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
                        crate::stat::spouge_evalue_with_gap_decay(
                            e_score_i32,
                            &scaled_kbp,
                            gbp,
                            prot.len() as i32,
                            subj_len as i32,
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
                if evalue > params.evalue_threshold {
                    continue;
                }
                let (query_start, query_end) = crate::util::protein_to_oriented_nuc_coords(
                    ph.query_start,
                    ph.query_end,
                    frame,
                );
                let (q_aln, s_aln) = if let (Some(ref qs), Some(ref ss)) = (&ph.qseq, &ph.sseq) {
                    (qs.as_bytes().to_vec(), ss.as_bytes().to_vec())
                } else {
                    (Vec::new(), Vec::new())
                };
                let midline: Vec<u8> = q_aln
                    .iter()
                    .zip(s_aln.iter())
                    .map(|(&q, &s)| {
                        if q == s {
                            q
                        } else if q != b'-' && s != b'-' {
                            b'+'
                        } else {
                            b' '
                        }
                    })
                    .collect();
                push_hsp_for_subject(
                    &mut results,
                    oid,
                    &title,
                    &accession,
                    subj_len,
                    &[],
                    Hsp {
                        score: ph.score,
                        bit_score,
                        evalue,
                        query_start,
                        query_end,
                        subject_start: ph.subject_start,
                        subject_end: ph.subject_end,
                        num_identities: ph.num_ident as usize,
                        num_gaps: ph.gap_opens as usize,
                        alignment_length: ph.align_length as usize,
                        query_aln: q_aln,
                        midline,
                        subject_aln: s_aln,
                        query_frame: frame,
                        subject_frame: 0,
                    },
                );
            }
        }
    }

    let mut results: Vec<SearchResult> = results.into_iter().flatten().collect();
    if params.sum_stats {
        apply_tblastn_linked_sum_stats(&mut results, &query_info, &prot_kbp);
    }
    for result in &mut results {
        prune_translated_hsp_variants(&mut result.hsps);
        result.hsps.sort_by(compare_hsps_by_evalue_then_coords);
        if let Some(max) = max_hsps {
            result.hsps.truncate(max);
        }
    }
    results.sort_by(compare_search_results);
    if results.len() > params.max_target_seqs {
        results.truncate(params.max_target_seqs);
    }
    results
}

/// Run a tblastn search (protein query vs translated nucleotide database).
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
    if crate::composition::read_composition(&query_aa, AA_SIZE).1 == 0 {
        return Vec::new();
    }
    let matrix = *get_matrix(params.matrix);
    let matrix_name = protein_matrix_name(params.matrix);
    let word_size = params.word_size.clamp(2, 6);
    let threshold = params
        .word_threshold
        .unwrap_or_else(|| suggested_word_threshold(params.matrix, crate::program::TBLASTN));

    let prot_kbp = protein_kbp_for_matrix(params.matrix, params.gap_open, params.gap_extend);
    let ln2 = crate::math::NCBIMATH_LN2;
    // tblastn intentionally keeps the IDEAL ungapped kbp here.
    // Switching to query-specific kbp (as we did for blastp in iter 99)
    // changes seed selection and re-routes the gapped DP, producing
    // higher scores than NCBI on real fixtures (e.g. NP_982592 vs seqn:
    // ideal=356 matches NCBI's 356; query-specific=403 diverges). Until
    // we trace the seed-selection sensitivity, leave tblastn on the ideal.
    let ungapped_kbp = crate::stat::protein_ungapped_kbp();
    let x_drop_ungapped = (params.x_drop_ungapped as f64 * ln2 / ungapped_kbp.lambda).ceil() as i32;
    let x_drop_gapped = (params.x_drop_gapped as f64 * ln2 / prot_kbp.lambda) as i32;
    let x_drop_final = (params.x_drop_final as f64 * ln2 / prot_kbp.lambda) as i32;
    let gap_trigger_raw = ((crate::stat::BLAST_GAP_TRIGGER_PROT * ln2 + ungapped_kbp.log_k)
        / ungapped_kbp.lambda) as i32;
    let total_subj_len: usize = (0..db.num_oids)
        .map(|oid| db.get_seq_len(oid) as usize)
        .sum();
    // NCBI uses the sequence-source minimum length for `cutoff_score_max`
    // (`blast_setup.c:970`). For translated-subject programs, the nucleotide
    // default BLAST_SEQSRC_MINLENGTH=10 is divided by three
    // (`blast_setup.c:970-973`), truncating to 3.
    const BLAST_SEQSRC_MINLENGTH_TBLASTN: i32 = 10;
    let min_subject_length: i32 = BLAST_SEQSRC_MINLENGTH_TBLASTN / 3;
    // BLAST_CalcEffLengths handles `db_length /= 3` internally for tblastn,
    // so we pass the raw total. For Gumbel-block lookup we still need the
    // translated total ourselves.
    let translated_total_subj_len = (total_subj_len / 3).max(1);
    let avg_subject_length = (translated_total_subj_len / db.num_oids.max(1) as usize).max(1);
    let mut query_info = crate::queryinfo::QueryInfo {
        num_queries: 1,
        contexts: vec![crate::queryinfo::ContextInfo {
            query_offset: 0,
            query_length: query_aa.len() as i32,
            eff_searchsp: 0,
            length_adjustment: 0,
            query_index: 0,
            frame: 0,
            is_valid: !query_aa.is_empty(),
        }],
        max_length: query_aa.len() as u32,
    };
    let scoring_options = crate::options::ScoringOptions {
        reward: 0,
        penalty: 0,
        gap_open: params.gap_open,
        gap_extend: params.gap_extend,
        gapped_calculation: true,
        matrix_name: Some(matrix_name.to_string()),
        is_ooframe: false,
    };
    let eff_params = crate::parameters::EffectiveLengthsParameters {
        options: crate::options::EffectiveLengthsOptions::default(),
        real_db_length: total_subj_len as i64,
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
    let gumbel_blk = protein_gumbel_for_matrix(
        params.matrix,
        params.gap_open,
        params.gap_extend,
        translated_total_subj_len as i64,
    );
    let max_hsps = params.max_hsps;

    // Build lookup table once per query.
    let lookup_table =
        crate::protein_lookup::ProteinLookupTable::build(&query_aa, word_size, &matrix, threshold);

    let results = map_database_oids(db, params, |oid| {
        let mut result: Option<SearchResult> = None;
        let subject_packed = db.get_sequence(oid);
        let subject_len = db.get_seq_len(oid) as usize;

        // Decode subject to BLASTNA (with ambiguity overlay when the .nsq
        // ambiguity table is present), then convert per-byte to NCBI4na for
        // `BLAST_GetAllTranslations`. The BLASTNA → NCBI4na step is the
        // existing `BLASTNA_TO_NCBI4NA` table from `blast_encoding.h:84`.
        let subj_blastna: Vec<u8> = match db.get_ambiguity_data(oid) {
            Some(amb) => crate::search::decode_packed_ncbi2na_with_ambiguity(
                subject_packed,
                subject_len,
                amb,
            ),
            None => crate::search::decode_packed_ncbi2na(subject_packed, subject_len),
        };
        let subj_ncbi4na: Vec<u8> = subj_blastna
            .iter()
            .map(|&b| crate::encoding::BLASTNA_TO_NCBI4NA[b as usize])
            .collect();

        let (translation_buffer, frame_offsets) = crate::util::blast_get_all_translations_ncbi4na(
            &subj_ncbi4na,
            subj_ncbi4na.len(),
            crate::util::lookup_genetic_code(params.db_gencode),
        );

        for ctx in 0..crate::util::NUM_FRAMES {
            let frame = crate::util::blast_context_to_frame_blastx(ctx as u32);
            let begin = (frame_offsets[ctx] + 1) as usize;
            let end = frame_offsets[ctx + 1] as usize;
            if begin >= end {
                continue;
            }
            let prot: &[u8] = &translation_buffer[begin..end];
            if prot.len() < word_size {
                continue;
            }
            if crate::composition::read_composition(prot, AA_SIZE).1 == 0 {
                continue;
            }
            let subj_prot_len = prot.len();
            let tblastn_seed_cutoff = protein_prelim_seed_cutoff(
                gap_trigger_raw,
                params.evalue_threshold,
                &prot_kbp,
                gumbel_blk.as_ref(),
                query_aa.len() as i32,
                min_subject_length,
                search_space,
            );
            let phits = protein_alignment_hits(
                &query_aa,
                prot,
                &matrix,
                &lookup_table,
                x_drop_ungapped,
                params.gap_open,
                params.gap_extend,
                x_drop_gapped,
                x_drop_final,
                tblastn_seed_cutoff,
            );
            if phits.is_empty() {
                continue;
            }
            let best_raw_ev = phits
                .iter()
                .map(|ph| {
                    if let Some(ref gbp) = gumbel_blk {
                        crate::stat::spouge_evalue_with_gap_decay(
                            ph.score,
                            &prot_kbp,
                            gbp,
                            query_aa.len() as i32,
                            subj_prot_len as i32,
                        )
                    } else {
                        prot_kbp.raw_to_evalue(ph.score, search_space)
                    }
                })
                .fold(f64::INFINITY, f64::min);
            if best_raw_ev > params.evalue_threshold {
                continue;
            }
            let accession = db
                .get_accession(oid)
                .unwrap_or_else(|| format!("oid_{}", oid));
            let title = String::from_utf8_lossy(db.get_header(oid)).to_string();
            let cutoff_s_tblastn = protein_eval_cutoff(
                params.evalue_threshold,
                &prot_kbp,
                gumbel_blk.as_ref(),
                query_aa.len() as i32,
                avg_subject_length as i32,
                search_space,
            );
            let (final_phits, use_adj_matrix, lambda_ratio_opt) =
                apply_compositional_adjustment_per_subject(
                    &query_aa,
                    prot,
                    &matrix,
                    phits,
                    params.gap_open,
                    params.gap_extend,
                    x_drop_final,
                    params.comp_adjust,
                    cutoff_s_tblastn,
                    prot_kbp.lambda,
                );
            for ph in final_phits {
                if ph.score < cutoff_s_tblastn {
                    continue;
                }
                let bit_score = prot_kbp.raw_to_bit(ph.score);
                let (e_score_i32, e_lambda) = match ph.scaled_score {
                    Some(s) if params.comp_adjust > 0 => {
                        (s, prot_kbp.lambda / COMPO_ADJUST_SCALE_FACTOR)
                    }
                    _ => (ph.score, prot_kbp.lambda),
                };
                let e_kbp = crate::stat::KarlinBlk {
                    lambda: e_lambda,
                    k: prot_kbp.k,
                    log_k: prot_kbp.log_k,
                    h: prot_kbp.h,
                    round_down: prot_kbp.round_down,
                };
                let evalue = if let Some(ref gbp) = gumbel_blk {
                    let base_ev = crate::stat::spouge_evalue_with_gap_decay(
                        e_score_i32,
                        &e_kbp,
                        gbp,
                        query_aa.len() as i32,
                        subj_prot_len as i32,
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
                        crate::stat::spouge_evalue_with_gap_decay(
                            e_score_i32,
                            &scaled_kbp,
                            gbp,
                            query_aa.len() as i32,
                            subj_prot_len as i32,
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
                if evalue > params.evalue_threshold {
                    continue;
                }
                let (subject_start, subject_end) = crate::util::protein_to_oriented_nuc_coords(
                    ph.subject_start,
                    ph.subject_end,
                    frame,
                );
                let (q_aln, s_aln) = if let (Some(ref qs), Some(ref ss)) = (&ph.qseq, &ph.sseq) {
                    (qs.as_bytes().to_vec(), ss.as_bytes().to_vec())
                } else {
                    (Vec::new(), Vec::new())
                };
                let midline: Vec<u8> = q_aln
                    .iter()
                    .zip(s_aln.iter())
                    .map(|(&q, &s)| {
                        if q == s {
                            q
                        } else if q != b'-' && s != b'-' {
                            b'+'
                        } else {
                            b' '
                        }
                    })
                    .collect();
                let hsp = Hsp {
                    score: ph.score,
                    bit_score,
                    evalue,
                    query_start: ph.query_start,
                    query_end: ph.query_end,
                    subject_start,
                    subject_end,
                    num_identities: ph.num_ident as usize,
                    num_gaps: ph.gap_opens as usize,
                    alignment_length: ph.align_length as usize,
                    query_aln: q_aln,
                    midline,
                    subject_aln: s_aln,
                    query_frame: 0,
                    subject_frame: frame,
                };
                match &mut result {
                    Some(existing) => existing.hsps.push(hsp),
                    None => {
                        result = Some(SearchResult {
                            subject_oid: oid,
                            subject_title: title.clone(),
                            subject_accession: accession.clone(),
                            subject_len,
                            hsps: vec![hsp],
                            taxids: vec![],
                        });
                    }
                }
            }
        }
        result
    });

    let mut results: Vec<SearchResult> = results.into_iter().flatten().collect();
    for result in &mut results {
        prune_translated_hsp_variants(&mut result.hsps);
        result.hsps.sort_by(compare_hsps_by_evalue_then_coords);
        if let Some(max) = max_hsps {
            result.hsps.truncate(max);
        }
    }
    results.sort_by(compare_search_results);
    if results.len() > params.max_target_seqs {
        results.truncate(params.max_target_seqs);
    }
    results
}

/// Run a tblastx search (translated nt query vs translated nt database).
pub fn tblastx(db: &BlastDb, query: &[u8], params: &SearchParams) -> Vec<SearchResult> {
    if query.len() < 3 {
        return Vec::new();
    }
    let query_ncbi4na = ascii_to_ncbi4na(query);
    let q_code = crate::util::lookup_genetic_code(params.query_gencode);
    let (mut query_translation, query_offsets) = crate::util::blast_get_all_translations_ncbi4na(
        &query_ncbi4na,
        query_ncbi4na.len(),
        q_code,
    );
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
    let word_size = params.word_size.clamp(2, 6);
    let threshold = params
        .word_threshold
        .unwrap_or_else(|| suggested_word_threshold(params.matrix, crate::program::TBLASTX));

    let ln2 = crate::math::NCBIMATH_LN2;
    let ungapped_kbp = crate::stat::protein_ungapped_kbp();
    let x_drop_ungapped = (params.x_drop_ungapped as f64 * ln2 / ungapped_kbp.lambda).ceil() as i32;
    let total_subj_len: usize = (0..db.num_oids)
        .map(|oid| db.get_seq_len(oid) as usize)
        .sum();
    let translated_total_subj_len = (total_subj_len / 3).max(1);

    // BLAST_CalcEffLengths once for all 6 query contexts; produces
    // per-context eff_searchsp + length_adjustment that the search loop
    // reads directly. Mirrors NCBI's setup-time invocation in
    // blast_engine.c. db_length /= 3 happens inside CalcEffLengths because
    // tblastx's subject_is_translated is true.
    let mut query_info_calc = crate::queryinfo::QueryInfo {
        num_queries: 1,
        contexts: (0..crate::util::NUM_FRAMES)
            .map(|ctx| {
                let begin = (query_offsets[ctx] + 1) as usize;
                let end = query_offsets[ctx + 1] as usize;
                let qlen = end.saturating_sub(begin) as i32;
                crate::queryinfo::ContextInfo {
                    query_offset: begin as i32,
                    query_length: qlen,
                    eff_searchsp: 0,
                    length_adjustment: 0,
                    query_index: 0,
                    frame: crate::util::blast_context_to_frame_blastx(ctx as u32),
                    is_valid: qlen > 0,
                }
            })
            .collect(),
        max_length: query_ncbi4na.len() as u32,
    };
    // tblastx is ungapped-only (`blast_options.c:869`). NCBI propagates
    // `gapped_calculation = FALSE` into `BLAST_CalcEffLengths`, which routes
    // alpha/beta through the ungapped lookup (`blast_setup.c:806`).
    let scoring_options = crate::options::ScoringOptions {
        reward: 0,
        penalty: 0,
        gap_open: params.gap_open,
        gap_extend: params.gap_extend,
        gapped_calculation: false,
        matrix_name: Some(matrix_name.to_string()),
        is_ooframe: false,
    };
    let eff_params = crate::parameters::EffectiveLengthsParameters {
        options: crate::options::EffectiveLengthsOptions::default(),
        real_db_length: total_subj_len as i64,
        real_num_seqs: db.num_oids as i32,
    };
    let kbp_array = vec![ungapped_kbp.clone(); crate::util::NUM_FRAMES];
    let kbp_std_array = vec![ungapped_kbp.clone(); crate::util::NUM_FRAMES];
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

    let mut query_contexts = Vec::new();
    let mut results = vec![None; db.num_oids as usize];
    for oid in 0..db.num_oids {
        let subject_packed = db.get_sequence(oid);
        let subject_len = db.get_seq_len(oid) as usize;

        // BLASTNA decode (with ambiguity overlay if available), then per-byte
        // BLASTNA → NCBI4na for `BLAST_GetAllTranslations`.
        let subj_blastna: Vec<u8> = match db.get_ambiguity_data(oid) {
            Some(amb) => crate::search::decode_packed_ncbi2na_with_ambiguity(
                subject_packed,
                subject_len,
                amb,
            ),
            None => crate::search::decode_packed_ncbi2na(subject_packed, subject_len),
        };
        let subj_ncbi4na: Vec<u8> = subj_blastna
            .iter()
            .map(|&b| crate::encoding::BLASTNA_TO_NCBI4NA[b as usize])
            .collect();

        let (subj_translation, subj_offsets) = crate::util::blast_get_all_translations_ncbi4na(
            &subj_ncbi4na,
            subj_ncbi4na.len(),
            crate::util::lookup_genetic_code(params.db_gencode),
        );

        for q_ctx in 0..crate::util::NUM_FRAMES {
            let qframe = crate::util::blast_context_to_frame_blastx(q_ctx as u32);
            let q_begin = (query_offsets[q_ctx] + 1) as usize;
            let q_end = query_offsets[q_ctx + 1] as usize;
            if q_begin >= q_end {
                continue;
            }
            let q_prot: &[u8] = &query_translation[q_begin..q_end];
            if q_prot.len() < word_size {
                continue;
            }
            if crate::composition::read_composition(q_prot, AA_SIZE).1 == 0 {
                continue;
            }
            let len_adj = query_info_calc.contexts[q_ctx].length_adjustment;
            let search_space = query_info_calc.contexts[q_ctx].eff_searchsp.max(1) as f64;

            // Per-context ungapped Karlin params from this query frame's
            // amino-acid composition. Mirrors `Blast_ScoreBlkKbpUngappedCalc`
            // (`blast_stat.c:2737`): `Blast_ResFreqString` →
            // `BlastScoreFreqCalc` → `Blast_KarlinBlkUngappedCalc`. For
            // translated queries (blastx/tblastx/rps-tblastn) NCBI also
            // applies an "ideal-Lambda cap" (line 2796) that we mirror here.
            let mut lo = i32::MAX;
            let mut hi = i32::MIN;
            for i in 0..AA_SIZE {
                for j in 0..AA_SIZE {
                    let s = matrix[i][j];
                    if s != 0 {
                        if s < lo {
                            lo = s;
                        }
                        if s > hi {
                            hi = s;
                        }
                    }
                }
            }
            if lo == i32::MAX {
                lo = -1;
            }
            if hi == i32::MIN {
                hi = 1;
            }
            let std_freq = crate::stat::protein_std_freq_ncbistdaa();
            let m = matrix;
            let matrix_fn = move |i: usize, j: usize| -> i32 { m[i][j] };
            let ctx_kbp_results = crate::stat::ungapped_kbp_calc_with_std(
                q_prot,
                &[crate::stat::UngappedKbpContext {
                    query_offset: 0,
                    query_length: q_prot.len() as i32,
                    is_valid: true,
                }],
                lo,
                hi,
                AA_SIZE,
                &[0u8, 2, 21, 23, 24, 25, 26, 27], // gap, B, X, Z, U, *, O, J
                &std_freq,
                &matrix_fn,
            );
            let ideal = crate::stat::protein_ungapped_kbp();
            let ctx_kbp = match ctx_kbp_results.into_iter().next().flatten() {
                Some(mut k) => {
                    // For translated queries: substitute ideal when computed
                    // Lambda >= ideal (more conservative). `blast_stat.c:2796`.
                    if k.lambda >= ideal.lambda {
                        k = ideal.clone();
                    }
                    k
                }
                None => ideal.clone(),
            };
            let ctx_gap_trigger_raw = ((crate::stat::BLAST_GAP_TRIGGER_PROT * ln2 + ctx_kbp.log_k)
                / ctx_kbp.lambda) as i32;

            query_contexts.push(TranslatedContextStats {
                frame: qframe,
                query_length: q_prot.len() as i32,
                eff_searchsp: search_space.max(1.0) as i64,
                length_adjustment: len_adj,
                kbp: ctx_kbp.clone(),
            });
            let lookup_table = crate::protein_lookup::ProteinLookupTable::build(
                q_prot, word_size, &matrix, threshold,
            );

            for s_ctx in 0..crate::util::NUM_FRAMES {
                let sframe = crate::util::blast_context_to_frame_blastx(s_ctx as u32);
                let s_begin = (subj_offsets[s_ctx] + 1) as usize;
                let s_end = subj_offsets[s_ctx + 1] as usize;
                if s_begin >= s_end {
                    continue;
                }
                let s_prot: &[u8] = &subj_translation[s_begin..s_end];
                if s_prot.len() < word_size {
                    continue;
                }
                if crate::composition::read_composition(s_prot, AA_SIZE).1 == 0 {
                    continue;
                }
                let tblastx_seed_cutoff = protein_prelim_seed_cutoff(
                    ctx_gap_trigger_raw,
                    crate::stat::CUTOFF_E_TBLASTX,
                    &ctx_kbp,
                    gumbel_blk.as_ref(),
                    q_prot.len() as i32,
                    1,
                    search_space,
                );
                let hit_cutoff = protein_eval_cutoff(
                    params.evalue_threshold,
                    &ctx_kbp,
                    gumbel_blk.as_ref(),
                    q_prot.len() as i32,
                    1,
                    search_space,
                );
                let save_cutoff = tblastx_seed_cutoff.min(hit_cutoff);
                let ungapped_hits = crate::protein_lookup::protein_scan_with_table(
                    q_prot,
                    s_prot,
                    &matrix,
                    &lookup_table,
                    x_drop_ungapped,
                );
                for ph in ungapped_hits {
                    if ph.score < save_cutoff {
                        continue;
                    }
                    // tblastx is ungapped-only (`blast_options.c:869`): gapped
                    // calculation is rejected at option-validation time. NCBI
                    // therefore saves ungapped HSPs and uses `sbp->kbp`
                    // (per-context ungapped Karlin params, computed from this
                    // query frame's amino-acid composition) for bit-score and
                    // e-value. See `blast_hits.c:1928`
                    // (`Blast_HSPListGetBitScores`): `kbp[hsp->context]->Lambda`.
                    let bit_score = ctx_kbp.raw_to_bit(ph.score);
                    let evalue = if let Some(ref gbp) = gumbel_blk {
                        crate::stat::spouge_evalue(
                            ph.score,
                            &ctx_kbp,
                            gbp,
                            q_prot.len() as i32,
                            s_prot.len() as i32,
                        )
                    } else {
                        // tblastx ungapped path: NCBI's `s_BlastEvenGapLinkHSPs`
                        // (single-HSP branch in `BLAST_LargeGapSumE`,
                        // `blast_stat.c:4557`) divides the raw KarlinStoE result
                        // by `BLAST_GapDecayDivisor(BLAST_GAP_DECAY_RATE, 1)`
                        // = 0.5, multiplying the e-value by 2.
                        let weight_divisor =
                            crate::stat::gap_decay_divisor(crate::stat::BLAST_GAP_DECAY_RATE, 1);
                        ctx_kbp.raw_to_evalue(ph.score, search_space) / weight_divisor
                    };
                    if evalue > params.evalue_threshold {
                        continue;
                    }
                    let accession = db
                        .get_accession(oid)
                        .unwrap_or_else(|| format!("oid_{}", oid));
                    let title = String::from_utf8_lossy(db.get_header(oid)).to_string();
                    let (query_start, query_end) = crate::util::protein_to_oriented_nuc_coords(
                        ph.query_start,
                        ph.query_end,
                        qframe,
                    );
                    let (subject_start, subject_end) = crate::util::protein_to_oriented_nuc_coords(
                        ph.subject_start,
                        ph.subject_end,
                        sframe,
                    );

                    let q_aln: Vec<u8> = q_prot[ph.query_start..ph.query_end]
                        .iter()
                        .map(|&aa| crate::protein::ncbistdaa_to_char(aa) as u8)
                        .collect();
                    let s_aln: Vec<u8> = s_prot[ph.subject_start..ph.subject_end]
                        .iter()
                        .map(|&aa| crate::protein::ncbistdaa_to_char(aa) as u8)
                        .collect();
                    let midline: Vec<u8> = q_aln
                        .iter()
                        .zip(s_aln.iter())
                        .map(|(&q, &s)| {
                            if q == s {
                                q
                            } else if q != b'-' && s != b'-' {
                                b'+'
                            } else {
                                b' '
                            }
                        })
                        .collect();
                    push_hsp_for_subject(
                        &mut results,
                        oid,
                        &title,
                        &accession,
                        subject_len,
                        &[],
                        Hsp {
                            score: ph.score,
                            bit_score,
                            evalue,
                            query_start,
                            query_end,
                            subject_start,
                            subject_end,
                            num_identities: ph.num_ident as usize,
                            num_gaps: 0,
                            alignment_length: ph.align_length as usize,
                            query_aln: q_aln,
                            midline,
                            subject_aln: s_aln,
                            query_frame: qframe,
                            subject_frame: sframe,
                        },
                    );
                }
            }
        }
    }

    apply_tblastx_linked_sum_stats(&mut results, &query_contexts);

    let mut results: Vec<SearchResult> = results.into_iter().flatten().collect();
    for result in &mut results {
        result.hsps.sort_by(compare_tblastx_hsps);
        if let Some(max) = max_hsps {
            result.hsps.truncate(max);
        }
    }
    results.sort_by(compare_search_results);
    if results.len() > params.max_target_seqs {
        results.truncate(params.max_target_seqs);
    }
    results
}

// ── Utility functions ───────────────────────────────────────────────────────

type ProteinCompositionAdjustment = Option<(Option<[[i32; AA_SIZE]; AA_SIZE]>, Option<f64>)>;

fn protein_composition_adjustment(
    query_aa: &[u8],
    subj_aa: &[u8],
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
    comp_mode: u8,
) -> ProteinCompositionAdjustment {
    if comp_mode == 0 {
        return None;
    }

    let (qcomp28, qn) = crate::composition::read_composition(query_aa, AA_SIZE);
    let (scomp28, sn) = crate::composition::read_composition(subj_aa, AA_SIZE);
    if qn == 0 || sn == 0 {
        return None;
    }

    let mut qp20 = [0.0f64; 20];
    let mut sp20 = [0.0f64; 20];
    crate::compo_mode_condition::gather_letter_probs(&qcomp28, &mut qp20);
    crate::compo_mode_condition::gather_letter_probs(&scomp28, &mut sp20);

    let rule = crate::compo_mode_condition::choose_matrix_adjust_rule(
        query_aa.len(),
        subj_aa.len(),
        &qp20,
        &sp20,
        comp_mode,
    );

    use crate::compo_mode_condition::MatrixAdjustRule;
    match rule {
        MatrixAdjustRule::DontAdjust => None,
        MatrixAdjustRule::ScaleOldMatrix => {
            let ungapped_lambda =
                crate::stat::protein_ungapped_kbp().lambda / COMPO_ADJUST_SCALE_FACTOR;
            let freq_ratios = crate::matrix::get_blosum62_freq_ratios();
            let start_matrix =
                crate::composition::matrix_from_freq_ratios(ungapped_lambda, &freq_ratios);
            crate::composition::composition_scale_matrix(
                &start_matrix,
                &qcomp28,
                &scomp28,
                ungapped_lambda,
                &freq_ratios,
            )
            .map(|adj_mat| (Some(adj_mat), None))
        }
        MatrixAdjustRule::UserSpecifiedRelEntropy
        | MatrixAdjustRule::UnconstrainedRelEntropy
        | MatrixAdjustRule::RelEntropyOldMatrixNewContext
        | MatrixAdjustRule::RelEntropyOldMatrixOldContext => {
            let (joint_probs, first_std, second_std) = crate::composition::blosum62_workspace();
            let mut adj_matrix = *matrix;
            let ungapped_lambda =
                crate::stat::protein_ungapped_kbp().lambda / COMPO_ADJUST_SCALE_FACTOR;
            let freq_ratios = crate::matrix::get_blosum62_freq_ratios();
            let start_matrix =
                crate::composition::matrix_from_freq_ratios(ungapped_lambda, &freq_ratios);
            let status = crate::composition::composition_matrix_adj(
                &mut adj_matrix,
                AA_SIZE,
                rule,
                qn,
                sn,
                &qcomp28,
                &scomp28,
                20,
                0.44,
                &joint_probs,
                &first_std,
                &second_std,
                ungapped_lambda,
                &start_matrix,
            );
            if status == 0 {
                Some((Some(adj_matrix), None))
            } else {
                // 1-1 with `Blast_AdjustScores` (composition_adjustment.c:1501–1530):
                // when `Blast_CompositionMatrixAdj` returns a non-fatal failure
                // (status > 0 — typically optimization didn't converge), NCBI
                // resets `*matrix_adjust_rule = eCompoScaleOldMatrix` and falls
                // through to `Blast_CompositionBasedStats`, which rescales the
                // matrix in place using `lambda_ratio`. We mirror that fallback
                // by invoking `composition_scale_matrix` (== ScaleOldMatrix path)
                // instead of returning a bare `lambda_ratio` for evalue scaling.
                crate::composition::composition_scale_matrix(
                    &start_matrix,
                    &qcomp28,
                    &scomp28,
                    ungapped_lambda,
                    &freq_ratios,
                )
                .map(|adj_mat| (Some(adj_mat), None))
            }
        }
    }
}

/// Apply NCBI-compatible composition-based adjustment to a per-subject
/// HSP list. Mirrors the matrix-adjust + re-align block used by blastp
/// (see lines 1411–1479) so blastx / tblastn / tblastx follow the same
/// `Blast_RedoAlignmentCore_MT` semantics: with the adjusted matrix
/// the alignment is re-run via `protein_gapped_align` at scaled gap
/// penalties, and the re-aligned coordinates / scores replace the
/// original ones.
///
/// Returns `(final_phits, use_adj_matrix, lambda_ratio_opt)` where
/// `lambda_ratio_opt` is `Some(lr)` only when the matrix-optimization
/// path fell back to lambda-only adjustment (status != 0 in NCBI's
/// `Blast_CompositionMatrixAdj`); callers use it to scale the e-value.
fn apply_compositional_adjustment_per_subject(
    query_aa: &[u8],
    subj_aa: &[u8],
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
    phits: Vec<crate::protein_lookup::ProteinHit>,
    gap_open: i32,
    gap_extend: i32,
    x_drop_final: i32,
    comp_mode: u8,
    cutoff_s: i32,
    gapped_lambda: f64,
) -> (Vec<crate::protein_lookup::ProteinHit>, bool, Option<f64>) {
    if comp_mode == 0 {
        return (phits, false, None);
    }
    let comp_scale = COMPO_ADJUST_SCALE_FACTOR;
    let scaled_gap_open = crate::math::nint(gap_open as f64 * comp_scale) as i32;
    let scaled_gap_extend = crate::math::nint(gap_extend as f64 * comp_scale) as i32;
    let scaled_x_drop_final = crate::math::nint(x_drop_final as f64 * comp_scale) as i32;

    let redo_subj_storage =
        kappa_redo_subject_sequence(query_aa.len(), subj_aa, &phits, gapped_lambda);
    let redo_subj_aa: &[u8] = redo_subj_storage.as_ref();
    let adj_result = protein_composition_adjustment(query_aa, redo_subj_aa, matrix, comp_mode);
    let lambda_ratio_opt = adj_result.as_ref().and_then(|(_, lr)| *lr);

    let _ = cutoff_s;

    if crate::composition::read_composition(redo_subj_aa, AA_SIZE).1 == 0 {
        return (Vec::new(), true, None);
    }

    let (final_phits, use_adj_matrix) = if let Some((None, _)) = adj_result {
        (Vec::new(), true)
    } else if let Some((Some(ref adj_mat), _)) = adj_result {
        // 1-1 with the single-HSP arm of NCBI's `Blast_RedoAlignmentCore_MT`
        // (redo_alignment.c:1430-1530): full SW under the adjusted matrix
        // → reverse SW for the start → forward-only X-drop within the
        // SW-derived bounds via `protein_sw_bounded_xdrop_align`.
        let sw_bounded = if phits.len() == 1 {
            let (sw_score, m_end_sw, q_end_sw) =
                crate::smith_waterman::blast_smith_waterman_score_only(
                    redo_subj_aa,
                    query_aa,
                    adj_mat,
                    scaled_gap_open,
                    scaled_gap_extend,
                );
            if sw_score > 0 {
                let (_, m_start_sw, q_start_sw) =
                    crate::smith_waterman::blast_smith_waterman_find_start(
                        redo_subj_aa,
                        query_aa,
                        adj_mat,
                        scaled_gap_open,
                        scaled_gap_extend,
                        m_end_sw,
                        q_end_sw,
                        sw_score,
                    );
                let q_extent = q_end_sw - q_start_sw;
                let s_extent = m_end_sw - m_start_sw;
                Some((q_start_sw, m_start_sw, q_extent, s_extent, sw_score))
            } else {
                None
            }
        } else {
            None
        };
        let mut new_phits = Vec::new();
        // Run BOTH bounded SW and bidirectional X-drop per phit; pick
        // the bounded result only when it strictly outscores
        // bidirectional (mirrors the inline blastp logic at api.rs:1483).
        // Bounded helper avoids alignment-end overshoot, but for
        // equal-score paths the bidirectional X-drop tends to land on
        // the NCBI-matching identity-rich alignment.
        let bounded_gr_opt =
            if let Some((q_start, m_start, q_extent, s_extent, target_score)) = sw_bounded {
                crate::protein::protein_sw_bounded_xdrop_align(
                    query_aa,
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
        {
            for ph in &phits {
                // NCBI's `s_RedoOneAlignment` uses `hsp->query.gapped_start`
                // (seed) — NOT alignment `query.offset` — when re-running
                // rescaled gapped DP under composition adjustment.
                // (`blast_kappa.c:1924-1926`).
                let bidir_gr = crate::protein::protein_gapped_align(
                    query_aa,
                    redo_subj_aa,
                    ph.gapped_start_q,
                    ph.gapped_start_s,
                    adj_mat,
                    scaled_gap_open,
                    scaled_gap_extend,
                    scaled_x_drop_final,
                );
                let gr_opt = match (&bounded_gr_opt, bidir_gr) {
                    (Some(b), Some(d)) => {
                        if b.score > d.score {
                            Some(b.clone())
                        } else {
                            Some(d)
                        }
                    }
                    (Some(b), None) => Some(b.clone()),
                    (None, d) => d,
                };
                if let Some(gr) = gr_opt {
                    let q_slice = &query_aa[gr.query_start..gr.query_end];
                    let s_slice = &redo_subj_aa[gr.subject_start..gr.subject_end];
                    let (qs, ss) = gr.edit_script.render_alignment(
                        q_slice,
                        s_slice,
                        crate::protein::ncbistdaa_to_char,
                    );
                    new_phits.push(crate::protein_lookup::ProteinHit {
                        query_start: gr.query_start,
                        query_end: gr.query_end,
                        subject_start: gr.subject_start,
                        subject_end: gr.subject_end,
                        score: crate::math::nint(gr.score as f64 / comp_scale) as i32,
                        num_ident: gr.num_ident,
                        align_length: gr.align_length,
                        mismatches: gr.mismatches,
                        gap_opens: gr.gap_opens,
                        qseq: Some(qs),
                        sseq: Some(ss),
                        scaled_score: Some(gr.score),
                        gapped_start_q: ph.gapped_start_q,
                        gapped_start_s: ph.gapped_start_s,
                    });
                }
            }
        }
        (new_phits, true)
    } else {
        (phits, false)
    };

    (final_phits, use_adj_matrix, lambda_ratio_opt)
}

/// Port of NCBI `Blast_HSPListPurgeHSPsWithCommonEndpoints`
/// (`blast_hits.c:2455`). Two-pass dedup:
/// 1. Sort by (query.offset, subject.offset) ASC, score DESC, query.end DESC,
///    subject.end DESC. Within consecutive HSPs sharing same query.offset and
///    same subject.offset: keep the FIRST (highest score), drop the rest.
/// 2. Sort by (query.end, subject.end) ASC, score DESC, query.offset DESC,
///    subject.offset DESC. Same rule on (query.end, subject.end).
///
/// For protein BLAST the engine calls this with `purge=TRUE` at
/// `blast_engine.c:544` after `BLAST_GetGappedScore`.
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
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'N' | b'n' => b'N',
            b'R' => b'Y',
            b'Y' => b'R',
            b'M' => b'K',
            b'K' => b'M',
            b'W' => b'W',
            b'S' => b'S',
            b'B' => b'V',
            b'V' => b'B',
            b'D' => b'H',
            b'H' => b'D',
            _ => b'N',
        })
        .collect()
}

/// Six-frame translation of an ASCII nucleotide sequence (public API).
/// Internally calls the 1-1 port of NCBI `BLAST_GetAllTranslations`
/// (NCBI4na path) and unpacks the resulting flat translation buffer
/// into per-frame `TranslatedFrame` structs for caller convenience.
pub fn six_frame_translate(nt_seq: &[u8]) -> [TranslatedFrame; 6] {
    six_frame_translate_with_table(nt_seq, &crate::util::STANDARD_GENETIC_CODE)
}

fn six_frame_translate_with_table(
    nt_seq: &[u8],
    genetic_code: &'static [u8; 64],
) -> [TranslatedFrame; 6] {
    let nt_ncbi4na = ascii_to_ncbi4na(nt_seq);
    let (translation_buffer, frame_offsets) = crate::util::blast_get_all_translations_ncbi4na(
        &nt_ncbi4na,
        nt_ncbi4na.len(),
        genetic_code,
    );
    let mut frames: [TranslatedFrame; 6] = std::array::from_fn(|ctx| {
        let frame = crate::util::blast_context_to_frame_blastx(ctx as u32);
        let begin = (frame_offsets[ctx] + 1) as usize;
        let end = frame_offsets[ctx + 1] as usize;
        let ascii: Vec<u8> = if begin < end {
            translation_buffer[begin..end]
                .iter()
                .map(|&b| ncbistdaa_to_ascii(b))
                .collect()
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
    // by `blast_context_to_frame_blastx`. Keep the array as-is.
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
pub fn apply_seg_ncbistdaa(seq: &mut [u8]) {
    apply_seg_ncbistdaa_with_options(seq, 12, 2.2, 2.5)
}

pub fn apply_seg_ncbistdaa_with_options(seq: &mut [u8], window: usize, locut: f64, hicut: f64) {
    let masked = AMINOACID_TO_NCBISTDAA[b'X' as usize & 0x7F];
    if is_single_residue_low_complexity(seq) {
        for aa in seq {
            if *aa != masked {
                *aa = masked;
            }
        }
        return;
    }
    let mask = crate::filter::seg_filter_ncbistdaa(seq, window, locut, hicut);
    for r in &mask.regions {
        let start = r.start.max(0) as usize;
        let end = (r.end as usize).min(seq.len());
        for aa in &mut seq[start..end] {
            if *aa != masked {
                *aa = masked;
            }
        }
    }
}

fn is_single_residue_low_complexity(seq: &[u8]) -> bool {
    if seq.len() < 12 {
        return false;
    }
    let mut residue = None;
    let mut true_count = 0usize;
    for &aa in seq {
        let is_true = crate::composition::TRUE_CHAR_POSITIONS.contains(&(aa as usize)) || aa == 24;
        if !is_true {
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
pub fn lowercase_mask(seq: &[u8]) -> Vec<bool> {
    seq.iter().map(|b| b.is_ascii_lowercase()).collect()
}

/// Apply repeat masking based on n-mer frequency.
pub fn apply_repeat_mask(seq: &mut [u8]) {
    let mask = repeat_mask(seq, 11, 2.0);
    for (i, masked) in mask.iter().enumerate() {
        if *masked {
            seq[i] = b'N';
        }
    }
}

/// Compute repeat mask: positions where n-mer frequency exceeds threshold.
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
    fn default() -> Self {
        Self::new()
    }
}

impl BlastnSearch {
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

    pub fn query(mut self, seq: &[u8]) -> Self {
        self.query_raw = seq.to_vec();
        self
    }
    pub fn subject(mut self, seq: &[u8]) -> Self {
        self.subject_raw = seq.to_vec();
        self
    }
    pub fn word_size(mut self, ws: usize) -> Self {
        self.word_size = ws;
        self
    }
    pub fn reward(mut self, r: i32) -> Self {
        self.reward = r;
        self
    }
    pub fn penalty(mut self, p: i32) -> Self {
        self.penalty = p;
        self
    }
    pub fn gap_open(mut self, g: i32) -> Self {
        self.gap_open = g;
        self
    }
    pub fn gap_extend(mut self, g: i32) -> Self {
        self.gap_extend = g;
        self
    }
    pub fn evalue(mut self, e: f64) -> Self {
        self.evalue = e;
        self
    }
    pub fn dust(mut self, d: bool) -> Self {
        self.dust = d;
        self
    }
    pub fn strand(mut self, s: Strand) -> Self {
        self.strand = s;
        self
    }

    pub fn run(&self) -> Vec<SearchHsp> {
        if self.query_raw.is_empty() || self.subject_raw.is_empty() {
            return Vec::new();
        }

        let mut query_plus: Vec<u8> = self
            .query_raw
            .iter()
            .map(|&b| blastn_iupac_to_blastna(b))
            .collect();

        if self.dust {
            let mask = crate::filter::dust_filter(&query_plus, 20, 64, 1);
            mask.apply(&mut query_plus, 14);
        }

        let query_plus_nomask: Vec<u8> = self
            .query_raw
            .iter()
            .map(|&b| blastn_iupac_to_blastna(b))
            .collect();
        let query_minus: Vec<u8> = query_plus
            .iter()
            .rev()
            .map(|&b| blastn_complement(b))
            .collect();
        let query_minus_nomask: Vec<u8> = query_plus_nomask
            .iter()
            .rev()
            .map(|&b| blastn_complement(b))
            .collect();

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

        let subject: Vec<u8> = self
            .subject_raw
            .iter()
            .map(|&b| blastn_iupac_to_blastna(b))
            .collect();

        let m = build_blastna_matrix(self.reward, self.penalty);
        let matrix_fn = |i: usize, j: usize| -> i32 { m[i][j] };
        let mut lo = i32::MAX;
        let mut hi = i32::MIN;
        for i in 0..16 {
            for j in 0..16 {
                let s = m[i][j];
                if s <= -100000000 || s >= 100000000 {
                    continue;
                }
                if s < lo {
                    lo = s;
                }
                if s > hi {
                    hi = s;
                }
            }
        }

        let ambig: &[u8] = &[14, 15];
        let ctx = UngappedKbpContext {
            query_offset: 0,
            query_length: query_plus.len() as i32,
            is_valid: true,
        };
        let kbp_results = ungapped_kbp_calc(&query_plus, &[ctx], lo, hi, 16, ambig, &matrix_fn);
        let ungapped_kbp = kbp_results[0].clone().unwrap_or(KarlinBlk {
            lambda: 1.374,
            k: 0.621,
            log_k: 0.621_f64.ln(),
            h: 1.286,
            round_down: false,
        });
        let (gapped_kbp, _) = nucl_gapped_kbp_lookup(
            self.gap_open,
            self.gap_extend,
            self.reward,
            self.penalty,
            &ungapped_kbp,
        )
        .unwrap_or((ungapped_kbp.clone(), false));

        let (alpha, beta) = nucl_alpha_beta(
            self.reward,
            self.penalty,
            self.gap_open,
            self.gap_extend,
            ungapped_kbp.lambda,
            ungapped_kbp.h,
            true,
        );
        let db_length = self.subject_raw.len() as i64;
        let (len_adj, _) = compute_length_adjustment_exact(
            gapped_kbp.k,
            gapped_kbp.log_k,
            alpha / gapped_kbp.lambda,
            beta,
            self.query_raw.len() as i32,
            db_length,
            1,
        );
        let eff_db = (db_length - len_adj as i64).max(1);
        let search_space = eff_db as f64 * (self.query_raw.len() as i32 - len_adj).max(1) as f64;
        let x_dropoff =
            (self.xdrop_gap_final * crate::math::NCBIMATH_LN2 / gapped_kbp.lambda) as i32;

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

fn blastn_iupac_to_blastna(b: u8) -> u8 {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        b'R' | b'r' => 4,
        b'Y' | b'y' => 5,
        b'M' | b'm' => 6,
        b'K' | b'k' => 7,
        b'W' | b'w' => 8,
        b'S' | b's' => 9,
        b'B' | b'b' => 10,
        b'D' | b'd' => 11,
        b'H' | b'h' => 12,
        b'V' | b'v' => 13,
        b'N' | b'n' => 14,
        _ => 15,
    }
}

fn blastn_complement(b: u8) -> u8 {
    match b {
        0 => 3,
        1 => 2,
        2 => 1,
        3 => 0,
        4 => 5,
        5 => 4,
        6 => 7,
        7 => 6,
        8 => 8,
        9 => 9,
        10 => 13,
        11 => 12,
        12 => 11,
        13 => 10,
        14 => 14,
        _ => 15,
    }
}

// ── Internal helpers ────────────────────────────────────────────────────────

/// Encode an ASCII IUPAC nucleotide query as NCBI4na via the
/// `IUPACNA_TO_NCBI4NA` table from `blast_encoding.h:84` (ported in
/// `crate::encoding`). `U`/`u` is mapped to T (NCBI's blast input layer
/// treats RNA bases as DNA before translation); table positions that
/// hold 0 (non-IUPAC characters) are remapped to N (15) so that
/// `codon_to_aa` produces NCBIstdaa `X` rather than translating a junk
/// byte. This mirrors what NCBI's CSeq_data conversion does on the way
/// into `BLAST_GetTranslation`.
fn ascii_to_ncbi4na(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .map(|&b| match b {
            b'U' | b'u' => 8,
            _ if (b as usize) < crate::encoding::IUPACNA_TO_NCBI4NA.len() => {
                let v = crate::encoding::IUPACNA_TO_NCBI4NA[b as usize];
                if v == 0 {
                    15
                } else {
                    v
                }
            }
            _ => 15,
        })
        .collect()
}

fn ncbistdaa_to_ascii(b: u8) -> u8 {
    if (b as usize) < AA_SIZE {
        NCBISTDAA_TO_AMINOACID[b as usize] as u8
    } else {
        b'X'
    }
}

fn encode_protein_query_nomask(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .map(|&b| AMINOACID_TO_NCBISTDAA[b as usize & 0x7F])
        .collect()
}

fn encode_protein_query(
    sequence: &[u8],
    filter_low_complexity: bool,
    seg_window: usize,
    seg_locut: f64,
    seg_hicut: f64,
) -> Vec<u8> {
    let mut query_aa = encode_protein_query_nomask(sequence);
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
        let mut seq = encode_protein_query_nomask(b"AAAAAAAAAAAAAAAAAAAA");
        apply_seg_ncbistdaa(&mut seq);
        assert_eq!(crate::composition::read_composition(&seq, AA_SIZE).1, 0);
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
}

fn build_midline(qseq: &str, sseq: &str) -> Vec<u8> {
    qseq.bytes()
        .zip(sseq.bytes())
        .map(|(q, s)| if q == s { b'|' } else { b' ' })
        .collect()
}

// ── Aliases for API compatibility ───────────────────────────────────────────

/// Get the scoring matrix for a given MatrixType.
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

fn protein_gapped_params_for_matrix(
    matrix: MatrixType,
    gap_open: i32,
    gap_extend: i32,
) -> Option<GappedParams> {
    crate::stat::lookup_matrix_params(protein_matrix_name(matrix), gap_open, gap_extend)
}

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

fn protein_gumbel_for_matrix(
    matrix: MatrixType,
    gap_open: i32,
    gap_extend: i32,
    db_length: i64,
) -> Option<GumbelBlk> {
    crate::stat::matrix_gumbel_blk(protein_matrix_name(matrix), gap_open, gap_extend, db_length)
}

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

fn suggested_word_threshold(matrix: MatrixType, program: crate::program::ProgramType) -> f64 {
    let mut threshold = match matrix {
        MatrixType::Blosum45 => 14.0,
        MatrixType::Blosum62 => 11.0,
        MatrixType::Blosum80 => 12.0,
        MatrixType::Pam30 => 16.0,
        MatrixType::Pam70 => 14.0,
        MatrixType::Identity => 27.0,
        _ => 11.0,
    };
    if crate::program::subject_is_translated(program) {
        threshold += 2.0;
    } else if matches!(program, crate::program::BLASTX | crate::program::TBLASTX) {
        threshold += 1.0;
    }
    threshold
}

/// Get a genetic code translation table by NCBI code number.
/// Returns a 64-byte table mapping codons to NCBIstdaa amino acid codes.
/// Supports all standard NCBI genetic codes (1-33).
pub fn get_codon_table(code: u8) -> &'static [u8; 64] {
    crate::util::lookup_genetic_code(code)
}

// ── Masking wrappers ────────────────────────────────────────────────────────

/// Apply DUST low-complexity masking in place (replaces masked positions with N).
pub fn apply_dust(seq: &mut [u8]) {
    let mask = crate::filter::dust_filter(seq, 20, 64, 1);
    for r in &mask.regions {
        let start = r.start.max(0) as usize;
        let end = (r.end as usize).min(seq.len());
        for i in start..end {
            seq[i] = b'N';
        }
    }
}

/// Apply SEG low-complexity masking in place (replaces masked positions with X).
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
pub fn apply_lowercase_mask_nucleotide(seq: &mut [u8]) {
    for b in seq.iter_mut() {
        if b.is_ascii_lowercase() {
            *b = b'N';
        }
    }
}

/// Replace lowercase characters with X (protein masking).
pub fn apply_lowercase_mask_protein(seq: &mut [u8]) {
    for b in seq.iter_mut() {
        if b.is_ascii_lowercase() {
            *b = b'X';
        }
    }
}

// ── Composition statistics ──────────────────────────────────────────────────

/// Compute amino acid composition (frequencies) from NCBIstdaa-encoded sequence.
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
                for i in 1..23usize {
                    for j in 1..23usize {
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
        for i in 1..23usize {
            for j in 1..23usize {
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

fn compo_expected_score(q: &[f64; 28], r: &[f64; 28], matrix: &ScoringMatrix) -> f64 {
    let mut mu = 0.0f64;
    for i in 1..23usize {
        for j in 1..23usize {
            mu += q[i] * r[j] * matrix.score(i as u8, j as u8) as f64;
        }
    }
    mu
}

fn expected_score_with_bg(matrix: &ScoringMatrix) -> f64 {
    let bg = &crate::matrix::AA_FREQUENCIES;
    // AA_FREQUENCIES is [f64; 20] in ACDEFGHIKLMNPQRSTVWY order, need to map to NCBIstdaa
    let ncbi_order = [
        1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22,
    ]; // NCBIstdaa indices
    let mut freq28 = [0.0f64; 28];
    for (k, &idx) in ncbi_order.iter().enumerate() {
        if k < bg.len() {
            freq28[idx] = bg[k];
        }
    }
    compo_expected_score(&freq28, &freq28, matrix)
}

// ── Low-level PSSM functions ────────────────────────────────────────────────

/// Build a PSSM from search results (for PSI-BLAST iteration).
pub fn build_pssm(
    query: &[u8],
    results: &[SearchResult],
    inclusion_evalue: f64,
    matrix_type: MatrixType,
    _lambda: f64,
) -> crate::pssm::Pssm {
    let query_aa: Vec<u8> = query
        .iter()
        .map(|&b| AMINOACID_TO_NCBISTDAA[b as usize & 0x7F])
        .collect();
    let matrix = get_matrix(matrix_type);
    let mut pssm = crate::pssm::Pssm::from_sequence(&query_aa, matrix);

    // Collect aligned subject sequences from results that pass inclusion threshold
    let aligned: Vec<Vec<u8>> = results
        .iter()
        .flat_map(|r| r.hsps.iter())
        .filter(|hsp| hsp.evalue <= inclusion_evalue)
        .filter_map(|hsp| {
            if !hsp.subject_aln.is_empty() {
                Some(
                    hsp.subject_aln
                        .iter()
                        .map(|&b| AMINOACID_TO_NCBISTDAA[b as usize & 0x7F])
                        .collect(),
                )
            } else {
                None
            }
        })
        .collect();

    if !aligned.is_empty() {
        pssm.update_from_alignment(&aligned, &crate::matrix::AA_FREQUENCIES, 0.5);
    }
    pssm
}

/// Search a database using a PSSM instead of a substitution matrix.
pub fn search_with_pssm(
    db: &BlastDb,
    _query: &[u8],
    pssm: &crate::pssm::Pssm,
    params: &SearchParams,
) -> Vec<SearchResult> {
    let prot_kbp = protein_kbp_for_matrix(params.matrix, params.gap_open, params.gap_extend);

    let total_subj_len: usize = (0..db.num_oids)
        .map(|oid| db.get_seq_len(oid) as usize)
        .sum();
    let search_space = (pssm.length as f64) * (total_subj_len as f64);

    let subj_pairs: Vec<(String, Vec<u8>)> = (0..db.num_oids)
        .map(|oid| {
            let acc = db
                .get_accession(oid)
                .unwrap_or_else(|| format!("oid_{}", oid));
            let raw = db.get_sequence(oid);
            let aa: Vec<u8> = raw.iter().filter(|&&b| b != 0).copied().collect();
            (acc, aa)
        })
        .collect();

    let hits = crate::pssm::psi_blast_iteration(
        pssm,
        &subj_pairs,
        params.evalue_threshold,
        search_space,
        prot_kbp.lambda,
        prot_kbp.k,
    );

    hits.iter()
        .map(|h| SearchResult {
            subject_oid: subj_pairs
                .iter()
                .position(|(id, _)| id == &h.subject_id)
                .unwrap_or(0) as u32,
            subject_title: h.subject_id.clone(),
            subject_accession: h.subject_id.clone(),
            subject_len: h.subject_len,
            hsps: vec![Hsp {
                score: h.score,
                bit_score: prot_kbp.raw_to_bit(h.score),
                evalue: h.evalue,
                query_start: 0,
                query_end: h.align_len,
                subject_start: h.subject_start,
                subject_end: h.subject_start + h.align_len,
                num_identities: 0,
                num_gaps: 0,
                alignment_length: h.align_len,
                query_aln: Vec::new(),
                midline: Vec::new(),
                subject_aln: Vec::new(),
                query_frame: 0,
                subject_frame: 0,
            }],
            taxids: vec![],
        })
        .collect()
}

/// Alias for `blastn_search` matching the old API.
/// Note: at the crate root this is available as `blastn_search` since the `blastn` name
/// is used by the builder-pattern module. Use `blast_rs::api::blastn()` for the function form.
pub fn blastn(db: &BlastDb, query: &[u8], params: &SearchParams) -> Vec<SearchResult> {
    blastn_search(db, query, params)
}

/// Alias for `blastp` — low-level protein search backend.
pub fn blast_search(db: &BlastDb, query: &[u8], params: &SearchParams) -> Vec<SearchResult> {
    blastp(db, query, params)
}

/// Low-level blastx search (alias for `blastx`).
pub fn blastx_search(db: &BlastDb, query: &[u8], params: &SearchParams) -> Vec<SearchResult> {
    blastx(db, query, params)
}

/// Low-level tblastn search (alias for `tblastn`).
pub fn tblastn_search(db: &BlastDb, query: &[u8], params: &SearchParams) -> Vec<SearchResult> {
    tblastn(db, query, params)
}

/// Low-level tblastx search (alias for `tblastx`).
pub fn tblastx_search(db: &BlastDb, query: &[u8], params: &SearchParams) -> Vec<SearchResult> {
    tblastx(db, query, params)
}

/// Low-level PSI-BLAST iteration (alias for `psiblast`).
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
    /// E-value threshold for including hits in PSSM construction (default: 0.001).
    pub inclusion_evalue: f64,
}

impl PsiblastParams {
    pub fn new(search: SearchParams) -> Self {
        PsiblastParams {
            search,
            num_iterations: 3,
            inclusion_evalue: 0.001,
        }
    }
    pub fn num_iterations(mut self, v: u32) -> Self {
        self.num_iterations = v;
        self
    }
    pub fn inclusion_evalue(mut self, v: f64) -> Self {
        self.inclusion_evalue = v;
        self
    }
}

/// Run iterative PSI-BLAST search.
/// Returns final-round hits and the resulting PSSM.
pub fn psiblast(
    db: &BlastDb,
    query: &[u8],
    params: &PsiblastParams,
) -> (Vec<SearchResult>, crate::pssm::Pssm) {
    let query_aa: Vec<u8> = query
        .iter()
        .map(|&b| AMINOACID_TO_NCBISTDAA[b as usize & 0x7F])
        .collect();
    let matrix = *get_matrix(params.search.matrix);

    // Initial PSSM from query
    let mut pssm = crate::pssm::Pssm::from_sequence(&query_aa, &matrix);

    let prot_kbp = protein_kbp_for_matrix(
        params.search.matrix,
        params.search.gap_open,
        params.search.gap_extend,
    );

    let total_subj_len: usize = (0..db.num_oids)
        .map(|oid| db.get_seq_len(oid) as usize)
        .sum();
    let search_space = (query_aa.len() as f64) * (total_subj_len as f64);

    let mut final_results = Vec::new();

    // Build subject pairs for PSI-BLAST iteration
    let subj_pairs: Vec<(String, Vec<u8>)> = (0..db.num_oids)
        .map(|oid| {
            let acc = db
                .get_accession(oid)
                .unwrap_or_else(|| format!("oid_{}", oid));
            let raw = db.get_sequence(oid);
            let aa: Vec<u8> = raw.iter().filter(|&&b| b != 0).copied().collect();
            (acc, aa)
        })
        .collect();

    for _iter in 0..params.num_iterations {
        let hits = crate::pssm::psi_blast_iteration(
            &pssm,
            &subj_pairs,
            params.inclusion_evalue,
            search_space,
            prot_kbp.lambda,
            prot_kbp.k,
        );

        if hits.is_empty() {
            break;
        }

        final_results = hits
            .iter()
            .map(|h| SearchResult {
                subject_oid: subj_pairs
                    .iter()
                    .position(|(id, _)| id == &h.subject_id)
                    .unwrap_or(0) as u32,
                subject_title: h.subject_id.clone(),
                subject_accession: h.subject_id.clone(),
                subject_len: h.subject_len,
                hsps: vec![Hsp {
                    score: h.score,
                    bit_score: prot_kbp.raw_to_bit(h.score),
                    evalue: h.evalue,
                    query_start: 0,
                    query_end: h.align_len,
                    subject_start: h.subject_start,
                    subject_end: h.subject_start + h.align_len,
                    num_identities: 0,
                    num_gaps: 0,
                    alignment_length: h.align_len,
                    query_aln: Vec::new(),
                    midline: Vec::new(),
                    subject_aln: Vec::new(),
                    query_frame: 0,
                    subject_frame: 0,
                }],
                taxids: vec![],
            })
            .collect();

        // Update PSSM from aligned sequences
        let aligned: Vec<Vec<u8>> = hits
            .iter()
            .filter_map(|h| {
                subj_pairs
                    .iter()
                    .find(|(id, _)| id == &h.subject_id)
                    .and_then(|(_, seq)| {
                        if h.subject_start >= seq.len() {
                            return None;
                        }
                        let end = (h.subject_start + h.align_len).min(seq.len());
                        Some(seq[h.subject_start..end].to_vec())
                    })
            })
            .collect();

        if !aligned.is_empty() {
            pssm.update_from_alignment(&aligned, &crate::matrix::AA_FREQUENCIES, 0.5);
        }
    }

    (final_results, pssm)
}

#[cfg(test)]
mod tests {
    use super::*;

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
}
