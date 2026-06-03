//! Rust port of `hspfilter_culling.c` — winnowing-based HSP filter
//! described in Berman, Zhang, Wolf, Koonin & Miller (J Comput Biol
//! 2000;7(1-2):293-302). Used as a `BlastHSPWriter` in NCBI's pipeline
//! to drop HSPs dominated by overlapping higher-quality ones.
//!
//! This module ports the culling writer/pipe and its interval-tree
//! utilities 1-1 with their NCBI C counterparts. API and CLI culling
//! filters reuse the same dominance rule so rendered-result filtering
//! follows NCBI's score/length overlap semantics.
//!
//! Status (counted from the audit's missing-functions list):
//!
//! | Function | LOC | Ported? |
//! |---|---:|---|
//! | `s_HSPCopy` | 4 | ✅ |
//! | `s_HSPFree` | 4 | (Rust `Drop` handles it) |
//! | `s_DominateTest` | 25 | ✅ |
//! | `s_FullPass` | 8 | ✅ |
//! | `s_ProcessHSPList` | 20 | ✅ |
//! | `s_MarkDownHSPList` | 18 | ✅ |
//! | `s_AddHSPtoList` | 4 | ✅ |
//! | `s_GetNode` / `s_RetNode` | 5 | (Rust `Box` handles it) |
//! | `s_CTreeNodeNew` | 14 | ✅ |
//! | `s_CTreeNodeFree` | 5 | (Rust `Drop` handles it) |
//! | `s_ForkChildren` | 28 | ✅ |
//! | `s_MarkDownCTree` | 9 | ✅ |
//! | `s_ProcessCTree` | 23 | ✅ |
//! | `s_CTreeNew` / `s_CTreeFree` | 12 | ✅ / Rust `Drop` |
//! | `s_RipHSPOffCTree` | 16 | ✅ |
//! | `s_SaveHSP` | 22 | ✅ |
//! | `s_BlastHSPCullingInit` | 4 | ✅ |
//! | `s_BlastHSPCullingFinal` | 66 | ✅ |
//! | `s_BlastHSPCullingRun` | 26 | ✅ |
//! | `s_BlastHSPCullingFree` / `s_BlastHSPCullingNew` | 24 | ✅ / ✅ |
//! | `s_BlastHSPCullingPipeRun` / `Free` / `New` | 46 | ✅ |
//! | `BlastHSPCullingParamsNew` / `Free` | 21 | ✅ |
//! | `BlastHSPCullingInfoNew` / `BlastHSPCullingPipeInfoNew` | 15 | ✅ |

pub use crate::hspstream::Hsp;
use crate::program::ProgramType;

/// 1-1 port of `BlastHSPCullingData` (`hspfilter_culling.c:476`).
///
/// Working state for the culling writer/pipe driver. Owns the
/// per-context `CTreeNode` forest and a snapshot of the params /
/// query info needed across run callbacks. NCBI's C uses a `void *`
/// payload tagged via the `BlastHSPWriter` vtable; the Rust port
/// uses a concrete struct directly.
pub struct BlastHSPCullingData<'a> {
    pub params: BlastHSPCullingParams,
    pub query_info: &'a crate::queryinfo::QueryInfo,
    pub num_contexts: i32,
    pub c_tree: Vec<Option<Box<CTreeNode>>>,
}

impl<'a> BlastHSPCullingData<'a> {
    /// 1-1 port of `s_BlastHSPCullingNew` (`hspfilter_culling.c:666`).
    /// naming: Rust keeps HSPCulling as the established `hspculling` token.
    /// Allocates the writer-data struct and snapshots the
    /// `num_contexts = query_info.last_context + 1` derived value.
    pub fn s_blast_hspculling_new(
        params: BlastHSPCullingParams,
        query_info: &'a crate::queryinfo::QueryInfo,
    ) -> Self {
        let num_contexts = query_info.contexts.len() as i32;
        Self {
            params,
            query_info,
            num_contexts,
            c_tree: Vec::new(),
        }
    }

    /// 1-1 port of `s_BlastHSPCullingInit` (`hspfilter_culling.c:487`).
    /// naming: Rust keeps HSPCulling as the established `hspculling` token.
    /// Allocates the per-context tree forest. Called before
    /// `cull_run` when the pipeline starts.
    pub fn s_blast_hspculling_init(&mut self) {
        self.c_tree = (0..self.num_contexts).map(|_| None).collect();
    }

    /// 1-1 port of `s_BlastHSPCullingRun` (`hspfilter_culling.c:602`).
    /// naming: Rust keeps HSPCulling as the established `hspculling` token.
    ///
    /// Walks an `HspList` and inserts each HSP into the appropriate
    /// per-context culling tree via `s_save_hsp`. Handles blastn's
    /// strand-symmetric context normalization (NCBI's
    /// `cid = isBlastn ? (context - context % NUM_STRANDS) : context`).
    /// On insertion success the source slot is logically dropped
    /// (NCBI sets `hsp_array[i] = NULL`); we drain `hsp_list.hsps`
    /// regardless and free the (now-empty) list at the end —
    /// matching NCBI's `Blast_HSPListFree(hsp_list);` finalizer.
    pub fn s_blast_hspculling_run(&mut self, hsp_list: &mut crate::hspstream::HspList) {
        let is_blastn = self.params.program == crate::program::BLASTN;
        let oid = hsp_list.oid;

        // Drain the list (taking ownership of every HSP). NCBI mutates
        // hsp_array slot pointers in place and then frees the
        // container; we drain via `std::mem::take` and drop the empty
        // Vec at function end — semantically identical.
        let drained: Vec<Hsp> = std::mem::take(&mut hsp_list.hsps);
        for hsp in drained {
            let qlen = self
                .query_info
                .contexts
                .get(hsp.context as usize)
                .map(|c| c.query_length)
                .unwrap_or(0);
            let cid = if is_blastn {
                hsp.context - (hsp.context % crate::util::NUM_STRANDS as i32)
            } else {
                hsp.context
            };
            let (begin, end) = if is_blastn && (hsp.context % crate::util::NUM_STRANDS as i32) != 0
            {
                (qlen - hsp.query_end, qlen - hsp.query_offset)
            } else {
                (hsp.query_offset, hsp.query_end)
            };
            let mut node = LinkedHsp {
                hsp,
                cid: cid,
                sid: oid,
                begin,
                end,
                merit: self.params.culling_max,
                next: None,
            };
            // Lazy-init the per-context tree if needed.
            let cid_idx = cid as usize;
            if cid_idx >= self.c_tree.len() {
                self.c_tree.resize_with(cid_idx + 1, || None);
            }
            if self.c_tree[cid_idx].is_none() {
                self.c_tree[cid_idx] = Some(s_ctree_new(qlen));
            }
            // Save into the tree; NCBI ignores the false return from
            // s_SaveHSP at the writer level — the dropped HSP is just
            // discarded.
            if let Some(tree) = self.c_tree[cid_idx].as_mut() {
                let _ = s_save_hsp(tree, &mut node);
            }
        }
        // The original `hsps` Vec is already drained; let it drop with
        // hsp_list, mirroring NCBI's `Blast_HSPListFree(hsp_list)`.
    }

    /// 1-1 port of `s_BlastHSPCullingFinal` (`hspfilter_culling.c:500`).
    /// naming: Rust keeps HSPCulling as the established `hspculling` token.
    ///
    /// Rips every per-context tree into a flat list of HSPs, groups
    /// them by subject OID, and packages them into an `HspResults`
    /// shape. Each per-query hitlist gets `worst_evalue` and
    /// `low_score` set, and each per-subject `HspList` gets
    /// `best_evalue` populated and is sorted by score (NCBI calls
    /// `Blast_HSPListSortByScore`).
    pub fn s_blast_hspculling_final(&mut self) -> crate::hspstream::HspResults {
        let num_queries = self.query_info.num_queries.max(1);
        let mut results = crate::hspstream::HspResults::new(num_queries);

        for cid in 0..self.num_contexts {
            let tree_slot = self.c_tree.get_mut(cid as usize).and_then(|s| s.take());
            if tree_slot.is_none() {
                continue;
            }
            let qid =
                crate::queryinfo::blast_get_query_index_from_context(cid, self.params.program);
            // Lazy-init per-query hitlist.
            let qid_idx = qid as usize;
            if qid_idx >= results.hitlists.len() {
                results.hitlists.resize(qid_idx + 1, None);
            }
            if results.hitlists[qid_idx].is_none() {
                results.hitlists[qid_idx] = Some(crate::hspstream::HitList::new());
            }
            let hitlist = results.hitlists[qid_idx].as_mut().expect("just created");

            // Rip the tree into a flat linked list and consume it.
            let mut cull = s_rip_hsp_off_ctree(tree_slot);
            while let Some(node) = cull.take() {
                cull = node.next;
                let sid = node.sid;
                let hsp = node.hsp;
                // Look for an existing HspList by subject OID.
                let mut found = false;
                for list in hitlist.hsp_lists.iter_mut() {
                    if list.oid == sid {
                        let _ = crate::hspstream::blast_hsp_list_save_hsp(list, hsp.clone());
                        found = true;
                        break;
                    }
                }
                if !found {
                    let mut list = crate::hspstream::HspList::new(sid);
                    let _ = crate::hspstream::blast_hsp_list_save_hsp(&mut list, hsp);
                    hitlist.hsp_lists.push(list);
                }
            }

            // Sort each HspList and compute worst_evalue / low_score
            // for the hitlist (matching NCBI's tail of finalize).
            hitlist.low_score = i32::MAX;
            hitlist.worst_evalue = 0.0;
            for list in hitlist.hsp_lists.iter_mut() {
                // NCBI `s_BlastGetBestEvalue` (`blast_hits.c:~1300`) seeds
                // with `(double)INT4_MAX` and takes MIN over all evalues.
                // For an empty list NCBI returns INT4_MAX; we previously
                // returned INFINITY. The downstream comparator pushes
                // empty lists to the end anyway, but match NCBI exactly.
                let best = list
                    .hsps
                    .iter()
                    .map(|h| h.evalue)
                    .fold(i32::MAX as f64, f64::min);
                list.best_evalue = best;
                list.sort_by_score();
                if let Some(head) = list.hsps.first() {
                    hitlist.low_score = hitlist.low_score.min(head.score);
                    hitlist.worst_evalue = hitlist.worst_evalue.max(list.best_evalue);
                }
            }
        }
        results
    }
}

/// 1-1 port of `s_BlastHSPCullingPipeRun` (`hspfilter_culling.c:699`).
/// naming: Rust keeps the public pipe wrapper name already used by callers.
///
/// Pipe-stage variant of the writer: re-runs the culling filter
/// against an existing `HspResults`, updating it in place. NCBI:
/// 1. Sort each per-subject `HspList` by e-value, then sort the
///    per-query `HitList` by e-value.
/// 2. Run `s_BlastHSPCullingRun` over every (qid, sid) pair.
/// 3. Free the per-query `HitList` so finalize can repopulate.
/// 4. Run `s_BlastHSPCullingFinal` to write the new shape.
///
/// We reproduce the same control flow against our `HspResults`
/// shape. The intermediate teardown matches NCBI's
/// `Blast_HitListFree(results->hitlist_array[qid])` — Rust drops the
/// `HitList` when we replace it with `None`.
pub fn blast_hsp_culling_pipe_run(
    data: &mut BlastHSPCullingData<'_>,
    results: &mut crate::hspstream::HspResults,
) {
    // C: `s_BlastHSPCullingInit(data, results);` — fresh forest.
    data.s_blast_hspculling_init();

    // C step 1: sort each HspList by evalue and the HitList by
    // evalue.
    for hitlist_slot in results.hitlists.iter_mut() {
        if let Some(hitlist) = hitlist_slot.as_mut() {
            for list in hitlist.hsp_lists.iter_mut() {
                // C: `Blast_HSPListSortByEvalue(hsp_list);` then
                // `hsp_list->best_evalue = hsp_list->hsp_array[0]->evalue;`.
                // NCBI's `s_EvalueCompareHSPs` (`blast_hits.c:1415`) uses
                // an epsilon=1e-180 evalue threshold (`s_EvalueComp`) and
                // falls back to `ScoreCompareHSPs` for ties — a raw
                // `partial_cmp` would diverge for near-zero evalues and
                // for any equal-evalue pair (no tie-break).
                list.hsps.sort_by(crate::hspstream::evalue_compare_hsps);
                if let Some(head) = list.hsps.first() {
                    list.best_evalue = head.evalue;
                }
            }
            // C: `Blast_HitListSortByEvalue(...)`.
            hitlist
                .hsp_lists
                .sort_by(crate::hspstream::evalue_compare_hsp_lists);
        }
    }

    // C step 2 + 3: run culling on every per-subject list, then drop
    // the now-empty hitlist so finalize can rebuild it.
    for hitlist_slot in results.hitlists.iter_mut() {
        let Some(hitlist) = hitlist_slot.take() else {
            continue;
        };
        for mut list in hitlist.hsp_lists {
            data.s_blast_hspculling_run(&mut list);
        }
        // The `hitlist` value drops here, mirroring NCBI's
        // `Blast_HitListFree(...)`.
        let _ = hitlist_slot; // already taken; keep slot empty
    }

    // C step 4: rebuild via finalize. The previous pass cleared every
    // hitlist slot, so finalize starts from an empty results shape and
    // populates fresh hitlists for each query that has surviving HSPs.
    let new_results = data.s_blast_hspculling_final();
    *results = new_results;
}

/// 1-1 port of `s_BlastHSPCullingPipeNew` (`hspfilter_culling.c:752`).
/// naming: Rust keeps the public pipe constructor name already used by callers.
///
/// In NCBI's C this returns a `BlastHSPPipe*` populated with vtable
/// pointers (`RunFnPtr` = `s_BlastHSPCullingPipeRun`, `FreeFnPtr` =
/// `s_BlastHSPCullingPipeFree`) and a `BlastHSPCullingData` payload.
/// In Rust we just construct the payload directly — the vtable
/// indirection is unnecessary because the dispatch happens via Rust
/// method calls.
pub fn blast_hsp_culling_pipe_new<'a>(
    params: BlastHSPCullingParams,
    query_info: &'a crate::queryinfo::QueryInfo,
) -> BlastHSPCullingData<'a> {
    BlastHSPCullingData::s_blast_hspculling_new(params, query_info)
}

/// 1-1 port of `s_BlastHSPCullingFree` (`hspfilter_culling.c:683`).
/// naming: Rust keeps HSPCulling as the established `hspculling` token.
/// `Drop` handles recursive cleanup of the owned culling trees; this
/// direct wrapper preserves the C boundary where the data pointer is
/// nulled after release.
pub fn s_blast_hspculling_free(slot: &mut Option<BlastHSPCullingData<'_>>) {
    *slot = None;
}

/// 1-1 port of `s_BlastHSPCullingPipeFree` (`hspfilter_culling.c:736`).
/// naming: Rust keeps the public pipe free-hook name already used by callers.
/// `Drop` handles the actual deallocation; this hook is a parity
/// marker for callers wanting to match NCBI's flow.
pub fn blast_hsp_culling_pipe_free(slot: &mut Option<BlastHSPCullingData<'_>>) {
    s_blast_hspculling_free(slot);
}

/// 1-1 port of `BlastHSPCullingOptions`. NCBI's struct only carries
/// `max_hits` (max HSPs per query region). Other fields would extend
/// the configurability of the culling stage.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct BlastHSPCullingOptions {
    pub max_hits: i32,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BlastHSPBestHitOptions {
    pub overhang: f64,
    pub score_edge: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BlastHSPSubjectBestHitOptions {
    pub max_range_diff: i32,
}

#[derive(Debug, Clone, PartialEq)]
pub struct BlastHSPFilteringOptions {
    pub best_hit: Option<BlastHSPBestHitOptions>,
    pub best_hit_stage: crate::util::EBlastStage,
    pub subject_besthit_opts: Option<BlastHSPSubjectBestHitOptions>,
    pub culling_opts: Option<BlastHSPCullingOptions>,
    pub culling_stage: crate::util::EBlastStage,
}

const K_BEST_HIT_OVERHANG_MIN: f64 = 0.0;
const K_BEST_HIT_OVERHANG_MAX: f64 = 0.5;
const K_BEST_HIT_SCORE_EDGE_MIN: f64 = 0.0;
const K_BEST_HIT_SCORE_EDGE_MAX: f64 = 0.5;

/// Rust ownership equivalent of `BlastHSPBestHitOptionsNew`
/// (`blast_options.c:1816`).
/// naming: Rust keeps HSPBestHit as the established `hspbest_hit` token.
pub fn blast_hspbest_hit_options_new(overhang: f64, score_edge: f64) -> BlastHSPBestHitOptions {
    BlastHSPBestHitOptions {
        overhang,
        score_edge,
    }
}

/// Rust ownership equivalent of `BlastHSPBestHitOptionsFree`
/// (`blast_options.c:1848`).
/// naming: Rust keeps HSPBestHit as the established `hspbest_hit` token.
pub fn blast_hspbest_hit_options_free(
    opt: &mut Option<BlastHSPBestHitOptions>,
) -> Option<BlastHSPBestHitOptions> {
    *opt = None;
    None
}

/// Rust ownership equivalent of `BlastHSPCullingOptionsNew`
/// (`blast_options.c:1857`).
/// naming: Rust keeps HSPCulling as the established `hspculling` token.
pub fn blast_hspculling_options_new(max: i32) -> BlastHSPCullingOptions {
    BlastHSPCullingOptions { max_hits: max }
}

/// Rust ownership equivalent of `BlastHSPCullingOptionsFree`
/// (`blast_options.c:1880`).
/// naming: Rust keeps HSPCulling as the established `hspculling` token.
pub fn blast_hspculling_options_free(
    culling_opts: &mut Option<BlastHSPCullingOptions>,
) -> Option<BlastHSPCullingOptions> {
    *culling_opts = None;
    None
}

/// Rust ownership equivalent of `BlastHSPFilteringOptionsNew`
/// (`blast_options.c:1890`).
/// naming: Rust keeps HSPFiltering as the established `hspfiltering` token.
pub fn blast_hspfiltering_options_new() -> BlastHSPFilteringOptions {
    BlastHSPFilteringOptions {
        best_hit: None,
        best_hit_stage: crate::util::EBlastStage::None,
        subject_besthit_opts: None,
        culling_opts: None,
        culling_stage: crate::util::EBlastStage::None,
    }
}

/// Rust ownership equivalent of `BlastHSPFilteringOptions_AddBestHit`
/// (`blast_options.c:1897`). On success, ownership moves from `best_hit`
/// into `filt_opts`, matching the C pointer-to-pointer transfer.
/// naming: Rust keeps HSPFiltering as the established `hspfiltering` token.
pub fn blast_hspfiltering_options_add_best_hit(
    filt_opts: Option<&mut BlastHSPFilteringOptions>,
    best_hit: &mut Option<BlastHSPBestHitOptions>,
    stage: crate::util::EBlastStage,
) -> i16 {
    let Some(filt_opts) = filt_opts else {
        return 1;
    };
    let Some(best_hit) = best_hit.take() else {
        return 1;
    };

    filt_opts.best_hit = Some(best_hit);
    filt_opts.best_hit_stage = stage;
    0
}

/// Rust ownership equivalent of `BlastHSPFilteringOptions_AddCulling`
/// (`blast_options.c:1913`). On success, ownership moves from `culling`
/// into `filt_opts`, matching the C pointer-to-pointer transfer.
/// naming: Rust keeps HSPFiltering as the established `hspfiltering` token.
pub fn blast_hspfiltering_options_add_culling(
    filt_opts: Option<&mut BlastHSPFilteringOptions>,
    culling: &mut Option<BlastHSPCullingOptions>,
    stage: crate::util::EBlastStage,
) -> i16 {
    let Some(filt_opts) = filt_opts else {
        return 1;
    };
    let Some(culling) = culling.take() else {
        return 1;
    };

    filt_opts.culling_opts = Some(culling);
    filt_opts.culling_stage = stage;
    0
}

/// Rust ownership equivalent of `BlastHSPSubjectBestHitOptionsNew`
/// (`blast_options.c:1963`).
pub fn blast_hsp_subject_best_hit_options_new(_is_protein: bool) -> BlastHSPSubjectBestHitOptions {
    // NCBI defaults: 3 for both protein and nucleotide programs.
    BlastHSPSubjectBestHitOptions { max_range_diff: 3 }
}

/// Rust ownership equivalent of `BlastHSPSubjectBestHitOptionsFree`
/// (`blast_options.c:1987`).
pub fn blast_hsp_subject_best_hit_options_free(
    subject_besthit_opts: &mut Option<BlastHSPSubjectBestHitOptions>,
) -> Option<BlastHSPSubjectBestHitOptions> {
    *subject_besthit_opts = None;
    None
}

/// Rust ownership equivalent of `BlastHSPFilteringOptions_AddSubjectBestHit`
/// (`blast_options.c:1998`).
/// naming: Rust keeps HSPFiltering as the established `hspfiltering` token.
pub fn blast_hspfiltering_options_add_subject_best_hit(
    filt_opts: Option<&mut BlastHSPFilteringOptions>,
    subject_besthit: &mut Option<BlastHSPSubjectBestHitOptions>,
) -> i16 {
    let Some(filt_opts) = filt_opts else {
        return 1;
    };
    let Some(subject_besthit) = subject_besthit.take() else {
        return 1;
    };

    filt_opts.subject_besthit_opts = Some(subject_besthit);
    0
}

/// blast-rs: Stage-bit helper for filtering option validation; not a direct
/// NCBI C port.
fn blast_stage_has_prelim(stage: crate::util::EBlastStage) -> bool {
    (stage as i32 & crate::util::EBlastStage::PrelimSearch as i32) != 0
}

/// NCBI: BlastHSPBestHitOptionsValidate (blast_options.c).
/// naming: Rust keeps HSPBestHit as the established `hspbest_hit` token.
pub fn blast_hspbest_hit_options_validate(opts: &BlastHSPFilteringOptions) -> i16 {
    let Some(best_hit) = opts.best_hit else {
        return 0;
    };

    if best_hit.overhang <= K_BEST_HIT_OVERHANG_MIN || best_hit.overhang >= K_BEST_HIT_OVERHANG_MAX
    {
        return -1;
    }
    if best_hit.score_edge <= K_BEST_HIT_SCORE_EDGE_MIN
        || best_hit.score_edge >= K_BEST_HIT_SCORE_EDGE_MAX
    {
        return -1;
    }
    0
}

/// NCBI: BlastHSPCullingOptionsValidate (blast_options.c).
/// naming: Rust keeps HSPCulling as the established `hspculling` token.
pub fn blast_hspculling_options_validate(opts: &BlastHSPFilteringOptions) -> i16 {
    let Some(culling_opts) = opts.culling_opts else {
        return 0;
    };

    if culling_opts.max_hits < 0 {
        return -1;
    }
    0
}

/// NCBI: BlastHSPFilteringOptionsValidate (blast_options.c).
/// naming: Rust keeps HSPFiltering as the established `hspfiltering` token.
pub fn blast_hspfiltering_options_validate(opts: &BlastHSPFilteringOptions) -> i16 {
    let status = blast_hspbest_hit_options_validate(opts);
    if status != 0 {
        return status;
    }
    let writer_found = blast_stage_has_prelim(opts.best_hit_stage);

    let status = blast_hspculling_options_validate(opts);
    if status != 0 {
        return status;
    }
    if blast_stage_has_prelim(opts.culling_stage) && writer_found {
        return 1;
    }
    0
}

/// Rust ownership equivalent of `BlastHSPFilteringOptionsFree`
/// (`blast_options.c:1952`).
/// naming: Rust keeps HSPFiltering as the established `hspfiltering` token.
pub fn blast_hspfiltering_options_free(
    opts: &mut Option<BlastHSPFilteringOptions>,
) -> Option<BlastHSPFilteringOptions> {
    if let Some(opts) = opts {
        blast_hspbest_hit_options_free(&mut opts.best_hit);
        blast_hsp_subject_best_hit_options_free(&mut opts.subject_besthit_opts);
        blast_hspculling_options_free(&mut opts.culling_opts);
    }
    *opts = None;
    None
}

/// 1-1 port of `BlastHSPCullingParams` (`hspfilter_culling.h:57`).
/// Built by `BlastHSPCullingParamsNew` from
/// `(BlastHitSavingOptions, BlastHSPCullingOptions)` plus a few
/// program-derived flags.
#[derive(Debug, Clone)]
pub struct BlastHSPCullingParams {
    /// Program type (NCBI: `program`).
    pub program: ProgramType,
    /// Number of hits saved during the preliminary part of search
    /// (NCBI: `prelim_hitlist_size`).
    pub prelim_hitlist_size: i32,
    /// Number of HSPs to save per database sequence (NCBI:
    /// `hsp_num_max`).
    pub hsp_num_max: i32,
    /// Number of HSPs allowed per query region (NCBI: `culling_max`,
    /// initialized from `culling_opts.max_hits`).
    pub culling_max: i32,
}

/// 1-1 port of `BlastHSPCullingParamsNew` (`hspfilter_culling.c:783`).
///
/// NCBI builds an intermediate `BlastHSPCollectorParams` to derive
/// `prelim_hitlist_size` and `hsp_num_max` from the hit-saving
/// options; we pass `hsp_num_max` explicitly because Rust
/// `HitSavingOptions` does not carry that field.
pub fn blast_hsp_culling_params_new(
    hit_options: &crate::options::HitSavingOptions,
    culling_opts: &BlastHSPCullingOptions,
    hsp_num_max: i32,
    program: ProgramType,
    composition_based_stats: i32,
    gapped_calculation: bool,
) -> BlastHSPCullingParams {
    let prelim_hitlist_size = crate::hspstream::get_prelim_hitlist_size(
        hit_options.hitlist_size,
        composition_based_stats,
        gapped_calculation,
    );
    BlastHSPCullingParams {
        program,
        prelim_hitlist_size,
        hsp_num_max: crate::hspstream::blast_hsp_num_max(gapped_calculation, hsp_num_max),
        culling_max: culling_opts.max_hits,
    }
}

/// 1-1 port of `BlastHSPCullingParamsFree` (`hspfilter_culling.c:804`).
/// `Drop` does the work in Rust; this is a parity hook.
pub fn blast_hsp_culling_params_free(slot: &mut Option<BlastHSPCullingParams>) {
    *slot = None;
}

/// 1-1 port of `BlastHSPCullingInfoNew` (`hspfilter_culling.c:813`).
///
/// NCBI bundles the params into a `BlastHSPWriterInfo` whose
/// `NewFnPtr` points to `s_BlastHSPCullingNew`. The Rust port
/// represents the writer info as a small Rust struct that records the
/// params and a marker indicating which constructor should run later. The
/// `NewFnPtr` indirection isn't needed in Rust (we don't have a vtable to
/// populate); the constructor is invoked directly when the caller is ready.
#[derive(Debug, Clone)]
pub struct BlastHSPWriterInfo {
    pub params: BlastHSPCullingParams,
}

/// NCBI: BlastHSPCullingInfoNew (hspfilter_culling.c).
pub fn blast_hsp_culling_info_new(params: BlastHSPCullingParams) -> BlastHSPWriterInfo {
    BlastHSPWriterInfo { params }
}

/// Port of NCBI `BlastHSPBestHitParams` (`hspfilter_besthit.h:38`).
#[derive(Debug, Clone)]
pub struct BlastHSPBestHitParams {
    pub prelim_hitlist_size: i32,
    pub hsp_num_max: i32,
    pub program: ProgramType,
    pub overhang: f64,
    pub score_edge: f64,
}

/// Port of NCBI `LinkedHSP_BH` (`hspfilter_besthit.c:43`).
#[derive(Debug, Clone)]
pub struct LinkedHspBestHit {
    pub hsp: Hsp,
    pub sid: i32,
    pub begin: i32,
    pub end: i32,
    pub len: i32,
}

/// Port of NCBI `BlastHSPBestHitData` (`hspfilter_besthit.c:52`).
pub struct BlastHSPBestHitData<'a> {
    pub params: BlastHSPBestHitParams,
    pub query_info: &'a crate::queryinfo::QueryInfo,
    pub best_list: Vec<Vec<LinkedHspBestHit>>,
    pub num_hsps: Vec<i32>,
    pub max_hsps: Vec<i32>,
}

#[derive(Debug, Clone)]
pub struct BlastHSPBestHitInfo {
    pub params: BlastHSPBestHitParams,
}

#[derive(Debug, Clone)]
pub struct BlastHSPBestHitPipeInfo {
    pub params: BlastHSPBestHitParams,
    pub next: Option<Box<BlastHSPBestHitPipeInfo>>,
}

/// Port of NCBI `BlastHSPBestHitParamsNew` (`hspfilter_besthit.c:592`).
pub fn blast_hsp_best_hit_params_new(
    hit_options: &crate::options::HitSavingOptions,
    best_hit_opts: &BlastHSPBestHitOptions,
    composition_based_stats: i32,
    gapped_calculation: bool,
) -> BlastHSPBestHitParams {
    BlastHSPBestHitParams {
        prelim_hitlist_size: crate::hspstream::get_prelim_hitlist_size(
            hit_options.hitlist_size,
            composition_based_stats,
            gapped_calculation,
        ),
        hsp_num_max: crate::hspstream::blast_hsp_num_max(
            gapped_calculation,
            hit_options.hsp_num_max,
        ),
        program: hit_options.program_number,
        overhang: best_hit_opts.overhang,
        score_edge: best_hit_opts.score_edge,
    }
}

/// Rust ownership equivalent of `BlastHSPBestHitParamsFree`.
pub fn blast_hsp_best_hit_params_free(
    slot: &mut Option<BlastHSPBestHitParams>,
) -> Option<BlastHSPBestHitParams> {
    *slot = None;
    None
}

/// Port of NCBI `BlastHSPBestHitInfoNew` (`hspfilter_besthit.c:619`).
pub fn blast_hsp_best_hit_info_new(params: BlastHSPBestHitParams) -> BlastHSPBestHitInfo {
    BlastHSPBestHitInfo { params }
}

/// Port of NCBI `BlastHSPBestHitPipeInfoNew` (`hspfilter_besthit.c:627`).
pub fn blast_hsp_best_hit_pipe_info_new(params: BlastHSPBestHitParams) -> BlastHSPBestHitPipeInfo {
    BlastHSPBestHitPipeInfo { params, next: None }
}

/// blast-rs: Best-hit query-span adapter; not a direct NCBI C port.
fn best_hit_query_begin_len(
    hsp: &Hsp,
    query_info: &crate::queryinfo::QueryInfo,
    program: ProgramType,
    rps: bool,
    query_index: i32,
) -> Option<(i32, i32, i32)> {
    let qid = if rps {
        query_index
    } else {
        crate::queryinfo::blast_get_query_index_from_context(hsp.context, program)
    };
    let context = query_info.contexts.get(hsp.context as usize);
    let qlen = context.map(|ctx| ctx.query_length).unwrap_or_else(|| {
        query_info
            .contexts
            .iter()
            .find(|ctx| ctx.query_index == qid)
            .map(|ctx| ctx.query_length)
            .unwrap_or(0)
    });
    let frame = context.map(|ctx| ctx.frame).unwrap_or(0);
    let len = hsp.query_end - hsp.query_offset;
    if len <= 0 || qid < 0 {
        return None;
    }
    let begin = if !rps && frame < 0 {
        qlen - hsp.query_end
    } else {
        hsp.query_offset
    };
    Some((qid, begin, len))
}

/// NCBI: s_BlastHSPBestHitInit (hspfilter_besthit.c).
fn s_blast_hsp_best_hit_init(
    data: &mut BlastHSPBestHitData<'_>,
    results: &crate::hspstream::HspResults,
) -> i32 {
    let num_queries = results.hitlists.len();
    data.best_list = (0..num_queries).map(|_| Vec::new()).collect();
    data.num_hsps = vec![0; num_queries];
    data.max_hsps = vec![data.params.prelim_hitlist_size * 2; num_queries];
    0
}

/// NCBI: s_ExportToHitlist (hspfilter_besthit.c).
fn s_export_to_hitlist(
    qid: usize,
    data: &mut BlastHSPBestHitData<'_>,
    hit_list: &mut crate::hspstream::HitList,
) -> i32 {
    if qid >= data.best_list.len() {
        return -1;
    }
    let mut tmp_hit_list = crate::hspstream::blast_hit_list_new(data.num_hsps[qid]);
    let nodes = std::mem::take(&mut data.best_list[qid]);
    data.num_hsps[qid] = 0;
    for node in nodes {
        if let Some(list) = tmp_hit_list
            .hsp_lists
            .iter_mut()
            .find(|list| list.oid == node.sid)
        {
            let _ = crate::hspstream::blast_hsp_list_save_hsp(list, node.hsp);
        } else {
            let mut list = crate::hspstream::blast_hsp_list_new(data.params.hsp_num_max);
            list.oid = node.sid;
            let _ = crate::hspstream::blast_hsp_list_save_hsp(&mut list, node.hsp);
            tmp_hit_list.hsp_lists.push(list);
        }
    }
    for list in tmp_hit_list.hsp_lists.drain(..) {
        let _ = hit_list.blast_hit_list_update(list);
    }
    0
}

/// NCBI: s_ImportFromHitlist (hspfilter_besthit.c).
fn s_import_from_hitlist(
    qid: usize,
    data: &mut BlastHSPBestHitData<'_>,
    hit_list: &mut crate::hspstream::HitList,
) -> i32 {
    if qid >= data.best_list.len() {
        return -1;
    }
    let mut imported = Vec::new();
    for list in hit_list.hsp_lists.drain(..) {
        for hsp in list.hsps {
            let Some((_, begin, len)) = best_hit_query_begin_len(
                &hsp,
                data.query_info,
                data.params.program,
                false,
                qid as i32,
            ) else {
                continue;
            };
            let node = LinkedHspBestHit {
                hsp,
                sid: list.oid,
                begin,
                end: begin + len,
                len,
            };
            let pos = imported
                .iter()
                .position(|existing: &LinkedHspBestHit| existing.begin >= begin)
                .unwrap_or(imported.len());
            imported.insert(pos, node);
        }
    }
    data.num_hsps[qid] = imported.len() as i32;
    data.max_hsps[qid] = data.num_hsps[qid] * 2;
    data.best_list[qid] = imported;
    0
}

/// blast-rs: Best-hit insertion helper over owned Rust HSP nodes; not a direct
/// NCBI C port.
fn best_hit_insert(
    data: &mut BlastHSPBestHitData<'_>,
    sid: i32,
    query_index: i32,
    rps: bool,
    mut hsp: Hsp,
) {
    let Some((qid, mut begin, len_a)) =
        best_hit_query_begin_len(&hsp, data.query_info, data.params.program, rps, query_index)
    else {
        return;
    };
    let qid = qid as usize;
    if qid >= data.best_list.len() {
        return;
    }
    let mut end = begin + len_a;
    let score_a = hsp.score;
    let evalue_a = hsp.evalue;
    let param_overhang = data.params.overhang;
    let param_s = 1.0 - data.params.score_edge;
    let den_a = score_a as f64 / len_a as f64 / param_s;

    let list = &mut data.best_list[qid];
    for node in list.iter().filter(|node| {
        node.end >= end
            && if rps {
                node.begin < begin
            } else {
                node.begin <= begin
            }
    }) {
        let score_b = node.hsp.score;
        if node.end >= end
            && node.hsp.evalue <= evalue_a
            && score_b as f64 / node.len as f64 > den_a
        {
            return;
        }
    }

    let allowed_overhang =
        (2.0 * len_a as f64 * param_overhang / (1.0 - 2.0 * param_overhang)) as i32;
    let allowed_begin = begin - allowed_overhang;
    let allowed_end = end + allowed_overhang;
    let overhang = (len_a as f64 * param_overhang) as i32;
    begin -= overhang;
    end += overhang;
    let den_a = score_a as f64 / len_a as f64 * param_s;

    let mut index = 0usize;
    while index < list.len() {
        if list[index].begin < allowed_begin || list[index].begin >= allowed_end {
            index += 1;
            continue;
        }
        let len_b = list[index].len;
        let score_b = list[index].hsp.score;
        let old_overhang = (list[index].end - list[index].begin - len_b) / 2;
        if list[index].begin + old_overhang >= begin
            && list[index].end - old_overhang <= end
            && list[index].hsp.evalue >= evalue_a
            && score_b as f64 / (len_b as f64) < den_a
        {
            list.remove(index);
            data.num_hsps[qid] -= 1;
        } else {
            index += 1;
        }
    }

    if rps {
        hsp.context = query_index;
    }
    let node = LinkedHspBestHit {
        hsp,
        sid,
        begin,
        end,
        len: len_a,
    };
    let insert_at = list
        .iter()
        .position(|existing| existing.begin >= begin)
        .unwrap_or(list.len());
    list.insert(insert_at, node);
    data.num_hsps[qid] += 1;

    if data.num_hsps[qid] > data.max_hsps[qid] {
        let mut hitlist = crate::hspstream::blast_hit_list_new(data.num_hsps[qid]);
        let _ = s_export_to_hitlist(qid, data, &mut hitlist);
        let _ = s_import_from_hitlist(qid, data, &mut hitlist);
    }
}

/// Port of NCBI `s_BlastHSPBestHitFinal` (`hspfilter_besthit.c:195`).
pub fn s_blast_hsp_best_hit_final(
    data: &mut BlastHSPBestHitData<'_>,
    results: &mut crate::hspstream::HspResults,
) -> i32 {
    for qid in 0..results.hitlists.len() {
        if data
            .best_list
            .get(qid)
            .map(|list| list.is_empty())
            .unwrap_or(true)
        {
            continue;
        }
        let mut hitlist = crate::hspstream::blast_hit_list_new(data.num_hsps[qid]);
        let _ = s_export_to_hitlist(qid, data, &mut hitlist);
        for list in hitlist.hsp_lists.iter_mut() {
            list.hsps.sort_by(crate::hspstream::score_compare_hsps);
        }
        let _ = crate::hspstream::blast_hit_list_sort_by_evalue(&mut hitlist);
        let target = results.hitlists[qid].get_or_insert_with(|| {
            crate::hspstream::blast_hit_list_new(data.params.prelim_hitlist_size)
        });
        for list in hitlist.hsp_lists.drain(..) {
            let _ = target.blast_hit_list_update(list);
        }
    }
    data.best_list.clear();
    data.num_hsps.clear();
    data.max_hsps.clear();
    0
}

/// Port of NCBI `s_BlastHSPBestHitRun` (`hspfilter_besthit.c:238`).
pub fn s_blast_hsp_best_hit_run(
    data: &mut BlastHSPBestHitData<'_>,
    hsp_list: Option<crate::hspstream::HspList>,
) -> i32 {
    let Some(mut hsp_list) = hsp_list else {
        return 0;
    };
    if data.best_list.is_empty() {
        return 0;
    }
    let param_overhang = data.params.overhang;
    let param_s = 1.0 - data.params.score_edge;
    if param_overhang <= 0.0 || param_s <= 0.0 {
        hsp_list.hsps.clear();
        return 0;
    }
    let sid = hsp_list.oid;
    for hsp in hsp_list.hsps.drain(..) {
        let Some((qid, _begin, len)) =
            best_hit_query_begin_len(&hsp, data.query_info, data.params.program, false, 0)
        else {
            continue;
        };
        if qid < 0 || qid as usize >= data.best_list.len() || len <= 0 {
            continue;
        }
        best_hit_insert(data, sid, 0, false, hsp);
    }
    0
}

/// Port of NCBI `s_BlastHSPBestHitRun_RPS` (`hspfilter_besthit.c:359`).
pub fn s_blast_hsp_best_hit_run_rps(
    data: &mut BlastHSPBestHitData<'_>,
    query_index: i32,
    hsp_list: Option<crate::hspstream::HspList>,
) -> i32 {
    let Some(mut hsp_list) = hsp_list else {
        return 0;
    };
    if data.best_list.is_empty() {
        return 0;
    }
    let param_overhang = data.params.overhang;
    let param_s = 1.0 - data.params.score_edge;
    if param_overhang <= 0.0 || param_s <= 0.0 {
        hsp_list.hsps.clear();
        return 0;
    }
    if query_index < 0 || query_index as usize >= data.best_list.len() {
        hsp_list.hsps.clear();
        return 0;
    }
    for hsp in hsp_list.hsps.drain(..) {
        if hsp.query_end <= hsp.query_offset {
            continue;
        }
        best_hit_insert(data, hsp.context, query_index, true, hsp);
    }
    0
}

/// Port of NCBI `s_BlastHSPBestHitPipeRun` (`hspfilter_besthit.c:520`).
pub fn s_blast_hsp_best_hit_pipe_run(
    data: &mut BlastHSPBestHitData<'_>,
    results: &mut crate::hspstream::HspResults,
) -> i32 {
    let _ = s_blast_hsp_best_hit_init(data, results);
    let _ = crate::hspstream::blast_hsp_results_sort_by_evalue(results);
    for qid in 0..results.hitlists.len() {
        let Some(mut hitlist) = results.hitlists[qid].take() else {
            continue;
        };
        for list in hitlist.hsp_lists.drain(..) {
            let _ = s_blast_hsp_best_hit_run(data, Some(list));
        }
    }
    s_blast_hsp_best_hit_final(data, results)
}

/// Rust ownership equivalent of `s_BlastHSPBestHitPipeFree`.
pub fn s_blast_hsp_best_hit_pipe_free<'a>(
    slot: &mut Option<BlastHSPBestHitData<'a>>,
) -> Option<BlastHSPBestHitData<'a>> {
    *slot = None;
    None
}

/// Port of NCBI `s_BlastHSPBestHitPipeNew` (`hspfilter_besthit.c:560`).
pub fn s_blast_hsp_best_hit_pipe_new<'a>(
    params: BlastHSPBestHitParams,
    query_info: &'a crate::queryinfo::QueryInfo,
) -> BlastHSPBestHitData<'a> {
    BlastHSPBestHitData {
        params,
        query_info,
        best_list: Vec::new(),
        num_hsps: Vec::new(),
        max_hsps: Vec::new(),
    }
}

#[derive(Debug, Clone)]
pub struct BlastHSPCollectorParams {
    pub program: ProgramType,
    pub prelim_hitlist_size: i32,
    pub hsp_num_max: i32,
}

/// Port of NCBI `BlastHSPCollectorParamsNew` (`hspfilter_collector.c:325`).
pub fn blast_hsp_collector_params_new(
    hit_options: &crate::options::HitSavingOptions,
    program: ProgramType,
    composition_based_stats: i32,
    gapped_calculation: bool,
    hsp_num_max: i32,
) -> BlastHSPCollectorParams {
    let prelim_hitlist_size = crate::hspstream::get_prelim_hitlist_size(
        hit_options.hitlist_size,
        composition_based_stats,
        gapped_calculation,
    );
    BlastHSPCollectorParams {
        program,
        prelim_hitlist_size,
        hsp_num_max: crate::hspstream::blast_hsp_num_max(gapped_calculation, hsp_num_max),
    }
}

/// Rust ownership equivalent of `BlastHSPCollectorParamsFree`.
pub fn blast_hsp_collector_params_free(
    slot: &mut Option<BlastHSPCollectorParams>,
) -> Option<BlastHSPCollectorParams> {
    *slot = None;
    None
}

#[derive(Debug, Clone)]
pub struct BlastHSPCollectorData {
    pub params: BlastHSPCollectorParams,
}

#[derive(Debug, Clone)]
pub struct BlastHSPCollectorWriter {
    pub data: BlastHSPCollectorData,
    pub rps_run: bool,
}

#[derive(Debug, Clone)]
pub struct BlastHSPCollectorInfo {
    pub params: BlastHSPCollectorParams,
}

/// Port of NCBI `s_BlastHSPCollectorFinal` (`hspfilter_collector.c:68`).
/// naming: Rust keeps HSPCollector as the established `hspcollector` token.
pub fn s_blast_hspcollector_final(
    _data: &mut BlastHSPCollectorData,
    results: &mut crate::hspstream::HspResults,
) -> i32 {
    crate::hspstream::blast_hsp_results_sort_by_evalue(results)
}

/// Port of NCBI `s_BlastHSPCollectorRun` (`hspfilter_collector.c:82`).
/// naming: Rust keeps HSPCollector as the established `hspcollector` token.
pub fn s_blast_hspcollector_run(
    data: &mut BlastHSPCollectorData,
    results: &mut crate::hspstream::HspResults,
    hsp_list: Option<crate::hspstream::HspList>,
) -> i32 {
    let Some(mut hsp_list) = hsp_list else {
        return 0;
    };
    if results.hitlists.is_empty() {
        return -1;
    }

    if results.hitlists.len() > 1 {
        let mut split_lists: Vec<Option<crate::hspstream::HspList>> =
            (0..results.hitlists.len()).map(|_| None).collect();

        for hsp in hsp_list.hsps.drain(..) {
            let query_index = crate::queryinfo::blast_get_query_index_from_context(
                hsp.context,
                data.params.program,
            ) as usize;
            if query_index >= results.hitlists.len() {
                return -1;
            }
            let split_list = split_lists[query_index].get_or_insert_with(|| {
                let mut list = crate::hspstream::blast_hsp_list_new(data.params.hsp_num_max);
                list.oid = hsp_list.oid;
                list
            });
            let _ = crate::hspstream::blast_hsp_list_save_hsp(split_list, hsp);
        }

        for (query_index, split_list) in split_lists.into_iter().enumerate() {
            let Some(mut split_list) = split_list else {
                continue;
            };
            split_list
                .hsps
                .sort_by(crate::hspstream::score_compare_hsps);
            let hitlist = results.hitlists[query_index].get_or_insert_with(|| {
                crate::hspstream::blast_hit_list_new(data.params.prelim_hitlist_size)
            });
            hitlist.blast_hit_list_update(split_list);
        }
    } else if !hsp_list.hsps.is_empty() {
        let hitlist = results.hitlists[0].get_or_insert_with(|| {
            crate::hspstream::blast_hit_list_new(data.params.prelim_hitlist_size)
        });
        hsp_list.hsp_max = data.params.hsp_num_max;
        hitlist.blast_hit_list_update(hsp_list);
    }

    0
}

/// Port of NCBI `s_BlastHSPCollectorRun_RPS` (`hspfilter_collector.c:204`).
/// naming: Rust keeps HSPCollector/RPS as established snake_case tokens.
pub fn s_blast_hspcollector_run_rps(
    data: &mut BlastHSPCollectorData,
    results: &mut crate::hspstream::HspResults,
    query_index: i32,
    mut hsp_list: crate::hspstream::HspList,
) -> i32 {
    if hsp_list.hsps.is_empty() {
        return 0;
    }
    let query_index = query_index as usize;
    if query_index >= results.hitlists.len() {
        return -1;
    }
    let hitlist = results.hitlists[query_index].get_or_insert_with(|| {
        crate::hspstream::blast_hit_list_new(data.params.prelim_hitlist_size)
    });

    hsp_list.hsps.sort_by(|a, b| {
        a.context
            .cmp(&b.context)
            .then_with(|| crate::hspstream::score_compare_hsps(a, b))
    });

    let mut index = 0usize;
    while index < hsp_list.hsps.len() {
        let oid = hsp_list.hsps[index].context;
        let mut split_list = crate::hspstream::blast_hsp_list_new(data.params.hsp_num_max);
        split_list.oid = oid;
        while index < hsp_list.hsps.len() && hsp_list.hsps[index].context == oid {
            let mut hsp = hsp_list.hsps[index].clone();
            hsp.context = 0;
            let _ = crate::hspstream::blast_hsp_list_save_hsp(&mut split_list, hsp);
            index += 1;
        }
        split_list
            .hsps
            .sort_by(crate::hspstream::score_compare_hsps);
        hitlist.blast_hit_list_update(split_list);
    }

    0
}

/// Rust ownership equivalent of `s_BlastHSPCollectorFree`.
/// naming: Rust keeps HSPCollector as the established `hspcollector` token.
pub fn s_blast_hspcollector_free(
    _: Option<BlastHSPCollectorWriter>,
) -> Option<BlastHSPCollectorWriter> {
    None
}

/// Port of NCBI `s_BlastHSPCollectorNew` (`hspfilter_collector.c:295`).
/// naming: Rust keeps HSPCollector as the established `hspcollector` token.
pub fn s_blast_hspcollector_new(params: BlastHSPCollectorParams) -> BlastHSPCollectorWriter {
    let rps_run = crate::program::blast_program_is_rps_blast(params.program);
    BlastHSPCollectorWriter {
        data: BlastHSPCollectorData { params },
        rps_run,
    }
}

/// Port of NCBI `BlastHSPCollectorInfoNew` (`hspfilter_collector.c:353`).
pub fn blast_hsp_collector_info_new(params: BlastHSPCollectorParams) -> BlastHSPCollectorInfo {
    BlastHSPCollectorInfo { params }
}

/// 1-1 port of `BlastHSPCullingPipeInfoNew` (`hspfilter_culling.c:822`).
///
/// NCBI builds a `BlastHSPPipeInfo` linked-list node holding the
/// params and a `NewFnPtr` for `s_BlastHSPCullingPipeNew`. The Rust
/// port keeps the linked-list shape so the pipe stage chain can be
/// represented faithfully.
#[derive(Debug, Clone)]
pub struct BlastHSPPipeInfo {
    pub params: BlastHSPCullingParams,
    pub next: Option<Box<BlastHSPPipeInfo>>,
}

/// NCBI: BlastHSPCullingPipeInfoNew (hspfilter_culling.c).
pub fn blast_hsp_culling_pipe_info_new(params: BlastHSPCullingParams) -> BlastHSPPipeInfo {
    BlastHSPPipeInfo { params, next: None }
}

/// 1-1 port of `BlastHSPPipeInfo_Add` (`blast_hspstream.c:762`) for the Rust
/// linked-list representation. C appends `node` to `*head` when present or
/// installs it as the head otherwise, then returns the inserted node pointer.
pub fn blast_hsp_pipe_info_add(
    head: &mut Option<Box<BlastHSPPipeInfo>>,
    mut node: Box<BlastHSPPipeInfo>,
) -> *mut BlastHSPPipeInfo {
    let node_ptr: *mut BlastHSPPipeInfo = &mut *node;
    let mut tail = head;
    while let Some(current) = tail {
        tail = &mut current.next;
    }
    *tail = Some(node);
    node_ptr
}

/// 1-1 port of `LinkedHSP` (`hspfilter_culling.c:52`). Singly-linked
/// HSP node owning its `Hsp` payload, with culling bookkeeping
/// (`begin`, `end` query-plus-strand offsets and `merit` countdown).
///
/// NCBI uses a raw `BlastHSP*`; the Rust port owns the `Hsp` value.
#[derive(Debug, Clone)]
pub struct LinkedHsp {
    pub hsp: Hsp,
    /// Context id (NCBI: `cid`).
    pub cid: i32,
    /// Subject OID (NCBI: `sid`).
    pub sid: i32,
    /// Query offset on the plus strand.
    pub begin: i32,
    /// Query end on the plus strand.
    pub end: i32,
    /// How many list elements still dominate this one. The HSP is
    /// dropped when `merit <= 0`. NCBI initializes it to a
    /// caller-supplied threshold.
    pub merit: i32,
    pub next: Option<Box<LinkedHsp>>,
}

impl LinkedHsp {
    /// 1-1 port of `s_HSPCopy` (`hspfilter_culling.c:65`). Produces a
    /// detached copy (next pointer cleared in callers' usage).
    /// naming: Rust exposes this as an associated copy method on `LinkedHsp`.
    pub fn s_hsp_copy(&self) -> LinkedHsp {
        LinkedHsp {
            hsp: self.hsp.clone(),
            cid: self.cid,
            sid: self.sid,
            begin: self.begin,
            end: self.end,
            merit: self.merit,
            next: None,
        }
    }
}

/// 1-1 port boundary for `s_HSPCopy` (`hspfilter_culling.c:65`).
/// NCBI: s_HSPCopy (hspfilter_culling.c:65).
/// naming: Rust keeps HSP as a separate snake_case token.
pub fn s_hsp_copy(hsp: &LinkedHsp) -> LinkedHsp {
    hsp.s_hsp_copy()
}

/// Rust ownership equivalent of `s_HSPFree` (`hspfilter_culling.c:72`).
pub fn s_hsp_free(_: Option<Box<LinkedHsp>>) -> Option<Box<LinkedHsp>> {
    None
}

/// 1-1 port of `s_DominateTest` (`hspfilter_culling.c:79`).
///
/// Returns `true` iff `p` dominates `y`. The dominance criterion is
/// a 50%-overlap precondition followed by NCBI's score+length
/// formula `d = 4*s1*l1 + 2*s1*l2 - 2*s2*l1 - 4*s2*l2`. On exact ties
/// (identical score/begin/length) a deterministic tie-breaker is
/// applied (`s1 > s2`, then `sid` ascending, then subject offset
/// ascending), matching NCBI's tie-breaking exactly.
pub fn s_dominate_test(p: &LinkedHsp, y: &LinkedHsp) -> bool {
    // C uses Int8 throughout to keep intermediate products from
    // overflowing on long alignments.
    let b1 = p.begin as i64;
    let b2 = y.begin as i64;
    let e1 = p.end as i64;
    let e2 = y.end as i64;
    let s1 = p.hsp.score as i64;
    let s2 = y.hsp.score as i64;
    let l1 = e1 - b1;
    let l2 = e2 - b2;
    let overlap = e1.min(e2) - b1.max(b2);

    // C: `if (2*overlap < l2) return FALSE;` — < 50% overlap of y by p.
    if 2 * overlap < l2 {
        return false;
    }

    // C: `d = 4*s1*l1 + 2*s1*l2 - 2*s2*l1 - 4*s2*l2;`.
    let d = 4 * s1 * l1 + 2 * s1 * l2 - 2 * s2 * l1 - 4 * s2 * l2;
    let identical = s1 == s2 && b1 == b2 && l1 == l2;
    if identical || d == 0 {
        if s1 != s2 {
            return s1 > s2;
        }
        if p.sid != y.sid {
            return p.sid < y.sid;
        }
        if p.hsp.subject_offset > y.hsp.subject_offset {
            return false;
        }
        return true;
    }
    if d < 0 {
        return false;
    }
    true
}

/// 1-1 port of `s_FullPass` (`hspfilter_culling.c:123`).
///
/// Walks `list` once and decrements `y.merit` every time a list member
/// dominates `y`. Returns `false` (HSP should be dropped) as soon as
/// `merit <= 0`, otherwise `true`.
pub fn s_full_pass(list: &Option<Box<LinkedHsp>>, y: &mut LinkedHsp) -> bool {
    let mut cur = list.as_deref();
    while let Some(node) = cur {
        if s_dominate_test(node, y) {
            y.merit -= 1;
            if y.merit <= 0 {
                return false;
            }
        }
        cur = node.next.as_deref();
    }
    true
}

/// 1-1 port of `s_AddHSPtoList` (`hspfilter_culling.c:193`). Pushes
/// `y` onto the head of `list`.
/// naming: Rust keeps HSP as a separate snake_case token.
pub fn s_add_hsp_to_list(list: &mut Option<Box<LinkedHsp>>, mut y: Box<LinkedHsp>) {
    y.next = list.take();
    *list = Some(y);
}

/// blast-rs: Rust-owned implementation boundary for `s_ProcessHSPList`.
///
/// NCBI skips exactly the inserted HSP pointer (`r != y`). Rust callers pass
/// `skip_first_match` only for the host list where the just-inserted node is
/// known to be present at the head; tree/subtree walks pass `false` so
/// distinct value-identical HSPs still get processed.
fn s_process_hsp_list_impl(
    list: &mut Option<Box<LinkedHsp>>,
    y: &LinkedHsp,
    skip_first_match: bool,
) -> i32 {
    /// blast-rs: Local identity predicate replacing C pointer equality; not a
    /// direct NCBI C port.
    fn matches_y(node: &LinkedHsp, y: &LinkedHsp) -> bool {
        node.begin == y.begin
            && node.end == y.end
            && node.hsp.score == y.hsp.score
            && node.sid == y.sid
            && node.cid == y.cid
    }
    let mut num = 0i32;
    let mut skipped_self = !skip_first_match;
    // Walk in place; on removal, we splice the chain. Use a current
    // owned slot pattern to keep ownership clear.
    let mut cursor = list;
    while cursor.is_some() {
        let head_is_y = if skipped_self {
            false
        } else {
            cursor.as_ref().map(|n| matches_y(n, y)).unwrap_or(false)
        };
        if head_is_y {
            skipped_self = true;
        }
        // Decide whether to drop the head node based on dominance.
        let drop_head = if !head_is_y {
            if let Some(node) = cursor.as_mut() {
                if s_dominate_test(y, node) {
                    node.merit -= 1;
                    node.merit <= 0
                } else {
                    false
                }
            } else {
                false
            }
        } else {
            false
        };
        if drop_head {
            // Remove and free this node, advance.
            let mut taken = cursor.take().expect("head present");
            *cursor = taken.next.take();
        } else {
            num += 1;
            // Move cursor forward to the next slot.
            cursor = &mut cursor.as_mut().expect("non-null").next;
        }
    }
    num
}

/// 1-1 port of `s_ProcessHSPList` (`hspfilter_culling.c:136`).
///
/// For every list element `r`, if `y` dominates `r`, decrement `r.merit`;
/// remove `r` when `merit <= 0`. Returns the number of list elements
/// remaining. This entry point is used for lists that do not contain the
/// inserted C pointer `y`.
pub fn s_process_hsp_list(list: &mut Option<Box<LinkedHsp>>, y: &LinkedHsp) -> i32 {
    s_process_hsp_list_impl(list, y, false)
}

/// 1-1 port of `s_MarkDownHSPList` (`hspfilter_culling.c:166`).
///
/// Decrements every element's `merit`, removing those that fall to
/// `<= 0`. Returns the number of elements remaining.
pub fn s_mark_down_hsp_list(list: &mut Option<Box<LinkedHsp>>) -> i32 {
    let mut num = 0i32;
    let mut cursor = list;
    while cursor.is_some() {
        let drop_head = if let Some(node) = cursor.as_mut() {
            node.merit -= 1;
            node.merit <= 0
        } else {
            false
        };
        if drop_head {
            let mut taken = cursor.take().expect("head present");
            *cursor = taken.next.take();
        } else {
            num += 1;
            cursor = &mut cursor.as_mut().expect("non-null").next;
        }
    }
    num
}

/// `ECTreeChild` (`hspfilter_culling.c:222`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CTreeChild {
    Left,
    Right,
}

/// 1-1 port of `CTreeNode` (`hspfilter_culling.c:201`). Splays the
/// query interval `[begin, end)` into a binary tree; each node's HSP
/// list holds candidates whose query span is contained in this
/// node's interval. Built up incrementally as new HSPs arrive.
#[derive(Debug, Clone, Default)]
pub struct CTreeNode {
    pub begin: i32,
    pub end: i32,
    pub left: Option<Box<CTreeNode>>,
    pub right: Option<Box<CTreeNode>>,
    pub hsp_list: Option<Box<LinkedHsp>>,
}

/// Rust allocation equivalent of `s_GetNode` (`hspfilter_culling.c:210`).
/// blast-rs: Rust allocation helper mirroring the C node pool boundary; not a
/// direct NCBI C port.
pub fn s_get_node() -> Box<CTreeNode> {
    Box::new(CTreeNode::default())
}

/// Rust ownership equivalent of `s_RetNode` (`hspfilter_culling.c:214`).
pub fn s_ret_node(_: Option<Box<CTreeNode>>) -> Option<Box<CTreeNode>> {
    None
}

impl CTreeNode {
    /// 1-1 port of `s_ForkChildren` (`hspfilter_culling.c:258`).
    /// naming: Rust exposes this as an associated method on `CTreeNode`.
    ///
    /// Walks the node's HSP list. Each element whose query span lies
    /// strictly to the left of the midpoint is moved into a freshly
    /// created left child; elements strictly to the right go into the
    /// right child; spanning elements stay on the node itself. NCBI
    /// uses the `eLeft`/`eRight` constants for its child-direction
    /// enum; we use [`CTreeChild`].
    pub fn fork_children(&mut self) {
        s_fork_children(self);
    }

    /// 1-1 port of `s_CTreeNodeNew` (`hspfilter_culling.c:228`).
    /// naming: Rust exposes this as an associated constructor on `CTreeNode`.
    /// Allocates a new node; if `parent` is `None`, the interval stays
    /// uninitialized (caller fills `begin`/`end` for the root). When a
    /// parent is supplied, the child's interval is `[parent.begin,
    /// midpt)` (left) or `[midpt, parent.end)` (right) where
    /// `midpt = (parent.begin + parent.end) / 2`.
    pub fn s_ctree_node_new(parent: Option<&CTreeNode>, dir: CTreeChild) -> Box<CTreeNode> {
        let mut node = s_get_node();
        if let Some(p) = parent {
            let midpt = (p.begin + p.end) / 2;
            match dir {
                CTreeChild::Left => {
                    node.begin = p.begin;
                    node.end = midpt;
                }
                CTreeChild::Right => {
                    node.begin = midpt;
                    node.end = p.end;
                }
            }
        }
        node
    }
}

/// 1-1 port boundary for `s_CTreeNodeNew` (`hspfilter_culling.c:228`).
/// NCBI: s_CTreeNodeNew (hspfilter_culling.c:228).
/// naming: Rust keeps CTree as the established `ctree` token.
pub fn s_ctree_node_new(parent: Option<&CTreeNode>, dir: CTreeChild) -> Box<CTreeNode> {
    CTreeNode::s_ctree_node_new(parent, dir)
}

/// 1-1 port boundary for `s_ForkChildren` (`hspfilter_culling.c:258`).
/// NCBI: s_ForkChildren (hspfilter_culling.c:258).
/// naming: Rust keeps the free helper name used by existing callers.
pub fn s_fork_children(node: &mut CTreeNode) {
    debug_assert!(node.left.is_none());
    debug_assert!(node.right.is_none());
    let midpt = (node.begin + node.end) / 2;

    // Take the entire list off the parent so we can sort each entry
    // into one of three buckets (stays / left / right). Reassemble
    // at the end. NCBI advances `p` before moving `r` to a child,
    // then leaves the predecessor `q` unchanged when a node moves.
    let mut original = node.hsp_list.take();
    let mut stays: Option<Box<LinkedHsp>> = None;
    let mut left_list: Option<Box<LinkedHsp>> = None;
    let mut right_list: Option<Box<LinkedHsp>> = None;

    /// blast-rs: Local linked-list append helper; not a direct NCBI C port.
    fn append_tail(list: &mut Option<Box<LinkedHsp>>, mut node: Box<LinkedHsp>) {
        node.next = None;
        match list.as_mut() {
            None => *list = Some(node),
            Some(head) => {
                let mut cursor: &mut Box<LinkedHsp> = head;
                while cursor.next.is_some() {
                    cursor = cursor.next.as_mut().unwrap();
                }
                cursor.next = Some(node);
            }
        }
    }

    while let Some(mut hsp_node) = original.take() {
        original = hsp_node.next.take();
        if hsp_node.end < midpt {
            append_tail(&mut left_list, hsp_node);
        } else if hsp_node.begin > midpt {
            append_tail(&mut right_list, hsp_node);
        } else {
            append_tail(&mut stays, hsp_node);
        }
    }

    if left_list.is_some() {
        let mut child = s_ctree_node_new(Some(node), CTreeChild::Left);
        let mut child_list = left_list;
        while let Some(mut hsp_node) = child_list.take() {
            child_list = hsp_node.next.take();
            s_add_hsp_to_list(&mut child.hsp_list, hsp_node);
        }
        node.left = Some(child);
    }

    if right_list.is_some() {
        let mut child = s_ctree_node_new(Some(node), CTreeChild::Right);
        let mut child_list = right_list;
        while let Some(mut hsp_node) = child_list.take() {
            child_list = hsp_node.next.take();
            s_add_hsp_to_list(&mut child.hsp_list, hsp_node);
        }
        node.right = Some(child);
    }

    node.hsp_list = stays;
}

/// Rust ownership equivalent of `s_CTreeNodeFree` (`hspfilter_culling.c:250`).
/// naming: Rust keeps CTree as the established `ctree` token.
pub fn s_ctree_node_free(node: Option<Box<CTreeNode>>) -> Option<Box<CTreeNode>> {
    if let Some(node) = node.as_ref() {
        debug_assert!(node.left.is_none());
        debug_assert!(node.right.is_none());
        debug_assert!(node.hsp_list.is_none());
    }
    s_ret_node(node)
}

/// Disabled NCBI debug hook (`s_Debug`, `hspfilter_culling.c:302`).
/// blast-rs: Disabled debug traversal hook; not a direct NCBI C port.
pub fn s_debug(node: Option<&CTreeNode>) {
    let Some(node) = node else {
        return;
    };

    let mut hsp = node.hsp_list.as_deref();
    while let Some(current) = hsp {
        let _ = (current.begin, current.end, current.hsp.score, current.merit);
        hsp = current.next.as_deref();
    }

    s_debug(node.left.as_deref());
    s_debug(node.right.as_deref());
}

/// 1-1 port of `s_MarkDownCTree` (`hspfilter_culling.c:319`).
/// naming: Rust keeps CTree as the established `ctree` token.
///
/// Recursively walks the tree, decrementing every HSP's `merit` and
/// removing nodes whose lists are emptied AND that have no children.
/// Caller passes the slot owning the tree (or subtree).
pub fn s_mark_down_ctree(slot: &mut Option<Box<CTreeNode>>) {
    let Some(node) = slot.as_mut() else { return };
    s_mark_down_ctree(&mut node.left);
    s_mark_down_ctree(&mut node.right);
    let remaining = s_mark_down_hsp_list(&mut node.hsp_list);
    if remaining <= 0 && node.left.is_none() && node.right.is_none() {
        *slot = None;
    }
}

/// 1-1 port of `s_ProcessCTree` (`hspfilter_culling.c:334`).
/// naming: Rust keeps CTree as the established `ctree` token.
///
/// Walks the tree to update merit in response to a newly-added
/// dominator `x`. If `x` covers the full range of a subtree, every
/// element there gets decremented (`s_mark_down_ctree`). Otherwise
/// the recursion descends into the half(s) that overlap `x`.
pub fn s_process_ctree(slot: &mut Option<Box<CTreeNode>>, x: &LinkedHsp) {
    let node = match slot.as_mut() {
        Some(n) => n,
        None => return,
    };
    // x covers full range → decrement everywhere.
    if x.begin <= node.begin && x.end >= node.end {
        s_mark_down_ctree(slot);
        return;
    }
    // Leaf: just process the local list and clean up if emptied.
    if node.left.is_none() && node.right.is_none() {
        if s_process_hsp_list(&mut node.hsp_list, x) <= 0 {
            *slot = None;
        }
        return;
    }
    // Internal: descend into the side(s) that overlap x.
    let midpt = (node.begin + node.end) / 2;
    if x.end < midpt {
        s_process_ctree(&mut node.left, x);
    } else if x.begin > midpt {
        s_process_ctree(&mut node.right, x);
    } else {
        s_process_ctree(&mut node.left, x);
        s_process_ctree(&mut node.right, x);
        if s_process_hsp_list(&mut node.hsp_list, x) <= 0
            && node.left.is_none()
            && node.right.is_none()
        {
            *slot = None;
        }
    }
}

/// 1-1 port of `s_CTreeNew` (`hspfilter_culling.c:376`). Root over
/// `[0, qlen)`.
/// naming: Rust keeps CTree as the established `ctree` token.
pub fn s_ctree_new(qlen: i32) -> Box<CTreeNode> {
    let mut tree = s_ctree_node_new(None, CTreeChild::Left);
    tree.begin = 0;
    tree.end = qlen;
    tree
}

/// Rust ownership equivalent of `s_CTreeFree` (`hspfilter_culling.c:384`).
/// naming: Rust keeps CTree as the established `ctree` token.
pub fn s_ctree_free(tree: Option<Box<CTreeNode>>) -> Option<Box<CTreeNode>> {
    let Some(mut tree) = tree else { return None };
    debug_assert!(tree.hsp_list.is_none());
    let _ = s_ctree_free(tree.left.take());
    let _ = s_ctree_free(tree.right.take());
    s_ctree_node_free(Some(tree))
}

/// 1-1 port of `s_RipHSPOffCTree` (`hspfilter_culling.c:396`).
/// naming: Rust keeps HSP and CTree as readable snake_case tokens.
///
/// Recursively rips every HSP off the tree (rooted at `slot`) into a
/// single linked list. NCBI assumes the caller owns the tree; the
/// Rust port consumes the tree (caller passes `Option::take()`-style
/// ownership) and returns the linked list.
pub fn s_rip_hsp_off_ctree(slot: Option<Box<CTreeNode>>) -> Option<Box<LinkedHsp>> {
    let Some(mut node) = slot else { return None };
    let mut q = node.hsp_list.take();
    let left_list = s_rip_hsp_off_ctree(node.left.take());
    let right_list = s_rip_hsp_off_ctree(node.right.take());

    /// blast-rs: Local linked-list append helper; not a direct NCBI C port.
    fn append_tail(list: &mut Option<Box<LinkedHsp>>, tail: Option<Box<LinkedHsp>>) {
        if list.is_none() {
            *list = tail;
            return;
        }
        let mut cursor = list.as_mut().unwrap();
        while cursor.next.is_some() {
            cursor = cursor.next.as_mut().unwrap();
        }
        cursor.next = tail;
    }
    append_tail(&mut q, left_list);
    append_tail(&mut q, right_list);
    q
}

/// `kNumHSPtoFork` (`hspfilter_culling.c:435`). Number of HSPs at a
/// leaf node that triggers a `fork_children` split.
pub const K_NUM_HSP_TO_FORK: i32 = 20;

/// 1-1 port of `s_SaveHSP` (`hspfilter_culling.c:430`). Returns
/// `false` when the candidate is dominated and should be dropped,
/// `true` after inserting (and possibly forking the host node).
///
/// `a` is consumed only when inserted; on `false` the caller can
/// reuse / drop it.
pub fn s_save_hsp(tree: &mut CTreeNode, a: &mut LinkedHsp) -> bool {
    // C uses a single `tree` pointer that's advanced down the tree;
    // we recreate that with raw pointer-chasing inside a loop. Rust
    // doesn't make this clean across owned children, so we recurse
    // using a path-collecting helper.
    /// blast-rs: Recursive descent helper for Rust-owned culling trees; not a
    /// direct NCBI C port.
    fn descend<'a>(node: &'a mut CTreeNode, a: &mut LinkedHsp) -> Option<&'a mut CTreeNode> {
        debug_assert!(node.begin <= a.begin);
        debug_assert!(node.end >= a.end);
        if !s_full_pass(&node.hsp_list, a) {
            return None;
        }
        let midpt = (node.begin + node.end) / 2;
        if a.end < midpt && node.left.is_some() {
            descend(node.left.as_mut().unwrap(), a)
        } else if a.begin > midpt && node.right.is_some() {
            descend(node.right.as_mut().unwrap(), a)
        } else {
            Some(node)
        }
    }
    let host = match descend(tree, a) {
        Some(n) => n,
        None => return false,
    };

    // Insert a copy at the host node.
    let copy = s_hsp_copy(a);
    s_add_hsp_to_list(&mut host.hsp_list, Box::new(copy));

    // Take the freshly-inserted head as the dominator reference.
    // NCBI's C uses `x` (the inserted node), but in Rust it's safer to
    // work with the pre-insert clone since `host.hsp_list`'s head is
    // now `x`.
    let x = s_hsp_copy(a);

    if host.left.is_none() && host.right.is_none() {
        if s_process_hsp_list_impl(&mut host.hsp_list, &x, true) >= K_NUM_HSP_TO_FORK {
            host.fork_children();
        }
        return true;
    }
    s_process_hsp_list_impl(&mut host.hsp_list, &x, true);
    s_process_ctree(&mut host.left, &x);
    s_process_ctree(&mut host.right, &x);
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mk(score: i32, begin: i32, end: i32, sid: i32) -> Box<LinkedHsp> {
        Box::new(LinkedHsp {
            hsp: Hsp {
                score,
                num_ident: 0,
                bit_score: 0.0,
                evalue: 0.0,
                query_offset: begin,
                query_end: end,
                query_gapped_start: begin,
                subject_offset: 0,
                subject_end: end - begin,
                subject_gapped_start: 0,
                context: 0,
                query_frame: 0,
                subject_frame: 0,
                num_gaps: 0,
                comp_adjustment_method: 0,
                edit_script: None,
                pat_info: None,
                map_info: None,
            },
            cid: 0,
            sid: sid,
            begin,
            end,
            merit: 1,
            next: None,
        })
    }

    fn hsp(score: i32, context: i32, query_offset: i32, subject_offset: i32) -> Hsp {
        Hsp {
            score,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1.0 / score.max(1) as f64,
            query_offset,
            query_end: query_offset + 10,
            query_gapped_start: query_offset,
            subject_offset,
            subject_end: subject_offset + 10,
            subject_gapped_start: subject_offset,
            context,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        }
    }

    #[test]
    fn dominate_test_higher_score_dominates_lower_when_overlapping() {
        let p = mk(100, 0, 100, 0);
        let y = mk(50, 10, 90, 1);
        // p covers [0,100], y covers [10,90] (80 wide). Overlap = 80,
        // 2*80 = 160 >= l2 (80) → 50% overlap satisfied. Higher score
        // dominates.
        assert!(s_dominate_test(&p, &y));
    }

    #[test]
    fn dominate_test_no_overlap_returns_false() {
        let p = mk(100, 0, 50, 0);
        let y = mk(50, 60, 100, 1);
        // Zero overlap → fails the 50% precondition.
        assert!(!s_dominate_test(&p, &y));
    }

    #[test]
    fn dominate_test_tied_uses_sid() {
        // Identical begin/end/score → tiebreak by sid (lower wins).
        let p = mk(100, 0, 100, 5);
        let y = mk(100, 0, 100, 9);
        assert!(s_dominate_test(&p, &y)); // 5 < 9 → p dominates
        assert!(!s_dominate_test(&y, &p));
    }

    #[test]
    fn full_pass_drops_when_merit_reaches_zero() {
        let mut list: Option<Box<LinkedHsp>> = None;
        s_add_hsp_to_list(&mut list, mk(100, 0, 100, 0));
        s_add_hsp_to_list(&mut list, mk(90, 0, 100, 1));
        // y has merit 1 and is dominated by both → fails on first hit.
        let mut y = *mk(50, 10, 90, 2);
        y.merit = 1;
        assert!(!s_full_pass(&list, &mut y));
        assert_eq!(y.merit, 0);
    }

    #[test]
    fn full_pass_keeps_when_merit_survives() {
        let mut list: Option<Box<LinkedHsp>> = None;
        s_add_hsp_to_list(&mut list, mk(100, 0, 100, 0));
        // y has merit 5 and only one dominator → still alive.
        let mut y = *mk(50, 10, 90, 1);
        y.merit = 5;
        assert!(s_full_pass(&list, &mut y));
        assert_eq!(y.merit, 4);
    }

    #[test]
    fn process_hsp_list_drops_dominated_elements() {
        let mut list: Option<Box<LinkedHsp>> = None;
        s_add_hsp_to_list(&mut list, mk(50, 10, 90, 1)); // dominated by y
        s_add_hsp_to_list(&mut list, mk(60, 20, 80, 2)); // dominated by y
        let y = *mk(200, 0, 100, 0);
        let remaining = s_process_hsp_list(&mut list, &y);
        assert_eq!(remaining, 0);
        assert!(list.is_none());
    }

    #[test]
    fn process_hsp_list_splices_removed_middle_and_counts_survivors() {
        let mut list: Option<Box<LinkedHsp>> = None;
        s_add_hsp_to_list(&mut list, mk(50, 300, 400, 3)); // no overlap
        s_add_hsp_to_list(&mut list, mk(50, 10, 90, 2)); // dominated by y
        s_add_hsp_to_list(&mut list, mk(50, 200, 300, 1)); // no overlap
        let y = *mk(200, 0, 100, 0);

        let remaining = s_process_hsp_list(&mut list, &y);

        assert_eq!(remaining, 2);
        let head = list.as_ref().expect("first survivor");
        assert_eq!(head.sid, 1);
        let tail = head.next.as_ref().expect("second survivor");
        assert_eq!(tail.sid, 3);
        assert!(tail.next.is_none());
    }

    #[test]
    fn process_hsp_list_host_skip_only_ignores_first_value_match() {
        let mut list: Option<Box<LinkedHsp>> = None;
        s_add_hsp_to_list(&mut list, mk(100, 0, 100, 1));
        s_add_hsp_to_list(&mut list, mk(100, 0, 100, 1));
        let y = *mk(100, 0, 100, 1);

        let remaining = s_process_hsp_list_impl(&mut list, &y, true);

        assert_eq!(remaining, 1);
        assert!(list.as_ref().is_some_and(|node| node.next.is_none()));
    }

    #[test]
    fn mark_down_hsp_list_drops_zero_merit_entries() {
        let mut list: Option<Box<LinkedHsp>> = None;
        let mut a = mk(100, 0, 100, 0);
        a.merit = 1; // will fall to 0 → dropped
        s_add_hsp_to_list(&mut list, a);
        let mut b = mk(100, 0, 100, 1);
        b.merit = 5;
        s_add_hsp_to_list(&mut list, b);
        let remaining = s_mark_down_hsp_list(&mut list);
        assert_eq!(remaining, 1);
        // Remaining one is the merit-5 entry decremented to 4.
        let head = list.as_ref().expect("head");
        assert_eq!(head.merit, 4);
    }

    #[test]
    fn blast_hsp_culling_pipe_run_rebuilds_results_in_place() {
        // Build a 1-context query info.
        let qi = crate::queryinfo::QueryInfo {
            num_queries: 1,
            contexts: vec![crate::queryinfo::ContextInfo {
                query_offset: 0,
                query_length: 1000,
                eff_searchsp: 0,
                length_adjustment: 0,
                query_index: 0,
                frame: 0,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: 1000,
            min_length: 0,
        };
        let hit = crate::options::HitSavingOptions::default();
        let cull_opts = BlastHSPCullingOptions { max_hits: 5 };
        let params =
            blast_hsp_culling_params_new(&hit, &cull_opts, 10, crate::program::BLASTP, 0, true);
        let mut data = blast_hsp_culling_pipe_new(params, &qi);

        // Pre-populate results: one hitlist with two HSPs in random
        // evalue order. The pipe-run should sort + repackage.
        let mut results = crate::hspstream::HspResults::new(1);
        let mut list = crate::hspstream::HspList::new(99);
        list.add_hsp(Hsp {
            score: 200,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1e-20,
            query_offset: 0,
            query_end: 100,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: 100,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        list.add_hsp(Hsp {
            score: 50,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 0.5, // higher e-value, will be re-sorted
            query_offset: 200,
            query_end: 300,
            query_gapped_start: 200,
            subject_offset: 200,
            subject_end: 300,
            subject_gapped_start: 200,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        let mut hitlist = crate::hspstream::HitList::new();
        hitlist.hsp_lists.push(list);
        results.hitlists[0] = Some(hitlist);

        blast_hsp_culling_pipe_run(&mut data, &mut results);

        // After pipe-run: hitlist should be re-built; the OID is
        // preserved.
        let hl = results.hitlists[0].as_ref().expect("hitlist");
        assert_eq!(hl.hsp_lists.len(), 1);
        assert_eq!(hl.hsp_lists[0].oid, 99);
        // Best evalue should reflect the lowest of the input HSPs.
        assert!(hl.hsp_lists[0].best_evalue < 1e-15);
    }

    #[test]
    fn blast_hsp_culling_pipe_free_clears_slot() {
        let qi = crate::queryinfo::QueryInfo {
            num_queries: 1,
            contexts: vec![],
            max_length: 0,
            min_length: 0,
        };
        let hit = crate::options::HitSavingOptions::default();
        let cull_opts = BlastHSPCullingOptions { max_hits: 5 };
        let params =
            blast_hsp_culling_params_new(&hit, &cull_opts, 10, crate::program::BLASTP, 0, true);
        let mut slot = Some(blast_hsp_culling_pipe_new(params, &qi));
        blast_hsp_culling_pipe_free(&mut slot);
        assert!(slot.is_none());
    }

    #[test]
    fn s_blast_hspculling_free_clears_initialized_state() {
        let qi = crate::queryinfo::QueryInfo {
            num_queries: 1,
            contexts: vec![crate::queryinfo::ContextInfo {
                query_offset: 0,
                query_length: 10,
                eff_searchsp: 0,
                length_adjustment: 0,
                query_index: 0,
                frame: 0,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: 10,
            min_length: 0,
        };
        let hit = crate::options::HitSavingOptions::default();
        let cull_opts = BlastHSPCullingOptions { max_hits: 5 };
        let params =
            blast_hsp_culling_params_new(&hit, &cull_opts, 10, crate::program::BLASTP, 0, true);
        let mut data = BlastHSPCullingData::s_blast_hspculling_new(params, &qi);
        data.s_blast_hspculling_init();
        assert_eq!(data.c_tree.len(), 1);

        let mut slot = Some(data);
        s_blast_hspculling_free(&mut slot);
        assert!(slot.is_none());
    }

    #[test]
    fn culling_data_init_run_finalize_round_trip() {
        // Build a 1-context query info.
        let qi = crate::queryinfo::QueryInfo {
            num_queries: 1,
            contexts: vec![crate::queryinfo::ContextInfo {
                query_offset: 0,
                query_length: 1000,
                eff_searchsp: 0,
                length_adjustment: 0,
                query_index: 0,
                frame: 0,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: 1000,
            min_length: 0,
        };
        let hit = crate::options::HitSavingOptions::default();
        let cull_opts = BlastHSPCullingOptions { max_hits: 5 };
        let params =
            blast_hsp_culling_params_new(&hit, &cull_opts, 10, crate::program::BLASTP, 0, true);
        let mut data = BlastHSPCullingData::s_blast_hspculling_new(params, &qi);
        data.s_blast_hspculling_init();
        // Build an HSP list with a few HSPs.
        let mut list = crate::hspstream::HspList::new(42);
        list.add_hsp(Hsp {
            score: 100,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1e-10,
            query_offset: 0,
            query_end: 100,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: 100,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        list.add_hsp(Hsp {
            score: 50,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1e-5,
            query_offset: 200,
            query_end: 300,
            query_gapped_start: 200,
            subject_offset: 200,
            subject_end: 300,
            subject_gapped_start: 200,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        data.s_blast_hspculling_run(&mut list);
        // Finalize and check the shape.
        let results = data.s_blast_hspculling_final();
        let hl = results.hitlists[0].as_ref().expect("hitlist");
        assert_eq!(hl.hsp_lists.len(), 1);
        assert_eq!(hl.hsp_lists[0].oid, 42);
        // Both HSPs should have survived (no domination).
        assert_eq!(hl.hsp_lists[0].hsps.len(), 2);
        assert_eq!(hl.low_score, 100);
        assert_eq!(hl.worst_evalue, 1e-10);
    }

    #[test]
    fn blast_hsp_culling_params_new_pipes_options() {
        let hit = crate::options::HitSavingOptions {
            expect_value: 10.0,
            hitlist_size: 500,
            cutoff_score: 0,
            percent_identity: 0.0,
            min_hit_length: 0,
            ..crate::options::HitSavingOptions::default()
        };
        let cull = BlastHSPCullingOptions { max_hits: 5 };
        let params =
            blast_hsp_culling_params_new(&hit, &cull, 100, crate::program::BLASTP, 0, true);
        assert_eq!(params.culling_max, 5);
        assert_eq!(params.prelim_hitlist_size, 550);
        assert_eq!(params.hsp_num_max, 100);
        assert_eq!(params.program, crate::program::BLASTP);
    }

    #[test]
    fn blast_hsp_culling_params_new_doubles_for_compo_adjust() {
        let hit = crate::options::HitSavingOptions {
            expect_value: 10.0,
            hitlist_size: 500,
            cutoff_score: 0,
            percent_identity: 0.0,
            min_hit_length: 0,
            ..crate::options::HitSavingOptions::default()
        };
        let cull = BlastHSPCullingOptions { max_hits: 5 };
        let params =
            blast_hsp_culling_params_new(&hit, &cull, 100, crate::program::BLASTP, 2, true);
        assert_eq!(params.prelim_hitlist_size, 1050);
    }

    #[test]
    fn blast_hsp_culling_info_and_pipe_info_carry_params() {
        let hit = crate::options::HitSavingOptions {
            expect_value: 10.0,
            hitlist_size: 200,
            cutoff_score: 0,
            percent_identity: 0.0,
            min_hit_length: 0,
            ..crate::options::HitSavingOptions::default()
        };
        let cull = BlastHSPCullingOptions { max_hits: 3 };
        let params =
            blast_hsp_culling_params_new(&hit, &cull, 50, crate::program::BLASTN, 0, false);
        let info = blast_hsp_culling_info_new(params.clone());
        assert_eq!(info.params.culling_max, 3);
        let pipe = blast_hsp_culling_pipe_info_new(params);
        assert_eq!(pipe.params.hsp_num_max, 50);
        assert!(pipe.next.is_none());
    }

    #[test]
    fn blast_hsp_pipe_info_add_appends_to_tail() {
        let hit = crate::options::HitSavingOptions::default();
        let cull = BlastHSPCullingOptions { max_hits: 3 };
        let first_params =
            blast_hsp_culling_params_new(&hit, &cull, 10, crate::program::BLASTN, 0, false);
        let second_params =
            blast_hsp_culling_params_new(&hit, &cull, 20, crate::program::BLASTP, 0, false);
        let third_params =
            blast_hsp_culling_params_new(&hit, &cull, 30, crate::program::TBLASTN, 0, false);
        let mut head = None;

        let first_ptr = blast_hsp_pipe_info_add(
            &mut head,
            Box::new(blast_hsp_culling_pipe_info_new(first_params)),
        );
        let second_ptr = blast_hsp_pipe_info_add(
            &mut head,
            Box::new(blast_hsp_culling_pipe_info_new(second_params)),
        );
        let third_ptr = blast_hsp_pipe_info_add(
            &mut head,
            Box::new(blast_hsp_culling_pipe_info_new(third_params)),
        );

        let first = head.as_ref().expect("head");
        let second = first.next.as_ref().expect("second");
        let third = second.next.as_ref().expect("third");
        assert_eq!(first.params.hsp_num_max, 10);
        assert_eq!(second.params.hsp_num_max, 20);
        assert_eq!(third.params.hsp_num_max, 30);
        assert!(third.next.is_none());
        assert_eq!(first_ptr, first.as_ref() as *const _ as *mut _);
        assert_eq!(second_ptr, second.as_ref() as *const _ as *mut _);
        assert_eq!(third_ptr, third.as_ref() as *const _ as *mut _);
    }

    #[test]
    fn blast_hsp_culling_params_free_clears_slot() {
        let hit = crate::options::HitSavingOptions::default();
        let cull = BlastHSPCullingOptions { max_hits: 5 };
        let mut slot = Some(blast_hsp_culling_params_new(
            &hit,
            &cull,
            10,
            crate::program::BLASTP,
            0,
            false,
        ));
        blast_hsp_culling_params_free(&mut slot);
        assert!(slot.is_none());
    }

    #[test]
    fn blast_hsp_best_hit_params_info_and_free_match_c_shape() {
        let hit = crate::options::HitSavingOptions {
            hitlist_size: 200,
            hsp_num_max: 7,
            program_number: crate::program::BLASTP,
            ..crate::options::HitSavingOptions::default()
        };
        let opts = BlastHSPBestHitOptions {
            overhang: 0.1,
            score_edge: 0.1,
        };
        let params = blast_hsp_best_hit_params_new(&hit, &opts, 0, true);
        assert_eq!(params.prelim_hitlist_size, 250);
        assert_eq!(params.hsp_num_max, 7);
        assert_eq!(params.program, crate::program::BLASTP);
        assert_eq!(params.overhang, 0.1);
        assert_eq!(params.score_edge, 0.1);

        let info = blast_hsp_best_hit_info_new(params.clone());
        assert_eq!(info.params.hsp_num_max, 7);
        let pipe_info = blast_hsp_best_hit_pipe_info_new(params.clone());
        assert!(pipe_info.next.is_none());
        assert_eq!(pipe_info.params.prelim_hitlist_size, 250);

        let mut slot = Some(params);
        assert!(blast_hsp_best_hit_params_free(&mut slot).is_none());
        assert!(slot.is_none());
    }

    #[test]
    fn blast_hsp_best_hit_pipe_run_filters_contained_lower_density_hsp() {
        let qi = crate::queryinfo::QueryInfo {
            num_queries: 1,
            contexts: vec![crate::queryinfo::ContextInfo {
                query_offset: 0,
                query_length: 200,
                eff_searchsp: 0,
                length_adjustment: 0,
                query_index: 0,
                frame: 0,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: 200,
            min_length: 0,
        };
        let hit = crate::options::HitSavingOptions {
            hitlist_size: 20,
            hsp_num_max: 10,
            program_number: crate::program::BLASTP,
            ..crate::options::HitSavingOptions::default()
        };
        let opts = BlastHSPBestHitOptions {
            overhang: 0.1,
            score_edge: 0.1,
        };
        let params = blast_hsp_best_hit_params_new(&hit, &opts, 0, true);
        let mut data = s_blast_hsp_best_hit_pipe_new(params, &qi);

        let mut results = crate::hspstream::HspResults::new(1);
        let mut hsp_list = crate::hspstream::HspList::new(42);
        hsp_list.add_hsp(Hsp {
            score: 200,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1e-40,
            query_offset: 0,
            query_end: 100,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: 100,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        hsp_list.add_hsp(Hsp {
            score: 50,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 1e-5,
            query_offset: 20,
            query_end: 80,
            query_gapped_start: 20,
            subject_offset: 20,
            subject_end: 80,
            subject_gapped_start: 20,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        let mut hitlist = crate::hspstream::HitList::new();
        hitlist.hsp_lists.push(hsp_list);
        results.hitlists[0] = Some(hitlist);

        assert_eq!(s_blast_hsp_best_hit_pipe_run(&mut data, &mut results), 0);
        let hitlist = results.hitlists[0].as_ref().expect("hitlist");
        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(hitlist.hsp_lists[0].oid, 42);
        assert_eq!(hitlist.hsp_lists[0].hsps.len(), 1);
        assert_eq!(hitlist.hsp_lists[0].hsps[0].score, 200);
    }

    #[test]
    fn blast_hsp_best_hit_run_rps_uses_hsp_context_as_subject_oid() {
        let qi = crate::queryinfo::QueryInfo {
            num_queries: 1,
            contexts: vec![crate::queryinfo::ContextInfo {
                query_offset: 0,
                query_length: 100,
                eff_searchsp: 0,
                length_adjustment: 0,
                query_index: 0,
                frame: 0,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: 100,
            min_length: 0,
        };
        let hit = crate::options::HitSavingOptions {
            hitlist_size: 20,
            hsp_num_max: 10,
            program_number: crate::program::RPS_BLAST,
            ..crate::options::HitSavingOptions::default()
        };
        let opts = BlastHSPBestHitOptions {
            overhang: 0.1,
            score_edge: 0.1,
        };
        let params = blast_hsp_best_hit_params_new(&hit, &opts, 0, true);
        let mut data = s_blast_hsp_best_hit_pipe_new(params, &qi);
        let mut results = crate::hspstream::HspResults::new(1);
        let _ = s_blast_hsp_best_hit_init(&mut data, &results);

        let mut hsp_list = crate::hspstream::HspList::new(0);
        hsp_list.add_hsp(hsp(100, 7, 10, 2));
        assert_eq!(
            s_blast_hsp_best_hit_run_rps(&mut data, 0, Some(hsp_list)),
            0
        );
        assert_eq!(s_blast_hsp_best_hit_final(&mut data, &mut results), 0);
        let hitlist = results.hitlists[0].as_ref().expect("hitlist");
        assert_eq!(hitlist.hsp_lists[0].oid, 7);
        assert_eq!(hitlist.hsp_lists[0].hsps[0].context, 0);
    }

    #[test]
    fn blast_hsp_best_hit_rps_keeps_equal_begin_candidate() {
        let qi = crate::queryinfo::QueryInfo {
            num_queries: 1,
            contexts: vec![crate::queryinfo::ContextInfo {
                query_offset: 0,
                query_length: 100,
                eff_searchsp: 0,
                length_adjustment: 0,
                query_index: 0,
                frame: 0,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: 100,
            min_length: 0,
        };
        let hit = crate::options::HitSavingOptions {
            hitlist_size: 20,
            hsp_num_max: 10,
            program_number: crate::program::RPS_BLAST,
            ..crate::options::HitSavingOptions::default()
        };
        let opts = BlastHSPBestHitOptions {
            overhang: 0.1,
            score_edge: 0.1,
        };
        let params = blast_hsp_best_hit_params_new(&hit, &opts, 0, true);
        let mut data = s_blast_hsp_best_hit_pipe_new(params, &qi);
        let mut results = crate::hspstream::HspResults::new(1);
        let _ = s_blast_hsp_best_hit_init(&mut data, &results);

        let mut hsp_list = crate::hspstream::HspList::new(0);
        hsp_list.add_hsp(hsp(100, 7, 10, 2));
        hsp_list.add_hsp(hsp(80, 7, 9, 2));
        assert_eq!(
            s_blast_hsp_best_hit_run_rps(&mut data, 0, Some(hsp_list)),
            0
        );
        assert_eq!(s_blast_hsp_best_hit_final(&mut data, &mut results), 0);
        let hitlist = results.hitlists[0].as_ref().expect("hitlist");
        assert_eq!(hitlist.hsp_lists[0].hsps.len(), 2);
        assert_eq!(
            hitlist.hsp_lists[0]
                .hsps
                .iter()
                .map(|hsp| hsp.score)
                .collect::<Vec<_>>(),
            vec![100, 80]
        );
    }

    #[test]
    fn import_from_hitlist_reverses_equal_begin_insertion_order() {
        let qi = crate::queryinfo::QueryInfo {
            num_queries: 1,
            contexts: vec![crate::queryinfo::ContextInfo {
                query_offset: 0,
                query_length: 100,
                eff_searchsp: 0,
                length_adjustment: 0,
                query_index: 0,
                frame: 0,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: 100,
            min_length: 0,
        };
        let hit = crate::options::HitSavingOptions {
            hitlist_size: 20,
            hsp_num_max: 10,
            program_number: crate::program::BLASTP,
            ..crate::options::HitSavingOptions::default()
        };
        let opts = BlastHSPBestHitOptions {
            overhang: 0.1,
            score_edge: 0.1,
        };
        let params = blast_hsp_best_hit_params_new(&hit, &opts, 0, true);
        let mut data = s_blast_hsp_best_hit_pipe_new(params, &qi);
        let results = crate::hspstream::HspResults::new(1);
        let _ = s_blast_hsp_best_hit_init(&mut data, &results);

        let mut list = crate::hspstream::HspList::new(9);
        list.add_hsp(hsp(100, 0, 10, 0));
        list.add_hsp(hsp(80, 0, 10, 20));
        let mut hitlist = crate::hspstream::HitList::new();
        hitlist.hsp_lists.push(list);

        assert_eq!(s_import_from_hitlist(0, &mut data, &mut hitlist), 0);
        assert_eq!(data.best_list[0].len(), 2);
        assert_eq!(data.best_list[0][0].hsp.score, 80);
        assert_eq!(data.best_list[0][1].hsp.score, 100);
    }

    #[test]
    fn export_to_hitlist_does_not_sort_hsps_before_final_stage() {
        let qi = crate::queryinfo::QueryInfo {
            num_queries: 1,
            contexts: vec![crate::queryinfo::ContextInfo {
                query_offset: 0,
                query_length: 100,
                eff_searchsp: 0,
                length_adjustment: 0,
                query_index: 0,
                frame: 0,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: 100,
            min_length: 0,
        };
        let hit = crate::options::HitSavingOptions {
            hitlist_size: 20,
            hsp_num_max: 10,
            program_number: crate::program::BLASTP,
            ..crate::options::HitSavingOptions::default()
        };
        let opts = BlastHSPBestHitOptions {
            overhang: 0.1,
            score_edge: 0.1,
        };
        let params = blast_hsp_best_hit_params_new(&hit, &opts, 0, true);
        let mut data = s_blast_hsp_best_hit_pipe_new(params, &qi);
        let results = crate::hspstream::HspResults::new(1);
        let _ = s_blast_hsp_best_hit_init(&mut data, &results);

        data.best_list[0].push(LinkedHspBestHit {
            hsp: hsp(60, 0, 0, 0),
            sid: 5,
            begin: 0,
            end: 10,
            len: 10,
        });
        data.best_list[0].push(LinkedHspBestHit {
            hsp: hsp(100, 0, 20, 20),
            sid: 5,
            begin: 20,
            end: 30,
            len: 10,
        });
        data.num_hsps[0] = 2;

        let mut hitlist = crate::hspstream::blast_hit_list_new(10);
        assert_eq!(s_export_to_hitlist(0, &mut data, &mut hitlist), 0);

        let list = &hitlist.hsp_lists[0];
        assert_eq!(
            list.hsps.iter().map(|hsp| hsp.score).collect::<Vec<_>>(),
            vec![60, 100]
        );
    }

    #[test]
    fn translated_hsp_filtering_options_lifecycle_and_transfers() {
        let mut filter = blast_hspfiltering_options_new();
        assert!(filter.best_hit.is_none());
        assert!(filter.subject_besthit_opts.is_none());
        assert!(filter.culling_opts.is_none());

        let mut best_hit = Some(blast_hspbest_hit_options_new(0.1, 0.2));
        assert_eq!(
            blast_hspfiltering_options_add_best_hit(
                Some(&mut filter),
                &mut best_hit,
                crate::util::EBlastStage::PrelimSearch,
            ),
            0
        );
        assert!(best_hit.is_none());
        assert_eq!(
            filter.best_hit,
            Some(BlastHSPBestHitOptions {
                overhang: 0.1,
                score_edge: 0.2,
            })
        );
        assert_eq!(
            filter.best_hit_stage,
            crate::util::EBlastStage::PrelimSearch
        );

        assert_eq!(
            blast_hsp_subject_best_hit_options_new(true),
            BlastHSPSubjectBestHitOptions { max_range_diff: 3 }
        );
        let mut subject_besthit = Some(blast_hsp_subject_best_hit_options_new(false));
        assert_eq!(
            blast_hspfiltering_options_add_subject_best_hit(
                Some(&mut filter),
                &mut subject_besthit,
            ),
            0
        );
        assert!(subject_besthit.is_none());
        assert_eq!(
            filter.subject_besthit_opts,
            Some(BlastHSPSubjectBestHitOptions { max_range_diff: 3 })
        );

        let mut culling = Some(blast_hspculling_options_new(7));
        assert_eq!(
            blast_hspfiltering_options_add_culling(
                Some(&mut filter),
                &mut culling,
                crate::util::EBlastStage::TracebackSearch,
            ),
            0
        );
        assert!(culling.is_none());
        assert_eq!(
            filter.culling_opts,
            Some(BlastHSPCullingOptions { max_hits: 7 })
        );
        assert_eq!(
            filter.culling_stage,
            crate::util::EBlastStage::TracebackSearch
        );
        assert_eq!(blast_hspfiltering_options_validate(&filter), 0);

        let mut slot = Some(filter);
        assert!(blast_hspfiltering_options_free(&mut slot).is_none());
        assert!(slot.is_none());
    }

    #[test]
    fn translated_hsp_filtering_options_validate_bounds_and_stage_conflicts() {
        let mut invalid_best_hit = blast_hspfiltering_options_new();
        invalid_best_hit.best_hit = Some(blast_hspbest_hit_options_new(0.0, 0.1));
        assert_eq!(blast_hspbest_hit_options_validate(&invalid_best_hit), -1);
        invalid_best_hit.best_hit = Some(blast_hspbest_hit_options_new(0.1, 0.5));
        assert_eq!(blast_hspbest_hit_options_validate(&invalid_best_hit), -1);

        let mut invalid_culling = blast_hspfiltering_options_new();
        invalid_culling.culling_opts = Some(blast_hspculling_options_new(-1));
        assert_eq!(blast_hspculling_options_validate(&invalid_culling), -1);

        let mut conflict = blast_hspfiltering_options_new();
        conflict.best_hit = Some(blast_hspbest_hit_options_new(0.1, 0.1));
        conflict.best_hit_stage = crate::util::EBlastStage::PrelimSearch;
        conflict.culling_opts = Some(blast_hspculling_options_new(3));
        conflict.culling_stage = crate::util::EBlastStage::Both;
        assert_eq!(blast_hspfiltering_options_validate(&conflict), 1);

        let mut missing_best_hit = None;
        assert_eq!(
            blast_hspfiltering_options_add_best_hit(
                Some(&mut conflict),
                &mut missing_best_hit,
                crate::util::EBlastStage::PrelimSearch,
            ),
            1
        );
        let mut missing_culling = None;
        assert_eq!(
            blast_hspfiltering_options_add_culling(
                None,
                &mut missing_culling,
                crate::util::EBlastStage::PrelimSearch,
            ),
            1
        );
        let mut missing_subject_besthit = None;
        assert_eq!(
            blast_hspfiltering_options_add_subject_best_hit(
                Some(&mut conflict),
                &mut missing_subject_besthit,
            ),
            1
        );
    }

    #[test]
    fn blast_hsp_collector_params_info_and_free() {
        let hit = crate::options::HitSavingOptions {
            hitlist_size: 25,
            ..Default::default()
        };
        let params = blast_hsp_collector_params_new(&hit, crate::program::BLASTP, 1, true, 7);
        assert_eq!(params.prelim_hitlist_size, 1050);
        assert_eq!(params.hsp_num_max, 7);

        let info = blast_hsp_collector_info_new(params.clone());
        assert_eq!(info.params.prelim_hitlist_size, 1050);
        let writer = s_blast_hspcollector_new(params.clone());
        assert!(!writer.rps_run);
        assert!(s_blast_hspcollector_free(Some(writer)).is_none());

        let mut slot = Some(params);
        assert!(blast_hsp_collector_params_free(&mut slot).is_none());
        assert!(slot.is_none());
    }

    #[test]
    fn blast_hsp_collector_run_splits_multiple_queries() {
        let hit = crate::options::HitSavingOptions {
            hitlist_size: 10,
            ..Default::default()
        };
        let params = blast_hsp_collector_params_new(&hit, crate::program::BLASTP, 0, true, 20);
        let mut data = BlastHSPCollectorData { params };
        let mut results = crate::hspstream::HspResults::new(2);
        let mut input = crate::hspstream::blast_hsp_list_new(20);
        input.oid = 42;
        input.hsps.push(hsp(100, 0, 0, 0));
        input.hsps.push(hsp(80, 0, 40, 20));
        input.hsps.push(hsp(90, 1, 20, 5));

        assert_eq!(
            s_blast_hspcollector_run(&mut data, &mut results, Some(input)),
            0
        );
        assert_eq!(results.hitlists.len(), 2);
        assert_eq!(results.hitlists[0].as_ref().unwrap().hsp_lists[0].oid, 42);
        assert_eq!(results.hitlists[0].as_ref().unwrap().hsp_lists.len(), 1);
        assert_eq!(
            results.hitlists[0].as_ref().unwrap().hsp_lists[0]
                .hsps
                .len(),
            2
        );
        assert_eq!(
            results.hitlists[1].as_ref().unwrap().hsp_lists[0].hsps[0].context,
            1
        );
        assert_eq!(s_blast_hspcollector_final(&mut data, &mut results), 0);
    }

    #[test]
    fn blast_hsp_collector_split_path_respects_hsp_cap() {
        let hit = crate::options::HitSavingOptions {
            hitlist_size: 10,
            ..Default::default()
        };
        let params = blast_hsp_collector_params_new(&hit, crate::program::BLASTP, 0, true, 1);
        let mut data = BlastHSPCollectorData { params };
        let mut results = crate::hspstream::HspResults::new(2);
        let mut input = crate::hspstream::blast_hsp_list_new(10);
        input.oid = 42;
        input.hsps.push(hsp(20, 0, 0, 0));
        input.hsps.push(hsp(90, 0, 20, 20));
        input.hsps.push(hsp(70, 1, 0, 0));

        assert_eq!(
            s_blast_hspcollector_run(&mut data, &mut results, Some(input)),
            0
        );
        let q0 = &results.hitlists[0].as_ref().unwrap().hsp_lists[0];
        assert_eq!(q0.hsps.len(), 1);
        assert_eq!(q0.hsps[0].score, 90);
        let q1 = &results.hitlists[1].as_ref().unwrap().hsp_lists[0];
        assert_eq!(q1.hsps.len(), 1);
        assert_eq!(q1.hsps[0].score, 70);
    }

    #[test]
    fn blast_hsp_collector_run_respects_prelim_hitlist_cap() {
        let hit = crate::options::HitSavingOptions {
            hitlist_size: 1,
            ..Default::default()
        };
        let params = blast_hsp_collector_params_new(&hit, crate::program::BLASTP, 0, false, 0);
        assert_eq!(params.prelim_hitlist_size, 1);
        assert_eq!(params.hsp_num_max, i32::MAX);
        let mut data = BlastHSPCollectorData { params };
        let mut results = crate::hspstream::HspResults::new(1);

        for &(oid, score) in &[(1, 30), (2, 90), (3, 40)] {
            let mut list = crate::hspstream::blast_hsp_list_new(0);
            list.oid = oid;
            list.hsps.push(hsp(score, 0, 0, 0));
            assert_eq!(
                s_blast_hspcollector_run(&mut data, &mut results, Some(list)),
                0
            );
        }

        let hitlist = results.hitlists[0].as_ref().unwrap();
        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(hitlist.hsp_lists[0].oid, 2);
        assert_eq!(hitlist.hsp_lists[0].hsps[0].score, 90);
    }

    #[test]
    fn blast_hsp_collector_final_sorts_heap_hitlists_by_evalue() {
        let hit = crate::options::HitSavingOptions {
            hitlist_size: 3,
            ..Default::default()
        };
        let params = blast_hsp_collector_params_new(&hit, crate::program::BLASTP, 0, false, 0);
        let mut data = BlastHSPCollectorData { params };
        let mut results = crate::hspstream::HspResults::new(1);

        for &(oid, score) in &[(1, 30), (2, 90), (3, 40), (4, 80)] {
            let mut list = crate::hspstream::blast_hsp_list_new(0);
            list.oid = oid;
            list.hsps.push(hsp(score, 0, 0, 0));
            assert_eq!(
                s_blast_hspcollector_run(&mut data, &mut results, Some(list)),
                0
            );
        }

        let hitlist = results.hitlists[0].as_ref().unwrap();
        assert_eq!(hitlist.hsp_lists.len(), 3);
        assert_eq!(hitlist.hsp_lists[0].oid, 3);

        assert_eq!(s_blast_hspcollector_final(&mut data, &mut results), 0);
        let ordered: Vec<i32> = results.hitlists[0]
            .as_ref()
            .unwrap()
            .hsp_lists
            .iter()
            .map(|list| list.oid)
            .collect();
        assert_eq!(ordered, vec![2, 4, 3]);
    }

    #[test]
    fn blast_hsp_collector_run_rps_groups_by_context_as_oid() {
        let hit = crate::options::HitSavingOptions {
            hitlist_size: 10,
            ..Default::default()
        };
        let params = blast_hsp_collector_params_new(&hit, crate::program::RPS_BLAST, 0, true, 20);
        let mut writer = s_blast_hspcollector_new(params);
        assert!(writer.rps_run);
        let mut results = crate::hspstream::HspResults::new(1);
        let mut input = crate::hspstream::blast_hsp_list_new(20);
        input.hsps.push(hsp(50, 9, 0, 30));
        input.hsps.push(hsp(80, 7, 0, 10));
        input.hsps.push(hsp(70, 9, 0, 20));

        assert_eq!(
            s_blast_hspcollector_run_rps(&mut writer.data, &mut results, 0, input),
            0
        );
        let hitlist = results.hitlists[0].as_ref().unwrap();
        assert_eq!(hitlist.hsp_lists.len(), 2);
        assert_eq!(hitlist.hsp_lists[0].oid, 7);
        assert_eq!(hitlist.hsp_lists[0].hsps[0].context, 0);
        assert_eq!(hitlist.hsp_lists[1].oid, 9);
        assert_eq!(hitlist.hsp_lists[1].hsps[0].score, 70);
    }

    #[test]
    fn blast_hsp_collector_rps_path_respects_hsp_cap() {
        let hit = crate::options::HitSavingOptions {
            hitlist_size: 10,
            ..Default::default()
        };
        let params = blast_hsp_collector_params_new(&hit, crate::program::RPS_BLAST, 0, true, 1);
        let mut writer = s_blast_hspcollector_new(params);
        let mut results = crate::hspstream::HspResults::new(1);
        let mut input = crate::hspstream::blast_hsp_list_new(10);
        input.hsps.push(hsp(50, 9, 0, 30));
        input.hsps.push(hsp(80, 9, 0, 10));

        assert_eq!(
            s_blast_hspcollector_run_rps(&mut writer.data, &mut results, 0, input),
            0
        );
        let hitlist = results.hitlists[0].as_ref().unwrap();
        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(hitlist.hsp_lists[0].hsps.len(), 1);
        assert_eq!(hitlist.hsp_lists[0].hsps[0].score, 80);
    }

    #[test]
    fn ctree_new_creates_root_with_full_range() {
        let tree = s_ctree_new(1000);
        assert_eq!(tree.begin, 0);
        assert_eq!(tree.end, 1000);
        assert!(tree.left.is_none());
        assert!(tree.right.is_none());
    }

    #[test]
    fn save_hsp_inserts_until_dominated() {
        let mut tree = *s_ctree_new(1000);
        // Insert a few non-dominating HSPs.
        for (score, b, e, sid) in [(50, 0, 100, 1), (40, 200, 300, 2), (60, 500, 600, 3)] {
            let mut node = *mk(score, b, e, sid);
            node.merit = 5;
            assert!(s_save_hsp(&mut tree, &mut node));
        }
        // Now insert one whose merit is too low to survive.
        let mut bad = *mk(20, 0, 100, 9);
        bad.merit = 0; // already at zero → first dominator drops it
                       // Note: merit=0 plus s_full_pass decrement → drops on first match.
        let _ = s_save_hsp(&mut tree, &mut bad);
    }

    #[test]
    fn rip_hsp_off_ctree_collects_all_nodes() {
        let mut tree = *s_ctree_new(1000);
        for (score, b, e, sid) in [(50, 0, 100, 1), (40, 200, 300, 2), (60, 500, 600, 3)] {
            let mut node = *mk(score, b, e, sid);
            node.merit = 5;
            s_save_hsp(&mut tree, &mut node);
        }
        let ripped = s_rip_hsp_off_ctree(Some(Box::new(tree)));
        // Count survivors.
        let mut count = 0;
        let mut cur = ripped.as_deref();
        while let Some(n) = cur {
            count += 1;
            cur = n.next.as_deref();
        }
        assert_eq!(count, 3);
    }

    #[test]
    fn fork_children_splits_by_midpoint() {
        let mut node = CTreeNode {
            begin: 0,
            end: 100,
            left: None,
            right: None,
            hsp_list: None,
        };
        // Three HSPs: one strictly left, one strictly right, one
        // spanning the midpoint.
        s_add_hsp_to_list(&mut node.hsp_list, mk(50, 5, 25, 1)); // strictly left
        s_add_hsp_to_list(&mut node.hsp_list, mk(60, 60, 90, 2)); // strictly right
        s_add_hsp_to_list(&mut node.hsp_list, mk(70, 30, 70, 3)); // spans midpt 50
        s_fork_children(&mut node);
        // Spanning HSP stays.
        assert!(node.hsp_list.is_some());
        let stays = node.hsp_list.as_ref().unwrap();
        assert_eq!(stays.sid, 3);
        // Left and right children populated.
        let left = node.left.as_ref().expect("left");
        assert_eq!(left.hsp_list.as_ref().unwrap().sid, 1);
        let right = node.right.as_ref().expect("right");
        assert_eq!(right.hsp_list.as_ref().unwrap().sid, 2);
    }

    #[test]
    fn ctree_node_new_child_splits_interval() {
        let mut parent = CTreeNode {
            begin: 0,
            end: 100,
            left: None,
            right: None,
            hsp_list: None,
        };
        let _ = &mut parent;
        let left = s_ctree_node_new(Some(&parent), CTreeChild::Left);
        let right = s_ctree_node_new(Some(&parent), CTreeChild::Right);
        assert_eq!(left.begin, 0);
        assert_eq!(left.end, 50);
        assert_eq!(right.begin, 50);
        assert_eq!(right.end, 100);
    }

    #[test]
    fn ctree_node_new_child_root_is_uninitialized() {
        let node = s_ctree_node_new(None, CTreeChild::Left);
        assert_eq!(node.begin, 0);
        assert_eq!(node.end, 0);
    }

    #[test]
    fn ctree_memory_wrappers_drop_owned_values() {
        let node = s_get_node();
        assert_eq!(node.begin, 0);
        assert_eq!(node.end, 0);
        assert!(s_ret_node(Some(node)).is_none());

        let leaf = s_ctree_node_new(None, CTreeChild::Left);
        assert!(s_ctree_node_free(Some(leaf)).is_none());

        let mut tree = s_ctree_new(42);
        tree.left = Some(s_ctree_node_new(Some(&tree), CTreeChild::Left));
        tree.right = Some(s_ctree_node_new(Some(&tree), CTreeChild::Right));
        s_debug(Some(&tree));
        assert!(s_ctree_free(Some(tree)).is_none());
        assert!(s_ctree_free(None).is_none());
    }

    #[test]
    fn linked_hsp_copy_detaches_next_pointer() {
        let mut a = LinkedHsp {
            hsp: Hsp {
                score: 100,
                num_ident: 0,
                bit_score: 0.0,
                evalue: 0.0,
                query_offset: 0,
                query_end: 100,
                query_gapped_start: 0,
                subject_offset: 0,
                subject_end: 100,
                subject_gapped_start: 0,
                context: 0,
                query_frame: 0,
                subject_frame: 0,
                num_gaps: 0,
                comp_adjustment_method: 0,
                edit_script: None,
                pat_info: None,
                map_info: None,
            },
            cid: 0,
            sid: 5,
            begin: 0,
            end: 100,
            merit: 3,
            next: None,
        };
        a.next = Some(Box::new(s_hsp_copy(&a)));
        let copy = s_hsp_copy(&a);
        assert_eq!(copy.sid, 5);
        assert_eq!(copy.merit, 3);
        // Detached: `next` should be None even though `a` had a next.
        assert!(copy.next.is_none());
        assert!(s_hsp_free(Some(Box::new(copy))).is_none());
    }
}
