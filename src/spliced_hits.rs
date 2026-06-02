//! Compatibility ports for `spliced_hits.c` and the small chain-list helpers
//! at the front of `hspfilter_mapper.c`.

use crate::filter::get_cutoff_score;
use crate::gapinfo::{
    GapAlignOpType, JumperEdit, MAPPER_ADAPTER, MAPPER_EXON, MAPPER_POLY_A, MAPPER_SPLICE_SIGNAL,
};
use crate::hspstream::{Hsp, HspList};
use crate::options::{HitSavingOptions, ScoringOptions};
use crate::program::{ProgramType, MAPPING};
use crate::queryinfo::{QueryInfo, E_FIRST_SEGMENT};
use crate::util::NUM_STRANDS;

#[derive(Debug, Clone)]
pub struct HSPContainer {
    pub hsp: Hsp,
    pub next: Option<Box<HSPContainer>>,
}

/// Port of NCBI `HSPContainerNew` (`spliced_hits.c:37`).
pub fn hsp_container_new(hsp: &mut Option<Hsp>) -> Option<Box<HSPContainer>> {
    hsp.take()
        .map(|hsp| Box::new(HSPContainer { hsp, next: None }))
}

/// Rust ownership equivalent of NCBI `HSPContainerFree` (`spliced_hits.c:50`).
pub fn hsp_container_free(_: Option<Box<HSPContainer>>) -> Option<Box<HSPContainer>> {
    None
}

/// Port of NCBI `HSPContainerDup` (`spliced_hits.c:66`).
pub fn hsp_container_dup(input: Option<&HSPContainer>) -> Option<Box<HSPContainer>> {
    let input = input?;
    Some(Box::new(HSPContainer {
        hsp: input.hsp.clone(),
        next: hsp_container_dup(input.next.as_deref()),
    }))
}

#[derive(Debug, Clone)]
pub struct HSPChain {
    pub context: i32,
    pub oid: i32,
    pub score: i32,
    pub hsps: Option<Box<HSPContainer>>,
    pub count: i32,
    pub pair: Option<usize>,
    pub pair_score: Option<i32>,
    pub pair_conf: u8,
    pub adapter: i32,
    pub poly_a: i32,
    pub next: Option<Box<HSPChain>>,
}

pub const PAIR_CONVERGENT: u8 = 0;
pub const PAIR_DIVERGENT: u8 = 1;
pub const PAIR_PARALLEL: u8 = 2;
pub const PAIR_NONE: u8 = 3;
const MAGICBLAST_MAX_INSERT_SIZE_SPLICED: i32 = 500_000;
const MAGICBLAST_MAX_INSERT_SIZE_NONSPLICED: i32 = 1_000;

/// Rust ownership equivalent of NCBI `HSPChainFree` (`spliced_hits.c:109`).
pub fn hsp_chain_free(chain_list: Option<Box<HSPChain>>) -> Option<Box<HSPChain>> {
    let mut chain = chain_list;
    while let Some(mut current) = chain {
        let next = current.next.take();
        current.pair = None;
        current.pair_score = None;
        current.hsps = hsp_container_free(current.hsps.take());
        chain = next;
    }
    None
}

/// Port of NCBI `HSPChainNew` (`spliced_hits.c:126`).
pub fn hsp_chain_new(context: i32) -> Box<HSPChain> {
    Box::new(HSPChain {
        context,
        oid: 0,
        score: 0,
        hsps: None,
        count: 0,
        pair: None,
        pair_score: None,
        pair_conf: PAIR_NONE,
        adapter: -1,
        poly_a: 0,
        next: None,
    })
}

/// Port of NCBI `CloneChain` (`spliced_hits.c:139`).
pub fn clone_chain(chain: Option<&HSPChain>) -> Option<Box<HSPChain>> {
    let chain = chain?;
    Some(Box::new(HSPChain {
        context: chain.context,
        oid: chain.oid,
        score: chain.score,
        hsps: hsp_container_dup(chain.hsps.as_deref()),
        count: chain.count,
        pair: None,
        pair_score: None,
        pair_conf: chain.pair_conf,
        adapter: chain.adapter,
        poly_a: chain.poly_a,
        next: None,
    }))
}

fn append_cloned_chain(list: &mut Option<Box<HSPChain>>, chain: &HSPChain) {
    let Some(cloned) = clone_chain(Some(chain)) else {
        return;
    };
    let mut cursor = list;
    loop {
        match cursor {
            Some(node) => cursor = &mut node.next,
            None => {
                *cursor = Some(cloned);
                break;
            }
        }
    }
}

/// Port of NCBI `FindPartialyCoveredQueries` (`hspfilter_mapper.c:4568`).
///
/// The upstream function name contains the same misspelling. It scans the
/// mapper's saved chains and clones chains for `oid` whose first HSP starts
/// after `word_size` or whose last HSP leaves more than `word_size` residues
/// uncovered at the query tail.
pub fn find_partialy_covered_queries(
    data: Option<&BlastHSPMapperData>,
    oid: i32,
    word_size: i32,
) -> Option<Box<HSPChain>> {
    let data = data?;
    let query_info = data.query_info.as_ref()?;
    let mut retval: Option<Box<HSPChain>> = None;

    for i in 0..query_info.num_queries.max(0) as usize {
        let mut chain = data.saved_chains.get(i).and_then(|chain| chain.as_deref());
        while let Some(current) = chain {
            if current.oid == oid && current.score >= 30 {
                if let Some(first) = current.hsps.as_deref() {
                    let mut should_clone = first.hsp.query_offset > word_size;
                    if !should_clone {
                        let mut last = first;
                        while let Some(next) = last.next.as_deref() {
                            last = next;
                        }
                        let context = last.hsp.context.max(0) as usize;
                        if let Some(ctx) = query_info.contexts.get(context) {
                            should_clone =
                                ctx.query_length.saturating_sub(last.hsp.query_end) > word_size;
                        }
                    }
                    if should_clone {
                        append_cloned_chain(&mut retval, current);
                    }
                }
            }
            chain = current.next.as_deref();
        }
    }

    retval
}

#[derive(Debug, Clone, Default)]
pub struct BlastMappingResults {
    pub num_queries: i32,
    pub chain_array: Vec<Option<Box<HSPChain>>>,
}

#[derive(Debug, Clone, Default)]
pub struct BlastHSPMapperData {
    pub params: Option<BlastHSPMapperParams>,
    pub query: Option<Vec<u8>>,
    pub query_info: Option<QueryInfo>,
    pub saved_chains: Vec<Option<Box<HSPChain>>>,
}

#[derive(Debug, Clone)]
pub struct BlastHSPMapperParams {
    pub program: ProgramType,
    pub hitlist_size: i32,
    pub paired: bool,
    pub splice: bool,
    pub longest_intron: i32,
    pub cutoff_score: i32,
    pub cutoff_score_fun: [i32; 2],
    pub cutoff_edit_dist: i32,
    pub scoring_options: ScoringOptions,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BlastHSPMapperCallback {
    PairedInit,
    Final,
    Free,
    SplicedPairedRun,
    PairedNew,
}

#[derive(Debug, Clone)]
pub struct BlastHSPMapperWriter {
    pub init_fn: BlastHSPMapperCallback,
    pub final_fn: BlastHSPMapperCallback,
    pub free_fn: BlastHSPMapperCallback,
    pub run_fn: BlastHSPMapperCallback,
    pub data: BlastHSPMapperData,
}

#[derive(Debug, Clone)]
pub struct BlastHSPMapperInfo {
    pub new_fn: BlastHSPMapperCallback,
    pub params: BlastHSPMapperParams,
}

pub const MAX_NUM_HSP_PATHS: usize = 40;

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct HSPNode {
    pub hsp_index: Option<usize>,
    pub best_score: i32,
    pub path_next: Option<usize>,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct HSPPath {
    pub start: Vec<Option<usize>>,
    pub num_paths: i32,
    pub score: i32,
}

/// Port of NCBI `Blast_MappingResultsNew` (`spliced_hits.c:166`).
pub fn blast_mapping_results_new() -> BlastMappingResults {
    BlastMappingResults::default()
}

/// Rust ownership equivalent of NCBI `Blast_MappingResultsFree`.
pub fn blast_mapping_results_free(_: Option<BlastMappingResults>) -> Option<BlastMappingResults> {
    None
}

/// Port of NCBI `BlastHSPMapperParamsNew` (`hspfilter_mapper.c:4897`).
///
/// The Rust `HitSavingOptions` type carries the fields used elsewhere in this
/// port. Mapper-only C fields that are not represented there are initialized to
/// the same neutral defaults used by the mapper when those features are absent.
pub fn blast_hsp_mapper_params_new(
    hit_options: Option<&HitSavingOptions>,
    scoring_options: Option<&ScoringOptions>,
) -> Option<BlastHSPMapperParams> {
    let hit_options = hit_options?;
    let scoring_options = scoring_options?;
    let mut mapper_scoring = scoring_options.clone();
    mapper_scoring.gap_open = -mapper_scoring.gap_open;
    mapper_scoring.gap_extend = -mapper_scoring.gap_extend;

    Some(BlastHSPMapperParams {
        program: MAPPING,
        hitlist_size: hit_options.hitlist_size.max(10),
        paired: true,
        splice: false,
        longest_intron: 0,
        cutoff_score: hit_options.cutoff_score,
        cutoff_score_fun: [0, 0],
        cutoff_edit_dist: -1,
        scoring_options: mapper_scoring,
    })
}

/// Rust ownership equivalent of NCBI `BlastHSPMapperParamsFree`.
pub fn blast_hsp_mapper_params_free(
    _: Option<BlastHSPMapperParams>,
) -> Option<BlastHSPMapperParams> {
    None
}

/// Port of NCBI `BlastHSPMapperInfoNew` (`hspfilter_mapper.c:4933`).
pub fn blast_hsp_mapper_info_new(params: BlastHSPMapperParams) -> BlastHSPMapperInfo {
    BlastHSPMapperInfo {
        new_fn: BlastHSPMapperCallback::PairedNew,
        params,
    }
}

pub fn s_get_chain_score(chain: Option<&HSPChain>) -> i32 {
    chain.map_or(0, |chain| chain.score + chain.pair_score.unwrap_or(0))
}

fn first_hsp(chain: &HSPChain) -> Option<&Hsp> {
    chain.hsps.as_deref().map(|container| &container.hsp)
}

fn duplicate_chain(a: &HSPChain, b: &HSPChain) -> bool {
    let (Some(a_hsp), Some(b_hsp)) = (first_hsp(a), first_hsp(b)) else {
        return false;
    };
    a.oid == b.oid
        && a.score == b.score
        && a_hsp.query_frame == b_hsp.query_frame
        && a_hsp.subject_offset == b_hsp.subject_offset
}

fn drain_chain_list(mut list: Option<Box<HSPChain>>) -> Vec<Box<HSPChain>> {
    let mut out = Vec::new();
    while let Some(mut chain) = list {
        list = chain.next.take();
        out.push(chain);
    }
    out
}

fn build_chain_list(mut chains: Vec<Box<HSPChain>>) -> Option<Box<HSPChain>> {
    let mut list = None;
    while let Some(mut chain) = chains.pop() {
        chain.next = list;
        list = Some(chain);
    }
    list
}

fn drain_hsp_container_list(mut list: Option<Box<HSPContainer>>) -> Vec<Box<HSPContainer>> {
    let mut out = Vec::new();
    while let Some(mut container) = list {
        list = container.next.take();
        out.push(container);
    }
    out
}

fn build_hsp_container_list(mut hsps: Vec<Box<HSPContainer>>) -> Option<Box<HSPContainer>> {
    let mut list = None;
    while let Some(mut container) = hsps.pop() {
        container.next = list;
        list = Some(container);
    }
    list
}

/// Score-sorted chain insertion used when no cutoff filtering is requested.
pub fn s_hsp_chain_list_insert_one_by_score(
    list: &mut Option<Box<HSPChain>>,
    chain: Option<Box<HSPChain>>,
    check_for_duplicates: bool,
) -> i32 {
    let Some(chain) = chain else {
        return -1;
    };

    let mut chains = drain_chain_list(list.take());
    if check_for_duplicates
        && chains
            .iter()
            .any(|existing| duplicate_chain(existing, &chain))
    {
        *list = build_chain_list(chains);
        return 0;
    }

    chains.push(chain);
    chains.sort_by(|a, b| b.score.cmp(&a.score));
    *list = build_chain_list(chains);
    0
}

/// Port of NCBI `s_TestCutoffs` (`hspfilter_mapper.c:124`).
pub fn s_test_cutoffs(chain: Option<&HSPChain>, cutoff_score: i32, cutoff_edit_dist: i32) -> bool {
    let Some(chain) = chain else {
        return false;
    };
    if chain.score < cutoff_score {
        return false;
    }
    if cutoff_edit_dist < 0 {
        return true;
    }

    let mut align_len: i32 = 0;
    let mut num_identical: i32 = 0;
    let mut container = chain.hsps.as_deref();
    while let Some(hsp_container) = container {
        let hsp = &hsp_container.hsp;
        let query_span = hsp.query_end.saturating_sub(hsp.query_offset).max(0);
        let subject_span = hsp.subject_end.saturating_sub(hsp.subject_offset).max(0);
        align_len = align_len.saturating_add(query_span.max(subject_span));
        num_identical = num_identical.saturating_add(hsp.num_ident);
        container = hsp_container.next.as_deref();
    }
    align_len.saturating_sub(num_identical) <= cutoff_edit_dist
}

/// Port of NCBI `s_HSPChainListInsertOne` (`hspfilter_mapper.c:165`).
pub fn s_hsp_chain_list_insert_one(
    list: &mut Option<Box<HSPChain>>,
    chain: Option<Box<HSPChain>>,
    check_for_duplicates: bool,
) -> i32 {
    let Some(chain) = chain else {
        return -1;
    };

    if list.is_none() {
        *list = Some(chain);
        return 0;
    }

    let list_score = s_get_chain_score(list.as_deref());
    let chain_score = s_get_chain_score(Some(&chain));
    if list_score > chain_score {
        return 0;
    }
    if list_score < chain_score {
        *list = Some(chain);
        return 0;
    }

    let mut chains = drain_chain_list(list.take());
    if check_for_duplicates
        && chains
            .iter()
            .any(|existing| duplicate_chain(existing, &chain))
    {
        *list = build_chain_list(chains);
        return 0;
    }

    chains.push(chain);
    *list = build_chain_list(chains);
    0
}

/// Port of NCBI `HSPChainListInsert` (`hspfilter_mapper.c:238`).
pub fn hsp_chain_list_insert(
    list: &mut Option<Box<HSPChain>>,
    chain: &mut Option<Box<HSPChain>>,
    cutoff_score: i32,
    cutoff_edit_dist: i32,
    check_for_duplicates: bool,
) -> i32 {
    let mut status = 0;
    for ch in drain_chain_list(chain.take()) {
        if cutoff_score <= 0 {
            status = s_hsp_chain_list_insert_one_by_score(list, Some(ch), check_for_duplicates);
        } else if s_test_cutoffs(Some(&ch), cutoff_score, cutoff_edit_dist) {
            status = s_hsp_chain_list_insert_one(list, Some(ch), check_for_duplicates);
        }
        if status != 0 {
            return status;
        }
    }
    0
}

/// Port of NCBI `HSPChainListTrim` (`hspfilter_mapper.c:277`).
pub fn hsp_chain_list_trim(list: &mut Option<Box<HSPChain>>, margin: i32) -> i32 {
    let Some(head) = list.as_ref() else {
        return -1;
    };
    let best_score = head.score;
    let mut chains = drain_chain_list(list.take());
    let keep_len = chains
        .iter()
        .take_while(|chain| best_score - chain.score <= margin)
        .count();
    chains.truncate(keep_len.max(1));
    *list = build_chain_list(chains);
    0
}

/// Port of NCBI `s_TestChains` (`hspfilter_mapper.c:301`).
pub fn s_test_chains(chain: Option<&HSPChain>) -> bool {
    let mut chain = chain;
    while let Some(current) = chain {
        let Some(mut container) = current.hsps.as_deref() else {
            return false;
        };
        if container.hsp.context != current.context {
            return false;
        }
        while let Some(next) = container.next.as_deref() {
            if next.hsp.context != current.context {
                return false;
            }
            if container.hsp.query_offset >= next.hsp.query_offset
                || container.hsp.subject_offset >= next.hsp.subject_offset
            {
                return false;
            }
            container = next;
        }
        chain = current.next.as_deref();
    }
    true
}

/// Port of NCBI debug helper `s_TestChainsSorted` (`hspfilter_mapper.c:366`).
pub fn s_test_chains_sorted(chain: Option<&HSPChain>) -> bool {
    if !s_test_chains(chain) {
        return false;
    }
    let mut prev = chain;
    let mut current = prev.and_then(|chain| chain.next.as_deref());
    while let (Some(p), Some(c)) = (prev, current) {
        if p.score < c.score {
            return false;
        }
        prev = current;
        current = c.next.as_deref();
    }
    true
}

/// Port of NCBI `s_BlastHSPMapperPairedInit` (`hspfilter_mapper.c:388`).
pub fn s_blast_hsp_mapper_paired_init(data: &mut BlastHSPMapperData, num_queries: i32) -> i32 {
    data.saved_chains.clear();
    data.saved_chains
        .extend((0..num_queries.max(0)).map(|_| None));
    0
}

/// Port of NCBI `s_BlastHSPMapperPairedNew` (`hspfilter_mapper.c:4871`).
pub fn s_blast_hsp_mapper_paired_new(
    params: Option<BlastHSPMapperParams>,
    query_info: Option<&QueryInfo>,
    query: Option<&[u8]>,
) -> BlastHSPMapperWriter {
    BlastHSPMapperWriter {
        init_fn: BlastHSPMapperCallback::PairedInit,
        final_fn: BlastHSPMapperCallback::Final,
        free_fn: BlastHSPMapperCallback::Free,
        run_fn: BlastHSPMapperCallback::SplicedPairedRun,
        data: BlastHSPMapperData {
            params,
            query: query.map(|query| query.to_vec()),
            query_info: query_info.cloned(),
            saved_chains: Vec::new(),
        },
    }
}

/// Rust ownership equivalent of NCBI `s_BlastHSPMapperFree`.
pub fn s_blast_hsp_mapper_free(_: Option<BlastHSPMapperWriter>) -> Option<BlastHSPMapperWriter> {
    None
}

/// Port of NCBI `s_GetOverlapCost` (`hspfilter_mapper.c:402`) for the Rust
/// `Hsp` shape.
pub fn s_get_overlap_cost(a: &Hsp, b: &Hsp, edit_penalty: i32) -> i32 {
    let a_edit_positions = mapper_edit_positions(a);
    let b_edit_positions = mapper_edit_positions(b);

    s_get_overlap_cost_with_edits(a, b, edit_penalty, &a_edit_positions, &b_edit_positions)
}

fn mapper_edit_positions(hsp: &Hsp) -> Vec<i32> {
    hsp.map_info
        .as_ref()
        .and_then(|info| info.edits.as_ref())
        .map(|edits| edits.edits.iter().map(|edit| edit.query_pos).collect())
        .unwrap_or_default()
}

pub fn s_get_overlap_cost_with_edits(
    a: &Hsp,
    b: &Hsp,
    edit_penalty: i32,
    a_edit_positions: &[i32],
    b_edit_positions: &[i32],
) -> i32 {
    if (a.query_offset <= b.query_offset && a.query_end >= b.query_end)
        || (b.query_offset <= a.query_offset && b.query_end >= a.query_end)
    {
        return a.score.min(b.score);
    }
    if a.query_end <= b.query_offset || b.query_end <= a.query_offset {
        return 0;
    }

    let (first, second, first_edits, second_edits) = if a.query_offset < b.query_offset {
        (a, b, a_edit_positions, b_edit_positions)
    } else {
        (b, a, b_edit_positions, a_edit_positions)
    };

    let mut first_overlap = first.query_end.saturating_sub(second.query_offset);
    let mut second_overlap = first_overlap;
    for &query_pos in first_edits {
        if query_pos >= second.query_offset {
            first_overlap = first_overlap.saturating_sub(edit_penalty);
        }
    }
    for &query_pos in second_edits {
        if query_pos < first.query_end {
            second_overlap = second_overlap.saturating_sub(edit_penalty);
        }
    }
    first_overlap.min(second_overlap)
}

fn s_compute_gap_score(length: i32, open_score: i32, extend_score: i32, seq_error: i32) -> i32 {
    if length < 4 {
        length * seq_error
    } else {
        open_score + length.min(4) * extend_score
    }
}

fn s_compute_alignment_score(hsp: &Hsp, score_options: &ScoringOptions, query_len: i32) -> i32 {
    let _ = query_len;
    if let Some(script) = hsp.edit_script.as_ref() {
        let mut score = 0;
        let edit_positions = hsp
            .map_info
            .as_ref()
            .and_then(|info| info.edits.as_ref())
            .map(|edits| edits.edits.as_slice())
            .unwrap_or(&[]);
        let mut query_pos = hsp.query_offset;
        for (op, count) in script.iter() {
            match op {
                GapAlignOpType::Sub => {
                    let query_end = query_pos.saturating_add(count);
                    let mismatches = edit_positions
                        .iter()
                        .filter(|edit| {
                            edit.query_pos >= query_pos
                                && edit.query_pos < query_end
                                && edit.query_base != 15
                                && edit.subject_base != 15
                                && edit.query_base != edit.subject_base
                        })
                        .count() as i32;
                    score += (count - mismatches) * score_options.reward
                        + mismatches * score_options.penalty;
                    query_pos = query_end;
                }
                GapAlignOpType::Ins | GapAlignOpType::Ins1 | GapAlignOpType::Ins2 => {
                    score -= score_options.gap_open + count * score_options.gap_extend;
                    query_pos = query_pos.saturating_add(count);
                }
                GapAlignOpType::Del | GapAlignOpType::Del1 | GapAlignOpType::Del2 => {
                    score -= score_options.gap_open + count * score_options.gap_extend;
                }
                GapAlignOpType::Decline => {}
            }
        }
        score
    } else {
        hsp.score
    }
}

/// Port of NCBI `s_ComputeChainScore` (`hspfilter_mapper.c:482`).
pub fn s_compute_chain_score(
    chain: Option<&mut HSPChain>,
    score_options: Option<&ScoringOptions>,
    query_len: i32,
    comp_hsp_score: bool,
) -> i32 {
    let Some(chain) = chain else {
        return -1;
    };
    let Some(score_options) = score_options else {
        return -1;
    };
    let Some(first) = chain.hsps.as_deref_mut() else {
        chain.score = 0;
        return 0;
    };

    if comp_hsp_score {
        first.hsp.score = s_compute_alignment_score(&first.hsp, score_options, query_len);
    }
    let mut score = first.hsp.score;
    let mut previous_query_end = first.hsp.query_end;
    let mut previous_subject_end = first.hsp.subject_end;
    let mut previous_right_splice = first
        .hsp
        .map_info
        .as_ref()
        .is_some_and(|info| info.right_edge & MAPPER_SPLICE_SIGNAL != 0);

    let mut current = first.next.as_deref_mut();
    while let Some(container) = current {
        if comp_hsp_score {
            container.hsp.score =
                s_compute_alignment_score(&container.hsp, score_options, query_len);
        }
        score += container.hsp.score;

        let current_left_splice = container
            .hsp
            .map_info
            .as_ref()
            .is_some_and(|info| info.left_edge & MAPPER_SPLICE_SIGNAL != 0);
        if !current_left_splice || !previous_right_splice {
            let query_gap = container
                .hsp
                .query_offset
                .saturating_sub(previous_query_end)
                .max(0);
            let subject_gap = container
                .hsp
                .subject_offset
                .saturating_sub(previous_subject_end)
                .max(0);
            score += s_compute_gap_score(query_gap, -12, -1, -4);
            score += s_compute_gap_score(subject_gap, -12, -1, -4);
        }

        previous_query_end = container.hsp.query_end;
        previous_subject_end = container.hsp.subject_end;
        previous_right_splice = container
            .hsp
            .map_info
            .as_ref()
            .is_some_and(|info| info.right_edge & MAPPER_SPLICE_SIGNAL != 0);
        current = container.next.as_deref_mut();
    }

    chain.score = score;
    score
}

/// Port of NCBI `s_FindFragmentStart` (`hspfilter_mapper.c:409`).
pub fn s_find_fragment_start(chain: Option<&HSPChain>) -> i32 {
    let Some(chain) = chain else {
        return -1;
    };
    let Some(first) = chain.hsps.as_deref() else {
        return -1;
    };
    if first.hsp.query_frame > 0 {
        return first.hsp.subject_offset;
    }

    let mut last = first;
    while let Some(next) = last.next.as_deref() {
        last = next;
    }
    last.hsp.subject_end.saturating_sub(1)
}

/// Port of NCBI debug helper `s_TestHSPRanges` (`hspfilter_mapper.c:577`).
pub fn s_test_hsp_ranges(hsp: &Hsp) -> bool {
    if hsp.query_offset < 0 || hsp.subject_offset < 0 {
        return false;
    }
    if hsp.query_end < hsp.query_offset || hsp.subject_end < hsp.subject_offset {
        return false;
    }
    if hsp.query_gapped_start > hsp.query_offset || hsp.subject_gapped_start > hsp.subject_offset {
        return false;
    }

    let query_span = i64::from(hsp.query_end.saturating_sub(hsp.query_offset));
    let subject_span = i64::from(hsp.subject_end.saturating_sub(hsp.subject_offset));
    if query_span == 0 && subject_span == 0 {
        return false;
    }

    if let Some(script) = hsp.edit_script.as_ref() {
        let mut query_used = 0i64;
        let mut subject_used = 0i64;
        for (op, count) in script.iter() {
            if count < 0 {
                return false;
            }
            let count = i64::from(count);
            match op {
                GapAlignOpType::Sub => {
                    query_used += count;
                    subject_used += count;
                }
                GapAlignOpType::Ins | GapAlignOpType::Ins1 | GapAlignOpType::Ins2 => {
                    query_used += count;
                }
                GapAlignOpType::Del | GapAlignOpType::Del1 | GapAlignOpType::Del2 => {
                    subject_used += count;
                }
                GapAlignOpType::Decline => {}
            }
        }
        if query_used > 0 && query_used != query_span {
            return false;
        }
        if subject_used > 0 && subject_used != subject_span {
            return false;
        }
    }

    true
}

fn mapper_trim_op_consumes_query(op: GapAlignOpType) -> bool {
    matches!(
        op,
        GapAlignOpType::Sub | GapAlignOpType::Ins | GapAlignOpType::Ins1 | GapAlignOpType::Ins2
    )
}

fn mapper_trim_op_consumes_subject(op: GapAlignOpType) -> bool {
    matches!(
        op,
        GapAlignOpType::Sub | GapAlignOpType::Del | GapAlignOpType::Del1 | GapAlignOpType::Del2
    )
}

fn mapper_trim_op_deltas(op: GapAlignOpType, count: i32) -> (i32, i32) {
    let query = if mapper_trim_op_consumes_query(op) {
        count
    } else {
        0
    };
    let subject = if mapper_trim_op_consumes_subject(op) {
        count
    } else {
        0
    };
    (query, subject)
}

fn mapper_fallback_score_from_script(
    hsp: &mut Hsp,
    mismatch_score: i32,
    gap_open_score: i32,
    gap_extend_score: i32,
) -> i32 {
    let Some(script) = hsp.edit_script.as_ref() else {
        hsp.num_ident = hsp
            .query_end
            .saturating_sub(hsp.query_offset)
            .min(hsp.subject_end.saturating_sub(hsp.subject_offset))
            .max(0);
        return hsp.num_ident;
    };

    let mut score: i32 = 0;
    let mut num_ident: i32 = 0;
    for (op, count) in script.iter() {
        match op {
            GapAlignOpType::Sub => {
                score = score.saturating_add(count);
                num_ident = num_ident.saturating_add(count);
            }
            GapAlignOpType::Ins
            | GapAlignOpType::Ins1
            | GapAlignOpType::Ins2
            | GapAlignOpType::Del
            | GapAlignOpType::Del1
            | GapAlignOpType::Del2 => {
                score = score.saturating_add(s_compute_gap_score(
                    count,
                    gap_open_score,
                    gap_extend_score,
                    mismatch_score,
                ));
            }
            GapAlignOpType::Decline => {}
        }
    }
    hsp.num_ident = num_ident;
    score
}

fn mapper_compute_alignment_score(
    hsp: &mut Hsp,
    mismatch_score: i32,
    gap_open_score: i32,
    gap_extend_score: i32,
) -> i32 {
    let Some(info) = hsp.map_info.as_ref() else {
        return mapper_fallback_score_from_script(
            hsp,
            mismatch_score,
            gap_open_score,
            gap_extend_score,
        );
    };
    let edit_slice = info
        .edits
        .as_ref()
        .map(|edits| edits.edits.as_slice())
        .unwrap_or(&[]);

    let gap_base = 15;
    let mut last_pos = hsp.query_offset;
    let mut score: i32 = 0;
    let mut num_identical: i32 = 0;
    let mut query_gap: i32 = 0;
    let mut subject_gap: i32 = 0;

    for edit in edit_slice {
        let num_matches = edit.query_pos.saturating_sub(last_pos);
        last_pos = edit.query_pos;
        score = score.saturating_add(num_matches);
        num_identical = num_identical.saturating_add(num_matches);

        if edit.query_base == gap_base {
            query_gap = query_gap.saturating_add(1);
            if subject_gap > 0 {
                score = score.saturating_add(s_compute_gap_score(subject_gap, -12, -1, -4));
                subject_gap = 0;
            }
        } else if edit.subject_base == gap_base {
            subject_gap = subject_gap.saturating_add(1);
            last_pos = last_pos.saturating_add(1);
            if query_gap > 0 {
                score = score.saturating_add(s_compute_gap_score(query_gap, -12, -1, -4));
                query_gap = 0;
            }
        } else {
            score = score.saturating_add(mismatch_score);
            last_pos = last_pos.saturating_add(1);
            if subject_gap > 0 {
                score = score.saturating_add(s_compute_gap_score(subject_gap, -12, -1, -4));
                subject_gap = 0;
            }
            if query_gap > 0 {
                score = score.saturating_add(s_compute_gap_score(query_gap, -12, -1, -4));
                query_gap = 0;
            }
        }
    }

    if subject_gap > 0 {
        score = score.saturating_add(s_compute_gap_score(subject_gap, -12, -1, -4));
    }
    if query_gap > 0 {
        score = score.saturating_add(s_compute_gap_score(query_gap, -12, -1, -4));
    }

    let tail_matches = hsp.query_end.saturating_sub(last_pos);
    score = score.saturating_add(tail_matches);
    num_identical = num_identical.saturating_add(tail_matches);
    hsp.num_ident = num_identical;
    score
}

/// Port of the edit-script/coordinate core of NCBI `s_TrimHSP`
/// (`hspfilter_mapper.c:741`).
///
/// Rewrites represented mapper `map_info` Jumper edits and subject-overhang
/// buffers alongside the BLAST gap script, matching NCBI's side effects.
pub fn s_trim_hsp(
    hsp: &mut Hsp,
    num: i32,
    is_query: bool,
    is_start: bool,
    mismatch_score: i32,
    gap_open_score: i32,
    gap_extend_score: i32,
    query_seq: Option<&[u8]>,
) -> i32 {
    if num <= 0 {
        return 0;
    }

    let query_span = hsp.query_end.saturating_sub(hsp.query_offset).max(0);
    let subject_span = hsp.subject_end.saturating_sub(hsp.subject_offset).max(0);
    if (is_query && num > query_span) || (!is_query && num > subject_span) {
        return -1;
    }

    let mut delta_query: i32 = 0;
    let mut delta_subject: i32 = 0;
    if let Some(script) = hsp.edit_script.as_mut() {
        let mut num_left = num;
        while num_left > 0 && !script.is_empty() {
            let idx = if is_start { 0 } else { script.len() - 1 };
            let (op, count) = script.get(idx).unwrap_or((GapAlignOpType::Sub, 0));
            if count <= 0 {
                script.remove(idx);
                continue;
            }

            let consumes_target = if is_query {
                mapper_trim_op_consumes_query(op)
            } else {
                mapper_trim_op_consumes_subject(op)
            };
            if consumes_target {
                let trim = count.min(num_left);
                let (dq, ds) = mapper_trim_op_deltas(op, trim);
                delta_query = delta_query.saturating_add(dq);
                delta_subject = delta_subject.saturating_add(ds);
                num_left = num_left.saturating_sub(trim);
                if trim == count {
                    script.remove(idx);
                } else {
                    script.set_num(idx, count - trim);
                }
            } else {
                let (dq, ds) = mapper_trim_op_deltas(op, count);
                delta_query = delta_query.saturating_add(dq);
                delta_subject = delta_subject.saturating_add(ds);
                script.remove(idx);
            }
        }
        if num_left > 0 {
            return -1;
        }
    } else {
        if is_query {
            delta_query = num;
        } else {
            delta_subject = num;
        }
    }

    mapper_trim_map_info(hsp, delta_query, delta_subject, is_start, query_seq);

    if is_start {
        hsp.query_offset = hsp.query_offset.saturating_add(delta_query);
        hsp.subject_offset = hsp.subject_offset.saturating_add(delta_subject);
        hsp.query_gapped_start = hsp.query_gapped_start.max(hsp.query_offset);
        hsp.subject_gapped_start = hsp.subject_gapped_start.max(hsp.subject_offset);
    } else {
        hsp.query_end = hsp.query_end.saturating_sub(delta_query);
        hsp.subject_end = hsp.subject_end.saturating_sub(delta_subject);
    }

    hsp.score =
        mapper_compute_alignment_score(hsp, mismatch_score, gap_open_score, gap_extend_score);
    0
}

fn mapper_trim_map_info(
    hsp: &mut Hsp,
    delta_query: i32,
    delta_subject: i32,
    is_start: bool,
    query_seq: Option<&[u8]>,
) {
    let Some(info) = hsp.map_info.as_mut() else {
        return;
    };
    let gap_base = 15;

    if let Some(query_seq) = query_seq {
        let overhangs = info.subject_overhangs.get_or_insert_with(Default::default);
        if is_start {
            if delta_subject > 0 {
                let old_len = overhangs.left.as_ref().map_or(0, Vec::len);
                let left = overhangs.left.get_or_insert_with(Vec::new);
                left.resize(old_len + delta_subject as usize, 0);
                if let Some(bases) = query_seq.get(
                    hsp.query_offset as usize
                        ..hsp
                            .query_offset
                            .saturating_add(delta_subject)
                            .max(hsp.query_offset) as usize,
                ) {
                    left[old_len..old_len + bases.len()].copy_from_slice(bases);
                }
                overhangs.left_len = left.len() as i32;

                let mut offset: i32 = 0;
                if let Some(edits) = info.edits.as_ref() {
                    for edit in edits
                        .edits
                        .iter()
                        .take_while(|edit| edit.query_pos < hsp.query_offset)
                    {
                        if edit.subject_base == gap_base {
                            offset = offset.saturating_sub(1);
                        } else {
                            let pos = edit.query_pos.saturating_add(offset);
                            if pos >= 0 {
                                if let Some(slot) = left.get_mut(old_len + pos as usize) {
                                    *slot = edit.subject_base;
                                }
                            }
                        }
                    }
                }
            }
        } else if delta_subject > 0 {
            let old_right = overhangs.right.take().unwrap_or_default();
            let mut right = vec![0; delta_subject as usize + old_right.len()];
            right[delta_subject as usize..].copy_from_slice(&old_right);
            if let Some(bases) = query_seq.get(
                hsp.query_end.saturating_sub(delta_subject).max(0) as usize
                    ..hsp.query_end.max(0) as usize,
            ) {
                right[..bases.len()].copy_from_slice(bases);
            }

            let mut k = 0;
            if let Some(edits) = info.edits.as_ref() {
                while k < edits.edits.len()
                    && edits.edits[k].query_pos < hsp.query_end.saturating_sub(delta_query)
                {
                    k += 1;
                }

                let mut offset = hsp.query_end.saturating_sub(delta_query).saturating_neg();
                for edit in edits.edits.iter().skip(k) {
                    if edit.query_pos >= hsp.query_end {
                        continue;
                    }
                    if edit.subject_base == gap_base {
                        offset = offset.saturating_sub(1);
                    } else {
                        let pos = edit.query_pos.saturating_add(offset);
                        if pos >= 0 {
                            if let Some(slot) = right.get_mut(pos as usize) {
                                *slot = edit.subject_base;
                            }
                        }
                    }
                }
            }
            overhangs.right = Some(right);
            overhangs.right_len = overhangs.right.as_ref().map_or(0, Vec::len) as i32;
        }
    }

    let Some(edits) = info.edits.as_mut() else {
        return;
    };
    let Some(script) = hsp.edit_script.as_ref() else {
        return;
    };

    if is_start {
        let mut k = edits.edits.len();
        let mut p = hsp.query_end.saturating_sub(1);
        for (op, count) in script.iter().rev() {
            if !matches!(
                op,
                GapAlignOpType::Del | GapAlignOpType::Del1 | GapAlignOpType::Del2
            ) {
                p = p.saturating_sub(count);
                while k > 0
                    && edits.edits[k - 1].query_pos > p
                    && edits.edits[k - 1].query_base != gap_base
                {
                    k -= 1;
                }
            } else {
                for _ in 0..count {
                    if k > 0 {
                        k -= 1;
                    }
                }
            }
        }
        if k > 0 {
            edits.edits.drain(0..k);
        }
    } else {
        let mut k = 0usize;
        let mut p = hsp.query_offset;
        for (op, count) in script.iter() {
            if !matches!(
                op,
                GapAlignOpType::Del | GapAlignOpType::Del1 | GapAlignOpType::Del2
            ) {
                p = p.saturating_add(count);
                while k < edits.edits.len()
                    && edits.edits[k].query_pos < p
                    && edits.edits[k].query_base != gap_base
                {
                    k += 1;
                }
            } else {
                for _ in 0..count {
                    if k < edits.edits.len() {
                        k += 1;
                    }
                }
            }
        }
        edits.edits.truncate(k);
        edits.num_edits = edits.edits.len() as i32;
    }
}

/// Port of NCBI `s_TrimOverlap` (`hspfilter_mapper.c:3407`).
///
/// Audit caveat: overlap removal delegates to `s_trim_hsp`, so it inherits its
/// simplified mapper edit/overhang trimming behavior.
pub fn s_trim_overlap(first: &mut Hsp, second: &mut Hsp, query: Option<&[u8]>) -> i32 {
    let query_overlap = first.query_end.saturating_sub(second.query_offset);
    if query_overlap > 0 {
        let second_span = second.query_end.saturating_sub(second.query_offset);
        let status = if second_span > query_overlap {
            s_trim_hsp(second, query_overlap, true, true, -4, -4, -4, query)
        } else {
            s_trim_hsp(first, query_overlap, true, false, -4, -4, -4, query)
        };
        if status != 0 {
            return status;
        }
    }

    let subject_overlap = first.subject_end.saturating_sub(second.subject_offset);
    if subject_overlap > 0 {
        let second_span = second.subject_end.saturating_sub(second.subject_offset);
        let status = if second_span > subject_overlap {
            s_trim_hsp(second, subject_overlap, false, true, -4, -4, -4, query)
        } else {
            s_trim_hsp(first, subject_overlap, false, false, -4, -4, -4, query)
        };
        if status != 0 {
            return status;
        }
    }

    if first.query_end > second.query_offset || first.subject_end > second.subject_offset {
        return -1;
    }
    0
}

/// Port of NCBI `s_TrimChainStartToSubjPos` (`hspfilter_mapper.c:4055`).
///
/// Rust trims/removes owned HSP containers, refreshes scores, and clears mapper
/// splice/exon bits on `map_info->left_edge`.
pub fn s_trim_chain_start_to_subj_pos(
    chain: Option<&mut HSPChain>,
    subj_pos: i32,
    mismatch_score: i32,
    gap_open_score: i32,
    gap_extend_score: i32,
    query: Option<&[u8]>,
) -> i32 {
    let Some(chain) = chain else {
        return 0;
    };
    if subj_pos < 0 {
        return 0;
    }

    while chain
        .hsps
        .as_ref()
        .is_some_and(|head| head.hsp.subject_end < subj_pos)
    {
        let Some(mut head) = chain.hsps.take() else {
            break;
        };
        chain.score -= head.hsp.score;
        chain.hsps = head.next.take();
    }

    let Some(head) = chain.hsps.as_mut() else {
        return -1;
    };
    if head.hsp.subject_offset >= subj_pos {
        return 0;
    }

    let num_left = subj_pos - head.hsp.subject_offset;
    let old_score = head.hsp.score;
    let status = s_trim_hsp(
        &mut head.hsp,
        num_left,
        false,
        true,
        mismatch_score,
        gap_open_score,
        gap_extend_score,
        query,
    );
    if status != 0 {
        return status;
    }
    chain.score -= old_score - head.hsp.score;
    mapper_clear_left_splice_exon(&mut head.hsp);

    let remove_head = chain
        .hsps
        .as_ref()
        .and_then(|head| head.next.as_ref().map(|next| (head, next)))
        .is_some_and(|(head, next)| head.hsp.query_offset >= next.hsp.query_offset);
    if remove_head {
        if let Some(mut head) = chain.hsps.take() {
            chain.hsps = head.next.take();
        }
    }

    0
}

/// Port of NCBI `s_TrimChainEndToSubjPos` (`hspfilter_mapper.c:4131`).
///
/// Rust trims/removes owned HSP containers, refreshes scores, and clears mapper
/// splice/exon bits on `map_info->right_edge`.
pub fn s_trim_chain_end_to_subj_pos(
    chain: Option<&mut HSPChain>,
    subj_pos: i32,
    mismatch_score: i32,
    gap_open_score: i32,
    gap_extend_score: i32,
    query: Option<&[u8]>,
) -> i32 {
    let Some(chain) = chain else {
        return -1;
    };
    if subj_pos <= 0 || query.is_none() {
        return -1;
    }
    if chain.hsps.is_none() {
        return -1;
    }

    let mut h_ptr = chain.hsps.as_deref_mut().unwrap() as *mut HSPContainer;
    unsafe {
        while let Some(next) = (*h_ptr).next.as_deref_mut() {
            if next.hsp.subject_end < subj_pos {
                h_ptr = next as *mut HSPContainer;
            } else {
                break;
            }
        }

        if let Some(next) = (*h_ptr).next.as_deref_mut() {
            if next.hsp.subject_offset < subj_pos {
                h_ptr = next as *mut HSPContainer;
            }
        }

        let h = &mut *h_ptr;
        let mut tail = h.next.take();
        while let Some(mut current) = tail {
            chain.score -= current.hsp.score;
            tail = current.next.take();
        }

        if h.hsp.subject_end <= subj_pos {
            return 0;
        }

        let num_left = h.hsp.subject_end.saturating_sub(subj_pos);
        let old_score = h.hsp.score;
        let status = s_trim_hsp(
            &mut h.hsp,
            num_left,
            false,
            false,
            mismatch_score,
            gap_open_score,
            gap_extend_score,
            query,
        );
        if status != 0 {
            return status;
        }
        chain.score -= old_score - h.hsp.score;
        mapper_clear_right_splice_exon(&mut h.hsp);

        if let Some(prev) = find_container_before(&mut chain.hsps, h_ptr) {
            if prev.hsp.query_end >= h.hsp.query_end {
                chain.score -= h.hsp.score;
                prev.next = None;
            }
        }
    }

    0
}

fn find_container_before(
    list: &mut Option<Box<HSPContainer>>,
    target: *const HSPContainer,
) -> Option<&mut HSPContainer> {
    let mut cursor = list.as_deref_mut();
    while let Some(current) = cursor {
        let is_previous = current
            .next
            .as_deref()
            .is_some_and(|next| std::ptr::eq(next as *const HSPContainer, target));
        if is_previous {
            return Some(current);
        }
        cursor = current.next.as_deref_mut();
    }
    None
}

/// Port of the chain/HSP trimming core of NCBI `s_SetAdapter`
/// (`hspfilter_mapper.c:1096`).
///
/// Marks mapper adapter/exon edge flags on `Hsp::map_info`, matching C's
/// `hsp->map_info` side effect.
pub fn s_set_adapter(
    chains: &mut Option<Box<HSPChain>>,
    adapter_pos: i32,
    query: Option<&[u8]>,
    query_len: i32,
    scores: Option<&ScoringOptions>,
) -> i16 {
    const MIN_ADAPTER_LEN: i32 = 3;
    if chains.is_none() || adapter_pos < 0 {
        return -1;
    }
    let Some(scores) = scores else {
        return -1;
    };

    let mut kept = Vec::new();
    for mut chain in drain_chain_list(chains.take()) {
        let Some(first) = chain.hsps.as_deref() else {
            continue;
        };

        if first.hsp.query_frame > 0 {
            let mut last = first;
            while let Some(next) = last.next.as_deref() {
                last = next;
            }
            if query_len.saturating_sub(last.hsp.query_end) < MIN_ADAPTER_LEN {
                kept.push(chain);
                continue;
            }

            chain.adapter = adapter_pos;
            if first.hsp.query_offset >= adapter_pos - 5 {
                continue;
            }

            let mut hsps = drain_hsp_container_list(chain.hsps.take());
            let mut recompute_score = false;
            if let Some(idx) = hsps
                .iter()
                .position(|container| container.hsp.query_end > adapter_pos)
            {
                if hsps[idx].hsp.query_offset < adapter_pos {
                    let trim_by = hsps[idx].hsp.query_end.saturating_sub(adapter_pos);
                    let old_score = hsps[idx].hsp.score;
                    let _ = s_trim_hsp(
                        &mut hsps[idx].hsp,
                        trim_by,
                        true,
                        false,
                        scores.penalty,
                        scores.gap_open,
                        scores.gap_extend,
                        query,
                    );
                    chain.score += hsps[idx].hsp.score - old_score;
                    hsps.truncate(idx + 1);
                } else {
                    hsps.truncate(idx);
                }
                recompute_score = true;
            }

            if let Some(last) = hsps.last_mut() {
                mapper_mark_right_edge(&mut last.hsp, MAPPER_ADAPTER | MAPPER_EXON);
            }
            chain.hsps = build_hsp_container_list(hsps);
            if recompute_score {
                let _ = s_compute_chain_score(Some(&mut chain), Some(scores), query_len, false);
            }

            if chain.hsps.is_some() {
                kept.push(chain);
            }
        } else {
            if first.hsp.query_offset < MIN_ADAPTER_LEN {
                kept.push(chain);
                continue;
            }

            chain.adapter = adapter_pos;
            let pos_minus = query_len.saturating_sub(adapter_pos).saturating_sub(1);
            let mut hsps = drain_hsp_container_list(chain.hsps.take());
            let Some(idx) = hsps
                .iter()
                .position(|container| container.hsp.query_end > pos_minus + 5)
            else {
                continue;
            };

            if idx > 0 {
                hsps.drain(0..idx);
            }
            if let Some(first_kept) = hsps.first_mut() {
                mapper_mark_left_edge(&mut first_kept.hsp, MAPPER_ADAPTER | MAPPER_EXON);
                if pos_minus >= first_kept.hsp.query_offset {
                    let trim_by = pos_minus
                        .saturating_sub(first_kept.hsp.query_offset)
                        .saturating_add(1);
                    let old_score = first_kept.hsp.score;
                    let _ = s_trim_hsp(
                        &mut first_kept.hsp,
                        trim_by,
                        true,
                        true,
                        scores.penalty,
                        scores.gap_open,
                        scores.gap_extend,
                        query,
                    );
                    chain.score += first_kept.hsp.score - old_score;
                }
            }
            chain.hsps = build_hsp_container_list(hsps);
            let _ = s_compute_chain_score(Some(&mut chain), Some(scores), query_len, false);
            if chain.hsps.is_some() {
                kept.push(chain);
            }
        }
    }

    *chains = build_chain_list(kept);
    0
}

/// Port of NCBI `s_FindAdapterInSequence` (`hspfilter_mapper.c:876`).
pub fn s_find_adapter_in_sequence(
    hsp_from: i32,
    hsp_to: i32,
    query: Option<&[u8]>,
    query_len: i32,
) -> i32 {
    const ADAPTERS: [&[u8]; 4] = [
        &[0, 2, 0, 3, 1, 2, 2, 0, 0, 2, 0, 2],
        &[0, 3, 2, 2, 0, 0, 3, 3, 1, 3, 1, 2],
        &[1, 3, 2, 3, 1, 3, 1, 3, 3, 0, 3, 0],
        &[2, 0, 3, 1, 2, 2, 0, 0, 2, 0, 2, 1, 0, 1, 0, 1, 2, 3, 1, 3],
    ];

    let Some(query) = query else {
        return -1;
    };
    let to = hsp_to.min(query_len).min(query.len() as i32);
    let query_limit = query_len.min(query.len() as i32);
    let from = hsp_from.max(0);
    if from >= to {
        return -1;
    }

    for adapter in ADAPTERS {
        let len = adapter.len() as i32;
        let mut q = (to - len).max(from);
        while q < query_limit - 4 {
            if query[q as usize..(q + 4) as usize] != adapter[..4] {
                q += 1;
                continue;
            }
            let mut errors = 2;
            let mut pos = q;
            while pos < query_limit && pos - q < len {
                let query_base = query[pos as usize];
                let adapter_base = adapter[(pos - q) as usize];
                if query_base != adapter_base {
                    errors -= 1;
                    if errors == 0 {
                        break;
                    }
                }
                pos += 1;
            }
            if pos == query_limit || pos - q == len {
                return q;
            }
            q += 1;
        }
    }

    -1
}

fn chain_query_coverage(chain: &HSPChain, query_len: i32) -> Option<(i32, i32)> {
    let first = chain.hsps.as_deref()?;
    if first.hsp.query_frame > 0 {
        let mut last = first;
        while let Some(next) = last.next.as_deref() {
            last = next;
        }
        Some((first.hsp.query_offset, last.hsp.query_end))
    } else {
        let mut last = first;
        while let Some(next) = last.next.as_deref() {
            last = next;
        }
        Some((
            query_len.saturating_sub(last.hsp.query_end),
            query_len.saturating_sub(first.hsp.query_offset),
        ))
    }
}

fn chain_adapter_overhang(chain: &HSPChain, query_len: i32) -> Option<i32> {
    let first = chain.hsps.as_deref()?;
    if first.hsp.query_frame > 0 {
        let mut last = first;
        while let Some(next) = last.next.as_deref() {
            last = next;
        }
        Some(query_len.saturating_sub(last.hsp.query_end))
    } else {
        Some(first.hsp.query_offset + 1)
    }
}

fn chain_adapter_search_bounds(chain: &HSPChain, query_len: i32) -> Option<(i32, i32)> {
    let mut hsp = chain.hsps.as_deref()?;
    if hsp.hsp.query_frame > 0 {
        while let Some(next) = hsp.next.as_deref() {
            hsp = next;
        }
        Some((hsp.hsp.query_offset, hsp.hsp.query_end))
    } else {
        Some((
            query_len.saturating_sub(hsp.hsp.query_end),
            query_len.saturating_sub(hsp.hsp.query_offset),
        ))
    }
}

/// Port of NCBI `s_FindAdapters` (`hspfilter_mapper.c:1291`) over the Rust
/// chain-list representation.
///
/// Adapter trimming delegates to [`s_set_adapter`], including C-shaped mapper
/// adapter/exon edge-bit annotations in `map_info`.
pub fn s_find_adapters(
    saved: &mut [Option<Box<HSPChain>>],
    query_sequence: &[u8],
    query_info: &QueryInfo,
    score_opts: &ScoringOptions,
) -> i32 {
    for query_idx in 0..query_info.num_queries.max(0) as usize {
        if query_idx >= saved.len() || saved[query_idx].is_none() {
            continue;
        }
        let plus_ctx = query_idx * NUM_STRANDS;
        let Some(context) = query_info.contexts.get(plus_ctx) else {
            continue;
        };
        let query_len = context.query_length.max(0);
        let query_start = context.query_offset.max(0) as usize;
        let query_end = query_start.saturating_add(query_len as usize);
        let Some(query) = query_sequence.get(query_start..query_end) else {
            continue;
        };

        let mut best = saved[query_idx].as_deref();
        let mut cursor = best.and_then(|chain| chain.next.as_deref());
        while let Some(chain) = cursor {
            if let Some(current_best) = best {
                if chain.score > current_best.score {
                    best = Some(chain);
                }
            }
            cursor = chain.next.as_deref();
        }

        let Some((from, to)) = best.and_then(|chain| chain_query_coverage(chain, query_len)) else {
            continue;
        };
        if from < 20 && to > query_len.saturating_sub(3) {
            continue;
        }

        let mut longest_overhang = 0;
        let mut search_bounds =
            best.and_then(|chain| chain_adapter_search_bounds(chain, query_len));
        let mut cursor = saved[query_idx]
            .as_deref()
            .and_then(|chain| chain.next.as_deref());
        while let Some(chain) = cursor {
            if let Some(overhang) = chain_adapter_overhang(chain, query_len) {
                if overhang > longest_overhang {
                    longest_overhang = overhang;
                    search_bounds = chain_adapter_search_bounds(chain, query_len);
                }
            }
            cursor = chain.next.as_deref();
        }

        let Some((from, to)) = search_bounds else {
            continue;
        };
        if to >= query_len.saturating_sub(3) {
            continue;
        }

        let adapter_pos = s_find_adapter_in_sequence(from, to, Some(query), query_len);
        if adapter_pos >= 0 {
            let _ = s_set_adapter(
                &mut saved[query_idx],
                adapter_pos,
                Some(query),
                query_len,
                Some(score_opts),
            );
        }
    }

    0
}

/// Port of the representable core of NCBI `s_MergeHSPs`
/// (`hspfilter_mapper.c:1580`).
///
/// Also merges represented mapper `map_info` edits, right-edge marker, and
/// right subject-overhang buffer, matching NCBI's side effects.
pub fn s_merge_hsps(
    first: &Hsp,
    second: &Hsp,
    query: Option<&[u8]>,
    score_opts: Option<&ScoringOptions>,
) -> Option<Hsp> {
    let score_opts = score_opts?;
    if second.query_offset < first.query_end || second.subject_offset < first.subject_end {
        return None;
    }

    let mut query_gap = second.subject_offset.saturating_sub(first.subject_end);
    let mut subject_gap = second.query_offset.saturating_sub(first.query_end);
    let mut mismatches = 0;
    if query_gap.max(subject_gap) < 4 {
        mismatches = query_gap.min(subject_gap);
        query_gap = query_gap.saturating_sub(mismatches);
        subject_gap = subject_gap.saturating_sub(mismatches);
    }

    let mut merged = first.clone();
    let mut script = first
        .edit_script
        .clone()
        .unwrap_or_else(crate::gapinfo::GapEditScript::new);
    if script.is_empty() && first.query_end > first.query_offset {
        script.push(
            GapAlignOpType::Sub,
            first
                .query_end
                .saturating_sub(first.query_offset)
                .min(first.subject_end.saturating_sub(first.subject_offset)),
        );
    }
    if mismatches > 0 {
        script.push(GapAlignOpType::Sub, mismatches);
    }
    if query_gap > 0 {
        script.push(GapAlignOpType::Del, query_gap);
    }
    if subject_gap > 0 {
        script.push(GapAlignOpType::Ins, subject_gap);
    }
    if let Some(second_script) = second.edit_script.as_ref() {
        for (op, count) in second_script.iter() {
            script.push(op, count);
        }
    } else if second.query_end > second.query_offset {
        script.push(
            GapAlignOpType::Sub,
            second
                .query_end
                .saturating_sub(second.query_offset)
                .min(second.subject_end.saturating_sub(second.subject_offset)),
        );
    }

    merge_hsp_map_info(
        &mut merged,
        second,
        query,
        mismatches,
        query_gap,
        subject_gap,
    );
    merged.query_end = second.query_end;
    merged.subject_end = second.subject_end;
    merged.edit_script = Some(script);
    merged.score = mapper_compute_alignment_score(
        &mut merged,
        score_opts.penalty,
        score_opts.gap_open,
        score_opts.gap_extend,
    );
    merged.num_gaps = merged
        .edit_script
        .as_ref()
        .map(|script| {
            script
                .iter()
                .filter(|(op, _)| !matches!(op, GapAlignOpType::Sub | GapAlignOpType::Decline))
                .map(|(_, count)| count)
                .sum()
        })
        .unwrap_or(0);
    Some(merged)
}

fn merge_hsp_map_info(
    merged: &mut Hsp,
    second: &Hsp,
    query: Option<&[u8]>,
    mismatches: i32,
    query_gap: i32,
    subject_gap: i32,
) {
    let Some(mut info) = merged.map_info.clone() else {
        return;
    };

    let mut edits = info.edits.take().unwrap_or_default();
    let first_edit_count = edits.edits.len();
    let gap_base = 15;
    let mut offset = merged.subject_offset.saturating_sub(merged.query_offset);
    for edit in edits.edits.iter().take(first_edit_count) {
        if edit.query_base == gap_base {
            offset = offset.saturating_add(1);
        } else if edit.subject_base == gap_base {
            offset = offset.saturating_sub(1);
        }
    }

    let right_overhang = info
        .subject_overhangs
        .as_ref()
        .and_then(|overhangs| overhangs.right.as_ref());

    for k in 0..mismatches {
        let query_pos = merged.query_end.saturating_add(k);
        let query_base = query
            .and_then(|query| query.get(query_pos as usize).copied())
            .unwrap_or(0);
        let subject_pos = query_pos
            .saturating_add(offset)
            .saturating_sub(merged.subject_end);
        let subject_base = right_overhang
            .and_then(|bases| bases.get(subject_pos as usize).copied())
            .unwrap_or(query_base);
        edits.edits.push(JumperEdit {
            query_pos,
            query_base,
            subject_base,
        });
    }

    for _ in 0..query_gap {
        let query_pos = merged.query_end.saturating_add(mismatches);
        let subject_pos = query_pos
            .saturating_add(offset)
            .saturating_sub(merged.subject_end);
        let subject_base = right_overhang
            .and_then(|bases| bases.get(subject_pos as usize).copied())
            .unwrap_or(0);
        edits.edits.push(JumperEdit {
            query_pos,
            query_base: gap_base,
            subject_base,
        });
        offset = offset.saturating_add(1);
    }

    for k in 0..subject_gap {
        let query_pos = merged
            .query_end
            .saturating_add(mismatches)
            .saturating_add(k);
        let query_base = query
            .and_then(|query| query.get(query_pos as usize).copied())
            .unwrap_or(0);
        edits.edits.push(JumperEdit {
            query_pos,
            query_base,
            subject_base: gap_base,
        });
    }

    if let Some(second_info) = second.map_info.as_ref() {
        if let Some(second_edits) = second_info.edits.as_ref() {
            edits.edits.extend_from_slice(&second_edits.edits);
        }
        edits.num_edits = edits.edits.len() as i32;
        info.right_edge = second_info.right_edge;
        if let Some(overhangs) = info.subject_overhangs.as_mut() {
            overhangs.right = second_info
                .subject_overhangs
                .as_ref()
                .and_then(|second_overhangs| second_overhangs.right.clone());
            overhangs.right_len = overhangs.right.as_ref().map_or(0, Vec::len) as i32;
        }
    }

    info.edits = Some(edits);
    merged.map_info = Some(info);
}

/// Port of NCBI `s_IntronToGap` (`hspfilter_mapper.c:3456`) for the Rust
/// owned linked-list representation.
///
/// Rust trims and merges the BLAST `Hsp` coordinates/gap script, relinks owned
/// `HSPContainer` nodes, and clears mapper splice edge bits on merge failure.
pub fn s_intron_to_gap(
    h: &mut HSPContainer,
    query: Option<&[u8]>,
    scoring_opts: Option<&ScoringOptions>,
) -> i32 {
    let Some(scoring_opts) = scoring_opts else {
        return -1;
    };
    if query.is_none() {
        return -1;
    }

    let Some(mut next) = h.next.take() else {
        return -1;
    };
    let following = next.next.take();

    let status = s_trim_overlap(&mut h.hsp, &mut next.hsp, query);
    if status != 0 {
        next.next = following;
        h.next = Some(next);
        return status;
    }

    if let Some(new_hsp) = s_merge_hsps(&h.hsp, &next.hsp, query, Some(scoring_opts)) {
        h.hsp = new_hsp;
        h.next = following;
    } else {
        mapper_clear_splice_signal(&mut h.hsp, &mut next.hsp);
        next.next = following;
        h.next = Some(next);
    }

    0
}

/// Conservative Rust port of NCBI `s_FindSpliceJunctions`
/// (`hspfilter_mapper.c:3518`) for owned `HSPChain` lists.
///
/// This keeps the same chain/pair traversal, overlap trimming,
/// short-intron-to-gap merge, chain score recomputation, and represented
/// `map_info`-backed overlap/gap splice discovery.
pub fn s_find_splice_junctions(
    chains: Option<&mut HSPChain>,
    query: Option<&[u8]>,
    query_len: i32,
    scoring_opts: Option<&ScoringOptions>,
) -> i32 {
    let Some(mut chain) = chains else {
        return -1;
    };
    let Some(query) = query else {
        return -1;
    };
    let Some(scoring_opts) = scoring_opts else {
        return -1;
    };

    loop {
        let Some(first_container) = chain.hsps.as_deref_mut() else {
            let _ = s_compute_chain_score(Some(chain), Some(scoring_opts), query_len, true);
            if let Some(next_chain) = chain.next.as_deref_mut() {
                chain = next_chain;
                continue;
            }
            break;
        };

        let mut h_ptr: *mut HSPContainer = first_container;
        // The C routine mutates a singly linked list while sometimes keeping
        // the current node after merging it with its successor. A raw pointer
        // keeps that control flow local without fighting Rust's borrow checker;
        // the pointer is always derived from the live owned chain list.
        unsafe {
            while (*h_ptr).next.is_some() {
                let h = &mut *h_ptr;
                let (query_gap, overlaps_query, can_short_merge, score_pair) = {
                    let next = h.next.as_ref().expect("checked above");
                    let query_gap = next.hsp.query_offset as i64 - h.hsp.query_end as i64;
                    let subject_gap = next.hsp.subject_offset as i64 - h.hsp.subject_end as i64;
                    let right_overhang_len = h
                        .hsp
                        .map_info
                        .as_ref()
                        .and_then(|info| info.subject_overhangs.as_ref())
                        .and_then(|overhangs| overhangs.right.as_ref())
                        .map_or(0, |right| right.len() as i64);
                    (
                        query_gap,
                        next.hsp.query_offset <= h.hsp.query_end
                            && next.hsp.query_offset > h.hsp.query_offset,
                        subject_gap - query_gap < 30 && subject_gap < right_overhang_len,
                        h.hsp.score > 50 && next.hsp.score > 50,
                    )
                };

                if can_short_merge {
                    if query_gap > 1 {
                        let status =
                            mapper_extend_hsp(&mut h.hsp, query_gap as i32, 0, false, scoring_opts);
                        if status != 0 {
                            return status;
                        }
                    }
                    let old_next_start = h.next.as_ref().map(|next| next.hsp.query_offset);
                    let status = s_intron_to_gap(h, Some(query), Some(scoring_opts));
                    if status != 0 {
                        return status;
                    }
                    if h.next.as_ref().map(|next| next.hsp.query_offset) == old_next_start {
                        h_ptr = h.next.as_deref_mut().expect("checked above") as *mut HSPContainer;
                    }
                    continue;
                }

                if overlaps_query {
                    let status = {
                        let next = h.next.as_deref_mut().expect("checked above");
                        s_find_splice_junctions_for_overlaps(
                            &mut h.hsp,
                            &mut next.hsp,
                            Some(query),
                            query_len,
                            !score_pair,
                        )
                    };
                    if status != 0 {
                        return status;
                    }
                    let has_splice = h
                        .hsp
                        .map_info
                        .as_ref()
                        .is_some_and(|info| info.right_edge & MAPPER_SPLICE_SIGNAL != 0);
                    if has_splice {
                        let should_merge = h
                            .next
                            .as_ref()
                            .is_some_and(|next| next.hsp.subject_offset - h.hsp.subject_end < 30);
                        if should_merge {
                            let old_next_start = h.next.as_ref().map(|next| next.hsp.query_offset);
                            let status = s_intron_to_gap(h, Some(query), Some(scoring_opts));
                            if status != 0 {
                                return status;
                            }
                            if h.next.as_ref().map(|next| next.hsp.query_offset) != old_next_start {
                                continue;
                            }
                        }
                    } else {
                        let status = {
                            let next = h.next.as_deref_mut().expect("checked above");
                            s_trim_overlap(&mut h.hsp, &mut next.hsp, Some(query))
                        };
                        if status != 0 {
                            return status;
                        }
                    }
                    h_ptr = h.next.as_deref_mut().expect("checked above") as *mut HSPContainer;
                    continue;
                }

                if query_gap > 0 {
                    let can_try_gap_splice = h
                        .hsp
                        .map_info
                        .as_ref()
                        .and_then(|info| info.subject_overhangs.as_ref())
                        .and_then(|overhangs| overhangs.right.as_ref())
                        .is_some_and(|right| query_gap < right.len() as i64)
                        && h.next.as_ref().is_some_and(|next| {
                            next.hsp
                                .map_info
                                .as_ref()
                                .and_then(|info| info.subject_overhangs.as_ref())
                                .and_then(|overhangs| overhangs.left.as_ref())
                                .is_some_and(|left| query_gap < left.len() as i64)
                        });
                    if can_try_gap_splice {
                        let status = {
                            let next = h.next.as_deref_mut().expect("checked above");
                            s_find_splice_junctions_for_gap_using_map_info(
                                &mut h.hsp,
                                &mut next.hsp,
                                Some(query),
                                query_len,
                                Some(scoring_opts),
                            )
                        };
                        if status != 0 {
                            return status;
                        }
                    }
                    let has_splice = h
                        .hsp
                        .map_info
                        .as_ref()
                        .is_some_and(|info| info.right_edge & MAPPER_SPLICE_SIGNAL != 0);
                    if has_splice {
                        let should_merge = h
                            .next
                            .as_ref()
                            .is_some_and(|next| next.hsp.subject_offset - h.hsp.subject_end < 30);
                        if should_merge {
                            let old_next_start = h.next.as_ref().map(|next| next.hsp.query_offset);
                            let status = s_intron_to_gap(h, Some(query), Some(scoring_opts));
                            if status != 0 {
                                return status;
                            }
                            if h.next.as_ref().map(|next| next.hsp.query_offset) != old_next_start {
                                continue;
                            }
                        }
                        h_ptr = h.next.as_deref_mut().expect("checked above") as *mut HSPContainer;
                        continue;
                    }
                    let status = {
                        let next = h.next.as_deref_mut().expect("checked above");
                        s_trim_overlap(&mut h.hsp, &mut next.hsp, Some(query))
                    };
                    if status != 0 {
                        return status;
                    }
                    h_ptr = h.next.as_deref_mut().expect("checked above") as *mut HSPContainer;
                    continue;
                }

                let status = {
                    let next = h.next.as_deref_mut().expect("checked above");
                    s_trim_overlap(&mut h.hsp, &mut next.hsp, Some(query))
                };
                if status != 0 {
                    return status;
                }
                h_ptr = h.next.as_deref_mut().expect("checked above") as *mut HSPContainer;
            }
        }

        let _ = s_compute_chain_score(Some(chain), Some(scoring_opts), query_len, true);
        if let Some(next_chain) = chain.next.as_deref_mut() {
            chain = next_chain;
        } else {
            break;
        }
    }

    0
}

fn hsp_frame_sign(chain: &HSPChain) -> i32 {
    chain
        .hsps
        .as_deref()
        .map(|container| container.hsp.query_frame.signum())
        .unwrap_or(0)
}

fn first_subject_offset(chain: &HSPChain) -> Option<i32> {
    chain
        .hsps
        .as_deref()
        .map(|container| container.hsp.subject_offset)
}

fn last_subject_end(chain: &HSPChain) -> Option<i32> {
    let mut current = chain.hsps.as_deref()?;
    while let Some(next) = current.next.as_deref() {
        current = next;
    }
    Some(current.hsp.subject_end)
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PairInfo {
    pub first: usize,
    pub second: usize,
    pub score: i32,
    pub trim_first: i32,
    pub trim_second: i32,
    pub valid_pair: bool,
    pub distance: i32,
    pub conf: u8,
}

/// Port of NCBI `s_ComparePairs` (`hspfilter_mapper.c`), matching the pair
/// ranking used by `s_FindBestPairs`: higher score first, then lower pair
/// configuration rank, then shorter distance.
pub fn s_compare_pairs(a: &PairInfo, b: &PairInfo) -> i32 {
    if a.score > b.score {
        -1
    } else if a.score < b.score {
        1
    } else if a.conf < b.conf {
        -1
    } else if a.conf > b.conf {
        1
    } else if a.distance < b.distance {
        -1
    } else if a.distance > b.distance {
        1
    } else {
        0
    }
}

/// Conservative Rust port of NCBI `s_FindBestPairs`
/// (`hspfilter_mapper.c:4214`).
///
/// Audit caveat: C stores raw `HSPChain*` pair links in `Pairinfo`, mutates the
/// two input lists in place, and may clone already-paired chains so several
/// chains can point at equivalent mate copies. Rust stores paired list indexes
/// plus cached partner scores, so this ports pair scoring, configuration
/// ranking, selected-pair marking, and overextension trimming for existing
/// owned chains. It does not clone additional mate chains for ambiguous
/// already-paired cases.
pub fn s_find_best_pairs(
    first_list: &mut Option<Box<HSPChain>>,
    second_list: &mut Option<Box<HSPChain>>,
    min_score: i32,
    pair_info: &mut Vec<PairInfo>,
    is_spliced: bool,
    scoring_options: Option<&ScoringOptions>,
    query: Option<&[u8]>,
) -> bool {
    let mut first_chains = drain_chain_list(first_list.take());
    let mut second_chains = drain_chain_list(second_list.take());
    pair_info.clear();

    let max_insert_size = if is_spliced {
        MAGICBLAST_MAX_INSERT_SIZE_SPLICED
    } else {
        MAGICBLAST_MAX_INSERT_SIZE_NONSPLICED
    };

    for (first_idx, first) in first_chains.iter().enumerate() {
        for (second_idx, second) in second_chains.iter().enumerate() {
            let first_frame = hsp_frame_sign(first);
            let second_frame = hsp_frame_sign(second);
            if first_frame == 0 || second_frame == 0 {
                continue;
            }

            let mut info = PairInfo {
                first: first_idx,
                second: second_idx,
                score: first.score + second.score,
                trim_first: 0,
                trim_second: 0,
                valid_pair: false,
                distance: 0,
                conf: PAIR_NONE,
            };

            if first_frame != second_frame {
                let (plus_is_first, plus, minus) = if first_frame > 0 {
                    (true, first.as_ref(), second.as_ref())
                } else {
                    (false, second.as_ref(), first.as_ref())
                };
                let Some(plus_start) = first_subject_offset(plus) else {
                    continue;
                };
                let Some(minus_start) = last_subject_end(minus) else {
                    continue;
                };
                let distance = minus_start - plus_start;
                info.distance = distance;

                if distance > 0 && distance < max_insert_size {
                    let Some(plus_end) = last_subject_end(plus) else {
                        continue;
                    };
                    let Some(minus_end) = first_subject_offset(minus) else {
                        continue;
                    };

                    if plus_end > minus_start {
                        info.score -= plus_end - minus_start;
                        if plus_is_first {
                            info.trim_first = minus_start;
                        } else {
                            info.trim_second = minus_start;
                        }
                    }
                    if minus_end < plus_start {
                        info.score -= plus_start - minus_end;
                        if plus_is_first {
                            info.trim_second = plus_start;
                        } else {
                            info.trim_first = plus_start;
                        }
                    }

                    info.conf = PAIR_CONVERGENT;
                } else {
                    info.conf = PAIR_DIVERGENT;
                }

                if info.score < min_score {
                    continue;
                }
            } else {
                info.conf = PAIR_PARALLEL;
                info.score -= 1;
            }

            pair_info.push(info);
        }
    }

    if pair_info.is_empty() {
        *first_list = build_chain_list(first_chains);
        *second_list = build_chain_list(second_chains);
        return false;
    }

    pair_info.sort_by(|a, b| s_compare_pairs(a, b).cmp(&0));

    let best_score = pair_info[0].score;
    let margin = 5;
    let mut first_used = vec![false; first_chains.len()];
    let mut second_used = vec![false; second_chains.len()];
    let mut convergent_found = false;
    let mut found = false;

    for info in pair_info.iter_mut() {
        if first_used[info.first]
            || second_used[info.second]
            || (convergent_found && info.conf != PAIR_CONVERGENT)
        {
            continue;
        }
        if info.conf == PAIR_CONVERGENT {
            convergent_found = true;
        }
        if best_score - info.score > margin {
            break;
        }

        first_used[info.first] = true;
        second_used[info.second] = true;
        info.valid_pair = true;
        let first_score = first_chains[info.first].score;
        let second_score = second_chains[info.second].score;
        first_chains[info.first].pair = Some(info.second);
        second_chains[info.second].pair = Some(info.first);
        first_chains[info.first].pair_score = Some(second_score);
        second_chains[info.second].pair_score = Some(first_score);
        first_chains[info.first].pair_conf = info.conf;
        second_chains[info.second].pair_conf = info.conf;
        found = true;
    }

    if found {
        if let (Some(scoring_options), Some(query)) = (scoring_options, query) {
            for info in pair_info.iter().filter(|info| info.valid_pair) {
                if info.trim_first > 0 {
                    if hsp_frame_sign(&first_chains[info.first]) > 0 {
                        let _ = s_trim_chain_end_to_subj_pos(
                            Some(&mut first_chains[info.first]),
                            info.trim_first,
                            scoring_options.penalty,
                            scoring_options.gap_open,
                            scoring_options.gap_extend,
                            Some(query),
                        );
                    } else {
                        let _ = s_trim_chain_start_to_subj_pos(
                            Some(&mut first_chains[info.first]),
                            info.trim_first,
                            scoring_options.penalty,
                            scoring_options.gap_open,
                            scoring_options.gap_extend,
                            Some(query),
                        );
                    }
                }
                if info.trim_second > 0 {
                    if hsp_frame_sign(&second_chains[info.second]) > 0 {
                        let _ = s_trim_chain_end_to_subj_pos(
                            Some(&mut second_chains[info.second]),
                            info.trim_second,
                            scoring_options.penalty,
                            scoring_options.gap_open,
                            scoring_options.gap_extend,
                            Some(query),
                        );
                    } else {
                        let _ = s_trim_chain_start_to_subj_pos(
                            Some(&mut second_chains[info.second]),
                            info.trim_second,
                            scoring_options.penalty,
                            scoring_options.gap_open,
                            scoring_options.gap_extend,
                            Some(query),
                        );
                    }
                }
            }
        }
    }

    *first_list = build_chain_list(first_chains);
    *second_list = build_chain_list(second_chains);
    found
}

fn s_mapper_cutoff_score(params: &BlastHSPMapperParams, query_len: i32) -> i32 {
    if params.cutoff_score_fun[1] != 0 {
        (params.cutoff_score_fun[0] + params.cutoff_score_fun[1] * query_len) / 100
    } else if params.cutoff_score == 0 {
        get_cutoff_score(query_len)
    } else {
        params.cutoff_score
    }
}

fn chain_from_hsp(mut hsp: Hsp, oid: i32, query_info: &QueryInfo) -> Option<Box<HSPChain>> {
    if hsp.context < 0 {
        return None;
    }
    if hsp.query_frame == 0 {
        hsp.query_frame = query_info
            .contexts
            .get(hsp.context as usize)
            .map_or(0, |context| context.frame);
    }

    let context = hsp.context;
    let score = hsp.score;
    let mut hsp_slot = Some(hsp);
    let mut chain = hsp_chain_new(context);
    chain.oid = oid;
    chain.score = score;
    chain.count = 1;
    chain.hsps = hsp_container_new(&mut hsp_slot);
    Some(chain)
}

fn hsp_chain_list_append_splice_candidate(
    list: &mut Option<Box<HSPChain>>,
    hsp: Hsp,
    oid: i32,
    query_info: &QueryInfo,
    cutoff_score: i32,
    cutoff_edit_dist: i32,
) -> i32 {
    let Some(mut chain) = chain_from_hsp(hsp, oid, query_info) else {
        return 0;
    };
    if cutoff_score > 0 && !s_test_cutoffs(Some(&chain), cutoff_score, cutoff_edit_dist) {
        return 0;
    }

    let Some(new_hsp) = first_hsp(&chain) else {
        return -1;
    };
    let new_query_offset = new_hsp.query_offset;
    let new_subject_offset = new_hsp.subject_offset;
    let new_query_frame = new_hsp.query_frame;
    let new_context = chain.context;
    let new_oid = chain.oid;

    let mut chains = drain_chain_list(list.take());
    for existing in chains.iter_mut().rev() {
        let Some(existing_first) = first_hsp(existing) else {
            continue;
        };
        let Some(existing_last) = last_container(existing) else {
            continue;
        };
        if existing.oid == new_oid
            && existing.context == new_context
            && existing_first.query_frame == new_query_frame
            && new_query_offset >= existing_last.hsp.query_end
            && new_subject_offset >= existing_last.hsp.subject_end
        {
            let Some(new_container) = chain.hsps.take() else {
                *list = build_chain_list(chains);
                return -1;
            };
            let Some(last) = last_container_mut(existing) else {
                *list = build_chain_list(chains);
                return -1;
            };
            last.next = Some(new_container);
            existing.score = existing.score.saturating_add(chain.score);
            *list = build_chain_list(chains);
            return 0;
        }
    }

    chains.push(chain);
    *list = build_chain_list(chains);
    0
}

/// Conservative Rust port of NCBI `s_BlastHSPMapperSplicedPairedRun`
/// (`hspfilter_mapper.c:3864`).
///
/// Audit caveat: C builds optimal paths with `s_FindBestPath` and mutates raw
/// `BlastHSPList`/`HSPChain` ownership. Rust receives an owned `HspList`,
/// groups same-query/same-subject splice candidates into representable
/// multi-HSP chains, uses `ContextInfo::segment_flags` to identify paired read
/// segments, and then reuses the audited chain insert, splice-junction, pairing,
/// and trim helpers.
pub fn s_blast_hsp_mapper_spliced_paired_run(
    data: &mut BlastHSPMapperData,
    hsp_list: Option<HspList>,
) -> i32 {
    let Some(mut hsp_list) = hsp_list else {
        return 0;
    };
    let Some(params) = data.params.clone() else {
        return -1;
    };
    let Some(query_info) = data.query_info.clone() else {
        return -1;
    };
    let query_owned = data.query.clone();
    let query = query_owned.as_deref();
    if params.splice && query.is_none() {
        return -1;
    }

    let num_queries = query_info.num_queries.max(0) as usize;
    if data.saved_chains.len() < num_queries {
        data.saved_chains.resize_with(num_queries, || None);
    }
    if hsp_list.hsps.is_empty() {
        return 0;
    }

    if params.splice {
        hsp_list.hsps.sort_by(|a, b| {
            a.context
                .cmp(&b.context)
                .then(a.subject_offset.cmp(&b.subject_offset))
        });
    } else {
        hsp_list
            .hsps
            .sort_by(|a, b| a.context.cmp(&b.context).then(b.score.cmp(&a.score)));
    }

    let mut per_query: Vec<Option<Box<HSPChain>>> = (0..num_queries).map(|_| None).collect();
    let mut previous_key: Option<(i32, i32, i32)> = None;
    for hsp in hsp_list.hsps {
        let key = (hsp.context, hsp.subject_offset, hsp.score);
        if previous_key == Some(key) {
            continue;
        }
        previous_key = Some(key);
        if hsp.context < 0 {
            continue;
        }

        let Some(context_info) = query_info.contexts.get(hsp.context as usize) else {
            continue;
        };
        let query_idx = (hsp.context as usize) / NUM_STRANDS;
        if query_idx >= num_queries {
            continue;
        }

        let cutoff_score = s_mapper_cutoff_score(&params, context_info.query_length);
        let status = if params.splice {
            hsp_chain_list_append_splice_candidate(
                &mut per_query[query_idx],
                hsp,
                hsp_list.oid,
                &query_info,
                cutoff_score,
                params.cutoff_edit_dist,
            )
        } else {
            let mut chain = chain_from_hsp(hsp, hsp_list.oid, &query_info);
            hsp_chain_list_insert(
                &mut per_query[query_idx],
                &mut chain,
                cutoff_score,
                params.cutoff_edit_dist,
                true,
            )
        };
        if status != 0 {
            return status;
        }
    }

    if params.splice {
        for list in per_query.iter_mut().filter_map(Option::as_deref_mut) {
            let context = list.context.max(0) as usize;
            let query_len = query_info
                .contexts
                .get(context)
                .map_or(0, |context| context.query_length);
            let status = s_find_splice_junctions(
                Some(list),
                query,
                query_len,
                Some(&params.scoring_options),
            );
            if status != 0 {
                return status;
            }
        }
    }

    let pair_bonus = if params.splice { 21 } else { 5 };
    let mut pair_info = Vec::new();
    let mut query_idx = 0;
    while query_idx < num_queries {
        let has_pair = query_info.query_segment_flags(query_idx) == E_FIRST_SEGMENT;
        let fragment_end = if has_pair {
            (query_idx + 2).min(num_queries)
        } else {
            query_idx + 1
        };

        if has_pair && query_idx + 1 < num_queries {
            let (left, right) = per_query.split_at_mut(query_idx + 1);
            let first = &mut left[query_idx];
            let second = &mut right[0];
            if first.is_some() && second.is_some() {
                let _ = s_find_best_pairs(
                    first,
                    second,
                    0,
                    &mut pair_info,
                    params.splice,
                    Some(&params.scoring_options),
                    query,
                );
            }
        }

        for idx in query_idx..fragment_end {
            if per_query[idx].is_none() {
                continue;
            }
            let status = hsp_chain_list_insert(
                &mut data.saved_chains[idx],
                &mut per_query[idx],
                0,
                params.cutoff_edit_dist,
                true,
            );
            if status != 0 {
                return status;
            }
            if data.saved_chains[idx].is_some() {
                let _ = hsp_chain_list_trim(&mut data.saved_chains[idx], pair_bonus);
            }
        }

        query_idx = fragment_end;
    }

    0
}

/// Port of NCBI `s_FindRearrangedPairs` (`hspfilter_mapper.c:1865`).
///
/// Audit caveat: C stores direct `HSPChain*` pair links. Linked Rust ownership
/// stores the paired query index and cached partner score instead.
pub fn s_find_rearranged_pairs(saved: &mut [Option<Box<HSPChain>>], query_info: &QueryInfo) -> i32 {
    let num_queries = query_info.num_queries.max(0) as usize;
    for query_idx in 0..num_queries.saturating_sub(1) {
        if query_idx + 1 >= saved.len() {
            break;
        }
        if query_info.query_segment_flags(query_idx) != E_FIRST_SEGMENT {
            continue;
        }

        let (left, right) = saved.split_at_mut(query_idx + 1);
        let Some(chain) = left[query_idx].as_deref_mut() else {
            continue;
        };
        let Some(thepair) = right[0].as_deref_mut() else {
            continue;
        };

        if chain.next.is_some() || chain.pair.is_some() || thepair.next.is_some() {
            continue;
        }
        if chain.oid == thepair.oid {
            continue;
        }
        if hsp_frame_sign(chain) == hsp_frame_sign(thepair) {
            continue;
        }

        chain.pair = Some(query_idx + 1);
        thepair.pair = Some(query_idx);
        chain.pair_score = Some(thepair.score);
        thepair.pair_score = Some(chain.score);
        chain.pair_conf = PAIR_PARALLEL;
        thepair.pair_conf = PAIR_PARALLEL;
    }

    0
}

fn chain_partner_score(chain: &HSPChain, head_scores: &[Option<i32>]) -> Option<i32> {
    chain.pair_score.or_else(|| {
        chain
            .pair
            .and_then(|idx| head_scores.get(idx).copied().flatten())
    })
}

fn chain_prune_effective_score(chain: &HSPChain, pair_bonus: i32) -> i32 {
    if chain.pair.is_some() || chain.pair_score.is_some() {
        chain.score + pair_bonus
    } else {
        chain.score
    }
}

fn clear_chain_pair(chain: &mut HSPChain) {
    chain.pair = None;
    chain.pair_score = None;
    chain.pair_conf = PAIR_NONE;
}

fn has_represented_mate(
    saved: &[Option<Box<HSPChain>>],
    self_query_idx: usize,
    score: i32,
    mate_query_idx: Option<usize>,
    mate_score: Option<i32>,
) -> bool {
    let Some(mate_score) = mate_score else {
        return mate_query_idx
            .and_then(|idx| saved.get(idx))
            .and_then(Option::as_deref)
            .is_some();
    };

    saved.iter().any(|list| {
        let mut current = list.as_deref();
        while let Some(chain) = current {
            if chain.score == mate_score
                && (chain.pair_score == Some(score) || chain.pair == Some(self_query_idx))
            {
                return true;
            }
            current = chain.next.as_deref();
        }
        false
    })
}

fn reconcile_pruned_pairs(saved: &mut [Option<Box<HSPChain>>], num_queries: usize) {
    let snapshot = saved.to_vec();
    for query_idx in 0..num_queries.min(saved.len()) {
        let mut current = saved[query_idx].as_deref_mut();
        while let Some(chain) = current {
            if chain.pair.is_some() || chain.pair_score.is_some() {
                let mate_query_idx = chain.pair.filter(|&idx| idx < snapshot.len());
                if !has_represented_mate(
                    &snapshot,
                    query_idx,
                    chain.score,
                    mate_query_idx,
                    chain.pair_score,
                ) {
                    clear_chain_pair(chain);
                }
            }
            current = chain.next.as_deref_mut();
        }
    }
}

/// Port of NCBI `s_PruneChains` (`hspfilter_mapper.c:1937`).
///
/// Audit caveat: C stores direct pair pointers to specific `HSPChain` nodes.
/// Rust stores a paired query index plus cached partner score and reconciles
/// those cached links after pruning, matching `HSPChainFree` nulling the mate
/// pointer in C for represented paired-chain cases.
pub fn s_prune_chains(
    saved: &mut [Option<Box<HSPChain>>],
    num_queries: i32,
    pair_bonus: i32,
) -> i32 {
    let num_queries = (num_queries.max(0) as usize).min(saved.len());
    for list in saved.iter_mut().take(num_queries) {
        let mut chains = drain_chain_list(list.take());
        if chains.len() <= 1 {
            *list = build_chain_list(chains);
            continue;
        }

        let best_score = chains.iter().map(|chain| chain.score).max().unwrap_or(0);
        chains.retain(|chain| chain_prune_effective_score(chain, pair_bonus) >= best_score);
        *list = build_chain_list(chains);
    }

    reconcile_pruned_pairs(saved, num_queries);

    let head_scores: Vec<Option<i32>> = saved
        .iter()
        .map(|list| list.as_ref().map(|chain| chain.score))
        .collect();

    for list in saved.iter_mut().take(num_queries) {
        let mut chains = drain_chain_list(list.take());
        if chains.is_empty() {
            *list = None;
            continue;
        }
        if chains.len() == 1 {
            chains[0].count = 1;
            *list = build_chain_list(chains);
            continue;
        }

        let best_score = chains
            .iter()
            .map(|chain| chain_prune_effective_score(chain, pair_bonus))
            .max()
            .unwrap_or(0);
        let best_pair_score = chains
            .iter()
            .filter_map(|chain| chain_partner_score(chain, &head_scores).map(|p| chain.score + p))
            .max()
            .unwrap_or(0);

        let mut score_order: Vec<(usize, i32)> = chains
            .iter()
            .enumerate()
            .map(|(idx, chain)| (idx, chain.score))
            .collect();
        score_order.sort_by(|a, b| b.1.cmp(&a.1));

        let mut count = 0;
        let mut i = 0;
        while i < score_order.len() {
            count += 1;
            let score = score_order[i].1;
            let mut k = i + 1;
            while k < score_order.len() && score_order[k].1 == score {
                count += 1;
                k += 1;
            }
            for &(idx, _) in &score_order[i..k] {
                chains[idx].count = count;
            }
            i = k;
        }

        chains.retain(|chain| {
            if let Some(partner_score) = chain_partner_score(chain, &head_scores) {
                chain.score + partner_score >= best_pair_score
            } else {
                chain_prune_effective_score(chain, pair_bonus) >= best_score
            }
        });
        *list = build_chain_list(chains);
    }

    0
}

/// Port of NCBI `s_RemoveOverlaps` (`hspfilter_mapper.c:2203`).
///
/// Audit caveat: C takes a raw `HSPChain*` and splices new chains directly
/// after it. Rust takes the list head containing that chain, splits adjacent
/// query-overlapping HSPs into sibling chains, and reattaches the original
/// tail. Mapper-only pair pointers are cleared on generated sibling chains.
pub fn s_remove_overlaps(
    chains: &mut Option<Box<HSPChain>>,
    score_opts: Option<&ScoringOptions>,
    query_len: i32,
) -> i32 {
    let Some(score_opts) = score_opts else {
        return -1;
    };
    let mut chain_list = drain_chain_list(chains.take());
    if chain_list.is_empty() {
        *chains = None;
        return -1;
    }

    let mut head = chain_list.remove(0);
    let tail = chain_list;
    let hsps = drain_hsp_container_list(head.hsps.take());
    if hsps.is_empty() {
        head.hsps = None;
        let mut rebuilt = vec![head];
        rebuilt.extend(tail);
        *chains = build_chain_list(rebuilt);
        return 0;
    }

    let mut chunks: Vec<Vec<Box<HSPContainer>>> = Vec::new();
    let mut current_chunk: Vec<Box<HSPContainer>> = Vec::new();
    for container in hsps {
        let starts_new_chunk = current_chunk
            .last()
            .map(|prev| prev.hsp.query_end > container.hsp.query_offset)
            .unwrap_or(false);
        if starts_new_chunk {
            chunks.push(current_chunk);
            current_chunk = Vec::new();
        }
        current_chunk.push(container);
    }
    if !current_chunk.is_empty() {
        chunks.push(current_chunk);
    }

    let rescore = chunks.len() > 1;
    let mut rebuilt = Vec::new();
    for (idx, chunk) in chunks.into_iter().enumerate() {
        let mut split_chain = if idx == 0 {
            (*head).clone()
        } else {
            let mut cloned = (*head).clone();
            cloned.pair = None;
            cloned.pair_score = None;
            cloned.next = None;
            cloned
        };
        split_chain.hsps = build_hsp_container_list(chunk);
        if rescore {
            let _ =
                s_compute_chain_score(Some(&mut split_chain), Some(score_opts), query_len, false);
        }
        rebuilt.push(Box::new(split_chain));
    }
    rebuilt.extend(tail);
    *chains = build_chain_list(rebuilt);
    0
}

fn sort_chain_list_by_oid(list: &mut Option<Box<HSPChain>>) {
    let mut chains = drain_chain_list(list.take());
    chains.sort_by(|a, b| s_compare_chains_by_oid(a, b).cmp(&0));
    *list = build_chain_list(chains);
}

fn remove_overlaps_from_chain_list(
    list: &mut Option<Box<HSPChain>>,
    score_opts: &ScoringOptions,
    query_len: i32,
) -> i32 {
    let chains = drain_chain_list(list.take());
    let mut rebuilt = Vec::new();
    for chain in chains {
        let mut single = Some(chain);
        let status = s_remove_overlaps(&mut single, Some(score_opts), query_len);
        if status != 0 {
            *list = build_chain_list(rebuilt);
            return status;
        }
        rebuilt.extend(drain_chain_list(single));
    }
    *list = build_chain_list(rebuilt);
    0
}

fn set_unique_mapping_counts(list: &mut Option<Box<HSPChain>>) {
    let mut chains = drain_chain_list(list.take());
    if chains.is_empty() {
        *list = None;
        return;
    }

    let mut num_unique = 1;
    for idx in 1..chains.len() {
        if chains[idx - 1].oid != chains[idx].oid
            || s_find_fragment_start(Some(&chains[idx - 1]))
                != s_find_fragment_start(Some(&chains[idx]))
        {
            num_unique += 1;
        }
    }
    for chain in &mut chains {
        chain.count = num_unique;
    }
    *list = build_chain_list(chains);
}

/// Port of NCBI `s_Finalize` (`hspfilter_mapper.c:2295`).
///
/// Audit caveat: C transfers the raw `saved` array into `results`. Rust moves
/// `saved` into `results.chain_array` with `mem::take`. Environment-gated C
/// behavior is preserved for pruning and overlap splitting; pair links use the
/// Rust paired-index/cache representation documented on the helper ports.
pub fn s_finalize(
    saved: &mut Vec<Option<Box<HSPChain>>>,
    results: &mut BlastMappingResults,
    query_info: &QueryInfo,
    query_sequence: &[u8],
    score_opts: &ScoringOptions,
    _is_paired: bool,
    cutoff_score: i32,
    cutoff_edit_dist: i32,
) -> i32 {
    const K_PAIR_BONUS: i32 = 21;
    let num_queries = (query_info.num_queries.max(0) as usize).min(saved.len());

    if std::env::var_os("MAPPER_NO_PRUNNING").is_none() {
        let status = s_prune_chains(saved, query_info.num_queries, K_PAIR_BONUS);
        if status != 0 {
            return status;
        }
    }

    let status = s_find_adapters(saved, query_sequence, query_info, score_opts);
    if status != 0 {
        return status;
    }
    let status = s_find_poly_a_tails(saved, query_sequence, query_info);
    if status != 0 {
        return status;
    }
    let status = s_find_rearranged_pairs(saved, query_info);
    if status != 0 {
        return status;
    }

    for list in saved.iter_mut().take(num_queries) {
        sort_chain_list_by_oid(list);
    }

    if std::env::var_os("MAPPER_NO_OVERLAPPED_HSP_MERGE").is_some() {
        for query_idx in 0..num_queries {
            let query_len = saved[query_idx]
                .as_ref()
                .and_then(|chain| query_info.contexts.get(chain.context.max(0) as usize))
                .map(|ctx| ctx.query_length)
                .unwrap_or(0);
            let status =
                remove_overlaps_from_chain_list(&mut saved[query_idx], score_opts, query_len);
            if status != 0 {
                return status;
            }
        }
    }

    for list in saved.iter_mut().take(num_queries) {
        let status = s_filter_chains(list, cutoff_score, cutoff_edit_dist);
        if status != 0 {
            return status;
        }
    }

    for list in saved.iter_mut().take(num_queries) {
        set_unique_mapping_counts(list);
    }

    results.num_queries = query_info.num_queries;
    results.chain_array = std::mem::take(saved);
    0
}

/// Port of NCBI `s_BlastHSPMapperFinal` (`hspfilter_mapper.c:2402`).
///
/// Audit caveat: C receives opaque `void*` writer data and mapping-result
/// pointers. Rust uses typed references and returns `-1` when required mapper
/// fields are absent instead of dereferencing null pointers.
pub fn s_blast_hsp_mapper_final(
    data: &mut BlastHSPMapperData,
    results: &mut BlastMappingResults,
) -> i32 {
    if data.saved_chains.is_empty() {
        return 0;
    }

    let Some(params) = data.params.clone() else {
        return -1;
    };
    let Some(query_info) = data.query_info.clone() else {
        return -1;
    };
    let Some(query) = data.query.as_deref() else {
        return -1;
    };

    s_finalize(
        &mut data.saved_chains,
        results,
        &query_info,
        query,
        &params.scoring_options,
        params.paired,
        params.cutoff_score,
        params.cutoff_edit_dist,
    )
}

fn mapper_splice_signal_at(subject: &[u8], i: usize) -> Option<u8> {
    let a = *subject.get(i)?;
    let b = *subject.get(i + 1)?;
    let c = *subject.get(i + 2)?;
    let d = *subject.get(i + 3)?;
    Some((a << 2) | b | (c << 6) | (d << 4))
}

fn mapper_set_splice_signal(first: &mut Hsp, second: &mut Hsp, signal: u8) {
    let first_info = first
        .map_info
        .get_or_insert_with(crate::hspstream::blast_hsp_mapping_info_new);
    first_info.right_edge = (signal & 0xf0) >> 4;
    first_info.right_edge |= MAPPER_SPLICE_SIGNAL;

    let second_info = second
        .map_info
        .get_or_insert_with(crate::hspstream::blast_hsp_mapping_info_new);
    second_info.left_edge = signal & 0x0f;
    second_info.left_edge |= MAPPER_SPLICE_SIGNAL;
}

fn mapper_clear_splice_signal(first: &mut Hsp, second: &mut Hsp) {
    if let Some(info) = first.map_info.as_mut() {
        info.right_edge &= !MAPPER_SPLICE_SIGNAL;
    }
    if let Some(info) = second.map_info.as_mut() {
        info.left_edge &= !MAPPER_SPLICE_SIGNAL;
    }
}

fn mapper_clear_left_splice_exon(hsp: &mut Hsp) {
    if let Some(info) = hsp.map_info.as_mut() {
        info.left_edge &= !MAPPER_SPLICE_SIGNAL;
        info.left_edge &= !MAPPER_EXON;
    }
}

fn mapper_clear_right_splice_exon(hsp: &mut Hsp) {
    if let Some(info) = hsp.map_info.as_mut() {
        info.right_edge &= !MAPPER_SPLICE_SIGNAL;
        info.right_edge &= !MAPPER_EXON;
    }
}

fn mapper_mark_left_edge(hsp: &mut Hsp, flags: u8) {
    let info = hsp
        .map_info
        .get_or_insert_with(crate::hspstream::blast_hsp_mapping_info_new);
    info.left_edge |= flags;
}

fn mapper_mark_right_edge(hsp: &mut Hsp, flags: u8) {
    let info = hsp
        .map_info
        .get_or_insert_with(crate::hspstream::blast_hsp_mapping_info_new);
    info.right_edge |= flags;
}

fn mapper_overlap_edit_counts(first: &Hsp, second: &Hsp) -> (usize, usize) {
    let first_count = first
        .map_info
        .as_ref()
        .and_then(|info| info.edits.as_ref())
        .map(|edits| {
            edits
                .edits
                .iter()
                .rev()
                .take_while(|edit| edit.query_pos >= second.query_offset)
                .count()
        })
        .unwrap_or(0);
    let second_count = second
        .map_info
        .as_ref()
        .and_then(|info| info.edits.as_ref())
        .map(|edits| {
            edits
                .edits
                .iter()
                .take_while(|edit| edit.query_pos < first.query_end)
                .count()
        })
        .unwrap_or(0);
    (first_count, second_count)
}

fn mapper_overlap_edits_remain(first: &Hsp, second: &Hsp) -> bool {
    let first_remains = first
        .map_info
        .as_ref()
        .and_then(|info| info.edits.as_ref())
        .and_then(|edits| edits.edits.last())
        .is_some_and(|edit| edit.query_pos >= second.query_offset);
    let second_remains = second
        .map_info
        .as_ref()
        .and_then(|info| info.edits.as_ref())
        .and_then(|edits| edits.edits.first())
        .is_some_and(|edit| edit.query_pos < first.query_end);
    first_remains || second_remains
}

fn mapper_trim_first_overlap_edits(
    first: &mut Hsp,
    second_query_offset: i32,
    query: &[u8],
    query_len: i32,
) -> i32 {
    const GAP_BASE: u8 = 15;
    loop {
        let edit = first
            .map_info
            .as_ref()
            .and_then(|info| info.edits.as_ref())
            .and_then(|edits| edits.edits.last())
            .copied();
        let Some(edit) = edit else {
            return 0;
        };
        if edit.query_pos < second_query_offset {
            return 0;
        }

        let mut edge = first
            .map_info
            .as_ref()
            .map(|info| info.right_edge)
            .unwrap_or_default();
        if edit.query_pos >= first.query_end.saturating_sub(1) {
            if edit.subject_base != GAP_BASE {
                edge >>= 2;
                edge |= (edit.subject_base as u8) << 2;
            }
        } else if edit.query_pos == query_len.saturating_sub(2) && edit.subject_base == GAP_BASE {
            if let Some(&base) = query.get(query_len.saturating_sub(1).max(0) as usize) {
                edge = (edge << 2) | base;
            }
        } else if edit.subject_base != GAP_BASE && edit.query_base != GAP_BASE {
            if let Some(&base) = query.get(edit.query_pos.saturating_add(1).max(0) as usize) {
                edge = ((edit.subject_base as u8) << 2) | base;
            }
        } else if edit.subject_base == GAP_BASE {
            if let (Some(&left), Some(&right)) = (
                query.get(edit.query_pos.saturating_add(1).max(0) as usize),
                query.get(edit.query_pos.saturating_add(2).max(0) as usize),
            ) {
                edge = (left << 2) | right;
            }
        } else {
            edge = ((edit.subject_base as u8) << 2)
                | query
                    .get(edit.query_pos.max(0) as usize)
                    .copied()
                    .unwrap_or_default();
        }

        let mut trim_by = first.query_end.saturating_sub(edit.query_pos);
        let status = if edit.query_base == GAP_BASE {
            trim_by = trim_by.saturating_add(1);
            s_trim_hsp(first, trim_by, false, false, -4, -4, -4, Some(query))
        } else {
            s_trim_hsp(first, trim_by, true, false, -4, -4, -4, Some(query))
        };
        if status != 0 {
            return status;
        }
        if let Some(info) = first.map_info.as_mut() {
            info.right_edge = edge;
        }
    }
}

fn mapper_trim_second_overlap_edits(second: &mut Hsp, first_query_end: i32, query: &[u8]) -> i32 {
    const GAP_BASE: u8 = 15;
    loop {
        let edit = second
            .map_info
            .as_ref()
            .and_then(|info| info.edits.as_ref())
            .and_then(|edits| edits.edits.first())
            .copied();
        let Some(edit) = edit else {
            return 0;
        };
        if edit.query_pos >= first_query_end {
            return 0;
        }

        let mut edge = second
            .map_info
            .as_ref()
            .map(|info| info.left_edge)
            .unwrap_or_default();
        if edit.query_pos == 0 {
            if edit.subject_base != GAP_BASE {
                edge <<= 2;
                edge |= edit.subject_base as u8;
            }
        } else if edit.query_pos == 1 && edit.subject_base == GAP_BASE {
            if let Some(&base) = query.first() {
                edge = (edge << 2) | base;
            }
        } else if edit.subject_base == GAP_BASE {
            if let (Some(&left), Some(&right)) = (
                query.get(edit.query_pos.saturating_sub(2).max(0) as usize),
                query.get(edit.query_pos.saturating_sub(1).max(0) as usize),
            ) {
                edge = (left << 2) | right;
            }
        } else {
            edge = (query
                .get(edit.query_pos.saturating_sub(1).max(0) as usize)
                .copied()
                .unwrap_or_default()
                << 2)
                | edit.subject_base as u8;
        }

        let trim_by = edit
            .query_pos
            .saturating_sub(second.query_offset)
            .saturating_add(1);
        let status = if edit.query_base == GAP_BASE {
            s_trim_hsp(second, trim_by, false, true, -4, -4, -4, Some(query))
        } else {
            s_trim_hsp(second, trim_by, true, true, -4, -4, -4, Some(query))
        };
        if status != 0 {
            return status;
        }
        if let Some(info) = second.map_info.as_mut() {
            info.left_edge = edge;
        }
    }
}

/// Conservative Rust port of NCBI `s_FindSpliceJunctionsForOverlaps`
/// (`hspfilter_mapper.c:2504`).
///
/// This ports C's coordinate guards, splice-signal search over the reconstructed
/// overlap when query flanks are available, mapper splice-edge updates, exact
/// split trimming at the signal, and the edit-aware overlap pre-split. The
/// caller performs overlap-removal fallback.
pub fn s_find_splice_junctions_for_overlaps(
    first: &mut Hsp,
    second: &mut Hsp,
    query: Option<&[u8]>,
    query_len: i32,
    consensus_only: bool,
) -> i32 {
    const SIGNALS: [u8; 23] = [
        0xb2, 0x71, 0x92, 0x79, 0x31, 0xb3, 0xbe, 0x41, 0xba, 0x51, 0xb0, 0xf1, 0x82, 0x7d, 0xf2,
        0x70, 0x32, 0x73, 0xa2, 0x75, 0x33, 0x30, 0xf3,
    ];

    if first.query_offset >= second.query_offset || first.subject_offset >= second.subject_offset {
        return -1;
    }

    let overlap_len = first.query_end.saturating_sub(second.query_offset);
    if overlap_len <= 0 {
        return 0;
    }

    let first_query_span = first.query_end.saturating_sub(first.query_offset);
    let second_query_span = second.query_end.saturating_sub(second.query_offset);
    if overlap_len >= first_query_span || overlap_len >= second_query_span {
        return -1;
    }

    if let Some(query) = query {
        let query_len = query_len.max(0).min(query.len() as i32);
        let (first_edits, second_edits) = mapper_overlap_edit_counts(first, second);
        if first_edits > 0 || second_edits > 0 {
            if first_edits > second_edits {
                let status =
                    mapper_trim_first_overlap_edits(first, second.query_offset, query, query_len);
                if status != 0 {
                    return status;
                }
            } else if second_edits > first_edits {
                let status = mapper_trim_second_overlap_edits(second, first.query_end, query);
                if status != 0 {
                    return status;
                }
            } else {
                mapper_clear_splice_signal(first, second);
                return 0;
            }

            if mapper_overlap_edits_remain(first, second) {
                mapper_clear_splice_signal(first, second);
                return 0;
            }
        }

        let overlap_len = first.query_end.saturating_sub(second.query_offset);
        if overlap_len <= 0 {
            mapper_clear_splice_signal(first, second);
            return 0;
        }
        if second.query_offset >= 0 && first.query_end <= query_len {
            if let (Some(first_right_edge), Some(second_left_edge)) = (
                first.map_info.as_ref().map(|info| info.right_edge & 0x0f),
                second.map_info.as_ref().map(|info| info.left_edge & 0x0f),
            ) {
                let mut subject = Vec::with_capacity(overlap_len as usize + 4);
                subject.push((second_left_edge >> 2) & 0x03);
                subject.push(second_left_edge & 0x03);
                subject.extend_from_slice(
                    &query[second.query_offset as usize..first.query_end as usize],
                );
                subject.push((first_right_edge >> 2) & 0x03);
                subject.push(first_right_edge & 0x03);

                let num_signals = if consensus_only { 2 } else { SIGNALS.len() };
                for &signal in SIGNALS.iter().take(num_signals) {
                    for i in 0..=overlap_len as usize {
                        if mapper_splice_signal_at(&subject, i) != Some(signal) {
                            continue;
                        }

                        let trim_first = first
                            .query_end
                            .saturating_sub(second.query_offset.saturating_add(i as i32));
                        if trim_first > 0 {
                            let status =
                                s_trim_hsp(first, trim_first, true, false, -4, -4, -4, Some(query));
                            if status != 0 {
                                return status;
                            }
                        }
                        if i > 0 {
                            let status =
                                s_trim_hsp(second, i as i32, true, true, -4, -4, -4, Some(query));
                            if status != 0 {
                                return status;
                            }
                        }
                        mapper_set_splice_signal(first, second, signal);
                        return 0;
                    }
                }
            }
        }
    }

    mapper_clear_splice_signal(first, second);
    0
}

/// Rust ownership equivalent of NCBI `s_ExtendAlignmentCleanup`
/// (`hspfilter_mapper.c:2779`).
///
/// Audit caveat: C frees a `BlastGapAlignStruct` and `JumperEditsBlock`; those
/// mapper-specific structs are not represented in `spliced_hits`, so callers
/// pass any owned Rust payloads that should be dropped at the cleanup point.
pub fn s_extend_alignment_cleanup<G, E>(
    subject: Option<Vec<u8>>,
    gap_align: Option<G>,
    edit_script: Option<crate::gapinfo::GapEditScript>,
    edits: Option<E>,
) {
    drop((subject, gap_align, edit_script, edits));
}

fn mapper_append_extension_ops(
    script: &mut crate::gapinfo::GapEditScript,
    query_len: i32,
    subject_len: i32,
) {
    let matches = query_len.min(subject_len);
    if matches > 0 {
        script.push(GapAlignOpType::Sub, matches);
    }
    if query_len > subject_len {
        script.push(GapAlignOpType::Ins, query_len.saturating_sub(subject_len));
    } else if subject_len > query_len {
        script.push(GapAlignOpType::Del, subject_len.saturating_sub(query_len));
    }
}

fn mapper_extend_hsp(
    hsp: &mut Hsp,
    query_len: i32,
    subject_len: i32,
    is_left: bool,
    score_opts: &ScoringOptions,
) -> i32 {
    if query_len < 0 || subject_len < 0 {
        return -1;
    }
    if query_len == 0 && subject_len == 0 {
        return 0;
    }

    let mut extension = crate::gapinfo::GapEditScript::new();
    mapper_append_extension_ops(&mut extension, query_len, subject_len);

    let mut script = hsp
        .edit_script
        .clone()
        .unwrap_or_else(crate::gapinfo::GapEditScript::new);
    if script.is_empty() && hsp.query_end > hsp.query_offset {
        script.push(
            GapAlignOpType::Sub,
            hsp.query_end
                .saturating_sub(hsp.query_offset)
                .min(hsp.subject_end.saturating_sub(hsp.subject_offset)),
        );
    }

    if is_left {
        for (op, count) in script.iter() {
            extension.push(op, count);
        }
        hsp.edit_script = Some(extension);
        hsp.query_offset = hsp.query_offset.saturating_sub(query_len);
        hsp.subject_offset = hsp.subject_offset.saturating_sub(subject_len);
        hsp.query_gapped_start = hsp.query_offset;
        hsp.subject_gapped_start = hsp.subject_offset;
    } else {
        for (op, count) in extension.iter() {
            script.push(op, count);
        }
        hsp.edit_script = Some(script);
        hsp.query_end = hsp.query_end.saturating_add(query_len);
        hsp.subject_end = hsp.subject_end.saturating_add(subject_len);
    }

    hsp.score = mapper_compute_alignment_score(
        hsp,
        score_opts.penalty,
        score_opts.gap_open,
        score_opts.gap_extend,
    );
    hsp.num_gaps = hsp
        .edit_script
        .as_ref()
        .map(|script| {
            script
                .iter()
                .filter(|(op, _)| !matches!(op, GapAlignOpType::Sub | GapAlignOpType::Decline))
                .map(|(_, count)| count)
                .sum()
        })
        .unwrap_or(0);
    0
}

/// NCBI: s_ExtendAlignment (`hspfilter_mapper.c:2795`).
#[allow(clippy::too_many_arguments)]
pub fn s_extend_alignment(
    hsp: &mut Hsp,
    query: Option<&[u8]>,
    query_from: i32,
    query_to: i32,
    subject_from: i32,
    subject_to: i32,
    score_options: Option<&ScoringOptions>,
    is_left: bool,
) -> i32 {
    let Some(query) = query else {
        return -1;
    };
    let Some(score_options) = score_options else {
        return -1;
    };
    let Some(map_info) = hsp.map_info.as_ref() else {
        return -1;
    };
    let Some(overhangs) = map_info.subject_overhangs.as_ref() else {
        return -1;
    };
    let (o_len, o_seq) = if is_left {
        (overhangs.left_len, overhangs.left.as_deref())
    } else {
        (overhangs.right_len, overhangs.right.as_deref())
    };
    let Some(o_seq) = o_seq else {
        return -1;
    };
    if o_len <= 0 || o_len as usize > o_seq.len() {
        return -1;
    }
    if query_from < 0
        || query_to < query_from - 1
        || subject_from < 0
        || subject_to < subject_from - 1
    {
        return -1;
    }
    let query_gap = if query_to >= query_from {
        query_to - query_from + 1
    } else {
        0
    };
    let subject_gap = if subject_to >= subject_from {
        subject_to - subject_from + 1
    } else {
        0
    };
    if query_gap < 0
        || subject_gap < 0
        || query_from as usize > query.len()
        || query_from.saturating_add(query_gap) as usize > query.len()
        || subject_from >= o_len
        || subject_from.saturating_add(subject_gap) > o_len
    {
        return -1;
    }
    if query_gap == 0 && subject_gap == 0 {
        return 0;
    }

    let subject = crate::encoding::pack_ncbi2na_bases(&o_seq[..o_len as usize]);
    let Some(mut gap_align) = crate::gapinfo::jumper_gap_align_new(o_len.saturating_mul(2)) else {
        return -1;
    };
    let Some(right_prelim) = gap_align.right_prelim_block.as_mut() else {
        return -1;
    };

    let query_start = query_from as usize;
    let query_end = query_start + query_gap as usize;
    let subject_start = subject_from as usize;
    let subject_end = subject_start + subject_gap as usize;
    let mut num_identical = 0;
    let mut ungapped_ext_len = 0;
    let (status, mut query_ext_len, mut subject_ext_len) = if query_gap > 0 && subject_gap > 0 {
        crate::gapinfo::jumper_extend_right_with_traceback(
            &query[query_start..query_end],
            &o_seq[subject_start..subject_end],
            1,
            0,
            0,
            0,
            20,
            20,
            right_prelim,
            &mut num_identical,
            false,
            &mut ungapped_ext_len,
        )
    } else {
        (0, 0, 0)
    };
    if status < 0 || query_ext_len > query_gap || subject_ext_len > subject_gap {
        return -1;
    }

    while query_ext_len < query_gap {
        if crate::gapinfo::jumper_prelim_edit_block_add(
            right_prelim,
            crate::gapinfo::JUMPER_INSERTION,
        ) != 0
        {
            return -1;
        }
        query_ext_len += 1;
    }
    while subject_ext_len < subject_gap {
        if crate::gapinfo::jumper_prelim_edit_block_add(
            right_prelim,
            crate::gapinfo::JUMPER_DELETION,
        ) != 0
        {
            return -1;
        }
        subject_ext_len += 1;
    }

    let Some(left_prelim) = gap_align.left_prelim_block.as_ref() else {
        return -1;
    };
    let Some(extension_script) =
        crate::gapinfo::jumper_prelim_edit_block_to_gap_edit_script(left_prelim, right_prelim)
    else {
        return -1;
    };
    let Some(extension_edits) = crate::gapinfo::jumper_find_edits(
        query,
        &subject,
        query_from,
        subject_from,
        query_from.saturating_add(query_ext_len),
        subject_from.saturating_add(subject_ext_len),
        left_prelim,
        right_prelim,
    ) else {
        return -1;
    };

    if is_left {
        let mut prefix = Some(extension_script);
        let mut current = hsp.edit_script.take();
        hsp.edit_script = crate::gapinfo::gap_edit_script_combine(&mut prefix, &mut current);
    } else {
        let mut current = hsp.edit_script.take();
        let mut suffix = Some(extension_script);
        hsp.edit_script = crate::gapinfo::gap_edit_script_combine(&mut current, &mut suffix);
    }
    if hsp.edit_script.is_none() {
        return -1;
    }

    let Some(map_info) = hsp.map_info.as_mut() else {
        return -1;
    };
    if is_left {
        let mut prefix = Some(extension_edits);
        let mut current = map_info.edits.take();
        map_info.edits = if current.is_some() {
            crate::gapinfo::jumper_edits_block_combine(&mut prefix, &mut current)
        } else {
            prefix
        };
    } else {
        let mut current = map_info.edits.take();
        let mut suffix = Some(extension_edits);
        map_info.edits = if current.is_some() {
            crate::gapinfo::jumper_edits_block_combine(&mut current, &mut suffix)
        } else {
            suffix
        };
    }
    if map_info.edits.is_none() {
        return -1;
    }

    if is_left {
        hsp.query_offset = hsp.query_offset.saturating_sub(query_ext_len);
        hsp.subject_offset = hsp.subject_offset.saturating_sub(subject_ext_len);
        hsp.query_gapped_start = hsp.query_offset;
        hsp.subject_gapped_start = hsp.subject_offset;
    } else {
        hsp.query_end = hsp.query_end.saturating_add(query_ext_len);
        hsp.subject_end = hsp.subject_end.saturating_add(subject_ext_len);
    }

    hsp.score = mapper_compute_alignment_score(
        hsp,
        score_options.penalty,
        score_options.gap_open,
        score_options.gap_extend,
    );
    hsp.num_gaps = hsp
        .edit_script
        .as_ref()
        .map(|script| {
            script
                .iter()
                .filter(|(op, _)| !matches!(op, GapAlignOpType::Sub | GapAlignOpType::Decline))
                .map(|(_, count)| count)
                .sum()
        })
        .unwrap_or(0);
    0
}

fn mapper_query_base(query: &[u8], pos: i32, query_len: i32) -> Option<u8> {
    if pos < 0 || pos >= query_len {
        return None;
    }
    query.get(pos as usize).copied()
}

/// Variant of [`s_find_splice_junctions_for_gap`] that reads overhangs from
/// the C-shaped `hsp->map_info` field now represented on Rust HSPs.
/// blast-rs: Native helper adapting Rust HSP map-info storage to the C-shaped port.
pub fn s_find_splice_junctions_for_gap_using_map_info(
    first: &mut Hsp,
    second: &mut Hsp,
    query: Option<&[u8]>,
    query_len: i32,
    score_opts: Option<&ScoringOptions>,
) -> i32 {
    let first_right = first
        .map_info
        .as_ref()
        .and_then(|info| info.subject_overhangs.as_ref())
        .and_then(|overhangs| overhangs.right.clone());
    let second_left = second
        .map_info
        .as_ref()
        .and_then(|info| info.subject_overhangs.as_ref())
        .and_then(|overhangs| overhangs.left.clone());

    s_find_splice_junctions_for_gap(
        first,
        second,
        query,
        query_len,
        score_opts,
        first_right.as_deref(),
        second_left.as_deref(),
    )
}

/// Conservative Rust port of NCBI `s_FindSpliceJunctionsForGap`
/// (`hspfilter_mapper.c:3010`).
///
/// This ports the canonical splice-signal search plus coordinate/gap-script
/// extension. When a signal is found it updates mapper edge bits on
/// `Hsp::map_info`, matching C's `hsp->map_info` side effect.
pub fn s_find_splice_junctions_for_gap(
    first: &mut Hsp,
    second: &mut Hsp,
    query: Option<&[u8]>,
    query_len: i32,
    score_opts: Option<&ScoringOptions>,
    first_right_overhang: Option<&[u8]>,
    second_left_overhang: Option<&[u8]>,
) -> i32 {
    let Some(query) = query else {
        return -1;
    };
    let Some(score_opts) = score_opts else {
        return -1;
    };
    let Some(first_right) = first_right_overhang else {
        return -1;
    };
    let Some(second_left) = second_left_overhang else {
        return -1;
    };

    let query_len = query_len.max(0).min(query.len() as i32);
    let query_gap = second.query_offset.saturating_sub(first.query_end);
    if query_gap < 0 {
        return 0;
    }
    if query_gap > first_right.len() as i32 - 2 || query_gap > second_left.len() as i32 - 2 {
        return 0;
    }

    const SIGNALS: [u8; 2] = [0xb2, 0x71];
    let first_len = first_right.len() as i32;
    let second_len = second_left.len() as i32;

    for q in 0..4 {
        if first.query_end.saturating_sub(q) <= first.query_offset {
            break;
        }
        let Some(high) = (match q {
            0 => first_right
                .first()
                .zip(first_right.get(1))
                .map(|(a, b)| (*a << 6) | (*b << 4)),
            1 => mapper_query_base(query, first.query_end.saturating_sub(1), query_len)
                .zip(first_right.first().copied())
                .map(|(a, b)| (a << 6) | (b << 4)),
            _ => mapper_query_base(query, first.query_end.saturating_sub(q), query_len)
                .zip(mapper_query_base(
                    query,
                    first.query_end.saturating_sub(q).saturating_add(1),
                    query_len,
                ))
                .map(|(a, b)| (a << 6) | (b << 4)),
        }) else {
            continue;
        };

        for &signal in &SIGNALS {
            if high != (signal & 0xf0) {
                continue;
            }
            let start = second_len
                .saturating_sub(1)
                .saturating_sub(query_gap)
                .saturating_sub(2)
                .saturating_sub(q);
            let lo = start.saturating_sub(1).max(0);
            let hi = start.saturating_add(1).min(second_len.saturating_sub(2));
            for i in lo..=hi {
                let low = (second_left[i as usize] << 2) | second_left[i as usize + 1];
                if high | low != signal {
                    continue;
                }
                let subject_gap = second_len
                    .saturating_sub(i.saturating_add(2))
                    .saturating_add(q);
                if query_gap.abs_diff(subject_gap) > 1 {
                    continue;
                }
                if q > 0 {
                    let status = s_trim_hsp(
                        first,
                        q,
                        true,
                        false,
                        score_opts.penalty,
                        score_opts.gap_open,
                        score_opts.gap_extend,
                        Some(query),
                    );
                    if status != 0 {
                        return status;
                    }
                }
                let query_ext = second.query_offset.saturating_sub(first.query_end);
                let subject_ext = second_len.saturating_sub(i.saturating_add(2));
                {
                    let info = second
                        .map_info
                        .get_or_insert_with(crate::hspstream::blast_hsp_mapping_info_new);
                    let overhangs = info.subject_overhangs.get_or_insert_with(Default::default);
                    overhangs.left_len = second_left.len() as i32;
                    overhangs.left = Some(second_left.to_vec());
                }
                let status = s_extend_alignment(
                    second,
                    Some(query),
                    second.query_offset.saturating_sub(query_ext),
                    second.query_offset.saturating_sub(1),
                    0,
                    subject_ext.saturating_sub(1),
                    Some(score_opts),
                    true,
                );
                if status == 0 {
                    mapper_set_splice_signal(first, second, signal);
                }
                return status;
            }
        }
    }

    for q in 0..4 {
        if second.query_offset.saturating_add(q) >= second.query_end {
            break;
        }
        let Some(low) = (match q {
            0 => second_left
                .get(second_left.len().saturating_sub(2))
                .zip(second_left.last())
                .map(|(a, b)| (*a << 2) | *b),
            1 => second_left
                .last()
                .copied()
                .zip(mapper_query_base(query, second.query_offset, query_len))
                .map(|(a, b)| (a << 2) | b),
            _ => mapper_query_base(
                query,
                second.query_offset.saturating_add(q).saturating_sub(2),
                query_len,
            )
            .zip(mapper_query_base(
                query,
                second.query_offset.saturating_add(q).saturating_sub(1),
                query_len,
            ))
            .map(|(a, b)| (a << 2) | b),
        }) else {
            continue;
        };

        for &signal in &SIGNALS {
            if low != (signal & 0x0f) {
                continue;
            }
            if q > 1
                && second
                    .map_info
                    .as_ref()
                    .and_then(|info| info.edits.as_ref())
                    .and_then(|edits| edits.edits.first())
                    .is_some_and(|edit| edit.query_pos < q - 2)
            {
                break;
            }
            let end = query_gap.saturating_add(q);
            let lo = end.saturating_sub(1).max(0);
            let hi = end.saturating_add(1).min(first_len.saturating_sub(2));
            for i in lo..=hi {
                let high = (first_right[i as usize] << 6) | (first_right[i as usize + 1] << 4);
                if high | low != signal {
                    continue;
                }
                let subject_gap = i.saturating_sub(q);
                if query_gap.abs_diff(subject_gap) > 1 {
                    continue;
                }
                if q > 0 {
                    let status = s_trim_hsp(
                        second,
                        q,
                        true,
                        true,
                        score_opts.penalty,
                        score_opts.gap_open,
                        score_opts.gap_extend,
                        Some(query),
                    );
                    if status != 0 {
                        return status;
                    }
                }
                let query_ext = second.query_offset.saturating_sub(first.query_end);
                {
                    let info = first
                        .map_info
                        .get_or_insert_with(crate::hspstream::blast_hsp_mapping_info_new);
                    let overhangs = info.subject_overhangs.get_or_insert_with(Default::default);
                    overhangs.right_len = first_right.len() as i32;
                    overhangs.right = Some(first_right.to_vec());
                }
                let status = s_extend_alignment(
                    first,
                    Some(query),
                    first.query_end,
                    first.query_end.saturating_add(query_ext).saturating_sub(1),
                    0,
                    i.saturating_sub(1),
                    Some(score_opts),
                    false,
                );
                if status == 0 {
                    mapper_set_splice_signal(first, second, signal);
                }
                return status;
            }
        }
    }

    mapper_clear_splice_signal(first, second);
    0
}

/// Port of NCBI `s_FindPolyAInSequence` (`hspfilter_mapper.c:1417`).
pub fn s_find_poly_a_in_sequence(sequence: Option<&[u8]>, length: i32) -> i32 {
    let Some(sequence) = sequence else {
        return -1;
    };
    let length = length.max(0).min(sequence.len() as i32);
    let k_base_a = 0;
    let k_max_errors = 3;
    let num_a;
    let mut err = 0;

    let mut i = length - 1;
    while i >= 0 && err < k_max_errors {
        if sequence[i as usize] != k_base_a {
            err += 1;
        }
        i -= 1;
    }
    i += 1;

    while i < length - 1
        && (sequence[i as usize] != k_base_a || sequence[(i + 1) as usize] != k_base_a)
    {
        if sequence[i as usize] != k_base_a {
            err -= 1;
        }
        i += 1;
    }

    num_a = length - i - err;
    if num_a < 3 || (num_a < 5 && err > 0) {
        return -1;
    }
    i
}

fn last_container(chain: &HSPChain) -> Option<&HSPContainer> {
    let mut h = chain.hsps.as_deref()?;
    while let Some(next) = h.next.as_deref() {
        h = next;
    }
    Some(h)
}

fn last_container_mut(chain: &mut HSPChain) -> Option<&mut HSPContainer> {
    let mut h = chain.hsps.as_deref_mut()?;
    loop {
        if h.next.is_none() {
            return Some(h);
        }
        h = h.next.as_deref_mut().unwrap();
    }
}

/// Port of NCBI `s_SetPolyATail` (`hspfilter_mapper.c:1462`).
pub fn s_set_poly_a_tail(
    chains: &mut Option<Box<HSPChain>>,
    positive_start: i32,
    negative_start: i32,
    query_len: i32,
) -> i32 {
    if chains.is_none() {
        return -1;
    }

    let mut ch = chains.as_deref_mut();
    while let Some(chain) = ch {
        let poly_a = if let Some(last) = last_container_mut(chain) {
            let hsp = &mut last.hsp;
            if query_len.saturating_sub(hsp.query_end) >= 5
                && ((hsp.query_frame < 0 && negative_start >= 0)
                    || (hsp.query_frame > 0 && positive_start >= 0))
            {
                mapper_mark_right_edge(hsp, MAPPER_POLY_A | MAPPER_EXON);
                Some(if hsp.query_frame > 0 {
                    positive_start.max(hsp.query_end)
                } else {
                    negative_start.max(hsp.query_end)
                })
            } else {
                None
            }
        } else {
            None
        };
        if let Some(poly_a) = poly_a {
            chain.poly_a = poly_a;
        }
        ch = chain.next.as_deref_mut();
    }

    0
}

/// Port of NCBI `s_FindPolyATails` (`hspfilter_mapper.c:1500`) for the Rust
/// chain/query-info model.
pub fn s_find_poly_a_tails(
    saved: &mut [Option<Box<HSPChain>>],
    query_sequence: &[u8],
    query_info: &QueryInfo,
) -> i32 {
    for query_idx in 0..query_info.num_queries.max(0) as usize {
        if query_idx >= saved.len() {
            break;
        }
        let Some(head) = saved[query_idx].as_ref() else {
            continue;
        };
        if head.adapter >= 0 {
            continue;
        }

        let plus_ctx = query_idx * NUM_STRANDS;
        let minus_ctx = plus_ctx + 1;
        let Some(plus) = query_info.contexts.get(plus_ctx) else {
            continue;
        };
        let query_len = plus.query_length.max(0);

        let mut best = head.as_ref();
        let mut cursor = head.next.as_deref();
        while let Some(ch) = cursor {
            if ch.score > best.score {
                best = ch;
            }
            cursor = ch.next.as_deref();
        }

        let Some(first) = best.hsps.as_deref() else {
            continue;
        };
        let Some(last) = last_container(best) else {
            continue;
        };
        let (from, to) = if first.hsp.query_frame > 0 {
            (first.hsp.query_offset, last.hsp.query_end)
        } else {
            (
                query_len.saturating_sub(last.hsp.query_end),
                query_len.saturating_sub(first.hsp.query_offset),
            )
        };
        if from < 4 && to > query_len.saturating_sub(3) {
            continue;
        }

        let plus_start = plus.query_offset.max(0) as usize;
        let plus_end = plus_start.saturating_add(query_len as usize);
        let positive_start =
            s_find_poly_a_in_sequence(query_sequence.get(plus_start..plus_end), query_len);

        let negative_start = if let Some(minus) = query_info.contexts.get(minus_ctx) {
            let minus_start = minus.query_offset.max(0) as usize;
            let minus_len = minus.query_length.max(0);
            let minus_end = minus_start.saturating_add(minus_len as usize);
            s_find_poly_a_in_sequence(query_sequence.get(minus_start..minus_end), minus_len)
        } else {
            -1
        };

        if positive_start >= 0 || negative_start >= 0 {
            s_set_poly_a_tail(
                &mut saved[query_idx],
                positive_start,
                negative_start,
                query_len,
            );
        }
    }

    0
}

/// Port of NCBI `s_FilterChains` (`hspfilter_mapper.c:2261`).
pub fn s_filter_chains(
    chains: &mut Option<Box<HSPChain>>,
    cutoff_score: i32,
    cutoff_edit_distance: i32,
) -> i32 {
    let mut kept = Vec::new();
    for chain in drain_chain_list(chains.take()) {
        if s_test_cutoffs(Some(&chain), cutoff_score, cutoff_edit_distance) {
            kept.push(chain);
        }
    }
    *chains = build_chain_list(kept);
    0
}

/// Port of NCBI `s_SortChains` (`hspfilter_mapper.c:2137`) for the Rust
/// linked-list representation. Chains are sorted in descending effective score.
pub fn s_sort_chains(saved: &mut [Option<Box<HSPChain>>]) -> i32 {
    for list in saved {
        if list.is_none() {
            continue;
        }
        let mut chains = drain_chain_list(list.take());
        chains.sort_by(|a, b| s_get_chain_score(Some(b)).cmp(&s_get_chain_score(Some(a))));
        *list = build_chain_list(chains);
    }
    0
}

/// Port of NCBI `s_CompareChainsByScore` (`hspfilter_mapper.c:1918`).
pub fn s_compare_chains_by_score(a: &HSPChain, b: &HSPChain) -> i32 {
    if a.score < b.score {
        1
    } else if a.score > b.score {
        -1
    } else {
        0
    }
}

/// Port of NCBI `s_CompareChainsByOid` (`hspfilter_mapper.c:2114`).
pub fn s_compare_chains_by_oid(a: &HSPChain, b: &HSPChain) -> i32 {
    if a.oid > b.oid {
        1
    } else if a.oid < b.oid {
        -1
    } else if s_find_fragment_start(Some(a)) > s_find_fragment_start(Some(b)) {
        1
    } else if s_find_fragment_start(Some(a)) < s_find_fragment_start(Some(b)) {
        -1
    } else {
        0
    }
}

/// Port of NCBI `s_HSPNodeArrayCopy` (`hspfilter_mapper.c:2431`).
pub fn s_hsp_node_array_copy(dest: &mut [HSPNode], source: &[HSPNode], num: i32) -> i32 {
    if num < 0 {
        return -1;
    }
    let num = num as usize;
    if dest.len() < num || source.len() < num {
        return -1;
    }

    for i in 0..num {
        dest[i] = source[i].clone();
        if let Some(next_idx) = source[i].path_next {
            if next_idx >= num {
                return -1;
            }
            dest[i].path_next = Some(next_idx);
        }
    }

    0
}

/// Rust ownership equivalent of NCBI `HSPPathFree` (`hspfilter_mapper.c:2454`).
pub fn hsp_path_free(_: Option<HSPPath>) -> Option<HSPPath> {
    None
}

/// Port of NCBI `HSPPathNew` (`hspfilter_mapper.c:2466`).
pub fn hsp_path_new() -> HSPPath {
    HSPPath {
        start: vec![None; MAX_NUM_HSP_PATHS],
        num_paths: 0,
        score: 0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hsp(score: i32, context: i32, q: i32, s: i32, ident: i32) -> Hsp {
        Hsp {
            score,
            num_ident: ident,
            bit_score: 0.0,
            evalue: 0.0,
            query_offset: q,
            query_end: q + 10,
            query_gapped_start: q,
            subject_offset: s,
            subject_end: s + 10,
            subject_gapped_start: s,
            context,
            query_frame: context,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        }
    }

    fn chain(score: i32, context: i32, oid: i32, q: i32, s: i32) -> Box<HSPChain> {
        let mut hsp_slot = Some(hsp(score, context, q, s, 9));
        let mut chain = hsp_chain_new(context);
        chain.oid = oid;
        chain.score = score;
        chain.hsps = hsp_container_new(&mut hsp_slot);
        chain
    }

    fn scores(list: &Option<Box<HSPChain>>) -> Vec<i32> {
        let mut out = Vec::new();
        let mut cursor = list.as_deref();
        while let Some(chain) = cursor {
            out.push(chain.score);
            cursor = chain.next.as_deref();
        }
        out
    }

    #[test]
    fn hsp_container_new_dup_free_take_ownership() {
        let mut missing_hsp = None;
        assert!(hsp_container_new(&mut missing_hsp).is_none());
        assert!(hsp_container_dup(None).is_none());

        let mut hsp_slot = Some(hsp(10, 0, 0, 0, 8));
        let container = hsp_container_new(&mut hsp_slot).unwrap();
        assert!(hsp_slot.is_none());
        let dup = hsp_container_dup(Some(&container)).unwrap();
        assert_eq!(dup.hsp.score, 10);
        assert!(hsp_container_free(Some(container)).is_none());
    }

    #[test]
    fn hsp_chain_new_clone_free_mapping_results() {
        let chain = chain(20, 2, 5, 0, 0);
        assert_eq!(chain.context, 2);
        assert_eq!(chain.adapter, -1);
        let cloned = clone_chain(Some(&chain)).unwrap();
        assert_eq!(cloned.context, 2);
        assert_eq!(cloned.oid, 5);
        assert!(cloned.next.is_none());
        assert!(clone_chain(None).is_none());
        assert!(hsp_chain_free(Some(chain)).is_none());
        assert!(hsp_chain_free(None).is_none());
        let results = blast_mapping_results_new();
        assert_eq!(results.num_queries, 0);
        assert!(blast_mapping_results_free(Some(results)).is_none());
        assert!(blast_mapping_results_free(None).is_none());
    }

    #[test]
    fn find_partialy_covered_queries_clones_matching_uncovered_chains() {
        let mut head = chain(35, 0, 7, 6, 0);
        head.next = Some(chain(40, 0, 8, 6, 0));
        let low_score = chain(20, 0, 7, 7, 0);
        let tail_uncovered = chain(32, 1, 7, 0, 0);
        let query_info = QueryInfo {
            num_queries: 2,
            contexts: vec![
                crate::queryinfo::ContextInfo {
                    query_offset: 0,
                    query_length: 30,
                    eff_searchsp: 0,
                    length_adjustment: 0,
                    query_index: 0,
                    frame: 0,
                    is_valid: true,
                    segment_flags: crate::queryinfo::E_NO_SEGMENTS,
                },
                crate::queryinfo::ContextInfo {
                    query_offset: 31,
                    query_length: 50,
                    eff_searchsp: 0,
                    length_adjustment: 0,
                    query_index: 1,
                    frame: 0,
                    is_valid: true,
                    segment_flags: crate::queryinfo::E_NO_SEGMENTS,
                },
            ],
            max_length: 50,
            min_length: 0,
        };
        let data = BlastHSPMapperData {
            query_info: Some(query_info),
            saved_chains: vec![Some(head), Some(tail_uncovered), Some(low_score)],
            ..Default::default()
        };

        let found = find_partialy_covered_queries(Some(&data), 7, 5);

        assert_eq!(scores(&found), vec![35, 32]);
        let first = found.as_ref().unwrap();
        assert_eq!(first.oid, 7);
        assert!(first.pair.is_none());
        assert_eq!(first.next.as_ref().unwrap().context, 1);
        assert!(find_partialy_covered_queries(None, 7, 5).is_none());
    }

    #[test]
    fn find_partialy_covered_queries_handles_extreme_trailing_endpoint() {
        let mut malformed = chain(35, 0, 7, 0, 0);
        malformed.hsps.as_mut().unwrap().hsp.query_end = i32::MIN;
        let query_info = QueryInfo {
            num_queries: 1,
            contexts: vec![crate::queryinfo::ContextInfo {
                query_offset: 0,
                query_length: 30,
                eff_searchsp: 0,
                length_adjustment: 0,
                query_index: 0,
                frame: 0,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: 30,
            min_length: 0,
        };
        let data = BlastHSPMapperData {
            query_info: Some(query_info),
            saved_chains: vec![Some(malformed)],
            ..Default::default()
        };

        let found = find_partialy_covered_queries(Some(&data), 7, 5);

        assert_eq!(scores(&found), vec![35]);
    }

    #[test]
    fn chain_list_insert_by_score_sorts_and_skips_duplicates() {
        let mut list = None;
        assert_eq!(
            s_hsp_chain_list_insert_one_by_score(&mut list, Some(chain(10, 0, 1, 0, 0)), true),
            0
        );
        assert_eq!(
            s_hsp_chain_list_insert_one_by_score(&mut list, Some(chain(30, 0, 2, 1, 1)), true),
            0
        );
        assert_eq!(scores(&list), vec![30, 10]);
        assert_eq!(
            s_hsp_chain_list_insert_one_by_score(&mut list, Some(chain(30, 0, 2, 1, 1)), true),
            0
        );
        assert_eq!(scores(&list), vec![30, 10]);
    }

    #[test]
    fn test_cutoffs_and_best_score_insert() {
        let mut list = Some(chain(20, 0, 1, 0, 0));
        assert_eq!(
            s_hsp_chain_list_insert_one(&mut list, Some(chain(10, 0, 2, 1, 1)), false),
            0
        );
        assert_eq!(scores(&list), vec![20]);
        assert_eq!(
            s_hsp_chain_list_insert_one(&mut list, Some(chain(40, 0, 3, 2, 2)), false),
            0
        );
        assert_eq!(scores(&list), vec![40]);

        assert!(s_test_cutoffs(list.as_deref(), 30, 1));
        assert!(!s_test_cutoffs(list.as_deref(), 30, 0));
        assert!(!s_test_cutoffs(None, 30, 0));

        let mut extreme = chain(100, 0, i32::MIN, 0, 0);
        extreme.hsps.as_mut().unwrap().hsp.query_end = i32::MAX;
        extreme.hsps.as_mut().unwrap().hsp.num_ident = 0;
        assert!(!s_test_cutoffs(Some(&extreme), 30, i32::MAX - 1));

        assert_eq!(
            s_hsp_chain_list_insert_one_by_score(&mut list, None, true),
            -1
        );
        assert_eq!(s_hsp_chain_list_insert_one(&mut list, None, true), -1);
        assert_eq!(scores(&list), vec![40]);
    }

    #[test]
    fn chain_list_insert_and_trim() {
        let mut incoming = Some(chain(50, 0, 1, 0, 0));
        incoming.as_mut().unwrap().next = Some(chain(45, 0, 2, 20, 20));
        incoming.as_mut().unwrap().next.as_mut().unwrap().next = Some(chain(20, 0, 3, 40, 40));

        let mut list = None;
        let mut missing_incoming = None;
        assert_eq!(
            hsp_chain_list_insert(&mut list, &mut missing_incoming, 0, -1, false),
            0
        );
        assert!(list.is_none());
        assert_eq!(
            hsp_chain_list_insert(&mut list, &mut incoming, 0, -1, false),
            0
        );
        assert!(incoming.is_none());
        assert_eq!(scores(&list), vec![50, 45, 20]);
        assert_eq!(hsp_chain_list_trim(&mut list, 5), 0);
        assert_eq!(scores(&list), vec![50, 45]);
        assert!(s_test_chains_sorted(list.as_deref()));

        let mut empty = None;
        assert_eq!(hsp_chain_list_trim(&mut empty, 5), -1);

        let mut below_cutoff = Some(chain(10, 0, 4, 60, 60));
        assert_eq!(
            hsp_chain_list_insert(&mut list, &mut below_cutoff, 30, -1, false),
            0
        );
        assert!(below_cutoff.is_none());
        assert_eq!(scores(&list), vec![50, 45]);

        let mut edit_distance_filtered = Some(chain(50, 0, 7, 80, 80));
        edit_distance_filtered
            .as_mut()
            .unwrap()
            .hsps
            .as_mut()
            .unwrap()
            .hsp
            .num_ident = 0;
        assert_eq!(
            hsp_chain_list_insert(&mut list, &mut edit_distance_filtered, 30, 1, false),
            0
        );
        assert!(edit_distance_filtered.is_none());
        assert_eq!(scores(&list), vec![50, 45]);

        let mut bad_context = chain(30, 0, 5, 0, 0);
        bad_context.hsps.as_mut().unwrap().hsp.context = 1;
        assert!(!s_test_chains(Some(&bad_context)));

        let mut bad_order = chain(30, 0, 6, 10, 10);
        let mut previous_hsp = Some(hsp(20, 0, 5, 5, 10));
        bad_order.hsps.as_mut().unwrap().next = hsp_container_new(&mut previous_hsp);
        assert!(!s_test_chains(Some(&bad_order)));

        let mut unsorted = chain(20, 0, 7, 0, 0);
        unsorted.next = Some(chain(30, 0, 8, 20, 20));
        assert!(!s_test_chains_sorted(Some(&unsorted)));
    }

    #[test]
    fn mapper_paired_init_allocates_saved_chain_slots() {
        let mut data = BlastHSPMapperData::default();
        assert_eq!(s_blast_hsp_mapper_paired_init(&mut data, 3), 0);
        assert_eq!(data.saved_chains.len(), 3);
        assert!(data.saved_chains.iter().all(Option::is_none));
        data.saved_chains[1] = Some(chain(20, 0, 1, 0, 0));
        assert_eq!(s_blast_hsp_mapper_paired_init(&mut data, 2), 0);
        assert_eq!(data.saved_chains.len(), 2);
        assert!(data.saved_chains.iter().all(Option::is_none));
        assert_eq!(s_blast_hsp_mapper_paired_init(&mut data, -1), 0);
        assert!(data.saved_chains.is_empty());
    }

    #[test]
    fn mapper_params_and_writer_lifecycle_match_c_setup() {
        let mut hit_options = HitSavingOptions::default();
        hit_options.hitlist_size = 4;
        hit_options.cutoff_score = 37;
        let scoring_options = ScoringOptions::new_blastn();

        let params = blast_hsp_mapper_params_new(Some(&hit_options), Some(&scoring_options))
            .expect("mapper params");
        assert_eq!(params.program, MAPPING);
        assert_eq!(params.hitlist_size, 10);
        assert_eq!(params.cutoff_score, 37);
        assert_eq!(params.cutoff_edit_dist, -1);
        assert_eq!(params.scoring_options.gap_open, -scoring_options.gap_open);
        assert_eq!(
            params.scoring_options.gap_extend,
            -scoring_options.gap_extend
        );
        assert!(blast_hsp_mapper_params_new(None, Some(&scoring_options)).is_none());
        assert!(blast_hsp_mapper_params_new(Some(&hit_options), None).is_none());

        let info = blast_hsp_mapper_info_new(params.clone());
        assert_eq!(info.new_fn, BlastHSPMapperCallback::PairedNew);

        let query_info = QueryInfo::new_blastn(&[12]);
        let query = [0, 1, 2, 3];
        let writer = s_blast_hsp_mapper_paired_new(Some(params), Some(&query_info), Some(&query));
        assert_eq!(writer.init_fn, BlastHSPMapperCallback::PairedInit);
        assert_eq!(writer.final_fn, BlastHSPMapperCallback::Final);
        assert_eq!(writer.free_fn, BlastHSPMapperCallback::Free);
        assert_eq!(writer.run_fn, BlastHSPMapperCallback::SplicedPairedRun);
        assert_eq!(writer.data.query.as_deref(), Some(query.as_slice()));
        assert_eq!(
            writer.data.query_info.as_ref().map(|info| info.num_queries),
            Some(1)
        );
        assert!(s_blast_hsp_mapper_free(Some(writer)).is_none());

        let missing_writer = s_blast_hsp_mapper_paired_new(None, None, None);
        assert!(missing_writer.data.params.is_none());
        assert!(missing_writer.data.query.is_none());
        assert!(missing_writer.data.query_info.is_none());
        assert!(missing_writer.data.saved_chains.is_empty());

        assert!(blast_hsp_mapper_params_free(None).is_none());
    }

    #[test]
    fn overlap_cost_handles_containment_disjoint_and_edit_penalties() {
        let a = hsp(50, 0, 0, 0, 10);
        let mut b = hsp(20, 0, 3, 3, 7);
        b.query_end = 8;
        assert_eq!(s_get_overlap_cost(&a, &b, 4), 20);

        let c = hsp(30, 0, 20, 20, 9);
        assert_eq!(s_get_overlap_cost(&a, &c, 4), 0);

        let mut d = hsp(30, 0, 6, 6, 9);
        d.query_end = 16;
        assert_eq!(s_get_overlap_cost_with_edits(&a, &d, 2, &[11], &[8]), 2);
        assert_eq!(s_get_overlap_cost_with_edits(&d, &a, 2, &[8], &[11]), 2);

        let mut edited_a = hsp(50, 0, 0, 0, 10);
        edited_a.map_info = Some(crate::hspstream::BlastHSPMappingInfo {
            edits: Some(crate::gapinfo::JumperEditsBlock {
                num_edits: 0,
                edits: vec![JumperEdit {
                    query_pos: 8,
                    query_base: 0,
                    subject_base: 1,
                }],
            }),
            subject_overhangs: None,
            left_edge: 0,
            right_edge: 0,

            flags: 0,
        });
        let mut edited_b = hsp(30, 0, 6, 6, 9);
        edited_b.query_end = 16;
        assert_eq!(s_get_overlap_cost(&edited_a, &edited_b, 2), 2);

        let mut reverse_containment = hsp(70, 0, -2, -2, 10);
        reverse_containment.query_end = 20;
        assert_eq!(s_get_overlap_cost(&a, &reverse_containment, 4), 50);

        let mut extreme_first = hsp(10, 0, i32::MIN, 0, 1);
        extreme_first.query_end = 0;
        let mut extreme_second = hsp(11, 0, -1, 0, 1);
        extreme_second.query_end = i32::MAX;
        assert_eq!(
            s_get_overlap_cost_with_edits(&extreme_first, &extreme_second, i32::MAX, &[0], &[-1],),
            i32::MIN + 2
        );
    }

    #[test]
    fn overlap_cost_wrapper_uses_mapper_edits_from_both_hsps() {
        let mut first = hsp(40, 0, 0, 0, 10);
        first.map_info = Some(crate::hspstream::BlastHSPMappingInfo {
            edits: Some(crate::gapinfo::JumperEditsBlock {
                num_edits: 0,
                edits: vec![
                    JumperEdit {
                        query_pos: 5,
                        query_base: 0,
                        subject_base: 1,
                    },
                    JumperEdit {
                        query_pos: 6,
                        query_base: 0,
                        subject_base: 1,
                    },
                ],
            }),
            subject_overhangs: None,
            left_edge: 0,
            right_edge: 0,

            flags: 0,
        });

        let mut second = hsp(35, 0, 5, 5, 10);
        second.query_end = 15;
        second.map_info = Some(crate::hspstream::BlastHSPMappingInfo {
            edits: Some(crate::gapinfo::JumperEditsBlock {
                num_edits: 0,
                edits: vec![
                    JumperEdit {
                        query_pos: 4,
                        query_base: 0,
                        subject_base: 1,
                    },
                    JumperEdit {
                        query_pos: 9,
                        query_base: 0,
                        subject_base: 1,
                    },
                ],
            }),
            subject_overhangs: None,
            left_edge: 0,
            right_edge: 0,

            flags: 0,
        });

        assert_eq!(
            s_get_overlap_cost(&first, &second, 2),
            s_get_overlap_cost_with_edits(&first, &second, 2, &[5, 6], &[4, 9])
        );
        assert_eq!(s_get_overlap_cost(&first, &second, 2), 1);
    }

    #[test]
    fn compute_chain_score_adds_hsps_and_gap_penalties() {
        let mut scored_chain = chain(30, 0, 1, 0, 0);
        let mut second_hsp = Some(hsp(20, 0, 15, 14, 8));
        scored_chain.hsps.as_mut().unwrap().next = hsp_container_new(&mut second_hsp);
        let opts = ScoringOptions::new_blastn();
        assert_eq!(
            s_compute_chain_score(Some(&mut scored_chain), Some(&opts), 100, false),
            18
        );
        assert_eq!(scored_chain.score, 18);

        assert_eq!(
            s_compute_chain_score(Some(&mut scored_chain), None, 100, false),
            -1
        );
        assert_eq!(scored_chain.score, 18);

        let mut empty_chain = hsp_chain_new(0);
        empty_chain.score = 27;
        assert_eq!(
            s_compute_chain_score(Some(&mut empty_chain), Some(&opts), 100, false),
            0
        );
        assert_eq!(empty_chain.score, 0);

        let mut scripted = chain(99, 0, 2, 0, 0);
        scripted.hsps.as_mut().unwrap().hsp.edit_script =
            Some(crate::gapinfo::GapEditScript::from_ops(vec![
                (GapAlignOpType::Sub, 5),
                (GapAlignOpType::Ins, 2),
            ]));
        assert_eq!(
            s_compute_chain_score(Some(&mut scripted), Some(&opts), 100, true),
            -4
        );
        assert_eq!(scripted.score, -4);

        let mut edited = chain(99, 0, 2, 0, 0);
        edited.hsps.as_mut().unwrap().hsp.edit_script = Some(
            crate::gapinfo::GapEditScript::from_ops(vec![(GapAlignOpType::Sub, 5)]),
        );
        edited.hsps.as_mut().unwrap().hsp.map_info = Some(crate::hspstream::BlastHSPMappingInfo {
            edits: Some(crate::gapinfo::JumperEditsBlock {
                num_edits: 0,
                edits: vec![JumperEdit {
                    query_pos: 4,
                    query_base: b'A',
                    subject_base: b'C',
                }],
            }),
            subject_overhangs: None,
            left_edge: 0,
            right_edge: 0,

            flags: 0,
        });
        assert_eq!(
            s_compute_chain_score(Some(&mut edited), Some(&opts), 100, true),
            1
        );
        assert_eq!(edited.score, 1);
        assert_eq!(edited.hsps.as_ref().unwrap().hsp.score, 1);

        let mut extreme_chain = chain(30, 0, 1, i32::MIN, i32::MIN);
        let mut extreme_second_hsp = hsp(20, 0, 0, 0, 1);
        extreme_second_hsp.query_offset = i32::MAX;
        extreme_second_hsp.subject_offset = i32::MAX;
        let mut extreme_second = Some(extreme_second_hsp);
        extreme_chain.hsps.as_mut().unwrap().next = hsp_container_new(&mut extreme_second);
        assert_eq!(
            s_compute_chain_score(Some(&mut extreme_chain), Some(&opts), 100, false),
            18
        );

        assert_eq!(s_compute_chain_score(None, Some(&opts), 100, false), -1);
    }

    #[test]
    fn compute_chain_score_skips_gap_penalty_between_splice_edges() {
        let mut scored_chain = chain(30, 0, 1, 0, 0);
        scored_chain.hsps.as_mut().unwrap().hsp.map_info =
            Some(crate::hspstream::BlastHSPMappingInfo {
                edits: None,
                subject_overhangs: None,
                left_edge: 0,
                right_edge: MAPPER_SPLICE_SIGNAL,

                flags: 0,
            });
        let mut second = hsp(20, 0, 15, 14, 8);
        second.map_info = Some(crate::hspstream::BlastHSPMappingInfo {
            edits: None,
            subject_overhangs: None,
            left_edge: MAPPER_SPLICE_SIGNAL,
            right_edge: 0,

            flags: 0,
        });
        let mut second_hsp = Some(second);
        scored_chain.hsps.as_mut().unwrap().next = hsp_container_new(&mut second_hsp);
        let opts = ScoringOptions::new_blastn();

        assert_eq!(
            s_compute_chain_score(Some(&mut scored_chain), Some(&opts), 100, false),
            50
        );
    }

    #[test]
    fn mapper_cutoff_and_chain_from_hsp_follow_context_rules() {
        let mut hit_options = HitSavingOptions::default();
        hit_options.cutoff_score = 0;
        let scores = ScoringOptions::new_blastn();
        let mut params = blast_hsp_mapper_params_new(Some(&hit_options), Some(&scores)).unwrap();
        params.cutoff_score_fun = [1000, 200];
        assert_eq!(s_mapper_cutoff_score(&params, 50), 110);

        params.cutoff_score_fun = [0, 0];
        params.cutoff_score = 37;
        assert_eq!(s_mapper_cutoff_score(&params, 50), 37);

        let query_info = QueryInfo::new_blastn(&[20]);
        let mut minus_frame_hsp = hsp(30, 1, 0, 0, 10);
        minus_frame_hsp.query_frame = 0;
        let chain = chain_from_hsp(minus_frame_hsp, 7, &query_info).unwrap();
        assert_eq!(chain.oid, 7);
        assert_eq!(chain.context, 1);
        assert_eq!(chain.count, 1);
        assert_eq!(chain.hsps.as_ref().unwrap().hsp.query_frame, -1);

        let mut preserved_frame = hsp(30, 0, 0, 0, 10);
        preserved_frame.query_frame = 2;
        let chain = chain_from_hsp(preserved_frame, 8, &query_info).unwrap();
        assert_eq!(chain.hsps.as_ref().unwrap().hsp.query_frame, 2);

        assert!(chain_from_hsp(hsp(30, -1, 0, 0, 10), 7, &query_info).is_none());
    }

    #[test]
    fn fragment_start_and_chain_comparators_match_c_ordering() {
        let mut plus = chain(30, 1, 4, 0, 11);
        let mut minus = chain(40, -1, 4, 0, 20);
        let mut tail = Some(hsp(5, -1, 20, 50, 5));
        minus.hsps.as_mut().unwrap().next = hsp_container_new(&mut tail);

        assert_eq!(s_find_fragment_start(Some(&plus)), 11);
        assert_eq!(s_find_fragment_start(Some(&minus)), 59);
        assert_eq!(s_find_fragment_start(None), -1);
        assert_eq!(s_find_fragment_start(Some(&hsp_chain_new(0))), -1);
        let mut extreme_minus = minus.clone();
        extreme_minus
            .hsps
            .as_mut()
            .unwrap()
            .next
            .as_mut()
            .unwrap()
            .hsp
            .subject_end = i32::MIN;
        assert_eq!(s_find_fragment_start(Some(&extreme_minus)), i32::MIN);

        assert_eq!(s_compare_chains_by_score(&plus, &minus), 1);
        plus.score = 40;
        assert_eq!(s_compare_chains_by_score(&plus, &minus), 0);
        plus.oid = 3;
        minus.oid = 4;
        assert_eq!(s_compare_chains_by_oid(&plus, &minus), -1);
        plus.oid = 4;
        assert_eq!(s_compare_chains_by_oid(&plus, &minus), -1);
        minus.hsps.as_mut().unwrap().next = None;
        minus.hsps.as_mut().unwrap().hsp.subject_end = 12;
        assert_eq!(s_compare_chains_by_oid(&plus, &minus), 0);
    }

    #[test]
    fn test_hsp_ranges_checks_spans_and_edit_script() {
        let mut valid = hsp(10, 0, 4, 8, 10);
        valid.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![(
            GapAlignOpType::Sub,
            10,
        )]));
        assert!(s_test_hsp_ranges(&valid));

        let mut invalid = valid.clone();
        invalid.query_end = invalid.query_offset + 9;
        assert!(!s_test_hsp_ranges(&invalid));

        let mut negative = valid.clone();
        negative.query_offset = -1;
        assert!(!s_test_hsp_ranges(&negative));

        let mut reversed_subject = valid.clone();
        reversed_subject.subject_end = reversed_subject.subject_offset;
        assert!(!s_test_hsp_ranges(&reversed_subject));

        let mut oversized_script = hsp(0, 0, 0, 0, 1);
        oversized_script.query_end = i32::MAX;
        oversized_script.subject_end = i32::MAX;
        oversized_script.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![
            (GapAlignOpType::Sub, i32::MAX),
            (GapAlignOpType::Ins, 1),
        ]));
        assert!(!s_test_hsp_ranges(&oversized_script));
    }

    #[test]
    fn trim_hsp_rejects_oversized_trims_and_preserves_noop_state() {
        let mut hsp = hsp(10, 0, 4, 8, 10);
        hsp.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![(
            GapAlignOpType::Sub,
            10,
        )]));
        hsp.map_info = Some(crate::hspstream::BlastHSPMappingInfo {
            edits: Some(crate::gapinfo::JumperEditsBlock {
                num_edits: 0,
                edits: vec![JumperEdit {
                    query_pos: 6,
                    query_base: 1,
                    subject_base: 2,
                }],
            }),
            subject_overhangs: Some(crate::gapinfo::SequenceOverhangs {
                left_len: 0,
                right_len: 0,
                left: Some(vec![1]),
                right: Some(vec![2]),
            }),
            left_edge: 3,
            right_edge: 4,

            flags: 0,
        });
        let original = hsp.clone();
        let assert_unchanged = |hsp: &Hsp| {
            assert_eq!(hsp.score, original.score);
            assert_eq!(hsp.num_ident, original.num_ident);
            assert_eq!(hsp.query_offset, original.query_offset);
            assert_eq!(hsp.query_end, original.query_end);
            assert_eq!(hsp.subject_offset, original.subject_offset);
            assert_eq!(hsp.subject_end, original.subject_end);
            assert_eq!(
                hsp.edit_script.as_ref().map(|script| script.ops_vec()),
                original.edit_script.as_ref().map(|script| script.ops_vec())
            );
            assert_eq!(
                hsp.map_info.as_ref().and_then(|info| info.edits.clone()),
                original
                    .map_info
                    .as_ref()
                    .and_then(|info| info.edits.clone())
            );
            assert_eq!(
                hsp.map_info
                    .as_ref()
                    .and_then(|info| info.subject_overhangs.clone()),
                original
                    .map_info
                    .as_ref()
                    .and_then(|info| info.subject_overhangs.clone())
            );
            assert_eq!(
                hsp.map_info.as_ref().map(|info| info.left_edge),
                original.map_info.as_ref().map(|info| info.left_edge)
            );
            assert_eq!(
                hsp.map_info.as_ref().map(|info| info.right_edge),
                original.map_info.as_ref().map(|info| info.right_edge)
            );
        };

        assert_eq!(s_trim_hsp(&mut hsp, 0, true, true, -4, -12, -1, None), 0);
        assert_unchanged(&hsp);

        assert_eq!(s_trim_hsp(&mut hsp, 11, true, true, -4, -12, -1, None), -1);
        assert_unchanged(&hsp);

        assert_eq!(
            s_trim_hsp(&mut hsp, 11, false, false, -4, -12, -1, None),
            -1
        );
        assert_unchanged(&hsp);

        let mut extreme_score = original.clone();
        extreme_score.query_offset = 0;
        extreme_score.subject_offset = 0;
        extreme_score.query_end = i32::MAX;
        extreme_score.subject_end = i32::MAX;
        extreme_score.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![
            (GapAlignOpType::Sub, i32::MAX),
            (GapAlignOpType::Sub, i32::MAX),
        ]));
        assert_eq!(
            s_trim_hsp(&mut extreme_score, 1, true, true, -4, -12, -1, None),
            0
        );
        assert_eq!(extreme_score.score, i32::MAX - 6);
        assert_eq!(extreme_score.num_ident, i32::MAX - 2);
    }

    #[test]
    fn trim_hsp_start_updates_script_and_coordinates() {
        let mut hsp = hsp(0, 0, 4, 8, 10);
        hsp.query_end = 18;
        hsp.subject_end = 22;
        hsp.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![
            (GapAlignOpType::Sub, 4),
            (GapAlignOpType::Ins, 2),
            (GapAlignOpType::Sub, 8),
        ]));
        hsp.map_info = Some(crate::hspstream::BlastHSPMappingInfo {
            edits: Some(crate::gapinfo::JumperEditsBlock {
                num_edits: 0,
                edits: vec![
                    JumperEdit {
                        query_pos: 2,
                        query_base: 0,
                        subject_base: 9,
                    },
                    JumperEdit {
                        query_pos: 6,
                        query_base: 1,
                        subject_base: 2,
                    },
                    JumperEdit {
                        query_pos: 10,
                        query_base: 3,
                        subject_base: 1,
                    },
                ],
            }),
            subject_overhangs: Some(crate::gapinfo::SequenceOverhangs {
                left_len: 0,
                right_len: 0,
                left: Some(vec![2]),
                right: None,
            }),
            left_edge: 0,
            right_edge: 0,

            flags: 0,
        });
        let query = vec![0, 1, 2, 3, 1, 1, 2, 3, 0, 1, 3, 3, 0, 1, 2, 3, 0, 1];

        assert_eq!(
            s_trim_hsp(&mut hsp, 5, true, true, -4, -12, -1, Some(&query)),
            0
        );
        assert_eq!(hsp.query_offset, 9);
        assert_eq!(hsp.subject_offset, 12);
        assert_eq!(
            hsp.edit_script.as_ref().unwrap().ops_vec(),
            vec![(GapAlignOpType::Ins, 1), (GapAlignOpType::Sub, 8)]
        );
        assert_eq!(hsp.num_ident, 8);
        let map_info = hsp.map_info.as_ref().unwrap();
        assert_eq!(
            map_info.subject_overhangs.as_ref().unwrap().left,
            Some(vec![2, 1, 1, 9, 3])
        );
        assert_eq!(
            map_info.edits.as_ref().unwrap().edits,
            vec![JumperEdit {
                query_pos: 10,
                query_base: 3,
                subject_base: 1,
            }]
        );
    }

    #[test]
    fn trim_hsp_end_can_trim_subject_through_gap_ops() {
        let mut hsp = hsp(0, 0, 4, 8, 10);
        hsp.query_end = 18;
        hsp.subject_end = 22;
        hsp.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![
            (GapAlignOpType::Sub, 8),
            (GapAlignOpType::Del, 3),
            (GapAlignOpType::Sub, 3),
        ]));
        hsp.map_info = Some(crate::hspstream::BlastHSPMappingInfo {
            edits: Some(crate::gapinfo::JumperEditsBlock {
                num_edits: 0,
                edits: vec![
                    JumperEdit {
                        query_pos: 6,
                        query_base: 1,
                        subject_base: 2,
                    },
                    JumperEdit {
                        query_pos: 12,
                        query_base: 15,
                        subject_base: 0,
                    },
                    JumperEdit {
                        query_pos: 12,
                        query_base: 15,
                        subject_base: 1,
                    },
                    JumperEdit {
                        query_pos: 16,
                        query_base: 2,
                        subject_base: 7,
                    },
                ],
            }),
            subject_overhangs: Some(crate::gapinfo::SequenceOverhangs {
                left_len: 0,
                right_len: 0,
                left: None,
                right: Some(vec![3]),
            }),
            left_edge: 0,
            right_edge: 0,

            flags: 0,
        });
        let query = vec![0, 1, 2, 3, 1, 2, 1, 0, 2, 3, 0, 1, 2, 0, 3, 3, 2, 1];

        assert_eq!(
            s_trim_hsp(&mut hsp, 4, false, false, -4, -12, -1, Some(&query)),
            0
        );
        assert_eq!(hsp.query_end, 15);
        assert_eq!(hsp.subject_end, 18);
        assert_eq!(
            hsp.edit_script.as_ref().unwrap().ops_vec(),
            vec![(GapAlignOpType::Sub, 8), (GapAlignOpType::Del, 2)]
        );
        assert_eq!(hsp.num_ident, 10);
        let map_info = hsp.map_info.as_ref().unwrap();
        assert_eq!(
            map_info.subject_overhangs.as_ref().unwrap().right,
            Some(vec![3, 7, 2, 1, 3])
        );
        assert_eq!(
            map_info.edits.as_ref().unwrap().edits,
            vec![
                JumperEdit {
                    query_pos: 6,
                    query_base: 1,
                    subject_base: 2,
                },
                JumperEdit {
                    query_pos: 12,
                    query_base: 15,
                    subject_base: 0,
                },
                JumperEdit {
                    query_pos: 12,
                    query_base: 15,
                    subject_base: 1,
                },
            ]
        );
    }

    #[test]
    fn trim_overlap_removes_query_overlap_from_second_when_it_survives() {
        let mut first = hsp(0, 0, 0, 20, 10);
        first.query_end = 10;
        first.subject_end = 30;
        first.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![(
            GapAlignOpType::Sub,
            10,
        )]));

        let mut second = hsp(0, 0, 7, 40, 10);
        second.query_end = 20;
        second.subject_end = 53;
        second.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![(
            GapAlignOpType::Sub,
            13,
        )]));

        assert_eq!(s_trim_overlap(&mut first, &mut second, None), 0);
        assert_eq!(first.query_end, 10);
        assert_eq!(second.query_offset, 10);
        assert_eq!(second.subject_offset, 43);
        assert_eq!(
            second.edit_script.as_ref().unwrap().ops_vec(),
            vec![(GapAlignOpType::Sub, 10)]
        );
    }

    #[test]
    fn trim_overlap_removes_subject_overlap_from_first_when_second_would_vanish() {
        let mut first = hsp(0, 0, 0, 0, 10);
        first.query_end = 10;
        first.subject_end = 20;
        first.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![(
            GapAlignOpType::Sub,
            20,
        )]));

        let mut second = hsp(0, 0, 12, 15, 5);
        second.query_end = 17;
        second.subject_end = 20;
        second.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![(
            GapAlignOpType::Sub,
            5,
        )]));

        assert_eq!(s_trim_overlap(&mut first, &mut second, None), 0);
        assert_eq!(first.subject_end, 15);
        assert_eq!(second.subject_offset, 15);
        assert_eq!(second.subject_end, 20);
        assert_eq!(
            first.edit_script.as_ref().unwrap().ops_vec(),
            vec![(GapAlignOpType::Sub, 15)]
        );
    }

    #[test]
    fn trim_overlap_rejects_extreme_malformed_overlap_without_overflow() {
        let mut first = hsp(0, 0, 0, 0, 10);
        first.query_end = i32::MAX;
        first.subject_end = i32::MAX;
        let mut second = hsp(0, 0, 0, 0, 10);
        second.query_offset = i32::MIN;
        second.subject_offset = i32::MIN;
        second.query_end = 0;
        second.subject_end = 0;

        assert_eq!(s_trim_overlap(&mut first, &mut second, None), -1);
    }

    #[test]
    fn trim_chain_start_to_subject_pos_drops_and_trims_boundary_hsp() {
        let mut first = hsp(8, 0, 0, 0, 8);
        first.query_end = 8;
        first.subject_end = 8;
        first.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![(
            GapAlignOpType::Sub,
            8,
        )]));
        let mut second = hsp(12, 0, 10, 10, 10);
        second.query_end = 20;
        second.subject_end = 20;
        second.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![(
            GapAlignOpType::Sub,
            10,
        )]));
        second.map_info = Some(crate::hspstream::BlastHSPMappingInfo {
            edits: None,
            subject_overhangs: None,
            left_edge: MAPPER_SPLICE_SIGNAL | MAPPER_EXON | 3,
            right_edge: MAPPER_SPLICE_SIGNAL | MAPPER_EXON | 2,

            flags: 0,
        });
        let mut first_slot = Some(first);
        let mut second_slot = Some(second);
        let mut chain = hsp_chain_new(0);
        chain.score = 20;
        chain.hsps = hsp_container_new(&mut first_slot);
        chain.hsps.as_mut().unwrap().next = hsp_container_new(&mut second_slot);

        assert_eq!(
            s_trim_chain_start_to_subj_pos(Some(&mut chain), 13, -4, -12, -1, Some(&[0; 40])),
            0
        );

        let head = chain.hsps.as_ref().unwrap();
        assert_eq!(head.hsp.query_offset, 13);
        assert_eq!(head.hsp.subject_offset, 13);
        assert_eq!(head.hsp.query_end, 20);
        assert_eq!(head.hsp.map_info.as_ref().unwrap().left_edge, 3);
        assert_ne!(
            head.hsp.map_info.as_ref().unwrap().right_edge & MAPPER_SPLICE_SIGNAL,
            0
        );
        assert_eq!(chain.score, 7);
        assert!(head.next.is_none());
    }

    #[test]
    fn trim_chain_end_to_subject_pos_drops_tail_and_trims_boundary_hsp() {
        let mut first = hsp(12, 0, 0, 0, 10);
        first.query_end = 10;
        first.subject_end = 10;
        first.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![(
            GapAlignOpType::Sub,
            10,
        )]));
        let mut second = hsp(14, 0, 12, 12, 12);
        second.query_end = 24;
        second.subject_end = 24;
        second.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![(
            GapAlignOpType::Sub,
            12,
        )]));
        second.map_info = Some(crate::hspstream::BlastHSPMappingInfo {
            edits: None,
            subject_overhangs: None,
            left_edge: MAPPER_SPLICE_SIGNAL | MAPPER_EXON | 1,
            right_edge: MAPPER_SPLICE_SIGNAL | MAPPER_EXON | 2,

            flags: 0,
        });
        let mut third = hsp(9, 0, 30, 30, 9);
        third.query_end = 39;
        third.subject_end = 39;
        third.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![(
            GapAlignOpType::Sub,
            9,
        )]));
        let mut first_slot = Some(first);
        let mut second_slot = Some(second);
        let mut third_slot = Some(third);
        let mut chain = hsp_chain_new(0);
        chain.score = 35;
        chain.hsps = hsp_container_new(&mut first_slot);
        chain.hsps.as_mut().unwrap().next = hsp_container_new(&mut second_slot);
        chain.hsps.as_mut().unwrap().next.as_mut().unwrap().next =
            hsp_container_new(&mut third_slot);

        assert_eq!(
            s_trim_chain_end_to_subj_pos(Some(&mut chain), 18, -4, -12, -1, Some(&[0; 40])),
            0
        );

        let head = chain.hsps.as_ref().unwrap();
        let second = head.next.as_ref().unwrap();
        assert_eq!(second.hsp.query_end, 18);
        assert_eq!(second.hsp.subject_end, 18);
        assert_ne!(
            second.hsp.map_info.as_ref().unwrap().left_edge & MAPPER_SPLICE_SIGNAL,
            0
        );
        assert_eq!(second.hsp.map_info.as_ref().unwrap().right_edge, 2);
        assert_eq!(
            second.hsp.edit_script.as_ref().unwrap().ops_vec(),
            vec![(GapAlignOpType::Sub, 6)]
        );
        assert_eq!(chain.score, 18);
        assert!(second.next.is_none());
    }

    #[test]
    fn set_adapter_trims_plus_chain_after_adapter_start() {
        let mut list = Some(chain(30, 1, 1, 0, 0));
        let mut second = hsp(20, 1, 12, 12, 13);
        second.query_end = 25;
        second.subject_end = 25;
        let mut second_slot = Some(second);
        list.as_mut().unwrap().hsps.as_mut().unwrap().next = hsp_container_new(&mut second_slot);

        let scores = ScoringOptions::new_blastn();
        assert_eq!(s_set_adapter(&mut list, 20, None, 40, Some(&scores)), 0);
        let chain = list.as_ref().unwrap();
        assert_eq!(chain.adapter, 20);
        let first = chain.hsps.as_ref().unwrap();
        let second = first.next.as_ref().unwrap();
        assert_eq!(second.hsp.query_end, 20);
        let edge = second.hsp.map_info.as_ref().unwrap().right_edge;
        assert_eq!(edge & MAPPER_ADAPTER, MAPPER_ADAPTER);
        assert_eq!(edge & MAPPER_EXON, MAPPER_EXON);
        assert!(second.next.is_none());
    }

    #[test]
    fn set_adapter_drops_minus_hsps_before_adapter_and_trims_start() {
        let mut list = Some(chain(30, -1, 1, 4, 4));
        list.as_mut().unwrap().hsps.as_mut().unwrap().hsp.query_end = 8;
        list.as_mut()
            .unwrap()
            .hsps
            .as_mut()
            .unwrap()
            .hsp
            .subject_end = 8;
        let mut second = hsp(20, -1, 10, 10, 15);
        second.query_end = 25;
        second.subject_end = 25;
        let mut second_slot = Some(second);
        list.as_mut().unwrap().hsps.as_mut().unwrap().next = hsp_container_new(&mut second_slot);

        let scores = ScoringOptions::new_blastn();
        assert_eq!(s_set_adapter(&mut list, 20, None, 40, Some(&scores)), 0);
        let chain = list.as_ref().unwrap();
        assert_eq!(chain.adapter, 20);
        let first = chain.hsps.as_ref().unwrap();
        assert_eq!(first.hsp.query_offset, 20);
        assert_eq!(first.hsp.query_end, 25);
        let edge = first.hsp.map_info.as_ref().unwrap().left_edge;
        assert_eq!(edge & MAPPER_ADAPTER, MAPPER_ADAPTER);
        assert_eq!(edge & MAPPER_EXON, MAPPER_EXON);
        assert!(first.next.is_none());
    }

    #[test]
    fn set_adapter_and_poly_a_tail_report_invalid_inputs_without_mutating_chain() {
        let mut missing_chain = None;
        let scores = ScoringOptions::new_blastn();
        assert_eq!(
            s_set_adapter(&mut missing_chain, 20, None, 40, Some(&scores)),
            -1
        );
        assert!(missing_chain.is_none());
        assert_eq!(s_set_poly_a_tail(&mut missing_chain, 12, -1, 20), -1);
        assert!(missing_chain.is_none());

        let mut list = Some(chain(30, 1, 1, 0, 0));
        let original_score = list.as_ref().unwrap().score;
        let original_adapter = list.as_ref().unwrap().adapter;
        let original_query_end = list.as_ref().unwrap().hsps.as_ref().unwrap().hsp.query_end;

        assert_eq!(s_set_adapter(&mut list, -1, None, 40, Some(&scores)), -1);
        assert_eq!(s_set_adapter(&mut list, 20, None, 40, None), -1);
        let chain = list.as_ref().unwrap();
        assert_eq!(chain.score, original_score);
        assert_eq!(chain.adapter, original_adapter);
        let hsp = &chain.hsps.as_ref().unwrap().hsp;
        assert_eq!(hsp.query_end, original_query_end);
        assert!(hsp.map_info.is_none());

        assert_eq!(
            s_set_adapter(&mut list, 20, None, i32::MIN, Some(&scores)),
            0
        );
        let chain = list.as_ref().unwrap();
        assert_eq!(chain.score, original_score);
        assert_eq!(chain.adapter, original_adapter);
        assert_eq!(
            chain.hsps.as_ref().unwrap().hsp.query_end,
            original_query_end
        );
        assert!(chain.hsps.as_ref().unwrap().hsp.map_info.is_none());
    }

    #[test]
    fn find_adapters_searches_longest_overhang_and_sets_adapter() {
        let mut query = vec![3; 122];
        query[30..42].copy_from_slice(&[0, 2, 0, 3, 1, 2, 2, 0, 0, 2, 0, 2]);
        let query_info = QueryInfo::new_blastn(&[60]);

        let mut list = Some(chain(30, 1, 1, 24, 24));
        {
            let hsp = &mut list.as_mut().unwrap().hsps.as_mut().unwrap().hsp;
            hsp.query_end = 42;
            hsp.subject_end = 42;
        }
        let mut saved = vec![list];
        let scores = ScoringOptions::new_blastn();

        assert_eq!(s_find_adapters(&mut saved, &query, &query_info, &scores), 0);
        let chain = saved[0].as_ref().unwrap();
        assert_eq!(chain.adapter, 30);
        let hsp = &chain.hsps.as_ref().unwrap().hsp;
        assert_eq!(hsp.query_offset, 24);
        assert_eq!(hsp.query_end, 30);
    }

    #[test]
    fn find_adapters_skips_missing_query_and_full_coverage() {
        let query_info = QueryInfo::new_blastn(&[60]);
        let scores = ScoringOptions::new_blastn();

        let mut short_query_saved = vec![Some(chain(30, 1, 1, 24, 24))];
        assert_eq!(
            s_find_adapters(&mut short_query_saved, &[3; 10], &query_info, &scores),
            0
        );
        assert_eq!(short_query_saved[0].as_ref().unwrap().adapter, -1);

        let mut full_coverage = Some(chain(40, 1, 2, 0, 0));
        full_coverage
            .as_mut()
            .unwrap()
            .hsps
            .as_mut()
            .unwrap()
            .hsp
            .query_end = 58;
        let mut query = vec![3; 122];
        query[30..42].copy_from_slice(&[0, 2, 0, 3, 1, 2, 2, 0, 0, 2, 0, 2]);
        let mut full_saved = vec![full_coverage];
        assert_eq!(
            s_find_adapters(&mut full_saved, &query, &query_info, &scores),
            0
        );
        assert_eq!(full_saved[0].as_ref().unwrap().adapter, -1);
    }

    #[test]
    fn merge_hsps_stitches_extents_and_gap_script() {
        let mut first = hsp(10, 1, 0, 0, 10);
        first.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![(
            GapAlignOpType::Sub,
            10,
        )]));
        let mut second = hsp(8, 1, 13, 12, 8);
        second.query_end = 21;
        second.subject_end = 20;
        second.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![(
            GapAlignOpType::Sub,
            8,
        )]));
        first.map_info = Some(crate::hspstream::BlastHSPMappingInfo {
            edits: Some(crate::gapinfo::JumperEditsBlock {
                num_edits: 0,
                edits: vec![JumperEdit {
                    query_pos: 4,
                    query_base: 1,
                    subject_base: 2,
                }],
            }),
            subject_overhangs: Some(crate::gapinfo::SequenceOverhangs {
                left_len: 0,
                right_len: 0,
                left: Some(vec![3, 2]),
                right: Some(vec![2, 3, 1]),
            }),
            left_edge: 7,
            right_edge: 8,

            flags: 0,
        });
        second.map_info = Some(crate::hspstream::BlastHSPMappingInfo {
            edits: Some(crate::gapinfo::JumperEditsBlock {
                num_edits: 0,
                edits: vec![JumperEdit {
                    query_pos: 16,
                    query_base: 0,
                    subject_base: 1,
                }],
            }),
            subject_overhangs: Some(crate::gapinfo::SequenceOverhangs {
                left_len: 0,
                right_len: 0,
                left: Some(vec![0]),
                right: Some(vec![1, 0]),
            }),
            left_edge: 9,
            right_edge: 10,

            flags: 0,
        });

        let scores = ScoringOptions::new_blastn();
        let query = vec![
            0, 1, 2, 3, 1, 1, 1, 1, 1, 1, 0, 1, 2, 3, 0, 1, 0, 1, 2, 3, 0,
        ];
        let merged = s_merge_hsps(&first, &second, Some(&query), Some(&scores)).expect("merged");
        assert_eq!(merged.query_offset, 0);
        assert_eq!(merged.query_end, 21);
        assert_eq!(merged.subject_offset, 0);
        assert_eq!(merged.subject_end, 20);
        assert_eq!(
            merged.edit_script.unwrap().ops_vec(),
            vec![
                (GapAlignOpType::Sub, 12),
                (GapAlignOpType::Ins, 1),
                (GapAlignOpType::Sub, 8),
            ]
        );
        assert_eq!(merged.num_ident, 16);
        let map_info = merged.map_info.as_ref().unwrap();
        assert_eq!(map_info.left_edge, 7);
        assert_eq!(map_info.right_edge, 10);
        assert_eq!(
            map_info.subject_overhangs.as_ref().unwrap().left,
            Some(vec![3, 2])
        );
        assert_eq!(
            map_info.subject_overhangs.as_ref().unwrap().right,
            Some(vec![1, 0])
        );
        assert_eq!(
            map_info.edits.as_ref().unwrap().edits,
            vec![
                JumperEdit {
                    query_pos: 4,
                    query_base: 1,
                    subject_base: 2,
                },
                JumperEdit {
                    query_pos: 10,
                    query_base: 0,
                    subject_base: 2,
                },
                JumperEdit {
                    query_pos: 11,
                    query_base: 1,
                    subject_base: 3,
                },
                JumperEdit {
                    query_pos: 12,
                    query_base: 2,
                    subject_base: 15,
                },
                JumperEdit {
                    query_pos: 16,
                    query_base: 0,
                    subject_base: 1,
                },
            ]
        );
    }

    #[test]
    fn merge_hsps_rejects_overlapping_extents() {
        let first = hsp(10, 1, 0, 0, 10);
        let second = hsp(8, 1, 5, 12, 8);
        let scores = ScoringOptions::new_blastn();
        assert!(s_merge_hsps(&first, &second, None, Some(&scores)).is_none());

        let non_overlapping = hsp(8, 1, 12, 12, 8);
        assert!(s_merge_hsps(&first, &non_overlapping, None, None).is_none());

        let mut extreme_first = hsp(10, 1, 0, 0, 10);
        extreme_first.query_end = i32::MAX - 2;
        extreme_first.subject_end = i32::MAX - 2;
        extreme_first.edit_script = None;
        let mut extreme_second = hsp(8, 1, 1, 1, 8);
        extreme_second.query_offset = i32::MAX - 1;
        extreme_second.subject_offset = i32::MAX - 1;
        extreme_second.query_end = i32::MAX;
        extreme_second.subject_end = i32::MAX;
        extreme_second.edit_script = None;
        let merged =
            s_merge_hsps(&extreme_first, &extreme_second, None, Some(&scores)).expect("merged");
        assert_eq!(merged.query_end, i32::MAX);
        assert_eq!(merged.subject_end, i32::MAX);
        assert_eq!(merged.num_ident, i32::MAX);
    }

    #[test]
    fn intron_to_gap_replaces_adjacent_hsp_pair_with_merged_hsp() {
        let scores = ScoringOptions::new_blastn();
        let mut first = hsp(10, 0, 0, 0, 10);
        first.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![(
            GapAlignOpType::Sub,
            10,
        )]));
        let mut second = hsp(8, 0, 12, 12, 8);
        second.query_end = 20;
        second.subject_end = 20;
        second.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![(
            GapAlignOpType::Sub,
            8,
        )]));
        let mut third_slot = Some(hsp(5, 0, 30, 30, 5));
        let mut second_slot = Some(second);
        let mut head = HSPContainer {
            hsp: first,
            next: hsp_container_new(&mut second_slot),
        };
        head.next.as_mut().unwrap().next = hsp_container_new(&mut third_slot);

        assert_eq!(s_intron_to_gap(&mut head, Some(&[0; 40]), Some(&scores)), 0);

        assert_eq!(head.hsp.query_offset, 0);
        assert_eq!(head.hsp.query_end, 20);
        assert_eq!(head.hsp.subject_end, 20);
        assert_eq!(
            head.hsp.edit_script.as_ref().unwrap().ops_vec(),
            vec![(GapAlignOpType::Sub, 20)]
        );
        assert_eq!(head.next.as_ref().unwrap().hsp.query_offset, 30);
        assert!(head.next.as_ref().unwrap().next.is_none());
    }

    #[test]
    fn intron_to_gap_trims_overlap_before_merging() {
        let scores = ScoringOptions::new_blastn();
        let mut first = hsp(10, 0, 0, 0, 10);
        first.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![(
            GapAlignOpType::Sub,
            10,
        )]));
        let mut second = hsp(8, 0, 8, 10, 8);
        second.query_end = 18;
        second.subject_end = 20;
        second.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![(
            GapAlignOpType::Sub,
            10,
        )]));
        let mut second_slot = Some(second);
        let mut head = HSPContainer {
            hsp: first,
            next: hsp_container_new(&mut second_slot),
        };

        assert_eq!(s_intron_to_gap(&mut head, Some(&[0; 40]), Some(&scores)), 0);

        assert_eq!(head.hsp.query_end, 18);
        assert_eq!(head.hsp.subject_end, 20);
        assert!(head.next.is_none());
        assert_eq!(
            head.hsp.edit_script.as_ref().unwrap().ops_vec(),
            vec![
                (GapAlignOpType::Sub, 10),
                (GapAlignOpType::Del, 2),
                (GapAlignOpType::Sub, 8)
            ]
        );
    }

    #[test]
    fn intron_to_gap_rejects_missing_state_without_unlinking_next_hsp() {
        let scores = ScoringOptions::new_blastn();
        let mut first = hsp(10, 0, 0, 0, 10);
        first.query_end = 10;
        first.subject_end = 10;
        let mut second = Some(hsp(10, 0, 14, 14, 10));
        let mut head = HSPContainer {
            hsp: first,
            next: hsp_container_new(&mut second),
        };

        assert_eq!(s_intron_to_gap(&mut head, None, Some(&scores)), -1);
        assert!(head.next.is_some());
        assert_eq!(head.next.as_ref().unwrap().hsp.query_offset, 14);

        assert_eq!(s_intron_to_gap(&mut head, Some(&[0; 40]), None), -1);
        assert!(head.next.is_some());

        let mut single = HSPContainer {
            hsp: hsp(10, 0, 0, 0, 10),
            next: None,
        };
        assert_eq!(
            s_intron_to_gap(&mut single, Some(&[0; 40]), Some(&scores)),
            -1
        );
        assert!(single.next.is_none());
    }

    #[test]
    fn find_splice_junctions_merges_short_gap_pair_and_recomputes_score() {
        let scores = ScoringOptions::new_blastn();
        let mut first = hsp(12, 0, 0, 0, 10);
        first.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![(
            GapAlignOpType::Sub,
            10,
        )]));
        first.map_info = Some(crate::hspstream::BlastHSPMappingInfo {
            edits: None,
            subject_overhangs: Some(crate::gapinfo::SequenceOverhangs {
                left_len: 0,
                right_len: 0,
                left: None,
                right: Some(vec![0, 0, 0]),
            }),
            left_edge: 0,
            right_edge: 0,

            flags: 0,
        });
        let mut second = hsp(10, 0, 12, 12, 8);
        second.query_end = 20;
        second.subject_end = 20;
        second.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![(
            GapAlignOpType::Sub,
            8,
        )]));
        let mut second_slot = Some(second);
        let mut chain = hsp_chain_new(0);
        chain.hsps = Some(Box::new(HSPContainer {
            hsp: first,
            next: hsp_container_new(&mut second_slot),
        }));
        chain.count = 2;

        assert_eq!(
            s_find_splice_junctions(Some(&mut chain), Some(&[0; 40]), 40, Some(&scores)),
            0
        );

        let head = chain.hsps.as_ref().unwrap();
        assert_eq!(head.hsp.query_offset, 0);
        assert_eq!(head.hsp.query_end, 20);
        assert_eq!(head.hsp.subject_end, 20);
        assert!(head.next.is_none());
        assert_eq!(
            head.hsp.edit_script.as_ref().unwrap().ops_vec(),
            vec![
                (GapAlignOpType::Sub, 10),
                (GapAlignOpType::Ins, 2),
                (GapAlignOpType::Del, 2),
                (GapAlignOpType::Sub, 8)
            ]
        );
        assert_eq!(chain.score, 0);
    }

    #[test]
    fn find_splice_junctions_trims_query_overlap_when_gap_state_is_unavailable() {
        let scores = ScoringOptions::new_blastn();
        let mut first = hsp(60, 0, 0, 0, 10);
        first.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![(
            GapAlignOpType::Sub,
            10,
        )]));
        let mut second = hsp(60, 0, 8, 70, 10);
        second.query_end = 18;
        second.subject_end = 80;
        second.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![(
            GapAlignOpType::Sub,
            10,
        )]));
        let mut second_slot = Some(second);
        let mut chain = hsp_chain_new(0);
        chain.hsps = Some(Box::new(HSPContainer {
            hsp: first,
            next: hsp_container_new(&mut second_slot),
        }));
        chain.count = 2;

        assert_eq!(
            s_find_splice_junctions(Some(&mut chain), Some(&[0; 40]), 40, Some(&scores)),
            0
        );

        let head = chain.hsps.as_ref().unwrap();
        let second = head.next.as_ref().unwrap();
        assert_eq!(head.hsp.query_end, 10);
        assert_eq!(second.hsp.query_offset, 10);
        assert_eq!(second.hsp.subject_offset, 72);
        assert_eq!(
            second.hsp.edit_script.as_ref().unwrap().ops_vec(),
            vec![(GapAlignOpType::Sub, 8)]
        );
    }

    #[test]
    fn find_splice_junctions_rejects_missing_state_and_scores_empty_chain() {
        let scores = ScoringOptions::new_blastn();
        let mut empty_chain = hsp_chain_new(0);
        empty_chain.score = 42;
        empty_chain.next = Some(chain(30, 0, 1, 0, 0));

        assert_eq!(
            s_find_splice_junctions(None, Some(&[0; 40]), 40, Some(&scores)),
            -1
        );
        assert_eq!(
            s_find_splice_junctions(Some(&mut empty_chain), None, 40, Some(&scores)),
            -1
        );
        assert_eq!(
            s_find_splice_junctions(Some(&mut empty_chain), Some(&[0; 40]), 40, None),
            -1
        );
        assert_eq!(
            s_find_splice_junctions(Some(&mut empty_chain), Some(&[0; 40]), 40, Some(&scores)),
            0
        );
        assert_eq!(empty_chain.score, 0);
        assert_eq!(empty_chain.next.as_ref().unwrap().score, 30);

        let mut extreme_first = hsp(60, 0, 0, 0, 1);
        extreme_first.query_end = i32::MAX;
        extreme_first.subject_end = i32::MAX;
        let mut extreme_second = hsp(60, 0, 0, 0, 1);
        extreme_second.query_offset = i32::MIN;
        extreme_second.subject_offset = i32::MIN;
        extreme_second.query_end = 0;
        extreme_second.subject_end = 0;
        let mut extreme_first_slot = Some(extreme_first);
        let mut extreme_second_slot = Some(extreme_second);
        let mut extreme_chain = hsp_chain_new(0);
        extreme_chain.hsps = hsp_container_new(&mut extreme_first_slot);
        extreme_chain.hsps.as_mut().unwrap().next = hsp_container_new(&mut extreme_second_slot);
        assert_eq!(
            s_find_splice_junctions(Some(&mut extreme_chain), Some(&[0; 40]), 40, Some(&scores)),
            -1
        );
    }

    #[test]
    fn find_best_pairs_marks_best_convergent_pair() {
        let mut first_list = Some(chain(30, 1, 1, 0, 100));
        let mut second_list = Some(chain(25, -1, 1, 0, 180));
        {
            let hsp = &mut second_list.as_mut().unwrap().hsps.as_mut().unwrap().hsp;
            hsp.subject_end = 200;
        }
        let mut pair_info = Vec::new();

        assert!(s_find_best_pairs(
            &mut first_list,
            &mut second_list,
            0,
            &mut pair_info,
            false,
            Some(&ScoringOptions::new_blastn()),
            Some(&[0; 240]),
        ));

        let first = first_list.as_ref().unwrap();
        let second = second_list.as_ref().unwrap();
        assert_eq!(first.pair, Some(0));
        assert_eq!(second.pair, Some(0));
        assert_eq!(first.pair_score, Some(25));
        assert_eq!(second.pair_score, Some(30));
        assert_eq!(first.pair_conf, PAIR_CONVERGENT);
        assert_eq!(second.pair_conf, PAIR_CONVERGENT);
        assert_eq!(pair_info[0].conf, PAIR_CONVERGENT);
        assert_eq!(pair_info[0].distance, 100);
    }

    #[test]
    fn compare_pairs_matches_find_best_pairs_sort_order() {
        let best = PairInfo {
            first: 0,
            second: 0,
            score: 60,
            trim_first: 0,
            trim_second: 0,
            valid_pair: false,
            distance: 120,
            conf: PAIR_CONVERGENT,
        };
        let lower_score = PairInfo {
            score: 59,
            ..best.clone()
        };
        let lower_conf = PairInfo {
            score: 60,
            conf: PAIR_PARALLEL,
            ..best.clone()
        };
        let farther = PairInfo {
            score: 60,
            distance: 150,
            ..best.clone()
        };

        assert_eq!(s_compare_pairs(&best, &best), 0);
        assert!(s_compare_pairs(&best, &lower_score) < 0);
        assert!(s_compare_pairs(&lower_score, &best) > 0);
        assert!(s_compare_pairs(&best, &lower_conf) < 0);
        assert!(s_compare_pairs(&best, &farther) < 0);

        let mut pairs = vec![farther, lower_conf, lower_score, best.clone()];
        pairs.sort_by(|a, b| s_compare_pairs(a, b).cmp(&0));
        assert_eq!(pairs[0], best);
    }

    #[test]
    fn find_best_pairs_trims_overextended_plus_chain_end() {
        let mut plus = hsp(120, 1, 0, 100, 120);
        plus.query_end = 120;
        plus.subject_end = 220;
        plus.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![(
            GapAlignOpType::Sub,
            120,
        )]));
        let mut first_slot = Some(plus);
        let mut first_list = Some(hsp_chain_new(1));
        first_list.as_mut().unwrap().score = 120;
        first_list.as_mut().unwrap().hsps = hsp_container_new(&mut first_slot);

        let mut minus = hsp(20, -1, 0, 180, 20);
        minus.query_end = 20;
        minus.subject_end = 200;
        minus.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![(
            GapAlignOpType::Sub,
            20,
        )]));
        let mut second_slot = Some(minus);
        let mut second_list = Some(hsp_chain_new(-1));
        second_list.as_mut().unwrap().score = 20;
        second_list.as_mut().unwrap().hsps = hsp_container_new(&mut second_slot);

        let mut pair_info = Vec::new();
        assert!(s_find_best_pairs(
            &mut first_list,
            &mut second_list,
            0,
            &mut pair_info,
            false,
            Some(&ScoringOptions::new_blastn()),
            Some(&[0; 240]),
        ));

        let first = first_list.as_ref().unwrap();
        let hsp = &first.hsps.as_ref().unwrap().hsp;
        assert_eq!(pair_info[0].trim_first, 200);
        assert_eq!(hsp.subject_end, 200);
        assert_eq!(hsp.query_end, 100);
        assert_eq!(
            hsp.edit_script.as_ref().unwrap().ops_vec(),
            vec![(GapAlignOpType::Sub, 100)]
        );
        assert_eq!(first.score, 100);
    }

    #[test]
    fn find_best_pairs_preserves_lists_when_no_candidate_passes_cutoff() {
        let mut first_list = Some(chain(20, 1, 1, 0, 100));
        let mut second_list = Some(chain(19, -1, 1, 0, 180));
        let mut pair_info = Vec::new();

        assert!(!s_find_best_pairs(
            &mut first_list,
            &mut second_list,
            100,
            &mut pair_info,
            false,
            None,
            None,
        ));

        assert!(pair_info.is_empty());
        assert_eq!(first_list.as_ref().unwrap().score, 20);
        assert_eq!(second_list.as_ref().unwrap().score, 19);
        assert_eq!(first_list.as_ref().unwrap().pair, None);
        assert_eq!(second_list.as_ref().unwrap().pair, None);
    }

    #[test]
    fn find_best_pairs_marks_parallel_same_strand_pair_without_trimming() {
        let mut first_list = Some(chain(30, 1, 1, 0, 100));
        let mut second_list = Some(chain(25, 1, 1, 0, 180));
        let mut pair_info = Vec::new();

        assert!(s_find_best_pairs(
            &mut first_list,
            &mut second_list,
            0,
            &mut pair_info,
            false,
            None,
            None,
        ));

        assert_eq!(pair_info[0].conf, PAIR_PARALLEL);
        assert_eq!(pair_info[0].score, 54);
        assert_eq!(pair_info[0].trim_first, 0);
        assert_eq!(pair_info[0].trim_second, 0);
        assert_eq!(first_list.as_ref().unwrap().pair_conf, PAIR_PARALLEL);
        assert_eq!(second_list.as_ref().unwrap().pair_conf, PAIR_PARALLEL);
    }

    #[test]
    fn spliced_paired_run_returns_zero_for_absent_hsp_list() {
        let mut data = BlastHSPMapperData::default();
        assert_eq!(s_blast_hsp_mapper_spliced_paired_run(&mut data, None), 0);
        assert!(data.saved_chains.is_empty());
    }

    #[test]
    fn spliced_paired_run_rejects_missing_required_state() {
        let mut hsp_list = HspList::new(7);
        hsp_list.hsps.push(hsp(30, 0, 0, 100, 10));

        let mut missing_params = BlastHSPMapperData {
            query_info: Some(QueryInfo::new_blastn(&[20])),
            ..Default::default()
        };
        assert_eq!(
            s_blast_hsp_mapper_spliced_paired_run(&mut missing_params, Some(hsp_list.clone())),
            -1
        );
        assert!(missing_params.saved_chains.is_empty());

        let scores = ScoringOptions::new_blastn();
        let mut params =
            blast_hsp_mapper_params_new(Some(&HitSavingOptions::default()), Some(&scores)).unwrap();
        params.splice = false;
        let mut missing_query_info = BlastHSPMapperData {
            params: Some(params.clone()),
            ..Default::default()
        };
        assert_eq!(
            s_blast_hsp_mapper_spliced_paired_run(&mut missing_query_info, Some(hsp_list.clone())),
            -1
        );
        assert!(missing_query_info.saved_chains.is_empty());

        params.splice = true;
        let mut missing_splice_query = BlastHSPMapperData {
            params: Some(params),
            query_info: Some(QueryInfo::new_blastn(&[20])),
            query: None,
            saved_chains: vec![None],
        };
        assert_eq!(
            s_blast_hsp_mapper_spliced_paired_run(&mut missing_splice_query, Some(hsp_list),),
            -1
        );
        assert!(missing_splice_query
            .saved_chains
            .iter()
            .all(Option::is_none));
    }

    #[test]
    fn spliced_paired_run_skips_invalid_contexts_duplicates_and_cutoff_failures() {
        let query_info = QueryInfo::new_blastn(&[50]);
        let scores = ScoringOptions::new_blastn();
        let mut params =
            blast_hsp_mapper_params_new(Some(&HitSavingOptions::default()), Some(&scores)).unwrap();
        params.cutoff_score = 25;
        params.cutoff_edit_dist = -1;

        let mut valid = hsp(40, 0, 0, 5, 10);
        valid.query_frame = 0;
        let duplicate = valid.clone();
        let low_score = hsp(20, 0, 10, 20, 10);
        let invalid_context = hsp(50, -1, 0, 30, 10);
        let out_of_range_context = hsp(60, 99, 0, 40, 10);

        let mut hsp_list = HspList::new(11);
        hsp_list.hsps.push(low_score);
        hsp_list.hsps.push(out_of_range_context);
        hsp_list.hsps.push(duplicate);
        hsp_list.hsps.push(invalid_context);
        hsp_list.hsps.push(valid);

        let mut data = BlastHSPMapperData {
            params: Some(params),
            query: Some(vec![0; 100]),
            query_info: Some(query_info),
            saved_chains: Vec::new(),
        };

        assert_eq!(
            s_blast_hsp_mapper_spliced_paired_run(&mut data, Some(hsp_list)),
            0
        );

        assert_eq!(data.saved_chains.len(), 1);
        let saved = data.saved_chains[0].as_ref().unwrap();
        assert_eq!(saved.score, 40);
        assert_eq!(saved.oid, 11);
        assert_eq!(saved.hsps.as_ref().unwrap().hsp.query_frame, 1);
        assert!(saved.next.is_none());
    }

    #[test]
    fn spliced_paired_run_uses_cutoff_function_when_present() {
        let query_info = QueryInfo::new_blastn(&[50]);
        let scores = ScoringOptions::new_blastn();
        let mut params =
            blast_hsp_mapper_params_new(Some(&HitSavingOptions::default()), Some(&scores)).unwrap();
        params.cutoff_score = 1;
        params.cutoff_score_fun = [1000, 200];
        params.cutoff_edit_dist = -1;

        let mut hsp_list = HspList::new(12);
        hsp_list.hsps.push(hsp(100, 0, 0, 5, 10));

        let mut data = BlastHSPMapperData {
            params: Some(params),
            query: Some(vec![0; 100]),
            query_info: Some(query_info),
            saved_chains: Vec::new(),
        };

        assert_eq!(
            s_blast_hsp_mapper_spliced_paired_run(&mut data, Some(hsp_list)),
            0
        );
        assert_eq!(data.saved_chains.len(), 1);
        assert!(data.saved_chains[0].is_none());
    }

    #[test]
    fn spliced_paired_run_saves_and_pairs_segment_flagged_query_chains() {
        let mut query_info = QueryInfo::new_blastn(&[80, 100, 100]);
        query_info.set_query_segment_flags(1, E_FIRST_SEGMENT);
        query_info.set_query_segment_flags(2, crate::queryinfo::E_LAST_SEGMENT);
        let scores = ScoringOptions::new_blastn();
        let mut params =
            blast_hsp_mapper_params_new(Some(&HitSavingOptions::default()), Some(&scores)).unwrap();
        params.cutoff_score = 1;
        params.cutoff_edit_dist = -1;

        let mut first = hsp(30, 2, 0, 100, 10);
        first.query_frame = 1;
        first.subject_end = 110;
        let mut second = hsp(25, 5, 0, 180, 10);
        second.query_frame = -1;
        second.subject_end = 200;

        let mut hsp_list = HspList::new(7);
        hsp_list.hsps.push(second);
        hsp_list.hsps.push(first);

        let mut data = BlastHSPMapperData {
            params: Some(params),
            query: Some(vec![0; 240]),
            query_info: Some(query_info),
            saved_chains: vec![None, None],
        };

        assert_eq!(
            s_blast_hsp_mapper_spliced_paired_run(&mut data, Some(hsp_list)),
            0
        );

        assert!(data.saved_chains[0].is_none());
        let first_chain = data.saved_chains[1].as_ref().unwrap();
        let second_chain = data.saved_chains[2].as_ref().unwrap();
        assert_eq!(first_chain.oid, 7);
        assert_eq!(second_chain.oid, 7);
        assert_eq!(first_chain.pair, Some(0));
        assert_eq!(second_chain.pair, Some(0));
        assert_eq!(first_chain.pair_score, Some(25));
        assert_eq!(second_chain.pair_score, Some(30));
        assert_eq!(first_chain.pair_conf, PAIR_CONVERGENT);
        assert_eq!(second_chain.pair_conf, PAIR_CONVERGENT);
    }

    #[test]
    fn spliced_paired_run_invokes_splice_junction_pass() {
        let query_info = QueryInfo::new_blastn(&[40]);
        let scores = ScoringOptions::new_blastn();
        let mut params =
            blast_hsp_mapper_params_new(Some(&HitSavingOptions::default()), Some(&scores)).unwrap();
        params.splice = true;
        params.cutoff_score = 1;
        params.cutoff_edit_dist = -1;

        let mut first = hsp(12, 0, 0, 0, 10);
        first.query_frame = 1;
        first.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![(
            GapAlignOpType::Sub,
            10,
        )]));
        first.map_info = Some(crate::hspstream::BlastHSPMappingInfo {
            edits: None,
            subject_overhangs: Some(crate::gapinfo::SequenceOverhangs {
                left_len: 0,
                right_len: 0,
                left: None,
                right: Some(vec![0, 0, 0]),
            }),
            left_edge: 0,
            right_edge: 0,

            flags: 0,
        });
        let mut second = hsp(10, 0, 12, 12, 8);
        second.query_frame = 1;
        second.query_end = 20;
        second.subject_end = 20;
        second.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![(
            GapAlignOpType::Sub,
            8,
        )]));

        let mut hsp_list = HspList::new(13);
        hsp_list.hsps.push(second);
        hsp_list.hsps.push(first);

        let mut data = BlastHSPMapperData {
            params: Some(params),
            query: Some(vec![0; 82]),
            query_info: Some(query_info),
            saved_chains: Vec::new(),
        };

        assert_eq!(
            s_blast_hsp_mapper_spliced_paired_run(&mut data, Some(hsp_list)),
            0
        );

        let saved = data.saved_chains[0].as_ref().unwrap();
        let head = saved.hsps.as_ref().unwrap();
        assert_eq!(saved.oid, 13);
        assert!(saved.score > 0);
        assert_eq!(head.hsp.query_offset, 0);
        assert_eq!(head.hsp.query_end, 20);
        assert_eq!(head.hsp.subject_end, 20);
        assert!(head.next.is_none());
        assert_eq!(
            head.hsp.edit_script.as_ref().unwrap().ops_vec(),
            vec![
                (GapAlignOpType::Sub, 10),
                (GapAlignOpType::Ins, 2),
                (GapAlignOpType::Del, 2),
                (GapAlignOpType::Sub, 8)
            ]
        );
    }

    #[test]
    fn find_rearranged_pairs_marks_parallel_adjacent_queries() {
        let mut query_info = QueryInfo::new_blastn(&[20, 20]);
        query_info.set_query_segment_flags(0, E_FIRST_SEGMENT);
        query_info.set_query_segment_flags(1, crate::queryinfo::E_LAST_SEGMENT);
        let mut saved = vec![
            Some(chain(30, 1, 100, 0, 0)),
            Some(chain(25, -1, 200, 0, 0)),
        ];

        assert_eq!(s_find_rearranged_pairs(&mut saved, &query_info), 0);
        assert_eq!(saved[0].as_ref().unwrap().pair, Some(1));
        assert_eq!(saved[1].as_ref().unwrap().pair, Some(0));
        assert_eq!(saved[0].as_ref().unwrap().pair_conf, PAIR_PARALLEL);
        assert_eq!(saved[1].as_ref().unwrap().pair_conf, PAIR_PARALLEL);
    }

    #[test]
    fn find_rearranged_pairs_uses_segment_flags_not_even_query_parity() {
        let mut query_info = QueryInfo::new_blastn(&[20, 20, 20]);
        query_info.set_query_segment_flags(1, E_FIRST_SEGMENT);
        query_info.set_query_segment_flags(2, crate::queryinfo::E_LAST_SEGMENT);
        let mut saved = vec![
            Some(chain(30, 1, 100, 0, 0)),
            Some(chain(31, 1, 101, 0, 0)),
            Some(chain(25, -1, 200, 0, 0)),
        ];

        assert_eq!(s_find_rearranged_pairs(&mut saved, &query_info), 0);
        assert_eq!(saved[0].as_ref().unwrap().pair, None);
        assert_eq!(saved[1].as_ref().unwrap().pair, Some(2));
        assert_eq!(saved[2].as_ref().unwrap().pair, Some(1));
    }

    #[test]
    fn find_rearranged_pairs_skips_same_oid_or_same_strand() {
        let mut query_info = QueryInfo::new_blastn(&[20, 20]);
        query_info.set_query_segment_flags(0, E_FIRST_SEGMENT);
        query_info.set_query_segment_flags(1, crate::queryinfo::E_LAST_SEGMENT);
        let mut same_oid = vec![
            Some(chain(30, 1, 100, 0, 0)),
            Some(chain(25, -1, 100, 0, 0)),
        ];
        assert_eq!(s_find_rearranged_pairs(&mut same_oid, &query_info), 0);
        assert_eq!(same_oid[0].as_ref().unwrap().pair, None);
        assert_eq!(same_oid[1].as_ref().unwrap().pair, None);

        let mut same_strand = vec![Some(chain(30, 1, 100, 0, 0)), Some(chain(25, 1, 200, 0, 0))];
        assert_eq!(s_find_rearranged_pairs(&mut same_strand, &query_info), 0);
        assert_eq!(same_strand[0].as_ref().unwrap().pair, None);
        assert_eq!(same_strand[1].as_ref().unwrap().pair, None);
    }

    #[test]
    fn prune_chains_keeps_best_single_and_sets_counts() {
        let mut list = Some(chain(30, 1, 1, 0, 0));
        let mut low = Some(chain(20, 1, 2, 0, 0));
        low.as_mut().unwrap().next = list.take();
        let mut saved = vec![low];

        assert_eq!(s_prune_chains(&mut saved, 1, 5), 0);
        let kept = saved[0].as_ref().unwrap();
        assert_eq!(kept.score, 30);
        assert_eq!(kept.count, 1);
        assert!(kept.next.is_none());
    }

    #[test]
    fn prune_chains_preserves_best_pair_with_bonus() {
        let mut paired = chain(28, 1, 1, 0, 0);
        paired.pair = Some(1);
        paired.pair_score = Some(25);
        let mut better_single = chain(30, 1, 2, 0, 0);
        better_single.next = Some(paired);

        let mut mate = chain(25, -1, 3, 0, 0);
        mate.pair = Some(0);
        mate.pair_score = Some(28);
        let mut saved = vec![Some(better_single), Some(mate)];

        assert_eq!(s_prune_chains(&mut saved, 2, 5), 0);
        let kept = saved[0].as_ref().unwrap();
        assert_eq!(kept.score, 28);
        assert_eq!(kept.pair, Some(1));
        assert_eq!(kept.count, 2);
        assert!(kept.next.is_none());
        assert_eq!(saved[1].as_ref().unwrap().count, 1);
    }

    #[test]
    fn prune_chains_clears_mate_when_partner_is_pruned_like_hsp_chain_free() {
        let mut pruned_partner = chain(10, 1, 1, 0, 0);
        pruned_partner.pair = Some(1);
        pruned_partner.pair_score = Some(100);
        let mut better_single = chain(30, 1, 2, 0, 0);
        better_single.next = Some(pruned_partner);

        let mut mate = chain(100, -1, 3, 0, 0);
        mate.pair = Some(0);
        mate.pair_score = Some(10);
        let mut saved = vec![Some(better_single), Some(mate)];

        assert_eq!(s_prune_chains(&mut saved, 2, 5), 0);
        let kept_single = saved[0].as_ref().unwrap();
        assert_eq!(kept_single.score, 30);
        assert!(kept_single.next.is_none());

        let orphaned_mate = saved[1].as_ref().unwrap();
        assert_eq!(orphaned_mate.score, 100);
        assert_eq!(orphaned_mate.pair, None);
        assert_eq!(orphaned_mate.pair_score, None);
        assert_eq!(orphaned_mate.pair_conf, PAIR_NONE);
    }

    #[test]
    fn prune_chains_tolerates_empty_and_out_of_range_query_counts() {
        let mut empty: Vec<Option<Box<HSPChain>>> = Vec::new();
        assert_eq!(s_prune_chains(&mut empty, 3, 5), 0);
        assert!(empty.is_empty());

        let mut saved = vec![Some(chain(22, 1, 1, 0, 0)), None];
        assert_eq!(s_prune_chains(&mut saved, -1, 5), 0);
        assert_eq!(saved[0].as_ref().unwrap().count, 0);

        assert_eq!(s_prune_chains(&mut saved, 10, 5), 0);
        assert_eq!(saved[0].as_ref().unwrap().count, 1);
        assert!(saved[1].is_none());
    }

    #[test]
    fn remove_overlaps_splits_adjacent_query_overlap_into_new_chain() {
        let mut list = Some(chain(30, 1, 1, 0, 0));
        let mut second = hsp(20, 1, 8, 20, 10);
        second.query_end = 18;
        second.subject_end = 30;
        let mut second_slot = Some(second);
        list.as_mut().unwrap().hsps.as_mut().unwrap().next = hsp_container_new(&mut second_slot);
        list.as_mut().unwrap().pair = Some(1);
        list.as_mut().unwrap().pair_score = Some(40);

        let scores = ScoringOptions::new_blastn();
        assert_eq!(s_remove_overlaps(&mut list, Some(&scores), 40), 0);

        let first = list.as_ref().unwrap();
        assert_eq!(first.hsps.as_ref().unwrap().hsp.query_offset, 0);
        assert!(first.hsps.as_ref().unwrap().next.is_none());
        assert_eq!(first.pair, Some(1));

        let second_chain = first.next.as_ref().unwrap();
        assert_eq!(second_chain.hsps.as_ref().unwrap().hsp.query_offset, 8);
        assert!(second_chain.hsps.as_ref().unwrap().next.is_none());
        assert_eq!(second_chain.pair, None);
        assert_eq!(second_chain.pair_score, None);
        assert!(second_chain.next.is_none());
    }

    #[test]
    fn remove_overlaps_leaves_non_overlapping_chain_intact() {
        let mut list = Some(chain(30, 1, 1, 0, 0));
        let mut second = hsp(20, 1, 12, 20, 10);
        second.query_end = 22;
        second.subject_end = 30;
        let mut second_slot = Some(second);
        list.as_mut().unwrap().hsps.as_mut().unwrap().next = hsp_container_new(&mut second_slot);

        let scores = ScoringOptions::new_blastn();
        assert_eq!(s_remove_overlaps(&mut list, Some(&scores), 40), 0);

        let first = list.as_ref().unwrap();
        assert!(first.hsps.as_ref().unwrap().next.is_some());
        assert!(first.next.is_none());
    }

    #[test]
    fn remove_overlaps_reports_invalid_inputs_without_mutating_chain() {
        let mut empty = None;
        let scores = ScoringOptions::new_blastn();
        assert_eq!(s_remove_overlaps(&mut empty, Some(&scores), 40), -1);
        assert!(empty.is_none());

        let mut list = Some(chain(30, 1, 1, 0, 0));
        let original_score = list.as_ref().unwrap().score;
        let original_pair = list.as_ref().unwrap().pair;
        let original_query_offset = list
            .as_ref()
            .unwrap()
            .hsps
            .as_ref()
            .unwrap()
            .hsp
            .query_offset;

        assert_eq!(s_remove_overlaps(&mut list, None, 40), -1);
        let chain = list.as_ref().unwrap();
        assert_eq!(chain.score, original_score);
        assert_eq!(chain.pair, original_pair);
        assert_eq!(
            chain.hsps.as_ref().unwrap().hsp.query_offset,
            original_query_offset
        );
        assert!(chain.next.is_none());
    }

    #[test]
    fn overlap_list_helpers_handle_empty_and_no_hsp_chains() {
        let scores = ScoringOptions::new_blastn();
        let mut empty = None;
        assert_eq!(remove_overlaps_from_chain_list(&mut empty, &scores, 40), 0);
        assert!(empty.is_none());

        let mut no_hsps = hsp_chain_new(1);
        no_hsps.oid = 9;
        no_hsps.score = 17;
        let mut list = Some(no_hsps);
        assert_eq!(remove_overlaps_from_chain_list(&mut list, &scores, 40), 0);
        let chain = list.as_ref().unwrap();
        assert_eq!(chain.oid, 9);
        assert_eq!(chain.score, 17);
        assert!(chain.hsps.is_none());
        assert!(chain.next.is_none());
    }

    #[test]
    fn finalize_sorts_filters_counts_and_moves_saved_chains_to_results() {
        let query_info = QueryInfo::new_blastn(&[12]);
        let query = vec![3; 100];
        let scores = ScoringOptions::new_blastn();

        let mut head = chain(30, 1, 2, 0, 0);
        head.next = Some(chain(30, 1, 1, 0, 0));
        let mut saved = vec![Some(head)];
        let mut results = blast_mapping_results_new();

        assert_eq!(
            s_finalize(
                &mut saved,
                &mut results,
                &query_info,
                &query,
                &scores,
                true,
                20,
                -1,
            ),
            0
        );

        assert!(saved.is_empty());
        assert_eq!(results.num_queries, 1);
        let first = results.chain_array[0].as_ref().unwrap();
        assert_eq!(first.oid, 1);
        assert_eq!(first.count, 2);
        let second = first.next.as_ref().unwrap();
        assert_eq!(second.oid, 2);
        assert_eq!(second.count, 2);
    }

    #[test]
    fn finalize_helpers_sort_by_oid_and_count_unique_fragments() {
        let mut list = Some(chain(30, 1, 3, 10, 30));
        list.as_mut().unwrap().next = Some(chain(30, 1, 1, 0, 10));
        list.as_mut().unwrap().next.as_mut().unwrap().next = Some(chain(30, 1, 1, 0, 10));
        list.as_mut()
            .unwrap()
            .next
            .as_mut()
            .unwrap()
            .next
            .as_mut()
            .unwrap()
            .next = Some(chain(30, 1, 2, 5, 20));

        sort_chain_list_by_oid(&mut list);
        let mut cursor = list.as_deref();
        let mut order = Vec::new();
        while let Some(chain) = cursor {
            order.push((chain.oid, s_find_fragment_start(Some(chain))));
            cursor = chain.next.as_deref();
        }
        assert_eq!(order, vec![(1, 10), (1, 10), (2, 20), (3, 30)]);

        set_unique_mapping_counts(&mut list);
        let mut cursor = list.as_deref();
        while let Some(chain) = cursor {
            assert_eq!(chain.count, 3);
            cursor = chain.next.as_deref();
        }

        let mut empty = None;
        set_unique_mapping_counts(&mut empty);
        assert!(empty.is_none());
    }

    #[test]
    fn blast_hsp_mapper_final_delegates_finalize_and_clears_saved_chains() {
        let query_info = QueryInfo::new_blastn(&[12]);
        let query = vec![3; 100];
        let scores = ScoringOptions::new_blastn();
        let mut params =
            blast_hsp_mapper_params_new(Some(&HitSavingOptions::default()), Some(&scores)).unwrap();
        params.cutoff_score = 20;
        params.cutoff_edit_dist = -1;

        let mut data = BlastHSPMapperData {
            params: Some(params),
            query: Some(query),
            query_info: Some(query_info),
            saved_chains: vec![Some(chain(30, 1, 1, 0, 0))],
        };
        let mut results = blast_mapping_results_new();

        assert_eq!(s_blast_hsp_mapper_final(&mut data, &mut results), 0);
        assert!(data.saved_chains.is_empty());
        assert_eq!(results.num_queries, 1);
        assert_eq!(results.chain_array[0].as_ref().unwrap().oid, 1);
    }

    #[test]
    fn blast_hsp_mapper_final_runs_adapter_and_poly_a_passes() {
        let query_info = QueryInfo::new_blastn(&[60, 20]);
        let mut query = vec![3; 164];
        query[30..42].copy_from_slice(&[0, 2, 0, 3, 1, 2, 2, 0, 0, 2, 0, 2]);
        query[143 + 14..143 + 20].copy_from_slice(&[0, 0, 0, 0, 0, 0]);

        let mut adapter_chain = chain(30, 1, 1, 24, 24);
        {
            let hsp = &mut adapter_chain.hsps.as_mut().unwrap().hsp;
            hsp.query_end = 42;
            hsp.subject_end = 42;
        }
        let poly_chain = chain(30, -1, 2, 2, 0);

        let scores = ScoringOptions::new_blastn();
        let mut params =
            blast_hsp_mapper_params_new(Some(&HitSavingOptions::default()), Some(&scores)).unwrap();
        params.cutoff_score = i32::MIN;
        params.cutoff_edit_dist = -1;

        let mut data = BlastHSPMapperData {
            params: Some(params),
            query: Some(query),
            query_info: Some(query_info),
            saved_chains: vec![Some(adapter_chain), Some(poly_chain)],
        };
        let mut results = blast_mapping_results_new();

        assert_eq!(s_blast_hsp_mapper_final(&mut data, &mut results), 0);
        assert!(data.saved_chains.is_empty());

        let adapter = results.chain_array[0].as_ref().unwrap();
        assert_eq!(adapter.adapter, 30);
        let adapter_hsp = &adapter.hsps.as_ref().unwrap().hsp;
        assert_eq!(adapter_hsp.query_offset, 24);
        assert_eq!(adapter_hsp.query_end, 30);
        let adapter_edge = adapter_hsp.map_info.as_ref().unwrap().right_edge;
        assert_eq!(adapter_edge & MAPPER_ADAPTER, MAPPER_ADAPTER);
        assert_eq!(adapter_edge & MAPPER_EXON, MAPPER_EXON);

        let poly = results.chain_array[1].as_ref().unwrap();
        assert_eq!(poly.poly_a, 14);
        let poly_edge = last_container(poly)
            .unwrap()
            .hsp
            .map_info
            .as_ref()
            .unwrap()
            .right_edge;
        assert_eq!(poly_edge & MAPPER_POLY_A, MAPPER_POLY_A);
        assert_eq!(poly_edge & MAPPER_EXON, MAPPER_EXON);
    }

    #[test]
    fn blast_hsp_mapper_final_rejects_missing_required_state() {
        let mut empty = BlastHSPMapperData::default();
        let mut results = blast_mapping_results_new();
        assert_eq!(s_blast_hsp_mapper_final(&mut empty, &mut results), 0);
        assert_eq!(results.num_queries, 0);

        let query_info = QueryInfo::new_blastn(&[12]);
        let query = vec![3; 100];
        let scores = ScoringOptions::new_blastn();
        let params =
            blast_hsp_mapper_params_new(Some(&HitSavingOptions::default()), Some(&scores)).unwrap();

        let mut missing_params = BlastHSPMapperData {
            params: None,
            query: Some(query.clone()),
            query_info: Some(query_info.clone()),
            saved_chains: vec![Some(chain(30, 1, 1, 0, 0))],
        };
        assert_eq!(
            s_blast_hsp_mapper_final(&mut missing_params, &mut results),
            -1
        );
        assert!(missing_params.saved_chains[0].is_some());

        let mut missing_query_info = BlastHSPMapperData {
            params: Some(params.clone()),
            query: Some(query.clone()),
            query_info: None,
            saved_chains: vec![Some(chain(30, 1, 1, 0, 0))],
        };
        assert_eq!(
            s_blast_hsp_mapper_final(&mut missing_query_info, &mut results),
            -1
        );
        assert!(missing_query_info.saved_chains[0].is_some());

        let mut missing_query = BlastHSPMapperData {
            params: Some(params),
            query: None,
            query_info: Some(query_info),
            saved_chains: vec![Some(chain(30, 1, 1, 0, 0))],
        };
        assert_eq!(
            s_blast_hsp_mapper_final(&mut missing_query, &mut results),
            -1
        );
        assert!(missing_query.saved_chains[0].is_some());
    }

    #[test]
    fn find_splice_junctions_for_overlaps_clears_splice_bits_without_trimming_second_hsp() {
        let mut first = hsp(30, 0, 0, 10, 20);
        first.query_end = 20;
        first.subject_end = 30;
        first.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![(
            GapAlignOpType::Sub,
            20,
        )]));

        let mut second = hsp(20, 0, 15, 40, 15);
        second.query_end = 30;
        second.subject_end = 55;
        second.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![(
            GapAlignOpType::Sub,
            15,
        )]));
        second.map_info = Some(crate::hspstream::BlastHSPMappingInfo {
            edits: None,
            subject_overhangs: None,
            left_edge: MAPPER_SPLICE_SIGNAL | 1,
            right_edge: 0,

            flags: 0,
        });

        assert_eq!(
            s_find_splice_junctions_for_overlaps(&mut first, &mut second, None, 30, false),
            0
        );
        assert_eq!(first.query_end, 20);
        assert_eq!(second.query_offset, 15);
        assert_eq!(second.subject_offset, 40);
        assert_eq!(
            second.edit_script.as_ref().unwrap().ops_vec(),
            vec![(GapAlignOpType::Sub, 15)]
        );
        assert_eq!(
            second.map_info.as_ref().unwrap().left_edge & MAPPER_SPLICE_SIGNAL,
            0
        );
    }

    #[test]
    fn find_splice_junctions_for_overlaps_clears_splice_bits_without_trimming_first_hsp() {
        let mut first = hsp(10, 0, 0, 10, 20);
        first.query_end = 20;
        first.subject_end = 30;
        first.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![(
            GapAlignOpType::Sub,
            20,
        )]));
        first.map_info = Some(crate::hspstream::BlastHSPMappingInfo {
            edits: None,
            subject_overhangs: None,
            left_edge: 0,
            right_edge: MAPPER_SPLICE_SIGNAL | 1,

            flags: 0,
        });

        let mut second = hsp(30, 0, 15, 40, 15);
        second.query_end = 30;
        second.subject_end = 55;
        second.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![(
            GapAlignOpType::Sub,
            15,
        )]));

        assert_eq!(
            s_find_splice_junctions_for_overlaps(&mut first, &mut second, None, 30, true),
            0
        );
        assert_eq!(first.query_end, 20);
        assert_eq!(first.subject_end, 30);
        assert_eq!(second.query_offset, 15);
        assert_eq!(
            first.edit_script.as_ref().unwrap().ops_vec(),
            vec![(GapAlignOpType::Sub, 20)]
        );
        assert_eq!(
            first.map_info.as_ref().unwrap().right_edge & MAPPER_SPLICE_SIGNAL,
            0
        );
    }

    #[test]
    fn find_splice_junctions_for_overlaps_rejects_extreme_malformed_overlap() {
        let mut first = hsp(10, 0, 0, 0, 20);
        first.query_offset = i32::MIN;
        first.query_end = i32::MAX;
        first.subject_end = 20;
        let mut second = hsp(20, 0, -1, 40, 15);
        second.query_end = i32::MAX;
        second.subject_end = 55;

        assert_eq!(
            s_find_splice_junctions_for_overlaps(&mut first, &mut second, None, 30, false),
            -1
        );
    }

    #[test]
    fn find_splice_junctions_for_overlaps_trims_first_when_it_has_more_overlap_edits() {
        let query = vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 3, 1, 2];
        let mut first = hsp(100, 0, 0, 0, 10);
        first.query_end = 10;
        first.subject_end = 10;
        first.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![(
            GapAlignOpType::Sub,
            10,
        )]));
        first.map_info = Some(crate::hspstream::BlastHSPMappingInfo {
            edits: Some(crate::gapinfo::JumperEditsBlock {
                num_edits: 0,
                edits: vec![
                    JumperEdit {
                        query_pos: 4,
                        query_base: 0,
                        subject_base: 1,
                    },
                    JumperEdit {
                        query_pos: 8,
                        query_base: 0,
                        subject_base: 2,
                    },
                    JumperEdit {
                        query_pos: 9,
                        query_base: 3,
                        subject_base: 1,
                    },
                ],
            }),
            subject_overhangs: None,
            left_edge: 0,
            right_edge: MAPPER_SPLICE_SIGNAL | 0x03,

            flags: 0,
        });

        let mut second = hsp(100, 0, 8, 20, 10);
        second.query_end = 18;
        second.subject_end = 30;
        second.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![(
            GapAlignOpType::Sub,
            10,
        )]));
        second.map_info = Some(crate::hspstream::BlastHSPMappingInfo {
            edits: Some(crate::gapinfo::JumperEditsBlock::default()),
            subject_overhangs: None,
            left_edge: MAPPER_SPLICE_SIGNAL | 0x02,
            right_edge: 0,

            flags: 0,
        });

        assert_eq!(
            s_find_splice_junctions_for_overlaps(&mut first, &mut second, Some(&query), 12, true),
            0
        );
        assert_eq!(first.query_end, 8);
        assert_eq!(first.subject_end, 8);
        assert_eq!(second.query_offset, 8);
        assert_eq!(
            first.edit_script.as_ref().unwrap().ops_vec(),
            vec![(GapAlignOpType::Sub, 8)]
        );
        assert_eq!(
            first
                .map_info
                .as_ref()
                .unwrap()
                .edits
                .as_ref()
                .unwrap()
                .edits,
            vec![JumperEdit {
                query_pos: 4,
                query_base: 0,
                subject_base: 1,
            }]
        );
        assert_eq!(
            first.map_info.as_ref().unwrap().right_edge & MAPPER_SPLICE_SIGNAL,
            0
        );
    }

    #[test]
    fn find_splice_junctions_for_overlaps_boundary_edit_does_not_block_signal_split() {
        let query = vec![0, 1, 2, 3, 0, 1, 2, 3, 2, 3, 0, 3, 1, 2];
        let mut first = hsp(100, 0, 0, 0, 12);
        first.query_end = 12;
        first.subject_end = 12;
        first.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![(
            GapAlignOpType::Sub,
            12,
        )]));
        first.map_info = Some(crate::hspstream::BlastHSPMappingInfo {
            edits: Some(crate::gapinfo::JumperEditsBlock {
                num_edits: 0,
                edits: vec![
                    JumperEdit {
                        query_pos: 10,
                        query_base: 0,
                        subject_base: 2,
                    },
                    JumperEdit {
                        query_pos: 11,
                        query_base: 3,
                        subject_base: 1,
                    },
                ],
            }),
            subject_overhangs: None,
            left_edge: 0,
            right_edge: 0,

            flags: 0,
        });

        let mut second = hsp(100, 0, 8, 20, 10);
        second.query_end = 18;
        second.subject_end = 30;
        second.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![(
            GapAlignOpType::Sub,
            10,
        )]));
        second.map_info = Some(crate::hspstream::BlastHSPMappingInfo {
            edits: Some(crate::gapinfo::JumperEditsBlock {
                num_edits: 0,
                edits: vec![JumperEdit {
                    query_pos: 10,
                    query_base: 0,
                    subject_base: 1,
                }],
            }),
            subject_overhangs: None,
            left_edge: 0x02,
            right_edge: 0,

            flags: 0,
        });

        assert_eq!(
            s_find_splice_junctions_for_overlaps(&mut first, &mut second, Some(&query), 14, true),
            0
        );
        assert_eq!(first.query_end, 8);
        assert_eq!(second.query_offset, 8);
        assert_ne!(
            first.map_info.as_ref().unwrap().right_edge & MAPPER_SPLICE_SIGNAL,
            0
        );
        assert_ne!(
            second.map_info.as_ref().unwrap().left_edge & MAPPER_SPLICE_SIGNAL,
            0
        );
    }

    #[test]
    fn find_splice_junctions_for_overlaps_trims_second_when_it_has_more_overlap_edits() {
        let query = vec![0, 1, 2, 3, 0, 1, 2, 3, 1, 2, 0, 3, 1, 2];
        let mut first = hsp(100, 0, 0, 0, 10);
        first.query_end = 10;
        first.subject_end = 10;
        first.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![(
            GapAlignOpType::Sub,
            10,
        )]));
        first.map_info = Some(crate::hspstream::BlastHSPMappingInfo {
            edits: Some(crate::gapinfo::JumperEditsBlock::default()),
            subject_overhangs: None,
            left_edge: 0,
            right_edge: MAPPER_SPLICE_SIGNAL | 0x03,

            flags: 0,
        });

        let mut second = hsp(100, 0, 8, 20, 10);
        second.query_end = 18;
        second.subject_end = 30;
        second.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![(
            GapAlignOpType::Sub,
            10,
        )]));
        second.map_info = Some(crate::hspstream::BlastHSPMappingInfo {
            edits: Some(crate::gapinfo::JumperEditsBlock {
                num_edits: 0,
                edits: vec![
                    JumperEdit {
                        query_pos: 8,
                        query_base: 1,
                        subject_base: 2,
                    },
                    JumperEdit {
                        query_pos: 9,
                        query_base: 2,
                        subject_base: 3,
                    },
                    JumperEdit {
                        query_pos: 12,
                        query_base: 1,
                        subject_base: 0,
                    },
                ],
            }),
            subject_overhangs: None,
            left_edge: MAPPER_SPLICE_SIGNAL | 0x02,
            right_edge: 0,

            flags: 0,
        });

        assert_eq!(
            s_find_splice_junctions_for_overlaps(&mut first, &mut second, Some(&query), 14, true),
            0
        );
        assert_eq!(first.query_end, 10);
        assert_eq!(second.query_offset, 10);
        assert_eq!(second.subject_offset, 22);
        assert_eq!(
            second.edit_script.as_ref().unwrap().ops_vec(),
            vec![(GapAlignOpType::Sub, 8)]
        );
        assert_eq!(
            second
                .map_info
                .as_ref()
                .unwrap()
                .edits
                .as_ref()
                .unwrap()
                .edits,
            vec![JumperEdit {
                query_pos: 12,
                query_base: 1,
                subject_base: 0,
            }]
        );
        assert_eq!(
            second.map_info.as_ref().unwrap().left_edge & MAPPER_SPLICE_SIGNAL,
            0
        );
    }

    #[test]
    fn find_splice_junctions_for_overlaps_equal_overlap_edits_clear_without_trimming() {
        let query = vec![0, 1, 2, 3, 0, 1, 2, 3, 1, 2, 0, 3];
        let mut first = hsp(100, 0, 0, 0, 10);
        first.query_end = 10;
        first.subject_end = 10;
        first.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![(
            GapAlignOpType::Sub,
            10,
        )]));
        first.map_info = Some(crate::hspstream::BlastHSPMappingInfo {
            edits: Some(crate::gapinfo::JumperEditsBlock {
                num_edits: 0,
                edits: vec![JumperEdit {
                    query_pos: 9,
                    query_base: 2,
                    subject_base: 3,
                }],
            }),
            subject_overhangs: None,
            left_edge: 0,
            right_edge: MAPPER_SPLICE_SIGNAL | 0x03,

            flags: 0,
        });

        let mut second = hsp(100, 0, 8, 20, 10);
        second.query_end = 18;
        second.subject_end = 30;
        second.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![(
            GapAlignOpType::Sub,
            10,
        )]));
        second.map_info = Some(crate::hspstream::BlastHSPMappingInfo {
            edits: Some(crate::gapinfo::JumperEditsBlock {
                num_edits: 0,
                edits: vec![JumperEdit {
                    query_pos: 8,
                    query_base: 1,
                    subject_base: 2,
                }],
            }),
            subject_overhangs: None,
            left_edge: MAPPER_SPLICE_SIGNAL | 0x02,
            right_edge: 0,

            flags: 0,
        });

        assert_eq!(
            s_find_splice_junctions_for_overlaps(&mut first, &mut second, Some(&query), 12, true),
            0
        );
        assert_eq!(first.query_end, 10);
        assert_eq!(second.query_offset, 8);
        assert_eq!(
            first.edit_script.as_ref().unwrap().ops_vec(),
            vec![(GapAlignOpType::Sub, 10)]
        );
        assert_eq!(
            second.edit_script.as_ref().unwrap().ops_vec(),
            vec![(GapAlignOpType::Sub, 10)]
        );
        assert_eq!(
            first.map_info.as_ref().unwrap().right_edge & MAPPER_SPLICE_SIGNAL,
            0
        );
        assert_eq!(
            second.map_info.as_ref().unwrap().left_edge & MAPPER_SPLICE_SIGNAL,
            0
        );
    }

    #[test]
    fn find_splice_junctions_for_overlaps_uses_signal_split_when_query_flanks_available() {
        let mut query = vec![1; 40];
        query[6..10].copy_from_slice(&[0, 2, 2, 3]);

        let mut first = hsp(100, 0, 0, 0, 20);
        first.query_end = 10;
        first.subject_end = 10;
        first.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![(
            GapAlignOpType::Sub,
            10,
        )]));
        first.map_info = Some(crate::hspstream::BlastHSPMappingInfo {
            edits: None,
            subject_overhangs: None,
            left_edge: 0,
            right_edge: 0x0b,

            flags: 0,
        });

        let mut second = hsp(100, 0, 8, 8, 20);
        second.query_offset = 8;
        second.query_end = 20;
        second.subject_offset = 8;
        second.subject_end = 20;
        second.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![(
            GapAlignOpType::Sub,
            10,
        )]));
        second.map_info = Some(crate::hspstream::BlastHSPMappingInfo {
            edits: None,
            subject_overhangs: None,
            left_edge: 0x02,
            right_edge: 0,

            flags: 0,
        });

        assert_eq!(
            s_find_splice_junctions_for_overlaps(&mut first, &mut second, Some(&query), 40, true),
            0
        );
        assert_eq!(first.query_end, 8);
        assert_eq!(first.subject_end, 8);
        assert_eq!(second.query_offset, 8);
        assert_eq!(second.subject_offset, 8);
        assert_ne!(
            first.map_info.as_ref().unwrap().right_edge & MAPPER_SPLICE_SIGNAL,
            0
        );
        assert_ne!(
            second.map_info.as_ref().unwrap().left_edge & MAPPER_SPLICE_SIGNAL,
            0
        );
        assert_eq!(
            first.edit_script.as_ref().unwrap().ops_vec(),
            vec![(GapAlignOpType::Sub, 8)]
        );
    }

    #[test]
    fn extend_alignment_cleanup_drops_owned_rust_payloads() {
        let script = crate::gapinfo::GapEditScript::from_ops(vec![(GapAlignOpType::Sub, 4)]);
        s_extend_alignment_cleanup(Some(vec![0, 1, 2, 3]), Some(()), Some(script), Some(()));
        s_extend_alignment_cleanup::<(), ()>(None, None, None, None);
    }

    #[test]
    fn find_splice_junctions_for_gap_extends_second_left_from_canonical_overhangs() {
        let scores = ScoringOptions::new_blastn();
        let query = vec![0; 40];

        let mut first = hsp(10, 0, 0, 0, 10);
        first.query_end = 10;
        first.subject_end = 10;
        first.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![(
            GapAlignOpType::Sub,
            10,
        )]));

        let mut second = hsp(10, 0, 14, 100, 10);
        second.query_end = 24;
        second.subject_end = 110;
        second.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![(
            GapAlignOpType::Sub,
            10,
        )]));

        let first_right = [2, 3, 0, 0, 0, 0];
        let second_left = [0, 2, 0, 0, 0, 0];
        assert_eq!(
            s_find_splice_junctions_for_gap(
                &mut first,
                &mut second,
                Some(&query),
                query.len() as i32,
                Some(&scores),
                Some(&first_right),
                Some(&second_left),
            ),
            0
        );

        assert_eq!(first.query_end, 10);
        assert_eq!(second.query_offset, 10);
        assert_eq!(second.subject_offset, 96);
        assert_ne!(
            first.map_info.as_ref().unwrap().right_edge & MAPPER_SPLICE_SIGNAL,
            0
        );
        assert_ne!(
            second.map_info.as_ref().unwrap().left_edge & MAPPER_SPLICE_SIGNAL,
            0
        );
        assert_eq!(
            second.edit_script.as_ref().unwrap().ops_vec(),
            vec![(GapAlignOpType::Sub, 14)]
        );
    }

    #[test]
    fn find_splice_junctions_for_gap_extends_first_right_from_canonical_overhangs() {
        let scores = ScoringOptions::new_blastn();
        let query = vec![0; 40];

        let mut first = hsp(10, 0, 0, 0, 10);
        first.query_end = 10;
        first.subject_end = 10;
        first.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![(
            GapAlignOpType::Sub,
            10,
        )]));

        let mut second = hsp(10, 0, 14, 100, 10);
        second.query_end = 24;
        second.subject_end = 110;
        second.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![(
            GapAlignOpType::Sub,
            10,
        )]));

        let first_right = [0, 0, 0, 0, 2, 3];
        let second_left = [0, 0, 0, 0, 0, 2];
        assert_eq!(
            s_find_splice_junctions_for_gap(
                &mut first,
                &mut second,
                Some(&query),
                query.len() as i32,
                Some(&scores),
                Some(&first_right),
                Some(&second_left),
            ),
            0
        );

        assert_eq!(first.query_end, 14);
        assert_eq!(first.subject_end, 14);
        assert_eq!(second.query_offset, 14);
        assert_ne!(
            first.map_info.as_ref().unwrap().right_edge & MAPPER_SPLICE_SIGNAL,
            0
        );
        assert_ne!(
            second.map_info.as_ref().unwrap().left_edge & MAPPER_SPLICE_SIGNAL,
            0
        );
        assert_eq!(
            first.edit_script.as_ref().unwrap().ops_vec(),
            vec![(GapAlignOpType::Sub, 14)]
        );
    }

    #[test]
    fn find_splice_junctions_for_gap_rejects_missing_state_and_clears_stale_edges() {
        let scores = ScoringOptions::new_blastn();
        let query = vec![1; 40];
        let mut first = hsp(10, 0, 0, 0, 10);
        first.query_end = 10;
        first.subject_end = 10;
        first.map_info = Some(crate::hspstream::BlastHSPMappingInfo {
            edits: None,
            subject_overhangs: None,
            left_edge: 0,
            right_edge: MAPPER_SPLICE_SIGNAL | MAPPER_EXON | 7,

            flags: 0,
        });
        let mut second = hsp(10, 0, 14, 100, 10);
        second.query_end = 24;
        second.subject_end = 110;
        second.map_info = Some(crate::hspstream::BlastHSPMappingInfo {
            edits: None,
            subject_overhangs: None,
            left_edge: MAPPER_SPLICE_SIGNAL | MAPPER_EXON | 3,
            right_edge: 0,

            flags: 0,
        });

        assert_eq!(
            s_find_splice_junctions_for_gap(
                &mut first,
                &mut second,
                None,
                query.len() as i32,
                Some(&scores),
                Some(&[1, 1, 1, 1, 1, 1]),
                Some(&[1, 1, 1, 1, 1, 1]),
            ),
            -1
        );
        assert_ne!(
            first.map_info.as_ref().unwrap().right_edge & MAPPER_SPLICE_SIGNAL,
            0
        );

        assert_eq!(
            s_find_splice_junctions_for_gap(
                &mut first,
                &mut second,
                Some(&query),
                query.len() as i32,
                None,
                Some(&[1, 1, 1, 1, 1, 1]),
                Some(&[1, 1, 1, 1, 1, 1]),
            ),
            -1
        );
        assert_ne!(
            second.map_info.as_ref().unwrap().left_edge & MAPPER_SPLICE_SIGNAL,
            0
        );

        assert_eq!(
            s_find_splice_junctions_for_gap(
                &mut first,
                &mut second,
                Some(&query),
                query.len() as i32,
                Some(&scores),
                Some(&[1, 1, 1, 1, 1, 1]),
                Some(&[1, 1, 1, 1, 1, 1]),
            ),
            0
        );
        assert_eq!(
            first.map_info.as_ref().unwrap().right_edge & MAPPER_SPLICE_SIGNAL,
            0
        );
        assert_eq!(
            second.map_info.as_ref().unwrap().left_edge & MAPPER_SPLICE_SIGNAL,
            0
        );

        let mut extreme_first = hsp(10, 0, 0, 0, 10);
        extreme_first.query_end = i32::MAX;
        extreme_first.subject_end = i32::MAX;
        let mut extreme_second = hsp(10, 0, 0, 0, 10);
        extreme_second.query_offset = i32::MIN;
        extreme_second.query_end = 0;
        extreme_second.subject_offset = i32::MIN;
        extreme_second.subject_end = 0;
        assert_eq!(
            s_find_splice_junctions_for_gap(
                &mut extreme_first,
                &mut extreme_second,
                Some(&query),
                query.len() as i32,
                Some(&scores),
                Some(&[1, 1, 1, 1, 1, 1]),
                Some(&[1, 1, 1, 1, 1, 1]),
            ),
            0
        );
    }

    #[test]
    fn find_splice_junctions_for_gap_reads_overhangs_from_hsp_map_info() {
        let scores = ScoringOptions::new_blastn();
        let query = vec![0; 40];

        let mut first = hsp(10, 0, 0, 0, 10);
        first.query_end = 10;
        first.subject_end = 10;
        first.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![(
            GapAlignOpType::Sub,
            10,
        )]));
        first.map_info = Some(crate::hspstream::BlastHSPMappingInfo {
            edits: Some(crate::gapinfo::JumperEditsBlock::default()),
            subject_overhangs: Some(crate::gapinfo::SequenceOverhangs {
                left_len: 0,
                right_len: 0,
                left: None,
                right: Some(vec![2, 3, 0, 0, 0, 0]),
            }),
            left_edge: 0,
            right_edge: 0,

            flags: 0,
        });

        let mut second = hsp(10, 0, 14, 100, 10);
        second.query_end = 24;
        second.subject_end = 110;
        second.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![(
            GapAlignOpType::Sub,
            10,
        )]));
        second.map_info = Some(crate::hspstream::BlastHSPMappingInfo {
            edits: Some(crate::gapinfo::JumperEditsBlock::default()),
            subject_overhangs: Some(crate::gapinfo::SequenceOverhangs {
                left_len: 0,
                right_len: 0,
                left: Some(vec![0, 2, 0, 0, 0, 0]),
                right: None,
            }),
            left_edge: 0,
            right_edge: 0,

            flags: 0,
        });

        assert_eq!(
            s_find_splice_junctions_for_gap_using_map_info(
                &mut first,
                &mut second,
                Some(&query),
                query.len() as i32,
                Some(&scores),
            ),
            0
        );

        assert_eq!(second.query_offset, 10);
        assert_eq!(second.subject_offset, 96);
        assert_ne!(
            first.map_info.as_ref().unwrap().right_edge & MAPPER_SPLICE_SIGNAL,
            0
        );
        assert_ne!(
            second.map_info.as_ref().unwrap().left_edge & MAPPER_SPLICE_SIGNAL,
            0
        );
    }

    #[test]
    fn find_adapter_in_sequence_allows_one_mismatch_near_hsp_end() {
        let mut query = vec![3; 30];
        query[10..22].copy_from_slice(&[0, 2, 0, 3, 1, 2, 2, 0, 0, 2, 0, 2]);
        assert_eq!(s_find_adapter_in_sequence(0, 22, Some(&query), 30), 10);
        query[15] = 0;
        assert_eq!(s_find_adapter_in_sequence(0, 22, Some(&query), 30), 10);
        query[10] = 1;
        assert_eq!(s_find_adapter_in_sequence(0, 22, Some(&query), 30), -1);
        assert_eq!(s_find_adapter_in_sequence(12, 12, Some(&query), 30), -1);
        assert_eq!(s_find_adapter_in_sequence(0, 30, None, 30), -1);
    }

    #[test]
    fn find_poly_a_in_sequence_matches_tail_rules() {
        assert_eq!(s_find_poly_a_in_sequence(Some(&[1, 2, 0, 0, 0, 0]), 6), 2);
        assert_eq!(
            s_find_poly_a_in_sequence(Some(&[1, 0, 0, 2, 0, 0, 0, 0, 0]), 9),
            1
        );
        assert_eq!(s_find_poly_a_in_sequence(Some(&[1, 2, 0, 0]), 4), -1);
        assert_eq!(s_find_poly_a_in_sequence(Some(&[0, 0, 0, 0]), -3), -1);
        assert_eq!(s_find_poly_a_in_sequence(None, 4), -1);
    }

    #[test]
    fn set_poly_a_tail_marks_eligible_chain_ends() {
        let mut list = Some(chain(30, 1, 1, 0, 0));
        assert_eq!(s_set_poly_a_tail(&mut list, 12, -1, 20), 0);
        assert_eq!(list.as_ref().unwrap().poly_a, 12);
        let edge = last_container(list.as_ref().unwrap())
            .unwrap()
            .hsp
            .map_info
            .as_ref()
            .unwrap()
            .right_edge;
        assert_eq!(edge & MAPPER_POLY_A, MAPPER_POLY_A);
        assert_eq!(edge & MAPPER_EXON, MAPPER_EXON);

        let mut too_close = Some(chain(30, 1, 1, 6, 0));
        assert_eq!(s_set_poly_a_tail(&mut too_close, 16, -1, 20), 0);
        assert_eq!(too_close.as_ref().unwrap().poly_a, 0);
        assert!(last_container(too_close.as_ref().unwrap())
            .unwrap()
            .hsp
            .map_info
            .is_none());

        let mut no_hsps = Some(hsp_chain_new(1));
        assert_eq!(s_set_poly_a_tail(&mut no_hsps, 12, -1, 20), 0);
        assert_eq!(no_hsps.as_ref().unwrap().poly_a, 0);

        let mut malformed_len = Some(chain(30, 1, 1, 0, 0));
        assert_eq!(s_set_poly_a_tail(&mut malformed_len, 12, -1, i32::MIN), 0);
        assert_eq!(malformed_len.as_ref().unwrap().poly_a, 0);
        assert!(last_container(malformed_len.as_ref().unwrap())
            .unwrap()
            .hsp
            .map_info
            .is_none());
    }

    #[test]
    fn find_poly_a_tails_scans_query_contexts() {
        let mut saved = vec![Some(chain(30, 1, 1, 2, 0))];
        let mut query = vec![3; 41];
        query[14..20].copy_from_slice(&[0, 0, 0, 0, 0, 0]);
        let query_info = QueryInfo::new_blastn(&[20]);
        assert_eq!(s_find_poly_a_tails(&mut saved, &query, &query_info), 0);
        assert_eq!(saved[0].as_ref().unwrap().poly_a, 14);
    }

    #[test]
    fn find_poly_a_tails_skips_adapter_full_coverage_and_missing_query_slices() {
        let query_info = QueryInfo::new_blastn(&[20]);

        let mut adapter_marked = Some(chain(30, 1, 1, 2, 0));
        adapter_marked.as_mut().unwrap().adapter = 12;
        let mut saved = vec![adapter_marked];
        let query = vec![0; 42];
        assert_eq!(s_find_poly_a_tails(&mut saved, &query, &query_info), 0);
        assert_eq!(saved[0].as_ref().unwrap().poly_a, 0);

        let mut full_coverage = Some(chain(30, 1, 2, 0, 0));
        full_coverage
            .as_mut()
            .unwrap()
            .hsps
            .as_mut()
            .unwrap()
            .hsp
            .query_end = 19;
        let mut saved = vec![full_coverage];
        assert_eq!(s_find_poly_a_tails(&mut saved, &query, &query_info), 0);
        assert_eq!(saved[0].as_ref().unwrap().poly_a, 0);

        let mut short_query_saved = vec![Some(chain(30, 1, 3, 2, 0))];
        assert_eq!(
            s_find_poly_a_tails(&mut short_query_saved, &[0; 10], &query_info),
            0
        );
        assert_eq!(short_query_saved[0].as_ref().unwrap().poly_a, 0);
    }

    #[test]
    fn sort_and_filter_chains_operate_on_linked_lists() {
        let mut list = Some(chain(10, 1, 1, 0, 0));
        list.as_mut().unwrap().next = Some(chain(50, 1, 2, 20, 20));
        list.as_mut().unwrap().next.as_mut().unwrap().next = Some(chain(30, 1, 3, 40, 40));

        let mut saved = vec![None, list];
        assert_eq!(s_sort_chains(&mut saved), 0);
        assert!(saved[0].is_none());
        assert_eq!(scores(&saved[1]), vec![50, 30, 10]);

        assert_eq!(s_filter_chains(&mut saved[1], 25, -1), 0);
        assert_eq!(scores(&saved[1]), vec![50, 30]);
    }

    #[test]
    fn filter_chains_applies_edit_distance_and_handles_empty_lists() {
        let mut empty = None;
        assert_eq!(s_filter_chains(&mut empty, 10, 0), 0);
        assert!(empty.is_none());

        let perfect = chain(30, 1, 1, 0, 0);
        let mut one_edit = chain(30, 1, 2, 20, 20);
        one_edit.hsps.as_mut().unwrap().hsp.num_ident = 9;
        let mut two_edits = chain(30, 1, 3, 40, 40);
        two_edits.hsps.as_mut().unwrap().hsp.num_ident = 8;

        let mut list = Some(two_edits);
        list.as_mut().unwrap().next = Some(one_edit);
        list.as_mut().unwrap().next.as_mut().unwrap().next = Some(perfect);

        assert_eq!(s_filter_chains(&mut list, 25, 1), 0);
        assert_eq!(scores(&list), vec![30, 30]);
        assert_eq!(list.as_ref().unwrap().oid, 2);
        assert_eq!(list.as_ref().unwrap().next.as_ref().unwrap().oid, 1);
    }

    #[test]
    fn hsp_node_array_copy_retargets_path_indexes() {
        let source = vec![
            HSPNode {
                hsp_index: Some(0),
                best_score: 50,
                path_next: Some(1),
            },
            HSPNode {
                hsp_index: Some(1),
                best_score: 20,
                path_next: None,
            },
        ];
        let mut dest = vec![HSPNode::default(); 2];
        assert_eq!(s_hsp_node_array_copy(&mut dest, &source, 2), 0);
        assert_eq!(dest, source);
        assert_eq!(s_hsp_node_array_copy(&mut dest, &source, 0), 0);
        assert_eq!(s_hsp_node_array_copy(&mut dest, &source, -1), -1);
        assert_eq!(s_hsp_node_array_copy(&mut dest[..1], &source, 2), -1);
        assert_eq!(s_hsp_node_array_copy(&mut dest, &source[..1], 2), -1);

        let invalid = vec![HSPNode {
            hsp_index: Some(0),
            best_score: 1,
            path_next: Some(3),
        }];
        assert_eq!(s_hsp_node_array_copy(&mut dest[..1], &invalid, 1), -1);
    }

    #[test]
    fn hsp_path_new_and_free_match_c_lifecycle() {
        let mut path = hsp_path_new();
        assert_eq!(path.start.len(), MAX_NUM_HSP_PATHS);
        assert_eq!(path.num_paths, 0);
        assert_eq!(path.score, 0);
        path.start[0] = Some(7);
        assert!(hsp_path_free(Some(path)).is_none());
        assert!(hsp_path_free(None).is_none());
    }
}
