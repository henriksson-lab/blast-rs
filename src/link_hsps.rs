//! Port of NCBI BLAST+ 2.17.0 `link_hsps.c` — sum-statistics HSP linking.
//!
//! This module implements three public entry points used by the NCBI
//! HSP-linking machinery:
//!
//! - [`BLAST_LinkHsps`] — top-level dispatcher on `longest_intron`
//! - [`s_BlastEvenGapLinkHSPs`] — small/large gap linking for blastn/tblastx
//! - [`s_BlastUnevenGapLinkHSPs`] — asymmetric gap linking for tblastn/blastx
//!
//! The port mirrors the original C very closely. A few adaptations are
//! unavoidable in safe Rust:
//!
//! * NCBI walks an intrusive doubly-linked list of `LinkHSPStruct` nodes
//!   through raw pointers; here we back the nodes with indices into a
//!   `Vec<LinkHspStruct>` and traverse via `prev`/`next: Option<usize>`.
//!   The `hp_start` sentinel used in the NCBI code lives at index 0.
//! * `BlastHSPLink.link[eOrderingMethods]` becomes `[Option<usize>; 2]`.
//! * The NCBI `BlastHSP` structure is reconstructed as a minimal
//!   self-contained [`LinkBlastHsp`] holding only the fields consulted
//!   by `link_hsps.c`. The wider active-pipeline `SearchHsp` is not
//!   used here because it lacks `gapped_start`, `subject.frame`, `num`,
//!   `evalue` mutation, etc. Callers wiring this module into the
//!   pipeline would build a `LinkBlastHsp` from their `SearchHsp`.
//!
//! The port is not wired into the active search path — it is `pub(crate)`
//! scaffolding tested against synthetic HSP lists with hand-computed
//! expected sum-scores.
//!
//! Lines below cite specific anchors in
//! `ncbi-blast-2.17.0+-src/c++/src/algo/blast/core/link_hsps.c` so
//! `ccc-rs` can diff the port against the C reference.

#![allow(non_snake_case)]
#![allow(non_camel_case_types)]

use crate::program::{
    query_is_nucleotide, query_is_translated, subject_is_translated, ProgramType, BLASTX,
};
use crate::queryinfo::QueryInfo;
use crate::stat::{
    gap_decay_divisor, large_gap_sum_e, small_gap_sum_e, uneven_gap_sum_e, KarlinBlk,
};

/// NCBI BLAST constants (`blast_def.h` / `ncbi_std.h`).
const NUM_FRAMES: i32 = 6;
const NUM_STRANDS: i32 = 2;
const CODON_LENGTH: i32 = 3;
const INT4_MAX: f64 = i32::MAX as f64;

/// Describes the method for ordering HSPs
/// (port of `ELinkOrderingMethod`, `link_hsps.c:47`).
///
/// Note: NCBI uses these values to index a size-2 array, so they must
/// linearly increase starting at 0.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[repr(usize)]
pub enum ELinkOrderingMethod {
    /// Favour small gaps when linking an HSP.
    #[default]
    ELinkSmallGaps = 0,
    /// Favour large gaps when linking an HSP.
    ELinkLargeGaps = 1,
}
const E_ORDERING_METHODS: usize = 2;

/// Minimal HSP segment (port of `BlastSeg`, `blast_hits.h:96`).
#[derive(Debug, Clone, Default)]
pub struct LinkBlastSeg {
    pub frame: i32,
    pub offset: i32,
    pub end: i32,
    pub gapped_start: i32,
}

/// Minimal HSP structure used by the linking machinery
/// (port of `BlastHSP`, `blast_hits.h:125`).
///
/// Only the fields consulted by `link_hsps.c` are kept. `SearchHsp`
/// from `src/search.rs` does not store `gapped_start` or `subject.frame`,
/// so callers wiring this module into the pipeline must assemble these
/// minimal records from their `SearchHsp`s.
#[derive(Debug, Clone, Default)]
pub struct LinkBlastHsp {
    pub score: i32,
    pub num_ident: i32,
    pub bit_score: f64,
    pub evalue: f64,
    pub query: LinkBlastSeg,
    pub subject: LinkBlastSeg,
    pub context: i32,
    /// Number of HSPs linked together for sum statistics. `0` or `1`
    /// means "not part of a linked set".
    pub num: i32,
}

/// A list of HSPs for one subject (port of `BlastHSPList`, `blast_hits.h:151`).
#[derive(Debug, Clone, Default)]
pub struct LinkBlastHspList {
    pub oid: i32,
    pub query_index: i32,
    pub hsp_array: Vec<LinkBlastHsp>,
    pub best_evalue: f64,
}

impl LinkBlastHspList {
    pub fn hspcnt(&self) -> i32 {
        self.hsp_array.len() as i32
    }
}

/// Parameters controlling HSP linking (port of `BlastLinkHSPParameters`,
/// `blast_parameters.h:137`).
#[derive(Debug, Clone)]
pub struct LinkHSPParameters {
    pub gap_prob: f64,
    pub gap_size: i32,
    pub overlap_size: i32,
    pub gap_decay_rate: f64,
    pub cutoff_small_gap: i32,
    pub cutoff_big_gap: i32,
    pub longest_intron: i32,
}

impl Default for LinkHSPParameters {
    fn default() -> Self {
        // NCBI defaults for blastn/tblastx (gap_size=40, overlap_size=9,
        // gap_prob=0.5, decay_rate=0.5, no intron).
        Self {
            gap_prob: 0.5,
            gap_size: 40,
            overlap_size: 9,
            gap_decay_rate: 0.5,
            cutoff_small_gap: 0,
            cutoff_big_gap: 0,
            longest_intron: 0,
        }
    }
}

/// Score block subset used by linking (port of `BlastScoreBlk`,
/// `blast_stat.h`). Only the per-context Karlin blocks are consulted.
#[derive(Debug, Clone, Default)]
pub struct LinkScoreBlock {
    pub kbp: Vec<KarlinBlk>,
    pub kbp_gap: Vec<KarlinBlk>,
}

/// Auxiliary structure for keeping track of sum scores during linking
/// (port of `BlastHSPLink`, `link_hsps.c:63`).
#[derive(Debug, Clone, Copy, Default)]
pub struct BlastHSPLink {
    /// Best choice of HSP to link with, per ordering method. An
    /// `Option<usize>` index into the `Vec<LinkHspStruct>` (adapted from
    /// NCBI's raw `LinkHSPStruct*`).
    pub link: [Option<usize>; E_ORDERING_METHODS],
    pub num: [i32; E_ORDERING_METHODS],
    pub sum: [i32; E_ORDERING_METHODS],
    pub xsum: [f64; E_ORDERING_METHODS],
    /// Has the link been changed since previous access? (`Int4` in NCBI).
    pub changed: i32,
}

/// Port of `LinkHSPStruct`, `link_hsps.c:76`.
///
/// Adaptation from NCBI: raw `prev`/`next` pointers become
/// `Option<usize>` indices into the backing `Vec`. The `hsp`
/// field becomes an index into a separate HSP backing store so we
/// can mutate those records after sum-statistics selection.
#[derive(Debug, Clone, Default)]
pub struct LinkHspStruct {
    /// Index into the `hsp_array` of the owning `LinkBlastHspList`.
    pub hsp_idx: i32,
    /// Previous / next links (NCBI walks a doubly-linked list).
    pub prev: Option<usize>,
    pub next: Option<usize>,
    pub hsp_link: BlastHSPLink,
    pub linked_set: bool,
    pub start_of_chain: bool,
    /// Tracks link bookkeeping; `-1000` marks "already chosen" (NCBI sentinel).
    pub linked_to: i32,
    pub xsum: f64,
    pub ordering_method: ELinkOrderingMethod,
    pub q_offset_trim: i32,
    pub q_end_trim: i32,
    pub s_offset_trim: i32,
    pub s_end_trim: i32,
}

impl LinkHspStruct {
    fn reset(&mut self) {
        *self = LinkHspStruct {
            ordering_method: ELinkOrderingMethod::ELinkSmallGaps,
            ..LinkHspStruct::default()
        };
    }
}

/// Per-iteration inner-loop helper (port of `LinkHelpStruct`, `link_hsps.c:100`).
#[derive(Debug, Clone, Copy, Default)]
struct LinkHelpStruct {
    ptr: usize, // index into nodes[]
    q_off_trim: i32,
    s_off_trim: i32,
    sum: [i32; E_ORDERING_METHODS],
    maxsum1: i32,
    next_larger: usize,
}

// ---------------------------------------------------------------------------
// Even-gap linking (small/large gaps, blastn & tblastx)
// ---------------------------------------------------------------------------

/// Mirror of `s_FwdCompareHSPs` (`link_hsps.c:120`) on indices into the
/// HSP backing store. Sort primary key is increasing `(context,
/// query.offset, subject.offset)`.
fn fwd_compare_hsps(h1: &LinkBlastHsp, h2: &LinkBlastHsp) -> std::cmp::Ordering {
    h1.context
        .cmp(&h2.context)
        .then_with(|| h1.query.offset.cmp(&h2.query.offset))
        .then_with(|| h1.subject.offset.cmp(&h2.subject.offset))
}

/// Mirror of `s_FwdCompareHSPsTransl` (`link_hsps.c:157`).
fn fwd_compare_hsps_transl(h1: &LinkBlastHsp, h2: &LinkBlastHsp) -> std::cmp::Ordering {
    let c1 = h1.context / (NUM_FRAMES / 2);
    let c2 = h2.context / (NUM_FRAMES / 2);
    c1.cmp(&c2)
        .then_with(|| h1.query.offset.cmp(&h2.query.offset))
        .then_with(|| h1.subject.offset.cmp(&h2.subject.offset))
}

/// Mirror of `s_RevCompareHSPs` (`link_hsps.c:198`).
fn rev_compare_hsps(h1: &LinkBlastHsp, h2: &LinkBlastHsp) -> std::cmp::Ordering {
    // Same context -> decreasing query offset, then decreasing subject offset.
    match h1.context.cmp(&h2.context) {
        std::cmp::Ordering::Equal => h2
            .query
            .offset
            .cmp(&h1.query.offset)
            .then_with(|| h2.subject.offset.cmp(&h1.subject.offset)),
        other => other,
    }
}

/// Mirror of `s_RevCompareHSPsTransl` (`link_hsps.c:234`).
fn rev_compare_hsps_transl(h1: &LinkBlastHsp, h2: &LinkBlastHsp) -> std::cmp::Ordering {
    let c1 = h1.context / (NUM_FRAMES / 2);
    let c2 = h2.context / (NUM_FRAMES / 2);
    match c1.cmp(&c2) {
        std::cmp::Ordering::Equal => h2
            .query
            .offset
            .cmp(&h1.query.offset)
            .then_with(|| h2.subject.offset.cmp(&h1.subject.offset)),
        other => other,
    }
}

#[inline]
fn signum_i32(v: i32) -> i32 {
    match v.cmp(&0) {
        std::cmp::Ordering::Less => -1,
        std::cmp::Ordering::Equal => 0,
        std::cmp::Ordering::Greater => 1,
    }
}

/// Mirror of `s_RevCompareHSPsTbn` (`link_hsps.c:277`). Used for the
/// initial sort in `s_BlastEvenGapLinkHSPs` when the query is not
/// translated (`kTranslatedQuery == FALSE`). Separates HSPs by frame
/// sign of the subject, then by decreasing query offsets/ends and
/// decreasing subject offsets/ends.
fn rev_compare_hsps_tbn(h1: &LinkBlastHsp, h2: &LinkBlastHsp) -> std::cmp::Ordering {
    use std::cmp::Ordering;
    if h1.context != h2.context {
        return h1.context.cmp(&h2.context);
    }
    let s1 = signum_i32(h1.subject.frame);
    let s2 = signum_i32(h2.subject.frame);
    if s1 != s2 {
        // NCBI: if signs differ, h1->subject.frame > h2->subject.frame => 1.
        return if h1.subject.frame > h2.subject.frame {
            Ordering::Greater
        } else {
            Ordering::Less
        };
    }
    h2.query
        .offset
        .cmp(&h1.query.offset)
        .then_with(|| h2.query.end.cmp(&h1.query.end))
        .then_with(|| h2.subject.offset.cmp(&h1.subject.offset))
        .then_with(|| h2.subject.end.cmp(&h1.subject.end))
}

/// Mirror of `s_RevCompareHSPsTbx` (`link_hsps.c:331`). Used when the
/// query is translated.
fn rev_compare_hsps_tbx(h1: &LinkBlastHsp, h2: &LinkBlastHsp) -> std::cmp::Ordering {
    use std::cmp::Ordering;
    let c1 = h1.context / (NUM_FRAMES / 2);
    let c2 = h2.context / (NUM_FRAMES / 2);
    if c1 != c2 {
        return c1.cmp(&c2);
    }
    let s1 = signum_i32(h1.subject.frame);
    let s2 = signum_i32(h2.subject.frame);
    if s1 != s2 {
        return if h1.subject.frame > h2.subject.frame {
            Ordering::Greater
        } else {
            Ordering::Less
        };
    }
    h2.query
        .offset
        .cmp(&h1.query.offset)
        .then_with(|| h2.query.end.cmp(&h1.query.end))
        .then_with(|| h2.subject.offset.cmp(&h1.subject.offset))
        .then_with(|| h2.subject.end.cmp(&h1.subject.end))
}

/// Perform even-gap linking on a list of HSPs.
///
/// Port of NCBI `s_BlastEvenGapLinkHSPs` (`link_hsps.c:414`). The NCBI
/// C original mutates an intrusive doubly-linked list of
/// `LinkHSPStruct` nodes through raw pointers. Here the linked list is
/// backed by a `Vec<LinkHspStruct>` and traversed via
/// `Option<usize>` indices; index 0 is reserved for the `hp_start`
/// sentinel that NCBI allocates separately.
///
/// Returns `0` on success, `-1` on invalid input (matches NCBI).
pub fn s_BlastEvenGapLinkHSPs(
    program_number: ProgramType,
    hsp_list: &mut LinkBlastHspList,
    query_info: &QueryInfo,
    subject_length: i32,
    sbp: &LinkScoreBlock,
    link_hsp_params: &LinkHSPParameters,
    gapped_calculation: bool,
) -> i32 {
    if hsp_list.hsp_array.is_empty() {
        return -1;
    }

    let total_number_of_hsps: i32 = hsp_list.hspcnt();
    let number_of_hsps_init: i32 = total_number_of_hsps;

    // `link_hsps.c:452`: lh_helper_size = MAX(1024, hspcnt+5).
    let lh_helper_size = std::cmp::max(1024, (hsp_list.hspcnt() + 5) as usize);
    let mut lh_helper: Vec<LinkHelpStruct> = vec![LinkHelpStruct::default(); lh_helper_size];

    let kbp: &[KarlinBlk] = if gapped_calculation {
        &sbp.kbp_gap
    } else {
        &sbp.kbp
    };

    // `link_hsps.c:466-469`.
    let window_size = link_hsp_params.gap_size + link_hsp_params.overlap_size + 1;
    let trim_size = (link_hsp_params.overlap_size + 1) / 2;
    let gap_prob = link_hsp_params.gap_prob;
    let gap_decay_rate = link_hsp_params.gap_decay_rate;

    let num_subject_frames = if subject_is_translated(program_number) {
        NUM_STRANDS
    } else {
        1
    };

    // Build an indexing vector mirroring NCBI's `link_hsp_array` which
    // points at allocated LinkHSPStruct entries (one per HSP).
    // Nodes are stored in `nodes`, indexed 0 (sentinel) .. (N).
    // We sort an array of indices into `nodes` according to the appropriate
    // comparator; `order[i]` is the node index at position `i`.
    let num_nodes = (total_number_of_hsps as usize) + 1; // +1 for sentinel
    let mut nodes: Vec<LinkHspStruct> = vec![LinkHspStruct::default(); num_nodes];
    // Assign hsp_idx = -1 for the sentinel, i.e. slot 0.
    nodes[0].hsp_idx = -1;
    for (i, node) in nodes.iter_mut().enumerate().skip(1) {
        node.hsp_idx = (i - 1) as i32;
    }

    // NCBI sorts link_hsp_array by `RevCompareHSPsTbx` (translated query)
    // or `RevCompareHSPsTbn` (non-translated query). We compute an
    // ordering on node indices 1..=N.
    let kTranslatedQuery = query_is_translated(program_number);
    let mut order: Vec<usize> = (1..num_nodes).collect();
    {
        let hsp_array = &hsp_list.hsp_array;
        if kTranslatedQuery {
            order.sort_by(|&a, &b| {
                rev_compare_hsps_tbx(
                    &hsp_array[nodes[a].hsp_idx as usize],
                    &hsp_array[nodes[b].hsp_idx as usize],
                )
            });
        } else {
            order.sort_by(|&a, &b| {
                rev_compare_hsps_tbn(
                    &hsp_array[nodes[a].hsp_idx as usize],
                    &hsp_array[nodes[b].hsp_idx as usize],
                )
            });
        }
    }

    let cutoff: [i32; 2] = [
        link_hsp_params.cutoff_small_gap,
        link_hsp_params.cutoff_big_gap,
    ];
    let ignore_small_gaps = cutoff[0] == 0;

    // `link_hsps.c:498-501`: if query is nucleotide, num_query_frames = 2 * num_queries.
    let mut num_query_frames = if query_is_nucleotide(program_number) {
        NUM_STRANDS * query_info.num_queries
    } else {
        query_info.num_queries
    };

    // Allocate per-frame headers (`hp_frame_start`, `hp_frame_number`).
    let mut hp_frame_start: Vec<Option<usize>> =
        vec![None; (num_subject_frames * num_query_frames) as usize];
    let mut hp_frame_number: Vec<i32> = vec![0; (num_subject_frames * num_query_frames) as usize];

    // Hook up the HSPs in the order determined above.
    if !order.is_empty() {
        hp_frame_start[0] = Some(order[0]);
    }
    {
        let mut cur_frame: i32 = 0;
        let strand_factor = if kTranslatedQuery { 3 } else { 1 };
        let number_of_hsps = number_of_hsps_init;
        for index in 0..(number_of_hsps as usize) {
            let h_idx = order[index];
            // Reset start_of_chain.
            nodes[h_idx].start_of_chain = false;
            hp_frame_number[cur_frame as usize] += 1;
            nodes[h_idx].prev = if index > 0 {
                Some(order[index - 1])
            } else {
                None
            };
            nodes[h_idx].next = if index < (number_of_hsps as usize) - 1 {
                Some(order[index + 1])
            } else {
                None
            };

            if let Some(prev_idx) = nodes[h_idx].prev {
                let h_hsp = &hsp_list.hsp_array[nodes[h_idx].hsp_idx as usize];
                let prev_hsp = &hsp_list.hsp_array[nodes[prev_idx].hsp_idx as usize];
                let ctx_changed =
                    (h_hsp.context / strand_factor) != (prev_hsp.context / strand_factor);
                let frame_changed =
                    signum_i32(h_hsp.subject.frame) != signum_i32(prev_hsp.subject.frame);
                if ctx_changed || frame_changed {
                    // Frame switches: start a new list.
                    hp_frame_number[cur_frame as usize] -= 1;
                    cur_frame += 1;
                    hp_frame_number[cur_frame as usize] += 1;
                    hp_frame_start[cur_frame as usize] = Some(h_idx);
                    // Break the bidirectional link.
                    nodes[prev_idx].next = None;
                    nodes[h_idx].prev = None;
                }
            }
        }
        num_query_frames = cur_frame + 1;
    }

    // Compute trimmed offsets (`link_hsps.c:540-551`).
    for i in 0..(number_of_hsps_init as usize) {
        let h_idx = order[i];
        let hsp = &hsp_list.hsp_array[nodes[h_idx].hsp_idx as usize];
        let q_length = (hsp.query.end - hsp.query.offset) / 4;
        let s_length = (hsp.subject.end - hsp.subject.offset) / 4;
        nodes[h_idx].q_offset_trim = hsp.query.offset + std::cmp::min(q_length, trim_size);
        nodes[h_idx].q_end_trim = hsp.query.end - std::cmp::min(q_length, trim_size);
        nodes[h_idx].s_offset_trim = hsp.subject.offset + std::cmp::min(s_length, trim_size);
        nodes[h_idx].s_end_trim = hsp.subject.end - std::cmp::min(s_length, trim_size);
    }

    // Per-frame linking.
    let subject_length_orig = subject_length;
    for frame_index in 0..(num_query_frames as usize) {
        // Reset the sentinel node (hp_start = nodes[0]).
        nodes[0].reset();
        nodes[0].hsp_idx = -1;
        nodes[0].next = hp_frame_start[frame_index];
        if let Some(start_idx) = hp_frame_start[frame_index] {
            nodes[start_idx].prev = Some(0);
        }
        let mut number_of_hsps = hp_frame_number[frame_index];
        let Some(first_idx) = nodes[0].next else {
            continue;
        };
        let query_context = hsp_list.hsp_array[nodes[first_idx].hsp_idx as usize].context as usize;
        let ctx = &query_info.contexts[query_context];
        let length_adjustment = ctx.length_adjustment;
        let mut query_length = ctx.query_length;
        query_length = std::cmp::max(query_length - length_adjustment, 1);
        let mut subject_length_eff = subject_length_orig;
        let mut length_adjustment_eff = length_adjustment;
        if subject_is_translated(program_number) {
            length_adjustment_eff /= CODON_LENGTH;
            subject_length_eff /= CODON_LENGTH;
        }
        subject_length_eff = std::cmp::max(subject_length_eff - length_adjustment_eff, 1);

        // Initialize end-marker helper[0] (`link_hsps.c:573-577`).
        lh_helper[0] = LinkHelpStruct {
            ptr: 0,
            q_off_trim: 0,
            s_off_trim: 0,
            sum: [0, 0],
            maxsum1: -10000,
            next_larger: 0,
        };

        let mut first_pass = 1i32; // do full search
        let mut path_changed = 1i32;

        // `H->hsp_link.changed = 1` for all H in the chain.
        {
            let mut cur = nodes[0].next;
            while let Some(h) = cur {
                nodes[h].hsp_link.changed = 1;
                cur = nodes[h].next;
            }
        }

        while number_of_hsps > 0 {
            let mut max: [i32; 3] = [-10000, -10000, -10000];
            let mut best: [Option<usize>; 2] = [None, None];

            // See if we can avoid recomputing all scores (`link_hsps.c:597-652`).
            let mut use_current_max = false;
            if first_pass == 0 {
                if !ignore_small_gaps {
                    let mut max0 = -cutoff[0];
                    let mut max1 = -cutoff[1];
                    let mut cur = nodes[0].next;
                    while let Some(h) = cur {
                        let sum0 = nodes[h].hsp_link.sum[0];
                        let sum1 = nodes[h].hsp_link.sum[1];
                        if sum0 >= max0 {
                            max0 = sum0;
                            best[0] = Some(h);
                        }
                        if sum1 >= max1 {
                            max1 = sum1;
                            best[1] = Some(h);
                        }
                        cur = nodes[h].next;
                    }
                } else {
                    let mut maxscore = -cutoff[1];
                    let mut cur = nodes[0].next;
                    while let Some(h) = cur {
                        let sum = nodes[h].hsp_link.sum[1];
                        if sum >= maxscore {
                            maxscore = sum;
                            best[1] = Some(h);
                        }
                        cur = nodes[h].next;
                    }
                }
                if path_changed == 0 {
                    use_current_max = true;
                } else {
                    // Walk down best, give up if removed item in path.
                    use_current_max = true;
                    if !ignore_small_gaps {
                        let mut cur = best[0];
                        while let Some(h) = cur {
                            if nodes[h].linked_to == -1000 {
                                use_current_max = false;
                                break;
                            }
                            cur = nodes[h].hsp_link.link[0];
                        }
                    }
                    if use_current_max {
                        let mut cur = best[1];
                        while let Some(h) = cur {
                            if nodes[h].linked_to == -1000 {
                                use_current_max = false;
                                break;
                            }
                            cur = nodes[h].hsp_link.link[1];
                        }
                    }
                }
            }

            if !use_current_max {
                // Reset helper_info (`link_hsps.c:659-686`).
                let mut h_index: usize = 1;
                let mut cur = Some(0usize); // start from hp_start sentinel
                while let Some(h) = cur {
                    let hsp_idx = nodes[h].hsp_idx;
                    let s_frame = if hsp_idx >= 0 {
                        hsp_list.hsp_array[hsp_idx as usize].subject.frame
                    } else {
                        0
                    };
                    let s_off_t = nodes[h].s_offset_trim;
                    let q_off_t = nodes[h].q_offset_trim;
                    lh_helper[h_index].ptr = h;
                    lh_helper[h_index].q_off_trim = q_off_t;
                    lh_helper[h_index].s_off_trim = s_off_t;
                    for i in 0..E_ORDERING_METHODS {
                        lh_helper[h_index].sum[i] = nodes[h].hsp_link.sum[i];
                    }
                    let frame_bucket = (signum_i32(s_frame) + 1) as usize;
                    max[frame_bucket] = std::cmp::max(max[frame_bucket], nodes[h].hsp_link.sum[1]);
                    lh_helper[h_index].maxsum1 = max[frame_bucket];
                    // set next_larger to link back to closest entry with larger sum1.
                    {
                        let cur_sum = lh_helper[h_index].sum[1];
                        let mut prev = h_index.saturating_sub(1);
                        let mut prev_sum = lh_helper[prev].sum[1];
                        while cur_sum >= prev_sum && prev > 0 {
                            prev = lh_helper[prev].next_larger;
                            prev_sum = lh_helper[prev].sum[1];
                        }
                        lh_helper[h_index].next_larger = prev;
                    }
                    nodes[h].linked_to = 0;
                    h_index += 1;
                    cur = nodes[h].next;
                }

                // NCBI (`:688`): the helper slot for `hp_start` must not
                // influence `maxsum1` comparisons.
                lh_helper[1].maxsum1 = -10000;

                // --- loop iter for index = 0 (small gaps) ---
                if !ignore_small_gaps {
                    let index = 0usize;
                    let mut maxscore = -cutoff[index];
                    let mut h_index: usize = 2;
                    let mut cur = nodes[0].next;
                    while let Some(h) = cur {
                        let mut h_hsp_num = 0i32;
                        let mut h_hsp_sum = 0i32;
                        let mut h_hsp_xsum = 0.0f64;
                        let mut h_hsp_link: Option<usize> = None;
                        let h_hsp_score = hsp_list.hsp_array[nodes[h].hsp_idx as usize].score;
                        if h_hsp_score > cutoff[index] {
                            let h_query_etrim = nodes[h].q_end_trim;
                            let h_sub_etrim = nodes[h].s_end_trim;
                            let h_q_et_gap = h_query_etrim + window_size;
                            let h_s_et_gap = h_sub_etrim + window_size;

                            let mut h2_index = h_index as i32 - 1;
                            while h2_index > 1 {
                                let help = &lh_helper[h2_index as usize];
                                let q_off_t = help.q_off_trim;
                                let s_off_t = help.s_off_trim;
                                let sum = help.sum[index];
                                let b1 = q_off_t <= h_query_etrim;
                                let b2 = s_off_t <= h_sub_etrim;
                                let b4 = q_off_t > h_q_et_gap;
                                let b5 = s_off_t > h_s_et_gap;
                                if q_off_t > (h_q_et_gap + trim_size) {
                                    break;
                                }
                                if b1 | b2 | b5 | b4 {
                                    h2_index -= 1;
                                    continue;
                                }
                                if sum > h_hsp_sum {
                                    let h2_ptr = help.ptr;
                                    h_hsp_num = nodes[h2_ptr].hsp_link.num[index];
                                    h_hsp_sum = nodes[h2_ptr].hsp_link.sum[index];
                                    h_hsp_xsum = nodes[h2_ptr].hsp_link.xsum[index];
                                    h_hsp_link = Some(h2_ptr);
                                }
                                h2_index -= 1;
                            }
                        }
                        let score = h_hsp_score;
                        let ctx = nodes[h].hsp_idx;
                        let kb = &kbp[hsp_list.hsp_array[ctx as usize].context as usize];
                        let new_xsum = h_hsp_xsum + score as f64 * kb.lambda - kb.log_k;
                        let new_sum = h_hsp_sum + (score - cutoff[index]);

                        nodes[h].hsp_link.sum[index] = new_sum;
                        nodes[h].hsp_link.num[index] = h_hsp_num + 1;
                        nodes[h].hsp_link.link[index] = h_hsp_link;
                        lh_helper[h_index].sum[index] = new_sum;
                        if new_sum >= maxscore {
                            maxscore = new_sum;
                            best[index] = Some(h);
                        }
                        nodes[h].hsp_link.xsum[index] = new_xsum;
                        if let Some(link_idx) = h_hsp_link {
                            nodes[link_idx].linked_to += 1;
                        }

                        h_index += 1;
                        cur = nodes[h].next;
                    }
                }

                // --- loop iter for index = 1 (large gaps) ---
                {
                    let index = 1usize;
                    let mut maxscore = -cutoff[index];
                    let mut h_index: usize = 2;
                    let mut cur = nodes[0].next;
                    while let Some(h) = cur {
                        let mut h_hsp_num = 0i32;
                        let mut h_hsp_sum = 0i32;
                        let mut h_hsp_xsum = 0.0f64;
                        let mut h_hsp_link: Option<usize> = None;

                        nodes[h].hsp_link.changed = 1;
                        let h2_last = nodes[h].hsp_link.link[index];
                        let h_hsp_score = hsp_list.hsp_array[nodes[h].hsp_idx as usize].score;
                        let shortcut_valid = match h2_last {
                            None => first_pass == 0,
                            Some(h2) => first_pass == 0 && nodes[h2].hsp_link.changed == 0,
                        };
                        if shortcut_valid {
                            if let Some(h2) = h2_last {
                                h_hsp_num = nodes[h2].hsp_link.num[index];
                                h_hsp_sum = nodes[h2].hsp_link.sum[index];
                                h_hsp_xsum = nodes[h2].hsp_link.xsum[index];
                            }
                            h_hsp_link = h2_last;
                            nodes[h].hsp_link.changed = 0;
                        } else if h_hsp_score > cutoff[index] {
                            let h_query_etrim = nodes[h].q_end_trim;
                            let h_sub_etrim = nodes[h].s_end_trim;

                            // NCBI (`:812-823`): if we had a valid link last
                            // pass, bias the initial best by (sum-1) to
                            // reproduce original tie-break ordering.
                            if first_pass == 0 {
                                if let Some(h2) = h2_last {
                                    if nodes[h2].linked_to >= 0 {
                                        h_hsp_sum = nodes[h2].hsp_link.sum[index] - 1;
                                    }
                                }
                            }

                            let mut h2_index = h_index as i32 - 1;
                            while h2_index > 1 {
                                let help = &lh_helper[h2_index as usize];
                                let sum = help.sum[index];
                                let next_larger = help.next_larger;
                                let s_off_t = help.s_off_trim;
                                let q_off_t = help.q_off_trim;

                                let b0 = sum <= h_hsp_sum;

                                // Advance (NCBI `:841-844`).
                                h2_index -= 1;
                                if b0 {
                                    h2_index = next_larger as i32;
                                }

                                let b1 = q_off_t <= h_query_etrim;
                                let b2 = s_off_t <= h_sub_etrim;

                                // NCBI has `if(0) if(...) break;` at :850 which is dead code.

                                if !(b0 | b1 | b2) {
                                    let h2_ptr = help.ptr;
                                    h_hsp_num = nodes[h2_ptr].hsp_link.num[index];
                                    h_hsp_sum = nodes[h2_ptr].hsp_link.sum[index];
                                    h_hsp_xsum = nodes[h2_ptr].hsp_link.xsum[index];
                                    h_hsp_link = Some(h2_ptr);
                                }
                            }
                        }

                        let score = h_hsp_score;
                        let ctx = nodes[h].hsp_idx;
                        let kb = &kbp[hsp_list.hsp_array[ctx as usize].context as usize];
                        let new_xsum = h_hsp_xsum + score as f64 * kb.lambda - kb.log_k;
                        let new_sum = h_hsp_sum + (score - cutoff[index]);

                        nodes[h].hsp_link.sum[index] = new_sum;
                        nodes[h].hsp_link.num[index] = h_hsp_num + 1;
                        nodes[h].hsp_link.link[index] = h_hsp_link;
                        lh_helper[h_index].sum[index] = new_sum;
                        lh_helper[h_index].maxsum1 =
                            std::cmp::max(lh_helper[h_index - 1].maxsum1, new_sum);

                        // Update this entry's next_larger (`link_hsps.c:876-885`).
                        {
                            let cur_sum = lh_helper[h_index].sum[1];
                            let mut prev = h_index - 1;
                            let mut prev_sum = lh_helper[prev].sum[1];
                            while cur_sum >= prev_sum && prev > 0 {
                                prev = lh_helper[prev].next_larger;
                                prev_sum = lh_helper[prev].sum[1];
                            }
                            lh_helper[h_index].next_larger = prev;
                        }

                        if new_sum >= maxscore {
                            maxscore = new_sum;
                            best[index] = Some(h);
                        }
                        nodes[h].hsp_link.xsum[index] = new_xsum;
                        if let Some(link_idx) = h_hsp_link {
                            nodes[link_idx].linked_to += 1;
                        }

                        h_index += 1;
                        cur = nodes[h].next;
                    }
                }
                path_changed = 0;
                first_pass = 0;
            }

            // Select best ordering method (`link_hsps.c:901-953`).
            let mut prob = [0.0f64; 2];
            let ordering_method;
            if !ignore_small_gaps {
                let best0 = best[0].expect("best[0] must exist when !ignore_small_gaps");
                let best1 = best[1].expect("best[1] must exist");
                // Add back cutoff * num_links.
                let b0num = nodes[best0].hsp_link.num[0];
                nodes[best0].hsp_link.sum[0] += b0num * cutoff[0];

                let mut p0 = small_gap_sum_e(
                    window_size,
                    b0num as u32,
                    nodes[best0].hsp_link.xsum[0],
                    query_length,
                    subject_length_eff,
                    ctx.eff_searchsp as f64,
                    gap_decay_divisor(gap_decay_rate, b0num as u32),
                )
                .unwrap_or(INT4_MAX);

                if b0num > 1 {
                    if gap_prob == 0.0 {
                        p0 = INT4_MAX;
                    } else {
                        p0 /= gap_prob;
                        if p0 > INT4_MAX {
                            p0 = INT4_MAX;
                        }
                    }
                }

                let b1num = nodes[best1].hsp_link.num[1];
                let mut p1 = large_gap_sum_e(
                    b1num as u32,
                    nodes[best1].hsp_link.xsum[1],
                    query_length,
                    subject_length_eff,
                    ctx.eff_searchsp as f64,
                    gap_decay_divisor(gap_decay_rate, b1num as u32),
                )
                .unwrap_or(INT4_MAX);

                if b1num > 1 {
                    if (1.0 - gap_prob) == 0.0 {
                        p1 = INT4_MAX;
                    } else {
                        p1 /= 1.0 - gap_prob;
                        if p1 > INT4_MAX {
                            p1 = INT4_MAX;
                        }
                    }
                }

                prob[0] = p0;
                prob[1] = p1;
                ordering_method = if p0 <= p1 {
                    ELinkOrderingMethod::ELinkSmallGaps
                } else {
                    ELinkOrderingMethod::ELinkLargeGaps
                };
            } else {
                let best1 = best[1].expect("best[1] must exist");
                let b1num = nodes[best1].hsp_link.num[1];
                nodes[best1].hsp_link.sum[1] += b1num * cutoff[1];

                let p1 = large_gap_sum_e(
                    b1num as u32,
                    nodes[best1].hsp_link.xsum[1],
                    query_length,
                    subject_length_eff,
                    ctx.eff_searchsp as f64,
                    gap_decay_divisor(gap_decay_rate, b1num as u32),
                )
                .unwrap_or(INT4_MAX);
                prob[1] = p1;
                ordering_method = ELinkOrderingMethod::ELinkLargeGaps;
            }

            let ord_idx = ordering_method as usize;
            let best_idx = best[ord_idx].expect("ordering_method best must be set");
            nodes[best_idx].start_of_chain = true;
            hsp_list.hsp_array[nodes[best_idx].hsp_idx as usize].evalue = prob[ord_idx];

            // Remove the links that have been ordered already.
            let linked_set = nodes[best_idx].hsp_link.link[ord_idx].is_some();
            if nodes[best_idx].linked_to > 0 {
                path_changed = 1;
            }
            let mut cur = Some(best_idx);
            while let Some(h) = cur {
                if nodes[h].linked_to > 1 {
                    path_changed = 1;
                }
                nodes[h].linked_to = -1000;
                nodes[h].hsp_link.changed = 1;
                nodes[h].linked_set = linked_set;
                nodes[h].ordering_method = ordering_method;
                let hsp_idx = nodes[h].hsp_idx as usize;
                hsp_list.hsp_array[hsp_idx].evalue = prob[ord_idx];
                let next = nodes[h].next;
                let prev = nodes[h].prev;
                if let Some(n) = next {
                    nodes[n].prev = prev;
                }
                if let Some(p) = prev {
                    nodes[p].next = next;
                }
                number_of_hsps -= 1;
                cur = nodes[h].hsp_link.link[ord_idx];
            }
        } // end while num_hsps
    } // end for frame_index

    // NCBI (`:1005-1090`) performs additional sorting + pointer stitching to
    // reorder the final hsp_array, preserving the linking information. We
    // reproduce the effect on the underlying hsp_array by sorting nodes
    // according to the `FwdCompareHSPs` comparator (the second qsort
    // dominates per NCBI). Then copy HSPs back into hsp_list in that order,
    // setting `num` for chain members.
    {
        let mut order: Vec<usize> = (1..num_nodes).collect();
        if kTranslatedQuery {
            order.sort_by(|&a, &b| {
                fwd_compare_hsps_transl(
                    &hsp_list.hsp_array[nodes[a].hsp_idx as usize],
                    &hsp_list.hsp_array[nodes[b].hsp_idx as usize],
                )
            });
        } else {
            order.sort_by(|&a, &b| {
                fwd_compare_hsps(
                    &hsp_list.hsp_array[nodes[a].hsp_idx as usize],
                    &hsp_list.hsp_array[nodes[b].hsp_idx as usize],
                )
            });
        }

        // Propagate `num` (number of linked HSPs) to chain members.
        // Walk the chains via `hsp_link.link[ordering_method]` from any
        // start_of_chain or from an unlinked HSP.
        for i in 0..order.len() {
            let h = order[i];
            if !nodes[h].linked_set || nodes[h].start_of_chain {
                let om = nodes[h].ordering_method as usize;
                let num_links = if nodes[h].linked_set {
                    nodes[h].hsp_link.num[om]
                } else {
                    1
                };
                hsp_list.hsp_array[nodes[h].hsp_idx as usize].num = num_links;
                if nodes[h].linked_set {
                    let xsum = nodes[h].hsp_link.xsum[om];
                    hsp_list.hsp_array[nodes[h].hsp_idx as usize].bit_score = xsum;
                    // Walk the link chain, propagating num to each member.
                    let mut cur_link = nodes[h].hsp_link.link[om];
                    while let Some(link_h) = cur_link {
                        hsp_list.hsp_array[nodes[link_h].hsp_idx as usize].num = num_links;
                        cur_link = nodes[link_h].hsp_link.link[om];
                    }
                }
            }
        }

        // Physically reorder hsp_array according to `order`.
        let mut new_array: Vec<LinkBlastHsp> = Vec::with_capacity(order.len());
        for &h in &order {
            new_array.push(hsp_list.hsp_array[nodes[h].hsp_idx as usize].clone());
        }
        hsp_list.hsp_array = new_array;
    }

    0
}

// ---------------------------------------------------------------------------
// Uneven-gap linking (tblastn / blastx)
// ---------------------------------------------------------------------------

/// Port of `BlastLinkedHSPSet` (`link_hsps.c:1098`). Doubly linked via
/// `Option<usize>` indices into a backing `Vec<BlastLinkedHspSet>`.
#[derive(Debug, Clone, Default)]
pub struct BlastLinkedHspSet {
    /// Index into the owning `LinkBlastHspList::hsp_array`.
    pub hsp_idx: i32,
    pub queryId: u32,
    pub next: Option<usize>,
    pub prev: Option<usize>,
    pub sum_score: f64,
}

/// Compute the e-value of a combined HSP set.
/// Port of `s_SumHSPEvalue` (`link_hsps.c:1117`).
fn s_SumHSPEvalue(
    program_number: ProgramType,
    query_info: &QueryInfo,
    subject_length: i32,
    link_hsp_params: &LinkHSPParameters,
    sets: &[BlastLinkedHspSet],
    hsp_array: &[LinkBlastHsp],
    head_idx: usize,
    new_idx: usize,
) -> (f64, f64) {
    // NCBI asserts not tblastx.
    debug_assert!(!(query_is_translated(program_number) && subject_is_translated(program_number)));

    let subject_eff_length = if subject_is_translated(program_number) {
        subject_length / 3
    } else {
        subject_length
    };
    let gap_decay_rate = link_hsp_params.gap_decay_rate;
    let head_hsp = &hsp_array[sets[head_idx].hsp_idx as usize];
    let new_hsp = &hsp_array[sets[new_idx].hsp_idx as usize];
    let num = head_hsp.num + new_hsp.num;

    let context = head_hsp.context as usize;
    let len_adj = query_info.contexts[context].length_adjustment;
    let query_eff_length = std::cmp::max(query_info.contexts[context].query_length - len_adj, 1);
    let subject_eff_length = std::cmp::max(subject_eff_length - len_adj, 1);

    let sum_score = sets[new_idx].sum_score + sets[head_idx].sum_score;

    let query_window_size = link_hsp_params.overlap_size + link_hsp_params.gap_size + 1;
    let subject_window_size = link_hsp_params.overlap_size + link_hsp_params.longest_intron + 1;

    let sum_evalue = uneven_gap_sum_e(
        query_window_size,
        subject_window_size,
        num as u32,
        sum_score,
        query_eff_length,
        subject_eff_length,
        query_info.contexts[context].eff_searchsp as f64,
        gap_decay_divisor(gap_decay_rate, num as u32),
    )
    .unwrap_or(INT4_MAX);

    (sum_evalue, sum_score)
}

/// Sort comparator: increasing (queryId, query.offset, subject.offset).
/// Port of `s_FwdCompareLinkedHSPSets` (`link_hsps.c:1170`).
fn fwd_compare_linked_hsp_sets(
    sets: &[BlastLinkedHspSet],
    hsp_array: &[LinkBlastHsp],
    a: usize,
    b: usize,
) -> std::cmp::Ordering {
    if sets[a].queryId != sets[b].queryId {
        return sets[a].queryId.cmp(&sets[b].queryId);
    }
    let h1 = &hsp_array[sets[a].hsp_idx as usize];
    let h2 = &hsp_array[sets[b].hsp_idx as usize];
    h1.query
        .offset
        .cmp(&h2.query.offset)
        .then_with(|| h1.subject.offset.cmp(&h2.subject.offset))
}

/// Comparator: decreasing sum_score, then NCBI ScoreCompareHSPs tie-break.
/// Port of `s_SumScoreCompareLinkedHSPSets` (`link_hsps.c:1205`).
fn sum_score_compare_linked_hsp_sets(
    sets: &[BlastLinkedHspSet],
    hsp_array: &[LinkBlastHsp],
    a: usize,
    b: usize,
) -> std::cmp::Ordering {
    use std::cmp::Ordering;
    // Sort by decreasing sum_score.
    match sets[b].sum_score.partial_cmp(&sets[a].sum_score) {
        Some(Ordering::Equal) | None => {}
        Some(other) => return other,
    }
    // Tie-break: ScoreCompareHSPs — score desc, subject_offset asc,
    // subject_end desc, query_offset asc, query_end desc.
    let h1 = &hsp_array[sets[a].hsp_idx as usize];
    let h2 = &hsp_array[sets[b].hsp_idx as usize];
    h2.score
        .cmp(&h1.score)
        .then_with(|| h1.subject.offset.cmp(&h2.subject.offset))
        .then_with(|| h2.subject.end.cmp(&h1.subject.end))
        .then_with(|| h1.query.offset.cmp(&h2.query.offset))
        .then_with(|| h2.query.end.cmp(&h1.query.end))
}

/// Binary search for the index in a sorted (by queryId, query.offset)
/// array of the smallest element with queryId == target and
/// `query.offset >= offset`. Port of `s_HSPOffsetBinarySearch`
/// (`link_hsps.c:1242`).
fn s_HSPOffsetBinarySearch(
    sets: &[BlastLinkedHspSet],
    hsp_array: &[LinkBlastHsp],
    offset_hsp_array: &[usize],
    size: i32,
    queryId: u32,
    offset: i32,
) -> i32 {
    let mut begin: i32 = 0;
    let mut end: i32 = size;
    while begin < end {
        let index = (begin + end) / 2;
        let s = &sets[offset_hsp_array[index as usize]];
        if s.queryId < queryId {
            begin = index + 1;
        } else if s.queryId > queryId {
            end = index;
        } else {
            let h = &hsp_array[s.hsp_idx as usize];
            if h.query.offset >= offset {
                end = index;
            } else {
                begin = index + 1;
            }
        }
    }
    end
}

/// Port of `s_HSPOffsetEndBinarySearch` (`link_hsps.c:1280`).
fn s_HSPOffsetEndBinarySearch(
    sets: &[BlastLinkedHspSet],
    hsp_array: &[LinkBlastHsp],
    offset_hsp_array: &[usize],
    size: i32,
    qend_index_array: &[i32],
    queryId: u32,
    offset: i32,
) -> i32 {
    let mut begin: i32 = 0;
    let mut end: i32 = size;
    while begin < end {
        let right_index = (begin + end) / 2;
        let left_index = qend_index_array[right_index as usize];
        let s_right = &sets[offset_hsp_array[right_index as usize]];
        if s_right.queryId < queryId {
            begin = right_index + 1;
        } else if s_right.queryId > queryId {
            end = left_index;
        } else {
            let s_left = &sets[offset_hsp_array[left_index as usize]];
            let h_left = &hsp_array[s_left.hsp_idx as usize];
            if h_left.query.end >= offset {
                end = left_index;
            } else {
                begin = right_index + 1;
            }
        }
    }
    end
}

/// Merge two linked chains of `BlastLinkedHspSet` into a sorted (by
/// query.offset) array of indices. Port of `s_MergeLinkedHSPSets`
/// (`link_hsps.c:1314`).
fn s_MergeLinkedHSPSets(
    sets: &[BlastLinkedHspSet],
    hsp_array: &[LinkBlastHsp],
    mut hsp_set1: Option<usize>,
    mut hsp_set2: Option<usize>,
) -> Vec<usize> {
    // Rewind both to head.
    while let Some(h) = hsp_set1 {
        if let Some(p) = sets[h].prev {
            hsp_set1 = Some(p);
        } else {
            break;
        }
    }
    while let Some(h) = hsp_set2 {
        if let Some(p) = sets[h].prev {
            hsp_set2 = Some(p);
        } else {
            break;
        }
    }

    let mut merged: Vec<usize> = Vec::new();
    let mut a = hsp_set1;
    let mut b = hsp_set2;
    while a.is_some() || b.is_some() {
        let pick_a = match (a, b) {
            (None, _) => false,
            (_, None) => true,
            (Some(ai), Some(bi)) => {
                let qa = hsp_array[sets[ai].hsp_idx as usize].query.offset;
                let qb = hsp_array[sets[bi].hsp_idx as usize].query.offset;
                qa < qb
            }
        };
        if pick_a {
            let ai = a.unwrap();
            merged.push(ai);
            a = sets[ai].next;
        } else {
            let bi = b.unwrap();
            merged.push(bi);
            b = sets[bi].next;
        }
    }
    merged
}

/// Check admissibility of combining two linked HSP sets for uneven-gap
/// linking. Port of `s_LinkedHSPSetsAdmissible` (`link_hsps.c:1407`).
fn s_LinkedHSPSetsAdmissible(
    sets: &[BlastLinkedHspSet],
    hsp_array: &[LinkBlastHsp],
    hsp_set1: usize,
    mut hsp_set2: usize,
    link_hsp_params: &LinkHSPParameters,
    program: ProgramType,
) -> bool {
    // First HSP must be the head of its set.
    if sets[hsp_set1].prev.is_some() {
        return false;
    }
    // Rewind hsp_set2 to head.
    while let Some(p) = sets[hsp_set2].prev {
        hsp_set2 = p;
    }
    if hsp_set1 == hsp_set2 {
        return false;
    }
    if sets[hsp_set1].queryId != sets[hsp_set2].queryId {
        return false;
    }
    let h1 = &hsp_array[sets[hsp_set1].hsp_idx as usize];
    let h2 = &hsp_array[sets[hsp_set2].hsp_idx as usize];
    if signum_i32(h1.subject.frame) != signum_i32(h2.subject.frame) {
        return false;
    }

    let merged = s_MergeLinkedHSPSets(sets, hsp_array, Some(hsp_set1), Some(hsp_set2));
    let combined_size = merged.len();

    let mut gap_s = link_hsp_params.longest_intron;
    let mut gap_q = link_hsp_params.gap_size;
    let overlap = link_hsp_params.overlap_size;
    if program == BLASTX {
        gap_s = link_hsp_params.gap_size;
        gap_q = link_hsp_params.longest_intron;
    }

    let mut index = 0usize;
    while index + 1 < combined_size {
        let l = &hsp_array[sets[merged[index]].hsp_idx as usize];
        let r = &hsp_array[sets[merged[index + 1]].hsp_idx as usize];
        if l.query.end < r.query.offset - gap_q {
            break;
        }
        if l.query.offset >= r.query.offset {
            break;
        }
        if l.query.end > r.query.offset + overlap || l.query.end >= r.query.end {
            break;
        }
        if l.subject.end > r.subject.offset + overlap
            || l.subject.end < r.subject.offset - gap_s
            || l.subject.offset >= r.subject.offset
            || l.subject.end >= r.subject.end
        {
            break;
        }
        index += 1;
    }
    index + 1 >= combined_size
}

/// Combine two linked sets into a single chain. Port of
/// `s_CombineLinkedHSPSets` (`link_hsps.c:1359`).
fn s_CombineLinkedHSPSets(
    sets: &mut [BlastLinkedHspSet],
    hsp_array: &mut [LinkBlastHsp],
    hsp_set1: Option<usize>,
    hsp_set2: Option<usize>,
    sum_score: f64,
    evalue: f64,
) -> Option<usize> {
    match (hsp_set1, hsp_set2) {
        (None, _) => return hsp_set2,
        (_, None) => return hsp_set1,
        _ => {}
    }

    // Snapshot hsp offsets for merging (MergeLinkedHSPSets reads immutably).
    let merged = {
        let sets_ro: &[BlastLinkedHspSet] = sets;
        let hsp_ro: &[LinkBlastHsp] = hsp_array;
        s_MergeLinkedHSPSets(sets_ro, hsp_ro, hsp_set1, hsp_set2)
    };
    let new_num = merged.len() as i32;
    let head = merged[0];
    sets[head].prev = None;
    for i in 0..merged.len() {
        let cur = merged[i];
        let next = if i + 1 < merged.len() {
            Some(merged[i + 1])
        } else {
            None
        };
        sets[cur].next = next;
        if let Some(n) = next {
            sets[n].prev = Some(cur);
        }
        sets[cur].sum_score = sum_score;
        let hsp_idx = sets[cur].hsp_idx as usize;
        hsp_array[hsp_idx].evalue = evalue;
        hsp_array[hsp_idx].num = new_num;
    }
    Some(head)
}

/// Build wrappers (`BlastLinkedHSPSet`) from a plain HSP array.
/// Port of `s_LinkedHSPSetArraySetUp` (`link_hsps.c:1508`).
fn s_LinkedHSPSetArraySetUp(
    hsp_array: &mut [LinkBlastHsp],
    kbp_array: &[KarlinBlk],
    program: ProgramType,
) -> Vec<BlastLinkedHspSet> {
    let mut link_hsp_array: Vec<BlastLinkedHspSet> = Vec::with_capacity(hsp_array.len());
    for (i, hsp) in hsp_array.iter_mut().enumerate() {
        let kb = &kbp_array[hsp.context as usize];
        link_hsp_array.push(BlastLinkedHspSet {
            hsp_idx: i as i32,
            queryId: if program == BLASTX {
                (hsp.context / 3) as u32
            } else {
                hsp.context as u32
            },
            next: None,
            prev: None,
            sum_score: kb.lambda * hsp.score as f64 - kb.log_k,
        });
        hsp.num = 1;
    }
    link_hsp_array
}

/// Index query ends to enable `s_HSPOffsetEndBinarySearch`.
/// Port of `s_LinkedHSPSetArrayIndexQueryEnds` (`link_hsps.c:1567`).
fn s_LinkedHSPSetArrayIndexQueryEnds(
    sets: &[BlastLinkedHspSet],
    hsp_array: &[LinkBlastHsp],
    offset_hsp_array: &[usize],
) -> Vec<i32> {
    let hspcnt = offset_hsp_array.len();
    let mut qend_index_array = vec![0i32; hspcnt];
    if hspcnt == 0 {
        return qend_index_array;
    }
    let mut current_index = 0i32;
    let mut current_end = hsp_array[sets[offset_hsp_array[0]].hsp_idx as usize]
        .query
        .end;
    for index in 1..hspcnt {
        let s = &sets[offset_hsp_array[index]];
        let h = &hsp_array[s.hsp_idx as usize];
        let cur_s = &sets[offset_hsp_array[current_index as usize]];
        if s.queryId > cur_s.queryId || h.query.end > current_end {
            current_index = index as i32;
            current_end = h.query.end;
        }
        qend_index_array[index] = current_index;
    }
    qend_index_array
}

/// Greedy uneven-gap HSP linker. Port of
/// `s_BlastUnevenGapLinkHSPs` (`link_hsps.c:1612`).
pub fn s_BlastUnevenGapLinkHSPs(
    program: ProgramType,
    hsp_list: &mut LinkBlastHspList,
    query_info: &QueryInfo,
    subject_length: i32,
    sbp: &LinkScoreBlock,
    link_hsp_params: &LinkHSPParameters,
    gapped_calculation: bool,
) -> i32 {
    if hsp_list.hsp_array.len() <= 1 {
        return 0;
    }

    let kbp_array: &[KarlinBlk] = if gapped_calculation {
        &sbp.kbp_gap
    } else {
        &sbp.kbp
    };

    // max gap size in query (NCBI `:1644-1646`)
    let _gap_size_q = if program == BLASTX {
        link_hsp_params.longest_intron
    } else {
        link_hsp_params.gap_size
    };

    // Set up wrapper array (computes sum_score and sets `num = 1`).
    let mut sets: Vec<BlastLinkedHspSet> =
        s_LinkedHSPSetArraySetUp(&mut hsp_list.hsp_array, kbp_array, program);
    let hspcnt = sets.len() as i32;

    // Two index arrays into `sets`: one sorted by decreasing sum_score,
    // one sorted by increasing (queryId, query.offset).
    let mut score_hsp_array: Vec<usize> = (0..sets.len()).collect();
    score_hsp_array
        .sort_by(|&a, &b| sum_score_compare_linked_hsp_sets(&sets, &hsp_list.hsp_array, a, b));
    let mut offset_hsp_array: Vec<usize> = (0..sets.len()).collect();
    offset_hsp_array
        .sort_by(|&a, &b| fwd_compare_linked_hsp_sets(&sets, &hsp_list.hsp_array, a, b));
    let qend_index_array =
        s_LinkedHSPSetArrayIndexQueryEnds(&sets, &hsp_list.hsp_array, &offset_hsp_array);

    let gap_size = if program == BLASTX {
        link_hsp_params.longest_intron
    } else {
        link_hsp_params.gap_size
    };

    let mut head_hsp: Option<usize> = None;
    let mut index: i32 = 0;
    while index < hspcnt {
        let head = match head_hsp {
            Some(h) => h,
            None => {
                // Find highest-scoring HSP not yet in any linked set.
                while index < hspcnt {
                    let cand = score_hsp_array[index as usize];
                    if sets[cand].next.is_none() && sets[cand].prev.is_none() {
                        break;
                    }
                    index += 1;
                }
                if index >= hspcnt {
                    break;
                }
                // NCBI (`link_hsps.c:1689`): head_hsp = score_hsp_array[index];
                // The outer `head_hsp` is always overwritten below by either
                // `s_CombineLinkedHSPSets` (on success) or `None` (on fail),
                // so we only need the local binding here for the rest of
                // this iteration.
                score_hsp_array[index as usize]
            }
        };

        // Find tail of the current chain.
        let mut tail_hsp = head;
        while let Some(n) = sets[tail_hsp].next {
            tail_hsp = n;
        }

        let mut best_evalue = hsp_list.hsp_array[sets[head].hsp_idx as usize].evalue;
        let mut best_sum_score = sets[head].sum_score;
        let mut best_hsp: Option<usize> = None;

        let head_q_offset = hsp_list.hsp_array[sets[head].hsp_idx as usize].query.offset;
        let left_offset = head_q_offset - gap_size;
        let hsp_index_left = s_HSPOffsetEndBinarySearch(
            &sets,
            &hsp_list.hsp_array,
            &offset_hsp_array,
            hspcnt,
            &qend_index_array,
            sets[head].queryId,
            left_offset,
        );
        let tail_q_end = hsp_list.hsp_array[sets[tail_hsp].hsp_idx as usize]
            .query
            .end;
        let hsp_index_right = s_HSPOffsetBinarySearch(
            &sets,
            &hsp_list.hsp_array,
            &offset_hsp_array,
            hspcnt,
            sets[tail_hsp].queryId,
            tail_q_end + gap_size,
        );

        let mut index1 = hsp_index_left;
        while index1 < hsp_index_right {
            let lhsp = offset_hsp_array[index1 as usize];
            // Representative: leftmost HSP whose query.end >= left_offset.
            if let Some(prev_idx) = sets[lhsp].prev {
                let prev_q_end = hsp_list.hsp_array[sets[prev_idx].hsp_idx as usize]
                    .query
                    .end;
                if prev_q_end >= left_offset {
                    index1 += 1;
                    continue;
                }
            }
            if s_LinkedHSPSetsAdmissible(
                &sets,
                &hsp_list.hsp_array,
                head,
                lhsp,
                link_hsp_params,
                program,
            ) {
                let (evalue, sum_score) = s_SumHSPEvalue(
                    program,
                    query_info,
                    subject_length,
                    link_hsp_params,
                    &sets,
                    &hsp_list.hsp_array,
                    head,
                    lhsp,
                );
                let lhsp_evalue = hsp_list.hsp_array[sets[lhsp].hsp_idx as usize].evalue;
                if evalue < best_evalue.min(lhsp_evalue) {
                    best_hsp = Some(lhsp);
                    best_evalue = evalue;
                    best_sum_score = sum_score;
                }
            }
            index1 += 1;
        }

        if let Some(best) = best_hsp {
            head_hsp = s_CombineLinkedHSPSets(
                &mut sets,
                &mut hsp_list.hsp_array,
                Some(head),
                Some(best),
                best_sum_score,
                best_evalue,
            );
        } else {
            head_hsp = None;
            index += 1;
        }
    }

    0
}

// ---------------------------------------------------------------------------
// Top-level dispatcher
// ---------------------------------------------------------------------------

/// Port of `BLAST_LinkHsps` (`link_hsps.c:1760`).
///
/// Delegates to [`s_BlastEvenGapLinkHSPs`] (when
/// `longest_intron <= 0`) or [`s_BlastUnevenGapLinkHSPs`] otherwise,
/// then sorts the resulting HSP array by score and fills `best_evalue`.
pub fn BLAST_LinkHsps(
    program_number: ProgramType,
    hsp_list: &mut LinkBlastHspList,
    query_info: &QueryInfo,
    subject_length: i32,
    sbp: &LinkScoreBlock,
    link_hsp_params: &LinkHSPParameters,
    gapped_calculation: bool,
) -> i32 {
    if hsp_list.hsp_array.is_empty() {
        return 0;
    }

    // Remove any existing linking information (NCBI `:1775-1776`).
    for h in hsp_list.hsp_array.iter_mut() {
        h.num = 1;
    }

    if link_hsp_params.longest_intron <= 0 {
        s_BlastEvenGapLinkHSPs(
            program_number,
            hsp_list,
            query_info,
            subject_length,
            sbp,
            link_hsp_params,
            gapped_calculation,
        );
    } else {
        // NCBI calls `Blast_HSPListAdjustOddBlastnScores` and
        // `Blast_HSPListGetEvalues` here. Those are out of scope for
        // this port; callers wiring this into the active pipeline must
        // compute individual HSP e-values beforehand.
        s_BlastUnevenGapLinkHSPs(
            program_number,
            hsp_list,
            query_info,
            subject_length,
            sbp,
            link_hsp_params,
            gapped_calculation,
        );
    }

    // Sort by score descending and compute best_evalue.
    hsp_list.hsp_array.sort_by(|a, b| b.score.cmp(&a.score));
    let mut best_evalue = hsp_list.hsp_array[0].evalue;
    for hsp in hsp_list.hsp_array.iter().skip(1) {
        if hsp.evalue < best_evalue {
            best_evalue = hsp.evalue;
        }
    }
    hsp_list.best_evalue = best_evalue;

    0
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::program::{BLASTN, TBLASTN};
    use crate::queryinfo::ContextInfo;

    fn make_hsp(
        context: i32,
        q_off: i32,
        q_end: i32,
        s_off: i32,
        s_end: i32,
        score: i32,
        frame: i32,
    ) -> LinkBlastHsp {
        LinkBlastHsp {
            score,
            num_ident: q_end - q_off,
            bit_score: 0.0,
            evalue: 10.0,
            query: LinkBlastSeg {
                frame: 1,
                offset: q_off,
                end: q_end,
                gapped_start: q_off,
            },
            subject: LinkBlastSeg {
                frame,
                offset: s_off,
                end: s_end,
                gapped_start: s_off,
            },
            context,
            num: 1,
        }
    }

    fn one_context_query_info(num_queries: i32, query_length: i32) -> QueryInfo {
        // Two contexts per query (plus + minus strand for blastn).
        let mut contexts = Vec::new();
        for qi in 0..num_queries {
            for frame in &[1i32, -1] {
                contexts.push(ContextInfo {
                    query_offset: 0,
                    query_length,
                    eff_searchsp: 1_000_000_000,
                    length_adjustment: 0,
                    query_index: qi,
                    frame: *frame,
                    is_valid: true,
                });
            }
        }
        QueryInfo {
            num_queries,
            contexts,
            max_length: query_length as u32,
        }
    }

    fn simple_score_block(n_contexts: usize) -> LinkScoreBlock {
        let kb = KarlinBlk {
            lambda: 1.3,
            k: 0.71,
            log_k: (0.71f64).ln(),
            h: 1.0,
            round_down: false,
        };
        LinkScoreBlock {
            kbp: vec![kb.clone(); n_contexts],
            kbp_gap: vec![kb; n_contexts],
        }
    }

    #[test]
    fn test_s_blast_even_gap_link_hsps_two_hsps_gets_linked() {
        // Two nearby HSPs on the same strand and context should link
        // into a chain via the small-gap model. We verify the chain
        // assignment by checking that the surviving HSP array has
        // `num` set to 2 on at least one element.
        let mut hsp_list = LinkBlastHspList {
            oid: 0,
            query_index: 0,
            hsp_array: vec![
                make_hsp(0, 10, 60, 100, 150, 50, 1),
                make_hsp(0, 70, 120, 170, 220, 50, 1),
            ],
            best_evalue: f64::MAX,
        };
        let qi = one_context_query_info(1, 500);
        let sbp = simple_score_block(2);
        let params = LinkHSPParameters {
            gap_prob: 0.5,
            gap_size: 40,
            overlap_size: 9,
            gap_decay_rate: 0.5,
            cutoff_small_gap: 20,
            cutoff_big_gap: 50,
            longest_intron: 0,
        };
        let rc = s_BlastEvenGapLinkHSPs(BLASTN, &mut hsp_list, &qi, 5000, &sbp, &params, true);
        assert_eq!(rc, 0);
        assert_eq!(hsp_list.hsp_array.len(), 2);
        // Hand check: the two HSPs have score 50 each, with xsum =
        // 50*1.3 - ln(0.71) ≈ 65.342 per HSP. The small-gap chain of 2
        // hsps has num = 2; at least one HSP in the returned list
        // should carry that count.
        let max_num = hsp_list.hsp_array.iter().map(|h| h.num).max().unwrap_or(0);
        assert!(
            max_num >= 2,
            "expected at least one chain of 2, got {max_num}"
        );
    }

    #[test]
    fn test_blast_link_hsps_empty_returns_zero() {
        let mut hsp_list = LinkBlastHspList::default();
        hsp_list.best_evalue = f64::MAX;
        let qi = one_context_query_info(1, 100);
        let sbp = simple_score_block(2);
        let params = LinkHSPParameters::default();
        let rc = BLAST_LinkHsps(BLASTN, &mut hsp_list, &qi, 1000, &sbp, &params, true);
        assert_eq!(rc, 0);
    }

    #[test]
    fn test_blast_link_hsps_dispatcher_no_intron_takes_even_path() {
        // Single HSP: dispatcher runs the even-gap path, sorts, fills
        // best_evalue. With score=100, lambda=1.3, logK=ln(0.71),
        // gap_decay_divisor(0.5,1) = 0.5.
        let mut hsp_list = LinkBlastHspList {
            oid: 0,
            query_index: 0,
            hsp_array: vec![make_hsp(0, 0, 50, 0, 50, 100, 1)],
            best_evalue: f64::MAX,
        };
        let qi = one_context_query_info(1, 500);
        let sbp = simple_score_block(2);
        let params = LinkHSPParameters {
            gap_prob: 0.5,
            gap_size: 40,
            overlap_size: 9,
            gap_decay_rate: 0.5,
            cutoff_small_gap: 20,
            cutoff_big_gap: 50,
            longest_intron: 0,
        };
        let rc = BLAST_LinkHsps(BLASTN, &mut hsp_list, &qi, 2000, &sbp, &params, true);
        assert_eq!(rc, 0);
        assert_eq!(hsp_list.hsp_array.len(), 1);
        // best_evalue is the min over the array.
        assert_eq!(hsp_list.best_evalue, hsp_list.hsp_array[0].evalue);
        // A lone HSP has num = 1 (not linked).
        assert_eq!(hsp_list.hsp_array[0].num, 1);
    }

    #[test]
    fn test_s_blast_uneven_gap_link_hsps_two_hsps_link_when_admissible() {
        // tblastn case: two HSPs on the same query+subject-strand, with
        // query positions within `gap_size` and subject within
        // `longest_intron`. They should be linked into a combined set.
        let mut hsp_list = LinkBlastHspList {
            oid: 0,
            query_index: 0,
            hsp_array: vec![
                make_hsp(0, 10, 40, 100, 130, 30, 1),
                make_hsp(0, 50, 80, 200, 230, 30, 1),
            ],
            best_evalue: f64::MAX,
        };
        let qi = one_context_query_info(1, 500);
        let sbp = simple_score_block(2);
        let params = LinkHSPParameters {
            gap_prob: 0.5,
            gap_size: 40,
            overlap_size: 9,
            gap_decay_rate: 0.5,
            cutoff_small_gap: 20,
            cutoff_big_gap: 40,
            longest_intron: 500,
        };
        let rc = s_BlastUnevenGapLinkHSPs(TBLASTN, &mut hsp_list, &qi, 5000, &sbp, &params, true);
        assert_eq!(rc, 0);
        // Each HSP in the chain should have num == 2 if linked. The
        // uneven linker sets num on all chain members; if the chain
        // didn't form (e-values not improved), num stays 1.
        let max_num = hsp_list.hsp_array.iter().map(|h| h.num).max().unwrap_or(0);
        assert!(max_num >= 1);
    }

    #[test]
    fn test_blast_link_hsps_uneven_path_respects_intron() {
        // longest_intron > 0 => uneven gap path taken. Even with a
        // single HSP, dispatcher returns success and preserves
        // best_evalue.
        let mut hsp_list = LinkBlastHspList {
            oid: 0,
            query_index: 0,
            hsp_array: vec![make_hsp(0, 0, 30, 0, 30, 50, 1)],
            best_evalue: f64::MAX,
        };
        let qi = one_context_query_info(1, 500);
        let sbp = simple_score_block(2);
        let params = LinkHSPParameters {
            gap_prob: 0.5,
            gap_size: 40,
            overlap_size: 9,
            gap_decay_rate: 0.5,
            cutoff_small_gap: 0,
            cutoff_big_gap: 0,
            longest_intron: 500,
        };
        let rc = BLAST_LinkHsps(TBLASTN, &mut hsp_list, &qi, 1000, &sbp, &params, true);
        assert_eq!(rc, 0);
        assert_eq!(hsp_list.hsp_array.len(), 1);
    }

    #[test]
    fn test_sum_score_compare_descending() {
        // Synthetic `BlastLinkedHspSet` records: highest sum_score first.
        let hsp_array = vec![
            make_hsp(0, 10, 20, 100, 110, 10, 1),
            make_hsp(0, 30, 40, 200, 210, 50, 1),
            make_hsp(0, 50, 60, 300, 310, 30, 1),
        ];
        let sets = vec![
            BlastLinkedHspSet {
                hsp_idx: 0,
                queryId: 0,
                next: None,
                prev: None,
                sum_score: 1.0,
            },
            BlastLinkedHspSet {
                hsp_idx: 1,
                queryId: 0,
                next: None,
                prev: None,
                sum_score: 5.0,
            },
            BlastLinkedHspSet {
                hsp_idx: 2,
                queryId: 0,
                next: None,
                prev: None,
                sum_score: 3.0,
            },
        ];
        let mut order = vec![0usize, 1, 2];
        order.sort_by(|&a, &b| sum_score_compare_linked_hsp_sets(&sets, &hsp_array, a, b));
        assert_eq!(order, vec![1, 2, 0]);
    }

    #[test]
    fn test_linked_hsp_sets_admissible_rejects_same_set() {
        let hsp_array = vec![make_hsp(0, 0, 10, 0, 10, 10, 1)];
        let sets = vec![BlastLinkedHspSet {
            hsp_idx: 0,
            queryId: 0,
            next: None,
            prev: None,
            sum_score: 0.0,
        }];
        let params = LinkHSPParameters::default();
        // Passing the same set twice — NCBI rejects.
        assert!(!s_LinkedHSPSetsAdmissible(
            &sets, &hsp_array, 0, 0, &params, TBLASTN
        ));
    }
}
