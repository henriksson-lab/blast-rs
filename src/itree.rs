//! Rust equivalent of blast_itree.c — interval tree for HSP containment.
//! Used to efficiently check if a new HSP is contained within existing ones.

use crate::hspstream::Hsp;

/// An interval [start, end) on one axis.
#[derive(Debug, Clone, Copy)]
pub struct Interval {
    pub start: i32,
    pub end: i32,
}

impl Interval {
    pub fn new(start: i32, end: i32) -> Self {
        Interval { start, end }
    }

    pub fn contains(&self, other: &Interval) -> bool {
        self.start <= other.start && self.end >= other.end
    }

    pub fn overlaps(&self, other: &Interval) -> bool {
        self.start < other.end && other.start < self.end
    }

    pub fn length(&self) -> i32 {
        self.end - self.start
    }
}

/// A 2D interval (query range + subject range) representing an HSP.
#[derive(Debug, Clone)]
pub struct Interval2D {
    pub query: Interval,
    pub subject: Interval,
    pub score: i32,
    pub query_index: i32,
    pub subject_frame: i32,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IntervalNodeDirection {
    Left,
    Right,
    Neither,
}

#[derive(Debug, Clone)]
pub struct SIntervalNode {
    pub leftend: i32,
    pub rightend: i32,
    pub leftptr: i32,
    pub midptr: i32,
    pub rightptr: i32,
    pub hsp: Option<Interval2D>,
}

impl Default for SIntervalNode {
    fn default() -> Self {
        Self {
            leftend: 0,
            rightend: 0,
            leftptr: 0,
            midptr: 0,
            rightptr: 0,
            hsp: None,
        }
    }
}

/// C-shaped interval tree for HSP containment checks. This mirrors NCBI's
/// `BlastIntervalTree`: nodes live in an index-addressed pool; non-leaf nodes
/// split a query or subject region into left/right children plus a midpoint
/// list, and leaves own the translated HSP interval.
pub struct IntervalTree {
    pub nodes: Vec<SIntervalNode>,
    pub num_alloc: i32,
    pub num_used: i32,
    pub s_min: i32,
    pub s_max: i32,
    hsp_count: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IntervalDirection {
    Left,
    Right,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CommonEndpointHsp {
    Input,
    Tree,
}

impl IntervalTree {
    pub fn new(q_max: i32, s_max: i32) -> Self {
        blast_interval_tree_init(0, q_max, 0, s_max)
    }

    /// Check if a new HSP is contained within any existing interval.
    /// Returns true if the new interval should be rejected (contained).
    pub fn is_contained(&self, query: Interval, subject: Interval, score: i32) -> bool {
        self.is_contained_with_metadata(query, subject, score, 0, 0)
    }

    /// Check containment with NCBI `s_HSPIsContained` metadata gates
    /// (`blast_itree.c:810`): the containing HSP must have the same query
    /// index, the same subject-frame sign, and score >= the candidate.
    pub fn is_contained_with_metadata(
        &self,
        query: Interval,
        subject: Interval,
        score: i32,
        query_index: i32,
        subject_frame: i32,
    ) -> bool {
        self.is_contained_with_metadata_and_min_diag_separation(
            query,
            subject,
            score,
            query_index,
            subject_frame,
            0,
        )
    }

    /// Check if a new HSP is contained within any existing interval, using
    /// NCBI megablast's diagonal-separation rule when nonzero.
    pub fn is_contained_with_min_diag_separation(
        &self,
        query: Interval,
        subject: Interval,
        score: i32,
        min_diag_separation: i32,
    ) -> bool {
        self.is_contained_with_metadata_and_min_diag_separation(
            query,
            subject,
            score,
            0,
            0,
            min_diag_separation,
        )
    }

    pub fn is_contained_with_metadata_and_min_diag_separation(
        &self,
        query: Interval,
        subject: Interval,
        score: i32,
        query_index: i32,
        subject_frame: i32,
        min_diag_separation: i32,
    ) -> bool {
        let hsp = Interval2D {
            query,
            subject,
            score,
            query_index,
            subject_frame,
        };
        blast_interval_tree_contains_interval(self, &hsp, min_diag_separation)
    }

    /// Add an interval to the tree, mirroring NCBI's `BlastIntervalTreeAddHSP`
    /// (`blast_itree.c:511`). Before inserting, NCBI checks whether an
    /// existing tree HSP shares the LEFT endpoint `(q.start, s.start)` or the
    /// RIGHT endpoint `(q.end, s.end)` via `s_IntervalTreeHasHSPEndpoint` →
    /// `s_HSPsHaveCommonEndpoint` (`blast_itree.c:242`). When such a sharing
    /// exists:
    ///  - if the tree HSP is higher-scoring, the new HSP is NOT inserted
    ///    (return early);
    ///  - if the new HSP is higher-scoring, the existing tree HSP is REMOVED
    ///    and the new HSP is inserted;
    ///  - on equal scores, the SHORTER HSP wins (by query length, then by
    ///    subject length); ties favor the HSP already in the tree.
    ///
    /// Without this dedup our tree accumulates multiple HSPs sharing
    /// endpoints, and the larger of them envelops legitimate seeds that
    /// NCBI's tree wouldn't envelop.
    pub fn insert(&mut self, query: Interval, subject: Interval, score: i32) {
        self.insert_with_metadata(query, subject, score, 0, 0)
    }

    pub fn insert_with_metadata(
        &mut self,
        query: Interval,
        subject: Interval,
        score: i32,
        query_index: i32,
        subject_frame: i32,
    ) {
        for which_end in [IntervalDirection::Left, IntervalDirection::Right] {
            let candidate = Interval2D {
                query,
                subject,
                score,
                query_index,
                subject_frame,
            };
            if s_interval_tree_has_hsp_endpoint_for_interval(self, &candidate, which_end) {
                return;
            }
        }
        let candidate = Interval2D {
            query,
            subject,
            score,
            query_index,
            subject_frame,
        };
        let _ = blast_interval_tree_add_interval(self, candidate, true);
    }

    pub fn len(&self) -> usize {
        self.hsp_count
    }

    pub fn reset(&mut self) {
        blast_interval_tree_reset(self);
    }
}

/// 1-1 translation of `Blast_IntervalTreeInit` (`blast_itree.c:145`) for the
/// Rust interval storage. The C tree records query and subject bounds in the
/// root node; the current Rust representation keeps only inserted intervals,
/// so the bounds are accepted for call-site parity.
pub fn blast_interval_tree_init(
    q_start: i32,
    q_end: i32,
    s_start: i32,
    s_end: i32,
) -> IntervalTree {
    let mut tree = IntervalTree {
        nodes: Vec::with_capacity(100),
        num_alloc: 100,
        num_used: 0,
        s_min: s_start,
        s_max: s_end,
        hsp_count: 0,
    };
    let mut retval = 0;
    let _ = s_interval_root_node_init(&mut tree, q_start, q_end, &mut retval);
    tree
}

/// 1-1 translation of `Blast_IntervalTreeFree` (`blast_itree.c:180`). Rust
/// ownership drops the tree; returning `None` matches the C function's NULL
/// return value.
pub fn blast_interval_tree_free(tree: &mut Option<IntervalTree>) -> Option<IntervalTree> {
    *tree = None;
    None
}

/// 1-1 translation of `Blast_IntervalTreeReset` (`blast_itree.c:193`).
pub fn blast_interval_tree_reset(tree: &mut IntervalTree) {
    tree.num_used = 1;
    tree.hsp_count = 0;
    tree.nodes.truncate(1);
    if let Some(root) = tree.nodes.get_mut(0) {
        root.leftptr = 0;
        root.midptr = 0;
        root.rightptr = 0;
        root.hsp = None;
    }
}

/// 1-1 translation of NCBI `s_IntervalNodeInit` (`blast_itree.c:58`) over
/// Rust's owned node pool.
pub fn s_interval_node_init(
    tree: &mut IntervalTree,
    parent_index: i32,
    dir: IntervalNodeDirection,
    ret_status: &mut i16,
) -> i32 {
    *ret_status = 0;
    if tree.num_used == tree.num_alloc {
        tree.num_alloc *= 2;
        let additional = (tree.num_alloc as usize).saturating_sub(tree.nodes.capacity());
        tree.nodes.reserve(additional);
    }

    let new_index = tree.num_used;
    tree.num_used += 1;
    tree.nodes.push(SIntervalNode::default());

    if dir == IntervalNodeDirection::Neither {
        return new_index;
    }

    let parent = tree.nodes[parent_index as usize].clone();
    let midpt = ((parent.leftend as i64 + parent.rightend as i64) / 2) as i32;
    let new_node = &mut tree.nodes[new_index as usize];
    new_node.leftptr = 0;
    new_node.midptr = 0;
    new_node.rightptr = 0;
    new_node.hsp = None;
    match dir {
        IntervalNodeDirection::Left => {
            new_node.leftend = parent.leftend;
            new_node.rightend = midpt;
        }
        IntervalNodeDirection::Right => {
            new_node.leftend = midpt + 1;
            new_node.rightend = parent.rightend;
        }
        IntervalNodeDirection::Neither => {}
    }
    new_index
}

/// 1-1 translation of NCBI `s_IntervalRootNodeInit` (`blast_itree.c:124`).
pub fn s_interval_root_node_init(
    tree: &mut IntervalTree,
    region_start: i32,
    region_end: i32,
    retval: &mut i16,
) -> i32 {
    let new_index = s_interval_node_init(tree, 0, IntervalNodeDirection::Neither, retval);
    if *retval != 0 {
        return 0;
    }
    let new_node = &mut tree.nodes[new_index as usize];
    new_node.leftptr = 0;
    new_node.midptr = 0;
    new_node.rightptr = 0;
    new_node.hsp = None;
    new_node.leftend = region_start;
    new_node.rightend = region_end;
    new_index
}

/// 1-1 translation of `BlastIntervalTreeAddHSP` (`blast_itree.c:511`) for the
/// local `Hsp` representation. `query_index` is the caller's translated
/// `s_GetQueryStrandOffset`/query-context key; the HSP's subject frame
/// preserves the NCBI same-strand metadata gate.
pub fn blast_interval_tree_add_hsp(hsp: &Hsp, tree: &mut IntervalTree, query_index: i32) -> i32 {
    let query = hsp_query_interval(hsp);
    let subject = hsp_subject_interval(hsp);
    debug_assert!(query.start <= query.end);
    debug_assert!(subject.start <= subject.end);

    let index_subject_range = tree.s_max > tree.s_min;
    if index_subject_range {
        debug_assert!(subject.start >= tree.s_min);
        debug_assert!(subject.end <= tree.s_max);
    }

    let candidate = Interval2D {
        query,
        subject,
        score: hsp.score,
        query_index,
        subject_frame: hsp.subject_frame,
    };

    if index_subject_range
        && s_interval_tree_has_hsp_endpoint_for_interval(tree, &candidate, IntervalDirection::Left)
    {
        return 0;
    }
    if index_subject_range
        && s_interval_tree_has_hsp_endpoint_for_interval(tree, &candidate, IntervalDirection::Right)
    {
        return 0;
    }

    let hsp = candidate;
    let index_subject = index_subject_range;
    let mut retval = 0;
    let mut root_index = 0i32;
    let new_index = s_interval_node_init(tree, 0, IntervalNodeDirection::Neither, &mut retval);
    if retval != 0 {
        return retval as i32;
    }
    tree.nodes[new_index as usize].leftptr = hsp.query_index;
    tree.nodes[new_index as usize].midptr = 0;
    tree.nodes[new_index as usize].hsp = Some(hsp.clone());
    tree.hsp_count += 1;

    let mut region_start = hsp.query.start;
    let mut region_end = hsp.query.end;
    let mut index_subject_range = false;

    loop {
        let middle = ((tree.nodes[root_index as usize].leftend as i64
            + tree.nodes[root_index as usize].rightend as i64)
            / 2) as i32;

        let (old_index, which_half) = if region_end < middle {
            if tree.nodes[root_index as usize].leftptr == 0 {
                tree.nodes[root_index as usize].leftptr = new_index;
                return 0;
            }
            let old_index = tree.nodes[root_index as usize].leftptr;
            if tree.nodes[old_index as usize].hsp.is_none() {
                root_index = old_index;
                continue;
            }
            (old_index, IntervalNodeDirection::Left)
        } else if region_start > middle {
            if tree.nodes[root_index as usize].rightptr == 0 {
                tree.nodes[root_index as usize].rightptr = new_index;
                return 0;
            }
            let old_index = tree.nodes[root_index as usize].rightptr;
            if tree.nodes[old_index as usize].hsp.is_none() {
                root_index = old_index;
                continue;
            }
            (old_index, IntervalNodeDirection::Right)
        } else if index_subject_range || !index_subject {
            tree.nodes[new_index as usize].midptr = tree.nodes[root_index as usize].midptr;
            tree.nodes[root_index as usize].midptr = new_index;
            return 0;
        } else {
            index_subject_range = true;
            if tree.nodes[root_index as usize].midptr == 0 {
                let mid_index =
                    s_interval_root_node_init(tree, tree.s_min, tree.s_max, &mut retval);
                if retval != 0 {
                    return retval as i32;
                }
                tree.nodes[root_index as usize].midptr = mid_index;
            }
            root_index = tree.nodes[root_index as usize].midptr;
            region_start = hsp.subject.start;
            region_end = hsp.subject.end;
            continue;
        };

        let mid_index = s_interval_node_init(tree, root_index, which_half, &mut retval);
        if retval != 0 {
            return retval as i32;
        }
        match which_half {
            IntervalNodeDirection::Left => tree.nodes[root_index as usize].leftptr = mid_index,
            IntervalNodeDirection::Right => tree.nodes[root_index as usize].rightptr = mid_index,
            IntervalNodeDirection::Neither => {}
        }

        let old_hsp = tree.nodes[old_index as usize]
            .hsp
            .as_ref()
            .expect("old leaf")
            .clone();
        let (old_region_start, old_region_end) = if index_subject_range {
            (old_hsp.subject.start, old_hsp.subject.end)
        } else {
            (old_hsp.query.start, old_hsp.query.end)
        };

        root_index = mid_index;
        let middle = ((tree.nodes[root_index as usize].leftend as i64
            + tree.nodes[root_index as usize].rightend as i64)
            / 2) as i32;
        if old_region_end < middle {
            tree.nodes[mid_index as usize].leftptr = old_index;
        } else if old_region_start > middle {
            tree.nodes[mid_index as usize].rightptr = old_index;
        } else if index_subject_range || !index_subject {
            tree.nodes[mid_index as usize].midptr = old_index;
        } else {
            let mid_index2 = s_interval_root_node_init(tree, tree.s_min, tree.s_max, &mut retval);
            if retval != 0 {
                return retval as i32;
            }
            tree.nodes[mid_index as usize].midptr = mid_index2;
            let old_region_start = old_hsp.subject.start;
            let old_region_end = old_hsp.subject.end;
            let middle = ((tree.nodes[mid_index2 as usize].leftend as i64
                + tree.nodes[mid_index2 as usize].rightend as i64)
                / 2) as i32;
            if old_region_end < middle {
                tree.nodes[mid_index2 as usize].leftptr = old_index;
            } else if old_region_start > middle {
                tree.nodes[mid_index2 as usize].rightptr = old_index;
            } else {
                tree.nodes[mid_index2 as usize].midptr = old_index;
            }
        }
    }
}

/// 1-1 translation of `BlastIntervalTreeContainsHSP` (`blast_itree.c:931`).
pub fn blast_interval_tree_contains_hsp(
    tree: &IntervalTree,
    hsp: &Hsp,
    query_index: i32,
    min_diag_separation: i32,
) -> bool {
    let candidate = Interval2D {
        query: hsp_query_interval(hsp),
        subject: hsp_subject_interval(hsp),
        score: hsp.score,
        query_index,
        subject_frame: hsp.subject_frame,
    };
    blast_interval_tree_contains_interval(tree, &candidate, min_diag_separation)
}

/// 1-1 translation of `BlastIntervalTreeNumRedundant` (`blast_itree.c:1196`).
/// NCBI counts only higher-scoring HSPs from the same query strand whose query
/// range contains the candidate range.
pub fn blast_interval_tree_num_redundant(tree: &IntervalTree, hsp: &Hsp, query_index: i32) -> i32 {
    let candidate = Interval2D {
        query: hsp_query_interval(hsp),
        subject: hsp_subject_interval(hsp),
        score: hsp.score,
        query_index,
        subject_frame: hsp.subject_frame,
    };
    blast_interval_tree_num_redundant_interval(tree, &candidate)
}

/// 1-1 translation of `BlastIntervalTreeMasksHSP` (`blast_itree.c:1107`) for
/// masklevel query-domain filtering over the node-pool interval tree.
pub fn blast_interval_tree_masks_hsp(
    tree: &IntervalTree,
    hsp: &Hsp,
    query_index: i32,
    masklevel: i32,
) -> bool {
    let candidate = Interval2D {
        query: hsp_query_interval(hsp),
        subject: hsp_subject_interval(hsp),
        score: hsp.score,
        query_index,
        subject_frame: hsp.subject_frame,
    };
    blast_interval_tree_masks_interval(tree, &candidate, 0, masklevel)
}

/// 1-1 translation of `s_HSPsHaveCommonEndpoint` (`blast_itree.c:274`) for
/// the flat Rust interval entry. The return value names which HSP C would
/// keep: input, tree, or no common endpoint.
pub fn s_hsps_have_common_endpoint(
    in_query: Interval,
    in_subject: Interval,
    in_score: i32,
    in_q_start: i32,
    in_subject_frame: i32,
    tree_hsp: &Interval2D,
    which_end: IntervalDirection,
) -> Option<CommonEndpointHsp> {
    if in_q_start != tree_hsp.query_index {
        return None;
    }
    if !same_subject_frame_sign(in_subject_frame, tree_hsp.subject_frame) {
        return None;
    }

    let same_endpoint = match which_end {
        IntervalDirection::Left => {
            in_query.start == tree_hsp.query.start && in_subject.start == tree_hsp.subject.start
        }
        IntervalDirection::Right => {
            in_query.end == tree_hsp.query.end && in_subject.end == tree_hsp.subject.end
        }
    };
    if !same_endpoint {
        return None;
    }

    if in_score > tree_hsp.score {
        return Some(CommonEndpointHsp::Input);
    }
    if in_score < tree_hsp.score {
        return Some(CommonEndpointHsp::Tree);
    }

    let in_q_length = in_query.end - in_query.start;
    let tree_q_length = tree_hsp.query.end - tree_hsp.query.start;
    if in_q_length > tree_q_length {
        return Some(CommonEndpointHsp::Tree);
    }
    if in_q_length < tree_q_length {
        return Some(CommonEndpointHsp::Input);
    }

    let in_s_length = in_subject.end - in_subject.start;
    let tree_s_length = tree_hsp.subject.end - tree_hsp.subject.start;
    if in_s_length > tree_s_length {
        return Some(CommonEndpointHsp::Tree);
    }
    if in_s_length < tree_s_length {
        return Some(CommonEndpointHsp::Input);
    }

    Some(CommonEndpointHsp::Tree)
}

fn s_interval_tree_has_hsp_endpoint_for_interval(
    tree: &mut IntervalTree,
    hsp: &Interval2D,
    which_end: IntervalDirection,
) -> bool {
    let mut root_index = 0usize;
    let target_offset = match which_end {
        IntervalDirection::Left => hsp.query.start,
        IntervalDirection::Right => hsp.query.end,
    };

    loop {
        let midptr = tree.nodes[root_index].midptr;
        if midptr != 0
            && s_midpoint_tree_has_hsp_endpoint_interval(
                tree,
                midptr,
                hsp,
                hsp.query_index,
                which_end,
            )
        {
            return true;
        }

        let middle = ((tree.nodes[root_index].leftend as i64
            + tree.nodes[root_index].rightend as i64)
            / 2) as i32;
        let tmp_index = if target_offset < middle {
            tree.nodes[root_index].leftptr
        } else if target_offset > middle {
            tree.nodes[root_index].rightptr
        } else {
            0
        };

        if tmp_index == 0 {
            return false;
        }

        let next_index = tmp_index as usize;
        if tree.nodes[next_index].hsp.is_some() {
            let best = {
                let next = &tree.nodes[next_index];
                s_hsps_have_common_endpoint(
                    hsp.query,
                    hsp.subject,
                    hsp.score,
                    hsp.query_index,
                    hsp.subject_frame,
                    next.hsp.as_ref().expect("leaf checked above"),
                    which_end,
                )
            };
            match best {
                Some(CommonEndpointHsp::Tree) => return true,
                Some(CommonEndpointHsp::Input) => {
                    if target_offset < middle {
                        tree.nodes[root_index].leftptr = 0;
                    } else if target_offset > middle {
                        tree.nodes[root_index].rightptr = 0;
                    }
                    tree.hsp_count = tree.hsp_count.saturating_sub(1);
                    return false;
                }
                None => return false,
            }
        }

        root_index = next_index;
    }
}

/// 1-1 translation of NCBI `s_MidpointTreeHasHSPEndpoint`
/// (`blast_itree.c:319`) over the subject-offset midpoint subtree.
pub fn s_midpoint_tree_has_hsp_endpoint(
    tree: &mut IntervalTree,
    root_index: i32,
    hsp: &Hsp,
    query_index: i32,
    which_end: IntervalDirection,
) -> bool {
    let candidate = Interval2D {
        query: hsp_query_interval(hsp),
        subject: hsp_subject_interval(hsp),
        score: hsp.score,
        query_index,
        subject_frame: hsp.subject_frame,
    };
    s_midpoint_tree_has_hsp_endpoint_interval(tree, root_index, &candidate, query_index, which_end)
}

/// 1-1 translation of NCBI `s_IntervalTreeHasHSPEndpoint`
/// (`blast_itree.c:429`) over the query-offset tree.
pub fn s_interval_tree_has_hsp_endpoint(
    tree: &mut IntervalTree,
    hsp: &Hsp,
    query_index: i32,
    which_end: IntervalDirection,
) -> bool {
    let candidate = Interval2D {
        query: hsp_query_interval(hsp),
        subject: hsp_subject_interval(hsp),
        score: hsp.score,
        query_index,
        subject_frame: hsp.subject_frame,
    };
    s_interval_tree_has_hsp_endpoint_for_interval(tree, &candidate, which_end)
}

/// 1-1 translation of NCBI `s_MidpointTreeContainsHSP`
/// (`blast_itree.c:865`) over the subject-offset midpoint subtree.
pub fn s_midpoint_tree_contains_hsp(
    tree: &IntervalTree,
    root_index: i32,
    hsp: &Hsp,
    query_index: i32,
    min_diag_separation: i32,
) -> bool {
    let candidate = Interval2D {
        query: hsp_query_interval(hsp),
        subject: hsp_subject_interval(hsp),
        score: hsp.score,
        query_index,
        subject_frame: hsp.subject_frame,
    };
    s_midpoint_tree_contains_interval(tree, root_index, &candidate, min_diag_separation)
}

fn s_midpoint_tree_has_hsp_endpoint_interval(
    tree: &mut IntervalTree,
    root_index: i32,
    hsp: &Interval2D,
    in_q_start: i32,
    which_end: IntervalDirection,
) -> bool {
    let mut root_index = root_index as usize;
    let target_offset = match which_end {
        IntervalDirection::Left => hsp.subject.start,
        IntervalDirection::Right => hsp.subject.end,
    };

    loop {
        let mut prev_index = root_index;
        let mut tmp_index = tree.nodes[root_index].midptr;
        while tmp_index != 0 {
            let next_index = tmp_index as usize;
            let next_midptr = tree.nodes[next_index].midptr;
            let best = {
                let next = &tree.nodes[next_index];
                s_hsps_have_common_endpoint(
                    hsp.query,
                    hsp.subject,
                    hsp.score,
                    in_q_start,
                    hsp.subject_frame,
                    next.hsp.as_ref().expect("midpoint list leaf"),
                    which_end,
                )
            };

            match best {
                Some(CommonEndpointHsp::Tree) => return true,
                Some(CommonEndpointHsp::Input) => {
                    tree.nodes[prev_index].midptr = next_midptr;
                    tree.hsp_count = tree.hsp_count.saturating_sub(1);
                    tmp_index = next_midptr;
                    continue;
                }
                None => {}
            }

            prev_index = next_index;
            tmp_index = next_midptr;
        }

        let middle = ((tree.nodes[root_index].leftend as i64
            + tree.nodes[root_index].rightend as i64)
            / 2) as i32;
        let tmp_index = if target_offset < middle {
            tree.nodes[root_index].leftptr
        } else if target_offset > middle {
            tree.nodes[root_index].rightptr
        } else {
            0
        };

        if tmp_index == 0 {
            return false;
        }

        let next_index = tmp_index as usize;
        if tree.nodes[next_index].hsp.is_some() {
            let best = {
                let next = &tree.nodes[next_index];
                s_hsps_have_common_endpoint(
                    hsp.query,
                    hsp.subject,
                    hsp.score,
                    in_q_start,
                    hsp.subject_frame,
                    next.hsp.as_ref().expect("leaf checked above"),
                    which_end,
                )
            };
            match best {
                Some(CommonEndpointHsp::Tree) => return true,
                Some(CommonEndpointHsp::Input) => {
                    if target_offset < middle {
                        tree.nodes[root_index].leftptr = 0;
                    } else if target_offset > middle {
                        tree.nodes[root_index].rightptr = 0;
                    }
                    tree.hsp_count = tree.hsp_count.saturating_sub(1);
                    return false;
                }
                None => return false,
            }
        }

        root_index = next_index;
    }
}

fn blast_interval_tree_add_interval(
    tree: &mut IntervalTree,
    hsp: Interval2D,
    index_subject: bool,
) -> i32 {
    if index_subject {
        if s_interval_tree_has_hsp_endpoint_for_interval(tree, &hsp, IntervalDirection::Left)
            || s_interval_tree_has_hsp_endpoint_for_interval(tree, &hsp, IntervalDirection::Right)
        {
            return 0;
        }
    }

    blast_interval_tree_add_prechecked_interval(tree, hsp, index_subject)
}

fn blast_interval_tree_add_prechecked_interval(
    tree: &mut IntervalTree,
    hsp: Interval2D,
    index_subject: bool,
) -> i32 {
    let mut retval = 0;
    let mut root_index = 0i32;
    let new_index = s_interval_node_init(tree, 0, IntervalNodeDirection::Neither, &mut retval);
    if retval != 0 {
        return retval as i32;
    }
    tree.nodes[new_index as usize].leftptr = hsp.query_index;
    tree.nodes[new_index as usize].midptr = 0;
    tree.nodes[new_index as usize].hsp = Some(hsp.clone());
    tree.hsp_count += 1;

    let mut region_start = hsp.query.start;
    let mut region_end = hsp.query.end;
    let mut index_subject_range = false;

    loop {
        let middle = ((tree.nodes[root_index as usize].leftend as i64
            + tree.nodes[root_index as usize].rightend as i64)
            / 2) as i32;

        let (old_index, which_half) = if region_end < middle {
            if tree.nodes[root_index as usize].leftptr == 0 {
                tree.nodes[root_index as usize].leftptr = new_index;
                return 0;
            }
            let old_index = tree.nodes[root_index as usize].leftptr;
            if tree.nodes[old_index as usize].hsp.is_none() {
                root_index = old_index;
                continue;
            }
            (old_index, IntervalNodeDirection::Left)
        } else if region_start > middle {
            if tree.nodes[root_index as usize].rightptr == 0 {
                tree.nodes[root_index as usize].rightptr = new_index;
                return 0;
            }
            let old_index = tree.nodes[root_index as usize].rightptr;
            if tree.nodes[old_index as usize].hsp.is_none() {
                root_index = old_index;
                continue;
            }
            (old_index, IntervalNodeDirection::Right)
        } else if index_subject_range || !index_subject {
            tree.nodes[new_index as usize].midptr = tree.nodes[root_index as usize].midptr;
            tree.nodes[root_index as usize].midptr = new_index;
            return 0;
        } else {
            index_subject_range = true;
            if tree.nodes[root_index as usize].midptr == 0 {
                let mid_index =
                    s_interval_root_node_init(tree, tree.s_min, tree.s_max, &mut retval);
                if retval != 0 {
                    return retval as i32;
                }
                tree.nodes[root_index as usize].midptr = mid_index;
            }
            root_index = tree.nodes[root_index as usize].midptr;
            region_start = hsp.subject.start;
            region_end = hsp.subject.end;
            continue;
        };

        let mid_index = s_interval_node_init(tree, root_index, which_half, &mut retval);
        if retval != 0 {
            return retval as i32;
        }
        match which_half {
            IntervalNodeDirection::Left => tree.nodes[root_index as usize].leftptr = mid_index,
            IntervalNodeDirection::Right => tree.nodes[root_index as usize].rightptr = mid_index,
            IntervalNodeDirection::Neither => {}
        }

        let old_hsp = tree.nodes[old_index as usize]
            .hsp
            .as_ref()
            .expect("old leaf")
            .clone();
        let (old_region_start, old_region_end) = if index_subject_range {
            (old_hsp.subject.start, old_hsp.subject.end)
        } else {
            (old_hsp.query.start, old_hsp.query.end)
        };

        root_index = mid_index;
        let middle = ((tree.nodes[root_index as usize].leftend as i64
            + tree.nodes[root_index as usize].rightend as i64)
            / 2) as i32;
        if old_region_end < middle {
            tree.nodes[mid_index as usize].leftptr = old_index;
        } else if old_region_start > middle {
            tree.nodes[mid_index as usize].rightptr = old_index;
        } else if index_subject_range || !index_subject {
            tree.nodes[mid_index as usize].midptr = old_index;
        } else {
            let mid_index2 = s_interval_root_node_init(tree, tree.s_min, tree.s_max, &mut retval);
            if retval != 0 {
                return retval as i32;
            }
            tree.nodes[mid_index as usize].midptr = mid_index2;
            let old_region_start = old_hsp.subject.start;
            let old_region_end = old_hsp.subject.end;
            let middle = ((tree.nodes[mid_index2 as usize].leftend as i64
                + tree.nodes[mid_index2 as usize].rightend as i64)
                / 2) as i32;
            if old_region_end < middle {
                tree.nodes[mid_index2 as usize].leftptr = old_index;
            } else if old_region_start > middle {
                tree.nodes[mid_index2 as usize].rightptr = old_index;
            } else {
                tree.nodes[mid_index2 as usize].midptr = old_index;
            }
        }
    }
}

/// 1-1 translation of `s_HSPIsContained` (`blast_itree.c:810`).
pub fn s_hsp_is_contained(
    in_hsp: &Interval2D,
    in_q_start: i32,
    tree_hsp: &Interval2D,
    tree_q_start: i32,
    min_diag_separation: i32,
) -> bool {
    if in_q_start != tree_q_start {
        return false;
    }
    in_hsp.score <= tree_hsp.score
        && same_subject_frame_sign(in_hsp.subject_frame, tree_hsp.subject_frame)
        && tree_hsp.query.contains(&in_hsp.query)
        && tree_hsp.subject.contains(&in_hsp.subject)
        && intervals_are_close_on_diagonal(
            tree_hsp.query,
            tree_hsp.subject,
            in_hsp.query,
            in_hsp.subject,
            min_diag_separation,
        )
}

fn s_midpoint_tree_contains_interval(
    tree: &IntervalTree,
    root_index: i32,
    hsp: &Interval2D,
    min_diag_separation: i32,
) -> bool {
    let mut node_index = root_index as usize;
    loop {
        if tree.nodes[node_index].hsp.is_some() {
            return s_hsp_is_contained(
                hsp,
                hsp.query_index,
                tree.nodes[node_index].hsp.as_ref().expect("leaf"),
                tree.nodes[node_index].leftptr,
                min_diag_separation,
            );
        }

        let mut tmp_index = tree.nodes[node_index].midptr;
        while tmp_index != 0 {
            let tmp_node = &tree.nodes[tmp_index as usize];
            if s_hsp_is_contained(
                hsp,
                hsp.query_index,
                tmp_node.hsp.as_ref().expect("midpoint list leaf"),
                tmp_node.leftptr,
                min_diag_separation,
            ) {
                return true;
            }
            tmp_index = tmp_node.midptr;
        }

        let middle = ((tree.nodes[node_index].leftend as i64
            + tree.nodes[node_index].rightend as i64)
            / 2) as i32;
        let tmp_index = if hsp.subject.end < middle {
            tree.nodes[node_index].leftptr
        } else if hsp.subject.start > middle {
            tree.nodes[node_index].rightptr
        } else {
            0
        };
        if tmp_index == 0 {
            return false;
        }
        node_index = tmp_index as usize;
    }
}

fn blast_interval_tree_contains_interval(
    tree: &IntervalTree,
    hsp: &Interval2D,
    min_diag_separation: i32,
) -> bool {
    let mut node_index = 0usize;
    loop {
        if tree.nodes[node_index].hsp.is_some() {
            return s_hsp_is_contained(
                hsp,
                hsp.query_index,
                tree.nodes[node_index].hsp.as_ref().expect("leaf"),
                tree.nodes[node_index].leftptr,
                min_diag_separation,
            );
        }

        let midptr = tree.nodes[node_index].midptr;
        if midptr > 0 && s_midpoint_tree_contains_interval(tree, midptr, hsp, min_diag_separation) {
            return true;
        }

        let middle = ((tree.nodes[node_index].leftend as i64
            + tree.nodes[node_index].rightend as i64)
            / 2) as i32;
        let tmp_index = if hsp.query.end < middle {
            tree.nodes[node_index].leftptr
        } else if hsp.query.start > middle {
            tree.nodes[node_index].rightptr
        } else {
            0
        };
        if tmp_index == 0 {
            return false;
        }
        node_index = tmp_index as usize;
    }
}

fn blast_interval_tree_num_redundant_interval(tree: &IntervalTree, hsp: &Interval2D) -> i32 {
    if tree.nodes[0].midptr != 0 && tree.nodes[tree.nodes[0].midptr as usize].hsp.is_none() {
        return num_redundant_all_nodes(tree, hsp);
    }

    let mut node_index = 0usize;
    let mut num_redundant = 0;
    loop {
        if tree.nodes[node_index].hsp.is_some() {
            let node = &tree.nodes[node_index];
            return num_redundant
                + s_hsp_query_range_is_contained(
                    hsp.query,
                    hsp.score,
                    hsp.query_index,
                    node.hsp.as_ref().expect("leaf").query,
                    node.hsp.as_ref().expect("leaf").score,
                    node.leftptr,
                );
        }

        let mut tmp_index = tree.nodes[node_index].midptr;
        while tmp_index != 0 {
            let tmp_node = &tree.nodes[tmp_index as usize];
            let tree_hsp = tmp_node.hsp.as_ref().expect("midpoint list leaf");
            num_redundant += s_hsp_query_range_is_contained(
                hsp.query,
                hsp.score,
                hsp.query_index,
                tree_hsp.query,
                tree_hsp.score,
                tmp_node.leftptr,
            );
            tmp_index = tmp_node.midptr;
        }

        let middle = ((tree.nodes[node_index].leftend as i64
            + tree.nodes[node_index].rightend as i64)
            / 2) as i32;
        let tmp_index = if hsp.query.end < middle {
            tree.nodes[node_index].leftptr
        } else if hsp.query.start > middle {
            tree.nodes[node_index].rightptr
        } else {
            0
        };
        if tmp_index == 0 {
            return num_redundant;
        }
        node_index = tmp_index as usize;
    }
}

fn num_redundant_all_nodes(tree: &IntervalTree, hsp: &Interval2D) -> i32 {
    tree.nodes
        .iter()
        .filter_map(|node| node.hsp.as_ref().map(|tree_hsp| (node, tree_hsp)))
        .map(|(node, tree_hsp)| {
            s_hsp_query_range_is_contained(
                hsp.query,
                hsp.score,
                hsp.query_index,
                tree_hsp.query,
                tree_hsp.score,
                node.leftptr,
            )
        })
        .sum()
}

fn blast_interval_tree_masks_interval(
    tree: &IntervalTree,
    hsp: &Interval2D,
    subtree_index: i32,
    masklevel: i32,
) -> bool {
    let mut node_index = subtree_index as usize;
    loop {
        if tree.nodes[node_index].hsp.is_some() {
            let node = &tree.nodes[node_index];
            return s_hsp_query_range_is_masklevel_contained(
                hsp.query,
                hsp.score,
                hsp.query_index,
                node.hsp.as_ref().expect("leaf").query,
                node.hsp.as_ref().expect("leaf").score,
                node.leftptr,
                masklevel,
            ) != 0;
        }

        let mut tmp_index = tree.nodes[node_index].midptr;
        while tmp_index != 0 {
            let tmp_node = &tree.nodes[tmp_index as usize];
            if let Some(tree_hsp) = &tmp_node.hsp {
                if s_hsp_query_range_is_masklevel_contained(
                    hsp.query,
                    hsp.score,
                    hsp.query_index,
                    tree_hsp.query,
                    tree_hsp.score,
                    tmp_node.leftptr,
                    masklevel,
                ) != 0
                {
                    return true;
                }
                tmp_index = tmp_node.midptr;
            } else {
                if masklevel_subtree_any_leaf(tree, tmp_index, hsp, masklevel) {
                    return true;
                }
                break;
            }
        }

        let middle = ((tree.nodes[node_index].leftend as i64
            + tree.nodes[node_index].rightend as i64)
            / 2) as i32;
        let tmp_index = if hsp.query.end < middle {
            tree.nodes[node_index].leftptr
        } else if hsp.query.start > middle {
            tree.nodes[node_index].rightptr
        } else {
            if tree.nodes[node_index].leftptr != 0
                && blast_interval_tree_masks_interval(
                    tree,
                    hsp,
                    tree.nodes[node_index].leftptr,
                    masklevel,
                )
            {
                return true;
            }
            if tree.nodes[node_index].rightptr != 0
                && blast_interval_tree_masks_interval(
                    tree,
                    hsp,
                    tree.nodes[node_index].rightptr,
                    masklevel,
                )
            {
                return true;
            }
            0
        };
        if tmp_index == 0 {
            return false;
        }
        node_index = tmp_index as usize;
    }
}

fn masklevel_subtree_any_leaf(
    tree: &IntervalTree,
    root_index: i32,
    hsp: &Interval2D,
    masklevel: i32,
) -> bool {
    let mut stack = vec![root_index];
    while let Some(index) = stack.pop() {
        let node = &tree.nodes[index as usize];
        if let Some(tree_hsp) = &node.hsp {
            if s_hsp_query_range_is_masklevel_contained(
                hsp.query,
                hsp.score,
                hsp.query_index,
                tree_hsp.query,
                tree_hsp.score,
                node.leftptr,
                masklevel,
            ) != 0
            {
                return true;
            }
        } else {
            if node.leftptr != 0 {
                stack.push(node.leftptr);
            }
            if node.midptr != 0 {
                stack.push(node.midptr);
            }
            if node.rightptr != 0 {
                stack.push(node.rightptr);
            }
        }
    }
    false
}

/// 1-1 translation of `s_HSPQueryRangeIsContained` (`blast_itree.c:1042`).
pub fn s_hsp_query_range_is_contained(
    in_query: Interval,
    in_score: i32,
    in_q_start: i32,
    tree_query: Interval,
    tree_score: i32,
    tree_q_start: i32,
) -> i32 {
    if in_q_start != tree_q_start || in_score >= tree_score {
        return 0;
    }
    if tree_query.start <= in_query.start && tree_query.end >= in_query.end {
        return 1;
    }
    0
}

/// 1-1 translation of `s_HSPQueryRangeIsMasklevelContained`
/// (`blast_itree.c:1086`) after caller-side query strand normalization.
pub fn s_hsp_query_range_is_masklevel_contained(
    in_query: Interval,
    in_score: i32,
    in_q_start: i32,
    tree_query: Interval,
    tree_score: i32,
    tree_q_start: i32,
    masklevel: i32,
) -> i32 {
    if in_q_start != tree_q_start || in_score > tree_score {
        return 0;
    }
    let in_len = in_query.end - in_query.start;
    if in_len == 0 {
        return 0;
    }

    let overlap_start = tree_query.start.max(in_query.start);
    let overlap_end = tree_query.end.min(in_query.end);
    let perc_overlap = ((100.0 * (overlap_end - overlap_start) as f64) / in_len as f64) as i32;
    if perc_overlap >= masklevel {
        return 1;
    }
    0
}

fn hsp_query_interval(hsp: &Hsp) -> Interval {
    Interval::new(hsp.query_offset, hsp.query_end)
}

fn hsp_subject_interval(hsp: &Hsp) -> Interval {
    Interval::new(hsp.subject_offset, hsp.subject_end)
}

fn same_subject_frame_sign(a: i32, b: i32) -> bool {
    (a == 0 && b == 0) || (a > 0 && b > 0) || (a < 0 && b < 0)
}

fn intervals_are_close_on_diagonal(
    existing_query: Interval,
    existing_subject: Interval,
    query: Interval,
    subject: Interval,
    min_diag_separation: i32,
) -> bool {
    if min_diag_separation == 0 {
        return true;
    }

    diagonal_distance(
        existing_query.start,
        existing_subject.start,
        query.start,
        subject.start,
    ) < min_diag_separation
        || diagonal_distance(
            existing_query.end,
            existing_subject.end,
            query.end,
            subject.end,
        ) < min_diag_separation
}

fn diagonal_distance(q1: i32, s1: i32, q2: i32, s2: i32) -> i32 {
    ((q1 - s1) - (q2 - s2)).abs()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_hsp(score: i32, q_start: i32, q_end: i32, s_start: i32, s_end: i32) -> Hsp {
        Hsp {
            score,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 0.0,
            query_offset: q_start,
            query_end: q_end,
            query_gapped_start: q_start,
            subject_offset: s_start,
            subject_end: s_end,
            subject_gapped_start: s_start,
            context: 0,
            query_frame: 0,
            subject_frame: 1,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        }
    }

    #[test]
    fn test_containment() {
        let mut tree = IntervalTree::new(100, 100);
        tree.insert(Interval::new(0, 50), Interval::new(0, 50), 100);

        assert!(tree.is_contained(Interval::new(10, 40), Interval::new(10, 40), 90));
        assert!(!tree.is_contained(Interval::new(10, 40), Interval::new(10, 40), 110));
        assert!(!tree.is_contained(Interval::new(0, 60), Interval::new(0, 50), 90));
    }

    #[test]
    fn test_megablast_diagonal_separation_keeps_distant_contained_hsp() {
        let mut tree = IntervalTree::new(100, 100);
        tree.insert(Interval::new(3, 40), Interval::new(3, 40), 37);

        assert!(!tree.is_contained_with_min_diag_separation(
            Interval::new(11, 40),
            Interval::new(3, 32),
            30,
            6
        ));
        assert!(tree.is_contained_with_min_diag_separation(
            Interval::new(5, 35),
            Interval::new(5, 35),
            30,
            6
        ));
    }

    #[test]
    fn test_containment_requires_same_query_index_and_subject_frame_sign() {
        let mut tree = IntervalTree::new(100, 100);
        tree.insert_with_metadata(Interval::new(0, 50), Interval::new(0, 50), 100, 0, 1);

        assert!(tree.is_contained_with_metadata(
            Interval::new(10, 40),
            Interval::new(10, 40),
            90,
            0,
            2
        ));
        assert!(!tree.is_contained_with_metadata(
            Interval::new(10, 40),
            Interval::new(10, 40),
            90,
            1,
            1
        ));
        assert!(!tree.is_contained_with_metadata(
            Interval::new(10, 40),
            Interval::new(10, 40),
            90,
            0,
            -1
        ));
    }

    #[test]
    fn translated_hsp_is_contained_matches_score_frame_query_and_diag_gates() {
        let tree_hsp = Interval2D {
            query: Interval::new(0, 50),
            subject: Interval::new(0, 50),
            score: 100,
            query_index: 0,
            subject_frame: 1,
        };
        let contained = Interval2D {
            query: Interval::new(10, 40),
            subject: Interval::new(10, 40),
            score: 90,
            query_index: 0,
            subject_frame: 2,
        };
        assert!(s_hsp_is_contained(&contained, 0, &tree_hsp, 0, 0));

        let mut too_good = contained.clone();
        too_good.score = 101;
        assert!(!s_hsp_is_contained(&too_good, 0, &tree_hsp, 0, 0));
        assert!(!s_hsp_is_contained(&contained, 1, &tree_hsp, 0, 0));

        let mut opposite_frame = contained.clone();
        opposite_frame.subject_frame = -1;
        assert!(!s_hsp_is_contained(&opposite_frame, 0, &tree_hsp, 0, 0));

        let shifted_diag = Interval2D {
            query: Interval::new(10, 40),
            subject: Interval::new(0, 30),
            score: 90,
            query_index: 0,
            subject_frame: 1,
        };
        assert!(s_hsp_is_contained(&shifted_diag, 0, &tree_hsp, 0, 0));
        assert!(!s_hsp_is_contained(&shifted_diag, 0, &tree_hsp, 0, 6));
    }

    #[test]
    fn test_common_endpoint_dedup_respects_query_index_and_frame_sign() {
        let mut tree = IntervalTree::new(100, 100);
        tree.insert_with_metadata(Interval::new(0, 50), Interval::new(0, 50), 100, 0, 1);
        tree.insert_with_metadata(Interval::new(0, 50), Interval::new(0, 50), 90, 1, 1);
        tree.insert_with_metadata(Interval::new(0, 50), Interval::new(0, 50), 80, 0, -1);

        assert_eq!(tree.len(), 3);
    }

    #[test]
    fn translated_common_endpoint_helper_uses_score_length_and_existing_tie() {
        let tree_hsp = Interval2D {
            query: Interval::new(0, 50),
            subject: Interval::new(10, 60),
            score: 100,
            query_index: 2,
            subject_frame: 1,
        };

        assert_eq!(
            s_hsps_have_common_endpoint(
                Interval::new(0, 40),
                Interval::new(10, 50),
                100,
                2,
                1,
                &tree_hsp,
                IntervalDirection::Left,
            ),
            Some(CommonEndpointHsp::Input)
        );
        assert_eq!(
            s_hsps_have_common_endpoint(
                Interval::new(0, 40),
                Interval::new(10, 50),
                90,
                2,
                1,
                &tree_hsp,
                IntervalDirection::Left,
            ),
            Some(CommonEndpointHsp::Tree)
        );
        assert_eq!(
            s_hsps_have_common_endpoint(
                tree_hsp.query,
                tree_hsp.subject,
                tree_hsp.score,
                2,
                1,
                &tree_hsp,
                IntervalDirection::Right,
            ),
            Some(CommonEndpointHsp::Tree)
        );
        assert_eq!(
            s_hsps_have_common_endpoint(
                tree_hsp.query,
                tree_hsp.subject,
                tree_hsp.score,
                3,
                1,
                &tree_hsp,
                IntervalDirection::Left,
            ),
            None
        );
    }

    #[test]
    fn translated_endpoint_tree_helpers_remove_worse_or_reject_input() {
        let mut tree = IntervalTree::new(100, 100);
        let existing = make_hsp(100, 0, 50, 10, 60);
        let stronger = make_hsp(110, 0, 40, 10, 50);
        let weaker = make_hsp(90, 0, 30, 10, 40);

        blast_interval_tree_add_hsp(&existing, &mut tree, 2);
        assert!(!s_interval_tree_has_hsp_endpoint(
            &mut tree,
            &stronger,
            2,
            IntervalDirection::Left
        ));
        assert_eq!(tree.len(), 0);

        blast_interval_tree_add_hsp(&existing, &mut tree, 2);
        let midpoint_root = tree.nodes[0].midptr;
        assert!(s_midpoint_tree_has_hsp_endpoint(
            &mut tree,
            midpoint_root,
            &weaker,
            2,
            IntervalDirection::Left
        ));
        assert_eq!(tree.len(), 1);
    }

    #[test]
    fn translated_midpoint_contains_hsp_matches_public_contains() {
        let mut tree = blast_interval_tree_init(0, 100, 0, 100);
        let outer = make_hsp(100, 0, 50, 0, 50);
        let inner = make_hsp(80, 10, 40, 10, 40);

        blast_interval_tree_add_hsp(&outer, &mut tree, 0);
        let midpoint_root = tree.nodes[0].midptr;

        assert!(s_midpoint_tree_contains_hsp(
            &tree,
            midpoint_root,
            &inner,
            0,
            0
        ));
        assert_eq!(
            s_midpoint_tree_contains_hsp(&tree, midpoint_root, &inner, 0, 6),
            blast_interval_tree_contains_hsp(&tree, &inner, 0, 6)
        );
    }

    #[test]
    fn translated_query_range_helpers_match_score_and_query_gates() {
        assert_eq!(
            s_hsp_query_range_is_contained(
                Interval::new(10, 30),
                80,
                0,
                Interval::new(0, 40),
                90,
                0,
            ),
            1
        );
        assert_eq!(
            s_hsp_query_range_is_contained(
                Interval::new(10, 30),
                90,
                0,
                Interval::new(0, 40),
                90,
                0,
            ),
            0
        );

        assert_eq!(
            s_hsp_query_range_is_masklevel_contained(
                Interval::new(20, 50),
                80,
                0,
                Interval::new(10, 40),
                80,
                0,
                60,
            ),
            1
        );
        assert_eq!(
            s_hsp_query_range_is_masklevel_contained(
                Interval::new(20, 50),
                90,
                0,
                Interval::new(10, 40),
                80,
                0,
                60,
            ),
            0
        );
    }

    #[test]
    fn translated_interval_tree_lifecycle_matches_c_shape() {
        let mut tree = Some(blast_interval_tree_init(0, 100, 0, 100));
        let hsp = make_hsp(50, 10, 30, 20, 40);
        assert_eq!(
            blast_interval_tree_add_hsp(&hsp, tree.as_mut().unwrap(), 0),
            0
        );
        assert_eq!(tree.as_ref().unwrap().len(), 1);

        blast_interval_tree_reset(tree.as_mut().unwrap());
        assert_eq!(tree.as_ref().unwrap().len(), 0);

        assert!(blast_interval_tree_free(&mut tree).is_none());
        assert!(tree.is_none());
    }

    #[test]
    fn translated_interval_tree_contains_hsp_uses_metadata_and_diagonal() {
        let mut tree = blast_interval_tree_init(0, 100, 0, 100);
        let outer = make_hsp(100, 0, 50, 0, 50);
        let inner = make_hsp(80, 10, 40, 10, 40);
        let distant_diag = make_hsp(80, 10, 40, 0, 30);

        blast_interval_tree_add_hsp(&outer, &mut tree, 2);

        assert!(blast_interval_tree_contains_hsp(&tree, &inner, 2, 0));
        assert!(!blast_interval_tree_contains_hsp(&tree, &inner, 3, 0));
        assert!(!blast_interval_tree_contains_hsp(
            &tree,
            &distant_diag,
            2,
            6
        ));
    }

    #[test]
    fn translated_interval_tree_counts_query_redundancy_only_for_higher_scores() {
        let mut tree = blast_interval_tree_init(0, 100, 0, 100);
        let stronger = make_hsp(100, 0, 50, 0, 50);
        let same_score = make_hsp(80, 0, 50, 10, 60);
        let candidate = make_hsp(80, 10, 40, 20, 50);

        blast_interval_tree_add_hsp(&stronger, &mut tree, 0);
        blast_interval_tree_add_hsp(&same_score, &mut tree, 0);

        assert_eq!(blast_interval_tree_num_redundant(&tree, &candidate, 0), 1);
        assert_eq!(blast_interval_tree_num_redundant(&tree, &candidate, 1), 0);
    }

    #[test]
    fn translated_interval_tree_masklevel_uses_query_overlap_percentage() {
        let mut tree = blast_interval_tree_init(0, 100, 0, 100);
        let existing = make_hsp(100, 10, 40, 0, 30);
        let candidate = make_hsp(80, 20, 50, 10, 40);
        let too_strong = make_hsp(120, 20, 50, 10, 40);

        blast_interval_tree_add_hsp(&existing, &mut tree, 0);

        assert!(blast_interval_tree_masks_hsp(&tree, &candidate, 0, 60));
        assert!(!blast_interval_tree_masks_hsp(&tree, &candidate, 0, 70));
        assert!(!blast_interval_tree_masks_hsp(&tree, &too_strong, 0, 60));
    }
}
