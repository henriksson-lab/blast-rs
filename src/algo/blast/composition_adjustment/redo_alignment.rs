use std::ffi::c_void;
use std::ptr;

use crate::algo::blast::composition_adjustment::compo_heap::{
    blast_compo_heap_filled_to_cutoff, blast_compo_heap_would_insert, BlastCompoHeap,
};
use crate::algo::blast::composition_adjustment::compo_mode_condition::MatrixAdjustRule;
use crate::algo::blast::composition_adjustment::composition_adjustment::{
    blast_adjust_scores, blast_get_composition_range, blast_read_aa_composition,
    BlastAminoAcidComposition, BlastCompositionWorkspace, BlastMatrixInfo, ECompoAdjustModes,
};
use crate::algo::blast::composition_adjustment::smith_waterman::{
    blast_forbidden_ranges_clear, blast_forbidden_ranges_push, blast_smith_waterman_find_start,
    blast_smith_waterman_find_start_pssm, blast_smith_waterman_score_only,
    blast_smith_waterman_score_only_pssm, BlastForbiddenRanges,
};
use crate::matrix::AA_SIZE;

pub const LOCAL_LN2: f64 = 0.69314718055994530941723212145818;
pub const EVALUE_STRETCH: i32 = 5;
pub const K_WINDOW_BORDER: i32 = 200;
pub const K_RE_MATRIX_ADJUSTMENT_PSEUDOCOUNTS: i32 = 20;
pub const KAPPA_BIT_TOL: f64 = 2.0;
pub const MINIMUM_LENGTH_NEAR_IDENTICAL: i32 = 50;
pub const COMPO_INTENSE_DEBUG: bool = false;

unsafe extern "C" {
    fn free(ptr: *mut c_void);
}

#[derive(Debug, Clone)]
pub struct BlastCompoAlignment {
    pub score: i32,
    pub matrix_adjust_rule: MatrixAdjustRule,
    pub query_index: i32,
    pub query_start: i32,
    pub query_end: i32,
    pub match_start: i32,
    pub match_end: i32,
    pub frame: i32,
    pub context: *mut c_void,
    pub next: Option<Box<BlastCompoAlignment>>,
}

#[allow(clippy::too_many_arguments)]
pub fn blast_compo_alignment_new(
    score: i32,
    matrix_adjust_rule: MatrixAdjustRule,
    query_start: i32,
    query_end: i32,
    query_index: i32,
    match_start: i32,
    match_end: i32,
    frame: i32,
    context: *mut c_void,
) -> Option<Box<BlastCompoAlignment>> {
    Some(Box::new(BlastCompoAlignment {
        score,
        matrix_adjust_rule,
        query_index,
        query_start,
        query_end,
        match_start,
        match_end,
        frame,
        context,
        next: None,
    }))
}

pub fn blast_compo_alignments_free(
    palign: &mut Option<Box<BlastCompoAlignment>>,
    free_context: Option<FreeAlignTracebackType>,
) {
    let mut align = palign.take();
    while let Some(mut current) = align {
        let align_next = current.next.take();
        if let Some(free_context) = free_context {
            if !current.context.is_null() {
                free_context(current.context);
                current.context = ptr::null_mut();
            }
        }
        align = align_next;
    }
}

pub fn s_alignments_rev(plist: &mut Option<Box<BlastCompoAlignment>>) {
    let mut list = plist.take();
    let mut new_list = None;
    while let Some(mut current) = list {
        let list_next = current.next.take();
        current.next = new_list;
        new_list = Some(current);
        list = list_next;
    }
    *plist = new_list;
}

pub fn s_alignment_cmp(a: &BlastCompoAlignment, b: &BlastCompoAlignment) -> i32 {
    let mut result = b.score.cmp(&a.score) as i32;
    if result == 0 {
        result = a.match_start.cmp(&b.match_start) as i32;
    }
    if result == 0 {
        result = b.match_end.cmp(&a.match_end) as i32;
    }
    if result == 0 {
        result = a.query_start.cmp(&b.query_start) as i32;
    }
    if result == 0 {
        result = b.query_end.cmp(&a.query_end) as i32;
    }
    result
}

pub fn s_alignments_are_sorted(alignments: Option<&BlastCompoAlignment>) -> i32 {
    let mut align = alignments;
    while let Some(current) = align {
        if current
            .next
            .as_deref()
            .is_some_and(|next| next.score > current.score)
        {
            return 0;
        }
        align = current.next.as_deref();
    }
    1
}

pub fn s_distinct_alignments_length(mut list: Option<&BlastCompoAlignment>) -> i32 {
    let mut length = 0;
    while let Some(current) = list {
        length += 1;
        list = current.next.as_deref();
    }
    length
}

pub fn s_distinct_alignments_sort(plist: &mut Option<Box<BlastCompoAlignment>>, hspcnt: i32) {
    if COMPO_INTENSE_DEBUG {
        assert_eq!(s_distinct_alignments_length(plist.as_deref()), hspcnt);
    }
    if hspcnt > 1 {
        let mut list = plist.take();
        let leftcnt = hspcnt / 2;
        let rightcnt = hspcnt - leftcnt;
        let mut leftlist = list.take();
        let mut cursor = &mut leftlist;
        for _ in 0..leftcnt - 1 {
            cursor = &mut cursor.as_mut().expect("non-null split list").next;
        }
        let mut rightlist = cursor.as_mut().expect("non-null split list").next.take();

        if COMPO_INTENSE_DEBUG {
            assert_eq!(s_distinct_alignments_length(rightlist.as_deref()), rightcnt);
            assert_eq!(s_distinct_alignments_length(leftlist.as_deref()), leftcnt);
        }

        if leftcnt > 1 {
            s_distinct_alignments_sort(&mut leftlist, leftcnt);
        }
        if rightcnt > 1 {
            s_distinct_alignments_sort(&mut rightlist, rightcnt);
        }

        let mut merged: Option<Box<BlastCompoAlignment>> = None;
        let mut tail = &mut merged;
        while leftlist.is_some() || rightlist.is_some() {
            if leftlist.is_none() {
                *tail = rightlist.take();
                break;
            } else if rightlist.is_none() {
                *tail = leftlist.take();
                break;
            } else {
                let take_left = s_alignment_cmp(
                    leftlist.as_deref().expect("left alignment"),
                    rightlist.as_deref().expect("right alignment"),
                ) < 0;
                let mut elt = if take_left {
                    let mut elt = leftlist.take().expect("left alignment");
                    leftlist = elt.next.take();
                    elt
                } else {
                    let mut elt = rightlist.take().expect("right alignment");
                    rightlist = elt.next.take();
                    elt
                };
                elt.next = None;
                *tail = Some(elt);
                tail = &mut tail.as_mut().expect("merged tail").next;
            }
        }
        *plist = merged;
        if COMPO_INTENSE_DEBUG {
            assert_eq!(s_distinct_alignments_length(plist.as_deref()), hspcnt);
            assert_eq!(s_alignments_are_sorted(plist.as_deref()), 1);
        }
    }
}

pub fn s_alignment_copy(align: &BlastCompoAlignment) -> Option<Box<BlastCompoAlignment>> {
    blast_compo_alignment_new(
        align.score,
        align.matrix_adjust_rule,
        align.query_start,
        align.query_end,
        align.query_index,
        align.match_start,
        align.match_end,
        align.frame,
        align.context,
    )
}

pub fn s_is_same_end_point(new_align: &BlastCompoAlignment, align: &BlastCompoAlignment) -> bool {
    assert_eq!(new_align.frame, align.frame);
    (align.query_start == new_align.query_start && align.match_start == new_align.match_start)
        || (align.query_end == new_align.query_end && align.match_end == new_align.match_end)
}

pub fn s_is_similar_end_point(
    new_align: &BlastCompoAlignment,
    align: &BlastCompoAlignment,
) -> bool {
    let start_contained = (align.query_start <= new_align.query_start
        && align.query_end >= new_align.query_start)
        && (align.match_start <= new_align.match_start && align.match_end >= new_align.match_start);
    let end_contained = (align.query_start <= new_align.query_end
        && align.query_end >= new_align.query_end)
        && (align.match_start <= new_align.match_end && align.match_end >= new_align.match_end);

    (start_contained
        && new_align.query_start - new_align.match_start == align.query_start - align.match_start)
        || (end_contained
            && new_align.query_end - new_align.match_end == align.query_end - align.match_end)
}

pub fn s_with_distinct_ends(
    p_new_align: &mut Option<Box<BlastCompoAlignment>>,
    p_old_alignments: &mut Option<Box<BlastCompoAlignment>>,
    free_align_context: Option<FreeAlignTracebackType>,
    is_same_adjustment: bool,
) {
    let Some(mut new_align) = p_new_align.take() else {
        return;
    };
    let mut include_new_align = true;
    let mut align = p_old_alignments.as_deref();
    while let Some(current) = align {
        if current.frame == new_align.frame
            && ((is_same_adjustment && s_is_same_end_point(&new_align, current))
                || (!is_same_adjustment && s_is_similar_end_point(&new_align, current)))
        {
            if new_align.score <= current.score {
                include_new_align = false;
                break;
            }
        }
        align = current.next.as_deref();
    }

    if include_new_align {
        let mut old_alignments = p_old_alignments.take();
        let new_frame = new_align.frame;
        let new_query_start = new_align.query_start;
        let new_match_start = new_align.match_start;
        let new_query_end = new_align.query_end;
        let new_match_end = new_align.match_end;
        new_align.next = None;
        let mut output = Some(new_align);
        let mut tail = &mut output.as_mut().expect("new alignment").next;

        while let Some(mut current) = old_alignments {
            let align_next = current.next.take();
            if current.frame == new_frame
                && ((current.query_start == new_query_start
                    && current.match_start == new_match_start)
                    || (current.query_end == new_query_end && current.match_end == new_match_end))
            {
                let mut single = Some(current);
                blast_compo_alignments_free(&mut single, free_align_context);
            } else {
                *tail = Some(current);
                tail = &mut tail.as_mut().expect("tail alignment").next;
            }
            old_alignments = align_next;
        }
        *p_old_alignments = output;
    } else {
        let mut rejected = Some(new_align);
        blast_compo_alignments_free(&mut rejected, free_align_context);
    }
}

#[derive(Debug, Clone)]
pub struct BlastCompoGappingParams {
    pub gap_open: i32,
    pub gap_extend: i32,
    pub decline_align: i32,
    pub x_dropoff: i32,
    pub context: *mut c_void,
}

#[derive(Debug, Clone, Copy)]
pub struct BlastCompoSequenceRange {
    pub begin: i32,
    pub end: i32,
    pub context: i32,
}

#[derive(Debug, Clone, Copy)]
pub struct BlastCompoSequenceData {
    pub data: *mut u8,
    pub length: i32,
    pub buffer: *mut u8,
}

#[derive(Debug, Clone, Copy)]
pub struct BlastCompoMatchingSequence {
    pub length: i32,
    pub index: i32,
    pub local_data: *mut c_void,
}

#[derive(Debug, Clone)]
pub struct BlastCompoQueryInfo {
    pub origin: i32,
    pub seq: BlastCompoSequenceData,
    pub composition: BlastAminoAcidComposition,
    pub eff_search_space: f64,
    pub words: *mut u64,
}

pub type CalcLambdaType = fn(probs: &[f64], min_score: i32, max_score: i32, lambda0: f64) -> f64;

pub type GetRangeType = fn(
    sequence: *const BlastCompoMatchingSequence,
    range: *const BlastCompoSequenceRange,
    data: *mut BlastCompoSequenceData,
    orig_query: *const BlastCompoSequenceData,
    q_range: *const BlastCompoSequenceRange,
    q_data: *mut BlastCompoSequenceData,
    query_words: *const u64,
    align: *const BlastCompoAlignment,
    should_test_identical: bool,
    compo_adjust_mode: ECompoAdjustModes,
    is_smith_waterman: bool,
    subject_maybe_biased: *mut bool,
) -> i32;

pub type RedoOneAlignmentType = fn(
    in_align: *mut BlastCompoAlignment,
    matrix_adjust_rule: MatrixAdjustRule,
    query_data: *mut BlastCompoSequenceData,
    query_range: *mut BlastCompoSequenceRange,
    ccat_query_length: i32,
    subject_data: *mut BlastCompoSequenceData,
    subject_range: *mut BlastCompoSequenceRange,
    full_subject_length: i32,
    gapping_params: *mut BlastCompoGappingParams,
) -> Option<Box<BlastCompoAlignment>>;

pub type NewXdropAlignType = fn(
    palign: *mut Option<Box<BlastCompoAlignment>>,
    pquery_end: *mut i32,
    pmatch_end: *mut i32,
    query_start: i32,
    match_start: i32,
    score: i32,
    query: *mut BlastCompoSequenceData,
    query_range: *mut BlastCompoSequenceRange,
    ccat_query_length: i32,
    subject: *mut BlastCompoSequenceData,
    subject_range: *mut BlastCompoSequenceRange,
    full_subject_length: i32,
    gapping_params: *mut BlastCompoGappingParams,
    matrix_adjust_rule: MatrixAdjustRule,
) -> i32;

pub type FreeAlignTracebackType = fn(traceback_data: *mut c_void);

#[derive(Debug, Clone, Copy)]
pub struct BlastRedoAlignCallbacks {
    pub calc_lambda: Option<CalcLambdaType>,
    pub get_range: Option<GetRangeType>,
    pub redo_one_alignment: Option<RedoOneAlignmentType>,
    pub new_xdrop_align: Option<NewXdropAlignType>,
    pub free_align_traceback: Option<FreeAlignTracebackType>,
}

#[derive(Debug, Clone)]
pub struct BlastRedoAlignParams {
    pub matrix_info: Option<Box<BlastMatrixInfo>>,
    pub gapping_params: Option<Box<BlastCompoGappingParams>>,
    pub compo_adjust_mode: ECompoAdjustModes,
    pub position_based: i32,
    pub re_pseudocounts: i32,
    pub subject_is_translated: i32,
    pub query_is_translated: i32,
    pub ccat_query_length: i32,
    pub cutoff_s: i32,
    pub cutoff_e: f64,
    pub do_link_hsps: i32,
    pub callbacks: *const BlastRedoAlignCallbacks,
    pub near_identical_cutoff: f64,
}

#[allow(clippy::too_many_arguments)]
pub fn blast_redo_align_params_new(
    pmatrix_info: &mut Option<Box<BlastMatrixInfo>>,
    pgapping_params: &mut Option<Box<BlastCompoGappingParams>>,
    compo_adjust_mode: ECompoAdjustModes,
    position_based: i32,
    query_is_translated: i32,
    subject_is_translated: i32,
    ccat_query_length: i32,
    cutoff_s: i32,
    cutoff_e: f64,
    do_link_hsps: i32,
    callbacks: *const BlastRedoAlignCallbacks,
    near_identical_cutoff: f64,
) -> Option<Box<BlastRedoAlignParams>> {
    let matrix_info = pmatrix_info.take();
    let gapping_params = pgapping_params.take();
    if matrix_info.is_none() || gapping_params.is_none() {
        return None;
    }

    Some(Box::new(BlastRedoAlignParams {
        matrix_info,
        gapping_params,
        compo_adjust_mode,
        position_based,
        re_pseudocounts: K_RE_MATRIX_ADJUSTMENT_PSEUDOCOUNTS,
        subject_is_translated,
        query_is_translated,
        ccat_query_length,
        cutoff_s,
        cutoff_e,
        do_link_hsps,
        callbacks,
        near_identical_cutoff,
    }))
}

pub fn blast_redo_align_params_free(pparams: &mut Option<Box<BlastRedoAlignParams>>) {
    if let Some(params) = pparams.as_mut() {
        params.matrix_info = None;
        params.gapping_params = None;
    }
    *pparams = None;
}

pub fn s_sequence_data_release(self_: &mut BlastCompoSequenceData) {
    if !self_.buffer.is_null() {
        unsafe {
            free(self_.buffer.cast::<c_void>());
        }
    }
    self_.data = ptr::null_mut();
    self_.buffer = ptr::null_mut();
}

#[derive(Debug, Clone)]
pub struct SWindowInfo {
    pub query_range: BlastCompoSequenceRange,
    pub subject_range: BlastCompoSequenceRange,
    pub align: Option<Box<BlastCompoAlignment>>,
    pub hspcnt: i32,
}

pub fn s_window_info_new(
    begin: i32,
    end: i32,
    context: i32,
    query_origin: i32,
    query_length: i32,
    query_index: i32,
    align: Option<Box<BlastCompoAlignment>>,
) -> Option<Box<SWindowInfo>> {
    let mut hspcnt = 0;
    let mut current = align.as_deref();
    while let Some(align) = current {
        hspcnt += 1;
        current = align.next.as_deref();
    }

    Some(Box::new(SWindowInfo {
        subject_range: BlastCompoSequenceRange {
            begin,
            end,
            context,
        },
        query_range: BlastCompoSequenceRange {
            begin: query_origin,
            end: query_origin + query_length,
            context: query_index,
        },
        align,
        hspcnt,
    }))
}

pub fn s_window_swap_range(self_: &mut SWindowInfo) {
    let range = BlastCompoSequenceRange {
        begin: self_.subject_range.begin,
        end: self_.subject_range.end,
        context: self_.subject_range.context,
    };
    self_.subject_range.begin = self_.query_range.begin;
    self_.subject_range.end = self_.query_range.end;
    self_.subject_range.context = self_.query_range.context;
    self_.query_range.begin = range.begin;
    self_.query_range.end = range.end;
    self_.query_range.context = range.context;
}

pub fn s_window_info_free(window: &mut Option<Box<SWindowInfo>>) {
    if let Some(window) = window.as_mut() {
        blast_compo_alignments_free(&mut window.align, None);
    }
    *window = None;
}

pub fn s_window_info_join(win1: &mut SWindowInfo, pwin2: &mut Option<Box<SWindowInfo>>) {
    let Some(mut win2) = pwin2.take() else {
        return;
    };

    assert_eq!(win1.subject_range.context, win2.subject_range.context);
    assert_eq!(win1.query_range.context, win2.query_range.context);

    win1.subject_range.begin = win1.subject_range.begin.min(win2.subject_range.begin);
    win1.subject_range.end = win1.subject_range.end.max(win2.subject_range.end);
    win1.hspcnt += win2.hspcnt;

    let mut tail = &mut win1.align;
    while tail.as_ref().is_some() {
        tail = &mut tail.as_mut().expect("alignment tail").next;
    }
    *tail = win2.align.take();

    let mut win2 = Some(win2);
    s_window_info_free(&mut win2);
    *pwin2 = None;
}

pub fn s_location_compare_windows(w1: &SWindowInfo, w2: &SWindowInfo) -> i32 {
    let sr1 = &w1.subject_range;
    let sr2 = &w2.subject_range;
    let qr1 = &w1.query_range;
    let qr2 = &w2.query_range;

    let mut result = qr1.context.cmp(&qr2.context) as i32;
    if result == 0 {
        result = sr1.context.cmp(&sr2.context) as i32;
    }
    if result == 0 {
        result = sr1.begin.cmp(&sr2.begin) as i32;
    }
    if result == 0 {
        result = sr1.end.cmp(&sr2.end) as i32;
    }
    if result == 0 {
        result = qr1.begin.cmp(&qr2.begin) as i32;
    }
    if result == 0 {
        result = qr1.end.cmp(&qr2.end) as i32;
    }
    result
}

pub fn s_subject_compare_windows(w1: &SWindowInfo, w2: &SWindowInfo) -> i32 {
    let sr1 = &w1.subject_range;
    let sr2 = &w2.subject_range;
    let qr1 = &w1.query_range;
    let qr2 = &w2.query_range;

    let mut result = sr1.begin.cmp(&sr2.begin) as i32;
    if result == 0 {
        result = sr1.end.cmp(&sr2.end) as i32;
    }
    if result == 0 {
        result = sr1.context.cmp(&sr2.context) as i32;
    }
    if result == 0 {
        result = qr1.begin.cmp(&qr2.begin) as i32;
    }
    if result == 0 {
        result = qr1.end.cmp(&qr2.end) as i32;
    }
    if result == 0 {
        result = qr1.context.cmp(&qr2.context) as i32;
    }
    result
}

pub fn s_get_translated_length(length: i32, frame: i32, is_pos_based: i32) -> i32 {
    if is_pos_based != 0 {
        let f = if frame < 0 {
            frame.abs() + 2
        } else {
            frame - 1
        };
        let nucl_length = (length - 5) / 2 + 2;
        return (nucl_length - f % 3) / 3;
    }

    (length - frame.abs() + 1) / 3
}

#[allow(clippy::too_many_arguments)]
pub fn s_windows_from_translated_aligns(
    alignments: Option<&BlastCompoAlignment>,
    query_info: &[BlastCompoQueryInfo],
    hspcnt: i32,
    border: i32,
    sequence_length: i32,
    pwindows: &mut Option<Vec<Option<Box<SWindowInfo>>>>,
    n_windows: &mut i32,
    subject_is_translated: i32,
    is_pos_based: i32,
) -> i32 {
    let mut windows: Vec<Option<Box<SWindowInfo>>> = vec![None; hspcnt.max(0) as usize];
    *n_windows = hspcnt;

    let mut align = alignments;
    let mut k = 0usize;
    while let Some(current) = align {
        let frame = current.frame;
        let query_index = current.query_index;
        let Some(query) = query_info.get(query_index as usize) else {
            *pwindows = None;
            return -1;
        };
        let query_length = query.seq.length;
        let translated_length = s_get_translated_length(sequence_length, frame, is_pos_based);

        let align_copy = s_alignment_copy(current);
        if align_copy.is_none() {
            *pwindows = None;
            return -1;
        }

        if subject_is_translated != 0 {
            let begin = 0.max(current.match_start - border);
            let end = translated_length.min(current.match_end + border);
            windows[k] =
                s_window_info_new(begin, end, frame, 0, query_length, query_index, align_copy);
        } else {
            let begin = 0.max(current.query_start - border);
            let end = query_length.min(current.query_end + border);
            windows[k] =
                s_window_info_new(begin, end, query_index, 0, sequence_length, 0, align_copy);
        }

        if windows[k].is_none() {
            *pwindows = None;
            return -1;
        }

        align = current.next.as_deref();
        k += 1;
    }

    windows.sort_by(|a, b| {
        let result = s_location_compare_windows(
            a.as_deref().expect("window"),
            b.as_deref().expect("window"),
        );
        result.cmp(&0)
    });

    let mut length_joined = 0usize;
    for k in 0..hspcnt.max(0) as usize {
        let mut window = windows[k].take();
        let should_join = if k + 1 < hspcnt.max(0) as usize {
            let next_window = windows[k + 1].as_deref().expect("next window");
            let window_ref = window.as_deref().expect("window");
            window_ref.subject_range.context == next_window.subject_range.context
                && window_ref.query_range.context == next_window.query_range.context
                && window_ref.subject_range.end >= next_window.subject_range.begin
        } else {
            false
        };

        if should_join {
            s_window_info_join(windows[k + 1].as_mut().expect("next window"), &mut window);
        } else {
            windows[length_joined] = window;
            length_joined += 1;
        }
    }
    *n_windows = length_joined as i32;

    for k in length_joined..hspcnt.max(0) as usize {
        windows[k] = None;
    }

    if subject_is_translated == 0 {
        for k in 0..length_joined {
            s_window_swap_range(windows[k].as_mut().expect("joined window"));
        }
    }

    for k in 0..length_joined {
        let window = windows[k].as_mut().expect("joined window");
        s_distinct_alignments_sort(&mut window.align, window.hspcnt);
    }
    windows.truncate(length_joined);
    windows.sort_by(|a, b| {
        let result =
            s_subject_compare_windows(a.as_deref().expect("window"), b.as_deref().expect("window"));
        result.cmp(&0)
    });
    *pwindows = Some(windows);
    0
}

pub fn s_windows_from_protein_aligns(
    alignments: Option<&BlastCompoAlignment>,
    query_info: &[BlastCompoQueryInfo],
    num_queries: i32,
    sequence_length: i32,
    pwindows: &mut Option<Vec<Option<Box<SWindowInfo>>>>,
    n_windows: &mut i32,
) -> i32 {
    let mut windows: Vec<Option<Box<SWindowInfo>>> = vec![None; num_queries.max(0) as usize];
    *n_windows = 0;
    *n_windows = num_queries;

    let mut align = alignments;
    while let Some(current) = align {
        let query_index = current.query_index;
        let Some(query) = query_info.get(query_index as usize) else {
            *pwindows = None;
            return -1;
        };
        let query_length = query.seq.length;

        if windows[query_index as usize].is_none() {
            windows[query_index as usize] =
                s_window_info_new(0, sequence_length, 0, 0, query_length, query_index, None);
            if windows[query_index as usize].is_none() {
                *pwindows = None;
                return -1;
            }
        }

        let Some(mut copied_align) = s_alignment_copy(current) else {
            *pwindows = None;
            return -1;
        };
        let window = windows[query_index as usize]
            .as_mut()
            .expect("query window");
        copied_align.next = window.align.take();
        window.align = Some(copied_align);
        window.hspcnt += 1;

        align = current.next.as_deref();
    }

    let mut window_index = 0usize;
    for query_index in 0..num_queries.max(0) as usize {
        if windows[query_index].is_some() {
            windows.swap(window_index, query_index);
            s_alignments_rev(&mut windows[window_index].as_mut().expect("window").align);
            window_index += 1;
        }
    }

    if window_index == 0 {
        *pwindows = None;
        return -1;
    }

    windows.truncate(window_index);
    *n_windows = window_index as i32;
    windows.sort_by(|a, b| {
        let result =
            s_subject_compare_windows(a.as_deref().expect("window"), b.as_deref().expect("window"));
        result.cmp(&0)
    });
    *pwindows = Some(windows);
    0
}

#[allow(clippy::too_many_arguments)]
pub fn s_windows_from_aligns(
    alignments: Option<&BlastCompoAlignment>,
    query_info: &[BlastCompoQueryInfo],
    hspcnt: i32,
    num_queries: i32,
    border: i32,
    sequence_length: i32,
    pwindows: &mut Option<Vec<Option<Box<SWindowInfo>>>>,
    n_windows: &mut i32,
    query_is_translated: i32,
    subject_is_translated: i32,
    is_pos_based: i32,
) -> i32 {
    if subject_is_translated != 0 || query_is_translated != 0 {
        s_windows_from_translated_aligns(
            alignments,
            query_info,
            hspcnt,
            border,
            sequence_length,
            pwindows,
            n_windows,
            subject_is_translated,
            is_pos_based,
        )
    } else {
        s_windows_from_protein_aligns(
            alignments,
            query_info,
            num_queries,
            sequence_length,
            pwindows,
            n_windows,
        )
    }
}

pub fn s_get_composition(
    composition: &mut BlastAminoAcidComposition,
    alphsize: i32,
    seq: &BlastCompoSequenceData,
    range: &BlastCompoSequenceRange,
    align: &BlastCompoAlignment,
    query_is_translated: bool,
    subject_is_translated: bool,
) {
    let data = seq.data;
    let length = range.end - range.begin;
    let sequence = if length > 0 && !data.is_null() {
        unsafe { std::slice::from_raw_parts(data, length as usize) }
    } else {
        &[]
    };

    let (left, right) = if query_is_translated || subject_is_translated {
        let start = if query_is_translated {
            align.query_start
        } else {
            align.match_start
        } - range.begin;
        let end = if query_is_translated {
            align.query_end
        } else {
            align.match_end
        } - range.begin;
        blast_get_composition_range(sequence, start, end)
    } else {
        (0, length)
    };

    let left = left.clamp(0, sequence.len() as i32) as usize;
    let right = right.clamp(left as i32, sequence.len() as i32) as usize;
    let alphsize = alphsize.max(0) as usize;
    let (prob, num_true_amino_acids) = blast_read_aa_composition(&sequence[left..right], alphsize);

    for i in 0..prob.len() {
        composition.prob[i] = prob[i];
    }
    composition.num_true_amino_acids = num_true_amino_acids as i32;
}

pub fn s_evalue_from_score(score: i32, lambda: f64, log_k: f64, searchsp: f64) -> f64 {
    searchsp * (-(lambda * score as f64) + log_k).exp()
}

pub fn s_is_contained(
    in_align: &BlastCompoAlignment,
    alignments: Option<&BlastCompoAlignment>,
    lambda: f64,
) -> bool {
    let query_offset = in_align.query_start;
    let query_end = in_align.query_end;
    let subject_offset = in_align.match_start;
    let subject_end = in_align.match_end;
    let score = in_align.score as f64;
    let score_thresh = score + KAPPA_BIT_TOL * LOCAL_LN2 / lambda;

    let mut align = alignments;
    while let Some(current) = align {
        let in_sign = if in_align.frame > 0 {
            1
        } else if in_align.frame < 0 {
            -1
        } else {
            0
        };
        let current_sign = if current.frame > 0 {
            1
        } else if current.frame < 0 {
            -1
        } else {
            0
        };
        if in_sign == current_sign
            && ((current.query_start <= query_offset && current.query_end >= query_offset)
                && (current.match_start <= subject_offset && current.match_end >= subject_offset))
            && ((current.query_start <= query_end && current.query_end >= query_end)
                && (current.match_start <= subject_end && current.match_end >= subject_end))
            && score_thresh <= current.score as f64
        {
            return true;
        }
        align = current.next.as_deref();
    }
    false
}

pub fn s_preliminary_test_near_identical(
    query_info: &[BlastCompoQueryInfo],
    window: &SWindowInfo,
    align: &BlastCompoAlignment,
    cutoff: f64,
) -> bool {
    let query_index = align.query_index;
    let Some(query_info) = query_info.get(query_index as usize) else {
        return false;
    };
    let query_length = query_info.seq.length;

    if cutoff > 0.0 {
        if align.match_end - align.match_start + 1 < query_length.min(MINIMUM_LENGTH_NEAR_IDENTICAL)
        {
            return false;
        }

        let align_len =
            (align.query_end - align.query_start).min(align.match_end - align.match_start);

        if (align.score as f64) / (align_len as f64) < cutoff {
            return false;
        }
    } else {
        if window.hspcnt > 1 || window.hspcnt < 1 {
            return false;
        }

        if align.query_end - align.query_start != align.match_end - align.match_start {
            return false;
        }

        if align.match_end - align.match_start + 1 < query_length.min(MINIMUM_LENGTH_NEAR_IDENTICAL)
        {
            return false;
        }
    }

    true
}

pub fn blast_compo_early_termination(
    evalue: f64,
    significant_matches: &[BlastCompoHeap],
    num_queries: i32,
) -> i32 {
    for i in 0..num_queries {
        let Some(heap) = significant_matches.get(i as usize) else {
            return 0;
        };
        if blast_compo_heap_filled_to_cutoff(heap) {
            let ecutoff = heap.ecutoff;
            if evalue <= EVALUE_STRETCH as f64 * ecutoff {
                return 0;
            }
        } else {
            return 0;
        }
    }
    1
}

#[allow(clippy::too_many_arguments)]
pub fn blast_redo_one_match(
    alignments: &mut [Option<Box<BlastCompoAlignment>>],
    params: &mut BlastRedoAlignParams,
    incoming_aligns: Option<&BlastCompoAlignment>,
    hspcnt: i32,
    lambda: f64,
    matching_seq: &mut BlastCompoMatchingSequence,
    ccat_query_length: i32,
    query_info: &mut [BlastCompoQueryInfo],
    num_queries: i32,
    matrix: &mut [Vec<i32>],
    alphsize: i32,
    nr_record: &mut BlastCompositionWorkspace,
    pvalue_for_this_pair: &mut f64,
    composition_test_index: i32,
    lambda_ratio: &mut f64,
) -> i32 {
    let mut status;
    let mut windows: Option<Vec<Option<Box<SWindowInfo>>>> = None;
    let mut n_windows = 0;
    let mut matrix_adjust_rule = MatrixAdjustRule::DontAdjust;

    let Some(scaled_matrix_info) = params.matrix_info.as_ref() else {
        return -1;
    };
    let compo_adjust_mode = params.compo_adjust_mode;
    let re_pseudocounts = params.re_pseudocounts;
    let query_is_translated = params.query_is_translated;
    let subject_is_translated = params.subject_is_translated;
    let Some(gapping_params) = params.gapping_params.as_mut() else {
        return -1;
    };
    assert!(!params.callbacks.is_null());
    let callbacks = unsafe { &*params.callbacks };
    let calc_lambda = callbacks
        .calc_lambda
        .expect("Blast_RedoOneMatch calc_lambda callback");
    let get_range = callbacks
        .get_range
        .expect("Blast_RedoOneMatch get_range callback");
    let redo_one_alignment = callbacks
        .redo_one_alignment
        .expect("Blast_RedoOneMatch redo_one_alignment callback");

    assert!((compo_adjust_mode as i32) < 2 || params.position_based == 0);

    for query_index in 0..num_queries.max(0) as usize {
        if let Some(slot) = alignments.get_mut(query_index) {
            blast_compo_alignments_free(slot, callbacks.free_align_traceback);
        }
    }

    status = s_windows_from_aligns(
        incoming_aligns,
        query_info,
        hspcnt,
        num_queries,
        K_WINDOW_BORDER,
        matching_seq.length,
        &mut windows,
        &mut n_windows,
        query_is_translated,
        subject_is_translated,
        params.position_based,
    );
    if status != 0 {
        return status;
    }

    let mut windows = windows.take().unwrap_or_default();
    for window_index in 0..n_windows.max(0) as usize {
        let Some(window) = windows.get_mut(window_index).and_then(Option::as_mut) else {
            continue;
        };
        let mut subject = BlastCompoSequenceData {
            data: ptr::null_mut(),
            length: 0,
            buffer: ptr::null_mut(),
        };
        let mut query = BlastCompoSequenceData {
            data: ptr::null_mut(),
            length: 0,
            buffer: ptr::null_mut(),
        };
        let mut near_identical_status = true;
        let mut old_near_identical_status = false;
        let mut subject_maybe_biased = true;
        let mut num_adjustments = 0;

        let query_index = window
            .align
            .as_deref()
            .map(|align| align.query_index)
            .unwrap_or(0);
        if query_index < 0 || query_index as usize >= query_info.len() {
            status = -1;
            break;
        }

        let mut hsp_index = 0;
        let mut in_align_ptr: *mut BlastCompoAlignment = window
            .align
            .as_deref_mut()
            .map(|align| align as *mut BlastCompoAlignment)
            .unwrap_or(ptr::null_mut());
        while !in_align_ptr.is_null() {
            let in_align = unsafe { &mut *in_align_ptr };

            if hsp_index == 0 || subject_maybe_biased {
                near_identical_status = s_preliminary_test_near_identical(
                    query_info,
                    window,
                    in_align,
                    params.near_identical_cutoff,
                );
            }

            if hsp_index == 0
                || (subject_maybe_biased && near_identical_status != old_near_identical_status)
            {
                s_sequence_data_release(&mut subject);
                s_sequence_data_release(&mut query);
                status = get_range(
                    matching_seq as *const BlastCompoMatchingSequence,
                    &window.subject_range as *const BlastCompoSequenceRange,
                    &mut subject as *mut BlastCompoSequenceData,
                    &query_info[query_index as usize].seq as *const BlastCompoSequenceData,
                    &window.query_range as *const BlastCompoSequenceRange,
                    &mut query as *mut BlastCompoSequenceData,
                    query_info[query_index as usize].words as *const u64,
                    in_align as *const BlastCompoAlignment,
                    near_identical_status,
                    compo_adjust_mode,
                    false,
                    &mut subject_maybe_biased as *mut bool,
                );
                if status != 0 {
                    break;
                }
            }

            if query_is_translated != 0 {
                s_get_composition(
                    &mut query_info[query_index as usize].composition,
                    alphsize,
                    &query,
                    &window.query_range,
                    in_align,
                    true,
                    false,
                );
            }

            let contained = alignments
                .get(query_index as usize)
                .and_then(Option::as_deref)
                .is_some_and(|existing| s_is_contained(in_align, Some(existing), lambda));

            if !contained {
                let mut adjust_search_failed = 0;
                if compo_adjust_mode != ECompoAdjustModes::NoCompositionBasedStats
                    && (subject_is_translated != 0
                        || hsp_index == 0
                        || near_identical_status != old_near_identical_status)
                {
                    let mut subject_composition = BlastAminoAcidComposition {
                        prob: [0.0; crate::algo::blast::composition_adjustment::composition_adjustment::COMPO_LARGEST_ALPHABET],
                        num_true_amino_acids: 0,
                    };
                    s_get_composition(
                        &mut subject_composition,
                        alphsize,
                        &subject,
                        &window.subject_range,
                        in_align,
                        false,
                        subject_is_translated != 0,
                    );
                    adjust_search_failed = blast_adjust_scores(
                        matrix,
                        &query_info[query_index as usize].composition,
                        query.length,
                        &subject_composition,
                        subject.length,
                        scaled_matrix_info,
                        compo_adjust_mode,
                        re_pseudocounts,
                        nr_record,
                        &mut matrix_adjust_rule,
                        calc_lambda,
                        pvalue_for_this_pair,
                        composition_test_index,
                        lambda_ratio,
                    );
                    if adjust_search_failed < 0 {
                        status = adjust_search_failed;
                        break;
                    }
                    num_adjustments += 1;
                }

                if adjust_search_failed == 0 {
                    let mut new_align = redo_one_alignment(
                        in_align as *mut BlastCompoAlignment,
                        matrix_adjust_rule,
                        &mut query as *mut BlastCompoSequenceData,
                        &mut window.query_range as *mut BlastCompoSequenceRange,
                        ccat_query_length,
                        &mut subject as *mut BlastCompoSequenceData,
                        &mut window.subject_range as *mut BlastCompoSequenceRange,
                        matching_seq.length,
                        gapping_params.as_mut() as *mut BlastCompoGappingParams,
                    );
                    if new_align
                        .as_ref()
                        .is_some_and(|align| align.score >= params.cutoff_s)
                    {
                        if let Some(slot) = alignments.get_mut(query_index as usize) {
                            s_with_distinct_ends(
                                &mut new_align,
                                slot,
                                callbacks.free_align_traceback,
                                num_adjustments == 1,
                            );
                        } else {
                            blast_compo_alignments_free(
                                &mut new_align,
                                callbacks.free_align_traceback,
                            );
                            status = -1;
                            break;
                        }
                    } else {
                        blast_compo_alignments_free(&mut new_align, callbacks.free_align_traceback);
                    }
                }
            }

            old_near_identical_status = near_identical_status;
            in_align_ptr = in_align
                .next
                .as_deref_mut()
                .map(|next| next as *mut BlastCompoAlignment)
                .unwrap_or(ptr::null_mut());
            hsp_index += 1;
        }

        if !subject.data.is_null() {
            s_sequence_data_release(&mut subject);
        }
        if !query.data.is_null() {
            s_sequence_data_release(&mut query);
        }
        if status != 0 {
            break;
        }
    }

    if status != 0 {
        for query_index in 0..num_queries.max(0) as usize {
            if let Some(slot) = alignments.get_mut(query_index) {
                blast_compo_alignments_free(slot, callbacks.free_align_traceback);
            }
        }
    }

    for window in windows.iter_mut() {
        s_window_info_free(window);
    }

    status
}

/// NCBI: `Blast_RedoOneMatchSmithWaterman` (`redo_alignment.c:1310`).
#[allow(clippy::too_many_arguments)]
pub fn blast_redo_one_match_smith_waterman(
    alignments: &mut [Option<Box<BlastCompoAlignment>>],
    params: &mut BlastRedoAlignParams,
    incoming_aligns: Option<&BlastCompoAlignment>,
    hspcnt: i32,
    lambda: f64,
    log_k: f64,
    matching_seq: &mut BlastCompoMatchingSequence,
    query_info: &mut [BlastCompoQueryInfo],
    num_queries: i32,
    matrix: &mut [Vec<i32>],
    alphsize: i32,
    nr_record: &mut BlastCompositionWorkspace,
    forbidden: &mut BlastForbiddenRanges,
    significant_matches: &mut [BlastCompoHeap],
    pvalue_for_this_pair: &mut f64,
    composition_test_index: i32,
    lambda_ratio: &mut f64,
) -> i32 {
    let mut status;
    let mut windows: Option<Vec<Option<Box<SWindowInfo>>>> = None;
    let mut n_windows = 0;
    let mut matrix_adjust_rule = MatrixAdjustRule::DontAdjust;

    let Some(scaled_matrix_info) = params.matrix_info.as_ref() else {
        return -1;
    };
    let compo_adjust_mode = params.compo_adjust_mode;
    let position_based = params.position_based;
    let re_pseudocounts = params.re_pseudocounts;
    let query_is_translated = params.query_is_translated;
    let subject_is_translated = params.subject_is_translated;
    let do_link_hsps = params.do_link_hsps;
    let ccat_query_length = params.ccat_query_length;
    let Some(gapping_params) = params.gapping_params.as_mut() else {
        return -1;
    };
    let gap_open = gapping_params.gap_open;
    let gap_extend = gapping_params.gap_extend;
    assert!(!params.callbacks.is_null());
    let callbacks = unsafe { &*params.callbacks };
    let calc_lambda = callbacks
        .calc_lambda
        .expect("Blast_RedoOneMatchSmithWaterman calc_lambda callback");
    let get_range = callbacks
        .get_range
        .expect("Blast_RedoOneMatchSmithWaterman get_range callback");
    let new_xdrop_align = callbacks
        .new_xdrop_align
        .expect("Blast_RedoOneMatchSmithWaterman new_xdrop_align callback");

    assert!((compo_adjust_mode as i32) < 2 || position_based == 0);

    for query_index in 0..num_queries.max(0) as usize {
        if let Some(slot) = alignments.get_mut(query_index) {
            blast_compo_alignments_free(slot, callbacks.free_align_traceback);
        }
    }

    status = s_windows_from_aligns(
        incoming_aligns,
        query_info,
        hspcnt,
        num_queries,
        K_WINDOW_BORDER,
        matching_seq.length,
        &mut windows,
        &mut n_windows,
        query_is_translated,
        subject_is_translated,
        position_based,
    );
    if status != 0 {
        return status;
    }

    let mut windows = windows.take().unwrap_or_default();
    for window_index in 0..n_windows.max(0) as usize {
        let Some(window) = windows.get_mut(window_index).and_then(Option::as_mut) else {
            continue;
        };
        let mut subject = BlastCompoSequenceData {
            data: ptr::null_mut(),
            length: 0,
            buffer: ptr::null_mut(),
        };
        let mut query = BlastCompoSequenceData {
            data: ptr::null_mut(),
            length: 0,
            buffer: ptr::null_mut(),
        };
        let mut adjust_search_failed = 0;

        let query_index = window.query_range.context;
        if query_index < 0 || query_index as usize >= query_info.len() {
            status = -1;
            break;
        }
        let Some(window_align) = window.align.as_deref() else {
            continue;
        };

        let near_identical_status = s_preliminary_test_near_identical(
            query_info,
            window,
            window_align,
            params.near_identical_cutoff,
        );

        status = get_range(
            matching_seq as *const BlastCompoMatchingSequence,
            &window.subject_range as *const BlastCompoSequenceRange,
            &mut subject as *mut BlastCompoSequenceData,
            &query_info[query_index as usize].seq as *const BlastCompoSequenceData,
            &window.query_range as *const BlastCompoSequenceRange,
            &mut query as *mut BlastCompoSequenceData,
            query_info[query_index as usize].words as *const u64,
            window_align as *const BlastCompoAlignment,
            near_identical_status,
            compo_adjust_mode,
            true,
            ptr::null_mut(),
        );
        if status != 0 {
            if !subject.data.is_null() {
                s_sequence_data_release(&mut subject);
            }
            if !query.data.is_null() {
                s_sequence_data_release(&mut query);
            }
            break;
        }

        if query_is_translated != 0 {
            s_get_composition(
                &mut query_info[query_index as usize].composition,
                alphsize,
                &query,
                &window.query_range,
                window_align,
                true,
                false,
            );
        }

        let searchsp = query_info[query_index as usize].eff_search_space;

        if compo_adjust_mode != ECompoAdjustModes::NoCompositionBasedStats {
            let mut subject_composition = BlastAminoAcidComposition {
                prob: [0.0; crate::algo::blast::composition_adjustment::composition_adjustment::COMPO_LARGEST_ALPHABET],
                num_true_amino_acids: 0,
            };
            s_get_composition(
                &mut subject_composition,
                alphsize,
                &subject,
                &window.subject_range,
                window_align,
                false,
                subject_is_translated != 0,
            );
            adjust_search_failed = blast_adjust_scores(
                matrix,
                &query_info[query_index as usize].composition,
                query.length,
                &subject_composition,
                subject.length,
                scaled_matrix_info,
                compo_adjust_mode,
                re_pseudocounts,
                nr_record,
                &mut matrix_adjust_rule,
                calc_lambda,
                pvalue_for_this_pair,
                composition_test_index,
                lambda_ratio,
            );
            if adjust_search_failed < 0 {
                status = adjust_search_failed;
            }
        }

        if status == 0 && adjust_search_failed == 0 {
            let mut square_matrix = [[0i32; AA_SIZE]; AA_SIZE];
            if position_based == 0 {
                if matrix.len() < AA_SIZE
                    || matrix.iter().take(AA_SIZE).any(|row| row.len() < AA_SIZE)
                {
                    status = -1;
                } else {
                    for i in 0..AA_SIZE {
                        square_matrix[i].copy_from_slice(&matrix[i][..AA_SIZE]);
                    }
                }
            }
            blast_forbidden_ranges_clear(forbidden);
            loop {
                let subject_slice: &[u8] = if subject.length <= 0 {
                    &[]
                } else if subject.data.is_null() {
                    status = -1;
                    &[]
                } else {
                    unsafe { std::slice::from_raw_parts(subject.data, subject.length as usize) }
                };
                let query_slice: &[u8] = if query.length <= 0 {
                    &[]
                } else if query.data.is_null() {
                    status = -1;
                    &[]
                } else {
                    unsafe { std::slice::from_raw_parts(query.data, query.length as usize) }
                };
                if status != 0 {
                    break;
                }

                let (a_sw_score, match_end_usize, query_end_usize) = if position_based != 0 {
                    blast_smith_waterman_score_only_pssm(
                        subject_slice,
                        query_slice,
                        matrix,
                        window.query_range.begin.max(0) as usize,
                        gap_open,
                        gap_extend,
                        Some(forbidden),
                    )
                } else {
                    blast_smith_waterman_score_only(
                        subject_slice,
                        query_slice,
                        &square_matrix,
                        gap_open,
                        gap_extend,
                        Some(forbidden),
                    )
                };
                let mut match_end = match_end_usize as i32;
                let mut query_end = query_end_usize as i32;

                let mut alignment_is_significant;
                if do_link_hsps != 0 {
                    alignment_is_significant = a_sw_score >= params.cutoff_s;
                } else {
                    let new_sw_evalue = s_evalue_from_score(a_sw_score, lambda, log_k, searchsp);
                    alignment_is_significant = new_sw_evalue < params.cutoff_e;
                    if alignments
                        .get(query_index as usize)
                        .and_then(Option::as_ref)
                        .is_none()
                    {
                        let Some(heap) = significant_matches.get_mut(query_index as usize) else {
                            status = -1;
                            break;
                        };
                        alignment_is_significant = alignment_is_significant
                            && blast_compo_heap_would_insert(
                                heap,
                                new_sw_evalue,
                                a_sw_score,
                                matching_seq.index,
                            );
                    }
                }

                if alignment_is_significant {
                    let (updated_score, match_start_usize, query_start_usize) =
                        if position_based != 0 {
                            blast_smith_waterman_find_start_pssm(
                                subject_slice,
                                query_slice,
                                matrix,
                                window.query_range.begin.max(0) as usize,
                                gap_open,
                                gap_extend,
                                match_end_usize,
                                query_end_usize,
                                a_sw_score,
                                Some(forbidden),
                            )
                        } else {
                            blast_smith_waterman_find_start(
                                subject_slice,
                                query_slice,
                                &square_matrix,
                                gap_open,
                                gap_extend,
                                match_end_usize,
                                query_end_usize,
                                a_sw_score,
                                Some(forbidden),
                            )
                        };
                    let _ = updated_score;
                    let match_start = match_start_usize as i32;
                    let query_start = query_start_usize as i32;

                    let mut new_align: Option<Box<BlastCompoAlignment>> = None;
                    status = new_xdrop_align(
                        &mut new_align as *mut Option<Box<BlastCompoAlignment>>,
                        &mut query_end as *mut i32,
                        &mut match_end as *mut i32,
                        query_start,
                        match_start,
                        a_sw_score,
                        &mut query as *mut BlastCompoSequenceData,
                        &mut window.query_range as *mut BlastCompoSequenceRange,
                        ccat_query_length,
                        &mut subject as *mut BlastCompoSequenceData,
                        &mut window.subject_range as *mut BlastCompoSequenceRange,
                        matching_seq.length,
                        gapping_params.as_mut() as *mut BlastCompoGappingParams,
                        matrix_adjust_rule,
                    );
                    if status != 0 {
                        blast_compo_alignments_free(&mut new_align, callbacks.free_align_traceback);
                        break;
                    }
                    if let Some(mut align) = new_align {
                        if let Some(slot) = alignments.get_mut(query_index as usize) {
                            align.next = slot.take();
                            *slot = Some(align);
                        } else {
                            let mut align = Some(align);
                            blast_compo_alignments_free(&mut align, callbacks.free_align_traceback);
                            status = -1;
                            break;
                        }
                    }

                    if window.hspcnt > 1 {
                        status = blast_forbidden_ranges_push(
                            forbidden,
                            query_start,
                            query_end,
                            match_start,
                            match_end,
                        );
                    }
                    if status != 0 {
                        break;
                    }
                }

                if !(alignment_is_significant && window.hspcnt > 1) {
                    break;
                }
            }
        }

        if !subject.data.is_null() {
            s_sequence_data_release(&mut subject);
        }
        if !query.data.is_null() {
            s_sequence_data_release(&mut query);
        }
        if status != 0 {
            break;
        }
    }

    if status != 0 {
        for query_index in 0..num_queries.max(0) as usize {
            if let Some(slot) = alignments.get_mut(query_index) {
                blast_compo_alignments_free(slot, callbacks.free_align_traceback);
            }
        }
    }

    for window in windows.iter_mut() {
        s_window_info_free(window);
    }

    status
}
