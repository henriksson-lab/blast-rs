//! Rust equivalent of blast_extend.c and na_ungapped.c
//! Ungapped and gapped extension structures.

use crate::parameters::InitialWordParameters;

pub const DIAGHASH_NUM_BUCKETS: usize = 512;
pub const DIAGHASH_CHAIN_LENGTH: usize = 256;
pub const MIN_INIT_HITLIST_SIZE: usize = 100;

/// Diagonal bookkeeping selector used by `BlastExtendWordNew`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DiagTableType {
    DiagArray,
    DiagHash,
}

/// Structure for keeping last hit information for a diagonal.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct DiagStruct {
    pub last_hit: i32,
    pub flag: bool,
}

/// Diagonal-array bookkeeping used by the ungapped extension word finders.
#[derive(Debug, Clone)]
pub struct BlastDiagTable {
    pub hit_level_array: Vec<DiagStruct>,
    pub hit_len_array: Option<Vec<u8>>,
    pub diag_array_length: i32,
    pub diag_mask: i32,
    pub offset: i32,
    pub window: i32,
    pub multiple_hits: bool,
    pub actual_window: i32,
}

/// One cell in the diagonal hash chain.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct DiagHashCell {
    pub diag: i32,
    pub level: i32,
    pub hit_saved: bool,
    pub hit_len: i32,
    pub next: u32,
}

/// Hash-table diagonal bookkeeping used by indexed nucleotide search.
#[derive(Debug, Clone)]
pub struct BlastDiagHash {
    pub num_buckets: u32,
    pub occupancy: u32,
    pub capacity: u32,
    pub backbone: Vec<u32>,
    pub chain: Vec<DiagHashCell>,
    pub offset: i32,
    pub window: i32,
}

/// Structure for keeping initial word extension information.
#[derive(Debug, Clone, Default)]
pub struct BlastExtendWord {
    pub diag_table: Option<BlastDiagTable>,
    pub hash_table: Option<BlastDiagHash>,
}

/// NCBI `BlastWordFinderType` callback selector.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum BlastWordFinderType {
    #[default]
    None,
    PhiBlast,
    Protein,
    Rps,
    Nucleotide,
    IndexedMegablast,
}

/// NCBI `BlastGetGappedScoreType` callback selector.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum BlastGetGappedScoreType {
    #[default]
    Standard,
    PhiBlast,
    SmithWaterman,
}

/// NCBI `JumperGappedType` callback selector.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum JumperGappedType {
    #[default]
    None,
    Jumper,
    ShortReadIndexed,
}

/// Rust owner for the NCBI `BlastCoreAuxStruct` search scratch bundle.
#[derive(Debug, Default)]
#[allow(non_snake_case)]
pub struct BlastCoreAuxStruct {
    pub WordFinder: BlastWordFinderType,
    pub GetGappedScore: BlastGetGappedScoreType,
    pub JumperGapped: JumperGappedType,
    pub ewp: Option<BlastExtendWord>,
    pub init_hitlist: Option<InitHitList>,
    pub offset_pairs: Vec<crate::lookup::OffsetPair>,
    pub mapper_wordhits: Option<crate::lookup::MapperWordHits>,
    pub translation_buffer: Option<Vec<u8>>,
    pub translation_table: Option<Vec<u8>>,
    pub translation_table_rc: Option<Vec<u8>>,
}

/// Initial hit from word finding (before ungapped extension).
#[derive(Debug, Clone)]
pub struct InitHsp {
    pub query_offset: i32,
    pub subject_offset: i32,
    pub ungapped_data: Option<UngappedData>,
}

/// PHI-BLAST initial hit from subject pattern scanning.
///
/// C stores this in `BlastOffsetPair::phi_offsets` and then calls
/// `BLAST_SaveInitialHit(init_hitlist, s_start, s_end, NULL)`. The generic
/// `InitHsp` storage is still used underneath, but this view keeps the PHI
/// names visible at call sites.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct PhiInitialHit {
    pub subject_start: i32,
    pub subject_end: i32,
}

/// Data from an ungapped extension.
#[derive(Debug, Clone)]
pub struct UngappedData {
    pub q_start: i32,
    pub s_start: i32,
    pub length: i32,
    pub score: i32,
}

/// Port of `s_GetUngappedHSPContext` (`blast_gapalign.c:2371`).
pub fn s_get_ungapped_hsp_context(
    query_info: &crate::queryinfo::QueryInfo,
    init_hsp: &InitHsp,
) -> i32 {
    crate::queryinfo::bsearch_context_info(init_hsp.query_offset, query_info)
}

/// Port of `s_AdjustInitialHSPOffsets` (`blast_gapalign.c:2384`).
pub fn s_adjust_initial_hsp_offsets(init_hsp: &mut InitHsp, query_start: i32) {
    init_hsp.query_offset -= query_start;
    if let Some(ungapped_data) = init_hsp.ungapped_data.as_mut() {
        ungapped_data.q_start -= query_start;
        debug_assert!(ungapped_data.q_start >= 0);
    }
}

/// Port of `s_SetUpLocalBlastSequenceBlk` (`blast_gapalign.c:2404`).
/// Returns the starting offset of the selected query context.
pub fn s_set_up_local_blast_sequence_blk(
    concatenated_query: &crate::util::BlastSequenceBlk,
    query_info: &crate::queryinfo::QueryInfo,
    context: i32,
    single_query: &mut crate::util::BlastSequenceBlk,
) -> i32 {
    let context = context.max(0) as usize;

    if let Some(oof_sequence) = concatenated_query.oof_sequence.as_ref() {
        let mixed_frame_context = context - context % crate::util::CODON_LENGTH;
        let Some(first_context) = query_info.contexts.get(mixed_frame_context) else {
            *single_query = crate::util::BlastSequenceBlk::default();
            return -1;
        };
        let Some(last_context) = query_info
            .contexts
            .get(mixed_frame_context + crate::util::CODON_LENGTH - 1)
        else {
            *single_query = crate::util::BlastSequenceBlk::default();
            return -1;
        };

        let query_start = first_context.query_offset;
        let query_end = last_context
            .query_offset
            .saturating_add(last_context.query_length);
        let query_length = query_end.saturating_sub(query_start);
        single_query.sequence = None;
        single_query.oof_sequence = slice_sequence(oof_sequence, query_start, query_length);
        single_query.length = query_length;
        return query_start;
    }

    let Some(context_info) = query_info.contexts.get(context) else {
        *single_query = crate::util::BlastSequenceBlk::default();
        return -1;
    };
    let query_start = context_info.query_offset;
    let query_length = context_info.query_length;
    single_query.sequence = concatenated_query
        .sequence
        .as_ref()
        .and_then(|sequence| slice_sequence(sequence, query_start, query_length));
    single_query.oof_sequence = None;
    single_query.length = query_length;
    query_start
}

/// Port of `s_AdjustHspOffsetsAndGetQueryData` (`blast_gapalign.c:2446`).
/// Returns the context containing `init_hsp` after populating `query_out`.
pub fn s_adjust_hsp_offsets_and_get_query_data(
    query: &crate::util::BlastSequenceBlk,
    query_info: &crate::queryinfo::QueryInfo,
    init_hsp: &mut InitHsp,
    query_out: &mut crate::util::BlastSequenceBlk,
) -> i32 {
    let context = s_get_ungapped_hsp_context(query_info, init_hsp);
    let query_start = s_set_up_local_blast_sequence_blk(query, query_info, context, query_out);
    if query_start >= 0 {
        s_adjust_initial_hsp_offsets(init_hsp, query_start);
    }
    context
}

/// Port of `BLAST_GetUngappedHSPList` (`blast_gapalign.c:4822`).
pub fn blast_get_ungapped_hsp_list(
    init_hitlist: Option<&mut InitHitList>,
    query_info: &crate::queryinfo::QueryInfo,
    subject_oid: i32,
    subject_frame: i32,
    hsp_num_max: i32,
    hsp_list: &mut Option<crate::hspstream::HspList>,
) -> i16 {
    let k_hsp_num_max = crate::hspstream::blast_hsp_num_max(false, hsp_num_max);

    let Some(init_hitlist) = init_hitlist else {
        if let Some(hsp_list) = hsp_list.as_mut() {
            hsp_list.hsps.clear();
            // NCBI `s_BlastGetBestEvalue` seeds with `(double)INT4_MAX`.
            hsp_list.best_evalue = i32::MAX as f64;
        }
        return 0;
    };

    for init_hsp in init_hitlist.hits.iter_mut() {
        if init_hsp.ungapped_data.is_none() {
            continue;
        }

        if hsp_list.is_none() {
            let mut new_list = crate::hspstream::blast_hsp_list_new(k_hsp_num_max);
            new_list.oid = subject_oid;
            *hsp_list = Some(new_list);
        }

        let context = s_get_ungapped_hsp_context(query_info, init_hsp);
        let Some(context_info) = query_info.contexts.get(context as usize) else {
            continue;
        };
        s_adjust_initial_hsp_offsets(init_hsp, context_info.query_offset);
        let ungapped_data = init_hsp
            .ungapped_data
            .as_ref()
            .expect("ungapped data checked before offset adjustment");

        let new_hsp = crate::hspstream::Hsp {
            score: ungapped_data.score,
            num_ident: 0,
            bit_score: 0.0,
            evalue: f64::MAX,
            query_offset: ungapped_data.q_start,
            query_end: ungapped_data.q_start + ungapped_data.length,
            query_gapped_start: init_hsp.query_offset,
            subject_offset: ungapped_data.s_start,
            subject_end: ungapped_data.s_start + ungapped_data.length,
            subject_gapped_start: init_hsp.subject_offset,
            context,
            query_frame: context_info.frame,
            subject_frame,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        };
        if let Some(hsp_list) = hsp_list.as_mut() {
            crate::hspstream::blast_hsp_list_save_hsp(hsp_list, new_hsp);
        }
    }

    if let Some(hsp_list) = hsp_list.as_mut() {
        hsp_list
            .hsps
            .sort_by(crate::hspstream::score_compare_hsps);
    }
    0
}

/// Port state for `BLAST_GetGappedScore` (`blast_gapalign.c:3739`).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct BlastGappedStats {
    pub extensions: i32,
}

/// NCBI: `s_BlastProtGappedAlignment` (`blast_gapalign.c:4298`).
fn s_blast_prot_gapped_alignment(
    query: &[u8],
    subject: &[u8],
    seed_q: usize,
    seed_s: usize,
    matrix: &[[i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE],
    gap_open: i32,
    gap_extend: i32,
    x_drop_gapped: i32,
) -> Option<crate::protein::ProteinGappedResult> {
    crate::protein::protein_gapped_align(
        query,
        subject,
        seed_q,
        seed_s,
        matrix,
        gap_open,
        gap_extend,
        x_drop_gapped,
    )
}

/// Port of `BLAST_GetGappedScore` (`blast_gapalign.c:3739`).
pub fn blast_get_gapped_score(
    program_number: crate::program::ProgramType,
    query: Option<&crate::util::BlastSequenceBlk>,
    query_info: Option<&crate::queryinfo::QueryInfo>,
    subject: Option<&crate::util::BlastSequenceBlk>,
    gap_align: Option<&mut crate::protein::BlastGapAlignStruct>,
    score_params: Option<&crate::parameters::ScoringParameters>,
    ext_params: Option<&crate::parameters::ExtensionParameters>,
    hit_params: Option<&crate::parameters::HitSavingParameters>,
    word_params: Option<&crate::parameters::InitialWordParameters>,
    init_hitlist: Option<&mut InitHitList>,
    hsp_list: &mut Option<crate::hspstream::HspList>,
    mut gapped_stats: Option<&mut BlastGappedStats>,
    _fence_hit: Option<&mut bool>,
) -> i16 {
    let (
        Some(query),
        Some(query_info),
        Some(subject),
        Some(gap_align),
        Some(score_params),
        Some(ext_params),
        Some(hit_params),
        Some(word_params),
        Some(init_hitlist),
    ) = (
        query,
        query_info,
        subject,
        gap_align,
        score_params,
        ext_params,
        hit_params,
        word_params,
        init_hitlist,
    )
    else {
        return 1;
    };

    if init_hitlist.total() == 0 {
        return 0;
    }

    let is_prot = program_number != crate::program::BLASTN
        && program_number != crate::program::PHI_BLASTN
        && program_number != crate::program::MAPPING;
    let is_greedy =
        ext_params.options.prelim_gap_ext == crate::options::PrelimGapExt::GreedyScoreOnly;

    let protein_matrix = if is_prot {
        if let Some(sbp) = gap_align.sbp.as_ref() {
            if sbp.matrix.data.len() >= crate::matrix::AA_SIZE
                && sbp
                    .matrix
                    .data
                    .iter()
                    .take(crate::matrix::AA_SIZE)
                    .all(|row| row.len() >= crate::matrix::AA_SIZE)
            {
                let mut matrix = [[0i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE];
                for (row_index, row) in matrix.iter_mut().enumerate() {
                    row.copy_from_slice(&sbp.matrix.data[row_index][..crate::matrix::AA_SIZE]);
                }
                Some(matrix)
            } else {
                return crate::util::BLASTERR_INVALIDPARAM as i16;
            }
        } else {
            match score_params.options.matrix_name.as_deref() {
                None | Some("BLOSUM62") => Some(crate::matrix::BLOSUM62),
                _ => return crate::util::BLASTERR_INVALIDPARAM as i16,
            }
        }
    } else {
        None
    };

    if hsp_list.is_none() {
        let k_hsp_num_max =
            crate::hspstream::blast_hsp_num_max(true, hit_params.options.hsp_num_max);
        let mut new_list = crate::hspstream::blast_hsp_list_new(k_hsp_num_max);
        new_list.oid = subject.oid;
        *hsp_list = Some(new_list);
    }

    blast_init_hit_list_sort_by_score(init_hitlist);
    let init_hsp_array = init_hitlist.hits.clone();
    let chained_kept = if ext_params.options.chaining
        && is_prot
        && score_params.options.matrix_name.as_deref() == Some("BLOSUM62")
    {
        let word_cutoff = word_params
            .cutoffs
            .first()
            .map(|cutoffs| cutoffs.cutoff_score)
            .unwrap_or(word_params.cutoff_score_min)
            .max(1);
        let hit_cutoff = hit_params.cutoff_score_min.max(1);
        s_chaining_alignment(
            &init_hsp_array,
            score_params.gap_open,
            score_params.gap_extend,
            word_cutoff,
            hit_cutoff,
        )
    } else {
        vec![true; init_hsp_array.len()]
    };
    let subject_tree_max = if crate::program::blast_subject_is_translated(program_number)
        && score_params.options.is_ooframe
    {
        2 * (subject
            .length
            .saturating_add(crate::util::CODON_LENGTH as i32))
        + 1
    } else {
        subject.length + 1
    };
    let mut interval_tree = Some(crate::itree::blast_interval_tree_init(
        0,
        query.length.saturating_add(1),
        0,
        subject_tree_max,
    ));
    let mut found_high_score = vec![false; query_info.num_queries.max(0) as usize];
    if !hit_params.low_score.is_empty() {
        for init_hsp in &init_hsp_array {
            let query_index = crate::queryinfo::blast_get_query_index_from_query_offset(
                init_hsp.query_offset,
                program_number,
                query_info,
            );
            if let Some(ungapped_data) = init_hsp.ungapped_data.as_ref() {
                let low_score = hit_params
                    .low_score
                    .get(query_index.max(0) as usize)
                    .copied()
                    .unwrap_or(i32::MAX);
                if ungapped_data.score > low_score {
                    if let Some(found) = found_high_score.get_mut(query_index.max(0) as usize) {
                        *found = true;
                    }
                }
            }
        }
    }

    for (mut init_hsp, keep_chained_hsp) in init_hsp_array.into_iter().zip(chained_kept.into_iter())
    {
        if !keep_chained_hsp {
            continue;
        }
        let mut query_tmp = crate::util::BlastSequenceBlk::default();
        let context = s_adjust_hsp_offsets_and_get_query_data(
            query,
            query_info,
            &mut init_hsp,
            &mut query_tmp,
        );
        let Some(context_info) = query_info.contexts.get(context.max(0) as usize) else {
            continue;
        };
        let query_index =
            crate::queryinfo::blast_get_query_index_from_context(context, program_number);

        if !hit_params.low_score.is_empty()
            && !found_high_score
                .get(query_index.max(0) as usize)
                .copied()
                .unwrap_or(false)
        {
            continue;
        }

        let (q_start, q_end, s_start, s_end, score) =
            if let Some(ungapped_data) = init_hsp.ungapped_data.as_ref() {
                (
                    ungapped_data.q_start,
                    ungapped_data.q_start + ungapped_data.length,
                    ungapped_data.s_start,
                    ungapped_data.s_start + ungapped_data.length,
                    ungapped_data.score,
                )
            } else {
                (
                    init_hsp.query_offset,
                    init_hsp.query_offset,
                    init_hsp.subject_offset,
                    init_hsp.subject_offset,
                    i32::MIN,
                )
            };

        let preliminary_hsp = crate::hspstream::Hsp {
            score,
            num_ident: 0,
            bit_score: 0.0,
            evalue: f64::MAX,
            query_offset: q_start,
            query_end: q_end,
            query_gapped_start: init_hsp.query_offset,
            subject_offset: s_start,
            subject_end: s_end,
            subject_gapped_start: init_hsp.subject_offset,
            context,
            query_frame: context_info.frame,
            subject_frame: subject.frame as i32,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        };
        if crate::itree::blast_interval_tree_contains_hsp(
            interval_tree.as_ref().expect("interval tree"),
            &preliminary_hsp,
            context,
            hit_params.options.min_diag_separation,
        ) {
            continue;
        }

        let cutoff_index = if crate::program::blast_program_is_rps_blast(program_number) {
            subject.oid
        } else {
            context
        };
        let cutoff = hit_params
            .cutoffs
            .get(cutoff_index.max(0) as usize)
            .map(|cutoffs| cutoffs.cutoff_score)
            .unwrap_or(hit_params.cutoff_score_min);
        if let Some(stats) = gapped_stats.as_deref_mut() {
            stats.extensions += 1;
        }

        if init_hsp.ungapped_data.is_none() {
            return crate::util::BLASTERR_INVALIDPARAM as i16;
        }

        let status = if is_prot {
            let Some(query_seq) = query_tmp.sequence.as_deref() else {
                return crate::util::BLASTERR_INVALIDPARAM as i16;
            };
            let Some(subject_seq) = subject.sequence.as_deref() else {
                return crate::util::BLASTERR_INVALIDPARAM as i16;
            };
            let Some(matrix) = protein_matrix.as_ref() else {
                return crate::util::BLASTERR_INVALIDPARAM as i16;
            };
            let q_start_usize = q_start.max(0) as usize;
            let s_start_usize = s_start.max(0) as usize;
            let q_len = q_end.saturating_sub(q_start).max(0) as usize;
            let s_len = s_end.saturating_sub(s_start).max(0) as usize;
            if q_len == 0
                || s_len == 0
                || q_start_usize >= query_seq.len()
                || s_start_usize >= subject_seq.len()
            {
                continue;
            }
            let (seed_q, seed_s) = crate::protein::blast_get_start_for_gapped_alignment(
                query_seq,
                subject_seq,
                q_start_usize,
                q_len,
                s_start_usize,
                s_len,
                matrix,
            );
            let Some(prelim) = s_blast_prot_gapped_alignment(
                query_seq,
                subject_seq,
                seed_q,
                seed_s,
                matrix,
                score_params.gap_open,
                score_params.gap_extend,
                ext_params.gap_x_dropoff,
            ) else {
                continue;
            };
            gap_align.query_start = prelim.query_start as i32;
            gap_align.query_stop = prelim.query_end as i32;
            gap_align.subject_start = prelim.subject_start as i32;
            gap_align.subject_stop = prelim.subject_end as i32;
            gap_align.score = prelim.score;
            gap_align.edit_script = None;
            0
        } else if is_greedy {
            let Some(query_sequence) = query_tmp.sequence.as_ref() else {
                return crate::util::BLASTERR_INVALIDPARAM as i16;
            };
            let Some(subject_sequence) = subject.sequence.as_deref() else {
                return crate::util::BLASTERR_INVALIDPARAM as i16;
            };
            if let Some(ungapped_data) = init_hsp.ungapped_data.as_ref() {
                init_hsp.query_offset = ungapped_data.q_start + ungapped_data.length / 2;
                init_hsp.subject_offset = ungapped_data.s_start + ungapped_data.length / 2;
            }
            let seed_q = init_hsp.query_offset.max(0) as usize;
            let seed_s = init_hsp.subject_offset.max(0) as usize;
            if seed_q >= query_sequence.len() || seed_s >= subject_sequence.len() {
                return crate::util::BLASTERR_INVALIDPARAM as i16;
            }
            match crate::greedy::blast_greedy_gapped_alignment(
                query_sequence,
                subject_sequence,
                seed_q,
                seed_s,
                score_params.reward,
                score_params.penalty,
                ext_params.gap_x_dropoff,
            ) {
                Some((score, query_start, query_stop, subject_start, subject_stop, edit_script)) => {
                    gap_align.score = score;
                    gap_align.query_start = query_start as i32;
                    gap_align.query_stop = query_stop as i32;
                    gap_align.subject_start = subject_start as i32;
                    gap_align.subject_stop = subject_stop as i32;
                    gap_align.greedy_query_seed_start = init_hsp.query_offset;
                    gap_align.greedy_subject_seed_start = init_hsp.subject_offset;
                    gap_align.edit_script = Some(edit_script);
                    init_hsp.query_offset = gap_align.greedy_query_seed_start;
                    init_hsp.subject_offset = gap_align.greedy_subject_seed_start;
                    0
                }
                None => {
                    gap_align.score = i32::MIN;
                    0
                }
            }
        } else {
            let Some(query_sequence) = query_tmp.sequence.as_ref() else {
                return 1;
            };
            let Some(subject_sequence) = subject.sequence.as_deref() else {
                return 1;
            };
            if s_end >= init_hsp.subject_offset.saturating_add(8) {
                init_hsp.subject_offset = init_hsp.subject_offset.saturating_add(3);
                init_hsp.query_offset = init_hsp.query_offset.saturating_add(3);
            }
            let seed_q = init_hsp.query_offset.max(0) as usize;
            let seed_s = init_hsp.subject_offset.max(0) as usize;
            if seed_q >= query_sequence.len() || seed_s >= subject_sequence.len() {
                return 1;
            }
            let (score, query_start, query_stop, subject_start, subject_stop) =
                crate::traceback::s_blast_dyn_prog_nt_gapped_alignment_decoded_extents(
                    query_sequence,
                    subject_sequence,
                    seed_q,
                    seed_s,
                    score_params.reward,
                    score_params.penalty,
                    score_params.gap_open,
                    score_params.gap_extend,
                    ext_params.gap_x_dropoff,
                    score,
                );
            gap_align.score = score;
            gap_align.query_start = query_start as i32;
            gap_align.query_stop = query_stop as i32;
            gap_align.subject_start = subject_start as i32;
            gap_align.subject_stop = subject_stop as i32;
            gap_align.edit_script = None;
            0
        };
        if status != 0 {
            return status;
        }

        if gap_align.score >= cutoff {
            let new_hsp = crate::hspstream::blast_hsp_init(
                gap_align.query_start,
                gap_align.query_stop,
                gap_align.subject_start,
                gap_align.subject_stop,
                init_hsp.query_offset,
                init_hsp.subject_offset,
                context,
                context_info.frame as i16,
                subject.frame,
                gap_align.score,
                gap_align.edit_script.take(),
            );
            if let Some(hsp_list) = hsp_list.as_mut() {
                let status = crate::hspstream::blast_hsp_list_save_hsp(hsp_list, new_hsp);
                if status != 0 {
                    return status as i16;
                }
                let Some(saved_hsp) = hsp_list.hsps.last() else {
                    return crate::diagnostics::BLASTERR_MEMORY;
                };
                let tree_status = crate::itree::blast_interval_tree_add_hsp(
                    saved_hsp,
                    interval_tree.as_mut().expect("interval tree"),
                    context,
                );
                if tree_status != 0 {
                    return tree_status as i16;
                }
            }
        } else if is_greedy {
            gap_align.edit_script = None;
        }
    }

    let _ = crate::itree::blast_interval_tree_free(&mut interval_tree);
    if let Some(hsp_list) = hsp_list.as_mut() {
        hsp_list
            .hsps
            .sort_by(crate::hspstream::score_compare_hsps);
    }
    0
}

fn slice_sequence(sequence: &[u8], start: i32, length: i32) -> Option<Vec<u8>> {
    let start = start.max(0) as usize;
    let length = length.max(0) as usize;
    let end = start.checked_add(length)?;
    if end <= sequence.len() {
        Some(sequence[start..end].to_vec())
    } else {
        None
    }
}

/// List of initial HSPs from the word-finding + ungapped extension phase.
#[derive(Debug)]
pub struct InitHitList {
    pub hits: Vec<InitHsp>,
    pub allocated: i32,
    pub do_not_reallocate: bool,
}

/// Node used by the gapped-alignment chaining workspace.
#[derive(Debug, Clone)]
pub struct BlastInitHspNode {
    pub init_hsp: Option<InitHsp>,
    pub best_score: i32,
    pub next: Option<usize>,
}

/// Rust model of NCBI `ChainingStruct` (`blast_gapalign_priv.h`).
#[derive(Debug, Clone, Default)]
pub struct ChainingStruct {
    pub nodes: Vec<BlastInitHspNode>,
    pub num_allocated: u32,
}

impl InitHitList {
    pub fn new() -> Self {
        InitHitList {
            hits: Vec::with_capacity(MIN_INIT_HITLIST_SIZE),
            allocated: MIN_INIT_HITLIST_SIZE as i32,
            do_not_reallocate: false,
        }
    }

    pub fn add(&mut self, hsp: InitHsp) {
        self.hits.push(hsp);
    }

    pub fn total(&self) -> usize {
        self.hits.len()
    }

    pub fn reset(&mut self) {
        self.hits.clear();
    }
}

/// Port of `s_BlastDiagTableNew` (`blast_extend.c:47`).
pub fn s_blast_diag_table_new(
    qlen: i32,
    multiple_hits: bool,
    window_size: i32,
) -> Option<BlastDiagTable> {
    let mut diag_array_length = 1i32;
    let target = qlen.saturating_add(window_size).max(1);
    while diag_array_length < target {
        diag_array_length = diag_array_length.checked_shl(1)?;
    }

    Some(BlastDiagTable {
        hit_level_array: Vec::new(),
        hit_len_array: None,
        diag_array_length,
        diag_mask: diag_array_length - 1,
        offset: window_size,
        window: window_size,
        multiple_hits,
        actual_window: 0,
    })
}

/// Port of `s_BlastDiagTableFree` (`blast_extend.c:79`).
pub fn s_blast_diag_table_free(diag_table: &mut Option<BlastDiagTable>) -> Option<BlastDiagTable> {
    if let Some(table) = diag_table.as_mut() {
        table.hit_level_array.clear();
        table.hit_len_array = None;
    }
    *diag_table = None;
    None
}

/// Port of `s_BlastDiagClear` (`blast_extend.c:92`).
pub fn s_blast_diag_clear(diag: Option<&mut BlastDiagTable>) -> i32 {
    let Some(diag) = diag else {
        return 0;
    };

    diag.offset = diag.window;
    for entry in &mut diag.hit_level_array {
        entry.flag = false;
        entry.last_hit = -diag.window;
    }
    if let Some(hit_len_array) = diag.hit_len_array.as_mut() {
        hit_len_array.fill(0);
    }
    0
}

/// Port of `BlastExtendWordNew` (`blast_extend.c:115`).
pub fn blast_extend_word_new(
    query_length: u32,
    word_params: &InitialWordParameters,
    diag_table_type: DiagTableType,
) -> Option<BlastExtendWord> {
    let window_size = word_params.options.window_size;
    match diag_table_type {
        DiagTableType::DiagArray => {
            let multiple_hits = window_size > 0;
            let mut diag_table =
                s_blast_diag_table_new(query_length as i32, multiple_hits, window_size)?;
            let diag_len = diag_table.diag_array_length as usize;

            diag_table.hit_level_array = vec![DiagStruct::default(); diag_len];
            if window_size > 0 {
                diag_table.hit_len_array = Some(vec![0; diag_len]);
            }

            Some(BlastExtendWord {
                diag_table: Some(diag_table),
                hash_table: None,
            })
        }
        DiagTableType::DiagHash => Some(BlastExtendWord {
            diag_table: None,
            hash_table: Some(BlastDiagHash {
                num_buckets: DIAGHASH_NUM_BUCKETS as u32,
                occupancy: 1,
                capacity: DIAGHASH_CHAIN_LENGTH as u32,
                backbone: vec![0; DIAGHASH_NUM_BUCKETS],
                chain: vec![DiagHashCell::default(); DIAGHASH_CHAIN_LENGTH],
                offset: window_size,
                window: window_size,
            }),
        }),
    }
}

/// Port of `s_BlastDiagHashFree` (`blast_extend.c:196`).
pub fn s_blast_diag_hash_free(hash_table: &mut Option<BlastDiagHash>) -> Option<BlastDiagHash> {
    if let Some(hash) = hash_table.as_mut() {
        hash.backbone.clear();
        hash.chain.clear();
        hash.occupancy = 0;
    }
    *hash_table = None;
    None
}

/// Port of `Blast_ExtendWordExit` (`blast_extend.c:167`).
pub fn blast_extend_word_exit(ewp: Option<&mut BlastExtendWord>, subject_length: i32) -> i16 {
    let Some(ewp) = ewp else {
        return -1;
    };

    if let Some(diag_table) = ewp.diag_table.as_mut() {
        if diag_table.offset >= i32::MAX / 4 {
            let _ = s_blast_diag_clear(Some(diag_table));
        } else {
            diag_table.offset = diag_table
                .offset
                .saturating_add(subject_length)
                .saturating_add(diag_table.window);
        }
    } else if let Some(hash_table) = ewp.hash_table.as_mut() {
        if hash_table.offset >= i32::MAX / 4 {
            hash_table.occupancy = 1;
            hash_table.offset = hash_table.window;
            hash_table.backbone.fill(0);
        } else {
            hash_table.offset = hash_table
                .offset
                .saturating_add(subject_length)
                .saturating_add(hash_table.window);
        }
    }

    0
}

/// Port of `BlastExtendWordFree` (`blast_extend.c:208`).
pub fn blast_extend_word_free(ewp: &mut Option<BlastExtendWord>) -> Option<BlastExtendWord> {
    if let Some(extend_word) = ewp.as_mut() {
        let _ = s_blast_diag_table_free(&mut extend_word.diag_table);
        let _ = s_blast_diag_hash_free(&mut extend_word.hash_table);
    }
    *ewp = None;
    None
}

/// Port of NCBI static `s_BlastCoreAuxStructFree`.
pub fn s_blast_core_aux_struct_free(
    aux_struct: &mut Option<BlastCoreAuxStruct>,
) -> Option<BlastCoreAuxStruct> {
    if let Some(aux) = aux_struct.as_mut() {
        let _ = blast_extend_word_free(&mut aux.ewp);
        let _ = blast_init_hit_list_free(&mut aux.init_hitlist);
        aux.offset_pairs.clear();
        let _ = crate::lookup::mapper_word_hits_free(&mut aux.mapper_wordhits);
        aux.translation_buffer = None;
        aux.translation_table = None;
        aux.translation_table_rc = None;
    }
    *aux_struct = None;
    None
}

/// Port of `BLAST_InitHitListNew` (`blast_extend.c:221`).
pub fn blast_init_hit_list_new() -> Option<InitHitList> {
    Some(InitHitList::new())
}

/// Port of `s_BlastInitHitListClean` (`blast_extend.c:247`).
pub fn s_blast_init_hit_list_clean(init_hitlist: &mut InitHitList) {
    blast_init_hit_list_reset(init_hitlist);
    init_hitlist.hits.shrink_to(0);
    init_hitlist.allocated = 0;
}

/// Port of `BlastInitHitListReset` (`blast_extend.c`).
pub fn blast_init_hit_list_reset(init_hitlist: &mut InitHitList) {
    init_hitlist.reset();
}

/// Port of `BlastInitHitListMove` (`blast_extend.c`).
pub fn blast_init_hit_list_move(dst: &mut InitHitList, src: &mut InitHitList) {
    s_blast_init_hit_list_clean(dst);
    dst.hits.append(&mut src.hits);
    dst.allocated = src.allocated;
    dst.do_not_reallocate = src.do_not_reallocate;
    src.allocated = 0;
    src.do_not_reallocate = false;
}

/// Port of `BLAST_InitHitListFree` (`blast_extend.c`).
pub fn blast_init_hit_list_free(init_hitlist: &mut Option<InitHitList>) -> Option<InitHitList> {
    if let Some(list) = init_hitlist.as_mut() {
        blast_init_hit_list_reset(list);
    }
    *init_hitlist = None;
    None
}

/// Port of `BLAST_SaveInitialHit` (`blast_extend.c:330`).
pub fn blast_save_initial_hit(
    init_hitlist: &mut InitHitList,
    q_off: i32,
    s_off: i32,
    ungapped_data: Option<UngappedData>,
) -> bool {
    if init_hitlist.hits.len() >= init_hitlist.allocated.max(0) as usize {
        if init_hitlist.do_not_reallocate {
            return false;
        }
        init_hitlist.allocated = if init_hitlist.allocated > 0 {
            init_hitlist.allocated.saturating_mul(2)
        } else {
            MIN_INIT_HITLIST_SIZE as i32
        };
    }
    init_hitlist.add(InitHsp {
        query_offset: q_off,
        subject_offset: s_off,
        ungapped_data,
    });
    true
}

/// Port of `BlastSaveInitHsp` (`blast_extend.c:367`).
pub fn blast_save_init_hsp(
    ungapped_hsps: &mut InitHitList,
    q_start: i32,
    s_start: i32,
    q_off: i32,
    s_off: i32,
    len: i32,
    score: i32,
) {
    let ungapped_data = UngappedData {
        q_start,
        s_start,
        length: len,
        score,
    };
    let _ = blast_save_initial_hit(ungapped_hsps, q_off, s_off, Some(ungapped_data));
}

/// Comparator used by the C hit-list sorter (`score_compare_match`).
pub fn score_compare_match(a: &InitHsp, b: &InitHsp) -> i32 {
    let Some(data_a) = a.ungapped_data.as_ref() else {
        return if b.ungapped_data.is_none() { 0 } else { 1 };
    };
    let Some(data_b) = b.ungapped_data.as_ref() else {
        return -1;
    };

    if data_a.score > data_b.score {
        -1
    } else if data_a.score < data_b.score {
        1
    } else if data_a.s_start < data_b.s_start {
        -1
    } else if data_a.s_start > data_b.s_start {
        1
    } else if data_a.length > data_b.length {
        -1
    } else if data_a.length < data_b.length {
        1
    } else if data_a.q_start < data_b.q_start {
        -1
    } else if data_a.q_start > data_b.q_start {
        1
    } else {
        0
    }
}

/// Port of `Blast_InitHitListSortByScore` (`blast_extend.c:311`).
pub fn blast_init_hit_list_sort_by_score(init_hitlist: &mut InitHitList) {
    init_hitlist
        .hits
        .sort_by(|a, b| score_compare_match(a, b).cmp(&0));
}

/// Port of `Blast_InitHitListIsSortedByScore` (`blast_extend.c:317`).
pub fn blast_init_hit_list_is_sorted_by_score(init_hitlist: &InitHitList) -> bool {
    init_hitlist
        .hits
        .windows(2)
        .all(|window| score_compare_match(&window[0], &window[1]) <= 0)
}

fn s_blast_init_hsp_gap_score(init_hsp: &InitHsp) -> Option<i32> {
    init_hsp.ungapped_data.as_ref().map(|data| data.score)
}

fn s_blast_init_hsp_query_start(init_hsp: &InitHsp) -> Option<i32> {
    init_hsp.ungapped_data.as_ref().map(|data| data.q_start)
}

fn s_blast_init_hsp_query_end(init_hsp: &InitHsp) -> Option<i32> {
    init_hsp
        .ungapped_data
        .as_ref()
        .map(|data| data.q_start.saturating_add(data.length))
}

fn s_blast_init_hsp_subject_start(init_hsp: &InitHsp) -> Option<i32> {
    init_hsp.ungapped_data.as_ref().map(|data| data.s_start)
}

fn s_blast_init_hsp_subject_end(init_hsp: &InitHsp) -> Option<i32> {
    init_hsp
        .ungapped_data
        .as_ref()
        .map(|data| data.s_start.saturating_add(data.length))
}

/// Port of `s_ChainingAlignment` (`blast_gapalign.c:3592`) over `BlastInitHSP`.
fn s_chaining_alignment(
    init_hsp_array: &[InitHsp],
    gap_open: i32,
    gap_extend: i32,
    word_cutoff: i32,
    hit_cutoff: i32,
) -> Vec<bool> {
    let mut chained_hsps: Vec<usize> = init_hsp_array
        .iter()
        .enumerate()
        .filter_map(|(index, init_hsp)| {
            init_hsp.ungapped_data.as_ref()?;
            Some(index)
        })
        .collect();
    chained_hsps.sort_by(|&a, &b| {
        s_compare_init_hsps_by_query_offset_score(init_hsp_array.get(a), init_hsp_array.get(b))
            .cmp(&0)
    });

    let gap_score = gap_open + gap_extend;
    let mut best_score: Vec<i32> = chained_hsps
        .iter()
        .map(|&index| s_blast_init_hsp_gap_score(&init_hsp_array[index]).unwrap_or(0))
        .collect();
    let mut kept = vec![true; init_hsp_array.len()];

    for k in (0..chained_hsps.len()).rev() {
        let k_hsp = &init_hsp_array[chained_hsps[k]];
        let Some(k_q_start) = s_blast_init_hsp_query_start(k_hsp) else {
            continue;
        };
        let Some(k_q_end) = s_blast_init_hsp_query_end(k_hsp) else {
            continue;
        };
        let Some(k_s_start) = s_blast_init_hsp_subject_start(k_hsp) else {
            continue;
        };
        let Some(k_s_end) = s_blast_init_hsp_subject_end(k_hsp) else {
            continue;
        };
        let self_score = best_score[k];

        for j in (k + 1)..chained_hsps.len() {
            let j_hsp = &init_hsp_array[chained_hsps[j]];
            let Some(j_q_start) = s_blast_init_hsp_query_start(j_hsp) else {
                continue;
            };
            let Some(j_s_start) = s_blast_init_hsp_subject_start(j_hsp) else {
                continue;
            };

            let q_diff = j_q_start - k_q_start + (k_q_end - k_q_start);
            let s_diff = j_s_start - k_s_start + (k_s_end - k_s_start);
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

    for (&index, &score) in chained_hsps.iter().zip(best_score.iter()) {
        if score - gap_score + word_cutoff - 1 < hit_cutoff {
            kept[index] = false;
        }
    }

    kept
}

/// Port of NCBI internal `s_TranslateHSPsToDNAPCoord`
/// (`blast_engine.c:324`).
///
/// Converts initial HSP offsets from translated protein coordinates back into
/// mixed DNA/protein coordinates for out-of-frame gapped extension, then
/// resorts the initial-hit list by score.
pub fn s_translate_hsps_to_dnapcoord(
    program: crate::program::ProgramType,
    init_hitlist: &mut InitHitList,
    query_info: &crate::queryinfo::QueryInfo,
    subject_frame: i16,
    subject_length: i32,
    offset: i32,
) {
    let codon_length = crate::util::CODON_LENGTH as i32;

    for init_hsp in &mut init_hitlist.hits {
        if program == crate::program::BLASTX {
            let contexts = &query_info.contexts;
            let context_idx =
                crate::queryinfo::bsearch_context_info(init_hsp.query_offset, query_info);
            if context_idx < 0 {
                continue;
            }
            let context_idx = context_idx as usize;
            let Some(context) = contexts.get(context_idx) else {
                continue;
            };

            let frame_idx = context_idx % crate::util::CODON_LENGTH as usize;
            let init_frame_idx = context_idx - frame_idx;
            let Some(init_context) = contexts.get(init_frame_idx) else {
                continue;
            };
            let frame_pos = init_context.query_offset + frame_idx as i32;

            init_hsp.query_offset =
                (init_hsp.query_offset - context.query_offset) * codon_length + frame_pos;
            if let Some(ungapped_data) = init_hsp.ungapped_data.as_mut() {
                ungapped_data.q_start =
                    (ungapped_data.q_start - context.query_offset) * codon_length + frame_pos;
            }
        } else {
            init_hsp.subject_offset += offset;
            if let Some(ungapped_data) = init_hsp.ungapped_data.as_mut() {
                ungapped_data.s_start += offset;
            }

            if subject_frame > 0 {
                let frame_shift = subject_frame as i32 - 1;
                init_hsp.subject_offset = init_hsp.subject_offset * codon_length + frame_shift;
                if let Some(ungapped_data) = init_hsp.ungapped_data.as_mut() {
                    ungapped_data.s_start = ungapped_data.s_start * codon_length + frame_shift;
                }
            } else {
                let frame_shift = subject_length - subject_frame as i32;
                init_hsp.subject_offset = init_hsp.subject_offset * codon_length + frame_shift;
                if let Some(ungapped_data) = init_hsp.ungapped_data.as_mut() {
                    ungapped_data.s_start = ungapped_data.s_start * codon_length + frame_shift;
                }
            }
        }
    }

    blast_init_hit_list_sort_by_score(init_hitlist);
}

/// Port of `s_CompareInitHSPsByQueryOffsetScore` (`blast_gapalign.c:3470`).
pub fn s_compare_init_hsps_by_query_offset_score(a: Option<&InitHsp>, b: Option<&InitHsp>) -> i32 {
    let Some(hsp_a) = a else {
        return if b.is_some() { -1 } else { 0 };
    };
    let Some(hsp_b) = b else {
        return 1;
    };
    let Some(data_a) = hsp_a.ungapped_data.as_ref() else {
        return if hsp_b.ungapped_data.is_some() { -1 } else { 0 };
    };
    let Some(data_b) = hsp_b.ungapped_data.as_ref() else {
        return 1;
    };

    if data_a.q_start < data_b.q_start {
        -1
    } else if data_a.q_start > data_b.q_start {
        1
    } else if data_a.score > data_b.score {
        -1
    } else if data_a.score < data_b.score {
        1
    } else {
        0
    }
}

/// Port of `ChainingStructNew` (`blast_gapalign.c:3455`).
pub fn chaining_struct_new() -> Option<ChainingStruct> {
    let mut chaining = ChainingStruct::default();
    chaining.nodes = Vec::new();
    chaining.num_allocated = 0;
    Some(chaining)
}

/// Port of `ChainingStructFree` (`blast_gapalign.c:3442`).
pub fn chaining_struct_free(ch: &mut Option<ChainingStruct>) -> Option<ChainingStruct> {
    if let Some(chaining) = ch.as_mut() {
        chaining.nodes.clear();
        chaining.num_allocated = 0;
    }
    *ch = None;
    None
}

/// Perform ungapped extension on packed nucleotide sequences.
///
/// Extends a seed hit in both directions until the score drops below
/// the x-dropoff threshold. Returns the ungapped alignment data.
pub fn na_ungapped_extend(
    query: &[u8],   // BLASTNA encoded
    subject: &[u8], // NCBI2na packed (4 bases/byte)
    q_offset: i32,  // seed position in query
    s_offset: i32,  // seed position in subject
    reward: i32,
    penalty: i32,
    x_dropoff: i32,
) -> Option<UngappedData> {
    na_ungapped_extend_len(
        query,
        subject,
        subject.len() * 4,
        q_offset,
        s_offset,
        reward,
        penalty,
        x_dropoff,
    )
}

/// Ungapped extension with explicit subject length in bases.
pub fn na_ungapped_extend_len(
    query: &[u8],       // BLASTNA encoded
    subject: &[u8],     // NCBI2na packed (4 bases/byte)
    subject_len: usize, // actual number of bases
    q_offset: i32,      // seed position in query
    s_offset: i32,      // seed position in subject
    reward: i32,
    penalty: i32,
    x_dropoff: i32,
) -> Option<UngappedData> {
    // Extend right from the seed
    let q_len = query.len() as i32;
    let s_len_bases = subject_len as i32;

    let mut score = 0i32;
    let mut sum = 0i32;
    let mut best_len_right = 0i32;
    let mut qi = q_offset;
    let mut si = s_offset;
    let x_dropoff_neg = -x_dropoff;
    let mut x_current = x_dropoff_neg;

    while qi < q_len && si < s_len_bases {
        let q_base = query[qi as usize];
        // Decode subject base from packed format
        let s_base = crate::encoding::ncbi2na_base_at(subject, si as usize);

        sum += blastna_score(q_base, s_base, reward, penalty);

        if sum > 0 {
            score += sum;
            best_len_right = qi - q_offset + 1;
            x_current = (-score).max(x_dropoff_neg);
            sum = 0;
        } else if sum < x_current {
            break;
        }
        qi += 1;
        si += 1;
    }

    // Extend left from the seed
    let mut score_left = 0i32;
    let mut sum_left = 0i32;
    let mut best_len_left = 0i32;
    qi = q_offset - 1;
    si = s_offset - 1;

    while qi >= 0 && si >= 0 {
        let q_base = query[qi as usize];
        let s_base = crate::encoding::ncbi2na_base_at(subject, si as usize);

        sum_left += blastna_score(q_base, s_base, reward, penalty);

        if sum_left > 0 {
            score_left += sum_left;
            best_len_left = q_offset - qi;
            sum_left = 0;
        } else if sum_left < x_dropoff_neg {
            break;
        }
        qi -= 1;
        si -= 1;
    }

    let total_score = score + score_left;
    if total_score <= 0 {
        return None;
    }

    Some(UngappedData {
        q_start: q_offset - best_len_left,
        s_start: s_offset - best_len_left,
        length: best_len_left + best_len_right,
        score: total_score,
    })
}

#[inline(always)]
fn blastna_score(a: u8, b: u8, reward: i32, penalty: i32) -> i32 {
    crate::encoding::blastna_pair_score(a, b, reward, penalty)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_chaining_struct_new_and_free_match_c_shape() {
        let mut chaining = chaining_struct_new();
        let ch = chaining.as_mut().expect("chaining workspace");
        assert!(ch.nodes.is_empty());
        assert_eq!(ch.num_allocated, 0);
        ch.nodes.push(BlastInitHspNode {
            init_hsp: Some(InitHsp {
                query_offset: 1,
                subject_offset: 2,
                ungapped_data: Some(UngappedData {
                    q_start: 1,
                    s_start: 2,
                    length: 3,
                    score: 11,
                }),
            }),
            best_score: 11,
            next: None,
        });
        ch.num_allocated = ch.nodes.capacity() as u32;

        assert!(chaining_struct_free(&mut chaining).is_none());
        assert!(chaining.is_none());
    }

    #[test]
    fn test_init_hitlist() {
        let mut list = InitHitList::new();
        list.add(InitHsp {
            query_offset: 0,
            subject_offset: 0,
            ungapped_data: None,
        });
        assert_eq!(list.total(), 1);
        list.reset();
        assert_eq!(list.total(), 0);
    }

    #[test]
    fn test_init_hitlist_c_lifecycle_helpers() {
        let fresh = blast_init_hit_list_new().expect("init hit list");
        assert_eq!(fresh.total(), 0);
        assert_eq!(fresh.allocated, MIN_INIT_HITLIST_SIZE as i32);
        assert!(!fresh.do_not_reallocate);

        let low = InitHsp {
            query_offset: 1,
            subject_offset: 2,
            ungapped_data: Some(UngappedData {
                q_start: 1,
                s_start: 2,
                length: 3,
                score: 7,
            }),
        };
        let high = InitHsp {
            query_offset: 4,
            subject_offset: 5,
            ungapped_data: Some(UngappedData {
                q_start: 4,
                s_start: 5,
                length: 6,
                score: 13,
            }),
        };

        assert_eq!(score_compare_match(&high, &low), -1);
        assert_eq!(score_compare_match(&low, &high), 1);
        assert_eq!(score_compare_match(&low, &low), 0);
        assert_eq!(
            s_compare_init_hsps_by_query_offset_score(Some(&low), Some(&high)),
            -1
        );
        assert_eq!(
            s_compare_init_hsps_by_query_offset_score(Some(&high), Some(&low)),
            1
        );
        assert_eq!(
            s_compare_init_hsps_by_query_offset_score(None, Some(&low)),
            -1
        );
        assert_eq!(
            score_compare_match(
                &InitHsp {
                    query_offset: 0,
                    subject_offset: 0,
                    ungapped_data: None,
                },
                &high
            ),
            1
        );

        let mut src = InitHitList::new();
        src.add(low.clone());
        src.add(high.clone());
        let mut dst = InitHitList::new();
        dst.add(InitHsp {
            query_offset: 9,
            subject_offset: 9,
            ungapped_data: None,
        });

        blast_init_hit_list_move(&mut dst, &mut src);
        assert_eq!(src.total(), 0);
        assert_eq!(src.allocated, 0);
        assert_eq!(dst.total(), 2);
        assert_eq!(dst.hits[0].query_offset, low.query_offset);
        assert_eq!(dst.hits[1].query_offset, high.query_offset);

        blast_init_hit_list_reset(&mut dst);
        assert_eq!(dst.total(), 0);

        let mut maybe_list = Some(InitHitList::new());
        maybe_list.as_mut().unwrap().add(high);
        assert!(blast_init_hit_list_free(&mut maybe_list).is_none());
        assert!(maybe_list.is_none());
    }

    #[test]
    fn test_translate_hsps_to_dnapcoord_blastx_uses_context_offsets() {
        let query_info = crate::queryinfo::QueryInfo {
            num_queries: 1,
            max_length: 10,
            min_length: 0,
            contexts: (0..6)
                .map(|idx| crate::queryinfo::ContextInfo {
                    query_offset: 100 + idx * 100,
                    query_length: 50,
                    eff_searchsp: 0,
                    length_adjustment: 0,
                    query_index: 0,
                    frame: idx + 1,
                    is_valid: true,
                    segment_flags: crate::queryinfo::E_NO_SEGMENTS,
                })
                .collect(),
        };
        let mut list = InitHitList::new();
        list.add(InitHsp {
            query_offset: 205,
            subject_offset: 7,
            ungapped_data: Some(UngappedData {
                q_start: 202,
                s_start: 7,
                length: 3,
                score: 11,
            }),
        });

        s_translate_hsps_to_dnapcoord(crate::program::BLASTX, &mut list, &query_info, 0, 0, 0);

        assert_eq!(list.hits[0].query_offset, 116);
        assert_eq!(list.hits[0].ungapped_data.as_ref().unwrap().q_start, 107);
        assert_eq!(list.hits[0].subject_offset, 7);
    }

    #[test]
    fn test_translate_hsps_to_dnapcoord_translated_subject_and_resort() {
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[20]);
        let mut list = InitHitList::new();
        list.add(InitHsp {
            query_offset: 5,
            subject_offset: 4,
            ungapped_data: Some(UngappedData {
                q_start: 5,
                s_start: 4,
                length: 2,
                score: 3,
            }),
        });
        list.add(InitHsp {
            query_offset: 2,
            subject_offset: 1,
            ungapped_data: Some(UngappedData {
                q_start: 2,
                s_start: 1,
                length: 2,
                score: 20,
            }),
        });

        s_translate_hsps_to_dnapcoord(crate::program::TBLASTN, &mut list, &query_info, 2, 0, 2);

        assert_eq!(list.hits[0].ungapped_data.as_ref().unwrap().score, 20);
        assert_eq!(list.hits[0].subject_offset, 10);
        assert_eq!(list.hits[0].ungapped_data.as_ref().unwrap().s_start, 10);
        assert_eq!(list.hits[1].subject_offset, 19);
        assert_eq!(list.hits[1].ungapped_data.as_ref().unwrap().s_start, 19);

        s_translate_hsps_to_dnapcoord(crate::program::TBLASTN, &mut list, &query_info, -1, 100, 0);
        assert_eq!(list.hits[0].subject_offset, 131);
        assert_eq!(list.hits[1].subject_offset, 158);
    }

    #[test]
    fn blast_core_aux_struct_free_clears_owned_scratch() {
        let mut aux = Some(BlastCoreAuxStruct {
            WordFinder: BlastWordFinderType::Nucleotide,
            GetGappedScore: BlastGetGappedScoreType::Standard,
            JumperGapped: JumperGappedType::None,
            ewp: Some(BlastExtendWord {
                diag_table: s_blast_diag_table_new(16, false, 8),
                hash_table: None,
            }),
            init_hitlist: blast_init_hit_list_new(),
            offset_pairs: vec![crate::lookup::OffsetPair {
                query_offset: 3,
                subject_offset: 5,
            }],
            mapper_wordhits: None,
            translation_buffer: Some(vec![1, 2, 3]),
            translation_table: Some(vec![4, 5, 6]),
            translation_table_rc: Some(vec![7, 8, 9]),
        });
        aux.as_mut()
            .unwrap()
            .init_hitlist
            .as_mut()
            .unwrap()
            .add(InitHsp {
                query_offset: 7,
                subject_offset: 11,
                ungapped_data: None,
            });

        assert!(s_blast_core_aux_struct_free(&mut aux).is_none());
        assert!(aux.is_none());
    }

    #[test]
    fn test_save_initial_hit_helpers_match_c_shape() {
        let mut list = InitHitList::new();
        assert!(blast_save_initial_hit(&mut list, 3, 5, None));
        assert_eq!(list.total(), 1);
        assert_eq!(list.hits[0].query_offset, 3);
        assert_eq!(list.hits[0].subject_offset, 5);
        assert!(list.hits[0].ungapped_data.is_none());

        blast_save_init_hsp(&mut list, 7, 11, 13, 17, 19, 23);
        assert_eq!(list.total(), 2);
        let hsp = &list.hits[1];
        assert_eq!(hsp.query_offset, 13);
        assert_eq!(hsp.subject_offset, 17);
        let data = hsp.ungapped_data.as_ref().expect("ungapped data");
        assert_eq!(data.q_start, 7);
        assert_eq!(data.s_start, 11);
        assert_eq!(data.length, 19);
        assert_eq!(data.score, 23);

        assert!(!blast_init_hit_list_is_sorted_by_score(&list));
        blast_init_hit_list_sort_by_score(&mut list);
        assert!(blast_init_hit_list_is_sorted_by_score(&list));
        assert_eq!(list.hits[0].ungapped_data.as_ref().unwrap().score, 23);

        list.allocated = list.total() as i32;
        list.do_not_reallocate = true;
        assert!(!blast_save_initial_hit(&mut list, 1, 2, None));
    }

    #[test]
    fn test_s_get_ungapped_hsp_context_uses_query_offset() {
        let query_info = crate::queryinfo::QueryInfo::new_blastn(&[10, 5]);
        let init_hsp = InitHsp {
            query_offset: 22,
            subject_offset: 7,
            ungapped_data: None,
        };

        assert_eq!(s_get_ungapped_hsp_context(&query_info, &init_hsp), 2);
    }

    #[test]
    fn test_s_adjust_initial_hsp_offsets_rebases_query_offsets() {
        let mut init_hsp = InitHsp {
            query_offset: 31,
            subject_offset: 7,
            ungapped_data: Some(UngappedData {
                q_start: 29,
                s_start: 5,
                length: 4,
                score: 12,
            }),
        };

        s_adjust_initial_hsp_offsets(&mut init_hsp, 20);

        assert_eq!(init_hsp.query_offset, 11);
        assert_eq!(init_hsp.subject_offset, 7);
        let ungapped_data = init_hsp.ungapped_data.unwrap();
        assert_eq!(ungapped_data.q_start, 9);
        assert_eq!(ungapped_data.s_start, 5);
        assert_eq!(ungapped_data.length, 4);
        assert_eq!(ungapped_data.score, 12);
    }

    #[test]
    fn test_blast_get_ungapped_hsp_list_matches_c_field_mapping() {
        let query_info = crate::queryinfo::QueryInfo::new_blastn(&[10, 8]);
        let mut init_hitlist = InitHitList::new();
        init_hitlist.add(InitHsp {
            query_offset: 23,
            subject_offset: 30,
            ungapped_data: Some(UngappedData {
                q_start: 22,
                s_start: 28,
                length: 4,
                score: 11,
            }),
        });
        init_hitlist.add(InitHsp {
            query_offset: 4,
            subject_offset: 18,
            ungapped_data: Some(UngappedData {
                q_start: 2,
                s_start: 16,
                length: 5,
                score: 23,
            }),
        });
        init_hitlist.add(InitHsp {
            query_offset: 7,
            subject_offset: 99,
            ungapped_data: None,
        });

        let mut hsp_list = None;
        assert_eq!(
            blast_get_ungapped_hsp_list(
                Some(&mut init_hitlist),
                &query_info,
                42,
                -1,
                0,
                &mut hsp_list,
            ),
            0
        );

        let hsp_list = hsp_list.expect("HSP list allocated for ungapped hits");
        assert_eq!(hsp_list.oid, 42);
        assert_eq!(hsp_list.hsps.len(), 2);
        assert_eq!(hsp_list.hsps[0].score, 23);
        assert_eq!(hsp_list.hsps[0].query_offset, 2);
        assert_eq!(hsp_list.hsps[0].query_end, 7);
        assert_eq!(hsp_list.hsps[0].query_gapped_start, 4);
        assert_eq!(hsp_list.hsps[0].subject_offset, 16);
        assert_eq!(hsp_list.hsps[0].subject_end, 21);
        assert_eq!(hsp_list.hsps[0].subject_gapped_start, 18);
        assert_eq!(hsp_list.hsps[0].context, 0);
        assert_eq!(hsp_list.hsps[0].query_frame, 1);
        assert_eq!(hsp_list.hsps[0].subject_frame, -1);

        assert_eq!(hsp_list.hsps[1].score, 11);
        assert_eq!(hsp_list.hsps[1].query_offset, 0);
        assert_eq!(hsp_list.hsps[1].query_gapped_start, 1);
        assert_eq!(hsp_list.hsps[1].context, 2);
    }

    #[test]
    fn test_blast_get_ungapped_hsp_list_clears_existing_list_without_hits() {
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[10]);
        let mut hsp_list = Some(crate::hspstream::blast_hsp_list_new(0));
        hsp_list.as_mut().unwrap().add_hsp(crate::hspstream::Hsp {
            score: 7,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 3.0,
            query_offset: 1,
            query_end: 2,
            query_gapped_start: 1,
            subject_offset: 4,
            subject_end: 5,
            subject_gapped_start: 4,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });

        assert_eq!(
            blast_get_ungapped_hsp_list(None, &query_info, 0, 0, 0, &mut hsp_list),
            0
        );
        assert!(hsp_list.as_ref().unwrap().hsps.is_empty());
        // NCBI `s_BlastGetBestEvalue` seeds with `(double)INT4_MAX`.
        assert_eq!(hsp_list.as_ref().unwrap().best_evalue, i32::MAX as f64);
    }

    #[test]
    fn test_blast_get_gapped_score_blastn_saves_hsp() {
        let query_info = crate::queryinfo::QueryInfo::new_blastn(&[8]);
        let query = crate::util::BlastSequenceBlk {
            sequence: Some(vec![0, 1, 2, 3, 0, 1, 2, 3]),
            length: 8,
            ..Default::default()
        };
        let subject = crate::util::BlastSequenceBlk {
            sequence: Some(vec![0, 1, 2, 3, 0, 1, 2, 3]),
            length: 8,
            frame: 1,
            oid: 17,
            ..Default::default()
        };
        let mut gap_align = crate::protein::BlastGapAlignStruct::default();
        let score_params = crate::parameters::ScoringParameters {
            options: crate::options::ScoringOptions::new_blastn(),
            reward: 2,
            penalty: -3,
            gap_open: 5,
            gap_extend: 2,
            shift_pen: i16::MAX as i32,
            scale_factor: 1.0,
        };
        let ext_params = crate::parameters::ExtensionParameters {
            options: crate::options::ExtensionOptions::new_blastn(),
            gap_x_dropoff: 20,
            gap_x_dropoff_final: 20,
            gap_trigger: 0,
        };
        let hit_params = crate::parameters::HitSavingParameters {
            cutoffs: vec![
                crate::parameters::BlastGappedCutoffs {
                    cutoff_score: 1,
                    cutoff_score_max: 1,
                };
                query_info.contexts.len()
            ],
            ..Default::default()
        };
        let word_params = InitialWordParameters {
            options: crate::options::InitialWordOptions::new_blastn(),
            x_dropoff_max: 20,
            cutoff_score_min: 1,
            cutoffs: Vec::new(),
            ungapped_extension: true,
            nucl_score_table: InitialWordParameters::build_nucl_score_table(2, -3),
        };
        let mut init_hitlist = InitHitList::new();
        init_hitlist.add(InitHsp {
            query_offset: 1,
            subject_offset: 1,
            ungapped_data: Some(UngappedData {
                q_start: 0,
                s_start: 0,
                length: 8,
                score: 16,
            }),
        });

        let mut hsp_list = None;
        let mut stats = BlastGappedStats::default();
        assert_eq!(
            blast_get_gapped_score(
                crate::program::BLASTN,
                Some(&query),
                Some(&query_info),
                Some(&subject),
                Some(&mut gap_align),
                Some(&score_params),
                Some(&ext_params),
                Some(&hit_params),
                Some(&word_params),
                Some(&mut init_hitlist),
                &mut hsp_list,
                Some(&mut stats),
                None,
            ),
            0
        );

        let hsp_list = hsp_list.expect("gapped HSP list");
        assert_eq!(stats.extensions, 1);
        assert_eq!(hsp_list.oid, 17);
        assert_eq!(hsp_list.hsps.len(), 1);
        assert_eq!(hsp_list.hsps[0].score, 16);
        assert_eq!(hsp_list.hsps[0].query_offset, 0);
        assert_eq!(hsp_list.hsps[0].subject_offset, 0);
        assert_eq!(hsp_list.hsps[0].query_gapped_start, 1);
        assert_eq!(hsp_list.hsps[0].subject_gapped_start, 1);
        assert_eq!(hsp_list.hsps[0].context, 0);
        assert_eq!(hsp_list.hsps[0].subject_frame, 1);
    }

    fn blastn_gapped_score_fixture() -> (
        crate::queryinfo::QueryInfo,
        crate::util::BlastSequenceBlk,
        crate::util::BlastSequenceBlk,
        crate::parameters::ScoringParameters,
        crate::parameters::ExtensionParameters,
        crate::parameters::HitSavingParameters,
        InitialWordParameters,
    ) {
        let query_info = crate::queryinfo::QueryInfo::new_blastn(&[8]);
        let query = crate::util::BlastSequenceBlk {
            sequence: Some(vec![0, 1, 2, 3, 0, 1, 2, 3]),
            length: 8,
            ..Default::default()
        };
        let subject = crate::util::BlastSequenceBlk {
            sequence: Some(vec![0, 1, 2, 3, 0, 1, 2, 3]),
            length: 8,
            frame: 1,
            oid: 17,
            ..Default::default()
        };
        let score_params = crate::parameters::ScoringParameters {
            options: crate::options::ScoringOptions::new_blastn(),
            reward: 2,
            penalty: -3,
            gap_open: 5,
            gap_extend: 2,
            shift_pen: i16::MAX as i32,
            scale_factor: 1.0,
        };
        let ext_params = crate::parameters::ExtensionParameters {
            options: crate::options::ExtensionOptions::new_blastn(),
            gap_x_dropoff: 20,
            gap_x_dropoff_final: 20,
            gap_trigger: 0,
        };
        let hit_params = crate::parameters::HitSavingParameters {
            cutoffs: vec![
                crate::parameters::BlastGappedCutoffs {
                    cutoff_score: 1,
                    cutoff_score_max: 1,
                };
                query_info.contexts.len()
            ],
            ..Default::default()
        };
        let word_params = InitialWordParameters {
            options: crate::options::InitialWordOptions::new_blastn(),
            x_dropoff_max: 20,
            cutoff_score_min: 1,
            cutoffs: Vec::new(),
            ungapped_extension: true,
            nucl_score_table: InitialWordParameters::build_nucl_score_table(2, -3),
        };
        (
            query_info,
            query,
            subject,
            score_params,
            ext_params,
            hit_params,
            word_params,
        )
    }

    #[test]
    fn test_blast_get_gapped_score_skips_interval_contained_duplicate() {
        let (query_info, query, subject, score_params, ext_params, hit_params, word_params) =
            blastn_gapped_score_fixture();
        let mut gap_align = crate::protein::BlastGapAlignStruct::default();
        let mut init_hitlist = InitHitList::new();
        init_hitlist.add(InitHsp {
            query_offset: 1,
            subject_offset: 1,
            ungapped_data: Some(UngappedData {
                q_start: 0,
                s_start: 0,
                length: 8,
                score: 16,
            }),
        });
        init_hitlist.add(InitHsp {
            query_offset: 3,
            subject_offset: 3,
            ungapped_data: Some(UngappedData {
                q_start: 2,
                s_start: 2,
                length: 3,
                score: 6,
            }),
        });

        let mut hsp_list = None;
        let mut stats = BlastGappedStats::default();
        assert_eq!(
            blast_get_gapped_score(
                crate::program::BLASTN,
                Some(&query),
                Some(&query_info),
                Some(&subject),
                Some(&mut gap_align),
                Some(&score_params),
                Some(&ext_params),
                Some(&hit_params),
                Some(&word_params),
                Some(&mut init_hitlist),
                &mut hsp_list,
                Some(&mut stats),
                None,
            ),
            0
        );

        assert_eq!(stats.extensions, 1);
        assert_eq!(hsp_list.expect("gapped HSP list").hsps.len(), 1);
    }

    #[test]
    fn test_blast_get_gapped_score_greedy_nucleotide_branch() {
        let (query_info, query, subject, score_params, mut ext_params, hit_params, word_params) =
            blastn_gapped_score_fixture();
        ext_params.options.prelim_gap_ext = crate::options::PrelimGapExt::GreedyScoreOnly;
        let mut gap_align = crate::protein::BlastGapAlignStruct::default();
        let mut init_hitlist = InitHitList::new();
        init_hitlist.add(InitHsp {
            query_offset: 1,
            subject_offset: 1,
            ungapped_data: Some(UngappedData {
                q_start: 0,
                s_start: 0,
                length: 8,
                score: 16,
            }),
        });

        let mut hsp_list = None;
        assert_eq!(
            blast_get_gapped_score(
                crate::program::BLASTN,
                Some(&query),
                Some(&query_info),
                Some(&subject),
                Some(&mut gap_align),
                Some(&score_params),
                Some(&ext_params),
                Some(&hit_params),
                Some(&word_params),
                Some(&mut init_hitlist),
                &mut hsp_list,
                None,
                None,
            ),
            0
        );

        let hsps = hsp_list.expect("gapped HSP list").hsps;
        assert_eq!(hsps.len(), 1);
        assert!(hsps[0].score >= 16);
        assert!(hsps[0].edit_script.is_some());
    }

    #[test]
    fn test_blast_get_gapped_score_blastn_restricted_align_does_not_block_nt() {
        let (query_info, query, subject, score_params, ext_params, mut hit_params, word_params) =
            blastn_gapped_score_fixture();
        hit_params.restricted_align = true;
        let mut gap_align = crate::protein::BlastGapAlignStruct::default();
        let mut init_hitlist = InitHitList::new();
        init_hitlist.add(InitHsp {
            query_offset: 1,
            subject_offset: 1,
            ungapped_data: Some(UngappedData {
                q_start: 0,
                s_start: 0,
                length: 8,
                score: 16,
            }),
        });

        let mut hsp_list = None;
        assert_eq!(
            blast_get_gapped_score(
                crate::program::BLASTN,
                Some(&query),
                Some(&query_info),
                Some(&subject),
                Some(&mut gap_align),
                Some(&score_params),
                Some(&ext_params),
                Some(&hit_params),
                Some(&word_params),
                Some(&mut init_hitlist),
                &mut hsp_list,
                None,
                None,
            ),
            0
        );
        assert_eq!(hsp_list.expect("gapped HSP list").hsps.len(), 1);
    }

    #[test]
    fn test_blast_get_gapped_score_protein_branch_saves_hsp() {
        let (query_info, query, subject, score_params, ext_params, hit_params, word_params) =
            blastn_gapped_score_fixture();
        let mut gap_align = crate::protein::BlastGapAlignStruct::default();
        let mut init_hitlist = InitHitList::new();
        init_hitlist.add(InitHsp {
            query_offset: 1,
            subject_offset: 1,
            ungapped_data: Some(UngappedData {
                q_start: 0,
                s_start: 0,
                length: 8,
                score: 16,
            }),
        });
        let mut hsp_list = None;

        assert_eq!(
            blast_get_gapped_score(
                crate::program::BLASTP,
                Some(&query),
                Some(&query_info),
                Some(&subject),
                Some(&mut gap_align),
                Some(&score_params),
                Some(&ext_params),
                Some(&hit_params),
                Some(&word_params),
                Some(&mut init_hitlist),
                &mut hsp_list,
                None,
                None,
            ),
            0
        );
        let hsps = hsp_list.expect("protein gapped HSP list").hsps;
        assert_eq!(hsps.len(), 1);
        assert!(hsps[0].score >= 16);
        assert_eq!(hsps[0].query_gapped_start, 1);
        assert_eq!(hsps[0].subject_gapped_start, 1);
    }

    #[test]
    fn test_blast_get_gapped_score_protein_chaining_filters_low_seed() {
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[16]);
        let query = crate::util::BlastSequenceBlk {
            sequence: Some(vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]),
            length: 16,
            ..Default::default()
        };
        let subject = crate::util::BlastSequenceBlk {
            sequence: query.sequence.clone(),
            length: 16,
            frame: 1,
            oid: 17,
            ..Default::default()
        };
        let mut gap_align = crate::protein::BlastGapAlignStruct::default();
        let score_params = crate::parameters::ScoringParameters {
            options: crate::options::ScoringOptions::new_blastp(),
            reward: 0,
            penalty: 0,
            gap_open: 5,
            gap_extend: 2,
            shift_pen: i16::MAX as i32,
            scale_factor: 1.0,
        };
        let mut ext_params = crate::parameters::ExtensionParameters {
            options: crate::options::ExtensionOptions::new_blastp(),
            gap_x_dropoff: 20,
            gap_x_dropoff_final: 20,
            gap_trigger: 0,
        };
        ext_params.options.chaining = true;
        let hit_params = crate::parameters::HitSavingParameters {
            cutoff_score_min: 16,
            cutoffs: vec![crate::parameters::BlastGappedCutoffs {
                cutoff_score: 16,
                cutoff_score_max: 16,
            }],
            ..Default::default()
        };
        let word_params = InitialWordParameters {
            options: crate::options::InitialWordOptions::new_blastp(),
            x_dropoff_max: 20,
            cutoff_score_min: 10,
            cutoffs: vec![crate::parameters::BlastUngappedCutoffs {
                x_dropoff_init: 20,
                cutoff_score: 10,
            }],
            ungapped_extension: true,
            nucl_score_table: InitialWordParameters::build_nucl_score_table(2, -3),
        };
        let mut init_hitlist = InitHitList::new();
        init_hitlist.add(InitHsp {
            query_offset: 1,
            subject_offset: 1,
            ungapped_data: Some(UngappedData {
                q_start: 0,
                s_start: 0,
                length: 4,
                score: 16,
            }),
        });
        init_hitlist.add(InitHsp {
            query_offset: 9,
            subject_offset: 9,
            ungapped_data: Some(UngappedData {
                q_start: 8,
                s_start: 8,
                length: 4,
                score: 1,
            }),
        });

        let mut hsp_list = None;
        let mut stats = BlastGappedStats::default();
        assert_eq!(
            blast_get_gapped_score(
                crate::program::BLASTP,
                Some(&query),
                Some(&query_info),
                Some(&subject),
                Some(&mut gap_align),
                Some(&score_params),
                Some(&ext_params),
                Some(&hit_params),
                Some(&word_params),
                Some(&mut init_hitlist),
                &mut hsp_list,
                Some(&mut stats),
                None,
            ),
            0
        );

        assert_eq!(stats.extensions, 1);
        assert_eq!(hsp_list.expect("protein gapped HSP list").hsps.len(), 1);
    }

    #[test]
    fn test_blast_get_gapped_score_protein_restricted_align_enters_branch() {
        let (query_info, query, subject, score_params, ext_params, mut hit_params, word_params) =
            blastn_gapped_score_fixture();
        hit_params.restricted_align = true;
        let mut gap_align = crate::protein::BlastGapAlignStruct::default();
        let mut init_hitlist = InitHitList::new();
        init_hitlist.add(InitHsp {
            query_offset: 1,
            subject_offset: 1,
            ungapped_data: Some(UngappedData {
                q_start: 0,
                s_start: 0,
                length: 8,
                score: 16,
            }),
        });
        let mut hsp_list = None;

        assert_eq!(
            blast_get_gapped_score(
                crate::program::BLASTP,
                Some(&query),
                Some(&query_info),
                Some(&subject),
                Some(&mut gap_align),
                Some(&score_params),
                Some(&ext_params),
                Some(&hit_params),
                Some(&word_params),
                Some(&mut init_hitlist),
                &mut hsp_list,
                None,
                None,
            ),
            0
        );
        assert_eq!(hsp_list.expect("protein gapped HSP list").hsps.len(), 1);
    }

    #[test]
    fn test_blast_get_gapped_score_requires_ungapped_data_for_nucleotide_prelim() {
        let (query_info, query, subject, score_params, ext_params, hit_params, word_params) =
            blastn_gapped_score_fixture();
        let mut gap_align = crate::protein::BlastGapAlignStruct::default();
        let mut init_hitlist = InitHitList::new();
        init_hitlist.add(InitHsp {
            query_offset: 1,
            subject_offset: 1,
            ungapped_data: None,
        });
        let mut hsp_list = None;

        assert_eq!(
            blast_get_gapped_score(
                crate::program::BLASTN,
                Some(&query),
                Some(&query_info),
                Some(&subject),
                Some(&mut gap_align),
                Some(&score_params),
                Some(&ext_params),
                Some(&hit_params),
                Some(&word_params),
                Some(&mut init_hitlist),
                &mut hsp_list,
                None,
                None,
            ),
            crate::util::BLASTERR_INVALIDPARAM as i16
        );
    }

    #[test]
    fn test_s_set_up_local_blast_sequence_blk_slices_context_sequence() {
        let query_info = crate::queryinfo::QueryInfo::new_blastn(&[4, 3]);
        let concatenated_query = crate::util::BlastSequenceBlk {
            sequence: Some(vec![
                1, 2, 3, 4, 15, 5, 6, 7, 8, 15, 9, 10, 11, 15, 12, 13, 14,
            ]),
            length: 17,
            ..Default::default()
        };
        let mut single_query = crate::util::BlastSequenceBlk::default();

        let query_start = s_set_up_local_blast_sequence_blk(
            &concatenated_query,
            &query_info,
            2,
            &mut single_query,
        );

        assert_eq!(query_start, 10);
        assert_eq!(single_query.length, 3);
        assert_eq!(single_query.sequence.as_deref(), Some(&[9, 10, 11][..]));
        assert!(single_query.oof_sequence.is_none());
    }

    #[test]
    fn test_s_set_up_local_blast_sequence_blk_handles_oof_context_group() {
        let query_info = crate::queryinfo::QueryInfo {
            num_queries: 1,
            contexts: vec![
                crate::queryinfo::ContextInfo {
                    query_offset: 0,
                    query_length: 2,
                    eff_searchsp: 0,
                    length_adjustment: 0,
                    query_index: 0,
                    frame: 1,
                    is_valid: true,
                    segment_flags: crate::queryinfo::E_NO_SEGMENTS,
                },
                crate::queryinfo::ContextInfo {
                    query_offset: 3,
                    query_length: 2,
                    eff_searchsp: 0,
                    length_adjustment: 0,
                    query_index: 0,
                    frame: 2,
                    is_valid: true,
                    segment_flags: crate::queryinfo::E_NO_SEGMENTS,
                },
                crate::queryinfo::ContextInfo {
                    query_offset: 6,
                    query_length: 2,
                    eff_searchsp: 0,
                    length_adjustment: 0,
                    query_index: 0,
                    frame: 3,
                    is_valid: true,
                    segment_flags: crate::queryinfo::E_NO_SEGMENTS,
                },
            ],
            max_length: 2,
            min_length: 0,
        };
        let concatenated_query = crate::util::BlastSequenceBlk {
            oof_sequence: Some(vec![1, 2, 15, 3, 4, 15, 5, 6]),
            length: 8,
            ..Default::default()
        };
        let mut single_query = crate::util::BlastSequenceBlk::default();

        let query_start = s_set_up_local_blast_sequence_blk(
            &concatenated_query,
            &query_info,
            1,
            &mut single_query,
        );

        assert_eq!(query_start, 0);
        assert_eq!(single_query.length, 8);
        assert!(single_query.sequence.is_none());
        assert_eq!(
            single_query.oof_sequence.as_deref(),
            Some(&[1, 2, 15, 3, 4, 15, 5, 6][..])
        );
    }

    #[test]
    fn test_s_adjust_hsp_offsets_and_get_query_data_rebases_hsp() {
        let query_info = crate::queryinfo::QueryInfo::new_blastn(&[4, 3]);
        let query = crate::util::BlastSequenceBlk {
            sequence: Some(vec![
                1, 2, 3, 4, 15, 5, 6, 7, 8, 15, 9, 10, 11, 15, 12, 13, 14,
            ]),
            length: 17,
            ..Default::default()
        };
        let mut init_hsp = InitHsp {
            query_offset: 12,
            subject_offset: 5,
            ungapped_data: Some(UngappedData {
                q_start: 11,
                s_start: 5,
                length: 2,
                score: 17,
            }),
        };
        let mut query_out = crate::util::BlastSequenceBlk::default();

        let context = s_adjust_hsp_offsets_and_get_query_data(
            &query,
            &query_info,
            &mut init_hsp,
            &mut query_out,
        );

        assert_eq!(context, 2);
        assert_eq!(query_out.sequence.as_deref(), Some(&[9, 10, 11][..]));
        assert_eq!(query_out.length, 3);
        assert_eq!(init_hsp.query_offset, 2);
        assert_eq!(init_hsp.ungapped_data.as_ref().unwrap().q_start, 1);
    }

    #[test]
    fn test_diag_table_and_extend_word_helpers_match_c_shape() {
        let mut diag = s_blast_diag_table_new(33, true, 40).expect("diag table");
        assert_eq!(diag.diag_array_length, 128);
        assert_eq!(diag.diag_mask, 127);
        assert_eq!(diag.offset, 40);
        assert_eq!(diag.window, 40);
        assert!(diag.multiple_hits);

        diag.hit_level_array = vec![
            DiagStruct {
                last_hit: 55,
                flag: true,
            };
            diag.diag_array_length as usize
        ];
        diag.hit_len_array = Some(vec![9; diag.diag_array_length as usize]);
        assert_eq!(s_blast_diag_clear(Some(&mut diag)), 0);
        assert_eq!(diag.offset, diag.window);
        assert!(diag
            .hit_level_array
            .iter()
            .all(|entry| !entry.flag && entry.last_hit == -40));
        assert!(diag
            .hit_len_array
            .as_ref()
            .unwrap()
            .iter()
            .all(|&len| len == 0));

        let mut options = crate::options::InitialWordOptions::new_blastn();
        options.window_size = 12;
        let word_params = InitialWordParameters {
            options,
            x_dropoff_max: 20,
            cutoff_score_min: 8,
            cutoffs: Vec::new(),
            ungapped_extension: true,
            nucl_score_table: InitialWordParameters::build_nucl_score_table(1, -3),
        };
        let extend_word =
            blast_extend_word_new(21, &word_params, DiagTableType::DiagArray).expect("extend word");
        let diag_table = extend_word.diag_table.as_ref().expect("diag table path");
        assert_eq!(diag_table.diag_array_length, 64);
        assert_eq!(diag_table.hit_level_array.len(), 64);
        assert_eq!(diag_table.hit_len_array.as_ref().unwrap().len(), 64);
        assert!(extend_word.hash_table.is_none());

        let mut extend_slot = Some(extend_word);
        assert_eq!(blast_extend_word_exit(extend_slot.as_mut(), 30), 0);
        assert_eq!(
            extend_slot
                .as_ref()
                .unwrap()
                .diag_table
                .as_ref()
                .unwrap()
                .offset,
            54
        );
        assert!(blast_extend_word_free(&mut extend_slot).is_none());
        assert!(extend_slot.is_none());

        let mut diag_slot = Some(diag);
        assert!(s_blast_diag_table_free(&mut diag_slot).is_none());
        assert!(diag_slot.is_none());

        let mut hash_extend =
            blast_extend_word_new(21, &word_params, DiagTableType::DiagHash).expect("hash path");
        assert!(hash_extend.diag_table.is_none());
        let hash_table = hash_extend.hash_table.as_ref().expect("hash table path");
        assert_eq!(hash_table.num_buckets, DIAGHASH_NUM_BUCKETS as u32);
        assert_eq!(hash_table.occupancy, 1);
        assert_eq!(hash_table.capacity, DIAGHASH_CHAIN_LENGTH as u32);
        assert_eq!(hash_table.backbone.len(), DIAGHASH_NUM_BUCKETS);
        assert_eq!(hash_table.chain.len(), DIAGHASH_CHAIN_LENGTH);
        assert_eq!(hash_table.offset, 12);
        assert_eq!(hash_table.window, 12);

        assert_eq!(blast_extend_word_exit(Some(&mut hash_extend), 20), 0);
        assert_eq!(hash_extend.hash_table.as_ref().unwrap().offset, 44);

        {
            let hash_table = hash_extend.hash_table.as_mut().unwrap();
            hash_table.offset = i32::MAX / 4;
            hash_table.occupancy = 77;
            hash_table.backbone.fill(9);
        }
        assert_eq!(blast_extend_word_exit(Some(&mut hash_extend), 20), 0);
        let hash_table = hash_extend.hash_table.as_ref().unwrap();
        assert_eq!(hash_table.offset, hash_table.window);
        assert_eq!(hash_table.occupancy, 1);
        assert!(hash_table.backbone.iter().all(|&bucket| bucket == 0));

        let mut hash_slot = hash_extend.hash_table.take();
        assert!(s_blast_diag_hash_free(&mut hash_slot).is_none());
        assert!(hash_slot.is_none());
    }

    #[test]
    fn test_ungapped_extend_perfect_match() {
        // Query: ACGT (BLASTNA: 0,1,2,3)
        let query = vec![0u8, 1, 2, 3];
        // Subject: ACGT packed = 0b00_01_10_11 = 0x1B
        let subject = vec![0x1Bu8];

        let result = na_ungapped_extend(&query, &subject, 0, 0, 2, -3, 20);
        assert!(result.is_some());
        let data = result.unwrap();
        assert_eq!(data.score, 8); // 4 matches * reward=2
        assert_eq!(data.length, 4);
    }

    #[test]
    fn test_ungapped_extend_with_xdrop() {
        // Build a sequence that starts matching then diverges.
        // Query:   A C G T  A A A A  (BLASTNA: 0,1,2,3,0,0,0,0)
        // Subject: A C G T  C C C C  (NCBI2na: 0,1,2,3,1,1,1,1)
        // First 4 positions match (+2 each = +8), then 4 mismatches (-3 each).
        // With x_dropoff=5 the extension should stop before consuming all mismatches.
        let query = vec![0u8, 1, 2, 3, 0, 0, 0, 0];
        let subject = crate::encoding::pack_ncbi2na_bases(&[0, 1, 2, 3, 1, 1, 1, 1]);

        let result = na_ungapped_extend(&query, &subject, 0, 0, 2, -3, 5);
        assert!(result.is_some());
        let data = result.unwrap();
        // Best score should be from the first 4 matching positions.
        assert_eq!(data.score, 8);
        // The alignment length should be 4 (only the matching region).
        assert_eq!(data.length, 4);
    }

    #[test]
    fn test_ungapped_extend_right_stops_when_total_score_goes_negative() {
        // NCBI's exact nucleotide extender uses a dynamic right-side x-drop:
        // after a +2 first base, a -3 mismatch makes the running total negative
        // and terminates the right extension even when x_dropoff is much larger.
        let query = vec![0u8, 0, 0, 0];
        let subject = crate::encoding::pack_ncbi2na_bases(&[0, 1, 0, 0]);

        let result = na_ungapped_extend_len(&query, &subject, 4, 0, 0, 2, -3, 20);
        assert!(result.is_some());
        let data = result.unwrap();
        assert_eq!(data.score, 2);
        assert_eq!(data.length, 1);
    }

    #[test]
    fn test_ungapped_extend_uses_blastna_ambiguity_score() {
        let query = vec![0u8, 0, 0, 0, 14, 0, 0, 0, 0, 0, 0, 0];
        let subject = crate::encoding::pack_ncbi2na_bases(&[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]);

        let result = na_ungapped_extend_len(&query, &subject, 12, 0, 0, 1, -5, 20);
        assert!(result.is_some());
        let data = result.unwrap();
        assert_eq!(data.score, 7);
        assert_eq!(data.length, 12);
    }

    #[test]
    fn test_ungapped_extend_at_boundary() {
        // Seed at position 0 — left extension has nothing to extend into.
        // Query and subject: ACGT (4 bases, all match)
        let query = vec![0u8, 1, 2, 3];
        let subject = crate::encoding::pack_ncbi2na_bases(&[0, 1, 2, 3]);

        let result = na_ungapped_extend(&query, &subject, 0, 0, 2, -3, 20);
        assert!(result.is_some());
        let data = result.unwrap();
        assert_eq!(data.q_start, 0);
        assert_eq!(data.s_start, 0);
        assert_eq!(data.score, 8);
        assert_eq!(data.length, 4);

        // Seed at the last position — right extension has nothing beyond it.
        let result2 = na_ungapped_extend(&query, &subject, 3, 3, 2, -3, 20);
        assert!(result2.is_some());
        let data2 = result2.unwrap();
        assert_eq!(data2.score, 8);
        assert_eq!(data2.length, 4);
        assert_eq!(data2.q_start, 0);
    }

    #[test]
    fn test_ungapped_extend_all_matches() {
        // 16 bases, all A, perfect match — extension should cover full length.
        let query = vec![0u8; 16];
        let subject = crate::encoding::pack_ncbi2na_bases(&[0u8; 16]);

        let result = na_ungapped_extend_len(&query, &subject, 16, 8, 8, 2, -3, 100);
        assert!(result.is_some());
        let data = result.unwrap();
        assert_eq!(data.score, 32); // 16 * 2
        assert_eq!(data.length, 16);
        assert_eq!(data.q_start, 0);
        assert_eq!(data.s_start, 0);
    }

    #[test]
    fn test_ungapped_extend_score_calculation() {
        // Manually verify score for a known pattern.
        // Query:   A C A C  G  (BLASTNA: 0,1,0,1,2)
        // Subject: A C G C  G  (NCBI2na: 0,1,2,1,2)
        // Position: 0:match(+2) 1:match(+2) 2:mismatch(-3) 3:match(+2) 4:match(+2)
        // Total = 2+2-3+2+2 = 5
        let query = vec![0u8, 1, 0, 1, 2];
        let subject = crate::encoding::pack_ncbi2na_bases(&[0, 1, 2, 1, 2]);

        let result = na_ungapped_extend_len(&query, &subject, 5, 0, 0, 2, -3, 20);
        assert!(result.is_some());
        let data = result.unwrap();
        assert_eq!(data.score, 5);
        assert_eq!(data.length, 5);
    }
}
