//! Rust equivalent of blast_hspstream.c — HSP collection and streaming.
//! Replaces the C vtable + mutex approach with Rust's type system.

use std::sync::Mutex;

use crate::gapinfo::{GapAlignOpType, JumperEditsBlock, SequenceOverhangs};
use crate::options::{HitSavingOptions, ScoringOptions};
use crate::split_query::{
    split_query_blk_allow_gap, split_query_blk_get_chunk_overlap_size,
    split_query_blk_get_context_offsets_for_chunk, split_query_blk_get_query_contexts_for_chunk,
    split_query_blk_get_query_indices_for_chunk, SSplitQueryBlk,
};

pub const K_BLAST_HSP_STREAM_ERROR: i32 = -1;
pub const K_BLAST_HSP_STREAM_SUCCESS: i32 = 0;
pub const K_BLAST_HSP_STREAM_EOF: i32 = 1;
pub const E_NONE: i32 = 0;
pub const E_PRELIM_SEARCH: i32 = 1 << 0;
pub const E_TRACEBACK_SEARCH: i32 = 1 << 1;
pub const E_BOTH: i32 = E_PRELIM_SEARCH | E_TRACEBACK_SEARCH;

/// A single High-Scoring Pair (alignment hit).
#[derive(Debug, Clone)]
pub struct Hsp {
    pub score: i32,
    pub num_ident: i32,
    pub bit_score: f64,
    pub evalue: f64,
    pub query_offset: i32,
    pub query_end: i32,
    pub query_gapped_start: i32,
    pub subject_offset: i32,
    pub subject_end: i32,
    pub subject_gapped_start: i32,
    pub context: i32,
    pub query_frame: i32,
    pub subject_frame: i32,
    pub num_gaps: i32,
    pub comp_adjustment_method: i32,
    pub edit_script: Option<crate::gapinfo::GapEditScript>,
    pub pat_info: Option<PhiPatInfo>,
    pub map_info: Option<BlastHSPMappingInfo>,
}

/// NCBI: `Blast_HSPInit` (`blast_hits.c:151`).
#[allow(clippy::too_many_arguments)]
pub fn blast_hsp_init(
    query_start: i32,
    query_end: i32,
    subject_start: i32,
    subject_end: i32,
    query_gapped_start: i32,
    subject_gapped_start: i32,
    query_context: i32,
    query_frame: i16,
    subject_frame: i16,
    score: i32,
    gap_edit: Option<crate::gapinfo::GapEditScript>,
) -> Hsp {
    let num_gaps = gap_edit
        .as_ref()
        .map(gap_edit_script_num_gap_opens)
        .unwrap_or(0);
    Hsp {
        score,
        num_ident: 0,
        bit_score: 0.0,
        evalue: 0.0,
        query_offset: query_start,
        query_end,
        query_gapped_start,
        subject_offset: subject_start,
        subject_end,
        subject_gapped_start,
        context: query_context,
        query_frame: query_frame as i32,
        subject_frame: subject_frame as i32,
        num_gaps,
        comp_adjustment_method: 0,
        edit_script: gap_edit,
        pat_info: None,
        map_info: None,
    }
}

/// NCBI: `BlastSeg` (`blast_hits.h:96`).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct BlastSeg {
    pub frame: i16,
    pub offset: i32,
    pub end: i32,
    pub gapped_start: i32,
}

impl BlastSeg {
    fn from_query_hsp(hsp: &Hsp) -> Self {
        Self {
            frame: hsp.query_frame as i16,
            offset: hsp.query_offset,
            end: hsp.query_end,
            gapped_start: hsp.query_gapped_start,
        }
    }

    fn from_subject_hsp(hsp: &Hsp) -> Self {
        Self {
            frame: hsp.subject_frame as i16,
            offset: hsp.subject_offset,
            end: hsp.subject_end,
            gapped_start: hsp.subject_gapped_start,
        }
    }
}

/// NCBI: `SPHIHspInfo` (`blast_hits.h:104`).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SPHIHspInfo {
    pub index: i32,
    pub length: i32,
}

impl From<PhiPatInfo> for SPHIHspInfo {
    fn from(value: PhiPatInfo) -> Self {
        Self {
            index: value.index as i32,
            length: value.length,
        }
    }
}

impl From<SPHIHspInfo> for PhiPatInfo {
    fn from(value: SPHIHspInfo) -> Self {
        Self {
            index: value.index.max(0) as usize,
            length: value.length,
        }
    }
}

/// NCBI: `BlastHSP` (`blast_hits.h:124`).
///
/// This is the upstream-shaped internal HSP. Existing Rust pipeline code still
/// uses [`Hsp`] in several places; conversion is intentionally explicit so the
/// translated traceback/CBS/evalue/link/reap path can migrate to `BlastHSP`
/// without carrying public API-shaped data through the core.
#[derive(Debug, Clone)]
pub struct BlastHSP {
    pub score: i32,
    pub num_ident: i32,
    pub bit_score: f64,
    pub evalue: f64,
    pub query: BlastSeg,
    pub subject: BlastSeg,
    pub context: i32,
    pub gap_info: Option<crate::gapinfo::GapEditScript>,
    pub num: i32,
    pub comp_adjustment_method: i16,
    pub pat_info: Option<SPHIHspInfo>,
    pub num_positives: i32,
    pub map_info: Option<BlastHSPMappingInfo>,
}

fn gap_edit_script_num_gap_opens(script: &crate::gapinfo::GapEditScript) -> i32 {
    script
        .iter()
        .filter(|(op, count)| {
            matches!(
                op,
                crate::gapinfo::GapAlignOpType::Ins | crate::gapinfo::GapAlignOpType::Del
            ) && *count > 0
        })
        .count() as i32
}

impl BlastHSP {
    pub fn from_legacy_hsp(hsp: Hsp) -> Self {
        Self {
            score: hsp.score,
            num_ident: hsp.num_ident,
            bit_score: hsp.bit_score,
            evalue: hsp.evalue,
            query: BlastSeg::from_query_hsp(&hsp),
            subject: BlastSeg::from_subject_hsp(&hsp),
            context: hsp.context,
            gap_info: hsp.edit_script,
            num: 0,
            comp_adjustment_method: hsp.comp_adjustment_method as i16,
            pat_info: hsp.pat_info.map(Into::into),
            num_positives: 0,
            map_info: hsp.map_info,
        }
    }

}

/// NCBI: `BlastHSPList` (`blast_hits.h:146`).
#[derive(Debug, Clone)]
pub struct BlastHSPList {
    pub oid: i32,
    pub query_index: i32,
    pub hsp_array: Vec<Option<BlastHSP>>,
    pub hspcnt: i32,
    pub allocated: i32,
    pub hsp_max: i32,
    pub do_not_reallocate: bool,
    pub best_evalue: f64,
}

impl BlastHSPList {
    pub fn from_legacy_hsp_list(list: HspList) -> Self {
        let hsp_array: Vec<Option<BlastHSP>> = list
            .hsps
            .into_iter()
            .map(BlastHSP::from_legacy_hsp)
            .map(Some)
            .collect();
        let hspcnt = hsp_array.len() as i32;
        Self {
            oid: list.oid,
            query_index: 0,
            allocated: hspcnt,
            hsp_array,
            hspcnt,
            hsp_max: list.hsp_max,
            do_not_reallocate: false,
            best_evalue: list.best_evalue,
        }
    }

    pub fn into_legacy_hsp_list(self) -> HspList {
        let hspcnt = self.hspcnt.max(0) as usize;
        HspList {
            oid: self.oid,
            hsps: self
                .hsp_array
                .into_iter()
                .take(hspcnt)
                .flatten()
                .map(|hsp| Hsp {
                    score: hsp.score,
                    num_ident: hsp.num_ident,
                    bit_score: hsp.bit_score,
                    evalue: hsp.evalue,
                    query_offset: hsp.query.offset,
                    query_end: hsp.query.end,
                    query_gapped_start: hsp.query.gapped_start,
                    subject_offset: hsp.subject.offset,
                    subject_end: hsp.subject.end,
                    subject_gapped_start: hsp.subject.gapped_start,
                    context: hsp.context,
                    query_frame: hsp.query.frame as i32,
                    subject_frame: hsp.subject.frame as i32,
                    num_gaps: hsp
                        .gap_info
                        .as_ref()
                        .map(gap_edit_script_num_gap_opens)
                        .unwrap_or(0),
                    comp_adjustment_method: hsp.comp_adjustment_method as i32,
                    edit_script: hsp.gap_info,
                    pat_info: hsp.pat_info.map(Into::into),
                    map_info: hsp.map_info,
                })
                .collect(),
            best_evalue: self.best_evalue,
            hsp_max: self.hsp_max,
        }
    }
}

/// PHI-BLAST pattern metadata stored on an HSP.
///
/// This mirrors the C `hsp->pat_info` fields used by
/// `PHIBlast_HSPResultsSplit` and PHI traceback without preserving the C raw
/// pointer shape.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct PhiPatInfo {
    pub index: usize,
    pub length: i32,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct BlastHSPMappingInfo {
    pub edits: Option<JumperEditsBlock>,
    pub left_edge: u8,
    pub right_edge: u8,
    pub flags: i32,
    pub subject_overhangs: Option<SequenceOverhangs>,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SBlastHitsParameters {
    pub prelim_hitlist_size: i32,
    pub hsp_num_max: i32,
}

#[derive(Debug, Clone, Copy)]
pub struct BlastTargetTranslationView<'a> {
    pub sequence: &'a [u8],
    pub pointer_offset: isize,
    pub translated_length: i32,
}

impl<'a> BlastTargetTranslationView<'a> {
    /// blast-rs: Bounds-checked view accessor for translated target slices; not a direct NCBI C port.
    pub fn get(&self, protein_offset: i32) -> Option<u8> {
        let index = protein_offset as isize + self.pointer_offset;
        if index < 0 {
            return None;
        }
        self.sequence.get(index as usize).copied()
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BlastPartialSubjectTranslation {
    pub translation_buffer: Vec<u8>,
    pub subject: Vec<u8>,
    pub subject_length: i32,
    pub start_shift: i32,
}

/// A list of HSPs for one query-subject pair.
#[derive(Debug, Clone)]
pub struct HspList {
    pub oid: i32,
    pub hsps: Vec<Hsp>,
    pub best_evalue: f64,
    pub hsp_max: i32,
}

impl HspList {
    /// blast-rs: Rust convenience constructor for owned HSP lists; not a direct NCBI C port.
    pub fn new(oid: i32) -> Self {
        HspList {
            oid,
            hsps: Vec::new(),
            // Match `blast_hsp_list_new` (and NCBI `s_BlastGetBestEvalue`'s
            // empty-list result): seed `best_evalue` with `(double)INT4_MAX`
            // so `blast_hsp_list_save_hsp`'s running-min invariant works.
            best_evalue: i32::MAX as f64,
            hsp_max: i32::MAX,
        }
    }

    /// blast-rs: Rust convenience wrapper around `blast_hsp_list_save_hsp`; not a direct NCBI C port.
    pub fn add_hsp(&mut self, hsp: Hsp) {
        let _ = blast_hsp_list_save_hsp(self, hsp);
    }

    /// blast-rs: Rust convenience wrapper around score sorting; not a direct NCBI C port.
    pub fn sort_by_score(&mut self) {
        self.hsps.sort_by(score_compare_hsps);
    }
}

/// blast-rs: Rust ownership equivalent of the C HSP free routine; not a direct NCBI C port.
pub fn blast_hsp_free(hsp: &mut Option<Hsp>) -> Option<Hsp> {
    if let Some(hsp) = hsp.as_mut() {
        hsp.edit_script = crate::gapinfo::gap_edit_script_delete(hsp.edit_script.take());
    }
    *hsp = None;
    None
}

/// NCBI: BlastHSPMappingInfoNew (blast_hits.c:207).
pub fn blast_hsp_mapping_info_new() -> BlastHSPMappingInfo {
    BlastHSPMappingInfo::default()
}

/// blast-rs: Rust ownership equivalent of the C mapping-info free routine; not a direct NCBI C port.
pub fn blast_hsp_mapping_info_free(_: Option<BlastHSPMappingInfo>) -> Option<BlastHSPMappingInfo> {
    None
}

/// NCBI: BlastHspNumMax (blast_hits.c:213).
pub fn blast_hsp_num_max(_gapped_calculation: bool, hsp_num_max: i32) -> i32 {
    if hsp_num_max > 0 {
        hsp_num_max
    } else {
        i32::MAX
    }
}

/// blast-rs: Plain edit-script wrapper around the shared identity counter; not
/// a direct NCBI C port.
pub fn blast_hsp_get_num_identities_plain(
    query: &[u8],
    subject: &[u8],
    hsp: &mut Hsp,
) -> (i16, i32) {
    let ops = hsp.edit_script.as_ref().map(|script| script.ops_vec());
    let (num_ident, align_length, _) =
        crate::blast_kappa::blast_hsp_get_num_identities(query, subject, hsp, ops.as_deref(), None);
    hsp.num_ident = num_ident;
    (0, align_length)
}

/// NCBI: s_Blast_HSPGetOOFNumIdentitiesAndPositives (blast_hits.c:850).
pub fn s_blast_hsp_get_oof_num_identities_and_positives(
    query: &[u8],
    subject: &[u8],
    hsp: &Hsp,
    program: crate::program::ProgramType,
    matrix: Option<&[Vec<i32>]>,
) -> (i16, i32, i32, i32) {
    let Some(edit_script) = hsp.edit_script.as_ref() else {
        return (-1, 0, 0, 0);
    };
    if query.is_empty() || subject.is_empty() {
        return (-1, 0, 0, 0);
    }

    let (q_seq, s_seq, mut q, mut s) =
        if program == crate::program::TBLASTN || program == crate::program::RPS_TBLASTN {
            (
                query,
                subject,
                hsp.query_offset as isize,
                hsp.subject_offset as isize,
            )
        } else {
            (
                subject,
                query,
                hsp.subject_offset as isize,
                hsp.query_offset as isize,
            )
        };

    let mut num_ident = 0;
    let mut num_pos = 0;
    let mut align_length = 0;

    for (op, count) in edit_script.iter() {
        match op {
            GapAlignOpType::Sub => {
                align_length += count;
                for _ in 0..count {
                    let q_res = q_seq.get(q as usize).copied().unwrap_or(u8::MAX);
                    let s_res = s_seq.get(s as usize).copied().unwrap_or(u8::MAX);
                    if q_res == s_res {
                        num_ident += 1;
                    } else if let Some(matrix) = matrix {
                        if matrix
                            .get(q_res as usize)
                            .and_then(|row| row.get(s_res as usize))
                            .copied()
                            .unwrap_or(0)
                            > 0
                        {
                            num_pos += 1;
                        }
                    }
                    q += 1;
                    s += crate::util::CODON_LENGTH as isize;
                }
            }
            GapAlignOpType::Ins => {
                align_length += count;
                s += count as isize * crate::util::CODON_LENGTH as isize;
            }
            GapAlignOpType::Del => {
                align_length += count;
                q += count as isize;
            }
            GapAlignOpType::Del2 => s -= 2,
            GapAlignOpType::Del1 => s -= 1,
            GapAlignOpType::Ins1 => s += 1,
            GapAlignOpType::Ins2 => s += 2,
            _ => {
                s += count as isize * crate::util::CODON_LENGTH as isize;
                q += count as isize;
            }
        }
    }

    let positives = if matrix.is_some() {
        num_pos + num_ident
    } else {
        0
    };
    (0, num_ident, align_length, positives)
}

/// NCBI: Blast_HSPGetTargetTranslation (blast_hits.c:1147).
///
/// C returns `translations[context] - range_start + 1`, which can point before
/// the allocation. Rust returns the allocation plus `pointer_offset` so callers
/// can perform the same coordinate conversion without exposing an invalid
/// slice start.
pub fn blast_hsp_get_target_translation<'a>(
    target_t: &'a mut crate::util::SBlastTargetTranslation,
    hsp: Option<&BlastHSP>,
) -> Option<BlastTargetTranslationView<'a>> {
    let hsp = hsp?;
    let context =
        crate::util::blast_frame_to_context(hsp.subject.frame as i32, target_t.program_number);
    if context < 0 {
        return None;
    }
    let context = context as usize;
    let range = target_t.range.as_mut()?;
    if context >= target_t.translations.len() || 2 * context + 1 >= range.len() {
        return None;
    }

    let mut start = range[2 * context];
    let mut stop = range[2 * context + 1];
    let subject_length = target_t
        .subject_blk
        .as_ref()
        .map(|blk| blk.length.max(0))
        .unwrap_or(0);
    let needs_partial_window = target_t.partial
        && (start != 0 || stop < subject_length / crate::util::CODON_LENGTH as i32 - 3);

    if needs_partial_window {
        const K_MAX_TRANSLATION: i32 = 99;
        let subject_blk = target_t.subject_blk.as_ref()?;
        let subject_seq = subject_blk
            .sequence
            .as_deref()
            .or_else(|| subject_blk.sequence_start.as_deref())?;

        let (nucl_start, mut nucl_end) = if hsp.subject.offset < 0 {
            (0, subject_length)
        } else {
            (
                (3 * hsp.subject.offset - K_MAX_TRANSLATION).max(0),
                subject_length.min(3 * hsp.subject.end + K_MAX_TRANSLATION),
            )
        };
        if subject_length - nucl_end <= 21 {
            nucl_end = subject_length;
        }

        let nucl_length = (nucl_end - nucl_start).max(0);
        let translation_length = 1 + nucl_length / crate::util::CODON_LENGTH as i32;
        let start_shift = nucl_start / crate::util::CODON_LENGTH as i32;
        let nucl_shift = if hsp.subject.frame < 0 {
            subject_length - nucl_start - nucl_length
        } else {
            nucl_start
        }
        .max(0) as usize;
        let nucl_length_usize = nucl_length as usize;
        let nucl_seq = subject_seq.get(nucl_shift..nucl_shift + nucl_length_usize)?;

        if start_shift < start || start_shift + translation_length > stop {
            let old_cached_width = stop - start;
            range[2 * context] = start_shift;
            start = start_shift;
            if translation_length > old_cached_width || target_t.translations[context].is_none() {
                target_t.translations[context] =
                    Some(vec![0; translation_length.max(0) as usize + 2]);
            }

            let rev = if hsp.subject.frame < 0 {
                crate::util::get_reverse_nucl_sequence(nucl_seq, nucl_length_usize)
            } else {
                Vec::new()
            };
            let translation = target_t.translations[context].as_mut()?;
            let length = crate::util::blast_get_translation(
                nucl_seq,
                &rev,
                nucl_length_usize,
                hsp.subject.frame as i32,
                translation,
                &target_t.gen_code_string,
            ) as i32;
            range[2 * context + 1] = start_shift + length;
            stop = range[2 * context + 1];

            if hsp.subject.offset >= 0 {
                translation[0] = crate::util::FENCE_SENTRY;
                if let Some(slot) = translation.get_mut(length as usize + 1) {
                    *slot = crate::util::FENCE_SENTRY;
                }
            }
        }
    }

    let sequence = target_t.translations[context].as_deref()?;
    Some(BlastTargetTranslationView {
        sequence,
        pointer_offset: 1 - start as isize,
        translated_length: stop,
    })
}

fn blast_seg_get_translated_offsets(
    offset: i32,
    end: i32,
    frame: i32,
    seq_length: i32,
) -> (i32, i32) {
    if frame > 0 {
        (
            offset * crate::util::CODON_LENGTH as i32 + frame - 1,
            end * crate::util::CODON_LENGTH as i32 + frame - 2,
        )
    } else if frame < 0 {
        (
            seq_length - offset * crate::util::CODON_LENGTH as i32 + frame,
            seq_length - end * crate::util::CODON_LENGTH as i32 + frame + 1,
        )
    } else {
        (offset + 1, end)
    }
}

/// NCBI: Blast_HSPGetAdjustedOffsets (blast_hits.c:1109).
///
/// Converts internal zero-based half-open HSP coordinates to the one-based
/// endpoint coordinates reported by BLAST, including translated query/subject
/// frames. Rust returns the four outputs as a tuple instead of out-pointers.
pub fn blast_hsp_get_adjusted_offsets(
    program: crate::program::ProgramType,
    hsp: &Hsp,
    query_length: i32,
    subject_length: i32,
) -> (i32, i32, i32, i32) {
    if hsp.edit_script.is_none() {
        return (
            hsp.query_offset + 1,
            hsp.query_end,
            hsp.subject_offset + 1,
            hsp.subject_end,
        );
    }

    if !crate::program::blast_query_is_translated(program)
        && !crate::program::blast_subject_is_translated(program)
    {
        if hsp.query_frame != hsp.subject_frame {
            let q_end = query_length - hsp.query_offset;
            let q_start = q_end - hsp.query_end + hsp.query_offset + 1;
            return (q_start, q_end, hsp.subject_end, hsp.subject_offset + 1);
        }
        return (
            hsp.query_offset + 1,
            hsp.query_end,
            hsp.subject_offset + 1,
            hsp.subject_end,
        );
    }

    let (q_start, q_end) = if crate::program::blast_query_is_translated(program) {
        blast_seg_get_translated_offsets(
            hsp.query_offset,
            hsp.query_end,
            hsp.query_frame,
            query_length,
        )
    } else {
        (hsp.query_offset + 1, hsp.query_end)
    };
    let (s_start, s_end) = if crate::program::blast_subject_is_translated(program) {
        let (_, s_end) = blast_seg_get_translated_offsets(
            hsp.subject_offset,
            hsp.subject_end,
            hsp.subject_frame,
            subject_length,
        );
        (q_start, s_end)
    } else {
        (hsp.subject_offset + 1, hsp.subject_end)
    };

    (q_start, q_end, s_start, s_end)
}

/// NCBI: Blast_HSPGetPartialSubjectTranslation (blast_hits.c:1230).
///
/// Performs the same conservative partial-window translation used by
/// `Blast_HSPGetTargetTranslation`, and shifts the HSP subject protein
/// coordinates so the returned translated subject starts at zero.
pub fn blast_hsp_get_partial_subject_translation(
    subject_blk: &crate::util::BlastSequenceBlk,
    hsp: &mut Hsp,
    is_ooframe: bool,
    gen_code_string: &[u8; 64],
) -> Result<BlastPartialSubjectTranslation, i16> {
    let subject_seq = subject_blk
        .sequence
        .as_deref()
        .or_else(|| subject_blk.sequence_start.as_deref())
        .ok_or(-1_i16)?;
    let subject_length = subject_blk.length.max(0);
    const MAX_FULL_TRANSLATION: i32 = 2100;

    let (nucl_start, nucl_length, start_shift) = if !is_ooframe {
        let start =
            (hsp.subject_offset * crate::util::CODON_LENGTH as i32 - MAX_FULL_TRANSLATION).max(0);
        let length = subject_length
            .min(hsp.subject_end * crate::util::CODON_LENGTH as i32 + MAX_FULL_TRANSLATION)
            - start;
        (
            start,
            length.max(0),
            start / crate::util::CODON_LENGTH as i32,
        )
    } else {
        let start = (hsp.subject_offset - MAX_FULL_TRANSLATION).max(0);
        let length = subject_length.min(hsp.subject_end + MAX_FULL_TRANSLATION) - start;
        (start, length.max(0), start)
    };
    let nucl_shift = if hsp.subject_frame < 0 {
        subject_length - nucl_start - nucl_length
    } else {
        nucl_start
    }
    .max(0) as usize;
    let nucl_length = nucl_length as usize;
    let nucl_seq = subject_seq
        .get(nucl_shift..nucl_shift + nucl_length)
        .ok_or(-1_i16)?;

    let partial = crate::util::blast_get_partial_translation(
        nucl_seq,
        nucl_length,
        hsp.subject_frame,
        gen_code_string,
        !is_ooframe,
        true,
        is_ooframe,
    )
    .map_err(|_| -1_i16)?;

    let (translation_buffer, subject, subject_length) = if is_ooframe {
        let mixed_seq = partial.mixed_seq.ok_or(-1_i16)?;
        let subject_length = partial.protein_length.unwrap_or(mixed_seq.len()) as i32;
        let subject_start = crate::util::CODON_LENGTH.min(mixed_seq.len());
        let subject = mixed_seq[subject_start..].to_vec();
        (Vec::new(), subject, subject_length)
    } else {
        let translation_buffer = partial.translation_buffer.ok_or(-1_i16)?;
        let subject_length = partial.protein_length.unwrap_or(0) as i32;
        let subject_end = 1 + subject_length.max(0) as usize;
        let subject = translation_buffer
            .get(1..subject_end)
            .ok_or(-1_i16)?
            .to_vec();
        (translation_buffer.clone(), subject, subject_length)
    };

    hsp.subject_offset -= start_shift;
    hsp.subject_end -= start_shift;
    hsp.subject_gapped_start -= start_shift;

    Ok(BlastPartialSubjectTranslation {
        translation_buffer,
        subject,
        subject_length,
        start_shift,
    })
}

/// NCBI: s_HSPTest (blast_hits.c:993).
pub fn s_hsp_test(hsp: &BlastHSP, hit_options: &HitSavingOptions, align_length: i32) -> bool {
    (hsp.num_ident as f64 * 100.0 < align_length as f64 * hit_options.percent_identity)
        || align_length < hit_options.min_hit_length
}

/// NCBI: Blast_HSPTest (blast_hits.c:1027).
pub fn blast_hsp_test(hsp: &Hsp, hit_options: &HitSavingOptions, align_length: i32) -> bool {
    let blast_hsp = BlastHSP::from_legacy_hsp(hsp.clone());
    s_hsp_test(&blast_hsp, hit_options, align_length)
}

/// NCBI: Blast_HSPGetQueryCoverage (blast_hits.c:1066).
pub fn blast_hsp_get_query_coverage(hsp: &Hsp, query_length: i32) -> f64 {
    let mut pct = 0.0;
    if query_length > 0 {
        pct = 100.0 * (hsp.query_end - hsp.query_offset) as f64 / query_length as f64;
        if pct < 99.0 {
            pct += 0.5;
        }
    }
    pct
}

/// NCBI: Blast_HSPQueryCoverageTest (blast_hits.c:1077).
pub fn blast_hsp_query_coverage_test(
    hsp: &Hsp,
    min_query_coverage_pct: f64,
    query_length: i32,
) -> bool {
    let hsp_coverage = blast_hsp_get_query_coverage(hsp, query_length);
    hsp_coverage < min_query_coverage_pct
}

/// NCBI: s_HSPStartDiag (blast_hits.c:1464).
pub fn s_hsp_start_diag(hsp: &Hsp) -> i32 {
    hsp.query_offset - hsp.subject_offset
}

/// NCBI: s_HSPEndDiag (blast_hits.c:1474).
pub fn s_hsp_end_diag(hsp: &Hsp) -> i32 {
    hsp.query_end - hsp.subject_end
}

/// NCBI: GetPrelimHitlistSize (blast_hits.c:43).
pub fn get_prelim_hitlist_size(
    hitlist_size: i32,
    composition_based_stats: i32,
    gapped_calculation: bool,
) -> i32 {
    let mut prelim_hitlist_size = hitlist_size;
    if composition_based_stats != 0 {
        if std::env::var_os("ADAPTIVE_CBS").is_some() {
            if hitlist_size < 1000 {
                prelim_hitlist_size = (prelim_hitlist_size + 1000).max(1500);
            } else {
                prelim_hitlist_size = prelim_hitlist_size.saturating_mul(2).saturating_add(50);
            }
        } else if hitlist_size <= 500 {
            prelim_hitlist_size = 1050;
        } else {
            prelim_hitlist_size = prelim_hitlist_size.saturating_mul(2).saturating_add(50);
        }
    } else if gapped_calculation {
        prelim_hitlist_size = prelim_hitlist_size
            .saturating_mul(2)
            .max(10)
            .min(prelim_hitlist_size.saturating_add(50));
    }
    prelim_hitlist_size
}

/// NCBI: SBlastHitsParametersNew (blast_hits.c:74).
/// naming: Historical Rust spelling keeps `sblast` as one token.
pub fn sblast_hits_parameters_new(
    hit_options: Option<&HitSavingOptions>,
    composition_based_stats: i32,
    scoring_options: Option<&ScoringOptions>,
    hsp_num_max: i32,
    retval: &mut Option<SBlastHitsParameters>,
) -> i16 {
    *retval = None;
    let (Some(hit_options), Some(scoring_options)) = (hit_options, scoring_options) else {
        return 1;
    };

    *retval = Some(SBlastHitsParameters {
        prelim_hitlist_size: get_prelim_hitlist_size(
            hit_options.hitlist_size,
            composition_based_stats,
            scoring_options.gapped_calculation,
        ),
        hsp_num_max: blast_hsp_num_max(scoring_options.gapped_calculation, hsp_num_max),
    });
    0
}

/// blast-rs: Rust ownership equivalent of the C hits-parameters free routine; not a direct NCBI C port.
pub fn sblast_hits_parameters_free(
    _: Option<SBlastHitsParameters>,
) -> Option<SBlastHitsParameters> {
    None
}

/// NCBI: Blast_HSPListNew (blast_hits.c:1558).
/// 1-1 translation for the Rust
/// `HspList` shape. `oid` is not initialized by the C allocator; callers
/// fill it later, so this translation uses the neutral sentinel `-1`.
pub fn blast_hsp_list_new(hsp_max: i32) -> HspList {
    let effective_max = if hsp_max > 0 { hsp_max } else { i32::MAX };
    HspList {
        oid: -1,
        hsps: Vec::with_capacity(100usize.min(effective_max as usize)),
        // NCBI `Blast_HSPListNew` uses `calloc` so `best_evalue` starts at 0.0;
        // however our Rust convenience treats `best_evalue` as the running
        // min, so we seed with `s_BlastGetBestEvalue`'s empty-list result
        // `(double)INT4_MAX`. `blast_hsp_list_save_hsp` then mins this against
        // each new HSP's evalue, yielding the same final value NCBI computes
        // explicitly via `s_BlastGetBestEvalue` later.
        best_evalue: i32::MAX as f64,
        hsp_max: effective_max,
    }
}

/// NCBI: Blast_HSPListFree (blast_hits.c:1542).
/// 1-1 translation. Rust drops
/// the list contents; clearing the option preserves the C return value shape.
pub fn blast_hsp_list_free(hsp_list: &mut Option<HspList>) -> Option<HspList> {
    *hsp_list = None;
    None
}

/// NCBI: BlastHSPListDup (blast_hits.c:1583).
/// 1-1 translation. `Clone`
/// performs the same detached copy for the Rust-owned HSP vector and edit
/// scripts.
pub fn blast_hsp_list_dup(hsp_list: Option<&HspList>) -> Option<HspList> {
    let Some(src) = hsp_list else {
        return None;
    };
    let mut dst = HspList::new(src.oid);
    dst.best_evalue = src.best_evalue;
    dst.hsp_max = src.hsp_max;
    dst.hsps.reserve(src.hsps.len());
    for hsp in src.hsps.iter() {
        dst.hsps.push(hsp.clone());
    }
    Some(dst)
}

/// NCBI: Blast_HSPListSwap (blast_hits.c:1614).
pub fn blast_hsp_list_swap(list1: &mut HspList, list2: &mut HspList) {
    std::mem::swap(list1, list2);
}

/// NCBI: Blast_HSPListSortByScore (blast_hits.c:1374).
pub fn blast_hsp_list_sort_by_score(hsp_list: Option<&mut BlastHSPList>) -> i32 {
    let Some(hsp_list) = hsp_list else {
        return 0;
    };
    if hsp_list.hsp_array.len() <= 1 {
        hsp_list.hspcnt = hsp_list.hsp_array.iter().filter(|hsp| hsp.is_some()).count() as i32;
        return 0;
    }
    hsp_list
        .hsp_array
        .sort_by(|a, b| score_compare_blast_hsps(a.as_ref(), b.as_ref()));
    hsp_list.hspcnt = hsp_list.hsp_array.iter().filter(|hsp| hsp.is_some()).count() as i32;
    0
}

/// NCBI: Blast_HSPUpdateWithTraceback (blast_traceback.c:78).
pub fn blast_hsp_update_with_traceback(
    gap_align: Option<&mut crate::blast_kappa::BlastGapAlignWorkspace>,
    hsp: Option<&mut Hsp>,
) -> i16 {
    let (Some(gap_align), Some(hsp)) = (gap_align, hsp) else {
        return -1;
    };

    hsp.score = gap_align.score;
    hsp.query_offset = gap_align.query_start;
    hsp.subject_offset = gap_align.subject_start;
    hsp.query_end = gap_align.query_stop;
    hsp.subject_end = gap_align.subject_stop;
    if gap_align.edit_script.is_some() {
        hsp.edit_script = gap_align.edit_script.take();
    }

    0
}

/// NCBI: Blast_HSPListSortByEvalue (blast_hits.c:1437).
///
/// NCBI checks `s_EvalueCompareHSPs` for an existing sort first and only
/// calls `qsort` if the list is not yet ordered. We mirror that.
pub fn blast_hsp_list_sort_by_evalue(hsp_list: Option<&mut HspList>) -> i32 {
    let Some(hsp_list) = hsp_list else {
        return 0;
    };
    if hsp_list.hsps.len() <= 1 {
        return 0;
    }
    let already_sorted = hsp_list
        .hsps
        .windows(2)
        .all(|pair| !evalue_compare_hsps(&pair[0], &pair[1]).is_gt());
    if already_sorted {
        return 0;
    }
    hsp_list.hsps.sort_by(evalue_compare_hsps);
    0
}

/// NCBI: s_Heapify (blast_hits.c:1627).
/// naming: Typed Rust helper keeps the C helper name while generalizing over element type.
///
/// C receives raw byte pointers plus element width; Rust receives the same
/// logical heap slice and element indices. The comparison polarity is the
/// same: elements that compare `Greater` are moved toward the heap root.
pub fn s_heapify<T>(
    heap: &mut [T],
    mut base: usize,
    lim: usize,
    last: usize,
    compar: impl Fn(&T, &T) -> std::cmp::Ordering,
) {
    if heap.is_empty() {
        return;
    }

    let mut left_son = 2 * base + 1;
    while base <= lim && left_son <= last {
        let large_son = if left_son == last {
            left_son
        } else if compar(&heap[left_son], &heap[left_son + 1]).is_ge() {
            left_son
        } else {
            left_son + 1
        };

        if compar(&heap[base], &heap[large_son]).is_lt() {
            heap.swap(base, large_son);
            base = large_son;
            left_son = 2 * base + 1;
        } else {
            break;
        }
    }
}

/// NCBI: s_CreateHeap (blast_hits.c:1660).
/// naming: Typed Rust helper keeps the C helper name while generalizing over element type.
pub fn s_create_heap<T>(heap: &mut [T], compar: impl Fn(&T, &T) -> std::cmp::Ordering + Copy) {
    let nel = heap.len();
    if nel <= 1 {
        return;
    }

    let last = nel - 1;
    let mut base = nel / 2;
    while base > 0 {
        base -= 1;
        let lim = nel / 2 - 1;
        s_heapify(heap, base, lim, last, compar);
    }
}

/// NCBI: s_BlastHSPListInsertHSPInHeap (blast_hits.c:1687).
/// naming: Rust spells HSP as a separate snake_case token.
/// Port for the Rust vector-backed HSP list. NCBI does
/// NOT touch `best_evalue` here — callers are responsible for setting it
/// when needed (typically via `Blast_HSPListGetEvalues` later).
pub fn s_blast_hsp_list_insert_hsp_in_heap(hsp_list: &mut HspList, hsp: Hsp) {
    if hsp_list.hsps.is_empty() {
        // Defensive: NCBI never calls this on an empty list (`Blast_HSPListSaveHSP`
        // only reaches the heap branch when full); track best_evalue here so
        // our Rust running-min invariant stays consistent.
        hsp_list.best_evalue = hsp_list.best_evalue.min(hsp.evalue);
        hsp_list.hsps.push(hsp);
        return;
    }

    if score_compare_hsps(&hsp, &hsp_list.hsps[0]).is_gt() {
        return;
    }

    hsp_list.hsps[0] = hsp;
    if hsp_list.hsps.len() >= 2 {
        let lim = hsp_list.hsps.len() / 2 - 1;
        let last = hsp_list.hsps.len() - 1;
        s_heapify(&mut hsp_list.hsps, 0, lim, last, score_compare_hsps);
    }
    // NCBI `s_BlastHSPListInsertHSPInHeap` (`blast_hits.c:1687`) does not
    // recompute `best_evalue`. Rust running-min invariant: best_evalue is
    // still valid (the replaced HSP at slot 0 might have been the min, but
    // since NCBI's design has best_evalue stale until explicit recompute,
    // we maintain the same level of staleness for parity).
}

/// NCBI: Blast_HSPListSaveHSP (blast_hits.c:1754).
/// naming: Rust spells HSP as a separate snake_case token.
/// 1-1 translation for the
/// vector-backed Rust list. When the list has reached `hsp_max`, the new HSP
/// replaces the current worst HSP only if it sorts better by
/// [`score_compare_hsps`].
pub fn blast_hsp_list_save_hsp(hsp_list: &mut HspList, hsp: Hsp) -> i32 {
    if hsp_list.hsps.len() < hsp_list.hsp_max.max(0) as usize {
        hsp_list.best_evalue = hsp_list.best_evalue.min(hsp.evalue);
        hsp_list.hsps.push(hsp);
        return 0;
    }

    if hsp_list.hsps.is_empty() {
        return 0;
    }

    s_create_heap(&mut hsp_list.hsps, score_compare_hsps);
    s_blast_hsp_list_insert_hsp_in_heap(hsp_list, hsp);
    0
}

/// blast-rs: Local interval containment helper extracted from HSP merge logic; not a direct NCBI C port.
fn contained_in_hsp(a: i32, b: i32, c: i32, d: i32, e: i32, f: i32) -> bool {
    a <= c && b >= c && d <= f && e >= f
}

/// NCBI: s_BlastMergeTwoHSPs (blast_hits.c:2672).
/// naming: Rust spells HSPs as a readable snake_case plural token.
fn s_blast_merge_two_hsps(hsp1: &mut Hsp, hsp2: &Hsp, allow_gap: bool) -> bool {
    debug_assert!(hsp1.edit_script.is_none() || hsp2.edit_script.is_none());

    if !allow_gap
        && hsp1.subject_offset - hsp2.subject_offset - hsp1.query_offset + hsp2.query_offset != 0
    {
        return false;
    }

    if hsp1.subject_frame != hsp2.subject_frame {
        return false;
    }

    if contained_in_hsp(
        hsp1.query_offset,
        hsp1.query_end,
        hsp2.query_offset,
        hsp1.subject_offset,
        hsp1.subject_end,
        hsp2.subject_offset,
    ) || contained_in_hsp(
        hsp1.query_offset,
        hsp1.query_end,
        hsp2.query_end,
        hsp1.subject_offset,
        hsp1.subject_end,
        hsp2.subject_end,
    ) {
        let hsp1_query_len = hsp1.query_end - hsp1.query_offset;
        let hsp2_query_len = hsp2.query_end - hsp2.query_offset;
        let denom = hsp1_query_len + hsp2_query_len;
        let score_density = if denom != 0 {
            (hsp1.score + hsp2.score) as f64 / denom as f64
        } else {
            0.0
        };

        hsp1.query_offset = hsp1.query_offset.min(hsp2.query_offset);
        hsp1.subject_offset = hsp1.subject_offset.min(hsp2.subject_offset);
        hsp1.query_end = hsp1.query_end.max(hsp2.query_end);
        hsp1.subject_end = hsp1.subject_end.max(hsp2.subject_end);
        if hsp2.score > hsp1.score {
            hsp1.query_gapped_start = hsp2.query_gapped_start;
            hsp1.subject_gapped_start = hsp2.subject_gapped_start;
            hsp1.score = hsp2.score;
        }

        let merged_len = hsp1.query_end - hsp1.query_offset;
        hsp1.score = ((score_density * merged_len as f64) as i32).max(hsp1.score);
        return true;
    }

    false
}

const OVERLAP_DIAG_CLOSE: i32 = 10;

/// blast-rs: Split-query context index helper for HSP-list merge; not a direct NCBI C port.
fn hsp_context_offset_index(context: i32, contexts_per_query: i32) -> Option<usize> {
    if contexts_per_query <= 0 {
        return None;
    }
    Some(context.rem_euclid(contexts_per_query) as usize)
}

/// NCBI: Blast_HSPListsMerge (blast_hits.c:2857).
/// naming: Rust spells HSP lists as separate snake_case tokens.
/// 1-1 translation for
/// Rust-owned HSP vectors. The split-offset arguments preserve C's overlap
/// strip merge for subject chunks (`contexts_per_query < 0`) and query chunks.
pub fn blast_hsp_lists_merge(
    hsp_list: &mut Option<HspList>,
    combined_hsp_list: &mut Option<HspList>,
    hsp_num_max: i32,
    split_offsets: Option<&[i32]>,
    contexts_per_query: i32,
    chunk_overlap_size: i32,
    allow_gap: bool,
    short_reads: bool,
) -> i32 {
    let Some(mut hsp_list_value) = hsp_list.take() else {
        return 0;
    };
    if hsp_list_value.hsps.is_empty() {
        *hsp_list = Some(hsp_list_value);
        return 0;
    }

    let Some(combined) = combined_hsp_list.as_mut() else {
        *combined_hsp_list = Some(hsp_list_value);
        return 0;
    };

    let mut hspcnt1 = 0usize;
    let mut hspcnt2 = 0usize;

    if let Some(split_offsets) = split_offsets {
        if contexts_per_query < 0 {
            let split = split_offsets.first().copied().unwrap_or(i32::MAX);
            for index1 in 0..combined.hsps.len() {
                if combined.hsps[index1].subject_end > split {
                    combined.hsps.swap(hspcnt1, index1);
                    hspcnt1 += 1;
                }
            }
            for index2 in 0..hsp_list_value.hsps.len() {
                if hsp_list_value.hsps[index2].subject_offset < split + chunk_overlap_size {
                    hsp_list_value.hsps.swap(hspcnt2, index2);
                    hspcnt2 += 1;
                }
            }
        } else {
            for index1 in 0..combined.hsps.len() {
                let Some(offset_idx) =
                    hsp_context_offset_index(combined.hsps[index1].context, contexts_per_query)
                else {
                    continue;
                };
                let Some(&split) = split_offsets.get(offset_idx) else {
                    continue;
                };
                if split < 0 {
                    continue;
                }
                if (combined.hsps[index1].query_frame >= 0
                    && combined.hsps[index1].query_end > split)
                    || (combined.hsps[index1].query_frame < 0
                        && combined.hsps[index1].query_offset < split + chunk_overlap_size)
                {
                    combined.hsps.swap(hspcnt1, index1);
                    hspcnt1 += 1;
                }
            }
            for index2 in 0..hsp_list_value.hsps.len() {
                let Some(offset_idx) = hsp_context_offset_index(
                    hsp_list_value.hsps[index2].context,
                    contexts_per_query,
                ) else {
                    continue;
                };
                let Some(&split) = split_offsets.get(offset_idx) else {
                    continue;
                };
                if split < 0 {
                    continue;
                }
                if (hsp_list_value.hsps[index2].query_frame < 0
                    && hsp_list_value.hsps[index2].query_end > split)
                    || (hsp_list_value.hsps[index2].query_frame >= 0
                        && hsp_list_value.hsps[index2].query_offset < split + chunk_overlap_size)
                {
                    hsp_list_value.hsps.swap(hspcnt2, index2);
                    hspcnt2 += 1;
                }
            }
        }
    }

    if hspcnt1 > 0 && hspcnt2 > 0 {
        let mut delete_hsp2 = vec![false; hsp_list_value.hsps.len()];
        for index1 in 0..hspcnt1 {
            for index2 in 0..hspcnt2 {
                if delete_hsp2[index2] {
                    continue;
                }
                if combined.hsps[index1].context != hsp_list_value.hsps[index2].context {
                    continue;
                }
                if short_reads {
                    delete_hsp2[index2] = true;
                    continue;
                }

                let (end_diag, start_diag) =
                    if contexts_per_query < 0 || combined.hsps[index1].query_frame >= 0 {
                        (
                            s_hsp_end_diag(&combined.hsps[index1]),
                            s_hsp_start_diag(&hsp_list_value.hsps[index2]),
                        )
                    } else {
                        (
                            s_hsp_end_diag(&hsp_list_value.hsps[index2]),
                            s_hsp_start_diag(&combined.hsps[index1]),
                        )
                    };

                if (end_diag - start_diag).abs() < OVERLAP_DIAG_CLOSE {
                    let hsp2 = hsp_list_value.hsps[index2].clone();
                    if s_blast_merge_two_hsps(&mut combined.hsps[index1], &hsp2, allow_gap) {
                        delete_hsp2[index2] = true;
                    }
                }
            }
        }

        let mut index = 0usize;
        hsp_list_value.hsps.retain(|_| {
            let keep = !delete_hsp2[index];
            index += 1;
            keep
        });
        hsp_list_value.best_evalue = s_blast_get_best_evalue(&hsp_list_value);
    }

    let new_hspcnt =
        (hsp_list_value.hsps.len() + combined.hsps.len()).min(hsp_num_max.max(0) as usize);
    s_blast_hsp_lists_combine_by_score(&mut hsp_list_value, combined, new_hspcnt as i32);
    0
}

/// blast-rs: Convenience wrapper around `blast_hsp_lists_merge`; not a direct NCBI C port.
pub fn blast_hsp_lists_merge_simple(
    hsp_list: &mut Option<HspList>,
    combined_hsp_list: &mut Option<HspList>,
    hsp_num_max: i32,
) -> i32 {
    blast_hsp_lists_merge(
        hsp_list,
        combined_hsp_list,
        hsp_num_max,
        None,
        0,
        0,
        true,
        false,
    )
}

/// NCBI: Blast_HSPListPHIGetBitScores (blast_hits.c:1934).
/// naming: Rust spells HSP/PHI as separate snake_case tokens.
/// NCBI reads `lambda` and `paramC` from
/// `sbp->kbp_gap[0]`; the Rust translation takes those score-block values
/// explicitly and updates every HSP in place.
pub fn blast_hsp_list_phi_get_bit_scores(hsp_list: &mut HspList, lambda: f64, param_c: f64) {
    let log_c = param_c.ln();
    for hsp in &mut hsp_list.hsps {
        hsp.bit_score =
            (hsp.score as f64 * lambda - log_c - (1.0 + hsp.score as f64 * lambda).ln())
                / crate::math::NCBIMATH_LN2;
    }
}

/// NCBI: PhiBlastGetEffectiveNumberOfPatterns (blast_hits.c:360).
///
/// C reads `occurrences` from `query_info->pattern_info` and reads the
/// minimum pattern length from `query_info->contexts[0].length_adjustment`.
/// Rust keeps those two inputs explicit so the counting rule remains local and
/// auditable.
pub fn phi_blast_get_effective_number_of_patterns(
    occurrence_offsets: &[i32],
    min_pattern_length: i32,
) -> i32 {
    if occurrence_offsets.len() <= 1 {
        return occurrence_offsets.len() as i32;
    }

    let mut count = 1;
    let mut last_effective_occurrence = occurrence_offsets[0];
    for &offset in occurrence_offsets.iter().skip(1) {
        if (offset - last_effective_occurrence) * 2 > min_pattern_length {
            last_effective_occurrence = offset;
            count += 1;
        }
    }
    count
}

/// NCBI: s_HSPPHIGetEvalue (blast_hits.c:399).
/// naming: Rust spells HSP/PHI as separate snake_case tokens.
///
/// The C helper reads `paramC`/`Lambda` from `sbp->kbp[0]`, the effective
/// query-pattern count from `PhiBlastGetEffectiveNumberOfPatterns`, and
/// `num_patterns_db` from `SPHIPatternSearchBlk`. Rust passes those values
/// explicitly while preserving the same formula.
pub fn s_hsp_phi_get_evalue(
    hsp: &mut Hsp,
    lambda: f64,
    param_c: f64,
    effective_num_patterns: i32,
    num_patterns_db: i32,
) {
    hsp.evalue = param_c
        * (1.0 + lambda * hsp.score as f64)
        * effective_num_patterns as f64
        * num_patterns_db as f64
        * (-lambda * hsp.score as f64).exp();
}

/// NCBI: Blast_HSPListPHIGetEvalues (blast_hits.c:1955).
/// naming: Rust spells HSP/PHI as separate snake_case tokens.
/// Port with the
/// C score-block and pattern-block fields passed explicitly.
pub fn blast_hsp_list_phi_get_evalues(
    hsp_list: &mut HspList,
    lambda: f64,
    param_c: f64,
    occurrence_offsets: &[i32],
    min_pattern_length: i32,
    num_patterns_db: i32,
) {
    if hsp_list.hsps.is_empty() {
        return;
    }
    let effective_num_patterns =
        phi_blast_get_effective_number_of_patterns(occurrence_offsets, min_pattern_length);
    for hsp in &mut hsp_list.hsps {
        s_hsp_phi_get_evalue(
            hsp,
            lambda,
            param_c,
            effective_num_patterns,
            num_patterns_db,
        );
    }
    hsp_list.best_evalue = s_blast_get_best_evalue(hsp_list);
}

/// NCBI: Blast_HSPListIsSortedByScore (blast_hits.c:1358).
pub fn blast_hsp_list_is_sorted_by_score(hsp_list: Option<&HspList>) -> bool {
    let Some(hsp_list) = hsp_list else {
        return true;
    };
    hsp_list
        .hsps
        .windows(2)
        .all(|pair| !score_compare_hsps(&pair[0], &pair[1]).is_gt())
}

/// NCBI: s_HSPListRescaleScores (blast_traceback.c:106).
///
/// Traceback scores may be held in a scaled integer space. C removes that
/// scaling with truncating integer division after adding half the divisor,
/// then sorts again because score ties can appear after rescaling.
pub fn s_hsp_list_rescale_scores(hsp_list: &mut HspList, scale_factor: f64) {
    for hsp in &mut hsp_list.hsps {
        hsp.score = ((hsp.score as f64 + 0.5 * scale_factor) / scale_factor) as i32;
    }
    hsp_list.hsps.sort_by(score_compare_hsps);
}

/// NCBI: s_BlastHSPRPSUpdate (blast_traceback.c:131).
/// naming: Rust spells HSP/RPS as separate snake_case tokens.
///
/// RPS traceback temporarily swaps query and subject; this helper flips
/// insertion/deletion edit operations so the script again describes the
/// restored query/subject orientation.
pub fn s_blast_hsp_rps_update(hsp: &mut Hsp) {
    let Some(edit_script) = hsp.edit_script.as_mut() else {
        return;
    };
    for index in 0..edit_script.len() {
        let (op, count) = edit_script.get(index).unwrap_or((GapAlignOpType::Sub, 0));
        let new_op = match op {
            GapAlignOpType::Ins => GapAlignOpType::Del,
            GapAlignOpType::Del => GapAlignOpType::Ins,
            other => other,
        };
        edit_script.set(index, new_op, count);
    }
}

/// NCBI: s_BlastHSPListRPSUpdate (blast_traceback.c:155).
///
/// For RPS programs, the traceback code has query/subject roles reversed.
/// This restores the HSP coordinates and frames, fixes the edit script, and
/// for RPS-tblastn derives the query context from the restored query frame.
pub fn s_blast_hsp_list_rps_update(program: crate::program::ProgramType, hsp_list: &mut HspList) {
    if !crate::program::blast_program_is_rps_blast(program) {
        return;
    }

    for hsp in &mut hsp_list.hsps {
        std::mem::swap(&mut hsp.query_offset, &mut hsp.subject_offset);
        std::mem::swap(&mut hsp.query_end, &mut hsp.subject_end);
        std::mem::swap(&mut hsp.query_gapped_start, &mut hsp.subject_gapped_start);
        std::mem::swap(&mut hsp.query_frame, &mut hsp.subject_frame);

        s_blast_hsp_rps_update(hsp);

        if program == crate::program::RPS_TBLASTN {
            hsp.context = crate::util::blast_frame_to_context(hsp.query_frame, program);
        }
    }

    hsp_list.hsps.sort_by(score_compare_hsps);
}

/// Native RPS profile database magic from NCBI `blast_rps.h`.
///
/// C uses this to distinguish the 26-column on-disk profile format from the
/// newer 28-column format.
pub const RPS_MAGIC_NUM: i32 = 0x1e16;
pub const RPS_MAGIC_NUM_28: i32 = 0x1e17;
pub const RPS_K_MULT: f64 = 1.2;
pub const RPS_FREQ_RATIO_SCALE: f64 = 1_000_000_000.0;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RpsProfileHeader {
    pub magic_number: i32,
    pub num_profiles: i32,
    pub start_offsets: Vec<i32>,
    pub pssm_values: Vec<i32>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RpsFreqRatiosHeader {
    pub start_offsets: Vec<i32>,
    pub freq_values: Vec<i32>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct RpsAuxInfo {
    pub matrix: String,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub ungapped_k: f64,
    pub ungapped_h: f64,
    pub max_db_seq_length: i32,
    pub db_length: i32,
    pub scale_factor: f64,
    pub profile_lengths: Vec<i32>,
    pub karlin_k: Vec<f64>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct RpsTracebackInfo {
    pub profile_header: RpsProfileHeader,
    pub freq_ratios_header: Option<RpsFreqRatiosHeader>,
    pub karlin_k: Vec<f64>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RpsNativeByteOrder {
    LittleEndian,
    BigEndian,
}

fn rps_native_magic_and_byte_order(bytes: &[u8]) -> Option<(i32, RpsNativeByteOrder)> {
    let magic_bytes: [u8; 4] = bytes.get(..4)?.try_into().ok()?;
    let little = i32::from_le_bytes(magic_bytes);
    if rps_alphabet_size(little).is_some() {
        return Some((little, RpsNativeByteOrder::LittleEndian));
    }
    let big = i32::from_be_bytes(magic_bytes);
    if rps_alphabet_size(big).is_some() {
        return Some((big, RpsNativeByteOrder::BigEndian));
    }
    None
}

fn rps_native_i32_values(bytes: &[u8], byte_order: RpsNativeByteOrder) -> Option<Vec<i32>> {
    if bytes.len() % std::mem::size_of::<i32>() != 0 {
        return None;
    }
    Some(
        bytes
            .chunks_exact(std::mem::size_of::<i32>())
            .map(|chunk| {
                let raw: [u8; 4] = chunk.try_into().expect("chunks_exact width");
                match byte_order {
                    RpsNativeByteOrder::LittleEndian => i32::from_le_bytes(raw),
                    RpsNativeByteOrder::BigEndian => i32::from_be_bytes(raw),
                }
            })
            .collect(),
    )
}

fn rps_native_profile_values(
    values: &[i32],
    rows: usize,
    alphabet_size: usize,
) -> Option<Vec<i32>> {
    let required_values = rows.checked_mul(alphabet_size)?;
    let mut pssm_values = values.get(..required_values)?.to_vec();

    // NCBI's pointer graph allocates one extra row pointer, but current native
    // `.rps`/`.freq` payloads carry rows through the final start offset. Keep
    // the owned traceback path's sentinel-row shape by appending a zero row
    // when the file does not contain one.
    let sentinel_values = alphabet_size;
    if values.len() >= required_values.checked_add(sentinel_values)? {
        pssm_values
            .extend_from_slice(values.get(required_values..required_values + sentinel_values)?);
    } else {
        pssm_values.extend(std::iter::repeat(0).take(sentinel_values));
    }
    Some(pssm_values)
}

/// blast-rs: native `.rps` profile-header parser for the owned RPS boundary;
/// not a direct NCBI C port.
///
/// The file layout is NCBI `BlastRPSProfileHeader`: magic, profile count,
/// `num_profiles + 1` start offsets, then flat `Int4` PSSM rows. Both little
/// and big-endian native files are accepted. The returned metadata keeps the
/// existing owned traceback invariant by carrying an extra sentinel row.
pub fn rps_profile_header_from_native_bytes(bytes: &[u8]) -> Result<RpsProfileHeader, i16> {
    let (magic_number, byte_order) = rps_native_magic_and_byte_order(bytes).ok_or(-1i16)?;
    let alphabet_size = rps_alphabet_size(magic_number).ok_or(-1i16)?;
    let values = rps_native_i32_values(bytes, byte_order).ok_or(-1i16)?;
    if values.len() < 3 {
        return Err(-1);
    }
    let num_profiles = values[1];
    if num_profiles < 0 {
        return Err(-1);
    }
    let offset_count = num_profiles as usize + 1;
    let offset_start = 2usize;
    let offset_end = offset_start.checked_add(offset_count).ok_or(-1i16)?;
    let start_offsets = values.get(offset_start..offset_end).ok_or(-1i16)?.to_vec();
    let num_rows = start_offsets.last().copied().ok_or(-1i16)?;
    if num_rows < 0 {
        return Err(-1);
    }
    let pssm_values = rps_native_profile_values(
        values.get(offset_end..).ok_or(-1i16)?,
        num_rows as usize,
        alphabet_size,
    )
    .ok_or(-1i16)?;

    Ok(RpsProfileHeader {
        magic_number,
        num_profiles,
        start_offsets,
        pssm_values,
    })
}

/// blast-rs: native `.freq` frequency-ratio parser for the owned RPS boundary;
/// not a direct NCBI C port.
pub fn rps_freq_ratios_header_from_native_bytes(bytes: &[u8]) -> Result<RpsFreqRatiosHeader, i16> {
    let (magic_number, byte_order) = rps_native_magic_and_byte_order(bytes).ok_or(-1i16)?;
    let alphabet_size = rps_alphabet_size(magic_number).ok_or(-1i16)?;
    let values = rps_native_i32_values(bytes, byte_order).ok_or(-1i16)?;
    if values.len() < 3 {
        return Err(-1);
    }
    let num_profiles = values[1];
    if num_profiles < 0 {
        return Err(-1);
    }
    let offset_count = num_profiles as usize + 1;
    let offset_start = 2usize;
    let offset_end = offset_start.checked_add(offset_count).ok_or(-1i16)?;
    let start_offsets = values.get(offset_start..offset_end).ok_or(-1i16)?.to_vec();
    let num_rows = start_offsets.last().copied().ok_or(-1i16)?;
    if num_rows < 0 {
        return Err(-1);
    }
    let freq_values = rps_native_profile_values(
        values.get(offset_end..).ok_or(-1i16)?,
        num_rows as usize,
        alphabet_size,
    )
    .ok_or(-1i16)?;

    Ok(RpsFreqRatiosHeader {
        start_offsets,
        freq_values,
    })
}

/// blast-rs: text `.aux` metadata parser for the owned RPS boundary;
/// not a direct NCBI C port.
///
/// NCBI `CRpsAuxFile::x_ReadFromFile` reads the scoring matrix, gap costs,
/// ungapped K/H, maximum DB sequence length, total DB length, scaling factor,
/// then `(profile_length, karlin_k)` pairs through EOF. Rust keeps those fields
/// as owned metadata and rejects incomplete or non-finite records.
pub fn rps_aux_info_from_bytes(bytes: &[u8]) -> Result<RpsAuxInfo, i16> {
    let text = std::str::from_utf8(bytes).map_err(|_| -1i16)?;
    let mut tokens = text.split_whitespace();
    let matrix = tokens.next().ok_or(-1i16)?.to_string();
    if matrix.is_empty() {
        return Err(-1);
    }

    let parse_i32 = |token: Option<&str>| -> Result<i32, i16> {
        token.ok_or(-1i16)?.parse::<i32>().map_err(|_| -1i16)
    };
    let parse_f64 = |token: Option<&str>| -> Result<f64, i16> {
        let value = token.ok_or(-1i16)?.parse::<f64>().map_err(|_| -1i16)?;
        if value.is_finite() {
            Ok(value)
        } else {
            Err(-1)
        }
    };

    let gap_open = parse_i32(tokens.next())?;
    let gap_extend = parse_i32(tokens.next())?;
    let ungapped_k = parse_f64(tokens.next())?;
    let ungapped_h = parse_f64(tokens.next())?;
    let max_db_seq_length = parse_i32(tokens.next())?;
    let db_length = parse_i32(tokens.next())?;
    let scale_factor = parse_f64(tokens.next())?;

    let mut profile_lengths = Vec::new();
    let mut karlin_k = Vec::new();
    while let Some(length_token) = tokens.next() {
        let profile_length = length_token.parse::<i32>().map_err(|_| -1i16)?;
        let k = parse_f64(tokens.next())?;
        if profile_length < 0 || k < 0.0 {
            return Err(-1);
        }
        profile_lengths.push(profile_length);
        karlin_k.push(k);
    }
    if karlin_k.is_empty() {
        return Err(-1);
    }

    Ok(RpsAuxInfo {
        matrix,
        gap_open,
        gap_extend,
        ungapped_k,
        ungapped_h,
        max_db_seq_length,
        db_length,
        scale_factor,
        profile_lengths,
        karlin_k,
    })
}

#[derive(Debug, Clone)]
pub struct RpsGapAlignData {
    pub concat_db_info: crate::queryinfo::QueryInfo,
    pub rps_pssm: Vec<Vec<i32>>,
    pub rps_freq: Option<Vec<Vec<i32>>>,
    pub rps_freq_ratios: Option<Vec<Vec<f64>>>,
    pub alphabet_size: usize,
    pub position_based: bool,
}

/// blast-rs: RPS profile magic-number decoder for owned traceback metadata; not a direct NCBI C port.
fn rps_alphabet_size(magic_number: i32) -> Option<usize> {
    match magic_number {
        RPS_MAGIC_NUM => Some(26),
        RPS_MAGIC_NUM_28 => Some(28),
        _ => None,
    }
}

/// blast-rs: RPS profile row decoder for owned integer tables; not a direct NCBI C port.
fn rps_rows_from_values(
    values: &[i32],
    num_rows: usize,
    alphabet_size: usize,
) -> Option<Vec<Vec<i32>>> {
    let total = num_rows.checked_add(1)?.checked_mul(alphabet_size)?;
    let rows = values.get(..total)?;
    Some(
        rows.chunks_exact(alphabet_size)
            .map(|row| row.to_vec())
            .collect(),
    )
}

/// NCBI: s_RPSFillFreqRatiosInPsiMatrix (blast_traceback.c:1001).
/// naming: Rust exposes the RPS acronym as a separate snake_case token.
///
/// C allocates `ncol x BLASTAA_SIZE` doubles, copies the on-disk integer
/// ratios divided by `FREQ_RATIO_SCALE` through the profile alphabet size, and
/// pads the remaining BLASTAA columns with zero.
pub fn s_rps_fill_freq_ratios_in_psi_matrix(
    freq: &[Vec<i32>],
    ncol: usize,
    alphabet_size: usize,
) -> Option<Vec<Vec<f64>>> {
    if alphabet_size > crate::encoding::BLASTAA_SIZE || freq.len() < ncol {
        return None;
    }
    let mut ratios = vec![vec![0.0; crate::encoding::BLASTAA_SIZE]; ncol];
    for (ic, row) in freq.iter().take(ncol).enumerate() {
        if row.len() < alphabet_size {
            return None;
        }
        for ir in 0..alphabet_size {
            ratios[ic][ir] = row[ir] as f64 / RPS_FREQ_RATIO_SCALE;
        }
    }
    Some(ratios)
}

/// blast-rs: Port-shaped fetch adapter for the `BlastSeqSrcGetSequence` call inside
/// `s_RPSComputeTraceback` (`blast_traceback.c:1182`).
/// not a direct NCBI C port.
///
/// The C loop chooses the traceback encoding from the program, fetches the RPS
/// consensus sequence for the current HSP-list OID, and skips that list if the
/// fetch fails. Rust returns `None` for the same skip condition.
pub fn s_rps_fetch_consensus_sequence(
    program_number: crate::program::ProgramType,
    seq_src: Option<&dyn crate::seqsrc::BlastSeqSource>,
    oid: i32,
) -> Option<crate::seqsrc::SeqData> {
    let mut arg = crate::seqsrc::BlastSeqSrcGetSeqArg {
        oid,
        encoding: blast_traceback_get_encoding(program_number).into(),
        check_oid_exclusion: true,
        ..crate::seqsrc::BlastSeqSrcGetSeqArg::default()
    };
    crate::seqsrc::blast_seq_src_get_sequence(seq_src, Some(&mut arg))
}

/// blast-rs: Port-shaped profile PSSM setup adapter for `s_RPSComputeTraceback`
/// (`blast_traceback.c:1196-1248`).
/// not a direct NCBI C port.
///
/// For CBS, C allocates a fresh `seq_length x BLASTAA_SIZE` PSSM scratch matrix
/// before redo alignment fills it. Without CBS, RPS-tblastn uses the selected
/// profile rows directly, while RPS-BLAST calls `RPSRescalePssm`.
pub fn s_rps_profile_pssm_for_traceback(
    program_number: crate::program::ProgramType,
    gap_data: &RpsGapAlignData,
    profile_index: usize,
    query_sequence: &[u8],
    consensus_length: usize,
    scale_factor: f64,
    sbp: &crate::stat::BlastScoreBlk,
    composition_based_stats: bool,
) -> Option<Vec<Vec<i32>>> {
    let profile_context = gap_data.concat_db_info.contexts.get(profile_index)?;
    if !profile_context.is_valid {
        return None;
    }
    let db_seq_start = profile_context.query_offset.max(0) as usize;
    let db_seq_end = db_seq_start.checked_add(consensus_length)?;
    if db_seq_end > gap_data.rps_pssm.len() {
        return None;
    }
    if composition_based_stats {
        return Some(vec![
            vec![0; crate::encoding::BLASTAA_SIZE];
            consensus_length
        ]);
    }
    if gap_data.alphabet_size > crate::encoding::BLASTAA_SIZE
        || !gap_data.rps_pssm[db_seq_start..db_seq_end]
            .iter()
            .all(|row| row.len() >= gap_data.alphabet_size)
    {
        return None;
    }
    if program_number == crate::program::RPS_TBLASTN {
        return Some(gap_data.rps_pssm[db_seq_start..db_seq_end].to_vec());
    }
    crate::stat::rps_rescale_pssm(
        scale_factor,
        query_sequence.len() as i32,
        Some(query_sequence),
        consensus_length as i32,
        Some(&gap_data.rps_pssm[db_seq_start..db_seq_end]),
        Some(sbp),
    )
}

/// blast-rs: Port-shaped Rust equivalent of NCBI `s_RPSGapAlignDataPrepare`
/// (`blast_traceback.c:944`).
/// not a direct NCBI C port.
///
/// C builds a `BlastQueryInfo` over the concatenated profile database and
/// attaches row pointers into the mmapped RPS PSSM/frequency-ratio blocks. Rust
/// keeps the same offsets and row layout in owned vectors so the traceback
/// driver can audit profile selection without pretending to own the full C
/// score block.
pub fn s_rps_gap_align_data_prepare(
    rps_info: Option<&RpsTracebackInfo>,
) -> Result<RpsGapAlignData, i16> {
    let Some(rps_info) = rps_info else {
        return Err(-1);
    };
    let profile_header = &rps_info.profile_header;
    if profile_header.num_profiles < 0 {
        return Err(-1);
    }
    let num_profiles = profile_header.num_profiles.max(0) as usize;
    if profile_header.start_offsets.len() < num_profiles + 1 {
        return Err(-1);
    }

    let Some(alphabet_size) = rps_alphabet_size(profile_header.magic_number) else {
        return Err(-1);
    };
    let num_pssm_rows = profile_header.start_offsets[num_profiles];
    if num_pssm_rows < 0 {
        return Err(-1);
    }
    let num_pssm_rows = num_pssm_rows as usize;
    let rps_pssm = rps_rows_from_values(&profile_header.pssm_values, num_pssm_rows, alphabet_size)
        .ok_or(-1i16)?;

    let rps_freq = if let Some(freq_header) = &rps_info.freq_ratios_header {
        if freq_header.start_offsets.len() < num_profiles + 1 {
            return Err(-1);
        }
        for profile_index in 0..num_profiles {
            let start = freq_header.start_offsets[profile_index];
            let end = freq_header.start_offsets[profile_index + 1];
            if start < 0 || end < start {
                return Err(-1);
            }
            if start != profile_header.start_offsets[profile_index]
                || end != profile_header.start_offsets[profile_index + 1]
            {
                return Err(-1);
            }
        }
        Some(
            rps_rows_from_values(&freq_header.freq_values, num_pssm_rows, alphabet_size)
                .ok_or(-1i16)?,
        )
    } else {
        None
    };
    let rps_freq_ratios = if let Some(freq) = &rps_freq {
        Some(
            s_rps_fill_freq_ratios_in_psi_matrix(freq, num_pssm_rows + 1, alphabet_size)
                .ok_or(-1i16)?,
        )
    } else {
        None
    };

    let mut contexts = Vec::with_capacity(num_profiles);
    let mut max_length = 0u32;
    let mut min_length = u32::MAX;
    for profile_index in 0..num_profiles {
        let start = profile_header.start_offsets[profile_index];
        let end = profile_header.start_offsets[profile_index + 1];
        if start < 0 || end < start {
            return Err(-1);
        }
        let query_length = end - start;
        let ql_u32 = query_length as u32;
        max_length = max_length.max(ql_u32);
        min_length = min_length.min(ql_u32);
        contexts.push(crate::queryinfo::ContextInfo {
            query_offset: start,
            query_length,
            eff_searchsp: 0,
            length_adjustment: 0,
            query_index: profile_index as i32,
            frame: 0,
            is_valid: query_length > 0,
            segment_flags: crate::queryinfo::E_NO_SEGMENTS,
        });
    }
    if min_length == u32::MAX {
        min_length = 0;
    }

    Ok(RpsGapAlignData {
        concat_db_info: crate::queryinfo::QueryInfo {
            num_queries: num_profiles as i32,
            contexts,
            max_length,
            min_length,
        },
        rps_pssm,
        rps_freq,
        rps_freq_ratios,
        alphabet_size,
        position_based: true,
    })
}

/// blast-rs: Port-shaped Rust equivalent of NCBI `s_RPSComputeTraceback`
/// (`blast_traceback.c:1058`) for the represented Rust stream/traceback layer.
/// not a direct NCBI C port.
///
/// The caller supplies the concrete traceback operation because Rust does not
/// yet expose the C `SeqSrc`/`BlastGapAlignStruct` pointer graph. This function
/// preserves the C orchestration: validate inputs, prepare concatenated RPS
/// profile contexts, drain the HSP stream by OID, select per-profile PSSM and
/// frequency rows, apply the RPS Karlin K multiplier, run traceback, restore
/// RPS HSP orientation, insert nonempty lists, close traceback pipes, and sort.
pub fn s_rps_compute_traceback<T>(
    program_number: crate::program::ProgramType,
    hsp_stream: Option<&HspStream>,
    query_info: Option<&crate::queryinfo::QueryInfo>,
    rps_info: Option<&RpsTracebackInfo>,
    hit_options: &HitSavingOptions,
    results: Option<&mut HspResults>,
    traceback_hsp_list: T,
    interrupt_search: Option<&mut dyn FnMut(&crate::util::SBlastProgress) -> bool>,
    progress_info: Option<&mut crate::util::SBlastProgress>,
) -> i16
where
    T: FnMut(&mut HspList, &RpsGapAlignData, usize, Option<f64>) -> i16,
{
    s_rps_compute_traceback_with_composition(
        program_number,
        hsp_stream,
        query_info,
        rps_info,
        false,
        hit_options,
        results,
        traceback_hsp_list,
        interrupt_search,
        progress_info,
    )
}

/// CBS-aware variant of [`s_rps_compute_traceback`].
/// blast-rs: Composition-aware Rust wrapper around RPS traceback orchestration; not a direct NCBI C port.
///
/// NCBI only installs `RPS_K_MULT * karlin_k[oid]` for RPS-tblastn inside the
/// composition-based branch; ordinary RPS-BLAST installs it for both CBS and
/// non-CBS paths.
pub fn s_rps_compute_traceback_with_composition<T>(
    program_number: crate::program::ProgramType,
    hsp_stream: Option<&HspStream>,
    query_info: Option<&crate::queryinfo::QueryInfo>,
    rps_info: Option<&RpsTracebackInfo>,
    composition_based_stats: bool,
    hit_options: &HitSavingOptions,
    results: Option<&mut HspResults>,
    mut traceback_hsp_list: T,
    mut interrupt_search: Option<&mut dyn FnMut(&crate::util::SBlastProgress) -> bool>,
    mut progress_info: Option<&mut crate::util::SBlastProgress>,
) -> i16
where
    T: FnMut(&mut HspList, &RpsGapAlignData, usize, Option<f64>) -> i16,
{
    let Some(stream) = hsp_stream else {
        return -1;
    };
    let Some(query_info) = query_info else {
        return -1;
    };
    let Some(results) = results else {
        return -1;
    };
    let Ok(gap_data) = s_rps_gap_align_data_prepare(rps_info) else {
        return -1;
    };

    if let Some(progress) = progress_info.as_deref_mut() {
        progress.stage = crate::util::EBlastStage::TracebackSearch;
    }

    loop {
        let (read_status, maybe_hsp_list) = blast_hsp_stream_read(Some(stream));
        if read_status == K_BLAST_HSP_STREAM_EOF {
            break;
        }
        if read_status != K_BLAST_HSP_STREAM_SUCCESS {
            BlastHSPStreamTBackClose(Some(stream), Some(&mut *results));
            return read_status as i16;
        }

        if let (Some(interrupt), Some(progress)) =
            (interrupt_search.as_deref_mut(), progress_info.as_deref())
        {
            if interrupt(progress) {
                BlastHSPStreamTBackClose(Some(stream), Some(&mut *results));
                let _ = blast_hsp_results_sort_by_evalue(results);
                return crate::diagnostics::BLASTERR_INTERRUPTED;
            }
        }

        let Some(mut hsp_list) = maybe_hsp_list else {
            continue;
        };
        if hsp_list.oid < 0 {
            continue;
        }
        let profile_index = hsp_list.oid as usize;
        let Some(profile_context) = gap_data.concat_db_info.contexts.get(profile_index) else {
            continue;
        };
        if !profile_context.is_valid {
            continue;
        }
        let db_seq_start = profile_context.query_offset.max(0) as usize;
        if db_seq_start >= gap_data.rps_pssm.len() {
            continue;
        }
        if let Some(freqs) = &gap_data.rps_freq {
            if db_seq_start >= freqs.len() {
                continue;
            }
        }
        let karlin_k = rps_info
            .and_then(|info| info.karlin_k.get(profile_index).copied())
            .filter(|_| program_number != crate::program::RPS_TBLASTN || composition_based_stats)
            .map(|k| RPS_K_MULT * k);

        let status = traceback_hsp_list(&mut hsp_list, &gap_data, profile_index, karlin_k);
        if status != 0 {
            BlastHSPStreamTBackClose(Some(stream), Some(&mut *results));
            return status;
        }

        s_blast_hsp_list_rps_update(program_number, &mut hsp_list);
        if hsp_list.hsps.is_empty() {
            continue;
        }
        let query_index = query_index_for_hsp_list(&hsp_list, query_info);
        let insert_status = blast_hsp_results_insert_hsp_list(
            results,
            query_index,
            hsp_list,
            hit_options.hitlist_size,
        );
        if insert_status != 0 {
            BlastHSPStreamTBackClose(Some(stream), Some(&mut *results));
            return insert_status as i16;
        }
    }

    BlastHSPStreamTBackClose(Some(stream), Some(&mut *results));
    let _ = blast_hsp_results_sort_by_evalue(results);
    0
}

/// blast-rs: Port-shaped coordinator for NCBI `s_HSPListPostTracebackUpdate`
/// (`blast_traceback.c:199`).
/// not a direct NCBI C port.
///
/// The C function is mostly orchestration: restore RPS HSP orientation, either
/// link HSPs or recompute ordinary e-values, reap by e-value, rescale raw
/// scores, and fill bit scores. Rust keeps the same sequencing while using the
/// already-translated HSP/link/stat helpers. Score-block data needed by the
/// linker is reconstructed from the Rust `BlastScoreBlk` view.
pub fn s_hsp_list_post_traceback_update(
    program_number: crate::program::ProgramType,
    hsp_list: &mut BlastHSPList,
    query_info: &crate::queryinfo::QueryInfo,
    score_params: &crate::parameters::ScoringParameters,
    hit_params: &crate::parameters::HitSavingParameters,
    sbp: &crate::stat::BlastScoreBlk,
    subject_length: i32,
) -> i16 {
    let k_gapped = score_params.options.gapped_calculation;
    if crate::program::blast_program_is_rps_blast(program_number) {
        for hsp in hsp_list.hsp_array.iter_mut().flatten() {
            std::mem::swap(&mut hsp.query, &mut hsp.subject);
            if let Some(edit_script) = hsp.gap_info.as_mut() {
                for index in 0..edit_script.len() {
                    let (op, count) = edit_script.get(index).unwrap_or((GapAlignOpType::Sub, 0));
                    let new_op = match op {
                        GapAlignOpType::Ins => GapAlignOpType::Del,
                        GapAlignOpType::Del => GapAlignOpType::Ins,
                        other => other,
                    };
                    edit_script.set(index, new_op, count);
                }
            }
            if program_number == crate::program::RPS_TBLASTN {
                hsp.context =
                    crate::util::blast_frame_to_context(hsp.query.frame as i32, program_number);
            }
        }
        hsp_list.hsp_array.sort_by(|a, b| score_compare_blast_hsps(a.as_ref(), b.as_ref()));
    }

    if let Some(link_hsp_params) = hit_params.link_hsp_params.as_ref() {
        let link_score_block = crate::link_hsps::LinkScoreBlock {
            kbp: sbp.kbp.clone(),
            kbp_gap: sbp.kbp_gap.clone(),
            gbp: None,
            link_gbp_db_length: None,
            recompute_evalues_before_uneven_linking: false,
        };
        let _ = crate::link_hsps::blast_link_hsp_list(
            program_number,
            hsp_list,
            query_info,
            subject_length,
            &link_score_block,
            link_hsp_params,
            k_gapped,
        );
    } else {
        let scale_factor = if crate::program::blast_program_is_rps_blast(program_number) {
            score_params.scale_factor
        } else {
            1.0
        };
        let source_query_index = hsp_list.query_index;
        let source_hsp_max = hsp_list.hsp_max;
        let mut legacy_list = std::mem::replace(
            hsp_list,
            BlastHSPList {
                oid: -1,
                query_index: source_query_index,
                hsp_array: Vec::new(),
                hspcnt: 0,
                allocated: 0,
                hsp_max: source_hsp_max,
                do_not_reallocate: false,
                best_evalue: i32::MAX as f64,
            },
        )
        .into_legacy_hsp_list();
        let _ = crate::blast_kappa::blast_hsp_list_get_evalues(
            program_number,
            query_info,
            subject_length,
            &mut legacy_list,
            k_gapped,
            false,
            sbp,
            0.0,
            scale_factor,
        );
        *hsp_list = BlastHSPList::from_legacy_hsp_list(legacy_list);
        hsp_list.query_index = source_query_index;
        hsp_list.hsp_max = source_hsp_max;
    }

    hsp_list
        .hsp_array
        .retain(|hsp| hsp.as_ref().is_some_and(|hsp| hsp.evalue <= hit_params.options.expect_value));
    hsp_list.hspcnt = hsp_list.hsp_array.len() as i32;
    hsp_list.allocated = hsp_list.hspcnt;
    for hsp in hsp_list.hsp_array.iter_mut().flatten() {
        hsp.score = ((hsp.score as f64 + 0.5 * score_params.scale_factor)
            / score_params.scale_factor) as i32;
    }
    hsp_list.hsp_array.sort_by(|a, b| score_compare_blast_hsps(a.as_ref(), b.as_ref()));

    let kbp_array = if k_gapped && !sbp.kbp_gap.is_empty() {
        &sbp.kbp_gap
    } else {
        &sbp.kbp
    };
    for hsp in hsp_list.hsp_array.iter_mut().flatten() {
        if let Some(ctx_kbp) = kbp_array.get(hsp.context.max(0) as usize) {
            hsp.bit_score =
                (hsp.score as f64 * ctx_kbp.lambda - ctx_kbp.log_k) / crate::math::NCBIMATH_LN2;
        }
    }
    hsp_list.best_evalue = hsp_list
        .hsp_array
        .iter()
        .filter_map(|hsp| hsp.as_ref().map(|hsp| hsp.evalue))
        .fold(i32::MAX as f64, f64::min);
    0
}

/// NCBI: Blast_HSPCalcLengthAndGaps (blast_hits.c:1055).
pub fn blast_hsp_calc_length_and_gaps(hsp: &Hsp) -> (i32, i32, i32) {
    let mut length = hsp.query_end - hsp.query_offset;
    let subject_length = hsp.subject_end - hsp.subject_offset;
    let mut gaps = 0;
    let mut gap_opens = 0;

    if let Some(edit_script) = &hsp.edit_script {
        for (op, count) in edit_script.iter() {
            // NCBI `Blast_HSPCalcLengthAndGaps` (`blast_hits.c:1065`) only
            // counts `eGapAlignDel` and `eGapAlignIns`; frame-shift variants
            // `Del1/Del2/Ins1/Ins2` are NOT counted as gaps (they're
            // out-of-frame nucleotide shifts, handled separately).
            match op {
                crate::gapinfo::GapAlignOpType::Del => {
                    length += count;
                    gaps += count;
                    gap_opens += 1;
                }
                crate::gapinfo::GapAlignOpType::Ins => {
                    gaps += count;
                    gap_opens += 1;
                }
                _ => {}
            }
        }
    } else if subject_length > length {
        length = subject_length;
    }

    (length, gaps, gap_opens)
}

/// NCBI: Blast_HSPListAdjustOffsets (blast_hits.c:3037).
pub fn blast_hsp_list_adjust_offsets(hsp_list: &mut HspList, offset: i32) {
    if offset == 0 {
        return;
    }
    for hsp in &mut hsp_list.hsps {
        hsp.subject_offset += offset;
        hsp.subject_end += offset;
        hsp.subject_gapped_start += offset;
    }
}

/// NCBI: s_AdjustSubjectForTranslatedSraSearch (blast_engine.c:1207).
pub fn s_adjust_subject_for_translated_sra_search(hsp_list: &mut HspList, offset: u8, length: i32) {
    for hsp in &mut hsp_list.hsps {
        if hsp.subject_frame > 0 {
            let mut fix_co = false;
            match offset {
                1 => match hsp.subject_frame {
                    1 => {
                        hsp.subject_frame = 3;
                        fix_co = true;
                    }
                    2 => hsp.subject_frame = 1,
                    3 => hsp.subject_frame = 2,
                    _ => {}
                },
                2 => match hsp.subject_frame {
                    1 => {
                        hsp.subject_frame = 2;
                        fix_co = true;
                    }
                    2 => {
                        hsp.subject_frame = 3;
                        fix_co = true;
                    }
                    3 => hsp.subject_frame = 1,
                    _ => {}
                },
                3 => fix_co = true,
                _ => {}
            }

            if fix_co {
                if hsp.subject_offset == 0 {
                    hsp.query_offset += 1;
                }
                if hsp.subject_offset > 0 {
                    hsp.subject_offset -= 1;
                }
                hsp.subject_end -= 1;
                hsp.subject_gapped_start -= 1;
            }
        } else {
            let max_end = length - offset as i32;
            let subj_end =
                (hsp.subject_end + 1) * crate::util::CODON_LENGTH as i32 - hsp.subject_frame - 2;
            if subj_end > max_end {
                hsp.subject_end -= 1;
                hsp.query_end -= 1;
            }
        }
    }
}

/// NCBI: s_BlastHSPListsCombineByScore (blast_hits.c:2749).
/// naming: Rust spells HSP lists as separate snake_case tokens.
/// Port for Rust-owned HSP vectors.
pub fn s_blast_hsp_lists_combine_by_score(
    hsp_list: &mut HspList,
    combined_hsp_list: &mut HspList,
    new_hspcnt: i32,
) {
    let new_hspcnt = new_hspcnt.max(0) as usize;
    let total_hsps = hsp_list.hsps.len() + combined_hsp_list.hsps.len();

    if new_hspcnt >= total_hsps {
        combined_hsp_list.hsps.append(&mut hsp_list.hsps);
        combined_hsp_list.hsps.sort_by(score_compare_hsps);
        combined_hsp_list.best_evalue = s_blast_get_best_evalue(combined_hsp_list);
        return;
    }

    combined_hsp_list.hsps.sort_by(score_compare_hsps);
    hsp_list.hsps.sort_by(score_compare_hsps);

    let old_hsps = std::mem::take(&mut combined_hsp_list.hsps);
    let new_hsps = std::mem::take(&mut hsp_list.hsps);
    let mut old_iter = old_hsps.into_iter().peekable();
    let mut new_iter = new_hsps.into_iter().peekable();
    let mut merged = Vec::with_capacity(new_hspcnt);

    while merged.len() < new_hspcnt {
        let take_old = match (old_iter.peek(), new_iter.peek()) {
            (Some(old), Some(new)) => score_compare_hsps(old, new).is_le(),
            (Some(_), None) => true,
            (None, Some(_)) => false,
            (None, None) => break,
        };

        if take_old {
            if let Some(hsp) = old_iter.next() {
                merged.push(hsp);
            }
        } else if let Some(hsp) = new_iter.next() {
            merged.push(hsp);
        }
    }

    combined_hsp_list.hsps = merged;
    combined_hsp_list.best_evalue = s_blast_get_best_evalue(combined_hsp_list);
}

/// NCBI: Blast_HSPListAppend (blast_hits.c:2809).
pub fn blast_hsp_list_append(
    old_hsp_list: &mut Option<HspList>,
    combined_hsp_list: &mut Option<HspList>,
    hsp_num_max: i32,
) -> i32 {
    let Some(mut old) = old_hsp_list.take() else {
        return 0;
    };
    if old.hsps.is_empty() {
        *old_hsp_list = Some(old);
        return 0;
    }

    let Some(combined) = combined_hsp_list.as_mut() else {
        *combined_hsp_list = Some(old);
        return 0;
    };

    let new_hspcnt = (combined.hsps.len() + old.hsps.len()).min(hsp_num_max.max(0) as usize);
    s_blast_hsp_lists_combine_by_score(&mut old, combined, new_hspcnt as i32);
    0
}

/// NCBI: ScoreCompareHSPs (blast_hits.c:1330).
/// Total ordering
/// used by `Blast_HSPListSortByScore`: primary key is score descending,
/// with tie-breakers on `(subject_offset asc, subject_end desc,
/// query_offset asc, query_end desc)`.
/// naming: Rust keeps HSPs as a readable snake_case plural token.
pub fn score_compare_hsps(a: &Hsp, b: &Hsp) -> std::cmp::Ordering {
    b.score
        .cmp(&a.score)
        .then_with(|| a.subject_offset.cmp(&b.subject_offset))
        .then_with(|| b.subject_end.cmp(&a.subject_end))
        .then_with(|| a.query_offset.cmp(&b.query_offset))
        .then_with(|| b.query_end.cmp(&a.query_end))
}

pub fn score_compare_blast_hsps(
    a: Option<&BlastHSP>,
    b: Option<&BlastHSP>,
) -> std::cmp::Ordering {
    match (a, b) {
        (Some(a), Some(b)) => b
            .score
            .cmp(&a.score)
            .then_with(|| a.subject.offset.cmp(&b.subject.offset))
            .then_with(|| b.subject.end.cmp(&a.subject.end))
            .then_with(|| a.query.offset.cmp(&b.query.offset))
            .then_with(|| b.query.end.cmp(&a.query.end)),
        (Some(_), None) => std::cmp::Ordering::Less,
        (None, Some(_)) => std::cmp::Ordering::Greater,
        (None, None) => std::cmp::Ordering::Equal,
    }
}

/// NCBI: s_QueryOffsetCompareHSPs (blast_hits.c:2268).
/// naming: Rust keeps HSPs as a readable snake_case plural token.
pub fn s_query_offset_compare_hsps(a: Option<&Hsp>, b: Option<&Hsp>) -> std::cmp::Ordering {
    let (Some(a), Some(b)) = (a, b) else {
        return match (a.is_some(), b.is_some()) {
            (false, false) => std::cmp::Ordering::Equal,
            (false, true) => std::cmp::Ordering::Greater,
            (true, false) => std::cmp::Ordering::Less,
            (true, true) => std::cmp::Ordering::Equal,
        };
    };

    a.context
        .cmp(&b.context)
        .then_with(|| a.query_offset.cmp(&b.query_offset))
        .then_with(|| a.subject_offset.cmp(&b.subject_offset))
        .then_with(|| b.score.cmp(&a.score))
        .then_with(|| b.query_end.cmp(&a.query_end))
        .then_with(|| b.subject_end.cmp(&a.subject_end))
}

/// NCBI: s_QueryEndCompareHSPs (blast_hits.c:2332).
/// naming: Rust keeps HSPs as a readable snake_case plural token.
pub fn s_query_end_compare_hsps(a: Option<&Hsp>, b: Option<&Hsp>) -> std::cmp::Ordering {
    let (Some(a), Some(b)) = (a, b) else {
        return match (a.is_some(), b.is_some()) {
            (false, false) => std::cmp::Ordering::Equal,
            (false, true) => std::cmp::Ordering::Greater,
            (true, false) => std::cmp::Ordering::Less,
            (true, true) => std::cmp::Ordering::Equal,
        };
    };

    a.context
        .cmp(&b.context)
        .then_with(|| a.query_end.cmp(&b.query_end))
        .then_with(|| a.subject_end.cmp(&b.subject_end))
        .then_with(|| b.score.cmp(&a.score))
        .then_with(|| b.query_offset.cmp(&a.query_offset))
        .then_with(|| b.subject_offset.cmp(&a.subject_offset))
}

/// NCBI: s_CutOffGapEditScript (blast_hits.c:2392).
pub fn s_cut_off_gap_edit_script(hsp: &mut Hsp, q_cut: i32, s_cut: i32, cut_begin: bool) {
    let Some(edit_script) = hsp.edit_script.as_mut() else {
        return;
    };

    let q_cut = q_cut - hsp.query_offset;
    let s_cut = s_cut - hsp.subject_offset;
    let mut qid = 0;
    let mut sid = 0;
    let mut found = false;
    let mut found_index = 0usize;
    let mut found_opid = 0i32;

    'outer: for (index, (op, num)) in edit_script.iter().enumerate() {
        let mut opid = 0;
        while opid < num {
            match op {
                GapAlignOpType::Sub => {
                    qid += 1;
                    sid += 1;
                    opid += 1;
                }
                GapAlignOpType::Del => {
                    sid += num;
                    opid += num;
                }
                GapAlignOpType::Ins => {
                    qid += num;
                    opid += num;
                }
                _ => {
                    opid += num;
                }
            }

            if qid >= q_cut && sid >= s_cut {
                found = true;
                found_index = index;
                found_opid = opid;
            }
            if found {
                break 'outer;
            }
        }
    }

    if !found {
        return;
    }

    let (_, found_num) = edit_script
        .get(found_index)
        .unwrap_or((GapAlignOpType::Sub, 0));
    if cut_begin {
        let mut new_ops = Vec::new();
        if found_opid < found_num {
            let (found_op, _) = edit_script
                .get(found_index)
                .unwrap_or((GapAlignOpType::Sub, 0));
            debug_assert_eq!(found_op, GapAlignOpType::Sub);
            new_ops.push((found_op, found_num - found_opid));
        }
        new_ops.extend(edit_script.iter().skip(found_index + 1));
        edit_script.replace_ops(new_ops);
        hsp.query_offset += qid;
        hsp.subject_offset += sid;
    } else {
        if found_opid < found_num {
            let (found_op, _) = edit_script
                .get(found_index)
                .unwrap_or((GapAlignOpType::Sub, 0));
            debug_assert_eq!(found_op, GapAlignOpType::Sub);
            edit_script.set_num(found_index, found_opid);
        }
        edit_script.truncate(found_index + 1);
        hsp.query_end = hsp.query_offset + qid;
        hsp.subject_end = hsp.subject_offset + sid;
    }
}

/// Results for one query: a collection of HspLists from different subjects.
#[derive(Debug, Clone)]
pub struct HitList {
    pub hsp_lists: Vec<HspList>,
    pub hsplist_max: i32,
    pub low_score: i32,
    pub worst_evalue: f64,
}

/// NCBI: `BlastHitList` (`blast_hits.h:161`).
#[derive(Debug, Clone)]
pub struct BlastHitList {
    pub hsplist_count: i32,
    pub hsplist_max: i32,
    pub worst_evalue: f64,
    pub low_score: i32,
    pub heapified: bool,
    pub hsplist_array: Vec<Option<BlastHSPList>>,
    pub hsplist_current: i32,
    pub num_hits: i32,
}

impl BlastHitList {
    pub fn from_legacy_hit_list(hitlist: HitList) -> Self {
        let hsplist_array: Vec<Option<BlastHSPList>> = hitlist
            .hsp_lists
            .into_iter()
            .map(BlastHSPList::from_legacy_hsp_list)
            .map(Some)
            .collect();
        let hsplist_count = hsplist_array.iter().filter(|list| list.is_some()).count() as i32;
        Self {
            hsplist_count,
            hsplist_max: hitlist.hsplist_max,
            worst_evalue: hitlist.worst_evalue,
            low_score: hitlist.low_score,
            heapified: false,
            hsplist_current: hsplist_array.len() as i32,
            hsplist_array,
            num_hits: 0,
        }
    }

    pub fn into_legacy_hit_list(self) -> HitList {
        let hsplist_count = self.hsplist_count.max(0) as usize;
        HitList {
            hsp_lists: self
                .hsplist_array
                .into_iter()
                .take(hsplist_count)
                .flatten()
                .map(BlastHSPList::into_legacy_hsp_list)
                .collect(),
            hsplist_max: self.hsplist_max,
            low_score: self.low_score,
            worst_evalue: self.worst_evalue,
        }
    }
}

impl HitList {
    /// blast-rs: Rust convenience constructor for owned hit lists; not a direct NCBI C port.
    pub fn new() -> Self {
        blast_hit_list_new(0)
    }

    /// NCBI: Blast_HitListUpdate (blast_hits.c:3243).
    /// Port for the
    /// Rust `HitList` shape.
    pub fn blast_hit_list_update(&mut self, mut hsp_list: HspList) -> i32 {
        hsp_list.best_evalue = s_blast_get_best_evalue(&hsp_list);
        if hsp_list.hsps.is_empty() {
            return 0;
        }
        let max = self.hsplist_max.max(0) as usize;
        if max == 0 || self.hsp_lists.len() < max {
            self.low_score = self.low_score.min(hsp_list.hsps[0].score);
            self.worst_evalue = self.worst_evalue.max(hsp_list.best_evalue);
            self.hsp_lists.push(hsp_list);
            if max != 0 && self.hsp_lists.len() == max {
                for list in &mut self.hsp_lists {
                    let _ = blast_hsp_list_sort_by_evalue(Some(list));
                    list.best_evalue = s_blast_get_best_evalue(list);
                }
                s_create_heap(&mut self.hsp_lists, evalue_compare_hsp_lists);
                if let Some(worst) = self.hsp_lists.first() {
                    self.worst_evalue = worst.best_evalue;
                    self.low_score = worst.hsps.first().map(|hsp| hsp.score).unwrap_or(i32::MAX);
                }
            }
            return 0;
        }

        let _ = blast_hsp_list_sort_by_evalue(Some(&mut hsp_list));
        hsp_list.best_evalue = s_blast_get_best_evalue(&hsp_list);
        if !evalue_compare_hsp_lists(&self.hsp_lists[0], &hsp_list).is_lt() {
            s_blast_hit_list_insert_hsp_list_in_heap(self, hsp_list);
        }
        0
    }

    /// blast-rs: Rust convenience wrapper around hit-list e-value sorting; not a direct NCBI C port.
    pub fn sort_by_evalue(&mut self) {
        let _ = blast_hit_list_sort_by_evalue(self);
    }
}

/// NCBI: Blast_HitListNew (blast_hits.c:3125).
pub fn blast_hit_list_new(hitlist_size: i32) -> HitList {
    HitList {
        hsp_lists: Vec::new(),
        hsplist_max: hitlist_size,
        low_score: i32::MAX,
        worst_evalue: 0.0,
    }
}

/// NCBI: Blast_HitListMerge (blast_hits.c:2119).
/// naming: Rust uses owned `Option<HitList>` parameters instead of C pointers.
pub fn blast_hit_list_merge(
    old_hit_list: &mut Option<HitList>,
    combined_hit_list: &mut Option<HitList>,
    contexts_per_query: i32,
    split_offsets: Option<&[i32]>,
    chunk_overlap_size: i32,
    allow_gap: bool,
) -> i32 {
    let Some(mut incoming) = old_hit_list.take() else {
        return 0;
    };
    if combined_hit_list.is_none() {
        *combined_hit_list = Some(incoming);
        return 0;
    }
    let Some(mut combined) = combined_hit_list.take() else {
        *combined_hit_list = Some(incoming);
        return 0;
    };

    let query_is_split = split_offsets.is_some_and(|offsets| {
        offsets
            .iter()
            .take(contexts_per_query.max(0) as usize)
            .any(|&offset| offset > 0)
    });
    incoming.hsp_lists.sort_by_key(|list| list.oid);
    combined.hsp_lists.sort_by_key(|list| list.oid);

    let mut merged_hitlist = blast_hit_list_new(incoming.hsplist_max);
    let mut incoming_iter = incoming.hsp_lists.into_iter().peekable();
    let mut combined_iter = combined.hsp_lists.into_iter().peekable();
    loop {
        match (
            incoming_iter.peek().map(|list| list.oid),
            combined_iter.peek().map(|list| list.oid),
        ) {
            (Some(incoming_oid), Some(combined_oid)) if incoming_oid == combined_oid => {
                let mut incoming_list = Some(incoming_iter.next().unwrap());
                let mut merged_list = Some(combined_iter.next().unwrap());
                let hsp_num_max = merged_list.as_ref().map(|list| list.hsp_max).unwrap_or(0);
                if query_is_split {
                    let _ = blast_hsp_lists_merge(
                        &mut incoming_list,
                        &mut merged_list,
                        hsp_num_max,
                        split_offsets,
                        contexts_per_query,
                        chunk_overlap_size,
                        allow_gap,
                        false,
                    );
                } else {
                    let _ =
                        blast_hsp_list_append(&mut incoming_list, &mut merged_list, hsp_num_max);
                }
                if let Some(merged) = merged_list {
                    let _ = merged_hitlist.blast_hit_list_update(merged);
                }
            }
            (Some(incoming_oid), Some(combined_oid)) if incoming_oid < combined_oid => {
                let _ = merged_hitlist.blast_hit_list_update(incoming_iter.next().unwrap());
            }
            (Some(_), Some(_)) => {
                let _ = merged_hitlist.blast_hit_list_update(combined_iter.next().unwrap());
            }
            (Some(_), None) => {
                let _ = merged_hitlist.blast_hit_list_update(incoming_iter.next().unwrap());
            }
            (None, Some(_)) => {
                let _ = merged_hitlist.blast_hit_list_update(combined_iter.next().unwrap());
            }
            (None, None) => break,
        }
    }

    // NCBI `Blast_HitListMerge` (`blast_hits.c:2216`) returns the merged
    // hitlist without a final evalue sort. Callers that need ordered output
    // (e.g. `Blast_HSPResultsSortByEvalue` before display) sort downstream.
    *combined_hit_list = Some(merged_hitlist);
    0
}

/// blast-rs: Convenience wrapper around `blast_hit_list_merge`; not a direct NCBI C port.
pub fn blast_hit_list_merge_simple(
    old_hit_list: &mut Option<HitList>,
    combined_hit_list: &mut Option<HitList>,
) -> i32 {
    blast_hit_list_merge(old_hit_list, combined_hit_list, 0, None, 0, true)
}

/// NCBI: Blast_HitListHSPListsFree (blast_hits.c:3155).
/// 1-1 translation for
/// the Rust vector-backed representation.
pub fn blast_hit_list_hsp_lists_free(hitlist: Option<&mut HitList>) -> i32 {
    let Some(hitlist) = hitlist else { return 0 };
    hitlist.hsp_lists.clear();
    0
}

/// NCBI: Blast_HitListFree (blast_hits.c:3139).
pub fn blast_hit_list_free(hitlist: &mut Option<HitList>) -> Option<HitList> {
    if let Some(hitlist) = hitlist {
        let _ = blast_hit_list_hsp_lists_free(Some(hitlist));
    }
    *hitlist = None;
    None
}

/// NCBI: s_BlastHitListPurge (blast_hits.c:3170).
/// 1-1 translation for the
/// Rust vector-backed hit list. C frees the first empty HSP list and every
/// list after it, even if a later list is non-empty.
pub fn s_blast_hit_list_purge(hitlist: Option<&mut HitList>) {
    let Some(hitlist) = hitlist else { return };
    if let Some(first_empty) = hitlist
        .hsp_lists
        .iter()
        .position(|hsp_list| hsp_list.hsps.is_empty())
    {
        hitlist.hsp_lists.truncate(first_empty);
    }
}

/// NCBI: Blast_HitListSortByEvalue (blast_hits.c:3330).
pub fn blast_hit_list_sort_by_evalue(hitlist: &mut HitList) -> i32 {
    if hitlist.hsp_lists.len() > 1 {
        hitlist.hsp_lists.sort_by(evalue_compare_hsp_lists);
    }
    s_blast_hit_list_purge(Some(hitlist));
    0
}

/// NCBI: s_BlastHitListInsertHSPListInHeap (blast_hits.c:3196).
/// Port for the Rust vector-backed hit list.
pub fn s_blast_hit_list_insert_hsp_list_in_heap(hit_list: &mut HitList, hsp_list: HspList) {
    if hit_list.hsp_lists.is_empty() {
        hit_list.hsp_lists.push(hsp_list);
    } else {
        hit_list.hsp_lists[0] = hsp_list;
        if hit_list.hsp_lists.len() >= 2 {
            let lim = hit_list.hsp_lists.len() / 2 - 1;
            let last = hit_list.hsp_lists.len() - 1;
            s_heapify(
                &mut hit_list.hsp_lists,
                0,
                lim,
                last,
                evalue_compare_hsp_lists,
            );
        }
    }

    if let Some(worst) = hit_list.hsp_lists.first() {
        hit_list.worst_evalue = worst.best_evalue;
        hit_list.low_score = worst.hsps.first().map(|hsp| hsp.score).unwrap_or(i32::MAX);
    }
}

/// NCBI: Blast_HSPListPurgeNullHSPs (blast_hits.c:2225).
/// naming: Rust keeps HSPs as a readable snake_case plural token.
/// `Vec<Hsp>` cannot contain null HSP pointers, so the Rust representation is
/// already compact and this is a no-op. NCBI does NOT recompute
/// `hsp_list->best_evalue` here either — callers are responsible, so this
/// intentionally leaves cached list state untouched.
pub fn blast_hsp_list_purge_null_hsps(_hsp_list: Option<&mut HspList>) -> i32 {
    0
}

/// NCBI: Blast_HitListPurgeNullHSPLists (blast_hits.c:3302).
/// `Vec<HspList>` cannot contain null HSP-list pointers, so the Rust
/// representation is already compact and this is a no-op. NCBI also does
/// no other work here (the C version just walks the array dropping nulls),
/// including no recomputation of `low_score` or `worst_evalue`.
pub fn blast_hit_list_purge_null_hsp_lists(_hitlist: Option<&mut HitList>) -> i32 {
    0
}

/// NCBI: s_BlastGetBestEvalue (blast_hits.c:1742).
fn s_blast_get_best_evalue(hsp_list: &HspList) -> f64 {
    // NCBI `s_BlastGetBestEvalue` (`blast_hits.c:1742`) seeds with `(double)INT4_MAX`.
    hsp_list
        .hsps
        .iter()
        .map(|hsp| hsp.evalue)
        .fold(i32::MAX as f64, f64::min)
}

/// NCBI: s_EvalueComp (blast_hits.c:1390).
/// naming: Rust drops the private `s_` prefix for this comparator.
/// Two e-values that
/// are both below `1.0e-180` compare equal; otherwise ordinary ordering.
pub fn evalue_comp(evalue1: f64, evalue2: f64) -> std::cmp::Ordering {
    const EPSILON: f64 = 1.0e-180;
    if evalue1 < EPSILON && evalue2 < EPSILON {
        return std::cmp::Ordering::Equal;
    }
    evalue1
        .partial_cmp(&evalue2)
        .unwrap_or(std::cmp::Ordering::Equal)
}

/// NCBI: s_EvalueCompareHSPLists (blast_hits.c:3078).
/// naming: Rust drops the private `s_` prefix and keeps HSP lists readable.
/// Sort
/// order for HSP lists: primary key is best_evalue (via `evalue_comp`),
/// then top-HSP score descending, then `oid` descending. Empty lists
/// sort to the end.
pub fn evalue_compare_hsp_lists(a: &HspList, b: &HspList) -> std::cmp::Ordering {
    use std::cmp::Ordering::*;
    match (a.hsps.is_empty(), b.hsps.is_empty()) {
        (true, true) => return Equal,
        (true, false) => return Greater,
        (false, true) => return Less,
        (false, false) => {}
    }
    let by_evalue = evalue_comp(a.best_evalue, b.best_evalue);
    if by_evalue != Equal {
        return by_evalue;
    }
    // NCBI `blast_hits.c:3099-3102`: top HSP score, descending.
    let top_a = a.hsps[0].score;
    let top_b = b.hsps[0].score;
    let by_score = top_b.cmp(&top_a);
    if by_score != Equal {
        return by_score;
    }
    // NCBI `blast_hits.c:3106`: `BLAST_CMP(h2->oid, h1->oid)` — oid descending.
    b.oid.cmp(&a.oid)
}

/// Complete search results for all queries.
#[derive(Debug, Clone)]
pub struct HspResults {
    pub hitlists: Vec<Option<HitList>>,
}

/// NCBI: `BlastHSPResults` (`blast_hits.h:173`).
#[derive(Debug, Clone)]
pub struct BlastHSPResults {
    pub num_queries: i32,
    pub hitlist_array: Vec<Option<BlastHitList>>,
}

impl BlastHSPResults {
    pub fn from_legacy_hsp_results(results: HspResults) -> Self {
        let hitlist_array: Vec<Option<BlastHitList>> = results
            .hitlists
            .into_iter()
            .map(|hitlist| hitlist.map(BlastHitList::from_legacy_hit_list))
            .collect();
        Self {
            num_queries: hitlist_array.len() as i32,
            hitlist_array,
        }
    }

    pub fn into_legacy_hsp_results(self) -> HspResults {
        let num_queries = self.num_queries.max(0) as usize;
        HspResults {
            hitlists: self
                .hitlist_array
                .into_iter()
                .take(num_queries)
                .map(|hitlist| hitlist.map(BlastHitList::into_legacy_hit_list))
                .collect(),
        }
    }
}

/// Rust ownership equivalent of traceback MT `SThreadLocalData`.
#[derive(Debug, Clone, Default)]
pub struct SThreadLocalData {
    pub gap_align: Option<crate::blast_kappa::BlastGapAlignWorkspace>,
    pub score_params: Option<crate::parameters::ScoringParameters>,
    pub ext_params: Option<crate::parameters::ExtensionParameters>,
    pub hit_params: Option<crate::parameters::HitSavingParameters>,
    pub eff_len_params: Option<crate::parameters::EffectiveLengthsParameters>,
    pub seqsrc: Option<*const dyn crate::seqsrc::BlastSeqSource>,
    pub query_info: Option<crate::queryinfo::QueryInfo>,
    pub results: Option<HspResults>,
}

/// Rust equivalent of traceback MT `SThreadLocalDataArray`.
#[derive(Debug, Clone, Default)]
pub struct SThreadLocalDataArray {
    pub tld: Vec<Option<SThreadLocalData>>,
    pub num_elems: u32,
}

#[derive(Debug, Clone, Default)]
pub struct BlastHspStreamResultBatch {
    pub hsplist_array: Vec<Option<HspList>>,
    pub num_hsplists: i32,
}

#[derive(Debug, Clone, Default)]
pub struct BlastHspStreamResultsBatchArray {
    pub array_of_batches: Vec<Option<BlastHspStreamResultBatch>>,
    pub num_batches: u32,
    pub num_allocated: u32,
}

impl HspResults {
    /// blast-rs: Rust convenience constructor for owned HSP results; not a direct NCBI C port.
    pub fn new(num_queries: i32) -> Self {
        blast_hsp_results_new(num_queries)
    }

    /// blast-rs: Rust convenience wrapper around result e-value sorting; not a direct NCBI C port.
    pub fn sort_by_evalue(&mut self) {
        let _ = blast_hsp_results_sort_by_evalue(self);
    }
}

/// NCBI: Blast_HSPResultsNew (blast_hits.c:3346).
pub fn blast_hsp_results_new(num_queries: i32) -> HspResults {
    HspResults {
        hitlists: (0..num_queries.max(0)).map(|_| None).collect(),
    }
}

/// NCBI: Blast_HSPResultsFree (blast_hits.c:3366).
pub fn blast_hsp_results_free(results: &mut Option<HspResults>) -> Option<HspResults> {
    if let Some(results) = results {
        for hitlist in &mut results.hitlists {
            let _ = blast_hit_list_free(hitlist);
        }
    }
    *results = None;
    None
}

/// NCBI: SThreadLocalDataNew (blast_traceback_mt_priv.c:37).
/// naming: Historical Rust spelling keeps `sthread` as one token.
pub fn sthread_local_data_new() -> SThreadLocalData {
    SThreadLocalData::default()
}

/// blast-rs: Rust ownership equivalent of the C thread-local-data free routine; not a direct NCBI C port.
pub fn sthread_local_data_free(tld: &mut Option<SThreadLocalData>) -> Option<SThreadLocalData> {
    if let Some(tld) = tld {
        tld.gap_align = None;
        tld.score_params = None;
        tld.ext_params = None;
        tld.hit_params = None;
        tld.eff_len_params = None;
        tld.seqsrc = None;
        let _ = blast_hsp_results_free(&mut tld.results);
        tld.query_info = None;
    }
    *tld = None;
    None
}

/// NCBI: SThreadLocalDataArrayNew (blast_traceback_mt_priv.c).
/// naming: Historical Rust spelling keeps `sthread` as one token.
pub fn sthread_local_data_array_new(num_threads: u32) -> SThreadLocalDataArray {
    SThreadLocalDataArray {
        tld: (0..num_threads)
            .map(|_| Some(sthread_local_data_new()))
            .collect(),
        num_elems: num_threads,
    }
}

/// blast-rs: Rust ownership equivalent of the C thread-local-data-array free routine; not a direct NCBI C port.
pub fn sthread_local_data_array_free(
    array: &mut Option<SThreadLocalDataArray>,
) -> Option<SThreadLocalDataArray> {
    if let Some(array) = array {
        for tld in array.tld.iter_mut().take(array.num_elems as usize) {
            let _ = sthread_local_data_free(tld);
        }
        array.tld.clear();
        array.num_elems = 0;
    }
    *array = None;
    None
}

/// NCBI: SThreadLocalDataArrayTrim (blast_traceback_mt_priv.c).
/// naming: Historical Rust spelling keeps `sthread` as one token.
pub fn sthread_local_data_array_trim(
    array: Option<&mut SThreadLocalDataArray>,
    actual_num_threads: u32,
) {
    let Some(array) = array else { return };
    if actual_num_threads >= array.num_elems {
        return;
    }
    for tld in array
        .tld
        .iter_mut()
        .take(array.num_elems as usize)
        .skip(actual_num_threads as usize)
    {
        let _ = sthread_local_data_free(tld);
    }
    array.tld.truncate(actual_num_threads as usize);
    array.num_elems = actual_num_threads;
}

/// blast-rs: Port-shaped setup for NCBI `SThreadLocalDataArraySetup`; not a direct NCBI C port.
/// naming: Historical Rust spelling keeps `sthread` as one token.
///
/// C also calls `BLAST_GapAlignSetUp` and copies a `BlastSeqSrc`; those raw
/// handles are represented elsewhere in Rust. This function mirrors the
/// per-thread query-info duplication and result allocation used by later
/// traceback consolidation.
pub fn sthread_local_data_array_setup(
    array: Option<&mut SThreadLocalDataArray>,
    query_info: &crate::queryinfo::QueryInfo,
    hitlist_size: i32,
) -> i16 {
    const BLASTERR_INVALIDPARAM: i16 = -1;
    let Some(array) = array else {
        return BLASTERR_INVALIDPARAM;
    };

    for tld in array.tld.iter_mut().take(array.num_elems as usize) {
        let Some(tld) = tld else {
            return BLASTERR_INVALIDPARAM;
        };
        let mut hit_params = crate::parameters::HitSavingParameters::default();
        hit_params.options.hitlist_size = hitlist_size;
        tld.query_info = Some(query_info.clone());
        tld.results = Some(blast_hsp_results_new(query_info.num_queries));
        tld.hit_params = Some(hit_params);
    }
    0
}

/// blast-rs: Thread-local consolidation sizing helper; not a direct NCBI C port.
fn s_count_hsp_lists_per_query(data: Option<&SThreadLocalDataArray>) -> Option<Vec<usize>> {
    let data = data?;
    let first_results = data.tld.first()?.as_ref()?.results.as_ref()?;
    let mut counts = vec![0; first_results.hitlists.len()];

    for tld in data.tld.iter().take(data.num_elems as usize).flatten() {
        let Some(results) = tld.results.as_ref() else {
            continue;
        };
        for (query_idx, hitlist) in results.hitlists.iter().enumerate().take(counts.len()) {
            if let Some(hitlist) = hitlist {
                counts[query_idx] += hitlist.hsp_lists.len();
            }
        }
    }
    Some(counts)
}

/// NCBI: SThreadLocalDataArrayConsolidateResults (blast_traceback_mt_priv.c).
/// naming: Historical Rust spelling keeps `sthread` as one token.
///
/// HSP lists are moved from each thread-local result into the returned
/// per-query hit lists. Empty lists are skipped, matching
/// `Blast_HSPList_IsEmpty`.
pub fn sthread_local_data_array_consolidate_results(
    array: Option<&mut SThreadLocalDataArray>,
) -> Option<HspResults> {
    let array = array?;
    let counts = s_count_hsp_lists_per_query(Some(array))?;
    let num_queries = counts.len() as i32;
    let hitlist_size = array
        .tld
        .first()
        .and_then(|tld| tld.as_ref())
        .and_then(|tld| tld.hit_params.as_ref())
        .map_or(0, |hit_params| hit_params.options.hitlist_size);
    let mut consolidated = blast_hsp_results_new(num_queries);

    for query_idx in 0..counts.len() {
        let mut hits4query = blast_hit_list_new(hitlist_size);
        hits4query.hsp_lists.reserve(counts[query_idx]);
        let mut stats_initialized = false;

        for tld in array
            .tld
            .iter_mut()
            .take(array.num_elems as usize)
            .flatten()
        {
            let Some(results) = tld.results.as_mut() else {
                continue;
            };
            let Some(Some(thread_hitlist)) = results.hitlists.get_mut(query_idx) else {
                continue;
            };

            for hsp_list in thread_hitlist.hsp_lists.drain(..) {
                if !hsp_list.hsps.is_empty() {
                    hits4query.hsp_lists.push(hsp_list);
                }
            }

            if !stats_initialized {
                hits4query.worst_evalue = thread_hitlist.worst_evalue;
                hits4query.low_score = thread_hitlist.low_score;
                stats_initialized = true;
            } else {
                hits4query.worst_evalue = hits4query.worst_evalue.max(thread_hitlist.worst_evalue);
                hits4query.low_score = hits4query.low_score.min(thread_hitlist.low_score);
            }
        }

        consolidated.hitlists[query_idx] = Some(hits4query);
    }
    Some(consolidated)
}

/// Port of NCBI `Blast_HSPStreamResultBatchInit` (`blast_hspstream.c:616`).
pub fn blast_hsp_stream_result_batch_init(num_hsplists: i32) -> BlastHspStreamResultBatch {
    BlastHspStreamResultBatch {
        hsplist_array: (0..num_hsplists.max(0)).map(|_| None).collect(),
        num_hsplists: 0,
    }
}

/// Rust ownership equivalent of NCBI `Blast_HSPStreamResultBatchFree`.
pub fn blast_hsp_stream_result_batch_free(
    mut batch: Option<BlastHspStreamResultBatch>,
) -> Option<BlastHspStreamResultBatch> {
    if let Some(batch) = batch.as_mut() {
        batch.hsplist_array.clear();
        batch.num_hsplists = 0;
    }
    None
}

/// Port of NCBI `Blast_HSPStreamResultBatchReset` (`blast_hspstream.c:639`).
#[allow(non_snake_case)]
pub fn Blast_HSPStreamResultBatchReset(
    batch: Option<&mut BlastHspStreamResultBatch>,
) -> Option<&mut BlastHspStreamResultBatch> {
    let batch = batch?;
    let num_hsplists = batch.num_hsplists.max(0) as usize;
    for slot in batch.hsplist_array.iter_mut().take(num_hsplists) {
        let _ = blast_hsp_list_free(slot);
    }
    batch.num_hsplists = 0;
    Some(batch)
}

/// Port of NCBI `s_BlastHSPStreamResultsBatchArrayReset`
/// (`blast_hspstream_mt_utils.c:42`).
pub fn s_blast_hsp_stream_results_batch_array_reset(
    batches: &mut Option<BlastHspStreamResultsBatchArray>,
) {
    let Some(batches) = batches.as_mut() else {
        return;
    };
    for slot in batches
        .array_of_batches
        .iter_mut()
        .take(batches.num_batches as usize)
    {
        let _ = Blast_HSPStreamResultBatchReset(slot.as_mut());
        let taken = slot.take();
        let _ = blast_hsp_stream_result_batch_free(taken);
    }
    batches.num_batches = 0;
}

/// Rust ownership equivalent of NCBI `BlastHSPStreamResultsBatchArrayFree`.
pub fn blast_hsp_stream_results_batch_array_free(
    mut batches: Option<BlastHspStreamResultsBatchArray>,
) -> Option<BlastHspStreamResultsBatchArray> {
    s_blast_hsp_stream_results_batch_array_reset(&mut batches);
    if let Some(batches) = batches.as_mut() {
        batches.array_of_batches.clear();
        batches.num_allocated = 0;
    }
    None
}

/// Port of NCBI `s_BlastHSPStreamResultsBatchArrayNew`
/// (`blast_hspstream_mt_utils.c:72`).
pub fn s_blast_hsp_stream_results_batch_array_new(
    num_elements: u32,
) -> BlastHspStreamResultsBatchArray {
    let num_elements = if num_elements == 0 { 10 } else { num_elements };
    BlastHspStreamResultsBatchArray {
        array_of_batches: (0..num_elements).map(|_| None).collect(),
        num_batches: 0,
        num_allocated: num_elements,
    }
}

/// Port of NCBI `BlastHSPStreamResultsBatchNew`.
pub fn blast_hsp_stream_results_batch_new() -> BlastHspStreamResultsBatchArray {
    s_blast_hsp_stream_results_batch_array_new(1)
}

/// Port of NCBI `s_BlastHSPStreamResultsBatchArrayAppend`
/// (`blast_hspstream_mt_utils.c:103`).
pub fn s_blast_hsp_stream_results_batch_array_append(
    batches: Option<&mut BlastHspStreamResultsBatchArray>,
    batch: Option<BlastHspStreamResultBatch>,
) -> i32 {
    const BLASTERR_INVALIDPARAM: i32 = -1;
    let Some(batches) = batches else {
        return BLASTERR_INVALIDPARAM;
    };
    let Some(batch) = batch else {
        return BLASTERR_INVALIDPARAM;
    };
    if batches.num_batches + 1 > batches.num_allocated {
        let new_allocated = batches.num_allocated.saturating_mul(2).max(1);
        batches
            .array_of_batches
            .resize(new_allocated as usize, None);
        batches.num_allocated = new_allocated;
    }
    let index = batches.num_batches as usize;
    batches.array_of_batches[index] = Some(batch);
    batches.num_batches += 1;
    0
}

/// 1-1 port of `Blast_HSPResultsInsertHSPList` (`blast_hits.c:3554`) for the
/// Rust `HspResults` shape.
///
/// The C function derives the query index from `hsp_list->query_index`; Rust
/// keeps that value at the call site, so it is supplied explicitly here.
pub fn blast_hsp_results_insert_hsp_list(
    results: &mut HspResults,
    query_index: i32,
    hsp_list: HspList,
    hitlist_size: i32,
) -> i32 {
    if hsp_list.hsps.is_empty() {
        return 0;
    }
    let idx = query_index as usize;
    if idx >= results.hitlists.len() {
        return 1;
    }
    let hitlist = results.hitlists[idx].get_or_insert_with(|| blast_hit_list_new(hitlist_size));
    hitlist.blast_hit_list_update(hsp_list)
}

/// 1-1 port of `Blast_HSPResultsSortByEvalue` (`blast_hits.c:3383`) for the
/// Rust `HspResults` shape.
pub fn blast_hsp_results_sort_by_evalue(results: &mut HspResults) -> i32 {
    for hitlist in results.hitlists.iter_mut().flatten() {
        let _ = blast_hit_list_sort_by_evalue(hitlist);
        s_blast_hit_list_purge(Some(hitlist));
    }
    0
}

/// NCBI: PHIBlast_HSPResultsSplit (blast_hits.c:3572).
/// naming: Historical Rust spelling keeps `phiblast` as one token.
///
/// C reads the pattern occurrence index from `hsp->pat_info->index`; Rust
/// stores the same PHI metadata in [`Hsp::pat_info`]. The split/clone/insert
/// and final sort mechanics mirror the C helper.
pub fn phiblast_hsp_results_split(
    results: Option<&HspResults>,
    num_patterns: i32,
) -> Option<Vec<Option<HspResults>>> {
    if num_patterns <= 0 {
        return None;
    }
    let num_patterns = num_patterns as usize;
    let mut phi_results = vec![None; num_patterns];

    let Some(results) = results else {
        return Some(phi_results);
    };
    let Some(Some(hit_list)) = results.hitlists.first() else {
        return Some(phi_results);
    };

    for hsp_list in &hit_list.hsp_lists {
        let mut hsplist_array: Vec<Option<HspList>> = vec![None; num_patterns];
        for hsp in &hsp_list.hsps {
            let Some(pattern_idx) = hsp.pat_info.map(|info| info.index) else {
                continue;
            };
            if pattern_idx >= num_patterns {
                continue;
            }
            let list = hsplist_array[pattern_idx].get_or_insert_with(|| blast_hsp_list_new(0));
            list.oid = hsp_list.oid;
            let _ = blast_hsp_list_save_hsp(list, hsp.clone());
        }

        for (pattern_idx, maybe_list) in hsplist_array.into_iter().enumerate() {
            if let Some(list) = maybe_list {
                let result =
                    phi_results[pattern_idx].get_or_insert_with(|| blast_hsp_results_new(1));
                let _ = blast_hsp_results_insert_hsp_list(result, 0, list, hit_list.hsplist_max);
            }
        }
    }

    for result in phi_results.iter_mut().flatten() {
        let _ = blast_hsp_results_sort_by_evalue(result);
    }

    Some(phi_results)
}

/// Port of NCBI `Blast_TracebackGetEncoding` (`blast_traceback.c:831`).
pub fn blast_traceback_get_encoding(
    program_number: crate::program::ProgramType,
) -> crate::seqsrc::SeqEncoding {
    if crate::program::blast_subject_is_protein(program_number)
        || crate::program::blast_subject_is_pssm(program_number)
    {
        crate::seqsrc::SeqEncoding::Protein
    } else if crate::program::blast_subject_is_translated(program_number) {
        crate::seqsrc::SeqEncoding::Ncbi4na
    } else {
        crate::seqsrc::SeqEncoding::Nucleotide
    }
}

/// 1-1 translation of `Blast_HSPResultsReverseSort` (`blast_hits.c:3404`).
pub fn blast_hsp_results_reverse_sort(results: &mut HspResults) -> i32 {
    for hitlist in results.hitlists.iter_mut().flatten() {
        if hitlist.hsp_lists.len() > 1 {
            hitlist
                .hsp_lists
                .sort_by(|a, b| evalue_compare_hsp_lists(b, a));
        }
        s_blast_hit_list_purge(Some(hitlist));
    }
    0
}

/// Port of NCBI `s_BlastPruneExtraHits` (`blast_traceback.c:890`).
pub fn s_blast_prune_extra_hits(results: &mut HspResults, hitlist_size: i32) {
    let hitlist_size = hitlist_size.max(0) as usize;
    for hitlist in results.hitlists.iter_mut().flatten() {
        if hitlist.hsp_lists.len() > hitlist_size {
            hitlist.hsp_lists.truncate(hitlist_size);
        }
    }
}

/// 1-1 translation of `Blast_HSPResultsReverseOrder` (`blast_hits.c:3420`).
pub fn blast_hsp_results_reverse_order(results: &mut HspResults) -> i32 {
    for hitlist in results.hitlists.iter_mut().flatten() {
        if hitlist.hsp_lists.len() > 1 {
            hitlist.hsp_lists.reverse();
        }
    }
    0
}

/// 1-1 translation of `Blast_HSPResultsApplyMasklevel` (`blast_hits.c:3467`)
/// for the Rust result and interval-tree representations. For each query,
/// all HSPs are drained into one array, sorted by decreasing raw score, then
/// reinserted only if their query range is not masklevel-covered by a
/// higher-scoring HSP already accepted for the same context.
pub fn blast_hsp_results_apply_masklevel(
    results: &mut HspResults,
    masklevel: i32,
    query_length: i32,
) -> i32 {
    let mut tree = crate::itree::blast_interval_tree_init(0, query_length + 1, 0, 0);

    for hitlist in results.hitlists.iter_mut().flatten() {
        let mut hsp_array: Vec<(usize, Hsp)> = Vec::new();
        for (list_idx, hsp_list) in hitlist.hsp_lists.iter_mut().enumerate() {
            for hsp in hsp_list.hsps.drain(..) {
                hsp_array.push((list_idx, hsp));
            }
            // NCBI `s_BlastGetBestEvalue` seeds with `(double)INT4_MAX`.
            hsp_list.best_evalue = i32::MAX as f64;
        }

        hsp_array.sort_by(|(_, a), (_, b)| b.score.cmp(&a.score));
        crate::itree::blast_interval_tree_reset(&mut tree);

        for (list_idx, hsp) in hsp_array {
            let query_index = hsp.context;
            if crate::itree::blast_interval_tree_masks_hsp(&tree, &hsp, query_index, masklevel) {
                continue;
            }
            crate::itree::blast_interval_tree_add_hsp(&hsp, &mut tree, query_index);
            let hsp_list = &mut hitlist.hsp_lists[list_idx];
            if hsp_list.hsps.is_empty() {
                hsp_list.best_evalue = hsp.evalue;
            }
            let _ = blast_hsp_list_save_hsp(hsp_list, hsp);
        }

        hitlist
            .hsp_lists
            .retain(|hsp_list| !hsp_list.hsps.is_empty());
        for hsp_list in &mut hitlist.hsp_lists {
            // NCBI `Blast_HSPResultsApplyMasklevel` (`blast_hits.c:3543`)
            // calls only `Blast_HSPListSortByScore` here, with no `best_evalue`
            // recompute — the field was already set on first re-insert
            // (line 3530: `hsplist->best_evalue = hsp->evalue;`).
            hsp_list.hsps.sort_by(score_compare_hsps);
        }
    }

    0
}

/// 1-1 translation of `s_TrimResultsByTotalHSPLimit`
/// (`blast_hits.c:3685`). The total HSP limit is applied independently per
/// query. HSP lists are considered in increasing HSP-count order, and each
/// list is truncated so that the first `i` lists keep at most `i*T/N` HSPs,
/// with at least one HSP per subject list.
pub fn s_trim_results_by_total_hsp_limit(results: &mut HspResults, total_hsp_limit: u32) -> bool {
    if total_hsp_limit == 0 {
        return false;
    }

    let mut hsp_limit_exceeded = false;
    for hitlist in results.hitlists.iter_mut().flatten() {
        let hsplist_count = hitlist.hsp_lists.len();
        if hsplist_count == 0 {
            continue;
        }

        let mut order: Vec<usize> = (0..hsplist_count).collect();
        order.sort_by_key(|&idx| hitlist.hsp_lists[idx].hsps.len());

        let mut tot_hsps = 0usize;
        let hsp_per_seq = (total_hsp_limit as usize / hsplist_count).max(1);
        for (subj_index, list_idx) in order.into_iter().enumerate() {
            let allowed_hsp_num = ((subj_index + 1) * hsp_per_seq).saturating_sub(tot_hsps);
            let hsp_list = &mut hitlist.hsp_lists[list_idx];
            if hsp_list.hsps.len() > allowed_hsp_num {
                // NCBI `s_TrimResultsByTotalHSPLimit` (`blast_hits.c:3722`)
                // only truncates the tail; no `best_evalue` recompute. The
                // head HSP (highest-score, lowest-evalue) is preserved.
                hsp_list.hsps.truncate(allowed_hsp_num);
                hsp_limit_exceeded = true;
            }
            tot_hsps += hsp_list.hsps.len();
        }
    }

    hsp_limit_exceeded
}

/// NCBI: s_EvalueCompareHSPs (blast_hits.c:1415).
/// naming: Rust keeps HSPs as a readable snake_case plural token.
pub fn evalue_compare_hsps(a: &Hsp, b: &Hsp) -> std::cmp::Ordering {
    let by_evalue = evalue_comp(a.evalue, b.evalue);
    if by_evalue != std::cmp::Ordering::Equal {
        return by_evalue;
    }
    score_compare_hsps(a, b)
}

/// 1-1 translation of `s_TrimResultsByTotalHSPLimitEx`
/// (`blast_hits.c:3767`) for Rust-owned results.
///
/// The extended variant differs from `s_TrimResultsByTotalHSPLimit`: when a
/// query exceeds the total HSP limit, C flattens all HSPs for that query,
/// keeps the best HSPs by e-value/score, then rebuilds subject lists grouped
/// by OID. `hsp_limit_exceeded`, when supplied, receives one flag per query.
pub fn s_trim_results_by_total_hsp_limit_ex(
    results: &mut HspResults,
    total_hsp_limit: u32,
    mut hsp_limit_exceeded: Option<&mut [bool]>,
) -> bool {
    if total_hsp_limit == 0 {
        return false;
    }

    let mut any_hsp_limit_exceeded = false;
    let keep_limit = total_hsp_limit as usize;

    for query_index in 0..results.hitlists.len() {
        if let Some(flags) = hsp_limit_exceeded.as_deref_mut() {
            if let Some(flag) = flags.get_mut(query_index) {
                *flag = false;
            }
        }

        let Some(hitlist) = results.hitlists[query_index].as_mut() else {
            continue;
        };
        let total_hsps: usize = hitlist
            .hsp_lists
            .iter()
            .map(|hsp_list| hsp_list.hsps.len())
            .sum();
        if total_hsps <= keep_limit {
            continue;
        }

        if let Some(flags) = hsp_limit_exceeded.as_deref_mut() {
            if let Some(flag) = flags.get_mut(query_index) {
                *flag = true;
            }
        }
        any_hsp_limit_exceeded = true;

        let max_hit_list_size = hitlist.hsplist_max;
        let mut everything: Vec<(i32, Hsp)> = Vec::with_capacity(total_hsps);
        for hsp_list in hitlist.hsp_lists.drain(..) {
            for hsp in hsp_list.hsps {
                everything.push((hsp_list.oid, hsp));
            }
        }
        everything.sort_by(|(_, a), (_, b)| evalue_compare_hsps(a, b));
        everything.truncate(keep_limit);
        everything.sort_by_key(|(oid, _)| *oid);

        let mut rebuilt = blast_hit_list_new(max_hit_list_size);
        let mut cursor = 0;
        while cursor < everything.len() {
            let oid = everything[cursor].0;
            let start = cursor;
            while cursor < everything.len() && everything[cursor].0 == oid {
                cursor += 1;
            }

            let mut list = blast_hsp_list_new((cursor - start) as i32);
            list.oid = oid;
            for (_, hsp) in everything[start..cursor].iter().cloned() {
                let _ = blast_hsp_list_save_hsp(&mut list, hsp);
            }
            let _ = rebuilt.blast_hit_list_update(list);
        }
        results.hitlists[query_index] = Some(rebuilt);
    }

    any_hsp_limit_exceeded
}

/// 1-1 translation of `Blast_TrimHSPListByMaxHsps` (`blast_hits.c:2049`)
/// for the Rust vector-backed HSP list.
pub fn blast_trim_hsp_list_by_max_hsps(
    hsp_list: Option<&mut HspList>,
    hit_options: &HitSavingOptions,
) -> i16 {
    let Some(hsp_list) = hsp_list else {
        return 0;
    };
    let hsp_max = hit_options.max_hsps_per_subject;
    if hsp_max == 0 || hsp_list.hsps.len() <= hsp_max as usize {
        return 0;
    }
    hsp_list.hsps.truncate(hsp_max.max(0) as usize);
    0
}

/// 1-1 translation of `Blast_HSPListReapByEvalue` (`blast_hits.c:1984`).
pub fn blast_hsp_list_reap_by_evalue(
    hsp_list: Option<&mut HspList>,
    hit_options: &HitSavingOptions,
) -> i16 {
    let Some(hsp_list) = hsp_list else {
        return 0;
    };
    // NCBI `Blast_HSPListReapByEvalue` (`blast_hits.c:1976`) does NOT
    // update `hsp_list->best_evalue` after filtering. Sister function
    // `Blast_HSPListReapByQueryCoverage` (`blast_hits.c:2010`) does
    // recompute. Match NCBI semantics here.
    hsp_list
        .hsps
        .retain(|hsp| hsp.evalue <= hit_options.expect_value);
    0
}

/// 1-1 translation of `Blast_HSPListPurgeHSPsWithCommonEndpoints`
/// (`blast_hits.c:2455`).
/// naming: Rust keeps HSPs as a readable snake_case plural token.
pub fn blast_hsp_list_purge_hsps_with_common_endpoints(
    program_number: crate::program::ProgramType,
    hsp_list: Option<&mut HspList>,
    purge: bool,
) -> i16 {
    let Some(hsp_list) = hsp_list else {
        return 0;
    };
    if hsp_list.hsps.len() <= 1 {
        return hsp_list.hsps.len() as i16;
    }
    if crate::program::blast_program_is_phi_blast(program_number) {
        return hsp_list.hsps.len() as i16;
    }
    let trim_blastn = program_number == crate::program::BLASTN && !purge;
    let mut c_return_count = hsp_list.hsps.len();

    fn trim_or_drop_common_endpoint(
        leader: &Hsp,
        hsp: &mut Hsp,
        cut_begin: bool,
        trim_blastn: bool,
    ) -> bool {
        if !trim_blastn || hsp.edit_script.is_none() {
            return false;
        }
        let before = (
            hsp.query_offset,
            hsp.query_end,
            hsp.subject_offset,
            hsp.subject_end,
        );
        let (q_cut, s_cut) = if cut_begin {
            (leader.query_end, leader.subject_end)
        } else {
            (leader.query_offset, leader.subject_offset)
        };
        s_cut_off_gap_edit_script(hsp, q_cut, s_cut, cut_begin);
        let after = (
            hsp.query_offset,
            hsp.query_end,
            hsp.subject_offset,
            hsp.subject_end,
        );
        before != after && hsp.query_offset < hsp.query_end && hsp.subject_offset < hsp.subject_end
    }

    // NCBI `s_QueryOffsetCompareHSPs` (`blast_hits.c:2268`) sorts by:
    //   context-asc, query.offset-asc, subject.offset-asc, score-desc,
    //   query.end-desc, subject.end-desc.
    // Missing `context` would cause HSPs from different query
    // contexts/strands to interleave, and the equality-check below would
    // then incorrectly merge cross-context pairs that share offsets.
    hsp_list
        .hsps
        .sort_by(|a, b| s_query_offset_compare_hsps(Some(a), Some(b)));
    // NCBI `blast_hits.c:2487-2492` requires matching context, query
    // offset, subject offset, AND subject frame before merging. Frame
    // matters for translated programs (blastx, tblastn, tblastx) where
    // different frames can share offset values.
    let mut active = std::mem::take(&mut hsp_list.hsps);
    let mut extras = Vec::new();
    let mut index = 0;
    while index < active.len() {
        while index + 1 < active.len()
            && active[index].context == active[index + 1].context
            && active[index].query_offset == active[index + 1].query_offset
            && active[index].subject_offset == active[index + 1].subject_offset
            && active[index].subject_frame == active[index + 1].subject_frame
        {
            c_return_count = c_return_count.saturating_sub(1);
            let leader = active[index].clone();
            let mut trailing = active.remove(index + 1);
            if trim_or_drop_common_endpoint(&leader, &mut trailing, true, trim_blastn) {
                extras.push(trailing);
            }
        }
        index += 1;
    }

    if active.len() > 1 {
        // NCBI `s_QueryEndCompareHSPs` (`blast_hits.c:2333`): context-asc,
        // query.end-asc, subject.end-asc, score-desc, query.offset-desc,
        // subject.offset-desc. Missing the context sort key would cause
        // cross-context interleaving and incorrect merges below.
        active.sort_by(|a, b| {
            a.context
                .cmp(&b.context)
                .then_with(|| a.query_end.cmp(&b.query_end))
                .then_with(|| a.subject_end.cmp(&b.subject_end))
                .then_with(|| b.score.cmp(&a.score))
                .then_with(|| b.query_offset.cmp(&a.query_offset))
                .then_with(|| b.subject_offset.cmp(&a.subject_offset))
        });
        // NCBI `blast_hits.c:2509-2514` matches on context, query.end,
        // subject.end, AND subject.frame for the second-pass merge.
        let mut index = 0;
        while index < active.len() {
            while index + 1 < active.len()
                && active[index].context == active[index + 1].context
                && active[index].query_end == active[index + 1].query_end
                && active[index].subject_end == active[index + 1].subject_end
                && active[index].subject_frame == active[index + 1].subject_frame
            {
                c_return_count = c_return_count.saturating_sub(1);
                let leader = active[index].clone();
                let mut trailing = active.remove(index + 1);
                if trim_or_drop_common_endpoint(&leader, &mut trailing, false, trim_blastn) {
                    extras.push(trailing);
                }
            }
            index += 1;
        }
    }

    active.append(&mut extras);
    hsp_list.hsps = active;
    hsp_list.best_evalue = s_blast_get_best_evalue(hsp_list);
    c_return_count as i16
}

/// `BlastHSPList`-shaped counterpart of
/// `Blast_HSPListPurgeHSPsWithCommonEndpoints` (`blast_hits.c:2455`).
///
/// This keeps the upstream `BlastHSP` fields intact, including `num`, so
/// translated traceback metadata survives the purge without converting through
/// the legacy Rust `HspList` adapter.
pub fn blast_hsp_list_purge_blast_hsps_with_common_endpoints(
    program_number: crate::program::ProgramType,
    hsp_list: Option<&mut BlastHSPList>,
    purge: bool,
) -> i16 {
    let Some(hsp_list) = hsp_list else {
        return 0;
    };
    let hspcnt = (hsp_list.hspcnt.max(0) as usize).min(hsp_list.hsp_array.len());
    if hspcnt == 0 {
        hsp_list.hsp_array.clear();
        hsp_list.hspcnt = 0;
        hsp_list.allocated = 0;
        hsp_list.best_evalue = i32::MAX as f64;
        return 0;
    }
    let active_count = hsp_list.hsp_array.iter().take(hspcnt).flatten().count();
    if active_count <= 1 {
        hsp_list.hsp_array = hsp_list
            .hsp_array
            .drain(..hspcnt)
            .flatten()
            .map(Some)
            .collect();
        hsp_list.hspcnt = active_count as i32;
        hsp_list.allocated = hsp_list.hspcnt;
        hsp_list.best_evalue = hsp_list
            .hsp_array
            .iter()
            .flatten()
            .map(|hsp| hsp.evalue)
            .fold(i32::MAX as f64, f64::min);
        return active_count as i16;
    }
    if crate::program::blast_program_is_phi_blast(program_number) {
        hsp_list.hsp_array = hsp_list
            .hsp_array
            .drain(..hspcnt)
            .flatten()
            .map(Some)
            .collect();
        hsp_list.hspcnt = active_count as i32;
        hsp_list.allocated = hsp_list.hspcnt;
        hsp_list.best_evalue = hsp_list
            .hsp_array
            .iter()
            .flatten()
            .map(|hsp| hsp.evalue)
            .fold(i32::MAX as f64, f64::min);
        return active_count as i16;
    }
    let trim_blastn = program_number == crate::program::BLASTN && !purge;
    debug_assert!(
        !trim_blastn,
        "BlastHSPList purge does not trim BLASTN traceback scripts"
    );
    let mut active: Vec<BlastHSP> = hsp_list.hsp_array.drain(..hspcnt).flatten().collect();
    let mut c_return_count = active.len();

    active.sort_by(|a, b| {
        a.context
            .cmp(&b.context)
            .then_with(|| a.query.offset.cmp(&b.query.offset))
            .then_with(|| a.subject.offset.cmp(&b.subject.offset))
            .then_with(|| b.score.cmp(&a.score))
            .then_with(|| b.query.end.cmp(&a.query.end))
            .then_with(|| b.subject.end.cmp(&a.subject.end))
    });
    let mut index = 0;
    while index < active.len() {
        while index + 1 < active.len()
            && active[index].context == active[index + 1].context
            && active[index].query.offset == active[index + 1].query.offset
            && active[index].subject.offset == active[index + 1].subject.offset
            && active[index].subject.frame == active[index + 1].subject.frame
        {
            c_return_count = c_return_count.saturating_sub(1);
            active.remove(index + 1);
        }
        index += 1;
    }

    if active.len() > 1 {
        active.sort_by(|a, b| {
            a.context
                .cmp(&b.context)
                .then_with(|| a.query.end.cmp(&b.query.end))
                .then_with(|| a.subject.end.cmp(&b.subject.end))
                .then_with(|| b.score.cmp(&a.score))
                .then_with(|| b.query.offset.cmp(&a.query.offset))
                .then_with(|| b.subject.offset.cmp(&a.subject.offset))
        });
        let mut index = 0;
        while index < active.len() {
            while index + 1 < active.len()
                && active[index].context == active[index + 1].context
                && active[index].query.end == active[index + 1].query.end
                && active[index].subject.end == active[index + 1].subject.end
                && active[index].subject.frame == active[index + 1].subject.frame
            {
                c_return_count = c_return_count.saturating_sub(1);
                active.remove(index + 1);
            }
            index += 1;
        }
    }

    hsp_list.best_evalue = active
        .iter()
        .map(|hsp| hsp.evalue)
        .fold(i32::MAX as f64, f64::min);
    hsp_list.hspcnt = active.len() as i32;
    hsp_list.allocated = hsp_list.hspcnt;
    hsp_list.hsp_array = active.into_iter().map(Some).collect();
    c_return_count as i16
}

/// 1-1 translation of `Blast_HSPListReapByQueryCoverage`
/// (`blast_hits.c:2010`) with query lengths supplied by context index.
pub fn blast_hsp_list_reap_by_query_coverage(
    hsp_list: Option<&mut HspList>,
    hit_options: &HitSavingOptions,
    query_lengths_by_context: &[i32],
) -> i16 {
    let Some(hsp_list) = hsp_list else {
        return 0;
    };
    if hsp_list.hsps.is_empty() || hit_options.query_cov_hsp_perc == 0.0 {
        return 0;
    }

    let old_len = hsp_list.hsps.len();
    hsp_list.hsps.retain(|hsp| {
        let query_length = query_lengths_by_context
            .get(hsp.context.max(0) as usize)
            .copied()
            .unwrap_or(0);
        !blast_hsp_query_coverage_test(hsp, hit_options.query_cov_hsp_perc, query_length)
    });
    if hsp_list.hsps.len() != old_len {
        hsp_list.best_evalue = s_blast_get_best_evalue(hsp_list);
    }
    0
}

/// 1-1 translation of `Blast_HSPListSubjectBestHit` (`blast_hits.c:2538`)
/// over Rust-owned HSP lists.
///
/// The C helper assumes HSPs are sorted by score and deletes later HSPs whose
/// query range is contained in a better HSP's query range, expanded by
/// `range_diff`. For BLASTN it repeats the check against the opposite query
/// strand context after flipping the better HSP's range through the query
/// length.
pub fn blast_hsp_list_subject_best_hit(
    program_number: crate::program::ProgramType,
    hsp_list: Option<&mut HspList>,
    range_diff: i32,
    query_lengths_by_context: &[i32],
) -> i32 {
    let Some(hsp_list) = hsp_list else {
        return 0;
    };
    if hsp_list.hsps.is_empty() {
        return 0;
    }
    if crate::program::blast_program_is_phi_blast(program_number) {
        return hsp_list.hsps.len() as i32;
    }

    let mut keep = vec![true; hsp_list.hsps.len()];
    for i in 0..hsp_list.hsps.len().saturating_sub(1) {
        if !keep[i] {
            continue;
        }
        let current = &hsp_list.hsps[i];
        let offset = (current.query_offset - range_diff).max(0);
        let mut end = current.query_end + range_diff;
        if end < 0 {
            end = current.query_end;
        }
        for j in (i + 1)..hsp_list.hsps.len() {
            let other = &hsp_list.hsps[j];
            if keep[j]
                && current.context == other.context
                && other.query_offset >= offset
                && other.query_end <= end
            {
                keep[j] = false;
            }
        }
    }

    hsp_list.hsps = hsp_list
        .hsps
        .drain(..)
        .zip(keep)
        .filter_map(|(hsp, keep)| keep.then_some(hsp))
        .collect();

    if program_number == crate::program::BLASTN {
        let mut keep = vec![true; hsp_list.hsps.len()];
        for i in 0..hsp_list.hsps.len().saturating_sub(1) {
            if !keep[i] {
                continue;
            }
            let current = &hsp_list.hsps[i];
            let current_context = current.context;
            let Some(query_len) = query_lengths_by_context
                .get(current_context.max(0) as usize)
                .copied()
            else {
                continue;
            };
            let target_context = if current.query_frame > 0 {
                current_context + 1
            } else {
                current_context - 1
            };
            let end = query_len - (current.query_offset - range_diff);
            let offset = query_len - (current.query_end + range_diff);
            for j in (i + 1)..hsp_list.hsps.len() {
                let other = &hsp_list.hsps[j];
                if keep[j]
                    && other.context == target_context
                    && other.query_offset >= offset
                    && other.query_end <= end
                {
                    keep[j] = false;
                }
            }
        }

        hsp_list.hsps = hsp_list
            .hsps
            .drain(..)
            .zip(keep)
            .filter_map(|(hsp, keep)| keep.then_some(hsp))
            .collect();
    }

    hsp_list.hsps.len() as i32
}

/// blast-rs: Port-shaped filter coordinator corresponding to `s_FilterBlastResults`; not a direct NCBI C port.
///
/// This applies the C max-HSP, query-coverage, and subject-best-hit branches
/// and purges empty subject lists after query-coverage filtering.
pub fn s_filter_blast_results(
    results: &mut HspResults,
    hit_options: &HitSavingOptions,
    query_lengths_by_context: &[i32],
    query_info: Option<&crate::queryinfo::QueryInfo>,
) {
    if let Some((best_hit_opts, query_info)) = hit_options
        .hsp_filt_opt
        .as_ref()
        .and_then(|opts| opts.best_hit.as_ref())
        .zip(query_info)
    {
        let params = crate::hspfilter_culling::blast_hsp_best_hit_params_new(
            hit_options,
            best_hit_opts,
            0,
            true,
        );
        let mut data = crate::hspfilter_culling::s_blast_hsp_best_hit_pipe_new(params, query_info);
        let _ = crate::hspfilter_culling::s_blast_hsp_best_hit_pipe_run(&mut data, results);
    }

    if let Some((culling_opts, query_info)) = hit_options
        .hsp_filt_opt
        .as_ref()
        .and_then(|opts| opts.culling_opts.as_ref())
        .zip(query_info)
    {
        let params = crate::hspfilter_culling::blast_hsp_culling_params_new(
            hit_options,
            culling_opts,
            0,
            hit_options.program_number,
            0,
            true,
        );
        let mut data = crate::hspfilter_culling::blast_hsp_culling_pipe_new(params, query_info);
        crate::hspfilter_culling::blast_hsp_culling_pipe_run(&mut data, results);
    }

    for hitlist in results.hitlists.iter_mut().flatten() {
        for hsp_list in &mut hitlist.hsp_lists {
            if hit_options.max_hsps_per_subject != 0 {
                let _ = blast_trim_hsp_list_by_max_hsps(Some(hsp_list), hit_options);
            }
            if hit_options.query_cov_hsp_perc != 0.0 {
                let _ = blast_hsp_list_reap_by_query_coverage(
                    Some(hsp_list),
                    hit_options,
                    query_lengths_by_context,
                );
            }
            if let Some(subject_besthit_opts) = hit_options
                .hsp_filt_opt
                .as_ref()
                .and_then(|opts| opts.subject_besthit_opts)
            {
                let _ = blast_hsp_list_subject_best_hit(
                    hit_options.program_number,
                    Some(hsp_list),
                    subject_besthit_opts.max_range_diff,
                    query_lengths_by_context,
                );
            }
        }
        if hit_options.query_cov_hsp_perc != 0.0 {
            s_blast_hit_list_purge(Some(hitlist));
        }
    }
}

/// Thread-safe HSP stream for collecting results during parallel search.
pub struct HspStream {
    program: crate::program::ProgramType,
    results: Mutex<HspResults>,
    sorted_hsplists: Mutex<Vec<SortedHspList>>,
    results_sorted: Mutex<bool>,
    sort_by_score: Mutex<Option<SortByScoreState>>,
    writer: Mutex<Option<BlastHSPWriter>>,
    writer_initialized: Mutex<bool>,
    writer_finalized: Mutex<bool>,
    closed: Mutex<bool>,
    pre_pipes: Mutex<Vec<BlastHSPPipe>>,
    tback_pipes: Mutex<Vec<BlastHSPPipe>>,
}

#[derive(Debug, Clone, Copy, Default)]
struct SortByScoreState {
    sort_on_read: bool,
    first_query_index: usize,
}

#[derive(Debug, Clone)]
struct SortedHspList {
    query_index: usize,
    hsp_list: HspList,
}

impl HspStream {
    /// blast-rs: Thread-safe Rust writer method for owned HSP lists; not a direct NCBI C port.
    pub fn blast_hspstream_write(&self, query_index: i32, hsp_list: HspList) -> i32 {
        if *self.results_sorted.lock().unwrap() {
            return K_BLAST_HSP_STREAM_ERROR;
        }

        let mut hsp_list = hsp_list;
        let mut writer = self.writer.lock().unwrap();
        if let Some(writer) = writer.as_mut() {
            if !*self.writer_initialized.lock().unwrap() {
                if let Some(init_fn) = writer.init_fn_ptr {
                    let _ = init_fn(writer.data.as_mut(), None);
                }
                *self.writer_initialized.lock().unwrap() = true;
            }
            if let Some(run_fn) = writer.run_fn_ptr {
                let mut blast_list = BlastHSPList::from_legacy_hsp_list(hsp_list);
                let status = run_fn(writer.data.as_mut(), &mut blast_list);
                if status != 0 {
                    return K_BLAST_HSP_STREAM_ERROR;
                }
                hsp_list = blast_list.into_legacy_hsp_list();
            }
        }
        drop(writer);

        let mut results = self.results.lock().unwrap();
        let insertion_status =
            blast_hsp_results_insert_hsp_list(&mut results, query_index, hsp_list, 0);
        if insertion_status != 0 {
            return insertion_status;
        }
        *self.closed.lock().unwrap() = false;
        insertion_status
    }

    /// blast-rs: Rust ownership extractor for collected stream results; not a direct NCBI C port.
    pub fn into_results(self) -> HspResults {
        let mut results = self.results.into_inner().unwrap();
        for sorted in self.sorted_hsplists.into_inner().unwrap() {
            if sorted.query_index >= results.hitlists.len() {
                continue;
            }
            let hitlist = results.hitlists[sorted.query_index].get_or_insert_with(HitList::new);
            hitlist.hsp_lists.push(sorted.hsp_list);
        }
        results
    }
}

/// NCBI: BlastHSPStreamNew (`blast_hspstream.c:653`).
pub fn blast_hsp_stream_new(
    program: crate::program::ProgramType,
    composition_based_stats: i32,
    sort_on_read: bool,
    num_queries: i32,
    writer: Option<BlastHSPWriter>,
) -> HspStream {
    let sort_by_score = if (crate::program::blast_query_is_protein(program)
        || crate::program::blast_query_is_pssm(program))
        && composition_based_stats != 0
    {
        Some(SortByScoreState {
            sort_on_read,
            first_query_index: 0,
        })
    } else {
        None
    };
    HspStream {
        program,
        results: Mutex::new(HspResults::new(num_queries)),
        sorted_hsplists: Mutex::new(Vec::with_capacity(100)),
        results_sorted: Mutex::new(false),
        sort_by_score: Mutex::new(sort_by_score),
        writer: Mutex::new(writer),
        writer_initialized: Mutex::new(false),
        writer_finalized: Mutex::new(false),
        closed: Mutex::new(false),
        pre_pipes: Mutex::new(Vec::new()),
        tback_pipes: Mutex::new(Vec::new()),
    }
}

/// blast-rs: Rust ownership equivalent of the C HSP-stream free routine; not a direct NCBI C port.
/// (`blast_hspstream.c:46`).
pub fn blast_hsp_stream_free(stream: &mut Option<HspStream>) -> Option<HspStream> {
    if let Some(stream) = stream.as_mut() {
        *stream.closed.lock().unwrap() = true;
        *stream.results_sorted.lock().unwrap() = true;
        stream.pre_pipes.lock().unwrap().clear();
        stream.tback_pipes.lock().unwrap().clear();
        stream.sorted_hsplists.lock().unwrap().clear();
        if let Some(writer) = stream.writer.lock().unwrap().take() {
            if let Some(free_fn) = writer.free_fn_ptr {
                let _ = free_fn(writer);
            }
        }
        let mut results = stream.results.lock().unwrap();
        results.hitlists.clear();
    }
    let _ = stream.take();
    None
}

/// NCBI: s_SortHSPListByOid (blast_hspstream.c:91).
/// naming: Rust spells HSP list as separate snake_case tokens.
pub fn s_sort_hsp_list_by_oid(left: &HspList, right: &HspList) -> i32 {
    right.oid - left.oid
}

fn s_apply_hsp_pipe_to_results(pipe: &mut BlastHSPPipe, results: &mut HspResults) -> i32 {
    let Some(run_fn) = pipe.run_fn_ptr else {
        return 0;
    };
    let mut blast_results = BlastHSPResults::from_legacy_hsp_results(results.clone());
    let status = run_fn(pipe.data.as_mut(), &mut blast_results);
    if status == 0 {
        *results = blast_results.into_legacy_hsp_results();
    }
    status
}

/// blast-rs: Rust lifecycle equivalent of static `s_FinalizeWriter`; not a direct NCBI C port.
pub fn s_finalize_writer(stream: Option<&HspStream>) {
    let Some(stream) = stream else {
        return;
    };
    if *stream.writer_finalized.lock().unwrap() {
        return;
    }
    {
        let mut writer = stream.writer.lock().unwrap();
        if let Some(writer) = writer.as_mut() {
            if !*stream.writer_initialized.lock().unwrap() {
                if let Some(init_fn) = writer.init_fn_ptr {
                    let _ = init_fn(writer.data.as_mut(), None);
                }
                *stream.writer_initialized.lock().unwrap() = true;
            }
            if let Some(final_fn) = writer.final_fn_ptr {
                let _ = final_fn(writer.data.as_mut(), None);
            }
        }
    }
    loop {
        let mut pipe = {
            let mut pipes = stream.pre_pipes.lock().unwrap();
            if pipes.is_empty() {
                break;
            }
            pipes.remove(0)
        };
        let mut results = stream.results.lock().unwrap();
        let _ = s_apply_hsp_pipe_to_results(&mut pipe, &mut results);
        drop(results);
        if let Some(free_fn) = pipe.free_fn_ptr {
            let _ = free_fn(pipe);
        }
    }
    *stream.writer_finalized.lock().unwrap() = true;
}

/// Port-shaped close corresponding to `BlastHSPStreamClose`.
pub fn blast_hsp_stream_close(stream: Option<&HspStream>) {
    let Some(stream) = stream else {
        return;
    };
    if *stream.results_sorted.lock().unwrap() {
        return;
    }
    s_finalize_writer(Some(stream));
    let mut results = stream.results.lock().unwrap();

    if let Some(sort_state) = *stream.sort_by_score.lock().unwrap() {
        if sort_state.sort_on_read {
            let _ = blast_hsp_results_reverse_sort(&mut results);
        } else {
            let _ = blast_hsp_results_reverse_order(&mut results);
        }
        *stream.results_sorted.lock().unwrap() = true;
        *stream.closed.lock().unwrap() = true;
        return;
    }

    let mut sorted_hsplists = stream.sorted_hsplists.lock().unwrap();
    for (query_index, hitlist) in results.hitlists.iter_mut().enumerate() {
        let Some(hitlist) = hitlist else {
            continue;
        };
        for hsp_list in hitlist.hsp_lists.drain(..) {
            sorted_hsplists.push(SortedHspList {
                query_index,
                hsp_list,
            });
        }
    }
    sorted_hsplists
        .sort_by(|left, right| s_sort_hsp_list_by_oid(&left.hsp_list, &right.hsp_list).cmp(&0));
    *stream.results_sorted.lock().unwrap() = true;
    *stream.closed.lock().unwrap() = true;
}

/// blast-rs: Port-shaped read corresponding to `BlastHSPStreamRead`; not a direct NCBI C port.
///
/// The C stream materializes a separate `sorted_hsplists` array during close
/// and then pops from its end. The Rust stream keeps ownership in
/// `HspResults`, so this removes the same lowest-OID list after the close-time
/// descending OID sort.
pub fn blast_hsp_stream_read(stream: Option<&HspStream>) -> (i32, Option<HspList>) {
    let Some(stream) = stream else {
        return (K_BLAST_HSP_STREAM_ERROR, None);
    };
    if !*stream.results_sorted.lock().unwrap() {
        blast_hsp_stream_close(Some(stream));
    }

    if stream.sort_by_score.lock().unwrap().is_some() {
        let mut sort_by_score = stream.sort_by_score.lock().unwrap();
        let Some(sort_state) = sort_by_score.as_mut() else {
            return (K_BLAST_HSP_STREAM_EOF, None);
        };
        let mut results = stream.results.lock().unwrap();
        let mut query_index = sort_state.first_query_index;
        while query_index < results.hitlists.len() {
            if results.hitlists[query_index]
                .as_ref()
                .is_some_and(|hitlist| !hitlist.hsp_lists.is_empty())
            {
                break;
            }
            query_index += 1;
        }
        if query_index >= results.hitlists.len() {
            return (K_BLAST_HSP_STREAM_EOF, None);
        }
        sort_state.first_query_index = query_index;
        let hitlist = results.hitlists[query_index]
            .as_mut()
            .expect("selected hitlist");
        let hsp_list = hitlist.hsp_lists.pop();
        if hitlist.hsp_lists.is_empty() {
            sort_state.first_query_index += 1;
            results.hitlists[query_index] = None;
        }
        return (K_BLAST_HSP_STREAM_SUCCESS, hsp_list);
    }

    let mut sorted_hsplists = stream.sorted_hsplists.lock().unwrap();
    let Some(sorted) = sorted_hsplists.pop() else {
        return (K_BLAST_HSP_STREAM_EOF, None);
    };
    (K_BLAST_HSP_STREAM_SUCCESS, Some(sorted.hsp_list))
}

fn blast_hsp_stream_batch_read_score_sorted(
    stream: &HspStream,
    batch: &mut BlastHspStreamResultBatch,
) -> i32 {
    let mut results = stream.results.lock().unwrap();
    let mut sort_by_score = stream.sort_by_score.lock().unwrap();
    let Some(sort_state) = sort_by_score.as_mut() else {
        return K_BLAST_HSP_STREAM_EOF;
    };
    let mut query_index = sort_state.first_query_index;
    while query_index < results.hitlists.len() {
        if results.hitlists[query_index]
            .as_ref()
            .is_some_and(|hitlist| !hitlist.hsp_lists.is_empty())
        {
            break;
        }
        query_index += 1;
    }
    if query_index >= results.hitlists.len() {
        return K_BLAST_HSP_STREAM_EOF;
    }
    sort_state.first_query_index = query_index;
    let hitlist = results.hitlists[query_index]
        .as_mut()
        .expect("selected hitlist");
    let Some(target_oid) = hitlist.hsp_lists.last().map(|list| list.oid) else {
        return K_BLAST_HSP_STREAM_EOF;
    };
    let mut out = Vec::new();
    while hitlist
        .hsp_lists
        .last()
        .is_some_and(|hsp_list| hsp_list.oid == target_oid)
    {
        out.push(hitlist.hsp_lists.pop().expect("checked last"));
    }
    out.reverse();
    if hitlist.hsp_lists.is_empty() {
        sort_state.first_query_index += 1;
        results.hitlists[query_index] = None;
    }
    fill_hsp_stream_batch(batch, out);
    K_BLAST_HSP_STREAM_SUCCESS
}

fn fill_hsp_stream_batch(batch: &mut BlastHspStreamResultBatch, out: Vec<HspList>) {
    if batch.hsplist_array.len() < out.len() {
        batch.hsplist_array.resize_with(out.len(), || None);
    }
    for (index, hsp_list) in out.into_iter().enumerate() {
        batch.hsplist_array[index] = Some(hsp_list);
    }
    batch.num_hsplists = batch
        .hsplist_array
        .iter()
        .take_while(|slot| slot.is_some())
        .count() as i32;
}

/// blast-rs: Port-shaped read corresponding to `BlastHSPStreamBatchRead`; not a direct NCBI C port.
/// (`blast_hspstream.c:568`).
///
/// Returns all HSP lists with the same target OID as the next stream read,
/// matching the C batching contract used by the MT traceback utilities.
pub fn blast_hsp_stream_batch_read(
    stream: Option<&HspStream>,
    batch: Option<&mut BlastHspStreamResultBatch>,
) -> i32 {
    let Some(stream) = stream else {
        return K_BLAST_HSP_STREAM_ERROR;
    };
    let Some(batch) = batch else {
        return K_BLAST_HSP_STREAM_ERROR;
    };
    if !*stream.results_sorted.lock().unwrap() {
        blast_hsp_stream_close(Some(stream));
    }
    let _ = Blast_HSPStreamResultBatchReset(Some(batch));

    if stream.sort_by_score.lock().unwrap().is_some() {
        return blast_hsp_stream_batch_read_score_sorted(stream, batch);
    }

    let mut sorted_hsplists = stream.sorted_hsplists.lock().unwrap();
    let target_oid = sorted_hsplists.last().map(|sorted| sorted.hsp_list.oid);
    let Some(target_oid) = target_oid else {
        return K_BLAST_HSP_STREAM_EOF;
    };

    let mut out = Vec::new();
    while sorted_hsplists
        .last()
        .is_some_and(|sorted| sorted.hsp_list.oid == target_oid)
    {
        out.push(sorted_hsplists.pop().expect("checked last").hsp_list);
    }

    out.reverse();
    fill_hsp_stream_batch(batch, out);
    K_BLAST_HSP_STREAM_SUCCESS
}

/// blast-rs: Rust set-based OID counter corresponding to `s_BlastHSPStreamCountNumOids`; not a direct NCBI C port.
pub fn s_blast_hsp_stream_count_num_oids(stream: Option<&HspStream>) -> u32 {
    let Some(stream) = stream else {
        return 0;
    };
    if !*stream.results_sorted.lock().unwrap() {
        blast_hsp_stream_close(Some(stream));
    }
    let mut oids = std::collections::BTreeSet::new();
    if stream.sort_by_score.lock().unwrap().is_some() {
        let results = stream.results.lock().unwrap();
        for hitlist in results.hitlists.iter().flatten() {
            for hsp_list in &hitlist.hsp_lists {
                oids.insert(hsp_list.oid);
            }
        }
    } else {
        let sorted_hsplists = stream.sorted_hsplists.lock().unwrap();
        for sorted in sorted_hsplists.iter() {
            oids.insert(sorted.hsp_list.oid);
        }
    }
    oids.len() as u32
}

/// blast-rs: Port-shaped batch conversion corresponding to `BlastHSPStreamToHSPStreamResultsBatch`; not a direct NCBI C port.
/// (`blast_hspstream_mt_utils.c:145`).
pub fn blast_hsp_stream_to_hsp_stream_results_batch(
    stream: Option<&HspStream>,
) -> (i32, Option<BlastHspStreamResultsBatchArray>) {
    const BLASTERR_INVALIDPARAM: i32 = -1;
    const BLASTERR_MEMORY: i32 = -2;
    let Some(stream) = stream else {
        return (BLASTERR_INVALIDPARAM, None);
    };

    let num_oids = s_blast_hsp_stream_count_num_oids(Some(stream));
    let mut batches = s_blast_hsp_stream_results_batch_array_new(num_oids);
    let num_queries = stream.results.lock().unwrap().hitlists.len() as i32;

    loop {
        let mut batch = blast_hsp_stream_result_batch_init(num_queries);
        match blast_hsp_stream_batch_read(Some(stream), Some(&mut batch)) {
            K_BLAST_HSP_STREAM_SUCCESS => {
                if s_blast_hsp_stream_results_batch_array_append(Some(&mut batches), Some(batch))
                    != 0
                {
                    let mut slot = Some(batches);
                    s_blast_hsp_stream_results_batch_array_reset(&mut slot);
                    let _ = blast_hsp_stream_results_batch_array_free(slot);
                    return (BLASTERR_MEMORY, None);
                }
            }
            K_BLAST_HSP_STREAM_EOF => break,
            status => return (status, None),
        }
    }

    (K_BLAST_HSP_STREAM_SUCCESS, Some(batches))
}

/// blast-rs: Port-shaped close corresponding to `BlastHSPStreamSimpleClose`; not a direct NCBI C port.
pub fn blast_hsp_stream_simple_close(stream: Option<&HspStream>) {
    blast_hsp_stream_close(stream);
}

/// blast-rs: Port-shaped close corresponding to `BlastHSPStreamMappingClose`; not a direct NCBI C port.
///
/// Audit note: Rust's stream does not carry C's mapper writer/vtable payload,
/// so this only coordinates the existing Rust stream close path.
pub fn blast_hsp_stream_mapping_close(stream: Option<&HspStream>) {
    let Some(stream) = stream else {
        return;
    };
    s_finalize_writer(Some(stream));
    stream.pre_pipes.lock().unwrap().clear();
    stream.tback_pipes.lock().unwrap().clear();
    blast_hsp_stream_close(Some(stream));
}

/// Port of NCBI `BlastHSPStreamTBackClose` (`blast_hspstream.c:241`).
/// naming: Rust keeps TBack as the established `tback` token.
///
/// Applies traceback-stage pipes to traceback results. Like the C function,
/// this does not close or sort the HSP stream.
#[allow(non_snake_case)]
pub fn BlastHSPStreamTBackClose(stream: Option<&HspStream>, results: Option<&mut HspResults>) {
    let Some(stream) = stream else {
        return;
    };
    let Some(results) = results else {
        return;
    };
    loop {
        let mut pipe = {
            let mut pipes = stream.tback_pipes.lock().unwrap();
            if pipes.is_empty() {
                break;
            }
            pipes.remove(0)
        };
        let _ = s_apply_hsp_pipe_to_results(&mut pipe, results);
        if let Some(free_fn) = pipe.free_fn_ptr {
            let _ = free_fn(pipe);
        }
    }
}

/// blast-rs: HSP-list to query-index mapping helper for Rust results insertion; not a direct NCBI C port.
fn query_index_for_hsp_list(hsp_list: &HspList, query_info: &crate::queryinfo::QueryInfo) -> i32 {
    hsp_list
        .hsps
        .first()
        .and_then(|hsp| query_info.contexts.get(hsp.context.max(0) as usize))
        .map_or(0, |ctx| ctx.query_index)
}

/// Top-level traceback route used by NCBI `BLAST_ComputeTraceback_MT`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TracebackDispatchKind {
    Rps,
    RedoAlignment,
    Phi,
    Ordinary,
}

/// blast-rs: Branch classifier extracted from traceback dispatch; not a direct
/// NCBI C port.
///
/// C checks RPS first, then composition-based/Smith-Waterman redo alignment,
/// then dispatches per-HSP-list PHI traceback inside the ordinary stream loop.
pub fn blast_traceback_dispatch_kind(
    program_number: crate::program::ProgramType,
    composition_based_stats: i32,
    smith_waterman_traceback: bool,
) -> TracebackDispatchKind {
    if crate::program::blast_program_is_rps_blast(program_number) {
        TracebackDispatchKind::Rps
    } else if composition_based_stats > 0 || smith_waterman_traceback {
        TracebackDispatchKind::RedoAlignment
    } else if crate::program::blast_program_is_phi_blast(program_number) {
        TracebackDispatchKind::Phi
    } else {
        TracebackDispatchKind::Ordinary
    }
}

/// blast-rs: Port-shaped Rust equivalent of `BLAST_ComputeTraceback`
/// (`blast_traceback.c:1390`) for the represented stream/list layer.
/// not a direct NCBI C port.
///
/// The C function owns SeqSrc fetching and program-specific traceback routing.
/// Rust keeps those lower-level translations as typed call sites, so this
/// public driver coordinates the shared stage behavior: mark progress, drain
/// the HSP stream by subject batches, run a caller-supplied traceback callback
/// for each HSP list, insert surviving lists into `HspResults`, apply
/// mask/filter/prune policy, sort database-search results, and honor an
/// interrupt callback.
pub fn blast_compute_traceback<T>(
    program_number: crate::program::ProgramType,
    hsp_stream: Option<&HspStream>,
    query_info: Option<&crate::queryinfo::QueryInfo>,
    hit_options: &HitSavingOptions,
    query_length: i32,
    database_search: bool,
    mut traceback_hsp_list: T,
    mut interrupt_search: Option<&mut dyn FnMut(&crate::util::SBlastProgress) -> bool>,
    mut progress_info: Option<&mut crate::util::SBlastProgress>,
) -> (i16, Option<HspResults>)
where
    T: FnMut(&mut HspList) -> i16,
{
    let Some(stream) = hsp_stream else {
        return (crate::util::BLASTERR_INVALIDPARAM, None);
    };
    let Some(query_info) = query_info else {
        return (crate::util::BLASTERR_INVALIDPARAM, None);
    };

    if let Some(progress) = progress_info.as_deref_mut() {
        progress.stage = crate::util::EBlastStage::TracebackSearch;
    }

    let mut results = HspResults::new(query_info.num_queries);
    let (batch_status, batches) = blast_hsp_stream_to_hsp_stream_results_batch(Some(stream));
    if batch_status != K_BLAST_HSP_STREAM_SUCCESS {
        BlastHSPStreamTBackClose(Some(stream), Some(&mut results));
        return (batch_status as i16, None);
    }

    let mut interrupted = false;

    if let Some(batches) = batches {
        for mut batch in batches.array_of_batches.into_iter().flatten() {
            if let (Some(interrupt), Some(progress)) =
                (interrupt_search.as_deref_mut(), progress_info.as_deref())
            {
                if interrupt(progress) {
                    interrupted = true;
                    break;
                }
            }

            for mut hsp_list in batch
                .hsplist_array
                .drain(..)
                .take(batch.num_hsplists.max(0) as usize)
                .flatten()
            {
                let status = traceback_hsp_list(&mut hsp_list);
                if status != 0 {
                    BlastHSPStreamTBackClose(Some(stream), Some(&mut results));
                    return (status, None);
                }
                if !hsp_list.hsps.is_empty() {
                    let query_index = query_index_for_hsp_list(&hsp_list, query_info);
                    let insert_status = blast_hsp_results_insert_hsp_list(
                        &mut results,
                        query_index,
                        hsp_list,
                        hit_options.hitlist_size,
                    );
                    if insert_status != 0 {
                        BlastHSPStreamTBackClose(Some(stream), Some(&mut results));
                        return (insert_status as i16, None);
                    }
                }
            }
        }
    }

    BlastHSPStreamTBackClose(Some(stream), Some(&mut results));

    if interrupted {
        return (crate::diagnostics::BLASTERR_INTERRUPTED, None);
    }

    if hit_options.mask_level < 101 {
        let _ =
            blast_hsp_results_apply_masklevel(&mut results, hit_options.mask_level, query_length);
    }
    if hit_options.query_cov_hsp_perc > 0.0
        || hit_options.max_hsps_per_subject > 0
        || hit_options
            .hsp_filt_opt
            .as_ref()
            .and_then(|opts| {
                opts.best_hit
                    .as_ref()
                    .map(|_| ())
                    .or(opts.culling_opts.map(|_| ()))
                    .or(opts.subject_besthit_opts.map(|_| ()))
            })
            .is_some()
    {
        let query_lengths: Vec<i32> = query_info
            .contexts
            .iter()
            .map(|ctx| ctx.query_length)
            .collect();
        s_filter_blast_results(&mut results, hit_options, &query_lengths, Some(query_info));
    }
    if database_search {
        let _ = blast_hsp_results_sort_by_evalue(&mut results);
    }
    s_blast_prune_extra_hits(&mut results, hit_options.hitlist_size);

    let _ = program_number;
    (0, Some(results))
}

/// NCBI: Blast_RunTracebackSearch (`blast_traceback.c:1801`).
#[allow(clippy::too_many_arguments)]
pub fn blast_run_traceback_search<T>(
    program_number: crate::program::ProgramType,
    hsp_stream: Option<&HspStream>,
    query_info: Option<&crate::queryinfo::QueryInfo>,
    hit_options: &HitSavingOptions,
    query_length: i32,
    database_search: bool,
    use_cbs_close: bool,
    traceback_hsp_list: T,
    num_threads: usize,
) -> (i16, Option<HspResults>)
where
    T: FnMut(&mut HspList) -> i16,
{
    blast_run_traceback_search_with_interrupt(
        program_number,
        hsp_stream,
        query_info,
        hit_options,
        query_length,
        database_search,
        use_cbs_close,
        traceback_hsp_list,
        None,
        None,
        num_threads,
    )
}

/// NCBI: Blast_RunTracebackSearchWithInterrupt (`blast_traceback.c:1821`).
/// Source boundary: blast_traceback.c:1821.
///
/// This wrapper closes the preliminary stream before delegating to
/// [`blast_compute_traceback`], matching the C entry point's ownership
/// transition from preliminary search to traceback search.
#[allow(clippy::too_many_arguments)]
pub fn blast_run_traceback_search_with_interrupt<T>(
    program_number: crate::program::ProgramType,
    hsp_stream: Option<&HspStream>,
    query_info: Option<&crate::queryinfo::QueryInfo>,
    hit_options: &HitSavingOptions,
    query_length: i32,
    database_search: bool,
    use_cbs_close: bool,
    mut traceback_hsp_list: T,
    interrupt_search: Option<&mut dyn FnMut(&crate::util::SBlastProgress) -> bool>,
    progress_info: Option<&mut crate::util::SBlastProgress>,
    num_threads: usize,
) -> (i16, Option<HspResults>)
where
    T: FnMut(&mut HspList) -> i16,
{
    let Some(stream) = hsp_stream else {
        return (crate::util::BLASTERR_INVALIDPARAM, None);
    };
    if query_info.is_none() {
        return (crate::util::BLASTERR_INVALIDPARAM, None);
    }
    let _num_threads = if num_threads == 0 { 1 } else { num_threads };
    if use_cbs_close {
        blast_hsp_cbs_stream_close(Some(stream), hit_options.hitlist_size);
    } else {
        blast_hsp_stream_close(Some(stream));
    }

    blast_compute_traceback(
        program_number,
        Some(stream),
        query_info,
        hit_options,
        query_length,
        database_search,
        |hsp_list| traceback_hsp_list(hsp_list),
        interrupt_search,
        progress_info,
    )
}

/// blast-rs: Context-to-frame helper for split-query HSP stream merge; not a direct NCBI C port.
fn blast_context_to_frame_for_context(context: i32, contexts_per_query: i32) -> i32 {
    match contexts_per_query {
        2 => {
            if context.rem_euclid(2) == 0 {
                1
            } else {
                -1
            }
        }
        6 => crate::util::blast_context_to_frame(context.rem_euclid(6) as u32),
        _ => 0,
    }
}

/// blast-rs: Port-shaped stream merge corresponding to `BlastHSPStreamMerge`; not a direct NCBI C port.
/// (`blast_hspstream.c:399`) for Rust-owned streams. `stream1` is the
/// local split-query chunk and `stream2` is the accumulated global stream.
pub fn blast_hsp_stream_merge(
    squery_blk: Option<&SSplitQueryBlk>,
    chunk_num: u32,
    stream1: Option<&HspStream>,
    stream2: Option<&HspStream>,
    contexts_per_query: i32,
) -> i32 {
    let (Some(stream1), Some(stream2)) = (stream1, stream2) else {
        return -1;
    };
    s_finalize_writer(Some(stream1));
    s_finalize_writer(Some(stream2));

    if squery_blk.is_none() {
        let src = stream1.results.lock().unwrap();
        let mut dst = stream2.results.lock().unwrap();
        if dst.hitlists.len() < src.hitlists.len() {
            dst.hitlists.resize_with(src.hitlists.len(), || None);
        }
        for (query_index, src_hitlist) in src.hitlists.iter().enumerate() {
            let Some(src_hitlist) = src_hitlist else {
                continue;
            };
            let mut incoming = Some(src_hitlist.clone());
            let _ = blast_hit_list_merge_simple(&mut incoming, &mut dst.hitlists[query_index]);
        }
        return 0;
    }

    if contexts_per_query <= 0 {
        return -1;
    }

    let (_, Some(query_list)) = split_query_blk_get_query_indices_for_chunk(squery_blk, chunk_num)
    else {
        return -1;
    };
    let (_, Some((context_list, num_contexts))) =
        split_query_blk_get_query_contexts_for_chunk(squery_blk, chunk_num)
    else {
        return -1;
    };
    let (_, Some(offset_list)) =
        split_query_blk_get_context_offsets_for_chunk(squery_blk, chunk_num)
    else {
        return -1;
    };
    let chunk_overlap_size = split_query_blk_get_chunk_overlap_size(squery_blk) as i32;
    let allow_gap = split_query_blk_allow_gap(squery_blk);

    let src = stream1.results.lock().unwrap();
    let mut dst = stream2.results.lock().unwrap();
    let global_query_max = query_list
        .iter()
        .copied()
        .take_while(|&query| query != u32::MAX)
        .max()
        .map(|query| query as usize + 1)
        .unwrap_or(0);
    if dst.hitlists.len() < global_query_max {
        dst.hitlists.resize_with(global_query_max, || None);
    }

    for (local_query, src_hitlist) in src.hitlists.iter().enumerate() {
        let Some(global_query) = query_list.get(local_query).copied() else {
            break;
        };
        if global_query == u32::MAX {
            break;
        }
        let Some(src_hitlist) = src_hitlist.as_ref() else {
            continue;
        };

        let mut split_points = vec![-1; contexts_per_query as usize];
        for frame_index in 0..contexts_per_query as usize {
            let local_context = local_query * contexts_per_query as usize + frame_index;
            let Some(&global_context) = context_list.get(local_context) else {
                continue;
            };
            if global_context < 0 {
                continue;
            }
            let offset_idx = global_context.rem_euclid(contexts_per_query) as usize;
            if let Some(&offset) = offset_list.get(local_context) {
                if offset != u32::MAX {
                    split_points[offset_idx] = offset as i32;
                }
            }
        }

        let mut adjusted_hitlist = src_hitlist.clone();
        for hsplist in &mut adjusted_hitlist.hsp_lists {
            for hsp in &mut hsplist.hsps {
                let local_context = hsp.context.max(0) as usize;
                if local_context >= num_contexts as usize {
                    continue;
                }
                let Some(&global_context) = context_list.get(local_context) else {
                    continue;
                };
                let Some(&offset) = offset_list.get(local_context) else {
                    continue;
                };
                if global_context < 0 || offset == u32::MAX {
                    continue;
                }
                let offset = offset as i32;
                hsp.context = global_context;
                hsp.query_offset += offset;
                hsp.query_end += offset;
                hsp.query_gapped_start += offset;
                hsp.query_frame =
                    blast_context_to_frame_for_context(hsp.context, contexts_per_query);
            }
        }

        let global_query = global_query as usize;
        if dst.hitlists.len() <= global_query {
            dst.hitlists.resize_with(global_query + 1, || None);
        }
        let mut incoming = Some(adjusted_hitlist);
        let _ = blast_hit_list_merge(
            &mut incoming,
            &mut dst.hitlists[global_query],
            contexts_per_query,
            Some(&split_points),
            chunk_overlap_size,
            allow_gap,
        );
    }

    for hitlist in dst.hitlists.iter_mut().flatten() {
        for hsplist in &mut hitlist.hsp_lists {
            hsplist.hsps.sort_by(score_compare_hsps);
        }
    }
    0
}

/// blast-rs: Convenience wrapper around `blast_hsp_stream_merge`; not a direct NCBI C port.
pub fn blast_hsp_stream_merge_simple(
    stream1: Option<&HspStream>,
    stream2: Option<&HspStream>,
) -> i32 {
    blast_hsp_stream_merge(None, 0, stream2, stream1, 0)
}

/// blast-rs: Port-shaped status helper corresponding to `BlastHSPStreamRegisterMTLock`; not a direct NCBI C port.
pub fn blast_hsp_stream_register_mt_lock(stream: Option<&HspStream>, has_lock: bool) -> i32 {
    if stream.is_none() || !has_lock {
        -1
    } else {
        0
    }
}

#[derive(Debug, Clone, Default)]
pub struct BlastHSPPipe {
    pub data: Option<BlastHSPPipeData>,
    pub run_fn_ptr: Option<BlastHSPPipeRunFn>,
    pub free_fn_ptr: Option<BlastHSPPipeFreeFn>,
    pub next: Option<Box<BlastHSPPipe>>,
}

#[derive(Debug, Clone, Default)]
pub struct BlastHSPWriter {
    pub data: Option<BlastHSPWriterData>,
    pub init_fn_ptr: Option<BlastHSPWriterInitFn>,
    pub run_fn_ptr: Option<BlastHSPWriterRunFn>,
    pub final_fn_ptr: Option<BlastHSPWriterFinalFn>,
    pub free_fn_ptr: Option<BlastHSPWriterFreeFn>,
}

#[derive(Debug, Clone, Default)]
pub struct BlastHSPWriterInfo {
    pub params: Option<BlastHSPWriterData>,
    pub new_fn_ptr: Option<BlastHSPWriterNewFn>,
}

#[derive(Debug, Clone, Default)]
pub struct BlastHSPPipeInfo {
    pub params: Option<BlastHSPPipeData>,
    pub new_fn_ptr: Option<BlastHSPPipeNewFn>,
    pub next: Option<Box<BlastHSPPipeInfo>>,
}

#[derive(Debug, Clone, Default)]
pub enum BlastHSPWriterData {
    #[default]
    None,
    HspResults(Box<BlastHSPResults>),
    Opaque,
}

#[derive(Debug, Clone, Default)]
pub enum BlastHSPPipeData {
    #[default]
    None,
    Opaque,
}

pub type BlastHSPWriterNewFn = fn(Option<BlastHSPWriterData>) -> BlastHSPWriter;
pub type BlastHSPWriterFreeFn = fn(BlastHSPWriter) -> Option<BlastHSPWriter>;
pub type BlastHSPWriterInitFn = fn(Option<&mut BlastHSPWriterData>, Option<&mut ()>) -> i32;
pub type BlastHSPWriterRunFn = fn(Option<&mut BlastHSPWriterData>, &mut BlastHSPList) -> i32;
pub type BlastHSPWriterFinalFn = fn(Option<&mut BlastHSPWriterData>, Option<&mut ()>) -> i32;
pub type BlastHSPPipeNewFn = fn(Option<BlastHSPPipeData>) -> BlastHSPPipe;
pub type BlastHSPPipeFreeFn = fn(BlastHSPPipe) -> Option<BlastHSPPipe>;
pub type BlastHSPPipeRunFn = fn(Option<&mut BlastHSPPipeData>, &mut BlastHSPResults) -> i32;

/// blast-rs: Port-shaped status helper corresponding to `BlastHSPStreamRegisterPipe`; not a direct NCBI C port.
pub fn blast_hsp_stream_register_pipe(
    stream: Option<&HspStream>,
    pipe: Option<BlastHSPPipe>,
    stage: i32,
) -> i32 {
    let Some(stream) = stream else {
        return -1;
    };
    let Some(pipe) = pipe else {
        return -1;
    };
    match stage {
        E_PRELIM_SEARCH => {
            stream.pre_pipes.lock().unwrap().push(pipe);
            0
        }
        E_TRACEBACK_SEARCH => {
            stream.tback_pipes.lock().unwrap().push(pipe);
            0
        }
        _ => -1,
    }
}

/// blast-rs: Port-shaped constructor corresponding to `BlastHSPWriterNew`; not a direct NCBI C port.
pub fn blast_hsp_writer_new() -> BlastHSPWriter {
    BlastHSPWriter::default()
}

/// blast-rs: Port-shaped constructor corresponding to `BlastHSPPipeNew`; not a direct NCBI C port.
pub fn blast_hsp_pipe_new() -> BlastHSPPipe {
    BlastHSPPipe::default()
}

/// blast-rs: Test/introspection helper for the Rust stream's pipe-chain equivalent; not a direct NCBI C port.
pub fn blast_hsp_stream_pipe_counts(stream: Option<&HspStream>) -> Option<(usize, usize)> {
    let stream = stream?;
    Some((
        stream.pre_pipes.lock().unwrap().len(),
        stream.tback_pipes.lock().unwrap().len(),
    ))
}

/// NCBI: s_TrimHitList (blast_hspstream.c:806).
pub fn s_trim_hit_list(hitlist: Option<&mut HitList>, count: i32) {
    let Some(hitlist) = hitlist else {
        return;
    };
    let count = count.max(0) as usize;
    if hitlist.hsp_lists.len() > count {
        hitlist.hsp_lists.truncate(count);
    }
}

/// Port of NCBI `BlastHSPCBSStreamClose` (`blast_hspstream.c:816`).
/// naming: Rust spells HSP/CBS as separate snake_case tokens.
pub fn blast_hsp_cbs_stream_close(stream: Option<&HspStream>, hitlist_size: i32) {
    let Some(stream) = stream else {
        return;
    };
    if *stream.closed.lock().unwrap() {
        return;
    }
    s_finalize_writer(Some(stream));
    {
        let mut results = stream.results.lock().unwrap();
        for hitlist in results.hitlists.iter_mut().flatten() {
            let ref_hit_num = 500.max(hitlist_size).max(0) as usize;
            let min_buf_size = ref_hit_num + 600;
            if min_buf_size + 100 >= hitlist.hsp_lists.len() {
                continue;
            }

            let _ = blast_hit_list_sort_by_evalue(hitlist);
            let best_evalue = hitlist.hsp_lists[ref_hit_num].best_evalue;
            let mut mag = -180_i32;
            if best_evalue != 0.0 {
                mag = best_evalue.log10() as i32;
            }
            if mag < -20 {
                mag = (mag * 90 / 100).max(mag + 10);
            } else {
                mag /= 2;
            }
            let evalue_limit = if mag >= 0 {
                best_evalue * 3.0
            } else {
                9.9 * 10_f64.powi(mag)
            };

            let max_index = hitlist.hsp_lists.len().saturating_sub(1);
            let mut index = min_buf_size;
            while index < max_index {
                let evalue = hitlist.hsp_lists[index].best_evalue;
                if evalue != 0.0 && evalue_limit < evalue {
                    s_trim_hit_list(Some(hitlist), index as i32);
                    break;
                }
                index += 100;
            }
        }
    }
    blast_hsp_stream_close(Some(stream));
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Mutex as TestMutex;

    static ENV_LOCK: TestMutex<()> = TestMutex::new(());

    struct RecordingSeqSrc {
        seqs: Vec<Vec<u8>>,
        encodings: TestMutex<Vec<crate::seqsrc::SeqEncoding>>,
    }

    impl crate::seqsrc::BlastSeqSource for RecordingSeqSrc {
        fn num_seqs(&self) -> i32 {
            self.seqs.len() as i32
        }

        fn total_length(&self) -> i64 {
            self.seqs.iter().map(|seq| seq.len() as i64).sum()
        }

        fn max_seq_len(&self) -> i32 {
            self.seqs
                .iter()
                .map(|seq| seq.len() as i32)
                .max()
                .unwrap_or(0)
        }

        fn avg_seq_len(&self) -> i32 {
            if self.seqs.is_empty() {
                0
            } else {
                (self.total_length() / self.seqs.len() as i64) as i32
            }
        }

        fn name(&self) -> &str {
            "recording"
        }

        fn is_protein(&self) -> bool {
            true
        }

        fn seq_len(&self, oid: i32) -> i32 {
            self.seqs
                .get(oid.max(0) as usize)
                .map_or(0, |seq| seq.len() as i32)
        }

        fn get_sequence(&self, arg: &crate::seqsrc::GetSeqArg) -> Option<crate::seqsrc::SeqData> {
            self.encodings.lock().unwrap().push(arg.encoding);
            let seq = self.seqs.get(arg.oid as usize)?;
            Some(crate::seqsrc::SeqData {
                sequence: seq.clone(),
                length: seq.len() as i32,
            })
        }

        fn iter_oids(&self) -> Box<dyn Iterator<Item = i32> + '_> {
            Box::new(0..self.num_seqs())
        }
    }

    #[test]
    fn test_hsp_stream() {
        let stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 2, None);

        let mut list1 = HspList::new(0);
        list1.add_hsp(Hsp {
            score: 100,
            num_ident: 50,
            bit_score: 91.5,
            evalue: 1e-20,
            query_offset: 0,
            query_end: 50,
            query_gapped_start: 0,
            subject_offset: 76,
            subject_end: 126,
            subject_gapped_start: 76,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        stream.blast_hspstream_write(0, list1);

        let mut list2 = HspList::new(5);
        list2.add_hsp(Hsp {
            score: 30,
            num_ident: 12,
            bit_score: 24.3,
            evalue: 1.7,
            query_offset: 10,
            query_end: 22,
            query_gapped_start: 10,
            subject_offset: 100,
            subject_end: 112,
            subject_gapped_start: 100,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        });
        stream.blast_hspstream_write(0, list2);

        let results = stream.into_results();
        assert!(results.hitlists[0].is_some());
        assert!(results.hitlists[1].is_none());
        let hl = results.hitlists[0].as_ref().unwrap();
        assert_eq!(hl.hsp_lists.len(), 2);
        assert_eq!(hl.hsp_lists[0].oid, 0);
        assert_eq!(hl.hsp_lists[1].oid, 5);
    }

    /// Helper to create an HSP with given score and evalue.
    fn make_hsp(score: i32, evalue: f64) -> Hsp {
        Hsp {
            score,
            num_ident: score / 2,
            bit_score: score as f64 * 0.9,
            evalue,
            query_offset: 0,
            query_end: score,
            query_gapped_start: 0,
            subject_offset: 0,
            subject_end: score,
            subject_gapped_start: 0,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        }
    }

    static PIPE_RUN_COUNT: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
    static PIPE_FREE_COUNT: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
    static WRITER_INIT_COUNT: std::sync::atomic::AtomicUsize =
        std::sync::atomic::AtomicUsize::new(0);
    static WRITER_FINAL_COUNT: std::sync::atomic::AtomicUsize =
        std::sync::atomic::AtomicUsize::new(0);
    static WRITER_FREE_COUNT: std::sync::atomic::AtomicUsize =
        std::sync::atomic::AtomicUsize::new(0);

    fn counting_pipe_run(
        _data: Option<&mut BlastHSPPipeData>,
        results: &mut BlastHSPResults,
    ) -> i32 {
        PIPE_RUN_COUNT.fetch_add(1, std::sync::atomic::Ordering::SeqCst);
        if let Some(Some(hitlist)) = results.hitlist_array.get_mut(0) {
            if let Some(Some(hsp_list)) = hitlist.hsplist_array.get_mut(0) {
                hsp_list.oid = 99;
            }
        }
        0
    }

    fn counting_pipe_free(pipe: BlastHSPPipe) -> Option<BlastHSPPipe> {
        PIPE_FREE_COUNT.fetch_add(1, std::sync::atomic::Ordering::SeqCst);
        let _ = pipe;
        None
    }

    fn counting_writer_init(
        _data: Option<&mut BlastHSPWriterData>,
        _results: Option<&mut ()>,
    ) -> i32 {
        WRITER_INIT_COUNT.fetch_add(1, std::sync::atomic::Ordering::SeqCst);
        0
    }

    fn counting_writer_final(
        _data: Option<&mut BlastHSPWriterData>,
        _results: Option<&mut ()>,
    ) -> i32 {
        WRITER_FINAL_COUNT.fetch_add(1, std::sync::atomic::Ordering::SeqCst);
        0
    }

    fn counting_writer_free(writer: BlastHSPWriter) -> Option<BlastHSPWriter> {
        WRITER_FREE_COUNT.fetch_add(1, std::sync::atomic::Ordering::SeqCst);
        let _ = writer;
        None
    }

    #[test]
    fn test_hsp_stream_ordering() {
        // Add HSP lists for multiple subjects in arbitrary order,
        // then sort by evalue and verify OID ordering matches evalue ordering.
        let stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);

        // Add subjects out of OID order with varying evalues
        let oids_and_evalues = vec![(7, 1e-5), (2, 1e-20), (5, 1e-10), (0, 1e-2), (9, 1e-30)];

        for &(oid, evalue) in &oids_and_evalues {
            let mut list = HspList::new(oid);
            list.add_hsp(make_hsp(50, evalue));
            stream.blast_hspstream_write(0, list);
        }

        let mut results = stream.into_results();
        let hitlist = results.hitlists[0].as_mut().unwrap();

        // Sort by OID to verify we can retrieve them in OID order
        hitlist.hsp_lists.sort_by_key(|l| l.oid);
        let oids: Vec<i32> = hitlist.hsp_lists.iter().map(|l| l.oid).collect();
        assert_eq!(oids, vec![0, 2, 5, 7, 9]);

        // Also verify sort_by_evalue produces the correct evalue ordering
        hitlist.sort_by_evalue();
        for i in 1..hitlist.hsp_lists.len() {
            assert!(
                hitlist.hsp_lists[i].best_evalue >= hitlist.hsp_lists[i - 1].best_evalue,
                "HspLists should be sorted by evalue: {} >= {}",
                hitlist.hsp_lists[i].best_evalue,
                hitlist.hsp_lists[i - 1].best_evalue
            );
        }
    }

    #[test]
    fn translated_hsp_stream_lifecycle_wrappers_close_merge_and_trim() {
        let stream1 = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);
        let stream2 = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);

        let mut list1 = HspList::new(1);
        list1.add_hsp(make_hsp(20, 1e-3));
        let mut list2 = HspList::new(3);
        list2.add_hsp(make_hsp(30, 1e-4));
        assert_eq!(stream1.blast_hspstream_write(0, list1), 0);
        assert_eq!(stream2.blast_hspstream_write(0, list2), 0);
        assert_eq!(
            blast_hsp_stream_merge_simple(Some(&stream1), Some(&stream2)),
            0
        );
        blast_hsp_stream_close(Some(&stream1));
        assert_eq!(stream1.blast_hspstream_write(0, HspList::new(9)), -1);

        let mut results = stream1.into_results();
        let hitlist = results.hitlists[0].as_mut().expect("hitlist");
        assert_eq!(hitlist.hsp_lists.len(), 2);
        assert_eq!(hitlist.hsp_lists[0].oid, 3);
        assert!(s_sort_hsp_list_by_oid(&hitlist.hsp_lists[0], &hitlist.hsp_lists[1]) < 0);
        s_trim_hit_list(Some(hitlist), 1);
        assert_eq!(hitlist.hsp_lists.len(), 1);

        let mut slot = Some(blast_hsp_stream_new(
            crate::program::UNDEFINED,
            0,
            false,
            1,
            None,
        ));
        assert!(blast_hsp_stream_free(&mut slot).is_none());
        assert!(slot.is_none());
    }

    #[test]
    fn translated_hsp_stream_merge_combines_same_subject_lists() {
        let stream1 = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);
        let stream2 = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);

        let mut dst_list = HspList::new(7);
        dst_list.add_hsp(make_hsp(30, 1e-8));
        let mut src_list = HspList::new(7);
        src_list.add_hsp(make_hsp(40, 1e-12));
        src_list.add_hsp(make_hsp(20, 1e-5));

        assert_eq!(stream1.blast_hspstream_write(0, dst_list), 0);
        assert_eq!(stream2.blast_hspstream_write(0, src_list), 0);
        assert_eq!(
            blast_hsp_stream_merge_simple(Some(&stream1), Some(&stream2)),
            0
        );

        let results = stream1.into_results();
        let hitlist = results.hitlists[0].as_ref().expect("hitlist");
        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(hitlist.hsp_lists[0].oid, 7);
        let scores: Vec<i32> = hitlist.hsp_lists[0]
            .hsps
            .iter()
            .map(|hsp| hsp.score)
            .collect();
        assert_eq!(scores, vec![40, 30, 20]);
        assert_eq!(hitlist.hsp_lists[0].best_evalue, 1e-12);
    }

    #[test]
    fn translated_hsp_stream_merge_maps_split_query_chunk_to_global_contexts() {
        let chunk_stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);
        let global_stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 3, None);

        let mut local_list = HspList::new(42);
        let mut local_hsp = make_hsp(50, 1e-12);
        local_hsp.context = 0;
        local_hsp.query_offset = 0;
        local_hsp.query_end = 50;
        local_hsp.subject_offset = 100;
        local_hsp.subject_end = 150;
        local_list.hsps = vec![local_hsp];
        assert_eq!(chunk_stream.blast_hspstream_write(0, local_list), 0);

        let mut global_hitlist = blast_hit_list_new(10);
        let mut global_list = HspList::new(42);
        let mut global_hsp = make_hsp(100, 1e-20);
        global_hsp.context = 2;
        global_hsp.query_frame = 1;
        global_hsp.query_offset = 0;
        global_hsp.query_end = 110;
        global_hsp.subject_offset = 0;
        global_hsp.subject_end = 110;
        global_list.hsps = vec![global_hsp];
        global_hitlist.hsp_lists.push(global_list);
        global_stream.results.lock().unwrap().hitlists[2] = Some(global_hitlist);

        let mut split = crate::split_query::split_query_blk_new(1, true).expect("split query");
        assert_eq!(
            crate::split_query::split_query_blk_add_query_to_chunk(Some(&mut split), 2, 0),
            0
        );
        assert_eq!(
            crate::split_query::split_query_blk_add_context_to_chunk(Some(&mut split), 2, 0),
            0
        );
        assert_eq!(
            crate::split_query::split_query_blk_add_context_to_chunk(Some(&mut split), 3, 0),
            0
        );
        assert_eq!(
            crate::split_query::split_query_blk_add_context_offset_to_chunk(
                Some(&mut split),
                100,
                0,
            ),
            0
        );
        assert_eq!(
            crate::split_query::split_query_blk_add_context_offset_to_chunk(
                Some(&mut split),
                200,
                0,
            ),
            0
        );
        assert_eq!(
            crate::split_query::split_query_blk_set_chunk_overlap_size(Some(&mut split), 40),
            0
        );

        assert_eq!(
            blast_hsp_stream_merge(
                Some(&split),
                0,
                Some(&chunk_stream),
                Some(&global_stream),
                2
            ),
            0
        );
        let results = global_stream.into_results();
        let hitlist = results.hitlists[2].as_ref().expect("global query hitlist");
        assert_eq!(hitlist.hsp_lists.len(), 1);
        let hsp = &hitlist.hsp_lists[0].hsps[0];
        assert_eq!(hsp.context, 2);
        assert_eq!(hsp.query_frame, 1);
        assert_eq!(hsp.query_offset, 0);
        assert_eq!(hsp.query_end, 150);
        assert_eq!(hsp.subject_end, 150);
        assert_eq!(hsp.score, 140);
    }

    #[test]
    fn translated_hsp_stream_tback_close_runs_traceback_pipes_only() {
        PIPE_RUN_COUNT.store(0, std::sync::atomic::Ordering::SeqCst);
        PIPE_FREE_COUNT.store(0, std::sync::atomic::Ordering::SeqCst);

        let stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);
        let pipe = BlastHSPPipe {
            data: Some(BlastHSPPipeData::Opaque),
            run_fn_ptr: Some(counting_pipe_run),
            free_fn_ptr: Some(counting_pipe_free),
            next: None,
        };
        assert_eq!(
            blast_hsp_stream_register_pipe(Some(&stream), Some(pipe), E_TRACEBACK_SEARCH),
            0
        );
        let mut oid5 = HspList::new(5);
        oid5.add_hsp(make_hsp(20, 1e-3));
        let mut oid1 = HspList::new(1);
        oid1.add_hsp(make_hsp(30, 1e-4));

        assert_eq!(stream.blast_hspstream_write(0, oid5), 0);
        assert_eq!(stream.blast_hspstream_write(0, oid1), 0);
        let mut results = HspResults::new(1);
        let mut traceback_list = HspList::new(7);
        traceback_list.add_hsp(make_hsp(40, 1e-8));
        assert_eq!(
            blast_hsp_results_insert_hsp_list(&mut results, 0, traceback_list, 0),
            0
        );
        BlastHSPStreamTBackClose(Some(&stream), Some(&mut results));

        assert_eq!(blast_hsp_stream_pipe_counts(Some(&stream)), Some((0, 0)));
        assert_eq!(PIPE_RUN_COUNT.load(std::sync::atomic::Ordering::SeqCst), 1);
        assert_eq!(PIPE_FREE_COUNT.load(std::sync::atomic::Ordering::SeqCst), 1);
        let hitlist = results.hitlists[0].as_ref().expect("traceback hitlist");
        assert_eq!(hitlist.hsp_lists[0].oid, 99);
        assert_eq!(stream.blast_hspstream_write(0, HspList::new(9)), 0);
        blast_hsp_stream_close(Some(&stream));
        let (status, list) = blast_hsp_stream_read(Some(&stream));
        assert_eq!(status, K_BLAST_HSP_STREAM_SUCCESS);
        assert_eq!(list.expect("first list").oid, 1);
    }

    #[test]
    fn translated_hsp_stream_close_preserves_hsp_order_within_subject() {
        let stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);
        let mut list = HspList::new(11);
        list.add_hsp(make_hsp(10, 1e-2));
        list.add_hsp(make_hsp(50, 1e-5));
        list.add_hsp(make_hsp(20, 1e-3));

        assert_eq!(stream.blast_hspstream_write(0, list), 0);
        blast_hsp_stream_close(Some(&stream));

        let (status, list) = blast_hsp_stream_read(Some(&stream));
        assert_eq!(status, K_BLAST_HSP_STREAM_SUCCESS);
        let scores: Vec<i32> = list
            .expect("closed stream should return written list")
            .hsps
            .iter()
            .map(|hsp| hsp.score)
            .collect();
        assert_eq!(scores, vec![10, 50, 20]);
    }

    #[test]
    fn translated_hsp_stream_read_and_batch_read_follow_closed_oid_order() {
        let stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 2, None);
        let mut q0_oid4 = HspList::new(4);
        q0_oid4.add_hsp(make_hsp(20, 1e-3));
        let mut q0_oid2 = HspList::new(2);
        q0_oid2.add_hsp(make_hsp(30, 1e-4));
        let mut q1_oid2 = HspList::new(2);
        q1_oid2.add_hsp(make_hsp(40, 1e-5));

        assert_eq!(stream.blast_hspstream_write(0, q0_oid4), 0);
        assert_eq!(stream.blast_hspstream_write(0, q0_oid2), 0);
        assert_eq!(stream.blast_hspstream_write(1, q1_oid2), 0);

        let mut batch = blast_hsp_stream_result_batch_init(1);
        assert_eq!(
            blast_hsp_stream_batch_read(Some(&stream), Some(&mut batch)),
            K_BLAST_HSP_STREAM_SUCCESS
        );
        assert_eq!(batch.num_hsplists, 2);
        assert!(batch
            .hsplist_array
            .iter()
            .take(2)
            .all(|slot| slot.as_ref().is_some_and(|list| list.oid == 2)));

        let (status, next) = blast_hsp_stream_read(Some(&stream));
        assert_eq!(status, K_BLAST_HSP_STREAM_SUCCESS);
        assert_eq!(next.expect("hsp list").oid, 4);

        let (status, next) = blast_hsp_stream_read(Some(&stream));
        assert_eq!(status, K_BLAST_HSP_STREAM_EOF);
        assert!(next.is_none());
    }

    #[test]
    fn translated_hsp_stream_batch_read_resets_reused_batch() {
        let stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 2, None);
        let mut q0_oid4 = HspList::new(4);
        q0_oid4.add_hsp(make_hsp(20, 1e-3));
        let mut q0_oid2 = HspList::new(2);
        q0_oid2.add_hsp(make_hsp(30, 1e-4));
        let mut q1_oid2 = HspList::new(2);
        q1_oid2.add_hsp(make_hsp(40, 1e-5));

        assert_eq!(stream.blast_hspstream_write(0, q0_oid4), 0);
        assert_eq!(stream.blast_hspstream_write(0, q0_oid2), 0);
        assert_eq!(stream.blast_hspstream_write(1, q1_oid2), 0);

        let mut batch = blast_hsp_stream_result_batch_init(2);
        assert_eq!(
            blast_hsp_stream_batch_read(Some(&stream), Some(&mut batch)),
            K_BLAST_HSP_STREAM_SUCCESS
        );
        assert_eq!(batch.num_hsplists, 2);
        assert!(batch.hsplist_array[1].is_some());

        assert_eq!(
            blast_hsp_stream_batch_read(Some(&stream), Some(&mut batch)),
            K_BLAST_HSP_STREAM_SUCCESS
        );
        assert_eq!(batch.num_hsplists, 1);
        assert_eq!(batch.hsplist_array[0].as_ref().expect("next list").oid, 4);
        assert!(batch.hsplist_array[1].is_none());

        assert_eq!(
            blast_hsp_stream_batch_read(Some(&stream), Some(&mut batch)),
            K_BLAST_HSP_STREAM_EOF
        );
        assert_eq!(batch.num_hsplists, 0);
        assert!(batch.hsplist_array[0].is_none());
    }

    #[test]
    fn translated_hsp_stream_close_uses_cbs_score_sorted_results_path() {
        let stream = blast_hsp_stream_new(crate::program::BLASTP, 1, true, 1, None);
        for &(oid, evalue) in &[(1, 1e-10), (2, 1e-2), (3, 1e-20)] {
            let mut list = HspList::new(oid);
            list.add_hsp(make_hsp(50, evalue));
            assert_eq!(stream.blast_hspstream_write(0, list), 0);
        }

        blast_hsp_stream_close(Some(&stream));

        assert!(stream.sorted_hsplists.lock().unwrap().is_empty());
        let hitlist_oids: Vec<i32> = stream.results.lock().unwrap().hitlists[0]
            .as_ref()
            .expect("hitlist")
            .hsp_lists
            .iter()
            .map(|list| list.oid)
            .collect();
        assert_eq!(hitlist_oids, vec![2, 1, 3]);

        let (status, first) = blast_hsp_stream_read(Some(&stream));
        assert_eq!(status, K_BLAST_HSP_STREAM_SUCCESS);
        assert_eq!(first.expect("first list").oid, 3);
    }

    #[test]
    fn translated_hsp_stream_writer_lifecycle_runs_once_and_frees() {
        WRITER_INIT_COUNT.store(0, std::sync::atomic::Ordering::SeqCst);
        WRITER_FINAL_COUNT.store(0, std::sync::atomic::Ordering::SeqCst);
        WRITER_FREE_COUNT.store(0, std::sync::atomic::Ordering::SeqCst);

        let writer = BlastHSPWriter {
            data: Some(BlastHSPWriterData::Opaque),
            init_fn_ptr: Some(counting_writer_init),
            run_fn_ptr: None,
            final_fn_ptr: Some(counting_writer_final),
            free_fn_ptr: Some(counting_writer_free),
        };
        let mut stream = Some(blast_hsp_stream_new(
            crate::program::BLASTN,
            0,
            false,
            1,
            Some(writer),
        ));
        let stream_ref = stream.as_ref().expect("stream");
        let mut list = HspList::new(8);
        list.add_hsp(make_hsp(40, 1e-8));
        assert_eq!(stream_ref.blast_hspstream_write(0, list), 0);

        blast_hsp_stream_close(Some(stream_ref));
        blast_hsp_stream_close(Some(stream_ref));
        assert_eq!(
            WRITER_INIT_COUNT.load(std::sync::atomic::Ordering::SeqCst),
            1
        );
        assert_eq!(
            WRITER_FINAL_COUNT.load(std::sync::atomic::Ordering::SeqCst),
            1
        );

        assert!(blast_hsp_stream_free(&mut stream).is_none());
        assert_eq!(
            WRITER_FREE_COUNT.load(std::sync::atomic::Ordering::SeqCst),
            1
        );
    }

    #[test]
    fn translated_hsp_stream_to_batches_groups_by_subject_oid() {
        let stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 2, None);
        let mut q0_oid4 = HspList::new(4);
        q0_oid4.add_hsp(make_hsp(20, 1e-3));
        let mut q0_oid2 = HspList::new(2);
        q0_oid2.add_hsp(make_hsp(30, 1e-4));
        let mut q1_oid2 = HspList::new(2);
        q1_oid2.add_hsp(make_hsp(40, 1e-5));

        assert_eq!(stream.blast_hspstream_write(0, q0_oid4), 0);
        assert_eq!(stream.blast_hspstream_write(0, q0_oid2), 0);
        assert_eq!(stream.blast_hspstream_write(1, q1_oid2), 0);
        assert_eq!(s_blast_hsp_stream_count_num_oids(Some(&stream)), 2);

        let (status, batches) = blast_hsp_stream_to_hsp_stream_results_batch(Some(&stream));
        assert_eq!(status, K_BLAST_HSP_STREAM_SUCCESS);
        let batches = batches.expect("batches");
        assert_eq!(batches.num_batches, 2);

        let first = batches.array_of_batches[0].as_ref().expect("first batch");
        let second = batches.array_of_batches[1].as_ref().expect("second batch");
        assert_eq!(first.num_hsplists, 2);
        assert_eq!(second.num_hsplists, 1);
        assert!(first
            .hsplist_array
            .iter()
            .take(2)
            .all(|slot| slot.as_ref().is_some_and(|list| list.oid == 2)));
        assert_eq!(second.hsplist_array[0].as_ref().unwrap().oid, 4);
    }

    #[test]
    fn translated_hsp_cbs_stream_close_trims_large_weak_tail() {
        let stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);
        let mut hitlist = HitList::new();
        for index in 0..1305 {
            let mut list = HspList::new(index);
            let evalue = if index <= 500 {
                1e-50
            } else {
                1e-5 + index as f64 * 1e-8
            };
            list.add_hsp(make_hsp(100 - (index as i32 % 10), evalue));
            list.best_evalue = s_blast_get_best_evalue(&list);
            hitlist.hsp_lists.push(list);
        }
        {
            let mut results = stream.results.lock().unwrap();
            results.hitlists[0] = Some(hitlist);
        }

        blast_hsp_cbs_stream_close(Some(&stream), 500);
        assert!(*stream.closed.lock().unwrap());
        let results = stream.into_results();
        let hitlist = results.hitlists[0].as_ref().expect("hitlist");
        assert!(hitlist.hsp_lists.len() < 1305);
    }

    #[test]
    fn translated_hsp_stream_pipe_and_writer_wrappers_match_status_shape() {
        let stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);
        assert_eq!(blast_hsp_stream_register_mt_lock(Some(&stream), true), 0);
        assert_eq!(blast_hsp_stream_register_mt_lock(None, true), -1);
        assert_eq!(
            blast_hsp_stream_register_pipe(
                Some(&stream),
                Some(blast_hsp_pipe_new()),
                E_PRELIM_SEARCH,
            ),
            0
        );
        assert_eq!(
            blast_hsp_stream_register_pipe(
                Some(&stream),
                Some(blast_hsp_pipe_new()),
                E_TRACEBACK_SEARCH,
            ),
            0
        );
        assert_eq!(blast_hsp_stream_pipe_counts(Some(&stream)), Some((1, 1)));
        assert_eq!(
            blast_hsp_stream_register_pipe(Some(&stream), Some(blast_hsp_pipe_new()), E_BOTH),
            -1
        );
        assert_eq!(
            blast_hsp_stream_register_pipe(Some(&stream), None, E_BOTH),
            -1
        );
        let writer = blast_hsp_writer_new();
        assert!(writer.data.is_none());
        assert!(writer.init_fn_ptr.is_none());
        assert!(writer.run_fn_ptr.is_none());
        assert!(writer.final_fn_ptr.is_none());
        assert!(writer.free_fn_ptr.is_none());
        blast_hsp_stream_simple_close(Some(&stream));
        assert_eq!(blast_hsp_stream_pipe_counts(Some(&stream)), Some((0, 1)));
        blast_hsp_stream_mapping_close(Some(&stream));
        let mut results = HspResults::new(1);
        BlastHSPStreamTBackClose(Some(&stream), Some(&mut results));
        assert_eq!(blast_hsp_stream_pipe_counts(Some(&stream)), Some((0, 0)));
    }

    #[test]
    fn translated_hsp_stream_close_runs_and_frees_prelim_pipes() {
        PIPE_RUN_COUNT.store(0, std::sync::atomic::Ordering::SeqCst);
        PIPE_FREE_COUNT.store(0, std::sync::atomic::Ordering::SeqCst);

        let stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);
        let mut list = HspList::new(7);
        list.add_hsp(make_hsp(40, 1e-8));
        assert_eq!(stream.blast_hspstream_write(0, list), 0);

        let pipe = BlastHSPPipe {
            data: Some(BlastHSPPipeData::Opaque),
            run_fn_ptr: Some(counting_pipe_run),
            free_fn_ptr: Some(counting_pipe_free),
            next: None,
        };
        assert_eq!(
            blast_hsp_stream_register_pipe(Some(&stream), Some(pipe), E_PRELIM_SEARCH),
            0
        );

        blast_hsp_stream_close(Some(&stream));

        assert_eq!(PIPE_RUN_COUNT.load(std::sync::atomic::Ordering::SeqCst), 1);
        assert_eq!(PIPE_FREE_COUNT.load(std::sync::atomic::Ordering::SeqCst), 1);
        assert_eq!(blast_hsp_stream_pipe_counts(Some(&stream)), Some((0, 0)));
        let (_, read_list) = blast_hsp_stream_read(Some(&stream));
        assert_eq!(read_list.unwrap().oid, 99);
    }

    #[test]
    fn hsp_and_hit_list_merge_append_non_split_lists() {
        let mut incoming = Some(HspList::new(7));
        incoming.as_mut().unwrap().add_hsp(make_hsp(20, 1e-5));
        let mut combined = Some(HspList::new(7));
        combined.as_mut().unwrap().add_hsp(make_hsp(30, 1e-8));

        assert_eq!(
            blast_hsp_lists_merge_simple(&mut incoming, &mut combined, 10),
            0
        );
        assert!(incoming.is_none());
        let merged = combined.as_ref().unwrap();
        assert_eq!(merged.hsps.len(), 2);
        assert_eq!(merged.best_evalue, 1e-8);

        let mut incoming = Some(HspList::new(7));
        incoming.as_mut().unwrap().add_hsp(make_hsp(50, 1e-20));
        incoming.as_mut().unwrap().add_hsp(make_hsp(15, 1e-3));
        let mut combined = Some(HspList::new(7));
        combined.as_mut().unwrap().add_hsp(make_hsp(40, 1e-12));
        combined.as_mut().unwrap().add_hsp(make_hsp(30, 1e-8));
        assert_eq!(
            blast_hsp_lists_merge_simple(&mut incoming, &mut combined, 2),
            0
        );
        assert!(incoming.is_none());
        let scores: Vec<i32> = combined
            .as_ref()
            .unwrap()
            .hsps
            .iter()
            .map(|hsp| hsp.score)
            .collect();
        assert_eq!(scores, vec![50, 40]);

        let mut empty_incoming = Some(HspList::new(7));
        let mut combined = Some(HspList::new(7));
        combined.as_mut().unwrap().add_hsp(make_hsp(30, 1e-8));
        assert_eq!(
            blast_hsp_lists_merge_simple(&mut empty_incoming, &mut combined, 10),
            0
        );
        assert!(empty_incoming.is_some());
        assert!(empty_incoming.as_ref().unwrap().hsps.is_empty());
        assert_eq!(combined.as_ref().unwrap().hsps.len(), 1);

        let mut incoming = Some(blast_hsp_list_new(5));
        incoming.as_mut().unwrap().oid = 7;
        incoming.as_mut().unwrap().add_hsp(make_hsp(50, 1e-20));
        let mut combined = Some(blast_hsp_list_new(5));
        combined.as_mut().unwrap().oid = 7;
        combined.as_mut().unwrap().add_hsp(make_hsp(40, 1e-12));
        assert_eq!(blast_hsp_list_append(&mut incoming, &mut combined, 1), 0);
        assert_eq!(combined.as_ref().unwrap().hsps.len(), 1);
        assert_eq!(combined.as_ref().unwrap().hsp_max, 5);

        let mut incoming = Some(HspList::new(9));
        let mut overlap_new = make_hsp(50, 1e-12);
        overlap_new.query_offset = 90;
        overlap_new.query_end = 150;
        overlap_new.subject_offset = 90;
        overlap_new.subject_end = 150;
        overlap_new.query_gapped_start = 90;
        overlap_new.subject_gapped_start = 90;
        let mut far_new = make_hsp(30, 1e-3);
        far_new.query_offset = 300;
        far_new.query_end = 340;
        far_new.subject_offset = 300;
        far_new.subject_end = 340;
        incoming.as_mut().unwrap().hsps = vec![overlap_new, far_new];

        let mut combined = Some(HspList::new(9));
        let mut overlap_old = make_hsp(100, 1e-20);
        overlap_old.query_offset = 0;
        overlap_old.query_end = 100;
        overlap_old.subject_offset = 0;
        overlap_old.subject_end = 100;
        combined.as_mut().unwrap().hsps = vec![overlap_old];

        assert_eq!(
            blast_hsp_lists_merge(
                &mut incoming,
                &mut combined,
                10,
                Some(&[80]),
                i32::MIN,
                40,
                true,
                false,
            ),
            0
        );
        let merged = combined.as_ref().unwrap();
        assert_eq!(merged.hsps.len(), 2);
        assert_eq!(merged.hsps[0].query_offset, 0);
        assert_eq!(merged.hsps[0].query_end, 150);
        assert_eq!(merged.hsps[0].subject_offset, 0);
        assert_eq!(merged.hsps[0].subject_end, 150);
        assert_eq!(merged.hsps[0].score, 140);

        let mut old_hit = Some(blast_hit_list_new(10));
        let mut oid7 = HspList::new(7);
        oid7.add_hsp(make_hsp(40, 1e-12));
        let mut oid3 = HspList::new(3);
        oid3.add_hsp(make_hsp(10, 1e-2));
        old_hit.as_mut().unwrap().blast_hit_list_update(oid7);
        old_hit.as_mut().unwrap().blast_hit_list_update(oid3);

        let mut combined_hit = Some(blast_hit_list_new(10));
        let mut existing = HspList::new(7);
        existing.add_hsp(make_hsp(25, 1e-6));
        combined_hit
            .as_mut()
            .unwrap()
            .blast_hit_list_update(existing);

        assert_eq!(
            blast_hit_list_merge_simple(&mut old_hit, &mut combined_hit),
            0
        );
        assert!(old_hit.is_none());
        let hit = combined_hit.as_ref().unwrap();
        // NCBI `Blast_HitListMerge` (`blast_hits.c:2216`) returns the merged
        // hitlist in oid-ascending order (the merge walks both inputs
        // sorted by oid). No final evalue sort here.
        assert_eq!(
            hit.hsp_lists.iter().map(|l| l.oid).collect::<Vec<_>>(),
            vec![3, 7]
        );
        assert_eq!(
            hit.hsp_lists
                .iter()
                .find(|l| l.oid == 7)
                .unwrap()
                .hsps
                .len(),
            2
        );
    }

    #[test]
    fn hit_list_merge_uses_append_when_split_offsets_are_not_active_like_c() {
        let mut old_hit = Some(blast_hit_list_new(10));
        let mut old_list = HspList::new(7);
        let mut old_hsp = make_hsp(50, 1e-20);
        old_hsp.query_offset = 0;
        old_hsp.query_end = 50;
        old_hsp.subject_offset = 0;
        old_hsp.subject_end = 50;
        old_list.add_hsp(old_hsp);
        old_hit.as_mut().unwrap().blast_hit_list_update(old_list);

        let mut combined_hit = Some(blast_hit_list_new(10));
        let mut combined_list = HspList::new(7);
        let mut combined_hsp = make_hsp(40, 1e-10);
        combined_hsp.query_offset = 0;
        combined_hsp.query_end = 40;
        combined_hsp.subject_offset = 0;
        combined_hsp.subject_end = 40;
        combined_list.add_hsp(combined_hsp);
        combined_hit
            .as_mut()
            .unwrap()
            .blast_hit_list_update(combined_list);

        assert_eq!(
            blast_hit_list_merge(&mut old_hit, &mut combined_hit, 1, Some(&[0]), 40, true,),
            0
        );

        let hit = combined_hit.as_ref().expect("merged hitlist");
        let list = hit
            .hsp_lists
            .iter()
            .find(|list| list.oid == 7)
            .expect("same oid list");
        assert_eq!(list.hsps.len(), 2);
        assert_eq!(
            list.hsps.iter().map(|hsp| hsp.score).collect::<Vec<_>>(),
            vec![50, 40]
        );
    }

    #[test]
    fn hit_list_merge_uses_old_hitlist_capacity_like_c() {
        let mut old_hit = Some(blast_hit_list_new(1));
        let mut old_list = HspList::new(1);
        old_list.add_hsp(make_hsp(30, 1e-5));
        old_hit.as_mut().unwrap().blast_hit_list_update(old_list);

        let mut combined_hit = Some(blast_hit_list_new(10));
        for (oid, score, evalue) in [(2, 40, 1e-10), (3, 50, 1e-20)] {
            let mut list = HspList::new(oid);
            list.add_hsp(make_hsp(score, evalue));
            combined_hit.as_mut().unwrap().blast_hit_list_update(list);
        }

        assert_eq!(
            blast_hit_list_merge_simple(&mut old_hit, &mut combined_hit),
            0
        );

        let hit = combined_hit.as_ref().expect("merged hitlist");
        assert_eq!(hit.hsplist_max, 1);
        assert_eq!(
            hit.hsp_lists
                .iter()
                .map(|list| list.oid)
                .collect::<Vec<_>>(),
            vec![3]
        );
    }

    #[test]
    fn test_hsp_stream_empty() {
        let stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 3, None);
        let results = stream.into_results();

        // All query slots should be None
        assert_eq!(results.hitlists.len(), 3);
        for hitlist in &results.hitlists {
            assert!(hitlist.is_none(), "Empty stream should have no results");
        }
    }

    #[test]
    fn test_hsp_stream_single_subject() {
        let stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);

        let mut list = HspList::new(42);
        list.add_hsp(make_hsp(75, 1e-15));
        list.add_hsp(make_hsp(30, 1e-3));
        stream.blast_hspstream_write(0, list);

        let results = stream.into_results();

        assert!(results.hitlists[0].is_some());
        let hitlist = results.hitlists[0].as_ref().unwrap();
        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(hitlist.hsp_lists[0].oid, 42);
        assert_eq!(hitlist.hsp_lists[0].hsps.len(), 2);
        // best_evalue should track the minimum
        assert_eq!(hitlist.hsp_lists[0].best_evalue, 1e-15);
    }

    #[test]
    fn test_blast_hsp_list_new_and_free() {
        let list = blast_hsp_list_new(5);
        assert_eq!(list.oid, -1);
        assert_eq!(list.hsp_max, 5);
        assert!(list.hsps.is_empty());
        assert!(list.hsps.capacity() <= 5);

        let unlimited = blast_hsp_list_new(0);
        assert_eq!(unlimited.hsp_max, i32::MAX);
        assert_eq!(unlimited.hsps.capacity(), 100);

        let mut slot = Some(list);
        assert!(blast_hsp_list_free(&mut slot).is_none());
        assert!(slot.is_none());
    }

    #[test]
    fn test_blast_hsp_lifecycle_and_mapping_info_helpers() {
        let mut hsp = make_hsp(30, 1e-4);
        hsp.edit_script = crate::gapinfo::gap_edit_script_new(2);
        if let Some(script) = hsp.edit_script.as_mut() {
            script.push(crate::gapinfo::GapAlignOpType::Sub, 10);
        }
        let mut slot = Some(hsp);
        assert!(blast_hsp_free(&mut slot).is_none());
        assert!(slot.is_none());
        let mut empty_slot = None;
        assert!(blast_hsp_free(&mut empty_slot).is_none());

        let mut info = blast_hsp_mapping_info_new();
        assert_eq!(info.left_edge, 0);
        assert_eq!(info.right_edge, 0);
        info.edits = crate::gapinfo::jumper_edits_block_new(1);
        info.subject_overhangs = Some(crate::gapinfo::SequenceOverhangs {
            left_len: 0,
            right_len: 0,
            left: Some(vec![1, 2]),
            right: Some(vec![3, 4]),
        });
        assert!(blast_hsp_mapping_info_free(Some(info)).is_none());

        assert_eq!(blast_hsp_num_max(false, 7), 7);
        assert_eq!(blast_hsp_num_max(true, 7), 7);
        assert_eq!(blast_hsp_num_max(false, 0), i32::MAX);
        assert_eq!(blast_hsp_num_max(true, 0), i32::MAX);
    }

    #[test]
    fn test_hit_parameter_helpers_match_c_sizing_rules() {
        let _guard = ENV_LOCK.lock().unwrap();
        let old_adaptive = std::env::var_os("ADAPTIVE_CBS");
        std::env::remove_var("ADAPTIVE_CBS");

        assert_eq!(get_prelim_hitlist_size(4, 0, true), 10);
        assert_eq!(get_prelim_hitlist_size(50, 0, true), 100);
        assert_eq!(get_prelim_hitlist_size(500, 1, true), 1050);
        assert_eq!(get_prelim_hitlist_size(501, 1, true), 1052);

        std::env::set_var("ADAPTIVE_CBS", "1");
        assert_eq!(get_prelim_hitlist_size(500, 1, true), 1500);
        assert_eq!(get_prelim_hitlist_size(1000, 1, true), 2050);

        if let Some(value) = old_adaptive {
            std::env::set_var("ADAPTIVE_CBS", value);
        } else {
            std::env::remove_var("ADAPTIVE_CBS");
        }

        let hit = HitSavingOptions {
            hitlist_size: 25,
            ..HitSavingOptions::default()
        };
        let scoring = ScoringOptions::new_blastn();
        let mut params = None;
        assert_eq!(
            sblast_hits_parameters_new(Some(&hit), 0, Some(&scoring), 7, &mut params),
            0
        );
        assert_eq!(
            params,
            Some(SBlastHitsParameters {
                prelim_hitlist_size: 50,
                hsp_num_max: 7,
            })
        );
        assert!(sblast_hits_parameters_free(params).is_none());

        let mut missing = Some(SBlastHitsParameters::default());
        assert_eq!(
            sblast_hits_parameters_new(None, 0, Some(&scoring), 7, &mut missing),
            1
        );
        assert!(missing.is_none());
    }

    #[test]
    fn test_hsp_identity_length_and_query_end_helpers_match_c_ordering() {
        let mut hsp = make_hsp(40, 1e-8);
        hsp.num_ident = 9;
        hsp.query_offset = 5;
        hsp.query_end = 20;
        hsp.subject_offset = 2;
        hsp.subject_end = 19;

        let hit_options = HitSavingOptions {
            percent_identity: 50.0,
            min_hit_length: 10,
            ..HitSavingOptions::default()
        };
        let blast_hsp = BlastHSP::from_legacy_hsp(hsp.clone());
        assert!(!s_hsp_test(&blast_hsp, &hit_options, 18));
        assert!(!blast_hsp_test(&hsp, &hit_options, 18));
        assert!(s_hsp_test(&blast_hsp, &hit_options, 9));

        let strict_identity = HitSavingOptions {
            percent_identity: 60.0,
            min_hit_length: 0,
            ..HitSavingOptions::default()
        };
        assert!(blast_hsp_test(&hsp, &strict_identity, 18));
        assert_eq!(s_hsp_start_diag(&hsp), 3);
        assert_eq!(s_hsp_end_diag(&hsp), 1);
        assert_eq!(blast_hsp_get_query_coverage(&hsp, 30), 50.5);
        assert!(blast_hsp_query_coverage_test(&hsp, 60.0, 30));
        assert!(!blast_hsp_query_coverage_test(&hsp, 50.0, 30));
        assert_eq!(blast_hsp_get_query_coverage(&hsp, 0), 0.0);

        let query = [0, 1, 2, 3, 0, 1];
        let subject = [0, 1, 3, 3, 0, 2];
        hsp.query_offset = 0;
        hsp.query_end = 6;
        hsp.subject_offset = 0;
        hsp.subject_end = 6;
        let (status, align_length) = blast_hsp_get_num_identities_plain(&query, &subject, &mut hsp);
        assert_eq!(status, 0);
        assert_eq!(align_length, 6);
        assert_eq!(hsp.num_ident, 4);

        let mut oof = make_hsp(40, 1e-8);
        oof.query_offset = 1;
        oof.subject_offset = 2;
        oof.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![
            (GapAlignOpType::Sub, 2),
            (GapAlignOpType::Ins, 1),
            (GapAlignOpType::Del, 1),
            (GapAlignOpType::Del2, 1),
            (GapAlignOpType::Ins2, 1),
            (GapAlignOpType::Sub, 1),
        ]));
        let oof_query = [0, 7, 8, 99, 10];
        let oof_subject = [0, 0, 7, 0, 0, 9, 0, 0, 0, 0, 0, 10];
        let mut matrix = vec![vec![0; 256]; 256];
        matrix[8][9] = 2;
        let (status, num_ident, align_length, num_pos) =
            s_blast_hsp_get_oof_num_identities_and_positives(
                &oof_query,
                &oof_subject,
                &oof,
                crate::program::TBLASTN,
                Some(&matrix),
            );
        assert_eq!(status, 0);
        assert_eq!(num_ident, 2);
        assert_eq!(align_length, 5);
        assert_eq!(num_pos, 3);

        let mut subject_blk = crate::util::BlastSequenceBlk::default();
        subject_blk.sequence = Some(vec![0, 1, 2, 4, 8, 1, 2, 4, 8, 1, 2, 4, 8, 1, 2, 4]);
        subject_blk.length = subject_blk.sequence.as_ref().unwrap().len() as i32;
        let mut target = crate::util::SBlastTargetTranslation {
            program_number: crate::program::TBLASTN,
            gen_code_string: crate::util::STANDARD_GENETIC_CODE,
            translations: vec![None; crate::util::NUM_FRAMES],
            partial: true,
            num_frames: crate::util::NUM_FRAMES as i32,
            range: Some(vec![0; 2 * crate::util::NUM_FRAMES]),
            subject_blk: Some(subject_blk),
        };
        let mut translated_hsp = BlastHSP::from_legacy_hsp(make_hsp(40, 1e-8));
        translated_hsp.subject.frame = 1;
        translated_hsp.subject.offset = 1;
        translated_hsp.subject.end = 3;
        let view = blast_hsp_get_target_translation(&mut target, Some(&translated_hsp))
            .expect("translation view");
        assert_eq!(view.pointer_offset, 1);
        assert_eq!(view.translated_length, 5);
        assert_eq!(view.sequence[0], crate::util::FENCE_SENTRY);
        assert_eq!(view.sequence[6], crate::util::FENCE_SENTRY);
        assert_eq!(view.get(0), Some(view.sequence[1]));

        let mut earlier_context = hsp.clone();
        earlier_context.context = -1;
        assert!(s_query_end_compare_hsps(Some(&earlier_context), Some(&hsp)).is_lt());
        assert!(s_query_end_compare_hsps(None, Some(&hsp)).is_gt());
        assert!(s_query_end_compare_hsps(Some(&hsp), None).is_lt());
        assert_eq!(
            s_query_end_compare_hsps(None, None),
            std::cmp::Ordering::Equal
        );

        let mut later_query = hsp.clone();
        later_query.query_end += 1;
        assert!(s_query_end_compare_hsps(Some(&hsp), Some(&later_query)).is_lt());

        let mut higher_score = hsp.clone();
        higher_score.score += 1;
        assert!(s_query_end_compare_hsps(Some(&higher_score), Some(&hsp)).is_lt());

        let mut shorter_query_range = hsp.clone();
        shorter_query_range.query_offset += 1;
        assert!(s_query_end_compare_hsps(Some(&shorter_query_range), Some(&hsp)).is_lt());
    }

    #[test]
    fn target_translation_view_handles_pretranslated_and_invalid_inputs() {
        let mut target = crate::util::SBlastTargetTranslation {
            program_number: crate::program::TBLASTN,
            gen_code_string: crate::util::STANDARD_GENETIC_CODE,
            translations: vec![None; crate::util::NUM_FRAMES],
            partial: false,
            num_frames: crate::util::NUM_FRAMES as i32,
            range: Some(vec![0; 2 * crate::util::NUM_FRAMES]),
            subject_blk: None,
        };
        target.translations[0] = Some(vec![10, 11, 12, 13, 14, 15]);
        target.range.as_mut().unwrap()[0] = 2;
        target.range.as_mut().unwrap()[1] = 5;
        let mut hsp = BlastHSP::from_legacy_hsp(make_hsp(30, 1.0e-6));
        hsp.subject.frame = 1;

        let view = blast_hsp_get_target_translation(&mut target, Some(&hsp)).expect("view");
        assert_eq!(view.pointer_offset, -1);
        assert_eq!(view.translated_length, 5);
        assert_eq!(view.get(2), Some(11));
        assert_eq!(view.get(3), Some(12));

        assert!(blast_hsp_get_target_translation(&mut target, None).is_none());

        let mut no_range = target.clone();
        no_range.range = None;
        assert!(blast_hsp_get_target_translation(&mut no_range, Some(&hsp)).is_none());

        let mut no_translation = target;
        no_translation.translations[0] = None;
        assert!(blast_hsp_get_target_translation(&mut no_translation, Some(&hsp)).is_none());

        let mut partial_without_subject = crate::util::SBlastTargetTranslation {
            program_number: crate::program::TBLASTN,
            gen_code_string: crate::util::STANDARD_GENETIC_CODE,
            translations: vec![None; crate::util::NUM_FRAMES],
            partial: true,
            num_frames: crate::util::NUM_FRAMES as i32,
            range: Some(vec![0; 2 * crate::util::NUM_FRAMES]),
            subject_blk: None,
        };
        assert!(
            blast_hsp_get_target_translation(&mut partial_without_subject, Some(&hsp)).is_none()
        );
    }

    #[test]
    fn target_translation_uses_program_aware_context_for_nucleotide_minus_strand() {
        let mut target = crate::util::SBlastTargetTranslation {
            program_number: crate::program::BLASTN,
            gen_code_string: crate::util::STANDARD_GENETIC_CODE,
            translations: vec![None; crate::util::NUM_FRAMES],
            partial: false,
            num_frames: crate::util::NUM_FRAMES as i32,
            range: Some(vec![0; 2 * crate::util::NUM_FRAMES]),
            subject_blk: None,
        };
        target.translations[1] = Some(vec![20, 21, 22, 23]);
        target.range.as_mut().unwrap()[2] = 1;
        target.range.as_mut().unwrap()[3] = 3;

        let mut hsp = BlastHSP::from_legacy_hsp(make_hsp(30, 1.0e-6));
        hsp.subject.frame = -1;

        let view = blast_hsp_get_target_translation(&mut target, Some(&hsp)).expect("view");
        assert_eq!(view.pointer_offset, 0);
        assert_eq!(view.translated_length, 3);
        assert_eq!(view.get(1), Some(21));
    }

    #[test]
    fn target_translation_reallocates_against_old_cached_width_like_c() {
        let sequence = (0..600).map(|i| [1u8, 2, 4, 8][i % 4]).collect::<Vec<_>>();
        let mut subject_blk = crate::util::BlastSequenceBlk::default();
        subject_blk.sequence = Some(sequence);
        subject_blk.length = 600;

        let mut target = crate::util::SBlastTargetTranslation {
            program_number: crate::program::TBLASTN,
            gen_code_string: crate::util::STANDARD_GENETIC_CODE,
            translations: vec![None; crate::util::NUM_FRAMES],
            partial: true,
            num_frames: crate::util::NUM_FRAMES as i32,
            range: Some(vec![0; 2 * crate::util::NUM_FRAMES]),
            subject_blk: Some(subject_blk),
        };
        target.range.as_mut().unwrap()[0] = 100;
        target.range.as_mut().unwrap()[1] = 200;
        target.translations[0] = Some(vec![0; 102]);

        let mut hsp = BlastHSP::from_legacy_hsp(make_hsp(30, 1.0e-6));
        hsp.subject.frame = 1;
        hsp.subject.offset = 20;
        hsp.subject.end = 80;

        let view = blast_hsp_get_target_translation(&mut target, Some(&hsp)).expect("view");
        assert_eq!(view.pointer_offset, 1);
        assert!(view.sequence.len() >= 116);
        assert_eq!(target.range.as_ref().unwrap()[0], 0);
        assert_eq!(target.range.as_ref().unwrap()[1], 113);
    }

    #[test]
    fn hsp_adjusted_offsets_handle_translated_frames() {
        let mut hsp = make_hsp(10, 1.0e-6);
        hsp.query_offset = 1;
        hsp.query_end = 4;
        hsp.query_frame = -2;
        hsp.subject_offset = 5;
        hsp.subject_end = 11;

        assert_eq!(
            blast_hsp_get_adjusted_offsets(crate::program::BLASTX, &hsp, 99, 100),
            (2, 4, 6, 11)
        );

        hsp.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![(
            GapAlignOpType::Sub,
            3,
        )]));
        assert_eq!(
            blast_hsp_get_adjusted_offsets(crate::program::BLASTX, &hsp, 99, 100),
            (94, 86, 6, 11)
        );

        hsp.query_frame = 0;
        hsp.subject_offset = 2;
        hsp.subject_end = 5;
        hsp.subject_frame = 2;
        assert_eq!(
            blast_hsp_get_adjusted_offsets(crate::program::TBLASTN, &hsp, 20, 120),
            (2, 4, 2, 15)
        );

        hsp.query_offset = 10;
        hsp.query_end = 20;
        hsp.query_frame = 1;
        hsp.subject_offset = 5;
        hsp.subject_end = 15;
        hsp.subject_frame = -1;
        assert_eq!(
            blast_hsp_get_adjusted_offsets(crate::program::BLASTN, &hsp, 100, 100),
            (81, 90, 15, 6)
        );
    }

    #[test]
    fn partial_subject_translation_shifts_hsp_and_returns_window() {
        let mut seq = Vec::new();
        for index in 0..240 {
            seq.push([1u8, 2, 4, 8][index % 4]);
        }
        let subject_blk = crate::util::BlastSequenceBlk {
            sequence: Some(seq),
            length: 240,
            ..Default::default()
        };
        let mut hsp = make_hsp(10, 1.0e-6);
        hsp.subject_frame = 1;
        hsp.subject_offset = 40;
        hsp.subject_end = 45;
        hsp.subject_gapped_start = 42;

        let partial = blast_hsp_get_partial_subject_translation(
            &subject_blk,
            &mut hsp,
            false,
            &crate::util::STANDARD_GENETIC_CODE,
        )
        .expect("partial translation");

        assert_eq!(partial.start_shift, 0);
        assert_eq!(
            (
                hsp.subject_offset,
                hsp.subject_end,
                hsp.subject_gapped_start
            ),
            (40, 45, 42)
        );
        assert_eq!(partial.subject_length, partial.subject.len() as i32);
        assert_eq!(
            partial.subject,
            partial.translation_buffer[1..1 + partial.subject_length as usize]
        );
    }

    #[test]
    fn query_offset_compare_hsps_matches_c_tie_breakers() {
        let hsp = make_hsp(30, 1.0e-6);
        let mut earlier_context = hsp.clone();
        earlier_context.context = -1;
        assert!(s_query_offset_compare_hsps(Some(&earlier_context), Some(&hsp)).is_lt());
        assert!(s_query_offset_compare_hsps(None, Some(&hsp)).is_gt());
        assert!(s_query_offset_compare_hsps(Some(&hsp), None).is_lt());
        assert_eq!(
            s_query_offset_compare_hsps(None, None),
            std::cmp::Ordering::Equal
        );

        let mut later_query = hsp.clone();
        later_query.query_offset += 1;
        assert!(s_query_offset_compare_hsps(Some(&hsp), Some(&later_query)).is_lt());

        let mut higher_score = hsp.clone();
        higher_score.score += 1;
        assert!(s_query_offset_compare_hsps(Some(&higher_score), Some(&hsp)).is_lt());

        let mut longer_query_range = hsp.clone();
        longer_query_range.query_end += 1;
        assert!(s_query_offset_compare_hsps(Some(&longer_query_range), Some(&hsp)).is_lt());
    }

    #[test]
    fn test_cut_off_gap_edit_script_matches_c_begin_and_end_cuts() {
        let mut begin = make_hsp(50, 1e-8);
        begin.query_offset = 10;
        begin.subject_offset = 100;
        begin.query_end = 20;
        begin.subject_end = 112;
        begin.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![
            (GapAlignOpType::Sub, 3),
            (GapAlignOpType::Del, 2),
            (GapAlignOpType::Sub, 5),
            (GapAlignOpType::Ins, 1),
        ]));

        s_cut_off_gap_edit_script(&mut begin, 15, 105, true);
        assert_eq!((begin.query_offset, begin.subject_offset), (15, 107));
        assert_eq!(
            begin.edit_script.as_ref().unwrap().ops_vec(),
            vec![(GapAlignOpType::Sub, 3), (GapAlignOpType::Ins, 1)]
        );

        let mut end = make_hsp(50, 1e-8);
        end.query_offset = 10;
        end.subject_offset = 100;
        end.query_end = 20;
        end.subject_end = 112;
        end.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![
            (GapAlignOpType::Sub, 3),
            (GapAlignOpType::Del, 2),
            (GapAlignOpType::Sub, 5),
            (GapAlignOpType::Ins, 1),
        ]));

        s_cut_off_gap_edit_script(&mut end, 15, 105, false);
        assert_eq!((end.query_end, end.subject_end), (15, 107));
        assert_eq!(
            end.edit_script.as_ref().unwrap().ops_vec(),
            vec![
                (GapAlignOpType::Sub, 3),
                (GapAlignOpType::Del, 2),
                (GapAlignOpType::Sub, 2)
            ]
        );

        let mut untouched = end.clone();
        s_cut_off_gap_edit_script(&mut untouched, 1000, 1000, false);
        assert_eq!(untouched.query_end, end.query_end);
        assert_eq!(
            untouched.edit_script.unwrap().ops_vec(),
            end.edit_script.unwrap().ops_vec()
        );
    }

    #[test]
    fn cut_off_gap_edit_script_does_not_advance_on_frameshift_ops() {
        let mut hsp = make_hsp(50, 1e-8);
        hsp.query_offset = 0;
        hsp.subject_offset = 0;
        hsp.query_end = 6;
        hsp.subject_end = 6;
        hsp.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![
            (GapAlignOpType::Sub, 1),
            (GapAlignOpType::Del2, 1),
            (GapAlignOpType::Ins1, 1),
            (GapAlignOpType::Sub, 5),
        ]));

        s_cut_off_gap_edit_script(&mut hsp, 2, 2, true);

        assert_eq!((hsp.query_offset, hsp.subject_offset), (2, 2));
        assert_eq!(
            hsp.edit_script.as_ref().unwrap().ops_vec(),
            vec![(GapAlignOpType::Sub, 4)]
        );
    }

    #[test]
    fn test_blast_hit_list_new_and_free() {
        let mut hitlist = blast_hit_list_new(7);
        assert_eq!(hitlist.hsplist_max, 7);
        assert_eq!(hitlist.low_score, i32::MAX);
        assert_eq!(hitlist.worst_evalue, 0.0);

        let mut list = HspList::new(42);
        list.add_hsp(make_hsp(31, 1e-8));
        assert_eq!(hitlist.blast_hit_list_update(list), 0);
        assert_eq!(hitlist.low_score, 31);
        assert_eq!(hitlist.worst_evalue, 1e-8);
        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(blast_hit_list_hsp_lists_free(Some(&mut hitlist)), 0);
        assert!(hitlist.hsp_lists.is_empty());

        let mut slot = Some(hitlist);
        assert!(blast_hit_list_free(&mut slot).is_none());
        assert!(slot.is_none());
    }

    #[test]
    fn test_blast_hit_list_sort_by_evalue_purges_empty_lists() {
        let mut hitlist = blast_hit_list_new(10);
        let mut worse = HspList::new(1);
        worse.add_hsp(make_hsp(20, 1e-3));
        let mut better = HspList::new(2);
        better.add_hsp(make_hsp(40, 1e-20));
        let empty = HspList::new(3);

        hitlist.hsp_lists.push(worse);
        hitlist.hsp_lists.push(empty);
        hitlist.hsp_lists.push(better);

        assert_eq!(blast_hit_list_sort_by_evalue(&mut hitlist), 0);
        let oids: Vec<i32> = hitlist.hsp_lists.iter().map(|list| list.oid).collect();
        assert_eq!(oids, vec![2, 1]);
    }

    #[test]
    fn test_blast_hit_list_update_respects_hsplist_max_heap() {
        let mut hitlist = blast_hit_list_new(3);
        for &(oid, score, evalue) in &[
            (1, 90, 1e-30),
            (2, 40, 1e-4),
            (3, 70, 1e-20),
            (4, 10, 1e-2),
            (5, 80, 1e-25),
        ] {
            let mut list = HspList::new(oid);
            list.add_hsp(make_hsp(score, evalue));
            assert_eq!(hitlist.blast_hit_list_update(list), 0);
        }

        let mut oids: Vec<i32> = hitlist.hsp_lists.iter().map(|list| list.oid).collect();
        oids.sort_unstable();
        assert_eq!(oids, vec![1, 3, 5]);
        assert_eq!(hitlist.hsp_lists[0].oid, 3);
        assert_eq!(hitlist.worst_evalue, hitlist.hsp_lists[0].best_evalue);
        assert_eq!(hitlist.low_score, hitlist.hsp_lists[0].hsps[0].score);

        assert_eq!(blast_hit_list_sort_by_evalue(&mut hitlist), 0);
        let ordered: Vec<i32> = hitlist.hsp_lists.iter().map(|list| list.oid).collect();
        assert_eq!(ordered, vec![1, 5, 3]);
    }

    #[test]
    fn test_blast_hit_list_update_sorts_full_heap_inputs_like_c() {
        let mut hitlist = blast_hit_list_new(2);
        for &(oid, score, evalue) in &[(1, 80, 1e-20), (2, 60, 1e-10)] {
            let mut list = HspList::new(oid);
            list.add_hsp(make_hsp(score, evalue));
            assert_eq!(hitlist.blast_hit_list_update(list), 0);
        }

        let mut unsorted_better = HspList::new(3);
        unsorted_better.add_hsp(make_hsp(10, 1e-2));
        unsorted_better.add_hsp(make_hsp(90, 1e-30));
        assert_eq!(hitlist.blast_hit_list_update(unsorted_better), 0);

        let oid3 = hitlist
            .hsp_lists
            .iter()
            .find(|list| list.oid == 3)
            .expect("replacement list");
        assert_eq!(oid3.hsps[0].score, 90);
        assert_eq!(oid3.best_evalue, 1e-30);

        let mut oids: Vec<i32> = hitlist.hsp_lists.iter().map(|list| list.oid).collect();
        oids.sort_unstable();
        assert_eq!(oids, vec![1, 3]);
    }

    #[test]
    fn test_blast_hit_list_update_replaces_tied_worst_like_c() {
        let mut hitlist = blast_hit_list_new(1);
        let mut saved = HspList::new(1);
        let mut saved_hsp = make_hsp(28, 1.0e-13);
        saved_hsp.query_offset = 5;
        saved.add_hsp(saved_hsp);
        assert_eq!(hitlist.blast_hit_list_update(saved), 0);

        let mut tied_newer = HspList::new(1);
        let mut tied_hsp = make_hsp(28, 1.0e-13);
        tied_hsp.query_offset = 99;
        tied_newer.add_hsp(tied_hsp);
        assert_eq!(hitlist.blast_hit_list_update(tied_newer), 0);

        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(hitlist.hsp_lists[0].hsps[0].query_offset, 99);
    }

    #[test]
    fn test_null_pointer_purge_translations_are_compact_in_rust() {
        let mut hsp_list = HspList::new(1);
        hsp_list.add_hsp(make_hsp(20, 1e-4));
        hsp_list.add_hsp(make_hsp(40, 1e-30));
        hsp_list.best_evalue = 9.9;
        assert_eq!(blast_hsp_list_purge_null_hsps(Some(&mut hsp_list)), 0);
        assert_eq!(hsp_list.hsps.len(), 2);
        assert_eq!(hsp_list.hsps[0].score, 20);
        assert_eq!(hsp_list.hsps[1].score, 40);
        assert_eq!(hsp_list.best_evalue, 9.9);

        let mut hitlist = blast_hit_list_new(5);
        hitlist.hsp_lists.push(hsp_list);
        let mut second_list = HspList::new(2);
        second_list.add_hsp(make_hsp(60, 1e-50));
        second_list.best_evalue = 1e-50;
        hitlist.hsp_lists.push(second_list);
        hitlist.low_score = 1234;
        hitlist.worst_evalue = 7.7;
        assert_eq!(blast_hit_list_purge_null_hsp_lists(Some(&mut hitlist)), 0);
        assert_eq!(hitlist.hsp_lists.len(), 2);
        assert_eq!(hitlist.hsp_lists[0].hsps[0].score, 20);
        assert_eq!(hitlist.hsp_lists[1].hsps[0].score, 60);
        assert_eq!(hitlist.low_score, 1234);
        assert_eq!(hitlist.worst_evalue, 7.7);

        assert_eq!(blast_hsp_list_purge_null_hsps(None), 0);
        assert_eq!(blast_hit_list_purge_null_hsp_lists(None), 0);
    }

    #[test]
    fn test_blast_hsp_results_new_and_free() {
        let results = blast_hsp_results_new(3);
        assert_eq!(results.hitlists.len(), 3);
        assert!(results.hitlists.iter().all(Option::is_none));

        let negative = blast_hsp_results_new(-1);
        assert!(negative.hitlists.is_empty());

        let mut slot = Some(results);
        assert!(blast_hsp_results_free(&mut slot).is_none());
        assert!(slot.is_none());
    }

    #[test]
    fn test_blast_hsp_results_sort_by_evalue_purges_first_empty_like_c() {
        let mut results = blast_hsp_results_new(1);
        let hitlist = results.hitlists[0].get_or_insert_with(|| blast_hit_list_new(10));

        let mut worse = HspList::new(1);
        worse.add_hsp(make_hsp(20, 1e-3));
        worse.best_evalue = s_blast_get_best_evalue(&worse);
        let empty = HspList::new(2);
        let mut better_after_empty = HspList::new(3);
        better_after_empty.add_hsp(make_hsp(80, 1e-40));
        better_after_empty.best_evalue = s_blast_get_best_evalue(&better_after_empty);

        hitlist.hsp_lists = vec![worse, empty, better_after_empty];

        assert_eq!(blast_hsp_results_sort_by_evalue(&mut results), 0);
        let hitlist = results.hitlists[0].as_ref().unwrap();
        assert_eq!(hitlist.hsp_lists.len(), 2);
        assert_eq!(hitlist.hsp_lists[0].oid, 3);
        assert_eq!(hitlist.hsp_lists[1].oid, 1);
    }

    #[test]
    fn test_blast_hsp_results_insert_hsp_list_edges() {
        let mut results = blast_hsp_results_new(1);

        assert_eq!(
            blast_hsp_results_insert_hsp_list(&mut results, 0, HspList::new(99), 10),
            0
        );
        assert!(results.hitlists[0].is_none());

        let mut negative_index = HspList::new(7);
        negative_index.add_hsp(make_hsp(50, 1e-6));
        assert_eq!(
            blast_hsp_results_insert_hsp_list(&mut results, -1, negative_index, 10),
            1
        );
        assert!(results.hitlists[0].is_none());

        let mut out_of_range = HspList::new(8);
        out_of_range.add_hsp(make_hsp(60, 1e-8));
        assert_eq!(
            blast_hsp_results_insert_hsp_list(&mut results, 2, out_of_range, 10),
            1
        );
        assert!(results.hitlists[0].is_none());

        let mut valid = HspList::new(9);
        valid.add_hsp(make_hsp(70, 1e-10));
        assert_eq!(
            blast_hsp_results_insert_hsp_list(&mut results, 0, valid, 10),
            0
        );
        let hitlist = results.hitlists[0].as_ref().expect("hitlist");
        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(hitlist.hsplist_max, 10);
        assert_eq!(hitlist.hsp_lists[0].oid, 9);
    }

    #[test]
    fn test_blast_hsp_results_insert_uses_hitlist_size_like_c() {
        let mut results = blast_hsp_results_new(1);
        for (oid, score, evalue) in [(1, 20, 1e-2), (2, 50, 1e-20), (3, 30, 1e-5)] {
            let mut list = HspList::new(oid);
            list.add_hsp(make_hsp(score, evalue));
            assert_eq!(
                blast_hsp_results_insert_hsp_list(&mut results, 0, list, 2),
                0
            );
        }

        let hitlist = results.hitlists[0].as_mut().expect("hitlist");
        let _ = blast_hit_list_sort_by_evalue(hitlist);
        assert_eq!(hitlist.hsplist_max, 2);
        assert_eq!(
            hitlist
                .hsp_lists
                .iter()
                .map(|list| list.oid)
                .collect::<Vec<_>>(),
            vec![2, 3]
        );
    }

    #[test]
    fn test_sthread_local_data_array_lifecycle_and_setup() {
        let mut array = sthread_local_data_array_new(2);
        assert_eq!(array.num_elems, 2);
        assert_eq!(array.tld.len(), 2);
        assert!(array.tld.iter().all(Option::is_some));

        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[11, 7]);
        assert_eq!(sthread_local_data_array_setup(None, &query_info, 9), -1);
        assert_eq!(
            sthread_local_data_array_setup(Some(&mut array), &query_info, 9),
            0
        );

        for tld in array.tld.iter().flatten() {
            let local_query = tld.query_info.as_ref().expect("thread query info");
            assert_eq!(local_query.num_queries, 2);
            assert_eq!(local_query.max_length, 11);
            assert_eq!(tld.hit_params.as_ref().unwrap().options.hitlist_size, 9);
            assert_eq!(
                tld.results.as_ref().expect("thread results").hitlists.len(),
                2
            );
        }

        let mut broken = sthread_local_data_array_new(2);
        broken.tld[1] = None;
        assert_eq!(
            sthread_local_data_array_setup(Some(&mut broken), &query_info, 9),
            -1
        );

        sthread_local_data_array_trim(Some(&mut array), 3);
        assert_eq!(array.num_elems, 2);
        sthread_local_data_array_trim(None, 1);
        sthread_local_data_array_trim(Some(&mut array), 1);
        assert_eq!(array.num_elems, 1);
        assert_eq!(array.tld.len(), 1);

        let mut slot = Some(array);
        assert!(sthread_local_data_array_free(&mut slot).is_none());
        assert!(slot.is_none());
    }

    #[test]
    fn test_sthread_local_data_array_consolidates_results_by_query() {
        let mut array = sthread_local_data_array_new(2);
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[10, 20]);
        assert_eq!(
            sthread_local_data_array_setup(Some(&mut array), &query_info, 12),
            0
        );

        let mut q0_thread0 = blast_hit_list_new(12);
        let mut list10 = HspList::new(10);
        list10.add_hsp(make_hsp(80, 1e-20));
        q0_thread0.hsp_lists.push(list10);
        q0_thread0.hsp_lists.push(HspList::new(999));
        q0_thread0.low_score = 80;
        q0_thread0.worst_evalue = 1e-20;
        array.tld[0]
            .as_mut()
            .unwrap()
            .results
            .as_mut()
            .unwrap()
            .hitlists[0] = Some(q0_thread0);

        let mut q0_thread1 = blast_hit_list_new(12);
        let mut list20 = HspList::new(20);
        list20.add_hsp(make_hsp(70, 1e-10));
        q0_thread1.hsp_lists.push(list20);
        q0_thread1.low_score = 70;
        q0_thread1.worst_evalue = 1e-10;
        array.tld[1]
            .as_mut()
            .unwrap()
            .results
            .as_mut()
            .unwrap()
            .hitlists[0] = Some(q0_thread1);

        let mut q1_thread1 = blast_hit_list_new(12);
        let mut list30 = HspList::new(30);
        list30.add_hsp(make_hsp(50, 1e-5));
        q1_thread1.hsp_lists.push(list30);
        array.tld[1]
            .as_mut()
            .unwrap()
            .results
            .as_mut()
            .unwrap()
            .hitlists[1] = Some(q1_thread1);

        let consolidated =
            sthread_local_data_array_consolidate_results(Some(&mut array)).expect("results");
        let q0 = consolidated.hitlists[0].as_ref().expect("query 0 hitlist");
        let q1 = consolidated.hitlists[1].as_ref().expect("query 1 hitlist");

        let q0_oids: Vec<i32> = q0.hsp_lists.iter().map(|list| list.oid).collect();
        assert_eq!(q0_oids, vec![10, 20]);
        assert_eq!(q0.low_score, 70);
        assert_eq!(q0.worst_evalue, 1e-10);
        assert_eq!(q1.hsp_lists.len(), 1);
        assert_eq!(q1.hsp_lists[0].oid, 30);

        for tld in array.tld.iter().flatten() {
            let results = tld.results.as_ref().unwrap();
            for hitlist in results.hitlists.iter().flatten() {
                assert!(hitlist.hsp_lists.is_empty());
            }
        }
    }

    #[test]
    fn test_hsp_stream_result_batch_lifecycle_and_append() {
        let mut batch = blast_hsp_stream_result_batch_init(2);
        assert_eq!(batch.hsplist_array.len(), 2);
        batch.hsplist_array[0] = Some(HspList::new(42));
        batch.num_hsplists = 1;
        let mut batch_slot = Some(batch);

        let reset = Blast_HSPStreamResultBatchReset(batch_slot.as_mut()).expect("reset batch");
        assert_eq!(reset.num_hsplists, 0);
        assert!(reset.hsplist_array[0].is_none());

        let mut batches = blast_hsp_stream_results_batch_new();
        assert_eq!(batches.num_allocated, 1);
        let default_batches = s_blast_hsp_stream_results_batch_array_new(0);
        assert_eq!(default_batches.num_allocated, 10);
        assert_eq!(default_batches.array_of_batches.len(), 10);
        assert_eq!(
            s_blast_hsp_stream_results_batch_array_append(
                Some(&mut batches),
                Some(blast_hsp_stream_result_batch_init(1)),
            ),
            0
        );
        assert_eq!(
            s_blast_hsp_stream_results_batch_array_append(
                Some(&mut batches),
                Some(blast_hsp_stream_result_batch_init(1)),
            ),
            0
        );
        assert_eq!(batches.num_batches, 2);
        assert_eq!(batches.num_allocated, 2);
        assert_eq!(
            s_blast_hsp_stream_results_batch_array_append(Some(&mut batches), None),
            -1
        );
        assert_eq!(
            s_blast_hsp_stream_results_batch_array_append(
                None,
                Some(blast_hsp_stream_result_batch_init(1)),
            ),
            -1
        );

        let mut batches_slot = Some(batches);
        s_blast_hsp_stream_results_batch_array_reset(&mut batches_slot);
        assert_eq!(batches_slot.as_ref().unwrap().num_batches, 0);
        assert!(blast_hsp_stream_results_batch_array_free(batches_slot).is_none());
    }

    #[test]
    fn test_blast_hsp_list_dup_and_swap() {
        let mut list1 = HspList::new(1);
        list1.add_hsp(make_hsp(20, 1e-4));
        let mut list2 = HspList::new(2);
        list2.add_hsp(make_hsp(40, 1e-8));

        let copy = blast_hsp_list_dup(Some(&list1)).expect("copy");
        assert_eq!(copy.oid, 1);
        assert_eq!(copy.hsps.len(), 1);
        assert_eq!(copy.hsps[0].score, 20);
        assert!(blast_hsp_list_dup(None).is_none());

        blast_hsp_list_swap(&mut list1, &mut list2);
        assert_eq!(list1.oid, 2);
        assert_eq!(list2.oid, 1);
    }

    #[test]
    fn test_blast_hsp_list_sort_by_score_matches_c_order() {
        let mut list = BlastHSPList {
            oid: 0,
            query_index: 0,
            hsp_array: Vec::new(),
            hspcnt: 0,
            allocated: 0,
            hsp_max: 0,
            do_not_reallocate: false,
            best_evalue: f64::MAX,
        };
        let mut hsp_a = make_hsp(40, 1e-4);
        hsp_a.subject_offset = 10;
        hsp_a.subject_end = 20;
        hsp_a.query_offset = 5;
        hsp_a.query_end = 15;
        let mut hsp_b = make_hsp(60, 1e-6);
        hsp_b.subject_offset = 30;
        hsp_b.subject_end = 40;
        hsp_b.query_offset = 7;
        hsp_b.query_end = 17;
        let mut hsp_c = make_hsp(40, 1e-5);
        hsp_c.subject_offset = 5;
        hsp_c.subject_end = 30;
        hsp_c.query_offset = 9;
        hsp_c.query_end = 19;
        list.hsp_array = vec![
            Some(BlastHSP::from_legacy_hsp(hsp_a)),
            Some(BlastHSP::from_legacy_hsp(hsp_b)),
            Some(BlastHSP::from_legacy_hsp(hsp_c)),
        ];
        list.hspcnt = list.hsp_array.len() as i32;
        list.allocated = list.hspcnt;
        list.best_evalue = 1e-6;

        assert_eq!(blast_hsp_list_sort_by_score(Some(&mut list)), 0);
        let ordered: Vec<(i32, i32, i32)> = list
            .hsp_array
            .iter()
            .filter_map(|hsp| hsp.as_ref())
            .map(|hsp| (hsp.score, hsp.subject.offset, hsp.subject.end))
            .collect();
        assert_eq!(ordered, vec![(60, 30, 40), (40, 5, 30), (40, 10, 20)]);
        assert_eq!(list.best_evalue, 1e-6);
        assert_eq!(blast_hsp_list_sort_by_score(None), 0);
    }

    #[test]
    fn test_blast_hsp_list_sort_by_evalue_matches_c_order() {
        let mut list = blast_hsp_list_new(0);
        let mut hsp_a = make_hsp(40, 1e-20);
        hsp_a.subject_offset = 20;
        hsp_a.subject_end = 50;
        let mut hsp_b = make_hsp(80, 1e-5);
        hsp_b.subject_offset = 10;
        hsp_b.subject_end = 30;
        let mut hsp_c = make_hsp(60, 1e-190);
        hsp_c.subject_offset = 30;
        hsp_c.subject_end = 40;
        let mut hsp_d = make_hsp(70, 1e-200);
        hsp_d.subject_offset = 5;
        hsp_d.subject_end = 25;
        list.hsps = vec![hsp_a, hsp_b, hsp_c, hsp_d];
        list.best_evalue = s_blast_get_best_evalue(&list);

        assert_eq!(blast_hsp_list_sort_by_evalue(Some(&mut list)), 0);
        let ordered: Vec<(i32, f64)> = list
            .hsps
            .iter()
            .map(|hsp| (hsp.score, hsp.evalue))
            .collect();
        assert_eq!(
            ordered,
            vec![(70, 1e-200), (60, 1e-190), (40, 1e-20), (80, 1e-5)]
        );
        assert_eq!(list.best_evalue, 1e-200);
        assert_eq!(blast_hsp_list_sort_by_evalue(None), 0);
    }

    #[test]
    fn test_blast_hsp_results_reverse_sort_and_order() {
        let mut results = blast_hsp_results_new(1);
        for &(oid, evalue) in &[(1, 1e-10), (2, 1e-2), (3, 1e-20)] {
            let mut list = HspList::new(oid);
            list.add_hsp(make_hsp(50, evalue));
            assert_eq!(
                blast_hsp_results_insert_hsp_list(&mut results, 0, list, 0),
                0
            );
        }

        assert_eq!(blast_hsp_results_reverse_sort(&mut results), 0);
        let hitlist = results.hitlists[0].as_ref().expect("hitlist");
        let reverse_sorted_oids: Vec<i32> = hitlist.hsp_lists.iter().map(|l| l.oid).collect();
        assert_eq!(reverse_sorted_oids, vec![2, 1, 3]);

        assert_eq!(blast_hsp_results_reverse_order(&mut results), 0);
        let hitlist = results.hitlists[0].as_ref().expect("hitlist");
        let reversed_oids: Vec<i32> = hitlist.hsp_lists.iter().map(|l| l.oid).collect();
        assert_eq!(reversed_oids, vec![3, 1, 2]);
    }

    #[test]
    fn test_blast_hsp_results_reverse_sort_purges_from_first_empty_like_c() {
        let mut results = blast_hsp_results_new(1);
        let mut hitlist = blast_hit_list_new(10);

        let mut non_empty = HspList::new(1);
        non_empty.add_hsp(make_hsp(50, 1e-20));
        hitlist.hsp_lists.push(non_empty);
        hitlist.hsp_lists.push(HspList::new(2));

        results.hitlists[0] = Some(hitlist);
        assert_eq!(blast_hsp_results_reverse_sort(&mut results), 0);

        let hitlist = results.hitlists[0].as_ref().expect("hitlist");
        assert!(hitlist.hsp_lists.is_empty());
    }

    #[test]
    fn traceback_get_encoding_matches_subject_type() {
        assert_eq!(
            blast_traceback_get_encoding(crate::program::BLASTP),
            crate::seqsrc::SeqEncoding::Protein
        );
        assert_eq!(
            blast_traceback_get_encoding(crate::program::BLASTX),
            crate::seqsrc::SeqEncoding::Protein
        );
        assert_eq!(
            blast_traceback_get_encoding(crate::program::RPS_TBLASTN),
            crate::seqsrc::SeqEncoding::Protein
        );
        assert_eq!(
            blast_traceback_get_encoding(crate::program::TBLASTN),
            crate::seqsrc::SeqEncoding::Ncbi4na
        );
        assert_eq!(
            blast_traceback_get_encoding(crate::program::TBLASTX),
            crate::seqsrc::SeqEncoding::Ncbi4na
        );
        assert_eq!(
            blast_traceback_get_encoding(crate::program::PSI_TBLASTN),
            crate::seqsrc::SeqEncoding::Ncbi4na
        );
        assert_eq!(
            blast_traceback_get_encoding(crate::program::BLASTN),
            crate::seqsrc::SeqEncoding::Nucleotide
        );
        assert_eq!(
            blast_traceback_get_encoding(crate::program::PHI_BLASTN),
            crate::seqsrc::SeqEncoding::Nucleotide
        );
    }

    #[test]
    fn blast_prune_extra_hits_truncates_each_query_hitlist() {
        let mut results = blast_hsp_results_new(2);
        for query_index in 0..2 {
            for oid in 0..4 {
                let mut list = HspList::new(10 * query_index + oid);
                list.add_hsp(make_hsp(50 - oid, 1e-10 * (oid + 1) as f64));
                assert_eq!(
                    blast_hsp_results_insert_hsp_list(&mut results, query_index, list, 0),
                    0
                );
            }
        }

        s_blast_prune_extra_hits(&mut results, 2);
        for hitlist in results.hitlists.iter().flatten() {
            assert_eq!(hitlist.hsp_lists.len(), 2);
        }
        assert_eq!(
            results.hitlists[0]
                .as_ref()
                .unwrap()
                .hsp_lists
                .iter()
                .map(|list| list.oid)
                .collect::<Vec<_>>(),
            vec![0, 1]
        );

        s_blast_prune_extra_hits(&mut results, -1);
        assert!(results
            .hitlists
            .iter()
            .flatten()
            .all(|hitlist| hitlist.hsp_lists.is_empty()));
    }

    #[test]
    fn test_blast_hsp_list_save_hsp_respects_hsp_max() {
        let mut list = blast_hsp_list_new(2);
        assert_eq!(blast_hsp_list_save_hsp(&mut list, make_hsp(10, 1e-1)), 0);
        assert_eq!(blast_hsp_list_save_hsp(&mut list, make_hsp(30, 1e-3)), 0);
        assert_eq!(blast_hsp_list_save_hsp(&mut list, make_hsp(20, 1e-2)), 0);

        list.hsps.sort_by(score_compare_hsps);
        let scores: Vec<i32> = list.hsps.iter().map(|h| h.score).collect();
        assert_eq!(scores, vec![30, 20]);
        assert_eq!(list.best_evalue, 1e-3);
    }

    #[test]
    fn test_heap_insert_helpers_keep_worst_at_root() {
        let mut hsps = vec![
            make_hsp(90, 1e-30),
            make_hsp(20, 1e-2),
            make_hsp(70, 1e-20),
            make_hsp(10, 1e-1),
        ];
        s_create_heap(&mut hsps, score_compare_hsps);
        assert_eq!(hsps[0].score, 10);

        let mut list = blast_hsp_list_new(3);
        list.hsps = vec![make_hsp(90, 1e-30), make_hsp(20, 1e-2), make_hsp(70, 1e-20)];
        list.best_evalue = s_blast_get_best_evalue(&list);
        s_create_heap(&mut list.hsps, score_compare_hsps);
        s_blast_hsp_list_insert_hsp_in_heap(&mut list, make_hsp(80, 1e-25));
        let mut scores: Vec<i32> = list.hsps.iter().map(|hsp| hsp.score).collect();
        scores.sort_unstable_by(|a, b| b.cmp(a));
        assert_eq!(scores, vec![90, 80, 70]);
        assert_eq!(list.best_evalue, 1e-30);
    }

    #[test]
    fn test_hit_list_heap_insert_updates_worst_summary() {
        let mut hitlist = blast_hit_list_new(3);
        for &(oid, score, evalue) in &[(1, 90, 1e-30), (2, 40, 1e-4), (3, 70, 1e-20)] {
            let mut list = HspList::new(oid);
            list.add_hsp(make_hsp(score, evalue));
            list.best_evalue = s_blast_get_best_evalue(&list);
            hitlist.hsp_lists.push(list);
        }
        s_create_heap(&mut hitlist.hsp_lists, evalue_compare_hsp_lists);
        assert_eq!(hitlist.hsp_lists[0].oid, 2);

        let mut replacement = HspList::new(4);
        replacement.add_hsp(make_hsp(80, 1e-25));
        replacement.best_evalue = s_blast_get_best_evalue(&replacement);
        s_blast_hit_list_insert_hsp_list_in_heap(&mut hitlist, replacement);

        let mut oids: Vec<i32> = hitlist.hsp_lists.iter().map(|list| list.oid).collect();
        oids.sort_unstable();
        assert_eq!(oids, vec![1, 3, 4]);
        assert_eq!(hitlist.worst_evalue, hitlist.hsp_lists[0].best_evalue);
        assert_eq!(hitlist.low_score, hitlist.hsp_lists[0].hsps[0].score);
    }

    #[test]
    fn test_blast_hsp_list_phi_get_bit_scores_uses_phi_formula() {
        let mut list = HspList::new(1);
        list.add_hsp(make_hsp(50, 1e-8));

        blast_hsp_list_phi_get_bit_scores(&mut list, 0.267, 0.041);

        let expected = (50.0_f64 * 0.267 - 0.041_f64.ln() - (1.0 + 50.0_f64 * 0.267).ln())
            / crate::math::NCBIMATH_LN2;
        assert!((list.hsps[0].bit_score - expected).abs() < 1e-12);
    }

    #[test]
    fn test_phi_effective_patterns_and_evalues_use_c_formula() {
        let occurrences = [0, 3, 9, 13, 22];
        assert_eq!(phi_blast_get_effective_number_of_patterns(&[], 8), 0);
        assert_eq!(phi_blast_get_effective_number_of_patterns(&[7], 8), 1);
        assert_eq!(
            phi_blast_get_effective_number_of_patterns(&occurrences, 8),
            3
        );

        let mut hsp = make_hsp(42, 0.0);
        s_hsp_phi_get_evalue(&mut hsp, 0.267, 0.041, 3, 17);
        let expected =
            0.041_f64 * (1.0 + 0.267_f64 * 42.0) * 3.0 * 17.0 * (-0.267_f64 * 42.0).exp();
        assert!((hsp.evalue - expected).abs() < 1e-15);

        let mut list = HspList::new(1);
        list.add_hsp(make_hsp(50, 0.0));
        list.add_hsp(make_hsp(30, 0.0));
        blast_hsp_list_phi_get_evalues(&mut list, 0.267, 0.041, &occurrences, 8, 17);
        assert_eq!(
            list.best_evalue,
            list.hsps[0].evalue.min(list.hsps[1].evalue)
        );
        assert!(list.hsps[0].evalue < list.hsps[1].evalue);
    }

    #[test]
    fn phiblast_hsp_results_split_groups_hsps_by_pattern_index() {
        let mut results = blast_hsp_results_new(1);
        let mut list = blast_hsp_list_new(0);
        list.oid = 42;
        let mut pattern_two = make_hsp(40, 1e-4);
        pattern_two.pat_info = Some(PhiPatInfo {
            index: 2,
            length: 6,
        });
        let mut pattern_zero = make_hsp(80, 1e-8);
        pattern_zero.pat_info = Some(PhiPatInfo {
            index: 0,
            length: 4,
        });
        let missing_pattern = make_hsp(90, 1e-9);
        let mut out_of_range_pattern = make_hsp(100, 1e-10);
        out_of_range_pattern.pat_info = Some(PhiPatInfo {
            index: 5,
            length: 8,
        });
        let _ = blast_hsp_list_save_hsp(&mut list, pattern_two);
        let _ = blast_hsp_list_save_hsp(&mut list, pattern_zero);
        let _ = blast_hsp_list_save_hsp(&mut list, missing_pattern);
        let _ = blast_hsp_list_save_hsp(&mut list, out_of_range_pattern);
        let _ = blast_hsp_results_insert_hsp_list(&mut results, 0, list, 100);

        let split = phiblast_hsp_results_split(Some(&results), 3).expect("split results");

        assert!(split[1].is_none());
        let p0 = split[0].as_ref().expect("pattern 0");
        let p2 = split[2].as_ref().expect("pattern 2");
        assert_eq!(
            p0.hitlists[0].as_ref().unwrap().hsp_lists[0].hsps[0].score,
            80
        );
        assert_eq!(
            p2.hitlists[0].as_ref().unwrap().hsp_lists[0].hsps[0].score,
            40
        );
        assert_eq!(p0.hitlists[0].as_ref().unwrap().hsp_lists[0].hsps.len(), 1);
        assert_eq!(p2.hitlists[0].as_ref().unwrap().hsp_lists[0].hsps.len(), 1);
        assert!(phiblast_hsp_results_split(Some(&results), 0).is_none());

        let empty_split = phiblast_hsp_results_split(None, 2).expect("empty split");
        assert_eq!(empty_split.len(), 2);
        assert!(empty_split.iter().all(Option::is_none));

        let no_hitlist = blast_hsp_results_new(1);
        let no_hitlist_split = phiblast_hsp_results_split(Some(&no_hitlist), 2).expect("split");
        assert_eq!(no_hitlist_split.len(), 2);
        assert!(no_hitlist_split.iter().all(Option::is_none));
    }

    #[test]
    fn hsp_list_post_traceback_update_recomputes_reaps_rescales_and_bits() {
        let mut scoring_options = crate::options::ScoringOptions::new_blastp();
        scoring_options.gapped_calculation = true;
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 2.0);
        let hit_params = crate::parameters::HitSavingParameters {
            options: crate::options::HitSavingOptions {
                expect_value: 1e-6,
                ..Default::default()
            },
            cutoff_score_min: 0,
            low_score: Vec::new(),
            cutoffs: Vec::new(),
            link_hsp_params: None,
            prelim_evalue: 1e-6,

            ..Default::default()
        };
        let mut query_info = crate::queryinfo::QueryInfo::new_blastp(&[100]);
        query_info.contexts[0].eff_searchsp = 10_000;
        let kbp = crate::stat::KarlinBlk {
            lambda: 0.5,
            k: 0.1,
            log_k: 0.1f64.ln(),
            h: 0.2,
            round_down: false,
        };
        let mut sbp =
            crate::stat::blast_score_blk_new(crate::encoding::BLASTAA_SEQ_CODE, 1).unwrap();
        sbp.kbp[0] = kbp.clone();
        sbp.kbp_gap = vec![kbp.clone()];

        let mut legacy = blast_hsp_list_new(0);
        legacy.oid = 7;
        let _ = blast_hsp_list_save_hsp(&mut legacy, make_hsp(50, f64::MAX));
        let _ = blast_hsp_list_save_hsp(&mut legacy, make_hsp(5, f64::MAX));
        let mut list = BlastHSPList::from_legacy_hsp_list(legacy);

        assert_eq!(
            s_hsp_list_post_traceback_update(
                crate::program::BLASTP,
                &mut list,
                &query_info,
                &score_params,
                &hit_params,
                &sbp,
                100,
            ),
            0
        );

        assert_eq!(list.hsp_array.len(), 1);
        let hsp = list.hsp_array[0].as_ref().expect("hsp");
        assert_eq!(hsp.score, 25);
        assert!((hsp.bit_score - kbp.raw_to_bit(25)).abs() < 1e-12);
        assert!(list.best_evalue <= 1e-6);
    }

    #[test]
    fn hsp_list_post_traceback_update_uses_each_hsp_context_statistics() {
        let mut scoring_options = crate::options::ScoringOptions::new_blastp();
        scoring_options.gapped_calculation = true;
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let hit_params = crate::parameters::HitSavingParameters {
            options: crate::options::HitSavingOptions {
                expect_value: 1000.0,
                ..Default::default()
            },
            cutoff_score_min: 0,
            low_score: Vec::new(),
            cutoffs: Vec::new(),
            link_hsp_params: None,
            prelim_evalue: 1000.0,

            ..Default::default()
        };
        let mut query_info = crate::queryinfo::QueryInfo::new_blastp(&[50, 50]);
        query_info.contexts[0].eff_searchsp = 1_000;
        query_info.contexts[1].eff_searchsp = 2_000;
        let kbp0 = crate::stat::KarlinBlk {
            lambda: 0.2,
            k: 0.5,
            log_k: 0.5f64.ln(),
            h: 0.2,
            round_down: false,
        };
        let kbp1 = crate::stat::KarlinBlk {
            lambda: 0.4,
            k: 0.25,
            log_k: 0.25f64.ln(),
            h: 0.2,
            round_down: false,
        };
        let mut sbp =
            crate::stat::blast_score_blk_new(crate::encoding::BLASTAA_SEQ_CODE, 2).unwrap();
        sbp.kbp = vec![kbp0.clone(), kbp1.clone()];
        sbp.kbp_gap = vec![kbp0.clone(), kbp1.clone()];

        let mut legacy = blast_hsp_list_new(0);
        let mut context0 = make_hsp(20, f64::MAX);
        context0.context = 0;
        let mut context1 = make_hsp(20, f64::MAX);
        context1.context = 1;
        let _ = blast_hsp_list_save_hsp(&mut legacy, context0);
        let _ = blast_hsp_list_save_hsp(&mut legacy, context1);
        let mut list = BlastHSPList::from_legacy_hsp_list(legacy);

        assert_eq!(
            s_hsp_list_post_traceback_update(
                crate::program::BLASTP,
                &mut list,
                &query_info,
                &score_params,
                &hit_params,
                &sbp,
                50,
            ),
            0
        );

        let by_context = |context| {
            list.hsp_array
                .iter()
                .flatten()
                .find(|hsp| hsp.context == context)
                .expect("context hsp")
        };
        assert!((by_context(0).evalue - kbp0.raw_to_evalue(20, 1_000.0)).abs() < 1e-12);
        assert!((by_context(1).evalue - kbp1.raw_to_evalue(20, 2_000.0)).abs() < 1e-12);
        assert!((by_context(0).bit_score - kbp0.raw_to_bit(20)).abs() < 1e-12);
        assert!((by_context(1).bit_score - kbp1.raw_to_bit(20)).abs() < 1e-12);
    }

    #[test]
    fn hsp_list_post_traceback_update_uses_get_evalues_and_bits_handoff() {
        let mut scoring_options = crate::options::ScoringOptions::new_blastp();
        scoring_options.gapped_calculation = true;
        let score_params =
            crate::parameters::ScoringParameters::from_options(&scoring_options, 1.0);
        let hit_params = crate::parameters::HitSavingParameters {
            options: crate::options::HitSavingOptions {
                expect_value: 1000.0,
                ..Default::default()
            },
            cutoff_score_min: 0,
            low_score: Vec::new(),
            cutoffs: Vec::new(),
            link_hsp_params: None,
            prelim_evalue: 1000.0,

            ..Default::default()
        };
        let mut query_info = crate::queryinfo::QueryInfo::new_blastp(&[50]);
        query_info.contexts[0].eff_searchsp = 75;
        let kbp = crate::stat::KarlinBlk {
            lambda: 0.3,
            k: 0.2,
            log_k: 0.2f64.ln(),
            h: 0.4,
            round_down: false,
        };
        let mut sbp =
            crate::stat::blast_score_blk_new(crate::encoding::BLASTAA_SEQ_CODE, 1).unwrap();
        sbp.kbp = vec![kbp.clone()];
        sbp.kbp_gap = vec![kbp.clone()];

        let mut legacy = blast_hsp_list_new(0);
        let mut hsp = make_hsp(12, f64::MAX);
        hsp.context = 4;
        let _ = blast_hsp_list_save_hsp(&mut legacy, hsp);
        let mut list = BlastHSPList::from_legacy_hsp_list(legacy);

        assert_eq!(
            s_hsp_list_post_traceback_update(
                crate::program::BLASTP,
                &mut list,
                &query_info,
                &score_params,
                &hit_params,
                &sbp,
                75,
            ),
            0
        );

        assert_eq!(list.hsp_array.len(), 1);
        let hsp = list.hsp_array[0].as_ref().expect("hsp");
        assert_eq!(hsp.context, 4);
        assert!((hsp.evalue - kbp.raw_to_evalue(12, 75.0)).abs() < 1e-12);
        assert_eq!(hsp.bit_score, 12.0 * 0.9);
    }

    #[test]
    fn test_blast_hsp_list_is_sorted_by_score() {
        let mut list = HspList::new(1);
        list.add_hsp(make_hsp(50, 1e-5));
        list.add_hsp(make_hsp(20, 1e-2));
        assert!(blast_hsp_list_is_sorted_by_score(Some(&list)));
        list.hsps.swap(0, 1);
        assert!(!blast_hsp_list_is_sorted_by_score(Some(&list)));
        assert!(blast_hsp_list_is_sorted_by_score(None));
    }

    #[test]
    fn test_traceback_rescale_scores_rounds_and_resorts() {
        let mut list = HspList::new(1);
        let mut later_subject = make_hsp(8, 1e-2);
        later_subject.subject_offset = 10;
        later_subject.subject_end = 18;
        let mut earlier_subject = make_hsp(7, 1e-3);
        earlier_subject.subject_offset = 2;
        earlier_subject.subject_end = 10;
        list.hsps = vec![later_subject, earlier_subject];

        s_hsp_list_rescale_scores(&mut list, 2.0);

        assert_eq!(list.hsps[0].score, 4);
        assert_eq!(list.hsps[1].score, 4);
        assert_eq!(list.hsps[0].subject_offset, 2);
        assert_eq!(list.hsps[1].subject_offset, 10);
    }

    #[test]
    fn test_rps_update_restores_segments_and_edit_script() {
        let mut list = HspList::new(1);
        let mut hsp = make_hsp(40, 1e-6);
        hsp.query_offset = 3;
        hsp.query_end = 13;
        hsp.query_gapped_start = 5;
        hsp.query_frame = -1;
        hsp.subject_offset = 30;
        hsp.subject_end = 44;
        hsp.subject_gapped_start = 31;
        hsp.subject_frame = 2;
        hsp.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![
            (GapAlignOpType::Sub, 4),
            (GapAlignOpType::Ins, 2),
            (GapAlignOpType::Del, 1),
        ]));
        list.hsps.push(hsp);

        s_blast_hsp_list_rps_update(crate::program::RPS_TBLASTN, &mut list);

        let hsp = &list.hsps[0];
        assert_eq!(
            (hsp.query_offset, hsp.query_end, hsp.query_gapped_start),
            (30, 44, 31)
        );
        assert_eq!(
            (
                hsp.subject_offset,
                hsp.subject_end,
                hsp.subject_gapped_start
            ),
            (3, 13, 5)
        );
        assert_eq!(
            (hsp.query_frame, hsp.subject_frame, hsp.context),
            (2, -1, 1)
        );
        assert_eq!(
            hsp.edit_script.as_ref().unwrap().ops_vec(),
            vec![
                (GapAlignOpType::Sub, 4),
                (GapAlignOpType::Del, 2),
                (GapAlignOpType::Ins, 1),
            ]
        );
    }

    #[test]
    fn test_blast_hsp_calc_length_and_gaps() {
        let mut hsp = make_hsp(50, 1e-8);
        hsp.query_offset = 10;
        hsp.query_end = 30;
        hsp.subject_offset = 100;
        hsp.subject_end = 125;
        assert_eq!(blast_hsp_calc_length_and_gaps(&hsp), (25, 0, 0));

        hsp.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![
            (crate::gapinfo::GapAlignOpType::Sub, 8),
            (crate::gapinfo::GapAlignOpType::Del, 3),
            (crate::gapinfo::GapAlignOpType::Ins, 2),
        ]));
        assert_eq!(blast_hsp_calc_length_and_gaps(&hsp), (23, 5, 2));
    }

    #[test]
    fn test_blast_hsp_list_adjust_offsets() {
        let mut list = HspList::new(1);
        list.add_hsp(make_hsp(40, 1e-6));
        blast_hsp_list_adjust_offsets(&mut list, 7);
        assert_eq!(list.hsps[0].subject_offset, 7);
        assert_eq!(list.hsps[0].subject_end, 47);
        assert_eq!(list.hsps[0].subject_gapped_start, 7);
    }

    #[test]
    fn test_adjust_subject_for_translated_sra_search_matches_frame_cases() {
        let mut list = HspList::new(1);
        let mut frame1_at_start = make_hsp(40, 1e-6);
        frame1_at_start.query_offset = 3;
        frame1_at_start.subject_offset = 0;
        frame1_at_start.subject_end = 12;
        frame1_at_start.subject_gapped_start = 0;
        frame1_at_start.subject_frame = 1;
        list.hsps.push(frame1_at_start);

        let mut frame2_no_coordinate_fix = make_hsp(30, 1e-5);
        frame2_no_coordinate_fix.subject_offset = 5;
        frame2_no_coordinate_fix.subject_end = 15;
        frame2_no_coordinate_fix.subject_gapped_start = 6;
        frame2_no_coordinate_fix.subject_frame = 2;
        list.hsps.push(frame2_no_coordinate_fix);

        s_adjust_subject_for_translated_sra_search(&mut list, 1, 100);

        assert_eq!(list.hsps[0].subject_frame, 3);
        assert_eq!(list.hsps[0].query_offset, 4);
        assert_eq!(list.hsps[0].subject_offset, 0);
        assert_eq!(list.hsps[0].subject_end, 11);
        assert_eq!(list.hsps[0].subject_gapped_start, -1);
        assert_eq!(list.hsps[1].subject_frame, 1);
        assert_eq!(list.hsps[1].subject_offset, 5);
        assert_eq!(list.hsps[1].subject_end, 15);

        let mut minus = HspList::new(2);
        let mut hsp = make_hsp(20, 1e-4);
        hsp.query_end = 14;
        hsp.subject_end = 10;
        hsp.subject_frame = -1;
        minus.hsps.push(hsp);

        s_adjust_subject_for_translated_sra_search(&mut minus, 2, 31);

        assert_eq!(minus.hsps[0].subject_end, 9);
        assert_eq!(minus.hsps[0].query_end, 13);
    }

    #[test]
    fn test_adjust_subject_for_translated_sra_search_covers_positive_offset_branches() {
        let mut offset2 = HspList::new(3);
        let mut frame1 = make_hsp(40, 1e-6);
        frame1.subject_offset = 4;
        frame1.subject_end = 12;
        frame1.subject_gapped_start = 5;
        frame1.subject_frame = 1;
        offset2.hsps.push(frame1);
        let mut frame2 = make_hsp(30, 1e-5);
        frame2.subject_offset = 0;
        frame2.subject_end = 10;
        frame2.subject_gapped_start = 0;
        frame2.subject_frame = 2;
        offset2.hsps.push(frame2);
        let mut frame3 = make_hsp(20, 1e-4);
        frame3.subject_offset = 7;
        frame3.subject_end = 16;
        frame3.subject_gapped_start = 8;
        frame3.subject_frame = 3;
        offset2.hsps.push(frame3);

        s_adjust_subject_for_translated_sra_search(&mut offset2, 2, 100);

        assert_eq!(offset2.hsps[0].subject_frame, 2);
        assert_eq!(offset2.hsps[0].subject_offset, 3);
        assert_eq!(offset2.hsps[0].subject_end, 11);
        assert_eq!(offset2.hsps[0].subject_gapped_start, 4);
        assert_eq!(offset2.hsps[1].subject_frame, 3);
        assert_eq!(offset2.hsps[1].query_offset, 1);
        assert_eq!(offset2.hsps[1].subject_offset, 0);
        assert_eq!(offset2.hsps[1].subject_end, 9);
        assert_eq!(offset2.hsps[1].subject_gapped_start, -1);
        assert_eq!(offset2.hsps[2].subject_frame, 1);
        assert_eq!(offset2.hsps[2].subject_offset, 7);
        assert_eq!(offset2.hsps[2].subject_end, 16);

        let mut offset3 = HspList::new(1);
        let mut hsp = make_hsp(10, 1e-3);
        hsp.query_offset = 5;
        hsp.subject_offset = 2;
        hsp.subject_end = 8;
        hsp.subject_gapped_start = 3;
        hsp.subject_frame = 3;
        offset3.hsps.push(hsp);

        s_adjust_subject_for_translated_sra_search(&mut offset3, 3, 100);

        assert_eq!(offset3.hsps[0].subject_frame, 3);
        assert_eq!(offset3.hsps[0].query_offset, 5);
        assert_eq!(offset3.hsps[0].subject_offset, 1);
        assert_eq!(offset3.hsps[0].subject_end, 7);
        assert_eq!(offset3.hsps[0].subject_gapped_start, 2);
    }

    #[test]
    fn test_blast_hsp_list_append_keeps_best_by_score() {
        let mut combined = HspList::new(1);
        combined.add_hsp(make_hsp(50, 1e-5));
        combined.add_hsp(make_hsp(10, 1e-1));
        let mut old = HspList::new(1);
        old.add_hsp(make_hsp(70, 1e-9));
        old.add_hsp(make_hsp(30, 1e-3));

        let mut old_slot = Some(old);
        let mut combined_slot = Some(combined);
        assert_eq!(
            blast_hsp_list_append(&mut old_slot, &mut combined_slot, 3),
            0
        );
        assert!(old_slot.is_none());
        let combined = combined_slot.as_ref().expect("combined");
        let scores: Vec<i32> = combined.hsps.iter().map(|h| h.score).collect();
        assert_eq!(scores, vec![70, 50, 30]);
        assert_eq!(combined.best_evalue, 1e-9);
    }

    #[test]
    fn test_hsp_lists_combine_by_score_merges_sorted_prefix() {
        let mut combined = HspList::new(1);
        combined.hsps = vec![make_hsp(40, 1e-4), make_hsp(90, 1e-30)];
        let mut new_list = HspList::new(1);
        new_list.hsps = vec![make_hsp(50, 1e-5), make_hsp(80, 1e-20), make_hsp(20, 1e-2)];

        s_blast_hsp_lists_combine_by_score(&mut new_list, &mut combined, 4);

        let scores: Vec<i32> = combined.hsps.iter().map(|hsp| hsp.score).collect();
        assert_eq!(scores, vec![90, 80, 50, 40]);
        assert!(new_list.hsps.is_empty());
        assert_eq!(combined.best_evalue, 1e-30);
    }

    #[test]
    fn test_blast_hsp_results_apply_masklevel_rebuilds_by_raw_score() {
        let mut results = blast_hsp_results_new(1);

        let mut list1 = HspList::new(10);
        let mut strong = make_hsp(100, 1e-20);
        strong.query_offset = 10;
        strong.query_end = 60;
        strong.subject_offset = 100;
        strong.subject_end = 150;
        strong.context = 0;
        list1.add_hsp(strong);

        let mut weak_masked = make_hsp(80, 1e-10);
        weak_masked.query_offset = 20;
        weak_masked.query_end = 50;
        weak_masked.subject_offset = 200;
        weak_masked.subject_end = 230;
        weak_masked.context = 0;
        list1.add_hsp(weak_masked);
        assert_eq!(
            blast_hsp_results_insert_hsp_list(&mut results, 0, list1, 0),
            0
        );

        let mut list2 = HspList::new(20);
        let mut different_context = make_hsp(70, 1e-6);
        different_context.query_offset = 20;
        different_context.query_end = 50;
        different_context.subject_offset = 300;
        different_context.subject_end = 330;
        different_context.context = 1;
        list2.add_hsp(different_context);
        assert_eq!(
            blast_hsp_results_insert_hsp_list(&mut results, 0, list2, 0),
            0
        );

        assert_eq!(blast_hsp_results_apply_masklevel(&mut results, 80, 100), 0);
        let hitlist = results.hitlists[0].as_ref().expect("hitlist");
        assert_eq!(hitlist.hsp_lists.len(), 2);
        let list1 = hitlist
            .hsp_lists
            .iter()
            .find(|list| list.oid == 10)
            .unwrap();
        let list2 = hitlist
            .hsp_lists
            .iter()
            .find(|list| list.oid == 20)
            .unwrap();
        assert_eq!(list1.hsps.len(), 1);
        assert_eq!(list1.hsps[0].score, 100);
        assert_eq!(list1.best_evalue, 1e-20);
        assert_eq!(list2.hsps.len(), 1);
        assert_eq!(list2.hsps[0].score, 70);
    }

    #[test]
    fn test_blast_hsp_results_apply_masklevel_removes_empty_lists() {
        let mut results = blast_hsp_results_new(1);

        let mut keeper = HspList::new(1);
        let mut strong = make_hsp(100, 1e-30);
        strong.query_offset = 0;
        strong.query_end = 100;
        strong.context = 0;
        keeper.add_hsp(strong);
        assert_eq!(
            blast_hsp_results_insert_hsp_list(&mut results, 0, keeper, 0),
            0
        );

        let mut removed = HspList::new(2);
        let mut weak = make_hsp(10, 1e-2);
        weak.query_offset = 10;
        weak.query_end = 90;
        weak.context = 0;
        removed.add_hsp(weak);
        assert_eq!(
            blast_hsp_results_insert_hsp_list(&mut results, 0, removed, 0),
            0
        );

        assert_eq!(blast_hsp_results_apply_masklevel(&mut results, 100, 100), 0);
        let hitlist = results.hitlists[0].as_ref().expect("hitlist");
        let oids: Vec<i32> = hitlist.hsp_lists.iter().map(|list| list.oid).collect();
        assert_eq!(oids, vec![1]);
    }

    #[test]
    fn test_s_trim_results_by_total_hsp_limit_uses_ncbi_allowance() {
        let mut results = blast_hsp_results_new(1);
        let mut hitlist = blast_hit_list_new(10);

        for &(oid, count) in &[(1, 1), (2, 3), (3, 5)] {
            let mut list = HspList::new(oid);
            for idx in 0..count {
                list.add_hsp(make_hsp(100 - idx, 1e-10 * (idx + 1) as f64));
            }
            hitlist.hsp_lists.push(list);
        }
        results.hitlists[0] = Some(hitlist);

        assert!(s_trim_results_by_total_hsp_limit(&mut results, 6));
        let hitlist = results.hitlists[0].as_ref().expect("hitlist");
        let counts: Vec<(i32, usize)> = hitlist
            .hsp_lists
            .iter()
            .map(|list| (list.oid, list.hsps.len()))
            .collect();
        assert_eq!(counts, vec![(1, 1), (2, 3), (3, 2)]);
    }

    #[test]
    fn test_s_trim_results_by_total_hsp_limit_zero_is_noop() {
        let mut results = blast_hsp_results_new(1);
        let mut hitlist = blast_hit_list_new(10);
        let mut list = HspList::new(1);
        list.add_hsp(make_hsp(10, 1e-2));
        hitlist.hsp_lists.push(list);
        results.hitlists[0] = Some(hitlist);

        assert!(!s_trim_results_by_total_hsp_limit(&mut results, 0));
        assert_eq!(
            results.hitlists[0].as_ref().unwrap().hsp_lists[0]
                .hsps
                .len(),
            1
        );
    }

    #[test]
    fn test_s_trim_results_by_total_hsp_limit_ex_rebuilds_best_hsps_by_oid() {
        let mut results = blast_hsp_results_new(2);
        let mut hitlist = blast_hit_list_new(9);

        let mut list_b = HspList::new(20);
        list_b.add_hsp(make_hsp(90, 1e-4));
        list_b.add_hsp(make_hsp(60, 1e-20));
        list_b.add_hsp(make_hsp(50, 1e-3));
        hitlist.hsp_lists.push(list_b);

        let mut list_a = HspList::new(10);
        list_a.add_hsp(make_hsp(80, 1e-10));
        list_a.add_hsp(make_hsp(70, 1e-30));
        hitlist.hsp_lists.push(list_a);

        let mut list_c = HspList::new(30);
        list_c.add_hsp(make_hsp(40, 1e-2));
        hitlist.hsp_lists.push(list_c);
        results.hitlists[0] = Some(hitlist);

        let mut flags = vec![true, true];
        assert!(s_trim_results_by_total_hsp_limit_ex(
            &mut results,
            3,
            Some(&mut flags)
        ));
        assert_eq!(flags, vec![true, false]);

        let hitlist = results.hitlists[0].as_ref().expect("hitlist");
        assert_eq!(hitlist.hsplist_max, 9);
        let summary: Vec<(i32, Vec<i32>)> = hitlist
            .hsp_lists
            .iter()
            .map(|list| {
                (
                    list.oid,
                    list.hsps.iter().map(|hsp| hsp.score).collect::<Vec<_>>(),
                )
            })
            .collect();
        assert_eq!(summary, vec![(10, vec![70, 80]), (20, vec![60])]);
    }

    #[test]
    fn test_s_trim_results_by_total_hsp_limit_ex_reports_without_flags() {
        let mut results = blast_hsp_results_new(1);
        let mut hitlist = blast_hit_list_new(9);
        let mut list = HspList::new(10);
        list.add_hsp(make_hsp(80, 1.0e-10));
        list.add_hsp(make_hsp(70, 1.0e-20));
        list.add_hsp(make_hsp(60, 1.0e-30));
        hitlist.hsp_lists.push(list);
        results.hitlists[0] = Some(hitlist);

        assert!(s_trim_results_by_total_hsp_limit_ex(&mut results, 2, None));

        let hitlist = results.hitlists[0].as_ref().expect("hitlist");
        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(
            hitlist.hsp_lists[0]
                .hsps
                .iter()
                .map(|hsp| hsp.score)
                .collect::<Vec<_>>(),
            vec![60, 70]
        );
    }

    #[test]
    fn test_s_trim_results_by_total_hsp_limit_ex_noop_clears_exceeded_flags() {
        let mut results = blast_hsp_results_new(2);
        let mut hitlist = blast_hit_list_new(9);
        let mut list = HspList::new(10);
        list.add_hsp(make_hsp(80, 1.0e-10));
        list.add_hsp(make_hsp(70, 1.0e-20));
        hitlist.hsp_lists.push(list);
        results.hitlists[0] = Some(hitlist);

        let mut flags = vec![true, true];
        assert!(!s_trim_results_by_total_hsp_limit_ex(
            &mut results,
            3,
            Some(&mut flags)
        ));

        assert_eq!(flags, vec![false, false]);
        let hitlist = results.hitlists[0].as_ref().expect("hitlist");
        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(hitlist.hsp_lists[0].hsps.len(), 2);
    }

    #[test]
    fn test_s_trim_results_by_total_hsp_limit_ex_tolerates_short_flag_slice() {
        let mut results = blast_hsp_results_new(2);
        for query_index in 0..2 {
            let mut hitlist = blast_hit_list_new(9);
            let mut list = HspList::new(query_index as i32);
            list.add_hsp(make_hsp(80, 1.0e-10));
            list.add_hsp(make_hsp(70, 1.0e-20));
            list.add_hsp(make_hsp(60, 1.0e-30));
            hitlist.hsp_lists.push(list);
            results.hitlists[query_index] = Some(hitlist);
        }

        let mut flags = vec![false];
        assert!(s_trim_results_by_total_hsp_limit_ex(
            &mut results,
            2,
            Some(&mut flags)
        ));

        assert_eq!(flags, vec![true]);
        for hitlist in results.hitlists.iter().flatten() {
            assert_eq!(hitlist.hsp_lists.len(), 1);
            assert_eq!(hitlist.hsp_lists[0].hsps.len(), 2);
        }
    }

    #[test]
    fn test_filter_blast_results_trims_and_reaps_query_coverage() {
        let mut results = blast_hsp_results_new(1);
        let mut hitlist = blast_hit_list_new(10);

        let mut list = HspList::new(11);
        let mut high_cov = make_hsp(100, 1e-20);
        high_cov.query_offset = 0;
        high_cov.query_end = 80;
        high_cov.context = 0;
        let mut low_cov = make_hsp(90, 1e-10);
        low_cov.query_offset = 0;
        low_cov.query_end = 10;
        low_cov.context = 0;
        let mut would_be_trimmed = make_hsp(80, 1e-5);
        would_be_trimmed.query_offset = 0;
        would_be_trimmed.query_end = 90;
        would_be_trimmed.context = 0;
        list.hsps = vec![high_cov, low_cov, would_be_trimmed];
        list.best_evalue = s_blast_get_best_evalue(&list);
        hitlist.hsp_lists.push(list);

        let mut emptied = HspList::new(22);
        let mut weak_coverage = make_hsp(50, 1e-2);
        weak_coverage.query_offset = 0;
        weak_coverage.query_end = 5;
        weak_coverage.context = 0;
        emptied.hsps.push(weak_coverage);
        emptied.best_evalue = s_blast_get_best_evalue(&emptied);
        hitlist.hsp_lists.push(emptied);
        results.hitlists[0] = Some(hitlist);

        let hit_options = HitSavingOptions {
            max_hsps_per_subject: 2,
            query_cov_hsp_perc: 50.0,
            ..HitSavingOptions::default()
        };
        s_filter_blast_results(&mut results, &hit_options, &[100], None);

        let hitlist = results.hitlists[0].as_ref().expect("hitlist");
        let summary: Vec<(i32, Vec<i32>)> = hitlist
            .hsp_lists
            .iter()
            .map(|list| {
                (
                    list.oid,
                    list.hsps.iter().map(|hsp| hsp.score).collect::<Vec<_>>(),
                )
            })
            .collect();
        assert_eq!(summary, vec![(11, vec![100])]);
    }

    #[test]
    fn test_blast_hsp_list_subject_best_hit_matches_c_ranges() {
        let mut list = blast_hsp_list_new(0);
        let mut best = make_hsp(100, 1e-30);
        best.query_offset = 10;
        best.query_end = 50;
        best.context = 0;
        let mut contained = make_hsp(80, 1e-20);
        contained.query_offset = 12;
        contained.query_end = 48;
        contained.context = 0;
        let mut outside = make_hsp(70, 1e-10);
        outside.query_offset = 52;
        outside.query_end = 70;
        outside.context = 0;
        list.hsps = vec![best, contained, outside];

        assert_eq!(
            blast_hsp_list_subject_best_hit(crate::program::BLASTP, Some(&mut list), 3, &[100]),
            2
        );
        assert_eq!(
            list.hsps.iter().map(|hsp| hsp.score).collect::<Vec<_>>(),
            vec![100, 70]
        );
    }

    #[test]
    fn test_blast_hsp_list_subject_best_hit_blastn_flips_opposite_strand() {
        let mut list = blast_hsp_list_new(0);
        let mut forward = make_hsp(100, 1e-30);
        forward.query_offset = 10;
        forward.query_end = 50;
        forward.context = 0;
        forward.query_frame = 1;
        let mut opposite_contained = make_hsp(80, 1e-20);
        opposite_contained.query_offset = 48;
        opposite_contained.query_end = 90;
        opposite_contained.context = 1;
        opposite_contained.query_frame = -1;
        let mut opposite_outside = make_hsp(70, 1e-10);
        opposite_outside.query_offset = 10;
        opposite_outside.query_end = 40;
        opposite_outside.context = 1;
        opposite_outside.query_frame = -1;
        list.hsps = vec![forward, opposite_contained, opposite_outside];

        assert_eq!(
            blast_hsp_list_subject_best_hit(
                crate::program::BLASTN,
                Some(&mut list),
                3,
                &[100, 100]
            ),
            2
        );
        assert_eq!(
            list.hsps.iter().map(|hsp| hsp.score).collect::<Vec<_>>(),
            vec![100, 70]
        );
    }

    #[test]
    fn test_filter_blast_results_applies_subject_best_hit_branch() {
        let mut results = blast_hsp_results_new(1);
        let mut hitlist = blast_hit_list_new(10);
        let mut list = HspList::new(11);
        let mut best = make_hsp(100, 1e-30);
        best.query_offset = 10;
        best.query_end = 50;
        best.context = 0;
        let mut contained = make_hsp(80, 1e-20);
        contained.query_offset = 12;
        contained.query_end = 48;
        contained.context = 0;
        list.hsps = vec![best, contained];
        list.best_evalue = s_blast_get_best_evalue(&list);
        hitlist.hsp_lists.push(list);
        results.hitlists[0] = Some(hitlist);

        let mut filter_opts = crate::hspfilter_culling::blast_hspfiltering_options_new();
        filter_opts.subject_besthit_opts =
            Some(crate::hspfilter_culling::BlastHSPSubjectBestHitOptions { max_range_diff: 3 });
        let hit_options = HitSavingOptions {
            program_number: crate::program::BLASTP,
            hsp_filt_opt: Some(filter_opts),
            ..HitSavingOptions::default()
        };

        s_filter_blast_results(&mut results, &hit_options, &[100], None);
        let hitlist = results.hitlists[0].as_ref().expect("hitlist");
        assert_eq!(hitlist.hsp_lists[0].hsps.len(), 1);
        assert_eq!(hitlist.hsp_lists[0].hsps[0].score, 100);
    }

    #[test]
    fn test_filter_blast_results_applies_best_hit_pipe() {
        let mut results = blast_hsp_results_new(1);
        let mut hitlist = blast_hit_list_new(10);
        let mut list = HspList::new(11);
        let mut best = make_hsp(100, 1e-30);
        best.query_offset = 0;
        best.query_end = 40;
        best.context = 0;
        let mut contained = make_hsp(40, 1e-20);
        contained.query_offset = 5;
        contained.query_end = 33;
        contained.context = 0;
        list.hsps = vec![best, contained];
        list.best_evalue = s_blast_get_best_evalue(&list);
        hitlist.hsp_lists.push(list);
        results.hitlists[0] = Some(hitlist);

        let mut filter_opts = crate::hspfilter_culling::blast_hspfiltering_options_new();
        filter_opts.best_hit = Some(crate::hspfilter_culling::blast_hspbest_hit_options_new(
            0.1, 0.1,
        ));
        let hit_options = HitSavingOptions {
            program_number: crate::program::BLASTP,
            hsp_filt_opt: Some(filter_opts),
            hitlist_size: 10,
            ..HitSavingOptions::default()
        };
        let query_info = crate::queryinfo::QueryInfo {
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
            min_length: 100,
        };

        s_filter_blast_results(&mut results, &hit_options, &[100], Some(&query_info));
        let hitlist = results.hitlists[0].as_ref().expect("hitlist");
        assert_eq!(hitlist.hsp_lists[0].hsps.len(), 1);
        assert_eq!(hitlist.hsp_lists[0].hsps[0].score, 100);
    }

    #[test]
    fn test_filter_blast_results_applies_culling_pipe_across_subjects() {
        let mut results = blast_hsp_results_new(1);
        let mut hitlist = blast_hit_list_new(10);

        let mut strong_list = HspList::new(11);
        let mut strong = make_hsp(100, 1e-30);
        strong.query_offset = 0;
        strong.query_end = 80;
        strong.subject_offset = 0;
        strong.subject_end = 80;
        strong.context = 0;
        strong_list.hsps.push(strong);
        strong_list.best_evalue = s_blast_get_best_evalue(&strong_list);

        let mut weak_list = HspList::new(12);
        let mut weak = make_hsp(40, 1e-20);
        weak.query_offset = 10;
        weak.query_end = 70;
        weak.subject_offset = 10;
        weak.subject_end = 70;
        weak.context = 0;
        weak_list.hsps.push(weak);
        weak_list.best_evalue = s_blast_get_best_evalue(&weak_list);

        hitlist.hsp_lists = vec![weak_list, strong_list];
        results.hitlists[0] = Some(hitlist);

        let mut filter_opts = crate::hspfilter_culling::blast_hspfiltering_options_new();
        filter_opts.culling_opts =
            Some(crate::hspfilter_culling::BlastHSPCullingOptions { max_hits: 1 });
        let hit_options = HitSavingOptions {
            program_number: crate::program::BLASTP,
            hsp_filt_opt: Some(filter_opts),
            hitlist_size: 10,
            ..HitSavingOptions::default()
        };
        let query_info = crate::queryinfo::QueryInfo {
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
            min_length: 100,
        };

        s_filter_blast_results(&mut results, &hit_options, &[100], Some(&query_info));

        let hitlist = results.hitlists[0].as_ref().expect("hitlist");
        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(hitlist.hsp_lists[0].oid, 11);
        assert_eq!(hitlist.hsp_lists[0].hsps.len(), 1);
        assert_eq!(hitlist.hsp_lists[0].hsps[0].score, 100);
    }

    #[test]
    fn test_filter_blast_results_skips_pipe_filters_without_query_info() {
        let mut results = blast_hsp_results_new(1);
        let mut hitlist = blast_hit_list_new(10);
        let mut list = HspList::new(11);
        let mut best = make_hsp(100, 1e-30);
        best.query_offset = 0;
        best.query_end = 40;
        best.context = 0;
        let mut contained = make_hsp(40, 1e-20);
        contained.query_offset = 5;
        contained.query_end = 33;
        contained.context = 0;
        list.hsps = vec![best, contained];
        list.best_evalue = s_blast_get_best_evalue(&list);
        hitlist.hsp_lists.push(list);
        results.hitlists[0] = Some(hitlist);

        let mut filter_opts = crate::hspfilter_culling::blast_hspfiltering_options_new();
        filter_opts.best_hit = Some(crate::hspfilter_culling::blast_hspbest_hit_options_new(
            0.1, 0.1,
        ));
        filter_opts.culling_opts =
            Some(crate::hspfilter_culling::BlastHSPCullingOptions { max_hits: 1 });
        let hit_options = HitSavingOptions {
            program_number: crate::program::BLASTP,
            hsp_filt_opt: Some(filter_opts),
            hitlist_size: 10,
            ..HitSavingOptions::default()
        };

        s_filter_blast_results(&mut results, &hit_options, &[100], None);

        let hitlist = results.hitlists[0].as_ref().expect("hitlist");
        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(
            hitlist.hsp_lists[0]
                .hsps
                .iter()
                .map(|hsp| hsp.score)
                .collect::<Vec<_>>(),
            vec![100, 40]
        );
    }

    #[test]
    fn test_direct_max_hsp_and_query_coverage_filter_edges() {
        let mut list = blast_hsp_list_new(0);
        list.hsps = vec![
            make_hsp(100, 1e-20),
            make_hsp(90, 1e-10),
            make_hsp(80, 1e-5),
        ];
        list.best_evalue = s_blast_get_best_evalue(&list);

        let negative_max = HitSavingOptions {
            max_hsps_per_subject: -1,
            ..HitSavingOptions::default()
        };
        assert_eq!(
            blast_trim_hsp_list_by_max_hsps(Some(&mut list), &negative_max),
            0
        );
        assert_eq!(list.hsps.len(), 3);
        assert_eq!(blast_trim_hsp_list_by_max_hsps(None, &negative_max), 0);

        let max_two = HitSavingOptions {
            max_hsps_per_subject: 2,
            ..HitSavingOptions::default()
        };
        assert_eq!(
            blast_trim_hsp_list_by_max_hsps(Some(&mut list), &max_two),
            0
        );
        assert_eq!(
            list.hsps.iter().map(|hsp| hsp.score).collect::<Vec<_>>(),
            vec![100, 90]
        );
        assert_eq!(list.best_evalue, 1e-20);

        let mut stale_best = blast_hsp_list_new(0);
        stale_best.hsps = vec![
            make_hsp(100, 1e-10),
            make_hsp(90, 1e-20),
            make_hsp(80, 1e-5),
        ];
        stale_best.best_evalue = s_blast_get_best_evalue(&stale_best);
        let max_one = HitSavingOptions {
            max_hsps_per_subject: 1,
            ..HitSavingOptions::default()
        };
        assert_eq!(
            blast_trim_hsp_list_by_max_hsps(Some(&mut stale_best), &max_one),
            0
        );
        assert_eq!(stale_best.hsps.len(), 1);
        assert_eq!(stale_best.hsps[0].evalue, 1e-10);
        assert_eq!(stale_best.best_evalue, 1e-20);

        let zero_cov = HitSavingOptions {
            query_cov_hsp_perc: 0.0,
            ..HitSavingOptions::default()
        };
        assert_eq!(
            blast_hsp_list_reap_by_query_coverage(Some(&mut list), &zero_cov, &[]),
            0
        );
        assert_eq!(list.hsps.len(), 2);

        list.hsps[0].context = 7;
        list.hsps[1].context = 8;
        let strict_cov = HitSavingOptions {
            query_cov_hsp_perc: 1.0,
            ..HitSavingOptions::default()
        };
        assert_eq!(
            blast_hsp_list_reap_by_query_coverage(Some(&mut list), &strict_cov, &[100]),
            0
        );
        assert!(list.hsps.is_empty());
        // NCBI `s_BlastGetBestEvalue` seeds with `(double)INT4_MAX`.
        assert_eq!(list.best_evalue, i32::MAX as f64);
        assert_eq!(
            blast_hsp_list_reap_by_query_coverage(None, &strict_cov, &[100]),
            0
        );
    }

    #[test]
    fn test_query_coverage_noop_preserves_cached_best_evalue() {
        let mut list = blast_hsp_list_new(0);
        let mut hsp = make_hsp(100, 1e-20);
        hsp.context = 0;
        hsp.query_offset = 0;
        hsp.query_end = 90;
        list.hsps = vec![hsp];
        list.best_evalue = 42.0;
        let permissive_cov = HitSavingOptions {
            query_cov_hsp_perc: 50.0,
            ..HitSavingOptions::default()
        };

        assert_eq!(
            blast_hsp_list_reap_by_query_coverage(Some(&mut list), &permissive_cov, &[100]),
            0
        );

        assert_eq!(list.hsps.len(), 1);
        assert_eq!(list.best_evalue, 42.0);
    }

    #[test]
    fn test_blast_hsp_list_reap_by_evalue_matches_c_threshold() {
        let mut list = blast_hsp_list_new(0);
        list.hsps = vec![
            make_hsp(10, 0.5),
            make_hsp(30, 1.0e-20),
            make_hsp(20, 2.0),
            make_hsp(40, 1.0),
        ];
        list.best_evalue = s_blast_get_best_evalue(&list);

        let hit_options = HitSavingOptions {
            expect_value: 1.0,
            ..HitSavingOptions::default()
        };
        assert_eq!(
            blast_hsp_list_reap_by_evalue(Some(&mut list), &hit_options),
            0
        );
        let remaining: Vec<(i32, f64)> = list
            .hsps
            .iter()
            .map(|hsp| (hsp.score, hsp.evalue))
            .collect();
        assert_eq!(remaining, vec![(10, 0.5), (30, 1.0e-20), (40, 1.0)]);
        assert_eq!(list.best_evalue, 1.0e-20);
        assert_eq!(blast_hsp_list_reap_by_evalue(None, &hit_options), 0);

        let strict = HitSavingOptions {
            expect_value: 0.0,
            ..HitSavingOptions::default()
        };
        assert_eq!(blast_hsp_list_reap_by_evalue(Some(&mut list), &strict), 0);
        assert!(list.hsps.is_empty());
        // NCBI `Blast_HSPListReapByEvalue` does NOT recompute `best_evalue`,
        // so it retains the pre-filter value of 1.0e-20.
        assert_eq!(list.best_evalue, 1.0e-20);
    }

    #[test]
    fn test_blast_hsp_list_purge_common_endpoints_matches_c_two_pass_rule() {
        let mut same_start_best_score = make_hsp(100, 1.0e-5);
        same_start_best_score.query_offset = 0;
        same_start_best_score.subject_offset = 0;
        same_start_best_score.query_end = 10;
        same_start_best_score.subject_end = 10;

        let mut same_start_worse_score = make_hsp(90, 1.0e-20);
        same_start_worse_score.query_offset = 0;
        same_start_worse_score.subject_offset = 0;
        same_start_worse_score.query_end = 12;
        same_start_worse_score.subject_end = 12;

        let mut same_end_lower_tiebreak = make_hsp(80, 1.0e-8);
        same_end_lower_tiebreak.query_offset = 2;
        same_end_lower_tiebreak.subject_offset = 2;
        same_end_lower_tiebreak.query_end = 20;
        same_end_lower_tiebreak.subject_end = 20;

        let mut same_end_higher_tiebreak = make_hsp(80, 1.0e-10);
        same_end_higher_tiebreak.query_offset = 4;
        same_end_higher_tiebreak.subject_offset = 4;
        same_end_higher_tiebreak.query_end = 20;
        same_end_higher_tiebreak.subject_end = 20;

        let mut same_end_worse_score = make_hsp(70, 1.0e-30);
        same_end_worse_score.query_offset = 3;
        same_end_worse_score.subject_offset = 3;
        same_end_worse_score.query_end = 20;
        same_end_worse_score.subject_end = 20;

        let mut list = blast_hsp_list_new(0);
        list.hsps = vec![
            same_end_worse_score,
            same_start_worse_score,
            same_end_lower_tiebreak,
            same_start_best_score,
            same_end_higher_tiebreak,
        ];
        list.best_evalue = s_blast_get_best_evalue(&list);

        assert_eq!(
            blast_hsp_list_purge_hsps_with_common_endpoints(
                crate::program::BLASTP,
                Some(&mut list),
                true,
            ),
            2
        );
        let remaining: Vec<(i32, i32, i32, i32, i32)> = list
            .hsps
            .iter()
            .map(|hsp| {
                (
                    hsp.score,
                    hsp.query_offset,
                    hsp.subject_offset,
                    hsp.query_end,
                    hsp.subject_end,
                )
            })
            .collect();
        assert_eq!(remaining, vec![(100, 0, 0, 10, 10), (80, 4, 4, 20, 20)]);
        assert_eq!(list.best_evalue, 1.0e-10);
        assert_eq!(
            blast_hsp_list_purge_hsps_with_common_endpoints(crate::program::BLASTP, None, true),
            0
        );
    }

    #[test]
    fn test_purge_common_endpoints_phi_blast_is_noop() {
        let mut first = make_hsp(100, 1.0e-5);
        first.query_offset = 0;
        first.subject_offset = 0;
        first.query_end = 10;
        first.subject_end = 10;

        let mut duplicate = make_hsp(90, 1.0e-20);
        duplicate.query_offset = 0;
        duplicate.subject_offset = 0;
        duplicate.query_end = 12;
        duplicate.subject_end = 12;

        let mut list = blast_hsp_list_new(0);
        list.hsps = vec![duplicate, first];
        list.best_evalue = s_blast_get_best_evalue(&list);

        assert_eq!(
            blast_hsp_list_purge_hsps_with_common_endpoints(
                crate::program::PHI_BLASTN,
                Some(&mut list),
                true,
            ),
            2
        );
        assert_eq!(list.hsps.len(), 2);
        assert_eq!(list.best_evalue, 1.0e-20);
    }

    #[test]
    fn test_purge_common_endpoints_blastn_can_trim_duplicate_start() {
        let mut leader = make_hsp(100, 1.0e-5);
        leader.query_offset = 0;
        leader.subject_offset = 0;
        leader.query_end = 4;
        leader.subject_end = 4;
        leader.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![(
            GapAlignOpType::Sub,
            4,
        )]));

        let mut duplicate = make_hsp(90, 1.0e-20);
        duplicate.query_offset = 0;
        duplicate.subject_offset = 0;
        duplicate.query_end = 8;
        duplicate.subject_end = 8;
        duplicate.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![(
            GapAlignOpType::Sub,
            8,
        )]));

        let mut list = blast_hsp_list_new(0);
        list.hsps = vec![duplicate, leader];
        list.best_evalue = s_blast_get_best_evalue(&list);

        let returned_count = blast_hsp_list_purge_hsps_with_common_endpoints(
            crate::program::BLASTN,
            Some(&mut list),
            false,
        );
        // C decrements its local `hsp_count` on the BLASTN trim path even
        // though the trimmed HSP remains in `hsp_array` beyond that sorted
        // prefix. Rust keeps the owned vector compact, but mirrors the return.
        assert_eq!(returned_count, 1);
        assert_eq!(list.hsps.len(), 2);
        let trimmed = list
            .hsps
            .iter()
            .find(|hsp| hsp.score == 90)
            .expect("trimmed hsp");
        assert_eq!((trimmed.query_offset, trimmed.subject_offset), (4, 4));
        assert_eq!((trimmed.query_end, trimmed.subject_end), (8, 8));
        assert_eq!(
            trimmed.edit_script.as_ref().unwrap().ops_vec(),
            vec![(GapAlignOpType::Sub, 4)]
        );
        assert_eq!(list.best_evalue, 1.0e-20);
    }

    #[test]
    fn test_purge_common_endpoints_blastn_excludes_first_pass_extras_from_second_pass() {
        let mut start_leader = make_hsp(100, 1.0e-5);
        start_leader.query_offset = 0;
        start_leader.subject_offset = 0;
        start_leader.query_end = 4;
        start_leader.subject_end = 4;
        start_leader.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![(
            GapAlignOpType::Sub,
            4,
        )]));

        let mut trimmed_extra = make_hsp(90, 1.0e-20);
        trimmed_extra.query_offset = 0;
        trimmed_extra.subject_offset = 0;
        trimmed_extra.query_end = 8;
        trimmed_extra.subject_end = 8;
        trimmed_extra.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![(
            GapAlignOpType::Sub,
            8,
        )]));

        let mut same_end_active = make_hsp(95, 1.0e-10);
        same_end_active.query_offset = 2;
        same_end_active.subject_offset = 2;
        same_end_active.query_end = 8;
        same_end_active.subject_end = 8;
        same_end_active.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![(
            GapAlignOpType::Sub,
            6,
        )]));

        let mut list = blast_hsp_list_new(0);
        list.hsps = vec![trimmed_extra, same_end_active, start_leader];
        list.best_evalue = s_blast_get_best_evalue(&list);

        let returned_count = blast_hsp_list_purge_hsps_with_common_endpoints(
            crate::program::BLASTN,
            Some(&mut list),
            false,
        );

        assert_eq!(returned_count, 2);
        assert_eq!(list.hsps.len(), 3);
        assert_eq!(
            list.hsps
                .iter()
                .map(|hsp| (
                    hsp.score,
                    hsp.query_offset,
                    hsp.subject_offset,
                    hsp.query_end,
                    hsp.subject_end,
                ))
                .collect::<Vec<_>>(),
            vec![(100, 0, 0, 4, 4), (95, 2, 2, 8, 8), (90, 4, 4, 8, 8)]
        );
    }

    #[test]
    fn blast_hsp_list_adapter_uses_hspcnt_not_allocated_slots() {
        let mut active = BlastHSP::from_legacy_hsp(make_hsp(50, 1.0e-20));
        active.num = 7;
        active.gap_info = Some(crate::gapinfo::GapEditScript::from_ops(vec![
            (GapAlignOpType::Sub, 2),
            (GapAlignOpType::Del2, 1),
            (GapAlignOpType::Ins1, 1),
            (GapAlignOpType::Del, 3),
            (GapAlignOpType::Ins, 2),
        ]));
        let mut inactive = BlastHSP::from_legacy_hsp(make_hsp(90, 1.0e-40));
        inactive.num = 99;

        let list = BlastHSPList {
            oid: 42,
            query_index: 0,
            hsp_array: vec![Some(active), Some(inactive)],
            hspcnt: 1,
            allocated: 2,
            hsp_max: 10,
            do_not_reallocate: false,
            best_evalue: 1.0e-20,
        };

        let legacy = list.into_legacy_hsp_list();

        assert_eq!(legacy.hsps.len(), 1);
        assert_eq!(legacy.hsps[0].score, 50);
        assert_eq!(legacy.hsps[0].num_gaps, 2);
    }

    #[test]
    fn blast_hit_list_adapter_uses_hsplist_count_not_allocated_slots() {
        let active = BlastHSPList::from_legacy_hsp_list({
            let mut list = HspList::new(11);
            list.add_hsp(make_hsp(50, 1.0e-20));
            list
        });
        let inactive = BlastHSPList::from_legacy_hsp_list({
            let mut list = HspList::new(99);
            list.add_hsp(make_hsp(90, 1.0e-40));
            list
        });

        let hit_list = BlastHitList {
            hsplist_count: 1,
            hsplist_max: 10,
            worst_evalue: 1.0e-20,
            low_score: 50,
            heapified: false,
            hsplist_array: vec![Some(active), Some(inactive)],
            hsplist_current: 2,
            num_hits: 0,
        };

        let legacy = hit_list.into_legacy_hit_list();

        assert_eq!(legacy.hsp_lists.len(), 1);
        assert_eq!(legacy.hsp_lists[0].oid, 11);
        assert_eq!(legacy.hsp_lists[0].hsps[0].score, 50);
    }

    #[test]
    fn blast_hsp_results_adapter_uses_num_queries_not_allocated_slots() {
        let active_hit_list = BlastHitList {
            hsplist_count: 0,
            hsplist_max: 10,
            worst_evalue: 0.0,
            low_score: i32::MAX,
            heapified: false,
            hsplist_array: Vec::new(),
            hsplist_current: 0,
            num_hits: 0,
        };
        let inactive_hit_list = BlastHitList {
            hsplist_count: 0,
            hsplist_max: 10,
            worst_evalue: 0.0,
            low_score: i32::MAX,
            heapified: false,
            hsplist_array: Vec::new(),
            hsplist_current: 0,
            num_hits: 0,
        };

        let results = BlastHSPResults {
            num_queries: 1,
            hitlist_array: vec![Some(active_hit_list), Some(inactive_hit_list)],
        };

        let legacy = results.into_legacy_hsp_results();

        assert_eq!(legacy.hitlists.len(), 1);
        assert!(legacy.hitlists[0].is_some());
    }

    #[test]
    fn blast_hsp_list_purge_ignores_inactive_allocated_slots_and_preserves_num() {
        let mut active = BlastHSP::from_legacy_hsp(make_hsp(50, 1.0e-20));
        active.query.offset = 0;
        active.subject.offset = 0;
        active.query.end = 10;
        active.subject.end = 10;
        active.subject.frame = 1;
        active.num = 7;

        let mut inactive_duplicate = active.clone();
        inactive_duplicate.score = 90;
        inactive_duplicate.evalue = 1.0e-40;
        inactive_duplicate.num = 99;

        let mut list = BlastHSPList {
            oid: 42,
            query_index: 0,
            hsp_array: vec![Some(active), Some(inactive_duplicate)],
            hspcnt: 1,
            allocated: 2,
            hsp_max: 10,
            do_not_reallocate: false,
            best_evalue: 1.0e-20,
        };

        assert_eq!(
            blast_hsp_list_purge_blast_hsps_with_common_endpoints(
                crate::program::TBLASTN,
                Some(&mut list),
                true,
            ),
            1
        );

        assert_eq!(list.hspcnt, 1);
        assert_eq!(list.hsp_array.len(), 1);
        let hsp = list.hsp_array[0].as_ref().expect("active hsp");
        assert_eq!(hsp.score, 50);
        assert_eq!(hsp.num, 7);
        assert_eq!(list.best_evalue, 1.0e-20);
    }

    #[test]
    fn blast_compute_traceback_drains_stream_runs_callback_and_collects_results() {
        let stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[30]);
        let mut list = HspList::new(42);
        let mut hsp = make_hsp(20, 1.0e-3);
        hsp.context = 0;
        list.add_hsp(hsp);
        assert_eq!(stream.blast_hspstream_write(0, list), 0);
        let mut progress = crate::util::s_blast_progress_new(None);
        let hit_options = HitSavingOptions {
            hitlist_size: 10,
            ..HitSavingOptions::default()
        };

        let (status, results) = blast_compute_traceback(
            crate::program::BLASTP,
            Some(&stream),
            Some(&query_info),
            &hit_options,
            30,
            true,
            |hsp_list| {
                hsp_list.hsps[0].score += 5;
                0
            },
            None,
            Some(&mut progress),
        );

        assert_eq!(status, 0);
        assert_eq!(progress.stage, crate::util::EBlastStage::TracebackSearch);
        let results = results.expect("results");
        let hitlist = results.hitlists[0].as_ref().expect("query hitlist");
        assert_eq!(hitlist.hsp_lists[0].oid, 42);
        assert_eq!(hitlist.hsp_lists[0].hsps[0].score, 25);
    }

    #[test]
    fn blast_compute_traceback_rejects_missing_required_inputs() {
        let stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[30]);

        let (status, results) = blast_compute_traceback(
            crate::program::BLASTP,
            None,
            Some(&query_info),
            &HitSavingOptions::default(),
            30,
            true,
            |_hsp_list| 0,
            None,
            None,
        );
        assert_eq!(status, crate::util::BLASTERR_INVALIDPARAM);
        assert!(results.is_none());

        let (status, results) = blast_compute_traceback(
            crate::program::BLASTP,
            Some(&stream),
            None,
            &HitSavingOptions::default(),
            30,
            true,
            |_hsp_list| 0,
            None,
            None,
        );
        assert_eq!(status, crate::util::BLASTERR_INVALIDPARAM);
        assert!(results.is_none());
        assert_eq!(stream.blast_hspstream_write(0, HspList::new(9)), 0);
    }

    #[test]
    fn blast_compute_traceback_accepts_empty_stream_without_callbacks() {
        let stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[30]);
        let mut callback_called = false;

        let (status, results) = blast_compute_traceback(
            crate::program::BLASTP,
            Some(&stream),
            Some(&query_info),
            &HitSavingOptions::default(),
            30,
            true,
            |_hsp_list| {
                callback_called = true;
                0
            },
            None,
            None,
        );

        assert_eq!(status, 0);
        assert!(!callback_called);
        let results = results.expect("empty traceback results");
        assert_eq!(results.hitlists.len(), 1);
        assert!(results.hitlists[0].is_none());
        assert_eq!(stream.blast_hspstream_write(0, HspList::new(9)), -1);
    }

    #[test]
    fn blast_compute_traceback_collates_same_oid_batch_by_query_context() {
        let stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 2, None);
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[30, 40]);

        let mut q0_list = HspList::new(17);
        let mut q0_hsp = make_hsp(20, 1.0e-3);
        q0_hsp.context = 0;
        q0_list.add_hsp(q0_hsp);
        assert_eq!(stream.blast_hspstream_write(0, q0_list), 0);

        let mut q1_list = HspList::new(17);
        let mut q1_hsp = make_hsp(30, 1.0e-4);
        q1_hsp.context = 1;
        q1_list.add_hsp(q1_hsp);
        assert_eq!(stream.blast_hspstream_write(1, q1_list), 0);

        let hit_options = HitSavingOptions {
            hitlist_size: 10,
            mask_level: 101,
            ..HitSavingOptions::default()
        };
        let mut seen_contexts = Vec::new();

        let (status, results) = blast_compute_traceback(
            crate::program::BLASTP,
            Some(&stream),
            Some(&query_info),
            &hit_options,
            40,
            true,
            |hsp_list| {
                seen_contexts.push(hsp_list.hsps[0].context);
                0
            },
            None,
            None,
        );

        assert_eq!(status, 0);
        assert_eq!(seen_contexts, vec![0, 1]);
        let results = results.expect("results");
        let q0 = results.hitlists[0].as_ref().expect("query 0 hitlist");
        let q1 = results.hitlists[1].as_ref().expect("query 1 hitlist");
        assert_eq!(q0.hsp_lists.len(), 1);
        assert_eq!(q1.hsp_lists.len(), 1);
        assert_eq!(q0.hsp_lists[0].oid, 17);
        assert_eq!(q0.hsp_lists[0].hsps[0].context, 0);
        assert_eq!(q1.hsp_lists[0].oid, 17);
        assert_eq!(q1.hsp_lists[0].hsps[0].context, 1);
    }

    #[test]
    fn blast_compute_traceback_skips_empty_lists_after_callback() {
        let stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[30]);
        let mut list = HspList::new(42);
        let mut hsp = make_hsp(20, 1.0e-3);
        hsp.context = 0;
        list.add_hsp(hsp);
        assert_eq!(stream.blast_hspstream_write(0, list), 0);
        let mut callback_called = false;

        let (status, results) = blast_compute_traceback(
            crate::program::BLASTP,
            Some(&stream),
            Some(&query_info),
            &HitSavingOptions::default(),
            30,
            true,
            |hsp_list| {
                callback_called = true;
                hsp_list.hsps.clear();
                0
            },
            None,
            None,
        );

        assert_eq!(status, 0);
        assert!(callback_called);
        let results = results.expect("results");
        assert!(results.hitlists[0].is_none());
        assert_eq!(stream.blast_hspstream_write(0, HspList::new(9)), -1);
    }

    #[test]
    fn blast_compute_traceback_skips_stream_empty_lists_without_callback() {
        let stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[30]);
        assert_eq!(stream.blast_hspstream_write(0, HspList::new(42)), 0);
        let mut callback_called = false;

        let (status, results) = blast_compute_traceback(
            crate::program::BLASTP,
            Some(&stream),
            Some(&query_info),
            &HitSavingOptions::default(),
            30,
            true,
            |hsp_list| {
                let _ = hsp_list;
                callback_called = true;
                0
            },
            None,
            None,
        );

        assert_eq!(status, 0);
        assert!(!callback_called);
        let results = results.expect("results");
        assert!(results.hitlists[0].is_none());
        assert_eq!(stream.blast_hspstream_write(0, HspList::new(9)), -1);
    }

    #[test]
    fn blast_compute_traceback_closes_stream_on_callback_error() {
        let stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[30]);
        let mut list = HspList::new(42);
        let mut hsp = make_hsp(20, 1.0e-3);
        hsp.context = 0;
        list.add_hsp(hsp);
        assert_eq!(stream.blast_hspstream_write(0, list), 0);

        let (status, results) = blast_compute_traceback(
            crate::program::BLASTP,
            Some(&stream),
            Some(&query_info),
            &HitSavingOptions::default(),
            30,
            true,
            |_hsp_list| 7,
            None,
            None,
        );

        assert_eq!(status, 7);
        assert!(results.is_none());
        assert_eq!(stream.blast_hspstream_write(0, HspList::new(9)), -1);
    }

    #[test]
    fn blast_compute_traceback_closes_stream_on_result_insert_error() {
        let stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);
        let query_info = crate::queryinfo::QueryInfo {
            num_queries: 0,
            contexts: Vec::new(),
            max_length: 0,
            min_length: 0,
        };
        let mut list = HspList::new(42);
        let mut hsp = make_hsp(20, 1.0e-3);
        hsp.context = 0;
        list.add_hsp(hsp);
        assert_eq!(stream.blast_hspstream_write(0, list), 0);

        let (status, results) = blast_compute_traceback(
            crate::program::BLASTP,
            Some(&stream),
            Some(&query_info),
            &HitSavingOptions::default(),
            30,
            true,
            |_hsp_list| 0,
            None,
            None,
        );

        assert_eq!(status, 1);
        assert!(results.is_none());
        assert_eq!(stream.blast_hspstream_write(0, HspList::new(9)), -1);
    }

    #[test]
    fn blast_compute_traceback_rejects_context_query_index_out_of_range() {
        let stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);
        let query_info = crate::queryinfo::QueryInfo {
            num_queries: 1,
            contexts: vec![crate::queryinfo::ContextInfo {
                query_offset: 0,
                query_length: 30,
                eff_searchsp: 0,
                length_adjustment: 0,
                query_index: 2,
                frame: 0,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: 30,
            min_length: 0,
        };
        let mut list = HspList::new(43);
        let mut hsp = make_hsp(20, 1.0e-3);
        hsp.context = 0;
        list.add_hsp(hsp);
        assert_eq!(stream.blast_hspstream_write(0, list), 0);

        let (status, results) = blast_compute_traceback(
            crate::program::BLASTP,
            Some(&stream),
            Some(&query_info),
            &HitSavingOptions::default(),
            30,
            true,
            |_hsp_list| 0,
            None,
            None,
        );

        assert_eq!(status, 1);
        assert!(results.is_none());
        assert_eq!(stream.blast_hspstream_write(0, HspList::new(9)), -1);
    }

    #[test]
    fn blast_compute_traceback_rejects_callback_updated_context_query_index() {
        let stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);
        let query_info = crate::queryinfo::QueryInfo {
            num_queries: 1,
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
                    query_offset: 30,
                    query_length: 30,
                    eff_searchsp: 0,
                    length_adjustment: 0,
                    query_index: 3,
                    frame: 0,
                    is_valid: true,
                    segment_flags: crate::queryinfo::E_NO_SEGMENTS,
                },
            ],
            max_length: 30,
            min_length: 0,
        };
        let mut list = HspList::new(43);
        let mut hsp = make_hsp(20, 1.0e-3);
        hsp.context = 0;
        list.add_hsp(hsp);
        assert_eq!(stream.blast_hspstream_write(0, list), 0);

        let (status, results) = blast_compute_traceback(
            crate::program::BLASTP,
            Some(&stream),
            Some(&query_info),
            &HitSavingOptions::default(),
            30,
            true,
            |hsp_list| {
                hsp_list.hsps[0].context = 1;
                0
            },
            None,
            None,
        );

        assert_eq!(status, 1);
        assert!(results.is_none());
        assert_eq!(stream.blast_hspstream_write(0, HspList::new(9)), -1);
    }

    #[test]
    fn blast_compute_traceback_interrupts_before_later_batch_callbacks() {
        let stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[30]);
        for oid in [1, 2] {
            let mut list = HspList::new(oid);
            let mut hsp = make_hsp(20 + oid, 1.0e-3);
            hsp.context = 0;
            list.add_hsp(hsp);
            assert_eq!(stream.blast_hspstream_write(0, list), 0);
        }
        let hit_options = HitSavingOptions {
            mask_level: 101,
            ..HitSavingOptions::default()
        };
        let mut progress = crate::util::s_blast_progress_new(None);
        let mut interrupt_calls = 0usize;
        let mut interrupt = |_progress: &crate::util::SBlastProgress| {
            interrupt_calls += 1;
            interrupt_calls == 2
        };
        let mut callback_oids = Vec::new();

        let (status, results) = blast_compute_traceback(
            crate::program::BLASTP,
            Some(&stream),
            Some(&query_info),
            &hit_options,
            30,
            true,
            |hsp_list| {
                callback_oids.push(hsp_list.oid);
                0
            },
            Some(&mut interrupt),
            Some(&mut progress),
        );

        assert_eq!(status, crate::diagnostics::BLASTERR_INTERRUPTED);
        assert!(results.is_none());
        assert_eq!(interrupt_calls, 2);
        assert_eq!(callback_oids, vec![1]);
        assert_eq!(progress.stage, crate::util::EBlastStage::TracebackSearch);
        assert_eq!(stream.blast_hspstream_write(0, HspList::new(9)), -1);
    }

    #[test]
    fn blast_compute_traceback_sorts_database_results_then_prunes_hitlist_size() {
        let stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[30]);
        for &(oid, evalue) in &[(7, 1.0e-3), (3, 1.0e-20)] {
            let mut list = HspList::new(oid);
            let mut hsp = make_hsp(20, evalue);
            hsp.context = 0;
            list.add_hsp(hsp);
            assert_eq!(stream.blast_hspstream_write(0, list), 0);
        }
        let hit_options = HitSavingOptions {
            hitlist_size: 1,
            ..HitSavingOptions::default()
        };

        let (status, results) = blast_compute_traceback(
            crate::program::BLASTP,
            Some(&stream),
            Some(&query_info),
            &hit_options,
            30,
            true,
            |_hsp_list| 0,
            None,
            None,
        );

        assert_eq!(status, 0);
        let results = results.expect("results");
        let hitlist = results.hitlists[0].as_ref().expect("query hitlist");
        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(hitlist.hsp_lists[0].oid, 3);
        assert_eq!(hitlist.hsp_lists[0].best_evalue, 1.0e-20);
    }

    #[test]
    fn blast_compute_traceback_skips_database_sort_for_non_database_search() {
        let stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[30]);
        for &(oid, evalue) in &[(3, 1.0e-2), (7, 1.0e-20)] {
            let mut list = HspList::new(oid);
            let mut hsp = make_hsp(20, evalue);
            hsp.context = 0;
            list.add_hsp(hsp);
            assert_eq!(stream.blast_hspstream_write(0, list), 0);
        }
        let hit_options = HitSavingOptions {
            hitlist_size: 2,
            mask_level: 101,
            ..HitSavingOptions::default()
        };

        let (status, results) = blast_compute_traceback(
            crate::program::BLASTP,
            Some(&stream),
            Some(&query_info),
            &hit_options,
            30,
            false,
            |_hsp_list| 0,
            None,
            None,
        );

        assert_eq!(status, 0);
        let results = results.expect("results");
        let hitlist = results.hitlists[0].as_ref().expect("query hitlist");
        assert_eq!(
            hitlist
                .hsp_lists
                .iter()
                .map(|list| (list.oid, list.best_evalue))
                .collect::<Vec<_>>(),
            vec![(3, 1.0e-2), (7, 1.0e-20)]
        );
    }

    #[test]
    fn blast_compute_traceback_applies_masklevel_after_result_collation() {
        let stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[100]);

        let mut strong_list = HspList::new(1);
        let mut strong = make_hsp(100, 1.0e-30);
        strong.context = 0;
        strong.query_offset = 0;
        strong.query_end = 100;
        strong.subject_offset = 0;
        strong.subject_end = 100;
        strong_list.add_hsp(strong);
        assert_eq!(stream.blast_hspstream_write(0, strong_list), 0);

        let mut weak_list = HspList::new(2);
        let mut weak = make_hsp(10, 1.0e-2);
        weak.context = 0;
        weak.query_offset = 10;
        weak.query_end = 90;
        weak.subject_offset = 10;
        weak.subject_end = 90;
        weak_list.add_hsp(weak);
        assert_eq!(stream.blast_hspstream_write(0, weak_list), 0);

        let hit_options = HitSavingOptions {
            hitlist_size: 10,
            mask_level: 100,
            ..HitSavingOptions::default()
        };

        let (status, results) = blast_compute_traceback(
            crate::program::BLASTP,
            Some(&stream),
            Some(&query_info),
            &hit_options,
            100,
            true,
            |_hsp_list| 0,
            None,
            None,
        );

        assert_eq!(status, 0);
        let results = results.expect("results");
        let hitlist = results.hitlists[0].as_ref().expect("query hitlist");
        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(hitlist.hsp_lists[0].oid, 1);
        assert_eq!(hitlist.hsp_lists[0].hsps[0].score, 100);
    }

    #[test]
    fn blast_compute_traceback_applies_max_hsps_and_query_coverage_filters() {
        let stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[100]);

        let mut kept_list = HspList::new(11);
        let mut high_cov = make_hsp(100, 1.0e-20);
        high_cov.context = 0;
        high_cov.query_offset = 0;
        high_cov.query_end = 80;
        kept_list.add_hsp(high_cov);
        let mut low_cov = make_hsp(90, 1.0e-10);
        low_cov.context = 0;
        low_cov.query_offset = 0;
        low_cov.query_end = 10;
        kept_list.add_hsp(low_cov);
        let mut trimmed_before_coverage = make_hsp(80, 1.0e-5);
        trimmed_before_coverage.context = 0;
        trimmed_before_coverage.query_offset = 0;
        trimmed_before_coverage.query_end = 90;
        kept_list.add_hsp(trimmed_before_coverage);
        assert_eq!(stream.blast_hspstream_write(0, kept_list), 0);

        let mut emptied_list = HspList::new(22);
        let mut weak_coverage = make_hsp(50, 1.0e-2);
        weak_coverage.context = 0;
        weak_coverage.query_offset = 0;
        weak_coverage.query_end = 5;
        emptied_list.add_hsp(weak_coverage);
        assert_eq!(stream.blast_hspstream_write(0, emptied_list), 0);

        let hit_options = HitSavingOptions {
            hitlist_size: 10,
            mask_level: 101,
            max_hsps_per_subject: 2,
            query_cov_hsp_perc: 50.0,
            ..HitSavingOptions::default()
        };

        let (status, results) = blast_compute_traceback(
            crate::program::BLASTP,
            Some(&stream),
            Some(&query_info),
            &hit_options,
            100,
            true,
            |_hsp_list| 0,
            None,
            None,
        );

        assert_eq!(status, 0);
        let results = results.expect("results");
        let hitlist = results.hitlists[0].as_ref().expect("query hitlist");
        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(hitlist.hsp_lists[0].oid, 11);
        assert_eq!(
            hitlist.hsp_lists[0]
                .hsps
                .iter()
                .map(|hsp| hsp.score)
                .collect::<Vec<_>>(),
            vec![100]
        );
    }

    #[test]
    fn blast_compute_traceback_applies_blastn_subject_best_hit_default_range() {
        let stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);
        let query_info = crate::queryinfo::QueryInfo::new_blastn(&[100]);

        let mut list = HspList::new(11);
        let mut forward = make_hsp(100, 1.0e-30);
        forward.context = 0;
        forward.query_frame = 1;
        forward.query_offset = 10;
        forward.query_end = 50;
        list.add_hsp(forward);

        let mut opposite_contained = make_hsp(80, 1.0e-20);
        opposite_contained.context = 1;
        opposite_contained.query_frame = -1;
        opposite_contained.query_offset = 48;
        opposite_contained.query_end = 90;
        list.add_hsp(opposite_contained);

        let mut opposite_outside = make_hsp(70, 1.0e-10);
        opposite_outside.context = 1;
        opposite_outside.query_frame = -1;
        opposite_outside.query_offset = 10;
        opposite_outside.query_end = 40;
        list.add_hsp(opposite_outside);
        assert_eq!(stream.blast_hspstream_write(0, list), 0);

        let mut filter_opts = crate::hspfilter_culling::blast_hspfiltering_options_new();
        filter_opts.subject_besthit_opts =
            Some(crate::hspfilter_culling::blast_hsp_subject_best_hit_options_new(false));
        let hit_options = HitSavingOptions {
            program_number: crate::program::BLASTN,
            hitlist_size: 10,
            mask_level: 101,
            hsp_filt_opt: Some(filter_opts),
            ..HitSavingOptions::default()
        };

        let (status, results) = blast_compute_traceback(
            crate::program::BLASTN,
            Some(&stream),
            Some(&query_info),
            &hit_options,
            100,
            true,
            |_hsp_list| 0,
            None,
            None,
        );

        assert_eq!(status, 0);
        let results = results.expect("results");
        let hitlist = results.hitlists[0].as_ref().expect("query hitlist");
        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(
            hitlist.hsp_lists[0]
                .hsps
                .iter()
                .map(|hsp| (hsp.score, hsp.context, hsp.query_offset, hsp.query_end))
                .collect::<Vec<_>>(),
            vec![(100, 0, 10, 50), (70, 1, 10, 40)]
        );
    }

    #[test]
    fn blast_compute_traceback_applies_culling_only_pipe() {
        let stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[100]);

        let mut weak_list = HspList::new(12);
        let mut weak = make_hsp(40, 1.0e-20);
        weak.context = 0;
        weak.query_offset = 10;
        weak.query_end = 70;
        weak.subject_offset = 10;
        weak.subject_end = 70;
        weak_list.add_hsp(weak);
        assert_eq!(stream.blast_hspstream_write(0, weak_list), 0);

        let mut strong_list = HspList::new(11);
        let mut strong = make_hsp(100, 1.0e-30);
        strong.context = 0;
        strong.query_offset = 0;
        strong.query_end = 80;
        strong.subject_offset = 0;
        strong.subject_end = 80;
        strong_list.add_hsp(strong);
        assert_eq!(stream.blast_hspstream_write(0, strong_list), 0);

        let mut filter_opts = crate::hspfilter_culling::blast_hspfiltering_options_new();
        filter_opts.culling_opts =
            Some(crate::hspfilter_culling::BlastHSPCullingOptions { max_hits: 1 });
        let hit_options = HitSavingOptions {
            hitlist_size: 10,
            mask_level: 101,
            hsp_filt_opt: Some(filter_opts),
            ..HitSavingOptions::default()
        };

        let (status, results) = blast_compute_traceback(
            crate::program::BLASTP,
            Some(&stream),
            Some(&query_info),
            &hit_options,
            100,
            true,
            |_hsp_list| 0,
            None,
            None,
        );

        assert_eq!(status, 0);
        let results = results.expect("results");
        let hitlist = results.hitlists[0].as_ref().expect("query hitlist");
        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(hitlist.hsp_lists[0].oid, 11);
        assert_eq!(hitlist.hsp_lists[0].hsps[0].score, 100);
    }

    #[test]
    fn blast_compute_traceback_query_coverage_reaps_missing_context_lengths() {
        let stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[100]);
        let mut list = HspList::new(44);
        let mut hsp = make_hsp(90, 1.0e-20);
        hsp.context = 7;
        hsp.query_offset = 0;
        hsp.query_end = 90;
        list.add_hsp(hsp);
        assert_eq!(stream.blast_hspstream_write(0, list), 0);
        let hit_options = HitSavingOptions {
            hitlist_size: 10,
            mask_level: 101,
            query_cov_hsp_perc: 1.0,
            ..HitSavingOptions::default()
        };

        let (status, results) = blast_compute_traceback(
            crate::program::BLASTP,
            Some(&stream),
            Some(&query_info),
            &hit_options,
            100,
            true,
            |_hsp_list| 0,
            None,
            None,
        );

        assert_eq!(status, 0);
        let results = results.expect("results");
        let hitlist = results.hitlists[0].as_ref().expect("query hitlist");
        assert!(hitlist.hsp_lists.is_empty());
    }

    #[test]
    fn blast_traceback_dispatch_kind_matches_c_branch_order() {
        assert_eq!(
            blast_traceback_dispatch_kind(crate::program::RPS_BLAST, 0, false),
            TracebackDispatchKind::Rps
        );
        assert_eq!(
            blast_traceback_dispatch_kind(crate::program::RPS_BLAST, 1, true),
            TracebackDispatchKind::Rps
        );
        assert_eq!(
            blast_traceback_dispatch_kind(crate::program::BLASTP, 1, false),
            TracebackDispatchKind::RedoAlignment
        );
        assert_eq!(
            blast_traceback_dispatch_kind(crate::program::BLASTP, 0, true),
            TracebackDispatchKind::RedoAlignment
        );
        assert_eq!(
            blast_traceback_dispatch_kind(crate::program::PHI_BLASTP, 0, false),
            TracebackDispatchKind::Phi
        );
        assert_eq!(
            blast_traceback_dispatch_kind(crate::program::PHI_BLASTP, 1, false),
            TracebackDispatchKind::RedoAlignment
        );
        assert_eq!(
            blast_traceback_dispatch_kind(crate::program::BLASTP, 0, false),
            TracebackDispatchKind::Ordinary
        );
    }

    #[test]
    fn blast_run_traceback_search_with_interrupt_closes_stream_and_honors_interrupt() {
        let stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[30]);
        let mut list = HspList::new(7);
        let mut hsp = make_hsp(20, 1.0e-3);
        hsp.context = 0;
        list.add_hsp(hsp);
        assert_eq!(stream.blast_hspstream_write(0, list), 0);
        let hit_options = HitSavingOptions::default();
        let mut progress = crate::util::s_blast_progress_new(None);
        let mut interrupt = |_progress: &crate::util::SBlastProgress| true;

        let (status, results) = blast_run_traceback_search_with_interrupt(
            crate::program::BLASTP,
            Some(&stream),
            Some(&query_info),
            &hit_options,
            30,
            true,
            false,
            |_hsp_list| 0,
            Some(&mut interrupt),
            Some(&mut progress),
            0,
        );

        assert_eq!(status, crate::diagnostics::BLASTERR_INTERRUPTED);
        assert!(results.is_none());
        assert_eq!(progress.stage, crate::util::EBlastStage::TracebackSearch);
        assert_eq!(stream.blast_hspstream_write(0, HspList::new(9)), -1);
    }

    #[test]
    fn blast_run_traceback_search_delegates_without_interrupt_state() {
        let stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[30]);
        let mut list = HspList::new(7);
        let mut hsp = make_hsp(20, 1.0e-3);
        hsp.context = 0;
        list.add_hsp(hsp);
        assert_eq!(stream.blast_hspstream_write(0, list), 0);
        let hit_options = HitSavingOptions {
            hitlist_size: 10,
            mask_level: 101,
            ..HitSavingOptions::default()
        };
        let mut callback_write_status = 0;

        let (status, results) = blast_run_traceback_search(
            crate::program::BLASTP,
            Some(&stream),
            Some(&query_info),
            &hit_options,
            30,
            true,
            false,
            |_hsp_list| {
                callback_write_status = stream.blast_hspstream_write(0, HspList::new(9));
                0
            },
            0,
        );

        assert_eq!(status, 0);
        let results = results.expect("results");
        let hitlist = results.hitlists[0].as_ref().expect("query hitlist");
        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(hitlist.hsp_lists[0].oid, 7);
        assert_eq!(callback_write_status, -1);
        assert_eq!(stream.blast_hspstream_write(0, HspList::new(10)), -1);
    }

    #[test]
    fn blast_run_traceback_search_with_interrupt_rejects_missing_stream() {
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[30]);

        let (status, results) = blast_run_traceback_search_with_interrupt(
            crate::program::BLASTP,
            None,
            Some(&query_info),
            &HitSavingOptions::default(),
            30,
            true,
            false,
            |_hsp_list| 0,
            None,
            None,
            0,
        );

        assert_eq!(status, crate::util::BLASTERR_INVALIDPARAM);
        assert!(results.is_none());
    }

    #[test]
    fn blast_run_traceback_search_with_interrupt_rejects_missing_query_before_close() {
        let stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);

        let (status, results) = blast_run_traceback_search_with_interrupt(
            crate::program::BLASTP,
            Some(&stream),
            None,
            &HitSavingOptions::default(),
            30,
            true,
            false,
            |_hsp_list| 0,
            None,
            None,
            0,
        );

        assert_eq!(status, crate::util::BLASTERR_INVALIDPARAM);
        assert!(results.is_none());
        assert_eq!(stream.blast_hspstream_write(0, HspList::new(9)), 0);
    }

    #[test]
    fn blast_run_traceback_search_with_interrupt_closes_ordinary_stream_before_callbacks() {
        let stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[30]);
        let mut list = HspList::new(7);
        let mut hsp = make_hsp(20, 1.0e-3);
        hsp.context = 0;
        list.add_hsp(hsp);
        assert_eq!(stream.blast_hspstream_write(0, list), 0);
        let hit_options = HitSavingOptions {
            hitlist_size: 10,
            mask_level: 101,
            ..HitSavingOptions::default()
        };
        let mut callback_write_status = 0;

        let (status, results) = blast_run_traceback_search_with_interrupt(
            crate::program::BLASTP,
            Some(&stream),
            Some(&query_info),
            &hit_options,
            30,
            true,
            false,
            |_hsp_list| {
                callback_write_status = stream.blast_hspstream_write(0, HspList::new(9));
                0
            },
            None,
            None,
            0,
        );

        assert_eq!(status, 0);
        assert!(results.is_some());
        assert_eq!(callback_write_status, -1);
        assert_eq!(stream.blast_hspstream_write(0, HspList::new(10)), -1);
    }

    #[test]
    fn blast_run_traceback_search_with_interrupt_uses_cbs_close_before_callbacks() {
        let stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[30]);
        for index in 0..1305 {
            let mut list = HspList::new(index);
            let evalue = if index <= 500 {
                1.0e-50
            } else {
                1.0e-5 + index as f64 * 1.0e-8
            };
            let mut hsp = make_hsp(100 - (index % 10), evalue);
            hsp.context = 0;
            list.add_hsp(hsp);
            assert_eq!(stream.blast_hspstream_write(0, list), 0);
        }
        let hit_options = HitSavingOptions {
            hitlist_size: 500,
            mask_level: 101,
            ..HitSavingOptions::default()
        };
        let mut callbacks = 0usize;

        let (status, results) = blast_run_traceback_search_with_interrupt(
            crate::program::BLASTP,
            Some(&stream),
            Some(&query_info),
            &hit_options,
            30,
            true,
            true,
            |_hsp_list| {
                callbacks += 1;
                0
            },
            None,
            None,
            0,
        );

        assert_eq!(status, 0);
        assert!(callbacks < 1305);
        assert_eq!(callbacks, 1100);
        let results = results.expect("results");
        let hitlist = results.hitlists[0].as_ref().expect("query hitlist");
        assert_eq!(hitlist.hsp_lists.len(), 500);
        assert_eq!(stream.blast_hspstream_write(0, HspList::new(9999)), -1);
    }

    fn make_rps_traceback_info() -> RpsTracebackInfo {
        let alphabet_size = 26;
        let num_rows = 5 + 1;
        RpsTracebackInfo {
            profile_header: RpsProfileHeader {
                magic_number: RPS_MAGIC_NUM,
                num_profiles: 2,
                start_offsets: vec![0, 2, 5],
                pssm_values: (0..(num_rows * alphabet_size)).map(|v| v as i32).collect(),
            },
            freq_ratios_header: Some(RpsFreqRatiosHeader {
                start_offsets: vec![0, 2, 5],
                freq_values: (0..(num_rows * alphabet_size))
                    .map(|v| 10_000 + v as i32)
                    .collect(),
            }),
            karlin_k: vec![0.1, 0.2],
        }
    }

    #[test]
    fn rps_gap_align_data_prepare_builds_concat_contexts_and_rows() {
        let rps_info = make_rps_traceback_info();

        let data = s_rps_gap_align_data_prepare(Some(&rps_info)).expect("RPS data");

        assert_eq!(data.alphabet_size, 26);
        assert!(data.position_based);
        assert_eq!(data.concat_db_info.num_queries, 2);
        assert_eq!(data.concat_db_info.contexts[0].query_offset, 0);
        assert_eq!(data.concat_db_info.contexts[0].query_length, 2);
        assert_eq!(data.concat_db_info.contexts[1].query_offset, 2);
        assert_eq!(data.concat_db_info.contexts[1].query_length, 3);
        assert_eq!(data.rps_pssm.len(), 6);
        assert_eq!(data.rps_pssm[2][0], 52);
        let freqs = data.rps_freq.as_ref().expect("freq rows");
        assert_eq!(freqs[5][25], 10_155);
        let ratios = data.rps_freq_ratios.as_ref().expect("freq ratio rows");
        // FREQ_RATIO_SCALE = 1e9 (NCBI blast_rps.h:83): 10155 / 1e9.
        assert!((ratios[5][25] - 1.0155e-5).abs() < 1e-12);
        assert!(s_rps_gap_align_data_prepare(None).is_err());
    }

    #[test]
    fn rps_gap_align_data_prepare_rejects_malformed_owned_profile_views() {
        let mut short_offsets = make_rps_traceback_info();
        short_offsets.profile_header.start_offsets = vec![0, 2];
        assert!(s_rps_gap_align_data_prepare(Some(&short_offsets)).is_err());

        let mut negative_profile_count = make_rps_traceback_info();
        negative_profile_count.profile_header.num_profiles = -1;
        negative_profile_count.profile_header.start_offsets = vec![0];
        negative_profile_count.profile_header.pssm_values = vec![0; 26];
        assert!(s_rps_gap_align_data_prepare(Some(&negative_profile_count)).is_err());

        let mut unknown_magic = make_rps_traceback_info();
        unknown_magic.profile_header.magic_number = 0x1234_5678;
        assert!(s_rps_gap_align_data_prepare(Some(&unknown_magic)).is_err());

        let mut negative_final_offset = make_rps_traceback_info();
        negative_final_offset.profile_header.start_offsets = vec![0, 2, -1];
        assert!(s_rps_gap_align_data_prepare(Some(&negative_final_offset)).is_err());

        let mut negative_profile_start = make_rps_traceback_info();
        negative_profile_start.profile_header.start_offsets = vec![-1, 2, 5];
        assert!(s_rps_gap_align_data_prepare(Some(&negative_profile_start)).is_err());

        let mut decreasing_offsets = make_rps_traceback_info();
        decreasing_offsets.profile_header.start_offsets = vec![0, 5, 2];
        assert!(s_rps_gap_align_data_prepare(Some(&decreasing_offsets)).is_err());

        let mut truncated_pssm = make_rps_traceback_info();
        truncated_pssm.profile_header.pssm_values.truncate(26);
        assert!(s_rps_gap_align_data_prepare(Some(&truncated_pssm)).is_err());

        let mut short_freq_offsets = make_rps_traceback_info();
        short_freq_offsets
            .freq_ratios_header
            .as_mut()
            .expect("freq header")
            .start_offsets = vec![0, 2];
        assert!(s_rps_gap_align_data_prepare(Some(&short_freq_offsets)).is_err());

        let mut decreasing_freq_offsets = make_rps_traceback_info();
        decreasing_freq_offsets
            .freq_ratios_header
            .as_mut()
            .expect("freq header")
            .start_offsets = vec![0, 5, 2];
        assert!(s_rps_gap_align_data_prepare(Some(&decreasing_freq_offsets)).is_err());

        let mut mismatched_freq_offsets = make_rps_traceback_info();
        mismatched_freq_offsets
            .freq_ratios_header
            .as_mut()
            .expect("freq header")
            .start_offsets = vec![0, 1, 5];
        assert!(s_rps_gap_align_data_prepare(Some(&mismatched_freq_offsets)).is_err());

        let mut truncated_freq = make_rps_traceback_info();
        truncated_freq
            .freq_ratios_header
            .as_mut()
            .expect("freq header")
            .freq_values
            .truncate(26);
        assert!(s_rps_gap_align_data_prepare(Some(&truncated_freq)).is_err());
    }

    #[test]
    fn rps_fill_freq_ratios_in_psi_matrix_scales_and_zero_pads() {
        let freq = vec![
            vec![1_000_000_000, 2_000_000_000, -500_000_000],
            vec![0, 1, 2],
        ];

        let ratios = s_rps_fill_freq_ratios_in_psi_matrix(&freq, 1, 3).expect("ratios");

        assert_eq!(ratios.len(), 1);
        assert_eq!(ratios[0][0], 1.0);
        assert_eq!(ratios[0][1], 2.0);
        assert_eq!(ratios[0][2], -0.5);
        assert!(ratios[0][3..].iter().all(|&value| value == 0.0));
        assert!(s_rps_fill_freq_ratios_in_psi_matrix(&freq, 3, 3).is_none());
        assert!(s_rps_fill_freq_ratios_in_psi_matrix(&freq, 1, 29).is_none());
    }

    #[test]
    fn rps_aux_info_from_bytes_parses_scoring_metadata_and_karlin_values() {
        let aux = b"BLOSUM62\n11 1\n0.041 0.14\n128 4096\n100.0\n8 0.21\n16 0.25\n";
        let info = rps_aux_info_from_bytes(aux).expect("RPS aux metadata");

        assert_eq!(info.matrix, "BLOSUM62");
        assert_eq!(info.gap_open, 11);
        assert_eq!(info.gap_extend, 1);
        assert!((info.ungapped_k - 0.041).abs() < 1e-12);
        assert!((info.ungapped_h - 0.14).abs() < 1e-12);
        assert_eq!(info.max_db_seq_length, 128);
        assert_eq!(info.db_length, 4096);
        assert!((info.scale_factor - 100.0).abs() < 1e-12);
        assert_eq!(info.profile_lengths, vec![8, 16]);
        assert_eq!(info.karlin_k, vec![0.21, 0.25]);

        assert!(rps_aux_info_from_bytes(b"BLOSUM62 11 1 0.1 0.2 8 16 100").is_err());
        assert!(rps_aux_info_from_bytes(b"BLOSUM62 11 1 0.1 0.2 8 16 100 4").is_err());
    }

    #[test]
    fn rps_gap_align_data_prepare_accepts_absent_frequency_ratios() {
        let mut rps_info = make_rps_traceback_info();
        rps_info.freq_ratios_header = None;

        let data = s_rps_gap_align_data_prepare(Some(&rps_info)).expect("RPS data");

        assert!(data.rps_freq.is_none());
        assert!(data.rps_freq_ratios.is_none());
        assert_eq!(data.rps_pssm.len(), 6);
        assert_eq!(data.concat_db_info.contexts[1].query_offset, 2);
        assert_eq!(data.concat_db_info.contexts[1].query_length, 3);
    }

    #[test]
    fn rps_fetch_consensus_sequence_uses_traceback_encoding_and_skips_missing_oid() {
        let seq_src = RecordingSeqSrc {
            seqs: vec![vec![1, 2, 3], vec![4, 5]],
            encodings: TestMutex::new(Vec::new()),
        };

        let data = s_rps_fetch_consensus_sequence(crate::program::RPS_BLAST, Some(&seq_src), 1)
            .expect("sequence");
        assert_eq!(data.sequence, vec![4, 5]);
        assert_eq!(data.length, 2);
        assert_eq!(
            seq_src.encodings.lock().unwrap().as_slice(),
            &[crate::seqsrc::SeqEncoding::Protein]
        );

        let data = s_rps_fetch_consensus_sequence(crate::program::RPS_TBLASTN, Some(&seq_src), 0)
            .expect("sequence");
        assert_eq!(data.sequence, vec![1, 2, 3]);
        assert_eq!(data.length, 3);

        assert!(
            s_rps_fetch_consensus_sequence(crate::program::RPS_TBLASTN, Some(&seq_src), 9)
                .is_none()
        );
        assert_eq!(
            seq_src.encodings.lock().unwrap().as_slice(),
            &[
                crate::seqsrc::SeqEncoding::Protein,
                crate::seqsrc::SeqEncoding::Protein,
                crate::seqsrc::SeqEncoding::Protein
            ]
        );
        assert!(s_rps_fetch_consensus_sequence(crate::program::RPS_BLAST, None, 0).is_none());
    }

    #[test]
    fn rps_profile_pssm_for_traceback_allocates_cbs_or_rescales_profile_rows() {
        let query = [
            crate::encoding::NCBISTDAA_A,
            crate::encoding::NCBISTDAA_A,
            crate::encoding::NCBISTDAA_C,
        ];
        let profile_rows: Vec<Vec<i32>> = crate::matrix::BLOSUM62
            .iter()
            .map(|row| row.to_vec())
            .collect();
        let profile_len = profile_rows.len();
        let gap_data = RpsGapAlignData {
            concat_db_info: crate::queryinfo::QueryInfo {
                num_queries: 1,
                contexts: vec![crate::queryinfo::ContextInfo {
                    query_offset: 0,
                    query_length: profile_len as i32,
                    eff_searchsp: 0,
                    length_adjustment: 0,
                    query_index: 0,
                    frame: 0,
                    is_valid: true,
                    segment_flags: crate::queryinfo::E_NO_SEGMENTS,
                }],
                max_length: profile_len as u32,
                min_length: 0,
            },
            rps_pssm: profile_rows,
            rps_freq: None,
            rps_freq_ratios: None,
            alphabet_size: crate::encoding::BLASTAA_SIZE,
            position_based: true,
        };
        let mut sbp =
            crate::stat::blast_score_blk_new(crate::encoding::BLASTAA_SEQ_CODE, 1).expect("sbp");
        sbp.name = Some("BLOSUM62".to_string());

        let scratch = s_rps_profile_pssm_for_traceback(
            crate::program::RPS_BLAST,
            &gap_data,
            0,
            &query,
            profile_len,
            2.0,
            &sbp,
            true,
        )
        .expect("cbs scratch");
        assert_eq!(
            scratch,
            vec![vec![0; crate::encoding::BLASTAA_SIZE]; profile_len]
        );

        let rescaled = s_rps_profile_pssm_for_traceback(
            crate::program::RPS_BLAST,
            &gap_data,
            0,
            &query,
            profile_len,
            2.0,
            &sbp,
            false,
        )
        .expect("rescaled pssm");
        assert_eq!(rescaled.len(), profile_len);
        assert_eq!(rescaled[0].len(), crate::encoding::BLASTAA_SIZE);
        assert_ne!(
            rescaled[0][crate::encoding::NCBISTDAA_A as usize],
            gap_data.rps_pssm[0][crate::encoding::NCBISTDAA_A as usize]
        );
        let direct = s_rps_profile_pssm_for_traceback(
            crate::program::RPS_TBLASTN,
            &gap_data,
            0,
            &query,
            profile_len,
            2.0,
            &sbp,
            false,
        )
        .expect("direct rps-tblastn pssm");
        assert_eq!(direct, gap_data.rps_pssm[..profile_len]);
        assert!(s_rps_profile_pssm_for_traceback(
            crate::program::RPS_BLAST,
            &gap_data,
            1,
            &query,
            profile_len,
            2.0,
            &sbp,
            false
        )
        .is_none());
        assert!(s_rps_profile_pssm_for_traceback(
            crate::program::RPS_BLAST,
            &gap_data,
            0,
            &query,
            profile_len + 1,
            2.0,
            &sbp,
            false
        )
        .is_none());

        let mut short_profile_row = gap_data.clone();
        short_profile_row.rps_pssm[0].truncate(crate::encoding::BLASTAA_SIZE - 1);
        assert!(s_rps_profile_pssm_for_traceback(
            crate::program::RPS_BLAST,
            &short_profile_row,
            0,
            &query,
            profile_len,
            2.0,
            &sbp,
            false
        )
        .is_none());
        assert!(s_rps_profile_pssm_for_traceback(
            crate::program::RPS_TBLASTN,
            &short_profile_row,
            0,
            &query,
            profile_len,
            2.0,
            &sbp,
            false
        )
        .is_none());

        let mut invalid_profile_context = gap_data.clone();
        invalid_profile_context.concat_db_info.contexts[0].is_valid = false;
        assert!(s_rps_profile_pssm_for_traceback(
            crate::program::RPS_BLAST,
            &invalid_profile_context,
            0,
            &query,
            profile_len,
            2.0,
            &sbp,
            true
        )
        .is_none());
    }

    #[test]
    fn rps_gap_align_data_prepare_uses_28_column_profile_format() {
        let alphabet_size = 28;
        let num_rows = 3 + 1;
        let rps_info = RpsTracebackInfo {
            profile_header: RpsProfileHeader {
                magic_number: RPS_MAGIC_NUM_28,
                num_profiles: 1,
                start_offsets: vec![0, 3],
                pssm_values: (0..(num_rows * alphabet_size)).map(|v| v as i32).collect(),
            },
            freq_ratios_header: None,
            karlin_k: vec![0.3],
        };

        let data = s_rps_gap_align_data_prepare(Some(&rps_info)).expect("RPS data");

        assert_eq!(data.alphabet_size, 28);
        assert_eq!(data.rps_pssm.len(), 4);
        assert_eq!(data.rps_pssm[3][27], 111);
        assert!(data.rps_freq.is_none());
    }

    #[test]
    fn rps_gap_align_data_prepare_marks_zero_length_profile_context_invalid() {
        let alphabet_size = 26;
        let num_rows = 2 + 1;
        let rps_info = RpsTracebackInfo {
            profile_header: RpsProfileHeader {
                magic_number: RPS_MAGIC_NUM,
                num_profiles: 2,
                start_offsets: vec![0, 0, 2],
                pssm_values: (0..(num_rows * alphabet_size)).map(|v| v as i32).collect(),
            },
            freq_ratios_header: None,
            karlin_k: vec![0.1, 0.2],
        };

        let data = s_rps_gap_align_data_prepare(Some(&rps_info)).expect("RPS data");

        assert_eq!(data.concat_db_info.num_queries, 2);
        assert_eq!(data.concat_db_info.max_length, 2);
        assert_eq!(data.concat_db_info.contexts[0].query_offset, 0);
        assert_eq!(data.concat_db_info.contexts[0].query_length, 0);
        assert!(!data.concat_db_info.contexts[0].is_valid);
        assert_eq!(data.concat_db_info.contexts[1].query_offset, 0);
        assert_eq!(data.concat_db_info.contexts[1].query_length, 2);
        assert!(data.concat_db_info.contexts[1].is_valid);
    }

    #[test]
    fn rps_compute_traceback_drains_stream_updates_orientation_and_inserts_results() {
        let stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[30]);
        let rps_info = make_rps_traceback_info();
        let mut list = HspList::new(1);
        let mut hsp = make_hsp(20, 1.0e-3);
        hsp.query_offset = 2;
        hsp.query_end = 7;
        hsp.query_gapped_start = 2;
        hsp.subject_offset = 10;
        hsp.subject_end = 16;
        hsp.subject_gapped_start = 11;
        hsp.edit_script = Some(crate::gapinfo::GapEditScript::from_ops(vec![
            (GapAlignOpType::Ins, 2),
            (GapAlignOpType::Del, 1),
        ]));
        list.add_hsp(hsp);
        assert_eq!(stream.blast_hspstream_write(0, list), 0);
        let mut results = HspResults::new(1);
        let hit_options = HitSavingOptions {
            hitlist_size: 10,
            ..HitSavingOptions::default()
        };
        let mut progress = crate::util::s_blast_progress_new(None);

        let status = s_rps_compute_traceback(
            crate::program::RPS_BLAST,
            Some(&stream),
            Some(&query_info),
            Some(&rps_info),
            &hit_options,
            Some(&mut results),
            |hsp_list, gap_data, profile_index, karlin_k| {
                assert_eq!(profile_index, 1);
                assert_eq!(gap_data.concat_db_info.contexts[1].query_offset, 2);
                assert_eq!(gap_data.rps_pssm[2][0], 52);
                let ratios = gap_data.rps_freq_ratios.as_ref().expect("freq ratios");
                // FREQ_RATIO_SCALE = 1e9 (NCBI blast_rps.h:83): 10052 / 1e9.
                assert!((ratios[2][0] - 1.0052e-5).abs() < 1e-12);
                assert!((karlin_k.expect("K") - 0.24).abs() < 1e-12);
                hsp_list.hsps[0].score += 5;
                0
            },
            None,
            Some(&mut progress),
        );

        assert_eq!(status, 0);
        assert_eq!(progress.stage, crate::util::EBlastStage::TracebackSearch);
        let hitlist = results.hitlists[0].as_ref().expect("query hitlist");
        assert_eq!(hitlist.hsp_lists[0].oid, 1);
        let hsp = &hitlist.hsp_lists[0].hsps[0];
        assert_eq!(hsp.score, 25);
        assert_eq!(hsp.query_offset, 10);
        assert_eq!(hsp.subject_offset, 2);
        let script = hsp.edit_script.as_ref().expect("edit script");
        assert_eq!(script.ops_vec()[0].0, GapAlignOpType::Del);
        assert_eq!(script.ops_vec()[1].0, GapAlignOpType::Ins);
    }

    #[test]
    fn rps_compute_traceback_rejects_missing_required_inputs() {
        let stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[30]);
        let rps_info = make_rps_traceback_info();
        let mut results = HspResults::new(1);

        assert_eq!(
            s_rps_compute_traceback(
                crate::program::RPS_BLAST,
                None,
                Some(&query_info),
                Some(&rps_info),
                &HitSavingOptions::default(),
                Some(&mut results),
                |_hsp_list, _gap_data, _profile_index, _karlin_k| 0,
                None,
                None,
            ),
            -1
        );
        assert_eq!(
            s_rps_compute_traceback(
                crate::program::RPS_BLAST,
                Some(&stream),
                None,
                Some(&rps_info),
                &HitSavingOptions::default(),
                Some(&mut results),
                |_hsp_list, _gap_data, _profile_index, _karlin_k| 0,
                None,
                None,
            ),
            -1
        );
        assert_eq!(
            s_rps_compute_traceback(
                crate::program::RPS_BLAST,
                Some(&stream),
                Some(&query_info),
                None,
                &HitSavingOptions::default(),
                Some(&mut results),
                |_hsp_list, _gap_data, _profile_index, _karlin_k| 0,
                None,
                None,
            ),
            -1
        );
        assert_eq!(
            s_rps_compute_traceback(
                crate::program::RPS_BLAST,
                Some(&stream),
                Some(&query_info),
                Some(&rps_info),
                &HitSavingOptions::default(),
                None,
                |_hsp_list, _gap_data, _profile_index, _karlin_k| 0,
                None,
                None,
            ),
            -1
        );
        assert_eq!(stream.blast_hspstream_write(0, HspList::new(9)), 0);
    }

    #[test]
    fn rps_compute_traceback_accepts_empty_stream_without_callbacks() {
        let stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[30]);
        let rps_info = make_rps_traceback_info();
        let mut results = HspResults::new(1);
        let mut callback_called = false;

        let status = s_rps_compute_traceback(
            crate::program::RPS_BLAST,
            Some(&stream),
            Some(&query_info),
            Some(&rps_info),
            &HitSavingOptions::default(),
            Some(&mut results),
            |_hsp_list, _gap_data, _profile_index, _karlin_k| {
                callback_called = true;
                0
            },
            None,
            None,
        );

        assert_eq!(status, 0);
        assert!(!callback_called);
        assert_eq!(results.hitlists.len(), 1);
        assert!(results.hitlists[0].is_none());
        assert_eq!(stream.blast_hspstream_write(0, HspList::new(9)), -1);
    }

    #[test]
    fn rps_compute_traceback_rejects_malformed_profile_before_stream_drain() {
        let stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[30]);
        let mut bad_rps_info = make_rps_traceback_info();
        bad_rps_info.profile_header.start_offsets = vec![0, 5, 2];
        let mut list = HspList::new(1);
        let mut hsp = make_hsp(20, 1.0e-3);
        hsp.context = 0;
        list.add_hsp(hsp);
        assert_eq!(stream.blast_hspstream_write(0, list), 0);
        let mut results = HspResults::new(1);
        let mut callback_called = false;

        assert_eq!(
            s_rps_compute_traceback(
                crate::program::RPS_BLAST,
                Some(&stream),
                Some(&query_info),
                Some(&bad_rps_info),
                &HitSavingOptions::default(),
                Some(&mut results),
                |_hsp_list, _gap_data, _profile_index, _karlin_k| {
                    callback_called = true;
                    0
                },
                None,
                None,
            ),
            -1
        );

        assert!(!callback_called);
        assert!(results.hitlists[0].is_none());
        assert_eq!(stream.blast_hspstream_write(0, HspList::new(9)), 0);
    }

    #[test]
    fn rps_compute_traceback_sorts_inserted_results_by_evalue() {
        let stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[30]);
        let rps_info = make_rps_traceback_info();
        for &(oid, evalue) in &[(0, 1.0e-2), (1, 1.0e-20)] {
            let mut list = HspList::new(oid);
            let mut hsp = make_hsp(20, evalue);
            hsp.context = 0;
            list.add_hsp(hsp);
            assert_eq!(stream.blast_hspstream_write(0, list), 0);
        }
        let mut results = HspResults::new(1);
        let hit_options = HitSavingOptions {
            hitlist_size: 10,
            ..HitSavingOptions::default()
        };

        let status = s_rps_compute_traceback(
            crate::program::RPS_BLAST,
            Some(&stream),
            Some(&query_info),
            Some(&rps_info),
            &hit_options,
            Some(&mut results),
            |_hsp_list, _gap_data, _profile_index, _karlin_k| 0,
            None,
            None,
        );

        assert_eq!(status, 0);
        let hitlist = results.hitlists[0].as_ref().expect("query hitlist");
        let summary: Vec<(i32, f64)> = hitlist
            .hsp_lists
            .iter()
            .map(|list| (list.oid, list.best_evalue))
            .collect();
        assert_eq!(summary, vec![(1, 1.0e-20), (0, 1.0e-2)]);
    }

    #[test]
    fn rps_compute_traceback_closes_stream_on_traceback_error() {
        let stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[30]);
        let rps_info = make_rps_traceback_info();
        let mut list = HspList::new(1);
        let mut hsp = make_hsp(20, 1.0e-3);
        hsp.context = 0;
        list.add_hsp(hsp);
        assert_eq!(stream.blast_hspstream_write(0, list), 0);
        let mut results = HspResults::new(1);

        let status = s_rps_compute_traceback(
            crate::program::RPS_BLAST,
            Some(&stream),
            Some(&query_info),
            Some(&rps_info),
            &HitSavingOptions::default(),
            Some(&mut results),
            |_hsp_list, _gap_data, _profile_index, _karlin_k| 7,
            None,
            None,
        );

        assert_eq!(status, 7);
        assert!(results.hitlists[0].is_none());
        assert_eq!(stream.blast_hspstream_write(0, HspList::new(9)), -1);
    }

    #[test]
    fn rps_compute_traceback_skips_empty_post_traceback_lists() {
        let stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[30]);
        let rps_info = make_rps_traceback_info();
        let mut list = HspList::new(1);
        let mut hsp = make_hsp(20, 1.0e-3);
        hsp.context = 0;
        list.add_hsp(hsp);
        assert_eq!(stream.blast_hspstream_write(0, list), 0);
        let mut results = HspResults::new(1);
        let mut callback_called = false;

        let status = s_rps_compute_traceback(
            crate::program::RPS_BLAST,
            Some(&stream),
            Some(&query_info),
            Some(&rps_info),
            &HitSavingOptions::default(),
            Some(&mut results),
            |hsp_list, _gap_data, _profile_index, _karlin_k| {
                callback_called = true;
                hsp_list.hsps.clear();
                0
            },
            None,
            None,
        );

        assert_eq!(status, 0);
        assert!(callback_called);
        assert!(results.hitlists[0].is_none());
        assert_eq!(stream.blast_hspstream_write(0, HspList::new(9)), -1);
    }

    #[test]
    fn rps_compute_traceback_closes_stream_on_result_insert_error() {
        let stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);
        let query_info = crate::queryinfo::QueryInfo {
            num_queries: 1,
            contexts: vec![crate::queryinfo::ContextInfo {
                query_offset: 0,
                query_length: 30,
                eff_searchsp: 0,
                length_adjustment: 0,
                query_index: 3,
                frame: 0,
                is_valid: true,
                segment_flags: crate::queryinfo::E_NO_SEGMENTS,
            }],
            max_length: 30,
            min_length: 0,
        };
        let rps_info = make_rps_traceback_info();
        let mut list = HspList::new(1);
        let mut hsp = make_hsp(20, 1.0e-3);
        hsp.context = 0;
        list.add_hsp(hsp);
        assert_eq!(stream.blast_hspstream_write(0, list), 0);
        let mut results = HspResults::new(1);
        let mut callback_called = false;

        let status = s_rps_compute_traceback(
            crate::program::RPS_BLAST,
            Some(&stream),
            Some(&query_info),
            Some(&rps_info),
            &HitSavingOptions::default(),
            Some(&mut results),
            |_hsp_list, _gap_data, _profile_index, _karlin_k| {
                callback_called = true;
                0
            },
            None,
            None,
        );

        assert_eq!(status, 1);
        assert!(callback_called);
        assert!(results.hitlists[0].is_none());
        assert_eq!(stream.blast_hspstream_write(0, HspList::new(9)), -1);
    }

    #[test]
    fn rps_compute_traceback_closes_stream_on_callback_updated_context_insert_error() {
        let stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);
        let query_info = crate::queryinfo::QueryInfo {
            num_queries: 1,
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
                    query_offset: 30,
                    query_length: 30,
                    eff_searchsp: 0,
                    length_adjustment: 0,
                    query_index: 3,
                    frame: 0,
                    is_valid: true,
                    segment_flags: crate::queryinfo::E_NO_SEGMENTS,
                },
            ],
            max_length: 30,
            min_length: 0,
        };
        let rps_info = make_rps_traceback_info();
        let mut list = HspList::new(1);
        let mut hsp = make_hsp(20, 1.0e-3);
        hsp.context = 0;
        list.add_hsp(hsp);
        assert_eq!(stream.blast_hspstream_write(0, list), 0);
        let mut results = HspResults::new(1);

        let status = s_rps_compute_traceback(
            crate::program::RPS_BLAST,
            Some(&stream),
            Some(&query_info),
            Some(&rps_info),
            &HitSavingOptions::default(),
            Some(&mut results),
            |hsp_list, _gap_data, profile_index, _karlin_k| {
                assert_eq!(profile_index, 1);
                hsp_list.hsps[0].context = 1;
                0
            },
            None,
            None,
        );

        assert_eq!(status, 1);
        assert!(results.hitlists[0].is_none());
        assert_eq!(stream.blast_hspstream_write(0, HspList::new(9)), -1);
    }

    #[test]
    fn rps_compute_traceback_closes_stream_on_interrupt_before_traceback() {
        let stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[30]);
        let rps_info = make_rps_traceback_info();
        let mut list = HspList::new(1);
        let mut hsp = make_hsp(20, 1.0e-3);
        hsp.context = 0;
        list.add_hsp(hsp);
        assert_eq!(stream.blast_hspstream_write(0, list), 0);
        let mut results = HspResults::new(1);
        let mut progress = crate::util::s_blast_progress_new(None);
        let mut interrupted = false;
        let mut interrupt = |_progress: &crate::util::SBlastProgress| {
            interrupted = true;
            true
        };

        let status = s_rps_compute_traceback(
            crate::program::RPS_BLAST,
            Some(&stream),
            Some(&query_info),
            Some(&rps_info),
            &HitSavingOptions::default(),
            Some(&mut results),
            |_hsp_list, _gap_data, _profile_index, _karlin_k| {
                panic!("traceback callback must not run after interrupt")
            },
            Some(&mut interrupt),
            Some(&mut progress),
        );

        assert_eq!(status, crate::diagnostics::BLASTERR_INTERRUPTED);
        assert!(interrupted);
        assert_eq!(progress.stage, crate::util::EBlastStage::TracebackSearch);
        assert!(results.hitlists[0].is_none());
        assert_eq!(stream.blast_hspstream_write(0, HspList::new(9)), -1);
    }

    #[test]
    fn rps_compute_traceback_sorts_partial_results_on_interrupt_like_c() {
        let stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[30]);
        let rps_info = make_rps_traceback_info();
        for &(oid, evalue) in &[(0, 1.0e-2), (1, 1.0e-20), (99, 1.0e-5)] {
            let mut list = HspList::new(oid);
            let mut hsp = make_hsp(20, evalue);
            hsp.context = 0;
            list.add_hsp(hsp);
            assert_eq!(stream.blast_hspstream_write(0, list), 0);
        }
        let mut results = HspResults::new(1);
        let mut interrupt_calls = 0usize;
        let mut interrupt = |_progress: &crate::util::SBlastProgress| {
            interrupt_calls += 1;
            interrupt_calls == 3
        };
        let mut progress = crate::util::s_blast_progress_new(None);
        let mut callbacks = 0usize;

        let status = s_rps_compute_traceback(
            crate::program::RPS_BLAST,
            Some(&stream),
            Some(&query_info),
            Some(&rps_info),
            &HitSavingOptions::default(),
            Some(&mut results),
            |_hsp_list, _gap_data, _profile_index, _karlin_k| {
                callbacks += 1;
                0
            },
            Some(&mut interrupt),
            Some(&mut progress),
        );

        assert_eq!(status, crate::diagnostics::BLASTERR_INTERRUPTED);
        assert_eq!(interrupt_calls, 3);
        assert_eq!(callbacks, 2);
        let hitlist = results.hitlists[0].as_ref().expect("partial results");
        let summary: Vec<(i32, f64)> = hitlist
            .hsp_lists
            .iter()
            .map(|list| (list.oid, list.best_evalue))
            .collect();
        assert_eq!(summary, vec![(1, 1.0e-20), (0, 1.0e-2)]);
        assert_eq!(stream.blast_hspstream_write(0, HspList::new(9)), -1);
    }

    #[test]
    fn rps_compute_traceback_skips_hsplists_without_profile_context() {
        let stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[30]);
        let rps_info = make_rps_traceback_info();
        let mut list = HspList::new(99);
        let mut hsp = make_hsp(20, 1.0e-3);
        hsp.context = 0;
        list.add_hsp(hsp);
        assert_eq!(stream.blast_hspstream_write(0, list), 0);
        let mut results = HspResults::new(1);
        let mut callback_called = false;

        let status = s_rps_compute_traceback(
            crate::program::RPS_BLAST,
            Some(&stream),
            Some(&query_info),
            Some(&rps_info),
            &HitSavingOptions::default(),
            Some(&mut results),
            |_hsp_list, _gap_data, _profile_index, _karlin_k| {
                callback_called = true;
                0
            },
            None,
            None,
        );

        assert_eq!(status, 0);
        assert!(!callback_called);
        assert!(results.hitlists[0].is_none());
        assert_eq!(stream.blast_hspstream_write(0, HspList::new(9)), -1);
    }

    #[test]
    fn rps_compute_traceback_skips_negative_oid_without_profile_context() {
        let stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[30]);
        let rps_info = make_rps_traceback_info();
        let mut list = HspList::new(-1);
        let mut hsp = make_hsp(20, 1.0e-3);
        hsp.context = 0;
        list.add_hsp(hsp);
        assert_eq!(stream.blast_hspstream_write(0, list), 0);
        let mut results = HspResults::new(1);
        let mut callback_called = false;

        let status = s_rps_compute_traceback(
            crate::program::RPS_BLAST,
            Some(&stream),
            Some(&query_info),
            Some(&rps_info),
            &HitSavingOptions::default(),
            Some(&mut results),
            |_hsp_list, _gap_data, _profile_index, _karlin_k| {
                callback_called = true;
                0
            },
            None,
            None,
        );

        assert_eq!(status, 0);
        assert!(!callback_called);
        assert!(results.hitlists[0].is_none());
        assert_eq!(stream.blast_hspstream_write(0, HspList::new(9)), -1);
    }

    #[test]
    fn rps_compute_traceback_skips_invalid_zero_length_profile_context() {
        let stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[30]);
        let alphabet_size = 26;
        let num_rows = 2 + 1;
        let rps_info = RpsTracebackInfo {
            profile_header: RpsProfileHeader {
                magic_number: RPS_MAGIC_NUM,
                num_profiles: 2,
                start_offsets: vec![0, 0, 2],
                pssm_values: (0..(num_rows * alphabet_size)).map(|v| v as i32).collect(),
            },
            freq_ratios_header: None,
            karlin_k: vec![0.1, 0.2],
        };
        let mut list = HspList::new(0);
        let mut hsp = make_hsp(20, 1.0e-3);
        hsp.context = 0;
        list.add_hsp(hsp);
        assert_eq!(stream.blast_hspstream_write(0, list), 0);
        let mut results = HspResults::new(1);
        let mut callback_called = false;

        let status = s_rps_compute_traceback(
            crate::program::RPS_BLAST,
            Some(&stream),
            Some(&query_info),
            Some(&rps_info),
            &HitSavingOptions::default(),
            Some(&mut results),
            |_hsp_list, _gap_data, _profile_index, _karlin_k| {
                callback_called = true;
                0
            },
            None,
            None,
        );

        assert_eq!(status, 0);
        assert!(!callback_called);
        assert!(results.hitlists[0].is_none());
        assert_eq!(stream.blast_hspstream_write(0, HspList::new(9)), -1);
    }

    #[test]
    fn rps_compute_traceback_continues_after_invalid_zero_length_profile_context() {
        let stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[30]);
        let alphabet_size = 26;
        let num_rows = 2 + 1;
        let rps_info = RpsTracebackInfo {
            profile_header: RpsProfileHeader {
                magic_number: RPS_MAGIC_NUM,
                num_profiles: 2,
                start_offsets: vec![0, 0, 2],
                pssm_values: (0..(num_rows * alphabet_size)).map(|v| v as i32).collect(),
            },
            freq_ratios_header: None,
            karlin_k: vec![0.1, 0.2],
        };
        let mut invalid_list = HspList::new(0);
        let mut invalid_hsp = make_hsp(20, 1.0e-3);
        invalid_hsp.context = 0;
        invalid_list.add_hsp(invalid_hsp);
        assert_eq!(stream.blast_hspstream_write(0, invalid_list), 0);

        let mut valid_list = HspList::new(1);
        let mut valid_hsp = make_hsp(30, 1.0e-12);
        valid_hsp.context = 0;
        valid_list.add_hsp(valid_hsp);
        assert_eq!(stream.blast_hspstream_write(0, valid_list), 0);
        let mut results = HspResults::new(1);
        let mut callback_profiles = Vec::new();

        let status = s_rps_compute_traceback(
            crate::program::RPS_BLAST,
            Some(&stream),
            Some(&query_info),
            Some(&rps_info),
            &HitSavingOptions::default(),
            Some(&mut results),
            |_hsp_list, gap_data, profile_index, karlin_k| {
                callback_profiles.push(profile_index);
                assert_eq!(profile_index, 1);
                assert_eq!(gap_data.concat_db_info.contexts[1].query_offset, 0);
                assert_eq!(gap_data.concat_db_info.contexts[1].query_length, 2);
                assert!((karlin_k.expect("RPS K") - 0.24).abs() < 1e-12);
                0
            },
            None,
            None,
        );

        assert_eq!(status, 0);
        assert_eq!(callback_profiles, vec![1]);
        let hitlist = results.hitlists[0].as_ref().expect("valid profile result");
        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(hitlist.hsp_lists[0].oid, 1);
        assert_eq!(stream.blast_hspstream_write(0, HspList::new(9)), -1);
    }

    #[test]
    fn rps_compute_traceback_rejects_malformed_frequency_payload_before_stream_drain() {
        let stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[30]);
        let mut rps_info = make_rps_traceback_info();
        let frequency_rows = rps_info
            .freq_ratios_header
            .as_mut()
            .expect("frequency ratios");
        frequency_rows
            .freq_values
            .truncate(2 * crate::encoding::BLASTAA_SIZE);
        let mut list = HspList::new(1);
        let mut hsp = make_hsp(20, 1.0e-3);
        hsp.context = 0;
        list.add_hsp(hsp);
        assert_eq!(stream.blast_hspstream_write(0, list), 0);
        let mut results = HspResults::new(1);
        let mut callback_called = false;

        let status = s_rps_compute_traceback(
            crate::program::RPS_BLAST,
            Some(&stream),
            Some(&query_info),
            Some(&rps_info),
            &HitSavingOptions::default(),
            Some(&mut results),
            |_hsp_list, _gap_data, _profile_index, _karlin_k| {
                callback_called = true;
                0
            },
            None,
            None,
        );

        assert_eq!(status, -1);
        assert!(!callback_called);
        assert!(results.hitlists[0].is_none());
        assert_eq!(stream.blast_hspstream_write(0, HspList::new(9)), 0);
    }

    #[test]
    fn rps_compute_traceback_allows_missing_karlin_k_metadata() {
        let stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[30]);
        let mut rps_info = make_rps_traceback_info();
        rps_info.karlin_k.clear();
        let mut list = HspList::new(1);
        let mut hsp = make_hsp(20, 1.0e-3);
        hsp.context = 0;
        list.add_hsp(hsp);
        assert_eq!(stream.blast_hspstream_write(0, list), 0);
        let mut results = HspResults::new(1);
        let mut callback_called = false;

        let status = s_rps_compute_traceback(
            crate::program::RPS_BLAST,
            Some(&stream),
            Some(&query_info),
            Some(&rps_info),
            &HitSavingOptions::default(),
            Some(&mut results),
            |_hsp_list, _gap_data, profile_index, karlin_k| {
                callback_called = true;
                assert_eq!(profile_index, 1);
                assert!(karlin_k.is_none());
                0
            },
            None,
            None,
        );

        assert_eq!(status, 0);
        assert!(callback_called);
        assert!(results.hitlists[0].is_some());
        assert_eq!(stream.blast_hspstream_write(0, HspList::new(9)), -1);
    }

    #[test]
    fn rps_compute_traceback_gates_rps_karlin_k_like_c() {
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[30]);
        let rps_info = make_rps_traceback_info();
        let hit_options = HitSavingOptions::default();

        let non_cbs_stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);
        let mut non_cbs_list = HspList::new(1);
        let mut non_cbs_hsp = make_hsp(20, 1.0e-3);
        non_cbs_hsp.query_frame = 1;
        non_cbs_list.add_hsp(non_cbs_hsp);
        assert_eq!(non_cbs_stream.blast_hspstream_write(0, non_cbs_list), 0);
        let mut non_cbs_results = HspResults::new(1);
        let mut non_cbs_called = false;

        let status = s_rps_compute_traceback(
            crate::program::RPS_TBLASTN,
            Some(&non_cbs_stream),
            Some(&query_info),
            Some(&rps_info),
            &hit_options,
            Some(&mut non_cbs_results),
            |_hsp_list, _gap_data, profile_index, karlin_k| {
                non_cbs_called = true;
                assert_eq!(profile_index, 1);
                assert!(karlin_k.is_none());
                0
            },
            None,
            None,
        );

        assert_eq!(status, 0);
        assert!(non_cbs_called);

        let cbs_stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);
        let mut cbs_list = HspList::new(1);
        let mut cbs_hsp = make_hsp(20, 1.0e-3);
        cbs_hsp.query_frame = 1;
        cbs_list.add_hsp(cbs_hsp);
        assert_eq!(cbs_stream.blast_hspstream_write(0, cbs_list), 0);
        let mut cbs_results = HspResults::new(1);
        let mut cbs_called = false;

        let status = s_rps_compute_traceback_with_composition(
            crate::program::RPS_TBLASTN,
            Some(&cbs_stream),
            Some(&query_info),
            Some(&rps_info),
            true,
            &hit_options,
            Some(&mut cbs_results),
            |_hsp_list, _gap_data, profile_index, karlin_k| {
                cbs_called = true;
                assert_eq!(profile_index, 1);
                assert!((karlin_k.expect("CBS K") - 0.24).abs() < 1e-12);
                0
            },
            None,
            None,
        );

        assert_eq!(status, 0);
        assert!(cbs_called);

        let rps_cbs_stream = blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);
        let mut rps_cbs_list = HspList::new(1);
        let rps_cbs_hsp = make_hsp(20, 1.0e-3);
        rps_cbs_list.add_hsp(rps_cbs_hsp);
        assert_eq!(rps_cbs_stream.blast_hspstream_write(0, rps_cbs_list), 0);
        let mut rps_cbs_results = HspResults::new(1);
        let mut rps_cbs_called = false;

        let status = s_rps_compute_traceback_with_composition(
            crate::program::RPS_BLAST,
            Some(&rps_cbs_stream),
            Some(&query_info),
            Some(&rps_info),
            true,
            &hit_options,
            Some(&mut rps_cbs_results),
            |_hsp_list, _gap_data, profile_index, karlin_k| {
                rps_cbs_called = true;
                assert_eq!(profile_index, 1);
                assert!((karlin_k.expect("RPS CBS K") - 0.24).abs() < 1e-12);
                0
            },
            None,
            None,
        );

        assert_eq!(status, 0);
        assert!(rps_cbs_called);
    }
}
