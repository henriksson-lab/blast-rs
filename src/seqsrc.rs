//! Rust equivalent of blast_seqsrc.c — sequence source trait.
//! Replaces the C vtable-based BlastSeqSrc with a Rust trait.

use std::any::Any;
use std::sync::Arc;

use crate::util::{BlastSequenceBlk, SSeqRange};

/// Direct Rust mirror of NCBI `EBlastEncoding` (`blast_encoding.h`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum EBlastEncoding {
    Protein = 0,
    Nucleotide = 1,
    Ncbi4na = 2,
    Ncbi2na = 3,
    Error = 255,
}

impl Default for EBlastEncoding {
    fn default() -> Self {
        Self::Protein
    }
}

/// Encoding for sequence retrieval.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SeqEncoding {
    /// NCBIstdaa for protein, NCBI2na packed for nucleotide (no sentinels)
    Protein,
    /// BLASTNA decoded with sentinel bytes
    Nucleotide,
    /// NCBI4na encoding
    Ncbi4na,
    /// NCBI2na encoding
    Ncbi2na,
    /// Invalid encoding sentinel.
    Error,
}

impl From<EBlastEncoding> for SeqEncoding {
    fn from(value: EBlastEncoding) -> Self {
        match value {
            EBlastEncoding::Protein => SeqEncoding::Protein,
            EBlastEncoding::Nucleotide => SeqEncoding::Nucleotide,
            EBlastEncoding::Ncbi4na => SeqEncoding::Ncbi4na,
            EBlastEncoding::Ncbi2na => SeqEncoding::Ncbi2na,
            EBlastEncoding::Error => SeqEncoding::Error,
        }
    }
}

impl From<SeqEncoding> for EBlastEncoding {
    fn from(value: SeqEncoding) -> Self {
        match value {
            SeqEncoding::Protein => EBlastEncoding::Protein,
            SeqEncoding::Nucleotide => EBlastEncoding::Nucleotide,
            SeqEncoding::Ncbi4na => EBlastEncoding::Ncbi4na,
            SeqEncoding::Ncbi2na => EBlastEncoding::Ncbi2na,
            SeqEncoding::Error => EBlastEncoding::Error,
        }
    }
}

/// Arguments for fetching a sequence.
#[derive(Debug, Clone, PartialEq)]
pub struct GetSeqArg {
    pub oid: i32,
    pub encoding: SeqEncoding,
    pub reset_ranges: bool,
    pub check_oid_exclusion: bool,
    pub ranges: Option<BlastSeqSrcSetRangesArg>,
}

impl Default for GetSeqArg {
    fn default() -> Self {
        Self {
            oid: 0,
            encoding: SeqEncoding::Protein,
            reset_ranges: false,
            check_oid_exclusion: false,
            ranges: None,
        }
    }
}

/// Direct Rust mirror of NCBI `BlastSeqSrcGetSeqArg` (`blast_seqsrc.h`).
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct BlastSeqSrcGetSeqArg {
    pub oid: i32,
    pub encoding: EBlastEncoding,
    pub reset_ranges: bool,
    pub check_oid_exclusion: bool,
    pub seq: Option<BlastSequenceBlk>,
    pub ranges: Option<BlastSeqSrcSetRangesArg>,
}

impl BlastSeqSrcGetSeqArg {
    pub fn as_get_seq_arg(&self) -> GetSeqArg {
        GetSeqArg {
            oid: self.oid,
            encoding: self.encoding.into(),
            reset_ranges: self.reset_ranges,
            check_oid_exclusion: self.check_oid_exclusion,
            ranges: self.ranges.clone(),
        }
    }
}

/// Fetched sequence data.
pub struct SeqData {
    pub sequence: Vec<u8>,
    pub length: i32,
}

pub const BLAST_SEQSRC_MINGAP: i32 = 1024;
pub const BLAST_SEQSRC_OVERHANG: i32 = 1024;
pub const BLAST_SEQSRC_MINLENGTH: i32 = 10;

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct BlastSeqSrcSetRangesArg {
    pub oid: i32,
    pub capacity: i32,
    pub num_ranges: i32,
    pub ranges: Vec<i32>,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct BlastGiList {
    pub data: Vec<i32>,
    pub num_allocated: usize,
    pub num_used: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BlastSeqSrcIterator {
    pub oid_list: Vec<i32>,
    pub chunk_sz: u32,
    pub current_pos: u32,
}

/// The sequence source trait — Rust replacement for BlastSeqSrc vtable.
pub trait BlastSeqSource: Send + Sync {
    /// Number of sequences in the database.
    fn num_seqs(&self) -> i32;

    /// Number of sequences from database statistics.
    fn num_seqs_stats(&self) -> i32 {
        self.num_seqs()
    }

    /// Total length of all sequences in bases/residues.
    fn total_length(&self) -> i64;

    /// Total length from database statistics.
    fn total_length_stats(&self) -> i64 {
        self.total_length()
    }

    /// Maximum sequence length.
    fn max_seq_len(&self) -> i32;

    /// Average sequence length.
    fn avg_seq_len(&self) -> i32;

    /// Minimum sequence length.
    fn min_seq_len(&self) -> i32 {
        BLAST_SEQSRC_MINLENGTH
    }

    /// Database name.
    fn name(&self) -> &str;

    /// Whether the database contains protein sequences.
    fn is_protein(&self) -> bool;

    /// Get the length of a specific sequence.
    fn seq_len(&self, oid: i32) -> i32;

    /// Fetch sequence data for the given OID with the specified encoding.
    fn get_sequence(&self, arg: &GetSeqArg) -> Option<SeqData>;

    /// Fetch sequence data through NCBI's `BlastSeqSrcGetSeqArg` object.
    fn get_sequence_blk(&self, arg: &mut BlastSeqSrcGetSeqArg) -> Option<SeqData> {
        if arg.reset_ranges && !self.is_protein() {
            arg.ranges = None;
        }
        let data = self.get_sequence(&arg.as_get_seq_arg())?;
        let copied = matches!(
            arg.encoding,
            EBlastEncoding::Nucleotide | EBlastEncoding::Ncbi4na
        );
        let has_sentinel_byte = arg.encoding == EBlastEncoding::Nucleotide;

        let mut seq = arg.seq.take().unwrap_or_default();
        seq.sequence_start = None;
        seq.sequence = None;
        seq.sequence_start_allocated = false;
        seq.sequence_allocated = false;

        if copied {
            seq.sequence_start_allocated = true;
            seq.sequence_start = Some(data.sequence.clone());
            seq.sequence = if has_sentinel_byte {
                Some(data.sequence.get(1..).unwrap_or(&[]).to_vec())
            } else {
                seq.sequence_start.clone()
            };
        } else {
            seq.sequence_allocated = true;
            seq.sequence = Some(data.sequence.clone());
        }

        seq.sequence_start_nomask = seq.sequence_start.clone();
        seq.sequence_nomask = seq.sequence.clone();
        seq.nomask_allocated = false;
        seq.length = data.length;
        seq.oid = arg.oid;
        seq.bases_offset = 0;
        arg.seq = Some(seq);
        Some(data)
    }

    /// Release a sequence obtained through `get_sequence_blk`.
    fn release_sequence(&self, arg: &mut BlastSeqSrcGetSeqArg) {
        if let Some(seq) = arg.seq.as_mut() {
            if seq.sequence_start_allocated {
                seq.sequence_start = None;
                seq.sequence_start_allocated = false;
                seq.sequence = None;
            }
            if seq.sequence_allocated {
                seq.sequence = None;
                seq.sequence_allocated = false;
            }
        }
        arg.ranges = None;
    }

    /// Iterator over OIDs.
    fn iter_oids(&self) -> Box<dyn Iterator<Item = i32> + '_>;

    /// Optional hook matching BlastSeqSrcSetNumberOfThreads.
    fn set_number_of_threads(&mut self, _n_threads: i32) {}

    /// Optional hook matching BlastSeqSrcGetSupportsPartialFetching.
    fn supports_partial_fetching(&self) -> bool {
        false
    }

    /// Optional hook matching BlastSeqSrcSetSeqRanges.
    fn set_seq_ranges(&mut self, _ranges: Option<&[SSeqRange]>) {}

    /// Optional hook matching BlastSeqSrcGetGis.
    fn get_gis(&self, _oid: i32) -> Option<BlastGiList> {
        None
    }

    /// Optional hook matching BlastSeqSrcIteratorNext.
    fn iter_next(&self, itr: &mut BlastSeqSrcIterator) -> i32 {
        if itr.current_pos == u32::MAX {
            itr.oid_list = self.iter_oids().collect();
            itr.current_pos = 0;
        }
        let Some(&oid) = itr.oid_list.get(itr.current_pos as usize) else {
            return SEQSRC_EOF;
        };
        itr.current_pos += 1;
        oid
    }

    /// Optional hook matching BlastSeqSrcResetChunkIterator.
    fn reset_chunk_iterator(&mut self) {}
}

/// EOF sentinel for sequence source iteration.
pub const SEQSRC_EOF: i32 = -1;

const INITIAL_GI_LIST_SIZE: usize = 10;
pub const K_BLAST_SEQ_SRC_DEFAULT_CHUNK_SIZE: u32 = 1024;

pub type BlastSeqSrcConstructor =
    fn(Option<Box<dyn Any + Send>>) -> Option<Arc<dyn BlastSeqSource>>;
pub type BlastSeqSrcDestructor = fn(&mut BlastSeqSrc);
pub type BlastSeqSrcCopier = fn(&BlastSeqSrc) -> Option<BlastSeqSrc>;

#[derive(Clone)]
pub struct BlastSeqSrc {
    pub source: Option<Arc<dyn BlastSeqSource>>,
    pub init_error: Option<String>,
    pub delete_fn: Option<BlastSeqSrcDestructor>,
    pub copy_fn: Option<BlastSeqSrcCopier>,
}

pub struct BlastSeqSrcNewInfo {
    pub constructor: Option<BlastSeqSrcConstructor>,
    pub ctor_argument: Option<Box<dyn Any + Send>>,
    pub delete_fn: Option<BlastSeqSrcDestructor>,
    pub copy_fn: Option<BlastSeqSrcCopier>,
}

/// Port of NCBI `BlastSeqSrcNew` (`blast_seqsrc.c:90`).
pub fn blast_seq_src_new(bssn_info: Option<BlastSeqSrcNewInfo>) -> Option<BlastSeqSrc> {
    let mut bssn_info = bssn_info?;
    let constructor = bssn_info.constructor?;
    let source = constructor(bssn_info.ctor_argument.take())?;
    Some(BlastSeqSrc {
        source: Some(source),
        init_error: None,
        delete_fn: bssn_info.delete_fn,
        copy_fn: bssn_info.copy_fn,
    })
}

/// Rust ownership equivalent of NCBI `BlastSeqSrcFree`.
pub fn blast_seq_src_free(mut seq_src: Option<BlastSeqSrc>) -> Option<BlastSeqSrc> {
    if let Some(seq_src) = seq_src.as_mut() {
        seq_src.init_error.take();
        if let Some(delete_fn) = seq_src.delete_fn {
            delete_fn(seq_src);
        }
        seq_src.source.take();
    }
    None
}

/// Port of NCBI `BlastSeqSrcCopy` (`blast_seqsrc.c:138`).
pub fn blast_seq_src_copy(seq_src: Option<&BlastSeqSrc>) -> Option<BlastSeqSrc> {
    let seq_src = seq_src?;
    if let Some(copy_fn) = seq_src.copy_fn {
        return copy_fn(seq_src);
    }
    Some(seq_src.clone())
}

/// Port of NCBI `BlastSeqSrcSetRangesArgNew` (`blast_seqsrc.c:476`).
pub fn blast_seq_src_set_ranges_arg_new(num_ranges: i32) -> BlastSeqSrcSetRangesArg {
    BlastSeqSrcSetRangesArg {
        oid: 0,
        capacity: num_ranges,
        num_ranges: 0,
        ranges: Vec::with_capacity((2 * num_ranges.max(0)) as usize),
    }
}

/// Port of NCBI `BlastSeqSrcSetRangesArgFree` (`blast_seqsrc.c:487`).
pub fn blast_seq_src_set_ranges_arg_free(
    arg: &mut Option<BlastSeqSrcSetRangesArg>,
) -> Option<BlastSeqSrcSetRangesArg> {
    *arg = None;
    None
}

/// Port of NCBI `BlastSeqSrcSetRangesArgAddRange` (`blast_seqsrc.c:499`).
pub fn blast_seq_src_set_ranges_arg_add_range(
    arg: &mut BlastSeqSrcSetRangesArg,
    begin: i32,
    end: i32,
    len: i32,
) -> i16 {
    if arg.num_ranges + 2 > arg.capacity {
        let new_size = arg.capacity.saturating_mul(2).max(arg.num_ranges + 2);
        arg.ranges
            .reserve((new_size.saturating_mul(2).max(0)) as usize);
        arg.capacity = new_size;
    }
    let begin = 0.max(begin - BLAST_SEQSRC_OVERHANG);
    let end = len.min(end + BLAST_SEQSRC_OVERHANG);
    arg.ranges.push(begin);
    arg.ranges.push(end);
    arg.num_ranges += 2;
    0
}

/// Port of NCBI `BlastSeqSrcSetRangesArgBuild` (`blast_seqsrc.c:520`).
pub fn blast_seq_src_set_ranges_arg_build(arg: &mut BlastSeqSrcSetRangesArg) {
    arg.num_ranges /= 2;
    if arg.num_ranges <= 1 {
        return;
    }

    let total = arg.num_ranges as usize;
    arg.ranges.truncate(total * 2);
    let mut pairs: Vec<[i32; 2]> = arg
        .ranges
        .chunks_exact(2)
        .map(|chunk| [chunk[0], chunk[1]])
        .collect();
    pairs.sort_by_key(|pair| pair[0]);

    let mut merged: Vec<[i32; 2]> = Vec::with_capacity(pairs.len());
    for pair in pairs {
        if let Some(last) = merged.last_mut() {
            if pair[0] > last[1] + BLAST_SEQSRC_MINGAP {
                merged.push(pair);
            } else if pair[1] > last[1] {
                last[1] = pair[1];
            }
        } else {
            merged.push(pair);
        }
    }

    arg.num_ranges = merged.len() as i32;
    arg.ranges.clear();
    for pair in merged {
        arg.ranges.extend_from_slice(&pair);
    }
}

/// Port of NCBI `BLAST_SetupPartialFetching` (`blast_traceback.c:1336`).
///
/// Builds the subject ranges needed to fetch only traceback neighborhoods for
/// all HSPs in a same-oid batch. Translated-subject programs expand protein
/// offsets back to nucleotide coordinates and reverse negative-frame ranges,
/// matching the C helper before range overhang/merge handling.
pub fn blast_setup_partial_fetching(
    program_number: crate::program::ProgramType,
    seq_src: &dyn BlastSeqSource,
    hsplist_array: &[&crate::hspstream::HspList],
) -> Option<BlastSeqSrcSetRangesArg> {
    let first = hsplist_array.first()?;
    let oid = first.oid;
    let len = seq_src.seq_len(oid);
    let num_hsps: i32 = hsplist_array
        .iter()
        .map(|hsp_list| hsp_list.hsps.len() as i32)
        .sum();
    let mut arg = blast_seq_src_set_ranges_arg_new(num_hsps);
    arg.oid = oid;

    for hsp_list in hsplist_array {
        for hsp in &hsp_list.hsps {
            let mut begin = hsp.subject_offset;
            let mut end = hsp.subject_end;
            if crate::program::blast_subject_is_translated(program_number) {
                begin = (begin - 2) * crate::util::CODON_LENGTH as i32;
                end = (end + 2) * crate::util::CODON_LENGTH as i32;
                if hsp.subject_frame < 0 {
                    let begin_new = len - end;
                    end = len - begin;
                    begin = begin_new;
                }
            }
            if blast_seq_src_set_ranges_arg_add_range(&mut arg, begin, end, len) != 0 {
                return None;
            }
        }
    }

    blast_seq_src_set_ranges_arg_build(&mut arg);
    Some(arg)
}

/// Port of NCBI `Blast_GiListNew` (`blast_seqsrc.c:335`).
pub fn blast_gi_list_new() -> BlastGiList {
    blast_gi_list_new_ex(INITIAL_GI_LIST_SIZE)
}

/// Port of NCBI `Blast_GiListNewEx` (`blast_seqsrc.c:341`).
pub fn blast_gi_list_new_ex(list_size: usize) -> BlastGiList {
    BlastGiList {
        data: Vec::with_capacity(list_size),
        num_allocated: list_size,
        num_used: 0,
    }
}

/// Port of NCBI `Blast_GiListFree` (`blast_seqsrc.c:356`).
pub fn blast_gi_list_free(gilist: &mut Option<BlastGiList>) -> Option<BlastGiList> {
    *gilist = None;
    None
}

/// Port of NCBI `Blast_GiList_ReallocIfNecessary` (`blast_seqsrc.c:370`).
pub fn blast_gi_list_realloc_if_necessary(gilist: &mut BlastGiList) -> i16 {
    if gilist.num_used + 1 > gilist.num_allocated {
        gilist.num_allocated = gilist
            .num_allocated
            .saturating_mul(2)
            .max(gilist.num_used + 1);
        gilist
            .data
            .reserve(gilist.num_allocated.saturating_sub(gilist.data.capacity()));
    }
    0
}

/// Port of NCBI `Blast_GiList_Append` (`blast_seqsrc.c:386`).
pub fn blast_gi_list_append(gilist: &mut BlastGiList, gi: i32) -> i16 {
    let retval = blast_gi_list_realloc_if_necessary(gilist);
    if retval != 0 {
        return retval;
    }
    gilist.data.push(gi);
    gilist.num_used += 1;
    retval
}

/// Port of NCBI `BlastSeqSrcIteratorNew` (`blast_seqsrc.c:414`).
pub fn blast_seq_src_iterator_new() -> BlastSeqSrcIterator {
    blast_seq_src_iterator_new_ex(0)
}

/// Port of NCBI `BlastSeqSrcIteratorNewEx` (`blast_seqsrc.c:421`).
pub fn blast_seq_src_iterator_new_ex(chunk_sz: u32) -> BlastSeqSrcIterator {
    let chunk_sz = if chunk_sz == 0 {
        K_BLAST_SEQ_SRC_DEFAULT_CHUNK_SIZE
    } else {
        chunk_sz
    };
    BlastSeqSrcIterator {
        oid_list: Vec::with_capacity(chunk_sz as usize),
        chunk_sz,
        current_pos: u32::MAX,
    }
}

/// Port of NCBI `BlastSeqSrcIteratorFree` (`blast_seqsrc.c:445`).
pub fn blast_seq_src_iterator_free(
    itr: &mut Option<BlastSeqSrcIterator>,
) -> Option<BlastSeqSrcIterator> {
    *itr = None;
    None
}

/// Port of NCBI `BlastSeqSrcIteratorNext` (`blast_seqsrc.c:457`).
pub fn blast_seq_src_iterator_next(
    seq_src: &dyn BlastSeqSource,
    itr: &mut BlastSeqSrcIterator,
) -> i32 {
    seq_src.iter_next(itr)
}

/// Port of NCBI `BlastSeqSrcResetChunkIterator` (`blast_seqsrc.c:466`).
pub fn blast_seq_src_reset_chunk_iterator(seq_src: &mut dyn BlastSeqSource) {
    seq_src.reset_chunk_iterator();
}

/// Port of NCBI `BlastSeqSrcSetNumberOfThreads` (`blast_seqsrc.c:200`).
pub fn blast_seq_src_set_number_of_threads(
    seq_src: Option<&mut dyn BlastSeqSource>,
    n_threads: i32,
) {
    if let Some(seq_src) = seq_src {
        seq_src.set_number_of_threads(n_threads);
    }
}

/// Port of NCBI `BlastSeqSrcGetTotLen` (`blast_seqsrc.c:251`).
pub fn blast_seq_src_get_tot_len(seq_src: &dyn BlastSeqSource) -> i64 {
    seq_src.total_length()
}

/// Port of NCBI `BlastSeqSrcGetTotLenStats` (`blast_seqsrc.c:259`).
pub fn blast_seq_src_get_tot_len_stats(seq_src: &dyn BlastSeqSource) -> i64 {
    seq_src.total_length_stats()
}

/// Port of NCBI `BlastSeqSrcGetNumSeqs` (`blast_seqsrc.c:177`).
pub fn blast_seq_src_get_num_seqs(seq_src: &dyn BlastSeqSource) -> i32 {
    seq_src.num_seqs()
}

/// Port of NCBI `BlastSeqSrcGetNumSeqsStats` (`blast_seqsrc.c:185`).
pub fn blast_seq_src_get_num_seqs_stats(seq_src: &dyn BlastSeqSource) -> i32 {
    seq_src.num_seqs_stats()
}

/// Port of NCBI `BlastSeqSrcGetMaxSeqLen` (`blast_seqsrc.c:193`).
pub fn blast_seq_src_get_max_seq_len(seq_src: &dyn BlastSeqSource) -> i32 {
    seq_src.max_seq_len()
}

/// Port of NCBI `BlastSeqSrcGetMinSeqLen` (`blast_seqsrc.c:201`).
pub fn blast_seq_src_get_min_seq_len(seq_src: &dyn BlastSeqSource) -> i32 {
    seq_src.min_seq_len()
}

/// Port of NCBI `BlastSeqSrcGetAvgSeqLen` (`blast_seqsrc.c:211`).
pub fn blast_seq_src_get_avg_seq_len(seq_src: &dyn BlastSeqSource) -> i32 {
    seq_src.avg_seq_len()
}

/// Port of NCBI `BlastSeqSrcGetName` (`blast_seqsrc.c:267`).
pub fn blast_seq_src_get_name(seq_src: &dyn BlastSeqSource) -> &str {
    seq_src.name()
}

/// Port of NCBI `BlastSeqSrcGetIsProt` (`blast_seqsrc.c:275`).
pub fn blast_seq_src_get_is_prot(seq_src: &dyn BlastSeqSource) -> bool {
    seq_src.is_protein()
}

/// Port of NCBI `BlastSeqSrcGetSequence` (`blast_seqsrc.c:271`).
pub fn blast_seq_src_get_sequence(
    seq_src: Option<&dyn BlastSeqSource>,
    getseq_arg: Option<&mut BlastSeqSrcGetSeqArg>,
) -> Option<SeqData> {
    seq_src?.get_sequence_blk(getseq_arg?)
}

/// Port of NCBI `BlastSeqSrcReleaseSequence` (`blast_seqsrc.c:290`).
pub fn blast_seq_src_release_sequence(
    seq_src: Option<&dyn BlastSeqSource>,
    getseq_arg: Option<&mut BlastSeqSrcGetSeqArg>,
) {
    if let (Some(seq_src), Some(getseq_arg)) = (seq_src, getseq_arg) {
        seq_src.release_sequence(getseq_arg);
    }
}

/// Port of NCBI `BlastSeqSrcGetSeqLen` (`blast_seqsrc.c:281`).
pub fn blast_seq_src_get_seq_len(seq_src: &dyn BlastSeqSource, oid: i32) -> i32 {
    seq_src.seq_len(oid)
}

/// Port of NCBI `BlastSeqSrcGetSupportsPartialFetching` (`blast_seqsrc.c:283`).
pub fn blast_seq_src_get_supports_partial_fetching(seq_src: &dyn BlastSeqSource) -> bool {
    seq_src.supports_partial_fetching()
}

/// Port of NCBI `BlastSeqSrcSetSeqRanges` (`blast_seqsrc.c:293`).
pub fn blast_seq_src_set_seq_ranges(
    seq_src: &mut dyn BlastSeqSource,
    ranges: Option<&[SSeqRange]>,
) {
    seq_src.set_seq_ranges(ranges);
}

/// Port of NCBI `BlastSeqSrcGetGis` (`blast_seqsrc.c:400`).
pub fn blast_seq_src_get_gis(seq_src: &dyn BlastSeqSource, oid: i32) -> Option<BlastGiList> {
    seq_src.get_gis(oid)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Mutex;

    struct MockSeqSrc {
        seqs: Vec<Vec<u8>>,
    }

    struct TrackingSeqSrc {
        seqs: Vec<Vec<u8>>,
        is_protein: bool,
        supports_partial: bool,
        n_threads: i32,
        ranges: Option<Vec<SSeqRange>>,
        gis: Vec<i32>,
        reset_count: i32,
    }

    struct RecordingSeqSrc {
        seqs: Vec<Vec<u8>>,
        is_protein: bool,
        last_arg: Mutex<Option<GetSeqArg>>,
    }

    impl BlastSeqSource for MockSeqSrc {
        fn num_seqs(&self) -> i32 {
            self.seqs.len() as i32
        }
        fn total_length(&self) -> i64 {
            self.seqs.iter().map(|s| s.len() as i64).sum()
        }
        fn max_seq_len(&self) -> i32 {
            self.seqs.iter().map(|s| s.len() as i32).max().unwrap_or(0)
        }
        fn avg_seq_len(&self) -> i32 {
            if self.seqs.is_empty() {
                0
            } else {
                (self.total_length() / self.num_seqs() as i64) as i32
            }
        }
        fn name(&self) -> &str {
            "mock"
        }
        fn is_protein(&self) -> bool {
            false
        }
        fn seq_len(&self, oid: i32) -> i32 {
            self.seqs[oid as usize].len() as i32
        }
        fn get_sequence(&self, arg: &GetSeqArg) -> Option<SeqData> {
            let seq = self.seqs.get(arg.oid as usize)?;
            Some(SeqData {
                sequence: seq.clone(),
                length: seq.len() as i32,
            })
        }
        fn iter_oids(&self) -> Box<dyn Iterator<Item = i32> + '_> {
            Box::new(0..self.num_seqs())
        }
    }

    impl BlastSeqSource for TrackingSeqSrc {
        fn num_seqs(&self) -> i32 {
            self.seqs.len() as i32
        }
        fn total_length(&self) -> i64 {
            self.seqs.iter().map(|s| s.len() as i64).sum()
        }
        fn max_seq_len(&self) -> i32 {
            self.seqs.iter().map(|s| s.len() as i32).max().unwrap_or(0)
        }
        fn avg_seq_len(&self) -> i32 {
            if self.seqs.is_empty() {
                0
            } else {
                (self.total_length() / self.num_seqs() as i64) as i32
            }
        }
        fn name(&self) -> &str {
            "tracking"
        }
        fn is_protein(&self) -> bool {
            self.is_protein
        }
        fn seq_len(&self, oid: i32) -> i32 {
            self.seqs[oid as usize].len() as i32
        }
        fn get_sequence(&self, arg: &GetSeqArg) -> Option<SeqData> {
            let seq = self.seqs.get(arg.oid as usize)?;
            Some(SeqData {
                sequence: seq.clone(),
                length: seq.len() as i32,
            })
        }
        fn iter_oids(&self) -> Box<dyn Iterator<Item = i32> + '_> {
            Box::new(0..self.num_seqs())
        }
        fn set_number_of_threads(&mut self, n_threads: i32) {
            self.n_threads = n_threads;
        }
        fn supports_partial_fetching(&self) -> bool {
            self.supports_partial
        }
        fn set_seq_ranges(&mut self, ranges: Option<&[SSeqRange]>) {
            self.ranges = ranges.map(|ranges| ranges.to_vec());
        }
        fn get_gis(&self, _oid: i32) -> Option<BlastGiList> {
            let mut gilist = blast_gi_list_new_ex(self.gis.len());
            for &gi in &self.gis {
                blast_gi_list_append(&mut gilist, gi);
            }
            Some(gilist)
        }
        fn iter_next(&self, itr: &mut BlastSeqSrcIterator) -> i32 {
            if itr.current_pos == u32::MAX {
                itr.oid_list = self.iter_oids().collect();
                itr.current_pos = 0;
            }
            let Some(&oid) = itr.oid_list.get(itr.current_pos as usize) else {
                return SEQSRC_EOF;
            };
            itr.current_pos += 1;
            oid
        }
        fn reset_chunk_iterator(&mut self) {
            self.reset_count += 1;
        }
    }

    impl BlastSeqSource for RecordingSeqSrc {
        fn num_seqs(&self) -> i32 {
            self.seqs.len() as i32
        }
        fn total_length(&self) -> i64 {
            self.seqs.iter().map(|s| s.len() as i64).sum()
        }
        fn max_seq_len(&self) -> i32 {
            self.seqs.iter().map(|s| s.len() as i32).max().unwrap_or(0)
        }
        fn avg_seq_len(&self) -> i32 {
            if self.seqs.is_empty() {
                0
            } else {
                (self.total_length() / self.num_seqs() as i64) as i32
            }
        }
        fn name(&self) -> &str {
            "recording"
        }
        fn is_protein(&self) -> bool {
            self.is_protein
        }
        fn seq_len(&self, oid: i32) -> i32 {
            self.seqs[oid as usize].len() as i32
        }
        fn get_sequence(&self, arg: &GetSeqArg) -> Option<SeqData> {
            *self.last_arg.lock().unwrap() = Some(arg.clone());
            let seq = self.seqs.get(arg.oid as usize)?;
            Some(SeqData {
                sequence: seq.clone(),
                length: seq.len() as i32,
            })
        }
        fn iter_oids(&self) -> Box<dyn Iterator<Item = i32> + '_> {
            Box::new(0..self.num_seqs())
        }
    }

    fn mock_seqsrc_constructor(_: Option<Box<dyn Any + Send>>) -> Option<Arc<dyn BlastSeqSource>> {
        Some(Arc::new(MockSeqSrc {
            seqs: vec![vec![1, 2, 3], vec![4, 5]],
        }))
    }

    fn mock_seqsrc_destructor(seq_src: &mut BlastSeqSrc) {
        seq_src.init_error = Some("freed".to_string());
    }

    fn mock_seqsrc_copier(seq_src: &BlastSeqSrc) -> Option<BlastSeqSrc> {
        let mut copy = seq_src.clone();
        copy.init_error = Some("copied".to_string());
        Some(copy)
    }

    #[test]
    fn test_blast_seq_src_new_copy_and_free_match_vtable_shape() {
        assert!(blast_seq_src_new(None).is_none());
        assert!(blast_seq_src_new(Some(BlastSeqSrcNewInfo {
            constructor: None,
            ctor_argument: None,
            delete_fn: None,
            copy_fn: None,
        }))
        .is_none());

        let seq_src = blast_seq_src_new(Some(BlastSeqSrcNewInfo {
            constructor: Some(mock_seqsrc_constructor),
            ctor_argument: None,
            delete_fn: Some(mock_seqsrc_destructor),
            copy_fn: Some(mock_seqsrc_copier),
        }))
        .expect("seqsrc");
        assert_eq!(seq_src.source.as_ref().unwrap().num_seqs(), 2);

        let copied = blast_seq_src_copy(Some(&seq_src)).expect("copy");
        assert_eq!(copied.init_error.as_deref(), Some("copied"));
        assert_eq!(copied.source.as_ref().unwrap().total_length(), 5);
        assert!(blast_seq_src_copy(None).is_none());
        assert!(blast_seq_src_free(Some(seq_src)).is_none());
    }

    #[test]
    fn test_mock_seqsrc() {
        let src = MockSeqSrc {
            seqs: vec![vec![0, 1, 2, 3], vec![0, 1]],
        };
        assert_eq!(src.num_seqs(), 2);
        assert_eq!(src.total_length(), 6);
        assert_eq!(src.max_seq_len(), 4);
        assert_eq!(src.seq_len(0), 4);
        assert_eq!(src.seq_len(1), 2);

        let data = src
            .get_sequence(&GetSeqArg {
                oid: 0,
                encoding: SeqEncoding::Protein,
                ..GetSeqArg::default()
            })
            .unwrap();
        assert_eq!(data.length, 4);
    }

    #[test]
    fn translated_seqsrc_getters_forward_to_trait_hooks() {
        let src = MockSeqSrc {
            seqs: vec![vec![0, 1, 2, 3], vec![0, 1]],
        };
        assert_eq!(blast_seq_src_get_num_seqs(&src), 2);
        assert_eq!(blast_seq_src_get_num_seqs_stats(&src), 2);
        assert_eq!(blast_seq_src_get_tot_len(&src), 6);
        assert_eq!(blast_seq_src_get_tot_len_stats(&src), 6);
        assert_eq!(blast_seq_src_get_max_seq_len(&src), 4);
        assert_eq!(blast_seq_src_get_min_seq_len(&src), BLAST_SEQSRC_MINLENGTH);
        assert_eq!(blast_seq_src_get_avg_seq_len(&src), 3);
        assert_eq!(blast_seq_src_get_name(&src), "mock");
        assert!(!blast_seq_src_get_is_prot(&src));
        assert_eq!(blast_seq_src_get_seq_len(&src, 1), 2);

        let mut arg = BlastSeqSrcGetSeqArg {
            oid: 0,
            encoding: EBlastEncoding::Protein,
            ..BlastSeqSrcGetSeqArg::default()
        };
        let data = blast_seq_src_get_sequence(Some(&src), Some(&mut arg)).expect("sequence");
        assert_eq!(data.sequence, vec![0, 1, 2, 3]);
        assert_eq!(arg.seq.as_ref().map(|seq| seq.length), Some(4));
        blast_seq_src_release_sequence(Some(&src), Some(&mut arg));
        assert!(arg.seq.is_some());
        assert!(arg.seq.as_ref().unwrap().sequence.is_none());
        assert!(blast_seq_src_get_sequence(None, Some(&mut arg)).is_none());
        assert!(blast_seq_src_get_sequence(Some(&src), None).is_none());
    }

    #[test]
    fn translated_seqsrc_get_sequence_blk_preserves_encoding_and_block_layout() {
        assert_eq!(SeqEncoding::from(EBlastEncoding::Error), SeqEncoding::Error);
        assert_eq!(
            EBlastEncoding::from(SeqEncoding::Error),
            EBlastEncoding::Error
        );

        let src = MockSeqSrc {
            seqs: vec![vec![99, 1, 2, 3], vec![4, 5, 6]],
        };
        let mut nucleotide_arg = BlastSeqSrcGetSeqArg {
            oid: 0,
            encoding: EBlastEncoding::Nucleotide,
            ..BlastSeqSrcGetSeqArg::default()
        };
        blast_seq_src_get_sequence(Some(&src), Some(&mut nucleotide_arg)).expect("sequence");
        let seq = nucleotide_arg.seq.as_ref().unwrap();
        assert_eq!(seq.sequence_start, Some(vec![99, 1, 2, 3]));
        assert_eq!(seq.sequence, Some(vec![1, 2, 3]));
        assert!(seq.sequence_start_allocated);
        assert!(!seq.sequence_allocated);

        let mut ncbi4na_arg = BlastSeqSrcGetSeqArg {
            oid: 1,
            encoding: EBlastEncoding::Ncbi4na,
            ..BlastSeqSrcGetSeqArg::default()
        };
        blast_seq_src_get_sequence(Some(&src), Some(&mut ncbi4na_arg)).expect("sequence");
        let seq = ncbi4na_arg.seq.as_ref().unwrap();
        assert_eq!(seq.sequence_start, Some(vec![4, 5, 6]));
        assert_eq!(seq.sequence, Some(vec![4, 5, 6]));
        assert!(seq.sequence_start_allocated);
        assert!(!seq.sequence_allocated);

        let mut protein_arg = BlastSeqSrcGetSeqArg {
            oid: 1,
            encoding: EBlastEncoding::Protein,
            ..BlastSeqSrcGetSeqArg::default()
        };
        blast_seq_src_get_sequence(Some(&src), Some(&mut protein_arg)).expect("sequence");
        let seq = protein_arg.seq.as_ref().unwrap();
        assert_eq!(seq.sequence_start, None);
        assert_eq!(seq.sequence, Some(vec![4, 5, 6]));
        assert!(!seq.sequence_start_allocated);
        assert!(seq.sequence_allocated);

        blast_seq_src_release_sequence(Some(&src), Some(&mut ncbi4na_arg));
        let seq = ncbi4na_arg.seq.as_ref().unwrap();
        assert!(seq.sequence_start.is_none());
        assert!(seq.sequence.is_none());
        assert!(!seq.sequence_start_allocated);
        assert!(ncbi4na_arg.ranges.is_none());
    }

    #[test]
    fn translated_seqsrc_get_sequence_blk_reset_ranges_before_fetch_for_nucleotide_sources() {
        let src = RecordingSeqSrc {
            seqs: vec![vec![0, 1, 2]],
            is_protein: false,
            last_arg: Mutex::new(None),
        };
        let mut arg = BlastSeqSrcGetSeqArg {
            oid: 0,
            encoding: EBlastEncoding::Nucleotide,
            reset_ranges: true,
            ranges: Some(BlastSeqSrcSetRangesArg {
                oid: 0,
                capacity: 1,
                num_ranges: 1,
                ranges: vec![0, 2],
            }),
            ..BlastSeqSrcGetSeqArg::default()
        };
        blast_seq_src_get_sequence(Some(&src), Some(&mut arg)).expect("sequence");
        assert!(arg.ranges.is_none());
        assert!(src
            .last_arg
            .lock()
            .unwrap()
            .as_ref()
            .unwrap()
            .ranges
            .is_none());

        let protein_src = RecordingSeqSrc {
            seqs: vec![vec![0, 1, 2]],
            is_protein: true,
            last_arg: Mutex::new(None),
        };
        let mut protein_arg = BlastSeqSrcGetSeqArg {
            oid: 0,
            encoding: EBlastEncoding::Protein,
            reset_ranges: true,
            ranges: Some(BlastSeqSrcSetRangesArg {
                oid: 0,
                capacity: 1,
                num_ranges: 1,
                ranges: vec![0, 2],
            }),
            ..BlastSeqSrcGetSeqArg::default()
        };
        blast_seq_src_get_sequence(Some(&protein_src), Some(&mut protein_arg)).expect("sequence");
        assert!(protein_arg.ranges.is_some());
        assert!(protein_src
            .last_arg
            .lock()
            .unwrap()
            .as_ref()
            .unwrap()
            .ranges
            .is_some());
    }

    fn test_hsp(
        subject_offset: i32,
        subject_end: i32,
        subject_frame: i32,
    ) -> crate::hspstream::Hsp {
        crate::hspstream::Hsp {
            score: 10,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 0.0,
            query_offset: 0,
            query_end: 0,
            query_gapped_start: 0,
            subject_offset,
            subject_end,
            subject_gapped_start: subject_offset,
            context: 0,
            query_frame: 0,
            subject_frame,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        }
    }

    #[test]
    fn setup_partial_fetching_builds_and_merges_subject_ranges() {
        let src = MockSeqSrc {
            seqs: vec![vec![0; 5000]],
        };
        let mut list = crate::hspstream::blast_hsp_list_new(0);
        list.oid = 0;
        let _ = crate::hspstream::blast_hsp_list_save_hsp(&mut list, test_hsp(1100, 1200, 0));
        let _ = crate::hspstream::blast_hsp_list_save_hsp(&mut list, test_hsp(1250, 1300, 0));

        let ranges =
            blast_setup_partial_fetching(crate::program::BLASTN, &src, &[&list]).expect("ranges");
        assert_eq!(ranges.oid, 0);
        assert_eq!(ranges.num_ranges, 1);
        assert_eq!(ranges.ranges, vec![76, 2324]);
    }

    #[test]
    fn setup_partial_fetching_translated_negative_frame_uses_nucleotide_reverse_range() {
        let src = MockSeqSrc {
            seqs: vec![vec![0; 9000]],
        };
        let mut list = crate::hspstream::blast_hsp_list_new(0);
        list.oid = 0;
        let _ = crate::hspstream::blast_hsp_list_save_hsp(&mut list, test_hsp(10, 20, -1));

        let ranges =
            blast_setup_partial_fetching(crate::program::TBLASTN, &src, &[&list]).expect("ranges");
        assert_eq!(ranges.num_ranges, 1);
        assert_eq!(ranges.ranges, vec![7910, 9000]);
    }

    /// Create a sequence source, iterate all sequences via iter_oids, verify count matches num_seqs.
    #[test]
    fn test_seqsrc_iteration() {
        let src = MockSeqSrc {
            seqs: vec![
                vec![0, 1, 2, 3],
                vec![4, 5],
                vec![6, 7, 8],
                vec![9, 10, 11, 12, 13],
            ],
        };
        let oids: Vec<i32> = src.iter_oids().collect();
        assert_eq!(oids.len() as i32, src.num_seqs());
        assert_eq!(oids, vec![0, 1, 2, 3]);

        // Verify each OID yields valid sequence data
        for oid in &oids {
            let data = src.get_sequence(&GetSeqArg {
                oid: *oid,
                encoding: SeqEncoding::Nucleotide,
                ..GetSeqArg::default()
            });
            assert!(data.is_some(), "Failed to get sequence for OID {}", oid);
        }

        // Empty source should iterate zero times
        let empty_src = MockSeqSrc { seqs: vec![] };
        let empty_oids: Vec<i32> = empty_src.iter_oids().collect();
        assert_eq!(empty_oids.len(), 0);
        assert_eq!(empty_src.num_seqs(), 0);
    }

    /// Verify sequence lengths from the source match the actual data lengths.
    #[test]
    fn test_seqsrc_length_queries() {
        let seqs = vec![
            vec![0u8; 100],
            vec![0u8; 1],
            vec![0u8; 50],
            vec![0u8; 200],
            vec![0u8; 75],
        ];
        let src = MockSeqSrc { seqs };

        assert_eq!(src.seq_len(0), 100);
        assert_eq!(src.seq_len(1), 1);
        assert_eq!(src.seq_len(2), 50);
        assert_eq!(src.seq_len(3), 200);
        assert_eq!(src.seq_len(4), 75);
        assert_eq!(src.max_seq_len(), 200);
        assert_eq!(src.min_seq_len(), BLAST_SEQSRC_MINLENGTH);
        // avg = (100+1+50+200+75)/5 = 85 (integer division of 426/5 = 85)
        assert_eq!(src.avg_seq_len(), 85);
    }

    /// Out-of-bounds OID access should return None from get_sequence.
    #[test]
    fn test_seqsrc_bounds_checking() {
        let src = MockSeqSrc {
            seqs: vec![vec![0, 1, 2], vec![3, 4]],
        };
        // Valid OIDs
        assert!(src
            .get_sequence(&GetSeqArg {
                oid: 0,
                encoding: SeqEncoding::Nucleotide,
                ..GetSeqArg::default()
            })
            .is_some());
        assert!(src
            .get_sequence(&GetSeqArg {
                oid: 1,
                encoding: SeqEncoding::Nucleotide,
                ..GetSeqArg::default()
            })
            .is_some());
        // Out of bounds
        assert!(src
            .get_sequence(&GetSeqArg {
                oid: 2,
                encoding: SeqEncoding::Nucleotide,
                ..GetSeqArg::default()
            })
            .is_none());
        assert!(src
            .get_sequence(&GetSeqArg {
                oid: 99,
                encoding: SeqEncoding::Nucleotide,
                ..GetSeqArg::default()
            })
            .is_none());
        // Negative OID (wraps to large usize, so should be None)
        assert!(src
            .get_sequence(&GetSeqArg {
                oid: -1,
                encoding: SeqEncoding::Nucleotide,
                ..GetSeqArg::default()
            })
            .is_none());
    }

    /// Sum of all sequence lengths should match total_length.
    #[test]
    fn test_seqsrc_total_length() {
        let src = MockSeqSrc {
            seqs: vec![vec![0; 10], vec![0; 20], vec![0; 30], vec![0; 40]],
        };
        // Total should be 10 + 20 + 30 + 40 = 100
        assert_eq!(src.total_length(), 100);

        // Verify by summing individual lengths
        let sum: i64 = src.iter_oids().map(|oid| src.seq_len(oid) as i64).sum();
        assert_eq!(sum, src.total_length());

        // Empty source should have total length 0
        let empty = MockSeqSrc { seqs: vec![] };
        assert_eq!(empty.total_length(), 0);
    }

    #[test]
    fn translated_seqsrc_vtable_wrappers_forward_to_trait_methods() {
        let mut src = TrackingSeqSrc {
            seqs: vec![vec![0; 10], vec![0; 20]],
            is_protein: true,
            supports_partial: true,
            n_threads: 1,
            ranges: None,
            gis: vec![42, 77],
            reset_count: 0,
        };

        assert_eq!(blast_seq_src_get_tot_len(&src), 30);
        assert!(blast_seq_src_get_is_prot(&src));
        assert!(blast_seq_src_get_supports_partial_fetching(&src));

        blast_seq_src_set_number_of_threads(Some(&mut src), 8);
        assert_eq!(src.n_threads, 8);

        let ranges = [
            SSeqRange { left: 2, right: 9 },
            SSeqRange {
                left: 14,
                right: 17,
            },
        ];
        blast_seq_src_set_seq_ranges(&mut src, Some(&ranges));
        assert_eq!(src.ranges, Some(ranges.to_vec()));

        blast_seq_src_set_number_of_threads(None, 99);

        let gis = blast_seq_src_get_gis(&src, 0).unwrap();
        assert_eq!(gis.data, vec![42, 77]);
        assert_eq!(gis.num_used, 2);

        let mut itr = blast_seq_src_iterator_new_ex(2);
        assert_eq!(itr.chunk_sz, 2);
        assert_eq!(itr.current_pos, u32::MAX);
        assert_eq!(blast_seq_src_iterator_next(&src, &mut itr), 0);
        assert_eq!(blast_seq_src_iterator_next(&src, &mut itr), 1);
        assert_eq!(blast_seq_src_iterator_next(&src, &mut itr), SEQSRC_EOF);

        blast_seq_src_reset_chunk_iterator(&mut src);
        assert_eq!(src.reset_count, 1);
    }

    #[test]
    fn translated_seqsrc_optional_hooks_default_to_c_absent_callback_shape() {
        let mut src = MockSeqSrc {
            seqs: vec![vec![0; 3]],
        };

        assert!(!blast_seq_src_get_supports_partial_fetching(&src));
        blast_seq_src_set_number_of_threads(Some(&mut src), 4);
        blast_seq_src_set_seq_ranges(&mut src, None);
        assert_eq!(blast_seq_src_get_tot_len(&src), 3);

        let mut itr = blast_seq_src_iterator_new();
        assert_eq!(itr.chunk_sz, K_BLAST_SEQ_SRC_DEFAULT_CHUNK_SIZE);
        assert_eq!(itr.current_pos, u32::MAX);
        assert_eq!(blast_seq_src_iterator_next(&src, &mut itr), 0);
        assert_eq!(blast_seq_src_iterator_next(&src, &mut itr), SEQSRC_EOF);
    }

    #[test]
    fn translated_set_ranges_arg_adds_overhang_and_merges_sorted_ranges() {
        let mut arg = blast_seq_src_set_ranges_arg_new(2);

        assert_eq!(
            blast_seq_src_set_ranges_arg_add_range(&mut arg, 3000, 3200, 10_000),
            0
        );
        assert_eq!(
            blast_seq_src_set_ranges_arg_add_range(&mut arg, 100, 200, 10_000),
            0
        );
        assert_eq!(
            blast_seq_src_set_ranges_arg_add_range(&mut arg, 4300, 4500, 10_000),
            0
        );
        assert!(arg.capacity >= 4);
        assert_eq!(arg.num_ranges, 6);
        assert_eq!(arg.ranges, vec![1976, 4224, 0, 1224, 3276, 5524]);

        blast_seq_src_set_ranges_arg_build(&mut arg);
        assert_eq!(arg.num_ranges, 1);
        assert_eq!(arg.ranges, vec![0, 5524]);
    }

    #[test]
    fn translated_set_ranges_arg_keeps_large_gaps_and_free_returns_none() {
        let mut arg = blast_seq_src_set_ranges_arg_new(4);

        blast_seq_src_set_ranges_arg_add_range(&mut arg, 10, 20, 20_000);
        blast_seq_src_set_ranges_arg_add_range(&mut arg, 5_000, 5_100, 20_000);
        blast_seq_src_set_ranges_arg_build(&mut arg);

        assert_eq!(arg.num_ranges, 2);
        assert_eq!(arg.ranges, vec![0, 1044, 3976, 6124]);

        let mut maybe_arg = Some(arg);
        assert_eq!(blast_seq_src_set_ranges_arg_free(&mut maybe_arg), None);
        assert_eq!(maybe_arg, None);
    }

    #[test]
    fn translated_gi_list_lifecycle_reallocates_and_appends() {
        let mut gilist = blast_gi_list_new_ex(1);

        assert_eq!(blast_gi_list_append(&mut gilist, 10), 0);
        assert_eq!(blast_gi_list_append(&mut gilist, 20), 0);
        assert_eq!(gilist.data, vec![10, 20]);
        assert_eq!(gilist.num_used, 2);
        assert!(gilist.num_allocated >= 2);

        let mut default_list = blast_gi_list_new();
        assert_eq!(default_list.num_allocated, 10);
        assert_eq!(blast_gi_list_realloc_if_necessary(&mut default_list), 0);
        assert_eq!(default_list.num_allocated, 10);

        let mut maybe_list = Some(gilist);
        assert_eq!(blast_gi_list_free(&mut maybe_list), None);
        assert_eq!(maybe_list, None);
    }

    #[test]
    fn translated_seqsrc_iterator_lifecycle_uses_default_chunk_and_frees() {
        let itr = blast_seq_src_iterator_new();
        assert_eq!(itr.chunk_sz, K_BLAST_SEQ_SRC_DEFAULT_CHUNK_SIZE);
        assert_eq!(itr.current_pos, u32::MAX);
        assert!(itr.oid_list.capacity() >= K_BLAST_SEQ_SRC_DEFAULT_CHUNK_SIZE as usize);

        let mut maybe_itr = Some(itr);
        assert_eq!(blast_seq_src_iterator_free(&mut maybe_itr), None);
        assert_eq!(maybe_itr, None);
    }
}
