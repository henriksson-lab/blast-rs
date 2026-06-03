//! Rust equivalent of blast_lookup.c, blast_nalookup.c, lookup_wrap.c
//! Lookup table structures for BLAST word finding.

/// Lookup table types. Variant names mirror NCBI `ELookupTableType`
/// (`lookup_wrap.h`) so they can be line-diffed against the C enum.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(clippy::enum_variant_names)]
pub enum LookupTableType {
    MegablastLookup,
    SmallNaLookup,
    NaLookup,
    AaLookup,
    CompressedAaLookup,
    RpsLookup,
    PhiLookup,
    PhiNaLookup,
    IndexedMbLookup,
    NaHashLookup,
}

/// General types of discontiguous word templates (`EDiscWordType`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DiscWordType {
    Coding = 0,
    Optimal = 1,
    TwoTemplates = 2,
}

/// Enumeration of all discontiguous word templates (`EDiscTemplateType`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DiscTemplateType {
    Contiguous = 0,
    Template11_16Coding = 1,
    Template11_16Optimal = 2,
    Template12_16Coding = 3,
    Template12_16Optimal = 4,
    Template11_18Coding = 5,
    Template11_18Optimal = 6,
    Template12_18Coding = 7,
    Template12_18Optimal = 8,
    Template11_21Coding = 9,
    Template11_21Optimal = 10,
    Template12_21Coding = 11,
    Template12_21Optimal = 12,
}

impl Default for DiscTemplateType {
    fn default() -> Self {
        Self::Contiguous
    }
}

/// A (query_offset, subject_offset) pair from the scan phase.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct OffsetPair {
    pub query_offset: i32,
    pub subject_offset: i32,
}

/// Seed pair accepted by the conservative RPS two-hit word finder.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RpsTwoHitSeed {
    pub first: OffsetPair,
    pub second: OffsetPair,
}

/// Rust-owned result for the represented RPS two-hit word-finder side effects.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct RpsWordFinderScan {
    pub seeds: Vec<RpsTwoHitSeed>,
    pub total_hits: i32,
    pub hits_extended: i32,
    init_hits: Vec<RpsInitHitRecord>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct RpsInitHitRecord {
    q_start: i32,
    s_start: i32,
    q_off: i32,
    s_off: i32,
    len: i32,
    score: i32,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct RpsTwoHitDiagState {
    last_hit: i32,
    flag: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum RpsTwoHitOutcome {
    NoReach {
        init_hit: Option<RpsInitHitRecord>,
    },
    Reached {
        init_hit: Option<RpsInitHitRecord>,
        s_last_off: i32,
    },
}

/// blast-rs: RPS scan/HSP bridge for the represented ungapped payload path;
/// not a direct NCBI C port.
///
/// The RPS traceback driver receives preliminary HSPs in the swapped
/// profile/query orientation used by NCBI before `s_BlastHSPListRPSUpdate`
/// restores the final user-visible coordinates. This helper preserves that
/// orientation while converting owned RPS word-finder payloads into the local
/// HSP-list representation.
pub fn blast_rps_word_finder_scan_to_hsp_list(
    scan: &RpsWordFinderScan,
    profile_oid: i32,
    hsp_max: i32,
) -> Option<crate::hspstream::HspList> {
    blast_rps_word_finder_scan_to_hsp_list_with_query_shift(scan, profile_oid, hsp_max, 0)
}

/// blast-rs: Driver-facing adapter from represented RPS word-finder payloads
/// into the traceback stream; not a direct NCBI C port.
///
/// Full RPS-BLAST database wiring is still outside this module, but this
/// preserves the boundary the driver needs after RPS ungapped extension: empty
/// scans are a no-op, and nonempty scans enter the regular HSP stream in the
/// pre-traceback RPS orientation consumed by `s_RPSComputeTraceback`.
pub fn blast_rps_word_finder_scan_write_hsp_stream(
    scan: &RpsWordFinderScan,
    stream: &crate::hspstream::HspStream,
    query_index: i32,
    profile_oid: i32,
    hsp_max: i32,
) -> i32 {
    let Some(hsp_list) = blast_rps_word_finder_scan_to_hsp_list(scan, profile_oid, hsp_max) else {
        return 0;
    };
    stream.blast_hspstream_write_blast_hsp_list(
        query_index,
        crate::hspstream::blast_hsp_list_from_legacy_hsp_list(hsp_list, query_index),
    )
}

/// blast-rs: bounded RPS driver adapter from an owned lookup table and subject
/// sequence into the traceback stream; not a direct NCBI C port.
///
/// The full CLI path still needs native CDD/RPS database loading, but this is
/// the driver boundary immediately below that reader: run the represented RPS
/// two-hit ungapped payload scan for one consensus/profile sequence and feed
/// the resulting preliminary HSPs into the same stream consumed by
/// `s_RPSComputeTraceback`.
pub fn blast_rps_subject_scan_write_hsp_stream(
    lookup: &mut BlastRpsLookupTable,
    subject: &[u8],
    stream: &crate::hspstream::HspStream,
    query_index: i32,
    profile_oid: i32,
    window_size: i32,
    x_dropoff: i32,
    score_cutoff: i32,
    hsp_max: i32,
    init_hitlist: Option<&mut crate::extend::InitHitList>,
    ungapped_stats: Option<&mut crate::diagnostics::UngappedStats>,
) -> i32 {
    let scan = s_blast_rps_word_finder_two_hit_with_extension_payloads(
        lookup,
        subject,
        window_size,
        x_dropoff,
        score_cutoff,
        init_hitlist,
        ungapped_stats,
    );
    blast_rps_word_finder_scan_write_hsp_stream(&scan, stream, query_index, profile_oid, hsp_max)
}

/// blast-rs: SeqSrc-backed RPS driver adapter into the traceback stream; not a
/// direct NCBI C port.
///
/// This is the owned Rust boundary corresponding to the native RPS driver
/// layer above `s_RPSComputeTraceback`: fetch each profile/consensus sequence
/// from `BlastSeqSrc` with the same encoding requested by traceback, run the
/// represented RPS ungapped payload scan, and write nonempty preliminary HSP
/// lists into the stream keyed by profile OID.
pub fn blast_rps_seqsrc_scan_write_hsp_stream(
    lookup: &mut BlastRpsLookupTable,
    seq_src: Option<&dyn crate::seqsrc::BlastSeqSource>,
    stream: &crate::hspstream::HspStream,
    program_number: crate::program::ProgramType,
    query_index: i32,
    window_size: i32,
    x_dropoff: i32,
    score_cutoff: i32,
    hsp_max: i32,
    mut init_hitlist: Option<&mut crate::extend::InitHitList>,
    mut ungapped_stats: Option<&mut crate::diagnostics::UngappedStats>,
) -> i32 {
    let Some(seq_src) = seq_src else {
        return -1;
    };
    let encoding = crate::hspstream::blast_traceback_get_encoding(program_number);
    let min_subject_len = lookup.wordsize.max(0);

    for oid in seq_src.iter_oids() {
        let Some(seq_data) = seq_src.get_sequence(&crate::seqsrc::GetSeqArg {
            oid,
            encoding,
            ..crate::seqsrc::GetSeqArg::default()
        }) else {
            continue;
        };
        if seq_data.length < min_subject_len || seq_data.sequence.len() < min_subject_len as usize {
            continue;
        }
        let status = blast_rps_subject_scan_write_hsp_stream(
            lookup,
            &seq_data.sequence,
            stream,
            query_index,
            oid,
            window_size,
            x_dropoff,
            score_cutoff,
            hsp_max,
            init_hitlist.as_deref_mut(),
            ungapped_stats.as_deref_mut(),
        );
        if status != 0 {
            return status;
        }
    }
    0
}

/// blast-rs: native owned RPS profile database scan into the traceback stream;
/// not a direct NCBI C port.
///
/// This is the driver boundary after owned profile database assembly: scan one
/// query sequence against the constructed RPS lookup, split represented
/// ungapped payloads by profile offsets, localize profile coordinates, and
/// write per-profile preliminary HSP lists for `s_RPSComputeTraceback`.
pub fn blast_rps_profile_database_scan_query_write_hsp_stream(
    profile_db: &mut OwnedRpsProfileDatabase,
    query: &[u8],
    stream: &crate::hspstream::HspStream,
    query_index: i32,
    window_size: i32,
    x_dropoff: i32,
    score_cutoff: i32,
    hsp_max: i32,
    init_hitlist: Option<&mut crate::extend::InitHitList>,
    ungapped_stats: Option<&mut crate::diagnostics::UngappedStats>,
) -> i32 {
    let scan = s_blast_rps_word_finder_two_hit_with_extension_payloads(
        &mut profile_db.lookup,
        query,
        window_size,
        x_dropoff,
        score_cutoff,
        init_hitlist,
        ungapped_stats,
    );
    if scan.init_hits.is_empty() {
        return 0;
    }

    let start_offsets = &profile_db.traceback_info.profile_header.start_offsets;
    let num_profiles = profile_db.traceback_info.profile_header.num_profiles.max(0) as usize;
    for profile_index in 0..num_profiles {
        let profile_start = start_offsets[profile_index];
        let profile_end = start_offsets[profile_index + 1];
        let mut profile_scan = RpsWordFinderScan {
            seeds: Vec::new(),
            total_hits: scan.total_hits,
            hits_extended: scan.hits_extended,
            init_hits: Vec::new(),
        };
        for hit in &scan.init_hits {
            let Some(hit_profile_index) = rps_profile_index_for_offset(start_offsets, hit.q_start)
            else {
                continue;
            };
            if hit_profile_index != profile_index || hit.q_start + hit.len > profile_end {
                continue;
            }
            profile_scan.init_hits.push(hit.clone());
        }
        if let Some(hsp_list) = blast_rps_word_finder_scan_to_hsp_list_with_query_shift(
            &profile_scan,
            profile_index as i32,
            hsp_max,
            profile_start,
        ) {
            let status = stream.blast_hspstream_write_blast_hsp_list(
                query_index,
                crate::hspstream::blast_hsp_list_from_legacy_hsp_list(hsp_list, query_index),
            );
            if status != 0 {
                return status;
            }
        }
    }
    0
}

/// Rust-owned equivalent of C `MapperWordHits`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MapperWordHits {
    pub pair_arrays: Vec<Vec<OffsetPair>>,
    pub num: Vec<i32>,
    pub num_arrays: i32,
    pub array_size: i32,
    pub divisor: i32,
    pub last_diag: Vec<i32>,
    pub last_pos: Vec<i32>,
}

/// Typed replacement for the nucleotide lookup callback pointer selected by
/// NCBI `BlastChooseNaExtend`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NaLookupCallback {
    Megablast,
    SmallNa,
    Na,
}

/// Typed replacement for the nucleotide extension callback pointer selected by
/// NCBI `BlastChooseNaExtend`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NaExtendCallback {
    Direct,
    Aligned,
    Generic,
    SmallAlignedOneByte,
    SmallGeneric,
}

/// Rust-owned view of the two callback fields C mutates on `LookupTableWrap`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct NaExtendChoice {
    pub lookup_callback: Option<NaLookupCallback>,
    pub extend_callback: Option<NaExtendCallback>,
}

pub const NA_HITS_PER_CELL: usize = 3;
pub const NA_WORDS_PER_HASH: usize = 3;
pub const NA_OFFSETS_PER_HASH: usize = 9;
const BLAST2NA_MASK: u8 = 0xfc;
const BITS_PER_NUC: u32 = 2;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NaLookupBackboneCell {
    pub num_used: i32,
    pub entries: [i32; NA_HITS_PER_CELL],
    pub overflow_cursor: i32,
}

impl Default for NaLookupBackboneCell {
    fn default() -> Self {
        Self {
            num_used: 0,
            entries: [0; NA_HITS_PER_CELL],
            overflow_cursor: 0,
        }
    }
}

/// Rust-owned representation of NCBI `BlastNaLookupTable` fields needed by
/// the standard nucleotide scanner helpers.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct BlastNaLookupTable {
    pub mask: i32,
    pub word_length: i32,
    pub lut_word_length: i32,
    pub scan_step: i32,
    pub backbone_size: i32,
    pub longest_chain: i32,
    pub thick_backbone: Vec<NaLookupBackboneCell>,
    pub overflow: Vec<i32>,
    pub overflow_size: i32,
    pub pv: Vec<u32>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NaHashLookupBackboneCell {
    pub num_words: i8,
    pub num_offsets: [i8; NA_WORDS_PER_HASH],
    pub words: [u32; NA_WORDS_PER_HASH],
    pub offsets: [i32; NA_OFFSETS_PER_HASH],
}

impl Default for NaHashLookupBackboneCell {
    fn default() -> Self {
        Self {
            num_words: 0,
            num_offsets: [0; NA_WORDS_PER_HASH],
            words: [0; NA_WORDS_PER_HASH],
            offsets: [0; NA_OFFSETS_PER_HASH],
        }
    }
}

/// Rust-owned representation of NCBI `BlastNaHashLookupTable`.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct BlastNaHashLookupTable {
    pub mask: i32,
    pub word_length: i32,
    pub lut_word_length: i32,
    pub scan_step: i32,
    pub backbone_size: i32,
    pub longest_chain: i32,
    pub thick_backbone: Vec<NaHashLookupBackboneCell>,
    pub overflow: Vec<i32>,
    pub offsets_size: i32,
    pub pv: Vec<u32>,
    pub pv_array_bts: i32,
}

/// Rust-owned equivalent of C `NaHashLookupThreadData`.
pub struct NaHashLookupThreadData {
    pub seq_arg: Vec<crate::seqsrc::GetSeqArg>,
    pub itr: Vec<crate::seqsrc::BlastSeqSrcIterator>,
    pub seq_src: Vec<crate::seqsrc::BlastSeqSrc>,
    pub word_counts: Vec<crate::util::BlastSparseUint1Array>,
    pub num_threads: i32,
}

#[derive(Debug, Clone)]
pub struct BackboneCell {
    pub word: u32,
    pub offset: i32,
    pub num_offsets: i32,
    pub next: Option<Box<BackboneCell>>,
}

impl Default for BackboneCell {
    fn default() -> Self {
        BackboneCell {
            word: 0,
            offset: 0,
            num_offsets: 0,
            next: None,
        }
    }
}

/// Nucleotide lookup table for standard blastn (word_size 4-12).
#[derive(Debug)]
pub struct SmallNaLookupTable {
    pub word_length: i32,
    pub backbone: Vec<i32>,
    pub overflow: Vec<i32>,
    pub pv_array: Vec<u32>, // presence vector for quick filtering
    pub longest_chain: i32,
    pub scan_step: i32,
}

/// Megablast lookup table (word_size >= 12, contiguous).
#[derive(Debug)]
pub struct MbLookupTable {
    pub word_length: i32,
    pub lut_word_length: i32, // may differ from word_length for discontiguous
    pub discontiguous: bool,
    pub template_length: i32,
    pub template_type: DiscTemplateType,
    pub two_templates: bool,
    pub second_template_type: DiscTemplateType,
    pub hashtable: Vec<i32>,
    pub hashtable2: Vec<i32>,
    pub next_pos: Vec<i32>,
    pub next_pos2: Vec<i32>,
    pub pv_array: Vec<u32>,
    pub pv_array_bts: i32,
    pub longest_chain: i32,
    pub scan_step: i32,
}

/// Protein lookup table.
#[derive(Debug)]
pub struct AaLookupTable {
    pub word_length: i32,
    pub threshold: f64,
    pub backbone: Vec<Vec<i32>>, // backbone[word_hash] = list of query offsets
    pub pv_array: Vec<u32>,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct RpsBackboneCell {
    pub num_used: i32,
    pub offset_pairs: Vec<OffsetPair>,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct RpsBucket {
    pub num_alloc: i32,
    pub num_used: i32,
    pub offset_pairs: Vec<OffsetPair>,
}

pub const RPS_BUCKET_SIZE: i32 = 2048;

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct BlastRpsInfo {
    pub alphabet_size: i32,
    pub wordsize: i32,
    pub rps_backbone: Vec<RpsBackboneCell>,
    pub rps_pssm: Vec<Vec<i32>>,
}

/// Rust-owned RPS profile database assembly boundary.
///
/// NCBI wires this through memory-mapped CDD/RPS files before lookup
/// construction. Rust keeps the same owned pieces together: traceback profile
/// metadata, consensus/profile sequences, the lookup construction input, and
/// the constructed RPS lookup table.
#[derive(Debug, Clone, PartialEq)]
pub struct OwnedRpsProfileDatabase {
    pub traceback_info: crate::hspstream::RpsTracebackInfo,
    pub consensus_sequences: Vec<Vec<u8>>,
    pub lookup_info: BlastRpsInfo,
    pub lookup: BlastRpsLookupTable,
}

/// File-backed inputs for the represented native RPS profile database boundary.
///
/// The `.rps` and optional `.freq` payloads use NCBI's native integer layout.
/// When `.aux` is supplied, its Karlin K values are wired into traceback
/// metadata; otherwise `karlin_k` is used as an already-decoded fallback.
/// Consensus sequences are still explicit because this boundary intentionally
/// avoids reaching into the ordinary BLAST database readers.
#[derive(Debug, Clone, PartialEq)]
pub struct NativeRpsProfileBundle {
    pub profile_bytes: Vec<u8>,
    pub freq_ratios_bytes: Option<Vec<u8>>,
    pub lookup_bytes: Option<Vec<u8>>,
    pub aux_bytes: Option<Vec<u8>>,
    pub consensus_sequences: Vec<Vec<u8>>,
    pub karlin_k: Vec<f64>,
    pub wordsize: i32,
}

/// Rust-owned representation of NCBI `BlastRPSLookupTable`.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct BlastRpsLookupTable {
    pub alphabet_size: i32,
    pub wordsize: i32,
    pub charsize: i32,
    pub backbone_size: i32,
    pub rps_backbone: Vec<RpsBackboneCell>,
    pub pv: Vec<u32>,
    pub rps_pssm: Vec<Vec<i32>>,
    pub num_buckets: i32,
    pub bucket_array: Vec<RpsBucket>,
}

const RPS_HITS_PER_CELL: usize = 3;
const RPS_LOOKUP_HEADER_WORDS: usize = 10;
const RPS_LOOKUP_BACKBONE_CELL_WORDS: usize = 1 + RPS_HITS_PER_CELL;
const RPS_LOOKUP_WORDSIZE: i32 = 3;

/// Generic wrapper around different lookup table types.
pub enum LookupTableWrap {
    Na(BlastNaLookupTable),
    NaHash(BlastNaHashLookupTable),
    SmallNa(SmallNaLookupTable),
    Megablast(MbLookupTable),
    Aa(AaLookupTable),
    Rps(BlastRpsLookupTable),
}

impl LookupTableWrap {
    pub fn table_type(&self) -> LookupTableType {
        match self {
            LookupTableWrap::Na(_) => LookupTableType::NaLookup,
            LookupTableWrap::NaHash(_) => LookupTableType::NaHashLookup,
            LookupTableWrap::SmallNa(_) => LookupTableType::SmallNaLookup,
            LookupTableWrap::Megablast(_) => LookupTableType::MegablastLookup,
            LookupTableWrap::Aa(_) => LookupTableType::AaLookup,
            LookupTableWrap::Rps(_) => LookupTableType::RpsLookup,
        }
    }

    pub fn word_length(&self) -> i32 {
        match self {
            LookupTableWrap::Na(t) => t.word_length,
            LookupTableWrap::NaHash(t) => t.word_length,
            LookupTableWrap::SmallNa(t) => t.word_length,
            LookupTableWrap::Megablast(t) => t.word_length,
            LookupTableWrap::Aa(t) => t.word_length,
            LookupTableWrap::Rps(t) => t.wordsize,
        }
    }
}

/// Port-shaped helper for NCBI static `s_GetMinimumSubjSeqLen`.
pub fn s_get_minimum_subj_seq_len(lookup_wrap: Option<&LookupTableWrap>) -> i32 {
    let Some(lookup_wrap) = lookup_wrap else {
        return 0;
    };
    match lookup_wrap {
        LookupTableWrap::Na(table) => table.word_length,
        LookupTableWrap::NaHash(table) => table.word_length,
        LookupTableWrap::SmallNa(table) => table.word_length,
        LookupTableWrap::Megablast(table) => table.word_length,
        LookupTableWrap::Aa(table) => table.word_length,
        LookupTableWrap::Rps(table) => table.wordsize,
    }
}

/// Port of NCBI inline `s_BlastLookupGetNumHits` (`blast_nascan.c:41`).
pub fn s_blast_lookup_get_num_hits(lookup: &BlastNaLookupTable, index: i32) -> i32 {
    if index < 0 {
        return 0;
    }
    let index = index as u32;
    let word = (index >> crate::stat::PV_ARRAY_BTS) as usize;
    let bit = index & crate::stat::PV_ARRAY_MASK;
    if lookup
        .pv
        .get(word)
        .is_some_and(|value| (value & (1u32 << bit)) != 0)
    {
        lookup
            .thick_backbone
            .get(index as usize)
            .map_or(0, |cell| cell.num_used)
    } else {
        0
    }
}

/// Port of NCBI inline `s_BlastLookupRetrieve` (`blast_nascan.c:58`).
pub fn s_blast_lookup_retrieve(
    lookup: &BlastNaLookupTable,
    index: i32,
    offset_pairs: &mut Vec<OffsetPair>,
    s_off: i32,
) {
    if index < 0 {
        return;
    }
    let Some(cell) = lookup.thick_backbone.get(index as usize) else {
        return;
    };
    let num_hits = cell.num_used.max(0) as usize;
    let lookup_pos: Vec<i32> = if num_hits <= NA_HITS_PER_CELL {
        cell.entries[..num_hits].to_vec()
    } else {
        let start = cell.overflow_cursor.max(0) as usize;
        lookup
            .overflow
            .get(start..start.saturating_add(num_hits))
            .unwrap_or(&[])
            .to_vec()
    };

    for q_off in lookup_pos {
        offset_pairs.push(OffsetPair {
            query_offset: q_off,
            subject_offset: s_off,
        });
    }
}

/// Port of NCBI inline `s_BlastMBLookupHasHits` (`blast_nascan.c:1393`).
pub fn s_blast_mb_lookup_has_hits(lookup: &MbLookupTable, index: i64) -> i32 {
    if index < 0 {
        return 0;
    }
    let index = index as u32;
    let pv_array_bts = lookup.pv_array_bts.max(0) as usize;
    let word = (index >> pv_array_bts) as usize;
    let bit = index & crate::stat::PV_ARRAY_MASK;
    if lookup
        .pv_array
        .get(word)
        .is_some_and(|value| (value & (1u32 << bit)) != 0)
    {
        1
    } else {
        0
    }
}

/// NCBI: s_BlastMBLookupRetrieve (blast_nascan.c:1413).
pub fn s_blast_mb_lookup_retrieve(
    lookup: &MbLookupTable,
    index: i64,
    offset_pairs: &mut Vec<OffsetPair>,
    s_off: i32,
) -> i32 {
    if index < 0 {
        return 0;
    }
    let mut q_off = lookup.hashtable.get(index as usize).copied().unwrap_or(0);
    let mut count = 0;
    while q_off != 0 {
        offset_pairs.push(OffsetPair {
            query_offset: q_off - 1,
            subject_offset: s_off,
        });
        count += 1;
        q_off = lookup.next_pos.get(q_off as usize).copied().unwrap_or(0);
    }
    count
}

/// NCBI: s_BlastMBLookupRetrieve2 (blast_nascan.c:1437).
pub fn s_blast_mb_lookup_retrieve2(
    lookup: &MbLookupTable,
    index: i64,
    offset_pairs: &mut Vec<OffsetPair>,
    s_off: i32,
) -> i32 {
    if index < 0 {
        return 0;
    }
    let mut q_off = lookup.hashtable2.get(index as usize).copied().unwrap_or(0);
    let mut count = 0;
    while q_off != 0 {
        offset_pairs.push(OffsetPair {
            query_offset: q_off - 1,
            subject_offset: s_off,
        });
        count += 1;
        q_off = lookup.next_pos2.get(q_off as usize).copied().unwrap_or(0);
    }
    count
}

/// Port of NCBI static `s_MBLookup` (`na_ungapped.c:49`).
pub fn s_mb_lookup(lookup_wrap: Option<&LookupTableWrap>, index: i32, q_pos: i32) -> bool {
    let Some(LookupTableWrap::Megablast(mb_lt)) = lookup_wrap else {
        return false;
    };
    if mb_lt.hashtable.is_empty() {
        return false;
    }
    let index = (index as usize) & (mb_lt.hashtable.len() - 1);
    let q_pos = q_pos + 1;
    if s_blast_mb_lookup_has_hits(mb_lt, index as i64) == 0 {
        return false;
    }

    let mut q_off = mb_lt.hashtable[index];
    while q_off != 0 {
        if q_off == q_pos {
            return true;
        }
        q_off = mb_lt.next_pos.get(q_off as usize).copied().unwrap_or(0);
    }
    false
}

/// Port of NCBI static `s_SmallNaLookup` (`na_ungapped.c:80`).
pub fn s_small_na_lookup(lookup_wrap: Option<&LookupTableWrap>, index: i32, q_pos: i32) -> bool {
    let Some(LookupTableWrap::SmallNa(lookup)) = lookup_wrap else {
        return false;
    };
    if lookup.backbone.is_empty() {
        return false;
    }
    let index = lookup.backbone[(index as usize) & (lookup.backbone.len() - 1)];

    if index == q_pos {
        return true;
    }
    if index == -1 || index >= 0 {
        return false;
    }

    let mut src_off = (-index) as usize;
    loop {
        let Some(&lookup_pos) = lookup.overflow.get(src_off) else {
            return false;
        };
        src_off += 1;
        if lookup_pos == q_pos {
            return true;
        }
        if lookup_pos < 0 {
            return false;
        }
    }
}

/// Retrieve all query offsets for a SmallNa lookup word, preserving C's compact
/// backbone/overflow layout used by `BlastSmallNaLookupTable`.
pub fn s_blast_small_na_lookup_retrieve(
    lookup: &SmallNaLookupTable,
    index: i32,
    offset_pairs: &mut Vec<OffsetPair>,
    s_off: i32,
) -> i32 {
    if lookup.backbone.is_empty() {
        return 0;
    }

    let index = lookup.backbone[(index as usize) & (lookup.backbone.len() - 1)];
    let before = offset_pairs.len();
    if index >= 0 {
        offset_pairs.push(OffsetPair {
            query_offset: index,
            subject_offset: s_off,
        });
    } else if index != -1 {
        let mut src_off = (-index) as usize;
        loop {
            let Some(&lookup_pos) = lookup.overflow.get(src_off) else {
                break;
            };
            src_off += 1;
            if lookup_pos < 0 {
                break;
            }
            offset_pairs.push(OffsetPair {
                query_offset: lookup_pos,
                subject_offset: s_off,
            });
        }
    }

    (offset_pairs.len() - before) as i32
}

/// Infer `BlastSmallNaLookupTable::lut_word_length` from the finalized
/// backbone size. Rust's compact SmallNa owner stores the same information as
/// C's `backbone_size = 1 << (BITS_PER_NUC * lut_word_length)`.
pub fn small_na_lut_word_length(lookup: &SmallNaLookupTable) -> i32 {
    let len = lookup.backbone.len();
    if len > 0 && len.is_power_of_two() {
        let bits = len.trailing_zeros();
        if bits % BITS_PER_NUC == 0 {
            let width = (bits / BITS_PER_NUC) as i32;
            if width > 0 {
                return width;
            }
        }
    }
    lookup.word_length
}

/// Port of NCBI static `s_NaLookup` (`na_ungapped.c:111`) over Rust's
/// `BlastNaLookupTable` owner.
pub fn s_na_lookup(lookup: Option<&BlastNaLookupTable>, index: i32, q_pos: i32) -> bool {
    let Some(lookup) = lookup else {
        return false;
    };
    if lookup.thick_backbone.is_empty() {
        return false;
    }
    let index = (index & lookup.mask) as usize;
    if s_blast_lookup_get_num_hits(lookup, index as i32) == 0 {
        return false;
    }
    let Some(cell) = lookup.thick_backbone.get(index) else {
        return false;
    };
    let num_hits = cell.num_used.max(0) as usize;
    let positions: Vec<i32> = if num_hits <= NA_HITS_PER_CELL {
        cell.entries.iter().take(num_hits).copied().collect()
    } else {
        let start = cell.overflow_cursor.max(0) as usize;
        lookup
            .overflow
            .get(start..start.saturating_add(num_hits))
            .unwrap_or(&[])
            .to_vec()
    };
    positions.into_iter().any(|pos| pos == q_pos)
}

/// Typed callback selection for NCBI `BlastChooseNaExtend`.
pub fn blast_choose_na_extend(lookup_wrap: Option<&LookupTableWrap>) -> Option<NaExtendChoice> {
    let lookup_wrap = lookup_wrap?;
    match lookup_wrap {
        LookupTableWrap::Megablast(lut) => {
            let extend_callback = if lut.lut_word_length == lut.word_length {
                NaExtendCallback::Direct
            } else if lut.lut_word_length % crate::stat::COMPRESSION_RATIO as i32 == 0
                && lut.scan_step % crate::stat::COMPRESSION_RATIO as i32 == 0
            {
                NaExtendCallback::Aligned
            } else {
                NaExtendCallback::Generic
            };
            Some(NaExtendChoice {
                lookup_callback: Some(NaLookupCallback::Megablast),
                extend_callback: Some(extend_callback),
            })
        }
        LookupTableWrap::SmallNa(lut) => {
            let lut_word_length = small_na_lut_word_length(lut);
            let extend_callback = if lut_word_length == lut.word_length {
                NaExtendCallback::Direct
            } else if lut_word_length % crate::stat::COMPRESSION_RATIO as i32 == 0
                && lut.scan_step % crate::stat::COMPRESSION_RATIO as i32 == 0
                && lut.word_length - lut_word_length <= 4
            {
                NaExtendCallback::SmallAlignedOneByte
            } else {
                NaExtendCallback::SmallGeneric
            };
            Some(NaExtendChoice {
                lookup_callback: Some(NaLookupCallback::SmallNa),
                extend_callback: Some(extend_callback),
            })
        }
        LookupTableWrap::Na(lut) => {
            let extend_callback = if lut.lut_word_length == lut.word_length {
                NaExtendCallback::Direct
            } else if lut.lut_word_length % crate::stat::COMPRESSION_RATIO as i32 == 0
                && lut.scan_step % crate::stat::COMPRESSION_RATIO as i32 == 0
            {
                NaExtendCallback::Aligned
            } else {
                NaExtendCallback::Generic
            };
            Some(NaExtendChoice {
                lookup_callback: Some(NaLookupCallback::Na),
                extend_callback: Some(extend_callback),
            })
        }
        LookupTableWrap::NaHash(_) => Some(NaExtendChoice {
            lookup_callback: None,
            extend_callback: None,
        }),
        LookupTableWrap::Aa(_) | LookupTableWrap::Rps(_) => None,
    }
}

/// Port of NCBI static `s_IsSeedMasked` (`na_ungapped.c:460`).
pub fn s_is_seed_masked(
    lookup_wrap: Option<&LookupTableWrap>,
    subject_packed: &[u8],
    s_off: i32,
    lut_word_length: i32,
    q_pos: i32,
) -> Option<bool> {
    let lookup_wrap = lookup_wrap?;
    if s_off < 0 || lut_word_length < 0 {
        return None;
    }
    let byte_offset = (s_off / 4) as usize;
    let shift = 2 * (16 - s_off % 4 - lut_word_length);
    if !(0..32).contains(&shift) {
        return None;
    }

    let b0 = u32::from(*subject_packed.get(byte_offset)?);
    let b1 = u32::from(*subject_packed.get(byte_offset + 1).unwrap_or(&0));
    let b2 = u32::from(*subject_packed.get(byte_offset + 2).unwrap_or(&0));
    let b3 = u32::from(*subject_packed.get(byte_offset + 3).unwrap_or(&0));
    let index = match shift {
        8 | 10 | 12 | 14 => ((b0 << 24) | (b1 << 16) | (b2 << 8)) >> shift,
        16 | 18 | 20 | 22 => ((b0 << 24) | (b1 << 16)) >> shift,
        24 => b0,
        _ => ((b0 << 24) | (b1 << 16) | (b2 << 8) | b3) >> shift,
    } as i32;

    let present = match blast_choose_na_extend(Some(lookup_wrap))?.lookup_callback? {
        NaLookupCallback::Megablast => s_mb_lookup(Some(lookup_wrap), index, q_pos),
        NaLookupCallback::SmallNa => s_small_na_lookup(Some(lookup_wrap), index, q_pos),
        NaLookupCallback::Na => {
            let LookupTableWrap::Na(lookup) = lookup_wrap else {
                return None;
            };
            s_na_lookup(Some(lookup), index, q_pos)
        }
    };
    Some(!present)
}

/// Port of NCBI static `s_TypeOfWord` (`na_ungapped.c:508`).
#[allow(clippy::too_many_arguments)]
pub fn s_type_of_word(
    subject_packed: &[u8],
    q_off: &mut i32,
    s_off: &mut i32,
    locations: Option<&[crate::util::SSeqRange]>,
    query_info: &crate::queryinfo::QueryInfo,
    s_range: u32,
    word_length: u32,
    lut_word_length: u32,
    lookup_wrap: Option<&LookupTableWrap>,
    check_double: bool,
    extended: &mut i32,
) -> i32 {
    *extended = 0;
    if word_length == lut_word_length {
        return 1;
    }
    if word_length < lut_word_length {
        return 0;
    }

    let word_length = word_length as i32;
    let lut_word_length = lut_word_length as i32;
    let mut q_end = *q_off + word_length;
    let mut s_end = *s_off + word_length;
    let context = crate::queryinfo::bsearch_context_info(q_end, query_info);
    let Some(context_info) = query_info.contexts.get(context.max(0) as usize) else {
        return 0;
    };
    let q_range = context_info.query_offset + context_info.query_length;

    if locations.is_some() {
        if s_is_seed_masked(
            lookup_wrap,
            subject_packed,
            s_end - lut_word_length,
            lut_word_length,
            q_end - lut_word_length,
        )
        .unwrap_or(true)
        {
            return 0;
        }
        loop {
            if !s_is_seed_masked(lookup_wrap, subject_packed, *s_off, lut_word_length, *q_off)
                .unwrap_or(true)
            {
                break;
            }
            *s_off += 1;
            *q_off += 1;
        }
    }

    let mut ext_to = word_length - (q_end - *q_off);
    let mut ext_max = (q_range - q_end).min(s_range as i32 - s_end);
    if ext_to > 0 || locations.is_some() {
        if ext_to > ext_max {
            return 0;
        }
        q_end += ext_to;
        s_end += ext_to;

        let mut s_pos = s_end - lut_word_length;
        let mut q_pos = q_end - lut_word_length;
        while s_pos > *s_off {
            if s_is_seed_masked(lookup_wrap, subject_packed, s_pos, lut_word_length, q_pos)
                .unwrap_or(true)
            {
                return 0;
            }
            s_pos -= lut_word_length;
            q_pos -= lut_word_length;
        }
        *extended = ext_to;
    }

    if !check_double {
        return 1;
    }

    ext_to += word_length;
    ext_max = ext_max.min(ext_to);
    let mut s_pos = s_end;
    let mut q_pos = q_end;
    while *extended + lut_word_length <= ext_max {
        if s_is_seed_masked(lookup_wrap, subject_packed, s_pos, lut_word_length, q_pos)
            .unwrap_or(true)
        {
            break;
        }
        s_pos += lut_word_length;
        q_pos += lut_word_length;
        *extended += lut_word_length;
    }

    s_pos -= lut_word_length - 1;
    q_pos -= lut_word_length - 1;
    while *extended < ext_max {
        if s_is_seed_masked(lookup_wrap, subject_packed, s_pos, lut_word_length, q_pos)
            .unwrap_or(true)
        {
            return 1;
        }
        *extended += 1;
        s_pos += 1;
        q_pos += 1;
    }

    if ext_max == ext_to {
        2
    } else {
        1
    }
}

/// Port of NCBI `BlastChooseNaLookupTable` (`blast_nalookup.c:46`).
pub fn blast_choose_na_lookup_table(
    lookup_options: Option<&crate::options::LookupTableOptions>,
    approx_table_entries: i32,
    max_q_off: i32,
    lut_width: &mut i32,
) -> Option<LookupTableType> {
    let lookup_options = lookup_options?;
    if lookup_options.mb_template_length > 0 {
        *lut_width = lookup_options.word_size;
        return Some(LookupTableType::MegablastLookup);
    }

    if crate::program::blast_program_is_mapping(lookup_options.program_number)
        && lookup_options.word_size >= 16
        && lookup_options.db_filter
    {
        *lut_width = 16;
        return Some(LookupTableType::NaHashLookup);
    }

    let mut lut_type;
    match lookup_options.word_size {
        4..=6 => {
            lut_type = LookupTableType::SmallNaLookup;
            *lut_width = lookup_options.word_size;
        }
        7 => {
            lut_type = LookupTableType::SmallNaLookup;
            *lut_width = if approx_table_entries < 250 { 6 } else { 7 };
        }
        8 => {
            lut_type = LookupTableType::SmallNaLookup;
            *lut_width = if approx_table_entries < 8500 { 7 } else { 8 };
        }
        9 => {
            if approx_table_entries < 1250 {
                *lut_width = 7;
                lut_type = LookupTableType::SmallNaLookup;
            } else if approx_table_entries < 21000 {
                *lut_width = 8;
                lut_type = LookupTableType::SmallNaLookup;
            } else {
                *lut_width = 9;
                lut_type = LookupTableType::MegablastLookup;
            }
        }
        10 => {
            if approx_table_entries < 1250 {
                *lut_width = 7;
                lut_type = LookupTableType::SmallNaLookup;
            } else if approx_table_entries < 8500 {
                *lut_width = 8;
                lut_type = LookupTableType::SmallNaLookup;
            } else if approx_table_entries < 18000 {
                *lut_width = 9;
                lut_type = LookupTableType::MegablastLookup;
            } else {
                *lut_width = 10;
                lut_type = LookupTableType::MegablastLookup;
            }
        }
        11 => {
            if approx_table_entries < 12000 {
                *lut_width = 8;
                lut_type = LookupTableType::SmallNaLookup;
            } else if approx_table_entries < 180000 {
                *lut_width = 10;
                lut_type = LookupTableType::MegablastLookup;
            } else {
                *lut_width = 11;
                lut_type = LookupTableType::MegablastLookup;
            }
        }
        12 => {
            if approx_table_entries < 8500 {
                *lut_width = 8;
                lut_type = LookupTableType::SmallNaLookup;
            } else if approx_table_entries < 18000 {
                *lut_width = 9;
                lut_type = LookupTableType::MegablastLookup;
            } else if approx_table_entries < 60000 {
                *lut_width = 10;
                lut_type = LookupTableType::MegablastLookup;
            } else if approx_table_entries < 900000 {
                *lut_width = 11;
                lut_type = LookupTableType::MegablastLookup;
            } else {
                *lut_width = 12;
                lut_type = LookupTableType::MegablastLookup;
            }
        }
        _ => {
            if approx_table_entries < 8500 {
                *lut_width = 8;
                lut_type = LookupTableType::SmallNaLookup;
            } else if approx_table_entries < 300000 {
                *lut_width = 11;
                lut_type = LookupTableType::MegablastLookup;
            } else {
                *lut_width = 12;
                lut_type = LookupTableType::MegablastLookup;
            }
        }
    }

    if lut_type == LookupTableType::SmallNaLookup
        && (approx_table_entries >= 32767 || max_q_off >= 32768)
    {
        lut_type = LookupTableType::NaLookup;
    }
    Some(lut_type)
}

/// Port of NCBI `EstimateNumTableEntries` (`lookup_util.c:190`).
pub fn estimate_num_table_entries(locations: &[crate::util::SSeqRange], max_off: &mut i32) -> i32 {
    let mut num_entries = 0i32;
    let mut curr_max = 0i32;
    for loc in locations {
        num_entries = num_entries.saturating_add(loc.right.saturating_sub(loc.left));
        curr_max = curr_max.max(loc.right);
    }
    *max_off = curr_max;
    num_entries
}

/// Rust ownership equivalent of NCBI `MapperWordHitsNew`.
pub fn mapper_word_hits_new(
    query: &crate::util::BlastSequenceBlk,
    query_info: &crate::queryinfo::QueryInfo,
) -> Option<MapperWordHits> {
    let num_arrays = (query_info.num_queries / 100).max(1);
    let array_size = 1000usize;
    let num_arrays_usize = num_arrays as usize;
    let context_slots = query_info
        .contexts
        .len()
        .max(query_info.num_queries.max(0) as usize);
    let mut last_pos = vec![0; context_slots];
    for slot in last_pos
        .iter_mut()
        .take(query_info.num_queries.max(0) as usize)
    {
        *slot = i32::MIN;
    }

    Some(MapperWordHits {
        pair_arrays: (0..num_arrays_usize)
            .map(|_| {
                vec![
                    OffsetPair {
                        query_offset: 0,
                        subject_offset: 0,
                    };
                    array_size
                ]
            })
            .collect(),
        num: vec![0; num_arrays_usize],
        num_arrays,
        array_size: array_size as i32,
        divisor: query.length / num_arrays + 1,
        last_diag: vec![0; context_slots],
        last_pos,
    })
}

/// Rust ownership equivalent of NCBI `MapperWordHitsFree`.
pub fn mapper_word_hits_free(word_hits: &mut Option<MapperWordHits>) -> Option<MapperWordHits> {
    *word_hits = None;
    None
}

/// 1-1 translation of `s_HasMaskAtHashEnabled` (`blast_nalookup.c:396`) over
/// the two C inputs that can enable mask-at-hash.
pub fn s_has_mask_at_hash_enabled(
    filtering_options_mask_at_hash: bool,
    filter_string: Option<&str>,
) -> bool {
    filtering_options_mask_at_hash || filter_string.is_some_and(|value| value.contains('m'))
}

/// 1-1 translation of `s_GetDiscTemplateType` (`blast_nalookup.c:638`).
pub fn s_get_disc_template_type(
    weight: i32,
    length: u8,
    word_type: DiscWordType,
) -> DiscTemplateType {
    match (weight, length, word_type) {
        (11, 16, DiscWordType::Coding | DiscWordType::TwoTemplates) => {
            DiscTemplateType::Template11_16Coding
        }
        (11, 16, DiscWordType::Optimal) => DiscTemplateType::Template11_16Optimal,
        (11, 18, DiscWordType::Coding | DiscWordType::TwoTemplates) => {
            DiscTemplateType::Template11_18Coding
        }
        (11, 18, DiscWordType::Optimal) => DiscTemplateType::Template11_18Optimal,
        (11, 21, DiscWordType::Coding | DiscWordType::TwoTemplates) => {
            DiscTemplateType::Template11_21Coding
        }
        (11, 21, DiscWordType::Optimal) => DiscTemplateType::Template11_21Optimal,
        (12, 16, DiscWordType::Coding | DiscWordType::TwoTemplates) => {
            DiscTemplateType::Template12_16Coding
        }
        (12, 16, DiscWordType::Optimal) => DiscTemplateType::Template12_16Optimal,
        (12, 18, DiscWordType::Coding | DiscWordType::TwoTemplates) => {
            DiscTemplateType::Template12_18Coding
        }
        (12, 18, DiscWordType::Optimal) => DiscTemplateType::Template12_18Optimal,
        (12, 21, DiscWordType::Coding | DiscWordType::TwoTemplates) => {
            DiscTemplateType::Template12_21Coding
        }
        (12, 21, DiscWordType::Optimal) => DiscTemplateType::Template12_21Optimal,
        _ => DiscTemplateType::Contiguous,
    }
}

/// Port of NCBI inline `ComputeDiscontiguousIndex` (`blast_nalookup.h:572`).
pub fn compute_discontiguous_index(accum: u64, template_type: DiscTemplateType) -> i32 {
    let lo = accum as u32;
    let hi = (accum >> 32) as u32;
    let index = match template_type {
        DiscTemplateType::Template11_16Coding => {
            (lo & 0x00000003)
                | ((lo & 0x000000f0) >> 2)
                | ((lo & 0x00003c00) >> 4)
                | ((lo & 0x000f0000) >> 6)
                | ((lo & 0x03c00000) >> 8)
                | ((lo & 0xf0000000) >> 10)
        }
        DiscTemplateType::Template11_16Optimal => {
            (lo & 0x0000003f)
                | ((lo & 0x00000f00) >> 2)
                | ((lo & 0x0003c000) >> 4)
                | ((lo & 0x00300000) >> 6)
                | ((lo & 0xfc000000) >> 10)
        }
        DiscTemplateType::Template11_18Coding => {
            (lo & 0x00000003)
                | ((lo & 0x000000f0) >> 2)
                | ((lo & 0x00003c00) >> 4)
                | ((lo & 0x00030000) >> 6)
                | ((lo & 0x03c00000) >> 10)
                | ((lo & 0xf0000000) >> 12)
                | ((hi & 0x0000000c) << 18)
        }
        DiscTemplateType::Template11_18Optimal => {
            (lo & 0x0000003f)
                | ((lo & 0x00000300) >> 2)
                | ((lo & 0x0003c000) >> 6)
                | ((lo & 0x00300000) >> 8)
                | ((lo & 0x0c000000) >> 12)
                | ((lo & 0xc0000000) >> 14)
                | ((hi & 0x0000000f) << 18)
        }
        DiscTemplateType::Template11_21Coding => {
            (lo & 0x00000003)
                | ((lo & 0x000000f0) >> 2)
                | ((lo & 0x00000c00) >> 4)
                | ((lo & 0x000f0000) >> 8)
                | ((lo & 0x00c00000) >> 10)
                | ((lo & 0xf0000000) >> 14)
                | ((hi & 0x0000000c) << 16)
                | ((hi & 0x00000300) << 12)
        }
        DiscTemplateType::Template11_21Optimal => {
            (lo & 0x0000003f)
                | ((lo & 0x00000300) >> 2)
                | ((lo & 0x0000c000) >> 6)
                | ((lo & 0x00c00000) >> 12)
                | ((lo & 0x0c000000) >> 14)
                | ((hi & 0x00000003) << 14)
                | ((hi & 0x000003f0) << 12)
        }
        DiscTemplateType::Template12_16Coding => {
            (lo & 0x00000003)
                | ((lo & 0x000000f0) >> 2)
                | ((lo & 0x00003c00) >> 4)
                | ((lo & 0x000f0000) >> 6)
                | ((lo & 0xffc00000) >> 8)
        }
        DiscTemplateType::Template12_16Optimal => {
            (lo & 0x0000003f)
                | ((lo & 0x00000f00) >> 2)
                | ((lo & 0x0003c000) >> 4)
                | ((lo & 0x00f00000) >> 6)
                | ((lo & 0xfc000000) >> 8)
        }
        DiscTemplateType::Template12_18Coding => {
            (lo & 0x00000003)
                | ((lo & 0x000000f0) >> 2)
                | ((lo & 0x00003c00) >> 4)
                | ((lo & 0x000f0000) >> 6)
                | ((lo & 0x03c00000) >> 8)
                | ((lo & 0xf0000000) >> 10)
                | ((hi & 0x0000000c) << 20)
        }
        DiscTemplateType::Template12_18Optimal => {
            (lo & 0x0000003f)
                | ((lo & 0x00000f00) >> 2)
                | ((lo & 0x0000c000) >> 4)
                | ((lo & 0x00f00000) >> 8)
                | ((lo & 0x0c000000) >> 10)
                | ((lo & 0xc0000000) >> 12)
                | ((hi & 0x0000000f) << 20)
        }
        DiscTemplateType::Template12_21Coding => {
            (lo & 0x00000003)
                | ((lo & 0x000000f0) >> 2)
                | ((lo & 0x00000c00) >> 4)
                | ((lo & 0x000f0000) >> 8)
                | ((lo & 0x03c00000) >> 10)
                | ((lo & 0xf0000000) >> 12)
                | ((hi & 0x0000000c) << 18)
                | ((hi & 0x00000300) << 14)
        }
        DiscTemplateType::Template12_21Optimal => {
            (lo & 0x0000003f)
                | ((lo & 0x00000300) >> 2)
                | ((lo & 0x0000c000) >> 6)
                | ((lo & 0x00f00000) >> 10)
                | ((lo & 0x0c000000) >> 12)
                | ((hi & 0x00000003) << 16)
                | ((hi & 0x000003f0) << 14)
        }
        DiscTemplateType::Contiguous => 0,
    };
    index as i32
}

/// 1-1 ownership translation of `BlastSmallNaLookupTableDestruct`.
pub fn blast_small_na_lookup_table_destruct(
    lookup: &mut Option<SmallNaLookupTable>,
) -> Option<SmallNaLookupTable> {
    *lookup = None;
    None
}

fn nucleotide_backbone_size(lut_width: i32) -> Option<usize> {
    if lut_width <= 0 {
        return None;
    }
    let bits = BITS_PER_NUC.checked_mul(lut_width as u32)?;
    let size = 1usize.checked_shl(bits)?;
    i32::try_from(size).ok()?;
    Some(size)
}

/// Rust ownership translation of NCBI `BlastSmallNaLookupTableNew`
/// (`blast_nalookup.c:378`).
pub fn blast_small_na_lookup_table_new(
    query_sequence: &[u8],
    locations: &[crate::util::SSeqRange],
    lut: &mut Option<SmallNaLookupTable>,
    lookup_options: &crate::options::LookupTableOptions,
    lut_width: i32,
) -> i16 {
    let Some(backbone_size) = nucleotide_backbone_size(lut_width) else {
        *lut = None;
        return -1;
    };
    if lookup_options.word_size < lut_width {
        *lut = None;
        return -1;
    }

    let mut lookup = SmallNaLookupTable {
        word_length: lookup_options.word_size,
        backbone: Vec::new(),
        overflow: Vec::new(),
        pv_array: Vec::new(),
        longest_chain: 0,
        scan_step: lookup_options.word_size - lut_width + 1,
    };
    if lookup.scan_step <= 0 {
        *lut = None;
        return -1;
    }

    let mut thin_backbone = vec![None; backbone_size];
    blast_lookup_index_query_exact_matches(
        &mut thin_backbone,
        lookup.word_length,
        BITS_PER_NUC as i32,
        lut_width,
        query_sequence,
        locations,
    );

    let status = s_blast_small_na_lookup_finalize(&mut thin_backbone, &mut lookup);
    if status != 0 {
        *lut = None;
        return status as i16;
    }

    *lut = Some(lookup);
    0
}

/// Port of NCBI `s_BlastSmallNaLookupFinalize` (`blast_nalookup.c:197`).
/// `thin_backbone` rows follow the C layout: slot 1 stores hit count and slots
/// 2.. hold query offsets.
pub fn s_blast_small_na_lookup_finalize(
    thin_backbone: &mut [Option<Vec<i32>>],
    lookup: &mut SmallNaLookupTable,
) -> i32 {
    let mut overflow_cells_needed = 2usize;
    let mut longest_chain = 0;

    for row in thin_backbone.iter().flatten() {
        let num_hits = row.get(1).copied().unwrap_or(0).max(0) as usize;
        if num_hits > 1 {
            overflow_cells_needed += num_hits + 1;
        }
        longest_chain = longest_chain.max(num_hits as i32);
    }

    if overflow_cells_needed >= 32768 {
        for row in thin_backbone {
            *row = None;
        }
        return -1;
    }

    lookup.backbone = vec![-1; thin_backbone.len()];
    lookup.overflow = vec![0; overflow_cells_needed];
    lookup.longest_chain = longest_chain;

    let mut overflow_cursor = 2usize;
    for (index, row) in thin_backbone.iter_mut().enumerate() {
        let Some(values) = row.take() else {
            continue;
        };
        let num_hits = values.get(1).copied().unwrap_or(0).max(0) as usize;
        if num_hits == 0 {
            continue;
        }
        if num_hits == 1 {
            lookup.backbone[index] = values.get(2).copied().unwrap_or(0);
        } else {
            lookup.backbone[index] = -(overflow_cursor as i32);
            for hit_index in 0..num_hits {
                lookup.overflow[overflow_cursor] = values.get(hit_index + 2).copied().unwrap_or(0);
                overflow_cursor += 1;
            }
            lookup.overflow[overflow_cursor] = -1;
            overflow_cursor += 1;
        }
    }
    lookup.overflow.truncate(overflow_cursor);
    0
}

/// NCBI: s_ComputeTableIndex (blast_lookup.c).
fn s_compute_table_index(wordsize: i32, charsize: i32, seq: &[u8]) -> usize {
    let mut index = 0usize;
    for &base in seq.iter().take(wordsize.max(0) as usize) {
        index = (index << charsize.max(0) as usize) | base as usize;
    }
    index
}

/// NCBI: BlastLookupAddWordHit (blast_lookup.c:33).
fn blast_lookup_add_word_hit(
    backbone: &mut [Option<Vec<i32>>],
    wordsize: i32,
    charsize: i32,
    seq: &[u8],
    query_offset: i32,
) {
    let index = s_compute_table_index(wordsize, charsize, seq);
    let Some(cell) = backbone.get_mut(index) else {
        return;
    };
    let chain = cell.get_or_insert_with(|| vec![8, 0]);
    chain.push(query_offset);
    chain[1] += 1;
    if chain[1] + 2 > chain[0] {
        chain[0] *= 2;
    }
}

/// Port of NCBI `BlastLookupIndexQueryExactMatches` (`blast_lookup.c:79`).
pub fn blast_lookup_index_query_exact_matches(
    backbone: &mut [Option<Vec<i32>>],
    word_length: i32,
    charsize: i32,
    lut_word_length: i32,
    query_sequence: &[u8],
    locations: &[crate::util::SSeqRange],
) {
    let word_length = word_length.max(0) as usize;
    let lut_word_length_usize = lut_word_length.max(0) as usize;
    let invalid_mask = ((0xffu16 << charsize.max(0) as u32) & 0xff) as u8;

    for loc in locations {
        if loc.left < 0 || loc.right < loc.left {
            continue;
        }
        let from = loc.left as usize;
        let to = loc.right as usize;
        if word_length > to - from + 1 || lut_word_length_usize == 0 {
            continue;
        }

        let last_start = to.saturating_sub(lut_word_length_usize - 1);
        for start in from..=last_start {
            let Some(word) = query_sequence.get(start..start + lut_word_length_usize) else {
                break;
            };
            if word.iter().any(|base| (base & invalid_mask) != 0) {
                continue;
            }
            blast_lookup_add_word_hit(backbone, lut_word_length, charsize, word, start as i32);
        }
    }
}

/// Port of NCBI `BlastHashLookupIndexQueryExactMatches` (`blast_lookup.c:237`).
pub fn blast_hash_lookup_index_query_exact_matches(
    backbone: &mut [BackboneCell],
    offsets: &mut [i32],
    word_length: i32,
    charsize: i32,
    lut_word_length: i32,
    query_sequence: &[u8],
    locations: &[crate::util::SSeqRange],
    hash_func: fn(&[u8], u32) -> u32,
    mask: u32,
    pv_array: Option<&[u32]>,
) -> i16 {
    let word_length_usize = word_length.max(0) as usize;
    let lut_word_length_usize = lut_word_length.max(0) as usize;
    if word_length_usize == 0 || lut_word_length_usize == 0 {
        return 0;
    }
    let invalid_mask = ((0xffu16 << charsize.max(0) as u32) & 0xff) as u8;

    for loc in locations {
        if loc.left < 0 || loc.right < loc.left {
            continue;
        }
        let from = loc.left as usize;
        let to = loc.right as usize;
        if word_length_usize > to - from + 1 || lut_word_length_usize > word_length_usize {
            continue;
        }

        let last_start = to + 1 - lut_word_length_usize;
        for start in from..=last_start {
            let Some(word) = query_sequence.get(start..start + lut_word_length_usize) else {
                break;
            };
            if word.iter().any(|base| (base & invalid_mask) != 0) {
                continue;
            }
            if s_add_word_hit(
                backbone,
                offsets,
                lut_word_length,
                charsize,
                word,
                start as i32,
                hash_func,
                mask,
                pv_array,
            ) != 0
            {
                return -1;
            }
        }
    }

    0
}

/// Port of NCBI `s_BlastNaLookupFinalize` (`blast_nalookup.c:443`).
pub fn s_blast_na_lookup_finalize(
    thin_backbone: &mut [Option<Vec<i32>>],
    lookup: &mut BlastNaLookupTable,
) {
    let backbone_size = thin_backbone.len();
    let mut overflow_cells_needed = 0usize;
    let mut longest_chain = 0;

    for row in thin_backbone.iter().flatten() {
        let num_hits = row.get(1).copied().unwrap_or(0).max(0) as usize;
        if num_hits > NA_HITS_PER_CELL {
            overflow_cells_needed += num_hits;
        }
        longest_chain = longest_chain.max(num_hits as i32);
    }

    lookup.backbone_size = backbone_size as i32;
    lookup.thick_backbone = vec![NaLookupBackboneCell::default(); backbone_size];
    lookup.pv = vec![0; (backbone_size >> crate::stat::PV_ARRAY_BTS) + 1];
    lookup.overflow = vec![0; overflow_cells_needed];
    lookup.longest_chain = longest_chain;

    let mut overflow_cursor = 0usize;
    for (index, row) in thin_backbone.iter_mut().enumerate() {
        let Some(values) = row.take() else {
            continue;
        };
        let num_hits = values.get(1).copied().unwrap_or(0).max(0) as usize;
        if num_hits == 0 {
            continue;
        }
        lookup.thick_backbone[index].num_used = num_hits as i32;
        if let Some(slot) = lookup.pv.get_mut(index >> crate::stat::PV_ARRAY_BTS) {
            *slot |= 1u32 << ((index as u32) & crate::stat::PV_ARRAY_MASK);
        }

        if num_hits <= NA_HITS_PER_CELL {
            for hit_index in 0..num_hits {
                lookup.thick_backbone[index].entries[hit_index] =
                    values.get(hit_index + 2).copied().unwrap_or(0);
            }
        } else {
            lookup.thick_backbone[index].overflow_cursor = overflow_cursor as i32;
            for hit_index in 0..num_hits {
                lookup.overflow[overflow_cursor] = values.get(hit_index + 2).copied().unwrap_or(0);
                overflow_cursor += 1;
            }
        }
    }
    lookup.overflow_size = overflow_cursor as i32;
    lookup.overflow.truncate(overflow_cursor);
}

/// 1-1 ownership translation of `BlastNaLookupTableDestruct`.
pub fn blast_na_lookup_table_destruct(
    lookup: &mut Option<BlastNaLookupTable>,
) -> Option<BlastNaLookupTable> {
    *lookup = None;
    None
}

/// Rust ownership translation of NCBI `BlastNaLookupTableNew`
/// (`blast_nalookup.c:548`).
pub fn blast_na_lookup_table_new(
    query_sequence: &[u8],
    locations: &[crate::util::SSeqRange],
    lut: &mut Option<BlastNaLookupTable>,
    lookup_options: &crate::options::LookupTableOptions,
    lut_width: i32,
) -> i16 {
    let Some(backbone_size) = nucleotide_backbone_size(lut_width) else {
        *lut = None;
        return -1;
    };
    if lookup_options.word_size < lut_width {
        *lut = None;
        return -1;
    }

    let mut lookup = BlastNaLookupTable {
        mask: backbone_size as i32 - 1,
        word_length: lookup_options.word_size,
        lut_word_length: lut_width,
        scan_step: lookup_options.word_size - lut_width + 1,
        backbone_size: backbone_size as i32,
        ..Default::default()
    };
    if lookup.scan_step <= 0 {
        *lut = None;
        return -1;
    }

    let mut thin_backbone = vec![None; backbone_size];
    blast_lookup_index_query_exact_matches(
        &mut thin_backbone,
        lookup.word_length,
        BITS_PER_NUC as i32,
        lookup.lut_word_length,
        query_sequence,
        locations,
    );

    s_blast_na_lookup_finalize(&mut thin_backbone, &mut lookup);
    *lut = Some(lookup);
    0
}

/// Port of NCBI `s_BlastNaHashLookupFinalize` (`blast_nalookup.c:2019`).
pub fn s_blast_na_hash_lookup_finalize(
    thin_backbone: &mut [BackboneCell],
    offsets: &[i32],
    lookup: &mut BlastNaHashLookupTable,
) {
    let backbone_size = thin_backbone.len();
    let pv_array_bts = if lookup.pv_array_bts > 0 {
        lookup.pv_array_bts as usize
    } else {
        crate::stat::PV_ARRAY_BTS
    };
    let mut overflow_cells_needed = 0usize;
    let mut longest_chain = 0;
    let mut max_word = 0usize;

    for head in thin_backbone.iter() {
        let mut num_hits = 0usize;
        let mut num_words = 0usize;
        if head.num_offsets > 0 {
            let mut cursor = Some(head);
            while let Some(cell) = cursor {
                num_hits += cell.num_offsets.max(0) as usize;
                num_words += 1;
                max_word = max_word.max(cell.word as usize);
                cursor = cell.next.as_deref();
            }
        }
        if num_words > NA_WORDS_PER_HASH || num_hits > NA_OFFSETS_PER_HASH {
            overflow_cells_needed += num_hits + (num_words * 2);
        }
        longest_chain = longest_chain.max(num_hits as i32);
    }

    lookup.backbone_size = backbone_size as i32;
    lookup.longest_chain = longest_chain;
    lookup.pv_array_bts = pv_array_bts as i32;
    lookup.thick_backbone = vec![NaHashLookupBackboneCell::default(); backbone_size];
    lookup.overflow = vec![0; overflow_cells_needed];
    lookup.pv.clear();
    lookup.pv.resize((max_word >> pv_array_bts) + 1, 0);

    let mut overflow_cursor = 0usize;
    for (index, head) in thin_backbone.iter_mut().enumerate() {
        if head.num_offsets == 0 {
            continue;
        }

        let mut num_words = 0usize;
        let mut num_offsets = 0usize;
        let mut cursor = Some(&*head);
        while let Some(cell) = cursor {
            num_words += 1;
            num_offsets += cell.num_offsets.max(0) as usize;
            cursor = cell.next.as_deref();
        }

        lookup.thick_backbone[index].num_words = num_words as i8;
        let is_overflow = num_words > NA_WORDS_PER_HASH || num_offsets > NA_OFFSETS_PER_HASH;

        if !is_overflow {
            let mut word_index = 0usize;
            let mut offset_index = 0usize;
            let mut cursor = Some(&*head);
            while let Some(cell) = cursor {
                lookup.thick_backbone[index].words[word_index] = cell.word;
                lookup.thick_backbone[index].num_offsets[word_index] = cell.num_offsets as i8;
                set_pv_bit(&mut lookup.pv, cell.word as usize, pv_array_bts);

                let mut j = cell.offset;
                while j != 0 && offset_index < NA_OFFSETS_PER_HASH {
                    lookup.thick_backbone[index].offsets[offset_index] = j - 1;
                    offset_index += 1;
                    j = offsets.get(j as usize).copied().unwrap_or(0);
                }

                word_index += 1;
                cursor = cell.next.as_deref();
            }
        } else {
            if num_words <= NA_WORDS_PER_HASH {
                let mut word_index = 0usize;
                let mut cursor = Some(&*head);
                while let Some(cell) = cursor {
                    lookup.thick_backbone[index].words[word_index] = cell.word;
                    word_index += 1;
                    cursor = cell.next.as_deref();
                }
            }

            lookup.thick_backbone[index].offsets[0] = overflow_cursor as i32;
            let mut cursor = Some(&*head);
            while let Some(cell) = cursor {
                lookup.overflow[overflow_cursor] = cell.word as i32;
                overflow_cursor += 1;
                lookup.overflow[overflow_cursor] = cell.num_offsets;
                overflow_cursor += 1;

                let mut j = cell.offset;
                while j != 0 {
                    lookup.overflow[overflow_cursor] = j - 1;
                    overflow_cursor += 1;
                    j = offsets.get(j as usize).copied().unwrap_or(0);
                }
                set_pv_bit(&mut lookup.pv, cell.word as usize, pv_array_bts);
                cursor = cell.next.as_deref();
            }
        }

        let mut next = head.next.take();
        let _ = backbone_cell_free(&mut next);
    }
    lookup.offsets_size = overflow_cursor as i32;
    lookup.overflow.truncate(overflow_cursor);
}

fn set_pv_bit(pv: &mut Vec<u32>, index: usize, pv_array_bts: usize) {
    let word = index >> pv_array_bts;
    if word >= pv.len() {
        pv.resize(word + 1, 0);
    }
    pv[word] |= 1u32 << ((index as u32) & crate::stat::PV_ARRAY_MASK);
}

fn clear_pv_bit(pv: &mut [u32], index: u32, pv_array_bts: usize) {
    let word = (index >> pv_array_bts) as usize;
    if let Some(slot) = pv.get_mut(word) {
        *slot &= !(1u32 << (index & crate::stat::PV_ARRAY_MASK));
    }
}

fn test_pv_bit(pv: &[u32], index: u64, pv_array_bts: usize) -> bool {
    let word = (index >> pv_array_bts) as usize;
    pv.get(word)
        .is_some_and(|slot| (*slot & (1u32 << ((index as u32) & crate::stat::PV_ARRAY_MASK))) != 0)
}

/// Port of NCBI static `FNV_hash` (`blast_nalookup.c:1390`).
pub fn fnv_hash(seq: &[u8], mask: u32) -> u32 {
    const FNV_PRIME: u32 = 16_777_619;
    const FNV_OFFSET_BASIS: u32 = 2_166_136_261;

    let mut hash = FNV_OFFSET_BASIS;
    for &byte in seq.iter().take(4) {
        hash = hash.wrapping_mul(FNV_PRIME);
        hash ^= u32::from(byte);
    }
    for _ in seq.len().min(4)..4 {
        hash = hash.wrapping_mul(FNV_PRIME);
    }
    hash & mask
}

/// Port of NCBI `s_NaHashLookupFillPV` (`blast_nalookup.c:1408`).
pub fn s_na_hash_lookup_fill_pv(
    query_sequence: &[u8],
    locations: &[crate::util::SSeqRange],
    lookup: &mut BlastNaHashLookupTable,
) -> i16 {
    let word_length = lookup.word_length.max(0) as usize;
    let lut_word_length = lookup.lut_word_length.max(0) as usize;
    let pv_array_bts = lookup.pv_array_bts.max(0) as usize;

    for loc in locations {
        if loc.left < 0 || loc.right < loc.left {
            continue;
        }
        let from = loc.left.max(0) as usize;
        let to = loc.right as usize;
        if word_length > to.saturating_sub(from) + 1 {
            continue;
        }
        let Some(window) = query_sequence.get(from..=to) else {
            continue;
        };

        let mut ecode = 0u32;
        let mut pos = lut_word_length.saturating_sub(1);
        for (seq_index, &base) in window.iter().enumerate() {
            if (base & BLAST2NA_MASK) != 0 {
                ecode = 0;
                pos = seq_index + lut_word_length;
                continue;
            }

            ecode = ecode.wrapping_shl(BITS_PER_NUC) | u32::from(base);
            if seq_index < pos {
                continue;
            }
            set_pv_bit(&mut lookup.pv, ecode as usize, pv_array_bts);
        }
    }

    0
}

/// Port of NCBI `s_NaHashLookupCountWordsInSubject_16_1`
/// (`blast_nalookup.c:1613`) over packed 2-bit subject bytes.
pub fn s_na_hash_lookup_count_words_in_subject_16_1(
    sequence: Option<&[u8]>,
    sequence_length: i32,
    lookup: Option<&BlastNaHashLookupTable>,
    counts: Option<&mut crate::util::BlastSparseUint1Array>,
    max_word_count: u8,
) -> i16 {
    let (Some(sequence), Some(lookup), Some(counts)) = (sequence, lookup, counts) else {
        return -1;
    };
    if lookup.pv.is_empty() || sequence_length < lookup.lut_word_length {
        return -1;
    }

    let k_num_words = sequence_length - lookup.lut_word_length;
    let max_read_index = if k_num_words > 0 {
        ((k_num_words - 1) / 4 + 4) as usize
    } else {
        3
    };
    if sequence.len() <= max_read_index {
        return -1;
    }

    let mask = (1u64 << (16 * BITS_PER_NUC)) - 1;
    let pv_array_bts = lookup.pv_array_bts.max(0) as usize;
    let mut shift = 8u32;
    let mut w = (u64::from(sequence[0]) << 24)
        | (u64::from(sequence[1]) << 16)
        | (u64::from(sequence[2]) << 8)
        | u64::from(sequence[3]);

    for i in 0..k_num_words {
        if i % 4 == 0 {
            shift = 8;
            w = (w << 8) | u64::from(sequence[(i / 4 + 4) as usize]);
        } else {
            shift -= 2;
        }

        let word = (w >> shift) & mask;
        if !test_pv_bit(&lookup.pv, word, pv_array_bts) {
            continue;
        }

        if let Some(value) =
            crate::util::blast_sparse_uint1_array_get_element(Some(&mut *counts), word as i64)
        {
            if *value < max_word_count {
                *value += 1;
            }
        }
    }

    0
}

/// Port-shaped equivalent of NCBI static
/// `s_NaHashLookupScanSubjectForWordCounts` (`blast_nalookup.c:1839`).
///
/// C allocates a sparse counter over the full 2^32 word space. Rust keeps the
/// same sparse bitfield/value logic but limits the backing bitfield to the PV
/// words actually present in the grow-on-set Rust PV representation.
pub fn s_na_hash_lookup_scan_subject_for_word_counts(
    seq_src: Option<&mut dyn crate::seqsrc::BlastSeqSource>,
    lookup: Option<&mut BlastNaHashLookupTable>,
    in_num_threads: u32,
    max_word_count: u8,
) -> i16 {
    let (Some(seq_src), Some(lookup)) = (seq_src, lookup) else {
        return -1;
    };
    if lookup.pv.is_empty() || lookup.lut_word_length != 16 {
        return -1;
    }

    let num_db_seqs = seq_src.num_seqs();
    if num_db_seqs <= 0 || in_num_threads == 0 {
        return -1;
    }
    let _num_threads = in_num_threads.min(num_db_seqs as u32);
    let sparse_len = (lookup.pv.len() as i64) << crate::stat::PV_ARRAY_BTS;
    let Some(mut word_counts) =
        crate::util::blast_sparse_uint1_array_new(Some(&lookup.pv), sparse_len)
    else {
        return -1;
    };

    crate::seqsrc::blast_seq_src_reset_chunk_iterator(seq_src);
    let oids: Vec<i32> = seq_src.iter_oids().collect();
    for oid in oids {
        let Some(seq_data) = seq_src.get_sequence(&crate::seqsrc::GetSeqArg {
            oid,
            encoding: crate::seqsrc::SeqEncoding::Protein,
            ..crate::seqsrc::GetSeqArg::default()
        }) else {
            continue;
        };
        let status = s_na_hash_lookup_count_words_in_subject_16_1(
            Some(&seq_data.sequence),
            seq_data.length,
            Some(&*lookup),
            Some(&mut word_counts),
            max_word_count,
        );
        if status != 0 {
            return status;
        }
    }

    let mut value_index = 0usize;
    for word_index in 0..word_counts.length as usize {
        let mut bit = 1u32;
        while bit != 0 {
            if (word_counts.bitfield[word_index] & bit) != 0 {
                let Some(&count) = word_counts.values.get(value_index) else {
                    return -1;
                };
                if count == 0 || count >= max_word_count {
                    word_counts.bitfield[word_index] &= !bit;
                }
                value_index += 1;
                if value_index >= word_counts.num_elements as usize {
                    break;
                }
            }
            bit = bit.wrapping_shl(1);
        }
    }

    lookup.pv = word_counts.bitfield;
    0
}

/// Port of NCBI `s_NaHashLookupRemovePolyAWords` (`blast_nalookup.c:1965`).
pub fn s_na_hash_lookup_remove_poly_a_words(lookup: Option<&mut BlastNaHashLookupTable>) -> i16 {
    let Some(lookup) = lookup else {
        return -1;
    };
    if lookup.lut_word_length != 16 || lookup.pv_array_bts != crate::stat::PV_ARRAY_BTS as i32 {
        return -1;
    }
    let pv_array_bts = lookup.pv_array_bts.max(0) as usize;
    let word_size = lookup.lut_word_length.max(0);

    clear_pv_bit(&mut lookup.pv, 0, pv_array_bts);
    clear_pv_bit(&mut lookup.pv, u32::MAX, pv_array_bts);

    for i in 1..4u32 {
        for k in 0..word_size.min(16) {
            let word = i << (k * 2);
            clear_pv_bit(&mut lookup.pv, word, pv_array_bts);
        }
    }

    for i in 0..3u32 {
        for k in 0..word_size.min(16) {
            let shift = (k * 2) as u32;
            let word = ((u32::MAX ^ (3u32 << shift)) | (i << shift)) & u32::MAX;
            clear_pv_bit(&mut lookup.pv, word, pv_array_bts);
        }
    }

    0
}

/// Rust ownership equivalent of NCBI `NaHashLookupThreadDataFree`.
pub fn na_hash_lookup_thread_data_free(
    thread_data: &mut Option<NaHashLookupThreadData>,
) -> Option<NaHashLookupThreadData> {
    if let Some(thread_data) = thread_data.as_mut() {
        for seq_arg in thread_data.seq_arg.iter_mut() {
            seq_arg.oid = 0;
        }
        thread_data.seq_arg.clear();

        for itr in thread_data.itr.drain(..) {
            let mut slot = Some(itr);
            let _ = crate::seqsrc::blast_seq_src_iterator_free(&mut slot);
        }

        for seq_src in thread_data.seq_src.drain(..) {
            let _ = crate::seqsrc::blast_seq_src_free(Some(seq_src));
        }

        for word_counts in thread_data.word_counts.drain(..) {
            let _ = crate::util::blast_sparse_uint1_array_free(Some(word_counts));
        }
        thread_data.num_threads = 0;
    }
    *thread_data = None;
    None
}

/// Port of NCBI static `NaHashLookupThreadDataNew`
/// (`blast_nalookup.c:1737`).
pub fn na_hash_lookup_thread_data_new(
    num_threads: i32,
    lookup: Option<&BlastNaHashLookupTable>,
    seq_src: Option<&crate::seqsrc::BlastSeqSrc>,
) -> Option<NaHashLookupThreadData> {
    let (Some(lookup), Some(seq_src)) = (lookup, seq_src) else {
        return None;
    };
    if num_threads < 1 {
        return None;
    }

    let sparse_len = 1i64.checked_shl((2 * lookup.lut_word_length).max(0) as u32)?;
    let first_counts = crate::util::blast_sparse_uint1_array_new(Some(&lookup.pv), sparse_len)?;
    let mut seq_arg = Vec::with_capacity(num_threads as usize);
    let mut itr = Vec::with_capacity(num_threads as usize);
    let mut seq_src_vec = Vec::with_capacity(num_threads as usize);
    let mut word_counts = Vec::with_capacity(num_threads as usize);

    for thread_index in 0..num_threads {
        seq_arg.push(crate::seqsrc::GetSeqArg {
            oid: 0,
            encoding: crate::seqsrc::SeqEncoding::Protein,
            ..crate::seqsrc::GetSeqArg::default()
        });
        seq_src_vec.push(crate::seqsrc::blast_seq_src_copy(Some(seq_src))?);
        itr.push(crate::seqsrc::blast_seq_src_iterator_new_ex(1));

        if thread_index == 0 {
            word_counts.push(first_counts.clone());
        } else {
            let mut counts = first_counts.clone();
            counts.values.fill(0);
            word_counts.push(counts);
        }
    }

    Some(NaHashLookupThreadData {
        seq_arg,
        itr,
        seq_src: seq_src_vec,
        word_counts,
        num_threads,
    })
}

/// 1-1 ownership translation of `BlastNaHashLookupTableDestruct`.
pub fn blast_na_hash_lookup_table_destruct(
    lookup: &mut Option<BlastNaHashLookupTable>,
) -> Option<BlastNaHashLookupTable> {
    *lookup = None;
    None
}

/// Rust ownership translation of NCBI `BlastNaHashLookupTableNew`
/// (`blast_nalookup.c:2222`).
pub fn blast_na_hash_lookup_table_new(
    query_sequence: &[u8],
    locations: &[crate::util::SSeqRange],
    lut: &mut Option<BlastNaHashLookupTable>,
    lookup_options: &crate::options::LookupTableOptions,
    seq_src: Option<&mut dyn crate::seqsrc::BlastSeqSource>,
    num_threads: u32,
) -> i16 {
    if lookup_options.db_filter && seq_src.is_none() {
        return -1;
    }

    let mut lookup = BlastNaHashLookupTable {
        word_length: lookup_options.word_size,
        lut_word_length: 16,
        scan_step: if lookup_options.db_filter {
            1
        } else {
            lookup_options.word_size - 16 + 1
        },
        pv_array_bts: crate::stat::PV_ARRAY_BTS as i32,
        ..Default::default()
    };
    if lookup.word_length < lookup.lut_word_length || lookup.scan_step <= 0 {
        *lut = None;
        return -1;
    }

    if s_na_hash_lookup_fill_pv(query_sequence, locations, &mut lookup) != 0 {
        *lut = None;
        return -1;
    }
    if s_na_hash_lookup_remove_poly_a_words(Some(&mut lookup)) != 0 {
        *lut = None;
        return -1;
    }

    if lookup_options.db_filter {
        let status = s_na_hash_lookup_scan_subject_for_word_counts(
            seq_src,
            Some(&mut lookup),
            num_threads,
            lookup_options.max_db_word_count,
        );
        if status != 0 {
            *lut = None;
            return status;
        }
    }

    let num_unique_words: i32 = lookup
        .pv
        .iter()
        .map(|&word| crate::util::s_popcount(word) as i32)
        .sum();
    let mut num_hash_bits = 8u32;
    while num_hash_bits < 32 && (1i64 << num_hash_bits) < i64::from(num_unique_words) {
        num_hash_bits += 1;
    }

    lookup.backbone_size = 1i32.checked_shl(num_hash_bits).unwrap_or(i32::MAX);
    lookup.mask = lookup.backbone_size.saturating_sub(1);
    let backbone_size = lookup.backbone_size.max(0) as usize;
    let mut thin_backbone = vec![BackboneCell::default(); backbone_size];
    let mut offsets = vec![0; query_sequence.len() + 1];

    let status = blast_hash_lookup_index_query_exact_matches(
        &mut thin_backbone,
        &mut offsets,
        lookup.word_length,
        BITS_PER_NUC as i32,
        lookup.lut_word_length,
        query_sequence,
        locations,
        fnv_hash,
        lookup.mask as u32,
        Some(&lookup.pv),
    );
    if status != 0 {
        *lut = None;
        return status;
    }

    s_blast_na_hash_lookup_finalize(&mut thin_backbone, &offsets, &mut lookup);
    *lut = Some(lookup);
    0
}

/// Port of NCBI static `s_BlastNaHashLookupRetieveHits`
/// (`blast_nascan.c:2701`) over the finalized Rust-owned hash table.
pub fn s_blast_na_hash_lookup_retrieve_hits(
    lookup: &BlastNaHashLookupTable,
    index: u32,
    s_off: i32,
    offset_pairs: &mut Vec<OffsetPair>,
) -> i32 {
    let pv_array_bts = if lookup.pv_array_bts > 0 {
        lookup.pv_array_bts as usize
    } else {
        0
    };
    if !test_pv_bit(&lookup.pv, u64::from(index), pv_array_bts) {
        return 0;
    }

    let hashed_index = fnv_hash(&index.to_ne_bytes(), lookup.mask as u32) as usize;
    if hashed_index >= lookup.thick_backbone.len() {
        return 0;
    }
    let cell = &lookup.thick_backbone[hashed_index];
    if cell.num_words <= 0 {
        return 0;
    }

    let mut num_hits = 0;
    let mut cursor = -1i32;
    if cell.num_words as usize <= NA_WORDS_PER_HASH {
        let mut start = 0usize;
        let mut word_index = 0usize;
        while word_index < cell.num_words as usize {
            let num_offsets = if cell.num_offsets[word_index] > 0 {
                cell.num_offsets[word_index] as usize
            } else {
                0
            };
            if cell.words[word_index] == index {
                if num_offsets > 0 {
                    let mut j = 0usize;
                    while j < num_offsets {
                        offset_pairs.push(OffsetPair {
                            query_offset: cell.offsets[start + j],
                            subject_offset: s_off,
                        });
                        num_hits += 1;
                        j += 1;
                    }
                } else {
                    cursor = cell.offsets[0];
                }
                break;
            }
            start += num_offsets;
            word_index += 1;
        }
    } else {
        cursor = cell.offsets[0];
    }

    if cursor >= 0 {
        let mut overflow_index = cursor as usize;
        let mut k = 0;
        while k < cell.num_words {
            if overflow_index + 1 >= lookup.overflow.len() {
                break;
            }
            let word = lookup.overflow[overflow_index] as u32;
            let num_offsets = lookup.overflow[overflow_index + 1];
            let num_offsets_usize = if num_offsets > 0 {
                num_offsets as usize
            } else {
                0
            };
            if word != index {
                overflow_index += num_offsets_usize + 2;
                k += 1;
                continue;
            }
            overflow_index += 2;
            let mut i = 0usize;
            while i < num_offsets_usize && overflow_index + i < lookup.overflow.len() {
                offset_pairs.push(OffsetPair {
                    query_offset: lookup.overflow[overflow_index + i],
                    subject_offset: s_off,
                });
                num_hits += 1;
                i += 1;
            }
            break;
        }
    }

    num_hits
}

/// Port of NCBI static `s_BlastNaHashScanSubject_Any`
/// (`blast_nascan.c:2860`) for the represented 16-base NaHash scan.
pub fn s_blast_na_hash_scan_subject_any(
    lookup: &BlastNaHashLookupTable,
    subject_packed: &[u8],
    subject_length: i32,
    scan_start: i32,
    scan_end: i32,
) -> Option<Vec<OffsetPair>> {
    if lookup.lut_word_length != 16 || lookup.scan_step <= 0 || subject_length < 16 {
        return None;
    }
    if scan_start < 0 || scan_end < scan_start || scan_end + lookup.lut_word_length > subject_length
    {
        return None;
    }
    if subject_length as usize > subject_packed.len().saturating_mul(4) {
        return None;
    }

    let mut offset_pairs = Vec::new();
    let mut s_off = scan_start;
    while s_off <= scan_end {
        let byte = (s_off / crate::stat::COMPRESSION_RATIO as i32) as usize;
        let phase = s_off % crate::stat::COMPRESSION_RATIO as i32;
        let index = if phase == 0 {
            let chunk = subject_packed.get(byte..byte + 4)?;
            u32::from_be_bytes([chunk[0], chunk[1], chunk[2], chunk[3]])
        } else {
            let chunk = subject_packed.get(byte..byte + 5)?;
            let w = (u64::from(chunk[0]) << 32)
                | (u64::from(chunk[1]) << 24)
                | (u64::from(chunk[2]) << 16)
                | (u64::from(chunk[3]) << 8)
                | u64::from(chunk[4]);
            let shift = 2 * (crate::stat::COMPRESSION_RATIO as i32 - phase);
            ((w >> shift) & u64::from(u32::MAX)) as u32
        };
        let _ = s_blast_na_hash_lookup_retrieve_hits(lookup, index, s_off, &mut offset_pairs);
        s_off = s_off.saturating_add(lookup.scan_step);
    }

    Some(offset_pairs)
}

/// Port of NCBI `s_RemovePolyAWords` (`blast_nalookup.c:903`).
pub fn s_remove_poly_a_words(mb_lt: &mut MbLookupTable) -> i16 {
    let word_size = mb_lt.lut_word_length.max(0);
    if let Some(slot) = mb_lt.hashtable.get_mut(0) {
        *slot = 0;
    }
    if word_size < 31 {
        let all_t = ((1u64 << (2 * word_size as u32)) - 1) as usize;
        if let Some(slot) = mb_lt.hashtable.get_mut(all_t) {
            *slot = 0;
        }
    }

    if word_size < 16 {
        return 0;
    }

    for i in 1..4u64 {
        for k in 0..word_size.min(16) {
            let word = (i << (k * 2)) as usize;
            if let Some(slot) = mb_lt.hashtable.get_mut(word) {
                *slot = 0;
            }
        }
    }

    for i in 0..3u64 {
        for k in 0..word_size.min(16) {
            let shift = (k * 2) as u32;
            let word = ((u32::MAX ^ (3u32 << shift)) | ((i as u32) << shift)) as usize;
            if let Some(slot) = mb_lt.hashtable.get_mut(word) {
                *slot = 0;
            }
        }
    }

    0
}

/// Port of NCBI `s_FillPV` (`blast_nalookup.c:812`) for contiguous
/// megablast database-word filtering.
pub fn s_fill_pv(
    query_sequence: &[u8],
    locations: &[crate::util::SSeqRange],
    mb_lt: &mut MbLookupTable,
    _lookup_options: Option<&crate::options::LookupTableOptions>,
) -> i16 {
    let lut_word_length = mb_lt.lut_word_length.max(0) as usize;
    let full_word_size = mb_lt.word_length.max(0) as usize;
    let hashsize = mb_lt.hashtable.len().max(
        1usize
            .checked_shl((2 * lut_word_length) as u32)
            .unwrap_or(0),
    );
    let lut_mask = (hashsize - 1) as u32;
    let pv_array_bts = mb_lt.pv_array_bts.max(0) as usize;

    for loc in locations {
        if loc.left < 0 || loc.right < loc.left {
            continue;
        }
        let from = loc.left as usize;
        let right = loc.right as usize;
        let loc_len = right - from + 1;
        if full_word_size > loc_len || lut_word_length == 0 {
            continue;
        }
        let Some(window) = query_sequence.get(from..=right) else {
            continue;
        };

        let mut ecode = 0u32;
        let mut valid_bases = 0usize;
        for &val in window {
            if (val & BLAST2NA_MASK) != 0 {
                ecode = 0;
                valid_bases = 0;
                continue;
            }

            ecode = ((ecode << BITS_PER_NUC) & lut_mask) + u32::from(val);
            valid_bases += 1;
            if valid_bases >= lut_word_length {
                set_pv_bit(&mut mb_lt.pv_array, ecode as usize, pv_array_bts);
            }
        }
    }

    0
}

/// Port of NCBI `s_FillContigMBTable` (`blast_nalookup.c:947`).
pub fn s_fill_contig_mb_table(
    query_sequence: &[u8],
    locations: &[crate::util::SSeqRange],
    mb_lt: &mut MbLookupTable,
    lookup_options: &crate::options::LookupTableOptions,
    counts: Option<&[u8]>,
) -> i16 {
    let lut_word_length = mb_lt.lut_word_length.max(0) as usize;
    let full_word_size = mb_lt.word_length.max(0) as usize;
    if lut_word_length == 0 || mb_lt.hashtable.is_empty() {
        return -1;
    }

    mb_lt.next_pos = vec![0; query_sequence.len() + 1];
    let hashsize = mb_lt.hashtable.len();
    let lut_mask = (hashsize - 1) as u32;
    let pv_array_bts = mb_lt.pv_array_bts.max(0) as usize;
    let compression_factor = 2048usize;
    let mut helper_array = vec![0u32; (hashsize / compression_factor).max(1)];

    if lookup_options.db_filter {
        mb_lt.pv_array.fill(0);
    }

    let scan_step = if lookup_options.stride > 0 {
        lookup_options.stride as usize
    } else {
        1
    };

    for loc in locations {
        if loc.left < 0 || loc.right < loc.left {
            continue;
        }
        let from = loc.left as usize;
        let right = loc.right as usize;
        let loc_len = right - from + 1;
        if full_word_size > loc_len || from >= query_sequence.len() {
            continue;
        }
        let last_start = right.saturating_sub(lut_word_length - 1);
        if last_start < from {
            continue;
        }

        let mut start = from;
        while start <= last_start {
            let Some(word_slice) = query_sequence.get(start..start + lut_word_length) else {
                break;
            };
            let mut ecode = 0u32;
            let mut ambiguous = false;
            for &val in word_slice {
                if (val & BLAST2NA_MASK) != 0 {
                    ambiguous = true;
                    break;
                }
                ecode = ((ecode << BITS_PER_NUC) & lut_mask) + u32::from(val);
            }
            if ambiguous {
                start += 1;
                continue;
            }

            if lookup_options.db_filter {
                let Some(counts) = counts else {
                    return -1;
                };
                let count_index = (ecode / 2) as usize;
                let Some(&packed_count) = counts.get(count_index) else {
                    return -1;
                };
                let count = if (ecode & 1) == 0 {
                    packed_count >> 4
                } else {
                    packed_count & 0x0f
                };
                if count >= lookup_options.max_db_word_count {
                    start += 1;
                    continue;
                }
            }

            let ecode_index = ecode as usize;
            let query_index = start + 1;
            if mb_lt.hashtable[ecode_index] == 0 {
                set_pv_bit(&mut mb_lt.pv_array, ecode_index, pv_array_bts);
            } else {
                let helper_index = (ecode_index / compression_factor).min(helper_array.len() - 1);
                helper_array[helper_index] += 1;
            }
            mb_lt.next_pos[query_index] = mb_lt.hashtable[ecode_index];
            mb_lt.hashtable[ecode_index] = query_index as i32;

            start += scan_step;
        }
    }

    if crate::program::blast_program_is_mapping(lookup_options.program_number) {
        s_remove_poly_a_words(mb_lt);
    }

    mb_lt.longest_chain = helper_array.into_iter().max().unwrap_or(2).max(2) as i32;
    0
}

fn s_insert_mb_word(
    mb_lt: &mut MbLookupTable,
    helper_array: &mut [u32],
    ecode: usize,
    query_index: usize,
    use_second_table: bool,
) -> i16 {
    let (hashtable, next_pos) = if use_second_table {
        (&mut mb_lt.hashtable2, &mut mb_lt.next_pos2)
    } else {
        (&mut mb_lt.hashtable, &mut mb_lt.next_pos)
    };
    if ecode >= hashtable.len() || query_index >= next_pos.len() {
        return -1;
    }

    if hashtable[ecode] == 0 {
        set_pv_bit(
            &mut mb_lt.pv_array,
            ecode,
            mb_lt.pv_array_bts.max(0) as usize,
        );
    } else if !helper_array.is_empty() {
        let helper_index = (ecode / 2048).min(helper_array.len() - 1);
        helper_array[helper_index] += 1;
    }
    next_pos[query_index] = hashtable[ecode];
    hashtable[ecode] = query_index as i32;
    0
}

/// Port of NCBI `s_FillDiscMBTable` (`blast_nalookup.c:656`).
pub fn s_fill_disc_mb_table(
    query_sequence: &[u8],
    locations: &[crate::util::SSeqRange],
    mb_lt: &mut MbLookupTable,
    lookup_options: &crate::options::LookupTableOptions,
) -> i16 {
    let word_size = lookup_options.word_size;
    let template_length = lookup_options.mb_template_length;
    if word_size <= 0 || template_length <= 0 || template_length < word_size {
        return -1;
    }

    let Ok(template_length_u8) = u8::try_from(template_length) else {
        return -1;
    };
    let word_type = match lookup_options.mb_template_type {
        0 => DiscWordType::Coding,
        1 => DiscWordType::Optimal,
        2 => DiscWordType::TwoTemplates,
        _ => return -1,
    };
    let template_type = s_get_disc_template_type(word_size, template_length_u8, word_type);
    if template_type == DiscTemplateType::Contiguous {
        return -1;
    }
    let second_template_type = if word_type == DiscWordType::TwoTemplates {
        s_get_disc_template_type(word_size, template_length_u8, DiscWordType::Optimal)
    } else {
        DiscTemplateType::Contiguous
    };

    mb_lt.lut_word_length = word_size;
    mb_lt.word_length = template_length;
    mb_lt.discontiguous = true;
    mb_lt.template_length = template_length;
    mb_lt.template_type = template_type;
    mb_lt.two_templates = word_type == DiscWordType::TwoTemplates;
    mb_lt.second_template_type = second_template_type;
    mb_lt.scan_step = 1;
    mb_lt.next_pos = vec![0; query_sequence.len() + 1];
    if word_type == DiscWordType::TwoTemplates {
        mb_lt.next_pos2 = vec![0; query_sequence.len() + 1];
        if mb_lt.hashtable2.len() != mb_lt.hashtable.len() {
            mb_lt.hashtable2 = vec![0; mb_lt.hashtable.len()];
        }
    }
    mb_lt.pv_array.fill(0);

    let mut helper_array = vec![0u32; (mb_lt.hashtable.len() / 2048).max(1)];
    let mut helper_array2 = if word_type == DiscWordType::TwoTemplates {
        vec![0u32; (mb_lt.hashtable.len() / 2048).max(1)]
    } else {
        Vec::new()
    };
    let template_length_usize = template_length as usize;
    let template_mask = if template_length >= 32 {
        u64::MAX
    } else {
        (1u64 << (2 * template_length as u32)) - 1
    };

    for loc in locations {
        if loc.left < 0 || loc.right < loc.left {
            continue;
        }
        let from = loc.left as usize;
        let right = loc.right as usize;
        if right >= query_sequence.len() || right + 1 - from < template_length_usize {
            continue;
        }

        let mut accum = 0u64;
        let mut valid_bases = 0usize;
        for pos in from..=right {
            let base = query_sequence[pos];
            if (base & BLAST2NA_MASK) != 0 {
                accum = 0;
                valid_bases = 0;
                continue;
            }
            accum = ((accum << BITS_PER_NUC) | u64::from(base)) & template_mask;
            valid_bases += 1;
            if valid_bases < template_length_usize {
                continue;
            }

            let query_index = pos + 2 - template_length_usize;
            let ecode1 = compute_discontiguous_index(accum, template_type) as usize;
            if s_insert_mb_word(mb_lt, &mut helper_array, ecode1, query_index, false) != 0 {
                return -1;
            }
            if word_type == DiscWordType::TwoTemplates {
                let ecode2 = compute_discontiguous_index(accum, second_template_type) as usize;
                if s_insert_mb_word(mb_lt, &mut helper_array2, ecode2, query_index, true) != 0 {
                    return -1;
                }
            }
        }
    }

    let longest_chain = helper_array.into_iter().max().unwrap_or(2).max(2) + 1;
    mb_lt.longest_chain = longest_chain as i32;
    if word_type == DiscWordType::TwoTemplates {
        let second_longest_chain = helper_array2.into_iter().max().unwrap_or(2).max(2) + 1;
        mb_lt.longest_chain += second_longest_chain as i32;
    }
    0
}

fn s_mb_disc_word_scan_subject_template(
    mb_lt: &MbLookupTable,
    subject: &[u8],
    subject_len: usize,
    max_hits: usize,
    scan_range: usize,
    template_length: usize,
    template_type: DiscTemplateType,
) -> Vec<OffsetPair> {
    if max_hits == 0
        || subject_len < template_length
        || template_type == DiscTemplateType::Contiguous
    {
        return Vec::new();
    }

    let last_start = scan_range
        .min(subject_len)
        .min(subject_len.saturating_sub(template_length) + 1);
    let template_mask = if template_length >= 32 {
        u64::MAX
    } else {
        (1u64 << (BITS_PER_NUC * template_length as u32)) - 1
    };
    let mut accum = 0u64;
    let mut hits = Vec::new();

    for pos in 0..subject_len {
        accum = ((accum << BITS_PER_NUC)
            | u64::from(crate::encoding::ncbi2na_base_at(subject, pos)))
            & template_mask;
        if pos + 1 < template_length {
            continue;
        }

        let s_off = pos + 1 - template_length;
        if s_off >= last_start {
            break;
        }

        let ecode = compute_discontiguous_index(accum, template_type) as i64;
        if s_blast_mb_lookup_has_hits(mb_lt, ecode) != 0 {
            let before = hits.len();
            s_blast_mb_lookup_retrieve(mb_lt, ecode, &mut hits, s_off as i32);
            if hits.len() >= max_hits {
                hits.truncate(max_hits);
                break;
            }
            if hits.len() == before {
                continue;
            }
        }
    }

    hits
}

/// Port-shaped dispatcher for NCBI discontiguous megablast subject scanners
/// (`s_MB_DiscWordScanSubject_*` in `blast_nascan.c`).
pub fn s_mb_disc_word_scan_subject(
    mb_lt: &MbLookupTable,
    subject: &[u8],
    subject_len: usize,
    scan_start: usize,
    scan_end_inclusive: usize,
) -> Vec<OffsetPair> {
    if !mb_lt.discontiguous
        || mb_lt.template_length <= 0
        || mb_lt.template_type == DiscTemplateType::Contiguous
        || subject_len < mb_lt.template_length as usize
        || scan_start > scan_end_inclusive
    {
        return Vec::new();
    }

    let template_length = mb_lt.template_length as usize;
    let last_start = scan_end_inclusive
        .min(subject_len.saturating_sub(template_length))
        .saturating_add(1);
    if scan_start >= last_start {
        return Vec::new();
    }

    let template_mask = if template_length >= 32 {
        u64::MAX
    } else {
        (1u64 << (BITS_PER_NUC * template_length as u32)) - 1
    };
    let mut accum = 0u64;
    let mut hits = Vec::new();

    for pos in 0..subject_len {
        accum = ((accum << BITS_PER_NUC)
            | u64::from(crate::encoding::ncbi2na_base_at(subject, pos)))
            & template_mask;
        if pos + 1 < template_length {
            continue;
        }

        let s_off = pos + 1 - template_length;
        if s_off < scan_start {
            continue;
        }
        if s_off >= last_start {
            break;
        }

        let ecode = compute_discontiguous_index(accum, mb_lt.template_type) as i64;
        crate::lookup::s_blast_mb_lookup_retrieve(mb_lt, ecode, &mut hits, s_off as i32);
        if mb_lt.two_templates && mb_lt.second_template_type != DiscTemplateType::Contiguous {
            let ecode2 = compute_discontiguous_index(accum, mb_lt.second_template_type) as i64;
            crate::lookup::s_blast_mb_lookup_retrieve2(mb_lt, ecode2, &mut hits, s_off as i32);
        }
    }

    hits
}

/// Port of NCBI static `s_MB_DiscWordScanSubject_11_18_1`
/// (`blast_nascan.c:2388`) over packed 2-bit subject bytes.
pub fn s_mb_disc_word_scan_subject_11_18_1(
    mb_lt: &MbLookupTable,
    subject: &[u8],
    subject_len: usize,
    max_hits: usize,
    scan_range: usize,
) -> Vec<OffsetPair> {
    debug_assert_eq!(mb_lt.lut_word_length, 11);
    debug_assert_eq!(mb_lt.word_length, 18);
    debug_assert_eq!(mb_lt.template_type, DiscTemplateType::Template11_18Coding);

    let mut hits = Vec::new();
    const TEMPLATE_LENGTH: usize = 18;
    if subject_len < TEMPLATE_LENGTH {
        return hits;
    }
    let last_start = scan_range
        .min(subject_len.saturating_sub(TEMPLATE_LENGTH))
        .saturating_add(1);
    let max_hits = max_hits.saturating_sub(mb_lt.longest_chain.max(0) as usize);
    let template_mask = (1u64 << (BITS_PER_NUC * TEMPLATE_LENGTH as u32)) - 1;
    let mut accum = 0u64;

    for pos in 0..subject_len {
        accum = ((accum << BITS_PER_NUC)
            | u64::from(crate::encoding::ncbi2na_base_at(subject, pos)))
            & template_mask;
        if pos + 1 < TEMPLATE_LENGTH {
            continue;
        }

        let s_off = pos + 1 - TEMPLATE_LENGTH;
        if s_off >= last_start {
            break;
        }
        let lo = accum as u32;
        let hi = (accum >> 32) as u32;
        let index = (lo & 0x00000003)
            | ((lo & 0x000000f0) >> 2)
            | ((lo & 0x00003c00) >> 4)
            | ((lo & 0x00030000) >> 6)
            | ((lo & 0x03c00000) >> 10)
            | ((lo & 0xf0000000) >> 12)
            | ((hi & 0x0000000c) << 18);
        if s_blast_mb_lookup_has_hits(mb_lt, i64::from(index)) != 0 {
            if hits.len() >= max_hits {
                break;
            }
            s_blast_mb_lookup_retrieve(mb_lt, i64::from(index), &mut hits, s_off as i32);
        }
    }

    hits
}

/// Port of NCBI static `s_MB_DiscWordScanSubject_11_21_1`
/// (`blast_nascan.c:2501`) over packed 2-bit subject bytes.
pub fn s_mb_disc_word_scan_subject_11_21_1(
    mb_lt: &MbLookupTable,
    subject: &[u8],
    subject_len: usize,
    max_hits: usize,
    scan_range: usize,
) -> Vec<OffsetPair> {
    debug_assert_eq!(mb_lt.lut_word_length, 11);
    debug_assert_eq!(mb_lt.word_length, 21);
    debug_assert_eq!(mb_lt.template_type, DiscTemplateType::Template11_21Coding);

    let mut hits = Vec::new();
    const TEMPLATE_LENGTH: usize = 21;
    if subject_len < TEMPLATE_LENGTH {
        return hits;
    }
    let last_start = scan_range
        .min(subject_len.saturating_sub(TEMPLATE_LENGTH))
        .saturating_add(1);
    let max_hits = max_hits.saturating_sub(mb_lt.longest_chain.max(0) as usize);
    let template_mask = (1u64 << (BITS_PER_NUC * TEMPLATE_LENGTH as u32)) - 1;
    let mut accum = 0u64;

    for pos in 0..subject_len {
        accum = ((accum << BITS_PER_NUC)
            | u64::from(crate::encoding::ncbi2na_base_at(subject, pos)))
            & template_mask;
        if pos + 1 < TEMPLATE_LENGTH {
            continue;
        }

        let s_off = pos + 1 - TEMPLATE_LENGTH;
        if s_off >= last_start {
            break;
        }
        let lo = accum as u32;
        let hi = (accum >> 32) as u32;
        let index = (lo & 0x00000003)
            | ((lo & 0x000000f0) >> 2)
            | ((lo & 0x00000c00) >> 4)
            | ((lo & 0x000f0000) >> 8)
            | ((lo & 0x00c00000) >> 10)
            | ((lo & 0xf0000000) >> 14)
            | ((hi & 0x0000000c) << 16)
            | ((hi & 0x00000300) << 12);
        if s_blast_mb_lookup_has_hits(mb_lt, i64::from(index)) != 0 {
            if hits.len() >= max_hits {
                break;
            }
            s_blast_mb_lookup_retrieve(mb_lt, i64::from(index), &mut hits, s_off as i32);
        }
    }

    hits
}

/// Port of NCBI `s_MBCountWordsInSubject_16_1` (`blast_nalookup.c:1127`)
/// over packed 2-bit subject bytes.
pub fn s_mb_count_words_in_subject_16_1(
    sequence: Option<&[u8]>,
    sequence_length: i32,
    mb_lt: Option<&MbLookupTable>,
    counts: Option<&mut [u8]>,
    max_word_count: u8,
) -> i16 {
    let (Some(sequence), Some(mb_lt), Some(counts)) = (sequence, mb_lt, counts) else {
        return -1;
    };
    if mb_lt.pv_array.is_empty() || sequence_length < mb_lt.lut_word_length {
        return -1;
    }

    let k_num_words = sequence_length - mb_lt.lut_word_length;
    let max_read_index = if k_num_words > 0 {
        ((k_num_words - 1) / 4 + 4) as usize
    } else {
        3
    };
    if sequence.len() <= max_read_index {
        return -1;
    }

    let mask = mb_lt.hashtable.len().saturating_sub(1) as u64;
    let pv_array_bts = mb_lt.pv_array_bts.max(0) as usize;
    let mut shift = 8u32;
    let mut w = (u64::from(sequence[0]) << 24)
        | (u64::from(sequence[1]) << 16)
        | (u64::from(sequence[2]) << 8)
        | u64::from(sequence[3]);

    for i in 0..k_num_words {
        if i % 4 == 0 {
            shift = 8;
            w = (w << 8) | u64::from(sequence[(i / 4 + 4) as usize]);
        } else {
            shift -= 2;
        }

        let word = (w >> shift) & mask;
        if !test_pv_bit(&mb_lt.pv_array, word, pv_array_bts) {
            continue;
        }

        let index = (word / 2) as usize;
        let Some(count) = counts.get_mut(index) else {
            return -1;
        };
        if (word & 1) != 0 {
            if (*count & 0x0f) < max_word_count {
                *count += 1;
            }
        } else if (*count >> 4) < max_word_count {
            *count += 1 << 4;
        }
    }

    0
}

/// Port of NCBI `s_ScanSubjectForWordCounts` (`blast_nalookup.c:1189`).
pub fn s_scan_subject_for_word_counts(
    seq_src: Option<&mut dyn crate::seqsrc::BlastSeqSource>,
    mb_lt: Option<&MbLookupTable>,
    counts: Option<&mut [u8]>,
    max_word_count: u8,
) -> i16 {
    let (Some(seq_src), Some(mb_lt), Some(counts)) = (seq_src, mb_lt, counts) else {
        return -1;
    };
    if mb_lt.pv_array.is_empty() {
        return -1;
    }

    crate::seqsrc::blast_seq_src_reset_chunk_iterator(seq_src);
    let oids: Vec<i32> = seq_src.iter_oids().collect();
    for oid in oids {
        if let Some(seq_data) = seq_src.get_sequence(&crate::seqsrc::GetSeqArg {
            oid,
            encoding: crate::seqsrc::SeqEncoding::Protein,
            ..crate::seqsrc::GetSeqArg::default()
        }) {
            let status = s_mb_count_words_in_subject_16_1(
                Some(&seq_data.sequence),
                seq_data.length,
                Some(mb_lt),
                Some(&mut *counts),
                max_word_count,
            );
            if status != 0 {
                return status;
            }
        }
    }

    0
}

/// Rust ownership translation of NCBI `BlastMBLookupTableNew`
/// (`blast_nalookup.c:1227`).
pub fn blast_mb_lookup_table_new(
    query_sequence: &[u8],
    locations: &[crate::util::SSeqRange],
    mb_lt_ptr: &mut Option<MbLookupTable>,
    lookup_options: &crate::options::LookupTableOptions,
    approx_table_entries: i32,
    lut_width: i32,
    seq_src: Option<&mut dyn crate::seqsrc::BlastSeqSource>,
) -> i16 {
    *mb_lt_ptr = None;
    let Some(hashsize) = nucleotide_backbone_size(lut_width) else {
        return -1;
    };
    if lut_width < 9 || lookup_options.word_size < lut_width {
        return -1;
    }

    let target_pv_size = 131_072usize;
    let small_query_cutoff = 15_000;
    let large_query_cutoff = 800_000;
    let mut pv_size = if lut_width <= 12 {
        if hashsize <= 8 * target_pv_size {
            hashsize >> crate::stat::PV_ARRAY_BTS
        } else {
            target_pv_size / crate::stat::PV_ARRAY_BYTES
        }
    } else {
        target_pv_size * 64 / crate::stat::PV_ARRAY_BYTES
    }
    .max(1);

    if !lookup_options.db_filter
        && (approx_table_entries <= small_query_cutoff
            || approx_table_entries >= large_query_cutoff)
    {
        pv_size = (pv_size / 2).max(1);
    }

    let pv_array_bts = {
        let entries_per_bit = (hashsize / pv_size.max(1)).max(1);
        entries_per_bit.trailing_zeros() as i32
    };

    let mut mb_lt = MbLookupTable {
        word_length: lookup_options.word_size,
        lut_word_length: lut_width,
        discontiguous: false,
        template_length: 0,
        template_type: DiscTemplateType::Contiguous,
        two_templates: false,
        second_template_type: DiscTemplateType::Contiguous,
        hashtable: vec![0; hashsize],
        hashtable2: Vec::new(),
        next_pos: Vec::new(),
        next_pos2: Vec::new(),
        pv_array: vec![0; pv_size],
        pv_array_bts,
        longest_chain: 0,
        scan_step: 0,
    };

    let mut counts = if lookup_options.db_filter {
        Some(vec![0u8; hashsize / 2])
    } else {
        None
    };

    if lookup_options.db_filter {
        let Some(seq_src) = seq_src else {
            return -1;
        };
        if s_fill_pv(query_sequence, locations, &mut mb_lt, Some(lookup_options)) != 0 {
            return -1;
        }
        let Some(counts_ref) = counts.as_mut() else {
            return -1;
        };
        if s_scan_subject_for_word_counts(
            Some(seq_src),
            Some(&mb_lt),
            Some(counts_ref.as_mut_slice()),
            lookup_options.max_db_word_count,
        ) != 0
        {
            return -1;
        }
    }

    let status = if lookup_options.mb_template_length > 0 {
        mb_lt.scan_step = 1;
        s_fill_disc_mb_table(query_sequence, locations, &mut mb_lt, lookup_options)
    } else {
        mb_lt.scan_step = lookup_options.word_size - lut_width + 1;
        if mb_lt.scan_step <= 0 {
            -1
        } else {
            s_fill_contig_mb_table(
                query_sequence,
                locations,
                &mut mb_lt,
                lookup_options,
                counts.as_deref(),
            )
        }
    };

    if status != 0 {
        return status;
    }

    *mb_lt_ptr = Some(mb_lt);
    0
}

/// 1-1 ownership translation of `BlastMBLookupTableDestruct`.
pub fn blast_mb_lookup_table_destruct(lookup: &mut Option<MbLookupTable>) -> Option<MbLookupTable> {
    *lookup = None;
    None
}

/// Conservative ownership port of NCBI `RPSLookupTableNew`
/// (`blast_aalookup.c:124`).
///
/// Audit caveat: the C function builds directly from memory-mapped RPS
/// database structures and preserves raw row pointers into the PSSM payload.
/// Rust models the represented ownership only: backbone cells, a presence
/// vector, cloned PSSM rows, and bucket storage for offset pairs.
pub fn rps_lookup_table_new(
    info: Option<&BlastRpsInfo>,
    lut: &mut Option<BlastRpsLookupTable>,
) -> i16 {
    *lut = None;
    let Some(info) = info else {
        return -1;
    };

    let alphabet_size = info.alphabet_size;
    if alphabet_size <= 0 || info.wordsize <= 0 {
        return -1;
    }

    let charsize = crate::util::ilog2(alphabet_size as i64) + 1;
    let backbone_size = info.rps_backbone.len() as i32;
    let mut pv = vec![0u32; ((backbone_size.max(0) as usize) >> crate::stat::PV_ARRAY_BTS) + 1];
    let mut bucket_array = Vec::new();

    for (index, cell) in info.rps_backbone.iter().enumerate() {
        if cell.num_used <= 0 {
            continue;
        }

        let word = index >> crate::stat::PV_ARRAY_BTS;
        let bit = (index as u32) & crate::stat::PV_ARRAY_MASK;
        if let Some(slot) = pv.get_mut(word) {
            *slot |= 1u32 << bit;
        }

        let num_used = cell.num_used.max(cell.offset_pairs.len() as i32);
        let num_alloc = num_used.max(1000);
        bucket_array.push(RpsBucket {
            num_alloc,
            num_used,
            offset_pairs: cell.offset_pairs.clone(),
        });
    }

    *lut = Some(BlastRpsLookupTable {
        alphabet_size,
        wordsize: info.wordsize,
        charsize,
        backbone_size,
        rps_backbone: info.rps_backbone.clone(),
        pv,
        rps_pssm: info.rps_pssm.clone(),
        num_buckets: bucket_array.len() as i32,
        bucket_array,
    });
    0
}

fn rps_profile_index_for_offset(start_offsets: &[i32], offset: i32) -> Option<usize> {
    if offset < 0 || start_offsets.len() < 2 {
        return None;
    }
    let offset = offset as usize;
    let starts: Vec<usize> = start_offsets
        .iter()
        .map(|&value| usize::try_from(value).ok())
        .collect::<Option<Vec<_>>>()?;
    let index = starts.partition_point(|&start| start <= offset);
    if index == 0 {
        return None;
    }
    let profile_index = index - 1;
    if profile_index + 1 >= starts.len() || offset >= starts[profile_index + 1] {
        return None;
    }
    Some(profile_index)
}

fn rps_lookup_magic_and_byte_order(
    bytes: &[u8],
) -> Option<(i32, crate::hspstream::RpsNativeByteOrder)> {
    let magic_bytes: [u8; 4] = bytes.get(..4)?.try_into().ok()?;
    let little = i32::from_le_bytes(magic_bytes);
    if matches!(
        little,
        crate::hspstream::RPS_MAGIC_NUM | crate::hspstream::RPS_MAGIC_NUM_28
    ) {
        return Some((little, crate::hspstream::RpsNativeByteOrder::LittleEndian));
    }
    let big = i32::from_be_bytes(magic_bytes);
    if matches!(
        big,
        crate::hspstream::RPS_MAGIC_NUM | crate::hspstream::RPS_MAGIC_NUM_28
    ) {
        return Some((big, crate::hspstream::RpsNativeByteOrder::BigEndian));
    }
    None
}

fn rps_lookup_native_i32_values(
    bytes: &[u8],
    byte_order: crate::hspstream::RpsNativeByteOrder,
) -> Option<Vec<i32>> {
    if bytes.len() % std::mem::size_of::<i32>() != 0 {
        return None;
    }
    Some(
        bytes
            .chunks_exact(std::mem::size_of::<i32>())
            .map(|chunk| {
                let raw: [u8; 4] = chunk.try_into().expect("chunks_exact width");
                match byte_order {
                    crate::hspstream::RpsNativeByteOrder::LittleEndian => i32::from_le_bytes(raw),
                    crate::hspstream::RpsNativeByteOrder::BigEndian => i32::from_be_bytes(raw),
                }
            })
            .collect(),
    )
}

fn rps_alphabet_size_from_magic(magic_number: i32) -> Option<i32> {
    match magic_number {
        crate::hspstream::RPS_MAGIC_NUM => Some(26),
        crate::hspstream::RPS_MAGIC_NUM_28 => Some(28),
        _ => None,
    }
}

/// blast-rs: native `.loo` lookup-table parser for the owned RPS boundary;
/// not a direct NCBI C port.
///
/// NCBI memory maps `BlastRPSLookupFileHeader`, followed by
/// `(backbone_size + 1)` compact backbone cells and an overflow `Int4` array.
/// Each compact cell stores profile offsets to the last residue in a word; the
/// existing Rust scanner applies the same first-residue correction during
/// lookup-table scans.
pub fn rps_lookup_info_from_native_lookup_bytes(
    bytes: &[u8],
    rps_pssm: Vec<Vec<i32>>,
) -> Result<BlastRpsInfo, i16> {
    let (magic_number, byte_order) = rps_lookup_magic_and_byte_order(bytes).ok_or(-1i16)?;
    let alphabet_size = rps_alphabet_size_from_magic(magic_number).ok_or(-1i16)?;
    let values = rps_lookup_native_i32_values(bytes, byte_order).ok_or(-1i16)?;
    if values.len() < RPS_LOOKUP_HEADER_WORDS {
        return Err(-1);
    }

    let num_lookup_tables = values[1];
    let num_hits = values[2];
    let num_filled_backbone_cells = values[3];
    let overflow_hits = values[4];
    let start_of_backbone = values[8];
    let end_of_overflow = values[9];
    if num_lookup_tables != 1
        || num_hits < 0
        || num_filled_backbone_cells < 0
        || overflow_hits < 0
        || start_of_backbone < 0
        || end_of_overflow < start_of_backbone
        || start_of_backbone as usize % std::mem::size_of::<i32>() != 0
        || end_of_overflow as usize % std::mem::size_of::<i32>() != 0
        || end_of_overflow as usize > bytes.len()
    {
        return Err(-1);
    }

    let charsize = crate::util::ilog2(alphabet_size as i64) + 1;
    let backbone_bits = RPS_LOOKUP_WORDSIZE.checked_mul(charsize).ok_or(-1i16)?;
    if !(1..usize::BITS as i32).contains(&backbone_bits) {
        return Err(-1);
    }
    let backbone_size = 1usize << backbone_bits as usize;
    let backbone_cells = backbone_size.checked_add(1).ok_or(-1i16)?;
    let backbone_words = backbone_cells
        .checked_mul(RPS_LOOKUP_BACKBONE_CELL_WORDS)
        .ok_or(-1i16)?;
    let backbone_start = start_of_backbone as usize / std::mem::size_of::<i32>();
    let overflow_start = backbone_start.checked_add(backbone_words).ok_or(-1i16)?;
    let overflow_end = overflow_start
        .checked_add(overflow_hits as usize)
        .ok_or(-1i16)?;
    if overflow_end > end_of_overflow as usize / std::mem::size_of::<i32>()
        || overflow_end > values.len()
    {
        return Err(-1);
    }

    let overflow_values = &values[overflow_start..overflow_end];
    let mut rps_backbone = vec![RpsBackboneCell::default(); backbone_size];
    let mut observed_hits = 0i32;
    let mut observed_filled = 0i32;
    for (index, cell) in rps_backbone.iter_mut().enumerate() {
        let cell_start = backbone_start + index * RPS_LOOKUP_BACKBONE_CELL_WORDS;
        let num_used = values[cell_start];
        if num_used < 0 {
            return Err(-1);
        }
        if num_used == 0 {
            continue;
        }

        observed_filled += 1;
        observed_hits = observed_hits.checked_add(num_used).ok_or(-1i16)?;
        let entries = &values[cell_start + 1..cell_start + 1 + RPS_HITS_PER_CELL];
        let mut offset_pairs = Vec::with_capacity(num_used as usize);
        if num_used as usize <= RPS_HITS_PER_CELL {
            offset_pairs.extend(entries.iter().take(num_used as usize).map(|&query_offset| {
                OffsetPair {
                    query_offset,
                    subject_offset: 0,
                }
            }));
        } else {
            if entries[1] < 0 || entries[1] as usize % std::mem::size_of::<i32>() != 0 {
                return Err(-1);
            }
            let overflow_index = entries[1] as usize / std::mem::size_of::<i32>();
            let overflow_count = num_used as usize - 1;
            let overflow_limit = overflow_index.checked_add(overflow_count).ok_or(-1i16)?;
            let overflow_slice = overflow_values
                .get(overflow_index..overflow_limit)
                .ok_or(-1i16)?;
            offset_pairs.push(OffsetPair {
                query_offset: entries[0],
                subject_offset: 0,
            });
            offset_pairs.extend(overflow_slice.iter().map(|&query_offset| OffsetPair {
                query_offset,
                subject_offset: 0,
            }));
        }
        cell.num_used = num_used;
        cell.offset_pairs = offset_pairs;
    }

    if observed_hits != num_hits || observed_filled != num_filled_backbone_cells {
        return Err(-1);
    }

    Ok(BlastRpsInfo {
        alphabet_size,
        wordsize: RPS_LOOKUP_WORDSIZE,
        rps_backbone,
        rps_pssm,
    })
}

fn blast_rps_word_finder_scan_to_hsp_list_with_query_shift(
    scan: &RpsWordFinderScan,
    profile_oid: i32,
    hsp_max: i32,
    query_shift: i32,
) -> Option<crate::hspstream::HspList> {
    if scan.init_hits.is_empty() {
        return None;
    }

    let mut hsp_list = crate::hspstream::blast_hsp_list_new(hsp_max);
    hsp_list.oid = profile_oid;
    for hit in &scan.init_hits {
        let q_start = hit.q_start - query_shift;
        let q_off = hit.q_off - query_shift;
        if q_start < 0 || q_off < 0 {
            continue;
        }
        let hsp = crate::hspstream::Hsp {
            score: hit.score,
            num_ident: 0,
            bit_score: 0.0,
            evalue: f64::MAX,
            query_offset: q_start,
            query_end: q_start + hit.len,
            query_gapped_start: q_off,
            subject_offset: hit.s_start,
            subject_end: hit.s_start + hit.len,
            subject_gapped_start: hit.s_off,
            context: 0,
            query_frame: 0,
            subject_frame: 0,
            num_gaps: 0,
            comp_adjustment_method: 0,
            edit_script: None,
            pat_info: None,
            map_info: None,
        };
        let _ = crate::hspstream::blast_hsp_list_save_hsp(&mut hsp_list, hsp);
    }
    if hsp_list.hsps.is_empty() {
        None
    } else {
        Some(hsp_list)
    }
}

/// blast-rs: assemble an owned RPS profile database from validated traceback
/// metadata plus consensus/profile sequences; not a direct NCBI C port.
///
/// This is the native boundary below a future CDD/RPS file reader. It validates
/// the profile offsets/PSSM rows through the traceback preparation path, builds
/// the RPS lookup backbone from consensus profile words, and constructs the
/// owned `BlastRpsLookupTable` used by the represented RPS word finder.
pub fn blast_rps_profile_database_new(
    traceback_info: crate::hspstream::RpsTracebackInfo,
    consensus_sequences: Vec<Vec<u8>>,
    wordsize: i32,
) -> Result<OwnedRpsProfileDatabase, i16> {
    if wordsize <= 0 {
        return Err(-1);
    }
    let gap_data = crate::hspstream::s_rps_gap_align_data_prepare(Some(&traceback_info))?;
    let alphabet_size = gap_data.alphabet_size as i32;
    validate_rps_consensus_sequences(&traceback_info, &consensus_sequences, alphabet_size)?;
    let charsize = crate::util::ilog2(alphabet_size as i64) + 1;
    let backbone_bits = charsize.checked_mul(wordsize).ok_or(-1i16)?;
    if !(1..usize::BITS as i32).contains(&backbone_bits) {
        return Err(-1);
    }
    let backbone_size = 1usize << backbone_bits as usize;
    let mut rps_backbone = vec![RpsBackboneCell::default(); backbone_size];

    for (profile_index, consensus) in consensus_sequences.iter().enumerate() {
        let profile_start = traceback_info.profile_header.start_offsets[profile_index];
        let profile_end = traceback_info.profile_header.start_offsets[profile_index + 1];
        if profile_start < 0 || profile_end < profile_start {
            return Err(-1);
        }
        if consensus.len() != (profile_end - profile_start) as usize {
            return Err(-1);
        }
        if consensus
            .iter()
            .any(|&residue| residue as i32 >= alphabet_size)
        {
            return Err(-1);
        }
        let wordsize = wordsize as usize;
        if consensus.len() < wordsize {
            continue;
        }
        for local_offset in 0..=consensus.len() - wordsize {
            let index =
                s_compute_table_index(wordsize as i32, charsize, &consensus[local_offset..]);
            let Some(cell) = rps_backbone.get_mut(index) else {
                return Err(-1);
            };
            cell.offset_pairs.push(OffsetPair {
                query_offset: profile_start + local_offset as i32 + wordsize as i32 - 1,
                subject_offset: 0,
            });
            cell.num_used = cell.offset_pairs.len() as i32;
        }
    }

    let lookup_info = BlastRpsInfo {
        alphabet_size,
        wordsize,
        rps_backbone,
        rps_pssm: gap_data.rps_pssm,
    };
    let mut lookup = None;
    if rps_lookup_table_new(Some(&lookup_info), &mut lookup) != 0 {
        return Err(-1);
    }
    Ok(OwnedRpsProfileDatabase {
        traceback_info,
        consensus_sequences,
        lookup_info,
        lookup: lookup.ok_or(-1i16)?,
    })
}

fn validate_rps_consensus_sequences(
    traceback_info: &crate::hspstream::RpsTracebackInfo,
    consensus_sequences: &[Vec<u8>],
    alphabet_size: i32,
) -> Result<(), i16> {
    let num_profiles = traceback_info.profile_header.num_profiles.max(0) as usize;
    if consensus_sequences.len() != num_profiles {
        return Err(-1);
    }
    for (profile_index, consensus) in consensus_sequences.iter().enumerate() {
        let profile_start = traceback_info.profile_header.start_offsets[profile_index];
        let profile_end = traceback_info.profile_header.start_offsets[profile_index + 1];
        if profile_start < 0 || profile_end < profile_start {
            return Err(-1);
        }
        if consensus.len() != (profile_end - profile_start) as usize {
            return Err(-1);
        }
        if consensus
            .iter()
            .any(|&residue| residue as i32 >= alphabet_size)
        {
            return Err(-1);
        }
    }
    Ok(())
}

fn blast_rps_profile_database_from_lookup_info(
    traceback_info: crate::hspstream::RpsTracebackInfo,
    consensus_sequences: Vec<Vec<u8>>,
    lookup_info: BlastRpsInfo,
) -> Result<OwnedRpsProfileDatabase, i16> {
    let mut lookup = None;
    if rps_lookup_table_new(Some(&lookup_info), &mut lookup) != 0 {
        return Err(-1);
    }
    Ok(OwnedRpsProfileDatabase {
        traceback_info,
        consensus_sequences,
        lookup_info,
        lookup: lookup.ok_or(-1i16)?,
    })
}

/// blast-rs: load an owned RPS profile database from native `.rps`/`.freq`
/// payloads plus explicit consensus sequences; not a direct NCBI C port.
///
/// This is the first file-backed boundary above `blast_rps_profile_database_new`:
/// it parses native-endian RPS profile metadata, validates optional frequency
/// ratios through the existing traceback preparation path, and then builds the
/// represented RPS lookup table from the supplied consensus profile sequences.
pub fn blast_rps_profile_database_from_native_bundle(
    bundle: NativeRpsProfileBundle,
) -> Result<OwnedRpsProfileDatabase, i16> {
    let profile_header =
        crate::hspstream::rps_profile_header_from_native_bytes(&bundle.profile_bytes)?;
    let num_profiles = profile_header.num_profiles;
    if num_profiles < 0 {
        return Err(-1);
    }
    let freq_ratios_header = bundle
        .freq_ratios_bytes
        .as_deref()
        .map(crate::hspstream::rps_freq_ratios_header_from_native_bytes)
        .transpose()?;
    let karlin_k = if let Some(aux_bytes) = bundle.aux_bytes.as_deref() {
        let aux_info = crate::hspstream::rps_aux_info_from_bytes(aux_bytes)?;
        if aux_info.karlin_k.len() < num_profiles as usize {
            return Err(-1);
        }
        if aux_info.profile_lengths.len() >= num_profiles as usize {
            for profile_index in 0..num_profiles as usize {
                let start = profile_header.start_offsets[profile_index];
                let end = profile_header.start_offsets[profile_index + 1];
                if start < 0
                    || end < start
                    || aux_info.profile_lengths[profile_index] != end - start
                {
                    return Err(-1);
                }
            }
        }
        aux_info.karlin_k
    } else {
        bundle.karlin_k
    };
    let traceback_info = crate::hspstream::RpsTracebackInfo {
        profile_header,
        freq_ratios_header,
        karlin_k,
    };
    let Some(lookup_bytes) = bundle.lookup_bytes else {
        return blast_rps_profile_database_new(
            traceback_info,
            bundle.consensus_sequences,
            bundle.wordsize,
        );
    };

    let gap_data = crate::hspstream::s_rps_gap_align_data_prepare(Some(&traceback_info))?;
    validate_rps_consensus_sequences(
        &traceback_info,
        &bundle.consensus_sequences,
        gap_data.alphabet_size as i32,
    )?;
    let lookup_info = rps_lookup_info_from_native_lookup_bytes(&lookup_bytes, gap_data.rps_pssm)?;
    if lookup_info.wordsize != bundle.wordsize {
        return Err(-1);
    }
    blast_rps_profile_database_from_lookup_info(
        traceback_info,
        bundle.consensus_sequences,
        lookup_info,
    )
}

fn rps_sidecar_path(base_path: &std::path::Path, extension: &str) -> std::path::PathBuf {
    base_path.with_extension(extension)
}

fn read_optional_rps_sidecar(
    base_path: &std::path::Path,
    extension: &str,
) -> Result<Option<Vec<u8>>, i16> {
    let path = rps_sidecar_path(base_path, extension);
    match std::fs::read(&path) {
        Ok(bytes) => Ok(Some(bytes)),
        Err(err) if err.kind() == std::io::ErrorKind::NotFound => Ok(None),
        Err(_) => Err(-1),
    }
}

/// blast-rs: assemble the native RPS profile bundle from database files;
/// not a direct NCBI C port.
///
/// NCBI resolves the RPS database through `CSeqDB::FindVolumePaths`, then maps
/// sidecar files with the same basename (`.rps`, optional `.freq`, optional
/// `.loo`, optional `.aux`) and fetches consensus profile sequences through
/// the protein BLAST database. Rust keeps that boundary reusable: this helper
/// reads the file payloads and materializes the consensus sequences, while
/// existing native parsers perform the format validation.
pub fn native_rps_profile_bundle_from_files<P: AsRef<std::path::Path>>(
    db_base_path: P,
    wordsize: i32,
    karlin_k: Vec<f64>,
) -> Result<NativeRpsProfileBundle, i16> {
    let db_base_path = db_base_path.as_ref();
    let profile_bytes = std::fs::read(rps_sidecar_path(db_base_path, "rps")).map_err(|_| -1i16)?;
    let freq_ratios_bytes = read_optional_rps_sidecar(db_base_path, "freq")?;
    let lookup_bytes = read_optional_rps_sidecar(db_base_path, "loo")?;
    let aux_bytes = read_optional_rps_sidecar(db_base_path, "aux")?;

    let consensus_db = crate::db::BlastDb::open(db_base_path).map_err(|_| -1i16)?;
    if consensus_db.db_type != crate::db::DbType::Protein {
        return Err(-1);
    }
    let consensus_sequences = (0..consensus_db.num_oids)
        .map(|oid| consensus_db.get_sequence(oid).to_vec())
        .collect();

    Ok(NativeRpsProfileBundle {
        profile_bytes,
        freq_ratios_bytes,
        lookup_bytes,
        aux_bytes,
        consensus_sequences,
        karlin_k,
        wordsize,
    })
}

/// blast-rs: convenience constructor for file-backed native RPS profile DBs;
/// not a direct NCBI C port.
pub fn blast_rps_profile_database_from_native_files<P: AsRef<std::path::Path>>(
    db_base_path: P,
    wordsize: i32,
    karlin_k: Vec<f64>,
) -> Result<OwnedRpsProfileDatabase, i16> {
    let bundle = native_rps_profile_bundle_from_files(db_base_path, wordsize, karlin_k)?;
    blast_rps_profile_database_from_native_bundle(bundle)
}

/// Rust ownership equivalent of NCBI `RPSLookupTableDestruct`
/// (`blast_aalookup.c:212`).
pub fn rps_lookup_table_destruct(
    lookup: &mut Option<BlastRpsLookupTable>,
) -> Option<BlastRpsLookupTable> {
    *lookup = None;
    None
}

/// Add one query/subject offset pair to an RPS bucket, matching
/// `s_AddToRPSBucket`'s grow-on-full behavior.
pub fn s_add_to_rps_bucket(bucket: &mut RpsBucket, q_off: u32, s_off: u32) {
    if bucket.num_alloc <= 0 {
        bucket.num_alloc = 1;
    }
    if bucket.num_used == bucket.num_alloc {
        bucket.num_alloc *= 2;
    }
    bucket.offset_pairs.push(OffsetPair {
        query_offset: q_off as i32,
        subject_offset: s_off as i32,
    });
    bucket.num_used += 1;
}

/// Scan a subject sequence against an RPS lookup table, matching
/// `BlastRPSScanSubject`'s bucketed offset collection.
///
/// Audit caveat: C returns hits through dynamically grown `RPSBucket` arrays
/// attached to a memory-mapped lookup table and can pause after four million
/// hits. Rust stores owned buckets and returns the count directly; the scan and
/// first-letter offset correction are kept in the same shape.
pub fn blast_rps_scan_subject(
    lookup: &mut BlastRpsLookupTable,
    sequence: &[u8],
    offset: &mut i32,
) -> i32 {
    for bucket in &mut lookup.bucket_array {
        bucket.num_used = 0;
        bucket.offset_pairs.clear();
    }

    if lookup.wordsize <= 0 || lookup.charsize < 0 || sequence.len() < lookup.wordsize as usize {
        *offset = sequence.len() as i32;
        return 0;
    }

    let start = (*offset).max(0) as usize;
    let wordsize = lookup.wordsize as usize;
    if start + wordsize > sequence.len() {
        *offset = sequence.len() as i32;
        return 0;
    }

    let table_correction = lookup.wordsize - 1;
    let last = sequence.len() - wordsize;
    let mut totalhits = 0;

    for subject_offset in start..=last {
        let index = s_compute_table_index(
            lookup.wordsize,
            lookup.charsize,
            &sequence[subject_offset..subject_offset + wordsize],
        );
        let Some(pv_word) = lookup.pv.get(index >> crate::stat::PV_ARRAY_BTS) else {
            continue;
        };
        if (pv_word & (1u32 << (index & crate::stat::PV_ARRAY_MASK as usize))) == 0 {
            continue;
        }

        let Some(cell) = lookup.rps_backbone.get(index) else {
            continue;
        };
        let numhits = cell.num_used.max(cell.offset_pairs.len() as i32);
        if numhits <= 0 {
            continue;
        }

        for pair in cell.offset_pairs.iter().take(numhits as usize) {
            let q_off = pair.query_offset - table_correction;
            if q_off < 0 {
                continue;
            }
            let bucket_index = (q_off / RPS_BUCKET_SIZE) as usize;
            if bucket_index >= lookup.bucket_array.len() {
                lookup
                    .bucket_array
                    .resize_with(bucket_index + 1, RpsBucket::default);
                lookup.num_buckets = lookup.bucket_array.len() as i32;
            }
            s_add_to_rps_bucket(
                &mut lookup.bucket_array[bucket_index],
                q_off as u32,
                subject_offset as u32,
            );
            totalhits += 1;
        }
    }

    *offset = (last + 1) as i32;
    totalhits
}

/// Name-matched scan wrapper for NCBI `BlastRPSWordFinder`.
///
/// Audit caveat: this is the pre-extension word-finding surface only. The C
/// function also dispatches one-hit/two-hit ungapped extension; those extension
/// states live outside the Rust RPS lookup module today. The returned
/// preliminary offset pairs are sorted by subject offset and query offset to
/// mirror the C word finder's sorted initial-hit-list boundary.
pub fn blast_rps_word_finder(lookup: &mut BlastRpsLookupTable, subject: &[u8]) -> Vec<OffsetPair> {
    let mut offset = 0;
    let _ = blast_rps_scan_subject(lookup, subject, &mut offset);
    let mut hits: Vec<OffsetPair> = lookup
        .bucket_array
        .iter()
        .flat_map(|bucket| bucket.offset_pairs.iter().copied())
        .collect();
    hits.sort_by_key(|hit| (hit.subject_offset, hit.query_offset));
    hits
}

/// NCBI: s_BlastRPSWordFinder_TwoHit (aa_ungapped.c).
///
/// Conservative C-shaped wrapper over the represented RPS lookup scan. NCBI's
/// static helper uses RPS lookup hits to drive ungapped extension and init-hit
/// list side effects; those RPS extension structures are not represented in
/// this module today. This wrapper preserves the two-hit gate that decides
/// which same-diagonal RPS word hits are eligible for that extension: the
/// second hit must be within `window_size` subject residues and must not
/// overlap the first word.
pub fn s_blast_rps_word_finder_two_hit(
    lookup: &mut BlastRpsLookupTable,
    subject: &[u8],
    window_size: i32,
) -> Vec<RpsTwoHitSeed> {
    s_blast_rps_word_finder_two_hit_scan(lookup, subject, window_size).seeds
}

fn s_blast_rps_word_finder_two_hit_scan(
    lookup: &mut BlastRpsLookupTable,
    subject: &[u8],
    window_size: i32,
) -> RpsWordFinderScan {
    s_blast_rps_word_finder_two_hit_scan_impl(lookup, subject, window_size, None)
}

#[derive(Debug, Clone, Copy)]
struct RpsTwoHitPayloadConfig {
    x_dropoff: i32,
    score_cutoff: i32,
}

fn s_blast_rps_word_finder_two_hit_scan_impl(
    lookup: &mut BlastRpsLookupTable,
    subject: &[u8],
    window_size: i32,
    payload_config: Option<RpsTwoHitPayloadConfig>,
) -> RpsWordFinderScan {
    if lookup.wordsize <= 0 || window_size <= lookup.wordsize {
        return RpsWordFinderScan::default();
    }

    let mut hits = blast_rps_word_finder(lookup, subject);
    let total_hits = hits.len() as i32;
    hits.retain(|hit| hit.query_offset >= 0 && hit.subject_offset >= 0);
    if hits.len() < 2 {
        return RpsWordFinderScan {
            total_hits,
            ..RpsWordFinderScan::default()
        };
    }

    let mut seeds = Vec::new();
    let mut init_hits = Vec::new();
    let mut diag_states: std::collections::HashMap<i32, RpsTwoHitDiagState> =
        std::collections::HashMap::new();
    let ws = lookup.wordsize;
    let mut hits_extended = 0i32;

    for second in hits {
        let diag = second.query_offset - second.subject_offset;
        let state = diag_states.entry(diag).or_insert(RpsTwoHitDiagState {
            last_hit: -window_size,
            flag: false,
        });
        let s_off = second.subject_offset;

        if state.flag {
            if s_off < state.last_hit {
                continue;
            }
            state.last_hit = s_off;
            state.flag = false;
            continue;
        }

        let diff = s_off - state.last_hit;
        if diff >= window_size {
            state.last_hit = s_off;
            state.flag = false;
            continue;
        }
        if diff < ws {
            continue;
        }

        let first = OffsetPair {
            query_offset: state.last_hit + diag,
            subject_offset: state.last_hit,
        };
        let seed = RpsTwoHitSeed { first, second };
        seeds.push(seed);
        hits_extended += 1;

        if let Some(config) = payload_config {
            match rps_extend_two_hit_payload(
                &lookup.rps_pssm,
                subject,
                ws as usize,
                config.x_dropoff,
                config.score_cutoff,
                seed,
            ) {
                RpsTwoHitOutcome::Reached {
                    init_hit,
                    s_last_off,
                } => {
                    if let Some(init_hit) = init_hit {
                        init_hits.push(init_hit);
                    }
                    state.last_hit = s_last_off - (ws - 1);
                    state.flag = true;
                }
                RpsTwoHitOutcome::NoReach { init_hit } => {
                    if let Some(init_hit) = init_hit {
                        init_hits.push(init_hit);
                    }
                    state.last_hit = s_off;
                    state.flag = false;
                }
            }
        } else {
            state.last_hit = s_off;
            state.flag = false;
        }
    }

    RpsWordFinderScan {
        seeds,
        total_hits,
        hits_extended,
        init_hits,
    }
}

fn rps_pssm_score(pssm: &[Vec<i32>], q_off: usize, residue: u8) -> Option<i32> {
    pssm.get(q_off)
        .and_then(|row| row.get(residue as usize))
        .copied()
}

fn rps_extend_left(
    pssm: &[Vec<i32>],
    subject: &[u8],
    q_start: usize,
    s_start: usize,
    x_dropoff: i32,
) -> (i32, i32) {
    let mut score = 0i32;
    let mut best = 0i32;
    let mut best_d = 0i32;
    if q_start == 0 || s_start == 0 {
        return (0, 0);
    }

    let mut qi = q_start - 1;
    let mut si = s_start - 1;
    let mut d = 1i32;
    loop {
        let Some(cell_score) = rps_pssm_score(pssm, qi, subject[si]) else {
            break;
        };
        score += cell_score;
        if score > best {
            best = score;
            best_d = d;
        }
        if best - score >= x_dropoff {
            break;
        }
        if qi == 0 || si == 0 {
            break;
        }
        qi -= 1;
        si -= 1;
        d += 1;
    }
    (best, best_d)
}

fn rps_extend_right(
    pssm: &[Vec<i32>],
    subject: &[u8],
    q_start: usize,
    s_start: usize,
    x_dropoff: i32,
    init_score: i32,
) -> (i32, i32, i32) {
    let mut score = init_score;
    let mut best = init_score;
    let mut best_d = 0i32;
    let mut qi = q_start;
    let mut si = s_start;
    let mut last_off_delta = 0i32;

    while qi < pssm.len() && si < subject.len() {
        let Some(cell_score) = rps_pssm_score(pssm, qi, subject[si]) else {
            break;
        };
        score += cell_score;
        if score > best {
            best = score;
            best_d = (qi - q_start + 1) as i32;
        }
        if score <= 0 || best - score >= x_dropoff {
            last_off_delta = (qi - q_start) as i32;
            break;
        }
        qi += 1;
        si += 1;
        last_off_delta = (qi - q_start) as i32;
    }

    (best, best_d, last_off_delta)
}

fn rps_extend_two_hit_payload(
    pssm: &[Vec<i32>],
    subject: &[u8],
    word_size: usize,
    x_dropoff: i32,
    score_cutoff: i32,
    seed: RpsTwoHitSeed,
) -> RpsTwoHitOutcome {
    let q_right_off = seed.second.query_offset as usize;
    let s_right_off = seed.second.subject_offset as usize;
    let s_left_off = seed.first.subject_offset as usize + word_size;

    let mut wscore = 0i32;
    let mut best_wscore = 0i32;
    let mut right_d = 0usize;
    for k in 0..word_size {
        let qi = q_right_off + k;
        let si = s_right_off + k;
        if si >= subject.len() {
            break;
        }
        let Some(cell_score) = rps_pssm_score(pssm, qi, subject[si]) else {
            break;
        };
        wscore += cell_score;
        if wscore > best_wscore {
            best_wscore = wscore;
            right_d = k + 1;
        }
    }

    let ext_q = q_right_off + right_d;
    let ext_s = s_right_off + right_d;
    let (left_score, left_d) = rps_extend_left(pssm, subject, ext_q, ext_s, x_dropoff);
    let reached_first = left_d >= (ext_s as i32 - s_left_off as i32);
    if !reached_first {
        let init_hit = (left_score >= score_cutoff && left_d > 0).then(|| {
            let q_start = ext_q as i32 - left_d;
            let s_start = ext_s as i32 - left_d;
            RpsInitHitRecord {
                q_start,
                s_start,
                q_off: seed.second.query_offset,
                s_off: seed.second.subject_offset,
                len: left_d,
                score: left_score,
            }
        });
        return RpsTwoHitOutcome::NoReach { init_hit };
    }

    let (right_score, right_d, s_last_off_delta) =
        rps_extend_right(pssm, subject, ext_q, ext_s, x_dropoff, left_score);
    let score = left_score.max(right_score);
    let s_last_off = ext_s as i32 + s_last_off_delta;
    let q_start = ext_q as i32 - left_d;
    let s_start = ext_s as i32 - left_d;
    let init_hit = (score >= score_cutoff).then_some(RpsInitHitRecord {
        q_start,
        s_start,
        q_off: seed.second.query_offset,
        s_off: seed.second.subject_offset,
        len: left_d + right_d,
        score,
    });
    RpsTwoHitOutcome::Reached {
        init_hit,
        s_last_off,
    }
}

/// blast-rs: side-channel preserving wrapper around
/// `s_BlastRPSWordFinder_TwoHit` (`aa_ungapped.c`).
///
/// The local RPS representation still stops at the two-hit seed gate, before
/// C's RPS ungapped extension. This wrapper nevertheless preserves the
/// C-shaped `BlastInitHitList`/`BlastUngappedStats` mutation boundary for the
/// accepted RPS seed pairs so callers can audit side-channel plumbing without
/// depending on the future extension implementation.
pub fn s_blast_rps_word_finder_two_hit_with_side_channels(
    lookup: &mut BlastRpsLookupTable,
    subject: &[u8],
    window_size: i32,
    init_hitlist: Option<&mut crate::extend::InitHitList>,
    ungapped_stats: Option<&mut crate::diagnostics::UngappedStats>,
) -> RpsWordFinderScan {
    let scan = s_blast_rps_word_finder_two_hit_scan(lookup, subject, window_size);
    let saved_hits = scan.seeds.len() as i32;
    if let Some(init_hitlist) = init_hitlist {
        for seed in &scan.seeds {
            crate::extend::blast_save_initial_hit(
                init_hitlist,
                seed.second.query_offset,
                seed.second.subject_offset,
                None,
            );
        }
    }
    crate::diagnostics::blast_ungapped_stats_update(
        ungapped_stats,
        scan.total_hits,
        scan.hits_extended,
        saved_hits,
    );
    scan
}

/// blast-rs: RPS two-hit side-channel wrapper with represented ungapped
/// extension payloads.
///
/// This keeps the existing conservative RPS seed gate, then uses the owned
/// RPS PSSM rows to fill `BlastInitHsp::ungapped_data` through
/// `BlastSaveInitHsp`, matching the payload boundary C reaches after
/// `s_BlastRPSWordFinder_TwoHit` drives ungapped extension.
pub fn s_blast_rps_word_finder_two_hit_with_extension_payloads(
    lookup: &mut BlastRpsLookupTable,
    subject: &[u8],
    window_size: i32,
    x_dropoff: i32,
    score_cutoff: i32,
    init_hitlist: Option<&mut crate::extend::InitHitList>,
    ungapped_stats: Option<&mut crate::diagnostics::UngappedStats>,
) -> RpsWordFinderScan {
    let scan = if lookup.wordsize > 0 && !lookup.rps_pssm.is_empty() {
        s_blast_rps_word_finder_two_hit_scan_impl(
            lookup,
            subject,
            window_size,
            Some(RpsTwoHitPayloadConfig {
                x_dropoff,
                score_cutoff,
            }),
        )
    } else {
        s_blast_rps_word_finder_two_hit_scan(lookup, subject, window_size)
    };

    let saved_hits = scan.init_hits.len() as i32;
    if let Some(init_hitlist) = init_hitlist {
        for hit in &scan.init_hits {
            crate::extend::blast_save_init_hsp(
                init_hitlist,
                hit.q_start,
                hit.s_start,
                hit.q_off,
                hit.s_off,
                hit.len,
                hit.score,
            );
        }
        crate::extend::blast_init_hit_list_sort_by_score(init_hitlist);
    }
    crate::diagnostics::blast_ungapped_stats_update(
        ungapped_stats,
        scan.total_hits,
        scan.hits_extended,
        saved_hits,
    );
    scan
}

pub const OFFSET_ARRAY_SIZE: i32 = 4096;

/// 1-1 ownership translation of `LookupTableWrapFree`.
pub fn lookup_table_wrap_free(lookup: &mut Option<LookupTableWrap>) -> Option<LookupTableWrap> {
    let Some(table) = lookup.take() else {
        return None;
    };

    match table {
        LookupTableWrap::Na(table) => {
            let mut inner = Some(table);
            let _ = blast_na_lookup_table_destruct(&mut inner);
        }
        LookupTableWrap::NaHash(table) => {
            let mut inner = Some(table);
            let _ = blast_na_hash_lookup_table_destruct(&mut inner);
        }
        LookupTableWrap::SmallNa(table) => {
            let mut inner = Some(table);
            let _ = blast_small_na_lookup_table_destruct(&mut inner);
        }
        LookupTableWrap::Megablast(table) => {
            let mut inner = Some(table);
            let _ = blast_mb_lookup_table_destruct(&mut inner);
        }
        LookupTableWrap::Aa(_table) => {}
        LookupTableWrap::Rps(table) => {
            let mut inner = Some(table);
            let _ = rps_lookup_table_destruct(&mut inner);
        }
    }

    *lookup = None;
    None
}

/// Rust ownership-shaped equivalent of NCBI `LookupTableWrapInit`
/// (`lookup_wrap.c:47`) once the concrete lookup table has been built.
pub fn lookup_table_wrap_init(
    lookup_wrap_ptr: &mut Option<LookupTableWrap>,
    lookup: LookupTableWrap,
) -> i16 {
    let _ = lookup_table_wrap_free(lookup_wrap_ptr);
    *lookup_wrap_ptr = Some(lookup);
    0
}

/// blast-rs: Nucleotide construction arm for represented Rust lookup variants;
/// not a direct NCBI C port.
pub fn lookup_table_wrap_init_nucleotide(
    query_sequence: &[u8],
    locations: &[crate::util::SSeqRange],
    lookup_options: &crate::options::LookupTableOptions,
    lookup_wrap_ptr: &mut Option<LookupTableWrap>,
    seq_src: Option<&mut dyn crate::seqsrc::BlastSeqSource>,
    num_threads: u32,
) -> i16 {
    match lookup_options.lut_type {
        LookupTableType::SmallNaLookup
        | LookupTableType::NaLookup
        | LookupTableType::MegablastLookup
        | LookupTableType::NaHashLookup => {}
        _ => {
            let _ = lookup_table_wrap_free(lookup_wrap_ptr);
            return -1;
        }
    }

    let _ = lookup_table_wrap_free(lookup_wrap_ptr);

    let mut max_q_off = 0;
    let num_table_entries = estimate_num_table_entries(locations, &mut max_q_off);
    let mut lut_width = 0;
    let Some(chosen_type) = blast_choose_na_lookup_table(
        Some(lookup_options),
        num_table_entries,
        max_q_off,
        &mut lut_width,
    ) else {
        return -1;
    };

    match chosen_type {
        LookupTableType::MegablastLookup => {
            let mut lookup = None;
            let status = blast_mb_lookup_table_new(
                query_sequence,
                locations,
                &mut lookup,
                lookup_options,
                num_table_entries,
                lut_width,
                seq_src,
            );
            if status != 0 {
                return status;
            }
            let Some(lookup) = lookup else {
                return -1;
            };
            *lookup_wrap_ptr = Some(LookupTableWrap::Megablast(lookup));
            0
        }
        LookupTableType::SmallNaLookup => {
            let mut small = None;
            let status = blast_small_na_lookup_table_new(
                query_sequence,
                locations,
                &mut small,
                lookup_options,
                lut_width,
            );
            if status == 0 {
                let Some(lookup) = small else {
                    return -1;
                };
                *lookup_wrap_ptr = Some(LookupTableWrap::SmallNa(lookup));
                return 0;
            }

            let mut standard = None;
            let fallback_status = blast_na_lookup_table_new(
                query_sequence,
                locations,
                &mut standard,
                lookup_options,
                lut_width,
            );
            if fallback_status != 0 {
                return fallback_status;
            }
            let Some(lookup) = standard else {
                return -1;
            };
            *lookup_wrap_ptr = Some(LookupTableWrap::Na(lookup));
            0
        }
        LookupTableType::NaHashLookup => {
            let mut lookup = None;
            let status = blast_na_hash_lookup_table_new(
                query_sequence,
                locations,
                &mut lookup,
                lookup_options,
                seq_src,
                num_threads,
            );
            if status != 0 {
                return status;
            }
            let Some(lookup) = lookup else {
                return -1;
            };
            *lookup_wrap_ptr = Some(LookupTableWrap::NaHash(lookup));
            0
        }
        LookupTableType::NaLookup => {
            let mut lookup = None;
            let status = blast_na_lookup_table_new(
                query_sequence,
                locations,
                &mut lookup,
                lookup_options,
                lut_width,
            );
            if status != 0 {
                return status;
            }
            let Some(lookup) = lookup else {
                return -1;
            };
            *lookup_wrap_ptr = Some(LookupTableWrap::Na(lookup));
            0
        }
        _ => -1,
    }
}

/// Rust ownership-shaped equivalent of NCBI `LookupTableWrapInit_MT`
/// (`lookup_wrap.c:61`) once the concrete lookup table has been built.
///
/// C constructs the table inside this function after switching on
/// `lookup_options->lut_type`. Rust also exposes
/// `lookup_table_wrap_init_nucleotide` for the represented nucleotide
/// construction branch; this lower-level helper preserves the ownership handoff
/// when callers already hold a constructed typed lookup table.
pub fn lookup_table_wrap_init_mt(
    lookup_wrap_ptr: &mut Option<LookupTableWrap>,
    lookup: LookupTableWrap,
    num_threads: u32,
) -> i16 {
    let _ = lookup_table_wrap_free(lookup_wrap_ptr);

    let initialized = match lookup {
        LookupTableWrap::Aa(table) => {
            // C eAaLookupTable: BlastAaLookupTableNew, index query, finalize.
            // Rust callers provide the finalized AaLookupTable value.
            let _bone_type_is_query_length_dependent = true;
            LookupTableWrap::Aa(table)
        }
        LookupTableWrap::Megablast(table) => {
            // C eMBLookupTable and chosen mixed/small-na fall-through paths:
            // allocation is completed by the nucleotide lookup builder before
            // ownership reaches this wrapper.
            LookupTableWrap::Megablast(table)
        }
        LookupTableWrap::Na(table) => {
            // C eNaLookupTable: standard nucleotide lookup ownership.
            LookupTableWrap::Na(table)
        }
        LookupTableWrap::NaHash(table) => {
            // C eNaHashLookupTable: construction receives num_threads before
            // ownership reaches the wrapper; scan callbacks are represented by
            // the typed NaHash scanner helpers.
            LookupTableWrap::NaHash(table)
        }
        LookupTableWrap::SmallNa(table) => {
            // C eSmallNaLookupTable can fall back to eNaLookupTable on
            // allocation failure. Rust represents the successfully constructed
            // small-table path here; standard NaLookup is held by scanner-side
            // structures outside LookupTableWrap.
            LookupTableWrap::SmallNa(table)
        }
        LookupTableWrap::Rps(table) => {
            // C eRPSLookupTable copies RPSInfo fields into BlastRPSLookupTable.
            // Rust callers provide the owned profile/backbone view directly.
            LookupTableWrap::Rps(table)
        }
    };

    // Keep num_threads observable at this boundary for audit: if/when the
    // represented NaHash lookup is added, this is the construction parameter
    // that must be threaded into that variant.
    let _nahash_num_threads = num_threads;
    *lookup_wrap_ptr = Some(initialized);
    0
}

/// 1-1 translation of `GetOffsetArraySize` for represented lookup variants.
pub fn get_offset_array_size(lookup: &LookupTableWrap) -> i32 {
    let longest_chain = match lookup {
        LookupTableWrap::Na(table) => table.longest_chain,
        LookupTableWrap::NaHash(table) => table.longest_chain,
        LookupTableWrap::SmallNa(table) => table.longest_chain,
        LookupTableWrap::Megablast(table) => table.longest_chain,
        LookupTableWrap::Aa(table) => table
            .backbone
            .iter()
            .map(|chain| chain.len() as i32)
            .max()
            .unwrap_or(0),
        LookupTableWrap::Rps(table) => table
            .bucket_array
            .iter()
            .map(|bucket| bucket.num_used)
            .max()
            .unwrap_or(0),
    };
    OFFSET_ARRAY_SIZE + longest_chain
}

/// 1-1 ownership translation of `BackboneCellFree`.
pub fn backbone_cell_free(cell: &mut Option<Box<BackboneCell>>) -> Option<Box<BackboneCell>> {
    let mut current = cell.take();
    while let Some(mut node) = current {
        current = node.next.take();
    }
    None
}

/// 1-1 translation of `BackboneCellNew`.
pub fn backbone_cell_new(word: u32, offset: i32) -> Option<Box<BackboneCell>> {
    let mut cell = Box::new(BackboneCell::default());
    if backbone_cell_init(Some(cell.as_mut()), word, offset) != 0 {
        let mut owned = Some(cell);
        return backbone_cell_free(&mut owned);
    }
    Some(cell)
}

/// 1-1 translation of `BackboneCellInit`.
pub fn backbone_cell_init(cell: Option<&mut BackboneCell>, word: u32, offset: i32) -> i32 {
    let Some(cell) = cell else {
        return -1;
    };
    cell.word = word;
    cell.offset = offset;
    cell.num_offsets = 1;
    0
}

/// 1-1 translation of `s_AddWordHit` over Rust-owned backbone cells.
pub fn s_add_word_hit(
    backbone: &mut [BackboneCell],
    offsets: &mut [i32],
    wordsize: i32,
    charsize: i32,
    seq: &[u8],
    mut offset: i32,
    hash_func: fn(&[u8], u32) -> u32,
    mask: u32,
    pv_array: Option<&[u32]>,
) -> i16 {
    let mut large_index = 0u32;
    for &base in seq.iter().take(wordsize as usize) {
        large_index = (large_index << charsize) | base as u32;
    }

    if let Some(pv_array) = pv_array {
        let word = (large_index >> crate::stat::PV_ARRAY_BTS) as usize;
        let bit = large_index & crate::stat::PV_ARRAY_MASK;
        if word >= pv_array.len() || (pv_array[word] & (1u32 << bit)) == 0 {
            return 0;
        }
    }

    let index = hash_func(&large_index.to_ne_bytes(), mask) as usize;
    if index >= backbone.len() {
        return -1;
    }
    offset += 1;
    if offset < 0 || offset as usize >= offsets.len() {
        return -1;
    }

    if backbone[index].num_offsets == 0 {
        return backbone_cell_init(Some(&mut backbone[index]), large_index, offset) as i16;
    }

    let mut current = &mut backbone[index];
    while current.next.is_some() && current.word != large_index {
        current = current.next.as_mut().expect("checked next");
    }

    if current.word == large_index {
        offsets[offset as usize] = current.offset;
        current.offset = offset;
        current.num_offsets += 1;
    } else {
        current.next = backbone_cell_new(large_index, offset);
        if current.next.is_none() {
            return -1;
        }
    }

    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_offset_pair() {
        let pair = OffsetPair {
            query_offset: 10,
            subject_offset: 20,
        };
        assert_eq!(pair.query_offset, 10);
        assert_eq!(pair.subject_offset, 20);
    }

    #[test]
    fn test_lookup_type() {
        let table = LookupTableWrap::Megablast(MbLookupTable {
            word_length: 28,
            lut_word_length: 12,
            discontiguous: false,
            template_length: 0,
            template_type: DiscTemplateType::Contiguous,
            two_templates: false,
            second_template_type: DiscTemplateType::Contiguous,
            hashtable: vec![],
            hashtable2: vec![],
            next_pos: vec![],
            next_pos2: vec![],
            pv_array: vec![],
            pv_array_bts: crate::stat::PV_ARRAY_BTS as i32,
            longest_chain: 0,
            scan_step: 1,
        });
        assert_eq!(table.table_type(), LookupTableType::MegablastLookup);
        assert_eq!(table.word_length(), 28);
    }

    #[test]
    fn translated_na_lookup_membership_helpers_match_c_shapes() {
        let mb = LookupTableWrap::Megablast(MbLookupTable {
            word_length: 12,
            lut_word_length: 12,
            discontiguous: false,
            template_length: 0,
            template_type: DiscTemplateType::Contiguous,
            two_templates: false,
            second_template_type: DiscTemplateType::Contiguous,
            hashtable: vec![0, 2, 0, 0],
            hashtable2: Vec::new(),
            next_pos: vec![0, 0, 0],
            next_pos2: Vec::new(),
            pv_array: vec![1 << 1],
            pv_array_bts: crate::stat::PV_ARRAY_BTS as i32,
            longest_chain: 1,
            scan_step: 1,
        });
        assert!(s_mb_lookup(Some(&mb), 1, 1));
        assert!(!s_mb_lookup(Some(&mb), 1, 0));
        assert_eq!(
            blast_choose_na_extend(Some(&mb)),
            Some(NaExtendChoice {
                lookup_callback: Some(NaLookupCallback::Megablast),
                extend_callback: Some(NaExtendCallback::Direct),
            })
        );
        assert_eq!(
            s_is_seed_masked(Some(&mb), &[1, 0, 0, 0], 0, 4, 1),
            Some(false)
        );
        assert_eq!(
            s_is_seed_masked(Some(&mb), &[1, 0, 0, 0], 0, 4, 0),
            Some(true)
        );

        let small = LookupTableWrap::SmallNa(SmallNaLookupTable {
            word_length: 4,
            backbone: vec![-2, 7],
            overflow: vec![0, 0, 3, 5, -1],
            pv_array: Vec::new(),
            longest_chain: 2,
            scan_step: 1,
        });
        assert!(s_small_na_lookup(Some(&small), 0, 3));
        assert!(s_small_na_lookup(Some(&small), 0, 5));
        assert!(s_small_na_lookup(Some(&small), 1, 7));
        assert!(!s_small_na_lookup(Some(&small), 0, 9));
        let LookupTableWrap::SmallNa(small_table) = &small else {
            unreachable!();
        };
        let mut pairs = Vec::new();
        assert_eq!(
            s_blast_small_na_lookup_retrieve(small_table, 0, &mut pairs, 11),
            2
        );
        assert_eq!(
            pairs,
            vec![
                OffsetPair {
                    query_offset: 3,
                    subject_offset: 11,
                },
                OffsetPair {
                    query_offset: 5,
                    subject_offset: 11,
                }
            ]
        );
        assert_eq!(
            blast_choose_na_extend(Some(&small)),
            Some(NaExtendChoice {
                lookup_callback: Some(NaLookupCallback::SmallNa),
                extend_callback: Some(NaExtendCallback::Direct),
            })
        );

        let mut na = BlastNaLookupTable {
            mask: 3,
            word_length: 4,
            lut_word_length: 4,
            thick_backbone: vec![NaLookupBackboneCell::default(); 4],
            pv: vec![0],
            ..BlastNaLookupTable::default()
        };
        na.thick_backbone[2].num_used = 2;
        na.thick_backbone[2].entries = [4, 8, 0];
        na.pv[0] = 1 << 2;
        assert!(s_na_lookup(Some(&na), 2, 8));
        assert!(!s_na_lookup(Some(&na), 2, 9));

        let na_wrap = LookupTableWrap::Na(na);
        assert_eq!(
            blast_choose_na_extend(Some(&na_wrap)),
            Some(NaExtendChoice {
                lookup_callback: Some(NaLookupCallback::Na),
                extend_callback: Some(NaExtendCallback::Direct),
            })
        );
        assert_eq!(
            s_is_seed_masked(Some(&na_wrap), &[2, 0, 0, 0], 0, 4, 8),
            Some(false)
        );

        let na_hash = LookupTableWrap::NaHash(BlastNaHashLookupTable::default());
        assert_eq!(
            blast_choose_na_extend(Some(&na_hash)),
            Some(NaExtendChoice {
                lookup_callback: None,
                extend_callback: None,
            })
        );
    }

    #[test]
    fn translated_blast_choose_na_extend_matches_c_extension_callback_rules() {
        let aligned = LookupTableWrap::Megablast(MbLookupTable {
            word_length: 28,
            lut_word_length: 12,
            discontiguous: false,
            template_length: 0,
            template_type: DiscTemplateType::Contiguous,
            two_templates: false,
            second_template_type: DiscTemplateType::Contiguous,
            hashtable: Vec::new(),
            hashtable2: Vec::new(),
            next_pos: Vec::new(),
            next_pos2: Vec::new(),
            pv_array: Vec::new(),
            pv_array_bts: crate::stat::PV_ARRAY_BTS as i32,
            longest_chain: 0,
            scan_step: 4,
        });
        assert_eq!(
            blast_choose_na_extend(Some(&aligned))
                .unwrap()
                .extend_callback,
            Some(NaExtendCallback::Aligned)
        );

        let generic = LookupTableWrap::Megablast(MbLookupTable {
            word_length: 28,
            lut_word_length: 12,
            discontiguous: false,
            template_length: 0,
            template_type: DiscTemplateType::Contiguous,
            two_templates: false,
            second_template_type: DiscTemplateType::Contiguous,
            hashtable: Vec::new(),
            hashtable2: Vec::new(),
            next_pos: Vec::new(),
            next_pos2: Vec::new(),
            pv_array: Vec::new(),
            pv_array_bts: crate::stat::PV_ARRAY_BTS as i32,
            longest_chain: 0,
            scan_step: 1,
        });
        assert_eq!(
            blast_choose_na_extend(Some(&generic))
                .unwrap()
                .extend_callback,
            Some(NaExtendCallback::Generic)
        );

        let small_generic = LookupTableWrap::SmallNa(SmallNaLookupTable {
            word_length: 9,
            backbone: vec![0; 1usize << (BITS_PER_NUC * 7)],
            overflow: Vec::new(),
            pv_array: Vec::new(),
            longest_chain: 0,
            scan_step: 3,
        });
        assert_eq!(
            small_na_lut_word_length(match &small_generic {
                LookupTableWrap::SmallNa(table) => table,
                _ => unreachable!(),
            }),
            7
        );
        assert_eq!(
            blast_choose_na_extend(Some(&small_generic))
                .unwrap()
                .extend_callback,
            Some(NaExtendCallback::SmallGeneric)
        );
    }

    fn all_a_mb_lookup(allowed_q_offsets: &[i32]) -> LookupTableWrap {
        let max_query_index = allowed_q_offsets
            .iter()
            .map(|q| q + 1)
            .max()
            .unwrap_or(0)
            .max(0) as usize;
        let mut hashtable = vec![0; 1];
        let mut next_pos = vec![0; max_query_index + 1];
        for &q_pos in allowed_q_offsets {
            let query_index = q_pos + 1;
            if query_index <= 0 {
                continue;
            }
            let query_index = query_index as usize;
            next_pos[query_index] = hashtable[0];
            hashtable[0] = query_index as i32;
        }

        LookupTableWrap::Megablast(MbLookupTable {
            word_length: 4,
            lut_word_length: 2,
            discontiguous: false,
            template_length: 0,
            template_type: DiscTemplateType::Contiguous,
            two_templates: false,
            second_template_type: DiscTemplateType::Contiguous,
            hashtable,
            hashtable2: Vec::new(),
            next_pos,
            next_pos2: Vec::new(),
            pv_array: vec![1],
            pv_array_bts: crate::stat::PV_ARRAY_BTS as i32,
            longest_chain: allowed_q_offsets.len() as i32,
            scan_step: 1,
        })
    }

    #[test]
    fn translated_type_of_word_matches_mask_and_double_word_paths() {
        let subject = crate::encoding::pack_ncbi2na_bases(&[0; 20]);
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[12]);
        let lookup = all_a_mb_lookup(&[1, 2, 3, 4, 6, 8, 10]);
        let mut q_off = 0;
        let mut s_off = 0;
        let mut extended = -1;

        assert_eq!(
            s_type_of_word(
                &subject,
                &mut q_off,
                &mut s_off,
                Some(&[crate::util::SSeqRange { left: 0, right: 0 }]),
                &query_info,
                20,
                4,
                2,
                Some(&lookup),
                false,
                &mut extended,
            ),
            1
        );
        assert_eq!((q_off, s_off, extended), (1, 1, 1));

        let lookup = all_a_mb_lookup(&[0, 2, 4, 6, 8, 10]);
        let mut q_off = 0;
        let mut s_off = 0;
        let mut extended = 0;
        assert_eq!(
            s_type_of_word(
                &subject,
                &mut q_off,
                &mut s_off,
                None,
                &query_info,
                20,
                4,
                2,
                Some(&lookup),
                true,
                &mut extended,
            ),
            2
        );
        assert_eq!((q_off, s_off, extended), (0, 0, 4));
    }

    #[test]
    fn translated_blast_choose_na_lookup_table_matches_c_thresholds() {
        let mut options = crate::options::LookupTableOptions {
            word_size: 9,
            program_number: crate::program::BLASTN,
            ..crate::options::LookupTableOptions::default()
        };
        let mut width = 0;

        assert_eq!(
            blast_choose_na_lookup_table(Some(&options), 1249, 100, &mut width),
            Some(LookupTableType::SmallNaLookup)
        );
        assert_eq!(width, 7);
        assert_eq!(
            blast_choose_na_lookup_table(Some(&options), 20_999, 100, &mut width),
            Some(LookupTableType::SmallNaLookup)
        );
        assert_eq!(width, 8);
        assert_eq!(
            blast_choose_na_lookup_table(Some(&options), 21_000, 100, &mut width),
            Some(LookupTableType::MegablastLookup)
        );
        assert_eq!(width, 9);

        options.word_size = 12;
        assert_eq!(
            blast_choose_na_lookup_table(Some(&options), 900_000, 100, &mut width),
            Some(LookupTableType::MegablastLookup)
        );
        assert_eq!(width, 12);

        options.word_size = 8;
        assert_eq!(
            blast_choose_na_lookup_table(Some(&options), 100, 32_768, &mut width),
            Some(LookupTableType::NaLookup)
        );
        assert_eq!(width, 7);

        options.mb_template_length = 18;
        options.word_size = 11;
        assert_eq!(
            blast_choose_na_lookup_table(Some(&options), 0, 0, &mut width),
            Some(LookupTableType::MegablastLookup)
        );
        assert_eq!(width, 11);

        options.mb_template_length = 0;
        options.program_number = crate::program::MAPPING;
        options.word_size = 16;
        options.db_filter = true;
        assert_eq!(
            blast_choose_na_lookup_table(Some(&options), 0, 0, &mut width),
            Some(LookupTableType::NaHashLookup)
        );
        assert_eq!(width, 16);
        assert!(blast_choose_na_lookup_table(None, 0, 0, &mut width).is_none());
    }

    #[test]
    fn estimate_num_table_entries_matches_c_right_minus_left_rule() {
        let ranges = [
            crate::util::SSeqRange { left: 0, right: 4 },
            crate::util::SSeqRange {
                left: 10,
                right: 15,
            },
        ];
        let mut max_off = 0;

        assert_eq!(estimate_num_table_entries(&ranges, &mut max_off), 9);
        assert_eq!(max_off, 15);
    }

    #[test]
    fn translated_mapper_word_hits_new_matches_c_layout_rules() {
        let query = crate::util::BlastSequenceBlk {
            length: 250,
            ..crate::util::BlastSequenceBlk::default()
        };
        let query_info = crate::queryinfo::QueryInfo::new_blastp(&vec![10; 250]);
        let word_hits = mapper_word_hits_new(&query, &query_info).expect("word hits");

        assert_eq!(word_hits.num_arrays, 2);
        assert_eq!(word_hits.array_size, 1000);
        assert_eq!(word_hits.pair_arrays.len(), 2);
        assert_eq!(word_hits.pair_arrays[0].len(), 1000);
        assert_eq!(word_hits.num, vec![0, 0]);
        assert_eq!(word_hits.divisor, 126);
        assert_eq!(word_hits.last_diag.len(), query_info.contexts.len());
        assert!(word_hits
            .last_pos
            .iter()
            .take(query_info.num_queries as usize)
            .all(|&pos| pos == i32::MIN));

        let mut slot = Some(word_hits);
        assert!(mapper_word_hits_free(&mut slot).is_none());
        assert!(slot.is_none());
    }

    #[test]
    fn test_small_na_lookup_table_word_length() {
        let table = SmallNaLookupTable {
            word_length: 8,
            backbone: vec![-1; 65536], // 4^8 = 65536 entries
            overflow: vec![],
            pv_array: vec![],
            longest_chain: 0,
            scan_step: 4,
        };
        assert_eq!(table.word_length, 8);
        assert_eq!(
            table.backbone.len(),
            65536,
            "8-mer backbone should have 4^8=65536 entries"
        );
        assert_eq!(table.scan_step, 4);
    }

    #[test]
    fn test_small_na_lookup_backbone_size() {
        // For word_size=7: 4^7 = 16384 backbone entries
        let table = SmallNaLookupTable {
            word_length: 7,
            backbone: vec![-1; 16384],
            overflow: vec![],
            pv_array: vec![],
            longest_chain: 0,
            scan_step: 4,
        };
        assert_eq!(table.backbone.len(), 16384);
    }

    #[test]
    fn test_megablast_lookup_table_properties() {
        let table = MbLookupTable {
            word_length: 28,
            lut_word_length: 12,
            discontiguous: false,
            template_length: 0,
            template_type: DiscTemplateType::Contiguous,
            two_templates: false,
            second_template_type: DiscTemplateType::Contiguous,
            hashtable: vec![0; 4_194_304], // 4^12 = 4194304
            hashtable2: vec![],
            next_pos: vec![],
            next_pos2: vec![],
            pv_array: vec![],
            pv_array_bts: crate::stat::PV_ARRAY_BTS as i32,
            longest_chain: 0,
            scan_step: 1,
        };
        assert_eq!(table.word_length, 28);
        assert_eq!(table.lut_word_length, 12);
        assert_eq!(
            table.hashtable.len(),
            4_194_304,
            "Megablast hash size should be 4^12"
        );
    }

    #[test]
    fn megablast_lookup_has_hits_and_retrieve_walk_chains() {
        let mut table = MbLookupTable {
            word_length: 28,
            lut_word_length: 12,
            discontiguous: false,
            template_length: 0,
            template_type: DiscTemplateType::Contiguous,
            two_templates: false,
            second_template_type: DiscTemplateType::Contiguous,
            hashtable: vec![0, 2],
            hashtable2: vec![0, 3],
            next_pos: vec![0, 0, 4, 0, 0],
            next_pos2: vec![0, 0, 0, 5, 0, 0],
            pv_array: vec![0],
            pv_array_bts: crate::stat::PV_ARRAY_BTS as i32,
            longest_chain: 2,
            scan_step: 1,
        };
        table.pv_array[0] = 1 << 1;

        assert_eq!(s_blast_mb_lookup_has_hits(&table, 1), 1);
        assert_eq!(s_blast_mb_lookup_has_hits(&table, 2), 0);

        let mut pairs = Vec::new();
        assert_eq!(s_blast_mb_lookup_retrieve(&table, 1, &mut pairs, 17), 2);
        assert_eq!(
            pairs,
            vec![
                OffsetPair {
                    query_offset: 1,
                    subject_offset: 17,
                },
                OffsetPair {
                    query_offset: 3,
                    subject_offset: 17,
                },
            ]
        );

        pairs.clear();
        assert_eq!(s_blast_mb_lookup_retrieve2(&table, 1, &mut pairs, 23), 2);
        assert_eq!(
            pairs,
            vec![
                OffsetPair {
                    query_offset: 2,
                    subject_offset: 23,
                },
                OffsetPair {
                    query_offset: 4,
                    subject_offset: 23,
                },
            ]
        );
    }

    #[test]
    fn standard_na_lookup_retrieve_uses_backbone_and_overflow_payloads() {
        let mut table = BlastNaLookupTable {
            thick_backbone: vec![
                NaLookupBackboneCell::default(),
                NaLookupBackboneCell {
                    num_used: 3,
                    entries: [4, 8, 12],
                    overflow_cursor: 0,
                },
                NaLookupBackboneCell {
                    num_used: 4,
                    entries: [0; NA_HITS_PER_CELL],
                    overflow_cursor: 1,
                },
            ],
            overflow: vec![99, 16, 20, 24, 28],
            pv: vec![0],
            ..Default::default()
        };
        table.pv[0] = (1 << 1) | (1 << 2);

        assert_eq!(s_blast_lookup_get_num_hits(&table, 1), 3);
        assert_eq!(s_blast_lookup_get_num_hits(&table, 2), 4);
        assert_eq!(s_blast_lookup_get_num_hits(&table, 3), 0);

        let mut pairs = Vec::new();
        s_blast_lookup_retrieve(&table, 1, &mut pairs, 30);
        assert_eq!(
            pairs,
            vec![
                OffsetPair {
                    query_offset: 4,
                    subject_offset: 30,
                },
                OffsetPair {
                    query_offset: 8,
                    subject_offset: 30,
                },
                OffsetPair {
                    query_offset: 12,
                    subject_offset: 30,
                },
            ]
        );

        pairs.clear();
        s_blast_lookup_retrieve(&table, 2, &mut pairs, 40);
        assert_eq!(
            pairs,
            vec![
                OffsetPair {
                    query_offset: 16,
                    subject_offset: 40,
                },
                OffsetPair {
                    query_offset: 20,
                    subject_offset: 40,
                },
                OffsetPair {
                    query_offset: 24,
                    subject_offset: 40,
                },
                OffsetPair {
                    query_offset: 28,
                    subject_offset: 40,
                },
            ]
        );
    }

    #[test]
    fn nucleotide_lookup_finalize_packs_thin_backbone_like_c() {
        let mut thin = vec![
            None,
            Some(vec![0, 1, 11]),
            Some(vec![0, 3, 21, 22, 23]),
            None,
        ];
        let mut small = SmallNaLookupTable {
            word_length: 8,
            backbone: Vec::new(),
            overflow: Vec::new(),
            pv_array: Vec::new(),
            longest_chain: 0,
            scan_step: 1,
        };
        assert_eq!(s_blast_small_na_lookup_finalize(&mut thin, &mut small), 0);
        assert_eq!(small.backbone, vec![-1, 11, -2, -1]);
        assert_eq!(small.overflow, vec![0, 0, 21, 22, 23, -1]);
        assert_eq!(small.longest_chain, 3);
        assert!(thin.iter().all(Option::is_none));

        let mut thin = vec![
            None,
            Some(vec![0, 2, 7, 8]),
            Some(vec![0, 4, 31, 32, 33, 34]),
        ];
        let mut standard = BlastNaLookupTable::default();
        s_blast_na_lookup_finalize(&mut thin, &mut standard);

        assert_eq!(standard.backbone_size, 3);
        assert_eq!(standard.longest_chain, 4);
        assert_eq!(standard.thick_backbone[1].num_used, 2);
        assert_eq!(standard.thick_backbone[1].entries, [7, 8, 0]);
        assert_eq!(standard.thick_backbone[2].num_used, 4);
        assert_eq!(standard.thick_backbone[2].overflow_cursor, 0);
        assert_eq!(standard.overflow, vec![31, 32, 33, 34]);
        assert_eq!(standard.overflow_size, 4);
        assert_eq!(s_blast_lookup_get_num_hits(&standard, 1), 2);
        assert!(thin.iter().all(Option::is_none));
    }

    #[test]
    fn lookup_index_query_exact_matches_builds_thin_backbone_like_c() {
        let mut thin = vec![None; 1 << 4];
        let ranges = [crate::util::SSeqRange { left: 0, right: 5 }];

        blast_lookup_index_query_exact_matches(
            &mut thin,
            3,
            BITS_PER_NUC as i32,
            2,
            &[0, 1, 0, 1, 4, 0],
            &ranges,
        );

        let ac = (0usize << 2) | 1;
        assert_eq!(thin[ac], Some(vec![8, 2, 0, 2]));
        assert!(thin[(1usize << 2) | 0].is_some());
        assert_eq!(
            thin.iter()
                .filter_map(|cell| cell.as_ref().map(|cell| cell[1]))
                .sum::<i32>(),
            3
        );
    }

    #[test]
    fn small_na_lookup_table_new_builds_backbone_and_scan_metadata() {
        let query = crate::encoding::encode_blastna_sequence(b"ACGTACGT");
        let range = [crate::util::SSeqRange {
            left: 0,
            right: query.len() as i32 - 1,
        }];
        let options = crate::options::LookupTableOptions {
            word_size: 4,
            ..Default::default()
        };
        let mut lookup = None;

        assert_eq!(
            blast_small_na_lookup_table_new(&query, &range, &mut lookup, &options, 4),
            0
        );
        let lookup = lookup.expect("SmallNa lookup");
        assert_eq!(lookup.word_length, 4);
        assert_eq!(lookup.scan_step, 1);
        assert_eq!(lookup.backbone.len(), 1usize << (BITS_PER_NUC * 4));
        assert!(lookup.longest_chain > 0);

        let index = s_compute_table_index(4, BITS_PER_NUC as i32, &query[0..4]) as i32;
        let wrap = LookupTableWrap::SmallNa(lookup);
        assert!(s_small_na_lookup(Some(&wrap), index, 0));
        assert!(s_small_na_lookup(Some(&wrap), index, 4));
        assert!(!s_small_na_lookup(Some(&wrap), index, 1));
    }

    #[test]
    fn standard_na_lookup_table_new_builds_thick_backbone_and_pv() {
        let query = crate::encoding::encode_blastna_sequence(b"ACGTACGTAC");
        let range = [crate::util::SSeqRange {
            left: 0,
            right: query.len() as i32 - 1,
        }];
        let options = crate::options::LookupTableOptions {
            word_size: 5,
            ..Default::default()
        };
        let mut lookup = None;

        assert_eq!(
            blast_na_lookup_table_new(&query, &range, &mut lookup, &options, 4),
            0
        );
        let lookup = lookup.expect("standard Na lookup");
        assert_eq!(lookup.word_length, 5);
        assert_eq!(lookup.lut_word_length, 4);
        assert_eq!(lookup.scan_step, 2);
        assert_eq!(lookup.backbone_size, 1i32 << (BITS_PER_NUC * 4));
        assert_eq!(lookup.mask, lookup.backbone_size - 1);
        assert!(lookup.longest_chain > 0);

        let index = s_compute_table_index(4, BITS_PER_NUC as i32, &query[0..4]) as i32;
        assert!(s_na_lookup(Some(&lookup), index, 0));
        assert!(s_na_lookup(Some(&lookup), index, 4));
        assert!(!s_na_lookup(Some(&lookup), index, 1));
        assert!(s_blast_lookup_get_num_hits(&lookup, index) >= 2);
        assert!(lookup.pv[(index as usize) >> crate::stat::PV_ARRAY_BTS] != 0);
    }

    #[test]
    fn megablast_lookup_table_new_builds_contiguous_lookup_like_c() {
        let query = crate::encoding::encode_blastna_sequence(b"ACGTACGTACGT");
        let range = [crate::util::SSeqRange {
            left: 0,
            right: query.len() as i32 - 1,
        }];
        let options = crate::options::LookupTableOptions {
            word_size: 9,
            program_number: crate::program::BLASTN,
            ..Default::default()
        };
        let mut lookup = None;

        assert_eq!(
            blast_mb_lookup_table_new(&query, &range, &mut lookup, &options, 21_000, 9, None),
            0
        );
        let lookup = lookup.expect("megablast lookup");
        assert_eq!(lookup.word_length, 9);
        assert_eq!(lookup.lut_word_length, 9);
        assert_eq!(lookup.scan_step, 1);
        assert_eq!(lookup.hashtable.len(), 1usize << (BITS_PER_NUC * 9));
        assert_eq!(lookup.next_pos.len(), query.len() + 1);
        assert!(!lookup.pv_array.is_empty());
        assert!(lookup.pv_array_bts >= crate::stat::PV_ARRAY_BTS as i32);

        let index = s_compute_table_index(9, BITS_PER_NUC as i32, &query[0..9]) as i32;
        let wrap = LookupTableWrap::Megablast(lookup);
        assert!(s_mb_lookup(Some(&wrap), index, 0));
        assert!(!s_mb_lookup(Some(&wrap), index, 1));
    }

    #[test]
    fn nucleotide_lookup_table_new_rejects_invalid_lut_width() {
        let query = crate::encoding::encode_blastna_sequence(b"ACGTACGT");
        let range = [crate::util::SSeqRange {
            left: 0,
            right: query.len() as i32 - 1,
        }];
        let options = crate::options::LookupTableOptions {
            word_size: 3,
            ..Default::default()
        };
        let mut small = Some(SmallNaLookupTable {
            word_length: 1,
            backbone: vec![0],
            overflow: Vec::new(),
            pv_array: Vec::new(),
            longest_chain: 0,
            scan_step: 1,
        });
        let mut standard = Some(BlastNaLookupTable::default());

        assert_eq!(
            blast_small_na_lookup_table_new(&query, &range, &mut small, &options, 4),
            -1
        );
        assert!(small.is_none());
        assert_eq!(
            blast_na_lookup_table_new(&query, &range, &mut standard, &options, 4),
            -1
        );
        assert!(standard.is_none());
    }

    #[test]
    fn translated_standard_na_lookup_destructor_drops_owned_table() {
        let mut lookup = Some(BlastNaLookupTable {
            thick_backbone: vec![NaLookupBackboneCell {
                num_used: 1,
                entries: [3, 0, 0],
                overflow_cursor: 0,
            }],
            overflow: vec![4],
            pv: vec![1],
            ..Default::default()
        });
        assert!(blast_na_lookup_table_destruct(&mut lookup).is_none());
        assert!(lookup.is_none());
    }

    #[test]
    fn nucleotide_hash_lookup_finalize_packs_inline_and_overflow_like_c() {
        let mask = 3;
        let bucket5 = fnv_hash(&5u32.to_ne_bytes(), mask) as usize;
        let bucket9 = fnv_hash(&6u32.to_ne_bytes(), mask) as usize;
        let bucket10 = fnv_hash(&10u32.to_ne_bytes(), mask) as usize;
        assert_ne!(bucket5, bucket9);
        assert_eq!(bucket9, bucket10);

        let mut thin = vec![BackboneCell::default(); 4];
        thin[bucket5] = BackboneCell {
            word: 5,
            offset: 1,
            num_offsets: 2,
            next: None,
        };
        thin[bucket9] = BackboneCell {
            word: 6,
            offset: 3,
            num_offsets: 4,
            next: Some(Box::new(BackboneCell {
                word: 10,
                offset: 7,
                num_offsets: 6,
                next: None,
            })),
        };
        let offsets = vec![0, 2, 0, 4, 5, 6, 0, 8, 9, 10, 11, 12, 0];
        let mut lookup = BlastNaHashLookupTable {
            pv_array_bts: crate::stat::PV_ARRAY_BTS as i32,
            ..Default::default()
        };

        s_blast_na_hash_lookup_finalize(&mut thin, &offsets, &mut lookup);

        assert_eq!(lookup.backbone_size, 4);
        assert_eq!(lookup.longest_chain, 10);
        assert_eq!(lookup.thick_backbone[bucket5].num_words, 1);
        assert_eq!(lookup.thick_backbone[bucket5].num_offsets, [2, 0, 0]);
        assert_eq!(lookup.thick_backbone[bucket5].words, [5, 0, 0]);
        assert_eq!(lookup.thick_backbone[bucket5].offsets[..2], [0, 1]);
        assert_eq!(lookup.thick_backbone[bucket9].num_words, 2);
        assert_eq!(lookup.thick_backbone[bucket9].words, [6, 10, 0]);
        assert_eq!(lookup.thick_backbone[bucket9].offsets[0], 0);
        assert_eq!(
            lookup.overflow,
            vec![6, 4, 2, 3, 4, 5, 10, 6, 6, 7, 8, 9, 10, 11]
        );
        assert_eq!(lookup.offsets_size, lookup.overflow.len() as i32);
        assert!(lookup.pv[0] & (1 << 5) != 0);
        assert!(lookup.pv[0] & (1 << 6) != 0);
        assert!(lookup.pv[0] & (1 << 10) != 0);
        assert!(thin[bucket9].next.is_none());

        lookup.mask = mask as i32;
        let mut pairs = Vec::new();
        assert_eq!(
            s_blast_na_hash_lookup_retrieve_hits(&lookup, 5, 21, &mut pairs),
            2
        );
        assert_eq!(
            pairs,
            vec![
                OffsetPair {
                    query_offset: 0,
                    subject_offset: 21,
                },
                OffsetPair {
                    query_offset: 1,
                    subject_offset: 21,
                },
            ]
        );

        pairs.clear();
        assert_eq!(
            s_blast_na_hash_lookup_retrieve_hits(&lookup, 10, 27, &mut pairs),
            6
        );
        assert_eq!(pairs[0].query_offset, 6);
        assert_eq!(pairs[5].query_offset, 11);
        assert_eq!(
            s_blast_na_hash_lookup_retrieve_hits(&lookup, 11, 27, &mut pairs),
            0
        );
    }

    #[test]
    fn na_hash_scan_subject_any_retrieves_16_base_hits() {
        let subject = crate::encoding::encode_blastna_sequence(b"ACGTACGTACGTACGTAC");
        let packed = crate::encoding::pack_ncbi2na_bases(&subject);
        let word = u32::from_be_bytes([packed[0], packed[1], packed[2], packed[3]]);
        let mut lookup = BlastNaHashLookupTable {
            mask: 15,
            word_length: 16,
            lut_word_length: 16,
            scan_step: 1,
            thick_backbone: vec![NaHashLookupBackboneCell::default(); 16],
            pv_array_bts: crate::stat::PV_ARRAY_BTS as i32,
            ..Default::default()
        };
        let hash_index = fnv_hash(&word.to_ne_bytes(), lookup.mask as u32) as usize;
        lookup.thick_backbone[hash_index] = NaHashLookupBackboneCell {
            num_words: 1,
            num_offsets: [1, 0, 0],
            words: [word, 0, 0],
            offsets: [4, 0, 0, 0, 0, 0, 0, 0, 0],
        };
        set_pv_bit(&mut lookup.pv, word as usize, crate::stat::PV_ARRAY_BTS);

        let hits = s_blast_na_hash_scan_subject_any(&lookup, &packed, subject.len() as i32, 0, 2)
            .expect("represented NaHash scan");
        assert_eq!(
            hits,
            vec![OffsetPair {
                query_offset: 4,
                subject_offset: 0,
            }]
        );
        assert!(
            s_blast_na_hash_scan_subject_any(&lookup, &packed, subject.len() as i32, 0, 4)
                .is_none()
        );
    }

    #[test]
    fn na_hash_lookup_table_new_builds_fnv_backbone_and_scan_metadata() {
        let query = crate::encoding::encode_blastna_sequence(b"ACGTACGTACGTACGTACGT");
        let range = [crate::util::SSeqRange {
            left: 0,
            right: query.len() as i32 - 1,
        }];
        let options = crate::options::LookupTableOptions {
            word_size: 16,
            db_filter: false,
            ..Default::default()
        };
        let mut lookup = None;

        assert_eq!(
            blast_na_hash_lookup_table_new(&query, &range, &mut lookup, &options, None, 1),
            0
        );
        let lookup = lookup.expect("NaHash lookup");
        assert_eq!(lookup.word_length, 16);
        assert_eq!(lookup.lut_word_length, 16);
        assert_eq!(lookup.scan_step, 1);
        assert_eq!(lookup.backbone_size, 256);
        assert_eq!(lookup.mask, 255);
        assert!(lookup.longest_chain > 0);

        let packed = crate::encoding::pack_ncbi2na_bases(&query);
        let hits = s_blast_na_hash_scan_subject_any(&lookup, &packed, query.len() as i32, 0, 4)
            .expect("NaHash scan");
        assert!(
            hits.iter()
                .any(|hit| hit.query_offset == 0 && hit.subject_offset == 0),
            "expected self hit from constructed NaHash lookup: {hits:?}"
        );
    }

    #[test]
    fn na_hash_lookup_table_new_indexes_leading_windows_for_long_words() {
        let query = crate::encoding::encode_blastna_sequence(b"ACGTACGTACGTACGTACGT");
        let range = [crate::util::SSeqRange {
            left: 0,
            right: query.len() as i32 - 1,
        }];
        let options = crate::options::LookupTableOptions {
            word_size: 20,
            db_filter: false,
            ..Default::default()
        };
        let mut lookup = None;

        assert_eq!(
            blast_na_hash_lookup_table_new(&query, &range, &mut lookup, &options, None, 1),
            0
        );
        let lookup = lookup.expect("long-word NaHash lookup");
        assert_eq!(lookup.word_length, 20);
        assert_eq!(lookup.lut_word_length, 16);
        assert_eq!(lookup.scan_step, 5);

        let packed = crate::encoding::pack_ncbi2na_bases(&query);
        let hits = s_blast_na_hash_scan_subject_any(&lookup, &packed, query.len() as i32, 0, 0)
            .expect("NaHash scan at leading word");
        assert!(
            hits.iter()
                .any(|hit| hit.query_offset == 0 && hit.subject_offset == 0),
            "C indexes the leading 16-base window for word_size > 16: {hits:?}"
        );
    }

    #[test]
    fn na_hash_lookup_table_new_applies_database_word_filter() {
        struct TestSeqSrc {
            seqs: Vec<Vec<u8>>,
            reset: bool,
        }

        impl crate::seqsrc::BlastSeqSource for TestSeqSrc {
            fn num_seqs(&self) -> i32 {
                self.seqs.len() as i32
            }

            fn total_length(&self) -> i64 {
                self.seqs.iter().map(|seq| seq.len() as i64).sum()
            }

            fn max_seq_len(&self) -> i32 {
                self.seqs.iter().map(Vec::len).max().unwrap_or(0) as i32
            }

            fn avg_seq_len(&self) -> i32 {
                if self.seqs.is_empty() {
                    0
                } else {
                    (self.total_length() / self.seqs.len() as i64) as i32
                }
            }

            fn name(&self) -> &str {
                "test"
            }

            fn is_protein(&self) -> bool {
                false
            }

            fn seq_len(&self, oid: i32) -> i32 {
                self.seqs[oid as usize].len() as i32
            }

            fn get_sequence(
                &self,
                arg: &crate::seqsrc::GetSeqArg,
            ) -> Option<crate::seqsrc::SeqData> {
                self.seqs
                    .get(arg.oid as usize)
                    .map(|seq| crate::seqsrc::SeqData {
                        sequence: seq.clone(),
                        length: seq.len() as i32 * 4,
                    })
            }

            fn iter_oids(&self) -> Box<dyn Iterator<Item = i32> + '_> {
                Box::new(0..self.seqs.len() as i32)
            }

            fn reset_chunk_iterator(&mut self) {
                self.reset = true;
            }
        }

        let query = crate::encoding::encode_blastna_sequence(b"ACGTACGTACGTACGT");
        let packed_query = crate::encoding::pack_ncbi2na_bases(&query);
        let range = [crate::util::SSeqRange {
            left: 0,
            right: query.len() as i32 - 1,
        }];
        let options = crate::options::LookupTableOptions {
            word_size: 16,
            db_filter: true,
            max_db_word_count: 2,
            ..Default::default()
        };
        let mut missing_seqsrc = None;
        assert_eq!(
            blast_na_hash_lookup_table_new(&query, &range, &mut missing_seqsrc, &options, None, 1),
            -1
        );

        let mut seqsrc = TestSeqSrc {
            seqs: vec![packed_query.clone(), packed_query],
            reset: false,
        };
        let mut lookup = None;
        assert_eq!(
            blast_na_hash_lookup_table_new(
                &query,
                &range,
                &mut lookup,
                &options,
                Some(&mut seqsrc),
                4,
            ),
            0
        );
        assert!(seqsrc.reset);
        let lookup = lookup.expect("filtered NaHash lookup");
        assert_eq!(lookup.scan_step, 1);
        assert_eq!(lookup.longest_chain, 0);
        assert!(lookup.thick_backbone.iter().all(|cell| cell.num_words == 0));
    }

    #[test]
    fn translated_hash_na_lookup_destructor_drops_owned_table() {
        let mut lookup = Some(BlastNaHashLookupTable {
            thick_backbone: vec![NaHashLookupBackboneCell {
                num_words: 1,
                num_offsets: [1, 0, 0],
                words: [7, 0, 0],
                offsets: [3, 0, 0, 0, 0, 0, 0, 0, 0],
            }],
            overflow: vec![7, 1, 3],
            pv: vec![1 << 7],
            ..Default::default()
        });
        assert!(blast_na_hash_lookup_table_destruct(&mut lookup).is_none());
        assert!(lookup.is_none());
    }

    #[test]
    fn na_hash_lookup_remove_poly_a_words_clears_c_poly_a_bits() {
        let mut lookup = BlastNaHashLookupTable {
            lut_word_length: 16,
            pv_array_bts: crate::stat::PV_ARRAY_BTS as i32,
            pv: vec![u32::MAX; 1],
            ..Default::default()
        };

        assert_eq!(s_na_hash_lookup_remove_poly_a_words(Some(&mut lookup)), 0);
        for index in [0u32, 1, 2, 3, 4, 8, 12] {
            assert_eq!(
                lookup.pv[(index >> crate::stat::PV_ARRAY_BTS) as usize]
                    & (1 << (index & crate::stat::PV_ARRAY_MASK)),
                0
            );
        }
        assert_eq!(s_na_hash_lookup_remove_poly_a_words(None), -1);

        let mut unsupported_width = BlastNaHashLookupTable {
            lut_word_length: 15,
            pv_array_bts: crate::stat::PV_ARRAY_BTS as i32,
            pv: vec![u32::MAX; 1],
            ..Default::default()
        };
        assert_eq!(
            s_na_hash_lookup_remove_poly_a_words(Some(&mut unsupported_width)),
            -1
        );

        let mut compressed_pv = BlastNaHashLookupTable {
            lut_word_length: 16,
            pv_array_bts: crate::stat::PV_ARRAY_BTS as i32 + 1,
            pv: vec![u32::MAX; 1],
            ..Default::default()
        };
        assert_eq!(
            s_na_hash_lookup_remove_poly_a_words(Some(&mut compressed_pv)),
            -1
        );
    }

    #[test]
    fn na_hash_lookup_fill_pv_sets_words_for_unmasked_ranges() {
        let mut lookup = BlastNaHashLookupTable {
            word_length: 3,
            lut_word_length: 3,
            pv_array_bts: crate::stat::PV_ARRAY_BTS as i32,
            pv: vec![0],
            ..Default::default()
        };
        let ranges = [crate::util::SSeqRange { left: 0, right: 6 }];

        assert_eq!(
            s_na_hash_lookup_fill_pv(&[0, 1, 2, 3, 4, 0, 1], &ranges, &mut lookup),
            0
        );

        let first_word = (0u32 << 4) | (1 << 2) | 2;
        let second_word = (1u32 << 4) | (2 << 2) | 3;
        let post_ambiguity_word = (0u32 << 4) | (0 << 2) | 1;
        assert!(lookup.pv[0] & (1 << first_word) != 0);
        assert!(lookup.pv[0] & (1 << second_word) != 0);
        assert_eq!(lookup.pv[0] & (1 << post_ambiguity_word), 0);
    }

    #[test]
    fn na_hash_lookup_count_words_in_subject_updates_sparse_counts_like_c() {
        let lookup = BlastNaHashLookupTable {
            lut_word_length: 16,
            pv_array_bts: crate::stat::PV_ARRAY_BTS as i32,
            pv: vec![1],
            ..Default::default()
        };
        let mut counts = crate::util::blast_sparse_uint1_array_new(Some(&[1]), 32)
            .expect("sparse counter for word zero");

        assert_eq!(
            s_na_hash_lookup_count_words_in_subject_16_1(
                Some(&[0, 0, 0, 0, 0]),
                18,
                Some(&lookup),
                Some(&mut counts),
                2,
            ),
            0
        );
        assert_eq!(counts.values, vec![2]);

        assert_eq!(
            s_na_hash_lookup_count_words_in_subject_16_1(
                Some(&[0, 0, 0, 0, 0]),
                18,
                Some(&lookup),
                Some(&mut counts),
                2,
            ),
            0
        );
        assert_eq!(counts.values, vec![2]);
        assert_eq!(
            s_na_hash_lookup_count_words_in_subject_16_1(
                None,
                18,
                Some(&lookup),
                Some(&mut counts),
                2
            ),
            -1
        );
    }

    #[test]
    fn megablast_remove_poly_a_words_clears_small_table_poly_a_entries() {
        let mut table = MbLookupTable {
            word_length: 4,
            lut_word_length: 4,
            discontiguous: false,
            template_length: 0,
            template_type: DiscTemplateType::Contiguous,
            two_templates: false,
            second_template_type: DiscTemplateType::Contiguous,
            hashtable: vec![7; 1 << 8],
            hashtable2: Vec::new(),
            next_pos: Vec::new(),
            next_pos2: Vec::new(),
            pv_array: Vec::new(),
            pv_array_bts: crate::stat::PV_ARRAY_BTS as i32,
            longest_chain: 0,
            scan_step: 1,
        };

        assert_eq!(s_remove_poly_a_words(&mut table), 0);
        assert_eq!(table.hashtable[0], 0);
        assert_eq!(table.hashtable[(1 << 8) - 1], 0);
        assert_eq!(table.hashtable[1], 7);
    }

    #[test]
    fn megablast_fill_pv_sets_contiguous_words_and_skips_ambiguities() {
        let mut table = MbLookupTable {
            word_length: 3,
            lut_word_length: 3,
            discontiguous: false,
            template_length: 0,
            template_type: DiscTemplateType::Contiguous,
            two_templates: false,
            second_template_type: DiscTemplateType::Contiguous,
            hashtable: vec![0; 1 << 6],
            hashtable2: Vec::new(),
            next_pos: Vec::new(),
            next_pos2: Vec::new(),
            pv_array: vec![0],
            pv_array_bts: crate::stat::PV_ARRAY_BTS as i32,
            longest_chain: 0,
            scan_step: 1,
        };
        let ranges = [crate::util::SSeqRange { left: 0, right: 6 }];

        assert_eq!(
            s_fill_pv(&[0, 1, 2, 3, 4, 0, 1], &ranges, &mut table, None),
            0
        );

        let first_word = (0u32 << 4) | (1 << 2) | 2;
        let second_word = (1u32 << 4) | (2 << 2) | 3;
        let post_ambiguity_word = (0u32 << 4) | (0 << 2) | 1;
        assert!(table.pv_array[0] & (1 << first_word) != 0);
        assert!(table.pv_array[0] & (1 << second_word) != 0);
        assert_eq!(table.pv_array[0] & (1 << post_ambiguity_word), 0);
    }

    #[test]
    fn megablast_fill_contig_table_chains_one_based_offsets_like_c() {
        let mut table = MbLookupTable {
            word_length: 3,
            lut_word_length: 3,
            discontiguous: false,
            template_length: 0,
            template_type: DiscTemplateType::Contiguous,
            two_templates: false,
            second_template_type: DiscTemplateType::Contiguous,
            hashtable: vec![0; 1 << 6],
            hashtable2: Vec::new(),
            next_pos: Vec::new(),
            next_pos2: Vec::new(),
            pv_array: vec![0],
            pv_array_bts: crate::stat::PV_ARRAY_BTS as i32,
            longest_chain: 0,
            scan_step: 1,
        };
        let ranges = [crate::util::SSeqRange { left: 0, right: 5 }];
        let options = crate::options::LookupTableOptions {
            word_size: 3,
            ..Default::default()
        };

        assert_eq!(
            s_fill_contig_mb_table(&[0, 1, 2, 0, 1, 2], &ranges, &mut table, &options, None),
            0
        );

        let word = (0usize << 4) | (1 << 2) | 2;
        assert_eq!(table.hashtable[word], 4);
        assert_eq!(table.next_pos[4], 1);
        assert!(table.pv_array[0] & (1 << word) != 0);
        assert_eq!(table.longest_chain, 2);
    }

    #[test]
    fn megablast_fill_contig_table_honors_db_filter_counts() {
        let mut table = MbLookupTable {
            word_length: 3,
            lut_word_length: 3,
            discontiguous: false,
            template_length: 0,
            template_type: DiscTemplateType::Contiguous,
            two_templates: false,
            second_template_type: DiscTemplateType::Contiguous,
            hashtable: vec![0; 1 << 6],
            hashtable2: Vec::new(),
            next_pos: Vec::new(),
            next_pos2: Vec::new(),
            pv_array: vec![u32::MAX],
            pv_array_bts: crate::stat::PV_ARRAY_BTS as i32,
            longest_chain: 0,
            scan_step: 1,
        };
        let ranges = [crate::util::SSeqRange { left: 0, right: 2 }];
        let options = crate::options::LookupTableOptions {
            word_size: 3,
            db_filter: true,
            max_db_word_count: 10,
            ..Default::default()
        };
        let word = (0usize << 4) | (1 << 2) | 2;
        let mut counts = vec![0u8; (1 << 6) / 2];
        counts[word / 2] = 10 << 4;

        assert_eq!(
            s_fill_contig_mb_table(&[0, 1, 2], &ranges, &mut table, &options, Some(&counts)),
            0
        );

        assert_eq!(table.hashtable[word], 0);
        assert_eq!(table.pv_array, vec![0]);
    }

    #[test]
    fn megablast_count_words_in_subject_updates_packed_nibble_counts() {
        let table = MbLookupTable {
            word_length: 16,
            lut_word_length: 16,
            discontiguous: false,
            template_length: 0,
            template_type: DiscTemplateType::Contiguous,
            two_templates: false,
            second_template_type: DiscTemplateType::Contiguous,
            hashtable: vec![0; 1 << 6],
            hashtable2: Vec::new(),
            next_pos: Vec::new(),
            next_pos2: Vec::new(),
            pv_array: vec![1],
            pv_array_bts: crate::stat::PV_ARRAY_BTS as i32,
            longest_chain: 0,
            scan_step: 1,
        };
        let mut counts = vec![0u8; table.hashtable.len() / 2];

        assert_eq!(
            s_mb_count_words_in_subject_16_1(
                Some(&[0, 0, 0, 0, 0]),
                18,
                Some(&table),
                Some(&mut counts),
                2
            ),
            0
        );
        assert_eq!(counts[0], 0x20);

        assert_eq!(
            s_mb_count_words_in_subject_16_1(
                Some(&[0, 0, 0, 0, 0]),
                18,
                Some(&table),
                Some(&mut counts),
                2
            ),
            0
        );
        assert_eq!(counts[0], 0x20);
    }

    #[test]
    fn na_hash_lookup_thread_data_new_matches_c_thread_layout() {
        struct EmptySeqSrc;

        impl crate::seqsrc::BlastSeqSource for EmptySeqSrc {
            fn num_seqs(&self) -> i32 {
                0
            }

            fn total_length(&self) -> i64 {
                0
            }

            fn max_seq_len(&self) -> i32 {
                0
            }

            fn avg_seq_len(&self) -> i32 {
                0
            }

            fn name(&self) -> &str {
                "empty"
            }

            fn is_protein(&self) -> bool {
                false
            }

            fn seq_len(&self, _oid: i32) -> i32 {
                0
            }

            fn get_sequence(
                &self,
                _arg: &crate::seqsrc::GetSeqArg,
            ) -> Option<crate::seqsrc::SeqData> {
                None
            }

            fn iter_oids(&self) -> Box<dyn Iterator<Item = i32> + '_> {
                Box::new(std::iter::empty())
            }
        }

        let lookup = BlastNaHashLookupTable {
            lut_word_length: 4,
            pv: {
                let mut pv = vec![0; 8];
                pv[0] = 0b1011;
                pv
            },
            pv_array_bts: crate::stat::PV_ARRAY_BTS as i32,
            ..BlastNaHashLookupTable::default()
        };
        let seq_src = crate::seqsrc::BlastSeqSrc {
            source: Some(std::sync::Arc::new(EmptySeqSrc)),
            init_error: None,
            delete_fn: None,
            copy_fn: None,
        };

        let mut thread_data =
            na_hash_lookup_thread_data_new(3, Some(&lookup), Some(&seq_src)).expect("thread data");
        assert_eq!(thread_data.num_threads, 3);
        assert_eq!(thread_data.seq_arg.len(), 3);
        assert!(thread_data
            .seq_arg
            .iter()
            .all(|arg| arg.encoding == crate::seqsrc::SeqEncoding::Protein));
        assert_eq!(
            thread_data
                .itr
                .iter()
                .map(|itr| itr.chunk_sz)
                .collect::<Vec<_>>(),
            vec![1, 1, 1]
        );
        assert_eq!(thread_data.seq_src.len(), 3);
        assert_eq!(thread_data.word_counts.len(), 3);
        assert_eq!(thread_data.word_counts[0].num_elements, 3);
        assert_eq!(
            thread_data.word_counts[1].counts,
            thread_data.word_counts[0].counts
        );
        thread_data.word_counts[1].values[0] = 7;
        assert_eq!(thread_data.word_counts[0].values[0], 0);

        let mut slot = Some(thread_data);
        assert!(na_hash_lookup_thread_data_free(&mut slot).is_none());
        assert!(slot.is_none());
        assert!(na_hash_lookup_thread_data_new(0, Some(&lookup), Some(&seq_src)).is_none());
        assert!(na_hash_lookup_thread_data_new(1, None, Some(&seq_src)).is_none());
        assert!(na_hash_lookup_thread_data_new(1, Some(&lookup), None).is_none());
    }

    #[test]
    fn megablast_scan_subject_for_word_counts_walks_seqsrc() {
        struct TestSeqSrc {
            seqs: Vec<Vec<u8>>,
            reset: bool,
        }

        impl crate::seqsrc::BlastSeqSource for TestSeqSrc {
            fn num_seqs(&self) -> i32 {
                self.seqs.len() as i32
            }

            fn total_length(&self) -> i64 {
                self.seqs.iter().map(|seq| seq.len() as i64).sum()
            }

            fn max_seq_len(&self) -> i32 {
                self.seqs.iter().map(Vec::len).max().unwrap_or(0) as i32
            }

            fn avg_seq_len(&self) -> i32 {
                if self.seqs.is_empty() {
                    0
                } else {
                    (self.total_length() / self.seqs.len() as i64) as i32
                }
            }

            fn name(&self) -> &str {
                "test"
            }

            fn is_protein(&self) -> bool {
                false
            }

            fn seq_len(&self, oid: i32) -> i32 {
                self.seqs[oid as usize].len() as i32
            }

            fn get_sequence(
                &self,
                arg: &crate::seqsrc::GetSeqArg,
            ) -> Option<crate::seqsrc::SeqData> {
                assert_eq!(arg.encoding, crate::seqsrc::SeqEncoding::Protein);
                self.seqs
                    .get(arg.oid as usize)
                    .map(|seq| crate::seqsrc::SeqData {
                        sequence: seq.clone(),
                        length: 18,
                    })
            }

            fn iter_oids(&self) -> Box<dyn Iterator<Item = i32> + '_> {
                Box::new(0..self.seqs.len() as i32)
            }

            fn reset_chunk_iterator(&mut self) {
                self.reset = true;
            }
        }

        let table = MbLookupTable {
            word_length: 16,
            lut_word_length: 16,
            discontiguous: false,
            template_length: 0,
            template_type: DiscTemplateType::Contiguous,
            two_templates: false,
            second_template_type: DiscTemplateType::Contiguous,
            hashtable: vec![0; 1 << 6],
            hashtable2: Vec::new(),
            next_pos: Vec::new(),
            next_pos2: Vec::new(),
            pv_array: vec![1],
            pv_array_bts: crate::stat::PV_ARRAY_BTS as i32,
            longest_chain: 0,
            scan_step: 1,
        };
        let mut counts = vec![0u8; table.hashtable.len() / 2];
        let mut seqsrc = TestSeqSrc {
            seqs: vec![vec![0, 0, 0, 0, 0], vec![0, 0, 0, 0, 0]],
            reset: false,
        };

        assert_eq!(
            s_scan_subject_for_word_counts(Some(&mut seqsrc), Some(&table), Some(&mut counts), 3),
            0
        );
        assert!(seqsrc.reset);
        assert_eq!(counts[0], 0x30);
    }

    #[test]
    fn test_aa_lookup_table_properties() {
        let table = AaLookupTable {
            word_length: 3,
            threshold: 11.0,
            backbone: vec![vec![]; 28 * 28 * 28], // 28^3 for NCBIstdaa
            pv_array: vec![],
        };
        assert_eq!(table.word_length, 3);
        assert_eq!(table.threshold, 11.0);
        assert_eq!(
            table.backbone.len(),
            21952,
            "3-mer AA backbone should have 28^3=21952 entries"
        );
    }

    #[test]
    fn rps_lookup_table_new_builds_presence_vector_and_buckets() {
        let info = BlastRpsInfo {
            alphabet_size: 26,
            wordsize: 3,
            rps_backbone: vec![
                RpsBackboneCell::default(),
                RpsBackboneCell {
                    num_used: 1,
                    offset_pairs: vec![OffsetPair {
                        query_offset: 2,
                        subject_offset: 4,
                    }],
                },
                RpsBackboneCell::default(),
                RpsBackboneCell {
                    num_used: 2,
                    offset_pairs: vec![
                        OffsetPair {
                            query_offset: 6,
                            subject_offset: 8,
                        },
                        OffsetPair {
                            query_offset: 10,
                            subject_offset: 12,
                        },
                    ],
                },
            ],
            rps_pssm: vec![vec![1, 2], vec![3, 4]],
        };

        let mut lookup = None;
        assert_eq!(rps_lookup_table_new(Some(&info), &mut lookup), 0);
        let lookup_ref = lookup.as_ref().unwrap();
        assert_eq!(lookup_ref.alphabet_size, 26);
        assert_eq!(lookup_ref.wordsize, 3);
        assert_eq!(lookup_ref.charsize, crate::util::ilog2(26) + 1);
        assert_eq!(lookup_ref.num_buckets, 2);
        assert_eq!(lookup_ref.bucket_array[0].num_alloc, 1000);
        assert_eq!(lookup_ref.bucket_array[1].num_used, 2);
        assert_eq!(lookup_ref.pv[0] & (1u32 << 1), 1u32 << 1);
        assert_eq!(lookup_ref.pv[0] & (1u32 << 3), 1u32 << 3);
        assert_eq!(lookup_ref.rps_pssm, info.rps_pssm);

        assert!(rps_lookup_table_destruct(&mut lookup).is_none());
        assert!(lookup.is_none());
    }

    #[test]
    fn rps_profile_database_new_builds_lookup_from_owned_consensus_profiles() {
        let mut pssm_values = vec![-8; 25 * 26];
        let mut query = vec![25u8; 24];
        query[3..6].copy_from_slice(&[1, 2, 3]);
        query[10..13].copy_from_slice(&[4, 5, 6]);
        for (subject_offset, residue) in query.iter().enumerate().take(13).skip(3) {
            let profile_query_offset = 8 + subject_offset as i32;
            pssm_values[profile_query_offset as usize * 26 + *residue as usize] = 5;
        }
        let rps_info = crate::hspstream::RpsTracebackInfo {
            profile_header: crate::hspstream::RpsProfileHeader {
                magic_number: crate::hspstream::RPS_MAGIC_NUM,
                num_profiles: 2,
                start_offsets: vec![0, 8, 24],
                pssm_values,
            },
            freq_ratios_header: None,
            karlin_k: vec![0.1, 0.25],
        };
        let consensus_sequences = vec![
            vec![7, 7, 7, 7, 7, 7, 7, 7],
            vec![7, 7, 7, 1, 2, 3, 7, 7, 7, 7, 4, 5, 6, 7, 7, 7],
        ];

        let mut profile_db =
            blast_rps_profile_database_new(rps_info.clone(), consensus_sequences, 3)
                .expect("owned RPS profile database");

        assert_eq!(profile_db.lookup.alphabet_size, 26);
        assert_eq!(profile_db.lookup.wordsize, 3);
        assert_eq!(profile_db.lookup.backbone_size, 32768);
        assert!(profile_db.lookup.num_buckets >= 2);
        let first_word_index = s_compute_table_index(3, 5, &[1, 2, 3]);
        assert_eq!(
            profile_db.lookup.rps_backbone[first_word_index].offset_pairs,
            vec![OffsetPair {
                query_offset: 13,
                subject_offset: 0,
            }]
        );

        let stream =
            crate::hspstream::blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);
        let mut init_hitlist = crate::extend::InitHitList::new();
        let mut stats = crate::diagnostics::UngappedStats::default();
        assert_eq!(
            blast_rps_profile_database_scan_query_write_hsp_stream(
                &mut profile_db,
                &query,
                &stream,
                0,
                40,
                6,
                1,
                4,
                Some(&mut init_hitlist),
                Some(&mut stats),
            ),
            0
        );
        assert_eq!(init_hitlist.total(), 1);
        assert_eq!(stats.lookup_hits, 2);
        assert_eq!(stats.init_extends, 1);
        assert_eq!(stats.good_init_extends, 1);

        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[32]);
        let mut results = crate::hspstream::HspResults::new(1);
        let mut callback_profiles = Vec::new();
        let status = crate::hspstream::s_rps_compute_traceback(
            crate::program::RPS_BLAST,
            Some(&stream),
            Some(&query_info),
            Some(&rps_info),
            &crate::options::HitSavingOptions {
                hitlist_size: 10,
                ..Default::default()
            },
            Some(&mut results),
            |hsp_list, _gap_data, profile_index, karlin_k| {
                callback_profiles.push(profile_index);
                assert_eq!(profile_index, 1);
                assert!((karlin_k.expect("RPS K") - 0.30).abs() < 1e-12);
                assert_eq!(hsp_list.oid, 1);
                assert_eq!(hsp_list.hsp_array[0].as_ref().unwrap().score, 50);
                assert_eq!(hsp_list.hsp_array[0].as_ref().unwrap().query.offset, 3);
                assert_eq!(hsp_list.hsp_array[0].as_ref().unwrap().subject.offset, 3);
                0
            },
            None,
            None,
        );

        assert_eq!(status, 0);
        assert_eq!(callback_profiles, vec![1]);
        let hitlist = results.hitlists[0].as_ref().expect("query hitlist");
        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(hitlist.hsp_lists[0].oid, 1);
        let hsp = &hitlist.hsp_lists[0].hsps[0];
        assert_eq!(hsp.query_offset, 3);
        assert_eq!(hsp.query_end, 13);
        assert_eq!(hsp.subject_offset, 3);
        assert_eq!(hsp.subject_end, 13);
    }

    #[test]
    fn rps_profile_database_from_native_bundle_parses_profile_file_and_scans() {
        fn native_profile_bytes(
            magic: i32,
            num_profiles: i32,
            start_offsets: &[i32],
            row_values: &[i32],
            big_endian: bool,
        ) -> Vec<u8> {
            let mut values = Vec::with_capacity(2 + start_offsets.len() + row_values.len());
            values.push(magic);
            values.push(num_profiles);
            values.extend_from_slice(start_offsets);
            values.extend_from_slice(row_values);

            let mut bytes = Vec::with_capacity(values.len() * std::mem::size_of::<i32>());
            for value in values {
                if big_endian {
                    bytes.extend_from_slice(&value.to_be_bytes());
                } else {
                    bytes.extend_from_slice(&value.to_le_bytes());
                }
            }
            bytes
        }

        let offsets = [0, 8, 24];
        let mut pssm_values = vec![-8; offsets[2] as usize * 26];
        let mut query = vec![25u8; 24];
        query[3..6].copy_from_slice(&[1, 2, 3]);
        query[10..13].copy_from_slice(&[4, 5, 6]);
        for (subject_offset, residue) in query.iter().enumerate().take(13).skip(3) {
            let profile_query_offset = 8 + subject_offset as i32;
            pssm_values[profile_query_offset as usize * 26 + *residue as usize] = 5;
        }

        let bundle = NativeRpsProfileBundle {
            profile_bytes: native_profile_bytes(
                crate::hspstream::RPS_MAGIC_NUM,
                2,
                &offsets,
                &pssm_values,
                true,
            ),
            freq_ratios_bytes: None,
            lookup_bytes: None,
            aux_bytes: Some(b"BLOSUM62 11 1 0.041 0.14 16 24 100.0 8 0.11 16 0.27\n".to_vec()),
            consensus_sequences: vec![
                vec![7, 7, 7, 7, 7, 7, 7, 7],
                vec![7, 7, 7, 1, 2, 3, 7, 7, 7, 7, 4, 5, 6, 7, 7, 7],
            ],
            karlin_k: vec![0.0, 0.0],
            wordsize: 3,
        };

        let mut profile_db =
            blast_rps_profile_database_from_native_bundle(bundle).expect("native RPS bundle");
        assert_eq!(
            profile_db.traceback_info.profile_header.start_offsets,
            offsets
        );
        assert_eq!(profile_db.traceback_info.karlin_k, vec![0.11, 0.27]);
        assert_eq!(profile_db.lookup.rps_pssm.len(), offsets[2] as usize + 1);
        assert_eq!(profile_db.lookup.rps_pssm[offsets[2] as usize], vec![0; 26]);

        let stream =
            crate::hspstream::blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);
        let mut init_hitlist = crate::extend::InitHitList::new();
        let mut stats = crate::diagnostics::UngappedStats::default();
        assert_eq!(
            blast_rps_profile_database_scan_query_write_hsp_stream(
                &mut profile_db,
                &query,
                &stream,
                0,
                40,
                6,
                1,
                4,
                Some(&mut init_hitlist),
                Some(&mut stats),
            ),
            0
        );
        assert_eq!(init_hitlist.total(), 1);
        assert_eq!(stats.lookup_hits, 2);
        assert_eq!(stats.good_init_extends, 1);
    }

    #[test]
    fn rps_profile_database_from_native_files_wires_consensus_db_sequences() {
        fn native_values(values: &[i32]) -> Vec<u8> {
            values
                .iter()
                .flat_map(|value| value.to_le_bytes())
                .collect()
        }

        let scratch = std::env::temp_dir().join(format!(
            "blast-rs-rps-files-{}-{}",
            std::process::id(),
            std::thread::current().name().unwrap_or("lookup")
        ));
        let _ = std::fs::remove_dir_all(&scratch);
        std::fs::create_dir_all(&scratch).expect("scratch dir");
        let db_base = scratch.join("rpsdb");

        let profile_ascii = b"GGGACDGGGHIKGGG";
        let profile_encoded = crate::encoding::encode_ncbistdaa_sequence(profile_ascii);
        let mut builder = crate::BlastDbBuilder::new(crate::DbType::Protein, "rps consensus");
        builder.add(crate::SequenceEntry {
            title: "profile0".to_string(),
            accession: "profile0".to_string(),
            sequence: profile_ascii.to_vec(),
            taxid: None,
        });
        builder.write(&db_base).expect("write protein consensus db");

        let offsets = [0, profile_encoded.len() as i32];
        let mut pssm_values = vec![-8; profile_encoded.len() * 26];
        for (offset, &residue) in profile_encoded.iter().enumerate() {
            pssm_values[offset * 26 + residue as usize] = 5;
        }
        let mut profile_file_values =
            vec![crate::hspstream::RPS_MAGIC_NUM, 1, offsets[0], offsets[1]];
        profile_file_values.extend_from_slice(&pssm_values);
        std::fs::write(
            rps_sidecar_path(&db_base, "rps"),
            native_values(&profile_file_values),
        )
        .expect("write rps sidecar");
        std::fs::write(
            rps_sidecar_path(&db_base, "aux"),
            format!(
                "BLOSUM62 11 1 0.041 0.14 {} {} 100.0 {} 0.21\n",
                profile_encoded.len(),
                profile_encoded.len(),
                profile_encoded.len()
            ),
        )
        .expect("write aux sidecar");

        let profile_db = blast_rps_profile_database_from_native_files(&db_base, 3, vec![0.0])
            .expect("file-backed native RPS database");
        assert_eq!(
            profile_db.consensus_sequences,
            vec![profile_encoded.clone()]
        );
        assert_eq!(profile_db.traceback_info.karlin_k, vec![0.21]);

        let acd_index =
            s_compute_table_index(3, profile_db.lookup.charsize, &profile_encoded[3..6]);
        assert_eq!(
            profile_db.lookup.rps_backbone[acd_index].offset_pairs,
            vec![OffsetPair {
                query_offset: 5,
                subject_offset: 0,
            }]
        );

        let _ = std::fs::remove_dir_all(&scratch);
    }

    #[test]
    fn rps_profile_database_from_native_bundle_reuses_native_lookup_file() {
        fn native_values(values: &[i32], big_endian: bool) -> Vec<u8> {
            values
                .iter()
                .flat_map(|value| {
                    if big_endian {
                        value.to_be_bytes()
                    } else {
                        value.to_le_bytes()
                    }
                })
                .collect()
        }

        fn native_profile_bytes(
            magic: i32,
            num_profiles: i32,
            start_offsets: &[i32],
            row_values: &[i32],
            big_endian: bool,
        ) -> Vec<u8> {
            let mut values = Vec::with_capacity(2 + start_offsets.len() + row_values.len());
            values.push(magic);
            values.push(num_profiles);
            values.extend_from_slice(start_offsets);
            values.extend_from_slice(row_values);
            native_values(&values, big_endian)
        }

        fn native_lookup_bytes(
            magic: i32,
            hits: &[(usize, Vec<i32>)],
            big_endian: bool,
        ) -> Vec<u8> {
            let alphabet_size = if magic == crate::hspstream::RPS_MAGIC_NUM {
                26
            } else {
                28
            };
            let charsize = crate::util::ilog2(alphabet_size) + 1;
            let backbone_size = 1usize << (RPS_LOOKUP_WORDSIZE * charsize) as usize;
            let start_of_backbone = (RPS_LOOKUP_HEADER_WORDS * std::mem::size_of::<i32>()) as i32;
            let mut backbone = vec![[0i32; RPS_LOOKUP_BACKBONE_CELL_WORDS]; backbone_size + 1];
            let mut overflow = Vec::new();
            let mut total_hits = 0i32;
            for &(index, ref offsets) in hits {
                total_hits += offsets.len() as i32;
                backbone[index][0] = offsets.len() as i32;
                if offsets.len() <= RPS_HITS_PER_CELL {
                    for (slot, &offset) in offsets.iter().enumerate() {
                        backbone[index][slot + 1] = offset;
                    }
                } else {
                    backbone[index][1] = offsets[0];
                    backbone[index][2] = (overflow.len() * std::mem::size_of::<i32>()) as i32;
                    overflow.extend_from_slice(&offsets[1..]);
                }
            }
            let end_of_overflow = start_of_backbone
                + (backbone.len() * RPS_LOOKUP_BACKBONE_CELL_WORDS * std::mem::size_of::<i32>())
                    as i32
                + (overflow.len() * std::mem::size_of::<i32>()) as i32;
            let mut values = vec![
                magic,
                1,
                total_hits,
                hits.len() as i32,
                overflow.len() as i32,
                0,
                0,
                0,
                start_of_backbone,
                end_of_overflow,
            ];
            for cell in backbone {
                values.extend_from_slice(&cell);
            }
            values.extend_from_slice(&overflow);
            native_values(&values, big_endian)
        }

        let offsets = [0, 8, 24];
        let mut pssm_values = vec![-8; offsets[2] as usize * 26];
        let mut query = vec![25u8; 24];
        query[3..6].copy_from_slice(&[1, 2, 3]);
        query[10..13].copy_from_slice(&[4, 5, 6]);
        for (subject_offset, residue) in query.iter().enumerate().take(13).skip(3) {
            let profile_query_offset = 8 + subject_offset as i32;
            pssm_values[profile_query_offset as usize * 26 + *residue as usize] = 5;
        }

        let first_word_index = s_compute_table_index(3, 5, &[1, 2, 3]);
        let second_word_index = s_compute_table_index(3, 5, &[4, 5, 6]);
        let bundle = NativeRpsProfileBundle {
            profile_bytes: native_profile_bytes(
                crate::hspstream::RPS_MAGIC_NUM,
                2,
                &offsets,
                &pssm_values,
                true,
            ),
            freq_ratios_bytes: None,
            lookup_bytes: Some(native_lookup_bytes(
                crate::hspstream::RPS_MAGIC_NUM,
                &[(first_word_index, vec![13]), (second_word_index, vec![20])],
                true,
            )),
            aux_bytes: None,
            consensus_sequences: vec![
                vec![7, 7, 7, 7, 7, 7, 7, 7],
                vec![7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
            ],
            karlin_k: vec![0.1, 0.25],
            wordsize: 3,
        };

        let mut profile_db =
            blast_rps_profile_database_from_native_bundle(bundle).expect("native RPS bundle");
        assert_eq!(
            profile_db.lookup.rps_backbone[first_word_index].offset_pairs,
            vec![OffsetPair {
                query_offset: 13,
                subject_offset: 0
            }]
        );

        let stream =
            crate::hspstream::blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);
        let mut init_hitlist = crate::extend::InitHitList::new();
        let mut stats = crate::diagnostics::UngappedStats::default();
        assert_eq!(
            blast_rps_profile_database_scan_query_write_hsp_stream(
                &mut profile_db,
                &query,
                &stream,
                0,
                40,
                6,
                1,
                4,
                Some(&mut init_hitlist),
                Some(&mut stats),
            ),
            0
        );
        assert_eq!(stats.lookup_hits, 2);
        assert_eq!(stats.good_init_extends, 1);
    }

    #[test]
    fn rps_lookup_info_from_native_lookup_bytes_decodes_overflow_hits() {
        fn native_values(values: &[i32]) -> Vec<u8> {
            values
                .iter()
                .flat_map(|value| value.to_le_bytes())
                .collect()
        }

        let charsize = crate::util::ilog2(26) + 1;
        let backbone_size = 1usize << (RPS_LOOKUP_WORDSIZE * charsize) as usize;
        let start_of_backbone = (RPS_LOOKUP_HEADER_WORDS * std::mem::size_of::<i32>()) as i32;
        let mut backbone = vec![[0i32; RPS_LOOKUP_BACKBONE_CELL_WORDS]; backbone_size + 1];
        let index = s_compute_table_index(3, 5, &[1, 2, 3]);
        backbone[index] = [4, 13, 0, 0];
        let overflow = [17, 20, 23];
        let end_of_overflow = start_of_backbone
            + (backbone.len() * RPS_LOOKUP_BACKBONE_CELL_WORDS * std::mem::size_of::<i32>()) as i32
            + (overflow.len() * std::mem::size_of::<i32>()) as i32;
        let mut values = vec![
            crate::hspstream::RPS_MAGIC_NUM,
            1,
            4,
            1,
            overflow.len() as i32,
            0,
            0,
            0,
            start_of_backbone,
            end_of_overflow,
        ];
        for cell in backbone {
            values.extend_from_slice(&cell);
        }
        values.extend_from_slice(&overflow);

        let info =
            rps_lookup_info_from_native_lookup_bytes(&native_values(&values), vec![vec![0; 26]])
                .expect("native lookup");
        assert_eq!(
            info.rps_backbone[index].offset_pairs,
            vec![
                OffsetPair {
                    query_offset: 13,
                    subject_offset: 0,
                },
                OffsetPair {
                    query_offset: 17,
                    subject_offset: 0,
                },
                OffsetPair {
                    query_offset: 20,
                    subject_offset: 0,
                },
                OffsetPair {
                    query_offset: 23,
                    subject_offset: 0,
                },
            ]
        );
    }

    #[test]
    fn rps_profile_database_new_rejects_malformed_owned_inputs() {
        let rps_info = crate::hspstream::RpsTracebackInfo {
            profile_header: crate::hspstream::RpsProfileHeader {
                magic_number: crate::hspstream::RPS_MAGIC_NUM,
                num_profiles: 1,
                start_offsets: vec![0, 4],
                pssm_values: vec![0; 5 * 26],
            },
            freq_ratios_header: None,
            karlin_k: vec![0.1],
        };

        assert!(blast_rps_profile_database_new(rps_info.clone(), vec![vec![1, 2, 3]], 3).is_err());
        assert!(
            blast_rps_profile_database_new(rps_info.clone(), vec![vec![1, 2, 3, 26]], 3).is_err()
        );
        assert!(blast_rps_profile_database_new(rps_info, vec![vec![1, 2, 3, 4]], 0).is_err());
    }

    #[test]
    fn rps_profile_database_from_native_bundle_rejects_bad_frequency_metadata() {
        fn native_values(values: &[i32]) -> Vec<u8> {
            values
                .iter()
                .flat_map(|value| value.to_le_bytes())
                .collect()
        }

        let profile_rows = vec![0; 4 * 26];
        let profile_bytes = {
            let mut values = vec![crate::hspstream::RPS_MAGIC_NUM, 1, 0, 4];
            values.extend_from_slice(&profile_rows);
            native_values(&values)
        };
        let bad_freq_bytes = {
            let mut values = vec![crate::hspstream::RPS_MAGIC_NUM, 1, 0, 3];
            values.extend_from_slice(&vec![0; 3 * 26]);
            native_values(&values)
        };

        let bundle = NativeRpsProfileBundle {
            profile_bytes,
            freq_ratios_bytes: Some(bad_freq_bytes),
            lookup_bytes: None,
            aux_bytes: None,
            consensus_sequences: vec![vec![1, 2, 3, 4]],
            karlin_k: vec![0.1],
            wordsize: 3,
        };

        assert!(blast_rps_profile_database_from_native_bundle(bundle).is_err());
    }

    #[test]
    fn rps_scan_subject_buckets_corrected_offsets() {
        let mut lookup = BlastRpsLookupTable {
            alphabet_size: 28,
            wordsize: 3,
            charsize: 5,
            backbone_size: 4096,
            rps_backbone: vec![RpsBackboneCell::default(); 4096],
            pv: vec![0; (4096 >> crate::stat::PV_ARRAY_BTS) + 1],
            rps_pssm: vec![],
            num_buckets: 1,
            bucket_array: vec![RpsBucket {
                num_alloc: 1,
                num_used: 1,
                offset_pairs: vec![OffsetPair {
                    query_offset: 99,
                    subject_offset: 99,
                }],
            }],
        };
        let index = s_compute_table_index(3, 5, &[1, 2, 3]);
        lookup.pv[index >> crate::stat::PV_ARRAY_BTS] |=
            1u32 << (index & crate::stat::PV_ARRAY_MASK as usize);
        lookup.rps_backbone[index] = RpsBackboneCell {
            num_used: 2,
            offset_pairs: vec![
                OffsetPair {
                    query_offset: RPS_BUCKET_SIZE + 2,
                    subject_offset: 0,
                },
                OffsetPair {
                    query_offset: 7,
                    subject_offset: 0,
                },
            ],
        };

        let mut offset = 0;
        assert_eq!(
            blast_rps_scan_subject(&mut lookup, &[9, 1, 2, 3, 8], &mut offset),
            2
        );
        assert_eq!(offset, 3);
        assert_eq!(lookup.num_buckets, 2);
        assert_eq!(
            lookup.bucket_array[0].offset_pairs,
            vec![OffsetPair {
                query_offset: 5,
                subject_offset: 1,
            }]
        );
        assert_eq!(
            lookup.bucket_array[1].offset_pairs,
            vec![OffsetPair {
                query_offset: RPS_BUCKET_SIZE,
                subject_offset: 1,
            }]
        );

        let hits = blast_rps_word_finder(&mut lookup, &[9, 1, 2, 3, 8, 1, 2, 3]);
        assert_eq!(
            hits,
            vec![
                OffsetPair {
                    query_offset: 5,
                    subject_offset: 1,
                },
                OffsetPair {
                    query_offset: RPS_BUCKET_SIZE,
                    subject_offset: 1,
                },
                OffsetPair {
                    query_offset: 5,
                    subject_offset: 5,
                },
                OffsetPair {
                    query_offset: RPS_BUCKET_SIZE,
                    subject_offset: 5,
                }
            ]
        );
    }

    #[test]
    fn rps_two_hit_word_finder_keeps_same_diagonal_in_window_pairs() {
        let mut lookup = BlastRpsLookupTable {
            alphabet_size: 28,
            wordsize: 3,
            charsize: 5,
            backbone_size: 32768,
            rps_backbone: vec![RpsBackboneCell::default(); 32768],
            pv: vec![0; (32768 >> crate::stat::PV_ARRAY_BTS) + 1],
            rps_pssm: vec![],
            num_buckets: 1,
            bucket_array: vec![RpsBucket {
                num_alloc: 1,
                num_used: 0,
                offset_pairs: vec![],
            }],
        };
        let words = [
            ([1u8, 2, 3], 9),
            ([4u8, 5, 6], 16),
            ([7u8, 8, 9], 24),
            ([10u8, 11, 12], 65),
        ];
        for (word, corrected_query_offset) in words {
            let index = s_compute_table_index(3, 5, &word);
            lookup.pv[index >> crate::stat::PV_ARRAY_BTS] |=
                1u32 << (index & crate::stat::PV_ARRAY_MASK as usize);
            lookup.rps_backbone[index] = RpsBackboneCell {
                num_used: 1,
                offset_pairs: vec![OffsetPair {
                    query_offset: corrected_query_offset + lookup.wordsize - 1,
                    subject_offset: 0,
                }],
            };
        }

        let mut subject = vec![25u8; 64];
        subject[5..8].copy_from_slice(&[1, 2, 3]);
        subject[12..15].copy_from_slice(&[4, 5, 6]);
        subject[20..23].copy_from_slice(&[7, 8, 9]);
        subject[61..64].copy_from_slice(&[10, 11, 12]);

        let seeds = s_blast_rps_word_finder_two_hit(&mut lookup, &subject, 40);
        assert_eq!(
            seeds,
            vec![
                RpsTwoHitSeed {
                    first: OffsetPair {
                        query_offset: 9,
                        subject_offset: 5,
                    },
                    second: OffsetPair {
                        query_offset: 16,
                        subject_offset: 12,
                    },
                },
                RpsTwoHitSeed {
                    first: OffsetPair {
                        query_offset: 16,
                        subject_offset: 12,
                    },
                    second: OffsetPair {
                        query_offset: 24,
                        subject_offset: 20,
                    },
                },
            ]
        );
    }

    #[test]
    fn rps_two_hit_word_finder_side_channels_fill_init_list_and_stats() {
        let mut lookup = BlastRpsLookupTable {
            alphabet_size: 28,
            wordsize: 3,
            charsize: 5,
            backbone_size: 32768,
            rps_backbone: vec![RpsBackboneCell::default(); 32768],
            pv: vec![0; (32768 >> crate::stat::PV_ARRAY_BTS) + 1],
            rps_pssm: vec![],
            num_buckets: 1,
            bucket_array: vec![RpsBucket::default()],
        };
        for (word, corrected_query_offset) in
            [([1u8, 2, 3], 9), ([4u8, 5, 6], 16), ([7u8, 8, 9], 24)]
        {
            let index = s_compute_table_index(3, 5, &word);
            lookup.pv[index >> crate::stat::PV_ARRAY_BTS] |=
                1u32 << (index & crate::stat::PV_ARRAY_MASK as usize);
            lookup.rps_backbone[index] = RpsBackboneCell {
                num_used: 1,
                offset_pairs: vec![OffsetPair {
                    query_offset: corrected_query_offset + lookup.wordsize - 1,
                    subject_offset: 0,
                }],
            };
        }

        let mut subject = vec![25u8; 32];
        subject[5..8].copy_from_slice(&[1, 2, 3]);
        subject[12..15].copy_from_slice(&[4, 5, 6]);
        subject[20..23].copy_from_slice(&[7, 8, 9]);
        let mut init_hitlist = crate::extend::InitHitList::new();
        let mut stats = crate::diagnostics::UngappedStats::default();

        let scan = s_blast_rps_word_finder_two_hit_with_side_channels(
            &mut lookup,
            &subject,
            40,
            Some(&mut init_hitlist),
            Some(&mut stats),
        );

        assert_eq!(scan.total_hits, 3);
        assert_eq!(scan.hits_extended, 2);
        assert_eq!(scan.seeds.len(), 2);
        assert_eq!(init_hitlist.total(), 2);
        assert_eq!(init_hitlist.hits[0].query_offset, 16);
        assert_eq!(init_hitlist.hits[0].subject_offset, 12);
        assert_eq!(init_hitlist.hits[1].query_offset, 24);
        assert_eq!(init_hitlist.hits[1].subject_offset, 20);
        assert_eq!(stats.lookup_hits, 3);
        assert_eq!(stats.num_seqs_lookup_hits, 1);
        assert_eq!(stats.init_extends, 2);
        assert_eq!(stats.good_init_extends, 2);
        assert_eq!(stats.num_seqs_passed, 1);
    }

    #[test]
    fn rps_two_hit_word_finder_extension_payloads_fill_ungapped_data() {
        let mut lookup = BlastRpsLookupTable {
            alphabet_size: 28,
            wordsize: 3,
            charsize: 5,
            backbone_size: 32768,
            rps_backbone: vec![RpsBackboneCell::default(); 32768],
            pv: vec![0; (32768 >> crate::stat::PV_ARRAY_BTS) + 1],
            rps_pssm: vec![vec![-8; 28]; 24],
            num_buckets: 1,
            bucket_array: vec![RpsBucket::default()],
        };
        for (word, corrected_query_offset) in [([1u8, 2, 3], 4), ([4u8, 5, 6], 11)] {
            let index = s_compute_table_index(3, 5, &word);
            lookup.pv[index >> crate::stat::PV_ARRAY_BTS] |=
                1u32 << (index & crate::stat::PV_ARRAY_MASK as usize);
            lookup.rps_backbone[index] = RpsBackboneCell {
                num_used: 1,
                offset_pairs: vec![OffsetPair {
                    query_offset: corrected_query_offset + lookup.wordsize - 1,
                    subject_offset: 0,
                }],
            };
        }

        let mut subject = vec![25u8; 24];
        subject[3..6].copy_from_slice(&[1, 2, 3]);
        subject[10..13].copy_from_slice(&[4, 5, 6]);
        for (subject_offset, residue) in subject.iter().enumerate().take(13).skip(3) {
            let query_offset = subject_offset + 1;
            lookup.rps_pssm[query_offset][*residue as usize] = 5;
        }

        let mut init_hitlist = crate::extend::InitHitList::new();
        let mut stats = crate::diagnostics::UngappedStats::default();

        let scan = s_blast_rps_word_finder_two_hit_with_extension_payloads(
            &mut lookup,
            &subject,
            40,
            6,
            1,
            Some(&mut init_hitlist),
            Some(&mut stats),
        );

        assert_eq!(scan.total_hits, 2);
        assert_eq!(scan.hits_extended, 1);
        assert_eq!(init_hitlist.total(), 1);
        let init = &init_hitlist.hits[0];
        assert_eq!(init.query_offset, 11);
        assert_eq!(init.subject_offset, 10);
        let data = init.ungapped_data.as_ref().expect("RPS ungapped payload");
        assert_eq!(data.q_start, 4);
        assert_eq!(data.s_start, 3);
        assert_eq!(data.length, 10);
        assert_eq!(data.score, 50);
        assert_eq!(stats.lookup_hits, 2);
        assert_eq!(stats.init_extends, 1);
        assert_eq!(stats.good_init_extends, 1);
    }

    #[test]
    fn rps_two_hit_payload_suppresses_extra_seed_inside_right_extension() {
        let mut lookup = BlastRpsLookupTable {
            alphabet_size: 28,
            wordsize: 3,
            charsize: 5,
            backbone_size: 32768,
            rps_backbone: vec![RpsBackboneCell::default(); 32768],
            pv: vec![0; (32768 >> crate::stat::PV_ARRAY_BTS) + 1],
            rps_pssm: vec![vec![-8; 28]; 24],
            num_buckets: 1,
            bucket_array: vec![RpsBucket::default()],
        };
        for (word, corrected_query_offset) in
            [([1u8, 2, 3], 4), ([4u8, 5, 6], 11), ([7u8, 8, 9], 15)]
        {
            let index = s_compute_table_index(3, 5, &word);
            lookup.pv[index >> crate::stat::PV_ARRAY_BTS] |=
                1u32 << (index & crate::stat::PV_ARRAY_MASK as usize);
            lookup.rps_backbone[index] = RpsBackboneCell {
                num_used: 1,
                offset_pairs: vec![OffsetPair {
                    query_offset: corrected_query_offset + lookup.wordsize - 1,
                    subject_offset: 0,
                }],
            };
        }

        let mut subject = vec![25u8; 20];
        subject[3..6].copy_from_slice(&[1, 2, 3]);
        subject[10..13].copy_from_slice(&[4, 5, 6]);
        subject[13] = 10;
        subject[14..17].copy_from_slice(&[7, 8, 9]);
        for (subject_offset, residue) in subject.iter().enumerate().take(20).skip(3) {
            let query_offset = subject_offset + 1;
            lookup.rps_pssm[query_offset][*residue as usize] = 5;
        }

        let mut init_hitlist = crate::extend::InitHitList::new();
        let mut stats = crate::diagnostics::UngappedStats::default();
        let scan = s_blast_rps_word_finder_two_hit_with_extension_payloads(
            &mut lookup,
            &subject,
            40,
            100,
            1,
            Some(&mut init_hitlist),
            Some(&mut stats),
        );

        assert_eq!(scan.total_hits, 3);
        assert_eq!(scan.hits_extended, 1);
        assert_eq!(
            scan.seeds,
            vec![RpsTwoHitSeed {
                first: OffsetPair {
                    query_offset: 4,
                    subject_offset: 3,
                },
                second: OffsetPair {
                    query_offset: 11,
                    subject_offset: 10,
                },
            }]
        );
        assert_eq!(init_hitlist.total(), 1);
        assert_eq!(stats.init_extends, 1);
        assert_eq!(stats.good_init_extends, 1);
    }

    #[test]
    fn rps_word_finder_payloads_convert_to_pre_traceback_hsp_list() {
        let mut lookup = BlastRpsLookupTable {
            alphabet_size: 28,
            wordsize: 3,
            charsize: 5,
            backbone_size: 32768,
            rps_backbone: vec![RpsBackboneCell::default(); 32768],
            pv: vec![0; (32768 >> crate::stat::PV_ARRAY_BTS) + 1],
            rps_pssm: vec![vec![-8; 28]; 24],
            num_buckets: 1,
            bucket_array: vec![RpsBucket::default()],
        };
        for (word, corrected_query_offset) in [([1u8, 2, 3], 4), ([4u8, 5, 6], 11)] {
            let index = s_compute_table_index(3, 5, &word);
            lookup.pv[index >> crate::stat::PV_ARRAY_BTS] |=
                1u32 << (index & crate::stat::PV_ARRAY_MASK as usize);
            lookup.rps_backbone[index] = RpsBackboneCell {
                num_used: 1,
                offset_pairs: vec![OffsetPair {
                    query_offset: corrected_query_offset + lookup.wordsize - 1,
                    subject_offset: 0,
                }],
            };
        }

        let mut subject = vec![25u8; 24];
        subject[3..6].copy_from_slice(&[1, 2, 3]);
        subject[10..13].copy_from_slice(&[4, 5, 6]);
        for (subject_offset, residue) in subject.iter().enumerate().take(13).skip(3) {
            let query_offset = subject_offset + 1;
            lookup.rps_pssm[query_offset][*residue as usize] = 5;
        }

        let scan = s_blast_rps_word_finder_two_hit_with_extension_payloads(
            &mut lookup,
            &subject,
            40,
            6,
            1,
            None,
            None,
        );
        let hsp_list = blast_rps_word_finder_scan_to_hsp_list(&scan, 7, 4).expect("RPS HSP list");

        assert_eq!(hsp_list.oid, 7);
        assert_eq!(hsp_list.hsps.len(), 1);
        let hsp = &hsp_list.hsps[0];
        assert_eq!(hsp.score, 50);
        assert_eq!(hsp.query_offset, 4);
        assert_eq!(hsp.query_end, 14);
        assert_eq!(hsp.query_gapped_start, 11);
        assert_eq!(hsp.subject_offset, 3);
        assert_eq!(hsp.subject_end, 13);
        assert_eq!(hsp.subject_gapped_start, 10);
        assert_eq!(hsp.evalue, f64::MAX);
        assert!(hsp.edit_script.is_none());

        let mut restored = hsp_list.clone();
        crate::hspstream::s_blast_hsp_list_rps_update(crate::program::RPS_BLAST, &mut restored);
        let restored_hsp = &restored.hsps[0];
        assert_eq!(restored_hsp.query_offset, 3);
        assert_eq!(restored_hsp.subject_offset, 4);
        assert_eq!(restored_hsp.query_gapped_start, 10);
        assert_eq!(restored_hsp.subject_gapped_start, 11);
    }

    #[test]
    fn rps_word_finder_payload_cutoff_blocks_hsp_bridge() {
        let mut lookup = BlastRpsLookupTable {
            alphabet_size: 28,
            wordsize: 3,
            charsize: 5,
            backbone_size: 32768,
            rps_backbone: vec![RpsBackboneCell::default(); 32768],
            pv: vec![0; (32768 >> crate::stat::PV_ARRAY_BTS) + 1],
            rps_pssm: vec![vec![-8; 28]; 24],
            num_buckets: 1,
            bucket_array: vec![RpsBucket::default()],
        };
        for (word, corrected_query_offset) in [([1u8, 2, 3], 4), ([4u8, 5, 6], 11)] {
            let index = s_compute_table_index(3, 5, &word);
            lookup.pv[index >> crate::stat::PV_ARRAY_BTS] |=
                1u32 << (index & crate::stat::PV_ARRAY_MASK as usize);
            lookup.rps_backbone[index] = RpsBackboneCell {
                num_used: 1,
                offset_pairs: vec![OffsetPair {
                    query_offset: corrected_query_offset + lookup.wordsize - 1,
                    subject_offset: 0,
                }],
            };
        }

        let mut subject = vec![25u8; 24];
        subject[3..6].copy_from_slice(&[1, 2, 3]);
        subject[10..13].copy_from_slice(&[4, 5, 6]);
        for (subject_offset, residue) in subject.iter().enumerate().take(13).skip(3) {
            let query_offset = subject_offset + 1;
            lookup.rps_pssm[query_offset][*residue as usize] = 5;
        }
        let mut stats = crate::diagnostics::UngappedStats::default();

        let scan = s_blast_rps_word_finder_two_hit_with_extension_payloads(
            &mut lookup,
            &subject,
            40,
            6,
            51,
            None,
            Some(&mut stats),
        );

        assert_eq!(scan.total_hits, 2);
        assert_eq!(scan.hits_extended, 1);
        assert_eq!(stats.lookup_hits, 2);
        assert_eq!(stats.init_extends, 1);
        assert_eq!(stats.good_init_extends, 0);
        assert!(blast_rps_word_finder_scan_to_hsp_list(&scan, 7, 4).is_none());
    }

    #[test]
    fn rps_word_finder_empty_payload_adapter_leaves_stream_writable() {
        let stream =
            crate::hspstream::blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);
        let scan = RpsWordFinderScan::default();

        assert_eq!(
            blast_rps_word_finder_scan_write_hsp_stream(&scan, &stream, 0, 7, 4),
            0
        );
        assert_eq!(
            stream.blast_hspstream_write_blast_hsp_list(
                0,
                crate::hspstream::blast_hsp_list_from_legacy_hsp_list(
                    crate::hspstream::HspList::new(7),
                    0,
                ),
            ),
            0
        );
    }

    #[test]
    fn rps_word_finder_payload_bridge_writes_traceback_stream() {
        let scan = RpsWordFinderScan {
            seeds: Vec::new(),
            total_hits: 4,
            hits_extended: 2,
            init_hits: vec![
                RpsInitHitRecord {
                    q_start: 6,
                    s_start: 2,
                    q_off: 8,
                    s_off: 4,
                    len: 4,
                    score: 45,
                },
                RpsInitHitRecord {
                    q_start: 3,
                    s_start: 1,
                    q_off: 5,
                    s_off: 3,
                    len: 4,
                    score: 20,
                },
            ],
        };
        let stream =
            crate::hspstream::blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);

        assert_eq!(
            blast_rps_word_finder_scan_write_hsp_stream(&scan, &stream, 0, 0, 1),
            0
        );

        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[32]);
        let rps_info = crate::hspstream::RpsTracebackInfo {
            profile_header: crate::hspstream::RpsProfileHeader {
                magic_number: crate::hspstream::RPS_MAGIC_NUM,
                num_profiles: 1,
                start_offsets: vec![0, 12],
                pssm_values: (0..(13 * 26)).collect(),
            },
            freq_ratios_header: None,
            karlin_k: vec![0.25],
        };
        let hit_options = crate::options::HitSavingOptions {
            hitlist_size: 10,
            ..Default::default()
        };
        let mut results = crate::hspstream::HspResults::new(1);

        let status = crate::hspstream::s_rps_compute_traceback(
            crate::program::RPS_BLAST,
            Some(&stream),
            Some(&query_info),
            Some(&rps_info),
            &hit_options,
            Some(&mut results),
            |hsp_list, _gap_data, profile_index, karlin_k| {
                assert_eq!(profile_index, 0);
                assert!((karlin_k.expect("RPS K") - 0.30).abs() < 1e-12);
                assert_eq!(hsp_list.hspcnt as usize, 1);
                let hsp = hsp_list.hsp_array[0].as_ref().unwrap();
                assert_eq!(hsp.score, 45);
                assert_eq!(hsp.query.offset, 6);
                assert_eq!(hsp.subject.offset, 2);
                0
            },
            None,
            None,
        );

        assert_eq!(status, 0);
        let hitlist = results.hitlists[0].as_ref().expect("query hitlist");
        assert_eq!(hitlist.hsp_lists.len(), 1);
        let hsp = &hitlist.hsp_lists[0].hsps[0];
        assert_eq!(hsp.score, 45);
        assert_eq!(hsp.query_offset, 2);
        assert_eq!(hsp.query_end, 6);
        assert_eq!(hsp.query_gapped_start, 4);
        assert_eq!(hsp.subject_offset, 6);
        assert_eq!(hsp.subject_end, 10);
        assert_eq!(hsp.subject_gapped_start, 8);
    }

    #[test]
    fn rps_subject_driver_adapter_scans_payloads_into_traceback_stream() {
        let mut lookup = BlastRpsLookupTable {
            alphabet_size: 28,
            wordsize: 3,
            charsize: 5,
            backbone_size: 32768,
            rps_backbone: vec![RpsBackboneCell::default(); 32768],
            pv: vec![0; (32768 >> crate::stat::PV_ARRAY_BTS) + 1],
            rps_pssm: vec![vec![-8; 28]; 24],
            num_buckets: 1,
            bucket_array: vec![RpsBucket::default()],
        };
        for (word, corrected_query_offset) in [([1u8, 2, 3], 4), ([4u8, 5, 6], 11)] {
            let index = s_compute_table_index(3, 5, &word);
            lookup.pv[index >> crate::stat::PV_ARRAY_BTS] |=
                1u32 << (index & crate::stat::PV_ARRAY_MASK as usize);
            lookup.rps_backbone[index] = RpsBackboneCell {
                num_used: 1,
                offset_pairs: vec![OffsetPair {
                    query_offset: corrected_query_offset + lookup.wordsize - 1,
                    subject_offset: 0,
                }],
            };
        }

        let mut subject = vec![25u8; 24];
        subject[3..6].copy_from_slice(&[1, 2, 3]);
        subject[10..13].copy_from_slice(&[4, 5, 6]);
        for (subject_offset, residue) in subject.iter().enumerate().take(13).skip(3) {
            let query_offset = subject_offset + 1;
            lookup.rps_pssm[query_offset][*residue as usize] = 5;
        }

        let stream =
            crate::hspstream::blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);
        let mut init_hitlist = crate::extend::InitHitList::new();
        let mut stats = crate::diagnostics::UngappedStats::default();

        assert_eq!(
            blast_rps_subject_scan_write_hsp_stream(
                &mut lookup,
                &subject,
                &stream,
                0,
                0,
                40,
                6,
                1,
                4,
                Some(&mut init_hitlist),
                Some(&mut stats),
            ),
            0
        );
        assert_eq!(init_hitlist.total(), 1);
        assert_eq!(stats.lookup_hits, 2);
        assert_eq!(stats.init_extends, 1);
        assert_eq!(stats.good_init_extends, 1);

        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[32]);
        let rps_info = crate::hspstream::RpsTracebackInfo {
            profile_header: crate::hspstream::RpsProfileHeader {
                magic_number: crate::hspstream::RPS_MAGIC_NUM,
                num_profiles: 1,
                start_offsets: vec![0, 12],
                pssm_values: (0..(13 * 26)).collect(),
            },
            freq_ratios_header: None,
            karlin_k: vec![0.25],
        };
        let hit_options = crate::options::HitSavingOptions {
            hitlist_size: 10,
            ..Default::default()
        };
        let mut results = crate::hspstream::HspResults::new(1);

        let status = crate::hspstream::s_rps_compute_traceback(
            crate::program::RPS_BLAST,
            Some(&stream),
            Some(&query_info),
            Some(&rps_info),
            &hit_options,
            Some(&mut results),
            |hsp_list, _gap_data, profile_index, _karlin_k| {
                assert_eq!(profile_index, 0);
                assert_eq!(hsp_list.hspcnt as usize, 1);
                assert_eq!(hsp_list.hsp_array[0].as_ref().unwrap().score, 50);
                0
            },
            None,
            None,
        );

        assert_eq!(status, 0);
        let hitlist = results.hitlists[0].as_ref().expect("query hitlist");
        assert_eq!(hitlist.hsp_lists.len(), 1);
        let hsp = &hitlist.hsp_lists[0].hsps[0];
        assert_eq!(hsp.query_offset, 3);
        assert_eq!(hsp.query_end, 13);
        assert_eq!(hsp.query_gapped_start, 10);
        assert_eq!(hsp.subject_offset, 4);
        assert_eq!(hsp.subject_end, 14);
        assert_eq!(hsp.subject_gapped_start, 11);
    }

    #[test]
    fn rps_subject_driver_adapter_no_hit_scan_leaves_traceback_stream_empty() {
        let mut lookup = BlastRpsLookupTable {
            alphabet_size: 28,
            wordsize: 3,
            charsize: 5,
            backbone_size: 32768,
            rps_backbone: vec![RpsBackboneCell::default(); 32768],
            pv: vec![0; (32768 >> crate::stat::PV_ARRAY_BTS) + 1],
            rps_pssm: vec![vec![-8; 28]; 8],
            num_buckets: 1,
            bucket_array: vec![RpsBucket::default()],
        };
        let subject = vec![25u8; 16];
        let stream =
            crate::hspstream::blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);
        let mut init_hitlist = crate::extend::InitHitList::new();
        let mut stats = crate::diagnostics::UngappedStats::default();

        assert_eq!(
            blast_rps_subject_scan_write_hsp_stream(
                &mut lookup,
                &subject,
                &stream,
                0,
                0,
                40,
                6,
                1,
                4,
                Some(&mut init_hitlist),
                Some(&mut stats),
            ),
            0
        );
        assert_eq!(init_hitlist.total(), 0);
        assert_eq!(stats.lookup_hits, 0);
        assert_eq!(stats.init_extends, 0);
        assert_eq!(stats.good_init_extends, 0);

        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[16]);
        let rps_info = crate::hspstream::RpsTracebackInfo {
            profile_header: crate::hspstream::RpsProfileHeader {
                magic_number: crate::hspstream::RPS_MAGIC_NUM,
                num_profiles: 1,
                start_offsets: vec![0, 8],
                pssm_values: (0..(9 * 26)).collect(),
            },
            freq_ratios_header: None,
            karlin_k: vec![0.25],
        };
        let mut results = crate::hspstream::HspResults::new(1);
        let mut callback_called = false;

        let status = crate::hspstream::s_rps_compute_traceback(
            crate::program::RPS_BLAST,
            Some(&stream),
            Some(&query_info),
            Some(&rps_info),
            &crate::options::HitSavingOptions::default(),
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
    }

    #[test]
    fn rps_seqsrc_driver_adapter_fetches_profiles_and_skips_missing_sequences() {
        struct TestRpsSeqSrc {
            entries: Vec<Option<Vec<u8>>>,
            encodings: std::sync::Mutex<Vec<crate::seqsrc::SeqEncoding>>,
        }

        impl crate::seqsrc::BlastSeqSource for TestRpsSeqSrc {
            fn num_seqs(&self) -> i32 {
                self.entries.len() as i32
            }

            fn total_length(&self) -> i64 {
                self.entries
                    .iter()
                    .flatten()
                    .map(|seq| seq.len() as i64)
                    .sum()
            }

            fn max_seq_len(&self) -> i32 {
                self.entries
                    .iter()
                    .flatten()
                    .map(Vec::len)
                    .max()
                    .unwrap_or(0) as i32
            }

            fn avg_seq_len(&self) -> i32 {
                if self.entries.is_empty() {
                    0
                } else {
                    (self.total_length() / self.entries.len() as i64) as i32
                }
            }

            fn name(&self) -> &str {
                "rps-test"
            }

            fn is_protein(&self) -> bool {
                true
            }

            fn seq_len(&self, oid: i32) -> i32 {
                self.entries
                    .get(oid as usize)
                    .and_then(|entry| entry.as_ref())
                    .map_or(0, |seq| seq.len() as i32)
            }

            fn get_sequence(
                &self,
                arg: &crate::seqsrc::GetSeqArg,
            ) -> Option<crate::seqsrc::SeqData> {
                self.encodings.lock().unwrap().push(arg.encoding);
                self.entries
                    .get(arg.oid as usize)
                    .and_then(|entry| entry.as_ref())
                    .map(|seq| crate::seqsrc::SeqData {
                        sequence: seq.clone(),
                        length: seq.len() as i32,
                    })
            }

            fn iter_oids(&self) -> Box<dyn Iterator<Item = i32> + '_> {
                Box::new(0..self.entries.len() as i32)
            }
        }

        let mut lookup = BlastRpsLookupTable {
            alphabet_size: 28,
            wordsize: 3,
            charsize: 5,
            backbone_size: 32768,
            rps_backbone: vec![RpsBackboneCell::default(); 32768],
            pv: vec![0; (32768 >> crate::stat::PV_ARRAY_BTS) + 1],
            rps_pssm: vec![vec![-8; 28]; 24],
            num_buckets: 1,
            bucket_array: vec![RpsBucket::default()],
        };
        for (word, corrected_query_offset) in [([1u8, 2, 3], 4), ([4u8, 5, 6], 11)] {
            let index = s_compute_table_index(3, 5, &word);
            lookup.pv[index >> crate::stat::PV_ARRAY_BTS] |=
                1u32 << (index & crate::stat::PV_ARRAY_MASK as usize);
            lookup.rps_backbone[index] = RpsBackboneCell {
                num_used: 1,
                offset_pairs: vec![OffsetPair {
                    query_offset: corrected_query_offset + lookup.wordsize - 1,
                    subject_offset: 0,
                }],
            };
        }

        let mut valid_profile = vec![25u8; 24];
        valid_profile[3..6].copy_from_slice(&[1, 2, 3]);
        valid_profile[10..13].copy_from_slice(&[4, 5, 6]);
        for (subject_offset, residue) in valid_profile.iter().enumerate().take(13).skip(3) {
            let query_offset = subject_offset + 1;
            lookup.rps_pssm[query_offset][*residue as usize] = 5;
        }

        let seq_src = TestRpsSeqSrc {
            entries: vec![Some(vec![1, 2]), None, Some(valid_profile)],
            encodings: std::sync::Mutex::new(Vec::new()),
        };
        let stream =
            crate::hspstream::blast_hsp_stream_new(crate::program::UNDEFINED, 0, false, 1, None);
        let mut init_hitlist = crate::extend::InitHitList::new();
        let mut stats = crate::diagnostics::UngappedStats::default();

        assert_eq!(
            blast_rps_seqsrc_scan_write_hsp_stream(
                &mut lookup,
                Some(&seq_src),
                &stream,
                crate::program::RPS_BLAST,
                0,
                40,
                6,
                1,
                4,
                Some(&mut init_hitlist),
                Some(&mut stats),
            ),
            0
        );
        assert_eq!(init_hitlist.total(), 1);
        assert_eq!(stats.lookup_hits, 2);
        assert_eq!(stats.init_extends, 1);
        assert_eq!(stats.good_init_extends, 1);
        assert_eq!(
            *seq_src.encodings.lock().unwrap(),
            vec![
                crate::seqsrc::SeqEncoding::Protein,
                crate::seqsrc::SeqEncoding::Protein,
                crate::seqsrc::SeqEncoding::Protein,
            ]
        );

        let query_info = crate::queryinfo::QueryInfo::new_blastp(&[32]);
        let rps_info = crate::hspstream::RpsTracebackInfo {
            profile_header: crate::hspstream::RpsProfileHeader {
                magic_number: crate::hspstream::RPS_MAGIC_NUM,
                num_profiles: 3,
                start_offsets: vec![0, 0, 0, 12],
                pssm_values: (0..(13 * 26)).collect(),
            },
            freq_ratios_header: None,
            karlin_k: vec![0.1, 0.2, 0.3],
        };
        let mut results = crate::hspstream::HspResults::new(1);
        let mut callback_profiles = Vec::new();

        let status = crate::hspstream::s_rps_compute_traceback(
            crate::program::RPS_BLAST,
            Some(&stream),
            Some(&query_info),
            Some(&rps_info),
            &crate::options::HitSavingOptions {
                hitlist_size: 10,
                ..Default::default()
            },
            Some(&mut results),
            |hsp_list, _gap_data, profile_index, karlin_k| {
                callback_profiles.push(profile_index);
                assert_eq!(profile_index, 2);
                assert!((karlin_k.expect("RPS K") - 0.36).abs() < 1e-12);
                assert_eq!(hsp_list.oid, 2);
                assert_eq!(hsp_list.hsp_array[0].as_ref().unwrap().score, 50);
                0
            },
            None,
            None,
        );

        assert_eq!(status, 0);
        assert_eq!(callback_profiles, vec![2]);
        let hitlist = results.hitlists[0].as_ref().expect("query hitlist");
        assert_eq!(hitlist.hsp_lists.len(), 1);
        assert_eq!(hitlist.hsp_lists[0].oid, 2);
        assert_eq!(hitlist.hsp_lists[0].hsps[0].query_offset, 3);
        assert_eq!(
            blast_rps_seqsrc_scan_write_hsp_stream(
                &mut lookup,
                None,
                &stream,
                crate::program::RPS_BLAST,
                0,
                40,
                6,
                1,
                4,
                None,
                None,
            ),
            -1
        );
    }

    #[test]
    fn test_lookup_table_wrap_variants() {
        let small = LookupTableWrap::SmallNa(SmallNaLookupTable {
            word_length: 11,
            backbone: vec![],
            overflow: vec![],
            pv_array: vec![],
            longest_chain: 0,
            scan_step: 4,
        });
        assert_eq!(small.table_type(), LookupTableType::SmallNaLookup);
        assert_eq!(small.word_length(), 11);

        let aa = LookupTableWrap::Aa(AaLookupTable {
            word_length: 3,
            threshold: 11.0,
            backbone: vec![],
            pv_array: vec![],
        });
        assert_eq!(aa.table_type(), LookupTableType::AaLookup);
        assert_eq!(aa.word_length(), 3);

        let rps = LookupTableWrap::Rps(BlastRpsLookupTable {
            wordsize: 4,
            ..BlastRpsLookupTable::default()
        });
        assert_eq!(rps.table_type(), LookupTableType::RpsLookup);
        assert_eq!(rps.word_length(), 4);

        let na = LookupTableWrap::Na(BlastNaLookupTable {
            word_length: 9,
            ..Default::default()
        });
        assert_eq!(na.table_type(), LookupTableType::NaLookup);
        assert_eq!(na.word_length(), 9);

        let na_hash = LookupTableWrap::NaHash(BlastNaHashLookupTable {
            word_length: 16,
            ..Default::default()
        });
        assert_eq!(na_hash.table_type(), LookupTableType::NaHashLookup);
        assert_eq!(na_hash.word_length(), 16);

        assert_eq!(s_get_minimum_subj_seq_len(Some(&small)), 11);
        assert_eq!(s_get_minimum_subj_seq_len(Some(&na)), 9);
        assert_eq!(s_get_minimum_subj_seq_len(Some(&na_hash)), 16);
        assert_eq!(s_get_minimum_subj_seq_len(Some(&aa)), 3);
        assert_eq!(s_get_minimum_subj_seq_len(Some(&rps)), 4);
        assert_eq!(s_get_minimum_subj_seq_len(None), 0);
    }

    #[test]
    fn translated_mask_at_hash_helper_matches_c_sources() {
        assert!(!s_has_mask_at_hash_enabled(false, None));
        assert!(s_has_mask_at_hash_enabled(true, None));
        assert!(s_has_mask_at_hash_enabled(false, Some("dust m")));
        assert!(!s_has_mask_at_hash_enabled(false, Some("dust")));
    }

    #[test]
    fn translated_disc_template_type_matches_ncbi_enum_values() {
        assert_eq!(
            s_get_disc_template_type(11, 16, DiscWordType::Coding),
            DiscTemplateType::Template11_16Coding
        );
        assert_eq!(
            s_get_disc_template_type(11, 16, DiscWordType::TwoTemplates),
            DiscTemplateType::Template11_16Coding
        );
        assert_eq!(
            s_get_disc_template_type(12, 21, DiscWordType::Optimal),
            DiscTemplateType::Template12_21Optimal
        );
        assert_eq!(
            s_get_disc_template_type(10, 16, DiscWordType::Coding),
            DiscTemplateType::Contiguous
        );
        assert_eq!(DiscTemplateType::Template12_21Optimal as i32, 12);
    }

    #[test]
    fn translated_disc_mb_table_fill_uses_discontiguous_index() {
        let query = vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
        let mut accum = 0u64;
        for &base in &query {
            accum = (accum << BITS_PER_NUC) | u64::from(base);
        }
        let ecode = compute_discontiguous_index(accum, DiscTemplateType::Template11_16Coding);
        assert_eq!(
            ecode,
            ((accum as u32 & 0x00000003)
                | ((accum as u32 & 0x000000f0) >> 2)
                | ((accum as u32 & 0x00003c00) >> 4)
                | ((accum as u32 & 0x000f0000) >> 6)
                | ((accum as u32 & 0x03c00000) >> 8)
                | ((accum as u32 & 0xf0000000) >> 10)) as i32
        );

        let mut table = MbLookupTable {
            word_length: 16,
            lut_word_length: 11,
            discontiguous: false,
            template_length: 0,
            template_type: DiscTemplateType::Contiguous,
            two_templates: false,
            second_template_type: DiscTemplateType::Contiguous,
            hashtable: vec![0; 1 << 22],
            hashtable2: Vec::new(),
            next_pos: Vec::new(),
            next_pos2: Vec::new(),
            pv_array: vec![0; (1 << 22) / 32],
            pv_array_bts: crate::stat::PV_ARRAY_BTS as i32,
            longest_chain: 0,
            scan_step: 0,
        };
        let options = crate::options::LookupTableOptions {
            word_size: 11,
            mb_template_length: 16,
            mb_template_type: DiscWordType::Coding as i32,
            ..crate::options::LookupTableOptions::default()
        };

        assert_eq!(
            s_fill_disc_mb_table(
                &query,
                &[crate::util::SSeqRange { left: 0, right: 15 }],
                &mut table,
                &options,
            ),
            0
        );
        assert_eq!(table.hashtable[ecode as usize], 1);
        assert_eq!(table.next_pos[1], 0);
        assert!(test_pv_bit(
            &table.pv_array,
            ecode as u64,
            crate::stat::PV_ARRAY_BTS
        ));
        assert_eq!(table.longest_chain, 3);
    }

    #[test]
    fn translated_disc_mb_two_template_longest_chain_uses_separate_c_helpers() {
        let query = vec![0; 18];
        let ranges = [crate::util::SSeqRange {
            left: 0,
            right: query.len() as i32 - 1,
        }];
        let options = crate::options::LookupTableOptions {
            word_size: 11,
            mb_template_length: 16,
            mb_template_type: DiscWordType::TwoTemplates as i32,
            ..Default::default()
        };
        let mut table = MbLookupTable {
            word_length: 16,
            lut_word_length: 11,
            discontiguous: false,
            template_length: 0,
            template_type: DiscTemplateType::Contiguous,
            two_templates: false,
            second_template_type: DiscTemplateType::Contiguous,
            hashtable: vec![0; 1 << 22],
            hashtable2: Vec::new(),
            next_pos: Vec::new(),
            next_pos2: Vec::new(),
            pv_array: vec![0; (1 << 22) / 32],
            pv_array_bts: crate::stat::PV_ARRAY_BTS as i32,
            longest_chain: 0,
            scan_step: 0,
        };

        assert_eq!(
            s_fill_disc_mb_table(&query, &ranges, &mut table, &options),
            0
        );
        assert_eq!(table.longest_chain, 6);
    }

    #[test]
    fn translated_disc_mb_subject_scans_retrieve_coding_template_hits() {
        for (template_length, scan) in [
            (
                18,
                s_mb_disc_word_scan_subject_11_18_1
                    as fn(&MbLookupTable, &[u8], usize, usize, usize) -> Vec<OffsetPair>,
            ),
            (21, s_mb_disc_word_scan_subject_11_21_1),
        ] {
            let query: Vec<u8> = (0..template_length).map(|i| (i % 4) as u8).collect();
            let mut table = MbLookupTable {
                word_length: template_length as i32,
                lut_word_length: 11,
                discontiguous: false,
                template_length: 0,
                template_type: DiscTemplateType::Contiguous,
                two_templates: false,
                second_template_type: DiscTemplateType::Contiguous,
                hashtable: vec![0; 1 << 22],
                hashtable2: Vec::new(),
                next_pos: Vec::new(),
                next_pos2: Vec::new(),
                pv_array: vec![0; (1 << 22) / 32],
                pv_array_bts: crate::stat::PV_ARRAY_BTS as i32,
                longest_chain: 0,
                scan_step: 0,
            };
            let options = crate::options::LookupTableOptions {
                word_size: 11,
                mb_template_length: template_length as i32,
                mb_template_type: DiscWordType::Coding as i32,
                ..crate::options::LookupTableOptions::default()
            };
            assert_eq!(
                s_fill_disc_mb_table(
                    &query,
                    &[crate::util::SSeqRange {
                        left: 0,
                        right: template_length as i32 - 1,
                    }],
                    &mut table,
                    &options,
                ),
                0
            );

            let packed_subject = crate::encoding::pack_ncbi2na_bases(&query);
            assert_eq!(
                scan(&table, &packed_subject, query.len(), 10, query.len()),
                vec![OffsetPair {
                    query_offset: 0,
                    subject_offset: 0,
                }]
            );
            assert!(scan(&table, &packed_subject, query.len(), 0, query.len()).is_empty());
        }
    }

    #[test]
    fn translated_disc_mb_subject_scan_retrieves_two_template_secondary_diagonals() {
        let query = crate::encoding::encode_blastna_sequence(
            b"ACGTACGTGACTTACCGTACGTACGTGACTTACCGTACGTACGTGACTTACCGTACGT",
        );
        let subject = crate::encoding::encode_blastna_sequence(
            b"GGGGACGTACGTGACTTACCGTACGTACGTGACTTACCGTACGTACGTGACTTACCGTACGTCCCC",
        );
        let packed_subject = crate::encoding::pack_ncbi2na_bases(&subject);
        let mut table = MbLookupTable {
            word_length: 16,
            lut_word_length: 11,
            discontiguous: false,
            template_length: 0,
            template_type: DiscTemplateType::Contiguous,
            two_templates: false,
            second_template_type: DiscTemplateType::Contiguous,
            hashtable: vec![0; 1 << 22],
            hashtable2: Vec::new(),
            next_pos: Vec::new(),
            next_pos2: Vec::new(),
            pv_array: vec![0; (1 << 22) / 32],
            pv_array_bts: crate::stat::PV_ARRAY_BTS as i32,
            longest_chain: 0,
            scan_step: 0,
        };
        let options = crate::options::LookupTableOptions {
            word_size: 11,
            mb_template_length: 16,
            mb_template_type: DiscWordType::TwoTemplates as i32,
            ..crate::options::LookupTableOptions::default()
        };
        assert_eq!(
            s_fill_disc_mb_table(
                &query,
                &[crate::util::SSeqRange {
                    left: 0,
                    right: query.len() as i32 - 1,
                }],
                &mut table,
                &options,
            ),
            0
        );

        let hits =
            s_mb_disc_word_scan_subject(&table, &packed_subject, subject.len(), 0, subject.len());
        assert!(
            hits.iter()
                .any(|pair| pair.query_offset == 0 && pair.subject_offset == 4),
            "primary full-length diagonal seed missing from two-template scan: {hits:?}"
        );
        assert!(
            hits.iter()
                .any(|pair| pair.query_offset == 18 && pair.subject_offset == 4),
            "secondary left-overlap diagonal seed missing from two-template scan: {hits:?}"
        );
        assert!(
            hits.iter()
                .any(|pair| pair.query_offset == 0 && pair.subject_offset == 22),
            "secondary right-overlap diagonal seed missing from two-template scan: {hits:?}"
        );
    }

    #[test]
    fn translated_lookup_destructors_drop_owned_tables() {
        let mut small = Some(SmallNaLookupTable {
            word_length: 8,
            backbone: vec![1, 2, 3],
            overflow: vec![4, 5],
            pv_array: vec![6],
            longest_chain: 2,
            scan_step: 4,
        });
        assert!(blast_small_na_lookup_table_destruct(&mut small).is_none());
        assert!(small.is_none());

        let mut mb = Some(MbLookupTable {
            word_length: 28,
            lut_word_length: 12,
            discontiguous: false,
            template_length: 0,
            template_type: DiscTemplateType::Contiguous,
            two_templates: false,
            second_template_type: DiscTemplateType::Contiguous,
            hashtable: vec![1],
            hashtable2: vec![0],
            next_pos: vec![2],
            next_pos2: vec![0],
            pv_array: vec![3],
            pv_array_bts: crate::stat::PV_ARRAY_BTS as i32,
            longest_chain: 1,
            scan_step: 1,
        });
        assert!(blast_mb_lookup_table_destruct(&mut mb).is_none());
        assert!(mb.is_none());
    }

    #[test]
    fn translated_lookup_wrap_free_and_offset_array_size_match_c_shape() {
        let small = LookupTableWrap::SmallNa(SmallNaLookupTable {
            word_length: 8,
            backbone: vec![],
            overflow: vec![],
            pv_array: vec![],
            longest_chain: 17,
            scan_step: 4,
        });
        assert_eq!(get_offset_array_size(&small), OFFSET_ARRAY_SIZE + 17);

        let na = LookupTableWrap::Na(BlastNaLookupTable {
            longest_chain: 5,
            ..Default::default()
        });
        assert_eq!(get_offset_array_size(&na), OFFSET_ARRAY_SIZE + 5);

        let na_hash = LookupTableWrap::NaHash(BlastNaHashLookupTable {
            longest_chain: 9,
            ..Default::default()
        });
        assert_eq!(get_offset_array_size(&na_hash), OFFSET_ARRAY_SIZE + 9);

        let aa = LookupTableWrap::Aa(AaLookupTable {
            word_length: 3,
            threshold: 11.0,
            backbone: vec![vec![1, 2], vec![3, 4, 5]],
            pv_array: vec![],
        });
        assert_eq!(get_offset_array_size(&aa), OFFSET_ARRAY_SIZE + 3);

        let mut wrapped = Some(small);
        assert!(lookup_table_wrap_free(&mut wrapped).is_none());
        assert!(wrapped.is_none());
    }

    #[test]
    fn translated_lookup_table_wrap_init_replaces_owned_table() {
        let mut wrapped = Some(LookupTableWrap::Megablast(MbLookupTable {
            word_length: 28,
            lut_word_length: 12,
            discontiguous: false,
            template_length: 0,
            template_type: DiscTemplateType::Contiguous,
            two_templates: false,
            second_template_type: DiscTemplateType::Contiguous,
            hashtable: vec![1],
            hashtable2: vec![0],
            next_pos: vec![2],
            next_pos2: vec![0],
            pv_array: vec![3],
            pv_array_bts: crate::stat::PV_ARRAY_BTS as i32,
            longest_chain: 1,
            scan_step: 1,
        }));

        let replacement = LookupTableWrap::SmallNa(SmallNaLookupTable {
            word_length: 11,
            backbone: vec![0; 16],
            overflow: vec![],
            pv_array: vec![],
            longest_chain: 4,
            scan_step: 4,
        });
        assert_eq!(lookup_table_wrap_init(&mut wrapped, replacement), 0);
        let wrapped_ref = wrapped.as_ref().expect("initialized lookup table");
        assert_eq!(wrapped_ref.table_type(), LookupTableType::SmallNaLookup);
        assert_eq!(wrapped_ref.word_length(), 11);
        assert_eq!(get_offset_array_size(wrapped_ref), OFFSET_ARRAY_SIZE + 4);

        let replacement = LookupTableWrap::Aa(AaLookupTable {
            word_length: 3,
            threshold: 11.0,
            backbone: vec![vec![0]],
            pv_array: vec![],
        });
        assert_eq!(lookup_table_wrap_init_mt(&mut wrapped, replacement, 4), 0);
        assert_eq!(
            wrapped.as_ref().expect("mt initialized table").table_type(),
            LookupTableType::AaLookup
        );
    }

    #[test]
    fn lookup_table_wrap_init_nucleotide_builds_chosen_variant() {
        let query = crate::encoding::encode_blastna_sequence(b"ACGTACGT");
        let ranges = [crate::util::SSeqRange {
            left: 0,
            right: query.len() as i32 - 1,
        }];
        let options = crate::options::LookupTableOptions {
            lut_type: LookupTableType::NaLookup,
            word_size: 4,
            program_number: crate::program::BLASTN,
            ..Default::default()
        };
        let mut wrapped = None;

        assert_eq!(
            lookup_table_wrap_init_nucleotide(&query, &ranges, &options, &mut wrapped, None, 1),
            0
        );
        let wrapped_ref = wrapped.as_ref().expect("nucleotide lookup");
        assert_eq!(wrapped_ref.table_type(), LookupTableType::SmallNaLookup);
        let index = s_compute_table_index(4, BITS_PER_NUC as i32, &query[0..4]) as i32;
        assert!(s_small_na_lookup(wrapped.as_ref(), index, 0));
    }

    #[test]
    fn lookup_table_wrap_init_nucleotide_falls_back_to_standard_na() {
        let query = crate::encoding::encode_blastna_sequence(b"ACGTACGT");
        let ranges = [crate::util::SSeqRange {
            left: 0,
            right: 32_768,
        }];
        let options = crate::options::LookupTableOptions {
            lut_type: LookupTableType::NaLookup,
            word_size: 8,
            program_number: crate::program::BLASTN,
            ..Default::default()
        };
        let mut wrapped = None;

        assert_eq!(
            lookup_table_wrap_init_nucleotide(&query, &ranges, &options, &mut wrapped, None, 1),
            0
        );
        let wrapped_ref = wrapped.as_ref().expect("standard nucleotide lookup");
        assert_eq!(wrapped_ref.table_type(), LookupTableType::NaLookup);
        assert_eq!(wrapped_ref.word_length(), 8);
    }

    fn identity_hash(bytes: &[u8], mask: u32) -> u32 {
        let mut raw = [0u8; 4];
        raw.copy_from_slice(&bytes[..4]);
        u32::from_ne_bytes(raw) & mask
    }

    #[test]
    fn translated_backbone_cell_lifecycle_matches_c_shape() {
        assert_eq!(backbone_cell_init(None, 1, 2), -1);

        let mut cell = backbone_cell_new(42, 7).expect("cell");
        assert_eq!(cell.word, 42);
        assert_eq!(cell.offset, 7);
        assert_eq!(cell.num_offsets, 1);

        cell.next = backbone_cell_new(43, 8);
        let mut owned = Some(cell);
        assert!(backbone_cell_free(&mut owned).is_none());
        assert!(owned.is_none());
    }

    #[test]
    fn translated_add_word_hit_chains_duplicate_offsets_and_filters_pv() {
        let mut backbone = vec![BackboneCell::default(); 16];
        let mut offsets = vec![0; 16];

        assert_eq!(
            s_add_word_hit(
                &mut backbone,
                &mut offsets,
                2,
                2,
                &[1, 2],
                3,
                identity_hash,
                15,
                None,
            ),
            0
        );
        assert_eq!(backbone[6].word, 6);
        assert_eq!(backbone[6].offset, 4);
        assert_eq!(backbone[6].num_offsets, 1);

        assert_eq!(
            s_add_word_hit(
                &mut backbone,
                &mut offsets,
                2,
                2,
                &[1, 2],
                5,
                identity_hash,
                15,
                None,
            ),
            0
        );
        assert_eq!(backbone[6].offset, 6);
        assert_eq!(backbone[6].num_offsets, 2);
        assert_eq!(offsets[6], 4);

        let pv = [0u32; 1];
        assert_eq!(
            s_add_word_hit(
                &mut backbone,
                &mut offsets,
                2,
                2,
                &[2, 1],
                7,
                identity_hash,
                15,
                Some(&pv),
            ),
            0
        );
        assert_eq!(backbone[9].num_offsets, 0);
    }

    #[test]
    fn translated_hash_lookup_index_query_exact_matches_builds_backbone_chains() {
        let query = vec![0, 1, 2, 3, 0, 1, 2, 3];
        let mut backbone = vec![BackboneCell::default(); 64];
        let mut offsets = vec![0; query.len() + 2];

        assert_eq!(
            blast_hash_lookup_index_query_exact_matches(
                &mut backbone,
                &mut offsets,
                4,
                BITS_PER_NUC as i32,
                3,
                &query,
                &[crate::util::SSeqRange { left: 0, right: 7 }],
                identity_hash,
                63,
                None,
            ),
            0
        );

        let duplicate_word = 0b0001_1011u32;
        let duplicate_index = identity_hash(&duplicate_word.to_ne_bytes(), 63);
        let duplicate_cell = &backbone[duplicate_index as usize];
        assert_eq!(duplicate_cell.word, duplicate_word);
        assert_eq!(duplicate_cell.offset, 6);
        assert_eq!(duplicate_cell.num_offsets, 2);
        assert_eq!(offsets[6], 2);

        let leading_word = 0b0001_10u32;
        let leading_index = identity_hash(&leading_word.to_ne_bytes(), 63);
        assert_eq!(backbone[leading_index as usize].word, leading_word);
        assert_eq!(
            backbone[leading_index as usize].num_offsets, 2,
            "C indexes leading lookup-width windows even when word_length > lut_word_length"
        );
        assert_eq!(backbone[leading_index as usize].offset, 5);
        assert_eq!(offsets[5], 1);
    }

    #[test]
    fn test_pv_array_bit_operations() {
        // Simulate PV array: set bit for index 100, check it
        let mut pv = [0u32; 8]; // 256 bits
        let idx = 100usize;
        let word = idx / 32;
        let bit = idx % 32;
        pv[word] |= 1 << bit;
        assert_ne!(pv[word] & (1 << bit), 0, "PV bit should be set");
        assert_eq!(pv[3] & (1 << 4), 1 << 4, "Bit 100 = word 3, bit 4");
        // Unset bit should be 0
        assert_eq!(pv[0] & (1 << 0), 0, "Bit 0 should not be set");
    }
}
