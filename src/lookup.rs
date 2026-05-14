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

/// Generic wrapper around different lookup table types.
pub enum LookupTableWrap {
    SmallNa(SmallNaLookupTable),
    Megablast(MbLookupTable),
    Aa(AaLookupTable),
    Rps(BlastRpsLookupTable),
}

impl LookupTableWrap {
    pub fn table_type(&self) -> LookupTableType {
        match self {
            LookupTableWrap::SmallNa(_) => LookupTableType::SmallNaLookup,
            LookupTableWrap::Megablast(_) => LookupTableType::MegablastLookup,
            LookupTableWrap::Aa(_) => LookupTableType::AaLookup,
            LookupTableWrap::Rps(_) => LookupTableType::RpsLookup,
        }
    }

    pub fn word_length(&self) -> i32 {
        match self {
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

fn blast_mb_lookup_retrieve_from_chain(
    hashtable: &[i32],
    next_pos: &[i32],
    index: i64,
    offset_pairs: &mut Vec<OffsetPair>,
    s_off: i32,
) -> i32 {
    if index < 0 {
        return 0;
    }
    let mut q_off = hashtable.get(index as usize).copied().unwrap_or(0);
    let mut count = 0;
    while q_off != 0 {
        offset_pairs.push(OffsetPair {
            query_offset: q_off - 1,
            subject_offset: s_off,
        });
        count += 1;
        q_off = next_pos.get(q_off as usize).copied().unwrap_or(0);
    }
    count
}

/// Port of NCBI inline `s_BlastMBLookupRetrieve` (`blast_nascan.c:1413`).
pub fn s_blast_mb_lookup_retrieve(
    lookup: &MbLookupTable,
    index: i64,
    offset_pairs: &mut Vec<OffsetPair>,
    s_off: i32,
) -> i32 {
    blast_mb_lookup_retrieve_from_chain(
        &lookup.hashtable,
        &lookup.next_pos,
        index,
        offset_pairs,
        s_off,
    )
}

/// Port of NCBI inline `s_BlastMBLookupRetrieve2` (`blast_nascan.c:1437`).
pub fn s_blast_mb_lookup_retrieve2(
    lookup: &MbLookupTable,
    index: i64,
    offset_pairs: &mut Vec<OffsetPair>,
    s_off: i32,
) -> i32 {
    blast_mb_lookup_retrieve_from_chain(
        &lookup.hashtable2,
        &lookup.next_pos2,
        index,
        offset_pairs,
        s_off,
    )
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
            let lut_word_length = lut.word_length;
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
        NaLookupCallback::Na => return None,
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

fn compute_table_index(wordsize: i32, charsize: i32, seq: &[u8]) -> usize {
    let mut index = 0usize;
    for &base in seq.iter().take(wordsize.max(0) as usize) {
        index = (index << charsize.max(0) as usize) | base as usize;
    }
    index
}

fn blast_lookup_add_word_hit_to_thin(
    backbone: &mut [Option<Vec<i32>>],
    wordsize: i32,
    charsize: i32,
    seq: &[u8],
    query_offset: i32,
) {
    let index = compute_table_index(wordsize, charsize, seq);
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
            blast_lookup_add_word_hit_to_thin(
                backbone,
                lut_word_length,
                charsize,
                word,
                start as i32,
            );
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

        let first_indexed_start = from + word_length_usize - lut_word_length_usize;
        let last_start = to + 1 - lut_word_length_usize;
        for start in first_indexed_start..=last_start {
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

/// Port of NCBI `s_NaHashLookupRemovePolyAWords` (`blast_nalookup.c:1965`).
pub fn s_na_hash_lookup_remove_poly_a_words(lookup: Option<&mut BlastNaHashLookupTable>) -> i16 {
    let Some(lookup) = lookup else {
        return -1;
    };
    let pv_array_bts = lookup.pv_array_bts.max(0) as usize;
    let word_size = lookup.lut_word_length.max(0);

    clear_pv_bit(&mut lookup.pv, 0, pv_array_bts);
    clear_pv_bit(&mut lookup.pv, u32::MAX, pv_array_bts);

    for i in 1..4u32 {
        let word = i;
        for _ in 0..word_size {
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
                    start += scan_step;
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
                if s_insert_mb_word(mb_lt, &mut helper_array, ecode2, query_index, true) != 0 {
                    return -1;
                }
            }
        }
    }

    mb_lt.longest_chain = helper_array.into_iter().max().unwrap_or(2).max(2) as i32;
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
    s_mb_disc_word_scan_subject_template(
        mb_lt,
        subject,
        subject_len,
        max_hits,
        scan_range,
        18,
        DiscTemplateType::Template11_18Coding,
    )
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
    s_mb_disc_word_scan_subject_template(
        mb_lt,
        subject,
        subject_len,
        max_hits,
        scan_range,
        21,
        DiscTemplateType::Template11_21Coding,
    )
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

    let mut alphabet_size = info.alphabet_size;
    if alphabet_size == 26 {
        alphabet_size = 28;
    }
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

fn rps_compute_table_index(wordsize: i32, charsize: i32, seq: &[u8]) -> usize {
    let mut index = 0usize;
    for &letter in seq.iter().take(wordsize.max(0) as usize) {
        index = (index << charsize.max(0) as usize) | letter as usize;
    }
    index
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
        let index = rps_compute_table_index(
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

pub const OFFSET_ARRAY_SIZE: i32 = 4096;

/// 1-1 ownership translation of `LookupTableWrapFree`.
pub fn lookup_table_wrap_free(lookup: &mut Option<LookupTableWrap>) -> Option<LookupTableWrap> {
    let Some(table) = lookup.take() else {
        return None;
    };

    match table {
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

/// Rust ownership-shaped equivalent of NCBI `LookupTableWrapInit_MT`
/// (`lookup_wrap.c:61`) once the concrete lookup table has been built.
///
/// C constructs the table inside this function after switching on
/// `lookup_options->lut_type`. In Rust, callers already construct the typed
/// table because the construction inputs live in program-specific modules. Keep
/// the represented switch here so the audited name still exposes the same
/// ownership and variant boundary as C. `num_threads` only affects C's
/// `eNaHashLookupTable` construction path, which is not represented by
/// `LookupTableWrap` today.
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
        let mut thin = vec![BackboneCell::default(); 3];
        thin[1] = BackboneCell {
            word: 5,
            offset: 1,
            num_offsets: 2,
            next: None,
        };
        thin[2] = BackboneCell {
            word: 9,
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

        assert_eq!(lookup.backbone_size, 3);
        assert_eq!(lookup.longest_chain, 10);
        assert_eq!(lookup.thick_backbone[1].num_words, 1);
        assert_eq!(lookup.thick_backbone[1].num_offsets, [2, 0, 0]);
        assert_eq!(lookup.thick_backbone[1].words, [5, 0, 0]);
        assert_eq!(lookup.thick_backbone[1].offsets[..2], [0, 1]);
        assert_eq!(lookup.thick_backbone[2].num_words, 2);
        assert_eq!(lookup.thick_backbone[2].words, [9, 10, 0]);
        assert_eq!(lookup.thick_backbone[2].offsets[0], 0);
        assert_eq!(
            lookup.overflow,
            vec![9, 4, 2, 3, 4, 5, 10, 6, 6, 7, 8, 9, 10, 11]
        );
        assert_eq!(lookup.offsets_size, lookup.overflow.len() as i32);
        assert!(lookup.pv[0] & (1 << 5) != 0);
        assert!(lookup.pv[0] & (1 << 9) != 0);
        assert!(lookup.pv[0] & (1 << 10) != 0);
        assert!(thin[2].next.is_none());
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
            lut_word_length: 4,
            pv_array_bts: crate::stat::PV_ARRAY_BTS as i32,
            pv: vec![u32::MAX; 1],
            ..Default::default()
        };

        assert_eq!(s_na_hash_lookup_remove_poly_a_words(Some(&mut lookup)), 0);
        for index in [0u32, 1, 2, 3] {
            assert_eq!(
                lookup.pv[(index >> crate::stat::PV_ARRAY_BTS) as usize]
                    & (1 << (index & crate::stat::PV_ARRAY_MASK)),
                0
            );
        }
        assert_eq!(s_na_hash_lookup_remove_poly_a_words(None), -1);
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
        assert_eq!(lookup_ref.alphabet_size, 28);
        assert_eq!(lookup_ref.wordsize, 3);
        assert_eq!(lookup_ref.charsize, crate::util::ilog2(28) + 1);
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
        let index = rps_compute_table_index(3, 5, &[1, 2, 3]);
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

        assert_eq!(s_get_minimum_subj_seq_len(Some(&small)), 11);
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
        assert_eq!(table.longest_chain, 2);
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
