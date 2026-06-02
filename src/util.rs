//! Rust equivalent of blast_util.c — BLAST utility functions.

use std::sync::{Mutex, OnceLock};

/// Translate a codon using a genetic code (NCBI2na fast path).
/// Takes 3 NCBI2na bases (2-bit unambiguous A=0/C=1/G=2/T=3) and
/// returns the amino acid in NCBIstdaa encoding.
///
/// For NCBI4na input with ambiguity resolution, use [`s_codon_to_aa`] —
/// the 1-1 port of NCBI's `s_CodonToAA`.
pub fn translate_codon(b1: u8, b2: u8, b3: u8, genetic_code: &[u8; 64]) -> u8 {
    let idx = ((b1 & 3) as usize) * 16 + ((b2 & 3) as usize) * 4 + (b3 & 3) as usize;
    genetic_code[idx]
}

/// Translate a codon given in NCBI4na (each byte is a 4-bit mask over
/// {A=1, C=2, G=4, T=8}) into an NCBIstdaa amino acid.
///
/// Port of NCBI `s_CodonToAA` (`blast_util.c:369`). When an input base
/// is an ambiguity code matching multiple unambiguous bases, enumerates
/// every compatible codon; returns the unique amino acid if they all
/// agree, otherwise returns NCBIstdaa `X`. `FENCE_SENTRY` is preserved.
///
/// NCBI's C iterates `mapping[4] = {T=8, C=2, A=1, G=4}` so its
/// `codes` table is TCAG-ordered. Our [`STANDARD_GENETIC_CODE`] is
/// ACGT-ordered, so `MAPPING` here is `[A=1, C=2, G=4, T=8]`.
pub fn s_codon_to_aa(codon: [u8; 3], genetic_code: &[u8; 64]) -> u8 {
    const MAPPING: [u8; 4] = [1, 2, 4, 8];

    let c0 = codon[0];
    let c1 = codon[1];
    let c2 = codon[2];
    if (c0 | c1 | c2) > 15 {
        if c0 == FENCE_SENTRY || c1 == FENCE_SENTRY || c2 == FENCE_SENTRY {
            return FENCE_SENTRY;
        }
    }

    // Fast path for unambiguous codons (each base is a single bit, i.e. one of
    // 1/2/4/8). This covers the vast majority of bytes in `blast_get_translation`
    // hot loops, eliminating the 4×4×4 ambiguity enumeration.
    let all_single = c0 <= 15
        && c1 <= 15
        && c2 <= 15
        && c0 != 0
        && (c0 & (c0 - 1)) == 0
        && c1 != 0
        && (c1 & (c1 - 1)) == 0
        && c2 != 0
        && (c2 & (c2 - 1)) == 0;
    if all_single {
        let i = c0.trailing_zeros() as usize;
        let j = c1.trailing_zeros() as usize;
        let k = c2.trailing_zeros() as usize;
        return genetic_code[i * 16 + j * 4 + k];
    }

    let mut aa: u8 = 0;
    for i in 0..4 {
        if c0 & MAPPING[i] != 0 {
            let index0 = i * 16;
            for j in 0..4 {
                if c1 & MAPPING[j] != 0 {
                    let index1 = index0 + j * 4;
                    for k in 0..4 {
                        if c2 & MAPPING[k] != 0 {
                            let taa = genetic_code[index1 + k];
                            if aa == 0 {
                                aa = taa;
                            } else if taa != aa {
                                return crate::encoding::NCBISTDAA_X;
                            }
                        }
                    }
                }
            }
        }
    }
    aa
}

pub const INIT_NUM_ELEMENTS: usize = 8;
pub const INIT_NUM_GEN_CODES: usize = 30;
pub const GENCODE_STRLEN: usize = 64;
pub const BLASTERR_INVALIDPARAM: i16 = 75;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EBlastStage {
    None = 0,
    PrelimSearch = 1,
    TracebackSearch = 2,
    Both = 3,
}

impl Default for EBlastStage {
    fn default() -> Self {
        Self::None
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SBlastProgress {
    pub stage: EBlastStage,
    pub user_data: Option<usize>,
}

/// Port of NCBI `SBlastProgressNew` (`blast_util.c:1387`).
pub fn s_blast_progress_new(user_data: Option<usize>) -> SBlastProgress {
    SBlastProgress {
        stage: EBlastStage::None,
        user_data,
    }
}

/// Rust equivalent of NCBI `SBlastProgressFree` (`blast_util.c:1397`).
pub fn s_blast_progress_free(_: Option<SBlastProgress>) -> Option<SBlastProgress> {
    None
}

/// Port of NCBI `SBlastProgressReset` (`blast_util.c:1406`).
pub fn s_blast_progress_reset(progress_info: Option<&mut SBlastProgress>) {
    if let Some(progress_info) = progress_info {
        progress_info.stage = EBlastStage::PrelimSearch;
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SDynamicUint4Array {
    pub num_used: u32,
    pub num_allocated: u32,
    pub data: Vec<u32>,
}

/// Port of NCBI `DynamicUint4ArrayNewEx` (`blast_dynarray.c:56`).
pub fn dynamic_uint4_array_new_ex(init_num_elements: u32) -> SDynamicUint4Array {
    SDynamicUint4Array {
        num_used: 0,
        num_allocated: init_num_elements,
        data: Vec::with_capacity(init_num_elements as usize),
    }
}

/// Port of NCBI `DynamicUint4ArrayNew` (`blast_dynarray.c:44`).
pub fn dynamic_uint4_array_new() -> SDynamicUint4Array {
    dynamic_uint4_array_new_ex(INIT_NUM_ELEMENTS as u32)
}

/// Rust equivalent of NCBI `DynamicUint4ArrayFree` (`blast_dynarray.c:76`).
pub fn dynamic_uint4_array_free(_: Option<SDynamicUint4Array>) -> Option<SDynamicUint4Array> {
    None
}

/// Port of NCBI `s_DynamicUint4Array_ReallocIfNecessary` (`blast_dynarray.c:94`).
pub fn s_dynamic_uint4_array_realloc_if_necessary(arr: &mut SDynamicUint4Array) -> i16 {
    if arr.num_used + 1 > arr.num_allocated {
        let new_allocated = if arr.num_allocated == 0 {
            1
        } else {
            arr.num_allocated * 2
        };
        arr.data
            .reserve((new_allocated as usize).saturating_sub(arr.data.capacity()));
        arr.num_allocated = new_allocated;
    }
    0
}

/// Port of NCBI `DynamicUint4Array_Append` (`blast_dynarray.c:113`).
pub fn dynamic_uint4_array_append(arr: &mut SDynamicUint4Array, element: u32) -> i16 {
    let retval = s_dynamic_uint4_array_realloc_if_necessary(arr);
    if retval != 0 {
        return retval;
    }
    arr.data.push(element);
    arr.num_used += 1;
    retval
}

/// Port of NCBI `DynamicUint4Array_Dup` (`blast_dynarray.c:129`).
pub fn dynamic_uint4_array_dup(src: Option<&SDynamicUint4Array>) -> Option<SDynamicUint4Array> {
    let src = src?;
    let mut retval = dynamic_uint4_array_new_ex(src.num_allocated);
    retval
        .data
        .extend_from_slice(&src.data[..src.num_used as usize]);
    retval.num_used = src.num_used;
    Some(retval)
}

/// Port of NCBI `DynamicUint4Array_Copy` (`blast_dynarray.c:146`).
pub fn dynamic_uint4_array_copy(dest: &mut SDynamicUint4Array, src: &SDynamicUint4Array) -> i32 {
    dest.data.clear();
    dest.data
        .extend_from_slice(&src.data[..src.num_used as usize]);
    dest.num_used = src.num_used;
    dest.num_allocated = src.num_allocated;
    0
}

/// Port of NCBI `DynamicUint4Array_AreEqual` (`blast_dynarray.c:182`).
pub fn dynamic_uint4_array_are_equal(
    a: Option<&SDynamicUint4Array>,
    b: Option<&SDynamicUint4Array>,
) -> bool {
    match (a, b) {
        (None, None) => true,
        (Some(_), None) | (None, Some(_)) => false,
        (Some(a), Some(b)) => {
            a.num_used == b.num_used
                && a.data[..a.num_used as usize] == b.data[..b.num_used as usize]
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SDynamicInt4Array {
    pub num_used: u32,
    pub num_allocated: u32,
    pub data: Vec<i32>,
}

/// Port of NCBI `DynamicInt4ArrayNew` (`blast_dynarray.c:213`).
pub fn dynamic_int4_array_new() -> SDynamicInt4Array {
    SDynamicInt4Array {
        num_used: 0,
        num_allocated: INIT_NUM_ELEMENTS as u32,
        data: Vec::with_capacity(INIT_NUM_ELEMENTS),
    }
}

/// Rust equivalent of NCBI `DynamicInt4ArrayFree` (`blast_dynarray.c:231`).
pub fn dynamic_int4_array_free(_: Option<SDynamicInt4Array>) -> Option<SDynamicInt4Array> {
    None
}

/// Port of NCBI `s_DynamicInt4Array_ReallocIfNecessary` (`blast_dynarray.c:249`).
pub fn s_dynamic_int4_array_realloc_if_necessary(arr: &mut SDynamicInt4Array) -> i16 {
    if arr.num_used + 1 > arr.num_allocated {
        let new_allocated = if arr.num_allocated == 0 {
            1
        } else {
            arr.num_allocated * 2
        };
        arr.data
            .reserve((new_allocated as usize).saturating_sub(arr.data.capacity()));
        arr.num_allocated = new_allocated;
    }
    0
}

/// Port of NCBI `DynamicInt4Array_Append` (`blast_dynarray.c:268`).
pub fn dynamic_int4_array_append(arr: &mut SDynamicInt4Array, element: i32) -> i16 {
    let retval = s_dynamic_int4_array_realloc_if_necessary(arr);
    if retval != 0 {
        return retval;
    }
    arr.data.push(element);
    arr.num_used += 1;
    retval
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SGenCodeNode {
    pub gc_id: u32,
    pub gc_str: Option<Vec<u8>>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SDynamicSGenCodeNodeArray {
    pub num_used: u32,
    pub num_allocated: u32,
    pub data: Vec<SGenCodeNode>,
}

/// Port of NCBI `DynamicSGenCodeNodeArrayNew` (`blast_dynarray.c:292`).
pub fn dynamic_s_gen_code_node_array_new() -> SDynamicSGenCodeNodeArray {
    SDynamicSGenCodeNodeArray {
        num_used: 0,
        num_allocated: INIT_NUM_GEN_CODES as u32,
        data: Vec::with_capacity(INIT_NUM_GEN_CODES),
    }
}

/// Rust equivalent of NCBI `DynamicSGenCodeNodeArrayFree` (`blast_dynarray.c:201`).
pub fn dynamic_s_gen_code_node_array_free(
    _: Option<SDynamicSGenCodeNodeArray>,
) -> Option<SDynamicSGenCodeNodeArray> {
    None
}

/// Port of NCBI `s_DynamicSGenCodeNodeArray_ReallocIfNecessary` (`blast_dynarray.c:313`).
pub fn s_dynamic_s_gen_code_node_array_realloc_if_necessary(
    arr: &mut SDynamicSGenCodeNodeArray,
) -> i16 {
    if arr.num_used + 1 > arr.num_allocated {
        let new_allocated = if arr.num_allocated == 0 {
            1
        } else {
            arr.num_allocated * 2
        };
        arr.data
            .reserve((new_allocated as usize).saturating_sub(arr.data.capacity()));
        arr.num_allocated = new_allocated;
    }
    0
}

/// Port of NCBI `s_DynamicSGenCodeNodeArray_IsSorted` (`blast_dynarray.c:335`).
pub fn s_dynamic_s_gen_code_node_array_is_sorted(arr: Option<&SDynamicSGenCodeNodeArray>) -> bool {
    let Some(arr) = arr else {
        return true;
    };
    if arr.num_used <= 1 {
        return true;
    }
    arr.data[..arr.num_used as usize]
        .windows(2)
        .all(|pair| pair[0].gc_id <= pair[1].gc_id)
}

/// Port of NCBI `s_SGenCodeNodeCompare` (`blast_dynarray.c:354`).
pub fn s_s_gen_code_node_compare(a: &SGenCodeNode, b: &SGenCodeNode) -> i32 {
    match a.gc_id.cmp(&b.gc_id) {
        std::cmp::Ordering::Less => -1,
        std::cmp::Ordering::Equal => 0,
        std::cmp::Ordering::Greater => 1,
    }
}

/// Port of NCBI `s_DynamicSGenCodeNodeArray_Sort` (`blast_dynarray.c:368`).
pub fn s_dynamic_s_gen_code_node_array_sort(arr: Option<&mut SDynamicSGenCodeNodeArray>) {
    let Some(arr) = arr else {
        return;
    };
    if arr.num_used <= 1 || s_dynamic_s_gen_code_node_array_is_sorted(Some(arr)) {
        return;
    }
    arr.data[..arr.num_used as usize].sort_by_key(|node| node.gc_id);
}

/// Port of NCBI `DynamicSGenCodeNodeArray_Append` (`blast_dynarray.c:384`).
pub fn dynamic_s_gen_code_node_array_append(
    arr: &mut SDynamicSGenCodeNodeArray,
    element: SGenCodeNode,
) -> i16 {
    let Some(gc_str) = element.gc_str.as_ref() else {
        return BLASTERR_INVALIDPARAM;
    };

    if dynamic_s_gen_code_node_array_find(arr, element.gc_id).is_some() {
        return 0;
    }

    let retval = s_dynamic_s_gen_code_node_array_realloc_if_necessary(arr);
    if retval != 0 {
        return retval;
    }

    let mut copied = vec![0; GENCODE_STRLEN];
    let len = gc_str.len().min(GENCODE_STRLEN);
    copied[..len].copy_from_slice(&gc_str[..len]);
    arr.data.push(SGenCodeNode {
        gc_id: element.gc_id,
        gc_str: Some(copied),
    });
    arr.num_used += 1;
    s_dynamic_s_gen_code_node_array_sort(Some(arr));
    retval
}

/// Port of NCBI `s_DynamicSGenCodeNodeArray_BinSearch` (`blast_dynarray.c:417`).
pub fn s_dynamic_s_gen_code_node_array_bin_search(
    arr: &SDynamicSGenCodeNodeArray,
    gen_code_id: u32,
) -> i32 {
    let mut b = 0i32;
    let mut e = arr.num_used as i32;
    while b < e - 1 {
        let m = (b + e) / 2;
        if arr.data[m as usize].gc_id > gen_code_id {
            e = m;
        } else {
            b = m;
        }
    }
    b
}

/// Port of NCBI `DynamicSGenCodeNodeArray_Find` (`blast_dynarray.c:438`).
pub fn dynamic_s_gen_code_node_array_find(
    arr: &SDynamicSGenCodeNodeArray,
    gen_code_id: u32,
) -> Option<&[u8]> {
    let index = s_dynamic_s_gen_code_node_array_bin_search(arr, gen_code_id) as usize;
    if index < arr.num_used as usize && arr.data[index].gc_id == gen_code_id {
        arr.data[index].gc_str.as_deref()
    } else {
        None
    }
}

static GEN_CODE_SINGLETON: OnceLock<Mutex<Option<SDynamicSGenCodeNodeArray>>> = OnceLock::new();

fn gen_code_singleton() -> &'static Mutex<Option<SDynamicSGenCodeNodeArray>> {
    GEN_CODE_SINGLETON.get_or_init(|| Mutex::new(None))
}

/// Port of NCBI `GenCodeSingletonInit`.
pub fn gen_code_singleton_init() {
    let mut singleton = gen_code_singleton()
        .lock()
        .expect("genetic-code singleton lock poisoned");
    if singleton.is_none() {
        *singleton = Some(dynamic_s_gen_code_node_array_new());
    }
}

/// Port of NCBI `GenCodeSingletonFini`.
pub fn gen_code_singleton_fini() {
    let mut singleton = gen_code_singleton()
        .lock()
        .expect("genetic-code singleton lock poisoned");
    let owned = singleton.take();
    *singleton = dynamic_s_gen_code_node_array_free(owned);
}

/// Port of NCBI `GenCodeSingletonAdd`.
pub fn gen_code_singleton_add(gen_code_id: u32, gen_code_str: &[u8]) -> i16 {
    let mut singleton = gen_code_singleton()
        .lock()
        .expect("genetic-code singleton lock poisoned");
    if singleton.is_none() {
        *singleton = Some(dynamic_s_gen_code_node_array_new());
    }
    let Some(arr) = singleton.as_mut() else {
        return BLASTERR_INVALIDPARAM;
    };
    dynamic_s_gen_code_node_array_append(
        arr,
        SGenCodeNode {
            gc_id: gen_code_id,
            gc_str: Some(gen_code_str.to_vec()),
        },
    )
}

/// Port-shaped Rust accessor for NCBI `GenCodeSingletonFind`.
pub fn gen_code_singleton_find(gen_code_id: u32) -> Option<Vec<u8>> {
    let singleton = gen_code_singleton()
        .lock()
        .expect("genetic-code singleton lock poisoned");
    let arr = singleton.as_ref()?;
    dynamic_s_gen_code_node_array_find(arr, gen_code_id).map(|code| code.to_vec())
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SSeqRange {
    pub left: i32,
    pub right: i32,
}

/// Port of NCBI `SSeqRangeNew` (`blast_util.c:47`).
pub fn s_seq_range_new(start: i32, stop: i32) -> SSeqRange {
    SSeqRange {
        left: start,
        right: stop,
    }
}

/// Port of NCBI `SSeqRangeArrayLessThanOrEqual` (`blast_util.c:56`).
pub fn s_seq_range_array_less_than_or_equal(ranges: Option<&[SSeqRange]>, target: i32) -> i32 {
    let Some(ranges) = ranges else {
        return -1;
    };
    if ranges.is_empty() {
        return -1;
    }

    let mut b = 0usize;
    let mut e = ranges.len();
    while b < e - 1 {
        let m = (b + e) / 2;
        if ranges[m].left > target {
            e = m;
        } else {
            b = m;
        }
    }

    if target > ranges[b].right && b < ranges.len() - 1 {
        (b + 1) as i32
    } else {
        b as i32
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum SubjectMaskingType {
    #[default]
    NoSubjMasking,
    SoftSubjMasking,
    HardSubjMasking,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct CompressedNucSeq {
    pub buffer: Vec<u8>,
    pub offset: usize,
}

impl CompressedNucSeq {
    /// Returns the logical `compressed_nuc_seq` view, equivalent to C's
    /// `compressed_nuc_seq_start + offset`.
    pub fn as_logical_slice(&self) -> &[u8] {
        self.buffer.get(self.offset..).unwrap_or(&[])
    }

    /// Returns the allocation start, equivalent to C's `compressed_nuc_seq_start`.
    pub fn as_start_slice(&self) -> &[u8] {
        &self.buffer
    }
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct BlastSequenceBlk {
    pub sequence: Option<Vec<u8>>,
    pub sequence_start: Option<Vec<u8>>,
    pub length: i32,
    pub frame: i16,
    pub subject_strand: i16,
    pub oid: i32,
    pub sequence_allocated: bool,
    pub sequence_start_allocated: bool,
    pub sequence_start_nomask: Option<Vec<u8>>,
    pub sequence_nomask: Option<Vec<u8>>,
    pub nomask_allocated: bool,
    pub oof_sequence: Option<Vec<u8>>,
    pub oof_sequence_allocated: bool,
    pub compressed_nuc_seq: Option<CompressedNucSeq>,
    pub compressed_nuc_seq_start: Option<Vec<u8>>,
    pub lcase_mask_allocated: bool,
    pub chunk: i32,
    pub gen_code_string: Option<[u8; 64]>,
    pub seq_ranges: Option<Vec<SSeqRange>>,
    pub num_seq_ranges: u32,
    pub seq_ranges_allocated: bool,
    pub mask_type: SubjectMaskingType,
    pub bases_offset: u8,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SBlastTargetTranslation {
    pub program_number: crate::program::ProgramType,
    pub gen_code_string: [u8; 64],
    pub translations: Vec<Option<Vec<u8>>>,
    pub partial: bool,
    pub num_frames: i32,
    pub range: Option<Vec<i32>>,
    pub subject_blk: Option<BlastSequenceBlk>,
}

/// Rust ownership equivalent of `BlastTargetTranslationFree` (`blast_util.c:1247`).
pub fn blast_target_translation_free(
    _: Option<SBlastTargetTranslation>,
) -> Option<SBlastTargetTranslation> {
    None
}

/// Port of NCBI `BlastTargetTranslationNew` (`blast_util.c:1268`).
pub fn blast_target_translation_new(
    subject_blk: &mut BlastSequenceBlk,
    gen_code_string: &[u8; 64],
    program_number: crate::program::ProgramType,
    is_ooframe: bool,
    target: &mut Option<SBlastTargetTranslation>,
) -> i16 {
    let partial = !is_ooframe;
    let mut retval = SBlastTargetTranslation {
        program_number,
        gen_code_string: *gen_code_string,
        translations: vec![None; NUM_FRAMES],
        partial,
        num_frames: NUM_FRAMES as i32,
        range: None,
        subject_blk: None,
    };

    if partial {
        retval.range = Some(vec![0; 2 * NUM_FRAMES]);
        retval.subject_blk = Some(subject_blk.clone());
    } else {
        // C passes subject_blk->sequence_start (the sentinel-prefixed buffer) to
        // both GetReverseNuclSequence and BLAST_GetTranslation (blast_util.c:1300).
        let seq = subject_blk.sequence_start.as_deref().unwrap_or(&[]);
        if is_ooframe {
            // C: BLAST_GetAllTranslations(sequence_start, eBlastEncodingNcbi4na,
            //    length, gen_code, NULL, NULL, &oof_sequence) — all 6 frames of
            //    both strands, interleaved into the mixed-frame sequence
            //    (blast_util.c:1289). The previous single-strand frame-1
            //    blast_get_partial_translation call produced a wrong (half-length,
            //    one-strand) sequence.
            let length = subject_blk.length.max(0) as usize;
            let (buffer, frame_offsets) = blast_get_all_translations(seq, length, gen_code_string);
            subject_blk.oof_sequence = Some(blast_all_translations_mixed_seq(
                &buffer,
                &frame_offsets,
                length,
            ));
            subject_blk.oof_sequence_allocated = true;
        } else {
            let rev = get_reverse_nucl_sequence(seq, subject_blk.length.max(0) as usize);
            for context in 0..NUM_FRAMES {
                let frame = blast_context_to_frame(context as u32);
                let mut translation = vec![0; subject_blk.length.max(0) as usize / 3 + 2];
                blast_get_translation(
                    seq,
                    &rev,
                    subject_blk.length.max(0) as usize,
                    frame,
                    &mut translation,
                    gen_code_string,
                );
                retval.translations[context] = Some(translation);
            }
        }
    }

    *target = Some(retval);
    0
}

/// Rust equivalent of `__sfree` (`blast_util.c:34`).
pub fn sfree<T>(slot: &mut Option<T>) {
    *slot = None;
}

/// Port of NCBI `iexp` (`lookup_util.c:46`).
pub fn iexp(x: i32, n: i32) -> i32 {
    let mut r = 1i32;
    let mut y = x;
    let mut n = n;

    if n == 0 {
        return 1;
    }
    if x == 0 {
        return 0;
    }

    while n > 1 {
        if n % 2 == 1 {
            r = r.wrapping_mul(y);
        }
        n >>= 1;
        y = y.wrapping_mul(y);
    }
    r.wrapping_mul(y)
}

/// Port of NCBI `ilog2` (`lookup_util.c:70`).
pub fn ilog2(mut x: i64) -> i32 {
    let mut lg = 0;
    if x == 0 {
        return 0;
    }

    while {
        x >>= 1;
        x != 0
    } {
        lg += 1;
    }
    lg
}

/// Port of NCBI internal `fkm_output` (`lookup_util.c:96`).
pub fn fkm_output(
    a: &[i32],
    n: i32,
    p: i32,
    output: &mut [u8],
    cursor: &mut i32,
    alphabet: Option<&[u8]>,
) {
    if n % p != 0 {
        return;
    }
    for i in 1..=p {
        let Some(slot) = output.get_mut(*cursor as usize) else {
            return;
        };
        let value = a.get(i as usize).copied().unwrap_or_default() as usize;
        *slot = alphabet
            .and_then(|abc| abc.get(value).copied())
            .unwrap_or(value as u8);
        *cursor += 1;
    }
}

/// Port of NCBI internal `fkm` (`lookup_util.c:140`).
pub fn fkm(
    a: &mut [i32],
    n: i32,
    k: i32,
    output: &mut [u8],
    cursor: &mut i32,
    alphabet: Option<&[u8]>,
) {
    if n <= 0 || k <= 0 || a.len() <= n as usize {
        return;
    }
    fkm_output(a, n, 1, output, cursor, alphabet);

    let mut i = n;
    loop {
        a[i as usize] += 1;

        for j in 1..=(n - i) {
            a[(j + i) as usize] = a[j as usize];
        }

        fkm_output(a, n, i, output, cursor, alphabet);
        i = n;

        while i > 0 && a[i as usize] == k - 1 {
            i -= 1;
        }
        if i <= 0 {
            break;
        }
    }
}

/// Port of NCBI internal `FNV_hash` (`blast_nalookup.c:1390`).
pub fn fnv_hash(seq: &[u8], mask: u32) -> u32 {
    const FNV_PRIME: u32 = 16_777_619;
    const FNV_OFFSET_BASIS: u32 = 2_166_136_261;
    let mut hash = FNV_OFFSET_BASIS;
    for i in 0..4 {
        hash = hash.wrapping_mul(FNV_PRIME);
        hash ^= seq.get(i).copied().unwrap_or(0) as u32;
    }
    hash & mask
}

/// Port of NCBI internal `s_Popcount` (`blast_nalookup.c:1476`).
pub fn s_popcount(mut v: u32) -> u32 {
    if v == 0 {
        return 0;
    }
    v = v - ((v >> 1) & 0x5555_5555);
    v = (v & 0x3333_3333) + ((v >> 2) & 0x3333_3333);
    v = (v + (v >> 4)) & 0x0f0f_0f0f;
    v = v.wrapping_mul(0x0101_0101);
    v >> 24
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BlastSparseUint1Array {
    pub bitfield: Vec<u32>,
    pub values: Vec<u8>,
    pub counts: Vec<i32>,
    pub num_elements: u32,
    pub length: u32,
}

/// Rust ownership equivalent of NCBI `BlastSparseUint1ArrayFree`.
pub fn blast_sparse_uint1_array_free(
    _: Option<BlastSparseUint1Array>,
) -> Option<BlastSparseUint1Array> {
    None
}

/// Port of NCBI `BlastSparseUint1ArrayNew` (`blast_nalookup.c:1520`).
pub fn blast_sparse_uint1_array_new(
    bitfield: Option<&[u32]>,
    len: i64,
) -> Option<BlastSparseUint1Array> {
    let bitfield = bitfield?;
    if len < 0 {
        return None;
    }
    let length = (len as u64 >> crate::stat::PV_ARRAY_BTS) as usize;
    if length == 0 || bitfield.len() < length {
        return None;
    }

    let mut counts = vec![0i32; length];
    counts[0] = s_popcount(bitfield[0]) as i32;
    for i in 1..length {
        counts[i] = counts[i - 1] + s_popcount(bitfield[i]) as i32;
    }

    let num_elements = counts[length - 1].max(0) as u32;
    Some(BlastSparseUint1Array {
        bitfield: bitfield[..length].to_vec(),
        values: vec![0; num_elements as usize],
        counts,
        num_elements,
        length: length as u32,
    })
}

/// Port of NCBI `BlastSparseUint1ArrayGetIndex` (`blast_nalookup.c:1556`).
pub fn blast_sparse_uint1_array_get_index(
    array: Option<&BlastSparseUint1Array>,
    index: i64,
) -> i32 {
    let Some(array) = array else {
        return -1;
    };
    if index < 0 {
        return -1;
    }

    let idx = (index as u64 >> crate::stat::PV_ARRAY_BTS) as usize;
    let bit_number = (index as u32) & crate::stat::PV_ARRAY_MASK;
    if idx >= array.length as usize {
        return -1;
    }

    let mut bit_count = if idx > 0 { array.counts[idx - 1] } else { 0 };
    bit_count += s_popcount(array.bitfield[idx] & ((1u32 << bit_number) - 1)) as i32;
    bit_count += 1;
    bit_count - 1
}

/// Port-shaped equivalent of NCBI `BlastSparseUint1ArrayGetElement`.
pub fn blast_sparse_uint1_array_get_element(
    array: Option<&mut BlastSparseUint1Array>,
    index: i64,
) -> Option<&mut u8> {
    let array = array?;
    let sparse_index = blast_sparse_uint1_array_get_index(Some(array), index);
    if sparse_index < 0 || sparse_index > array.num_elements as i32 {
        return None;
    }
    array.values.get_mut(sparse_index as usize)
}

/// Port of NCBI `s_BlastSequenceBlkFreeSeqRanges` (`blast_util.c:90`).
pub fn s_blast_sequence_blk_free_seq_ranges(seq_blk: &mut BlastSequenceBlk) {
    if seq_blk.seq_ranges_allocated {
        seq_blk.seq_ranges = None;
        seq_blk.num_seq_ranges = 0;
        seq_blk.seq_ranges_allocated = false;
    }
}

/// Port of NCBI `BlastSeqBlkNew` (`blast_util.c:133`).
pub fn blast_seq_blk_new(retval: Option<&mut Option<BlastSequenceBlk>>) -> i16 {
    let Some(retval) = retval else {
        return -1;
    };
    *retval = Some(BlastSequenceBlk::default());
    0
}

/// Port of NCBI `BlastSetUp_SeqBlkNew` (`blast_util.c:101`).
pub fn blast_set_up_seq_blk_new(
    buffer: &[u8],
    length: i32,
    seq_blk: &mut Option<BlastSequenceBlk>,
    buffer_allocated: bool,
) -> i16 {
    if seq_blk.is_none() && blast_seq_blk_new(Some(seq_blk)) != 0 {
        return -1;
    }
    let seq_blk = seq_blk.as_mut().expect("sequence block allocated above");

    if buffer_allocated {
        seq_blk.sequence_start_allocated = true;
        seq_blk.sequence_start = Some(buffer.to_vec());
        seq_blk.sequence = Some(buffer.get(1..).unwrap_or(&[]).to_vec());
    } else {
        seq_blk.sequence = Some(buffer.to_vec());
        seq_blk.sequence_start = None;
    }

    seq_blk.sequence_start_nomask = seq_blk.sequence_start.clone();
    seq_blk.sequence_nomask = seq_blk.sequence.clone();
    seq_blk.nomask_allocated = false;
    seq_blk.length = length;
    seq_blk.bases_offset = 0;
    0
}

/// Port of NCBI `BlastSeqBlkSetSequence` (`blast_util.c:147`).
pub fn blast_seq_blk_set_sequence(
    seq_blk: Option<&mut BlastSequenceBlk>,
    sequence: &[u8],
    seqlen: i32,
) -> i16 {
    let Some(seq_blk) = seq_blk else {
        return -1;
    };

    seq_blk.sequence_start_allocated = true;
    seq_blk.sequence_start = Some(sequence.to_vec());
    seq_blk.sequence = Some(sequence.get(1..).unwrap_or(&[]).to_vec());
    seq_blk.sequence_start_nomask = seq_blk.sequence_start.clone();
    seq_blk.sequence_nomask = seq_blk.sequence.clone();
    seq_blk.nomask_allocated = false;
    seq_blk.length = seqlen;
    seq_blk.oof_sequence = None;
    0
}

/// Port of NCBI `BlastSeqBlkSetCompressedSequence` (`blast_util.c:167`).
pub fn blast_seq_blk_set_compressed_sequence(
    seq_blk: Option<&mut BlastSequenceBlk>,
    sequence: &[u8],
) -> i16 {
    let Some(seq_blk) = seq_blk else {
        return -1;
    };

    seq_blk.sequence_allocated = true;
    seq_blk.sequence = Some(sequence.to_vec());
    seq_blk.oof_sequence = None;
    0
}

/// Port of NCBI `BlastSeqBlkSetSeqRanges` (`blast_util.c:182`).
pub fn blast_seq_blk_set_seq_ranges(
    seq_blk: Option<&mut BlastSequenceBlk>,
    seq_ranges: Option<&[SSeqRange]>,
    num_seq_ranges: u32,
    copy_seq_ranges: bool,
    mask_type: SubjectMaskingType,
) -> i16 {
    let Some(seq_blk) = seq_blk else {
        return -1;
    };
    let Some(seq_ranges) = seq_ranges else {
        return -1;
    };
    if num_seq_ranges == 0 || seq_ranges.len() < num_seq_ranges as usize {
        return -1;
    }

    s_blast_sequence_blk_free_seq_ranges(seq_blk);
    let mut tmp = seq_ranges[..num_seq_ranges as usize].to_vec();
    tmp[0].left = 0;
    tmp[num_seq_ranges as usize - 1].right = seq_blk.length;
    seq_blk.seq_ranges = Some(tmp);
    seq_blk.num_seq_ranges = num_seq_ranges;
    seq_blk.seq_ranges_allocated = copy_seq_ranges;
    seq_blk.mask_type = mask_type;
    0
}

/// Port of NCBI `BlastSequenceBlkClean` (`blast_util.c:220`).
pub fn blast_sequence_blk_clean(seq_blk: Option<&mut BlastSequenceBlk>) {
    let Some(seq_blk) = seq_blk else {
        return;
    };

    seq_blk.sequence = None;
    seq_blk.sequence_start = None;
    seq_blk.sequence_start_nomask = None;
    seq_blk.sequence_nomask = None;
    seq_blk.oof_sequence = None;
    seq_blk.compressed_nuc_seq = None;
    seq_blk.compressed_nuc_seq_start = None;
    seq_blk.sequence_allocated = false;
    seq_blk.sequence_start_allocated = false;
    seq_blk.nomask_allocated = false;
    seq_blk.oof_sequence_allocated = false;
    seq_blk.lcase_mask_allocated = false;
    s_blast_sequence_blk_free_seq_ranges(seq_blk);
}

/// Port of NCBI `BlastSequenceBlkFree` (`blast_util.c:245`).
pub fn blast_sequence_blk_free(seq_blk: Option<BlastSequenceBlk>) -> Option<BlastSequenceBlk> {
    drop(seq_blk);
    None
}

/// Port of NCBI `BlastSequenceBlkCopy` (`blast_util.c:259`).
pub fn blast_sequence_blk_copy(copy: &mut Option<BlastSequenceBlk>, src: &BlastSequenceBlk) {
    let mut copied = src.clone();
    copied.sequence_allocated = false;
    copied.sequence_start_allocated = false;
    copied.oof_sequence_allocated = false;
    copied.lcase_mask_allocated = false;
    copied.seq_ranges_allocated = false;
    *copy = Some(copied);
}

/// Port of NCBI `BlastCompressBlastnaSequence` (`blast_util.c:459`).
pub fn blast_compress_blastna_sequence(seq_blk: Option<&mut BlastSequenceBlk>) -> i16 {
    let Some(seq_blk) = seq_blk else {
        return -1;
    };
    let Some(old_seq) = seq_blk.sequence.as_ref() else {
        return -1;
    };
    let len = seq_blk.length.max(0) as usize;
    if old_seq.len() < len {
        return -1;
    }

    let mut new_seq = vec![0u8; len + 3];
    let max_start = 3usize.min(len);
    let mut curr_letter = 0u8;
    for i in 0..max_start {
        curr_letter = (curr_letter << 2) | (old_seq[i] & 3);
        new_seq[i + 3 - max_start] = curr_letter;
    }

    for i in max_start..len {
        curr_letter = (curr_letter << 2) | (old_seq[i] & 3);
        new_seq[i + 3 - max_start] = curr_letter;
    }

    for i in 0..max_start {
        curr_letter <<= 2;
        new_seq[3 + len - (max_start - i)] = curr_letter;
    }

    seq_blk.compressed_nuc_seq_start = Some(new_seq.clone());
    seq_blk.compressed_nuc_seq = Some(CompressedNucSeq {
        buffer: new_seq,
        offset: 3,
    });
    0
}

/// Port of NCBI `BSearchInt4` (`blast_util.c:1231`).
/// Returns the lower bracket index from a sorted ascending array. Like the C
/// helper, values smaller than the first element and empty arrays return `0`.
pub fn b_search_int4(n: i32, a: &[i32]) -> i32 {
    let mut b = 0usize;
    let mut e = a.len();

    while b < e.saturating_sub(1) {
        let m = (b + e) / 2;
        if a[m] > n {
            e = m;
        } else {
            b = m;
        }
    }

    b as i32
}

/// Port of NCBI `BLAST_StrToUpper` (`blast_util.c:1352`).
/// Returns `None` for C's NULL input; otherwise duplicates the byte string and
/// uppercases ASCII letters in the duplicate.
pub fn blast_str_to_upper(string: Option<&[u8]>) -> Option<Vec<u8>> {
    let string = string?;
    let mut retval = string.to_vec();

    for p in &mut retval {
        *p = p.to_ascii_uppercase();
    }

    Some(retval)
}

/// Rust equivalent of `BlastMemDup` (`ncbi_std.c:35`).
pub fn blast_mem_dup<T: Clone>(orig: Option<&T>) -> Option<T> {
    orig.cloned()
}

/// `eBlastEncodingNucleotide` (`blast_encoding.h:55`).
pub const E_BLAST_ENCODING_NUCLEOTIDE: i32 = 1;
/// `eBlastEncodingNcbi4na` (`blast_encoding.h:59`).
pub const E_BLAST_ENCODING_NCBI4NA: i32 = 2;
/// `eBlastEncodingNcbi2na` (`blast_encoding.h:60`).
pub const E_BLAST_ENCODING_NCBI2NA: i32 = 3;

/// Port of NCBI `BLAST_PackDNA` (`blast_util.c:870`).
/// Returns `(status, packed_seq)`, where the packed sequence includes the
/// final BLAST remainder byte.
pub fn blast_pack_dna(buffer: &[u8], length: usize, encoding: i32) -> (i16, Option<Vec<u8>>) {
    let new_length = length / 4 + 1;
    let mut new_buffer = vec![0u8; new_length];
    let mut index = 0usize;
    let mut new_index = 0usize;

    while new_index < new_length - 1 {
        new_buffer[new_index] = pack_dna_base(buffer[index], encoding) << 6
            | pack_dna_base(buffer[index + 1], encoding) << 4
            | pack_dna_base(buffer[index + 2], encoding) << 2
            | pack_dna_base(buffer[index + 3], encoding);
        new_index += 1;
        index += 4;
    }

    new_buffer[new_index] = (length % 4) as u8;
    while index < length {
        let shift = match index % 4 {
            0 => 6,
            1 => 4,
            2 => 2,
            _ => unreachable!(),
        };
        new_buffer[new_index] |= pack_dna_base(buffer[index], encoding) << shift;
        index += 1;
    }

    (0, Some(new_buffer))
}

fn pack_dna_base(base: u8, encoding: i32) -> u8 {
    if encoding == E_BLAST_ENCODING_NUCLEOTIDE {
        base & 0x03
    } else {
        crate::encoding::NCBI4NA_TO_BLASTNA[base as usize] & 0x03
    }
}

/// Port of NCBI `s_BlastGetTranslationTable` (`blast_util.c:990`).
/// The returned table is indexed by three packed NCBI2na bases.
pub fn s_blast_get_translation_table(
    genetic_code: Option<&[u8; 64]>,
    reverse_complement: bool,
) -> Option<[u8; 64]> {
    let genetic_code = genetic_code?;
    let mut translation = [0u8; 64];

    for index1 in 0..4 {
        for index2 in 0..4 {
            for index3 in 0..4 {
                if reverse_complement {
                    translation[(index3 << 4) + (index2 << 2) + index1] = translate_codon(
                        3 - index1 as u8,
                        3 - index2 as u8,
                        3 - index3 as u8,
                        genetic_code,
                    );
                } else {
                    translation[(index1 << 4) + (index2 << 2) + index3] =
                        translate_codon(index1 as u8, index2 as u8, index3 as u8, genetic_code);
                }
            }
        }
    }

    Some(translation)
}

/// Port of NCBI `BLAST_TranslateCompressedSequence` (`blast_util.c:508`).
/// Translates packed NCBI2na sequence bytes with NULLB sentinels around the
/// protein output. Negative frames expect the reverse-complement translation
/// table from [`s_blast_get_translation_table`], matching the C caller.
pub fn blast_translate_compressed_sequence(
    translation: &[u8; 64],
    length: usize,
    nt_seq: Option<&[u8]>,
    frame: i32,
    prot_seq: Option<&mut [u8]>,
) -> usize {
    let Some(nt_seq) = nt_seq else {
        return 0;
    };
    let Some(prot_seq) = prot_seq else {
        return 0;
    };
    if frame == 0 || length + 1 < frame.unsigned_abs() as usize + CODON_LENGTH {
        return 0;
    }
    if nt_seq.is_empty() {
        return 0;
    }

    prot_seq[0] = NULLB;
    let prot_seq_start = 1usize;
    let mut out = prot_seq_start;
    let mut codon = 0usize;
    let remainder = length % 4;

    if frame > 0 {
        let full_bytes = length / 4usize;
        let nt_seq_end = full_bytes.saturating_sub(1);
        let last_remainder = (4 * (length / 4) + 1 - frame as usize) % CODON_LENGTH;
        let total_remainder = last_remainder + remainder;

        let mut state = frame as usize - 1;
        let mut nt = 0usize;
        let mut byte_value = nt_seq[nt] as usize;

        while nt < nt_seq_end {
            match state {
                0 => {
                    codon = byte_value >> 2;
                    prot_seq[out] = translation[codon];
                    out += 1;

                    codon = (byte_value & 3) << 4;
                    nt += 1;
                    byte_value = nt_seq[nt] as usize;
                    codon += byte_value >> 4;
                    prot_seq[out] = translation[codon];
                    out += 1;
                    if nt >= nt_seq_end {
                        state = 2;
                    } else {
                        codon = (byte_value & 15) << 2;
                        nt += 1;
                        byte_value = nt_seq[nt] as usize;
                        codon += byte_value >> 6;
                        prot_seq[out] = translation[codon];
                        out += 1;
                        if nt >= nt_seq_end {
                            state = 1;
                        } else {
                            codon = byte_value & 63;
                            prot_seq[out] = translation[codon];
                            out += 1;
                            nt += 1;
                            byte_value = nt_seq[nt] as usize;
                            state = 0;
                        }
                    }
                }
                3 => {
                    codon = (byte_value & 3) << 4;
                    nt += 1;
                    byte_value = nt_seq[nt] as usize;
                    codon += byte_value >> 4;
                    prot_seq[out] = translation[codon];
                    out += 1;
                    if nt >= nt_seq_end {
                        state = 2;
                    } else {
                        codon = (byte_value & 15) << 2;
                        nt += 1;
                        byte_value = nt_seq[nt] as usize;
                        codon += byte_value >> 6;
                        prot_seq[out] = translation[codon];
                        out += 1;
                        if nt >= nt_seq_end {
                            state = 1;
                        } else {
                            codon = byte_value & 63;
                            prot_seq[out] = translation[codon];
                            out += 1;
                            nt += 1;
                            byte_value = nt_seq[nt] as usize;
                            state = 0;
                        }
                    }
                }
                2 => {
                    codon = (byte_value & 15) << 2;
                    nt += 1;
                    byte_value = nt_seq[nt] as usize;
                    codon += byte_value >> 6;
                    prot_seq[out] = translation[codon];
                    out += 1;
                    if nt >= nt_seq_end {
                        state = 1;
                    } else {
                        codon = byte_value & 63;
                        prot_seq[out] = translation[codon];
                        out += 1;
                        nt += 1;
                        byte_value = nt_seq[nt] as usize;
                        state = 0;
                    }
                }
                1 => {
                    codon = byte_value & 63;
                    prot_seq[out] = translation[codon];
                    out += 1;
                    nt += 1;
                    byte_value = nt_seq[nt] as usize;
                    state = 0;
                }
                _ => unreachable!(),
            }

            while nt + 10 < nt_seq_end {
                nt += 1;
                let byte_value1 = nt_seq[nt] as usize;
                nt += 1;
                let byte_value2 = nt_seq[nt] as usize;
                nt += 1;
                let byte_value3 = nt_seq[nt] as usize;

                codon = byte_value >> 2;
                prot_seq[out] = translation[codon];
                out += 1;

                codon = (byte_value & 3) << 4;
                codon += byte_value1 >> 4;
                prot_seq[out] = translation[codon];
                out += 1;

                nt += 1;
                let byte_value4 = nt_seq[nt] as usize;
                codon = (byte_value1 & 15) << 2;
                codon += byte_value2 >> 6;
                prot_seq[out] = translation[codon];
                out += 1;

                codon = byte_value2 & 63;
                nt += 1;
                let byte_value5 = nt_seq[nt] as usize;
                prot_seq[out] = translation[codon];
                out += 1;

                codon = byte_value3 >> 2;
                prot_seq[out] = translation[codon];
                out += 1;

                nt += 1;
                byte_value = nt_seq[nt] as usize;
                codon = (byte_value3 & 3) << 4;
                codon += byte_value4 >> 4;
                prot_seq[out] = translation[codon];
                out += 1;

                codon = (byte_value4 & 15) << 2;
                codon += byte_value5 >> 6;
                prot_seq[out] = translation[codon];
                out += 1;

                codon = byte_value5 & 63;
                prot_seq[out] = translation[codon];
                out += 1;
                state = 0;
            }
        }

        if state == 1 {
            byte_value = nt_seq[nt] as usize;
            codon = byte_value & 63;
            state = 0;
            prot_seq[out] = translation[codon];
            out += 1;
        } else if state == 0 {
            byte_value = nt_seq[nt] as usize;
            codon = byte_value >> 2;
            state = 3;
            prot_seq[out] = translation[codon];
            out += 1;
        }

        if full_bytes > 0 && total_remainder >= CODON_LENGTH {
            if nt_seq.len() <= nt_seq_end + 1 {
                prot_seq[out] = NULLB;
                return out - prot_seq_start;
            }
            byte_value = nt_seq[nt_seq_end] as usize;
            let last_byte = nt_seq[nt_seq_end + 1] as usize;
            if state == 0 {
                codon = last_byte >> 2;
            } else if state == 2 {
                codon = (byte_value & 15) << 2;
                codon += last_byte >> 6;
            } else if state == 3 {
                codon = (byte_value & 3) << 4;
                codon += last_byte >> 4;
            }
            prot_seq[out] = translation[codon];
            out += 1;
        }
    } else {
        let nt_seq_start = 0usize;
        let full_bytes = length / 4usize;
        let mut nt = full_bytes;
        let mut state = remainder as i32 + frame;
        let mut emitted_from_partial_tail = false;

        if state >= 0 {
            if nt >= nt_seq.len() {
                prot_seq[out] = NULLB;
                return out - prot_seq_start;
            }
            let last_byte = nt_seq[nt] as usize;
            nt = nt.saturating_sub(1);
            emitted_from_partial_tail = true;
            if state == 0 {
                codon = last_byte >> 6;
                let byte_value = nt_seq[nt] as usize;
                codon += (byte_value & 15) << 2;
                state = 1;
            } else if state == 1 {
                codon = last_byte >> 4;
                let byte_value = nt_seq[nt] as usize;
                codon += (byte_value & 3) << 4;
                state = 2;
            } else if state == 2 {
                codon = last_byte >> 2;
                state = 3;
            }
            prot_seq[out] = translation[codon];
            out += 1;
        } else {
            state = 3 + remainder as i32 + frame + 1;
            nt = nt.saturating_sub(1);
        }

        let mut byte_value = nt_seq[nt] as usize;
        while nt > nt_seq_start {
            match state {
                3 => {
                    codon = byte_value & 63;
                    prot_seq[out] = translation[codon];
                    out += 1;

                    codon = byte_value >> 6;
                    nt -= 1;
                    byte_value = nt_seq[nt] as usize;
                    codon += (byte_value & 15) << 2;
                    prot_seq[out] = translation[codon];
                    out += 1;
                    if nt <= nt_seq_start {
                        state = 1;
                    } else {
                        codon = byte_value >> 4;
                        nt -= 1;
                        byte_value = nt_seq[nt] as usize;
                        codon += (byte_value & 3) << 4;
                        prot_seq[out] = translation[codon];
                        out += 1;
                        if nt <= nt_seq_start {
                            state = 2;
                        } else {
                            codon = byte_value >> 2;
                            prot_seq[out] = translation[codon];
                            out += 1;
                            nt -= 1;
                            byte_value = nt_seq[nt] as usize;
                            state = 3;
                        }
                    }
                }
                0 => {
                    codon = byte_value >> 6;
                    nt -= 1;
                    byte_value = nt_seq[nt] as usize;
                    codon += (byte_value & 15) << 2;
                    prot_seq[out] = translation[codon];
                    out += 1;
                    if nt <= nt_seq_start {
                        state = 1;
                    } else {
                        codon = byte_value >> 4;
                        nt -= 1;
                        byte_value = nt_seq[nt] as usize;
                        codon += (byte_value & 3) << 4;
                        prot_seq[out] = translation[codon];
                        out += 1;
                        if nt <= nt_seq_start {
                            state = 2;
                        } else {
                            codon = byte_value >> 2;
                            prot_seq[out] = translation[codon];
                            out += 1;
                            nt -= 1;
                            byte_value = nt_seq[nt] as usize;
                            state = 3;
                        }
                    }
                }
                1 => {
                    codon = byte_value >> 4;
                    nt -= 1;
                    byte_value = nt_seq[nt] as usize;
                    codon += (byte_value & 3) << 4;
                    prot_seq[out] = translation[codon];
                    out += 1;
                    if nt <= nt_seq_start {
                        state = 2;
                    } else {
                        codon = byte_value >> 2;
                        prot_seq[out] = translation[codon];
                        out += 1;
                        nt -= 1;
                        byte_value = nt_seq[nt] as usize;
                        state = 3;
                    }
                }
                2 => {
                    codon = byte_value >> 2;
                    prot_seq[out] = translation[codon];
                    out += 1;
                    nt -= 1;
                    byte_value = nt_seq[nt] as usize;
                    state = 3;
                }
                _ => unreachable!(),
            }

            while nt > nt_seq_start + 10 {
                nt -= 1;
                let byte_value1 = nt_seq[nt] as usize;
                nt -= 1;
                let byte_value2 = nt_seq[nt] as usize;
                nt -= 1;
                let byte_value3 = nt_seq[nt] as usize;

                codon = byte_value & 63;
                prot_seq[out] = translation[codon];
                out += 1;

                codon = byte_value >> 6;
                codon += (byte_value1 & 15) << 2;
                prot_seq[out] = translation[codon];
                out += 1;

                nt -= 1;
                let byte_value4 = nt_seq[nt] as usize;
                codon = byte_value1 >> 4;
                codon += (byte_value2 & 3) << 4;
                prot_seq[out] = translation[codon];
                out += 1;

                codon = byte_value2 >> 2;
                prot_seq[out] = translation[codon];
                out += 1;

                nt -= 1;
                let byte_value5 = nt_seq[nt] as usize;

                codon = byte_value3 & 63;
                prot_seq[out] = translation[codon];
                out += 1;

                nt -= 1;
                byte_value = nt_seq[nt] as usize;
                codon = byte_value3 >> 6;
                codon += (byte_value4 & 15) << 2;
                prot_seq[out] = translation[codon];
                out += 1;

                codon = byte_value4 >> 4;
                codon += (byte_value5 & 3) << 4;
                prot_seq[out] = translation[codon];
                out += 1;

                codon = byte_value5 >> 2;
                prot_seq[out] = translation[codon];
                out += 1;
            }
        }

        byte_value = nt_seq[nt] as usize;
        if full_bytes == 0 && emitted_from_partial_tail {
        } else if state == 3 {
            codon = byte_value & 63;
            prot_seq[out] = translation[codon];
            out += 1;
        } else if state == 2 {
            codon = byte_value >> 2;
            prot_seq[out] = translation[codon];
            out += 1;
        }
    }

    prot_seq[out] = NULLB;
    out - prot_seq_start
}

/// Result bundle for [`blast_get_partial_translation`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PartialTranslation {
    pub translation_buffer: Option<Vec<u8>>,
    pub protein_length: Option<usize>,
    pub mixed_seq: Option<Vec<u8>>,
}

/// Port of NCBI `Blast_GetPartialTranslation` (`blast_util.c:1141`).
/// Returns translated data for one frame, or an out-of-frame mixed sequence
/// when `want_mixed_seq` is true.
pub fn blast_get_partial_translation(
    nucl_seq: &[u8],
    nucl_length: usize,
    frame: i32,
    genetic_code: &[u8; 64],
    want_translation_buffer: bool,
    want_protein_length: bool,
    want_mixed_seq: bool,
) -> Result<PartialTranslation, ()> {
    let nucl_seq_rev = if frame < 0 {
        get_reverse_nucl_sequence(nucl_seq, nucl_length)
    } else {
        Vec::new()
    };

    if !want_mixed_seq {
        let mut translation_buffer = vec![NULLB; nucl_length / CODON_LENGTH + 2];
        let length = blast_get_translation(
            nucl_seq,
            &nucl_seq_rev,
            nucl_length,
            frame,
            &mut translation_buffer,
            genetic_code,
        );
        return Ok(PartialTranslation {
            translation_buffer: want_translation_buffer.then_some(translation_buffer),
            protein_length: want_protein_length.then_some(length),
            mixed_seq: None,
        });
    }

    let frame_sign = if frame < 0 { -1 } else { 1 };
    let mut translation_buffer = vec![NULLB; nucl_length + 2];
    let mut frame_offsets = [0usize; CODON_LENGTH];
    let mut offset = 0usize;

    for index in 1..=CODON_LENGTH {
        let length = blast_get_translation(
            nucl_seq,
            &nucl_seq_rev,
            nucl_length,
            frame_sign * index as i32,
            &mut translation_buffer[offset..],
            genetic_code,
        );
        frame_offsets[index - 1] = offset;
        offset += length + 1;
    }

    let mut mixed_seq = vec![NULLB; nucl_length + 1];
    for (index, seq) in mixed_seq.iter_mut().enumerate() {
        *seq = translation_buffer[frame_offsets[index % CODON_LENGTH] + (index / CODON_LENGTH)];
    }

    Ok(PartialTranslation {
        translation_buffer: want_translation_buffer.then_some(translation_buffer),
        protein_length: want_protein_length.then_some(nucl_length),
        mixed_seq: Some(mixed_seq),
    })
}

/// Port of NCBI `BLAST_GetTranslatedProteinLength` (`blast_util.c:922`).
pub fn blast_get_translated_protein_length(nucleotide_length: usize, context: u32) -> usize {
    let offset = (context as usize) % CODON_LENGTH;
    if nucleotide_length == 0 || nucleotide_length <= offset {
        0
    } else {
        (nucleotide_length - offset) / CODON_LENGTH
    }
}

// Macros from `blast_util.h` and `blast_def.h` ported as `const`/`fn`.

/// `NULLB` (`ncbi_std.h:181`).
pub const NULLB: u8 = 0;
/// `FENCE_SENTRY` (`blast_util.h:364`).
pub const FENCE_SENTRY: u8 = 201;
/// `CODON_LENGTH` (`blast_def.h:63`).
pub const CODON_LENGTH: usize = 3;
/// `NUM_FRAMES` (`blast_def.h:88`).
pub const NUM_FRAMES: usize = 6;
/// `NUM_STRANDS` (`blast_def.h:93`).
pub const NUM_STRANDS: usize = 2;

/// Port of `IS_residue(x)` (`blast_util.h:48`).
#[inline]
pub fn is_residue(x: u8) -> bool {
    x <= 250
}

/// Port of `BLAST_ContextToFrame` (`blast_util.c:839`) for blastx-style
/// programs (eBlastTypeBlastx/Tblastx/RpsTblastn). Maps the 6 contexts
/// 0..5 to frames 1, 2, 3, -1, -2, -3.
pub fn blast_context_to_frame(context_number: u32) -> i32 {
    match (context_number as usize) % NUM_FRAMES {
        0 => 1,
        1 => 2,
        2 => 3,
        3 => -1,
        4 => -2,
        5 => -3,
        _ => unreachable!(),
    }
}

/// Port of `BLAST_FrameToContext` (`blast_util.c:1211`). Maps a reading frame
/// to a context index for any program class: translated programs (frame ±1..±3
/// → 0..5 via `frame-1` / `2-frame`), nucleotide programs (frame ±1 → 0/1), and
/// protein programs (frame 0 → 0). Unlike [`blast_context_to_frame`] this covers
/// all three branches, matching the C dispatch.
pub fn blast_frame_to_context(frame: i32, program: crate::program::ProgramType) -> i32 {
    use crate::program::{
        blast_query_is_nucleotide, blast_query_is_translated, blast_subject_is_nucleotide,
        blast_subject_is_translated,
    };
    if blast_query_is_translated(program) || blast_subject_is_translated(program) {
        debug_assert!((-3..=3).contains(&frame));
        if frame > 0 {
            frame - 1
        } else {
            2 - frame
        }
    } else if blast_query_is_nucleotide(program) || blast_subject_is_nucleotide(program) {
        debug_assert!(frame == 1 || frame == -1);
        if frame == 1 {
            0
        } else {
            1
        }
    } else {
        debug_assert!(frame == 0);
        0
    }
}

/// Port of `GetReverseNuclSequence` (`blast_util.c:807`). Builds the
/// reverse-complement of an NCBI4na buffer with NULLB sentinels at
/// positions 0 and length+1, matching the C layout used by
/// `BLAST_GetTranslation` (`query_seq_rev + 1` skips the leading NULLB).
pub fn get_reverse_nucl_sequence(sequence: &[u8], length: usize) -> Vec<u8> {
    // C: rev_sequence = malloc(length + 2);
    //    rev_sequence[0] = rev_sequence[length+1] = NULLB;
    let mut rev = vec![NULLB; length + 2];
    for index in 0..length {
        let b = sequence[index];
        rev[length - index] = if b == FENCE_SENTRY {
            FENCE_SENTRY
        } else {
            crate::encoding::complement_ncbi4na_base(b)
        };
    }
    rev
}

/// Port of `BLAST_GetTranslation` (`blast_util.c:427`). Translates one
/// frame of an NCBI4na nucleotide sequence into NCBIstdaa. Writes
/// `prot_seq[0] = NULLB` then residues at indices 1..=k then a trailing
/// NULLB; returns `k` (the number of residues written).
///
/// `query_seq_rev` must be the buffer returned by [`get_reverse_nucl_sequence`]
/// (length+2 bytes, leading NULLB at [0]). It is dereferenced only when
/// `frame < 0`. Frames in {1,2,3} read from `query_seq` directly; frames
/// in {-1,-2,-3} read from `&query_seq_rev[1..]`.
pub fn blast_get_translation(
    query_seq: &[u8],
    query_seq_rev: &[u8],
    nt_length: usize,
    frame: i32,
    prot_seq: &mut [u8],
    genetic_code: &[u8; 64],
) -> usize {
    // C: nucl_seq = (frame >= 0 ? query_seq : query_seq_rev + 1);
    let nucl_seq: &[u8] = if frame >= 0 {
        query_seq
    } else {
        &query_seq_rev[1..]
    };

    // C: prot_seq[0] = NULLB; index_prot = 1;
    prot_seq[0] = NULLB;
    let mut index_prot: usize = 1;

    // C: for (index = ABS(frame) - 1; index < nt_length - 2; index += CODON_LENGTH)
    let start = (frame.unsigned_abs() as usize) - 1;
    if nt_length >= 2 && nucl_seq.len() >= nt_length {
        let nucl_ptr = nucl_seq.as_ptr();
        let prot_ptr = prot_seq.as_mut_ptr();
        let end = nt_length - 2;
        let mut index = start;
        while index < end {
            // SAFETY: `index < end = nt_length - 2`, plus nt_length <= nucl_seq.len()
            // and prot_seq must have at least `(end - start) / 3 + 2` slots which
            // callers guarantee. Avoids bounds checks in the hot translation loop.
            let codon = unsafe {
                [
                    *nucl_ptr.add(index),
                    *nucl_ptr.add(index + 1),
                    *nucl_ptr.add(index + 2),
                ]
            };
            let residue = s_codon_to_aa(codon, genetic_code);
            // C: if (IS_residue(residue) || residue == FENCE_SENTRY)
            if is_residue(residue) || residue == FENCE_SENTRY {
                unsafe {
                    *prot_ptr.add(index_prot) = residue;
                }
                index_prot += 1;
            }
            index += CODON_LENGTH;
        }
    }
    prot_seq[index_prot] = NULLB;
    index_prot - 1
}

/// Port of `BLAST_GetAllTranslations` (`blast_util.c:1045`) for the
/// `eBlastEncodingNcbi4na` path only. Translates a nucleotide sequence
/// in all 6 frames into a single concatenated buffer separated by NULLB
/// sentinels, returning `(translation_buffer, frame_offsets)`.
///
/// Frame `ctx`'s residues are at
/// `&translation_buffer[frame_offsets[ctx] + 1 .. frame_offsets[ctx+1]]`
/// (matching C's `translation_buffer + frame_offsets[ctx] + 1`).
pub fn blast_get_all_translations(
    nucl_seq: &[u8],
    nucl_length: usize,
    genetic_code: &[u8; 64],
) -> (Vec<u8>, [u32; NUM_FRAMES + 1]) {
    // C: buffer_length = 2*(nucl_length+1)+2;
    let buffer_length = 2 * (nucl_length + 1) + 2;
    let mut translation_buffer = vec![0u8; buffer_length];
    let nucl_seq_rev = get_reverse_nucl_sequence(nucl_seq, nucl_length);

    let mut offset: u32 = 0;
    let mut frame_offsets = [0u32; NUM_FRAMES + 1];
    frame_offsets[0] = 0;

    for context in 0..NUM_FRAMES {
        let frame = blast_context_to_frame(context as u32);
        let length = blast_get_translation(
            nucl_seq,
            &nucl_seq_rev,
            nucl_length,
            frame,
            &mut translation_buffer[offset as usize..],
            genetic_code,
        );
        // C: offset += length + 1;  (1 extra byte for the inter-frame NULLB)
        offset += (length as u32) + 1;
        frame_offsets[context + 1] = offset;
    }

    (translation_buffer, frame_offsets)
}

/// Port of `BLAST_GetAllTranslations` for the `eBlastEncodingNcbi2na` path:
/// translate all six frames DIRECTLY from the packed NCBI2na subject, exactly
/// as NCBI's preliminary engine does (`blast_util.c:1067-1091`) — a forward and
/// a reverse-complement translation table feeding `BLAST_TranslateCompressedSequence`
/// per frame. This avoids the per-OID decode→BLASTNA→NCBI4na buffers the
/// `blast_get_all_translations` (ncbi4na) path requires.
///
/// Produces a byte-identical `translation_buffer` / `frame_offsets` to the
/// ncbi4na path for ambiguity-free subjects. Subjects carrying ambiguity codes
/// MUST stay on the ncbi4na path (`BLAST_TranslateCompressedSequence` cannot
/// represent ambiguities; NCBI's packed path likewise ignores them).
pub fn blast_get_all_translations_packed(
    packed_nucl: &[u8],
    nucl_length: usize,
    genetic_code: &[u8; 64],
) -> (Vec<u8>, [u32; NUM_FRAMES + 1]) {
    // C: buffer_length = 2*(nucl_length+1)+2;
    let buffer_length = 2 * (nucl_length + 1) + 2;
    let mut translation_buffer = vec![0u8; buffer_length];
    let table_fwd = s_blast_get_translation_table(Some(genetic_code), false)
        .expect("forward translation table");
    let table_rc = s_blast_get_translation_table(Some(genetic_code), true)
        .expect("reverse-complement translation table");

    let mut offset: u32 = 0;
    let mut frame_offsets = [0u32; NUM_FRAMES + 1];
    frame_offsets[0] = 0;

    for context in 0..NUM_FRAMES {
        let frame = blast_context_to_frame(context as u32);
        let table = if frame > 0 { &table_fwd } else { &table_rc };
        let length = blast_translate_compressed_sequence(
            table,
            nucl_length,
            Some(packed_nucl),
            frame,
            Some(&mut translation_buffer[offset as usize..]),
        );
        // C: offset += length + 1;  (1 extra byte for the inter-frame NULLB)
        offset += (length as u32) + 1;
        frame_offsets[context + 1] = offset;
    }

    (translation_buffer, frame_offsets)
}

/// Port of the out-of-frame mixed-sequence construction block of
/// `BLAST_GetAllTranslations` (`blast_util.c:1112-1126`). Given the 6-frame
/// concatenated `translation_buffer` and its `frame_offsets` (as produced by
/// [`blast_get_all_translations`]), builds the OOF mixed-frame sequence used by
/// out-of-frame gapping: length `2*nucl_length+3`, interleaving the three
/// frames of each strand base-by-base (`+1,+2,+3` then `-1,-2,-3`), with a
/// trailing NULLB sentinel. Split out as a separate helper so the common
/// `blast_get_all_translations` callers (which pass C's `mixed_seq_ptr == NULL`)
/// keep their `(buffer, frame_offsets)` signature.
pub fn blast_all_translations_mixed_seq(
    translation_buffer: &[u8],
    frame_offsets: &[u32; NUM_FRAMES + 1],
    nucl_length: usize,
) -> Vec<u8> {
    // C: malloc(2*nucl_length+3); the trailing byte is the final NULLB sentinel.
    let mut mixed_seq = vec![NULLB; 2 * nucl_length + 3];
    let mut out = 0usize;
    let mut index = 0usize;
    while index < NUM_FRAMES {
        for i in 0..=nucl_length {
            let context = i % CODON_LENGTH;
            let offset = i / CODON_LENGTH;
            mixed_seq[out] = translation_buffer[frame_offsets[index + context] as usize + offset];
            out += 1;
        }
        index += CODON_LENGTH;
    }
    // `out` is now 2*(nucl_length+1); mixed_seq[out] stays NULLB (C's `*seq = NULLB`).
    mixed_seq
}

/// Port of `BLAST_CreateMixedFrameDNATranslation` (`blast_util.c:931`). Builds
/// the query-side out-of-frame mixed-frame DNA sequence: for each query, the
/// three forward (then three reverse) reading frames are interleaved
/// base-by-base, prefixed by three NULLB sentinels, and written at the
/// first-frame context's `query_offset`. The result is stored in
/// `query_blk.oof_sequence`. Used only by out-of-frame (frameshift) gapping in
/// blastx; returns -1 if the query block has no `sequence`.
pub fn blast_create_mixed_frame_dna_translation(
    query_blk: &mut BlastSequenceBlk,
    query_info: &crate::queryinfo::QueryInfo,
) -> i16 {
    let total_length = query_info.seq_buf_len();
    // C: malloc(total_length+1) — 1 extra byte for the final sentinel.
    let mut buffer = vec![NULLB; total_length + 1];
    let contexts = &query_info.contexts;
    if contexts.is_empty() {
        return -1;
    }
    let last_context = contexts.len() - 1;
    let sequence = match query_blk.sequence.as_deref() {
        Some(s) => s,
        None => return -1,
    };

    let mut seq = 0usize; // running write cursor; mirrors C's post-incremented `seq`
    let mut wrote_any = false;
    let mut index = 0usize;
    while index <= last_context {
        // query_length == 0 indicates this context is not searched.
        if contexts[index].query_length == 0 {
            index += CODON_LENGTH;
            continue;
        }
        seq = contexts[index].query_offset.max(0) as usize;
        let mut length = [0i32; CODON_LENGTH];
        for (i, len) in length.iter_mut().enumerate() {
            buffer[seq] = NULLB;
            seq += 1;
            *len = contexts[index + i].query_length;
        }
        let mut i = 0usize;
        loop {
            let context = i % CODON_LENGTH;
            let offset = i / CODON_LENGTH;
            // Once one frame is past its end, we are done.
            if offset >= length[context].max(0) as usize {
                break;
            }
            let src = contexts[index + context].query_offset.max(0) as usize + offset;
            buffer[seq] = sequence[src];
            seq += 1;
            i += 1;
        }
        wrote_any = true;
        index += CODON_LENGTH;
    }
    // C: if (seq) *seq = NULLB; — trailing sentinel after the last written base.
    if wrote_any {
        buffer[seq] = NULLB;
    }

    query_blk.oof_sequence = Some(buffer);
    query_blk.oof_sequence_allocated = true;
    0
}

/// Standard genetic code (NCBI translation table 1).
/// Maps 64 codons (in NCBI2na order: A=0,C=1,G=2,T=3) to NCBIstdaa amino acids.
pub static STANDARD_GENETIC_CODE: [u8; 64] = [
    // AAA AAC AAG AAT ACA ACC ACG ACT AGA AGC AGG AGT ATA ATC ATG ATT
    10, 13, 10, 13, 18, 18, 18, 18, 16, 17, 16, 17, 9, 9, 12, 9,
    // CAA CAC CAG CAT CCA CCC CCG CCT CGA CGC CGG CGT CTA CTC CTG CTT
    15, 8, 15, 8, 14, 14, 14, 14, 16, 16, 16, 16, 11, 11, 11, 11,
    // GAA GAC GAG GAT GCA GCC GCG GCT GGA GGC GGG GGT GTA GTC GTG GTT
    5, 4, 5, 4, 1, 1, 1, 1, 7, 7, 7, 7, 19, 19, 19, 19,
    // TAA TAC TAG TAT TCA TCC TCG TCT TGA TGC TGG TGT TTA TTC TTG TTT
    25, 22, 25, 22, 17, 17, 17, 17, 25, 3, 20, 3, 11, 6, 11, 6,
];

/// Build a genetic code table by modifying the standard code.
const fn make_code(changes: &[(usize, u8)]) -> [u8; 64] {
    let mut t = STANDARD_GENETIC_CODE;
    let mut i = 0;
    while i < changes.len() {
        t[changes[i].0] = changes[i].1;
        i += 1;
    }
    t
}

/// Code 2: Vertebrate Mitochondrial (AGA=*, AGG=*, ATA=M, TGA=W)
pub static GENETIC_CODE_2: [u8; 64] = make_code(&[(8, 25), (10, 25), (12, 12), (56, 20)]);
/// Code 3: Yeast Mitochondrial (ATA=M, CTN=T, TGA=W)
pub static GENETIC_CODE_3: [u8; 64] =
    make_code(&[(12, 12), (28, 18), (29, 18), (30, 18), (31, 18), (56, 20)]);
/// Code 4: Mold/Protozoan Mitochondrial (TGA=W)
pub static GENETIC_CODE_4: [u8; 64] = make_code(&[(56, 20)]);
/// Code 5: Invertebrate Mitochondrial (AGA=S, AGG=S, ATA=M, TGA=W)
pub static GENETIC_CODE_5: [u8; 64] = make_code(&[(8, 17), (10, 17), (12, 12), (56, 20)]);
/// Code 6: Ciliate Nuclear (TAA=Q, TAG=Q)
pub static GENETIC_CODE_6: [u8; 64] = make_code(&[(48, 15), (50, 15)]);
/// Code 9: Echinoderm Mitochondrial (AAA=N, AGA=S, AGG=S, TGA=W)
pub static GENETIC_CODE_9: [u8; 64] = make_code(&[(0, 13), (8, 17), (10, 17), (56, 20)]);
/// Code 10: Euplotid Nuclear (TGA=C)
pub static GENETIC_CODE_10: [u8; 64] = make_code(&[(56, 3)]);
/// Code 12: Alternative Yeast Nuclear (CTG=S)
pub static GENETIC_CODE_12: [u8; 64] = make_code(&[(30, 17)]);
/// Code 13: Ascidian Mitochondrial (AGA=G, AGG=G, ATA=M, TGA=W)
pub static GENETIC_CODE_13: [u8; 64] = make_code(&[(8, 7), (10, 7), (12, 12), (56, 20)]);
/// Code 14: Alternative Flatworm Mitochondrial (AAA=N, AGA=S, AGG=S, TAA=Y, TGA=W)
pub static GENETIC_CODE_14: [u8; 64] = make_code(&[(0, 13), (8, 17), (10, 17), (48, 22), (56, 20)]);
/// Code 15: Blepharisma Nuclear (TAG=Q)
pub static GENETIC_CODE_15: [u8; 64] = make_code(&[(50, 15)]);
/// Code 16: Chlorophycean Mitochondrial (TAG=L)
pub static GENETIC_CODE_16: [u8; 64] = make_code(&[(50, 11)]);
/// Code 21: Trematode Mitochondrial (AAA=N, AGA=S, AGG=S, ATA=M, TGA=W)
pub static GENETIC_CODE_21: [u8; 64] = make_code(&[(0, 13), (8, 17), (10, 17), (12, 12), (56, 20)]);
/// Code 22: Scenedesmus obliquus Mitochondrial (TCA=*, TAG=L)
pub static GENETIC_CODE_22: [u8; 64] = make_code(&[(52, 25), (50, 11)]);
/// Code 23: Thraustochytrium Mitochondrial (TTA=*)
pub static GENETIC_CODE_23: [u8; 64] = make_code(&[(60, 25)]);
/// Code 24: Rhabdopleuridae Mitochondrial (AGA=S, AGG=K, TGA=W)
pub static GENETIC_CODE_24: [u8; 64] = make_code(&[(8, 17), (10, 10), (56, 20)]);
/// Code 25: Candidate Division SR1 (TGA=G)
pub static GENETIC_CODE_25: [u8; 64] = make_code(&[(56, 7)]);
/// Code 26: Pachysolen tannophilus Nuclear (CTG=A)
pub static GENETIC_CODE_26: [u8; 64] = make_code(&[(30, 1)]);
/// Code 27: Karyorelictea Nuclear (TAA=Q, TAG=Q, TGA=W)
pub static GENETIC_CODE_27: [u8; 64] = make_code(&[(48, 15), (50, 15), (56, 20)]);
/// Code 28: Condylostoma Nuclear (TAA=Q, TAG=Q, TGA=W)
pub static GENETIC_CODE_28: [u8; 64] = make_code(&[(48, 15), (50, 15), (56, 20)]);
/// Code 29: Mesodinium Nuclear (TAA=Y, TAG=Y)
pub static GENETIC_CODE_29: [u8; 64] = make_code(&[(48, 22), (50, 22)]);
/// Code 30: Peritrich Nuclear (TAA=E, TAG=E)
pub static GENETIC_CODE_30: [u8; 64] = make_code(&[(48, 5), (50, 5)]);
/// Code 31: Blastocrithidia Nuclear (TAA=E, TAG=E, TGA=W)
pub static GENETIC_CODE_31: [u8; 64] = make_code(&[(48, 5), (50, 5), (56, 20)]);
/// Code 33: Cephalodiscidae Mitochondrial (AGA=S, AGG=K, TAA=Y, TGA=W)
pub static GENETIC_CODE_33: [u8; 64] = make_code(&[(8, 17), (10, 10), (48, 22), (56, 20)]);

/// Look up a genetic code table by NCBI code number.
pub fn lookup_genetic_code(code: u8) -> &'static [u8; 64] {
    match code {
        1 | 11 => &STANDARD_GENETIC_CODE,
        2 => &GENETIC_CODE_2,
        3 => &GENETIC_CODE_3,
        4 => &GENETIC_CODE_4,
        5 => &GENETIC_CODE_5,
        6 => &GENETIC_CODE_6,
        9 => &GENETIC_CODE_9,
        10 => &GENETIC_CODE_10,
        12 => &GENETIC_CODE_12,
        13 => &GENETIC_CODE_13,
        14 => &GENETIC_CODE_14,
        15 => &GENETIC_CODE_15,
        16 => &GENETIC_CODE_16,
        21 => &GENETIC_CODE_21,
        22 => &GENETIC_CODE_22,
        23 => &GENETIC_CODE_23,
        24 => &GENETIC_CODE_24,
        25 => &GENETIC_CODE_25,
        26 => &GENETIC_CODE_26,
        27 => &GENETIC_CODE_27,
        28 => &GENETIC_CODE_28,
        29 => &GENETIC_CODE_29,
        30 => &GENETIC_CODE_30,
        31 => &GENETIC_CODE_31,
        33 => &GENETIC_CODE_33,
        _ => &STANDARD_GENETIC_CODE,
    }
}

/// Map protein coordinates back to nucleotide coordinates given a reading frame.
///
/// For blastx: maps query protein coords → query nucleotide coords.
/// For tblastn: maps subject protein coords → subject nucleotide coords.
///
/// - `prot_start`, `prot_end`: 0-based protein coordinates from ungapped extension
/// - `frame`: reading frame (1,2,3 for forward, -1,-2,-3 for reverse)
/// - `nuc_len`: total nucleotide sequence length
///
/// Returns 1-based nucleotide coordinates (start, end) as BLAST reports them.
pub fn protein_to_nuc_coords(
    prot_start: i32,
    prot_end: i32,
    frame: i32,
    nuc_len: i32,
) -> (i32, i32) {
    if frame > 0 {
        let offset = frame - 1;
        (prot_start * 3 + offset + 1, prot_end * 3 + offset)
    } else {
        // Reverse strand: BLAST reports minus-frame coordinates in strand
        // orientation, so the displayed interval runs high→low.
        let offset = (-frame) - 1;
        let rc_start = prot_start * 3 + offset;
        let rc_end = prot_end * 3 + offset - 1;
        (nuc_len - rc_start, nuc_len - rc_end)
    }
}

/// Map protein coordinates onto nucleotide coordinates within the oriented
/// frame-specific sequence used during translated search.
///
/// The returned interval is 0-based with an exclusive end. For negative frames
/// the coordinates are relative to the reverse-complement orientation that was
/// translated, matching the crate's internal convention for minus-strand
/// blastn HSP coordinates before output formatting maps them back to BLAST's
/// displayed coordinates.
pub fn protein_to_oriented_nuc_coords(
    prot_start: usize,
    prot_end: usize,
    frame: i32,
) -> (usize, usize) {
    let offset = frame.unsigned_abs() as usize - 1;
    (prot_start * 3 + offset, prot_end * 3 + offset)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_blast_progress_new_free_reset() {
        let mut progress = s_blast_progress_new(Some(0xfeed));
        assert_eq!(progress.stage, EBlastStage::None);
        assert_eq!(progress.user_data, Some(0xfeed));

        progress.stage = EBlastStage::TracebackSearch;
        s_blast_progress_reset(Some(&mut progress));
        assert_eq!(progress.stage, EBlastStage::PrelimSearch);
        assert_eq!(progress.user_data, Some(0xfeed));
        s_blast_progress_reset(None);
        assert_eq!(s_blast_progress_free(Some(progress)), None);
    }

    #[test]
    fn test_dynamic_uint4_array_append_dup_copy_equal() {
        let mut arr = dynamic_uint4_array_new();
        assert_eq!(arr.num_allocated, INIT_NUM_ELEMENTS as u32);
        for value in 0..9 {
            assert_eq!(dynamic_uint4_array_append(&mut arr, value), 0);
        }
        assert_eq!(arr.num_used, 9);
        assert_eq!(arr.num_allocated, 16);
        assert_eq!(&arr.data[..9], &[0, 1, 2, 3, 4, 5, 6, 7, 8]);

        let dup = dynamic_uint4_array_dup(Some(&arr)).unwrap();
        assert!(dynamic_uint4_array_are_equal(Some(&arr), Some(&dup)));
        assert!(!dynamic_uint4_array_are_equal(Some(&arr), None));
        assert!(dynamic_uint4_array_are_equal(None, None));

        let mut copied = dynamic_uint4_array_new_ex(1);
        assert_eq!(dynamic_uint4_array_copy(&mut copied, &arr), 0);
        assert!(dynamic_uint4_array_are_equal(Some(&arr), Some(&copied)));
        assert_eq!(dynamic_uint4_array_free(Some(copied)), None);
    }

    #[test]
    fn test_dynamic_int4_array_append_grows() {
        let mut arr = dynamic_int4_array_new();
        for value in -4..5 {
            assert_eq!(dynamic_int4_array_append(&mut arr, value), 0);
        }
        assert_eq!(arr.num_used, 9);
        assert_eq!(arr.num_allocated, 16);
        assert_eq!(&arr.data[..9], &[-4, -3, -2, -1, 0, 1, 2, 3, 4]);
        assert_eq!(dynamic_int4_array_free(Some(arr)), None);
    }

    #[test]
    fn test_dynamic_s_gen_code_node_array_append_sort_find() {
        let mut arr = dynamic_s_gen_code_node_array_new();
        let code_a = vec![b'A'; GENCODE_STRLEN];
        let code_b = vec![b'B'; GENCODE_STRLEN];

        assert_eq!(
            dynamic_s_gen_code_node_array_append(
                &mut arr,
                SGenCodeNode {
                    gc_id: 11,
                    gc_str: Some(code_b.clone()),
                },
            ),
            0
        );
        assert_eq!(
            dynamic_s_gen_code_node_array_append(
                &mut arr,
                SGenCodeNode {
                    gc_id: 1,
                    gc_str: Some(code_a.clone()),
                },
            ),
            0
        );
        assert!(s_dynamic_s_gen_code_node_array_is_sorted(Some(&arr)));
        assert_eq!(arr.data[0].gc_id, 1);
        assert_eq!(
            dynamic_s_gen_code_node_array_find(&arr, 1),
            Some(code_a.as_slice())
        );
        assert_eq!(
            dynamic_s_gen_code_node_array_find(&arr, 11),
            Some(code_b.as_slice())
        );
        assert_eq!(dynamic_s_gen_code_node_array_find(&arr, 99), None);

        assert_eq!(
            dynamic_s_gen_code_node_array_append(
                &mut arr,
                SGenCodeNode {
                    gc_id: 1,
                    gc_str: Some(vec![b'Z'; GENCODE_STRLEN]),
                },
            ),
            0
        );
        assert_eq!(arr.num_used, 2);
        assert_eq!(
            dynamic_s_gen_code_node_array_append(
                &mut arr,
                SGenCodeNode {
                    gc_id: 2,
                    gc_str: None,
                },
            ),
            BLASTERR_INVALIDPARAM
        );
        assert_eq!(s_s_gen_code_node_compare(&arr.data[0], &arr.data[1]), -1);
        assert_eq!(dynamic_s_gen_code_node_array_free(Some(arr)), None);
    }

    #[test]
    fn test_gen_code_singleton_init_add_find_fini() {
        gen_code_singleton_fini();
        assert_eq!(gen_code_singleton_find(42), None);

        gen_code_singleton_init();
        let code_a = vec![b'A'; GENCODE_STRLEN];
        let code_b = vec![b'B'; GENCODE_STRLEN];
        assert_eq!(gen_code_singleton_add(42, &code_a), 0);
        assert_eq!(gen_code_singleton_add(7, &code_b), 0);
        assert_eq!(gen_code_singleton_find(42), Some(code_a.clone()));
        assert_eq!(gen_code_singleton_find(7), Some(code_b));

        assert_eq!(gen_code_singleton_add(42, &[b'Z'; GENCODE_STRLEN]), 0);
        assert_eq!(gen_code_singleton_find(42), Some(code_a));

        gen_code_singleton_fini();
        assert_eq!(gen_code_singleton_find(42), None);
    }

    #[test]
    fn test_s_seq_range_helpers_match_ncbi() {
        let ranges = [
            s_seq_range_new(10, 20),
            s_seq_range_new(30, 40),
            s_seq_range_new(50, 60),
        ];
        assert_eq!(s_seq_range_array_less_than_or_equal(None, 10), -1);
        assert_eq!(s_seq_range_array_less_than_or_equal(Some(&[]), 10), -1);
        assert_eq!(s_seq_range_array_less_than_or_equal(Some(&ranges), 5), 0);
        assert_eq!(s_seq_range_array_less_than_or_equal(Some(&ranges), 35), 1);
        assert_eq!(s_seq_range_array_less_than_or_equal(Some(&ranges), 45), 2);
        assert_eq!(s_seq_range_array_less_than_or_equal(Some(&ranges), 80), 2);
    }

    #[test]
    fn test_blast_sequence_block_setup_and_ranges() {
        let mut slot = None;
        assert_eq!(blast_seq_blk_new(Some(&mut slot)), 0);
        let mut blk = slot.unwrap();
        assert_eq!(
            blast_seq_blk_set_sequence(Some(&mut blk), &[NULLB, 1, 2, 3], 3),
            0
        );
        assert_eq!(blk.sequence.as_deref(), Some(&[1, 2, 3][..]));
        assert_eq!(blk.sequence_nomask.as_deref(), Some(&[1, 2, 3][..]));
        assert!(blk.sequence_start_allocated);

        let ranges = [s_seq_range_new(7, 8), s_seq_range_new(9, 10)];
        assert_eq!(
            blast_seq_blk_set_seq_ranges(
                Some(&mut blk),
                Some(&ranges),
                2,
                true,
                SubjectMaskingType::SoftSubjMasking,
            ),
            0
        );
        assert_eq!(
            blk.seq_ranges.as_deref(),
            Some(&[s_seq_range_new(0, 8), s_seq_range_new(9, 3)][..])
        );
        assert_eq!(blk.num_seq_ranges, 2);
        assert!(blk.seq_ranges_allocated);

        s_blast_sequence_blk_free_seq_ranges(&mut blk);
        assert_eq!(blk.seq_ranges, None);
        assert_eq!(blk.num_seq_ranges, 0);
        assert!(!blk.seq_ranges_allocated);
    }

    #[test]
    fn test_blast_sequence_block_copy_clean_free_and_sfree() {
        let mut blk = None;
        assert_eq!(
            blast_set_up_seq_blk_new(&[NULLB, 4, 5], 2, &mut blk, true),
            0
        );
        let mut blk = blk.unwrap();
        blk.oof_sequence = Some(vec![1, 2, 3]);
        blk.oof_sequence_allocated = true;
        blk.seq_ranges = Some(vec![s_seq_range_new(0, 2)]);
        blk.seq_ranges_allocated = true;
        blk.num_seq_ranges = 1;

        let mut copy = None;
        blast_sequence_blk_copy(&mut copy, &blk);
        let copy = copy.unwrap();
        assert_eq!(copy.sequence, blk.sequence);
        assert!(!copy.sequence_start_allocated);
        assert!(!copy.oof_sequence_allocated);
        assert!(!copy.seq_ranges_allocated);

        blast_sequence_blk_clean(Some(&mut blk));
        assert_eq!(blk.sequence, None);
        assert_eq!(blk.oof_sequence, None);
        assert_eq!(blk.seq_ranges, None);
        assert_eq!(blast_sequence_blk_free(Some(blk)), None);

        let mut owned = Some(vec![1u8, 2]);
        sfree(&mut owned);
        assert_eq!(owned, None);
    }

    #[test]
    fn test_blast_target_translation_new_partial_and_oof() {
        let mut blk = None;
        assert_eq!(
            blast_set_up_seq_blk_new(&[NULLB, 1, 8, 4, 1, 2, 4], 6, &mut blk, true),
            0
        );
        let mut blk = blk.unwrap();
        let mut target = None;
        assert_eq!(
            blast_target_translation_new(
                &mut blk,
                &STANDARD_GENETIC_CODE,
                crate::program::TBLASTN,
                false,
                &mut target,
            ),
            0
        );
        let target = target.unwrap();
        assert!(target.partial);
        assert_eq!(target.num_frames, NUM_FRAMES as i32);
        assert_eq!(target.range.as_deref(), Some(&[0; 12][..]));
        assert!(target.translations.iter().all(Option::is_none));
        assert_eq!(blast_target_translation_free(Some(target)), None);

        let mut target = None;
        assert_eq!(
            blast_target_translation_new(
                &mut blk,
                &STANDARD_GENETIC_CODE,
                crate::program::TBLASTN,
                true,
                &mut target,
            ),
            0
        );
        let target = target.unwrap();
        assert!(!target.partial);
        assert_eq!(target.range, None);
        assert!(blk.oof_sequence_allocated);
        // C BLAST_GetAllTranslations builds a mixed-frame sequence over both
        // strands: 2*nucl_length+3 bytes (blast_util.c:1115).
        assert_eq!(
            blk.oof_sequence.as_ref().unwrap().len(),
            2 * blk.length as usize + 3
        );
    }

    #[test]
    fn test_blast_seq_blk_set_compressed_and_compress_blastna() {
        let mut blk = BlastSequenceBlk::default();
        assert_eq!(
            blast_seq_blk_set_compressed_sequence(Some(&mut blk), &[0b0001_1011, 0]),
            0
        );
        assert_eq!(blk.sequence.as_deref(), Some(&[0b0001_1011, 0][..]));
        assert!(blk.sequence_allocated);

        blk.sequence = Some(vec![0, 1, 2, 3, 0]);
        blk.length = 5;
        assert_eq!(blast_compress_blastna_sequence(Some(&mut blk)), 0);
        let compressed = blk.compressed_nuc_seq.unwrap();
        assert_eq!(compressed.offset, 3);
        assert_eq!(
            blk.compressed_nuc_seq_start.as_deref(),
            Some(compressed.as_start_slice())
        );
        assert_eq!(compressed.as_start_slice()[0], 0b00);
        assert_eq!(compressed.as_start_slice()[1], 0b00_01);
        assert_eq!(compressed.as_start_slice()[2], 0b00_01_10);

        let logical = compressed.as_logical_slice();
        assert_eq!(logical[0], 0b00_01_10_11);
        assert_eq!(logical[1], 0b01_10_11_00);
        assert_eq!(logical[2], 0b10_11_00_00);
        assert_eq!(logical[4], 0);
    }

    #[test]
    fn test_b_search_int4_matches_ncbi_lower_bracket() {
        let values = [2, 4, 8, 16, 32];
        assert_eq!(b_search_int4(1, &values), 0);
        assert_eq!(b_search_int4(2, &values), 0);
        assert_eq!(b_search_int4(7, &values), 1);
        assert_eq!(b_search_int4(8, &values), 2);
        assert_eq!(b_search_int4(100, &values), 4);
        assert_eq!(b_search_int4(100, &[]), 0);
    }

    #[test]
    fn test_blast_str_to_upper_matches_ncbi_duplicate() {
        let original = b"Blastp-62*";
        assert_eq!(
            blast_str_to_upper(Some(original)),
            Some(b"BLASTP-62*".to_vec())
        );
        assert_eq!(blast_str_to_upper(Some(b"")), Some(Vec::new()));
        assert_eq!(blast_str_to_upper(None), None);
    }

    #[test]
    fn test_blast_mem_dup_clones_value() {
        let original = vec![1u8, 2, 3];
        let duplicate = blast_mem_dup(Some(&original)).unwrap();
        assert_eq!(duplicate, original);
        assert_ne!(duplicate.as_ptr(), original.as_ptr());
        assert_eq!(blast_mem_dup::<Vec<u8>>(None), None);
    }

    #[test]
    fn test_blast_pack_dna_matches_ncbi_layout() {
        assert_eq!(
            blast_pack_dna(&[0, 1, 2, 3], 4, E_BLAST_ENCODING_NUCLEOTIDE),
            (0, Some(vec![0b0001_1011, 0]))
        );
        assert_eq!(
            blast_pack_dna(&[3, 2, 1], 3, E_BLAST_ENCODING_NUCLEOTIDE),
            (0, Some(vec![0b1110_0100 | 3]))
        );
        assert_eq!(
            blast_pack_dna(&[1, 2, 4, 8], 4, E_BLAST_ENCODING_NCBI4NA),
            (0, Some(vec![0b0001_1011, 0]))
        );
    }

    #[test]
    fn test_s_blast_get_translation_table_matches_ncbi2na_codons() {
        let table = s_blast_get_translation_table(Some(&STANDARD_GENETIC_CODE), false).unwrap();
        let rc_table = s_blast_get_translation_table(Some(&STANDARD_GENETIC_CODE), true).unwrap();

        assert_eq!(
            table[0b00_11_10],
            translate_codon(0, 3, 2, &STANDARD_GENETIC_CODE)
        );
        assert_eq!(
            rc_table[0b10_11_00],
            translate_codon(3, 0, 1, &STANDARD_GENETIC_CODE)
        );
        assert_eq!(s_blast_get_translation_table(None, false), None);
    }

    #[test]
    fn test_blast_translate_compressed_sequence_matches_ncbi4na_path() {
        let ncbi4na = vec![1u8, 2, 4, 8, 1, 2, 4, 8, 1, 2, 4];
        let (packed_status, packed) =
            blast_pack_dna(&ncbi4na, ncbi4na.len(), E_BLAST_ENCODING_NCBI4NA);
        assert_eq!(packed_status, 0);
        let packed = packed.unwrap();
        let table = s_blast_get_translation_table(Some(&STANDARD_GENETIC_CODE), false).unwrap();
        let rc_table = s_blast_get_translation_table(Some(&STANDARD_GENETIC_CODE), true).unwrap();
        let (flat, offsets) =
            blast_get_all_translations(&ncbi4na, ncbi4na.len(), &STANDARD_GENETIC_CODE);

        for (context, frame) in [1, 2, 3, -1, -2, -3].into_iter().enumerate() {
            let mut prot = vec![0u8; ncbi4na.len() / CODON_LENGTH + 3];
            let translation = if frame > 0 { &table } else { &rc_table };
            let len = blast_translate_compressed_sequence(
                translation,
                ncbi4na.len(),
                Some(&packed),
                frame,
                Some(&mut prot),
            );
            let begin = offsets[context] as usize;
            let end = offsets[context + 1] as usize;
            assert_eq!(&prot[..=len], &flat[begin..end], "frame {frame}");
        }
        let mut prot = vec![99u8; 4];
        assert_eq!(
            blast_translate_compressed_sequence(&table, ncbi4na.len(), None, 1, Some(&mut prot)),
            0
        );
        assert_eq!(
            blast_translate_compressed_sequence(&table, ncbi4na.len(), Some(&packed), 1, None),
            0
        );
    }

    #[test]
    fn test_blast_translate_compressed_sequence_matches_ncbi4na_all_tails_and_fast_loop() {
        let table = s_blast_get_translation_table(Some(&STANDARD_GENETIC_CODE), false).unwrap();
        let rc_table = s_blast_get_translation_table(Some(&STANDARD_GENETIC_CODE), true).unwrap();
        let bases = [1u8, 2, 4, 8];

        for len in 3..96usize {
            let ncbi4na: Vec<u8> = (0..len).map(|i| bases[(i * 7 + len) & 3]).collect();
            let (packed_status, packed) =
                blast_pack_dna(&ncbi4na, ncbi4na.len(), E_BLAST_ENCODING_NCBI4NA);
            assert_eq!(packed_status, 0);
            let packed = packed.unwrap();
            let (flat, offsets) =
                blast_get_all_translations(&ncbi4na, ncbi4na.len(), &STANDARD_GENETIC_CODE);

            for (context, frame) in [1, 2, 3, -1, -2, -3].into_iter().enumerate() {
                let mut prot = vec![0u8; ncbi4na.len() + 3];
                let translation = if frame > 0 { &table } else { &rc_table };
                let got_len = blast_translate_compressed_sequence(
                    translation,
                    ncbi4na.len(),
                    Some(&packed),
                    frame,
                    Some(&mut prot),
                );
                let begin = offsets[context] as usize;
                let end = offsets[context + 1] as usize;
                assert_eq!(
                    &prot[..=got_len],
                    &flat[begin..end],
                    "len {len} frame {frame}"
                );
            }
        }
    }

    #[test]
    fn test_blast_get_partial_translation_single_frame() {
        let seq = vec![1u8, 8, 4, 1, 2, 4, 8];
        let partial = blast_get_partial_translation(
            &seq,
            seq.len(),
            1,
            &STANDARD_GENETIC_CODE,
            true,
            true,
            false,
        )
        .unwrap();
        let mut expected = vec![NULLB; seq.len() / CODON_LENGTH + 2];
        let rev = Vec::new();
        let expected_len = blast_get_translation(
            &seq,
            &rev,
            seq.len(),
            1,
            &mut expected,
            &STANDARD_GENETIC_CODE,
        );
        assert_eq!(partial.protein_length, Some(expected_len));
        assert_eq!(
            partial.translation_buffer.as_deref(),
            Some(expected.as_slice())
        );
        assert_eq!(partial.mixed_seq, None);
    }

    #[test]
    fn test_blast_get_partial_translation_mixed_frames() {
        let seq = vec![1u8, 2, 4, 8, 1, 2, 4, 8, 1];
        let partial = blast_get_partial_translation(
            &seq,
            seq.len(),
            -1,
            &STANDARD_GENETIC_CODE,
            true,
            true,
            true,
        )
        .unwrap();
        assert_eq!(partial.protein_length, Some(seq.len()));
        let mixed = partial.mixed_seq.unwrap();
        assert_eq!(mixed.len(), seq.len() + 1);

        let rev = get_reverse_nucl_sequence(&seq, seq.len());
        let mut frame_buffers = Vec::new();
        for frame in [-1, -2, -3] {
            let mut buffer = vec![NULLB; seq.len() / CODON_LENGTH + 2];
            blast_get_translation(
                &seq,
                &rev,
                seq.len(),
                frame,
                &mut buffer,
                &STANDARD_GENETIC_CODE,
            );
            frame_buffers.push(buffer);
        }
        for (index, &residue) in mixed.iter().enumerate() {
            assert_eq!(
                residue,
                frame_buffers[index % CODON_LENGTH][index / CODON_LENGTH]
            );
        }
    }

    #[test]
    fn test_blast_get_translated_protein_length_matches_ncbi() {
        assert_eq!(blast_get_translated_protein_length(0, 0), 0);
        assert_eq!(blast_get_translated_protein_length(1, 1), 0);
        assert_eq!(blast_get_translated_protein_length(3, 0), 1);
        assert_eq!(blast_get_translated_protein_length(3, 1), 0);
        assert_eq!(blast_get_translated_protein_length(10, 5), 2);
    }

    #[test]
    fn test_translate_codon_met() {
        // ATG = methionine (M = 12 in NCBIstdaa)
        // A=0, T=3, G=2 → idx = 0*16 + 3*4 + 2 = 14
        assert_eq!(translate_codon(0, 3, 2, &STANDARD_GENETIC_CODE), 12);
    }

    #[test]
    fn test_translate_codon_stop() {
        // TAA = stop.
        // T=3, A=0, A=0 → idx = 3*16 + 0*4 + 0 = 48
        assert_eq!(
            translate_codon(3, 0, 0, &STANDARD_GENETIC_CODE),
            crate::encoding::NCBISTDAA_STOP
        );
    }

    #[test]
    fn test_six_frame() {
        // ACGTACGTAC (10 bases) in NCBI4na: A=1, C=2, G=4, T=8
        let seq = vec![1u8, 2, 4, 8, 1, 2, 4, 8, 1, 2];
        let (_buf, offsets) = blast_get_all_translations(&seq, seq.len(), &STANDARD_GENETIC_CODE);
        // Frame +1 covers ACG TAC GTA → 3 residues.
        assert_eq!(offsets[1] - offsets[0] - 1, 3);
        // Frame names are produced in order 1, 2, 3, -1, -2, -3.
        assert_eq!(blast_context_to_frame(0), 1);
        assert_eq!(blast_context_to_frame(3), -1);
    }

    #[test]
    fn test_protein_to_nuc_coords_forward_frame1() {
        // 30bp sequence, frame +1, protein pos 2..5 (0-based)
        // nuc_start = 2*3 + 0 + 1 = 7, nuc_end = 5*3 + 0 = 15
        let (start, end) = protein_to_nuc_coords(2, 5, 1, 30);
        assert_eq!(start, 7);
        assert_eq!(end, 15);
    }

    #[test]
    fn test_protein_to_nuc_coords_forward_frame2() {
        // frame +2: offset = 1
        // nuc_start = 0*3 + 1 + 1 = 2, nuc_end = 3*3 + 1 = 10
        let (start, end) = protein_to_nuc_coords(0, 3, 2, 30);
        assert_eq!(start, 2);
        assert_eq!(end, 10);
    }

    #[test]
    fn test_protein_to_nuc_coords_forward_frame3() {
        // frame +3: offset = 2
        // nuc_start = 1*3 + 2 + 1 = 6, nuc_end = 4*3 + 2 = 14
        let (start, end) = protein_to_nuc_coords(1, 4, 3, 30);
        assert_eq!(start, 6);
        assert_eq!(end, 14);
    }

    #[test]
    fn test_protein_to_nuc_coords_reverse_frame1() {
        // 30bp sequence, frame -1, protein pos 0..3 (0-based on rev-comp)
        // BLAST displays minus-strand intervals as high→low.
        // rc_start=0, rc_end=8, so the forward coordinates are 30..22.
        let (start, end) = protein_to_nuc_coords(0, 3, -1, 30);
        assert_eq!(start, 30);
        assert_eq!(end, 22);
    }

    #[test]
    fn test_protein_to_nuc_coords_reverse_frame2() {
        // 30bp, frame -2, protein pos 1..4
        // rc_start = 4, rc_end = 12, so the forward coordinates are 26..18.
        let (start, end) = protein_to_nuc_coords(1, 4, -2, 30);
        assert_eq!(start, 26);
        assert_eq!(end, 18);
    }

    #[test]
    fn test_protein_to_oriented_nuc_coords_forward_frame2() {
        let (start, end) = protein_to_oriented_nuc_coords(0, 3, 2);
        assert_eq!((start, end), (1, 10));
    }

    #[test]
    fn test_protein_to_oriented_nuc_coords_reverse_frame3() {
        let (start, end) = protein_to_oriented_nuc_coords(2, 5, -3);
        assert_eq!((start, end), (8, 17));
    }

    #[test]
    fn test_translate_codon_all_standard_amino_acids() {
        // Verify a selection of codons from the standard genetic code
        // GCN = Ala (1): GCA = idx 2*16+1*4+0 = 36
        assert_eq!(translate_codon(2, 1, 0, &STANDARD_GENETIC_CODE), 1); // GCA -> A
                                                                         // TGC = Cys (3): idx = 3*16+2*4+1 = 57
        assert_eq!(translate_codon(3, 2, 1, &STANDARD_GENETIC_CODE), 3); // TGC -> C
                                                                         // GAT = Asp (4): idx = 2*16+0*4+3 = 35
        assert_eq!(translate_codon(2, 0, 3, &STANDARD_GENETIC_CODE), 4); // GAT -> D
                                                                         // TTT = Phe (6): idx = 3*16+3*4+3 = 63
        assert_eq!(translate_codon(3, 3, 3, &STANDARD_GENETIC_CODE), 6); // TTT -> F
                                                                         // GGT = Gly (7): idx = 2*16+2*4+3 = 43
        assert_eq!(translate_codon(2, 2, 3, &STANDARD_GENETIC_CODE), 7); // GGT -> G
                                                                         // TGG = Trp (20): idx = 3*16+2*4+2 = 58
        assert_eq!(translate_codon(3, 2, 2, &STANDARD_GENETIC_CODE), 20); // TGG -> W
    }

    #[test]
    fn test_translate_codon_stop_codons() {
        // TAA: T=3, A=0, A=0 -> idx=48
        assert_eq!(
            translate_codon(3, 0, 0, &STANDARD_GENETIC_CODE),
            crate::encoding::NCBISTDAA_STOP
        );
        // TAG: T=3, A=0, G=2 -> idx=50
        assert_eq!(
            translate_codon(3, 0, 2, &STANDARD_GENETIC_CODE),
            crate::encoding::NCBISTDAA_STOP
        );
        // TGA: T=3, G=2, A=0 -> idx=56
        assert_eq!(
            translate_codon(3, 2, 0, &STANDARD_GENETIC_CODE),
            crate::encoding::NCBISTDAA_STOP
        );
    }

    #[test]
    fn test_s_codon_to_aa_unambiguous() {
        // ATG -> M (12): NCBI4na A=1, T=8, G=4
        assert_eq!(s_codon_to_aa([1, 8, 4], &STANDARD_GENETIC_CODE), 12);
        // TAA -> *.
        assert_eq!(
            s_codon_to_aa([8, 1, 1], &STANDARD_GENETIC_CODE),
            crate::encoding::NCBISTDAA_STOP
        );
        // GCT -> A (1)
        assert_eq!(s_codon_to_aa([4, 2, 8], &STANDARD_GENETIC_CODE), 1);
    }

    #[test]
    fn test_s_codon_to_aa_ambiguity_resolves_to_x() {
        // CTN: CTA=L, CTC=L, CTG=L, CTT=L — all leucine, so resolves to L (11).
        // CTN with N=15
        assert_eq!(s_codon_to_aa([2, 8, 15], &STANDARD_GENETIC_CODE), 11);

        // ATN: ATA=I, ATC=I, ATG=M, ATT=I — mixed, must return X.
        assert_eq!(
            s_codon_to_aa([1, 8, 15], &STANDARD_GENETIC_CODE),
            crate::encoding::NCBISTDAA_X
        );

        // YTA (Y=C|T=10): CTA=L, TTA=L — both L, so resolves to L (11).
        assert_eq!(s_codon_to_aa([10, 8, 1], &STANDARD_GENETIC_CODE), 11);

        // RTA (R=A|G=5): ATA=I, GTA=V — mixed, X.
        assert_eq!(
            s_codon_to_aa([5, 8, 1], &STANDARD_GENETIC_CODE),
            crate::encoding::NCBISTDAA_X
        );
    }

    #[test]
    fn test_s_codon_to_aa_n_in_first_position_mixed() {
        // NCT (N=15): all four first-base options. ACT=T, CCT=P, GCT=A, TCT=S.
        // Mixed → X.
        assert_eq!(
            s_codon_to_aa([15, 2, 8], &STANDARD_GENETIC_CODE),
            crate::encoding::NCBISTDAA_X
        );
    }

    #[test]
    fn test_s_codon_to_aa_malformed_falls_through_like_c() {
        // NCBI `s_CodonToAA` only special-cases FENCE_SENTRY when a byte is
        // >15; otherwise the mapping loop still consumes whatever low ncbi4na
        // bits are set. 200 has only the T bit set among the low four bits,
        // so AT(200) follows ATT -> I.
        assert_eq!(s_codon_to_aa([1, 8, 200], &STANDARD_GENETIC_CODE), 9);
        // 16 has no low ncbi4na bits set, so the mapping loop never assigns
        // an amino acid and returns the C initial value, NULLB.
        assert_eq!(s_codon_to_aa([16, 8, 4], &STANDARD_GENETIC_CODE), NULLB);
        // 24 has the low T bit set in addition to high malformed bits, so it
        // falls through as T.
        assert_eq!(
            s_codon_to_aa([24, 1, 1], &STANDARD_GENETIC_CODE),
            crate::encoding::NCBISTDAA_STOP
        );
    }

    #[test]
    fn test_s_codon_to_aa_preserves_fence_sentry() {
        assert_eq!(
            s_codon_to_aa([1, FENCE_SENTRY, 4], &STANDARD_GENETIC_CODE),
            FENCE_SENTRY
        );
    }

    #[test]
    fn test_get_reverse_nucl_sequence_layout() {
        // NCBI4na: A=1, C=2, G=4, T=8 → ACGT
        let seq = vec![1u8, 2, 4, 8];
        let rev = get_reverse_nucl_sequence(&seq, 4);
        // Length+2 with NULLB sentinels at [0] and [length+1]; reverse-complement
        // of ACGT is ACGT (palindrome), so rev[1..=4] should be 1, 2, 4, 8.
        assert_eq!(rev.len(), 6);
        assert_eq!(rev[0], NULLB);
        assert_eq!(rev[5], NULLB);
        assert_eq!(&rev[1..=4], &[1u8, 2, 4, 8]);

        // Non-palindrome: AAAC → GTTT in reverse-complement.
        let seq = vec![1u8, 1, 1, 2];
        let rev = get_reverse_nucl_sequence(&seq, 4);
        assert_eq!(&rev[1..=4], &[4u8, 8, 8, 8]);
    }

    #[test]
    fn test_blast_get_translation_forward_and_reverse() {
        // ATG GCT (NCBI4na) → M A in frame +1.
        let seq: Vec<u8> = vec![1, 8, 4, 4, 2, 8];
        let rev = get_reverse_nucl_sequence(&seq, seq.len());
        let mut prot = vec![0u8; seq.len() / CODON_LENGTH + 2];
        let n = blast_get_translation(&seq, &rev, seq.len(), 1, &mut prot, &STANDARD_GENETIC_CODE);
        assert_eq!(n, 2);
        assert_eq!(prot[0], NULLB);
        assert_eq!(prot[1], 12); // Met
        assert_eq!(prot[2], 1); // Ala
        assert_eq!(prot[3], NULLB);

        // Frame -1 reads the reverse-complement of ATGGCT = AGCCAT, codon AGC = S(17)
        // and CAT = H(8).
        let n_rev =
            blast_get_translation(&seq, &rev, seq.len(), -1, &mut prot, &STANDARD_GENETIC_CODE);
        assert_eq!(n_rev, 2);
        assert_eq!(prot[1], 17); // Ser
        assert_eq!(prot[2], 8); // His
    }

    #[test]
    fn test_blast_get_all_translations_offsets_and_frames() {
        let seq: Vec<u8> = vec![1, 8, 4, 4, 2, 8]; // ATGGCT
        let (buf, offsets) = blast_get_all_translations(&seq, seq.len(), &STANDARD_GENETIC_CODE);
        // Six contexts, frame ∈ {1,2,3,-1,-2,-3}; lengths in residues are
        // 2, 1, 1, 2, 1, 1 for ATGGCT (frame +1: ATG GCT → M A; frame +2:
        // TGG → W; frame +3: GGC → G; reverse-complement is AGCCAT,
        // similarly).
        for ctx in 0..NUM_FRAMES {
            let begin = (offsets[ctx] + 1) as usize; // skip leading NULLB sentinel
            let end = offsets[ctx + 1] as usize;
            assert!(begin <= end);
            // Trailing NULLB at offsets[ctx+1] - 1.
            assert!(end == 0 || buf[end] == NULLB || buf[end] == 0);
        }
        // Frame +1 has 2 residues.
        assert_eq!(offsets[1] - offsets[0] - 1, 2);
    }

    #[test]
    fn test_blast_frame_to_context_all_branches() {
        use crate::program::{BLASTN, BLASTP, BLASTX, TBLASTN};
        // Translated: ±1..±3 → 0..5 (frame-1 / 2-frame).
        for (frame, ctx) in [(1, 0), (2, 1), (3, 2), (-1, 3), (-2, 4), (-3, 5)] {
            assert_eq!(blast_frame_to_context(frame, BLASTX), ctx);
            assert_eq!(blast_frame_to_context(frame, TBLASTN), ctx);
        }
        // Nucleotide: +1 → 0, -1 → 1.
        assert_eq!(blast_frame_to_context(1, BLASTN), 0);
        assert_eq!(blast_frame_to_context(-1, BLASTN), 1);
        // Protein: frame 0 → 0.
        assert_eq!(blast_frame_to_context(0, BLASTP), 0);
    }

    #[test]
    fn test_blast_all_translations_mixed_seq_matches_c_layout() {
        let seq: Vec<u8> = vec![1, 8, 4, 4, 2, 8]; // ATGGCT
        let nucl_length = seq.len();
        let (buf, offsets) = blast_get_all_translations(&seq, nucl_length, &STANDARD_GENETIC_CODE);
        let mixed = blast_all_translations_mixed_seq(&buf, &offsets, nucl_length);
        // C: malloc(2*nucl_length+3), final byte is the trailing NULLB.
        assert_eq!(mixed.len(), 2 * nucl_length + 3);
        assert_eq!(*mixed.last().unwrap(), NULLB);
        // Reproduce C's interleave directly and compare.
        let mut expected = vec![NULLB; 2 * nucl_length + 3];
        let mut out = 0usize;
        let mut index = 0usize;
        while index < NUM_FRAMES {
            for i in 0..=nucl_length {
                expected[out] =
                    buf[offsets[index + (i % CODON_LENGTH)] as usize + (i / CODON_LENGTH)];
                out += 1;
            }
            index += CODON_LENGTH;
        }
        assert_eq!(mixed, expected);
    }

    #[test]
    fn test_blast_create_mixed_frame_dna_translation_interleaves_frames() {
        use crate::queryinfo::{ContextInfo, QueryInfo};
        let ctx = |query_offset: i32, query_length: i32| ContextInfo {
            query_offset,
            query_length,
            eff_searchsp: 0,
            length_adjustment: 0,
            query_index: 0,
            frame: 0,
            is_valid: query_length > 0,
            segment_flags: 0,
        };
        // Three forward-frame contexts at offsets 0/3/6 with lengths 2/2/1.
        let qinfo = QueryInfo {
            num_queries: 1,
            contexts: vec![ctx(0, 2), ctx(3, 2), ctx(6, 1)],
            max_length: 2,
            min_length: 1,
        };
        // seq_buf_len: last ctx offset 6 + len 1 + 2 = 9.
        assert_eq!(qinfo.seq_buf_len(), 9);

        let mut blk = BlastSequenceBlk::default();
        blk.sequence = Some(vec![10, 11, 12, 13, 14, 15, 16, 17, 18, 19]);
        assert_eq!(
            blast_create_mixed_frame_dna_translation(&mut blk, &qinfo),
            0
        );
        assert!(blk.oof_sequence_allocated);
        // 3 leading NULLBs, then interleave f0[0],f1[0],f2[0],f0[1],f1[1] (f2
        // exhausted at offset 1), then a trailing NULLB. Buffer is seq_buf_len+1.
        assert_eq!(
            blk.oof_sequence.as_deref(),
            Some(&[NULLB, NULLB, NULLB, 10, 13, 16, 11, 14, NULLB, NULLB][..])
        );
    }

    #[test]
    fn test_genetic_code_2_vertebrate_mito() {
        // AGA -> * instead of R
        assert_eq!(
            translate_codon(0, 2, 0, &GENETIC_CODE_2),
            crate::encoding::NCBISTDAA_STOP
        );
        // AGG -> * instead of R
        assert_eq!(
            translate_codon(0, 2, 2, &GENETIC_CODE_2),
            crate::encoding::NCBISTDAA_STOP
        );
        // ATA -> M (12) instead of I
        assert_eq!(translate_codon(0, 3, 0, &GENETIC_CODE_2), 12);
        // TGA -> W (20) instead of * (stop)
        assert_eq!(translate_codon(3, 2, 0, &GENETIC_CODE_2), 20);
    }

    #[test]
    fn test_genetic_code_6_ciliate() {
        // TAA -> Q (15) instead of * (stop)
        assert_eq!(translate_codon(3, 0, 0, &GENETIC_CODE_6), 15);
        // TAG -> Q (15) instead of * (stop)
        assert_eq!(translate_codon(3, 0, 2, &GENETIC_CODE_6), 15);
        // TGA should still be stop
        assert_eq!(
            translate_codon(3, 2, 0, &GENETIC_CODE_6),
            crate::encoding::NCBISTDAA_STOP
        );
    }

    #[test]
    fn test_lookup_genetic_code_standard() {
        let code = lookup_genetic_code(1);
        assert_eq!(code as *const _, &STANDARD_GENETIC_CODE as *const _);
        // Code 11 (bacterial) uses same table as standard
        let code11 = lookup_genetic_code(11);
        assert_eq!(code11 as *const _, &STANDARD_GENETIC_CODE as *const _);
    }

    #[test]
    fn test_lookup_genetic_code_unknown_returns_standard() {
        let code = lookup_genetic_code(99);
        assert_eq!(code as *const _, &STANDARD_GENETIC_CODE as *const _);
    }

    #[test]
    fn test_lookup_utility_math_and_hash_helpers_match_c_rules() {
        assert_eq!(iexp(3, 0), 1);
        assert_eq!(iexp(0, 5), 0);
        assert_eq!(iexp(2, 10), 1024);
        assert_eq!(iexp(-2, 3), -8);

        assert_eq!(ilog2(0), 0);
        assert_eq!(ilog2(1), 0);
        assert_eq!(ilog2(2), 1);
        assert_eq!(ilog2(3), 1);
        assert_eq!(ilog2(16), 4);

        let mut a = vec![0; 4];
        let mut output = vec![0; 8];
        let mut cursor = 0;
        fkm(&mut a, 3, 2, &mut output, &mut cursor, None);
        assert_eq!(cursor, 8);
        assert_eq!(output, vec![0, 0, 0, 1, 0, 1, 1, 1]);

        let mut translated = vec![0; 4];
        let mut cursor = 0;
        fkm_output(&[0, 0, 1], 2, 2, &mut translated, &mut cursor, Some(b"AC"));
        assert_eq!(&translated[..2], b"AC");

        assert_eq!(s_popcount(0), 0);
        assert_eq!(s_popcount(0xffff_ffff), 32);
        assert_eq!(s_popcount(0x1010_1010), 4);

        let seq = [1u8, 2, 3, 4];
        let expected = {
            let mut hash = 2_166_136_261u32;
            for base in seq {
                hash = hash.wrapping_mul(16_777_619);
                hash ^= base as u32;
            }
            hash & 0xffff
        };
        assert_eq!(fnv_hash(&seq, 0xffff), expected);
        assert_eq!(fnv_hash(&seq[..2], 0), 0);
    }

    #[test]
    fn test_sparse_uint1_array_indexes_values_and_free() {
        let bitfield = [0b1011u32, 0b1000u32];
        let mut array = blast_sparse_uint1_array_new(Some(&bitfield), 64).expect("sparse array");
        assert_eq!(array.length, 2);
        assert_eq!(array.counts, vec![3, 4]);
        assert_eq!(array.num_elements, 4);

        assert_eq!(blast_sparse_uint1_array_get_index(Some(&array), 0), 0);
        assert_eq!(blast_sparse_uint1_array_get_index(Some(&array), 1), 1);
        assert_eq!(blast_sparse_uint1_array_get_index(Some(&array), 3), 2);
        assert_eq!(blast_sparse_uint1_array_get_index(Some(&array), 35), 3);
        assert_eq!(blast_sparse_uint1_array_get_index(Some(&array), 64), -1);
        assert_eq!(blast_sparse_uint1_array_get_index(None, 0), -1);

        *blast_sparse_uint1_array_get_element(Some(&mut array), 1).expect("element") = 7;
        assert_eq!(array.values[1], 7);
        assert!(blast_sparse_uint1_array_get_element(Some(&mut array), 64).is_none());
        assert!(blast_sparse_uint1_array_new(None, 64).is_none());
        assert!(blast_sparse_uint1_array_free(Some(array)).is_none());
    }

    #[test]
    fn test_six_frame_translation_short_seq() {
        // Sequence shorter than 3 bases (NCBI4na: A=1, C=2): no codons in any frame.
        let seq = vec![1u8, 2]; // AC
        let (_buf, offsets) = blast_get_all_translations(&seq, seq.len(), &STANDARD_GENETIC_CODE);
        for ctx in 0..NUM_FRAMES {
            assert_eq!(
                offsets[ctx + 1] - offsets[ctx],
                1,
                "expected zero residues in ctx {ctx}"
            );
        }
    }

    #[test]
    fn test_six_frame_translation_exact_codon() {
        // Exactly 3 bases (NCBI4na): ATG = Met. A=1, T=8, G=4.
        let seq = vec![1u8, 8, 4];
        let (buf, offsets) = blast_get_all_translations(&seq, seq.len(), &STANDARD_GENETIC_CODE);
        // Frame +1 has 1 residue (Met = 12).
        assert_eq!(offsets[1] - offsets[0] - 1, 1);
        assert_eq!(buf[(offsets[0] + 1) as usize], 12);
        // Frames +2 and +3 have 0 residues.
        assert_eq!(offsets[2] - offsets[1] - 1, 0);
        assert_eq!(offsets[3] - offsets[2] - 1, 0);
    }

    #[test]
    fn test_six_frame_frame_numbers() {
        // Frame numbering follows BLAST_ContextToFrame for blastx-style programs:
        // contexts 0..5 → frames 1, 2, 3, -1, -2, -3.
        assert_eq!(blast_context_to_frame(0), 1);
        assert_eq!(blast_context_to_frame(1), 2);
        assert_eq!(blast_context_to_frame(2), 3);
        assert_eq!(blast_context_to_frame(3), -1);
        assert_eq!(blast_context_to_frame(4), -2);
        assert_eq!(blast_context_to_frame(5), -3);
    }

    #[test]
    fn test_protein_to_nuc_coords_forward_frame1_start() {
        // First protein position, frame +1
        let (start, end) = protein_to_nuc_coords(0, 1, 1, 30);
        // nuc_start = 0*3 + 0 + 1 = 1, nuc_end = 1*3 + 0 = 3
        assert_eq!(start, 1);
        assert_eq!(end, 3);
    }

    #[test]
    fn test_protein_to_nuc_coords_reverse_frame3() {
        // 30bp, frame -3, protein pos 0..2
        // rc_start = 2, rc_end = 7, so the forward coordinates are 28..23.
        let (start, end) = protein_to_nuc_coords(0, 2, -3, 30);
        assert_eq!(start, 28);
        assert_eq!(end, 23);
    }

    #[test]
    fn test_genetic_code_4_mold_mito() {
        // TGA -> W (20) instead of stop
        assert_eq!(translate_codon(3, 2, 0, &GENETIC_CODE_4), 20);
        // Everything else same as standard - ATG still Met
        assert_eq!(translate_codon(0, 3, 2, &GENETIC_CODE_4), 12);
    }

    #[test]
    fn test_genetic_code_10_euplotid() {
        // TGA -> C (3) instead of stop
        assert_eq!(translate_codon(3, 2, 0, &GENETIC_CODE_10), 3);
    }

    #[test]
    fn test_genetic_code_25_sr1() {
        // TGA -> G (7) instead of stop
        assert_eq!(translate_codon(3, 2, 0, &GENETIC_CODE_25), 7);
    }

    #[test]
    fn test_genetic_code_28_condylostoma() {
        // TAA/TAG -> Q and TGA -> W.
        assert_eq!(translate_codon(3, 0, 0, lookup_genetic_code(28)), 15);
        assert_eq!(translate_codon(3, 0, 2, lookup_genetic_code(28)), 15);
        assert_eq!(translate_codon(3, 2, 0, lookup_genetic_code(28)), 20);
    }
}
