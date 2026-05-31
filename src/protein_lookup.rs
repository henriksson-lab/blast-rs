//! Protein word lookup table and scan for blastp/blastx/tblastn.
//!
//! Replaces the O(n^2) brute-force approach with the standard BLAST
//! neighborhood-word lookup table: for each w-mer in the query, all
//! amino-acid words scoring >= threshold (against the scoring matrix)
//! are hashed into a table keyed by word index. During scanning, each
//! subject w-mer is hashed once, checked against a presence vector for
//! fast rejection, and only on a PV hit are backbone entries examined
//! and ungapped extensions triggered.

#[cfg(test)]
use crate::encoding::encode_ncbistdaa_sequence;
use crate::encoding::ncbistdaa_to_aminoacid_char;
use crate::matrix::AA_SIZE;
use crate::protein::{blast_get_start_for_gapped_alignment, protein_gapped_align};
use crate::pssm::Pssm;

/// Result of a protein hit after extension.
#[derive(Debug, Clone)]
pub struct ProteinHit {
    pub query_start: usize,
    pub query_end: usize,
    pub subject_start: usize,
    pub subject_end: usize,
    pub score: i32,
    pub num_ident: i32,
    pub align_length: i32,
    pub mismatches: i32,
    pub gap_opens: i32,
    pub qseq: Option<String>,
    pub sseq: Option<String>,
    /// Pre-rounding scaled DP score from comp_adjust path.
    /// `None` outside comp_adjust. NCBI feeds this scaled value into Spouge with a
    /// matching scaled Lambda; we use it so the rounding-to-display-units doesn't
    /// distort the e-value (`s_HSPListNormalizeScores` happens *after* e-value calc).
    pub scaled_score: Option<i32>,
    /// Mirror of NCBI's `BlastHSP.query.gapped_start` and `subject.gapped_start`
    /// (`blast_hits.h:62`): the seed (q, s) used when running the gapped DP for
    /// this HSP. NCBI's `s_RedoOneAlignment` (`blast_kappa.c:1924-1939`) uses
    /// this seed — NOT the alignment's q_start — when re-running the rescaled
    /// gapped DP under composition adjustment. Default 0 keeps backward
    /// compatibility for callers that haven't computed this yet.
    pub gapped_start_q: usize,
    pub gapped_start_s: usize,
}

/// blast-rs: ProteinHit comparator modeled on score-ordered HSP sorting; not a
/// direct NCBI C port.
/// Same tie-breaker as `hspstream::score_compare_hsps` / `search::score_compare_search_hsps`
/// but over `ProteinHit`'s `usize` offset fields.
fn score_compare_protein_hits(a: &ProteinHit, b: &ProteinHit) -> std::cmp::Ordering {
    b.score
        .cmp(&a.score)
        .then_with(|| a.subject_start.cmp(&b.subject_start))
        .then_with(|| b.subject_end.cmp(&a.subject_end))
        .then_with(|| a.query_start.cmp(&b.query_start))
        .then_with(|| b.query_end.cmp(&a.query_end))
}

fn blast_init_protein_hits_sort_by_score(hits: &mut [ProteinHit]) {
    hits.sort_unstable_by(score_compare_protein_hits);
}

/// Protein word lookup table.
///
/// Uses a CSR (compressed sparse row) layout for cache-friendly access:
/// `data[offsets[hash]..offsets[hash+1]]` holds query offsets whose
/// neighborhood contains the word that hashes to `hash`.
/// `pv` is a presence-vector bit array for fast rejection.
/// Maximum hits stored inline per backbone cell (matches NCBI AA_HITS_PER_CELL).
const HITS_PER_CELL: usize = 3;

/// One cell of the thick backbone. Stores up to HITS_PER_CELL query offsets
/// inline (no pointer chase). Overflows to a separate array.
#[derive(Clone, Copy)]
struct BackboneCell {
    num_used: u16,
    /// If num_used <= HITS_PER_CELL: inline entries.
    /// If num_used > HITS_PER_CELL: entries[0] is overflow cursor into `overflow`.
    entries: [i32; HITS_PER_CELL],
}

pub struct ProteinLookupTable {
    word_size: usize,
    /// Thick backbone: each cell stores up to 3 hits inline.
    backbone: Vec<BackboneCell>,
    /// Overflow array for cells with > HITS_PER_CELL hits.
    overflow: Vec<i32>,
    pv: Vec<u64>,
    /// Compressed-alphabet lookup table, populated only for `word_size >= 5`.
    /// NCBI's `BLAST_FillLookupTableOptions` (`blast_options.c`) switches to
    /// `eCompressedAaLookupTable` for `word_size > 5`; the dense neighborhood
    /// backbone (2^(word_size*5) cells with 28-letter neighbor enumeration)
    /// is infeasible for `word_size >= 5`. When present, the two-hit/one-hit
    /// scanners drive the identical extension logic from the compressed scan's
    /// `(query_offset, subject_offset)` stream instead of the dense backbone.
    compressed: Option<BlastCompressedAaLookupTable>,
}

const COMPRESSED_HITS_PER_BACKBONE_CELL: usize = 4;
const COMPRESSED_HITS_PER_OVERFLOW_CELL: usize = 4;
const COMPRESSED_HITS_CELL_MASK: i32 = 0x03;
const COMPRESSED_OVERFLOW_CELLS_IN_BANK: i32 = 209_710;
const COMPRESSED_OVERFLOW_MAX_BANKS: i32 = 1024;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CompressedOverflowCell {
    pub next: Option<usize>,
    pub query_offsets: [i32; COMPRESSED_HITS_PER_OVERFLOW_CELL],
}

impl Default for CompressedOverflowCell {
    fn default() -> Self {
        Self {
            next: None,
            query_offsets: [0; COMPRESSED_HITS_PER_OVERFLOW_CELL],
        }
    }
}

/// Backbone cell for the compressed-alphabet lookup table.
///
/// NCBI `CompressedLookupBackboneCell` (`blast_aalookup.h:221-231`) is a 24-byte
/// struct: `Int4 num_used; Int4 query_offset; union { Int4 query_offsets[4];
/// struct { Int4 query_offsets[2]; CompressedOverflowCell* head; } overflow_list;
/// } payload;`. The two union arms share the same 16-byte `payload`, selected by
/// `num_used`:
///   * arm A (`num_used <= 5`): `payload` holds up to 4 inline query offsets
///     (`payload.query_offsets[0..4]`, for entries 2..5);
///   * arm B (`num_used > 5`): `payload[0..2]` hold `overflow_list.query_offsets`
///     and `payload[2..4]` hold the 8-byte `head` pointer.
///
/// blast-rs mirrors this with a single 16-byte `payload: [i32; 4]`, packing the
/// overflow-cell index (instead of a raw pointer) into `payload[2..4]` for arm B.
/// This collapses the previous 48-byte representation to NCBI's 24 bytes, halving
/// the per-query backbone calloc/zero cost (~11.39M cells per query).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CompressedLookupBackboneCell {
    pub num_used: i32,
    pub query_offset: i32,
    /// The 16-byte union payload; interpretation selected by `num_used`
    /// exactly as NCBI's `union payload` (see struct docs).
    pub payload: [i32; COMPRESSED_HITS_PER_BACKBONE_CELL],
}

impl Default for CompressedLookupBackboneCell {
    fn default() -> Self {
        Self {
            num_used: 0,
            query_offset: 0,
            payload: [0; COMPRESSED_HITS_PER_BACKBONE_CELL],
        }
    }
}

impl CompressedLookupBackboneCell {
    /// Arm A: the inline query offsets (`payload.query_offsets`, NCBI
    /// `blast_aalookup.h:227`). Valid when `num_used <= 5`; entry slot `i`
    /// (0-based) holds the offset for the `(i+2)`-th hit.
    #[inline]
    fn inline_offsets(&self) -> &[i32; COMPRESSED_HITS_PER_BACKBONE_CELL] {
        &self.payload
    }

    #[inline]
    fn inline_offsets_mut(&mut self) -> &mut [i32; COMPRESSED_HITS_PER_BACKBONE_CELL] {
        &mut self.payload
    }

    /// Arm B: the two query offsets spilled to the backbone when the overflow
    /// list is created (NCBI `overflow_list.query_offsets`, `blast_aalookup.h:213`).
    #[inline]
    fn overflow_offsets(&self) -> [i32; COMPRESSED_HITS_PER_BACKBONE_CELL - 2] {
        [self.payload[0], self.payload[1]]
    }

    #[inline]
    fn set_overflow_offsets(&mut self, a: i32, b: i32) {
        self.payload[0] = a;
        self.payload[1] = b;
    }

    /// Arm B: the overflow-list head index, packed into the 8 bytes that hold
    /// the `CompressedOverflowCell* head` pointer in NCBI (`blast_aalookup.h:217`).
    /// We store `index + 1` so that all-zero payload (the calloc'd state) reads
    /// back as `None`, matching NCBI's NULL head.
    #[inline]
    fn overflow_head(&self) -> Option<usize> {
        let lo = self.payload[2] as u32 as u64;
        let hi = self.payload[3] as u32 as u64;
        let packed = lo | (hi << 32);
        if packed == 0 {
            None
        } else {
            Some((packed - 1) as usize)
        }
    }

    #[inline]
    fn set_overflow_head(&mut self, head: Option<usize>) {
        let packed: u64 = match head {
            None => 0,
            Some(idx) => (idx as u64) + 1,
        };
        self.payload[2] = packed as u32 as i32;
        self.payload[3] = (packed >> 32) as u32 as i32;
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BlastCompressedAaLookupTable {
    pub threshold: i32,
    pub word_length: usize,
    pub alphabet_size: usize,
    pub compressed_alphabet_size: usize,
    pub longest_chain: i32,
    pub backbone: Vec<CompressedLookupBackboneCell>,
    pub overflow_cells: Vec<CompressedOverflowCell>,
    pub curr_overflow_cell: i32,
    pub curr_overflow_bank: i32,
    pub pv: Vec<u32>,
    pub pv_array_bts: usize,
    pub reciprocal_alphabet_size: u64,
    pub compress_table: Vec<u8>,
    pub scaled_compress_table: Vec<i32>,
    pub neighbor_matches: i32,
    pub exact_matches: i32,
}

impl BlastCompressedAaLookupTable {
    pub fn new(word_length: usize, compressed_alphabet_size: usize, backbone_size: usize) -> Self {
        let mut compress_table = vec![u8::MAX; AA_SIZE];
        for (idx, value) in compress_table.iter_mut().enumerate() {
            *value = idx as u8;
        }
        let scaled_compress_table =
            build_scaled_compress_table(word_length, compressed_alphabet_size, &compress_table);
        Self {
            threshold: 0,
            word_length,
            alphabet_size: AA_SIZE,
            compressed_alphabet_size,
            longest_chain: 0,
            backbone: vec![CompressedLookupBackboneCell::default(); backbone_size],
            overflow_cells: Vec::new(),
            curr_overflow_cell: COMPRESSED_OVERFLOW_CELLS_IN_BANK,
            curr_overflow_bank: -1,
            pv: Vec::new(),
            pv_array_bts: crate::stat::PV_ARRAY_BTS,
            reciprocal_alphabet_size: compressed_reciprocal_alphabet_size(compressed_alphabet_size),
            compress_table,
            scaled_compress_table,
            neighbor_matches: 0,
            exact_matches: 0,
        }
    }
}

fn compressed_reciprocal_alphabet_size(compressed_alphabet_size: usize) -> u64 {
    if compressed_alphabet_size == 0 {
        0
    } else {
        (1u64 << 32).div_ceil(compressed_alphabet_size as u64)
    }
}

fn compressed_scale(word_length: usize, compressed_alphabet_size: usize) -> i32 {
    let mut scale = 1i32;
    for _ in 1..word_length {
        scale = scale.saturating_mul(compressed_alphabet_size as i32);
    }
    scale
}

fn build_scaled_compress_table(
    word_length: usize,
    compressed_alphabet_size: usize,
    compress_table: &[u8],
) -> Vec<i32> {
    let scale = compressed_scale(word_length, compressed_alphabet_size);
    compress_table
        .iter()
        .map(|&compressed| {
            if compressed as usize >= compressed_alphabet_size {
                -1
            } else {
                compressed as i32 * scale
            }
        })
        .collect()
}

fn compressed_reciprocal_preshift(index: usize, reciprocal: u64) -> usize {
    (((index as u64) * reciprocal) >> 32) as usize
}

fn compressed_prefix_index(
    lookup: &BlastCompressedAaLookupTable,
    subject: &[u8],
    start: usize,
    width: usize,
) -> Option<usize> {
    let mut index = 0usize;
    let mut factor = 1usize;
    for &letter in subject[start..start + width].iter() {
        let compressed = *lookup.compress_table.get(letter as usize)?;
        if compressed as usize >= lookup.compressed_alphabet_size {
            return None;
        }
        index += compressed as usize * factor;
        factor *= lookup.compressed_alphabet_size;
    }
    Some(index)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AaLookupBoneType {
    Backbone,
    Smallbone,
}

struct NeighborInfo<'a> {
    query_word: &'a [u8],
    subject_word: &'a mut [u8],
    alphabet_size: usize,
    word_size: usize,
    matrix: &'a [[i32; AA_SIZE]; AA_SIZE],
    row_max: &'a [i32; AA_SIZE],
    offset_list: &'a [i32],
    threshold: i32,
}

struct PssmNeighborInfo<'a> {
    pssm: &'a Pssm,
    subject_word: &'a mut [u8],
    alphabet_size: usize,
    word_size: usize,
    row_max: &'a [i32],
    offset: i32,
    threshold: i32,
    query_offset: usize,
}

pub struct CompressedNeighborInfo<'a> {
    pub compressed_alphabet_size: usize,
    pub wordsize: usize,
    pub matrix: &'a [[i32; AA_SIZE]; AA_SIZE],
    pub row_max: [i32; AA_SIZE],
    pub threshold: i32,
    pub matrix_sorted: [[i32; AA_SIZE]; AA_SIZE],
    pub matrix_sorted_char: [[u8; AA_SIZE]; AA_SIZE],
}

impl ProteinLookupTable {
    /// Build the lookup table from a query sequence.
    ///
    /// * `query` - NCBIstdaa-encoded query
    /// * `word_size` - typically 3 for blastp
    /// * `matrix` - scoring matrix (AA_SIZE x AA_SIZE)
    /// * `threshold` - minimum neighborhood word score (typically 11)
    pub fn build(
        query: &[u8],
        word_size: usize,
        matrix: &[[i32; AA_SIZE]; AA_SIZE],
        threshold: f64,
    ) -> Self {
        // NCBI switches to the compressed-alphabet lookup table for large word
        // sizes (`eCompressedAaLookupTable`), where the dense 28-letter
        // neighborhood backbone explodes (ws6 -> ~22 GB / billions of
        // insertions). The dense path below is left BYTE-FOR-BYTE unchanged for
        // `word_size <= 4`; only `word_size >= 5` takes the compressed branch.
        // NCBI BLAST_FillLookupTableOptions (blast_options.c:1165) routes to the
        // compressed table only for word_size > 5; ws5 stays on the dense
        // 28-letter table. (See the build below.)
        if word_size > 5 {
            return Self::build_compressed(query, word_size, matrix, threshold);
        }
        let alphabet_size = AA_SIZE;
        // NCBI sizes the AA backbone to the highest valid BLASTAA index plus
        // one; the rolling scan mask is still the full charsize bit mask.
        let table_size = aa_lookup_backbone_size(word_size);
        let mut backbone: Vec<Vec<i32>> = vec![Vec::new(); table_size];
        let mut exact_backbone: Vec<Vec<i32>> = vec![Vec::new(); table_size];

        // Precompute per-row maximums for branch-and-bound pruning.
        // row_max[aa] = max score achievable when the query letter is `aa`.
        let mut row_max = [0i32; AA_SIZE];
        for q in 0..AA_SIZE {
            let mut mx = i32::MIN;
            for s in 0..AA_SIZE {
                if matrix[q][s] > mx {
                    mx = matrix[q][s];
                }
            }
            row_max[q] = mx;
        }

        let thresh_i = threshold as i32; // integer threshold for comparison
        if query.len() >= word_size {
            for i in 0..=(query.len() - word_size) {
                let query_word = &query[i..i + word_size];
                if !protein_lookup_word_is_valid(query_word.iter().copied()) {
                    continue;
                }
                let hash = word_hash(query_word, alphabet_size);
                exact_backbone[hash].push(i as i32);
            }

            let mut word_buf = vec![0u8; word_size];
            for offsets in &exact_backbone {
                if offsets.is_empty() {
                    continue;
                }
                let query_offset = offsets[0] as usize;
                let query_word = &query[query_offset..query_offset + word_size];
                s_add_word_hits(
                    &mut backbone,
                    matrix,
                    query_word,
                    offsets,
                    thresh_i,
                    &row_max,
                    &mut word_buf,
                    alphabet_size,
                    word_size,
                );
            }
        }

        // Build presence vector.
        let pv_len = table_size.div_ceil(64);
        let mut pv = vec![0u64; pv_len];
        for (idx, entries) in backbone.iter().enumerate() {
            if !entries.is_empty() {
                pv[idx >> 6] |= 1u64 << (idx & 63);
            }
        }

        // Convert to thick backbone (inline ≤3 hits, overflow for rest).
        // Matches NCBI BLAST+ AaLookupBackboneCell layout.
        let empty_cell = BackboneCell {
            num_used: 0,
            entries: [0; HITS_PER_CELL],
        };
        let mut thick: Vec<BackboneCell> = vec![empty_cell; table_size];
        let mut overflow: Vec<i32> = Vec::new();

        for (idx, entries) in backbone.iter().enumerate() {
            let n = entries.len();
            if n == 0 {
                continue;
            }
            thick[idx].num_used = n as u16;
            if n <= HITS_PER_CELL {
                // Store inline — no pointer chase needed at scan time
                for (i, &val) in entries.iter().enumerate() {
                    thick[idx].entries[i] = val;
                }
            } else {
                // Store overflow cursor in entries[0]
                thick[idx].entries[0] = overflow.len() as i32;
                overflow.extend_from_slice(entries);
            }
        }

        ProteinLookupTable {
            word_size,
            backbone: thick,
            overflow,
            pv,
            compressed: None,
        }
    }

    /// Build the compressed-alphabet lookup table for `word_size` 5/6/7.
    ///
    /// NCBI: `BlastCompressedAaLookupTableNew` (`blast_aalookup.c:1270`). The
    /// 28-letter protein alphabet is collapsed to a 15-letter alphabet
    /// (`s_alphabet15`, for ws 5/6) or 10-letter (`s_alphabet10`, ws 7) so
    /// neighborhood enumeration stays bounded. A custom non-square
    /// `BLASTAA_SIZE x compressed_alphabet_size` score matrix is built, scaled
    /// by `kMatrixScale / ungapped_lambda` (kMatrixScale = 100), and the
    /// neighbor threshold is correspondingly scaled by `kMatrixScale`.
    fn build_compressed(
        query: &[u8],
        word_size: usize,
        matrix: &[[i32; AA_SIZE]; AA_SIZE],
        threshold: f64,
    ) -> Self {
        const K_MATRIX_SCALE: f64 = 100.0;
        let compressed_alphabet_size = if word_size == 7 { 10 } else { 15 };
        let backbone_size = (compressed_alphabet_size as f64).powi(word_size as i32) as usize + 1;

        let mut lookup =
            BlastCompressedAaLookupTable::new(word_size, compressed_alphabet_size, backbone_size);
        lookup.threshold = (K_MATRIX_SCALE * threshold) as i32;
        lookup.alphabet_size = AA_SIZE;
        lookup.reciprocal_alphabet_size =
            compressed_reciprocal_alphabet_size(compressed_alphabet_size);

        // Build the compressed translation table (protein letter -> compressed
        // letter), the reverse table, and the scaled non-square score matrix.
        let matrix_name = matrix_name_for_scoring(matrix);
        let (compress_table, compressed_matrix) =
            build_compressed_alphabet(matrix_name, compressed_alphabet_size, K_MATRIX_SCALE);
        lookup.compress_table = compress_table;
        lookup.scaled_compress_table = build_scaled_compress_table(
            word_size,
            compressed_alphabet_size,
            &lookup.compress_table,
        );

        // Index the query's neighboring words into the compressed backbone.
        s_compressed_add_neighboring_words(&mut lookup, &compressed_matrix, query, None);
        s_compressed_lookup_finalize(&mut lookup);

        ProteinLookupTable {
            word_size,
            backbone: Vec::new(),
            overflow: Vec::new(),
            pv: Vec::new(),
            compressed: Some(lookup),
        }
    }
}

/// blast-rs: identify the named scoring matrix backing a raw score matrix so the
/// compressed lookup can fetch the matrix-specific ungapped lambda and frequency
/// ratios (NCBI threads the matrix name via `sbp->name`). Falls back to BLOSUM62.
fn matrix_name_for_scoring(matrix: &[[i32; AA_SIZE]; AA_SIZE]) -> &'static str {
    let candidates: [(&'static str, &'static [[i32; AA_SIZE]; AA_SIZE]); 8] = [
        ("BLOSUM62", &crate::matrix::BLOSUM62),
        ("BLOSUM45", &crate::matrix::BLOSUM45),
        ("BLOSUM50", &crate::matrix::BLOSUM50),
        ("BLOSUM80", &crate::matrix::BLOSUM80),
        ("BLOSUM90", &crate::matrix::BLOSUM90),
        ("PAM30", &crate::matrix::PAM30),
        ("PAM70", &crate::matrix::PAM70),
        ("PAM250", &crate::matrix::PAM250),
    ];
    for (name, candidate) in candidates {
        if candidate == matrix {
            return name;
        }
    }
    "BLOSUM62"
}

/// NCBI: `SCompressedAlphabetNew` + `s_BuildCompressedTranslation` +
/// `s_BuildCompressedScoreMatrix` + `s_GetCompressedProbs` (`blast_stat.c`).
///
/// Returns `(compress_table, compressed_matrix)` where `compress_table[aa]` maps
/// a BLASTAA letter to its compressed letter (or `compressed_alphabet_size` for
/// letters with no group), and `compressed_matrix[q][s]` is the score for query
/// BLASTAA letter `q` against compressed letter `s`, scaled by
/// `matrix_scale_factor / ungapped_lambda` and rounded — exactly as NCBI's
/// lookup-only score matrix. Only the first `compressed_alphabet_size` columns
/// of `compressed_matrix` are meaningful.
fn build_compressed_alphabet(
    matrix_name: &str,
    compressed_alphabet_size: usize,
    matrix_scale_factor: f64,
) -> (Vec<u8>, [[i32; AA_SIZE]; AA_SIZE]) {
    // 23-to-15 (SE_B(14)) and 23-to-10 (SE-V(10)) compressed alphabets, letter
    // groups separated by spaces, in decreasing combined residue frequency
    // (NCBI `s_alphabet15` / `s_alphabet10`).
    const ALPHABET15: &str = "ST IJV LM KR EQZ A G BD P N F Y H C W";
    const ALPHABET10: &str = "IJLMV AST BDENZ KQR G FY P H C W";
    let alphabet_string = if compressed_alphabet_size == 10 {
        ALPHABET10
    } else {
        ALPHABET15
    };

    // s_BuildCompressedTranslation: parse the grouping string into a translation
    // table and a reverse (compressed letter -> list of BLASTAA letters) table.
    let mut compress_table = vec![compressed_alphabet_size as u8; AA_SIZE];
    let mut rev_table: Vec<Vec<usize>> = vec![Vec::new(); compressed_alphabet_size];
    let mut compressed_letter = 0usize;
    for ch in alphabet_string.chars() {
        if ch.is_whitespace() {
            compressed_letter += 1;
        } else if ch.is_ascii_alphabetic() {
            let aa = crate::encoding::AMINOACID_TO_NCBISTDAA[ch as usize & 0x7f] as usize;
            if compressed_letter < compressed_alphabet_size && aa < AA_SIZE {
                compress_table[aa] = compressed_letter as u8;
                rev_table[compressed_letter].push(aa);
            }
        }
    }

    // s_BuildCompressedScoreMatrix.
    let lambda = crate::stat::protein_ungapped_kbp_for_matrix(matrix_name).lambda;
    let scale = if lambda > 0.0 {
        matrix_scale_factor / lambda
    } else {
        matrix_scale_factor
    };
    let std_freqs = crate::matrix::get_matrix_freq_ratios(matrix_name)
        .unwrap_or_else(|| crate::matrix::get_matrix_freq_ratios("BLOSUM62").unwrap());
    let prob = crate::stat::protein_std_freq_ncbistdaa();

    // s_GetCompressedProbs: P(aa | compressed letter) = prob[aa] / sum of group.
    let mut compressed_prob = [0.0f64; AA_SIZE];
    for group in rev_table.iter() {
        let prob_sum: f64 = group.iter().map(|&aa| prob[aa]).sum();
        if prob_sum > 0.0 {
            for &aa in group {
                compressed_prob[aa] = prob[aa] / prob_sum;
            }
        }
    }

    // BLAST_SCORE_MIN (NCBI) is -32768.
    const BLAST_SCORE_MIN: f64 = -32768.0;
    let min_freq = BLAST_SCORE_MIN / scale;
    let mut compressed_matrix = [[0i32; AA_SIZE]; AA_SIZE];
    for q in 0..AA_SIZE {
        for (s, group) in rev_table.iter().enumerate() {
            let mut val = 0.0f64;
            for &aa in group {
                val += std_freqs[q][aa] * compressed_prob[aa];
            }
            let val = if val < 1e-8 { min_freq } else { val.ln() };
            compressed_matrix[q][s] = crate::math::blast_nint(val * scale) as i32;
        }
    }

    (compress_table, compressed_matrix)
}

/// NCBI: `BlastAaLookupTableNew` backbone sizing (`blast_aalookup.c`).
/// naming: Local helper for one constructor substep, not the full C constructor.
fn aa_lookup_backbone_size(word_size: usize) -> usize {
    let mut backbone_size = 0usize;
    for i in 0..word_size {
        backbone_size |= (AA_SIZE - 1) << (i * CHARSIZE);
    }
    backbone_size + 1
}

/// NCBI: BlastAaLookupFinalize (blast_aalookup.c).
/// Conservative Rust port of NCBI `BlastAaLookupFinalize`
/// (`blast_aalookup.c:267`).
///
/// Audit caveat: C finalizes from a temporary `thin_backbone` allocation into
/// either `AaLookupBackboneCell` or `AaLookupSmallboneCell` storage, then frees
/// the thin rows. Rust `ProteinLookupTable::build` already stores the finalized
/// thick-backbone representation, so this routine recomputes the presence
/// vector and compacts overflow storage in place. `Smallbone` is accepted for
/// C-name parity, but uses the same `u16` hit-count representation as the
/// existing Rust scanner.
pub fn blast_aa_lookup_finalize(
    lookup: Option<&mut ProteinLookupTable>,
    bone_type: AaLookupBoneType,
) -> i32 {
    let Some(lookup) = lookup else {
        return -1;
    };
    let _smallbone = matches!(bone_type, AaLookupBoneType::Smallbone);

    let pv_len = lookup.backbone.len().div_ceil(64);
    let mut pv = vec![0u64; pv_len];
    let mut compact_overflow = Vec::new();

    for (index, cell) in lookup.backbone.iter_mut().enumerate() {
        let num_used = cell.num_used as usize;
        if num_used == 0 {
            cell.entries = [0; HITS_PER_CELL];
            continue;
        }

        pv[index >> 6] |= 1u64 << (index & 63);
        if num_used <= HITS_PER_CELL {
            for slot in cell.entries.iter_mut().skip(num_used) {
                *slot = 0;
            }
            continue;
        }

        let old_start = cell.entries[0].max(0) as usize;
        let old_end = old_start
            .saturating_add(num_used)
            .min(lookup.overflow.len());
        let new_start = compact_overflow.len();
        compact_overflow.extend_from_slice(&lookup.overflow[old_start..old_end]);
        cell.entries = [0; HITS_PER_CELL];
        cell.entries[0] = new_start as i32;
    }

    lookup.pv = pv;
    lookup.overflow = compact_overflow;
    0
}

/// NCBI: s_AddWordHits (blast_aalookup.c).
fn s_add_word_hits(
    backbone: &mut [Vec<i32>],
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
    query_word: &[u8],
    offset_list: &[i32],
    threshold: i32,
    row_max: &[i32; AA_SIZE],
    subject_word: &mut [u8],
    alphabet_size: usize,
    word_size: usize,
) {
    let self_score: i32 = query_word
        .iter()
        .map(|&aa| matrix[aa as usize][aa as usize])
        .sum();

    if threshold == 0 || self_score < threshold {
        for &query_offset in offset_list {
            blast_lookup_add_word_hit(backbone, word_size, CHARSIZE, query_word, query_offset);
        }
    }

    if threshold == 0 {
        return;
    }

    let mut info = NeighborInfo {
        query_word,
        subject_word,
        alphabet_size,
        word_size,
        matrix,
        row_max,
        offset_list,
        threshold,
    };

    let mut score = row_max[query_word[0] as usize];
    for &aa in &query_word[1..] {
        score += row_max[aa as usize];
    }

    s_add_word_hits_core(backbone, &mut info, score, 0);
}

/// NCBI: s_AddWordHitsCore (blast_aalookup.c).
fn s_add_word_hits_core(
    backbone: &mut [Vec<i32>],
    info: &mut NeighborInfo<'_>,
    mut score: i32,
    current_pos: usize,
) {
    score -= info.row_max[info.query_word[current_pos] as usize];
    let row = &info.matrix[info.query_word[current_pos] as usize];

    if current_pos == info.word_size - 1 {
        for (aa, &cell_score) in row.iter().take(info.alphabet_size).enumerate() {
            if score + cell_score >= info.threshold {
                info.subject_word[current_pos] = aa as u8;
                for &query_offset in info.offset_list {
                    blast_lookup_add_word_hit(
                        backbone,
                        info.word_size,
                        CHARSIZE,
                        info.subject_word,
                        query_offset,
                    );
                }
            }
        }
        return;
    }

    for (aa, &cell_score) in row.iter().take(info.alphabet_size).enumerate() {
        if score + cell_score >= info.threshold {
            info.subject_word[current_pos] = aa as u8;
            s_add_word_hits_core(backbone, info, score + cell_score, current_pos + 1);
        }
    }
}

/// NCBI: s_AddPSSMNeighboringWords (blast_aalookup.c).
/// Add neighboring words from position-specific rows, matching
/// `s_AddPSSMNeighboringWords` in NCBI's protein lookup builder.
///
/// Audit caveat: C receives a raw `BlastSeqLoc` list and `Int4**` row pointers.
/// The Rust API receives a `Pssm` plus optional half-open query ranges, then
/// follows the same row-max pruning and word insertion logic over Vec-backed
/// thin-backbone cells.
pub fn s_add_pssm_neighboring_words(
    backbone: &mut [Vec<i32>],
    pssm: &Pssm,
    word_size: usize,
    threshold: i32,
    query_bias: i32,
    ranges: Option<&[(usize, usize)]>,
) -> i32 {
    if word_size == 0 || pssm.length < word_size {
        return 0;
    }

    let table_size = aa_lookup_backbone_size(word_size);
    if backbone.len() < table_size {
        return -1;
    }

    let mut total = 0;
    let mut subject_word = vec![0u8; word_size];

    if let Some(ranges) = ranges {
        for &(start, end) in ranges {
            total += s_add_pssm_neighboring_range(
                backbone,
                pssm,
                word_size,
                threshold,
                query_bias,
                start,
                end,
                &mut subject_word,
            );
        }
    } else {
        total += s_add_pssm_neighboring_range(
            backbone,
            pssm,
            word_size,
            threshold,
            query_bias,
            0,
            pssm.length,
            &mut subject_word,
        );
    }

    total
}

/// blast-rs: Range wrapper around the C s_AddPSSMNeighboringWords helper; not
/// a direct NCBI C port.
fn s_add_pssm_neighboring_range(
    backbone: &mut [Vec<i32>],
    pssm: &Pssm,
    word_size: usize,
    threshold: i32,
    query_bias: i32,
    start: usize,
    end: usize,
    subject_word: &mut [u8],
) -> i32 {
    let start = start.min(pssm.length);
    let end = end.min(pssm.length);
    if end < start || end - start < word_size {
        return 0;
    }

    let last_offset = end - word_size;
    let mut total = 0;
    for query_offset in start..=last_offset {
        let Some(query_offset_i32) = i32::try_from(query_offset).ok() else {
            return -1;
        };
        let Some(offset) = query_offset_i32.checked_add(query_bias) else {
            return -1;
        };

        let mut row_max = vec![i32::MIN; word_size];
        for (i, max_score) in row_max.iter_mut().enumerate() {
            *max_score = pssm.scores[query_offset + i]
                .iter()
                .take(AA_SIZE)
                .copied()
                .max()
                .unwrap_or(i32::MIN);
        }

        total += s_add_pssm_word_hits(
            backbone,
            pssm,
            query_offset,
            word_size,
            threshold,
            offset,
            &row_max,
            subject_word,
        );
    }

    total
}

/// NCBI: s_AddPSSMWordHits (blast_aalookup.c).
fn s_add_pssm_word_hits(
    backbone: &mut [Vec<i32>],
    pssm: &Pssm,
    query_offset: usize,
    word_size: usize,
    threshold: i32,
    offset: i32,
    row_max: &[i32],
    subject_word: &mut [u8],
) -> i32 {
    let score = row_max.iter().take(word_size).sum();
    let mut info = PssmNeighborInfo {
        pssm,
        subject_word,
        alphabet_size: AA_SIZE,
        word_size,
        row_max,
        offset,
        threshold,
        query_offset,
    };

    s_add_pssm_word_hits_core(backbone, &mut info, score, 0)
}

/// NCBI: s_AddPSSMWordHitsCore (blast_aalookup.c).
fn s_add_pssm_word_hits_core(
    backbone: &mut [Vec<i32>],
    info: &mut PssmNeighborInfo<'_>,
    mut score: i32,
    current_pos: usize,
) -> i32 {
    score -= info.row_max[current_pos];
    let row = &info.pssm.scores[info.query_offset + current_pos];

    if current_pos == info.word_size - 1 {
        let mut added = 0;
        for (aa, &cell_score) in row.iter().take(info.alphabet_size).enumerate() {
            if score + cell_score >= info.threshold {
                info.subject_word[current_pos] = aa as u8;
                blast_lookup_add_word_hit(
                    backbone,
                    info.word_size,
                    CHARSIZE,
                    info.subject_word,
                    info.offset,
                );
                added += 1;
            }
        }
        return added;
    }

    let mut added = 0;
    for (aa, &cell_score) in row.iter().take(info.alphabet_size).enumerate() {
        if score + cell_score >= info.threshold {
            info.subject_word[current_pos] = aa as u8;
            added += s_add_pssm_word_hits_core(backbone, info, score + cell_score, current_pos + 1);
        }
    }
    added
}

/// NCBI: s_CompressedListGetNewCell (blast_aalookup.c).
/// Allocate the next compressed-overflow cell, matching
/// `s_CompressedListGetNewCell`'s bank/cursor state.
///
/// Audit caveat: C stores banks of raw cells and returns a raw pointer. Rust
/// stores the same cells in one Vec and returns the stable cell index used by
/// the compressed backbone's linked-list fields.
pub fn s_compressed_list_get_new_cell(lookup: &mut BlastCompressedAaLookupTable) -> Option<usize> {
    if lookup.curr_overflow_cell == COMPRESSED_OVERFLOW_CELLS_IN_BANK {
        let bank_idx = lookup.curr_overflow_bank + 1;
        if bank_idx >= COMPRESSED_OVERFLOW_MAX_BANKS {
            return None;
        }
        lookup.curr_overflow_bank = bank_idx;
        lookup.curr_overflow_cell = 0;
    }

    let cell_index = lookup.overflow_cells.len();
    lookup
        .overflow_cells
        .push(CompressedOverflowCell::default());
    lookup.curr_overflow_cell += 1;
    Some(cell_index)
}

/// NCBI: s_CompressedLookupAddWordHit (blast_aalookup.c).
/// Add one query offset to a compressed-alphabet lookup cell, following
/// `s_CompressedLookupAddWordHit`'s inline-to-overflow transition.
///
/// The first five hits are stored in the backbone representation. On the sixth
/// hit, the last two inline payload offsets plus the new offset are moved to an
/// overflow cell; subsequent hits fill the overflow head and allocate another
/// head cell whenever the current one reaches the C cell capacity of four.
pub fn s_compressed_lookup_add_word_hit(
    lookup: &mut BlastCompressedAaLookupTable,
    index: usize,
    query_offset: i32,
) -> i32 {
    if index >= lookup.backbone.len() {
        return -1;
    }

    // NCBI: s_CompressedLookupAddWordHit (blast_aalookup.c:790-847).
    let num_entries = lookup.backbone[index].num_used;
    match num_entries {
        // C: case 0 — backbone_cell->query_offset = query_offset; (:792)
        0 => {
            lookup.backbone[index].query_offset = query_offset;
        }
        // C: case 1..4 — payload.query_offsets[num_entries-1] = query_offset; (:798)
        1..=4 => {
            lookup.backbone[index].inline_offsets_mut()[(num_entries - 1) as usize] = query_offset;
        }
        // C: case 5 — spill last two inline offsets + new offset into a fresh
        // overflow cell; move first two inline offsets into the overflow arm. (:800-822)
        5 => {
            let payload = *lookup.backbone[index].inline_offsets();
            let tmp_offset_0 = payload[0];
            let tmp_offset_1 = payload[1];
            let Some(new_cell) = s_compressed_list_get_new_cell(lookup) else {
                return -1;
            };
            lookup.overflow_cells[new_cell].next = None;
            lookup.overflow_cells[new_cell].query_offsets[0] = payload[2];
            lookup.overflow_cells[new_cell].query_offsets[1] = payload[3];
            lookup.overflow_cells[new_cell].query_offsets[2] = query_offset;
            lookup.backbone[index].set_overflow_offsets(tmp_offset_0, tmp_offset_1);
            lookup.backbone[index].set_overflow_head(Some(new_cell));
        }
        // C: default — continue existing overflow list, allocating a new head
        // cell whenever the current one is full. (:825-845)
        _ => {
            let cell_offset = ((num_entries - 3) & COMPRESSED_HITS_CELL_MASK) as usize;
            if cell_offset == 0 {
                let old_head = lookup.backbone[index].overflow_head();
                let Some(new_cell) = s_compressed_list_get_new_cell(lookup) else {
                    return -1;
                };
                lookup.overflow_cells[new_cell].next = old_head;
                lookup.backbone[index].set_overflow_head(Some(new_cell));
            }

            let Some(head) = lookup.backbone[index].overflow_head() else {
                return -1;
            };
            lookup.overflow_cells[head].query_offsets[cell_offset] = query_offset;
        }
    }

    lookup.backbone[index].num_used += 1;
    0
}

/// NCBI: s_CompressedLookupAddEncoded (blast_aalookup.c).
/// Add one already-compressed word to the compressed lookup table, preserving
/// the fixed-radix index formulas and packed backbone/overflow insertion
/// layout used by the C implementation.
pub fn s_compressed_lookup_add_encoded(
    lookup: &mut BlastCompressedAaLookupTable,
    word: &[u8],
    query_offset: i32,
) -> i32 {
    const W7P1: [usize; 10] = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90];
    const W7P2: [usize; 10] = [0, 100, 200, 300, 400, 500, 600, 700, 800, 900];
    const W7P3: [usize; 10] = [0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000];
    const W7P4: [usize; 10] = [
        0, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000,
    ];
    const W7P5: [usize; 10] = [
        0, 100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000,
    ];
    const W7P6: [usize; 10] = [
        0, 1000000, 2000000, 3000000, 4000000, 5000000, 6000000, 7000000, 8000000, 9000000,
    ];
    const W6P1: [usize; 15] = [
        0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180, 195, 210,
    ];
    const W6P2: [usize; 15] = [
        0, 225, 450, 675, 900, 1125, 1350, 1575, 1800, 2025, 2250, 2475, 2700, 2925, 3150,
    ];
    const W6P3: [usize; 15] = [
        0, 3375, 6750, 10125, 13500, 16875, 20250, 23625, 27000, 30375, 33750, 37125, 40500, 43875,
        47250,
    ];
    const W6P4: [usize; 15] = [
        0, 50625, 101250, 151875, 202500, 253125, 303750, 354375, 405000, 455625, 506250, 556875,
        607500, 658125, 708750,
    ];
    const W6P5: [usize; 15] = [
        0, 759375, 1518750, 2278125, 3037500, 3796875, 4556250, 5315625, 6075000, 6834375, 7593750,
        8353125, 9112500, 9871875, 10631250,
    ];

    let index = match lookup.word_length {
        5 => {
            if word.len() < 5 || word[..5].iter().any(|&c| c as usize >= 15) {
                return -1;
            }
            word[0] as usize
                + W6P1[word[1] as usize]
                + W6P2[word[2] as usize]
                + W6P3[word[3] as usize]
                + W6P4[word[4] as usize]
        }
        6 => {
            if word.len() < 6 || word[..6].iter().any(|&c| c as usize >= 15) {
                return -1;
            }
            word[0] as usize
                + W6P1[word[1] as usize]
                + W6P2[word[2] as usize]
                + W6P3[word[3] as usize]
                + W6P4[word[4] as usize]
                + W6P5[word[5] as usize]
        }
        7 => {
            if word.len() < 7 || word[..7].iter().any(|&c| c as usize >= 10) {
                return -1;
            }
            word[0] as usize
                + W7P1[word[1] as usize]
                + W7P2[word[2] as usize]
                + W7P3[word[3] as usize]
                + W7P4[word[4] as usize]
                + W7P5[word[5] as usize]
                + W7P6[word[6] as usize]
        }
        _ => return -1,
    };

    s_compressed_lookup_add_word_hit(lookup, index, query_offset)
}

/// blast-rs: Unencoded-input adapter for the C s_CompressedLookupAddEncoded
/// helper; not a direct NCBI C port.
fn s_compressed_lookup_add_unencoded(
    lookup: &mut BlastCompressedAaLookupTable,
    word: &[u8],
    query_offset: i32,
) -> i32 {
    if word.len() < lookup.word_length {
        return -1;
    }

    let mut encoded = vec![0u8; lookup.word_length];
    for (dst, &letter) in encoded.iter_mut().zip(word.iter()) {
        let Some(&compressed) = lookup.compress_table.get(letter as usize) else {
            return 0;
        };
        if compressed as usize >= lookup.compressed_alphabet_size {
            return 0;
        }
        *dst = compressed;
    }

    s_compressed_lookup_add_encoded(lookup, &encoded, query_offset)
}

/// NCBI: s_LoadSortedMatrix (blast_aalookup.c).
pub fn s_load_sorted_matrix(info: &mut CompressedNeighborInfo<'_>) {
    for long_char in 0..AA_SIZE {
        let mut pairs: Vec<(i32, u8)> = (0..info.compressed_alphabet_size)
            .map(|short_char| {
                (
                    info.row_max[long_char] - info.matrix[long_char][short_char],
                    short_char as u8,
                )
            })
            .collect();
        pairs.sort_by_key(|&(diff, letter)| (diff, letter));

        for (i, &(_, letter)) in pairs.iter().enumerate() {
            info.matrix_sorted[long_char][i] = info.matrix[long_char][letter as usize];
            info.matrix_sorted_char[long_char][i] = letter;
        }
    }
}

/// blast-rs: Owned setup for NCBI `CompressedNeighborInfo`.
/// Mirrors the row-max and sorted-matrix preparation that the recursive
/// compressed-neighbor enumerator consumes, but returns a Rust value instead
/// of filling caller-owned C storage.
fn compressed_neighbor_info<'a>(
    lookup: &BlastCompressedAaLookupTable,
    matrix: &'a [[i32; AA_SIZE]; AA_SIZE],
) -> CompressedNeighborInfo<'a> {
    let mut info = CompressedNeighborInfo {
        compressed_alphabet_size: lookup.compressed_alphabet_size,
        wordsize: lookup.word_length,
        matrix,
        row_max: [0; AA_SIZE],
        threshold: lookup.threshold,
        matrix_sorted: [[0; AA_SIZE]; AA_SIZE],
        matrix_sorted_char: [[0; AA_SIZE]; AA_SIZE],
    };

    for row in 0..lookup.alphabet_size.min(AA_SIZE) {
        info.row_max[row] = matrix[row][0];
        for col in 1..lookup.compressed_alphabet_size.min(AA_SIZE) {
            info.row_max[row] = info.row_max[row].max(matrix[row][col]);
        }
    }
    s_load_sorted_matrix(&mut info);
    info
}

/// NCBI: s_CompressedAddWordHitsCore (blast_aalookup.c).
/// Recursive compressed-neighbor enumeration, matching
/// `s_CompressedAddWordHitsCore`.
///
/// Audit caveat: C keeps `query_word`, `subject_word`, and lookup pointers
/// inside `CompressedNeighborInfo`; Rust passes those slices explicitly to keep
/// mutable lookup access separate from immutable scoring state.
pub fn s_compressed_add_word_hits_core(
    lookup: &mut BlastCompressedAaLookupTable,
    info: &CompressedNeighborInfo<'_>,
    query_word: &[u8],
    subject_word: &mut [u8],
    query_offset: i32,
    mut score: i32,
    current_pos: usize,
) -> i32 {
    let curr_query_char = query_word[current_pos] as usize;
    score -= info.row_max[curr_query_char];
    let row_sorted = &info.matrix_sorted[curr_query_char];
    let char_sorted = &info.matrix_sorted_char[curr_query_char];

    if current_pos == info.wordsize - 1 {
        let mut added = 0;
        for i in 0..info.compressed_alphabet_size {
            if score + row_sorted[i] < info.threshold {
                break;
            }
            subject_word[current_pos] = char_sorted[i];
            if s_compressed_lookup_add_encoded(lookup, subject_word, query_offset) == 0 {
                lookup.neighbor_matches += 1;
                added += 1;
            }
        }
        return added;
    }

    let mut added = 0;
    for i in 0..info.compressed_alphabet_size {
        if score + row_sorted[i] < info.threshold {
            break;
        }
        subject_word[current_pos] = char_sorted[i];
        added += s_compressed_add_word_hits_core(
            lookup,
            info,
            query_word,
            subject_word,
            query_offset,
            score + row_sorted[i],
            current_pos + 1,
        );
    }
    added
}

/// NCBI: s_CompressedAddWordHits (blast_aalookup.c).
pub fn s_compressed_add_word_hits(
    lookup: &mut BlastCompressedAaLookupTable,
    info: &CompressedNeighborInfo<'_>,
    query: &[u8],
    query_offset: usize,
) -> i32 {
    if lookup.word_length == 0 || query_offset + lookup.word_length > query.len() {
        return -1;
    }

    lookup.exact_matches += 1;
    let query_word = &query[query_offset..query_offset + lookup.word_length];
    let mut score = 0;
    for &letter in query_word {
        let Some(&compressed) = lookup.compress_table.get(letter as usize) else {
            return 0;
        };
        if compressed as usize >= lookup.compressed_alphabet_size {
            return 0;
        }
        score += info.matrix[letter as usize][compressed as usize];
    }

    let mut added = 0;
    if lookup.threshold == 0 || score < lookup.threshold {
        added += s_compressed_lookup_add_unencoded(lookup, query_word, query_offset as i32);
    } else {
        lookup.neighbor_matches -= 1;
    }

    if lookup.threshold == 0 {
        return added;
    }

    let mut subject_word = vec![0u8; lookup.word_length];
    score = query_word
        .iter()
        .map(|&letter| info.row_max[letter as usize])
        .sum();

    added
        + s_compressed_add_word_hits_core(
            lookup,
            info,
            query_word,
            &mut subject_word,
            query_offset as i32,
            score,
            0,
        )
}

/// NCBI: s_CompressedAddNeighboringWords (blast_aalookup.c).
/// Index compressed-alphabet neighboring words over optional half-open query
/// ranges, matching `s_CompressedAddNeighboringWords`'s setup and loop shape.
pub fn s_compressed_add_neighboring_words(
    lookup: &mut BlastCompressedAaLookupTable,
    compressed_matrix: &[[i32; AA_SIZE]; AA_SIZE],
    query: &[u8],
    ranges: Option<&[(usize, usize)]>,
) -> i32 {
    if lookup.word_length == 0 || query.len() < lookup.word_length {
        return 0;
    }

    let info = compressed_neighbor_info(lookup, compressed_matrix);
    let mut total = 0;
    if let Some(ranges) = ranges {
        for &(start, end) in ranges {
            let start = start.min(query.len());
            let end = end.min(query.len());
            if end < start || end - start < lookup.word_length {
                continue;
            }
            for offset in start..=end - lookup.word_length {
                total += s_compressed_add_word_hits(lookup, &info, query, offset);
            }
        }
    } else {
        for offset in 0..=query.len() - lookup.word_length {
            total += s_compressed_add_word_hits(lookup, &info, query, offset);
        }
    }
    total
}

/// NCBI: s_CompressedLookupFinalize (blast_aalookup.c).
/// Complete compressed lookup construction by building the presence vector and
/// longest-chain metadata, matching `s_CompressedLookupFinalize`.
pub fn s_compressed_lookup_finalize(lookup: &mut BlastCompressedAaLookupTable) -> i32 {
    const TARGET_PV_BYTES: usize = 262_144;

    let backbone_size = lookup.backbone.len();
    let occupied = lookup
        .backbone
        .iter()
        .filter(|cell| cell.num_used > 0)
        .count();

    let mut pv_array_bts = crate::stat::PV_ARRAY_BTS;
    if occupied as f64 <= 0.01 * backbone_size as f64 {
        pv_array_bts +=
            crate::util::ilog2((backbone_size / (8 * TARGET_PV_BYTES)) as i64).max(0) as usize;
    }

    let mut pv = vec![0u32; (backbone_size >> pv_array_bts) + 1];
    let mut longest_chain = 0;
    for (index, cell) in lookup.backbone.iter().enumerate() {
        let count = cell.num_used;
        if count > 0 {
            pv[index >> pv_array_bts] |= 1u32 << (index & crate::stat::PV_ARRAY_MASK as usize);
            longest_chain = longest_chain.max(count);
        }
    }

    lookup.pv = pv;
    lookup.pv_array_bts = pv_array_bts;
    lookup.longest_chain = longest_chain;
    0
}

thread_local! {
    /// FIX C: per-thread reusable scan buffer for the compressed word finders.
    /// Mirrors NCBI's single `offset_pairs` array, allocated once per worker and
    /// reused (capacity retained) across every subject (`blast_engine.c:1039`).
    static COMPRESSED_SCAN_BUF: std::cell::Cell<Vec<crate::lookup::OffsetPair>> =
        const { std::cell::Cell::new(Vec::new()) };
}

/// blast-rs: inline-read helper for FIX B. Mirrors the offset-copy block of
/// NCBI's `s_BlastCompressedAaScanSubject` (`blast_aascan.c:345-396`), pushing
/// `(query_offset, s_off)` pairs for one backbone cell directly into `out`
/// without an intermediate per-cell `Vec`. Returns the number of hits pushed.
#[inline]
fn compressed_push_cell_hits(
    lookup: &BlastCompressedAaLookupTable,
    cell: &CompressedLookupBackboneCell,
    s_off: i32,
    out: &mut Vec<crate::lookup::OffsetPair>,
) -> usize {
    let numhits = cell.num_used.max(0) as usize;
    if numhits == 0 {
        return 0;
    }

    // C: dest[0] = backbone_cell->query_offset (:348-349)
    out.push(crate::lookup::OffsetPair {
        query_offset: cell.query_offset,
        subject_offset: s_off,
    });

    if numhits <= COMPRESSED_HITS_PER_BACKBONE_CELL + 1 {
        // C: hits all live in the backbone (:351-357)
        for &q in &cell.inline_offsets()[..numhits - 1] {
            out.push(crate::lookup::OffsetPair {
                query_offset: q,
                subject_offset: s_off,
            });
        }
        return numhits;
    }

    // C: hits are in the backbone cell and in the overflow list (:359-395)
    let overflow = cell.overflow_offsets();
    out.push(crate::lookup::OffsetPair {
        query_offset: overflow[0],
        subject_offset: s_off,
    });
    out.push(crate::lookup::OffsetPair {
        query_offset: overflow[1],
        subject_offset: s_off,
    });

    let first_cell_entries = ((numhits as i32 - 3) & COMPRESSED_HITS_CELL_MASK).max(0) as usize;
    let mut current = cell.overflow_head();
    if let Some(head) = current {
        for &q in &lookup.overflow_cells[head].query_offsets[..first_cell_entries] {
            out.push(crate::lookup::OffsetPair {
                query_offset: q,
                subject_offset: s_off,
            });
        }
        if first_cell_entries != 0 {
            current = lookup.overflow_cells[head].next;
        }
    }
    while let Some(cell_index) = current {
        for &q in &lookup.overflow_cells[cell_index].query_offsets {
            out.push(crate::lookup::OffsetPair {
                query_offset: q,
                subject_offset: s_off,
            });
        }
        current = lookup.overflow_cells[cell_index].next;
    }
    numhits
}

/// NCBI: s_BlastCompressedAaScanSubject (blast_aalookup.c).
/// Scan a subject sequence with a compressed amino-acid lookup table, matching
/// `s_BlastCompressedAaScanSubject`'s PV test and offset-pair copy behavior.
/// Hits are appended into the caller-owned `hits` buffer (FIX C: the buffer is
/// reused across subjects — NCBI reuses one `offset_pairs` array sized once in
/// `blast_engine.c:1039-1040`). Returns the number of pairs written.
///
/// Audit caveat: C uses a reciprocal-multiply rolling index over
/// `scaled_compress_table`; Rust mirrors that rolling index, with slice bounds
/// checks and a rebuilt scaled table if callers changed `compress_table`.
pub fn s_blast_compressed_aa_scan_subject_into(
    lookup: &BlastCompressedAaLookupTable,
    subject: &[u8],
    array_size: usize,
    scan_range: Option<(usize, usize)>,
    hits: &mut Vec<crate::lookup::OffsetPair>,
) -> usize {
    hits.clear();
    if lookup.word_length == 0
        || lookup.compressed_alphabet_size == 0
        || subject.len() < lookup.word_length
        || array_size == 0
    {
        return 0;
    }

    let start = scan_range.map(|range| range.0).unwrap_or(0);
    let end = scan_range
        .map(|range| range.1)
        .unwrap_or(subject.len() - lookup.word_length)
        .min(subject.len() - lookup.word_length);
    if start > end {
        return 0;
    }

    let scale = compressed_scale(lookup.word_length, lookup.compressed_alphabet_size);
    let scaled_table_is_current = lookup.scaled_compress_table.len() == lookup.compress_table.len()
        && lookup
            .compress_table
            .iter()
            .zip(&lookup.scaled_compress_table)
            .all(|(&compressed, &scaled)| {
                let expected = if compressed as usize >= lookup.compressed_alphabet_size {
                    -1
                } else {
                    compressed as i32 * scale
                };
                scaled == expected
            });
    let rebuilt_scaled_compress_table;
    let scaled_compress_table = if scaled_table_is_current {
        &lookup.scaled_compress_table
    } else {
        rebuilt_scaled_compress_table = build_scaled_compress_table(
            lookup.word_length,
            lookup.compressed_alphabet_size,
            &lookup.compress_table,
        );
        &rebuilt_scaled_compress_table
    };
    let reciprocal = if lookup.reciprocal_alphabet_size != 0 {
        lookup.reciprocal_alphabet_size
    } else {
        compressed_reciprocal_alphabet_size(lookup.compressed_alphabet_size)
    };
    let prefix_width = lookup.word_length - 1;
    let mut subject_offset = start;

    'prime: while subject_offset <= end {
        let Some(prefix_index) =
            compressed_prefix_index(lookup, subject, subject_offset, prefix_width)
        else {
            subject_offset += 1;
            continue;
        };
        let mut shifted_index = prefix_index;

        while subject_offset <= end {
            let next_letter = subject[subject_offset + prefix_width] as usize;
            let scaled_compressed = scaled_compress_table
                .get(next_letter)
                .copied()
                .unwrap_or(-1);
            if scaled_compressed < 0 {
                subject_offset += 1;
                continue 'prime;
            }

            let index = shifted_index + scaled_compressed as usize;
            shifted_index = compressed_reciprocal_preshift(index, reciprocal);

            let Some(pv_word) = lookup.pv.get(index >> lookup.pv_array_bts) else {
                subject_offset += 1;
                continue;
            };
            if (pv_word & (1u32 << (index & crate::stat::PV_ARRAY_MASK as usize))) == 0 {
                subject_offset += 1;
                continue;
            }

            let Some(cell) = lookup.backbone.get(index) else {
                subject_offset += 1;
                continue;
            };
            if cell.num_used <= 0 {
                subject_offset += 1;
                continue;
            }
            // C: copy hits only if `numhits <= array_size - totalhits`; on
            // overflow set s_range[1] and return without partial copy
            // (blast_aascan.c:340-405). FIX B: read offsets inline into `hits`,
            // no per-cell Vec.
            let numhits = cell.num_used.max(0) as usize;
            if numhits > array_size.saturating_sub(hits.len()) {
                return hits.len();
            }
            compressed_push_cell_hits(lookup, cell, subject_offset as i32, hits);
            subject_offset += 1;
        }
    }
    hits.len()
}

fn s_blast_compressed_aa_scan_subject_batched<F>(
    lookup: &BlastCompressedAaLookupTable,
    subject: &[u8],
    array_size: usize,
    hits: &mut Vec<crate::lookup::OffsetPair>,
    mut consume_batch: F,
) -> usize
where
    F: FnMut(&[crate::lookup::OffsetPair]),
{
    if lookup.word_length == 0 || subject.len() < lookup.word_length {
        hits.clear();
        return 0;
    }

    let last_pos = subject.len() - lookup.word_length;
    let longest_chain = lookup.longest_chain.max(1) as usize;
    let batch_capacity = array_size.max(longest_chain).max(1);
    let scan_span = (batch_capacity / longest_chain).max(1);
    let mut total_hits = 0usize;
    let mut start = 0usize;

    while start <= last_pos {
        let end = (start + scan_span - 1).min(last_pos);
        let written = s_blast_compressed_aa_scan_subject_into(
            lookup,
            subject,
            batch_capacity,
            Some((start, end)),
            hits,
        );
        total_hits += written;
        consume_batch(hits);
        start = end + 1;
    }

    total_hits
}

/// Bits per residue for hashing (ceil(log2(alphabet_size))).
/// AA_SIZE=28, charsize=5 (rounds up to 32). Matches NCBI BLAST+ `charsize`.
const CHARSIZE: usize = 5;

/// NCBI: BlastLookupAddWordHit (blast_lookup.c).
/// 1-1 port of `BlastLookupAddWordHit` (`blast_lookup.c:33`) over the Rust
/// vector-backed thin backbone.
///
/// NCBI stores each cell as a manually grown `[capacity, used, hits...]` chain.
/// The Rust representation uses `Vec<i32>` per cell; `push` is the exact
/// ownership-safe equivalent of appending `query_offset` after growing the C
/// chain when needed.
pub(crate) fn blast_lookup_add_word_hit(
    backbone: &mut [Vec<i32>],
    wordsize: usize,
    charsize: usize,
    seq: &[u8],
    query_offset: i32,
) {
    if !protein_lookup_word_is_valid(seq.iter().take(wordsize).copied()) {
        return;
    }
    let index = s_compute_table_index(wordsize, charsize, seq);
    if let Some(cell) = backbone.get_mut(index) {
        cell.push(query_offset);
    }
}

/// NCBI: s_ComputeTableIndex (lookup_util.c).
fn s_compute_table_index(wordsize: usize, charsize: usize, word: &[u8]) -> usize {
    let mut index = 0usize;
    for &letter in word.iter().take(wordsize) {
        index = (index << charsize) | letter as usize;
    }
    index
}

fn protein_lookup_word_is_valid<I>(word: I) -> bool
where
    I: IntoIterator<Item = u8>,
{
    word.into_iter().all(|letter| (letter as usize) < AA_SIZE)
}

/// Hash a word using NCBI-style shift+or (matching ComputeTableIndex).
/// hash = (w[0] << (n-1)*charsize) | (w[1] << (n-2)*charsize) | ... | w[n-1]
#[inline]
fn word_hash(word: &[u8], _alphabet_size: usize) -> usize {
    s_compute_table_index(word.len(), CHARSIZE, word)
}

/// NCBI: s_ComputeTableIndexIncremental (lookup_util.c).
#[inline]
#[allow(dead_code)]
fn s_compute_table_index_incremental(prev_hash: usize, new_char: u8, mask: usize) -> usize {
    ((prev_hash << CHARSIZE) | new_char as usize) & mask
}

/// Build a merged presence vector from multiple lookup tables.
///
/// The merged PV is the bitwise OR of all individual PVs. A subject word
/// only needs to be checked against individual queries if the merged PV bit
/// is set — this skips ~90% of positions when queries are diverse.
pub fn merge_pv(tables: &[&ProteinLookupTable]) -> Vec<u64> {
    if tables.is_empty() {
        return Vec::new();
    }
    let pv_len = tables[0].pv.len();
    let mut merged = vec![0u64; pv_len];
    for table in tables {
        for (i, &bits) in table.pv.iter().enumerate() {
            merged[i] |= bits;
        }
    }
    merged
}

fn merged_pv_is_superset(tables: &[&ProteinLookupTable], merged_pv: &[u64]) -> bool {
    if tables.is_empty() {
        return false;
    }
    let pv_len = tables[0].pv.len();
    if merged_pv.len() != pv_len || tables.iter().any(|table| table.pv.len() != pv_len) {
        return false;
    }
    tables.iter().all(|table| {
        table
            .pv
            .iter()
            .zip(merged_pv.iter())
            .all(|(&table_bits, &merged_bits)| table_bits & !merged_bits == 0)
    })
}

fn subject_has_merged_pv_hit(subject: &[u8], word_size: usize, merged_pv: &[u64]) -> bool {
    if subject.len() < word_size || word_size == 0 {
        return false;
    }
    let mask = (1usize << (word_size * CHARSIZE)) - 1;
    let last_pos = subject.len() - word_size;
    let mut hash = word_hash(&subject[0..word_size], 0);
    for s_pos in 0..=last_pos {
        if s_pos > 0 {
            hash = s_compute_table_index_incremental(hash, subject[s_pos + word_size - 1], mask);
        }
        if merged_pv
            .get(hash >> 6)
            .is_some_and(|bits| bits & (1u64 << (hash & 63)) != 0)
        {
            return true;
        }
    }
    false
}

/// blast-rs: Public batch-scan wrapper for one subject against multiple
/// protein queries. This is not an NCBI C entry point; it combines a merged PV
/// prefilter with the translated single-query two-hit scanner below.
///
/// The merged PV is used only as a conservative no-hit prefilter. Once any
/// possible word is present, each query still runs through
/// [`protein_scan_with_table_reuse`] so the NCBI two-hit diagonal state,
/// overlap checks, and extension gates stay identical to the single-query path.
///
/// Returns a `Vec` of `(query_index, Vec<ProteinHit>)` for queries that had hits.
pub fn s_blast_aa_scan_subject(
    queries: &[&[u8]],
    tables: &[&ProteinLookupTable],
    merged_pv: &[u64],
    subject: &[u8],
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
    x_dropoff: i32,
) -> Vec<(usize, Vec<ProteinHit>)> {
    if tables.is_empty() || subject.is_empty() {
        return Vec::new();
    }
    let word_size = tables[0].word_size;
    if merged_pv_is_superset(tables, merged_pv)
        && !subject_has_merged_pv_hit(subject, word_size, merged_pv)
    {
        return Vec::new();
    }

    // Keep batch BLASTP byte-for-byte aligned with the single-query protein
    // scanner. Earlier code used the merged PV but extended every accepted word
    // hit with a simpler one-hit diagonal filter; NCBI's protein path goes
    // through `s_BlastAaWordFinder_TwoHit`, including the diagonal-offset trick
    // and the left-must-reach-first-hit check. Reusing the same translated
    // scanner here trades some batching speed for faithful behavior.
    let mut results = Vec::new();
    let mut diag_buf = Vec::new();
    for (qi, (&query, &table)) in queries.iter().zip(tables.iter()).enumerate() {
        let hits = protein_scan_with_table_reuse(
            query,
            subject,
            matrix,
            table,
            x_dropoff,
            TWO_HIT_WINDOW,
            1,
            &mut diag_buf,
        );
        if !hits.is_empty() {
            results.push((qi, hits));
        }
    }
    results
}

/// blast-rs: Public convenience wrapper that builds a protein lookup table and
/// scans one subject. Keep this Rust API name stable; the translated NCBI
/// scanner is reached through `protein_scan_with_table`.
///
/// Returns a list of `ProteinHit` sorted by descending score.
pub fn blast_choose_protein_scan_subject(
    query: &[u8],
    subject: &[u8],
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
    word_size: usize,
    threshold: f64,
    x_dropoff: i32,
) -> Vec<ProteinHit> {
    if query.len() < word_size || subject.len() < word_size {
        return Vec::new();
    }
    let table = ProteinLookupTable::build(query, word_size, matrix, threshold);
    protein_scan_with_table(query, subject, matrix, &table, x_dropoff)
}

/// blast-rs: Rust convenience wrapper around `s_blast_aa_word_finder_two_hit`
/// that owns the diagonal scratch buffer for callers that do not batch scans.
/// Keep this public API name stable for Rust call sites.
pub fn protein_scan_with_table(
    query: &[u8],
    subject: &[u8],
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
    table: &ProteinLookupTable,
    x_dropoff: i32,
) -> Vec<ProteinHit> {
    let mut diag_buf = Vec::new();
    protein_scan_with_table_reuse(
        query,
        subject,
        matrix,
        table,
        x_dropoff,
        TWO_HIT_WINDOW,
        1,
        &mut diag_buf,
    )
}

/// Two-hit window size (matches NCBI BLAST+ default for protein).
const TWO_HIT_WINDOW: i32 = 40;

/// Inter-query sentinel byte for the concatenated multi-query buffer (see
/// `blastp_batch`). Chosen `>= AA_SIZE` so the lookup builder's
/// `protein_lookup_word_is_valid` skips any word spanning a query boundary, and
/// so ungapped extension stops here (the value never appears in an encoded
/// query, where `encode_ncbistdaa_sequence` only emits 0..=27). This mirrors
/// NCBI's NULLB inter-context sentinel in the concatenated query (NCBI uses 0;
/// 255 is used here to guarantee no interference with the single-query path,
/// whose extension boundary check would otherwise need a never-in-query value).
pub(crate) const PROTEIN_CONCAT_SENTINEL: u8 = 255;

/// Unsupported branches for the local `BlastAaWordFinder` representation.
///
/// Retained for API stability. The `WindowSize` variant is no longer produced:
/// the two-hit scanner now threads an arbitrary two-hit window through
/// (`-window_size`), so non-default windows are represented. No
/// `BlastAaWordFinder` dispatch branch is currently unsupported.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BlastAaWordFinderUnsupported {
    /// Historically reported when a non-default two-hit window was requested.
    /// Now unused — non-default windows are honored by the scanner.
    WindowSize { requested: i32, supported: i32 },
}

/// blast-rs: owned result for the represented protein word-finder side
/// channels that NCBI fills through `BlastInitHitList`/`BlastUngappedStats`.
#[derive(Debug, Clone, Default)]
pub struct BlastAaWordFinderScan {
    pub hits: Vec<ProteinHit>,
    pub total_hits: i32,
    pub hits_extended: i32,
}

#[derive(Debug, Clone)]
struct BlastAaInitHitRecord {
    q_start: i32,
    s_start: i32,
    q_off: i32,
    s_off: i32,
    len: i32,
    score: i32,
}

struct BlastAaWordFinderSideChannels<'a> {
    init_hitlist: Option<&'a mut crate::extend::InitHitList>,
    ungapped_stats: Option<&'a mut crate::diagnostics::UngappedStats>,
}

/// NCBI: s_BlastAaExtendLeft (aa_ungapped.c).
/// Extend left from position (q_start-1, s_start-1) with x-dropoff.
/// Returns (best_score, left_displacement, num_identities).
/// Matches NCBI `s_BlastAaExtendLeft`.
#[inline]
fn s_blast_aa_extend_left(
    query: &[u8],
    subject: &[u8],
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
    q_start: usize,
    s_start: usize,
    x_dropoff: i32,
) -> (i32, i32, i32) {
    let mut score = 0i32;
    let mut best = 0i32;
    let mut best_d = 0i32;
    let mut ident = 0i32;
    let mut best_ident = 0i32;
    if q_start == 0 || s_start == 0 {
        return (0, 0, 0);
    }
    let mut qi = q_start - 1;
    let mut si = s_start - 1;
    let mut d = 1i32;
    loop {
        unsafe {
            let q = *query.as_ptr().add(qi);
            // Concatenated multi-query buffers place a sentinel between queries;
            // stop extension at the query boundary (NCBI stops at the NULLB
            // inter-context sentinel). Never triggers for single-query scans,
            // whose encoded query contains no sentinel byte.
            if q == PROTEIN_CONCAT_SENTINEL {
                break;
            }
            let s = *subject.as_ptr().add(si);
            score += *matrix.get_unchecked(q as usize).get_unchecked(s as usize);
            if q == s {
                ident += 1;
            }
        }
        if score > best {
            best = score;
            best_d = d;
            best_ident = ident;
        }
        // NCBI `s_BlastAaExtendLeft` (`aa_ungapped.c:915`) uses `>=` here
        // — the inline comment in the C source explicitly notes this is
        // intentional and differs from older code that used `>`.
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
    (best, best_d, best_ident)
}

/// NCBI: s_BlastAaExtendRight (aa_ungapped.c).
/// Extend right from position (q_start, s_start) with x-dropoff.
/// `init_score` is the cumulative score from left extension.
/// Returns (best_score, right_displacement, num_identities, s_last_off_delta).
///
/// `s_last_off_delta` is the rightmost subject offset EXAMINED by the loop
/// (one past the last cell processed) — used by NCBI's
/// `s_BlastAaWordFinder_TwoHit` to set `last_hit = s_last_off - (wordsize - 1)`
/// for diagonal bookkeeping (`aa_ungapped.c:599`). Distinct from `best_d`,
/// which is the displacement to the maximum-score cell. NCBI tracks this
/// via `s_BlastAaExtendRight`'s `s_last_off` out-parameter
/// (`aa_ungapped.c:864`: `*s_last_off = s_off + i;`).
#[inline]
fn s_blast_aa_extend_right(
    query: &[u8],
    subject: &[u8],
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
    q_start: usize,
    s_start: usize,
    x_dropoff: i32,
    init_score: i32,
) -> (i32, i32, i32, i32) {
    let mut score = init_score;
    let mut best = init_score;
    let mut best_d = 0i32;
    let mut ident = 0i32;
    let mut best_ident = 0i32;
    let mut qi = q_start;
    let mut si = s_start;
    let qlen = query.len();
    let slen = subject.len();
    let mut last_off_delta: i32 = 0;
    while qi < qlen && si < slen {
        unsafe {
            let q = *query.as_ptr().add(qi);
            // Concatenated multi-query buffers place a sentinel between queries;
            // stop extension at the query boundary (NCBI stops at the NULLB
            // inter-context sentinel). Never triggers for single-query scans,
            // whose encoded query contains no sentinel byte.
            if q == PROTEIN_CONCAT_SENTINEL {
                break;
            }
            let s = *subject.as_ptr().add(si);
            score += *matrix.get_unchecked(q as usize).get_unchecked(s as usize);
            if q == s {
                ident += 1;
            }
        }
        if score > best {
            best = score;
            best_d = (qi - q_start + 1) as i32;
            best_ident = ident;
        }
        // NCBI `s_BlastAaExtendRight` (`aa_ungapped.c:859`):
        //   `if (score <= 0 || (maxscore - score) >= dropoff) break;`
        // The right extension has BOTH a `score <= 0` early termination AND
        // the `>=` X-drop test (the inline comment notes the `>=` is intentional).
        // Left extension uses only the `>=` X-drop test.
        if score <= 0 || best - score >= x_dropoff {
            // NCBI's `*s_last_off = s_off + i` where `i` is the for-loop
            // counter value at break — i.e., the index of the cell that
            // triggered the break (NOT one past). With our loop using a
            // pre-increment pattern that DOESN'T advance qi/si on break,
            // `qi - q_start` matches NCBI's `i` exactly.
            last_off_delta = (qi - q_start) as i32;
            break;
        }
        qi += 1;
        si += 1;
        // After the increment, `qi - q_start` equals the number of
        // completed iterations. If we exit the while via the bounds
        // check, this matches NCBI's `i = n` after loop completion.
        last_off_delta = (qi - q_start) as i32;
    }
    (best, best_d, best_ident, last_off_delta)
}

/// NCBI: s_BlastAaExtendOneHit (aa_ungapped.c).
/// Port of NCBI `s_BlastAaExtendOneHit` (`aa_ungapped.c:1019`), reduced to
/// ordinary substitution-matrix scoring over the Rust slice representation.
fn s_blast_aa_extend_one_hit(
    query: &[u8],
    subject: &[u8],
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
    word_size: usize,
    x_dropoff: i32,
    q_off: usize,
    s_off: usize,
) -> (Option<ProteinHit>, i32) {
    let mut score = 0i32;
    let mut sum = 0i32;
    let mut q_left_off = q_off;
    let mut q_best_left_off = q_off;
    // NCBI s_BlastAaExtendOneHit (aa_ungapped.c:1030): init q_off + word_size
    // (persists only if no in-word position scores above 0).
    let mut q_right_off = q_off + word_size;

    for i in 0..word_size {
        let qi = q_off + i;
        let si = s_off + i;
        if qi >= query.len() || si >= subject.len() {
            break;
        }
        sum += matrix[query[qi] as usize][subject[si] as usize];
        if sum > score {
            score = sum;
            q_best_left_off = q_left_off;
            q_right_off = qi;
        } else if sum <= 0 {
            sum = 0;
            q_left_off = qi + 1;
        }
    }

    let q_left_off = q_best_left_off;
    let init_hit_width = q_right_off.saturating_sub(q_left_off) + 1;
    let diag_delta = s_off as isize - q_off as isize;
    let s_left_off = (q_left_off as isize + diag_delta) as usize;
    let s_right_off = (q_right_off as isize + diag_delta) as usize;
    let init_ident = (q_left_off..=q_right_off)
        .filter(|&qi| {
            let si = (qi as isize + diag_delta) as usize;
            query.get(qi) == subject.get(si)
        })
        .count() as i32;

    let (left_score, left_d, left_ident) = s_blast_aa_extend_left_with_score(
        query, subject, matrix, q_left_off, s_left_off, x_dropoff, score,
    );
    let (total_score, right_d, right_ident, s_last_off_delta) = s_blast_aa_extend_right(
        query,
        subject,
        matrix,
        q_right_off + 1,
        s_right_off + 1,
        x_dropoff,
        left_score,
    );
    let s_last_off = s_right_off as i32 + 1 + s_last_off_delta;
    if total_score <= 0 {
        return (None, s_last_off);
    }

    let qs = q_left_off - left_d as usize;
    let ss = s_left_off - left_d as usize;
    let alen = left_d + right_d + init_hit_width as i32;
    let ident = left_ident + init_ident + right_ident;
    (
        Some(ProteinHit {
            query_start: qs,
            query_end: qs + alen as usize,
            subject_start: ss,
            subject_end: ss + alen as usize,
            score: total_score,
            num_ident: ident,
            align_length: alen,
            mismatches: alen - ident,
            gap_opens: 0,
            qseq: None,
            sseq: None,
            scaled_score: None,
            gapped_start_q: 0,
            gapped_start_s: 0,
        }),
        s_last_off,
    )
}

#[inline]
fn s_blast_aa_extend_left_with_score(
    query: &[u8],
    subject: &[u8],
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
    q_start: usize,
    s_start: usize,
    x_dropoff: i32,
    init_score: i32,
) -> (i32, i32, i32) {
    let mut score = init_score;
    let mut best = init_score;
    let mut best_d = 0i32;
    let mut ident = 0i32;
    let mut best_ident = 0i32;
    if q_start == 0 || s_start == 0 {
        return (best, 0, 0);
    }
    let mut qi = q_start - 1;
    let mut si = s_start - 1;
    let mut d = 1i32;
    loop {
        unsafe {
            let q = *query.as_ptr().add(qi);
            // Concatenated multi-query buffers place a sentinel between queries;
            // stop extension at the query boundary (NCBI stops at the NULLB
            // inter-context sentinel). Never triggers for single-query scans,
            // whose encoded query contains no sentinel byte.
            if q == PROTEIN_CONCAT_SENTINEL {
                break;
            }
            let s = *subject.as_ptr().add(si);
            score += *matrix.get_unchecked(q as usize).get_unchecked(s as usize);
            if q == s {
                ident += 1;
            }
        }
        if score > best {
            best = score;
            best_d = d;
            best_ident = ident;
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
    (best, best_d, best_ident)
}

/// Result of `s_blast_aa_extend_two_hit`. NCBI distinguishes:
///   - left extension did NOT reach the first hit (`right_extend == FALSE`),
///   - left reached first AND right extended (`right_extend == TRUE`), with the
///     HSP only saved if `score >= cutoff`.
/// We mirror that here: `s_last_off` is `Some` exactly when right_extend was
/// true; `hit` is `Some` only when `score > 0` (caller-side cutoff).
pub enum TwoHitOutcome {
    /// Left extension did not reach first hit. NCBI's caller updates
    /// `diag.last_hit` from the current subject offset, not `s_last_off`.
    /// The left-only extension may still be saved before that diagonal update.
    NoReach { hit: Option<ProteinHit> },
    /// Left reached first hit (right extension ran). NCBI's caller updates
    /// `diag.last_hit` from `s_last_off - (wordsize - 1)`. HSP is only saved
    /// when present (`score > 0`).
    Reached {
        hit: Option<ProteinHit>,
        s_last_off: i32,
    },
}

/// NCBI: s_BlastAaExtendTwoHit (aa_ungapped.c).
/// Port of NCBI `s_BlastAaExtendTwoHit` (`aa_ungapped.c:1089`), reduced to
/// the local slice-based scanner state. Returns NCBI's `right_extend` /
/// `s_last_off` plus the optional ungapped HSP.
fn s_blast_aa_extend_two_hit(
    query: &[u8],
    subject: &[u8],
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
    word_size: usize,
    x_dropoff: i32,
    q_right_off: usize,
    s_right_off: usize,
    s_left_off: usize,
) -> TwoHitOutcome {
    let qlen = query.len();

    // Find best start within the word.
    let mut wscore = 0i32;
    let mut best_wscore = 0i32;
    let mut right_d = 0usize;
    for k in 0..word_size {
        let qi = q_right_off + k;
        let si = s_right_off + k;
        if qi < qlen && si < subject.len() {
            wscore += matrix[query[qi] as usize][subject[si] as usize];
        }
        if wscore > best_wscore {
            best_wscore = wscore;
            right_d = k + 1;
        }
    }
    let ext_q = q_right_off + right_d;
    let ext_s = s_right_off + right_d;

    let (left_score, left_d, left_ident) =
        s_blast_aa_extend_left(query, subject, matrix, ext_q, ext_s, x_dropoff);
    let reached_first = left_d >= (ext_s as i32 - s_left_off as i32);
    if !reached_first {
        // NCBI `s_BlastAaExtendTwoHit` leaves `right_extend = FALSE`, but it
        // still returns `left_score` plus hsp_q/hsp_s/hsp_len. The caller saves
        // that left-only HSP before updating diag.last_hit from subject_offset.
        let hit = (left_score > 0 && left_d > 0).then(|| {
            let qs = ext_q - left_d as usize;
            let qe = ext_q;
            let ss = ext_s - left_d as usize;
            let se = ext_s;
            let alen = left_d;
            ProteinHit {
                query_start: qs,
                query_end: qe,
                subject_start: ss,
                subject_end: se,
                score: left_score,
                num_ident: left_ident,
                align_length: alen,
                mismatches: alen - left_ident,
                gap_opens: 0,
                qseq: None,
                sseq: None,
                scaled_score: None,
                gapped_start_q: 0,
                gapped_start_s: 0,
            }
        });
        return TwoHitOutcome::NoReach { hit };
    }

    let (right_score, right_d_r, right_ident, s_last_off_delta) =
        s_blast_aa_extend_right(query, subject, matrix, ext_q, ext_s, x_dropoff, left_score);
    let total_score = left_score.max(right_score);
    let s_last_off = ext_s as i32 + s_last_off_delta;
    if total_score <= 0 {
        // NCBI still sets `*right_extend = TRUE` here — only HSP saving is
        // gated on `score >= cutoff`. Surface `s_last_off` so the caller can
        // update `diag.last_hit` via `s_last_off - (ws - 1)` exactly like NCBI.
        return TwoHitOutcome::Reached {
            hit: None,
            s_last_off,
        };
    }

    let qs = ext_q - left_d as usize;
    let qe = ext_q + right_d_r as usize;
    let ss = ext_s - left_d as usize;
    let se = ext_s + right_d_r as usize;
    let alen = (qe - qs) as i32;
    let ident = left_ident + right_ident;

    TwoHitOutcome::Reached {
        hit: Some(ProteinHit {
            query_start: qs,
            query_end: qe,
            subject_start: ss,
            subject_end: se,
            score: total_score,
            num_ident: ident,
            align_length: alen,
            mismatches: alen - ident,
            gap_opens: 0,
            qseq: None,
            sseq: None,
            scaled_score: None,
            gapped_start_q: 0,
            gapped_start_s: 0,
        }),
        s_last_off,
    }
}

/// NCBI: s_BlastAaWordFinder_TwoHit (aa_ungapped.c).
/// Port of NCBI `s_BlastAaWordFinder_TwoHit` (`aa_ungapped.c:440`), over the
/// Rust lookup-table representation.
///
/// Uses the NCBI BLAST+ two-hit algorithm: a word hit on a diagonal only
/// triggers ungapped extension if a second hit was seen on the same diagonal
/// within `TWO_HIT_WINDOW` positions. Use
/// [`blast_aa_word_finder_with_side_channels`] for the C-shaped
/// `BlastInitHitList`/stats effects.
pub fn s_blast_aa_word_finder_two_hit(
    query: &[u8],
    subject: &[u8],
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
    table: &ProteinLookupTable,
    x_dropoff: i32,
    window: i32,
    diag_buf: &mut Vec<(i32, bool)>,
) -> Vec<ProteinHit> {
    // Convenience wrapper: save every positively-scoring ungapped HSP (NCBI's
    // minimum `cutoff_score` is 1, so `>= 1` reproduces the historic `> 0`
    // gate). Callers that have a real per-context ungapped cutoff should go
    // through `protein_scan_with_table_reuse`/`blast_aa_word_finder`.
    s_blast_aa_word_finder_two_hit_scan(
        query, subject, matrix, table, x_dropoff, window, 1, diag_buf, None, None,
    )
    .hits
}

/// Compressed-alphabet two-hit word finder.
///
/// NCBI: `s_BlastAaWordFinder_TwoHit` (`aa_ungapped.c:440`) driven by the
/// compressed `scansub` callback (`s_BlastCompressedAaScanSubject`). The
/// per-hit diagonal/two-hit gating and extension are byte-identical to the
/// dense path above; only the hit-stream source differs (compressed scan
/// emits `(query_offset, subject_offset)` pairs ordered by subject offset).
/// The real (full) scoring matrix is used for extension, exactly like NCBI.
#[allow(clippy::too_many_arguments)]
fn s_blast_aa_word_finder_two_hit_compressed(
    query: &[u8],
    subject: &[u8],
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
    table: &ProteinLookupTable,
    x_dropoff: i32,
    window: i32,
    cutoff_score: i32,
    diag_buf: &mut Vec<(i32, bool)>,
    context_starts: Option<&[usize]>,
    mut side_channels: Option<&mut BlastAaWordFinderSideChannels<'_>>,
) -> BlastAaWordFinderScan {
    let word_size = table.word_size;
    let lookup = table.compressed.as_ref().expect("compressed table present");

    let mut diag_count = 1usize;
    while diag_count < query.len() + window as usize {
        diag_count <<= 1;
    }
    let diag_mask = diag_count - 1;
    diag_buf.clear();
    diag_buf.resize(diag_count, (-window, false));
    let diag_array = diag_buf;
    let diag_offset = window;

    let mut hits: Vec<ProteinHit> = Vec::new();
    let ws = word_size as i32;
    let mut total_hits = 0i32;
    let mut hits_extended = 0i32;

    // Collect compressed word hits in subject-offset-ordered batches, mirroring
    // NCBI's repeated `scansub` calls. Dense/repetitive subjects can exceed one
    // offset-pair buffer, so drive the existing scan_range cursor until the
    // subject is exhausted.
    // FIX C: scan into a per-thread reused buffer (NCBI reuses one
    // `offset_pairs` array across all subjects, sized once in
    // `blast_engine.c:1039-1040`); rayon gives each worker its own thread, so a
    // thread-local exactly mirrors that per-thread ownership.
    let array_size = subject.len().saturating_add(1).saturating_mul(8).max(1024);
    let mut scan_buf = COMPRESSED_SCAN_BUF.take();
    total_hits += s_blast_compressed_aa_scan_subject_batched(
        lookup,
        subject,
        array_size,
        &mut scan_buf,
        |scan_batch| {
            for pair in scan_batch {
                let q_pos = pair.query_offset as usize;
                let s_pos = pair.subject_offset as usize;
                let s_off = s_pos as i32;
                let diag = q_pos.wrapping_sub(s_pos) & diag_mask;
                let (last_hit, flag) = diag_array[diag];

                if flag {
                    if s_off + diag_offset < last_hit {
                        continue;
                    }
                    diag_array[diag] = (s_off + diag_offset, false);
                    continue;
                }
                let diff = s_off - (last_hit - diag_offset);
                if diff >= window {
                    diag_array[diag] = (s_off + diag_offset, false);
                    continue;
                }
                if diff < ws {
                    continue;
                }
                if let Some(starts) = context_starts {
                    let ctx_start = context_start_for_offset(starts, q_pos);
                    if (diff as usize) > q_pos || q_pos - (diff as usize) < ctx_start {
                        diag_array[diag] = (s_off + diag_offset, false);
                        continue;
                    }
                }

                let s_left_off = (last_hit - diag_offset + ws) as usize; // end of first hit
                let s_right_off = s_pos;
                let q_right_off = q_pos;
                hits_extended += 1;

                match s_blast_aa_extend_two_hit(
                    query,
                    subject,
                    matrix,
                    word_size,
                    x_dropoff,
                    q_right_off,
                    s_right_off,
                    s_left_off,
                ) {
                    TwoHitOutcome::Reached { hit, s_last_off } => {
                        // NCBI `s_BlastAaWordFinder_TwoHit` (`aa_ungapped.c:588`) saves
                        // the init HSP only when `score >= cutoffs->cutoff_score`; the
                        // diagonal `last_hit`/`flag` update below is INDEPENDENT of the
                        // cutoff and happens regardless of whether the HSP was saved.
                        diag_array[diag] = (s_last_off - (ws - 1) + diag_offset, true);
                        if let Some(hit) = hit {
                            if hit.score >= cutoff_score {
                                if let Some(channels) = side_channels.as_deref_mut() {
                                    if let Some(init_hitlist) = channels.init_hitlist.as_deref_mut()
                                    {
                                        crate::extend::blast_save_init_hsp(
                                            init_hitlist,
                                            hit.query_start as i32,
                                            hit.subject_start as i32,
                                            q_right_off as i32,
                                            s_right_off as i32,
                                            hit.align_length,
                                            hit.score,
                                        );
                                    }
                                }
                                hits.push(hit);
                            }
                        }
                    }
                    TwoHitOutcome::NoReach { hit } => {
                        if let Some(hit) = hit {
                            if hit.score >= cutoff_score {
                                if let Some(channels) = side_channels.as_deref_mut() {
                                    if let Some(init_hitlist) = channels.init_hitlist.as_deref_mut()
                                    {
                                        crate::extend::blast_save_init_hsp(
                                            init_hitlist,
                                            hit.query_start as i32,
                                            hit.subject_start as i32,
                                            q_right_off as i32,
                                            s_right_off as i32,
                                            hit.align_length,
                                            hit.score,
                                        );
                                    }
                                }
                                hits.push(hit);
                            }
                        }
                        diag_array[diag] = (s_off + diag_offset, false);
                    }
                }
            }
        },
    ) as i32;
    // FIX C: return the (now-populated, retained-capacity) buffer to the
    // thread-local for the next subject.
    COMPRESSED_SCAN_BUF.set(scan_buf);

    let scan = BlastAaWordFinderScan {
        hits,
        total_hits,
        hits_extended,
    };
    if let Some(channels) = side_channels {
        let saved_hits = channels
            .init_hitlist
            .as_deref()
            .map_or(scan.hits.len() as i32, |init_hitlist| {
                init_hitlist.total() as i32
            });
        crate::diagnostics::blast_ungapped_stats_update(
            channels.ungapped_stats.as_deref_mut(),
            scan.total_hits,
            scan.hits_extended,
            saved_hits,
        );
    }
    scan
}

#[allow(clippy::too_many_arguments)]
fn s_blast_aa_word_finder_two_hit_scan(
    query: &[u8],
    subject: &[u8],
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
    table: &ProteinLookupTable,
    x_dropoff: i32,
    window: i32,
    cutoff_score: i32,
    diag_buf: &mut Vec<(i32, bool)>,
    context_starts: Option<&[usize]>,
    mut side_channels: Option<&mut BlastAaWordFinderSideChannels<'_>>,
) -> BlastAaWordFinderScan {
    let word_size = table.word_size;
    if query.len() < word_size || subject.len() < word_size {
        return BlastAaWordFinderScan::default();
    }

    if table.compressed.is_some() {
        return s_blast_aa_word_finder_two_hit_compressed(
            query,
            subject,
            matrix,
            table,
            x_dropoff,
            window,
            cutoff_score,
            diag_buf,
            context_starts,
            side_channels,
        );
    }

    let mut diag_count = 1usize;
    while diag_count < query.len() + window as usize {
        diag_count <<= 1;
    }
    let diag_mask = diag_count - 1;
    diag_buf.clear();
    diag_buf.resize(diag_count, (-window, false));
    let diag_array = diag_buf;

    let mut hits: Vec<ProteinHit> = Vec::new();
    let ws = word_size as i32;
    let mask = (1usize << (word_size * CHARSIZE)) - 1;

    let last_pos = subject.len() - word_size;
    let mut hash = word_hash(&subject[0..word_size], 0);

    // SAFETY: All indices are bounds-checked by loop structure:
    // - s_pos in 0..=last_pos where last_pos = subject.len() - word_size
    // - hash < table_size (= alphabet_size^word_size), PV/backbone sized to table_size
    // - diag < diag_count = query.len() + subject.len()
    // Removing bounds checks in this hot loop matches NCBI C approach.
    let pv = table.pv.as_ptr();
    let backbone = table.backbone.as_ptr();
    let overflow = table.overflow.as_ptr();
    let subj = subject.as_ptr();
    let diag_ptr = diag_array.as_mut_ptr();
    let diag_offset = window;
    let mut total_hits = 0i32;
    let mut hits_extended = 0i32;

    for s_pos in 0..=last_pos {
        if s_pos > 0 {
            unsafe {
                // NCBI ComputeTableIndexIncremental: shift + or + mask (no multiply)
                let new = *subj.add(s_pos + word_size - 1) as usize;
                hash = ((hash << CHARSIZE) | new) & mask;
            }
        }

        unsafe {
            if *pv.add(hash >> 6) & (1u64 << (hash & 63)) == 0 {
                continue;
            }

            let cell = &*backbone.add(hash);
            let num = cell.num_used as usize;
            if num == 0 {
                continue;
            }
            total_hits += num as i32;

            let (hit_ptr, hit_len) = if num <= HITS_PER_CELL {
                (cell.entries.as_ptr(), num)
            } else {
                let cursor = cell.entries[0] as usize;
                (overflow.add(cursor), num)
            };

            for i in 0..hit_len {
                let q_pos = *hit_ptr.add(i) as usize;
                let s_off = s_pos as i32;
                let diag = q_pos.wrapping_sub(s_pos) & diag_mask;
                let (last_hit, flag) = *diag_ptr.add(diag);

                if flag {
                    if s_off + diag_offset < last_hit {
                        continue;
                    }
                    *diag_ptr.add(diag) = (s_off + diag_offset, false);
                    continue;
                }
                let diff = s_off - (last_hit - diag_offset);
                if diff >= window {
                    *diag_ptr.add(diag) = (s_off + diag_offset, false);
                    continue;
                }
                if diff < ws {
                    continue;
                }
                if let Some(starts) = context_starts {
                    let ctx_start = context_start_for_offset(starts, q_pos);
                    if (diff as usize) > q_pos || q_pos - (diff as usize) < ctx_start {
                        *diag_ptr.add(diag) = (s_off + diag_offset, false);
                        continue;
                    }
                }

                let s_left_off = (last_hit - diag_offset + ws) as usize; // end of first hit
                let s_right_off = s_pos;
                let q_right_off = q_pos;
                hits_extended += 1;

                match s_blast_aa_extend_two_hit(
                    query,
                    subject,
                    matrix,
                    word_size,
                    x_dropoff,
                    q_right_off,
                    s_right_off,
                    s_left_off,
                ) {
                    TwoHitOutcome::Reached { hit, s_last_off } => {
                        // NCBI `s_BlastAaWordFinder_TwoHit` (`aa_ungapped.c:412`)
                        // sets diag.last_hit from `s_last_off` whenever
                        // `right_extend == TRUE`, even if the HSP score was
                        // below the save cutoff. The diagonal update is
                        // INDEPENDENT of the cutoff; only HSP saving is gated on
                        // `score >= cutoffs->cutoff_score` (`aa_ungapped.c:588`).
                        *diag_ptr.add(diag) = (s_last_off - (ws - 1) + diag_offset, true);
                        if let Some(hit) = hit {
                            if hit.score >= cutoff_score {
                                let record = BlastAaInitHitRecord {
                                    q_start: hit.query_start as i32,
                                    s_start: hit.subject_start as i32,
                                    q_off: q_right_off as i32,
                                    s_off: s_right_off as i32,
                                    len: hit.align_length,
                                    score: hit.score,
                                };
                                if let Some(channels) = side_channels.as_deref_mut() {
                                    if let Some(init_hitlist) = channels.init_hitlist.as_deref_mut()
                                    {
                                        crate::extend::blast_save_init_hsp(
                                            init_hitlist,
                                            record.q_start,
                                            record.s_start,
                                            record.q_off,
                                            record.s_off,
                                            record.len,
                                            record.score,
                                        );
                                    }
                                }
                                hits.push(hit);
                            }
                        }
                    }
                    TwoHitOutcome::NoReach { hit } => {
                        if let Some(hit) = hit {
                            if hit.score >= cutoff_score {
                                let record = BlastAaInitHitRecord {
                                    q_start: hit.query_start as i32,
                                    s_start: hit.subject_start as i32,
                                    q_off: q_right_off as i32,
                                    s_off: s_right_off as i32,
                                    len: hit.align_length,
                                    score: hit.score,
                                };
                                if let Some(channels) = side_channels.as_deref_mut() {
                                    if let Some(init_hitlist) = channels.init_hitlist.as_deref_mut()
                                    {
                                        crate::extend::blast_save_init_hsp(
                                            init_hitlist,
                                            record.q_start,
                                            record.s_start,
                                            record.q_off,
                                            record.s_off,
                                            record.len,
                                            record.score,
                                        );
                                    }
                                }
                                hits.push(hit);
                            }
                        }
                        *diag_ptr.add(diag) = (s_off + diag_offset, false);
                    }
                }
            }
        }
    }

    let scan = BlastAaWordFinderScan {
        hits,
        total_hits,
        hits_extended,
    };
    if let Some(channels) = side_channels {
        let saved_hits = channels
            .init_hitlist
            .as_deref()
            .map_or(scan.hits.len() as i32, |init_hitlist| {
                init_hitlist.total() as i32
            });
        crate::diagnostics::blast_ungapped_stats_update(
            channels.ungapped_stats.as_deref_mut(),
            scan.total_hits,
            scan.hits_extended,
            saved_hits,
        );
    }
    scan
}

#[inline]
fn context_start_for_offset(context_starts: &[usize], query_offset: usize) -> usize {
    match context_starts.binary_search(&query_offset) {
        Ok(i) => context_starts[i],
        Err(0) => 0,
        Err(i) => context_starts[i - 1],
    }
}

/// NCBI: s_BlastAaWordFinder_OneHit (aa_ungapped.c).
/// Port of NCBI `s_BlastAaWordFinder_OneHit` (`aa_ungapped.c:713`) over the
/// local lookup-table and diagonal-scratch representations.
pub fn s_blast_aa_word_finder_one_hit(
    query: &[u8],
    subject: &[u8],
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
    table: &ProteinLookupTable,
    x_dropoff: i32,
    diag_buf: &mut Vec<(i32, bool)>,
) -> Vec<ProteinHit> {
    s_blast_aa_word_finder_one_hit_scan(query, subject, matrix, table, x_dropoff, 1, diag_buf, None)
        .hits
}

/// Compressed-alphabet one-hit word finder. Mirrors the dense one-hit gating
/// (`s_BlastAaWordFinder_OneHit`) but sources hits from the compressed scan.
fn s_blast_aa_word_finder_one_hit_compressed(
    query: &[u8],
    subject: &[u8],
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
    table: &ProteinLookupTable,
    x_dropoff: i32,
    cutoff_score: i32,
    diag_buf: &mut Vec<(i32, bool)>,
    mut side_channels: Option<&mut BlastAaWordFinderSideChannels<'_>>,
) -> BlastAaWordFinderScan {
    let word_size = table.word_size;
    let lookup = table.compressed.as_ref().expect("compressed table present");

    // NCBI s_BlastDiagTableNew (blast_extend.c:52) sizes the one-hit diagonal
    // table to the query length only (window_size 0); the modulo aliasing is
    // intentional. (Two-hit paths size to query.len()+window.)
    let mut diag_count = 1usize;
    while diag_count < query.len() {
        diag_count <<= 1;
    }
    let diag_mask = diag_count - 1;
    diag_buf.clear();
    diag_buf.resize(diag_count, (0, false));
    let diag_array = diag_buf;

    let mut hits = Vec::new();
    let ws = word_size as i32;
    let mut total_hits = 0i32;
    let mut hits_extended = 0i32;

    // FIX C: reuse the per-thread scan buffer (see two-hit path above).
    let array_size = subject.len().saturating_add(1).saturating_mul(8).max(1024);
    let mut scan_buf = COMPRESSED_SCAN_BUF.take();
    total_hits += s_blast_compressed_aa_scan_subject_batched(
        lookup,
        subject,
        array_size,
        &mut scan_buf,
        |scan_batch| {
            for pair in scan_batch {
                let q_pos = pair.query_offset as usize;
                let s_pos = pair.subject_offset as usize;
                let diag = s_pos.wrapping_sub(q_pos) & diag_mask;
                let last_hit = diag_array[diag].0;
                let diff = s_pos as i32 - last_hit;
                if diff < 0 {
                    continue;
                }
                hits_extended += 1;

                let (hit, s_last_off) = s_blast_aa_extend_one_hit(
                    query, subject, matrix, word_size, x_dropoff, q_pos, s_pos,
                );
                diag_array[diag] = (s_last_off - (ws - 1), false);
                if let Some(hit) = hit {
                    if hit.score >= cutoff_score {
                        if let Some(channels) = side_channels.as_deref_mut() {
                            if let Some(init_hitlist) = channels.init_hitlist.as_deref_mut() {
                                crate::extend::blast_save_init_hsp(
                                    init_hitlist,
                                    hit.query_start as i32,
                                    hit.subject_start as i32,
                                    q_pos as i32,
                                    s_pos as i32,
                                    hit.align_length,
                                    hit.score,
                                );
                            }
                        }
                        hits.push(hit);
                    }
                }
            }
        },
    ) as i32;
    // FIX C: return the buffer to the thread-local for the next subject.
    COMPRESSED_SCAN_BUF.set(scan_buf);

    let scan = BlastAaWordFinderScan {
        hits,
        total_hits,
        hits_extended,
    };
    if let Some(channels) = side_channels {
        let saved_hits = channels
            .init_hitlist
            .as_deref()
            .map_or(scan.hits.len() as i32, |init_hitlist| {
                init_hitlist.total() as i32
            });
        crate::diagnostics::blast_ungapped_stats_update(
            channels.ungapped_stats.as_deref_mut(),
            scan.total_hits,
            scan.hits_extended,
            saved_hits,
        );
    }
    scan
}

fn s_blast_aa_word_finder_one_hit_scan(
    query: &[u8],
    subject: &[u8],
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
    table: &ProteinLookupTable,
    x_dropoff: i32,
    cutoff_score: i32,
    diag_buf: &mut Vec<(i32, bool)>,
    mut side_channels: Option<&mut BlastAaWordFinderSideChannels<'_>>,
) -> BlastAaWordFinderScan {
    let word_size = table.word_size;
    if query.len() < word_size || subject.len() < word_size {
        return BlastAaWordFinderScan::default();
    }

    if table.compressed.is_some() {
        return s_blast_aa_word_finder_one_hit_compressed(
            query,
            subject,
            matrix,
            table,
            x_dropoff,
            cutoff_score,
            diag_buf,
            side_channels,
        );
    }

    // NCBI s_BlastDiagTableNew (blast_extend.c:52) sizes the one-hit diagonal
    // table to the query length only (window_size 0); the modulo aliasing is
    // intentional. (Two-hit paths size to query.len()+window.)
    let mut diag_count = 1usize;
    while diag_count < query.len() {
        diag_count <<= 1;
    }
    let diag_mask = diag_count - 1;
    diag_buf.clear();
    diag_buf.resize(diag_count, (0, false));

    let mut hits = Vec::new();
    let ws = word_size as i32;
    let mask = (1usize << (word_size * CHARSIZE)) - 1;
    let last_pos = subject.len() - word_size;
    let mut hash = word_hash(&subject[0..word_size], 0);

    let pv = table.pv.as_ptr();
    let backbone = table.backbone.as_ptr();
    let overflow = table.overflow.as_ptr();
    let subj = subject.as_ptr();
    let diag_ptr = diag_buf.as_mut_ptr();
    let mut total_hits = 0i32;
    let mut hits_extended = 0i32;

    for s_pos in 0..=last_pos {
        if s_pos > 0 {
            unsafe {
                let new = *subj.add(s_pos + word_size - 1) as usize;
                hash = ((hash << CHARSIZE) | new) & mask;
            }
        }

        unsafe {
            if *pv.add(hash >> 6) & (1u64 << (hash & 63)) == 0 {
                continue;
            }

            let cell = &*backbone.add(hash);
            let num = cell.num_used as usize;
            if num == 0 {
                continue;
            }
            total_hits += num as i32;

            let (hit_ptr, hit_len) = if num <= HITS_PER_CELL {
                (cell.entries.as_ptr(), num)
            } else {
                let cursor = cell.entries[0] as usize;
                (overflow.add(cursor), num)
            };

            for i in 0..hit_len {
                let q_pos = *hit_ptr.add(i) as usize;
                let diag = s_pos.wrapping_sub(q_pos) & diag_mask;
                let last_hit = (*diag_ptr.add(diag)).0;
                let diff = s_pos as i32 - last_hit;
                if diff < 0 {
                    continue;
                }
                hits_extended += 1;

                let (hit, s_last_off) = s_blast_aa_extend_one_hit(
                    query, subject, matrix, word_size, x_dropoff, q_pos, s_pos,
                );
                *diag_ptr.add(diag) = (s_last_off - (ws - 1), false);
                if let Some(hit) = hit {
                    if hit.score >= cutoff_score {
                        let record = BlastAaInitHitRecord {
                            q_start: hit.query_start as i32,
                            s_start: hit.subject_start as i32,
                            q_off: q_pos as i32,
                            s_off: s_pos as i32,
                            len: hit.align_length,
                            score: hit.score,
                        };
                        if let Some(channels) = side_channels.as_deref_mut() {
                            if let Some(init_hitlist) = channels.init_hitlist.as_deref_mut() {
                                crate::extend::blast_save_init_hsp(
                                    init_hitlist,
                                    record.q_start,
                                    record.s_start,
                                    record.q_off,
                                    record.s_off,
                                    record.len,
                                    record.score,
                                );
                            }
                        }
                        hits.push(hit);
                    }
                }
            }
        }
    }

    let scan = BlastAaWordFinderScan {
        hits,
        total_hits,
        hits_extended,
    };
    if let Some(channels) = side_channels {
        let saved_hits = channels
            .init_hitlist
            .as_deref()
            .map_or(scan.hits.len() as i32, |init_hitlist| {
                init_hitlist.total() as i32
            });
        crate::diagnostics::blast_ungapped_stats_update(
            channels.ungapped_stats.as_deref_mut(),
            scan.total_hits,
            scan.hits_extended,
            saved_hits,
        );
    }
    scan
}

/// NCBI: BlastAaWordFinder (aa_ungapped.c).
/// C-shaped dispatch wrapper over the locally represented protein word-finder
/// branch. NCBI switches between one-hit and two-hit scanning from the lookup
/// options window.
#[allow(clippy::too_many_arguments)]
pub fn blast_aa_word_finder(
    query: &[u8],
    subject: &[u8],
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
    table: &ProteinLookupTable,
    x_dropoff: i32,
    window_size: i32,
    cutoff_score: i32,
    diag_buf: &mut Vec<(i32, bool)>,
) -> Result<Vec<ProteinHit>, BlastAaWordFinderUnsupported> {
    let mut scan = if window_size <= 0 {
        s_blast_aa_word_finder_one_hit_scan(
            query,
            subject,
            matrix,
            table,
            x_dropoff,
            cutoff_score,
            diag_buf,
            None,
        )
    } else {
        s_blast_aa_word_finder_two_hit_scan(
            query,
            subject,
            matrix,
            table,
            x_dropoff,
            window_size,
            cutoff_score,
            diag_buf,
            None,
            None,
        )
    };
    blast_init_protein_hits_sort_by_score(&mut scan.hits);
    Ok(scan.hits)
}

/// blast-rs: side-channel preserving wrapper around `BlastAaWordFinder`
/// (`aa_ungapped.c`).
///
/// Side-channel preserving wrapper for callers that want the C-shaped
/// `BlastInitHitList` and `BlastUngappedStats` effects. The scan records the
/// exact seed offsets that C passes to `BlastSaveInitHsp`, then sorts the
/// resulting init-hit list by score just like `BlastAaWordFinder`.
#[allow(clippy::too_many_arguments)]
pub fn blast_aa_word_finder_with_side_channels(
    query: &[u8],
    subject: &[u8],
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
    table: &ProteinLookupTable,
    x_dropoff: i32,
    window_size: i32,
    cutoff_score: i32,
    diag_buf: &mut Vec<(i32, bool)>,
    init_hitlist: Option<&mut crate::extend::InitHitList>,
    ungapped_stats: Option<&mut crate::diagnostics::UngappedStats>,
) -> Result<BlastAaWordFinderScan, BlastAaWordFinderUnsupported> {
    let mut side_channels = BlastAaWordFinderSideChannels {
        init_hitlist,
        ungapped_stats,
    };
    let mut scan = if window_size <= 0 {
        s_blast_aa_word_finder_one_hit_scan(
            query,
            subject,
            matrix,
            table,
            x_dropoff,
            cutoff_score,
            diag_buf,
            Some(&mut side_channels),
        )
    } else {
        s_blast_aa_word_finder_two_hit_scan(
            query,
            subject,
            matrix,
            table,
            x_dropoff,
            window_size,
            cutoff_score,
            diag_buf,
            None,
            Some(&mut side_channels),
        )
    };
    blast_init_protein_hits_sort_by_score(&mut scan.hits);

    if let Some(init_hitlist) = side_channels.init_hitlist {
        crate::extend::blast_init_hit_list_sort_by_score(init_hitlist);
    }

    Ok(scan)
}

/// blast-rs: Rust convenience wrapper around `s_blast_aa_word_finder_two_hit`
/// that lets higher-level search code reuse the diagonal scratch buffer across
/// subjects.
#[allow(clippy::too_many_arguments)]
pub fn protein_scan_with_table_reuse(
    query: &[u8],
    subject: &[u8],
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
    table: &ProteinLookupTable,
    x_dropoff: i32,
    window: i32,
    cutoff_score: i32,
    diag_buf: &mut Vec<(i32, bool)>,
) -> Vec<ProteinHit> {
    // Keep the C-shaped dispatch edge visible: NCBI routes callers through
    // `BlastAaWordFinder`, which then selects the one-hit or two-hit scanner.
    // `cutoff_score` is NCBI's per-context ungapped `cutoffs->cutoff_score`,
    // gating which ungapped HSPs are saved (`aa_ungapped.c:588`).
    blast_aa_word_finder(
        query,
        subject,
        matrix,
        table,
        x_dropoff,
        window,
        cutoff_score,
        diag_buf,
    )
    .expect("local protein lookup wrapper supports one-hit and two-hit scans")
}

/// Scan a concatenated query buffer while preserving NCBI's same-context
/// two-hit diagonal rule. `context_starts` must be sorted query offsets.
#[allow(clippy::too_many_arguments)]
pub fn protein_scan_with_table_reuse_contexts(
    query: &[u8],
    subject: &[u8],
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
    table: &ProteinLookupTable,
    x_dropoff: i32,
    window: i32,
    cutoff_score: i32,
    diag_buf: &mut Vec<(i32, bool)>,
    context_starts: &[usize],
) -> Vec<ProteinHit> {
    if window <= 0 {
        let mut scan = s_blast_aa_word_finder_one_hit_scan(
            query,
            subject,
            matrix,
            table,
            x_dropoff,
            cutoff_score,
            diag_buf,
            None,
        )
        .hits;
        blast_init_protein_hits_sort_by_score(&mut scan);
        scan
    } else {
        let mut scan = s_blast_aa_word_finder_two_hit_scan(
            query,
            subject,
            matrix,
            table,
            x_dropoff,
            window,
            cutoff_score,
            diag_buf,
            Some(context_starts),
            None,
        )
        .hits;
        // This path scans concatenated Rust query buffers and must pass the
        // same-context guard directly into the two-hit body. NCBI still sorts at
        // `BlastAaWordFinder` after scanning, so do the equivalent here.
        blast_init_protein_hits_sort_by_score(&mut scan);
        scan
    }
}

/// blast-rs: Rust convenience wrapper that builds a protein lookup table, scans
/// ungapped seeds, then runs local gapped extension on accepted seeds.
///
/// This is the standard two-phase BLAST approach:
/// 1. Find seeds via lookup table + ungapped extension
/// 2. For seeds above a cutoff, perform gapped alignment with X-dropoff DP
pub fn protein_gapped_scan(
    query: &[u8],
    subject: &[u8],
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
    word_size: usize,
    threshold: f64,
    ungap_x_dropoff: i32,
    gap_open: i32,
    gap_extend: i32,
    gap_x_dropoff: i32,
    ungap_cutoff: i32,
) -> Vec<ProteinHit> {
    let table = ProteinLookupTable::build(query, word_size, matrix, threshold);
    protein_gapped_scan_with_table(
        query,
        subject,
        matrix,
        &table,
        ungap_x_dropoff,
        gap_open,
        gap_extend,
        gap_x_dropoff,
        ungap_cutoff,
    )
}

/// blast-rs: Rust convenience wrapper for the two-phase protein scan when the
/// caller already owns a pre-built lookup table.
pub fn protein_gapped_scan_with_table(
    query: &[u8],
    subject: &[u8],
    matrix: &[[i32; AA_SIZE]; AA_SIZE],
    table: &ProteinLookupTable,
    ungap_x_dropoff: i32,
    gap_open: i32,
    gap_extend: i32,
    gap_x_dropoff: i32,
    ungap_cutoff: i32,
) -> Vec<ProteinHit> {
    // Phase 1: ungapped seeds
    let ungapped = protein_scan_with_table(query, subject, matrix, table, ungap_x_dropoff);

    // Phase 2: gapped extension on seeds above cutoff
    let mut gapped_hits = Vec::new();
    // Track which diagonals we've already done gapped extension on
    let diag_count = query.len() + subject.len();
    let mut gapped_diag: Vec<bool> = vec![false; diag_count];

    for uh in &ungapped {
        if uh.score < ungap_cutoff {
            continue;
        }

        let (seed_q, seed_s) = blast_get_start_for_gapped_alignment(
            query,
            subject,
            uh.query_start,
            uh.query_end.saturating_sub(uh.query_start),
            uh.subject_start,
            uh.subject_end.saturating_sub(uh.subject_start),
            matrix,
        );
        let diag = seed_s + query.len() - seed_q;
        if diag < diag_count && gapped_diag[diag] {
            continue;
        }
        if diag < diag_count {
            gapped_diag[diag] = true;
        }

        if let Some(gr) = protein_gapped_align(
            query,
            subject,
            seed_q,
            seed_s,
            matrix,
            gap_open,
            gap_extend,
            gap_x_dropoff,
        ) {
            let q_slice = &query[gr.query_start..gr.query_end];
            let s_slice = &subject[gr.subject_start..gr.subject_end];
            let (qseq, sseq) =
                gr.edit_script
                    .render_alignment(q_slice, s_slice, ncbistdaa_to_aminoacid_char);
            gapped_hits.push(ProteinHit {
                query_start: gr.query_start,
                query_end: gr.query_end,
                subject_start: gr.subject_start,
                subject_end: gr.subject_end,
                score: gr.score,
                num_ident: gr.num_ident,
                align_length: gr.align_length,
                mismatches: gr.mismatches,
                gap_opens: gr.gap_opens,
                qseq: Some(qseq),
                sseq: Some(sseq),
                scaled_score: None,
                gapped_start_q: seed_q,
                gapped_start_s: seed_s,
            });
        }
    }

    // If no seeds passed the cutoff, fall back to ungapped hits
    if gapped_hits.is_empty() {
        return ungapped;
    }

    gapped_hits.sort_by(score_compare_protein_hits);
    gapped_hits
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrix::AA_SIZE;

    fn simple_matrix() -> [[i32; AA_SIZE]; AA_SIZE] {
        let mut m = [[0i32; AA_SIZE]; AA_SIZE];
        for i in 1..21 {
            m[i][i] = 4;
        }
        for i in 1..21 {
            for j in 1..21 {
                if i != j {
                    m[i][j] = -1;
                }
            }
        }
        m
    }

    #[test]
    fn compressed_reciprocal_alphabet_size_matches_ncbi_ceil_constants() {
        assert_eq!(compressed_reciprocal_alphabet_size(15), 286_331_154);
        assert_eq!(compressed_reciprocal_alphabet_size(10), 429_496_730);
        assert_eq!(compressed_reciprocal_alphabet_size(0), 0);
    }

    #[test]
    fn test_word_hash() {
        // Shift-based hash: h = (h << 5) | b for each byte
        // [1,2,3] = (1 << 10) | (2 << 5) | 3 = 1024 + 64 + 3 = 1091
        assert_eq!(word_hash(&[1, 2, 3], 28), (1 << 10) | (2 << 5) | 3);
        assert_eq!(word_hash(&[0, 0, 0], 28), 0);
        assert_eq!(word_hash(&[0, 0, 1], 28), 1);
    }

    #[test]
    fn aa_lookup_backbone_size_matches_ncbi_valid_index_range() {
        let c_sized = aa_lookup_backbone_size(3);
        let full_charsize_space = 1usize << (3 * CHARSIZE);

        assert_eq!(c_sized, word_hash(&[27, 27, 27], AA_SIZE) + 1);
        assert_eq!(c_sized, 28_540);
        assert!(c_sized < full_charsize_space);
    }

    #[test]
    fn test_lookup_table_build() {
        let m = simple_matrix();
        // Query: 3 amino acids, word_size=3, threshold=12 (exact match only with score 4*3=12)
        let query = vec![1u8, 2, 3];
        let table = ProteinLookupTable::build(&query, 3, &m, 12.0);
        assert_eq!(table.backbone.len(), aa_lookup_backbone_size(3));
        let hash = word_hash(&[1, 2, 3], AA_SIZE);
        let cell = &table.backbone[hash];
        assert!(cell.num_used > 0, "Cell should have hits");
        let hits: &[i32] = if (cell.num_used as usize) <= HITS_PER_CELL {
            &cell.entries[..cell.num_used as usize]
        } else {
            let c = cell.entries[0] as usize;
            &table.overflow[c..c + cell.num_used as usize]
        };
        assert!(
            hits.contains(&0),
            "Exact match word should have query offset 0"
        );
        // PV bit should be set.
        assert_ne!(table.pv[hash >> 6] & (1u64 << (hash & 63)), 0);
    }

    #[test]
    fn blast_aa_lookup_finalize_rebuilds_pv_and_compacts_overflow() {
        let mut table = ProteinLookupTable {
            word_size: 3,
            backbone: vec![
                BackboneCell {
                    num_used: 0,
                    entries: [99; HITS_PER_CELL],
                },
                BackboneCell {
                    num_used: 2,
                    entries: [7, 8, 99],
                },
                BackboneCell {
                    num_used: 4,
                    entries: [2, 0, 0],
                },
            ],
            overflow: vec![100, 101, 9, 10, 11, 12, 200],
            pv: vec![u64::MAX],
            compressed: None,
        };

        assert_eq!(
            blast_aa_lookup_finalize(Some(&mut table), AaLookupBoneType::Backbone),
            0
        );
        assert_eq!(table.pv[0] & 1, 0);
        assert_ne!(table.pv[0] & (1u64 << 1), 0);
        assert_ne!(table.pv[0] & (1u64 << 2), 0);
        assert_eq!(table.backbone[0].entries, [0; HITS_PER_CELL]);
        assert_eq!(table.backbone[1].entries, [7, 8, 0]);
        assert_eq!(table.backbone[2].entries[0], 0);
        assert_eq!(table.overflow, vec![9, 10, 11, 12]);
        assert_eq!(
            blast_aa_lookup_finalize(None, AaLookupBoneType::Smallbone),
            -1
        );
    }

    #[test]
    fn add_pssm_neighboring_words_indexes_threshold_words() {
        let mut scores = vec![[-10i32; AA_SIZE]; 3];
        scores[0][1] = 6;
        scores[1][2] = 5;
        scores[1][3] = 7;
        scores[2][4] = 6;
        let pssm = Pssm {
            scores,
            length: 3,
            info_content: vec![0.0; 3],
            start_numerator: None,
            ancillary_gap_kbp: None,
        };
        let table_size = 1usize << (2 * CHARSIZE);
        let mut backbone = vec![Vec::new(); table_size];

        assert_eq!(
            s_add_pssm_neighboring_words(&mut backbone, &pssm, 2, 11, 5, None),
            4
        );

        assert_eq!(backbone[word_hash(&[1, 2], AA_SIZE)], vec![5]);
        assert_eq!(backbone[word_hash(&[1, 3], AA_SIZE)], vec![5]);
        assert_eq!(backbone[word_hash(&[2, 4], AA_SIZE)], vec![6]);
        assert_eq!(backbone[word_hash(&[3, 4], AA_SIZE)], vec![6]);
    }

    #[test]
    fn compressed_lookup_add_word_hit_transitions_to_overflow_cells() {
        let mut lookup = BlastCompressedAaLookupTable::new(3, 15, 1);
        for offset in 10..18 {
            assert_eq!(s_compressed_lookup_add_word_hit(&mut lookup, 0, offset), 0);
        }

        let cell = lookup.backbone[0];
        assert_eq!(cell.num_used, 8);
        assert_eq!(cell.query_offset, 10);
        assert_eq!(cell.overflow_offsets(), [11, 12]);
        let head = cell.overflow_head().expect("overflow head");
        assert_eq!(lookup.overflow_cells[head].query_offsets, [17, 0, 0, 0]);
        let next = lookup.overflow_cells[head]
            .next
            .expect("older overflow cell");
        assert_eq!(lookup.overflow_cells[next].query_offsets, [13, 14, 15, 16]);
        assert_eq!(lookup.curr_overflow_bank, 0);
        assert_eq!(lookup.curr_overflow_cell, 2);
    }

    #[test]
    fn compressed_lookup_add_encoded_uses_c_radices() {
        let index5 = 1 + 15 * 2 + 225 * 3 + 3_375 * 4 + 50_625 * 5;
        let mut lookup5 = BlastCompressedAaLookupTable::new(5, 15, index5 + 1);
        assert_eq!(
            s_compressed_lookup_add_encoded(&mut lookup5, &[1, 2, 3, 4, 5], 42),
            0
        );
        assert_eq!(lookup5.backbone[index5].query_offset, 42);

        let index6 = index5 + 759_375 * 6;
        let mut lookup6 = BlastCompressedAaLookupTable::new(6, 15, index6 + 1);
        assert_eq!(
            s_compressed_lookup_add_encoded(&mut lookup6, &[1, 2, 3, 4, 5, 6], 43),
            0
        );
        assert_eq!(lookup6.backbone[index6].query_offset, 43);

        let index7 = 1 + 10 * 2 + 100 * 3 + 1_000 * 4 + 10_000 * 5 + 100_000 * 6;
        let index7 = index7 + 1_000_000 * 7;
        let mut lookup7 = BlastCompressedAaLookupTable::new(7, 10, index7 + 1);
        assert_eq!(
            s_compressed_lookup_add_encoded(&mut lookup7, &[1, 2, 3, 4, 5, 6, 7], 44),
            0
        );
        assert_eq!(lookup7.backbone[index7].query_offset, 44);
    }

    #[test]
    fn compressed_lookup_add_encoded_rejects_unsupported_word_lengths() {
        let mut lookup = BlastCompressedAaLookupTable::new(4, 15, 1);

        assert_eq!(
            s_compressed_lookup_add_encoded(&mut lookup, &[1, 2, 3, 4], 42),
            -1
        );
        assert_eq!(lookup.backbone[0].num_used, 0);
    }

    #[test]
    fn compressed_add_neighboring_words_uses_sorted_matrix_pruning() {
        let index = 1 + 15 + 225 + 3_375 + 50_625;
        let mut lookup = BlastCompressedAaLookupTable::new(5, 2, index + 1);
        lookup.threshold = 25;
        lookup.compress_table = vec![u8::MAX; AA_SIZE];
        lookup.compress_table[1] = 1;

        let mut matrix = [[-10i32; AA_SIZE]; AA_SIZE];
        matrix[1][1] = 5;
        matrix[1][0] = -10;
        assert_eq!(
            s_compressed_add_neighboring_words(&mut lookup, &matrix, &[1, 1, 1, 1, 1], None),
            1
        );

        assert_eq!(lookup.backbone[index].query_offset, 0);
        assert_eq!(lookup.backbone[index].num_used, 1);
        assert_eq!(lookup.exact_matches, 1);
        assert_eq!(lookup.neighbor_matches, 0);
    }

    #[test]
    fn compressed_lookup_finalize_builds_pv_and_longest_chain() {
        let mut lookup = BlastCompressedAaLookupTable::new(5, 15, 130);
        lookup.backbone[3].num_used = 2;
        lookup.backbone[64].num_used = 5;
        lookup.backbone[129].num_used = 1;

        assert_eq!(s_compressed_lookup_finalize(&mut lookup), 0);

        assert_eq!(lookup.longest_chain, 5);
        assert_eq!(lookup.pv_array_bts, crate::stat::PV_ARRAY_BTS);
        assert_ne!(lookup.pv[3 >> lookup.pv_array_bts] & (1u32 << (3 & 31)), 0);
        assert_ne!(
            lookup.pv[64 >> lookup.pv_array_bts] & (1u32 << (64 & 31)),
            0
        );
        assert_ne!(
            lookup.pv[129 >> lookup.pv_array_bts] & (1u32 << (129 & 31)),
            0
        );
        assert_eq!(lookup.pv[0] & (1u32 << 4), 0);
    }

    #[test]
    fn compressed_scan_subject_copies_inline_hits_after_pv_test() {
        let index = 1 + 15 * 2 + 225 * 3 + 3_375 * 4 + 50_625 * 5;
        let mut lookup = BlastCompressedAaLookupTable::new(5, 15, index + 1);
        assert_eq!(
            s_compressed_lookup_add_encoded(&mut lookup, &[1, 2, 3, 4, 5], 11),
            0
        );
        assert_eq!(
            s_compressed_lookup_add_encoded(&mut lookup, &[1, 2, 3, 4, 5], 17),
            0
        );
        assert_eq!(s_compressed_lookup_finalize(&mut lookup), 0);

        let mut hits = Vec::new();
        s_blast_compressed_aa_scan_subject_into(
            &lookup,
            &[9, 1, 2, 3, 4, 5, 8],
            8,
            None,
            &mut hits,
        );
        assert_eq!(
            hits,
            vec![
                crate::lookup::OffsetPair {
                    query_offset: 11,
                    subject_offset: 1,
                },
                crate::lookup::OffsetPair {
                    query_offset: 17,
                    subject_offset: 1,
                },
            ]
        );
    }

    #[test]
    fn compressed_scan_subject_pauses_when_cell_exceeds_capacity() {
        // NCBI s_BlastCompressedAaScanSubject (blast_aascan.c:340-405) copies a
        // cell's hits only when they all fit; on overflow it sets the resume
        // cursor and returns WITHOUT a partial copy. With a 2-hit cell and
        // array_size=1, no hits are emitted.
        let index = 1 + 15 * 2 + 225 * 3 + 3_375 * 4 + 50_625 * 5;
        let mut lookup = BlastCompressedAaLookupTable::new(5, 15, index + 1);
        for query_offset in [11, 17] {
            assert_eq!(
                s_compressed_lookup_add_encoded(&mut lookup, &[1, 2, 3, 4, 5], query_offset),
                0
            );
        }
        assert_eq!(s_compressed_lookup_finalize(&mut lookup), 0);

        let mut hits = Vec::new();
        let n = s_blast_compressed_aa_scan_subject_into(
            &lookup,
            &[9, 1, 2, 3, 4, 5, 8],
            1,
            None,
            &mut hits,
        );
        assert_eq!(n, 0);
        assert!(hits.is_empty());
    }

    #[test]
    fn compressed_scan_subject_batched_completes_dense_repetitive_hits() {
        let index = 1 + 15 * 2 + 225 * 3 + 3_375 * 4 + 50_625 * 5;
        let mut lookup = BlastCompressedAaLookupTable::new(5, 15, index + 1);
        assert_eq!(
            s_compressed_lookup_add_encoded(&mut lookup, &[1, 2, 3, 4, 5], 42),
            0
        );
        assert_eq!(s_compressed_lookup_finalize(&mut lookup), 0);

        let mut subject = Vec::new();
        for _ in 0..20 {
            subject.extend_from_slice(&[1, 2, 3, 4, 5]);
        }

        let mut single_batch = Vec::new();
        let truncated =
            s_blast_compressed_aa_scan_subject_into(&lookup, &subject, 4, None, &mut single_batch);
        assert_eq!(truncated, 4);

        let mut batched = Vec::new();
        let mut collected = Vec::new();
        let total = s_blast_compressed_aa_scan_subject_batched(
            &lookup,
            &subject,
            4,
            &mut batched,
            |batch| {
                collected.extend_from_slice(batch);
            },
        );

        assert_eq!(total, 20);
        assert_eq!(collected.len(), 20);
        assert_eq!(collected.first().map(|p| p.subject_offset), Some(0));
        assert_eq!(collected.last().map(|p| p.subject_offset), Some(95));
        assert!(collected.iter().all(|p| p.query_offset == 42));
    }

    #[test]
    fn test_protein_scan_identical() {
        let m = simple_matrix();
        let query = vec![1u8, 2, 3, 4, 5, 6, 7, 8];
        let subject = query.clone();
        let hits = blast_choose_protein_scan_subject(&query, &subject, &m, 3, 11.0, 20);
        assert!(!hits.is_empty(), "Should find hits for identical sequences");
        // Best hit should cover full length.
        let best = &hits[0];
        assert_eq!(best.score, 32); // 8 * 4
    }

    #[test]
    fn test_protein_scan_no_match() {
        let m = simple_matrix();
        // Query all 1s, subject all 2s — word score = -1*3 = -3 < threshold 11
        let query = vec![1u8, 1, 1, 1, 1];
        let subject = vec![2u8, 2, 2, 2, 2];
        let hits = blast_choose_protein_scan_subject(&query, &subject, &m, 3, 11.0, 20);
        assert!(
            hits.is_empty(),
            "Should find no hits for unrelated sequences"
        );
    }

    #[test]
    fn test_protein_scan_short_sequences() {
        let m = simple_matrix();
        let query = vec![1u8, 2];
        let subject = vec![1u8, 2, 3];
        let hits = blast_choose_protein_scan_subject(&query, &subject, &m, 3, 11.0, 20);
        assert!(
            hits.is_empty(),
            "Sequences shorter than word_size yield no hits"
        );
    }

    #[test]
    fn test_batch_scan_matches_single_two_hit_scanner() {
        let m = simple_matrix();
        let q1 = vec![1u8, 2, 3, 4, 5, 6, 7, 8, 9];
        let q2 = vec![9u8, 8, 7, 6, 5, 4, 3, 2, 1];
        let subject = vec![0u8, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 9, 8, 7, 6, 5, 4, 3, 2, 1];
        let t1 = ProteinLookupTable::build(&q1, 3, &m, 11.0);
        let t2 = ProteinLookupTable::build(&q2, 3, &m, 11.0);
        let tables = vec![&t1, &t2];
        let queries = vec![q1.as_slice(), q2.as_slice()];
        let merged = merge_pv(&tables);

        let batch = s_blast_aa_scan_subject(&queries, &tables, &merged, &subject, &m, 20);
        let singles = [
            protein_scan_with_table(&q1, &subject, &m, &t1, 20),
            protein_scan_with_table(&q2, &subject, &m, &t2, 20),
        ];

        assert_eq!(batch.len(), 2);
        for (query_idx, hits) in batch {
            let expected = &singles[query_idx];
            assert_eq!(hits.len(), expected.len());
            for (hit, exp) in hits.iter().zip(expected.iter()) {
                assert_eq!(
                    (
                        hit.query_start,
                        hit.query_end,
                        hit.subject_start,
                        hit.subject_end,
                        hit.score,
                        hit.num_ident,
                    ),
                    (
                        exp.query_start,
                        exp.query_end,
                        exp.subject_start,
                        exp.subject_end,
                        exp.score,
                        exp.num_ident,
                    )
                );
            }
        }
    }

    #[test]
    fn batch_scan_merged_pv_no_hit_prefilter_preserves_empty_result() {
        let m = simple_matrix();
        let q1 = vec![1u8, 2, 3, 4, 5, 6, 7, 8, 9];
        let q2 = vec![9u8, 8, 7, 6, 5, 4, 3, 2, 1];
        let subject = vec![27u8; 32];
        let t1 = ProteinLookupTable::build(&q1, 3, &m, 12.0);
        let t2 = ProteinLookupTable::build(&q2, 3, &m, 12.0);
        let tables = vec![&t1, &t2];
        let queries = vec![q1.as_slice(), q2.as_slice()];
        let merged = merge_pv(&tables);

        let batch = s_blast_aa_scan_subject(&queries, &tables, &merged, &subject, &m, 20);
        let singles = [
            protein_scan_with_table(&q1, &subject, &m, &t1, 20),
            protein_scan_with_table(&q2, &subject, &m, &t2, 20),
        ];

        assert!(batch.is_empty());
        assert!(singles.iter().all(Vec::is_empty));
    }

    #[test]
    fn batch_scan_ignores_merged_pv_that_is_not_a_superset() {
        let m = simple_matrix();
        let query = vec![1u8, 2, 3, 4, 5, 6, 7, 8, 9];
        let subject = query.clone();
        let table = ProteinLookupTable::build(&query, 3, &m, 11.0);
        let tables = vec![&table];
        let queries = vec![query.as_slice()];
        let invalid_merged = vec![0u64; table.pv.len()];

        let batch = s_blast_aa_scan_subject(&queries, &tables, &invalid_merged, &subject, &m, 20);
        let single = protein_scan_with_table(&query, &subject, &m, &table, 20);

        assert!(!single.is_empty());
        assert_eq!(batch.len(), 1);
        assert_eq!(batch[0].1.len(), single.len());
    }

    #[test]
    fn test_lookup_table_pv_array_exact_match() {
        let m = simple_matrix();
        // Query: [1, 2, 3, 4, 5], word_size=3, threshold=12 (exact match only)
        let query = vec![1u8, 2, 3, 4, 5];
        let table = ProteinLookupTable::build(&query, 3, &m, 12.0);
        // Should have entries for words [1,2,3], [2,3,4], [3,4,5]
        let hashes = [
            word_hash(&[1, 2, 3], AA_SIZE),
            word_hash(&[2, 3, 4], AA_SIZE),
            word_hash(&[3, 4, 5], AA_SIZE),
        ];
        for (i, &hash) in hashes.iter().enumerate() {
            let cell = &table.backbone[hash];
            assert!(
                cell.num_used > 0,
                "Word at position {} should have entries",
                i
            );
            let hits: &[i32] = if (cell.num_used as usize) <= HITS_PER_CELL {
                &cell.entries[..cell.num_used as usize]
            } else {
                let c = cell.entries[0] as usize;
                &table.overflow[c..c + cell.num_used as usize]
            };
            assert!(
                hits.contains(&(i as i32)),
                "Word at position {} should map to query offset {}",
                i,
                i
            );
            // PV bit should be set
            assert_ne!(
                table.pv[hash >> 6] & (1u64 << (hash & 63)),
                0,
                "PV bit should be set for word at position {}",
                i
            );
        }
    }

    #[test]
    fn test_lookup_table_no_neighbors_high_threshold() {
        let m = simple_matrix();
        // threshold=13 > max possible score (4*3=12), so no NEIGHBORHOOD
        // word can meet the threshold. NCBI's `s_AddWordHits`
        // (`blast_aalookup.c:504-509`) explicitly inserts the EXACT query
        // word when its self-score is below the threshold (otherwise the
        // exact match would be missed entirely). With 3 query positions
        // (5 - 3 + 1 = 3 windows), we expect exactly 3 exact-word entries.
        let query = vec![1u8, 2, 3, 4, 5];
        let table = ProteinLookupTable::build(&query, 3, &m, 13.0);
        let total_entries: usize = table.backbone.iter().map(|c| c.num_used as usize).sum();
        assert_eq!(
            total_entries, 3,
            "Only the 3 exact query words should be inserted (no neighbors), got {}",
            total_entries
        );
    }

    #[test]
    fn test_lookup_table_neighbors_low_threshold() {
        let m = simple_matrix();
        // threshold=7 allows words with 2 matches + 1 mismatch (4+4-1=7)
        let query = vec![1u8, 2, 3];
        let table = ProteinLookupTable::build(&query, 3, &m, 7.0);
        let total_entries: usize = table.backbone.iter().map(|c| c.num_used as usize).sum();
        assert!(
            total_entries > 1,
            "Low threshold should produce neighboring words (got {})",
            total_entries
        );
    }

    #[test]
    fn test_protein_scan_with_table_reuse() {
        let m = crate::matrix::BLOSUM62;
        let query = encode_ncbistdaa_sequence(b"AGIMVAGI");
        let table = ProteinLookupTable::build(&query, 3, &m, 11.0);
        // Scan two different subjects with same table
        let subject1 = query.clone();
        let subject2 = vec![2u8, 3, 4, 5, 6, 7]; // completely different
        let hits1 = protein_scan_with_table(&query, &subject1, &m, &table, 40);
        let _hits2 = protein_scan_with_table(&query, &subject2, &m, &table, 40);
        assert!(
            !hits1.is_empty(),
            "Should find hits for identical sequences"
        );
        // hits2 may or may not be empty depending on BLOSUM62 neighborhood
    }

    #[test]
    fn two_hit_cutoff_gate_drops_sub_cutoff_hsps() {
        // NCBI `s_BlastAaWordFinder_TwoHit` (`aa_ungapped.c:588`) only saves an
        // init HSP when `score >= cutoffs->cutoff_score`. With cutoff 1 every
        // positive HSP is saved (historic `> 0` behavior); raising the cutoff
        // above the best HSP score must drop all of them, while the gate strictly
        // shrinks (never grows) the saved set.
        let m = crate::matrix::BLOSUM62;
        let query = encode_ncbistdaa_sequence(b"AGIMVAGIKLMNPQRS");
        let subject = query.clone();
        let table = ProteinLookupTable::build(&query, 3, &m, 11.0);

        let mut diag = Vec::new();
        let all = blast_aa_word_finder(
            &query,
            &subject,
            &m,
            &table,
            40,
            TWO_HIT_WINDOW,
            1,
            &mut diag,
        )
        .expect("two-hit scan");
        assert!(
            !all.is_empty(),
            "expected at least one positive ungapped HSP"
        );
        let best = all.iter().map(|h| h.score).max().unwrap();

        // A cutoff equal to the best score keeps only the top HSP(s); a cutoff
        // above it drops everything.
        let mut diag2 = Vec::new();
        let at_best = blast_aa_word_finder(
            &query,
            &subject,
            &m,
            &table,
            40,
            TWO_HIT_WINDOW,
            best,
            &mut diag2,
        )
        .expect("two-hit scan");
        assert!(at_best.iter().all(|h| h.score >= best));
        assert!(at_best.len() <= all.len());

        let mut diag3 = Vec::new();
        let above = blast_aa_word_finder(
            &query,
            &subject,
            &m,
            &table,
            40,
            TWO_HIT_WINDOW,
            best + 1,
            &mut diag3,
        )
        .expect("two-hit scan");
        assert!(
            above.is_empty(),
            "cutoff above best score must drop all HSPs"
        );
    }

    #[test]
    fn one_hit_cutoff_gate_drops_sub_cutoff_hsps() {
        let m = simple_matrix();
        let query = vec![1u8, 2, 3, 4, 5];
        let subject = query.clone();
        let table = ProteinLookupTable::build(&query, 3, &m, 12.0);

        let mut diag = Vec::new();
        let all = blast_aa_word_finder(&query, &subject, &m, &table, 20, 0, 1, &mut diag)
            .expect("one-hit scan");
        assert_eq!(all.len(), 1);
        let best = all[0].score;

        let mut diag2 = Vec::new();
        let at_best = blast_aa_word_finder(&query, &subject, &m, &table, 20, 0, best, &mut diag2)
            .expect("one-hit scan");
        assert_eq!(at_best.len(), 1);
        assert_eq!(at_best[0].score, best);

        let mut diag3 = Vec::new();
        let above = blast_aa_word_finder(&query, &subject, &m, &table, 20, 0, best + 1, &mut diag3)
            .expect("one-hit scan");
        assert!(
            above.is_empty(),
            "cutoff above best score must drop all HSPs"
        );
    }

    #[test]
    fn compressed_one_hit_cutoff_gate_drops_sub_cutoff_hsps() {
        let m = crate::matrix::BLOSUM62;
        let query = encode_ncbistdaa_sequence(b"AGIMVA");
        let subject = query.clone();
        let table = ProteinLookupTable::build(&query, 6, &m, 11.0);
        assert!(table.compressed.is_some());

        let mut diag = Vec::new();
        let all = blast_aa_word_finder(&query, &subject, &m, &table, 40, 0, 1, &mut diag)
            .expect("compressed one-hit scan");
        assert!(!all.is_empty());
        let best = all.iter().map(|hit| hit.score).max().unwrap();

        let mut diag2 = Vec::new();
        let above = blast_aa_word_finder(&query, &subject, &m, &table, 40, 0, best + 1, &mut diag2)
            .expect("compressed one-hit scan");
        assert!(
            above.is_empty(),
            "cutoff above best score must drop all HSPs"
        );
    }

    #[test]
    fn blast_aa_word_finder_dispatches_default_two_hit_window() {
        let m = simple_matrix();
        let query = vec![1u8, 2, 3, 4, 5, 6, 7, 8, 9];
        let subject = vec![0u8, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0];
        let table = ProteinLookupTable::build(&query, 3, &m, 12.0);
        let mut direct_diag = Vec::new();
        let mut dispatch_diag = Vec::new();

        let direct = s_blast_aa_word_finder_two_hit(
            &query,
            &subject,
            &m,
            &table,
            20,
            TWO_HIT_WINDOW,
            &mut direct_diag,
        );
        let dispatch = blast_aa_word_finder(
            &query,
            &subject,
            &m,
            &table,
            20,
            TWO_HIT_WINDOW,
            1,
            &mut dispatch_diag,
        )
        .expect("default protein window should use the represented two-hit scanner");

        assert_eq!(dispatch.len(), direct.len());
        for (hit, exp) in dispatch.iter().zip(direct.iter()) {
            assert_eq!(
                (
                    hit.query_start,
                    hit.query_end,
                    hit.subject_start,
                    hit.subject_end,
                    hit.score,
                    hit.num_ident,
                ),
                (
                    exp.query_start,
                    exp.query_end,
                    exp.subject_start,
                    exp.subject_end,
                    exp.score,
                    exp.num_ident,
                )
            );
        }
    }

    #[test]
    fn blast_aa_word_finder_honors_non_default_window() {
        // Non-default two-hit windows are now threaded through to the scanner
        // (previously rejected as BlastAaWordFinderUnsupported::WindowSize). The
        // dispatch result must match a direct two-hit scan using the same window.
        let m = simple_matrix();
        let query = vec![1u8, 2, 3, 4, 5, 6, 7, 8, 9];
        let subject = vec![0u8, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0];
        let table = ProteinLookupTable::build(&query, 3, &m, 12.0);
        let mut direct_diag = Vec::new();
        let mut dispatch_diag = Vec::new();

        let direct =
            s_blast_aa_word_finder_two_hit(&query, &subject, &m, &table, 20, 12, &mut direct_diag);
        let dispatch =
            blast_aa_word_finder(&query, &subject, &m, &table, 20, 12, 1, &mut dispatch_diag)
                .expect("non-default window is now supported");
        assert_eq!(dispatch.len(), direct.len());
        for (hit, exp) in dispatch.iter().zip(direct.iter()) {
            assert_eq!(
                (hit.query_start, hit.subject_start, hit.score),
                (exp.query_start, exp.subject_start, exp.score)
            );
        }
    }

    #[test]
    fn blast_aa_word_finder_dispatches_one_hit_window() {
        let m = simple_matrix();
        let query = vec![1u8, 2, 3, 4, 5];
        let subject = vec![9u8, 1, 2, 3, 4, 5, 9];
        let table = ProteinLookupTable::build(&query, 3, &m, 12.0);
        let mut direct_diag = Vec::new();
        let mut dispatch_diag = Vec::new();

        let direct =
            s_blast_aa_word_finder_one_hit(&query, &subject, &m, &table, 20, &mut direct_diag);
        let dispatch =
            blast_aa_word_finder(&query, &subject, &m, &table, 20, 0, 1, &mut dispatch_diag)
                .expect("zero window should use the represented one-hit scanner");

        assert!(!dispatch.is_empty());
        assert_eq!(dispatch.len(), direct.len());
        for (hit, exp) in dispatch.iter().zip(direct.iter()) {
            assert_eq!(
                (
                    hit.query_start,
                    hit.query_end,
                    hit.subject_start,
                    hit.subject_end,
                    hit.score,
                    hit.num_ident,
                ),
                (
                    exp.query_start,
                    exp.query_end,
                    exp.subject_start,
                    exp.subject_end,
                    exp.score,
                    exp.num_ident,
                )
            );
        }
    }

    #[test]
    fn blast_aa_word_finder_with_side_channels_fills_init_list_and_stats() {
        let m = simple_matrix();
        let query = vec![1u8, 2, 3, 4, 5, 6];
        let subject = query.clone();
        let table = ProteinLookupTable::build(&query, 3, &m, 12.0);
        let mut diag = Vec::new();
        let mut init_hitlist = crate::extend::InitHitList::new();
        let mut stats = crate::diagnostics::UngappedStats::default();

        let scan = blast_aa_word_finder_with_side_channels(
            &query,
            &subject,
            &m,
            &table,
            20,
            0,
            1,
            &mut diag,
            Some(&mut init_hitlist),
            Some(&mut stats),
        )
        .expect("one-hit side-channel scan");

        assert_eq!(scan.hits.len(), 1);
        assert!(scan.total_hits >= scan.hits_extended);
        assert_eq!(scan.hits_extended, 1);
        assert_eq!(init_hitlist.total(), scan.hits.len());
        assert_eq!(init_hitlist.hits[0].query_offset, 0);
        assert_eq!(init_hitlist.hits[0].subject_offset, 0);
        assert_eq!(
            init_hitlist.hits[0]
                .ungapped_data
                .as_ref()
                .map(|data| { (data.q_start, data.s_start, data.length, data.score) }),
            Some((
                scan.hits[0].query_start as i32,
                scan.hits[0].subject_start as i32,
                scan.hits[0].align_length,
                scan.hits[0].score,
            ))
        );
        assert_eq!(stats.lookup_hits, scan.total_hits as i64);
        assert_eq!(stats.num_seqs_lookup_hits, 1);
        assert_eq!(stats.init_extends, scan.hits_extended);
        assert_eq!(stats.good_init_extends, scan.hits.len() as i32);
        assert_eq!(stats.num_seqs_passed, 1);
    }

    #[test]
    fn blast_aa_word_finder_side_channels_preserve_two_hit_seed_offsets() {
        let m = simple_matrix();
        let query = vec![1u8, 2, 3, 4, 5, 6, 7, 8, 9];
        let subject = vec![0u8, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0];
        let table = ProteinLookupTable::build(&query, 3, &m, 12.0);
        let mut diag = Vec::new();
        let mut init_hitlist = crate::extend::InitHitList::new();
        let mut stats = crate::diagnostics::UngappedStats::default();

        let scan = blast_aa_word_finder_with_side_channels(
            &query,
            &subject,
            &m,
            &table,
            20,
            TWO_HIT_WINDOW,
            1,
            &mut diag,
            Some(&mut init_hitlist),
            Some(&mut stats),
        )
        .expect("two-hit side-channel scan");

        assert!(!scan.hits.is_empty());
        assert_eq!(init_hitlist.total(), scan.hits.len());
        assert!(
            init_hitlist.hits.iter().any(|init| {
                let data = init.ungapped_data.as_ref().expect("ungapped data");
                init.query_offset > data.q_start && init.subject_offset > data.s_start
            }),
            "two-hit side channels should preserve triggering seed offsets, not only HSP starts"
        );
        assert_eq!(stats.good_init_extends, scan.hits.len() as i32);
    }

    #[test]
    fn one_hit_scanner_suppresses_covered_diagonal_words() {
        let m = simple_matrix();
        let query = vec![1u8, 2, 3, 4, 5, 6];
        let subject = query.clone();
        let table = ProteinLookupTable::build(&query, 3, &m, 12.0);
        let mut diag = Vec::new();

        let hits = s_blast_aa_word_finder_one_hit(&query, &subject, &m, &table, 20, &mut diag);

        assert_eq!(hits.len(), 1);
        assert_eq!(
            (
                hits[0].query_start,
                hits[0].query_end,
                hits[0].subject_start,
                hits[0].subject_end,
            ),
            (0, 6, 0, 6)
        );
    }

    #[test]
    fn two_hit_scanner_saves_left_only_hsp_when_first_hit_not_reached() {
        let m = simple_matrix();
        let query = vec![1u8, 2, 3, 4, 4, 4, 4, 4, 4, 4, 5, 6, 7];
        let subject = vec![1u8, 2, 3, 8, 8, 8, 8, 8, 8, 8, 5, 6, 7];
        let table = ProteinLookupTable::build(&query, 3, &m, 12.0);

        let hits = protein_scan_with_table(&query, &subject, &m, &table, 1);

        assert!(
            hits.iter().any(|hit| {
                (
                    hit.query_start,
                    hit.query_end,
                    hit.subject_start,
                    hit.subject_end,
                    hit.score,
                    hit.num_ident,
                    hit.align_length,
                ) == (10, 13, 10, 13, 12, 3, 3)
            }),
            "left-only HSP from the second word hit should be saved: {hits:?}"
        );
    }

    #[test]
    fn test_word_hash_deterministic() {
        // Same input should always produce same hash
        let h1 = word_hash(&[5, 10, 15], AA_SIZE);
        let h2 = word_hash(&[5, 10, 15], AA_SIZE);
        assert_eq!(h1, h2);
        // Different input should (generally) produce different hash
        let h3 = word_hash(&[5, 10, 16], AA_SIZE);
        assert_ne!(h1, h3);
    }

    #[test]
    fn test_word_hash_range() {
        // All hashes should be in [0, 2^(word_size*CHARSIZE))
        let max_hash = 1usize << (3 * CHARSIZE); // 32768
        for a in 0..AA_SIZE {
            for b in 0..AA_SIZE {
                for c in 0..AA_SIZE {
                    let h = word_hash(&[a as u8, b as u8, c as u8], AA_SIZE);
                    assert!(
                        h < max_hash,
                        "Hash {} out of range for [{},{},{}]",
                        h,
                        a,
                        b,
                        c
                    );
                }
            }
        }
    }

    #[test]
    fn test_protein_gapped_scan_with_table_finds_alignment() {
        let m = crate::matrix::BLOSUM62;
        // Two related protein sequences
        let query = encode_ncbistdaa_sequence(b"AGIMVAGIMV");
        let subject = query.clone();
        let table = ProteinLookupTable::build(&query, 3, &m, 11.0);
        let hits = protein_gapped_scan_with_table(
            &query, &subject, &m, &table, 40, // ungap_x_dropoff
            11, 1,   // gap_open, gap_extend
            260, // gap_x_dropoff
            40,  // ungap_cutoff
        );
        assert!(
            !hits.is_empty(),
            "Gapped scan should find alignment for identical sequences"
        );
        let best = &hits[0];
        assert!(best.score > 0);
        assert_eq!(best.query_start, 0);
    }

    #[test]
    fn test_protein_gapped_scan_uses_ncbi_window_start_not_midpoint() {
        let m = simple_matrix();
        let query = vec![1u8, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14];
        let subject = query.clone();
        let table = ProteinLookupTable::build(&query, 3, &m, 12.0);
        let hits = protein_gapped_scan_with_table(
            &query, &subject, &m, &table, 20, // ungap_x_dropoff
            5, 1,   // gap_open, gap_extend
            100, // gap_x_dropoff
            1,   // ungap_cutoff
        );

        assert!(!hits.is_empty());
        let expected = blast_get_start_for_gapped_alignment(
            &query,
            &subject,
            0,
            query.len(),
            0,
            subject.len(),
            &m,
        );
        assert_eq!((hits[0].gapped_start_q, hits[0].gapped_start_s), expected);
        assert_ne!(
            (hits[0].gapped_start_q, hits[0].gapped_start_s),
            (query.len() / 2, subject.len() / 2)
        );
    }
}
