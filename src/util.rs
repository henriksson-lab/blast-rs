//! Rust equivalent of blast_util.c — BLAST utility functions.

/// Translate a codon using a genetic code (NCBI2na fast path).
/// Takes 3 NCBI2na bases (2-bit unambiguous A=0/C=1/G=2/T=3) and
/// returns the amino acid in NCBIstdaa encoding.
///
/// For NCBI4na input with ambiguity resolution, use [`codon_to_aa`] —
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
/// agree, otherwise returns NCBIstdaa `X` (21). Bytes > 15 (malformed
/// or FENCE_SENTRY) also yield `X`.
///
/// NCBI's C iterates `mapping[4] = {T=8, C=2, A=1, G=4}` so its
/// `codes` table is TCAG-ordered. Our [`STANDARD_GENETIC_CODE`] is
/// ACGT-ordered, so `MAPPING` here is `[A=1, C=2, G=4, T=8]`.
pub fn codon_to_aa(codon: [u8; 3], genetic_code: &[u8; 64]) -> u8 {
    const X_RESIDUE: u8 = 21;
    const MAPPING: [u8; 4] = [1, 2, 4, 8];

    if (codon[0] | codon[1] | codon[2]) > 15 {
        return X_RESIDUE;
    }

    let mut aa: u8 = 0;
    for i in 0..4 {
        if codon[0] & MAPPING[i] != 0 {
            let index0 = i * 16;
            for j in 0..4 {
                if codon[1] & MAPPING[j] != 0 {
                    let index1 = index0 + j * 4;
                    for k in 0..4 {
                        if codon[2] & MAPPING[k] != 0 {
                            let taa = genetic_code[index1 + k];
                            if aa == 0 {
                                aa = taa;
                            } else if taa != aa {
                                return X_RESIDUE;
                            }
                        }
                    }
                }
            }
        }
    }
    aa
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
pub fn blast_context_to_frame_blastx(context_number: u32) -> i32 {
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

/// Port of `GetReverseNuclSequence` (`blast_util.c:807`). Builds the
/// reverse-complement of an NCBI4na buffer with NULLB sentinels at
/// positions 0 and length+1, matching the C layout used by
/// `BLAST_GetTranslation` (`query_seq_rev + 1` skips the leading NULLB).
pub fn get_reverse_nucl_sequence(sequence: &[u8], length: usize) -> Vec<u8> {
    // C: rev_sequence = malloc(length + 2);
    //    rev_sequence[0] = rev_sequence[length+1] = NULLB;
    let mut rev = vec![NULLB; length + 2];
    // Conversion table from forward to reverse strand residue (NCBI4na).
    // (`blast_util.c:814-819`. The C source comment says "blastna
    // encoding" but the table values are NCBI4na bit masks; see the
    // identities A(1)↔T(8), C(2)↔G(4), N(15)↔N(15).)
    const CONVERSION_TABLE: [u8; 16] = [0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15];
    for index in 0..length {
        let b = sequence[index];
        rev[length - index] = if b == FENCE_SENTRY {
            FENCE_SENTRY
        } else {
            CONVERSION_TABLE[b as usize & 0x0f]
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
    if nt_length >= 2 {
        let mut index = start;
        while index < nt_length - 2 {
            let codon = [nucl_seq[index], nucl_seq[index + 1], nucl_seq[index + 2]];
            let residue = codon_to_aa(codon, genetic_code);
            // C: if (IS_residue(residue) || residue == FENCE_SENTRY)
            if is_residue(residue) || residue == FENCE_SENTRY {
                prot_seq[index_prot] = residue;
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
pub fn blast_get_all_translations_ncbi4na(
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
        let frame = blast_context_to_frame_blastx(context as u32);
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
    fn test_translate_codon_met() {
        // ATG = methionine (M = 12 in NCBIstdaa)
        // A=0, T=3, G=2 → idx = 0*16 + 3*4 + 2 = 14
        assert_eq!(translate_codon(0, 3, 2, &STANDARD_GENETIC_CODE), 12);
    }

    #[test]
    fn test_translate_codon_stop() {
        // TAA = stop (* = 25 in NCBIstdaa)
        // T=3, A=0, A=0 → idx = 3*16 + 0*4 + 0 = 48
        assert_eq!(translate_codon(3, 0, 0, &STANDARD_GENETIC_CODE), 25);
    }

    #[test]
    fn test_six_frame() {
        // ACGTACGTAC (10 bases) in NCBI4na: A=1, C=2, G=4, T=8
        let seq = vec![1u8, 2, 4, 8, 1, 2, 4, 8, 1, 2];
        let (_buf, offsets) =
            blast_get_all_translations_ncbi4na(&seq, seq.len(), &STANDARD_GENETIC_CODE);
        // Frame +1 covers ACG TAC GTA → 3 residues.
        assert_eq!(offsets[1] - offsets[0] - 1, 3);
        // Frame names are produced in order 1, 2, 3, -1, -2, -3.
        assert_eq!(blast_context_to_frame_blastx(0), 1);
        assert_eq!(blast_context_to_frame_blastx(3), -1);
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
        // TAA (25): T=3, A=0, A=0 -> idx=48
        assert_eq!(translate_codon(3, 0, 0, &STANDARD_GENETIC_CODE), 25);
        // TAG (25): T=3, A=0, G=2 -> idx=50
        assert_eq!(translate_codon(3, 0, 2, &STANDARD_GENETIC_CODE), 25);
        // TGA (25): T=3, G=2, A=0 -> idx=56
        assert_eq!(translate_codon(3, 2, 0, &STANDARD_GENETIC_CODE), 25);
    }

    #[test]
    fn test_codon_to_aa_unambiguous() {
        // ATG -> M (12): NCBI4na A=1, T=8, G=4
        assert_eq!(codon_to_aa([1, 8, 4], &STANDARD_GENETIC_CODE), 12);
        // TAA -> * (25)
        assert_eq!(codon_to_aa([8, 1, 1], &STANDARD_GENETIC_CODE), 25);
        // GCT -> A (1)
        assert_eq!(codon_to_aa([4, 2, 8], &STANDARD_GENETIC_CODE), 1);
    }

    #[test]
    fn test_codon_to_aa_ambiguity_resolves_to_x() {
        // CTN: CTA=L, CTC=L, CTG=L, CTT=L — all leucine, so resolves to L (11).
        // CTN with N=15
        assert_eq!(codon_to_aa([2, 8, 15], &STANDARD_GENETIC_CODE), 11);

        // ATN: ATA=I, ATC=I, ATG=M, ATT=I — mixed, must return X (21).
        assert_eq!(codon_to_aa([1, 8, 15], &STANDARD_GENETIC_CODE), 21);

        // YTA (Y=C|T=10): CTA=L, TTA=L — both L, so resolves to L (11).
        assert_eq!(codon_to_aa([10, 8, 1], &STANDARD_GENETIC_CODE), 11);

        // RTA (R=A|G=5): ATA=I, GTA=V — mixed, X.
        assert_eq!(codon_to_aa([5, 8, 1], &STANDARD_GENETIC_CODE), 21);
    }

    #[test]
    fn test_codon_to_aa_n_in_first_position_mixed() {
        // NCT (N=15): all four first-base options. ACT=T, CCT=P, GCT=A, TCT=S.
        // Mixed → X.
        assert_eq!(codon_to_aa([15, 2, 8], &STANDARD_GENETIC_CODE), 21);
    }

    #[test]
    fn test_codon_to_aa_malformed_returns_x() {
        // Any byte > 15 (e.g. FENCE_SENTRY) → X.
        assert_eq!(codon_to_aa([1, 8, 200], &STANDARD_GENETIC_CODE), 21);
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
    fn test_blast_get_all_translations_ncbi4na_offsets_and_frames() {
        let seq: Vec<u8> = vec![1, 8, 4, 4, 2, 8]; // ATGGCT
        let (buf, offsets) =
            blast_get_all_translations_ncbi4na(&seq, seq.len(), &STANDARD_GENETIC_CODE);
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
    fn test_genetic_code_2_vertebrate_mito() {
        // AGA -> * (stop, 25) instead of R
        assert_eq!(translate_codon(0, 2, 0, &GENETIC_CODE_2), 25);
        // AGG -> * (stop, 25) instead of R
        assert_eq!(translate_codon(0, 2, 2, &GENETIC_CODE_2), 25);
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
        assert_eq!(translate_codon(3, 2, 0, &GENETIC_CODE_6), 25);
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
    fn test_six_frame_translation_short_seq() {
        // Sequence shorter than 3 bases (NCBI4na: A=1, C=2): no codons in any frame.
        let seq = vec![1u8, 2]; // AC
        let (_buf, offsets) =
            blast_get_all_translations_ncbi4na(&seq, seq.len(), &STANDARD_GENETIC_CODE);
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
        let (buf, offsets) =
            blast_get_all_translations_ncbi4na(&seq, seq.len(), &STANDARD_GENETIC_CODE);
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
        assert_eq!(blast_context_to_frame_blastx(0), 1);
        assert_eq!(blast_context_to_frame_blastx(1), 2);
        assert_eq!(blast_context_to_frame_blastx(2), 3);
        assert_eq!(blast_context_to_frame_blastx(3), -1);
        assert_eq!(blast_context_to_frame_blastx(4), -2);
        assert_eq!(blast_context_to_frame_blastx(5), -3);
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
}
