//! Rust equivalent of blast_util.c — BLAST utility functions.

/// Translate a codon using a genetic code.
/// Takes 3 NCBI2na bases (2-bit unambiguous A=0/C=1/G=2/T=3) and
/// returns the amino acid in NCBIstdaa encoding.
///
/// Companion to NCBI's `s_CodonToAA` (`blast_util.c:369`), but simpler:
/// NCBI's version takes NCBI4na bases and iterates over all possible
/// combinations (returning X when a codon with an ambiguity code
/// translates to multiple amino acids). This Rust port assumes the
/// caller has already resolved ambiguities or is using NCBI2na input.
pub fn translate_codon(b1: u8, b2: u8, b3: u8, genetic_code: &[u8; 64]) -> u8 {
    let idx = ((b1 & 3) as usize) * 16 + ((b2 & 3) as usize) * 4 + (b3 & 3) as usize;
    genetic_code[idx]
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

/// Translate a nucleotide sequence to protein in all 6 reading frames.
/// Returns Vec of (frame, protein_sequence) where frame is 1,2,3,-1,-2,-3.
pub fn six_frame_translation(nuc_seq: &[u8], genetic_code: &[u8; 64]) -> Vec<(i32, Vec<u8>)> {
    let mut results = Vec::new();
    let len = nuc_seq.len();

    // Forward frames 1, 2, 3
    for frame_offset in 0..3 {
        let mut protein = Vec::new();
        let mut i = frame_offset;
        while i + 2 < len {
            protein.push(translate_codon(
                nuc_seq[i],
                nuc_seq[i + 1],
                nuc_seq[i + 2],
                genetic_code,
            ));
            i += 3;
        }
        results.push((frame_offset as i32 + 1, protein));
    }

    // Reverse complement for frames -1, -2, -3
    let rc: Vec<u8> = nuc_seq.iter().rev().map(|&b| 3 - (b & 3)).collect();
    for frame_offset in 0..3 {
        let mut protein = Vec::new();
        let mut i = frame_offset;
        while i + 2 < rc.len() {
            protein.push(translate_codon(rc[i], rc[i + 1], rc[i + 2], genetic_code));
            i += 3;
        }
        results.push((-(frame_offset as i32 + 1), protein));
    }

    results
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
        // Reverse strand: map RC protein coords to forward strand nucleotide coords.
        // Returns (low, high) on forward strand; caller uses qframe to indicate direction.
        let offset = (-frame) - 1;
        let rc_start = prot_start * 3 + offset;
        let rc_end = prot_end * 3 + offset - 1;
        (nuc_len - rc_end, nuc_len - rc_start)
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
        // ACGTACGTAC (10 bases)
        let seq = vec![0u8, 1, 2, 3, 0, 1, 2, 3, 0, 1];
        let frames = six_frame_translation(&seq, &STANDARD_GENETIC_CODE);
        assert_eq!(frames.len(), 6);
        assert_eq!(frames[0].0, 1); // frame +1
        assert_eq!(frames[3].0, -1); // frame -1
                                     // Frame +1: ACG TAC GTA → 3 codons
        assert_eq!(frames[0].1.len(), 3);
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
        // offset = 0
        // rc_start = 0*3 + 0 = 0, rc_end = 3*3 + 0 - 1 = 8
        // nuc: (30 - 8, 30 - 0 + 1) = (22, 31) — but 31 > 30, so:
        // Actually: start on forward = nuc_len - rc_end = 30 - 8 = 22
        //           end on forward = nuc_len - rc_start + 1 = 30 - 0 + 1 = 31
        // Hmm, that's > nuc_len. Let me reconsider.
        // For BLAST output on minus strand, coordinates go high→low.
        // A protein alignment at positions 0..3 on the reverse complement
        // means nucleotides rc[0..9] which maps to forward strand [21..30].
        // 30bp sequence, frame -1, protein pos 0..3 (0-based on rev-comp)
        // rc_start=0, rc_end=8, fwd: (30-8, 30-0) = (22, 30)
        let (start, end) = protein_to_nuc_coords(0, 3, -1, 30);
        assert_eq!(start, 22);
        assert_eq!(end, 30);
    }

    #[test]
    fn test_protein_to_nuc_coords_reverse_frame2() {
        // 30bp, frame -2, protein pos 1..4
        // offset = 1
        // rc_start = 1*3 + 1 = 4, rc_end = 4*3 + 1 - 1 = 12
        // nuc: (30 - 12, 30 - 4) = (18, 26)
        let (start, end) = protein_to_nuc_coords(1, 4, -2, 30);
        assert_eq!(start, 18);
        assert_eq!(end, 26);
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
        // Sequence shorter than 3 bases: no codons in any frame
        let seq = vec![0u8, 1]; // AC
        let frames = six_frame_translation(&seq, &STANDARD_GENETIC_CODE);
        assert_eq!(frames.len(), 6);
        for (_, protein) in &frames {
            assert_eq!(protein.len(), 0);
        }
    }

    #[test]
    fn test_six_frame_translation_exact_codon() {
        // Exactly 3 bases: one codon in frame +1, zero in +2 and +3
        let seq = vec![0u8, 3, 2]; // ATG = Met
        let frames = six_frame_translation(&seq, &STANDARD_GENETIC_CODE);
        assert_eq!(frames[0].1.len(), 1); // frame +1
        assert_eq!(frames[0].1[0], 12); // Met = 12
        assert_eq!(frames[1].1.len(), 0); // frame +2
        assert_eq!(frames[2].1.len(), 0); // frame +3
    }

    #[test]
    fn test_six_frame_frame_numbers() {
        let seq = vec![0u8, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3]; // 12 bases
        let frames = six_frame_translation(&seq, &STANDARD_GENETIC_CODE);
        assert_eq!(frames[0].0, 1);
        assert_eq!(frames[1].0, 2);
        assert_eq!(frames[2].0, 3);
        assert_eq!(frames[3].0, -1);
        assert_eq!(frames[4].0, -2);
        assert_eq!(frames[5].0, -3);
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
        // offset = 2
        // rc_start = 0*3 + 2 = 2, rc_end = 2*3 + 2 - 1 = 7
        // nuc: (30 - 7, 30 - 2) = (23, 28)
        let (start, end) = protein_to_nuc_coords(0, 2, -3, 30);
        assert_eq!(start, 23);
        assert_eq!(end, 28);
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
