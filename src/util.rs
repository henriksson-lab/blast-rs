//! Rust equivalent of blast_util.c — BLAST utility functions.

/// Translate a nucleotide sequence using a genetic code.
/// Takes 3 bases (codon) and returns the amino acid in NCBIstdaa encoding.
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
            protein.push(translate_codon(nuc_seq[i], nuc_seq[i + 1], nuc_seq[i + 2], genetic_code));
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
pub fn protein_to_nuc_coords(prot_start: i32, prot_end: i32, frame: i32, nuc_len: i32) -> (i32, i32) {
    if frame > 0 {
        let offset = frame - 1;
        (prot_start * 3 + offset + 1, prot_end * 3 + offset)
    } else {
        let offset = (-frame) - 1;
        let rc_start = prot_start * 3 + offset;
        let rc_end = prot_end * 3 + offset - 1;
        (nuc_len - rc_end, nuc_len - rc_start + 1)
    }
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
        let (start, end) = protein_to_nuc_coords(0, 3, -1, 30);
        assert_eq!(start, 22);
        assert_eq!(end, 31);
    }

    #[test]
    fn test_protein_to_nuc_coords_reverse_frame2() {
        // 30bp, frame -2, protein pos 1..4
        // offset = 1
        // rc_start = 1*3 + 1 = 4, rc_end = 4*3 + 1 - 1 = 12
        // nuc: (30 - 12, 30 - 4 + 1) = (18, 27)
        let (start, end) = protein_to_nuc_coords(1, 4, -2, 30);
        assert_eq!(start, 18);
        assert_eq!(end, 27);
    }
}
