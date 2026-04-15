//! Rust equivalent of blast_encoding.c
//! These tables will replace the C versions when the FFI layer is removed.

pub const BLASTNA_SIZE: usize = 16;
pub const BLASTAA_SIZE: usize = 28;

pub static NCBI4NA_TO_BLASTNA: [u8; BLASTNA_SIZE] =
    [15, 0, 1, 6, 2, 4, 9, 13, 3, 8, 5, 12, 7, 11, 10, 14];

pub static BLASTNA_TO_NCBI4NA: [u8; BLASTNA_SIZE] =
    [1, 2, 4, 8, 5, 10, 3, 12, 9, 6, 14, 13, 11, 7, 15, 0];

pub static BLASTNA_TO_IUPACNA: [i8; BLASTNA_SIZE] = [
    b'A' as i8, b'C' as i8, b'G' as i8, b'T' as i8, b'R' as i8, b'Y' as i8, b'M' as i8, b'K' as i8,
    b'W' as i8, b'S' as i8, b'B' as i8, b'D' as i8, b'H' as i8, b'V' as i8, b'N' as i8, b'-' as i8,
];

pub static NCBI4NA_TO_IUPACNA: [i8; BLASTNA_SIZE] = [
    b'-' as i8, b'A' as i8, b'C' as i8, b'M' as i8, b'G' as i8, b'R' as i8, b'S' as i8, b'V' as i8,
    b'T' as i8, b'W' as i8, b'Y' as i8, b'H' as i8, b'K' as i8, b'D' as i8, b'B' as i8, b'N' as i8,
];

pub static IUPACNA_TO_BLASTNA: [u8; 128] = [
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 0, 10, 1, 11, 15, 15, 2,
    12, 15, 15, 7, 15, 6, 14, 15, // @ABCDEFGHIJKLMNO
    15, 15, 4, 9, 3, 15, 13, 8, 15, 5, 15, 15, 15, 15, 15, 15, // PQRSTUVWXYZ[\]^_
    15, 0, 10, 1, 11, 15, 15, 2, 12, 15, 15, 7, 15, 6, 14, 15, // `abcdefghijklmno
    15, 15, 4, 9, 3, 15, 13, 8, 15, 5, 15, 15, 15, 15, 15, 15, // pqrstuvwxyz{|}~DEL
];

pub static IUPACNA_TO_NCBI4NA: [u8; 128] = [
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 14, 2, 13, 0, 0, 4, 11, 0, 0, 12, 0, 3, 15, 0, // @ABCDEFGHIJKLMNO
    0, 0, 5, 6, 8, 0, 7, 9, 0, 10, 0, 0, 0, 0, 0, 0, // PQRSTUVWXYZ[\]^_
    0, 1, 14, 2, 13, 0, 0, 4, 11, 0, 0, 12, 0, 3, 15, 0, // `abcdefghijklmno
    0, 0, 5, 6, 8, 0, 7, 9, 0, 10, 0, 0, 0, 0, 0, 0, // pqrstuvwxyz{|}~DEL
];

pub static AMINOACID_TO_NCBISTDAA: [u8; 128] = [
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 27, 10, 11, 12, 13, 26, // @ABCDEFGHIJKLMNO
    14, 15, 16, 17, 18, 24, 19, 20, 21, 22, 23, 0, 0, 0, 0, 0, // PQRSTUVWXYZ[\]^_
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 27, 10, 11, 12, 13, 26, // `abcdefghijklmno
    14, 15, 16, 17, 18, 24, 19, 20, 21, 22, 23, 0, 0, 0, 0, 0, // pqrstuvwxyz{|}~DEL
];

pub static NCBISTDAA_TO_AMINOACID: [i8; BLASTAA_SIZE] = [
    b'-' as i8, b'A' as i8, b'B' as i8, b'C' as i8, b'D' as i8, b'E' as i8, b'F' as i8, b'G' as i8,
    b'H' as i8, b'I' as i8, b'K' as i8, b'L' as i8, b'M' as i8, b'N' as i8, b'P' as i8, b'Q' as i8,
    b'R' as i8, b'S' as i8, b'T' as i8, b'V' as i8, b'W' as i8, b'X' as i8, b'Y' as i8, b'Z' as i8,
    b'U' as i8, b'*' as i8, b'O' as i8, b'J' as i8,
];

pub static PROT_SENTINEL: u8 = 0;

pub static NUCL_SENTINEL: u8 = 0xF;

#[cfg(test)]
mod tests {
    use super::*;

    // -- test_iupacna_to_blastna_all_bases --
    #[test]
    fn test_iupacna_to_blastna_all_bases() {
        // Standard bases uppercase
        assert_eq!(IUPACNA_TO_BLASTNA[b'A' as usize], 0);
        assert_eq!(IUPACNA_TO_BLASTNA[b'C' as usize], 1);
        assert_eq!(IUPACNA_TO_BLASTNA[b'G' as usize], 2);
        assert_eq!(IUPACNA_TO_BLASTNA[b'T' as usize], 3);
        // Standard bases lowercase
        assert_eq!(IUPACNA_TO_BLASTNA[b'a' as usize], 0);
        assert_eq!(IUPACNA_TO_BLASTNA[b'c' as usize], 1);
        assert_eq!(IUPACNA_TO_BLASTNA[b'g' as usize], 2);
        assert_eq!(IUPACNA_TO_BLASTNA[b't' as usize], 3);
    }

    // -- test_blastna_to_iupacna_roundtrip --
    #[test]
    fn test_blastna_to_iupacna_roundtrip() {
        // For each standard base: IUPAC -> BLASTNA -> IUPAC should be identity
        for &ch in &[b'A', b'C', b'G', b'T'] {
            let blastna = IUPACNA_TO_BLASTNA[ch as usize];
            let back = BLASTNA_TO_IUPACNA[blastna as usize] as u8;
            assert_eq!(back, ch, "Roundtrip failed for {}", ch as char);
        }
    }

    // -- test_aminoacid_to_ncbistdaa --
    #[test]
    fn test_aminoacid_to_ncbistdaa() {
        // Verify standard amino acids
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'A' as usize], 1);
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'C' as usize], 3);
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'D' as usize], 4);
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'E' as usize], 5);
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'F' as usize], 6);
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'G' as usize], 7);
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'H' as usize], 8);
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'I' as usize], 9);
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'K' as usize], 10);
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'L' as usize], 11);
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'M' as usize], 12);
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'N' as usize], 13);
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'P' as usize], 14);
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'Q' as usize], 15);
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'R' as usize], 16);
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'S' as usize], 17);
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'T' as usize], 18);
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'V' as usize], 19);
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'W' as usize], 20);
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'Y' as usize], 22);
        // Stop codon
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'*' as usize], 25);
        // Lowercase matches uppercase
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'a' as usize], 1);
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'm' as usize], 12);
    }

    // -- test_ncbistdaa_to_aminoacid_roundtrip --
    #[test]
    fn test_ncbistdaa_to_aminoacid_roundtrip() {
        // For the 20 standard amino acids, encoding then decoding should round-trip
        for &ch in b"ACDEFGHIKLMNPQRSTVWY" {
            let stdaa = AMINOACID_TO_NCBISTDAA[ch as usize];
            let back = NCBISTDAA_TO_AMINOACID[stdaa as usize] as u8;
            assert_eq!(back, ch, "Roundtrip failed for amino acid {}", ch as char);
        }
    }

    // -- test_ambiguity_codes --
    #[test]
    fn test_ambiguity_codes() {
        // N = any base
        assert_eq!(IUPACNA_TO_BLASTNA[b'N' as usize], 14);
        assert_eq!(IUPACNA_TO_BLASTNA[b'n' as usize], 14);
        // R = purine (A or G)
        assert_eq!(IUPACNA_TO_BLASTNA[b'R' as usize], 4);
        assert_eq!(IUPACNA_TO_BLASTNA[b'r' as usize], 4);
        // Y = pyrimidine (C or T)
        assert_eq!(IUPACNA_TO_BLASTNA[b'Y' as usize], 5);
        assert_eq!(IUPACNA_TO_BLASTNA[b'y' as usize], 5);
        // M = amino (A or C)
        assert_eq!(IUPACNA_TO_BLASTNA[b'M' as usize], 6);
        // K = keto (G or T)
        assert_eq!(IUPACNA_TO_BLASTNA[b'K' as usize], 7);
        // S = strong (G or C)
        assert_eq!(IUPACNA_TO_BLASTNA[b'S' as usize], 9);
        // W = weak (A or T)
        assert_eq!(IUPACNA_TO_BLASTNA[b'W' as usize], 8);
        // B = not A (C, G, T)
        assert_eq!(IUPACNA_TO_BLASTNA[b'B' as usize], 10);
        // D = not C (A, G, T)
        assert_eq!(IUPACNA_TO_BLASTNA[b'D' as usize], 11);
        // H = not G (A, C, T)
        assert_eq!(IUPACNA_TO_BLASTNA[b'H' as usize], 12);
        // V = not T (A, C, G)
        assert_eq!(IUPACNA_TO_BLASTNA[b'V' as usize], 13);
    }

    #[test]
    fn test_ncbi4na_blastna_roundtrip() {
        // Converting NCBI4NA -> BLASTNA -> NCBI4NA should be identity for all 16 values
        for i in 0u8..16 {
            let blastna = NCBI4NA_TO_BLASTNA[i as usize];
            let back = BLASTNA_TO_NCBI4NA[blastna as usize];
            assert_eq!(back, i, "NCBI4NA roundtrip failed for {}", i);
        }
    }

    #[test]
    fn test_blastna_ncbi4na_roundtrip() {
        // Converting BLASTNA -> NCBI4NA -> BLASTNA should be identity for all 16 values
        for i in 0u8..16 {
            let ncbi4na = BLASTNA_TO_NCBI4NA[i as usize];
            let back = NCBI4NA_TO_BLASTNA[ncbi4na as usize];
            assert_eq!(back, i, "BLASTNA roundtrip failed for {}", i);
        }
    }

    #[test]
    fn test_sentinel_values() {
        assert_eq!(PROT_SENTINEL, 0);
        assert_eq!(NUCL_SENTINEL, 0xF);
    }

    #[test]
    fn test_ncbi4na_to_iupacna_standard_bases() {
        // In NCBI4NA: A=1, C=2, G=4, T=8
        assert_eq!(NCBI4NA_TO_IUPACNA[1] as u8, b'A');
        assert_eq!(NCBI4NA_TO_IUPACNA[2] as u8, b'C');
        assert_eq!(NCBI4NA_TO_IUPACNA[4] as u8, b'G');
        assert_eq!(NCBI4NA_TO_IUPACNA[8] as u8, b'T');
        // Gap = 0
        assert_eq!(NCBI4NA_TO_IUPACNA[0] as u8, b'-');
        // N = 15
        assert_eq!(NCBI4NA_TO_IUPACNA[15] as u8, b'N');
    }

    #[test]
    fn test_special_amino_acids() {
        // U = Selenocysteine
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'U' as usize], 24);
        assert_eq!(NCBISTDAA_TO_AMINOACID[24] as u8, b'U');
        // O = Pyrrolysine
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'O' as usize], 26);
        assert_eq!(NCBISTDAA_TO_AMINOACID[26] as u8, b'O');
        // J = Leucine or Isoleucine
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'J' as usize], 27);
        assert_eq!(NCBISTDAA_TO_AMINOACID[27] as u8, b'J');
        // X = unknown
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'X' as usize], 21);
        assert_eq!(NCBISTDAA_TO_AMINOACID[21] as u8, b'X');
        // B = Asp or Asn
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'B' as usize], 2);
        assert_eq!(NCBISTDAA_TO_AMINOACID[2] as u8, b'B');
        // Z = Glu or Gln
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'Z' as usize], 23);
        assert_eq!(NCBISTDAA_TO_AMINOACID[23] as u8, b'Z');
    }

    #[test]
    fn test_iupacna_to_ncbi4na_standard_bases() {
        assert_eq!(IUPACNA_TO_NCBI4NA[b'A' as usize], 1);
        assert_eq!(IUPACNA_TO_NCBI4NA[b'C' as usize], 2);
        assert_eq!(IUPACNA_TO_NCBI4NA[b'G' as usize], 4);
        assert_eq!(IUPACNA_TO_NCBI4NA[b'T' as usize], 8);
        // Lowercase
        assert_eq!(IUPACNA_TO_NCBI4NA[b'a' as usize], 1);
        assert_eq!(IUPACNA_TO_NCBI4NA[b'c' as usize], 2);
        assert_eq!(IUPACNA_TO_NCBI4NA[b'g' as usize], 4);
        assert_eq!(IUPACNA_TO_NCBI4NA[b't' as usize], 8);
        // N = 15 (all bits set)
        assert_eq!(IUPACNA_TO_NCBI4NA[b'N' as usize], 15);
    }
}
