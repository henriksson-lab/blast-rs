//! Sequence encoding conversions matching the C BLAST core.

/// IUPAC nucleotide character to BLASTNA encoding.
/// BLASTNA is a permutation of NCBI4na where the first 4 values (0-3)
/// match NCBI2na: A=0, C=1, G=2, T=3.
/// Ambiguity codes use values 4-15, with 15 = N (any base).
pub fn iupacna_to_blastna(c: u8) -> u8 {
    match c {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        b'R' | b'r' => 4, // A or G
        b'Y' | b'y' => 5, // C or T
        b'M' | b'm' => 6, // A or C
        b'K' | b'k' => 7, // G or T
        b'W' | b'w' => 8, // A or T
        b'S' | b's' => 9, // C or G
        b'B' | b'b' => 10, // C, G, or T
        b'D' | b'd' => 11, // A, G, or T
        b'H' | b'h' => 12, // A, C, or T
        b'V' | b'v' => 13, // A, C, or G
        b'N' | b'n' => 14, // any base
        _ => 15,           // gap or unknown
    }
}

/// IUPAC amino acid character to NCBIstdaa encoding.
pub fn aminoacid_to_ncbistdaa(c: u8) -> u8 {
    match c {
        b'-' => 0,  // Gap
        b'A' | b'a' => 1,
        b'B' | b'b' => 2,  // Asp or Asn
        b'C' | b'c' => 3,
        b'D' | b'd' => 4,
        b'E' | b'e' => 5,
        b'F' | b'f' => 6,
        b'G' | b'g' => 7,
        b'H' | b'h' => 8,
        b'I' | b'i' => 9,
        b'K' | b'k' => 10,
        b'L' | b'l' => 11,
        b'M' | b'm' => 12,
        b'N' | b'n' => 13,
        b'P' | b'p' => 14,
        b'Q' | b'q' => 15,
        b'R' | b'r' => 16,
        b'S' | b's' => 17,
        b'T' | b't' => 18,
        b'V' | b'v' => 19,
        b'W' | b'w' => 20,
        b'X' | b'x' => 21, // Unknown
        b'Y' | b'y' => 22,
        b'Z' | b'z' => 23, // Glu or Gln
        b'U' | b'u' => 24, // Selenocysteine
        b'*' => 25,         // Stop
        b'O' | b'o' => 26,  // Pyrrolysine
        b'J' | b'j' => 27,  // Leu or Ile
        _ => 21,             // Unknown -> X
    }
}

/// BLASTNA sentinel byte value (used to mark query boundaries).
pub const BLASTNA_SENTINEL: u8 = 15;
