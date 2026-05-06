//! Sequence encoding conversions matching the C BLAST core.

/// IUPAC nucleotide character to BLASTNA encoding.
/// BLASTNA is a permutation of NCBI4na where the first 4 values (0-3)
/// match NCBI2na: A=0, C=1, G=2, T=3.
/// Ambiguity codes use values 4-15, with 15 = N (any base).
pub fn iupacna_to_blastna(c: u8) -> u8 {
    crate::encoding::iupacna_to_blastna_base(c)
}

/// IUPAC amino acid character to NCBIstdaa encoding.
pub fn aminoacid_to_ncbistdaa(c: u8) -> u8 {
    crate::encoding::aminoacid_to_ncbistdaa_base(c)
}

/// BLASTNA sentinel byte value (used to mark query boundaries).
/// Alias for `crate::encoding::NUCL_SENTINEL` (both = 0xF = 15).
#[allow(dead_code)]
pub const BLASTNA_SENTINEL: u8 = crate::encoding::NUCL_SENTINEL;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn compatibility_encoders_delegate_to_canonical_tables() {
        for b in 0u8..=127 {
            assert_eq!(
                iupacna_to_blastna(b),
                crate::encoding::iupacna_to_blastna_base(b),
                "IUPACNA byte {b}"
            );
            assert_eq!(
                aminoacid_to_ncbistdaa(b),
                crate::encoding::aminoacid_to_ncbistdaa_base(b),
                "amino-acid byte {b}"
            );
        }
        assert_eq!(iupacna_to_blastna(b'U'), 3);
        assert_eq!(iupacna_to_blastna(b'u'), 3);
    }
}
