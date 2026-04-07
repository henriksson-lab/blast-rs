//! Rust equivalent of blast_program.c — program type classification.

use blast_core_sys as ffi;

pub type ProgramType = ffi::EBlastProgramType;

// Bitmask constants from blast_program.h
const PROTEIN_QUERY_MASK: u32 = 1 << 0;
const PROTEIN_SUBJECT_MASK: u32 = 1 << 1;
const NUCLEOTIDE_QUERY_MASK: u32 = 1 << 2;
const NUCLEOTIDE_SUBJECT_MASK: u32 = 1 << 3;
const TRANSLATED_QUERY_MASK: u32 = 1 << 4;
const TRANSLATED_SUBJECT_MASK: u32 = 1 << 5;
const PSSM_QUERY_MASK: u32 = 1 << 6;
const PSSM_SUBJECT_MASK: u32 = 1 << 7;
const PATTERN_QUERY_MASK: u32 = 1 << 8;
const MAPPING_MASK: u32 = 1 << 9;

fn has_flag(p: ProgramType, mask: u32) -> bool {
    (p & mask) != 0
}

pub fn query_is_protein(p: ProgramType) -> bool { has_flag(p, PROTEIN_QUERY_MASK) }
pub fn query_is_nucleotide(p: ProgramType) -> bool { has_flag(p, NUCLEOTIDE_QUERY_MASK) }
pub fn query_is_pssm(p: ProgramType) -> bool { has_flag(p, PSSM_QUERY_MASK) }
pub fn subject_is_protein(p: ProgramType) -> bool { has_flag(p, PROTEIN_SUBJECT_MASK) }
pub fn subject_is_nucleotide(p: ProgramType) -> bool { has_flag(p, NUCLEOTIDE_SUBJECT_MASK) }
pub fn subject_is_translated(p: ProgramType) -> bool { has_flag(p, TRANSLATED_SUBJECT_MASK) }
pub fn query_is_translated(p: ProgramType) -> bool { has_flag(p, TRANSLATED_QUERY_MASK) }
pub fn query_is_pattern(p: ProgramType) -> bool { has_flag(p, PATTERN_QUERY_MASK) }

pub fn is_psi_blast(p: ProgramType) -> bool { has_flag(p, PSSM_QUERY_MASK) }
pub fn is_phi_blast(p: ProgramType) -> bool { has_flag(p, PATTERN_QUERY_MASK) }
pub fn is_rps_blast(p: ProgramType) -> bool { has_flag(p, PSSM_SUBJECT_MASK) }
pub fn is_mapping(p: ProgramType) -> bool { has_flag(p, MAPPING_MASK) }

pub fn is_nucleotide(p: ProgramType) -> bool {
    query_is_nucleotide(p) && subject_is_nucleotide(p)
        && !query_is_translated(p) && !subject_is_translated(p)
}

pub fn is_valid(p: ProgramType) -> bool {
    matches!(p,
        ffi::EBlastProgramType_eBlastTypeBlastp
        | ffi::EBlastProgramType_eBlastTypeBlastn
        | ffi::EBlastProgramType_eBlastTypeBlastx
        | ffi::EBlastProgramType_eBlastTypeTblastn
        | ffi::EBlastProgramType_eBlastTypeTblastx
        | ffi::EBlastProgramType_eBlastTypePsiBlast
        | ffi::EBlastProgramType_eBlastTypePsiTblastn
        | ffi::EBlastProgramType_eBlastTypeRpsBlast
        | ffi::EBlastProgramType_eBlastTypeRpsTblastn
        | ffi::EBlastProgramType_eBlastTypePhiBlastp
        | ffi::EBlastProgramType_eBlastTypePhiBlastn
        | ffi::EBlastProgramType_eBlastTypeMapping
    )
}

/// Number of contexts for a given program type.
pub fn num_contexts(p: ProgramType) -> u32 {
    if is_nucleotide(p) || is_mapping(p) {
        2 // plus and minus strand
    } else if query_is_translated(p) {
        6 // 6 reading frames
    } else {
        1 // protein
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_blastn() {
        let p = ffi::EBlastProgramType_eBlastTypeBlastn;
        assert!(is_nucleotide(p));
        assert!(query_is_nucleotide(p));
        assert!(subject_is_nucleotide(p));
        assert!(!query_is_protein(p));
        assert!(!subject_is_translated(p));
        assert!(is_valid(p));
        assert_eq!(num_contexts(p), 2);
    }

    #[test]
    fn test_blastp() {
        let p = ffi::EBlastProgramType_eBlastTypeBlastp;
        assert!(!is_nucleotide(p));
        assert!(query_is_protein(p));
        assert!(subject_is_protein(p));
        assert!(is_valid(p));
        assert_eq!(num_contexts(p), 1);
    }

    #[test]
    fn test_blastx() {
        let p = ffi::EBlastProgramType_eBlastTypeBlastx;
        assert!(query_is_nucleotide(p));
        assert!(query_is_translated(p));
        assert!(subject_is_protein(p));
        assert!(is_valid(p));
        assert_eq!(num_contexts(p), 6);
    }
}
