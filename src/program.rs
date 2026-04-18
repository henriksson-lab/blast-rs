//! Program type classification (pure Rust, no FFI).

pub type ProgramType = u32;

// Program type constants matching NCBI BLAST+ EBlastProgramType
pub const BLASTP: ProgramType = 0x01 | 0x02; // protein query + protein subject
pub const BLASTN: ProgramType = 0x04 | 0x08; // nucleotide query + nucleotide subject
pub const BLASTX: ProgramType = 0x04 | 0x02 | 0x10; // nucleotide query (translated) + protein subject
pub const TBLASTN: ProgramType = 0x01 | 0x08 | 0x20; // protein query + nucleotide subject (translated)
pub const TBLASTX: ProgramType = 0x04 | 0x08 | 0x10 | 0x20; // both translated
pub const PSI_BLAST: ProgramType = 0x01 | 0x02 | 0x40; // PSSM query
pub const RPS_BLAST: ProgramType = 0x01 | 0x80; // PSSM subject
pub const MAPPING: ProgramType = 0x04 | 0x08 | 0x200; // mapping

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

pub fn query_is_protein(p: ProgramType) -> bool {
    has_flag(p, PROTEIN_QUERY_MASK)
}
pub fn query_is_nucleotide(p: ProgramType) -> bool {
    has_flag(p, NUCLEOTIDE_QUERY_MASK)
}
pub fn query_is_pssm(p: ProgramType) -> bool {
    has_flag(p, PSSM_QUERY_MASK)
}
pub fn subject_is_protein(p: ProgramType) -> bool {
    has_flag(p, PROTEIN_SUBJECT_MASK)
}
pub fn subject_is_nucleotide(p: ProgramType) -> bool {
    has_flag(p, NUCLEOTIDE_SUBJECT_MASK)
}
pub fn subject_is_translated(p: ProgramType) -> bool {
    has_flag(p, TRANSLATED_SUBJECT_MASK)
}
pub fn query_is_translated(p: ProgramType) -> bool {
    has_flag(p, TRANSLATED_QUERY_MASK)
}
pub fn query_is_pattern(p: ProgramType) -> bool {
    has_flag(p, PATTERN_QUERY_MASK)
}

pub fn is_psi_blast(p: ProgramType) -> bool {
    has_flag(p, PSSM_QUERY_MASK)
}
pub fn is_phi_blast(p: ProgramType) -> bool {
    has_flag(p, PATTERN_QUERY_MASK)
}
pub fn is_rps_blast(p: ProgramType) -> bool {
    has_flag(p, PSSM_SUBJECT_MASK)
}
pub fn is_mapping(p: ProgramType) -> bool {
    has_flag(p, MAPPING_MASK)
}

pub fn is_nucleotide(p: ProgramType) -> bool {
    query_is_nucleotide(p)
        && subject_is_nucleotide(p)
        && !query_is_translated(p)
        && !subject_is_translated(p)
}

pub fn is_valid(p: ProgramType) -> bool {
    matches!(
        p,
        BLASTP | BLASTN | BLASTX | TBLASTN | TBLASTX | PSI_BLAST | RPS_BLAST | MAPPING
    )
}

/// Port of NCBI `BLAST_GetNumberOfContexts` (`blast_util.c:1373`).
/// Returns `NUM_FRAMES` (6) for translated queries, `NUM_STRANDS` (2)
/// for nucleotide queries, `1` for protein queries on valid programs,
/// and `0` for invalid programs — matching NCBI's dispatch order.
pub fn num_contexts(p: ProgramType) -> u32 {
    if query_is_translated(p) {
        6
    } else if query_is_nucleotide(p) {
        2
    } else if is_valid(p) {
        1
    } else {
        0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_blastn() {
        assert!(is_nucleotide(BLASTN));
        assert!(query_is_nucleotide(BLASTN));
        assert!(subject_is_nucleotide(BLASTN));
        assert!(!query_is_protein(BLASTN));
        assert!(is_valid(BLASTN));
        assert_eq!(num_contexts(BLASTN), 2);
    }

    #[test]
    fn test_blastp() {
        assert!(!is_nucleotide(BLASTP));
        assert!(query_is_protein(BLASTP));
        assert!(subject_is_protein(BLASTP));
        assert!(is_valid(BLASTP));
        assert_eq!(num_contexts(BLASTP), 1);
    }

    #[test]
    fn test_blastx() {
        assert!(query_is_nucleotide(BLASTX));
        assert!(query_is_translated(BLASTX));
        assert!(subject_is_protein(BLASTX));
        assert!(is_valid(BLASTX));
        assert_eq!(num_contexts(BLASTX), 6);
    }
}
