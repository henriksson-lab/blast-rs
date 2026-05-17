//! Program type classification (pure Rust, no FFI).

pub type ProgramType = u32;

// Program type constants matching NCBI BLAST+ EBlastProgramType
pub const BLASTP: ProgramType = 0x01 | 0x02; // protein query + protein subject
pub const BLASTN: ProgramType = 0x04 | 0x08; // nucleotide query + nucleotide subject
pub const BLASTX: ProgramType = 0x04 | 0x02 | 0x10; // nucleotide query (translated) + protein subject
pub const TBLASTN: ProgramType = 0x01 | 0x08 | 0x20; // protein query + nucleotide subject (translated)
pub const TBLASTX: ProgramType = 0x04 | 0x08 | 0x10 | 0x20; // both translated
pub const PSI_BLAST: ProgramType = 0x01 | 0x02 | 0x40; // PSSM query
pub const PSI_TBLASTN: ProgramType = 0x01 | 0x08 | 0x20 | 0x40; // PSSM query + translated subject
pub const RPS_BLAST: ProgramType = 0x01 | 0x80; // PSSM subject
pub const RPS_TBLASTN: ProgramType = 0x04 | 0x02 | 0x10 | 0x80; // translated query + PSSM subject
pub const PHI_BLASTP: ProgramType = 0x01 | 0x02 | 0x100; // pattern query + protein subject
pub const PHI_BLASTN: ProgramType = 0x04 | 0x08 | 0x100; // pattern query + nucleotide subject
pub const MAPPING: ProgramType = 0x04 | 0x08 | 0x200; // mapping
pub const UNDEFINED: ProgramType = 0x0;

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

pub fn blast_query_is_protein(p: ProgramType) -> bool {
    has_flag(p, PROTEIN_QUERY_MASK)
}
pub fn blast_query_is_nucleotide(p: ProgramType) -> bool {
    has_flag(p, NUCLEOTIDE_QUERY_MASK)
}
pub fn blast_query_is_pssm(p: ProgramType) -> bool {
    has_flag(p, PSSM_QUERY_MASK)
}
pub fn blast_subject_is_pssm(p: ProgramType) -> bool {
    has_flag(p, PSSM_SUBJECT_MASK)
}
pub fn blast_subject_is_protein(p: ProgramType) -> bool {
    has_flag(p, PROTEIN_SUBJECT_MASK)
}
pub fn blast_subject_is_nucleotide(p: ProgramType) -> bool {
    has_flag(p, NUCLEOTIDE_SUBJECT_MASK)
}
pub fn blast_subject_is_translated(p: ProgramType) -> bool {
    has_flag(p, TRANSLATED_SUBJECT_MASK)
}
pub fn blast_query_is_translated(p: ProgramType) -> bool {
    has_flag(p, TRANSLATED_QUERY_MASK)
}
pub fn blast_query_is_pattern(p: ProgramType) -> bool {
    has_flag(p, PATTERN_QUERY_MASK)
}

pub fn blast_program_is_psi_blast(p: ProgramType) -> bool {
    has_flag(p, PSSM_QUERY_MASK)
}
pub fn blast_program_is_phi_blast(p: ProgramType) -> bool {
    has_flag(p, PATTERN_QUERY_MASK)
}
pub fn blast_program_is_rps_blast(p: ProgramType) -> bool {
    has_flag(p, PSSM_SUBJECT_MASK)
}
pub fn blast_program_is_mapping(p: ProgramType) -> bool {
    has_flag(p, MAPPING_MASK)
}

pub fn blast_program_is_nucleotide(p: ProgramType) -> bool {
    blast_query_is_nucleotide(p)
        && blast_subject_is_nucleotide(p)
        && !blast_query_is_translated(p)
        && !blast_subject_is_translated(p)
}

pub fn blast_program_is_valid(p: ProgramType) -> bool {
    match p {
        BLASTP => true,
        BLASTN => true,
        BLASTX => true,
        TBLASTN => true,
        TBLASTX => true,
        PSI_BLAST => true,
        PSI_TBLASTN => true,
        RPS_BLAST => true,
        RPS_TBLASTN => true,
        PHI_BLASTP => true,
        PHI_BLASTN => true,
        MAPPING => true,
        _ => false,
    }
}

/// Port of NCBI `BlastProgram2Number` (`blast_util.c:278`).
/// naming: Rust preserves the C `2` spelling used in the source symbol.
/// Returns `(status, program)`, where status is `1` only for a NULL input in C.
/// Unknown non-null names return success with `eBlastTypeUndefined`.
pub fn blast_program_2_number(program: Option<&str>) -> (i16, ProgramType) {
    let Some(program) = program else {
        return (1, UNDEFINED);
    };

    let number = if program.eq_ignore_ascii_case("blastn") {
        BLASTN
    } else if program.eq_ignore_ascii_case("blastp") {
        BLASTP
    } else if program.eq_ignore_ascii_case("blastx") {
        BLASTX
    } else if program.eq_ignore_ascii_case("tblastn") {
        TBLASTN
    } else if program.eq_ignore_ascii_case("tblastx") {
        TBLASTX
    } else if program.eq_ignore_ascii_case("rpsblast") {
        RPS_BLAST
    } else if program.eq_ignore_ascii_case("rpstblastn") {
        RPS_TBLASTN
    } else if program.eq_ignore_ascii_case("psiblast") {
        PSI_BLAST
    } else if program.eq_ignore_ascii_case("psitblastn") {
        PSI_TBLASTN
    } else if program.eq_ignore_ascii_case("phiblastn") {
        PHI_BLASTN
    } else if program.eq_ignore_ascii_case("phiblastp") {
        PHI_BLASTP
    } else if program.eq_ignore_ascii_case("mapper") {
        MAPPING
    } else {
        UNDEFINED
    };

    (0, number)
}

/// Port of NCBI `BlastNumber2Program` (`blast_util.c:312`).
/// naming: Rust preserves the C `2` spelling used in the source symbol.
/// Rust returns an owned canonical program name instead of allocating through
/// the caller's `char**`; `None` models the C NULL output pointer error.
pub fn blast_number_2_program(number: ProgramType, program: Option<()>) -> (i16, Option<String>) {
    if program.is_none() {
        return (1, None);
    }

    let name = match number {
        BLASTN => "blastn",
        BLASTP => "blastp",
        BLASTX => "blastx",
        TBLASTN => "tblastn",
        TBLASTX => "tblastx",
        RPS_BLAST => "rpsblast",
        RPS_TBLASTN => "rpstblastn",
        PSI_BLAST => "psiblast",
        PSI_TBLASTN => "psitblastn",
        PHI_BLASTP => "phiblastp",
        PHI_BLASTN => "phiblastn",
        MAPPING => "mapper",
        _ => "unknown",
    };

    (0, Some(name.to_owned()))
}

/// Port of NCBI `BLAST_GetNumberOfContexts` (`blast_util.c:1373`).
/// Returns `NUM_FRAMES` (6) for translated queries, `NUM_STRANDS` (2)
/// for nucleotide queries, `1` for protein queries on valid programs,
/// and `0` for invalid programs — matching NCBI's dispatch order.
pub fn num_contexts(p: ProgramType) -> u32 {
    if blast_query_is_translated(p) {
        6
    } else if blast_query_is_nucleotide(p) {
        2
    } else if blast_program_is_valid(p) {
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
        assert!(blast_program_is_nucleotide(BLASTN));
        assert!(blast_query_is_nucleotide(BLASTN));
        assert!(blast_subject_is_nucleotide(BLASTN));
        assert!(!blast_query_is_protein(BLASTN));
        assert!(blast_program_is_valid(BLASTN));
        assert_eq!(num_contexts(BLASTN), 2);
    }

    #[test]
    fn test_blastp() {
        assert!(!blast_program_is_nucleotide(BLASTP));
        assert!(blast_query_is_protein(BLASTP));
        assert!(blast_subject_is_protein(BLASTP));
        assert!(blast_program_is_valid(BLASTP));
        assert_eq!(num_contexts(BLASTP), 1);
    }

    #[test]
    fn test_blastx() {
        assert!(blast_query_is_nucleotide(BLASTX));
        assert!(blast_query_is_translated(BLASTX));
        assert!(blast_subject_is_protein(BLASTX));
        assert!(blast_program_is_valid(BLASTX));
        assert_eq!(num_contexts(BLASTX), 6);
    }

    #[test]
    fn program_is_valid_matches_ncbi_extended_programs() {
        assert!(blast_program_is_valid(PSI_TBLASTN));
        assert!(blast_program_is_valid(RPS_TBLASTN));
        assert!(blast_program_is_valid(PHI_BLASTP));
        assert!(blast_program_is_valid(PHI_BLASTN));
        assert!(!blast_program_is_valid(0));
    }

    #[test]
    fn blast_program_2_number_matches_ncbi_names() {
        let cases = [
            ("blastn", BLASTN),
            ("blastp", BLASTP),
            ("blastx", BLASTX),
            ("tblastn", TBLASTN),
            ("tblastx", TBLASTX),
            ("rpsblast", RPS_BLAST),
            ("rpstblastn", RPS_TBLASTN),
            ("psiblast", PSI_BLAST),
            ("psitblastn", PSI_TBLASTN),
            ("phiblastn", PHI_BLASTN),
            ("phiblastp", PHI_BLASTP),
            ("mapper", MAPPING),
        ];

        for (name, program) in cases {
            assert_eq!(blast_program_2_number(Some(name)), (0, program));
            assert_eq!(
                blast_program_2_number(Some(&name.to_ascii_uppercase())),
                (0, program)
            );
        }
        assert_eq!(blast_program_2_number(Some("junk")), (0, UNDEFINED));
        assert_eq!(blast_program_2_number(None), (1, UNDEFINED));
    }

    #[test]
    fn blast_number_2_program_matches_ncbi_names() {
        let cases = [
            (BLASTN, "blastn"),
            (BLASTP, "blastp"),
            (BLASTX, "blastx"),
            (TBLASTN, "tblastn"),
            (TBLASTX, "tblastx"),
            (RPS_BLAST, "rpsblast"),
            (RPS_TBLASTN, "rpstblastn"),
            (PSI_BLAST, "psiblast"),
            (PSI_TBLASTN, "psitblastn"),
            (PHI_BLASTP, "phiblastp"),
            (PHI_BLASTN, "phiblastn"),
            (MAPPING, "mapper"),
        ];

        for (program, name) in cases {
            assert_eq!(
                blast_number_2_program(program, Some(())),
                (0, Some(name.to_owned()))
            );
        }
        assert_eq!(
            blast_number_2_program(UNDEFINED, Some(())),
            (0, Some("unknown".to_owned()))
        );
        assert_eq!(blast_number_2_program(BLASTP, None), (1, None));
    }
}
