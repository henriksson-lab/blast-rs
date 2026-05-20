//! Rust equivalent of blast_encoding.c
//! These tables will replace the C versions when the FFI layer is removed.

/// NCBI `BLAST2NA_SIZE` (`blast_encoding.h:91`): compressed 2-bit
/// nucleotide alphabet — 4 symbols (A, C, G, T).
pub const BLAST2NA_SIZE: usize = 4;
/// NCBI `BLASTNA_SIZE` (`blast_encoding.h:92`): full 16-symbol nucleotide
/// alphabet (A, C, G, T, ambiguity codes, gap).
pub const BLASTNA_SIZE: usize = 16;
/// NCBI `BLASTAA_SIZE` (`blast_encoding.h:93`): 28-symbol amino-acid
/// alphabet (NCBIstdaa).
pub const BLASTAA_SIZE: usize = 28;

/// NCBI `BLASTNA_SEQ_CODE` (`blast_encoding.h:96`): sequence-code
/// identifier for BLASTNA (BLAST extension, not a registered ASN.1 code).
pub const BLASTNA_SEQ_CODE: i32 = 99;
/// NCBI `BLASTAA_SEQ_CODE` (`blast_encoding.h:98`): NCBIstdaa = ASN.1 code 11.
pub const BLASTAA_SEQ_CODE: i32 = 11;
/// NCBI `NCBI4NA_SEQ_CODE` (`blast_encoding.h:99`): NCBI4na = ASN.1 code 4.
pub const NCBI4NA_SEQ_CODE: i32 = 4;

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

/// NCBIstdaa gap/sentinel character (`eGapChar`).
pub const NCBISTDAA_GAP: u8 = 0;
/// NCBIstdaa alanine character (`A`).
pub const NCBISTDAA_A: u8 = 1;
/// NCBIstdaa aspartate/asparagine ambiguity character (`B`).
pub const NCBISTDAA_B: u8 = 2;
/// NCBIstdaa cysteine character (`C`).
pub const NCBISTDAA_C: u8 = 3;
/// NCBIstdaa aspartate character (`D`).
pub const NCBISTDAA_D: u8 = 4;
/// NCBIstdaa glutamate character (`E`).
pub const NCBISTDAA_E: u8 = 5;
/// NCBIstdaa isoleucine character (`I`).
pub const NCBISTDAA_I: u8 = 9;
/// NCBIstdaa leucine character (`L`).
pub const NCBISTDAA_L: u8 = 11;
/// NCBIstdaa asparagine character (`N`).
pub const NCBISTDAA_N: u8 = 13;
/// NCBIstdaa glutamine character (`Q`).
pub const NCBISTDAA_Q: u8 = 15;
/// NCBIstdaa unknown amino-acid character (`X`).
pub const NCBISTDAA_X: u8 = 21;
/// NCBIstdaa glutamate/glutamine ambiguity character (`Z`).
pub const NCBISTDAA_Z: u8 = 23;
/// NCBIstdaa selenocysteine character (`U`).
pub const NCBISTDAA_U: u8 = 24;
/// NCBIstdaa stop character (`*`).
pub const NCBISTDAA_STOP: u8 = 25;
/// NCBIstdaa pyrrolysine character (`O`).
pub const NCBISTDAA_O: u8 = 26;
/// NCBIstdaa isoleucine/leucine ambiguity character (`J`).
pub const NCBISTDAA_J: u8 = 27;
/// The 20 standard amino acids in NCBIstdaa code order.
///
/// NCBIstdaa reserves 21 for `X`, so tyrosine (`Y`) is 22 rather than part of
/// a contiguous 1..=20 standard-residue range.
pub const NCBISTDAA_STANDARD_RESIDUES: [u8; 20] = [
    1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22,
];

/// Bit i is set iff i is one of the 20 standard NCBIstdaa residue codes.
/// Used for the hot-path `is_ncbistdaa_standard_residue` check — replaces
/// `.contains(&b)` on `NCBISTDAA_STANDARD_RESIDUES`, which was emitting
/// memchr calls for every byte in `blast_read_aa_composition` (~50% of tblastx
/// per-subject runtime).
const STANDARD_RESIDUE_MASK: u32 = {
    let mut mask: u32 = 0;
    let mut i = 0;
    while i < NCBISTDAA_STANDARD_RESIDUES.len() {
        mask |= 1 << NCBISTDAA_STANDARD_RESIDUES[i];
        i += 1;
    }
    mask
};

/// Return whether a code is one of the 20 standard NCBIstdaa residues.
#[inline]
pub fn is_ncbistdaa_standard_residue(b: u8) -> bool {
    b < 32 && (STANDARD_RESIDUE_MASK & (1u32 << b)) != 0
}

/// Return whether a code contributes to amino-acid composition.
///
/// NCBI's composition path treats selenocysteine (`U`) as a true residue and
/// folds it into cysteine after counting.
pub fn is_ncbistdaa_composition_residue(b: u8) -> bool {
    is_ncbistdaa_standard_residue(b) || b == NCBISTDAA_U
}

const OLD_AMBIGUITY_MAX_POSITION: u32 = 0x00ff_ffff;
const OLD_AMBIGUITY_MAX_RUN: u32 = 16;
const NEW_AMBIGUITY_MAX_RUN: u32 = 4096;

/// Convert one ASCII amino-acid byte to NCBIstdaa.
pub fn aminoacid_to_ncbistdaa_base(b: u8) -> u8 {
    AMINOACID_TO_NCBISTDAA[b as usize & 0x7f]
}

/// Convert an ASCII amino-acid sequence to NCBIstdaa.
pub fn encode_ncbistdaa_sequence(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .map(|&b| aminoacid_to_ncbistdaa_base(b))
        .collect()
}

/// Convert one NCBIstdaa byte to an ASCII amino-acid byte, using `X` for
/// out-of-range values.
pub fn ncbistdaa_to_aminoacid_base(b: u8) -> u8 {
    if (b as usize) < NCBISTDAA_TO_AMINOACID.len() {
        NCBISTDAA_TO_AMINOACID[b as usize] as u8
    } else {
        b'X'
    }
}

/// Convert one NCBIstdaa byte to an ASCII amino-acid display character.
pub fn ncbistdaa_to_aminoacid_char(b: u8) -> char {
    ncbistdaa_to_aminoacid_base(b) as char
}

/// Convert an NCBIstdaa sequence to ASCII amino-acid bytes.
pub fn ncbistdaa_to_aminoacid_sequence(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .map(|&b| ncbistdaa_to_aminoacid_base(b))
        .collect()
}

/// Convert an NCBIstdaa sequence to an ASCII amino-acid display string.
pub fn ncbistdaa_to_aminoacid_string(seq: &[u8]) -> String {
    // SAFETY: the canonical mapping emits ASCII amino-acid bytes only.
    unsafe { String::from_utf8_unchecked(ncbistdaa_to_aminoacid_sequence(seq)) }
}

/// Convert one ASCII IUPAC nucleotide to NCBI4na, mapping RNA `U` to DNA `T`.
/// Non-IUPAC bytes become `N`, matching BLAST's conservative input handling.
pub fn iupacna_to_ncbi4na_base(b: u8) -> u8 {
    match b {
        b'U' | b'u' => 8,
        _ if (b as usize) < IUPACNA_TO_NCBI4NA.len() => {
            let v = IUPACNA_TO_NCBI4NA[b as usize];
            if v == 0 {
                15
            } else {
                v
            }
        }
        _ => 15,
    }
}

/// Convert an ASCII IUPAC nucleotide sequence to NCBI4na.
///
/// RNA `U`/`u` is mapped to DNA `T`; non-IUPAC bytes become `N` so
/// downstream translation produces NCBIstdaa `X` instead of translating a
/// junk byte.
pub fn encode_ncbi4na_sequence(seq: &[u8]) -> Vec<u8> {
    seq.iter().map(|&b| iupacna_to_ncbi4na_base(b)).collect()
}

/// Convert one ASCII IUPAC nucleotide to BLASTNA, mapping RNA `U` to DNA `T`.
pub fn iupacna_to_blastna_base(b: u8) -> u8 {
    match b {
        b'U' | b'u' => 3,
        _ if (b as usize) < IUPACNA_TO_BLASTNA.len() => IUPACNA_TO_BLASTNA[b as usize],
        _ => 15,
    }
}

/// Convert an ASCII IUPAC nucleotide sequence to BLASTNA.
pub fn encode_blastna_sequence(seq: &[u8]) -> Vec<u8> {
    seq.iter().map(|&b| iupacna_to_blastna_base(b)).collect()
}

pub fn is_unambiguous_ncbi4na(b: u8) -> bool {
    matches!(b, 1 | 2 | 4 | 8)
}

pub fn blastna_to_ncbi4na_base(b: u8) -> u8 {
    BLASTNA_TO_NCBI4NA.get(b as usize).copied().unwrap_or(0)
}

/// Convert a BLASTNA sequence to NCBI4na bytes.
pub fn blastna_to_ncbi4na_sequence(seq: &[u8]) -> Vec<u8> {
    seq.iter().map(|&b| blastna_to_ncbi4na_base(b)).collect()
}

pub fn ncbi4na_to_blastna_base(b: u8) -> u8 {
    NCBI4NA_TO_BLASTNA.get(b as usize).copied().unwrap_or(15)
}

pub fn ncbi4na_to_iupacna_base(b: u8) -> u8 {
    NCBI4NA_TO_IUPACNA
        .get(b as usize)
        .copied()
        .unwrap_or(b'N' as i8) as u8
}

/// Convert one NCBI4na byte to an ASCII IUPAC/gap display character.
pub fn ncbi4na_to_iupacna_char(b: u8) -> char {
    ncbi4na_to_iupacna_base(b) as char
}

/// Convert an NCBI4na sequence to ASCII IUPAC/gap display bytes.
pub fn ncbi4na_to_iupacna_sequence(seq: &[u8]) -> Vec<u8> {
    seq.iter().map(|&b| ncbi4na_to_iupacna_base(b)).collect()
}

/// Convert an NCBI4na sequence to an ASCII IUPAC/gap display string.
pub fn ncbi4na_to_iupacna_string(seq: &[u8]) -> String {
    // SAFETY: the canonical mapping emits ASCII IUPAC/gap bytes only.
    unsafe { String::from_utf8_unchecked(ncbi4na_to_iupacna_sequence(seq)) }
}

pub fn blastna_bases_overlap(a: u8, b: u8) -> bool {
    blastna_to_ncbi4na_base(a) & blastna_to_ncbi4na_base(b) != 0
}

pub fn blastna_degeneracy(b: u8) -> usize {
    (0..4)
        .filter(|&base| blastna_bases_overlap(b, base))
        .count()
}

pub fn blastna_pair_score(a: u8, b: u8, reward: i32, penalty: i32) -> i32 {
    if a < 4 && b < 4 {
        return if a == b { reward } else { penalty };
    }
    if a >= NUCL_SENTINEL || b >= NUCL_SENTINEL {
        return i32::MIN / 2;
    }
    if !blastna_bases_overlap(a, b) {
        return penalty;
    }
    let degen = blastna_degeneracy(a.max(b)) as i32;
    // NCBI `blast_stat.c:1106`:
    //   (Int4)BLAST_Nint((double)((degeneracy[i2]-1)*penalty + reward) / degeneracy[i2])
    crate::math::blast_nint(((degen - 1) * penalty + reward) as f64 / degen as f64) as i32
}

pub fn complement_blastna_base(b: u8) -> u8 {
    const TABLE: [u8; BLASTNA_SIZE] = [3, 2, 1, 0, 5, 4, 7, 6, 8, 9, 13, 12, 11, 10, 14, 15];
    TABLE.get(b as usize).copied().unwrap_or(15)
}

/// Reverse-complement a BLASTNA sequence.
pub fn reverse_complement_blastna_sequence(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| complement_blastna_base(b))
        .collect()
}

pub fn blastna_to_iupacna_base(b: u8) -> u8 {
    BLASTNA_TO_IUPACNA
        .get(b as usize)
        .copied()
        .unwrap_or(b'-' as i8) as u8
}

/// Convert one BLASTNA byte to an ASCII IUPAC/gap display character.
pub fn blastna_to_iupacna_char(b: u8) -> char {
    blastna_to_iupacna_base(b) as char
}

/// Convert a BLASTNA sequence to ASCII IUPAC/gap display bytes.
pub fn blastna_to_iupacna_sequence(seq: &[u8]) -> Vec<u8> {
    seq.iter().map(|&b| blastna_to_iupacna_base(b)).collect()
}

/// Convert a BLASTNA sequence to an ASCII IUPAC/gap display string.
pub fn blastna_to_iupacna_string(seq: &[u8]) -> String {
    // SAFETY: the canonical mapping emits ASCII IUPAC/gap bytes only.
    unsafe { String::from_utf8_unchecked(blastna_to_iupacna_sequence(seq)) }
}

pub fn complement_iupacna_ascii_base(b: u8) -> u8 {
    let blastna = iupacna_to_blastna_base(b);
    if blastna == NUCL_SENTINEL {
        b'N'
    } else {
        blastna_to_iupacna_base(complement_blastna_base(blastna))
    }
}

/// Reverse-complement an ASCII IUPAC nucleotide sequence.
pub fn reverse_complement_iupacna_sequence(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| complement_iupacna_ascii_base(b))
        .collect()
}

pub fn complement_ncbi4na_base(b: u8) -> u8 {
    const TABLE: [u8; BLASTNA_SIZE] = [0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15];
    TABLE[(b as usize) & 0x0f]
}

pub fn ncbi2na_base_at(packed: &[u8], pos: usize) -> u8 {
    let byte = packed[pos >> 2];
    (byte >> (6 - 2 * (pos & 3))) & 3
}

/// Convert either an ASCII IUPAC nucleotide or a BLASTNA code to the 2-bit
/// NCBI2na base used by DUST and packed nucleotide subjects.
pub fn blastna_or_iupac_to_ncbi2na_base(base: u8) -> u8 {
    match base {
        b'C' | b'c' | 1 => 1,
        b'G' | b'g' | 2 => 2,
        b'T' | b't' | b'U' | b'u' | 3 => 3,
        b'N' | b'n' => 0,
        _ if base <= 15 => base & 0x3,
        _ => 0,
    }
}

/// Pack pre-encoded NCBI2na/BLASTNA bases into 2-bit bytes without the BLAST
/// database remainder sentinel byte.
pub fn pack_ncbi2na_bases(bases: &[u8]) -> Vec<u8> {
    let mut packed = vec![0u8; bases.len().div_ceil(4)];
    for (i, &base) in bases.iter().enumerate() {
        packed[i >> 2] |= (base & 0x03) << (6 - 2 * (i & 3));
    }
    packed
}

/// Pack an ASCII IUPAC nucleotide sequence as BLAST DB NCBI2na.
///
/// Ambiguous bases are encoded as A in the packed stream; callers that need
/// exact ambiguity preservation should write [`encode_ncbi2na_ambiguity_data`]
/// after the packed bytes.
pub fn encode_ncbi2na_sequence(seq: &[u8]) -> Vec<u8> {
    let mut packed = Vec::with_capacity(seq.len() / 4 + 1);
    let full_bytes = seq.len() / 4;
    let remainder = seq.len() % 4;
    for i in 0..full_bytes {
        let b = (iupacna_to_ncbi2na_base(seq[i * 4]) << 6)
            | (iupacna_to_ncbi2na_base(seq[i * 4 + 1]) << 4)
            | (iupacna_to_ncbi2na_base(seq[i * 4 + 2]) << 2)
            | iupacna_to_ncbi2na_base(seq[i * 4 + 3]);
        packed.push(b);
    }
    if remainder > 0 {
        let mut last = 0u8;
        for j in 0..remainder {
            last |= iupacna_to_ncbi2na_base(seq[full_bytes * 4 + j]) << (6 - 2 * j);
        }
        last |= remainder as u8;
        packed.push(last);
    } else {
        packed.push(0);
    }
    packed
}

fn iupacna_to_ncbi2na_base(b: u8) -> u8 {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' | b'U' | b'u' => 3,
        _ => 0,
    }
}

/// Encode the BLAST v4 nucleotide DB ambiguity sidecar.
///
/// BLAST stores subject sequence bases packed as 2-bit NCBI2na and records
/// ambiguous IUPAC symbols separately as NCBI4na runs. The compact one-word
/// format is used when possible; the extended two-word format is used when a
/// run cannot fit the old 24-bit-position/4-bit-length layout.
pub fn encode_ncbi2na_ambiguity_data(seq: &[u8]) -> Vec<u8> {
    let mut runs = Vec::new();
    let mut i = 0usize;
    while i < seq.len() {
        let ncbi4na = iupacna_to_ncbi4na_base(seq[i]);
        if is_unambiguous_ncbi4na(ncbi4na) {
            i += 1;
            continue;
        }

        let start = i;
        i += 1;
        while i < seq.len() && iupacna_to_ncbi4na_base(seq[i]) == ncbi4na {
            i += 1;
        }
        runs.push((start as u32, (i - start) as u32, ncbi4na as u32));
    }

    encode_ncbi2na_ambiguity_runs(&runs)
}

fn encode_ncbi2na_ambiguity_runs(runs: &[(u32, u32, u32)]) -> Vec<u8> {
    if runs.is_empty() {
        return Vec::new();
    }

    let new_format = runs.iter().any(|&(position, run_len, _)| {
        position > OLD_AMBIGUITY_MAX_POSITION || run_len > OLD_AMBIGUITY_MAX_RUN
    });
    let max_run = if new_format {
        NEW_AMBIGUITY_MAX_RUN
    } else {
        OLD_AMBIGUITY_MAX_RUN
    };

    let record_count: u32 = runs
        .iter()
        .map(|&(_, run_len, _)| (run_len + max_run - 1) / max_run)
        .sum();
    let mut out = Vec::with_capacity(4 + record_count as usize * if new_format { 8 } else { 4 });
    let header = if new_format {
        0x8000_0000 | record_count
    } else {
        record_count
    };
    out.extend_from_slice(&header.to_be_bytes());

    for &(mut position, mut run_len, ncbi4na) in runs {
        while run_len > 0 {
            let chunk_len = run_len.min(max_run);
            let run_len_minus_one = chunk_len - 1;
            if new_format {
                let word = (ncbi4na << 28) | (run_len_minus_one << 16);
                out.extend_from_slice(&word.to_be_bytes());
                out.extend_from_slice(&position.to_be_bytes());
            } else {
                let word = (ncbi4na << 28) | (run_len_minus_one << 24) | position;
                out.extend_from_slice(&word.to_be_bytes());
            }
            position = position.saturating_add(chunk_len);
            run_len -= chunk_len;
        }
    }

    out
}

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
        for &ch in b"ACGT" {
            let blastna = IUPACNA_TO_BLASTNA[ch as usize];
            let back = BLASTNA_TO_IUPACNA[blastna as usize] as u8;
            assert_eq!(back, ch, "Roundtrip failed for {}", ch as char);
        }
    }

    // -- test_aminoacid_to_ncbistdaa --
    #[test]
    fn test_aminoacid_to_ncbistdaa() {
        // Verify standard amino acids
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'A' as usize], NCBISTDAA_A);
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'C' as usize], NCBISTDAA_C);
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'D' as usize], NCBISTDAA_D);
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'E' as usize], NCBISTDAA_E);
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'F' as usize], 6);
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'G' as usize], 7);
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'H' as usize], 8);
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'I' as usize], NCBISTDAA_I);
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'K' as usize], 10);
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'L' as usize], NCBISTDAA_L);
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'M' as usize], 12);
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'N' as usize], NCBISTDAA_N);
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'P' as usize], 14);
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'Q' as usize], NCBISTDAA_Q);
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'R' as usize], 16);
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'S' as usize], 17);
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'T' as usize], 18);
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'V' as usize], 19);
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'W' as usize], 20);
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'Y' as usize], 22);
        // Stop codon
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'*' as usize], NCBISTDAA_STOP);
        // Lowercase matches uppercase
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'a' as usize], NCBISTDAA_A);
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'm' as usize], 12);
    }

    #[test]
    fn test_aminoacid_to_ncbistdaa_base_and_sequence() {
        assert_eq!(aminoacid_to_ncbistdaa_base(b'A'), 1);
        assert_eq!(aminoacid_to_ncbistdaa_base(b'm'), 12);
        assert_eq!(aminoacid_to_ncbistdaa_base(b'*'), NCBISTDAA_STOP);
        assert_eq!(
            encode_ncbistdaa_sequence(b"Am*"),
            vec![1, 12, NCBISTDAA_STOP]
        );
    }

    #[test]
    fn test_ncbistdaa_residue_classification_helpers() {
        assert!(is_ncbistdaa_standard_residue(1)); // A
        assert!(is_ncbistdaa_standard_residue(22)); // Y
        assert!(!is_ncbistdaa_standard_residue(2)); // B
        assert!(!is_ncbistdaa_standard_residue(NCBISTDAA_X));
        assert!(!is_ncbistdaa_standard_residue(NCBISTDAA_U));

        assert!(is_ncbistdaa_composition_residue(1));
        assert!(is_ncbistdaa_composition_residue(NCBISTDAA_U));
        assert!(!is_ncbistdaa_composition_residue(NCBISTDAA_X));
    }

    #[test]
    fn test_ncbistdaa_to_aminoacid_base_maps_out_of_range_to_x() {
        assert_eq!(ncbistdaa_to_aminoacid_base(1), b'A');
        assert_eq!(ncbistdaa_to_aminoacid_base(NCBISTDAA_STOP), b'*');
        assert_eq!(ncbistdaa_to_aminoacid_base(255), b'X');
    }

    #[test]
    fn test_ncbistdaa_to_aminoacid_char_helper() {
        assert_eq!(ncbistdaa_to_aminoacid_char(1), 'A');
        assert_eq!(ncbistdaa_to_aminoacid_char(NCBISTDAA_STOP), '*');
        assert_eq!(ncbistdaa_to_aminoacid_char(255), 'X');
    }

    #[test]
    fn test_ncbistdaa_to_aminoacid_sequence_helper() {
        assert_eq!(
            ncbistdaa_to_aminoacid_sequence(&[1, 12, NCBISTDAA_STOP, 28]),
            b"AM*X"
        );
    }

    #[test]
    fn test_ncbistdaa_to_aminoacid_string_helper() {
        assert_eq!(
            ncbistdaa_to_aminoacid_string(&[1, 12, NCBISTDAA_STOP, 28]),
            "AM*X"
        );
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
        assert_eq!(ncbi4na_to_iupacna_base(1), b'A');
        assert_eq!(ncbi4na_to_iupacna_base(2), b'C');
        assert_eq!(ncbi4na_to_iupacna_base(4), b'G');
        assert_eq!(ncbi4na_to_iupacna_base(8), b'T');
        // Gap = 0
        assert_eq!(ncbi4na_to_iupacna_base(0), b'-');
        // N = 15
        assert_eq!(ncbi4na_to_iupacna_base(15), b'N');
        assert_eq!(ncbi4na_to_iupacna_base(99), b'N');
        assert_eq!(ncbi4na_to_iupacna_char(3), 'M');
        assert_eq!(
            ncbi4na_to_iupacna_sequence(&[0, 1, 2, 4, 8, 15, 99]),
            b"-ACGTNN"
        );
        assert_eq!(
            ncbi4na_to_iupacna_string(&[0, 1, 2, 4, 8, 15, 99]),
            "-ACGTNN"
        );
    }

    #[test]
    fn test_special_amino_acids() {
        // U = Selenocysteine
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'U' as usize], NCBISTDAA_U);
        assert_eq!(NCBISTDAA_TO_AMINOACID[NCBISTDAA_U as usize] as u8, b'U');
        // O = Pyrrolysine
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'O' as usize], NCBISTDAA_O);
        assert_eq!(NCBISTDAA_TO_AMINOACID[NCBISTDAA_O as usize] as u8, b'O');
        // J = Leucine or Isoleucine
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'J' as usize], NCBISTDAA_J);
        assert_eq!(NCBISTDAA_TO_AMINOACID[NCBISTDAA_J as usize] as u8, b'J');
        // X = unknown
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'X' as usize], NCBISTDAA_X);
        assert_eq!(NCBISTDAA_TO_AMINOACID[NCBISTDAA_X as usize] as u8, b'X');
        // B = Asp or Asn
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'B' as usize], NCBISTDAA_B);
        assert_eq!(NCBISTDAA_TO_AMINOACID[NCBISTDAA_B as usize] as u8, b'B');
        // Z = Glu or Gln
        assert_eq!(AMINOACID_TO_NCBISTDAA[b'Z' as usize], NCBISTDAA_Z);
        assert_eq!(NCBISTDAA_TO_AMINOACID[NCBISTDAA_Z as usize] as u8, b'Z');
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

    #[test]
    fn test_iupacna_to_ncbi4na_base_maps_u_and_invalid_to_blast_input_values() {
        assert_eq!(iupacna_to_ncbi4na_base(b'U'), 8);
        assert_eq!(iupacna_to_ncbi4na_base(b'u'), 8);
        assert_eq!(iupacna_to_ncbi4na_base(b'?'), 15);
    }

    #[test]
    fn test_encode_ncbi4na_sequence_helper() {
        assert_eq!(encode_ncbi4na_sequence(b"ACGTU?"), vec![1, 2, 4, 8, 8, 15]);
    }

    #[test]
    fn test_iupacna_to_blastna_base_maps_u_to_t() {
        assert_eq!(iupacna_to_blastna_base(b'U'), 3);
        assert_eq!(iupacna_to_blastna_base(b'u'), 3);
        assert_eq!(iupacna_to_blastna_base(b'R'), 4);
        assert_eq!(iupacna_to_blastna_base(b'?'), 15);
    }

    #[test]
    fn test_encode_blastna_sequence_helper() {
        assert_eq!(
            encode_blastna_sequence(b"ACGTURYN?"),
            vec![0, 1, 2, 3, 3, 4, 5, 14, 15]
        );
    }

    #[test]
    fn test_blastna_ncbi4na_base_helpers() {
        assert_eq!(blastna_to_ncbi4na_base(0), 1);
        assert_eq!(blastna_to_ncbi4na_base(14), 15);
        assert_eq!(blastna_to_ncbi4na_base(99), 0);
        assert_eq!(ncbi4na_to_blastna_base(1), 0);
        assert_eq!(ncbi4na_to_blastna_base(15), 14);
        assert_eq!(ncbi4na_to_blastna_base(99), 15);
    }

    #[test]
    fn test_blastna_to_ncbi4na_sequence_helper() {
        assert_eq!(
            blastna_to_ncbi4na_sequence(&[0, 1, 2, 3, 4, 14, 99]),
            vec![1, 2, 4, 8, 5, 15, 0]
        );
    }

    #[test]
    fn test_blastna_overlap_and_degeneracy_helpers() {
        assert!(blastna_bases_overlap(4, 0)); // R overlaps A.
        assert!(blastna_bases_overlap(4, 2)); // R overlaps G.
        assert!(!blastna_bases_overlap(4, 1)); // R does not overlap C.
        assert_eq!(blastna_degeneracy(0), 1);
        assert_eq!(blastna_degeneracy(4), 2);
        assert_eq!(blastna_degeneracy(14), 4);
        assert_eq!(blastna_degeneracy(15), 0);
    }

    #[test]
    fn test_blastna_pair_score_helper() {
        assert_eq!(blastna_pair_score(0, 0, 1, -3), 1);
        assert_eq!(blastna_pair_score(0, 1, 1, -3), -3);
        assert_eq!(blastna_pair_score(4, 0, 1, -3), -1);
        assert_eq!(blastna_pair_score(0, 4, 1, -3), -1);
        assert_eq!(blastna_pair_score(4, 1, 1, -3), -3);
        assert_eq!(blastna_pair_score(15, 0, 1, -3), i32::MIN / 2);
    }

    #[test]
    fn test_complement_blastna_base_helper() {
        assert_eq!(complement_blastna_base(0), 3);
        assert_eq!(complement_blastna_base(1), 2);
        assert_eq!(complement_blastna_base(4), 5);
        assert_eq!(complement_blastna_base(10), 13);
        assert_eq!(complement_blastna_base(14), 14);
        assert_eq!(complement_blastna_base(99), 15);
    }

    #[test]
    fn test_reverse_complement_blastna_sequence_helper() {
        assert_eq!(
            reverse_complement_blastna_sequence(&[0, 1, 2, 3, 4, 14]),
            vec![14, 5, 0, 1, 2, 3]
        );
    }

    #[test]
    fn test_blastna_to_iupacna_base_helper() {
        assert_eq!(blastna_to_iupacna_base(0), b'A');
        assert_eq!(blastna_to_iupacna_base(4), b'R');
        assert_eq!(blastna_to_iupacna_base(10), b'B');
        assert_eq!(blastna_to_iupacna_base(14), b'N');
        assert_eq!(blastna_to_iupacna_base(15), b'-');
        assert_eq!(blastna_to_iupacna_base(99), b'-');
    }

    #[test]
    fn test_blastna_to_iupacna_char_helper() {
        assert_eq!(blastna_to_iupacna_char(0), 'A');
        assert_eq!(blastna_to_iupacna_char(15), '-');
    }

    #[test]
    fn test_blastna_to_iupacna_sequence_helper() {
        assert_eq!(
            blastna_to_iupacna_sequence(&[0, 1, 2, 3, 4, 5, 14, 15, 99]),
            b"ACGTRYN--"
        );
    }

    #[test]
    fn test_blastna_to_iupacna_string_helper() {
        assert_eq!(
            blastna_to_iupacna_string(&[0, 1, 2, 3, 4, 5, 14, 15, 99]),
            "ACGTRYN--"
        );
    }

    #[test]
    fn test_complement_iupacna_ascii_base_helper() {
        assert_eq!(complement_iupacna_ascii_base(b'A'), b'T');
        assert_eq!(complement_iupacna_ascii_base(b'a'), b'T');
        assert_eq!(complement_iupacna_ascii_base(b'U'), b'A');
        assert_eq!(complement_iupacna_ascii_base(b'u'), b'A');
        assert_eq!(complement_iupacna_ascii_base(b'R'), b'Y');
        assert_eq!(complement_iupacna_ascii_base(b'r'), b'Y');
        assert_eq!(complement_iupacna_ascii_base(b'B'), b'V');
        assert_eq!(complement_iupacna_ascii_base(b'b'), b'V');
        assert_eq!(complement_iupacna_ascii_base(b'?'), b'N');
    }

    #[test]
    fn test_reverse_complement_iupacna_sequence_helper() {
        assert_eq!(
            reverse_complement_iupacna_sequence(b"rYmKwsBdHvN"),
            b"NBDHVSWMKRY"
        );
        assert_eq!(reverse_complement_iupacna_sequence(b"UuTt"), b"AAAA");
    }

    #[test]
    fn test_complement_ncbi4na_base_helper() {
        assert_eq!(complement_ncbi4na_base(1), 8);
        assert_eq!(complement_ncbi4na_base(8), 1);
        assert_eq!(complement_ncbi4na_base(2), 4);
        assert_eq!(complement_ncbi4na_base(4), 2);
        assert_eq!(complement_ncbi4na_base(15), 15);
        assert_eq!(complement_ncbi4na_base(0xff), 15);
    }

    #[test]
    fn test_ncbi2na_base_at_helper() {
        let packed = vec![0b0001_1011, 0b1110_0100];
        assert_eq!(ncbi2na_base_at(&packed, 0), 0);
        assert_eq!(ncbi2na_base_at(&packed, 1), 1);
        assert_eq!(ncbi2na_base_at(&packed, 2), 2);
        assert_eq!(ncbi2na_base_at(&packed, 3), 3);
        assert_eq!(ncbi2na_base_at(&packed, 4), 3);
        assert_eq!(ncbi2na_base_at(&packed, 5), 2);
        assert_eq!(ncbi2na_base_at(&packed, 6), 1);
        assert_eq!(ncbi2na_base_at(&packed, 7), 0);
    }

    #[test]
    fn test_blastna_or_iupac_to_ncbi2na_base_preserves_legacy_dust_mapping() {
        assert_eq!(blastna_or_iupac_to_ncbi2na_base(b'A'), 0);
        assert_eq!(blastna_or_iupac_to_ncbi2na_base(b'c'), 1);
        assert_eq!(blastna_or_iupac_to_ncbi2na_base(b'U'), 3);
        assert_eq!(blastna_or_iupac_to_ncbi2na_base(6), 2);
        assert_eq!(blastna_or_iupac_to_ncbi2na_base(b'?'), 0);
    }

    #[test]
    fn test_pack_ncbi2na_bases_has_no_db_remainder_sentinel() {
        assert_eq!(pack_ncbi2na_bases(&[0, 1, 2, 3]), vec![0b0001_1011]);
        assert_eq!(
            pack_ncbi2na_bases(&[3, 2, 1, 0, 1]),
            vec![0b1110_0100, 0b0100_0000]
        );
        assert_eq!(pack_ncbi2na_bases(&[]), Vec::<u8>::new());
    }

    #[test]
    fn test_encode_ncbi2na_sequence_packs_bases_and_remainder() {
        assert_eq!(encode_ncbi2na_sequence(b"ACGT"), vec![0b0001_1011, 0]);
        assert_eq!(encode_ncbi2na_sequence(b"ACG"), vec![0b0001_1000 | 3]);
        assert_eq!(encode_ncbi2na_sequence(b"AU"), vec![0b0011_0000 | 2]);
    }

    #[test]
    fn test_encode_ncbi2na_ambiguity_data_uses_old_format_when_possible() {
        let encoded = encode_ncbi2na_ambiguity_data(b"ACGTRRNNY");
        let words: Vec<u32> = encoded
            .chunks_exact(4)
            .map(|chunk| u32::from_be_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]))
            .collect();

        assert_eq!(words[0], 3);
        assert_eq!(words[1], (5u32 << 28) | (1u32 << 24) | 4u32);
        assert_eq!(words[2], (15u32 << 28) | (1u32 << 24) | 6u32);
        assert_eq!(words[3], (10u32 << 28) | 8u32);
    }

    #[test]
    fn test_encode_ncbi2na_ambiguity_runs_uses_new_format_for_large_position() {
        let encoded = encode_ncbi2na_ambiguity_runs(&[(0x0100_0000, 2, 15)]);
        let words: Vec<u32> = encoded
            .chunks_exact(4)
            .map(|chunk| u32::from_be_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]))
            .collect();

        assert_eq!(words[0], 0x8000_0001);
        assert_eq!(words[1], (15u32 << 28) | (1u32 << 16));
        assert_eq!(words[2], 0x0100_0000);
    }

    #[test]
    fn test_encode_ncbi2na_ambiguity_runs_uses_new_format_for_long_run() {
        let encoded = encode_ncbi2na_ambiguity_runs(&[(2, 4097, 15)]);
        let words: Vec<u32> = encoded
            .chunks_exact(4)
            .map(|chunk| u32::from_be_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]))
            .collect();

        assert_eq!(words[0], 0x8000_0002);
        assert_eq!(words[1], (15u32 << 28) | (4095u32 << 16));
        assert_eq!(words[2], 2);
        assert_eq!(words[3], 15u32 << 28);
        assert_eq!(words[4], 4098);
    }
}
