//! Rust equivalent of blast_encoding.c
//! These tables will replace the C versions when the FFI layer is removed.

pub const BLASTNA_SIZE: usize = 16;
pub const BLASTAA_SIZE: usize = 28;

pub static NCBI4NA_TO_BLASTNA: [u8; BLASTNA_SIZE] = [
    15, 0, 1, 6, 2, 4, 9, 13, 3, 8, 5, 12, 7, 11, 10, 14,
];


pub static BLASTNA_TO_NCBI4NA: [u8; BLASTNA_SIZE] = [
    1, 2, 4, 8, 5, 10, 3, 12, 9, 6, 14, 13, 11, 7, 15, 0,
];


pub static BLASTNA_TO_IUPACNA: [i8; BLASTNA_SIZE] = [
    b'A' as i8, b'C' as i8, b'G' as i8, b'T' as i8,
    b'R' as i8, b'Y' as i8, b'M' as i8, b'K' as i8,
    b'W' as i8, b'S' as i8, b'B' as i8, b'D' as i8,
    b'H' as i8, b'V' as i8, b'N' as i8, b'-' as i8,
];


pub static NCBI4NA_TO_IUPACNA: [i8; BLASTNA_SIZE] = [
    b'-' as i8, b'A' as i8, b'C' as i8, b'M' as i8,
    b'G' as i8, b'R' as i8, b'S' as i8, b'V' as i8,
    b'T' as i8, b'W' as i8, b'Y' as i8, b'H' as i8,
    b'K' as i8, b'D' as i8, b'B' as i8, b'N' as i8,
];


pub static IUPACNA_TO_BLASTNA: [u8; 128] = [
    15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
    15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
    15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
    15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
    15, 0,10, 1,11,15,15, 2,12,15,15, 7,15, 6,14,15,
    15,15, 4, 9, 3,15,13, 8,15, 5,15,15,15,15,15,15,
    15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
    15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
];


pub static IUPACNA_TO_NCBI4NA: [u8; 128] = [
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1,14, 2,13, 0, 0, 4,11, 0, 0,12, 0, 3,15, 0,
    0, 0, 5, 6, 8, 0, 7, 9, 0,10, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
];


pub static AMINOACID_TO_NCBISTDAA: [u8; 128] = [
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,25, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9,27,10,11,12,13,26,
    14,15,16,17,18,24,19,20,21,22,23, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
];


pub static NCBISTDAA_TO_AMINOACID: [i8; BLASTAA_SIZE] = [
    b'-' as i8, b'A' as i8, b'B' as i8, b'C' as i8,
    b'D' as i8, b'E' as i8, b'F' as i8, b'G' as i8,
    b'H' as i8, b'I' as i8, b'K' as i8, b'L' as i8,
    b'M' as i8, b'N' as i8, b'P' as i8, b'Q' as i8,
    b'R' as i8, b'S' as i8, b'T' as i8, b'V' as i8,
    b'W' as i8, b'X' as i8, b'Y' as i8, b'Z' as i8,
    b'U' as i8, b'*' as i8, b'O' as i8, b'J' as i8,
];


pub static kProtSentinel: u8 = 0;


pub static kNuclSentinel: u8 = 0xF;
