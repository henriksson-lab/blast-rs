mod encoding;
mod fasta;

pub use encoding::{iupacna_to_blastna, aminoacid_to_ncbistdaa};
pub use fasta::{FastaRecord, parse_fasta};
