mod encoding;
mod fasta;

pub use encoding::{aminoacid_to_ncbistdaa, iupacna_to_blastna};
pub use fasta::{parse_fasta, parse_fasta_with_default_id, FastaRecord};
