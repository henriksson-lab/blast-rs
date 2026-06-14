// Upstream sources:
// - ncbi-blast-2.17.0+-src/c++/include/algo/blast/blastinput/blast_fasta_input.hpp
// - ncbi-blast-2.17.0+-src/c++/src/algo/blast/blastinput/blast_fasta_input.cpp

use super::blast_input::BlastInputSourceConfig;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BlastFastaInputSource {
    pub config: BlastInputSourceConfig,
    pub line_reader: Option<String>,
    pub input_reader: Option<String>,
    pub read_proteins: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ShortReadInputFormat {
    Fasta,
    Fastc,
    Fastq,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ShortReadFastaInputSource {
    pub bases_added: u64,
    pub seq_buffer_len: u64,
    pub line_reader: Option<String>,
    pub second_line_reader: Option<String>,
    pub sequence: String,
    pub is_paired: bool,
    pub format: ShortReadInputFormat,
    pub id: u32,
    pub parse_seq_ids: bool,
}
