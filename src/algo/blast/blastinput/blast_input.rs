// Upstream sources:
// - ncbi-blast-2.17.0+-src/c++/include/algo/blast/blastinput/blast_input.hpp
// - ncbi-blast-2.17.0+-src/c++/src/algo/blast/blastinput/blast_input.cpp

use super::blast_scope_src::DataLoaderConfig;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NaStrand {
    Unknown,
    Plus,
    Minus,
    Both,
    Other,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SeqRange {
    pub from: u64,
    pub to: u64,
    pub empty: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BlastInputSourceConfig {
    pub strand: NaStrand,
    pub lower_case_mask: bool,
    pub believe_deflines: bool,
    pub skip_seq_check: bool,
    pub range: SeqRange,
    pub data_loader_config: DataLoaderConfig,
    pub retrieve_seq_data: bool,
    pub local_id_counter: i32,
    pub seq_len_threshold_to_guess: u32,
    pub local_id_prefix: String,
    pub gaps_to_ns: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InputExceptionCode {
    InvalidInput,
    EmptyInput,
    MissingProteinId,
    MissingNucleotideId,
    MissingSequenceData,
    InvalidSequenceData,
    InvalidRange,
    SequenceDataError,
    Unknown,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BlastInputSource {
    pub source_name: String,
    pub exhausted: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BlastInput {
    pub source: Option<Box<BlastInputSource>>,
    pub batch_size: u64,
    pub num_seqs: i64,
    pub total_length: i64,
}
