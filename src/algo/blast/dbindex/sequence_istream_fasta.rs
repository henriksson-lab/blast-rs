// Upstream sources:
// - ncbi-blast-2.17.0+-src/c++/include/algo/blast/dbindex/sequence_istream.hpp
// - ncbi-blast-2.17.0+-src/c++/include/algo/blast/dbindex/sequence_istream_fasta.hpp
// - ncbi-blast-2.17.0+-src/c++/src/algo/blast/dbindex/sequence_istream_fasta.cpp

use super::sequence_istream_bdb::SequenceStreamData;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SequenceIStreamFasta {
    pub stream_allocated: bool,
    pub istream: Option<String>,
    pub fasta_reader: Option<String>,
    pub seq_positions: Vec<u64>,
    pub name: String,
    pub cache: Option<SequenceStreamData>,
    pub use_cache: bool,
}
