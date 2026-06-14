// Upstream sources:
// - ncbi-blast-2.17.0+-src/c++/include/algo/blast/dbindex/sequence_istream.hpp
// - ncbi-blast-2.17.0+-src/c++/include/algo/blast/dbindex/sequence_istream_bdb.hpp
// - ncbi-blast-2.17.0+-src/c++/src/algo/blast/dbindex/sequence_istream_bdb.cpp

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SequenceStreamData {
    pub seq_entry: Option<String>,
    pub mask_locs: Vec<String>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SequenceIStreamExceptionCode {
    Eof,
    BadSequence,
    Io,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SequenceIStreamBlastDb {
    pub seqdb: Option<String>,
    pub oid: i32,
    pub filter_algo_id: i32,
    pub use_filter: bool,
}
