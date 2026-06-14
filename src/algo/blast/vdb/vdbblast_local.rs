// Upstream sources:
// - ncbi-blast-2.17.0+-src/c++/include/algo/blast/vdb/vdbblast_local.hpp
// - ncbi-blast-2.17.0+-src/c++/src/algo/blast/vdb/vdbblast_local.cpp

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LocalVdbStruct {
    pub chunks_for_thread: Vec<Vec<String>>,
    pub total_num_seqs: u64,
    pub total_length: u64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SraSearchMode {
    SraSearch,
    CsraSearch,
    WgsSearch,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LocalVdbBlast {
    pub query_vector: Option<String>,
    pub options_handle: Option<String>,
    pub total_num_seqs: u64,
    pub total_length: u64,
    pub chunks_for_thread: Vec<Vec<String>>,
    pub num_threads: u32,
    pub num_extensions: i32,
    pub include_filtered_reads: bool,
    pub pssm: Option<String>,
}
