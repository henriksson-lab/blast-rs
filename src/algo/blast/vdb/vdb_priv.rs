// Upstream sources:
// - ncbi-blast-2.17.0+-src/c++/src/algo/blast/vdb/vdb_priv.h
// - ncbi-blast-2.17.0+-src/c++/src/algo/blast/vdb/vdb_priv.c

use super::error_priv::VdbSrcErrMsg;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Vdb2naIcReader {
    pub reader: Option<String>,
    pub row_id: i64,
    pub start: u64,
    pub len: u64,
    pub buffer: Vec<u8>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VdbPartialFetchingRanges {
    pub ranges: Vec<i32>,
    pub num_ranges: i32,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VdbData {
    pub vdb_manager: Option<String>,
    pub run_set: Option<String>,
    pub ref_set: Option<String>,
    pub reader_2na: Option<String>,
    pub reader_4na: Option<String>,
    pub reader_cache_2na: Option<Vdb2naIcReader>,
    pub run_set_name: String,
    pub num_seqs: u64,
    pub total_length: u64,
    pub status: VdbSrcErrMsg,
    pub include_filtered_reads: bool,
    pub range_list: Option<VdbPartialFetchingRanges>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VdbNewArgs {
    pub vdb_accessions: String,
    pub is_protein: bool,
    pub use_cache: bool,
    pub include_filtered_reads: bool,
    pub range_list: Option<VdbPartialFetchingRanges>,
}
