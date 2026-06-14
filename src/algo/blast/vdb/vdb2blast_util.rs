// Upstream sources:
// - ncbi-blast-2.17.0+-src/c++/include/algo/blast/vdb/vdb2blast_util.hpp
// - ncbi-blast-2.17.0+-src/c++/src/algo/blast/vdb/vdb2blast_util.cpp

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VdbIdType {
    SeqId,
    Accession,
    SpotId,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VdbBlastUtil {
    pub own_seq_src: bool,
    pub all_runs: String,
    pub seq_src: Option<String>,
    pub is_csra_util: bool,
    pub include_filtered_reads: bool,
}
