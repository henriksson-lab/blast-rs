// Upstream sources:
// - ncbi-blast-2.17.0+-src/c++/include/algo/blast/format/build_archive.hpp
// - ncbi-blast-2.17.0+-src/c++/src/algo/blast/format/build_archive.cpp

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BlastArchiveBuildInputs {
    pub queries: Option<String>,
    pub options: Option<String>,
    pub alignments: Option<String>,
    pub masks: Vec<String>,
    pub messages: Vec<String>,
    pub pssm: Option<String>,
    pub num_iters: u32,
}
