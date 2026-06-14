// Upstream sources:
// - ncbi-blast-2.17.0+-src/c++/include/algo/blast/format/blastxml2_format.hpp
// - ncbi-blast-2.17.0+-src/c++/src/algo/blast/format/blastxml2_format.cpp

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BlastXml2ReportData {
    pub program: String,
    pub database_name: String,
    pub query_id: String,
    pub messages: Vec<String>,
}
