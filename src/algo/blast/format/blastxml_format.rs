// Upstream sources:
// - ncbi-blast-2.17.0+-src/c++/include/algo/blast/format/blastxml_format.hpp
// - ncbi-blast-2.17.0+-src/c++/src/algo/blast/format/blastxml_format.cpp

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BlastXmlReportData {
    pub program: String,
    pub database_name: String,
    pub query_ids: Vec<String>,
    pub messages: Vec<String>,
}
