// Upstream sources:
// - ncbi-blast-2.17.0+-src/c++/include/algo/blast/format/data4xmlformat.hpp
// - ncbi-blast-2.17.0+-src/c++/src/algo/blast/format/data4xmlformat.cpp

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CmdLineBlastXmlReportData {
    pub queries: Option<String>,
    pub options: Option<String>,
    pub db_name: String,
    pub query_genetic_code: i32,
    pub db_genetic_code: i32,
    pub ancillary_data: Vec<String>,
    pub alignments: Vec<String>,
    pub masks: Vec<String>,
    pub no_hits_found: bool,
    pub errors: Vec<String>,
    pub matrix: Vec<Vec<i32>>,
    pub num_sequences: i32,
    pub num_bases: i64,
}
