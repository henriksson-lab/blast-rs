// Upstream sources:
// - ncbi-blast-2.17.0+-src/c++/include/algo/blast/format/data4xml2format.hpp
// - ncbi-blast-2.17.0+-src/c++/src/algo/blast/format/data4xml2format.cpp

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CmdLineBlastXml2ReportData {
    pub query: Option<String>,
    pub options: Option<String>,
    pub scope: Option<String>,
    pub db_name: String,
    pub num_sequences: i64,
    pub num_bases: i64,
    pub tax_db_found: bool,
    pub is_bl2seq: bool,
    pub is_iterative: bool,
    pub ancillary_data: Vec<String>,
    pub alignments: Vec<String>,
    pub errors: Vec<String>,
    pub matrix: Option<String>,
    pub subject_ids: Vec<String>,
    pub query_masks: Vec<String>,
}
