// Upstream source:
// - ncbi-blast-2.17.0+-src/c++/src/algo/blast/dbindex/dbindex_search.cpp

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DbIndexSearchContext {
    pub index_name: String,
    pub word_size: u64,
    pub lookup_table: Option<String>,
    pub sequence_blocks: Vec<String>,
    pub diagnostics: Option<String>,
}
