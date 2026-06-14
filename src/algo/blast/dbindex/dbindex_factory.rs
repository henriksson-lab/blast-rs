// Upstream source:
// - ncbi-blast-2.17.0+-src/c++/src/algo/blast/dbindex/dbindex_factory.cpp

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DbIndexFactoryRequest {
    pub index_name: String,
    pub old_style_index: bool,
    pub use_index: bool,
    pub molecule_type: Option<String>,
}
