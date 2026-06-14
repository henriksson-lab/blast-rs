// Upstream sources:
// - ncbi-blast-2.17.0+-src/c++/include/algo/blast/blastinput/blast_scope_src.hpp
// - ncbi-blast-2.17.0+-src/c++/src/algo/blast/blastinput/blast_scope_src.cpp

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DataLoaderConfigOpt {
    UseBlastDbDataLoader,
    UseGenbankDataLoader,
    UseNoDataLoaders,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DataLoaderConfig {
    pub use_blast_dbs: bool,
    pub blast_db_name: String,
    pub is_loading_proteins: bool,
    pub use_genbank: bool,
    pub use_fixed_size_slices: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BlastScopeSource {
    pub object_manager: Option<String>,
    pub config: DataLoaderConfig,
    pub blast_db_loader_name: String,
    pub genbank_loader_name: String,
}
