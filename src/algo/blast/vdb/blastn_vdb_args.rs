// Upstream sources:
// - ncbi-blast-2.17.0+-src/c++/include/algo/blast/vdb/blastn_vdb_args.hpp
// - ncbi-blast-2.17.0+-src/c++/src/algo/blast/vdb/blastn_vdb_args.cpp

use crate::algo::blast::blastinput::blast_args::{BlastAppArgs, BlastDatabaseArgs};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BlastVDatabaseArgs {
    pub base: BlastDatabaseArgs,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SraSearchModeArgs {
    pub search_mode: super::vdbblast_local::SraSearchMode,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BlastnVdbAppArgs {
    pub base: BlastAppArgs,
    pub vdb_search_mode_args: Option<SraSearchModeArgs>,
}
