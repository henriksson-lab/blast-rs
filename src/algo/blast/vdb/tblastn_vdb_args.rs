// Upstream sources:
// - ncbi-blast-2.17.0+-src/c++/include/algo/blast/vdb/tblastn_vdb_args.hpp
// - ncbi-blast-2.17.0+-src/c++/src/algo/blast/vdb/tblastn_vdb_args.cpp

use crate::algo::blast::blastinput::blast_args::{BlastAppArgs, PsiBlastArgs};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TblastnVdbAppArgs {
    pub base: BlastAppArgs,
    pub vdb_search_mode_args: Option<super::blastn_vdb_args::SraSearchModeArgs>,
    pub psi_blast_args: Option<PsiBlastArgs>,
}
