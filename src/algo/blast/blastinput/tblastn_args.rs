// Upstream sources:
// - ncbi-blast-2.17.0+-src/c++/include/algo/blast/blastinput/tblastn_args.hpp
// - ncbi-blast-2.17.0+-src/c++/src/algo/blast/blastinput/tblastn_args.cpp

use super::blast_args::{BlastAppArgs, PsiBlastArgs};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TblastnAppArgs {
    pub base: BlastAppArgs,
    pub psi_blast_args: Option<PsiBlastArgs>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TblastnNodeArgs {
    pub base: TblastnAppArgs,
    pub output_stream: String,
    pub input_stream: Option<String>,
}
