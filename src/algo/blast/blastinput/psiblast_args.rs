// Upstream sources:
// - ncbi-blast-2.17.0+-src/c++/include/algo/blast/blastinput/psiblast_args.hpp
// - ncbi-blast-2.17.0+-src/c++/src/algo/blast/blastinput/psiblast_args.cpp

use super::blast_args::{BlastAppArgs, PsiBlastArgs};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PsiBlastAppArgs {
    pub base: BlastAppArgs,
    pub psi_blast_args: Option<PsiBlastArgs>,
}
