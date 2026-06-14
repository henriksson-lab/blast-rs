// Upstream sources:
// - ncbi-blast-2.17.0+-src/c++/include/algo/blast/blastinput/deltablast_args.hpp
// - ncbi-blast-2.17.0+-src/c++/src/algo/blast/blastinput/deltablast_args.cpp

use super::blast_args::{BlastAppArgs, DeltaBlastArgs, PsiBlastArgs};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DeltaBlastAppArgs {
    pub base: BlastAppArgs,
    pub delta_blast_args: Option<DeltaBlastArgs>,
    pub psi_blast_args: Option<PsiBlastArgs>,
}
