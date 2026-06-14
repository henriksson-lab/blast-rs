// Upstream sources:
// - ncbi-blast-2.17.0+-src/c++/include/algo/blast/blastinput/igblastn_args.hpp
// - ncbi-blast-2.17.0+-src/c++/src/algo/blast/blastinput/igblastn_args.cpp

use super::blast_args::{BlastAppArgs, IgBlastArgs};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IgBlastnAppArgs {
    pub base: BlastAppArgs,
    pub ig_blast_args: Option<IgBlastArgs>,
}
