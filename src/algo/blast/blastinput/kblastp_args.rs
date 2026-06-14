// Upstream sources:
// - ncbi-blast-2.17.0+-src/c++/include/algo/blast/blastinput/kblastp_args.hpp
// - ncbi-blast-2.17.0+-src/c++/src/algo/blast/blastinput/kblastp_args.cpp

use super::blast_args::{BlastAppArgs, KBlastpArgs};

#[derive(Debug, Clone, PartialEq)]
pub struct KBlastpAppArgs {
    pub base: BlastAppArgs,
    pub kblastp_args: Option<KBlastpArgs>,
}
