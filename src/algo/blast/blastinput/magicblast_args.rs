// Upstream sources:
// - ncbi-blast-2.17.0+-src/c++/include/algo/blast/blastinput/magicblast_args.hpp
// - ncbi-blast-2.17.0+-src/c++/src/algo/blast/blastinput/magicblast_args.cpp

use super::blast_args::BlastAppArgs;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MagicBlastAppArgs {
    pub base: BlastAppArgs,
}
