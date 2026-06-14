// Upstream sources:
// - ncbi-blast-2.17.0+-src/c++/include/algo/blast/blastinput/rmblastn_args.hpp
// - ncbi-blast-2.17.0+-src/c++/src/algo/blast/blastinput/rmblastn_args.cpp

use super::blast_args::BlastAppArgs;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RmBlastnAppArgs {
    pub base: BlastAppArgs,
}
