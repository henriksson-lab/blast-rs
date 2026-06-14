// Upstream sources:
// - ncbi-blast-2.17.0+-src/c++/include/algo/blast/blastinput/tblastx_args.hpp
// - ncbi-blast-2.17.0+-src/c++/src/algo/blast/blastinput/tblastx_args.cpp

use super::blast_args::BlastAppArgs;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TblastxAppArgs {
    pub base: BlastAppArgs,
}
