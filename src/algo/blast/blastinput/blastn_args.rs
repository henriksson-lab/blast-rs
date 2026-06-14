// Upstream sources:
// - ncbi-blast-2.17.0+-src/c++/include/algo/blast/blastinput/blastn_args.hpp
// - ncbi-blast-2.17.0+-src/c++/src/algo/blast/blastinput/blastn_args.cpp

use super::blast_args::BlastAppArgs;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BlastnAppArgs {
    pub base: BlastAppArgs,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BlastnNodeArgs {
    pub base: BlastnAppArgs,
    pub output_stream: String,
    pub input_stream: Option<String>,
}
