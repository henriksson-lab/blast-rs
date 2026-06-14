// Upstream sources:
// - ncbi-blast-2.17.0+-src/c++/include/algo/blast/blastinput/blastp_args.hpp
// - ncbi-blast-2.17.0+-src/c++/src/algo/blast/blastinput/blastp_args.cpp

use super::blast_args::BlastAppArgs;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BlastpAppArgs {
    pub base: BlastAppArgs,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BlastpNodeArgs {
    pub base: BlastpAppArgs,
    pub output_stream: String,
    pub input_stream: Option<String>,
}
