// Upstream sources:
// - ncbi-blast-2.17.0+-src/c++/include/algo/blast/blastinput/blastx_args.hpp
// - ncbi-blast-2.17.0+-src/c++/src/algo/blast/blastinput/blastx_args.cpp

use super::blast_args::BlastAppArgs;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BlastxAppArgs {
    pub base: BlastAppArgs,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BlastxNodeArgs {
    pub base: BlastxAppArgs,
    pub output_stream: String,
    pub input_stream: Option<String>,
}
