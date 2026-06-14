// Upstream sources:
// - ncbi-blast-2.17.0+-src/c++/include/algo/blast/blastinput/rpstblastn_args.hpp
// - ncbi-blast-2.17.0+-src/c++/src/algo/blast/blastinput/rpstblastn_args.cpp

use super::blast_args::BlastAppArgs;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RpstBlastnAppArgs {
    pub base: BlastAppArgs,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RpstBlastnNodeArgs {
    pub base: RpstBlastnAppArgs,
    pub output_stream: String,
    pub input_stream: Option<String>,
}
