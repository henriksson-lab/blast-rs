// Upstream sources:
// - ncbi-blast-2.17.0+-src/c++/include/algo/blast/blastinput/rpsblast_args.hpp
// - ncbi-blast-2.17.0+-src/c++/src/algo/blast/blastinput/rpsblast_args.cpp

use super::blast_args::{BlastAppArgs, MtArgs};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RpsBlastMtArgs {
    pub base: MtArgs,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RpsBlastAppArgs {
    pub base: BlastAppArgs,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RpsBlastNodeArgs {
    pub base: RpsBlastAppArgs,
    pub output_stream: String,
    pub input_stream: Option<String>,
}
