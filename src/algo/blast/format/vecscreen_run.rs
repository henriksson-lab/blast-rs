// Upstream sources:
// - ncbi-blast-2.17.0+-src/c++/include/algo/blast/format/vecscreen_run.hpp
// - ncbi-blast-2.17.0+-src/c++/src/algo/blast/format/vecscreen_run.cpp

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VecScreenVersion {
    pub version: String,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VecscreenSummary {
    pub seqid: Option<String>,
    pub range_from: u64,
    pub range_to: u64,
    pub match_type: String,
    pub aligns: Vec<String>,
    pub drops: Vec<String>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VecscreenOutputFormat {
    Text,
    Html,
    Tabular,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VecscreenFormatter {
    pub screener: Option<String>,
    pub scope: Option<String>,
    pub outfmt: VecscreenOutputFormat,
    pub html_output: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VecscreenRun {
    pub seq_loc: Option<String>,
    pub scope: Option<String>,
    pub db: String,
    pub terminal_flexibility: u64,
    pub vecscreen: Option<String>,
    pub queries: Option<String>,
    pub seqalign_set: Option<String>,
    pub raw_blast_results: Option<String>,
}
