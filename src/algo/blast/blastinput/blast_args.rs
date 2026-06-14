// Upstream sources:
// - ncbi-blast-2.17.0+-src/c++/include/algo/blast/blastinput/blast_args.hpp
// - ncbi-blast-2.17.0+-src/c++/src/algo/blast/blastinput/blast_args.cpp

use super::blast_input::SeqRange;
use std::collections::BTreeSet;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BlastCmdLineArgs {
    pub arg_descriptions: Vec<String>,
    pub parsed_args: Vec<String>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StdCmdLineArgs {
    pub input_stream: Option<String>,
    pub output_stream: Option<String>,
    pub decompress_input_stream: Option<String>,
    pub compress_output_stream: Option<String>,
    pub query_tmp_input_file: Option<String>,
    pub gzip_enabled: bool,
    pub sra_accession_enabled: bool,
    pub unaligned_output_stream: Option<String>,
    pub unaligned_compress_output_stream: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProgramDescriptionArgs {
    pub program_name: String,
    pub program_description: String,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TaskCmdLineArgs {
    pub supported_tasks: BTreeSet<String>,
    pub default_task: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GenericSearchArgs {
    pub query_is_protein: bool,
    pub is_rps_blast: bool,
    pub show_percent_identity: bool,
    pub is_tblastx: bool,
    pub is_ig_blast: bool,
    pub suppress_sum_stats: bool,
    pub is_blastn: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FilteringArgs {
    pub query_is_protein: bool,
    pub filter_by_default: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CompositionBasedStatsArgs {
    pub is_two_and_three_supported: bool,
    pub default_option: String,
    pub zero_option_description: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GeneticCodeTarget {
    Query,
    Database,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GeneticCodeArgs {
    pub target: GeneticCodeTarget,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GapTriggerArgs {
    pub query_is_protein: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PssmEngineArgs {
    pub is_delta_blast: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PsiTargetDatabase {
    BlastDb,
    CddDb,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PsiBlastArgs {
    pub database_target: PsiTargetDatabase,
    pub num_iterations: usize,
    pub checkpoint_output: Option<String>,
    pub ascii_matrix_output: Option<String>,
    pub pssm: Option<String>,
    pub is_delta_blast: bool,
    pub save_last_pssm: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub struct KBlastpArgs {
    pub jaccard_distance: f64,
    pub min_hits: i32,
    pub db_index: String,
    pub candidate_seqs: i32,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DeltaBlastArgs {
    pub domain_db: Option<String>,
    pub show_domain_hits: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct QueryOptionsArgs {
    pub strand: super::blast_input::NaStrand,
    pub range: SeqRange,
    pub use_lowercase_mask: bool,
    pub parse_deflines: bool,
    pub query_cannot_be_nucl: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MapperInputFormat {
    Fasta,
    Fastc,
    Fastq,
    Sra,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MapperQueryOptionsArgs {
    pub base: QueryOptionsArgs,
    pub is_paired: bool,
    pub input_format: MapperInputFormat,
    pub sra_accessions: Vec<String>,
    pub mate_input_stream: Option<String>,
    pub decompress_input_stream: Option<String>,
    pub enable_sra_cache: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BlastDatabaseArgs {
    pub search_db: Option<String>,
    pub request_molecule_type: bool,
    pub is_rps_blast: bool,
    pub is_ig_blast: bool,
    pub is_protein: bool,
    pub is_mapper: bool,
    pub is_kblast: bool,
    pub subjects: Option<String>,
    pub scope: Option<String>,
    pub supports_database_masking: bool,
    pub support_ipg_filtering: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IgBlastArgs {
    pub is_protein: bool,
    pub ig_options: Option<String>,
    pub scope: Option<String>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutputFormat {
    Pairwise,
    QueryAnchored,
    QueryAnchoredNoIdentities,
    FlatQueryAnchored,
    FlatQueryAnchoredNoIdentities,
    Xml,
    Tabular,
    TabularWithComments,
    TextASN1,
    BinaryASN1,
    Csv,
    BlastArchive,
    SeqalignJson,
    MultipleFileXml2,
    MultipleFileJson,
    Sam,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FormatFlag {
    Default,
    Archive,
    AlignmentOnly,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FormattingArgs {
    pub output_format: OutputFormat,
    pub show_gis: bool,
    pub num_descriptions: u64,
    pub num_alignments: u64,
    pub default_num_descriptions: u64,
    pub default_num_alignments: u64,
    pub html: bool,
    pub is_ig_blast: bool,
    pub custom_output_format_spec: String,
    pub line_length: usize,
    pub format_flags: FormatFlag,
    pub hits_sort_option: i32,
    pub hsps_sort_option: i32,
    pub custom_delimiter: String,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MapperFormattingArgs {
    pub base: FormattingArgs,
    pub trim_read_ids: bool,
    pub print_unaligned: bool,
    pub no_discordant: bool,
    pub fwd_rev: bool,
    pub rev_fwd: bool,
    pub fwd_only: bool,
    pub rev_only: bool,
    pub only_strand_specific: bool,
    pub print_md_tag: bool,
    pub unaligned_output_format: OutputFormat,
    pub user_tag: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MtMode {
    NotSupported,
    Automatic,
    Forced,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MtArgs {
    pub num_threads: usize,
    pub mt_mode: MtMode,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RemoteArgs {
    pub is_remote: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DebugArgs {
    pub debug_output: bool,
    pub remote_debug_output: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BlastAppArgs {
    pub args: Vec<BlastCmdLineArgs>,
    pub query_options_args: Option<QueryOptionsArgs>,
    pub blast_database_args: Option<BlastDatabaseArgs>,
    pub formatting_args: Option<FormattingArgs>,
    pub mt_args: Option<MtArgs>,
    pub remote_args: Option<RemoteArgs>,
    pub std_cmd_line_args: Option<StdCmdLineArgs>,
    pub search_strategy_args: Option<String>,
    pub debug_args: Option<DebugArgs>,
    pub hsp_filtering_args: Option<String>,
    pub options_handle: Option<String>,
    pub task: String,
    pub client_id: String,
    pub is_ungapped: bool,
}
