// Upstream sources:
// - ncbi-blast-2.17.0+-src/c++/include/algo/blast/format/blast_format.hpp
// - ncbi-blast-2.17.0+-src/c++/src/algo/blast/format/blast_format.cpp

use crate::algo::blast::blastinput::blast_args::OutputFormat;
use crate::algo::blast::blastinput::blast_input::SeqRange;
use crate::algo::blast::format::blastfmtutil::{BlastFormatUtilDbInfo, BlastXmlIncremental};

#[derive(Debug, Clone, PartialEq)]
pub struct BlastFormatClone {
    pub na: String,
    pub chain_type: String,
    pub aa: String,
    pub v_gene: String,
    pub d_gene: String,
    pub j_gene: String,
    pub c_gene: String,
    pub seqid: String,
    pub identity: f64,
    pub productive: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DisplayOption {
    Descriptions,
    Alignments,
    Metadata,
    DescriptionsWithTemplates,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BlastFormat {
    pub format_type: OutputFormat,
    pub is_html: bool,
    pub db_is_aa: bool,
    pub believe_query: bool,
    pub outfile: String,
    pub num_summary: i32,
    pub num_alignments: i32,
    pub hitlist_size: i32,
    pub program: String,
    pub db_name: String,
    pub query_gen_code: i32,
    pub db_gen_code: i32,
    pub show_gi: bool,
    pub show_linked_set_size: bool,
    pub is_ungapped_search: bool,
    pub matrix_name: Option<String>,
    pub scoring_matrix: Vec<Vec<i32>>,
    pub scope: Option<String>,
    pub is_bl2seq: bool,
    pub is_db_scan: bool,
    pub subject_tag: String,
    pub is_remote_search: bool,
    pub queries_formatted: i32,
    pub megablast: bool,
    pub indexed_megablast: bool,
    pub seq_info_src: Option<String>,
    pub db_info: Vec<BlastFormatUtilDbInfo>,
    pub search_db: Option<String>,
    pub accumulated_queries: Option<String>,
    pub accumulated_results: Option<String>,
    pub disable_ka_stats: bool,
    pub custom_output_format_spec: String,
    pub blast_xml_incremental: Option<BlastXmlIncremental>,
    pub ig_options: Option<String>,
    pub domain_db_info: Vec<BlastFormatUtilDbInfo>,
    pub options: Option<String>,
    pub is_vdb: bool,
    pub query_range: SeqRange,
    pub is_iterative: bool,
    pub base_file: String,
    pub xml_file_count: i32,
    pub line_length: usize,
    pub original_exception_mask: u32,
    pub sam_formatter: Option<String>,
    pub cmdline: String,
    pub long_seq_id: bool,
    pub defline_templates: Option<String>,
    pub align_templates: Option<String>,
    pub align_seq_list: String,
    pub hits_sort_option: i32,
    pub hsps_sort_option: i32,
    pub custom_delimiter: String,
}
