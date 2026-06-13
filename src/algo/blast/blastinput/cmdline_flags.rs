// Upstream sources:
// - ncbi-blast-2.17.0+-src/c++/include/algo/blast/blastinput/cmdline_flags.hpp
// - ncbi-blast-2.17.0+-src/c++/src/algo/blast/blastinput/cmdline_flags.cpp

pub const K_ARG_QUERY: &str = "query";
pub const K_DFLT_ARG_QUERY: &str = "-";
pub const K_ARG_OUTPUT: &str = "out";

pub const K_ARG_DB: &str = "db";
pub const K_ARG_SUBJECT: &str = "subject";

pub const K_ARG_DB_SIZE: &str = "dbsize";
pub const K_ARG_DB_TYPE: &str = "dbtype";

pub const K_ARG_GI_LIST: &str = "gilist";
pub const K_ARG_DB_TITLE: &str = "title";
pub const K_ARG_SEQ_ID_LIST: &str = "seqidlist";
pub const K_ARG_NEGATIVE_GI_LIST: &str = "negative_gilist";
pub const K_ARG_NEGATIVE_SEQID_LIST: &str = "negative_seqidlist";
pub const K_ARG_DB_SOFT_MASK: &str = "db_soft_mask";
pub const K_ARG_DB_HARD_MASK: &str = "db_hard_mask";

pub const K_ARG_IPG_LIST: &str = "ipglist";
pub const K_ARG_NEGATIVE_IPG_LIST: &str = "negative_ipglist";

pub const K_TASK: &str = "task";

pub const K_ARG_QUERY_GENETIC_CODE: &str = "query_gencode";
pub const K_ARG_DB_GENETIC_CODE: &str = "db_gencode";

pub const K_ARG_REMOTE: &str = "remote";
pub const K_ARG_NUM_THREADS: &str = "num_threads";
pub const K_DFLT_IG_BLAST_NUM_THREADS: usize = 4;

pub const K_ARG_MATRIX_NAME: &str = "matrix";

pub const K_ARG_COMPLEXITY_ADJ: &str = "complexity_adjust";
pub const K_ARG_MASK_LEVEL: &str = "mask_level";
pub const K_DFLT_ARG_MASK_LEVEL: &str = "-1";

pub const K_ARG_EVALUE: &str = "evalue";
pub const K_ARG_MIN_RAW_GAPPED_SCORE: &str = "min_raw_gapped_score";

pub const K_ARG_MAX_TARGET_SEQUENCES: &str = "max_target_seqs";
pub const K_DFLT_ARG_MAX_TARGET_SEQUENCES: u32 = 100;

pub const K_ARG_GAP_OPEN: &str = "gapopen";
pub const K_ARG_GAP_EXTEND: &str = "gapextend";

pub const K_ARG_MISMATCH: &str = "penalty";
pub const K_ARG_MATCH: &str = "reward";

pub const K_ARG_UNGAPPED_X_DROPOFF: &str = "xdrop_ungap";
pub const K_ARG_GAPPED_X_DROPOFF: &str = "xdrop_gap";
pub const K_ARG_FINAL_GAPPED_X_DROPOFF: &str = "xdrop_gap_final";

pub const K_ARG_WINDOW_SIZE: &str = "window_size";
pub const K_ARG_OFF_DIAGONAL_RANGE: &str = "off_diagonal_range";
pub const K_DFLT_OFF_DIAGONAL_RANGE: i32 = 0;

pub const K_ARG_WORD_SIZE: &str = "word_size";
pub const K_ARG_WORD_SCORE_THRESHOLD: &str = "threshold";
pub const K_ARG_EFF_SEARCH_SPACE: &str = "searchsp";
pub const K_ARG_USE_SW_TRACEBACK: &str = "use_sw_tback";

pub const K_ARG_USE_LCASE_MASKING: &str = "lcase_masking";
pub const K_DFLT_ARG_USE_LCASE_MASKING: bool = false;
pub const K_ARG_STRAND: &str = "strand";
pub const K_DFLT_ARG_STRAND: &str = "both";
pub const K_ARG_QUERY_LOCATION: &str = "query_loc";
pub const K_ARG_SUBJECT_LOCATION: &str = "subject_loc";
pub const K_ARG_PARSE_DEFLINES: &str = "parse_deflines";
pub const K_DFLT_ARG_PARSE_DEFLINES: bool = false;

pub const K_ARG_MAX_INTRON_LENGTH: &str = "max_intron_length";
pub const K_DFLT_ARG_MAX_INTRON_LENGTH: i32 = 0;

pub const K_ARG_CULLING_LIMIT: &str = "culling_limit";
pub const K_DFLT_ARG_CULLING_LIMIT: i32 = 0;

pub const K_ARG_BEST_HIT_OVERHANG: &str = "best_hit_overhang";
pub const K_DFLT_ARG_BEST_HIT_OVERHANG: f64 = 0.1;
pub const K_ARG_BEST_HIT_SCORE_EDGE: &str = "best_hit_score_edge";
pub const K_DFLT_ARG_BEST_HIT_SCORE_EDGE: f64 = 0.1;

pub const K_ARG_SUBJECT_BEST_HIT: &str = "subject_besthit";
pub const K_ARG_FRAME_SHIFT_PENALTY: &str = "frame_shift_penalty";
pub const K_ARG_GAP_TRIGGER: &str = "gap_trigger";
pub const K_ARG_UNGAPPED: &str = "ungapped";

pub const K_ARG_COMP_BASED_STATS: &str = "comp_based_stats";
pub const K_DFLT_ARG_COMP_BASED_STATS: &str = "2";
pub const K_DFLT_ARG_COMP_BASED_STATS_DELTA: &str = "1";
pub const K_DFLT_ARG_COMP_BASED_STATS_RPS: &str = "1";

pub const K_DFLT_ARG_NO_FILTERING: &str = "no";
pub const K_DFLT_ARG_APPLY_FILTERING: &str = "yes";

pub const K_ARG_SEG_FILTERING: &str = "seg";
pub const K_DFLT_ARG_SEG_FILTERING: &str = "12 2.2 2.5";

pub const K_ARG_DUST_FILTERING: &str = "dust";
pub const K_DFLT_ARG_DUST_FILTERING: &str = "20 64 1";

pub const K_ARG_FILTERING_DB: &str = "filtering_db";
pub const K_ARG_WINDOW_MASKER_TAX_ID: &str = "window_masker_taxid";
pub const K_ARG_WINDOW_MASKER_DATABASE: &str = "window_masker_db";
pub const K_ARG_LOOKUP_TABLE_MASKING_ONLY: &str = "soft_masking";
pub const K_DFLT_ARG_LOOKUP_TABLE_MASKING_ONLY_PROT: &str = "false";
pub const K_DFLT_ARG_LOOKUP_TABLE_MASKING_ONLY_NUCL: &str = "true";

pub const K_ARG_PSI_NUM_ITERATIONS: &str = "num_iterations";
pub const K_DFLT_ARG_PSI_NUM_ITERATIONS: u32 = 1;
pub const K_ARG_PSI_INPUT_CHKPNT_FILE: &str = "in_pssm";
pub const K_ARG_PSI_OUTPUT_CHKPNT_FILE: &str = "out_pssm";
pub const K_ARG_ASCII_PSSM_OUTPUT_FILE: &str = "out_ascii_pssm";
pub const K_ARG_MSA_INPUT_FILE: &str = "in_msa";
pub const K_ARG_IGNORE_MSA_MASTER: &str = "ignore_msa_master";
pub const K_ARG_MSA_MASTER_INDEX: &str = "msa_master_idx";
pub const K_ARG_PSI_PSEUDOCOUNT: &str = "pseudocount";
pub const K_ARG_PSI_INCLUSION_E_THRESHOLD: &str = "inclusion_ethresh";
pub const K_ARG_SAVE_LAST_PSSM: &str = "save_pssm_after_last_round";
pub const K_ARG_SAVE_ALL_PSSMS: &str = "save_each_pssm";
pub const K_ARG_PHI_PATTERN_FILE: &str = "phi_pattern";

pub const K_ARG_GL_SUBJECT: &str = "germline_sequence_";
pub const K_ARG_GL_DATABASE: &str = "germline_db_";
pub const K_ARG_C_REGION_DATABASE: &str = "c_region_db";
pub const K_ARG_GL_NUM_ALIGN: &str = "num_alignments_";
pub const K_ARG_C_REGION_NUM_ALIGN: &str = "num_alignments_C";
pub const K_ARG_GL_CHAIN_TYPE: &str = "auxiliary_data";
pub const K_ARG_CUSTOM_INTERNAL_DATA: &str = "custom_internal_data";
pub const K_ARG_D_FRAME_DEFINITION_FILE: &str = "d_frame_data";
pub const K_ARG_GL_ORIGIN: &str = "organism";
pub const K_ARG_GL_DOMAIN_SYSTEM: &str = "domain_system";
pub const K_ARG_GL_FOCUS_V: &str = "focus_on_V_segment";
pub const K_ARG_EXTEND_ALIGN_5END: &str = "extend_align5end";
pub const K_ARG_EXTEND_ALIGN_3END: &str = "extend_align3end";
pub const K_ARG_DETECT_OVERLAP: &str = "allow_vdj_overlap";
pub const K_ARG_MIN_V_LENGTH: &str = "min_V_length";
pub const K_ARG_MIN_J_LENGTH: &str = "min_J_length";
pub const K_ARG_NUM_CLONOTYPE: &str = "num_clonotype";
pub const K_ARG_CLONOTYPE_FILE: &str = "clonotype_out";
pub const K_ARG_TRANSLATE: &str = "show_translation";
pub const K_ARG_MIN_D_MATCH: &str = "min_D_match";
pub const K_ARG_D_PENALTY: &str = "D_penalty";
pub const K_ARG_V_PENALTY: &str = "V_penalty";
pub const K_ARG_J_PENALTY: &str = "J_penalty";
pub const K_ARG_IG_SEQ_TYPE: &str = "ig_seqtype";
pub const K_ARG_MAX_HSPS_PER_SUBJECT: &str = "max_hsps";
pub const K_ARG_SUM_STATS: &str = "sum_stats";
pub const K_ARG_PERCENT_IDENTITY: &str = "perc_identity";
pub const K_ARG_NO_GREEDY_EXTENSION: &str = "no_greedy";
pub const K_ARG_DMB_TEMPLATE_TYPE: &str = "template_type";
pub const K_ARG_DMB_TEMPLATE_LENGTH: &str = "template_length";
pub const K_ARG_NO_DISCORDANT: &str = "no_discordant";
pub const K_ARG_FWD_REV: &str = "fr";
pub const K_ARG_REV_FWD: &str = "rf";
pub const K_ARG_FWD_ONLY: &str = "f";
pub const K_ARG_REV_ONLY: &str = "r";
pub const K_ARG_ONLY_STRAND_SPECIFIC: &str = "only_strand_specific";

pub const K_ARG_INPUT_SEARCH_STRATEGY: &str = "import_search_strategy";
pub const K_ARG_OUTPUT_SEARCH_STRATEGY: &str = "export_search_strategy";

pub const K_ARG_USE_INDEX: &str = "use_index";
pub const K_ARG_INDEX_NAME: &str = "index_name";
pub const K_ARG_OLD_STYLE_INDEX: &str = "old_style_index";
pub const K_DFLT_ARG_OLD_STYLE_INDEX: &str = "true";
pub const K_DFLT_ARG_USE_INDEX: bool = false;

pub const K_ARG_ENTREZ_QUERY: &str = "entrez_query";

pub const K_ARG_RPS_DB: &str = "rpsdb";
pub const K_DFLT_ARG_RPS_DB: &str = "cdd_delta";
pub const K_ARG_DOMAIN_INCLUSION_E_THRESHOLD: &str = "domain_inclusion_ethresh";
pub const K_ARG_SHOW_DOMAIN_HITS: &str = "show_domain_hits";

pub const K_ARG_J_DISTANCE: &str = "thresh";
pub const K_DFLT_ARG_J_DISTANCE: &str = "0.1";
pub const K_ARG_MIN_HITS: &str = "min_hits";
pub const K_DFLT_ARG_MIN_HITS: &str = "0";
pub const K_ARG_CANDIDATE_SEQS: &str = "candidates";
pub const K_DFLT_ARG_CANDIDATE_SEQS: &str = "1000";

pub const K_ARG_RID: &str = "rid";
pub const K_ARG_ARCHIVE: &str = "archive";

pub const K_ARG_QUERY_COV_HSP_PERC: &str = "qcov_hsp_perc";
pub const K_ARG_LINE_LENGTH: &str = "line_length";
pub const K_ARG_TAX_ID_LIST: &str = "taxids";
pub const K_ARG_TAX_ID_LIST_FILE: &str = "taxidlist";
pub const K_ARG_NEGATIVE_TAX_ID_LIST: &str = "negative_taxids";
pub const K_ARG_NEGATIVE_TAX_ID_LIST_FILE: &str = "negative_taxidlist";
pub const K_ARG_NO_TAX_ID_EXPANSION: &str = "no_taxid_expansion";

pub const K_ARG_PAIRED: &str = "paired";
pub const K_ARG_SCORE: &str = "score";
pub const K_ARG_LIMIT_LOOKUP: &str = "limit_lookup";
pub const K_ARG_SPLICE: &str = "splice";
pub const K_ARG_LOOKUP_STRIDE: &str = "lookup_stride";
pub const K_ARG_INPUT_FORMAT: &str = "infmt";
pub const K_ARG_QUALITY_FILTER: &str = "validate_seqs";
pub const K_ARG_QUERY_MATE: &str = "query_mate";
pub const K_ARG_REF_TYPE: &str = "reftype";
pub const K_ARG_OUTPUT_GZIP: &str = "gzo";
pub const K_ARG_SRA_ACCESSION: &str = "sra";
pub const K_ARG_SRA_ACCESSION_BATCH: &str = "sra_batch";
pub const K_ARG_NO_READ_ID_TRIM: &str = "no_query_id_trim";
pub const K_ARG_NO_UNALIGNED: &str = "no_unaligned";
pub const K_ARG_ENABLE_SRA_CACHE: &str = "sra_cache";
pub const K_ARG_MAX_EDIT_DIST: &str = "max_edit_dist";
pub const K_ARG_MAX_DB_WORD_COUNT: &str = "max_db_word_count";
pub const K_ARG_PRINT_MD_TAG: &str = "md_tag";
pub const K_ARG_UNALIGNED_OUTPUT: &str = "out_unaligned";
pub const K_ARG_UNALIGNED_FORMAT: &str = "unaligned_fmt";
pub const K_ARG_USER_TAG: &str = "tag";

pub const K_ARG_MT_MODE: &str = "mt_mode";
