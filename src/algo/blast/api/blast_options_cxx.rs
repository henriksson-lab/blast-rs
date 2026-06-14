/// NCBI C++: `CBlastOptions::EAPILocality` (`blast_options.hpp`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum EApiLocality {
    ELocal,
    ERemote,
    EBoth,
}

/// NCBI C++: `CBlastOptions` (`blast_options.hpp`).
#[derive(Clone, Debug)]
pub struct CBlastOptions {
    pub m_Local: Option<Box<CBlastOptionsLocal>>,
    pub m_Remote: Option<Box<CBlastOptionsRemote>>,
    pub m_ProgramName: String,
    pub m_ServiceName: String,
    pub m_DefaultsMode: bool,
    pub m_GenCodeSingletonVar: CAutomaticGenCodeSingleton,
}

#[derive(Clone, Debug)]
pub struct CBlastOptionsLocal {
    pub m_QueryOpts: CQuerySetUpOptions,
    pub m_LutOpts: CLookupTableOptions,
    pub m_InitWordOpts: CBlastInitialWordOptions,
    pub m_ExtnOpts: CBlastExtensionOptions,
    pub m_HitSaveOpts: CBlastHitSavingOptions,
    pub m_PSIBlastOpts: CPSIBlastOptions,
    pub m_DeltaBlastOpts: CPSIBlastOptions,
    pub m_DbOpts: CBlastDatabaseOptions,
    pub m_ScoringOpts: CBlastScoringOptions,
    pub m_EffLenOpts: CBlastEffectiveLengthsOptions,
    pub m_Program: EProgram,
    pub m_UseMBIndex: bool,
    pub m_ForceMBIndex: bool,
    pub m_OldStyleMBIndex: bool,
    pub m_MBIndexLoaded: bool,
    pub m_MBIndexName: String,
}

#[derive(Clone, Debug)]
pub struct CBlastOptionsRemote {
    pub m_Blast4Opts: Option<CBlast4Parameters>,
}

#[derive(Clone, Debug)]
pub struct CAutomaticGenCodeSingleton;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum EProgram {
    EBlastn,
    EBlastp,
    EBlastx,
    ETblastn,
    ETblastx,
    ERpsBlast,
    ERpsTblastn,
    EPhiBlastn,
    EPhiBlastp,
    EPsiBlast,
    EPsiTblastn,
    EDiscMegablast,
    EMagicBlast,
    EMapper,
    EUnknown,
}

#[derive(Clone, Debug)]
pub struct CQuerySetUpOptions {
    pub filter_string: Option<String>,
    pub strand_option: i32,
    pub genetic_code: i32,
    pub read_quality_filtering: bool,
}

#[derive(Clone, Debug)]
pub struct CLookupTableOptions {
    pub threshold: f64,
    pub lut_type: i32,
    pub word_size: i32,
    pub stride: u32,
    pub mb_template_length: u8,
    pub mb_template_type: u8,
    pub db_filter: bool,
    pub max_db_word_count: u8,
    pub phi_pattern: Option<String>,
}

#[derive(Clone, Debug)]
pub struct CBlastInitialWordOptions {
    pub window_size: i32,
    pub scan_range: i32,
    pub x_dropoff: f64,
    pub gap_trigger: f64,
}

#[derive(Clone, Debug)]
pub struct CBlastExtensionOptions {
    pub gap_x_dropoff: f64,
    pub gap_x_dropoff_final: f64,
    pub prelim_gap_ext: i32,
    pub traceback_ext: i32,
    pub composition_based_stats: i32,
    pub unified_p: i32,
    pub max_mismatches: i32,
    pub mismatch_window: i32,
    pub chaining: bool,
}

#[derive(Clone, Debug)]
pub struct CBlastHitSavingOptions {
    pub hitlist_size: i32,
    pub hsp_num_max: i32,
    pub max_hsps_per_subject: i32,
    pub culling_limit: i32,
    pub expect_value: f64,
    pub cutoff_score: i32,
    pub cutoff_score_fun: [i32; 2],
    pub percent_identity: f64,
    pub max_edit_distance: i32,
    pub query_cov_hsp_perc: f64,
    pub min_diag_separation: i32,
    pub do_sum_stats: bool,
    pub longest_intron: i32,
    pub mask_level: i32,
    pub low_score_perc: f64,
    pub paired: bool,
    pub splice: bool,
}

#[derive(Clone, Debug)]
pub struct CPSIBlastOptions {
    pub inclusion_ethresh: f64,
    pub pseudo_count: i32,
    pub nsg_compatibility_mode: bool,
}

#[derive(Clone, Debug)]
pub struct CBlastDatabaseOptions {
    pub genetic_code: i32,
}

#[derive(Clone, Debug)]
pub struct CBlastScoringOptions {
    pub program_number: u32,
    pub matrix: Option<String>,
    pub gapped_calculation: bool,
    pub complexity_adjusted_scoring: bool,
    pub reward: i32,
    pub penalty: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub shift_pen: i32,
    pub is_ooframe: bool,
}

#[derive(Clone, Debug)]
pub struct CBlastEffectiveLengthsOptions {
    pub db_length: i64,
    pub dbseq_num: i32,
    pub searchsp_eff: Vec<i64>,
}

#[derive(Clone, Debug)]
pub struct CBlast4Parameters {
    pub values: Vec<String>,
}
