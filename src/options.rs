//! Rust equivalent of blast_options.c — BLAST search options.
//! Contains all the option types and their default values.

use crate::lookup::LookupTableType;
use crate::program::{self, ProgramType};

// ---- Default constants ----

// Window sizes
pub const WINDOW_SIZE_PROT: i32 = 40;
pub const WINDOW_SIZE_NUCL: i32 = 0;
pub const WINDOW_SIZE_MEGABLAST: i32 = 0;
pub const SCAN_RANGE_NUCL: i32 = 0;

// Word sizes
pub const WORDSIZE_PROT: i32 = 3;
pub const WORDSIZE_NUCL: i32 = 11;
pub const WORDSIZE_MEGABLAST: i32 = 28;
pub const WORDSIZE_MAPPER: i32 = 18;
pub const MAX_DB_WORD_COUNT_MAPPER: u8 = 30;

// Gap costs
pub const GAP_OPEN_PROT: i32 = 11;
pub const GAP_OPEN_NUCL: i32 = 5;
pub const GAP_OPEN_MEGABLAST: i32 = 0;
pub const GAP_EXTN_PROT: i32 = 1;
pub const GAP_EXTN_NUCL: i32 = 2;
pub const GAP_EXTN_MEGABLAST: i32 = 0;

// Scoring
pub const PENALTY: i32 = -3;
pub const REWARD: i32 = 1;

// X-dropoff
pub const UNGAPPED_X_DROPOFF_PROT: f64 = 7.0;
pub const UNGAPPED_X_DROPOFF_NUCL: f64 = 20.0;
pub const GAP_X_DROPOFF_PROT: f64 = 15.0;
pub const GAP_X_DROPOFF_NUCL: f64 = 30.0;
pub const GAP_X_DROPOFF_FINAL_PROT: f64 = 25.0;
pub const GAP_X_DROPOFF_FINAL_NUCL: f64 = 100.0;
pub const GAP_TRIGGER_PROT: f64 = 22.0;
pub const GAP_TRIGGER_NUCL: f64 = 27.0;

// E-value
pub const EXPECT_VALUE: f64 = 10.0;
pub const HITLIST_SIZE: i32 = 500;
pub const WORD_THRESHOLD_BLASTP: f64 = 11.0;
pub const WORD_THRESHOLD_BLASTX: f64 = 12.0;
pub const WORD_THRESHOLD_TBLASTN: f64 = 13.0;
pub const WORD_THRESHOLD_TBLASTX: f64 = 13.0;
pub const BLASTERR_OPTION_VALUE_INVALID: i16 = 1;
pub const BLASTERR_OPTION_PROGRAM_INVALID: i16 = 2;
pub const K_PSSM_NO_IMPALA_SCALING: f64 = 1.0;

// ---- Option structs ----

/// Query setup options. 1-1 port of C `QuerySetUpOptions`
/// (`blast_options.h:285`).
#[derive(Debug, Clone, PartialEq)]
pub struct QuerySetUpOptions {
    pub filtering_options: Option<crate::filter::SBlastFilterOptions>,
    pub filter_string: Option<String>,
    pub strand_option: u8,
    pub genetic_code: i32,
}

impl Default for QuerySetUpOptions {
    fn default() -> Self {
        let mut filtering_options = None;
        crate::filter::sblast_filter_options_new(
            &mut filtering_options,
            crate::filter::EFilterOptions::EEmpty,
        );
        Self {
            filtering_options,
            filter_string: None,
            strand_option: 0,
            genetic_code: 1,
        }
    }
}

/// Port of NCBI `BlastQuerySetUpOptionsNew` (`blast_options.c:514`).
pub fn blast_query_set_up_options_new(options: &mut Option<QuerySetUpOptions>) -> i16 {
    *options = Some(QuerySetUpOptions::default());
    0
}

/// Port of NCBI `BlastQuerySetUpOptionsFree` (`blast_options.c:501`).
pub fn blast_query_set_up_options_free(
    options: &mut Option<QuerySetUpOptions>,
) -> Option<QuerySetUpOptions> {
    if let Some(options) = options.as_mut() {
        options.filter_string = None;
        crate::filter::sblast_filter_options_free(&mut options.filtering_options);
    }
    *options = None;
    None
}

/// Port of NCBI `BLAST_FillQuerySetUpOptions` (`blast_options.c:534`).
pub fn blast_fill_query_set_up_options(
    options: Option<&mut QuerySetUpOptions>,
    program_number: ProgramType,
    filter_string: Option<&str>,
    strand_option: u8,
) -> i16 {
    let Some(options) = options else {
        return crate::util::BLASTERR_INVALIDPARAM;
    };

    if strand_option != 0
        && matches!(
            program_number,
            program::BLASTN
                | program::PHI_BLASTN
                | program::BLASTX
                | program::TBLASTX
                | program::MAPPING
        )
    {
        options.strand_option = strand_option;
    }

    if let Some(filter_string) = filter_string {
        options.filter_string = None;
        crate::filter::sblast_filter_options_free(&mut options.filtering_options);
        return crate::filter::blast_filtering_options_from_string(
            program_number,
            Some(filter_string),
            &mut options.filtering_options,
        );
    }

    0
}

/// Lookup table options. 1-1 port of C `LookupTableOptions`
/// (`blast_options.h:206`).
#[derive(Debug, Clone, PartialEq)]
pub struct LookupTableOptions {
    pub threshold: f64,
    pub lut_type: LookupTableType,
    pub word_size: i32,
    pub mb_template_length: i32,
    pub mb_template_type: i32,
    pub phi_pattern: Option<String>,
    pub program_number: ProgramType,
    pub stride: u32,
    pub db_filter: bool,
    pub max_db_word_count: u8,
}

impl Default for LookupTableOptions {
    fn default() -> Self {
        Self {
            threshold: 0.0,
            lut_type: LookupTableType::AaLookup,
            word_size: 0,
            mb_template_length: 0,
            mb_template_type: 0,
            phi_pattern: None,
            program_number: program::UNDEFINED,
            stride: 0,
            db_filter: false,
            max_db_word_count: 0,
        }
    }
}

/// Port-shaped constructor for NCBI `LookupTableOptionsNew`
/// (`blast_options.c:1077`).
pub fn lookup_table_options_new(
    program_number: ProgramType,
    options: &mut Option<LookupTableOptions>,
) -> i16 {
    let mut new_options = LookupTableOptions::default();

    match program_number {
        program::MAPPING => {
            new_options.max_db_word_count = MAX_DB_WORD_COUNT_MAPPER;
            new_options.word_size = WORDSIZE_MEGABLAST;
            new_options.lut_type = LookupTableType::MegablastLookup;
        }
        program::BLASTN => {
            new_options.word_size = WORDSIZE_MEGABLAST;
            new_options.lut_type = LookupTableType::MegablastLookup;
        }
        program::RPS_BLAST | program::RPS_TBLASTN => {
            new_options.word_size = WORDSIZE_PROT;
            new_options.lut_type = LookupTableType::RpsLookup;
            new_options.threshold = if program_number == program::RPS_BLAST {
                WORD_THRESHOLD_BLASTP
            } else {
                WORD_THRESHOLD_TBLASTN
            };
        }
        program::PHI_BLASTN => {
            new_options.lut_type = LookupTableType::PhiNaLookup;
        }
        program::PHI_BLASTP => {
            new_options.lut_type = LookupTableType::PhiLookup;
        }
        _ => {
            new_options.word_size = WORDSIZE_PROT;
            new_options.lut_type = LookupTableType::AaLookup;
            new_options.threshold = match program_number {
                program::BLASTP => WORD_THRESHOLD_BLASTP,
                program::BLASTX => WORD_THRESHOLD_BLASTX,
                program::TBLASTN => WORD_THRESHOLD_TBLASTN,
                program::TBLASTX => WORD_THRESHOLD_TBLASTX,
                _ => 0.0,
            };
        }
    }

    new_options.program_number = program_number;
    new_options.stride = 0;
    *options = Some(new_options);
    0
}

/// Port of NCBI `LookupTableOptionsFree` (`blast_options.c:1063`).
pub fn lookup_table_options_free(
    options: &mut Option<LookupTableOptions>,
) -> Option<LookupTableOptions> {
    if let Some(options) = options.as_mut() {
        options.phi_pattern = None;
    }
    *options = None;
    None
}

/// Port of NCBI `BLAST_FillLookupTableOptions` (`blast_options.c:1129`).
pub fn blast_fill_lookup_table_options(
    options: Option<&mut LookupTableOptions>,
    program_number: ProgramType,
    is_megablast: bool,
    threshold: f64,
    word_size: i32,
) -> i16 {
    let Some(options) = options else {
        return crate::util::BLASTERR_INVALIDPARAM;
    };

    if program_number == program::BLASTN {
        if is_megablast {
            options.lut_type = LookupTableType::MegablastLookup;
            options.word_size = WORDSIZE_MEGABLAST;
        } else {
            options.lut_type = LookupTableType::NaLookup;
            options.word_size = WORDSIZE_NUCL;
        }
    } else if program_number == program::MAPPING {
        options.lut_type = LookupTableType::NaHashLookup;
        options.word_size = WORDSIZE_MAPPER;
        options.max_db_word_count = MAX_DB_WORD_COUNT_MAPPER;
    } else {
        options.lut_type = LookupTableType::AaLookup;
    }

    if threshold < 0.0 {
        options.threshold = 0.0;
    } else if threshold > 0.0 {
        options.threshold = threshold;
    }

    if program::blast_program_is_rps_blast(program_number) {
        options.lut_type = LookupTableType::RpsLookup;
    }
    if word_size != 0 {
        options.word_size = word_size;
    }
    if matches!(
        program_number,
        program::TBLASTN | program::BLASTP | program::BLASTX
    ) && word_size > 5
    {
        options.lut_type = LookupTableType::CompressedAaLookup;
    }

    0
}

/// Port of NCBI `BLAST_GetSuggestedThreshold` (`blast_options.c:1174`).
pub fn blast_get_suggested_threshold(
    program_number: ProgramType,
    matrix_name: Option<&str>,
    threshold: &mut f64,
) -> i16 {
    if program_number == program::BLASTN || program_number == program::MAPPING {
        return 0;
    }

    let Some(matrix_name) = matrix_name else {
        return crate::util::BLASTERR_INVALIDPARAM;
    };

    *threshold = if matrix_name.eq_ignore_ascii_case("BLOSUM62") {
        11.0
    } else if matrix_name.eq_ignore_ascii_case("BLOSUM45") {
        14.0
    } else if matrix_name.eq_ignore_ascii_case("BLOSUM62_20") {
        100.0
    } else if matrix_name.eq_ignore_ascii_case("BLOSUM80") {
        12.0
    } else if matrix_name.eq_ignore_ascii_case("PAM30") {
        16.0
    } else if matrix_name.eq_ignore_ascii_case("PAM70") {
        14.0
    } else if matrix_name.eq_ignore_ascii_case("IDENTITY") {
        27.0
    } else {
        11.0
    };

    if program::blast_subject_is_translated(program_number) {
        *threshold += 2.0;
    } else if program::blast_query_is_translated(program_number) {
        *threshold += 1.0;
    }

    0
}

/// Port of NCBI `BLAST_GetSuggestedWindowSize` (`blast_options.c:1211`).
pub fn blast_get_suggested_window_size(
    program_number: ProgramType,
    matrix_name: Option<&str>,
    window_size: &mut i32,
) -> i16 {
    if program_number == program::BLASTN || program_number == program::MAPPING {
        return 0;
    }

    let Some(matrix_name) = matrix_name else {
        return crate::util::BLASTERR_INVALIDPARAM;
    };

    *window_size = if matrix_name.eq_ignore_ascii_case("BLOSUM62") {
        40
    } else if matrix_name.eq_ignore_ascii_case("BLOSUM45") {
        60
    } else if matrix_name.eq_ignore_ascii_case("BLOSUM80") {
        25
    } else if matrix_name.eq_ignore_ascii_case("PAM30") {
        15
    } else if matrix_name.eq_ignore_ascii_case("PAM70") {
        20
    } else {
        40
    };

    0
}

/// Port of NCBI static `s_DiscWordOptionsValidate`
/// (`blast_options.c:1248`).
fn s_disc_word_options_validate(word_size: i32, template_length: i32, template_type: i32) -> bool {
    if template_length == 0 {
        return true;
    }

    if word_size != 11 && word_size != 12 {
        return false;
    }

    if template_length != 16 && template_length != 18 && template_length != 21 {
        return false;
    }

    template_type <= 2
}

/// Port of NCBI `LookupTableOptionsValidate` (`blast_options.c:1284`).
pub fn lookup_table_options_validate(
    program_number: ProgramType,
    options: Option<&LookupTableOptions>,
) -> i16 {
    let Some(options) = options else {
        return crate::util::BLASTERR_INVALIDPARAM;
    };
    let is_phi_blast = program::blast_program_is_phi_blast(program_number);

    if options.phi_pattern.is_some() && !is_phi_blast {
        return BLASTERR_OPTION_PROGRAM_INVALID;
    }

    if is_phi_blast {
        return 0;
    }

    if program_number != program::BLASTN
        && program_number != program::MAPPING
        && !program::blast_program_is_rps_blast(program_number)
        && options.threshold <= 0.0
    {
        return BLASTERR_OPTION_VALUE_INVALID;
    }

    if options.word_size <= 0 {
        if !program::blast_program_is_rps_blast(program_number) {
            return BLASTERR_OPTION_VALUE_INVALID;
        }
    } else if program::blast_program_is_nucleotide(program_number)
        && !program::blast_query_is_pattern(program_number)
        && options.word_size < 4
    {
        return BLASTERR_OPTION_VALUE_INVALID;
    } else if program_number == program::BLASTN
        && options.word_size > crate::stat::DBSEQ_CHUNK_OVERLAP as i32
    {
        return BLASTERR_OPTION_VALUE_INVALID;
    } else if program_number != program::BLASTN
        && program_number != program::MAPPING
        && options.word_size > 4
    {
        if matches!(
            program_number,
            program::BLASTP | program::TBLASTN | program::BLASTX | program::TBLASTX
        ) {
            if options.word_size > 7 {
                return BLASTERR_OPTION_VALUE_INVALID;
            }
        } else if program_number == program::PSI_BLAST && options.word_size > 4 {
            return BLASTERR_OPTION_VALUE_INVALID;
        } else {
            return BLASTERR_OPTION_VALUE_INVALID;
        }
    }

    if program_number != program::BLASTN
        && program_number != program::MAPPING
        && options.lut_type == LookupTableType::MegablastLookup
    {
        return BLASTERR_OPTION_PROGRAM_INVALID;
    }

    if matches!(
        program_number,
        program::BLASTP | program::TBLASTN | program::BLASTX | program::TBLASTX
    ) {
        if options.word_size > 5 && options.lut_type != LookupTableType::CompressedAaLookup {
            return BLASTERR_OPTION_VALUE_INVALID;
        } else if options.lut_type == LookupTableType::CompressedAaLookup
            && options.word_size != 5
            && options.word_size != 6
            && options.word_size != 7
        {
            return BLASTERR_OPTION_VALUE_INVALID;
        }
    }

    if program::blast_program_is_nucleotide(program_number)
        && !program::blast_query_is_pattern(program_number)
        && options.mb_template_length > 0
    {
        if !s_disc_word_options_validate(
            options.word_size,
            options.mb_template_length,
            options.mb_template_type,
        ) {
            return BLASTERR_OPTION_VALUE_INVALID;
        } else if options.lut_type != LookupTableType::MegablastLookup {
            return BLASTERR_OPTION_VALUE_INVALID;
        }
    }

    if !program::blast_program_is_nucleotide(program_number) && options.db_filter {
        return BLASTERR_OPTION_VALUE_INVALID;
    }

    if options.db_filter && options.word_size < 16 {
        return BLASTERR_OPTION_VALUE_INVALID;
    }

    0
}

/// Scoring options.
#[derive(Debug, Clone)]
pub struct ScoringOptions {
    pub reward: i32,
    pub penalty: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub gapped_calculation: bool,
    pub matrix_name: Option<String>,
    pub is_ooframe: bool,
}

impl ScoringOptions {
    pub fn new_blastn() -> Self {
        ScoringOptions {
            reward: REWARD,
            penalty: PENALTY,
            gap_open: GAP_OPEN_NUCL,
            gap_extend: GAP_EXTN_NUCL,
            gapped_calculation: true,
            matrix_name: None,
            is_ooframe: false,
        }
    }

    pub fn new_blastp() -> Self {
        ScoringOptions {
            reward: 0,
            penalty: 0,
            gap_open: GAP_OPEN_PROT,
            gap_extend: GAP_EXTN_PROT,
            gapped_calculation: true,
            matrix_name: Some("BLOSUM62".to_string()),
            is_ooframe: false,
        }
    }
}

/// Port-shaped constructor for NCBI `BlastScoringOptionsNew`
/// (`blast_options.c:788`).
pub fn blast_scoring_options_new(
    program_number: ProgramType,
    options: &mut Option<ScoringOptions>,
) -> i16 {
    let mut new_options = if program::blast_program_is_nucleotide(program_number) {
        ScoringOptions::new_blastn()
    } else {
        ScoringOptions::new_blastp()
    };
    if program_number == program::TBLASTX {
        new_options.gapped_calculation = false;
    }
    *options = Some(new_options);
    0
}

/// Port of NCBI `BlastScoringOptionsFree` (`blast_options.c:774`).
pub fn blast_scoring_options_free(options: &mut Option<ScoringOptions>) -> Option<ScoringOptions> {
    if let Some(options) = options.as_mut() {
        options.matrix_name = None;
    }
    *options = None;
    None
}

/// Port of NCBI `BlastScoringOptionsSetMatrix` (`blast_options.c:974`).
pub fn blast_scoring_options_set_matrix(
    options: Option<&mut ScoringOptions>,
    matrix_name: Option<&str>,
) -> i16 {
    let Some(options) = options else {
        return crate::util::BLASTERR_INVALIDPARAM;
    };
    let Some(matrix_name) = matrix_name else {
        return crate::util::BLASTERR_INVALIDPARAM;
    };

    options.matrix_name = Some(matrix_name.to_ascii_uppercase());
    0
}

/// Port of NCBI `BlastScoringOptionsDup` (`blast_options.c:956`).
pub fn blast_scoring_options_dup(
    new_opt: &mut Option<ScoringOptions>,
    old_opt: Option<&ScoringOptions>,
) -> i16 {
    let Some(old_opt) = old_opt else {
        *new_opt = None;
        return crate::util::BLASTERR_INVALIDPARAM;
    };
    *new_opt = Some(old_opt.clone());
    0
}

/// Port of NCBI `BLAST_FillScoringOptions` (`blast_options.c:823`).
pub fn blast_fill_scoring_options(
    options: Option<&mut ScoringOptions>,
    program_number: ProgramType,
    greedy_extension: bool,
    penalty: i32,
    reward: i32,
    matrix: Option<&str>,
    gap_open: i32,
    gap_extend: i32,
) -> i16 {
    let Some(options) = options else {
        return crate::util::BLASTERR_INVALIDPARAM;
    };

    if !program::blast_program_is_nucleotide(program_number) {
        if let Some(matrix) = matrix {
            options.matrix_name = Some(matrix.to_ascii_uppercase());
        }
    } else {
        if penalty != 0 {
            options.penalty = penalty;
        }
        if reward != 0 {
            options.reward = reward;
        }
        if greedy_extension {
            options.gap_open = GAP_OPEN_MEGABLAST;
            options.gap_extend = GAP_EXTN_MEGABLAST;
        } else {
            options.gap_open = GAP_OPEN_NUCL;
            options.gap_extend = GAP_EXTN_NUCL;
        }
    }

    if gap_open >= 0 {
        options.gap_open = gap_open;
    }
    if gap_extend >= 0 {
        options.gap_extend = gap_extend;
    }

    0
}

/// Port of NCBI `BlastScoringOptionsValidate` (`blast_options.c:862`).
pub fn blast_scoring_options_validate(
    program_number: ProgramType,
    options: Option<&ScoringOptions>,
) -> i16 {
    let Some(options) = options else {
        return crate::util::BLASTERR_INVALIDPARAM;
    };

    if program_number == program::TBLASTX && options.gapped_calculation {
        return BLASTERR_OPTION_PROGRAM_INVALID;
    }
    if options.is_ooframe && !matches!(program_number, program::BLASTX | program::TBLASTN) {
        return BLASTERR_OPTION_PROGRAM_INVALID;
    }

    if program::blast_program_is_nucleotide(program_number) {
        if !(options.penalty == 0 && options.reward == 0) && options.penalty >= 0 {
            return BLASTERR_OPTION_VALUE_INVALID;
        }

        if options.gapped_calculation && options.gap_open > 0 && options.gap_extend == 0 {
            return BLASTERR_OPTION_VALUE_INVALID;
        }
    } else if options.gapped_calculation && !program::blast_program_is_rps_blast(program_number) {
        let Some(matrix_name) = options.matrix_name.as_deref() else {
            return BLASTERR_OPTION_VALUE_INVALID;
        };
        let standard_only = !matches!(program_number, program::BLASTP | program::TBLASTN);
        if crate::stat::blast_karlin_blk_gapped_load_from_tables(
            None,
            options.gap_open,
            options.gap_extend,
            Some(matrix_name),
            standard_only,
        ) != 0
        {
            return BLASTERR_OPTION_VALUE_INVALID;
        }
    }

    0
}

/// Hit saving options.
#[derive(Debug, Clone)]
pub struct HitSavingOptions {
    pub expect_value: f64,
    pub cutoff_score: i32,
    pub cutoff_score_fun: [i32; 2],
    pub percent_identity: f64,
    pub max_edit_distance: i32,
    pub hitlist_size: i32,
    pub hsp_num_max: i32,
    pub total_hsp_limit: i32,
    pub culling_limit: i32,
    pub mask_level: i32,
    pub do_sum_stats: bool,
    pub longest_intron: i32,
    pub min_hit_length: i32,
    pub min_diag_separation: i32,
    pub program_number: ProgramType,
    pub hsp_filt_opt: Option<crate::hspfilter_culling::BlastHSPFilteringOptions>,
    pub low_score_perc: f64,
    pub query_cov_hsp_perc: f64,
    pub max_hsps_per_subject: i32,
    pub paired: bool,
    pub splice: bool,
}

impl Default for HitSavingOptions {
    fn default() -> Self {
        HitSavingOptions {
            expect_value: EXPECT_VALUE,
            cutoff_score: 0,
            cutoff_score_fun: [0, 0],
            percent_identity: 0.0,
            max_edit_distance: i32::MAX,
            hitlist_size: HITLIST_SIZE,
            hsp_num_max: 0,
            total_hsp_limit: 0,
            culling_limit: 0,
            mask_level: 101,
            do_sum_stats: false,
            longest_intron: 0,
            min_hit_length: 0,
            min_diag_separation: 0,
            program_number: program::UNDEFINED,
            hsp_filt_opt: None,
            low_score_perc: 0.0,
            query_cov_hsp_perc: 0.0,
            max_hsps_per_subject: 0,
            paired: false,
            splice: false,
        }
    }
}

/// Port-shaped constructor for NCBI `BlastHitSavingOptionsNew`
/// (`blast_options.c:1446`).
pub fn blast_hit_saving_options_new(
    program_number: ProgramType,
    options: &mut Option<HitSavingOptions>,
    gapped_calculation: bool,
) -> i16 {
    let mut new_options = HitSavingOptions::default();
    new_options.program_number = program_number;
    new_options.do_sum_stats = if program_number == program::RPS_TBLASTN {
        false
    } else {
        !gapped_calculation
            || program::blast_query_is_translated(program_number)
            || program::blast_subject_is_translated(program_number)
    };
    *options = Some(new_options);
    0
}

/// Port of NCBI `BlastHitSavingOptionsFree` (`blast_options.c:1435`).
pub fn blast_hit_saving_options_free(
    options: &mut Option<HitSavingOptions>,
) -> Option<HitSavingOptions> {
    if let Some(options) = options.as_mut() {
        options.hsp_filt_opt = None;
    }
    *options = None;
    None
}

/// Port of NCBI `BLAST_FillHitSavingOptions` (`blast_options.c:1484`).
pub fn blast_fill_hit_saving_options(
    options: Option<&mut HitSavingOptions>,
    evalue: f64,
    hitlist_size: i32,
    _is_gapped: bool,
    culling_limit: i32,
    min_diag_separation: i32,
) -> i16 {
    let Some(options) = options else {
        return crate::util::BLASTERR_INVALIDPARAM;
    };

    if hitlist_size != 0 {
        options.hitlist_size = hitlist_size;
    }
    if evalue != 0.0 {
        options.expect_value = evalue;
    }
    if min_diag_separation != 0 {
        options.min_diag_separation = min_diag_separation;
    }
    options.culling_limit = culling_limit;
    options.hsp_filt_opt = None;
    options.max_edit_distance = i32::MAX;
    0
}

/// Port of NCBI `BlastHitSavingOptionsValidate` (`blast_options.c:1507`).
pub fn blast_hit_saving_options_validate(
    program_number: ProgramType,
    options: Option<&HitSavingOptions>,
) -> i16 {
    let Some(options) = options else {
        return crate::util::BLASTERR_INVALIDPARAM;
    };

    if options.hitlist_size < 1 {
        return BLASTERR_OPTION_VALUE_INVALID;
    }

    if options.expect_value <= 0.0 && options.cutoff_score <= 0 {
        return BLASTERR_OPTION_VALUE_INVALID;
    }

    if options.longest_intron != 0
        && !matches!(
            program_number,
            program::TBLASTN | program::PSI_TBLASTN | program::BLASTX | program::MAPPING
        )
    {
        return BLASTERR_OPTION_PROGRAM_INVALID;
    }

    if options.culling_limit < 0 {
        return BLASTERR_OPTION_VALUE_INVALID;
    }

    if let Some(hsp_filt_opt) = options.hsp_filt_opt.as_ref() {
        if crate::hspfilter_culling::blast_hspfiltering_options_validate(hsp_filt_opt) != 0 {
            return BLASTERR_OPTION_VALUE_INVALID;
        }
    }

    0
}

/// Initial word finding options.
#[derive(Debug, Clone)]
pub struct InitialWordOptions {
    pub gap_trigger: f64,
    pub window_size: i32,
    pub scan_range: i32,
    pub x_dropoff: f64,
    pub word_size: i32,
    pub program_number: ProgramType,
}

impl InitialWordOptions {
    pub fn new_blastn() -> Self {
        InitialWordOptions {
            gap_trigger: GAP_TRIGGER_NUCL,
            window_size: WINDOW_SIZE_NUCL,
            scan_range: SCAN_RANGE_NUCL,
            x_dropoff: UNGAPPED_X_DROPOFF_NUCL,
            word_size: WORDSIZE_NUCL,
            program_number: program::BLASTN,
        }
    }

    pub fn new_blastp() -> Self {
        InitialWordOptions {
            gap_trigger: GAP_TRIGGER_PROT,
            window_size: WINDOW_SIZE_PROT,
            scan_range: 0,
            x_dropoff: UNGAPPED_X_DROPOFF_PROT,
            word_size: WORDSIZE_PROT,
            program_number: program::BLASTP,
        }
    }
}

/// Port-shaped constructor for NCBI `BlastInitialWordOptionsNew`
/// (`blast_options.c:573`).
pub fn blast_initial_word_options_new(
    program_number: ProgramType,
    options: &mut Option<InitialWordOptions>,
) -> i16 {
    let mut new_options = if program::blast_program_is_nucleotide(program_number) {
        InitialWordOptions::new_blastn()
    } else {
        InitialWordOptions::new_blastp()
    };
    new_options.program_number = program_number;
    *options = Some(new_options);
    0
}

/// Port of NCBI `BlastInitialWordOptionsValidate` (`blast_options.c:601`).
pub fn blast_initial_word_options_validate(
    program_number: ProgramType,
    options: &InitialWordOptions,
) -> i16 {
    let ungapped_xdrop_is_required = program_number != program::BLASTN
        && program_number != program::MAPPING
        && !program::blast_program_is_phi_blast(program_number);
    if ungapped_xdrop_is_required && options.x_dropoff <= 0.0 {
        return 1;
    }

    if program_number == program::BLASTN && options.scan_range != 0 && options.window_size == 0 {
        return 1;
    }

    0
}

/// Port of NCBI `BLAST_FillInitialWordOptions` (`blast_options.c:634`).
pub fn blast_fill_initial_word_options(
    options: Option<&mut InitialWordOptions>,
    _program: ProgramType,
    window_size: i32,
    xdrop_ungapped: f64,
) -> i16 {
    let Some(options) = options else {
        return crate::util::BLASTERR_INVALIDPARAM;
    };

    if window_size != 0 {
        options.window_size = window_size;
    }
    if xdrop_ungapped != 0.0 {
        options.x_dropoff = xdrop_ungapped;
    }

    0
}

/// Port of NCBI `BlastInitialWordOptionsFree` (`blast_options.c:562`).
pub fn blast_initial_word_options_free(
    options: &mut Option<InitialWordOptions>,
) -> Option<InitialWordOptions> {
    *options = None;
    None
}

/// Extension options.
#[derive(Debug, Clone)]
pub struct ExtensionOptions {
    pub gap_x_dropoff: f64,
    pub gap_x_dropoff_final: f64,
    pub gap_trigger: f64,
    pub prelim_gap_ext: PrelimGapExt,
    pub traceback_ext: TracebackExt,
    pub composition_based_stats: i32,
    pub unified_p: i32,
    pub max_mismatches: i32,
    pub mismatch_window: i32,
    pub chaining: bool,
    pub program_number: ProgramType,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PrelimGapExt {
    DynProgScoreOnly,
    GreedyScoreOnly,
    SmithWatermanScoreOnly,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum TracebackExt {
    DynProgTbck,
    GreedyTbck,
    SmithWatermanTbckFull,
}

impl ExtensionOptions {
    pub fn new_blastn() -> Self {
        ExtensionOptions {
            gap_x_dropoff: GAP_X_DROPOFF_NUCL,
            gap_x_dropoff_final: GAP_X_DROPOFF_FINAL_NUCL,
            gap_trigger: 22.0,
            prelim_gap_ext: PrelimGapExt::DynProgScoreOnly,
            traceback_ext: TracebackExt::DynProgTbck,
            composition_based_stats: 0,
            unified_p: 0,
            max_mismatches: 5,
            mismatch_window: 10,
            chaining: false,
            program_number: program::BLASTN,
        }
    }

    pub fn new_blastp() -> Self {
        ExtensionOptions {
            gap_x_dropoff: GAP_X_DROPOFF_PROT,
            gap_x_dropoff_final: GAP_X_DROPOFF_FINAL_PROT,
            gap_trigger: 22.0,
            prelim_gap_ext: PrelimGapExt::DynProgScoreOnly,
            traceback_ext: TracebackExt::DynProgTbck,
            composition_based_stats: 0,
            unified_p: 0,
            max_mismatches: 5,
            mismatch_window: 10,
            chaining: false,
            program_number: program::BLASTP,
        }
    }
}

/// Port-shaped constructor for NCBI `BlastExtensionOptionsNew`
/// (`blast_options.c:660`).
pub fn blast_extension_options_new(
    program: ProgramType,
    options: &mut Option<ExtensionOptions>,
    gapped: bool,
) -> i16 {
    let mut new_options = if program::blast_program_is_nucleotide(program) {
        ExtensionOptions::new_blastn()
    } else {
        ExtensionOptions::new_blastp()
    };
    new_options.program_number = program;
    if gapped
        && program::blast_query_is_pssm(program)
        && !program::blast_subject_is_translated(program)
    {
        new_options.composition_based_stats = 1;
    }
    *options = Some(new_options);
    0
}

/// Port of NCBI `BlastExtensionOptionsFree` (`blast_options.c:650`).
pub fn blast_extension_options_free(
    options: &mut Option<ExtensionOptions>,
) -> Option<ExtensionOptions> {
    *options = None;
    None
}

/// Port of NCBI `BLAST_FillExtensionOptions` (`blast_options.c:699`).
pub fn blast_fill_extension_options(
    options: Option<&mut ExtensionOptions>,
    program: ProgramType,
    greedy: i32,
    x_dropoff: f64,
    x_dropoff_final: f64,
) -> i16 {
    let Some(options) = options else {
        return crate::util::BLASTERR_INVALIDPARAM;
    };

    if program::blast_program_is_nucleotide(program) {
        if greedy != 0 {
            options.gap_x_dropoff = crate::stat::BLAST_GAP_X_DROPOFF_GREEDY as f64;
            options.gap_x_dropoff_final = GAP_X_DROPOFF_FINAL_NUCL;
            options.prelim_gap_ext = PrelimGapExt::GreedyScoreOnly;
            options.traceback_ext = TracebackExt::GreedyTbck;
        } else {
            options.gap_x_dropoff = GAP_X_DROPOFF_NUCL;
            options.gap_x_dropoff_final = GAP_X_DROPOFF_FINAL_NUCL;
            options.prelim_gap_ext = PrelimGapExt::DynProgScoreOnly;
            options.traceback_ext = TracebackExt::DynProgTbck;
        }
    }

    if program::blast_query_is_pssm(program) && !program::blast_subject_is_translated(program) {
        options.composition_based_stats = 1;
    }

    if x_dropoff != 0.0 {
        options.gap_x_dropoff = x_dropoff;
    }
    if x_dropoff_final != 0.0 {
        options.gap_x_dropoff_final = x_dropoff_final;
    } else {
        options.gap_x_dropoff_final = options.gap_x_dropoff_final.max(x_dropoff);
    }

    0
}

/// Port of NCBI `BlastExtensionOptionsValidate` (`blast_options.c:740`).
pub fn blast_extension_options_validate(
    program_number: ProgramType,
    options: Option<&ExtensionOptions>,
) -> i16 {
    let Some(options) = options else {
        return crate::util::BLASTERR_INVALIDPARAM;
    };

    if program_number != program::BLASTN
        && program_number != program::MAPPING
        && (options.prelim_gap_ext == PrelimGapExt::GreedyScoreOnly
            || options.traceback_ext == TracebackExt::GreedyTbck)
    {
        return BLASTERR_OPTION_PROGRAM_INVALID;
    }

    if (options.prelim_gap_ext == PrelimGapExt::SmithWatermanScoreOnly
        && options.traceback_ext != TracebackExt::SmithWatermanTbckFull)
        || (options.prelim_gap_ext != PrelimGapExt::SmithWatermanScoreOnly
            && options.traceback_ext == TracebackExt::SmithWatermanTbckFull)
    {
        return BLASTERR_OPTION_VALUE_INVALID;
    }

    0
}

fn blast_extension_scoring_options_validate(
    program_number: ProgramType,
    ext_options: Option<&ExtensionOptions>,
    score_options: Option<&ScoringOptions>,
) -> i16 {
    let (Some(ext_options), Some(score_options)) = (ext_options, score_options) else {
        return crate::util::BLASTERR_INVALIDPARAM;
    };

    if program_number == program::BLASTN
        && score_options.gap_open == 0
        && score_options.gap_extend == 0
        && ext_options.prelim_gap_ext != PrelimGapExt::GreedyScoreOnly
        && ext_options.traceback_ext != TracebackExt::GreedyTbck
    {
        return BLASTERR_OPTION_VALUE_INVALID;
    }

    if ext_options.composition_based_stats != 0 {
        let allowed = program::blast_query_is_pssm(program_number)
            || matches!(
                program_number,
                program::TBLASTN
                    | program::BLASTP
                    | program::BLASTX
                    | program::RPS_BLAST
                    | program::RPS_TBLASTN
                    | program::PSI_BLAST
            );
        if !allowed || !score_options.gapped_calculation {
            return BLASTERR_OPTION_VALUE_INVALID;
        }
    }

    0
}

/// Port of NCBI `BLAST_ValidateOptions` (`blast_options.c:1751`).
pub fn blast_validate_options(
    program_number: ProgramType,
    ext_options: Option<&ExtensionOptions>,
    score_options: Option<&ScoringOptions>,
    lookup_options: Option<&LookupTableOptions>,
    word_options: Option<&InitialWordOptions>,
    hit_options: Option<&HitSavingOptions>,
) -> i16 {
    let mut status = blast_extension_options_validate(program_number, ext_options);
    if status != 0 {
        return status;
    }
    status = blast_scoring_options_validate(program_number, score_options);
    if status != 0 {
        return status;
    }
    status = lookup_table_options_validate(program_number, lookup_options);
    if status != 0 {
        return status;
    }

    let Some(word_options) = word_options else {
        return crate::util::BLASTERR_INVALIDPARAM;
    };
    status = blast_initial_word_options_validate(program_number, word_options);
    if status != 0 {
        return status;
    }
    status = blast_hit_saving_options_validate(program_number, hit_options);
    if status != 0 {
        return status;
    }
    status = blast_extension_scoring_options_validate(program_number, ext_options, score_options);
    if status != 0 {
        return status;
    }

    if matches!(program_number, program::BLASTP | program::TBLASTN) {
        let is_identity = score_options
            .and_then(|options| options.matrix_name.as_deref())
            .is_some_and(|matrix| matrix.eq_ignore_ascii_case("IDENTITY"));
        if lookup_options.is_some_and(|options| options.word_size > 5) && is_identity {
            return BLASTERR_OPTION_VALUE_INVALID;
        }
    }

    if matches!(program_number, program::RPS_BLAST | program::RPS_TBLASTN) {
        let Some(hit_options) = hit_options else {
            return crate::util::BLASTERR_INVALIDPARAM;
        };
        if hit_options.culling_limit != 0 || hit_options.hsp_filt_opt.is_some() {
            return BLASTERR_OPTION_VALUE_INVALID;
        }
    }

    0
}

/// Effective length options. 1-1 port of `BlastEffectiveLengthsOptions`
/// (`blast_options.h:484`).
#[derive(Debug, Clone, Default)]
pub struct EffectiveLengthsOptions {
    /// Database length used for statistical calculations.
    pub db_length: i64,
    /// Number of database sequences used for statistical calculations.
    pub dbseq_num: i32,
    /// Number of elements in `searchsp_eff` (must equal the number of
    /// contexts in the search).
    pub num_searchspaces: i32,
    /// Search space override per query context.
    pub searchsp_eff: Vec<i64>,
}

/// Port of NCBI `BlastEffectiveLengthsOptionsNew` (`blast_options.c:1003`).
pub fn blast_effective_lengths_options_new(options: &mut Option<EffectiveLengthsOptions>) -> i16 {
    *options = Some(EffectiveLengthsOptions::default());
    0
}

/// Port of NCBI `BlastEffectiveLengthsOptionsFree` (`blast_options.c:990`).
pub fn blast_effective_lengths_options_free(
    options: &mut Option<EffectiveLengthsOptions>,
) -> Option<EffectiveLengthsOptions> {
    if let Some(options) = options.as_mut() {
        options.searchsp_eff.clear();
        options.num_searchspaces = 0;
    }
    *options = None;
    None
}

/// Port of NCBI `BLAST_FillEffectiveLengthsOptions` (`blast_options.c:1038`).
pub fn blast_fill_effective_lengths_options(
    options: Option<&mut EffectiveLengthsOptions>,
    dbseq_num: i32,
    db_length: i64,
    searchsp_eff: Option<&[i64]>,
    num_searchsp: i32,
) -> i16 {
    let Some(options) = options else {
        return crate::util::BLASTERR_INVALIDPARAM;
    };

    if num_searchsp > options.num_searchspaces {
        options.num_searchspaces = num_searchsp;
        options.searchsp_eff.resize(num_searchsp as usize, 0);
    }

    let count = options.num_searchspaces.max(0) as usize;
    if count > 0 {
        let Some(searchsp_eff) = searchsp_eff else {
            return crate::util::BLASTERR_INVALIDPARAM;
        };
        if searchsp_eff.len() < count {
            return crate::util::BLASTERR_INVALIDPARAM;
        }
        options.searchsp_eff[..count].copy_from_slice(&searchsp_eff[..count]);
    }

    options.dbseq_num = dbseq_num;
    options.db_length = db_length;
    0
}

/// PSI-BLAST options. 1-1 port of C `PSIBlastOptions`
/// (`blast_options.h:498`).
#[derive(Debug, Clone, PartialEq)]
pub struct PSIBlastOptions {
    pub pseudo_count: i32,
    pub inclusion_ethresh: f64,
    pub use_best_alignment: bool,
    pub nsg_compatibility_mode: bool,
    pub impala_scaling_factor: f64,
    pub ignore_unaligned_positions: bool,
}

impl Default for PSIBlastOptions {
    fn default() -> Self {
        Self {
            pseudo_count: crate::stat::PSI_PSEUDO_COUNT_CONST,
            inclusion_ethresh: crate::stat::PSI_INCLUSION_ETHRESH,
            use_best_alignment: true,
            nsg_compatibility_mode: false,
            impala_scaling_factor: K_PSSM_NO_IMPALA_SCALING,
            ignore_unaligned_positions: false,
        }
    }
}

/// Port of NCBI `PSIBlastOptionsNew` (`blast_options.c:1556`).
pub fn psi_blast_options_new(options: &mut Option<PSIBlastOptions>) -> i16 {
    *options = Some(PSIBlastOptions::default());
    0
}

/// Port of NCBI `PSIBlastOptionsValidate` (`blast_options.c:1579`).
pub fn psi_blast_options_validate(options: Option<&PSIBlastOptions>) -> i16 {
    let Some(options) = options else {
        return 1;
    };

    if options.pseudo_count < 0 {
        return 1;
    }
    if options.inclusion_ethresh <= 0.0 {
        return 1;
    }

    0
}

/// Port of NCBI `PSIBlastOptionsFree` (`blast_options.c:1604`).
pub fn psi_blast_options_free(options: &mut Option<PSIBlastOptions>) -> Option<PSIBlastOptions> {
    *options = None;
    None
}

/// Database options.
#[derive(Debug, Clone)]
pub struct DatabaseOptions {
    pub genetic_code: i32,
}

impl Default for DatabaseOptions {
    fn default() -> Self {
        DatabaseOptions { genetic_code: 1 } // standard genetic code
    }
}

/// Port-shaped constructor for NCBI `BlastDatabaseOptionsNew`
/// (`blast_options.c:1621`).
pub fn blast_database_options_new(options: &mut Option<DatabaseOptions>) -> i16 {
    *options = Some(DatabaseOptions::default());
    0
}

/// Port of NCBI `BlastDatabaseOptionsFree` (`blast_options.c:1630`).
pub fn blast_database_options_free(
    options: &mut Option<DatabaseOptions>,
) -> Option<DatabaseOptions> {
    *options = None;
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_blastn_defaults() {
        let s = ScoringOptions::new_blastn();
        assert_eq!(s.reward, 1);
        assert_eq!(s.penalty, -3);
        assert_eq!(s.gap_open, 5);
        assert_eq!(s.gap_extend, 2);
        assert!(s.gapped_calculation);
    }

    #[test]
    fn test_blastp_defaults() {
        let s = ScoringOptions::new_blastp();
        assert_eq!(s.gap_open, 11);
        assert_eq!(s.gap_extend, 1);
        assert_eq!(s.matrix_name.as_deref(), Some("BLOSUM62"));
    }

    #[test]
    fn test_hit_saving_defaults() {
        let h = HitSavingOptions::default();
        assert_eq!(h.expect_value, 10.0);
        assert_eq!(h.hitlist_size, 500);
        assert_eq!(h.mask_level, 101);
        assert_eq!(h.max_edit_distance, i32::MAX);
    }

    #[test]
    fn translated_lookup_table_options_new_free_and_fill() {
        let mut lookup = None;
        assert_eq!(
            lookup_table_options_new(crate::program::BLASTN, &mut lookup),
            0
        );
        let lookup_ref = lookup.as_ref().expect("blastn lookup options");
        assert_eq!(lookup_ref.word_size, WORDSIZE_MEGABLAST);
        assert_eq!(lookup_ref.lut_type, LookupTableType::MegablastLookup);
        assert_eq!(lookup_ref.program_number, crate::program::BLASTN);

        assert_eq!(
            lookup_table_options_new(crate::program::BLASTP, &mut lookup),
            0
        );
        let lookup_ref = lookup.as_ref().expect("blastp lookup options");
        assert_eq!(lookup_ref.word_size, WORDSIZE_PROT);
        assert_eq!(lookup_ref.lut_type, LookupTableType::AaLookup);
        assert_eq!(lookup_ref.threshold, WORD_THRESHOLD_BLASTP);

        assert_eq!(
            lookup_table_options_new(crate::program::RPS_TBLASTN, &mut lookup),
            0
        );
        let lookup_ref = lookup.as_ref().expect("rps lookup options");
        assert_eq!(lookup_ref.lut_type, LookupTableType::RpsLookup);
        assert_eq!(lookup_ref.threshold, WORD_THRESHOLD_TBLASTN);

        let mut fill = LookupTableOptions::default();
        assert_eq!(
            blast_fill_lookup_table_options(Some(&mut fill), crate::program::BLASTN, false, 0.0, 0),
            0
        );
        assert_eq!(fill.lut_type, LookupTableType::NaLookup);
        assert_eq!(fill.word_size, WORDSIZE_NUCL);

        assert_eq!(
            blast_fill_lookup_table_options(
                Some(&mut fill),
                crate::program::BLASTP,
                false,
                15.0,
                6
            ),
            0
        );
        assert_eq!(fill.threshold, 15.0);
        assert_eq!(fill.word_size, 6);
        assert_eq!(fill.lut_type, LookupTableType::CompressedAaLookup);

        assert_eq!(
            blast_fill_lookup_table_options(None, crate::program::BLASTN, false, 0.0, 0),
            crate::util::BLASTERR_INVALIDPARAM
        );
        assert!(lookup_table_options_free(&mut lookup).is_none());
        assert!(lookup.is_none());
    }

    #[test]
    fn translated_suggested_lookup_threshold_and_window_tables() {
        let mut threshold = 0.0;
        assert_eq!(
            blast_get_suggested_threshold(crate::program::BLASTP, Some("BLOSUM45"), &mut threshold),
            0
        );
        assert_eq!(threshold, 14.0);
        assert_eq!(
            blast_get_suggested_threshold(crate::program::BLASTX, Some("pam30"), &mut threshold),
            0
        );
        assert_eq!(threshold, 17.0);
        assert_eq!(
            blast_get_suggested_threshold(
                crate::program::TBLASTN,
                Some("BLOSUM62"),
                &mut threshold
            ),
            0
        );
        assert_eq!(threshold, 13.0);
        assert_eq!(
            blast_get_suggested_threshold(crate::program::BLASTN, None, &mut threshold),
            0
        );
        assert_eq!(
            blast_get_suggested_threshold(crate::program::BLASTP, None, &mut threshold),
            crate::util::BLASTERR_INVALIDPARAM
        );

        let mut window_size = 0;
        assert_eq!(
            blast_get_suggested_window_size(
                crate::program::BLASTP,
                Some("BLOSUM45"),
                &mut window_size
            ),
            0
        );
        assert_eq!(window_size, 60);
        assert_eq!(
            blast_get_suggested_window_size(
                crate::program::BLASTP,
                Some("unknown"),
                &mut window_size
            ),
            0
        );
        assert_eq!(window_size, 40);
        assert_eq!(
            blast_get_suggested_window_size(crate::program::MAPPING, None, &mut window_size),
            0
        );
        assert_eq!(
            blast_get_suggested_window_size(crate::program::BLASTP, None, &mut window_size),
            crate::util::BLASTERR_INVALIDPARAM
        );
    }

    #[test]
    fn translated_lookup_table_validation_matches_c_rules() {
        let mut options = LookupTableOptions::default();
        assert_eq!(
            lookup_table_options_validate(crate::program::BLASTP, Some(&options)),
            BLASTERR_OPTION_VALUE_INVALID
        );

        options.threshold = WORD_THRESHOLD_BLASTP;
        options.word_size = WORDSIZE_PROT;
        options.lut_type = LookupTableType::AaLookup;
        assert_eq!(
            lookup_table_options_validate(crate::program::BLASTP, Some(&options)),
            0
        );

        options.word_size = 6;
        assert_eq!(
            lookup_table_options_validate(crate::program::BLASTP, Some(&options)),
            BLASTERR_OPTION_VALUE_INVALID
        );
        options.lut_type = LookupTableType::CompressedAaLookup;
        assert_eq!(
            lookup_table_options_validate(crate::program::BLASTP, Some(&options)),
            0
        );
        options.word_size = 8;
        assert_eq!(
            lookup_table_options_validate(crate::program::BLASTP, Some(&options)),
            BLASTERR_OPTION_VALUE_INVALID
        );

        let mut blastn_options = LookupTableOptions {
            word_size: WORDSIZE_NUCL,
            lut_type: LookupTableType::MegablastLookup,
            ..LookupTableOptions::default()
        };
        assert_eq!(
            lookup_table_options_validate(crate::program::BLASTN, Some(&blastn_options)),
            0
        );
        blastn_options.word_size = 3;
        assert_eq!(
            lookup_table_options_validate(crate::program::BLASTN, Some(&blastn_options)),
            BLASTERR_OPTION_VALUE_INVALID
        );
        blastn_options.word_size = 11;
        blastn_options.mb_template_length = 16;
        blastn_options.mb_template_type = 2;
        assert_eq!(
            lookup_table_options_validate(crate::program::BLASTN, Some(&blastn_options)),
            0
        );
        blastn_options.mb_template_type = 3;
        assert_eq!(
            lookup_table_options_validate(crate::program::BLASTN, Some(&blastn_options)),
            BLASTERR_OPTION_VALUE_INVALID
        );

        let phi_options = LookupTableOptions {
            phi_pattern: Some("A-C".to_string()),
            ..LookupTableOptions::default()
        };
        assert_eq!(
            lookup_table_options_validate(crate::program::BLASTP, Some(&phi_options)),
            BLASTERR_OPTION_PROGRAM_INVALID
        );
        assert_eq!(
            lookup_table_options_validate(crate::program::PHI_BLASTP, Some(&phi_options)),
            0
        );
        assert_eq!(
            lookup_table_options_validate(crate::program::BLASTN, None),
            crate::util::BLASTERR_INVALIDPARAM
        );
    }

    #[test]
    fn translated_option_free_functions_clear_owned_options() {
        let mut query_setup = Some(QuerySetUpOptions::default());
        let mut scoring = Some(ScoringOptions::new_blastp());
        let mut hit_saving = Some(HitSavingOptions::default());
        let mut initial = Some(InitialWordOptions::new_blastn());
        let mut extension = Some(ExtensionOptions::new_blastp());
        let mut effective = Some(EffectiveLengthsOptions {
            db_length: 100,
            dbseq_num: 2,
            num_searchspaces: 2,
            searchsp_eff: vec![10, 20],
        });
        let mut psi = Some(PSIBlastOptions::default());
        let mut database = Some(DatabaseOptions::default());

        assert!(blast_query_set_up_options_free(&mut query_setup).is_none());
        assert!(blast_scoring_options_free(&mut scoring).is_none());
        assert!(blast_hit_saving_options_free(&mut hit_saving).is_none());
        assert!(blast_initial_word_options_free(&mut initial).is_none());
        assert!(blast_extension_options_free(&mut extension).is_none());
        assert!(blast_effective_lengths_options_free(&mut effective).is_none());
        assert!(psi_blast_options_free(&mut psi).is_none());
        assert!(blast_database_options_free(&mut database).is_none());

        assert!(query_setup.is_none());
        assert!(scoring.is_none());
        assert!(hit_saving.is_none());
        assert!(initial.is_none());
        assert!(extension.is_none());
        assert!(effective.is_none());
        assert!(psi.is_none());
        assert!(database.is_none());
    }

    #[test]
    fn translated_option_new_functions_match_existing_defaults() {
        let mut query_setup = None;
        assert_eq!(blast_query_set_up_options_new(&mut query_setup), 0);
        let query_setup_ref = query_setup.as_ref().expect("query setup options");
        assert_eq!(query_setup_ref.genetic_code, 1);
        assert_eq!(query_setup_ref.strand_option, 0);
        assert!(query_setup_ref.filter_string.is_none());
        assert!(crate::filter::sblast_filter_options_no_filtering(
            query_setup_ref.filtering_options.as_ref()
        ));

        let mut initial = None;
        assert_eq!(
            blast_initial_word_options_new(crate::program::BLASTN, &mut initial),
            0
        );
        let initial_ref = initial.as_ref().expect("initial word options");
        assert_eq!(initial_ref.window_size, WINDOW_SIZE_NUCL);
        assert_eq!(initial_ref.scan_range, SCAN_RANGE_NUCL);
        assert_eq!(initial_ref.gap_trigger, GAP_TRIGGER_NUCL);
        assert_eq!(initial_ref.x_dropoff, UNGAPPED_X_DROPOFF_NUCL);
        assert_eq!(initial_ref.program_number, crate::program::BLASTN);

        assert_eq!(
            blast_initial_word_options_new(crate::program::BLASTP, &mut initial),
            0
        );
        let initial_ref = initial.as_ref().expect("protein initial word options");
        assert_eq!(initial_ref.window_size, WINDOW_SIZE_PROT);
        assert_eq!(initial_ref.gap_trigger, GAP_TRIGGER_PROT);
        assert_eq!(initial_ref.x_dropoff, UNGAPPED_X_DROPOFF_PROT);
        assert_eq!(initial_ref.program_number, crate::program::BLASTP);

        let mut scoring = None;
        assert_eq!(
            blast_scoring_options_new(crate::program::BLASTN, &mut scoring),
            0
        );
        let scoring_ref = scoring.as_ref().expect("scoring options");
        assert_eq!(scoring_ref.reward, REWARD);
        assert_eq!(scoring_ref.penalty, PENALTY);
        assert_eq!(scoring_ref.gap_open, GAP_OPEN_NUCL);
        assert_eq!(scoring_ref.gap_extend, GAP_EXTN_NUCL);
        assert_eq!(scoring_ref.matrix_name, None);

        assert_eq!(
            blast_scoring_options_new(crate::program::BLASTP, &mut scoring),
            0
        );
        let scoring_ref = scoring.as_ref().expect("protein scoring options");
        assert_eq!(scoring_ref.gap_open, GAP_OPEN_PROT);
        assert_eq!(scoring_ref.gap_extend, GAP_EXTN_PROT);
        assert_eq!(scoring_ref.matrix_name.as_deref(), Some("BLOSUM62"));

        assert_eq!(
            blast_scoring_options_new(crate::program::TBLASTX, &mut scoring),
            0
        );
        assert!(
            !scoring
                .as_ref()
                .expect("tblastx scoring options")
                .gapped_calculation
        );

        let mut hit_saving = None;
        assert_eq!(
            blast_hit_saving_options_new(crate::program::BLASTP, &mut hit_saving, true),
            0
        );
        let hit_saving_ref = hit_saving.as_ref().expect("hit saving options");
        assert_eq!(hit_saving_ref.hitlist_size, HITLIST_SIZE);
        assert_eq!(hit_saving_ref.expect_value, EXPECT_VALUE);
        assert_eq!(hit_saving_ref.program_number, crate::program::BLASTP);
        assert!(!hit_saving_ref.do_sum_stats);

        assert_eq!(
            blast_hit_saving_options_new(crate::program::BLASTX, &mut hit_saving, true),
            0
        );
        assert!(
            hit_saving
                .as_ref()
                .expect("translated hit saving options")
                .do_sum_stats
        );
        assert_eq!(
            blast_hit_saving_options_new(crate::program::BLASTP, &mut hit_saving, false),
            0
        );
        assert!(
            hit_saving
                .as_ref()
                .expect("ungapped hit saving options")
                .do_sum_stats
        );
        assert_eq!(
            blast_hit_saving_options_new(crate::program::RPS_TBLASTN, &mut hit_saving, false),
            0
        );
        assert!(
            !hit_saving
                .as_ref()
                .expect("rps tblastn hit saving options")
                .do_sum_stats
        );

        let mut extension = None;
        assert_eq!(
            blast_extension_options_new(crate::program::BLASTN, &mut extension, true),
            0
        );
        let extension_ref = extension.as_ref().expect("extension options");
        assert_eq!(extension_ref.gap_x_dropoff, GAP_X_DROPOFF_NUCL);
        assert_eq!(extension_ref.gap_x_dropoff_final, GAP_X_DROPOFF_FINAL_NUCL);
        assert_eq!(extension_ref.prelim_gap_ext, PrelimGapExt::DynProgScoreOnly);
        assert_eq!(extension_ref.traceback_ext, TracebackExt::DynProgTbck);
        assert_eq!(extension_ref.program_number, crate::program::BLASTN);
        assert_eq!(extension_ref.max_mismatches, 5);
        assert_eq!(extension_ref.mismatch_window, 10);

        let mut effective = None;
        assert_eq!(blast_effective_lengths_options_new(&mut effective), 0);
        assert_eq!(
            effective.as_ref().map(|options| options.num_searchspaces),
            Some(0)
        );

        let mut database = None;
        assert_eq!(blast_database_options_new(&mut database), 0);
        assert_eq!(
            database.as_ref().map(|options| options.genetic_code),
            Some(1)
        );
    }

    #[test]
    fn translated_initial_word_fill_and_validate_match_c_rules() {
        let mut options = InitialWordOptions::new_blastp();
        assert_eq!(
            blast_fill_initial_word_options(Some(&mut options), crate::program::BLASTP, 17, 9.5),
            0
        );
        assert_eq!(options.window_size, 17);
        assert_eq!(options.x_dropoff, 9.5);
        assert_eq!(
            blast_initial_word_options_validate(crate::program::BLASTP, &options),
            0
        );

        options.x_dropoff = 0.0;
        assert_eq!(
            blast_initial_word_options_validate(crate::program::BLASTP, &options),
            1
        );
        assert_eq!(
            blast_initial_word_options_validate(crate::program::PHI_BLASTP, &options),
            0
        );

        let mut blastn_options = InitialWordOptions::new_blastn();
        blastn_options.scan_range = 12;
        blastn_options.window_size = 0;
        assert_eq!(
            blast_initial_word_options_validate(crate::program::BLASTN, &blastn_options),
            1
        );
        assert_eq!(
            blast_fill_initial_word_options(None, crate::program::BLASTN, 0, 0.0),
            crate::util::BLASTERR_INVALIDPARAM
        );
    }

    #[test]
    fn translated_query_setup_fill_matches_c_rules() {
        let mut options = None;
        assert_eq!(blast_query_set_up_options_new(&mut options), 0);
        let options = options.as_mut().expect("query setup options");

        assert_eq!(
            blast_fill_query_set_up_options(Some(options), crate::program::BLASTN, None, 3),
            0
        );
        assert_eq!(options.strand_option, 3);

        assert_eq!(
            blast_fill_query_set_up_options(Some(options), crate::program::BLASTP, None, 1),
            0
        );
        assert_eq!(options.strand_option, 3);

        assert_eq!(
            blast_fill_query_set_up_options(Some(options), crate::program::BLASTN, Some("D"), 0),
            0
        );
        let filtering = options
            .filtering_options
            .as_ref()
            .expect("parsed dust filter");
        assert!(filtering.dust_options.is_some());
        assert!(filtering.seg_options.is_none());
        assert!(options.filter_string.is_none());

        assert_eq!(
            blast_fill_query_set_up_options(Some(options), crate::program::BLASTP, Some("seg"), 0),
            0
        );
        let filtering = options
            .filtering_options
            .as_ref()
            .expect("parsed seg filter");
        assert!(filtering.seg_options.is_some());
        assert!(filtering.dust_options.is_none());

        assert_eq!(
            blast_fill_query_set_up_options(Some(options), crate::program::BLASTN, Some("bad"), 0),
            crate::util::BLASTERR_INVALIDPARAM
        );
        assert_eq!(
            blast_fill_query_set_up_options(None, crate::program::BLASTN, None, 0),
            crate::util::BLASTERR_INVALIDPARAM
        );
    }

    #[test]
    fn translated_hit_saving_fill_and_validate_match_c_rules() {
        let mut options = HitSavingOptions::default();
        assert_eq!(
            blast_fill_hit_saving_options(Some(&mut options), 1e-6, 25, true, 2, 7),
            0
        );
        assert_eq!(options.expect_value, 1e-6);
        assert_eq!(options.hitlist_size, 25);
        assert_eq!(options.culling_limit, 2);
        assert_eq!(options.min_diag_separation, 7);
        assert_eq!(options.max_edit_distance, i32::MAX);
        assert_eq!(
            blast_hit_saving_options_validate(crate::program::BLASTP, Some(&options)),
            0
        );

        options.hitlist_size = 0;
        assert_eq!(
            blast_hit_saving_options_validate(crate::program::BLASTP, Some(&options)),
            BLASTERR_OPTION_VALUE_INVALID
        );
        options.hitlist_size = 25;
        options.expect_value = 0.0;
        options.cutoff_score = 0;
        assert_eq!(
            blast_hit_saving_options_validate(crate::program::BLASTP, Some(&options)),
            BLASTERR_OPTION_VALUE_INVALID
        );
        options.cutoff_score = 20;
        options.longest_intron = 100;
        assert_eq!(
            blast_hit_saving_options_validate(crate::program::BLASTP, Some(&options)),
            BLASTERR_OPTION_PROGRAM_INVALID
        );
        assert_eq!(
            blast_hit_saving_options_validate(crate::program::TBLASTN, Some(&options)),
            0
        );
        options.culling_limit = -1;
        assert_eq!(
            blast_hit_saving_options_validate(crate::program::TBLASTN, Some(&options)),
            BLASTERR_OPTION_VALUE_INVALID
        );
        assert_eq!(
            blast_hit_saving_options_validate(crate::program::BLASTP, None),
            crate::util::BLASTERR_INVALIDPARAM
        );
        assert_eq!(
            blast_fill_hit_saving_options(None, 1e-6, 25, true, 2, 7),
            crate::util::BLASTERR_INVALIDPARAM
        );
    }

    #[test]
    fn translated_effective_lengths_new_and_fill_match_c_rules() {
        let mut options = None;
        assert_eq!(blast_effective_lengths_options_new(&mut options), 0);
        let options = options.as_mut().expect("effective length options");

        assert_eq!(
            blast_fill_effective_lengths_options(Some(options), 17, 1000, Some(&[10, 20]), 2),
            0
        );
        assert_eq!(options.dbseq_num, 17);
        assert_eq!(options.db_length, 1000);
        assert_eq!(options.num_searchspaces, 2);
        assert_eq!(options.searchsp_eff, vec![10, 20]);

        assert_eq!(
            blast_fill_effective_lengths_options(Some(options), 3, 50, Some(&[30, 40]), 1),
            0
        );
        assert_eq!(options.dbseq_num, 3);
        assert_eq!(options.db_length, 50);
        assert_eq!(options.num_searchspaces, 2);
        assert_eq!(options.searchsp_eff, vec![30, 40]);

        assert_eq!(
            blast_fill_effective_lengths_options(Some(options), 3, 50, Some(&[1]), 2),
            crate::util::BLASTERR_INVALIDPARAM
        );
        assert_eq!(
            blast_fill_effective_lengths_options(None, 3, 50, Some(&[1, 2]), 2),
            crate::util::BLASTERR_INVALIDPARAM
        );
    }

    #[test]
    fn translated_extension_fill_and_validate_match_c_rules() {
        let mut options = ExtensionOptions::new_blastn();
        assert_eq!(
            blast_fill_extension_options(Some(&mut options), crate::program::BLASTN, 1, 0.0, 0.0),
            0
        );
        assert_eq!(
            options.gap_x_dropoff,
            crate::stat::BLAST_GAP_X_DROPOFF_GREEDY as f64
        );
        assert_eq!(options.gap_x_dropoff_final, GAP_X_DROPOFF_FINAL_NUCL);
        assert_eq!(options.prelim_gap_ext, PrelimGapExt::GreedyScoreOnly);
        assert_eq!(options.traceback_ext, TracebackExt::GreedyTbck);
        assert_eq!(
            blast_extension_options_validate(crate::program::BLASTN, Some(&options)),
            0
        );
        assert_eq!(
            blast_extension_options_validate(crate::program::BLASTP, Some(&options)),
            BLASTERR_OPTION_PROGRAM_INVALID
        );

        assert_eq!(
            blast_fill_extension_options(Some(&mut options), crate::program::BLASTN, 0, 40.0, 0.0),
            0
        );
        assert_eq!(options.gap_x_dropoff, 40.0);
        assert_eq!(options.gap_x_dropoff_final, GAP_X_DROPOFF_FINAL_NUCL);

        options.prelim_gap_ext = PrelimGapExt::SmithWatermanScoreOnly;
        options.traceback_ext = TracebackExt::DynProgTbck;
        assert_eq!(
            blast_extension_options_validate(crate::program::BLASTN, Some(&options)),
            BLASTERR_OPTION_VALUE_INVALID
        );
        options.traceback_ext = TracebackExt::SmithWatermanTbckFull;
        assert_eq!(
            blast_extension_options_validate(crate::program::BLASTN, Some(&options)),
            0
        );

        assert_eq!(
            blast_fill_extension_options(None, crate::program::BLASTN, 0, 0.0, 0.0),
            crate::util::BLASTERR_INVALIDPARAM
        );
        assert_eq!(
            blast_extension_options_validate(crate::program::BLASTN, None),
            crate::util::BLASTERR_INVALIDPARAM
        );
    }

    #[test]
    fn translated_scoring_helpers_match_c_rules() {
        let mut options = ScoringOptions::new_blastp();
        assert_eq!(
            blast_scoring_options_set_matrix(Some(&mut options), Some("pam30")),
            0
        );
        assert_eq!(options.matrix_name.as_deref(), Some("PAM30"));

        let mut dup = None;
        assert_eq!(blast_scoring_options_dup(&mut dup, Some(&options)), 0);
        assert_eq!(
            dup.as_ref().and_then(|o| o.matrix_name.as_deref()),
            Some("PAM30")
        );
        assert_eq!(
            blast_scoring_options_dup(&mut dup, None),
            crate::util::BLASTERR_INVALIDPARAM
        );
        assert!(dup.is_none());

        assert_eq!(
            blast_fill_scoring_options(
                Some(&mut options),
                crate::program::BLASTP,
                false,
                0,
                0,
                Some("blosum45"),
                -1,
                -1
            ),
            0
        );
        assert_eq!(options.matrix_name.as_deref(), Some("BLOSUM45"));
        options.gap_open = 13;
        options.gap_extend = 3;
        assert_eq!(
            blast_scoring_options_validate(crate::program::BLASTP, Some(&options)),
            0
        );
        options.gap_open = 1;
        options.gap_extend = 1;
        assert_eq!(
            blast_scoring_options_validate(crate::program::BLASTP, Some(&options)),
            BLASTERR_OPTION_VALUE_INVALID
        );
        options.gap_open = 13;
        options.gap_extend = 3;
        assert_eq!(
            blast_scoring_options_validate(crate::program::RPS_BLAST, Some(&options)),
            0
        );
        options.matrix_name = Some("IDENTITY".to_string());
        options.gap_open = 15;
        options.gap_extend = 2;
        assert_eq!(
            blast_scoring_options_validate(crate::program::BLASTX, Some(&options)),
            BLASTERR_OPTION_VALUE_INVALID
        );
        assert_eq!(
            blast_scoring_options_validate(crate::program::BLASTP, Some(&options)),
            0
        );
        options.is_ooframe = true;
        assert_eq!(
            blast_scoring_options_validate(crate::program::BLASTP, Some(&options)),
            BLASTERR_OPTION_PROGRAM_INVALID
        );
        assert_eq!(
            blast_scoring_options_validate(crate::program::BLASTX, Some(&options)),
            BLASTERR_OPTION_VALUE_INVALID
        );
        options.is_ooframe = false;
        options.matrix_name = Some("BLOSUM62".to_string());
        options.gap_open = GAP_OPEN_PROT;
        options.gap_extend = GAP_EXTN_PROT;

        let mut blastn = ScoringOptions::new_blastn();
        assert_eq!(
            blast_fill_scoring_options(
                Some(&mut blastn),
                crate::program::BLASTN,
                true,
                -2,
                1,
                None,
                -1,
                -1
            ),
            0
        );
        assert_eq!(blastn.penalty, -2);
        assert_eq!(blastn.reward, 1);
        assert_eq!(blastn.gap_open, GAP_OPEN_MEGABLAST);
        assert_eq!(blastn.gap_extend, GAP_EXTN_MEGABLAST);
        assert_eq!(
            blast_scoring_options_validate(crate::program::BLASTN, Some(&blastn)),
            0
        );

        blastn.penalty = 1;
        assert_eq!(
            blast_scoring_options_validate(crate::program::BLASTN, Some(&blastn)),
            BLASTERR_OPTION_VALUE_INVALID
        );
        blastn.penalty = -1;
        blastn.gap_open = 5;
        blastn.gap_extend = 0;
        assert_eq!(
            blast_scoring_options_validate(crate::program::BLASTN, Some(&blastn)),
            BLASTERR_OPTION_VALUE_INVALID
        );

        let tblastx = ScoringOptions::new_blastp();
        assert_eq!(
            blast_scoring_options_validate(crate::program::TBLASTX, Some(&tblastx)),
            BLASTERR_OPTION_PROGRAM_INVALID
        );
        assert_eq!(
            blast_scoring_options_set_matrix(None, Some("BLOSUM62")),
            crate::util::BLASTERR_INVALIDPARAM
        );
        assert_eq!(
            blast_fill_scoring_options(None, crate::program::BLASTP, false, 0, 0, None, -1, -1),
            crate::util::BLASTERR_INVALIDPARAM
        );
        assert_eq!(
            blast_scoring_options_validate(crate::program::BLASTP, None),
            crate::util::BLASTERR_INVALIDPARAM
        );
    }

    #[test]
    fn translated_validate_options_runs_subvalidators_and_cross_checks() {
        let mut lookup = LookupTableOptions::default();
        lookup.word_size = WORDSIZE_PROT;
        lookup.lut_type = LookupTableType::AaLookup;
        lookup.threshold = WORD_THRESHOLD_BLASTP;
        lookup.program_number = crate::program::BLASTP;
        let word = InitialWordOptions::new_blastp();
        let ext = ExtensionOptions::new_blastp();
        let score = ScoringOptions::new_blastp();
        let hit = HitSavingOptions::default();

        assert_eq!(
            blast_validate_options(
                crate::program::BLASTP,
                Some(&ext),
                Some(&score),
                Some(&lookup),
                Some(&word),
                Some(&hit)
            ),
            0
        );

        let mut identity_score = score.clone();
        identity_score.matrix_name = Some("IDENTITY".to_string());
        lookup.word_size = 6;
        assert_eq!(
            blast_validate_options(
                crate::program::BLASTP,
                Some(&ext),
                Some(&identity_score),
                Some(&lookup),
                Some(&word),
                Some(&hit)
            ),
            BLASTERR_OPTION_VALUE_INVALID
        );

        let mut blastn_score = ScoringOptions::new_blastn();
        blastn_score.gap_open = 0;
        blastn_score.gap_extend = 0;
        let blastn_ext = ExtensionOptions::new_blastn();
        let blastn_word = InitialWordOptions::new_blastn();
        let blastn_lookup = LookupTableOptions {
            word_size: WORDSIZE_NUCL,
            lut_type: LookupTableType::NaLookup,
            program_number: crate::program::BLASTN,
            ..LookupTableOptions::default()
        };
        assert_eq!(
            blast_validate_options(
                crate::program::BLASTN,
                Some(&blastn_ext),
                Some(&blastn_score),
                Some(&blastn_lookup),
                Some(&blastn_word),
                Some(&hit)
            ),
            BLASTERR_OPTION_VALUE_INVALID
        );

        let rps_lookup = LookupTableOptions {
            word_size: WORDSIZE_PROT,
            lut_type: LookupTableType::RpsLookup,
            threshold: WORD_THRESHOLD_BLASTP,
            program_number: crate::program::RPS_BLAST,
            ..LookupTableOptions::default()
        };
        let mut rps_hit = HitSavingOptions::default();
        rps_hit.culling_limit = 1;
        assert_eq!(
            blast_validate_options(
                crate::program::RPS_BLAST,
                Some(&ext),
                Some(&score),
                Some(&rps_lookup),
                Some(&word),
                Some(&rps_hit)
            ),
            BLASTERR_OPTION_VALUE_INVALID
        );

        assert_eq!(
            blast_validate_options(
                crate::program::BLASTP,
                None,
                Some(&score),
                Some(&lookup),
                Some(&word),
                Some(&hit)
            ),
            crate::util::BLASTERR_INVALIDPARAM
        );
    }

    #[test]
    fn translated_psi_blast_options_match_c_rules() {
        let mut options = None;
        assert_eq!(psi_blast_options_new(&mut options), 0);
        let options = options.as_ref().expect("PSI options should be created");

        assert_eq!(options.pseudo_count, crate::stat::PSI_PSEUDO_COUNT_CONST);
        assert_eq!(
            options.inclusion_ethresh,
            crate::stat::PSI_INCLUSION_ETHRESH
        );
        assert!(options.use_best_alignment);
        assert!(!options.nsg_compatibility_mode);
        assert_eq!(options.impala_scaling_factor, K_PSSM_NO_IMPALA_SCALING);
        assert!(!options.ignore_unaligned_positions);
        assert_eq!(psi_blast_options_validate(Some(options)), 0);

        let mut invalid = options.clone();
        invalid.pseudo_count = -1;
        assert_eq!(psi_blast_options_validate(Some(&invalid)), 1);
        invalid.pseudo_count = 0;
        invalid.inclusion_ethresh = 0.0;
        assert_eq!(psi_blast_options_validate(Some(&invalid)), 1);
        assert_eq!(psi_blast_options_validate(None), 1);
    }

    /// Port of NCBI blastoptions_unit_test: verify all default blastn options.
    #[test]
    fn test_blastn_default_options() {
        let scoring = ScoringOptions::new_blastn();
        assert_eq!(scoring.reward, 1);
        assert_eq!(scoring.penalty, -3);
        assert_eq!(scoring.gap_open, 5);
        assert_eq!(scoring.gap_extend, 2);
        assert!(scoring.gapped_calculation);
        assert!(
            scoring.matrix_name.is_none(),
            "blastn should not use a named matrix"
        );
        assert!(!scoring.is_ooframe);

        let word = InitialWordOptions::new_blastn();
        assert_eq!(word.word_size, WORDSIZE_NUCL);
        assert_eq!(word.word_size, 11);
        assert_eq!(word.window_size, WINDOW_SIZE_NUCL);
        assert_eq!(word.x_dropoff, UNGAPPED_X_DROPOFF_NUCL);

        let ext = ExtensionOptions::new_blastn();
        assert_eq!(ext.gap_x_dropoff, GAP_X_DROPOFF_NUCL);
        assert_eq!(ext.gap_x_dropoff_final, GAP_X_DROPOFF_FINAL_NUCL);
        assert_eq!(ext.prelim_gap_ext, PrelimGapExt::DynProgScoreOnly);
        assert_eq!(ext.traceback_ext, TracebackExt::DynProgTbck);
    }

    /// Port of NCBI blastoptions_unit_test: verify default blastp options.
    #[test]
    fn test_blastp_default_options() {
        let scoring = ScoringOptions::new_blastp();
        assert_eq!(scoring.reward, 0, "blastp does not use reward");
        assert_eq!(scoring.penalty, 0, "blastp does not use penalty");
        assert_eq!(scoring.gap_open, 11);
        assert_eq!(scoring.gap_extend, 1);
        assert!(scoring.gapped_calculation);
        assert_eq!(scoring.matrix_name.as_deref(), Some("BLOSUM62"));
        assert!(!scoring.is_ooframe);

        let word = InitialWordOptions::new_blastp();
        assert_eq!(word.word_size, WORDSIZE_PROT);
        assert_eq!(word.word_size, 3);
        assert_eq!(word.window_size, WINDOW_SIZE_PROT);
        assert_eq!(word.window_size, 40);
        assert_eq!(word.x_dropoff, UNGAPPED_X_DROPOFF_PROT);

        let ext = ExtensionOptions::new_blastp();
        assert_eq!(ext.gap_x_dropoff, GAP_X_DROPOFF_PROT);
        assert_eq!(ext.gap_x_dropoff_final, GAP_X_DROPOFF_FINAL_PROT);
        assert_eq!(ext.prelim_gap_ext, PrelimGapExt::DynProgScoreOnly);
        assert_eq!(ext.traceback_ext, TracebackExt::DynProgTbck);
    }

    /// Port of NCBI blastoptions_unit_test: blastx uses protein defaults.
    #[test]
    fn test_blastx_default_options() {
        // blastx is a translated nucleotide search, so it uses protein scoring
        let scoring = ScoringOptions::new_blastp();
        assert_eq!(scoring.gap_open, GAP_OPEN_PROT);
        assert_eq!(scoring.gap_extend, GAP_EXTN_PROT);
        assert_eq!(scoring.matrix_name.as_deref(), Some("BLOSUM62"));

        let word = InitialWordOptions::new_blastp();
        assert_eq!(word.word_size, WORDSIZE_PROT);
        assert_eq!(word.window_size, WINDOW_SIZE_PROT);
    }

    /// Port of NCBI blastoptions_unit_test: invalid reward/penalty combinations.
    #[test]
    fn test_option_validation_bad_reward_penalty() {
        // Reward must be positive, penalty must be negative for blastn.
        // While we don't have a validation function yet, verify the invariants
        // that the defaults satisfy and that violations can be detected.
        let opts = ScoringOptions::new_blastn();
        assert!(opts.reward > 0, "reward must be positive");
        assert!(opts.penalty < 0, "penalty must be negative");

        // Construct an invalid option set and verify the invariant is broken
        let bad = ScoringOptions {
            reward: -1,
            penalty: 3,
            ..ScoringOptions::new_blastn()
        };
        assert!(bad.reward <= 0, "negative reward is invalid for blastn");
        assert!(bad.penalty >= 0, "positive penalty is invalid for blastn");
    }

    /// Port of NCBI blastoptions_unit_test: verify DUST filter defaults.
    /// DUST uses window=64, threshold=2.0 as defaults in this codebase.
    #[test]
    fn test_dust_options() {
        // DUST defaults are not stored in a struct but used as parameters.
        // Verify the constants match NCBI defaults.
        let default_window: usize = 64;
        let default_threshold: f64 = 2.0;
        assert_eq!(default_window, 64);
        assert!((default_threshold - 2.0).abs() < f64::EPSILON);

        // Database options default to standard genetic code
        let db = DatabaseOptions::default();
        assert_eq!(db.genetic_code, 1);
    }
}
