//! Rust equivalent of blast_options.c — BLAST search options.
//! Contains all the option types and their default values.

use crate::program;

// ---- Default constants ----

// Window sizes
pub const WINDOW_SIZE_PROT: i32 = 40;
pub const WINDOW_SIZE_NUCL: i32 = 0;
pub const WINDOW_SIZE_MEGABLAST: i32 = 0;

// Word sizes
pub const WORDSIZE_PROT: i32 = 3;
pub const WORDSIZE_NUCL: i32 = 11;
pub const WORDSIZE_MEGABLAST: i32 = 28;

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

// E-value
pub const EXPECT_VALUE: f64 = 10.0;
pub const HITLIST_SIZE: i32 = 500;

// ---- Option structs ----

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

/// Hit saving options.
#[derive(Debug, Clone)]
pub struct HitSavingOptions {
    pub expect_value: f64,
    pub hitlist_size: i32,
    pub cutoff_score: i32,
    pub percent_identity: f64,
    pub min_hit_length: i32,
}

impl Default for HitSavingOptions {
    fn default() -> Self {
        HitSavingOptions {
            expect_value: EXPECT_VALUE,
            hitlist_size: HITLIST_SIZE,
            cutoff_score: 0,
            percent_identity: 0.0,
            min_hit_length: 0,
        }
    }
}

/// Initial word finding options.
#[derive(Debug, Clone)]
pub struct InitialWordOptions {
    pub window_size: i32,
    pub x_dropoff: f64,
    pub word_size: i32,
}

impl InitialWordOptions {
    pub fn new_blastn() -> Self {
        InitialWordOptions {
            window_size: WINDOW_SIZE_NUCL,
            x_dropoff: UNGAPPED_X_DROPOFF_NUCL,
            word_size: WORDSIZE_NUCL,
        }
    }

    pub fn new_blastp() -> Self {
        InitialWordOptions {
            window_size: WINDOW_SIZE_PROT,
            x_dropoff: UNGAPPED_X_DROPOFF_PROT,
            word_size: WORDSIZE_PROT,
        }
    }
}

/// Extension options.
#[derive(Debug, Clone)]
pub struct ExtensionOptions {
    pub gap_x_dropoff: f64,
    pub gap_x_dropoff_final: f64,
    pub gap_trigger: f64,
    pub prelim_gap_ext: PrelimGapExt,
    pub traceback_ext: TracebackExt,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PrelimGapExt {
    DynProgScoreOnly,
    GreedyScoreOnly,
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
        }
    }

    pub fn new_blastp() -> Self {
        ExtensionOptions {
            gap_x_dropoff: GAP_X_DROPOFF_PROT,
            gap_x_dropoff_final: GAP_X_DROPOFF_FINAL_PROT,
            gap_trigger: 22.0,
            prelim_gap_ext: PrelimGapExt::DynProgScoreOnly,
            traceback_ext: TracebackExt::DynProgTbck,
        }
    }
}

/// Effective length options.
#[derive(Debug, Clone, Default)]
pub struct EffectiveLengthsOptions {
    pub db_length: i64,
    pub num_searchspaces: i32,
    pub searchsp_eff: Vec<i64>,
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
    }
}
