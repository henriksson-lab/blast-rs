//! Rust equivalent of blast_parameters.c — computed search parameters.
//! These are derived from options + scoring blocks at search time.

use crate::options::*;

/// Scoring parameters computed from options and score block.
#[derive(Debug, Clone)]
pub struct ScoringParameters {
    pub options: ScoringOptions,
    pub reward: i32,
    pub penalty: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub scale_factor: f64,
}

impl ScoringParameters {
    pub fn from_options(opts: &ScoringOptions, scale_factor: f64) -> Self {
        ScoringParameters {
            options: opts.clone(),
            reward: opts.reward,
            penalty: opts.penalty,
            gap_open: opts.gap_open,
            gap_extend: opts.gap_extend,
            scale_factor,
        }
    }
}

/// Extension parameters computed from options and score block.
#[derive(Debug, Clone)]
pub struct ExtensionParameters {
    pub options: ExtensionOptions,
    pub gap_x_dropoff: i32,       // raw score x-dropoff for preliminary
    pub gap_x_dropoff_final: i32, // raw score x-dropoff for traceback
    pub gap_trigger: i32,         // raw score trigger for gapped extension
}

/// Hit saving parameters computed from options.
#[derive(Debug, Clone)]
pub struct HitSavingParameters {
    pub options: HitSavingOptions,
    pub cutoff_score_min: i32,
    pub low_score: Vec<i32>,
}

/// Initial word parameters computed from options and score block.
#[derive(Debug, Clone)]
pub struct InitialWordParameters {
    pub options: InitialWordOptions,
    pub x_dropoff_max: i32,
    pub cutoff_score_min: i32,
    /// Scoring table for 4-base nucleotide words (256 entries).
    pub nucl_score_table: [i32; 256],
}

impl InitialWordParameters {
    /// Build the nucleotide score table for 4-base packed words.
    /// Each byte encodes 4 bases (2 bits each), and the table stores
    /// the combined match/mismatch score for all 256 possible byte values.
    /// Build the nucleotide score table for 4-base packed words.
    /// Indexed by XOR of query and subject packed bytes.
    /// Each bit pair in the XOR result: 00 = match (reward), != 00 = mismatch (penalty).
    pub fn build_nucl_score_table(reward: i32, penalty: i32) -> [i32; 256] {
        let mut table = [0i32; 256];
        for xor_val in 0..256u32 {
            let mut score = 0i32;
            for pos in 0..4 {
                let bits = (xor_val >> (6 - 2 * pos)) & 3;
                score += if bits == 0 { reward } else { penalty };
            }
            table[xor_val as usize] = score;
        }
        table
    }
}

/// Effective length parameters.
#[derive(Debug, Clone)]
pub struct EffectiveLengthsParameters {
    pub options: EffectiveLengthsOptions,
    pub real_db_length: i64,
    pub real_num_seqs: i32,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_scoring_params_from_options() {
        let opts = ScoringOptions::new_blastn();
        let params = ScoringParameters::from_options(&opts, 1.0);
        assert_eq!(params.reward, 1);
        assert_eq!(params.penalty, -3);
        assert_eq!(params.gap_open, 5);
        assert_eq!(params.gap_extend, 2);
        assert_eq!(params.scale_factor, 1.0);
    }

    /// Port of NCBI blastsetup_unit_test: derive parameters from blastn options.
    #[test]
    fn test_parameters_from_blastn_options() {
        let scoring_opts = ScoringOptions::new_blastn();
        let params = ScoringParameters::from_options(&scoring_opts, 1.0);

        // Parameters should faithfully carry over scoring options
        assert_eq!(params.reward, REWARD);
        assert_eq!(params.penalty, PENALTY);
        assert_eq!(params.gap_open, GAP_OPEN_NUCL);
        assert_eq!(params.gap_extend, GAP_EXTN_NUCL);
        assert_eq!(params.scale_factor, 1.0);

        // The stored options should match
        assert_eq!(params.options.reward, params.reward);
        assert_eq!(params.options.penalty, params.penalty);
        assert_eq!(params.options.gap_open, params.gap_open);
        assert_eq!(params.options.gap_extend, params.gap_extend);

        // Verify nucl score table for reward=1, penalty=-3
        let table = InitialWordParameters::build_nucl_score_table(
            scoring_opts.reward,
            scoring_opts.penalty,
        );
        // XOR 0x00 = all 4 bases match -> 4 * reward = 4
        assert_eq!(table[0x00], 4 * scoring_opts.reward);
        // XOR 0xFF = all 4 bases mismatch -> 4 * penalty = -12
        assert_eq!(table[0xFF], 4 * scoring_opts.penalty);
        // XOR 0x01 = 3 match + 1 mismatch -> 3*1 + 1*(-3) = 0
        assert_eq!(table[0x01], 3 * scoring_opts.reward + scoring_opts.penalty);
    }

    /// Port of NCBI blastsetup_unit_test: derive parameters from blastp options.
    #[test]
    fn test_parameters_from_blastp_options() {
        let scoring_opts = ScoringOptions::new_blastp();
        let params = ScoringParameters::from_options(&scoring_opts, 1.0);

        assert_eq!(params.gap_open, GAP_OPEN_PROT);
        assert_eq!(params.gap_extend, GAP_EXTN_PROT);
        assert_eq!(params.options.matrix_name.as_deref(), Some("BLOSUM62"));
        assert!(params.options.gapped_calculation);

        // With a non-default scale factor
        let scaled = ScoringParameters::from_options(&scoring_opts, 2.5);
        assert_eq!(scaled.scale_factor, 2.5);
        // Gap costs remain unscaled in the parameters struct
        assert_eq!(scaled.gap_open, GAP_OPEN_PROT);
        assert_eq!(scaled.gap_extend, GAP_EXTN_PROT);
    }
}
