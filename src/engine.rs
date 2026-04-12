//! Rust equivalent of blast_engine.c — top-level search orchestration.
//! This module will eventually replace Blast_RunFullSearch.

use crate::options::*;

/// Search configuration combining all options.
pub struct SearchConfig {
    pub scoring: ScoringOptions,
    pub hit_saving: HitSavingOptions,
    pub initial_word: InitialWordOptions,
    pub extension: ExtensionOptions,
    pub eff_lengths: EffectiveLengthsOptions,
    pub db: DatabaseOptions,
}

impl SearchConfig {
    /// Create default configuration for blastn.
    pub fn default_blastn() -> Self {
        SearchConfig {
            scoring: ScoringOptions::new_blastn(),
            hit_saving: HitSavingOptions::default(),
            initial_word: InitialWordOptions::new_blastn(),
            extension: ExtensionOptions::new_blastn(),
            eff_lengths: EffectiveLengthsOptions::default(),
            db: DatabaseOptions::default(),
        }
    }
}

/// Result of a BLAST search for one query.
#[derive(Debug)]
pub struct QueryResult {
    pub query_id: String,
    pub hits: Vec<AlignmentHit>,
}

/// A single alignment hit.
#[derive(Debug)]
pub struct AlignmentHit {
    pub subject_oid: i32,
    pub score: i32,
    pub bit_score: f64,
    pub evalue: f64,
    pub query_start: i32,
    pub query_end: i32,
    pub subject_start: i32,
    pub subject_end: i32,
    pub num_ident: i32,
    pub align_length: i32,
    pub num_gaps: i32,
    pub context: i32, // 0=plus, 1=minus
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_config() {
        let config = SearchConfig::default_blastn();
        assert_eq!(config.scoring.reward, 1);
        assert_eq!(config.scoring.penalty, -3);
        assert_eq!(config.initial_word.window_size, 0);
        assert_eq!(config.hit_saving.expect_value, 10.0);
    }

    #[test]
    fn test_default_config_extension_options() {
        let config = SearchConfig::default_blastn();
        // Extension options should have blastn defaults
        assert_eq!(config.extension.gap_x_dropoff, 30.0);
        assert_eq!(config.extension.gap_x_dropoff_final, 100.0);
    }

    #[test]
    fn test_default_config_db_options() {
        let config = SearchConfig::default_blastn();
        assert_eq!(config.db.genetic_code, 1);
    }

    #[test]
    fn test_default_config_eff_lengths() {
        let config = SearchConfig::default_blastn();
        assert_eq!(config.eff_lengths.db_length, 0);
    }

    #[test]
    fn test_alignment_hit_debug() {
        let hit = AlignmentHit {
            subject_oid: 0,
            score: 100,
            bit_score: 50.5,
            evalue: 1e-10,
            query_start: 1,
            query_end: 50,
            subject_start: 10,
            subject_end: 60,
            num_ident: 45,
            align_length: 50,
            num_gaps: 2,
            context: 0,
        };
        let debug_str = format!("{:?}", hit);
        assert!(debug_str.contains("score: 100"));
        assert!(debug_str.contains("evalue"));
    }

    #[test]
    fn test_query_result_debug() {
        let result = QueryResult {
            query_id: "test_query".to_string(),
            hits: vec![],
        };
        let debug_str = format!("{:?}", result);
        assert!(debug_str.contains("test_query"));
        assert!(debug_str.contains("hits"));
    }

    #[test]
    fn test_alignment_hit_context_plus_minus() {
        let plus_hit = AlignmentHit {
            subject_oid: 1,
            score: 50,
            bit_score: 25.0,
            evalue: 0.01,
            query_start: 1,
            query_end: 20,
            subject_start: 1,
            subject_end: 20,
            num_ident: 18,
            align_length: 20,
            num_gaps: 0,
            context: 0,
        };
        let minus_hit = AlignmentHit {
            subject_oid: 1,
            score: 50,
            bit_score: 25.0,
            evalue: 0.01,
            query_start: 1,
            query_end: 20,
            subject_start: 1,
            subject_end: 20,
            num_ident: 18,
            align_length: 20,
            num_gaps: 0,
            context: 1,
        };
        assert_eq!(plus_hit.context, 0);
        assert_eq!(minus_hit.context, 1);
    }

    #[test]
    fn test_query_result_with_hits() {
        let result = QueryResult {
            query_id: "seq1".to_string(),
            hits: vec![
                AlignmentHit {
                    subject_oid: 0,
                    score: 200,
                    bit_score: 100.0,
                    evalue: 1e-30,
                    query_start: 1,
                    query_end: 100,
                    subject_start: 50,
                    subject_end: 150,
                    num_ident: 90,
                    align_length: 100,
                    num_gaps: 5,
                    context: 0,
                },
                AlignmentHit {
                    subject_oid: 1,
                    score: 150,
                    bit_score: 75.0,
                    evalue: 1e-20,
                    query_start: 10,
                    query_end: 80,
                    subject_start: 200,
                    subject_end: 270,
                    num_ident: 60,
                    align_length: 70,
                    num_gaps: 3,
                    context: 1,
                },
            ],
        };
        assert_eq!(result.hits.len(), 2);
        assert!(result.hits[0].score > result.hits[1].score);
        assert!(result.hits[0].evalue < result.hits[1].evalue);
    }
}
