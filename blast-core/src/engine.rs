//! Rust equivalent of blast_engine.c — top-level search orchestration.
//! This module will eventually replace Blast_RunFullSearch.

use crate::hspstream::{HspStream, HspResults, HspList, Hsp};
use crate::seqsrc::BlastSeqSource;
use crate::options::*;
use crate::stat::KarlinBlk;

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
}
