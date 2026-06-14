use crate::algo::blast::core::blast_lookup::LookupTableWrap;
use crate::algo::blast::core::blast_stat::BlastScoreBlk;

pub struct SPatternUnit {
    pub allowed_letters: String,
    pub disallowed_letters: String,
    pub at_least: usize,
    pub at_most: usize,
    pub is_x: bool,
}

pub struct CSeedTop {
    pub pattern: String,
    pub lookup: LookupTableWrap,
    pub score_blk: BlastScoreBlk,
    pub units: Vec<SPatternUnit>,
}
