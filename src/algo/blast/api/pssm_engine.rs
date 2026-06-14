use crate::algo::blast::core::blast_stat::BlastScoreBlk;

pub struct IPssmInputData;
pub struct IPssmInputFreqRatios;
pub struct IPssmInputCdd;

pub enum CPssmEngineErrCode {
    NullInputData,
    InvalidInputData,
}

pub struct CPssmEngineException {
    pub err_code: CPssmEngineErrCode,
}

pub struct CPssmEngine {
    pub pssm_input: *mut IPssmInputData,
    pub pssm_input_freq_ratios: *mut IPssmInputFreqRatios,
    pub score_blk: BlastScoreBlk,
    pub pssm_input_cdd: *mut IPssmInputCdd,
}

pub struct CScorematPssmConverter;
