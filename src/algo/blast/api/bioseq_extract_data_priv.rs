use std::sync::Arc;

use crate::algo::blast::api::query_data::{
    CBioseq, CBioseqSet, CSeqData, CSeqId, CSeqLoc, EBlastProgramType, TMaskedQueryRegions,
};
use crate::algo::blast::core::blast_encoding::EBlastEncoding;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CSeqUtilCoding {
    NotSet,
    Iupacna,
    Iupacaa,
    Ncbi2na,
    Ncbi4na,
    Ncbieaa,
    Ncbistdaa,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CSeqDataChoice {
    NotSet,
    Iupacna,
    Iupacaa,
    Ncbi2na,
    Ncbi4na,
    Ncbieaa,
    Ncbistdaa,
}

#[derive(Debug, Clone)]
pub struct SBlastSequence {
    pub data: Vec<u8>,
    pub length: u32,
    pub encoding: EBlastEncoding,
}

#[derive(Debug, Clone)]
pub struct CBlastSeqVectorFromCSeqData {
    pub seq_data: Option<Arc<CSeqData>>,
    pub length: u32,
    pub sequence_data: Vec<i8>,
    pub encoding: CSeqUtilCoding,
    pub coding: CSeqDataChoice,
    pub plus_strand: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ENaStrand {
    Unknown,
    Plus,
    Minus,
    Both,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ESentinelType {
    NoSentinels,
    Sentinels,
}

#[derive(Debug, Clone)]
pub struct CBlastQuerySourceBioseqSet {
    pub is_prot: bool,
    pub bioseqs: Vec<Arc<CBioseq>>,
    pub bioseq_set: Option<Arc<CBioseqSet>>,
    pub masks: Vec<Option<Arc<CSeqLoc>>>,
    pub masked_regions: Vec<TMaskedQueryRegions>,
    pub seq_locs: Vec<Arc<CSeqLoc>>,
    pub seq_ids: Vec<Arc<CSeqId>>,
    pub genetic_code_ids: Vec<u32>,
    pub strands: Vec<ENaStrand>,
    pub titles: Vec<String>,
    pub segment_info: Vec<i32>,
    pub program: EBlastProgramType,
}
