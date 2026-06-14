use std::sync::Arc;

use crate::algo::blast::api::query_data::{
    CSeqAlignSet, CSeqId, EBlastProgramType, TMaskedQueryRegions, TQueryMessages, TSearchMessages,
    TSeqLocInfoVector,
};
use crate::algo::blast::core::blast_query_info::SPHIQueryInfo;
use crate::algo::blast::core::blast_stat::{GumbelBlk, KarlinBlk};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EResultType {
    DatabaseSearch,
    SequenceComparison,
}

pub type TSeqAlignVector = Vec<Option<Arc<CSeqAlignSet>>>;

#[derive(Debug, Clone)]
pub struct CBlastAncillaryData {
    pub program_type: EBlastProgramType,
    pub gumbel_blk: *mut GumbelBlk,
    pub ungapped_karlin_blk: *mut KarlinBlk,
    pub gapped_karlin_blk: *mut KarlinBlk,
    pub psi_ungapped_karlin_blk: *mut KarlinBlk,
    pub psi_gapped_karlin_blk: *mut KarlinBlk,
    pub search_space: i64,
    pub length_adjustment: i64,
}

#[derive(Debug, Clone)]
pub struct CSearchResults {
    pub query_id: Option<Arc<CSeqId>>,
    pub alignment: Option<Arc<CSeqAlignSet>>,
    pub errors: TQueryMessages,
    pub masks: TMaskedQueryRegions,
    pub subject_masks: TSeqLocInfoVector,
    pub ancillary_data: Option<Arc<CBlastAncillaryData>>,
    pub rid: String,
    pub phi_query_info: *mut SPHIQueryInfo,
}

pub type TQueryIdVector = Vec<Arc<CSeqId>>;
pub type TAncillaryVector = Vec<Arc<CBlastAncillaryData>>;

#[derive(Debug, Clone)]
pub struct CSearchResultSet {
    pub result_type: EResultType,
    pub num_queries: usize,
    pub results: Vec<Arc<CSearchResults>>,
    pub is_phi_blast: bool,
    pub query_masks: TSeqLocInfoVector,
    pub alignments: TSeqAlignVector,
    pub messages: TSearchMessages,
    pub query_ids: TQueryIdVector,
    pub ancillary_data: TAncillaryVector,
}
