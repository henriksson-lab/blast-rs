use std::sync::Arc;

use crate::algo::blast::api::bioseq_extract_data_priv::CBlastQuerySourceBioseqSet;
use crate::algo::blast::api::query_data::{
    CBioseq, CBioseqSet, CSeqId, CSeqLoc, TMaskedQueryRegions,
};
use crate::algo::blast::core::blast_util::SSeqRange;

pub type TMaskedSubjRegions = TMaskedQueryRegions;

#[derive(Debug, Clone)]
pub struct IBlastSeqInfoSrc {
    pub ids: Vec<Vec<Arc<CSeqId>>>,
    pub seq_locs: Vec<Option<Arc<CSeqLoc>>>,
    pub lengths: Vec<u32>,
    pub has_gi_list: bool,
    pub masks: Vec<TMaskedSubjRegions>,
    pub can_return_partial_sequence: bool,
}

#[derive(Debug, Clone)]
pub struct CBioseqSeqInfoSrc {
    pub data_source: CBlastQuerySourceBioseqSet,
    pub bioseq: Option<Arc<CBioseq>>,
    pub bioseq_set: Option<Arc<CBioseqSet>>,
    pub is_prot: bool,
    pub target_ranges: Vec<SSeqRange>,
}
