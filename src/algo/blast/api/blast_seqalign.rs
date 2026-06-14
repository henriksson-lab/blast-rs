use std::sync::Arc;

use crate::algo::blast::api::blast_results::{EResultType, TSeqAlignVector};
use crate::algo::blast::api::query_data::{
    CSeqAlign, CSeqAlignSet, CSeqId, CSeqLoc, CSplicedSeg, CStdSeg, EBlastProgramType,
    ILocalQueryData, TSeqLocInfoVector,
};
use crate::algo::blast::api::seqinfosrc_bioseq::IBlastSeqInfoSrc;
use crate::algo::blast::core::blast_hits::{BlastHSPList, BlastHSPResults, BlastHitList};
use crate::algo::blast::core::spliced_hits::HSPChain;

#[derive(Debug, Clone)]
pub struct SeqAlignConversionContext {
    pub program: EBlastProgramType,
    pub hsp_results: *mut BlastHSPResults,
    pub hsp_list: *mut BlastHSPList,
    pub hit_list: *mut BlastHitList,
    pub query_id: Option<Arc<CSeqId>>,
    pub subject_id: Option<Arc<CSeqId>>,
    pub query_loc: Option<Arc<CSeqLoc>>,
    pub query_length: i32,
    pub subject_length: i32,
    pub is_ooframe: bool,
    pub gapped: bool,
    pub oof_mode: bool,
    pub result_type: EResultType,
    pub seqid_list: Vec<String>,
    pub seq_aligns: Vec<Arc<CSeqAlign>>,
    pub seq_align_set: Option<Arc<CSeqAlignSet>>,
    pub std_seg_list: Vec<Arc<CStdSeg>>,
    pub subj_masks: Vec<TSeqLocInfoVector>,
    pub local_data: *mut ILocalQueryData,
    pub seqinfo_src: *const IBlastSeqInfoSrc,
    pub spliced_seg: Option<Arc<CSplicedSeg>>,
    pub hsp_chain: *const HSPChain,
}

pub type TLocalBlastSeqAlignVector = TSeqAlignVector;
