use std::sync::Arc;

use crate::algo::blast::api::blast_options_builder::CBlastOptionsBuilder;
use crate::algo::blast::api::blast_options_handle::CBlastOptionsHandle;
use crate::algo::blast::api::query_data::{
    CBioseqSet, CPssmWithParameters, CSeqAlignSet, CSeqId, CSeqLoc,
};
use crate::algo::blast::api::remote_blast::CBlast4QueueSearchRequest;
use crate::algo::blast::core::blast_util::{SSeqRange, SubjectMaskingType};

#[derive(Clone, Debug)]
pub struct CImportStrategyData {
    pub options_handle: Option<Arc<CBlastOptionsHandle>>,
    pub filtering_id: i32,
    pub query_range: SSeqRange,
    pub task: String,
    pub psi_num_of_iterations: u32,
    pub filtering_key: String,
    pub subject_masking_type: SubjectMaskingType,
    pub pssm: Option<Arc<CPssmWithParameters>>,
    pub bioseq_set: Option<Arc<CBioseqSet>>,
    pub seq_loc: Option<Arc<CSeqLoc>>,
    pub seq_id: Option<Arc<CSeqId>>,
    pub seq_align_set: Option<Arc<CSeqAlignSet>>,
}

#[derive(Clone, Debug)]
pub struct CExportStrategy {
    pub data: Option<Arc<CImportStrategyData>>,
    pub queue_search_request: Option<Arc<CBlast4QueueSearchRequest>>,
    pub client_id: String,
    pub psi_num_iterations: u32,
}

#[derive(Clone, Debug)]
pub struct CImportStrategy {
    pub data: Option<Arc<CImportStrategyData>>,
    pub request: Option<Arc<CBlast4Request>>,
    pub service: String,
    pub options_builder: Option<Arc<CBlastOptionsBuilder>>,
    pub ignore_unsupported_options: bool,
}

#[derive(Clone, Debug)]
pub struct CBlast4Request;
