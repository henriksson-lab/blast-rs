use std::sync::Arc;

use crate::algo::blast::api::blast_results::{CSearchResultSet, TAncillaryVector, TSeqAlignVector};
use crate::algo::blast::api::query_data::{
    CBlastOptions, CPssmWithParameters, CSearchMessage, CSeqId, CSeqLoc, EBlastProgramType,
    IQueryFactory, TSearchMessages, TSeqLocInfoVector,
};
use crate::algo::blast::api::setup_factory::SInternalData;
use crate::algo::blast::api::split_query_cxx::CQuerySplitter;
use crate::algo::blast::core::blast_query_info::BlastQueryInfo;
use crate::algo::blast::core::blast_seqsrc::BlastSeqSrc;

#[derive(Debug, Clone)]
pub struct BlastMessage;

#[derive(Debug, Clone)]
pub struct SBlastSetupData {
    pub internal_data: Option<Arc<SInternalData>>,
    pub query_splitter: Option<Arc<CQuerySplitter>>,
    pub masks: TSeqLocInfoVector,
    pub messages: TSearchMessages,
    pub query_factory: Option<Arc<IQueryFactory>>,
    pub options: Option<Arc<CBlastOptions>>,
    pub pssm: Option<Arc<CPssmWithParameters>>,
    pub seqsrc: *mut BlastSeqSrc,
    pub num_threads: usize,
}

#[derive(Debug, Clone)]
pub struct BlastAuxPrivState {
    pub seqids: Vec<Arc<CSeqId>>,
    pub whole_seq_loc: Option<Arc<CSeqLoc>>,
    pub search_messages: TSearchMessages,
    pub blast_messages: Vec<Arc<BlastMessage>>,
    pub query_info: *const BlastQueryInfo,
    pub program: EBlastProgramType,
    pub alignments: TSeqAlignVector,
    pub ancillary_data: TAncillaryVector,
    pub result_set: Option<Arc<CSearchResultSet>>,
    pub subject_masks: Vec<TSeqLocInfoVector>,
    pub query_masks: TSeqLocInfoVector,
    pub raw_messages: Vec<Arc<CSearchMessage>>,
}
