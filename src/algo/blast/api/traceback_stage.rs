use std::sync::Arc;

use super::blast_options_cxx::{CBlastOptions, CBlastOptionsMemento};
use super::query_data::{IQueryFactory, TSearchMessages};
use super::seqinfosrc_bioseq::IBlastSeqInfoSrc;
use super::setup_factory::{SDatabaseScanData, SInternalData};

/// NCBI C++: `CBlastTracebackSearch` (`traceback_stage.hpp`).
#[derive(Clone, Debug)]
pub struct CBlastTracebackSearch {
    pub m_QueryFactory: Option<Arc<IQueryFactory>>,
    pub m_Options: Option<Arc<CBlastOptions>>,
    pub m_InternalData: Option<Arc<SInternalData>>,
    pub m_OptsMemento: Option<Arc<CBlastOptionsMemento>>,
    pub m_Messages: TSearchMessages,
    pub m_SeqInfoSrc: Option<Arc<IBlastSeqInfoSrc>>,
    pub m_ResultType: EResultType,
    pub m_DBscanInfo: Option<Arc<SDatabaseScanData>>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum EResultType {
    EDatabaseSearch,
    EBl2Seq,
}
