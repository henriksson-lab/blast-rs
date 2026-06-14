use std::sync::Arc;

use super::blast_options_cxx::CBlastOptions;
use super::local_db_adapter::CLocalDbAdapter;
use super::prelim_stage::CBlastPrelimSearch;
use super::query_data::{IQueryFactory, TSearchMessages};
use super::seqinfosrc_bioseq::IBlastSeqInfoSrc;
use super::setup_factory::SInternalData;
use super::traceback_stage::CBlastTracebackSearch;

/// NCBI C++: `CLocalBlast` (`local_blast.hpp`).
#[derive(Clone, Debug)]
pub struct CLocalBlast {
    pub m_QueryFactory: Option<Arc<IQueryFactory>>,
    pub m_Opts: Option<Arc<CBlastOptions>>,
    pub m_InternalData: Option<Arc<SInternalData>>,
    pub m_PrelimSearch: Option<Arc<CBlastPrelimSearch>>,
    pub m_TbackSearch: Option<Arc<CBlastTracebackSearch>>,
    pub m_LocalDbAdapter: Option<Arc<CLocalDbAdapter>>,
    pub m_SeqInfoSrc: Option<Arc<IBlastSeqInfoSrc>>,
    pub m_Messages: TSearchMessages,
    pub m_batch_num_str: String,
}
