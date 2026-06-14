use std::sync::Arc;

use super::blast_options_cxx::CBlastOptions;
use super::local_db_adapter::CLocalDbAdapter;
use super::query_data::{IQueryFactory, TSearchMessages, TSeqLocInfoVector};
use super::setup_factory::SInternalData;
use super::uniform_search::CSearchDatabase;

/// NCBI C++: `CBlastPrelimSearch` (`prelim_stage.hpp`).
#[derive(Clone, Debug)]
pub struct CBlastPrelimSearch {
    pub m_QueryFactory: Option<Arc<IQueryFactory>>,
    pub m_InternalData: Option<Arc<SInternalData>>,
    pub m_Options: Option<Arc<CBlastOptions>>,
    pub m_DbAdapter: Option<Arc<CLocalDbAdapter>>,
    pub m_DbInfo: Option<Arc<CSearchDatabase>>,
    pub m_Messages: TSearchMessages,
    pub m_MasksForAllQueries: TSeqLocInfoVector,
}
