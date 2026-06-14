use std::sync::Arc;

use super::blast_options_handle::CBlastOptionsHandle;
use super::query_data::{IQueryFactory, TSeqLocVector};
use super::seqinfosrc_bioseq::IBlastSeqInfoSrc;
use super::uniform_search::CSearchDatabase;
use crate::algo::blast::core::blast_seqsrc::BlastSeqSrc;

/// NCBI C++: `CLocalDbAdapter` (`local_db_adapter.hpp`).
#[derive(Clone, Debug)]
pub struct CLocalDbAdapter {
    pub m_SeqSrc: Option<Arc<BlastSeqSrc>>,
    pub m_SeqInfoSrc: Option<Arc<IBlastSeqInfoSrc>>,
    pub m_DbInfo: Option<Arc<CSearchDatabase>>,
    pub m_SubjectFactory: Option<Arc<IQueryFactory>>,
    pub m_OptsHandle: Option<Arc<CBlastOptionsHandle>>,
    pub m_Subjects: TSeqLocVector,
    pub m_DbName: String,
    pub m_DbScanMode: bool,
}
