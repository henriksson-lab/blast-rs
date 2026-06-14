use std::sync::Arc;

use super::blast_options_handle::CBlastOptionsHandle;
use super::query_data::CPssmWithParameters;
use super::query_data::{IQueryFactory, IRemoteQueryData};
use super::remote_blast::CRemoteBlast;
use super::uniform_search::CSearchDatabase;

/// NCBI C++: `CRemoteSeqSearch` (`remote_search.hpp`).
#[derive(Clone, Debug)]
pub struct CRemoteSeqSearch {
    pub m_SearchOpts: Option<Arc<CBlastOptionsHandle>>,
    pub m_RemoteBlast: Option<Arc<CRemoteBlast>>,
    pub m_Queries: Option<Arc<IRemoteQueryData>>,
    pub m_Subject: Option<Arc<CSearchDatabase>>,
    pub m_Warnings: Vec<String>,
}

/// NCBI C++: `CRemotePssmSearch` (`remote_search.hpp`).
#[derive(Clone, Debug)]
pub struct CRemotePssmSearch {
    pub m_SearchOpts: Option<Arc<CBlastOptionsHandle>>,
    pub m_RemoteBlast: Option<Arc<CRemoteBlast>>,
    pub m_Pssm: Option<Arc<CPssmWithParameters>>,
    pub m_Subject: Option<Arc<CSearchDatabase>>,
    pub m_Warnings: Vec<String>,
}

/// NCBI C++: `CRemoteSearchFactory` (`remote_search.hpp`).
#[derive(Clone, Debug)]
pub struct CRemoteSearchFactory;

#[derive(Clone, Debug)]
pub struct CRemoteSearchQueryFactory {
    pub m_QueryFactory: Option<Arc<IQueryFactory>>,
}
