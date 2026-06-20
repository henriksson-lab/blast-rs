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

impl CRemoteSearchFactory {
    pub fn get_seq_search(&self) -> CRemoteSeqSearch {
        CRemoteSeqSearch {
            m_SearchOpts: None,
            m_RemoteBlast: None,
            m_Queries: None,
            m_Subject: None,
            m_Warnings: Vec::new(),
        }
    }

    pub fn get_pssm_search(&self) -> CRemotePssmSearch {
        CRemotePssmSearch {
            m_SearchOpts: None,
            m_RemoteBlast: None,
            m_Pssm: None,
            m_Subject: None,
            m_Warnings: Vec::new(),
        }
    }
}

#[derive(Clone, Debug)]
pub struct CRemoteSearchQueryFactory {
    pub m_QueryFactory: Option<Arc<IQueryFactory>>,
}
