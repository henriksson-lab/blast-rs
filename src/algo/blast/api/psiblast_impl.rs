use crate::algo::blast::api::blast_options_handle::CBlastOptionsHandle;
use crate::algo::blast::api::blast_results::CSearchResultSet;
use crate::algo::blast::api::local_db_adapter::CLocalDbAdapter;
use crate::algo::blast::api::query_data::{CPssmWithParameters, IQueryFactory};

pub enum EResultType {
    Alignment,
    SearchResultSet,
}

pub struct CPsiBlastImpl {
    pub pssm: Option<Box<CPssmWithParameters>>,
    pub query: Option<Box<IQueryFactory>>,
    pub subject: Option<Box<CLocalDbAdapter>>,
    pub opts_handle: Option<Box<CBlastOptionsHandle>>,
    pub results: Option<Box<CSearchResultSet>>,
    pub result_type: EResultType,
}
