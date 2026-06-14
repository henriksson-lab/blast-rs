use crate::algo::blast::api::blast_results::CSearchResultSet;
use crate::algo::blast::api::blast_rps_options::CBlastRPSOptionsHandle;
use crate::algo::blast::api::deltablast_options::CDeltaBlastOptionsHandle;
use crate::algo::blast::api::local_db_adapter::CLocalDbAdapter;
use crate::algo::blast::api::query_data::{CPssmWithParameters, IQueryFactory};

pub struct CDeltaBlast {
    pub queries: Option<Box<IQueryFactory>>,
    pub subject: Option<Box<CLocalDbAdapter>>,
    pub domain_db: Option<Box<CLocalDbAdapter>>,
    pub options: Option<Box<CDeltaBlastOptionsHandle>>,
    pub rps_options: Option<Box<CBlastRPSOptionsHandle>>,
    pub pssm: Vec<Box<CPssmWithParameters>>,
    pub domain_results: Option<Box<CSearchResultSet>>,
    pub results: Option<Box<CSearchResultSet>>,
}
