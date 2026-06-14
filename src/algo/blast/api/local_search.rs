use std::sync::Arc;

use crate::algo::blast::api::blast_options_handle::CBlastOptionsHandle;
use crate::algo::blast::api::local_blast::CLocalBlast;
use crate::algo::blast::api::psiblast::CPsiBlast;
use crate::algo::blast::api::query_data::{CPssmWithParameters, IQueryFactory};
use crate::algo::blast::api::uniform_search::{CSearchDatabase, ISearch};

#[derive(Clone, Debug)]
pub struct CLocalSearch {
    pub base: ISearch,
    pub search_opts: Option<Arc<CBlastOptionsHandle>>,
    pub local_blast: Option<Arc<CLocalBlast>>,
    pub database: Option<Arc<CSearchDatabase>>,
    pub query_factory: Option<Arc<IQueryFactory>>,
    pub warnings: Vec<String>,
}

#[derive(Clone, Debug)]
pub struct CPsiSearch {
    pub base: ISearch,
    pub search_opts: Option<Arc<CBlastOptionsHandle>>,
    pub psi_blast: Option<Arc<CPsiBlast>>,
    pub pssm: Option<Arc<CPssmWithParameters>>,
    pub subject: Option<Arc<CSearchDatabase>>,
}
