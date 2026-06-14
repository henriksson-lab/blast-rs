use std::sync::Arc;

use crate::algo::blast::api::blast_results::CSearchResultSet;
use crate::algo::blast::api::magicblast_options::CMagicBlastOptionsHandle;
use crate::algo::blast::api::query_data::{IQueryFactory, TMaskedQueryRegions};
use crate::algo::blast::api::uniform_search::CSearchDatabase;

#[derive(Clone, Debug)]
pub struct CMagicBlast {
    pub query_factory: Option<Arc<IQueryFactory>>,
    pub database: Option<Arc<CSearchDatabase>>,
    pub options: Option<Arc<CMagicBlastOptionsHandle>>,
    pub results: Option<Arc<CSearchResultSet>>,
    pub splice_junctions: Vec<String>,
    pub query_masks: Vec<TMaskedQueryRegions>,
}
