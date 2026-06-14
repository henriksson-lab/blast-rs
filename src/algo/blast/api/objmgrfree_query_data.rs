use std::sync::Arc;

use crate::algo::blast::api::query_data::{CBioseq, CBioseqSet, IQueryFactory};

#[derive(Debug, Clone)]
pub struct CObjMgrFreeQueryFactory {
    pub base: IQueryFactory,
    pub bioseq: Option<Arc<CBioseq>>,
    pub bioseqs: Option<Arc<CBioseqSet>>,
}
