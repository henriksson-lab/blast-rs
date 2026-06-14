use std::sync::Arc;

use crate::algo::blast::api::query_data::{
    CBlastQueryVector, CScope, IQueryFactory, TSeqLocInfoVector, TSeqLocVector,
};

#[derive(Debug, Clone)]
pub struct CObjMgrQueryFactory {
    pub base: IQueryFactory,
    pub seq_loc_vector: TSeqLocVector,
    pub query_vector: Option<Arc<CBlastQueryVector>>,
    pub extracted_scopes: Vec<Arc<CScope>>,
    pub user_specified_masks: TSeqLocInfoVector,
}
