use std::sync::Arc;

use crate::algo::blast::api::query_data::{
    CBlastOptions, CBlastQueryVector, CScope, ILocalQueryData, IQueryFactory, TSeqLocInfoVector,
};
use crate::algo::blast::api::split_query_blk::CSplitQueryBlk;

#[derive(Debug, Clone)]
pub struct CQuerySplitter {
    pub query_factory: Option<Arc<IQueryFactory>>,
    pub options: *const CBlastOptions,
    pub num_chunks: u32,
    pub split_blk: Option<Arc<CSplitQueryBlk>>,
    pub query_chunk_factories: Vec<Arc<IQueryFactory>>,
    pub local_query_data: Option<Arc<ILocalQueryData>>,
    pub total_query_length: usize,
    pub chunk_size: usize,
    pub scopes: Vec<Arc<CScope>>,
    pub user_specified_masks: TSeqLocInfoVector,
    pub split_queries_in_chunk: Vec<Arc<CBlastQueryVector>>,
}
