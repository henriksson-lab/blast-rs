use std::sync::Arc;

use crate::algo::blast::api::query_data::{CBlastOptions, EBlastProgramType, ILocalQueryData};
use crate::algo::blast::api::setup_factory::SInternalData;
use crate::algo::blast::api::split_query_blk::CSplitQueryBlk;

#[derive(Debug, Clone)]
pub struct CContextTranslator {
    pub contexts_per_chunk: Vec<Vec<i32>>,
    pub starting_chunks: Vec<Vec<i32>>,
    pub absolute_contexts: Vec<Vec<i32>>,
    pub split_query_blk: Option<Arc<CSplitQueryBlk>>,
}

#[derive(Debug, Clone)]
pub struct CQueryDataPerChunk {
    pub program: EBlastProgramType,
    pub query_indices_per_chunk: Vec<Vec<usize>>,
    pub query_lengths: Vec<usize>,
    pub last_chunk_for_query_cache: Vec<i32>,
    pub local_query_data: Option<Arc<ILocalQueryData>>,
}

#[derive(Debug, Clone)]
pub struct SplitQueryAuxState {
    pub options: Option<Arc<CBlastOptions>>,
    pub full_data: Option<Arc<SInternalData>>,
    pub num_threaded: usize,
}
