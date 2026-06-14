use std::sync::Arc;

use crate::algo::blast::api::query_data::{EBlastProgramType, IQueryFactory};
use crate::algo::blast::core::blast_query_info::BlastQueryInfo;

#[derive(Clone, Debug)]
pub struct CEffectiveSearchSpaceCalculator {
    pub query_factory: Option<Arc<IQueryFactory>>,
    pub program: EBlastProgramType,
    pub query_info: *mut BlastQueryInfo,
    pub db_num_seqs: i32,
    pub db_num_bases: i64,
}
