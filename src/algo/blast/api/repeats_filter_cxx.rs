use std::sync::Arc;

use crate::algo::blast::api::query_data::{CSeqLoc, TMaskedQueryRegions};

#[derive(Clone, Debug)]
pub struct CRepeatsFilter {
    pub query_loc: Option<Arc<CSeqLoc>>,
    pub database: String,
    pub masks: TMaskedQueryRegions,
}
