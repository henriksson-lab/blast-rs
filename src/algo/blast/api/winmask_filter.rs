use std::sync::Arc;

use crate::algo::blast::api::query_data::{CSeqLoc, TMaskedQueryRegions};

#[derive(Clone, Debug)]
pub struct CWindowMaskerFilter {
    pub query_loc: Option<Arc<CSeqLoc>>,
    pub taxid: i32,
    pub database: String,
    pub masks: TMaskedQueryRegions,
}
