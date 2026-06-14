use std::sync::Arc;

use crate::algo::blast::api::query_data::{CSeqLoc, EBlastProgramType, TMaskedQueryRegions};

#[derive(Clone, Debug)]
pub struct CDustFilter {
    pub query_loc: Option<Arc<CSeqLoc>>,
    pub program: EBlastProgramType,
    pub level: i32,
    pub window: i32,
    pub linker: i32,
    pub masks: TMaskedQueryRegions,
}
