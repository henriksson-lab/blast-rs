use std::sync::Arc;

use crate::algo::blast::api::query_data::{
    CSeqId, CSeqLoc, CSeqLocInfo, EBlastProgramType, TMaskedQueryRegions,
};
use crate::algo::blast::core::blast_filter::BlastMaskLoc;

#[derive(Clone, Debug)]
pub struct CFrameFinder {
    pub frame: i32,
}

#[derive(Clone, Debug)]
pub struct CAutomaticGenCodeSingleton {
    pub ref_counter: u32,
    pub genetic_codes: Vec<i32>,
}

#[derive(Clone, Debug)]
pub struct CBlastAppDiagHandler {
    pub saved_messages: Vec<String>,
    pub save: bool,
}

#[derive(Clone, Debug)]
pub struct CConstRefCSeqIdLessThan {
    pub lhs: Option<Arc<CSeqId>>,
    pub rhs: Option<Arc<CSeqId>>,
}

#[derive(Clone, Debug)]
pub struct BlastAuxConversionState {
    pub program: EBlastProgramType,
    pub seq_loc: Option<Arc<CSeqLoc>>,
    pub seq_ids: Vec<Arc<CSeqId>>,
    pub masked_regions: TMaskedQueryRegions,
    pub seq_loc_infos: Vec<Arc<CSeqLocInfo>>,
    pub mask: *const BlastMaskLoc,
}
