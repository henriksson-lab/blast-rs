use std::sync::Arc;

use crate::algo::blast::api::query_data::{
    CBlastQueryVector, CScope, CSeqLoc, EBlastProgramType, TMaskedQueryRegions, TSeqLocVector,
};
use crate::algo::blast::core::blast_query_info::BlastQueryInfo;
use crate::algo::blast::core::blast_util::BLAST_SequenceBlk;

#[derive(Clone, Debug)]
pub struct BlastObjMgrToolsState {
    pub seq_locs: TSeqLocVector,
    pub query_vector: Option<Arc<CBlastQueryVector>>,
    pub program: EBlastProgramType,
    pub scope: Option<Arc<CScope>>,
    pub seq_loc: Option<Arc<CSeqLoc>>,
    pub masks: TMaskedQueryRegions,
    pub query_info: *mut BlastQueryInfo,
    pub sequence_blk: *mut BLAST_SequenceBlk,
}
