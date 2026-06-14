use std::sync::Arc;

use crate::algo::blast::api::bioseq_extract_data_priv::{ENaStrand, ESentinelType};
use crate::algo::blast::api::query_data::{
    CSeqLoc, EBlastProgramType, IBlastQuerySource, TSearchMessages,
};
use crate::algo::blast::core::blast_encoding::EBlastEncoding;
use crate::algo::blast::core::blast_query_info::BlastQueryInfo;
use crate::algo::blast::core::blast_util::BLAST_SequenceBlk;

#[derive(Clone, Debug)]
pub struct CBlastSetupCxxState {
    pub query_source: Option<Arc<IBlastQuerySource>>,
    pub program: EBlastProgramType,
    pub strand: ENaStrand,
    pub encoding: EBlastEncoding,
    pub sentinel: ESentinelType,
    pub query_info: *mut BlastQueryInfo,
    pub sequence_blk: *mut BLAST_SequenceBlk,
    pub messages: TSearchMessages,
    pub seq_locs: Vec<Arc<CSeqLoc>>,
}
