use std::sync::Arc;

use crate::algo::blast::api::bioseq_extract_data_priv::{ENaStrand, ESentinelType};
use crate::algo::blast::api::query_data::{
    CBlastOptions, CBlastQueryVector, CScope, CSeqLoc, EBlastProgramType, IBlastQuerySource,
    TMaskedQueryRegions, TSeqLocVector,
};
use crate::algo::blast::core::blast_encoding::EBlastEncoding;
use crate::algo::blast::core::blast_query_info::BlastQueryInfo;
use crate::algo::blast::core::blast_util::BLAST_SequenceBlk;

#[derive(Debug, Clone)]
pub struct CBlastQuerySourceOM {
    pub base: IBlastQuerySource,
    pub query_vector: Option<Arc<CBlastQueryVector>>,
    pub seq_loc_vector: Option<TSeqLocVector>,
    pub own_seq_loc_vector: bool,
    pub options: *const CBlastOptions,
    pub calculated_masks: bool,
    pub program: EBlastProgramType,
    pub strand: ENaStrand,
    pub encoding: EBlastEncoding,
    pub sentinel: ESentinelType,
    pub masks: Vec<TMaskedQueryRegions>,
    pub query_info: *mut BlastQueryInfo,
    pub sequence_blk: *mut BLAST_SequenceBlk,
    pub packed_seqints: Vec<Arc<CSeqLoc>>,
    pub scopes: Vec<Arc<CScope>>,
}
