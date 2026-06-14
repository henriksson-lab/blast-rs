use std::sync::Arc;

use crate::algo::blast::api::query_data::CSeqId;
use crate::algo::blast::api::seqinfosrc_bioseq::IBlastSeqInfoSrc;

pub type TGi = i64;

#[derive(Clone, Debug)]
pub struct BlastSeqInfoSrcAuxState {
    pub seqinfo_src: Option<Arc<IBlastSeqInfoSrc>>,
    pub oid: i32,
    pub seqid: Option<Arc<CSeqId>>,
    pub length: u32,
    pub gis: Vec<TGi>,
    pub seqids: Vec<String>,
    pub use_gis: bool,
}
