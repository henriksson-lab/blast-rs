use std::sync::Arc;

use crate::algo::blast::api::seqinfosrc_bioseq::IBlastSeqInfoSrc;
use crate::algo::blast::api::seqsrc_seqdb::CSeqDB;

#[derive(Debug, Clone)]
pub struct CSeqDbSeqInfoSrc {
    pub base: IBlastSeqInfoSrc,
    pub seqdb: Option<Arc<CSeqDB>>,
    pub filtering_algo_id: i32,
    pub dbname: String,
    pub is_protein: bool,
}
