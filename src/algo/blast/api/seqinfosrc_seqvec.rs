use crate::algo::blast::api::query_data::TSeqLocVector;
use crate::algo::blast::api::seqinfosrc_bioseq::IBlastSeqInfoSrc;

#[derive(Debug, Clone)]
pub struct CSeqVecSeqInfoSrc {
    pub base: IBlastSeqInfoSrc,
    pub seq_vec: TSeqLocVector,
}
