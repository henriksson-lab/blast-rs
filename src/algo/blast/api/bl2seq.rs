use std::sync::Arc;

use crate::algo::blast::api::blast_options_handle::CBlastOptionsHandle;
use crate::algo::blast::api::blast_results::{CSearchResultSet, TAncillaryVector, TSeqAlignVector};
use crate::algo::blast::api::local_blast::CLocalBlast;
use crate::algo::blast::api::query_data::{TSearchMessages, TSeqLocVector};
use crate::algo::blast::api::setup_factory::TInterruptFnPtr;

#[derive(Clone, Debug)]
pub struct CBl2Seq {
    pub queries: TSeqLocVector,
    pub subjects: TSeqLocVector,
    pub opts_handle: Option<Arc<CBlastOptionsHandle>>,
    pub blast: Option<Arc<CLocalBlast>>,
    pub dbscan_mode: bool,
    pub messages: TSearchMessages,
    pub ancillary_data: TAncillaryVector,
    pub results: Option<Arc<CSearchResultSet>>,
    pub interrupt_fn: Option<TInterruptFnPtr>,
    pub interrupt_user_data: usize,
    pub seq_aligns: TSeqAlignVector,
}

#[derive(Clone, Debug)]
pub struct CBlastFilterTest {
    pub messages: TSearchMessages,
}
