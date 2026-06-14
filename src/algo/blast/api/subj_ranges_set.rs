use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

use crate::algo::blast::api::seqsrc_seqdb::CSeqDBRangeList;

pub const K_HSP_EXPAND_SIZE: i32 = 1024;
pub const K_MIN_GAP: i32 = 1024;

#[derive(Debug, Clone)]
pub struct CSubjectRanges {
    pub query_oids: BTreeSet<i32>,
    pub offsets: CSeqDBRangeList,
}

#[derive(Debug, Clone)]
pub struct CSubjectRangesSet {
    pub subj_ranges: BTreeMap<i32, Arc<CSubjectRanges>>,
    pub expand_hsp: i32,
    pub min_gap: i32,
}
