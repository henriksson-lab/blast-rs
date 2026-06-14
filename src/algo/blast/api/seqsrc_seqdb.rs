use std::any::Any;
use std::collections::BTreeSet;
use std::sync::Arc;

use crate::algo::blast::core::blast_seqsrc::BlastSeqSrcGetSeqArg;
use crate::algo::blast::core::blast_util::SubjectMaskingType;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CSeqDBOidListType {
    OidList,
    OidRange,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CSeqDBSeqType {
    Protein,
    Nucleotide,
    Unknown,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CSeqDBSummaryType {
    UnfilteredAll,
    FilteredAll,
    FilteredRange,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CSeqDBMmapFileType {
    IndexFile,
    SequenceFile,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CSeqDBMmapStrategy {
    Normal,
    Sequential,
    WillNeed,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CSeqDBOffsetPair {
    pub first: u32,
    pub second: u32,
}

#[derive(Debug, Clone)]
pub struct CSeqDBSequenceRanges {
    pub size: usize,
    pub capacity: usize,
    pub data: Vec<CSeqDBOffsetPair>,
}

pub type CSeqDBRangeList = BTreeSet<(i32, i32)>;

#[derive(Debug, Clone)]
pub struct CSeqDB {
    pub impl_ptr: *mut CSeqDBImpl,
}

#[derive(Debug, Clone)]
pub struct CSeqDBExpert {
    pub base: CSeqDB,
}

#[derive(Debug, Clone)]
pub struct CSeqDBImpl;

#[derive(Debug, Clone)]
pub struct SSeqDBSeqSrcData {
    pub seqdb: Option<Arc<CSeqDBExpert>>,
    pub mask_algo_id: i32,
    pub mask_type: SubjectMaskingType,
    pub copied: bool,
    pub is_protein: bool,
    pub seq_ranges: CSeqDBSequenceRanges,
}

pub type TSeqDBData = SSeqDBSeqSrcData;

pub fn s_seq_db_get_is_prot(seqdb_handle: &dyn Any, _args: Option<&mut dyn Any>) -> bool {
    seqdb_handle
        .downcast_ref::<TSeqDBData>()
        .map(|datap| datap.is_protein)
        .unwrap_or(false)
}

pub fn s_seq_db_release_sequence(seqdb_handle: &dyn Any, args: &mut BlastSeqSrcGetSeqArg) {
    let Some(datap) = seqdb_handle.downcast_ref::<TSeqDBData>() else {
        return;
    };
    if args.seq.is_null() {
        return;
    }

    unsafe {
        if (*args.seq).sequence_start_allocated {
            if datap.copied {
                (*args.seq).sequence_start = None;
            }
            (*args.seq).sequence_start_allocated = false;
        }
        if (*args.seq).sequence_allocated {
            (*args.seq).sequence = None;
            (*args.seq).sequence_allocated = false;
        }
    }

    args.ranges = None;
}
