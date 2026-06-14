use std::any::Any;
use std::cell::Cell;
use std::ffi::c_char;
use std::ffi::CStr;

use crate::algo::blast::core::blast_seqsrc::{BlastSeqSrcGetSeqArg, BlastSeqSrcSetRangesArg};

pub const K_EMPTY_VDB_NAME: &str = "";

#[derive(Debug, Clone)]
pub struct VdbBlast2naReader;

#[derive(Debug, Clone)]
pub struct Packed2naRead;

#[derive(Debug, Clone)]
pub struct VdbBlastMgr;

#[derive(Debug, Clone)]
pub struct VdbBlastRunSet;

#[derive(Debug, Clone)]
pub struct VdbBlastReferenceSet;

#[derive(Debug, Clone)]
pub struct VdbBlast4naReader;

#[derive(Debug, Clone)]
pub struct SVDBSRCSeqSrc2naIterReader {
    pub reader: *mut VdbBlast2naReader,
    pub buffer: *mut Packed2naRead,
    pub max_index: u32,
    pub current_index: u32,
    pub chunked_buf: *mut u8,
    pub chunked_buf_size: u32,
}

pub type TVDB2naICReader = SVDBSRCSeqSrc2naIterReader;

#[derive(Debug, Clone)]
pub struct SVDBSRCPartialFetchingRanges {
    pub oid: i32,
    pub num_ranges: i32,
    pub ranges: *mut i32,
}

pub type TVDBPartialFetchingRanges = SVDBSRCPartialFetchingRanges;

#[derive(Debug, Clone)]
pub struct SVBDSRCSeqSrcData {
    pub hdl_mgr: *mut VdbBlastMgr,
    pub num_runs: u32,
    pub num_seqs: u64,
    pub run_set: *mut VdbBlastRunSet,
    pub ref_set: *mut VdbBlastReferenceSet,
    pub names: *mut c_char,
    pub reader_2na: *mut TVDB2naICReader,
    pub reader_4na: *mut VdbBlast4naReader,
    pub range_list: Cell<*mut TVDBPartialFetchingRanges>,
    pub is_initialized: bool,
}

pub type TVDBData = SVBDSRCSeqSrcData;

#[derive(Debug, Clone)]
pub struct SVDBSRCNewArgs {
    pub vdb_run_accessions: *const *const c_char,
    pub num_runs: u32,
    pub is_protein: bool,
    pub is_run_excluded: *mut bool,
    pub status: u32,
    pub is_csra: bool,
    pub include_filtered_reads: bool,
}

pub type TVDBNewArgs = SVDBSRCNewArgs;

pub fn s_vdbsrc_get_num_seqs(vdb_data_handle: &dyn Any, _dummy: Option<&mut dyn Any>) -> i32 {
    let Some(vdb_data) = vdb_data_handle.downcast_ref::<TVDBData>() else {
        return 0;
    };

    let retval = if vdb_data.is_initialized {
        vdb_data.num_seqs
    } else {
        0
    };

    if retval > i32::MAX as u64 {
        return crate::algo::blast::vdb::common_priv::VDBSRC_OVERFLOW_RV;
    }

    retval as i32
}

pub fn s_vdbsrc_get_tot_len_stats(_vdb_data_handle: &dyn Any, _dummy: Option<&mut dyn Any>) -> i64 {
    0
}

pub fn s_vdbsrc_get_num_seqs_stats(
    _vdb_data_handle: &dyn Any,
    _dummy: Option<&mut dyn Any>,
) -> i32 {
    0
}

pub fn s_vdbsrc_get_name(
    vdb_data_handle: &dyn Any,
    _dummy: Option<&mut dyn Any>,
) -> Option<String> {
    let Some(vdb_data) = vdb_data_handle.downcast_ref::<TVDBData>() else {
        return Some(K_EMPTY_VDB_NAME.to_string());
    };

    if !vdb_data.is_initialized || vdb_data.names.is_null() {
        return Some(K_EMPTY_VDB_NAME.to_string());
    }

    Some(
        unsafe { CStr::from_ptr(vdb_data.names) }
            .to_string_lossy()
            .into_owned(),
    )
}

pub fn s_vdbsrc_release_sequence(vdb_data_handle: &dyn Any, args: &mut BlastSeqSrcGetSeqArg) {
    if args.seq.is_null() {
        return;
    }

    unsafe {
        if (*args.seq).sequence_start_allocated {
            (*args.seq).sequence_start = None;
            (*args.seq).sequence_start_allocated = false;
        }
        if (*args.seq).sequence_allocated {
            (*args.seq).sequence = None;
            (*args.seq).sequence_allocated = false;
        }
    }

    let Some(vdb_data) = vdb_data_handle.downcast_ref::<TVDBData>() else {
        return;
    };
    let range_list = vdb_data.range_list.get();
    if range_list.is_null() {
        return;
    }

    unsafe {
        if (*range_list).oid == args.oid {
            if !(*range_list).ranges.is_null() && (*range_list).num_ranges > 0 {
                let len = (*range_list).num_ranges as usize * 2;
                drop(Box::from_raw(std::ptr::slice_from_raw_parts_mut(
                    (*range_list).ranges,
                    len,
                )));
            }
            drop(Box::from_raw(range_list));
            vdb_data.range_list.set(std::ptr::null_mut());
        }
    }
}

pub fn s_vdbsrc_reset_chunk_iterator(_vdb_data_handle: &dyn Any) {}

pub fn s_vdbsrc_get_supports_partial_fetching(
    vdb_data_handle: &dyn Any,
    _dummy: Option<&mut dyn Any>,
) -> bool {
    vdb_data_handle
        .downcast_ref::<TVDBData>()
        .map(|vdb_data| !vdb_data.ref_set.is_null())
        .unwrap_or(false)
}

pub fn s_vdbsrc_set_seq_range(vdb_data_handle: &dyn Any, args: &mut BlastSeqSrcSetRangesArg) {
    let Some(vdb_data) = vdb_data_handle.downcast_ref::<TVDBData>() else {
        return;
    };
    if args.num_ranges < 0 {
        return;
    }
    let range_len = args.num_ranges as usize * 2;
    if args.ranges.len() < range_len {
        return;
    }

    unsafe {
        let old_range_list = vdb_data.range_list.get();
        if !old_range_list.is_null() {
            let range_list = old_range_list;
            if !(*range_list).ranges.is_null() && (*range_list).num_ranges > 0 {
                let old_len = (*range_list).num_ranges as usize * 2;
                drop(Box::from_raw(std::ptr::slice_from_raw_parts_mut(
                    (*range_list).ranges,
                    old_len,
                )));
            }
            drop(Box::from_raw(range_list));
            vdb_data.range_list.set(std::ptr::null_mut());
        }

        let ranges_ptr = if range_len > 0 {
            let mut ranges = args.ranges[..range_len].to_vec().into_boxed_slice();
            let ranges_ptr = ranges.as_mut_ptr();
            let _ = Box::into_raw(ranges);
            ranges_ptr
        } else {
            std::ptr::null_mut()
        };

        vdb_data
            .range_list
            .set(Box::into_raw(Box::new(SVDBSRCPartialFetchingRanges {
                oid: args.oid,
                num_ranges: args.num_ranges,
                ranges: ranges_ptr,
            })));
    }
}
