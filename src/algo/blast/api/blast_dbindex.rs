use std::sync::Arc;
use std::sync::{Mutex, OnceLock};

use crate::algo::blast::api::query_data::ILocalQueryData;
use crate::algo::blast::api::seqsrc_seqdb::CSeqDB;
use crate::algo::blast::core::blast_extend::InitHitList;
use crate::algo::blast::core::blast_lookup::LookupTableWrap;
use crate::algo::blast::core::blast_options::{InitialWordOptions, LookupTableOptions};
use crate::algo::blast::core::blast_util::BLAST_SequenceBlk;

#[derive(Clone, Debug)]
pub struct CBlastDbIndex {
    pub seqdb: Option<Arc<CSeqDB>>,
    pub local_query_data: Option<Arc<ILocalQueryData>>,
    pub index_name: String,
    pub loaded: bool,
    pub old_style: bool,
}

#[derive(Clone, Debug)]
pub struct CBlastSeqLocWrap;

pub trait CIndexedDb: Send {
    /// NCBI C++: `CIndexedDb::CheckOid`.
    fn check_oid(&mut self, oid: i32, last_vol_id: *mut i32) -> i32;

    /// NCBI C++: `CIndexedDb::EndSearchIndication`.
    fn end_search_indication(&mut self, last_vol_id: i32);

    /// NCBI C++: `CIndexedDb::DoPreSearch`.
    fn do_pre_search(
        &mut self,
        queries: Option<&mut BLAST_SequenceBlk>,
        lut_options: Option<&LookupTableOptions>,
        word_options: Option<&InitialWordOptions>,
    );

    /// NCBI C++: `CIndexedDb::SetQueryInfo`.
    fn set_query_info(&mut self, locs_wrap: Option<Arc<CBlastSeqLocWrap>>);

    /// NCBI C++: `CIndexedDb::GetResults`.
    fn get_results(&mut self, oid: u64, chunk: u64, init_hitlist: *mut InitHitList) -> u64;

    /// NCBI C++: `CIndexedDb::MinIndexWordSize`.
    fn min_index_word_size(&mut self) -> i32;

    /// NCBI C++: `CIndexedDb_New::SetMultipleThreads`.
    fn set_multiple_threads(&mut self, _multiple_threads: bool) {}

    /// NCBI C++: `CIndexedDb_New::SetNumThreads`.
    fn set_num_threads(&mut self, _n_threads: usize) {}
}

pub type DbIndexSetUsingThreadsFn = fn(bool);
pub type DbIndexSetNumThreadsFn = fn(usize);
pub type DbIndexSetQueryInfoFn = fn(&mut LookupTableWrap, Option<Arc<CBlastSeqLocWrap>>);
pub type DbIndexRunSearchFn =
    fn(Option<&mut BLAST_SequenceBlk>, Option<&LookupTableOptions>, Option<&InitialWordOptions>);

fn index_set_instance() -> &'static Mutex<Option<Box<dyn CIndexedDb>>> {
    static INDEX_SET_INSTANCE: OnceLock<Mutex<Option<Box<dyn CIndexedDb>>>> = OnceLock::new();
    INDEX_SET_INSTANCE.get_or_init(|| Mutex::new(None))
}

fn set_query_info_fn() -> &'static Mutex<DbIndexSetQueryInfoFn> {
    static SET_QUERY_INFO_FN: OnceLock<Mutex<DbIndexSetQueryInfoFn>> = OnceLock::new();
    SET_QUERY_INFO_FN.get_or_init(|| Mutex::new(null_set_query_info))
}

fn set_using_threads_fn() -> &'static Mutex<DbIndexSetUsingThreadsFn> {
    static SET_USING_THREADS_FN: OnceLock<Mutex<DbIndexSetUsingThreadsFn>> = OnceLock::new();
    SET_USING_THREADS_FN.get_or_init(|| Mutex::new(null_set_using_threads))
}

fn set_num_threads_fn() -> &'static Mutex<DbIndexSetNumThreadsFn> {
    static SET_NUM_THREADS_FN: OnceLock<Mutex<DbIndexSetNumThreadsFn>> = OnceLock::new();
    SET_NUM_THREADS_FN.get_or_init(|| Mutex::new(null_set_num_threads))
}

fn run_search_fn() -> &'static Mutex<DbIndexRunSearchFn> {
    static RUN_SEARCH_FN: OnceLock<Mutex<DbIndexRunSearchFn>> = OnceLock::new();
    RUN_SEARCH_FN.get_or_init(|| Mutex::new(null_run_search))
}

/// NCBI C++: `CIndexedDb::Index_Set_Instance`.
pub fn cindexed_db_set_instance(index_set_instance: Option<Box<dyn CIndexedDb>>) {
    if let Ok(mut current) = self::index_set_instance().lock() {
        *current = index_set_instance;
    }
}

/// NCBI C++: `NullSetUsingThreads`.
fn null_set_using_threads(_: bool) {}

/// NCBI C++: `NullSetNumThreads`.
fn null_set_num_threads(_: usize) {}

/// NCBI C++: `NullSetQueryInfo`.
fn null_set_query_info(_lt_wrap: &mut LookupTableWrap, _locs_wrap: Option<Arc<CBlastSeqLocWrap>>) {}

/// NCBI C++: `NullRunSearch`.
fn null_run_search(
    _: Option<&mut BLAST_SequenceBlk>,
    _: Option<&LookupTableOptions>,
    _: Option<&InitialWordOptions>,
) {
}

/// NCBI C++: `IndexedDbRunSearch`.
fn indexed_db_run_search(
    queries: Option<&mut BLAST_SequenceBlk>,
    lut_options: Option<&LookupTableOptions>,
    word_options: Option<&InitialWordOptions>,
) {
    let Ok(mut index_set_instance) = self::index_set_instance().lock() else {
        return;
    };
    let Some(idb) = index_set_instance.as_mut() else {
        return;
    };
    idb.do_pre_search(queries, lut_options, word_options);
}

/// NCBI C++: `IndexedDbSetUsingThreads`.
fn indexed_db_set_using_threads(multiple_threads: bool) {
    let Ok(mut index_set_instance) = self::index_set_instance().lock() else {
        return;
    };
    let Some(idb) = index_set_instance.as_mut() else {
        return;
    };
    idb.set_multiple_threads(multiple_threads);
}

/// NCBI C++: `IndexedDbSetNumThreads`.
fn indexed_db_set_num_threads(n_threads: usize) {
    let Ok(mut index_set_instance) = self::index_set_instance().lock() else {
        return;
    };
    let Some(idb) = index_set_instance.as_mut() else {
        return;
    };
    idb.set_num_threads(n_threads);
}

/// NCBI C++: `IndexedDbSetQueryInfo`.
fn indexed_db_set_query_info(
    lt_wrap: &mut LookupTableWrap,
    locs_wrap: Option<Arc<CBlastSeqLocWrap>>,
) {
    let Ok(mut index_set_instance) = self::index_set_instance().lock() else {
        return;
    };
    let Some(idb) = index_set_instance.as_mut() else {
        return;
    };

    lt_wrap.read_indexed_db = Some(s_mb_idb_get_results);
    lt_wrap.check_index_oid = Some(s_mb_idb_check_oid);
    lt_wrap.end_search_indication = Some(s_mb_idx_end_search_indication);
    idb.set_query_info(locs_wrap);
}

/// NCBI C++: `SetUpDbIndexCallbacks`.
pub fn set_up_db_index_callbacks() {
    if let Ok(mut callback) = self::set_using_threads_fn().lock() {
        *callback = indexed_db_set_using_threads;
    }
    if let Ok(mut callback) = self::set_num_threads_fn().lock() {
        *callback = indexed_db_set_num_threads;
    }
    if let Ok(mut callback) = self::set_query_info_fn().lock() {
        *callback = indexed_db_set_query_info;
    }
    if let Ok(mut callback) = self::run_search_fn().lock() {
        *callback = indexed_db_run_search;
    }
}

/// NCBI C++: `ClearDbIndexCallbacks`.
pub fn clear_db_index_callbacks() {
    if let Ok(mut callback) = self::set_using_threads_fn().lock() {
        *callback = null_set_using_threads;
    }
    if let Ok(mut callback) = self::set_num_threads_fn().lock() {
        *callback = null_set_num_threads;
    }
    if let Ok(mut callback) = self::set_query_info_fn().lock() {
        *callback = null_set_query_info;
    }
    if let Ok(mut callback) = self::run_search_fn().lock() {
        *callback = null_run_search;
    }
}

/// NCBI C++: `GetDbIndexSetUsingThreadsFn`.
pub fn get_db_index_set_using_threads_fn() -> DbIndexSetUsingThreadsFn {
    self::set_using_threads_fn()
        .lock()
        .map(|callback| *callback)
        .unwrap_or(null_set_using_threads)
}

/// NCBI C++: `GetDbIndexSetNumThreadsFn`.
pub fn get_db_index_set_num_threads_fn() -> DbIndexSetNumThreadsFn {
    self::set_num_threads_fn()
        .lock()
        .map(|callback| *callback)
        .unwrap_or(null_set_num_threads)
}

/// NCBI C++: `GetDbIndexSetQueryInfoFn`.
pub fn get_db_index_set_query_info_fn() -> DbIndexSetQueryInfoFn {
    self::set_query_info_fn()
        .lock()
        .map(|callback| *callback)
        .unwrap_or(null_set_query_info)
}

/// NCBI C++: `GetDbIndexRunSearchFn`.
pub fn get_db_index_run_search_fn() -> DbIndexRunSearchFn {
    self::run_search_fn()
        .lock()
        .map(|callback| *callback)
        .unwrap_or(null_run_search)
}

/// NCBI C++: `s_MB_IdxEndSearchIndication`.
unsafe extern "C" fn s_mb_idx_end_search_indication(last_vol_id: i32) {
    let Ok(mut index_set_instance) = self::index_set_instance().lock() else {
        return;
    };
    let Some(idb) = index_set_instance.as_mut() else {
        return;
    };
    idb.end_search_indication(last_vol_id);
}

/// NCBI C++: `s_MB_IdbCheckOid`.
unsafe extern "C" fn s_mb_idb_check_oid(oid: i32, last_vol_id: *mut i32) -> i32 {
    debug_assert!(oid >= 0);
    let Ok(mut index_set_instance) = self::index_set_instance().lock() else {
        return 2;
    };
    let Some(idb) = index_set_instance.as_mut() else {
        return 2;
    };
    idb.check_oid(oid, last_vol_id)
}

/// NCBI C++: `s_MB_IdbGetResults`.
unsafe extern "C" fn s_mb_idb_get_results(
    oid_i: i32,
    chunk_i: i32,
    init_hitlist: *mut InitHitList,
) -> u64 {
    debug_assert!(oid_i >= 0);
    debug_assert!(chunk_i >= 0);
    debug_assert!(!init_hitlist.is_null());

    let Ok(mut index_set_instance) = self::index_set_instance().lock() else {
        return 0;
    };
    let Some(idb) = index_set_instance.as_mut() else {
        return 0;
    };

    idb.get_results(oid_i as u64, chunk_i as u64, init_hitlist)
}
