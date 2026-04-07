//! Implement the C BlastSeqSrc vtable backed by our Rust BlastDb reader.

use crate::index::BlastDb;
use blast_core_sys as ffi;
use std::ffi::CString;
use std::os::raw::{c_char, c_int, c_void};
use std::ptr;
use std::sync::{Arc, Mutex};

/// Internal data stored in the BlastSeqSrc's DataStructure pointer.
struct RustSeqSrcData {
    db: BlastDb,
    name: CString,
    /// Mutex-protected counter for thread-safe iteration.
    next_oid: Mutex<i32>,
}

/// Create a `BlastSeqSrc` backed by a Rust `BlastDb`.
///
/// The returned pointer is owned by the caller and must be freed with
/// `BlastSeqSrcFree`. The `BlastDb` is moved into the SeqSrc.
pub fn create_seq_src(db: BlastDb) -> *mut ffi::BlastSeqSrc {
    let name = CString::new(db.title.as_str()).unwrap_or_else(|_| CString::new("unknown").unwrap());
    let data = Box::new(RustSeqSrcData {
        db,
        name,
        next_oid: Mutex::new(0),
    });
    let data_ptr = Box::into_raw(data) as *mut c_void;

    let new_info = ffi::BlastSeqSrcNewInfo {
        constructor: Some(rust_seqsrc_constructor),
        ctor_argument: data_ptr,
    };

    unsafe { ffi::BlastSeqSrcNew(&new_info) }
}

// ---------------------------------------------------------------------------
// BlastSeqSrc vtable callbacks
// ---------------------------------------------------------------------------

unsafe extern "C" fn rust_seqsrc_constructor(
    seq_src: *mut ffi::BlastSeqSrc,
    arg: *mut c_void,
) -> *mut ffi::BlastSeqSrc {
    // `arg` is our RustSeqSrcData pointer
    unsafe {
        ffi::_BlastSeqSrcImpl_SetDataStructure(seq_src, arg);
        ffi::_BlastSeqSrcImpl_SetDeleteFnPtr(seq_src, Some(rust_seqsrc_destructor));
        ffi::_BlastSeqSrcImpl_SetCopyFnPtr(seq_src, None);
        ffi::_BlastSeqSrcImpl_SetGetNumSeqs(seq_src, Some(rust_get_num_seqs));
        ffi::_BlastSeqSrcImpl_SetGetNumSeqsStats(seq_src, Some(rust_get_num_seqs));
        ffi::_BlastSeqSrcImpl_SetGetMaxSeqLen(seq_src, Some(rust_get_max_seq_len));
        ffi::_BlastSeqSrcImpl_SetGetMinSeqLen(seq_src, Some(rust_get_min_seq_len));
        ffi::_BlastSeqSrcImpl_SetGetAvgSeqLen(seq_src, Some(rust_get_avg_seq_len));
        ffi::_BlastSeqSrcImpl_SetGetTotLen(seq_src, Some(rust_get_tot_len));
        ffi::_BlastSeqSrcImpl_SetGetTotLenStats(seq_src, Some(rust_get_tot_len));
        ffi::_BlastSeqSrcImpl_SetGetName(seq_src, Some(rust_get_name));
        ffi::_BlastSeqSrcImpl_SetGetIsProt(seq_src, Some(rust_get_is_prot));
        ffi::_BlastSeqSrcImpl_SetGetSequence(seq_src, Some(rust_get_sequence));
        ffi::_BlastSeqSrcImpl_SetGetSeqLen(seq_src, Some(rust_get_seq_len));
        ffi::_BlastSeqSrcImpl_SetReleaseSequence(seq_src, Some(rust_release_sequence));
        ffi::_BlastSeqSrcImpl_SetIterNext(seq_src, Some(rust_iter_next));
        ffi::_BlastSeqSrcImpl_SetResetChunkIterator(seq_src, Some(rust_reset_chunk_iter));
        ffi::_BlastSeqSrcImpl_SetSetNumberOfThreads(seq_src, Some(rust_set_num_threads));
        ffi::_BlastSeqSrcImpl_SetGetSupportsPartialFetching(seq_src, None);
        ffi::_BlastSeqSrcImpl_SetSetSeqRange(seq_src, None);
    }
    seq_src
}

unsafe extern "C" fn rust_seqsrc_destructor(
    seq_src: *mut ffi::BlastSeqSrc,
) -> *mut ffi::BlastSeqSrc {
    let data_ptr = unsafe { ffi::_BlastSeqSrcImpl_GetDataStructure(seq_src) };
    if !data_ptr.is_null() {
        unsafe {
            drop(Box::from_raw(data_ptr as *mut RustSeqSrcData));
            ffi::_BlastSeqSrcImpl_SetDataStructure(seq_src, ptr::null_mut());
        }
    }
    // Must return NULL — BlastSeqSrcFree asserts this
    ptr::null_mut()
}

#[inline]
unsafe fn get_data<'a>(impl_ptr: *mut c_void) -> &'a RustSeqSrcData {
    unsafe { &*(impl_ptr as *const RustSeqSrcData) }
}

unsafe extern "C" fn rust_get_num_seqs(impl_ptr: *mut c_void, _arg: *mut c_void) -> ffi::Int4 {
    let data = unsafe { get_data(impl_ptr) };
    data.db.num_oids as ffi::Int4
}

unsafe extern "C" fn rust_get_max_seq_len(impl_ptr: *mut c_void, _arg: *mut c_void) -> ffi::Int4 {
    let data = unsafe { get_data(impl_ptr) };
    data.db.max_seq_len as ffi::Int4
}

unsafe extern "C" fn rust_get_min_seq_len(impl_ptr: *mut c_void, _arg: *mut c_void) -> ffi::Int4 {
    // Not tracked by BlastDb currently; return 1 as minimum
    1
}

unsafe extern "C" fn rust_get_avg_seq_len(impl_ptr: *mut c_void, _arg: *mut c_void) -> ffi::Int4 {
    let data = unsafe { get_data(impl_ptr) };
    if data.db.num_oids > 0 {
        (data.db.total_length / data.db.num_oids as u64) as ffi::Int4
    } else {
        0
    }
}

unsafe extern "C" fn rust_get_tot_len(impl_ptr: *mut c_void, _arg: *mut c_void) -> ffi::Int8 {
    let data = unsafe { get_data(impl_ptr) };
    data.db.total_length as ffi::Int8
}

unsafe extern "C" fn rust_get_name(
    impl_ptr: *mut c_void,
    _arg: *mut c_void,
) -> *const c_char {
    let data = unsafe { get_data(impl_ptr) };
    data.name.as_ptr()
}

unsafe extern "C" fn rust_get_is_prot(impl_ptr: *mut c_void, _arg: *mut c_void) -> ffi::Boolean {
    let data = unsafe { get_data(impl_ptr) };
    match data.db.db_type {
        crate::index::DbType::Protein => 1,
        crate::index::DbType::Nucleotide => 0,
    }
}

unsafe extern "C" fn rust_get_seq_len(impl_ptr: *mut c_void, arg: *mut c_void) -> ffi::Int4 {
    let data = unsafe { get_data(impl_ptr) };
    if arg.is_null() {
        return -1;
    }
    let oid = unsafe { *(arg as *const i32) };
    if oid < 0 || oid as u32 >= data.db.num_oids {
        return -1;
    }
    data.db.get_seq_len(oid as u32) as ffi::Int4
}

/// NCBI4NA to BLASTNA conversion table.
/// BLASTNA is a permutation where A=0,C=1,G=2,T=3 (matching NCBI2na).
const NCBI4NA_TO_BLASTNA: [u8; 16] = [
    15, // 0 = gap -> sentinel
    0,  // 1 = A
    1,  // 2 = C
    6,  // 3 = M (AC)
    2,  // 4 = G
    4,  // 5 = R (AG)
    9,  // 6 = S (CG)
    13, // 7 = V (ACG)
    3,  // 8 = T
    8,  // 9 = W (AT)
    5,  // 10 = Y (CT)
    12, // 11 = H (ACT)
    7,  // 12 = K (GT)
    11, // 13 = D (AGT)
    10, // 14 = B (CGT)
    14, // 15 = N (ACGT)
];

/// Decode a packed 2-bit nucleotide sequence to BLASTNA encoding.
/// The packed format stores 4 bases per byte: bits 7-6 = base0, 5-4 = base1, etc.
/// NCBI2na: A=0, C=1, G=2, T=3 — these map directly to BLASTNA 0-3.
fn decode_ncbi2na_to_blastna(packed: &[u8], seq_len: u32) -> Vec<u8> {
    let mut result = Vec::with_capacity(seq_len as usize + 2);
    // Leading sentinel byte
    result.push(15); // BLASTNA sentinel

    let full_bytes = seq_len as usize / 4;
    let remainder = (seq_len % 4) as usize;

    for i in 0..full_bytes {
        let byte = packed[i];
        result.push((byte >> 6) & 3); // base 0
        result.push((byte >> 4) & 3); // base 1
        result.push((byte >> 2) & 3); // base 2
        result.push(byte & 3);        // base 3
    }
    if remainder > 0 && full_bytes < packed.len() {
        let byte = packed[full_bytes];
        for j in 0..remainder {
            result.push((byte >> (6 - 2 * j)) & 3);
        }
    }

    // Trailing sentinel byte
    result.push(15);
    result
}

unsafe extern "C" fn rust_get_sequence(
    impl_ptr: *mut c_void,
    arg: *mut ffi::BlastSeqSrcGetSeqArg,
) -> ffi::Int2 {
    if arg.is_null() {
        return -1;
    }
    let data = unsafe { get_data(impl_ptr) };
    let get_arg = unsafe { &mut *arg };
    let oid = get_arg.oid as u32;

    if oid >= data.db.num_oids {
        return -1;
    }

    let seq_data = data.db.get_sequence(oid);
    let seq_len = data.db.get_seq_len(oid);
    let encoding = get_arg.encoding;

    let is_nuc = data.db.db_type == crate::index::DbType::Nucleotide;
    let needs_decoded = encoding == ffi::EBlastEncoding_eBlastEncodingNucleotide
        || encoding == ffi::EBlastEncoding_eBlastEncodingNcbi4na;

    if is_nuc && !needs_decoded {
        // Preliminary search: provide raw NCBI2na packed data, NO sentinel bytes
        // This matches how seqsrc_seqdb.cpp works with eBlastEncodingProtein
        let buf_len = seq_data.len();
        let buf = unsafe { libc::malloc(buf_len) as *mut u8 };
        if buf.is_null() {
            return -1;
        }
        unsafe {
            ptr::copy_nonoverlapping(seq_data.as_ptr(), buf, buf_len);
        }

        // Allocate sequence block if needed
        if get_arg.seq.is_null() {
            let mut seq_blk: *mut ffi::BLAST_SequenceBlk = ptr::null_mut();
            let rc = unsafe { ffi::BlastSeqBlkNew(&mut seq_blk) };
            if rc != 0 || seq_blk.is_null() {
                unsafe { libc::free(buf as *mut c_void) };
                return -1;
            }
            get_arg.seq = seq_blk;
        }

        let seq_blk = unsafe { &mut *get_arg.seq };
        seq_blk.sequence = buf;
        seq_blk.sequence_start = ptr::null_mut();
        seq_blk.length = seq_len as i32;
        seq_blk.sequence_allocated = 1;
        seq_blk.sequence_start_allocated = 0;
        seq_blk.oid = oid as i32;
    } else if is_nuc {
        // Traceback: provide decoded BLASTNA with sentinel bytes
        let decoded = decode_ncbi2na_to_blastna(seq_data, seq_len);
        let buf_len = decoded.len();
        let buf = unsafe { libc::malloc(buf_len) as *mut u8 };
        if buf.is_null() {
            return -1;
        }
        unsafe {
            ptr::copy_nonoverlapping(decoded.as_ptr(), buf, buf_len);
        }

        let rc = unsafe {
            ffi::BlastSetUp_SeqBlkNew(buf, seq_len as i32, &mut get_arg.seq, 1)
        };
        if rc != 0 {
            unsafe { libc::free(buf as *mut c_void) };
            return -1;
        }
        let seq_blk = unsafe { &mut *get_arg.seq };
        seq_blk.oid = oid as i32;
    } else {
        // Protein
        let rc = unsafe {
            ffi::BlastSetUp_SeqBlkNew(
                seq_data.as_ptr() as *const u8,
                seq_len as i32,
                &mut get_arg.seq,
                0, // not allocated — just point to mmap'd data
            )
        };
        if rc != 0 {
            return -1;
        }
        let seq_blk = unsafe { &mut *get_arg.seq };
        seq_blk.oid = oid as i32;
    }

    0
}

unsafe extern "C" fn rust_release_sequence(
    _impl_ptr: *mut c_void,
    arg: *mut ffi::BlastSeqSrcGetSeqArg,
) {
    if arg.is_null() {
        return;
    }
    let get_arg = unsafe { &mut *arg };
    if !get_arg.seq.is_null() {
        unsafe {
            ffi::BlastSequenceBlkClean(get_arg.seq);
        }
    }
}

unsafe extern "C" fn rust_iter_next(
    impl_ptr: *mut c_void,
    itr: *mut ffi::BlastSeqSrcIterator,
) -> ffi::Int4 {
    let data = unsafe { get_data(impl_ptr) };
    let itr = unsafe { &mut *itr };

    // Simple sequential iteration using oid_range for chunk-based access
    if itr.current_pos >= itr.oid_range[1] as u32 {
        // Need a new chunk
        let mut next = data.next_oid.lock().unwrap();
        let start = *next;
        if start >= data.db.num_oids as i32 {
            return ffi::BLAST_SEQSRC_EOF; // -1
        }
        let chunk = if itr.chunk_sz > 0 { itr.chunk_sz as i32 } else { 1024 };
        let end = std::cmp::min(start + chunk, data.db.num_oids as i32);
        *next = end;
        drop(next);

        itr.oid_range[0] = start;
        itr.oid_range[1] = end;
        itr.current_pos = start as u32;
    }

    let oid = itr.current_pos as ffi::Int4;
    itr.current_pos += 1;
    oid
}

unsafe extern "C" fn rust_reset_chunk_iter(impl_ptr: *mut c_void) {
    let data = unsafe { get_data(impl_ptr) };
    *data.next_oid.lock().unwrap() = 0;
}

unsafe extern "C" fn rust_set_num_threads(
    _impl_ptr: *mut c_void,
    _num_threads: c_int,
) {
    // No-op for now; Rayon handles threading
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::index::BlastDb;
    use std::path::PathBuf;

    fn test_db_path() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .join("ncbi-blast-2.17.0+-src/c++/src/algo/blast/unit_tests/api/data/seqn")
    }

    #[test]
    fn test_create_seq_src() {
        let db = BlastDb::open(&test_db_path()).unwrap();
        let num_oids = db.num_oids;
        let seq_src = create_seq_src(db);
        assert!(!seq_src.is_null());

        // Test that we can query it through the C API
        let n = unsafe { ffi::BlastSeqSrcGetNumSeqs(seq_src) };
        assert_eq!(n, num_oids as i32);

        // Clean up
        unsafe {
            ffi::BlastSeqSrcFree(seq_src);
        }
    }
}
