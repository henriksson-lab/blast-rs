//! Safe Rust wrappers around the BLAST C core engine.

pub use blast_core_sys as ffi;

pub mod blastn;
pub mod diagnostics;
pub mod encoding;
pub mod engine;
pub mod extend;
pub mod filter;
pub mod gapalign;
pub mod gapinfo;
pub mod greedy;
pub mod hits;
pub mod hspstream;
pub mod itree;
pub mod listnode;
pub mod lookup;
pub mod math;
pub mod options;
pub mod parameters;
pub mod program;
pub mod protein;
pub mod pssm;
pub mod rps;
pub mod queryinfo;
pub mod search;
pub mod sequence;
pub mod seqsrc;
pub mod stat;
pub mod traceback;
pub mod matrix;
pub mod util;

use std::ptr;

/// Result type for BLAST operations.
pub type Result<T> = std::result::Result<T, BlastError>;

/// BLAST error type.
#[derive(Debug)]
pub struct BlastError {
    pub code: i32,
    pub message: String,
}

impl std::fmt::Display for BlastError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "BLAST error {}: {}", self.code, self.message)
    }
}

impl std::error::Error for BlastError {}

// ---------------------------------------------------------------------------
// RAII wrappers for C core structures
// ---------------------------------------------------------------------------

/// Owned wrapper around `BLAST_SequenceBlk`.
pub struct SequenceBlk {
    pub ptr: *mut ffi::BLAST_SequenceBlk,
}

impl Drop for SequenceBlk {
    fn drop(&mut self) {
        if !self.ptr.is_null() {
            unsafe {
                ffi::BlastSequenceBlkFree(self.ptr);
            }
        }
    }
}

/// Owned wrapper around `BlastQueryInfo`.
pub struct QueryInfo {
    pub ptr: *mut ffi::BlastQueryInfo,
}

impl Drop for QueryInfo {
    fn drop(&mut self) {
        if !self.ptr.is_null() {
            unsafe {
                ffi::BlastQueryInfoFree(self.ptr);
            }
        }
    }
}

impl QueryInfo {
    /// Create a new BlastQueryInfo for the given program and number of queries.
    pub fn new(program: ffi::EBlastProgramType, num_queries: i32) -> Result<Self> {
        let ptr = unsafe { ffi::BlastQueryInfoNew(program, num_queries) };
        if ptr.is_null() {
            return Err(BlastError {
                code: -1,
                message: "Failed to allocate BlastQueryInfo".into(),
            });
        }
        Ok(QueryInfo { ptr })
    }
}

/// Owned wrapper around `BlastScoreBlk`.
pub struct ScoreBlk {
    pub ptr: *mut ffi::BlastScoreBlk,
}

impl Drop for ScoreBlk {
    fn drop(&mut self) {
        if !self.ptr.is_null() {
            unsafe {
                ffi::BlastScoreBlkFree(self.ptr);
            }
        }
    }
}

impl ScoreBlk {
    pub fn new(alphabet: u8, num_contexts: i32) -> Result<Self> {
        let ptr = unsafe { ffi::BlastScoreBlkNew(alphabet, num_contexts) };
        if ptr.is_null() {
            return Err(BlastError {
                code: -1,
                message: "Failed to allocate BlastScoreBlk".into(),
            });
        }
        Ok(ScoreBlk { ptr })
    }
}

/// Owned wrapper around `BlastHSPResults`.
pub struct HspResults {
    pub ptr: *mut ffi::BlastHSPResults,
}

impl Drop for HspResults {
    fn drop(&mut self) {
        if !self.ptr.is_null() {
            unsafe {
                ffi::Blast_HSPResultsFree(self.ptr);
            }
        }
    }
}

impl HspResults {
    pub fn new(num_queries: i32) -> Result<Self> {
        let ptr = unsafe { ffi::Blast_HSPResultsNew(num_queries) };
        if ptr.is_null() {
            return Err(BlastError {
                code: -1,
                message: "Failed to allocate BlastHSPResults".into(),
            });
        }
        Ok(HspResults { ptr })
    }

    /// Number of queries in results.
    pub fn num_queries(&self) -> i32 {
        unsafe { (*self.ptr).num_queries }
    }
}

/// Owned wrapper around `BlastHSPStream`.
pub struct HspStream {
    pub ptr: *mut ffi::BlastHSPStream,
}

impl Drop for HspStream {
    fn drop(&mut self) {
        if !self.ptr.is_null() {
            unsafe {
                ffi::BlastHSPStreamFree(self.ptr);
            }
        }
    }
}

/// Owned wrapper around `BlastDiagnostics`.
pub struct Diagnostics {
    pub ptr: *mut ffi::BlastDiagnostics,
}

impl Drop for Diagnostics {
    fn drop(&mut self) {
        if !self.ptr.is_null() {
            unsafe {
                ffi::Blast_DiagnosticsFree(self.ptr);
            }
        }
    }
}

impl Diagnostics {
    pub fn new() -> Result<Self> {
        let ptr = unsafe { ffi::Blast_DiagnosticsInit() };
        if ptr.is_null() {
            return Err(BlastError {
                code: -1,
                message: "Failed to allocate BlastDiagnostics".into(),
            });
        }
        Ok(Diagnostics { ptr })
    }
}

/// Owned wrapper around `LookupTableWrap`.
pub struct LookupTable {
    pub ptr: *mut ffi::LookupTableWrap,
}

impl Drop for LookupTable {
    fn drop(&mut self) {
        if !self.ptr.is_null() {
            unsafe {
                ffi::LookupTableWrapFree(self.ptr);
            }
        }
    }
}

/// Owned wrapper around `BlastSeqSrc`.
pub struct SeqSrc {
    pub ptr: *mut ffi::BlastSeqSrc,
}

impl Drop for SeqSrc {
    fn drop(&mut self) {
        if !self.ptr.is_null() {
            unsafe {
                ffi::BlastSeqSrcFree(self.ptr);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Options wrappers
// ---------------------------------------------------------------------------

/// Helper macro for simple options with (program, &mut ptr) constructor.
macro_rules! options_wrapper_simple {
    ($name:ident, $ffi_type:ty, $free_fn:path, $new_fn:path) => {
        pub struct $name {
            pub ptr: *mut $ffi_type,
        }

        impl Drop for $name {
            fn drop(&mut self) {
                if !self.ptr.is_null() {
                    unsafe {
                        $free_fn(self.ptr);
                    }
                }
            }
        }

        impl $name {
            pub fn new(program: ffi::EBlastProgramType) -> Result<Self> {
                let mut ptr: *mut $ffi_type = ptr::null_mut();
                let status = unsafe { $new_fn(program, &mut ptr) };
                if status != 0 || ptr.is_null() {
                    return Err(BlastError {
                        code: status as i32,
                        message: format!("Failed to create {}", stringify!($name)),
                    });
                }
                Ok($name { ptr })
            }
        }
    };
}

options_wrapper_simple!(
    ScoringOptions,
    ffi::BlastScoringOptions,
    ffi::BlastScoringOptionsFree,
    ffi::BlastScoringOptionsNew
);

options_wrapper_simple!(
    InitialWordOptions,
    ffi::BlastInitialWordOptions,
    ffi::BlastInitialWordOptionsFree,
    ffi::BlastInitialWordOptionsNew
);

// ExtensionOptions takes an extra `gapped` parameter
pub struct ExtensionOptions {
    pub ptr: *mut ffi::BlastExtensionOptions,
}

impl Drop for ExtensionOptions {
    fn drop(&mut self) {
        if !self.ptr.is_null() {
            unsafe {
                ffi::BlastExtensionOptionsFree(self.ptr);
            }
        }
    }
}

impl ExtensionOptions {
    pub fn new(program: ffi::EBlastProgramType, gapped: bool) -> Result<Self> {
        let mut ptr: *mut ffi::BlastExtensionOptions = ptr::null_mut();
        let status =
            unsafe { ffi::BlastExtensionOptionsNew(program, &mut ptr, gapped as u8) };
        if status != 0 || ptr.is_null() {
            return Err(BlastError {
                code: status as i32,
                message: "Failed to create ExtensionOptions".into(),
            });
        }
        Ok(ExtensionOptions { ptr })
    }
}

// HitSavingOptions takes an extra `gapped_calculation` parameter
pub struct HitSavingOptions {
    pub ptr: *mut ffi::BlastHitSavingOptions,
}

impl Drop for HitSavingOptions {
    fn drop(&mut self) {
        if !self.ptr.is_null() {
            unsafe {
                ffi::BlastHitSavingOptionsFree(self.ptr);
            }
        }
    }
}

impl HitSavingOptions {
    pub fn new(program: ffi::EBlastProgramType, gapped: bool) -> Result<Self> {
        let mut ptr: *mut ffi::BlastHitSavingOptions = ptr::null_mut();
        let status =
            unsafe { ffi::BlastHitSavingOptionsNew(program, &mut ptr, gapped as u8) };
        if status != 0 || ptr.is_null() {
            return Err(BlastError {
                code: status as i32,
                message: "Failed to create HitSavingOptions".into(),
            });
        }
        Ok(HitSavingOptions { ptr })
    }
}

/// EffectiveLengthsOptions has a different constructor signature (no program arg).
pub struct EffectiveLengthsOptions {
    pub ptr: *mut ffi::BlastEffectiveLengthsOptions,
}

impl Drop for EffectiveLengthsOptions {
    fn drop(&mut self) {
        if !self.ptr.is_null() {
            unsafe {
                ffi::BlastEffectiveLengthsOptionsFree(self.ptr);
            }
        }
    }
}

impl EffectiveLengthsOptions {
    pub fn new() -> Result<Self> {
        let mut ptr: *mut ffi::BlastEffectiveLengthsOptions = ptr::null_mut();
        let status = unsafe { ffi::BlastEffectiveLengthsOptionsNew(&mut ptr) };
        if status != 0 || ptr.is_null() {
            return Err(BlastError {
                code: status as i32,
                message: "Failed to create EffectiveLengthsOptions".into(),
            });
        }
        Ok(EffectiveLengthsOptions { ptr })
    }
}

/// Database options.
pub struct DatabaseOptions {
    pub ptr: *mut ffi::BlastDatabaseOptions,
}

impl Drop for DatabaseOptions {
    fn drop(&mut self) {
        if !self.ptr.is_null() {
            unsafe {
                ffi::BlastDatabaseOptionsFree(self.ptr);
            }
        }
    }
}

impl DatabaseOptions {
    pub fn new() -> Result<Self> {
        let mut ptr: *mut ffi::BlastDatabaseOptions = ptr::null_mut();
        let status = unsafe { ffi::BlastDatabaseOptionsNew(&mut ptr) };
        if status != 0 || ptr.is_null() {
            return Err(BlastError {
                code: status as i32,
                message: "Failed to create DatabaseOptions".into(),
            });
        }
        Ok(DatabaseOptions { ptr })
    }
}

// ---------------------------------------------------------------------------
// High-level search function
// ---------------------------------------------------------------------------

/// Run a full BLAST search using the C core engine.
///
/// This wraps `Blast_RunFullSearch`, which performs preliminary search + traceback.
pub fn run_full_search(
    program: ffi::EBlastProgramType,
    query: &mut SequenceBlk,
    query_info: &mut QueryInfo,
    seq_src: &SeqSrc,
    score_blk: &mut ScoreBlk,
    scoring_opts: &ScoringOptions,
    lookup: &mut LookupTable,
    word_opts: &InitialWordOptions,
    ext_opts: &ExtensionOptions,
    hit_opts: &HitSavingOptions,
    eff_len_opts: &EffectiveLengthsOptions,
    db_opts: &DatabaseOptions,
    hsp_stream: &mut HspStream,
    diagnostics: &mut Diagnostics,
) -> Result<HspResults> {
    let mut results_ptr: *mut ffi::BlastHSPResults = ptr::null_mut();

    let status = unsafe {
        ffi::Blast_RunFullSearch(
            program,
            query.ptr,
            query_info.ptr,
            seq_src.ptr,
            score_blk.ptr,
            scoring_opts.ptr,
            lookup.ptr,
            word_opts.ptr,
            ext_opts.ptr,
            hit_opts.ptr,
            eff_len_opts.ptr,
            ptr::null(), // psi_options
            db_opts.ptr,
            hsp_stream.ptr,
            ptr::null(), // rps_info
            diagnostics.ptr,
            &mut results_ptr,
            None,        // interrupt callback
            ptr::null_mut(), // progress info
        )
    };

    if status != 0 {
        return Err(BlastError {
            code: status,
            message: "Blast_RunFullSearch failed".into(),
        });
    }

    if results_ptr.is_null() {
        return Err(BlastError {
            code: -1,
            message: "Blast_RunFullSearch returned null results".into(),
        });
    }

    Ok(HspResults { ptr: results_ptr })
}
