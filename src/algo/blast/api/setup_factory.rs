use std::fmt;

use crate::algo::blast::api::rps_aux::CBlastRPSInfo;
use crate::algo::blast::core::blast_hspstream::HspStream;
use crate::algo::blast::core::blast_lookup::LookupTableWrap;
use crate::algo::blast::core::blast_query_info::BlastQueryInfo;
use crate::algo::blast::core::blast_seqsrc::BlastSeqSrc;
use crate::algo::blast::core::blast_stat::BlastScoreBlk;
use crate::algo::blast::core::blast_util::{BLAST_SequenceBlk, SBlastProgress};
use crate::algo::blast::core::pattern::PhiPatternSearchBlk;

/// NCBI C++: `CStructWrapper<TData>` (`setup_factory.hpp`).
pub struct CStructWrapper<TData> {
    pub m_Data: Option<TData>,
    pub m_DeleteFunction: Option<fn(TData) -> Option<TData>>,
}

impl<TData> CStructWrapper<TData> {
    /// NCBI C++: `CStructWrapper::CStructWrapper`.
    pub fn new(obj: Option<TData>, dfun: Option<fn(TData) -> Option<TData>>) -> Self {
        Self {
            m_Data: obj,
            m_DeleteFunction: dfun,
        }
    }

    /// NCBI C++: `CStructWrapper::GetPointer`.
    pub fn get_pointer(&mut self) -> Option<&mut TData> {
        self.m_Data.as_mut()
    }

    /// NCBI C++: `CStructWrapper::operator*`.
    pub fn deref_mut(&mut self) -> Option<&mut TData> {
        self.m_Data.as_mut()
    }

    /// NCBI C++: `CStructWrapper::operator->`.
    pub fn arrow(&mut self) -> Option<&mut TData> {
        self.m_Data.as_mut()
    }
}

impl<TData> Drop for CStructWrapper<TData> {
    fn drop(&mut self) {
        if let Some(data) = self.m_Data.take() {
            if let Some(delete_function) = self.m_DeleteFunction {
                let _ = delete_function(data);
            }
        }
    }
}

/// NCBI C++: `WrapStruct` (`setup_factory.hpp`).
pub fn wrap_struct<TData>(
    obj: Option<TData>,
    del: Option<fn(TData) -> Option<TData>>,
) -> CStructWrapper<TData> {
    CStructWrapper::new(obj, del)
}

pub type TBlastScoreBlk = CStructWrapper<BlastScoreBlk>;
pub type TLookupTableWrap = CStructWrapper<LookupTableWrap>;
pub type TBlastDiagnostics = CStructWrapper<BlastDiagnostics>;
pub type TBlastHSPStream = CStructWrapper<HspStream>;
pub type TBlastSeqSrc = CStructWrapper<BlastSeqSrc>;
pub type TSPHIPatternSearchBlk = CStructWrapper<PhiPatternSearchBlk>;
pub type TInterruptFnPtr = usize;

/// NCBI C++: `SInternalData` (`setup_factory.hpp`).
pub struct SInternalData {
    pub m_Queries: Option<BLAST_SequenceBlk>,
    pub m_QueryInfo: Option<BlastQueryInfo>,
    pub m_ScoreBlk: Option<TBlastScoreBlk>,
    pub m_LookupTable: Option<TLookupTableWrap>,
    pub m_Diagnostics: Option<TBlastDiagnostics>,
    pub m_HspStream: Option<TBlastHSPStream>,
    pub m_SeqSrc: Option<TBlastSeqSrc>,
    pub m_RpsData: Option<CBlastRPSInfo>,
    pub m_FnInterrupt: Option<TInterruptFnPtr>,
    pub m_ProgressMonitor: Option<SBlastProgress>,
}

impl SInternalData {
    /// NCBI C++: `SInternalData::SInternalData`.
    pub fn new() -> Self {
        Self {
            m_Queries: None,
            m_QueryInfo: None,
            m_ScoreBlk: None,
            m_LookupTable: None,
            m_Diagnostics: None,
            m_HspStream: None,
            m_SeqSrc: None,
            m_RpsData: None,
            m_FnInterrupt: None,
            m_ProgressMonitor: None,
        }
    }
}

impl fmt::Debug for SInternalData {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("SInternalData")
            .field("m_Queries", &self.m_Queries.is_some())
            .field("m_QueryInfo", &self.m_QueryInfo.is_some())
            .field("m_ScoreBlk", &self.m_ScoreBlk.is_some())
            .field("m_LookupTable", &self.m_LookupTable.is_some())
            .field("m_Diagnostics", &self.m_Diagnostics.is_some())
            .field("m_HspStream", &self.m_HspStream.is_some())
            .field("m_SeqSrc", &self.m_SeqSrc.is_some())
            .field("m_RpsData", &self.m_RpsData.is_some())
            .field("m_FnInterrupt", &self.m_FnInterrupt.is_some())
            .field("m_ProgressMonitor", &self.m_ProgressMonitor.is_some())
            .finish()
    }
}

#[derive(Clone, Debug)]
pub struct BlastDiagnostics;

/// NCBI C++: `CThreadable` (`setup_factory.hpp`).
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CThreadable {
    pub m_NumThreads: usize,
}

impl CThreadable {
    pub const K_MIN_NUM_THREADS: usize = 1;

    /// NCBI C++: `CThreadable::CThreadable`.
    pub fn new() -> Self {
        Self {
            m_NumThreads: Self::K_MIN_NUM_THREADS,
        }
    }

    /// NCBI C++: `CThreadable::SetNumberOfThreads`.
    pub fn set_number_of_threads(&mut self, nthreads: usize) {
        self.m_NumThreads = if nthreads == 0 {
            Self::K_MIN_NUM_THREADS
        } else {
            nthreads
        };
    }

    /// NCBI C++: `CThreadable::GetNumberOfThreads`.
    pub fn get_number_of_threads(&self) -> usize {
        debug_assert!(self.m_NumThreads >= Self::K_MIN_NUM_THREADS);
        self.m_NumThreads
    }

    /// NCBI C++: `CThreadable::IsMultiThreaded`.
    pub fn is_multi_threaded(&self) -> bool {
        self.m_NumThreads > Self::K_MIN_NUM_THREADS
    }
}

/// NCBI C++: `SDatabaseScanData` (`setup_factory.hpp/.cpp`).
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SDatabaseScanData {
    pub kNoPhiBlastPattern: i32,
    pub m_NumPatOccurInDB: i32,
}

impl SDatabaseScanData {
    /// NCBI C++: `SDatabaseScanData::SDatabaseScanData`.
    pub fn new() -> Self {
        Self {
            kNoPhiBlastPattern: -1,
            m_NumPatOccurInDB: -1,
        }
    }
}
