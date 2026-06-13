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
