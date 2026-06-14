use std::sync::Arc;

use super::blast_options_cxx::CBlastOptionsMemento;
use super::setup_factory::SInternalData;

/// NCBI C++: `CPrelimSearchRunner` (`prelim_search_runner.hpp`).
#[derive(Clone, Debug)]
pub struct CPrelimSearchRunner {
    pub m_InternalData: Option<Arc<SInternalData>>,
    pub m_OptsMemento: Option<Arc<CBlastOptionsMemento>>,
}

/// NCBI C++: `CPrelimSearchThread` (`prelim_search_runner.hpp`).
#[derive(Debug)]
pub struct CPrelimSearchThread {
    pub m_InternalData: SInternalData,
    pub m_OptsMemento: Option<Arc<CBlastOptionsMemento>>,
}
