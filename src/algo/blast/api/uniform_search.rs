use std::sync::Arc;

pub use crate::algo::blast::api::blast_options_cxx::EProgram;
use crate::algo::blast::api::blast_results::CSearchResultSet;
use crate::algo::blast::api::query_data::CSeqId;
use crate::algo::blast::api::seqsrc_seqdb::CSeqDB;

/// NCBI C++: `CSearchException::EErrCode` (`uniform_search.hpp`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum CSearchExceptionErrCode {
    EConfigErr,
    EMemErr,
    EInternal,
}

/// NCBI C++: `CSearchDatabase::EMoleculeType` (`uniform_search.hpp`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum EMoleculeType {
    EBlastDbIsProtein,
    EBlastDbIsNucleotide,
}

/// NCBI C++: `ESubjectMaskingType` users in `uniform_search.hpp`.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ESubjectMaskingType {
    ESoftSubjMasking,
    EHardSubjMasking,
}

/// NCBI C++: `CSearchDatabase` (`uniform_search.hpp`).
#[derive(Clone, Debug)]
pub struct CSearchDatabase {
    pub m_DbName: String,
    pub m_MolType: EMoleculeType,
    pub m_EntrezQueryLimitation: String,
    pub m_GiList: Option<Arc<CSeqDBGiList>>,
    pub m_NegativeGiList: Option<Arc<CSeqDBGiList>>,
    pub m_GiListSet: bool,
    pub m_FilteringAlgorithmString: String,
    pub m_FilteringAlgorithmId: i32,
    pub m_MaskType: ESubjectMaskingType,
    pub m_NeedsFilteringTranslation: bool,
    pub m_DbInitialized: bool,
    pub m_SeqDb: Option<Arc<CSeqDB>>,
}

/// NCBI C++: `ISearch` (`uniform_search.hpp`).
#[derive(Clone, Debug)]
pub struct ISearch;

/// NCBI C++: `ISeqSearch` (`uniform_search.hpp`).
#[derive(Clone, Debug)]
pub struct ISeqSearch;

/// NCBI C++: `IPssmSearch` (`uniform_search.hpp`).
#[derive(Clone, Debug)]
pub struct IPssmSearch;

/// NCBI C++: `ISearchFactory` (`uniform_search.hpp`).
#[derive(Clone, Debug)]
pub struct ISearchFactory;

#[derive(Clone, Debug)]
pub struct CSeqDBGiList;
