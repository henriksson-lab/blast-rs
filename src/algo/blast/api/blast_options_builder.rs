use std::marker::PhantomData;

use super::blast_options_cxx::EApiLocality;
use super::uniform_search::ESubjectMaskingType;

/// NCBI C++: `CBlastOptionsBuilder::SOptional<T>` (`blast_options_builder.hpp`).
#[derive(Clone, Debug)]
pub struct SOptional<T> {
    pub m_IsSet: bool,
    pub m_Value: T,
}

/// NCBI C++: `CBlastOptionsBuilder` (`blast_options_builder.hpp`).
#[derive(Clone, Debug)]
pub struct CBlastOptionsBuilder {
    pub m_Program: String,
    pub m_Service: String,
    pub m_PerformCulling: bool,
    pub m_HspRangeMax: i32,
    pub m_EntrezQuery: SOptional<String>,
    pub m_FirstDbSeq: SOptional<i32>,
    pub m_FinalDbSeq: SOptional<i32>,
    pub m_GiList: SOptional<Vec<TGi>>,
    pub m_NegativeGiList: SOptional<Vec<TGi>>,
    pub m_DbFilteringAlgorithmId: SOptional<i32>,
    pub m_DbFilteringAlgorithmKey: SOptional<String>,
    pub m_SubjectMaskingType: SOptional<ESubjectMaskingType>,
    pub m_QueryMasks: SOptional<TMaskList>,
    pub m_IgnoreQueryMasks: bool,
    pub m_QueryRange: TSeqRange,
    pub m_Locality: EApiLocality,
    pub m_IgnoreUnsupportedOptions: bool,
    pub m_ForceMbIndex: bool,
    pub m_MbIndexName: String,
    pub m_TaxidList: SOptional<Vec<TTaxId>>,
    pub m_NegativeTaxidList: SOptional<Vec<TTaxId>>,
}

pub type TValueList = Vec<CBlast4Parameter>;
pub type TMaskList = Vec<CBlast4Mask>;
pub type TGi = i64;
pub type TTaxId = i32;

#[derive(Clone, Debug)]
pub struct TSeqRange {
    pub from: usize,
    pub to: usize,
}

#[derive(Clone, Debug)]
pub struct CBlast4Parameter {
    pub _opaque: PhantomData<()>,
}

#[derive(Clone, Debug)]
pub struct CBlast4Mask {
    pub _opaque: PhantomData<()>,
}
