use std::collections::BTreeSet;
use std::sync::Arc;

use super::blast_options_builder::CBlast4Parameter;
use super::blast_options_handle::CBlastOptionsHandle;
use super::query_data::{CBioseq, CPssmWithParameters, CSeqData, CSeqLoc};
use super::query_data::{EBlastProgramType, NUCLEOTIDE_QUERY_MASK, TRANSLATED_QUERY_MASK};

/// NCBI C++: `CRemoteBlastException::EErrCode` (`remote_blast.hpp`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum CRemoteBlastExceptionErrCode {
    EServiceNotAvailable,
    EIncompleteConfig,
}

/// NCBI C++: `CRemoteBlast::ESearchStatus` (`remote_blast.hpp`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ESearchStatus {
    EStatusUnknown,
    EStatusDone,
    EStatusPending,
    EStatusFailed,
}

/// NCBI C++: `CRemoteBlast::EDebugMode` (`remote_blast.hpp`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum EDebugMode {
    EDebug,
    ESilent,
}

/// NCBI C++: `CRemoteBlast::EState` (`remote_blast.hpp`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum EState {
    EStart,
    EFailed,
    EWait,
    EDone,
}

/// NCBI C++: `CRemoteBlast::EImmediacy` (`remote_blast.hpp`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum EImmediacy {
    EPollAsync,
    EPollImmed,
}

/// NCBI C++: `CRemoteBlast::ENeedConfig` (`remote_blast.hpp`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ENeedConfig {
    ENoConfig,
    EProgram,
    EService,
    EQueries,
    ESubject,
    ENeedAll,
}

/// NCBI C++: `CRemoteBlast` (`remote_blast.hpp`).
#[derive(Clone, Debug)]
pub struct CRemoteBlast {
    pub m_CBOH: Option<Arc<CBlastOptionsHandle>>,
    pub m_QSR: Option<Arc<CBlast4QueueSearchRequest>>,
    pub m_Reply: Option<Arc<CBlast4Reply>>,
    pub m_Archive: Option<Arc<CBlast4Archive>>,
    pub m_ReadFile: bool,
    pub m_ObjectStream: Option<Arc<CObjectIStream>>,
    pub m_ObjectType: EFormat,
    pub m_Errs: Vec<String>,
    pub m_Warn: Vec<String>,
    pub m_RID: String,
    pub m_ErrIgn: i32,
    pub m_Pending: bool,
    pub m_Verbose: EDebugMode,
    pub m_NeedConfig: ENeedConfig,
    pub m_Dbs: Option<Arc<CBlast4Database>>,
    pub m_SubjectSequences: Vec<Arc<CBioseq>>,
    pub m_SubjectSeqLocs: TSeqLocList,
    pub m_DbFilteringAlgorithmId: i32,
    pub m_DbFilteringAlgorithmKey: String,
    pub m_Task: String,
    pub m_ClientId: String,
    pub m_use_disk_cache: bool,
    pub m_disk_cache_error_flag: bool,
    pub m_disk_cache_error_msg: String,
    pub m_TaxidList: BTreeSet<TTaxId>,
    pub m_NegativeTaxidList: BTreeSet<TTaxId>,
}

pub type TSeqLocList = Vec<Arc<CSeqLoc>>;
pub type TKarlinAltschulBlocks = Vec<Arc<CBlast4KaBlock>>;
pub type TSeqIntervalVector = Vec<Arc<CSeqInterval>>;
pub type TSeqDataVector = Vec<Arc<CSeqData>>;
pub type TValueList = Vec<Arc<CBlast4Parameter>>;
pub type TTaxId = i32;

#[derive(Clone, Debug)]
pub struct CBlast4QueueSearchRequest;

#[derive(Clone, Debug)]
pub struct CBlast4Reply;

#[derive(Clone, Debug)]
pub struct CBlast4Archive;

#[derive(Clone, Debug)]
pub struct CObjectIStream;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum EFormat {
    EUnknown,
    EBinary,
    EText,
    EXml,
}

#[derive(Clone, Debug)]
pub struct CBlast4Database;

#[derive(Clone, Debug)]
pub struct CBlast4KaBlock;

#[derive(Clone, Debug)]
pub struct CSeqInterval;

#[derive(Clone, Debug)]
pub struct CRemoteBlastPssmQuery {
    pub m_Pssm: Option<Arc<CPssmWithParameters>>,
}

/// NCBI C++: `objects::EBlast4_frame_type`.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum EBlast4FrameType {
    NotSet,
    Plus1,
    Plus2,
    Plus3,
    Minus1,
    Minus2,
    Minus3,
}

/// NCBI C++: `CSeqLocInfo::ETranslationFrame`.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ETranslationFrame {
    NotSet,
    Plus1,
    Plus2,
    Plus3,
    Minus1,
    Minus2,
    Minus3,
}

/// NCBI C++: `FrameNumber2NetworkFrame`.
pub fn frame_number_2_network_frame(frame: i32, program: EBlastProgramType) -> EBlast4FrameType {
    let program = program as u32;
    if (program & TRANSLATED_QUERY_MASK) != 0 {
        return match frame {
            1 => EBlast4FrameType::Plus1,
            2 => EBlast4FrameType::Plus2,
            3 => EBlast4FrameType::Plus3,
            -1 => EBlast4FrameType::Minus1,
            -2 => EBlast4FrameType::Minus2,
            -3 => EBlast4FrameType::Minus3,
            _ => panic!("invalid translated query frame"),
        };
    }

    if (program & NUCLEOTIDE_QUERY_MASK) != 0 {
        debug_assert!(frame == -1 || frame == 1);
        return EBlast4FrameType::NotSet;
    }

    EBlast4FrameType::NotSet
}

/// NCBI C++: `NetworkFrame2FrameNumber`.
pub fn network_frame_2_frame_number(
    frame: EBlast4FrameType,
    program: EBlastProgramType,
) -> ETranslationFrame {
    let program = program as u32;
    if (program & TRANSLATED_QUERY_MASK) != 0 {
        return match frame {
            EBlast4FrameType::Plus1 => ETranslationFrame::Plus1,
            EBlast4FrameType::Plus2 => ETranslationFrame::Plus2,
            EBlast4FrameType::Plus3 => ETranslationFrame::Plus3,
            EBlast4FrameType::Minus1 => ETranslationFrame::Minus1,
            EBlast4FrameType::Minus2 => ETranslationFrame::Minus2,
            EBlast4FrameType::Minus3 => ETranslationFrame::Minus3,
            EBlast4FrameType::NotSet => panic!("invalid translated network frame"),
        };
    }

    ETranslationFrame::NotSet
}
