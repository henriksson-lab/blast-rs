// Upstream sources:
// - ncbi-blast-2.17.0+-src/c++/src/algo/blast/vdb/error_priv.h
// - ncbi-blast-2.17.0+-src/c++/src/algo/blast/vdb/error_priv.c

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VdbSrcErrCode {
    NoError,
    ManagerLoadError,
    ManagerFreeError,
    RunSetLoadError,
    RefSetLoadError,
    RunSetFreeError,
    UninitVdbDataError,
    AddRunsError,
    GetRunSetNameError,
    GetRunSetNumSeqError,
    Reader2naError,
    Reader4naError,
    Read2naCacheError,
    Read2naCopyError,
    Read4naCacheError,
    Read4naCopyError,
    NoMemForVdbData,
    NoMemForRuns,
    NoMemForChunkSeq,
    SeqLengthError,
    ReadIdMismatch,
    Seq4naStringError,
    Seq2naStringError,
    NumSeqOverflowError,
    FilteredRead,
    IdOutOfRange,
    Seq4naRefSeqBufOverflow,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VdbSrcErrMsg {
    pub is_error: bool,
    pub result_code: u32,
    pub local_code: VdbSrcErrCode,
    pub msg_context: Option<String>,
}
