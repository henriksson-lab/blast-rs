use std::sync::Arc;

use crate::algo::blast::core::blast_query_info::BlastQueryInfo;
use crate::algo::blast::core::blast_util::BLAST_SequenceBlk;

pub use crate::algo::blast::api::blast_options_cxx::CBlastOptions;

pub const BLAST_GENETIC_CODE: u32 = 1;

pub const PROTEIN_QUERY_MASK: u32 = 0x1 << 0;
pub const PROTEIN_SUBJECT_MASK: u32 = 0x1 << 1;
pub const NUCLEOTIDE_QUERY_MASK: u32 = 0x1 << 2;
pub const NUCLEOTIDE_SUBJECT_MASK: u32 = 0x1 << 3;
pub const TRANSLATED_QUERY_MASK: u32 = 0x1 << 4;
pub const TRANSLATED_SUBJECT_MASK: u32 = 0x1 << 5;
pub const PSSM_QUERY_MASK: u32 = 0x1 << 6;
pub const PSSM_SUBJECT_MASK: u32 = 0x1 << 7;
pub const PATTERN_QUERY_MASK: u32 = 0x1 << 8;
pub const MAPPING_MASK: u32 = 0x1 << 9;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u32)]
pub enum EBlastProgramType {
    Blastp = PROTEIN_QUERY_MASK | PROTEIN_SUBJECT_MASK,
    Blastn = NUCLEOTIDE_QUERY_MASK | NUCLEOTIDE_SUBJECT_MASK,
    Blastx = NUCLEOTIDE_QUERY_MASK | PROTEIN_SUBJECT_MASK | TRANSLATED_QUERY_MASK,
    Tblastn = PROTEIN_QUERY_MASK | NUCLEOTIDE_SUBJECT_MASK | TRANSLATED_SUBJECT_MASK,
    Tblastx = NUCLEOTIDE_QUERY_MASK
        | NUCLEOTIDE_SUBJECT_MASK
        | TRANSLATED_QUERY_MASK
        | TRANSLATED_SUBJECT_MASK,
    PsiBlast = PSSM_QUERY_MASK | PROTEIN_QUERY_MASK | PROTEIN_SUBJECT_MASK,
    PsiTblastn =
        PSSM_QUERY_MASK | PROTEIN_QUERY_MASK | NUCLEOTIDE_SUBJECT_MASK | TRANSLATED_SUBJECT_MASK,
    RpsBlast = PSSM_SUBJECT_MASK | PROTEIN_QUERY_MASK | PROTEIN_SUBJECT_MASK,
    RpsTblastn =
        PSSM_SUBJECT_MASK | NUCLEOTIDE_QUERY_MASK | PROTEIN_SUBJECT_MASK | TRANSLATED_QUERY_MASK,
    PhiBlastp = PATTERN_QUERY_MASK | PROTEIN_QUERY_MASK | PROTEIN_SUBJECT_MASK,
    PhiBlastn = PATTERN_QUERY_MASK | NUCLEOTIDE_QUERY_MASK | NUCLEOTIDE_SUBJECT_MASK,
    Mapping = NUCLEOTIDE_QUERY_MASK | NUCLEOTIDE_SUBJECT_MASK | MAPPING_MASK,
    Undefined = 0,
}

#[derive(Debug, Clone)]
pub struct CSeqLoc;

#[derive(Debug, Clone)]
pub struct CScope;

#[derive(Debug, Clone)]
pub struct CSeqLocInfo;

#[derive(Debug, Clone)]
pub struct CBioseqSet;

#[derive(Debug, Clone)]
pub struct CBioseq;

#[derive(Debug, Clone)]
pub struct CSeqData;

#[derive(Debug, Clone)]
pub struct CSeqId;

#[derive(Debug, Clone)]
pub struct CSeqAlign;

#[derive(Debug, Clone)]
pub struct CSeqAlignSet;

#[derive(Debug, Clone)]
pub struct CStdSeg;

#[derive(Debug, Clone)]
pub struct CSplicedSeg;

#[derive(Debug, Clone)]
pub struct CPssmWithParameters;

#[derive(Debug, Clone)]
pub struct CSeqVector;

pub type TMaskedQueryRegions = Vec<Arc<CSeqLocInfo>>;
pub type TSeqLocInfoVector = Vec<TMaskedQueryRegions>;
pub type TSeqLocVector = Vec<SSeqLoc>;

#[derive(Debug, Clone)]
pub struct SSeqLoc {
    pub seqloc: Option<Arc<CSeqLoc>>,
    pub scope: Option<Arc<CScope>>,
    pub mask: Option<Arc<CSeqLoc>>,
    pub ignore_strand_in_mask: bool,
    pub genetic_code_id: u32,
}

#[derive(Debug, Clone)]
pub struct CBlastSearchQuery {
    pub seqloc: Option<Arc<CSeqLoc>>,
    pub scope: Option<Arc<CScope>>,
    pub mask: TMaskedQueryRegions,
    pub genetic_code_id: u32,
}

#[derive(Debug, Clone)]
pub struct CBlastQueryVector {
    pub queries: Vec<Arc<CBlastSearchQuery>>,
}

#[derive(Debug, Clone)]
pub struct ILocalQueryData {
    pub seq_blk: Option<BLAST_SequenceBlk>,
    pub query_info: Option<BlastQueryInfo>,
    pub messages: TSearchMessages,
    pub sum_of_sequence_lengths: usize,
}

#[derive(Debug, Clone)]
pub struct IRemoteQueryData {
    pub bioseqs: Option<Arc<CBioseqSet>>,
    pub seq_locs: TRemoteSeqLocs,
}

#[derive(Debug, Clone)]
pub struct IQueryFactory {
    pub local_query_data: Option<Arc<ILocalQueryData>>,
    pub remote_query_data: Option<Arc<IRemoteQueryData>>,
}

#[derive(Debug, Clone)]
pub struct IBlastQuerySource {
    pub bioseq_set_source_present: bool,
}

pub type TRemoteSeqLocs = Vec<Arc<CSeqLoc>>;
pub type TSearchMessages = Vec<TQueryMessages>;
pub type TQueryMessages = Vec<CSearchMessage>;

#[derive(Debug, Clone)]
pub struct CSearchMessage {
    pub query_index: usize,
    pub severity: i32,
    pub message: String,
}
