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

impl ILocalQueryData {
    pub fn x_validate_index(&self, index: usize) {
        let num_queries = self
            .query_info
            .as_ref()
            .map(|query_info| query_info.num_queries.max(0) as usize)
            .unwrap_or(0);
        if index > num_queries {
            panic!("Index {index} out of range ({num_queries} max)");
        }
    }

    pub fn is_valid_query(&self, index: usize) -> bool {
        self.x_validate_index(index);

        let query_info = self.query_info.as_ref().expect("BlastQueryInfo");
        let mut all_contexts_valid = true;
        for context_index in query_info.first_context..=query_info.last_context {
            let Some(context) = query_info.contexts.get(context_index as usize) else {
                continue;
            };
            if context.query_index == index as i32 && !context.is_valid {
                all_contexts_valid = false;
                break;
            }
        }
        all_contexts_valid
    }

    pub fn is_at_least_one_query_valid(&self) -> bool {
        let num_queries = self
            .query_info
            .as_ref()
            .map(|query_info| query_info.num_queries.max(0) as usize)
            .unwrap_or(0);
        for index in 0..num_queries {
            if self.is_valid_query(index) {
                return true;
            }
        }
        false
    }

    pub fn get_query_messages(&self, index: usize, qmsgs: &mut TQueryMessages) {
        self.x_validate_index(index);
        *qmsgs = self.messages[index].clone();
    }

    pub fn get_messages(&self, messages: &mut TSearchMessages) {
        *messages = self.messages.clone();
    }

    pub fn flush_sequence_data(&mut self) {
        self.seq_blk = None;
    }
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
