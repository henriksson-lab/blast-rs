use crate::algo::blast::api::query_data::{CBioseq, CSeqAlignSet, CSeqId};
pub struct PsiCdMsa;
pub struct PsiCdMsaCell;
pub struct PsiCdMsaCellData;
pub struct PsiMsaDimensions;
pub struct PsiBlastOptions;
pub struct PsiDiagnosticsRequest;

pub struct IPssmInputCdd {
    pub base: IPssmInputBase,
}

pub struct IPssmInputBase {
    pub gap_existence: i32,
    pub gap_extension: i32,
}

pub struct CCddInputData {
    pub query_data: Vec<u8>,
    pub query_title: String,
    pub db_name: String,
    pub seqalign_set: Option<Box<CSeqAlignSet>>,
    pub hits: Vec<Box<CddHit>>,
    pub cdd_data: PsiCdMsa,
    pub msa_dimensions: PsiMsaDimensions,
    pub msa_data: Vec<PsiCdMsaCell>,
    pub msa: *mut *mut PsiCdMsaCell,
    pub opts: PsiBlastOptions,
    pub matrix_name: String,
    pub diagnostics_request: *mut PsiDiagnosticsRequest,
    pub min_evalue: f64,
    pub query_bioseq: Option<Box<CBioseq>>,
    pub gap_existence: i32,
    pub gap_extension: i32,
}

pub struct CddRange {
    pub from: i32,
    pub to: i32,
}

pub struct CddHitSegment {
    pub query_range: CddRange,
    pub subject_range: CddRange,
    pub msa_data: Vec<PsiCdMsaCellData>,
    pub w_freqs_data: Vec<f64>,
}

pub enum CddHitApplyTo {
    Query,
    Subject,
}

pub struct CddHit {
    pub subject_id: Option<Box<CSeqId>>,
    pub evalue: f64,
    pub msa_idx: i32,
    pub segment_list: Vec<Box<CddHitSegment>>,
}
