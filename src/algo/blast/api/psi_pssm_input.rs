use crate::algo::blast::api::query_data::{CBioseq, CScope, CSeqAlignSet};

pub struct PsiMsa;
pub struct PsiMsaDimensions;
pub struct PsiBlastOptions;
pub struct PsiDiagnosticsRequest;

pub struct CPsiBlastInputFreqRatios {
    pub query: *mut u8,
    pub query_length: u32,
    pub matrix_name: Option<String>,
    pub gap_existence: i32,
    pub gap_extension: i32,
    pub freq_ratios: Vec<Vec<f64>>,
    pub impala_scale_factor: f64,
}

pub struct CPsiBlastInputData {
    pub query: *mut u8,
    pub query_title: String,
    pub scope: Option<Box<CScope>>,
    pub msa: *mut PsiMsa,
    pub msa_dimensions: PsiMsaDimensions,
    pub seq_align_set: Option<Box<CSeqAlignSet>>,
    pub opts: PsiBlastOptions,
    pub diagnostics_request: *mut PsiDiagnosticsRequest,
    pub matrix_name: String,
    pub gap_existence: i32,
    pub gap_extension: i32,
    pub query_bioseq: Option<Box<CBioseq>>,
}
