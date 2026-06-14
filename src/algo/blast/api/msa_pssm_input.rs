pub struct CSeqEntry;
use crate::algo::blast::api::query_data::CBioseq;

pub struct PsiMsa;
pub struct PsiMsaDimensions;
pub struct PsiBlastOptions;
pub struct PsiDiagnosticsRequest;

pub struct CPsiBlastInputClustalW {
    pub query: Vec<u8>,
    pub ascii_msa: Vec<String>,
    pub msa: *mut PsiMsa,
    pub msa_dimensions: PsiMsaDimensions,
    pub opts: PsiBlastOptions,
    pub diagnostics_request: *mut PsiDiagnosticsRequest,
    pub matrix_name: String,
    pub gap_existence: i32,
    pub gap_extension: i32,
    pub seq_entry: Option<Box<CSeqEntry>>,
    pub query_bioseq: Option<Box<CBioseq>>,
}
