use crate::algo::blast::api::blast_options_cxx::CBlastOptions;

pub struct CRPSTBlastnOptionsHandle {
    pub opts: Option<Box<CBlastOptions>>,
}
