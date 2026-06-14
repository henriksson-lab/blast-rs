use crate::algo::blast::api::blast_options_cxx::CBlastOptions;

pub struct CPHIBlastProtOptionsHandle {
    pub opts: Option<Box<CBlastOptions>>,
}
