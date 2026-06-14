use crate::algo::blast::api::blast_options_cxx::CBlastOptions;

pub struct CPHIBlastNuclOptionsHandle {
    pub opts: Option<Box<CBlastOptions>>,
}
