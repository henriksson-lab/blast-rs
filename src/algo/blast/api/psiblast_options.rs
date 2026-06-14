use crate::algo::blast::api::blast_options_cxx::CBlastOptions;

pub struct CPSIBlastOptionsHandle {
    pub opts: Option<Box<CBlastOptions>>,
}
