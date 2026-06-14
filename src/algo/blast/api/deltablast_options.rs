use crate::algo::blast::api::blast_options_cxx::CBlastOptions;

pub struct CDeltaBlastOptionsHandle {
    pub opts: Option<Box<CBlastOptions>>,
}
