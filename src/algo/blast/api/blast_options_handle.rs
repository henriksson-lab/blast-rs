use std::sync::Arc;

use crate::algo::blast::api::blast_options_cxx::{CBlastOptions, EApiLocality, EProgram};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ETaskSet {
    NuclNucl,
    ProtProt,
    Mapping,
    All,
}

#[derive(Clone, Debug)]
pub struct CBlastOptionsFactory {
    pub tasks: Vec<String>,
    pub documentation: Vec<(String, String)>,
    pub default_locality: EApiLocality,
}

#[derive(Clone, Debug)]
pub struct CBlastOptionsHandle {
    pub opts: Option<Arc<CBlastOptions>>,
    pub defaults_mode: bool,
    pub locality: EApiLocality,
    pub program: EProgram,
}
