use crate::algo::blast::api::blast_options_cxx::EProgram;
use crate::algo::blast::api::query_data::EBlastProgramType;

#[derive(Clone, Debug)]
pub struct CBlastOptionsMemento {
    pub program: EProgram,
    pub program_type: EBlastProgramType,
    pub program_name: String,
    pub service_name: String,
    pub defaults_mode: bool,
}
