use crate::algo::blast::api::blast_options_handle::CBlastOptionsHandle;

#[derive(Clone, Debug)]
pub struct CTBlastxOptionsHandle {
    pub base: CBlastOptionsHandle,
    pub strand_option: i32,
    pub db_genetic_code: i32,
    pub query_genetic_code: i32,
}
