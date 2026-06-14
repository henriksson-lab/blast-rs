use crate::algo::blast::api::blast_advprot_options::CBlastAdvancedProteinOptionsHandle;

#[derive(Clone, Debug)]
pub struct CBlastxOptionsHandle {
    pub base: CBlastAdvancedProteinOptionsHandle,
    pub strand_option: i32,
    pub query_genetic_code: i32,
    pub out_of_frame_mode: bool,
    pub frame_shift_penalty: i32,
    pub longest_intron_length: i32,
}
