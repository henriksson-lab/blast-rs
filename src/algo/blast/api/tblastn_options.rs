use crate::algo::blast::api::blast_advprot_options::CBlastAdvancedProteinOptionsHandle;

#[derive(Clone, Debug)]
pub struct CTBlastnOptionsHandle {
    pub base: CBlastAdvancedProteinOptionsHandle,
    pub out_of_frame_mode: bool,
    pub frame_shift_penalty: i32,
    pub longest_intron_length: i32,
    pub db_genetic_code: i32,
}
