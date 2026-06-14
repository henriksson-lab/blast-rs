use crate::algo::blast::api::blast_nucl_options::CBlastNucleotideOptionsHandle;

#[derive(Clone, Debug)]
pub struct CDiscNucleotideOptionsHandle {
    pub base: CBlastNucleotideOptionsHandle,
    pub mb_template_length: u8,
    pub mb_template_type: u8,
    pub word_size: i32,
}
