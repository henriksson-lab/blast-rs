use crate::algo::blast::api::blast_prot_options::CBlastProteinOptionsHandle;

#[derive(Clone, Debug)]
pub struct CBlastAdvancedProteinOptionsHandle {
    pub base: CBlastProteinOptionsHandle,
    pub composition_based_stats: i32,
    pub smith_waterman_mode: bool,
}
