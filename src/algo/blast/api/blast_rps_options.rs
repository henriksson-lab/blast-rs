use crate::algo::blast::api::blast_prot_options::CBlastProteinOptionsHandle;

#[derive(Clone, Debug)]
pub struct CBlastRPSOptionsHandle {
    pub base: CBlastProteinOptionsHandle,
    pub word_threshold: f64,
    pub word_size: i32,
    pub composition_based_stats: bool,
    pub gap_opening_cost: i32,
    pub gap_extension_cost: i32,
}
