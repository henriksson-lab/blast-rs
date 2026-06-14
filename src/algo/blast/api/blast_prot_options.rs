use crate::algo::blast::api::blast_options_handle::CBlastOptionsHandle;

#[derive(Clone, Debug)]
pub struct CBlastProteinOptionsHandle {
    pub base: CBlastOptionsHandle,
    pub word_threshold: f64,
    pub word_size: i32,
    pub x_dropoff: f64,
    pub seg_filtering: bool,
    pub matrix_name: String,
    pub gap_opening_cost: i32,
    pub gap_extension_cost: i32,
    pub chaining: bool,
}
