use crate::algo::blast::api::blast_options_handle::CBlastOptionsHandle;

#[derive(Clone, Debug)]
pub struct CBlastNucleotideOptionsHandle {
    pub base: CBlastOptionsHandle,
    pub lookup_table_type: i32,
    pub word_size: i32,
    pub strand_option: i32,
    pub dust_filtering: bool,
    pub repeat_filtering: bool,
    pub window_masker_taxid: i32,
    pub window_masker_database: String,
    pub x_dropoff: f64,
    pub gap_extn_algorithm: i32,
    pub gap_traceback_algorithm: i32,
    pub match_reward: i32,
    pub mismatch_penalty: i32,
    pub matrix_name: String,
    pub gap_opening_cost: i32,
    pub gap_extension_cost: i32,
}
