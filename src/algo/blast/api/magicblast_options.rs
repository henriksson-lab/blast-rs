use crate::algo::blast::api::blast_nucl_options::CBlastNucleotideOptionsHandle;

#[derive(Clone, Debug)]
pub struct CMagicBlastOptionsHandle {
    pub base: CBlastNucleotideOptionsHandle,
    pub lookup_db_filter: bool,
    pub max_db_word_count: u8,
    pub lookup_table_stride: i32,
    pub read_quality_filtering: bool,
    pub read_max_fraction_ambiguous: f64,
    pub min_dimer_entropy: i32,
    pub cutoff_score_coeffs: Vec<f64>,
    pub max_edit_distance: i32,
    pub paired: bool,
    pub splice_alignments: bool,
    pub longest_intron_length: i32,
}
