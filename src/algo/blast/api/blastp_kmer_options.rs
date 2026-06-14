use crate::algo::blast::api::blast_advprot_options::CBlastAdvancedProteinOptionsHandle;

#[derive(Clone, Debug)]
pub struct CBlastpKmerOptionsHandle {
    pub base: CBlastAdvancedProteinOptionsHandle,
    pub kmer_length: i32,
    pub kmer_count: i32,
    pub jaccard_threshold: f64,
}
