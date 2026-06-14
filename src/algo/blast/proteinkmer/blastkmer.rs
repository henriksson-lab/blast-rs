#[derive(Clone, Debug)]
pub struct BlastKmer {
    pub word_size: usize,
    pub alphabet_size: usize,
    pub kmers: Vec<Vec<u8>>,
}
