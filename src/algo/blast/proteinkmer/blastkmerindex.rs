#[derive(Clone, Debug)]
pub struct BlastKmerIndex {
    pub word_size: usize,
    pub offsets: Vec<usize>,
    pub positions: Vec<u32>,
}
