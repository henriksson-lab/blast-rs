#[derive(Clone, Debug)]
pub struct BlastKmerResults {
    pub query_id: String,
    pub subject_ids: Vec<String>,
    pub scores: Vec<f64>,
}
