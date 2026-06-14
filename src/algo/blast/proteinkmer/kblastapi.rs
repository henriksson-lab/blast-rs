use std::sync::Arc;

use crate::algo::blast::proteinkmer::blastkmerresults::BlastKmerResults;

#[derive(Clone, Debug)]
pub struct KBlastApi {
    pub database: String,
    pub query_file: String,
    pub results: Option<Arc<BlastKmerResults>>,
}
