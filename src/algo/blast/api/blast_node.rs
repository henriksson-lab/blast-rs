use std::sync::Arc;

use crate::algo::blast::api::blast_results::CSearchResultSet;

#[derive(Clone, Debug)]
pub struct CBlastNode {
    pub id: String,
    pub children: Vec<Arc<CBlastNode>>,
    pub results: Option<Arc<CSearchResultSet>>,
    pub messages: Vec<String>,
}
