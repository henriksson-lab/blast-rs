use std::sync::Arc;

#[derive(Clone, Debug)]
pub struct CIgBlastOptions {
    pub organism: String,
    pub domain_system: String,
    pub germline_db: String,
    pub auxiliary_data: String,
}

#[derive(Clone, Debug)]
pub struct CIgBlast {
    pub options: Option<Arc<CIgBlastOptions>>,
    pub rearrangement_results: Vec<String>,
}
