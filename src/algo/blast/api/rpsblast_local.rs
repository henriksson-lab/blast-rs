use crate::algo::blast::api::blast_options_handle::CBlastOptionsHandle;
use crate::algo::blast::api::query_data::CBlastQueryVector;

pub struct CLocalRPSBlast {
    pub num_of_threads: u32,
    pub db_name: String,
    pub opt_handle: Option<Box<CBlastOptionsHandle>>,
    pub query_vector: Option<Box<CBlastQueryVector>>,
    pub num_of_dbs: u32,
    pub rps_databases: Vec<String>,
}
