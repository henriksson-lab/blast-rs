#[derive(Clone, Debug)]
pub struct CBlastUsageReport {
    pub app_name: String,
    pub program: String,
    pub task: String,
    pub database: String,
    pub query_count: usize,
    pub enabled: bool,
    pub sent: bool,
}
