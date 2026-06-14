// Upstream sources:
// - ncbi-blast-2.17.0+-src/c++/include/algo/blast/format/blast_async_format.hpp
// - ncbi-blast-2.17.0+-src/c++/src/algo/blast/format/blast_async_format.cpp

use std::collections::BTreeMap;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FormatResultValues {
    pub query_vector: Option<String>,
    pub blast_results: Option<String>,
    pub formatter: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BlastAsyncFormatThread {
    pub results_map: BTreeMap<i32, Vec<FormatResultValues>>,
    pub done: bool,
    pub semaphore: Option<String>,
}
