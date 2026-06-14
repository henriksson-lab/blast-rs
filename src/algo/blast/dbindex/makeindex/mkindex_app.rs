// Upstream sources:
// - ncbi-blast-2.17.0+-src/c++/src/algo/blast/dbindex/makeindex/mkindex_app.hpp
// - ncbi-blast-2.17.0+-src/c++/src/algo/blast/dbindex/makeindex/mkindex_app.cpp

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MkIndexApplication {
    pub usage_line: String,
}
