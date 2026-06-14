// Upstream source:
// - ncbi-blast-2.17.0+-src/c++/src/algo/blast/dbindex/makeindex/main.cpp

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MkIndexMain {
    pub application: Option<super::mkindex_app::MkIndexApplication>,
}
