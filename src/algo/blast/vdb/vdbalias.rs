// Upstream sources:
// - ncbi-blast-2.17.0+-src/c++/include/algo/blast/vdb/vdbalias.hpp
// - ncbi-blast-2.17.0+-src/c++/src/algo/blast/vdb/vdbalias.cpp

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VdbAliasNode {
    pub name: String,
    pub children: Vec<VdbAliasNode>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VdbAliasUtil {
    pub aliases: Vec<VdbAliasNode>,
}
