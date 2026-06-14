// Upstream sources:
// - ncbi-blast-2.17.0+-src/c++/include/algo/blast/blastinput/blast_asn1_input.hpp
// - ncbi-blast-2.17.0+-src/c++/src/algo/blast/blastinput/blast_asn1_input.cpp

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Asn1InputSourceOmf {
    pub bases_added: u64,
    pub input_stream: Option<String>,
    pub second_input_stream: Option<String>,
    pub is_paired: bool,
    pub is_binary: bool,
}
