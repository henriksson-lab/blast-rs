// Upstream sources:
// - ncbi-blast-2.17.0+-src/c++/include/algo/blast/format/blastfmtutil.hpp
// - ncbi-blast-2.17.0+-src/c++/src/algo/blast/format/blastfmtutil.cpp

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SubjectScores {
    None,
    Score,
    BitScore,
    Evalue,
    Identity,
    Coverage,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BlastFormatUtilDbInfo {
    pub name: String,
    pub definition: String,
    pub date: String,
    pub num_sequences: i64,
    pub total_length: i64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BlastFormattingMatrix {
    pub matrix_name: String,
    pub scores: Vec<Vec<i32>>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BlastXmlIncremental {
    pub iteration_num: i32,
    pub serial_xml_end: String,
}
