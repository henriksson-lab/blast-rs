// Upstream sources:
// - ncbi-blast-2.17.0+-src/c++/include/algo/blast/dbindex/dbindex.hpp
// - ncbi-blast-2.17.0+-src/c++/src/algo/blast/dbindex/dbindex.cpp

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IndexSuperHeaderErrCode {
    File,
    Read,
    Write,
    Endian,
    Version,
    Size,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Endianness {
    LittleEndian,
    BigEndian,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct IndexSuperHeaderBase {
    pub actual_size: usize,
    pub endianness: u32,
    pub version: u32,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct IndexSuperHeaderV1 {
    pub base: IndexSuperHeaderBase,
    pub num_seq: u32,
    pub num_vol: u32,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct IndexHeader {
    pub legacy: bool,
    pub hkey_width: u64,
    pub stride: u64,
    pub ws_hint: u64,
    pub max_chunk_size: u64,
    pub chunk_overlap: u64,
    pub start: u64,
    pub start_chunk: u64,
    pub stop: u64,
    pub stop_chunk: u64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VectorWrap<T> {
    pub base: Option<usize>,
    pub data: Vec<T>,
    pub vec: bool,
    pub size: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DbIndexErrCode {
    BadOption,
    BadSequence,
    BadVersion,
    BadData,
    Io,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DbIndexOptions {
    pub idmap: bool,
    pub legacy: bool,
    pub stride: u64,
    pub ws_hint: u64,
    pub hkey_width: u64,
    pub chunk_size: u64,
    pub chunk_overlap: u64,
    pub report_level: u64,
    pub max_index_size: u64,
    pub stat_file_name: String,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DbIndexSearchResults {
    pub word_size: u64,
    pub start: u64,
    pub results: Vec<Option<String>>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DbIndexVolInfo {
    pub mapfile: Option<String>,
    pub map: Option<usize>,
    pub map_start: Option<usize>,
    pub offset_data: Option<String>,
    pub subject_map_offset: usize,
    pub version: u64,
    pub stride: u64,
}
