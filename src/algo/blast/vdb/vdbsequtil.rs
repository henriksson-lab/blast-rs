// Upstream sources:
// - ncbi-blast-2.17.0+-src/c++/src/algo/blast/vdb/vdbsequtil.h
// - ncbi-blast-2.17.0+-src/c++/src/algo/blast/vdb/vdbsequtil.c

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VdbByteArray {
    pub data: Vec<u8>,
    pub size: usize,
    pub capacity: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VdbNuclDataFormat {
    Ncbistdaa,
    Ncbi2na,
    Ncbi4na,
    Iupacna,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VdbNuclDataRequest {
    pub format: VdbNuclDataFormat,
    pub oid: i32,
    pub offset: i64,
    pub length: i64,
    pub reverse: bool,
    pub buffer: VdbByteArray,
}
