//! Rust equivalent of blast_lookup.c, blast_nalookup.c, lookup_wrap.c
//! Lookup table structures for BLAST word finding.

/// Lookup table types.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LookupTableType {
    MegablastLookup,
    SmallNaLookup,
    NaLookup,
    AaLookup,
    CompressedAaLookup,
    RpsLookup,
    PhiLookup,
    PhiNaLookup,
    IndexedMbLookup,
    NaHashLookup,
}

/// A (query_offset, subject_offset) pair from the scan phase.
#[derive(Debug, Clone, Copy)]
pub struct OffsetPair {
    pub query_offset: i32,
    pub subject_offset: i32,
}

/// Nucleotide lookup table for standard blastn (word_size 4-12).
#[derive(Debug)]
pub struct SmallNaLookupTable {
    pub word_length: i32,
    pub backbone: Vec<i32>,
    pub overflow: Vec<i32>,
    pub pv_array: Vec<u32>, // presence vector for quick filtering
    pub scan_step: i32,
}

/// Megablast lookup table (word_size >= 12, contiguous).
#[derive(Debug)]
pub struct MbLookupTable {
    pub word_length: i32,
    pub lut_word_length: i32, // may differ from word_length for discontiguous
    pub hashtable: Vec<i32>,
    pub next_pos: Vec<i32>,
    pub pv_array: Vec<u32>,
    pub longest_chain: i32,
    pub scan_step: i32,
}

/// Protein lookup table.
#[derive(Debug)]
pub struct AaLookupTable {
    pub word_length: i32,
    pub threshold: f64,
    pub backbone: Vec<Vec<i32>>, // backbone[word_hash] = list of query offsets
    pub pv_array: Vec<u32>,
}

/// Generic wrapper around different lookup table types.
pub enum LookupTableWrap {
    SmallNa(SmallNaLookupTable),
    Megablast(MbLookupTable),
    Aa(AaLookupTable),
}

impl LookupTableWrap {
    pub fn table_type(&self) -> LookupTableType {
        match self {
            LookupTableWrap::SmallNa(_) => LookupTableType::SmallNaLookup,
            LookupTableWrap::Megablast(_) => LookupTableType::MegablastLookup,
            LookupTableWrap::Aa(_) => LookupTableType::AaLookup,
        }
    }

    pub fn word_length(&self) -> i32 {
        match self {
            LookupTableWrap::SmallNa(t) => t.word_length,
            LookupTableWrap::Megablast(t) => t.word_length,
            LookupTableWrap::Aa(t) => t.word_length,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_offset_pair() {
        let pair = OffsetPair {
            query_offset: 10,
            subject_offset: 20,
        };
        assert_eq!(pair.query_offset, 10);
        assert_eq!(pair.subject_offset, 20);
    }

    #[test]
    fn test_lookup_type() {
        let table = LookupTableWrap::Megablast(MbLookupTable {
            word_length: 28,
            lut_word_length: 12,
            hashtable: vec![],
            next_pos: vec![],
            pv_array: vec![],
            longest_chain: 0,
            scan_step: 1,
        });
        assert_eq!(table.table_type(), LookupTableType::MegablastLookup);
        assert_eq!(table.word_length(), 28);
    }
}
