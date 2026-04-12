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

    #[test]
    fn test_small_na_lookup_table_word_length() {
        let table = SmallNaLookupTable {
            word_length: 8,
            backbone: vec![-1; 65536], // 4^8 = 65536 entries
            overflow: vec![],
            pv_array: vec![],
            scan_step: 4,
        };
        assert_eq!(table.word_length, 8);
        assert_eq!(table.backbone.len(), 65536, "8-mer backbone should have 4^8=65536 entries");
        assert_eq!(table.scan_step, 4);
    }

    #[test]
    fn test_small_na_lookup_backbone_size() {
        // For word_size=7: 4^7 = 16384 backbone entries
        let table = SmallNaLookupTable {
            word_length: 7,
            backbone: vec![-1; 16384],
            overflow: vec![],
            pv_array: vec![],
            scan_step: 4,
        };
        assert_eq!(table.backbone.len(), 16384);
    }

    #[test]
    fn test_megablast_lookup_table_properties() {
        let table = MbLookupTable {
            word_length: 28,
            lut_word_length: 12,
            hashtable: vec![0; 4_194_304], // 4^12 = 4194304
            next_pos: vec![],
            pv_array: vec![],
            longest_chain: 0,
            scan_step: 1,
        };
        assert_eq!(table.word_length, 28);
        assert_eq!(table.lut_word_length, 12);
        assert_eq!(table.hashtable.len(), 4_194_304, "Megablast hash size should be 4^12");
    }

    #[test]
    fn test_aa_lookup_table_properties() {
        let table = AaLookupTable {
            word_length: 3,
            threshold: 11.0,
            backbone: vec![vec![]; 28 * 28 * 28], // 28^3 for NCBIstdaa
            pv_array: vec![],
        };
        assert_eq!(table.word_length, 3);
        assert_eq!(table.threshold, 11.0);
        assert_eq!(table.backbone.len(), 21952, "3-mer AA backbone should have 28^3=21952 entries");
    }

    #[test]
    fn test_lookup_table_wrap_variants() {
        let small = LookupTableWrap::SmallNa(SmallNaLookupTable {
            word_length: 11, backbone: vec![], overflow: vec![], pv_array: vec![], scan_step: 4,
        });
        assert_eq!(small.table_type(), LookupTableType::SmallNaLookup);
        assert_eq!(small.word_length(), 11);

        let aa = LookupTableWrap::Aa(AaLookupTable {
            word_length: 3, threshold: 11.0, backbone: vec![], pv_array: vec![],
        });
        assert_eq!(aa.table_type(), LookupTableType::AaLookup);
        assert_eq!(aa.word_length(), 3);
    }

    #[test]
    fn test_pv_array_bit_operations() {
        // Simulate PV array: set bit for index 100, check it
        let mut pv = vec![0u32; 8]; // 256 bits
        let idx = 100usize;
        let word = idx / 32;
        let bit = idx % 32;
        pv[word] |= 1 << bit;
        assert_ne!(pv[word] & (1 << bit), 0, "PV bit should be set");
        assert_eq!(pv[3] & (1 << 4), 1 << 4, "Bit 100 = word 3, bit 4");
        // Unset bit should be 0
        assert_eq!(pv[0] & (1 << 0), 0, "Bit 0 should not be set");
    }
}
