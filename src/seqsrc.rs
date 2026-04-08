//! Rust equivalent of blast_seqsrc.c — sequence source trait.
//! Replaces the C vtable-based BlastSeqSrc with a Rust trait.

/// Encoding for sequence retrieval.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SeqEncoding {
    /// NCBIstdaa for protein, NCBI2na packed for nucleotide (no sentinels)
    Protein,
    /// BLASTNA decoded with sentinel bytes
    Nucleotide,
    /// NCBI4na encoding
    Ncbi4na,
    /// NCBI2na encoding
    Ncbi2na,
}

/// Arguments for fetching a sequence.
pub struct GetSeqArg {
    pub oid: i32,
    pub encoding: SeqEncoding,
}

/// Fetched sequence data.
pub struct SeqData {
    pub sequence: Vec<u8>,
    pub length: i32,
}

/// The sequence source trait — Rust replacement for BlastSeqSrc vtable.
pub trait BlastSeqSource: Send + Sync {
    /// Number of sequences in the database.
    fn num_seqs(&self) -> i32;

    /// Total length of all sequences in bases/residues.
    fn total_length(&self) -> i64;

    /// Maximum sequence length.
    fn max_seq_len(&self) -> i32;

    /// Average sequence length.
    fn avg_seq_len(&self) -> i32;

    /// Minimum sequence length.
    fn min_seq_len(&self) -> i32 { 1 }

    /// Database name.
    fn name(&self) -> &str;

    /// Whether the database contains protein sequences.
    fn is_protein(&self) -> bool;

    /// Get the length of a specific sequence.
    fn seq_len(&self, oid: i32) -> i32;

    /// Fetch sequence data for the given OID with the specified encoding.
    fn get_sequence(&self, arg: &GetSeqArg) -> Option<SeqData>;

    /// Iterator over OIDs.
    fn iter_oids(&self) -> Box<dyn Iterator<Item = i32> + '_>;
}

/// EOF sentinel for sequence source iteration.
pub const SEQSRC_EOF: i32 = -1;

#[cfg(test)]
mod tests {
    use super::*;

    struct MockSeqSrc {
        seqs: Vec<Vec<u8>>,
    }

    impl BlastSeqSource for MockSeqSrc {
        fn num_seqs(&self) -> i32 { self.seqs.len() as i32 }
        fn total_length(&self) -> i64 { self.seqs.iter().map(|s| s.len() as i64).sum() }
        fn max_seq_len(&self) -> i32 { self.seqs.iter().map(|s| s.len() as i32).max().unwrap_or(0) }
        fn avg_seq_len(&self) -> i32 {
            if self.seqs.is_empty() { 0 } else { (self.total_length() / self.num_seqs() as i64) as i32 }
        }
        fn name(&self) -> &str { "mock" }
        fn is_protein(&self) -> bool { false }
        fn seq_len(&self, oid: i32) -> i32 { self.seqs[oid as usize].len() as i32 }
        fn get_sequence(&self, arg: &GetSeqArg) -> Option<SeqData> {
            let seq = self.seqs.get(arg.oid as usize)?;
            Some(SeqData { sequence: seq.clone(), length: seq.len() as i32 })
        }
        fn iter_oids(&self) -> Box<dyn Iterator<Item = i32> + '_> {
            Box::new(0..self.num_seqs())
        }
    }

    #[test]
    fn test_mock_seqsrc() {
        let src = MockSeqSrc {
            seqs: vec![vec![0, 1, 2, 3], vec![0, 1]],
        };
        assert_eq!(src.num_seqs(), 2);
        assert_eq!(src.total_length(), 6);
        assert_eq!(src.max_seq_len(), 4);
        assert_eq!(src.seq_len(0), 4);
        assert_eq!(src.seq_len(1), 2);

        let data = src.get_sequence(&GetSeqArg { oid: 0, encoding: SeqEncoding::Protein }).unwrap();
        assert_eq!(data.length, 4);
    }
}
