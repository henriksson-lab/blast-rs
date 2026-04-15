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
    fn min_seq_len(&self) -> i32 {
        1
    }

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
        fn num_seqs(&self) -> i32 {
            self.seqs.len() as i32
        }
        fn total_length(&self) -> i64 {
            self.seqs.iter().map(|s| s.len() as i64).sum()
        }
        fn max_seq_len(&self) -> i32 {
            self.seqs.iter().map(|s| s.len() as i32).max().unwrap_or(0)
        }
        fn avg_seq_len(&self) -> i32 {
            if self.seqs.is_empty() {
                0
            } else {
                (self.total_length() / self.num_seqs() as i64) as i32
            }
        }
        fn name(&self) -> &str {
            "mock"
        }
        fn is_protein(&self) -> bool {
            false
        }
        fn seq_len(&self, oid: i32) -> i32 {
            self.seqs[oid as usize].len() as i32
        }
        fn get_sequence(&self, arg: &GetSeqArg) -> Option<SeqData> {
            let seq = self.seqs.get(arg.oid as usize)?;
            Some(SeqData {
                sequence: seq.clone(),
                length: seq.len() as i32,
            })
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

        let data = src
            .get_sequence(&GetSeqArg {
                oid: 0,
                encoding: SeqEncoding::Protein,
            })
            .unwrap();
        assert_eq!(data.length, 4);
    }

    /// Create a sequence source, iterate all sequences via iter_oids, verify count matches num_seqs.
    #[test]
    fn test_seqsrc_iteration() {
        let src = MockSeqSrc {
            seqs: vec![
                vec![0, 1, 2, 3],
                vec![4, 5],
                vec![6, 7, 8],
                vec![9, 10, 11, 12, 13],
            ],
        };
        let oids: Vec<i32> = src.iter_oids().collect();
        assert_eq!(oids.len() as i32, src.num_seqs());
        assert_eq!(oids, vec![0, 1, 2, 3]);

        // Verify each OID yields valid sequence data
        for oid in &oids {
            let data = src.get_sequence(&GetSeqArg {
                oid: *oid,
                encoding: SeqEncoding::Nucleotide,
            });
            assert!(data.is_some(), "Failed to get sequence for OID {}", oid);
        }

        // Empty source should iterate zero times
        let empty_src = MockSeqSrc { seqs: vec![] };
        let empty_oids: Vec<i32> = empty_src.iter_oids().collect();
        assert_eq!(empty_oids.len(), 0);
        assert_eq!(empty_src.num_seqs(), 0);
    }

    /// Verify sequence lengths from the source match the actual data lengths.
    #[test]
    fn test_seqsrc_length_queries() {
        let seqs = vec![
            vec![0u8; 100],
            vec![0u8; 1],
            vec![0u8; 50],
            vec![0u8; 200],
            vec![0u8; 75],
        ];
        let src = MockSeqSrc { seqs };

        assert_eq!(src.seq_len(0), 100);
        assert_eq!(src.seq_len(1), 1);
        assert_eq!(src.seq_len(2), 50);
        assert_eq!(src.seq_len(3), 200);
        assert_eq!(src.seq_len(4), 75);
        assert_eq!(src.max_seq_len(), 200);
        assert_eq!(src.min_seq_len(), 1);
        // avg = (100+1+50+200+75)/5 = 85 (integer division of 426/5 = 85)
        assert_eq!(src.avg_seq_len(), 85);
    }

    /// Out-of-bounds OID access should return None from get_sequence.
    #[test]
    fn test_seqsrc_bounds_checking() {
        let src = MockSeqSrc {
            seqs: vec![vec![0, 1, 2], vec![3, 4]],
        };
        // Valid OIDs
        assert!(src
            .get_sequence(&GetSeqArg {
                oid: 0,
                encoding: SeqEncoding::Nucleotide
            })
            .is_some());
        assert!(src
            .get_sequence(&GetSeqArg {
                oid: 1,
                encoding: SeqEncoding::Nucleotide
            })
            .is_some());
        // Out of bounds
        assert!(src
            .get_sequence(&GetSeqArg {
                oid: 2,
                encoding: SeqEncoding::Nucleotide
            })
            .is_none());
        assert!(src
            .get_sequence(&GetSeqArg {
                oid: 99,
                encoding: SeqEncoding::Nucleotide
            })
            .is_none());
        // Negative OID (wraps to large usize, so should be None)
        assert!(src
            .get_sequence(&GetSeqArg {
                oid: -1,
                encoding: SeqEncoding::Nucleotide
            })
            .is_none());
    }

    /// Sum of all sequence lengths should match total_length.
    #[test]
    fn test_seqsrc_total_length() {
        let src = MockSeqSrc {
            seqs: vec![vec![0; 10], vec![0; 20], vec![0; 30], vec![0; 40]],
        };
        // Total should be 10 + 20 + 30 + 40 = 100
        assert_eq!(src.total_length(), 100);

        // Verify by summing individual lengths
        let sum: i64 = src.iter_oids().map(|oid| src.seq_len(oid) as i64).sum();
        assert_eq!(sum, src.total_length());

        // Empty source should have total length 0
        let empty = MockSeqSrc { seqs: vec![] };
        assert_eq!(empty.total_length(), 0);
    }
}
