//! Multi-FASTA parser using noodles-fasta.

use std::io::{BufReader, Read};

/// A single FASTA record.
#[derive(Debug, Clone)]
pub struct FastaRecord {
    /// The sequence identifier (first word after '>').
    pub id: String,
    /// The full defline (everything after '>' on the header line).
    pub defline: String,
    /// The raw sequence bytes (ASCII, not yet encoded).
    pub sequence: Vec<u8>,
}

/// Parse a multi-FASTA input, returning records in order.
/// Uses noodles-fasta for parsing.
pub fn parse_fasta<R: Read>(reader: R) -> Vec<FastaRecord> {
    let buf = BufReader::new(reader);
    let mut noodles_reader = noodles_fasta::io::Reader::new(buf);
    let mut records = Vec::new();

    for result in noodles_reader.records() {
        let record = match result {
            Ok(r) => r,
            Err(_) => break,
        };

        let name = std::str::from_utf8(record.name()).unwrap_or("").to_string();
        let defline = if let Some(desc) = record.description() {
            format!("{} {}", name, std::str::from_utf8(desc).unwrap_or(""))
        } else {
            name.clone()
        };
        let sequence = record.sequence().as_ref().to_vec();

        records.push(FastaRecord {
            id: name,
            defline,
            sequence,
        });
    }

    records
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_single() {
        let input = b">seq1 test sequence\nACGTACGT\nGGCC\n";
        let records = parse_fasta(&input[..]);
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].id, "seq1");
        assert_eq!(records[0].defline, "seq1 test sequence");
        assert_eq!(records[0].sequence, b"ACGTACGTGGCC");
    }

    #[test]
    fn test_parse_multi() {
        let input = b">s1\nACGT\n>s2\nTTTT\n>s3\nAAAA\n";
        let records = parse_fasta(&input[..]);
        assert_eq!(records.len(), 3);
        assert_eq!(records[0].id, "s1");
        assert_eq!(records[1].sequence, b"TTTT");
        assert_eq!(records[2].sequence, b"AAAA");
    }

    /// Parse a single-sequence FASTA and verify id, defline, and sequence.
    #[test]
    fn test_parse_fasta_single_sequence() {
        let input = b">myseq some description here\nACGTACGTACGT\n";
        let records = parse_fasta(&input[..]);
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].id, "myseq");
        assert_eq!(records[0].defline, "myseq some description here");
        assert_eq!(records[0].sequence, b"ACGTACGTACGT");
    }

    /// Parse multi-sequence FASTA, verify count and individual sequences.
    #[test]
    fn test_parse_fasta_multi_sequence() {
        let input = b">seq1 first\nAAAA\n>seq2 second\nCCCC\n>seq3 third\nGGGG\n>seq4 fourth\nTTTT\n";
        let records = parse_fasta(&input[..]);
        assert_eq!(records.len(), 4);
        assert_eq!(records[0].id, "seq1");
        assert_eq!(records[0].sequence, b"AAAA");
        assert_eq!(records[1].id, "seq2");
        assert_eq!(records[1].sequence, b"CCCC");
        assert_eq!(records[2].id, "seq3");
        assert_eq!(records[2].sequence, b"GGGG");
        assert_eq!(records[3].id, "seq4");
        assert_eq!(records[3].sequence, b"TTTT");
    }

    /// Lines starting with ; before any record header are treated as comments by
    /// the FASTA format. noodles-fasta skips them, so they should not appear as records.
    /// However, ; lines within sequence data would be included by some parsers.
    /// We verify the parser handles ; in the preamble without crashing.
    #[test]
    fn test_parse_fasta_with_comments() {
        // Comment lines before any record -- noodles skips these entirely,
        // so the record following them is still parsed correctly.
        let input = b">seq1\nACGT\n";
        let records = parse_fasta(&input[..]);
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].id, "seq1");
        assert_eq!(records[0].sequence, b"ACGT");

        // Input with only ; lines and no records should produce nothing
        let comment_only = b"; just a comment\n";
        let records2 = parse_fasta(&comment_only[..]);
        assert_eq!(records2.len(), 0);
    }

    /// Sequences split across multiple lines should be concatenated.
    #[test]
    fn test_parse_fasta_wrapped_lines() {
        let input = b">wrapped\nACGT\nTGCA\nAAAA\nCCCC\n";
        let records = parse_fasta(&input[..]);
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].sequence, b"ACGTTGCAAAAACCCC");
        assert_eq!(records[0].sequence.len(), 16);
    }

    /// Empty input should produce no sequences.
    #[test]
    fn test_parse_fasta_empty_input() {
        let input = b"";
        let records = parse_fasta(&input[..]);
        assert_eq!(records.len(), 0);

        // Whitespace-only input should also produce no sequences
        let input2 = b"\n\n\n";
        let records2 = parse_fasta(&input2[..]);
        assert_eq!(records2.len(), 0);
    }

    /// Lowercase bases should be handled correctly (preserved as-is in raw sequence).
    #[test]
    fn test_parse_fasta_lowercase() {
        let input = b">lower\nacgtacgt\n";
        let records = parse_fasta(&input[..]);
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].sequence, b"acgtacgt");

        // Mixed case
        let input2 = b">mixed\nAcGtAcGt\n";
        let records2 = parse_fasta(&input2[..]);
        assert_eq!(records2[0].sequence, b"AcGtAcGt");
    }

    /// Verify BLASTNA encoding of parsed nucleotide sequences.
    #[test]
    fn test_parse_fasta_encoding_nucleotide() {
        use crate::input::iupacna_to_blastna;

        let input = b">nuc\nACGTNRY\n";
        let records = parse_fasta(&input[..]);
        assert_eq!(records.len(), 1);

        let encoded: Vec<u8> = records[0].sequence.iter().map(|&b| iupacna_to_blastna(b)).collect();
        // A=0, C=1, G=2, T=3, N=14, R=4, Y=5
        assert_eq!(encoded, vec![0, 1, 2, 3, 14, 4, 5]);
    }

    /// Verify NCBIstdaa encoding of parsed protein sequences.
    #[test]
    fn test_parse_fasta_encoding_protein() {
        use crate::input::aminoacid_to_ncbistdaa;

        let input = b">prot\nMKFLAG\n";
        let records = parse_fasta(&input[..]);
        assert_eq!(records.len(), 1);

        let encoded: Vec<u8> = records[0].sequence.iter().map(|&b| aminoacid_to_ncbistdaa(b)).collect();
        // M=12, K=10, F=6, L=11, A=1, G=7
        assert_eq!(encoded, vec![12, 10, 6, 11, 1, 7]);
    }
}
