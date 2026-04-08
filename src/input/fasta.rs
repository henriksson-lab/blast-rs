//! Multi-FASTA parser using noodles-fasta.

use std::io::{BufRead, BufReader, Read};

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
}
