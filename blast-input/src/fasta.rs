//! Multi-FASTA parser.

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
pub fn parse_fasta<R: Read>(reader: R) -> Vec<FastaRecord> {
    let buf = BufReader::new(reader);
    let mut records = Vec::new();
    let mut current_id = String::new();
    let mut current_defline = String::new();
    let mut current_seq = Vec::new();

    for line in buf.lines() {
        let line = match line {
            Ok(l) => l,
            Err(_) => break,
        };
        let line = line.trim_end();
        if line.is_empty() {
            continue;
        }
        if let Some(header) = line.strip_prefix('>') {
            // Save previous record if any
            if !current_defline.is_empty() || !current_seq.is_empty() {
                records.push(FastaRecord {
                    id: current_id,
                    defline: current_defline,
                    sequence: current_seq,
                });
            }
            current_defline = header.to_string();
            current_id = header
                .split_whitespace()
                .next()
                .unwrap_or("")
                .to_string();
            current_seq = Vec::new();
        } else {
            // Sequence line: append bytes, skip whitespace
            for &b in line.as_bytes() {
                if !b.is_ascii_whitespace() {
                    current_seq.push(b);
                }
            }
        }
    }

    // Don't forget the last record
    if !current_defline.is_empty() || !current_seq.is_empty() {
        records.push(FastaRecord {
            id: current_id,
            defline: current_defline,
            sequence: current_seq,
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
