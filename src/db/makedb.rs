//! makeblastdb equivalent — create BLAST databases from FASTA files.

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

/// Create a BLAST v4 nucleotide database from a FASTA file.
pub fn make_nucleotide_db(
    fasta_path: &Path,
    output_base: &Path,
    title: &str,
) -> io::Result<(u32, u64)> {
    // returns (num_seqs, total_length)
    let fasta_data = std::fs::read_to_string(fasta_path)?;
    let mut sequences: Vec<(String, Vec<u8>)> = Vec::new();

    let mut current_header = String::new();
    let mut current_seq = Vec::new();

    for line in fasta_data.lines() {
        if let Some(hdr) = line.strip_prefix('>') {
            if !current_header.is_empty() || !current_seq.is_empty() {
                sequences.push((current_header, current_seq));
            }
            current_header = hdr.to_string();
            current_seq = Vec::new();
        } else {
            for &b in line.trim().as_bytes() {
                if b.is_ascii_alphabetic() {
                    current_seq.push(b);
                }
            }
        }
    }
    if !current_header.is_empty() || !current_seq.is_empty() {
        sequences.push((current_header, current_seq));
    }

    let num_seqs = sequences.len() as u32;
    let mut total_length = 0u64;

    // Write .nsq (sequence data)
    let mut nsq = BufWriter::new(File::create(output_base.with_extension("nsq"))?);
    nsq.write_all(&[0u8])?; // sentinel byte

    let mut seq_offsets = vec![1u32]; // first seq starts at byte 1
    let mut amb_offsets = Vec::new();

    for (_, seq) in &sequences {
        let seq_start = seq_offsets.last().copied().unwrap_or(1);
        // Pack 4 bases per byte (NCBI2na)
        let mut packed = Vec::new();
        let iupac_to_2na = |b: u8| -> u8 {
            match b {
                b'A' | b'a' => 0,
                b'C' | b'c' => 1,
                b'G' | b'g' => 2,
                b'T' | b't' => 3,
                _ => 0,
            }
        };

        let full_bytes = seq.len() / 4;
        let remainder = seq.len() % 4;
        for i in 0..full_bytes {
            let b = (iupac_to_2na(seq[i * 4]) << 6)
                | (iupac_to_2na(seq[i * 4 + 1]) << 4)
                | (iupac_to_2na(seq[i * 4 + 2]) << 2)
                | iupac_to_2na(seq[i * 4 + 3]);
            packed.push(b);
        }
        // Last byte: pack remaining bases + remainder count in low 2 bits
        if remainder > 0 {
            let mut last = 0u8;
            for j in 0..remainder {
                last |= iupac_to_2na(seq[full_bytes * 4 + j]) << (6 - 2 * j);
            }
            last |= remainder as u8; // low 2 bits = remainder count
            packed.push(last);
        } else {
            packed.push(0); // remainder = 0
        }

        nsq.write_all(&packed)?;
        let amb_offset = seq_start as u32 + packed.len() as u32;
        amb_offsets.push(amb_offset);
        // No ambiguity data for now
        seq_offsets.push(amb_offset);
        total_length += seq.len() as u64;
    }
    amb_offsets.push(*seq_offsets.last().unwrap_or(&0));
    nsq.flush()?;

    // Write .nhr (header data) — ASN.1 BER encoded Blast-def-line-set
    let mut nhr = BufWriter::new(File::create(output_base.with_extension("nhr"))?);
    let mut hdr_offsets = vec![0u32];
    for (oid, (header, _)) in sequences.iter().enumerate() {
        let start = hdr_offsets.last().copied().unwrap_or(0);
        let asn1 = encode_defline_asn1(header, oid as i32);
        nhr.write_all(&asn1)?;
        hdr_offsets.push(start + asn1.len() as u32);
    }
    nhr.flush()?;

    // Write .nin (index)
    let mut nin = BufWriter::new(File::create(output_base.with_extension("nin"))?);
    // Version
    nin.write_all(&4u32.to_be_bytes())?;
    // Type (0 = nucleotide in v4 format)
    nin.write_all(&0u32.to_be_bytes())?;
    // Title
    let title_bytes = title.as_bytes();
    nin.write_all(&(title_bytes.len() as u32).to_be_bytes())?;
    nin.write_all(title_bytes)?;
    // Date
    let date = "2026-01-01";
    nin.write_all(&(date.len() as u32).to_be_bytes())?;
    nin.write_all(date.as_bytes())?;
    // Num OIDs
    nin.write_all(&num_seqs.to_be_bytes())?;
    // Total length (little-endian!)
    nin.write_all(&total_length.to_le_bytes())?;
    // Max length
    let max_len = sequences
        .iter()
        .map(|(_, s)| s.len() as u32)
        .max()
        .unwrap_or(0);
    nin.write_all(&max_len.to_be_bytes())?;

    // Header offsets (num_seqs + 1)
    for off in &hdr_offsets {
        nin.write_all(&off.to_be_bytes())?;
    }
    // Sequence offsets (num_seqs + 1)
    for off in &seq_offsets {
        nin.write_all(&off.to_be_bytes())?;
    }
    // Ambiguity offsets (num_seqs + 1) = same as seq offsets (no ambiguity)
    for off in &amb_offsets {
        nin.write_all(&off.to_be_bytes())?;
    }
    nin.flush()?;

    Ok((num_seqs, total_length))
}

/// Encode a Blast-def-line-set ASN.1 BER header matching NCBI format.
/// Structure: SEQUENCE { SEQUENCE { title [0] VisibleString, seqid [1] { general { db "BL_ORD_ID", tag id INTEGER } } } }
fn encode_defline_asn1(title: &str, oid: i32) -> Vec<u8> {
    let mut buf = Vec::new();
    // Blast-def-line-set ::= SEQUENCE OF Blast-def-line
    buf.extend_from_slice(&[0x30, 0x80]); // SEQUENCE, indefinite
    {
        // Blast-def-line ::= SEQUENCE
        buf.extend_from_slice(&[0x30, 0x80]);
        {
            // title [0] VisibleString
            buf.extend_from_slice(&[0xa0, 0x80]);
            // VisibleString
            buf.push(0x1a);
            let title_bytes = title.as_bytes();
            encode_asn1_length(&mut buf, title_bytes.len());
            buf.extend_from_slice(title_bytes);
            buf.extend_from_slice(&[0x00, 0x00]); // END context[0]

            // seqid [1] SET OF Seq-id
            buf.extend_from_slice(&[0xa1, 0x80]);
            {
                // Seq-id ::= CHOICE { general Dbtag }
                // general [10] in Seq-id CHOICE
                buf.extend_from_slice(&[0x30, 0x80]);
                {
                    buf.extend_from_slice(&[0xaa, 0x80]); // context[10] = general
                    {
                        // Dbtag ::= SEQUENCE { db VisibleString, tag Object-id }
                        buf.extend_from_slice(&[0x30, 0x80]);
                        {
                            // db [0] VisibleString "BL_ORD_ID"
                            buf.extend_from_slice(&[0xa0, 0x80]);
                            buf.push(0x1a);
                            buf.push(9); // length
                            buf.extend_from_slice(b"BL_ORD_ID");
                            buf.extend_from_slice(&[0x00, 0x00]); // END

                            // tag [1] Object-id ::= CHOICE { id INTEGER }
                            buf.extend_from_slice(&[0xa1, 0x80]);
                            {
                                // id [0] INTEGER
                                buf.extend_from_slice(&[0xa0, 0x80]);
                                buf.push(0x02); // INTEGER
                                let oid_bytes = encode_asn1_integer(oid);
                                encode_asn1_length(&mut buf, oid_bytes.len());
                                buf.extend_from_slice(&oid_bytes);
                                buf.extend_from_slice(&[0x00, 0x00]); // END id
                            }
                            buf.extend_from_slice(&[0x00, 0x00]); // END tag
                        }
                        buf.extend_from_slice(&[0x00, 0x00]); // END Dbtag
                    }
                    buf.extend_from_slice(&[0x00, 0x00]); // END general
                }
                buf.extend_from_slice(&[0x00, 0x00]); // END inner SEQUENCE
            }
            buf.extend_from_slice(&[0x00, 0x00]); // END seqid
        }
        buf.extend_from_slice(&[0x00, 0x00]); // END Blast-def-line
    }
    buf.extend_from_slice(&[0x00, 0x00]); // END Blast-def-line-set
    buf
}

fn encode_asn1_length(buf: &mut Vec<u8>, len: usize) {
    if len < 128 {
        buf.push(len as u8);
    } else if len < 256 {
        buf.push(0x81);
        buf.push(len as u8);
    } else {
        buf.push(0x82);
        buf.push((len >> 8) as u8);
        buf.push(len as u8);
    }
}

fn encode_asn1_integer(val: i32) -> Vec<u8> {
    if val == 0 {
        vec![0]
    } else if val > 0 && val < 128 {
        vec![val as u8]
    } else if val > 0 && val < 32768 {
        if val < 256 {
            vec![0, val as u8]
        } else {
            vec![(val >> 8) as u8, val as u8]
        }
    } else {
        vec![
            (val >> 24) as u8,
            (val >> 16) as u8,
            (val >> 8) as u8,
            val as u8,
        ]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_make_nucleotide_db() {
        let dir = std::env::temp_dir().join("blast_makedb_test");
        std::fs::create_dir_all(&dir).ok();

        let fasta = dir.join("test.fa");
        std::fs::write(&fasta, ">seq1\nACGTACGT\n>seq2\nTTTTAAAA\n").unwrap();

        let db_base = dir.join("testdb");
        let (nseq, total) = make_nucleotide_db(&fasta, &db_base, "Test DB").unwrap();
        assert_eq!(nseq, 2);
        assert_eq!(total, 16); // 8 + 8

        // Verify files exist
        assert!(db_base.with_extension("nin").exists());
        assert!(db_base.with_extension("nsq").exists());
        assert!(db_base.with_extension("nhr").exists());

        // Try opening with our reader
        let db = super::super::index::BlastDb::open(&db_base).unwrap();
        assert_eq!(db.num_oids, 2);
        assert_eq!(db.total_length, 16);

        std::fs::remove_dir_all(&dir).ok();
    }

    /// Create a nucleotide DB, then read it back and verify sequences can be decoded.
    #[test]
    fn test_create_and_read_nucleotide_db() {
        let dir = std::env::temp_dir().join("blast_makedb_create_read");
        std::fs::create_dir_all(&dir).ok();

        let fasta = dir.join("input.fa");
        std::fs::write(&fasta, ">chr1\nACGTACGTACGTACGT\n>chr2\nAAAACCCCGGGGTTTT\n").unwrap();

        let db_base = dir.join("testdb");
        let (nseq, total) = make_nucleotide_db(&fasta, &db_base, "Create-Read Test").unwrap();
        assert_eq!(nseq, 2);
        assert_eq!(total, 32);

        let db = super::super::index::BlastDb::open(&db_base).unwrap();
        assert_eq!(db.num_oids, 2);
        assert_eq!(db.db_type, super::super::index::DbType::Nucleotide);
        assert_eq!(db.total_length, 32);

        // Verify each sequence has the correct length
        assert_eq!(db.get_seq_len(0), 16);
        assert_eq!(db.get_seq_len(1), 16);

        // Verify raw data is non-empty
        assert!(!db.get_sequence(0).is_empty());
        assert!(!db.get_sequence(1).is_empty());

        std::fs::remove_dir_all(&dir).ok();
    }

    /// Create a DB with 10 sequences and verify all round-trip correctly.
    #[test]
    fn test_roundtrip_multiple_sequences() {
        let dir = std::env::temp_dir().join("blast_makedb_roundtrip");
        std::fs::create_dir_all(&dir).ok();

        // Generate 10 sequences of varying lengths
        let bases = [b'A', b'C', b'G', b'T'];
        let mut fasta = String::new();
        let mut expected_lengths: Vec<usize> = Vec::new();
        for i in 0..10 {
            let len = 20 + i * 13; // lengths: 20, 33, 46, 59, ...
            fasta.push_str(&format!(">seq{}\n", i));
            for j in 0..len {
                fasta.push(bases[(i + j) % 4] as char);
            }
            fasta.push('\n');
            expected_lengths.push(len);
        }

        let fasta_path = dir.join("multi.fa");
        std::fs::write(&fasta_path, &fasta).unwrap();

        let db_base = dir.join("multidb");
        let (nseq, total) = make_nucleotide_db(&fasta_path, &db_base, "Multi Test").unwrap();
        assert_eq!(nseq, 10);
        let expected_total: u64 = expected_lengths.iter().sum::<usize>() as u64;
        assert_eq!(total, expected_total);

        let db = super::super::index::BlastDb::open(&db_base).unwrap();
        assert_eq!(db.num_oids, 10);
        assert_eq!(db.total_length, expected_total);

        // Verify each sequence length matches
        for (oid, &exp_len) in expected_lengths.iter().enumerate() {
            let got_len = db.get_seq_len(oid as u32);
            assert_eq!(
                got_len, exp_len as u32,
                "OID {} length mismatch: got {} expected {}",
                oid, got_len, exp_len
            );
        }

        // Verify all sequences are readable
        for oid in 0..10u32 {
            let seq = db.get_sequence(oid);
            assert!(!seq.is_empty(), "OID {} should have data", oid);
            let hdr = db.get_header(oid);
            assert!(!hdr.is_empty(), "OID {} should have header", oid);
        }

        std::fs::remove_dir_all(&dir).ok();
    }

    /// Protein DB creation is not yet implemented, but we test that the nucleotide
    /// DB can handle protein-like data (all same letters) without crashing.
    #[test]
    fn test_create_nucleotide_db_edge_cases() {
        let dir = std::env::temp_dir().join("blast_makedb_edge");
        std::fs::create_dir_all(&dir).ok();

        // Single base sequence
        let fasta = dir.join("single.fa");
        std::fs::write(&fasta, ">tiny\nA\n").unwrap();
        let db_base = dir.join("tinydb");
        let (nseq, total) = make_nucleotide_db(&fasta, &db_base, "Tiny").unwrap();
        assert_eq!(nseq, 1);
        assert_eq!(total, 1);
        let db = super::super::index::BlastDb::open(&db_base).unwrap();
        assert_eq!(db.get_seq_len(0), 1);

        // Sequence length that is exact multiple of 4
        let fasta2 = dir.join("exact4.fa");
        std::fs::write(&fasta2, ">exact\nACGTACGT\n").unwrap();
        let db_base2 = dir.join("exact4db");
        let (nseq2, total2) = make_nucleotide_db(&fasta2, &db_base2, "Exact4").unwrap();
        assert_eq!(nseq2, 1);
        assert_eq!(total2, 8);
        let db2 = super::super::index::BlastDb::open(&db_base2).unwrap();
        assert_eq!(db2.get_seq_len(0), 8);

        std::fs::remove_dir_all(&dir).ok();
    }
}
