//! makeblastdb equivalent — create BLAST databases from FASTA files.

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

use crate::db::defline::encode_defline_asn1;
use crate::db::index_writer::write_index_file;
use crate::encoding::{
    encode_ncbi2na_ambiguity_data, encode_ncbi2na_sequence, encode_ncbistdaa_sequence,
};

fn db_component_path(base_path: &Path, ext: &str) -> std::path::PathBuf {
    let mut path = base_path.as_os_str().to_os_string();
    path.push(".");
    path.push(ext);
    path.into()
}

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
    let mut nsq = BufWriter::new(File::create(db_component_path(output_base, "nsq"))?);
    nsq.write_all(&[0u8])?; // sentinel byte

    let mut seq_offsets = vec![1u32]; // first seq starts at byte 1
    let mut amb_offsets = Vec::new();

    for (_, seq) in &sequences {
        let seq_start = seq_offsets.last().copied().unwrap_or(1);
        let packed = encode_ncbi2na_sequence(seq);

        nsq.write_all(&packed)?;
        let amb_offset = seq_start as u32 + packed.len() as u32;
        amb_offsets.push(amb_offset);
        let ambiguity_data = encode_ncbi2na_ambiguity_data(seq);
        nsq.write_all(&ambiguity_data)?;
        seq_offsets.push(amb_offset + ambiguity_data.len() as u32);
        total_length += seq.len() as u64;
    }
    amb_offsets.push(*seq_offsets.last().unwrap_or(&0));
    nsq.flush()?;

    // Write .nhr (header data) — ASN.1 BER encoded Blast-def-line-set
    let mut nhr = BufWriter::new(File::create(db_component_path(output_base, "nhr"))?);
    let mut hdr_offsets = vec![0u32];
    for (oid, (header, _)) in sequences.iter().enumerate() {
        let start = hdr_offsets.last().copied().unwrap_or(0);
        let asn1 = encode_defline_asn1(header, oid as i32);
        nhr.write_all(&asn1)?;
        hdr_offsets.push(start + asn1.len() as u32);
    }
    nhr.flush()?;

    let max_len = sequences
        .iter()
        .map(|(_, s)| s.len() as u32)
        .max()
        .unwrap_or(0);
    write_index_file(
        &db_component_path(output_base, "nin"),
        4,
        crate::db::DbType::Nucleotide,
        title,
        num_seqs,
        total_length,
        max_len,
        &hdr_offsets,
        &seq_offsets,
        Some(&amb_offsets),
    )?;

    Ok((num_seqs, total_length))
}

/// Create a BLAST v4 protein database from a FASTA file.
/// blast-rs: Standalone Rust makeblastdb-style protein writer; not a direct
/// NCBI C port.
pub fn make_protein_db(
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
                if b.is_ascii_alphabetic() || b == b'*' {
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

    // Write .psq (protein sequence data in NCBIstdaa)
    let mut psq = BufWriter::new(File::create(db_component_path(output_base, "psq"))?);
    psq.write_all(&[0u8])?; // sentinel byte

    let mut seq_offsets = vec![1u32]; // first seq starts at byte 1
    for (_, seq) in &sequences {
        let encoded = encode_ncbistdaa_sequence(seq);
        psq.write_all(&encoded)?;
        psq.write_all(&[0u8])?; // sentinel between sequences
        seq_offsets.push(seq_offsets.last().copied().unwrap_or(1) + encoded.len() as u32 + 1);
        total_length += seq.len() as u64;
    }
    psq.flush()?;

    // Write .phr (header data) — ASN.1 BER encoded Blast-def-line-set
    let mut phr = BufWriter::new(File::create(db_component_path(output_base, "phr"))?);
    let mut hdr_offsets = vec![0u32];
    for (oid, (header, _)) in sequences.iter().enumerate() {
        let start = hdr_offsets.last().copied().unwrap_or(0);
        let asn1 = encode_defline_asn1(header, oid as i32);
        phr.write_all(&asn1)?;
        hdr_offsets.push(start + asn1.len() as u32);
    }
    phr.flush()?;

    let max_len = sequences
        .iter()
        .map(|(_, s)| s.len() as u32)
        .max()
        .unwrap_or(0);
    write_index_file(
        &db_component_path(output_base, "pin"),
        4,
        crate::db::DbType::Protein,
        title,
        num_seqs,
        total_length,
        max_len,
        &hdr_offsets,
        &seq_offsets,
        None,
    )?;

    Ok((num_seqs, total_length))
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
        assert!(db_component_path(&db_base, "nin").exists());
        assert!(db_component_path(&db_base, "nsq").exists());
        assert!(db_component_path(&db_base, "nhr").exists());

        // Try opening with our reader
        let db = super::super::index::BlastDb::open(&db_base).unwrap();
        assert_eq!(db.num_oids, 2);
        assert_eq!(db.total_length, 16);
        assert!(is_blastdb_date(&db.date));
        assert_eq!(db.get_accession(0).as_deref(), Some("seq1"));
        assert_eq!(db.get_defline(0).as_deref(), Some("seq1"));
        let raw_header = std::fs::read(db_component_path(&db_base, "nhr")).unwrap();
        assert!(raw_header
            .windows(b"BL_ORD_ID".len())
            .any(|w| w == b"BL_ORD_ID"));

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn test_make_nucleotide_db_appends_component_extensions() {
        let dir = std::env::temp_dir().join("blast_makedb_append_ext_test");
        std::fs::create_dir_all(&dir).ok();

        let fasta = dir.join("input.fa");
        std::fs::write(&fasta, ">seq1\nACGT\n").unwrap();

        let db_base = dir.join("testdb.00");
        make_nucleotide_db(&fasta, &db_base, "Append Ext Test").unwrap();
        assert!(db_component_path(&db_base, "nin").exists());
        assert!(db_component_path(&db_base, "nsq").exists());
        assert!(db_component_path(&db_base, "nhr").exists());
        assert!(!db_base.with_extension("nin").exists());

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
        assert_eq!(db.get_accession(0).as_deref(), Some("chr1"));
        assert_eq!(db.get_accession(1).as_deref(), Some("chr2"));

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

    /// Nucleotide DB creation handles small edge-case sequences without crashing.
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

    #[test]
    fn test_make_protein_db_roundtrips_with_reader() {
        let dir = std::env::temp_dir().join("blast_makedb_protein");
        std::fs::create_dir_all(&dir).ok();

        let fasta = dir.join("protein.fa");
        std::fs::write(&fasta, ">prot1 desc\nMTEYK\n>prot2\nACDEFGHIK\n").unwrap();

        let db_base = dir.join("protdb");
        let (nseq, total) = make_protein_db(&fasta, &db_base, "Protein Test").unwrap();
        assert_eq!(nseq, 2);
        assert_eq!(total, 14);

        assert!(db_component_path(&db_base, "pin").exists());
        assert!(db_component_path(&db_base, "psq").exists());
        assert!(db_component_path(&db_base, "phr").exists());

        let db = super::super::index::BlastDb::open(&db_base).unwrap();
        assert_eq!(db.num_oids, 2);
        assert_eq!(db.db_type, super::super::index::DbType::Protein);
        assert_eq!(db.total_length, 14);
        assert_eq!(db.get_seq_len(0), 5);
        assert_eq!(db.get_seq_len(1), 9);
        assert_eq!(db.get_accession(0).as_deref(), Some("prot1"));
        assert_eq!(db.get_defline(0).as_deref(), Some("prot1 desc"));
        assert_eq!(
            db.get_sequence(0),
            crate::encoding::encode_ncbistdaa_sequence(b"MTEYK")
        );
        assert_eq!(
            db.get_sequence(1),
            crate::encoding::encode_ncbistdaa_sequence(b"ACDEFGHIK")
        );

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn test_make_nucleotide_db_writes_ambiguity_data() {
        let dir = std::env::temp_dir().join("blast_makedb_ambiguity");
        std::fs::create_dir_all(&dir).ok();

        let fasta = dir.join("amb.fa");
        std::fs::write(&fasta, ">amb1\nACGTRRNNY\n").unwrap();

        let db_base = dir.join("ambdb");
        let (nseq, total) = make_nucleotide_db(&fasta, &db_base, "Ambiguity Test").unwrap();
        assert_eq!(nseq, 1);
        assert_eq!(total, 9);

        let db = super::super::index::BlastDb::open(&db_base).unwrap();
        let ambiguity = db
            .get_ambiguity_data(0)
            .expect("ambiguous FASTA bases should be encoded");
        assert!(!ambiguity.is_empty());

        let decoded = crate::search::decode_packed_ncbi2na_with_ambiguity(
            db.get_sequence(0),
            db.get_seq_len(0) as usize,
            ambiguity,
        );
        assert_eq!(
            decoded,
            vec![0, 1, 2, 3, 4, 4, 14, 14, 5],
            "decoded BLASTNA sequence should preserve R/N/Y ambiguity"
        );

        std::fs::remove_dir_all(&dir).ok();
    }

    fn is_blastdb_date(date: &str) -> bool {
        date.len() == 10
            && date.as_bytes()[4] == b'-'
            && date.as_bytes()[7] == b'-'
            && date
                .bytes()
                .enumerate()
                .all(|(i, b)| matches!(i, 4 | 7) || b.is_ascii_digit())
    }
}
