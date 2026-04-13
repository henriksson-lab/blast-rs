//! Stress tests for blastn-short with short primer queries against larger databases.
//! These tests verify that the search doesn't stack-overflow with many sequences
//! and multithreaded execution.

use tempfile::TempDir;
use blast_rs::db::DbType;
use blast_rs::{BlastDbBuilder, SequenceEntry, SearchParams, blastn};
use std::path::Path;

fn nt_entry(acc: &str, title: &str, seq: &str) -> SequenceEntry {
    SequenceEntry {
        title: title.to_string(),
        accession: acc.to_string(),
        sequence: seq.as_bytes().to_vec(),
        taxid: None,
    }
}

fn build_nucleotide_db(entries: Vec<SequenceEntry>) -> (TempDir, blast_rs::db::BlastDb) {
    let tmp = TempDir::new().unwrap();
    let base = tmp.path().join("testdb");
    let mut builder = BlastDbBuilder::new(DbType::Nucleotide, "stress test db");
    for e in entries { builder.add(e); }
    builder.write(&base).unwrap();
    let db = blast_rs::db::BlastDb::open(&base).unwrap();
    (tmp, db)
}

/// Generate a pseudo-random nucleotide sequence.
fn random_nt_seq(len: usize, seed: u64) -> Vec<u8> {
    let bases = b"ACGT";
    let mut state = seed;
    (0..len).map(|_| {
        state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        bases[((state >> 33) % 4) as usize]
    }).collect()
}

fn blastn_short_params() -> SearchParams {
    let mut params = SearchParams::blastn();
    params.word_size = 7;
    params.match_score = 1;
    params.mismatch = -3;
    params.gap_open = 5;
    params.gap_extend = 2;
    params.evalue_threshold = 5.0;
    params.max_target_seqs = 10000;
    params.filter_low_complexity = false;
    params
}

/// Stress test: 2000 subjects, 8 threads, short primer.
/// This simulates the conditions that caused the original stack overflow.
#[test]
fn stress_blastn_short_2000_subjects_8_threads() {
    let primer = b"ATCGATCGATCGATCGATCGATCGATCG"; // 28bp primer

    let mut entries = Vec::new();
    for i in 0..2000u64 {
        let mut seq = random_nt_seq(5000, i * 31337 + 1);
        // Embed primer in ~20% of subjects
        if i % 5 == 0 {
            let pos = (i as usize * 47) % (seq.len() - primer.len());
            seq[pos..pos + primer.len()].copy_from_slice(primer);
        }
        entries.push(nt_entry(
            &format!("stress_{}", i),
            &format!("stress test sequence {}", i),
            &String::from_utf8(seq).unwrap(),
        ));
    }

    let (_tmp, db) = build_nucleotide_db(entries);

    let mut params = blastn_short_params();
    params.num_threads = 8;

    let results = blastn(&db, primer, &params);
    assert!(!results.is_empty(), "Should find hits in stress test");
    eprintln!("Stress test found {} subject hits", results.len());
}

/// Stress test: many subjects with long sequences.
#[test]
fn stress_blastn_short_long_subjects() {
    let primer = b"GCTAGCTAGCTAGCTAGCTAGCTA"; // 24bp

    let mut entries = Vec::new();
    for i in 0..100u64 {
        let mut seq = random_nt_seq(50_000, i * 42 + 99);
        // Embed primer
        let pos = (i as usize * 123) % (seq.len() - primer.len());
        seq[pos..pos + primer.len()].copy_from_slice(primer);
        entries.push(nt_entry(
            &format!("longseq_{}", i),
            &format!("long sequence {}", i),
            &String::from_utf8(seq).unwrap(),
        ));
    }

    let (_tmp, db) = build_nucleotide_db(entries);

    let mut params = blastn_short_params();
    params.num_threads = 4;

    let results = blastn(&db, primer, &params);
    assert!(!results.is_empty(), "Should find hits in long-subject stress test");
    eprintln!("Long-subject stress test found {} subject hits", results.len());
}

/// Stress test: very short primer (15bp) which is common for PCR primers.
#[test]
fn stress_blastn_short_15bp_primer() {
    let primer = b"ATCGATCGATCGATC"; // 15bp

    let mut entries = Vec::new();
    for i in 0..500u64 {
        let mut seq = random_nt_seq(10_000, i * 9999);
        let pos = (i as usize * 67) % (seq.len() - primer.len());
        seq[pos..pos + primer.len()].copy_from_slice(primer);
        entries.push(nt_entry(
            &format!("s15_{}", i),
            &format!("15bp test seq {}", i),
            &String::from_utf8(seq).unwrap(),
        ));
    }

    let (_tmp, db) = build_nucleotide_db(entries);

    let mut params = blastn_short_params();
    params.num_threads = 8;
    params.max_hsps = Some(1);

    let results = blastn(&db, primer, &params);
    assert!(!results.is_empty(), "Should find 15bp primer hits");
    eprintln!("15bp primer stress test found {} subject hits", results.len());
}

/// Test with real V5 NCBI database (16S_ribosomal_RNA) if available.
/// This is the closest we can get to the user's reported core_nt scenario.
#[test]
fn stress_real_v5_database_if_available() {
    let db_path = Path::new("/husky/henriksson/for_claude/blast/16S_ribosomal_RNA");
    if !db_path.with_extension("nin").exists() {
        eprintln!("Skipping: 16S_ribosomal_RNA database not found at {:?}", db_path);
        return;
    }

    let db = blast_rs::db::BlastDb::open(db_path).expect("Failed to open 16S database");
    assert_eq!(db.db_type, DbType::Nucleotide);
    assert!(db.num_oids > 1000, "16S DB should have thousands of sequences, got {}", db.num_oids);
    eprintln!("16S DB: {} OIDs, version {}", db.num_oids, db.version);

    // 19bp primer that matches many 16S sequences
    let primer = b"GTGCCAGCAGCCGCGGTAA";

    let mut params = blastn_short_params();
    params.num_threads = 8;
    params.max_hsps = Some(1);

    let results = blastn(&db, primer, &params);
    assert!(!results.is_empty(), "Should find primer hits in 16S database");
    eprintln!("Real V5 DB stress test found {} subject hits", results.len());

    // Verify taxonomy data is populated
    let with_taxids: usize = results.iter().filter(|r| !r.taxids.is_empty()).count();
    eprintln!("  {} hits with taxonomy IDs", with_taxids);
    assert!(with_taxids > 0, "V5 database should have taxonomy IDs");
}
