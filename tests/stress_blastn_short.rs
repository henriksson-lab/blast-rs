//! Stress tests for blastn-short with short primer queries against larger databases.
//! These tests verify that the search doesn't stack-overflow with many sequences
//! and multithreaded execution.

use blast_rs::db::DbType;
use blast_rs::{blastn, BlastDbBuilder, SearchParams, SequenceEntry};
use std::path::Path;
use tempfile::TempDir;

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
    for e in entries {
        builder.add(e);
    }
    builder.write(&base).unwrap();
    let db = blast_rs::db::BlastDb::open(&base).unwrap();
    (tmp, db)
}

/// Generate a pseudo-random nucleotide sequence.
fn random_nt_seq(len: usize, seed: u64) -> Vec<u8> {
    let bases = b"ACGT";
    let mut state = seed;
    (0..len)
        .map(|_| {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            bases[((state >> 33) % 4) as usize]
        })
        .collect()
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
