//! Comparison tests: run both Rust and reference BLAST, diff output.

use std::path::{Path, PathBuf};
use std::process::Command;

fn project_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .to_path_buf()
}

fn rust_blastn() -> PathBuf {
    project_root().join("target/debug/blast-cli")
}

fn ref_blastn() -> PathBuf {
    project_root().join("ncbi-blast-2.17.0+-src/c++/ReleaseMT/bin/blastn")
}

fn test_db(name: &str) -> PathBuf {
    project_root().join(format!(
        "ncbi-blast-2.17.0+-src/c++/src/algo/blast/unit_tests/api/data/{}",
        name
    ))
}

fn run_rust(query: &Path, db: &Path, word_size: i32) -> String {
    let output = Command::new(rust_blastn())
        .args([
            "--query", query.to_str().unwrap(),
            "--db", db.to_str().unwrap(),
            "--word_size", &word_size.to_string(),
        ])
        .output()
        .expect("Failed to run Rust blastn");
    String::from_utf8_lossy(&output.stdout).to_string()
}

fn run_ref(query: &Path, db: &Path, word_size: i32) -> String {
    let output = Command::new(ref_blastn())
        .args([
            "-query", query.to_str().unwrap(),
            "-db", db.to_str().unwrap(),
            "-outfmt", "6",
            "-word_size", &word_size.to_string(),
            "-reward", "1",
            "-penalty", "-3",
            "-gapopen", "5",
            "-gapextend", "2",
        ])
        .output()
        .expect("Failed to run reference blastn");
    String::from_utf8_lossy(&output.stdout).to_string()
}

#[test]
fn test_blastn_seqn_word11() {
    let query = project_root().join("tests/fixtures/test_query.fa");
    let db = test_db("seqn");

    if !ref_blastn().exists() {
        eprintln!("Skipping: reference blastn not built");
        return;
    }

    let rust_out = run_rust(&query, &db, 11);
    let ref_out = run_ref(&query, &db, 11);

    assert_eq!(
        rust_out, ref_out,
        "\n=== RUST ===\n{}\n=== REF ===\n{}",
        rust_out, ref_out
    );
}

#[test]
fn test_blastn_seqn_word28_megablast() {
    let query = project_root().join("tests/fixtures/test_query.fa");
    let db = test_db("seqn");

    if !ref_blastn().exists() {
        eprintln!("Skipping: reference blastn not built");
        return;
    }

    let rust_out = run_rust(&query, &db, 28);
    let ref_out = run_ref(&query, &db, 28);

    assert_eq!(
        rust_out, ref_out,
        "\n=== RUST ===\n{}\n=== REF ===\n{}",
        rust_out, ref_out
    );
}

#[test]
fn test_blastn_no_hits() {
    // Random query that shouldn't match
    let query = project_root().join("tests/fixtures/test_query_random.fa");
    let db = test_db("seqn");

    if !query.exists() || !ref_blastn().exists() {
        eprintln!("Skipping: missing fixtures or reference");
        return;
    }

    let rust_out = run_rust(&query, &db, 28);
    let ref_out = run_ref(&query, &db, 28);

    assert_eq!(
        rust_out, ref_out,
        "No-hit search should produce identical (empty) output"
    );
}
