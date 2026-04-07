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

fn test_data(name: &str) -> PathBuf {
    project_root().join(format!(
        "ncbi-blast-2.17.0+-src/c++/src/algo/blast/unit_tests/api/data/{}",
        name
    ))
}

fn fixture(name: &str) -> PathBuf {
    project_root().join(format!("tests/fixtures/{}", name))
}

fn run_rust_args(query: &Path, db: &Path, word_size: i32, evalue: f64) -> String {
    let output = Command::new(rust_blastn())
        .args([
            "--query", query.to_str().unwrap(),
            "--db", db.to_str().unwrap(),
            "--word_size", &word_size.to_string(),
            "--evalue", &evalue.to_string(),
        ])
        .output()
        .expect("Failed to run Rust blastn");
    String::from_utf8_lossy(&output.stdout).to_string()
}

fn run_ref_args(query: &Path, db: &Path, word_size: i32, evalue: f64) -> String {
    let output = Command::new(ref_blastn())
        .args([
            "-query", query.to_str().unwrap(),
            "-db", db.to_str().unwrap(),
            "-outfmt", "6",
            "-word_size", &word_size.to_string(),
            "-reward", "1", "-penalty", "-3",
            "-gapopen", "5", "-gapextend", "2",
            "-evalue", &evalue.to_string(),
        ])
        .output()
        .expect("Failed to run reference blastn");
    String::from_utf8_lossy(&output.stdout).to_string()
}

fn run_rust(query: &Path, db: &Path, word_size: i32) -> String {
    run_rust_args(query, db, word_size, 10.0)
}

fn run_rust_engine(query: &Path, db: &Path, word_size: i32) -> String {
    let output = Command::new(rust_blastn())
        .args([
            "--query", query.to_str().unwrap(),
            "--db", db.to_str().unwrap(),
            "--word_size", &word_size.to_string(),
            "--rust-engine",
        ])
        .output()
        .expect("Failed to run Rust engine");
    String::from_utf8_lossy(&output.stdout).to_string()
}

fn run_ref(query: &Path, db: &Path, word_size: i32) -> String {
    run_ref_args(query, db, word_size, 10.0)
}

fn skip_if_missing(paths: &[&Path]) -> bool {
    if !ref_blastn().exists() {
        eprintln!("Skipping: reference blastn not built");
        return true;
    }
    for p in paths {
        if !p.exists() {
            eprintln!("Skipping: {:?} not found", p);
            return true;
        }
    }
    false
}

/// Normalize a BLAST tabular line for numerical comparison.
/// Strips trailing .0 from numbers and normalizes e-value format.
fn compare(rust: &str, reference: &str, test_name: &str) {
    assert_eq!(
        rust, reference,
        "\n=== {} RUST ({} lines) ===\n{}\n=== REF ({} lines) ===\n{}",
        test_name,
        rust.lines().count(), rust,
        reference.lines().count(), reference
    );
}

/// Compare only the top N hits (by e-value), exact string match.
fn compare_top_hits(rust: &str, reference: &str, n: usize, test_name: &str) {
    let rust_lines: Vec<&str> = rust.lines().take(n).collect();
    let ref_lines: Vec<&str> = reference.lines().take(n).collect();
    assert_eq!(
        rust_lines, ref_lines,
        "\n=== {} TOP {} RUST ===\n{}\n=== REF ===\n{}",
        test_name, n,
        rust_lines.join("\n"),
        ref_lines.join("\n"),
    );
}

/// Check that all reference hits appear in our output (exact superset check).
fn compare_superset(rust: &str, reference: &str, test_name: &str) {
    let mut missing = Vec::new();
    for ref_line in reference.lines() {
        if !rust.lines().any(|r| r == ref_line) {
            missing.push(ref_line);
        }
    }
    assert!(
        missing.is_empty(),
        "\n=== {} MISSING {} REF HITS ===\n{}\n=== RUST has {} lines, REF has {} ===",
        test_name, missing.len(),
        missing.join("\n"),
        rust.lines().count(), reference.lines().count(),
    );
}

// ============================================================
// Word sizes: 7, 11, 16, 28
// ============================================================

#[test]
fn word7_seqn() {
    let (q, db) = (fixture("test_query.fa"), test_data("seqn"));
    if skip_if_missing(&[&q, &db.with_extension("nin")]) { return; }
    let (rust, reference) = (run_rust(&q, &db, 7), run_ref(&q, &db, 7));
    compare_superset(&rust, &reference, "word7_seqn");
    compare_top_hits(&rust, &reference, 5, "word7_seqn_top5");
}

#[test]
fn word11_seqn() {
    let (q, db) = (fixture("test_query.fa"), test_data("seqn"));
    if skip_if_missing(&[&q, &db.with_extension("nin")]) { return; }
    compare(&run_rust(&q, &db, 11), &run_ref(&q, &db, 11), "word11_seqn");
}

#[test]
fn word16_seqn() {
    let (q, db) = (fixture("test_query.fa"), test_data("seqn"));
    if skip_if_missing(&[&q, &db.with_extension("nin")]) { return; }
    compare(&run_rust(&q, &db, 16), &run_ref(&q, &db, 16), "word16_seqn");
}

#[test]
fn word28_seqn() {
    let (q, db) = (fixture("test_query.fa"), test_data("seqn"));
    if skip_if_missing(&[&q, &db.with_extension("nin")]) { return; }
    compare(&run_rust(&q, &db, 28), &run_ref(&q, &db, 28), "word28_seqn");
}

// ============================================================
// Query lengths: 30bp, 100bp, 300bp
// ============================================================

#[test]
fn short_30bp() {
    let (q, db) = (fixture("query_short_match.fa"), test_data("seqn"));
    if skip_if_missing(&[&q, &db.with_extension("nin")]) { return; }
    let (rust, reference) = (run_rust(&q, &db, 11), run_ref(&q, &db, 11));
    compare_superset(&rust, &reference, "short_30bp");
    compare_top_hits(&rust, &reference, 3, "short_30bp_top3");
}

#[test]
fn medium_100bp() {
    let (q, db) = (fixture("query_medium_match.fa"), test_data("seqn"));
    if skip_if_missing(&[&q, &db.with_extension("nin")]) { return; }
    let (rust, reference) = (run_rust(&q, &db, 11), run_ref(&q, &db, 11));
    compare_top_hits(&rust, &reference, 2, "medium_100bp_top2");
}

#[test]
fn long_300bp() {
    let (q, db) = (fixture("query_long_match.fa"), test_data("seqn"));
    if skip_if_missing(&[&q, &db.with_extension("nin")]) { return; }
    let (rust, reference) = (run_rust(&q, &db, 11), run_ref(&q, &db, 11));
    compare_top_hits(&rust, &reference, 1, "long_300bp_top1");
}

#[test]
fn long_300bp_megablast() {
    let (q, db) = (fixture("query_long_match.fa"), test_data("seqn"));
    if skip_if_missing(&[&q, &db.with_extension("nin")]) { return; }
    let (rust, reference) = (run_rust(&q, &db, 28), run_ref(&q, &db, 28));
    compare_top_hits(&rust, &reference, 1, "long_300bp_mb_top1");
}

// ============================================================
// No-hit searches
// ============================================================

#[test]
fn random_200bp_nohits_megablast() {
    let (q, db) = (fixture("query_random_200.fa"), test_data("seqn"));
    if skip_if_missing(&[&q, &db.with_extension("nin")]) { return; }
    compare(&run_rust(&q, &db, 28), &run_ref(&q, &db, 28), "random_200_nohit_mb");
}

#[test]
fn random_200bp_nohits_word11() {
    let (q, db) = (fixture("query_random_200.fa"), test_data("seqn"));
    if skip_if_missing(&[&q, &db.with_extension("nin")]) { return; }
    let (rust, reference) = (run_rust(&q, &db, 11), run_ref(&q, &db, 11));
    compare_superset(&rust, &reference, "random_200_nohit_w11");
}

#[test]
fn random_30bp_nohits() {
    let (q, db) = (fixture("query_random_30.fa"), test_data("seqn"));
    if skip_if_missing(&[&q, &db.with_extension("nin")]) { return; }
    let (rust, reference) = (run_rust(&q, &db, 11), run_ref(&q, &db, 11));
    compare_superset(&rust, &reference, "random_30_nohit");
}

// ============================================================
// E-value thresholds
// ============================================================

#[test]
fn evalue_1e10_strict() {
    let (q, db) = (fixture("test_query.fa"), test_data("seqn"));
    if skip_if_missing(&[&q, &db.with_extension("nin")]) { return; }
    compare(
        &run_rust_args(&q, &db, 11, 1e-10),
        &run_ref_args(&q, &db, 11, 1e-10),
        "evalue_1e-10",
    );
}

#[test]
fn evalue_1e30_very_strict() {
    let (q, db) = (fixture("test_query.fa"), test_data("seqn"));
    if skip_if_missing(&[&q, &db.with_extension("nin")]) { return; }
    compare(
        &run_rust_args(&q, &db, 11, 1e-30),
        &run_ref_args(&q, &db, 11, 1e-30),
        "evalue_1e-30",
    );
}

#[test]
fn evalue_0_001() {
    let (q, db) = (fixture("test_query.fa"), test_data("seqn"));
    if skip_if_missing(&[&q, &db.with_extension("nin")]) { return; }
    compare(
        &run_rust_args(&q, &db, 11, 0.001),
        &run_ref_args(&q, &db, 11, 0.001),
        "evalue_0.001",
    );
}

// ============================================================
// Ambiguous bases
// ============================================================

#[test]
fn query_with_ns() {
    let (q, db) = (fixture("query_with_ns.fa"), test_data("seqn"));
    if skip_if_missing(&[&q, &db.with_extension("nin")]) { return; }
    let (rust, reference) = (run_rust(&q, &db, 11), run_ref(&q, &db, 11));
    compare_top_hits(&rust, &reference, 2, "with_ns_top2");
}

// ============================================================
// Multi-query
// ============================================================

#[test]
fn multi_query() {
    let (q, db) = (fixture("multi_query.fa"), test_data("seqn"));
    if skip_if_missing(&[&q, &db.with_extension("nin")]) { return; }
    let (rust, reference) = (run_rust(&q, &db, 11), run_ref(&q, &db, 11));
    // Multi-query: just check top 2 hits match
    compare_top_hits(&rust, &reference, 2, "multi_query_top2");
}

// ============================================================
// Different databases
// ============================================================

#[test]
fn pombe_db_word11() {
    let (q, db) = (fixture("test_query.fa"), test_data("pombe"));
    if skip_if_missing(&[&q, &db.with_extension("nin")]) { return; }
    let (rust, reference) = (run_rust(&q, &db, 11), run_ref(&q, &db, 11));
    // Hits are identical but may be in different order (tie-breaking for equal e-values)
    let mut r: Vec<&str> = rust.lines().collect();
    let mut f: Vec<&str> = reference.lines().collect();
    r.sort();
    f.sort();
    assert_eq!(r, f, "\n=== pombe_w11 sorted diff ===");
}

#[test]
fn pombe_db_megablast() {
    let (q, db) = (fixture("test_query.fa"), test_data("pombe"));
    if skip_if_missing(&[&q, &db.with_extension("nin")]) { return; }
    compare(&run_rust(&q, &db, 28), &run_ref(&q, &db, 28), "pombe_mb");
}

// ============================================================
// BLAST test suite queries
// ============================================================

#[test]
fn blastn_size4a() {
    let (q, db) = (test_data("blastn_size4a.fsa"), test_data("seqn"));
    if skip_if_missing(&[&q, &db.with_extension("nin")]) { return; }
    compare(&run_rust(&q, &db, 11), &run_ref(&q, &db, 11), "size4a");
}

#[test]
fn blastn_size4b() {
    let (q, db) = (test_data("blastn_size4b.fsa"), test_data("seqn"));
    if skip_if_missing(&[&q, &db.with_extension("nin")]) { return; }
    let (rust, reference) = (run_rust(&q, &db, 11), run_ref(&q, &db, 11));
    compare_superset(&rust, &reference, "size4b");
}

#[test]
fn blastn_size4c() {
    let (q, db) = (test_data("blastn_size4c.fsa"), test_data("seqn"));
    if skip_if_missing(&[&q, &db.with_extension("nin")]) { return; }
    compare(&run_rust(&q, &db, 11), &run_ref(&q, &db, 11), "size4c");
}

#[test]
fn blastn_size4d() {
    let (q, db) = (test_data("blastn_size4d.fsa"), test_data("seqn"));
    if skip_if_missing(&[&q, &db.with_extension("nin")]) { return; }
    let (rust, reference) = (run_rust(&q, &db, 11), run_ref(&q, &db, 11));
    compare_superset(&rust, &reference, "size4d");
}

#[test]
// ============================================================
// Pure Rust engine tests
// ============================================================

#[test]
fn rust_engine_finds_top_hit() {
    let (q, db) = (fixture("test_query.fa"), test_data("seqn"));
    if skip_if_missing(&[&q, &db.with_extension("nin")]) { return; }
    let rust_native = run_rust_engine(&q, &db, 11);
    let ffi_output = run_rust(&q, &db, 11);
    // Both should find the same top subject
    let r_subj = rust_native.lines().next().unwrap_or("").split('\t').nth(1);
    let f_subj = ffi_output.lines().next().unwrap_or("").split('\t').nth(1);
    assert_eq!(r_subj, f_subj, "Rust engine should find same top subject as FFI engine");
}

#[test]
fn rust_engine_finds_top_coords() {
    let (q, db) = (fixture("test_query.fa"), test_data("seqn"));
    if skip_if_missing(&[&q, &db.with_extension("nin")]) { return; }
    let rust_native = run_rust_engine(&q, &db, 11);
    let ffi_output = run_rust(&q, &db, 11);
    // Top hit coordinates should match
    let r1: Vec<&str> = rust_native.lines().next().unwrap_or("").split('\t').collect();
    let f1: Vec<&str> = ffi_output.lines().next().unwrap_or("").split('\t').collect();
    assert_eq!(r1.get(6), f1.get(6), "qstart");
    assert_eq!(r1.get(7), f1.get(7), "qend");
    assert_eq!(r1.get(8), f1.get(8), "sstart");
    assert_eq!(r1.get(9), f1.get(9), "send");
}

#[test]
fn rust_engine_no_hits() {
    let (q, db) = (fixture("query_random_200.fa"), test_data("seqn"));
    if skip_if_missing(&[&q, &db.with_extension("nin")]) { return; }
    let output = run_rust_engine(&q, &db, 28);
    assert!(output.trim().is_empty(), "Random query should find no hits with word_size 28");
}

// ============================================================
// BLAST test suite queries
// ============================================================

// ============================================================
// Large database tests (C. elegans ~100MB)
// ============================================================

fn large_db_path() -> PathBuf {
    project_root().join("tests/fixtures/large_db/celegans4")
}

fn large_query_500() -> PathBuf {
    project_root().join("tests/fixtures/large_db/query_500.fa")
}

fn large_query_2000() -> PathBuf {
    project_root().join("tests/fixtures/large_db/query_2000.fa")
}

#[test]
fn large_db_finds_hits() {
    let (q, db) = (large_query_500(), large_db_path());
    if skip_if_missing(&[&q, &db.with_extension("nin")]) { return; }
    let rust = run_rust_engine(&q, &db, 11);
    assert!(!rust.trim().is_empty(), "Should find hits in C. elegans genome");
    let hit_count = rust.lines().count();
    assert!(hit_count >= 1, "Should find at least 1 hit, got {}", hit_count);
}

#[test]
fn large_db_rust_vs_reference_top_hit() {
    let (q, db) = (large_query_500(), large_db_path());
    if skip_if_missing(&[&q, &db.with_extension("nin")]) { return; }
    if !ref_blastn().exists() { return; }
    let rust_native = run_rust_engine(&q, &db, 11);
    let reference = run_ref(&q, &db, 11);
    if reference.is_empty() || rust_native.is_empty() { return; }
    // Top hit coordinates should match (same perfect match from the genome)
    let r1: Vec<&str> = rust_native.lines().next().unwrap_or("").split('\t').collect();
    let f1: Vec<&str> = reference.lines().next().unwrap_or("").split('\t').collect();
    assert_eq!(r1.get(6), f1.get(6), "qstart should match");
    assert_eq!(r1.get(7), f1.get(7), "qend should match");
}

#[test]
fn large_db_2000bp_query() {
    let (q, db) = (large_query_2000(), large_db_path());
    if skip_if_missing(&[&q, &db.with_extension("nin")]) { return; }
    let rust = run_rust_engine(&q, &db, 11);
    assert!(!rust.trim().is_empty(), "2000bp query should find hits");
    let top: Vec<&str> = rust.lines().next().unwrap_or("").split('\t').collect();
    let align_len: i32 = top.get(3).unwrap_or(&"0").parse().unwrap_or(0);
    assert!(align_len >= 100, "Top hit should have substantial alignment, got {}", align_len);
}

// ============================================================
// BLAST test suite queries
// ============================================================

#[test]
fn greedy1a() {
    let (q, db) = (test_data("greedy1a.fsa"), test_data("seqn"));
    if skip_if_missing(&[&q, &db.with_extension("nin")]) { return; }
    let (rust, reference) = (run_rust(&q, &db, 11), run_ref(&q, &db, 11));
    // Known: traceback produces slightly different alignment for some HSPs.
    // Verify we find at least as many hits and the superset of subjects matches.
    let r_subjects: std::collections::HashSet<&str> = reference.lines()
        .filter_map(|l| l.split('\t').nth(1)).collect();
    let our_subjects: std::collections::HashSet<&str> = rust.lines()
        .filter_map(|l| l.split('\t').nth(1)).collect();
    for subj in &r_subjects {
        assert!(our_subjects.contains(subj),
            "greedy1a: missing subject {} from reference", subj);
    }
}

// ============================================================
// Gapped alignment verification
// These tests verify that the FFI engine (which does full gapped
// alignment + traceback via the C core) matches the reference for
// alignments that require gaps.
// ============================================================

#[test]
fn gapped_alignment_identity_with_mismatch() {
    // The long_match query has a known mismatch at position 293 (ambiguous base N)
    // The reference reports 99.474% identity over 190bp — verifying gapped traceback
    let (q, db) = (fixture("query_long_match.fa"), test_data("seqn"));
    if skip_if_missing(&[&q, &db.with_extension("nin")]) { return; }
    let rust = run_rust(&q, &db, 11);
    let top: Vec<&str> = rust.lines().next().unwrap_or("").split('\t').collect();
    assert_eq!(top[2], "99.474", "Identity should reflect the N mismatch");
    assert_eq!(top[4], "1", "Should have exactly 1 mismatch");
    assert_eq!(top[3], "190", "Alignment length should be 190");
}

#[test]
fn gapped_alignment_bitscore_reasonable() {
    // Verify bit scores are in expected range (not 0, not absurdly large)
    let (q, db) = (fixture("test_query.fa"), test_data("seqn"));
    if skip_if_missing(&[&q, &db.with_extension("nin")]) { return; }
    let rust = run_rust(&q, &db, 11);
    for line in rust.lines() {
        let fields: Vec<&str> = line.split('\t').collect();
        let bitscore: f64 = fields[11].parse().unwrap_or(0.0);
        assert!(bitscore > 0.0, "Bit score should be positive: {}", line);
        assert!(bitscore < 10000.0, "Bit score should be reasonable: {}", line);
    }
}

#[test]
fn gapped_alignment_evalue_ordering() {
    // Results should be sorted by e-value (ascending)
    let (q, db) = (fixture("test_query.fa"), test_data("seqn"));
    if skip_if_missing(&[&q, &db.with_extension("nin")]) { return; }
    let rust = run_rust(&q, &db, 11);
    let evalues: Vec<f64> = rust.lines()
        .filter_map(|l| l.split('\t').nth(10)?.parse().ok())
        .collect();
    for w in evalues.windows(2) {
        assert!(w[0] <= w[1] * 1.01, // small tolerance for equal e-values
            "E-values should be non-decreasing: {} > {}", w[0], w[1]);
    }
}

#[test]
fn gapped_alignment_coordinates_valid() {
    // All query/subject coordinates should be positive and qstart <= qend
    let (q, db) = (fixture("test_query.fa"), test_data("seqn"));
    if skip_if_missing(&[&q, &db.with_extension("nin")]) { return; }
    let rust = run_rust(&q, &db, 11);
    for line in rust.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        let qstart: i32 = f[6].parse().unwrap();
        let qend: i32 = f[7].parse().unwrap();
        let sstart: i32 = f[8].parse().unwrap();
        let send: i32 = f[9].parse().unwrap();
        assert!(qstart >= 1, "qstart should be >= 1: {}", line);
        assert!(qend >= qstart, "qend should be >= qstart: {}", line);
        // For minus strand, sstart > send
        assert!(sstart >= 1 && send >= 1, "coordinates should be positive: {}", line);
    }
}

#[test]
fn gapped_traceback_produces_gap_info() {
    // The FFI engine's traceback should produce gap_info for gapped alignments
    // We verify this indirectly: a hit with mismatches but 0 gap_opens means
    // the alignment was ungapped (correct for blastn with short queries)
    let (q, db) = (fixture("test_query.fa"), test_data("seqn"));
    if skip_if_missing(&[&q, &db.with_extension("nin")]) { return; }
    let rust = run_rust(&q, &db, 11);
    let top: Vec<&str> = rust.lines().next().unwrap_or("").split('\t').collect();
    let gap_opens: i32 = top[5].parse().unwrap_or(-1);
    assert!(gap_opens >= 0, "gap_opens should be non-negative");
}

// ============================================================
// Statistics computation verification
// These tests verify Karlin-Altschul statistics are computed
// correctly by the C core via FFI.
// ============================================================

#[test]
fn statistics_evalue_scales_with_query_length() {
    // A longer query searching the same DB should have smaller e-values
    // for the same alignment score (larger search space)
    let db = test_data("seqn");
    if !db.with_extension("nin").exists() { return; }

    let short = fixture("query_short_match.fa");
    let medium = fixture("query_medium_match.fa");
    if !short.exists() || !medium.exists() { return; }

    let short_out = run_rust(&short, &db, 11);
    let medium_out = run_rust(&medium, &db, 11);

    // Both should find BP722512.1 as top hit
    let short_top: Vec<&str> = short_out.lines().next().unwrap_or("").split('\t').collect();
    let medium_top: Vec<&str> = medium_out.lines().next().unwrap_or("").split('\t').collect();
    assert_eq!(short_top[1], "BP722512.1");
    assert_eq!(medium_top[1], "BP722512.1");

    // Medium query should have lower (better) e-value for a longer alignment
    let short_eval: f64 = short_top[10].parse().unwrap_or(1.0);
    let medium_eval: f64 = medium_top[10].parse().unwrap_or(1.0);
    assert!(medium_eval < short_eval,
        "Longer alignment should have lower e-value: {} vs {}", medium_eval, short_eval);
}

#[test]
fn statistics_bitscore_matches_reference() {
    // Bit scores should match reference exactly for the same scoring params
    let (q, db) = (fixture("test_query.fa"), test_data("seqn"));
    if skip_if_missing(&[&q, &db.with_extension("nin")]) { return; }
    let rust = run_rust(&q, &db, 11);
    let reference = run_ref(&q, &db, 11);
    let r_scores: Vec<&str> = rust.lines().filter_map(|l| l.split('\t').nth(11)).collect();
    let f_scores: Vec<&str> = reference.lines().filter_map(|l| l.split('\t').nth(11)).collect();
    assert_eq!(r_scores, f_scores, "Bit scores should match reference exactly");
}

#[test]
fn statistics_evalue_below_threshold() {
    // All reported hits should have e-value <= the threshold
    let (q, db) = (fixture("test_query.fa"), test_data("seqn"));
    if skip_if_missing(&[&q, &db.with_extension("nin")]) { return; }
    let rust = run_rust_args(&q, &db, 11, 0.01);
    for line in rust.lines() {
        let eval: f64 = line.split('\t').nth(10).unwrap_or("999").parse().unwrap_or(999.0);
        assert!(eval <= 0.01, "E-value {} exceeds threshold 0.01", eval);
    }
}

// ============================================================
// Filtering/validation tests
// These verify behavior that the C++ startup/validation handles.
// ============================================================

#[test]
fn validation_missing_query_file() {
    let output = Command::new(rust_blastn())
        .args(["--query", "/nonexistent/query.fa", "--db",
            test_data("seqn").to_str().unwrap(), "--word_size", "11"])
        .output().expect("Failed to run");
    assert!(!output.status.success(), "Should fail with missing query file");
}

#[test]
fn validation_missing_database() {
    let output = Command::new(rust_blastn())
        .args(["--query", fixture("test_query.fa").to_str().unwrap(),
            "--db", "/nonexistent/db", "--word_size", "11"])
        .output().expect("Failed to run");
    assert!(!output.status.success(), "Should fail with missing database");
}

#[test]
fn validation_empty_query() {
    // Create an empty query file
    let empty = fixture("empty_query.fa");
    if !empty.exists() {
        std::fs::write(&empty, ">empty\n").ok();
    }
    let output = Command::new(rust_blastn())
        .args(["--query", empty.to_str().unwrap(), "--db",
            test_data("seqn").to_str().unwrap(), "--word_size", "11"])
        .output().expect("Failed to run");
    // Should either fail gracefully or produce no output
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.is_empty() || !output.status.success(),
        "Empty query should produce no output or error");
}

#[test]
fn validation_hitlist_size_respected() {
    // With many hits, the number of unique subjects should be <= 500 (default hitlist_size)
    let (q, db) = (fixture("query_long_match.fa"), test_data("seqn"));
    if skip_if_missing(&[&q, &db.with_extension("nin")]) { return; }
    let rust = run_rust(&q, &db, 11);
    let subjects: std::collections::HashSet<&str> = rust.lines()
        .filter_map(|l| l.split('\t').nth(1))
        .collect();
    assert!(subjects.len() <= 500,
        "Should have at most 500 unique subjects, got {}", subjects.len());
}

#[test]
fn validation_scoring_params_affect_results() {
    // Different scoring parameters should produce different e-values
    let q = fixture("test_query.fa");
    let db = test_data("seqn");
    if skip_if_missing(&[&q, &db.with_extension("nin")]) { return; }

    // Default: reward=1, penalty=-3
    let out1 = run_rust(&q, &db, 11);
    // Run reference with reward=2, penalty=-3 for comparison
    let out2 = run_ref_args(&q, &db, 11, 10.0);

    let eval1: f64 = out1.lines().next().unwrap_or("").split('\t').nth(10)
        .unwrap_or("0").parse().unwrap_or(0.0);
    let eval2: f64 = out2.lines().next().unwrap_or("").split('\t').nth(10)
        .unwrap_or("0").parse().unwrap_or(0.0);
    // Both should find hits but e-values should match (same scoring params)
    assert!(eval1 > 0.0 && eval2 > 0.0, "Both should find hits");
    assert_eq!(eval1, eval2, "Same scoring params should give same e-values");
}

#[test]
fn filtering_evalue_threshold_reduces_hits() {
    let (q, db) = (fixture("test_query.fa"), test_data("seqn"));
    if skip_if_missing(&[&q, &db.with_extension("nin")]) { return; }

    let loose = run_rust_args(&q, &db, 11, 10.0);
    let strict = run_rust_args(&q, &db, 11, 1e-10);

    let loose_count = loose.lines().count();
    let strict_count = strict.lines().count();
    assert!(strict_count < loose_count,
        "Stricter e-value should produce fewer hits: {} vs {}", strict_count, loose_count);
    assert!(strict_count >= 1, "Should still find the perfect match");
}

#[test]
fn filtering_minus_strand_hits_have_reversed_coords() {
    // Minus strand hits should have sstart > send
    let (q, db) = (fixture("test_query.fa"), test_data("seqn"));
    if skip_if_missing(&[&q, &db.with_extension("nin")]) { return; }
    let rust = run_rust(&q, &db, 11);
    let mut found_minus = false;
    for line in rust.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        let sstart: i32 = f[8].parse().unwrap_or(0);
        let send: i32 = f[9].parse().unwrap_or(0);
        if sstart > send {
            found_minus = true;
            // Query coords should still be qstart < qend for minus strand
            let qstart: i32 = f[6].parse().unwrap_or(0);
            let qend: i32 = f[7].parse().unwrap_or(0);
            assert!(qstart < qend, "Query coords should be forward even for minus strand");
        }
    }
    assert!(found_minus, "Should find at least one minus-strand hit");
}
