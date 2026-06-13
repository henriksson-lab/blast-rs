//! Integration tests ported from the previous blast-rs implementation.
//!
//! Tests build in-memory BLAST databases, run searches, and validate results.
#![allow(clippy::approx_constant)]

use tempfile::TempDir;

use blast_rs::db::DbType;
use blast_rs::{
    blastn, blastp, blastx, parse_fasta, reverse_complement, six_frame_translate, tblastn, tblastx,
    BlastDbBuilder, BlastnSearch, SearchParams, SearchResult, SequenceEntry, Strand,
};

// ── Helpers ──────────────────────────────────────────────────────────────────

fn blast_cli_bin_for_tests() -> Option<std::path::PathBuf> {
    if let Some(path) = std::env::var_os("BLAST_RS_CLI_BIN")
        .or_else(|| std::env::var_os("CARGO_BIN_EXE_blast-cli"))
        .map(std::path::PathBuf::from)
    {
        return Some(path);
    }

    let exe = std::env::current_exe().ok()?;
    let deps_dir = exe.parent()?;
    let debug_dir = if deps_dir.file_name().is_some_and(|name| name == "deps") {
        deps_dir.parent()?
    } else {
        deps_dir
    };
    let profile_bin = debug_dir.join("blast-cli");
    if profile_bin.exists() {
        return Some(profile_bin);
    }
    debug_dir
        .parent()
        .map(|target_dir| target_dir.join("release").join("blast-cli"))
        .filter(|path| path.exists())
}

fn ncbi_bin(name: &str) -> std::path::PathBuf {
    std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("ncbi-blast-2.17.0+-src/c++/ReleaseMT/bin")
        .join(name)
}

fn run_blastn_subject_cli(args: &[&str], stdin: Option<&[u8]>) -> std::process::Output {
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        panic!("build blast-cli first or set BLAST_RS_CLI_BIN");
    };
    let mut cmd = std::process::Command::new(blast_cli);
    cmd.arg("blastn").args(args);
    if stdin.is_some() {
        cmd.stdin(std::process::Stdio::piped());
    }
    let mut child = cmd.spawn().expect("spawn blast-cli blastn");
    if let Some(input) = stdin {
        use std::io::Write;
        child
            .stdin
            .as_mut()
            .expect("child stdin")
            .write_all(input)
            .expect("write child stdin");
    }
    child.wait_with_output().expect("wait for blast-cli blastn")
}

fn ascii_reverse_complement(seq: &str) -> String {
    String::from_utf8(blast_rs::algo::blast::api::reverse_complement(
        seq.as_bytes(),
    ))
    .unwrap()
}

fn build_protein_db(entries: Vec<SequenceEntry>) -> (TempDir, blast_rs::db::BlastDb) {
    let tmp = TempDir::new().unwrap();
    let base = tmp.path().join("testdb");
    let mut builder = BlastDbBuilder::new(DbType::Protein, "test protein db");
    for e in entries {
        builder.add(e);
    }
    builder.write(&base).unwrap();
    let db = blast_rs::db::BlastDb::open(&base).unwrap();
    (tmp, db)
}

fn build_nucleotide_db(entries: Vec<SequenceEntry>) -> (TempDir, blast_rs::db::BlastDb) {
    let tmp = TempDir::new().unwrap();
    let base = tmp.path().join("testdb");
    let mut builder = BlastDbBuilder::new(DbType::Nucleotide, "test nt db");
    for e in entries {
        builder.add(e);
    }
    builder.write(&base).unwrap();
    let db = blast_rs::db::BlastDb::open(&base).unwrap();
    (tmp, db)
}

fn protein_entry(acc: &str, title: &str, seq: &str) -> SequenceEntry {
    SequenceEntry {
        title: title.to_string(),
        accession: acc.to_string(),
        sequence: seq.as_bytes().to_vec(),
        taxid: None,
    }
}

fn nt_entry(acc: &str, title: &str, seq: &str) -> SequenceEntry {
    SequenceEntry {
        title: title.to_string(),
        accession: acc.to_string(),
        sequence: seq.as_bytes().to_vec(),
        taxid: None,
    }
}

fn assert_blastn_subject_matches_ncbi(
    query_fasta: &str,
    subject_fasta: &str,
    rust_extra_args: &[&str],
    ncbi_extra_args: &[&str],
) {
    assert_blastn_subject_outfmt_matches_ncbi(
        query_fasta,
        subject_fasta,
        "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen",
        rust_extra_args,
        ncbi_extra_args,
    );
}

fn assert_blastn_subject_outfmt_matches_ncbi(
    query_fasta: &str,
    subject_fasta: &str,
    outfmt: &str,
    rust_extra_args: &[&str],
    ncbi_extra_args: &[&str],
) {
    assert_blastn_subject_task_outfmt_matches_ncbi(
        query_fasta,
        subject_fasta,
        "blastn-short",
        outfmt,
        rust_extra_args,
        ncbi_extra_args,
    );
}

fn assert_blastn_subject_task_outfmt_matches_ncbi(
    query_fasta: &str,
    subject_fasta: &str,
    task: &str,
    outfmt: &str,
    rust_extra_args: &[&str],
    ncbi_extra_args: &[&str],
) {
    if !ncbi_bin("blastn").exists() {
        eprintln!("Skipping: blastn not found");
        return;
    }
    let Some(blast_cli) = std::env::var_os("BLAST_RS_CLI_BIN")
        .or_else(|| std::env::var_os("CARGO_BIN_EXE_blast-cli"))
        .map(std::path::PathBuf::from)
    else {
        eprintln!("Skipping: set BLAST_RS_CLI_BIN or CARGO_BIN_EXE_blast-cli to run CLI parity");
        return;
    };

    let tmp = TempDir::new().expect("tempdir");
    let query = tmp.path().join("query.fa");
    let subject = tmp.path().join("subject.fa");
    let rust_out = tmp.path().join("rust.tsv");
    let ncbi_out = tmp.path().join("ncbi.tsv");
    std::fs::write(&query, query_fasta).expect("write query FASTA");
    std::fs::write(&subject, subject_fasta).expect("write subject FASTA");

    let mut rust_cmd = std::process::Command::new(blast_cli);
    rust_cmd
        .arg("blastn")
        .arg("--query")
        .arg(&query)
        .arg("--subject")
        .arg(&subject)
        .arg("--task")
        .arg(task)
        .arg("--outfmt")
        .arg(outfmt)
        .arg("--num_threads")
        .arg("1")
        .arg("--out")
        .arg(&rust_out);
    for arg in rust_extra_args {
        rust_cmd.arg(arg);
    }
    let rust_status = rust_cmd.status().expect("run blast-cli subject parity");
    assert!(
        rust_status.success(),
        "blast-cli exited with {}",
        rust_status
    );

    let mut ncbi_cmd = std::process::Command::new(ncbi_bin("blastn"));
    ncbi_cmd
        .arg("-query")
        .arg(&query)
        .arg("-subject")
        .arg(&subject)
        .arg("-task")
        .arg(task)
        .arg("-outfmt")
        .arg(outfmt)
        .arg("-num_threads")
        .arg("1")
        .arg("-out")
        .arg(&ncbi_out);
    for arg in ncbi_extra_args {
        ncbi_cmd.arg(arg);
    }
    let ncbi_status = ncbi_cmd.status().expect("run NCBI blastn subject parity");
    assert!(
        ncbi_status.success(),
        "NCBI blastn exited with {}",
        ncbi_status
    );

    let rust = std::fs::read(&rust_out).expect("read rust output");
    let ncbi = std::fs::read(&ncbi_out).expect("read ncbi output");
    assert_eq!(
        rust, ncbi,
        "Rust --subject output differs from NCBI\nRust: {:?}\nNCBI: {:?}",
        rust_out, ncbi_out
    );
}

fn assert_blastn_subject_task_outfmt_matches_ncbi_sorted_lines(
    query_fasta: &str,
    subject_fasta: &str,
    task: &str,
    outfmt: &str,
    rust_extra_args: &[&str],
    ncbi_extra_args: &[&str],
) {
    if !ncbi_bin("blastn").exists() {
        eprintln!("Skipping: blastn not found");
        return;
    }
    let Some(blast_cli) = std::env::var_os("BLAST_RS_CLI_BIN")
        .or_else(|| std::env::var_os("CARGO_BIN_EXE_blast-cli"))
        .map(std::path::PathBuf::from)
    else {
        eprintln!("Skipping: set BLAST_RS_CLI_BIN or CARGO_BIN_EXE_blast-cli to run CLI parity");
        return;
    };

    let tmp = TempDir::new().expect("tempdir");
    let query = tmp.path().join("query.fa");
    let subject = tmp.path().join("subject.fa");
    let rust_out = tmp.path().join("rust.tsv");
    let ncbi_out = tmp.path().join("ncbi.tsv");
    std::fs::write(&query, query_fasta).expect("write query FASTA");
    std::fs::write(&subject, subject_fasta).expect("write subject FASTA");

    let mut rust_cmd = std::process::Command::new(blast_cli);
    rust_cmd
        .arg("blastn")
        .arg("--query")
        .arg(&query)
        .arg("--subject")
        .arg(&subject)
        .arg("--task")
        .arg(task)
        .arg("--outfmt")
        .arg(outfmt)
        .arg("--num_threads")
        .arg("1")
        .arg("--out")
        .arg(&rust_out);
    for arg in rust_extra_args {
        rust_cmd.arg(arg);
    }
    let rust_status = rust_cmd.status().expect("run blast-cli subject parity");
    assert!(
        rust_status.success(),
        "blast-cli exited with {}",
        rust_status
    );

    let mut ncbi_cmd = std::process::Command::new(ncbi_bin("blastn"));
    ncbi_cmd
        .arg("-query")
        .arg(&query)
        .arg("-subject")
        .arg(&subject)
        .arg("-task")
        .arg(task)
        .arg("-outfmt")
        .arg(outfmt)
        .arg("-num_threads")
        .arg("1")
        .arg("-out")
        .arg(&ncbi_out);
    for arg in ncbi_extra_args {
        ncbi_cmd.arg(arg);
    }
    let ncbi_status = ncbi_cmd.status().expect("run NCBI blastn subject parity");
    assert!(
        ncbi_status.success(),
        "NCBI blastn exited with {}",
        ncbi_status
    );

    let sort_lines = |bytes: Vec<u8>| {
        let mut lines: Vec<String> = String::from_utf8(bytes)
            .expect("UTF-8 tabular output")
            .lines()
            .map(ToOwned::to_owned)
            .collect();
        lines.sort();
        lines
    };
    let rust = sort_lines(std::fs::read(&rust_out).expect("read rust output"));
    let ncbi = sort_lines(std::fs::read(&ncbi_out).expect("read ncbi output"));
    assert_eq!(
        rust, ncbi,
        "Rust --subject sorted output differs from NCBI\nRust: {:?}\nNCBI: {:?}",
        rust_out, ncbi_out
    );
}

fn assert_blastp_subject_outfmt_matches_ncbi(
    query_fasta: &str,
    subject_fasta: &str,
    outfmt: &str,
    rust_extra_args: &[&str],
    ncbi_extra_args: &[&str],
) {
    if !ncbi_bin("blastp").exists() {
        eprintln!("Skipping: blastp not found");
        return;
    }
    let Some(blast_cli) = std::env::var_os("BLAST_RS_CLI_BIN")
        .or_else(|| std::env::var_os("CARGO_BIN_EXE_blast-cli"))
        .map(std::path::PathBuf::from)
    else {
        eprintln!("Skipping: set BLAST_RS_CLI_BIN or CARGO_BIN_EXE_blast-cli to run CLI parity");
        return;
    };

    let tmp = TempDir::new().expect("tempdir");
    let query = tmp.path().join("query.fa");
    let subject = tmp.path().join("subject.fa");
    let rust_out = tmp.path().join("rust.tsv");
    let ncbi_out = tmp.path().join("ncbi.tsv");
    std::fs::write(&query, query_fasta).expect("write query FASTA");
    std::fs::write(&subject, subject_fasta).expect("write subject FASTA");

    let mut rust_cmd = std::process::Command::new(blast_cli);
    rust_cmd
        .arg("blastp")
        .arg("--query")
        .arg(&query)
        .arg("--subject")
        .arg(&subject)
        .arg("--outfmt")
        .arg(outfmt)
        .arg("--num_threads")
        .arg("1")
        .arg("--out")
        .arg(&rust_out);
    for arg in rust_extra_args {
        rust_cmd.arg(arg);
    }
    let rust_status = rust_cmd
        .status()
        .expect("run blast-cli blastp subject parity");
    assert!(
        rust_status.success(),
        "blast-cli blastp exited with {}",
        rust_status
    );

    let mut ncbi_cmd = std::process::Command::new(ncbi_bin("blastp"));
    ncbi_cmd
        .arg("-query")
        .arg(&query)
        .arg("-subject")
        .arg(&subject)
        .arg("-outfmt")
        .arg(outfmt)
        .arg("-num_threads")
        .arg("1")
        .arg("-out")
        .arg(&ncbi_out);
    for arg in ncbi_extra_args {
        ncbi_cmd.arg(arg);
    }
    let ncbi_status = ncbi_cmd.status().expect("run NCBI blastp subject parity");
    assert!(
        ncbi_status.success(),
        "NCBI blastp exited with {}",
        ncbi_status
    );

    let rust = std::fs::read(&rust_out).expect("read rust output");
    let ncbi = std::fs::read(&ncbi_out).expect("read ncbi output");
    assert_eq!(
        rust, ncbi,
        "Rust blastp --subject output differs from NCBI\nRust: {:?}\nNCBI: {:?}",
        rust_out, ncbi_out
    );
}

fn assert_blastx_subject_outfmt_matches_ncbi(
    query_fasta: &str,
    subject_fasta: &str,
    outfmt: &str,
    rust_extra_args: &[&str],
    ncbi_extra_args: &[&str],
) {
    if !ncbi_bin("blastx").exists() {
        eprintln!("Skipping: blastx not found");
        return;
    }
    let Some(blast_cli) = std::env::var_os("BLAST_RS_CLI_BIN")
        .or_else(|| std::env::var_os("CARGO_BIN_EXE_blast-cli"))
        .map(std::path::PathBuf::from)
    else {
        eprintln!("Skipping: set BLAST_RS_CLI_BIN or CARGO_BIN_EXE_blast-cli to run CLI parity");
        return;
    };

    let tmp = TempDir::new().expect("tempdir");
    let query = tmp.path().join("query.fa");
    let subject = tmp.path().join("subject.fa");
    let rust_out = tmp.path().join("rust.tsv");
    let ncbi_out = tmp.path().join("ncbi.tsv");
    std::fs::write(&query, query_fasta).expect("write query FASTA");
    std::fs::write(&subject, subject_fasta).expect("write subject FASTA");

    let mut rust_cmd = std::process::Command::new(blast_cli);
    rust_cmd
        .arg("blastx")
        .arg("--query")
        .arg(&query)
        .arg("--subject")
        .arg(&subject)
        .arg("--outfmt")
        .arg(outfmt)
        .arg("--num_threads")
        .arg("1")
        .arg("--out")
        .arg(&rust_out);
    for arg in rust_extra_args {
        rust_cmd.arg(arg);
    }
    let rust_status = rust_cmd
        .status()
        .expect("run blast-cli blastx subject parity");
    assert!(
        rust_status.success(),
        "blast-cli blastx exited with {}",
        rust_status
    );

    let mut ncbi_cmd = std::process::Command::new(ncbi_bin("blastx"));
    ncbi_cmd
        .arg("-query")
        .arg(&query)
        .arg("-subject")
        .arg(&subject)
        .arg("-outfmt")
        .arg(outfmt)
        .arg("-num_threads")
        .arg("1")
        .arg("-out")
        .arg(&ncbi_out);
    for arg in ncbi_extra_args {
        ncbi_cmd.arg(arg);
    }
    let ncbi_status = ncbi_cmd.status().expect("run NCBI blastx subject parity");
    assert!(
        ncbi_status.success(),
        "NCBI blastx exited with {}",
        ncbi_status
    );

    let rust = std::fs::read(&rust_out).expect("read rust output");
    let ncbi = std::fs::read(&ncbi_out).expect("read ncbi output");
    assert_eq!(
        rust, ncbi,
        "Rust blastx --subject output differs from NCBI\nRust: {:?}\nNCBI: {:?}",
        rust_out, ncbi_out
    );
}

fn assert_blastp_db_outfmt_matches_ncbi(
    query_fasta: &str,
    db_fasta: &str,
    outfmt: &str,
    rust_extra_args: &[&str],
    ncbi_extra_args: &[&str],
) {
    if !ncbi_bin("blastp").exists() || !ncbi_bin("makeblastdb").exists() {
        eprintln!("Skipping: blastp or makeblastdb not found");
        return;
    }
    let Some(blast_cli) = std::env::var_os("BLAST_RS_CLI_BIN")
        .or_else(|| std::env::var_os("CARGO_BIN_EXE_blast-cli"))
        .map(std::path::PathBuf::from)
    else {
        eprintln!("Skipping: set BLAST_RS_CLI_BIN or CARGO_BIN_EXE_blast-cli to run CLI parity");
        return;
    };

    let tmp = TempDir::new().expect("tempdir");
    let query = tmp.path().join("query.fa");
    let db_fasta_path = tmp.path().join("db.fa");
    let db = tmp.path().join("testdb");
    let rust_out = tmp.path().join("rust.txt");
    let ncbi_out = tmp.path().join("ncbi.txt");
    std::fs::write(&query, query_fasta).expect("write query FASTA");
    std::fs::write(&db_fasta_path, db_fasta).expect("write db FASTA");

    let make_status = std::process::Command::new(ncbi_bin("makeblastdb"))
        .arg("-in")
        .arg(&db_fasta_path)
        .arg("-dbtype")
        .arg("prot")
        .arg("-out")
        .arg(&db)
        .stdout(std::process::Stdio::null())
        .status()
        .expect("run makeblastdb");
    assert!(
        make_status.success(),
        "makeblastdb exited with {make_status}"
    );

    let mut rust_cmd = std::process::Command::new(blast_cli);
    rust_cmd
        .arg("blastp")
        .arg("--query")
        .arg(&query)
        .arg("--db")
        .arg(&db)
        .arg("--outfmt")
        .arg(outfmt)
        .arg("--num_threads")
        .arg("1")
        .arg("--out")
        .arg(&rust_out);
    for arg in rust_extra_args {
        rust_cmd.arg(arg);
    }
    let rust_status = rust_cmd.status().expect("run blast-cli blastp DB parity");
    assert!(
        rust_status.success(),
        "blast-cli blastp exited with {}",
        rust_status
    );

    let mut ncbi_cmd = std::process::Command::new(ncbi_bin("blastp"));
    ncbi_cmd
        .arg("-query")
        .arg(&query)
        .arg("-db")
        .arg(&db)
        .arg("-outfmt")
        .arg(outfmt)
        .arg("-num_threads")
        .arg("1")
        .arg("-out")
        .arg(&ncbi_out);
    for arg in ncbi_extra_args {
        ncbi_cmd.arg(arg);
    }
    let ncbi_status = ncbi_cmd.status().expect("run NCBI blastp DB parity");
    assert!(
        ncbi_status.success(),
        "NCBI blastp exited with {}",
        ncbi_status
    );

    let rust = std::fs::read(&rust_out).expect("read rust output");
    let ncbi = std::fs::read(&ncbi_out).expect("read ncbi output");
    assert_eq!(
        rust, ncbi,
        "Rust blastp DB output differs from NCBI\nRust: {:?}\nNCBI: {:?}",
        rust_out, ncbi_out
    );
}

fn assert_translated_subject_outfmt_matches_ncbi(
    program: &str,
    ncbi_program: &str,
    query_fasta: &str,
    subject_fasta: &str,
    outfmt: &str,
    rust_extra_args: &[&str],
    ncbi_extra_args: &[&str],
) {
    if !std::path::Path::new(ncbi_program).exists() {
        eprintln!("Skipping: {ncbi_program} not found");
        return;
    }
    let Some(blast_cli) = std::env::var_os("BLAST_RS_CLI_BIN")
        .or_else(|| std::env::var_os("CARGO_BIN_EXE_blast-cli"))
        .map(std::path::PathBuf::from)
    else {
        eprintln!("Skipping: set BLAST_RS_CLI_BIN or CARGO_BIN_EXE_blast-cli to run CLI parity");
        return;
    };

    let tmp = TempDir::new().expect("tempdir");
    let query = tmp.path().join("query.fa");
    let subject = tmp.path().join("subject.fa");
    let rust_out = tmp.path().join("rust.tsv");
    let ncbi_out = tmp.path().join("ncbi.tsv");
    std::fs::write(&query, query_fasta).expect("write query FASTA");
    std::fs::write(&subject, subject_fasta).expect("write subject FASTA");

    let mut rust_cmd = std::process::Command::new(blast_cli);
    rust_cmd
        .arg(program)
        .arg("--query")
        .arg(&query)
        .arg("--subject")
        .arg(&subject)
        .arg("--outfmt")
        .arg(outfmt)
        .arg("--num_threads")
        .arg("1")
        .arg("--out")
        .arg(&rust_out);
    for arg in rust_extra_args {
        rust_cmd.arg(arg);
    }
    let rust_status = rust_cmd.status().expect("run blast-cli translated parity");
    assert!(
        rust_status.success(),
        "blast-cli {program} exited with {}",
        rust_status
    );

    let mut ncbi_cmd = std::process::Command::new(ncbi_program);
    ncbi_cmd
        .arg("-query")
        .arg(&query)
        .arg("-subject")
        .arg(&subject)
        .arg("-outfmt")
        .arg(outfmt)
        .arg("-num_threads")
        .arg("1")
        .arg("-out")
        .arg(&ncbi_out);
    for arg in ncbi_extra_args {
        ncbi_cmd.arg(arg);
    }
    let ncbi_status = ncbi_cmd.status().expect("run NCBI translated parity");
    assert!(
        ncbi_status.success(),
        "NCBI {program} exited with {}",
        ncbi_status
    );

    let rust = std::fs::read(&rust_out).expect("read rust output");
    let ncbi = std::fs::read(&ncbi_out).expect("read ncbi output");
    assert_eq!(
        rust, ncbi,
        "Rust {program} --subject output differs from NCBI\nRust: {:?}\nNCBI: {:?}",
        rust_out, ncbi_out
    );
}

fn assert_translated_subject_outfmt_matches_ncbi_sorted_lines(
    program: &str,
    ncbi_program: &str,
    query_fasta: &str,
    subject_fasta: &str,
    outfmt: &str,
    rust_extra_args: &[&str],
    ncbi_extra_args: &[&str],
) {
    if !std::path::Path::new(ncbi_program).exists() {
        eprintln!("Skipping: {ncbi_program} not found");
        return;
    }
    let Some(blast_cli) = std::env::var_os("BLAST_RS_CLI_BIN")
        .or_else(|| std::env::var_os("CARGO_BIN_EXE_blast-cli"))
        .map(std::path::PathBuf::from)
    else {
        eprintln!("Skipping: set BLAST_RS_CLI_BIN or CARGO_BIN_EXE_blast-cli to run CLI parity");
        return;
    };

    let tmp = TempDir::new().expect("tempdir");
    let query = tmp.path().join("query.fa");
    let subject = tmp.path().join("subject.fa");
    let rust_out = tmp.path().join("rust.tsv");
    let ncbi_out = tmp.path().join("ncbi.tsv");
    std::fs::write(&query, query_fasta).expect("write query FASTA");
    std::fs::write(&subject, subject_fasta).expect("write subject FASTA");

    let mut rust_cmd = std::process::Command::new(blast_cli);
    rust_cmd
        .arg(program)
        .arg("--query")
        .arg(&query)
        .arg("--subject")
        .arg(&subject)
        .arg("--outfmt")
        .arg(outfmt)
        .arg("--num_threads")
        .arg("1")
        .arg("--out")
        .arg(&rust_out);
    for arg in rust_extra_args {
        rust_cmd.arg(arg);
    }
    let rust_status = rust_cmd.status().expect("run blast-cli translated parity");
    assert!(
        rust_status.success(),
        "blast-cli {program} exited with {}",
        rust_status
    );

    let mut ncbi_cmd = std::process::Command::new(ncbi_program);
    ncbi_cmd
        .arg("-query")
        .arg(&query)
        .arg("-subject")
        .arg(&subject)
        .arg("-outfmt")
        .arg(outfmt)
        .arg("-num_threads")
        .arg("1")
        .arg("-out")
        .arg(&ncbi_out);
    for arg in ncbi_extra_args {
        ncbi_cmd.arg(arg);
    }
    let ncbi_status = ncbi_cmd.status().expect("run NCBI translated parity");
    assert!(
        ncbi_status.success(),
        "NCBI {program} exited with {}",
        ncbi_status
    );

    let sort_lines = |bytes: Vec<u8>| {
        let mut lines: Vec<String> = String::from_utf8(bytes)
            .expect("UTF-8 tabular output")
            .lines()
            .map(ToOwned::to_owned)
            .collect();
        lines.sort();
        lines
    };

    let rust = sort_lines(std::fs::read(&rust_out).expect("read rust output"));
    let ncbi = sort_lines(std::fs::read(&ncbi_out).expect("read ncbi output"));
    assert_eq!(
        rust, ncbi,
        "Rust {program} --subject sorted output differs from NCBI\nRust: {:?}\nNCBI: {:?}",
        rust_out, ncbi_out
    );
}

fn assert_translated_db_outfmt_matches_ncbi_sorted_lines(
    program: &str,
    ncbi_program: &str,
    dbtype: &str,
    query_fasta: &str,
    db_fasta: &str,
    outfmt: &str,
    rust_extra_args: &[&str],
    ncbi_extra_args: &[&str],
) {
    if !std::path::Path::new(ncbi_program).exists() || !ncbi_bin("makeblastdb").exists() {
        eprintln!("Skipping: {ncbi_program} or makeblastdb not found");
        return;
    }
    let Some(blast_cli) = std::env::var_os("BLAST_RS_CLI_BIN")
        .or_else(|| std::env::var_os("CARGO_BIN_EXE_blast-cli"))
        .map(std::path::PathBuf::from)
    else {
        eprintln!("Skipping: set BLAST_RS_CLI_BIN or CARGO_BIN_EXE_blast-cli to run CLI parity");
        return;
    };

    let tmp = TempDir::new().expect("tempdir");
    let query = tmp.path().join("query.fa");
    let db_fasta_path = tmp.path().join("db.fa");
    let db = tmp.path().join("testdb");
    let rust_out = tmp.path().join("rust.tsv");
    let ncbi_out = tmp.path().join("ncbi.tsv");
    std::fs::write(&query, query_fasta).expect("write query FASTA");
    std::fs::write(&db_fasta_path, db_fasta).expect("write db FASTA");

    let make_status = std::process::Command::new(ncbi_bin("makeblastdb"))
        .arg("-in")
        .arg(&db_fasta_path)
        .arg("-dbtype")
        .arg(dbtype)
        .arg("-out")
        .arg(&db)
        .stdout(std::process::Stdio::null())
        .status()
        .expect("run makeblastdb");
    assert!(
        make_status.success(),
        "makeblastdb exited with {make_status}"
    );

    let mut rust_cmd = std::process::Command::new(blast_cli);
    rust_cmd
        .arg(program)
        .arg("--query")
        .arg(&query)
        .arg("--db")
        .arg(&db)
        .arg("--outfmt")
        .arg(outfmt)
        .arg("--num_threads")
        .arg("1")
        .arg("--out")
        .arg(&rust_out);
    for arg in rust_extra_args {
        rust_cmd.arg(arg);
    }
    let rust_status = rust_cmd
        .status()
        .expect("run blast-cli translated DB parity");
    assert!(
        rust_status.success(),
        "blast-cli {program} exited with {}",
        rust_status
    );

    let mut ncbi_cmd = std::process::Command::new(ncbi_program);
    ncbi_cmd
        .arg("-query")
        .arg(&query)
        .arg("-db")
        .arg(&db)
        .arg("-outfmt")
        .arg(outfmt)
        .arg("-num_threads")
        .arg("1")
        .arg("-out")
        .arg(&ncbi_out);
    for arg in ncbi_extra_args {
        ncbi_cmd.arg(arg);
    }
    let ncbi_status = ncbi_cmd.status().expect("run NCBI translated DB parity");
    assert!(
        ncbi_status.success(),
        "NCBI {program} exited with {}",
        ncbi_status
    );

    let sort_lines = |bytes: Vec<u8>| {
        let mut lines: Vec<String> = String::from_utf8(bytes)
            .expect("UTF-8 tabular output")
            .lines()
            .map(ToOwned::to_owned)
            .collect();
        lines.sort();
        lines
    };

    let rust = sort_lines(std::fs::read(&rust_out).expect("read rust output"));
    let ncbi = sort_lines(std::fs::read(&ncbi_out).expect("read ncbi output"));
    assert_eq!(
        rust, ncbi,
        "Rust {program} DB sorted output differs from NCBI\nRust: {:?}\nNCBI: {:?}",
        rust_out, ncbi_out
    );
}

fn assert_translated_db_outfmt_matches_ncbi(
    program: &str,
    ncbi_program: &str,
    dbtype: &str,
    query_fasta: &str,
    db_fasta: &str,
    outfmt: &str,
    rust_extra_args: &[&str],
    ncbi_extra_args: &[&str],
) {
    assert_translated_db_outfmt_matches_ncbi_with_num_threads(
        program,
        ncbi_program,
        dbtype,
        query_fasta,
        db_fasta,
        outfmt,
        rust_extra_args,
        ncbi_extra_args,
        "1",
    );
}

fn assert_translated_db_outfmt_matches_ncbi_with_num_threads(
    program: &str,
    ncbi_program: &str,
    dbtype: &str,
    query_fasta: &str,
    db_fasta: &str,
    outfmt: &str,
    rust_extra_args: &[&str],
    ncbi_extra_args: &[&str],
    num_threads: &str,
) {
    if !std::path::Path::new(ncbi_program).exists() || !ncbi_bin("makeblastdb").exists() {
        eprintln!("Skipping: {ncbi_program} or makeblastdb not found");
        return;
    }
    let Some(blast_cli) = std::env::var_os("BLAST_RS_CLI_BIN")
        .or_else(|| std::env::var_os("CARGO_BIN_EXE_blast-cli"))
        .map(std::path::PathBuf::from)
    else {
        eprintln!("Skipping: set BLAST_RS_CLI_BIN or CARGO_BIN_EXE_blast-cli to run CLI parity");
        return;
    };

    let tmp = TempDir::new().expect("tempdir");
    let query = tmp.path().join("query.fa");
    let db_fasta_path = tmp.path().join("db.fa");
    let db = tmp.path().join("testdb");
    let rust_out = tmp.path().join("rust.tsv");
    let ncbi_out = tmp.path().join("ncbi.tsv");
    std::fs::write(&query, query_fasta).expect("write query FASTA");
    std::fs::write(&db_fasta_path, db_fasta).expect("write db FASTA");

    let make_status = std::process::Command::new(ncbi_bin("makeblastdb"))
        .arg("-in")
        .arg(&db_fasta_path)
        .arg("-dbtype")
        .arg(dbtype)
        .arg("-out")
        .arg(&db)
        .stdout(std::process::Stdio::null())
        .status()
        .expect("run makeblastdb");
    assert!(
        make_status.success(),
        "makeblastdb exited with {make_status}"
    );

    let mut rust_cmd = std::process::Command::new(blast_cli);
    rust_cmd
        .arg(program)
        .arg("--query")
        .arg(&query)
        .arg("--db")
        .arg(&db)
        .arg("--outfmt")
        .arg(outfmt)
        .arg("--num_threads")
        .arg(num_threads)
        .arg("--out")
        .arg(&rust_out);
    for arg in rust_extra_args {
        rust_cmd.arg(arg);
    }
    let rust_status = rust_cmd
        .status()
        .expect("run blast-cli translated DB parity");
    assert!(
        rust_status.success(),
        "blast-cli {program} exited with {}",
        rust_status
    );

    let mut ncbi_cmd = std::process::Command::new(ncbi_program);
    ncbi_cmd
        .arg("-query")
        .arg(&query)
        .arg("-db")
        .arg(&db)
        .arg("-outfmt")
        .arg(outfmt)
        .arg("-num_threads")
        .arg(num_threads)
        .arg("-out")
        .arg(&ncbi_out);
    for arg in ncbi_extra_args {
        ncbi_cmd.arg(arg);
    }
    let ncbi_status = ncbi_cmd.status().expect("run NCBI translated DB parity");
    assert!(
        ncbi_status.success(),
        "NCBI {program} exited with {}",
        ncbi_status
    );

    let rust = std::fs::read(&rust_out).expect("read rust output");
    let ncbi = std::fs::read(&ncbi_out).expect("read ncbi output");
    assert_eq!(
        rust, ncbi,
        "Rust {program} DB output differs from NCBI\nRust: {:?}\nNCBI: {:?}",
        rust_out, ncbi_out
    );
}

fn normalize_sam_for_cli_parity(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes)
        .lines()
        .map(|line| {
            if line.starts_with("@PG\t") {
                "@PG\tID:0\tVN:2.17.0+\tCL:<normalized>\tPN:blastn".to_string()
            } else {
                line.to_string()
            }
        })
        .collect::<Vec<_>>()
        .join("\n")
        + "\n"
}

fn assert_blastn_subject_sam_matches_ncbi(
    query_fasta: &str,
    subject_fasta: &str,
    rust_extra_args: &[&str],
    ncbi_extra_args: &[&str],
) {
    if !ncbi_bin("blastn").exists() {
        eprintln!("Skipping: blastn not found");
        return;
    }
    let Some(blast_cli) = std::env::var_os("BLAST_RS_CLI_BIN")
        .or_else(|| std::env::var_os("CARGO_BIN_EXE_blast-cli"))
        .map(std::path::PathBuf::from)
    else {
        eprintln!("Skipping: set BLAST_RS_CLI_BIN or CARGO_BIN_EXE_blast-cli to run CLI parity");
        return;
    };

    let tmp = TempDir::new().expect("tempdir");
    let query = tmp.path().join("query.fa");
    let subject = tmp.path().join("subject.fa");
    let rust_out = tmp.path().join("rust.sam");
    let ncbi_out = tmp.path().join("ncbi.sam");
    std::fs::write(&query, query_fasta).expect("write query FASTA");
    std::fs::write(&subject, subject_fasta).expect("write subject FASTA");

    let mut rust_cmd = std::process::Command::new(blast_cli);
    rust_cmd
        .arg("blastn")
        .arg("--query")
        .arg(&query)
        .arg("--subject")
        .arg(&subject)
        .arg("--task")
        .arg("blastn-short")
        .arg("--outfmt")
        .arg("17")
        .arg("--num_threads")
        .arg("1")
        .arg("--out")
        .arg(&rust_out);
    for arg in rust_extra_args {
        rust_cmd.arg(arg);
    }
    let rust_status = rust_cmd.status().expect("run blast-cli SAM parity");
    assert!(
        rust_status.success(),
        "blast-cli exited with {}",
        rust_status
    );

    let mut ncbi_cmd = std::process::Command::new(ncbi_bin("blastn"));
    ncbi_cmd
        .arg("-query")
        .arg(&query)
        .arg("-subject")
        .arg(&subject)
        .arg("-task")
        .arg("blastn-short")
        .arg("-outfmt")
        .arg("17")
        .arg("-num_threads")
        .arg("1")
        .arg("-out")
        .arg(&ncbi_out);
    for arg in ncbi_extra_args {
        ncbi_cmd.arg(arg);
    }
    let ncbi_status = ncbi_cmd.status().expect("run NCBI blastn SAM parity");
    assert!(
        ncbi_status.success(),
        "NCBI blastn exited with {}",
        ncbi_status
    );

    let rust = std::fs::read(&rust_out).expect("read rust SAM output");
    let ncbi = std::fs::read(&ncbi_out).expect("read ncbi SAM output");
    assert_eq!(
        normalize_sam_for_cli_parity(&rust),
        normalize_sam_for_cli_parity(&ncbi),
        "Rust --subject SAM output differs from NCBI\nRust: {:?}\nNCBI: {:?}",
        rust_out,
        ncbi_out
    );
}

fn assert_blastn_db_sam_matches_ncbi(
    query_fasta: &str,
    db_fasta: &str,
    rust_extra_args: &[&str],
    ncbi_extra_args: &[&str],
) {
    if !ncbi_bin("blastn").exists() || !ncbi_bin("makeblastdb").exists() {
        eprintln!("Skipping: blastn or makeblastdb not found");
        return;
    }
    let Some(blast_cli) = std::env::var_os("BLAST_RS_CLI_BIN")
        .or_else(|| std::env::var_os("CARGO_BIN_EXE_blast-cli"))
        .map(std::path::PathBuf::from)
    else {
        eprintln!("Skipping: set BLAST_RS_CLI_BIN or CARGO_BIN_EXE_blast-cli to run CLI parity");
        return;
    };

    let tmp = TempDir::new().expect("tempdir");
    let query = tmp.path().join("query.fa");
    let db_fasta_path = tmp.path().join("db.fa");
    let db = tmp.path().join("testdb");
    let rust_out = tmp.path().join("rust.sam");
    let ncbi_out = tmp.path().join("ncbi.sam");
    std::fs::write(&query, query_fasta).expect("write query FASTA");
    std::fs::write(&db_fasta_path, db_fasta).expect("write db FASTA");

    let make_status = std::process::Command::new(ncbi_bin("makeblastdb"))
        .arg("-in")
        .arg(&db_fasta_path)
        .arg("-dbtype")
        .arg("nucl")
        .arg("-out")
        .arg(&db)
        .status()
        .expect("run makeblastdb");
    assert!(
        make_status.success(),
        "makeblastdb exited with {make_status}"
    );

    let mut rust_cmd = std::process::Command::new(blast_cli);
    rust_cmd
        .arg("blastn")
        .arg("--query")
        .arg(&query)
        .arg("--db")
        .arg(&db)
        .arg("--task")
        .arg("blastn-short")
        .arg("--outfmt")
        .arg("17")
        .arg("--num_threads")
        .arg("1")
        .arg("--out")
        .arg(&rust_out);
    for arg in rust_extra_args {
        rust_cmd.arg(arg);
    }
    let rust_status = rust_cmd.status().expect("run blast-cli DB SAM parity");
    assert!(
        rust_status.success(),
        "blast-cli exited with {}",
        rust_status
    );

    let mut ncbi_cmd = std::process::Command::new(ncbi_bin("blastn"));
    ncbi_cmd
        .arg("-query")
        .arg(&query)
        .arg("-db")
        .arg(&db)
        .arg("-task")
        .arg("blastn-short")
        .arg("-outfmt")
        .arg("17")
        .arg("-num_threads")
        .arg("1")
        .arg("-out")
        .arg(&ncbi_out);
    for arg in ncbi_extra_args {
        ncbi_cmd.arg(arg);
    }
    let ncbi_status = ncbi_cmd.status().expect("run NCBI blastn DB SAM parity");
    assert!(
        ncbi_status.success(),
        "NCBI blastn exited with {}",
        ncbi_status
    );

    let rust = std::fs::read(&rust_out).expect("read rust SAM output");
    let ncbi = std::fs::read(&ncbi_out).expect("read ncbi SAM output");
    assert_eq!(
        normalize_sam_for_cli_parity(&rust),
        normalize_sam_for_cli_parity(&ncbi),
        "Rust DB SAM output differs from NCBI\nRust: {:?}\nNCBI: {:?}",
        rust_out,
        ncbi_out
    );
}

fn assert_blastn_db_xml_matches_ncbi(
    query_fasta: &str,
    db_fasta: &str,
    rust_extra_args: &[&str],
    ncbi_extra_args: &[&str],
) {
    if !ncbi_bin("blastn").exists() || !ncbi_bin("makeblastdb").exists() {
        eprintln!("Skipping: blastn or makeblastdb not found");
        return;
    }
    let Some(blast_cli) = std::env::var_os("BLAST_RS_CLI_BIN")
        .or_else(|| std::env::var_os("CARGO_BIN_EXE_blast-cli"))
        .map(std::path::PathBuf::from)
    else {
        eprintln!("Skipping: set BLAST_RS_CLI_BIN or CARGO_BIN_EXE_blast-cli to run CLI parity");
        return;
    };

    let tmp = TempDir::new().expect("tempdir");
    let query = tmp.path().join("query.fa");
    let db_fasta_path = tmp.path().join("db.fa");
    let db = tmp.path().join("testdb");
    let rust_out = tmp.path().join("rust.xml");
    let ncbi_out = tmp.path().join("ncbi.xml");
    std::fs::write(&query, query_fasta).expect("write query FASTA");
    std::fs::write(&db_fasta_path, db_fasta).expect("write db FASTA");

    let make_status = std::process::Command::new(ncbi_bin("makeblastdb"))
        .arg("-in")
        .arg(&db_fasta_path)
        .arg("-dbtype")
        .arg("nucl")
        .arg("-out")
        .arg(&db)
        .stdout(std::process::Stdio::null())
        .status()
        .expect("run makeblastdb");
    assert!(
        make_status.success(),
        "makeblastdb exited with {make_status}"
    );

    let mut rust_cmd = std::process::Command::new(blast_cli);
    rust_cmd
        .arg("blastn")
        .arg("--query")
        .arg(&query)
        .arg("--db")
        .arg(&db)
        .arg("--task")
        .arg("blastn-short")
        .arg("--outfmt")
        .arg("5")
        .arg("--num_threads")
        .arg("1")
        .arg("--out")
        .arg(&rust_out);
    for arg in rust_extra_args {
        rust_cmd.arg(arg);
    }
    let rust_status = rust_cmd.status().expect("run blast-cli DB XML parity");
    assert!(
        rust_status.success(),
        "blast-cli exited with {}",
        rust_status
    );

    let mut ncbi_cmd = std::process::Command::new(ncbi_bin("blastn"));
    ncbi_cmd
        .arg("-query")
        .arg(&query)
        .arg("-db")
        .arg(&db)
        .arg("-task")
        .arg("blastn-short")
        .arg("-outfmt")
        .arg("5")
        .arg("-num_threads")
        .arg("1")
        .arg("-out")
        .arg(&ncbi_out);
    for arg in ncbi_extra_args {
        ncbi_cmd.arg(arg);
    }
    let ncbi_status = ncbi_cmd.status().expect("run NCBI blastn DB XML parity");
    assert!(
        ncbi_status.success(),
        "NCBI blastn exited with {}",
        ncbi_status
    );

    let rust = std::fs::read(&rust_out).expect("read rust XML output");
    let ncbi = std::fs::read(&ncbi_out).expect("read ncbi XML output");
    assert_eq!(
        rust, ncbi,
        "Rust DB XML output differs from NCBI\nRust: {:?}\nNCBI: {:?}",
        rust_out, ncbi_out
    );
}

fn assert_blastn_db_outfmt_matches_ncbi(
    query_fasta: &str,
    db_fasta: &str,
    outfmt: &str,
    makeblastdb_extra_args: &[&str],
    rust_extra_args: &[&str],
    ncbi_extra_args: &[&str],
) {
    assert_blastn_db_outfmt_matches_ncbi_with_num_threads(
        query_fasta,
        db_fasta,
        outfmt,
        makeblastdb_extra_args,
        rust_extra_args,
        ncbi_extra_args,
        "1",
    );
}

fn assert_blastn_db_outfmt_matches_ncbi_with_num_threads(
    query_fasta: &str,
    db_fasta: &str,
    outfmt: &str,
    makeblastdb_extra_args: &[&str],
    rust_extra_args: &[&str],
    ncbi_extra_args: &[&str],
    num_threads: &str,
) {
    if !ncbi_bin("blastn").exists() || !ncbi_bin("makeblastdb").exists() {
        eprintln!("Skipping: blastn or makeblastdb not found");
        return;
    }
    let Some(blast_cli) = std::env::var_os("BLAST_RS_CLI_BIN")
        .or_else(|| std::env::var_os("CARGO_BIN_EXE_blast-cli"))
        .map(std::path::PathBuf::from)
    else {
        eprintln!("Skipping: set BLAST_RS_CLI_BIN or CARGO_BIN_EXE_blast-cli to run CLI parity");
        return;
    };

    let tmp = TempDir::new().expect("tempdir");
    let query = tmp.path().join("query.fa");
    let db_fasta_path = tmp.path().join("db.fa");
    let db = tmp.path().join("testdb");
    let rust_out = tmp.path().join("rust.tsv");
    let ncbi_out = tmp.path().join("ncbi.tsv");
    let taxids_file = tmp.path().join("taxids.txt");
    let seqids_file = tmp.path().join("seqids.txt");
    std::fs::write(&query, query_fasta).expect("write query FASTA");
    std::fs::write(&db_fasta_path, db_fasta).expect("write db FASTA");
    std::fs::write(&taxids_file, "9606\n").expect("write taxid list");
    std::fs::write(&seqids_file, "s2\n").expect("write seqid list");
    let rust_has_task = rust_extra_args
        .iter()
        .any(|arg| matches!(*arg, "--task" | "-task"));
    let ncbi_has_task = ncbi_extra_args
        .iter()
        .any(|arg| matches!(*arg, "--task" | "-task"));

    let mut make_cmd = std::process::Command::new(ncbi_bin("makeblastdb"));
    make_cmd
        .arg("-in")
        .arg(&db_fasta_path)
        .arg("-dbtype")
        .arg("nucl")
        .arg("-out")
        .arg(&db)
        .stdout(std::process::Stdio::null());
    for arg in makeblastdb_extra_args {
        make_cmd.arg(arg);
    }
    let make_status = make_cmd.status().expect("run makeblastdb");
    assert!(
        make_status.success(),
        "makeblastdb exited with {make_status}"
    );

    let mut rust_cmd = std::process::Command::new(blast_cli);
    rust_cmd
        .arg("blastn")
        .arg("--query")
        .arg(&query)
        .arg("--db")
        .arg(&db);
    if !rust_has_task {
        rust_cmd.arg("--task").arg("blastn-short");
    }
    rust_cmd
        .arg("--outfmt")
        .arg(outfmt)
        .arg("--num_threads")
        .arg(num_threads)
        .arg("--out")
        .arg(&rust_out);
    for arg in rust_extra_args {
        if *arg == "{taxids_file}" {
            rust_cmd.arg(&taxids_file);
        } else if *arg == "{seqids_file}" {
            rust_cmd.arg(&seqids_file);
        } else {
            rust_cmd.arg(arg);
        }
    }
    let rust_status = rust_cmd.status().expect("run blast-cli DB parity");
    assert!(
        rust_status.success(),
        "blast-cli exited with {}",
        rust_status
    );

    let mut ncbi_cmd = std::process::Command::new(ncbi_bin("blastn"));
    ncbi_cmd.arg("-query").arg(&query).arg("-db").arg(&db);
    if !ncbi_has_task {
        ncbi_cmd.arg("-task").arg("blastn-short");
    }
    ncbi_cmd
        .arg("-outfmt")
        .arg(outfmt)
        .arg("-num_threads")
        .arg(num_threads)
        .arg("-out")
        .arg(&ncbi_out);
    for arg in ncbi_extra_args {
        if *arg == "{taxids_file}" {
            ncbi_cmd.arg(&taxids_file);
        } else if *arg == "{seqids_file}" {
            ncbi_cmd.arg(&seqids_file);
        } else {
            ncbi_cmd.arg(arg);
        }
    }
    let ncbi_status = ncbi_cmd.status().expect("run NCBI blastn DB parity");
    assert!(
        ncbi_status.success(),
        "NCBI blastn exited with {}",
        ncbi_status
    );

    let rust = std::fs::read(&rust_out).expect("read rust output");
    let ncbi = std::fs::read(&ncbi_out).expect("read ncbi output");
    assert_eq!(
        rust, ncbi,
        "Rust DB output differs from NCBI\nRust: {:?}\nNCBI: {:?}",
        rust_out, ncbi_out
    );
}

fn assert_blastn_db_outfmt_matches_ncbi_sorted_lines(
    query_fasta: &str,
    db_fasta: &str,
    outfmt: &str,
    makeblastdb_extra_args: &[&str],
    rust_extra_args: &[&str],
    ncbi_extra_args: &[&str],
) {
    if !ncbi_bin("blastn").exists() || !ncbi_bin("makeblastdb").exists() {
        eprintln!("Skipping: blastn or makeblastdb not found");
        return;
    }
    let Some(blast_cli) = std::env::var_os("BLAST_RS_CLI_BIN")
        .or_else(|| std::env::var_os("CARGO_BIN_EXE_blast-cli"))
        .map(std::path::PathBuf::from)
    else {
        eprintln!("Skipping: set BLAST_RS_CLI_BIN or CARGO_BIN_EXE_blast-cli to run CLI parity");
        return;
    };

    let tmp = TempDir::new().expect("tempdir");
    let query = tmp.path().join("query.fa");
    let db_fasta_path = tmp.path().join("db.fa");
    let db = tmp.path().join("testdb");
    let rust_out = tmp.path().join("rust.tsv");
    let ncbi_out = tmp.path().join("ncbi.tsv");
    std::fs::write(&query, query_fasta).expect("write query FASTA");
    std::fs::write(&db_fasta_path, db_fasta).expect("write db FASTA");
    let rust_has_task = rust_extra_args
        .iter()
        .any(|arg| matches!(*arg, "--task" | "-task"));
    let ncbi_has_task = ncbi_extra_args
        .iter()
        .any(|arg| matches!(*arg, "--task" | "-task"));

    let mut make_cmd = std::process::Command::new(ncbi_bin("makeblastdb"));
    make_cmd
        .arg("-in")
        .arg(&db_fasta_path)
        .arg("-dbtype")
        .arg("nucl")
        .arg("-out")
        .arg(&db)
        .stdout(std::process::Stdio::null());
    for arg in makeblastdb_extra_args {
        make_cmd.arg(arg);
    }
    let make_status = make_cmd.status().expect("run makeblastdb");
    assert!(
        make_status.success(),
        "makeblastdb exited with {make_status}"
    );

    let mut rust_cmd = std::process::Command::new(blast_cli);
    rust_cmd
        .arg("blastn")
        .arg("--query")
        .arg(&query)
        .arg("--db")
        .arg(&db);
    if !rust_has_task {
        rust_cmd.arg("--task").arg("blastn-short");
    }
    rust_cmd
        .arg("--outfmt")
        .arg(outfmt)
        .arg("--num_threads")
        .arg("1")
        .arg("--out")
        .arg(&rust_out);
    for arg in rust_extra_args {
        rust_cmd.arg(arg);
    }
    let rust_status = rust_cmd.status().expect("run blast-cli DB parity");
    assert!(
        rust_status.success(),
        "blast-cli exited with {}",
        rust_status
    );

    let mut ncbi_cmd = std::process::Command::new(ncbi_bin("blastn"));
    ncbi_cmd.arg("-query").arg(&query).arg("-db").arg(&db);
    if !ncbi_has_task {
        ncbi_cmd.arg("-task").arg("blastn-short");
    }
    ncbi_cmd
        .arg("-outfmt")
        .arg(outfmt)
        .arg("-num_threads")
        .arg("1")
        .arg("-out")
        .arg(&ncbi_out);
    for arg in ncbi_extra_args {
        ncbi_cmd.arg(arg);
    }
    let ncbi_status = ncbi_cmd.status().expect("run NCBI blastn DB parity");
    assert!(
        ncbi_status.success(),
        "NCBI blastn exited with {}",
        ncbi_status
    );

    let sort_lines = |bytes: Vec<u8>| {
        let mut lines: Vec<String> = String::from_utf8(bytes)
            .expect("UTF-8 tabular output")
            .lines()
            .map(ToOwned::to_owned)
            .collect();
        lines.sort();
        lines
    };
    let rust = sort_lines(std::fs::read(&rust_out).expect("read rust output"));
    let ncbi = sort_lines(std::fs::read(&ncbi_out).expect("read ncbi output"));
    assert_eq!(
        rust, ncbi,
        "Rust DB sorted output differs from NCBI\nRust: {:?}\nNCBI: {:?}",
        rust_out, ncbi_out
    );
}

fn assert_blastn_indexed_db_outfmt_matches_ncbi(
    query_fasta: &str,
    db_fasta: &str,
    outfmt: &str,
    rust_extra_args: &[&str],
    ncbi_extra_args: &[&str],
) {
    if !ncbi_bin("blastn").exists()
        || !ncbi_bin("makeblastdb").exists()
        || !ncbi_bin("makembindex").exists()
    {
        eprintln!("Skipping: blastn, makeblastdb, or makembindex not found");
        return;
    }
    let Some(blast_cli) = std::env::var_os("BLAST_RS_CLI_BIN")
        .or_else(|| std::env::var_os("CARGO_BIN_EXE_blast-cli"))
        .map(std::path::PathBuf::from)
    else {
        eprintln!("Skipping: set BLAST_RS_CLI_BIN or CARGO_BIN_EXE_blast-cli to run CLI parity");
        return;
    };

    let tmp = TempDir::new().expect("tempdir");
    let query = tmp.path().join("query.fa");
    let db_fasta_path = tmp.path().join("db.fa");
    let db = tmp.path().join("testdb");
    let rust_out = tmp.path().join("rust.tsv");
    let ncbi_out = tmp.path().join("ncbi.tsv");
    std::fs::write(&query, query_fasta).expect("write query FASTA");
    std::fs::write(&db_fasta_path, db_fasta).expect("write db FASTA");

    let make_status = std::process::Command::new(ncbi_bin("makeblastdb"))
        .arg("-in")
        .arg(&db_fasta_path)
        .arg("-dbtype")
        .arg("nucl")
        .arg("-out")
        .arg(&db)
        .stdout(std::process::Stdio::null())
        .status()
        .expect("run makeblastdb");
    assert!(
        make_status.success(),
        "makeblastdb exited with {make_status}"
    );

    let index_status = std::process::Command::new(ncbi_bin("makembindex"))
        .arg("-input")
        .arg(&db)
        .arg("-iformat")
        .arg("blastdb")
        .arg("-verbosity")
        .arg("quiet")
        .stdout(std::process::Stdio::null())
        .status()
        .expect("run makembindex");
    assert!(
        index_status.success(),
        "makembindex exited with {index_status}"
    );

    let mut rust_cmd = std::process::Command::new(blast_cli);
    rust_cmd
        .arg("blastn")
        .arg("--query")
        .arg(&query)
        .arg("--db")
        .arg(&db)
        .arg("--outfmt")
        .arg(outfmt)
        .arg("--num_threads")
        .arg("1")
        .arg("--out")
        .arg(&rust_out);
    for arg in rust_extra_args {
        if *arg == "{db}" {
            rust_cmd.arg(&db);
        } else {
            rust_cmd.arg(arg);
        }
    }
    let rust_status = rust_cmd.status().expect("run blast-cli indexed DB parity");
    assert!(
        rust_status.success(),
        "blast-cli exited with {}",
        rust_status
    );

    let mut ncbi_cmd = std::process::Command::new(ncbi_bin("blastn"));
    ncbi_cmd
        .arg("-query")
        .arg(&query)
        .arg("-db")
        .arg(&db)
        .arg("-outfmt")
        .arg(outfmt)
        .arg("-num_threads")
        .arg("1")
        .arg("-out")
        .arg(&ncbi_out);
    for arg in ncbi_extra_args {
        if *arg == "{db}" {
            ncbi_cmd.arg(&db);
        } else {
            ncbi_cmd.arg(arg);
        }
    }
    let ncbi_status = ncbi_cmd.status().expect("run NCBI indexed DB parity");
    assert!(
        ncbi_status.success(),
        "NCBI blastn exited with {}",
        ncbi_status
    );

    let rust = std::fs::read(&rust_out).expect("read rust output");
    let ncbi = std::fs::read(&ncbi_out).expect("read ncbi output");
    assert_eq!(
        rust, ncbi,
        "Rust indexed DB output differs from NCBI\nRust: {:?}\nNCBI: {:?}",
        rust_out, ncbi_out
    );
}

fn large_db_fixture_paths(query_name: &str) -> Option<(std::path::PathBuf, std::path::PathBuf)> {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/large_db");
    let query = root.join(query_name);
    let db = root.join("celegans");
    if query.exists() && db.with_extension("nin").exists() {
        Some((query, db))
    } else {
        eprintln!(
            "Skipping: large_db fixture not present under {}",
            root.display()
        );
        None
    }
}

fn run_large_db_blastn(
    query_name: &str,
    max_hsps: Option<&str>,
    extra_env: Option<(&str, &str)>,
) -> Option<(TempDir, std::path::PathBuf)> {
    let (query, db) = large_db_fixture_paths(query_name)?;
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli or set BLAST_RS_CLI_BIN to run CLI parity");
        return None;
    };

    let tmp = TempDir::new().expect("tempdir");
    let out = tmp.path().join("out.tsv");
    let mut cmd = std::process::Command::new(blast_cli);
    cmd.arg("blastn")
        .arg("--task")
        .arg("blastn-short")
        .arg("--dust")
        .arg("no")
        .arg("--evalue")
        .arg("10")
        .arg("--query")
        .arg(query)
        .arg("--db")
        .arg(db)
        .arg("--outfmt")
        .arg("6 qseqid sseqid pident length qstart qend sstart send bitscore evalue")
        .arg("--num_threads")
        .arg("1")
        .arg("--out")
        .arg(&out);
    if let Some(max_hsps) = max_hsps {
        cmd.arg("--max_hsps").arg(max_hsps);
    }
    if let Some((key, value)) = extra_env {
        cmd.env(key, value);
    }
    let status = cmd.status().expect("run blast-cli large_db parity");
    assert!(status.success(), "blast-cli exited with {status}");
    Some((tmp, out))
}

fn assert_large_db_blastn_matches_ncbi(query_name: &str, max_hsps: Option<&str>) {
    if !ncbi_bin("blastn").exists() {
        eprintln!("Skipping: blastn not found");
        return;
    }
    let Some((query, db)) = large_db_fixture_paths(query_name) else {
        return;
    };
    let Some((rust_tmp, rust_out)) = run_large_db_blastn(query_name, max_hsps, None) else {
        return;
    };

    let ncbi_tmp = TempDir::new().expect("tempdir");
    let ncbi_out = ncbi_tmp.path().join("ncbi.tsv");
    let mut ncbi_cmd = std::process::Command::new(ncbi_bin("blastn"));
    ncbi_cmd
        .arg("-task")
        .arg("blastn-short")
        .arg("-dust")
        .arg("no")
        .arg("-evalue")
        .arg("10")
        .arg("-query")
        .arg(query)
        .arg("-db")
        .arg(db)
        .arg("-outfmt")
        .arg("6 qseqid sseqid pident length qstart qend sstart send bitscore evalue")
        .arg("-num_threads")
        .arg("1")
        .arg("-out")
        .arg(&ncbi_out);
    if let Some(max_hsps) = max_hsps {
        ncbi_cmd.arg("-max_hsps").arg(max_hsps);
    }
    let status = ncbi_cmd.status().expect("run NCBI blastn large_db parity");
    assert!(status.success(), "NCBI blastn exited with {status}");

    let rust = std::fs::read(&rust_out).expect("read rust output");
    let ncbi = std::fs::read(&ncbi_out).expect("read ncbi output");
    assert_eq!(
        rust, ncbi,
        "Rust large_db output differs from NCBI\nRust: {:?}\nNCBI: {:?}",
        rust_out, ncbi_out
    );
    drop(rust_tmp);
}

fn patch_blastdb_index_total_length(index_path: &std::path::Path, total_length: u64) {
    fn read_u32_be(data: &[u8], offset: &mut usize) -> u32 {
        let value = u32::from_be_bytes([
            data[*offset],
            data[*offset + 1],
            data[*offset + 2],
            data[*offset + 3],
        ]);
        *offset += 4;
        value
    }
    fn skip_blast_string(data: &[u8], offset: &mut usize) {
        let len = read_u32_be(data, offset) as usize;
        *offset += len;
    }

    let mut data = std::fs::read(index_path).expect("read BLAST index");
    let mut offset = 0usize;
    let version = read_u32_be(&data, &mut offset);
    let _seq_type = read_u32_be(&data, &mut offset);
    if version == 5 {
        let _volume_number = read_u32_be(&data, &mut offset);
    }
    skip_blast_string(&data, &mut offset);
    if version == 5 {
        skip_blast_string(&data, &mut offset);
    }
    skip_blast_string(&data, &mut offset);
    let _num_oids = read_u32_be(&data, &mut offset);
    let total_offset = offset;
    data[total_offset..total_offset + 8].copy_from_slice(&total_length.to_le_bytes());
    std::fs::write(index_path, data).expect("write patched BLAST index");
}

fn write_taxonomy4blast_sqlite(dir: &std::path::Path) {
    let db = dir.join("taxonomy4blast.sqlite3");
    let conn = rusqlite::Connection::open(&db).expect("open taxonomy4blast.sqlite3");
    conn.execute(
        "CREATE TABLE TaxidInfo (taxid INTEGER PRIMARY KEY, parent INTEGER)",
        [],
    )
    .expect("create TaxidInfo");
    conn.execute(
        "CREATE INDEX TaxidInfoCompositeIdx_parent ON TaxidInfo(parent,taxid)",
        [],
    )
    .expect("create TaxidInfo parent index");
    for (taxid, parent) in [(1, 1), (9605, 1), (9606, 9605), (63221, 9606)] {
        conn.execute(
            "INSERT INTO TaxidInfo(taxid,parent) VALUES (?1,?2)",
            (taxid, parent),
        )
        .expect("insert taxonomy row");
    }
}

fn run_rust_taxonomy_filter_case(extra_args: &[&str]) -> Vec<u8> {
    if !ncbi_bin("makeblastdb").exists() {
        eprintln!("Skipping: makeblastdb not found");
        return Vec::new();
    }
    let Some(blast_cli) = std::env::var_os("BLAST_RS_CLI_BIN")
        .or_else(|| std::env::var_os("CARGO_BIN_EXE_blast-cli"))
        .map(std::path::PathBuf::from)
    else {
        eprintln!("Skipping: set BLAST_RS_CLI_BIN or CARGO_BIN_EXE_blast-cli to run CLI parity");
        return Vec::new();
    };

    let tmp = TempDir::new().expect("tempdir");
    let query = tmp.path().join("query.fa");
    let db_fasta_path = tmp.path().join("db.fa");
    let taxid_map = tmp.path().join("taxid_map.txt");
    let db = tmp.path().join("testdb");
    let out = tmp.path().join("rust.tsv");
    std::fs::write(&query, ">q1\nACGTACGTACGTACGTACGT\n").expect("write query FASTA");
    std::fs::write(
        &db_fasta_path,
        ">s9606\nACGTACGTACGTACGTACGT\n>s63221\nACGTACGTACGTACGTACGT\n",
    )
    .expect("write db FASTA");
    std::fs::write(&taxid_map, "s9606 9606\ns63221 63221\n").expect("write taxid map");
    write_taxonomy4blast_sqlite(tmp.path());

    let make_status = std::process::Command::new(ncbi_bin("makeblastdb"))
        .arg("-in")
        .arg(&db_fasta_path)
        .arg("-dbtype")
        .arg("nucl")
        .arg("-out")
        .arg(&db)
        .arg("-parse_seqids")
        .arg("-taxid_map")
        .arg(&taxid_map)
        .stdout(std::process::Stdio::null())
        .status()
        .expect("run makeblastdb");
    assert!(
        make_status.success(),
        "makeblastdb exited with {make_status}"
    );

    let mut cmd = std::process::Command::new(blast_cli);
    cmd.arg("blastn")
        .arg("--query")
        .arg(&query)
        .arg("--db")
        .arg(&db)
        .arg("--task")
        .arg("blastn-short")
        .arg("--outfmt")
        .arg("6 sseqid staxid")
        .arg("--num_threads")
        .arg("1")
        .arg("--dust")
        .arg("no")
        .arg("--max_target_seqs")
        .arg("10")
        .arg("--max_hsps")
        .arg("1")
        .arg("--out")
        .arg(&out)
        .env("BLASTDB", tmp.path());
    for arg in extra_args {
        cmd.arg(arg);
    }
    let status = cmd.status().expect("run blast-cli taxonomy expansion case");
    assert!(status.success(), "blast-cli exited with {status}");
    std::fs::read(&out).expect("read taxonomy expansion output")
}

fn pssm_round_output_path(base: &std::path::Path, round: usize) -> std::path::PathBuf {
    let mut path = base.as_os_str().to_os_string();
    path.push(format!(".{round}"));
    std::path::PathBuf::from(path)
}

fn assert_psiblast_subject_outfmt_matches_ncbi(
    query_fasta: &str,
    subject_fasta: &str,
    outfmt: &str,
    rust_extra: &[&str],
    ncbi_extra: &[&str],
) {
    if !ncbi_bin("psiblast").exists() {
        eprintln!("Skipping: psiblast not found");
        return;
    }
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    let tmp = TempDir::new().expect("tempdir");
    let query = tmp.path().join("query.fa");
    let subject = tmp.path().join("subject.fa");
    std::fs::write(&query, query_fasta).expect("write query");
    std::fs::write(&subject, subject_fasta).expect("write subject");

    let mut rust_cmd = std::process::Command::new(&blast_cli);
    rust_cmd
        .arg("psiblast")
        .arg("--query")
        .arg(&query)
        .arg("--subject")
        .arg(&subject)
        .arg("--outfmt")
        .arg(outfmt);
    rust_cmd.args(rust_extra);
    let rust = rust_cmd
        .output()
        .unwrap_or_else(|err| panic!("run blast-cli psiblast {outfmt}: {err}"));

    let mut ncbi_cmd = std::process::Command::new(ncbi_bin("psiblast"));
    ncbi_cmd
        .arg("-query")
        .arg(&query)
        .arg("-subject")
        .arg(&subject)
        .arg("-outfmt")
        .arg(outfmt);
    ncbi_cmd.args(ncbi_extra);
    let ncbi = ncbi_cmd
        .output()
        .unwrap_or_else(|err| panic!("run NCBI psiblast {outfmt}: {err}"));

    assert!(rust.status.success(), "blast-cli psiblast failed");
    assert!(ncbi.status.success(), "NCBI psiblast failed");
    assert_eq!(rust.stdout, ncbi.stdout, "psiblast {outfmt} stdout differs");
    assert_eq!(
        String::from_utf8_lossy(&rust.stderr),
        String::from_utf8_lossy(&ncbi.stderr),
        "psiblast {outfmt} stderr differs"
    );
}

// ── BLASTP tests ─────────────────────────────────────────────────────────────

// ── BLASTN tests ─────────────────────────────────────────────────────────────

// ── BLASTX test ─────────────────────────────────────────────────────────────

// ── TBLASTN test ────────────────────────────────────────────────────────────

// ── TBLASTX test ────────────────────────────────────────────────────────────

// ── Database round-trip tests ────────────────────────────────────────────────

// ── Reverse complement test ─────────────────────────────────────────────────

// ── Six-frame translation test ──────────────────────────────────────────────

// ── FASTA parsing edge cases ─────────────────────────────────────────────────

// ── Multi-query search ──────────────────────────────────────────────────────

// ── Edge cases ──────────────────────────────────────────────────────────────

// ── Large query test ────────────────────────────────────────────────────────

// ── Blastn with multiple mismatches ─────────────────────────────────────────

// ── Multi-subject index tests ───────────────────────────────────────────────

// ── Prokka-style annotation performance test ─────────────────────────────────

/// Performance test: build a protein DB from prokka's sprot FASTA (~25K entries),
/// then search 5 query proteins against it. This mimics the prokka annotation
/// workload. Run with: cargo test --release -- --ignored test_blastp_prokka_sprot
///
/// Target: should complete in under 30 seconds for 5 queries (Perl Prokka does
/// 63 queries against the same DB in ~12 seconds using NCBI BLAST+).
/// Sensitivity test: verify that blastp finds the same hits as NCBI BLAST+.
///
/// These are real CDS proteins from E. faecium plasmid AUS0004_p1 that
/// NCBI BLAST+ (via Perl Prokka) annotates against the sprot database.
/// blast-rs must find the same top hit and hit count for each of these queries.
///
/// Run with: cargo test --release -- --ignored test_blastp_sensitivity
/// Test composition-based matrix adjustment on the srtA/Sortase A fixture.
/// The query at position 00009 in the E. faecium plasmid is a legitimate
/// P0DPQ5 Sortase A hit; composition adjustment should retain it while
/// reshaping the alignment and score consistently with the NCBI reference.
///
/// Run with: cargo test --release -- --ignored test_comp_adjust_srta
/// Real-life test: compare blast-rs against NCBI BLAST+ output for known prokka queries.
/// Tests that top hit, score, alignment length, and e-value match NCBI BLAST+.
/// Each query was run through NCBI blastp 2.12.0+ with comp_based_stats=0 against
/// the Bacteria/sprot database and the top hit recorded.
///
/// Run with: cargo test --release -- --ignored test_blastp_vs_ncbi
/// Real-life test: verify hit count matches NCBI for the full plasmid annotation.
/// Runs all 63 CDS from the test plasmid against sprot and counts annotated hits.
/// NCBI BLAST+ with comp_based_stats=0, evalue=1e-9 finds ~11 hits for this dataset.
///
/// Run with: cargo test --release -- --ignored test_blastp_plasmid_annotation_count
/// Stress test for composition-based statistics lambda ratio.
/// Uses compositionally biased queries (Pro-rich, Glu-rich, Ala-rich) that
/// produce large mode 0 vs mode 1 score differences in NCBI BLAST+.
/// This test exposes bugs in lambda ratio computation — the score should
/// DECREASE (or stay similar) with comp_adjust=1 vs comp_adjust=0 for biased sequences.
///
/// NCBI BLAST+ 2.12.0 reference values (comp_based_stats=0 → comp_based_stats=1):
///   pro_rich: score 164→105 (top hit changes from P9WJC5→Q70XJ9)
///   glu_rich: score 97→56  (36% reduction)
///   ala_rich: score 79→58  (27% reduction)
///   srtA:     score 234→221, alignment 225→141
///   normal:   score 162→166 (slight increase, balanced composition)
///
/// Run with: cargo test --release -- --ignored test_lambda_ratio_stress
/// External sensitivity fixture for srtA and umuC against Prokka sprot.
// ── End-to-end API tests (ported from NCBI bl2seq / blastengine / traceback) ─

/// Search a nucleotide sequence against itself. Should produce a single perfect
/// alignment covering the full length with 100% identity.
/// Search completely unrelated sequences and verify no hits at strict evalue.
/// Subject contains two separate regions matching the query. Verify multiple
/// HSPs are returned.
/// Protein self-search: search a protein against itself and verify perfect alignment.
/// Search two related but not identical proteins. Verify a hit with positive
/// score and reasonable e-value.
/// BLASTX: nucleotide query encoding a known protein searched against a protein
/// database. Should find the translated match.
/// TBLASTN: protein query against nucleotide subject that encodes the protein.
/// Should find the translated nucleotide match.
/// Query whose reverse complement matches the subject. Verify hit on minus strand.
/// Very short query (15bp) with word_size=7. Should still find a match.
/// Set a strict evalue threshold and verify only significant hits pass.
/// Search against a BLAST database and against a subject FASTA should produce
/// equivalent results for the same sequences.
/// Use BlastnSearch builder API directly for a seq-vs-seq search.
/// External srtA/P0DPQ5 composition-adjustment internals.
// ── Stress tests for short-primer / large-DB scenarios (stack overflow regression) ──

/// Generate a random-ish nucleotide sequence of given length using a simple LCG.
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

/// Full Swiss-Prot blastp benchmark: build a protein DB from UniProt Swiss-Prot
/// (~570K entries), search 100 query sequences, and compare results with NCBI BLAST+.
///
/// Requires Swiss-Prot FASTA at one of the checked paths.
/// Run with: cargo test --release -- --ignored test_blastp_swissprot
/// External lambda-ratio reference values for biased sequences.
/// Run with: cargo test --release -- --ignored test_lambda_ratio_biased_sprot_reference --nocapture
/// Verify that blast_karlin_lambda_nr with Robinson&Robinson frequencies gives the
/// correct standard BLOSUM62 ungapped lambda (0.3176).
// ── Real core_nt parity tests ────────────────────────────────────────────────

const CORE_NT_PRIMER_ID: &str = "realistic_primer";
const CORE_NT_PRIMER_SEQ: &str = "GTCTCCTCTGACTTCAACAGCG";
const CORE_NT_BASE: &str = "/husky/henriksson/for_claude/blast/core_nt/core_nt";

fn blast_cli_bin() -> std::path::PathBuf {
    std::env::var_os("CARGO_BIN_EXE_blast-cli")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|| {
            let profile = if cfg!(debug_assertions) {
                "debug"
            } else {
                "release"
            };
            std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                .join("target")
                .join(profile)
                .join("blast-cli")
        })
}

fn write_core_nt_primer(tmp: &TempDir) -> std::path::PathBuf {
    let query = tmp.path().join("core_nt_realistic_primer.fa");
    std::fs::write(
        &query,
        format!(">{}\n{}\n", CORE_NT_PRIMER_ID, CORE_NT_PRIMER_SEQ),
    )
    .expect("write primer query");
    query
}

fn run_core_nt_rust(query: &std::path::Path, db: &std::path::Path, out: &std::path::Path) {
    let status = std::process::Command::new(blast_cli_bin())
        .arg("blastn")
        .arg("--query")
        .arg(query)
        .arg("--db")
        .arg(db)
        .arg("--task")
        .arg("blastn-short")
        .arg("--outfmt")
        .arg("6")
        .arg("--max_target_seqs")
        .arg("500")
        .arg("--num_threads")
        .arg("8")
        .arg("--out")
        .arg(out)
        .status()
        .expect("run blast-cli");
    assert!(status.success(), "blast-cli exited with {}", status);
}

fn run_core_nt_ncbi(query: &std::path::Path, db: &std::path::Path, out: &std::path::Path) {
    let status = std::process::Command::new(ncbi_bin("blastn"))
        .arg("-query")
        .arg(query)
        .arg("-db")
        .arg(db)
        .arg("-task")
        .arg("blastn-short")
        .arg("-outfmt")
        .arg("6")
        .arg("-max_target_seqs")
        .arg("500")
        .arg("-num_threads")
        .arg("8")
        .arg("-out")
        .arg(out)
        .status()
        .expect("run blastn");
    assert!(status.success(), "NCBI blastn exited with {}", status);
}

fn assert_core_nt_outfmt6_matches_ncbi(db_suffix: &str) {
    let db = std::path::PathBuf::from(format!("{}{}", CORE_NT_BASE, db_suffix));
    let db_index = {
        let mut p = db.as_os_str().to_os_string();
        p.push(".nin");
        std::path::PathBuf::from(p)
    };
    if !db_index.exists() && !db.with_extension("nal").exists() {
        eprintln!("Skipping: core_nt database not found at {:?}", db);
        return;
    }
    if !ncbi_bin("blastn").exists() {
        eprintln!("Skipping: blastn not found");
        return;
    }

    let tmp = TempDir::new().expect("tempdir");
    let query = write_core_nt_primer(&tmp);
    let rust_out = tmp.path().join("rust.tsv");
    let ncbi_out = tmp.path().join("ncbi.tsv");

    run_core_nt_rust(&query, &db, &rust_out);
    run_core_nt_ncbi(&query, &db, &ncbi_out);

    let rust = std::fs::read(&rust_out).expect("read rust output");
    let ncbi = std::fs::read(&ncbi_out).expect("read ncbi output");
    assert_eq!(
        rust, ncbi,
        "Rust outfmt 6 output differs from NCBI for {:?}\nRust: {:?}\nNCBI: {:?}",
        db, rust_out, ncbi_out
    );
}

/// Real core_nt parity check for the first volume.
///
/// Run with: cargo test --release -- --ignored test_core_nt00_primer_outfmt6_matches_ncbi
fn assert_core_nt_taxonomy_outfmt_matches_ncbi(db_suffix: &str) {
    assert_core_nt_taxonomy_outfmt_matches_ncbi_with_blastdb(db_suffix, None);
}

fn assert_core_nt_taxonomy_outfmt_matches_ncbi_with_blastdb(
    db_suffix: &str,
    blastdb: Option<&std::path::Path>,
) {
    let db = std::path::PathBuf::from(format!("{}{}", CORE_NT_BASE, db_suffix));
    let db_index = {
        let mut p = db.as_os_str().to_os_string();
        p.push(".nin");
        std::path::PathBuf::from(p)
    };
    if !db_index.exists() && !db.with_extension("nal").exists() {
        eprintln!("Skipping: core_nt database not found at {:?}", db);
        return;
    }
    if !ncbi_bin("blastn").exists() {
        eprintln!("Skipping: blastn not found");
        return;
    }

    let tmp = TempDir::new().expect("tempdir");
    let query = write_core_nt_primer(&tmp);
    let rust_out = tmp.path().join("rust_tax.tsv");
    let ncbi_out = tmp.path().join("ncbi_tax.tsv");
    let outfmt = "6 qseqid qlen qstart qend saccver slen sstart send bitscore staxid ssciname sskingdom length pident";

    let mut rust_cmd = std::process::Command::new(blast_cli_bin());
    if let Some(blastdb) = blastdb {
        rust_cmd.env("BLASTDB", blastdb);
    }
    let rust_status = rust_cmd
        .arg("blastn")
        .arg("--task")
        .arg("blastn-short")
        .arg("--evalue")
        .arg("5")
        .arg("--max_target_seqs")
        .arg("10000")
        .arg("--max_hsps")
        .arg("1")
        .arg("--num_threads")
        .arg("8")
        .arg("--perc_identity")
        .arg("90")
        .arg("--db")
        .arg(&db)
        .arg("--outfmt")
        .arg(outfmt)
        .arg("--query")
        .arg(&query)
        .arg("--out")
        .arg(&rust_out)
        .status()
        .expect("run blast-cli taxonomy regression");
    assert!(
        rust_status.success(),
        "blast-cli exited with {}",
        rust_status
    );

    let mut ncbi_cmd = std::process::Command::new(ncbi_bin("blastn"));
    if let Some(blastdb) = blastdb {
        ncbi_cmd.env("BLASTDB", blastdb);
    }
    let ncbi_status = ncbi_cmd
        .arg("-task")
        .arg("blastn-short")
        .arg("-evalue")
        .arg("5")
        .arg("-max_target_seqs")
        .arg("10000")
        .arg("-max_hsps")
        .arg("1")
        .arg("-num_threads")
        .arg("8")
        .arg("-perc_identity")
        .arg("90")
        .arg("-db")
        .arg(&db)
        .arg("-outfmt")
        .arg(outfmt)
        .arg("-query")
        .arg(&query)
        .arg("-out")
        .arg(&ncbi_out)
        .status()
        .expect("run NCBI blastn taxonomy regression");
    assert!(
        ncbi_status.success(),
        "NCBI blastn exited with {}",
        ncbi_status
    );

    let rust = std::fs::read(&rust_out).expect("read rust taxonomy output");
    let ncbi = std::fs::read(&ncbi_out).expect("read ncbi taxonomy output");
    assert_eq!(
        rust, ncbi,
        "Rust taxonomy outfmt output differs from NCBI\nRust: {:?}\nNCBI: {:?}",
        rust_out, ncbi_out
    );
}

/// Regression test for the user-reported taxonomy-heavy primer command.
///
/// Run with: cargo test --release -- --ignored test_core_nt00_primer_taxonomy_outfmt_matches_ncbi
fn assert_core_nt_title_outfmt_matches_ncbi(db_suffix: &str) {
    let db = std::path::PathBuf::from(format!("{}{}", CORE_NT_BASE, db_suffix));
    let db_index = {
        let mut p = db.as_os_str().to_os_string();
        p.push(".nin");
        std::path::PathBuf::from(p)
    };
    if !db_index.exists() && !db.with_extension("nal").exists() {
        eprintln!("Skipping: core_nt database not found at {:?}", db);
        return;
    }
    if !ncbi_bin("blastn").exists() {
        eprintln!("Skipping: blastn not found");
        return;
    }

    let tmp = TempDir::new().expect("tempdir");
    let query = write_core_nt_primer(&tmp);
    let rust_out = tmp.path().join("rust_title.tsv");
    let ncbi_out = tmp.path().join("ncbi_title.tsv");
    let outfmt = "6 qseqid saccver stitle salltitles bitscore";

    let rust_status = std::process::Command::new(blast_cli_bin())
        .arg("blastn")
        .arg("--query")
        .arg(&query)
        .arg("--db")
        .arg(&db)
        .arg("--task")
        .arg("blastn-short")
        .arg("--outfmt")
        .arg(outfmt)
        .arg("--max_target_seqs")
        .arg("20")
        .arg("--num_threads")
        .arg("8")
        .arg("--out")
        .arg(&rust_out)
        .status()
        .expect("run blast-cli title parity");
    assert!(
        rust_status.success(),
        "blast-cli exited with {}",
        rust_status
    );

    let ncbi_status = std::process::Command::new(ncbi_bin("blastn"))
        .arg("-query")
        .arg(&query)
        .arg("-db")
        .arg(&db)
        .arg("-task")
        .arg("blastn-short")
        .arg("-outfmt")
        .arg(outfmt)
        .arg("-max_target_seqs")
        .arg("20")
        .arg("-num_threads")
        .arg("8")
        .arg("-out")
        .arg(&ncbi_out)
        .status()
        .expect("run blastn title parity");
    assert!(
        ncbi_status.success(),
        "NCBI blastn exited with {}",
        ncbi_status
    );

    let rust = std::fs::read(&rust_out).expect("read rust title output");
    let ncbi = std::fs::read(&ncbi_out).expect("read ncbi title output");
    assert_eq!(
        rust, ncbi,
        "Rust title outfmt output differs from NCBI\nRust: {:?}\nNCBI: {:?}",
        rust_out, ncbi_out
    );
}

/// Real core_nt title-field parity check for the first volume.
///
/// Run with: cargo test --release -- --ignored test_core_nt00_primer_title_outfmt_matches_ncbi
fn assert_core_nt_outfmt7_matches_ncbi(db_suffix: &str) {
    let db = std::path::PathBuf::from(format!("{}{}", CORE_NT_BASE, db_suffix));
    let db_index = {
        let mut p = db.as_os_str().to_os_string();
        p.push(".nin");
        std::path::PathBuf::from(p)
    };
    if !db_index.exists() && !db.with_extension("nal").exists() {
        eprintln!("Skipping: core_nt database not found at {:?}", db);
        return;
    }
    if !ncbi_bin("blastn").exists() {
        eprintln!("Skipping: blastn not found");
        return;
    }

    let tmp = TempDir::new().expect("tempdir");
    let query = write_core_nt_primer(&tmp);
    let rust_out = tmp.path().join("rust_outfmt7.tsv");
    let ncbi_out = tmp.path().join("ncbi_outfmt7.tsv");

    let rust_status = std::process::Command::new(blast_cli_bin())
        .arg("blastn")
        .arg("--query")
        .arg(&query)
        .arg("--db")
        .arg(&db)
        .arg("--task")
        .arg("blastn-short")
        .arg("--outfmt")
        .arg("7")
        .arg("--max_target_seqs")
        .arg("5")
        .arg("--num_threads")
        .arg("8")
        .arg("--out")
        .arg(&rust_out)
        .status()
        .expect("run blast-cli outfmt7 parity");
    assert!(
        rust_status.success(),
        "blast-cli exited with {}",
        rust_status
    );

    let ncbi_status = std::process::Command::new(ncbi_bin("blastn"))
        .arg("-query")
        .arg(&query)
        .arg("-db")
        .arg(&db)
        .arg("-task")
        .arg("blastn-short")
        .arg("-outfmt")
        .arg("7")
        .arg("-max_target_seqs")
        .arg("5")
        .arg("-num_threads")
        .arg("8")
        .arg("-out")
        .arg(&ncbi_out)
        .status()
        .expect("run blastn outfmt7 parity");
    assert!(
        ncbi_status.success(),
        "NCBI blastn exited with {}",
        ncbi_status
    );

    let rust = std::fs::read(&rust_out).expect("read rust outfmt7 output");
    let ncbi = std::fs::read(&ncbi_out).expect("read ncbi outfmt7 output");
    assert_eq!(
        rust, ncbi,
        "Rust outfmt 7 output differs from NCBI\nRust: {:?}\nNCBI: {:?}",
        rust_out, ncbi_out
    );
}
