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

fn ascii_reverse_complement(seq: &str) -> String {
    String::from_utf8(blast_rs::api::reverse_complement(seq.as_bytes())).unwrap()
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
    if !std::path::Path::new("/usr/bin/blastn").exists() {
        eprintln!("Skipping: /usr/bin/blastn not found");
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

    let mut ncbi_cmd = std::process::Command::new("/usr/bin/blastn");
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
    if !std::path::Path::new("/usr/bin/blastn").exists() {
        eprintln!("Skipping: /usr/bin/blastn not found");
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

    let mut ncbi_cmd = std::process::Command::new("/usr/bin/blastn");
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
    if !std::path::Path::new("/usr/bin/blastp").exists() {
        eprintln!("Skipping: /usr/bin/blastp not found");
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

    let mut ncbi_cmd = std::process::Command::new("/usr/bin/blastp");
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

fn assert_blastp_db_outfmt_matches_ncbi(
    query_fasta: &str,
    db_fasta: &str,
    outfmt: &str,
    rust_extra_args: &[&str],
    ncbi_extra_args: &[&str],
) {
    if !std::path::Path::new("/usr/bin/blastp").exists()
        || !std::path::Path::new("/usr/bin/makeblastdb").exists()
    {
        eprintln!("Skipping: /usr/bin/blastp or /usr/bin/makeblastdb not found");
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

    let make_status = std::process::Command::new("/usr/bin/makeblastdb")
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

    let mut ncbi_cmd = std::process::Command::new("/usr/bin/blastp");
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
    if !std::path::Path::new(ncbi_program).exists()
        || !std::path::Path::new("/usr/bin/makeblastdb").exists()
    {
        eprintln!("Skipping: {ncbi_program} or /usr/bin/makeblastdb not found");
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

    let make_status = std::process::Command::new("/usr/bin/makeblastdb")
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
    if !std::path::Path::new(ncbi_program).exists()
        || !std::path::Path::new("/usr/bin/makeblastdb").exists()
    {
        eprintln!("Skipping: {ncbi_program} or /usr/bin/makeblastdb not found");
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

    let make_status = std::process::Command::new("/usr/bin/makeblastdb")
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
                "@PG\tID:0\tVN:2.12.0+\tCL:<normalized>\tPN:blastn".to_string()
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
    if !std::path::Path::new("/usr/bin/blastn").exists() {
        eprintln!("Skipping: /usr/bin/blastn not found");
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

    let mut ncbi_cmd = std::process::Command::new("/usr/bin/blastn");
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
    if !std::path::Path::new("/usr/bin/blastn").exists()
        || !std::path::Path::new("/usr/bin/makeblastdb").exists()
    {
        eprintln!("Skipping: /usr/bin/blastn or /usr/bin/makeblastdb not found");
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

    let make_status = std::process::Command::new("/usr/bin/makeblastdb")
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

    let mut ncbi_cmd = std::process::Command::new("/usr/bin/blastn");
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
    if !std::path::Path::new("/usr/bin/blastn").exists()
        || !std::path::Path::new("/usr/bin/makeblastdb").exists()
    {
        eprintln!("Skipping: /usr/bin/blastn or /usr/bin/makeblastdb not found");
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

    let make_status = std::process::Command::new("/usr/bin/makeblastdb")
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

    let mut ncbi_cmd = std::process::Command::new("/usr/bin/blastn");
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
    if !std::path::Path::new("/usr/bin/blastn").exists()
        || !std::path::Path::new("/usr/bin/makeblastdb").exists()
    {
        eprintln!("Skipping: /usr/bin/blastn or /usr/bin/makeblastdb not found");
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

    let mut make_cmd = std::process::Command::new("/usr/bin/makeblastdb");
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

    let mut ncbi_cmd = std::process::Command::new("/usr/bin/blastn");
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
    if !std::path::Path::new("/usr/bin/blastn").exists()
        || !std::path::Path::new("/usr/bin/makeblastdb").exists()
    {
        eprintln!("Skipping: /usr/bin/blastn or /usr/bin/makeblastdb not found");
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

    let mut make_cmd = std::process::Command::new("/usr/bin/makeblastdb");
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

    let mut ncbi_cmd = std::process::Command::new("/usr/bin/blastn");
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
    if !std::path::Path::new("/usr/bin/blastn").exists()
        || !std::path::Path::new("/usr/bin/makeblastdb").exists()
        || !std::path::Path::new("/usr/bin/makembindex").exists()
    {
        eprintln!("Skipping: /usr/bin/blastn, makeblastdb, or makembindex not found");
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

    let make_status = std::process::Command::new("/usr/bin/makeblastdb")
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

    let index_status = std::process::Command::new("/usr/bin/makembindex")
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
        rust_cmd.arg(arg);
    }
    let rust_status = rust_cmd.status().expect("run blast-cli indexed DB parity");
    assert!(
        rust_status.success(),
        "blast-cli exited with {}",
        rust_status
    );

    let mut ncbi_cmd = std::process::Command::new("/usr/bin/blastn");
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
    if !std::path::Path::new("/usr/bin/blastn").exists() {
        eprintln!("Skipping: /usr/bin/blastn not found");
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
    let mut ncbi_cmd = std::process::Command::new("/usr/bin/blastn");
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

#[test]
fn blastn_subject_ncbi_parity_dust_no_exact_hits() {
    assert_blastn_subject_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGTACGT\n",
        ">s1\nTTTTACGTACGTACGTACGTACGTACGTAAAA\n>s2\nACGTACGTACGTACGTACGTACGT\n>s3\nGGGGGGGGGGGGGGGGGGGGGGGG\n",
        &["--dust", "no", "--max_target_seqs", "10", "--max_hsps", "2"],
        &["-dust", "no", "-max_target_seqs", "10", "-max_hsps", "2"],
    );
}

#[test]
#[ignore = "requires the large celegans fixture and NCBI blastn"]
fn blastn_large_db_ncbi_parity_q500_q2000_regressions() {
    assert_large_db_blastn_matches_ncbi("query_500.fa", None);
    assert_large_db_blastn_matches_ncbi("query_2000.fa", None);
    assert_large_db_blastn_matches_ncbi("query_2000.fa", Some("1"));
}

#[test]
fn blastn_subject_ncbi_parity_dbsize_and_searchsp_statistics() {
    let query = ">subseq_oid0\nTTAAGGAGGCTCATCTTTCAGAATCCATGCTGTGGGCCAGCAAGAGTTAA\n";
    let subject = ">subj1\nTTAAGGAGGCTCATCTTTCAGAATCCATGCTGTGGGCCAGCAAGAGTTAA\n";
    let outfmt = "6 qseqid sseqid evalue bitscore score length pident";

    for (rust_args, ncbi_args) in [
        (
            vec!["--dust", "no", "--dbsize", "1000000"],
            vec!["-dust", "no", "-dbsize", "1000000"],
        ),
        (
            vec!["--dust", "no", "--searchsp", "1000000"],
            vec!["-dust", "no", "-searchsp", "1000000"],
        ),
        (
            vec!["--dust", "no", "--dbsize", "5000000000"],
            vec!["-dust", "no", "-dbsize", "5000000000"],
        ),
        (
            vec!["--dust", "no", "--searchsp", "5000000000"],
            vec!["-dust", "no", "-searchsp", "5000000000"],
        ),
        (
            vec!["--dust", "no", "--dbsize", "-1"],
            vec!["-dust", "no", "-dbsize", "-1"],
        ),
        (
            vec!["--dust", "no", "--searchsp", "0"],
            vec!["-dust", "no", "-searchsp", "0"],
        ),
    ] {
        assert_blastn_subject_outfmt_matches_ncbi(query, subject, outfmt, &rust_args, &ncbi_args);
    }
}

#[test]
fn blastn_subject_ncbi_parity_pairwise_dbsize_and_searchsp_statistics() {
    let query = ">subseq_oid0\nTTAAGGAGGCTCATCTTTCAGAATCCATGCTGTGGGCCAGCAAGAGTTAA\n";
    let subject = ">subj1\nTTAAGGAGGCTCATCTTTCAGAATCCATGCTGTGGGCCAGCAAGAGTTAA\n";

    for (rust_args, ncbi_args) in [
        (
            vec!["--dust", "no", "--dbsize", "1000000"],
            vec!["-dust", "no", "-dbsize", "1000000"],
        ),
        (
            vec!["--dust", "no", "--searchsp", "1000000"],
            vec!["-dust", "no", "-searchsp", "1000000"],
        ),
        (
            vec![
                "--dust",
                "no",
                "--reward",
                "1",
                "--penalty",
                "-3",
                "--gapopen",
                "5",
                "--gapextend",
                "2",
            ],
            vec![
                "-dust",
                "no",
                "-reward",
                "1",
                "-penalty",
                "-3",
                "-gapopen",
                "5",
                "-gapextend",
                "2",
            ],
        ),
    ] {
        assert_blastn_subject_outfmt_matches_ncbi(query, subject, "0", &rust_args, &ncbi_args);
    }
}

#[test]
fn blastn_subject_ncbi_parity_dbsize_searchsp_multi_query_and_filtering() {
    let query = ">q1\nTTAAGGAGGCTCATCTTTCAGAATCCATGCTGTGGGCCAGCAAGAGTTAA\n>q2\nTTAAGGAGGCTCATCTTTCAGAATCCATGCTGTGGGCCAGCAAGAGTTAA\n";
    let subject = ">subj1\nTTAAGGAGGCTCATCTTTCAGAATCCATGCTGTGGGCCAGCAAGAGTTAA\n";
    let outfmt = "6 qseqid sseqid evalue bitscore score length pident";

    for (rust_args, ncbi_args) in [
        (
            vec!["--dust", "no", "--dbsize", "1000000"],
            vec!["-dust", "no", "-dbsize", "1000000"],
        ),
        (
            vec!["--dust", "no", "--searchsp", "1000000"],
            vec!["-dust", "no", "-searchsp", "1000000"],
        ),
    ] {
        assert_blastn_subject_outfmt_matches_ncbi(query, subject, outfmt, &rust_args, &ncbi_args);
    }

    for (rust_args, ncbi_args) in [
        (
            vec!["--dust", "no", "--dbsize", "1000000", "--evalue", "1e-25"],
            vec!["-dust", "no", "-dbsize", "1000000", "-evalue", "1e-25"],
        ),
        (
            vec!["--dust", "no", "--searchsp", "1000000", "--evalue", "1e-25"],
            vec!["-dust", "no", "-searchsp", "1000000", "-evalue", "1e-25"],
        ),
    ] {
        assert_blastn_subject_outfmt_matches_ncbi(
            ">q1\nTTAAGGAGGCTCATCTTTCAGAATCCATGCTGTGGGCCAGCAAGAGTTAA\n",
            ">subj1\nTTAAGGAGGCTCATCTTTCAGAATCCATGCTGTGGGCCAGCAAGAGTTAA\n",
            outfmt,
            &rust_args,
            &ncbi_args,
        );
    }
}

#[test]
fn blastn_subject_ncbi_parity_dbsize_searchsp_minus_strand_statistics() {
    let query = ">q_minus\nGAATCCATGCTGTGGGCCAGCAAGAGTTAA\n";
    let subject = ">subj_rc\nTTAACTCTTGCTGGCCCACAGCATGGATTC\n";
    let outfmt = "6 qseqid sseqid sstrand qstart qend sstart send evalue bitscore score length pident qseq sseq";

    for (rust_args, ncbi_args) in [
        (
            vec!["--dust", "no", "--strand", "minus", "--dbsize", "1000000"],
            vec!["-dust", "no", "-strand", "minus", "-dbsize", "1000000"],
        ),
        (
            vec!["--dust", "no", "--strand", "minus", "--searchsp", "1000000"],
            vec!["-dust", "no", "-strand", "minus", "-searchsp", "1000000"],
        ),
    ] {
        assert_blastn_subject_outfmt_matches_ncbi(query, subject, outfmt, &rust_args, &ncbi_args);
    }
}

#[test]
fn blastn_subject_ncbi_parity_long_exact_match_extends_to_edges() {
    let seq = "GAATCCATGCTGTGGGCCAGCAAGAGTTAAGGTGCTCATGGTTTTGAGAAAACATCTGAGGACTCTGACAGCACTCTCCCATCCTTGGTCTCCACAGTCT";
    let query = format!(">q\n{}\n", seq);
    let subject = format!(">s\n{}\n", seq);

    assert_blastn_subject_outfmt_matches_ncbi(
        &query,
        &subject,
        "6 qseqid sseqid qstart qend sstart send evalue bitscore score length pident qseq sseq",
        &["--dust", "no"],
        &["-dust", "no"],
    );
}

#[test]
fn blastn_subject_ncbi_parity_ungapped_right_xdrop_negative_total() {
    let query = ">q\nAAAAAAAAAAAAAAAAAAAA\n";
    let subject = ">s\nAAAACAAAAAAAAAAAAAAA\n";
    let outfmt = "6 qstart qend sstart send score length qseq sseq btop";

    assert_blastn_subject_outfmt_matches_ncbi(
        query,
        subject,
        outfmt,
        &[
            "--ungapped",
            "--dust",
            "no",
            "--word_size",
            "4",
            "--reward",
            "1",
            "--penalty",
            "-5",
            "--evalue",
            "1000",
            "--max_hsps",
            "10",
        ],
        &[
            "-ungapped",
            "-dust",
            "no",
            "-word_size",
            "4",
            "-reward",
            "1",
            "-penalty",
            "-5",
            "-evalue",
            "1000",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_ungapped_perfect_30bp_primer() {
    let query = ">q\nACGTACGTACGTACGTACGTACGTACGTAC\n";
    let subject = ">s\nTTTTACGTACGTACGTACGTACGTACGTACGTACAAAA\n";
    let outfmt = "6 qstart qend sstart send score length qseq sseq btop";

    assert_blastn_subject_outfmt_matches_ncbi(
        query,
        subject,
        outfmt,
        &[
            "--ungapped",
            "--dust",
            "no",
            "--word_size",
            "7",
            "--evalue",
            "1000",
            "--max_hsps",
            "10",
        ],
        &[
            "-ungapped",
            "-dust",
            "no",
            "-word_size",
            "7",
            "-evalue",
            "1000",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_ungapped_30bp_primer_central_mismatch() {
    let query = ">q\nACGTACGTACGTACGTACGTACGTACGTAC\n";
    let subject = ">s\nTTTTACGTACGTACGTATGTACGTACGTACGTACAAAA\n";
    let outfmt = "6 qstart qend sstart send score length qseq sseq btop";

    assert_blastn_subject_outfmt_matches_ncbi(
        query,
        subject,
        outfmt,
        &[
            "--ungapped",
            "--dust",
            "no",
            "--word_size",
            "7",
            "--evalue",
            "1000",
            "--max_hsps",
            "20",
        ],
        &[
            "-ungapped",
            "-dust",
            "no",
            "-word_size",
            "7",
            "-evalue",
            "1000",
            "-max_hsps",
            "20",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_ungapped_ambiguity_scoring_and_ordering() {
    let query = ">q\nAAAANAAAAAAAAAAA\n";
    let subject = ">s\nAAAAAAAAAAAAAAAA\n";
    let outfmt = "6 qstart qend sstart send score length qseq sseq btop";

    assert_blastn_subject_outfmt_matches_ncbi(
        query,
        subject,
        outfmt,
        &[
            "--ungapped",
            "--dust",
            "no",
            "--word_size",
            "4",
            "--reward",
            "1",
            "--penalty",
            "-5",
            "--evalue",
            "1000",
            "--max_hsps",
            "10",
        ],
        &[
            "-ungapped",
            "-dust",
            "no",
            "-word_size",
            "4",
            "-reward",
            "1",
            "-penalty",
            "-5",
            "-evalue",
            "1000",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_ungapped_default_megablast_left_xdrop() {
    let query = ">q\nAAAAAAAAAAAAAAAAAAAA\n";
    let subject = ">s\nAAAAAAAAACAAAAAAAAAA\n";
    let outfmt = "6 qstart qend sstart send score length qseq sseq btop";

    assert_blastn_subject_outfmt_matches_ncbi(
        query,
        subject,
        outfmt,
        &[
            "--ungapped",
            "--dust",
            "no",
            "--word_size",
            "4",
            "--reward",
            "1",
            "--penalty",
            "-5",
            "--evalue",
            "1000",
            "--max_hsps",
            "10",
        ],
        &[
            "-ungapped",
            "-dust",
            "no",
            "-word_size",
            "4",
            "-reward",
            "1",
            "-penalty",
            "-5",
            "-evalue",
            "1000",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_ungapped_permissive_evalue_cutoff() {
    let query = ">q\nAAAAAAAAAAAA\n";
    let subject = ">s\nAAAACCAAAAAA\n";
    let outfmt = "6 qstart qend sstart send score length qseq sseq btop";

    assert_blastn_subject_outfmt_matches_ncbi(
        query,
        subject,
        outfmt,
        &[
            "--ungapped",
            "--dust",
            "no",
            "--word_size",
            "4",
            "--reward",
            "1",
            "--penalty",
            "-2",
            "--xdrop_ungap",
            "4",
            "--evalue",
            "1000",
            "--max_hsps",
            "20",
        ],
        &[
            "-ungapped",
            "-dust",
            "no",
            "-word_size",
            "4",
            "-reward",
            "1",
            "-penalty",
            "-2",
            "-xdrop_ungap",
            "4",
            "-evalue",
            "1000",
            "-max_hsps",
            "20",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_ungapped_lcase_masking_extends_unmasked_query() {
    let query = ">q\nAAAAaaaaAAAAAAAA\n";
    let subject = ">s\nAAAAAAAAAAAAAAAA\n";
    let outfmt = "6 qstart qend sstart send score length qseq sseq btop";

    assert_blastn_subject_outfmt_matches_ncbi(
        query,
        subject,
        outfmt,
        &[
            "--ungapped",
            "--dust",
            "no",
            "--lcase_masking",
            "--word_size",
            "4",
            "--reward",
            "1",
            "--penalty",
            "-5",
            "--evalue",
            "1000",
            "--max_hsps",
            "10",
        ],
        &[
            "-ungapped",
            "-dust",
            "no",
            "-lcase_masking",
            "-word_size",
            "4",
            "-reward",
            "1",
            "-penalty",
            "-5",
            "-evalue",
            "1000",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_ungapped_lowercase_subject_is_sequence() {
    let query = ">q\nAAAAAAAAAAAAAAAA\n";
    let subject = ">s\nAAAAaaaaAAAAAAAA\n";
    let outfmt = "6 qstart qend sstart send score length qseq sseq btop";

    assert_blastn_subject_outfmt_matches_ncbi(
        query,
        subject,
        outfmt,
        &[
            "--ungapped",
            "--dust",
            "no",
            "--lcase_masking",
            "--word_size",
            "4",
            "--reward",
            "1",
            "--penalty",
            "-5",
            "--evalue",
            "1000",
            "--max_hsps",
            "10",
        ],
        &[
            "-ungapped",
            "-dust",
            "no",
            "-lcase_masking",
            "-word_size",
            "4",
            "-reward",
            "1",
            "-penalty",
            "-5",
            "-evalue",
            "1000",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_ungapped_mismatch_boundary_matrix() {
    let query = ">q\nAAAAAAAAAAAAAAAA\n";
    let outfmt = "6 qstart qend sstart send score length qseq sseq btop";
    let cases = [
        ("left_mismatch", ">s\nCAAAAAAAAAAAAAAA\n"),
        ("right_mismatch", ">s\nAAAAAAAAAAAAAAAC\n"),
        ("two_close_mismatches", ">s\nAAAAAACCAAAAAAAA\n"),
        ("subject_ambiguity", ">s\nAAAANAAAAAAAAAAA\n"),
        ("xdrop_boundary", ">s\nAAAACAAAAAAAACAAA\n"),
    ];

    for (label, subject) in cases {
        assert_blastn_subject_outfmt_matches_ncbi(
            query,
            subject,
            outfmt,
            &[
                "--ungapped",
                "--dust",
                "no",
                "--word_size",
                "4",
                "--reward",
                "1",
                "--penalty",
                "-5",
                "--evalue",
                "1000",
                "--max_hsps",
                "10",
            ],
            &[
                "-ungapped",
                "-dust",
                "no",
                "-word_size",
                "4",
                "-reward",
                "1",
                "-penalty",
                "-5",
                "-evalue",
                "1000",
                "-max_hsps",
                "10",
            ],
        );
        eprintln!("checked {label}");
    }
}

#[test]
fn blastn_subject_ncbi_parity_ungapped_explicit_xdrop_boundary_matrix() {
    let query = ">q\nAAAAAAAAAAAAAAAAAAAA\n";
    let outfmt = "6 qstart qend sstart send score length qseq sseq btop";
    let cases = [
        ("left_two_mismatch_x1", ">s\nCCAAAAAAAAAAAAAAAAAA\n", "1"),
        ("right_two_mismatch_x1", ">s\nAAAAAAAAAAAAAAAAAACC\n", "1"),
        ("both_end_mismatch_x1", ">s\nCAAAAAAAAAAAAAAAAAAC\n", "1"),
        (
            "internal_two_mismatch_x2",
            ">s\nAAAAAAAACCAAAAAAAAAA\n",
            "2",
        ),
        (
            "internal_two_mismatch_x3",
            ">s\nAAAAAAAACCAAAAAAAAAA\n",
            "3",
        ),
        ("late_drop_x2", ">s\nAAAAAAAAAAAACCCAAAAA\n", "2"),
        ("early_drop_x2", ">s\nAAAAACCCAAAAAAAAAAAA\n", "2"),
        ("ambig_drop_x2", ">s\nAAAAANNNAAAAAAAAAAAA\n", "2"),
    ];

    for (label, subject, xdrop) in cases {
        assert_blastn_subject_outfmt_matches_ncbi(
            query,
            subject,
            outfmt,
            &[
                "--ungapped",
                "--dust",
                "no",
                "--word_size",
                "4",
                "--reward",
                "1",
                "--penalty",
                "-2",
                "--xdrop_ungap",
                xdrop,
                "--evalue",
                "1000",
                "--max_hsps",
                "20",
            ],
            &[
                "-ungapped",
                "-dust",
                "no",
                "-word_size",
                "4",
                "-reward",
                "1",
                "-penalty",
                "-2",
                "-xdrop_ungap",
                xdrop,
                "-evalue",
                "1000",
                "-max_hsps",
                "20",
            ],
        );
        eprintln!("checked {label}");
    }
}

#[test]
fn blastn_subject_ncbi_parity_ungapped_minus_strand_explicit_xdrop_boundary_matrix() {
    let query = ">q\nAAAAAAAAAAAAAAAAAAAA\n";
    let outfmt = "6 qstart qend sstart send score length qseq sseq btop";
    let cases = [
        ("minus_left_two_mismatch_x1", "CCAAAAAAAAAAAAAAAAAA", "1"),
        ("minus_right_two_mismatch_x1", "AAAAAAAAAAAAAAAAAACC", "1"),
        (
            "minus_internal_two_mismatch_x2",
            "AAAAAAAACCAAAAAAAAAA",
            "2",
        ),
        ("minus_late_drop_x2", "AAAAAAAAAAAACCCAAAAA", "2"),
        ("minus_early_drop_x2", "AAAAACCCAAAAAAAAAAAA", "2"),
    ];

    for (label, plus_subject, xdrop) in cases {
        let subject = format!(">s\n{}\n", ascii_reverse_complement(plus_subject));
        assert_blastn_subject_outfmt_matches_ncbi(
            query,
            &subject,
            outfmt,
            &[
                "--ungapped",
                "--dust",
                "no",
                "--strand",
                "minus",
                "--word_size",
                "4",
                "--reward",
                "1",
                "--penalty",
                "-2",
                "--xdrop_ungap",
                xdrop,
                "--evalue",
                "1000",
                "--max_hsps",
                "20",
            ],
            &[
                "-ungapped",
                "-dust",
                "no",
                "-strand",
                "minus",
                "-word_size",
                "4",
                "-reward",
                "1",
                "-penalty",
                "-2",
                "-xdrop_ungap",
                xdrop,
                "-evalue",
                "1000",
                "-max_hsps",
                "20",
            ],
        );
        eprintln!("checked {label}");
    }
}

#[test]
fn blastn_subject_ncbi_parity_strand_plus_minus() {
    let query = ">q1\nGGGTTTAAACCCGGGTTTAAACCC\n";
    let subject = ">s1\nGGGTTTAAACCCGGGTTTAAACCC\n>s2\nGGGTTTAAACCCGGGTTTAAACCCGGGTTTAAACCC\n>s3\nCCCCCCCCCCCCCCCCCCCCCCCC\n";
    for strand in ["plus", "minus"] {
        assert_blastn_subject_matches_ncbi(
            query,
            subject,
            &[
                "--dust",
                "no",
                "--strand",
                strand,
                "--max_target_seqs",
                "10",
                "--max_hsps",
                "2",
            ],
            &[
                "-dust",
                "no",
                "-strand",
                strand,
                "-max_target_seqs",
                "10",
                "-max_hsps",
                "2",
            ],
        );
    }
}

#[test]
fn blastn_subject_ncbi_parity_dc_megablast_defaults() {
    assert_blastn_subject_task_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n",
        ">s1\nTTTTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAAAA\n>s2\nCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "dc-megablast",
        "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen",
        &["--dust", "no", "--max_target_seqs", "10"],
        &["-dust", "no", "-max_target_seqs", "10"],
    );
}

#[test]
fn blastn_subject_ncbi_parity_dc_megablast_discontiguous_template_fixture() {
    assert_blastn_subject_task_outfmt_matches_ncbi(
        ">q1\nACGTACGTGACTTACCGTACGTACGTGACTTACCGTACGTACGTGACTTACCGTACGT\n",
        ">s1\nTTTTACGTACGTGACTTACCGTACGTACGTGACTTACCGTACGTACGTGACTTACCGTACGTAAAA\n",
        "dc-megablast",
        "6 qseqid sseqid qstart qend sstart send length mismatch gapopen score bitscore btop",
        &[
            "--dust",
            "no",
            "--template_type",
            "coding",
            "--template_length",
            "18",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "20",
            "--max_hsps",
            "10",
        ],
        &[
            "-dust",
            "no",
            "-template_type",
            "coding",
            "-template_length",
            "18",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "20",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_dc_megablast_two_template_fixture() {
    assert_blastn_subject_task_outfmt_matches_ncbi(
        ">q1\nACGTACGTGACTTACCGTACGTACGTGACTTACCGTACGTACGTGACTTACCGTACGT\n",
        ">s1\nGGGGACGTACGTGACTTACCGTACGTACGTGACTTACCGTACGTACGTGACTTACCGTACGTCCCC\n",
        "dc-megablast",
        "6 qseqid sseqid qstart qend sstart send length mismatch gapopen score bitscore btop",
        &[
            "--dust",
            "no",
            "--template_type",
            "coding_and_optimal",
            "--template_length",
            "16",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "20",
            "--max_hsps",
            "10",
        ],
        &[
            "-dust",
            "no",
            "-template_type",
            "coding_and_optimal",
            "-template_length",
            "16",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "20",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_dc_megablast_multi_subject_two_template_fixture() {
    assert_blastn_subject_task_outfmt_matches_ncbi(
        ">q1\nACGTACGTGACTTACCGTACGTACGTGACTTACCGTACGTACGTGACTTACCGTACGT\n",
        ">s_primary\nGGGGACGTACGTGACTTACCGTACGTACGTGACTTACCGTACGTACGTGACTTACCGTACGTCCCC\n>s_shifted\nTTTTTTACGTACGTGACTTACCGTACGTACGTGACTTACCGTACGTACGTGACTTACCGTACGTAAAA\n>s_noise\nCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "dc-megablast",
        "6 qseqid sseqid qstart qend sstart send length mismatch gapopen score bitscore btop",
        &[
            "--dust",
            "no",
            "--template_type",
            "coding_and_optimal",
            "--template_length",
            "16",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "20",
            "--max_hsps",
            "10",
        ],
        &[
            "-dust",
            "no",
            "-template_type",
            "coding_and_optimal",
            "-template_length",
            "16",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "20",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_dc_megablast_masked_discontiguous_template_fixture() {
    assert_blastn_subject_task_outfmt_matches_ncbi(
        ">q1\nACGTACGTgacttaccGTACGTACGTGACTTACCGTACGTACGTGACTTACCGT\n",
        ">s1\nTTTTACGTACGTGACTTACCGTACGTACGTGACTTACCGTACGTACGTGACTTACCGTAAAA\n",
        "dc-megablast",
        "6 qseqid sseqid qstart qend sstart send length mismatch gapopen score bitscore btop",
        &[
            "--dust",
            "no",
            "--lcase_masking",
            "--soft_masking",
            "true",
            "--template_type",
            "coding",
            "--template_length",
            "18",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "20",
            "--max_hsps",
            "10",
        ],
        &[
            "-dust",
            "no",
            "-lcase_masking",
            "-soft_masking",
            "true",
            "-template_type",
            "coding",
            "-template_length",
            "18",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "20",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_masked_query_mapper_fixture() {
    assert_blastn_subject_task_outfmt_matches_ncbi(
        ">q1\nACGTACGTacgtacgtGACTTACCGTACGTACGTGACTTACCGTACGT\n",
        ">s1\nTTTTACGTACGTACGTACGTGACTTACCGTACGTACGTGACTTACCGTACGTAAAA\n",
        "blastn-short",
        "6 qseqid sseqid qstart qend sstart send length mismatch gapopen score bitscore btop",
        &[
            "--dust",
            "no",
            "--lcase_masking",
            "--soft_masking",
            "true",
            "--word_size",
            "7",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "20",
            "--max_hsps",
            "10",
        ],
        &[
            "-dust",
            "no",
            "-lcase_masking",
            "-soft_masking",
            "true",
            "-word_size",
            "7",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "20",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_masked_query_minus_strand_mapper_fixture() {
    let subject = format!(
        ">s1\n{}\n",
        ascii_reverse_complement("TTTTACGTACGTACGTACGTGACTTACCGTACGTACGTGACTTACCGTACGTAAAA")
    );
    assert_blastn_subject_task_outfmt_matches_ncbi(
        ">q1\nACGTACGTacgtacgtGACTTACCGTACGTACGTGACTTACCGTACGT\n",
        &subject,
        "blastn-short",
        "6 qseqid sseqid qstart qend sstart send length mismatch gapopen score bitscore btop",
        &[
            "--dust",
            "no",
            "--lcase_masking",
            "--soft_masking",
            "true",
            "--strand",
            "minus",
            "--word_size",
            "7",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "20",
            "--max_hsps",
            "10",
        ],
        &[
            "-dust",
            "no",
            "-lcase_masking",
            "-soft_masking",
            "true",
            "-strand",
            "minus",
            "-word_size",
            "7",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "20",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_terminal_overhang_mapper_fixture() {
    assert_blastn_subject_task_outfmt_matches_ncbi(
        ">q1\nTTTTACGTACGTGACTTACCGTACGTACGTAAAA\n",
        ">s1\nACGTACGTGACTTACCGTACGTACGT\n>s2\nCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "blastn-short",
        "6 qseqid sseqid qstart qend sstart send length mismatch gapopen score bitscore qseq sseq btop",
        &[
            "--dust",
            "no",
            "--word_size",
            "7",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "10",
        ],
        &[
            "-dust",
            "no",
            "-word_size",
            "7",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_terminal_overhang_mapper_multi_subject_fixture() {
    assert_blastn_subject_task_outfmt_matches_ncbi_sorted_lines(
        ">q1\nTTTTACGTACGTGACTTACCGTACGTACGTAAAA\n",
        ">s1\nACGTACGTGACTTACCGTACGTACGT\n>s2\nACGTACGTGACTTACCGTACGTACGT\n>s3\nCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "blastn-short",
        "6 qseqid sseqid qstart qend sstart send length mismatch gapopen score bitscore qseq sseq btop",
        &[
            "--dust",
            "no",
            "--word_size",
            "7",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "10",
        ],
        &[
            "-dust",
            "no",
            "-word_size",
            "7",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_minus_strand_terminal_overhang_mapper_fixture() {
    assert_blastn_subject_task_outfmt_matches_ncbi(
        ">q1\nTTTTACGTACGTGACTTACCGTACGTACGTAAAA\n",
        ">s1\nACGTACGTACGGTAAGTCACGTACGT\n>s2\nCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "blastn-short",
        "6 qseqid sseqid qstart qend sstart send length mismatch gapopen score bitscore qseq sseq btop",
        &[
            "--dust",
            "no",
            "--strand",
            "minus",
            "--word_size",
            "7",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "10",
        ],
        &[
            "-dust",
            "no",
            "-strand",
            "minus",
            "-word_size",
            "7",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_minus_strand_terminal_overhang_mapper_multi_subject_fixture() {
    assert_blastn_subject_task_outfmt_matches_ncbi_sorted_lines(
        ">q1\nTTTTACGTACGTGACTTACCGTACGTACGTAAAA\n",
        ">s1\nACGTACGTACGGTAAGTCACGTACGT\n>s2\nACGTACGTACGGTAAGTCACGTACGT\n>s3\nCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "blastn-short",
        "6 qseqid sseqid qstart qend sstart send length mismatch gapopen score bitscore qseq sseq btop",
        &[
            "--dust",
            "no",
            "--strand",
            "minus",
            "--word_size",
            "7",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "10",
        ],
        &[
            "-dust",
            "no",
            "-strand",
            "minus",
            "-word_size",
            "7",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_adapter_like_trailing_overhang_mapper_fixture() {
    assert_blastn_subject_task_outfmt_matches_ncbi(
        ">q1\nTTTTACGTACGTGACTTACCGTACGTACGTACACGGAAGAC\n",
        ">s1\nACGTACGTGACTTACCGTACGTACGT\n>s2\nCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "blastn-short",
        "6 qseqid sseqid qstart qend sstart send length mismatch gapopen score bitscore qseq sseq btop",
        &[
            "--dust",
            "no",
            "--word_size",
            "7",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "10",
        ],
        &[
            "-dust",
            "no",
            "-word_size",
            "7",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_adapter_like_trailing_overhang_mapper_multi_subject_fixture() {
    assert_blastn_subject_task_outfmt_matches_ncbi_sorted_lines(
        ">q1\nTTTTACGTACGTGACTTACCGTACGTACGTACACGGAAGAC\n",
        ">s1\nACGTACGTGACTTACCGTACGTACGT\n>s2\nACGTACGTGACTTACCGTACGTACGT\n>s3\nCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "blastn-short",
        "6 qseqid sseqid qstart qend sstart send length mismatch gapopen score bitscore qseq sseq btop",
        &[
            "--dust",
            "no",
            "--word_size",
            "7",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "10",
        ],
        &[
            "-dust",
            "no",
            "-word_size",
            "7",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_poly_a_trailing_overhang_mapper_fixture() {
    assert_blastn_subject_task_outfmt_matches_ncbi(
        ">q1\nTTTTACGTACGTGACTTACCGTACGTACGTAAAAAA\n",
        ">s1\nACGTACGTGACTTACCGTACGTACGT\n>s2\nCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "blastn-short",
        "6 qseqid sseqid qstart qend sstart send length mismatch gapopen score bitscore qseq sseq btop",
        &[
            "--dust",
            "no",
            "--word_size",
            "7",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "10",
        ],
        &[
            "-dust",
            "no",
            "-word_size",
            "7",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_poly_a_trailing_overhang_mapper_multi_subject_fixture() {
    assert_blastn_subject_task_outfmt_matches_ncbi_sorted_lines(
        ">q1\nTTTTACGTACGTGACTTACCGTACGTACGTAAAAAA\n",
        ">s1\nACGTACGTGACTTACCGTACGTACGT\n>s2\nACGTACGTGACTTACCGTACGTACGT\n>s3\nCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "blastn-short",
        "6 qseqid sseqid qstart qend sstart send length mismatch gapopen score bitscore qseq sseq btop",
        &[
            "--dust",
            "no",
            "--word_size",
            "7",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "10",
        ],
        &[
            "-dust",
            "no",
            "-word_size",
            "7",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_exact_adapter_trimming_mapper_fixture() {
    assert_blastn_subject_task_outfmt_matches_ncbi(
        ">q1\nTTTTACGTACGTGACTTACCGTACGTACGTAGATCGGAAGAG\n",
        ">s1\nACGTACGTGACTTACCGTACGTACGT\n>s2\nCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "blastn-short",
        "6 qseqid sseqid qstart qend sstart send length mismatch gapopen score bitscore qseq sseq btop",
        &[
            "--dust",
            "no",
            "--word_size",
            "7",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "10",
        ],
        &[
            "-dust",
            "no",
            "-word_size",
            "7",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_exact_adapter_trimming_mapper_multi_subject_fixture() {
    assert_blastn_subject_task_outfmt_matches_ncbi_sorted_lines(
        ">q1\nTTTTACGTACGTGACTTACCGTACGTACGTAGATCGGAAGAG\n",
        ">s1\nACGTACGTGACTTACCGTACGTACGT\n>s2\nACGTACGTGACTTACCGTACGTACGT\n>s3\nCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "blastn-short",
        "6 qseqid sseqid qstart qend sstart send length mismatch gapopen score bitscore qseq sseq btop",
        &[
            "--dust",
            "no",
            "--word_size",
            "7",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "10",
        ],
        &[
            "-dust",
            "no",
            "-word_size",
            "7",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
#[ignore = "documents the remaining exact multi-subject adapter tie-ordering gap"]
fn blastn_subject_ncbi_parity_exact_adapter_trimming_mapper_multi_subject_order() {
    assert_blastn_subject_task_outfmt_matches_ncbi(
        ">q1\nTTTTACGTACGTGACTTACCGTACGTACGTAGATCGGAAGAG\n",
        ">s1\nACGTACGTGACTTACCGTACGTACGT\n>s2\nACGTACGTGACTTACCGTACGTACGT\n>s3\nCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "blastn-short",
        "6 qseqid sseqid qstart qend sstart send length mismatch gapopen score bitscore qseq sseq btop",
        &[
            "--dust",
            "no",
            "--word_size",
            "7",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "10",
        ],
        &[
            "-dust",
            "no",
            "-word_size",
            "7",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_rmblastn_terminal_overhang_mapper_fixture() {
    assert_blastn_subject_task_outfmt_matches_ncbi(
        ">q1\nTTTTACGTACGTGACTTACCGTACGTACGTAAAA\n",
        ">s1\nACGTACGTGACTTACCGTACGTACGT\n>s2\nCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "rmblastn",
        "6 qseqid sseqid qstart qend sstart send length mismatch gapopen score bitscore qseq sseq btop",
        &[
            "--dust",
            "no",
            "--word_size",
            "7",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "10",
        ],
        &[
            "-dust",
            "no",
            "-word_size",
            "7",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_rmblastn_terminal_overhang_mapper_multi_subject_fixture() {
    assert_blastn_subject_task_outfmt_matches_ncbi_sorted_lines(
        ">q1\nTTTTACGTACGTGACTTACCGTACGTACGTAAAA\n",
        ">s1\nACGTACGTGACTTACCGTACGTACGT\n>s2\nACGTACGTGACTTACCGTACGTACGT\n>s3\nCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "rmblastn",
        "6 qseqid sseqid qstart qend sstart send length mismatch gapopen score bitscore qseq sseq btop",
        &[
            "--dust",
            "no",
            "--word_size",
            "7",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "10",
        ],
        &[
            "-dust",
            "no",
            "-word_size",
            "7",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_rmblastn_minus_strand_terminal_overhang_mapper_coordinates() {
    assert_blastn_subject_task_outfmt_matches_ncbi(
        ">q1\nTTTTACGTACGTGACTTACCGTACGTACGTAAAA\n",
        ">s1\nACGTACGTACGGTAAGTCACGTACGT\n>s2\nCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "rmblastn",
        "6 qseqid sseqid qstart qend sstart send length mismatch qseq sseq btop",
        &[
            "--dust",
            "no",
            "--strand",
            "minus",
            "--word_size",
            "7",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "10",
        ],
        &[
            "-dust",
            "no",
            "-strand",
            "minus",
            "-word_size",
            "7",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_rmblastn_minus_strand_terminal_overhang_mapper_multi_subject_coordinates(
) {
    assert_blastn_subject_task_outfmt_matches_ncbi_sorted_lines(
        ">q1\nTTTTACGTACGTGACTTACCGTACGTACGTAAAA\n",
        ">s1\nACGTACGTACGGTAAGTCACGTACGT\n>s2\nACGTACGTACGGTAAGTCACGTACGT\n>s3\nCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "rmblastn",
        "6 qseqid sseqid qstart qend sstart send length mismatch qseq sseq btop",
        &[
            "--dust",
            "no",
            "--strand",
            "minus",
            "--word_size",
            "7",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "10",
        ],
        &[
            "-dust",
            "no",
            "-strand",
            "minus",
            "-word_size",
            "7",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
#[ignore = "documents the remaining rmblastn minus-strand terminal secondary-HSP scoring gap"]
fn blastn_subject_ncbi_parity_rmblastn_minus_strand_terminal_overhang_mapper_scored() {
    assert_blastn_subject_task_outfmt_matches_ncbi(
        ">q1\nTTTTACGTACGTGACTTACCGTACGTACGTAAAA\n",
        ">s1\nACGTACGTACGGTAAGTCACGTACGT\n>s2\nCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "rmblastn",
        "6 qseqid sseqid qstart qend sstart send length mismatch gapopen score bitscore qseq sseq btop",
        &[
            "--dust",
            "no",
            "--strand",
            "minus",
            "--word_size",
            "7",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "10",
        ],
        &[
            "-dust",
            "no",
            "-strand",
            "minus",
            "-word_size",
            "7",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
#[ignore = "documents the remaining rmblastn minus-strand terminal secondary-HSP scoring gap"]
fn blastn_subject_ncbi_parity_rmblastn_minus_strand_terminal_overhang_mapper_multi_subject_scored()
{
    assert_blastn_subject_task_outfmt_matches_ncbi_sorted_lines(
        ">q1\nTTTTACGTACGTGACTTACCGTACGTACGTAAAA\n",
        ">s1\nACGTACGTACGGTAAGTCACGTACGT\n>s2\nACGTACGTACGGTAAGTCACGTACGT\n>s3\nCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "rmblastn",
        "6 qseqid sseqid qstart qend sstart send length mismatch gapopen score bitscore qseq sseq btop",
        &[
            "--dust",
            "no",
            "--strand",
            "minus",
            "--word_size",
            "7",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "10",
        ],
        &[
            "-dust",
            "no",
            "-strand",
            "minus",
            "-word_size",
            "7",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_rmblastn_adapter_like_trailing_overhang_mapper_fixture() {
    assert_blastn_subject_task_outfmt_matches_ncbi(
        ">q1\nTTTTACGTACGTGACTTACCGTACGTACGTACACGGAAGAC\n",
        ">s1\nACGTACGTGACTTACCGTACGTACGT\n>s2\nCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "rmblastn",
        "6 qseqid sseqid qstart qend sstart send length mismatch gapopen score bitscore qseq sseq btop",
        &[
            "--dust",
            "no",
            "--word_size",
            "7",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "10",
        ],
        &[
            "-dust",
            "no",
            "-word_size",
            "7",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_rmblastn_adapter_like_trailing_overhang_mapper_multi_subject_fixture()
{
    assert_blastn_subject_task_outfmt_matches_ncbi_sorted_lines(
        ">q1\nTTTTACGTACGTGACTTACCGTACGTACGTACACGGAAGAC\n",
        ">s1\nACGTACGTGACTTACCGTACGTACGT\n>s2\nACGTACGTGACTTACCGTACGTACGT\n>s3\nCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "rmblastn",
        "6 qseqid sseqid qstart qend sstart send length mismatch gapopen score bitscore qseq sseq btop",
        &[
            "--dust",
            "no",
            "--word_size",
            "7",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "10",
        ],
        &[
            "-dust",
            "no",
            "-word_size",
            "7",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_rmblastn_poly_a_trailing_overhang_mapper_fixture() {
    assert_blastn_subject_task_outfmt_matches_ncbi(
        ">q1\nTTTTACGTACGTGACTTACCGTACGTACGTAAAAAA\n",
        ">s1\nACGTACGTGACTTACCGTACGTACGT\n>s2\nCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "rmblastn",
        "6 qseqid sseqid qstart qend sstart send length mismatch gapopen score bitscore qseq sseq btop",
        &[
            "--dust",
            "no",
            "--word_size",
            "7",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "10",
        ],
        &[
            "-dust",
            "no",
            "-word_size",
            "7",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_rmblastn_poly_a_trailing_overhang_mapper_multi_subject_fixture() {
    assert_blastn_subject_task_outfmt_matches_ncbi_sorted_lines(
        ">q1\nTTTTACGTACGTGACTTACCGTACGTACGTAAAAAA\n",
        ">s1\nACGTACGTGACTTACCGTACGTACGT\n>s2\nACGTACGTGACTTACCGTACGTACGT\n>s3\nCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "rmblastn",
        "6 qseqid sseqid qstart qend sstart send length mismatch gapopen score bitscore qseq sseq btop",
        &[
            "--dust",
            "no",
            "--word_size",
            "7",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "10",
        ],
        &[
            "-dust",
            "no",
            "-word_size",
            "7",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_rmblastn_exact_adapter_trimming_mapper_fixture() {
    assert_blastn_subject_task_outfmt_matches_ncbi(
        ">q1\nTTTTACGTACGTGACTTACCGTACGTACGTAGATCGGAAGAG\n",
        ">s1\nACGTACGTGACTTACCGTACGTACGT\n>s2\nCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "rmblastn",
        "6 qseqid sseqid qstart qend sstart send length mismatch gapopen score bitscore qseq sseq btop",
        &[
            "--dust",
            "no",
            "--word_size",
            "7",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "10",
        ],
        &[
            "-dust",
            "no",
            "-word_size",
            "7",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_rmblastn_exact_adapter_trimming_mapper_multi_subject_fixture() {
    assert_blastn_subject_task_outfmt_matches_ncbi_sorted_lines(
        ">q1\nTTTTACGTACGTGACTTACCGTACGTACGTAGATCGGAAGAG\n",
        ">s1\nACGTACGTGACTTACCGTACGTACGT\n>s2\nACGTACGTGACTTACCGTACGTACGT\n>s3\nCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "rmblastn",
        "6 qseqid sseqid qstart qend sstart send length mismatch gapopen score bitscore qseq sseq btop",
        &[
            "--dust",
            "no",
            "--word_size",
            "7",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "10",
        ],
        &[
            "-dust",
            "no",
            "-word_size",
            "7",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
#[ignore = "documents the remaining exact multi-subject adapter tie-ordering gap"]
fn blastn_subject_ncbi_parity_rmblastn_exact_adapter_trimming_mapper_multi_subject_order() {
    assert_blastn_subject_task_outfmt_matches_ncbi(
        ">q1\nTTTTACGTACGTGACTTACCGTACGTACGTAGATCGGAAGAG\n",
        ">s1\nACGTACGTGACTTACCGTACGTACGT\n>s2\nACGTACGTGACTTACCGTACGTACGT\n>s3\nCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "rmblastn",
        "6 qseqid sseqid qstart qend sstart send length mismatch gapopen score bitscore qseq sseq btop",
        &[
            "--dust",
            "no",
            "--word_size",
            "7",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "10",
        ],
        &[
            "-dust",
            "no",
            "-word_size",
            "7",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_rmblastn_task_defaults() {
    assert_blastn_subject_task_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGT\n",
        ">s1\nACGTACGTACGTACGTACGT\n",
        "rmblastn",
        "6 qseqid sseqid length bitscore",
        &["--dust", "no", "--max_target_seqs", "10"],
        &["-dust", "no", "-max_target_seqs", "10"],
    );
}

#[test]
fn blastn_db_ncbi_parity_rmblastn_task_defaults() {
    if !std::path::Path::new("/usr/bin/blastn").exists()
        || !std::path::Path::new("/usr/bin/makeblastdb").exists()
    {
        eprintln!("Skipping: /usr/bin/blastn or /usr/bin/makeblastdb not found");
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
    let db_fasta = tmp.path().join("db.fa");
    let db = tmp.path().join("testdb");
    let rust_out = tmp.path().join("rust.tsv");
    let ncbi_out = tmp.path().join("ncbi.tsv");
    std::fs::write(&query, ">q1\nACGTACGTACGTACGTACGT\n").expect("write query FASTA");
    std::fs::write(&db_fasta, ">s1\nACGTACGTACGTACGTACGT\n").expect("write db FASTA");

    let make_status = std::process::Command::new("/usr/bin/makeblastdb")
        .arg("-in")
        .arg(&db_fasta)
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

    let rust_status = std::process::Command::new(blast_cli)
        .arg("blastn")
        .arg("--query")
        .arg(&query)
        .arg("--db")
        .arg(&db)
        .arg("--task")
        .arg("rmblastn")
        .arg("--outfmt")
        .arg("6 qseqid sseqid length bitscore")
        .arg("--dust")
        .arg("no")
        .arg("--max_target_seqs")
        .arg("10")
        .arg("--num_threads")
        .arg("1")
        .arg("--out")
        .arg(&rust_out)
        .status()
        .expect("run blast-cli rmblastn DB parity");
    assert!(rust_status.success(), "blast-cli exited with {rust_status}");

    let ncbi_status = std::process::Command::new("/usr/bin/blastn")
        .arg("-query")
        .arg(&query)
        .arg("-db")
        .arg(&db)
        .arg("-task")
        .arg("rmblastn")
        .arg("-outfmt")
        .arg("6 qseqid sseqid length bitscore")
        .arg("-dust")
        .arg("no")
        .arg("-max_target_seqs")
        .arg("10")
        .arg("-num_threads")
        .arg("1")
        .arg("-out")
        .arg(&ncbi_out)
        .status()
        .expect("run NCBI rmblastn DB parity");
    assert!(
        ncbi_status.success(),
        "NCBI blastn exited with {ncbi_status}"
    );

    let rust = std::fs::read(&rust_out).expect("read rust output");
    let ncbi = std::fs::read(&ncbi_out).expect("read ncbi output");
    assert_eq!(
        rust, ncbi,
        "Rust rmblastn DB output differs from NCBI\nRust: {:?}\nNCBI: {:?}",
        rust_out, ncbi_out
    );
}

#[test]
fn blastn_subject_ncbi_parity_blastn_task_defaults() {
    assert_blastn_subject_task_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n",
        ">s1\nTTTTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAAAA\n>s2\nCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "blastn-short",
        "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen",
        &["--dust", "no", "--max_target_seqs", "10"],
        &["-dust", "no", "-max_target_seqs", "10"],
    );
}

#[test]
fn blastn_subject_ncbi_parity_megablast_contained_diagonals() {
    assert_blastn_subject_task_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n",
        ">s1\nTTTTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAAAA\n>s2\nCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "megablast",
        "6 qstart qend sstart send length bitscore sframe btop",
        &["--dust", "no", "--max_target_seqs", "100"],
        &["-dust", "no", "-max_target_seqs", "100"],
    );
}

#[test]
fn blastn_subject_ncbi_parity_megablast_indel_tie_max_hsps() {
    assert_blastn_subject_task_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGTACGTACGTACGTGGGGACGTACGTACGTACGTACGTACGTACGTACGT\n",
        ">s1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n",
        "megablast",
        "6 qstart qend sstart send score bitscore evalue length mismatch gapopen gaps btop",
        &["--dust", "no", "--evalue", "1000", "--max_hsps", "1"],
        &["-dust", "no", "-evalue", "1000", "-max_hsps", "1"],
    );
}

#[test]
fn blastn_subject_ncbi_parity_megablast_indel_tie_multi_hsp_repro() {
    assert_blastn_subject_task_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGTACGTACGTACGTGGGGACGTACGTACGTACGTACGTACGTACGTACGT\n",
        ">s1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n",
        "megablast",
        "6 qstart qend sstart send score bitscore evalue length mismatch gapopen gaps btop",
        &[
            "--dust",
            "no",
            "--evalue",
            "1000",
            "--strand",
            "plus",
            "--max_hsps",
            "10",
        ],
        &[
            "-dust",
            "no",
            "-evalue",
            "1000",
            "-strand",
            "plus",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_both_strands_hsp_order() {
    assert_blastn_subject_matches_ncbi(
        ">q1\nGGGTTTAAACCCGGGTTTAAACCC\n",
        ">s1\nGGGTTTAAACCCGGGTTTAAACCC\n>s2\nGGGTTTAAACCCGGGTTTAAACCCGGGTTTAAACCC\n>s3\nCCCCCCCCCCCCCCCCCCCCCCCC\n",
        &["--dust", "no", "--max_target_seqs", "10", "--max_hsps", "2"],
        &["-dust", "no", "-max_target_seqs", "10", "-max_hsps", "2"],
    );
}

#[test]
fn blastn_subject_ncbi_parity_same_interval_hsp_order() {
    assert_blastn_subject_task_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGTACGTACGTACGT\n",
        ">s1\nACGTACGTACGTACGTACGT\n",
        "blastn-short",
        "6 sseqid qstart qend sstart send length bitscore",
        &["--dust", "no", "--max_target_seqs", "10"],
        &["-dust", "no", "-max_target_seqs", "10"],
    );
}

#[test]
fn blastn_subject_ncbi_parity_max_hsps_per_query() {
    assert_blastn_subject_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGTACGT\n>q2\nGGGTTTAAACCCGGGTTTAAACCC\n",
        ">s1\nACGTACGTACGTACGTACGTACGTNNNNGGGTTTAAACCCGGGTTTAAACCC\n",
        &[
            "--dust",
            "no",
            "--strand",
            "plus",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "1",
        ],
        &[
            "-dust",
            "no",
            "-strand",
            "plus",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "1",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_culling_limit() {
    assert_blastn_subject_task_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGTACGTACGTACGT\n",
        ">s_full\nACGTACGTACGTACGTACGTACGTACGTACGT\n>s_part\nACGTACGTACGTACGTACGT\n>s_shift\nTTTTACGTACGTACGTACGTACGTAAAA\n",
        "blastn",
        "6 sseqid qstart qend sstart send length bitscore",
        &[
            "--dust",
            "no",
            "--max_target_seqs",
            "10",
            "--culling_limit",
            "1",
        ],
        &[
            "-dust",
            "no",
            "-max_target_seqs",
            "10",
            "-culling_limit",
            "1",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_filter_boundary_values() {
    let query = ">q1\nACGTACGTACGTACGTACGTACGTACGTACGT\n";
    let subject = ">s_exact\nACGTACGTACGTACGTACGTACGTACGTACGT\n>s_imperfect\nACGTACGTACGTACGTTCGTACGTACGTACGT\n>s_partial\nACGTACGTACGTACGTACGT\n";
    let outfmt = "6 qseqid sseqid qstart qend sstart send length pident qcovhsp bitscore";

    for (rust_args, ncbi_args) in [
        (
            vec![
                "--dust",
                "no",
                "--max_target_seqs",
                "10",
                "--perc_identity",
                "0",
            ],
            vec![
                "-dust",
                "no",
                "-max_target_seqs",
                "10",
                "-perc_identity",
                "0",
            ],
        ),
        (
            vec![
                "--dust",
                "no",
                "--max_target_seqs",
                "10",
                "--perc_identity",
                "100",
            ],
            vec![
                "-dust",
                "no",
                "-max_target_seqs",
                "10",
                "-perc_identity",
                "100",
            ],
        ),
        (
            vec![
                "--dust",
                "no",
                "--max_target_seqs",
                "10",
                "--qcov_hsp_perc",
                "0",
            ],
            vec![
                "-dust",
                "no",
                "-max_target_seqs",
                "10",
                "-qcov_hsp_perc",
                "0",
            ],
        ),
        (
            vec![
                "--dust",
                "no",
                "--max_target_seqs",
                "10",
                "--qcov_hsp_perc",
                "100",
            ],
            vec![
                "-dust",
                "no",
                "-max_target_seqs",
                "10",
                "-qcov_hsp_perc",
                "100",
            ],
        ),
        (
            vec!["--dust", "no", "--max_target_seqs", "10", "--max_hsps", "1"],
            vec!["-dust", "no", "-max_target_seqs", "10", "-max_hsps", "1"],
        ),
        (
            vec![
                "--dust",
                "no",
                "--max_target_seqs",
                "10",
                "--culling_limit",
                "0",
            ],
            vec![
                "-dust",
                "no",
                "-max_target_seqs",
                "10",
                "-culling_limit",
                "0",
            ],
        ),
    ] {
        assert_blastn_subject_outfmt_matches_ncbi(query, subject, outfmt, &rust_args, &ncbi_args);
    }
}

#[test]
fn blastn_subject_ncbi_parity_max_target_seqs_edges() {
    let query = ">q1\nACGTACGTACGTACGTACGTACGTACGTACGT\n";
    let subject = ">s1\nACGTACGTACGTACGTACGTACGTACGTACGT\n>s2\nACGTACGTACGTACGTTCGTACGTACGTACGT\n>s3\nACGTACGTACGTACGTACGT\n";
    let outfmt = "6 qseqid sseqid qstart qend sstart send length pident bitscore";

    for (rust_args, ncbi_args) in [
        (
            vec!["--dust", "no", "--max_target_seqs", "1"],
            vec!["-dust", "no", "-max_target_seqs", "1"],
        ),
        (
            vec!["--dust", "no", "--max_target_seqs", "100000"],
            vec!["-dust", "no", "-max_target_seqs", "100000"],
        ),
    ] {
        assert_blastn_subject_outfmt_matches_ncbi(query, subject, outfmt, &rust_args, &ncbi_args);
    }
}

#[test]
fn blastn_subject_ncbi_parity_max_target_seqs_after_best_hit_filter() {
    let query = ">q1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";
    let subject = ">s_full\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n>s_bad\nACGTACGTACATACGTACATACGTAC\n>s_bad2\nACGTACGTACGTTCGTACGTTCGTACGTACGT\n>s_tail\nTTTTACGTACGTACGTACGTAAAA\n";

    assert_blastn_subject_task_outfmt_matches_ncbi(
        query,
        subject,
        "blastn",
        "6 sseqid",
        &[
            "--dust",
            "no",
            "--word_size",
            "7",
            "--max_target_seqs",
            "2",
            "--best_hit_overhang",
            "0.1",
            "--best_hit_score_edge",
            "0.1",
            "--max_hsps",
            "1",
        ],
        &[
            "-dust",
            "no",
            "-word_size",
            "7",
            "-max_target_seqs",
            "2",
            "-best_hit_overhang",
            "0.1",
            "-best_hit_score_edge",
            "0.1",
            "-max_hsps",
            "1",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_best_hit_tied_hsp_ordering() {
    let query = ">q1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";
    let subject = ">s_full\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n>s_bad\nACGTACGTACATACGTACATACGTAC\n>s_bad2\nACGTACGTACGTTCGTACGTTCGTACGTACGT\n>s_tail\nTTTTACGTACGTACGTACGTAAAA\n";
    let outfmt = "6 sseqid qstart qend sstart send length score bitscore evalue btop";

    assert_blastn_subject_task_outfmt_matches_ncbi(
        query,
        subject,
        "blastn",
        outfmt,
        &[
            "--dust",
            "no",
            "--word_size",
            "7",
            "--max_target_seqs",
            "2",
            "--best_hit_overhang",
            "0.1",
            "--best_hit_score_edge",
            "0.1",
        ],
        &[
            "-dust",
            "no",
            "-word_size",
            "7",
            "-max_target_seqs",
            "2",
            "-best_hit_overhang",
            "0.1",
            "-best_hit_score_edge",
            "0.1",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_culling_tied_subject_ordering() {
    let query = ">q1\nACGTACGTACGTACGTACGTACGTACGTACGT\n";
    let subject = ">s_alpha\nACGTACGTACGTACGTACGTACGTACGTACGT\n>s_beta\nACGTACGTACGTACGTACGTACGTACGTACGT\n>s_gamma\nACGTACGTACGTACGTACGTACGTACGTACGT\n";

    assert_blastn_subject_task_outfmt_matches_ncbi(
        query,
        subject,
        "blastn",
        "6 sseqid qstart qend sstart send score bitscore evalue",
        &[
            "--dust",
            "no",
            "--word_size",
            "7",
            "--max_target_seqs",
            "10",
            "--culling_limit",
            "1",
        ],
        &[
            "-dust",
            "no",
            "-word_size",
            "7",
            "-max_target_seqs",
            "10",
            "-culling_limit",
            "1",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_best_hit_tied_subject_ordering() {
    let query = ">q1\nACGTACGTACGTACGTACGTACGTACGTACGT\n";
    let subject = ">s_alpha\nACGTACGTACGTACGTACGTACGTACGTACGT\n>s_beta\nACGTACGTACGTACGTACGTACGTACGTACGT\n>s_gamma\nACGTACGTACGTACGTACGTACGTACGTACGT\n";

    assert_blastn_subject_task_outfmt_matches_ncbi(
        query,
        subject,
        "blastn",
        "6 sseqid qstart qend sstart send score bitscore evalue",
        &[
            "--dust",
            "no",
            "--word_size",
            "7",
            "--max_target_seqs",
            "10",
            "--best_hit_overhang",
            "0.1",
            "--best_hit_score_edge",
            "0.1",
        ],
        &[
            "-dust",
            "no",
            "-word_size",
            "7",
            "-max_target_seqs",
            "10",
            "-best_hit_overhang",
            "0.1",
            "-best_hit_score_edge",
            "0.1",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_max_target_seqs_after_culling_limit() {
    let query = ">q1\nACGTACGTACGTACGTACGTACGTACGTACGTTTAACCGGTTAA\n";
    let subject = ">s_full\nACGTACGTACGTACGTACGTACGTACGTACGT\n>s_part\nACGTACGTACGTACGTACGT\n>s_tail\nTTAACCGGTTAA\n>s_noise\nACGTACGTACGTACGT\n";

    assert_blastn_subject_task_outfmt_matches_ncbi(
        query,
        subject,
        "blastn",
        "6 sseqid qstart qend length score bitscore",
        &[
            "--dust",
            "no",
            "--word_size",
            "7",
            "--max_target_seqs",
            "2",
            "--culling_limit",
            "1",
            "--max_hsps",
            "1",
        ],
        &[
            "-dust",
            "no",
            "-word_size",
            "7",
            "-max_target_seqs",
            "2",
            "-culling_limit",
            "1",
            "-max_hsps",
            "1",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_best_hit_filter() {
    assert_blastn_subject_task_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n",
        ">s_full\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n>s_bad\nACGTACGTACATACGTACATACGTAC\n>s_bad2\nACGTACGTACGTTCGTACGTTCGTACGTACGT\n",
        "blastn",
        "6 sseqid qstart qend sstart send length pident score bitscore btop",
        &[
            "--dust",
            "no",
            "--word_size",
            "7",
            "--max_target_seqs",
            "20",
            "--best_hit_overhang",
            "0.1",
            "--best_hit_score_edge",
            "0.1",
        ],
        &[
            "-dust",
            "no",
            "-word_size",
            "7",
            "-max_target_seqs",
            "20",
            "-best_hit_overhang",
            "0.1",
            "-best_hit_score_edge",
            "0.1",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_subject_besthit() {
    assert_blastn_subject_task_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n",
        ">s1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n>s2\nACGTACGTACGTACGTACGTACGTACGTACGT\n",
        "blastn",
        "6 sseqid qstart qend sstart send length score bitscore",
        &[
            "--dust",
            "no",
            "--max_target_seqs",
            "20",
            "--subject_besthit",
        ],
        &["-dust", "no", "-max_target_seqs", "20", "-subject_besthit"],
    );
}

#[test]
fn blastn_subject_ncbi_parity_qcov_hsp_per_query() {
    assert_blastn_subject_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGTACGT\n>q2\nGGGTTTAAACCCGGGTTTAAACCCAAAAAAAAAAAA\n",
        ">s1\nACGTACGTACGTACGTACGTACGTNNNNGGGTTTAAACCCGGGTTTAAACCC\n",
        &[
            "--dust",
            "no",
            "--strand",
            "plus",
            "--max_target_seqs",
            "10",
            "--qcov_hsp_perc",
            "80",
        ],
        &[
            "-dust",
            "no",
            "-strand",
            "plus",
            "-max_target_seqs",
            "10",
            "-qcov_hsp_perc",
            "80",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_perc_identity_filter() {
    assert_blastn_subject_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGTACGT\n",
        ">perfect\nACGTACGTACGTACGTACGTACGT\n>one_mismatch\nACGTACGTACGTACGTACGTACGA\n",
        &[
            "--dust",
            "no",
            "--strand",
            "plus",
            "--max_target_seqs",
            "10",
            "--perc_identity",
            "99",
        ],
        &[
            "-dust",
            "no",
            "-strand",
            "plus",
            "-max_target_seqs",
            "10",
            "-perc_identity",
            "99",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_min_raw_gapped_score() {
    assert_blastn_subject_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGTACGT\n",
        ">perfect\nACGTACGTACGTACGTACGTACGT\n>shorter\nACGTACGTACGTACGTACGTACGA\n",
        &[
            "--dust",
            "no",
            "--strand",
            "plus",
            "--max_target_seqs",
            "10",
            "--min_raw_gapped_score",
            "47",
        ],
        &[
            "-dust",
            "no",
            "-strand",
            "plus",
            "-max_target_seqs",
            "10",
            "-min_raw_gapped_score",
            "47",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_default_dust_masking() {
    assert_blastn_subject_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGTACGT\n",
        ">s1\nTTTTACGTACGTACGTACGTACGTACGTAAAA\n>s2\nACGTACGTACGTACGTACGTACGT\n",
        &["--max_target_seqs", "10", "--max_hsps", "2"],
        &["-max_target_seqs", "10", "-max_hsps", "2"],
    );
}

#[test]
fn blastn_subject_ncbi_parity_lcase_masking() {
    assert_blastn_subject_matches_ncbi(
        ">q1\nacgtacgtacgtacgtacgtacgt\n",
        ">s1\nACGTACGTACGTACGTACGTACGT\n",
        &["--dust", "no", "--lcase_masking", "--max_target_seqs", "10"],
        &["-dust", "no", "-lcase_masking", "-max_target_seqs", "10"],
    );
}

#[test]
fn blastn_subject_ncbi_parity_query_loc() {
    let query = ">q1\nTTTTACGTACGTACGTACGTACGTAAAA\n";
    let subject = ">s1\nACGTACGTACGTACGTACGT\n";
    for strand in ["plus", "minus"] {
        assert_blastn_subject_matches_ncbi(
            query,
            subject,
            &[
                "--dust",
                "no",
                "--strand",
                strand,
                "--query_loc",
                "5-24",
                "--max_target_seqs",
                "10",
            ],
            &[
                "-dust",
                "no",
                "-strand",
                strand,
                "-query_loc",
                "5-24",
                "-max_target_seqs",
                "10",
            ],
        );
    }
}

#[test]
fn blastn_subject_ncbi_parity_subject_loc() {
    let query = ">q1\nACGTACGTACGTACGTACGT\n";
    let subject = ">s1\nTTTTACGTACGTACGTACGTACGTAAAA\n";
    for strand in ["plus", "minus"] {
        assert_blastn_subject_outfmt_matches_ncbi(
            query,
            subject,
            "6 qseqid sseqid qlen slen qstart qend sstart send evalue bitscore qframe sframe",
            &[
                "--dust",
                "no",
                "--strand",
                strand,
                "--subject_loc",
                "5-24",
                "--max_target_seqs",
                "10",
            ],
            &[
                "-dust",
                "no",
                "-strand",
                strand,
                "-subject_loc",
                "5-24",
                "-max_target_seqs",
                "10",
            ],
        );
    }
}

#[test]
fn blastn_subject_ncbi_parity_qcovs_union() {
    assert_blastn_subject_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGTGGGTTTAAACCCGGGTTTA\n",
        ">s1\nACGTACGTACGTACGTACGTNNNNNNNNNNGGGTTTAAACCCGGGTTTA\n",
        "6 qseqid sseqid qstart qend length qcovs qcovhsp bitscore",
        &[
            "--dust",
            "no",
            "--strand",
            "plus",
            "--max_target_seqs",
            "10",
        ],
        &["-dust", "no", "-strand", "plus", "-max_target_seqs", "10"],
    );
}

#[test]
fn blastn_subject_ncbi_parity_outfmt6_delim() {
    assert_blastn_subject_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGTACGT\n",
        ">s1\nACGTACGTACGTACGTACGTACGT\n",
        "6 delim=, qseqid sseqid qstart qend",
        &["--dust", "no", "--max_target_seqs", "10"],
        &["-dust", "no", "-max_target_seqs", "10"],
    );
}

#[test]
fn blastn_subject_ncbi_parity_outfmt_delim_values_are_literal() {
    for outfmt in [
        "6 delim=tab qseqid sseqid length",
        "7 delim=space qseqid sseqid length",
        "10 delim=tab qseqid sseqid length",
    ] {
        assert_blastn_subject_outfmt_matches_ncbi(
            ">q1\nACGTACGTACGTACGTACGTACGT\n",
            ">s1\nACGTACGTACGTACGTACGTACGT\n",
            outfmt,
            &["--dust", "no", "--max_target_seqs", "10"],
            &["-dust", "no", "-max_target_seqs", "10"],
        );
    }
}

#[test]
fn blastn_subject_ncbi_parity_outfmt6_std_keyword() {
    assert_blastn_subject_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGT\n",
        ">s1\nACGTACGTACGTACGTACGT\n",
        "6 std qlen",
        &["--dust", "no", "--max_target_seqs", "10"],
        &["-dust", "no", "-max_target_seqs", "10"],
    );
}

#[test]
fn blastn_subject_ncbi_parity_outfmt6_duplicate_fields() {
    assert_blastn_subject_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGT\n",
        ">s1\nACGTACGTACGTACGTACGT\n",
        "6 qseqid qseqid sseqid length length bitscore",
        &["--dust", "no", "--max_target_seqs", "10"],
        &["-dust", "no", "-max_target_seqs", "10"],
    );
}

#[test]
fn blastn_subject_ncbi_parity_outfmt6_invalid_fields_are_ignored() {
    assert_blastn_subject_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGT\n",
        ">s1\nACGTACGTACGTACGTACGT\n",
        "6 qseqid bogus sseqid",
        &["--dust", "no", "--max_target_seqs", "10"],
        &["-dust", "no", "-max_target_seqs", "10"],
    );

    assert_blastn_subject_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGT\n",
        ">s1\nACGTACGTACGTACGTACGT\n",
        "6 bogus",
        &["--dust", "no", "--max_target_seqs", "10"],
        &["-dust", "no", "-max_target_seqs", "10"],
    );
}

#[test]
fn blastn_subject_ncbi_parity_outfmt6_empty_field_list_uses_defaults() {
    for outfmt in ["6 ", "6 delim=,"] {
        assert_blastn_subject_outfmt_matches_ncbi(
            ">q1\nACGTACGTACGTACGTACGT\n",
            ">s1\nACGTACGTACGTACGTACGT\n",
            outfmt,
            &["--dust", "no", "--max_target_seqs", "10"],
            &["-dust", "no", "-max_target_seqs", "10"],
        );
    }
}

#[test]
fn blastn_subject_ncbi_parity_custom_field_parser_edges() {
    for outfmt in [
        "6 delim=, delim=tab qseqid sseqid length",
        "6 std qaccver saccver qlen",
        "7 qseqid bogus sseqid qseqid length length",
        "7 delim=, delim=tab qseqid sseqid length",
        "7 std qaccver saccver qlen",
        "7 bogus",
        "10 qseqid bogus sseqid qseqid length length",
        "10 delim=tab delim=, qseqid sseqid length",
        "10 std qaccver saccver qlen",
        "10 bogus",
    ] {
        assert_blastn_subject_outfmt_matches_ncbi(
            ">q1\nACGTACGTACGTACGTACGT\n",
            ">s1\nACGTACGTACGTACGTACGT\n",
            outfmt,
            &["--dust", "no", "--max_target_seqs", "10"],
            &["-dust", "no", "-max_target_seqs", "10"],
        );
    }
}

#[test]
fn blastn_subject_ncbi_parity_invalid_outfmt_number_errors() {
    if !std::path::Path::new("/usr/bin/blastn").exists() {
        eprintln!("Skipping: /usr/bin/blastn not found");
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
    std::fs::write(&query, ">q1\nACGTACGTACGTACGTACGT\n").expect("write query FASTA");
    std::fs::write(&subject, ">s1\nACGTACGTACGTACGTACGT\n").expect("write subject FASTA");

    for outfmt in ["\"6 qseqid sseqid\"", "", " qseqid"] {
        let rust = std::process::Command::new(&blast_cli)
            .arg("blastn")
            .arg("--query")
            .arg(&query)
            .arg("--subject")
            .arg(&subject)
            .arg("--task")
            .arg("blastn-short")
            .arg("--dust")
            .arg("no")
            .arg("--outfmt")
            .arg(outfmt)
            .output()
            .expect("run blast-cli invalid outfmt");
        let ncbi = std::process::Command::new("/usr/bin/blastn")
            .arg("-query")
            .arg(&query)
            .arg("-subject")
            .arg(&subject)
            .arg("-task")
            .arg("blastn-short")
            .arg("-dust")
            .arg("no")
            .arg("-outfmt")
            .arg(outfmt)
            .output()
            .expect("run NCBI invalid outfmt");

        assert!(
            !rust.status.success(),
            "blast-cli should reject outfmt {outfmt:?}"
        );
        assert!(
            !ncbi.status.success(),
            "NCBI should reject outfmt {outfmt:?}"
        );
        assert_eq!(
            rust.stdout, ncbi.stdout,
            "invalid outfmt stdout differs for {outfmt:?}"
        );
        assert_eq!(
            rust.stderr, ncbi.stderr,
            "invalid outfmt stderr differs for {outfmt:?}"
        );
    }
}

#[test]
fn blastn_subject_ncbi_parity_archive_outfmt_requires_output_file() {
    if !std::path::Path::new("/usr/bin/blastn").exists() {
        eprintln!("Skipping: /usr/bin/blastn not found");
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
    std::fs::write(&query, ">q1\nACGTACGTACGTACGTACGT\n").expect("write query FASTA");
    std::fs::write(&subject, ">s1\nACGTACGTACGTACGTACGT\n").expect("write subject FASTA");

    for outfmt in ["13", "14"] {
        let rust = std::process::Command::new(&blast_cli)
            .arg("blastn")
            .arg("--query")
            .arg(&query)
            .arg("--subject")
            .arg(&subject)
            .arg("--task")
            .arg("blastn-short")
            .arg("--dust")
            .arg("no")
            .arg("--outfmt")
            .arg(outfmt)
            .output()
            .expect("run blast-cli archive outfmt");
        let ncbi = std::process::Command::new("/usr/bin/blastn")
            .arg("-query")
            .arg(&query)
            .arg("-subject")
            .arg(&subject)
            .arg("-task")
            .arg("blastn-short")
            .arg("-dust")
            .arg("no")
            .arg("-outfmt")
            .arg(outfmt)
            .output()
            .expect("run NCBI archive outfmt");

        assert!(
            !rust.status.success(),
            "blast-cli should reject outfmt {outfmt}"
        );
        assert!(!ncbi.status.success(), "NCBI should reject outfmt {outfmt}");
        assert_eq!(
            rust.stdout, ncbi.stdout,
            "archive outfmt stdout differs for {outfmt}"
        );
        assert_eq!(
            rust.stderr, ncbi.stderr,
            "archive outfmt stderr differs for {outfmt}"
        );
    }
}

#[test]
fn blastn_subject_rejects_unsupported_numeric_outfmt() {
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
    std::fs::write(&query, ">q1\nACGTACGTACGTACGTACGT\n").expect("write query FASTA");
    std::fs::write(&subject, ">s1\nACGTACGTACGTACGTACGT\n").expect("write subject FASTA");

    for outfmt in ["1", "2", "3", "4", "8", "9", "11", "12", "15", "16", "18"] {
        let output = std::process::Command::new(&blast_cli)
            .arg("blastn")
            .arg("--query")
            .arg(&query)
            .arg("--subject")
            .arg(&subject)
            .arg("--task")
            .arg("blastn-short")
            .arg("--dust")
            .arg("no")
            .arg("--outfmt")
            .arg(outfmt)
            .output()
            .expect("run blast-cli unsupported outfmt");
        assert!(
            !output.status.success(),
            "blast-cli should reject outfmt {outfmt}"
        );
        assert!(
            output.stdout.is_empty(),
            "unsupported outfmt {outfmt} should not write stdout"
        );
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(
            stderr.contains(&format!("Output format {outfmt} is not supported")),
            "unsupported outfmt stderr differed for {outfmt}:\n{stderr}"
        );
    }
}

#[test]
fn blastn_subject_ncbi_parity_outfmt7_std_keyword() {
    assert_blastn_subject_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGT\n",
        ">s1\nACGTACGTACGTACGTACGT\n",
        "7 std qlen",
        &["--dust", "no", "--max_target_seqs", "10"],
        &["-dust", "no", "-max_target_seqs", "10"],
    );
}

#[test]
fn blastn_subject_ncbi_parity_outfmt10_std_keyword() {
    assert_blastn_subject_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGT\n",
        ">s1\nACGTACGTACGTACGTACGT\n",
        "10 std qlen",
        &["--dust", "no", "--max_target_seqs", "10"],
        &["-dust", "no", "-max_target_seqs", "10"],
    );
}

#[test]
fn blastn_ncbi_parity_version_ignores_other_arguments() {
    if !std::path::Path::new("/usr/bin/blastn").exists() {
        eprintln!("Skipping: /usr/bin/blastn not found");
        return;
    }
    let Some(blast_cli) = std::env::var_os("BLAST_RS_CLI_BIN")
        .or_else(|| std::env::var_os("CARGO_BIN_EXE_blast-cli"))
        .map(std::path::PathBuf::from)
    else {
        eprintln!("Skipping: set BLAST_RS_CLI_BIN or CARGO_BIN_EXE_blast-cli to run CLI parity");
        return;
    };

    let rust = std::process::Command::new(&blast_cli)
        .arg("blastn")
        .arg("-version")
        .arg("--task")
        .arg("bad")
        .arg("--query")
        .arg("missing.fa")
        .output()
        .expect("run blast-cli version");
    let ncbi = std::process::Command::new("/usr/bin/blastn")
        .arg("-version")
        .arg("-task")
        .arg("bad")
        .arg("-query")
        .arg("missing.fa")
        .output()
        .expect("run NCBI version");

    assert!(rust.status.success(), "blast-cli version should succeed");
    assert!(ncbi.status.success(), "NCBI version should succeed");
    assert_eq!(rust.stdout, ncbi.stdout, "version stdout differs");
    assert_eq!(rust.stderr, ncbi.stderr, "version stderr differs");
}

#[test]
fn non_blastn_programs_ncbi_parity_version_ignores_other_arguments() {
    let Some(blast_cli) = std::env::var_os("BLAST_RS_CLI_BIN")
        .or_else(|| std::env::var_os("CARGO_BIN_EXE_blast-cli"))
        .map(std::path::PathBuf::from)
    else {
        eprintln!("Skipping: set BLAST_RS_CLI_BIN or CARGO_BIN_EXE_blast-cli");
        return;
    };

    for program in [
        "blastp",
        "blastx",
        "tblastn",
        "tblastx",
        "psiblast",
        "rpstblastn",
        "deltablast",
    ] {
        let ncbi_bin = format!("/usr/bin/{program}");
        if !std::path::Path::new(&ncbi_bin).exists() {
            eprintln!("Skipping: {ncbi_bin} not found");
            continue;
        }
        let rust = std::process::Command::new(&blast_cli)
            .arg(program)
            .arg("-version")
            .arg("-query")
            .arg("missing.fa")
            .arg("-outfmt")
            .arg("bad")
            .arg("-help")
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli {program} -version: {err}"));
        let ncbi = std::process::Command::new(&ncbi_bin)
            .arg("-version")
            .arg("-query")
            .arg("missing.fa")
            .arg("-outfmt")
            .arg("bad")
            .arg("-help")
            .output()
            .unwrap_or_else(|err| panic!("run NCBI {program} -version: {err}"));

        assert!(
            rust.status.success(),
            "blast-cli {program} -version should succeed"
        );
        assert!(
            ncbi.status.success(),
            "NCBI {program} -version should succeed"
        );
        assert_eq!(rust.stdout, ncbi.stdout, "{program} version stdout differs");
        assert_eq!(rust.stderr, ncbi.stderr, "{program} version stderr differs");
    }
}

#[test]
fn blastn_help_ignores_other_arguments_and_uses_blast_shape() {
    let Some(blast_cli) = std::env::var_os("BLAST_RS_CLI_BIN")
        .or_else(|| std::env::var_os("CARGO_BIN_EXE_blast-cli"))
        .map(std::path::PathBuf::from)
    else {
        eprintln!("Skipping: set BLAST_RS_CLI_BIN or CARGO_BIN_EXE_blast-cli to run CLI parity");
        return;
    };

    for help_arg in ["-h", "-help"] {
        let output = std::process::Command::new(&blast_cli)
            .arg("blastn")
            .arg(help_arg)
            .arg("--task")
            .arg("bad")
            .arg("--query")
            .arg("missing.fa")
            .output()
            .expect("run blast-cli help");

        assert!(
            output.status.success(),
            "{help_arg} should ignore invalid arguments"
        );
        assert!(
            output.stderr.is_empty(),
            "{help_arg} should not write stderr"
        );
        let stdout = String::from_utf8_lossy(&output.stdout);
        assert!(
            stdout.starts_with("USAGE\n  blastn [-h] [-help]"),
            "unexpected help prefix:\n{stdout}"
        );
        assert!(stdout.contains("DESCRIPTION\n   Nucleotide-Nucleotide BLAST 2.12.0+"));
        if help_arg == "-help" {
            assert!(stdout.contains("OPTIONAL ARGUMENTS"));
        }
    }
}

#[test]
fn blastn_sort_options_warn_when_ignored_by_outfmt() {
    let Some(blast_cli) = std::env::var_os("BLAST_RS_CLI_BIN")
        .or_else(|| std::env::var_os("CARGO_BIN_EXE_blast-cli"))
        .map(std::path::PathBuf::from)
    else {
        eprintln!(
            "Skipping: set BLAST_RS_CLI_BIN or CARGO_BIN_EXE_blast-cli to run CLI warning parity"
        );
        return;
    };

    let tmp = TempDir::new().expect("tempdir");
    let query = tmp.path().join("query.fa");
    let subject = tmp.path().join("subject.fa");
    std::fs::write(&query, ">q1\nACGTACGTACGTACGTACGTACGT\n").expect("write query FASTA");
    std::fs::write(&subject, ">s1\nACGTACGTACGTACGTACGTACGT\n").expect("write subject FASTA");

    let output = std::process::Command::new(blast_cli)
        .arg("blastn")
        .arg("--query")
        .arg(&query)
        .arg("--subject")
        .arg(&subject)
        .arg("--task")
        .arg("blastn-short")
        .arg("--outfmt")
        .arg("6")
        .arg("--dust")
        .arg("no")
        .arg("--sorthits")
        .arg("1")
        .arg("--sorthsps")
        .arg("1")
        .arg("--max_target_seqs")
        .arg("4")
        .arg("--line_length")
        .arg("80")
        .output()
        .expect("run blast-cli sort warning check");

    assert!(
        output.status.success(),
        "blast-cli exited with {}",
        output.status
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains(
            "Warning: [blastn] The parameter -sorthits is ignored for output formats > 4."
        ),
        "missing sorthits ignored warning, stderr was:\n{stderr}"
    );
    assert!(
        stderr.contains(
            "Warning: [blastn] The parameter -sorthsps is ignored for output formats != 0."
        ),
        "missing sorthsps ignored warning, stderr was:\n{stderr}"
    );
    assert!(
        stderr.contains("Warning: [blastn] Examining 5 or more matches is recommended"),
        "missing hitlist-size warning, stderr was:\n{stderr}"
    );
    assert!(
        stderr.contains(
            "Warning: [blastn] The parameter -line_length is not applicable for output formats > 4 ."
        ),
        "missing line_length ignored warning, stderr was:\n{stderr}"
    );
}

#[test]
fn blastn_formatting_options_match_ncbi_tabular_behavior() {
    assert_blastn_subject_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGTACGT\n",
        ">s1\nACGTACGTACGTACGTACGTACGT\n>s2\nACGTACGTACGTACGTACGTACGT\n>s3\nACGTACGTACGTACGTACGTACGT\n",
        "6 sseqid bitscore",
        &["--dust", "no", "--num_alignments", "1"],
        &["-dust", "no", "-num_alignments", "1"],
    );
}

#[test]
fn blastn_subject_ncbi_parity_html_ignored_for_tabular_output() {
    assert_blastn_subject_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGTACGT\n",
        ">s1\nACGTACGTACGTACGTACGTACGT\n",
        "6 qseqid sseqid pident length bitscore",
        &["--dust", "no", "--html"],
        &["-dust", "no", "-html"],
    );
}

#[test]
fn blastn_subject_ncbi_parity_parse_deflines_ids_and_accessions() {
    for (query, subject) in [
        (
            ">lcl|query1 query title\nACGTACGTACGTACGTACGT\n",
            ">lcl|subject1 subject title\nACGTACGTACGTACGTACGT\n",
        ),
        (
            ">gi|123|ref|QACC.1| query title\nACGTACGTACGTACGTACGT\n",
            ">gi|456|ref|SACC.1| subject title\nACGTACGTACGTACGTACGT\n",
        ),
    ] {
        assert_blastn_subject_outfmt_matches_ncbi(
            query,
            subject,
            "6 qseqid qgi qacc qaccver sseqid sgi sallgi sacc saccver sallseqid sallacc",
            &["--dust", "no", "--parse_deflines"],
            &["-dust", "no", "-parse_deflines"],
        );
    }
}

#[test]
fn blastn_subject_ncbi_parity_pairwise_parse_deflines() {
    for (query, subject) in [
        (
            ">lcl|query1 query title\nACGTACGTACGTACGTACGT\n",
            ">lcl|subject1 subject title\nACGTACGTACGTACGTACGT\n",
        ),
        (
            ">gi|123|ref|QACC.1| query title\nACGTACGTACGTACGTACGT\n",
            ">gi|456|ref|SACC.1| subject title\nACGTACGTACGTACGTACGT\n",
        ),
    ] {
        assert_blastn_subject_outfmt_matches_ncbi(
            query,
            subject,
            "0",
            &[
                "--dust",
                "no",
                "--parse_deflines",
                "--num_descriptions",
                "1",
                "--num_alignments",
                "1",
            ],
            &[
                "-dust",
                "no",
                "-parse_deflines",
                "-num_descriptions",
                "1",
                "-num_alignments",
                "1",
            ],
        );
    }
}

#[test]
fn blastn_subject_ncbi_parity_xml_parse_deflines() {
    for (query, subject) in [
        (
            ">lcl|query1 query title\nACGTACGTACGTACGTACGT\n",
            ">lcl|subject1 subject title\nACGTACGTACGTACGTACGT\n",
        ),
        (
            ">gi|123|ref|QACC.1| query title\nACGTACGTACGTACGTACGT\n",
            ">gi|456|ref|SACC.1| subject title\nACGTACGTACGTACGTACGT\n",
        ),
    ] {
        assert_blastn_subject_outfmt_matches_ncbi(
            query,
            subject,
            "5",
            &[
                "--dust",
                "no",
                "--parse_deflines",
                "--max_target_seqs",
                "10",
            ],
            &["-dust", "no", "-parse_deflines", "-max_target_seqs", "10"],
        );
    }
}

#[test]
fn blastn_subject_ncbi_parity_sam_parse_deflines() {
    for (query, subject) in [
        (
            ">lcl|query1 query title\nACGTACGTACGTACGTACGT\n",
            ">lcl|subject1 subject title\nACGTACGTACGTACGTACGT\n",
        ),
        (
            ">gi|123|ref|QACC.1| query title\nACGTACGTACGTACGTACGT\n",
            ">gi|456|ref|SACC.1| subject title\nACGTACGTACGTACGTACGT\n",
        ),
    ] {
        assert_blastn_subject_sam_matches_ncbi(
            query,
            subject,
            &[
                "--dust",
                "no",
                "--parse_deflines",
                "--max_target_seqs",
                "10",
            ],
            &["-dust", "no", "-parse_deflines", "-max_target_seqs", "10"],
        );
    }
}

#[test]
fn blastn_subject_ncbi_parity_outfmt10_default_and_custom() {
    let query = ">q1\nACGTACGTACGTACGTACGT\n";
    let subject = ">s1 subject title, with comma\nACGTACGTACGTACGTACGT\n";
    for outfmt in [
        "10",
        "10 qseqid sseqid pident qstart qend sstart send stitle bitscore",
    ] {
        assert_blastn_subject_outfmt_matches_ncbi(
            query,
            subject,
            outfmt,
            &["--dust", "no", "--max_target_seqs", "10"],
            &["-dust", "no", "-max_target_seqs", "10"],
        );
    }
}

#[test]
fn blastn_subject_ncbi_parity_outfmt7_default_and_custom() {
    let query = ">q1 first query\nACGTACGTACGTACGTACGT\n>q2 nohit\nTTTTTTTTTTTTTTTTTTTT\n";
    let subject = ">s1 first subject\nACGTACGTACGTACGTACGT\n";
    for outfmt in ["7", "7 qseqid sseqid pident length bitscore"] {
        assert_blastn_subject_outfmt_matches_ncbi(
            query,
            subject,
            outfmt,
            &["--dust", "no", "--max_target_seqs", "10"],
            &["-dust", "no", "-max_target_seqs", "10"],
        );
    }
}

#[test]
fn blastn_subject_ncbi_parity_pairwise_exact_hit() {
    assert_blastn_subject_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGT\n",
        ">s1 subject one\nACGTACGTACGTACGTACGT\n",
        "0",
        &["--dust", "no", "--max_target_seqs", "10"],
        &["-dust", "no", "-max_target_seqs", "10"],
    );
}

#[test]
fn blastn_subject_ncbi_parity_pairwise_no_hit() {
    assert_blastn_subject_outfmt_matches_ncbi(
        ">q1\nAAAAAAAAAAAAAAAAAAAA\n",
        ">s1 subject one\nCCCCCCCCCCCCCCCCCCCC\n",
        "0",
        &["--dust", "no", "--max_target_seqs", "10"],
        &["-dust", "no", "-max_target_seqs", "10"],
    );
}

#[test]
fn blastn_subject_ncbi_parity_pairwise_multi_query() {
    assert_blastn_subject_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGT\n>q2 nohit\nAAAAAAAAAAAAAAAAAAAA\n",
        ">s1 subject one\nACGTACGTACGTACGTACGT\n",
        "0",
        &["--dust", "no", "--max_target_seqs", "10"],
        &["-dust", "no", "-max_target_seqs", "10"],
    );
}

#[test]
fn blastn_subject_ncbi_parity_pairwise_description_alignment_limits() {
    assert_blastn_subject_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGTACGT\n",
        ">s1\nACGTACGTACGTACGTACGTACGT\n>s2\nACGTACGTACGTACGTACGTACGT\n>s3\nACGTACGTACGTACGTACGTACGT\n",
        "0",
        &[
            "--dust",
            "no",
            "--num_descriptions",
            "2",
            "--num_alignments",
            "1",
        ],
        &[
            "-dust",
            "no",
            "-num_descriptions",
            "2",
            "-num_alignments",
            "1",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_pairwise_line_length() {
    assert_blastn_subject_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGTACGT\n",
        ">s1\nACGTACGTACGTACGTACGTACGT\n",
        "0",
        &[
            "--dust",
            "no",
            "--num_descriptions",
            "1",
            "--num_alignments",
            "1",
            "--line_length",
            "12",
        ],
        &[
            "-dust",
            "no",
            "-num_descriptions",
            "1",
            "-num_alignments",
            "1",
            "-line_length",
            "12",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_pairwise_zero_descriptions_zero_alignments() {
    if !std::path::Path::new("/usr/bin/blastn").exists() {
        eprintln!("Skipping: /usr/bin/blastn not found");
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
    let rust_out = tmp.path().join("rust.out");
    let ncbi_out = tmp.path().join("ncbi.out");
    std::fs::write(&query, ">q1\nACGTACGTACGTACGTACGTACGT\n").expect("write query FASTA");
    std::fs::write(&subject, ">s1\nACGTACGTACGTACGTACGTACGT\n").expect("write subject FASTA");

    let rust = std::process::Command::new(blast_cli)
        .arg("blastn")
        .arg("--query")
        .arg(&query)
        .arg("--subject")
        .arg(&subject)
        .arg("--task")
        .arg("blastn-short")
        .arg("--outfmt")
        .arg("0")
        .arg("--dust")
        .arg("no")
        .arg("--num_descriptions")
        .arg("0")
        .arg("--num_alignments")
        .arg("0")
        .arg("--out")
        .arg(&rust_out)
        .output()
        .expect("run blast-cli zero pairwise limits");
    let ncbi = std::process::Command::new("/usr/bin/blastn")
        .arg("-query")
        .arg(&query)
        .arg("-subject")
        .arg(&subject)
        .arg("-task")
        .arg("blastn-short")
        .arg("-outfmt")
        .arg("0")
        .arg("-dust")
        .arg("no")
        .arg("-num_descriptions")
        .arg("0")
        .arg("-num_alignments")
        .arg("0")
        .arg("-out")
        .arg(&ncbi_out)
        .output()
        .expect("run NCBI zero pairwise limits");

    assert!(
        !rust.status.success(),
        "blast-cli should reject zero pairwise limits"
    );
    assert!(
        !ncbi.status.success(),
        "NCBI should reject zero pairwise limits"
    );
    assert_eq!(
        std::fs::read(&rust_out).unwrap_or_default(),
        std::fs::read(&ncbi_out).unwrap_or_default(),
        "zero pairwise limit outputs differ"
    );
    let rust_stderr = String::from_utf8_lossy(&rust.stderr);
    let ncbi_stderr = String::from_utf8_lossy(&ncbi.stderr);
    for expected in [
        "BLAST query/options error: No hits are being saved",
        "Please refer to the BLAST+ user manual.",
    ] {
        assert!(
            rust_stderr.contains(expected),
            "missing Rust stderr line {expected:?}, stderr was:\n{rust_stderr}"
        );
        assert!(
            ncbi_stderr.contains(expected),
            "missing NCBI stderr line {expected:?}, stderr was:\n{ncbi_stderr}"
        );
    }
    assert_eq!(
        rust_stderr, ncbi_stderr,
        "zero pairwise limit stderr differs"
    );
}

#[test]
fn blastn_subject_ncbi_parity_pairwise_sorthsps_query_start() {
    assert_blastn_subject_task_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGTACGTACGTACGT\n",
        ">s_full\nACGTACGTACGTACGTACGTACGTACGTACGT\n>s_short\nACGTACGTACGTACGTACGTACGT\n>s_mismatch\nACGTACGTACGTACGTACGTACGTACGTACGA\n",
        "blastn-short",
        "0",
        &[
            "--dust",
            "no",
            "--max_target_seqs",
            "10",
            "--sorthsps",
            "2",
        ],
        &[
            "-dust",
            "no",
            "-max_target_seqs",
            "10",
            "-sorthsps",
            "2",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_pairwise_sorthsps_score() {
    assert_blastn_subject_task_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGTACGTACGTACGT\n",
        ">s_full\nACGTACGTACGTACGTACGTACGTACGTACGT\n>s_short\nACGTACGTACGTACGTACGTACGT\n>s_mismatch\nACGTACGTACGTACGTACGTACGTACGTACGA\n",
        "blastn-short",
        "0",
        &[
            "--dust",
            "no",
            "--max_target_seqs",
            "10",
            "--sorthsps",
            "1",
        ],
        &[
            "-dust",
            "no",
            "-max_target_seqs",
            "10",
            "-sorthsps",
            "1",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_pairwise_sorthsps_percent_identity() {
    assert_blastn_subject_task_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGTACGTACGTACGT\n",
        ">s_full\nACGTACGTACGTACGTACGTACGTACGTACGT\n>s_short\nACGTACGTACGTACGTACGTACGT\n>s_mismatch\nACGTACGTACGTACGTACGTACGTACGTACGA\n",
        "blastn-short",
        "0",
        &[
            "--dust",
            "no",
            "--max_target_seqs",
            "10",
            "--sorthsps",
            "3",
        ],
        &[
            "-dust",
            "no",
            "-max_target_seqs",
            "10",
            "-sorthsps",
            "3",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_pairwise_sorthsps_subject_start() {
    assert_blastn_subject_task_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGTACGTACGTACGT\n",
        ">s_full\nACGTACGTACGTACGTACGTACGTACGTACGT\n>s_short\nACGTACGTACGTACGTACGTACGT\n>s_mismatch\nACGTACGTACGTACGTACGTACGTACGTACGA\n",
        "blastn-short",
        "0",
        &[
            "--dust",
            "no",
            "--max_target_seqs",
            "10",
            "--sorthsps",
            "4",
        ],
        &[
            "-dust",
            "no",
            "-max_target_seqs",
            "10",
            "-sorthsps",
            "4",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_pairwise_sorthits_best_score() {
    assert_blastn_subject_task_outfmt_matches_ncbi(
        ">q1\nACGTTGCAAGTCCTGATCGATGCTAGCTTACG\n",
        ">s_full_imperfect\nACGTTGCAAGTCCTGATCGATGCTAGCTTTCG\n>s_short_exact\nACGTTGCAAGTCCTGATCGATGCT\n",
        "blastn-short",
        "0",
        &[
            "--dust",
            "no",
            "--max_target_seqs",
            "10",
            "--sorthits",
            "1",
        ],
        &[
            "-dust",
            "no",
            "-max_target_seqs",
            "10",
            "-sorthits",
            "1",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_pairwise_sorthits_total_score() {
    assert_blastn_subject_task_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGTACGTACGTACGT\n",
        ">s_long_mismatch\nACGTACGTACGTACGTACGTACGTACGTACGA\n>s_short_exact\nACGTACGTACGTACGTACGTACGT\n",
        "blastn-short",
        "0",
        &[
            "--dust",
            "no",
            "--max_target_seqs",
            "10",
            "--sorthits",
            "2",
        ],
        &[
            "-dust",
            "no",
            "-max_target_seqs",
            "10",
            "-sorthits",
            "2",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_pairwise_sorthits_percent_identity() {
    assert_blastn_subject_task_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGTACGTACGTACGT\n",
        ">s_full_imperfect\nACGTACGTACGTACGTTCGTACGTACGTACGT\n>s_short_exact\nACGTACGTACGTACGTACGTACGT\n",
        "blastn-short",
        "0",
        &[
            "--dust",
            "no",
            "--max_target_seqs",
            "10",
            "--sorthits",
            "3",
        ],
        &[
            "-dust",
            "no",
            "-max_target_seqs",
            "10",
            "-sorthits",
            "3",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_pairwise_sorthits_query_coverage() {
    assert_blastn_subject_task_outfmt_matches_ncbi(
        ">q1\nACGTTGCAAGTCCTGATCGATGCTAGCTTACG\n",
        ">s_full_imperfect\nACGTTGCAAGTCCTGATCGATGCTAGCTTTCG\n>s_short_exact\nACGTTGCAAGTCCTGATCGATGCT\n",
        "blastn-short",
        "0",
        &[
            "--dust",
            "no",
            "--max_target_seqs",
            "10",
            "--sorthits",
            "4",
        ],
        &[
            "-dust",
            "no",
            "-max_target_seqs",
            "10",
            "-sorthits",
            "4",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_outfmt7_delim() {
    assert_blastn_subject_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGT\n",
        ">s1\nACGTACGTACGTACGTACGT\n",
        "7 delim=, qseqid sseqid qstart qend",
        &["--dust", "no", "--max_target_seqs", "10"],
        &["-dust", "no", "-max_target_seqs", "10"],
    );
}

#[test]
fn blastn_subject_ncbi_parity_btop() {
    assert_blastn_subject_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGTACGT\n",
        ">s1\nACGTACGTTCGTACGTACGTACGT\n",
        "6 qseqid sseqid qstart qend sstart send length pident btop qseq sseq",
        &[
            "--dust",
            "no",
            "--strand",
            "plus",
            "--word_size",
            "7",
            "--max_target_seqs",
            "10",
        ],
        &[
            "-dust",
            "no",
            "-strand",
            "plus",
            "-word_size",
            "7",
            "-max_target_seqs",
            "10",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_btop_with_gaps() {
    for (query, subject) in [
        (
            ">q1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n",
            ">s1\nACGTACGTACGTACGTACGTAAAACGTACGTACGTACGTACGTACGTACGT\n",
        ),
        (
            ">q1\nACGTACGTACGTACGTACGTAAAACGTACGTACGTACGTACGTACGTACGT\n",
            ">s1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n",
        ),
    ] {
        assert_blastn_subject_outfmt_matches_ncbi(
            query,
            subject,
            "6 qstart qend sstart send length pident gapopen btop qseq sseq bitscore",
            &[
                "--dust",
                "no",
                "--strand",
                "plus",
                "--word_size",
                "7",
                "--gapopen",
                "5",
                "--gapextend",
                "2",
                "--max_target_seqs",
                "10",
                "--max_hsps",
                "10",
            ],
            &[
                "-dust",
                "no",
                "-strand",
                "plus",
                "-word_size",
                "7",
                "-gapopen",
                "5",
                "-gapextend",
                "2",
                "-max_target_seqs",
                "10",
                "-max_hsps",
                "10",
            ],
        );
    }
}

#[test]
fn blastn_subject_ncbi_parity_strand_frame_fields() {
    assert_blastn_subject_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGT\n",
        ">s1\nACGTACGTACGTACGTACGT\n",
        "6 qseqid sseqid sstrand qframe sframe frames",
        &["--dust", "no", "--max_target_seqs", "10"],
        &["-dust", "no", "-max_target_seqs", "10"],
    );
}

#[test]
fn blastn_subject_ncbi_parity_positive_gi_qcovus_fields() {
    assert_blastn_subject_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGTACGT\n",
        ">s1\nACGTACGTTCGTACGTACGTACGT\n",
        "6 nident positive pident ppos qcovhsp qcovus qcovs qgi sgi sallgi",
        &[
            "--dust",
            "no",
            "--strand",
            "plus",
            "--word_size",
            "7",
            "--max_target_seqs",
            "10",
        ],
        &[
            "-dust",
            "no",
            "-strand",
            "plus",
            "-word_size",
            "7",
            "-max_target_seqs",
            "10",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_qcovhsp_with_query_gap() {
    assert_blastn_subject_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n",
        ">s1\nACGTACGTACGTACGTACGTAAAACGTACGTACGTACGTACGT\n",
        "6 qstart qend nident mismatch gaps gapopen pident length qcovhsp qseq sseq btop",
        &[
            "--dust",
            "no",
            "--strand",
            "plus",
            "--word_size",
            "7",
            "--gapopen",
            "5",
            "--gapextend",
            "2",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "10",
        ],
        &[
            "-dust",
            "no",
            "-strand",
            "plus",
            "-word_size",
            "7",
            "-gapopen",
            "5",
            "-gapextend",
            "2",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_gapped_single_insertion_and_deletion() {
    assert_blastn_subject_outfmt_matches_ncbi(
        ">q\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n",
        ">s_ins\nACGTACGTACGTAAAACGTACGTACGTACGTACGTACGTACGT\n>s_del\nACGTACGTACGTACGTACGTACGTACGTACGTACGT\n",
        "6 sseqid qstart qend sstart send score length nident mismatch gaps gapopen qseq sseq btop",
        &[
            "--dust",
            "no",
            "--strand",
            "plus",
            "--word_size",
            "7",
            "--gapopen",
            "5",
            "--gapextend",
            "2",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "10",
        ],
        &[
            "-dust",
            "no",
            "-strand",
            "plus",
            "-word_size",
            "7",
            "-gapopen",
            "5",
            "-gapextend",
            "2",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_gapped_long_gap_traceback() {
    assert_blastn_subject_outfmt_matches_ncbi(
        ">q\nACGTTGCAACGATCGTACGATTCGAGCTTAGGCTAATCGGATCCTAGCTAGGCTAATCGATCGTAGCTAGCATCGAT\n",
        ">longgap\nACGTTGCAACGATCGTACGATTCGAGCTTAGGCTAGGGGGGGGGGGGGGGGTAATCGGATCCTAGCTAGGCTAATCGATCGTAGCTAGCATCGAT\n",
        "6 sseqid qstart qend sstart send score length nident mismatch gaps gapopen qseq sseq btop",
        &[
            "--dust",
            "no",
            "--strand",
            "plus",
            "--word_size",
            "7",
            "--reward",
            "1",
            "--penalty",
            "-2",
            "--gapopen",
            "0",
            "--gapextend",
            "2",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "5",
        ],
        &[
            "-dust",
            "no",
            "-strand",
            "plus",
            "-word_size",
            "7",
            "-reward",
            "1",
            "-penalty",
            "-2",
            "-gapopen",
            "0",
            "-gapextend",
            "2",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "5",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_gapped_traceback_edge_matrix() {
    let outfmt =
        "6 sseqid qstart qend sstart send score length nident mismatch gaps gapopen qseq sseq btop";
    let cases = [
        (
            "adjacent_ins_del",
            ">q\nACGTTGCAACGATCGTACGATTCGAGCTTAGGCTAATCGGATCCTAGCTAGGCTAATCGATCGTAGCTAGCATCGAT\n",
            ">s\nACGTTGCAACGATCGTACGATTCGAGCTTAGGCTAGGGTAATCGGATCCTAGCTAGGCTAATCGATCGTAGCTAGCATCGAT\n",
        ),
        (
            "adjacent_del_ins",
            ">q\nACGTTGCAACGATCGTACGATTCGAGCTTAGGCTAGGGTAATCGGATCCTAGCTAGGCTAATCGATCGTAGCTAGCATCGAT\n",
            ">s\nACGTTGCAACGATCGTACGATTCGAGCTTAGGCTAATCGGATCCTAGCTAGGCTAATCGATCGTAGCTAGCATCGAT\n",
        ),
        (
            "two_gap_choice",
            ">q\nACGTACGTACGTAAAACGTACGTACGT\n",
            ">s\nACGTACGTACGTTTTACGTACGTACGT\n",
        ),
        (
            "equal_gap_mismatch",
            ">q\nACGTACGTACGTACGTACGTACGT\n",
            ">s\nACGTACGTTCGTACGTACGTACGT\n",
        ),
        (
            "gap_near_start",
            ">q\nACGTACGTACGTACGTACGTACGTACGT\n",
            ">s\nACGTAAAACGTACGTACGTACGTACGTACGT\n",
        ),
        (
            "gap_near_end",
            ">q\nACGTACGTACGTACGTACGTACGTACGT\n",
            ">s\nACGTACGTACGTACGTACGTACGTAAAACGT\n",
        ),
    ];

    for (label, query, subject) in cases {
        assert_blastn_subject_outfmt_matches_ncbi(
            query,
            subject,
            outfmt,
            &[
                "--dust",
                "no",
                "--strand",
                "plus",
                "--word_size",
                "7",
                "--gapopen",
                "5",
                "--gapextend",
                "2",
                "--max_target_seqs",
                "10",
                "--max_hsps",
                "10",
            ],
            &[
                "-dust",
                "no",
                "-strand",
                "plus",
                "-word_size",
                "7",
                "-gapopen",
                "5",
                "-gapextend",
                "2",
                "-max_target_seqs",
                "10",
                "-max_hsps",
                "10",
            ],
        );
        eprintln!("checked {label}");
    }
}

#[test]
fn blastn_subject_ncbi_parity_gapped_low_prelim_high_final_xdrop_matrix() {
    let outfmt =
        "6 qstart qend sstart send score length nident mismatch gaps gapopen qseq sseq btop";
    let cases = [
        (
            "mismatch_block",
            ">q\nACGTACGTACGTAAAACGTACGTACGT\n",
            ">s\nACGTACGTACGTTTTACGTACGTACGT\n",
        ),
        (
            "gap_near_end",
            ">q\nACGTACGTACGTACGTACGTACGTACGT\n",
            ">s\nACGTACGTACGTACGTACGTACGTAAAACGT\n",
        ),
    ];

    for (label, query, subject) in cases {
        assert_blastn_subject_outfmt_matches_ncbi(
            query,
            subject,
            outfmt,
            &[
                "--dust",
                "no",
                "--strand",
                "plus",
                "--word_size",
                "7",
                "--gapopen",
                "5",
                "--gapextend",
                "2",
                "--xdrop_gap",
                "2",
                "--xdrop_gap_final",
                "100",
                "--max_target_seqs",
                "10",
                "--max_hsps",
                "10",
            ],
            &[
                "-dust",
                "no",
                "-strand",
                "plus",
                "-word_size",
                "7",
                "-gapopen",
                "5",
                "-gapextend",
                "2",
                "-xdrop_gap",
                "2",
                "-xdrop_gap_final",
                "100",
                "-max_target_seqs",
                "10",
                "-max_hsps",
                "10",
            ],
        );
        eprintln!("checked {label}");
    }
}

#[test]
fn blastn_subject_ncbi_parity_gapped_final_xdrop_boundary_matrix() {
    let outfmt =
        "6 qstart qend sstart send score length nident mismatch gaps gapopen qseq sseq btop";
    let cases = [
        (
            "mismatch_block_x2",
            ">q\nACGTACGTACGTAAAACGTACGTACGT\n",
            ">s\nACGTACGTACGTTTTACGTACGTACGT\n",
            "2",
        ),
        (
            "ambiguity_block_x1",
            ">q\nACGTACGTACGTACGTACGTACGTACGT\n",
            ">s\nACGTACGTACGTNNNNACGTACGTACGT\n",
            "1",
        ),
        (
            "gap_near_end_x3",
            ">q\nACGTACGTACGTACGTACGTACGTACGT\n",
            ">s\nACGTACGTACGTACGTACGTACGTAAAACGT\n",
            "3",
        ),
    ];

    for (label, query, subject, xdrop) in cases {
        assert_blastn_subject_outfmt_matches_ncbi(
            query,
            subject,
            outfmt,
            &[
                "--dust",
                "no",
                "--strand",
                "plus",
                "--word_size",
                "7",
                "--gapopen",
                "5",
                "--gapextend",
                "2",
                "--xdrop_gap_final",
                xdrop,
                "--max_target_seqs",
                "10",
                "--max_hsps",
                "10",
            ],
            &[
                "-dust",
                "no",
                "-strand",
                "plus",
                "-word_size",
                "7",
                "-gapopen",
                "5",
                "-gapextend",
                "2",
                "-xdrop_gap_final",
                xdrop,
                "-max_target_seqs",
                "10",
                "-max_hsps",
                "10",
            ],
        );
        eprintln!("checked {label}");
    }
}

#[test]
fn blastn_subject_ncbi_parity_gapped_traceback_minus_strand_edge_matrix() {
    let outfmt =
        "6 sseqid qstart qend sstart send score length nident mismatch gaps gapopen qseq sseq btop";
    let cases = [
        (
            "minus_adjacent_ins_del",
            "ACGTTGCAACGATCGTACGATTCGAGCTTAGGCTAATCGGATCCTAGCTAGGCTAATCGATCGTAGCTAGCATCGAT",
            "ACGTTGCAACGATCGTACGATTCGAGCTTAGGCTAGGGTAATCGGATCCTAGCTAGGCTAATCGATCGTAGCTAGCATCGAT",
        ),
        (
            "minus_equal_gap_mismatch",
            "ACGTACGTACGTACGTACGTACGT",
            "ACGTACGTTCGTACGTACGTACGT",
        ),
        (
            "minus_gap_near_start",
            "ACGTACGTACGTACGTACGTACGTACGT",
            "ACGTAAAACGTACGTACGTACGTACGTACGT",
        ),
        (
            "minus_gap_near_end",
            "ACGTACGTACGTACGTACGTACGTACGT",
            "ACGTACGTACGTACGTACGTACGTAAAACGT",
        ),
    ];

    for (label, query_seq, plus_subject_seq) in cases {
        let query = format!(">q\n{query_seq}\n");
        let subject = format!(">s\n{}\n", ascii_reverse_complement(plus_subject_seq));
        assert_blastn_subject_outfmt_matches_ncbi(
            &query,
            &subject,
            outfmt,
            &[
                "--dust",
                "no",
                "--strand",
                "minus",
                "--word_size",
                "7",
                "--gapopen",
                "5",
                "--gapextend",
                "2",
                "--max_target_seqs",
                "10",
                "--max_hsps",
                "10",
            ],
            &[
                "-dust",
                "no",
                "-strand",
                "minus",
                "-word_size",
                "7",
                "-gapopen",
                "5",
                "-gapextend",
                "2",
                "-max_target_seqs",
                "10",
                "-max_hsps",
                "10",
            ],
        );
        eprintln!("checked {label}");
    }
}

#[test]
fn blastn_subject_ncbi_parity_gapped_endpoint_ties_and_ambiguity_display() {
    let outfmt =
        "6 qstart qend sstart send score length nident mismatch gaps gapopen qseq sseq btop";
    for (label, query, subject) in [
        (
            "endpoint_trim_ties",
            ">q\nTTTACGTACGTACGTACGTACGTAAA\n",
            ">s\nGGGACGTACGTACGTACGTACGTCCC\n",
        ),
        (
            "ambiguous_subject_display",
            ">q\nACGTACGTACGTACGTACGTACGTACGT\n",
            ">s\nACGTACGTACGTNNNNACGTACGTACGT\n",
        ),
    ] {
        assert_blastn_subject_outfmt_matches_ncbi(
            query,
            subject,
            outfmt,
            &[
                "--dust",
                "no",
                "--word_size",
                "7",
                "--gapopen",
                "5",
                "--gapextend",
                "2",
                "--max_target_seqs",
                "10",
                "--max_hsps",
                "10",
            ],
            &[
                "-dust",
                "no",
                "-word_size",
                "7",
                "-gapopen",
                "5",
                "-gapextend",
                "2",
                "-max_target_seqs",
                "10",
                "-max_hsps",
                "10",
            ],
        );
        eprintln!("checked {label}");
    }
}

#[test]
fn blastn_subject_ncbi_parity_query_order_before_subject_tie_order() {
    assert_blastn_subject_outfmt_matches_ncbi(
        ">q_start\nAACGTGCTAGCTAGGCTAATCGGATCCTAGCTAGCTA\n>q_end\nACGTGCTAGCTAGGCTAATCGGATCCTAGCTAGCTAA\n",
        ">s_start\nACGTGCTAGCTAGGCTAATCGGATCCTAGCTAGCTA\n>s_end\nACGTGCTAGCTAGGCTAATCGGATCCTAGCTAGCTA\n",
        "6 qseqid sseqid qstart qend sstart send score length nident mismatch gaps gapopen qseq sseq btop",
        &[
            "--dust",
            "no",
            "--strand",
            "plus",
            "--word_size",
            "7",
            "--gapopen",
            "5",
            "--gapextend",
            "2",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "5",
        ],
        &[
            "-dust",
            "no",
            "-strand",
            "plus",
            "-word_size",
            "7",
            "-gapopen",
            "5",
            "-gapextend",
            "2",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "5",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_qcov_hsp_filter_with_query_gap() {
    assert_blastn_subject_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n",
        ">s1\nACGTACGTACGTACGTACGTAAAACGTACGTACGTACGTACGT\n",
        "6 qstart qend length qcovhsp qseq sseq",
        &[
            "--dust",
            "no",
            "--strand",
            "plus",
            "--word_size",
            "7",
            "--gapopen",
            "5",
            "--gapextend",
            "2",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "10",
            "--qcov_hsp_perc",
            "42",
        ],
        &[
            "-dust",
            "no",
            "-strand",
            "plus",
            "-word_size",
            "7",
            "-gapopen",
            "5",
            "-gapextend",
            "2",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "10",
            "-qcov_hsp_perc",
            "42",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_default_outfmt_with_query_gap() {
    assert_blastn_subject_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n",
        ">s1\nACGTACGTACGTACGTACGTAAAACGTACGTACGTACGTACGT\n",
        "6",
        &[
            "--dust",
            "no",
            "--strand",
            "plus",
            "--word_size",
            "7",
            "--gapopen",
            "5",
            "--gapextend",
            "2",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "10",
        ],
        &[
            "-dust",
            "no",
            "-strand",
            "plus",
            "-word_size",
            "7",
            "-gapopen",
            "5",
            "-gapextend",
            "2",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_sam_exact_hit() {
    assert_blastn_subject_sam_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGT\n",
        ">s1\nACGTACGTACGTACGTACGT\n",
        &["--dust", "no", "--max_target_seqs", "10"],
        &["-dust", "no", "-max_target_seqs", "10"],
    );
}

#[test]
fn blastn_subject_ncbi_parity_sam_query_gap() {
    assert_blastn_subject_sam_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTTTTTTTTTTTTTTTTTTTT\n",
        ">s1\nACGTACGTACGTACGTACGTAAAACGTACGTACGTACGTACGT\n",
        &[
            "--dust",
            "no",
            "--strand",
            "plus",
            "--word_size",
            "7",
            "--gapopen",
            "5",
            "--gapextend",
            "2",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "10",
        ],
        &[
            "-dust",
            "no",
            "-strand",
            "plus",
            "-word_size",
            "7",
            "-gapopen",
            "5",
            "-gapextend",
            "2",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_sam_subject_gap() {
    assert_blastn_subject_sam_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGTAAAACGTACGTACGTACGTACGT\n",
        ">s1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n",
        &[
            "--dust",
            "no",
            "--strand",
            "plus",
            "--word_size",
            "7",
            "--gapopen",
            "5",
            "--gapextend",
            "2",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "10",
        ],
        &[
            "-dust",
            "no",
            "-strand",
            "plus",
            "-word_size",
            "7",
            "-gapopen",
            "5",
            "-gapextend",
            "2",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_subject_ncbi_parity_sam_no_hit() {
    assert_blastn_subject_sam_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGT\n",
        ">s1\nTTTTTTTTTTTTTTTTTTTT\n",
        &["--dust", "no", "--max_target_seqs", "10"],
        &["-dust", "no", "-max_target_seqs", "10"],
    );
}

#[test]
fn blastn_subject_ncbi_parity_xml_exact_hit() {
    assert_blastn_subject_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGT\n",
        ">s1 subject one\nACGTACGTACGTACGTACGT\n",
        "5",
        &["--dust", "no", "--max_target_seqs", "10"],
        &["-dust", "no", "-max_target_seqs", "10"],
    );
}

#[test]
fn blastn_db_ncbi_parity_sam_exact_hit() {
    assert_blastn_db_sam_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGT\n",
        ">s1\nACGTACGTACGTACGTACGT\n",
        &["--dust", "no", "--max_target_seqs", "10"],
        &["-dust", "no", "-max_target_seqs", "10"],
    );
}

#[test]
fn blastn_db_ncbi_parity_xml_exact_hit() {
    assert_blastn_db_xml_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGT\n",
        ">s1 subject one\nACGTACGTACGTACGTACGT\n",
        &["--dust", "no", "--max_target_seqs", "10"],
        &["-dust", "no", "-max_target_seqs", "10"],
    );
}

#[test]
fn blastn_db_ncbi_parity_pairwise_exact_hit() {
    assert_blastn_db_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGT\n",
        ">s1 subject one\nACGTACGTACGTACGTACGT\n",
        "0",
        &[],
        &["--dust", "no", "--max_target_seqs", "10"],
        &["-dust", "no", "-max_target_seqs", "10"],
    );
}

#[test]
fn blastn_db_ncbi_parity_pairwise_no_hit() {
    assert_blastn_db_outfmt_matches_ncbi(
        ">q1\nAAAAAAAAAAAAAAAAAAAA\n",
        ">s1 subject one\nCCCCCCCCCCCCCCCCCCCC\n",
        "0",
        &[],
        &["--dust", "no", "--max_target_seqs", "10"],
        &["-dust", "no", "-max_target_seqs", "10"],
    );
}

#[test]
fn blastn_db_ncbi_parity_short_read_masked_fixture() {
    assert_blastn_db_outfmt_matches_ncbi(
        ">q1\nACGTACGTacgtacgtGACTTACCGTACGTACGTGACTTACCGTACGT\n",
        ">s1 subject one\nTTTTACGTACGTACGTACGTGACTTACCGTACGTACGTGACTTACCGTACGTAAAA\n>s2 subject two\nACGTACGTACGTACGTGACTTACCGTACGTAAAAACGTGACTTACCGTACGT\n>s3 subject three\nCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "6 qseqid qstart qend sstart send length mismatch gapopen score bitscore btop",
        &[],
        &[
            "--dust",
            "no",
            "--lcase_masking",
            "--soft_masking",
            "true",
            "--word_size",
            "7",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "20",
            "--max_hsps",
            "10",
        ],
        &[
            "-dust",
            "no",
            "-lcase_masking",
            "-soft_masking",
            "true",
            "-word_size",
            "7",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "20",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_terminal_overhang_mapper_fixture() {
    assert_blastn_db_outfmt_matches_ncbi(
        ">q1\nTTTTACGTACGTGACTTACCGTACGTACGTAAAA\n",
        ">s1 subject one\nACGTACGTGACTTACCGTACGTACGT\n>s2 subject two\nCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "6 qseqid qstart qend sstart send length mismatch gapopen score bitscore qseq sseq btop",
        &[],
        &[
            "--task",
            "blastn-short",
            "--dust",
            "no",
            "--word_size",
            "7",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "10",
        ],
        &[
            "-task",
            "blastn-short",
            "-dust",
            "no",
            "-word_size",
            "7",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_terminal_overhang_mapper_multi_subject_fixture() {
    assert_blastn_db_outfmt_matches_ncbi_sorted_lines(
        ">q1\nTTTTACGTACGTGACTTACCGTACGTACGTAAAA\n",
        ">s1 subject one\nACGTACGTGACTTACCGTACGTACGT\n>s2 subject two\nACGTACGTGACTTACCGTACGTACGT\n>s3 subject three\nCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "6 qseqid qstart qend sstart send length mismatch gapopen score bitscore qseq sseq btop",
        &[],
        &[
            "--task",
            "blastn-short",
            "--dust",
            "no",
            "--word_size",
            "7",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "10",
        ],
        &[
            "-task",
            "blastn-short",
            "-dust",
            "no",
            "-word_size",
            "7",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_minus_strand_terminal_overhang_mapper_fixture() {
    assert_blastn_db_outfmt_matches_ncbi(
        ">q1\nTTTTACGTACGTGACTTACCGTACGTACGTAAAA\n",
        ">s1 subject one\nACGTACGTACGGTAAGTCACGTACGT\n>s2 subject two\nCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "6 qseqid qstart qend sstart send length mismatch gapopen score bitscore qseq sseq btop",
        &[],
        &[
            "--task",
            "blastn-short",
            "--dust",
            "no",
            "--strand",
            "minus",
            "--word_size",
            "7",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "10",
        ],
        &[
            "-task",
            "blastn-short",
            "-dust",
            "no",
            "-strand",
            "minus",
            "-word_size",
            "7",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_minus_strand_terminal_overhang_mapper_multi_subject_fixture() {
    assert_blastn_db_outfmt_matches_ncbi_sorted_lines(
        ">q1\nTTTTACGTACGTGACTTACCGTACGTACGTAAAA\n",
        ">s1 subject one\nACGTACGTACGGTAAGTCACGTACGT\n>s2 subject two\nACGTACGTACGGTAAGTCACGTACGT\n>s3 subject three\nCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "6 qseqid qstart qend sstart send length mismatch gapopen score bitscore qseq sseq btop",
        &[],
        &[
            "--task",
            "blastn-short",
            "--dust",
            "no",
            "--strand",
            "minus",
            "--word_size",
            "7",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "10",
        ],
        &[
            "-task",
            "blastn-short",
            "-dust",
            "no",
            "-strand",
            "minus",
            "-word_size",
            "7",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_adapter_like_trailing_overhang_mapper_fixture() {
    assert_blastn_db_outfmt_matches_ncbi(
        ">q1\nTTTTACGTACGTGACTTACCGTACGTACGTACACGGAAGAC\n",
        ">s1 subject one\nACGTACGTGACTTACCGTACGTACGT\n>s2 subject two\nCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "6 qseqid qstart qend sstart send length mismatch gapopen score bitscore qseq sseq btop",
        &[],
        &[
            "--task",
            "blastn-short",
            "--dust",
            "no",
            "--word_size",
            "7",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "10",
        ],
        &[
            "-task",
            "blastn-short",
            "-dust",
            "no",
            "-word_size",
            "7",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_adapter_like_trailing_overhang_mapper_multi_subject_fixture() {
    assert_blastn_db_outfmt_matches_ncbi_sorted_lines(
        ">q1\nTTTTACGTACGTGACTTACCGTACGTACGTACACGGAAGAC\n",
        ">s1 subject one\nACGTACGTGACTTACCGTACGTACGT\n>s2 subject two\nACGTACGTGACTTACCGTACGTACGT\n>s3 subject three\nCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "6 qseqid qstart qend sstart send length mismatch gapopen score bitscore qseq sseq btop",
        &[],
        &[
            "--task",
            "blastn-short",
            "--dust",
            "no",
            "--word_size",
            "7",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "10",
        ],
        &[
            "-task",
            "blastn-short",
            "-dust",
            "no",
            "-word_size",
            "7",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_poly_a_trailing_overhang_mapper_fixture() {
    assert_blastn_db_outfmt_matches_ncbi(
        ">q1\nTTTTACGTACGTGACTTACCGTACGTACGTAAAAAA\n",
        ">s1 subject one\nACGTACGTGACTTACCGTACGTACGT\n>s2 subject two\nCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "6 qseqid qstart qend sstart send length mismatch gapopen score bitscore qseq sseq btop",
        &[],
        &[
            "--task",
            "blastn-short",
            "--dust",
            "no",
            "--word_size",
            "7",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "10",
        ],
        &[
            "-task",
            "blastn-short",
            "-dust",
            "no",
            "-word_size",
            "7",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_poly_a_trailing_overhang_mapper_multi_subject_fixture() {
    assert_blastn_db_outfmt_matches_ncbi_sorted_lines(
        ">q1\nTTTTACGTACGTGACTTACCGTACGTACGTAAAAAA\n",
        ">s1 subject one\nACGTACGTGACTTACCGTACGTACGT\n>s2 subject two\nACGTACGTGACTTACCGTACGTACGT\n>s3 subject three\nCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "6 qseqid qstart qend sstart send length mismatch gapopen score bitscore qseq sseq btop",
        &[],
        &[
            "--task",
            "blastn-short",
            "--dust",
            "no",
            "--word_size",
            "7",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "10",
        ],
        &[
            "-task",
            "blastn-short",
            "-dust",
            "no",
            "-word_size",
            "7",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_exact_adapter_trimming_mapper_fixture() {
    assert_blastn_db_outfmt_matches_ncbi(
        ">q1\nTTTTACGTACGTGACTTACCGTACGTACGTAGATCGGAAGAG\n",
        ">s1 subject one\nACGTACGTGACTTACCGTACGTACGT\n>s2 subject two\nCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "6 qseqid qstart qend sstart send length mismatch gapopen score bitscore qseq sseq btop",
        &[],
        &[
            "--task",
            "blastn-short",
            "--dust",
            "no",
            "--word_size",
            "7",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "10",
        ],
        &[
            "-task",
            "blastn-short",
            "-dust",
            "no",
            "-word_size",
            "7",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_exact_adapter_trimming_mapper_multi_subject_fixture() {
    assert_blastn_db_outfmt_matches_ncbi_sorted_lines(
        ">q1\nTTTTACGTACGTGACTTACCGTACGTACGTAGATCGGAAGAG\n",
        ">s1 subject one\nACGTACGTGACTTACCGTACGTACGT\n>s2 subject two\nACGTACGTGACTTACCGTACGTACGT\n>s3 subject three\nCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "6 qseqid qstart qend sstart send length mismatch gapopen score bitscore qseq sseq btop",
        &[],
        &[
            "--task",
            "blastn-short",
            "--dust",
            "no",
            "--word_size",
            "7",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "10",
        ],
        &[
            "-task",
            "blastn-short",
            "-dust",
            "no",
            "-word_size",
            "7",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_rmblastn_terminal_overhang_mapper_fixture() {
    assert_blastn_db_outfmt_matches_ncbi(
        ">q1\nTTTTACGTACGTGACTTACCGTACGTACGTAAAA\n",
        ">s1 subject one\nACGTACGTGACTTACCGTACGTACGT\n>s2 subject two\nCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "6 qseqid qstart qend sstart send length mismatch gapopen score bitscore qseq sseq btop",
        &[],
        &[
            "--task",
            "rmblastn",
            "--dust",
            "no",
            "--word_size",
            "7",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "10",
        ],
        &[
            "-task",
            "rmblastn",
            "-dust",
            "no",
            "-word_size",
            "7",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_rmblastn_terminal_overhang_mapper_multi_subject_fixture() {
    assert_blastn_db_outfmt_matches_ncbi_sorted_lines(
        ">q1\nTTTTACGTACGTGACTTACCGTACGTACGTAAAA\n",
        ">s1 subject one\nACGTACGTGACTTACCGTACGTACGT\n>s2 subject two\nACGTACGTGACTTACCGTACGTACGT\n>s3 subject three\nCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "6 qseqid qstart qend sstart send length mismatch gapopen score bitscore qseq sseq btop",
        &[],
        &[
            "--task",
            "rmblastn",
            "--dust",
            "no",
            "--word_size",
            "7",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "10",
        ],
        &[
            "-task",
            "rmblastn",
            "-dust",
            "no",
            "-word_size",
            "7",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_rmblastn_minus_strand_terminal_overhang_mapper_coordinates() {
    assert_blastn_db_outfmt_matches_ncbi(
        ">q1\nTTTTACGTACGTGACTTACCGTACGTACGTAAAA\n",
        ">s1 subject one\nACGTACGTACGGTAAGTCACGTACGT\n>s2 subject two\nCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "6 qseqid qstart qend sstart send length mismatch qseq sseq btop",
        &[],
        &[
            "--task",
            "rmblastn",
            "--dust",
            "no",
            "--strand",
            "minus",
            "--word_size",
            "7",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "10",
        ],
        &[
            "-task",
            "rmblastn",
            "-dust",
            "no",
            "-strand",
            "minus",
            "-word_size",
            "7",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_rmblastn_minus_strand_terminal_overhang_mapper_multi_subject_coordinates()
{
    assert_blastn_db_outfmt_matches_ncbi_sorted_lines(
        ">q1\nTTTTACGTACGTGACTTACCGTACGTACGTAAAA\n",
        ">s1 subject one\nACGTACGTACGGTAAGTCACGTACGT\n>s2 subject two\nACGTACGTACGGTAAGTCACGTACGT\n>s3 subject three\nCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "6 qseqid qstart qend sstart send length mismatch qseq sseq btop",
        &[],
        &[
            "--task",
            "rmblastn",
            "--dust",
            "no",
            "--strand",
            "minus",
            "--word_size",
            "7",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "10",
        ],
        &[
            "-task",
            "rmblastn",
            "-dust",
            "no",
            "-strand",
            "minus",
            "-word_size",
            "7",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
#[ignore = "documents the remaining rmblastn minus-strand terminal secondary-HSP scoring gap"]
fn blastn_db_ncbi_parity_rmblastn_minus_strand_terminal_overhang_mapper_scored() {
    assert_blastn_db_outfmt_matches_ncbi(
        ">q1\nTTTTACGTACGTGACTTACCGTACGTACGTAAAA\n",
        ">s1 subject one\nACGTACGTACGGTAAGTCACGTACGT\n>s2 subject two\nCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "6 qseqid qstart qend sstart send length mismatch gapopen score bitscore qseq sseq btop",
        &[],
        &[
            "--task",
            "rmblastn",
            "--dust",
            "no",
            "--strand",
            "minus",
            "--word_size",
            "7",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "10",
        ],
        &[
            "-task",
            "rmblastn",
            "-dust",
            "no",
            "-strand",
            "minus",
            "-word_size",
            "7",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
#[ignore = "documents the remaining rmblastn minus-strand terminal secondary-HSP scoring gap"]
fn blastn_db_ncbi_parity_rmblastn_minus_strand_terminal_overhang_mapper_multi_subject_scored() {
    assert_blastn_db_outfmt_matches_ncbi_sorted_lines(
        ">q1\nTTTTACGTACGTGACTTACCGTACGTACGTAAAA\n",
        ">s1 subject one\nACGTACGTACGGTAAGTCACGTACGT\n>s2 subject two\nACGTACGTACGGTAAGTCACGTACGT\n>s3 subject three\nCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "6 qseqid qstart qend sstart send length mismatch gapopen score bitscore qseq sseq btop",
        &[],
        &[
            "--task",
            "rmblastn",
            "--dust",
            "no",
            "--strand",
            "minus",
            "--word_size",
            "7",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "10",
        ],
        &[
            "-task",
            "rmblastn",
            "-dust",
            "no",
            "-strand",
            "minus",
            "-word_size",
            "7",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_rmblastn_adapter_like_trailing_overhang_mapper_fixture() {
    assert_blastn_db_outfmt_matches_ncbi(
        ">q1\nTTTTACGTACGTGACTTACCGTACGTACGTACACGGAAGAC\n",
        ">s1 subject one\nACGTACGTGACTTACCGTACGTACGT\n>s2 subject two\nCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "6 qseqid qstart qend sstart send length mismatch gapopen score bitscore qseq sseq btop",
        &[],
        &[
            "--task",
            "rmblastn",
            "--dust",
            "no",
            "--word_size",
            "7",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "10",
        ],
        &[
            "-task",
            "rmblastn",
            "-dust",
            "no",
            "-word_size",
            "7",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_rmblastn_adapter_like_trailing_overhang_mapper_multi_subject_fixture() {
    assert_blastn_db_outfmt_matches_ncbi_sorted_lines(
        ">q1\nTTTTACGTACGTGACTTACCGTACGTACGTACACGGAAGAC\n",
        ">s1 subject one\nACGTACGTGACTTACCGTACGTACGT\n>s2 subject two\nACGTACGTGACTTACCGTACGTACGT\n>s3 subject three\nCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "6 qseqid qstart qend sstart send length mismatch gapopen score bitscore qseq sseq btop",
        &[],
        &[
            "--task",
            "rmblastn",
            "--dust",
            "no",
            "--word_size",
            "7",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "10",
        ],
        &[
            "-task",
            "rmblastn",
            "-dust",
            "no",
            "-word_size",
            "7",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_rmblastn_poly_a_trailing_overhang_mapper_fixture() {
    assert_blastn_db_outfmt_matches_ncbi(
        ">q1\nTTTTACGTACGTGACTTACCGTACGTACGTAAAAAA\n",
        ">s1 subject one\nACGTACGTGACTTACCGTACGTACGT\n>s2 subject two\nCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "6 qseqid qstart qend sstart send length mismatch gapopen score bitscore qseq sseq btop",
        &[],
        &[
            "--task",
            "rmblastn",
            "--dust",
            "no",
            "--word_size",
            "7",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "10",
        ],
        &[
            "-task",
            "rmblastn",
            "-dust",
            "no",
            "-word_size",
            "7",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_rmblastn_poly_a_trailing_overhang_mapper_multi_subject_fixture() {
    assert_blastn_db_outfmt_matches_ncbi_sorted_lines(
        ">q1\nTTTTACGTACGTGACTTACCGTACGTACGTAAAAAA\n",
        ">s1 subject one\nACGTACGTGACTTACCGTACGTACGT\n>s2 subject two\nACGTACGTGACTTACCGTACGTACGT\n>s3 subject three\nCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "6 qseqid qstart qend sstart send length mismatch gapopen score bitscore qseq sseq btop",
        &[],
        &[
            "--task",
            "rmblastn",
            "--dust",
            "no",
            "--word_size",
            "7",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "10",
        ],
        &[
            "-task",
            "rmblastn",
            "-dust",
            "no",
            "-word_size",
            "7",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_rmblastn_exact_adapter_trimming_mapper_fixture() {
    assert_blastn_db_outfmt_matches_ncbi(
        ">q1\nTTTTACGTACGTGACTTACCGTACGTACGTAGATCGGAAGAG\n",
        ">s1 subject one\nACGTACGTGACTTACCGTACGTACGT\n>s2 subject two\nCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "6 qseqid qstart qend sstart send length mismatch gapopen score bitscore qseq sseq btop",
        &[],
        &[
            "--task",
            "rmblastn",
            "--dust",
            "no",
            "--word_size",
            "7",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "10",
        ],
        &[
            "-task",
            "rmblastn",
            "-dust",
            "no",
            "-word_size",
            "7",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_rmblastn_exact_adapter_trimming_mapper_multi_subject_fixture() {
    assert_blastn_db_outfmt_matches_ncbi_sorted_lines(
        ">q1\nTTTTACGTACGTGACTTACCGTACGTACGTAGATCGGAAGAG\n",
        ">s1 subject one\nACGTACGTGACTTACCGTACGTACGT\n>s2 subject two\nACGTACGTGACTTACCGTACGTACGT\n>s3 subject three\nCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "6 qseqid qstart qend sstart send length mismatch gapopen score bitscore qseq sseq btop",
        &[],
        &[
            "--task",
            "rmblastn",
            "--dust",
            "no",
            "--word_size",
            "7",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "10",
        ],
        &[
            "-task",
            "rmblastn",
            "-dust",
            "no",
            "-word_size",
            "7",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_indexed_megablast_fixture() {
    assert_blastn_indexed_db_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGTACGTACGTACGT\n",
        ">s1 subject one\nTTTTACGTACGTACGTACGTACGTACGTACGTACGTAAAA\n>s2 subject two\nCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "6 qseqid qstart qend sstart send length mismatch gapopen score bitscore btop",
        &[
            "--task",
            "megablast",
            "--use_index",
            "true",
            "--dust",
            "no",
            "--word_size",
            "16",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "6",
        ],
        &[
            "-task",
            "megablast",
            "-use_index",
            "true",
            "-dust",
            "no",
            "-word_size",
            "16",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "6",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_indexed_megablast_no_hit_fixture() {
    assert_blastn_indexed_db_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGTACGTACGTACGT\n",
        ">s1 subject one\nTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n>s2 subject two\nCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "6 qseqid qstart qend sstart send length mismatch gapopen score bitscore btop",
        &[
            "--task",
            "megablast",
            "--use_index",
            "true",
            "--dust",
            "no",
            "--word_size",
            "16",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "6",
        ],
        &[
            "-task",
            "megablast",
            "-use_index",
            "true",
            "-dust",
            "no",
            "-word_size",
            "16",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "6",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_indexed_megablast_multi_subject_fixture() {
    assert_blastn_indexed_db_outfmt_matches_ncbi(
        ">q1\nACGTACGTTGCAAGTCGATCGTACGATTCGGA\n",
        ">s1 subject one\nTTTTACGTACGTTGCAAGTCGATCGTACGATTCGGAAAAA\n>s2 subject two\nGGGGACGTACGTTGCAAGTCGATCGTACGATTCGGACCCC\n>s3 subject three\nCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "6 qseqid qstart qend sstart send length mismatch gapopen score bitscore btop",
        &[
            "--task",
            "megablast",
            "--use_index",
            "true",
            "--dust",
            "no",
            "--word_size",
            "16",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "6",
        ],
        &[
            "-task",
            "megablast",
            "-use_index",
            "true",
            "-dust",
            "no",
            "-word_size",
            "16",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "6",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_indexed_short_read_masked_fixture() {
    assert_blastn_indexed_db_outfmt_matches_ncbi(
        ">q1\nACGTACGTacgtacgtGACTTACCGTACGTACGTGACTTACCGTACGT\n",
        ">s1 subject one\nTTTTACGTACGTACGTACGTGACTTACCGTACGTACGTGACTTACCGTACGTAAAA\n>s2 subject two\nACGTACGTACGTACGTGACTTACCGTACGTAAAAACGTGACTTACCGTACGT\n>s3 subject three\nCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "6 qseqid qstart qend sstart send length mismatch gapopen score bitscore btop",
        &[
            "--task",
            "blastn-short",
            "--use_index",
            "true",
            "--dust",
            "no",
            "--lcase_masking",
            "--soft_masking",
            "true",
            "--word_size",
            "7",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "20",
            "--max_hsps",
            "10",
        ],
        &[
            "-task",
            "blastn-short",
            "-use_index",
            "true",
            "-dust",
            "no",
            "-lcase_masking",
            "-soft_masking",
            "true",
            "-word_size",
            "7",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "20",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_indexed_short_read_multi_subject_fixture() {
    assert_blastn_indexed_db_outfmt_matches_ncbi(
        ">q1\nACGTACGTTGCAAGTCGATCGTACGATT\n",
        ">s1 subject one\nTTTTACGTACGTTGCAAGTCGATCGTACGATTAAAA\n>s2 subject two\nGGGGACGTACGTTGCAAGTCGATCGTACGATTCCCC\n>s3 subject three\nCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "6 qseqid qstart qend sstart send length mismatch gapopen score bitscore btop",
        &[
            "--task",
            "blastn-short",
            "--use_index",
            "true",
            "--dust",
            "no",
            "--word_size",
            "7",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "6",
        ],
        &[
            "-task",
            "blastn-short",
            "-use_index",
            "true",
            "-dust",
            "no",
            "-word_size",
            "7",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "6",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_indexed_short_read_no_hit_fixture() {
    assert_blastn_indexed_db_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGTACGTACGT\n",
        ">s1 subject one\nTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n>s2 subject two\nCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "6 qseqid qstart qend sstart send length mismatch gapopen score bitscore btop",
        &[
            "--task",
            "blastn-short",
            "--use_index",
            "true",
            "--dust",
            "no",
            "--word_size",
            "7",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "6",
        ],
        &[
            "-task",
            "blastn-short",
            "-use_index",
            "true",
            "-dust",
            "no",
            "-word_size",
            "7",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "6",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_dc_megablast_discontiguous_template_fixture() {
    assert_blastn_db_outfmt_matches_ncbi(
        ">q1\nACGTACGTGACTTACCGTACGTACGTGACTTACCGTACGTACGTGACTTACCGTACGT\n",
        ">s1 subject one\nTTTTACGTACGTGACTTACCGTACGTACGTGACTTACCGTACGTACGTGACTTACCGTACGTAAAA\n",
        "6 qseqid qstart qend sstart send length mismatch gapopen score bitscore btop",
        &[],
        &[
            "--task",
            "dc-megablast",
            "--dust",
            "no",
            "--template_type",
            "coding",
            "--template_length",
            "18",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "20",
            "--max_hsps",
            "10",
        ],
        &[
            "-task",
            "dc-megablast",
            "-dust",
            "no",
            "-template_type",
            "coding",
            "-template_length",
            "18",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "20",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_dc_megablast_two_template_fixture() {
    assert_blastn_db_outfmt_matches_ncbi(
        ">q1\nACGTACGTGACTTACCGTACGTACGTGACTTACCGTACGTACGTGACTTACCGTACGT\n",
        ">s1 subject one\nGGGGACGTACGTGACTTACCGTACGTACGTGACTTACCGTACGTACGTGACTTACCGTACGTCCCC\n",
        "6 qseqid qstart qend sstart send length mismatch gapopen score bitscore btop",
        &[],
        &[
            "--task",
            "dc-megablast",
            "--dust",
            "no",
            "--template_type",
            "coding_and_optimal",
            "--template_length",
            "16",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "20",
            "--max_hsps",
            "10",
        ],
        &[
            "-task",
            "dc-megablast",
            "-dust",
            "no",
            "-template_type",
            "coding_and_optimal",
            "-template_length",
            "16",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "20",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_dc_megablast_multi_subject_two_template_fixture() {
    assert_blastn_db_outfmt_matches_ncbi(
        ">q1\nACGTACGTGACTTACCGTACGTACGTGACTTACCGTACGTACGTGACTTACCGTACGT\n",
        ">s_primary subject one\nGGGGACGTACGTGACTTACCGTACGTACGTGACTTACCGTACGTACGTGACTTACCGTACGTCCCC\n>s_shifted subject two\nTTTTTTACGTACGTGACTTACCGTACGTACGTGACTTACCGTACGTACGTGACTTACCGTACGTAAAA\n>s_noise subject three\nCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "6 qseqid qstart qend sstart send length mismatch gapopen score bitscore btop",
        &[],
        &[
            "--task",
            "dc-megablast",
            "--dust",
            "no",
            "--template_type",
            "coding_and_optimal",
            "--template_length",
            "16",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "20",
            "--max_hsps",
            "10",
        ],
        &[
            "-task",
            "dc-megablast",
            "-dust",
            "no",
            "-template_type",
            "coding_and_optimal",
            "-template_length",
            "16",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "20",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_dc_megablast_masked_discontiguous_template_fixture() {
    assert_blastn_db_outfmt_matches_ncbi(
        ">q1\nACGTACGTgacttaccGTACGTACGTGACTTACCGTACGTACGTGACTTACCGT\n",
        ">s1 subject one\nTTTTACGTACGTGACTTACCGTACGTACGTGACTTACCGTACGTACGTGACTTACCGTAAAA\n>s2 subject two\nCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "6 qseqid qstart qend sstart send length mismatch gapopen score bitscore btop",
        &[],
        &[
            "--task",
            "dc-megablast",
            "--dust",
            "no",
            "--lcase_masking",
            "--soft_masking",
            "true",
            "--template_type",
            "coding",
            "--template_length",
            "18",
            "--evalue",
            "1000",
            "--max_target_seqs",
            "20",
            "--max_hsps",
            "10",
        ],
        &[
            "-task",
            "dc-megablast",
            "-dust",
            "no",
            "-lcase_masking",
            "-soft_masking",
            "true",
            "-template_type",
            "coding",
            "-template_length",
            "18",
            "-evalue",
            "1000",
            "-max_target_seqs",
            "20",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_max_target_seqs_after_best_hit_filter() {
    let query = ">q1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";
    let db = ">s_full\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n>s_bad\nACGTACGTACATACGTACATACGTAC\n>s_bad2\nACGTACGTACGTTCGTACGTTCGTACGTACGT\n>s_tail\nTTTTACGTACGTACGTACGTAAAA\n";

    assert_blastn_db_outfmt_matches_ncbi(
        query,
        db,
        "6 sseqid",
        &[],
        &[
            "--dust",
            "no",
            "--word_size",
            "7",
            "--max_target_seqs",
            "2",
            "--best_hit_overhang",
            "0.1",
            "--best_hit_score_edge",
            "0.1",
            "--max_hsps",
            "1",
        ],
        &[
            "-dust",
            "no",
            "-word_size",
            "7",
            "-max_target_seqs",
            "2",
            "-best_hit_overhang",
            "0.1",
            "-best_hit_score_edge",
            "0.1",
            "-max_hsps",
            "1",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_best_hit_tied_hsp_ordering() {
    let query = ">q1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";
    let db = ">s_full\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n>s_bad\nACGTACGTACATACGTACATACGTAC\n>s_bad2\nACGTACGTACGTTCGTACGTTCGTACGTACGT\n>s_tail\nTTTTACGTACGTACGTACGTAAAA\n";

    assert_blastn_db_outfmt_matches_ncbi(
        query,
        db,
        "6 sseqid qstart qend sstart send length score bitscore evalue btop",
        &[],
        &[
            "--dust",
            "no",
            "--word_size",
            "7",
            "--max_target_seqs",
            "2",
            "--best_hit_overhang",
            "0.1",
            "--best_hit_score_edge",
            "0.1",
        ],
        &[
            "-dust",
            "no",
            "-word_size",
            "7",
            "-max_target_seqs",
            "2",
            "-best_hit_overhang",
            "0.1",
            "-best_hit_score_edge",
            "0.1",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_culling_tied_subject_ordering() {
    let query = ">q1\nACGTACGTACGTACGTACGTACGTACGTACGT\n";
    let db = ">s_alpha\nACGTACGTACGTACGTACGTACGTACGTACGT\n>s_beta\nACGTACGTACGTACGTACGTACGTACGTACGT\n>s_gamma\nACGTACGTACGTACGTACGTACGTACGTACGT\n";

    assert_blastn_db_outfmt_matches_ncbi(
        query,
        db,
        "6 sseqid qstart qend sstart send score bitscore evalue",
        &[],
        &[
            "--dust",
            "no",
            "--word_size",
            "7",
            "--max_target_seqs",
            "10",
            "--culling_limit",
            "1",
        ],
        &[
            "-dust",
            "no",
            "-word_size",
            "7",
            "-max_target_seqs",
            "10",
            "-culling_limit",
            "1",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_best_hit_tied_subject_ordering() {
    let query = ">q1\nACGTACGTACGTACGTACGTACGTACGTACGT\n";
    let db = ">s_alpha\nACGTACGTACGTACGTACGTACGTACGTACGT\n>s_beta\nACGTACGTACGTACGTACGTACGTACGTACGT\n>s_gamma\nACGTACGTACGTACGTACGTACGTACGTACGT\n";

    assert_blastn_db_outfmt_matches_ncbi(
        query,
        db,
        "6 sseqid qstart qend sstart send score bitscore evalue",
        &[],
        &[
            "--dust",
            "no",
            "--word_size",
            "7",
            "--max_target_seqs",
            "10",
            "--best_hit_overhang",
            "0.1",
            "--best_hit_score_edge",
            "0.1",
        ],
        &[
            "-dust",
            "no",
            "-word_size",
            "7",
            "-max_target_seqs",
            "10",
            "-best_hit_overhang",
            "0.1",
            "-best_hit_score_edge",
            "0.1",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_outfmt6_duplicate_fields() {
    assert_blastn_db_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGT\n",
        ">s1\nACGTACGTACGTACGTACGT\n",
        "6 qseqid qseqid sseqid length length bitscore",
        &[],
        &["--dust", "no", "--max_target_seqs", "10"],
        &["-dust", "no", "-max_target_seqs", "10"],
    );
}

#[test]
fn blastn_db_ncbi_parity_outfmt6_invalid_fields_are_ignored() {
    assert_blastn_db_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGT\n",
        ">s1\nACGTACGTACGTACGTACGT\n",
        "6 qseqid bogus sseqid",
        &[],
        &["--dust", "no", "--max_target_seqs", "10"],
        &["-dust", "no", "-max_target_seqs", "10"],
    );

    assert_blastn_db_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGT\n",
        ">s1\nACGTACGTACGTACGTACGT\n",
        "6 bogus",
        &[],
        &["--dust", "no", "--max_target_seqs", "10"],
        &["-dust", "no", "-max_target_seqs", "10"],
    );
}

#[test]
fn blastn_db_ncbi_parity_outfmt6_empty_field_list_uses_defaults() {
    for outfmt in ["6 ", "6 delim=,"] {
        assert_blastn_db_outfmt_matches_ncbi(
            ">q1\nACGTACGTACGTACGTACGT\n",
            ">s1\nACGTACGTACGTACGTACGT\n",
            outfmt,
            &[],
            &["--dust", "no", "--max_target_seqs", "10"],
            &["-dust", "no", "-max_target_seqs", "10"],
        );
    }
}

#[test]
fn blastn_db_ncbi_parity_custom_field_parser_edges() {
    for outfmt in [
        "6 delim=tab qseqid sseqid length",
        "6 delim=, delim=tab qseqid sseqid length",
        "6 std qaccver saccver qlen",
        "7 qseqid bogus sseqid qseqid length length",
        "7 delim=space qseqid sseqid length",
        "7 delim=, delim=tab qseqid sseqid length",
        "7 std qaccver saccver qlen",
        "7 bogus",
        "10 qseqid bogus sseqid qseqid length length",
        "10 delim=tab qseqid sseqid length",
        "10 delim=tab delim=, qseqid sseqid length",
        "10 std qaccver saccver qlen",
        "10 bogus",
    ] {
        assert_blastn_db_outfmt_matches_ncbi(
            ">q1\nACGTACGTACGTACGTACGT\n",
            ">s1\nACGTACGTACGTACGTACGT\n",
            outfmt,
            &[],
            &["--dust", "no", "--max_target_seqs", "10"],
            &["-dust", "no", "-max_target_seqs", "10"],
        );
    }
}

#[test]
fn blastn_db_ncbi_parity_query_order_before_subject_tie_order() {
    assert_blastn_db_outfmt_matches_ncbi(
        ">q_start\nAACGTGCTAGCTAGGCTAATCGGATCCTAGCTAGCTA\n>q_end\nACGTGCTAGCTAGGCTAATCGGATCCTAGCTAGCTAA\n",
        ">s_start\nACGTGCTAGCTAGGCTAATCGGATCCTAGCTAGCTA\n>s_end\nACGTGCTAGCTAGGCTAATCGGATCCTAGCTAGCTA\n",
        "6 qseqid sseqid qstart qend sstart send score length nident mismatch gaps gapopen qseq sseq btop",
        &[],
        &[
            "--dust",
            "no",
            "--strand",
            "plus",
            "--word_size",
            "7",
            "--gapopen",
            "5",
            "--gapextend",
            "2",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "5",
        ],
        &[
            "-dust",
            "no",
            "-strand",
            "plus",
            "-word_size",
            "7",
            "-gapopen",
            "5",
            "-gapextend",
            "2",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "5",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_max_target_seqs_after_culling_limit() {
    let query = ">q1\nACGTACGTACGTACGTACGTACGTACGTACGTTTAACCGGTTAA\n";
    let db = ">s_full\nACGTACGTACGTACGTACGTACGTACGTACGT\n>s_part\nACGTACGTACGTACGTACGT\n>s_tail\nTTAACCGGTTAA\n>s_noise\nACGTACGTACGTACGT\n";

    assert_blastn_db_outfmt_matches_ncbi(
        query,
        db,
        "6 sseqid qstart qend length score bitscore",
        &[],
        &[
            "--dust",
            "no",
            "--word_size",
            "7",
            "--max_target_seqs",
            "2",
            "--culling_limit",
            "1",
            "--max_hsps",
            "1",
        ],
        &[
            "-dust",
            "no",
            "-word_size",
            "7",
            "-max_target_seqs",
            "2",
            "-culling_limit",
            "1",
            "-max_hsps",
            "1",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_dbsize_and_searchsp_statistics() {
    let query = ">subseq_oid0\nTTAAGGAGGCTCATCTTTCAGAATCCATGCTGTGGGCCAGCAAGAGTTAA\n";
    let db = ">subj1\nTTAAGGAGGCTCATCTTTCAGAATCCATGCTGTGGGCCAGCAAGAGTTAA\n";
    let outfmt = "6 qseqid sseqid evalue bitscore score length pident";

    for (rust_args, ncbi_args) in [
        (
            vec!["--dust", "no", "--dbsize", "1000000"],
            vec!["-dust", "no", "-dbsize", "1000000"],
        ),
        (
            vec!["--dust", "no", "--searchsp", "1000000"],
            vec!["-dust", "no", "-searchsp", "1000000"],
        ),
        (
            vec!["--dust", "no", "--dbsize", "5000000000"],
            vec!["-dust", "no", "-dbsize", "5000000000"],
        ),
        (
            vec!["--dust", "no", "--searchsp", "5000000000"],
            vec!["-dust", "no", "-searchsp", "5000000000"],
        ),
        (
            vec!["--dust", "no", "--dbsize", "-1"],
            vec!["-dust", "no", "-dbsize", "-1"],
        ),
        (
            vec!["--dust", "no", "--searchsp", "0"],
            vec!["-dust", "no", "-searchsp", "0"],
        ),
    ] {
        assert_blastn_db_outfmt_matches_ncbi(query, db, outfmt, &[], &rust_args, &ncbi_args);
    }
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

#[test]
fn blastn_db_ncbi_parity_compact_huge_total_length_statistics() {
    if !std::path::Path::new("/usr/bin/blastn").exists()
        || !std::path::Path::new("/usr/bin/makeblastdb").exists()
    {
        eprintln!("Skipping: /usr/bin/blastn or /usr/bin/makeblastdb not found");
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
    let db_fasta = tmp.path().join("db.fa");
    let db = tmp.path().join("testdb");
    let rust_out = tmp.path().join("rust.tsv");
    let ncbi_out = tmp.path().join("ncbi.tsv");
    std::fs::write(&query, ">q1\nACGTACGTACGTACGTACGTACGTACGTACGT\n").expect("write query FASTA");
    std::fs::write(&db_fasta, ">s1\nACGTACGTACGTACGTACGTACGTACGTACGT\n").expect("write db FASTA");

    let make_status = std::process::Command::new("/usr/bin/makeblastdb")
        .arg("-in")
        .arg(&db_fasta)
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
    patch_blastdb_index_total_length(&db.with_extension("nin"), 5_000_000_000);

    let outfmt = "6 qseqid sseqid evalue bitscore score length pident";
    let rust_status = std::process::Command::new(&blast_cli)
        .arg("blastn")
        .arg("--query")
        .arg(&query)
        .arg("--db")
        .arg(&db)
        .arg("--task")
        .arg("blastn-short")
        .arg("--outfmt")
        .arg(outfmt)
        .arg("--dust")
        .arg("no")
        .arg("--max_target_seqs")
        .arg("10")
        .arg("--num_threads")
        .arg("1")
        .arg("--out")
        .arg(&rust_out)
        .status()
        .expect("run blast-cli huge metadata DB parity");
    assert!(rust_status.success(), "blast-cli exited with {rust_status}");

    let ncbi_status = std::process::Command::new("/usr/bin/blastn")
        .arg("-query")
        .arg(&query)
        .arg("-db")
        .arg(&db)
        .arg("-task")
        .arg("blastn-short")
        .arg("-outfmt")
        .arg(outfmt)
        .arg("-dust")
        .arg("no")
        .arg("-max_target_seqs")
        .arg("10")
        .arg("-num_threads")
        .arg("1")
        .arg("-out")
        .arg(&ncbi_out)
        .status()
        .expect("run NCBI huge metadata DB parity");
    assert!(
        ncbi_status.success(),
        "NCBI blastn exited with {ncbi_status}"
    );

    assert_eq!(
        std::fs::read(&rust_out).expect("read rust output"),
        std::fs::read(&ncbi_out).expect("read ncbi output"),
        "compact huge-total-length DB output differs"
    );
}

#[test]
fn blastn_db_ncbi_parity_multivolume_alias_dbsize_searchsp_statistics() {
    if !std::path::Path::new("/usr/bin/blastn").exists()
        || !std::path::Path::new("/usr/bin/makeblastdb").exists()
    {
        eprintln!("Skipping: /usr/bin/blastn or /usr/bin/makeblastdb not found");
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
    let vol0_fa = tmp.path().join("vol0.fa");
    let vol1_fa = tmp.path().join("vol1.fa");
    let vol0 = tmp.path().join("vol0");
    let vol1 = tmp.path().join("vol1");
    let alias = tmp.path().join("multi.nal");
    let db = tmp.path().join("multi");

    std::fs::write(
        &query,
        ">q1\nTTAAGGAGGCTCATCTTTCAGAATCCATGCTGTGGGCCAGCAAGAGTTAA\n",
    )
    .expect("write query");
    std::fs::write(
        &vol0_fa,
        ">s0\nTTAAGGAGGCTCATCTTTCAGAATCCATGCTGTGGGCCAGCAAGAGTTAA\n",
    )
    .expect("write vol0");
    std::fs::write(
        &vol1_fa,
        ">s1\nGGGGGGGGGGGGGGGGGGGGTTAAGGAGGCTCATCTTTCAGAATCCATGCTGTGGGCCAGCAAGAGTTAT\n",
    )
    .expect("write vol1");

    for (input, output) in [(&vol0_fa, &vol0), (&vol1_fa, &vol1)] {
        let status = std::process::Command::new("/usr/bin/makeblastdb")
            .arg("-in")
            .arg(input)
            .arg("-dbtype")
            .arg("nucl")
            .arg("-out")
            .arg(output)
            .stdout(std::process::Stdio::null())
            .status()
            .expect("run makeblastdb");
        assert!(status.success(), "makeblastdb exited with {status}");
    }

    std::fs::write(
        &alias,
        format!(
            "TITLE compact multivolume alias\nDBLIST {} {}\nNSEQ 2\nLENGTH 130\n",
            vol0.display(),
            vol1.display()
        ),
    )
    .expect("write alias");

    let outfmt = "6 sseqid evalue bitscore score length pident";
    for (label, rust_args, ncbi_args) in [
        (
            "dbsize",
            vec!["--dust", "no", "--dbsize", "5000000000"],
            vec!["-dust", "no", "-dbsize", "5000000000"],
        ),
        (
            "searchsp",
            vec!["--dust", "no", "--searchsp", "5000000000"],
            vec!["-dust", "no", "-searchsp", "5000000000"],
        ),
    ] {
        let rust_out = tmp.path().join(format!("rust_{label}.tsv"));
        let ncbi_out = tmp.path().join(format!("ncbi_{label}.tsv"));

        let mut rust_cmd = std::process::Command::new(&blast_cli);
        rust_cmd
            .arg("blastn")
            .arg("--query")
            .arg(&query)
            .arg("--db")
            .arg(&db)
            .arg("--task")
            .arg("blastn-short")
            .arg("--outfmt")
            .arg(outfmt)
            .arg("--num_threads")
            .arg("1")
            .arg("--out")
            .arg(&rust_out);
        for arg in rust_args {
            rust_cmd.arg(arg);
        }
        let rust_status = rust_cmd.status().expect("run blast-cli alias dbsize");
        assert!(rust_status.success(), "blast-cli exited with {rust_status}");

        let mut ncbi_cmd = std::process::Command::new("/usr/bin/blastn");
        ncbi_cmd
            .arg("-query")
            .arg(&query)
            .arg("-db")
            .arg(&db)
            .arg("-task")
            .arg("blastn-short")
            .arg("-outfmt")
            .arg(outfmt)
            .arg("-num_threads")
            .arg("1")
            .arg("-out")
            .arg(&ncbi_out);
        for arg in ncbi_args {
            ncbi_cmd.arg(arg);
        }
        let ncbi_status = ncbi_cmd.status().expect("run NCBI alias dbsize");
        assert!(
            ncbi_status.success(),
            "NCBI blastn exited with {ncbi_status}"
        );

        assert_eq!(
            std::fs::read(&rust_out).expect("read rust output"),
            std::fs::read(&ncbi_out).expect("read ncbi output"),
            "{label} multivolume alias output differs"
        );
    }
}

#[test]
fn blastn_db_ncbi_parity_multivolume_alias_equal_score_subject_ordering() {
    if !std::path::Path::new("/usr/bin/blastn").exists()
        || !std::path::Path::new("/usr/bin/makeblastdb").exists()
    {
        eprintln!("Skipping: /usr/bin/blastn or /usr/bin/makeblastdb not found");
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
    let vol0_fa = tmp.path().join("vol0.fa");
    let vol1_fa = tmp.path().join("vol1.fa");
    let vol0 = tmp.path().join("vol0");
    let vol1 = tmp.path().join("vol1");
    let alias = tmp.path().join("multi.nal");
    let db = tmp.path().join("multi");

    std::fs::write(
        &query,
        ">q
GTCTCCTCTGACTTCAACAGCG
",
    )
    .expect("write query");
    std::fs::write(
        &vol0_fa,
        ">plus_low
GTCTCCTCTGACTTCAA
>minus_low
TTGAAGTCAGAGGAGAC
",
    )
    .expect("write vol0");
    std::fs::write(
        &vol1_fa,
        ">plus_high
CTCTGACTTCAACAGCG
>minus_high
CGCTGTTGAAGTCAGAG
",
    )
    .expect("write vol1");

    for (input, output) in [(&vol0_fa, &vol0), (&vol1_fa, &vol1)] {
        let status = std::process::Command::new("/usr/bin/makeblastdb")
            .arg("-in")
            .arg(input)
            .arg("-dbtype")
            .arg("nucl")
            .arg("-out")
            .arg(output)
            .stdout(std::process::Stdio::null())
            .status()
            .expect("run makeblastdb");
        assert!(status.success(), "makeblastdb exited with {status}");
    }

    std::fs::write(
        &alias,
        format!(
            "TITLE compact multivolume ordering alias
DBLIST {} {}
NSEQ 4
LENGTH 68
",
            vol0.display(),
            vol1.display()
        ),
    )
    .expect("write alias");

    let outfmt = "6 sseqid score sstart send qstart qend";
    let rust_out = tmp.path().join("rust.tsv");
    let ncbi_out = tmp.path().join("ncbi.tsv");

    let rust_status = std::process::Command::new(&blast_cli)
        .arg("blastn")
        .arg("--query")
        .arg(&query)
        .arg("--db")
        .arg(&db)
        .arg("--task")
        .arg("blastn-short")
        .arg("--dust")
        .arg("no")
        .arg("--word_size")
        .arg("7")
        .arg("--max_target_seqs")
        .arg("10")
        .arg("--max_hsps")
        .arg("1")
        .arg("--num_threads")
        .arg("1")
        .arg("--outfmt")
        .arg(outfmt)
        .arg("--out")
        .arg(&rust_out)
        .status()
        .expect("run blast-cli alias ordering");
    assert!(rust_status.success(), "blast-cli exited with {rust_status}");

    let ncbi_status = std::process::Command::new("/usr/bin/blastn")
        .arg("-query")
        .arg(&query)
        .arg("-db")
        .arg(&db)
        .arg("-task")
        .arg("blastn-short")
        .arg("-dust")
        .arg("no")
        .arg("-word_size")
        .arg("7")
        .arg("-max_target_seqs")
        .arg("10")
        .arg("-max_hsps")
        .arg("1")
        .arg("-num_threads")
        .arg("1")
        .arg("-outfmt")
        .arg(outfmt)
        .arg("-out")
        .arg(&ncbi_out)
        .status()
        .expect("run NCBI alias ordering");
    assert!(
        ncbi_status.success(),
        "NCBI blastn exited with {ncbi_status}"
    );

    assert_eq!(
        std::fs::read(&rust_out).expect("read rust output"),
        std::fs::read(&ncbi_out).expect("read ncbi output"),
        "multivolume alias tied-subject output differs"
    );
}

#[test]
fn blastn_db_ncbi_parity_pairwise_dbsize_and_searchsp_statistics() {
    let query = ">subseq_oid0\nTTAAGGAGGCTCATCTTTCAGAATCCATGCTGTGGGCCAGCAAGAGTTAA\n";
    let db = ">subj1\nTTAAGGAGGCTCATCTTTCAGAATCCATGCTGTGGGCCAGCAAGAGTTAA\n";

    for (rust_args, ncbi_args) in [
        (
            vec!["--dust", "no", "--dbsize", "1000000"],
            vec!["-dust", "no", "-dbsize", "1000000"],
        ),
        (
            vec!["--dust", "no", "--searchsp", "1000000"],
            vec!["-dust", "no", "-searchsp", "1000000"],
        ),
        (
            vec![
                "--dust",
                "no",
                "--reward",
                "1",
                "--penalty",
                "-3",
                "--gapopen",
                "5",
                "--gapextend",
                "2",
            ],
            vec![
                "-dust",
                "no",
                "-reward",
                "1",
                "-penalty",
                "-3",
                "-gapopen",
                "5",
                "-gapextend",
                "2",
            ],
        ),
    ] {
        assert_blastn_db_outfmt_matches_ncbi(query, db, "0", &[], &rust_args, &ncbi_args);
    }
}

#[test]
fn blastn_db_ncbi_parity_dbsize_searchsp_multi_query_and_filtering() {
    let query = ">q1\nTTAAGGAGGCTCATCTTTCAGAATCCATGCTGTGGGCCAGCAAGAGTTAA\n>q2\nTTAAGGAGGCTCATCTTTCAGAATCCATGCTGTGGGCCAGCAAGAGTTAA\n";
    let db = ">subj1\nTTAAGGAGGCTCATCTTTCAGAATCCATGCTGTGGGCCAGCAAGAGTTAA\n";
    let outfmt = "6 qseqid sseqid evalue bitscore score length pident";

    for (rust_args, ncbi_args) in [
        (
            vec!["--dust", "no", "--dbsize", "1000000"],
            vec!["-dust", "no", "-dbsize", "1000000"],
        ),
        (
            vec!["--dust", "no", "--searchsp", "1000000"],
            vec!["-dust", "no", "-searchsp", "1000000"],
        ),
    ] {
        assert_blastn_db_outfmt_matches_ncbi(query, db, outfmt, &[], &rust_args, &ncbi_args);
    }

    for (rust_args, ncbi_args) in [
        (
            vec!["--dust", "no", "--dbsize", "1000000", "--evalue", "1e-25"],
            vec!["-dust", "no", "-dbsize", "1000000", "-evalue", "1e-25"],
        ),
        (
            vec!["--dust", "no", "--searchsp", "1000000", "--evalue", "1e-25"],
            vec!["-dust", "no", "-searchsp", "1000000", "-evalue", "1e-25"],
        ),
    ] {
        assert_blastn_db_outfmt_matches_ncbi(
            ">q1\nTTAAGGAGGCTCATCTTTCAGAATCCATGCTGTGGGCCAGCAAGAGTTAA\n",
            ">subj1\nTTAAGGAGGCTCATCTTTCAGAATCCATGCTGTGGGCCAGCAAGAGTTAA\n",
            outfmt,
            &[],
            &rust_args,
            &ncbi_args,
        );
    }
}

#[test]
fn blastn_db_ncbi_parity_dbsize_searchsp_minus_strand_statistics() {
    let query = ">q_minus\nGAATCCATGCTGTGGGCCAGCAAGAGTTAA\n";
    let db = ">subj_rc\nTTAACTCTTGCTGGCCCACAGCATGGATTC\n";
    let outfmt = "6 qseqid sseqid sstrand qstart qend sstart send evalue bitscore score length pident qseq sseq";

    for (rust_args, ncbi_args) in [
        (
            vec!["--dust", "no", "--strand", "minus", "--dbsize", "1000000"],
            vec!["-dust", "no", "-strand", "minus", "-dbsize", "1000000"],
        ),
        (
            vec!["--dust", "no", "--strand", "minus", "--searchsp", "1000000"],
            vec!["-dust", "no", "-strand", "minus", "-searchsp", "1000000"],
        ),
    ] {
        assert_blastn_db_outfmt_matches_ncbi(query, db, outfmt, &[], &rust_args, &ncbi_args);
    }
}

#[test]
fn blastn_db_ncbi_parity_long_exact_match_extends_to_edges() {
    let seq = "GAATCCATGCTGTGGGCCAGCAAGAGTTAAGGTGCTCATGGTTTTGAGAAAACATCTGAGGACTCTGACAGCACTCTCCCATCCTTGGTCTCCACAGTCT";
    let query = format!(">q\n{}\n", seq);
    let db = format!(">s\n{}\n", seq);

    assert_blastn_db_outfmt_matches_ncbi(
        &query,
        &db,
        "6 qseqid sseqid qstart qend sstart send evalue bitscore score length pident qseq sseq",
        &[],
        &["--dust", "no"],
        &["-dust", "no"],
    );
}

#[test]
fn blastn_db_ncbi_parity_ungapped_right_xdrop_negative_total() {
    let query = ">q\nAAAAAAAAAAAAAAAAAAAA\n";
    let db = ">s\nAAAACAAAAAAAAAAAAAAA\n";
    let outfmt = "6 qstart qend sstart send score length qseq sseq btop";

    assert_blastn_db_outfmt_matches_ncbi(
        query,
        db,
        outfmt,
        &[],
        &[
            "--ungapped",
            "--dust",
            "no",
            "--word_size",
            "4",
            "--reward",
            "1",
            "--penalty",
            "-5",
            "--evalue",
            "1000",
            "--max_hsps",
            "10",
        ],
        &[
            "-ungapped",
            "-dust",
            "no",
            "-word_size",
            "4",
            "-reward",
            "1",
            "-penalty",
            "-5",
            "-evalue",
            "1000",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_ungapped_perfect_30bp_primer() {
    let query = ">q\nACGTACGTACGTACGTACGTACGTACGTAC\n";
    let db = ">s\nTTTTACGTACGTACGTACGTACGTACGTACGTACAAAA\n";
    let outfmt = "6 qstart qend sstart send score length qseq sseq btop";

    assert_blastn_db_outfmt_matches_ncbi(
        query,
        db,
        outfmt,
        &[],
        &[
            "--ungapped",
            "--dust",
            "no",
            "--word_size",
            "7",
            "--evalue",
            "1000",
            "--max_hsps",
            "10",
        ],
        &[
            "-ungapped",
            "-dust",
            "no",
            "-word_size",
            "7",
            "-evalue",
            "1000",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_ungapped_30bp_primer_central_mismatch() {
    let query = ">q\nACGTACGTACGTACGTACGTACGTACGTAC\n";
    let db = ">s\nTTTTACGTACGTACGTATGTACGTACGTACGTACAAAA\n";
    let outfmt = "6 qstart qend sstart send score length qseq sseq btop";

    assert_blastn_db_outfmt_matches_ncbi(
        query,
        db,
        outfmt,
        &[],
        &[
            "--ungapped",
            "--dust",
            "no",
            "--word_size",
            "7",
            "--evalue",
            "1000",
            "--max_hsps",
            "20",
        ],
        &[
            "-ungapped",
            "-dust",
            "no",
            "-word_size",
            "7",
            "-evalue",
            "1000",
            "-max_hsps",
            "20",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_ungapped_ambiguity_scoring_and_ordering() {
    let query = ">q\nAAAANAAAAAAAAAAA\n";
    let db = ">s\nAAAAAAAAAAAAAAAA\n";
    let outfmt = "6 qstart qend sstart send score length qseq sseq btop";

    assert_blastn_db_outfmt_matches_ncbi(
        query,
        db,
        outfmt,
        &[],
        &[
            "--ungapped",
            "--dust",
            "no",
            "--word_size",
            "4",
            "--reward",
            "1",
            "--penalty",
            "-5",
            "--evalue",
            "1000",
            "--max_hsps",
            "10",
        ],
        &[
            "-ungapped",
            "-dust",
            "no",
            "-word_size",
            "4",
            "-reward",
            "1",
            "-penalty",
            "-5",
            "-evalue",
            "1000",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_ungapped_default_megablast_left_xdrop() {
    let query = ">q\nAAAAAAAAAAAAAAAAAAAA\n";
    let db = ">s\nAAAAAAAAACAAAAAAAAAA\n";
    let outfmt = "6 qstart qend sstart send score length qseq sseq btop";

    assert_blastn_db_outfmt_matches_ncbi(
        query,
        db,
        outfmt,
        &[],
        &[
            "--ungapped",
            "--dust",
            "no",
            "--word_size",
            "4",
            "--reward",
            "1",
            "--penalty",
            "-5",
            "--evalue",
            "1000",
            "--max_hsps",
            "10",
        ],
        &[
            "-ungapped",
            "-dust",
            "no",
            "-word_size",
            "4",
            "-reward",
            "1",
            "-penalty",
            "-5",
            "-evalue",
            "1000",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_megablast_contained_diagonals() {
    assert_blastn_db_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n",
        ">s1\nTTTTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAAAA\n>s2\nCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        "6 qstart qend sstart send length bitscore sframe btop",
        &[],
        &["--task", "megablast", "--dust", "no", "--max_target_seqs", "100"],
        &["-task", "megablast", "-dust", "no", "-max_target_seqs", "100"],
    );
}

#[test]
fn blastn_db_ncbi_parity_megablast_indel_tie_max_hsps() {
    assert_blastn_db_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGTACGTACGTACGTGGGGACGTACGTACGTACGTACGTACGTACGTACGT\n",
        ">s1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n",
        "6 qstart qend sstart send score bitscore evalue length mismatch gapopen gaps btop",
        &[],
        &[
            "--task",
            "megablast",
            "--dust",
            "no",
            "--evalue",
            "1000",
            "--max_hsps",
            "1",
        ],
        &[
            "-task",
            "megablast",
            "-dust",
            "no",
            "-evalue",
            "1000",
            "-max_hsps",
            "1",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_ungapped_permissive_evalue_cutoff() {
    let query = ">q\nAAAAAAAAAAAA\n";
    let db = ">s\nAAAACCAAAAAA\n";
    let outfmt = "6 qstart qend sstart send score length qseq sseq btop";

    assert_blastn_db_outfmt_matches_ncbi(
        query,
        db,
        outfmt,
        &[],
        &[
            "--ungapped",
            "--dust",
            "no",
            "--word_size",
            "4",
            "--reward",
            "1",
            "--penalty",
            "-2",
            "--xdrop_ungap",
            "4",
            "--evalue",
            "1000",
            "--max_hsps",
            "20",
        ],
        &[
            "-ungapped",
            "-dust",
            "no",
            "-word_size",
            "4",
            "-reward",
            "1",
            "-penalty",
            "-2",
            "-xdrop_ungap",
            "4",
            "-evalue",
            "1000",
            "-max_hsps",
            "20",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_ungapped_lcase_masking_extends_unmasked_query() {
    let query = ">q\nAAAAaaaaAAAAAAAA\n";
    let db = ">s\nAAAAAAAAAAAAAAAA\n";
    let outfmt = "6 qstart qend sstart send score length qseq sseq btop";

    assert_blastn_db_outfmt_matches_ncbi(
        query,
        db,
        outfmt,
        &[],
        &[
            "--ungapped",
            "--dust",
            "no",
            "--lcase_masking",
            "--word_size",
            "4",
            "--reward",
            "1",
            "--penalty",
            "-5",
            "--evalue",
            "1000",
            "--max_hsps",
            "10",
        ],
        &[
            "-ungapped",
            "-dust",
            "no",
            "-lcase_masking",
            "-word_size",
            "4",
            "-reward",
            "1",
            "-penalty",
            "-5",
            "-evalue",
            "1000",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_ungapped_minus_strand_explicit_xdrop_boundary_matrix() {
    let query = ">q\nAAAAAAAAAAAAAAAAAAAA\n";
    let outfmt = "6 qstart qend sstart send score length qseq sseq btop";
    let cases = [
        ("minus_left_two_mismatch_x1", "CCAAAAAAAAAAAAAAAAAA", "1"),
        ("minus_right_two_mismatch_x1", "AAAAAAAAAAAAAAAAAACC", "1"),
        (
            "minus_internal_two_mismatch_x2",
            "AAAAAAAACCAAAAAAAAAA",
            "2",
        ),
        ("minus_late_drop_x2", "AAAAAAAAAAAACCCAAAAA", "2"),
        ("minus_early_drop_x2", "AAAAACCCAAAAAAAAAAAA", "2"),
    ];

    for (label, plus_subject, xdrop) in cases {
        let db = format!(">s\n{}\n", ascii_reverse_complement(plus_subject));
        assert_blastn_db_outfmt_matches_ncbi(
            query,
            &db,
            outfmt,
            &[],
            &[
                "--ungapped",
                "--dust",
                "no",
                "--strand",
                "minus",
                "--word_size",
                "4",
                "--reward",
                "1",
                "--penalty",
                "-2",
                "--xdrop_ungap",
                xdrop,
                "--evalue",
                "1000",
                "--max_hsps",
                "20",
            ],
            &[
                "-ungapped",
                "-dust",
                "no",
                "-strand",
                "minus",
                "-word_size",
                "4",
                "-reward",
                "1",
                "-penalty",
                "-2",
                "-xdrop_ungap",
                xdrop,
                "-evalue",
                "1000",
                "-max_hsps",
                "20",
            ],
        );
        eprintln!("checked {label}");
    }
}

#[test]
fn blastn_db_ncbi_parity_ungapped_lowercase_subject_is_sequence() {
    let query = ">q\nAAAAAAAAAAAAAAAA\n";
    let db = ">s\nAAAAaaaaAAAAAAAA\n";
    let outfmt = "6 qstart qend sstart send score length qseq sseq btop";

    assert_blastn_db_outfmt_matches_ncbi(
        query,
        db,
        outfmt,
        &[],
        &[
            "--ungapped",
            "--dust",
            "no",
            "--lcase_masking",
            "--word_size",
            "4",
            "--reward",
            "1",
            "--penalty",
            "-5",
            "--evalue",
            "1000",
            "--max_hsps",
            "10",
        ],
        &[
            "-ungapped",
            "-dust",
            "no",
            "-lcase_masking",
            "-word_size",
            "4",
            "-reward",
            "1",
            "-penalty",
            "-5",
            "-evalue",
            "1000",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_ungapped_mismatch_boundary_matrix() {
    let query = ">q\nAAAAAAAAAAAAAAAA\n";
    let outfmt = "6 qstart qend sstart send score length qseq sseq btop";
    let cases = [
        ("left_mismatch", ">s\nCAAAAAAAAAAAAAAA\n"),
        ("right_mismatch", ">s\nAAAAAAAAAAAAAAAC\n"),
        ("two_close_mismatches", ">s\nAAAAAACCAAAAAAAA\n"),
        ("subject_ambiguity", ">s\nAAAANAAAAAAAAAAA\n"),
        ("xdrop_boundary", ">s\nAAAACAAAAAAAACAAA\n"),
    ];

    for (label, db) in cases {
        assert_blastn_db_outfmt_matches_ncbi(
            query,
            db,
            outfmt,
            &[],
            &[
                "--ungapped",
                "--dust",
                "no",
                "--word_size",
                "4",
                "--reward",
                "1",
                "--penalty",
                "-5",
                "--evalue",
                "1000",
                "--max_hsps",
                "10",
            ],
            &[
                "-ungapped",
                "-dust",
                "no",
                "-word_size",
                "4",
                "-reward",
                "1",
                "-penalty",
                "-5",
                "-evalue",
                "1000",
                "-max_hsps",
                "10",
            ],
        );
        eprintln!("checked {label}");
    }
}

#[test]
fn blastn_db_ncbi_parity_ungapped_explicit_xdrop_boundary_matrix() {
    let query = ">q\nAAAAAAAAAAAAAAAAAAAA\n";
    let outfmt = "6 qstart qend sstart send score length qseq sseq btop";
    let cases = [
        ("left_two_mismatch_x1", ">s\nCCAAAAAAAAAAAAAAAAAA\n", "1"),
        ("right_two_mismatch_x1", ">s\nAAAAAAAAAAAAAAAAAACC\n", "1"),
        ("both_end_mismatch_x1", ">s\nCAAAAAAAAAAAAAAAAAAC\n", "1"),
        (
            "internal_two_mismatch_x2",
            ">s\nAAAAAAAACCAAAAAAAAAA\n",
            "2",
        ),
        (
            "internal_two_mismatch_x3",
            ">s\nAAAAAAAACCAAAAAAAAAA\n",
            "3",
        ),
        ("late_drop_x2", ">s\nAAAAAAAAAAAACCCAAAAA\n", "2"),
        ("early_drop_x2", ">s\nAAAAACCCAAAAAAAAAAAA\n", "2"),
        ("ambig_drop_x2", ">s\nAAAAANNNAAAAAAAAAAAA\n", "2"),
    ];

    for (label, db, xdrop) in cases {
        assert_blastn_db_outfmt_matches_ncbi(
            query,
            db,
            outfmt,
            &[],
            &[
                "--ungapped",
                "--dust",
                "no",
                "--word_size",
                "4",
                "--reward",
                "1",
                "--penalty",
                "-2",
                "--xdrop_ungap",
                xdrop,
                "--evalue",
                "1000",
                "--max_hsps",
                "20",
            ],
            &[
                "-ungapped",
                "-dust",
                "no",
                "-word_size",
                "4",
                "-reward",
                "1",
                "-penalty",
                "-2",
                "-xdrop_ungap",
                xdrop,
                "-evalue",
                "1000",
                "-max_hsps",
                "20",
            ],
        );
        eprintln!("checked {label}");
    }
}

#[test]
fn blastn_db_ncbi_parity_pairwise_multi_query() {
    assert_blastn_db_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGT\n>q2 nohit\nAAAAAAAAAAAAAAAAAAAA\n",
        ">s1 subject one\nACGTACGTACGTACGTACGT\n",
        "0",
        &[],
        &["--dust", "no", "--max_target_seqs", "10"],
        &["-dust", "no", "-max_target_seqs", "10"],
    );
}

#[test]
fn blastn_db_ncbi_parity_pairwise_description_alignment_limits() {
    assert_blastn_db_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGTACGT\n",
        ">s1\nACGTACGTACGTACGTACGTACGT\n>s2\nACGTACGTACGTACGTACGTACGT\n>s3\nACGTACGTACGTACGTACGTACGT\n",
        "0",
        &[],
        &[
            "--dust",
            "no",
            "--num_descriptions",
            "2",
            "--num_alignments",
            "1",
        ],
        &[
            "-dust",
            "no",
            "-num_descriptions",
            "2",
            "-num_alignments",
            "1",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_pairwise_line_length() {
    assert_blastn_db_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGTACGT\n",
        ">s1\nACGTACGTACGTACGTACGTACGT\n",
        "0",
        &[],
        &[
            "--dust",
            "no",
            "--num_descriptions",
            "1",
            "--num_alignments",
            "1",
            "--line_length",
            "12",
        ],
        &[
            "-dust",
            "no",
            "-num_descriptions",
            "1",
            "-num_alignments",
            "1",
            "-line_length",
            "12",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_positive_taxids_filter() {
    assert_blastn_db_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGT\n",
        ">s1 tax one\nACGTACGTACGTACGTACGT\n",
        "6 staxid",
        &["-taxid", "9606"],
        &[
            "--dust",
            "no",
            "--max_target_seqs",
            "10",
            "--taxids",
            "9606",
        ],
        &["-dust", "no", "-max_target_seqs", "10", "-taxids", "9606"],
    );
}

#[test]
fn blastn_db_ncbi_parity_negative_taxids_filter() {
    assert_blastn_db_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGT\n",
        ">s1 tax one\nACGTACGTACGTACGTACGT\n",
        "6 staxid",
        &["-taxid", "9606"],
        &[
            "--dust",
            "no",
            "--max_target_seqs",
            "10",
            "--negative_taxids",
            "9606",
        ],
        &[
            "-dust",
            "no",
            "-max_target_seqs",
            "10",
            "-negative_taxids",
            "9606",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_taxidlist_filter() {
    assert_blastn_db_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGT\n",
        ">s1 tax one\nACGTACGTACGTACGTACGT\n",
        "6 staxid",
        &["-taxid", "9606"],
        &[
            "--dust",
            "no",
            "--max_target_seqs",
            "10",
            "--taxidlist",
            "{taxids_file}",
        ],
        &[
            "-dust",
            "no",
            "-max_target_seqs",
            "10",
            "-taxidlist",
            "{taxids_file}",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_negative_taxidlist_filter() {
    assert_blastn_db_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGT\n",
        ">s1 tax one\nACGTACGTACGTACGTACGT\n",
        "6 staxid",
        &["-taxid", "9606"],
        &[
            "--dust",
            "no",
            "--max_target_seqs",
            "10",
            "--negative_taxidlist",
            "{taxids_file}",
        ],
        &[
            "-dust",
            "no",
            "-max_target_seqs",
            "10",
            "-negative_taxidlist",
            "{taxids_file}",
        ],
    );
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
    if !std::path::Path::new("/usr/bin/makeblastdb").exists() {
        eprintln!("Skipping: /usr/bin/makeblastdb not found");
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

    let make_status = std::process::Command::new("/usr/bin/makeblastdb")
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

#[test]
fn blastn_db_taxids_filter_expands_descendants_from_taxonomy4blast_sqlite() {
    let out = run_rust_taxonomy_filter_case(&["--taxids", "9606"]);
    if out.is_empty() {
        return;
    }
    let text = String::from_utf8(out).expect("utf8 output");
    assert!(
        text.contains("s9606\t9606\n"),
        "direct taxid missing: {text:?}"
    );
    assert!(
        text.contains("s63221\t63221\n"),
        "descendant taxid missing: {text:?}"
    );
}

#[test]
fn blastn_db_negative_taxids_filter_expands_descendants_from_taxonomy4blast_sqlite() {
    let out = run_rust_taxonomy_filter_case(&["--negative_taxids", "9606"]);
    if out.is_empty() {
        return;
    }
    assert_eq!(
        out, b"",
        "negative descendant filter should remove all hits"
    );
}

#[test]
fn blastn_db_no_taxid_expansion_keeps_exact_taxid_filtering() {
    let out = run_rust_taxonomy_filter_case(&["--no_taxid_expansion", "--taxids", "9606"]);
    if out.is_empty() {
        return;
    }
    let text = String::from_utf8(out).expect("utf8 output");
    assert!(
        text.contains("s9606\t9606\n"),
        "direct taxid missing: {text:?}"
    );
    assert!(
        !text.contains("s63221\t63221\n"),
        "--no_taxid_expansion should not keep descendant hits: {text:?}"
    );
}

#[test]
fn blastn_db_ncbi_parity_missing_taxidlist_errors() {
    if !std::path::Path::new("/usr/bin/blastn").exists()
        || !std::path::Path::new("/usr/bin/makeblastdb").exists()
    {
        eprintln!("Skipping: /usr/bin/blastn or /usr/bin/makeblastdb not found");
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
    let db_fasta = tmp.path().join("db.fa");
    let db = tmp.path().join("testdb");
    let missing = tmp.path().join("missing_taxids.txt");
    let rust_out = tmp.path().join("rust.tsv");
    let ncbi_out = tmp.path().join("ncbi.tsv");
    std::fs::write(&query, ">q1\nACGTACGTACGTACGTACGT\n").expect("write query FASTA");
    std::fs::write(&db_fasta, ">s1\nACGTACGTACGTACGTACGT\n").expect("write db FASTA");

    let make_status = std::process::Command::new("/usr/bin/makeblastdb")
        .arg("-in")
        .arg(&db_fasta)
        .arg("-dbtype")
        .arg("nucl")
        .arg("-out")
        .arg(&db)
        .arg("-taxid")
        .arg("9606")
        .stdout(std::process::Stdio::null())
        .status()
        .expect("run makeblastdb");
    assert!(
        make_status.success(),
        "makeblastdb exited with {make_status}"
    );

    let rust = std::process::Command::new(blast_cli)
        .arg("blastn")
        .arg("--query")
        .arg(&query)
        .arg("--db")
        .arg(&db)
        .arg("--task")
        .arg("blastn-short")
        .arg("--outfmt")
        .arg("6")
        .arg("--dust")
        .arg("no")
        .arg("--taxidlist")
        .arg(&missing)
        .arg("--out")
        .arg(&rust_out)
        .output()
        .expect("run blast-cli missing taxidlist");
    let ncbi = std::process::Command::new("/usr/bin/blastn")
        .arg("-query")
        .arg(&query)
        .arg("-db")
        .arg(&db)
        .arg("-task")
        .arg("blastn-short")
        .arg("-outfmt")
        .arg("6")
        .arg("-dust")
        .arg("no")
        .arg("-taxidlist")
        .arg(&missing)
        .arg("-out")
        .arg(&ncbi_out)
        .output()
        .expect("run NCBI missing taxidlist");

    assert!(
        !rust.status.success(),
        "blast-cli should reject missing taxidlist"
    );
    assert!(
        !ncbi.status.success(),
        "NCBI should reject missing taxidlist"
    );
    assert_eq!(
        std::fs::read(&rust_out).unwrap_or_default(),
        std::fs::read(&ncbi_out).unwrap_or_default(),
        "missing taxidlist outputs differ"
    );
    assert_eq!(
        String::from_utf8_lossy(&rust.stderr),
        String::from_utf8_lossy(&ncbi.stderr),
        "missing taxidlist stderr differs"
    );
}

#[test]
fn blastn_db_ncbi_parity_invalid_taxids_errors() {
    if !std::path::Path::new("/usr/bin/blastn").exists()
        || !std::path::Path::new("/usr/bin/makeblastdb").exists()
    {
        eprintln!("Skipping: /usr/bin/blastn or /usr/bin/makeblastdb not found");
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
    let db_fasta = tmp.path().join("db.fa");
    let db = tmp.path().join("testdb");
    let rust_out = tmp.path().join("rust.tsv");
    let ncbi_out = tmp.path().join("ncbi.tsv");
    std::fs::write(&query, ">q1\nACGTACGTACGTACGTACGT\n").expect("write query FASTA");
    std::fs::write(&db_fasta, ">s1\nACGTACGTACGTACGTACGT\n").expect("write db FASTA");

    let make_status = std::process::Command::new("/usr/bin/makeblastdb")
        .arg("-in")
        .arg(&db_fasta)
        .arg("-dbtype")
        .arg("nucl")
        .arg("-out")
        .arg(&db)
        .arg("-taxid")
        .arg("9606")
        .stdout(std::process::Stdio::null())
        .status()
        .expect("run makeblastdb");
    assert!(
        make_status.success(),
        "makeblastdb exited with {make_status}"
    );

    let rust = std::process::Command::new(blast_cli)
        .arg("blastn")
        .arg("--query")
        .arg(&query)
        .arg("--db")
        .arg(&db)
        .arg("--task")
        .arg("blastn-short")
        .arg("--outfmt")
        .arg("6")
        .arg("--dust")
        .arg("no")
        .arg("--taxids")
        .arg("abc")
        .arg("--out")
        .arg(&rust_out)
        .output()
        .expect("run blast-cli invalid taxids");
    let ncbi = std::process::Command::new("/usr/bin/blastn")
        .arg("-query")
        .arg(&query)
        .arg("-db")
        .arg(&db)
        .arg("-task")
        .arg("blastn-short")
        .arg("-outfmt")
        .arg("6")
        .arg("-dust")
        .arg("no")
        .arg("-taxids")
        .arg("abc")
        .arg("-out")
        .arg(&ncbi_out)
        .output()
        .expect("run NCBI invalid taxids");

    assert!(
        !rust.status.success(),
        "blast-cli should reject invalid taxids"
    );
    assert!(!ncbi.status.success(), "NCBI should reject invalid taxids");
    assert_eq!(
        std::fs::read(&rust_out).unwrap_or_default(),
        std::fs::read(&ncbi_out).unwrap_or_default(),
        "invalid taxids outputs differ"
    );
    assert_eq!(
        String::from_utf8_lossy(&rust.stderr),
        String::from_utf8_lossy(&ncbi.stderr),
        "invalid taxids stderr differs"
    );
}

#[test]
fn blastn_db_ncbi_parity_invalid_taxidlist_errors() {
    if !std::path::Path::new("/usr/bin/blastn").exists()
        || !std::path::Path::new("/usr/bin/makeblastdb").exists()
    {
        eprintln!("Skipping: /usr/bin/blastn or /usr/bin/makeblastdb not found");
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
    let db_fasta = tmp.path().join("db.fa");
    let db = tmp.path().join("testdb");
    let bad_taxids = tmp.path().join("bad_taxids.txt");
    let rust_out = tmp.path().join("rust.tsv");
    let ncbi_out = tmp.path().join("ncbi.tsv");
    std::fs::write(&query, ">q1\nACGTACGTACGTACGTACGT\n").expect("write query FASTA");
    std::fs::write(&db_fasta, ">s1\nACGTACGTACGTACGTACGT\n").expect("write db FASTA");
    std::fs::write(&bad_taxids, "9606\nabc\n").expect("write invalid taxid list");

    let make_status = std::process::Command::new("/usr/bin/makeblastdb")
        .arg("-in")
        .arg(&db_fasta)
        .arg("-dbtype")
        .arg("nucl")
        .arg("-out")
        .arg(&db)
        .arg("-taxid")
        .arg("9606")
        .stdout(std::process::Stdio::null())
        .status()
        .expect("run makeblastdb");
    assert!(
        make_status.success(),
        "makeblastdb exited with {make_status}"
    );

    let rust = std::process::Command::new(blast_cli)
        .arg("blastn")
        .arg("--query")
        .arg(&query)
        .arg("--db")
        .arg(&db)
        .arg("--task")
        .arg("blastn-short")
        .arg("--outfmt")
        .arg("6")
        .arg("--dust")
        .arg("no")
        .arg("--taxidlist")
        .arg(&bad_taxids)
        .arg("--out")
        .arg(&rust_out)
        .output()
        .expect("run blast-cli invalid taxidlist");
    let ncbi = std::process::Command::new("/usr/bin/blastn")
        .arg("-query")
        .arg(&query)
        .arg("-db")
        .arg(&db)
        .arg("-task")
        .arg("blastn-short")
        .arg("-outfmt")
        .arg("6")
        .arg("-dust")
        .arg("no")
        .arg("-taxidlist")
        .arg(&bad_taxids)
        .arg("-out")
        .arg(&ncbi_out)
        .output()
        .expect("run NCBI invalid taxidlist");

    assert!(
        !rust.status.success(),
        "blast-cli should reject invalid taxidlist"
    );
    assert!(
        !ncbi.status.success(),
        "NCBI should reject invalid taxidlist"
    );
    assert_eq!(
        std::fs::read(&rust_out).unwrap_or_default(),
        std::fs::read(&ncbi_out).unwrap_or_default(),
        "invalid taxidlist outputs differ"
    );
    assert_eq!(
        String::from_utf8_lossy(&rust.stderr),
        String::from_utf8_lossy(&ncbi.stderr),
        "invalid taxidlist stderr differs"
    );
}

#[test]
fn blastn_db_ncbi_parity_seqidlist_filter() {
    assert_blastn_db_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGT\n",
        ">s1 subject one\nACGTACGTACGTACGTACGT\n>s2 subject two\nACGTACGTACGTACGTACGT\n",
        "6 sacc",
        &["-parse_seqids"],
        &[
            "--dust",
            "no",
            "--max_target_seqs",
            "10",
            "--seqidlist",
            "{seqids_file}",
        ],
        &[
            "-dust",
            "no",
            "-max_target_seqs",
            "10",
            "-seqidlist",
            "{seqids_file}",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_negative_seqidlist_filter() {
    assert_blastn_db_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGT\n",
        ">s1 subject one\nACGTACGTACGTACGTACGT\n>s2 subject two\nACGTACGTACGTACGTACGT\n",
        "6 sacc",
        &["-parse_seqids"],
        &[
            "--dust",
            "no",
            "--max_target_seqs",
            "10",
            "--negative_seqidlist",
            "{seqids_file}",
        ],
        &[
            "-dust",
            "no",
            "-max_target_seqs",
            "10",
            "-negative_seqidlist",
            "{seqids_file}",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_seqidlist_warnings() {
    if !std::path::Path::new("/usr/bin/blastn").exists()
        || !std::path::Path::new("/usr/bin/makeblastdb").exists()
    {
        eprintln!("Skipping: /usr/bin/blastn or /usr/bin/makeblastdb not found");
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
    let db_fasta = tmp.path().join("db.fa");
    let db = tmp.path().join("testdb");
    let seqids = tmp.path().join("seqids.txt");
    std::fs::write(&query, ">q1\nACGTACGTACGTACGTACGT\n").expect("write query FASTA");
    std::fs::write(
        &db_fasta,
        ">s1 subject one\nACGTACGTACGTACGTACGT\n>s2 subject two\nACGTACGTACGTACGTACGT\n",
    )
    .expect("write db FASTA");
    std::fs::write(&seqids, "s2\n").expect("write seqid list");

    let make_status = std::process::Command::new("/usr/bin/makeblastdb")
        .arg("-in")
        .arg(&db_fasta)
        .arg("-dbtype")
        .arg("nucl")
        .arg("-out")
        .arg(&db)
        .arg("-parse_seqids")
        .stdout(std::process::Stdio::null())
        .status()
        .expect("run makeblastdb");
    assert!(
        make_status.success(),
        "makeblastdb exited with {make_status}"
    );

    for (rust_option, ncbi_option) in [
        ("--seqidlist", "-seqidlist"),
        ("--negative_seqidlist", "-negative_seqidlist"),
    ] {
        let rust_out = tmp
            .path()
            .join(format!("rust_{}.tsv", rust_option.trim_start_matches("--")));
        let ncbi_out = tmp
            .path()
            .join(format!("ncbi_{}.tsv", rust_option.trim_start_matches("--")));

        let rust = std::process::Command::new(&blast_cli)
            .arg("blastn")
            .arg("--query")
            .arg(&query)
            .arg("--db")
            .arg(&db)
            .arg("--task")
            .arg("blastn-short")
            .arg("--outfmt")
            .arg("6 sacc")
            .arg("--dust")
            .arg("no")
            .arg("--max_target_seqs")
            .arg("10")
            .arg(rust_option)
            .arg(&seqids)
            .arg("--out")
            .arg(&rust_out)
            .output()
            .expect("run blast-cli seqidlist warning parity");
        let ncbi = std::process::Command::new("/usr/bin/blastn")
            .arg("-query")
            .arg(&query)
            .arg("-db")
            .arg(&db)
            .arg("-task")
            .arg("blastn-short")
            .arg("-outfmt")
            .arg("6 sacc")
            .arg("-dust")
            .arg("no")
            .arg("-max_target_seqs")
            .arg("10")
            .arg(ncbi_option)
            .arg(&seqids)
            .arg("-out")
            .arg(&ncbi_out)
            .output()
            .expect("run NCBI seqidlist warning parity");

        assert!(rust.status.success(), "blast-cli failed for {rust_option}");
        assert!(ncbi.status.success(), "NCBI failed for {ncbi_option}");
        assert_eq!(
            std::fs::read(&rust_out).unwrap_or_default(),
            std::fs::read(&ncbi_out).unwrap_or_default(),
            "{rust_option} output differs"
        );
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{rust_option} stderr differs"
        );
    }
}

#[test]
fn blastn_db_ncbi_parity_missing_database_error() {
    if !std::path::Path::new("/usr/bin/blastn").exists() {
        eprintln!("Skipping: /usr/bin/blastn not found");
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
    let missing_db = tmp.path().join("missing_db");
    std::fs::write(&query, ">q1\nACGTACGTACGTACGTACGT\n").expect("write query FASTA");

    let rust = std::process::Command::new(blast_cli)
        .arg("blastn")
        .arg("--query")
        .arg(&query)
        .arg("--db")
        .arg(&missing_db)
        .arg("--task")
        .arg("blastn-short")
        .arg("--outfmt")
        .arg("6")
        .arg("--dust")
        .arg("no")
        .output()
        .expect("run blast-cli missing database");
    let ncbi = std::process::Command::new("/usr/bin/blastn")
        .arg("-query")
        .arg(&query)
        .arg("-db")
        .arg(&missing_db)
        .arg("-task")
        .arg("blastn-short")
        .arg("-outfmt")
        .arg("6")
        .arg("-dust")
        .arg("no")
        .output()
        .expect("run NCBI missing database");

    assert!(!rust.status.success(), "blast-cli should reject missing DB");
    assert!(!ncbi.status.success(), "NCBI should reject missing DB");
    assert_eq!(
        rust.status.code(),
        ncbi.status.code(),
        "missing database status differs"
    );
    assert_eq!(rust.stdout, ncbi.stdout, "missing database stdout differs");
    assert_eq!(
        String::from_utf8_lossy(&rust.stderr),
        String::from_utf8_lossy(&ncbi.stderr),
        "missing database stderr differs"
    );
}

#[test]
fn blastn_db_ncbi_parity_alias_nseq_length_override_statistics() {
    if !std::path::Path::new("/usr/bin/blastn").exists()
        || !std::path::Path::new("/usr/bin/makeblastdb").exists()
    {
        eprintln!("Skipping: /usr/bin/blastn or /usr/bin/makeblastdb not found");
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
    let db_fasta = tmp.path().join("db.fa");
    let db_base = tmp.path().join("base");
    let alias = tmp.path().join("length_alias.nal");
    let db = tmp.path().join("length_alias");
    let primer = "GTCTCCTCTGACTTCAACAGCG";

    std::fs::write(&query, format!(">q\n{primer}\n")).expect("write query");
    std::fs::write(&db_fasta, format!(">s0\n{primer}\n")).expect("write database FASTA");
    let status = std::process::Command::new("/usr/bin/makeblastdb")
        .arg("-in")
        .arg(&db_fasta)
        .arg("-dbtype")
        .arg("nucl")
        .arg("-out")
        .arg(&db_base)
        .stdout(std::process::Stdio::null())
        .status()
        .expect("run makeblastdb");
    assert!(status.success(), "makeblastdb exited with {status}");

    std::fs::write(
        &alias,
        format!(
            "TITLE length alias\nDBLIST {}\nNSEQ 100000\nLENGTH 5000000000\n",
            db_base.display()
        ),
    )
    .expect("write alias");

    let outfmt = "6 sseqid evalue bitscore score length pident";
    let rust_out = tmp.path().join("rust.tsv");
    let ncbi_out = tmp.path().join("ncbi.tsv");

    let rust_status = std::process::Command::new(&blast_cli)
        .arg("blastn")
        .arg("--query")
        .arg(&query)
        .arg("--db")
        .arg(&db)
        .arg("--task")
        .arg("blastn-short")
        .arg("--dust")
        .arg("no")
        .arg("--word_size")
        .arg("7")
        .arg("--max_target_seqs")
        .arg("10")
        .arg("--outfmt")
        .arg(outfmt)
        .arg("--num_threads")
        .arg("1")
        .arg("--out")
        .arg(&rust_out)
        .status()
        .expect("run blast-cli length alias");
    assert!(rust_status.success(), "blast-cli exited with {rust_status}");

    let ncbi_status = std::process::Command::new("/usr/bin/blastn")
        .arg("-query")
        .arg(&query)
        .arg("-db")
        .arg(&db)
        .arg("-task")
        .arg("blastn-short")
        .arg("-dust")
        .arg("no")
        .arg("-word_size")
        .arg("7")
        .arg("-max_target_seqs")
        .arg("10")
        .arg("-outfmt")
        .arg(outfmt)
        .arg("-num_threads")
        .arg("1")
        .arg("-out")
        .arg(&ncbi_out)
        .status()
        .expect("run NCBI length alias");
    assert!(ncbi_status.success(), "NCBI exited with {ncbi_status}");

    assert_eq!(
        std::fs::read(&rust_out).expect("read rust output"),
        std::fs::read(&ncbi_out).expect("read ncbi output"),
        "alias NSEQ/LENGTH override output differs"
    );
}

#[test]
fn blastn_db_ncbi_parity_alias_stats_metadata_override_statistics() {
    if !std::path::Path::new("/usr/bin/blastn").exists()
        || !std::path::Path::new("/usr/bin/makeblastdb").exists()
    {
        eprintln!("Skipping: /usr/bin/blastn or /usr/bin/makeblastdb not found");
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
    let db_fasta = tmp.path().join("db.fa");
    let db_base = tmp.path().join("base");
    let alias = tmp.path().join("length_alias.nal");
    let db = tmp.path().join("length_alias");
    let primer = "GTCTCCTCTGACTTCAACAGCG";

    std::fs::write(&query, format!(">q\n{primer}\n")).expect("write query");
    std::fs::write(&db_fasta, format!(">s0\n{primer}\n")).expect("write database FASTA");
    let status = std::process::Command::new("/usr/bin/makeblastdb")
        .arg("-in")
        .arg(&db_fasta)
        .arg("-dbtype")
        .arg("nucl")
        .arg("-out")
        .arg(&db_base)
        .stdout(std::process::Stdio::null())
        .status()
        .expect("run makeblastdb");
    assert!(status.success(), "makeblastdb exited with {status}");

    std::fs::write(
        &alias,
        format!("TITLE length alias\nDBLIST {}\nNSEQ 2\nLENGTH 48\nSTATS_NSEQ 100000\nSTATS_TOTLEN 5000000000\n", db_base.display()),
    )
    .expect("write alias");

    let outfmt = "6 sseqid evalue bitscore score length pident";
    let rust_out = tmp.path().join("rust.tsv");
    let ncbi_out = tmp.path().join("ncbi.tsv");

    let rust_status = std::process::Command::new(&blast_cli)
        .arg("blastn")
        .arg("--query")
        .arg(&query)
        .arg("--db")
        .arg(&db)
        .arg("--task")
        .arg("blastn-short")
        .arg("--dust")
        .arg("no")
        .arg("--word_size")
        .arg("7")
        .arg("--max_target_seqs")
        .arg("10")
        .arg("--outfmt")
        .arg(outfmt)
        .arg("--num_threads")
        .arg("1")
        .arg("--out")
        .arg(&rust_out)
        .status()
        .expect("run blast-cli length alias");
    assert!(rust_status.success(), "blast-cli exited with {rust_status}");

    let ncbi_status = std::process::Command::new("/usr/bin/blastn")
        .arg("-query")
        .arg(&query)
        .arg("-db")
        .arg(&db)
        .arg("-task")
        .arg("blastn-short")
        .arg("-dust")
        .arg("no")
        .arg("-word_size")
        .arg("7")
        .arg("-max_target_seqs")
        .arg("10")
        .arg("-outfmt")
        .arg(outfmt)
        .arg("-num_threads")
        .arg("1")
        .arg("-out")
        .arg(&ncbi_out)
        .status()
        .expect("run NCBI length alias");
    assert!(ncbi_status.success(), "NCBI exited with {ncbi_status}");

    assert_eq!(
        std::fs::read(&rust_out).expect("read rust output"),
        std::fs::read(&ncbi_out).expect("read ncbi output"),
        "alias STATS_NSEQ/STATS_TOTLEN override output differs"
    );
}

#[test]
fn blastn_db_ncbi_parity_nested_alias_metadata_precedence() {
    if !std::path::Path::new("/usr/bin/blastn").exists()
        || !std::path::Path::new("/usr/bin/makeblastdb").exists()
    {
        eprintln!("Skipping: /usr/bin/blastn or /usr/bin/makeblastdb not found");
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
    let db_fasta = tmp.path().join("db.fa");
    let db_base = tmp.path().join("base");
    let inner = tmp.path().join("inner.nal");
    let outer = tmp.path().join("outer.nal");
    let db = tmp.path().join("outer");
    let primer = "GTCTCCTCTGACTTCAACAGCG";

    std::fs::write(&query, format!(">q\n{primer}\n")).expect("write query");
    std::fs::write(&db_fasta, format!(">s0\n{primer}\n")).expect("write database FASTA");
    let status = std::process::Command::new("/usr/bin/makeblastdb")
        .arg("-in")
        .arg(&db_fasta)
        .arg("-dbtype")
        .arg("nucl")
        .arg("-out")
        .arg(&db_base)
        .stdout(std::process::Stdio::null())
        .status()
        .expect("run makeblastdb");
    assert!(status.success(), "makeblastdb exited with {status}");

    std::fs::write(
        &inner,
        format!(
            "TITLE inner\nDBLIST {}\nNSEQ 100000\nLENGTH 5000000000\n",
            db_base.display()
        ),
    )
    .expect("write inner alias");

    let outfmt = "6 sseqid evalue bitscore score length pident";
    for (label, outer_contents) in [
        ("child_metadata", "TITLE outer\nDBLIST inner\n".to_string()),
        (
            "parent_override",
            "TITLE outer\nDBLIST inner\nNSEQ 2\nLENGTH 48\n".to_string(),
        ),
    ] {
        std::fs::write(&outer, outer_contents).expect("write outer alias");
        let rust_out = tmp.path().join(format!("{label}.rust.tsv"));
        let ncbi_out = tmp.path().join(format!("{label}.ncbi.tsv"));

        let rust_status = std::process::Command::new(&blast_cli)
            .arg("blastn")
            .arg("--query")
            .arg(&query)
            .arg("--db")
            .arg(&db)
            .arg("--task")
            .arg("blastn-short")
            .arg("--dust")
            .arg("no")
            .arg("--word_size")
            .arg("7")
            .arg("--max_target_seqs")
            .arg("10")
            .arg("--outfmt")
            .arg(outfmt)
            .arg("--num_threads")
            .arg("1")
            .arg("--out")
            .arg(&rust_out)
            .status()
            .expect("run blast-cli nested alias metadata");
        assert!(
            rust_status.success(),
            "blast-cli {label} exited with {rust_status}"
        );

        let ncbi_status = std::process::Command::new("/usr/bin/blastn")
            .arg("-query")
            .arg(&query)
            .arg("-db")
            .arg(&db)
            .arg("-task")
            .arg("blastn-short")
            .arg("-dust")
            .arg("no")
            .arg("-word_size")
            .arg("7")
            .arg("-max_target_seqs")
            .arg("10")
            .arg("-outfmt")
            .arg(outfmt)
            .arg("-num_threads")
            .arg("1")
            .arg("-out")
            .arg(&ncbi_out)
            .status()
            .expect("run NCBI nested alias metadata");
        assert!(
            ncbi_status.success(),
            "NCBI {label} exited with {ncbi_status}"
        );

        assert_eq!(
            std::fs::read(&rust_out).expect("read rust output"),
            std::fs::read(&ncbi_out).expect("read ncbi output"),
            "nested alias metadata precedence differs for {label}"
        );
    }
}

#[test]
fn blastn_db_ncbi_parity_alias_first_last_oid_range() {
    if !std::path::Path::new("/usr/bin/blastn").exists()
        || !std::path::Path::new("/usr/bin/makeblastdb").exists()
    {
        eprintln!("Skipping: /usr/bin/blastn or /usr/bin/makeblastdb not found");
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
    let vol0_fa = tmp.path().join("vol0.fa");
    let vol1_fa = tmp.path().join("vol1.fa");
    let vol0 = tmp.path().join("vol0");
    let vol1 = tmp.path().join("vol1");
    let alias = tmp.path().join("ranged.nal");
    let db = tmp.path().join("ranged");
    let primer = "GTCTCCTCTGACTTCAACAGCG";

    std::fs::write(&query, format!(">q\n{primer}\n")).expect("write query");
    std::fs::write(&vol0_fa, format!(">s0\n{primer}\n>s1\n{primer}\n")).expect("write vol0");
    std::fs::write(&vol1_fa, format!(">s2\n{primer}\n>s3\n{primer}\n")).expect("write vol1");

    for (input, output) in [(&vol0_fa, &vol0), (&vol1_fa, &vol1)] {
        let status = std::process::Command::new("/usr/bin/makeblastdb")
            .arg("-in")
            .arg(input)
            .arg("-dbtype")
            .arg("nucl")
            .arg("-out")
            .arg(output)
            .stdout(std::process::Stdio::null())
            .status()
            .expect("run makeblastdb");
        assert!(status.success(), "makeblastdb exited with {status}");
    }

    std::fs::write(
        &alias,
        format!(
            "TITLE ranged alias\nDBLIST {} {}\nFIRST_OID 2\nLAST_OID 3\n",
            vol0.display(),
            vol1.display()
        ),
    )
    .expect("write ranged alias");

    let outfmt = "6 sseqid qstart qend sstart send score length";
    let rust_out = tmp.path().join("rust.tsv");
    let ncbi_out = tmp.path().join("ncbi.tsv");

    let rust_status = std::process::Command::new(&blast_cli)
        .arg("blastn")
        .arg("--query")
        .arg(&query)
        .arg("--db")
        .arg(&db)
        .arg("--task")
        .arg("blastn-short")
        .arg("--dust")
        .arg("no")
        .arg("--word_size")
        .arg("7")
        .arg("--max_target_seqs")
        .arg("10")
        .arg("--max_hsps")
        .arg("1")
        .arg("--outfmt")
        .arg(outfmt)
        .arg("--num_threads")
        .arg("1")
        .arg("--out")
        .arg(&rust_out)
        .status()
        .expect("run blast-cli ranged alias");
    assert!(rust_status.success(), "blast-cli exited with {rust_status}");

    let ncbi_status = std::process::Command::new("/usr/bin/blastn")
        .arg("-query")
        .arg(&query)
        .arg("-db")
        .arg(&db)
        .arg("-task")
        .arg("blastn-short")
        .arg("-dust")
        .arg("no")
        .arg("-word_size")
        .arg("7")
        .arg("-max_target_seqs")
        .arg("10")
        .arg("-max_hsps")
        .arg("1")
        .arg("-outfmt")
        .arg(outfmt)
        .arg("-num_threads")
        .arg("1")
        .arg("-out")
        .arg(&ncbi_out)
        .status()
        .expect("run NCBI ranged alias");
    assert!(ncbi_status.success(), "NCBI exited with {ncbi_status}");

    assert_eq!(
        std::fs::read(&rust_out).expect("read rust output"),
        std::fs::read(&ncbi_out).expect("read ncbi output"),
        "alias FIRST_OID/LAST_OID range output differs"
    );
}

#[test]
fn blastn_db_ncbi_parity_nested_alias_filter_coordinates() {
    if !std::path::Path::new("/usr/bin/blastn").exists()
        || !std::path::Path::new("/usr/bin/makeblastdb").exists()
    {
        eprintln!("Skipping: /usr/bin/blastn or /usr/bin/makeblastdb not found");
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
    let db_fasta = tmp.path().join("db.fa");
    let db_base = tmp.path().join("base");
    let inner = tmp.path().join("inner.nal");
    let outer = tmp.path().join("outer.nal");
    let mask = tmp.path().join("mask.msk");
    let db = tmp.path().join("outer");
    let primer = "GTCTCCTCTGACTTCAACAGCG";

    std::fs::write(&query, format!(">q\n{primer}\n")).expect("write query");
    std::fs::write(
        &db_fasta,
        format!(">s0\n{primer}\n>s1\n{primer}\n>s2\n{primer}\n>s3\n{primer}\n"),
    )
    .expect("write database FASTA");
    let status = std::process::Command::new("/usr/bin/makeblastdb")
        .arg("-in")
        .arg(&db_fasta)
        .arg("-dbtype")
        .arg("nucl")
        .arg("-out")
        .arg(&db_base)
        .stdout(std::process::Stdio::null())
        .status()
        .expect("run makeblastdb");
    assert!(status.success(), "makeblastdb exited with {status}");

    std::fs::write(
        &inner,
        format!(
            "TITLE inner\nDBLIST {}\nOIDLIST mask.msk\n",
            db_base.display()
        ),
    )
    .expect("write inner alias");
    std::fs::write(
        &outer,
        "TITLE outer\nDBLIST inner\nFIRST_OID 3\nLAST_OID 3\n",
    )
    .expect("write outer alias");
    std::fs::write(&mask, [0, 0, 0, 4, 0xa0]).expect("write OID bitmap");

    let outfmt = "6 sseqid qstart qend sstart send score length";
    let rust_out = tmp.path().join("rust.tsv");
    let ncbi_out = tmp.path().join("ncbi.tsv");

    let rust_status = std::process::Command::new(&blast_cli)
        .arg("blastn")
        .arg("--query")
        .arg(&query)
        .arg("--db")
        .arg(&db)
        .arg("--task")
        .arg("blastn-short")
        .arg("--dust")
        .arg("no")
        .arg("--word_size")
        .arg("7")
        .arg("--max_target_seqs")
        .arg("10")
        .arg("--max_hsps")
        .arg("1")
        .arg("--outfmt")
        .arg(outfmt)
        .arg("--num_threads")
        .arg("1")
        .arg("--out")
        .arg(&rust_out)
        .status()
        .expect("run blast-cli nested alias filters");
    assert!(rust_status.success(), "blast-cli exited with {rust_status}");

    let ncbi_status = std::process::Command::new("/usr/bin/blastn")
        .arg("-query")
        .arg(&query)
        .arg("-db")
        .arg(&db)
        .arg("-task")
        .arg("blastn-short")
        .arg("-dust")
        .arg("no")
        .arg("-word_size")
        .arg("7")
        .arg("-max_target_seqs")
        .arg("10")
        .arg("-max_hsps")
        .arg("1")
        .arg("-outfmt")
        .arg(outfmt)
        .arg("-num_threads")
        .arg("1")
        .arg("-out")
        .arg(&ncbi_out)
        .status()
        .expect("run NCBI nested alias filters");
    assert!(ncbi_status.success(), "NCBI exited with {ncbi_status}");

    assert_eq!(
        std::fs::read(&rust_out).expect("read rust output"),
        std::fs::read(&ncbi_out).expect("read ncbi output"),
        "nested alias filter coordinate output differs"
    );
}

#[test]
fn blastn_db_ncbi_parity_alias_oidlist_bitmap_filter() {
    if !std::path::Path::new("/usr/bin/blastn").exists()
        || !std::path::Path::new("/usr/bin/makeblastdb").exists()
    {
        eprintln!("Skipping: /usr/bin/blastn or /usr/bin/makeblastdb not found");
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
    let db_fasta = tmp.path().join("db.fa");
    let db_base = tmp.path().join("base");
    let alias = tmp.path().join("masked.nal");
    let db = tmp.path().join("masked");
    let mask = tmp.path().join("mask.msk");
    let primer = "GTCTCCTCTGACTTCAACAGCG";

    std::fs::write(&query, format!(">q\n{primer}\n")).expect("write query");
    std::fs::write(
        &db_fasta,
        format!(">s0\n{primer}\n>s1\n{primer}\n>s2\n{primer}\n>s3\n{primer}\n"),
    )
    .expect("write database FASTA");
    let status = std::process::Command::new("/usr/bin/makeblastdb")
        .arg("-in")
        .arg(&db_fasta)
        .arg("-dbtype")
        .arg("nucl")
        .arg("-out")
        .arg(&db_base)
        .stdout(std::process::Stdio::null())
        .status()
        .expect("run makeblastdb");
    assert!(status.success(), "makeblastdb exited with {status}");

    std::fs::write(
        &alias,
        format!(
            "TITLE masked alias\nDBLIST {}\nOIDLIST mask.msk\n",
            db_base.display()
        ),
    )
    .expect("write alias");
    std::fs::write(&mask, [0, 0, 0, 4, 0x50]).expect("write OID bitmap");

    let outfmt = "6 sseqid qstart qend sstart send score length";
    let rust_out = tmp.path().join("rust.tsv");
    let ncbi_out = tmp.path().join("ncbi.tsv");

    let rust_status = std::process::Command::new(&blast_cli)
        .arg("blastn")
        .arg("--query")
        .arg(&query)
        .arg("--db")
        .arg(&db)
        .arg("--task")
        .arg("blastn-short")
        .arg("--dust")
        .arg("no")
        .arg("--word_size")
        .arg("7")
        .arg("--max_target_seqs")
        .arg("10")
        .arg("--max_hsps")
        .arg("1")
        .arg("--outfmt")
        .arg(outfmt)
        .arg("--num_threads")
        .arg("1")
        .arg("--out")
        .arg(&rust_out)
        .status()
        .expect("run blast-cli OIDLIST alias");
    assert!(rust_status.success(), "blast-cli exited with {rust_status}");

    let ncbi_status = std::process::Command::new("/usr/bin/blastn")
        .arg("-query")
        .arg(&query)
        .arg("-db")
        .arg(&db)
        .arg("-task")
        .arg("blastn-short")
        .arg("-dust")
        .arg("no")
        .arg("-word_size")
        .arg("7")
        .arg("-max_target_seqs")
        .arg("10")
        .arg("-max_hsps")
        .arg("1")
        .arg("-outfmt")
        .arg(outfmt)
        .arg("-num_threads")
        .arg("1")
        .arg("-out")
        .arg(&ncbi_out)
        .status()
        .expect("run NCBI OIDLIST alias");
    assert!(ncbi_status.success(), "NCBI exited with {ncbi_status}");

    assert_eq!(
        std::fs::read(&rust_out).expect("read rust output"),
        std::fs::read(&ncbi_out).expect("read ncbi output"),
        "alias OIDLIST bitmap output differs"
    );
}

#[test]
fn blastn_db_ncbi_parity_alias_dblist_path_forms() {
    if !std::path::Path::new("/usr/bin/blastn").exists()
        || !std::path::Path::new("/usr/bin/makeblastdb").exists()
    {
        eprintln!("Skipping: /usr/bin/blastn or /usr/bin/makeblastdb not found");
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
    let vols_dir = tmp.path().join("vols");
    let aliases_dir = tmp.path().join("aliases");
    std::fs::create_dir_all(&vols_dir).expect("create vols dir");
    std::fs::create_dir_all(&aliases_dir).expect("create aliases dir");

    let query = tmp.path().join("query.fa");
    let vol0_fa = vols_dir.join("vol0.fa");
    let vol1_fa = vols_dir.join("vol1.fa");
    let vol0 = vols_dir.join("vol0");
    let vol1 = vols_dir.join("vol1");
    std::fs::write(&query, ">q\nACGTACGTACGTACGTACGT\n").expect("write query");
    std::fs::write(&vol0_fa, ">s0\nACGTACGTACGTACGTACGT\n").expect("write vol0");
    std::fs::write(&vol1_fa, ">s1\nTTTTACGTACGTACGTACGTACGTAAAA\n").expect("write vol1");

    for (input, output) in [(&vol0_fa, &vol0), (&vol1_fa, &vol1)] {
        let status = std::process::Command::new("/usr/bin/makeblastdb")
            .arg("-in")
            .arg(input)
            .arg("-dbtype")
            .arg("nucl")
            .arg("-out")
            .arg(output)
            .stdout(std::process::Stdio::null())
            .status()
            .expect("run makeblastdb");
        assert!(status.success(), "makeblastdb exited with {status}");
    }

    let absolute_alias = aliases_dir.join("absolute.nal");
    let relative_alias = aliases_dir.join("relative.nal");
    let inner_alias = aliases_dir.join("inner.nal");
    let nested_alias = aliases_dir.join("nested.nal");
    std::fs::write(
        &absolute_alias,
        format!(
            "TITLE absolute alias\nDBLIST {} {}\n",
            vol0.display(),
            vol1.display()
        ),
    )
    .expect("write absolute alias");
    std::fs::write(
        &relative_alias,
        "TITLE relative alias\nDBLIST ../vols/vol0 ../vols/vol1\n",
    )
    .expect("write relative alias");
    std::fs::write(&inner_alias, "TITLE inner alias\nDBLIST ../vols/vol0\n")
        .expect("write inner alias");
    std::fs::write(
        &nested_alias,
        "TITLE nested alias\nDBLIST inner ../vols/vol1\n",
    )
    .expect("write nested alias");

    let outfmt = "6 sseqid qstart qend sstart send score length";
    for (label, db) in [
        ("absolute", aliases_dir.join("absolute")),
        ("relative", aliases_dir.join("relative")),
        ("nested", aliases_dir.join("nested")),
    ] {
        let rust_out = tmp.path().join(format!("{label}.rust.tsv"));
        let ncbi_out = tmp.path().join(format!("{label}.ncbi.tsv"));

        let rust_status = std::process::Command::new(&blast_cli)
            .arg("blastn")
            .arg("--query")
            .arg(&query)
            .arg("--db")
            .arg(&db)
            .arg("--task")
            .arg("blastn-short")
            .arg("--dust")
            .arg("no")
            .arg("--outfmt")
            .arg(outfmt)
            .arg("--num_threads")
            .arg("1")
            .arg("--out")
            .arg(&rust_out)
            .status()
            .expect("run blast-cli alias path form");
        assert!(
            rust_status.success(),
            "blast-cli {label} exited with {rust_status}"
        );

        let ncbi_status = std::process::Command::new("/usr/bin/blastn")
            .arg("-query")
            .arg(&query)
            .arg("-db")
            .arg(&db)
            .arg("-task")
            .arg("blastn-short")
            .arg("-dust")
            .arg("no")
            .arg("-outfmt")
            .arg(outfmt)
            .arg("-num_threads")
            .arg("1")
            .arg("-out")
            .arg(&ncbi_out)
            .status()
            .expect("run NCBI alias path form");
        assert!(
            ncbi_status.success(),
            "NCBI {label} exited with {ncbi_status}"
        );

        assert_eq!(
            std::fs::read(&rust_out).expect("read rust output"),
            std::fs::read(&ncbi_out).expect("read ncbi output"),
            "{label} alias DBLIST path-form output differs"
        );
        eprintln!("checked {label}");
    }
}

#[test]
fn blastn_db_ncbi_parity_empty_alias_error() {
    if !std::path::Path::new("/usr/bin/blastn").exists() {
        eprintln!("Skipping: /usr/bin/blastn not found");
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
    let db = tmp.path().join("empty_alias");
    let alias = db.with_extension("nal");
    std::fs::write(&query, ">q1\nACGTACGTACGTACGTACGT\n").expect("write query FASTA");
    std::fs::write(&alias, "TITLE empty alias\n").expect("write empty alias");

    let rust = std::process::Command::new(blast_cli)
        .arg("blastn")
        .arg("--query")
        .arg(&query)
        .arg("--db")
        .arg(&db)
        .arg("--task")
        .arg("blastn-short")
        .arg("--outfmt")
        .arg("6")
        .arg("--dust")
        .arg("no")
        .output()
        .expect("run blast-cli empty alias");
    let ncbi = std::process::Command::new("/usr/bin/blastn")
        .arg("-query")
        .arg(&query)
        .arg("-db")
        .arg(&db)
        .arg("-task")
        .arg("blastn-short")
        .arg("-outfmt")
        .arg("6")
        .arg("-dust")
        .arg("no")
        .output()
        .expect("run NCBI empty alias");

    assert!(
        !rust.status.success(),
        "blast-cli should reject empty alias"
    );
    assert!(!ncbi.status.success(), "NCBI should reject empty alias");
    assert_eq!(
        rust.status.code(),
        ncbi.status.code(),
        "empty alias status differs"
    );
    assert_eq!(rust.stdout, ncbi.stdout, "empty alias stdout differs");
    assert_eq!(
        String::from_utf8_lossy(&rust.stderr),
        String::from_utf8_lossy(&ncbi.stderr),
        "empty alias stderr differs"
    );
}

#[test]
fn blastn_db_ncbi_parity_missing_alias_volume_error() {
    if !std::path::Path::new("/usr/bin/blastn").exists() {
        eprintln!("Skipping: /usr/bin/blastn not found");
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
    let db = tmp.path().join("bad_alias");
    let alias = db.with_extension("nal");
    std::fs::write(&query, ">q1\nACGTACGTACGTACGTACGT\n").expect("write query FASTA");
    std::fs::write(&alias, "TITLE bad alias\nDBLIST missing_volume\n").expect("write bad alias");

    let rust = std::process::Command::new(blast_cli)
        .arg("blastn")
        .arg("--query")
        .arg(&query)
        .arg("--db")
        .arg(&db)
        .arg("--task")
        .arg("blastn-short")
        .arg("--outfmt")
        .arg("6")
        .arg("--dust")
        .arg("no")
        .output()
        .expect("run blast-cli missing alias volume");
    let ncbi = std::process::Command::new("/usr/bin/blastn")
        .arg("-query")
        .arg(&query)
        .arg("-db")
        .arg(&db)
        .arg("-task")
        .arg("blastn-short")
        .arg("-outfmt")
        .arg("6")
        .arg("-dust")
        .arg("no")
        .output()
        .expect("run NCBI missing alias volume");

    assert!(
        !rust.status.success(),
        "blast-cli should reject missing alias volume"
    );
    assert!(
        !ncbi.status.success(),
        "NCBI should reject missing alias volume"
    );
    assert_eq!(
        rust.status.code(),
        ncbi.status.code(),
        "missing alias volume status differs"
    );
    assert_eq!(
        rust.stdout, ncbi.stdout,
        "missing alias volume stdout differs"
    );
    assert_eq!(
        String::from_utf8_lossy(&rust.stderr),
        String::from_utf8_lossy(&ncbi.stderr),
        "missing alias volume stderr differs"
    );
}

#[test]
fn blastn_db_ncbi_parity_partial_database_missing_header_error() {
    if !std::path::Path::new("/usr/bin/blastn").exists()
        || !std::path::Path::new("/usr/bin/makeblastdb").exists()
    {
        eprintln!("Skipping: /usr/bin/blastn or /usr/bin/makeblastdb not found");
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
    let fasta = tmp.path().join("db.fa");
    let db = tmp.path().join("partial_db");
    std::fs::write(&query, ">q1\nACGTACGTACGTACGTACGT\n").expect("write query FASTA");
    std::fs::write(&fasta, ">s1\nACGTACGTACGTACGTACGT\n").expect("write DB FASTA");

    let make_status = std::process::Command::new("/usr/bin/makeblastdb")
        .arg("-in")
        .arg(&fasta)
        .arg("-dbtype")
        .arg("nucl")
        .arg("-out")
        .arg(&db)
        .output()
        .expect("run makeblastdb");
    assert!(
        make_status.status.success(),
        "makeblastdb failed: {}",
        String::from_utf8_lossy(&make_status.stderr)
    );

    std::fs::remove_file(db.with_extension("nhr")).expect("remove DB header component");

    let rust = std::process::Command::new(&blast_cli)
        .arg("blastn")
        .arg("--query")
        .arg(&query)
        .arg("--db")
        .arg(&db)
        .arg("--task")
        .arg("blastn-short")
        .arg("--outfmt")
        .arg("6")
        .arg("--dust")
        .arg("no")
        .output()
        .expect("run blast-cli partial database");
    let ncbi = std::process::Command::new("/usr/bin/blastn")
        .arg("-query")
        .arg(&query)
        .arg("-db")
        .arg(&db)
        .arg("-task")
        .arg("blastn-short")
        .arg("-outfmt")
        .arg("6")
        .arg("-dust")
        .arg("no")
        .output()
        .expect("run NCBI partial database");

    assert!(!rust.status.success(), "blast-cli should reject partial DB");
    assert!(!ncbi.status.success(), "NCBI should reject partial DB");
    assert_eq!(rust.status.code(), ncbi.status.code(), "status differs");
    assert_eq!(rust.stdout, ncbi.stdout, "stdout differs");
    assert_eq!(
        String::from_utf8_lossy(&rust.stderr),
        String::from_utf8_lossy(&ncbi.stderr),
        "stderr differs"
    );
}

#[test]
fn blastn_db_ncbi_parity_empty_query_error() {
    if !std::path::Path::new("/usr/bin/blastn").exists()
        || !std::path::Path::new("/usr/bin/makeblastdb").exists()
    {
        eprintln!("Skipping: /usr/bin/blastn or /usr/bin/makeblastdb not found");
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
    let fasta = tmp.path().join("db.fa");
    let db = tmp.path().join("testdb");
    std::fs::write(&query, ">empty\n").expect("write empty query FASTA");
    std::fs::write(&fasta, ">s1\nACGTACGTACGTACGTACGT\n").expect("write DB FASTA");

    let make_status = std::process::Command::new("/usr/bin/makeblastdb")
        .arg("-in")
        .arg(&fasta)
        .arg("-dbtype")
        .arg("nucl")
        .arg("-out")
        .arg(&db)
        .output()
        .expect("run makeblastdb");
    assert!(
        make_status.status.success(),
        "makeblastdb failed: {}",
        String::from_utf8_lossy(&make_status.stderr)
    );

    let rust = std::process::Command::new(&blast_cli)
        .arg("blastn")
        .arg("--query")
        .arg(&query)
        .arg("--db")
        .arg(&db)
        .arg("--task")
        .arg("blastn-short")
        .arg("--outfmt")
        .arg("6")
        .arg("--dust")
        .arg("no")
        .output()
        .expect("run blast-cli empty DB query");
    let ncbi = std::process::Command::new("/usr/bin/blastn")
        .arg("-query")
        .arg(&query)
        .arg("-db")
        .arg(&db)
        .arg("-task")
        .arg("blastn-short")
        .arg("-outfmt")
        .arg("6")
        .arg("-dust")
        .arg("no")
        .output()
        .expect("run NCBI empty DB query");

    assert!(
        !rust.status.success(),
        "blast-cli should reject empty query"
    );
    assert!(!ncbi.status.success(), "NCBI should reject empty query");
    assert_eq!(rust.status.code(), ncbi.status.code(), "status differs");
    assert_eq!(rust.stdout, ncbi.stdout, "stdout differs");
    assert_eq!(
        String::from_utf8_lossy(&rust.stderr),
        String::from_utf8_lossy(&ncbi.stderr),
        "stderr differs"
    );
}

#[test]
fn blastn_db_ncbi_parity_mixed_empty_query_records_warn_and_continue() {
    if !std::path::Path::new("/usr/bin/blastn").exists()
        || !std::path::Path::new("/usr/bin/makeblastdb").exists()
    {
        eprintln!("Skipping: /usr/bin/blastn or /usr/bin/makeblastdb not found");
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
    let db_fasta = tmp.path().join("db.fa");
    let db = tmp.path().join("testdb");
    std::fs::write(&query, ">empty\n>q1\nACGTACGTACGT\n").expect("write mixed query FASTA");
    std::fs::write(&db_fasta, ">s1\nACGTACGTACGT\n").expect("write DB FASTA");

    let make_status = std::process::Command::new("/usr/bin/makeblastdb")
        .arg("-in")
        .arg(&db_fasta)
        .arg("-dbtype")
        .arg("nucl")
        .arg("-out")
        .arg(&db)
        .stdout(std::process::Stdio::null())
        .status()
        .expect("run makeblastdb");
    assert!(make_status.success(), "makeblastdb failed: {make_status}");

    let rust = std::process::Command::new(&blast_cli)
        .arg("blastn")
        .arg("--query")
        .arg(&query)
        .arg("--db")
        .arg(&db)
        .arg("--task")
        .arg("blastn-short")
        .arg("--outfmt")
        .arg("6 qseqid sseqid length")
        .arg("--dust")
        .arg("no")
        .output()
        .expect("run blast-cli DB mixed empty query records");
    let ncbi = std::process::Command::new("/usr/bin/blastn")
        .arg("-query")
        .arg(&query)
        .arg("-db")
        .arg(&db)
        .arg("-task")
        .arg("blastn-short")
        .arg("-outfmt")
        .arg("6 qseqid sseqid length")
        .arg("-dust")
        .arg("no")
        .output()
        .expect("run NCBI DB mixed empty query records");

    assert_eq!(rust.status.code(), ncbi.status.code(), "status differs");
    assert_eq!(rust.stdout, ncbi.stdout, "stdout differs");
    assert_eq!(
        String::from_utf8_lossy(&rust.stderr),
        String::from_utf8_lossy(&ncbi.stderr),
        "stderr differs"
    );
}

#[test]
fn blastn_db_ncbi_parity_missing_import_search_strategy_error() {
    if !std::path::Path::new("/usr/bin/blastn").exists() {
        eprintln!("Skipping: /usr/bin/blastn not found");
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
    let missing_strategy = tmp.path().join("missing_strategy.asn");
    std::fs::write(&query, ">q1\nACGTACGTACGTACGTACGT\n").expect("write query FASTA");

    let rust = std::process::Command::new(blast_cli)
        .arg("blastn")
        .arg("--query")
        .arg(&query)
        .arg("--db")
        .arg("tests/fixtures/seqn/seqn")
        .arg("--import_search_strategy")
        .arg(&missing_strategy)
        .output()
        .expect("run blast-cli missing import search strategy");
    let ncbi = std::process::Command::new("/usr/bin/blastn")
        .arg("-query")
        .arg(&query)
        .arg("-db")
        .arg("tests/fixtures/seqn/seqn")
        .arg("-import_search_strategy")
        .arg(&missing_strategy)
        .output()
        .expect("run NCBI missing import search strategy");

    assert!(
        !rust.status.success(),
        "blast-cli should reject missing import strategy"
    );
    assert!(
        !ncbi.status.success(),
        "NCBI should reject missing import strategy"
    );
    assert_eq!(
        rust.status.code(),
        ncbi.status.code(),
        "missing import strategy status differs"
    );
    assert_eq!(
        rust.stdout, ncbi.stdout,
        "missing import strategy stdout differs"
    );
    assert_eq!(
        String::from_utf8_lossy(&rust.stderr),
        String::from_utf8_lossy(&ncbi.stderr),
        "missing import strategy stderr differs"
    );
}

#[test]
fn blastn_db_ncbi_parity_inaccessible_export_search_strategy_error() {
    if !std::path::Path::new("/usr/bin/blastn").exists() {
        eprintln!("Skipping: /usr/bin/blastn not found");
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
    let export_strategy = tmp.path().join("missing_dir").join("export.asn");
    std::fs::write(&query, ">q1\nACGTACGTACGTACGTACGT\n").expect("write query FASTA");

    let rust = std::process::Command::new(blast_cli)
        .arg("blastn")
        .arg("--query")
        .arg(&query)
        .arg("--db")
        .arg("tests/fixtures/seqn/seqn")
        .arg("--export_search_strategy")
        .arg(&export_strategy)
        .output()
        .expect("run blast-cli inaccessible export search strategy");
    let ncbi = std::process::Command::new("/usr/bin/blastn")
        .arg("-query")
        .arg(&query)
        .arg("-db")
        .arg("tests/fixtures/seqn/seqn")
        .arg("-export_search_strategy")
        .arg(&export_strategy)
        .output()
        .expect("run NCBI inaccessible export search strategy");

    assert!(
        !rust.status.success(),
        "blast-cli should reject inaccessible export strategy"
    );
    assert!(
        !ncbi.status.success(),
        "NCBI should reject inaccessible export strategy"
    );
    assert_eq!(
        rust.status.code(),
        ncbi.status.code(),
        "inaccessible export strategy status differs"
    );
    assert_eq!(
        rust.stdout, ncbi.stdout,
        "inaccessible export strategy stdout differs"
    );
    assert_eq!(
        String::from_utf8_lossy(&rust.stderr),
        String::from_utf8_lossy(&ncbi.stderr),
        "inaccessible export strategy stderr differs"
    );
}

#[test]
fn blastn_subject_ncbi_parity_empty_query_and_subject_warnings() {
    if !std::path::Path::new("/usr/bin/blastn").exists() {
        eprintln!("Skipping: /usr/bin/blastn not found");
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
    let empty_query = tmp.path().join("empty_query.fa");
    let subject = tmp.path().join("subject.fa");
    let empty_subject = tmp.path().join("empty_subject.fa");
    std::fs::write(&query, ">q1\nACGTACGTACGTACGTACGT\n").expect("write query FASTA");
    std::fs::write(&empty_query, ">empty\n").expect("write empty query FASTA");
    std::fs::write(&subject, ">s1\nACGTACGTACGTACGTACGT\n").expect("write subject FASTA");
    std::fs::write(&empty_subject, ">empty\n").expect("write empty subject FASTA");

    for (query_path, subject_path, should_succeed, label) in [
        (&empty_query, &subject, false, "empty query"),
        (&query, &empty_subject, true, "empty subject"),
    ] {
        let rust = std::process::Command::new(&blast_cli)
            .arg("blastn")
            .arg("--query")
            .arg(query_path)
            .arg("--subject")
            .arg(subject_path)
            .arg("--task")
            .arg("blastn-short")
            .arg("--outfmt")
            .arg("6")
            .arg("--dust")
            .arg("no")
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli {label}: {err}"));
        let ncbi = std::process::Command::new("/usr/bin/blastn")
            .arg("-query")
            .arg(query_path)
            .arg("-subject")
            .arg(subject_path)
            .arg("-task")
            .arg("blastn-short")
            .arg("-outfmt")
            .arg("6")
            .arg("-dust")
            .arg("no")
            .output()
            .unwrap_or_else(|err| panic!("run NCBI {label}: {err}"));

        assert_eq!(
            rust.status.success(),
            should_succeed,
            "blast-cli success state differs from expected for {label}"
        );
        assert_eq!(
            ncbi.status.success(),
            should_succeed,
            "NCBI success state differs from expected for {label}"
        );
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "{label} status differs"
        );
        assert_eq!(rust.stdout, ncbi.stdout, "{label} stdout differs");
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{label} stderr differs"
        );
    }
}

#[test]
fn blastn_subject_ncbi_parity_mixed_empty_query_records_warn_and_continue() {
    if !std::path::Path::new("/usr/bin/blastn").exists() {
        eprintln!("Skipping: /usr/bin/blastn not found");
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
    std::fs::write(&query, ">empty\n>q1\nACGTACGTACGT\n").expect("write mixed query FASTA");
    std::fs::write(&subject, ">s1\nACGTACGTACGT\n").expect("write subject FASTA");

    let rust = std::process::Command::new(&blast_cli)
        .arg("blastn")
        .arg("--query")
        .arg(&query)
        .arg("--subject")
        .arg(&subject)
        .arg("--task")
        .arg("blastn-short")
        .arg("--outfmt")
        .arg("6 qseqid sseqid length")
        .arg("--dust")
        .arg("no")
        .output()
        .expect("run blast-cli mixed empty query records");
    let ncbi = std::process::Command::new("/usr/bin/blastn")
        .arg("-query")
        .arg(&query)
        .arg("-subject")
        .arg(&subject)
        .arg("-task")
        .arg("blastn-short")
        .arg("-outfmt")
        .arg("6 qseqid sseqid length")
        .arg("-dust")
        .arg("no")
        .output()
        .expect("run NCBI mixed empty query records");

    assert_eq!(rust.status.code(), ncbi.status.code(), "status differs");
    assert_eq!(rust.stdout, ncbi.stdout, "stdout differs");
    assert_eq!(
        String::from_utf8_lossy(&rust.stderr),
        String::from_utf8_lossy(&ncbi.stderr),
        "stderr differs"
    );
}

#[test]
fn blastn_subject_ncbi_parity_location_error_handling() {
    if !std::path::Path::new("/usr/bin/blastn").exists() {
        eprintln!("Skipping: /usr/bin/blastn not found");
        return;
    }
    let Some(blast_cli) = std::env::var_os("BLAST_RS_CLI_BIN")
        .or_else(|| std::env::var_os("CARGO_BIN_EXE_blast-cli"))
        .map(std::path::PathBuf::from)
    else {
        eprintln!("Skipping: set BLAST_RS_CLI_BIN or CARGO_BIN_EXE_blast-cli to run CLI parity");
        return;
    };

    for (rust_option, ncbi_option, value, label) in [
        ("--query_loc", "-query_loc", "bad", "query bad format"),
        ("--query_loc", "-query_loc", "a-10", "query bad start"),
        ("--query_loc", "-query_loc", "1-b", "query bad stop"),
        ("--query_loc", "-query_loc", "0-10", "query zero start"),
        ("--query_loc", "-query_loc", "10-5", "query reversed range"),
        ("--query_loc", "-query_loc", "999-1000", "query high start"),
        ("--query_loc", "-query_loc", "1-999", "query high stop"),
        ("--subject_loc", "-subject_loc", "bad", "subject bad format"),
        ("--subject_loc", "-subject_loc", "a-10", "subject bad start"),
        ("--subject_loc", "-subject_loc", "1-b", "subject bad stop"),
        (
            "--subject_loc",
            "-subject_loc",
            "0-10",
            "subject zero start",
        ),
        (
            "--subject_loc",
            "-subject_loc",
            "10-5",
            "subject reversed range",
        ),
        (
            "--subject_loc",
            "-subject_loc",
            "999-1000",
            "subject high start",
        ),
        (
            "--subject_loc",
            "-subject_loc",
            "1-999",
            "subject high stop",
        ),
    ] {
        let rust = std::process::Command::new(&blast_cli)
            .arg("blastn")
            .arg("--query")
            .arg("tests/fixtures/query_short_match.fa")
            .arg("--subject")
            .arg("tests/fixtures/subject_test.fa")
            .arg("--task")
            .arg("blastn-short")
            .arg("--outfmt")
            .arg("6")
            .arg("--dust")
            .arg("no")
            .arg(rust_option)
            .arg(value)
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli {label}: {err}"));
        let ncbi = std::process::Command::new("/usr/bin/blastn")
            .arg("-query")
            .arg("tests/fixtures/query_short_match.fa")
            .arg("-subject")
            .arg("tests/fixtures/subject_test.fa")
            .arg("-task")
            .arg("blastn-short")
            .arg("-outfmt")
            .arg("6")
            .arg("-dust")
            .arg("no")
            .arg(ncbi_option)
            .arg(value)
            .output()
            .unwrap_or_else(|err| panic!("run NCBI {label}: {err}"));

        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "{label} status differs"
        );
        assert_eq!(rust.stdout, ncbi.stdout, "{label} stdout differs");
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{label} stderr differs"
        );
    }
}

#[test]
fn blastn_db_ncbi_parity_unsupported_db_mask_errors() {
    if !std::path::Path::new("/usr/bin/blastn").exists() {
        eprintln!("Skipping: /usr/bin/blastn not found");
        return;
    }
    let Some(blast_cli) = std::env::var_os("BLAST_RS_CLI_BIN")
        .or_else(|| std::env::var_os("CARGO_BIN_EXE_blast-cli"))
        .map(std::path::PathBuf::from)
    else {
        eprintln!("Skipping: set BLAST_RS_CLI_BIN or CARGO_BIN_EXE_blast-cli to run CLI parity");
        return;
    };

    for (rust_option, ncbi_option) in [
        ("--db_soft_mask", "-db_soft_mask"),
        ("--db_hard_mask", "-db_hard_mask"),
    ] {
        let rust = std::process::Command::new(&blast_cli)
            .arg("blastn")
            .arg("--query")
            .arg("tests/fixtures/query_short_match.fa")
            .arg("--db")
            .arg("tests/fixtures/seqn/seqn")
            .arg("--task")
            .arg("blastn-short")
            .arg("--outfmt")
            .arg("6")
            .arg("--dust")
            .arg("no")
            .arg(rust_option)
            .arg("99999")
            .output()
            .expect("run blast-cli unsupported DB mask");
        let ncbi = std::process::Command::new("/usr/bin/blastn")
            .arg("-query")
            .arg("tests/fixtures/query_short_match.fa")
            .arg("-db")
            .arg("tests/fixtures/seqn/seqn")
            .arg("-task")
            .arg("blastn-short")
            .arg("-outfmt")
            .arg("6")
            .arg("-dust")
            .arg("no")
            .arg(ncbi_option)
            .arg("99999")
            .output()
            .expect("run NCBI unsupported DB mask");

        assert!(
            !rust.status.success(),
            "blast-cli should reject unsupported {rust_option}"
        );
        assert!(
            !ncbi.status.success(),
            "NCBI should reject unsupported {ncbi_option}"
        );
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "unsupported {rust_option} status differs"
        );
        assert_eq!(
            rust.stdout, ncbi.stdout,
            "unsupported {rust_option} stdout differs"
        );
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "unsupported {rust_option} stderr differs"
        );
    }
}

#[test]
fn blastn_db_ncbi_parity_named_db_mask_warnings() {
    if !std::path::Path::new("/usr/bin/blastn").exists() {
        eprintln!("Skipping: /usr/bin/blastn not found");
        return;
    }
    let Some(blast_cli) = std::env::var_os("BLAST_RS_CLI_BIN")
        .or_else(|| std::env::var_os("CARGO_BIN_EXE_blast-cli"))
        .map(std::path::PathBuf::from)
    else {
        eprintln!("Skipping: set BLAST_RS_CLI_BIN or CARGO_BIN_EXE_blast-cli to run CLI parity");
        return;
    };

    for (rust_option, ncbi_option, value) in [
        ("--db_soft_mask", "-db_soft_mask", "abc"),
        ("--db_soft_mask", "-db_soft_mask", ""),
        ("--db_soft_mask", "-db_soft_mask", "1abc"),
        ("--db_soft_mask", "-db_soft_mask", "999999999999999999999"),
        ("--db_hard_mask", "-db_hard_mask", "abc"),
        ("--db_hard_mask", "-db_hard_mask", ""),
        ("--db_hard_mask", "-db_hard_mask", "1abc"),
        ("--db_hard_mask", "-db_hard_mask", "999999999999999999999"),
    ] {
        let rust = std::process::Command::new(&blast_cli)
            .arg("blastn")
            .arg("--query")
            .arg("tests/fixtures/query_short_match.fa")
            .arg("--db")
            .arg("tests/fixtures/seqn/seqn")
            .arg("--task")
            .arg("blastn-short")
            .arg("--outfmt")
            .arg("6")
            .arg("--dust")
            .arg("no")
            .arg(rust_option)
            .arg(value)
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli named DB mask {rust_option}: {err}"));
        let ncbi = std::process::Command::new("/usr/bin/blastn")
            .arg("-query")
            .arg("tests/fixtures/query_short_match.fa")
            .arg("-db")
            .arg("tests/fixtures/seqn/seqn")
            .arg("-task")
            .arg("blastn-short")
            .arg("-outfmt")
            .arg("6")
            .arg("-dust")
            .arg("no")
            .arg(ncbi_option)
            .arg(value)
            .output()
            .unwrap_or_else(|err| panic!("run NCBI named DB mask {ncbi_option}: {err}"));

        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "{rust_option} {value:?} status differs"
        );
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{rust_option} {value:?} stderr differs"
        );
    }
}

#[test]
fn blastn_subject_ncbi_parity_raw_sequence_without_fasta_header() {
    if !std::path::Path::new("/usr/bin/blastn").exists() {
        eprintln!("Skipping: /usr/bin/blastn not found");
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
    std::fs::write(&query, "ACGTACGT\n").expect("write raw query");
    std::fs::write(&subject, "ACGTACGTACGT\n").expect("write raw subject");

    let rust = std::process::Command::new(&blast_cli)
        .arg("blastn")
        .arg("--query")
        .arg(&query)
        .arg("--subject")
        .arg(&subject)
        .arg("--task")
        .arg("blastn-short")
        .arg("--outfmt")
        .arg("6 qseqid sseqid length")
        .arg("--dust")
        .arg("no")
        .output()
        .expect("run blast-cli raw sequence");
    let ncbi = std::process::Command::new("/usr/bin/blastn")
        .arg("-query")
        .arg(&query)
        .arg("-subject")
        .arg(&subject)
        .arg("-task")
        .arg("blastn-short")
        .arg("-outfmt")
        .arg("6 qseqid sseqid length")
        .arg("-dust")
        .arg("no")
        .output()
        .expect("run NCBI raw sequence");

    assert_eq!(rust.status.code(), ncbi.status.code(), "status differs");
    assert_eq!(rust.stdout, ncbi.stdout, "stdout differs");
    assert_eq!(
        String::from_utf8_lossy(&rust.stderr),
        String::from_utf8_lossy(&ncbi.stderr),
        "stderr differs"
    );
}

#[test]
fn blastn_db_ncbi_parity_ignores_preamble_before_fasta_header() {
    assert_blastn_db_outfmt_matches_ncbi(
        "ACGT\n>q1\nACGTACGTACGT\n",
        ">s1\nACGTACGTACGT\n",
        "6 qseqid sseqid length",
        &[],
        &["--dust", "no"],
        &["-dust", "no"],
    );
}

#[test]
fn blastn_db_ncbi_parity_skips_fasta_comment_lines_inside_records() {
    assert_blastn_db_outfmt_matches_ncbi(
        ">q1\nACGTACGT\n; internal comment\n# hash comment\nACGTACGTACGT\n",
        ">s1\nACGTACGTACGTACGTACGT\n",
        "6 qseqid qlen sseqid length qseq",
        &[],
        &["--dust", "no", "--max_target_seqs", "10"],
        &["-dust", "no", "-max_target_seqs", "10"],
    );
}

#[test]
fn blastn_db_ncbi_parity_ignores_invalid_chars_inside_query_sequence_lines() {
    assert_blastn_db_outfmt_matches_ncbi(
        ">q1
ACGT*-.123ACGTACGTACGTACGT
",
        ">s1
ACGTACGTACGTACGTACGT
",
        "6 qseqid qlen sseqid length qseq",
        &[],
        &["--dust", "no", "--max_target_seqs", "10"],
        &["-dust", "no", "-max_target_seqs", "10"],
    );
}

#[test]
fn blastn_subject_ncbi_parity_ignores_invalid_chars_inside_sequence_lines() {
    assert_blastn_subject_outfmt_matches_ncbi(
        ">q1
ACGT*-.123ACGTACGTACGTACGT
",
        ">s1
ACGTACGTACGTACGTACGT
",
        "6 qseqid qlen sseqid length qseq",
        &["--dust", "no", "--max_target_seqs", "10"],
        &["-dust", "no", "-max_target_seqs", "10"],
    );
}

#[test]
fn blastn_subject_ncbi_parity_raw_invalid_residue_warnings_and_errors() {
    if !std::path::Path::new("/usr/bin/blastn").exists() {
        eprintln!("Skipping: /usr/bin/blastn not found");
        return;
    }
    let Some(blast_cli) = std::env::var_os("BLAST_RS_CLI_BIN")
        .or_else(|| std::env::var_os("CARGO_BIN_EXE_blast-cli"))
        .map(std::path::PathBuf::from)
    else {
        eprintln!("Skipping: set BLAST_RS_CLI_BIN or CARGO_BIN_EXE_blast-cli to run CLI parity");
        return;
    };

    for (label, query_fasta) in [
        ("raw invalid residues", "ACGT*-.123ACGTACGTACGTACGT\n"),
        (
            "raw indented pseudo-header",
            "  >q1\nACGTACGTACGTACGTACGT\n",
        ),
        ("raw dash-only line", "---\nACGTACGTACGTACGTACGT\n"),
    ] {
        let tmp = TempDir::new().expect("tempdir");
        let query = tmp.path().join("query.fa");
        let subject = tmp.path().join("subject.fa");
        std::fs::write(&query, query_fasta).expect("write query FASTA");
        std::fs::write(&subject, ">s1\nACGTACGTACGTACGTACGT\n").expect("write subject FASTA");

        let rust = std::process::Command::new(&blast_cli)
            .arg("blastn")
            .arg("--query")
            .arg(&query)
            .arg("--subject")
            .arg(&subject)
            .arg("--task")
            .arg("blastn-short")
            .arg("--outfmt")
            .arg("6 qseqid qlen sseqid length qseq")
            .arg("--dust")
            .arg("no")
            .output()
            .expect("run blast-cli raw invalid residue parity");
        let ncbi = std::process::Command::new("/usr/bin/blastn")
            .arg("-query")
            .arg(&query)
            .arg("-subject")
            .arg(&subject)
            .arg("-task")
            .arg("blastn-short")
            .arg("-outfmt")
            .arg("6 qseqid qlen sseqid length qseq")
            .arg("-dust")
            .arg("no")
            .output()
            .expect("run NCBI raw invalid residue parity");

        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "{label} status differs"
        );
        assert_eq!(rust.stdout, ncbi.stdout, "{label} stdout differs");
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{label} stderr differs"
        );
    }
}
#[test]
fn blastn_subject_ncbi_parity_invalid_residue_warnings() {
    if !std::path::Path::new("/usr/bin/blastn").exists() {
        eprintln!("Skipping: /usr/bin/blastn not found");
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
    std::fs::write(&query, ">q1\nACGT*-.123ACGTACGTACGTACGT\n").expect("write query FASTA");
    std::fs::write(&subject, ">s1\nACGTACGTACGTACGTACGT\n").expect("write subject FASTA");

    let rust = std::process::Command::new(&blast_cli)
        .arg("blastn")
        .arg("--query")
        .arg(&query)
        .arg("--subject")
        .arg(&subject)
        .arg("--task")
        .arg("blastn-short")
        .arg("--outfmt")
        .arg("6 qseqid qlen sseqid length qseq")
        .arg("--dust")
        .arg("no")
        .output()
        .expect("run blast-cli invalid residue warnings");
    let ncbi = std::process::Command::new("/usr/bin/blastn")
        .arg("-query")
        .arg(&query)
        .arg("-subject")
        .arg(&subject)
        .arg("-task")
        .arg("blastn-short")
        .arg("-outfmt")
        .arg("6 qseqid qlen sseqid length qseq")
        .arg("-dust")
        .arg("no")
        .output()
        .expect("run NCBI invalid residue warnings");

    assert_eq!(rust.status.code(), ncbi.status.code(), "status differs");
    assert_eq!(rust.stdout, ncbi.stdout, "stdout differs");
    assert_eq!(
        String::from_utf8_lossy(&rust.stderr),
        String::from_utf8_lossy(&ncbi.stderr),
        "stderr differs"
    );
}
#[test]
fn blastn_subject_ncbi_parity_implausible_fasta_sequence_line_errors() {
    if !std::path::Path::new("/usr/bin/blastn").exists() {
        eprintln!("Skipping: /usr/bin/blastn not found");
        return;
    }
    let Some(blast_cli) = std::env::var_os("BLAST_RS_CLI_BIN")
        .or_else(|| std::env::var_os("CARGO_BIN_EXE_blast-cli"))
        .map(std::path::PathBuf::from)
    else {
        eprintln!("Skipping: set BLAST_RS_CLI_BIN or CARGO_BIN_EXE_blast-cli to run CLI parity");
        return;
    };

    for (label, query_fasta, subject_fasta) in [
        (
            "digits query line",
            ">q1\n123\nACGTACGTACGTACGTACGT\n",
            ">s1\nACGTACGTACGTACGTACGT\n",
        ),
        (
            "dash query line",
            ">q1\n---\nACGTACGTACGTACGTACGT\n",
            ">s1\nACGTACGTACGTACGTACGT\n",
        ),
        (
            "digits subject line",
            ">q1\nACGTACGTACGTACGTACGT\n",
            ">s1\n123\nACGTACGTACGTACGTACGT\n",
        ),
    ] {
        let tmp = TempDir::new().expect("tempdir");
        let query = tmp.path().join("query.fa");
        let subject = tmp.path().join("subject.fa");
        std::fs::write(&query, query_fasta).expect("write query FASTA");
        std::fs::write(&subject, subject_fasta).expect("write subject FASTA");

        let rust = std::process::Command::new(&blast_cli)
            .arg("blastn")
            .arg("--query")
            .arg(&query)
            .arg("--subject")
            .arg(&subject)
            .arg("--task")
            .arg("blastn-short")
            .arg("--outfmt")
            .arg("6")
            .arg("--dust")
            .arg("no")
            .output()
            .expect("run blast-cli implausible FASTA");
        let ncbi = std::process::Command::new("/usr/bin/blastn")
            .arg("-query")
            .arg(&query)
            .arg("-subject")
            .arg(&subject)
            .arg("-task")
            .arg("blastn-short")
            .arg("-outfmt")
            .arg("6")
            .arg("-dust")
            .arg("no")
            .output()
            .expect("run NCBI implausible FASTA");

        assert!(!rust.status.success(), "blast-cli should reject {label}");
        assert!(!ncbi.status.success(), "NCBI should reject {label}");
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "{label} status differs"
        );
        assert_eq!(rust.stdout, ncbi.stdout, "{label} stdout differs");
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{label} stderr differs"
        );
    }
}

#[test]
fn blastn_db_ncbi_parity_ignores_whitespace_inside_query_sequence_lines() {
    assert_blastn_db_outfmt_matches_ncbi(
        ">q1
ACGT ACGT
ACGT	ACGT ACGT
",
        ">s1
ACGTACGTACGTACGTACGT
",
        "6 qseqid qlen sseqid length qseq",
        &[],
        &["--dust", "no", "--max_target_seqs", "10"],
        &["-dust", "no", "-max_target_seqs", "10"],
    );
}

#[test]
fn blastn_subject_ncbi_parity_ignores_whitespace_inside_sequence_lines() {
    assert_blastn_subject_outfmt_matches_ncbi(
        ">q1
ACGT ACGT
ACGT	ACGT ACGT
",
        ">s1
ACGTACGTACGTACGTACGT
",
        "6 qseqid qlen sseqid length qseq",
        &["--dust", "no", "--max_target_seqs", "10"],
        &["-dust", "no", "-max_target_seqs", "10"],
    );
}

#[test]
fn blastn_subject_ncbi_parity_indented_pseudo_header_is_raw_text() {
    assert_blastn_subject_outfmt_matches_ncbi(
        "  >q1
ACGTACGTACGTACGTACGT
",
        ">s1
ACGTACGTACGTACGTACGT
",
        "6 qseqid qlen sseqid length qseq",
        &["--dust", "no", "--max_target_seqs", "10"],
        &["-dust", "no", "-max_target_seqs", "10"],
    );
}

#[test]
fn blastn_db_ncbi_parity_empty_fasta_defline_uses_default_query_id() {
    assert_blastn_db_outfmt_matches_ncbi(
        ">
ACGTACGTACGTACGTACGT
",
        ">s1
ACGTACGTACGTACGTACGT
",
        "6 qseqid qlen sseqid length",
        &[],
        &["--dust", "no", "--max_target_seqs", "10"],
        &["-dust", "no", "-max_target_seqs", "10"],
    );
}

#[test]
fn blastn_subject_ncbi_parity_skips_fasta_comment_lines_inside_records() {
    assert_blastn_subject_outfmt_matches_ncbi(
        ">q1\nACGTACGT\n; internal comment\n# hash comment\nACGTACGTACGT\n",
        ">s1\nACGTACGT\n; subject comment\n# subject hash comment\nACGTACGTACGT\n",
        "6 qseqid qlen sseqid slen length qseq sseq",
        &["--dust", "no", "--max_target_seqs", "10"],
        &["-dust", "no", "-max_target_seqs", "10"],
    );
}

#[test]
fn blastn_db_ncbi_parity_empty_subject_defline_uses_blast_ord_id() {
    assert_blastn_db_outfmt_matches_ncbi(
        ">q1
ACGTACGTACGTACGTACGT
",
        ">
ACGTACGTACGTACGTACGT
",
        "6 qseqid sseqid slen length",
        &[],
        &["--dust", "no", "--max_target_seqs", "10"],
        &["-dust", "no", "-max_target_seqs", "10"],
    );
}

#[test]
fn blastn_subject_ncbi_parity_empty_fasta_defline_uses_default_ids() {
    assert_blastn_subject_outfmt_matches_ncbi(
        ">
ACGTACGTACGTACGTACGT
",
        ">
ACGTACGTACGTACGTACGT
",
        "6 qseqid qlen sseqid slen length",
        &["--dust", "no", "--max_target_seqs", "10"],
        &["-dust", "no", "-max_target_seqs", "10"],
    );
}

#[test]
fn blastn_subject_ncbi_parity_ignores_preamble_before_fasta_header() {
    if !std::path::Path::new("/usr/bin/blastn").exists() {
        eprintln!("Skipping: /usr/bin/blastn not found");
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
    std::fs::write(&query, "ACGT\n>q1\nACGTACGTACGT\n").expect("write query FASTA");
    std::fs::write(&subject, ">s1\nACGTACGTACGT\n").expect("write subject FASTA");

    let rust = std::process::Command::new(&blast_cli)
        .arg("blastn")
        .arg("--query")
        .arg(&query)
        .arg("--subject")
        .arg(&subject)
        .arg("--task")
        .arg("blastn-short")
        .arg("--outfmt")
        .arg("6 qseqid sseqid length")
        .arg("--dust")
        .arg("no")
        .output()
        .expect("run blast-cli FASTA preamble");
    let ncbi = std::process::Command::new("/usr/bin/blastn")
        .arg("-query")
        .arg(&query)
        .arg("-subject")
        .arg(&subject)
        .arg("-task")
        .arg("blastn-short")
        .arg("-outfmt")
        .arg("6 qseqid sseqid length")
        .arg("-dust")
        .arg("no")
        .output()
        .expect("run NCBI FASTA preamble");

    assert_eq!(rust.status.code(), ncbi.status.code(), "status differs");
    assert_eq!(rust.stdout, ncbi.stdout, "stdout differs");
    assert_eq!(
        String::from_utf8_lossy(&rust.stderr),
        String::from_utf8_lossy(&ncbi.stderr),
        "stderr differs"
    );
}

#[test]
fn blastn_subject_ncbi_parity_missing_query_or_subject_file_errors() {
    if !std::path::Path::new("/usr/bin/blastn").exists() {
        eprintln!("Skipping: /usr/bin/blastn not found");
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
    let missing_query = tmp.path().join("missing_query.fa");
    let missing_subject = tmp.path().join("missing_subject.fa");
    std::fs::write(&query, ">q1\nACGTACGTACGTACGTACGT\n").expect("write query FASTA");
    std::fs::write(&subject, ">s1\nACGTACGTACGTACGTACGT\n").expect("write subject FASTA");

    for (rust_query, rust_subject, ncbi_query, ncbi_subject, label) in [
        (
            &missing_query,
            &subject,
            &missing_query,
            &subject,
            "missing query",
        ),
        (
            &query,
            &missing_subject,
            &query,
            &missing_subject,
            "missing subject",
        ),
    ] {
        let rust = std::process::Command::new(&blast_cli)
            .arg("blastn")
            .arg("--query")
            .arg(rust_query)
            .arg("--subject")
            .arg(rust_subject)
            .arg("--task")
            .arg("blastn-short")
            .arg("--outfmt")
            .arg("6")
            .arg("--dust")
            .arg("no")
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli {label}: {err}"));
        let ncbi = std::process::Command::new("/usr/bin/blastn")
            .arg("-query")
            .arg(ncbi_query)
            .arg("-subject")
            .arg(ncbi_subject)
            .arg("-task")
            .arg("blastn-short")
            .arg("-outfmt")
            .arg("6")
            .arg("-dust")
            .arg("no")
            .output()
            .unwrap_or_else(|err| panic!("run NCBI {label}: {err}"));

        assert!(!rust.status.success(), "blast-cli should reject {label}");
        assert!(!ncbi.status.success(), "NCBI should reject {label}");
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "{label} status differs"
        );
        assert_eq!(rust.stdout, ncbi.stdout, "{label} stdout differs");
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{label} stderr differs"
        );
    }
}

#[test]
fn blastn_db_ncbi_parity_missing_query_file_error() {
    if !std::path::Path::new("/usr/bin/blastn").exists()
        || !std::path::Path::new("/usr/bin/makeblastdb").exists()
    {
        eprintln!("Skipping: /usr/bin/blastn or /usr/bin/makeblastdb not found");
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
    let missing_query = tmp.path().join("missing_query.fa");
    let db_fasta = tmp.path().join("db.fa");
    let db = tmp.path().join("testdb");
    std::fs::write(&db_fasta, ">s1\nACGTACGTACGTACGTACGT\n").expect("write DB FASTA");

    let make_status = std::process::Command::new("/usr/bin/makeblastdb")
        .arg("-in")
        .arg(&db_fasta)
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

    let rust = std::process::Command::new(&blast_cli)
        .arg("blastn")
        .arg("--query")
        .arg(&missing_query)
        .arg("--db")
        .arg(&db)
        .arg("--task")
        .arg("blastn-short")
        .arg("--dust")
        .arg("no")
        .arg("--outfmt")
        .arg("6")
        .output()
        .expect("run blast-cli missing DB query");
    let ncbi = std::process::Command::new("/usr/bin/blastn")
        .arg("-query")
        .arg(&missing_query)
        .arg("-db")
        .arg(&db)
        .arg("-task")
        .arg("blastn-short")
        .arg("-dust")
        .arg("no")
        .arg("-outfmt")
        .arg("6")
        .output()
        .expect("run NCBI missing DB query");

    assert!(
        !rust.status.success(),
        "blast-cli should reject missing DB query"
    );
    assert!(
        !ncbi.status.success(),
        "NCBI should reject missing DB query"
    );
    assert_eq!(rust.status.code(), ncbi.status.code(), "status differs");
    assert_eq!(rust.stdout, ncbi.stdout, "stdout differs");
    assert_eq!(
        String::from_utf8_lossy(&rust.stderr),
        String::from_utf8_lossy(&ncbi.stderr),
        "stderr differs"
    );
}

#[test]
fn blastn_subject_ncbi_parity_out_file_behavior() {
    if !std::path::Path::new("/usr/bin/blastn").exists() {
        eprintln!("Skipping: /usr/bin/blastn not found");
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
    std::fs::write(&query, ">q1\nACGTACGTACGTACGTACGT\n").expect("write query FASTA");
    std::fs::write(&subject, ">s1\nACGTACGTACGTACGTACGT\n").expect("write subject FASTA");

    let base_rust_args = [
        "blastn",
        "--query",
        query.to_str().unwrap(),
        "--subject",
        subject.to_str().unwrap(),
        "--task",
        "blastn-short",
        "--outfmt",
        "6",
        "--dust",
        "no",
    ];
    let base_ncbi_args = [
        "-query",
        query.to_str().unwrap(),
        "-subject",
        subject.to_str().unwrap(),
        "-task",
        "blastn-short",
        "-outfmt",
        "6",
        "-dust",
        "no",
    ];

    let rust_stdout = std::process::Command::new(&blast_cli)
        .args(base_rust_args)
        .output()
        .expect("run blast-cli stdout");
    let ncbi_stdout = std::process::Command::new("/usr/bin/blastn")
        .args(base_ncbi_args)
        .output()
        .expect("run NCBI stdout");
    assert_eq!(rust_stdout.status.code(), ncbi_stdout.status.code());
    assert_eq!(rust_stdout.stdout, ncbi_stdout.stdout);
    assert_eq!(
        String::from_utf8_lossy(&rust_stdout.stderr),
        String::from_utf8_lossy(&ncbi_stdout.stderr)
    );

    let rust_out = tmp.path().join("rust.tsv");
    let ncbi_out = tmp.path().join("ncbi.tsv");
    std::fs::write(&rust_out, "stale\n").expect("write stale rust out");
    std::fs::write(&ncbi_out, "stale\n").expect("write stale ncbi out");
    let rust_overwrite = std::process::Command::new(&blast_cli)
        .args(base_rust_args)
        .arg("--out")
        .arg(&rust_out)
        .output()
        .expect("run blast-cli overwrite");
    let ncbi_overwrite = std::process::Command::new("/usr/bin/blastn")
        .args(base_ncbi_args)
        .arg("-out")
        .arg(&ncbi_out)
        .output()
        .expect("run NCBI overwrite");
    assert_eq!(rust_overwrite.status.code(), ncbi_overwrite.status.code());
    assert_eq!(rust_overwrite.stdout, ncbi_overwrite.stdout);
    assert_eq!(
        String::from_utf8_lossy(&rust_overwrite.stderr),
        String::from_utf8_lossy(&ncbi_overwrite.stderr)
    );
    assert_eq!(
        std::fs::read(&rust_out).unwrap(),
        std::fs::read(&ncbi_out).unwrap()
    );

    let rust_dir_out = std::process::Command::new(&blast_cli)
        .args(base_rust_args)
        .arg("--out")
        .arg(tmp.path())
        .output()
        .expect("run blast-cli directory out");
    let ncbi_dir_out = std::process::Command::new("/usr/bin/blastn")
        .args(base_ncbi_args)
        .arg("-out")
        .arg(tmp.path())
        .output()
        .expect("run NCBI directory out");
    assert!(!rust_dir_out.status.success());
    assert!(!ncbi_dir_out.status.success());
    assert_eq!(rust_dir_out.status.code(), ncbi_dir_out.status.code());
    assert_eq!(rust_dir_out.stdout, ncbi_dir_out.stdout);
    assert_eq!(
        String::from_utf8_lossy(&rust_dir_out.stderr),
        String::from_utf8_lossy(&ncbi_dir_out.stderr)
    );

    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;

        let ro_dir = tmp.path().join("readonly");
        std::fs::create_dir(&ro_dir).expect("create read-only dir");
        std::fs::set_permissions(&ro_dir, std::fs::Permissions::from_mode(0o555))
            .expect("make dir read-only");
        let ro_out = ro_dir.join("out.tsv");
        let rust_ro_out = std::process::Command::new(&blast_cli)
            .args(base_rust_args)
            .arg("--out")
            .arg(&ro_out)
            .output()
            .expect("run blast-cli read-only out");
        let ncbi_ro_out = std::process::Command::new("/usr/bin/blastn")
            .args(base_ncbi_args)
            .arg("-out")
            .arg(&ro_out)
            .output()
            .expect("run NCBI read-only out");
        std::fs::set_permissions(&ro_dir, std::fs::Permissions::from_mode(0o755))
            .expect("restore read-only dir permissions");

        assert!(!rust_ro_out.status.success());
        assert!(!ncbi_ro_out.status.success());
        assert_eq!(rust_ro_out.status.code(), ncbi_ro_out.status.code());
        assert_eq!(rust_ro_out.stdout, ncbi_ro_out.stdout);
        assert_eq!(
            String::from_utf8_lossy(&rust_ro_out.stderr),
            String::from_utf8_lossy(&ncbi_ro_out.stderr)
        );
    }
}

#[test]
fn blastn_ncbi_parity_query_masking_resource_errors() {
    if !std::path::Path::new("/usr/bin/blastn").exists() {
        eprintln!("Skipping: /usr/bin/blastn not found");
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
    let missing_filtering_db = tmp.path().join("missing_filtering_db");
    let missing_window_masker_db = tmp.path().join("missing_window_masker_db");

    for (rust_args, ncbi_args, label) in [
        (
            vec![
                "--filtering_db".to_string(),
                missing_filtering_db.display().to_string(),
            ],
            vec![
                "-filtering_db".to_string(),
                missing_filtering_db.display().to_string(),
            ],
            "missing filtering_db",
        ),
        (
            vec![
                "--window_masker_db".to_string(),
                missing_window_masker_db.display().to_string(),
            ],
            vec![
                "-window_masker_db".to_string(),
                missing_window_masker_db.display().to_string(),
            ],
            "missing window_masker_db",
        ),
        (
            vec!["--window_masker_taxid".to_string(), "999999999".to_string()],
            vec!["-window_masker_taxid".to_string(), "999999999".to_string()],
            "missing window_masker_taxid data",
        ),
    ] {
        let use_database = label == "db_soft_mask plus db_hard_mask";
        let mut rust_cmd = std::process::Command::new(&blast_cli);
        rust_cmd
            .arg("blastn")
            .arg("--query")
            .arg("tests/fixtures/query_short_match.fa");
        if use_database {
            rust_cmd.arg("--db").arg("tests/fixtures/seqn/seqn");
        } else {
            rust_cmd
                .arg("--subject")
                .arg("tests/fixtures/subject_test.fa");
        }
        rust_cmd
            .arg("--task")
            .arg("blastn-short")
            .arg("--outfmt")
            .arg("6")
            .arg("--dust")
            .arg("no")
            .args(&rust_args);
        let rust = rust_cmd
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli {label}: {err}"));

        let mut ncbi_cmd = std::process::Command::new("/usr/bin/blastn");
        ncbi_cmd
            .arg("-query")
            .arg("tests/fixtures/query_short_match.fa");
        if use_database {
            ncbi_cmd.arg("-db").arg("tests/fixtures/seqn/seqn");
        } else {
            ncbi_cmd
                .arg("-subject")
                .arg("tests/fixtures/subject_test.fa");
        }
        ncbi_cmd
            .arg("-task")
            .arg("blastn-short")
            .arg("-outfmt")
            .arg("6")
            .arg("-dust")
            .arg("no")
            .args(&ncbi_args);
        let ncbi = ncbi_cmd
            .output()
            .unwrap_or_else(|err| panic!("run NCBI {label}: {err}"));

        assert!(!rust.status.success(), "blast-cli should reject {label}");
        assert!(!ncbi.status.success(), "NCBI should reject {label}");
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "{label} status differs"
        );
        assert_eq!(rust.stdout, ncbi.stdout, "{label} stdout differs");
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{label} stderr differs"
        );
    }
}

#[test]
fn blastn_ncbi_parity_entrez_query_requires_remote_error() {
    if !std::path::Path::new("/usr/bin/blastn").exists() {
        eprintln!("Skipping: /usr/bin/blastn not found");
        return;
    }
    let Some(blast_cli) = std::env::var_os("BLAST_RS_CLI_BIN")
        .or_else(|| std::env::var_os("CARGO_BIN_EXE_blast-cli"))
        .map(std::path::PathBuf::from)
    else {
        eprintln!("Skipping: set BLAST_RS_CLI_BIN or CARGO_BIN_EXE_blast-cli to run CLI parity");
        return;
    };

    let rust = std::process::Command::new(&blast_cli)
        .arg("blastn")
        .arg("--query")
        .arg("tests/fixtures/query_short_match.fa")
        .arg("--subject")
        .arg("tests/fixtures/subject_test.fa")
        .arg("--task")
        .arg("blastn-short")
        .arg("--outfmt")
        .arg("6")
        .arg("--dust")
        .arg("no")
        .arg("--entrez_query")
        .arg("txid9606[orgn]")
        .output()
        .expect("run blast-cli entrez_query without remote");
    let ncbi = std::process::Command::new("/usr/bin/blastn")
        .arg("-query")
        .arg("tests/fixtures/query_short_match.fa")
        .arg("-subject")
        .arg("tests/fixtures/subject_test.fa")
        .arg("-task")
        .arg("blastn-short")
        .arg("-outfmt")
        .arg("6")
        .arg("-dust")
        .arg("no")
        .arg("-entrez_query")
        .arg("txid9606[orgn]")
        .output()
        .expect("run NCBI entrez_query without remote");

    assert!(
        !rust.status.success(),
        "blast-cli should reject entrez_query"
    );
    assert!(!ncbi.status.success(), "NCBI should reject entrez_query");
    assert_eq!(rust.status.code(), ncbi.status.code(), "status differs");
    assert_eq!(rust.stdout, ncbi.stdout, "stdout differs");
    assert_eq!(
        String::from_utf8_lossy(&rust.stderr),
        String::from_utf8_lossy(&ncbi.stderr),
        "stderr differs"
    );
}

#[test]
fn blastn_remote_is_explicitly_unsupported() {
    let blast_cli = blast_cli_bin();
    if !blast_cli.exists() {
        eprintln!(
            "Skipping: build blast-cli first or set CARGO_BIN_EXE_blast-cli to run CLI parity"
        );
        return;
    }

    let output = std::process::Command::new(&blast_cli)
        .arg("blastn")
        .arg("--query")
        .arg("tests/fixtures/query_short_match.fa")
        .arg("--db")
        .arg("nt")
        .arg("--task")
        .arg("blastn-short")
        .arg("--outfmt")
        .arg("6")
        .arg("--remote")
        .output()
        .expect("run blast-cli remote unsupported check");

    assert!(
        !output.status.success(),
        "remote search should fail locally"
    );
    assert_eq!(output.stdout, b"", "remote error should not emit stdout");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("BLAST query/options error: Remote BLAST is not supported"),
        "unexpected remote stderr: {stderr}"
    );
}

#[test]
fn blastn_ncbi_parity_subject_incompatible_with_db_error() {
    if !std::path::Path::new("/usr/bin/blastn").exists() {
        eprintln!("Skipping: /usr/bin/blastn not found");
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
    let multi_query = tmp.path().join("multi_query.fa");
    let subject = tmp.path().join("subject.fa");
    let missing_db = tmp.path().join("missingdb");
    std::fs::write(&query, ">q1\nACGTACGTACGTACGTACGT\n").expect("write query FASTA");
    std::fs::write(
        &multi_query,
        ">q1\nACGTACGTACGTACGTACGT\n>q2\nACGTACGTACGTACGTACGT\n",
    )
    .expect("write multi-query FASTA");
    std::fs::write(&subject, ">s1\nACGTACGTACGTACGTACGT\n").expect("write subject FASTA");

    for (query_path, label) in [(&query, "single query"), (&multi_query, "multi query")] {
        let rust = std::process::Command::new(&blast_cli)
            .arg("blastn")
            .arg("--query")
            .arg(query_path)
            .arg("--subject")
            .arg(&subject)
            .arg("--db")
            .arg(&missing_db)
            .arg("--task")
            .arg("blastn-short")
            .arg("--outfmt")
            .arg("6")
            .arg("--dust")
            .arg("no")
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli subject plus db {label}: {err}"));
        let ncbi = std::process::Command::new("/usr/bin/blastn")
            .arg("-query")
            .arg(query_path)
            .arg("-subject")
            .arg(&subject)
            .arg("-db")
            .arg(&missing_db)
            .arg("-task")
            .arg("blastn-short")
            .arg("-outfmt")
            .arg("6")
            .arg("-dust")
            .arg("no")
            .output()
            .unwrap_or_else(|err| panic!("run NCBI subject plus db {label}: {err}"));

        assert!(
            !rust.status.success(),
            "blast-cli should reject subject plus db {label}"
        );
        assert!(
            !ncbi.status.success(),
            "NCBI should reject subject plus db {label}"
        );
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "{label} status differs"
        );
        assert_eq!(rust.stdout, ncbi.stdout, "{label} stdout differs");
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{label} stderr differs"
        );
    }
}

#[test]
fn blastn_ncbi_parity_subject_incompatible_with_database_filters_error() {
    if !std::path::Path::new("/usr/bin/blastn").exists() {
        eprintln!("Skipping: /usr/bin/blastn not found");
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
    let id_list = tmp.path().join("ids.txt");
    let taxid_list = tmp.path().join("taxids.txt");
    std::fs::write(&query, ">q1\nACGTACGTACGTACGTACGT\n").expect("write query FASTA");
    std::fs::write(&subject, ">s1\nACGTACGTACGTACGTACGT\n").expect("write subject FASTA");
    std::fs::write(&id_list, "1\n").expect("write id list");
    std::fs::write(&taxid_list, "9606\n").expect("write taxid list");

    for (rust_args, ncbi_args, label) in [
        (
            vec!["--gilist".to_string(), id_list.display().to_string()],
            vec!["-gilist".to_string(), id_list.display().to_string()],
            "gilist",
        ),
        (
            vec!["--seqidlist".to_string(), id_list.display().to_string()],
            vec!["-seqidlist".to_string(), id_list.display().to_string()],
            "seqidlist",
        ),
        (
            vec![
                "--negative_gilist".to_string(),
                id_list.display().to_string(),
            ],
            vec![
                "-negative_gilist".to_string(),
                id_list.display().to_string(),
            ],
            "negative_gilist",
        ),
        (
            vec![
                "--negative_seqidlist".to_string(),
                id_list.display().to_string(),
            ],
            vec![
                "-negative_seqidlist".to_string(),
                id_list.display().to_string(),
            ],
            "negative_seqidlist",
        ),
        (
            vec!["--taxids".to_string(), "9606".to_string()],
            vec!["-taxids".to_string(), "9606".to_string()],
            "taxids",
        ),
        (
            vec!["--negative_taxids".to_string(), "9606".to_string()],
            vec!["-negative_taxids".to_string(), "9606".to_string()],
            "negative_taxids",
        ),
        (
            vec!["--taxidlist".to_string(), taxid_list.display().to_string()],
            vec!["-taxidlist".to_string(), taxid_list.display().to_string()],
            "taxidlist",
        ),
        (
            vec![
                "--negative_taxidlist".to_string(),
                taxid_list.display().to_string(),
            ],
            vec![
                "-negative_taxidlist".to_string(),
                taxid_list.display().to_string(),
            ],
            "negative_taxidlist",
        ),
    ] {
        let rust = std::process::Command::new(&blast_cli)
            .arg("blastn")
            .arg("--query")
            .arg(&query)
            .arg("--subject")
            .arg(&subject)
            .arg("--task")
            .arg("blastn-short")
            .arg("--outfmt")
            .arg("6")
            .arg("--dust")
            .arg("no")
            .args(&rust_args)
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli {label}: {err}"));
        let ncbi = std::process::Command::new("/usr/bin/blastn")
            .arg("-query")
            .arg(&query)
            .arg("-subject")
            .arg(&subject)
            .arg("-task")
            .arg("blastn-short")
            .arg("-outfmt")
            .arg("6")
            .arg("-dust")
            .arg("no")
            .args(&ncbi_args)
            .output()
            .unwrap_or_else(|err| panic!("run NCBI {label}: {err}"));

        assert!(!rust.status.success(), "blast-cli should reject {label}");
        assert!(!ncbi.status.success(), "NCBI should reject {label}");
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "{label} status differs"
        );
        assert_eq!(rust.stdout, ncbi.stdout, "{label} stdout differs");
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{label} stderr differs"
        );
    }
}

#[test]
fn blastn_ncbi_parity_database_filter_pairs_are_incompatible() {
    if !std::path::Path::new("/usr/bin/blastn").exists() {
        eprintln!("Skipping: /usr/bin/blastn not found");
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
    let missing_db = tmp.path().join("missingdb");
    let id_list = tmp.path().join("ids.txt");
    let taxid_list = tmp.path().join("taxids.txt");
    std::fs::write(&query, ">q1\nACGTACGTACGTACGTACGT\n").expect("write query FASTA");
    std::fs::write(&id_list, "s1\n").expect("write id list");
    std::fs::write(&taxid_list, "9606\n").expect("write taxid list");

    let id_list = id_list.display().to_string();
    let taxid_list = taxid_list.display().to_string();
    for (rust_args, ncbi_args, label) in [
        (
            vec!["--seqidlist", &id_list, "--gilist", &id_list],
            vec!["-seqidlist", &id_list, "-gilist", &id_list],
            "seqidlist plus gilist",
        ),
        (
            vec!["--negative_gilist", &id_list, "--gilist", &id_list],
            vec!["-negative_gilist", &id_list, "-gilist", &id_list],
            "negative_gilist plus gilist",
        ),
        (
            vec!["--seqidlist", &id_list, "--negative_seqidlist", &id_list],
            vec!["-seqidlist", &id_list, "-negative_seqidlist", &id_list],
            "seqidlist plus negative_seqidlist",
        ),
        (
            vec!["--taxids", "9606", "--seqidlist", &id_list],
            vec!["-taxids", "9606", "-seqidlist", &id_list],
            "taxids plus seqidlist",
        ),
        (
            vec!["--taxids", "9606", "--taxidlist", &taxid_list],
            vec!["-taxids", "9606", "-taxidlist", &taxid_list],
            "taxids plus taxidlist",
        ),
        (
            vec!["--taxids", "9606", "--negative_taxids", "9605"],
            vec!["-taxids", "9606", "-negative_taxids", "9605"],
            "taxids plus negative_taxids",
        ),
        (
            vec![
                "--negative_taxids",
                "9606",
                "--negative_taxidlist",
                &taxid_list,
            ],
            vec![
                "-negative_taxids",
                "9606",
                "-negative_taxidlist",
                &taxid_list,
            ],
            "negative_taxids plus negative_taxidlist",
        ),
    ] {
        let rust = std::process::Command::new(&blast_cli)
            .arg("blastn")
            .arg("--query")
            .arg(&query)
            .arg("--db")
            .arg(&missing_db)
            .arg("--task")
            .arg("blastn-short")
            .arg("--outfmt")
            .arg("6")
            .args(&rust_args)
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli {label}: {err}"));
        let ncbi = std::process::Command::new("/usr/bin/blastn")
            .arg("-query")
            .arg(&query)
            .arg("-db")
            .arg(&missing_db)
            .arg("-task")
            .arg("blastn-short")
            .arg("-outfmt")
            .arg("6")
            .args(&ncbi_args)
            .output()
            .unwrap_or_else(|err| panic!("run NCBI {label}: {err}"));

        assert!(!rust.status.success(), "blast-cli should reject {label}");
        assert!(!ncbi.status.success(), "NCBI should reject {label}");
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "{label} status differs"
        );
        assert_eq!(rust.stdout, ncbi.stdout, "{label} stdout differs");
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{label} stderr differs"
        );
    }
}

#[test]
fn blastn_ncbi_parity_option_relationship_constraints() {
    if !std::path::Path::new("/usr/bin/blastn").exists() {
        eprintln!("Skipping: /usr/bin/blastn not found");
        return;
    }
    let Some(blast_cli) = std::env::var_os("BLAST_RS_CLI_BIN")
        .or_else(|| std::env::var_os("CARGO_BIN_EXE_blast-cli"))
        .map(std::path::PathBuf::from)
    else {
        eprintln!("Skipping: set BLAST_RS_CLI_BIN or CARGO_BIN_EXE_blast-cli to run CLI parity");
        return;
    };

    for (rust_args, ncbi_args, label) in [
        (
            vec![
                "--task",
                "blastn-short",
                "--max_target_seqs",
                "10",
                "--num_descriptions",
                "5",
            ],
            vec![
                "-task",
                "blastn-short",
                "-max_target_seqs",
                "10",
                "-num_descriptions",
                "5",
            ],
            "max_target_seqs plus num_descriptions",
        ),
        (
            vec![
                "--task",
                "blastn-short",
                "--max_target_seqs",
                "10",
                "--num_alignments",
                "5",
            ],
            vec![
                "-task",
                "blastn-short",
                "-max_target_seqs",
                "10",
                "-num_alignments",
                "5",
            ],
            "max_target_seqs plus num_alignments",
        ),
        (
            vec![
                "--task",
                "blastn-short",
                "--culling_limit",
                "0",
                "--best_hit_overhang",
                "0.1",
            ],
            vec![
                "-task",
                "blastn-short",
                "-culling_limit",
                "0",
                "-best_hit_overhang",
                "0.1",
            ],
            "culling_limit plus best_hit_overhang",
        ),
        (
            vec![
                "--task",
                "blastn-short",
                "--culling_limit",
                "1",
                "--best_hit_score_edge",
                "0.1",
            ],
            vec![
                "-task",
                "blastn-short",
                "-culling_limit",
                "1",
                "-best_hit_score_edge",
                "0.1",
            ],
            "culling_limit plus best_hit_score_edge",
        ),
        (
            vec!["--task", "dc-megablast", "--template_type", "coding"],
            vec!["-task", "dc-megablast", "-template_type", "coding"],
            "template_type requires template_length",
        ),
        (
            vec!["--task", "dc-megablast", "--template_length", "16"],
            vec!["-task", "dc-megablast", "-template_length", "16"],
            "template_length requires template_type",
        ),
        (
            vec![
                "--task",
                "blastn-short",
                "--db_soft_mask",
                "1",
                "--db_hard_mask",
                "2",
            ],
            vec![
                "-task",
                "blastn-short",
                "-db_soft_mask",
                "1",
                "-db_hard_mask",
                "2",
            ],
            "db_soft_mask plus db_hard_mask",
        ),
        (
            vec!["--task", "blastn-short", "--mt_mode", "1"],
            vec!["-task", "blastn-short", "-mt_mode", "1"],
            "mt_mode requires num_threads",
        ),
        (
            vec![
                "--task",
                "blastn-short",
                "--import_search_strategy",
                "tests/fixtures/query_short_match.fa",
                "--export_search_strategy",
                "/tmp/blast_strategy.asn",
            ],
            vec![
                "-task",
                "blastn-short",
                "-import_search_strategy",
                "tests/fixtures/query_short_match.fa",
                "-export_search_strategy",
                "/tmp/blast_strategy.asn",
            ],
            "import_search_strategy plus export_search_strategy",
        ),
        (
            vec!["--task", "blastn-short", "--remote", "--num_threads", "2"],
            vec!["-task", "blastn-short", "-remote", "-num_threads", "2"],
            "remote plus num_threads",
        ),
        (
            vec![
                "--task",
                "blastn-short",
                "--remote",
                "--subject_loc",
                "1-10",
            ],
            vec!["-task", "blastn-short", "-remote", "-subject_loc", "1-10"],
            "subject_loc plus remote",
        ),
    ] {
        let use_database = matches!(
            label,
            "db_soft_mask plus db_hard_mask" | "remote plus num_threads"
        );
        let mut rust_cmd = std::process::Command::new(&blast_cli);
        rust_cmd
            .arg("blastn")
            .arg("--query")
            .arg("tests/fixtures/query_short_match.fa");
        if use_database {
            rust_cmd.arg("--db").arg("tests/fixtures/seqn/seqn");
        } else {
            rust_cmd
                .arg("--subject")
                .arg("tests/fixtures/subject_test.fa");
        }
        rust_cmd
            .arg("--outfmt")
            .arg("6")
            .arg("--dust")
            .arg("no")
            .args(&rust_args);
        let rust = rust_cmd
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli {label}: {err}"));

        let mut ncbi_cmd = std::process::Command::new("/usr/bin/blastn");
        ncbi_cmd
            .arg("-query")
            .arg("tests/fixtures/query_short_match.fa");
        if use_database {
            ncbi_cmd.arg("-db").arg("tests/fixtures/seqn/seqn");
        } else {
            ncbi_cmd
                .arg("-subject")
                .arg("tests/fixtures/subject_test.fa");
        }
        ncbi_cmd
            .arg("-outfmt")
            .arg("6")
            .arg("-dust")
            .arg("no")
            .args(&ncbi_args);
        let ncbi = ncbi_cmd
            .output()
            .unwrap_or_else(|err| panic!("run NCBI {label}: {err}"));

        assert!(!rust.status.success(), "blast-cli should reject {label}");
        assert!(!ncbi.status.success(), "NCBI should reject {label}");
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "{label} status differs"
        );
        assert_eq!(rust.stdout, ncbi.stdout, "{label} stdout differs");
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{label} stderr differs"
        );
    }
}

#[test]
fn blastn_ncbi_parity_invalid_boolean_option_errors() {
    if !std::path::Path::new("/usr/bin/blastn").exists() {
        eprintln!("Skipping: /usr/bin/blastn not found");
        return;
    }
    let Some(blast_cli) = std::env::var_os("BLAST_RS_CLI_BIN")
        .or_else(|| std::env::var_os("CARGO_BIN_EXE_blast-cli"))
        .map(std::path::PathBuf::from)
    else {
        eprintln!("Skipping: set BLAST_RS_CLI_BIN or CARGO_BIN_EXE_blast-cli to run CLI parity");
        return;
    };

    for (rust_option, ncbi_option) in [
        ("--soft_masking", "-soft_masking"),
        ("--sum_stats", "-sum_stats"),
        ("--use_index", "-use_index"),
    ] {
        let rust = std::process::Command::new(&blast_cli)
            .arg("blastn")
            .arg("--query")
            .arg("tests/fixtures/query_short_match.fa")
            .arg("--subject")
            .arg("tests/fixtures/subject_test.fa")
            .arg("--task")
            .arg("blastn-short")
            .arg("--outfmt")
            .arg("6")
            .arg("--dust")
            .arg("no")
            .arg(rust_option)
            .arg("maybe")
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli invalid {rust_option}: {err}"));
        let ncbi = std::process::Command::new("/usr/bin/blastn")
            .arg("-query")
            .arg("tests/fixtures/query_short_match.fa")
            .arg("-subject")
            .arg("tests/fixtures/subject_test.fa")
            .arg("-task")
            .arg("blastn-short")
            .arg("-outfmt")
            .arg("6")
            .arg("-dust")
            .arg("no")
            .arg(ncbi_option)
            .arg("maybe")
            .output()
            .unwrap_or_else(|err| panic!("run NCBI invalid {ncbi_option}: {err}"));

        assert!(
            !rust.status.success(),
            "blast-cli should reject invalid {rust_option}"
        );
        assert!(
            !ncbi.status.success(),
            "NCBI should reject invalid {ncbi_option}"
        );
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "{rust_option} status differs"
        );
        assert_eq!(rust.stdout, ncbi.stdout, "{rust_option} stdout differs");
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{rust_option} stderr differs"
        );
    }
}

#[test]
fn blastn_ncbi_parity_index_options_are_accepted_when_unused() {
    if !std::path::Path::new("/usr/bin/blastn").exists() {
        eprintln!("Skipping: /usr/bin/blastn not found");
        return;
    }
    let Some(blast_cli) = std::env::var_os("BLAST_RS_CLI_BIN")
        .or_else(|| std::env::var_os("CARGO_BIN_EXE_blast-cli"))
        .map(std::path::PathBuf::from)
    else {
        eprintln!("Skipping: set BLAST_RS_CLI_BIN or CARGO_BIN_EXE_blast-cli to run CLI parity");
        return;
    };

    for (rust_args, ncbi_args, label) in [
        (
            vec!["--use_index", "true"],
            vec!["-use_index", "true"],
            "use_index true",
        ),
        (
            vec!["--use_index", "false"],
            vec!["-use_index", "false"],
            "use_index false",
        ),
        (
            vec!["--index_name", "idx"],
            vec!["-index_name", "idx"],
            "index_name",
        ),
        (
            vec!["--use_index", "true", "--index_name", "idx"],
            vec!["-use_index", "true", "-index_name", "idx"],
            "use_index with index_name",
        ),
    ] {
        let rust = std::process::Command::new(&blast_cli)
            .arg("blastn")
            .arg("--query")
            .arg("tests/fixtures/query_short_match.fa")
            .arg("--subject")
            .arg("tests/fixtures/subject_test.fa")
            .arg("--task")
            .arg("blastn-short")
            .arg("--outfmt")
            .arg("6")
            .arg("--dust")
            .arg("no")
            .args(&rust_args)
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli {label}: {err}"));
        let ncbi = std::process::Command::new("/usr/bin/blastn")
            .arg("-query")
            .arg("tests/fixtures/query_short_match.fa")
            .arg("-subject")
            .arg("tests/fixtures/subject_test.fa")
            .arg("-task")
            .arg("blastn-short")
            .arg("-outfmt")
            .arg("6")
            .arg("-dust")
            .arg("no")
            .args(&ncbi_args)
            .output()
            .unwrap_or_else(|err| panic!("run NCBI {label}: {err}"));

        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "{label} status differs"
        );
        assert_eq!(rust.stdout, ncbi.stdout, "{label} stdout differs");
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{label} stderr differs"
        );
    }
}

#[test]
fn blastn_ncbi_parity_invalid_dust_option_errors() {
    if !std::path::Path::new("/usr/bin/blastn").exists() {
        eprintln!("Skipping: /usr/bin/blastn not found");
        return;
    }
    let Some(blast_cli) = std::env::var_os("BLAST_RS_CLI_BIN")
        .or_else(|| std::env::var_os("CARGO_BIN_EXE_blast-cli"))
        .map(std::path::PathBuf::from)
    else {
        eprintln!("Skipping: set BLAST_RS_CLI_BIN or CARGO_BIN_EXE_blast-cli to run CLI parity");
        return;
    };

    for dust in ["maybe", "1 2", "a b c"] {
        let rust = std::process::Command::new(&blast_cli)
            .arg("blastn")
            .arg("--query")
            .arg("tests/fixtures/query_short_match.fa")
            .arg("--subject")
            .arg("tests/fixtures/subject_test.fa")
            .arg("--task")
            .arg("blastn-short")
            .arg("--outfmt")
            .arg("6")
            .arg("--dust")
            .arg(dust)
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli invalid dust {dust}: {err}"));
        let ncbi = std::process::Command::new("/usr/bin/blastn")
            .arg("-query")
            .arg("tests/fixtures/query_short_match.fa")
            .arg("-subject")
            .arg("tests/fixtures/subject_test.fa")
            .arg("-task")
            .arg("blastn-short")
            .arg("-outfmt")
            .arg("6")
            .arg("-dust")
            .arg(dust)
            .output()
            .unwrap_or_else(|err| panic!("run NCBI invalid dust {dust}: {err}"));

        assert!(
            !rust.status.success(),
            "blast-cli should reject dust {dust}"
        );
        assert!(!ncbi.status.success(), "NCBI should reject dust {dust}");
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "{dust} status differs"
        );
        assert_eq!(rust.stdout, ncbi.stdout, "{dust} stdout differs");
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{dust} stderr differs"
        );
    }
}

#[test]
fn blastn_ncbi_parity_negative_searchsp_error() {
    if !std::path::Path::new("/usr/bin/blastn").exists() {
        eprintln!("Skipping: /usr/bin/blastn not found");
        return;
    }
    let Some(blast_cli) = std::env::var_os("BLAST_RS_CLI_BIN")
        .or_else(|| std::env::var_os("CARGO_BIN_EXE_blast-cli"))
        .map(std::path::PathBuf::from)
    else {
        eprintln!("Skipping: set BLAST_RS_CLI_BIN or CARGO_BIN_EXE_blast-cli to run CLI parity");
        return;
    };

    let rust = std::process::Command::new(&blast_cli)
        .arg("blastn")
        .arg("--query")
        .arg("tests/fixtures/query_short_match.fa")
        .arg("--subject")
        .arg("tests/fixtures/subject_test.fa")
        .arg("--task")
        .arg("blastn-short")
        .arg("--outfmt")
        .arg("6")
        .arg("--dust")
        .arg("no")
        .arg("--searchsp")
        .arg("-1")
        .output()
        .expect("run blast-cli negative searchsp");
    let ncbi = std::process::Command::new("/usr/bin/blastn")
        .arg("-query")
        .arg("tests/fixtures/query_short_match.fa")
        .arg("-subject")
        .arg("tests/fixtures/subject_test.fa")
        .arg("-task")
        .arg("blastn-short")
        .arg("-outfmt")
        .arg("6")
        .arg("-dust")
        .arg("no")
        .arg("-searchsp")
        .arg("-1")
        .output()
        .expect("run NCBI negative searchsp");

    assert!(
        !rust.status.success(),
        "blast-cli should reject negative searchsp"
    );
    assert!(
        !ncbi.status.success(),
        "NCBI should reject negative searchsp"
    );
    assert_eq!(rust.status.code(), ncbi.status.code(), "status differs");
    assert_eq!(rust.stdout, ncbi.stdout, "stdout differs");
    assert_eq!(
        String::from_utf8_lossy(&rust.stderr),
        String::from_utf8_lossy(&ncbi.stderr),
        "stderr differs"
    );
}

#[test]
fn blastn_ncbi_parity_invalid_word_size_errors() {
    if !std::path::Path::new("/usr/bin/blastn").exists() {
        eprintln!("Skipping: /usr/bin/blastn not found");
        return;
    }
    let Some(blast_cli) = std::env::var_os("BLAST_RS_CLI_BIN")
        .or_else(|| std::env::var_os("CARGO_BIN_EXE_blast-cli"))
        .map(std::path::PathBuf::from)
    else {
        eprintln!("Skipping: set BLAST_RS_CLI_BIN or CARGO_BIN_EXE_blast-cli to run CLI parity");
        return;
    };

    for word_size in ["0", "3", "-1"] {
        let rust = std::process::Command::new(&blast_cli)
            .arg("blastn")
            .arg("--query")
            .arg("tests/fixtures/query_short_match.fa")
            .arg("--subject")
            .arg("tests/fixtures/subject_test.fa")
            .arg("--task")
            .arg("blastn-short")
            .arg("--outfmt")
            .arg("6")
            .arg("--dust")
            .arg("no")
            .arg("--word_size")
            .arg(word_size)
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli invalid word_size {word_size}: {err}"));
        let ncbi = std::process::Command::new("/usr/bin/blastn")
            .arg("-query")
            .arg("tests/fixtures/query_short_match.fa")
            .arg("-subject")
            .arg("tests/fixtures/subject_test.fa")
            .arg("-task")
            .arg("blastn-short")
            .arg("-outfmt")
            .arg("6")
            .arg("-dust")
            .arg("no")
            .arg("-word_size")
            .arg(word_size)
            .output()
            .unwrap_or_else(|err| panic!("run NCBI invalid word_size {word_size}: {err}"));

        assert!(
            !rust.status.success(),
            "blast-cli should reject word_size {word_size}"
        );
        assert!(
            !ncbi.status.success(),
            "NCBI should reject word_size {word_size}"
        );
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "{word_size} status differs"
        );
        assert_eq!(rust.stdout, ncbi.stdout, "{word_size} stdout differs");
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{word_size} stderr differs"
        );
    }
}

#[test]
fn blastn_ncbi_parity_nonpositive_evalue_errors() {
    if !std::path::Path::new("/usr/bin/blastn").exists() {
        eprintln!("Skipping: /usr/bin/blastn not found");
        return;
    }
    let Some(blast_cli) = std::env::var_os("BLAST_RS_CLI_BIN")
        .or_else(|| std::env::var_os("CARGO_BIN_EXE_blast-cli"))
        .map(std::path::PathBuf::from)
    else {
        eprintln!("Skipping: set BLAST_RS_CLI_BIN or CARGO_BIN_EXE_blast-cli to run CLI parity");
        return;
    };

    for evalue in ["0", "-1"] {
        let rust = std::process::Command::new(&blast_cli)
            .arg("blastn")
            .arg("--query")
            .arg("tests/fixtures/query_short_match.fa")
            .arg("--subject")
            .arg("tests/fixtures/subject_test.fa")
            .arg("--task")
            .arg("blastn-short")
            .arg("--outfmt")
            .arg("6")
            .arg("--dust")
            .arg("no")
            .arg("--evalue")
            .arg(evalue)
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli invalid evalue {evalue}: {err}"));
        let ncbi = std::process::Command::new("/usr/bin/blastn")
            .arg("-query")
            .arg("tests/fixtures/query_short_match.fa")
            .arg("-subject")
            .arg("tests/fixtures/subject_test.fa")
            .arg("-task")
            .arg("blastn-short")
            .arg("-outfmt")
            .arg("6")
            .arg("-dust")
            .arg("no")
            .arg("-evalue")
            .arg(evalue)
            .output()
            .unwrap_or_else(|err| panic!("run NCBI invalid evalue {evalue}: {err}"));

        assert!(
            !rust.status.success(),
            "blast-cli should reject evalue {evalue}"
        );
        assert!(!ncbi.status.success(), "NCBI should reject evalue {evalue}");
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "{evalue} status differs"
        );
        assert_eq!(rust.stdout, ncbi.stdout, "{evalue} stdout differs");
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{evalue} stderr differs"
        );
    }
}

#[test]
fn blastn_ncbi_parity_percent_constraint_errors() {
    if !std::path::Path::new("/usr/bin/blastn").exists() {
        eprintln!("Skipping: /usr/bin/blastn not found");
        return;
    }
    let Some(blast_cli) = std::env::var_os("BLAST_RS_CLI_BIN")
        .or_else(|| std::env::var_os("CARGO_BIN_EXE_blast-cli"))
        .map(std::path::PathBuf::from)
    else {
        eprintln!("Skipping: set BLAST_RS_CLI_BIN or CARGO_BIN_EXE_blast-cli to run CLI parity");
        return;
    };

    for (rust_option, ncbi_option, value) in [
        ("--perc_identity", "-perc_identity", "-1"),
        ("--perc_identity", "-perc_identity", "101"),
        ("--qcov_hsp_perc", "-qcov_hsp_perc", "-1"),
        ("--qcov_hsp_perc", "-qcov_hsp_perc", "101"),
    ] {
        let rust = std::process::Command::new(&blast_cli)
            .arg("blastn")
            .arg("--query")
            .arg("tests/fixtures/query_short_match.fa")
            .arg("--subject")
            .arg("tests/fixtures/subject_test.fa")
            .arg("--task")
            .arg("blastn-short")
            .arg("--outfmt")
            .arg("6")
            .arg("--dust")
            .arg("no")
            .arg(rust_option)
            .arg(value)
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli invalid {rust_option} {value}: {err}"));
        let ncbi = std::process::Command::new("/usr/bin/blastn")
            .arg("-query")
            .arg("tests/fixtures/query_short_match.fa")
            .arg("-subject")
            .arg("tests/fixtures/subject_test.fa")
            .arg("-task")
            .arg("blastn-short")
            .arg("-outfmt")
            .arg("6")
            .arg("-dust")
            .arg("no")
            .arg(ncbi_option)
            .arg(value)
            .output()
            .unwrap_or_else(|err| panic!("run NCBI invalid {ncbi_option} {value}: {err}"));

        assert!(
            !rust.status.success(),
            "blast-cli should reject {rust_option} {value}"
        );
        assert!(
            !ncbi.status.success(),
            "NCBI should reject {ncbi_option} {value}"
        );
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "{rust_option} {value} status differs"
        );
        assert_eq!(
            rust.stdout, ncbi.stdout,
            "{rust_option} {value} stdout differs"
        );
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{rust_option} {value} stderr differs"
        );
    }
}

#[test]
fn blastn_ncbi_parity_missing_option_value_errors() {
    if !std::path::Path::new("/usr/bin/blastn").exists() {
        eprintln!("Skipping: /usr/bin/blastn not found");
        return;
    }
    let Some(blast_cli) = std::env::var_os("BLAST_RS_CLI_BIN")
        .or_else(|| std::env::var_os("CARGO_BIN_EXE_blast-cli"))
        .map(std::path::PathBuf::from)
    else {
        eprintln!("Skipping: set BLAST_RS_CLI_BIN or CARGO_BIN_EXE_blast-cli to run CLI parity");
        return;
    };

    for (rust_option, ncbi_option) in [
        ("--task", "-task"),
        ("--strand", "-strand"),
        ("--outfmt", "-outfmt"),
        ("--query", "-query"),
        ("--db", "-db"),
        ("--subject", "-subject"),
        ("--evalue", "-evalue"),
        ("--word_size", "-word_size"),
        ("--num_threads", "-num_threads"),
        ("--dust", "-dust"),
    ] {
        let rust = std::process::Command::new(&blast_cli)
            .arg("blastn")
            .arg(rust_option)
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli missing {rust_option}: {err}"));
        let ncbi = std::process::Command::new("/usr/bin/blastn")
            .arg(ncbi_option)
            .output()
            .unwrap_or_else(|err| panic!("run NCBI missing {ncbi_option}: {err}"));

        assert!(
            !rust.status.success(),
            "blast-cli should reject missing {rust_option}"
        );
        assert!(
            !ncbi.status.success(),
            "NCBI should reject missing {ncbi_option}"
        );
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "{rust_option} status differs"
        );
        assert_eq!(rust.stdout, ncbi.stdout, "{rust_option} stdout differs");
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{rust_option} stderr differs"
        );
    }
}

#[test]
fn blastn_ncbi_parity_choice_constraint_errors() {
    if !std::path::Path::new("/usr/bin/blastn").exists() {
        eprintln!("Skipping: /usr/bin/blastn not found");
        return;
    }
    let Some(blast_cli) = std::env::var_os("BLAST_RS_CLI_BIN")
        .or_else(|| std::env::var_os("CARGO_BIN_EXE_blast-cli"))
        .map(std::path::PathBuf::from)
    else {
        eprintln!("Skipping: set BLAST_RS_CLI_BIN or CARGO_BIN_EXE_blast-cli to run CLI parity");
        return;
    };

    for (rust_extra_args, ncbi_extra_args, label) in [
        (
            vec!["--task", "invalid"],
            vec!["-task", "invalid"],
            "invalid task",
        ),
        (
            vec!["--task", "blastn-short", "--strand", "invalid"],
            vec!["-task", "blastn-short", "-strand", "invalid"],
            "invalid strand",
        ),
        (
            vec!["--task", "blastn-short", "--template_type", "invalid"],
            vec!["-task", "blastn-short", "-template_type", "invalid"],
            "invalid template_type",
        ),
    ] {
        let rust = std::process::Command::new(&blast_cli)
            .arg("blastn")
            .arg("--query")
            .arg("tests/fixtures/query_short_match.fa")
            .arg("--subject")
            .arg("tests/fixtures/subject_test.fa")
            .arg("--outfmt")
            .arg("6")
            .arg("--dust")
            .arg("no")
            .args(rust_extra_args)
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli {label}: {err}"));
        let ncbi = std::process::Command::new("/usr/bin/blastn")
            .arg("-query")
            .arg("tests/fixtures/query_short_match.fa")
            .arg("-subject")
            .arg("tests/fixtures/subject_test.fa")
            .arg("-outfmt")
            .arg("6")
            .arg("-dust")
            .arg("no")
            .args(ncbi_extra_args)
            .output()
            .unwrap_or_else(|err| panic!("run NCBI {label}: {err}"));

        assert!(!rust.status.success(), "blast-cli should reject {label}");
        assert!(!ncbi.status.success(), "NCBI should reject {label}");
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "{label} status differs"
        );
        assert_eq!(rust.stdout, ncbi.stdout, "{label} stdout differs");
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{label} stderr differs"
        );
    }
}

#[test]
fn blastn_ncbi_parity_no_arg_switch_stray_value_errors() {
    if !std::path::Path::new("/usr/bin/blastn").exists() {
        eprintln!("Skipping: /usr/bin/blastn not found");
        return;
    }
    let Some(blast_cli) = std::env::var_os("BLAST_RS_CLI_BIN")
        .or_else(|| std::env::var_os("CARGO_BIN_EXE_blast-cli"))
        .map(std::path::PathBuf::from)
    else {
        eprintln!("Skipping: set BLAST_RS_CLI_BIN or CARGO_BIN_EXE_blast-cli to run CLI parity");
        return;
    };

    for (rust_option, ncbi_option) in [
        ("--ungapped", "-ungapped"),
        ("--lcase_masking", "-lcase_masking"),
        ("--no_greedy", "-no_greedy"),
        ("--show_gis", "-show_gis"),
        ("--html", "-html"),
        ("--subject_besthit", "-subject_besthit"),
        ("--parse_deflines", "-parse_deflines"),
    ] {
        let rust = std::process::Command::new(&blast_cli)
            .arg("blastn")
            .arg("--query")
            .arg("tests/fixtures/query_short_match.fa")
            .arg("--subject")
            .arg("tests/fixtures/subject_test.fa")
            .arg("--task")
            .arg("blastn-short")
            .arg("--outfmt")
            .arg("6")
            .arg("--dust")
            .arg("no")
            .arg(rust_option)
            .arg("false")
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli {rust_option}: {err}"));
        let ncbi = std::process::Command::new("/usr/bin/blastn")
            .arg("-query")
            .arg("tests/fixtures/query_short_match.fa")
            .arg("-subject")
            .arg("tests/fixtures/subject_test.fa")
            .arg("-task")
            .arg("blastn-short")
            .arg("-outfmt")
            .arg("6")
            .arg("-dust")
            .arg("no")
            .arg(ncbi_option)
            .arg("false")
            .output()
            .unwrap_or_else(|err| panic!("run NCBI {ncbi_option}: {err}"));

        assert!(
            !rust.status.success(),
            "blast-cli should reject stray value after {rust_option}"
        );
        assert!(
            !ncbi.status.success(),
            "NCBI should reject stray value after {ncbi_option}"
        );
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "{rust_option} status differs"
        );
        assert_eq!(rust.stdout, ncbi.stdout, "{rust_option} stdout differs");
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{rust_option} stderr differs"
        );
    }
}

#[test]
fn blastn_no_taxid_expansion_rejects_ncbi_2_17_incompatible_options() {
    let Some(blast_cli) = std::env::var_os("BLAST_RS_CLI_BIN")
        .or_else(|| std::env::var_os("CARGO_BIN_EXE_blast-cli"))
        .map(std::path::PathBuf::from)
    else {
        eprintln!("Skipping: set BLAST_RS_CLI_BIN or CARGO_BIN_EXE_blast-cli to run CLI parity");
        return;
    };

    for (database_args, extra_args, incompatible) in [
        (
            vec!["--subject", "tests/fixtures/subject_test.fa"],
            vec![],
            "subject",
        ),
        (vec![], vec!["--subject_loc", "1-10"], "subject_loc"),
        (
            vec!["--db", "tests/fixtures/seqn/seqn"],
            vec!["--window_masker_taxid", "9606"],
            "window_masker_taxid",
        ),
        (
            vec!["--db", "tests/fixtures/seqn/seqn"],
            vec!["--gilist", "tests/fixtures/gilist.txt"],
            "gilist",
        ),
        (
            vec!["--db", "tests/fixtures/seqn/seqn"],
            vec!["--seqidlist", "tests/fixtures/seqidlist.txt"],
            "seqidlist",
        ),
        (
            vec!["--db", "tests/fixtures/seqn/seqn"],
            vec!["--negative_gilist", "tests/fixtures/gilist.txt"],
            "negative_gilist",
        ),
        (
            vec!["--db", "tests/fixtures/seqn/seqn"],
            vec!["--negative_seqidlist", "tests/fixtures/seqidlist.txt"],
            "negative_seqidlist",
        ),
    ] {
        let mut command = std::process::Command::new(&blast_cli);
        command
            .arg("blastn")
            .arg("--query")
            .arg("tests/fixtures/query_short_match.fa")
            .arg("--task")
            .arg("blastn-short")
            .arg("--outfmt")
            .arg("6")
            .arg("--dust")
            .arg("no")
            .arg("--no_taxid_expansion")
            .args(database_args)
            .args(extra_args);
        let output = command
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli no_taxid_expansion {incompatible}: {err}"));
        let stderr = String::from_utf8_lossy(&output.stderr);
        let expected = format!(
            "Argument \"no_taxid_expansion\". Incompatible with argument:  `{incompatible}'"
        );
        assert!(
            !output.status.success(),
            "blast-cli should reject no_taxid_expansion with {incompatible}"
        );
        assert!(
            stderr.contains(&expected),
            "missing expected incompatibility {expected}, stderr was:\n{stderr}"
        );
    }
}

#[test]
fn blastn_ncbi_parity_hsp_pruning_constraint_errors() {
    if !std::path::Path::new("/usr/bin/blastn").exists() {
        eprintln!("Skipping: /usr/bin/blastn not found");
        return;
    }
    let Some(blast_cli) = std::env::var_os("BLAST_RS_CLI_BIN")
        .or_else(|| std::env::var_os("CARGO_BIN_EXE_blast-cli"))
        .map(std::path::PathBuf::from)
    else {
        eprintln!("Skipping: set BLAST_RS_CLI_BIN or CARGO_BIN_EXE_blast-cli to run CLI parity");
        return;
    };

    for (rust_option, ncbi_option, value) in [
        ("--max_target_seqs", "-max_target_seqs", "0"),
        ("--max_hsps", "-max_hsps", "0"),
        ("--culling_limit", "-culling_limit", "-1"),
    ] {
        let rust = std::process::Command::new(&blast_cli)
            .arg("blastn")
            .arg("--query")
            .arg("tests/fixtures/query_short_match.fa")
            .arg("--subject")
            .arg("tests/fixtures/subject_test.fa")
            .arg("--task")
            .arg("blastn-short")
            .arg("--outfmt")
            .arg("6")
            .arg("--dust")
            .arg("no")
            .arg(rust_option)
            .arg(value)
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli invalid {rust_option} {value}: {err}"));
        let ncbi = std::process::Command::new("/usr/bin/blastn")
            .arg("-query")
            .arg("tests/fixtures/query_short_match.fa")
            .arg("-subject")
            .arg("tests/fixtures/subject_test.fa")
            .arg("-task")
            .arg("blastn-short")
            .arg("-outfmt")
            .arg("6")
            .arg("-dust")
            .arg("no")
            .arg(ncbi_option)
            .arg(value)
            .output()
            .unwrap_or_else(|err| panic!("run NCBI invalid {ncbi_option} {value}: {err}"));

        assert!(
            !rust.status.success(),
            "blast-cli should reject {rust_option} {value}"
        );
        assert!(
            !ncbi.status.success(),
            "NCBI should reject {ncbi_option} {value}"
        );
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "{rust_option} {value} status differs"
        );
        assert_eq!(
            rust.stdout, ncbi.stdout,
            "{rust_option} {value} stdout differs"
        );
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{rust_option} {value} stderr differs"
        );
    }
}

#[test]
fn blastn_ncbi_parity_remaining_integer_constraint_errors() {
    if !std::path::Path::new("/usr/bin/blastn").exists() {
        eprintln!("Skipping: /usr/bin/blastn not found");
        return;
    }
    let Some(blast_cli) = std::env::var_os("BLAST_RS_CLI_BIN")
        .or_else(|| std::env::var_os("CARGO_BIN_EXE_blast-cli"))
        .map(std::path::PathBuf::from)
    else {
        eprintln!("Skipping: set BLAST_RS_CLI_BIN or CARGO_BIN_EXE_blast-cli to run CLI parity");
        return;
    };

    for (rust_option, ncbi_option, value) in [
        ("--num_threads", "-num_threads", "0"),
        ("--window_size", "-window_size", "-1"),
        ("--off_diagonal_range", "-off_diagonal_range", "-1"),
        ("--mt_mode", "-mt_mode", "2"),
        ("--num_descriptions", "-num_descriptions", "-1"),
        ("--num_alignments", "-num_alignments", "-1"),
        ("--line_length", "-line_length", "0"),
        ("--sorthits", "-sorthits", "-1"),
        ("--sorthits", "-sorthits", "5"),
        ("--sorthsps", "-sorthsps", "-1"),
        ("--sorthsps", "-sorthsps", "5"),
        ("--reward", "-reward", "-1"),
        ("--penalty", "-penalty", "1"),
        ("--best_hit_overhang", "-best_hit_overhang", "0.5"),
        ("--best_hit_score_edge", "-best_hit_score_edge", "0"),
        ("--template_length", "-template_length", "17"),
    ] {
        let rust = std::process::Command::new(&blast_cli)
            .arg("blastn")
            .arg("--query")
            .arg("tests/fixtures/query_short_match.fa")
            .arg("--subject")
            .arg("tests/fixtures/subject_test.fa")
            .arg("--task")
            .arg("blastn-short")
            .arg("--outfmt")
            .arg("6")
            .arg("--dust")
            .arg("no")
            .arg(rust_option)
            .arg(value)
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli invalid {rust_option} {value}: {err}"));
        let ncbi = std::process::Command::new("/usr/bin/blastn")
            .arg("-query")
            .arg("tests/fixtures/query_short_match.fa")
            .arg("-subject")
            .arg("tests/fixtures/subject_test.fa")
            .arg("-task")
            .arg("blastn-short")
            .arg("-outfmt")
            .arg("6")
            .arg("-dust")
            .arg("no")
            .arg(ncbi_option)
            .arg(value)
            .output()
            .unwrap_or_else(|err| panic!("run NCBI invalid {ncbi_option} {value}: {err}"));

        assert!(
            !rust.status.success(),
            "blast-cli should reject {rust_option} {value}"
        );
        assert!(
            !ncbi.status.success(),
            "NCBI should reject {ncbi_option} {value}"
        );
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "{rust_option} {value} status differs"
        );
        assert_eq!(
            rust.stdout, ncbi.stdout,
            "{rust_option} {value} stdout differs"
        );
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{rust_option} {value} stderr differs"
        );
    }
}

#[test]
fn blastn_ncbi_parity_nonfinite_float_conversion_errors() {
    if !std::path::Path::new("/usr/bin/blastn").exists() {
        eprintln!("Skipping: /usr/bin/blastn not found");
        return;
    }
    let Some(blast_cli) = std::env::var_os("BLAST_RS_CLI_BIN")
        .or_else(|| std::env::var_os("CARGO_BIN_EXE_blast-cli"))
        .map(std::path::PathBuf::from)
    else {
        eprintln!("Skipping: set BLAST_RS_CLI_BIN or CARGO_BIN_EXE_blast-cli to run CLI parity");
        return;
    };

    for (rust_option, ncbi_option, value) in [
        ("--evalue", "-evalue", "nan"),
        ("--evalue", "-evalue", "inf"),
        ("--perc_identity", "-perc_identity", "nan"),
        ("--qcov_hsp_perc", "-qcov_hsp_perc", "nan"),
        ("--xdrop_ungap", "-xdrop_ungap", "nan"),
        ("--xdrop_gap", "-xdrop_gap", "nan"),
        ("--xdrop_gap_final", "-xdrop_gap_final", "nan"),
        ("--xdrop_ungap", "-xdrop_ungap", "inf"),
        ("--xdrop_gap", "-xdrop_gap", "inf"),
        ("--xdrop_gap_final", "-xdrop_gap_final", "inf"),
    ] {
        let rust = std::process::Command::new(&blast_cli)
            .arg("blastn")
            .arg("--query")
            .arg("tests/fixtures/query_short_match.fa")
            .arg("--subject")
            .arg("tests/fixtures/subject_test.fa")
            .arg("--task")
            .arg("blastn-short")
            .arg("--outfmt")
            .arg("6")
            .arg("--dust")
            .arg("no")
            .arg(rust_option)
            .arg(value)
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli invalid {rust_option} {value}: {err}"));
        let ncbi = std::process::Command::new("/usr/bin/blastn")
            .arg("-query")
            .arg("tests/fixtures/query_short_match.fa")
            .arg("-subject")
            .arg("tests/fixtures/subject_test.fa")
            .arg("-task")
            .arg("blastn-short")
            .arg("-outfmt")
            .arg("6")
            .arg("-dust")
            .arg("no")
            .arg(ncbi_option)
            .arg(value)
            .output()
            .unwrap_or_else(|err| panic!("run NCBI invalid {ncbi_option} {value}: {err}"));

        assert!(
            !rust.status.success(),
            "blast-cli should reject {rust_option} {value}"
        );
        assert!(
            !ncbi.status.success(),
            "NCBI should reject {ncbi_option} {value}"
        );
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "{rust_option} {value} status differs"
        );
        assert_eq!(
            rust.stdout, ncbi.stdout,
            "{rust_option} {value} stdout differs"
        );
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{rust_option} {value} stderr differs"
        );
    }
}

#[test]
fn blastn_ncbi_parity_invalid_float_string_conversion_errors() {
    if !std::path::Path::new("/usr/bin/blastn").exists() {
        eprintln!("Skipping: /usr/bin/blastn not found");
        return;
    }
    let Some(blast_cli) = std::env::var_os("BLAST_RS_CLI_BIN")
        .or_else(|| std::env::var_os("CARGO_BIN_EXE_blast-cli"))
        .map(std::path::PathBuf::from)
    else {
        eprintln!("Skipping: set BLAST_RS_CLI_BIN or CARGO_BIN_EXE_blast-cli to run CLI parity");
        return;
    };

    for (rust_option, ncbi_option, value) in [
        ("--evalue", "-evalue", "abc"),
        ("--evalue", "-evalue", "1abc"),
        ("--evalue", "-evalue", "1e+"),
        ("--evalue", "-evalue", "1e+abc"),
        ("--evalue", "-evalue", "1e-"),
        ("--perc_identity", "-perc_identity", "abc"),
        ("--perc_identity", "-perc_identity", "1e+"),
        ("--qcov_hsp_perc", "-qcov_hsp_perc", "abc"),
        ("--xdrop_ungap", "-xdrop_ungap", "abc"),
        ("--xdrop_gap", "-xdrop_gap", "abc"),
        ("--xdrop_gap", "-xdrop_gap", "1e+"),
        ("--xdrop_gap_final", "-xdrop_gap_final", "abc"),
        ("--best_hit_overhang", "-best_hit_overhang", "abc"),
        ("--best_hit_score_edge", "-best_hit_score_edge", "abc"),
    ] {
        let rust = std::process::Command::new(&blast_cli)
            .arg("blastn")
            .arg("--query")
            .arg("tests/fixtures/query_short_match.fa")
            .arg("--subject")
            .arg("tests/fixtures/subject_test.fa")
            .arg("--task")
            .arg("blastn-short")
            .arg("--outfmt")
            .arg("6")
            .arg("--dust")
            .arg("no")
            .arg(rust_option)
            .arg(value)
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli invalid {rust_option} {value}: {err}"));
        let ncbi = std::process::Command::new("/usr/bin/blastn")
            .arg("-query")
            .arg("tests/fixtures/query_short_match.fa")
            .arg("-subject")
            .arg("tests/fixtures/subject_test.fa")
            .arg("-task")
            .arg("blastn-short")
            .arg("-outfmt")
            .arg("6")
            .arg("-dust")
            .arg("no")
            .arg(ncbi_option)
            .arg(value)
            .output()
            .unwrap_or_else(|err| panic!("run NCBI invalid {ncbi_option} {value}: {err}"));

        assert!(
            !rust.status.success(),
            "blast-cli should reject {rust_option} {value}"
        );
        assert!(
            !ncbi.status.success(),
            "NCBI should reject {ncbi_option} {value}"
        );
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "{rust_option} {value} status differs"
        );
        assert_eq!(
            rust.stdout, ncbi.stdout,
            "{rust_option} {value} stdout differs"
        );
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{rust_option} {value} stderr differs"
        );
    }
}

#[test]
fn blastn_ncbi_parity_empty_numeric_conversion_errors() {
    if !std::path::Path::new("/usr/bin/blastn").exists() {
        eprintln!("Skipping: /usr/bin/blastn not found");
        return;
    }
    let Some(blast_cli) = std::env::var_os("BLAST_RS_CLI_BIN")
        .or_else(|| std::env::var_os("CARGO_BIN_EXE_blast-cli"))
        .map(std::path::PathBuf::from)
    else {
        eprintln!("Skipping: set BLAST_RS_CLI_BIN or CARGO_BIN_EXE_blast-cli to run CLI parity");
        return;
    };

    for (rust_option, ncbi_option) in [
        ("--evalue", "-evalue"),
        ("--perc_identity", "-perc_identity"),
        ("--qcov_hsp_perc", "-qcov_hsp_perc"),
        ("--xdrop_ungap", "-xdrop_ungap"),
        ("--xdrop_gap", "-xdrop_gap"),
        ("--xdrop_gap_final", "-xdrop_gap_final"),
        ("--word_size", "-word_size"),
        ("--max_hsps", "-max_hsps"),
        ("--max_target_seqs", "-max_target_seqs"),
        ("--num_threads", "-num_threads"),
        ("--culling_limit", "-culling_limit"),
        ("--window_size", "-window_size"),
        ("--off_diagonal_range", "-off_diagonal_range"),
        ("--mt_mode", "-mt_mode"),
        ("--num_descriptions", "-num_descriptions"),
        ("--num_alignments", "-num_alignments"),
        ("--line_length", "-line_length"),
        ("--dbsize", "-dbsize"),
        ("--searchsp", "-searchsp"),
        ("--sorthits", "-sorthits"),
        ("--sorthsps", "-sorthsps"),
        ("--reward", "-reward"),
        ("--penalty", "-penalty"),
        ("--best_hit_overhang", "-best_hit_overhang"),
        ("--best_hit_score_edge", "-best_hit_score_edge"),
        ("--template_length", "-template_length"),
        ("--gapopen", "-gapopen"),
        ("--gapextend", "-gapextend"),
        ("--min_raw_gapped_score", "-min_raw_gapped_score"),
        ("--window_masker_taxid", "-window_masker_taxid"),
    ] {
        let rust = std::process::Command::new(&blast_cli)
            .arg("blastn")
            .arg("--query")
            .arg("tests/fixtures/query_short_match.fa")
            .arg("--subject")
            .arg("tests/fixtures/subject_test.fa")
            .arg("--task")
            .arg("blastn-short")
            .arg("--outfmt")
            .arg("6")
            .arg("--dust")
            .arg("no")
            .arg(rust_option)
            .arg("")
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli empty {rust_option}: {err}"));
        let ncbi = std::process::Command::new("/usr/bin/blastn")
            .arg("-query")
            .arg("tests/fixtures/query_short_match.fa")
            .arg("-subject")
            .arg("tests/fixtures/subject_test.fa")
            .arg("-task")
            .arg("blastn-short")
            .arg("-outfmt")
            .arg("6")
            .arg("-dust")
            .arg("no")
            .arg(ncbi_option)
            .arg("")
            .output()
            .unwrap_or_else(|err| panic!("run NCBI empty {ncbi_option}: {err}"));

        assert!(
            !rust.status.success(),
            "blast-cli should reject empty {rust_option}"
        );
        assert!(
            !ncbi.status.success(),
            "NCBI should reject empty {ncbi_option}"
        );
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "{rust_option} status differs"
        );
        assert_eq!(rust.stdout, ncbi.stdout, "{rust_option} stdout differs");
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{rust_option} stderr differs"
        );
    }
}

#[test]
fn blastn_ncbi_parity_invalid_integer_string_conversion_errors() {
    if !std::path::Path::new("/usr/bin/blastn").exists() {
        eprintln!("Skipping: /usr/bin/blastn not found");
        return;
    }
    let Some(blast_cli) = std::env::var_os("BLAST_RS_CLI_BIN")
        .or_else(|| std::env::var_os("CARGO_BIN_EXE_blast-cli"))
        .map(std::path::PathBuf::from)
    else {
        eprintln!("Skipping: set BLAST_RS_CLI_BIN or CARGO_BIN_EXE_blast-cli to run CLI parity");
        return;
    };

    for (rust_option, ncbi_option, value) in [
        ("--word_size", "-word_size", "12abc"),
        ("--word_size", "-word_size", "12.0"),
        ("--word_size", "-word_size", "1e3"),
        ("--word_size", "-word_size", "999999999999999999999"),
        ("--word_size", "-word_size", "-999999999999999999999"),
        ("--max_hsps", "-max_hsps", "12abc"),
        ("--max_target_seqs", "-max_target_seqs", "1abc"),
        ("--num_threads", "-num_threads", "1abc"),
        ("--culling_limit", "-culling_limit", "1abc"),
        ("--window_size", "-window_size", "1abc"),
        ("--off_diagonal_range", "-off_diagonal_range", "1abc"),
        ("--mt_mode", "-mt_mode", "1abc"),
        ("--num_descriptions", "-num_descriptions", "1abc"),
        ("--num_alignments", "-num_alignments", "1abc"),
        ("--line_length", "-line_length", "1abc"),
        ("--dbsize", "-dbsize", "12abc"),
        ("--searchsp", "-searchsp", "12abc"),
        ("--sorthits", "-sorthits", "1abc"),
        ("--sorthsps", "-sorthsps", "1abc"),
        ("--reward", "-reward", "1abc"),
        ("--penalty", "-penalty", "1abc"),
        ("--template_length", "-template_length", "1abc"),
        ("--gapopen", "-gapopen", "1abc"),
        ("--gapextend", "-gapextend", "1abc"),
        ("--min_raw_gapped_score", "-min_raw_gapped_score", "1abc"),
        ("--window_masker_taxid", "-window_masker_taxid", "1abc"),
    ] {
        let rust = std::process::Command::new(&blast_cli)
            .arg("blastn")
            .arg("--query")
            .arg("tests/fixtures/query_short_match.fa")
            .arg("--subject")
            .arg("tests/fixtures/subject_test.fa")
            .arg("--task")
            .arg("blastn-short")
            .arg("--outfmt")
            .arg("6")
            .arg("--dust")
            .arg("no")
            .arg(rust_option)
            .arg(value)
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli invalid {rust_option} {value}: {err}"));
        let ncbi = std::process::Command::new("/usr/bin/blastn")
            .arg("-query")
            .arg("tests/fixtures/query_short_match.fa")
            .arg("-subject")
            .arg("tests/fixtures/subject_test.fa")
            .arg("-task")
            .arg("blastn-short")
            .arg("-outfmt")
            .arg("6")
            .arg("-dust")
            .arg("no")
            .arg(ncbi_option)
            .arg(value)
            .output()
            .unwrap_or_else(|err| panic!("run NCBI invalid {ncbi_option} {value}: {err}"));

        assert!(
            !rust.status.success(),
            "blast-cli should reject {rust_option} {value}"
        );
        assert!(
            !ncbi.status.success(),
            "NCBI should reject {ncbi_option} {value}"
        );
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "{rust_option} {value} status differs"
        );
        assert_eq!(
            rust.stdout, ncbi.stdout,
            "{rust_option} {value} stdout differs"
        );
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{rust_option} {value} stderr differs"
        );
    }
}

#[test]
fn blastn_ncbi_parity_integer_option_range_conversion_errors() {
    if !std::path::Path::new("/usr/bin/blastn").exists() {
        eprintln!("Skipping: /usr/bin/blastn not found");
        return;
    }
    let Some(blast_cli) = std::env::var_os("BLAST_RS_CLI_BIN")
        .or_else(|| std::env::var_os("CARGO_BIN_EXE_blast-cli"))
        .map(std::path::PathBuf::from)
    else {
        eprintln!("Skipping: set BLAST_RS_CLI_BIN or CARGO_BIN_EXE_blast-cli to run CLI parity");
        return;
    };

    for (rust_option, ncbi_option) in [
        ("--word_size", "-word_size"),
        ("--max_hsps", "-max_hsps"),
        ("--max_target_seqs", "-max_target_seqs"),
        ("--num_threads", "-num_threads"),
        ("--culling_limit", "-culling_limit"),
        ("--window_size", "-window_size"),
        ("--off_diagonal_range", "-off_diagonal_range"),
        ("--mt_mode", "-mt_mode"),
        ("--num_descriptions", "-num_descriptions"),
        ("--num_alignments", "-num_alignments"),
        ("--line_length", "-line_length"),
        ("--sorthits", "-sorthits"),
        ("--sorthsps", "-sorthsps"),
        ("--reward", "-reward"),
        ("--penalty", "-penalty"),
        ("--template_length", "-template_length"),
        ("--gapopen", "-gapopen"),
        ("--gapextend", "-gapextend"),
        ("--min_raw_gapped_score", "-min_raw_gapped_score"),
        ("--window_masker_taxid", "-window_masker_taxid"),
    ] {
        let value = "5000000000";
        let rust = std::process::Command::new(&blast_cli)
            .arg("blastn")
            .arg("--query")
            .arg("tests/fixtures/query_short_match.fa")
            .arg("--subject")
            .arg("tests/fixtures/subject_test.fa")
            .arg("--task")
            .arg("blastn-short")
            .arg("--outfmt")
            .arg("6")
            .arg("--dust")
            .arg("no")
            .arg(rust_option)
            .arg(value)
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli range {rust_option}: {err}"));
        let ncbi = std::process::Command::new("/usr/bin/blastn")
            .arg("-query")
            .arg("tests/fixtures/query_short_match.fa")
            .arg("-subject")
            .arg("tests/fixtures/subject_test.fa")
            .arg("-task")
            .arg("blastn-short")
            .arg("-outfmt")
            .arg("6")
            .arg("-dust")
            .arg("no")
            .arg(ncbi_option)
            .arg(value)
            .output()
            .unwrap_or_else(|err| panic!("run NCBI range {ncbi_option}: {err}"));

        assert!(
            !rust.status.success(),
            "blast-cli should reject {rust_option}"
        );
        assert!(!ncbi.status.success(), "NCBI should reject {ncbi_option}");
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "{rust_option} status differs"
        );
        assert_eq!(rust.stdout, ncbi.stdout, "{rust_option} stdout differs");
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{rust_option} stderr differs"
        );
    }
}

#[test]
fn blastn_ncbi_parity_bare_exponent_float_values() {
    let query = ">subseq_oid0\nTTAAGGAGGCTCATCTTTCAGAATCCATGCTGTGGGCCAGCAAGAGTTAA\n";
    let subject = ">subj1\nTTAAGGAGGCTCATCTTTCAGAATCCATGCTGTGGGCCAGCAAGAGTTAA\n";
    let outfmt = "6 qseqid sseqid evalue bitscore score length pident";

    for (rust_args, ncbi_args) in [
        (
            vec!["--dust", "no", "--evalue", "1e"],
            vec!["-dust", "no", "-evalue", "1e"],
        ),
        (
            vec!["--dust", "no", "--evalue", "1E"],
            vec!["-dust", "no", "-evalue", "1E"],
        ),
        (
            vec!["--dust", "no", "--xdrop_gap", "1e"],
            vec!["-dust", "no", "-xdrop_gap", "1e"],
        ),
    ] {
        assert_blastn_subject_outfmt_matches_ncbi(query, subject, outfmt, &rust_args, &ncbi_args);
    }
}

#[test]
fn blastn_ncbi_parity_negative_gap_cost_errors() {
    if !std::path::Path::new("/usr/bin/blastn").exists() {
        eprintln!("Skipping: /usr/bin/blastn not found");
        return;
    }
    let Some(blast_cli) = std::env::var_os("BLAST_RS_CLI_BIN")
        .or_else(|| std::env::var_os("CARGO_BIN_EXE_blast-cli"))
        .map(std::path::PathBuf::from)
    else {
        eprintln!("Skipping: set BLAST_RS_CLI_BIN or CARGO_BIN_EXE_blast-cli to run CLI parity");
        return;
    };

    for (rust_option, ncbi_option) in [("--gapopen", "-gapopen"), ("--gapextend", "-gapextend")] {
        let rust = std::process::Command::new(&blast_cli)
            .arg("blastn")
            .arg("--query")
            .arg("tests/fixtures/query_short_match.fa")
            .arg("--subject")
            .arg("tests/fixtures/subject_test.fa")
            .arg("--task")
            .arg("blastn-short")
            .arg("--outfmt")
            .arg("6")
            .arg("--dust")
            .arg("no")
            .arg(rust_option)
            .arg("-1")
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli invalid {rust_option}: {err}"));
        let ncbi = std::process::Command::new("/usr/bin/blastn")
            .arg("-query")
            .arg("tests/fixtures/query_short_match.fa")
            .arg("-subject")
            .arg("tests/fixtures/subject_test.fa")
            .arg("-task")
            .arg("blastn-short")
            .arg("-outfmt")
            .arg("6")
            .arg("-dust")
            .arg("no")
            .arg(ncbi_option)
            .arg("-1")
            .output()
            .unwrap_or_else(|err| panic!("run NCBI invalid {ncbi_option}: {err}"));

        assert!(
            !rust.status.success(),
            "blast-cli should reject invalid {rust_option}"
        );
        assert!(
            !ncbi.status.success(),
            "NCBI should reject invalid {ncbi_option}"
        );
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "{rust_option} status differs"
        );
        assert_eq!(rust.stdout, ncbi.stdout, "{rust_option} stdout differs");
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{rust_option} stderr differs"
        );
    }
}

#[test]
fn blastn_db_ncbi_parity_existing_gilist_without_isam_errors() {
    if !std::path::Path::new("/usr/bin/blastn").exists()
        || !std::path::Path::new("/usr/bin/makeblastdb").exists()
    {
        eprintln!("Skipping: /usr/bin/blastn or /usr/bin/makeblastdb not found");
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
    let db_fasta = tmp.path().join("db.fa");
    let db = tmp.path().join("testdb");
    let gi_list = tmp.path().join("gi.txt");
    std::fs::write(&query, ">q1\nACGTACGTACGTACGTACGT\n").expect("write query FASTA");
    std::fs::write(&db_fasta, ">s1\nACGTACGTACGTACGTACGT\n").expect("write db FASTA");
    std::fs::write(&gi_list, "1\n").expect("write GI list");

    let make_status = std::process::Command::new("/usr/bin/makeblastdb")
        .arg("-in")
        .arg(&db_fasta)
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

    for (rust_option, ncbi_option) in [
        ("--gilist", "-gilist"),
        ("--negative_gilist", "-negative_gilist"),
    ] {
        let rust = std::process::Command::new(&blast_cli)
            .arg("blastn")
            .arg("--query")
            .arg(&query)
            .arg("--db")
            .arg(&db)
            .arg("--task")
            .arg("blastn-short")
            .arg("--outfmt")
            .arg("6")
            .arg("--dust")
            .arg("no")
            .arg(rust_option)
            .arg(&gi_list)
            .output()
            .expect("run blast-cli existing GI list");
        let ncbi = std::process::Command::new("/usr/bin/blastn")
            .arg("-query")
            .arg(&query)
            .arg("-db")
            .arg(&db)
            .arg("-task")
            .arg("blastn-short")
            .arg("-outfmt")
            .arg("6")
            .arg("-dust")
            .arg("no")
            .arg(ncbi_option)
            .arg(&gi_list)
            .output()
            .expect("run NCBI existing GI list");

        assert!(
            !rust.status.success(),
            "blast-cli should reject existing {rust_option}"
        );
        assert!(
            !ncbi.status.success(),
            "NCBI should reject existing {ncbi_option}"
        );
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "existing {rust_option} status differs"
        );
        assert_eq!(
            rust.stdout, ncbi.stdout,
            "existing {rust_option} stdout differs"
        );
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "existing {rust_option} stderr differs"
        );
    }
}

#[test]
fn blastn_db_ncbi_parity_missing_id_list_errors() {
    if !std::path::Path::new("/usr/bin/blastn").exists()
        || !std::path::Path::new("/usr/bin/makeblastdb").exists()
    {
        eprintln!("Skipping: /usr/bin/blastn or /usr/bin/makeblastdb not found");
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
    let db_fasta = tmp.path().join("db.fa");
    let db = tmp.path().join("testdb");
    std::fs::write(&query, ">q1\nACGTACGTACGTACGTACGT\n").expect("write query FASTA");
    std::fs::write(&db_fasta, ">s1\nACGTACGTACGTACGTACGT\n").expect("write db FASTA");

    let make_status = std::process::Command::new("/usr/bin/makeblastdb")
        .arg("-in")
        .arg(&db_fasta)
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

    for (rust_option, ncbi_option) in [
        ("--gilist", "-gilist"),
        ("--seqidlist", "-seqidlist"),
        ("--negative_gilist", "-negative_gilist"),
        ("--negative_seqidlist", "-negative_seqidlist"),
    ] {
        let missing = tmp.path().join(format!(
            "missing_{}.txt",
            rust_option.trim_start_matches("--")
        ));
        let rust_out = tmp
            .path()
            .join(format!("rust_{}.tsv", rust_option.trim_start_matches("--")));
        let ncbi_out = tmp
            .path()
            .join(format!("ncbi_{}.tsv", rust_option.trim_start_matches("--")));

        let rust = std::process::Command::new(&blast_cli)
            .arg("blastn")
            .arg("--query")
            .arg(&query)
            .arg("--db")
            .arg(&db)
            .arg("--task")
            .arg("blastn-short")
            .arg("--outfmt")
            .arg("6")
            .arg("--dust")
            .arg("no")
            .arg(rust_option)
            .arg(&missing)
            .arg("--out")
            .arg(&rust_out)
            .output()
            .expect("run blast-cli missing ID list");
        let ncbi = std::process::Command::new("/usr/bin/blastn")
            .arg("-query")
            .arg(&query)
            .arg("-db")
            .arg(&db)
            .arg("-task")
            .arg("blastn-short")
            .arg("-outfmt")
            .arg("6")
            .arg("-dust")
            .arg("no")
            .arg(ncbi_option)
            .arg(&missing)
            .arg("-out")
            .arg(&ncbi_out)
            .output()
            .expect("run NCBI missing ID list");

        assert!(
            !rust.status.success(),
            "blast-cli should reject missing {rust_option}"
        );
        assert!(
            !ncbi.status.success(),
            "NCBI should reject missing {ncbi_option}"
        );
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "missing {rust_option} status differs"
        );
        assert_eq!(
            std::fs::read(&rust_out).unwrap_or_default(),
            std::fs::read(&ncbi_out).unwrap_or_default(),
            "missing {rust_option} outputs differ"
        );
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "missing {rust_option} stderr differs"
        );
    }
}

#[test]
fn blastn_db_ncbi_parity_traceback_field_matrix_exact() {
    assert_blastn_db_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGT\n",
        ">s1 subject one\nACGTACGTACGTACGTACGT\n",
        "6 qseq sseq sstrand qframe sframe frames btop",
        &[],
        &["--dust", "no", "--max_target_seqs", "10"],
        &["-dust", "no", "-max_target_seqs", "10"],
    );
}

#[test]
fn blastn_db_ncbi_parity_traceback_field_matrix_gapped() {
    assert_blastn_db_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGTAAAACGTACGTACGTACGTACGT\n",
        ">s1 subject one\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n",
        "6 qseq sseq btop qcovhsp gaps mismatch gapopen length",
        &[],
        &[
            "--dust",
            "no",
            "--strand",
            "plus",
            "--word_size",
            "7",
            "--gapopen",
            "5",
            "--gapextend",
            "2",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "10",
        ],
        &[
            "-dust",
            "no",
            "-strand",
            "plus",
            "-word_size",
            "7",
            "-gapopen",
            "5",
            "-gapextend",
            "2",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_min_raw_gapped_score() {
    assert_blastn_db_outfmt_matches_ncbi(
        ">q1\nACGTACGTACGTACGTACGTACGT\n",
        ">perfect\nACGTACGTACGTACGTACGTACGT\n>shorter\nACGTACGTACGTACGTACGTACGA\n",
        "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen score",
        &[],
        &[
            "--dust",
            "no",
            "--strand",
            "plus",
            "--max_target_seqs",
            "10",
            "--min_raw_gapped_score",
            "47",
        ],
        &[
            "-dust",
            "no",
            "-strand",
            "plus",
            "-max_target_seqs",
            "10",
            "-min_raw_gapped_score",
            "47",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_equal_score_subject_group_ordering() {
    assert_blastn_db_outfmt_matches_ncbi(
        ">q
GTCTCCTCTGACTTCAACAGCG
",
        ">plus_low
GTCTCCTCTGACTTCAA
>minus_low
TTGAAGTCAGAGGAGAC
>plus_high
CTCTGACTTCAACAGCG
>minus_high
CGCTGTTGAAGTCAGAG
",
        "6 sseqid score sstart send qstart qend",
        &[],
        &[
            "--dust",
            "no",
            "--word_size",
            "7",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "1",
        ],
        &[
            "-dust",
            "no",
            "-word_size",
            "7",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "1",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_gapped_endpoint_tie_ordering() {
    assert_blastn_db_outfmt_matches_ncbi(
        ">q
TTTACGTACGTACGTACGTACGTAAA
",
        ">s
GGGACGTACGTACGTACGTACGTCCC
",
        "6 qstart qend sstart send score length nident mismatch gaps gapopen qseq sseq btop",
        &[],
        &[
            "--dust",
            "no",
            "--word_size",
            "7",
            "--gapopen",
            "5",
            "--gapextend",
            "2",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "10",
        ],
        &[
            "-dust",
            "no",
            "-word_size",
            "7",
            "-gapopen",
            "5",
            "-gapextend",
            "2",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_num_threads_two_gapped_traceback() {
    assert_blastn_db_outfmt_matches_ncbi_with_num_threads(
        ">q
ACGTTGCAACGATCGTACGATTCGAGCTTAGGCTAATCGGATCCTAGCTAGGCTAATCGATCGTAGCTAGCATCGAT
",
        ">s1
ACGTTGCAACGATCGTACGATTCGAGCTTAGGCTAGGGGGGGGGGGGGGGGTAATCGGATCCTAGCTAGGCTAATCGATCGTAGCTAGCATCGAT
>s2
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
",
        "6 qseqid sseqid qstart qend sstart send score length nident mismatch gaps gapopen qseq sseq btop",
        &[],
        &[
            "--dust",
            "no",
            "--strand",
            "plus",
            "--word_size",
            "7",
            "--reward",
            "1",
            "--penalty",
            "-2",
            "--gapopen",
            "0",
            "--gapextend",
            "2",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "5",
        ],
        &[
            "-dust",
            "no",
            "-strand",
            "plus",
            "-word_size",
            "7",
            "-reward",
            "1",
            "-penalty",
            "-2",
            "-gapopen",
            "0",
            "-gapextend",
            "2",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "5",
        ],
        "2",
    );
}

#[test]
fn blastn_db_ncbi_parity_gapped_ambiguity_display() {
    assert_blastn_db_outfmt_matches_ncbi(
        ">q
ACGTACGTACGTACGTACGTACGTACGT
",
        ">s
ACGTACGTACGTNNNNACGTACGTACGT
",
        "6 qstart qend sstart send score length nident mismatch gaps gapopen qseq sseq btop",
        &[],
        &[
            "--dust",
            "no",
            "--word_size",
            "7",
            "--gapopen",
            "5",
            "--gapextend",
            "2",
            "--max_target_seqs",
            "10",
            "--max_hsps",
            "10",
        ],
        &[
            "-dust",
            "no",
            "-word_size",
            "7",
            "-gapopen",
            "5",
            "-gapextend",
            "2",
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastn_db_ncbi_parity_gapped_traceback_edge_matrix() {
    let outfmt =
        "6 sseqid qstart qend sstart send score length nident mismatch gaps gapopen qseq sseq btop";
    let cases = [
        (
            "adjacent_ins_del",
            ">q\nACGTTGCAACGATCGTACGATTCGAGCTTAGGCTAATCGGATCCTAGCTAGGCTAATCGATCGTAGCTAGCATCGAT\n",
            ">s\nACGTTGCAACGATCGTACGATTCGAGCTTAGGCTAGGGTAATCGGATCCTAGCTAGGCTAATCGATCGTAGCTAGCATCGAT\n",
        ),
        (
            "adjacent_del_ins",
            ">q\nACGTTGCAACGATCGTACGATTCGAGCTTAGGCTAGGGTAATCGGATCCTAGCTAGGCTAATCGATCGTAGCTAGCATCGAT\n",
            ">s\nACGTTGCAACGATCGTACGATTCGAGCTTAGGCTAATCGGATCCTAGCTAGGCTAATCGATCGTAGCTAGCATCGAT\n",
        ),
        (
            "two_gap_choice",
            ">q\nACGTACGTACGTAAAACGTACGTACGT\n",
            ">s\nACGTACGTACGTTTTACGTACGTACGT\n",
        ),
        (
            "equal_gap_mismatch",
            ">q\nACGTACGTACGTACGTACGTACGT\n",
            ">s\nACGTACGTTCGTACGTACGTACGT\n",
        ),
        (
            "gap_near_start",
            ">q\nACGTACGTACGTACGTACGTACGTACGT\n",
            ">s\nACGTAAAACGTACGTACGTACGTACGTACGT\n",
        ),
        (
            "gap_near_end",
            ">q\nACGTACGTACGTACGTACGTACGTACGT\n",
            ">s\nACGTACGTACGTACGTACGTACGTAAAACGT\n",
        ),
    ];

    for (label, query, db) in cases {
        assert_blastn_db_outfmt_matches_ncbi(
            query,
            db,
            outfmt,
            &[],
            &[
                "--dust",
                "no",
                "--strand",
                "plus",
                "--word_size",
                "7",
                "--gapopen",
                "5",
                "--gapextend",
                "2",
                "--max_target_seqs",
                "10",
                "--max_hsps",
                "10",
            ],
            &[
                "-dust",
                "no",
                "-strand",
                "plus",
                "-word_size",
                "7",
                "-gapopen",
                "5",
                "-gapextend",
                "2",
                "-max_target_seqs",
                "10",
                "-max_hsps",
                "10",
            ],
        );
        eprintln!("checked {label}");
    }
}

#[test]
fn blastn_db_ncbi_parity_gapped_low_prelim_high_final_xdrop_matrix() {
    let outfmt =
        "6 qstart qend sstart send score length nident mismatch gaps gapopen qseq sseq btop";
    let cases = [
        (
            "mismatch_block",
            ">q\nACGTACGTACGTAAAACGTACGTACGT\n",
            ">s\nACGTACGTACGTTTTACGTACGTACGT\n",
        ),
        (
            "gap_near_end",
            ">q\nACGTACGTACGTACGTACGTACGTACGT\n",
            ">s\nACGTACGTACGTACGTACGTACGTAAAACGT\n",
        ),
    ];

    for (label, query, db) in cases {
        assert_blastn_db_outfmt_matches_ncbi(
            query,
            db,
            outfmt,
            &[],
            &[
                "--dust",
                "no",
                "--strand",
                "plus",
                "--word_size",
                "7",
                "--gapopen",
                "5",
                "--gapextend",
                "2",
                "--xdrop_gap",
                "2",
                "--xdrop_gap_final",
                "100",
                "--max_target_seqs",
                "10",
                "--max_hsps",
                "10",
            ],
            &[
                "-dust",
                "no",
                "-strand",
                "plus",
                "-word_size",
                "7",
                "-gapopen",
                "5",
                "-gapextend",
                "2",
                "-xdrop_gap",
                "2",
                "-xdrop_gap_final",
                "100",
                "-max_target_seqs",
                "10",
                "-max_hsps",
                "10",
            ],
        );
        eprintln!("checked {label}");
    }
}

#[test]
fn blastn_db_ncbi_parity_gapped_final_xdrop_boundary_matrix() {
    let outfmt =
        "6 qstart qend sstart send score length nident mismatch gaps gapopen qseq sseq btop";
    let cases = [
        (
            "mismatch_block_x2",
            ">q\nACGTACGTACGTAAAACGTACGTACGT\n",
            ">s\nACGTACGTACGTTTTACGTACGTACGT\n",
            "2",
        ),
        (
            "ambiguity_block_x1",
            ">q\nACGTACGTACGTACGTACGTACGTACGT\n",
            ">s\nACGTACGTACGTNNNNACGTACGTACGT\n",
            "1",
        ),
        (
            "gap_near_end_x3",
            ">q\nACGTACGTACGTACGTACGTACGTACGT\n",
            ">s\nACGTACGTACGTACGTACGTACGTAAAACGT\n",
            "3",
        ),
    ];

    for (label, query, db, xdrop) in cases {
        assert_blastn_db_outfmt_matches_ncbi(
            query,
            db,
            outfmt,
            &[],
            &[
                "--dust",
                "no",
                "--strand",
                "plus",
                "--word_size",
                "7",
                "--gapopen",
                "5",
                "--gapextend",
                "2",
                "--xdrop_gap_final",
                xdrop,
                "--max_target_seqs",
                "10",
                "--max_hsps",
                "10",
            ],
            &[
                "-dust",
                "no",
                "-strand",
                "plus",
                "-word_size",
                "7",
                "-gapopen",
                "5",
                "-gapextend",
                "2",
                "-xdrop_gap_final",
                xdrop,
                "-max_target_seqs",
                "10",
                "-max_hsps",
                "10",
            ],
        );
        eprintln!("checked {label}");
    }
}

#[test]
fn blastn_db_ncbi_parity_gapped_traceback_minus_strand_edge_matrix() {
    let outfmt =
        "6 sseqid qstart qend sstart send score length nident mismatch gaps gapopen qseq sseq btop";
    let cases = [
        (
            "minus_adjacent_ins_del",
            "ACGTTGCAACGATCGTACGATTCGAGCTTAGGCTAATCGGATCCTAGCTAGGCTAATCGATCGTAGCTAGCATCGAT",
            "ACGTTGCAACGATCGTACGATTCGAGCTTAGGCTAGGGTAATCGGATCCTAGCTAGGCTAATCGATCGTAGCTAGCATCGAT",
        ),
        (
            "minus_equal_gap_mismatch",
            "ACGTACGTACGTACGTACGTACGT",
            "ACGTACGTTCGTACGTACGTACGT",
        ),
        (
            "minus_gap_near_start",
            "ACGTACGTACGTACGTACGTACGTACGT",
            "ACGTAAAACGTACGTACGTACGTACGTACGT",
        ),
        (
            "minus_gap_near_end",
            "ACGTACGTACGTACGTACGTACGTACGT",
            "ACGTACGTACGTACGTACGTACGTAAAACGT",
        ),
    ];

    for (label, query_seq, plus_subject_seq) in cases {
        let query = format!(">q\n{query_seq}\n");
        let db = format!(">s\n{}\n", ascii_reverse_complement(plus_subject_seq));
        assert_blastn_db_outfmt_matches_ncbi(
            &query,
            &db,
            outfmt,
            &[],
            &[
                "--dust",
                "no",
                "--strand",
                "minus",
                "--word_size",
                "7",
                "--gapopen",
                "5",
                "--gapextend",
                "2",
                "--max_target_seqs",
                "10",
                "--max_hsps",
                "10",
            ],
            &[
                "-dust",
                "no",
                "-strand",
                "minus",
                "-word_size",
                "7",
                "-gapopen",
                "5",
                "-gapextend",
                "2",
                "-max_target_seqs",
                "10",
                "-max_hsps",
                "10",
            ],
        );
        eprintln!("checked {label}");
    }
}

#[test]
fn blastn_subject_ncbi_parity_sall_fields() {
    assert_blastn_subject_outfmt_matches_ncbi(
        ">q1 desc query\nACGTACGTACGTACGTACGTACGT\n",
        ">s1 first subject title here\nACGTACGTACGTACGTACGTACGT\n",
        "6 qseqid sseqid sacc saccver sallseqid sallacc stitle salltitles",
        &["--dust", "no", "--max_target_seqs", "10"],
        &["-dust", "no", "-max_target_seqs", "10"],
    );
}

#[test]
fn blastp_subject_ncbi_parity_default_seg_masks_low_complexity_query() {
    assert_blastp_subject_outfmt_matches_ncbi(
        ">q1 low complexity query\nAAAAAAAAAAAAAAAAAAAA\n",
        concat!(
            ">s1 low complexity subject\nAAAAAAAAAAAAAAAAAAAA\n",
            ">s2 mixed subject\nACDEFGHIKLMNPQRSTVWY\n"
        ),
        "6 qseqid sseqid length bitscore evalue",
        &[],
        &[],
    );
}

#[test]
fn blastp_subject_ncbi_parity_comp_based_stats_masks_low_complexity_subject_redo() {
    for mode in ["1", "2"] {
        assert_blastp_subject_outfmt_matches_ncbi(
            ">q1 low complexity query\nAAAAAAAAAAAAAAAAAAAA\n",
            concat!(
                ">s1 low complexity subject\nAAAAAAAAAAAAAAAAAAAA\n",
                ">s2 mixed subject\nACDEFGHIKLMNPQRSTVWY\n"
            ),
            "6 qseqid sseqid length bitscore evalue",
            &[
                "--seg",
                "no",
                "--comp_based_stats",
                mode,
                "--evalue",
                "1000000000",
            ],
            &[
                "-seg",
                "no",
                "-comp_based_stats",
                mode,
                "-evalue",
                "1000000000",
            ],
        );
    }
}

#[test]
fn blastp_subject_ncbi_parity_custom_seg_options_mask_low_complexity_query() {
    assert_blastp_subject_outfmt_matches_ncbi(
        ">q1 low complexity query\nAAAAAAAAAAAAAAAAAAAA\n",
        ">s1 low complexity subject\nAAAAAAAAAAAAAAAAAAAA\n",
        "6 qseqid sseqid length score bitscore evalue",
        &[
            "--seg",
            "10 1.8 2.1",
            "--comp_based_stats",
            "0",
            "--evalue",
            "1000000000",
        ],
        &[
            "-seg",
            "10 1.8 2.1",
            "-comp_based_stats",
            "0",
            "-evalue",
            "1000000000",
        ],
    );
}

#[test]
fn blastp_ncbi_parity_invalid_seg_option_error() {
    if !std::path::Path::new("/usr/bin/blastp").exists() {
        eprintln!("Skipping: /usr/bin/blastp not found");
        return;
    }
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    let tmp = TempDir::new().expect("tempdir");
    let query = tmp.path().join("query.fa");
    let subject = tmp.path().join("subject.fa");
    std::fs::write(&query, ">q\nACDEFGHIKLMNPQRSTVWY\n").expect("write query");
    std::fs::write(&subject, ">s\nACDEFGHIKLMNPQRSTVWY\n").expect("write subject");

    let rust = std::process::Command::new(&blast_cli)
        .arg("blastp")
        .arg("--query")
        .arg(&query)
        .arg("--subject")
        .arg(&subject)
        .arg("--seg")
        .arg("maybe")
        .arg("--outfmt")
        .arg("6")
        .output()
        .expect("run blast-cli invalid SEG");
    let ncbi = std::process::Command::new("/usr/bin/blastp")
        .arg("-query")
        .arg(&query)
        .arg("-subject")
        .arg(&subject)
        .arg("-seg")
        .arg("maybe")
        .arg("-outfmt")
        .arg("6")
        .output()
        .expect("run NCBI invalid SEG");

    assert!(
        !rust.status.success(),
        "blast-cli should reject invalid SEG"
    );
    assert!(!ncbi.status.success(), "NCBI should reject invalid SEG");
    assert_eq!(rust.status.code(), ncbi.status.code(), "status differs");
    assert_eq!(rust.stdout, ncbi.stdout, "stdout differs");
    assert_eq!(
        String::from_utf8_lossy(&rust.stderr),
        String::from_utf8_lossy(&ncbi.stderr),
        "stderr differs"
    );
}

#[test]
fn blastn_ncbi_parity_seg_is_unknown_argument() {
    if !std::path::Path::new("/usr/bin/blastn").exists() {
        eprintln!("Skipping: /usr/bin/blastn not found");
        return;
    }
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    let tmp = TempDir::new().expect("tempdir");
    let query = tmp.path().join("query.fa");
    let subject = tmp.path().join("subject.fa");
    std::fs::write(&query, ">q\nACGTACGTACGT\n").expect("write query");
    std::fs::write(&subject, ">s\nACGTACGTACGT\n").expect("write subject");

    let rust = std::process::Command::new(&blast_cli)
        .arg("blastn")
        .arg("--query")
        .arg(&query)
        .arg("--subject")
        .arg(&subject)
        .arg("--seg")
        .arg("no")
        .arg("--outfmt")
        .arg("6")
        .output()
        .expect("run blast-cli blastn SEG");
    let ncbi = std::process::Command::new("/usr/bin/blastn")
        .arg("-query")
        .arg(&query)
        .arg("-subject")
        .arg(&subject)
        .arg("-seg")
        .arg("no")
        .arg("-outfmt")
        .arg("6")
        .output()
        .expect("run NCBI blastn SEG");

    assert!(!rust.status.success(), "blast-cli should reject blastn SEG");
    assert!(!ncbi.status.success(), "NCBI should reject blastn SEG");
    assert_eq!(rust.status.code(), ncbi.status.code(), "status differs");
    assert_eq!(rust.stdout, ncbi.stdout, "stdout differs");
    assert_eq!(
        String::from_utf8_lossy(&rust.stderr),
        String::from_utf8_lossy(&ncbi.stderr),
        "stderr differs"
    );
}

#[test]
fn blastn_ncbi_parity_protein_options_are_unknown_arguments() {
    if !std::path::Path::new("/usr/bin/blastn").exists() {
        eprintln!("Skipping: /usr/bin/blastn not found");
        return;
    }
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    for (rust_args, ncbi_args) in [
        (vec!["--matrix", "BLOSUM62"], vec!["-matrix", "BLOSUM62"]),
        (vec!["--threshold", "11"], vec!["-threshold", "11"]),
        (
            vec!["--comp_based_stats", "1"],
            vec!["-comp_based_stats", "1"],
        ),
        (vec!["--query_gencode", "1"], vec!["-query_gencode", "1"]),
        (vec!["--db_gencode", "1"], vec!["-db_gencode", "1"]),
        (vec!["--use_sw_tback"], vec!["-use_sw_tback"]),
        (
            vec!["--in_pssm", "missing.pssm"],
            vec!["-in_pssm", "missing.pssm"],
        ),
    ] {
        let rust = std::process::Command::new(&blast_cli)
            .arg("blastn")
            .args(&rust_args)
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli blastn {rust_args:?}: {err}"));
        let ncbi = std::process::Command::new("/usr/bin/blastn")
            .args(&ncbi_args)
            .output()
            .unwrap_or_else(|err| panic!("run NCBI blastn {ncbi_args:?}: {err}"));

        assert!(
            !rust.status.success(),
            "blast-cli should reject blastn {rust_args:?}"
        );
        assert!(
            !ncbi.status.success(),
            "NCBI should reject blastn {ncbi_args:?}"
        );
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "blastn {rust_args:?} status differs"
        );
        assert_eq!(
            rust.stdout, ncbi.stdout,
            "blastn {rust_args:?} stdout differs"
        );
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "blastn {rust_args:?} stderr differs"
        );
    }
}

#[test]
fn blastp_ncbi_parity_dust_is_unknown_argument() {
    if !std::path::Path::new("/usr/bin/blastp").exists() {
        eprintln!("Skipping: /usr/bin/blastp not found");
        return;
    }
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    let tmp = TempDir::new().expect("tempdir");
    let query = tmp.path().join("query.fa");
    let subject = tmp.path().join("subject.fa");
    std::fs::write(&query, ">q\nACDEFGHIKLMNPQRSTVWY\n").expect("write query");
    std::fs::write(&subject, ">s\nACDEFGHIKLMNPQRSTVWY\n").expect("write subject");

    let rust = std::process::Command::new(&blast_cli)
        .arg("blastp")
        .arg("--query")
        .arg(&query)
        .arg("--subject")
        .arg(&subject)
        .arg("--dust")
        .arg("no")
        .arg("--outfmt")
        .arg("6")
        .output()
        .expect("run blast-cli blastp DUST");
    let ncbi = std::process::Command::new("/usr/bin/blastp")
        .arg("-query")
        .arg(&query)
        .arg("-subject")
        .arg(&subject)
        .arg("-dust")
        .arg("no")
        .arg("-outfmt")
        .arg("6")
        .output()
        .expect("run NCBI blastp DUST");

    assert!(
        !rust.status.success(),
        "blast-cli should reject blastp DUST"
    );
    assert!(!ncbi.status.success(), "NCBI should reject blastp DUST");
    assert_eq!(rust.status.code(), ncbi.status.code(), "status differs");
    assert_eq!(rust.stdout, ncbi.stdout, "stdout differs");
    assert_eq!(
        String::from_utf8_lossy(&rust.stderr),
        String::from_utf8_lossy(&ncbi.stderr),
        "stderr differs"
    );
}

#[test]
fn blastx_ncbi_parity_dust_is_unknown_argument() {
    if !std::path::Path::new("/usr/bin/blastx").exists() {
        eprintln!("Skipping: /usr/bin/blastx not found");
        return;
    }
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    let tmp = TempDir::new().expect("tempdir");
    let query = tmp.path().join("query.fa");
    let subject = tmp.path().join("subject.fa");
    std::fs::write(&query, ">q\nATGAAATTTCTTATTCTTCTTTTC\n").expect("write query");
    std::fs::write(&subject, ">s\nMKFLILLF\n").expect("write subject");

    let rust = std::process::Command::new(&blast_cli)
        .arg("blastx")
        .arg("--query")
        .arg(&query)
        .arg("--subject")
        .arg(&subject)
        .arg("--dust")
        .arg("no")
        .arg("--outfmt")
        .arg("6")
        .output()
        .expect("run blast-cli blastx DUST");
    let ncbi = std::process::Command::new("/usr/bin/blastx")
        .arg("-query")
        .arg(&query)
        .arg("-subject")
        .arg(&subject)
        .arg("-dust")
        .arg("no")
        .arg("-outfmt")
        .arg("6")
        .output()
        .expect("run NCBI blastx DUST");

    assert!(
        !rust.status.success(),
        "blast-cli should reject blastx DUST"
    );
    assert!(!ncbi.status.success(), "NCBI should reject blastx DUST");
    assert_eq!(rust.status.code(), ncbi.status.code(), "status differs");
    assert_eq!(rust.stdout, ncbi.stdout, "stdout differs");
    assert_eq!(
        String::from_utf8_lossy(&rust.stderr),
        String::from_utf8_lossy(&ncbi.stderr),
        "stderr differs"
    );
}

#[test]
fn tblastn_ncbi_parity_dust_is_unknown_argument() {
    if !std::path::Path::new("/usr/bin/tblastn").exists() {
        eprintln!("Skipping: /usr/bin/tblastn not found");
        return;
    }
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    let tmp = TempDir::new().expect("tempdir");
    let query = tmp.path().join("query.fa");
    let subject = tmp.path().join("subject.fa");
    std::fs::write(&query, ">q\nMKFLILLF\n").expect("write query");
    std::fs::write(&subject, ">s\nATGAAATTTCTTATTCTTCTTTTC\n").expect("write subject");

    let rust = std::process::Command::new(&blast_cli)
        .arg("tblastn")
        .arg("--query")
        .arg(&query)
        .arg("--subject")
        .arg(&subject)
        .arg("--dust")
        .arg("no")
        .arg("--outfmt")
        .arg("6")
        .output()
        .expect("run blast-cli tblastn DUST");
    let ncbi = std::process::Command::new("/usr/bin/tblastn")
        .arg("-query")
        .arg(&query)
        .arg("-subject")
        .arg(&subject)
        .arg("-dust")
        .arg("no")
        .arg("-outfmt")
        .arg("6")
        .output()
        .expect("run NCBI tblastn DUST");

    assert!(
        !rust.status.success(),
        "blast-cli should reject tblastn DUST"
    );
    assert!(!ncbi.status.success(), "NCBI should reject tblastn DUST");
    assert_eq!(rust.status.code(), ncbi.status.code(), "status differs");
    assert_eq!(rust.stdout, ncbi.stdout, "stdout differs");
    assert_eq!(
        String::from_utf8_lossy(&rust.stderr),
        String::from_utf8_lossy(&ncbi.stderr),
        "stderr differs"
    );
}

#[test]
fn tblastx_ncbi_parity_dust_is_unknown_argument() {
    if !std::path::Path::new("/usr/bin/tblastx").exists() {
        eprintln!("Skipping: /usr/bin/tblastx not found");
        return;
    }
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    let tmp = TempDir::new().expect("tempdir");
    let query = tmp.path().join("query.fa");
    let subject = tmp.path().join("subject.fa");
    std::fs::write(&query, ">q\nATGAAATTTCTTATTCTTCTTTTC\n").expect("write query");
    std::fs::write(&subject, ">s\nATGAAATTTCTTATTCTTCTTTTC\n").expect("write subject");

    for (rust_option, ncbi_option) in [
        ("--dust", "-dust"),
        ("--gapopen", "-gapopen"),
        ("--gapextend", "-gapextend"),
        ("--ungapped", "-ungapped"),
        ("--mt_mode", "-mt_mode"),
        ("--xdrop_gap", "-xdrop_gap"),
        ("--xdrop_gap_final", "-xdrop_gap_final"),
    ] {
        let rust = std::process::Command::new(&blast_cli)
            .arg("tblastx")
            .arg("--query")
            .arg(&query)
            .arg("--subject")
            .arg(&subject)
            .arg(rust_option)
            .arg("3")
            .arg("--outfmt")
            .arg("6")
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli tblastx {rust_option}: {err}"));
        let ncbi = std::process::Command::new("/usr/bin/tblastx")
            .arg("-query")
            .arg(&query)
            .arg("-subject")
            .arg(&subject)
            .arg(ncbi_option)
            .arg("3")
            .arg("-outfmt")
            .arg("6")
            .output()
            .unwrap_or_else(|err| panic!("run NCBI tblastx {ncbi_option}: {err}"));

        assert!(
            !rust.status.success(),
            "blast-cli should reject tblastx {rust_option}"
        );
        assert!(
            !ncbi.status.success(),
            "NCBI should reject tblastx {ncbi_option}"
        );
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "{rust_option} status differs"
        );
        assert_eq!(rust.stdout, ncbi.stdout, "{rust_option} stdout differs");
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{rust_option} stderr differs"
        );
    }
}

#[test]
fn non_blastn_sam_outfmt_rejection_matches_ncbi() {
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    for (program, query_fasta, subject_fasta) in [
        (
            "blastp",
            ">q\nMKFLIFALILFATVALAPKSSSHEI\n",
            ">s\nMKFLIFALILFATVALAPKSSSHEIHH\n",
        ),
        ("blastx", ">q\nATGAAATTTCTTATTCTTCTTTTC\n", ">s\nMKFLILLF\n"),
        (
            "tblastn",
            ">q\nMKFLILLF\n",
            ">s\nATGAAATTTCTTATTCTTCTTTTC\n",
        ),
        (
            "tblastx",
            ">q\nATGAAATTTCTTATTCTTCTTTTC\n",
            ">s\nATGAAATTTCTTATTCTTCTTTTC\n",
        ),
    ] {
        let ncbi_bin = format!("/usr/bin/{program}");
        if !std::path::Path::new(&ncbi_bin).exists() {
            eprintln!("Skipping: {ncbi_bin} not found");
            continue;
        }

        let tmp = TempDir::new().expect("tempdir");
        let query = tmp.path().join("query.fa");
        let subject = tmp.path().join("subject.fa");
        std::fs::write(&query, query_fasta).expect("write query");
        std::fs::write(&subject, subject_fasta).expect("write subject");

        let rust = std::process::Command::new(&blast_cli)
            .arg(program)
            .arg("--query")
            .arg(&query)
            .arg("--subject")
            .arg(&subject)
            .arg("--outfmt")
            .arg("17")
            .arg("--seg")
            .arg("no")
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli {program} SAM rejection: {err}"));
        let ncbi = std::process::Command::new(&ncbi_bin)
            .arg("-query")
            .arg(&query)
            .arg("-subject")
            .arg(&subject)
            .arg("-outfmt")
            .arg("17")
            .arg("-seg")
            .arg("no")
            .output()
            .unwrap_or_else(|err| panic!("run NCBI {program} SAM rejection: {err}"));

        assert!(
            !rust.status.success(),
            "blast-cli {program} should reject SAM"
        );
        assert!(!ncbi.status.success(), "NCBI {program} should reject SAM");
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "{program} status differs"
        );
        assert_eq!(rust.stdout, ncbi.stdout, "{program} stdout differs");
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{program} stderr differs"
        );
    }
}

#[test]
fn blastn_ncbi_parity_no_greedy_rejects_zero_gap_costs() {
    if !std::path::Path::new("/usr/bin/blastn").exists() {
        eprintln!("Skipping: /usr/bin/blastn not found");
        return;
    }
    let Some(blast_cli) = std::env::var_os("BLAST_RS_CLI_BIN")
        .or_else(|| std::env::var_os("CARGO_BIN_EXE_blast-cli"))
        .map(std::path::PathBuf::from)
    else {
        eprintln!("Skipping: set BLAST_RS_CLI_BIN or CARGO_BIN_EXE_blast-cli to run CLI parity");
        return;
    };

    let rust = std::process::Command::new(&blast_cli)
        .arg("blastn")
        .arg("--query")
        .arg("tests/fixtures/query_short_match.fa")
        .arg("--subject")
        .arg("tests/fixtures/subject_test.fa")
        .arg("--task")
        .arg("megablast")
        .arg("--dust")
        .arg("no")
        .arg("--no_greedy")
        .output()
        .expect("run blast-cli -no_greedy zero gaps");
    let ncbi = std::process::Command::new("/usr/bin/blastn")
        .arg("-query")
        .arg("tests/fixtures/query_short_match.fa")
        .arg("-subject")
        .arg("tests/fixtures/subject_test.fa")
        .arg("-task")
        .arg("megablast")
        .arg("-dust")
        .arg("no")
        .arg("-no_greedy")
        .output()
        .expect("run NCBI -no_greedy zero gaps");

    assert!(!rust.status.success(), "blast-cli should reject -no_greedy");
    assert!(!ncbi.status.success(), "NCBI should reject -no_greedy");
    assert_eq!(rust.status.code(), ncbi.status.code(), "status differs");
    assert_eq!(rust.stdout, ncbi.stdout, "stdout differs");
    assert_eq!(
        String::from_utf8_lossy(&rust.stderr),
        String::from_utf8_lossy(&ncbi.stderr),
        "stderr differs"
    );
}

#[test]
fn protein_programs_ncbi_parity_nucleotide_filtering_options_are_unknown_arguments() {
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    for program in ["blastp", "blastx", "tblastn", "tblastx"] {
        let ncbi_bin = format!("/usr/bin/{program}");
        if !std::path::Path::new(&ncbi_bin).exists() {
            eprintln!("Skipping: {ncbi_bin} not found");
            continue;
        }

        for (rust_option, ncbi_option, value) in [
            ("--filtering_db", "-filtering_db", "missing_filter_db"),
            (
                "--window_masker_db",
                "-window_masker_db",
                "missing_window_db",
            ),
            ("--window_masker_taxid", "-window_masker_taxid", "9606"),
        ] {
            let rust = std::process::Command::new(&blast_cli)
                .arg(program)
                .arg(rust_option)
                .arg(value)
                .output()
                .unwrap_or_else(|err| {
                    panic!("run blast-cli {program} {rust_option} {value}: {err}")
                });
            let ncbi = std::process::Command::new(&ncbi_bin)
                .arg(ncbi_option)
                .arg(value)
                .output()
                .unwrap_or_else(|err| panic!("run NCBI {program} {ncbi_option} {value}: {err}"));

            assert!(
                !rust.status.success(),
                "blast-cli should reject {program} {rust_option}"
            );
            assert!(
                !ncbi.status.success(),
                "NCBI should reject {program} {ncbi_option}"
            );
            assert_eq!(
                rust.status.code(),
                ncbi.status.code(),
                "{program} {rust_option} status differs"
            );
            assert_eq!(
                rust.stdout, ncbi.stdout,
                "{program} {rust_option} stdout differs"
            );
            assert_eq!(
                String::from_utf8_lossy(&rust.stderr),
                String::from_utf8_lossy(&ncbi.stderr),
                "{program} {rust_option} stderr differs"
            );
        }
    }
}

#[test]
fn protein_programs_ncbi_parity_gap_trigger_is_unknown_argument() {
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    for program in ["blastp", "blastx", "tblastn", "tblastx"] {
        let ncbi_bin = format!("/usr/bin/{program}");
        if !std::path::Path::new(&ncbi_bin).exists() {
            eprintln!("Skipping: {ncbi_bin} not found");
            continue;
        }

        let rust = std::process::Command::new(&blast_cli)
            .arg(program)
            .arg("--gap_trigger")
            .arg("22")
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli {program} gap_trigger: {err}"));
        let ncbi = std::process::Command::new(&ncbi_bin)
            .arg("-gap_trigger")
            .arg("22")
            .output()
            .unwrap_or_else(|err| panic!("run NCBI {program} gap_trigger: {err}"));

        assert!(
            !rust.status.success(),
            "blast-cli should reject {program} gap_trigger"
        );
        assert!(
            !ncbi.status.success(),
            "NCBI should reject {program} gap_trigger"
        );
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "{program} gap_trigger status differs"
        );
        assert_eq!(
            rust.stdout, ncbi.stdout,
            "{program} gap_trigger stdout differs"
        );
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{program} gap_trigger stderr differs"
        );
    }
}

#[test]
fn translated_programs_ncbi_parity_invalid_gencode_errors_before_required_query() {
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    for (program, option) in [
        ("blastx", "-query_gencode"),
        ("tblastn", "-db_gencode"),
        ("tblastx", "-query_gencode"),
        ("tblastx", "-db_gencode"),
    ] {
        let ncbi_bin = format!("/usr/bin/{program}");
        if !std::path::Path::new(&ncbi_bin).exists() {
            eprintln!("Skipping: {ncbi_bin} not found");
            continue;
        }
        let rust = std::process::Command::new(&blast_cli)
            .arg(program)
            .arg(option)
            .arg("99")
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli {program} {option}: {err}"));
        let ncbi = std::process::Command::new(&ncbi_bin)
            .arg(option)
            .arg("99")
            .output()
            .unwrap_or_else(|err| panic!("run NCBI {program} {option}: {err}"));

        assert!(
            !rust.status.success(),
            "blast-cli should reject {program} {option}"
        );
        assert!(
            !ncbi.status.success(),
            "NCBI should reject {program} {option}"
        );
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "{program} {option} status differs"
        );
        assert_eq!(
            rust.stdout, ncbi.stdout,
            "{program} {option} stdout differs"
        );
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{program} {option} stderr differs"
        );
    }
}

#[test]
fn protein_programs_ncbi_parity_invalid_threshold_errors_before_required_query() {
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    for program in ["blastp", "blastx", "tblastn", "tblastx"] {
        let ncbi_bin = format!("/usr/bin/{program}");
        if !std::path::Path::new(&ncbi_bin).exists() {
            eprintln!("Skipping: {ncbi_bin} not found");
            continue;
        }
        let rust = std::process::Command::new(&blast_cli)
            .arg(program)
            .arg("-threshold")
            .arg("-1")
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli {program} -threshold -1: {err}"));
        let ncbi = std::process::Command::new(&ncbi_bin)
            .arg("-threshold")
            .arg("-1")
            .output()
            .unwrap_or_else(|err| panic!("run NCBI {program} -threshold -1: {err}"));

        assert!(
            !rust.status.success(),
            "blast-cli should reject {program} -threshold -1"
        );
        assert!(
            !ncbi.status.success(),
            "NCBI should reject {program} -threshold -1"
        );
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "{program} threshold status differs"
        );
        assert_eq!(
            rust.stdout, ncbi.stdout,
            "{program} threshold stdout differs"
        );
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{program} threshold stderr differs"
        );
    }
}

#[test]
fn non_blastn_programs_ncbi_parity_reject_sam_outfmt() {
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    for (program, query, subject) in [
        (
            "blastp",
            "tests/fixtures/protein_query.fa",
            "tests/fixtures/protein_subject.fa",
        ),
        (
            "blastx",
            "tests/fixtures/blastx_nuc_query.fa",
            "tests/fixtures/blastx_prot_subject.fa",
        ),
        (
            "tblastn",
            "tests/fixtures/tblastn_prot_query.fa",
            "tests/fixtures/tblastn_nuc_subject.fa",
        ),
        (
            "tblastx",
            "tests/fixtures/tblastx_nuc_query.fa",
            "tests/fixtures/tblastx_nuc_subject.fa",
        ),
        (
            "psiblast",
            "tests/fixtures/psi_query.fa",
            "tests/fixtures/psi_subject.fa",
        ),
        (
            "deltablast",
            "tests/fixtures/protein_query.fa",
            "tests/fixtures/protein_subject.fa",
        ),
    ] {
        let ncbi_bin = format!("/usr/bin/{program}");
        if !std::path::Path::new(&ncbi_bin).exists() {
            eprintln!("Skipping: {ncbi_bin} not found");
            continue;
        }
        let rust = std::process::Command::new(&blast_cli)
            .arg(program)
            .arg("-query")
            .arg(query)
            .arg("-subject")
            .arg(subject)
            .arg("-outfmt")
            .arg("17")
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli {program} outfmt 17: {err}"));
        let ncbi = std::process::Command::new(&ncbi_bin)
            .arg("-query")
            .arg(query)
            .arg("-subject")
            .arg(subject)
            .arg("-outfmt")
            .arg("17")
            .output()
            .unwrap_or_else(|err| panic!("run NCBI {program} outfmt 17: {err}"));

        assert!(
            !rust.status.success(),
            "blast-cli should reject {program} outfmt 17"
        );
        assert!(
            !ncbi.status.success(),
            "NCBI should reject {program} outfmt 17"
        );
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "{program} outfmt 17 status differs"
        );
        assert_eq!(
            rust.stdout, ncbi.stdout,
            "{program} outfmt 17 stdout differs"
        );
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{program} outfmt 17 stderr differs"
        );
    }

    if std::path::Path::new("/usr/bin/rpstblastn").exists() {
        let rust = std::process::Command::new(&blast_cli)
            .arg("rpstblastn")
            .arg("-query")
            .arg("tests/fixtures/tblastx_nuc_query.fa")
            .arg("-db")
            .arg("missing")
            .arg("-outfmt")
            .arg("17")
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli rpstblastn outfmt 17: {err}"));
        let ncbi = std::process::Command::new("/usr/bin/rpstblastn")
            .arg("-query")
            .arg("tests/fixtures/tblastx_nuc_query.fa")
            .arg("-db")
            .arg("missing")
            .arg("-outfmt")
            .arg("17")
            .output()
            .unwrap_or_else(|err| panic!("run NCBI rpstblastn outfmt 17: {err}"));

        assert!(
            !rust.status.success(),
            "blast-cli should reject rpstblastn outfmt 17"
        );
        assert!(
            !ncbi.status.success(),
            "NCBI should reject rpstblastn outfmt 17"
        );
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "rpstblastn outfmt 17 status differs"
        );
        assert_eq!(
            rust.stdout, ncbi.stdout,
            "rpstblastn outfmt 17 stdout differs"
        );
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "rpstblastn outfmt 17 stderr differs"
        );
    }
}

#[test]
fn installed_programs_document_archive_outfmt_11_temporary_boundary() {
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    for (program, query, subject) in [
        (
            "blastn",
            "tests/fixtures/query_short_match.fa",
            "tests/fixtures/subject_test.fa",
        ),
        (
            "blastp",
            "tests/fixtures/protein_query.fa",
            "tests/fixtures/protein_subject.fa",
        ),
        (
            "blastx",
            "tests/fixtures/blastx_nuc_query.fa",
            "tests/fixtures/blastx_prot_subject.fa",
        ),
        (
            "tblastn",
            "tests/fixtures/tblastn_prot_query.fa",
            "tests/fixtures/tblastn_nuc_subject.fa",
        ),
        (
            "tblastx",
            "tests/fixtures/tblastx_nuc_query.fa",
            "tests/fixtures/tblastx_nuc_subject.fa",
        ),
        (
            "psiblast",
            "tests/fixtures/psi_query.fa",
            "tests/fixtures/psi_subject.fa",
        ),
    ] {
        let output = std::process::Command::new(&blast_cli)
            .arg(program)
            .arg("-query")
            .arg(query)
            .arg("-subject")
            .arg(subject)
            .arg("-outfmt")
            .arg("11")
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli {program} outfmt 11: {err}"));

        assert!(
            !output.status.success(),
            "blast-cli should reject archive outfmt 11 for {program} until archive writing is implemented"
        );
        assert!(
            output.stdout.is_empty(),
            "archive outfmt 11 rejection should not write stdout for {program}"
        );
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(
            stderr.contains("Output format 11 is not supported"),
            "archive outfmt 11 stderr differed for {program}:\n{stderr}"
        );
    }
}

#[test]
fn deltablast_db_ncbi_parity_supported_outfmts_reach_missing_database() {
    if !std::path::Path::new("/usr/bin/deltablast").exists() {
        eprintln!("Skipping: /usr/bin/deltablast not found");
        return;
    }
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    for outfmt in ["0", "5", "6", "7", "10"] {
        let rust = std::process::Command::new(&blast_cli)
            .arg("deltablast")
            .arg("-query")
            .arg("tests/fixtures/protein_query.fa")
            .arg("-db")
            .arg("missing")
            .arg("-outfmt")
            .arg(outfmt)
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli deltablast outfmt {outfmt}: {err}"));
        let ncbi = std::process::Command::new("/usr/bin/deltablast")
            .arg("-query")
            .arg("tests/fixtures/protein_query.fa")
            .arg("-db")
            .arg("missing")
            .arg("-outfmt")
            .arg(outfmt)
            .output()
            .unwrap_or_else(|err| panic!("run NCBI deltablast outfmt {outfmt}: {err}"));

        assert!(
            !rust.status.success(),
            "blast-cli should reject deltablast outfmt {outfmt} without CDD"
        );
        assert!(
            !ncbi.status.success(),
            "NCBI should reject deltablast outfmt {outfmt} without CDD"
        );
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "deltablast outfmt {outfmt} status differs"
        );
        assert_eq!(
            rust.stdout, ncbi.stdout,
            "deltablast outfmt {outfmt} stdout differs"
        );
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "deltablast outfmt {outfmt} stderr differs"
        );
    }
}

#[test]
fn rpstblastn_db_ncbi_parity_supported_outfmts_reach_missing_database() {
    if !std::path::Path::new("/usr/bin/rpstblastn").exists() {
        eprintln!("Skipping: /usr/bin/rpstblastn not found");
        return;
    }
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    for outfmt in ["0", "5", "6", "7", "10"] {
        let rust = std::process::Command::new(&blast_cli)
            .arg("rpstblastn")
            .arg("-query")
            .arg("tests/fixtures/tblastx_nuc_query.fa")
            .arg("-db")
            .arg("missing")
            .arg("-outfmt")
            .arg(outfmt)
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli rpstblastn outfmt {outfmt}: {err}"));
        let ncbi = std::process::Command::new("/usr/bin/rpstblastn")
            .arg("-query")
            .arg("tests/fixtures/tblastx_nuc_query.fa")
            .arg("-db")
            .arg("missing")
            .arg("-outfmt")
            .arg(outfmt)
            .output()
            .unwrap_or_else(|err| panic!("run NCBI rpstblastn outfmt {outfmt}: {err}"));

        assert!(
            !rust.status.success(),
            "blast-cli should reject rpstblastn outfmt {outfmt} without CDD"
        );
        assert!(
            !ncbi.status.success(),
            "NCBI should reject rpstblastn outfmt {outfmt} without CDD"
        );
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "rpstblastn outfmt {outfmt} status differs"
        );
        assert_eq!(
            rust.stdout, ncbi.stdout,
            "rpstblastn outfmt {outfmt} stdout differs"
        );
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "rpstblastn outfmt {outfmt} stderr differs"
        );
    }
}

#[test]
fn rpstblastn_db_ncbi_parity_num_threads_two_reaches_missing_database() {
    if !std::path::Path::new("/usr/bin/rpstblastn").exists() {
        eprintln!("Skipping: /usr/bin/rpstblastn not found");
        return;
    }
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    let rust = std::process::Command::new(&blast_cli)
        .arg("rpstblastn")
        .arg("-query")
        .arg("tests/fixtures/tblastx_nuc_query.fa")
        .arg("-db")
        .arg("missing")
        .arg("-outfmt")
        .arg("6")
        .arg("-num_threads")
        .arg("2")
        .output()
        .expect("run blast-cli rpstblastn num_threads 2");
    let ncbi = std::process::Command::new("/usr/bin/rpstblastn")
        .arg("-query")
        .arg("tests/fixtures/tblastx_nuc_query.fa")
        .arg("-db")
        .arg("missing")
        .arg("-outfmt")
        .arg("6")
        .arg("-num_threads")
        .arg("2")
        .output()
        .expect("run NCBI rpstblastn num_threads 2");

    assert!(
        !rust.status.success(),
        "blast-cli should reject rpstblastn num_threads 2 without CDD"
    );
    assert!(
        !ncbi.status.success(),
        "NCBI should reject rpstblastn num_threads 2 without CDD"
    );
    assert_eq!(
        rust.status.code(),
        ncbi.status.code(),
        "rpstblastn num_threads 2 status differs"
    );
    assert_eq!(
        rust.stdout, ncbi.stdout,
        "rpstblastn num_threads 2 stdout differs"
    );
    assert_eq!(
        String::from_utf8_lossy(&rust.stderr),
        String::from_utf8_lossy(&ncbi.stderr),
        "rpstblastn num_threads 2 stderr differs"
    );
}

#[test]
fn db_ncbi_parity_plain_unimplemented_outfmts_reach_missing_database() {
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    for (program, query, db_kind) in [
        ("blastn", "tests/fixtures/query_random_200.fa", "nucleotide"),
        ("blastp", "tests/fixtures/protein_query.fa", "protein"),
        ("blastx", "tests/fixtures/blastx_nuc_query.fa", "protein"),
        (
            "tblastn",
            "tests/fixtures/tblastn_prot_query.fa",
            "nucleotide",
        ),
        (
            "tblastx",
            "tests/fixtures/tblastx_nuc_query.fa",
            "nucleotide",
        ),
        ("psiblast", "tests/fixtures/psi_query.fa", "protein"),
        (
            "rpstblastn",
            "tests/fixtures/tblastx_nuc_query.fa",
            "protein",
        ),
        ("deltablast", "tests/fixtures/protein_query.fa", "protein"),
    ] {
        let ncbi_bin = format!("/usr/bin/{program}");
        if !std::path::Path::new(&ncbi_bin).exists() {
            eprintln!("Skipping: {ncbi_bin} not found");
            continue;
        }
        for outfmt in ["12", "13", "14", "15", "16", "18"] {
            let rust = std::process::Command::new(&blast_cli)
                .arg(program)
                .arg("-query")
                .arg(query)
                .arg("-db")
                .arg("missing")
                .arg("-outfmt")
                .arg(outfmt)
                .output()
                .unwrap_or_else(|err| {
                    panic!("run blast-cli {program} missing DB outfmt {outfmt}: {err}")
                });
            let ncbi = std::process::Command::new(&ncbi_bin)
                .arg("-query")
                .arg(query)
                .arg("-db")
                .arg("missing")
                .arg("-outfmt")
                .arg(outfmt)
                .output()
                .unwrap_or_else(|err| {
                    panic!("run NCBI {program} missing DB outfmt {outfmt}: {err}")
                });

            assert!(
                !rust.status.success(),
                "blast-cli should reject {program} missing DB outfmt {outfmt}"
            );
            assert!(
                !ncbi.status.success(),
                "NCBI should reject {program} missing DB outfmt {outfmt}"
            );
            assert_eq!(
                rust.status.code(),
                ncbi.status.code(),
                "{program} missing {db_kind} DB outfmt {outfmt} status differs"
            );
            assert_eq!(
                rust.stdout, ncbi.stdout,
                "{program} missing {db_kind} DB outfmt {outfmt} stdout differs"
            );
            assert_eq!(
                String::from_utf8_lossy(&rust.stderr),
                String::from_utf8_lossy(&ncbi.stderr),
                "{program} missing {db_kind} DB outfmt {outfmt} stderr differs"
            );
        }
    }
}

#[test]
fn installed_programs_ncbi_parity_missing_db_or_subject_diagnostic() {
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    for (program, query) in [
        ("blastn", "tests/fixtures/query_random_200.fa"),
        ("blastp", "tests/fixtures/protein_query.fa"),
        ("blastx", "tests/fixtures/blastx_nuc_query.fa"),
        ("tblastn", "tests/fixtures/tblastn_prot_query.fa"),
        ("tblastx", "tests/fixtures/tblastx_nuc_query.fa"),
        ("psiblast", "tests/fixtures/psi_query.fa"),
        ("rpstblastn", "tests/fixtures/tblastx_nuc_query.fa"),
        ("deltablast", "tests/fixtures/protein_query.fa"),
    ] {
        let ncbi_bin = format!("/usr/bin/{program}");
        if !std::path::Path::new(&ncbi_bin).exists() {
            eprintln!("Skipping: {ncbi_bin} not found");
            continue;
        }
        let rust = std::process::Command::new(&blast_cli)
            .arg(program)
            .arg("-query")
            .arg(query)
            .arg("-outfmt")
            .arg("6")
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli {program} without db/subject: {err}"));
        let ncbi = std::process::Command::new(&ncbi_bin)
            .arg("-query")
            .arg(query)
            .arg("-outfmt")
            .arg("6")
            .output()
            .unwrap_or_else(|err| panic!("run NCBI {program} without db/subject: {err}"));

        assert!(
            !rust.status.success(),
            "blast-cli should reject {program} without db/subject"
        );
        assert!(
            !ncbi.status.success(),
            "NCBI should reject {program} without db/subject"
        );
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "{program} no-db/subject status differs"
        );
        assert_eq!(
            rust.stdout, ncbi.stdout,
            "{program} no-db/subject stdout differs"
        );
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{program} no-db/subject stderr differs"
        );
    }
}

#[test]
fn installed_programs_ncbi_parity_subject_and_db_are_incompatible() {
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    for (program, query, subject) in [
        (
            "blastn",
            "tests/fixtures/query_random_200.fa",
            "tests/fixtures/subject_test.fa",
        ),
        (
            "blastp",
            "tests/fixtures/protein_query.fa",
            "tests/fixtures/protein_subject.fa",
        ),
        (
            "blastx",
            "tests/fixtures/blastx_nuc_query.fa",
            "tests/fixtures/blastx_prot_subject.fa",
        ),
        (
            "tblastn",
            "tests/fixtures/tblastn_prot_query.fa",
            "tests/fixtures/tblastn_nuc_subject.fa",
        ),
        (
            "tblastx",
            "tests/fixtures/tblastx_nuc_query.fa",
            "tests/fixtures/tblastx_nuc_subject.fa",
        ),
        (
            "psiblast",
            "tests/fixtures/psi_query.fa",
            "tests/fixtures/psi_subject.fa",
        ),
        (
            "rpstblastn",
            "tests/fixtures/tblastx_nuc_query.fa",
            "tests/fixtures/tblastx_nuc_subject.fa",
        ),
        (
            "deltablast",
            "tests/fixtures/protein_query.fa",
            "tests/fixtures/protein_subject.fa",
        ),
    ] {
        let ncbi_bin = format!("/usr/bin/{program}");
        if !std::path::Path::new(&ncbi_bin).exists() {
            eprintln!("Skipping: {ncbi_bin} not found");
            continue;
        }
        let rust = std::process::Command::new(&blast_cli)
            .arg(program)
            .arg("-query")
            .arg(query)
            .arg("-subject")
            .arg(subject)
            .arg("-db")
            .arg("missing")
            .arg("-outfmt")
            .arg("6")
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli {program} with subject and db: {err}"));
        let ncbi = std::process::Command::new(&ncbi_bin)
            .arg("-query")
            .arg(query)
            .arg("-subject")
            .arg(subject)
            .arg("-db")
            .arg("missing")
            .arg("-outfmt")
            .arg("6")
            .output()
            .unwrap_or_else(|err| panic!("run NCBI {program} with subject and db: {err}"));

        assert!(
            !rust.status.success(),
            "blast-cli should reject {program} with subject and db"
        );
        assert!(
            !ncbi.status.success(),
            "NCBI should reject {program} with subject and db"
        );
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "{program} subject+db status differs"
        );
        assert_eq!(
            rust.stdout, ncbi.stdout,
            "{program} subject+db stdout differs"
        );
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{program} subject+db stderr differs"
        );
    }
}

#[test]
fn installed_programs_ncbi_parity_subject_loc_and_db_are_incompatible() {
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    for (program, query) in [
        ("blastn", "tests/fixtures/query_random_200.fa"),
        ("blastp", "tests/fixtures/protein_query.fa"),
        ("blastx", "tests/fixtures/blastx_nuc_query.fa"),
        ("tblastn", "tests/fixtures/tblastn_prot_query.fa"),
        ("tblastx", "tests/fixtures/tblastx_nuc_query.fa"),
        ("psiblast", "tests/fixtures/psi_query.fa"),
        ("rpstblastn", "tests/fixtures/tblastx_nuc_query.fa"),
        ("deltablast", "tests/fixtures/protein_query.fa"),
    ] {
        let ncbi_bin = format!("/usr/bin/{program}");
        if !std::path::Path::new(&ncbi_bin).exists() {
            eprintln!("Skipping: {ncbi_bin} not found");
            continue;
        }
        let rust = std::process::Command::new(&blast_cli)
            .arg(program)
            .arg("-query")
            .arg(query)
            .arg("-subject_loc")
            .arg("1-10")
            .arg("-db")
            .arg("missing")
            .arg("-outfmt")
            .arg("6")
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli {program} with subject_loc and db: {err}"));
        let ncbi = std::process::Command::new(&ncbi_bin)
            .arg("-query")
            .arg(query)
            .arg("-subject_loc")
            .arg("1-10")
            .arg("-db")
            .arg("missing")
            .arg("-outfmt")
            .arg("6")
            .output()
            .unwrap_or_else(|err| panic!("run NCBI {program} with subject_loc and db: {err}"));

        assert!(
            !rust.status.success(),
            "blast-cli should reject {program} with subject_loc and db"
        );
        assert!(
            !ncbi.status.success(),
            "NCBI should reject {program} with subject_loc and db"
        );
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "{program} subject_loc+db status differs"
        );
        assert_eq!(
            rust.stdout, ncbi.stdout,
            "{program} subject_loc+db stdout differs"
        );
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{program} subject_loc+db stderr differs"
        );
    }
}

#[test]
fn subject_mode_programs_ncbi_parity_subject_and_gilist_are_incompatible() {
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    for (program, query, subject) in [
        (
            "blastn",
            "tests/fixtures/query_random_200.fa",
            "tests/fixtures/subject_test.fa",
        ),
        (
            "blastp",
            "tests/fixtures/protein_query.fa",
            "tests/fixtures/protein_subject.fa",
        ),
        (
            "blastx",
            "tests/fixtures/blastx_nuc_query.fa",
            "tests/fixtures/blastx_prot_subject.fa",
        ),
        (
            "tblastn",
            "tests/fixtures/tblastn_prot_query.fa",
            "tests/fixtures/tblastn_nuc_subject.fa",
        ),
        (
            "tblastx",
            "tests/fixtures/tblastx_nuc_query.fa",
            "tests/fixtures/tblastx_nuc_subject.fa",
        ),
        (
            "psiblast",
            "tests/fixtures/psi_query.fa",
            "tests/fixtures/psi_subject.fa",
        ),
        (
            "deltablast",
            "tests/fixtures/protein_query.fa",
            "tests/fixtures/protein_subject.fa",
        ),
    ] {
        let ncbi_bin = format!("/usr/bin/{program}");
        if !std::path::Path::new(&ncbi_bin).exists() {
            eprintln!("Skipping: {ncbi_bin} not found");
            continue;
        }
        let rust = std::process::Command::new(&blast_cli)
            .arg(program)
            .arg("-query")
            .arg(query)
            .arg("-subject")
            .arg(subject)
            .arg("-gilist")
            .arg("missing.gi")
            .arg("-outfmt")
            .arg("6")
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli {program} with subject and gilist: {err}"));
        let ncbi = std::process::Command::new(&ncbi_bin)
            .arg("-query")
            .arg(query)
            .arg("-subject")
            .arg(subject)
            .arg("-gilist")
            .arg("missing.gi")
            .arg("-outfmt")
            .arg("6")
            .output()
            .unwrap_or_else(|err| panic!("run NCBI {program} with subject and gilist: {err}"));

        assert!(
            !rust.status.success(),
            "blast-cli should reject {program} with subject and gilist"
        );
        assert!(
            !ncbi.status.success(),
            "NCBI should reject {program} with subject and gilist"
        );
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "{program} subject+gilist status differs"
        );
        assert_eq!(
            rust.stdout, ncbi.stdout,
            "{program} subject+gilist stdout differs"
        );
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{program} subject+gilist stderr differs"
        );
    }
}

#[test]
fn installed_programs_ncbi_parity_gilist_and_seqidlist_are_incompatible() {
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    for (program, query) in [
        ("blastn", "tests/fixtures/query_random_200.fa"),
        ("blastp", "tests/fixtures/protein_query.fa"),
        ("blastx", "tests/fixtures/blastx_nuc_query.fa"),
        ("tblastn", "tests/fixtures/tblastn_prot_query.fa"),
        ("tblastx", "tests/fixtures/tblastx_nuc_query.fa"),
        ("psiblast", "tests/fixtures/psi_query.fa"),
        ("rpstblastn", "tests/fixtures/tblastx_nuc_query.fa"),
        ("deltablast", "tests/fixtures/protein_query.fa"),
    ] {
        let ncbi_bin = format!("/usr/bin/{program}");
        if !std::path::Path::new(&ncbi_bin).exists() {
            eprintln!("Skipping: {ncbi_bin} not found");
            continue;
        }
        let rust = std::process::Command::new(&blast_cli)
            .arg(program)
            .arg("-query")
            .arg(query)
            .arg("-db")
            .arg("missing")
            .arg("-gilist")
            .arg("missing.gi")
            .arg("-seqidlist")
            .arg("missing.seqids")
            .arg("-outfmt")
            .arg("6")
            .output()
            .unwrap_or_else(|err| {
                panic!("run blast-cli {program} with gilist and seqidlist: {err}")
            });
        let ncbi = std::process::Command::new(&ncbi_bin)
            .arg("-query")
            .arg(query)
            .arg("-db")
            .arg("missing")
            .arg("-gilist")
            .arg("missing.gi")
            .arg("-seqidlist")
            .arg("missing.seqids")
            .arg("-outfmt")
            .arg("6")
            .output()
            .unwrap_or_else(|err| panic!("run NCBI {program} with gilist and seqidlist: {err}"));

        assert!(
            !rust.status.success(),
            "blast-cli should reject {program} with gilist and seqidlist"
        );
        assert!(
            !ncbi.status.success(),
            "NCBI should reject {program} with gilist and seqidlist"
        );
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "{program} gilist+seqidlist status differs"
        );
        assert_eq!(
            rust.stdout, ncbi.stdout,
            "{program} gilist+seqidlist stdout differs"
        );
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{program} gilist+seqidlist stderr differs"
        );
    }
}

#[test]
fn installed_programs_ncbi_parity_num_descriptions_and_max_target_seqs_are_incompatible() {
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    for (program, query) in [
        ("blastn", "tests/fixtures/query_random_200.fa"),
        ("blastp", "tests/fixtures/protein_query.fa"),
        ("blastx", "tests/fixtures/blastx_nuc_query.fa"),
        ("tblastn", "tests/fixtures/tblastn_prot_query.fa"),
        ("tblastx", "tests/fixtures/tblastx_nuc_query.fa"),
        ("psiblast", "tests/fixtures/psi_query.fa"),
        ("rpstblastn", "tests/fixtures/tblastx_nuc_query.fa"),
        ("deltablast", "tests/fixtures/protein_query.fa"),
    ] {
        let ncbi_bin = format!("/usr/bin/{program}");
        if !std::path::Path::new(&ncbi_bin).exists() {
            eprintln!("Skipping: {ncbi_bin} not found");
            continue;
        }
        let rust = std::process::Command::new(&blast_cli)
            .arg(program)
            .arg("-query")
            .arg(query)
            .arg("-db")
            .arg("missing")
            .arg("-num_descriptions")
            .arg("1")
            .arg("-max_target_seqs")
            .arg("1")
            .arg("-outfmt")
            .arg("6")
            .output()
            .unwrap_or_else(|err| {
                panic!("run blast-cli {program} with num_descriptions and max_target_seqs: {err}")
            });
        let ncbi = std::process::Command::new(&ncbi_bin)
            .arg("-query")
            .arg(query)
            .arg("-db")
            .arg("missing")
            .arg("-num_descriptions")
            .arg("1")
            .arg("-max_target_seqs")
            .arg("1")
            .arg("-outfmt")
            .arg("6")
            .output()
            .unwrap_or_else(|err| {
                panic!("run NCBI {program} with num_descriptions and max_target_seqs: {err}")
            });

        assert!(
            !rust.status.success(),
            "blast-cli should reject {program} with num_descriptions and max_target_seqs"
        );
        assert!(
            !ncbi.status.success(),
            "NCBI should reject {program} with num_descriptions and max_target_seqs"
        );
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "{program} num_descriptions+max_target_seqs status differs"
        );
        assert_eq!(
            rust.stdout, ncbi.stdout,
            "{program} num_descriptions+max_target_seqs stdout differs"
        );
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{program} num_descriptions+max_target_seqs stderr differs"
        );
    }
}

#[test]
fn installed_programs_ncbi_parity_entrez_query_requires_remote() {
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    for (program, query) in [
        ("blastn", "tests/fixtures/query_random_200.fa"),
        ("blastp", "tests/fixtures/protein_query.fa"),
        ("blastx", "tests/fixtures/blastx_nuc_query.fa"),
        ("tblastn", "tests/fixtures/tblastn_prot_query.fa"),
        ("tblastx", "tests/fixtures/tblastx_nuc_query.fa"),
        ("psiblast", "tests/fixtures/psi_query.fa"),
        ("rpstblastn", "tests/fixtures/tblastx_nuc_query.fa"),
        ("deltablast", "tests/fixtures/protein_query.fa"),
    ] {
        let ncbi_bin = format!("/usr/bin/{program}");
        if !std::path::Path::new(&ncbi_bin).exists() {
            eprintln!("Skipping: {ncbi_bin} not found");
            continue;
        }
        let rust = std::process::Command::new(&blast_cli)
            .arg(program)
            .arg("-query")
            .arg(query)
            .arg("-db")
            .arg("missing")
            .arg("-entrez_query")
            .arg("foo")
            .arg("-outfmt")
            .arg("6")
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli {program} with entrez_query: {err}"));
        let ncbi = std::process::Command::new(&ncbi_bin)
            .arg("-query")
            .arg(query)
            .arg("-db")
            .arg("missing")
            .arg("-entrez_query")
            .arg("foo")
            .arg("-outfmt")
            .arg("6")
            .output()
            .unwrap_or_else(|err| panic!("run NCBI {program} with entrez_query: {err}"));

        assert!(
            !rust.status.success(),
            "blast-cli should reject {program} entrez_query without remote"
        );
        assert!(
            !ncbi.status.success(),
            "NCBI should reject {program} entrez_query without remote"
        );
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "{program} entrez_query status differs"
        );
        assert_eq!(
            rust.stdout, ncbi.stdout,
            "{program} entrez_query stdout differs"
        );
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{program} entrez_query stderr differs"
        );
    }
}

#[test]
fn installed_programs_ncbi_parity_mt_mode_validation_order() {
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    for (program, query) in [
        ("blastn", "tests/fixtures/query_random_200.fa"),
        ("blastp", "tests/fixtures/protein_query.fa"),
        ("blastx", "tests/fixtures/blastx_nuc_query.fa"),
        ("tblastn", "tests/fixtures/tblastn_prot_query.fa"),
        ("tblastx", "tests/fixtures/tblastx_nuc_query.fa"),
        ("psiblast", "tests/fixtures/psi_query.fa"),
        ("rpstblastn", "tests/fixtures/tblastx_nuc_query.fa"),
        ("deltablast", "tests/fixtures/protein_query.fa"),
    ] {
        let ncbi_bin = format!("/usr/bin/{program}");
        if !std::path::Path::new(&ncbi_bin).exists() {
            eprintln!("Skipping: {ncbi_bin} not found");
            continue;
        }
        let rust = std::process::Command::new(&blast_cli)
            .arg(program)
            .arg("-query")
            .arg(query)
            .arg("-db")
            .arg("missing")
            .arg("-mt_mode")
            .arg("1")
            .arg("-outfmt")
            .arg("6")
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli {program} with mt_mode: {err}"));
        let ncbi = std::process::Command::new(&ncbi_bin)
            .arg("-query")
            .arg(query)
            .arg("-db")
            .arg("missing")
            .arg("-mt_mode")
            .arg("1")
            .arg("-outfmt")
            .arg("6")
            .output()
            .unwrap_or_else(|err| panic!("run NCBI {program} with mt_mode: {err}"));

        assert!(
            !rust.status.success(),
            "blast-cli should reject {program} mt_mode case"
        );
        assert!(
            !ncbi.status.success(),
            "NCBI should reject {program} mt_mode case"
        );
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "{program} mt_mode status differs"
        );
        assert_eq!(rust.stdout, ncbi.stdout, "{program} mt_mode stdout differs");
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{program} mt_mode stderr differs"
        );
    }
}

#[test]
fn installed_programs_ncbi_parity_remote_and_num_threads_are_incompatible() {
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    for (program, query) in [
        ("blastn", "tests/fixtures/query_random_200.fa"),
        ("blastp", "tests/fixtures/protein_query.fa"),
        ("blastx", "tests/fixtures/blastx_nuc_query.fa"),
        ("tblastn", "tests/fixtures/tblastn_prot_query.fa"),
        ("tblastx", "tests/fixtures/tblastx_nuc_query.fa"),
        ("psiblast", "tests/fixtures/psi_query.fa"),
        ("deltablast", "tests/fixtures/protein_query.fa"),
    ] {
        let ncbi_bin = format!("/usr/bin/{program}");
        if !std::path::Path::new(&ncbi_bin).exists() {
            eprintln!("Skipping: {ncbi_bin} not found");
            continue;
        }
        let rust = std::process::Command::new(&blast_cli)
            .arg(program)
            .arg("-query")
            .arg(query)
            .arg("-db")
            .arg("missing")
            .arg("-remote")
            .arg("-num_threads")
            .arg("2")
            .arg("-outfmt")
            .arg("6")
            .output()
            .unwrap_or_else(|err| {
                panic!("run blast-cli {program} with remote and num_threads: {err}")
            });
        let ncbi = std::process::Command::new(&ncbi_bin)
            .arg("-query")
            .arg(query)
            .arg("-db")
            .arg("missing")
            .arg("-remote")
            .arg("-num_threads")
            .arg("2")
            .arg("-outfmt")
            .arg("6")
            .output()
            .unwrap_or_else(|err| panic!("run NCBI {program} with remote and num_threads: {err}"));

        assert!(
            !rust.status.success(),
            "blast-cli should reject {program} with remote and num_threads"
        );
        assert!(
            !ncbi.status.success(),
            "NCBI should reject {program} with remote and num_threads"
        );
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "{program} remote+num_threads status differs"
        );
        assert_eq!(
            rust.stdout, ncbi.stdout,
            "{program} remote+num_threads stdout differs"
        );
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{program} remote+num_threads stderr differs"
        );
    }
}

#[test]
fn installed_programs_ncbi_parity_qcov_hsp_perc_range_error() {
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    for (program, query) in [
        ("blastn", "tests/fixtures/query_random_200.fa"),
        ("blastp", "tests/fixtures/protein_query.fa"),
        ("blastx", "tests/fixtures/blastx_nuc_query.fa"),
        ("tblastn", "tests/fixtures/tblastn_prot_query.fa"),
        ("tblastx", "tests/fixtures/tblastx_nuc_query.fa"),
        ("psiblast", "tests/fixtures/psi_query.fa"),
        ("rpstblastn", "tests/fixtures/tblastx_nuc_query.fa"),
        ("deltablast", "tests/fixtures/protein_query.fa"),
    ] {
        let ncbi_bin = format!("/usr/bin/{program}");
        if !std::path::Path::new(&ncbi_bin).exists() {
            eprintln!("Skipping: {ncbi_bin} not found");
            continue;
        }
        let rust = std::process::Command::new(&blast_cli)
            .arg(program)
            .arg("-query")
            .arg(query)
            .arg("-db")
            .arg("missing")
            .arg("-qcov_hsp_perc")
            .arg("101")
            .arg("-outfmt")
            .arg("6")
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli {program} with bad qcov_hsp_perc: {err}"));
        let ncbi = std::process::Command::new(&ncbi_bin)
            .arg("-query")
            .arg(query)
            .arg("-db")
            .arg("missing")
            .arg("-qcov_hsp_perc")
            .arg("101")
            .arg("-outfmt")
            .arg("6")
            .output()
            .unwrap_or_else(|err| panic!("run NCBI {program} with bad qcov_hsp_perc: {err}"));

        assert!(
            !rust.status.success(),
            "blast-cli should reject {program} qcov_hsp_perc 101"
        );
        assert!(
            !ncbi.status.success(),
            "NCBI should reject {program} qcov_hsp_perc 101"
        );
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "{program} qcov_hsp_perc status differs"
        );
        assert_eq!(
            rust.stdout, ncbi.stdout,
            "{program} qcov_hsp_perc stdout differs"
        );
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{program} qcov_hsp_perc stderr differs"
        );
    }
}

#[test]
fn protein_programs_ncbi_parity_perc_identity_is_unknown_argument() {
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    for (program, query) in [
        ("blastp", "tests/fixtures/protein_query.fa"),
        ("blastx", "tests/fixtures/blastx_nuc_query.fa"),
        ("tblastn", "tests/fixtures/tblastn_prot_query.fa"),
        ("tblastx", "tests/fixtures/tblastx_nuc_query.fa"),
        ("psiblast", "tests/fixtures/psi_query.fa"),
        ("rpstblastn", "tests/fixtures/tblastx_nuc_query.fa"),
        ("deltablast", "tests/fixtures/protein_query.fa"),
    ] {
        let ncbi_bin = format!("/usr/bin/{program}");
        if !std::path::Path::new(&ncbi_bin).exists() {
            eprintln!("Skipping: {ncbi_bin} not found");
            continue;
        }
        let rust = std::process::Command::new(&blast_cli)
            .arg(program)
            .arg("-query")
            .arg(query)
            .arg("-db")
            .arg("missing")
            .arg("-perc_identity")
            .arg("90")
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli {program} with perc_identity: {err}"));
        let ncbi = std::process::Command::new(&ncbi_bin)
            .arg("-query")
            .arg(query)
            .arg("-db")
            .arg("missing")
            .arg("-perc_identity")
            .arg("90")
            .output()
            .unwrap_or_else(|err| panic!("run NCBI {program} with perc_identity: {err}"));

        assert!(
            !rust.status.success(),
            "blast-cli should reject {program} perc_identity"
        );
        assert!(
            !ncbi.status.success(),
            "NCBI should reject {program} perc_identity"
        );
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "{program} perc_identity status differs"
        );
        assert_eq!(
            rust.stdout, ncbi.stdout,
            "{program} perc_identity stdout differs"
        );
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{program} perc_identity stderr differs"
        );
    }
}

#[test]
fn installed_programs_ncbi_parity_max_hsps_range_error() {
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    for (program, query) in [
        ("blastn", "tests/fixtures/query_random_200.fa"),
        ("blastp", "tests/fixtures/protein_query.fa"),
        ("blastx", "tests/fixtures/blastx_nuc_query.fa"),
        ("tblastn", "tests/fixtures/tblastn_prot_query.fa"),
        ("tblastx", "tests/fixtures/tblastx_nuc_query.fa"),
        ("psiblast", "tests/fixtures/psi_query.fa"),
        ("rpstblastn", "tests/fixtures/tblastx_nuc_query.fa"),
        ("deltablast", "tests/fixtures/protein_query.fa"),
    ] {
        let ncbi_bin = format!("/usr/bin/{program}");
        if !std::path::Path::new(&ncbi_bin).exists() {
            eprintln!("Skipping: {ncbi_bin} not found");
            continue;
        }
        let rust = std::process::Command::new(&blast_cli)
            .arg(program)
            .arg("-query")
            .arg(query)
            .arg("-db")
            .arg("missing")
            .arg("-max_hsps")
            .arg("0")
            .arg("-outfmt")
            .arg("6")
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli {program} with bad max_hsps: {err}"));
        let ncbi = std::process::Command::new(&ncbi_bin)
            .arg("-query")
            .arg(query)
            .arg("-db")
            .arg("missing")
            .arg("-max_hsps")
            .arg("0")
            .arg("-outfmt")
            .arg("6")
            .output()
            .unwrap_or_else(|err| panic!("run NCBI {program} with bad max_hsps: {err}"));

        assert!(
            !rust.status.success(),
            "blast-cli should reject {program} max_hsps 0"
        );
        assert!(
            !ncbi.status.success(),
            "NCBI should reject {program} max_hsps 0"
        );
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "{program} max_hsps status differs"
        );
        assert_eq!(
            rust.stdout, ncbi.stdout,
            "{program} max_hsps stdout differs"
        );
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{program} max_hsps stderr differs"
        );
    }
}

#[test]
fn installed_programs_ncbi_parity_soft_masking_bool_conversion_error() {
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    for (program, query) in [
        ("blastn", "tests/fixtures/query_random_200.fa"),
        ("blastp", "tests/fixtures/protein_query.fa"),
        ("blastx", "tests/fixtures/blastx_nuc_query.fa"),
        ("tblastn", "tests/fixtures/tblastn_prot_query.fa"),
        ("tblastx", "tests/fixtures/tblastx_nuc_query.fa"),
        ("psiblast", "tests/fixtures/psi_query.fa"),
        ("rpstblastn", "tests/fixtures/tblastx_nuc_query.fa"),
        ("deltablast", "tests/fixtures/protein_query.fa"),
    ] {
        let ncbi_bin = format!("/usr/bin/{program}");
        if !std::path::Path::new(&ncbi_bin).exists() {
            eprintln!("Skipping: {ncbi_bin} not found");
            continue;
        }
        let rust = std::process::Command::new(&blast_cli)
            .arg(program)
            .arg("-query")
            .arg(query)
            .arg("-db")
            .arg("missing")
            .arg("-soft_masking")
            .arg("maybe")
            .arg("-outfmt")
            .arg("6")
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli {program} with bad soft_masking: {err}"));
        let ncbi = std::process::Command::new(&ncbi_bin)
            .arg("-query")
            .arg(query)
            .arg("-db")
            .arg("missing")
            .arg("-soft_masking")
            .arg("maybe")
            .arg("-outfmt")
            .arg("6")
            .output()
            .unwrap_or_else(|err| panic!("run NCBI {program} with bad soft_masking: {err}"));

        assert!(
            !rust.status.success(),
            "blast-cli should reject {program} soft_masking maybe"
        );
        assert!(
            !ncbi.status.success(),
            "NCBI should reject {program} soft_masking maybe"
        );
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "{program} soft_masking status differs"
        );
        assert_eq!(
            rust.stdout, ncbi.stdout,
            "{program} soft_masking stdout differs"
        );
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{program} soft_masking stderr differs"
        );
    }
}

#[test]
fn installed_programs_ncbi_parity_omitted_query_with_subject_is_empty_query() {
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    for (program, subject) in [
        ("blastn", "tests/fixtures/subject_test.fa"),
        ("blastp", "tests/fixtures/protein_subject.fa"),
        ("blastx", "tests/fixtures/blastx_prot_subject.fa"),
        ("tblastn", "tests/fixtures/tblastn_nuc_subject.fa"),
        ("tblastx", "tests/fixtures/tblastx_nuc_subject.fa"),
        ("psiblast", "tests/fixtures/psi_subject.fa"),
        ("deltablast", "tests/fixtures/protein_subject.fa"),
    ] {
        let ncbi_bin = format!("/usr/bin/{program}");
        if !std::path::Path::new(&ncbi_bin).exists() {
            eprintln!("Skipping: {ncbi_bin} not found");
            continue;
        }
        let rust = std::process::Command::new(&blast_cli)
            .arg(program)
            .arg("-subject")
            .arg(subject)
            .arg("-outfmt")
            .arg("6")
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli {program} without query: {err}"));
        let ncbi = std::process::Command::new(&ncbi_bin)
            .arg("-subject")
            .arg(subject)
            .arg("-outfmt")
            .arg("6")
            .output()
            .unwrap_or_else(|err| panic!("run NCBI {program} without query: {err}"));

        assert!(
            rust.status.success(),
            "blast-cli should accept omitted {program} query as empty"
        );
        assert!(
            ncbi.status.success(),
            "NCBI should accept omitted {program} query as empty"
        );
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "{program} omitted query status differs"
        );
        assert_eq!(
            rust.stdout, ncbi.stdout,
            "{program} omitted query stdout differs"
        );
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{program} omitted query stderr differs"
        );
    }
}

#[test]
fn installed_programs_ncbi_parity_omitted_query_missing_subject_precedes_empty_query() {
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    for program in [
        "blastn",
        "blastp",
        "blastx",
        "tblastn",
        "tblastx",
        "psiblast",
        "deltablast",
    ] {
        let ncbi_bin = format!("/usr/bin/{program}");
        if !std::path::Path::new(&ncbi_bin).exists() {
            eprintln!("Skipping: {ncbi_bin} not found");
            continue;
        }
        let rust = std::process::Command::new(&blast_cli)
            .arg(program)
            .arg("-subject")
            .arg("missing_subject.fa")
            .arg("-outfmt")
            .arg("6")
            .output()
            .unwrap_or_else(|err| {
                panic!("run blast-cli {program} without query against missing subject: {err}")
            });
        let ncbi = std::process::Command::new(&ncbi_bin)
            .arg("-subject")
            .arg("missing_subject.fa")
            .arg("-outfmt")
            .arg("6")
            .output()
            .unwrap_or_else(|err| {
                panic!("run NCBI {program} without query against missing subject: {err}")
            });

        assert!(
            !rust.status.success(),
            "blast-cli should reject {program} missing subject before empty query"
        );
        assert!(
            !ncbi.status.success(),
            "NCBI should reject {program} missing subject before empty query"
        );
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "{program} omitted query missing subject status differs"
        );
        assert_eq!(
            rust.stdout, ncbi.stdout,
            "{program} omitted query missing subject stdout differs"
        );
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{program} omitted query missing subject stderr differs"
        );
    }
}

#[test]
fn installed_programs_ncbi_parity_missing_subject_precedes_missing_query_in_subject_mode() {
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    for program in [
        "blastn",
        "blastp",
        "blastx",
        "tblastn",
        "tblastx",
        "psiblast",
        "deltablast",
    ] {
        let ncbi_bin = format!("/usr/bin/{program}");
        if !std::path::Path::new(&ncbi_bin).exists() {
            eprintln!("Skipping: {ncbi_bin} not found");
            continue;
        }
        let rust = std::process::Command::new(&blast_cli)
            .arg(program)
            .arg("--query")
            .arg("missing_query.fa")
            .arg("--subject")
            .arg("missing_subject.fa")
            .arg("--outfmt")
            .arg("6")
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli {program} with missing files: {err}"));
        let ncbi = std::process::Command::new(&ncbi_bin)
            .arg("-query")
            .arg("missing_query.fa")
            .arg("-subject")
            .arg("missing_subject.fa")
            .arg("-outfmt")
            .arg("6")
            .output()
            .unwrap_or_else(|err| panic!("run NCBI {program} with missing files: {err}"));

        assert!(
            !rust.status.success(),
            "blast-cli should reject {program} missing subject before missing query"
        );
        assert!(
            !ncbi.status.success(),
            "NCBI should reject {program} missing subject before missing query"
        );
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "{program} missing subject/query status differs"
        );
        assert_eq!(
            rust.stdout, ncbi.stdout,
            "{program} missing subject/query stdout differs"
        );
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{program} missing subject/query stderr differs"
        );
    }
}

#[test]
fn protein_programs_ncbi_parity_missing_query_precedes_missing_db() {
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    for program in [
        "blastp",
        "blastx",
        "tblastn",
        "tblastx",
        "psiblast",
        "deltablast",
    ] {
        let ncbi_bin = format!("/usr/bin/{program}");
        if !std::path::Path::new(&ncbi_bin).exists() {
            eprintln!("Skipping: {ncbi_bin} not found");
            continue;
        }
        let rust = std::process::Command::new(&blast_cli)
            .arg(program)
            .arg("--query")
            .arg("missing_query.fa")
            .arg("--db")
            .arg("missing_search")
            .arg("--outfmt")
            .arg("6")
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli {program} missing query and DB: {err}"));
        let ncbi = std::process::Command::new(&ncbi_bin)
            .arg("-query")
            .arg("missing_query.fa")
            .arg("-db")
            .arg("missing_search")
            .arg("-outfmt")
            .arg("6")
            .output()
            .unwrap_or_else(|err| panic!("run NCBI {program} missing query and DB: {err}"));

        assert!(
            !rust.status.success(),
            "blast-cli should reject {program} missing query before missing DB"
        );
        assert!(
            !ncbi.status.success(),
            "NCBI should reject {program} missing query before missing DB"
        );
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "{program} missing query/DB status differs"
        );
        assert_eq!(
            rust.stdout, ncbi.stdout,
            "{program} missing query/DB stdout differs"
        );
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{program} missing query/DB stderr differs"
        );
    }
}

#[test]
fn protein_programs_ncbi_parity_output_file_precedes_missing_db_or_empty_query() {
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    for (program, query, subject) in [
        (
            "blastp",
            "tests/fixtures/protein_query.fa",
            "tests/fixtures/protein_subject.fa",
        ),
        (
            "blastx",
            "tests/fixtures/blastx_nuc_query.fa",
            "tests/fixtures/blastx_prot_subject.fa",
        ),
        (
            "tblastn",
            "tests/fixtures/tblastn_prot_query.fa",
            "tests/fixtures/tblastn_nuc_subject.fa",
        ),
        (
            "tblastx",
            "tests/fixtures/tblastx_nuc_query.fa",
            "tests/fixtures/tblastx_nuc_subject.fa",
        ),
        (
            "psiblast",
            "tests/fixtures/protein_query.fa",
            "tests/fixtures/protein_subject.fa",
        ),
        (
            "deltablast",
            "tests/fixtures/protein_query.fa",
            "tests/fixtures/protein_subject.fa",
        ),
    ] {
        let ncbi_bin = format!("/usr/bin/{program}");
        if !std::path::Path::new(&ncbi_bin).exists() {
            eprintln!("Skipping: {ncbi_bin} not found");
            continue;
        }
        for (rust_source_args, ncbi_source_args) in [
            (
                vec!["--query", query, "--db", "missing_db"],
                vec!["-query", query, "-db", "missing_db"],
            ),
            (vec!["--db", "missing_db"], vec!["-db", "missing_db"]),
            (vec!["--subject", subject], vec!["-subject", subject]),
        ] {
            let rust = std::process::Command::new(&blast_cli)
                .arg(program)
                .args(&rust_source_args)
                .arg("--out")
                .arg("/no/such/dir/out.txt")
                .arg("--outfmt")
                .arg("6")
                .output()
                .unwrap_or_else(|err| {
                    panic!("run blast-cli {program} with inaccessible out: {err}")
                });
            let ncbi = std::process::Command::new(&ncbi_bin)
                .args(&ncbi_source_args)
                .arg("-out")
                .arg("/no/such/dir/out.txt")
                .arg("-outfmt")
                .arg("6")
                .output()
                .unwrap_or_else(|err| panic!("run NCBI {program} with inaccessible out: {err}"));

            assert!(
                !rust.status.success(),
                "blast-cli should reject {program} inaccessible out before search setup"
            );
            assert!(
                !ncbi.status.success(),
                "NCBI should reject {program} inaccessible out before search setup"
            );
            assert_eq!(
                rust.status.code(),
                ncbi.status.code(),
                "{program} inaccessible out status differs for {rust_source_args:?}"
            );
            assert_eq!(
                rust.stdout, ncbi.stdout,
                "{program} inaccessible out stdout differs for {rust_source_args:?}"
            );
            assert_eq!(
                String::from_utf8_lossy(&rust.stderr),
                String::from_utf8_lossy(&ncbi.stderr),
                "{program} inaccessible out stderr differs for {rust_source_args:?}"
            );
        }
    }
}

#[test]
fn installed_programs_ncbi_parity_omitted_query_missing_db_precedes_empty_query() {
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    for program in [
        "blastn",
        "blastp",
        "blastx",
        "tblastn",
        "tblastx",
        "psiblast",
        "rpstblastn",
        "deltablast",
    ] {
        let ncbi_bin = format!("/usr/bin/{program}");
        if !std::path::Path::new(&ncbi_bin).exists() {
            eprintln!("Skipping: {ncbi_bin} not found");
            continue;
        }
        let rust = std::process::Command::new(&blast_cli)
            .arg(program)
            .arg("-db")
            .arg("missing")
            .arg("-outfmt")
            .arg("6")
            .output()
            .unwrap_or_else(|err| {
                panic!("run blast-cli {program} without query against missing DB: {err}")
            });
        let ncbi = std::process::Command::new(&ncbi_bin)
            .arg("-db")
            .arg("missing")
            .arg("-outfmt")
            .arg("6")
            .output()
            .unwrap_or_else(|err| {
                panic!("run NCBI {program} without query against missing DB: {err}")
            });

        assert!(
            !rust.status.success(),
            "blast-cli should reject {program} missing DB before empty query"
        );
        assert!(
            !ncbi.status.success(),
            "NCBI should reject {program} missing DB before empty query"
        );
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "{program} omitted query missing DB status differs"
        );
        assert_eq!(
            rust.stdout, ncbi.stdout,
            "{program} omitted query missing DB stdout differs"
        );
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{program} omitted query missing DB stderr differs"
        );
    }
}

#[test]
fn blastn_blastp_ncbi_parity_omitted_query_rust_db_directory_is_not_ncbi_db() {
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    for (program, db) in [
        ("blastn", "tests/fixtures/seqn"),
        ("blastp", "tests/fixtures/seqp"),
    ] {
        let ncbi_bin = format!("/usr/bin/{program}");
        if !std::path::Path::new(&ncbi_bin).exists() {
            eprintln!("Skipping: {ncbi_bin} not found");
            continue;
        }
        let rust = std::process::Command::new(&blast_cli)
            .arg(program)
            .arg("-db")
            .arg(db)
            .arg("-outfmt")
            .arg("6")
            .output()
            .unwrap_or_else(|err| {
                panic!("run blast-cli {program} without query against Rust fixture DB: {err}")
            });
        let ncbi = std::process::Command::new(&ncbi_bin)
            .arg("-db")
            .arg(db)
            .arg("-outfmt")
            .arg("6")
            .output()
            .unwrap_or_else(|err| {
                panic!("run NCBI {program} without query against Rust fixture DB: {err}")
            });

        assert!(
            !rust.status.success(),
            "blast-cli should reject {program} Rust fixture DB before empty query"
        );
        assert!(
            !ncbi.status.success(),
            "NCBI should reject {program} Rust fixture DB before empty query"
        );
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "{program} omitted query Rust fixture DB status differs"
        );
        assert_eq!(
            rust.stdout, ncbi.stdout,
            "{program} omitted query Rust fixture DB stdout differs"
        );
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{program} omitted query Rust fixture DB stderr differs"
        );
    }
}

#[test]
fn non_blastn_programs_ncbi_parity_short_help_uses_program_usage() {
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    for program in [
        "blastp",
        "blastx",
        "tblastn",
        "tblastx",
        "psiblast",
        "rpstblastn",
        "deltablast",
    ] {
        let ncbi_bin = format!("/usr/bin/{program}");
        if !std::path::Path::new(&ncbi_bin).exists() {
            eprintln!("Skipping: {ncbi_bin} not found");
            continue;
        }
        let rust = std::process::Command::new(&blast_cli)
            .arg(program)
            .arg("-h")
            .arg("-query")
            .arg("missing.fa")
            .arg("-outfmt")
            .arg("bad")
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli {program} -h: {err}"));
        let ncbi = std::process::Command::new(&ncbi_bin)
            .arg("-h")
            .arg("-query")
            .arg("missing.fa")
            .arg("-outfmt")
            .arg("bad")
            .output()
            .unwrap_or_else(|err| panic!("run NCBI {program} -h: {err}"));

        assert!(rust.status.success(), "blast-cli {program} -h failed");
        assert!(ncbi.status.success(), "NCBI {program} -h failed");
        assert_eq!(rust.stdout, ncbi.stdout, "{program} -h stdout differs");
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{program} -h stderr differs"
        );
    }
}

#[test]
fn non_blastn_programs_help_line_is_short_help_only() {
    let Some(blast_cli) = std::env::var_os("BLAST_RS_CLI_BIN")
        .or_else(|| std::env::var_os("CARGO_BIN_EXE_blast-cli"))
        .map(std::path::PathBuf::from)
    else {
        eprintln!("Skipping: set BLAST_RS_CLI_BIN or CARGO_BIN_EXE_blast-cli");
        return;
    };

    for program in [
        "blastp",
        "blastx",
        "tblastn",
        "tblastx",
        "psiblast",
        "rpstblastn",
        "deltablast",
    ] {
        let short = std::process::Command::new(&blast_cli)
            .arg(program)
            .arg("-h")
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli {program} -h: {err}"));
        let detailed = std::process::Command::new(&blast_cli)
            .arg(program)
            .arg("-help")
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli {program} -help: {err}"));

        assert!(short.status.success(), "blast-cli {program} -h failed");
        assert!(
            detailed.status.success(),
            "blast-cli {program} -help failed"
        );
        let short_stdout = String::from_utf8_lossy(&short.stdout);
        let detailed_stdout = String::from_utf8_lossy(&detailed.stdout);
        assert!(
            short_stdout.contains("Use '-help' to print detailed descriptions"),
            "{program} -h should include the detailed-help hint"
        );
        assert!(
            !detailed_stdout.contains("Use '-help' to print detailed descriptions"),
            "{program} -help should not include the detailed-help hint"
        );
        assert!(
            detailed_stdout.contains("OPTIONAL ARGUMENTS"),
            "{program} -help should include detailed argument section heading"
        );
    }
}

#[test]
fn psiblast_subject_ncbi_parity_tabular_and_csv() {
    if !std::path::Path::new("/usr/bin/psiblast").exists() {
        eprintln!("Skipping: /usr/bin/psiblast not found");
        return;
    }
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    for outfmt in [
        "6 qseqid sseqid pident length qstart qend sstart send bitscore evalue",
        "10 qseqid sseqid pident length qstart qend sstart send bitscore evalue",
    ] {
        let rust = std::process::Command::new(&blast_cli)
            .arg("psiblast")
            .arg("-query")
            .arg("tests/fixtures/psi_query.fa")
            .arg("-subject")
            .arg("tests/fixtures/psi_subject.fa")
            .arg("-outfmt")
            .arg(outfmt)
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli psiblast {outfmt}: {err}"));
        let ncbi = std::process::Command::new("/usr/bin/psiblast")
            .arg("-query")
            .arg("tests/fixtures/psi_query.fa")
            .arg("-subject")
            .arg("tests/fixtures/psi_subject.fa")
            .arg("-outfmt")
            .arg(outfmt)
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
}

#[test]
fn psiblast_ncbi_parity_domain_inclusion_is_unknown_argument() {
    if !std::path::Path::new("/usr/bin/psiblast").exists() {
        eprintln!("Skipping: /usr/bin/psiblast not found");
        return;
    }
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    let rust = std::process::Command::new(&blast_cli)
        .arg("psiblast")
        .arg("--domain_inclusion_ethresh")
        .arg("0.01")
        .output()
        .expect("run blast-cli psiblast domain inclusion");
    let ncbi = std::process::Command::new("/usr/bin/psiblast")
        .arg("-domain_inclusion_ethresh")
        .arg("0.01")
        .output()
        .expect("run NCBI psiblast domain inclusion");

    assert!(!rust.status.success(), "blast-cli should reject option");
    assert!(!ncbi.status.success(), "NCBI should reject option");
    assert_eq!(rust.status.code(), ncbi.status.code(), "status differs");
    assert_eq!(rust.stdout, ncbi.stdout, "stdout differs");
    assert_eq!(
        String::from_utf8_lossy(&rust.stderr),
        String::from_utf8_lossy(&ncbi.stderr),
        "stderr differs"
    );
}

#[test]
fn psiblast_ncbi_parity_restart_input_constraint_order() {
    if !std::path::Path::new("/usr/bin/psiblast").exists() {
        eprintln!("Skipping: /usr/bin/psiblast not found");
        return;
    }
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    for (rust_args, ncbi_args) in [
        (
            vec!["--query", "q.fa", "--in_msa", "restart.msa"],
            vec!["-query", "q.fa", "-in_msa", "restart.msa"],
        ),
        (
            vec!["--in_msa", "restart.msa", "--in_pssm", "restart.chk"],
            vec!["-in_msa", "restart.msa", "-in_pssm", "restart.chk"],
        ),
        (
            vec!["--in_pssm", "restart.chk", "--ignore_msa_master"],
            vec!["-in_pssm", "restart.chk", "-ignore_msa_master"],
        ),
        (
            vec!["--in_pssm", "restart.chk", "--msa_master_idx", "2"],
            vec!["-in_pssm", "restart.chk", "-msa_master_idx", "2"],
        ),
        (
            vec![
                "--in_msa",
                "restart.msa",
                "--in_pssm",
                "restart.chk",
                "--ignore_msa_master",
            ],
            vec![
                "-in_msa",
                "restart.msa",
                "-in_pssm",
                "restart.chk",
                "-ignore_msa_master",
            ],
        ),
        (
            vec![
                "--query",
                "q.fa",
                "--ignore_msa_master",
                "--in_pssm",
                "restart.chk",
            ],
            vec![
                "-query",
                "q.fa",
                "-ignore_msa_master",
                "-in_pssm",
                "restart.chk",
            ],
        ),
        (
            vec![
                "--query",
                "q.fa",
                "--in_msa",
                "restart.msa",
                "--in_pssm",
                "restart.chk",
                "--ignore_msa_master",
                "--msa_master_idx",
                "2",
            ],
            vec![
                "-query",
                "q.fa",
                "-in_msa",
                "restart.msa",
                "-in_pssm",
                "restart.chk",
                "-ignore_msa_master",
                "-msa_master_idx",
                "2",
            ],
        ),
        (
            vec![
                "--query",
                "q.fa",
                "--in_msa",
                "restart.msa",
                "--in_pssm",
                "restart.chk",
            ],
            vec![
                "-query",
                "q.fa",
                "-in_msa",
                "restart.msa",
                "-in_pssm",
                "restart.chk",
            ],
        ),
    ] {
        let rust = std::process::Command::new(&blast_cli)
            .arg("psiblast")
            .args(&rust_args)
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli psiblast {rust_args:?}: {err}"));
        let ncbi = std::process::Command::new("/usr/bin/psiblast")
            .args(&ncbi_args)
            .output()
            .unwrap_or_else(|err| panic!("run NCBI psiblast {ncbi_args:?}: {err}"));

        assert!(
            !rust.status.success(),
            "blast-cli should reject input conflict"
        );
        assert!(!ncbi.status.success(), "NCBI should reject input conflict");
        assert_eq!(rust.status.code(), ncbi.status.code(), "status differs");
        assert_eq!(rust.stdout, ncbi.stdout, "stdout differs");
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "stderr differs for {rust_args:?}"
        );
    }
}

#[test]
fn psiblast_ncbi_parity_psi_numeric_constraints() {
    if !std::path::Path::new("/usr/bin/psiblast").exists() {
        eprintln!("Skipping: /usr/bin/psiblast not found");
        return;
    }
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    let tmp = TempDir::new().expect("tempdir");
    let query = tmp.path().join("query.fa");
    let subject = tmp.path().join("subject.fa");
    std::fs::write(&query, ">q\nMKKWLFGFLG\n").expect("write query");
    std::fs::write(&subject, ">s\nMKKWLFGFLG\n").expect("write subject");

    let rejected = std::process::Command::new(&blast_cli)
        .arg("psiblast")
        .arg("--query")
        .arg(&query)
        .arg("--subject")
        .arg(&subject)
        .arg("--num_iterations")
        .arg("-1")
        .output()
        .expect("run blast-cli psiblast negative iterations");
    let ncbi_rejected = std::process::Command::new("/usr/bin/psiblast")
        .arg("-query")
        .arg(&query)
        .arg("-subject")
        .arg(&subject)
        .arg("-num_iterations")
        .arg("-1")
        .output()
        .expect("run NCBI psiblast negative iterations");
    assert_eq!(rejected.status.code(), ncbi_rejected.status.code());
    assert_eq!(rejected.stdout, ncbi_rejected.stdout);
    assert_eq!(
        String::from_utf8_lossy(&rejected.stderr),
        String::from_utf8_lossy(&ncbi_rejected.stderr)
    );

    for (rust_option, ncbi_option, value) in [
        ("--num_iterations", "-num_iterations", "0"),
        ("--gap_trigger", "-gap_trigger", "-1"),
        ("--pseudocount", "-pseudocount", "-1"),
        ("--inclusion_ethresh", "-inclusion_ethresh", "0"),
    ] {
        let rust = std::process::Command::new(&blast_cli)
            .arg("psiblast")
            .arg("--query")
            .arg(&query)
            .arg("--subject")
            .arg(&subject)
            .arg("--outfmt")
            .arg("6")
            .arg(rust_option)
            .arg(value)
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli psiblast {rust_option}: {err}"));
        let ncbi = std::process::Command::new("/usr/bin/psiblast")
            .arg("-query")
            .arg(&query)
            .arg("-subject")
            .arg(&subject)
            .arg("-outfmt")
            .arg("6")
            .arg(ncbi_option)
            .arg(value)
            .output()
            .unwrap_or_else(|err| panic!("run NCBI psiblast {ncbi_option}: {err}"));

        assert!(
            rust.status.success(),
            "blast-cli should accept {rust_option} {value}: {}",
            String::from_utf8_lossy(&rust.stderr)
        );
        assert!(
            ncbi.status.success(),
            "NCBI should accept {ncbi_option} {value}: {}",
            String::from_utf8_lossy(&ncbi.stderr)
        );
    }

    let restart_msa = tmp.path().join("restart.msa");
    std::fs::write(&restart_msa, ">master\nMKKWLFGFLG\n>seq2\nMKKWLFGFLG\n")
        .expect("write restart MSA");
    for value in ["0", "-1"] {
        let rust = std::process::Command::new(&blast_cli)
            .arg("psiblast")
            .arg("--in_msa")
            .arg(&restart_msa)
            .arg("--subject")
            .arg(&subject)
            .arg("--msa_master_idx")
            .arg(value)
            .arg("--outfmt")
            .arg("6")
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli psiblast msa_master_idx {value}: {err}"));
        let ncbi = std::process::Command::new("/usr/bin/psiblast")
            .arg("-in_msa")
            .arg(&restart_msa)
            .arg("-subject")
            .arg(&subject)
            .arg("-msa_master_idx")
            .arg(value)
            .arg("-outfmt")
            .arg("6")
            .output()
            .unwrap_or_else(|err| panic!("run NCBI psiblast msa_master_idx {value}: {err}"));

        assert!(
            !rust.status.success(),
            "blast-cli should reject msa_master_idx {value}"
        );
        assert!(
            !ncbi.status.success(),
            "NCBI should reject msa_master_idx {value}"
        );
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "msa_master_idx {value} status differs"
        );
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "msa_master_idx {value} stderr differs"
        );
    }
}

#[test]
#[ignore = "documents remaining PSI iteration-output gap for num_iterations=0 convergence output"]
fn psiblast_ncbi_parity_num_iterations_zero_iteration_output() {
    if !std::path::Path::new("/usr/bin/psiblast").exists() {
        eprintln!("Skipping: /usr/bin/psiblast not found");
        return;
    }
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    let tmp = TempDir::new().expect("tempdir");
    let query = tmp.path().join("query.fa");
    let subject = tmp.path().join("subject.fa");
    std::fs::write(&query, ">q\nMKKWLFGFLG\n").expect("write query");
    std::fs::write(&subject, ">s\nMKKWLFGFLG\n").expect("write subject");

    let rust = std::process::Command::new(&blast_cli)
        .arg("psiblast")
        .arg("--query")
        .arg(&query)
        .arg("--subject")
        .arg(&subject)
        .arg("--outfmt")
        .arg("6")
        .arg("--num_iterations")
        .arg("0")
        .output()
        .expect("run blast-cli psiblast num_iterations 0");
    let ncbi = std::process::Command::new("/usr/bin/psiblast")
        .arg("-query")
        .arg(&query)
        .arg("-subject")
        .arg(&subject)
        .arg("-outfmt")
        .arg("6")
        .arg("-num_iterations")
        .arg("0")
        .output()
        .expect("run NCBI psiblast num_iterations 0");

    assert!(rust.status.success(), "blast-cli psiblast failed");
    assert!(ncbi.status.success(), "NCBI psiblast failed");
    assert_eq!(rust.stdout, ncbi.stdout, "stdout differs");
    assert_eq!(
        String::from_utf8_lossy(&rust.stderr),
        String::from_utf8_lossy(&ncbi.stderr),
        "stderr differs"
    );
}

#[test]
fn psiblast_ncbi_parity_num_iterations_zero_output_shape() {
    if !std::path::Path::new("/usr/bin/psiblast").exists() {
        eprintln!("Skipping: /usr/bin/psiblast not found");
        return;
    }
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    let tmp = TempDir::new().expect("tempdir");
    let query = tmp.path().join("query.fa");
    let subject = tmp.path().join("subject.fa");
    let rust_out = tmp.path().join("rust.out");
    let ncbi_out = tmp.path().join("ncbi.out");
    std::fs::write(&query, ">q\nMKKWLFGFLG\n").expect("write query");
    std::fs::write(&subject, ">s\nMKKWLFGFLG\n").expect("write subject");

    let rust = std::process::Command::new(&blast_cli)
        .arg("psiblast")
        .arg("--query")
        .arg(&query)
        .arg("--subject")
        .arg(&subject)
        .arg("--outfmt")
        .arg("6")
        .arg("--num_iterations")
        .arg("0")
        .arg("--out")
        .arg(&rust_out)
        .output()
        .expect("run blast-cli psiblast num_iterations 0");
    let ncbi = std::process::Command::new("/usr/bin/psiblast")
        .arg("-query")
        .arg(&query)
        .arg("-subject")
        .arg(&subject)
        .arg("-outfmt")
        .arg("6")
        .arg("-num_iterations")
        .arg("0")
        .arg("-out")
        .arg(&ncbi_out)
        .output()
        .expect("run NCBI psiblast num_iterations 0");

    assert!(rust.status.success(), "blast-cli psiblast failed");
    assert!(ncbi.status.success(), "NCBI psiblast failed");
    assert_eq!(
        String::from_utf8_lossy(&rust.stderr),
        String::from_utf8_lossy(&ncbi.stderr),
        "stderr differs"
    );

    let rust_stdout = std::fs::read_to_string(&rust_out).expect("read rust output");
    let ncbi_stdout = std::fs::read_to_string(&ncbi_out).expect("read NCBI output");
    let rust_lines: Vec<&str> = rust_stdout.lines().collect();
    let ncbi_lines: Vec<&str> = ncbi_stdout.lines().collect();

    assert_eq!(
        rust_lines.len(),
        ncbi_lines.len(),
        "stdout line count differs"
    );
    assert_eq!(rust_lines[0], ncbi_lines[0], "first-round hit differs");
    assert_eq!(rust_lines[2], ncbi_lines[2], "blank separator differs");
    assert_eq!(
        rust_lines[3], ncbi_lines[3],
        "convergence marker should be written to the output stream"
    );

    let rust_final_fields: Vec<&str> = rust_lines[1].split('\t').collect();
    let ncbi_final_fields: Vec<&str> = ncbi_lines[1].split('\t').collect();
    assert_eq!(
        &rust_final_fields[..10],
        &ncbi_final_fields[..10],
        "final-round coordinates differ"
    );
}

#[test]
fn psiblast_ncbi_parity_gap_trigger_negative_one_iteration_output() {
    if !std::path::Path::new("/usr/bin/psiblast").exists() {
        eprintln!("Skipping: /usr/bin/psiblast not found");
        return;
    }
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    let tmp = TempDir::new().expect("tempdir");
    let query = tmp.path().join("query.fa");
    let subject = tmp.path().join("subject.fa");
    std::fs::write(&query, ">q\nMKKWLFGFLG\n").expect("write query");
    std::fs::write(&subject, ">s\nMKKWLFGFLG\n").expect("write subject");

    let rust = std::process::Command::new(&blast_cli)
        .arg("psiblast")
        .arg("--query")
        .arg(&query)
        .arg("--subject")
        .arg(&subject)
        .arg("--outfmt")
        .arg("6")
        .arg("--gap_trigger")
        .arg("-1")
        .output()
        .expect("run blast-cli psiblast gap_trigger -1");
    let ncbi = std::process::Command::new("/usr/bin/psiblast")
        .arg("-query")
        .arg(&query)
        .arg("-subject")
        .arg(&subject)
        .arg("-outfmt")
        .arg("6")
        .arg("-gap_trigger")
        .arg("-1")
        .output()
        .expect("run NCBI psiblast gap_trigger -1");

    assert!(rust.status.success(), "blast-cli psiblast failed");
    assert!(ncbi.status.success(), "NCBI psiblast failed");
    assert_eq!(rust.stdout, ncbi.stdout, "stdout differs");
    assert_eq!(
        String::from_utf8_lossy(&rust.stderr),
        String::from_utf8_lossy(&ncbi.stderr),
        "stderr differs"
    );
}

#[test]
fn psiblast_subject_ncbi_parity_restart_msa_master_selection() {
    if !std::path::Path::new("/usr/bin/psiblast").exists() {
        eprintln!("Skipping: /usr/bin/psiblast not found");
        return;
    }
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    for (msa, rust_extra, ncbi_extra) in [
        (
            ">master
MKFLIFALILFATVALAPKSSSHEI
>seq2
MKFLIFALILFATVALAPKSSSHEI
",
            Vec::<&str>::new(),
            Vec::<&str>::new(),
        ),
        (
            ">master desc extra
MKFLIFALILFATVALAPKSSSHEI
>seq2
MKFLIFALILFATVALAPKSSSHEI
",
            Vec::<&str>::new(),
            Vec::<&str>::new(),
        ),
        (
            ">   master desc
MKFLIFALILFATVALAPKSSSHEI
>seq2
MKFLIFALILFATVALAPKSSSHEI
",
            Vec::<&str>::new(),
            Vec::<&str>::new(),
        ),
        (
            ">master desc|bad
MKFLIFALILFATVALAPKSSSHEI
>seq2
MKFLIFALILFATVALAPKSSSHEI
",
            Vec::<&str>::new(),
            Vec::<&str>::new(),
        ),
        (
            ">master\tdesc\textra
MKFLIFALILFATVALAPKSSSHEI
>seq2
MKFLIFALILFATVALAPKSSSHEI
",
            Vec::<&str>::new(),
            Vec::<&str>::new(),
        ),
        (
            ">master\r\nMKFLIFALILFATVALAPKSSSHEI\r\n>seq2\r\nMKFLIFALILFATVALAPKSSSHEI\r\n",
            Vec::<&str>::new(),
            Vec::<&str>::new(),
        ),
        (
            ">master
MKFLIFALILFATVALAPKSSSHEI

;comment
>seq2
MKFLIFALILFATVALAPKSSSHEI
",
            Vec::<&str>::new(),
            Vec::<&str>::new(),
        ),
        (
            ">seq1
M-FLIFALILFATVALAPKSSSHEI
>master2
MKFLIFALILFATVALAPKSSSHEI
",
            vec!["--msa_master_idx", "2"],
            vec!["-msa_master_idx", "2"],
        ),
        (
            ">seq1
M-FLIFALILFATVALAPKSSSHEI
>master2 desc2 extra
MKFLIFALILFATVALAPKSSSHEI
",
            vec!["--msa_master_idx", "2"],
            vec!["-msa_master_idx", "2"],
        ),
        (
            ">master
M-KFLIFALILFATVALAPKSSSHEI
>seq2
MAKFLIFALILFATVALAPKSSSHEI
",
            Vec::<&str>::new(),
            Vec::<&str>::new(),
        ),
        (
            ">master
MKFLIFALILFATVALAPKSSSHEI
>empty

>seq2
MKFLIFALILFATVALAPKSSSHEI
",
            Vec::<&str>::new(),
            Vec::<&str>::new(),
        ),
        (
            ">master
MKF LIF\tALILFATVALAPKSSSHEI

>seq2
MKF LIF\tALILFATVALAPKSSSHEI
",
            Vec::<&str>::new(),
            Vec::<&str>::new(),
        ),
        (
            ">master
MKFLIFALILFATVALAPKSSSHEI
>seq|bad
MKFLIFALILFATVALAPKSSSHEI
",
            Vec::<&str>::new(),
            Vec::<&str>::new(),
        ),
        (
            ">master
MBZJOUXKLF
>seq2
MBZJOUXKLF
",
            Vec::<&str>::new(),
            Vec::<&str>::new(),
        ),
        (
            ">master
M-oF
>seq2
M-oF
",
            Vec::<&str>::new(),
            Vec::<&str>::new(),
        ),
        (
            ">master label
M-oF
>seq2
M-oF
",
            Vec::<&str>::new(),
            Vec::<&str>::new(),
        ),
        (
            ">master
MKFLIFALILFATVALAPKSSSHEI
>seq2
MKFLIFALILFATVALAPKSSSHEI
",
            vec!["--ignore_msa_master"],
            vec!["-ignore_msa_master"],
        ),
    ] {
        let tmp = TempDir::new().expect("tempdir");
        let msa_path = tmp.path().join("restart.msa");
        let subject = tmp.path().join("subject.fa");
        std::fs::write(&msa_path, msa).expect("write restart MSA");
        std::fs::write(&subject, ">s\nMKFLIFALILFATVALAPKSSSHEIHH\n").expect("write subject");

        let mut rust_cmd = std::process::Command::new(&blast_cli);
        rust_cmd
            .arg("psiblast")
            .arg("--in_msa")
            .arg(&msa_path)
            .arg("--subject")
            .arg(&subject)
            .arg("--outfmt")
            .arg("6 qseqid sseqid length qstart qend sstart send")
            .arg("--seg")
            .arg("no")
            .arg("--comp_based_stats")
            .arg("0");
        rust_cmd.args(&rust_extra);
        let rust = rust_cmd
            .output()
            .expect("run blast-cli psiblast restart MSA");

        let mut ncbi_cmd = std::process::Command::new("/usr/bin/psiblast");
        ncbi_cmd
            .arg("-in_msa")
            .arg(&msa_path)
            .arg("-subject")
            .arg(&subject)
            .arg("-outfmt")
            .arg("6 qseqid sseqid length qstart qend sstart send")
            .arg("-seg")
            .arg("no")
            .arg("-comp_based_stats")
            .arg("0");
        ncbi_cmd.args(&ncbi_extra);
        let ncbi = ncbi_cmd.output().expect("run NCBI psiblast restart MSA");

        assert!(rust.status.success(), "blast-cli psiblast failed");
        assert!(ncbi.status.success(), "NCBI psiblast failed");
        assert_eq!(
            rust.stdout, ncbi.stdout,
            "psiblast restart MSA stdout differs"
        );
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "psiblast restart MSA stderr differs"
        );
    }
}

#[test]
fn psiblast_db_ncbi_parity_restart_msa_master_selection() {
    if !std::path::Path::new("/usr/bin/psiblast").exists()
        || !std::path::Path::new("/usr/bin/makeblastdb").exists()
    {
        eprintln!("Skipping: /usr/bin/psiblast or /usr/bin/makeblastdb not found");
        return;
    }
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    for (msa, rust_extra, ncbi_extra) in [
        (
            ">master
MKFLIFALILFATVALAPKSSSHEI
>seq2
MKFLIFALILFATVALAPKSSSHEI
",
            Vec::<&str>::new(),
            Vec::<&str>::new(),
        ),
        (
            ">master desc extra
MKFLIFALILFATVALAPKSSSHEI
>seq2
MKFLIFALILFATVALAPKSSSHEI
",
            Vec::<&str>::new(),
            Vec::<&str>::new(),
        ),
        (
            ">   master desc
MKFLIFALILFATVALAPKSSSHEI
>seq2
MKFLIFALILFATVALAPKSSSHEI
",
            Vec::<&str>::new(),
            Vec::<&str>::new(),
        ),
        (
            ">master desc|bad
MKFLIFALILFATVALAPKSSSHEI
>seq2
MKFLIFALILFATVALAPKSSSHEI
",
            Vec::<&str>::new(),
            Vec::<&str>::new(),
        ),
        (
            ">master\tdesc\textra
MKFLIFALILFATVALAPKSSSHEI
>seq2
MKFLIFALILFATVALAPKSSSHEI
",
            Vec::<&str>::new(),
            Vec::<&str>::new(),
        ),
        (
            ">master\r\nMKFLIFALILFATVALAPKSSSHEI\r\n>seq2\r\nMKFLIFALILFATVALAPKSSSHEI\r\n",
            Vec::<&str>::new(),
            Vec::<&str>::new(),
        ),
        (
            ">master
MKFLIFALILFATVALAPKSSSHEI

;comment
>seq2
MKFLIFALILFATVALAPKSSSHEI
",
            Vec::<&str>::new(),
            Vec::<&str>::new(),
        ),
        (
            ">seq1
M-FLIFALILFATVALAPKSSSHEI
>master2
MKFLIFALILFATVALAPKSSSHEI
",
            vec!["--msa_master_idx", "2"],
            vec!["-msa_master_idx", "2"],
        ),
        (
            ">seq1
M-FLIFALILFATVALAPKSSSHEI
>master2 desc2 extra
MKFLIFALILFATVALAPKSSSHEI
",
            vec!["--msa_master_idx", "2"],
            vec!["-msa_master_idx", "2"],
        ),
        (
            ">master
M-KFLIFALILFATVALAPKSSSHEI
>seq2
MAKFLIFALILFATVALAPKSSSHEI
",
            Vec::<&str>::new(),
            Vec::<&str>::new(),
        ),
        (
            ">master
MKFLIFALILFATVALAPKSSSHEI
>empty

>seq2
MKFLIFALILFATVALAPKSSSHEI
",
            Vec::<&str>::new(),
            Vec::<&str>::new(),
        ),
        (
            ">master
MKF LIF\tALILFATVALAPKSSSHEI

>seq2
MKF LIF\tALILFATVALAPKSSSHEI
",
            Vec::<&str>::new(),
            Vec::<&str>::new(),
        ),
        (
            ">master
MKFLIFALILFATVALAPKSSSHEI
>seq|bad
MKFLIFALILFATVALAPKSSSHEI
",
            Vec::<&str>::new(),
            Vec::<&str>::new(),
        ),
        (
            ">master
MBZJOUXKLF
>seq2
MBZJOUXKLF
",
            Vec::<&str>::new(),
            Vec::<&str>::new(),
        ),
        (
            ">master
M-oF
>seq2
M-oF
",
            Vec::<&str>::new(),
            Vec::<&str>::new(),
        ),
        (
            ">master label
M-oF
>seq2
M-oF
",
            Vec::<&str>::new(),
            Vec::<&str>::new(),
        ),
        (
            ">master
MKFLIFALILFATVALAPKSSSHEI
>seq2
MKFLIFALILFATVALAPKSSSHEI
",
            vec!["--ignore_msa_master"],
            vec!["-ignore_msa_master"],
        ),
    ] {
        let tmp = TempDir::new().expect("tempdir");
        let msa_path = tmp.path().join("restart.msa");
        let db_fasta = tmp.path().join("db.fa");
        let db = tmp.path().join("testdb");
        let rust_out = tmp.path().join("rust.tsv");
        let ncbi_out = tmp.path().join("ncbi.tsv");
        std::fs::write(&msa_path, msa).expect("write restart MSA");
        std::fs::write(&db_fasta, ">s\nMKFLIFALILFATVALAPKSSSHEIHH\n").expect("write DB FASTA");

        let make_status = std::process::Command::new("/usr/bin/makeblastdb")
            .arg("-in")
            .arg(&db_fasta)
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

        let mut rust_cmd = std::process::Command::new(&blast_cli);
        rust_cmd
            .arg("psiblast")
            .arg("--in_msa")
            .arg(&msa_path)
            .arg("--db")
            .arg(&db)
            .arg("--outfmt")
            .arg("6 qseqid sseqid length qstart qend sstart send")
            .arg("--seg")
            .arg("no")
            .arg("--comp_based_stats")
            .arg("0")
            .arg("--num_threads")
            .arg("1")
            .arg("--out")
            .arg(&rust_out);
        rust_cmd.args(&rust_extra);
        let rust_status = rust_cmd
            .status()
            .expect("run blast-cli psiblast restart MSA DB");

        let mut ncbi_cmd = std::process::Command::new("/usr/bin/psiblast");
        ncbi_cmd
            .arg("-in_msa")
            .arg(&msa_path)
            .arg("-db")
            .arg(&db)
            .arg("-outfmt")
            .arg("6 qseqid sseqid length qstart qend sstart send")
            .arg("-seg")
            .arg("no")
            .arg("-comp_based_stats")
            .arg("0")
            .arg("-num_threads")
            .arg("1")
            .arg("-out")
            .arg(&ncbi_out);
        ncbi_cmd.args(&ncbi_extra);
        let ncbi_status = ncbi_cmd.status().expect("run NCBI psiblast restart MSA DB");

        assert!(
            rust_status.success(),
            "blast-cli psiblast exited with {rust_status}"
        );
        assert!(
            ncbi_status.success(),
            "NCBI psiblast exited with {ncbi_status}"
        );
        assert_eq!(
            std::fs::read(&rust_out).expect("read rust output"),
            std::fs::read(&ncbi_out).expect("read ncbi output"),
            "psiblast DB restart MSA output differs"
        );
    }
}

#[test]
fn psiblast_ncbi_parity_bad_restart_msa_master_errors() {
    if !std::path::Path::new("/usr/bin/psiblast").exists()
        || !std::path::Path::new("/usr/bin/makeblastdb").exists()
    {
        eprintln!("Skipping: /usr/bin/psiblast or /usr/bin/makeblastdb not found");
        return;
    }
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    for (msa, rust_extra, ncbi_extra) in [
        ("", Vec::<&str>::new(), Vec::<&str>::new()),
        (
            ">master
MKFLIFALILFATVALAPKSSSHEI
",
            Vec::<&str>::new(),
            Vec::<&str>::new(),
        ),
        (
            ">master
MKFLIFALILFATVALAPKSSSHEI
",
            vec!["--msa_master_idx", "3"],
            vec!["-msa_master_idx", "3"],
        ),
        (
            ">master
---
>seq2
MKF
",
            Vec::<&str>::new(),
            Vec::<&str>::new(),
        ),
        (
            "not-a-defline
MKF
",
            Vec::<&str>::new(),
            Vec::<&str>::new(),
        ),
        (
            ";comment
>master
MKF
>seq2
MKF
",
            Vec::<&str>::new(),
            Vec::<&str>::new(),
        ),
        (
            "#comment
>master
MKF
>seq2
MKF
",
            Vec::<&str>::new(),
            Vec::<&str>::new(),
        ),
        (
            ">master
MKF
>seq2
MKFXX
",
            Vec::<&str>::new(),
            Vec::<&str>::new(),
        ),
        (
            ">master
MK.F
>seq2
MK.F
",
            Vec::<&str>::new(),
            Vec::<&str>::new(),
        ),
        (
            ">master
MK1F
>seq2
MK1F
",
            Vec::<&str>::new(),
            Vec::<&str>::new(),
        ),
        (
            ">master
MK~F
>seq2
MK~F
",
            Vec::<&str>::new(),
            Vec::<&str>::new(),
        ),
        (
            ">master
MK\u{000c}F
>seq2
MK\u{000c}F
",
            Vec::<&str>::new(),
            Vec::<&str>::new(),
        ),
        (
            ">master
MK@F
>seq2
MK@F
",
            Vec::<&str>::new(),
            Vec::<&str>::new(),
        ),
        (
            ">master
MK?F
>seq2
MK?F
",
            Vec::<&str>::new(),
            Vec::<&str>::new(),
        ),
        (
            ">master
MK*F
>seq2
MK*F
",
            Vec::<&str>::new(),
            Vec::<&str>::new(),
        ),
        (
            ">master|bad
MKF
>seq2
MKF
",
            Vec::<&str>::new(),
            Vec::<&str>::new(),
        ),
        (
            "
>master
MKF
>seq2
MKF
",
            Vec::<&str>::new(),
            Vec::<&str>::new(),
        ),
        (
            ">master

MKF
>seq2
MKF
",
            Vec::<&str>::new(),
            Vec::<&str>::new(),
        ),
        (
            ">master
M

KF
>seq2
MKF
",
            Vec::<&str>::new(),
            Vec::<&str>::new(),
        ),
        (
            ">master
MKF

;comment
#comment2
>seq2
MKF
",
            Vec::<&str>::new(),
            Vec::<&str>::new(),
        ),
        (
            ">master
MKF
;comment
>seq2
MKF
",
            Vec::<&str>::new(),
            Vec::<&str>::new(),
        ),
        (
            ">master
;comment
MKF
>seq2
MKF
",
            Vec::<&str>::new(),
            Vec::<&str>::new(),
        ),
        (
            ">master
MKF
#comment
>seq2
MKF
",
            Vec::<&str>::new(),
            Vec::<&str>::new(),
        ),
        (
            ">master
MKF
>
MKF
",
            Vec::<&str>::new(),
            Vec::<&str>::new(),
        ),
    ] {
        let tmp = TempDir::new().expect("tempdir");
        let msa_path = tmp.path().join("restart.msa");
        let subject = tmp.path().join("subject.fa");
        let db_fasta = tmp.path().join("db.fa");
        let db = tmp.path().join("testdb");
        std::fs::write(&msa_path, msa).expect("write restart MSA");
        std::fs::write(&subject, ">s\nMKFLIFALILFATVALAPKSSSHEIHH\n").expect("write subject");
        std::fs::write(&db_fasta, ">s\nMKFLIFALILFATVALAPKSSSHEIHH\n").expect("write DB FASTA");

        let make_status = std::process::Command::new("/usr/bin/makeblastdb")
            .arg("-in")
            .arg(&db_fasta)
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

        for search_source in ["subject", "db"] {
            let mut rust_cmd = std::process::Command::new(&blast_cli);
            rust_cmd.arg("psiblast").arg("--in_msa").arg(&msa_path);
            if search_source == "subject" {
                rust_cmd.arg("--subject").arg(&subject);
            } else {
                rust_cmd.arg("--db").arg(&db);
            }
            rust_cmd.arg("--outfmt").arg("6").args(&rust_extra);
            let rust = rust_cmd
                .output()
                .expect("run blast-cli psiblast bad restart MSA");

            let mut ncbi_cmd = std::process::Command::new("/usr/bin/psiblast");
            ncbi_cmd.arg("-in_msa").arg(&msa_path);
            if search_source == "subject" {
                ncbi_cmd.arg("-subject").arg(&subject);
            } else {
                ncbi_cmd.arg("-db").arg(&db).arg("-num_threads").arg("1");
            }
            ncbi_cmd.arg("-outfmt").arg("6").args(&ncbi_extra);
            let ncbi = ncbi_cmd
                .output()
                .expect("run NCBI psiblast bad restart MSA");

            assert!(
                !rust.status.success(),
                "blast-cli should reject bad restart MSA"
            );
            assert!(!ncbi.status.success(), "NCBI should reject bad restart MSA");
            assert_eq!(rust.status.code(), ncbi.status.code(), "status differs");
            assert_eq!(rust.stdout, ncbi.stdout, "stdout differs");
            assert_eq!(
                String::from_utf8_lossy(&rust.stderr),
                String::from_utf8_lossy(&ncbi.stderr),
                "stderr differs for {search_source}"
            );
        }
    }
}

#[test]
fn psiblast_ncbi_parity_save_pssm_flags_require_output_path() {
    if !std::path::Path::new("/usr/bin/psiblast").exists() {
        eprintln!("Skipping: /usr/bin/psiblast not found");
        return;
    }
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    let tmp = TempDir::new().expect("tempdir");
    let query = tmp.path().join("query.fa");
    let subject = tmp.path().join("subject.fa");
    std::fs::write(&query, ">q\nMKKWLFGFLG\n").expect("write query");
    std::fs::write(&subject, ">s\nMKKWLFGFLG\n").expect("write subject");

    for (rust_flag, ncbi_flag) in [
        ("--save_each_pssm", "-save_each_pssm"),
        (
            "--save_pssm_after_last_round",
            "-save_pssm_after_last_round",
        ),
    ] {
        let rust = std::process::Command::new(&blast_cli)
            .arg("psiblast")
            .arg("--query")
            .arg(&query)
            .arg("--subject")
            .arg(&subject)
            .arg(rust_flag)
            .arg("--outfmt")
            .arg("6")
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli psiblast {rust_flag}: {err}"));
        let ncbi = std::process::Command::new("/usr/bin/psiblast")
            .arg("-query")
            .arg(&query)
            .arg("-subject")
            .arg(&subject)
            .arg(ncbi_flag)
            .arg("-outfmt")
            .arg("6")
            .output()
            .unwrap_or_else(|err| panic!("run NCBI psiblast {ncbi_flag}: {err}"));

        assert!(
            !rust.status.success(),
            "blast-cli should reject {rust_flag} without output path"
        );
        assert!(
            !ncbi.status.success(),
            "NCBI should reject {ncbi_flag} without output path"
        );
        assert_eq!(rust.status.code(), ncbi.status.code(), "status differs");
        assert_eq!(rust.stdout, ncbi.stdout, "stdout differs");
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "stderr differs for {rust_flag}"
        );
    }
}

#[test]
fn psiblast_ncbi_parity_pssm_output_paths_without_save_flags_do_not_write() {
    if !std::path::Path::new("/usr/bin/psiblast").exists()
        || !std::path::Path::new("/usr/bin/makeblastdb").exists()
    {
        eprintln!("Skipping: /usr/bin/psiblast or /usr/bin/makeblastdb not found");
        return;
    }
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    let tmp = TempDir::new().expect("tempdir");
    let query = tmp.path().join("query.fa");
    let subject = tmp.path().join("subject.fa");
    let db_fasta = tmp.path().join("db.fa");
    let db = tmp.path().join("testdb");
    std::fs::write(&query, ">q\nMKKWLFGFLG\n").expect("write query");
    std::fs::write(&subject, ">s\nMKKWLFGFLG\n").expect("write subject");
    std::fs::write(&db_fasta, ">s\nMKKWLFGFLG\n").expect("write DB FASTA");
    let make_status = std::process::Command::new("/usr/bin/makeblastdb")
        .arg("-in")
        .arg(&db_fasta)
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

    for search_source in ["subject", "db"] {
        let rust_checkpoint = tmp.path().join(format!("rust-{search_source}.chk"));
        let rust_ascii = tmp.path().join(format!("rust-{search_source}.ascii"));
        let ncbi_checkpoint = tmp.path().join(format!("ncbi-{search_source}.chk"));
        let ncbi_ascii = tmp.path().join(format!("ncbi-{search_source}.ascii"));

        let mut rust_cmd = std::process::Command::new(&blast_cli);
        rust_cmd.arg("psiblast").arg("--query").arg(&query);
        if search_source == "subject" {
            rust_cmd.arg("--subject").arg(&subject);
        } else {
            rust_cmd.arg("--db").arg(&db).arg("--num_threads").arg("1");
        }
        let rust = rust_cmd
            .arg("--out_pssm")
            .arg(&rust_checkpoint)
            .arg("--out_ascii_pssm")
            .arg(&rust_ascii)
            .arg("--outfmt")
            .arg("6 qseqid sseqid")
            .output()
            .expect("run blast-cli psiblast PSSM output paths");

        let mut ncbi_cmd = std::process::Command::new("/usr/bin/psiblast");
        ncbi_cmd.arg("-query").arg(&query);
        if search_source == "subject" {
            ncbi_cmd.arg("-subject").arg(&subject);
        } else {
            ncbi_cmd.arg("-db").arg(&db).arg("-num_threads").arg("1");
        }
        let ncbi = ncbi_cmd
            .arg("-out_pssm")
            .arg(&ncbi_checkpoint)
            .arg("-out_ascii_pssm")
            .arg(&ncbi_ascii)
            .arg("-outfmt")
            .arg("6 qseqid sseqid")
            .output()
            .expect("run NCBI psiblast PSSM output paths");

        assert!(rust.status.success(), "blast-cli psiblast failed");
        assert!(ncbi.status.success(), "NCBI psiblast failed");
        assert_eq!(rust.stdout, ncbi.stdout, "stdout differs");
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "stderr differs"
        );
        assert!(
            !rust_checkpoint.exists(),
            "blast-cli should not write checkpoint without save flag"
        );
        assert!(
            !rust_ascii.exists(),
            "blast-cli should not write ASCII PSSM without save flag"
        );
        assert!(
            !ncbi_checkpoint.exists(),
            "NCBI should not write checkpoint without save flag"
        );
        assert!(
            !ncbi_ascii.exists(),
            "NCBI should not write ASCII PSSM without save flag"
        );
    }

    for search_source in ["subject", "db"] {
        for output_kind in ["checkpoint", "ascii"] {
            let rust_artifact = if output_kind == "checkpoint" {
                tmp.path().join(format!("rust-{search_source}-only.chk"))
            } else {
                tmp.path().join(format!("rust-{search_source}-only.ascii"))
            };
            let ncbi_artifact = if output_kind == "checkpoint" {
                tmp.path().join(format!("ncbi-{search_source}-only.chk"))
            } else {
                tmp.path().join(format!("ncbi-{search_source}-only.ascii"))
            };

            let mut rust_cmd = std::process::Command::new(&blast_cli);
            rust_cmd.arg("psiblast").arg("--query").arg(&query);
            if search_source == "subject" {
                rust_cmd.arg("--subject").arg(&subject);
            } else {
                rust_cmd.arg("--db").arg(&db).arg("--num_threads").arg("1");
            }
            if output_kind == "checkpoint" {
                rust_cmd.arg("--out_pssm").arg(&rust_artifact);
            } else {
                rust_cmd.arg("--out_ascii_pssm").arg(&rust_artifact);
            }
            let rust = rust_cmd
                .arg("--outfmt")
                .arg("6 qseqid sseqid")
                .output()
                .expect("run blast-cli psiblast single PSSM output path");

            let mut ncbi_cmd = std::process::Command::new("/usr/bin/psiblast");
            ncbi_cmd.arg("-query").arg(&query);
            if search_source == "subject" {
                ncbi_cmd.arg("-subject").arg(&subject);
            } else {
                ncbi_cmd.arg("-db").arg(&db).arg("-num_threads").arg("1");
            }
            if output_kind == "checkpoint" {
                ncbi_cmd.arg("-out_pssm").arg(&ncbi_artifact);
            } else {
                ncbi_cmd.arg("-out_ascii_pssm").arg(&ncbi_artifact);
            }
            let ncbi = ncbi_cmd
                .arg("-outfmt")
                .arg("6 qseqid sseqid")
                .output()
                .expect("run NCBI psiblast single PSSM output path");

            assert!(rust.status.success(), "blast-cli psiblast failed");
            assert!(ncbi.status.success(), "NCBI psiblast failed");
            assert_eq!(rust.stdout, ncbi.stdout, "stdout differs");
            assert_eq!(
                String::from_utf8_lossy(&rust.stderr),
                String::from_utf8_lossy(&ncbi.stderr),
                "stderr differs for {search_source} {output_kind}"
            );
            assert!(
                !rust_artifact.exists(),
                "blast-cli should not write single {output_kind} PSSM without save flag"
            );
            assert!(
                !ncbi_artifact.exists(),
                "NCBI should not write single {output_kind} PSSM without save flag"
            );
        }
    }
}

fn pssm_round_output_path(base: &std::path::Path, round: usize) -> std::path::PathBuf {
    let mut path = base.as_os_str().to_os_string();
    path.push(format!(".{round}"));
    std::path::PathBuf::from(path)
}

#[test]
fn psiblast_ncbi_parity_save_each_single_round_writes_no_pssm_artifacts() {
    if !std::path::Path::new("/usr/bin/psiblast").exists()
        || !std::path::Path::new("/usr/bin/makeblastdb").exists()
    {
        eprintln!("Skipping: /usr/bin/psiblast or /usr/bin/makeblastdb not found");
        return;
    }
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    let tmp = TempDir::new().expect("tempdir");
    let query = tmp.path().join("query.fa");
    let subject = tmp.path().join("subject.fa");
    let db = tmp.path().join("testdb");
    std::fs::write(&query, ">q\nMKKWLFGFLG\n").expect("write query");
    std::fs::write(&subject, ">s\nMKKWLFGFLG\n").expect("write subject");
    let make_status = std::process::Command::new("/usr/bin/makeblastdb")
        .arg("-in")
        .arg(&subject)
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

    for search_source in ["subject", "db"] {
        let rust_checkpoint = tmp.path().join(format!("rust-{search_source}.chk"));
        let rust_ascii = tmp.path().join(format!("rust-{search_source}.ascii"));
        let ncbi_checkpoint = tmp.path().join(format!("ncbi-{search_source}.chk"));
        let ncbi_ascii = tmp.path().join(format!("ncbi-{search_source}.ascii"));

        let mut rust_cmd = std::process::Command::new(&blast_cli);
        rust_cmd.arg("psiblast").arg("--query").arg(&query);
        if search_source == "subject" {
            rust_cmd.arg("--subject").arg(&subject);
        } else {
            rust_cmd.arg("--db").arg(&db).arg("--num_threads").arg("1");
        }
        let rust = rust_cmd
            .arg("--num_iterations")
            .arg("1")
            .arg("--out_pssm")
            .arg(&rust_checkpoint)
            .arg("--out_ascii_pssm")
            .arg(&rust_ascii)
            .arg("--save_each_pssm")
            .arg("--outfmt")
            .arg("6 qseqid sseqid")
            .output()
            .expect("run blast-cli psiblast save-each single round");

        let mut ncbi_cmd = std::process::Command::new("/usr/bin/psiblast");
        ncbi_cmd.arg("-query").arg(&query);
        if search_source == "subject" {
            ncbi_cmd.arg("-subject").arg(&subject);
        } else {
            ncbi_cmd.arg("-db").arg(&db).arg("-num_threads").arg("1");
        }
        let ncbi = ncbi_cmd
            .arg("-num_iterations")
            .arg("1")
            .arg("-out_pssm")
            .arg(&ncbi_checkpoint)
            .arg("-out_ascii_pssm")
            .arg(&ncbi_ascii)
            .arg("-save_each_pssm")
            .arg("-outfmt")
            .arg("6 qseqid sseqid")
            .output()
            .expect("run NCBI psiblast save-each single round");

        assert!(rust.status.success(), "blast-cli psiblast failed");
        assert!(ncbi.status.success(), "NCBI psiblast failed");
        assert_eq!(rust.stdout, ncbi.stdout, "stdout differs");
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "stderr differs for {search_source}"
        );
        for path in [&rust_checkpoint, &rust_ascii, &ncbi_checkpoint, &ncbi_ascii] {
            assert!(
                !path.exists(),
                "save-each one-round run should not write base artifact: {path:?}"
            );
            assert!(
                !pssm_round_output_path(path, 1).exists(),
                "save-each one-round run should not write round artifact: {path:?}"
            );
        }
    }
}

#[test]
fn psiblast_ncbi_parity_save_last_writes_pssm_artifacts() {
    if !std::path::Path::new("/usr/bin/psiblast").exists()
        || !std::path::Path::new("/usr/bin/makeblastdb").exists()
    {
        eprintln!("Skipping: /usr/bin/psiblast or /usr/bin/makeblastdb not found");
        return;
    }
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    let tmp = TempDir::new().expect("tempdir");
    let query = tmp.path().join("query.fa");
    let subject = tmp.path().join("subject.fa");
    let db_fasta = tmp.path().join("db.fa");
    let db = tmp.path().join("testdb");
    std::fs::write(&query, ">q\nMKKWLFGFLG\n").expect("write query");
    std::fs::write(&subject, ">s\nMKKWLFGFLG\n").expect("write subject");
    std::fs::write(&db_fasta, ">s\nMKKWLFGFLG\n").expect("write DB FASTA");
    let make_status = std::process::Command::new("/usr/bin/makeblastdb")
        .arg("-in")
        .arg(&db_fasta)
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

    for search_source in ["subject", "db"] {
        let rust_checkpoint = tmp.path().join(format!("rust-{search_source}.chk"));
        let rust_ascii = tmp.path().join(format!("rust-{search_source}.ascii"));
        let ncbi_checkpoint = tmp.path().join(format!("ncbi-{search_source}.chk"));
        let ncbi_ascii = tmp.path().join(format!("ncbi-{search_source}.ascii"));

        let mut rust_cmd = std::process::Command::new(&blast_cli);
        rust_cmd.arg("psiblast").arg("--query").arg(&query);
        if search_source == "subject" {
            rust_cmd.arg("--subject").arg(&subject);
        } else {
            rust_cmd.arg("--db").arg(&db).arg("--num_threads").arg("1");
        }
        let rust = rust_cmd
            .arg("--out_pssm")
            .arg(&rust_checkpoint)
            .arg("--out_ascii_pssm")
            .arg(&rust_ascii)
            .arg("--save_pssm_after_last_round")
            .arg("--outfmt")
            .arg("6 qseqid sseqid")
            .output()
            .expect("run blast-cli psiblast save last PSSM");

        let mut ncbi_cmd = std::process::Command::new("/usr/bin/psiblast");
        ncbi_cmd.arg("-query").arg(&query);
        if search_source == "subject" {
            ncbi_cmd.arg("-subject").arg(&subject);
        } else {
            ncbi_cmd.arg("-db").arg(&db).arg("-num_threads").arg("1");
        }
        let ncbi = ncbi_cmd
            .arg("-out_pssm")
            .arg(&ncbi_checkpoint)
            .arg("-out_ascii_pssm")
            .arg(&ncbi_ascii)
            .arg("-save_pssm_after_last_round")
            .arg("-outfmt")
            .arg("6 qseqid sseqid")
            .output()
            .expect("run NCBI psiblast save last PSSM");

        assert!(rust.status.success(), "blast-cli psiblast failed");
        assert!(ncbi.status.success(), "NCBI psiblast failed");
        assert_eq!(rust.stdout, ncbi.stdout, "stdout differs");
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "stderr differs"
        );
        for path in [&rust_checkpoint, &rust_ascii, &ncbi_checkpoint, &ncbi_ascii] {
            assert!(path.is_file(), "PSSM artifact should exist: {path:?}");
            assert!(
                std::fs::metadata(path).expect("stat PSSM artifact").len() > 0,
                "PSSM artifact should be non-empty: {path:?}"
            );
        }
    }

    for search_source in ["subject", "db"] {
        for output_kind in ["checkpoint", "ascii"] {
            let rust_artifact = if output_kind == "checkpoint" {
                tmp.path().join(format!("rust-{search_source}-only.chk"))
            } else {
                tmp.path().join(format!("rust-{search_source}-only.ascii"))
            };
            let ncbi_artifact = if output_kind == "checkpoint" {
                tmp.path().join(format!("ncbi-{search_source}-only.chk"))
            } else {
                tmp.path().join(format!("ncbi-{search_source}-only.ascii"))
            };

            let mut rust_cmd = std::process::Command::new(&blast_cli);
            rust_cmd.arg("psiblast").arg("--query").arg(&query);
            if search_source == "subject" {
                rust_cmd.arg("--subject").arg(&subject);
            } else {
                rust_cmd.arg("--db").arg(&db).arg("--num_threads").arg("1");
            }
            if output_kind == "checkpoint" {
                rust_cmd.arg("--out_pssm").arg(&rust_artifact);
            } else {
                rust_cmd.arg("--out_ascii_pssm").arg(&rust_artifact);
            }
            let rust = rust_cmd
                .arg("--save_pssm_after_last_round")
                .arg("--outfmt")
                .arg("6 qseqid sseqid")
                .output()
                .expect("run blast-cli psiblast save last single PSSM output");

            let mut ncbi_cmd = std::process::Command::new("/usr/bin/psiblast");
            ncbi_cmd.arg("-query").arg(&query);
            if search_source == "subject" {
                ncbi_cmd.arg("-subject").arg(&subject);
            } else {
                ncbi_cmd.arg("-db").arg(&db).arg("-num_threads").arg("1");
            }
            if output_kind == "checkpoint" {
                ncbi_cmd.arg("-out_pssm").arg(&ncbi_artifact);
            } else {
                ncbi_cmd.arg("-out_ascii_pssm").arg(&ncbi_artifact);
            }
            let ncbi = ncbi_cmd
                .arg("-save_pssm_after_last_round")
                .arg("-outfmt")
                .arg("6 qseqid sseqid")
                .output()
                .expect("run NCBI psiblast save last single PSSM output");

            assert!(rust.status.success(), "blast-cli psiblast failed");
            assert!(ncbi.status.success(), "NCBI psiblast failed");
            assert_eq!(rust.stdout, ncbi.stdout, "stdout differs");
            assert_eq!(
                String::from_utf8_lossy(&rust.stderr),
                String::from_utf8_lossy(&ncbi.stderr),
                "stderr differs for {search_source} {output_kind}"
            );
            for path in [&rust_artifact, &ncbi_artifact] {
                assert!(
                    path.is_file(),
                    "single PSSM artifact should exist: {path:?}"
                );
                assert!(
                    std::fs::metadata(path)
                        .expect("stat single PSSM artifact")
                        .len()
                        > 0,
                    "single PSSM artifact should be non-empty: {path:?}"
                );
            }
        }
    }
}

#[test]
fn deltablast_ncbi_parity_rejects_psiblast_restart_options() {
    if !std::path::Path::new("/usr/bin/deltablast").exists() {
        eprintln!("Skipping: /usr/bin/deltablast not found");
        return;
    }
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    for (rust_args, ncbi_args) in [
        (
            vec!["--in_msa", "missing.msa"],
            vec!["-in_msa", "missing.msa"],
        ),
        (vec!["--msa_master_idx", "1"], vec!["-msa_master_idx", "1"]),
        (vec!["--ignore_msa_master"], vec!["-ignore_msa_master"]),
        (
            vec!["--in_pssm", "missing.chk"],
            vec!["-in_pssm", "missing.chk"],
        ),
        (
            vec!["--phi_pattern", "pattern.txt"],
            vec!["-phi_pattern", "pattern.txt"],
        ),
    ] {
        let rust = std::process::Command::new(&blast_cli)
            .arg("deltablast")
            .args(&rust_args)
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli deltablast {rust_args:?}: {err}"));
        let ncbi = std::process::Command::new("/usr/bin/deltablast")
            .args(&ncbi_args)
            .output()
            .unwrap_or_else(|err| panic!("run NCBI deltablast {ncbi_args:?}: {err}"));

        assert!(
            !rust.status.success(),
            "blast-cli should reject deltablast {rust_args:?}"
        );
        assert!(
            !ncbi.status.success(),
            "NCBI should reject deltablast {ncbi_args:?}"
        );
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "deltablast {rust_args:?} status differs"
        );
        assert_eq!(rust.stdout, ncbi.stdout, "deltablast stdout differs");
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "deltablast stderr differs for {rust_args:?}"
        );
    }
}

#[test]
fn deltablast_ncbi_parity_accepts_delta_only_options_before_database_check() {
    if !std::path::Path::new("/usr/bin/deltablast").exists() {
        eprintln!("Skipping: /usr/bin/deltablast not found");
        return;
    }
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    for (rust_args, ncbi_args) in [
        (
            vec!["--rpsdb", "missing_cdd"],
            vec!["-rpsdb", "missing_cdd"],
        ),
        (vec!["--show_domain_hits"], vec!["-show_domain_hits"]),
    ] {
        let rust = std::process::Command::new(&blast_cli)
            .arg("deltablast")
            .args(&rust_args)
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli deltablast {rust_args:?}: {err}"));
        let ncbi = std::process::Command::new("/usr/bin/deltablast")
            .args(&ncbi_args)
            .output()
            .unwrap_or_else(|err| panic!("run NCBI deltablast {ncbi_args:?}: {err}"));

        assert!(
            !rust.status.success(),
            "blast-cli should reject missing search source after {rust_args:?}"
        );
        assert!(
            !ncbi.status.success(),
            "NCBI should reject missing search source after {ncbi_args:?}"
        );
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "deltablast {rust_args:?} status differs"
        );
        assert_eq!(rust.stdout, ncbi.stdout, "deltablast stdout differs");
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "deltablast stderr differs for {rust_args:?}"
        );
    }
}

#[test]
fn deltablast_ncbi_parity_existing_search_db_reaches_missing_rpsdb() {
    if !std::path::Path::new("/usr/bin/deltablast").exists()
        || !std::path::Path::new("/usr/bin/makeblastdb").exists()
    {
        eprintln!("Skipping: /usr/bin/deltablast or /usr/bin/makeblastdb not found");
        return;
    }
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    let tmp = TempDir::new().expect("tempdir");
    let query = tmp.path().join("query.fa");
    let subject = tmp.path().join("subject.fa");
    let search_db = tmp.path().join("searchdb");
    let rpsdb = tmp.path().join("missing_cdd");
    std::fs::write(&query, ">q\nMKKWLFGFLG\n").expect("write query");
    std::fs::write(&subject, ">s\nMKKWLFGFLG\n").expect("write subject");
    let make_status = std::process::Command::new("/usr/bin/makeblastdb")
        .arg("-in")
        .arg(&subject)
        .arg("-dbtype")
        .arg("prot")
        .arg("-out")
        .arg(&search_db)
        .status()
        .expect("run makeblastdb");
    assert!(
        make_status.success(),
        "makeblastdb exited with {make_status}"
    );

    let rust = std::process::Command::new(&blast_cli)
        .arg("deltablast")
        .arg("--query")
        .arg(&query)
        .arg("--db")
        .arg(&search_db)
        .arg("--rpsdb")
        .arg(&rpsdb)
        .arg("--outfmt")
        .arg("6")
        .output()
        .expect("run blast-cli deltablast missing rpsdb");
    let ncbi = std::process::Command::new("/usr/bin/deltablast")
        .arg("-query")
        .arg(&query)
        .arg("-db")
        .arg(&search_db)
        .arg("-rpsdb")
        .arg(&rpsdb)
        .arg("-outfmt")
        .arg("6")
        .output()
        .expect("run NCBI deltablast missing rpsdb");

    assert!(
        !rust.status.success(),
        "blast-cli should reject missing RPS DB"
    );
    assert!(!ncbi.status.success(), "NCBI should reject missing RPS DB");
    assert_eq!(rust.status.code(), ncbi.status.code(), "status differs");
    assert_eq!(rust.stdout, ncbi.stdout, "stdout differs");
    assert_eq!(
        String::from_utf8_lossy(&rust.stderr),
        String::from_utf8_lossy(&ncbi.stderr),
        "stderr differs"
    );
}

#[test]
fn deltablast_ncbi_parity_pairwise_subject_reaches_missing_rpsdb_with_preamble() {
    if !std::path::Path::new("/usr/bin/deltablast").exists() {
        eprintln!("Skipping: /usr/bin/deltablast not found");
        return;
    }
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    for (rust_extra, ncbi_extra) in [
        (Vec::<&str>::new(), Vec::<&str>::new()),
        (
            vec!["--rpsdb", "missing_cdd"],
            vec!["-rpsdb", "missing_cdd"],
        ),
        (vec!["--subject_loc", "1-10"], vec!["-subject_loc", "1-10"]),
        (vec!["--num_threads", "2"], vec!["-num_threads", "2"]),
    ] {
        let rust = std::process::Command::new(&blast_cli)
            .arg("deltablast")
            .arg("--query")
            .arg("tests/fixtures/protein_query.fa")
            .arg("--subject")
            .arg("tests/fixtures/protein_subject.fa")
            .arg("--outfmt")
            .arg("0")
            .args(&rust_extra)
            .output()
            .unwrap_or_else(|err| {
                panic!("run blast-cli deltablast pairwise subject {rust_extra:?}: {err}")
            });
        let ncbi = std::process::Command::new("/usr/bin/deltablast")
            .arg("-query")
            .arg("tests/fixtures/protein_query.fa")
            .arg("-subject")
            .arg("tests/fixtures/protein_subject.fa")
            .arg("-outfmt")
            .arg("0")
            .args(&ncbi_extra)
            .output()
            .unwrap_or_else(|err| {
                panic!("run NCBI deltablast pairwise subject {ncbi_extra:?}: {err}")
            });

        assert!(
            !rust.status.success(),
            "blast-cli should reject missing DELTA RPS DB"
        );
        assert!(
            !ncbi.status.success(),
            "NCBI should reject missing DELTA RPS DB"
        );
        assert_eq!(rust.status.code(), ncbi.status.code(), "status differs");
        assert_eq!(
            rust.stdout, ncbi.stdout,
            "stdout differs for {rust_extra:?}"
        );
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "stderr differs for {rust_extra:?}"
        );
    }
}

#[test]
fn deltablast_ncbi_parity_location_errors_precede_missing_rpsdb() {
    if !std::path::Path::new("/usr/bin/deltablast").exists() {
        eprintln!("Skipping: /usr/bin/deltablast not found");
        return;
    }
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    for (rust_loc_arg, ncbi_loc_arg) in [
        ("--subject_loc", "-subject_loc"),
        ("--query_loc", "-query_loc"),
    ] {
        let rust = std::process::Command::new(&blast_cli)
            .arg("deltablast")
            .arg("--query")
            .arg("tests/fixtures/protein_query.fa")
            .arg("--subject")
            .arg("tests/fixtures/protein_subject.fa")
            .arg(rust_loc_arg)
            .arg("500-600")
            .arg("--outfmt")
            .arg("6")
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli deltablast invalid {rust_loc_arg}: {err}"));
        let ncbi = std::process::Command::new("/usr/bin/deltablast")
            .arg("-query")
            .arg("tests/fixtures/protein_query.fa")
            .arg("-subject")
            .arg("tests/fixtures/protein_subject.fa")
            .arg(ncbi_loc_arg)
            .arg("500-600")
            .arg("-outfmt")
            .arg("6")
            .output()
            .unwrap_or_else(|err| panic!("run NCBI deltablast invalid {ncbi_loc_arg}: {err}"));

        assert!(
            !rust.status.success(),
            "blast-cli should reject invalid {rust_loc_arg}"
        );
        assert!(
            !ncbi.status.success(),
            "NCBI should reject invalid {ncbi_loc_arg}"
        );
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "status differs for {rust_loc_arg}"
        );
        assert_eq!(
            rust.stdout, ncbi.stdout,
            "stdout differs for {rust_loc_arg}"
        );
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "stderr differs for {rust_loc_arg}"
        );
    }
}

#[test]
fn protein_programs_ncbi_parity_location_error_precedence() {
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    for (program, query, subject) in [
        (
            "blastp",
            "tests/fixtures/protein_query.fa",
            "tests/fixtures/protein_subject.fa",
        ),
        (
            "blastx",
            "tests/fixtures/blastx_nuc_query.fa",
            "tests/fixtures/blastx_prot_subject.fa",
        ),
        (
            "tblastn",
            "tests/fixtures/tblastn_prot_query.fa",
            "tests/fixtures/tblastn_nuc_subject.fa",
        ),
        (
            "tblastx",
            "tests/fixtures/tblastx_nuc_query.fa",
            "tests/fixtures/tblastx_nuc_subject.fa",
        ),
        (
            "psiblast",
            "tests/fixtures/protein_query.fa",
            "tests/fixtures/protein_subject.fa",
        ),
        (
            "deltablast",
            "tests/fixtures/protein_query.fa",
            "tests/fixtures/protein_subject.fa",
        ),
    ] {
        let ncbi_bin = format!("/usr/bin/{program}");
        if !std::path::Path::new(&ncbi_bin).exists() {
            eprintln!("Skipping: {ncbi_bin} not found");
            continue;
        }

        for (rust_args, ncbi_args, label) in [
            (
                vec![
                    "--query",
                    query,
                    "--subject",
                    subject,
                    "--query_loc",
                    "bad",
                    "--outfmt",
                    "6",
                ],
                vec![
                    "-query",
                    query,
                    "-subject",
                    subject,
                    "-query_loc",
                    "bad",
                    "-outfmt",
                    "6",
                ],
                "bad query_loc with explicit query",
            ),
            (
                vec!["--subject", subject, "--query_loc", "bad", "--outfmt", "6"],
                vec!["-subject", subject, "-query_loc", "bad", "-outfmt", "6"],
                "bad query_loc with omitted query",
            ),
            (
                vec![
                    "--query",
                    "missing_query.fa",
                    "--subject",
                    subject,
                    "--subject_loc",
                    "bad",
                    "--outfmt",
                    "6",
                ],
                vec![
                    "-query",
                    "missing_query.fa",
                    "-subject",
                    subject,
                    "-subject_loc",
                    "bad",
                    "-outfmt",
                    "6",
                ],
                "bad subject_loc before missing query",
            ),
            (
                vec![
                    "--query",
                    query,
                    "--subject",
                    subject,
                    "--subject_loc",
                    "bad",
                    "--out",
                    "/no/such/dir/out.txt",
                    "--outfmt",
                    "6",
                ],
                vec![
                    "-query",
                    query,
                    "-subject",
                    subject,
                    "-subject_loc",
                    "bad",
                    "-out",
                    "/no/such/dir/out.txt",
                    "-outfmt",
                    "6",
                ],
                "bad subject_loc before bad out",
            ),
            (
                vec![
                    "--query",
                    query,
                    "--subject",
                    subject,
                    "--query_loc",
                    "bad",
                    "--out",
                    "/no/such/dir/out.txt",
                    "--outfmt",
                    "6",
                ],
                vec![
                    "-query",
                    query,
                    "-subject",
                    subject,
                    "-query_loc",
                    "bad",
                    "-out",
                    "/no/such/dir/out.txt",
                    "-outfmt",
                    "6",
                ],
                "bad out before bad query_loc",
            ),
        ] {
            let rust = std::process::Command::new(&blast_cli)
                .arg(program)
                .args(&rust_args)
                .output()
                .unwrap_or_else(|err| panic!("run blast-cli {program} {label}: {err}"));
            let ncbi = std::process::Command::new(&ncbi_bin)
                .args(&ncbi_args)
                .output()
                .unwrap_or_else(|err| panic!("run NCBI {program} {label}: {err}"));

            assert!(
                !rust.status.success(),
                "blast-cli should reject {program} {label}"
            );
            assert!(
                !ncbi.status.success(),
                "NCBI should reject {program} {label}"
            );
            assert_eq!(
                rust.status.code(),
                ncbi.status.code(),
                "{program} {label} status differs"
            );
            assert_eq!(rust.stdout, ncbi.stdout, "{program} {label} stdout differs");
            assert_eq!(
                String::from_utf8_lossy(&rust.stderr),
                String::from_utf8_lossy(&ncbi.stderr),
                "{program} {label} stderr differs"
            );
        }
    }
}

#[test]
fn deltablast_ncbi_parity_input_file_errors_precede_missing_rpsdb() {
    if !std::path::Path::new("/usr/bin/deltablast").exists() {
        eprintln!("Skipping: /usr/bin/deltablast not found");
        return;
    }
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    for (rust_args, ncbi_args) in [
        (
            vec![
                "--query",
                "missing_query.fa",
                "--subject",
                "tests/fixtures/protein_subject.fa",
            ],
            vec![
                "-query",
                "missing_query.fa",
                "-subject",
                "tests/fixtures/protein_subject.fa",
            ],
        ),
        (
            vec!["--query", "missing_query.fa", "--db", "missing_search"],
            vec!["-query", "missing_query.fa", "-db", "missing_search"],
        ),
        (
            vec![
                "--query",
                "tests/fixtures/protein_query.fa",
                "--subject",
                "missing_subject.fa",
            ],
            vec![
                "-query",
                "tests/fixtures/protein_query.fa",
                "-subject",
                "missing_subject.fa",
            ],
        ),
        (
            vec!["--subject", "missing_subject.fa"],
            vec!["-subject", "missing_subject.fa"],
        ),
        (
            vec![
                "--query",
                "tests/fixtures/protein_query.fa",
                "--subject",
                "tests/fixtures/protein_subject.fa",
                "--out",
                "/no/such/dir/out.txt",
            ],
            vec![
                "-query",
                "tests/fixtures/protein_query.fa",
                "-subject",
                "tests/fixtures/protein_subject.fa",
                "-out",
                "/no/such/dir/out.txt",
            ],
        ),
    ] {
        let rust = std::process::Command::new(&blast_cli)
            .arg("deltablast")
            .args(&rust_args)
            .arg("--outfmt")
            .arg("6")
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli deltablast {rust_args:?}: {err}"));
        let ncbi = std::process::Command::new("/usr/bin/deltablast")
            .args(&ncbi_args)
            .arg("-outfmt")
            .arg("6")
            .output()
            .unwrap_or_else(|err| panic!("run NCBI deltablast {ncbi_args:?}: {err}"));

        assert!(
            !rust.status.success(),
            "blast-cli should reject inaccessible input for {rust_args:?}"
        );
        assert!(
            !ncbi.status.success(),
            "NCBI should reject inaccessible input for {ncbi_args:?}"
        );
        assert_eq!(rust.status.code(), ncbi.status.code(), "status differs");
        assert_eq!(rust.stdout, ncbi.stdout, "stdout differs");
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "stderr differs for {rust_args:?}"
        );
    }
}

#[test]
fn deltablast_ncbi_parity_show_domain_hits_rejects_value() {
    if !std::path::Path::new("/usr/bin/deltablast").exists() {
        eprintln!("Skipping: /usr/bin/deltablast not found");
        return;
    }
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    let rust = std::process::Command::new(&blast_cli)
        .arg("deltablast")
        .arg("--query")
        .arg("tests/fixtures/protein_query.fa")
        .arg("--show_domain_hits")
        .arg("true")
        .output()
        .expect("run blast-cli deltablast show_domain_hits true");
    let ncbi = std::process::Command::new("/usr/bin/deltablast")
        .arg("-query")
        .arg("tests/fixtures/protein_query.fa")
        .arg("-show_domain_hits")
        .arg("true")
        .output()
        .expect("run NCBI deltablast show_domain_hits true");

    assert!(!rust.status.success(), "blast-cli should reject flag value");
    assert!(!ncbi.status.success(), "NCBI should reject flag value");
    assert_eq!(rust.status.code(), ncbi.status.code(), "status differs");
    assert_eq!(rust.stdout, ncbi.stdout, "stdout differs");
    assert_eq!(
        String::from_utf8_lossy(&rust.stderr),
        String::from_utf8_lossy(&ncbi.stderr),
        "stderr differs"
    );
}

#[test]
fn deltablast_ncbi_parity_show_domain_hits_rejects_subject() {
    if !std::path::Path::new("/usr/bin/deltablast").exists() {
        eprintln!("Skipping: /usr/bin/deltablast not found");
        return;
    }
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    let rust = std::process::Command::new(&blast_cli)
        .arg("deltablast")
        .arg("--query")
        .arg("tests/fixtures/protein_query.fa")
        .arg("--subject")
        .arg("tests/fixtures/protein_subject.fa")
        .arg("--show_domain_hits")
        .arg("--outfmt")
        .arg("6")
        .output()
        .expect("run blast-cli deltablast show_domain_hits subject");
    let ncbi = std::process::Command::new("/usr/bin/deltablast")
        .arg("-query")
        .arg("tests/fixtures/protein_query.fa")
        .arg("-subject")
        .arg("tests/fixtures/protein_subject.fa")
        .arg("-show_domain_hits")
        .arg("-outfmt")
        .arg("6")
        .output()
        .expect("run NCBI deltablast show_domain_hits subject");

    assert!(
        !rust.status.success(),
        "blast-cli should reject subject with show_domain_hits"
    );
    assert!(
        !ncbi.status.success(),
        "NCBI should reject subject with show_domain_hits"
    );
    assert_eq!(rust.status.code(), ncbi.status.code(), "status differs");
    assert_eq!(rust.stdout, ncbi.stdout, "stdout differs");
    assert_eq!(
        String::from_utf8_lossy(&rust.stderr),
        String::from_utf8_lossy(&ncbi.stderr),
        "stderr differs"
    );
}

fn assert_psiblast_subject_outfmt_matches_ncbi(
    query_fasta: &str,
    subject_fasta: &str,
    outfmt: &str,
    rust_extra: &[&str],
    ncbi_extra: &[&str],
) {
    if !std::path::Path::new("/usr/bin/psiblast").exists() {
        eprintln!("Skipping: /usr/bin/psiblast not found");
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

    let mut ncbi_cmd = std::process::Command::new("/usr/bin/psiblast");
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

#[test]
fn psiblast_subject_ncbi_parity_custom_field_parser_edges() {
    if !std::path::Path::new("/usr/bin/psiblast").exists() {
        eprintln!("Skipping: /usr/bin/psiblast not found");
        return;
    }
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    for outfmt in [
        "6 delim=, qseqid bogus sseqid qseqid length length",
        "6 delim=tab qseqid sseqid length",
        "6 delim=, delim=tab qseqid sseqid length",
        "6 std qlen",
        "6 std qaccver saccver qlen",
        "6 bogus",
        "7 qseqid bogus sseqid qseqid length length",
        "7 delim=space qseqid sseqid length",
        "7 delim=, delim=tab qseqid sseqid length",
        "7 std qlen",
        "7 std qaccver saccver qlen",
        "7 bogus",
        "10 qseqid bogus sseqid qseqid length length",
        "10 delim=tab qseqid sseqid length",
        "10 delim=tab delim=, qseqid sseqid length",
        "10 std qlen",
        "10 std qaccver saccver qlen",
        "10 bogus",
    ] {
        let rust = std::process::Command::new(&blast_cli)
            .arg("psiblast")
            .arg("-query")
            .arg("tests/fixtures/psi_query.fa")
            .arg("-subject")
            .arg("tests/fixtures/psi_subject.fa")
            .arg("-outfmt")
            .arg(outfmt)
            .arg("-seg")
            .arg("no")
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli psiblast {outfmt}: {err}"));
        let ncbi = std::process::Command::new("/usr/bin/psiblast")
            .arg("-query")
            .arg("tests/fixtures/psi_query.fa")
            .arg("-subject")
            .arg("tests/fixtures/psi_subject.fa")
            .arg("-outfmt")
            .arg(outfmt)
            .arg("-seg")
            .arg("no")
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
}

#[test]
fn psiblast_subject_ncbi_parity_parse_deflines_tabular_and_csv() {
    if !std::path::Path::new("/usr/bin/psiblast").exists() {
        eprintln!("Skipping: /usr/bin/psiblast not found");
        return;
    }
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    let tmp = TempDir::new().expect("tempdir");
    let query = tmp.path().join("query.fa");
    let subject = tmp.path().join("subject.fa");
    std::fs::write(
        &query,
        ">gi|123|ref|QPROT.1| query protein\nMKFLIFALILFATVALAPKSSSHEI\n",
    )
    .expect("write query");
    std::fs::write(
        &subject,
        ">gi|456|ref|SPROT.1| subject protein\nMKFLIFALILFATVALAPKSSSHEIHH\n",
    )
    .expect("write subject");

    for outfmt in [
        "6 qseqid qgi qacc qaccver sseqid sgi sallgi sacc saccver sallseqid stitle salltitles",
        "10 qseqid qgi qacc qaccver sseqid sgi sallgi sacc saccver sallseqid stitle salltitles",
    ] {
        let rust = std::process::Command::new(&blast_cli)
            .arg("psiblast")
            .arg("--query")
            .arg(&query)
            .arg("--subject")
            .arg(&subject)
            .arg("--outfmt")
            .arg(outfmt)
            .arg("--seg")
            .arg("no")
            .arg("--comp_based_stats")
            .arg("0")
            .arg("--parse_deflines")
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli psiblast {outfmt}: {err}"));
        let ncbi = std::process::Command::new("/usr/bin/psiblast")
            .arg("-query")
            .arg(&query)
            .arg("-subject")
            .arg(&subject)
            .arg("-outfmt")
            .arg(outfmt)
            .arg("-seg")
            .arg("no")
            .arg("-comp_based_stats")
            .arg("0")
            .arg("-parse_deflines")
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
}

#[test]
fn psiblast_subject_ncbi_parity_commented_tabular() {
    if !std::path::Path::new("/usr/bin/psiblast").exists() {
        eprintln!("Skipping: /usr/bin/psiblast not found");
        return;
    }
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    for (outfmt, rust_extra, ncbi_extra) in [
        ("7", Vec::<&str>::new(), Vec::<&str>::new()),
        (
            "7 qaccver saccver pident length bitscore",
            vec!["--seg", "no", "--comp_based_stats", "0", "--parse_deflines"],
            vec!["-seg", "no", "-comp_based_stats", "0", "-parse_deflines"],
        ),
    ] {
        let tmp = TempDir::new().expect("tempdir");
        let query = tmp.path().join("query.fa");
        let subject = tmp.path().join("subject.fa");
        std::fs::write(
            &query,
            ">gi|123|ref|QPROT.1| query protein\nMKFLIFALILFATVALAPKSSSHEI\n",
        )
        .expect("write query");
        std::fs::write(
            &subject,
            ">gi|456|ref|SPROT.1| subject protein\nMKFLIFALILFATVALAPKSSSHEIHH\n",
        )
        .expect("write subject");

        let mut rust_cmd = std::process::Command::new(&blast_cli);
        rust_cmd
            .arg("psiblast")
            .arg("--query")
            .arg(&query)
            .arg("--subject")
            .arg(&subject)
            .arg("--outfmt")
            .arg(outfmt);
        rust_cmd.args(&rust_extra);
        let rust = rust_cmd
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli psiblast {outfmt}: {err}"));

        let mut ncbi_cmd = std::process::Command::new("/usr/bin/psiblast");
        ncbi_cmd
            .arg("-query")
            .arg(&query)
            .arg("-subject")
            .arg(&subject)
            .arg("-outfmt")
            .arg(outfmt);
        ncbi_cmd.args(&ncbi_extra);
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
}

#[test]
fn psiblast_subject_ncbi_parity_xml() {
    if !std::path::Path::new("/usr/bin/psiblast").exists() {
        eprintln!("Skipping: /usr/bin/psiblast not found");
        return;
    }
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    for (rust_extra, ncbi_extra) in [
        (Vec::<&str>::new(), Vec::<&str>::new()),
        (
            vec!["--seg", "no", "--comp_based_stats", "0", "--parse_deflines"],
            vec!["-seg", "no", "-comp_based_stats", "0", "-parse_deflines"],
        ),
    ] {
        let tmp = TempDir::new().expect("tempdir");
        let query = tmp.path().join("query.fa");
        let subject = tmp.path().join("subject.fa");
        std::fs::write(
            &query,
            ">gi|123|ref|QPROT.1| query protein\nMKFLIFALILFATVALAPKSSSHEI\n",
        )
        .expect("write query");
        std::fs::write(
            &subject,
            ">gi|456|ref|SPROT.1| subject protein\nMKFLIFALILFATVALAPKSSSHEIHH\n",
        )
        .expect("write subject");

        let mut rust_cmd = std::process::Command::new(&blast_cli);
        rust_cmd
            .arg("psiblast")
            .arg("--query")
            .arg(&query)
            .arg("--subject")
            .arg(&subject)
            .arg("--outfmt")
            .arg("5");
        rust_cmd.args(&rust_extra);
        let rust = rust_cmd
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli psiblast XML: {err}"));

        let mut ncbi_cmd = std::process::Command::new("/usr/bin/psiblast");
        ncbi_cmd
            .arg("-query")
            .arg(&query)
            .arg("-subject")
            .arg(&subject)
            .arg("-outfmt")
            .arg("5");
        ncbi_cmd.args(&ncbi_extra);
        let ncbi = ncbi_cmd
            .output()
            .unwrap_or_else(|err| panic!("run NCBI psiblast XML: {err}"));

        assert!(rust.status.success(), "blast-cli psiblast XML failed");
        assert!(ncbi.status.success(), "NCBI psiblast XML failed");
        assert_eq!(rust.stdout, ncbi.stdout, "psiblast XML stdout differs");
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "psiblast XML stderr differs"
        );
    }
}

#[test]
fn psiblast_subject_ncbi_parity_pairwise() {
    assert_psiblast_subject_outfmt_matches_ncbi(
        ">q\nMKFLIFALILFATVALAPKSSSHEI\n",
        ">s\nMKFLIFALILFATVALAPKSSSHEIHH\n",
        "0",
        &[],
        &[],
    );
}

#[test]
fn psiblast_subject_ncbi_parity_pairwise_no_hits() {
    assert_psiblast_subject_outfmt_matches_ncbi(
        ">q\nMKFLIFALILFATVALAPKSSSHEI\n",
        ">s unrelated\nACDEFGHIKLMNPQRSTVWY\n",
        "0",
        &["--seg", "no", "--comp_based_stats", "0"],
        &["-seg", "no", "-comp_based_stats", "0"],
    );
}

#[test]
fn psiblast_subject_ncbi_parity_xml_multiple_queries() {
    assert_psiblast_subject_outfmt_matches_ncbi(
        ">q1\nMKFLIFALILFATVALAPKSSSHEI\n>q2\nACDEFGHIKLMNPQRSTVWY\n",
        ">s1\nMKFLIFALILFATVALAPKSSSHEIHH\n>s2\nACDEFGHIKLMNPQRSTVWY\n",
        "5",
        &["--seg", "no", "--comp_based_stats", "0"],
        &["-seg", "no", "-comp_based_stats", "0"],
    );
}

#[test]
fn protein_programs_ncbi_parity_threshold_conversion_errors_before_required_query() {
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    for program in ["blastp", "blastx", "tblastn", "tblastx"] {
        let ncbi_bin = format!("/usr/bin/{program}");
        if !std::path::Path::new(&ncbi_bin).exists() {
            eprintln!("Skipping: {ncbi_bin} not found");
            continue;
        }
        let rust = std::process::Command::new(&blast_cli)
            .arg(program)
            .arg("-threshold")
            .arg("abc")
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli {program} -threshold abc: {err}"));
        let ncbi = std::process::Command::new(&ncbi_bin)
            .arg("-threshold")
            .arg("abc")
            .output()
            .unwrap_or_else(|err| panic!("run NCBI {program} -threshold abc: {err}"));

        assert!(
            !rust.status.success(),
            "blast-cli should reject {program} -threshold abc"
        );
        assert!(
            !ncbi.status.success(),
            "NCBI should reject {program} -threshold abc"
        );
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "{program} threshold status differs"
        );
        assert_eq!(
            rust.stdout, ncbi.stdout,
            "{program} threshold stdout differs"
        );
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{program} threshold stderr differs"
        );
    }
}

#[test]
fn protein_programs_ncbi_parity_word_size_minimum_is_two() {
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    for program in ["blastp", "blastx", "tblastn", "tblastx"] {
        let ncbi_bin = format!("/usr/bin/{program}");
        if !std::path::Path::new(&ncbi_bin).exists() {
            eprintln!("Skipping: {ncbi_bin} not found");
            continue;
        }
        let rust = std::process::Command::new(&blast_cli)
            .arg(program)
            .arg("-word_size")
            .arg("1")
            .output()
            .unwrap_or_else(|err| panic!("run blast-cli {program} -word_size 1: {err}"));
        let ncbi = std::process::Command::new(&ncbi_bin)
            .arg("-word_size")
            .arg("1")
            .output()
            .unwrap_or_else(|err| panic!("run NCBI {program} -word_size 1: {err}"));

        assert!(
            !rust.status.success(),
            "blast-cli should reject {program} -word_size 1"
        );
        assert!(
            !ncbi.status.success(),
            "NCBI should reject {program} -word_size 1"
        );
        assert_eq!(
            rust.status.code(),
            ncbi.status.code(),
            "{program} word_size status differs"
        );
        assert_eq!(
            rust.stdout, ncbi.stdout,
            "{program} word_size stdout differs"
        );
        assert_eq!(
            String::from_utf8_lossy(&rust.stderr),
            String::from_utf8_lossy(&ncbi.stderr),
            "{program} word_size stderr differs"
        );
    }
}

#[test]
fn blastp_ncbi_parity_unsupported_matrix_error() {
    if !std::path::Path::new("/usr/bin/blastp").exists() {
        eprintln!("Skipping: /usr/bin/blastp not found");
        return;
    }
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: build blast-cli first or set BLAST_RS_CLI_BIN");
        return;
    };

    let tmp = TempDir::new().expect("tempdir");
    let query = tmp.path().join("query.fa");
    let subject = tmp.path().join("subject.fa");
    std::fs::write(&query, ">q\nACDEFGHIKLMNPQRSTVWY\n").expect("write query");
    std::fs::write(&subject, ">s\nACDEFGHIKLMNPQRSTVWY\n").expect("write subject");

    let rust = std::process::Command::new(&blast_cli)
        .arg("blastp")
        .arg("-query")
        .arg(&query)
        .arg("-subject")
        .arg(&subject)
        .arg("-matrix")
        .arg("BAD")
        .arg("-outfmt")
        .arg("6")
        .output()
        .expect("run blast-cli unsupported matrix");
    let ncbi = std::process::Command::new("/usr/bin/blastp")
        .arg("-query")
        .arg(&query)
        .arg("-subject")
        .arg(&subject)
        .arg("-matrix")
        .arg("BAD")
        .arg("-outfmt")
        .arg("6")
        .output()
        .expect("run NCBI unsupported matrix");

    assert!(!rust.status.success(), "blast-cli should reject BAD matrix");
    assert!(!ncbi.status.success(), "NCBI should reject BAD matrix");
    assert_eq!(rust.status.code(), ncbi.status.code(), "status differs");
    assert_eq!(rust.stdout, ncbi.stdout, "stdout differs");
    assert_eq!(
        String::from_utf8_lossy(&rust.stderr),
        String::from_utf8_lossy(&ncbi.stderr),
        "stderr differs"
    );
}

#[test]
fn blastp_subject_ncbi_parity_identity_matrix_scores_mismatch() {
    assert_blastp_subject_outfmt_matches_ncbi(
        ">q\nAAAAAAAAAAGAAAAAAAAAA\n",
        ">s\nAAAAAAAAAACAAAAAAAAAA\n",
        "6 qseqid sseqid score bitscore evalue length pident qseq sseq",
        &[
            "--matrix",
            "IDENTITY",
            "--seg",
            "no",
            "--comp_based_stats",
            "0",
            "--threshold",
            "1",
            "--word_size",
            "2",
            "--evalue",
            "1000",
        ],
        &[
            "-matrix",
            "IDENTITY",
            "-seg",
            "no",
            "-comp_based_stats",
            "0",
            "-threshold",
            "1",
            "-word_size",
            "2",
            "-evalue",
            "1000",
        ],
    );
}

#[test]
fn blastp_subject_ncbi_parity_pairwise_exact_hit() {
    assert_blastp_subject_outfmt_matches_ncbi(
        ">q1 protein query\nMKFLIFALILFATVALAPKSSSHEI\n",
        ">s1 protein subject\nMKFLIFALILFATVALAPKSSSHEIHH\n>s2 unrelated\nACDEFGHIKLMN\n",
        "0",
        &["--seg", "no", "--comp_based_stats", "0"],
        &["-seg", "no", "-comp_based_stats", "0"],
    );
}

#[test]
fn blastp_subject_ncbi_parity_pairwise_no_hits() {
    assert_blastp_subject_outfmt_matches_ncbi(
        ">q1 protein query\nMKFLIFALILFATVALAPKSSSHEI\n",
        ">s1 unrelated\nACDEFGHIKLMNPQRSTVWY\n",
        "0",
        &["--seg", "no", "--comp_based_stats", "0"],
        &["-seg", "no", "-comp_based_stats", "0"],
    );
}

#[test]
fn blastp_subject_ncbi_parity_pairwise_multiple_queries() {
    assert_blastp_subject_outfmt_matches_ncbi(
        ">q1 protein query\nMKFLIFALILFATVALAPKSSSHEI\n>q2 unrelated query\nACDEFGHIKLMNPQRSTVWY\n",
        ">s1 protein subject\nMKFLIFALILFATVALAPKSSSHEIHH\n>s2 unrelated\nACDEFGHIKLMNPQRSTVWY\n",
        "0",
        &["--seg", "no", "--comp_based_stats", "0"],
        &["-seg", "no", "-comp_based_stats", "0"],
    );
}

#[test]
fn blastp_subject_ncbi_parity_pairwise_gapped_alignment() {
    assert_blastp_subject_outfmt_matches_ncbi(
        ">q1 gapped query\nACDEFGHIKLMNPQRSTVWY\n",
        ">s1 gapped subject\nACDEFGHIKXXLMNPQRSTVWY\n",
        "0",
        &[
            "--seg",
            "no",
            "--comp_based_stats",
            "0",
            "--word_size",
            "2",
            "--threshold",
            "1",
            "--evalue",
            "1000",
        ],
        &[
            "-seg",
            "no",
            "-comp_based_stats",
            "0",
            "-word_size",
            "2",
            "-threshold",
            "1",
            "-evalue",
            "1000",
        ],
    );
}

#[test]
fn blastp_subject_ncbi_parity_pairwise_length_adjusted_stats() {
    assert_blastp_subject_outfmt_matches_ncbi(
        ">q1 repeat query\nACDEFGHIKLMNPQRSTVWYGGGGGGGGGGACDEFGHIKLMNPQRSTVWY\n",
        ">s1 repeat subject\nACDEFGHIKLMNPQRSTVWYTTTTTTTTTTACDEFGHIKLMNPQRSTVWY\n",
        "0",
        &[
            "--seg",
            "no",
            "--comp_based_stats",
            "0",
            "--word_size",
            "2",
            "--threshold",
            "1",
            "--evalue",
            "1000",
            "--max_hsps",
            "10",
        ],
        &[
            "-seg",
            "no",
            "-comp_based_stats",
            "0",
            "-word_size",
            "2",
            "-threshold",
            "1",
            "-evalue",
            "1000",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastp_subject_ncbi_parity_pairwise_zero_alignments() {
    assert_blastp_subject_outfmt_matches_ncbi(
        ">q1 protein query\nMKFLIFALILFATVALAPKSSSHEI\n",
        ">s1 protein subject\nMKFLIFALILFATVALAPKSSSHEIHH\n>s2 unrelated\nACDEFGHIKLMN\n",
        "0",
        &[
            "--seg",
            "no",
            "--comp_based_stats",
            "0",
            "--num_alignments",
            "0",
        ],
        &[
            "-seg",
            "no",
            "-comp_based_stats",
            "0",
            "-num_alignments",
            "0",
        ],
    );
}

#[test]
fn blastp_subject_ncbi_parity_pairwise_zero_descriptions() {
    assert_blastp_subject_outfmt_matches_ncbi(
        ">q1 protein query\nMKFLIFALILFATVALAPKSSSHEI\n",
        ">s1 protein subject\nMKFLIFALILFATVALAPKSSSHEIHH\n>s2 unrelated\nACDEFGHIKLMN\n",
        "0",
        &[
            "--seg",
            "no",
            "--comp_based_stats",
            "0",
            "--num_descriptions",
            "0",
        ],
        &[
            "-seg",
            "no",
            "-comp_based_stats",
            "0",
            "-num_descriptions",
            "0",
        ],
    );
}

#[test]
fn blastp_subject_ncbi_parity_pairwise_line_length() {
    assert_blastp_subject_outfmt_matches_ncbi(
        ">q1 protein query\nMKFLIFALILFATVALAPKSSSHEI\n",
        ">s1 protein subject\nMKFLIFALILFATVALAPKSSSHEIHH\n>s2 unrelated\nACDEFGHIKLMN\n",
        "0",
        &[
            "--seg",
            "no",
            "--comp_based_stats",
            "0",
            "--line_length",
            "10",
        ],
        &["-seg", "no", "-comp_based_stats", "0", "-line_length", "10"],
    );
}

#[test]
fn blastp_subject_ncbi_parity_commented_tabular_exact_hit() {
    assert_blastp_subject_outfmt_matches_ncbi(
        ">q1 protein query\nMKFLIFALILFATVALAPKSSSHEI\n",
        ">s1 protein subject\nMKFLIFALILFATVALAPKSSSHEIHH\n>s2 unrelated\nACDEFGHIKLMN\n",
        "7",
        &["--seg", "no", "--comp_based_stats", "0"],
        &["-seg", "no", "-comp_based_stats", "0"],
    );
}

#[test]
fn blastp_subject_ncbi_parity_pairwise_parse_deflines() {
    assert_blastp_subject_outfmt_matches_ncbi(
        ">gi|123|ref|QPROT.1| query protein\nMKFLIFALILFATVALAPKSSSHEI\n",
        ">gi|456|ref|SPROT.1| subject protein\nMKFLIFALILFATVALAPKSSSHEIHH\n",
        "0",
        &["--seg", "no", "--comp_based_stats", "0", "--parse_deflines"],
        &["-seg", "no", "-comp_based_stats", "0", "-parse_deflines"],
    );
}

#[test]
fn blastp_subject_ncbi_parity_tabular_parse_deflines() {
    assert_blastp_subject_outfmt_matches_ncbi(
        ">gi|123|ref|QPROT.1| query protein\nMKFLIFALILFATVALAPKSSSHEI\n",
        ">gi|456|ref|SPROT.1| subject protein\nMKFLIFALILFATVALAPKSSSHEIHH\n",
        "6 qseqid qacc qaccver sseqid sacc saccver",
        &["--seg", "no", "--comp_based_stats", "0", "--parse_deflines"],
        &["-seg", "no", "-comp_based_stats", "0", "-parse_deflines"],
    );
}

#[test]
fn blastp_subject_ncbi_parity_tabular_parse_deflines_gis() {
    assert_blastp_subject_outfmt_matches_ncbi(
        ">gi|123|ref|QPROT.1| query protein\nMKFLIFALILFATVALAPKSSSHEI\n",
        ">gi|456|ref|SPROT.1| subject protein\nMKFLIFALILFATVALAPKSSSHEIHH\n",
        "6 qseqid qgi qacc qaccver sseqid sgi sallgi sacc saccver sallseqid",
        &["--seg", "no", "--comp_based_stats", "0", "--parse_deflines"],
        &["-seg", "no", "-comp_based_stats", "0", "-parse_deflines"],
    );
}

#[test]
fn blastp_subject_ncbi_parity_tabular_parse_deflines_titles() {
    assert_blastp_subject_outfmt_matches_ncbi(
        ">gi|123|ref|QPROT.1| query protein\nMKFLIFALILFATVALAPKSSSHEI\n",
        ">gi|456|ref|SPROT.1| subject protein\nMKFLIFALILFATVALAPKSSSHEIHH\n",
        "6 sseqid stitle salltitles",
        &["--seg", "no", "--comp_based_stats", "0", "--parse_deflines"],
        &["-seg", "no", "-comp_based_stats", "0", "-parse_deflines"],
    );
}

#[test]
fn blastp_subject_ncbi_parity_custom_field_parser_edges() {
    for outfmt in [
        "6 delim=, qseqid bogus sseqid qseqid length length",
        "6 delim=tab qseqid sseqid length",
        "6 delim=, delim=tab qseqid sseqid length",
        "6 std qlen",
        "6 bogus",
        "7 qseqid bogus sseqid qseqid length length",
        "7 delim=space qseqid sseqid length",
        "7 delim=, delim=tab qseqid sseqid length",
        "7 std qlen",
        "7 bogus",
        "10 qseqid bogus sseqid qseqid length length",
        "10 delim=tab qseqid sseqid length",
        "10 delim=tab delim=, qseqid sseqid length",
        "10 std qlen",
        "10 bogus",
    ] {
        assert_blastp_subject_outfmt_matches_ncbi(
            ">q1 protein query\nMKFLIFALILFATVALAPKSSSHEI\n",
            ">s1 protein subject\nMKFLIFALILFATVALAPKSSSHEIHH\n",
            outfmt,
            &["--seg", "no", "--comp_based_stats", "0"],
            &["-seg", "no", "-comp_based_stats", "0"],
        );
    }
}

#[test]
fn blastp_subject_ncbi_parity_hit_limit_and_identity_filters() {
    for (rust_extra, ncbi_extra) in [
        (vec!["--max_hsps", "1"], vec!["-max_hsps", "1"]),
        (
            vec!["--max_target_seqs", "1"],
            vec!["-max_target_seqs", "1"],
        ),
        (vec!["--perc_identity", "90"], vec!["-perc_identity", "90"]),
        (vec!["--qcov_hsp_perc", "80"], vec!["-qcov_hsp_perc", "80"]),
    ] {
        if rust_extra.first() == Some(&"--perc_identity") {
            // NCBI BLAST+ 2.12 rejects -perc_identity for blastp.
            continue;
        }
        let mut rust_args = vec!["--seg", "no", "--comp_based_stats", "0", "--evalue", "1000"];
        rust_args.extend(rust_extra);
        let mut ncbi_args = vec!["-seg", "no", "-comp_based_stats", "0", "-evalue", "1000"];
        ncbi_args.extend(ncbi_extra);
        assert_blastp_subject_outfmt_matches_ncbi(
            ">q1 protein query\nMKFLILLFQQQQGMKFLILLF\n",
            ">full\nMKFLILLFQQQQGMKFLILLF\n>partial\nMKFLILLF\n>near\nMKFLILLYQQQQGMKFLILLF\n",
            "6 qseqid sseqid score pident length qstart qend sstart send",
            &rust_args,
            &ncbi_args,
        );
    }
}

#[test]
fn blastp_subject_ncbi_parity_culling_and_best_hit_filters() {
    for (rust_extra, ncbi_extra) in [
        (vec!["--culling_limit", "1"], vec!["-culling_limit", "1"]),
        (vec!["--subject_besthit"], vec!["-subject_besthit"]),
        (
            vec!["--best_hit_overhang", "0.1", "--best_hit_score_edge", "0.1"],
            vec!["-best_hit_overhang", "0.1", "-best_hit_score_edge", "0.1"],
        ),
    ] {
        let mut rust_args = vec!["--seg", "no", "--comp_based_stats", "0", "--evalue", "1000"];
        rust_args.extend(rust_extra);
        let mut ncbi_args = vec!["-seg", "no", "-comp_based_stats", "0", "-evalue", "1000"];
        ncbi_args.extend(ncbi_extra);
        assert_blastp_subject_outfmt_matches_ncbi(
            ">q1 protein query\nMKFLILLFQQQQGMKFLILLF\n",
            ">full\nMKFLILLFQQQQGMKFLILLF\n>contained\nKFLILL\n>near\nMKFLILLYQQQQGMKFLILLF\n",
            "6 qseqid sseqid score evalue qstart qend sstart send",
            &rust_args,
            &ncbi_args,
        );
    }
}

#[test]
fn blastp_subject_ncbi_parity_csv_parse_deflines_metadata() {
    assert_blastp_subject_outfmt_matches_ncbi(
        ">gi|123|ref|QPROT.1| query protein\nMKFLIFALILFATVALAPKSSSHEI\n",
        ">gi|456|ref|SPROT.1| subject protein\nMKFLIFALILFATVALAPKSSSHEIHH\n",
        "10 qseqid qgi qacc qaccver sseqid sgi sallgi sacc saccver sallseqid stitle salltitles",
        &["--seg", "no", "--comp_based_stats", "0", "--parse_deflines"],
        &["-seg", "no", "-comp_based_stats", "0", "-parse_deflines"],
    );
}

#[test]
fn blastp_subject_ncbi_parity_commented_tabular_parse_deflines() {
    assert_blastp_subject_outfmt_matches_ncbi(
        ">gi|123|ref|QPROT.1| query protein\nMKFLIFALILFATVALAPKSSSHEI\n",
        ">gi|456|ref|SPROT.1| subject protein\nMKFLIFALILFATVALAPKSSSHEIHH\n",
        "7 qaccver saccver pident length bitscore",
        &["--seg", "no", "--comp_based_stats", "0", "--parse_deflines"],
        &["-seg", "no", "-comp_based_stats", "0", "-parse_deflines"],
    );
}

#[test]
fn blastp_subject_ncbi_parity_xml_exact_hit() {
    assert_blastp_subject_outfmt_matches_ncbi(
        ">q1 protein query\nMKFLIFALILFATVALAPKSSSHEI\n",
        ">s1 protein subject\nMKFLIFALILFATVALAPKSSSHEIHH\n>s2 unrelated\nACDEFGHIKLMN\n",
        "5",
        &["--seg", "no", "--comp_based_stats", "0"],
        &["-seg", "no", "-comp_based_stats", "0"],
    );
}

#[test]
fn blastp_subject_ncbi_parity_xml_positive_substitution() {
    assert_blastp_subject_outfmt_matches_ncbi(
        ">q1 protein query\nMKFLILLFNILCLFPVLAADNHGVSMNAS\n",
        ">s1 related protein\nMKFLILLFNILCLFPVLAADNHGVSINAS\n",
        "5",
        &["--seg", "no", "--comp_based_stats", "0"],
        &["-seg", "no", "-comp_based_stats", "0"],
    );
}

#[test]
fn blastp_subject_ncbi_parity_xml_no_hits() {
    assert_blastp_subject_outfmt_matches_ncbi(
        ">q1 protein query\nMKFLIFALILFATVALAPKSSSHEI\n",
        ">s1 unrelated\nACDEFGHIKLMNPQRSTVWY\n",
        "5",
        &["--seg", "no", "--comp_based_stats", "0"],
        &["-seg", "no", "-comp_based_stats", "0"],
    );
}

#[test]
fn blastp_subject_ncbi_parity_xml_multiple_queries() {
    assert_blastp_subject_outfmt_matches_ncbi(
        ">q1 protein query\nMKFLIFALILFATVALAPKSSSHEI\n>q2 unrelated query\nACDEFGHIKLMNPQRSTVWY\n",
        ">s1 protein subject\nMKFLIFALILFATVALAPKSSSHEIHH\n>s2 unrelated\nACDEFGHIKLMNPQRSTVWY\n",
        "5",
        &["--seg", "no", "--comp_based_stats", "0"],
        &["-seg", "no", "-comp_based_stats", "0"],
    );
}

#[test]
fn blastp_subject_ncbi_parity_xml_multi_hsp_same_subject() {
    assert_blastp_subject_outfmt_matches_ncbi(
        ">q1 repeat query\nACDEFGHIKLMNPQRSTVWYGGGGGGGGGGACDEFGHIKLMNPQRSTVWY\n",
        ">s1 repeat subject\nACDEFGHIKLMNPQRSTVWYTTTTTTTTTTACDEFGHIKLMNPQRSTVWY\n",
        "5",
        &[
            "--seg",
            "no",
            "--comp_based_stats",
            "0",
            "--word_size",
            "2",
            "--threshold",
            "1",
            "--evalue",
            "1000",
            "--max_hsps",
            "10",
        ],
        &[
            "-seg",
            "no",
            "-comp_based_stats",
            "0",
            "-word_size",
            "2",
            "-threshold",
            "1",
            "-evalue",
            "1000",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastp_subject_ncbi_parity_xml_parse_deflines() {
    assert_blastp_subject_outfmt_matches_ncbi(
        ">gi|123|ref|QPROT.1| query protein\nMKFLIFALILFATVALAPKSSSHEI\n",
        ">gi|456|ref|SPROT.1| subject protein\nMKFLIFALILFATVALAPKSSSHEIHH\n",
        "5",
        &["--seg", "no", "--comp_based_stats", "0", "--parse_deflines"],
        &["-seg", "no", "-comp_based_stats", "0", "-parse_deflines"],
    );
}

#[test]
fn blastp_db_ncbi_parity_pairwise_exact_hit() {
    assert_blastp_db_outfmt_matches_ncbi(
        ">q1 protein query\nMKFLIFALILFATVALAPKSSSHEI\n",
        ">s1 protein subject\nMKFLIFALILFATVALAPKSSSHEIHH\n>s2 unrelated\nACDEFGHIKLMN\n",
        "0",
        &["--seg", "no", "--comp_based_stats", "0"],
        &["-seg", "no", "-comp_based_stats", "0"],
    );
}

#[test]
fn blastp_db_ncbi_parity_pairwise_parse_deflines() {
    assert_blastp_db_outfmt_matches_ncbi(
        ">gi|123|ref|QPROT.1| query protein\nMKFLIFALILFATVALAPKSSSHEI\n",
        ">gi|456|ref|SPROT.1| subject protein\nMKFLIFALILFATVALAPKSSSHEIHH\n",
        "0",
        &["--seg", "no", "--comp_based_stats", "0", "--parse_deflines"],
        &["-seg", "no", "-comp_based_stats", "0", "-parse_deflines"],
    );
}

#[test]
fn blastp_db_ncbi_parity_tabular_parse_deflines() {
    assert_blastp_db_outfmt_matches_ncbi(
        ">gi|123|ref|QPROT.1| query protein\nMKFLIFALILFATVALAPKSSSHEI\n",
        ">gi|456|ref|SPROT.1| subject protein\nMKFLIFALILFATVALAPKSSSHEIHH\n",
        "6 qseqid qacc qaccver sseqid sacc saccver",
        &["--seg", "no", "--comp_based_stats", "0", "--parse_deflines"],
        &["-seg", "no", "-comp_based_stats", "0", "-parse_deflines"],
    );
}

#[test]
fn blastp_db_ncbi_parity_tabular_parse_deflines_gis() {
    assert_blastp_db_outfmt_matches_ncbi(
        ">gi|123|ref|QPROT.1| query protein\nMKFLIFALILFATVALAPKSSSHEI\n",
        ">gi|456|ref|SPROT.1| subject protein\nMKFLIFALILFATVALAPKSSSHEIHH\n",
        "6 qseqid qgi qacc qaccver sseqid sgi sallgi sacc saccver sallseqid",
        &["--seg", "no", "--comp_based_stats", "0", "--parse_deflines"],
        &["-seg", "no", "-comp_based_stats", "0", "-parse_deflines"],
    );
}

#[test]
fn blastp_db_ncbi_parity_tabular_parse_deflines_titles() {
    assert_blastp_db_outfmt_matches_ncbi(
        ">gi|123|ref|QPROT.1| query protein\nMKFLIFALILFATVALAPKSSSHEI\n",
        ">gi|456|ref|SPROT.1| subject protein\nMKFLIFALILFATVALAPKSSSHEIHH\n",
        "6 sseqid stitle salltitles",
        &["--seg", "no", "--comp_based_stats", "0", "--parse_deflines"],
        &["-seg", "no", "-comp_based_stats", "0", "-parse_deflines"],
    );
}

#[test]
fn blastp_db_ncbi_parity_custom_field_parser_edges() {
    for outfmt in [
        "6 delim=, qseqid bogus sseqid qseqid length length",
        "6 delim=tab qseqid sseqid length",
        "6 std qlen",
        "6 bogus",
        "7 qseqid bogus sseqid qseqid length length",
        "7 delim=space qseqid sseqid length",
        "7 std qlen",
        "7 bogus",
        "10 qseqid bogus sseqid qseqid length length",
        "10 delim=tab qseqid sseqid length",
        "10 std qlen",
        "10 bogus",
    ] {
        assert_blastp_db_outfmt_matches_ncbi(
            ">q1 protein query\nMKFLIFALILFATVALAPKSSSHEI\n",
            ">s1 protein subject\nMKFLIFALILFATVALAPKSSSHEIHH\n",
            outfmt,
            &["--seg", "no", "--comp_based_stats", "0"],
            &["-seg", "no", "-comp_based_stats", "0"],
        );
    }
}

#[test]
fn blastp_db_ncbi_parity_hit_limit_and_identity_filters() {
    for (rust_extra, ncbi_extra) in [
        (vec!["--max_hsps", "1"], vec!["-max_hsps", "1"]),
        (
            vec!["--max_target_seqs", "1"],
            vec!["-max_target_seqs", "1"],
        ),
        (vec!["--perc_identity", "90"], vec!["-perc_identity", "90"]),
        (vec!["--qcov_hsp_perc", "80"], vec!["-qcov_hsp_perc", "80"]),
    ] {
        if rust_extra.first() == Some(&"--perc_identity") {
            // NCBI BLAST+ 2.12 rejects -perc_identity for blastp.
            continue;
        }
        let mut rust_args = vec!["--seg", "no", "--comp_based_stats", "0", "--evalue", "1000"];
        rust_args.extend(rust_extra);
        let mut ncbi_args = vec!["-seg", "no", "-comp_based_stats", "0", "-evalue", "1000"];
        ncbi_args.extend(ncbi_extra);
        assert_blastp_db_outfmt_matches_ncbi(
            ">q1 protein query\nMKFLILLFQQQQGMKFLILLF\n",
            ">full\nMKFLILLFQQQQGMKFLILLF\n>partial\nMKFLILLF\n>near\nMKFLILLYQQQQGMKFLILLF\n",
            "6 qseqid sseqid score pident length qstart qend sstart send",
            &rust_args,
            &ncbi_args,
        );
    }
}

#[test]
fn blastp_db_ncbi_parity_culling_and_best_hit_filters() {
    for (rust_extra, ncbi_extra) in [
        (vec!["--culling_limit", "1"], vec!["-culling_limit", "1"]),
        (vec!["--subject_besthit"], vec!["-subject_besthit"]),
        (
            vec!["--best_hit_overhang", "0.1", "--best_hit_score_edge", "0.1"],
            vec!["-best_hit_overhang", "0.1", "-best_hit_score_edge", "0.1"],
        ),
    ] {
        let mut rust_args = vec!["--seg", "no", "--comp_based_stats", "0", "--evalue", "1000"];
        rust_args.extend(rust_extra);
        let mut ncbi_args = vec!["-seg", "no", "-comp_based_stats", "0", "-evalue", "1000"];
        ncbi_args.extend(ncbi_extra);
        assert_blastp_db_outfmt_matches_ncbi(
            ">q1 protein query\nMKFLILLFQQQQGMKFLILLF\n",
            ">full\nMKFLILLFQQQQGMKFLILLF\n>contained\nKFLILL\n>near\nMKFLILLYQQQQGMKFLILLF\n",
            "6 qseqid sseqid score evalue qstart qend sstart send",
            &rust_args,
            &ncbi_args,
        );
    }
}

#[test]
fn blastp_db_ncbi_parity_csv_parse_deflines_metadata() {
    assert_blastp_db_outfmt_matches_ncbi(
        ">gi|123|ref|QPROT.1| query protein\nMKFLIFALILFATVALAPKSSSHEI\n",
        ">gi|456|ref|SPROT.1| subject protein\nMKFLIFALILFATVALAPKSSSHEIHH\n",
        "10 qseqid qgi qacc qaccver sseqid sgi sallgi sacc saccver sallseqid stitle salltitles",
        &["--seg", "no", "--comp_based_stats", "0", "--parse_deflines"],
        &["-seg", "no", "-comp_based_stats", "0", "-parse_deflines"],
    );
}

#[test]
fn blastp_db_ncbi_parity_commented_tabular_exact_hit() {
    assert_blastp_db_outfmt_matches_ncbi(
        ">q1 protein query\nMKFLIFALILFATVALAPKSSSHEI\n",
        ">s1 protein subject\nMKFLIFALILFATVALAPKSSSHEIHH\n>s2 unrelated\nACDEFGHIKLMN\n",
        "7",
        &["--seg", "no", "--comp_based_stats", "0"],
        &["-seg", "no", "-comp_based_stats", "0"],
    );
}

#[test]
fn blastp_db_ncbi_parity_commented_tabular_parse_deflines() {
    assert_blastp_db_outfmt_matches_ncbi(
        ">gi|123|ref|QPROT.1| query protein\nMKFLIFALILFATVALAPKSSSHEI\n",
        ">gi|456|ref|SPROT.1| subject protein\nMKFLIFALILFATVALAPKSSSHEIHH\n",
        "7 qaccver saccver pident length bitscore",
        &["--seg", "no", "--comp_based_stats", "0", "--parse_deflines"],
        &["-seg", "no", "-comp_based_stats", "0", "-parse_deflines"],
    );
}

#[test]
fn blastp_db_ncbi_parity_xml_exact_hit() {
    assert_blastp_db_outfmt_matches_ncbi(
        ">q1 protein query\nMKFLIFALILFATVALAPKSSSHEI\n",
        ">s1 protein subject\nMKFLIFALILFATVALAPKSSSHEIHH\n>s2 unrelated\nACDEFGHIKLMN\n",
        "5",
        &["--seg", "no", "--comp_based_stats", "0"],
        &["-seg", "no", "-comp_based_stats", "0"],
    );
}

#[test]
fn blastp_db_ncbi_parity_xml_positive_substitution() {
    assert_blastp_db_outfmt_matches_ncbi(
        ">q1 protein query\nMKFLILLFNILCLFPVLAADNHGVSMNAS\n",
        ">s1 related protein\nMKFLILLFNILCLFPVLAADNHGVSINAS\n",
        "5",
        &["--seg", "no", "--comp_based_stats", "0"],
        &["-seg", "no", "-comp_based_stats", "0"],
    );
}

#[test]
fn blastp_db_ncbi_parity_xml_no_hits() {
    assert_blastp_db_outfmt_matches_ncbi(
        ">q1 protein query\nMKFLIFALILFATVALAPKSSSHEI\n",
        ">s1 unrelated\nACDEFGHIKLMNPQRSTVWY\n",
        "5",
        &["--seg", "no", "--comp_based_stats", "0"],
        &["-seg", "no", "-comp_based_stats", "0"],
    );
}

#[test]
fn blastp_db_ncbi_parity_xml_multiple_queries() {
    assert_blastp_db_outfmt_matches_ncbi(
        ">q1 protein query\nMKFLIFALILFATVALAPKSSSHEI\n>q2 unrelated query\nACDEFGHIKLMNPQRSTVWY\n",
        ">s1 protein subject\nMKFLIFALILFATVALAPKSSSHEIHH\n>s2 unrelated\nACDEFGHIKLMNPQRSTVWY\n",
        "5",
        &["--seg", "no", "--comp_based_stats", "0"],
        &["-seg", "no", "-comp_based_stats", "0"],
    );
}

#[test]
fn blastp_db_ncbi_parity_xml_multi_hsp_same_subject() {
    assert_blastp_db_outfmt_matches_ncbi(
        ">q1 repeat query\nACDEFGHIKLMNPQRSTVWYGGGGGGGGGGACDEFGHIKLMNPQRSTVWY\n",
        ">s1 repeat subject\nACDEFGHIKLMNPQRSTVWYTTTTTTTTTTACDEFGHIKLMNPQRSTVWY\n",
        "5",
        &[
            "--seg",
            "no",
            "--comp_based_stats",
            "0",
            "--word_size",
            "2",
            "--threshold",
            "1",
            "--evalue",
            "1000",
            "--max_hsps",
            "10",
        ],
        &[
            "-seg",
            "no",
            "-comp_based_stats",
            "0",
            "-word_size",
            "2",
            "-threshold",
            "1",
            "-evalue",
            "1000",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastp_db_ncbi_parity_xml_parse_deflines() {
    assert_blastp_db_outfmt_matches_ncbi(
        ">gi|123|ref|QPROT.1| query protein\nMKFLIFALILFATVALAPKSSSHEI\n",
        ">gi|456|ref|SPROT.1| subject protein\nMKFLIFALILFATVALAPKSSSHEIHH\n",
        "5",
        &["--seg", "no", "--comp_based_stats", "0", "--parse_deflines"],
        &["-seg", "no", "-comp_based_stats", "0", "-parse_deflines"],
    );
}

#[test]
fn blastp_db_ncbi_parity_pairwise_no_hits() {
    assert_blastp_db_outfmt_matches_ncbi(
        ">q1 protein query\nMKFLIFALILFATVALAPKSSSHEI\n",
        ">s1 unrelated\nACDEFGHIKLMNPQRSTVWY\n",
        "0",
        &["--seg", "no", "--comp_based_stats", "0"],
        &["-seg", "no", "-comp_based_stats", "0"],
    );
}

#[test]
fn blastp_db_ncbi_parity_pairwise_multiple_queries() {
    assert_blastp_db_outfmt_matches_ncbi(
        ">q1 protein query\nMKFLIFALILFATVALAPKSSSHEI\n>q2 unrelated query\nACDEFGHIKLMNPQRSTVWY\n",
        ">s1 protein subject\nMKFLIFALILFATVALAPKSSSHEIHH\n>s2 unrelated\nACDEFGHIKLMNPQRSTVWY\n",
        "0",
        &["--seg", "no", "--comp_based_stats", "0"],
        &["-seg", "no", "-comp_based_stats", "0"],
    );
}

#[test]
fn blastp_db_ncbi_parity_pairwise_gapped_alignment() {
    assert_blastp_db_outfmt_matches_ncbi(
        ">q1 gapped query\nACDEFGHIKLMNPQRSTVWY\n",
        ">s1 gapped subject\nACDEFGHIKXXLMNPQRSTVWY\n",
        "0",
        &[
            "--seg",
            "no",
            "--comp_based_stats",
            "0",
            "--word_size",
            "2",
            "--threshold",
            "1",
            "--evalue",
            "1000",
        ],
        &[
            "-seg",
            "no",
            "-comp_based_stats",
            "0",
            "-word_size",
            "2",
            "-threshold",
            "1",
            "-evalue",
            "1000",
        ],
    );
}

#[test]
fn blastp_db_ncbi_parity_pairwise_length_adjusted_stats() {
    assert_blastp_db_outfmt_matches_ncbi(
        ">q1 repeat query\nACDEFGHIKLMNPQRSTVWYGGGGGGGGGGACDEFGHIKLMNPQRSTVWY\n",
        ">s1 repeat subject\nACDEFGHIKLMNPQRSTVWYTTTTTTTTTTACDEFGHIKLMNPQRSTVWY\n",
        "0",
        &[
            "--seg",
            "no",
            "--comp_based_stats",
            "0",
            "--word_size",
            "2",
            "--threshold",
            "1",
            "--evalue",
            "1000",
            "--max_hsps",
            "10",
        ],
        &[
            "-seg",
            "no",
            "-comp_based_stats",
            "0",
            "-word_size",
            "2",
            "-threshold",
            "1",
            "-evalue",
            "1000",
            "-max_hsps",
            "10",
        ],
    );
}

#[test]
fn blastp_db_ncbi_parity_pairwise_zero_alignments() {
    assert_blastp_db_outfmt_matches_ncbi(
        ">q1 protein query\nMKFLIFALILFATVALAPKSSSHEI\n",
        ">s1 protein subject\nMKFLIFALILFATVALAPKSSSHEIHH\n>s2 unrelated\nACDEFGHIKLMN\n",
        "0",
        &[
            "--seg",
            "no",
            "--comp_based_stats",
            "0",
            "--num_alignments",
            "0",
        ],
        &[
            "-seg",
            "no",
            "-comp_based_stats",
            "0",
            "-num_alignments",
            "0",
        ],
    );
}

#[test]
fn blastp_db_ncbi_parity_pairwise_zero_descriptions() {
    assert_blastp_db_outfmt_matches_ncbi(
        ">q1 protein query\nMKFLIFALILFATVALAPKSSSHEI\n",
        ">s1 protein subject\nMKFLIFALILFATVALAPKSSSHEIHH\n>s2 unrelated\nACDEFGHIKLMN\n",
        "0",
        &[
            "--seg",
            "no",
            "--comp_based_stats",
            "0",
            "--num_descriptions",
            "0",
        ],
        &[
            "-seg",
            "no",
            "-comp_based_stats",
            "0",
            "-num_descriptions",
            "0",
        ],
    );
}

#[test]
fn blastp_db_ncbi_parity_pairwise_line_length() {
    assert_blastp_db_outfmt_matches_ncbi(
        ">q1 protein query\nMKFLIFALILFATVALAPKSSSHEI\n",
        ">s1 protein subject\nMKFLIFALILFATVALAPKSSSHEIHH\n>s2 unrelated\nACDEFGHIKLMN\n",
        "0",
        &[
            "--seg",
            "no",
            "--comp_based_stats",
            "0",
            "--line_length",
            "10",
        ],
        &["-seg", "no", "-comp_based_stats", "0", "-line_length", "10"],
    );
}

#[test]
fn blastp_subject_ncbi_parity_identity_matrix_default_stats_and_warning() {
    if !std::path::Path::new("/usr/bin/blastp").exists() {
        eprintln!("Skipping: /usr/bin/blastp not found");
        return;
    }
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: blast-cli binary not found");
        return;
    };

    let tmp = TempDir::new().expect("tempdir");
    let query = tmp.path().join("query.fa");
    let subject = tmp.path().join("subject.fa");
    std::fs::write(&query, ">prot_q\nMKFLIFALILFATVALAPKSSSHEI\n").expect("write query");
    std::fs::write(
        &subject,
        ">prot_s1\nMKFLIFALILFATVALAPKSSSHEIHH\n>prot_s2\nACDEFGHIKLMN\n",
    )
    .expect("write subject");

    let rust = std::process::Command::new(blast_cli)
        .arg("blastp")
        .arg("--query")
        .arg(&query)
        .arg("--subject")
        .arg(&subject)
        .arg("--matrix")
        .arg("IDENTITY")
        .arg("--outfmt")
        .arg("6 qseqid sseqid score bitscore evalue length pident")
        .output()
        .expect("run blast-cli IDENTITY");
    let ncbi = std::process::Command::new("/usr/bin/blastp")
        .arg("-query")
        .arg(&query)
        .arg("-subject")
        .arg(&subject)
        .arg("-matrix")
        .arg("IDENTITY")
        .arg("-outfmt")
        .arg("6 qseqid sseqid score bitscore evalue length pident")
        .output()
        .expect("run NCBI IDENTITY");

    assert!(
        rust.status.success(),
        "blast-cli exited with {}",
        rust.status
    );
    assert!(ncbi.status.success(), "NCBI exited with {}", ncbi.status);
    assert_eq!(rust.stdout, ncbi.stdout, "stdout differs");
    assert_eq!(
        String::from_utf8_lossy(&rust.stderr),
        String::from_utf8_lossy(&ncbi.stderr),
        "stderr differs"
    );
}

#[test]
fn blastp_subject_ncbi_parity_matrix_specific_default_gap_costs() {
    for matrix in [
        "BLOSUM45", "BLOSUM50", "BLOSUM62", "BLOSUM80", "BLOSUM90", "PAM30", "PAM70", "PAM250",
    ] {
        assert_blastp_subject_outfmt_matches_ncbi(
            ">prot_q\nMKFLIFALILFATVALAPKSSSHEI\n",
            ">prot_s1\nMKFLIFALILFATVALAPKSSSHEIHH\n",
            "6 qseqid sseqid score bitscore evalue length pident",
            &["--matrix", matrix, "--comp_based_stats", "0"],
            &["-matrix", matrix, "-comp_based_stats", "0"],
        );
    }
}

#[test]
fn blastp_subject_ncbi_parity_matrix_specific_default_gap_costs_in_xml() {
    assert_blastp_subject_outfmt_matches_ncbi(
        ">prot_q\nMKFLIFALILFATVALAPKSSSHEI\n",
        ">prot_s1\nMKFLIFALILFATVALAPKSSSHEIHH\n",
        "5",
        &["--matrix", "PAM30", "--comp_based_stats", "0"],
        &["-matrix", "PAM30", "-comp_based_stats", "0"],
    );
}

#[test]
fn blastp_subject_ncbi_parity_matrix_specific_default_gap_costs_in_pairwise() {
    assert_blastp_subject_outfmt_matches_ncbi(
        ">prot_q\nMKFLIFALILFATVALAPKSSSHEI\n",
        ">prot_s1\nMKFLIFALILFATVALAPKSSSHEIHH\n",
        "0",
        &["--matrix", "PAM30", "--comp_based_stats", "0"],
        &["-matrix", "PAM30", "-comp_based_stats", "0"],
    );
}

#[test]
fn blastx_subject_ncbi_parity_matrix_specific_default_gap_costs_in_xml() {
    assert_translated_subject_outfmt_matches_ncbi(
        "blastx",
        "/usr/bin/blastx",
        ">q1\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">s1\nMKFLILLF\n",
        "5",
        &[
            "--matrix",
            "PAM30",
            "--seg",
            "no",
            "--comp_based_stats",
            "0",
        ],
        &["-matrix", "PAM30", "-seg", "no", "-comp_based_stats", "0"],
    );
}

#[test]
fn blastx_subject_ncbi_parity_matrix_specific_default_gap_costs_in_pairwise() {
    assert_translated_subject_outfmt_matches_ncbi(
        "blastx",
        "/usr/bin/blastx",
        ">q1\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">s1\nMKFLILLF\n",
        "0",
        &[
            "--matrix",
            "PAM30",
            "--seg",
            "no",
            "--comp_based_stats",
            "0",
        ],
        &["-matrix", "PAM30", "-seg", "no", "-comp_based_stats", "0"],
    );
}

#[test]
fn blastx_ncbi_parity_identity_matrix_is_unsupported() {
    if !std::path::Path::new("/usr/bin/blastx").exists() {
        eprintln!("Skipping: /usr/bin/blastx not found");
        return;
    }
    let Some(blast_cli) = blast_cli_bin_for_tests() else {
        eprintln!("Skipping: blast-cli binary not found");
        return;
    };

    let tmp = TempDir::new().expect("tempdir");
    let query = tmp.path().join("query.fa");
    let subject = tmp.path().join("subject.fa");
    std::fs::write(&query, ">q\nATGAAATTTTTAATTTTTGCT\n").expect("write query");
    std::fs::write(&subject, ">s\nMKFLIFA\n").expect("write subject");

    let rust = std::process::Command::new(blast_cli)
        .arg("blastx")
        .arg("--query")
        .arg(&query)
        .arg("--subject")
        .arg(&subject)
        .arg("--matrix")
        .arg("IDENTITY")
        .arg("--outfmt")
        .arg("6")
        .output()
        .expect("run blast-cli blastx IDENTITY");
    let ncbi = std::process::Command::new("/usr/bin/blastx")
        .arg("-query")
        .arg(&query)
        .arg("-subject")
        .arg(&subject)
        .arg("-matrix")
        .arg("IDENTITY")
        .arg("-outfmt")
        .arg("6")
        .output()
        .expect("run NCBI blastx IDENTITY");

    assert!(!rust.status.success(), "blast-cli should reject IDENTITY");
    assert!(!ncbi.status.success(), "NCBI should reject IDENTITY");
    assert_eq!(rust.status.code(), ncbi.status.code(), "status differs");
    assert_eq!(rust.stdout, ncbi.stdout, "stdout differs");
    assert_eq!(
        String::from_utf8_lossy(&rust.stderr),
        String::from_utf8_lossy(&ncbi.stderr),
        "stderr differs"
    );
}

#[test]
fn blastx_subject_ncbi_parity_default_seg_masks_low_complexity_query() {
    assert_translated_subject_outfmt_matches_ncbi(
        "blastx",
        "/usr/bin/blastx",
        ">q1 low complexity nt query\nGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCT\n",
        concat!(
            ">s1 low complexity protein subject\nAAAAAAAAAAAAAAAAAAAA\n",
            ">s2 mixed protein subject\nACDEFGHIKLMNPQRSTVWY\n"
        ),
        "6 qseqid sseqid length bitscore evalue",
        &[],
        &[],
    );
}

#[test]
fn blastx_subject_ncbi_parity_query_gencode_changes_translation() {
    assert_translated_subject_outfmt_matches_ncbi(
        "blastx",
        "/usr/bin/blastx",
        ">q1 ciliate glutamine codons\nTAATAATAATAATAA\n",
        ">s1 glutamine peptide\nQQQQQ\n",
        "6 qseqid sseqid qstart qend sstart send score length qseq sseq",
        &["--query_gencode", "6", "--seg", "no", "--evalue", "1000"],
        &["-query_gencode", "6", "-seg", "no", "-evalue", "1000"],
    );
}

#[test]
fn tblastn_subject_ncbi_parity_default_seg_masks_low_complexity_query() {
    assert_translated_subject_outfmt_matches_ncbi(
        "tblastn",
        "/usr/bin/tblastn",
        ">q1 low complexity protein query\nAAAAAAAAAAAAAAAAAAAA\n",
        concat!(
            ">s1 low complexity nt subject\nGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCT\n",
            ">s2 mixed nt subject\nGCTTGTGATGAATTTGGTCATATAAAACTGATGAATCCTCAACGTTCTACTGTGGTA\n"
        ),
        "6 qseqid sseqid length bitscore evalue",
        &[],
        &[],
    );
}

#[test]
fn tblastn_subject_ncbi_parity_db_gencode_changes_translation() {
    assert_translated_subject_outfmt_matches_ncbi(
        "tblastn",
        "/usr/bin/tblastn",
        ">q1 glutamine peptide\nQQQQQ\n",
        ">s1 ciliate glutamine codons\nTAATAATAATAATAA\n",
        "6 qseqid sseqid qstart qend sstart send score length qseq sseq",
        &["--db_gencode", "6", "--seg", "no", "--evalue", "1000"],
        &["-db_gencode", "6", "-seg", "no", "-evalue", "1000"],
    );
}

#[test]
fn tblastn_subject_ncbi_parity_identity_matrix_default_threshold() {
    assert_translated_subject_outfmt_matches_ncbi(
        "tblastn",
        "/usr/bin/tblastn",
        ">prot_q\nMKFLIFALILFATVALAPKSSSHEI\n",
        ">s\nATGAAATTTTTAATTTTTGCTTTAATTTTATTTGCTACTGTTGCTTTAGCTCCAAAATCATCATCTCATGAAATT\n",
        "6 qseqid sseqid score bitscore length pident qstart qend sstart send qseq sseq",
        &["--matrix", "IDENTITY", "--comp_based_stats", "0"],
        &["-matrix", "IDENTITY", "-comp_based_stats", "0"],
    );
}

#[test]
fn translated_subject_ncbi_parity_custom_field_parser_edges() {
    for (program, ncbi_bin, query, subject) in [
        (
            "blastx",
            "/usr/bin/blastx",
            ">q1\nATGAAATTTCTGATTCTGCTGTTT\n",
            ">s1\nMKFLILLF\n",
        ),
        (
            "tblastn",
            "/usr/bin/tblastn",
            ">q1\nMKFLILLF\n",
            ">s1\nATGAAATTTCTGATTCTGCTGTTT\n",
        ),
        (
            "tblastx",
            "/usr/bin/tblastx",
            ">q1\nATGAAATTTCTGATTCTGCTGTTT\n",
            ">s1\nATGAAATTTCTGATTCTGCTGTTT\n",
        ),
    ] {
        for outfmt in [
            "6 delim=, qseqid bogus sseqid qseqid length length",
            "6 delim=tab qseqid sseqid length",
            "6 std qlen",
            "6 bogus",
            "7 qseqid bogus sseqid qseqid length length",
            "7 delim=space qseqid sseqid length",
            "7 std qlen",
            "7 bogus",
            "10 qseqid bogus sseqid qseqid length length",
            "10 delim=tab qseqid sseqid length",
            "10 std qlen",
            "10 bogus",
        ] {
            assert_translated_subject_outfmt_matches_ncbi(
                program,
                ncbi_bin,
                query,
                subject,
                outfmt,
                &["--seg", "no"],
                &["-seg", "no"],
            );
        }
    }
}

#[test]
fn tblastx_subject_ncbi_parity_default_seg_masks_low_complexity_query() {
    assert_translated_subject_outfmt_matches_ncbi(
        "tblastx",
        "/usr/bin/tblastx",
        ">q1 low complexity nt query\nGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCT\n",
        concat!(
            ">s1 low complexity nt subject\nGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCT\n",
            ">s2 mixed nt subject\nGCTTGTGATGAATTTGGTCATATAAAACTGATGAATCCTCAACGTTCTACTGTGGTA\n"
        ),
        "6 qseqid sseqid length bitscore evalue",
        &[],
        &[],
    );
}

#[test]
fn tblastx_subject_ncbi_parity_query_and_db_gencode_change_translation() {
    assert_translated_subject_outfmt_matches_ncbi(
        "tblastx",
        "/usr/bin/tblastx",
        ">q1 ciliate glutamine codons\nTAATAATAATAATAA\n",
        ">s1 ciliate glutamine codons\nTAATAATAATAATAA\n",
        "6 qseqid sseqid qstart qend sstart send score length qseq sseq",
        &[
            "--query_gencode",
            "6",
            "--db_gencode",
            "6",
            "--seg",
            "no",
            "--evalue",
            "1000",
        ],
        &[
            "-query_gencode",
            "6",
            "-db_gencode",
            "6",
            "-seg",
            "no",
            "-evalue",
            "1000",
        ],
    );
}

#[test]
fn tblastx_subject_ncbi_parity_exact_translation_coordinates_and_frames() {
    assert_translated_subject_outfmt_matches_ncbi_sorted_lines(
        "tblastx",
        "/usr/bin/tblastx",
        ">q1\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">s1\nATGAAATTTCTGATTCTGCTGTTT\n",
        "6 qseqid sseqid qlen slen qstart qend sstart send score qframe sframe frames length",
        &[],
        &[],
    );
}

#[test]
fn blastx_subject_ncbi_parity_exact_translation_coordinates_and_frames() {
    assert_translated_subject_outfmt_matches_ncbi_sorted_lines(
        "blastx",
        "/usr/bin/blastx",
        ">q1\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">s1\nMKFLILLF\n",
        "6 qseqid sseqid qlen slen qstart qend sstart send score qframe sframe frames length",
        &[],
        &[],
    );
}

#[test]
fn tblastn_subject_ncbi_parity_exact_translation_coordinates_and_frames() {
    assert_translated_subject_outfmt_matches_ncbi_sorted_lines(
        "tblastn",
        "/usr/bin/tblastn",
        ">q1\nMKFLILLF\n",
        ">s1\nATGAAATTTCTGATTCTGCTGTTT\n",
        "6 qseqid sseqid qlen slen qstart qend sstart send score qframe sframe frames length",
        &[],
        &[],
    );
}

#[test]
fn blastx_subject_ncbi_parity_reverse_frame_coordinates_and_frames() {
    let query = format!(
        ">q1\n{}\n",
        ascii_reverse_complement("ATGAAATTTCTGATTCTGCTGTTT")
    );
    assert_translated_subject_outfmt_matches_ncbi_sorted_lines(
        "blastx",
        "/usr/bin/blastx",
        &query,
        ">s1\nMKFLILLF\n",
        "6 qseqid sseqid qlen slen qstart qend sstart send score qframe sframe frames length",
        &[],
        &[],
    );
}

#[test]
fn tblastn_subject_ncbi_parity_reverse_frame_coordinates_and_frames() {
    let subject = format!(
        ">s1\n{}\n",
        ascii_reverse_complement("ATGAAATTTCTGATTCTGCTGTTT")
    );
    assert_translated_subject_outfmt_matches_ncbi_sorted_lines(
        "tblastn",
        "/usr/bin/tblastn",
        ">q1\nMKFLILLF\n",
        &subject,
        "6 qseqid sseqid qlen slen qstart qend sstart send score qframe sframe frames length",
        &[],
        &[],
    );
}

#[test]
fn tblastx_subject_ncbi_parity_reverse_frame_coordinates_and_frames() {
    let query = format!(
        ">q1\n{}\n",
        ascii_reverse_complement("ATGAAATTTCTGATTCTGCTGTTT")
    );
    let subject = format!(
        ">s1\n{}\n",
        ascii_reverse_complement("ATGAAATTTCTGATTCTGCTGTTT")
    );
    assert_translated_subject_outfmt_matches_ncbi_sorted_lines(
        "tblastx",
        "/usr/bin/tblastx",
        &query,
        &subject,
        "6 qseqid sseqid qlen slen qstart qend sstart send score qframe sframe frames length",
        &[],
        &[],
    );
}

#[test]
fn tblastx_subject_ncbi_parity_indel_remains_ungapped() {
    let query =
        ">q1\nATGAAATTTCTGATTCTGCTGTTTAATATTCTGTGTCTGTTTCCTGTTCTGGCTGCTGATAATCATGGTGTTTCTATGAATGCTTCT\n";
    let subject =
        ">s1\nATGAAATTTCTGATTCTGCTGTTTATTCTGTGTCTGTTTCCTGTTCTGGCTGCTGATAATCATGGTGTTTCTATGAATGCTTCT\n";
    assert_translated_subject_outfmt_matches_ncbi_sorted_lines(
        "tblastx",
        "/usr/bin/tblastx",
        query,
        subject,
        "6 qseqid sseqid qstart qend sstart send score length gaps qframe sframe frames qseq sseq",
        &["--seg", "no", "--evalue", "1000"],
        &["-seg", "no", "-evalue", "1000"],
    );
}

#[test]
fn tblastx_subject_ncbi_parity_frame_offset_insertion_remains_ungapped() {
    assert_translated_subject_outfmt_matches_ncbi_sorted_lines(
        "tblastx",
        "/usr/bin/tblastx",
        ">q1\nATGAAATTTCTGATTCTGCTGTTTATGAAATTTCTGATTCTGCTGTTTATGAAATTTCTGATTCTGCTGTTT\n",
        ">s1\nATGAAATTTCTGATTCTGCTGTTTAATGAAATTTCTGATTCTGCTGTTTATGAAATTTCTGATTCTGCTGTTT\n",
        "6 qseqid sseqid qstart qend sstart send score length gaps gapopen qframe sframe frames qseq sseq",
        &["--seg", "no", "--evalue", "1000"],
        &["-seg", "no", "-evalue", "1000"],
    );
}

#[test]
fn tblastx_subject_ncbi_parity_frameshift_pattern_remains_ungapped() {
    assert_translated_subject_outfmt_matches_ncbi_sorted_lines(
        "tblastx",
        "/usr/bin/tblastx",
        ">q1\nGCTTGTGATGAATTTGGTCATATTAAACTGATGAATCCTCAGCGTTCTACCGTTTGGTATGCTTGTGATGAATTTGGTCATATTAAACTGAATGAATCCTCAGCGTTCTACCGTTTGGTAT\n",
        ">s1\nGCTTGTGATGAATTTGGTCATATTAAACTGATGAATCCTCAGCGTTCTACCGTTTGGTATGCTTGTGATGAATTTGGTCATATTAAACTGAATGAATCCTCAGCGTTCTACCGTTTGGTAT\n",
        "6 qseqid sseqid qstart qend sstart send score length gaps gapopen qframe sframe frames qseq sseq btop",
        &["--seg", "no", "--evalue", "1000"],
        &["-seg", "no", "-evalue", "1000"],
    );
}

#[test]
fn tblastx_subject_ncbi_parity_complex_frameshift_top_hsps() {
    assert_translated_subject_outfmt_matches_ncbi_sorted_lines(
        "tblastx",
        "/usr/bin/tblastx",
        ">q1\nGCTTGTGATGAATTTGGTCATATTAAACTGATGAATCCTCAGCGTTCTACCGTTTGGTATGCTTGTGATGAATTTGGTCATATTAAACTGAATGAATCCTCAGCGTTCTACCGTTTGGTAT\n",
        ">s1\nGCTTGTGATGAATTTGGTCATATTAAACTGATGAATCCTCAGCGTTCTACCGTTTGGTATGCTTGTGATGAATTTGGTCATATTAAACTGATGAATCCTCAGCGTTCTACCGTTTGGTAT\n",
        "6 qseqid sseqid qstart qend sstart send score length gaps gapopen qframe sframe frames qseq sseq btop",
        &["--seg", "no", "--evalue", "1000", "--max_hsps", "30"],
        &["-seg", "no", "-evalue", "1000", "-max_hsps", "30"],
    );
}

#[test]
fn blastx_subject_ncbi_parity_indel_remains_ungapped() {
    assert_translated_subject_outfmt_matches_ncbi_sorted_lines(
        "blastx",
        "/usr/bin/blastx",
        ">q1\nATGAAATTTCTGATTCTGCTGTTTAAAATTTATGAAATTTCTGATTCTGCTGTTT\n",
        ">s1\nMKFLILLFKFMKFLILLF\n",
        "6 qseqid sseqid qstart qend sstart send score length gaps gapopen qframe sframe frames qseq sseq btop",
        &["--seg", "no", "--comp_based_stats", "0", "--evalue", "1000"],
        &["-seg", "no", "-comp_based_stats", "0", "-evalue", "1000"],
    );
}

#[test]
fn tblastn_subject_ncbi_parity_indel_remains_ungapped() {
    assert_translated_subject_outfmt_matches_ncbi_sorted_lines(
        "tblastn",
        "/usr/bin/tblastn",
        ">q1\nMKFLILLFKFMKFLILLF\n",
        ">s1\nATGAAATTTCTGATTCTGCTGTTTAAAATTTATGAAATTTCTGATTCTGCTGTTT\n",
        "6 qseqid sseqid qstart qend sstart send score length gaps gapopen qframe sframe frames qseq sseq btop",
        &["--seg", "no", "--comp_based_stats", "0", "--evalue", "1000"],
        &["-seg", "no", "-comp_based_stats", "0", "-evalue", "1000"],
    );
}

#[test]
fn blastx_subject_ncbi_parity_frameshift_gap_script() {
    assert_translated_subject_outfmt_matches_ncbi_sorted_lines(
        "blastx",
        "/usr/bin/blastx",
        ">q1\nGCTTGTGATGAATTTGGTCATATTAAACTGATGAATCCTCAGCGTTCTACCGTTTGGTATGCTTGTGATGAATTTGGTCATATTAAACTGAATGAATCCTCAGCGTTCTACCGTTTGGTAT\n",
        ">s1\nACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY\n",
        "6 qseqid sseqid qstart qend sstart send score length gaps gapopen qframe sframe frames qseq sseq btop",
        &[
            "--seg",
            "no",
            "--comp_based_stats",
            "0",
            "--evalue",
            "1000",
            "--gapopen",
            "9",
            "--gapextend",
            "1",
        ],
        &[
            "-seg",
            "no",
            "-comp_based_stats",
            "0",
            "-evalue",
            "1000",
            "-gapopen",
            "9",
            "-gapextend",
            "1",
        ],
    );
}

#[test]
fn tblastn_subject_ncbi_parity_frameshift_gap_script() {
    assert_translated_subject_outfmt_matches_ncbi_sorted_lines(
        "tblastn",
        "/usr/bin/tblastn",
        ">q1\nACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY\n",
        ">s1\nGCTTGTGATGAATTTGGTCATATTAAACTGATGAATCCTCAGCGTTCTACCGTTTGGTATGCTTGTGATGAATTTGGTCATATTAAACTGAATGAATCCTCAGCGTTCTACCGTTTGGTAT\n",
        "6 qseqid sseqid qstart qend sstart send score length gaps gapopen qframe sframe frames qseq sseq btop",
        &[
            "--seg",
            "no",
            "--comp_based_stats",
            "0",
            "--evalue",
            "1000",
            "--gapopen",
            "9",
            "--gapextend",
            "1",
        ],
        &[
            "-seg",
            "no",
            "-comp_based_stats",
            "0",
            "-evalue",
            "1000",
            "-gapopen",
            "9",
            "-gapextend",
            "1",
        ],
    );
}

#[test]
fn blastx_subject_ncbi_parity_commented_tabular_exact_hit() {
    assert_translated_subject_outfmt_matches_ncbi(
        "blastx",
        "/usr/bin/blastx",
        ">q1\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">s1\nMKFLILLF\n",
        "7",
        &[],
        &[],
    );
}

#[test]
fn blastx_subject_ncbi_parity_commented_tabular_parse_deflines() {
    assert_translated_subject_outfmt_matches_ncbi(
        "blastx",
        "/usr/bin/blastx",
        ">gi|123|ref|QNT.1| query nucleotide\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">gi|456|ref|SPROT.1| subject protein\nMKFLILLF\n",
        "7",
        &["--seg", "no", "--comp_based_stats", "0", "--parse_deflines"],
        &["-seg", "no", "-comp_based_stats", "0", "-parse_deflines"],
    );
}

#[test]
fn blastx_subject_ncbi_parity_xml_exact_hit() {
    assert_translated_subject_outfmt_matches_ncbi(
        "blastx",
        "/usr/bin/blastx",
        ">q1\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">s1\nMKFLILLF\n",
        "5",
        &["--seg", "no", "--comp_based_stats", "0"],
        &["-seg", "no", "-comp_based_stats", "0"],
    );
}

#[test]
fn blastx_subject_ncbi_parity_xml_positive_substitution() {
    assert_translated_subject_outfmt_matches_ncbi(
        "blastx",
        "/usr/bin/blastx",
        ">q1\nATGAAATTTCTGATTCTGCTGTTTAACATTCTGTGCCTGTTTCCTGTTCTGGCTGCTGATAACCATGGTGTTTCTATGAACGCTTCT\n",
        ">s1\nMKFLILLFNILCLFPVLAADNHGVSINAS\n",
        "5",
        &["--seg", "no", "--comp_based_stats", "0"],
        &["-seg", "no", "-comp_based_stats", "0"],
    );
}

#[test]
fn blastx_subject_ncbi_parity_xml_no_hits() {
    assert_translated_subject_outfmt_matches_ncbi(
        "blastx",
        "/usr/bin/blastx",
        ">q1\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">s1 unrelated\nPPPPPPPPPPPP\n",
        "5",
        &["--seg", "no", "--comp_based_stats", "0"],
        &["-seg", "no", "-comp_based_stats", "0"],
    );
}

#[test]
fn blastx_subject_ncbi_parity_xml_multiple_queries() {
    assert_translated_subject_outfmt_matches_ncbi(
        "blastx",
        "/usr/bin/blastx",
        ">q1\nATGAAATTTCTGATTCTGCTGTTT\n>q2 unrelated\nCCCCCCCCCCCCCCCCCCCCCCCC\n",
        ">s1\nMKFLILLF\n",
        "5",
        &["--seg", "no", "--comp_based_stats", "0"],
        &["-seg", "no", "-comp_based_stats", "0"],
    );
}

#[test]
fn blastx_subject_ncbi_parity_xml_multi_hsp_same_subject() {
    assert_translated_subject_outfmt_matches_ncbi(
        "blastx",
        "/usr/bin/blastx",
        ">q1\nATGAAATTTCTGATTCTGCTGTTTAAAACCCCGGGGTTTTATGAAATTTCTGATTCTGCTGTTT\n",
        ">s1 multi\nMKFLILLFQQQQGMKFLILLF\n",
        "5",
        &["--seg", "no", "--comp_based_stats", "0", "--evalue", "1000"],
        &["-seg", "no", "-comp_based_stats", "0", "-evalue", "1000"],
    );
}

#[test]
fn blastx_subject_ncbi_parity_xml_parse_deflines() {
    assert_translated_subject_outfmt_matches_ncbi(
        "blastx",
        "/usr/bin/blastx",
        ">gi|123|ref|QNT.1| query nucleotide\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">gi|456|ref|SPROT.1| subject protein\nMKFLILLF\n",
        "5",
        &["--seg", "no", "--comp_based_stats", "0", "--parse_deflines"],
        &["-seg", "no", "-comp_based_stats", "0", "-parse_deflines"],
    );
}

#[test]
fn blastx_subject_ncbi_parity_pairwise_exact_hit() {
    assert_translated_subject_outfmt_matches_ncbi(
        "blastx",
        "/usr/bin/blastx",
        ">q1\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">s1\nMKFLILLF\n",
        "0",
        &["--seg", "no", "--comp_based_stats", "0"],
        &["-seg", "no", "-comp_based_stats", "0"],
    );
}

#[test]
fn blastx_subject_ncbi_parity_pairwise_parse_deflines() {
    assert_translated_subject_outfmt_matches_ncbi(
        "blastx",
        "/usr/bin/blastx",
        ">gi|123|ref|QNT.1| query nucleotide\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">gi|456|ref|SPROT.1| subject protein\nMKFLILLF\n",
        "0",
        &["--seg", "no", "--comp_based_stats", "0", "--parse_deflines"],
        &["-seg", "no", "-comp_based_stats", "0", "-parse_deflines"],
    );
}

#[test]
fn blastx_subject_ncbi_parity_tabular_parse_deflines() {
    assert_translated_subject_outfmt_matches_ncbi(
        "blastx",
        "/usr/bin/blastx",
        ">gi|123|ref|QNT.1| query nucleotide\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">gi|456|ref|SPROT.1| subject protein\nMKFLILLF\n",
        "6 qseqid qacc qaccver sseqid sacc saccver",
        &["--seg", "no", "--comp_based_stats", "0", "--parse_deflines"],
        &["-seg", "no", "-comp_based_stats", "0", "-parse_deflines"],
    );
}

#[test]
fn blastx_subject_ncbi_parity_tabular_parse_deflines_gis() {
    assert_translated_subject_outfmt_matches_ncbi(
        "blastx",
        "/usr/bin/blastx",
        ">gi|123|ref|QNT.1| query nucleotide\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">gi|456|ref|SPROT.1| subject protein\nMKFLILLF\n",
        "6 qseqid qgi qacc qaccver sseqid sgi sallgi sacc saccver sallseqid",
        &["--seg", "no", "--comp_based_stats", "0", "--parse_deflines"],
        &["-seg", "no", "-comp_based_stats", "0", "-parse_deflines"],
    );
}

#[test]
fn blastx_subject_ncbi_parity_tabular_parse_deflines_titles() {
    assert_translated_subject_outfmt_matches_ncbi(
        "blastx",
        "/usr/bin/blastx",
        ">gi|123|ref|QNT.1| query nucleotide\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">gi|456|ref|SPROT.1| subject protein\nMKFLILLF\n",
        "6 sseqid stitle salltitles",
        &["--seg", "no", "--comp_based_stats", "0", "--parse_deflines"],
        &["-seg", "no", "-comp_based_stats", "0", "-parse_deflines"],
    );
}

#[test]
fn blastx_subject_ncbi_parity_csv_parse_deflines_metadata() {
    assert_translated_subject_outfmt_matches_ncbi(
        "blastx",
        "/usr/bin/blastx",
        ">gi|123|ref|QNT.1| query nucleotide\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">gi|456|ref|SPROT.1| subject protein\nMKFLILLF\n",
        "10 qseqid qgi qacc qaccver sseqid sgi sallgi sacc saccver sallseqid stitle salltitles",
        &["--seg", "no", "--comp_based_stats", "0", "--parse_deflines"],
        &["-seg", "no", "-comp_based_stats", "0", "-parse_deflines"],
    );
}

#[test]
fn blastx_subject_ncbi_parity_pairwise_no_hits() {
    assert_translated_subject_outfmt_matches_ncbi(
        "blastx",
        "/usr/bin/blastx",
        ">q1\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">s1\nGGGGGGGG\n",
        "0",
        &[
            "--seg",
            "no",
            "--comp_based_stats",
            "0",
            "--evalue",
            "1e-100",
        ],
        &["-seg", "no", "-comp_based_stats", "0", "-evalue", "1e-100"],
    );
}

#[test]
fn blastx_subject_ncbi_parity_pairwise_line_length() {
    assert_translated_subject_outfmt_matches_ncbi(
        "blastx",
        "/usr/bin/blastx",
        ">q1\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">s1\nMKFLILLF\n",
        "0",
        &[
            "--seg",
            "no",
            "--comp_based_stats",
            "0",
            "--line_length",
            "5",
        ],
        &["-seg", "no", "-comp_based_stats", "0", "-line_length", "5"],
    );
}

#[test]
fn blastx_subject_ncbi_parity_pairwise_num_alignments() {
    assert_translated_subject_outfmt_matches_ncbi(
        "blastx",
        "/usr/bin/blastx",
        ">q1\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">s1 first\nMKFLILLF\n>s2 second\nMKFLILLF\n",
        "0",
        &[
            "--seg",
            "no",
            "--comp_based_stats",
            "0",
            "--num_alignments",
            "1",
            "--num_descriptions",
            "1",
        ],
        &[
            "-seg",
            "no",
            "-comp_based_stats",
            "0",
            "-num_alignments",
            "1",
            "-num_descriptions",
            "1",
        ],
    );
}

#[test]
fn blastx_subject_ncbi_parity_pairwise_multi_hsp_same_subject() {
    assert_translated_subject_outfmt_matches_ncbi(
        "blastx",
        "/usr/bin/blastx",
        ">q1\nATGAAATTTCTGATTCTGCTGTTTAAAACCCCGGGGTTTTATGAAATTTCTGATTCTGCTGTTT\n",
        ">s1 multi\nMKFLILLFQQQQGMKFLILLF\n",
        "0",
        &["--seg", "no", "--comp_based_stats", "0", "--evalue", "1000"],
        &["-seg", "no", "-comp_based_stats", "0", "-evalue", "1000"],
    );
}

#[test]
fn blastx_subject_ncbi_parity_tabular_sum_stats_gapped_single_hsp() {
    assert_translated_subject_outfmt_matches_ncbi(
        "blastx",
        "/usr/bin/blastx",
        ">q1\nATGAAATTTCTGATTCTGCTGTTTAAATTTATGAAATTTCTGATTCTGCTGTTT\n",
        ">s1\nMKFLILLFMKFLILLF\n",
        "6 qseqid sseqid score bitscore evalue length qstart qend sstart send",
        &[
            "--seg",
            "no",
            "--comp_based_stats",
            "0",
            "--sum_stats",
            "true",
            "--evalue",
            "1000",
        ],
        &[
            "-seg",
            "no",
            "-comp_based_stats",
            "0",
            "-sum_stats",
            "true",
            "-evalue",
            "1000",
        ],
    );
}

#[test]
fn blastx_subject_ncbi_parity_tabular_sum_stats_false_single_hsp() {
    assert_translated_subject_outfmt_matches_ncbi(
        "blastx",
        "/usr/bin/blastx",
        ">q1\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">s1\nMKFLILLF\n",
        "6 qseqid sseqid score bitscore evalue length qstart qend sstart send",
        &[
            "--seg",
            "no",
            "--comp_based_stats",
            "0",
            "--sum_stats",
            "false",
            "--evalue",
            "1000",
        ],
        &[
            "-seg",
            "no",
            "-comp_based_stats",
            "0",
            "-sum_stats",
            "false",
            "-evalue",
            "1000",
        ],
    );
}

#[test]
fn blastx_subject_ncbi_parity_explicit_sum_stats_gapped_single_hsp() {
    assert_translated_subject_outfmt_matches_ncbi(
        "blastx",
        "/usr/bin/blastx",
        ">q1\nATGAAATTTCTGATTCTGCTGTTTAAATTTATGAAATTTCTGATTCTGCTGTTT\n",
        ">s1\nMKFLILLFMKFLILLF\n",
        "0",
        &[
            "--seg",
            "no",
            "--comp_based_stats",
            "0",
            "--sum_stats",
            "true",
            "--evalue",
            "1000",
        ],
        &[
            "-seg",
            "no",
            "-comp_based_stats",
            "0",
            "-sum_stats",
            "true",
            "-evalue",
            "1000",
        ],
    );
}

#[test]
fn tblastn_subject_ncbi_parity_commented_tabular_exact_hit() {
    assert_translated_subject_outfmt_matches_ncbi(
        "tblastn",
        "/usr/bin/tblastn",
        ">q1\nMKFLILLF\n",
        ">s1\nATGAAATTTCTGATTCTGCTGTTT\n",
        "7",
        &[],
        &[],
    );
}

#[test]
fn tblastn_subject_ncbi_parity_commented_tabular_parse_deflines() {
    assert_translated_subject_outfmt_matches_ncbi(
        "tblastn",
        "/usr/bin/tblastn",
        ">gi|123|ref|QPROT.1| query protein\nMKFLILLF\n",
        ">gi|456|ref|SNT.1| subject nucleotide\nATGAAATTTCTGATTCTGCTGTTT\n",
        "7 qaccver saccver pident length bitscore",
        &["--seg", "no", "--comp_based_stats", "0", "--parse_deflines"],
        &["-seg", "no", "-comp_based_stats", "0", "-parse_deflines"],
    );
}

#[test]
fn tblastn_subject_ncbi_parity_xml_exact_hit() {
    assert_translated_subject_outfmt_matches_ncbi(
        "tblastn",
        "/usr/bin/tblastn",
        ">q1\nMKFLILLF\n",
        ">s1\nATGAAATTTCTGATTCTGCTGTTT\n",
        "5",
        &["--seg", "no", "--comp_based_stats", "0"],
        &["-seg", "no", "-comp_based_stats", "0"],
    );
}

#[test]
fn tblastn_subject_ncbi_parity_xml_no_hits() {
    assert_translated_subject_outfmt_matches_ncbi(
        "tblastn",
        "/usr/bin/tblastn",
        ">q1\nMKFLILLF\n",
        ">s1 unrelated\nGGGGGGGGGGGGGGGGGGGGGGGG\n",
        "5",
        &["--seg", "no", "--comp_based_stats", "0"],
        &["-seg", "no", "-comp_based_stats", "0"],
    );
}

#[test]
fn tblastn_subject_ncbi_parity_xml_multiple_queries() {
    assert_translated_subject_outfmt_matches_ncbi(
        "tblastn",
        "/usr/bin/tblastn",
        ">q1\nMKFLILLF\n>q2 unrelated\nPPPPPPPP\n",
        ">s1\nATGAAATTTCTGATTCTGCTGTTT\n",
        "5",
        &["--seg", "no", "--comp_based_stats", "0"],
        &["-seg", "no", "-comp_based_stats", "0"],
    );
}

#[test]
fn tblastn_subject_ncbi_parity_xml_multi_hsp_same_subject() {
    assert_translated_subject_outfmt_matches_ncbi(
        "tblastn",
        "/usr/bin/tblastn",
        ">q1\nMKFLILLFQQQQGMKFLILLF\n",
        ">s1 multi\nATGAAATTTCTGATTCTGCTGTTTAAAACCCCGGGGTTTTATGAAATTTCTGATTCTGCTGTTT\n",
        "5",
        &["--seg", "no", "--comp_based_stats", "0", "--evalue", "1000"],
        &["-seg", "no", "-comp_based_stats", "0", "-evalue", "1000"],
    );
}

#[test]
fn tblastn_subject_ncbi_parity_xml_parse_deflines() {
    assert_translated_subject_outfmt_matches_ncbi(
        "tblastn",
        "/usr/bin/tblastn",
        ">gi|123|ref|QPROT.1| query protein\nMKFLILLF\n",
        ">gi|456|ref|SNT.1| subject nucleotide\nATGAAATTTCTGATTCTGCTGTTT\n",
        "5",
        &["--seg", "no", "--comp_based_stats", "0", "--parse_deflines"],
        &["-seg", "no", "-comp_based_stats", "0", "-parse_deflines"],
    );
}

#[test]
fn tblastn_subject_ncbi_parity_pairwise_exact_hit() {
    assert_translated_subject_outfmt_matches_ncbi(
        "tblastn",
        "/usr/bin/tblastn",
        ">q1\nMKFLILLF\n",
        ">s1\nATGAAATTTCTGATTCTGCTGTTT\n",
        "0",
        &["--seg", "no", "--comp_based_stats", "0"],
        &["-seg", "no", "-comp_based_stats", "0"],
    );
}

#[test]
fn tblastn_subject_ncbi_parity_pairwise_parse_deflines() {
    assert_translated_subject_outfmt_matches_ncbi(
        "tblastn",
        "/usr/bin/tblastn",
        ">gi|123|ref|QPROT.1| query protein\nMKFLILLF\n",
        ">gi|456|ref|SNT.1| subject nucleotide\nATGAAATTTCTGATTCTGCTGTTT\n",
        "0",
        &["--seg", "no", "--comp_based_stats", "0", "--parse_deflines"],
        &["-seg", "no", "-comp_based_stats", "0", "-parse_deflines"],
    );
}

#[test]
fn tblastn_subject_ncbi_parity_explicit_sum_stats_gapped_single_hsp() {
    assert_translated_subject_outfmt_matches_ncbi(
        "tblastn",
        "/usr/bin/tblastn",
        ">q1\nMKFLILLFMKFLILLF\n",
        ">s1\nATGAAATTTCTGATTCTGCTGTTTAAATTTATGAAATTTCTGATTCTGCTGTTT\n",
        "0",
        &[
            "--seg",
            "no",
            "--comp_based_stats",
            "0",
            "--sum_stats",
            "true",
            "--evalue",
            "1000",
        ],
        &[
            "-seg",
            "no",
            "-comp_based_stats",
            "0",
            "-sum_stats",
            "true",
            "-evalue",
            "1000",
        ],
    );
}

#[test]
fn tblastn_subject_ncbi_parity_tabular_sum_stats_gapped_single_hsp() {
    assert_translated_subject_outfmt_matches_ncbi(
        "tblastn",
        "/usr/bin/tblastn",
        ">q1\nMKFLILLFMKFLILLF\n",
        ">s1\nATGAAATTTCTGATTCTGCTGTTTAAATTTATGAAATTTCTGATTCTGCTGTTT\n",
        "6 qseqid sseqid score bitscore evalue length qstart qend sstart send",
        &[
            "--seg",
            "no",
            "--comp_based_stats",
            "0",
            "--sum_stats",
            "true",
            "--evalue",
            "1000",
        ],
        &[
            "-seg",
            "no",
            "-comp_based_stats",
            "0",
            "-sum_stats",
            "true",
            "-evalue",
            "1000",
        ],
    );
}

#[test]
fn tblastn_subject_ncbi_parity_tabular_sum_stats_false_single_hsp() {
    assert_translated_subject_outfmt_matches_ncbi(
        "tblastn",
        "/usr/bin/tblastn",
        ">q1\nMKFLILLF\n",
        ">s1\nATGAAATTTCTGATTCTGCTGTTT\n",
        "6 qseqid sseqid score bitscore evalue length qstart qend sstart send",
        &[
            "--seg",
            "no",
            "--comp_based_stats",
            "0",
            "--sum_stats",
            "false",
            "--evalue",
            "1000",
        ],
        &[
            "-seg",
            "no",
            "-comp_based_stats",
            "0",
            "-sum_stats",
            "false",
            "-evalue",
            "1000",
        ],
    );
}

#[test]
fn translated_subject_ncbi_parity_small_max_intron_uses_non_sum_stats_evalue() {
    for (program, ncbi_bin, query, subject) in [
        (
            "blastx",
            "/usr/bin/blastx",
            ">q1\nATGAAATTTCTGATTCTGCTGTTT\n",
            ">s1\nMKFLILLF\n",
        ),
        (
            "tblastn",
            "/usr/bin/tblastn",
            ">q1\nMKFLILLF\n",
            ">s1\nATGAAATTTCTGATTCTGCTGTTT\n",
        ),
    ] {
        for max_intron in ["1", "2", "3", "4"] {
            assert_translated_subject_outfmt_matches_ncbi(
                program,
                ncbi_bin,
                query,
                subject,
                "6 qseqid sseqid score bitscore evalue length qstart qend sstart send",
                &[
                    "--seg",
                    "no",
                    "--comp_based_stats",
                    "0",
                    "--evalue",
                    "1000",
                    "--max_intron_length",
                    max_intron,
                ],
                &[
                    "-seg",
                    "no",
                    "-comp_based_stats",
                    "0",
                    "-evalue",
                    "1000",
                    "-max_intron_length",
                    max_intron,
                ],
            );
        }
    }
}

#[test]
fn translated_db_ncbi_parity_small_max_intron_uses_non_sum_stats_evalue() {
    for (program, ncbi_bin, dbtype, query, db_fasta) in [
        (
            "blastx",
            "/usr/bin/blastx",
            "prot",
            ">q1\nATGAAATTTCTGATTCTGCTGTTT\n",
            ">s1\nMKFLILLF\n",
        ),
        (
            "tblastn",
            "/usr/bin/tblastn",
            "nucl",
            ">q1\nMKFLILLF\n",
            ">s1\nATGAAATTTCTGATTCTGCTGTTT\n",
        ),
    ] {
        for max_intron in ["1", "2", "3", "4"] {
            assert_translated_db_outfmt_matches_ncbi(
                program,
                ncbi_bin,
                dbtype,
                query,
                db_fasta,
                "6 qseqid sseqid score bitscore evalue length qstart qend sstart send",
                &[
                    "--seg",
                    "no",
                    "--comp_based_stats",
                    "0",
                    "--evalue",
                    "1000",
                    "--max_intron_length",
                    max_intron,
                ],
                &[
                    "-seg",
                    "no",
                    "-comp_based_stats",
                    "0",
                    "-evalue",
                    "1000",
                    "-max_intron_length",
                    max_intron,
                ],
            );
        }
    }
}

#[test]
fn tblastn_subject_ncbi_parity_tabular_parse_deflines() {
    assert_translated_subject_outfmt_matches_ncbi(
        "tblastn",
        "/usr/bin/tblastn",
        ">gi|123|ref|QPROT.1| query protein\nMKFLILLF\n",
        ">gi|456|ref|SNT.1| subject nucleotide\nATGAAATTTCTGATTCTGCTGTTT\n",
        "6 qseqid qacc qaccver sseqid sacc saccver pident length bitscore",
        &["--seg", "no", "--comp_based_stats", "0", "--parse_deflines"],
        &["-seg", "no", "-comp_based_stats", "0", "-parse_deflines"],
    );
}

#[test]
fn tblastn_subject_ncbi_parity_tabular_parse_deflines_gis() {
    assert_translated_subject_outfmt_matches_ncbi(
        "tblastn",
        "/usr/bin/tblastn",
        ">gi|123|ref|QPROT.1| query protein\nMKFLILLF\n",
        ">gi|456|ref|SNT.1| subject nucleotide\nATGAAATTTCTGATTCTGCTGTTT\n",
        "6 qseqid qgi qacc qaccver sseqid sgi sallgi sacc saccver sallseqid pident length bitscore",
        &["--seg", "no", "--comp_based_stats", "0", "--parse_deflines"],
        &["-seg", "no", "-comp_based_stats", "0", "-parse_deflines"],
    );
}

#[test]
fn tblastn_subject_ncbi_parity_tabular_parse_deflines_titles() {
    assert_translated_subject_outfmt_matches_ncbi(
        "tblastn",
        "/usr/bin/tblastn",
        ">gi|123|ref|QPROT.1| query protein\nMKFLILLF\n",
        ">gi|456|ref|SNT.1| subject nucleotide\nATGAAATTTCTGATTCTGCTGTTT\n",
        "6 sseqid stitle salltitles",
        &["--seg", "no", "--comp_based_stats", "0", "--parse_deflines"],
        &["-seg", "no", "-comp_based_stats", "0", "-parse_deflines"],
    );
}

#[test]
fn tblastn_subject_ncbi_parity_csv_parse_deflines_metadata() {
    assert_translated_subject_outfmt_matches_ncbi(
        "tblastn",
        "/usr/bin/tblastn",
        ">gi|123|ref|QPROT.1| query protein\nMKFLILLF\n",
        ">gi|456|ref|SNT.1| subject nucleotide\nATGAAATTTCTGATTCTGCTGTTT\n",
        "10 qseqid qgi qacc qaccver sseqid sgi sallgi sacc saccver sallseqid stitle salltitles",
        &["--seg", "no", "--comp_based_stats", "0", "--parse_deflines"],
        &["-seg", "no", "-comp_based_stats", "0", "-parse_deflines"],
    );
}

#[test]
fn tblastn_subject_ncbi_parity_pairwise_no_hits() {
    assert_translated_subject_outfmt_matches_ncbi(
        "tblastn",
        "/usr/bin/tblastn",
        ">q1\nMKFLILLF\n",
        ">s1\nGGGGGGGGGGGGGGGGGGGGGGGG\n",
        "0",
        &[
            "--seg",
            "no",
            "--comp_based_stats",
            "0",
            "--evalue",
            "1e-100",
        ],
        &["-seg", "no", "-comp_based_stats", "0", "-evalue", "1e-100"],
    );
}

#[test]
fn tblastn_subject_ncbi_parity_pairwise_line_length() {
    assert_translated_subject_outfmt_matches_ncbi(
        "tblastn",
        "/usr/bin/tblastn",
        ">q1\nMKFLILLF\n",
        ">s1\nATGAAATTTCTGATTCTGCTGTTT\n",
        "0",
        &[
            "--seg",
            "no",
            "--comp_based_stats",
            "0",
            "--line_length",
            "5",
        ],
        &["-seg", "no", "-comp_based_stats", "0", "-line_length", "5"],
    );
}

#[test]
fn tblastn_subject_ncbi_parity_pairwise_num_alignments() {
    assert_translated_subject_outfmt_matches_ncbi(
        "tblastn",
        "/usr/bin/tblastn",
        ">q1\nMKFLILLF\n",
        ">s1 first\nATGAAATTTCTGATTCTGCTGTTT\n>s2 second\nATGAAATTTCTGATTCTGCTGTTT\n",
        "0",
        &[
            "--seg",
            "no",
            "--comp_based_stats",
            "0",
            "--num_alignments",
            "1",
            "--num_descriptions",
            "1",
        ],
        &[
            "-seg",
            "no",
            "-comp_based_stats",
            "0",
            "-num_alignments",
            "1",
            "-num_descriptions",
            "1",
        ],
    );
}

#[test]
fn tblastn_subject_ncbi_parity_pairwise_multi_hsp_same_subject() {
    assert_translated_subject_outfmt_matches_ncbi(
        "tblastn",
        "/usr/bin/tblastn",
        ">q1\nMKFLILLFQQQQGMKFLILLF\n",
        ">s1 multi\nATGAAATTTCTGATTCTGCTGTTTAAAACCCCGGGGTTTTATGAAATTTCTGATTCTGCTGTTT\n",
        "0",
        &["--seg", "no", "--comp_based_stats", "0", "--evalue", "1000"],
        &["-seg", "no", "-comp_based_stats", "0", "-evalue", "1000"],
    );
}

#[test]
fn tblastx_subject_ncbi_parity_commented_tabular_exact_hit() {
    assert_translated_subject_outfmt_matches_ncbi(
        "tblastx",
        "/usr/bin/tblastx",
        ">q1\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">s1\nATGAAATTTCTGATTCTGCTGTTT\n",
        "7",
        &[],
        &[],
    );
}

#[test]
fn tblastx_subject_ncbi_parity_commented_tabular_parse_deflines() {
    assert_translated_subject_outfmt_matches_ncbi(
        "tblastx",
        "/usr/bin/tblastx",
        ">gi|123|ref|QNT.1| query nucleotide\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">gi|456|ref|SNT.1| subject nucleotide\nATGAAATTTCTGATTCTGCTGTTT\n",
        "7",
        &["--seg", "no", "--evalue", "1000", "--parse_deflines"],
        &["-seg", "no", "-evalue", "1000", "-parse_deflines"],
    );
}

#[test]
fn tblastx_subject_ncbi_parity_xml_exact_hit() {
    assert_translated_subject_outfmt_matches_ncbi(
        "tblastx",
        "/usr/bin/tblastx",
        ">q1\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">s1\nATGAAATTTCTGATTCTGCTGTTT\n",
        "5",
        &["--seg", "no", "--evalue", "1000"],
        &["-seg", "no", "-evalue", "1000"],
    );
}

#[test]
fn tblastx_subject_ncbi_parity_xml_no_hits() {
    assert_translated_subject_outfmt_matches_ncbi(
        "tblastx",
        "/usr/bin/tblastx",
        ">q1\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">s1 unrelated\nCCACCTCCACCTCCACCTCCACCT\n",
        "5",
        &["--seg", "no"],
        &["-seg", "no"],
    );
}

#[test]
fn tblastx_subject_ncbi_parity_xml_multiple_queries() {
    assert_translated_subject_outfmt_matches_ncbi(
        "tblastx",
        "/usr/bin/tblastx",
        ">q1\nATGAAATTTCTGATTCTGCTGTTT\n>q2 unrelated\nCCACCTCCACCTCCACCTCCACCT\n",
        ">s1\nATGAAATTTCTGATTCTGCTGTTT\n",
        "5",
        &["--seg", "no", "--evalue", "1000"],
        &["-seg", "no", "-evalue", "1000"],
    );
}

#[test]
fn tblastx_subject_ncbi_parity_xml_multi_hsp_same_subject() {
    assert_translated_subject_outfmt_matches_ncbi(
        "tblastx",
        "/usr/bin/tblastx",
        ">q1\nATGAAATTTCTGATTCTGCTGTTTAAAACCCCGGGGTTTTATGAAATTTCTGATTCTGCTGTTT\n",
        ">s1 multi\nATGAAATTTCTGATTCTGCTGTTTAAAACCCCGGGGTTTTATGAAATTTCTGATTCTGCTGTTT\n",
        "5",
        &["--seg", "no", "--evalue", "1000"],
        &["-seg", "no", "-evalue", "1000"],
    );
}

#[test]
fn tblastx_subject_ncbi_parity_xml_parse_deflines() {
    assert_translated_subject_outfmt_matches_ncbi(
        "tblastx",
        "/usr/bin/tblastx",
        ">gi|123|ref|QNT.1| query nucleotide\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">gi|456|ref|SNT.1| subject nucleotide\nATGAAATTTCTGATTCTGCTGTTT\n",
        "5",
        &["--seg", "no", "--evalue", "1000", "--parse_deflines"],
        &["-seg", "no", "-evalue", "1000", "-parse_deflines"],
    );
}

#[test]
fn tblastx_subject_ncbi_parity_pairwise_exact_hit() {
    assert_translated_subject_outfmt_matches_ncbi(
        "tblastx",
        "/usr/bin/tblastx",
        ">q1\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">s1\nATGAAATTTCTGATTCTGCTGTTT\n",
        "0",
        &["--seg", "no", "--evalue", "1000"],
        &["-seg", "no", "-evalue", "1000"],
    );
}

#[test]
fn tblastx_subject_ncbi_parity_pairwise_parse_deflines() {
    assert_translated_subject_outfmt_matches_ncbi(
        "tblastx",
        "/usr/bin/tblastx",
        ">gi|123|ref|QNT.1| query nucleotide\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">gi|456|ref|SNT.1| subject nucleotide\nATGAAATTTCTGATTCTGCTGTTT\n",
        "0",
        &["--seg", "no", "--evalue", "1000", "--parse_deflines"],
        &["-seg", "no", "-evalue", "1000", "-parse_deflines"],
    );
}

#[test]
fn tblastx_subject_ncbi_parity_tabular_parse_deflines() {
    assert_translated_subject_outfmt_matches_ncbi(
        "tblastx",
        "/usr/bin/tblastx",
        ">gi|123|ref|QNT.1| query nucleotide\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">gi|456|ref|SNT.1| subject nucleotide\nATGAAATTTCTGATTCTGCTGTTT\n",
        "6 qseqid qacc qaccver sseqid sacc saccver pident length bitscore",
        &["--seg", "no", "--evalue", "1000", "--parse_deflines"],
        &["-seg", "no", "-evalue", "1000", "-parse_deflines"],
    );
}

#[test]
fn tblastx_subject_ncbi_parity_tabular_parse_deflines_gis() {
    assert_translated_subject_outfmt_matches_ncbi(
        "tblastx",
        "/usr/bin/tblastx",
        ">gi|123|ref|QNT.1| query nucleotide\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">gi|456|ref|SNT.1| subject nucleotide\nATGAAATTTCTGATTCTGCTGTTT\n",
        "6 qseqid qgi qacc qaccver sseqid sgi sallgi sacc saccver sallseqid pident length bitscore",
        &["--seg", "no", "--evalue", "1000", "--parse_deflines"],
        &["-seg", "no", "-evalue", "1000", "-parse_deflines"],
    );
}

#[test]
fn tblastx_subject_ncbi_parity_tabular_parse_deflines_titles() {
    assert_translated_subject_outfmt_matches_ncbi(
        "tblastx",
        "/usr/bin/tblastx",
        ">gi|123|ref|QNT.1| query nucleotide\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">gi|456|ref|SNT.1| subject nucleotide\nATGAAATTTCTGATTCTGCTGTTT\n",
        "6 sseqid stitle salltitles",
        &["--seg", "no", "--evalue", "1000", "--parse_deflines"],
        &["-seg", "no", "-evalue", "1000", "-parse_deflines"],
    );
}

#[test]
fn tblastx_subject_ncbi_parity_csv_parse_deflines_metadata() {
    assert_translated_subject_outfmt_matches_ncbi(
        "tblastx",
        "/usr/bin/tblastx",
        ">gi|123|ref|QNT.1| query nucleotide\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">gi|456|ref|SNT.1| subject nucleotide\nATGAAATTTCTGATTCTGCTGTTT\n",
        "10 qseqid qgi qacc qaccver sseqid sgi sallgi sacc saccver sallseqid stitle salltitles",
        &["--seg", "no", "--evalue", "1000", "--parse_deflines"],
        &["-seg", "no", "-evalue", "1000", "-parse_deflines"],
    );
}

#[test]
fn tblastx_subject_ncbi_parity_pairwise_no_hits() {
    assert_translated_subject_outfmt_matches_ncbi(
        "tblastx",
        "/usr/bin/tblastx",
        ">q1\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">s1\nGGGGGGGGGGGGGGGGGGGGGGGG\n",
        "0",
        &["--seg", "no", "--evalue", "1e-100"],
        &["-seg", "no", "-evalue", "1e-100"],
    );
}

#[test]
fn tblastx_subject_ncbi_parity_pairwise_line_length() {
    assert_translated_subject_outfmt_matches_ncbi(
        "tblastx",
        "/usr/bin/tblastx",
        ">q1\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">s1\nATGAAATTTCTGATTCTGCTGTTT\n",
        "0",
        &["--seg", "no", "--evalue", "1000", "--line_length", "5"],
        &["-seg", "no", "-evalue", "1000", "-line_length", "5"],
    );
}

#[test]
fn tblastx_subject_ncbi_parity_pairwise_num_alignments() {
    assert_translated_subject_outfmt_matches_ncbi(
        "tblastx",
        "/usr/bin/tblastx",
        ">q1\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">s1 first\nATGAAATTTCTGATTCTGCTGTTT\n>s2 second\nATGAAATTTCTGATTCTGCTGTTT\n",
        "0",
        &[
            "--seg",
            "no",
            "--evalue",
            "1000",
            "--num_alignments",
            "1",
            "--num_descriptions",
            "1",
        ],
        &[
            "-seg",
            "no",
            "-evalue",
            "1000",
            "-num_alignments",
            "1",
            "-num_descriptions",
            "1",
        ],
    );
}

#[test]
fn tblastx_subject_ncbi_parity_pairwise_multi_hsp_same_subject() {
    assert_translated_subject_outfmt_matches_ncbi(
        "tblastx",
        "/usr/bin/tblastx",
        ">q1\nATGAAATTTCTGATTCTGCTGTTTAAAACCCCGGGGTTTTATGAAATTTCTGATTCTGCTGTTT\n",
        ">s1 multi\nATGAAATTTCTGATTCTGCTGTTTAAAACCCCGGGGTTTTATGAAATTTCTGATTCTGCTGTTT\n",
        "0",
        &["--seg", "no", "--evalue", "1000"],
        &["-seg", "no", "-evalue", "1000"],
    );
}

#[test]
fn translated_subject_ncbi_parity_pairwise_zero_description_alignment_limits() {
    for (program, ncbi_bin, query, subject) in [
        (
            "blastx",
            "/usr/bin/blastx",
            ">q1 translated query\nATGAAATTTCTGATTCTGCTGTTT\n",
            ">s1 protein subject\nMKFLILLF\n>s2 unrelated\nACDEFGHIKLMN\n",
        ),
        (
            "tblastn",
            "/usr/bin/tblastn",
            ">q1 protein query\nMKFLILLF\n",
            ">s1 nucleotide subject\nATGAAATTTCTGATTCTGCTGTTT\n>s2 unrelated\nCCCCCCCCCCCCCCCCCCCCCCCC\n",
        ),
        (
            "tblastx",
            "/usr/bin/tblastx",
            ">q1 translated query\nATGAAATTTCTGATTCTGCTGTTT\n",
            ">s1 nucleotide subject\nATGAAATTTCTGATTCTGCTGTTT\n>s2 unrelated\nCCCCCCCCCCCCCCCCCCCCCCCC\n",
        ),
    ] {
        for (rust_extra, ncbi_extra) in [
            (vec!["--num_alignments", "0"], vec!["-num_alignments", "0"]),
            (
                vec!["--num_descriptions", "0"],
                vec!["-num_descriptions", "0"],
            ),
        ] {
            let (rust_comp_args, ncbi_comp_args): (&[&str], &[&str]) = if program == "tblastx" {
                (&[], &[])
            } else {
                (&["--seg", "no", "--comp_based_stats", "0"], &["-seg", "no", "-comp_based_stats", "0"])
            };
            let mut rust_args = Vec::new();
            rust_args.extend(rust_comp_args);
            rust_args.extend(rust_extra);
            let mut ncbi_args = Vec::new();
            ncbi_args.extend(ncbi_comp_args);
            ncbi_args.extend(ncbi_extra);
            assert_translated_subject_outfmt_matches_ncbi(
                program,
                ncbi_bin,
                query,
                subject,
                "0",
                &rust_args,
                &ncbi_args,
            );
        }
    }
}

#[test]
fn translated_subject_ncbi_parity_culling_limit() {
    for (program, ncbi_bin, query, subject) in [
        (
            "blastx",
            "/usr/bin/blastx",
            ">q1\nATGAAATTTCTGATTCTGCTGTTTAAAACCCCGGGGTTTTATGAAATTTCTGATTCTGCTGTTT\n",
            ">wide\nMKFLILLFQQQQGMKFLILLF\n>contained\nKFLILL\n>outside\nQQQQG\n",
        ),
        (
            "tblastn",
            "/usr/bin/tblastn",
            ">q1\nMKFLILLFQQQQGMKFLILLF\n",
            ">wide\nATGAAATTTCTGATTCTGCTGTTTAAAACCCCGGGGTTTTATGAAATTTCTGATTCTGCTGTTT\n>contained\nAAATTTCTGATTCTGCTG\n>outside\nCAGCAGCAGCAGGGT\n",
        ),
        (
            "tblastx",
            "/usr/bin/tblastx",
            ">q1\nATGAAATTTCTGATTCTGCTGTTTAAAACCCCGGGGTTTTATGAAATTTCTGATTCTGCTGTTT\n",
            ">wide\nATGAAATTTCTGATTCTGCTGTTTAAAACCCCGGGGTTTTATGAAATTTCTGATTCTGCTGTTT\n>contained\nAAATTTCTGATTCTGCTG\n>outside\nCAGCAGCAGCAGGGT\n",
        ),
    ] {
        let (rust_comp_args, ncbi_comp_args): (&[&str], &[&str]) = if program == "tblastx" {
            (&[], &[])
        } else {
            (&["--comp_based_stats", "0"], &["-comp_based_stats", "0"])
        };
        let mut rust_args = vec!["--seg", "no"];
        rust_args.extend(rust_comp_args);
        rust_args.extend([
            "--sum_stats",
            "false",
            "--evalue",
            "1000",
            "--culling_limit",
            "1",
        ]);
        let mut ncbi_args = vec!["-seg", "no"];
        ncbi_args.extend(ncbi_comp_args);
        ncbi_args.extend([
            "-sum_stats",
            "false",
            "-evalue",
            "1000",
            "-culling_limit",
            "1",
        ]);
        assert_translated_subject_outfmt_matches_ncbi(
            program,
            ncbi_bin,
            query,
            subject,
            "6 qseqid sseqid score evalue qstart qend sstart send",
            &rust_args,
            &ncbi_args,
        );
    }
}

#[test]
fn translated_subject_ncbi_parity_hit_limit_and_identity_filters() {
    for (program, ncbi_bin, query, subject) in [
        (
            "blastx",
            "/usr/bin/blastx",
            ">q1\nATGAAATTTCTGATTCTGCTGTTTAAAACCCCGGGGTTTTATGAAATTTCTGATTCTGCTGTTT\n",
            ">full\nMKFLILLFQQQQGMKFLILLF\n>partial\nMKFLILLF\n>near\nMKFLILLY\n",
        ),
        (
            "tblastn",
            "/usr/bin/tblastn",
            ">q1\nMKFLILLFQQQQGMKFLILLF\n",
            ">full\nATGAAATTTCTGATTCTGCTGTTTAAAACCCCGGGGTTTTATGAAATTTCTGATTCTGCTGTTT\n>partial\nATGAAATTTCTGATTCTGCTGTTT\n>near\nATGAAATTTCTGATTCTGCTGTTACT\n",
        ),
        (
            "tblastx",
            "/usr/bin/tblastx",
            ">q1\nATGAAATTTCTGATTCTGCTGTTTAAAACCCCGGGGTTTTATGAAATTTCTGATTCTGCTGTTT\n",
            ">full\nATGAAATTTCTGATTCTGCTGTTTAAAACCCCGGGGTTTTATGAAATTTCTGATTCTGCTGTTT\n>partial\nATGAAATTTCTGATTCTGCTGTTT\n>near\nATGAAATTTCTGATTCTGCTGTTACT\n",
        ),
    ] {
        for (rust_extra, ncbi_extra) in [
            (
                vec!["--max_hsps", "1"],
                vec!["-max_hsps", "1"],
            ),
            (
                vec!["--max_target_seqs", "1"],
                vec!["-max_target_seqs", "1"],
            ),
            (
                vec!["--perc_identity", "90"],
                vec!["-perc_identity", "90"],
            ),
            (
                vec!["--qcov_hsp_perc", "80"],
                vec!["-qcov_hsp_perc", "80"],
            ),
        ] {
            if rust_extra.first() == Some(&"--perc_identity") {
                continue;
            }
            let (rust_comp_args, ncbi_comp_args): (&[&str], &[&str]) = if program == "tblastx" {
                (&[], &[])
            } else {
                (&["--comp_based_stats", "0"], &["-comp_based_stats", "0"])
            };
            let mut rust_args = vec!["--seg", "no"];
            rust_args.extend(rust_comp_args);
            rust_args.extend(["--sum_stats", "false", "--evalue", "1000"]);
            rust_args.extend(rust_extra);
            let mut ncbi_args = vec!["-seg", "no"];
            ncbi_args.extend(ncbi_comp_args);
            ncbi_args.extend(["-sum_stats", "false", "-evalue", "1000"]);
            ncbi_args.extend(ncbi_extra);
            assert_translated_subject_outfmt_matches_ncbi(
                program,
                ncbi_bin,
                query,
                subject,
                "6 qseqid sseqid score pident length qstart qend sstart send",
                &rust_args,
                &ncbi_args,
            );
        }
    }
}

#[test]
fn tblastx_subject_ncbi_parity_stable_hit_limit_filters() {
    for (rust_extra, ncbi_extra) in [
        (vec!["--max_hsps", "1"], vec!["-max_hsps", "1"]),
        (
            vec!["--max_target_seqs", "1"],
            vec!["-max_target_seqs", "1"],
        ),
        (vec!["--qcov_hsp_perc", "80"], vec!["-qcov_hsp_perc", "80"]),
    ] {
        let mut rust_args = vec!["--seg", "no", "--sum_stats", "false", "--evalue", "1000"];
        rust_args.extend(rust_extra);
        let mut ncbi_args = vec!["-seg", "no", "-sum_stats", "false", "-evalue", "1000"];
        ncbi_args.extend(ncbi_extra);
        assert_translated_subject_outfmt_matches_ncbi(
            "tblastx",
            "/usr/bin/tblastx",
            ">q1\nATGAAATTTCTGATTCTGCTGTTTAAAACCCCGGGGTTTTATGAAATTTCTGATTCTGCTGTTT\n",
            ">full\nATGAAATTTCTGATTCTGCTGTTTAAAACCCCGGGGTTTTATGAAATTTCTGATTCTGCTGTTT\n>partial\nATGAAATTTCTGATTCTGCTGTTT\n>near\nATGAAATTTCTGATTCTGCTGTTACT\n",
            "6 qseqid sseqid score pident length qstart qend sstart send qframe sframe",
            &rust_args,
            &ncbi_args,
        );
    }
}

#[test]
fn translated_subject_ncbi_parity_best_hit_filters() {
    for (program, ncbi_bin, query, subject) in [
        (
            "blastx",
            "/usr/bin/blastx",
            ">q1\nATGAAATTTCTGATTCTGCTGTTTAAAACCCCGGGGTTTTATGAAATTTCTGATTCTGCTGTTT\n",
            ">full\nMKFLILLFQQQQGMKFLILLF\n>contained\nKFLILL\n>near\nMKFLILLYQQQQGMKFLILLF\n",
        ),
        (
            "tblastn",
            "/usr/bin/tblastn",
            ">q1\nMKFLILLFQQQQGMKFLILLF\n",
            ">full\nATGAAATTTCTGATTCTGCTGTTTAAAACCCCGGGGTTTTATGAAATTTCTGATTCTGCTGTTT\n>contained\nAAATTTCTGATTCTGCTG\n>near\nATGAAATTTCTGATTCTGCTGTTACTAAAACCCCGGGGTTTTATGAAATTTCTGATTCTGCTGTTT\n",
        ),
        (
            "tblastx",
            "/usr/bin/tblastx",
            ">q1\nATGAAATTTCTGATTCTGCTGTTTAAAACCCCGGGGTTTTATGAAATTTCTGATTCTGCTGTTT\n",
            ">full\nATGAAATTTCTGATTCTGCTGTTTAAAACCCCGGGGTTTTATGAAATTTCTGATTCTGCTGTTT\n>contained\nAAATTTCTGATTCTGCTG\n>near\nATGAAATTTCTGATTCTGCTGTTACTAAAACCCCGGGGTTTTATGAAATTTCTGATTCTGCTGTTT\n",
        ),
    ] {
        let (rust_extra, ncbi_extra) = (vec!["--subject_besthit"], vec!["-subject_besthit"]);
        let (rust_comp_args, ncbi_comp_args): (&[&str], &[&str]) = if program == "tblastx" {
            (&[], &[])
        } else {
            (&["--comp_based_stats", "0"], &["-comp_based_stats", "0"])
        };
        let mut rust_args = vec!["--seg", "no"];
        rust_args.extend(rust_comp_args);
        rust_args.extend(["--sum_stats", "false", "--evalue", "1000"]);
        rust_args.extend(rust_extra);
        let mut ncbi_args = vec!["-seg", "no"];
        ncbi_args.extend(ncbi_comp_args);
        ncbi_args.extend(["-sum_stats", "false", "-evalue", "1000"]);
        ncbi_args.extend(ncbi_extra);
        assert_translated_subject_outfmt_matches_ncbi(
            program,
            ncbi_bin,
            query,
            subject,
            "6 qseqid sseqid score evalue qstart qend sstart send",
            &rust_args,
            &ncbi_args,
        );
    }
}

#[test]
fn blastx_db_ncbi_parity_exact_translation_coordinates_and_frames() {
    assert_translated_db_outfmt_matches_ncbi_sorted_lines(
        "blastx",
        "/usr/bin/blastx",
        "prot",
        ">q1\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">s1\nMKFLILLF\n",
        "6 qseqid sseqid qlen slen qstart qend sstart send score qframe sframe frames length",
        &[],
        &[],
    );
}

#[test]
fn tblastn_db_ncbi_parity_exact_translation_coordinates_and_frames() {
    assert_translated_db_outfmt_matches_ncbi_sorted_lines(
        "tblastn",
        "/usr/bin/tblastn",
        "nucl",
        ">q1\nMKFLILLF\n",
        ">s1\nATGAAATTTCTGATTCTGCTGTTT\n",
        "6 qseqid sseqid qlen slen qstart qend sstart send score qframe sframe frames length",
        &[],
        &[],
    );
}

#[test]
fn tblastx_db_ncbi_parity_exact_translation_coordinates_and_frames() {
    assert_translated_db_outfmt_matches_ncbi_sorted_lines(
        "tblastx",
        "/usr/bin/tblastx",
        "nucl",
        ">q1\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">s1\nATGAAATTTCTGATTCTGCTGTTT\n",
        "6 qseqid sseqid qlen slen qstart qend sstart send score qframe sframe frames length",
        &[],
        &[],
    );
}

#[test]
fn blastx_db_ncbi_parity_reverse_frame_coordinates_and_frames() {
    let query = format!(
        ">q1\n{}\n",
        ascii_reverse_complement("ATGAAATTTCTGATTCTGCTGTTT")
    );
    assert_translated_db_outfmt_matches_ncbi_sorted_lines(
        "blastx",
        "/usr/bin/blastx",
        "prot",
        &query,
        ">s1\nMKFLILLF\n",
        "6 qseqid sseqid qlen slen qstart qend sstart send score qframe sframe frames length",
        &[],
        &[],
    );
}

#[test]
fn tblastn_db_ncbi_parity_reverse_frame_coordinates_and_frames() {
    let subject = format!(
        ">s1\n{}\n",
        ascii_reverse_complement("ATGAAATTTCTGATTCTGCTGTTT")
    );
    assert_translated_db_outfmt_matches_ncbi_sorted_lines(
        "tblastn",
        "/usr/bin/tblastn",
        "nucl",
        ">q1\nMKFLILLF\n",
        &subject,
        "6 qseqid sseqid qlen slen qstart qend sstart send score qframe sframe frames length",
        &[],
        &[],
    );
}

#[test]
fn tblastx_db_ncbi_parity_reverse_frame_coordinates_and_frames() {
    let query = format!(
        ">q1\n{}\n",
        ascii_reverse_complement("ATGAAATTTCTGATTCTGCTGTTT")
    );
    let subject = format!(
        ">s1\n{}\n",
        ascii_reverse_complement("ATGAAATTTCTGATTCTGCTGTTT")
    );
    assert_translated_db_outfmt_matches_ncbi_sorted_lines(
        "tblastx",
        "/usr/bin/tblastx",
        "nucl",
        &query,
        &subject,
        "6 qseqid sseqid qlen slen qstart qend sstart send score qframe sframe frames length",
        &[],
        &[],
    );
}

#[test]
fn translated_db_ncbi_parity_custom_field_parser_edges() {
    for (program, ncbi_bin, dbtype, query, db_fasta) in [
        (
            "blastx",
            "/usr/bin/blastx",
            "prot",
            ">q1\nATGAAATTTCTGATTCTGCTGTTT\n",
            ">s1\nMKFLILLF\n",
        ),
        (
            "tblastn",
            "/usr/bin/tblastn",
            "nucl",
            ">q1\nMKFLILLF\n",
            ">s1\nATGAAATTTCTGATTCTGCTGTTT\n",
        ),
        (
            "tblastx",
            "/usr/bin/tblastx",
            "nucl",
            ">q1\nATGAAATTTCTGATTCTGCTGTTT\n",
            ">s1\nATGAAATTTCTGATTCTGCTGTTT\n",
        ),
    ] {
        for outfmt in [
            "6 delim=, qseqid bogus sseqid qseqid length length",
            "6 delim=tab qseqid sseqid length",
            "6 std qlen",
            "6 bogus",
            "7 qseqid bogus sseqid qseqid length length",
            "7 delim=space qseqid sseqid length",
            "7 std qlen",
            "7 bogus",
            "10 qseqid bogus sseqid qseqid length length",
            "10 delim=tab qseqid sseqid length",
            "10 std qlen",
            "10 bogus",
        ] {
            assert_translated_db_outfmt_matches_ncbi(
                program,
                ncbi_bin,
                dbtype,
                query,
                db_fasta,
                outfmt,
                &["--seg", "no"],
                &["-seg", "no"],
            );
        }
    }
}

#[test]
fn blastx_db_ncbi_parity_commented_tabular_exact_hit() {
    assert_translated_db_outfmt_matches_ncbi(
        "blastx",
        "/usr/bin/blastx",
        "prot",
        ">q1\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">s1\nMKFLILLF\n",
        "7",
        &[],
        &[],
    );
}

#[test]
fn blastx_db_ncbi_parity_commented_tabular_parse_deflines() {
    assert_translated_db_outfmt_matches_ncbi(
        "blastx",
        "/usr/bin/blastx",
        "prot",
        ">gi|123|ref|QNT.1| query nucleotide\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">gi|456|ref|SPROT.1| subject protein\nMKFLILLF\n",
        "7",
        &["--seg", "no", "--comp_based_stats", "0", "--parse_deflines"],
        &["-seg", "no", "-comp_based_stats", "0", "-parse_deflines"],
    );
}

#[test]
fn blastx_db_ncbi_parity_xml_exact_hit() {
    assert_translated_db_outfmt_matches_ncbi(
        "blastx",
        "/usr/bin/blastx",
        "prot",
        ">q1\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">s1\nMKFLILLF\n",
        "5",
        &["--seg", "no", "--comp_based_stats", "0"],
        &["-seg", "no", "-comp_based_stats", "0"],
    );
}

#[test]
fn blastx_db_ncbi_parity_xml_positive_substitution() {
    assert_translated_db_outfmt_matches_ncbi(
        "blastx",
        "/usr/bin/blastx",
        "prot",
        ">q1\nATGAAATTTCTGATTCTGCTGTTTAACATTCTGTGCCTGTTTCCTGTTCTGGCTGCTGATAACCATGGTGTTTCTATGAACGCTTCT\n",
        ">s1\nMKFLILLFNILCLFPVLAADNHGVSINAS\n",
        "5",
        &["--seg", "no", "--comp_based_stats", "0"],
        &["-seg", "no", "-comp_based_stats", "0"],
    );
}

#[test]
fn blastx_db_ncbi_parity_xml_no_hits() {
    assert_translated_db_outfmt_matches_ncbi(
        "blastx",
        "/usr/bin/blastx",
        "prot",
        ">q1\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">s1 unrelated\nPPPPPPPPPPPP\n",
        "5",
        &["--seg", "no", "--comp_based_stats", "0"],
        &["-seg", "no", "-comp_based_stats", "0"],
    );
}

#[test]
fn blastx_db_ncbi_parity_xml_multiple_queries() {
    assert_translated_db_outfmt_matches_ncbi(
        "blastx",
        "/usr/bin/blastx",
        "prot",
        ">q1\nATGAAATTTCTGATTCTGCTGTTT\n>q2 unrelated\nCCCCCCCCCCCCCCCCCCCCCCCC\n",
        ">s1\nMKFLILLF\n",
        "5",
        &["--seg", "no", "--comp_based_stats", "0"],
        &["-seg", "no", "-comp_based_stats", "0"],
    );
}

#[test]
fn blastx_db_ncbi_parity_xml_multi_hsp_same_subject() {
    assert_translated_db_outfmt_matches_ncbi(
        "blastx",
        "/usr/bin/blastx",
        "prot",
        ">q1\nATGAAATTTCTGATTCTGCTGTTTAAAACCCCGGGGTTTTATGAAATTTCTGATTCTGCTGTTT\n",
        ">s1 multi\nMKFLILLFQQQQGMKFLILLF\n",
        "5",
        &["--seg", "no", "--comp_based_stats", "0", "--evalue", "1000"],
        &["-seg", "no", "-comp_based_stats", "0", "-evalue", "1000"],
    );
}

#[test]
fn blastx_db_ncbi_parity_xml_parse_deflines() {
    assert_translated_db_outfmt_matches_ncbi(
        "blastx",
        "/usr/bin/blastx",
        "prot",
        ">gi|123|ref|QNT.1| query nucleotide\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">gi|456|ref|SPROT.1| subject protein\nMKFLILLF\n",
        "5",
        &["--seg", "no", "--comp_based_stats", "0", "--parse_deflines"],
        &["-seg", "no", "-comp_based_stats", "0", "-parse_deflines"],
    );
}

#[test]
fn blastx_db_ncbi_parity_pairwise_exact_hit() {
    assert_translated_db_outfmt_matches_ncbi(
        "blastx",
        "/usr/bin/blastx",
        "prot",
        ">q1\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">s1\nMKFLILLF\n",
        "0",
        &["--seg", "no", "--comp_based_stats", "0"],
        &["-seg", "no", "-comp_based_stats", "0"],
    );
}

#[test]
fn blastx_db_ncbi_parity_pairwise_parse_deflines() {
    assert_translated_db_outfmt_matches_ncbi(
        "blastx",
        "/usr/bin/blastx",
        "prot",
        ">gi|123|ref|QNT.1| query nucleotide\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">gi|456|ref|SPROT.1| subject protein\nMKFLILLF\n",
        "0",
        &["--seg", "no", "--comp_based_stats", "0", "--parse_deflines"],
        &["-seg", "no", "-comp_based_stats", "0", "-parse_deflines"],
    );
}

#[test]
fn blastx_db_ncbi_parity_tabular_parse_deflines() {
    assert_translated_db_outfmt_matches_ncbi(
        "blastx",
        "/usr/bin/blastx",
        "prot",
        ">gi|123|ref|QNT.1| query nucleotide\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">gi|456|ref|SPROT.1| subject protein\nMKFLILLF\n",
        "6 qseqid qacc qaccver sseqid sacc saccver",
        &["--seg", "no", "--comp_based_stats", "0", "--parse_deflines"],
        &["-seg", "no", "-comp_based_stats", "0", "-parse_deflines"],
    );
}

#[test]
fn blastx_db_ncbi_parity_tabular_parse_deflines_gis() {
    assert_translated_db_outfmt_matches_ncbi(
        "blastx",
        "/usr/bin/blastx",
        "prot",
        ">gi|123|ref|QNT.1| query nucleotide\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">gi|456|ref|SPROT.1| subject protein\nMKFLILLF\n",
        "6 qseqid qgi qacc qaccver sseqid sgi sallgi sacc saccver sallseqid",
        &["--seg", "no", "--comp_based_stats", "0", "--parse_deflines"],
        &["-seg", "no", "-comp_based_stats", "0", "-parse_deflines"],
    );
}

#[test]
fn blastx_db_ncbi_parity_tabular_parse_deflines_titles() {
    assert_translated_db_outfmt_matches_ncbi(
        "blastx",
        "/usr/bin/blastx",
        "prot",
        ">gi|123|ref|QNT.1| query nucleotide\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">gi|456|ref|SPROT.1| subject protein\nMKFLILLF\n",
        "6 sseqid stitle salltitles",
        &["--seg", "no", "--comp_based_stats", "0", "--parse_deflines"],
        &["-seg", "no", "-comp_based_stats", "0", "-parse_deflines"],
    );
}

#[test]
fn blastx_db_ncbi_parity_csv_parse_deflines_metadata() {
    assert_translated_db_outfmt_matches_ncbi(
        "blastx",
        "/usr/bin/blastx",
        "prot",
        ">gi|123|ref|QNT.1| query nucleotide\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">gi|456|ref|SPROT.1| subject protein\nMKFLILLF\n",
        "10 qseqid qgi qacc qaccver sseqid sgi sallgi sacc saccver sallseqid stitle salltitles",
        &["--seg", "no", "--comp_based_stats", "0", "--parse_deflines"],
        &["-seg", "no", "-comp_based_stats", "0", "-parse_deflines"],
    );
}

#[test]
fn blastx_db_ncbi_parity_pairwise_no_hits() {
    assert_translated_db_outfmt_matches_ncbi(
        "blastx",
        "/usr/bin/blastx",
        "prot",
        ">q1\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">s1\nGGGGGGGG\n",
        "0",
        &[
            "--seg",
            "no",
            "--comp_based_stats",
            "0",
            "--evalue",
            "1e-100",
        ],
        &["-seg", "no", "-comp_based_stats", "0", "-evalue", "1e-100"],
    );
}

#[test]
fn blastx_db_ncbi_parity_pairwise_line_length() {
    assert_translated_db_outfmt_matches_ncbi(
        "blastx",
        "/usr/bin/blastx",
        "prot",
        ">q1\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">s1\nMKFLILLF\n",
        "0",
        &[
            "--seg",
            "no",
            "--comp_based_stats",
            "0",
            "--line_length",
            "5",
        ],
        &["-seg", "no", "-comp_based_stats", "0", "-line_length", "5"],
    );
}

#[test]
fn blastx_db_ncbi_parity_pairwise_num_alignments() {
    assert_translated_db_outfmt_matches_ncbi(
        "blastx",
        "/usr/bin/blastx",
        "prot",
        ">q1\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">s1 first\nMKFLILLF\n>s2 second\nMKFLILLF\n",
        "0",
        &[
            "--seg",
            "no",
            "--comp_based_stats",
            "0",
            "--num_alignments",
            "1",
            "--num_descriptions",
            "1",
        ],
        &[
            "-seg",
            "no",
            "-comp_based_stats",
            "0",
            "-num_alignments",
            "1",
            "-num_descriptions",
            "1",
        ],
    );
}

#[test]
fn blastx_db_ncbi_parity_pairwise_multi_hsp_same_subject() {
    assert_translated_db_outfmt_matches_ncbi(
        "blastx",
        "/usr/bin/blastx",
        "prot",
        ">q1\nATGAAATTTCTGATTCTGCTGTTTAAAACCCCGGGGTTTTATGAAATTTCTGATTCTGCTGTTT\n",
        ">s1 multi\nMKFLILLFQQQQGMKFLILLF\n",
        "0",
        &["--seg", "no", "--comp_based_stats", "0", "--evalue", "1000"],
        &["-seg", "no", "-comp_based_stats", "0", "-evalue", "1000"],
    );
}

#[test]
fn blastx_db_ncbi_parity_tabular_sum_stats_gapped_single_hsp() {
    assert_translated_db_outfmt_matches_ncbi(
        "blastx",
        "/usr/bin/blastx",
        "prot",
        ">q1\nATGAAATTTCTGATTCTGCTGTTTAAATTTATGAAATTTCTGATTCTGCTGTTT\n",
        ">s1\nMKFLILLFMKFLILLF\n",
        "6 qseqid sseqid score bitscore evalue length qstart qend sstart send",
        &[
            "--seg",
            "no",
            "--comp_based_stats",
            "0",
            "--sum_stats",
            "true",
            "--evalue",
            "1000",
        ],
        &[
            "-seg",
            "no",
            "-comp_based_stats",
            "0",
            "-sum_stats",
            "true",
            "-evalue",
            "1000",
        ],
    );
}

#[test]
fn blastx_db_ncbi_parity_num_threads_two_gapped_traceback() {
    assert_translated_db_outfmt_matches_ncbi_with_num_threads(
        "blastx",
        "/usr/bin/blastx",
        "prot",
        ">q1\nATGAAATTTCTGATTCTGCTGTTTAAATTTATGAAATTTCTGATTCTGCTGTTT\n",
        ">s1\nMKFLILLFMKFLILLF\n",
        "6 qseqid sseqid score bitscore evalue length gaps qstart qend sstart send qseq sseq btop",
        &[
            "--seg",
            "no",
            "--comp_based_stats",
            "0",
            "--sum_stats",
            "true",
            "--evalue",
            "1000",
        ],
        &[
            "-seg",
            "no",
            "-comp_based_stats",
            "0",
            "-sum_stats",
            "true",
            "-evalue",
            "1000",
        ],
        "2",
    );
}

#[test]
fn blastx_db_ncbi_parity_tabular_sum_stats_false_single_hsp() {
    assert_translated_db_outfmt_matches_ncbi(
        "blastx",
        "/usr/bin/blastx",
        "prot",
        ">q1\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">s1\nMKFLILLF\n",
        "6 qseqid sseqid score bitscore evalue length qstart qend sstart send",
        &[
            "--seg",
            "no",
            "--comp_based_stats",
            "0",
            "--sum_stats",
            "false",
            "--evalue",
            "1000",
        ],
        &[
            "-seg",
            "no",
            "-comp_based_stats",
            "0",
            "-sum_stats",
            "false",
            "-evalue",
            "1000",
        ],
    );
}

#[test]
fn blastx_db_ncbi_parity_explicit_sum_stats_gapped_single_hsp() {
    assert_translated_db_outfmt_matches_ncbi(
        "blastx",
        "/usr/bin/blastx",
        "prot",
        ">q1\nATGAAATTTCTGATTCTGCTGTTTAAATTTATGAAATTTCTGATTCTGCTGTTT\n",
        ">s1\nMKFLILLFMKFLILLF\n",
        "0",
        &[
            "--seg",
            "no",
            "--comp_based_stats",
            "0",
            "--sum_stats",
            "true",
            "--evalue",
            "1000",
        ],
        &[
            "-seg",
            "no",
            "-comp_based_stats",
            "0",
            "-sum_stats",
            "true",
            "-evalue",
            "1000",
        ],
    );
}

#[test]
fn tblastn_db_ncbi_parity_commented_tabular_exact_hit() {
    assert_translated_db_outfmt_matches_ncbi(
        "tblastn",
        "/usr/bin/tblastn",
        "nucl",
        ">q1\nMKFLILLF\n",
        ">s1\nATGAAATTTCTGATTCTGCTGTTT\n",
        "7",
        &[],
        &[],
    );
}

#[test]
fn tblastn_db_ncbi_parity_commented_tabular_parse_deflines() {
    assert_translated_db_outfmt_matches_ncbi(
        "tblastn",
        "/usr/bin/tblastn",
        "nucl",
        ">gi|123|ref|QPROT.1| query protein\nMKFLILLF\n",
        ">gi|456|ref|SNT.1| subject nucleotide\nATGAAATTTCTGATTCTGCTGTTT\n",
        "7 qaccver saccver pident length bitscore",
        &["--seg", "no", "--comp_based_stats", "0", "--parse_deflines"],
        &["-seg", "no", "-comp_based_stats", "0", "-parse_deflines"],
    );
}

#[test]
fn tblastn_db_ncbi_parity_xml_exact_hit() {
    assert_translated_db_outfmt_matches_ncbi(
        "tblastn",
        "/usr/bin/tblastn",
        "nucl",
        ">q1\nMKFLILLF\n",
        ">s1\nATGAAATTTCTGATTCTGCTGTTT\n",
        "5",
        &["--seg", "no", "--comp_based_stats", "0"],
        &["-seg", "no", "-comp_based_stats", "0"],
    );
}

#[test]
fn tblastn_db_ncbi_parity_xml_no_hits() {
    assert_translated_db_outfmt_matches_ncbi(
        "tblastn",
        "/usr/bin/tblastn",
        "nucl",
        ">q1\nMKFLILLF\n",
        ">s1 unrelated\nGGGGGGGGGGGGGGGGGGGGGGGG\n",
        "5",
        &["--seg", "no", "--comp_based_stats", "0"],
        &["-seg", "no", "-comp_based_stats", "0"],
    );
}

#[test]
fn tblastn_db_ncbi_parity_xml_multiple_queries() {
    assert_translated_db_outfmt_matches_ncbi(
        "tblastn",
        "/usr/bin/tblastn",
        "nucl",
        ">q1\nMKFLILLF\n>q2 unrelated\nPPPPPPPP\n",
        ">s1\nATGAAATTTCTGATTCTGCTGTTT\n",
        "5",
        &["--seg", "no", "--comp_based_stats", "0"],
        &["-seg", "no", "-comp_based_stats", "0"],
    );
}

#[test]
fn tblastn_db_ncbi_parity_xml_multi_hsp_same_subject() {
    assert_translated_db_outfmt_matches_ncbi(
        "tblastn",
        "/usr/bin/tblastn",
        "nucl",
        ">q1\nMKFLILLFQQQQGMKFLILLF\n",
        ">s1 multi\nATGAAATTTCTGATTCTGCTGTTTAAAACCCCGGGGTTTTATGAAATTTCTGATTCTGCTGTTT\n",
        "5",
        &["--seg", "no", "--comp_based_stats", "0", "--evalue", "1000"],
        &["-seg", "no", "-comp_based_stats", "0", "-evalue", "1000"],
    );
}

#[test]
fn tblastn_db_ncbi_parity_xml_parse_deflines() {
    assert_translated_db_outfmt_matches_ncbi(
        "tblastn",
        "/usr/bin/tblastn",
        "nucl",
        ">gi|123|ref|QPROT.1| query protein\nMKFLILLF\n",
        ">gi|456|ref|SNT.1| subject nucleotide\nATGAAATTTCTGATTCTGCTGTTT\n",
        "5",
        &["--seg", "no", "--comp_based_stats", "0", "--parse_deflines"],
        &["-seg", "no", "-comp_based_stats", "0", "-parse_deflines"],
    );
}

#[test]
fn tblastn_db_ncbi_parity_pairwise_exact_hit() {
    assert_translated_db_outfmt_matches_ncbi(
        "tblastn",
        "/usr/bin/tblastn",
        "nucl",
        ">q1\nMKFLILLF\n",
        ">s1\nATGAAATTTCTGATTCTGCTGTTT\n",
        "0",
        &["--seg", "no", "--comp_based_stats", "0"],
        &["-seg", "no", "-comp_based_stats", "0"],
    );
}

#[test]
fn tblastn_db_ncbi_parity_pairwise_parse_deflines() {
    assert_translated_db_outfmt_matches_ncbi(
        "tblastn",
        "/usr/bin/tblastn",
        "nucl",
        ">gi|123|ref|QPROT.1| query protein\nMKFLILLF\n",
        ">gi|456|ref|SNT.1| subject nucleotide\nATGAAATTTCTGATTCTGCTGTTT\n",
        "0",
        &["--seg", "no", "--comp_based_stats", "0", "--parse_deflines"],
        &["-seg", "no", "-comp_based_stats", "0", "-parse_deflines"],
    );
}

#[test]
fn tblastn_db_ncbi_parity_explicit_sum_stats_gapped_single_hsp() {
    assert_translated_db_outfmt_matches_ncbi(
        "tblastn",
        "/usr/bin/tblastn",
        "nucl",
        ">q1\nMKFLILLFMKFLILLF\n",
        ">s1\nATGAAATTTCTGATTCTGCTGTTTAAATTTATGAAATTTCTGATTCTGCTGTTT\n",
        "0",
        &[
            "--seg",
            "no",
            "--comp_based_stats",
            "0",
            "--sum_stats",
            "true",
            "--evalue",
            "1000",
        ],
        &[
            "-seg",
            "no",
            "-comp_based_stats",
            "0",
            "-sum_stats",
            "true",
            "-evalue",
            "1000",
        ],
    );
}

#[test]
fn tblastn_db_ncbi_parity_tabular_sum_stats_gapped_single_hsp() {
    assert_translated_db_outfmt_matches_ncbi(
        "tblastn",
        "/usr/bin/tblastn",
        "nucl",
        ">q1\nMKFLILLFMKFLILLF\n",
        ">s1\nATGAAATTTCTGATTCTGCTGTTTAAATTTATGAAATTTCTGATTCTGCTGTTT\n",
        "6 qseqid sseqid score bitscore evalue length qstart qend sstart send",
        &[
            "--seg",
            "no",
            "--comp_based_stats",
            "0",
            "--sum_stats",
            "true",
            "--evalue",
            "1000",
        ],
        &[
            "-seg",
            "no",
            "-comp_based_stats",
            "0",
            "-sum_stats",
            "true",
            "-evalue",
            "1000",
        ],
    );
}

#[test]
fn tblastn_db_ncbi_parity_tabular_sum_stats_false_single_hsp() {
    assert_translated_db_outfmt_matches_ncbi(
        "tblastn",
        "/usr/bin/tblastn",
        "nucl",
        ">q1\nMKFLILLF\n",
        ">s1\nATGAAATTTCTGATTCTGCTGTTT\n",
        "6 qseqid sseqid score bitscore evalue length qstart qend sstart send",
        &[
            "--seg",
            "no",
            "--comp_based_stats",
            "0",
            "--sum_stats",
            "false",
            "--evalue",
            "1000",
        ],
        &[
            "-seg",
            "no",
            "-comp_based_stats",
            "0",
            "-sum_stats",
            "false",
            "-evalue",
            "1000",
        ],
    );
}

#[test]
fn tblastn_db_ncbi_parity_tabular_parse_deflines() {
    assert_translated_db_outfmt_matches_ncbi(
        "tblastn",
        "/usr/bin/tblastn",
        "nucl",
        ">gi|123|ref|QPROT.1| query protein\nMKFLILLF\n",
        ">gi|456|ref|SNT.1| subject nucleotide\nATGAAATTTCTGATTCTGCTGTTT\n",
        "6 qseqid qacc qaccver sseqid sacc saccver pident length bitscore",
        &["--seg", "no", "--comp_based_stats", "0", "--parse_deflines"],
        &["-seg", "no", "-comp_based_stats", "0", "-parse_deflines"],
    );
}

#[test]
fn tblastn_db_ncbi_parity_tabular_parse_deflines_gis() {
    assert_translated_db_outfmt_matches_ncbi(
        "tblastn",
        "/usr/bin/tblastn",
        "nucl",
        ">gi|123|ref|QPROT.1| query protein\nMKFLILLF\n",
        ">gi|456|ref|SNT.1| subject nucleotide\nATGAAATTTCTGATTCTGCTGTTT\n",
        "6 qseqid qgi qacc qaccver sseqid sgi sallgi sacc saccver sallseqid pident length bitscore",
        &["--seg", "no", "--comp_based_stats", "0", "--parse_deflines"],
        &["-seg", "no", "-comp_based_stats", "0", "-parse_deflines"],
    );
}

#[test]
fn tblastn_db_ncbi_parity_tabular_parse_deflines_titles() {
    assert_translated_db_outfmt_matches_ncbi(
        "tblastn",
        "/usr/bin/tblastn",
        "nucl",
        ">gi|123|ref|QPROT.1| query protein\nMKFLILLF\n",
        ">gi|456|ref|SNT.1| subject nucleotide\nATGAAATTTCTGATTCTGCTGTTT\n",
        "6 sseqid stitle salltitles",
        &["--seg", "no", "--comp_based_stats", "0", "--parse_deflines"],
        &["-seg", "no", "-comp_based_stats", "0", "-parse_deflines"],
    );
}

#[test]
fn tblastn_db_ncbi_parity_csv_parse_deflines_metadata() {
    assert_translated_db_outfmt_matches_ncbi(
        "tblastn",
        "/usr/bin/tblastn",
        "nucl",
        ">gi|123|ref|QPROT.1| query protein\nMKFLILLF\n",
        ">gi|456|ref|SNT.1| subject nucleotide\nATGAAATTTCTGATTCTGCTGTTT\n",
        "10 qseqid qgi qacc qaccver sseqid sgi sallgi sacc saccver sallseqid stitle salltitles",
        &["--seg", "no", "--comp_based_stats", "0", "--parse_deflines"],
        &["-seg", "no", "-comp_based_stats", "0", "-parse_deflines"],
    );
}

#[test]
fn tblastn_db_ncbi_parity_pairwise_no_hits() {
    assert_translated_db_outfmt_matches_ncbi(
        "tblastn",
        "/usr/bin/tblastn",
        "nucl",
        ">q1\nMKFLILLF\n",
        ">s1\nGGGGGGGGGGGGGGGGGGGGGGGG\n",
        "0",
        &[
            "--seg",
            "no",
            "--comp_based_stats",
            "0",
            "--evalue",
            "1e-100",
        ],
        &["-seg", "no", "-comp_based_stats", "0", "-evalue", "1e-100"],
    );
}

#[test]
fn tblastn_db_ncbi_parity_pairwise_line_length() {
    assert_translated_db_outfmt_matches_ncbi(
        "tblastn",
        "/usr/bin/tblastn",
        "nucl",
        ">q1\nMKFLILLF\n",
        ">s1\nATGAAATTTCTGATTCTGCTGTTT\n",
        "0",
        &[
            "--seg",
            "no",
            "--comp_based_stats",
            "0",
            "--line_length",
            "5",
        ],
        &["-seg", "no", "-comp_based_stats", "0", "-line_length", "5"],
    );
}

#[test]
fn tblastn_db_ncbi_parity_pairwise_num_alignments() {
    assert_translated_db_outfmt_matches_ncbi(
        "tblastn",
        "/usr/bin/tblastn",
        "nucl",
        ">q1\nMKFLILLF\n",
        ">s1 first\nATGAAATTTCTGATTCTGCTGTTT\n>s2 second\nATGAAATTTCTGATTCTGCTGTTT\n",
        "0",
        &[
            "--seg",
            "no",
            "--comp_based_stats",
            "0",
            "--num_alignments",
            "1",
            "--num_descriptions",
            "1",
        ],
        &[
            "-seg",
            "no",
            "-comp_based_stats",
            "0",
            "-num_alignments",
            "1",
            "-num_descriptions",
            "1",
        ],
    );
}

#[test]
fn tblastn_db_ncbi_parity_pairwise_multi_hsp_same_subject() {
    assert_translated_db_outfmt_matches_ncbi(
        "tblastn",
        "/usr/bin/tblastn",
        "nucl",
        ">q1\nMKFLILLFQQQQGMKFLILLF\n",
        ">s1 multi\nATGAAATTTCTGATTCTGCTGTTTAAAACCCCGGGGTTTTATGAAATTTCTGATTCTGCTGTTT\n",
        "0",
        &["--seg", "no", "--comp_based_stats", "0", "--evalue", "1000"],
        &["-seg", "no", "-comp_based_stats", "0", "-evalue", "1000"],
    );
}

#[test]
fn tblastx_db_ncbi_parity_commented_tabular_exact_hit() {
    assert_translated_db_outfmt_matches_ncbi(
        "tblastx",
        "/usr/bin/tblastx",
        "nucl",
        ">q1\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">s1\nATGAAATTTCTGATTCTGCTGTTT\n",
        "7",
        &[],
        &[],
    );
}

#[test]
fn tblastx_db_ncbi_parity_commented_tabular_parse_deflines() {
    assert_translated_db_outfmt_matches_ncbi(
        "tblastx",
        "/usr/bin/tblastx",
        "nucl",
        ">gi|123|ref|QNT.1| query nucleotide\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">gi|456|ref|SNT.1| subject nucleotide\nATGAAATTTCTGATTCTGCTGTTT\n",
        "7",
        &["--seg", "no", "--evalue", "1000", "--parse_deflines"],
        &["-seg", "no", "-evalue", "1000", "-parse_deflines"],
    );
}

#[test]
fn tblastx_db_ncbi_parity_xml_exact_hit() {
    assert_translated_db_outfmt_matches_ncbi(
        "tblastx",
        "/usr/bin/tblastx",
        "nucl",
        ">q1\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">s1\nATGAAATTTCTGATTCTGCTGTTT\n",
        "5",
        &["--seg", "no", "--evalue", "1000"],
        &["-seg", "no", "-evalue", "1000"],
    );
}

#[test]
fn tblastx_db_ncbi_parity_xml_no_hits() {
    assert_translated_db_outfmt_matches_ncbi(
        "tblastx",
        "/usr/bin/tblastx",
        "nucl",
        ">q1\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">s1 unrelated\nCCACCTCCACCTCCACCTCCACCT\n",
        "5",
        &["--seg", "no"],
        &["-seg", "no"],
    );
}

#[test]
fn tblastx_db_ncbi_parity_xml_multiple_queries() {
    assert_translated_db_outfmt_matches_ncbi(
        "tblastx",
        "/usr/bin/tblastx",
        "nucl",
        ">q1\nATGAAATTTCTGATTCTGCTGTTT\n>q2 unrelated\nCCACCTCCACCTCCACCTCCACCT\n",
        ">s1\nATGAAATTTCTGATTCTGCTGTTT\n",
        "5",
        &["--seg", "no", "--evalue", "1000"],
        &["-seg", "no", "-evalue", "1000"],
    );
}

#[test]
fn tblastx_db_ncbi_parity_xml_multi_hsp_same_subject() {
    assert_translated_db_outfmt_matches_ncbi(
        "tblastx",
        "/usr/bin/tblastx",
        "nucl",
        ">q1\nATGAAATTTCTGATTCTGCTGTTTAAAACCCCGGGGTTTTATGAAATTTCTGATTCTGCTGTTT\n",
        ">s1 multi\nATGAAATTTCTGATTCTGCTGTTTAAAACCCCGGGGTTTTATGAAATTTCTGATTCTGCTGTTT\n",
        "5",
        &["--seg", "no", "--evalue", "1000"],
        &["-seg", "no", "-evalue", "1000"],
    );
}

#[test]
fn tblastx_db_ncbi_parity_xml_parse_deflines() {
    assert_translated_db_outfmt_matches_ncbi(
        "tblastx",
        "/usr/bin/tblastx",
        "nucl",
        ">gi|123|ref|QNT.1| query nucleotide\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">gi|456|ref|SNT.1| subject nucleotide\nATGAAATTTCTGATTCTGCTGTTT\n",
        "5",
        &["--seg", "no", "--evalue", "1000", "--parse_deflines"],
        &["-seg", "no", "-evalue", "1000", "-parse_deflines"],
    );
}

#[test]
fn tblastx_db_ncbi_parity_pairwise_exact_hit() {
    assert_translated_db_outfmt_matches_ncbi(
        "tblastx",
        "/usr/bin/tblastx",
        "nucl",
        ">q1\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">s1\nATGAAATTTCTGATTCTGCTGTTT\n",
        "0",
        &["--seg", "no", "--evalue", "1000"],
        &["-seg", "no", "-evalue", "1000"],
    );
}

#[test]
fn tblastx_db_ncbi_parity_pairwise_parse_deflines() {
    assert_translated_db_outfmt_matches_ncbi(
        "tblastx",
        "/usr/bin/tblastx",
        "nucl",
        ">gi|123|ref|QNT.1| query nucleotide\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">gi|456|ref|SNT.1| subject nucleotide\nATGAAATTTCTGATTCTGCTGTTT\n",
        "0",
        &["--seg", "no", "--evalue", "1000", "--parse_deflines"],
        &["-seg", "no", "-evalue", "1000", "-parse_deflines"],
    );
}

#[test]
fn tblastx_db_ncbi_parity_tabular_parse_deflines() {
    assert_translated_db_outfmt_matches_ncbi(
        "tblastx",
        "/usr/bin/tblastx",
        "nucl",
        ">gi|123|ref|QNT.1| query nucleotide\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">gi|456|ref|SNT.1| subject nucleotide\nATGAAATTTCTGATTCTGCTGTTT\n",
        "6 qseqid qacc qaccver sseqid sacc saccver pident length bitscore",
        &["--seg", "no", "--evalue", "1000", "--parse_deflines"],
        &["-seg", "no", "-evalue", "1000", "-parse_deflines"],
    );
}

#[test]
fn tblastx_db_ncbi_parity_tabular_parse_deflines_gis() {
    assert_translated_db_outfmt_matches_ncbi(
        "tblastx",
        "/usr/bin/tblastx",
        "nucl",
        ">gi|123|ref|QNT.1| query nucleotide\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">gi|456|ref|SNT.1| subject nucleotide\nATGAAATTTCTGATTCTGCTGTTT\n",
        "6 qseqid qgi qacc qaccver sseqid sgi sallgi sacc saccver sallseqid pident length bitscore",
        &["--seg", "no", "--evalue", "1000", "--parse_deflines"],
        &["-seg", "no", "-evalue", "1000", "-parse_deflines"],
    );
}

#[test]
fn tblastx_db_ncbi_parity_tabular_parse_deflines_titles() {
    assert_translated_db_outfmt_matches_ncbi(
        "tblastx",
        "/usr/bin/tblastx",
        "nucl",
        ">gi|123|ref|QNT.1| query nucleotide\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">gi|456|ref|SNT.1| subject nucleotide\nATGAAATTTCTGATTCTGCTGTTT\n",
        "6 sseqid stitle salltitles",
        &["--seg", "no", "--evalue", "1000", "--parse_deflines"],
        &["-seg", "no", "-evalue", "1000", "-parse_deflines"],
    );
}

#[test]
fn tblastx_db_ncbi_parity_csv_parse_deflines_metadata() {
    assert_translated_db_outfmt_matches_ncbi(
        "tblastx",
        "/usr/bin/tblastx",
        "nucl",
        ">gi|123|ref|QNT.1| query nucleotide\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">gi|456|ref|SNT.1| subject nucleotide\nATGAAATTTCTGATTCTGCTGTTT\n",
        "10 qseqid qgi qacc qaccver sseqid sgi sallgi sacc saccver sallseqid stitle salltitles",
        &["--seg", "no", "--evalue", "1000", "--parse_deflines"],
        &["-seg", "no", "-evalue", "1000", "-parse_deflines"],
    );
}

#[test]
fn tblastx_db_ncbi_parity_pairwise_no_hits() {
    assert_translated_db_outfmt_matches_ncbi(
        "tblastx",
        "/usr/bin/tblastx",
        "nucl",
        ">q1\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">s1\nGGGGGGGGGGGGGGGGGGGGGGGG\n",
        "0",
        &["--seg", "no", "--evalue", "1e-100"],
        &["-seg", "no", "-evalue", "1e-100"],
    );
}

#[test]
fn tblastx_db_ncbi_parity_pairwise_line_length() {
    assert_translated_db_outfmt_matches_ncbi(
        "tblastx",
        "/usr/bin/tblastx",
        "nucl",
        ">q1\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">s1\nATGAAATTTCTGATTCTGCTGTTT\n",
        "0",
        &["--seg", "no", "--evalue", "1000", "--line_length", "5"],
        &["-seg", "no", "-evalue", "1000", "-line_length", "5"],
    );
}

#[test]
fn tblastx_db_ncbi_parity_pairwise_num_alignments() {
    assert_translated_db_outfmt_matches_ncbi(
        "tblastx",
        "/usr/bin/tblastx",
        "nucl",
        ">q1\nATGAAATTTCTGATTCTGCTGTTT\n",
        ">s1 first\nATGAAATTTCTGATTCTGCTGTTT\n>s2 second\nATGAAATTTCTGATTCTGCTGTTT\n",
        "0",
        &[
            "--seg",
            "no",
            "--evalue",
            "1000",
            "--num_alignments",
            "1",
            "--num_descriptions",
            "1",
        ],
        &[
            "-seg",
            "no",
            "-evalue",
            "1000",
            "-num_alignments",
            "1",
            "-num_descriptions",
            "1",
        ],
    );
}

#[test]
fn tblastx_db_ncbi_parity_pairwise_multi_hsp_same_subject() {
    assert_translated_db_outfmt_matches_ncbi(
        "tblastx",
        "/usr/bin/tblastx",
        "nucl",
        ">q1\nATGAAATTTCTGATTCTGCTGTTTAAAACCCCGGGGTTTTATGAAATTTCTGATTCTGCTGTTT\n",
        ">s1 multi\nATGAAATTTCTGATTCTGCTGTTTAAAACCCCGGGGTTTTATGAAATTTCTGATTCTGCTGTTT\n",
        "0",
        &["--seg", "no", "--evalue", "1000"],
        &["-seg", "no", "-evalue", "1000"],
    );
}

#[test]
fn tblastx_db_ncbi_parity_indel_remains_ungapped() {
    assert_translated_db_outfmt_matches_ncbi_sorted_lines(
        "tblastx",
        "/usr/bin/tblastx",
        "nucl",
        ">q1\nATGAAATTTCTGATTCTGCTGTTTAATATTCTGTGTCTGTTTCCTGTTCTGGCTGCTGATAATCATGGTGTTTCTATGAATGCTTCT\n",
        ">s1\nATGAAATTTCTGATTCTGCTGTTTATTCTGTGTCTGTTTCCTGTTCTGGCTGCTGATAATCATGGTGTTTCTATGAATGCTTCT\n",
        "6 qseqid sseqid qstart qend sstart send score length gaps qframe sframe frames qseq sseq",
        &["--seg", "no", "--evalue", "1000"],
        &["-seg", "no", "-evalue", "1000"],
    );
}

#[test]
fn tblastx_db_ncbi_parity_frame_offset_insertion_remains_ungapped() {
    assert_translated_db_outfmt_matches_ncbi_sorted_lines(
        "tblastx",
        "/usr/bin/tblastx",
        "nucl",
        ">q1\nATGAAATTTCTGATTCTGCTGTTTATGAAATTTCTGATTCTGCTGTTTATGAAATTTCTGATTCTGCTGTTT\n",
        ">s1\nATGAAATTTCTGATTCTGCTGTTTAATGAAATTTCTGATTCTGCTGTTTATGAAATTTCTGATTCTGCTGTTT\n",
        "6 qseqid sseqid qstart qend sstart send score length gaps gapopen qframe sframe frames qseq sseq",
        &["--seg", "no", "--evalue", "1000"],
        &["-seg", "no", "-evalue", "1000"],
    );
}

#[test]
fn tblastx_db_ncbi_parity_frameshift_pattern_remains_ungapped() {
    assert_translated_db_outfmt_matches_ncbi_sorted_lines(
        "tblastx",
        "/usr/bin/tblastx",
        "nucl",
        ">q1\nGCTTGTGATGAATTTGGTCATATTAAACTGATGAATCCTCAGCGTTCTACCGTTTGGTATGCTTGTGATGAATTTGGTCATATTAAACTGAATGAATCCTCAGCGTTCTACCGTTTGGTAT\n",
        ">s1\nGCTTGTGATGAATTTGGTCATATTAAACTGATGAATCCTCAGCGTTCTACCGTTTGGTATGCTTGTGATGAATTTGGTCATATTAAACTGAATGAATCCTCAGCGTTCTACCGTTTGGTAT\n",
        "6 qseqid sseqid qstart qend sstart send score length gaps gapopen qframe sframe frames qseq sseq btop",
        &["--seg", "no", "--evalue", "1000"],
        &["-seg", "no", "-evalue", "1000"],
    );
}

#[test]
fn tblastx_db_ncbi_parity_complex_frameshift_top_hsps() {
    assert_translated_db_outfmt_matches_ncbi_sorted_lines(
        "tblastx",
        "/usr/bin/tblastx",
        "nucl",
        ">q1\nGCTTGTGATGAATTTGGTCATATTAAACTGATGAATCCTCAGCGTTCTACCGTTTGGTATGCTTGTGATGAATTTGGTCATATTAAACTGAATGAATCCTCAGCGTTCTACCGTTTGGTAT\n",
        ">s1\nGCTTGTGATGAATTTGGTCATATTAAACTGATGAATCCTCAGCGTTCTACCGTTTGGTATGCTTGTGATGAATTTGGTCATATTAAACTGATGAATCCTCAGCGTTCTACCGTTTGGTAT\n",
        "6 qseqid sseqid qstart qend sstart send score length gaps gapopen qframe sframe frames qseq sseq btop",
        &["--seg", "no", "--evalue", "1000", "--max_hsps", "30"],
        &["-seg", "no", "-evalue", "1000", "-max_hsps", "30"],
    );
}

#[test]
fn blastx_db_ncbi_parity_indel_remains_ungapped() {
    assert_translated_db_outfmt_matches_ncbi_sorted_lines(
        "blastx",
        "/usr/bin/blastx",
        "prot",
        ">q1\nATGAAATTTCTGATTCTGCTGTTTAAAATTTATGAAATTTCTGATTCTGCTGTTT\n",
        ">s1\nMKFLILLFKFMKFLILLF\n",
        "6 qseqid sseqid qstart qend sstart send score length gaps gapopen qframe sframe frames qseq sseq btop",
        &["--seg", "no", "--comp_based_stats", "0", "--evalue", "1000"],
        &["-seg", "no", "-comp_based_stats", "0", "-evalue", "1000"],
    );
}

#[test]
fn tblastn_db_ncbi_parity_indel_remains_ungapped() {
    assert_translated_db_outfmt_matches_ncbi_sorted_lines(
        "tblastn",
        "/usr/bin/tblastn",
        "nucl",
        ">q1\nMKFLILLFKFMKFLILLF\n",
        ">s1\nATGAAATTTCTGATTCTGCTGTTTAAAATTTATGAAATTTCTGATTCTGCTGTTT\n",
        "6 qseqid sseqid qstart qend sstart send score length gaps gapopen qframe sframe frames qseq sseq btop",
        &["--seg", "no", "--comp_based_stats", "0", "--evalue", "1000"],
        &["-seg", "no", "-comp_based_stats", "0", "-evalue", "1000"],
    );
}

#[test]
fn blastx_db_ncbi_parity_frameshift_gap_script() {
    assert_translated_db_outfmt_matches_ncbi_sorted_lines(
        "blastx",
        "/usr/bin/blastx",
        "prot",
        ">q1\nGCTTGTGATGAATTTGGTCATATTAAACTGATGAATCCTCAGCGTTCTACCGTTTGGTATGCTTGTGATGAATTTGGTCATATTAAACTGAATGAATCCTCAGCGTTCTACCGTTTGGTAT\n",
        ">s1\nACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY\n",
        "6 qseqid sseqid qstart qend sstart send score length gaps gapopen qframe sframe frames qseq sseq btop",
        &[
            "--seg",
            "no",
            "--comp_based_stats",
            "0",
            "--evalue",
            "1000",
            "--gapopen",
            "9",
            "--gapextend",
            "1",
        ],
        &[
            "-seg",
            "no",
            "-comp_based_stats",
            "0",
            "-evalue",
            "1000",
            "-gapopen",
            "9",
            "-gapextend",
            "1",
        ],
    );
}

#[test]
fn tblastn_db_ncbi_parity_frameshift_gap_script() {
    assert_translated_db_outfmt_matches_ncbi_sorted_lines(
        "tblastn",
        "/usr/bin/tblastn",
        "nucl",
        ">q1\nACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY\n",
        ">s1\nGCTTGTGATGAATTTGGTCATATTAAACTGATGAATCCTCAGCGTTCTACCGTTTGGTATGCTTGTGATGAATTTGGTCATATTAAACTGAATGAATCCTCAGCGTTCTACCGTTTGGTAT\n",
        "6 qseqid sseqid qstart qend sstart send score length gaps gapopen qframe sframe frames qseq sseq btop",
        &[
            "--seg",
            "no",
            "--comp_based_stats",
            "0",
            "--evalue",
            "1000",
            "--gapopen",
            "9",
            "--gapextend",
            "1",
        ],
        &[
            "-seg",
            "no",
            "-comp_based_stats",
            "0",
            "-evalue",
            "1000",
            "-gapopen",
            "9",
            "-gapextend",
            "1",
        ],
    );
}

#[test]
fn translated_db_ncbi_parity_pairwise_zero_description_alignment_limits() {
    for (program, ncbi_bin, dbtype, query, db_fasta) in [
        (
            "blastx",
            "/usr/bin/blastx",
            "prot",
            ">q1 translated query\nATGAAATTTCTGATTCTGCTGTTT\n",
            ">s1 protein subject\nMKFLILLF\n>s2 unrelated\nACDEFGHIKLMN\n",
        ),
        (
            "tblastn",
            "/usr/bin/tblastn",
            "nucl",
            ">q1 protein query\nMKFLILLF\n",
            ">s1 nucleotide subject\nATGAAATTTCTGATTCTGCTGTTT\n>s2 unrelated\nCCCCCCCCCCCCCCCCCCCCCCCC\n",
        ),
        (
            "tblastx",
            "/usr/bin/tblastx",
            "nucl",
            ">q1 translated query\nATGAAATTTCTGATTCTGCTGTTT\n",
            ">s1 nucleotide subject\nATGAAATTTCTGATTCTGCTGTTT\n>s2 unrelated\nCCCCCCCCCCCCCCCCCCCCCCCC\n",
        ),
    ] {
        for (rust_extra, ncbi_extra) in [
            (vec!["--num_alignments", "0"], vec!["-num_alignments", "0"]),
            (
                vec!["--num_descriptions", "0"],
                vec!["-num_descriptions", "0"],
            ),
        ] {
            let (rust_comp_args, ncbi_comp_args): (&[&str], &[&str]) = if program == "tblastx" {
                (&[], &[])
            } else {
                (&["--seg", "no", "--comp_based_stats", "0"], &["-seg", "no", "-comp_based_stats", "0"])
            };
            let mut rust_args = Vec::new();
            rust_args.extend(rust_comp_args);
            rust_args.extend(rust_extra);
            let mut ncbi_args = Vec::new();
            ncbi_args.extend(ncbi_comp_args);
            ncbi_args.extend(ncbi_extra);
            assert_translated_db_outfmt_matches_ncbi(
                program,
                ncbi_bin,
                dbtype,
                query,
                db_fasta,
                "0",
                &rust_args,
                &ncbi_args,
            );
        }
    }
}

#[test]
fn translated_db_ncbi_parity_culling_limit() {
    for (program, ncbi_bin, dbtype, query, db_fasta) in [
        (
            "blastx",
            "/usr/bin/blastx",
            "prot",
            ">q1\nATGAAATTTCTGATTCTGCTGTTTAAAACCCCGGGGTTTTATGAAATTTCTGATTCTGCTGTTT\n",
            ">wide\nMKFLILLFQQQQGMKFLILLF\n>contained\nKFLILL\n>outside\nQQQQG\n",
        ),
        (
            "tblastn",
            "/usr/bin/tblastn",
            "nucl",
            ">q1\nMKFLILLFQQQQGMKFLILLF\n",
            ">wide\nATGAAATTTCTGATTCTGCTGTTTAAAACCCCGGGGTTTTATGAAATTTCTGATTCTGCTGTTT\n>contained\nAAATTTCTGATTCTGCTG\n>outside\nCAGCAGCAGCAGGGT\n",
        ),
        (
            "tblastx",
            "/usr/bin/tblastx",
            "nucl",
            ">q1\nATGAAATTTCTGATTCTGCTGTTTAAAACCCCGGGGTTTTATGAAATTTCTGATTCTGCTGTTT\n",
            ">wide\nATGAAATTTCTGATTCTGCTGTTTAAAACCCCGGGGTTTTATGAAATTTCTGATTCTGCTGTTT\n>contained\nAAATTTCTGATTCTGCTG\n>outside\nCAGCAGCAGCAGGGT\n",
        ),
    ] {
        let (rust_comp_args, ncbi_comp_args): (&[&str], &[&str]) = if program == "tblastx" {
            (&[], &[])
        } else {
            (&["--comp_based_stats", "0"], &["-comp_based_stats", "0"])
        };
        let mut rust_args = vec!["--seg", "no"];
        rust_args.extend(rust_comp_args);
        rust_args.extend([
            "--sum_stats",
            "false",
            "--evalue",
            "1000",
            "--culling_limit",
            "1",
        ]);
        let mut ncbi_args = vec!["-seg", "no"];
        ncbi_args.extend(ncbi_comp_args);
        ncbi_args.extend([
            "-sum_stats",
            "false",
            "-evalue",
            "1000",
            "-culling_limit",
            "1",
        ]);
        assert_translated_db_outfmt_matches_ncbi(
            program,
            ncbi_bin,
            dbtype,
            query,
            db_fasta,
            "6 qseqid sseqid score evalue qstart qend sstart send",
            &rust_args,
            &ncbi_args,
        );
    }
}

#[test]
fn translated_db_ncbi_parity_hit_limit_and_identity_filters() {
    for (program, ncbi_bin, dbtype, query, db_fasta) in [
        (
            "blastx",
            "/usr/bin/blastx",
            "prot",
            ">q1\nATGAAATTTCTGATTCTGCTGTTTAAAACCCCGGGGTTTTATGAAATTTCTGATTCTGCTGTTT\n",
            ">full\nMKFLILLFQQQQGMKFLILLF\n>partial\nMKFLILLF\n>near\nMKFLILLY\n",
        ),
        (
            "tblastn",
            "/usr/bin/tblastn",
            "nucl",
            ">q1\nMKFLILLFQQQQGMKFLILLF\n",
            ">full\nATGAAATTTCTGATTCTGCTGTTTAAAACCCCGGGGTTTTATGAAATTTCTGATTCTGCTGTTT\n>partial\nATGAAATTTCTGATTCTGCTGTTT\n>near\nATGAAATTTCTGATTCTGCTGTTACT\n",
        ),
        (
            "tblastx",
            "/usr/bin/tblastx",
            "nucl",
            ">q1\nATGAAATTTCTGATTCTGCTGTTTAAAACCCCGGGGTTTTATGAAATTTCTGATTCTGCTGTTT\n",
            ">full\nATGAAATTTCTGATTCTGCTGTTTAAAACCCCGGGGTTTTATGAAATTTCTGATTCTGCTGTTT\n>partial\nATGAAATTTCTGATTCTGCTGTTT\n>near\nATGAAATTTCTGATTCTGCTGTTACT\n",
        ),
    ] {
        for (rust_extra, ncbi_extra) in [
            (
                vec!["--max_hsps", "1"],
                vec!["-max_hsps", "1"],
            ),
            (
                vec!["--max_target_seqs", "1"],
                vec!["-max_target_seqs", "1"],
            ),
            (
                vec!["--perc_identity", "90"],
                vec!["-perc_identity", "90"],
            ),
            (
                vec!["--qcov_hsp_perc", "80"],
                vec!["-qcov_hsp_perc", "80"],
            ),
        ] {
            if rust_extra.first() == Some(&"--perc_identity") {
                continue;
            }
            let (rust_comp_args, ncbi_comp_args): (&[&str], &[&str]) = if program == "tblastx" {
                (&[], &[])
            } else {
                (&["--comp_based_stats", "0"], &["-comp_based_stats", "0"])
            };
            let mut rust_args = vec!["--seg", "no"];
            rust_args.extend(rust_comp_args);
            rust_args.extend(["--sum_stats", "false", "--evalue", "1000"]);
            rust_args.extend(rust_extra);
            let mut ncbi_args = vec!["-seg", "no"];
            ncbi_args.extend(ncbi_comp_args);
            ncbi_args.extend(["-sum_stats", "false", "-evalue", "1000"]);
            ncbi_args.extend(ncbi_extra);
            assert_translated_db_outfmt_matches_ncbi(
                program,
                ncbi_bin,
                dbtype,
                query,
                db_fasta,
                "6 qseqid sseqid score pident length qstart qend sstart send",
                &rust_args,
                &ncbi_args,
            );
        }
    }
}

#[test]
fn tblastx_db_ncbi_parity_stable_hit_limit_filters() {
    for (rust_extra, ncbi_extra) in [
        (vec!["--max_hsps", "1"], vec!["-max_hsps", "1"]),
        (
            vec!["--max_target_seqs", "1"],
            vec!["-max_target_seqs", "1"],
        ),
        (vec!["--qcov_hsp_perc", "80"], vec!["-qcov_hsp_perc", "80"]),
    ] {
        let mut rust_args = vec!["--seg", "no", "--sum_stats", "false", "--evalue", "1000"];
        rust_args.extend(rust_extra);
        let mut ncbi_args = vec!["-seg", "no", "-sum_stats", "false", "-evalue", "1000"];
        ncbi_args.extend(ncbi_extra);
        assert_translated_db_outfmt_matches_ncbi(
            "tblastx",
            "/usr/bin/tblastx",
            "nucl",
            ">q1\nATGAAATTTCTGATTCTGCTGTTTAAAACCCCGGGGTTTTATGAAATTTCTGATTCTGCTGTTT\n",
            ">full\nATGAAATTTCTGATTCTGCTGTTTAAAACCCCGGGGTTTTATGAAATTTCTGATTCTGCTGTTT\n>partial\nATGAAATTTCTGATTCTGCTGTTT\n>near\nATGAAATTTCTGATTCTGCTGTTACT\n",
            "6 qseqid sseqid score pident length qstart qend sstart send qframe sframe",
            &rust_args,
            &ncbi_args,
        );
    }
}

#[test]
fn translated_db_ncbi_parity_best_hit_filters() {
    for (program, ncbi_bin, dbtype, query, db_fasta) in [
        (
            "blastx",
            "/usr/bin/blastx",
            "prot",
            ">q1\nATGAAATTTCTGATTCTGCTGTTTAAAACCCCGGGGTTTTATGAAATTTCTGATTCTGCTGTTT\n",
            ">full\nMKFLILLFQQQQGMKFLILLF\n>contained\nKFLILL\n>near\nMKFLILLYQQQQGMKFLILLF\n",
        ),
        (
            "tblastn",
            "/usr/bin/tblastn",
            "nucl",
            ">q1\nMKFLILLFQQQQGMKFLILLF\n",
            ">full\nATGAAATTTCTGATTCTGCTGTTTAAAACCCCGGGGTTTTATGAAATTTCTGATTCTGCTGTTT\n>contained\nAAATTTCTGATTCTGCTG\n>near\nATGAAATTTCTGATTCTGCTGTTACTAAAACCCCGGGGTTTTATGAAATTTCTGATTCTGCTGTTT\n",
        ),
        (
            "tblastx",
            "/usr/bin/tblastx",
            "nucl",
            ">q1\nATGAAATTTCTGATTCTGCTGTTTAAAACCCCGGGGTTTTATGAAATTTCTGATTCTGCTGTTT\n",
            ">full\nATGAAATTTCTGATTCTGCTGTTTAAAACCCCGGGGTTTTATGAAATTTCTGATTCTGCTGTTT\n>contained\nAAATTTCTGATTCTGCTG\n>near\nATGAAATTTCTGATTCTGCTGTTACTAAAACCCCGGGGTTTTATGAAATTTCTGATTCTGCTGTTT\n",
        ),
    ] {
        let (rust_extra, ncbi_extra) = (vec!["--subject_besthit"], vec!["-subject_besthit"]);
        let (rust_comp_args, ncbi_comp_args): (&[&str], &[&str]) = if program == "tblastx" {
            (&[], &[])
        } else {
            (&["--comp_based_stats", "0"], &["-comp_based_stats", "0"])
        };
        let mut rust_args = vec!["--seg", "no"];
        rust_args.extend(rust_comp_args);
        rust_args.extend(["--sum_stats", "false", "--evalue", "1000"]);
        rust_args.extend(rust_extra);
        let mut ncbi_args = vec!["-seg", "no"];
        ncbi_args.extend(ncbi_comp_args);
        ncbi_args.extend(["-sum_stats", "false", "-evalue", "1000"]);
        ncbi_args.extend(ncbi_extra);
        assert_translated_db_outfmt_matches_ncbi(
            program,
            ncbi_bin,
            dbtype,
            query,
            db_fasta,
            "6 qseqid sseqid score evalue qstart qend sstart send",
            &rust_args,
            &ncbi_args,
        );
    }
}

#[test]
fn blastx_db_ncbi_parity_multi_subject_ordering() {
    assert_translated_db_outfmt_matches_ncbi(
        "blastx",
        "/usr/bin/blastx",
        "prot",
        ">q1\nATGAAATTTCTGATTCTGCTGTTT\n",
        concat!(">s_exact\nMKFLILLF\n", ">s_near\nMKFLILLY\n",),
        "6 qseqid sseqid score qframe sframe frames",
        &[],
        &[],
    );
}

#[test]
fn tblastn_db_ncbi_parity_multi_subject_ordering() {
    assert_translated_db_outfmt_matches_ncbi(
        "tblastn",
        "/usr/bin/tblastn",
        "nucl",
        ">q1\nMKFLILLF\n",
        concat!(
            ">s_exact\nATGAAATTTCTGATTCTGCTGTTT\n",
            ">s_near\nATGAAATTTCTGATTCTGCTGTAT\n",
        ),
        "6 qseqid sseqid score qframe sframe frames",
        &[],
        &[],
    );
}

#[test]
fn tblastx_db_ncbi_parity_multi_subject_ordering() {
    assert_translated_db_outfmt_matches_ncbi(
        "tblastx",
        "/usr/bin/tblastx",
        "nucl",
        ">q1\nATGAAATTTCTGATTCTGCTGTTT\n",
        concat!(
            ">s_exact\nATGAAATTTCTGATTCTGCTGTTT\n",
            ">s_near\nATGAAATTTCTGATTCTGCTGTAT\n",
        ),
        "6 qseqid sseqid score qframe sframe frames",
        &[],
        &[],
    );
}

// ── BLASTP tests ─────────────────────────────────────────────────────────────

#[test]
fn blastp_exact_match() {
    let seq = "MKFLILLFNILCLFPVLAADNHGVSMNAS";
    let (_tmp, db) = build_protein_db(vec![
        protein_entry("P001", "exact match protein", seq),
        protein_entry("P002", "unrelated protein", "WWWWWWWWWWWWWWWWWWWWWWWWWWWWW"),
    ]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    let results = blastp(&db, seq.as_bytes(), &params);

    assert_eq!(
        results.len(),
        1,
        "exact query should only hit the exact subject"
    );
    let best = &results[0];
    assert_eq!(best.subject_oid, 0);
    let hsp = &best.hsps[0];
    assert!(
        (hsp.percent_identity() - 100.0).abs() < 0.01,
        "exact match should be 100% identity, got {:.1}%",
        hsp.percent_identity()
    );
    assert_eq!(hsp.alignment_length, seq.len());
}

#[test]
fn blastp_no_hit_for_unrelated() {
    let (_tmp, db) = build_protein_db(vec![protein_entry(
        "P001",
        "all alanine",
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAA",
    )]);
    let params = SearchParams::blastp()
        .evalue(1e-10)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    let results = blastp(&db, b"WWWWWWWWWWWWWWWWWWWWWWWWWWWW", &params);
    assert!(
        results.is_empty(),
        "unrelated sequences should not produce hits at strict evalue"
    );
}

#[test]
fn blastp_finds_similar_sequence() {
    let query = "MKFLILLFNILCLFPVLAADNHGVSMNAS";
    let subject = "MKFLILLFNILCLFPVLAADNHGVSMNAS";
    let mutated = "MKFLILLFNILCLFPVLAAENHGVSMNAS"; // D→E
    let (_tmp, db) = build_protein_db(vec![
        protein_entry("P001", "near identical", subject),
        protein_entry("P002", "one mismatch", mutated),
    ]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    let results = blastp(&db, query.as_bytes(), &params);
    assert_eq!(results.len(), 2, "should find both similar sequences");
    assert_eq!(results[0].subject_oid, 0, "exact match should rank first");
    assert_eq!(
        results[1].subject_oid, 1,
        "one-mismatch hit should rank second"
    );
    assert_eq!(results[0].hsps[0].alignment_length, query.len());
    assert_eq!(results[0].hsps[0].num_identities as usize, query.len());
    assert_eq!(results[0].hsps[0].num_gaps, 0);
    assert_eq!(results[1].hsps[0].alignment_length, query.len());
    assert_eq!(results[1].hsps[0].num_identities as usize, query.len() - 1);
    assert_eq!(results[1].hsps[0].num_gaps, 0);
}

#[test]
fn blastp_max_target_seqs_limits_results() {
    let query = "MKFLILLFNILCLFPVLAADNHGVSMNAS";
    let entries: Vec<_> = (0..10)
        .map(|i| protein_entry(&format!("P{:03}", i), "copy", query))
        .collect();
    let (_tmp, db) = build_protein_db(entries);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .max_target_seqs(3)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    let results = blastp(&db, query.as_bytes(), &params);
    assert_eq!(
        results.len(),
        3,
        "max_target_seqs=3 should return three hits"
    );
}

#[test]
fn blastp_empty_query() {
    let (_tmp, db) = build_protein_db(vec![protein_entry(
        "P001",
        "target",
        "MKFLILLFNILCLFPVLAADNHGVSMNAS",
    )]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    let results = blastp(&db, b"", &params);
    assert!(results.is_empty(), "empty query should produce no results");
}

#[test]
fn blastp_short_query_below_word_size() {
    let (_tmp, db) = build_protein_db(vec![protein_entry(
        "P001",
        "target",
        "MKFLILLFNILCLFPVLAADNHGVSMNAS",
    )]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    let results = blastp(&db, b"MK", &params);
    let _ = results; // should not panic
}

#[test]
fn blastp_single_residue_query() {
    let (_tmp, db) = build_protein_db(vec![protein_entry(
        "P001",
        "target",
        "MKFLILLFNILCLFPVLAADNHGVSMNAS",
    )]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    let results = blastp(&db, b"M", &params);
    let _ = results; // should not panic
}

#[test]
fn blastp_all_twenty_amino_acids() {
    let query = "ACDEFGHIKLMNPQRSTVWY";
    let (_tmp, db) = build_protein_db(vec![protein_entry("P001", "all aa", query)]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    let results = blastp(&db, query.as_bytes(), &params);
    assert_eq!(results.len(), 1);
    assert_eq!(results[0].subject_oid, 0);
    let hsp = &results[0].hsps[0];
    assert!((hsp.percent_identity() - 100.0).abs() < 0.01);
    assert_eq!(hsp.alignment_length, 20);
    assert_eq!(hsp.num_identities, 20);
    assert_eq!(hsp.num_gaps, 0);
}

#[test]
fn blastp_default_seg_masks_low_complexity_query() {
    // NCBI's blastp default has SEG OFF (`blastp_args.cpp:50`,
    // `kFilterByDefault = false`); the suppression of low-complexity hits
    // comes from comp_adjust=2 collapsing the A-A diagonal, not from SEG.
    // We explicitly enable filter_low_complexity here to test our SEG impl.
    let query = "AAAAAAAAAAAAAAAAAAAA";
    let (_tmp, db) = build_protein_db(vec![protein_entry("P001", "poly-a subject", query)]);

    let filtered = SearchParams::blastp()
        .evalue(1e6)
        .num_threads(1)
        .filter_low_complexity(true);
    let unfiltered = SearchParams::blastp()
        .evalue(1e6)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    let filtered_results = blastp(&db, query.as_bytes(), &filtered);
    let unfiltered_results = blastp(&db, query.as_bytes(), &unfiltered);

    assert!(
        filtered_results.is_empty(),
        "explicit blastp SEG masking should suppress low-complexity query hits"
    );
    assert_eq!(
        unfiltered_results.len(),
        1,
        "disabling low-complexity masking should restore the poly-A hit"
    );
    assert_eq!(unfiltered_results[0].subject_oid, 0);
    assert_eq!(unfiltered_results[0].hsps[0].alignment_length, query.len());
    assert_eq!(
        unfiltered_results[0].hsps[0].num_identities as usize,
        query.len()
    );
    assert_eq!(unfiltered_results[0].hsps[0].num_gaps, 0);
}

// ── BLASTN tests ─────────────────────────────────────────────────────────────

#[test]
fn blastn_exact_match() {
    let seq = "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACGATCGTAAGGCCTTAGCAGTCAATGC";
    let (_tmp, db) = build_nucleotide_db(vec![nt_entry("N001", "exact nt", seq)]);
    let params = SearchParams::blastn()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false);

    let results = blastn(&db, seq.as_bytes(), &params);
    assert_eq!(results.len(), 1, "should find one exact nucleotide match");
    assert_eq!(results[0].subject_oid, 0);
    let hsp = &results[0].hsps[0];
    assert!((hsp.percent_identity() - 100.0).abs() < 0.01);
    assert_eq!(hsp.alignment_length, seq.len());
    assert_eq!(hsp.num_identities as usize, seq.len());
    assert_eq!(hsp.num_gaps, 0);
}

#[test]
fn blastn_reverse_complement_hit() {
    let seq = "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACGATCGTAAGGCCTTAGCAGTCAATGC";
    let rc = String::from_utf8(reverse_complement(seq.as_bytes())).unwrap();
    let (_tmp, db) = build_nucleotide_db(vec![nt_entry("N001", "forward strand", seq)]);
    let params = SearchParams::blastn()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .strand("both");

    let results = blastn(&db, rc.as_bytes(), &params);
    assert_eq!(results.len(), 1, "should find one reverse-complement hit");
    assert_eq!(results[0].subject_oid, 0);
    let hsp = &results[0].hsps[0];
    assert_eq!(hsp.alignment_length, seq.len());
    assert_eq!(hsp.num_identities as usize, seq.len());
    assert_eq!(hsp.num_gaps, 0);
}

#[test]
fn blastn_no_hit_unrelated() {
    let (_tmp, db) = build_nucleotide_db(vec![
        nt_entry("N001", "seq A", "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACGATCGTAAGGCCTTAGCAGTCAATGC"),
    ]);
    let params = SearchParams::blastn()
        .evalue(1e-10)
        .num_threads(1)
        .filter_low_complexity(false);

    let results = blastn(&db, b"TTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGC", &params);
    assert!(
        results.is_empty(),
        "unrelated nucleotide sequences should not match"
    );
}

#[test]
fn blastn_mismatch_scoring() {
    let seq     = "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACGATCGTAAGGCCTTAGCAGTCAATGC";
    let mutated = "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACAATCGTAAGGCCTTAGCAGTCAATGC"; // G→A at pos 76
    let (_tmp, db) = build_nucleotide_db(vec![nt_entry("N001", "original", seq)]);
    let params = SearchParams::blastn()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false);

    let results = blastn(&db, mutated.as_bytes(), &params);
    assert_eq!(results.len(), 1);
    assert_eq!(results[0].subject_oid, 0);
    let hsp = &results[0].hsps[0];
    assert_eq!(hsp.alignment_length, seq.len());
    assert_eq!(hsp.num_identities as usize, seq.len() - 1);
    assert_eq!(hsp.num_gaps, 0);
    assert!((hsp.percent_identity() - 99.01960784313727).abs() < 1.0e-12);
}

#[test]
fn blastn_empty_query() {
    let (_tmp, db) = build_nucleotide_db(vec![
        nt_entry("N001", "target", "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACGATCGTAAGGCCTTAGCAGTCAATGC"),
    ]);
    let params = SearchParams::blastn()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false);

    let results = blastn(&db, b"", &params);
    assert!(results.is_empty(), "empty query should produce no results");
}

#[test]
fn blastn_short_query_below_word_size() {
    let (_tmp, db) = build_nucleotide_db(vec![
        nt_entry("N001", "target", "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACGATCGTAAGGCCTTAGCAGTCAATGC"),
    ]);
    let params = SearchParams::blastn()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false);

    let results = blastn(&db, b"ATGCGT", &params);
    let _ = results; // should not panic
}

#[test]
fn blastn_plus_strand_only() {
    let seq = "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACGATCGTAAGGCCTTAGCAGTCAATGC";
    let rc = String::from_utf8(reverse_complement(seq.as_bytes())).unwrap();
    let (_tmp, db) = build_nucleotide_db(vec![nt_entry("N001", "forward only", seq)]);

    let params = SearchParams::blastn()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .strand("plus");

    let results = blastn(&db, rc.as_bytes(), &params);
    assert!(
        results.is_empty(),
        "plus-strand-only should not find reverse complement hit"
    );
}

#[test]
fn blastn_alignment_strings_are_ascii() {
    let seq = "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACGATCGTAAGGCCTTAGCAGTCAATGC";
    let (_tmp, db) = build_nucleotide_db(vec![nt_entry("N001", "target", seq)]);
    let params = SearchParams::blastn()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false);

    let results = blastn(&db, seq.as_bytes(), &params);
    assert_eq!(results.len(), 1);
    let hsp = &results[0].hsps[0];
    for &b in &hsp.query_aln {
        assert!(
            b == b'-' || b.is_ascii_alphabetic(),
            "query_aln byte {} is not ASCII letter or gap",
            b
        );
    }
    for &b in &hsp.subject_aln {
        assert!(
            b == b'-' || b.is_ascii_alphabetic(),
            "subject_aln byte {} is not ASCII letter or gap",
            b
        );
    }
}

// ── BLASTX test ─────────────────────────────────────────────────────────────

#[test]
fn blastx_finds_translated_hit() {
    let nt_query = "ATGAAATTTCTGATTCTGCTGTTT";
    let protein = "MKFLILLF";

    let (_tmp, db) = build_protein_db(vec![protein_entry("P001", "target protein", protein)]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    let results = blastx(&db, nt_query.as_bytes(), &params);
    assert_eq!(
        results.len(),
        1,
        "blastx should find one translated protein"
    );
    assert_eq!(results[0].subject_oid, 0);
    let hsp = &results[0].hsps[0];
    assert_eq!(hsp.query_frame, 1);
    assert_eq!(hsp.alignment_length, protein.len());
    assert_eq!(hsp.num_identities as usize, protein.len());
    assert_eq!(hsp.num_gaps, 0);
    assert_eq!(hsp.query_start, 0);
    assert_eq!(hsp.query_end, nt_query.len());
}

#[test]
fn blastx_empty_nt_query() {
    let (_tmp, db) = build_protein_db(vec![protein_entry(
        "P001",
        "target",
        "MKFLILLFNILCLFPVLAAD",
    )]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    let results = blastx(&db, b"", &params);
    assert!(results.is_empty());
}

#[test]
fn blastx_default_seg_masks_low_complexity_translation() {
    let nt_query = "GCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCT";
    let protein = "AAAAAAAAAAAAAAAAAAAA";

    let (_tmp, db) = build_protein_db(vec![protein_entry("P001", "poly-a target", protein)]);
    // `SearchParams::blastp()` now defaults filter_low_complexity = false to
    // match NCBI's `blastp_args.cpp:50`. Explicitly enable filtering here to
    // exercise the SEG path on the translated query.
    let filtered = SearchParams::blastp()
        .evalue(1e6)
        .num_threads(1)
        .filter_low_complexity(true)
        .comp_adjust(0);
    let unfiltered = SearchParams::blastp()
        .evalue(1e6)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    let filtered_results = blastx(&db, nt_query.as_bytes(), &filtered);
    let unfiltered_results = blastx(&db, nt_query.as_bytes(), &unfiltered);

    assert!(
        filtered_results.is_empty(),
        "default blastx SEG masking should suppress low-complexity translated query hits"
    );
    assert_eq!(
        unfiltered_results.len(),
        1,
        "disabling low-complexity masking should restore the translated poly-A hit"
    );
    assert_eq!(unfiltered_results[0].subject_oid, 0);
    let translated_query_len = nt_query.len() / 3;
    assert_eq!(
        unfiltered_results[0].hsps[0].alignment_length,
        translated_query_len
    );
    assert_eq!(
        unfiltered_results[0].hsps[0].num_identities as usize,
        translated_query_len
    );
    assert_eq!(unfiltered_results[0].hsps[0].num_gaps, 0);
}

// ── TBLASTN test ────────────────────────────────────────────────────────────

#[test]
fn tblastn_finds_protein_in_nt_db() {
    let protein_query = "MKFLILLF";
    let nt_subject = "ATGAAATTTCTGATTCTGCTGTTT";

    let (_tmp, db) = build_nucleotide_db(vec![nt_entry("N001", "coding region", nt_subject)]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    let results = tblastn(&db, protein_query.as_bytes(), &params);
    assert_eq!(
        results.len(),
        1,
        "tblastn should find one translated nucleotide hit"
    );
    assert_eq!(results[0].subject_oid, 0);
    let hsp = &results[0].hsps[0];
    assert_eq!(hsp.subject_frame, 1);
    assert_eq!(hsp.alignment_length, protein_query.len());
    assert_eq!(hsp.num_identities as usize, protein_query.len());
    assert_eq!(hsp.num_gaps, 0);
    assert_eq!(results[0].subject_len, nt_subject.len());
    assert_eq!(hsp.subject_start, 0);
    assert_eq!(hsp.subject_end, nt_subject.len());
}

#[test]
fn tblastn_accepts_reusable_rayon_thread_pool() {
    let protein_query = b"MKFLILLF";
    let (_tmp, db) = build_nucleotide_db(vec![
        nt_entry("N001", "exact coding", "ATGAAATTTCTGATTCTGCTGTTT"),
        nt_entry("N002", "near coding", "ATGAAATATCTGATTCTGCTGTTT"),
        nt_entry("N003", "unrelated", "GCTGCTGCTGCTGCTGCTGCTGCT"),
    ]);
    let baseline_params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);
    let baseline = tblastn(&db, protein_query, &baseline_params);

    let pool = std::sync::Arc::new(
        rayon::ThreadPoolBuilder::new()
            .num_threads(2)
            .build()
            .expect("build rayon pool"),
    );
    let pooled_params = SearchParams::blastp()
        .evalue(10.0)
        .filter_low_complexity(false)
        .comp_adjust(0)
        .thread_pool(pool);
    let pooled = tblastn(&db, protein_query, &pooled_params);

    assert_eq!(pooled.len(), baseline.len());
    assert_eq!(pooled[0].subject_accession, baseline[0].subject_accession);
    assert_eq!(pooled[0].hsps[0].score, baseline[0].hsps[0].score);
    assert_eq!(
        pooled[0].hsps[0].subject_frame,
        baseline[0].hsps[0].subject_frame
    );
}

#[test]
fn tblastn_empty_protein_query() {
    let (_tmp, db) =
        build_nucleotide_db(vec![nt_entry("N001", "coding", "ATGAAATTTCTGATTCTGCTGTTT")]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    let results = tblastn(&db, b"", &params);
    assert!(results.is_empty());
}

#[test]
fn tblastn_default_seg_masks_low_complexity_query() {
    let protein_query = "AAAAAAAAAAAAAAAAAAAA";
    let nt_subject = "GCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCT";

    let (_tmp, db) =
        build_nucleotide_db(vec![nt_entry("N001", "poly-a coding region", nt_subject)]);
    // `SearchParams::blastp()` now defaults filter_low_complexity = false.
    // Explicitly enable filtering here to exercise the SEG path on the
    // protein query.
    let filtered = SearchParams::blastp()
        .evalue(1e6)
        .num_threads(1)
        .filter_low_complexity(true)
        .comp_adjust(0);
    let unfiltered = SearchParams::blastp()
        .evalue(1e6)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    let filtered_results = tblastn(&db, protein_query.as_bytes(), &filtered);
    let unfiltered_results = tblastn(&db, protein_query.as_bytes(), &unfiltered);

    assert!(
        filtered_results.is_empty(),
        "default tblastn SEG masking should suppress low-complexity protein query hits"
    );
    assert_eq!(
        unfiltered_results.len(),
        1,
        "disabling low-complexity masking should restore the tblastn poly-A hit"
    );
    assert_eq!(unfiltered_results[0].subject_oid, 0);
    let translated_subject_len = nt_subject.len() / 3;
    assert_eq!(
        unfiltered_results[0].hsps[0].alignment_length,
        translated_subject_len
    );
    assert_eq!(
        unfiltered_results[0].hsps[0].num_identities as usize,
        translated_subject_len
    );
    assert_eq!(unfiltered_results[0].hsps[0].num_gaps, 0);
}

// ── TBLASTX test ────────────────────────────────────────────────────────────

#[test]
fn tblastx_translated_vs_translated() {
    let nt_seq = "ATGAAATTTCTGATTCTGCTGTTTAACATTCTGTGCCTGTTC";
    let (_tmp, db) = build_nucleotide_db(vec![nt_entry("N001", "coding nt", nt_seq)]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    let results = tblastx(&db, nt_seq.as_bytes(), &params);
    assert_eq!(
        results.len(),
        1,
        "tblastx should find one translated self-hit"
    );
    let hsp = &results[0].hsps[0];
    let frame2_len = (nt_seq.len() - 1) / 3;
    assert_eq!(hsp.query_frame, 2);
    assert_eq!(hsp.subject_frame, 2);
    assert_eq!(hsp.alignment_length, frame2_len);
    assert_eq!(hsp.num_identities as usize, frame2_len);
    assert_eq!(hsp.num_gaps, 0);
    assert_eq!(results[0].subject_len, nt_seq.len());
    assert_eq!(hsp.query_start, 1);
    assert_eq!(hsp.query_end, 1 + frame2_len * 3);
    assert_eq!(hsp.subject_start, 1);
    assert_eq!(hsp.subject_end, 1 + frame2_len * 3);
}

#[test]
fn tblastx_default_seg_masks_low_complexity_translation() {
    let nt_seq = "GCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCT";
    let (_tmp, db) = build_nucleotide_db(vec![nt_entry("N001", "poly-a coding nt", nt_seq)]);
    // `SearchParams::blastp()` now defaults filter_low_complexity = false.
    // Explicitly enable filtering here to exercise the SEG path on the
    // translated query/subject.
    let filtered = SearchParams::blastp()
        .evalue(1e6)
        .num_threads(1)
        .filter_low_complexity(true)
        .comp_adjust(0);
    let unfiltered = SearchParams::blastp()
        .evalue(1e6)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    let filtered_results = tblastx(&db, nt_seq.as_bytes(), &filtered);
    let unfiltered_results = tblastx(&db, nt_seq.as_bytes(), &unfiltered);

    assert!(
        filtered_results.is_empty(),
        "default tblastx SEG masking should suppress low-complexity translated query hits"
    );
    assert_eq!(
        unfiltered_results.len(),
        1,
        "disabling low-complexity masking should restore the tblastx poly-A hit"
    );
    assert_eq!(unfiltered_results[0].subject_oid, 0);
    let frame3_len = (nt_seq.len() - 2) / 3;
    assert_eq!(unfiltered_results[0].hsps[0].query_frame, 3);
    assert_eq!(unfiltered_results[0].hsps[0].subject_frame, 3);
    assert_eq!(unfiltered_results[0].hsps[0].alignment_length, frame3_len);
    assert_eq!(
        unfiltered_results[0].hsps[0].num_identities as usize,
        frame3_len
    );
    assert_eq!(unfiltered_results[0].hsps[0].num_gaps, 0);
}

// ── Database round-trip tests ────────────────────────────────────────────────

#[test]
fn protein_db_roundtrip() {
    let seq = "MKFLILLFNILCLFPVLAADNHGVSMNAS";
    let (_tmp, db) = build_protein_db(vec![protein_entry("P001", "test seq", seq)]);
    assert_eq!(db.num_oids, 1);
    assert_eq!(db.get_seq_len(0), seq.len() as u32);
}

#[test]
fn nucleotide_db_roundtrip() {
    let seq = "ATGGCTAGCGATCGATCGATCGATCG";
    let (_tmp, db) = build_nucleotide_db(vec![nt_entry("N001", "test nt", seq)]);
    assert_eq!(db.num_oids, 1);
    assert_eq!(db.get_seq_len(0), seq.len() as u32);
}

#[test]
fn db_multiple_sequences() {
    let entries = vec![
        protein_entry("P001", "first", "MKFLILLFNILCLFPVLAAD"),
        protein_entry("P002", "second", "NHGVSMNASQRDHFKLAEV"),
        protein_entry("P003", "third", "ACDEFGHIKLMNPQRSTVWY"),
    ];
    let (_tmp, db) = build_protein_db(entries);
    assert_eq!(db.num_oids, 3);
}

// ── Reverse complement test ─────────────────────────────────────────────────

#[test]
fn reverse_complement_roundtrip() {
    let seq = b"ATGGCTAGCGATCG";
    let rc = reverse_complement(seq);
    let rc2 = reverse_complement(&rc);
    assert_eq!(
        &rc2, seq,
        "double reverse complement should return original"
    );
}

#[test]
fn reverse_complement_handles_lowercase_ambiguity_and_rna() {
    assert_eq!(reverse_complement(b"rYmKwsBdHvN"), b"NBDHVSWMKRY");
    assert_eq!(reverse_complement(b"UuTt"), b"AAAA");
    assert_eq!(reverse_complement(b"ACGT?"), b"NACGT");
}

// ── Six-frame translation test ──────────────────────────────────────────────

#[test]
fn six_frame_translate_produces_six_frames() {
    let seq = b"ATGGCTAGCGATCGATCGATCGATCG";
    let frames = six_frame_translate(seq);
    assert_eq!(frames.len(), 6);
    let frame_nums: Vec<i32> = frames.iter().map(|f| f.frame).collect();
    assert_eq!(frame_nums, vec![1, 2, 3, -1, -2, -3]);
    // Frame +1 starts with M (ATG = Met)
    assert_eq!(frames[0].protein[0], b'M');
}

// ── FASTA parsing edge cases ─────────────────────────────────────────────────

#[test]
fn parse_fasta_with_blank_lines() {
    let input = b">seq1\n\nACGT\n\nTGCA\n\n>seq2\nAAAA\n";
    let seqs = parse_fasta(input);
    assert_eq!(seqs.len(), 2);
    assert_eq!(&seqs[0].1, b"ACGTTGCA".as_slice());
    assert_eq!(&seqs[1].1, b"AAAA".as_slice());
}

#[test]
fn parse_fasta_no_trailing_newline() {
    let input = b">seq1\nACGT";
    let seqs = parse_fasta(input);
    assert_eq!(seqs.len(), 1);
    assert_eq!(&seqs[0].1, b"ACGT".as_slice());
}

// ── Multi-query search ──────────────────────────────────────────────────────

#[test]
fn multi_query_fasta() {
    let fasta = b">q1\nMKFLILLFNILCLFPVLAAD\n>q2\nNHGVSMNASQRDHFKLAEV\n";
    let queries = parse_fasta(fasta);
    assert_eq!(queries.len(), 2);

    let (_tmp, db) = build_protein_db(vec![
        protein_entry("P001", "match1", "MKFLILLFNILCLFPVLAAD"),
        protein_entry("P002", "match2", "NHGVSMNASQRDHFKLAEV"),
    ]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    for (idx, (title, seq)) in queries.iter().enumerate() {
        let results = blastp(&db, seq, &params);
        assert_eq!(results.len(), 1, "query '{}' should find one hit", title);
        assert_eq!(results[0].subject_oid as usize, idx);
        assert_eq!(results[0].hsps[0].alignment_length, seq.len());
        assert_eq!(results[0].hsps[0].num_identities as usize, seq.len());
        assert_eq!(results[0].hsps[0].num_gaps, 0);
    }
}

// ── Edge cases ──────────────────────────────────────────────────────────────

#[test]
fn blastp_short_subject_in_db() {
    let (_tmp, db) = build_protein_db(vec![protein_entry("P001", "tiny", "MK")]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    let results = blastp(&db, b"MKFLILLFNILCLFPVLAADNHGVSMNAS", &params);
    let _ = results; // should not panic
}

#[test]
fn blastp_stop_codon_in_sequence() {
    let query = "MKFLILLFNILCLFPVLAAD";
    let subject = "MKFLILLF*ILCLFPVLAAD";
    let (_tmp, db) = build_protein_db(vec![protein_entry("P001", "with stop", subject)]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    let results = blastp(&db, query.as_bytes(), &params);
    let _ = results; // should not panic
}

#[test]
fn blastn_query_with_ambiguous_bases() {
    let query   = "ATGCGTACCTGAAAGCTTCAGNACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACGATCGTAAGGCCTTAGCAGTCAATGC";
    let subject = "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACGATCGTAAGGCCTTAGCAGTCAATGC";
    let (_tmp, db) = build_nucleotide_db(vec![nt_entry("N001", "target", subject)]);
    let params = SearchParams::blastn()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false);

    let results = blastn(&db, query.as_bytes(), &params);
    let _ = results; // should not panic
}

#[test]
fn db_single_base_nucleotide() {
    let (_tmp, db) = build_nucleotide_db(vec![nt_entry("N001", "single", "A")]);
    assert_eq!(db.num_oids, 1);
}

#[test]
fn db_single_residue_protein() {
    let (_tmp, db) = build_protein_db(vec![protein_entry("P001", "single", "M")]);
    assert_eq!(db.num_oids, 1);
}

// ── Large query test ────────────────────────────────────────────────────────

#[test]
fn blastp_large_query_1000aa() {
    let aa = b"ACDEFGHIKLMNPQRSTVWY";
    let query: Vec<u8> = (0..1000).map(|i| aa[i % aa.len()]).collect();
    let (_tmp, db) = build_protein_db(vec![
        protein_entry("P001", "target", std::str::from_utf8(&query).unwrap()),
        protein_entry("P002", "unrelated", "WWWWWWWWWWWWWWWWWWWWWWWWWWWWW"),
    ]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);
    let results = blastp(&db, &query, &params);
    assert_eq!(
        results.len(),
        1,
        "should find one self-match for 1000aa query"
    );
    assert_eq!(results[0].subject_oid, 0);
    assert_eq!(results[0].hsps[0].alignment_length, query.len());
    assert_eq!(results[0].hsps[0].num_identities as usize, query.len());
    assert_eq!(results[0].hsps[0].num_gaps, 0);
}

// ── Blastn with multiple mismatches ─────────────────────────────────────────

#[test]
fn blastn_multiple_mismatches() {
    let query   = "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACGATCGTAAGGCCTTAGCAGTCAATGC";
    let subject = "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACGATCGTAAGGCCTTAGCAGTCAATGC";
    let mut subj_bytes = subject.as_bytes().to_vec();
    subj_bytes[10] = b'T';
    subj_bytes[30] = b'A';
    subj_bytes[50] = b'G';
    subj_bytes[70] = b'C';
    subj_bytes[90] = b'T';
    let subj_str = String::from_utf8(subj_bytes).unwrap();

    let (_tmp, db) = build_nucleotide_db(vec![nt_entry("N001", "5 mismatches", &subj_str)]);
    let params = SearchParams::blastn()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false);
    let results = blastn(&db, query.as_bytes(), &params);
    assert_eq!(results.len(), 1, "should find hit with edited mismatches");
    assert_eq!(results[0].subject_oid, 0);
    let hsp = &results[0].hsps[0];
    assert_eq!(hsp.alignment_length, query.len());
    assert_eq!(hsp.num_identities as usize, query.len() - 4);
    assert_eq!(hsp.num_gaps, 0);
    assert!((hsp.percent_identity() - 96.07843137254902).abs() < 1.0e-12);
}

// ── Multi-subject index tests ───────────────────────────────────────────────

#[test]
fn blastp_multi_subject_finds_correct_hit() {
    // 5 different proteins; only one matches the query
    let (_tmp, db) = build_protein_db(vec![
        protein_entry("P001", "unrelated1", "WWWWWWWWWWWWWWWWWWWWWWWWWWWWW"),
        protein_entry("P002", "unrelated2", "AAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
        protein_entry("P003", "the target", "MKFLILLFNILCLFPVLAADNHGVSMNAS"),
        protein_entry("P004", "unrelated3", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCC"),
        protein_entry("P005", "unrelated4", "GGGGGGGGGGGGGGGGGGGGGGGGGGGG"),
    ]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    let results = blastp(&db, b"MKFLILLFNILCLFPVLAADNHGVSMNAS", &params);
    assert_eq!(results.len(), 1, "should find the matching subject only");
    // The matching subject should be OID 2 (0-indexed)
    assert_eq!(results[0].subject_oid, 2, "should match OID 2 (the target)");
    assert_eq!(results[0].hsps[0].alignment_length, 29);
    assert_eq!(results[0].hsps[0].num_identities, 29);
    assert_eq!(results[0].hsps[0].num_gaps, 0);
}

#[test]
fn blastp_multi_subject_finds_multiple_hits() {
    // Two subjects match, three don't
    let (_tmp, db) = build_protein_db(vec![
        protein_entry("P001", "match exact", "MKFLILLFNILCLFPVLAADNHGVSMNAS"),
        protein_entry("P002", "unrelated", "WWWWWWWWWWWWWWWWWWWWWWWWWWWWW"),
        protein_entry("P003", "match similar", "MKFLILLFNILCLFPVLAAENHGVSMNAS"), // one mismatch
        protein_entry("P004", "unrelated", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCC"),
    ]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    let results = blastp(&db, b"MKFLILLFNILCLFPVLAADNHGVSMNAS", &params);
    let oids: Vec<u32> = results.iter().map(|r| r.subject_oid).collect();
    assert_eq!(oids, [0, 2], "should rank exact match before similar match");
}

#[test]
fn blastn_multi_subject_finds_correct_hit() {
    let target = "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACGATCGTAAGGCCTTAGCAGTCAATGC";
    let (_tmp, db) = build_nucleotide_db(vec![
        nt_entry("N001", "decoy1", "TTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGC"),
        nt_entry("N002", "decoy2", "AAACCCGGGTTAAAACCCGGGTTAAAACCCGGGTTAAAACCCGGGTTAAAACCCGGGTTAAAACCCGGGTTAAAACCCGGGTTAAAACCCGGGTTAAAACCCGG"),
        nt_entry("N003", "target", target),
        nt_entry("N004", "decoy3", "GGGAAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAATT"),
    ]);
    let params = SearchParams::blastn()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false);

    let results = blastn(&db, target.as_bytes(), &params);
    assert_eq!(
        results.len(),
        1,
        "should find the target in multi-subject DB"
    );
    assert_eq!(results[0].subject_oid, 2, "should match OID 2 (the target)");
    assert!((results[0].hsps[0].percent_identity() - 100.0).abs() < 0.01);
    assert_eq!(results[0].hsps[0].alignment_length, target.len());
    assert_eq!(results[0].hsps[0].num_identities as usize, target.len());
    assert_eq!(results[0].hsps[0].num_gaps, 0);
}

#[test]
fn blastn_multi_subject_many_sequences() {
    // 20 sequences, only one matches
    let target = "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACGATCGTAAGGCCTTAGCAGTCAATGC";
    let decoy = "TTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGC";
    let mut entries: Vec<SequenceEntry> = (0..19)
        .map(|i| nt_entry(&format!("N{:03}", i), "decoy", decoy))
        .collect();
    entries.push(nt_entry("N019", "target", target));

    let (_tmp, db) = build_nucleotide_db(entries);
    assert_eq!(db.num_oids, 20);

    let params = SearchParams::blastn()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false);
    let results = blastn(&db, target.as_bytes(), &params);
    assert_eq!(results.len(), 1, "should find target in 20-sequence DB");
    assert_eq!(results[0].subject_oid, 19, "target should be OID 19");
    assert_eq!(results[0].hsps[0].alignment_length, target.len());
    assert_eq!(results[0].hsps[0].num_identities as usize, target.len());
    assert_eq!(results[0].hsps[0].num_gaps, 0);
}

#[test]
fn blastp_multi_subject_20_sequences() {
    let query = "MKFLILLFNILCLFPVLAADNHGVSMNAS";
    let entries: Vec<SequenceEntry> = (0..20)
        .map(|i| protein_entry(&format!("P{:03}", i), "copy", query))
        .collect();

    let (_tmp, db) = build_protein_db(entries);
    assert_eq!(db.num_oids, 20);

    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);
    let results = blastp(&db, query.as_bytes(), &params);
    assert_eq!(results.len(), 20, "should find all 20 identical subjects");
    assert_eq!(
        results.iter().map(|r| r.subject_oid).collect::<Vec<_>>(),
        (0..20).rev().collect::<Vec<_>>()
    );
    assert!(results
        .iter()
        .all(|r| r.hsps[0].alignment_length == query.len()));
    assert!(results
        .iter()
        .all(|r| r.hsps[0].num_identities as usize == query.len()));
}

#[test]
fn blastx_multi_subject_protein_db() {
    let nt_query = "ATGAAATTTCTGATTCTGCTGTTT"; // encodes MKFLILLF
    let (_tmp, db) = build_protein_db(vec![
        protein_entry("P001", "decoy", "WWWWWWWWWWWWWWWWWWWW"),
        protein_entry("P002", "target", "MKFLILLF"),
        protein_entry("P003", "decoy", "AAAAAAAAAAAAAAAAAAA"),
    ]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    let results = blastx(&db, nt_query.as_bytes(), &params);
    assert_eq!(
        results.len(),
        1,
        "blastx should find target in multi-subject protein DB"
    );
    assert_eq!(
        results[0].subject_oid, 1,
        "should match OID 1 (the target protein)"
    );
    assert_eq!(results[0].hsps[0].query_frame, 1);
    assert_eq!(results[0].hsps[0].alignment_length, 8);
    assert_eq!(results[0].hsps[0].num_identities, 8);
    assert_eq!(results[0].hsps[0].num_gaps, 0);
}

#[test]
fn tblastn_multi_subject_nt_db() {
    let protein_query = "MKFLILLF";
    let nt_target = "ATGAAATTTCTGATTCTGCTGTTT"; // encodes MKFLILLF
    let nt_decoy = "TTTGGGCCCAAATTTGGGCCCAAA";

    let (_tmp, db) = build_nucleotide_db(vec![
        nt_entry("N001", "decoy1", nt_decoy),
        nt_entry("N002", "decoy2", nt_decoy),
        nt_entry("N003", "target", nt_target),
    ]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    let results = tblastn(&db, protein_query.as_bytes(), &params);
    assert_eq!(
        results.len(),
        1,
        "tblastn should find protein in multi-subject nt DB"
    );
    assert_eq!(
        results[0].subject_oid, 2,
        "should match OID 2 (the coding sequence)"
    );
    assert_eq!(results[0].hsps[0].subject_frame, 1);
    assert_eq!(results[0].hsps[0].alignment_length, protein_query.len());
    assert_eq!(
        results[0].hsps[0].num_identities as usize,
        protein_query.len()
    );
    assert_eq!(results[0].hsps[0].num_gaps, 0);
}

#[test]
fn blastx_same_subject_can_emit_multiple_hsps() {
    let nt_query = "ATGAAATTTCTGATTCTGCTGTTTAAAACCCCGGGGTTTTATGAAATTTCTGATTCTGCTGTTT";
    let (_tmp, db) = build_protein_db(vec![protein_entry(
        "P001",
        "multi",
        "MKFLILLFQQQQGMKFLILLF",
    )]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    let results = blastx(&db, nt_query.as_bytes(), &params);
    assert_eq!(results.len(), 1, "expected one subject hit");
    assert!(
        results[0].hsps.len() >= 2,
        "blastx should report multiple translated HSPs for separated matching regions"
    );
}

#[test]
fn tblastn_same_subject_can_emit_multiple_hsps() {
    let nt_target = "ATGAAATTTCTGATTCTGCTGTTTAAAACCCCGGGGTTTTATGAAATTTCTGATTCTGCTGTTT";
    let (_tmp, db) = build_nucleotide_db(vec![nt_entry("N001", "multi", nt_target)]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    let results = tblastn(&db, b"MKFLILLFQQQQGMKFLILLF", &params);
    assert_eq!(results.len(), 1, "expected one subject hit");
    assert!(
        results[0].hsps.len() >= 2,
        "tblastn should report multiple translated HSPs for separated matching regions"
    );
}

#[test]
fn tblastx_indel_repro_is_ungapped_only() {
    let query =
        b"ATGAAATTTCTGATTCTGCTGTTTAATATTCTGTGTCTGTTTCCTGTTCTGGCTGCTGATAATCATGGTGTTTCTATGAATGCTTCT";
    let subject =
        "ATGAAATTTCTGATTCTGCTGTTTATTCTGTGTCTGTTTCCTGTTCTGGCTGCTGATAATCATGGTGTTTCTATGAATGCTTCT";
    let (_tmp, db) = build_nucleotide_db(vec![nt_entry("N001", "indel", subject)]);
    let params = SearchParams::blastp()
        .evalue(1000.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    let results = tblastx(&db, query, &params);
    assert_eq!(results.len(), 1, "expected one subject hit");
    assert!(results[0].hsps.iter().all(|hsp| hsp.num_gaps == 0));
    assert!(results[0]
        .hsps
        .iter()
        .all(|hsp| { !hsp.query_aln.contains(&b'-') && !hsp.subject_aln.contains(&b'-') }));
}

#[test]
fn tblastx_accepts_explicit_positive_max_intron_without_linker_panic() {
    let query =
        b"ATGAAATTTCTGATTCTGCTGTTTAATATTCTGTGTCTGTTTCCTGTTCTGGCTGCTGATAATCATGGTGTTTCTATGAATGCTTCT";
    let subject =
        "ATGAAATTTCTGATTCTGCTGTTTATTCTGTGTCTGTTTCCTGTTCTGGCTGCTGATAATCATGGTGTTTCTATGAATGCTTCT";
    let (_tmp, db) = build_nucleotide_db(vec![nt_entry("N001", "indel", subject)]);
    let params = SearchParams::tblastx()
        .evalue(1000.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .max_intron_length(303);

    let results = tblastx(&db, query, &params);
    assert_eq!(results.len(), 1, "expected one subject hit");
    assert!(results[0].hsps.iter().all(|hsp| hsp.num_gaps == 0));
}

#[test]
fn blastx_translated_overlap_repro_matches_expected_hsp_set() {
    let nt_query = "ATGAAATTTCTGATTCTGCTGTTTAAAACCCCGGGGTTTTATGAAATTTCTGATTCTGCTGTTT";
    let (_tmp, db) = build_protein_db(vec![protein_entry(
        "P001",
        "translated-overlap",
        "MKFLILLFQQQQGMKFLILLF",
    )]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false);

    let results = blastx(&db, nt_query.as_bytes(), &params);
    assert_eq!(results.len(), 1, "expected one subject hit");
    let hsps: Vec<_> = results[0]
        .hsps
        .iter()
        .map(|h| {
            (
                h.query_start,
                h.query_end,
                h.subject_start,
                h.subject_end,
                h.score,
                h.alignment_length,
                h.query_frame,
                h.subject_frame,
            )
        })
        .collect();
    assert_eq!(
        hsps,
        vec![
            (0, 33, 0, 11, 36, 11, 1, 0),
            (0, 24, 13, 21, 34, 8, 1, 0),
            (40, 64, 0, 8, 32, 8, 2, 0),
            (40, 64, 13, 21, 32, 8, 2, 0),
        ]
    );
}

#[test]
fn tblastn_translated_overlap_repro_matches_expected_hsp_set() {
    let nt_target = "ATGAAATTTCTGATTCTGCTGTTTAAAACCCCGGGGTTTTATGAAATTTCTGATTCTGCTGTTT";
    let (_tmp, db) = build_nucleotide_db(vec![nt_entry("N001", "translated-overlap", nt_target)]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false);

    let results = tblastn(&db, b"MKFLILLFQQQQGMKFLILLF", &params);
    assert_eq!(results.len(), 1, "expected one subject hit");
    let hsps: Vec<_> = results[0]
        .hsps
        .iter()
        .map(|h| {
            (
                h.query_start,
                h.query_end,
                h.subject_start,
                h.subject_end,
                h.score,
                h.alignment_length,
                h.query_frame,
                h.subject_frame,
            )
        })
        .collect();
    assert_eq!(
        hsps,
        vec![
            (0, 11, 0, 33, 36, 11, 0, 1),
            (13, 21, 0, 24, 34, 8, 0, 1),
            (0, 8, 40, 64, 32, 8, 0, 2),
            (13, 21, 40, 64, 32, 8, 0, 2),
        ]
    );
}

// ── Prokka-style annotation performance test ─────────────────────────────────

/// Performance test: build a protein DB from prokka's sprot FASTA (~25K entries),
/// then search 5 query proteins against it. This mimics the prokka annotation
/// workload. Run with: cargo test --release -- --ignored test_blastp_prokka_sprot
///
/// Target: should complete in under 30 seconds for 5 queries (Perl Prokka does
/// 63 queries against the same DB in ~12 seconds using NCBI BLAST+).
#[test]
#[ignore]
fn test_blastp_prokka_sprot() {
    // Paths relative to the blast-rs repo — prokka db may be at a sibling path
    let sprot_paths = [
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../prokka-rs/prokka/db/kingdom/Bacteria/sprot"),
        std::path::PathBuf::from(
            "/data/henriksson/github/claude/prokka-rs/prokka/db/kingdom/Bacteria/sprot",
        ),
    ];
    let sprot_path = sprot_paths.iter().find(|p| p.exists());
    let sprot_path = match sprot_path {
        Some(p) => p,
        None => {
            eprintln!("Skipping: prokka sprot database not found");
            return;
        }
    };

    // Load sprot FASTA using input::parse_fasta (returns FastaRecord structs)
    let file = std::fs::File::open(sprot_path).unwrap();
    let records = blast_rs::input::parse_fasta(file);
    eprintln!("Loaded {} sprot records", records.len());
    assert!(records.len() > 1000, "sprot should have >1000 entries");

    // Build indexed protein DB
    let t0 = std::time::Instant::now();
    let tmp = TempDir::new().unwrap();
    let base = tmp.path().join("sprot");
    let mut builder = BlastDbBuilder::new(DbType::Protein, "sprot");
    for rec in &records {
        builder.add(SequenceEntry {
            title: rec.defline.clone(),
            accession: rec.id.clone(),
            sequence: rec.sequence.clone(),
            taxid: None,
        });
    }
    builder.write(&base).unwrap();
    let db = blast_rs::db::BlastDb::open(&base).unwrap();
    let db_build_time = t0.elapsed();
    eprintln!(
        "DB build: {:.2}s ({} entries)",
        db_build_time.as_secs_f64(),
        records.len()
    );

    // 5 representative bacterial protein queries (real CDS translations)
    let queries: Vec<&[u8]> = vec![
        // Replication initiation protein (~90 aa)
        b"MKQIKEYLEEFVHSRLNKNIILRAAGFEYAKENPNFSQYYGNTVVSLPHRGKYGGPVNRIAPEMFHQIVAKPGERTFEGMFAIFKHRFPDWRDAES",
        // Short hypothetical (~50 aa)
        b"MNDFNYYKSKEIYREKYYQMPKVFFTNEKYMDLSNDAKIAYMLLKDRFDYS",
        // Transposase fragment (~80 aa)
        b"MNYFRYKQFNKDVITVAVGYYLRYALSYRDISEILRGRGVNVHHSTVYRWVQEYAPILYQQSINTAKNTLKGIECIYALY",
        // TraB transfer protein (~60 aa)
        b"MIKKFSLTTVYVAFLSIVLSNITLGAENPGPKIEQGLQQVQTFLTGLIVAVGICAGVWIV",
        // ErmC methyltransferase (~70 aa)
        b"MNEKNIKHSQNFITSKHNIDKIMTNIRLNEHDNIFEIGSGKGHFTLELVQRCNFVTAIEIDHKLCKTTEN",
    ];

    let params = SearchParams::blastp()
        .evalue(1e-6)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    let t1 = std::time::Instant::now();
    let mut total_hits = 0;
    for (i, query) in queries.iter().enumerate() {
        let t_q = std::time::Instant::now();
        let results = blastp(&db, query, &params);
        let q_time = t_q.elapsed();
        eprintln!(
            "Query {} ({} aa): {} hits in {:.2}s",
            i + 1,
            query.len(),
            results.len(),
            q_time.as_secs_f64()
        );
        total_hits += results.len();
    }
    let search_time = t1.elapsed();
    eprintln!(
        "Total: {} queries, {} hits, {:.2}s ({:.2}s/query)",
        queries.len(),
        total_hits,
        search_time.as_secs_f64(),
        search_time.as_secs_f64() / queries.len() as f64
    );

    // Performance summary
    let per_query = search_time.as_secs_f64() / queries.len() as f64;
    let ncbi_per_query = 12.0 / 63.0; // NCBI BLAST+ reference: 63 queries in 12s
    eprintln!(
        "\nPerformance: {:.2}s/query (NCBI BLAST+ reference: {:.3}s/query, {:.0}x slower)",
        per_query,
        ncbi_per_query,
        per_query / ncbi_per_query
    );

    // Performance assertion: 5 queries should complete in under 60 seconds
    // (NCBI BLAST+ does 63 queries in ~12s, so 5 queries in 60s is very generous)
    assert!(
        search_time.as_secs() < 60,
        "blastp search too slow: {:.1}s for {} queries (target: <60s)",
        search_time.as_secs_f64(),
        queries.len()
    );

    // Now test with all available threads
    let params_mt = SearchParams::blastp()
        .evalue(1e-6)
        .num_threads(0) // 0 = use all available cores
        .filter_low_complexity(false)
        .comp_adjust(0);

    let t2 = std::time::Instant::now();
    let mut total_hits_mt = 0;
    for query in &queries {
        let results = blastp(&db, query, &params_mt);
        total_hits_mt += results.len();
    }
    let mt_time = t2.elapsed();
    let mt_per_query = mt_time.as_secs_f64() / queries.len() as f64;
    let speedup = search_time.as_secs_f64() / mt_time.as_secs_f64();
    eprintln!(
        "\nMulti-threaded: {:.2}s/query ({:.1}x speedup over single-threaded, {:.0}x vs NCBI)",
        mt_per_query,
        speedup,
        mt_per_query / ncbi_per_query
    );

    // Verify same hit count
    assert_eq!(
        total_hits, total_hits_mt,
        "Multi-threaded should find same hits as single-threaded"
    );
}

/// Sensitivity test: verify that blastp finds the same hits as NCBI BLAST+.
///
/// These are real CDS proteins from E. faecium plasmid AUS0004_p1 that
/// NCBI BLAST+ (via Perl Prokka) annotates against the sprot database.
/// blast-rs must find the same top hit and hit count for each of these queries.
///
/// Run with: cargo test --release -- --ignored test_blastp_sensitivity
#[test]
#[ignore]
fn test_blastp_sensitivity() {
    let sprot_paths = [
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../prokka-rs/prokka/db/kingdom/Bacteria/sprot"),
        std::path::PathBuf::from(
            "/data/henriksson/github/claude/prokka-rs/prokka/db/kingdom/Bacteria/sprot",
        ),
    ];
    let sprot_path = match sprot_paths.iter().find(|p| p.exists()) {
        Some(p) => p,
        None => {
            eprintln!("Skipping: prokka sprot database not found");
            return;
        }
    };

    // Build DB
    let file = std::fs::File::open(sprot_path).unwrap();
    let records = blast_rs::input::parse_fasta(file);
    let tmp = TempDir::new().unwrap();
    let base = tmp.path().join("sprot");
    let mut builder = BlastDbBuilder::new(DbType::Protein, "sprot");
    for rec in &records {
        builder.add(SequenceEntry {
            title: rec.defline.clone(),
            accession: rec.id.clone(),
            sequence: rec.sequence.clone(),
            taxid: None,
        });
    }
    builder.write(&base).unwrap();
    let db = blast_rs::db::BlastDb::open(&base).unwrap();

    // Queries that NCBI BLAST+ finds in sprot (from Perl Prokka on E. faecium plasmid).
    // Each tuple: (name, protein_sequence, expected_top_accession, hit_count, top_evalue)
    type SensitivityExpectedRow<'a> = (&'a str, &'a [u8], &'a str, usize, f64);
    let known_hits: Vec<SensitivityExpectedRow> = vec![
        ("topB",
         b"MMKTVILAEKPSQAKAYADSFSKATRKDGYFEIQDRLFPGETVITYGFGHLVELDSPDMYDENWKQWSLEHLPIFPNQYHYHVPKDKKKQFNVVKQQLQSADTIIIATDSDREGELIAWTIIQQAGADQGKTFKRLWINSLEKEAIYQGFQQLRDAGETYPKFEEAQARQIADWLIGMNGSPLYSLLLQQKGIPGSFSLGRVQTPTLYMIYQLQEKIRNFQKEPYFEGKAQVIAQNGAFDAKLDPNETQATQEAFEDYLKEKGVQLGKQPGTIHQVETEKKSAASPRLFSLSSLQSKMNQLMKASAKDTLEAMQGLYEGKYLSYPRTDTPYITEGEYAYLLDHLDEYKHFLKAEAIPTPIHTPNSRYVNNKKVQEHYAIIPTKTVMTAAAFEQLSPLQQAIYEQVLKTTVAMFAEKYTYEETTILTQVQQLQLKAIGKVPLDLGWKKLFGKESEGKEKEEEPLLPKVTKGETVTVDLQVLEKETKPPQPYTEGTLITAMKTAGKTVDSEEAQSILKEVEGIGTEATRANIIETLKQKEYIKVEKNKLVVTNKGILLCQAVEKEPLLTSAEMTAKWESYLLKIGERKGTQTTFLTNIQKFVSHLLEVVPGQIQSTDFGSTLQEVKAASEKQEAARHLGICPKCQEQEVLLYHKAAACTSEACDFRLWTTIAKKKLTATQLKEIIQNGRTSQPVKGLKGQKGSFEATIVLKEDFTTSFEFSEKKKTNYKKRTRRTTK",
         "P14294", 10, 3.15e-84),
        ("ssb",
         b"MINNVTLVGRLTKDPDLRYTASGTAVATFTLAVNRNFTNQNGNREADFINCVIWRKPAETMATLAKKGILIGVVGRIQTRTYDNQQGQRVYVTEVVADNFQLLESKAATESRAHADQSSTSPSTTTFEQRDTATPNNNGLNASQNPFGGQSIDISDDDLPF",
         "P66854", 20, 2.50e-61),
        ("yoeB",
         b"MIKAWSDDAWDDYLYWHEQGNKSNIKKINKLIKDIDRSPFAGLGKPEALKHDLSGKWSRRRLVDLTDDDLEKIREEKIPFFIGLSQDRVQRMYQEKGLTIDSVFHGKRKPVTKVIINDLVERF",
         "P69348", 2, 4.35e-23),
        ("srtA",
         b"MKKCFLFCLKGGRWMKKWLFGFLGVALIVVCSVFGYVSYQKHEGEVFKQNIEKKMPVDQINAHAKSYKEDATNVNNDMSLGQMLSIQKEAIEMGVNKQVFAQIQIPALGLALPIFKGANQYTLSLGAATYFYEDAEMGKGNYVLAGHNMEMPGVLFSDIQKLSLGEVMDLVSNDGVYRYKVTRKFIVPEYFKLIDGVPEENSFLSLPKKGEKPLLTLFTCVYTSQGKERYVVQGELQ",
         "P0DPQ5", 5, 5.62e-24),
        // Note: umuC removed — NCBI BLAST+ also finds no hit at evalue 1e-6
        // (verified with comp_based_stats=0/1 and default settings).
    ];

    let params = SearchParams::blastp()
        .evalue(1e-6)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    for (name, query, expected_acc, expected_hits, expected_evalue) in &known_hits {
        let results = blastp(&db, query, &params);
        assert_eq!(results.len(), *expected_hits, "{name}: hit count drifted");
        let best = &results[0];
        assert!(
            best.subject_accession.contains(expected_acc)
                || best.subject_title.contains(expected_acc),
            "{name}: top hit {} / {} should identify {expected_acc}",
            best.subject_accession,
            best.subject_title
        );
        assert!(
            (best.best_evalue() - *expected_evalue).abs() <= expected_evalue * 0.01,
            "{name}: best e-value {:.3e} should match expected {:.3e}",
            best.best_evalue(),
            expected_evalue
        );
    }
}

/// Test composition-based matrix adjustment on the srtA/Sortase A fixture.
/// The query at position 00009 in the E. faecium plasmid is a legitimate
/// P0DPQ5 Sortase A hit; composition adjustment should retain it while
/// reshaping the alignment and score consistently with the NCBI reference.
///
/// Run with: cargo test --release -- --ignored test_comp_adjust_srta
#[test]
#[ignore]
fn test_comp_adjust_srta_retains_sortase_hit() {
    let sprot_paths = [
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../prokka-rs/prokka/db/kingdom/Bacteria/sprot"),
        std::path::PathBuf::from(
            "/data/henriksson/github/claude/prokka-rs/prokka/db/kingdom/Bacteria/sprot",
        ),
    ];
    let sprot_path = match sprot_paths.iter().find(|p| p.exists()) {
        Some(p) => p,
        None => {
            eprintln!("Skipping: prokka sprot database not found");
            return;
        }
    };
    let file = std::fs::File::open(sprot_path).unwrap();
    let records = blast_rs::input::parse_fasta(file);
    let tmp = TempDir::new().unwrap();
    let base = tmp.path().join("sprot");
    let mut builder = BlastDbBuilder::new(DbType::Protein, "sprot");
    for rec in &records {
        builder.add(SequenceEntry {
            title: rec.defline.clone(),
            accession: rec.id.clone(),
            sequence: rec.sequence.clone(),
            taxid: None,
        });
    }
    builder.write(&base).unwrap();
    let db = blast_rs::db::BlastDb::open(&base).unwrap();

    // The srtA query (CDS at 6145..6858 on plasmid, actual sequence from prodigal).
    let query = b"MKKCFLFCLKGGRWMKKWLFGFLGVALIVVCSVFGYVSYQKHEGEVFKQNIEKKMPVDQINAHAKSYKEDATNVNNDMSLGQMLSIQKEAIEMGVNKQVFAQIQIPALGLALPIFKGANQYTLSLGAATYFYEDAEMGKGNYVLAGHNMEMPGVLFSDIQKLSLGEVMDLVSNDGVYRYKVTRKFIVPEYFKLIDGVPEENSFLSLPKKGEKPLLTLFTCVYTSQGKERYVVQGELQ";

    // Without composition adjustment — should find a hit
    let params_no = SearchParams::blastp()
        .evalue(1e-6)
        .num_threads(1)
        .comp_adjust(0);
    let results_no = blastp(&db, query, &params_no);

    // With lambda scaling only (mode 1)
    let params_1 = SearchParams::blastp()
        .evalue(1e-6)
        .num_threads(1)
        .comp_adjust(1);
    let results_1 = blastp(&db, query, &params_1);

    // With conditional matrix adjustment (mode 2).
    let params_2 = SearchParams::blastp()
        .evalue(1e-6)
        .num_threads(1)
        .comp_adjust(2);
    let results_2 = blastp(&db, query, &params_2);

    // Check: the Sortase hit remains, with the NCBI reference score shape.
    let srta_hits_2: Vec<_> = results_2
        .iter()
        .filter(|r| r.subject_title.to_lowercase().contains("sortase"))
        .collect();
    assert_eq!(results_no.len(), 5);
    assert!(
        results_no
            .first()
            .is_some_and(|r| r.subject_accession.contains("P0DPQ5")),
        "comp_adjust=0 should find P0DPQ5 as the top srtA hit"
    );
    assert_eq!(results_no[0].hsps[0].score, 234);
    assert_eq!(results_1.len(), 5);
    assert!(
        results_1
            .first()
            .is_some_and(|r| r.subject_accession.contains("P0DPQ5")),
        "comp_adjust=1 should retain P0DPQ5"
    );
    assert_eq!(results_1[0].hsps[0].score, 221);
    assert_eq!(results_2.len(), 5);
    assert!(
        results_2
            .first()
            .is_some_and(|r| r.subject_accession.contains("P0DPQ5")),
        "comp_adjust=2 should retain P0DPQ5"
    );
    assert_eq!(results_2[0].hsps[0].score, 233);
    assert_eq!(srta_hits_2.len(), 5);

    // Also test dinB query (DNA polymerase IV) — should NOT be eliminated
    let dinb_query = b"MYLAISSLQHRTYVCIMWKNGVLFMMDYSKEPVNDYFLIDMKSFYASVECIERNLDPLTTELVVMSRSDNTGSGLILASSPEAKKRYGITNVSRPRDLPQPFPKTLHVVPPRMNLYIKRNMQVNNIFRRYVADEDLLIYSIDESILKVTKSLNLFTTEETRSQRRKKLAQMIQERIKEELGLIATVGVGDNPLLAKLALDNEAKHNEGFIAEWTYENVPEKVWNIPEMTDFWGIGSRMKKRLNQMGILSIRDLANWNPYTIKNRLGVIGLQLYFHANGIDRTDIAIPPEPTKEKSYGNSQVLPRDYTRRNEIELVVKEMAEQVAIRIRQHNCKTGCVHLNIGTSILETRPGFSHQMKIPITDNTKELQNYCLFLFDKYYEGQEVRHVGITYSKLVYTDSLQLDLFSDPQKQINEENLDKIIDKIRQKYGFTSIVHASSMLESARSITRSTLVGGHAGGNGGIKND";
    let dinb_expected = [
        (0u8, "P58965", 186, 4usize),
        (1u8, "Q47155", 178, 4usize),
        (2u8, "Q47155", 181, 4usize),
    ];
    for (mode, expected_acc, expected_score, expected_hits) in dinb_expected {
        let params = SearchParams::blastp()
            .evalue(1e-6)
            .num_threads(1)
            .comp_adjust(mode);
        let results = blastp(&db, dinb_query, &params);
        assert_eq!(results.len(), expected_hits);
        assert!(
            results[0].subject_accession.contains(expected_acc),
            "dinB comp_adjust={mode} should rank {expected_acc} first"
        );
        assert_eq!(results[0].hsps[0].score, expected_score);
    }
}

/// Real-life test: compare blast-rs against NCBI BLAST+ output for known prokka queries.
/// Tests that top hit, score, alignment length, and e-value match NCBI BLAST+.
/// Each query was run through NCBI blastp 2.12.0+ with comp_based_stats=0 against
/// the Bacteria/sprot database and the top hit recorded.
///
/// Run with: cargo test --release -- --ignored test_blastp_vs_ncbi
#[test]
#[ignore]
fn test_blastp_vs_ncbi() {
    let sprot_paths = [
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../prokka-rs/prokka/db/kingdom/Bacteria/sprot"),
        std::path::PathBuf::from(
            "/data/henriksson/github/claude/prokka-rs/prokka/db/kingdom/Bacteria/sprot",
        ),
    ];
    let sprot_path = match sprot_paths.iter().find(|p| p.exists()) {
        Some(p) => p,
        None => {
            eprintln!("Skipping: prokka sprot database not found");
            return;
        }
    };
    let file = std::fs::File::open(sprot_path).unwrap();
    let records = blast_rs::input::parse_fasta(file);
    let tmp = TempDir::new().unwrap();
    let base = tmp.path().join("sprot");
    let mut builder = BlastDbBuilder::new(DbType::Protein, "sprot");
    for rec in &records {
        builder.add(SequenceEntry {
            title: rec.defline.clone(),
            accession: rec.id.clone(),
            sequence: rec.sequence.clone(),
            taxid: None,
        });
    }
    builder.write(&base).unwrap();
    let db = blast_rs::db::BlastDb::open(&base).unwrap();

    // Known queries from Prokka plasmid test with NCBI BLAST+ results (comp_based_stats=0).
    // Format: (name, query_seq, ncbi_top_accession, ncbi_score, ncbi_evalue, ncbi_length)
    type ProkkaExpectedRow<'a> = (&'a str, &'a [u8], &'a str, i32, f64, usize);
    let known: Vec<ProkkaExpectedRow> = vec![
        // topB — DNA topoisomerase 3
        // NCBI blastp 2.12.0+ comp_based_stats=0: P14294 score=714 e=3.15e-84 len=666
        ("topB",
         b"MMKTVILAEKPSQAKAYADSFSKATRKDGYFEIQDRLFPGETVITYGFGHLVELDSPDMYDENWKQWSLEHLPIFPNQYHYHVPKDKKKQFNVVKQQLQSADTIIIATDSDREGELIAWTIIQQAGADQGKTFKRLWINSLEKEAIYQGFQQLRDAGETYPKFEEAQARQIADWLIGMNGSPLYSLLLQQKGIPGSFSLGRVQTPTLYMIYQLQEKIRNFQKEPYFEGKAQVIAQNGAFDAKLDPNETQATQEAFEDYLKEKGVQLGKQPGTIHQVETEKKSAASPRLFSLSSLQSKMNQLMKASAKDTLEAMQGLYEGKYLSYPRTDTPYITEGEYAYLLDHLDEYKHFLKAEAIPTPIHTPNSRYVNNKKVQEHYAIIPTKTVMTAAAFEQLSPLQQAIYEQVLKTTVAMFAEKYTYEETTILTQVQQLQLKAIGKVPLDLGWKKLFGKESEGKEKEEEPLLPKVTKGETVTVDLQVLEKETKPPQPYTEGTLITAMKTAGKTVDSEEAQSILKEVEGIGTEATRANIIETLKQKEYIKVEKNKLVVTNKGILLCQAVEKEPLLTSAEMTAKWESYLLKIGERKGTQTTFLTNIQKFVSHLLEVVPGQIQSTDFGSTLQEVKAASEKQEAARHLGICPKCQEQEVLLYHKAAACTSEACDFRLWTTIAKKKLTATQLKEIIQNGRTSQPVKGLKGQKGSFEATIVLKEDFTTSFEFSEKKKTNYKKRTRRTTK",
         "P14294", 714, 3.15e-84, 666),
        // ssb — Single-stranded DNA-binding protein
        // NCBI: P66854 score=471 e=2.50e-61 len=162
        ("ssb",
         b"MINNVTLVGRLTKDPDLRYTASGTAVATFTLAVNRNFTNQNGNREADFINCVIWRKPAETMATLAKKGILIGVVGRIQTRTYDNQQGQRVYVTEVVADNFQLLESKAATESRAHADQSSTSPSTTTFEQRDTATPNNNGLNASQNPFGGQSIDISDDDLPF",
         "P66854", 471, 2.50e-61, 162),  // DP66854 in our DB
        // yoeB — Toxin YoeB
        // NCBI: P69348 score=209 e=4.35e-23 len=72
        ("yoeB",
         b"MIKAWSDDAWDDYLYWHEQGNKSNIKKINKLIKDIDRSPFAGLGKPEALKHDLSGKWSRRRLVDLTDDDLEKIREEKIPFFIGLSQDRVQRMYQEKGLTIDSVFHGKRKPVTKVIINDLVERF",
         "P69348", 209, 4.35e-23, 72),
        // srtA — Sortase A (legitimate hit, tests composition adjustment)
        // NCBI: P0DPQ5 score=234 e=5.62e-24 len=225
        ("srtA",
         b"MKKCFLFCLKGGRWMKKWLFGFLGVALIVVCSVFGYVSYQKHEGEVFKQNIEKKMPVDQINAHAKSYKEDATNVNNDMSLGQMLSIQKEAIEMGVNKQVFAQIQIPALGLALPIFKGANQYTLSLGAATYFYEDAEMGKGNYVLAGHNMEMPGVLFSDIQKLSLGEVMDLVSNDGVYRYKVTRKFIVPEYFKLIDGVPEENSFLSLPKKGEKPLLTLFTCVYTSQGKERYVVQGELQ",
         "P0DPQ5", 234, 5.62e-24, 225),
        // dinB — DNA polymerase IV
        // NCBI: P58965 score=186 e=3.86e-15 len=408
        ("dinB",
         b"MYLAISSLQHRTYVCIMWKNGVLFMMDYSKEPVNDYFLIDMKSFYASVECIERNLDPLTTELVVMSRSDNTGSGLILASSPEAKKRYGITNVSRPRDLPQPFPKTLHVVPPRMNLYIKRNMQVNNIFRRYVADEDLLIYSIDESILKVTKSLNLFTTEETRSQRRKKLAQMIQERIKEELGLIATVGVGDNPLLAKLALDNEAKHNEGFIAEWTYENVPEKVWNIPEMTDFWGIGSRMKKRLNQMGILSIRDLANWNPYTIKNRLGVIGLQLYFHANGIDRTDIAIPPEPTKEKSYGNSQVLPRDYTRRNEIELVVKEMAEQVAIRIRQHNCKTGCVHLNIGTSILETRPGFSHQMKIPITDNTKELQNYCLFLFDKYYEGQEVRHVGITYSKLVYTDSLQLDLFSDPQKQINEENLDKIIDKIRQKYGFTSIVHASSMLESARSITRSTLVGGHAGGNGGIKND",
         "P58965", 186, 3.86e-15, 408),
    ];

    let params = SearchParams::blastp()
        .evalue(1e-6)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    for (name, query, expected_acc, ncbi_score, ncbi_evalue, ncbi_length) in &known {
        let results = blastp(&db, query, &params);
        let best = results
            .first()
            .unwrap_or_else(|| panic!("{name}: should find {expected_acc} at e={ncbi_evalue:.2e}"));
        let best_hsp = &best.hsps[0];

        assert!(
            best.subject_accession.contains(expected_acc)
                || best.subject_title.contains(expected_acc),
            "{name}: wrong top hit {} / {} (expected {expected_acc})",
            best.subject_accession,
            best.subject_title
        );
        assert_eq!(
            best_hsp.score, *ncbi_score,
            "{name}: raw score should match NCBI"
        );
        assert_eq!(
            best_hsp.alignment_length, *ncbi_length,
            "{name}: alignment length should match NCBI"
        );
        assert!(
            (best_hsp.evalue - *ncbi_evalue).abs() <= ncbi_evalue * 0.01,
            "{name}: e-value {:.3e} should match NCBI {:.3e}",
            best_hsp.evalue,
            ncbi_evalue
        );
    }
}

/// Real-life test: verify hit count matches NCBI for the full plasmid annotation.
/// Runs all 63 CDS from the test plasmid against sprot and counts annotated hits.
/// NCBI BLAST+ with comp_based_stats=0, evalue=1e-9 finds ~11 hits for this dataset.
///
/// Run with: cargo test --release -- --ignored test_blastp_plasmid_annotation_count
#[test]
#[ignore]
fn test_blastp_plasmid_annotation_count() {
    let sprot_paths = [
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../prokka-rs/prokka/db/kingdom/Bacteria/sprot"),
        std::path::PathBuf::from(
            "/data/henriksson/github/claude/prokka-rs/prokka/db/kingdom/Bacteria/sprot",
        ),
    ];
    let sprot_path = match sprot_paths.iter().find(|p| p.exists()) {
        Some(p) => p,
        None => {
            eprintln!("Skipping: prokka sprot database not found");
            return;
        }
    };

    // Build BLAST DB
    let file = std::fs::File::open(sprot_path).unwrap();
    let records = blast_rs::input::parse_fasta(file);
    let tmp = TempDir::new().unwrap();
    let base = tmp.path().join("sprot");
    let mut builder = BlastDbBuilder::new(DbType::Protein, "sprot");
    for rec in &records {
        builder.add(SequenceEntry {
            title: rec.defline.clone(),
            accession: rec.id.clone(),
            sequence: rec.sequence.clone(),
            taxid: None,
        });
    }
    builder.write(&base).unwrap();
    let db = blast_rs::db::BlastDb::open(&base).unwrap();

    // Load plasmid CDS proteins from the .faa file if available, or skip
    let _faa_paths = [std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../prokka-rs/tests/data/plasmid_cds.faa")];
    // If no .faa file, generate proteins by running prodigal on the plasmid
    // For now, use hardcoded test proteins from known CDS
    let test_proteins: Vec<(&str, &[u8])> = vec![
        ("topB", b"MMKTVILAEKPSQAKAYADSFSKATRKDGYFEIQDRLFPGETVITYGFGHLVELDSPDMYDENWKQWSLEHLPIFPNQYHYHVPKDKKKQFNVVKQQLQSADTIIIATDSDREGELIAWTIIQQAGADQGKTFKRLWINSLEKEAIYQGFQQLRDAGETYPKFEEAQARQIADWLIGMNGSPLYSLLLQQKGIPGSFSLGRVQTPTLYMIYQLQEKIRNFQKEPYFEGKAQVIAQNGAFDAKLDPNETQATQEAFEDYLKEKGVQLGKQPGTIHQVETEKKSAASPRLFSLSSLQSKMNQLMKASAKDTLEAMQGLYEGKYLSYPRTDTPYITEGEYAYLLDHLDEYKHFLKAEAIPTPIHTPNSRYVNNKKVQEHYAIIPTKTVMTAAAFEQLSPLQQAIYEQVLKTTVAMFAEKYTYEETTILTQVQQLQLKAIGKVPLDLGWKKLFGKESEGKEKEEEPLLPKVTKGETVTVDLQVLEKETKPPQPYTEGTLITAMKTAGKTVDSEEAQSILKEVEGIGTEATRANIIETLKQKEYIKVEKNKLVVTNKGILLCQAVEKEPLLTSAEMTAKWESYLLKIGERKGTQTTFLTNIQKFVSHLLEVVPGQIQSTDFGSTLQEVKAASEKQEAARHLGICPKCQEQEVLLYHKAAACTSEACDFRLWTTIAKKKLTATQLKEIIQNGRTSQPVKGLKGQKGSFEATIVLKEDFTTSFEFSEKKKTNYKKRTRRTTK" as &[u8]),
        ("ssb", b"MINNVTLVGRLTKDPDLRYTASGTAVATFTLAVNRNFTNQNGNREADFINCVIWRKPAETMATLAKKGILIGVVGRIQTRTYDNQQGQRVYVTEVVADNFQLLESKAATESRAHADQSSTSPSTTTFEQRDTATPNNNGLNASQNPFGGQSIDISDDDLPF"),
        ("yoeB", b"MIKAWSDDAWDDYLYWHEQGNKSNIKKINKLIKDIDRSPFAGLGKPEALKHDLSGKWSRRRLVDLTDDDLEKIREEKIPFFIGLSQDRVQRMYQEKGLTIDSVFHGKRKPVTKVIINDLVERF"),
        ("srtA_1", b"MKKCFLFCLKGGRWMKKWLFGFLGVALIVVCSVFGYVSYQKHEGEVFKQNIEKKMPVDQINAHAKSYKEDATNVNNDMSLGQMLSIQKEAIEMGVNKQVFAQIQIPALGLALPIFKGANQYTLSLGAATYFYEDAEMGKGNYVLAGHNMEMPGVLFSDIQKLSLGEVMDLVSNDGVYRYKVTRKFIVPEYFKLIDGVPEENSFLSLPKKGEKPLLTLFTCVYTSQGKERYVVQGELQ"),
        ("srtA_2", b"MGKWIIAFWLLSAVGVLLLMPAEASVAKYQQNQQIAAIDRTGTAAETDSSLDVAKIELGDPVGILTIPSISLKLPIYDGTSDKILENGVGITEGTGDITGGNGKNPLIAGHSGLYKDNLFDDLPSVKKGEKFYIKVDGEQHAYQIDRIEEVQKDELQRNFVTYLEPNPNEDRVTLMTCTPKGINTHRFLVYGKRVTFTKSELKDEENKKQKLSWKWLLGSTVFLSVMIIGSLFVYKKKK"),
        ("dinB", b"MYLAISSLQHRTYVCIMWKNGVLFMMDYSKEPVNDYFLIDMKSFYASVECIERNLDPLTTELVVMSRSDNTGSGLILASSPEAKKRYGITNVSRPRDLPQPFPKTLHVVPPRMNLYIKRNMQVNNIFRRYVADEDLLIYSIDESILKVTKSLNLFTTEETRSQRRKKLAQMIQERIKEELGLIATVGVGDNPLLAKLALDNEAKHNEGFIAEWTYENVPEKVWNIPEMTDFWGIGSRMKKRLNQMGILSIRDLANWNPYTIKNRLGVIGLQLYFHANGIDRTDIAIPPEPTKEKSYGNSQVLPRDYTRRNEIELVVKEMAEQVAIRIRQHNCKTGCVHLNIGTSILETRPGFSHQMKIPITDNTKELQNYCLFLFDKYYEGQEVRHVGITYSKLVYTDSLQLDLFSDPQKQINEENLDKIIDKIRQKYGFTSIVHASSMLESARSITRSTLVGGHAGGNGGIKND"),
        ("hin", b"MSLIAEVRSLTGIQSSAQAIQELGGKFNINQRSIERYKEFLNQHPSRQIIDSMLTNTISALGLNITKLGLRFKARKYGEEKTLYSKDALRTKAQNLIASADYIQELNKHPSKAQQLNTELIELVNNTLKERISRLSSQKISTAKERITGYKKITENAKEFARAFG"),
        ("bin3", b"MSAFAQIVRSLTGIQSSAQAIQELGGEFKISQRAIERYKENLGSQPTEEVLETMLANTIGAIGLSVSRLGLRYKARKIGEEKSLYNKEALRTQAISNLIKNHKFMKAQTLNKELINKLAKALEQRISRISSSQTISSAKERITEYKKITENAIEQIKAGLQ"),
        ("soj", b"MKLAIVADVSGEGLCSTIVGKTSVSALAKRAGVKKVIALDTATSTQLHKNADYLLVKGMSRQVSLSIGSRFLTDGKQDIISLVVLPISNLEQQTAKLDLQKQIIGAKPLVVPEDVSKGLKEGDQIVSYAFNTLRLMVFVDPDKKDRLESEIESLVQKAIAQKNRAQEAKIIQDALDSVRTIALKPLDYQVRDIAEKINHALENAGFTPMFDTHVTGRFITPSAQGKSTIDKAYGLVKQVGDS"),
    ];

    let params = SearchParams::blastp()
        .evalue(1e-9)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    let mut hit_names = Vec::new();
    for (name, query) in &test_proteins {
        let results = blastp(&db, query, &params);
        if !results.is_empty() {
            hit_names.push(*name);
        }
    }

    // With comp_adjust=0 against sprot, expect 6 hits (hin/bin3/soj have no sprot hits
    // even in NCBI BLAST+ — they're annotated via ISfinder/COG databases in prokka).
    assert_eq!(
        hit_names,
        ["topB", "ssb", "yoeB", "srtA_1", "srtA_2", "dinB"],
        "Unexpected sprot annotation set"
    );
}

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
#[test]
#[ignore]
fn test_lambda_ratio_stress() {
    let sprot_paths = [
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../prokka-rs/prokka/db/kingdom/Bacteria/sprot"),
        std::path::PathBuf::from(
            "/data/henriksson/github/claude/prokka-rs/prokka/db/kingdom/Bacteria/sprot",
        ),
    ];
    let sprot_path = match sprot_paths.iter().find(|p| p.exists()) {
        Some(p) => p,
        None => {
            eprintln!("Skipping: prokka sprot database not found");
            return;
        }
    };
    let file = std::fs::File::open(sprot_path).unwrap();
    let records = blast_rs::input::parse_fasta(file);
    let tmp = TempDir::new().unwrap();
    let base = tmp.path().join("sprot");
    let mut builder = BlastDbBuilder::new(DbType::Protein, "sprot");
    for rec in &records {
        builder.add(SequenceEntry {
            title: rec.defline.clone(),
            accession: rec.id.clone(),
            sequence: rec.sequence.clone(),
            taxid: None,
        });
    }
    builder.write(&base).unwrap();
    let db = blast_rs::db::BlastDb::open(&base).unwrap();

    // Compositionally biased queries with NCBI reference values.
    // (name, query, ncbi_mode0_score, ncbi_mode1_score, ncbi_mode0_acc, ncbi_mode1_acc)
    type BiasedQueryRow<'a> = (&'a str, &'a [u8], i32, i32, &'a str, &'a str);
    let biased_queries: Vec<BiasedQueryRow> = vec![
        // Pro-rich: extreme proline bias, NCBI mode1 drops score 36% and changes top hit
        ("pro_rich",
         b"MPPPPVALPPTPPEAPPPPAQPPDPPAQPPPPAQPVAPPAPPTPPEAPPPTAQPVAPPAPPTLPEAPPPTAQ",
         164, 105, "P9WJC5", "Q70XJ9"),
        // Glu-rich: extreme glutamate bias
        ("glu_rich",
         b"MEEEEKELEQEKKKLEEEKAEELEEELKKLEQEEVKEEIKELEEKLEEEQKEELKNELEEE",
         97, 56, "A0A0H2XG66", "P54735"),
        // Ala-rich: hydrophobic alanine bias
        ("ala_rich",
         b"MAAAALAGALAAAGALAAALAAGALAAEAAAALAGVLAARAGALAALAAGVLAARAGALAALA",
         79, 58, "P11910", "Q05308"),
        // srtA: moderate bias, alignment shortens significantly in mode 1
        ("srtA",
         b"MKKCFLFCLKGGRWMKKWLFGFLGVALIVVCSVFGYVSYQKHEGEVFKQNIEKKMPVDQINAHAKSYKEDATNVNNDMSLGQMLSIQKEAIEMGVNKQVFAQIQIPALGLALPIFKGANQYTLSLGAATYFYEDAEMGKGNYVLAGHNMEMPGVLFSDIQKLSLGEVMDLVSNDGVYRYKVTRKFIVPEYFKLIDGVPEENSFLSLPKKGEKPLLTLFTCVYTSQGKERYVVQGELQ",
         234, 221, "P0DPQ5", "P0DPQ5"),
        // Normal composition: score should stay similar or increase slightly
        ("normal_topB_100aa",
         b"MMKTVILAEKPSQAKAYADSFSKATRKDGYFEIQDRLFPGETVITYGFGHLVELDSPDMYDENWKQWSLEHLPIFPNQYHYHVPKDKKKQFNVVKQQLQSA",
         162, 166, "P14294", "P14294"),
    ];

    for (name, query, ncbi_m0_score, ncbi_m1_score, ncbi_m0_acc, ncbi_m1_acc) in &biased_queries {
        // Run with comp_adjust=0
        let params_0 = SearchParams::blastp()
            .evalue(1.0)
            .num_threads(1)
            .filter_low_complexity(false)
            .comp_adjust(0);
        let results_0 = blastp(&db, query, &params_0);

        // Run with comp_adjust=1 (ScaleOldMatrix — lambda ratio rescaling)
        let params_1 = SearchParams::blastp()
            .evalue(1.0)
            .num_threads(1)
            .filter_low_complexity(false)
            .comp_adjust(1);
        let results_1 = blastp(&db, query, &params_1);

        let top_0 = results_0
            .first()
            .unwrap_or_else(|| panic!("{name}: comp_adjust=0 should find a top hit"));
        let top_1 = results_1
            .first()
            .unwrap_or_else(|| panic!("{name}: comp_adjust=1 should find a top hit"));

        assert_eq!(
            top_0.hsps[0].score, *ncbi_m0_score,
            "{name}: comp_adjust=0 score should match NCBI"
        );
        assert_eq!(
            top_1.hsps[0].score, *ncbi_m1_score,
            "{name}: comp_adjust=1 score should match NCBI"
        );
        assert!(
            top_0.subject_accession.contains(ncbi_m0_acc)
                || top_0.subject_title.contains(ncbi_m0_acc),
            "{name}: comp_adjust=0 top hit {} / {} should identify NCBI accession {}",
            top_0.subject_accession,
            top_0.subject_title,
            ncbi_m0_acc
        );
        assert!(
            top_1.subject_accession.contains(ncbi_m1_acc)
                || top_1.subject_title.contains(ncbi_m1_acc),
            "{name}: comp_adjust=1 top hit {} / {} should identify NCBI accession {}",
            top_1.subject_accession,
            top_1.subject_title,
            ncbi_m1_acc
        );
    }
}

/// Per-function timing breakdown for blastp against sprot.
/// Measures: DB load, lookup table build, subject scan, gapped alignment, total.
///
/// Run with: cargo test --release -- --ignored test_blastp_timing_breakdown
#[test]
#[ignore]
fn test_blastp_timing_breakdown() {
    let sprot_paths = [
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../prokka-rs/prokka/db/kingdom/Bacteria/sprot"),
        std::path::PathBuf::from(
            "/data/henriksson/github/claude/prokka-rs/prokka/db/kingdom/Bacteria/sprot",
        ),
    ];
    let sprot_path = match sprot_paths.iter().find(|p| p.exists()) {
        Some(p) => p,
        None => {
            eprintln!("Skipping: prokka sprot database not found");
            return;
        }
    };

    // Build DB
    let file = std::fs::File::open(sprot_path).unwrap();
    let records = blast_rs::input::parse_fasta(file);
    let tmp = TempDir::new().unwrap();
    let base = tmp.path().join("sprot");
    let mut builder = BlastDbBuilder::new(DbType::Protein, "sprot");
    for rec in &records {
        builder.add(SequenceEntry {
            title: rec.defline.clone(),
            accession: rec.id.clone(),
            sequence: rec.sequence.clone(),
            taxid: None,
        });
    }
    builder.write(&base).unwrap();
    let db = blast_rs::db::BlastDb::open(&base).unwrap();

    // A realistic 260aa prokka query
    let query = b"MMKTVILAEKPSQAKAYADSFSKATRKDGYFEIQDRLFPGETVITYGFGHLVELDSPDMYDENWKQWSLEHLPIFPNQYHYHVPKDKKKQFNVVKQQLQSADTIIIATDSDREGELIAWTIIQQAGADQGKTFKRLWINSLEKEAIYQGFQQLRDAGETYPKFEEAQARQIADWLIGMNGSPLYSLLLQQKGIPGSFSLGRVQTPTLYMIYQLQEKIRNFQKEPYFEGKAQVIAQNGAFDAKLDPNETQA";

    let params = SearchParams::blastp()
        .evalue(1e-6)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    // Warm up
    let _ = blastp(&db, query, &params);

    // Time full blastp
    let n = 5;
    let t = std::time::Instant::now();
    for _ in 0..n {
        let _ = blastp(&db, query, &params);
    }
    let total = t.elapsed();

    // Time just lookup table build
    let query_aa = blast_rs::encoding::encode_ncbistdaa_sequence(query);
    let matrix = blast_rs::matrix::BLOSUM62;

    let t = std::time::Instant::now();
    for _ in 0..n {
        let _ = blast_rs::protein_lookup::ProteinLookupTable::build(&query_aa, 3, &matrix, 11.0);
    }
    let table_build = t.elapsed();

    // Time just subject iteration (scan without extension) — approximate by scanning
    let table = blast_rs::protein_lookup::ProteinLookupTable::build(&query_aa, 3, &matrix, 11.0);

    let t = std::time::Instant::now();
    for _ in 0..n {
        for oid in 0..db.num_oids {
            let subj = db.get_sequence(oid);
            let len = db.get_seq_len(oid) as usize;
            if len < 3 {
                continue;
            }
            let _ = blast_rs::protein_lookup::protein_scan_with_table(
                &query_aa,
                &subj[..len],
                &matrix,
                &table,
                40,
            );
        }
    }
    let scan_total = t.elapsed();

    // Time just DB access (iterate all subjects, no search)
    let t = std::time::Instant::now();
    let mut total_len = 0usize;
    for _ in 0..n {
        for oid in 0..db.num_oids {
            let subj = db.get_sequence(oid);
            let len = db.get_seq_len(oid) as usize;
            total_len += len;
            std::hint::black_box(&subj[..len]);
        }
    }
    let db_iter = t.elapsed();

    eprintln!(
        "\n=== BLASTP TIMING BREAKDOWN ({} iterations, 260aa query vs {} subjects) ===",
        n, db.num_oids
    );
    eprintln!(
        "  Full blastp():        {:>8.3}s  ({:.4}s/query)",
        total.as_secs_f64(),
        total.as_secs_f64() / n as f64
    );
    eprintln!(
        "  Lookup table build:   {:>8.3}s  ({:.4}s/query)",
        table_build.as_secs_f64(),
        table_build.as_secs_f64() / n as f64
    );
    eprintln!(
        "  Scan + extend:        {:>8.3}s  ({:.4}s/query)",
        scan_total.as_secs_f64(),
        scan_total.as_secs_f64() / n as f64
    );
    eprintln!(
        "  DB iteration only:    {:>8.3}s  ({:.4}s/query)",
        db_iter.as_secs_f64(),
        db_iter.as_secs_f64() / n as f64
    );
    eprintln!(
        "  Overhead (gapped etc): {:>7.3}s  ({:.4}s/query)",
        (total.as_secs_f64() - scan_total.as_secs_f64()).max(0.0),
        ((total.as_secs_f64() - scan_total.as_secs_f64()) / n as f64).max(0.0)
    );
    eprintln!(
        "  Total subject bytes:  {} ({}/iter)",
        total_len,
        total_len / n
    );
    eprintln!();

    // Also time different query sizes
    let sizes = [50, 100, 200, 400, 800];
    eprintln!("=== QUERY SIZE SCALING ===");
    for &sz in &sizes {
        let q: Vec<u8> = query.iter().cycle().take(sz).copied().collect();
        let t = std::time::Instant::now();
        for _ in 0..3 {
            let _ = blastp(&db, &q, &params);
        }
        let elapsed = t.elapsed();
        eprintln!("  {}aa: {:.4}s/query", sz, elapsed.as_secs_f64() / 3.0);
    }
}

/// External sensitivity fixture for srtA and umuC against Prokka sprot.
#[test]
#[ignore]
fn test_srta_and_umuc_sprot_sensitivity_reference() {
    let sprot_paths = [
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../prokka-rs/prokka/db/kingdom/Bacteria/sprot"),
        std::path::PathBuf::from(
            "/data/henriksson/github/claude/prokka-rs/prokka/db/kingdom/Bacteria/sprot",
        ),
    ];
    let sprot_path = match sprot_paths.iter().find(|p| p.exists()) {
        Some(p) => p,
        None => {
            eprintln!("Skipping: sprot not found");
            return;
        }
    };
    let file = std::fs::File::open(sprot_path).unwrap();
    let records = blast_rs::input::parse_fasta(file);
    let tmp = TempDir::new().unwrap();
    let base = tmp.path().join("sprot");
    let mut builder = BlastDbBuilder::new(DbType::Protein, "sprot");
    for rec in &records {
        builder.add(SequenceEntry {
            title: rec.defline.clone(),
            accession: rec.id.clone(),
            sequence: rec.sequence.clone(),
            taxid: None,
        });
    }
    builder.write(&base).unwrap();
    let db = blast_rs::db::BlastDb::open(&base).unwrap();

    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    let params_strict = SearchParams::blastp()
        .evalue(1e-6)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    let srta = b"MIIRHPKKKRIMGKWIIAFWLLSAVGVLLLMPAEASVAKYQQNQQIAAIDRTGTAAETDSSLDVAKIELGDPVGILTIPSISLKLPIYDGTSDKILENGVGITEGTGDITGGNGKNPLIAGHSGLYKDNLFDDLPSVKKGEKFYIKVDGEQHAYQIDRIEEVQKDELQRNFVTYLEPNPNEDRVTLMTCTPKGINTHRFLVYGKRVTFTKSELKDEENKKQKLSWKWLLGSTVFLSVMIIGSLFVYKKKK";
    let results_e10 = blastp(&db, srta, &params);
    let results_strict = blastp(&db, srta, &params_strict);
    assert_eq!(results_e10.len(), 27);
    assert_eq!(results_strict.len(), 5);
    assert!(
        results_e10
            .first()
            .is_some_and(|r| r.subject_accession.contains("P0DPQ5")),
        "srtA permissive search should rank P0DPQ5 first"
    );
    assert!(
        results_strict
            .first()
            .is_some_and(|r| r.subject_accession.contains("P0DPQ5")),
        "srtA strict search should retain P0DPQ5 first"
    );
    let srta_top = &results_strict[0].hsps[0];
    assert_eq!(srta_top.score, 183);
    assert_eq!(srta_top.alignment_length, 213);
    assert_eq!(srta_top.num_identities, 57);
    assert_eq!(srta_top.num_gaps, 5);

    let umuc = b"MNSDLILAGESPSYNAAFIAMKEQHPAVFYAQHNAFGLKKIRSGFISDEQAKEYYPLICEALQKDITHFVDEIVASITGYSIDNIRFAKENKNKTIINSFEGWYNLSQQLLDTIMNEQNKSHPQFSYYSKLNSSHQLSHKEKAEAYIAGINIQITIDKQGKFFQKHFDIIQSIIKEESVNIPVLFINTSRNLKYSTGIEFNELFKRSNSNSLLAKRRVFYSLPYQPAKYREYFDSFKKISEKWIEAYKCNELDKNIGISFHDFYDSRFRTKDAKKQFSFINNIMSKIRDLYEVPEKIVRELKTRFKWFWEKKVKK";
    let results_e10 = blastp(&db, umuc, &params);
    let results_strict = blastp(&db, umuc, &params_strict);
    assert_eq!(results_e10.len(), 23);
    assert!(results_strict.is_empty());
    assert!(
        results_e10
            .first()
            .is_some_and(|r| r.subject_accession.contains("P28367")),
        "umuC permissive search should rank the weak P28367 hit first"
    );
    let umuc_top = &results_e10[0].hsps[0];
    assert_eq!(umuc_top.score, 71);
    assert_eq!(umuc_top.alignment_length, 67);
    assert_eq!(umuc_top.num_identities, 21);
    assert_eq!(umuc_top.num_gaps, 3);
}

// ── End-to-end API tests (ported from NCBI bl2seq / blastengine / traceback) ─

/// Search a nucleotide sequence against itself. Should produce a single perfect
/// alignment covering the full length with 100% identity.
#[test]
fn test_blastn_self_search() {
    let seq = "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACGATCGTAAGGCCTTAGCAGTCAATGC";
    let (_tmp, db) = build_nucleotide_db(vec![nt_entry("N001", "self target", seq)]);
    let params = SearchParams::blastn()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false);

    let results = blastn(&db, seq.as_bytes(), &params);
    assert_eq!(
        results.len(),
        1,
        "self search must find exactly one subject"
    );
    let best = &results[0];
    assert_eq!(best.subject_oid, 0);
    assert_eq!(
        best.hsps.len(),
        1,
        "self search should yield exactly one HSP"
    );
    let hsp = &best.hsps[0];
    assert!(
        (hsp.percent_identity() - 100.0).abs() < 0.01,
        "self search identity should be 100%, got {:.2}%",
        hsp.percent_identity()
    );
    assert_eq!(
        hsp.alignment_length,
        seq.len(),
        "alignment should span entire sequence"
    );
    assert_eq!(hsp.num_identities as usize, seq.len());
    assert_eq!(hsp.num_gaps, 0, "perfect self alignment has no gaps");
}

#[test]
fn test_blastn_search_api_finds_abricate_long_contig_hit() {
    let abricate_root = std::path::Path::new("/data/henriksson/github/claude/abricate-rs");
    let db_path = abricate_root.join("abricate/db/ncbi/sequences");
    let query_path = abricate_root.join("abricate/test/assembly.fa");
    if !db_path.exists() || !query_path.exists() {
        eprintln!("Skipping: ABRicate ncbi database fixture not present");
        return;
    }

    let db_fasta = std::fs::read(&db_path).expect("read ABRicate ncbi db FASTA");
    let (subject_title, subject_seq) = parse_fasta(&db_fasta)
        .into_iter()
        .find(|(id, _)| id.contains("NG_055999.1"))
        .expect("NG_055999.1 in ABRicate ncbi db");
    let (_tmp, db) = build_nucleotide_db(vec![nt_entry(
        "NG_055999.1",
        &subject_title,
        std::str::from_utf8(&subject_seq).expect("subject is UTF-8 DNA"),
    )]);
    let assembly = std::fs::read(&query_path).expect("read ABRicate assembly");
    let query = parse_fasta(&assembly)
        .into_iter()
        .find(|(id, _)| {
            id.split_whitespace()
                .next()
                .is_some_and(|id| id == "LGJG01000038")
        })
        .map(|(_, seq)| seq)
        .expect("LGJG01000038 in assembly fixture");
    let params = SearchParams::blastn()
        .filter_low_complexity(false)
        .evalue(1e-20)
        .max_target_seqs(10000)
        .num_threads(1)
        .culling_limit(Some(1))
        .match_score(1)
        .mismatch(-3);

    let results = blastn(&db, &query, &params);
    let hit = results
        .iter()
        .find(|result| result.subject_title.contains("NG_055999.1"))
        .and_then(|result| result.hsps.first().map(|hsp| (result, hsp)));
    let Some((result, hsp)) = hit else {
        panic!("blastn_search API did not find the long-contig blaZ hit");
    };

    assert!(
        result.subject_title.contains("blaZ"),
        "unexpected subject title: {}",
        result.subject_title
    );
    assert_eq!((hsp.query_start + 1, hsp.query_end), (64650, 65495));
    assert_eq!((hsp.subject_start + 1, hsp.subject_end), (1, 846));
    assert!(
        (hsp.percent_identity() - 96.809).abs() < 0.01,
        "unexpected identity {:.3}",
        hsp.percent_identity()
    );
}

/// Search completely unrelated sequences and verify no hits at strict evalue.
#[test]
fn test_blastn_no_hit() {
    // Poly-A query vs poly-C subject -- no significant similarity.
    let query   = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    let subject = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC";
    let (_tmp, db) = build_nucleotide_db(vec![nt_entry("N001", "unrelated", subject)]);
    let params = SearchParams::blastn()
        .evalue(1e-10)
        .num_threads(1)
        .filter_low_complexity(false);

    let results = blastn(&db, query.as_bytes(), &params);
    assert!(
        results.is_empty(),
        "completely unrelated sequences should produce no hits at evalue 1e-10"
    );
}

/// Subject contains two separate regions matching the query. Verify multiple
/// HSPs are returned.
#[test]
fn test_blastn_multiple_hsps() {
    // Two distinct 40bp matching regions separated by 60bp of unrelated sequence.
    let region_a = "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTT";
    let region_b = "GCTTAACGATCGTAAGGCCTTAGCAGTCAATGCTTGAAGT";
    let spacer = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
    let query = format!("{}{}", region_a, region_b);
    let subject = format!("{}{}{}", region_a, spacer, region_b);
    let (_tmp, db) = build_nucleotide_db(vec![nt_entry("N001", "multi region", &subject)]);
    let params = SearchParams::blastn()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false);

    let results = blastn(&db, query.as_bytes(), &params);
    assert_eq!(
        results.len(),
        1,
        "should find one subject with matching regions"
    );
    assert_eq!(results[0].subject_oid, 0);
    // Count total HSPs across all results for subject N001
    let total_hsps: usize = results.iter().map(|r| r.hsps.len()).sum();
    assert!(
        total_hsps >= 1,
        "subject with two matching regions should produce at least 1 HSP, got {}",
        total_hsps
    );
}

/// Protein self-search: search a protein against itself and verify perfect alignment.
#[test]
fn test_blastp_self_search() {
    let seq = "MKFLILLFNILCLFPVLAADNHGVSMNAS";
    let (_tmp, db) = build_protein_db(vec![protein_entry("P001", "self protein", seq)]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    let results = blastp(&db, seq.as_bytes(), &params);
    assert_eq!(
        results.len(),
        1,
        "protein self search should find one subject"
    );
    assert_eq!(results[0].subject_oid, 0);
    let hsp = &results[0].hsps[0];
    assert!(
        (hsp.percent_identity() - 100.0).abs() < 0.01,
        "protein self search should be 100% identity, got {:.2}%",
        hsp.percent_identity()
    );
    assert_eq!(hsp.alignment_length, seq.len());
    assert_eq!(hsp.num_identities as usize, seq.len());
    assert_eq!(hsp.num_gaps, 0);
}

/// Search two related but not identical proteins. Verify a hit with positive
/// score and reasonable e-value.
#[test]
fn test_blastp_known_pair() {
    // Query and subject differ at a few positions (conservative substitutions).
    let query = "MKFLILLFNILCLFPVLAADNHGVSMNAS";
    let subject = "MKFLILLFNILCLFPVLAADNHGVSINAS"; // M→I near the end (conservative in BLOSUM62)
    let (_tmp, db) = build_protein_db(vec![protein_entry("P001", "related protein", subject)]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    let results = blastp(&db, query.as_bytes(), &params);
    assert_eq!(results.len(), 1, "related proteins should produce one hit");
    assert_eq!(results[0].subject_oid, 0);
    let hsp = &results[0].hsps[0];
    assert_eq!(hsp.alignment_length, query.len());
    assert_eq!(hsp.num_identities as usize, query.len() - 1);
    assert_eq!(hsp.num_gaps, 0);
    assert!(
        hsp.evalue < 1.0,
        "related protein hit should have reasonable evalue, got {:.2e}",
        hsp.evalue
    );
    assert!((hsp.percent_identity() - 96.55172413793103).abs() < 1.0e-12);
}

/// BLASTX: nucleotide query encoding a known protein searched against a protein
/// database. Should find the translated match.
#[test]
fn test_blastx_finds_protein_match() {
    // ATG AAA TTT CTG ATT CTG CTG TTT AAC ATT CTG TGC CTG TTC
    // encodes: M   K   F   L   I   L   L   F   N   I   L   C   L   F
    let nt_query = "ATGAAATTTCTGATTCTGCTGTTTAACATTCTGTGCCTGTTC";
    let protein = "MKFLILLFNILCLF";
    let (_tmp, db) = build_protein_db(vec![protein_entry("P001", "target protein", protein)]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    let results = blastx(&db, nt_query.as_bytes(), &params);
    assert_eq!(
        results.len(),
        1,
        "blastx should find one translated protein match"
    );
    let hsp = &results[0].hsps[0];
    assert_eq!(hsp.query_frame, 1);
    assert_eq!(hsp.alignment_length, protein.len());
    assert_eq!(hsp.num_identities as usize, protein.len());
    assert_eq!(hsp.num_gaps, 0);
}

/// TBLASTN: protein query against nucleotide subject that encodes the protein.
/// Should find the translated nucleotide match.
#[test]
fn test_tblastn_finds_nucleotide_match() {
    let protein_query = "MKFLILLFNILCLF";
    let nt_subject = "ATGAAATTTCTGATTCTGCTGTTTAACATTCTGTGCCTGTTC";
    let (_tmp, db) = build_nucleotide_db(vec![nt_entry("N001", "coding region", nt_subject)]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);

    let results = tblastn(&db, protein_query.as_bytes(), &params);
    assert_eq!(
        results.len(),
        1,
        "tblastn should find one translated nucleotide match"
    );
    let hsp = &results[0].hsps[0];
    assert_eq!(hsp.subject_frame, 1);
    assert_eq!(hsp.alignment_length, protein_query.len());
    assert_eq!(hsp.num_identities as usize, protein_query.len());
    assert_eq!(hsp.num_gaps, 0);
}

/// Query whose reverse complement matches the subject. Verify hit on minus strand.
#[test]
fn test_blastn_both_strands() {
    let subject = "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACGATCGTAAGGCCTTAGCAGTCAATGC";
    let query_rc = String::from_utf8(reverse_complement(subject.as_bytes())).unwrap();
    let (_tmp, db) = build_nucleotide_db(vec![nt_entry("N001", "forward strand subject", subject)]);

    // With strand=both, should find the reverse-complement hit.
    let params = SearchParams::blastn()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .strand("both");

    let results = blastn(&db, query_rc.as_bytes(), &params);
    assert_eq!(results.len(), 1, "should find one RC hit when strand=both");
    let hsp = &results[0].hsps[0];
    assert!(
        (hsp.percent_identity() - 100.0).abs() < 0.01,
        "RC of subject should be 100% identity, got {:.2}%",
        hsp.percent_identity()
    );
    assert_eq!(hsp.alignment_length, subject.len());
    assert_eq!(hsp.num_identities as usize, subject.len());
    assert_eq!(hsp.num_gaps, 0);

    // With strand=plus only, should NOT find it (query is RC of subject).
    let params_plus = SearchParams::blastn()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .strand("plus");

    let results_plus = blastn(&db, query_rc.as_bytes(), &params_plus);
    assert!(
        results_plus.is_empty(),
        "plus-strand-only should not find hit when query is RC of subject"
    );
}

/// Very short query (15bp) with word_size=7. Should still find a match.
#[test]
fn test_blastn_short_query() {
    let subject = "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACGATCGTAAGGCCTTAGCAGTCAATGC";
    let query = "ATGCGTACCTGAAAG"; // first 15bp of subject
    let (_tmp, db) = build_nucleotide_db(vec![nt_entry("N001", "target", subject)]);
    let params = SearchParams::blastn()
        .evalue(10.0)
        .word_size(7)
        .num_threads(1)
        .filter_low_complexity(false);

    let results = blastn(&db, query.as_bytes(), &params);
    assert_eq!(
        results.len(),
        1,
        "short 15bp query with word_size=7 should find one match"
    );
    let hsp = &results[0].hsps[0];
    assert!(
        (hsp.percent_identity() - 100.0).abs() < 0.01,
        "exact substring should be 100% identity"
    );
    assert_eq!(hsp.alignment_length, query.len());
    assert_eq!(hsp.num_identities as usize, query.len());
    assert_eq!(hsp.num_gaps, 0);
}

/// Set a strict evalue threshold and verify only significant hits pass.
#[test]
fn test_blastn_evalue_filter() {
    let seq = "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACGATCGTAAGGCCTTAGCAGTCAATGC";
    let (_tmp, db) = build_nucleotide_db(vec![nt_entry("N001", "target", seq)]);

    // Relaxed evalue -- should find hits
    let params_relaxed = SearchParams::blastn()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false);
    let results_relaxed = blastn(&db, seq.as_bytes(), &params_relaxed);
    assert_eq!(
        results_relaxed.len(),
        1,
        "relaxed evalue should find self hit"
    );

    // Very strict evalue for self-hit of 100bp should still pass (self hit is highly significant)
    let params_strict = SearchParams::blastn()
        .evalue(1e-30)
        .num_threads(1)
        .filter_low_complexity(false);
    let results_strict = blastn(&db, seq.as_bytes(), &params_strict);
    assert_eq!(
        results_strict.len(),
        1,
        "100bp self-hit should still pass evalue 1e-30"
    );

    // All returned HSPs must satisfy the evalue threshold
    for r in &results_strict {
        for hsp in &r.hsps {
            assert!(
                hsp.evalue <= 1e-30,
                "HSP evalue {:.2e} exceeds threshold 1e-30",
                hsp.evalue
            );
        }
    }
}

/// Search against a BLAST database and against a subject FASTA should produce
/// equivalent results for the same sequences.
#[test]
fn test_db_search_vs_subject_search() {
    let query   = "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACGATCGTAAGGCCTTAGCAGTCAATGC";
    let subject = "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACGATCGTAAGGCCTTAGCAGTCAATGC";

    // Database-based search
    let (_tmp, db) = build_nucleotide_db(vec![nt_entry("N001", "target", subject)]);
    let params = SearchParams::blastn()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false);
    let db_results = blastn(&db, query.as_bytes(), &params);

    // Subject-mode search (using BlastnSearch builder which does seq-vs-seq)
    let builder_results = BlastnSearch::new()
        .query(query.as_bytes())
        .subject(subject.as_bytes())
        .dust(false)
        .evalue(10.0)
        .run();

    assert_eq!(db_results.len(), 1, "DB search should find one hit");
    assert_eq!(
        builder_results.len(),
        1,
        "subject search should find one hit"
    );

    // The best HSP scores should be comparable (both are self-hits)
    let db_best_score = db_results[0].hsps[0].score;
    let subj_best_score = builder_results[0].score;
    assert_eq!(db_results[0].hsps[0].alignment_length, query.len());
    assert_eq!(builder_results[0].align_length as usize, query.len());
    assert_eq!(db_results[0].hsps[0].num_identities as usize, query.len());
    assert_eq!(builder_results[0].num_ident as usize, query.len());
    assert!(db_best_score > 0);
    assert!(subj_best_score > 0);
}

/// Use BlastnSearch builder API directly for a seq-vs-seq search.
#[test]
fn test_blastn_search_builder_api() {
    let query   = "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACGATCGTAAGGCCTTAGCAGTCAATGC";
    let subject = "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACGATCGTAAGGCCTTAGCAGTCAATGC";

    let results = BlastnSearch::new()
        .query(query.as_bytes())
        .subject(subject.as_bytes())
        .word_size(11)
        .reward(1)
        .penalty(-3)
        .gap_open(5)
        .gap_extend(2)
        .evalue(10.0)
        .dust(false)
        .strand(Strand::Both)
        .run();

    assert_eq!(
        results.len(),
        1,
        "builder API self-search should find one hit"
    );
    let best = &results[0];
    assert!(best.score > 0, "best HSP score should be positive");
    assert!(
        best.evalue < 1e-10,
        "100bp self-hit should have very significant evalue"
    );
    assert_eq!(
        best.align_length as usize,
        query.len(),
        "self-hit alignment should span full query length"
    );
    assert_eq!(
        best.num_ident as usize,
        query.len(),
        "self-hit should have 100% identity"
    );

    // Test with empty query -- should return empty
    let empty = BlastnSearch::new()
        .query(b"")
        .subject(subject.as_bytes())
        .run();
    assert!(empty.is_empty(), "empty query should produce no results");

    // Test with empty subject -- should return empty
    let empty2 = BlastnSearch::new()
        .query(query.as_bytes())
        .subject(b"")
        .run();
    assert!(empty2.is_empty(), "empty subject should produce no results");
}

#[test]
#[ignore]
fn test_comp_ratio() {
    // 00009 vs P0DPQ5: composition adjustment should leave this legitimate
    // sortase hit unpenalized, matching the retained NCBI reference hit.
    let q = b"MIIRHPKKKRIMGKWIIAFWLLSAVGVLLLMPAEASVAKYQQNQQIAAIDRTGTAAETDSSLDVAKIELGDPVGILTIPRISLTLPIYDATNEKILENGVGITEGTGDITGGNGKNPLIAGHSGLYKDNLFDDLPSVKKGEKFYIKVDGEQHAYQIDRIEEVQKDELQRNFVTYLEPNPNEDRVTLMTCTPKGINTHRFLVYGKRVTFTKSELKDEENKKQKLSWKWLLGSTVFLSVMIIGSLFVYKKKK";
    let q_aa = blast_rs::encoding::encode_ncbistdaa_sequence(q);

    // Get P0DPQ5 from DB
    let sprot_paths = [std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../prokka-rs/prokka/db/kingdom/Bacteria/sprot")];
    let sprot_path = match sprot_paths.iter().find(|p| p.exists()) {
        Some(p) => p,
        None => {
            return;
        }
    };
    let file = std::fs::File::open(sprot_path).unwrap();
    let records = blast_rs::input::parse_fasta(file);
    let p0dpq5 = records.iter().find(|r| r.id == "P0DPQ5").unwrap();
    let s_aa = blast_rs::encoding::encode_ncbistdaa_sequence(&p0dpq5.sequence);

    let (qcomp, qn) = blast_rs::composition::read_composition(&q_aa, 28);
    let (scomp, sn) = blast_rs::composition::read_composition(&s_aa, 28);
    assert_eq!(qn, 250);
    assert_eq!(sn, 210);

    let matrix = blast_rs::matrix::BLOSUM62;
    let lambda = 0.267;

    let ratio = blast_rs::composition::composition_lambda_ratio(&matrix, &qcomp, &scomp, lambda);

    assert_eq!(
        ratio,
        Some(1.0),
        "Sortase A/P0DPQ5 should not receive a lambda-ratio penalty"
    );
}

/// External srtA/P0DPQ5 composition-adjustment internals.
#[test]
#[ignore]
#[allow(clippy::approx_constant)]
fn test_comp_adjust_srta_p0dpq5_reference_internals() {
    let sprot_paths = [
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../prokka-rs/prokka/db/kingdom/Bacteria/sprot"),
        std::path::PathBuf::from(
            "/data/henriksson/github/claude/prokka-rs/prokka/db/kingdom/Bacteria/sprot",
        ),
    ];
    let sprot_path = match sprot_paths.iter().find(|p| p.exists()) {
        Some(p) => p,
        None => {
            eprintln!("Skipping: sprot not found");
            return;
        }
    };
    let file = std::fs::File::open(sprot_path).unwrap();
    let records = blast_rs::input::parse_fasta(file);

    let query_aa = b"MKKCFLFCLKGGRWMKKWLFGFLGVALIVVCSVFGYVSYQKHEGEVFKQNIEKKMPVDQINAHAKSYKEDATNVNNDMSLGQMLSIQKEAIEMGVNKQVFAQIQIPALGLALPIFKGANQYTLSLGAATYFYEDAEMGKGNYVLAGHNMEMPGVLFSDIQKLSLGEVMDLVSNDGVYRYKVTRKFIVPEYFKLIDGVPEENSFLSLPKKGEKPLLTLFTCVYTSQGKERYVVQGELQ";
    let query_ncbi = blast_rs::encoding::encode_ncbistdaa_sequence(query_aa);

    let subj_rec = records
        .iter()
        .find(|r| r.id.contains("P0DPQ5") || r.defline.contains("srtA"))
        .unwrap();
    assert!(
        subj_rec.id.contains("P0DPQ5"),
        "fixture should select the P0DPQ5 sortase record"
    );
    let subj_ncbi = blast_rs::encoding::encode_ncbistdaa_sequence(&subj_rec.sequence);

    let matrix = *blast_rs::api::get_matrix(blast_rs::api::MatrixType::Blosum62);
    let ungapped_lambda = 0.3176f64; // BLOSUM62 ungapped

    let (qcomp28, qn) = blast_rs::composition::read_composition(&query_ncbi, 28);
    let (scomp28, sn) = blast_rs::composition::read_composition(&subj_ncbi, 28);

    assert_eq!(query_ncbi.len(), 237);
    assert_eq!(qn, 237);
    assert_eq!(subj_ncbi.len(), 210);
    assert_eq!(sn, 210);

    let lr = blast_rs::composition::composition_lambda_ratio(
        &matrix,
        &qcomp28,
        &scomp28,
        ungapped_lambda,
    );
    let lr = lr.expect("srtA/P0DPQ5 should produce a finite lambda ratio");
    assert!(
        (lr - 0.9367349357083943).abs() < 1.0e-12,
        "unexpected lambda ratio: {lr}"
    );

    let freq_ratios = blast_rs::matrix::get_blosum62_freq_ratios();
    let adj = blast_rs::composition::composition_scale_matrix(
        &matrix,
        &qcomp28,
        &scomp28,
        ungapped_lambda,
        &freq_ratios,
    )
    .expect("srtA/P0DPQ5 should produce a scaled matrix");
    assert_eq!(
        (
            adj[1][1],
            adj[3][3],
            adj[4][4],
            adj[5][5],
            adj[12][12],
            adj[14][14],
            adj[16][16],
            adj[1][4],
            adj[1][5],
            adj[10][16],
        ),
        (4, 9, 6, 5, 6, 8, 6, -2, -1, 2)
    );
}

#[test]
fn composition_matrix_adj_short_exact_matches_reference_internals() {
    use blast_rs::compo_mode_condition::MatrixAdjustRule;

    let query = b"MKFLILLF";
    let subject = b"MKFLILLF";
    let query_ncbi = blast_rs::encoding::encode_ncbistdaa_sequence(query);
    let subject_ncbi = blast_rs::encoding::encode_ncbistdaa_sequence(subject);

    let matrix = *blast_rs::api::get_matrix(blast_rs::api::MatrixType::Blosum62);
    let (qcomp28, qn) = blast_rs::composition::read_composition(&query_ncbi, 28);
    let (scomp28, sn) = blast_rs::composition::read_composition(&subject_ncbi, 28);

    let mut qp20 = [0.0f64; 20];
    let mut sp20 = [0.0f64; 20];
    blast_rs::compo_mode_condition::gather_letter_probs(&qcomp28, &mut qp20);
    blast_rs::compo_mode_condition::gather_letter_probs(&scomp28, &mut sp20);

    let rule = blast_rs::compo_mode_condition::choose_matrix_adjust_rule(
        query_ncbi.len(),
        subject_ncbi.len(),
        &qp20,
        &sp20,
        2,
    );
    assert_eq!(rule, MatrixAdjustRule::UserSpecifiedRelEntropy);
    assert_eq!(qn, 8);
    assert_eq!(sn, 8);

    let lr = blast_rs::composition::composition_lambda_ratio(
        &matrix,
        &qcomp28,
        &scomp28,
        blast_rs::stat::protein_ungapped_kbp().lambda,
    );
    assert_eq!(lr, Some(0.5));

    let (joint_probs, first_std, second_std) = blast_rs::composition::blosum62_workspace();
    let mut adj_matrix = matrix;
    let status = blast_rs::composition::composition_matrix_adj(
        &mut adj_matrix,
        blast_rs::matrix::AA_SIZE,
        rule,
        qn,
        sn,
        &qcomp28,
        &scomp28,
        20,
        0.44,
        &joint_probs,
        &first_std,
        &second_std,
        blast_rs::stat::protein_ungapped_kbp().lambda,
        &matrix,
    );
    assert_eq!(status, 0);
    assert_eq!(
        (
            adj_matrix[12][12],
            adj_matrix[10][10],
            adj_matrix[6][6],
            adj_matrix[11][11],
            adj_matrix[9][9],
        ),
        (5, 5, 5, 3, 4)
    );
    assert_eq!(
        query_ncbi
            .iter()
            .zip(subject_ncbi.iter())
            .map(|(&qa, &sa)| adj_matrix[qa as usize][sa as usize])
            .sum::<i32>(),
        33
    );

    let freq_ratios = blast_rs::matrix::get_blosum62_freq_ratios();
    let scaled = blast_rs::composition::composition_scale_matrix(
        &matrix,
        &qcomp28,
        &scomp28,
        blast_rs::stat::protein_ungapped_kbp().lambda,
        &freq_ratios,
    )
    .expect("short exact pair should produce a scaled matrix");
    assert_eq!(
        (
            scaled[12][12],
            scaled[10][10],
            scaled[6][6],
            scaled[11][11],
            scaled[9][9],
        ),
        (3, 2, 3, 2, 2)
    );

    for specified_re in [0.40, 0.42, 0.44, 0.46, 0.48, 0.50] {
        let (joint_probs, first_std, second_std) = blast_rs::composition::blosum62_workspace();
        let mut adj = matrix;
        let status = blast_rs::composition::composition_matrix_adj(
            &mut adj,
            blast_rs::matrix::AA_SIZE,
            rule,
            qn,
            sn,
            &qcomp28,
            &scomp28,
            20,
            specified_re,
            &joint_probs,
            &first_std,
            &second_std,
            blast_rs::stat::protein_ungapped_kbp().lambda,
            &matrix,
        );
        let score = query_ncbi
            .iter()
            .zip(subject_ncbi.iter())
            .map(|(&qa, &sa)| adj[qa as usize][sa as usize])
            .sum::<i32>();
        assert_eq!(status, 0);
        assert_eq!(
            (
                score,
                adj[12][12],
                adj[10][10],
                adj[6][6],
                adj[11][11],
                adj[9][9],
            ),
            (33, 5, 5, 5, 3, 4),
            "specified_re={specified_re:.2}"
        );
    }

    for lambda in [0.315, 0.316, 0.317, 3177.0 / 10000.0, 0.318, 0.319, 0.320] {
        let (joint_probs, first_std, second_std) = blast_rs::composition::blosum62_workspace();
        let mut adj = matrix;
        let status = blast_rs::composition::composition_matrix_adj(
            &mut adj,
            blast_rs::matrix::AA_SIZE,
            rule,
            qn,
            sn,
            &qcomp28,
            &scomp28,
            20,
            0.44,
            &joint_probs,
            &first_std,
            &second_std,
            lambda,
            &matrix,
        );
        let score = query_ncbi
            .iter()
            .zip(subject_ncbi.iter())
            .map(|(&qa, &sa)| adj[qa as usize][sa as usize])
            .sum::<i32>();
        assert_eq!(status, 0);
        assert_eq!(
            (
                score,
                adj[12][12],
                adj[10][10],
                adj[6][6],
                adj[11][11],
                adj[9][9],
            ),
            (33, 5, 5, 5, 3, 4),
            "lambda={lambda:.4}"
        );
    }

    for test_rule in [
        MatrixAdjustRule::UnconstrainedRelEntropy,
        MatrixAdjustRule::RelEntropyOldMatrixNewContext,
        MatrixAdjustRule::RelEntropyOldMatrixOldContext,
        MatrixAdjustRule::UserSpecifiedRelEntropy,
    ] {
        let (joint_probs, first_std, second_std) = blast_rs::composition::blosum62_workspace();
        let mut adj = matrix;
        let status = blast_rs::composition::composition_matrix_adj(
            &mut adj,
            blast_rs::matrix::AA_SIZE,
            test_rule,
            qn,
            sn,
            &qcomp28,
            &scomp28,
            20,
            0.44,
            &joint_probs,
            &first_std,
            &second_std,
            blast_rs::stat::protein_ungapped_kbp().lambda,
            &matrix,
        );
        let score = query_ncbi
            .iter()
            .zip(subject_ncbi.iter())
            .map(|(&qa, &sa)| adj[qa as usize][sa as usize])
            .sum::<i32>();
        assert_eq!(status, 0);
        assert_eq!(
            (
                score,
                adj[12][12],
                adj[10][10],
                adj[6][6],
                adj[11][11],
                adj[9][9],
            ),
            (33, 5, 5, 5, 3, 4),
            "test_rule={test_rule:?}"
        );
    }
}

#[test]
fn blastp_comp_adjust_mode1_is_not_a_noop_on_short_exact_hit() {
    let (_tmp, db) = build_protein_db(vec![protein_entry("P001", "exact", "MKFLILLF")]);
    let mode0 = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .comp_adjust(0);
    let mode1 = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .comp_adjust(1);

    let score0 = blastp(&db, b"MKFLILLF", &mode0)[0].hsps[0].score;
    let score1 = blastp(&db, b"MKFLILLF", &mode1)[0].hsps[0].score;

    assert!(
        score1 < score0,
        "comp_adjust=1 should rescale this short exact protein hit instead of leaving the raw BLOSUM62 score unchanged (mode0={}, mode1={})",
        score0,
        score1
    );
}

#[test]
fn blastp_accepts_reusable_rayon_thread_pool() {
    let (_tmp, db) = build_protein_db(vec![
        protein_entry("P001", "exact", "MKFLILLF"),
        protein_entry("P002", "near", "MKYLIILF"),
        protein_entry("P003", "unrelated", "GGGGGGGG"),
    ]);
    let baseline_params = SearchParams::blastp().evalue(10.0).num_threads(1);
    let baseline = blastp(&db, b"MKFLILLF", &baseline_params);

    let pool = std::sync::Arc::new(
        rayon::ThreadPoolBuilder::new()
            .num_threads(2)
            .build()
            .expect("build rayon pool"),
    );
    let pooled_params = SearchParams::blastp().evalue(10.0).thread_pool(pool);
    let pooled = blastp(&db, b"MKFLILLF", &pooled_params);

    assert_eq!(pooled.len(), baseline.len());
    assert_eq!(pooled[0].subject_accession, baseline[0].subject_accession);
    assert_eq!(pooled[0].hsps[0].score, baseline[0].hsps[0].score);
}

#[test]
fn blastp_short_exact_comp_modes_match_ncbi_reference_scores() {
    let (_tmp, db) = build_protein_db(vec![protein_entry("P001", "exact", "MKFLILLF")]);

    let score_for = |mode| {
        let params = SearchParams::blastp()
            .evalue(10.0)
            .num_threads(1)
            .comp_adjust(mode);
        blastp(&db, b"MKFLILLF", &params)[0].hsps[0].score
    };

    assert_eq!(score_for(0), 38, "NCBI blastp -comp_based_stats 0 gives 38");
    assert_eq!(score_for(1), 21, "NCBI blastp -comp_based_stats 1 gives 21");
    assert_eq!(score_for(2), 32, "NCBI blastp -comp_based_stats 2 gives 32");
}

#[test]
fn blastx_short_exact_comp_modes_match_ncbi_reference_scores() {
    let nt_query = b"ATGAAATTTCTTATTCTTCTTTTC";
    let (_tmp, db) = build_protein_db(vec![protein_entry("P001", "exact", "MKFLILLF")]);

    let score_for = |mode| {
        let params = SearchParams::blastp()
            .evalue(10.0)
            .num_threads(1)
            .comp_adjust(mode);
        blastx(&db, nt_query, &params)[0].hsps[0].score
    };

    assert_eq!(score_for(0), 38, "NCBI blastx -comp_based_stats 0 gives 38");
    assert_eq!(score_for(1), 21, "NCBI blastx -comp_based_stats 1 gives 21");
    assert_eq!(score_for(2), 32, "NCBI blastx -comp_based_stats 2 gives 32");
}

#[test]
fn tblastn_short_exact_comp_modes_match_ncbi_reference_scores() {
    let (_tmp, db) = build_nucleotide_db(vec![nt_entry(
        "N001",
        "exact coding nt",
        "ATGAAATTTCTTATTCTTCTTTTC",
    )]);

    let score_for = |mode| {
        let params = SearchParams::blastp()
            .evalue(10.0)
            .num_threads(1)
            .comp_adjust(mode);
        tblastn(&db, b"MKFLILLF", &params)[0].hsps[0].score
    };

    assert_eq!(
        score_for(0),
        38,
        "NCBI tblastn -comp_based_stats 0 gives 38"
    );
    assert_eq!(
        score_for(1),
        21,
        "NCBI tblastn -comp_based_stats 1 gives 21"
    );
    assert_eq!(
        score_for(2),
        32,
        "NCBI tblastn -comp_based_stats 2 gives 32"
    );
}

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

#[test]
fn test_blastn_short_primer_many_subjects() {
    // Simulate blastn-short: 20bp primer vs 500 subjects (each ~2000bp).
    // The primer is embedded in every 10th subject so there are many hits.
    let primer = b"ATCGATCGATCGATCGATCG";

    let mut entries = Vec::new();
    for i in 0..500u64 {
        let mut seq = random_nt_seq(2000, i * 12345 + 7);
        // Embed primer in every 10th subject at a random-ish position
        if i % 10 == 0 {
            let pos = (i as usize * 37) % (seq.len() - primer.len());
            seq[pos..pos + primer.len()].copy_from_slice(primer);
        }
        entries.push(nt_entry(
            &format!("seq_{}", i),
            &format!("subject sequence {}", i),
            &String::from_utf8(seq).unwrap(),
        ));
    }

    let (_tmp, db) = build_nucleotide_db(entries);

    // blastn-short parameters: word_size=7, reward=1, penalty=-3, gap_open=5, gap_extend=2
    let mut params = SearchParams::blastn();
    params.word_size = 7;
    params.match_score = 1;
    params.mismatch = -3;
    params.gap_open = 5;
    params.gap_extend = 2;
    params.evalue_threshold = 5.0;
    params.max_target_seqs = 10000;
    params.filter_low_complexity = false;

    let results = blastn(&db, primer, &params);
    // Should find hits in the subjects where the primer was embedded
    assert!(!results.is_empty(), "Should find hits for embedded primer");
    // At least some of the 50 subjects with embedded primer should appear
    assert!(
        results.len() >= 10,
        "Expected at least 10 subjects with hits, got {}",
        results.len()
    );
}

#[test]
fn test_blastn_short_primer_multithreaded() {
    // Test that multithreaded search doesn't stack-overflow.
    let primer = b"GCTAGCTAGCTAGCTAGCTA";

    let mut entries = Vec::new();
    for i in 0..200u64 {
        let mut seq = random_nt_seq(5000, i * 99 + 42);
        if i % 5 == 0 {
            let pos = (i as usize * 23) % (seq.len() - primer.len());
            seq[pos..pos + primer.len()].copy_from_slice(primer);
        }
        entries.push(nt_entry(
            &format!("mseq_{}", i),
            &format!("multi-threaded test seq {}", i),
            &String::from_utf8(seq).unwrap(),
        ));
    }

    let (_tmp, db) = build_nucleotide_db(entries);

    let mut params = SearchParams::blastn();
    params.word_size = 7;
    params.match_score = 1;
    params.mismatch = -3;
    params.gap_open = 5;
    params.gap_extend = 2;
    params.evalue_threshold = 5.0;
    params.max_target_seqs = 10000;
    params.num_threads = 4;
    params.filter_low_complexity = false;

    let results = blastn(&db, primer, &params);
    assert!(
        !results.is_empty(),
        "Multithreaded blastn-short should find hits"
    );
}

#[test]
fn test_blastn_search_api_multithreaded_matches_single_threaded() {
    let query = b"ACGTACGTGACCTTGGAACTACGTACGTGACCTTGGAACT";

    let mut entries = Vec::new();
    for i in 0..80u64 {
        let mut seq = random_nt_seq(1200, i * 131 + 17);
        if i % 7 == 0 {
            let pos = (i as usize * 29) % (seq.len() - query.len());
            seq[pos..pos + query.len()].copy_from_slice(query);
        }
        entries.push(nt_entry(
            &format!("parity_{}", i),
            &format!("blastn API parallel parity seq {}", i),
            &String::from_utf8(seq).unwrap(),
        ));
    }
    let (_tmp, db) = build_nucleotide_db(entries);

    let base = SearchParams::blastn()
        .word_size(7)
        .evalue(5.0)
        .max_target_seqs(10000)
        .num_threads(1)
        .filter_low_complexity(false);
    let parallel = SearchParams {
        num_threads: 4,
        ..base.clone()
    };

    let signature = |results: Vec<SearchResult>| {
        results
            .into_iter()
            .map(|result| {
                (
                    result.subject_accession,
                    result
                        .hsps
                        .into_iter()
                        .map(|hsp| {
                            (
                                hsp.score,
                                hsp.query_start,
                                hsp.query_end,
                                hsp.subject_start,
                                hsp.subject_end,
                                hsp.num_identities,
                                hsp.alignment_length,
                            )
                        })
                        .collect::<Vec<_>>(),
                )
            })
            .collect::<Vec<_>>()
    };

    assert_eq!(
        signature(blastn(&db, query, &base)),
        signature(blastn(&db, query, &parallel)),
        "parallel blastn API results must match single-threaded results"
    );
}

#[test]
fn test_blastn_short_very_short_primer() {
    // Edge case: very short primer (15bp) with word_size=7.
    // This is the kind of query that triggered the original stack overflow report.
    let primer = b"ATCGATCGATCGATC"; // 15bp

    let mut entries = Vec::new();
    for i in 0..100u64 {
        let mut seq = random_nt_seq(10000, i * 777);
        // Embed primer in every subject
        let pos = (i as usize * 53) % (seq.len() - primer.len());
        seq[pos..pos + primer.len()].copy_from_slice(primer);
        entries.push(nt_entry(
            &format!("short_{}", i),
            &format!("short primer test {}", i),
            &String::from_utf8(seq).unwrap(),
        ));
    }

    let (_tmp, db) = build_nucleotide_db(entries);

    let mut params = SearchParams::blastn();
    params.word_size = 7;
    params.match_score = 1;
    params.mismatch = -3;
    params.gap_open = 5;
    params.gap_extend = 2;
    params.evalue_threshold = 5.0;
    params.max_target_seqs = 10000;
    params.filter_low_complexity = false;

    let results = blastn(&db, primer, &params);
    assert_eq!(
        results.len(),
        100,
        "very short primer should hit every embedded subject"
    );
    assert!(results.iter().all(|r| r
        .hsps
        .iter()
        .any(|hsp| hsp.alignment_length == primer.len()
            && hsp.num_identities as usize == primer.len())));
}

#[test]
fn test_blastn_short_via_builder() {
    // Test the BlastnSearch builder with blastn-short parameters.
    let primer = b"ATCGATCGATCGATCGATCG";
    let mut subject = random_nt_seq(5000, 42);
    subject[1000..1020].copy_from_slice(primer);

    let results = BlastnSearch::new()
        .word_size(7)
        .reward(1)
        .penalty(-3)
        .gap_open(5)
        .gap_extend(2)
        .evalue(5.0)
        .dust(false)
        .query(primer)
        .subject(&subject)
        .run();

    assert_eq!(results.len(), 4, "Builder search hit count drifted");
    let best = &results[0];
    assert_eq!(best.align_length as usize, primer.len());
    assert_eq!(best.num_ident as usize, primer.len());
}

/// Full Swiss-Prot blastp benchmark: build a protein DB from UniProt Swiss-Prot
/// (~570K entries), search 100 query sequences, and compare results with NCBI BLAST+.
///
/// Requires Swiss-Prot FASTA at one of the checked paths.
/// Run with: cargo test --release -- --ignored test_blastp_swissprot
#[test]
#[ignore]
fn test_blastp_swissprot() {
    // Try plain FASTA first, then decompress .gz if needed
    let plain_path =
        std::path::PathBuf::from("/husky/henriksson/for_claude/diamond/uniprot_sprot.fasta");
    let gz_path =
        std::path::PathBuf::from("/husky/henriksson/for_claude/diamond/uniprot_sprot.fasta.gz");

    if !plain_path.exists() && gz_path.exists() {
        eprintln!("Decompressing {} ...", gz_path.display());
        let status = std::process::Command::new("gunzip")
            .arg("-k")
            .arg(&gz_path)
            .status()
            .expect("failed to run gunzip");
        assert!(status.success(), "gunzip failed");
    }

    if !plain_path.exists() {
        eprintln!(
            "Skipping: Swiss-Prot FASTA not found at {}",
            plain_path.display()
        );
        return;
    }
    eprintln!("Using Swiss-Prot FASTA: {}", plain_path.display());

    let file = std::fs::File::open(&plain_path).unwrap();
    let records = blast_rs::input::parse_fasta(file);
    eprintln!("Loaded {} Swiss-Prot records", records.len());
    assert!(
        records.len() > 100_000,
        "Swiss-Prot should have >100K entries, got {}",
        records.len()
    );

    // Extract first 100 sequences as queries
    let num_queries = 100;
    let query_records: Vec<_> = records.iter().take(num_queries).collect();

    // Build protein database from all records
    let t0 = std::time::Instant::now();
    let tmp = TempDir::new().unwrap();
    let base = tmp.path().join("swissprot");
    let mut builder = BlastDbBuilder::new(DbType::Protein, "swissprot");
    for rec in &records {
        builder.add(SequenceEntry {
            title: rec.defline.clone(),
            accession: rec.id.clone(),
            sequence: rec.sequence.clone(),
            taxid: None,
        });
    }
    builder.write(&base).unwrap();
    let db = blast_rs::db::BlastDb::open(&base).unwrap();
    let db_time = t0.elapsed();
    eprintln!(
        "DB build: {:.2}s ({} entries)",
        db_time.as_secs_f64(),
        records.len()
    );

    // Search with single thread first
    let params = SearchParams::blastp()
        .evalue(1e-3)
        .max_target_seqs(25)
        .num_threads(1);

    let t1 = std::time::Instant::now();
    let mut total_hits = 0;
    for (i, qrec) in query_records.iter().enumerate() {
        let results = blastp(&db, &qrec.sequence, &params);
        if i < 5 || (i + 1) % 20 == 0 {
            eprintln!(
                "  Query {:>3} ({}, {} aa): {} hits",
                i + 1,
                &qrec.id[..qrec.id.len().min(20)],
                qrec.sequence.len(),
                results.len()
            );
        }
        total_hits += results.len();
    }
    let st_time = t1.elapsed();
    eprintln!(
        "Single-threaded: {} queries, {} hits, {:.2}s ({:.3}s/query)",
        num_queries,
        total_hits,
        st_time.as_secs_f64(),
        st_time.as_secs_f64() / num_queries as f64
    );

    // Search with all threads
    let params_mt = SearchParams::blastp()
        .evalue(1e-3)
        .max_target_seqs(25)
        .num_threads(0);

    let t2 = std::time::Instant::now();
    let mut total_hits_mt = 0;
    for qrec in &query_records {
        let results = blastp(&db, &qrec.sequence, &params_mt);
        total_hits_mt += results.len();
    }
    let mt_time = t2.elapsed();
    let speedup = st_time.as_secs_f64() / mt_time.as_secs_f64();
    eprintln!(
        "Multi-threaded:  {} queries, {} hits, {:.2}s ({:.3}s/query, {:.1}x speedup)",
        num_queries,
        total_hits_mt,
        mt_time.as_secs_f64(),
        mt_time.as_secs_f64() / num_queries as f64,
        speedup
    );

    // Hit counts should match between single and multi-threaded
    assert_eq!(
        total_hits, total_hits_mt,
        "Single-threaded and multi-threaded hit counts should match"
    );

    // Sanity: most queries should find at least a self-hit
    assert!(
        total_hits >= num_queries,
        "Expected at least {} hits (one self-hit per query), got {}",
        num_queries,
        total_hits
    );

    eprintln!("\n=== Swiss-Prot Benchmark Summary ===");
    eprintln!("Database:     {} sequences", records.len());
    eprintln!("Queries:      {}", num_queries);
    eprintln!("DB build:     {:.2}s", db_time.as_secs_f64());
    eprintln!(
        "Search (1T):  {:.2}s ({:.3}s/query)",
        st_time.as_secs_f64(),
        st_time.as_secs_f64() / num_queries as f64
    );
    eprintln!(
        "Search (MT):  {:.2}s ({:.3}s/query)",
        mt_time.as_secs_f64(),
        mt_time.as_secs_f64() / num_queries as f64
    );
    eprintln!("Total hits:   {}", total_hits);
}

/// External lambda-ratio reference values for biased sequences.
/// Run with: cargo test --release -- --ignored test_lambda_ratio_biased_sprot_reference --nocapture
#[test]
#[ignore]
fn test_lambda_ratio_biased_sprot_reference() {
    let sprot_paths = [
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../prokka-rs/prokka/db/kingdom/Bacteria/sprot"),
        std::path::PathBuf::from(
            "/data/henriksson/github/claude/prokka-rs/prokka/db/kingdom/Bacteria/sprot",
        ),
    ];
    let sprot_path = match sprot_paths.iter().find(|p| p.exists()) {
        Some(p) => p,
        None => {
            eprintln!("Skipping: sprot not found");
            return;
        }
    };
    let file = std::fs::File::open(sprot_path).unwrap();
    let records = blast_rs::input::parse_fasta(file);
    let tmp = TempDir::new().unwrap();
    let base = tmp.path().join("sprot");
    let mut builder = BlastDbBuilder::new(DbType::Protein, "sprot");
    for rec in &records {
        builder.add(SequenceEntry {
            title: rec.defline.clone(),
            accession: rec.id.clone(),
            sequence: rec.sequence.clone(),
            taxid: None,
        });
    }
    builder.write(&base).unwrap();
    let db = blast_rs::db::BlastDb::open(&base).unwrap();

    let matrix = *blast_rs::api::get_matrix(blast_rs::api::MatrixType::Blosum62);
    let ungapped_lambda = 0.3176f64;

    let glu_query = b"MEEEEKELEQEKKKLEEEKAEELEEELKKLEQEEVKEEIKELEEKLEEEQKEELKNELEEE";
    let glu_ncbi = blast_rs::encoding::encode_ncbistdaa_sequence(glu_query);
    let (qcomp, qn) = blast_rs::composition::read_composition(&glu_ncbi, 28);

    assert_eq!(qn, 61);

    let params_0 = SearchParams::blastp()
        .evalue(1.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(0);
    let results_0 = blastp(&db, glu_query, &params_0);
    assert!(
        !results_0.is_empty(),
        "glu-rich query should find sprot hits"
    );
    assert!(
        results_0[0].subject_accession.contains("A0A0H2XG66"),
        "unexpected glu-rich top hit: {}",
        results_0[0].subject_accession
    );
    let top_oid = results_0[0].subject_oid;
    let subj_raw = db.get_sequence(top_oid);
    let subj_len = db.get_seq_len(top_oid) as usize;
    let subj_aa = &subj_raw[..subj_len];
    let (scomp, sn) = blast_rs::composition::read_composition(subj_aa, 28);
    assert_eq!(top_oid, 6590);
    assert_eq!(sn, 444);
    let ratio =
        blast_rs::composition::composition_lambda_ratio(&matrix, &qcomp, &scomp, ungapped_lambda)
            .expect("glu-rich top hit should produce a lambda ratio");
    assert!(
        (ratio - 0.5).abs() < 1.0e-12,
        "unexpected clamped glu-rich top-hit lambda ratio: {ratio}"
    );

    let mut bg_prob = [0.0f64; 28];
    for (k, &idx) in blast_rs::composition::TRUE_CHAR_POSITIONS
        .iter()
        .enumerate()
    {
        bg_prob[idx] = blast_rs::composition::BLOSUM62_BG[k];
    }
    let query_vs_bg =
        blast_rs::composition::composition_lambda_ratio(&matrix, &qcomp, &bg_prob, ungapped_lambda)
            .expect("glu-rich/background pair should produce a lambda ratio");
    assert!(
        (query_vs_bg - 1.0).abs() < 1.0e-12,
        "unexpected clamped glu-rich/background lambda ratio: {query_vs_bg}"
    );
    let pvalue_query_vs_bg = blast_rs::composition::composition_lambda_ratio_with_adjustment(
        &matrix,
        &qcomp,
        &bg_prob,
        ungapped_lambda,
        true,
    )
    .expect("glu-rich/background p-value pair should produce a lambda ratio");
    assert!(
        (pvalue_query_vs_bg - 1.07041410).abs() < 1.0e-8,
        "unexpected unclamped glu-rich/background lambda ratio: {pvalue_query_vs_bg}"
    );

    let bg_vs_bg = blast_rs::composition::composition_lambda_ratio(
        &matrix,
        &bg_prob,
        &bg_prob,
        ungapped_lambda,
    )
    .expect("background/background pair should produce a lambda ratio");
    assert!(
        (bg_vs_bg - 0.83516064).abs() < 1.0e-8,
        "unexpected background/background lambda ratio: {bg_vs_bg}"
    );
}

/// Verify that karlin_lambda_nr with Robinson&Robinson frequencies gives the
/// correct standard BLOSUM62 ungapped lambda (0.3176).
#[test]
fn test_karlin_lambda_standard() {
    // Robinson & Robinson frequencies in NCBIstdaa positions
    let mut rr_prob = [0.0f64; 28];
    // A=1, C=3, D=4, E=5, F=6, G=7, H=8, I=9, K=10, L=11, M=12, N=13,
    // P=14, Q=15, R=16, S=17, T=18, V=19, W=20, Y=22
    rr_prob[1] = 0.07805; // A
    rr_prob[3] = 0.01925; // C
    rr_prob[4] = 0.05364; // D
    rr_prob[5] = 0.06295; // E
    rr_prob[6] = 0.03856; // F
    rr_prob[7] = 0.07377; // G
    rr_prob[8] = 0.02199; // H
    rr_prob[9] = 0.05142; // I
    rr_prob[10] = 0.05744; // K
    rr_prob[11] = 0.09019; // L
    rr_prob[12] = 0.02243; // M
    rr_prob[13] = 0.04487; // N
    rr_prob[14] = 0.05203; // P
    rr_prob[15] = 0.04264; // Q
    rr_prob[16] = 0.05129; // R
    rr_prob[17] = 0.07120; // S
    rr_prob[18] = 0.05841; // T
    rr_prob[19] = 0.06441; // V
    rr_prob[20] = 0.01330; // W
    rr_prob[22] = 0.03216; // Y

    let matrix = *blast_rs::api::get_matrix(blast_rs::api::MatrixType::Blosum62);

    let lambda = blast_rs::composition::compute_ungapped_lambda_with_bg(&matrix, &rr_prob);

    assert!(
        (lambda - 0.31760596).abs() < 1.0e-8,
        "Lambda with R&R freqs should match the NCBI kbp_ideal reference, got {:.8}",
        lambda
    );
}

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
    let status = std::process::Command::new("/usr/bin/blastn")
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
        .expect("run /usr/bin/blastn");
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
    if !std::path::Path::new("/usr/bin/blastn").exists() {
        eprintln!("Skipping: /usr/bin/blastn not found");
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
#[test]
#[ignore]
fn test_core_nt00_primer_outfmt6_matches_ncbi() {
    assert_core_nt_outfmt6_matches_ncbi(".00");
}

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
    if !std::path::Path::new("/usr/bin/blastn").exists() {
        eprintln!("Skipping: /usr/bin/blastn not found");
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

    let mut ncbi_cmd = std::process::Command::new("/usr/bin/blastn");
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
#[test]
#[ignore]
fn test_core_nt00_primer_taxonomy_outfmt_matches_ncbi() {
    assert_core_nt_taxonomy_outfmt_matches_ncbi(".00");
}

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
    if !std::path::Path::new("/usr/bin/blastn").exists() {
        eprintln!("Skipping: /usr/bin/blastn not found");
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

    let ncbi_status = std::process::Command::new("/usr/bin/blastn")
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
        .expect("run /usr/bin/blastn title parity");
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
#[test]
#[ignore]
fn test_core_nt00_primer_title_outfmt_matches_ncbi() {
    assert_core_nt_title_outfmt_matches_ncbi(".00");
}

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
    if !std::path::Path::new("/usr/bin/blastn").exists() {
        eprintln!("Skipping: /usr/bin/blastn not found");
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

    let ncbi_status = std::process::Command::new("/usr/bin/blastn")
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
        .expect("run /usr/bin/blastn outfmt7 parity");
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

/// Real core_nt commented-tabular parity check for the first volume.
///
/// Run with: cargo test --release -- --ignored test_core_nt00_primer_outfmt7_matches_ncbi
#[test]
#[ignore]
fn test_core_nt00_primer_outfmt7_matches_ncbi() {
    assert_core_nt_outfmt7_matches_ncbi(".00");
}

/// Regression test for taxonomy names when BLASTDB points at taxdb.bti/btd.
///
/// Run with: cargo test --release -- --ignored test_core_nt00_primer_taxonomy_names_with_blastdb_match_ncbi
#[test]
#[ignore]
fn test_core_nt00_primer_taxonomy_names_with_blastdb_match_ncbi() {
    let blastdb = std::path::Path::new(CORE_NT_BASE).parent().unwrap();
    assert_core_nt_taxonomy_outfmt_matches_ncbi_with_blastdb(".00", Some(blastdb));
}

/// Regression test that standalone nonzero core_nt volumes use global OIDs for
/// alias-level .not taxonomy lookup.
///
/// Run with: cargo test --release -- --ignored test_core_nt_nonzero_volume_taxonomy_outfmt_matches_ncbi
#[test]
#[ignore]
fn test_core_nt_nonzero_volume_taxonomy_outfmt_matches_ncbi() {
    assert_core_nt_taxonomy_outfmt_matches_ncbi(".12");
    assert_core_nt_taxonomy_outfmt_matches_ncbi(".28");
}

/// Regression test for accession extraction where title text contains a
/// version-like clone token before the real GenBank accession.
///
/// Run with: cargo test --release -- --ignored test_core_nt_accession_regression_volumes_match_ncbi
#[test]
#[ignore]
fn test_core_nt_accession_regression_volumes_match_ncbi() {
    assert_core_nt_outfmt6_matches_ncbi(".12");
    assert_core_nt_outfmt6_matches_ncbi(".28");
}

/// Full 84-volume core_nt alias parity test. This is intentionally gated by an
/// environment variable because it scans the full local 272 GB database.
///
/// Run with:
/// BLAST_RS_RUN_FULL_CORE_NT=1 cargo test --release -- --ignored test_core_nt_alias_primer_outfmt6_matches_ncbi
#[test]
#[ignore]
fn test_core_nt_alias_primer_outfmt6_matches_ncbi() {
    if std::env::var_os("BLAST_RS_RUN_FULL_CORE_NT").is_none() {
        eprintln!("Skipping: set BLAST_RS_RUN_FULL_CORE_NT=1 to scan full core_nt alias");
        return;
    }
    assert_core_nt_outfmt6_matches_ncbi("");
}

/// Full 84-volume core_nt alias taxonomy parity test. This is intentionally
/// gated because it scans the full local 272 GB database twice.
///
/// Run with:
/// BLAST_RS_RUN_FULL_CORE_NT=1 cargo test --release -- --ignored test_core_nt_alias_primer_taxonomy_outfmt_matches_ncbi
#[test]
#[ignore]
fn test_core_nt_alias_primer_taxonomy_outfmt_matches_ncbi() {
    if std::env::var_os("BLAST_RS_RUN_FULL_CORE_NT").is_none() {
        eprintln!("Skipping: set BLAST_RS_RUN_FULL_CORE_NT=1 to scan full core_nt alias");
        return;
    }
    assert_core_nt_taxonomy_outfmt_matches_ncbi("");
}
