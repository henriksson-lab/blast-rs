use std::collections::BTreeSet;
use std::path::{Path, PathBuf};
use std::process::Command;

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn require_file(path: &Path) -> bool {
    if path.exists() {
        true
    } else {
        eprintln!("SKIP: missing {}", path.display());
        false
    }
}

fn ncbi_tool(name: &str) -> Option<PathBuf> {
    let path = repo_root().join(format!("ncbi-blast-2.17.0+-src/c++/ReleaseMT/bin/{name}"));
    if require_file(&path) {
        Some(path)
    } else {
        None
    }
}

fn blast_cli() -> Option<PathBuf> {
    if let Ok(path) = std::env::var("CARGO_BIN_EXE_blast-cli") {
        let path = PathBuf::from(path);
        if require_file(&path) {
            return Some(path);
        }
    }
    let path = repo_root().join("target/release/blast-cli");
    if require_file(&path) {
        Some(path)
    } else {
        None
    }
}

fn run_lines(program: &Path, args: &[&str]) -> Vec<String> {
    let output = Command::new(program)
        .args(args)
        .current_dir(repo_root())
        .output()
        .unwrap_or_else(|err| panic!("failed to run {}: {err}", program.display()));
    assert!(
        output.status.success(),
        "{} failed with status {}\nstderr:\n{}",
        program.display(),
        output.status,
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8_lossy(&output.stdout)
        .lines()
        .filter(|line| !line.trim().is_empty() && !line.starts_with('#'))
        .map(str::to_owned)
        .collect()
}

fn set(lines: Vec<String>) -> BTreeSet<String> {
    lines.into_iter().collect()
}

fn diff_counts(a: &BTreeSet<String>, b: &BTreeSet<String>) -> (usize, usize) {
    (a.difference(b).count(), b.difference(a).count())
}

fn query_coordinate_set(lines: Vec<String>) -> BTreeSet<String> {
    lines
        .into_iter()
        .map(|line| {
            let cols: Vec<&str> = line.split('\t').collect();
            assert!(cols.len() >= 6, "outfmt row has too few columns: {line}");
            // Drop sseqid. blast-rs and NCBI can format equivalent protein IDs
            // differently for the same sequence; this test isolates hit geometry.
            format!(
                "{}\t{}\t{}\t{}\t{}",
                cols[0], cols[2], cols[3], cols[4], cols[5]
            )
        })
        .collect()
}

fn prefix_set(lines: Vec<String>, ncols: usize) -> BTreeSet<String> {
    lines
        .into_iter()
        .map(|line| {
            let cols: Vec<&str> = line.split('\t').collect();
            assert!(
                cols.len() >= ncols,
                "outfmt row has too few columns for {ncols}-column key: {line}"
            );
            cols[..ncols].join("\t")
        })
        .collect()
}

#[test]
#[ignore = "requires local NCBI 2.17 binaries and /tmp/bench real-data fixtures"]
fn real_data_megablast_nt20_coordinate_parity_guard() {
    let Some(ncbi) = ncbi_tool("blastn") else {
        return;
    };
    let Some(rs) = blast_cli() else {
        return;
    };
    if !require_file(Path::new("/tmp/bench/nt20.fa")) {
        return;
    }

    let outfmt = "6 qseqid sseqid qstart qend sstart send";
    let common = [
        "-query",
        "/tmp/bench/nt20.fa",
        "-db",
        "tests/fixtures/large_db/celegans",
        "-num_threads",
        "1",
        "-outfmt",
        outfmt,
        "-evalue",
        "10",
    ];
    let ncbi_set = set(run_lines(&ncbi, &common));
    let rs_set = set(run_lines(
        &rs,
        &["blastn"].into_iter().chain(common).collect::<Vec<_>>(),
    ));

    assert_eq!(
        ncbi_set, rs_set,
        "default megablast coordinate parity regressed"
    );
}

#[test]
#[ignore = "requires local NCBI 2.17 binaries and /tmp/bench real-data fixtures"]
fn real_data_regular_blastn_nt3_coordinate_parity_guard() {
    let Some(ncbi) = ncbi_tool("blastn") else {
        return;
    };
    let Some(rs) = blast_cli() else {
        return;
    };
    if !require_file(Path::new("/tmp/bench/nt3.fa")) {
        return;
    }

    let outfmt = "6 qseqid sseqid sstart send qstart qend";
    let ncbi_args = [
        "-task",
        "blastn",
        "-query",
        "/tmp/bench/nt3.fa",
        "-db",
        "tests/fixtures/large_db/celegans",
        "-num_threads",
        "1",
        "-outfmt",
        outfmt,
        "-evalue",
        "10",
    ];
    let rs_args = [
        "blastn",
        "-task",
        "blastn",
        "-query",
        "/tmp/bench/nt3.fa",
        "-db",
        "tests/fixtures/large_db/celegans",
        "-num_threads",
        "1",
        "-outfmt",
        outfmt,
        "-evalue",
        "10",
    ];
    let ncbi_set = set(run_lines(&ncbi, &ncbi_args));
    let rs_set = set(run_lines(&rs, &rs_args));

    assert_eq!(ncbi_set.len(), 3646);
    assert_eq!(rs_set.len(), 3646);
    assert_eq!(
        ncbi_set, rs_set,
        "regular blastn nt3 coordinate parity regressed"
    );
}

#[test]
#[ignore = "requires local NCBI 2.17 binaries and /tmp/bench real-data fixtures"]
fn real_data_blastp_cbs0_query_coordinate_parity_guard() {
    let Some(ncbi) = ncbi_tool("blastp") else {
        return;
    };
    let Some(rs) = blast_cli() else {
        return;
    };
    if !require_file(Path::new("/tmp/bench/prot30.fa")) {
        return;
    }

    let outfmt = "6 qseqid sseqid qstart qend sstart send";
    let common = [
        "-query",
        "/tmp/bench/prot30.fa",
        "-db",
        "tests/fixtures/seqp/seqp",
        "-num_threads",
        "1",
        "-comp_based_stats",
        "0",
        "-outfmt",
        outfmt,
        "-evalue",
        "10",
    ];
    let ncbi_set = query_coordinate_set(run_lines(&ncbi, &common));
    let rs_set = query_coordinate_set(run_lines(
        &rs,
        &["blastp"].into_iter().chain(common).collect::<Vec<_>>(),
    ));

    assert_eq!(
        ncbi_set, rs_set,
        "blastp cbs0 query-coordinate parity regressed"
    );
}

#[test]
#[ignore = "requires local NCBI 2.17 binaries and /tmp/bench real-data fixtures"]
fn real_data_blastp_cbs2_documents_current_composition_gap() {
    let Some(ncbi) = ncbi_tool("blastp") else {
        return;
    };
    let Some(rs) = blast_cli() else {
        return;
    };
    if !require_file(Path::new("/tmp/bench/prot30.fa")) {
        return;
    }

    let outfmt = "6 qseqid sseqid qstart qend sstart send";
    let common = [
        "-query",
        "/tmp/bench/prot30.fa",
        "-db",
        "tests/fixtures/seqp/seqp",
        "-num_threads",
        "1",
        "-comp_based_stats",
        "2",
        "-outfmt",
        outfmt,
        "-evalue",
        "10",
    ];
    let ncbi_set = query_coordinate_set(run_lines(&ncbi, &common));
    let rs_set = query_coordinate_set(run_lines(
        &rs,
        &["blastp"].into_iter().chain(common).collect::<Vec<_>>(),
    ));

    assert_eq!(ncbi_set.len(), 334);
    assert_eq!(rs_set.len(), 342);
    assert_eq!(diff_counts(&ncbi_set, &rs_set), (53, 61));
}

#[test]
#[ignore = "requires local NCBI 2.17 binaries and /tmp/bench real-data fixtures"]
fn real_data_blastx_nt20_seqp_cbs0_documents_current_coordinate_gap() {
    let Some(ncbi) = ncbi_tool("blastx") else {
        return;
    };
    let Some(rs) = blast_cli() else {
        return;
    };
    if !require_file(Path::new("/tmp/bench/nt20.fa")) {
        return;
    }

    let outfmt = "6 qseqid sseqid qstart qend sstart send qframe sframe";
    let common = [
        "-query",
        "/tmp/bench/nt20.fa",
        "-db",
        "tests/fixtures/seqp/seqp",
        "-num_threads",
        "1",
        "-comp_based_stats",
        "0",
        "-outfmt",
        outfmt,
        "-evalue",
        "10",
    ];
    let ncbi_set = set(run_lines(&ncbi, &common));
    let rs_set = set(run_lines(
        &rs,
        &["blastx"].into_iter().chain(common).collect::<Vec<_>>(),
    ));

    assert_eq!(ncbi_set.len(), 113);
    assert_eq!(rs_set.len(), 107);
    assert_eq!(diff_counts(&ncbi_set, &rs_set), (6, 0));
}

#[test]
#[ignore = "requires local NCBI 2.17 binaries and /tmp/bench real-data fixtures"]
fn real_data_blastx_nt20_seqp_cbs2_documents_current_composition_gap() {
    let Some(ncbi) = ncbi_tool("blastx") else {
        return;
    };
    let Some(rs) = blast_cli() else {
        return;
    };
    if !require_file(Path::new("/tmp/bench/nt20.fa")) {
        return;
    }

    let outfmt = "6 qseqid sseqid qstart qend sstart send qframe sframe";
    let common = [
        "-query",
        "/tmp/bench/nt20.fa",
        "-db",
        "tests/fixtures/seqp/seqp",
        "-num_threads",
        "1",
        "-comp_based_stats",
        "2",
        "-outfmt",
        outfmt,
        "-evalue",
        "10",
    ];
    let ncbi_set = set(run_lines(&ncbi, &common));
    let rs_set = set(run_lines(
        &rs,
        &["blastx"].into_iter().chain(common).collect::<Vec<_>>(),
    ));

    assert_eq!(ncbi_set.len(), 73);
    assert_eq!(rs_set.len(), 82);
    assert_eq!(diff_counts(&ncbi_set, &rs_set), (26, 35));
}

#[test]
#[ignore = "requires local NCBI 2.17 binaries and /tmp/bench real-data fixtures"]
fn real_data_tblastn_prot30_seqn_cbs0_documents_current_coordinate_gap() {
    let Some(ncbi) = ncbi_tool("tblastn") else {
        return;
    };
    let Some(rs) = blast_cli() else {
        return;
    };
    if !require_file(Path::new("/tmp/bench/prot30.fa")) {
        return;
    }

    let outfmt = "6 qseqid sseqid qstart qend sstart send qframe sframe";
    let common = [
        "-query",
        "/tmp/bench/prot30.fa",
        "-db",
        "tests/fixtures/seqn/seqn",
        "-num_threads",
        "1",
        "-comp_based_stats",
        "0",
        "-outfmt",
        outfmt,
        "-evalue",
        "10",
    ];
    let ncbi_set = set(run_lines(&ncbi, &common));
    let rs_set = set(run_lines(
        &rs,
        &["tblastn"].into_iter().chain(common).collect::<Vec<_>>(),
    ));

    assert_eq!(ncbi_set.len(), 88);
    assert_eq!(rs_set.len(), 82);
    assert_eq!(diff_counts(&ncbi_set, &rs_set), (6, 0));
}

#[test]
#[ignore = "requires local NCBI 2.17 binaries and /tmp/bench real-data fixtures"]
fn real_data_tblastn_prot30_seqn_cbs2_documents_current_composition_gap() {
    let Some(ncbi) = ncbi_tool("tblastn") else {
        return;
    };
    let Some(rs) = blast_cli() else {
        return;
    };
    if !require_file(Path::new("/tmp/bench/prot30.fa")) {
        return;
    }

    let outfmt = "6 qseqid sseqid qstart qend sstart send qframe sframe";
    let common = [
        "-query",
        "/tmp/bench/prot30.fa",
        "-db",
        "tests/fixtures/seqn/seqn",
        "-num_threads",
        "1",
        "-comp_based_stats",
        "2",
        "-outfmt",
        outfmt,
        "-evalue",
        "10",
    ];
    let ncbi_set = set(run_lines(&ncbi, &common));
    let rs_set = set(run_lines(
        &rs,
        &["tblastn"].into_iter().chain(common).collect::<Vec<_>>(),
    ));

    assert_eq!(ncbi_set.len(), 71);
    assert_eq!(rs_set.len(), 79);
    assert_eq!(diff_counts(&ncbi_set, &rs_set), (28, 36));
}

#[test]
#[ignore = "requires local NCBI 2.17 binaries and /tmp/bench real-data fixtures"]
fn real_data_tblastx_nt20_seqn_documents_current_coordinate_gap() {
    let Some(ncbi) = ncbi_tool("tblastx") else {
        return;
    };
    let Some(rs) = blast_cli() else {
        return;
    };
    if !require_file(Path::new("/tmp/bench/nt20.fa")) {
        return;
    }

    let outfmt = "6 qseqid sseqid qstart qend sstart send qframe sframe";
    let common = [
        "-query",
        "/tmp/bench/nt20.fa",
        "-db",
        "tests/fixtures/seqn/seqn",
        "-num_threads",
        "1",
        "-outfmt",
        outfmt,
        "-evalue",
        "10",
    ];
    let ncbi_set = set(run_lines(&ncbi, &common));
    let rs_set = set(run_lines(
        &rs,
        &["tblastx"].into_iter().chain(common).collect::<Vec<_>>(),
    ));

    assert_eq!(ncbi_set.len(), 458);
    assert_eq!(rs_set.len(), 465);
    assert_eq!(diff_counts(&ncbi_set, &rs_set), (8, 15));
}

#[test]
#[ignore = "requires local NCBI 2.17 binaries and /tmp/bench real-data fixtures"]
fn real_data_tblastx_nt10_seqn_documents_current_coordinate_gap() {
    let Some(ncbi) = ncbi_tool("tblastx") else {
        return;
    };
    let Some(rs) = blast_cli() else {
        return;
    };
    if !require_file(Path::new("/tmp/bench/nt10.fa")) {
        return;
    }

    let outfmt = "6 qseqid sseqid qstart qend sstart send qframe sframe";
    let common = [
        "-query",
        "/tmp/bench/nt10.fa",
        "-db",
        "tests/fixtures/seqn/seqn",
        "-num_threads",
        "1",
        "-outfmt",
        outfmt,
        "-evalue",
        "10",
    ];
    let ncbi_set = set(run_lines(&ncbi, &common));
    let rs_set = set(run_lines(
        &rs,
        &["tblastx"].into_iter().chain(common).collect::<Vec<_>>(),
    ));

    assert_eq!(ncbi_set.len(), 263);
    assert_eq!(rs_set.len(), 266);
    assert_eq!(diff_counts(&ncbi_set, &rs_set), (7, 10));
}

#[test]
#[ignore = "requires local NCBI 2.17 binaries and /tmp/bench real-data fixtures"]
fn real_data_regular_blastn_nt10_yeast_documents_current_overreport_gap() {
    let Some(ncbi) = ncbi_tool("blastn") else {
        return;
    };
    let Some(rs) = blast_cli() else {
        return;
    };
    if !require_file(Path::new("/tmp/bench/nt10.fa")) {
        return;
    }

    let outfmt = "6 qseqid sseqid sstart send qstart qend";
    let common = [
        "-task",
        "blastn",
        "-query",
        "/tmp/bench/nt10.fa",
        "-db",
        "tests/fixtures/large_db/yeast",
        "-num_threads",
        "1",
        "-outfmt",
        outfmt,
        "-evalue",
        "10",
    ];
    let ncbi_set = set(run_lines(&ncbi, &common));
    let rs_set = set(run_lines(
        &rs,
        &["blastn"].into_iter().chain(common).collect::<Vec<_>>(),
    ));

    assert_eq!(ncbi_set.len(), 418);
    assert_eq!(rs_set.len(), 418);
    assert_eq!(diff_counts(&ncbi_set, &rs_set), (0, 0));
}

#[test]
#[ignore = "requires local NCBI 2.17 binaries and /tmp/bench real-data fixtures"]
fn real_data_regular_blastn_nt20_yeast_documents_current_coordinate_gap() {
    let Some(ncbi) = ncbi_tool("blastn") else {
        return;
    };
    let Some(rs) = blast_cli() else {
        return;
    };
    if !require_file(Path::new("/tmp/bench/nt20.fa")) {
        return;
    }

    let outfmt = "6 qseqid sseqid sstart send qstart qend";
    let common = [
        "-task",
        "blastn",
        "-query",
        "/tmp/bench/nt20.fa",
        "-db",
        "tests/fixtures/large_db/yeast",
        "-num_threads",
        "1",
        "-outfmt",
        outfmt,
        "-evalue",
        "10",
    ];
    let ncbi_set = set(run_lines(&ncbi, &common));
    let rs_set = set(run_lines(
        &rs,
        &["blastn"].into_iter().chain(common).collect::<Vec<_>>(),
    ));

    assert_eq!(ncbi_set.len(), 769);
    assert_eq!(rs_set.len(), 769);
    assert_eq!(diff_counts(&ncbi_set, &rs_set), (0, 0));
}

#[test]
#[ignore = "requires local NCBI 2.17 binaries and /tmp/bench real-data fixtures"]
fn real_data_regular_blastn_nt20_pombe_documents_current_overreport_gap() {
    let Some(ncbi) = ncbi_tool("blastn") else {
        return;
    };
    let Some(rs) = blast_cli() else {
        return;
    };
    if !require_file(Path::new("/tmp/bench/nt20.fa")) {
        return;
    }

    let outfmt = "6 qseqid sseqid sstart send qstart qend";
    let common = [
        "-task",
        "blastn",
        "-query",
        "/tmp/bench/nt20.fa",
        "-db",
        "tests/fixtures/pombe/pombe",
        "-num_threads",
        "1",
        "-outfmt",
        outfmt,
        "-evalue",
        "10",
    ];
    let ncbi_set = set(run_lines(&ncbi, &common));
    let rs_set = set(run_lines(
        &rs,
        &["blastn"].into_iter().chain(common).collect::<Vec<_>>(),
    ));

    assert_eq!(ncbi_set.len(), 473);
    assert_eq!(rs_set.len(), 473);
    assert_eq!(diff_counts(&ncbi_set, &rs_set), (0, 0));
}

#[test]
#[ignore = "requires local NCBI 2.17 binaries and scratch real-data fixtures"]
fn real_data_blastp_caa85678_cbs2_documents_score_drift_without_coordinate_gap() {
    let Some(ncbi) = ncbi_tool("blastp") else {
        return;
    };
    let Some(rs) = blast_cli() else {
        return;
    };
    if !require_file(Path::new("scratch/parity15/CAA85678.1.fa")) {
        return;
    }

    let outfmt = "6 qseqid sseqid qstart qend sstart send bitscore evalue length pident";
    let common = [
        "-query",
        "scratch/parity15/CAA85678.1.fa",
        "-db",
        "tests/fixtures/seqp/seqp",
        "-num_threads",
        "1",
        "-comp_based_stats",
        "2",
        "-outfmt",
        outfmt,
        "-evalue",
        "10",
    ];
    let ncbi_lines = run_lines(&ncbi, &common);
    let rs_lines = run_lines(
        &rs,
        &["blastp"].into_iter().chain(common).collect::<Vec<_>>(),
    );
    let ncbi_coords = prefix_set(ncbi_lines.clone(), 6);
    let rs_coords = prefix_set(rs_lines.clone(), 6);
    let ncbi_full = set(ncbi_lines);
    let rs_full = set(rs_lines);

    assert_eq!(ncbi_coords, rs_coords);
    assert_eq!(ncbi_full.len(), 13);
    assert_eq!(rs_full.len(), 13);
    assert_eq!(diff_counts(&ncbi_full, &rs_full), (4, 4));
}

// ---------------------------------------------------------------------------
// Regression tests for non-parity cases found in the 2026-05-29 real-data
// sweep of all subtools. Each documents current behaviour against NCBI 2.17;
// when a gap is closed the corresponding assertion should be tightened.
// ---------------------------------------------------------------------------

/// Run a program and return (exit_code, non-comment stdout lines) WITHOUT
/// asserting success — for error/exit-code parity checks.
fn run_status_and_lines(program: &Path, args: &[&str]) -> (i32, BTreeSet<String>) {
    let output = Command::new(program)
        .args(args)
        .current_dir(repo_root())
        .output()
        .unwrap_or_else(|err| panic!("failed to run {}: {err}", program.display()));
    let code = output.status.code().unwrap_or(-1);
    let lines = String::from_utf8_lossy(&output.stdout)
        .lines()
        .filter(|line| !line.trim().is_empty() && !line.starts_with('#'))
        .map(str::to_owned)
        .collect();
    (code, lines)
}

/// Returns true if `program args` exits on its own within `secs`, false if it
/// had to be killed (i.e. it hung).
fn run_completes_within(program: &Path, args: &[&str], secs: u64) -> bool {
    let mut child = Command::new(program)
        .args(args)
        .current_dir(repo_root())
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null())
        .spawn()
        .unwrap_or_else(|err| panic!("failed to spawn {}: {err}", program.display()));
    let start = std::time::Instant::now();
    loop {
        if child.try_wait().expect("try_wait").is_some() {
            return true;
        }
        if start.elapsed() > std::time::Duration::from_secs(secs) {
            let _ = child.kill();
            let _ = child.wait();
            return false;
        }
        std::thread::sleep(std::time::Duration::from_millis(200));
    }
}

/// Map coordinate key (qseqid, qstart, qend, sstart, send) -> (bitscore, evalue)
/// from an `-outfmt "6 qseqid sseqid qstart qend sstart send bitscore evalue"`.
fn coord_bitscore_evalue(lines: Vec<String>) -> std::collections::BTreeMap<String, (String, String)> {
    lines
        .into_iter()
        .filter_map(|line| {
            let c: Vec<&str> = line.split('\t').collect();
            if c.len() < 8 {
                return None;
            }
            Some((
                format!("{}\t{}\t{}\t{}\t{}", c[0], c[2], c[3], c[4], c[5]),
                (c[6].to_string(), c[7].to_string()),
            ))
        })
        .collect()
}

/// blastp `-matrix PAM250` defaults to gapopen/gapextend 15/2 (was wrongly 14/2,
/// which produced incorrect Karlin-Altschul lambda/K — self-hit bitscore 1528 vs
/// NCBI 1603). After the fix, cbs0 coordinate parity collapses to the borderline
/// score-floor residual. Guards the fix.
#[test]
#[ignore = "requires local NCBI 2.17 binaries and /tmp/bench real-data fixtures"]
fn real_data_blastp_pam250_default_gap_parity_guard() {
    let Some(ncbi) = ncbi_tool("blastp") else {
        return;
    };
    let Some(rs) = blast_cli() else {
        return;
    };
    if !require_file(Path::new("/tmp/bench/prot30.fa")) {
        return;
    }

    let outfmt = "6 qseqid sseqid qstart qend sstart send";
    let common = [
        "-query",
        "/tmp/bench/prot30.fa",
        "-db",
        "tests/fixtures/seqp/seqp",
        "-num_threads",
        "1",
        "-matrix",
        "PAM250",
        "-comp_based_stats",
        "0",
        "-outfmt",
        outfmt,
        "-evalue",
        "10",
    ];
    let ncbi_set = query_coordinate_set(run_lines(&ncbi, &common));
    let rs_set = query_coordinate_set(run_lines(
        &rs,
        &["blastp"].into_iter().chain(common).collect::<Vec<_>>(),
    ));

    // PAM250 self-hit bitscore must match NCBI (1603) — the cleanest K-A probe.
    let rs_self = run_lines(
        &rs,
        &["blastp"]
            .into_iter()
            .chain([
                "-query",
                "/tmp/bench/prot30.fa",
                "-db",
                "tests/fixtures/seqp/seqp",
                "-num_threads",
                "1",
                "-matrix",
                "PAM250",
                "-comp_based_stats",
                "0",
                "-outfmt",
                "6 qseqid bitscore",
                "-evalue",
                "10",
            ])
            .collect::<Vec<_>>(),
    );
    let max_bits = rs_self
        .iter()
        .filter_map(|l| l.split('\t').nth(1))
        .filter_map(|b| b.parse::<f64>().ok())
        .fold(0.0_f64, f64::max);
    assert!(
        (max_bits - 1603.0).abs() < 1.0,
        "PAM250 default-gap self-hit bitscore drifted: {max_bits} (expected ~1603)"
    );
    assert_eq!(diff_counts(&ncbi_set, &rs_set), (1, 3));
}

/// blastn `-ungapped` near-parity. Was a massive over-report (≈3× on yeast/nt10,
/// ≈47× on celegans) because DUST query-masking was disabled in `-ungapped` mode.
/// Fixed 2026-05-29 (`src/main.rs`: DUST stays on as soft masking for ungapped,
/// matching NCBI). Now full recall; tiny residual at the e-value cutoff edge.
#[test]
#[ignore = "requires local NCBI 2.17 binaries and /tmp/bench real-data fixtures"]
fn real_data_blastn_ungapped_near_parity() {
    let Some(ncbi) = ncbi_tool("blastn") else {
        return;
    };
    let Some(rs) = blast_cli() else {
        return;
    };
    if !require_file(Path::new("/tmp/bench/nt10.fa")) {
        return;
    }

    let outfmt = "6 qseqid sseqid qstart qend sstart send";
    let common = [
        "-task",
        "blastn",
        "-ungapped",
        "-query",
        "/tmp/bench/nt10.fa",
        "-db",
        "tests/fixtures/large_db/yeast",
        "-num_threads",
        "1",
        "-outfmt",
        outfmt,
        "-evalue",
        "10",
    ];
    let ncbi_set = set(run_lines(&ncbi, &common));
    let rs_set = set(run_lines(
        &rs,
        &["blastn"].into_iter().chain(common).collect::<Vec<_>>(),
    ));

    assert_eq!(ncbi_set.len(), 550);
    let (only_ncbi, only_rs) = diff_counts(&ncbi_set, &rs_set);
    assert_eq!(only_ncbi, 0, "rs must recall every NCBI ungapped HSP");
    // Residual: 2 borderline 24-bp HSPs at the evalue=10 cutoff edge.
    assert!(
        only_rs <= 2,
        "blastn -ungapped over-report regressed ({only_rs}); expected <=2 borderline"
    );
}

/// blastn `-task dc-megablast` (discontiguous seeding). Was (2,0) recall+over-report
/// on yeast/nt20; fixed 2026-05-29 (`src/search.rs`: two-hit extension with the
/// disc window size 40, and removal of an over-strict score floor) → exact parity.
/// (celegans still has a deeper residual recall gap on the large genome.)
#[test]
#[ignore = "requires local NCBI 2.17 binaries and /tmp/bench real-data fixtures"]
fn real_data_blastn_dc_megablast_yeast_parity() {
    let Some(ncbi) = ncbi_tool("blastn") else {
        return;
    };
    let Some(rs) = blast_cli() else {
        return;
    };
    if !require_file(Path::new("/tmp/bench/nt20.fa")) {
        return;
    }

    let outfmt = "6 qseqid sseqid qstart qend sstart send";
    let common = [
        "-task",
        "dc-megablast",
        "-query",
        "/tmp/bench/nt20.fa",
        "-db",
        "tests/fixtures/large_db/yeast",
        "-num_threads",
        "1",
        "-outfmt",
        outfmt,
        "-evalue",
        "10",
    ];
    let ncbi_set = set(run_lines(&ncbi, &common));
    let rs_set = set(run_lines(
        &rs,
        &["blastn"].into_iter().chain(common).collect::<Vec<_>>(),
    ));

    assert_eq!(ncbi_set.len(), 7);
    assert_eq!(rs_set.len(), 7);
    assert_eq!(diff_counts(&ncbi_set, &rs_set), (0, 0));
}

/// Translated-search e-value drift: on tblastn HSPs with IDENTICAL coordinates
/// and bitscore, rs computes a different (lower) e-value than NCBI. Bitscore
/// parity proves lambda/K are right; the divergence is in the per-query
/// effective-length / sum-statistics term. Documents the systemic bug.
#[test]
#[ignore = "requires local NCBI 2.17 binaries and /tmp/bench real-data fixtures"]
fn real_data_tblastn_evalue_drift_on_shared_hsps() {
    let Some(ncbi) = ncbi_tool("tblastn") else {
        return;
    };
    let Some(rs) = blast_cli() else {
        return;
    };
    if !require_file(Path::new("/tmp/bench/prot30.fa")) {
        return;
    }

    let outfmt = "6 qseqid sseqid qstart qend sstart send bitscore evalue";
    let common = [
        "-query",
        "/tmp/bench/prot30.fa",
        "-db",
        "tests/fixtures/seqn/seqn",
        "-num_threads",
        "1",
        "-comp_based_stats",
        "0",
        "-outfmt",
        outfmt,
        "-evalue",
        "10",
    ];
    let ncbi_map = coord_bitscore_evalue(run_lines(&ncbi, &common));
    let rs_map = coord_bitscore_evalue(run_lines(
        &rs,
        &["tblastn"].into_iter().chain(common).collect::<Vec<_>>(),
    ));

    let mut shared = 0usize;
    let mut bit_mismatch = 0usize;
    let mut eval_mismatch = 0usize;
    for (key, (nb, ne)) in &ncbi_map {
        if let Some((rb, re)) = rs_map.get(key) {
            shared += 1;
            if nb != rb {
                bit_mismatch += 1;
            }
            if ne != re {
                eval_mismatch += 1;
            }
        }
    }
    assert!(shared > 10, "need a meaningful shared-HSP sample, got {shared}");
    assert_eq!(bit_mismatch, 0, "bitscores must match on shared HSPs");
    // Fixed 2026-05-29: tblastn now feeds the per-frame translated protein
    // subject length into the Spouge e-value (api.rs), matching NCBI exactly.
    assert_eq!(
        eval_mismatch, 0,
        "tblastn e-values must match NCBI on shared HSPs ({eval_mismatch}/{shared} drift)"
    );
}

/// blastx exhibits the same translated-search e-value drift as tblastn.
#[test]
#[ignore = "requires local NCBI 2.17 binaries and /tmp/bench real-data fixtures"]
fn real_data_blastx_evalue_drift_on_shared_hsps() {
    let Some(ncbi) = ncbi_tool("blastx") else {
        return;
    };
    let Some(rs) = blast_cli() else {
        return;
    };
    if !require_file(Path::new("/tmp/bench/nt10.fa")) {
        return;
    }

    let outfmt = "6 qseqid sseqid qstart qend sstart send bitscore evalue";
    let common = [
        "-query",
        "/tmp/bench/nt10.fa",
        "-db",
        "tests/fixtures/seqp/seqp",
        "-num_threads",
        "1",
        "-comp_based_stats",
        "0",
        "-outfmt",
        outfmt,
        "-evalue",
        "10",
    ];
    let ncbi_map = coord_bitscore_evalue(run_lines(&ncbi, &common));
    let rs_map = coord_bitscore_evalue(run_lines(
        &rs,
        &["blastx"].into_iter().chain(common).collect::<Vec<_>>(),
    ));

    let mut shared = 0usize;
    let mut eval_mismatch = 0usize;
    for (key, (_nb, ne)) in &ncbi_map {
        if let Some((_rb, re)) = rs_map.get(key) {
            shared += 1;
            if ne != re {
                eval_mismatch += 1;
            }
        }
    }
    assert!(shared > 10, "need a meaningful shared-HSP sample, got {shared}");
    // Fixed 2026-05-29: blastx singleton e-values now use the per-frame Spouge
    // value (api.rs num>1 gating), matching NCBI exactly.
    assert_eq!(
        eval_mismatch, 0,
        "blastx e-values must match NCBI on shared HSPs ({eval_mismatch}/{shared} drift)"
    );
}

/// NCBI rejects `-max_intron_length` on tblastx ("Uneven gap linking of HSPs is
/// allowed for blastx, tblastn, and psitblastn only"); blast-rs currently
/// accepts it (and even alters output via HSP linking). Documents the
/// error-parity gap.
#[test]
#[ignore = "requires local NCBI 2.17 binaries and /tmp/bench real-data fixtures"]
fn real_data_tblastx_should_reject_max_intron_length_like_ncbi() {
    let Some(ncbi) = ncbi_tool("tblastx") else {
        return;
    };
    let Some(rs) = blast_cli() else {
        return;
    };
    if !require_file(Path::new("/tmp/bench/nt10.fa")) {
        return;
    }

    let args = [
        "-query",
        "/tmp/bench/nt10.fa",
        "-db",
        "tests/fixtures/seqn/seqn",
        "-num_threads",
        "1",
        "-max_intron_length",
        "100",
        "-outfmt",
        "6",
    ];
    let (ncbi_code, _) = run_status_and_lines(&ncbi, &args);
    let (rs_code, _) = run_status_and_lines(
        &rs,
        &["tblastx"].into_iter().chain(args).collect::<Vec<_>>(),
    );

    assert_ne!(ncbi_code, 0, "NCBI tblastx must reject -max_intron_length");
    // Fixed 2026-05-29: blast-rs now rejects -max_intron_length on tblastx too.
    assert_ne!(
        rs_code, 0,
        "blast-rs must reject -max_intron_length on tblastx (matches NCBI)"
    );
}

/// `-outfmt 4` (flat query-anchored). Implemented 2026-05-29
/// (`src/format/flat_query_anchored.rs`); blastn output is byte-identical to NCBI
/// 2.17. (blastp/tblastn match in layout modulo the pre-existing protein score
/// drift; blastx/tblastx outfmt-4 are not yet supported — translated-query
/// per-frame grouping.)
#[test]
#[ignore = "requires local NCBI 2.17 binaries and /tmp/bench real-data fixtures"]
fn real_data_blastn_outfmt4_matches_ncbi() {
    let Some(ncbi) = ncbi_tool("blastn") else {
        return;
    };
    let Some(rs) = blast_cli() else {
        return;
    };
    if !require_file(Path::new("/tmp/bench/nt20.fa")) {
        return;
    }

    let args = [
        "-query",
        "/tmp/bench/nt20.fa",
        "-db",
        "tests/fixtures/large_db/yeast",
        "-num_threads",
        "1",
        "-outfmt",
        "4",
        "-evalue",
        "10",
    ];
    let ncbi_out = Command::new(&ncbi)
        .args(args)
        .current_dir(repo_root())
        .output()
        .expect("run ncbi blastn -outfmt 4");
    let rs_out = Command::new(&rs)
        .args(["blastn"].into_iter().chain(args))
        .current_dir(repo_root())
        .output()
        .expect("run rs blastn -outfmt 4");

    assert!(ncbi_out.status.success(), "NCBI supports -outfmt 4");
    assert!(rs_out.status.success(), "blast-rs now implements -outfmt 4");
    assert_eq!(
        String::from_utf8_lossy(&rs_out.stdout),
        String::from_utf8_lossy(&ncbi_out.stdout),
        "blastn -outfmt 4 must be byte-identical to NCBI 2.17"
    );
}

/// blastp `-word_size 6` used to HANG (the dense 28-letter neighborhood table
/// explodes for word_size>=5). Fixed 2026-05-29 by routing word_size>=5 through
/// NCBI's compressed-alphabet lookup (protein_lookup.rs). Now completes with full
/// recall vs NCBI; a residual over-report remains (weak hits the compressed
/// seeding keeps — same family as the protein score-floor over-report at ws3).
#[test]
#[ignore = "requires local NCBI 2.17 binaries and /tmp/bench real-data fixtures"]
fn real_data_blastp_word_size_6_completes_with_full_recall() {
    let Some(ncbi) = ncbi_tool("blastp") else {
        return;
    };
    let Some(rs) = blast_cli() else {
        return;
    };
    if !require_file(Path::new("/tmp/bench/prot30.fa")) {
        return;
    }

    let outfmt = "6 qseqid sseqid qstart qend sstart send";
    let common = [
        "-query",
        "/tmp/bench/prot30.fa",
        "-db",
        "tests/fixtures/seqp/seqp",
        "-num_threads",
        "1",
        "-word_size",
        "6",
        "-outfmt",
        outfmt,
        "-evalue",
        "10",
    ];
    // Completes (run_lines would block forever if it still hung).
    let ncbi_set = set(run_lines(&ncbi, &common));
    let rs_set = set(run_lines(
        &rs,
        &["blastp"].into_iter().chain(common).collect::<Vec<_>>(),
    ));
    let (only_ncbi, only_rs) = diff_counts(&ncbi_set, &rs_set);
    assert_eq!(only_ncbi, 0, "blastp -word_size 6 must recall every NCBI HSP");
    // Residual over-report (weak hits from the lossy compressed seeding).
    assert!(
        only_rs > 0,
        "over-report gone ({only_rs}) — tighten this test toward exact parity"
    );
}
