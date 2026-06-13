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

/// Map coordinate key (qseqid, qstart, qend, sstart, send) -> (bitscore, evalue)
/// from an `-outfmt "6 qseqid sseqid qstart qend sstart send bitscore evalue"`.
fn coord_bitscore_evalue(
    lines: Vec<String>,
) -> std::collections::BTreeMap<String, (String, String)> {
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
