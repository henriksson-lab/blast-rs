//! Performance tests for blast-rs on the real AMRFinderPlus search workloads.
//!
//! These compare `blast-cli` (blast-rs) against the matching system NCBI tool on
//! the actual AMR searches that `amrfinder` runs — a representative real-world
//! workload (large-ish protein DB, many hits) that this crate's own fixtures
//! don't cover. The workload data (the `AMRProt.fa` database and the query FASTAs)
//! lives in the sibling `ncbi-amrfinderplus-rs` repo and is NOT vendored here
//! (it would bloat the published crate), so these tests skip unless that repo is
//! present next to this one (or pointed at via env vars).
//!
//! They are `#[ignore]`d so a normal `cargo test` skips them. Run explicitly:
//!
//! ```text
//! cargo test --test amr_workload_perf -- --ignored --nocapture
//! ```
//!
//! Each skips gracefully (passing) if `blast-cli`, the NCBI tool, `/usr/bin/time`,
//! or the AMR workload are unavailable. Override via env vars:
//!   BLAST_RS_CLI_BIN=<path to blast-cli>  (default: CARGO_BIN_EXE_blast-cli, then target/release/blast-cli)
//!   BLAST_RS_PERF_AMR_DIR=<amrfinder-rs repo root>  (default: ../ncbi-amrfinderplus-rs)
//!   BLAST_RS_PERF_DB=<path to AMRProt.fa>            (overrides DB discovery)
//!   BLAST_RS_PERF_QUERY_PROT=<protein FASTA>         (default: <amr_dir>/tests/golden/test_prot.fa)
//!   BLAST_RS_PERF_QUERY_DNA=<nucleotide FASTA>       (default: <amr_dir>/tests/golden/test_dna.fa)
//!   NCBI_BLASTP / NCBI_BLASTX=<path>                 (default: /usr/bin/<tool>, then PATH)
//!   BLAST_RS_PERF_REPS=<n>                (default: 3; set 1 for a quick blastx check)
//!   BLAST_RS_PERF_THREADS=<n>             (default: 4)
//!   BLAST_RS_PERF_GENCODE=<n>             (default: 11, amrfinder's translation table)
//!   BLAST_RS_PERF_MIN_RECALL=<0..1>       (default: 0.95)

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::Instant;

const OUTFMT: &str = "6 sseqid qseqid sstart send slen qstart qend qlen sseq qseq";

fn env_or<T: std::str::FromStr>(key: &str, default: T) -> T {
    std::env::var(key)
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(default)
}

fn manifest_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// Root of the sibling amrfinder-rs repo that carries the AMR workload.
fn amr_dir() -> PathBuf {
    std::env::var_os("BLAST_RS_PERF_AMR_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|| manifest_dir().join("../ncbi-amrfinderplus-rs"))
}

/// Locate the AMR protein BLAST database (the directory holding `AMRProt.fa.pin`).
fn find_amr_db() -> Option<PathBuf> {
    if let Some(p) = std::env::var_os("BLAST_RS_PERF_DB") {
        let p = PathBuf::from(p);
        return p.exists().then_some(p);
    }
    let db_root = amr_dir().join("amrfinder_db");
    for e in std::fs::read_dir(&db_root).ok()?.flatten() {
        if e.path().join("AMRProt.fa.pin").exists() {
            return Some(e.path().join("AMRProt.fa"));
        }
    }
    None
}

fn find_blast_cli() -> Option<PathBuf> {
    if let Some(p) =
        std::env::var_os("BLAST_RS_CLI_BIN").or_else(|| std::env::var_os("CARGO_BIN_EXE_blast-cli"))
    {
        let p = PathBuf::from(p);
        if p.exists() {
            return Some(p);
        }
    }
    let rel = manifest_dir().join("target/release/blast-cli");
    rel.exists().then_some(rel)
}

/// Find an NCBI tool (`blastp`, `blastx`, ...): env override `NCBI_<TOOL>`, then
/// `/usr/bin/<tool>`, then a PATH lookup.
fn find_ncbi(tool: &str) -> Option<PathBuf> {
    let env_key = format!("NCBI_{}", tool.to_uppercase());
    if let Ok(p) = std::env::var(&env_key) {
        let p = PathBuf::from(p);
        return p.exists().then_some(p);
    }
    let usr = PathBuf::from(format!("/usr/bin/{tool}"));
    if usr.exists() {
        return Some(usr);
    }
    let out = Command::new("which").arg(tool).output().ok()?;
    if out.status.success() {
        let p = PathBuf::from(String::from_utf8_lossy(&out.stdout).trim());
        if p.exists() {
            return Some(p);
        }
    }
    None
}

struct Measure {
    best_secs: f64,
    peak_rss_kb: Option<u64>,
}

/// Run `prog args...` `reps` times under `/usr/bin/time -v`, returning the best
/// (minimum) wall time and the maximum observed peak RSS. Wall time is measured
/// in-process with `Instant`; RSS is parsed from `/usr/bin/time`'s report.
fn run_timed(prog: &Path, args: &[String], reps: usize) -> Measure {
    let have_time = Path::new("/usr/bin/time").exists();
    let mut best = f64::INFINITY;
    let mut peak: Option<u64> = None;

    for _ in 0..reps {
        let start = Instant::now();
        let output = if have_time {
            Command::new("/usr/bin/time")
                .arg("-v")
                .arg(prog)
                .args(args)
                .output()
        } else {
            Command::new(prog).args(args).output()
        }
        .expect("failed to spawn blast process");
        let secs = start.elapsed().as_secs_f64();

        assert!(
            output.status.success(),
            "blast process failed: {} {:?}\nstderr:\n{}",
            prog.display(),
            args,
            String::from_utf8_lossy(&output.stderr)
        );

        best = best.min(secs);
        if have_time {
            if let Some(rss) = parse_max_rss_kb(&String::from_utf8_lossy(&output.stderr)) {
                peak = Some(peak.map_or(rss, |p| p.max(rss)));
            }
        }
    }
    Measure {
        best_secs: best,
        peak_rss_kb: peak,
    }
}

fn parse_max_rss_kb(time_stderr: &str) -> Option<u64> {
    for line in time_stderr.lines() {
        if let Some(rest) = line.split("Maximum resident set size (kbytes):").nth(1) {
            return rest.trim().parse().ok();
        }
    }
    None
}

/// Parse `outfmt 6` output into the set of (query, subject-accession) hit pairs.
/// Column layout is `sseqid qseqid sstart send slen qstart qend qlen sseq qseq`
/// (amrfinder's order). The subject accession is the leading `|`-delimited token,
/// normalizing NCBI's full pipe-delimited defline to match blast-rs's accession.
fn hit_pairs(path: &Path) -> BTreeSet<(String, String)> {
    let text = std::fs::read_to_string(path).unwrap_or_default();
    text.lines()
        .filter_map(|line| {
            let mut cols = line.split('\t');
            let sseqid = cols.next()?;
            let qseqid = cols.next()?;
            let accession = sseqid.split('|').next().unwrap_or(sseqid);
            Some((qseqid.to_string(), accession.to_string()))
        })
        .collect()
}

/// Run both tools, report timing/RSS/hit-agreement, and assert the recall of NCBI
/// hits by blast-rs meets `min_recall`.
#[allow(clippy::too_many_arguments)]
fn perf_compare(
    label: &str,
    blast_cli: &Path,
    brs_args: &[String],
    brs_out: &Path,
    ncbi: &Path,
    ncbi_args: &[String],
    ncbi_out: &Path,
    reps: usize,
    min_recall: f64,
) {
    let brs = run_timed(blast_cli, brs_args, reps);
    let ncbi_m = run_timed(ncbi, ncbi_args, reps);

    let brs_pairs = hit_pairs(brs_out);
    let ncbi_pairs = hit_pairs(ncbi_out);
    let shared = brs_pairs.intersection(&ncbi_pairs).count();
    let only_ncbi = ncbi_pairs.len() - shared;
    let only_brs = brs_pairs.len() - shared;
    let union = brs_pairs.len() + ncbi_pairs.len() - shared;
    let recall = if ncbi_pairs.is_empty() {
        1.0
    } else {
        shared as f64 / ncbi_pairs.len() as f64
    };
    let precision = if brs_pairs.is_empty() {
        0.0
    } else {
        shared as f64 / brs_pairs.len() as f64
    };
    let jaccard = if union == 0 {
        1.0
    } else {
        shared as f64 / union as f64
    };

    let rss = |m: &Measure| match m.peak_rss_kb {
        Some(kb) => format!("{:8.1} MB", kb as f64 / 1024.0),
        None => "      n/a".to_string(),
    };
    let speedup = ncbi_m.best_secs / brs.best_secs;

    eprintln!("\n[amr_workload_perf] ===== {label} performance (best of {reps}) =====");
    eprintln!("  tool       best_time      peak_rss     hit_pairs");
    eprintln!(
        "  blast-rs  {:9.3}s   {}   {}",
        brs.best_secs,
        rss(&brs),
        brs_pairs.len()
    );
    eprintln!(
        "  NCBI      {:9.3}s   {}   {}",
        ncbi_m.best_secs,
        rss(&ncbi_m),
        ncbi_pairs.len()
    );
    eprintln!(
        "  -> blast-rs is {:.2}x {} than NCBI",
        if speedup >= 1.0 {
            speedup
        } else {
            1.0 / speedup
        },
        if speedup >= 1.0 { "faster" } else { "slower" }
    );
    eprintln!(
        "[amr_workload_perf] hit agreement: shared={shared} only_ncbi={only_ncbi} only_blast_rs={only_brs}"
    );
    eprintln!("  recall(of NCBI)={recall:.4}  precision={precision:.4}  jaccard={jaccard:.4}");

    assert!(
        recall >= min_recall,
        "{label}: blast-rs recall of NCBI hits {recall:.4} < required {min_recall:.4} \
         (shared={shared}, ncbi={}). blast-rs may be missing real hits.",
        ncbi_pairs.len()
    );
}

/// Shared inputs/binaries for the perf tests. Returns `None` (after logging a SKIP)
/// when a required input is missing, so the calling test passes as skipped.
struct Common {
    blast_cli: PathBuf,
    db: PathBuf,
    reps: usize,
    threads: usize,
    min_recall: f64,
}

fn common() -> Option<Common> {
    let blast_cli = match find_blast_cli() {
        Some(p) => p,
        None => {
            eprintln!("[amr_workload_perf] SKIP: blast-cli not found (set BLAST_RS_CLI_BIN or build --release)");
            return None;
        }
    };
    let db = match find_amr_db() {
        Some(d) => d,
        None => {
            eprintln!(
                "[amr_workload_perf] SKIP: AMR database not found under {} (set BLAST_RS_PERF_DB or BLAST_RS_PERF_AMR_DIR)",
                amr_dir().display()
            );
            return None;
        }
    };
    Some(Common {
        blast_cli,
        db,
        reps: env_or("BLAST_RS_PERF_REPS", 3),
        threads: env_or("BLAST_RS_PERF_THREADS", 4),
        min_recall: env_or("BLAST_RS_PERF_MIN_RECALL", 0.95),
    })
}

fn query_path(env_key: &str, default_name: &str) -> PathBuf {
    std::env::var_os(env_key)
        .map(PathBuf::from)
        .unwrap_or_else(|| amr_dir().join("tests/golden").join(default_name))
}

#[test]
#[ignore = "perf benchmark; run with --ignored --nocapture"]
fn blastp_amr_workload_perf_vs_ncbi() {
    let Some(c) = common() else { return };
    let query = query_path("BLAST_RS_PERF_QUERY_PROT", "test_prot.fa");
    if !query.exists() {
        eprintln!(
            "[amr_workload_perf] SKIP: protein query {} not found",
            query.display()
        );
        return;
    }
    let ncbi = match find_ncbi("blastp") {
        Some(p) => p,
        None => {
            eprintln!("[amr_workload_perf] SKIP: NCBI blastp not found (set NCBI_BLASTP or install BLAST+)");
            return;
        }
    };

    let tmp = tempfile::tempdir().expect("tempdir");
    let brs_out = tmp.path().join("blast_rs.tsv");
    let ncbi_out = tmp.path().join("ncbi.tsv");
    let q = query.display().to_string();
    let d = c.db.display().to_string();
    let th = c.threads.to_string();

    // blast-rs (blast-cli blastp) — long-flag interface; no `-task` (blastn-only).
    let brs_args: Vec<String> = vec![
        "blastp".into(),
        "-q".into(),
        q.clone(),
        "-d".into(),
        d.clone(),
        "--comp_based_stats".into(),
        "0".into(),
        "--seg".into(),
        "no".into(),
        "--max_target_seqs".into(),
        "10000".into(),
        "--dbsize".into(),
        "10000".into(),
        "--evalue".into(),
        "1e-10".into(),
        "--word_size".into(),
        "5".into(),
        "--num_threads".into(),
        th.clone(),
        "--outfmt".into(),
        OUTFMT.into(),
        "-o".into(),
        brs_out.display().to_string(),
    ];
    // NCBI blastp — matches amrfinder's blastp-fast invocation exactly.
    let ncbi_args: Vec<String> = vec![
        "-query".into(),
        q,
        "-db".into(),
        d,
        "-comp_based_stats".into(),
        "0".into(),
        "-seg".into(),
        "no".into(),
        "-max_target_seqs".into(),
        "10000".into(),
        "-dbsize".into(),
        "10000".into(),
        "-evalue".into(),
        "1e-10".into(),
        "-word_size".into(),
        "5".into(),
        "-task".into(),
        "blastp-fast".into(),
        "-num_threads".into(),
        th,
        "-outfmt".into(),
        OUTFMT.into(),
        "-out".into(),
        ncbi_out.display().to_string(),
    ];

    eprintln!(
        "[amr_workload_perf] workload: blastp {} vs {} | threads={} reps={}",
        query.display(),
        c.db.display(),
        c.threads,
        c.reps
    );
    perf_compare(
        "blastp",
        &c.blast_cli,
        &brs_args,
        &brs_out,
        &ncbi,
        &ncbi_args,
        &ncbi_out,
        c.reps,
        c.min_recall,
    );
}

#[test]
#[ignore = "perf benchmark; run with --ignored --nocapture"]
fn blastx_amr_workload_perf_vs_ncbi() {
    let Some(c) = common() else { return };
    let query = query_path("BLAST_RS_PERF_QUERY_DNA", "test_dna.fa");
    if !query.exists() {
        eprintln!(
            "[amr_workload_perf] SKIP: nucleotide query {} not found",
            query.display()
        );
        return;
    }
    let ncbi = match find_ncbi("blastx") {
        Some(p) => p,
        None => {
            eprintln!("[amr_workload_perf] SKIP: NCBI blastx not found (set NCBI_BLASTX or install BLAST+)");
            return;
        }
    };
    let gencode: u32 = env_or("BLAST_RS_PERF_GENCODE", 11);

    let tmp = tempfile::tempdir().expect("tempdir");
    let brs_out = tmp.path().join("blast_rs.tsv");
    let ncbi_out = tmp.path().join("ncbi.tsv");
    let q = query.display().to_string();
    let d = c.db.display().to_string();
    let th = c.threads.to_string();
    let gc = gencode.to_string();

    // blast-rs (blast-cli blastx): translated nucleotide query vs AMRProt protein DB.
    let brs_args: Vec<String> = vec![
        "blastx".into(),
        "-q".into(),
        q.clone(),
        "-d".into(),
        d.clone(),
        "--comp_based_stats".into(),
        "0".into(),
        "--seg".into(),
        "no".into(),
        "--max_target_seqs".into(),
        "10000".into(),
        "--dbsize".into(),
        "10000".into(),
        "--evalue".into(),
        "1e-10".into(),
        "--word_size".into(),
        "5".into(),
        "--query_gencode".into(),
        gc.clone(),
        "--num_threads".into(),
        th.clone(),
        "--outfmt".into(),
        OUTFMT.into(),
        "-o".into(),
        brs_out.display().to_string(),
    ];
    // NCBI blastx — matches amrfinder's blastx invocation (amrfinder-rs src/pipeline.rs).
    let ncbi_args: Vec<String> = vec![
        "-query".into(),
        q,
        "-db".into(),
        d,
        "-comp_based_stats".into(),
        "0".into(),
        "-seg".into(),
        "no".into(),
        "-max_target_seqs".into(),
        "10000".into(),
        "-dbsize".into(),
        "10000".into(),
        "-evalue".into(),
        "1e-10".into(),
        "-word_size".into(),
        "5".into(),
        "-query_gencode".into(),
        gc,
        "-num_threads".into(),
        th,
        "-outfmt".into(),
        OUTFMT.into(),
        "-out".into(),
        ncbi_out.display().to_string(),
    ];

    eprintln!(
        "[amr_workload_perf] workload: blastx {} vs {} | gencode={} threads={} reps={}",
        query.display(),
        c.db.display(),
        gencode,
        c.threads,
        c.reps
    );
    perf_compare(
        "blastx",
        &c.blast_cli,
        &brs_args,
        &brs_out,
        &ncbi,
        &ncbi_args,
        &ncbi_out,
        c.reps,
        c.min_recall,
    );
}
