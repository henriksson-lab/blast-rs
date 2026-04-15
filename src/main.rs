//! BLAST command-line interface.

use blast_rs::db::{BlastDb, DbType};
use blast_rs::format::{format_tabular, TabularHit};
use blast_rs::input::{iupacna_to_blastna, parse_fasta};
use clap::{Parser, Subcommand};
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::PathBuf;

#[derive(Parser)]
#[command(
    name = "blast-cli",
    version,
    about = "BLAST sequence search (Rust implementation)",
    allow_negative_numbers = true
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

fn parse_percent_0_to_100(value: &str) -> Result<f64, String> {
    let parsed = value
        .parse::<f64>()
        .map_err(|_| format!("invalid floating-point value: `{value}`"))?;
    if (0.0..=100.0).contains(&parsed) {
        Ok(parsed)
    } else {
        Err(format!("value must be between 0 and 100: `{value}`"))
    }
}

fn parse_best_hit_fraction(value: &str) -> Result<f64, String> {
    let parsed = value
        .parse::<f64>()
        .map_err(|_| format!("invalid floating-point value: `{value}`"))?;
    if parsed > 0.0 && parsed < 0.5 {
        Ok(parsed)
    } else {
        Err(format!("value must be greater than 0 and less than 0.5: `{value}`"))
    }
}

fn parse_template_length(value: &str) -> Result<i32, String> {
    let parsed = value
        .parse::<i32>()
        .map_err(|_| format!("invalid integer value: `{value}`"))?;
    if matches!(parsed, 16 | 18 | 21) {
        Ok(parsed)
    } else {
        Err(format!("value must be one of 16, 18, or 21: `{value}`"))
    }
}

#[derive(Subcommand)]
enum Commands {
    /// Nucleotide-nucleotide BLAST
    Blastn(BlastnArgs),
    /// Protein-protein BLAST
    Blastp(BlastnArgs),
    /// Translated nucleotide query vs protein subject
    Blastx(BlastnArgs),
    /// Protein query vs translated nucleotide subject
    Tblastn(BlastnArgs),
    /// Translated query vs translated subject
    Tblastx(BlastnArgs),
    /// Position-specific iterated BLAST
    Psiblast(BlastnArgs),
    /// Reverse position-specific BLAST
    Rpsblast(BlastnArgs),
    /// Translated RPS-BLAST
    Rpstblastn(BlastnArgs),
    /// Domain enhanced lookup time accelerated BLAST
    Deltablast(BlastnArgs),
}

#[derive(Parser, Clone)]
#[command(allow_negative_numbers = true)]
struct BlastnArgs {
    /// Query file in FASTA format
    #[arg(short, long)]
    query: PathBuf,

    /// BLAST database name (path without extension)
    #[arg(short, long)]
    db: Option<PathBuf>,

    /// Subject FASTA file (alternative to --db)
    #[arg(short, long)]
    subject: Option<PathBuf>,

    /// Output file (default: stdout)
    #[arg(short, long)]
    out: Option<PathBuf>,

    /// Expectation value threshold
    #[arg(short, long, default_value = "10.0")]
    evalue: f64,

    /// Number of threads
    #[arg(long = "num_threads", default_value = "1", value_parser = clap::value_parser!(i32).range(1..))]
    num_threads: i32,

    /// Output format: 0=pairwise, 5=XML, 6=tabular, 17=SAM, or '6 col1 col2...'
    #[arg(long = "outfmt", default_value = "6")]
    outfmt: String,

    /// Word size for initial seed
    #[arg(short = 'W', long = "word_size")]
    word_size: Option<i32>,

    /// Reward for nucleotide match
    #[arg(long, value_parser = clap::value_parser!(i32).range(0..))]
    reward: Option<i32>,

    /// Penalty for nucleotide mismatch
    #[arg(long, value_parser = clap::value_parser!(i32).range(..=0))]
    penalty: Option<i32>,

    /// Cost to open a gap
    #[arg(long)]
    gapopen: Option<i32>,

    /// Cost to extend a gap
    #[arg(long)]
    gapextend: Option<i32>,

    /// Query strand(s) to search: both, plus, minus
    #[arg(long, default_value = "both", value_parser = ["both", "plus", "minus"])]
    strand: String,

    /// Maximum number of target sequences to report
    #[arg(
        long = "max_target_seqs",
        alias = "max-target-seqs",
        conflicts_with_all = ["num_descriptions", "num_alignments"],
        value_parser = clap::value_parser!(i32).range(1..)
    )]
    max_target_seqs: Option<i32>,

    /// DUST low-complexity filtering: yes/no
    #[arg(long = "dust", default_value = "yes")]
    dust: String,

    /// X-dropoff for ungapped extensions
    #[arg(long = "xdrop_ungap", default_value = "20.0")]
    xdrop_ungap: f64,

    /// X-dropoff for preliminary gapped extensions
    #[arg(long = "xdrop_gap", default_value = "30.0")]
    xdrop_gap: f64,

    /// X-dropoff for final gapped extensions
    #[arg(long = "xdrop_gap_final", default_value = "100.0")]
    xdrop_gap_final: f64,

    /// Perform ungapped alignment only
    #[arg(long, default_value = "false")]
    ungapped: bool,

    /// Minimum percent identity to report
    #[arg(long = "perc_identity", default_value = "0.0", value_parser = parse_percent_0_to_100)]
    perc_identity: f64,

    /// Minimum query coverage per HSP (percent)
    #[arg(long = "qcov_hsp_perc", default_value = "0.0", value_parser = parse_percent_0_to_100)]
    qcov_hsp_perc: f64,

    /// Maximum number of HSPs per subject
    #[arg(long = "max_hsps", alias = "max-hsps", value_parser = clap::value_parser!(i32).range(1..))]
    max_hsps: Option<i32>,

    /// Culling limit: delete hits enveloped by higher-scoring hits
    #[arg(long = "culling_limit", default_value = "0", value_parser = clap::value_parser!(i32).range(0..))]
    culling_limit: i32,

    /// Window size for multiple hit algorithm
    #[arg(long = "window_size", default_value = "0", value_parser = clap::value_parser!(i32).range(0..))]
    window_size: i32,

    /// Soft masking (filter query but use for lookup)
    #[arg(long = "soft_masking", default_value = "true")]
    soft_masking: String,

    /// Use lowercase masking in query
    #[arg(long = "lcase_masking", default_value = "false")]
    lcase_masking: bool,

    /// Restrict search to query location (format: start-stop)
    #[arg(long = "query_loc", alias = "query-loc")]
    query_loc: Option<String>,

    /// Effective database size
    #[arg(long = "dbsize", default_value = "0")]
    dbsize: i64,

    /// Effective search space
    #[arg(long = "searchsp", default_value = "0", value_parser = clap::value_parser!(i64).range(0..))]
    searchsp: i64,

    /// Task: megablast, blastn, blastn-short, dc-megablast
    #[arg(long = "task", value_parser = ["blastn", "blastn-short", "dc-megablast", "megablast"])]
    task: Option<String>,

    /// Perform no-greedy extension (standard DP instead of greedy)
    #[arg(long = "no_greedy", default_value = "false")]
    no_greedy: bool,

    /// GI list to restrict database sequences
    #[arg(long)]
    gilist: Option<PathBuf>,

    /// SeqID list to restrict database sequences
    #[arg(long)]
    seqidlist: Option<PathBuf>,

    /// Negative GI list (exclude these)
    #[arg(long = "negative_gilist", alias = "negative-gilist")]
    negative_gilist: Option<PathBuf>,

    /// Negative SeqID list (exclude these)
    #[arg(long = "negative_seqidlist", alias = "negative-seqidlist")]
    negative_seqidlist: Option<PathBuf>,

    /// TaxID list to restrict database sequences
    #[arg(long)]
    taxidlist: Option<PathBuf>,

    /// TaxIDs to restrict (comma-separated)
    #[arg(long)]
    taxids: Option<String>,

    /// Negative TaxID list (exclude these)
    #[arg(long = "negative_taxidlist", alias = "negative-taxidlist")]
    negative_taxidlist: Option<PathBuf>,

    /// Negative TaxIDs (comma-separated)
    #[arg(long = "negative_taxids", alias = "negative-taxids")]
    negative_taxids: Option<String>,

    /// Import search strategy from file
    #[arg(long = "import_search_strategy", alias = "import-search-strategy")]
    import_search_strategy: Option<PathBuf>,

    /// Export search strategy to file
    #[arg(long = "export_search_strategy", alias = "export-search-strategy")]
    export_search_strategy: Option<PathBuf>,

    /// Show GIs in deflines
    #[arg(long = "show_gis", alias = "show-gis", default_value = "false")]
    show_gis: bool,

    /// Number of descriptions to show (pairwise output)
    #[arg(long = "num_descriptions", alias = "num-descriptions", value_parser = clap::value_parser!(i32).range(0..))]
    num_descriptions: Option<i32>,

    /// Number of alignments to show (pairwise output)
    #[arg(long = "num_alignments", alias = "num-alignments", value_parser = clap::value_parser!(i32).range(0..))]
    num_alignments: Option<i32>,

    /// Line length for pairwise output
    #[arg(long = "line_length", alias = "line-length", value_parser = clap::value_parser!(i32).range(1..))]
    line_length: Option<i32>,

    /// Produce HTML output
    #[arg(long, default_value = "false")]
    html: bool,

    /// Sort hits: 0=evalue, 1=bitscore, 2=total_score, 3=pct_identity, 4=query_coverage
    #[arg(long = "sorthits", alias = "sort_hits", default_value = "0", value_parser = clap::value_parser!(i32).range(0..=4))]
    sorthits: i32,

    /// Sort HSPs: 0=evalue, 1=score, 2=query_start, 3=pct_identity, 4=subject_start
    #[arg(long = "sorthsps", alias = "sort_hsps", default_value = "0", value_parser = clap::value_parser!(i32).range(0..=4))]
    sorthsps: i32,

    /// Database for filtering (e.g. repeats)
    #[arg(long)]
    filtering_db: Option<PathBuf>,

    /// WindowMasker database
    #[arg(long)]
    window_masker_db: Option<PathBuf>,

    /// WindowMasker TaxID
    #[arg(long)]
    window_masker_taxid: Option<i32>,

    /// Database soft mask algorithm ID
    #[arg(long)]
    db_soft_mask: Option<i32>,

    /// Database hard mask algorithm ID
    #[arg(long)]
    db_hard_mask: Option<i32>,

    /// Subject best hit
    #[arg(
        long = "subject_besthit",
        alias = "subject-besthit",
        default_value = "false"
    )]
    subject_besthit: bool,

    /// Best hit overhang
    #[arg(long = "best_hit_overhang", alias = "best-hit-overhang", value_parser = parse_best_hit_fraction)]
    best_hit_overhang: Option<f64>,

    /// Best hit score edge
    #[arg(long = "best_hit_score_edge", alias = "best-hit-score-edge", value_parser = parse_best_hit_fraction)]
    best_hit_score_edge: Option<f64>,

    /// Minimum raw gapped score
    #[arg(long = "min_raw_gapped_score", alias = "min-raw-gapped-score")]
    min_raw_gapped_score: Option<i32>,

    /// Template type for discontiguous megablast: coding, optimal, coding_and_optimal
    #[arg(long = "template_type", alias = "template-type", value_parser = ["coding", "optimal", "coding_and_optimal"])]
    template_type: Option<String>,

    /// Template length for discontiguous megablast: 16, 18, 21
    #[arg(long = "template_length", alias = "template-length", value_parser = parse_template_length)]
    template_length: Option<i32>,

    /// Parse deflines in FASTA input
    #[arg(
        long = "parse_deflines",
        alias = "parse-deflines",
        default_value = "false"
    )]
    parse_deflines: bool,

    /// Restrict search to subject location (format: start-stop)
    #[arg(long = "subject_loc", alias = "subject-loc")]
    subject_loc: Option<String>,

    /// Entrez query to restrict database (remote searches)
    #[arg(long = "entrez_query", alias = "entrez-query")]
    entrez_query: Option<String>,

    /// Remote search at NCBI
    #[arg(long, default_value = "false")]
    remote: bool,

    /// Use index for search
    #[arg(long = "use_index", alias = "use-index", default_value = "false")]
    use_index: bool,

    /// Index name
    #[arg(long = "index_name", alias = "index-name")]
    index_name: Option<String>,

    /// Multithreading mode: 0=split by DB, 1=split by query
    #[arg(long = "mt_mode", alias = "mt-mode", default_value = "0", value_parser = clap::value_parser!(i32).range(0..=1))]
    mt_mode: i32,

    /// Off-diagonal range for multi-hit extension
    #[arg(
        long = "off_diagonal_range",
        alias = "off-diagonal-range",
        default_value = "0",
        value_parser = clap::value_parser!(i32).range(0..)
    )]
    off_diagonal_range: i32,

    /// Disable TaxID expansion
    #[arg(
        long = "no_taxid_expansion",
        alias = "no-taxid-expansion",
        default_value = "false"
    )]
    no_taxid_expansion: bool,
}

impl BlastnArgs {
    fn effective_max_target_seqs(&self) -> i32 {
        self.max_target_seqs
            .or(self.num_alignments.filter(|_| outfmt_number(&self.outfmt) > 4))
            .unwrap_or(500)
    }

    /// Apply task-specific defaults for parameters not explicitly set by the user.
    fn apply_task_defaults(&mut self) {
        if let Some(ref task) = self.task {
            match task.as_str() {
                "blastn-short" => {
                    // NCBI defaults for blastn-short
                    if self.word_size.is_none() {
                        self.word_size = Some(7);
                    }
                    if self.reward.is_none() {
                        self.reward = Some(1);
                    }
                    if self.penalty.is_none() {
                        self.penalty = Some(-3);
                    }
                    if self.gapopen.is_none() {
                        self.gapopen = Some(5);
                    }
                    if self.gapextend.is_none() {
                        self.gapextend = Some(2);
                    }
                }
                "blastn" => {
                    if self.word_size.is_none() {
                        self.word_size = Some(11);
                    }
                    if self.reward.is_none() {
                        self.reward = Some(2);
                    }
                    if self.penalty.is_none() {
                        self.penalty = Some(-3);
                    }
                    if self.gapopen.is_none() {
                        self.gapopen = Some(5);
                    }
                    if self.gapextend.is_none() {
                        self.gapextend = Some(2);
                    }
                }
                "megablast" => {
                    if self.word_size.is_none() {
                        self.word_size = Some(28);
                    }
                    if self.reward.is_none() {
                        self.reward = Some(1);
                    }
                    if self.penalty.is_none() {
                        self.penalty = Some(-2);
                    }
                    if self.gapopen.is_none() {
                        self.gapopen = Some(0);
                    }
                    if self.gapextend.is_none() {
                        self.gapextend = Some(0);
                    }
                }
                "dc-megablast" => {
                    if self.word_size.is_none() {
                        self.word_size = Some(28);
                    }
                    if self.reward.is_none() {
                        self.reward = Some(2);
                    }
                    if self.penalty.is_none() {
                        self.penalty = Some(-3);
                    }
                    if self.gapopen.is_none() {
                        self.gapopen = Some(5);
                    }
                    if self.gapextend.is_none() {
                        self.gapextend = Some(2);
                    }
                }
                _ => {}
            }
        }
        // NCBI blastn defaults to the megablast task when -task is omitted.
        if self.word_size.is_none() {
            self.word_size = Some(28);
        }
        if self.reward.is_none() {
            self.reward = Some(1);
        }
        if self.penalty.is_none() {
            self.penalty = Some(-2);
        }
        if self.gapopen.is_none() {
            self.gapopen = Some(0);
        }
        if self.gapextend.is_none() {
            self.gapextend = Some(0);
        }
    }

    fn word_size(&self) -> i32 {
        self.word_size.unwrap_or(11)
    }
    fn reward(&self) -> i32 {
        self.reward.unwrap_or(1)
    }
    fn penalty(&self) -> i32 {
        self.penalty.unwrap_or(-3)
    }
    fn gapopen(&self) -> i32 {
        self.gapopen.unwrap_or(5)
    }
    fn gapextend(&self) -> i32 {
        self.gapextend.unwrap_or(2)
    }
}

fn main() {
    // Run on a thread with 64 MB stack to avoid stack overflow on large databases.
    // The main thread's default 8 MB stack is insufficient for large DB operations
    // (e.g., core_nt with 1.5M OIDs per volume).
    let builder = std::thread::Builder::new().stack_size(256 * 1024 * 1024);
    let handler = builder
        .spawn(main_inner)
        .expect("Failed to spawn main thread");
    if let Err(e) = handler.join() {
        eprintln!("Fatal error: {:?}", e);
        std::process::exit(1);
    }
}

fn main_inner() {
    let cli = Cli::parse();

    let (program, mut args) = match cli.command {
        Commands::Blastn(a) => ("blastn", a),
        Commands::Blastp(a) => ("blastp", a),
        Commands::Blastx(a) => ("blastx", a),
        Commands::Tblastn(a) => ("tblastn", a),
        Commands::Tblastx(a) => ("tblastx", a),
        Commands::Psiblast(a) => ("psiblast", a),
        Commands::Rpsblast(a) => ("rpsblast", a),
        Commands::Rpstblastn(a) => ("rpstblastn", a),
        Commands::Deltablast(a) => ("deltablast", a),
    };

    args.apply_task_defaults();

    // Warn about unsupported options that are silently ignored
    if args.gilist.is_some() {
        eprintln!("Warning: --gilist is not yet supported, ignoring");
    }
    if args.seqidlist.is_some() {
        eprintln!("Warning: --seqidlist is not yet supported, ignoring");
    }
    if args.negative_gilist.is_some() {
        eprintln!("Warning: --negative-gilist is not yet supported, ignoring");
    }
    if args.negative_seqidlist.is_some() {
        eprintln!("Warning: --negative-seqidlist is not yet supported, ignoring");
    }
    if args.remote {
        eprintln!("Warning: --remote is not supported, ignoring");
    }
    emit_sort_option_warnings(program, &args);
    emit_hitlist_size_warnings(program, &args);
    emit_formatting_option_warnings(program, &args);

    if args.db.is_none() && args.subject.is_none() {
        eprintln!("Error: Either --db or --subject is required");
        std::process::exit(1);
    }

    let result = match program {
        "blastn" => run_blastn(&args),
        "blastp" => run_blastp(&args),
        "blastx" => run_blastx(&args),
        "tblastn" => run_tblastn(&args),
        "tblastx" => run_tblastx(&args),
        "psiblast" => run_psiblast(&args),
        "rpsblast" => run_rpsblast(&args),
        "rpstblastn" => run_rpstblastn(&args),
        "deltablast" => run_deltablast(&args),
        _ => unreachable!(),
    };

    if let Err(e) = result {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    }
}

fn outfmt_number(outfmt: &str) -> i32 {
    outfmt
        .split_whitespace()
        .next()
        .and_then(|part| part.parse().ok())
        .unwrap_or(6)
}

fn emit_sort_option_warnings(program: &str, args: &BlastnArgs) {
    let outfmt_num = outfmt_number(&args.outfmt);
    if args.sorthits != 0 && outfmt_num > 4 {
        eprintln!(
            "Warning: [{}] The parameter -sorthits is ignored for output formats > 4.",
            program
        );
    }
    if args.sorthsps != 0 && outfmt_num != 0 {
        eprintln!(
            "Warning: [{}] The parameter -sorthsps is ignored for output formats != 0.",
            program
        );
    }
}

fn emit_hitlist_size_warnings(program: &str, args: &BlastnArgs) {
    if args.effective_max_target_seqs() < 5 {
        eprintln!(
            "Warning: [{}] Examining 5 or more matches is recommended",
            program
        );
    }
}

fn emit_formatting_option_warnings(program: &str, args: &BlastnArgs) {
    let outfmt_num = outfmt_number(&args.outfmt);
    if outfmt_num <= 4 {
        return;
    }
    if args.num_descriptions.is_some() {
        eprintln!(
            "Warning: [{}] The parameter -num_descriptions is ignored for output formats > 4 . Use -max_target_seqs to control output",
            program
        );
    }
    if args.line_length.is_some() {
        eprintln!(
            "Warning: [{}] The parameter -line_length is not applicable for output formats > 4 .",
            program
        );
    }
}

fn run_blastp(args: &BlastnArgs) -> Result<(), Box<dyn std::error::Error>> {
    let query_file = File::open(&args.query)?;
    let queries = parse_fasta(query_file);
    if queries.is_empty() {
        return Err("No sequences found in query file".into());
    }

    if let Some(ref subject_path) = args.subject {
        let subject_file = File::open(subject_path)?;
        let subjects = parse_fasta(subject_file);

        let matrix = blast_rs::matrix::BLOSUM62;

        let prot_kbp = blast_rs::stat::lookup_protein_params(args.gapopen(), args.gapextend())
            .map(|p| blast_rs::stat::KarlinBlk {
                lambda: p.lambda,
                k: p.k,
                log_k: p.k.ln(),
                h: p.h,
            })
            .unwrap_or_else(blast_rs::stat::protein_ungapped_kbp);

        let total_subj_len: usize = subjects.iter().map(|s| s.sequence.len()).sum();
        let avg_q_len =
            queries.iter().map(|q| q.sequence.len()).sum::<usize>() / queries.len().max(1);
        let len_adj = blast_rs::stat::compute_length_adjustment(
            avg_q_len as i32,
            total_subj_len as i64,
            subjects.len() as i32,
            &prot_kbp,
        );
        let search_space = blast_rs::stat::compute_search_space(
            avg_q_len as i64,
            total_subj_len as i64,
            subjects.len() as i32,
            len_adj,
        );

        // Protein default word_size=3 (NCBI default), cap at 6 to keep table size manageable
        let word_size = if args.word_size.is_none() {
            3
        } else {
            args.word_size().max(2).min(6) as usize
        };
        let threshold = 11.0; // NCBI default for word_size=3

        let mut hits = Vec::new();
        for qrec in &queries {
            let query_aa: Vec<u8> = qrec
                .sequence
                .iter()
                .map(|&b| blast_rs::input::aminoacid_to_ncbistdaa(b))
                .collect();

            for srec in &subjects {
                let subj_aa: Vec<u8> = srec
                    .sequence
                    .iter()
                    .map(|&b| blast_rs::input::aminoacid_to_ncbistdaa(b))
                    .collect();

                let phits = blast_rs::protein_lookup::protein_gapped_scan(
                    &query_aa,
                    &subj_aa,
                    &matrix,
                    word_size,
                    threshold,
                    30,
                    args.gapopen(),
                    args.gapextend(),
                    50,
                    0,
                );
                for ph in phits {
                    let evalue = prot_kbp.raw_to_evalue(ph.score, search_space);
                    if evalue > args.evalue {
                        continue;
                    }
                    hits.push(TabularHit {
                        query_id: qrec.id.clone(),
                        subject_id: srec.id.clone(),
                        subject_title: String::new(),
                        pct_identity: if ph.align_length > 0 {
                            100.0 * ph.num_ident as f64 / ph.align_length as f64
                        } else {
                            0.0
                        },
                        align_len: ph.align_length,
                        mismatches: ph.mismatches,
                        gap_opens: ph.gap_opens,
                        query_start: ph.query_start as i32 + 1,
                        query_end: ph.query_end as i32,
                        subject_start: ph.subject_start as i32 + 1,
                        subject_end: ph.subject_end as i32,
                        evalue,
                        bit_score: prot_kbp.raw_to_bit(ph.score),
                        query_len: query_aa.len() as i32,
                        subject_len: subj_aa.len() as i32,
                        raw_score: ph.score,
                        qseq: ph.qseq.clone(),
                        sseq: ph.sseq.clone(),
                        qframe: 0,
                        sframe: 0,
                        subject_taxids: vec![],
                        subject_sci_name: String::new(),
                        subject_common_name: String::new(),
                        subject_blast_name: String::new(),
                        subject_kingdom: String::new(),
                        num_ident: ph.num_ident,
                    });
                }
            }
        }

        hits.sort_by(|a, b| {
            b.bit_score
                .partial_cmp(&a.bit_score)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let mut seen = std::collections::HashSet::new();
        hits.retain(|h| {
            seen.insert((
                h.query_id.clone(),
                h.subject_id.clone(),
                h.query_start,
                h.subject_start,
            ))
        });

        let stdout = io::stdout();
        let mut writer: Box<dyn Write> = if let Some(ref path) = args.out {
            Box::new(BufWriter::new(File::create(path)?))
        } else {
            Box::new(BufWriter::new(stdout.lock()))
        };
        write_tabular_output(&mut writer, &hits, &args.outfmt)?;
        writer.flush()?;
        return Ok(());
    }

    // Database mode
    let db_path = args
        .db
        .as_ref()
        .ok_or("blastp requires --db or --subject")?;
    let db = BlastDb::open(db_path)?;
    if db.db_type != DbType::Protein {
        return Err("blastp requires a protein database".into());
    }

    let mut params = blast_rs::api::SearchParams::blastp()
        .evalue(args.evalue)
        .num_threads(if args.num_threads <= 0 {
            1
        } else {
            args.num_threads as usize
        });
    // Don't override gap costs with blastn defaults (5/2).
    // SearchParams::blastp() already sets the correct protein defaults (11/1).
    // Only override if user explicitly set non-default values on the command line.
    if args.gapopen() != 5 || args.gapextend() != 2 {
        params.gap_open = args.gapopen();
        params.gap_extend = args.gapextend();
    }
    if let Some(ws) = args.word_size {
        params.word_size = ws as usize;
    }
    params.max_target_seqs = args.effective_max_target_seqs() as usize;
    params.max_hsps = args.max_hsps.map(|max| max as usize);
    // TODO: composition adjustment is very slow, disable for now
    params.comp_adjust = 0;

    let mut all_hits = Vec::new();
    for qrec in &queries {
        eprintln!(
            "DEBUG run_blastp: query='{}' len={} gap_open={} gap_extend={}",
            qrec.id,
            qrec.sequence.len(),
            params.gap_open,
            params.gap_extend
        );
        let results = blast_rs::api::blastp(&db, &qrec.sequence, &params);
        eprintln!("DEBUG run_blastp: got {} results", results.len());
        for sr in results {
            for hsp in sr.hsps {
                all_hits.push(TabularHit {
                    query_id: qrec.id.clone(),
                    subject_id: sr.subject_accession.clone(),
                    subject_title: String::new(),
                    pct_identity: hsp.percent_identity(),
                    align_len: hsp.alignment_length as i32,
                    mismatches: (hsp.alignment_length - hsp.num_identities - hsp.num_gaps) as i32,
                    gap_opens: hsp.num_gaps as i32,
                    query_start: hsp.query_start as i32 + 1,
                    query_end: hsp.query_end as i32,
                    subject_start: hsp.subject_start as i32 + 1,
                    subject_end: hsp.subject_end as i32,
                    evalue: hsp.evalue,
                    bit_score: hsp.bit_score,
                    query_len: qrec.sequence.len() as i32,
                    subject_len: sr.subject_len as i32,
                    raw_score: hsp.score,
                    qseq: if hsp.query_aln.is_empty() {
                        None
                    } else {
                        Some(String::from_utf8_lossy(&hsp.query_aln).into_owned())
                    },
                    sseq: if hsp.subject_aln.is_empty() {
                        None
                    } else {
                        Some(String::from_utf8_lossy(&hsp.subject_aln).into_owned())
                    },
                    qframe: 0,
                    sframe: 0,
                    subject_taxids: sr.taxids.clone(),
                    subject_sci_name: String::new(),
                    subject_common_name: String::new(),
                    subject_blast_name: String::new(),
                    subject_kingdom: String::new(),
                    num_ident: hsp.num_identities as i32,
                });
            }
        }
    }

    all_hits.sort_by(|a, b| {
        a.evalue
            .partial_cmp(&b.evalue)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let stdout = io::stdout();
    let mut writer: Box<dyn Write> = if let Some(ref path) = args.out {
        Box::new(BufWriter::new(File::create(path)?))
    } else {
        Box::new(BufWriter::new(stdout.lock()))
    };
    write_tabular_output(&mut writer, &all_hits, &args.outfmt)?;
    writer.flush()?;
    Ok(())
}

fn run_blastx(args: &BlastnArgs) -> Result<(), Box<dyn std::error::Error>> {
    use blast_rs::util::{six_frame_translation, STANDARD_GENETIC_CODE};

    let query_file = File::open(&args.query)?;
    let queries = parse_fasta(query_file);
    if queries.is_empty() {
        return Err("No sequences found in query file".into());
    }

    let subject_path = args
        .subject
        .as_ref()
        .ok_or("blastx requires --subject (protein FASTA)")?;
    let subject_file = File::open(subject_path)?;
    let subjects = parse_fasta(subject_file);

    let matrix = blast_rs::matrix::BLOSUM62;

    let prot_kbp = blast_rs::stat::protein_ungapped_kbp();
    let total_subj_len: usize = subjects.iter().map(|s| s.sequence.len()).sum();

    let mut all_hits = Vec::new();

    for qrec in &queries {
        // Encode nucleotide query
        let nuc_seq: Vec<u8> = qrec
            .sequence
            .iter()
            .map(|&b| match b {
                b'A' | b'a' => 0,
                b'C' | b'c' => 1,
                b'G' | b'g' => 2,
                b'T' | b't' => 3,
                _ => 0,
            })
            .collect();

        // Translate in 6 frames
        let frames = six_frame_translation(&nuc_seq, &STANDARD_GENETIC_CODE);

        for (frame, protein_query) in &frames {
            if protein_query.len() < 3 {
                continue;
            }

            let search_space = blast_rs::stat::compute_search_space(
                protein_query.len() as i64,
                total_subj_len as i64,
                subjects.len() as i32,
                0,
            );

            for srec in &subjects {
                let subj_aa: Vec<u8> = srec
                    .sequence
                    .iter()
                    .map(|&b| blast_rs::input::aminoacid_to_ncbistdaa(b))
                    .collect();

                let phits = blast_rs::protein_lookup::protein_gapped_scan(
                    &protein_query,
                    &subj_aa,
                    &matrix,
                    3,
                    11.0,
                    30,
                    args.gapopen(),
                    args.gapextend(),
                    50,
                    0,
                );
                let nuc_len = nuc_seq.len() as i32;

                for ph in phits {
                    let evalue = prot_kbp.raw_to_evalue(ph.score, search_space);
                    if evalue > args.evalue {
                        continue;
                    }

                    let (q_nuc_start, q_nuc_end) = blast_rs::util::protein_to_nuc_coords(
                        ph.query_start as i32,
                        ph.query_end as i32,
                        *frame,
                        nuc_len,
                    );

                    all_hits.push(TabularHit {
                        query_id: qrec.id.clone(),
                        subject_id: srec.id.clone(),
                        subject_title: String::new(),
                        pct_identity: if ph.align_length > 0 {
                            100.0 * ph.num_ident as f64 / ph.align_length as f64
                        } else {
                            0.0
                        },
                        align_len: ph.align_length,
                        mismatches: ph.mismatches,
                        gap_opens: ph.gap_opens,
                        query_start: q_nuc_start,
                        query_end: q_nuc_end,
                        subject_start: ph.subject_start as i32 + 1,
                        subject_end: ph.subject_end as i32,
                        evalue,
                        bit_score: prot_kbp.raw_to_bit(ph.score),
                        query_len: nuc_len,
                        subject_len: subj_aa.len() as i32,
                        raw_score: ph.score,
                        qseq: ph.qseq.clone(),
                        sseq: ph.sseq.clone(),
                        qframe: *frame,
                        sframe: 0,
                        subject_taxids: vec![],
                        subject_sci_name: String::new(),
                        subject_common_name: String::new(),
                        subject_blast_name: String::new(),
                        subject_kingdom: String::new(),
                        num_ident: ph.num_ident,
                    });
                }
            }
        }
    }

    all_hits.sort_by(|a, b| {
        a.evalue
            .partial_cmp(&b.evalue)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let mut seen = std::collections::HashSet::new();
    all_hits.retain(|h| seen.insert((h.query_id.clone(), h.query_start, h.subject_start)));

    let stdout = io::stdout();
    let mut writer: Box<dyn Write> = if let Some(ref path) = args.out {
        Box::new(BufWriter::new(File::create(path)?))
    } else {
        Box::new(BufWriter::new(stdout.lock()))
    };
    write_tabular_output(&mut writer, &all_hits, &args.outfmt)?;
    writer.flush()?;
    Ok(())
}

fn run_tblastn(args: &BlastnArgs) -> Result<(), Box<dyn std::error::Error>> {
    use blast_rs::util::{six_frame_translation, STANDARD_GENETIC_CODE};

    // tblastn: protein query vs translated nucleotide subject
    let query_file = File::open(&args.query)?;
    let queries = parse_fasta(query_file);
    let subject_path = args
        .subject
        .as_ref()
        .ok_or("tblastn requires --subject (nucleotide FASTA)")?;
    let subject_file = File::open(subject_path)?;
    let subjects = parse_fasta(subject_file);

    let matrix = blast_rs::matrix::BLOSUM62;
    let prot_kbp = blast_rs::stat::protein_ungapped_kbp();

    let mut all_hits = Vec::new();
    for qrec in &queries {
        let query_aa: Vec<u8> = qrec
            .sequence
            .iter()
            .map(|&b| blast_rs::input::aminoacid_to_ncbistdaa(b))
            .collect();

        for srec in &subjects {
            let nuc_seq: Vec<u8> = srec
                .sequence
                .iter()
                .map(|&b| match b {
                    b'A' | b'a' => 0,
                    b'C' | b'c' => 1,
                    b'G' | b'g' => 2,
                    b'T' | b't' => 3,
                    _ => 0,
                })
                .collect();
            let frames = six_frame_translation(&nuc_seq, &STANDARD_GENETIC_CODE);

            for (frame, subj_protein) in &frames {
                if subj_protein.len() < 3 || query_aa.len() < 3 {
                    continue;
                }
                let search_space = (query_aa.len() * subj_protein.len()) as f64;
                let subj_nuc_len = nuc_seq.len() as i32;

                let phits = blast_rs::protein_lookup::protein_gapped_scan(
                    &query_aa,
                    &subj_protein,
                    &matrix,
                    3,
                    11.0,
                    30,
                    args.gapopen(),
                    args.gapextend(),
                    50,
                    0,
                );

                for ph in phits {
                    let evalue = prot_kbp.raw_to_evalue(ph.score, search_space);
                    if evalue > args.evalue {
                        continue;
                    }

                    let (s_nuc_start, s_nuc_end) = blast_rs::util::protein_to_nuc_coords(
                        ph.subject_start as i32,
                        ph.subject_end as i32,
                        *frame,
                        subj_nuc_len,
                    );

                    all_hits.push(TabularHit {
                        query_id: qrec.id.clone(),
                        subject_id: srec.id.clone(),
                        subject_title: String::new(),
                        pct_identity: if ph.align_length > 0 {
                            100.0 * ph.num_ident as f64 / ph.align_length as f64
                        } else {
                            0.0
                        },
                        align_len: ph.align_length,
                        mismatches: ph.mismatches,
                        gap_opens: ph.gap_opens,
                        query_start: ph.query_start as i32 + 1,
                        query_end: ph.query_end as i32,
                        subject_start: s_nuc_start,
                        subject_end: s_nuc_end,
                        evalue,
                        bit_score: prot_kbp.raw_to_bit(ph.score),
                        query_len: query_aa.len() as i32,
                        subject_len: subj_nuc_len,
                        raw_score: ph.score,
                        qseq: ph.qseq.clone(),
                        sseq: ph.sseq.clone(),
                        qframe: 0,
                        sframe: *frame,
                        subject_taxids: vec![],
                        subject_sci_name: String::new(),
                        subject_common_name: String::new(),
                        subject_blast_name: String::new(),
                        subject_kingdom: String::new(),
                        num_ident: ph.num_ident,
                    });
                }
            }
        }
    }
    all_hits.sort_by(|a, b| {
        a.evalue
            .partial_cmp(&b.evalue)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let mut seen = std::collections::HashSet::new();
    all_hits.retain(|h| {
        seen.insert((
            h.query_id.clone(),
            h.subject_id.clone(),
            h.query_start,
            h.subject_start,
        ))
    });

    let stdout = io::stdout();
    let mut writer: Box<dyn Write> = if let Some(ref path) = args.out {
        Box::new(BufWriter::new(File::create(path)?))
    } else {
        Box::new(BufWriter::new(stdout.lock()))
    };
    write_tabular_output(&mut writer, &all_hits, &args.outfmt)?;
    writer.flush()?;
    Ok(())
}

fn run_psiblast(args: &BlastnArgs) -> Result<(), Box<dyn std::error::Error>> {
    use blast_rs::pssm::Pssm;

    let query_file = File::open(&args.query)?;
    let queries = parse_fasta(query_file);
    if queries.is_empty() {
        return Err("No query sequences".into());
    }

    let subject_path = args.subject.as_ref().ok_or("psiblast requires --subject")?;
    let subject_file = File::open(subject_path)?;
    let subjects = parse_fasta(subject_file);

    let matrix = blast_rs::matrix::BLOSUM62;

    let prot_kbp = blast_rs::stat::protein_ungapped_kbp();
    let total_subj_len: usize = subjects.iter().map(|s| s.sequence.len()).sum();
    let search_space = (queries[0].sequence.len() * total_subj_len) as f64;

    // Build initial PSSM from query
    let query_aa: Vec<u8> = queries[0]
        .sequence
        .iter()
        .map(|&b| blast_rs::input::aminoacid_to_ncbistdaa(b))
        .collect();
    let mut pssm = Pssm::from_sequence(&query_aa, &matrix);

    // Prepare subjects as (id, aa_sequence) pairs
    let subj_pairs: Vec<(String, Vec<u8>)> = subjects
        .iter()
        .map(|s| {
            (
                s.id.clone(),
                s.sequence
                    .iter()
                    .map(|&b| blast_rs::input::aminoacid_to_ncbistdaa(b))
                    .collect(),
            )
        })
        .collect();

    // Run 3 iterations
    let mut all_hits = Vec::new();
    for _iter in 0..3 {
        let results = blast_rs::pssm::psi_blast_iteration(
            &pssm,
            &subj_pairs,
            args.evalue,
            search_space,
            prot_kbp.lambda,
            prot_kbp.k,
        );

        if results.is_empty() {
            break;
        }

        all_hits.clear();
        for hit in &results {
            // Extract aligned sequences (ungapped PSSM scan)
            let subj_seq = subj_pairs
                .iter()
                .find(|(sid, _)| sid == &hit.subject_id)
                .map(|(_, s)| s.as_slice());
            let qseq: String = query_aa[..hit.align_len]
                .iter()
                .map(|&b| blast_rs::protein::ncbistdaa_to_char(b))
                .collect();
            let sseq: String = subj_seq
                .map(|s| {
                    s[hit.subject_start..hit.subject_start + hit.align_len]
                        .iter()
                        .map(|&b| blast_rs::protein::ncbistdaa_to_char(b))
                        .collect()
                })
                .unwrap_or_default();
            // Compute identity
            let num_ident = query_aa[..hit.align_len]
                .iter()
                .zip(
                    subj_seq
                        .map(|s| &s[hit.subject_start..hit.subject_start + hit.align_len])
                        .unwrap_or(&[]),
                )
                .filter(|(a, b)| a == b)
                .count() as i32;
            let pct_identity = if hit.align_len > 0 {
                100.0 * num_ident as f64 / hit.align_len as f64
            } else {
                0.0
            };

            all_hits.push(TabularHit {
                query_id: queries[0].id.clone(),
                subject_id: hit.subject_id.clone(),
                subject_title: String::new(),
                pct_identity,
                align_len: hit.align_len as i32,
                mismatches: hit.align_len as i32 - num_ident,
                gap_opens: 0,
                query_start: 1,
                query_end: hit.align_len as i32,
                subject_start: hit.subject_start as i32 + 1,
                subject_end: (hit.subject_start + hit.align_len) as i32,
                evalue: hit.evalue,
                bit_score: prot_kbp.raw_to_bit(hit.score),
                query_len: pssm.length as i32,
                subject_len: hit.subject_len as i32,
                raw_score: hit.score,
                qseq: Some(qseq),
                sseq: Some(sseq),
                qframe: 0,
                sframe: 0,
                subject_taxids: vec![],
                subject_sci_name: String::new(),
                subject_common_name: String::new(),
                subject_blast_name: String::new(),
                subject_kingdom: String::new(),
                num_ident,
            });
        }

        // Update PSSM from aligned sequences (simplified)
        let aligned: Vec<Vec<u8>> = results
            .iter()
            .filter_map(|h| {
                subj_pairs
                    .iter()
                    .find(|(sid, _)| sid == &h.subject_id)
                    .map(|(_, s)| s.clone())
            })
            .collect();
        pssm.update_from_alignment(&aligned, &blast_rs::matrix::AA_FREQUENCIES, 10.0);
    }

    let stdout = io::stdout();
    let mut writer: Box<dyn Write> = if let Some(ref path) = args.out {
        Box::new(BufWriter::new(File::create(path)?))
    } else {
        Box::new(BufWriter::new(stdout.lock()))
    };
    write_tabular_output(&mut writer, &all_hits, &args.outfmt)?;
    writer.flush()?;
    Ok(())
}

fn run_rpsblast(_args: &BlastnArgs) -> Result<(), Box<dyn std::error::Error>> {
    Err("rpsblast requires a pre-built PSSM database (CDD). Use --program blastp for protein search.".into())
}

fn run_rpstblastn(_args: &BlastnArgs) -> Result<(), Box<dyn std::error::Error>> {
    // rpstblastn: translated nucleotide query vs RPS (PSSM) database
    // Like rpsblast but translates the nucleotide query first
    Err("rpstblastn requires a pre-built PSSM database (CDD). Use --program blastx for translated search.".into())
}

fn run_deltablast(args: &BlastnArgs) -> Result<(), Box<dyn std::error::Error>> {
    // deltablast: domain-enhanced lookup time accelerated BLAST
    // Uses CDD domains to construct initial PSSM, then runs PSI-BLAST
    // For --subject mode, fall back to PSI-BLAST behavior
    if args.subject.is_some() {
        eprintln!("Note: deltablast falling back to psiblast (no CDD database)");
        return run_psiblast(args);
    }
    Err("deltablast requires CDD database. Use --program psiblast with --subject for iterative search.".into())
}

fn run_tblastx(args: &BlastnArgs) -> Result<(), Box<dyn std::error::Error>> {
    use blast_rs::util::{six_frame_translation, STANDARD_GENETIC_CODE};

    // tblastx: translated nucleotide query vs translated nucleotide subject
    let query_file = File::open(&args.query)?;
    let queries = parse_fasta(query_file);
    let subject_path = args
        .subject
        .as_ref()
        .ok_or("tblastx requires --subject (nucleotide FASTA)")?;
    let subject_file = File::open(subject_path)?;
    let subjects = parse_fasta(subject_file);

    let matrix = blast_rs::matrix::BLOSUM62;
    let prot_kbp = blast_rs::stat::protein_ungapped_kbp();

    let mut all_hits = Vec::new();
    for qrec in &queries {
        let q_nuc: Vec<u8> = qrec
            .sequence
            .iter()
            .map(|&b| match b {
                b'A' | b'a' => 0,
                b'C' | b'c' => 1,
                b'G' | b'g' => 2,
                b'T' | b't' => 3,
                _ => 0,
            })
            .collect();
        let q_frames = six_frame_translation(&q_nuc, &STANDARD_GENETIC_CODE);

        for srec in &subjects {
            let s_nuc: Vec<u8> = srec
                .sequence
                .iter()
                .map(|&b| match b {
                    b'A' | b'a' => 0,
                    b'C' | b'c' => 1,
                    b'G' | b'g' => 2,
                    b'T' | b't' => 3,
                    _ => 0,
                })
                .collect();
            let s_frames = six_frame_translation(&s_nuc, &STANDARD_GENETIC_CODE);

            for (qframe, q_prot) in &q_frames {
                for (sframe, s_prot) in &s_frames {
                    if q_prot.len() < 3 || s_prot.len() < 3 {
                        continue;
                    }
                    let search_space = (q_prot.len() * s_prot.len()) as f64;

                    let phits = blast_rs::protein_lookup::protein_gapped_scan(
                        q_prot,
                        s_prot,
                        &matrix,
                        3,
                        13.0,
                        25,
                        args.gapopen(),
                        args.gapextend(),
                        50,
                        0,
                    );

                    for ph in phits {
                        let evalue = prot_kbp.raw_to_evalue(ph.score, search_space);
                        if evalue > args.evalue {
                            continue;
                        }
                        all_hits.push(TabularHit {
                            query_id: qrec.id.clone(),
                            subject_id: srec.id.clone(),
                            subject_title: String::new(),
                            pct_identity: if ph.align_length > 0 {
                                100.0 * ph.num_ident as f64 / ph.align_length as f64
                            } else {
                                0.0
                            },
                            align_len: ph.align_length,
                            mismatches: ph.mismatches,
                            gap_opens: ph.gap_opens,
                            query_start: ph.query_start as i32 + 1,
                            query_end: ph.query_end as i32,
                            subject_start: ph.subject_start as i32 + 1,
                            subject_end: ph.subject_end as i32,
                            evalue,
                            bit_score: prot_kbp.raw_to_bit(ph.score),
                            query_len: q_prot.len() as i32,
                            subject_len: s_prot.len() as i32,
                            raw_score: ph.score,
                            qseq: ph.qseq.clone(),
                            sseq: ph.sseq.clone(),
                            qframe: *qframe,
                            sframe: *sframe,
                            subject_taxids: vec![],
                            subject_sci_name: String::new(),
                            subject_common_name: String::new(),
                            subject_blast_name: String::new(),
                            subject_kingdom: String::new(),
                            num_ident: ph.num_ident,
                        });
                    }
                }
            }
        }
    }
    all_hits.sort_by(|a, b| {
        a.evalue
            .partial_cmp(&b.evalue)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let mut seen = std::collections::HashSet::new();
    all_hits.retain(|h| {
        seen.insert((
            h.query_id.clone(),
            h.subject_id.clone(),
            h.query_start,
            h.subject_start,
        ))
    });

    let stdout = io::stdout();
    let mut writer: Box<dyn Write> = if let Some(ref path) = args.out {
        Box::new(BufWriter::new(File::create(path)?))
    } else {
        Box::new(BufWriter::new(stdout.lock()))
    };
    write_tabular_output(&mut writer, &all_hits, &args.outfmt)?;
    writer.flush()?;
    Ok(())
}

fn run_blastn(args: &BlastnArgs) -> Result<(), Box<dyn std::error::Error>> {
    let query_file = File::open(&args.query)?;
    let records = parse_fasta(query_file);
    if records.is_empty() {
        return Err("No query sequences found".into());
    }

    // Subject mode (FASTA vs FASTA)
    if let Some(ref subject_path) = args.subject {
        let subject_records = parse_fasta(File::open(subject_path)?);
        return run_blastn_subject(args, &records, &subject_records);
    }

    let db_path = args
        .db
        .as_ref()
        .ok_or("Either --db or --subject is required")?;
    #[cfg(not(test))]
    let mut db = BlastDb::open(db_path)?;
    #[cfg(test)]
    let db = BlastDb::open(db_path)?;
    if db.db_type != DbType::Nucleotide {
        return Err("blastn requires a nucleotide database".into());
    }
    #[cfg(not(test))]
    if outfmt_requests_taxonomy(&args.outfmt) {
        db.load_tax_lookup_from_base_path(db_path)?;
    }

    run_blastn_rust(args, &records, db)
}

#[cfg_attr(test, allow(dead_code))]
fn outfmt_requests_taxonomy(outfmt: &str) -> bool {
    let mut parts = outfmt.split_whitespace();
    if parts
        .next()
        .and_then(|s| s.parse::<i32>().ok())
        .unwrap_or(6)
        != 6
    {
        return false;
    }
    parts.any(|field| {
        matches!(
            field,
            "staxid"
                | "staxids"
                | "ssciname"
                | "scomname"
                | "sblastname"
                | "sskingdom"
                | "sskingdoms"
        )
    })
}

/// Pure Rust blastn search — no FFI calls.
fn run_blastn_rust(
    args: &BlastnArgs,
    records: &[blast_rs::input::FastaRecord],
    db: BlastDb,
) -> Result<(), Box<dyn std::error::Error>> {
    use rayon::prelude::*;

    // NCBI BLAST resolves taxonomy names through taxdb on BLASTDB. A taxdb
    // sitting next to an explicitly path-qualified database is not enough.
    let tax_name_db = std::env::var_os("BLASTDB").and_then(|paths| {
        std::env::split_paths(&paths).find_map(|path| blast_rs::db::TaxNameDb::open(&path).ok())
    });

    // Pure Rust KBP computation — per-context (plus + minus) for exact e-value matching
    // The C engine computes separate KBP for each context based on strand composition
    let (kbp_plus, kbp_minus, searchsp_plus, searchsp_minus) = {
        let rec = &records[0];
        let (loc_start, loc_end) = query_loc_bounds(args, rec.sequence.len())?;
        let query_plus: Vec<u8> = rec.sequence[loc_start..loc_end]
            .iter()
            .map(|&b| iupacna_to_blastna(b))
            .collect();
        let query_minus: Vec<u8> = query_plus
            .iter()
            .rev()
            .map(|&b| complement_blastna(b))
            .collect();

        const BLASTNA_SIZE: usize = 16;
        const BLASTNA_TO_NCBI4NA: [u8; BLASTNA_SIZE] =
            [1, 2, 4, 8, 5, 10, 3, 12, 9, 6, 14, 13, 11, 7, 15, 0];
        let reward = args.reward();
        let penalty = args.penalty();
        let matrix_fn = |i: usize, j: usize| -> i32 {
            if i >= BLASTNA_SIZE || j >= BLASTNA_SIZE {
                return penalty;
            }
            if i == BLASTNA_SIZE - 1 || j == BLASTNA_SIZE - 1 {
                return i32::MIN / 2;
            }
            if BLASTNA_TO_NCBI4NA[i] & BLASTNA_TO_NCBI4NA[j] != 0 {
                let d = (0..4)
                    .filter(|&k| BLASTNA_TO_NCBI4NA[j] & BLASTNA_TO_NCBI4NA[k] != 0)
                    .count() as i32;
                (((d - 1) as f64 * penalty as f64 + reward as f64) / d as f64).round() as i32
            } else {
                penalty
            }
        };
        let mut lo = i32::MAX;
        let mut hi = i32::MIN;
        for i in 0..BLASTNA_SIZE {
            for j in 0..BLASTNA_SIZE {
                let s = matrix_fn(i, j);
                if s <= -100000000 || s >= 100000000 {
                    continue;
                }
                if s < lo {
                    lo = s;
                }
                if s > hi {
                    hi = s;
                }
            }
        }

        let ambig: &[u8] = &[14, 15];
        let qlen = query_plus.len() as i32;

        // Compute KBP for both contexts
        let ctxs = [blast_rs::stat::UngappedKbpContext {
            query_offset: 0,
            query_length: qlen,
            is_valid: true,
        }];
        let plus_kbp_results = blast_rs::stat::ungapped_kbp_calc(
            &query_plus,
            &ctxs,
            lo,
            hi,
            BLASTNA_SIZE,
            ambig,
            &matrix_fn,
        );
        let minus_kbp_results = blast_rs::stat::ungapped_kbp_calc(
            &query_minus,
            &ctxs,
            lo,
            hi,
            BLASTNA_SIZE,
            ambig,
            &matrix_fn,
        );

        let default_kbp = blast_rs::stat::KarlinBlk {
            lambda: 1.374,
            k: 0.621,
            log_k: 0.621_f64.ln(),
            h: 1.286,
        };
        let ungapped_plus = plus_kbp_results[0].clone().unwrap_or(default_kbp.clone());
        let ungapped_minus = minus_kbp_results[0].clone().unwrap_or(default_kbp.clone());

        // Gapped KBP (from table — may fall back to ungapped for large gap costs)
        let (gkbp_plus, _) = blast_rs::stat::nucl_gapped_kbp_lookup(
            args.gapopen(),
            args.gapextend(),
            reward,
            penalty,
            &ungapped_plus,
        )
        .unwrap_or((ungapped_plus.clone(), false));
        let (gkbp_minus, _) = blast_rs::stat::nucl_gapped_kbp_lookup(
            args.gapopen(),
            args.gapextend(),
            reward,
            penalty,
            &ungapped_minus,
        )
        .unwrap_or((ungapped_minus.clone(), false));

        // Per-context search space
        let compute_searchsp =
            |kbp: &blast_rs::stat::KarlinBlk, ukbp: &blast_rs::stat::KarlinBlk| -> f64 {
                let (alpha, beta) = blast_rs::stat::nucl_alpha_beta(
                    reward,
                    penalty,
                    args.gapopen(),
                    args.gapextend(),
                    ukbp.lambda,
                    ukbp.h,
                    true,
                );
                let (len_adj, _) = blast_rs::stat::compute_length_adjustment_exact(
                    kbp.k,
                    kbp.log_k,
                    alpha / kbp.lambda,
                    beta,
                    qlen,
                    db.total_length as i64,
                    db.num_oids as i32,
                );
                let eff_db = std::cmp::max(
                    db.total_length as i64 - db.num_oids as i64 * len_adj as i64,
                    1,
                );
                eff_db as f64 * (qlen - len_adj).max(1) as f64
            };
        let sp_plus = compute_searchsp(&gkbp_plus, &ungapped_plus);
        let sp_minus = compute_searchsp(&gkbp_minus, &ungapped_minus);

        (gkbp_plus, gkbp_minus, sp_plus, sp_minus)
    };
    // Use plus-strand KBP as the primary (for seed finding threshold)
    let kbp = &kbp_plus;
    let search_space = searchsp_plus;
    let word_size = args.word_size() as usize;
    #[cfg(not(test))]
    let use_contiguous_megablast_lookup = args
        .task
        .as_deref()
        .map_or(true, |task| task == "megablast");
    let reward = args.reward();
    let penalty = args.penalty();
    let gapopen = args.gapopen();
    let gapextend = args.gapextend();
    let evalue = args.evalue;

    let mut all_hits: Vec<(u32, TabularHit)> = Vec::new();

    let apply_dust = args.dust != "no";

    // Pre-compute invariant values outside the per-query loop
    let cutoff_score = {
        let e = evalue.max(1.0e-297);
        ((kbp.k * search_space / e).ln() / kbp.lambda).ceil() as i32
    };
    let gapped_x_dropoff = (args.xdrop_gap_final * std::f64::consts::LN_2 / kbp.lambda) as i32;
    let num_threads = if args.num_threads <= 0 {
        rayon::current_num_threads()
    } else {
        args.num_threads as usize
    };
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .stack_size(64 * 1024 * 1024) // 64 MB per thread to avoid stack overflow on large DBs
        .build()
        .unwrap_or_else(|_| rayon::ThreadPoolBuilder::new().build().unwrap());

    let search_plus = args.strand != "minus";
    let search_minus = args.strand != "plus";

    // Pre-encode all queries (masked + unmasked, plus + minus strands)
    struct EncodedQuery {
        id: String,
        seq_len: i32,
        query_offset: i32,
        plus_masked: Vec<u8>,
        minus_masked: Vec<u8>,
        plus_nomask: Vec<u8>,
        minus_nomask: Vec<u8>,
    }
    let encoded_queries: Vec<EncodedQuery> = records
        .iter()
        .map(|rec| -> Result<EncodedQuery, Box<dyn std::error::Error>> {
            let (loc_start, loc_end) = query_loc_bounds(args, rec.sequence.len())?;
            let raw_query = &rec.sequence[loc_start..loc_end];
            let plus_nomask: Vec<u8> = raw_query.iter().map(|&b| iupacna_to_blastna(b)).collect();
            let minus_nomask: Vec<u8> = plus_nomask
                .iter()
                .rev()
                .map(|&b| complement_blastna(b))
                .collect();
            let mut plus_masked = plus_nomask.clone();
            if apply_dust {
                apply_blastn_dust_mask(&mut plus_masked);
            }
            if args.lcase_masking {
                apply_lowercase_mask(raw_query, &mut plus_masked);
            }
            let minus_masked: Vec<u8> = plus_masked
                .iter()
                .rev()
                .map(|&b| complement_blastna(b))
                .collect();
            Ok(EncodedQuery {
                id: rec.id.clone(),
                seq_len: plus_nomask.len() as i32,
                query_offset: loc_start as i32,
                plus_masked: if search_plus { plus_masked } else { Vec::new() },
                minus_masked: if search_minus {
                    minus_masked
                } else {
                    Vec::new()
                },
                plus_nomask: if search_plus { plus_nomask } else { Vec::new() },
                minus_nomask: if search_minus {
                    minus_nomask
                } else {
                    Vec::new()
                },
            })
        })
        .collect::<Result<_, _>>()?;

    // Scan subjects in parallel. For each subject, check ALL queries.
    // This keeps subject data in cache across queries (subject-major order).
    let hitlist_size = args.effective_max_target_seqs();
    let prelim_hitlist_size =
        std::cmp::min(std::cmp::max(2 * hitlist_size, 10), hitlist_size + 50) as usize;

    // Collect hits: (query_idx, oid, hsps)
    #[cfg(not(test))]
    let per_subject_hits: Vec<Vec<(usize, u32, Vec<blast_rs::search::SearchHsp>)>> = if num_threads
        <= 1
    {
        let prepared_queries: Vec<_> = encoded_queries
            .iter()
            .map(|eq| {
                if use_contiguous_megablast_lookup {
                    blast_rs::search::PreparedBlastnQuery::new_megablast(
                        &eq.plus_masked,
                        &eq.minus_masked,
                        word_size,
                    )
                } else {
                    blast_rs::search::PreparedBlastnQuery::new(
                        &eq.plus_masked,
                        &eq.minus_masked,
                        word_size,
                    )
                }
            })
            .collect();
        let mut scratch: Vec<_> = prepared_queries
            .iter()
            .map(|prepared| prepared.last_hit_scratch())
            .collect();

        let mut collected = Vec::new();
        for (volume_idx, (start_oid, end_oid)) in db.volume_oid_ranges().into_iter().enumerate() {
            db.advise_volume_sequential(volume_idx);
            for oid in start_oid..end_oid {
                let local_oid = oid - start_oid;
                let (packed, seq_len) = db.get_volume_sequence_and_len(volume_idx, local_oid);
                let seq_len = seq_len as usize;
                let mut results = Vec::new();
                for (qi, eq) in encoded_queries.iter().enumerate() {
                    let mut hsps = blast_rs::search::blastn_gapped_search_packed_prepared(
                        &prepared_queries[qi],
                        &eq.plus_nomask,
                        &eq.minus_nomask,
                        packed,
                        seq_len,
                        reward,
                        penalty,
                        gapopen,
                        gapextend,
                        gapped_x_dropoff,
                        &kbp,
                        search_space,
                        evalue,
                        &mut scratch[qi],
                    );
                    hsps.retain(|h| h.score >= cutoff_score);
                    if !hsps.is_empty() {
                        results.push((qi, oid, hsps));
                    }
                }
                if !results.is_empty() {
                    collected.push(results);
                }
            }
            db.advise_volume_dontneed(volume_idx);
        }
        collected
    } else {
        let prepared_queries: Vec<_> = encoded_queries
            .iter()
            .map(|eq| {
                if use_contiguous_megablast_lookup {
                    blast_rs::search::PreparedBlastnQuery::new_megablast(
                        &eq.plus_masked,
                        &eq.minus_masked,
                        word_size,
                    )
                } else {
                    blast_rs::search::PreparedBlastnQuery::new(
                        &eq.plus_masked,
                        &eq.minus_masked,
                        word_size,
                    )
                }
            })
            .collect();

        let mut collected = Vec::new();
        for (volume_idx, (start_oid, end_oid)) in db.volume_oid_ranges().into_iter().enumerate() {
            db.advise_volume_sequential(volume_idx);
            let mut volume_hits: Vec<_> = pool.install(|| {
                (start_oid..end_oid)
                    .into_par_iter()
                    .map_init(
                        || {
                            prepared_queries
                                .iter()
                                .map(|prepared| prepared.last_hit_scratch())
                                .collect::<Vec<_>>()
                        },
                        |scratch, oid| {
                            let local_oid = oid - start_oid;
                            let (packed, seq_len) =
                                db.get_volume_sequence_and_len(volume_idx, local_oid);
                            let seq_len = seq_len as usize;
                            let mut results = Vec::new();
                            for (qi, eq) in encoded_queries.iter().enumerate() {
                                let mut hsps =
                                    blast_rs::search::blastn_gapped_search_packed_prepared(
                                        &prepared_queries[qi],
                                        &eq.plus_nomask,
                                        &eq.minus_nomask,
                                        packed,
                                        seq_len,
                                        reward,
                                        penalty,
                                        gapopen,
                                        gapextend,
                                        gapped_x_dropoff,
                                        &kbp,
                                        search_space,
                                        evalue,
                                        &mut scratch[qi],
                                    );
                                hsps.retain(|h| h.score >= cutoff_score);
                                if !hsps.is_empty() {
                                    results.push((qi, oid, hsps));
                                }
                            }
                            if results.is_empty() {
                                None
                            } else {
                                Some(results)
                            }
                        },
                    )
                    .filter_map(|hits| hits)
                    .collect()
            });
            collected.append(&mut volume_hits);
            db.advise_volume_dontneed(volume_idx);
        }
        collected
    };

    #[cfg(test)]
    let per_subject_hits: Vec<Vec<(usize, u32, Vec<blast_rs::search::SearchHsp>)>> =
        pool.install(|| {
            (0..db.num_oids)
                .into_par_iter()
                .filter_map(|oid| {
                    let packed = db.get_sequence(oid);
                    let seq_len = db.get_seq_len(oid) as usize;
                    let mut results = Vec::new();
                    for (qi, eq) in encoded_queries.iter().enumerate() {
                        let mut hsps = blast_rs::search::blastn_gapped_search_packed(
                            &eq.plus_masked,
                            &eq.minus_masked,
                            &eq.plus_nomask,
                            &eq.minus_nomask,
                            packed,
                            seq_len,
                            word_size,
                            reward,
                            penalty,
                            gapopen,
                            gapextend,
                            gapped_x_dropoff,
                            &kbp,
                            search_space,
                            evalue,
                        );
                        hsps.retain(|h| h.score >= cutoff_score);
                        if !hsps.is_empty() {
                            results.push((qi, oid, hsps));
                        }
                    }
                    if results.is_empty() {
                        None
                    } else {
                        Some(results)
                    }
                })
                .collect()
        });

    // Flatten and group by query
    let mut per_query_oid_hits: Vec<Vec<(u32, Vec<blast_rs::search::SearchHsp>)>> =
        vec![Vec::new(); encoded_queries.len()];
    for subject_results in per_subject_hits {
        for (qi, oid, hsps) in subject_results {
            per_query_oid_hits[qi].push((oid, hsps));
        }
    }

    // Process each query's results
    for (qi, mut oid_hits) in per_query_oid_hits.into_iter().enumerate() {
        let eq = &encoded_queries[qi];

        // Sort by OID for deterministic output
        oid_hits.sort_by_key(|&(oid, _)| oid);

        // Hitlist pruning
        if oid_hits.len() > prelim_hitlist_size {
            oid_hits.sort_by(|(a_oid, a_hsps), (b_oid, b_hsps)| {
                compare_oid_hsps_for_hitlist(*a_oid, a_hsps, *b_oid, b_hsps)
            });
            oid_hits.truncate(prelim_hitlist_size);
            oid_hits.sort_by_key(|&(oid, _)| oid);
        }

        for (oid, hsps) in oid_hits {
            for hsp in hsps {
                let subject_id = db
                    .get_accession(oid)
                    .unwrap_or_else(|| format!("gnl|BL_ORD_ID|{}", oid));
                let query_len = eq.seq_len;
                let query_offset = eq.query_offset;

                let (q_start, q_end) = if hsp.context == 1 {
                    (
                        query_offset + query_len - hsp.query_end + 1,
                        query_offset + query_len - hsp.query_start,
                    )
                } else {
                    (
                        query_offset + hsp.query_start + 1,
                        query_offset + hsp.query_end,
                    )
                };
                let (s_start, s_end) = if hsp.context == 1 {
                    (hsp.subject_end, hsp.subject_start + 1)
                } else {
                    (hsp.subject_start + 1, hsp.subject_end)
                };

                all_hits.push((
                    oid,
                    TabularHit {
                        query_id: eq.id.clone(),
                        subject_id,
                        subject_title: String::new(),
                        pct_identity: if hsp.align_length > 0 {
                            100.0 * hsp.num_ident as f64 / hsp.align_length as f64
                        } else {
                            0.0
                        },
                        align_len: hsp.align_length,
                        mismatches: hsp.mismatches,
                        gap_opens: hsp.gap_opens,
                        query_start: q_start,
                        query_end: q_end,
                        subject_start: s_start,
                        subject_end: s_end,
                        evalue: if hsp.context == 1 {
                            kbp_minus.raw_to_evalue(hsp.score, searchsp_minus)
                        } else {
                            kbp_plus.raw_to_evalue(hsp.score, searchsp_plus)
                        },
                        bit_score: if hsp.context == 1 {
                            kbp_minus.raw_to_bit(hsp.score)
                        } else {
                            kbp_plus.raw_to_bit(hsp.score)
                        },
                        query_len,
                        subject_len: db.get_seq_len(oid) as i32,
                        raw_score: hsp.score,
                        qseq: hsp.qseq.clone(),
                        sseq: hsp.sseq.clone(),
                        qframe: 1,
                        sframe: if hsp.context == 1 { -1 } else { 1 },
                        subject_taxids: db.get_taxids(oid),
                        subject_sci_name: String::new(),
                        subject_common_name: String::new(),
                        subject_blast_name: String::new(),
                        subject_kingdom: String::new(),
                        num_ident: hsp.num_ident,
                    },
                ));
            }
        }
    }

    // Sort to match C reference output order.
    // The C engine's hitlist preserves the order subjects were inserted during
    // traceback. The traceback processes HSPLists from sorted_hsplists
    // (sorted by descending OID, read from the end = ascending OID order).
    // So hitlist order = ascending OID = insertion order.
    // The output then iterates the hitlist sequentially.
    // Within a subject, HSPs are sorted by score descending.
    //
    // But with max_target_seqs limiting: the hitlist uses a heap that keeps
    // the best subjects. When the heap is full, new subjects replace the worst.
    // The final order after heap extraction is NOT simple ascending OID —
    // it depends on insertion/eviction patterns.
    //
    // For exact matching, we'd need to replicate the heap behavior. For now,
    // sort by: best score descending (best subject first), then ascending OID,
    // then within subject by score descending.
    // Sort to match C engine's Blast_HSPResultsSortByEvalue:
    // Group hits by OID, sort groups by (ascending e-value, descending score, descending OID).
    // Within each group: sort HSPs by descending score.
    // Use C's libc qsort for exact platform-matching tie-breaking behavior.
    {
        // Build per-OID groups
        let mut groups: std::collections::BTreeMap<u32, Vec<TabularHit>> =
            std::collections::BTreeMap::new();
        let mut best_evalue: std::collections::HashMap<u32, f64> = std::collections::HashMap::new();
        let mut best_score: std::collections::HashMap<u32, i32> = std::collections::HashMap::new();
        for (oid, hit) in all_hits {
            let ee = best_evalue.entry(oid).or_insert(hit.evalue);
            if hit.evalue < *ee {
                *ee = hit.evalue;
            }
            let es = best_score
                .entry(oid)
                .or_insert((hit.bit_score * 1000.0) as i32);
            let raw = (hit.bit_score * 1000.0) as i32;
            if raw > *es {
                *es = raw;
            }
            groups.entry(oid).or_default().push(hit);
        }

        // Sort HSPs within each group by descending score
        for hsps in groups.values_mut() {
            hsps.sort_by(|a, b| {
                b.bit_score
                    .partial_cmp(&a.bit_score)
                    .unwrap_or(std::cmp::Ordering::Equal)
                    .then(a.query_start.cmp(&b.query_start))
            });
        }

        // Sort groups using C's qsort comparator logic with libc qsort for exact match
        let mut sorted_oids: Vec<u32> = groups.keys().copied().collect();

        // Use libc qsort to match C's exact tie-breaking behavior
        sorted_oids.sort_by(|&a, &b| {
            let ev_a = best_evalue[&a];
            let ev_b = best_evalue[&b];
            ev_a.partial_cmp(&ev_b)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| best_score[&b].cmp(&best_score[&a]))
                .then_with(|| b.cmp(&a))
        });

        // Reconstruct all_hits in sorted order
        all_hits = Vec::new();
        for oid in sorted_oids {
            if let Some(hsps) = groups.remove(&oid) {
                for hit in hsps {
                    all_hits.push((oid, hit));
                }
            }
        }
    }
    // Filter by e-value using per-context KBP (hits may have been admitted by
    // the search using plus-strand KBP but have worse e-value with their actual context KBP)
    all_hits.retain(|(_, hit)| hit.evalue <= evalue);

    // Apply max_target_seqs per query (matches C engine behavior)
    let max_subjects = args.effective_max_target_seqs() as usize;
    let mut seen_per_query: std::collections::HashMap<String, std::collections::HashSet<String>> =
        std::collections::HashMap::new();
    all_hits.retain(|(_, hit)| {
        let seen = seen_per_query.entry(hit.query_id.clone()).or_default();
        if seen.len() >= max_subjects && !seen.contains(&hit.subject_id) {
            false
        } else {
            seen.insert(hit.subject_id.clone());
            true
        }
    });

    // Apply post-search filters
    let query_len = records
        .first()
        .map(|r| r.sequence.len() as i32)
        .unwrap_or(0);
    let subject_deflines: std::collections::HashMap<String, String> = all_hits
        .iter()
        .filter_map(|(oid, hit)| {
            db_pairwise_subject_defline(&db, *oid, &hit.subject_id)
                .map(|defline| (hit.subject_id.clone(), defline))
        })
        .collect();
    for (oid, hit) in &mut all_hits {
        if let Some(title) = db_subject_title(&db, *oid, &hit.subject_id) {
            hit.subject_title = title;
        }
    }
    let db_sam_labels: std::collections::HashMap<_, _> = all_hits
        .iter()
        .map(|(oid, hit)| (sam_hit_key(hit), format!("BL_ORD_ID:{}", oid)))
        .collect();
    let db_xml_groups: std::collections::HashMap<_, _> = all_hits
        .iter()
        .map(|(oid, hit)| (sam_hit_key(hit), format!("BL_ORD_ID:{}", oid)))
        .collect();
    let db_xml_metadata: std::collections::HashMap<String, (String, String, String, i32)> =
        all_hits
            .iter()
            .filter_map(|(oid, hit)| {
                let group = format!("BL_ORD_ID:{}", oid);
                let defline = db_subject_title(&db, *oid, &hit.subject_id)
                    .or_else(|| db_subject_defline(&db, *oid, &hit.subject_id))?;
                Some((
                    group,
                    (
                        format!("gnl|BL_ORD_ID|{}", oid),
                        defline,
                        oid.to_string(),
                        hit.subject_len,
                    ),
                ))
            })
            .collect();
    let mut hits: Vec<TabularHit> = all_hits.into_iter().map(|(_, h)| h).collect();
    apply_filters(&mut hits, args, query_len);
    let db_sam_labels: std::collections::HashMap<String, String> = hits
        .iter()
        .filter_map(|hit| {
            db_sam_labels
                .get(&sam_hit_key(hit))
                .map(|label| (sam_hit_label_key(hit), label.clone()))
        })
        .collect();
    let db_xml_groups: std::collections::HashMap<String, String> = hits
        .iter()
        .filter_map(|hit| {
            db_xml_groups
                .get(&sam_hit_key(hit))
                .map(|label| (sam_hit_label_key(hit), label.clone()))
        })
        .collect();

    // Enrich hits with taxonomy names from taxdb (if available)
    if let Some(ref tdb) = tax_name_db {
        for hit in &mut hits {
            if let Some(&taxid) = hit.subject_taxids.first() {
                if let Some(info) = tdb.get_info(taxid) {
                    hit.subject_sci_name = info.scientific_name;
                    hit.subject_common_name = info.common_name;
                    hit.subject_blast_name = info.blast_name;
                    hit.subject_kingdom = if info.kingdom == "-" {
                        String::new()
                    } else {
                        info.kingdom
                    };
                }
            }
        }
    }

    // Output
    let stdout = io::stdout();
    let mut writer: Box<dyn Write> = if let Some(ref path) = args.out {
        Box::new(BufWriter::new(File::create(path)?))
    } else {
        Box::new(BufWriter::new(stdout.lock()))
    };
    let outfmt_parts: Vec<&str> = args.outfmt.split_whitespace().collect();
    let outfmt_num: i32 = outfmt_parts[0].parse().unwrap_or(6);

    if outfmt_num == 0 {
        write_pairwise_db_report_preamble(&mut writer, &db)?;
        for query in records {
            let query_hits = pairwise_query_hits(&hits, &query.id, args.sorthits, args.sorthsps);
            write_pairwise_db_query_header(
                &mut writer,
                query,
                &query_hits,
                &subject_deflines,
                args,
            )?;
            let mut last_pairwise_subject: Option<&str> = None;
            for hit in query_hits {
                let mut query_aln = hit
                    .qseq
                    .as_deref()
                    .map(alignment_string_to_blastna)
                    .unwrap_or_default();
                let mut subject_aln = hit
                    .sseq
                    .as_deref()
                    .map(alignment_string_to_blastna)
                    .unwrap_or_default();
                if hit.subject_start > hit.subject_end {
                    query_aln = reverse_complement_alignment_blastna(&query_aln);
                    subject_aln = reverse_complement_alignment_blastna(&subject_aln);
                }
                let has_previous_subject = last_pairwise_subject.is_some();
                let show_subject_header = last_pairwise_subject != Some(hit.subject_id.as_str());
                last_pairwise_subject = Some(hit.subject_id.as_str());
                if show_subject_header && has_previous_subject {
                    writeln!(writer)?;
                }
                write_pairwise_alignment(
                    &mut writer,
                    hit,
                    &query_aln,
                    &subject_aln,
                    show_subject_header,
                    subject_deflines.get(&hit.subject_id).map(String::as_str),
                )?;
            }
            write_pairwise_db_query_stats(&mut writer, search_space)?;
        }
        write_pairwise_db_database_footer(&mut writer, &db, args)?;
    } else if outfmt_num == 5 {
        let database_label = args
            .db
            .as_ref()
            .map(|db| db.display().to_string())
            .unwrap_or_default();
        write_blastn_db_xml_output(
            &mut writer,
            &hits,
            records,
            args,
            &db,
            &database_label,
            &db_xml_groups,
            &db_xml_metadata,
        )?;
    } else if outfmt_num == 17 {
        write_blastn_sam_output(&mut writer, &hits, records, &db_sam_labels)?;
    } else {
        let database_label = args
            .db
            .as_ref()
            .map(|db| db.display().to_string())
            .unwrap_or_else(|| "N/A".to_string());
        write_blastn_tabular_output(&mut writer, &hits, &args.outfmt, records, &database_label)?;
    }
    writer.flush()?;

    Ok(())
}

fn write_tabular_output<W: Write>(
    writer: &mut W,
    hits: &[TabularHit],
    outfmt: &str,
) -> std::io::Result<()> {
    let outfmt_parts: Vec<&str> = outfmt.split_whitespace().collect();
    let outfmt_num: i32 = outfmt_parts
        .first()
        .and_then(|s| s.parse().ok())
        .unwrap_or(6);
    if outfmt_num == 10 {
        let cols = if outfmt_parts.len() > 1 {
            outfmt_parts[1..].join(" ")
        } else {
            blast_rs::format::DEFAULT_TABULAR_COLUMNS.to_string()
        };
        blast_rs::format::format_tabular_custom_with_delimiter(writer, hits, &cols, ",")
    } else if outfmt_parts.len() > 1 {
        let cols = outfmt_parts[1..].join(" ");
        blast_rs::format::format_tabular_custom(writer, hits, &cols)
    } else {
        format_tabular(writer, hits)
    }
}

fn write_blastn_tabular_output<W: Write>(
    writer: &mut W,
    hits: &[TabularHit],
    outfmt: &str,
    queries: &[blast_rs::input::FastaRecord],
    database_label: &str,
) -> std::io::Result<()> {
    let outfmt_parts: Vec<&str> = outfmt.split_whitespace().collect();
    let outfmt_num: i32 = outfmt_parts
        .first()
        .and_then(|s| s.parse().ok())
        .unwrap_or(6);
    if outfmt_num != 7 {
        return write_tabular_output(writer, hits, outfmt);
    }

    let cols = if outfmt_parts.len() > 1 {
        outfmt_parts[1..].join(" ")
    } else {
        blast_rs::format::DEFAULT_TABULAR_COLUMNS.to_string()
    };
    let cols_for_header = blast_rs::format::expanded_column_tokens(&cols);
    let field_names = cols_for_header
        .iter()
        .map(|c| blast_rs::format::field_display_name(c))
        .collect::<Vec<_>>()
        .join(", ");

    for query in queries {
        let query_hits: Vec<TabularHit> = hits
            .iter()
            .filter(|hit| hit.query_id == query.id)
            .cloned()
            .collect();
        writeln!(writer, "# BLASTN 2.12.0+")?;
        writeln!(writer, "# Query: {}", query.defline)?;
        writeln!(writer, "# Database: {}", database_label)?;
        if !query_hits.is_empty() {
            writeln!(writer, "# Fields: {}", field_names)?;
        }
        writeln!(writer, "# {} hits found", query_hits.len())?;
        if !query_hits.is_empty() {
            blast_rs::format::format_tabular_custom(writer, &query_hits, &cols)?;
        }
    }
    writeln!(writer, "# BLAST processed {} queries", queries.len())?;
    Ok(())
}

fn apply_blastn_dust_mask(seq: &mut [u8]) {
    let window = seq.len().min(64);
    if window >= 3 {
        let mask = blast_rs::filter::dust_filter(seq, window, 2.0);
        mask.apply(seq, 14);
    }
}

fn apply_lowercase_mask(raw_query: &[u8], encoded_query: &mut [u8]) {
    for (raw, encoded) in raw_query.iter().zip(encoded_query.iter_mut()) {
        if raw.is_ascii_lowercase() {
            *encoded = 14;
        }
    }
}

fn query_loc_bounds(
    args: &BlastnArgs,
    query_len: usize,
) -> Result<(usize, usize), Box<dyn std::error::Error>> {
    let Some(loc) = args.query_loc.as_deref() else {
        return Ok((0, query_len));
    };
    loc_bounds(loc, query_len, "query_loc")
}

fn subject_loc_bounds(
    args: &BlastnArgs,
    subject_len: usize,
) -> Result<(usize, usize), Box<dyn std::error::Error>> {
    let Some(loc) = args.subject_loc.as_deref() else {
        return Ok((0, subject_len));
    };
    loc_bounds(loc, subject_len, "subject_loc")
}

fn loc_bounds(
    loc: &str,
    seq_len: usize,
    arg_name: &str,
) -> Result<(usize, usize), Box<dyn std::error::Error>> {
    let Some((start, end)) = loc.split_once('-') else {
        return Err(format!("Invalid {} {:?}: expected start-end", arg_name, loc).into());
    };
    let start: usize = start
        .parse()
        .map_err(|_| format!("Invalid {} {:?}: bad start", arg_name, loc))?;
    let end: usize = end
        .parse()
        .map_err(|_| format!("Invalid {} {:?}: bad end", arg_name, loc))?;
    if start == 0 || end < start || end > seq_len {
        return Err(format!(
            "Invalid {} {:?} for sequence length {}",
            arg_name, loc, seq_len
        )
        .into());
    }
    Ok((start - 1, end))
}

/// Search query against subject FASTA sequences (no database needed).
fn run_blastn_subject(
    args: &BlastnArgs,
    queries: &[blast_rs::input::FastaRecord],
    subjects: &[blast_rs::input::FastaRecord],
) -> Result<(), Box<dyn std::error::Error>> {
    use blast_rs::search::blastn_gapped_search_nomask;

    let total_subj_len: usize = subjects
        .iter()
        .map(|s| subject_loc_bounds(args, s.sequence.len()).map(|(start, end)| end - start))
        .collect::<Result<Vec<_>, _>>()?
        .into_iter()
        .sum();
    let word_size = args.word_size() as usize;
    let search_plus = args.strand != "minus";
    let search_minus = args.strand != "plus";

    let mut all_hits = Vec::new();

    for query_rec in queries {
        let (loc_start, loc_end) = query_loc_bounds(args, query_rec.sequence.len())?;
        let raw_query = &query_rec.sequence[loc_start..loc_end];
        let query_plus_nomask: Vec<u8> = raw_query.iter().map(|&b| iupacna_to_blastna(b)).collect();
        let query_minus_nomask: Vec<u8> = query_plus_nomask
            .iter()
            .rev()
            .map(|&b| complement_blastna(b))
            .collect();
        let mut query_plus = query_plus_nomask.clone();
        if args.dust != "no" {
            apply_blastn_dust_mask(&mut query_plus);
        }
        if args.lcase_masking {
            apply_lowercase_mask(raw_query, &mut query_plus);
        }
        let query_minus: Vec<u8> = query_plus
            .iter()
            .rev()
            .map(|&b| complement_blastna(b))
            .collect();
        let (kbp, search_space, _) = blastn_subject_stats(
            args,
            &query_plus_nomask,
            total_subj_len as i64,
            subjects.len() as i32,
        );

        for subj_rec in subjects {
            let (subj_loc_start, subj_loc_end) = subject_loc_bounds(args, subj_rec.sequence.len())?;
            let raw_subject = &subj_rec.sequence[subj_loc_start..subj_loc_end];
            let subject_offset = subj_loc_start as i32;
            let subject: Vec<u8> = raw_subject.iter().map(|&b| iupacna_to_blastna(b)).collect();

            let hsps = blastn_gapped_search_nomask(
                &query_plus,
                &query_minus,
                &query_plus_nomask,
                &query_minus_nomask,
                &subject,
                word_size,
                args.reward(),
                args.penalty(),
                args.gapopen(),
                args.gapextend(),
                20,
                &kbp,
                search_space,
                args.evalue,
            );

            for hsp in hsps {
                if (hsp.context == 0 && !search_plus) || (hsp.context == 1 && !search_minus) {
                    continue;
                }
                let query_len = query_plus_nomask.len() as i32;
                let query_offset = loc_start as i32;
                let (q_start, q_end) = if hsp.context == 1 {
                    (
                        query_offset + query_len - hsp.query_end + 1,
                        query_offset + query_len - hsp.query_start,
                    )
                } else {
                    (
                        query_offset + hsp.query_start + 1,
                        query_offset + hsp.query_end,
                    )
                };
                let (s_start, s_end) = if hsp.context == 1 {
                    (
                        subject_offset + hsp.subject_end,
                        subject_offset + hsp.subject_start + 1,
                    )
                } else {
                    (
                        subject_offset + hsp.subject_start + 1,
                        subject_offset + hsp.subject_end,
                    )
                };

                all_hits.push(TabularHit {
                    query_id: query_rec.id.clone(),
                    subject_id: subj_rec.id.clone(),
                    subject_title: String::new(),
                    pct_identity: if hsp.align_length > 0 {
                        100.0 * hsp.num_ident as f64 / hsp.align_length as f64
                    } else {
                        0.0
                    },
                    align_len: hsp.align_length,
                    mismatches: hsp.mismatches,
                    gap_opens: hsp.gap_opens,
                    query_start: q_start,
                    query_end: q_end,
                    subject_start: s_start,
                    subject_end: s_end,
                    evalue: hsp.evalue,
                    bit_score: hsp.bit_score,
                    query_len,
                    subject_len: subj_rec.sequence.len() as i32,
                    raw_score: hsp.score,
                    qseq: hsp.qseq.clone(),
                    sseq: hsp.sseq.clone(),
                    qframe: 1,
                    sframe: if hsp.context == 1 { -1 } else { 1 },
                    subject_taxids: vec![],
                    subject_sci_name: String::new(),
                    subject_common_name: String::new(),
                    subject_blast_name: String::new(),
                    subject_kingdom: String::new(),
                    num_ident: hsp.num_ident,
                });
            }
        }
    }

    all_hits.sort_by(|a, b| {
        let a_subject_lo = a.subject_start.min(a.subject_end);
        let b_subject_lo = b.subject_start.min(b.subject_end);
        let a_subject_hi = a.subject_start.max(a.subject_end);
        let b_subject_hi = b.subject_start.max(b.subject_end);

        a.evalue
            .partial_cmp(&b.evalue)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| b.raw_score.cmp(&a.raw_score))
            .then_with(|| b.subject_id.cmp(&a.subject_id))
            .then_with(|| a_subject_lo.cmp(&b_subject_lo))
            .then_with(|| a_subject_hi.cmp(&b_subject_hi))
            .then_with(|| hsp_query_order_start(a).cmp(&hsp_query_order_start(b)))
            .then_with(|| b.sframe.cmp(&a.sframe))
    });
    let max_subjects = args.effective_max_target_seqs() as usize;
    let mut seen_per_query: std::collections::HashMap<String, std::collections::HashSet<String>> =
        std::collections::HashMap::new();
    all_hits.retain(|hit| {
        let seen = seen_per_query.entry(hit.query_id.clone()).or_default();
        if seen.len() >= max_subjects && !seen.contains(&hit.subject_id) {
            false
        } else {
            seen.insert(hit.subject_id.clone());
            true
        }
    });
    let query_len = queries
        .first()
        .map(|r| r.sequence.len() as i32)
        .unwrap_or(0);
    apply_filters(&mut all_hits, args, query_len);
    let subject_deflines: std::collections::HashMap<String, String> = subjects
        .iter()
        .map(|rec| (rec.id.clone(), rec.defline.clone()))
        .collect();

    let stdout = io::stdout();
    let mut writer: Box<dyn Write> = if let Some(ref path) = args.out {
        Box::new(BufWriter::new(File::create(path)?))
    } else {
        Box::new(BufWriter::new(stdout.lock()))
    };
    let outfmt_parts: Vec<&str> = args.outfmt.split_whitespace().collect();
    let outfmt_num: i32 = outfmt_parts[0].parse().unwrap_or(6);
    if outfmt_num == 0 {
        write_pairwise_subject_report_preamble(&mut writer, subjects, args, total_subj_len as i64)?;
        for query in queries {
            let query_hits =
                pairwise_query_hits(&all_hits, &query.id, args.sorthits, args.sorthsps);
            write_pairwise_subject_query_header(&mut writer, query, subjects, &query_hits, args)?;
            let mut last_pairwise_subject: Option<&str> = None;
            for hit in query_hits {
                let mut query_aln = hit
                    .qseq
                    .as_deref()
                    .map(alignment_string_to_blastna)
                    .unwrap_or_default();
                let mut subject_aln = hit
                    .sseq
                    .as_deref()
                    .map(alignment_string_to_blastna)
                    .unwrap_or_default();
                if hit.subject_start > hit.subject_end {
                    query_aln = reverse_complement_alignment_blastna(&query_aln);
                    subject_aln = reverse_complement_alignment_blastna(&subject_aln);
                }
                let has_previous_subject = last_pairwise_subject.is_some();
                let show_subject_header = last_pairwise_subject != Some(hit.subject_id.as_str());
                last_pairwise_subject = Some(hit.subject_id.as_str());
                if show_subject_header && has_previous_subject {
                    writeln!(writer)?;
                }
                let subject_display = subject_deflines
                    .get(&hit.subject_id)
                    .map(|defline| format!(" {}", defline));
                write_pairwise_alignment(
                    &mut writer,
                    hit,
                    &query_aln,
                    &subject_aln,
                    show_subject_header,
                    subject_display.as_deref(),
                )?;
            }
            write_pairwise_subject_query_stats(
                &mut writer,
                args,
                query,
                total_subj_len as i64,
                subjects.len() as i32,
            )?;
        }
        write_pairwise_subject_database_footer(&mut writer, subjects, args, total_subj_len as i64)?;
    } else if outfmt_num == 5 {
        write_blastn_subject_xml_output(
            &mut writer,
            &all_hits,
            queries,
            subjects,
            args,
            total_subj_len as i64,
        )?;
    } else if outfmt_num == 17 {
        write_blastn_subject_sam_output(&mut writer, &all_hits, queries, subjects)?;
    } else {
        let database_label = args
            .subject
            .as_ref()
            .map(|subject| format!("User specified sequence set (Input: {})", subject.display()))
            .unwrap_or_else(|| "User specified sequence set".to_string());
        write_blastn_tabular_output(
            &mut writer,
            &all_hits,
            &args.outfmt,
            queries,
            &database_label,
        )?;
    }
    writer.flush()?;
    Ok(())
}

fn format_sam_float(value: f64) -> String {
    if value != 0.0 && value.abs() < 0.0001 {
        format!("{:.5e}", value)
    } else {
        let abs = value.abs();
        let digits_before_decimal = if abs >= 1.0 {
            abs.log10().floor() as i32 + 1
        } else {
            0
        };
        let decimals = (6 - digits_before_decimal).max(0) as usize;
        let mut s = format!("{value:.decimals$}");
        while s.contains('.') && s.ends_with('0') {
            s.pop();
        }
        if s.ends_with('.') {
            s.pop();
        }
        s
    }
}

fn xml_escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for ch in s.chars() {
        match ch {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            '\'' => out.push_str("&apos;"),
            _ => out.push(ch),
        }
    }
    out
}

fn format_xml_evalue(value: f64) -> String {
    if value == 0.0 {
        "0".to_string()
    } else {
        let s = format!("{value:.5e}");
        if let Some((mantissa, exponent)) = s.split_once('e') {
            let exponent = exponent
                .trim_start_matches('+')
                .trim_start_matches('0')
                .replace("-0", "-");
            format!("{mantissa}e{exponent}")
        } else {
            s
        }
    }
}

fn format_xml_stat_float(value: f64) -> String {
    if value.abs() < 1.0 {
        format!("{value:.15}")
    } else {
        format!("{value:.14}")
    }
}

fn blastn_xml_midline(hit: &TabularHit) -> String {
    let (Some(qseq), Some(sseq)) = (hit.qseq.as_deref(), hit.sseq.as_deref()) else {
        return "|".repeat(hit.num_ident.max(0) as usize);
    };
    qseq.bytes()
        .zip(sseq.bytes())
        .map(|(q, s)| if q != b'-' && q == s { '|' } else { ' ' })
        .collect()
}

fn hsp_query_order_start(hit: &TabularHit) -> i32 {
    if hit.sframe < 0 {
        hit.query_len - hit.query_start.max(hit.query_end) + 1
    } else {
        hit.query_start.min(hit.query_end)
    }
}

fn pairwise_query_hits<'a>(
    hits: &'a [TabularHit],
    query_id: &str,
    sorthits: i32,
    sorthsps: i32,
) -> Vec<&'a TabularHit> {
    let query_hits: Vec<&TabularHit> = hits.iter().filter(|hit| hit.query_id == query_id).collect();
    if query_hits.len() <= 1 {
        return query_hits;
    }

    let mut groups: Vec<Vec<&TabularHit>> = Vec::new();
    let mut start = 0;
    while start < query_hits.len() {
        let subject_id = query_hits[start].subject_id.as_str();
        let mut end = start + 1;
        while end < query_hits.len() && query_hits[end].subject_id == subject_id {
            end += 1;
        }
        let mut group = query_hits[start..end].to_vec();
        if sorthsps != 0 {
            group.sort_by(|a, b| compare_pairwise_hsps(a, b, sorthsps));
        }
        groups.push(group);
        start = end;
    }

    if sorthits != 0 {
        groups.sort_by(|a, b| compare_pairwise_hit_groups(a, b, sorthits));
    }

    groups.into_iter().flatten().collect()
}

fn compare_pairwise_hsps(a: &TabularHit, b: &TabularHit, sorthsps: i32) -> std::cmp::Ordering {
    let ord = match sorthsps {
        1 => b
            .raw_score
            .cmp(&a.raw_score)
            .then_with(|| b.bit_score.partial_cmp(&a.bit_score).unwrap_or(std::cmp::Ordering::Equal)),
        2 => a
            .query_start
            .min(a.query_end)
            .cmp(&b.query_start.min(b.query_end)),
        3 => b
            .pct_identity
            .partial_cmp(&a.pct_identity)
            .unwrap_or(std::cmp::Ordering::Equal),
        4 => a
            .subject_start
            .min(a.subject_end)
            .cmp(&b.subject_start.min(b.subject_end)),
        _ => a
            .evalue
            .partial_cmp(&b.evalue)
            .unwrap_or(std::cmp::Ordering::Equal),
    };

    ord.then_with(|| a.evalue.partial_cmp(&b.evalue).unwrap_or(std::cmp::Ordering::Equal))
        .then_with(|| b.raw_score.cmp(&a.raw_score))
        .then_with(|| hsp_query_order_start(a).cmp(&hsp_query_order_start(b)))
        .then_with(|| a.subject_start.min(a.subject_end).cmp(&b.subject_start.min(b.subject_end)))
        .then_with(|| b.sframe.cmp(&a.sframe))
}

fn pairwise_best_hit<'a>(hits: &'a [&'a TabularHit]) -> &'a TabularHit {
    hits.iter()
        .copied()
        .min_by(|a, b| {
            a.evalue
                .partial_cmp(&b.evalue)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| b.raw_score.cmp(&a.raw_score))
        })
        .expect("pairwise hit group should not be empty")
}

fn pairwise_total_bit_score(hits: &[&TabularHit]) -> f64 {
    hits.iter().map(|hit| hit.bit_score).sum()
}

fn format_pairwise_total_bit_score(hits: &[&TabularHit]) -> String {
    blast_rs::format::format_bitscore(pairwise_total_bit_score(hits))
}

fn pairwise_query_coverage(hits: &[&TabularHit]) -> i32 {
    let Some(first) = hits.first() else {
        return 0;
    };
    if first.query_len <= 0 {
        return 0;
    }

    let mut ranges: Vec<(i32, i32)> = hits
        .iter()
        .map(|hit| {
            (
                hit.query_start.min(hit.query_end),
                hit.query_start.max(hit.query_end),
            )
        })
        .collect();
    ranges.sort_unstable_by_key(|&(start, end)| (start, end));

    let mut covered = 0;
    let mut current = ranges[0];
    for &(start, end) in &ranges[1..] {
        if start <= current.1 + 1 {
            current.1 = current.1.max(end);
        } else {
            covered += current.1 - current.0 + 1;
            current = (start, end);
        }
    }
    covered += current.1 - current.0 + 1;

    let percent = (100.0 * covered as f64 / first.query_len as f64) as i32;
    if percent > 100 {
        99
    } else {
        percent
    }
}

fn pairwise_max_identity(hits: &[&TabularHit]) -> i32 {
    hits.iter()
        .map(|hit| hit.pct_identity.round() as i32)
        .max()
        .unwrap_or(0)
}

fn compare_pairwise_hit_groups(
    a: &[&TabularHit],
    b: &[&TabularHit],
    sorthits: i32,
) -> std::cmp::Ordering {
    let a_best = pairwise_best_hit(a);
    let b_best = pairwise_best_hit(b);

    let ord = match sorthits {
        1 => b_best
            .bit_score
            .partial_cmp(&a_best.bit_score)
            .unwrap_or(std::cmp::Ordering::Equal),
        2 => pairwise_total_bit_score(b)
            .partial_cmp(&pairwise_total_bit_score(a))
            .unwrap_or(std::cmp::Ordering::Equal),
        3 => pairwise_max_identity(b).cmp(&pairwise_max_identity(a)),
        4 => pairwise_query_coverage(b).cmp(&pairwise_query_coverage(a)),
        _ => a_best
            .evalue
            .partial_cmp(&b_best.evalue)
            .unwrap_or(std::cmp::Ordering::Equal),
    };

    ord.then_with(|| {
        a_best
            .evalue
            .partial_cmp(&b_best.evalue)
            .unwrap_or(std::cmp::Ordering::Equal)
    })
    .then_with(|| b_best.raw_score.cmp(&a_best.raw_score))
    .then_with(|| a_best.subject_id.cmp(&b_best.subject_id))
}

fn write_pairwise_hit_summary_header<W: Write>(writer: &mut W, sorthits: i32) -> io::Result<()> {
    if sorthits == 0 {
        writeln!(
            writer,
            "                                                                      Score     E"
        )?;
        writeln!(
            writer,
            "Sequences producing significant alignments:                          (Bits)  Value"
        )?;
    } else {
        writeln!(
            writer,
            "                                                                      Score   Total  Query    E     Max"
        )?;
        writeln!(
            writer,
            "Sequences producing significant alignments:                          (Bits)   Score  cover   Value  Ident"
        )?;
    }
    Ok(())
}

fn write_pairwise_hit_summary_row<W: Write>(
    writer: &mut W,
    desc: &str,
    hits: &[&TabularHit],
    sorthits: i32,
) -> io::Result<()> {
    let best = pairwise_best_hit(hits);
    if sorthits == 0 {
        writeln!(
            writer,
            "{:<68}  {:>4}    {}",
            desc,
            blast_rs::format::format_bitscore(best.bit_score),
            blast_rs::format::format_pairwise_evalue(best.evalue),
        )
    } else {
        fn pad_to_column(line: &mut String, column: usize) {
            if line.len() < column {
                line.push_str(&" ".repeat(column - line.len()));
            }
        }

        let max_identity = format!("{}%", pairwise_max_identity(hits));
        let query_coverage = format!("{}%", pairwise_query_coverage(hits));
        let mut line = format!(
            "{:<68}  {}",
            desc,
            blast_rs::format::format_bitscore(best.bit_score)
        );
        pad_to_column(&mut line, 79);
        line.push_str(&format_pairwise_total_bit_score(hits));
        pad_to_column(&mut line, 86);
        line.push_str(&query_coverage);
        pad_to_column(&mut line, 93);
        line.push_str(&blast_rs::format::format_pairwise_evalue(best.evalue));
        pad_to_column(&mut line, 100);
        line.push_str(&max_identity);
        pad_to_column(&mut line, 110);
        writeln!(writer, "{line}")
    }
}

fn write_blastn_subject_xml_output<W: Write>(
    writer: &mut W,
    hits: &[TabularHit],
    queries: &[blast_rs::input::FastaRecord],
    subjects: &[blast_rs::input::FastaRecord],
    args: &BlastnArgs,
    total_subject_len: i64,
) -> std::io::Result<()> {
    let subject_deflines: std::collections::HashMap<&str, &str> = subjects
        .iter()
        .map(|rec| (rec.id.as_str(), rec.defline.as_str()))
        .collect();
    let subject_accessions: std::collections::HashMap<&str, String> = subjects
        .iter()
        .enumerate()
        .map(|(i, rec)| (rec.id.as_str(), format!("Subject_{}", i + 1)))
        .collect();

    writeln!(writer, "<?xml version=\"1.0\"?>")?;
    writeln!(writer, "<!DOCTYPE BlastOutput PUBLIC \"-//NCBI//NCBI BlastOutput/EN\" \"http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd\">")?;
    writeln!(writer, "<BlastOutput>")?;
    writeln!(
        writer,
        "  <BlastOutput_program>blastn</BlastOutput_program>"
    )?;
    writeln!(
        writer,
        "  <BlastOutput_version>BLASTN 2.12.0+</BlastOutput_version>"
    )?;
    writeln!(writer, "  <BlastOutput_reference>Stephen F. Altschul, Thomas L. Madden, Alejandro A. Sch&amp;auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), &quot;Gapped BLAST and PSI-BLAST: a new generation of protein database search programs&quot;, Nucleic Acids Res. 25:3389-3402.</BlastOutput_reference>")?;
    writeln!(writer, "  <BlastOutput_db></BlastOutput_db>")?;
    if let Some(query) = queries.first() {
        writeln!(
            writer,
            "  <BlastOutput_query-ID>Query_1</BlastOutput_query-ID>"
        )?;
        writeln!(
            writer,
            "  <BlastOutput_query-def>{}</BlastOutput_query-def>",
            xml_escape(&query.defline)
        )?;
        writeln!(
            writer,
            "  <BlastOutput_query-len>{}</BlastOutput_query-len>",
            query.sequence.len()
        )?;
    }
    writeln!(writer, "  <BlastOutput_param>")?;
    writeln!(writer, "    <Parameters>")?;
    writeln!(
        writer,
        "      <Parameters_expect>{}</Parameters_expect>",
        format_sam_float(args.evalue)
    )?;
    writeln!(
        writer,
        "      <Parameters_sc-match>{}</Parameters_sc-match>",
        args.reward()
    )?;
    writeln!(
        writer,
        "      <Parameters_sc-mismatch>{}</Parameters_sc-mismatch>",
        args.penalty()
    )?;
    writeln!(
        writer,
        "      <Parameters_gap-open>{}</Parameters_gap-open>",
        args.gapopen()
    )?;
    writeln!(
        writer,
        "      <Parameters_gap-extend>{}</Parameters_gap-extend>",
        args.gapextend()
    )?;
    writeln!(writer, "      <Parameters_filter>m;</Parameters_filter>")?;
    writeln!(writer, "    </Parameters>")?;
    writeln!(writer, "  </BlastOutput_param>")?;
    writeln!(writer, "<BlastOutput_iterations>")?;

    for (query_index, query) in queries.iter().enumerate() {
        let (loc_start, loc_end) = query_loc_bounds(args, query.sequence.len())
            .map_err(|err| io::Error::new(io::ErrorKind::InvalidInput, err.to_string()))?;
        let query_plus_nomask: Vec<u8> = query.sequence[loc_start..loc_end]
            .iter()
            .map(|&b| iupacna_to_blastna(b))
            .collect();
        let (kbp, search_space, len_adj) = blastn_subject_stats(
            args,
            &query_plus_nomask,
            total_subject_len,
            subjects.len() as i32,
        );

        writeln!(writer, "<Iteration>")?;
        writeln!(
            writer,
            "  <Iteration_iter-num>{}</Iteration_iter-num>",
            query_index + 1
        )?;
        writeln!(
            writer,
            "  <Iteration_query-ID>Query_{}</Iteration_query-ID>",
            query_index + 1
        )?;
        writeln!(
            writer,
            "  <Iteration_query-def>{}</Iteration_query-def>",
            xml_escape(&query.defline)
        )?;
        writeln!(
            writer,
            "  <Iteration_query-len>{}</Iteration_query-len>",
            query.sequence.len()
        )?;
        writeln!(writer, "<Iteration_hits>")?;

        let mut subject_ord = std::collections::BTreeMap::<&str, Vec<&TabularHit>>::new();
        for hit in hits.iter().filter(|hit| hit.query_id == query.id) {
            subject_ord
                .entry(hit.subject_id.as_str())
                .or_default()
                .push(hit);
        }
        for (hit_num, (subject_id, hsps)) in subject_ord.iter().enumerate() {
            let Some(first) = hsps.first() else {
                continue;
            };
            writeln!(writer, "<Hit>")?;
            writeln!(writer, "  <Hit_num>{}</Hit_num>", hit_num + 1)?;
            writeln!(writer, "  <Hit_id>{}</Hit_id>", xml_escape(subject_id))?;
            writeln!(
                writer,
                "  <Hit_def>{}</Hit_def>",
                xml_escape(
                    subject_deflines
                        .get(subject_id)
                        .copied()
                        .unwrap_or(subject_id)
                )
            )?;
            writeln!(
                writer,
                "  <Hit_accession>{}</Hit_accession>",
                subject_accessions
                    .get(subject_id)
                    .map(String::as_str)
                    .unwrap_or(subject_id)
            )?;
            writeln!(writer, "  <Hit_len>{}</Hit_len>", first.subject_len)?;
            writeln!(writer, "  <Hit_hsps>")?;
            for (hsp_num, hit) in hsps.iter().enumerate() {
                writeln!(writer, "    <Hsp>")?;
                writeln!(writer, "      <Hsp_num>{}</Hsp_num>", hsp_num + 1)?;
                writeln!(
                    writer,
                    "      <Hsp_bit-score>{:.2}</Hsp_bit-score>",
                    hit.bit_score
                )?;
                writeln!(writer, "      <Hsp_score>{}</Hsp_score>", hit.raw_score)?;
                writeln!(
                    writer,
                    "      <Hsp_evalue>{}</Hsp_evalue>",
                    format_xml_evalue(hit.evalue)
                )?;
                writeln!(
                    writer,
                    "      <Hsp_query-from>{}</Hsp_query-from>",
                    hit.query_start
                )?;
                writeln!(
                    writer,
                    "      <Hsp_query-to>{}</Hsp_query-to>",
                    hit.query_end
                )?;
                writeln!(
                    writer,
                    "      <Hsp_hit-from>{}</Hsp_hit-from>",
                    hit.subject_start
                )?;
                writeln!(writer, "      <Hsp_hit-to>{}</Hsp_hit-to>", hit.subject_end)?;
                writeln!(
                    writer,
                    "      <Hsp_query-frame>{}</Hsp_query-frame>",
                    hit.qframe
                )?;
                writeln!(
                    writer,
                    "      <Hsp_hit-frame>{}</Hsp_hit-frame>",
                    hit.sframe
                )?;
                writeln!(
                    writer,
                    "      <Hsp_identity>{}</Hsp_identity>",
                    hit.num_ident
                )?;
                writeln!(
                    writer,
                    "      <Hsp_positive>{}</Hsp_positive>",
                    hit.num_ident
                )?;
                writeln!(writer, "      <Hsp_gaps>{}</Hsp_gaps>", sam_gap_count(hit))?;
                writeln!(
                    writer,
                    "      <Hsp_align-len>{}</Hsp_align-len>",
                    hit.align_len
                )?;
                writeln!(
                    writer,
                    "      <Hsp_qseq>{}</Hsp_qseq>",
                    xml_escape(hit.qseq.as_deref().unwrap_or(""))
                )?;
                writeln!(
                    writer,
                    "      <Hsp_hseq>{}</Hsp_hseq>",
                    xml_escape(hit.sseq.as_deref().unwrap_or(""))
                )?;
                writeln!(
                    writer,
                    "      <Hsp_midline>{}</Hsp_midline>",
                    xml_escape(&blastn_xml_midline(hit))
                )?;
                writeln!(writer, "    </Hsp>")?;
            }
            writeln!(writer, "  </Hit_hsps>")?;
            writeln!(writer, "</Hit>")?;
        }

        writeln!(writer, "</Iteration_hits>")?;
        writeln!(writer, "  <Iteration_stat>")?;
        writeln!(writer, "    <Statistics>")?;
        writeln!(writer, "      <Statistics_db-num>0</Statistics_db-num>")?;
        writeln!(writer, "      <Statistics_db-len>0</Statistics_db-len>")?;
        writeln!(
            writer,
            "      <Statistics_hsp-len>{}</Statistics_hsp-len>",
            len_adj
        )?;
        writeln!(
            writer,
            "      <Statistics_eff-space>{}</Statistics_eff-space>",
            search_space.round() as u64
        )?;
        writeln!(
            writer,
            "      <Statistics_kappa>{}</Statistics_kappa>",
            format_xml_stat_float(kbp.k)
        )?;
        writeln!(
            writer,
            "      <Statistics_lambda>{}</Statistics_lambda>",
            format_xml_stat_float(kbp.lambda)
        )?;
        writeln!(
            writer,
            "      <Statistics_entropy>{}</Statistics_entropy>",
            format_xml_stat_float(kbp.h)
        )?;
        writeln!(writer, "    </Statistics>")?;
        writeln!(writer, "  </Iteration_stat>")?;
        writeln!(writer, "</Iteration>")?;
    }
    writeln!(writer, "</BlastOutput_iterations>")?;
    writeln!(writer, "</BlastOutput>")?;
    writeln!(writer)?;
    Ok(())
}

fn write_blastn_db_xml_output<W: Write>(
    writer: &mut W,
    hits: &[TabularHit],
    queries: &[blast_rs::input::FastaRecord],
    args: &BlastnArgs,
    db: &BlastDb,
    database_label: &str,
    hit_groups: &std::collections::HashMap<String, String>,
    hit_metadata: &std::collections::HashMap<String, (String, String, String, i32)>,
) -> std::io::Result<()> {
    writeln!(writer, "<?xml version=\"1.0\"?>")?;
    writeln!(writer, "<!DOCTYPE BlastOutput PUBLIC \"-//NCBI//NCBI BlastOutput/EN\" \"http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd\">")?;
    writeln!(writer, "<BlastOutput>")?;
    writeln!(
        writer,
        "  <BlastOutput_program>blastn</BlastOutput_program>"
    )?;
    writeln!(
        writer,
        "  <BlastOutput_version>BLASTN 2.12.0+</BlastOutput_version>"
    )?;
    writeln!(writer, "  <BlastOutput_reference>Stephen F. Altschul, Thomas L. Madden, Alejandro A. Sch&amp;auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), &quot;Gapped BLAST and PSI-BLAST: a new generation of protein database search programs&quot;, Nucleic Acids Res. 25:3389-3402.</BlastOutput_reference>")?;
    writeln!(
        writer,
        "  <BlastOutput_db>{}</BlastOutput_db>",
        xml_escape(database_label)
    )?;
    if let Some(query) = queries.first() {
        writeln!(
            writer,
            "  <BlastOutput_query-ID>Query_1</BlastOutput_query-ID>"
        )?;
        writeln!(
            writer,
            "  <BlastOutput_query-def>{}</BlastOutput_query-def>",
            xml_escape(&query.defline)
        )?;
        writeln!(
            writer,
            "  <BlastOutput_query-len>{}</BlastOutput_query-len>",
            query.sequence.len()
        )?;
    }
    writeln!(writer, "  <BlastOutput_param>")?;
    writeln!(writer, "    <Parameters>")?;
    writeln!(
        writer,
        "      <Parameters_expect>{}</Parameters_expect>",
        format_sam_float(args.evalue)
    )?;
    writeln!(
        writer,
        "      <Parameters_sc-match>{}</Parameters_sc-match>",
        args.reward()
    )?;
    writeln!(
        writer,
        "      <Parameters_sc-mismatch>{}</Parameters_sc-mismatch>",
        args.penalty()
    )?;
    writeln!(
        writer,
        "      <Parameters_gap-open>{}</Parameters_gap-open>",
        args.gapopen()
    )?;
    writeln!(
        writer,
        "      <Parameters_gap-extend>{}</Parameters_gap-extend>",
        args.gapextend()
    )?;
    writeln!(writer, "      <Parameters_filter>m;</Parameters_filter>")?;
    writeln!(writer, "    </Parameters>")?;
    writeln!(writer, "  </BlastOutput_param>")?;
    writeln!(writer, "<BlastOutput_iterations>")?;

    for (query_index, query) in queries.iter().enumerate() {
        let (loc_start, loc_end) = query_loc_bounds(args, query.sequence.len())
            .map_err(|err| io::Error::new(io::ErrorKind::InvalidInput, err.to_string()))?;
        let query_plus_nomask: Vec<u8> = query.sequence[loc_start..loc_end]
            .iter()
            .map(|&b| iupacna_to_blastna(b))
            .collect();
        let (kbp, search_space, len_adj) = blastn_subject_stats(
            args,
            &query_plus_nomask,
            db.total_length as i64,
            db.num_oids as i32,
        );

        writeln!(writer, "<Iteration>")?;
        writeln!(
            writer,
            "  <Iteration_iter-num>{}</Iteration_iter-num>",
            query_index + 1
        )?;
        writeln!(
            writer,
            "  <Iteration_query-ID>Query_{}</Iteration_query-ID>",
            query_index + 1
        )?;
        writeln!(
            writer,
            "  <Iteration_query-def>{}</Iteration_query-def>",
            xml_escape(&query.defline)
        )?;
        writeln!(
            writer,
            "  <Iteration_query-len>{}</Iteration_query-len>",
            query.sequence.len()
        )?;
        writeln!(writer, "<Iteration_hits>")?;

        let mut group_order = Vec::<String>::new();
        let mut grouped: std::collections::HashMap<String, Vec<&TabularHit>> =
            std::collections::HashMap::new();
        for hit in hits.iter().filter(|hit| hit.query_id == query.id) {
            let group = hit_groups
                .get(&sam_hit_label_key(hit))
                .cloned()
                .unwrap_or_else(|| hit.subject_id.clone());
            if !grouped.contains_key(&group) {
                group_order.push(group.clone());
            }
            grouped.entry(group).or_default().push(hit);
        }

        for (hit_num, group) in group_order.iter().enumerate() {
            let Some(hsps) = grouped.get(group) else {
                continue;
            };
            let Some(first) = hsps.first() else {
                continue;
            };
            let fallback = (
                first.subject_id.clone(),
                first.subject_id.clone(),
                first.subject_id.clone(),
                first.subject_len,
            );
            let (hit_id, hit_def, hit_accession, hit_len) =
                hit_metadata.get(group).unwrap_or(&fallback);
            writeln!(writer, "<Hit>")?;
            writeln!(writer, "  <Hit_num>{}</Hit_num>", hit_num + 1)?;
            writeln!(writer, "  <Hit_id>{}</Hit_id>", xml_escape(hit_id))?;
            writeln!(writer, "  <Hit_def>{}</Hit_def>", xml_escape(hit_def))?;
            writeln!(
                writer,
                "  <Hit_accession>{}</Hit_accession>",
                xml_escape(hit_accession)
            )?;
            writeln!(writer, "  <Hit_len>{}</Hit_len>", hit_len)?;
            writeln!(writer, "  <Hit_hsps>")?;
            for (hsp_num, hit) in hsps.iter().enumerate() {
                writeln!(writer, "    <Hsp>")?;
                writeln!(writer, "      <Hsp_num>{}</Hsp_num>", hsp_num + 1)?;
                writeln!(
                    writer,
                    "      <Hsp_bit-score>{:.2}</Hsp_bit-score>",
                    hit.bit_score
                )?;
                writeln!(writer, "      <Hsp_score>{}</Hsp_score>", hit.raw_score)?;
                writeln!(
                    writer,
                    "      <Hsp_evalue>{}</Hsp_evalue>",
                    format_xml_evalue(hit.evalue)
                )?;
                writeln!(
                    writer,
                    "      <Hsp_query-from>{}</Hsp_query-from>",
                    hit.query_start
                )?;
                writeln!(
                    writer,
                    "      <Hsp_query-to>{}</Hsp_query-to>",
                    hit.query_end
                )?;
                writeln!(
                    writer,
                    "      <Hsp_hit-from>{}</Hsp_hit-from>",
                    hit.subject_start
                )?;
                writeln!(writer, "      <Hsp_hit-to>{}</Hsp_hit-to>", hit.subject_end)?;
                writeln!(
                    writer,
                    "      <Hsp_query-frame>{}</Hsp_query-frame>",
                    hit.qframe
                )?;
                writeln!(
                    writer,
                    "      <Hsp_hit-frame>{}</Hsp_hit-frame>",
                    hit.sframe
                )?;
                writeln!(
                    writer,
                    "      <Hsp_identity>{}</Hsp_identity>",
                    hit.num_ident
                )?;
                writeln!(
                    writer,
                    "      <Hsp_positive>{}</Hsp_positive>",
                    hit.num_ident
                )?;
                writeln!(writer, "      <Hsp_gaps>{}</Hsp_gaps>", sam_gap_count(hit))?;
                writeln!(
                    writer,
                    "      <Hsp_align-len>{}</Hsp_align-len>",
                    hit.align_len
                )?;
                writeln!(
                    writer,
                    "      <Hsp_qseq>{}</Hsp_qseq>",
                    xml_escape(hit.qseq.as_deref().unwrap_or(""))
                )?;
                writeln!(
                    writer,
                    "      <Hsp_hseq>{}</Hsp_hseq>",
                    xml_escape(hit.sseq.as_deref().unwrap_or(""))
                )?;
                writeln!(
                    writer,
                    "      <Hsp_midline>{}</Hsp_midline>",
                    xml_escape(&blastn_xml_midline(hit))
                )?;
                writeln!(writer, "    </Hsp>")?;
            }
            writeln!(writer, "  </Hit_hsps>")?;
            writeln!(writer, "</Hit>")?;
        }

        writeln!(writer, "</Iteration_hits>")?;
        writeln!(writer, "  <Iteration_stat>")?;
        writeln!(writer, "    <Statistics>")?;
        writeln!(
            writer,
            "      <Statistics_db-num>{}</Statistics_db-num>",
            db.num_oids
        )?;
        writeln!(
            writer,
            "      <Statistics_db-len>{}</Statistics_db-len>",
            db.total_length
        )?;
        writeln!(
            writer,
            "      <Statistics_hsp-len>{}</Statistics_hsp-len>",
            len_adj
        )?;
        writeln!(
            writer,
            "      <Statistics_eff-space>{}</Statistics_eff-space>",
            search_space.round() as u64
        )?;
        writeln!(
            writer,
            "      <Statistics_kappa>{}</Statistics_kappa>",
            format_xml_stat_float(kbp.k)
        )?;
        writeln!(
            writer,
            "      <Statistics_lambda>{}</Statistics_lambda>",
            format_xml_stat_float(kbp.lambda)
        )?;
        writeln!(
            writer,
            "      <Statistics_entropy>{}</Statistics_entropy>",
            format_xml_stat_float(kbp.h)
        )?;
        writeln!(writer, "    </Statistics>")?;
        writeln!(writer, "  </Iteration_stat>")?;
        writeln!(writer, "</Iteration>")?;
    }
    writeln!(writer, "</BlastOutput_iterations>")?;
    writeln!(writer, "</BlastOutput>")?;
    writeln!(writer)?;
    Ok(())
}

fn sam_cigar_subject_as_read(hit: &TabularHit) -> String {
    let (Some(qseq), Some(sseq)) = (hit.qseq.as_deref(), hit.sseq.as_deref()) else {
        return format!("{}M", hit.align_len);
    };

    let mut cigar = String::new();
    let mut current_op = '\0';
    let mut current_len = 0usize;
    for (q, s) in qseq.bytes().zip(sseq.bytes()) {
        let op = if q == b'-' {
            'I'
        } else if s == b'-' {
            'D'
        } else {
            'M'
        };
        if op == current_op {
            current_len += 1;
        } else {
            if current_len > 0 {
                cigar.push_str(&current_len.to_string());
                cigar.push(current_op);
            }
            current_op = op;
            current_len = 1;
        }
    }
    if current_len > 0 {
        cigar.push_str(&current_len.to_string());
        cigar.push(current_op);
    }
    if cigar.is_empty() {
        format!("{}M", hit.align_len)
    } else {
        cigar
    }
}

fn sam_gap_count(hit: &TabularHit) -> i32 {
    let (Some(qseq), Some(sseq)) = (hit.qseq.as_deref(), hit.sseq.as_deref()) else {
        return (hit.align_len - hit.num_ident - hit.mismatches).max(0);
    };
    qseq.bytes()
        .chain(sseq.bytes())
        .filter(|&base| base == b'-')
        .count() as i32
}

fn sam_mismatch_count(hit: &TabularHit) -> i32 {
    let (Some(qseq), Some(sseq)) = (hit.qseq.as_deref(), hit.sseq.as_deref()) else {
        return hit.mismatches;
    };
    qseq.bytes()
        .zip(sseq.bytes())
        .filter(|&(q, s)| q != b'-' && s != b'-' && q != s)
        .count() as i32
}

fn sam_pairwise_identity(hit: &TabularHit) -> f64 {
    let (Some(qseq), Some(sseq)) = (hit.qseq.as_deref(), hit.sseq.as_deref()) else {
        return hit.pct_identity;
    };
    let mut compared = 0usize;
    let mut identical = 0usize;
    for (q, s) in qseq.bytes().zip(sseq.bytes()) {
        if q != b'-' && s != b'-' {
            compared += 1;
            if q == s {
                identical += 1;
            }
        }
    }
    if compared == 0 {
        0.0
    } else {
        100.0 * identical as f64 / compared as f64
    }
}

fn sam_hit_key(hit: &TabularHit) -> (String, String, i32, i32, i32, i32, i32) {
    (
        hit.query_id.clone(),
        hit.subject_id.clone(),
        hit.query_start,
        hit.query_end,
        hit.subject_start,
        hit.subject_end,
        hit.raw_score,
    )
}

fn sam_hit_label_key(hit: &TabularHit) -> String {
    format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}",
        hit.query_id,
        hit.subject_id,
        hit.query_start,
        hit.query_end,
        hit.subject_start,
        hit.subject_end,
        hit.raw_score
    )
}

fn write_blastn_subject_sam_output<W: Write>(
    writer: &mut W,
    hits: &[TabularHit],
    queries: &[blast_rs::input::FastaRecord],
    subjects: &[blast_rs::input::FastaRecord],
) -> std::io::Result<()> {
    let subject_labels: std::collections::HashMap<&str, String> = subjects
        .iter()
        .enumerate()
        .map(|(i, rec)| (rec.id.as_str(), format!("Subject_{}", i + 1)))
        .collect();
    let hit_labels: std::collections::HashMap<String, String> = hits
        .iter()
        .filter_map(|hit| {
            subject_labels
                .get(hit.subject_id.as_str())
                .map(|label| (sam_hit_label_key(hit), label.clone()))
        })
        .collect();

    write_blastn_sam_output(writer, hits, queries, &hit_labels)
}

fn write_blastn_sam_output<W: Write>(
    writer: &mut W,
    hits: &[TabularHit],
    queries: &[blast_rs::input::FastaRecord],
    hit_labels: &std::collections::HashMap<String, String>,
) -> std::io::Result<()> {
    if hits.is_empty() {
        return Ok(());
    }

    let query_labels: std::collections::HashMap<&str, String> = queries
        .iter()
        .enumerate()
        .map(|(i, rec)| (rec.id.as_str(), format!("Query_{}", i + 1)))
        .collect();

    writeln!(writer, "@HD\tVN:1.2\tSO:coordinate\tGO:reference")?;
    for (i, query) in queries.iter().enumerate() {
        writeln!(
            writer,
            "@SQ\tSN:Query_{}\tLN:{}",
            i + 1,
            query.sequence.len()
        )?;
    }
    writeln!(
        writer,
        "@PG\tID:0\tVN:2.12.0+\tCL:blast-cli blastn -outfmt 17 \tPN:blastn"
    )?;

    for hit in hits {
        let subject_label = hit_labels
            .get(&sam_hit_label_key(hit))
            .map(String::as_str)
            .unwrap_or(hit.subject_id.as_str());
        let query_label = query_labels
            .get(hit.query_id.as_str())
            .map(String::as_str)
            .unwrap_or(hit.query_id.as_str());
        let flag = if hit.subject_start > hit.subject_end {
            16
        } else {
            0
        };
        let pos = hit.query_start.min(hit.query_end);
        let cigar = sam_cigar_subject_as_read(hit);
        let edit_distance = sam_gap_count(hit) + sam_mismatch_count(hit);
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t255\t{}\t*\t0\t0\t*\t*\tAS:i:{}\tEV:f:{}\tNM:i:{}\tPI:f:{:.2}\tBS:f:{}",
            subject_label,
            flag,
            query_label,
            pos,
            cigar,
            hit.raw_score,
            format_sam_float(hit.evalue),
            edit_distance,
            sam_pairwise_identity(hit),
            format_sam_float(hit.bit_score),
        )?;
    }
    Ok(())
}

fn compare_oid_hsps_for_hitlist(
    a_oid: u32,
    a_hsps: &[blast_rs::search::SearchHsp],
    b_oid: u32,
    b_hsps: &[blast_rs::search::SearchHsp],
) -> std::cmp::Ordering {
    let a_best = a_hsps.iter().map(|h| h.evalue).fold(f64::MAX, f64::min);
    let b_best = b_hsps.iter().map(|h| h.evalue).fold(f64::MAX, f64::min);
    let a_score = a_hsps.iter().map(|h| h.score).max().unwrap_or(i32::MIN);
    let b_score = b_hsps.iter().map(|h| h.score).max().unwrap_or(i32::MIN);

    a_best
        .partial_cmp(&b_best)
        .unwrap_or(std::cmp::Ordering::Equal)
        .then_with(|| b_score.cmp(&a_score))
        .then_with(|| b_oid.cmp(&a_oid))
}

fn blastn_subject_stats(
    args: &BlastnArgs,
    query_plus: &[u8],
    total_subject_len: i64,
    num_subjects: i32,
) -> (blast_rs::stat::KarlinBlk, f64, i32) {
    const BLASTNA_SIZE: usize = 16;
    const BLASTNA_TO_NCBI4NA: [u8; BLASTNA_SIZE] =
        [1, 2, 4, 8, 5, 10, 3, 12, 9, 6, 14, 13, 11, 7, 15, 0];

    let reward = args.reward();
    let penalty = args.penalty();
    let matrix_fn = |i: usize, j: usize| -> i32 {
        if i >= BLASTNA_SIZE || j >= BLASTNA_SIZE {
            return penalty;
        }
        if i == BLASTNA_SIZE - 1 || j == BLASTNA_SIZE - 1 {
            return i32::MIN / 2;
        }
        if BLASTNA_TO_NCBI4NA[i] & BLASTNA_TO_NCBI4NA[j] != 0 {
            let d = (0..4)
                .filter(|&k| BLASTNA_TO_NCBI4NA[j] & BLASTNA_TO_NCBI4NA[k] != 0)
                .count() as i32;
            (((d - 1) as f64 * penalty as f64 + reward as f64) / d as f64).round() as i32
        } else {
            penalty
        }
    };

    let mut lo = i32::MAX;
    let mut hi = i32::MIN;
    for i in 0..BLASTNA_SIZE {
        for j in 0..BLASTNA_SIZE {
            let s = matrix_fn(i, j);
            if s <= -100000000 || s >= 100000000 {
                continue;
            }
            lo = lo.min(s);
            hi = hi.max(s);
        }
    }

    let qlen = query_plus.len() as i32;
    let ctxs = [blast_rs::stat::UngappedKbpContext {
        query_offset: 0,
        query_length: qlen,
        is_valid: true,
    }];
    let ambig: &[u8] = &[14, 15];
    let default_kbp = blast_rs::stat::KarlinBlk {
        lambda: 1.374,
        k: 0.621,
        log_k: 0.621_f64.ln(),
        h: 1.286,
    };
    let ungapped = blast_rs::stat::ungapped_kbp_calc(
        query_plus,
        &ctxs,
        lo,
        hi,
        BLASTNA_SIZE,
        ambig,
        &matrix_fn,
    )[0]
    .clone()
    .unwrap_or(default_kbp);
    let (kbp, _) = blast_rs::stat::nucl_gapped_kbp_lookup(
        args.gapopen(),
        args.gapextend(),
        reward,
        penalty,
        &ungapped,
    )
    .unwrap_or((ungapped.clone(), false));

    let (alpha, beta) = blast_rs::stat::nucl_alpha_beta(
        reward,
        penalty,
        args.gapopen(),
        args.gapextend(),
        ungapped.lambda,
        ungapped.h,
        true,
    );
    let (len_adj, _) = blast_rs::stat::compute_length_adjustment_exact(
        kbp.k,
        kbp.log_k,
        alpha / kbp.lambda,
        beta,
        qlen,
        total_subject_len,
        num_subjects,
    );
    let eff_db = std::cmp::max(total_subject_len - num_subjects as i64 * len_adj as i64, 1);
    let search_space = eff_db as f64 * (qlen - len_adj).max(1) as f64;
    (kbp, search_space, len_adj)
}

fn db_subject_defline(db: &BlastDb, oid: u32, subject_id: &str) -> Option<String> {
    let title = extract_header_title(db.get_header(oid))?;
    if title == subject_id || title.starts_with(&format!("{} ", subject_id)) {
        Some(title)
    } else {
        Some(format!("{} {}", subject_id, title))
    }
}

fn db_subject_title(db: &BlastDb, oid: u32, subject_id: &str) -> Option<String> {
    let title = extract_header_title(db.get_header(oid))?;
    if title == subject_id {
        return Some(String::new());
    }
    if let Some(rest) = title.strip_prefix(subject_id) {
        if rest.starts_with(char::is_whitespace) {
            return Some(rest.trim_start().to_string());
        }
    }
    Some(title)
}

fn db_pairwise_subject_defline(db: &BlastDb, oid: u32, subject_id: &str) -> Option<String> {
    if subject_id == "BL_ORD_ID" || subject_id.starts_with("gnl|BL_ORD_ID|") {
        db_subject_title(db, oid, subject_id)
    } else {
        db_subject_defline(db, oid, subject_id)
    }
}

fn extract_header_title(hdr: &[u8]) -> Option<String> {
    let mut i = 0;
    while i + 1 < hdr.len() {
        if matches!(hdr[i], 0x1a | 0x0c) {
            if let Some((len, len_len)) = read_ber_len(hdr, i + 1) {
                let start = i + 1 + len_len;
                let end = start.saturating_add(len);
                if len > 0 && end <= hdr.len() {
                    let bytes = &hdr[start..end];
                    if bytes
                        .iter()
                        .all(|&b| b == b'\t' || b == b' ' || b.is_ascii_graphic())
                    {
                        let s = String::from_utf8_lossy(bytes).trim().to_string();
                        if looks_like_header_title(&s) {
                            return Some(s);
                        }
                    }
                    i = end;
                    continue;
                }
            }
        }
        i += 1;
    }
    None
}

fn read_ber_len(buf: &[u8], pos: usize) -> Option<(usize, usize)> {
    let first = *buf.get(pos)?;
    if first & 0x80 == 0 {
        return Some((first as usize, 1));
    }
    let count = (first & 0x7f) as usize;
    if count == 0 || count > std::mem::size_of::<usize>() || pos + 1 + count > buf.len() {
        return None;
    }
    let mut len = 0usize;
    for &b in &buf[pos + 1..pos + 1 + count] {
        len = (len << 8) | b as usize;
    }
    Some((len, 1 + count))
}

fn looks_like_header_title(s: &str) -> bool {
    if s.is_empty() {
        return false;
    }
    let has_word_separator = s.bytes().any(|b| b == b' ' || b == b'\t');
    let has_lowercase = s.bytes().any(|b| b.is_ascii_lowercase());
    let seqid_only = s
        .bytes()
        .all(|b| b.is_ascii_alphanumeric() || matches!(b, b'_' | b'.' | b'-' | b'|'));
    (has_word_separator || has_lowercase) && !seqid_only
}

fn write_pairwise_subject_report_preamble<W: Write>(
    writer: &mut W,
    subjects: &[blast_rs::input::FastaRecord],
    args: &BlastnArgs,
    total_subject_len: i64,
) -> io::Result<()> {
    writeln!(writer, "BLASTN 2.12.0+")?;
    writeln!(writer)?;
    writeln!(writer)?;
    writeln!(
        writer,
        "Reference: Stephen F. Altschul, Thomas L. Madden, Alejandro A."
    )?;
    writeln!(
        writer,
        "Schaffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J."
    )?;
    writeln!(
        writer,
        "Lipman (1997), \"Gapped BLAST and PSI-BLAST: a new generation of"
    )?;
    writeln!(
        writer,
        "protein database search programs\", Nucleic Acids Res. 25:3389-3402."
    )?;
    writeln!(writer)?;
    writeln!(writer)?;
    writeln!(writer)?;
    write_pairwise_subject_database_line(writer, "", args.subject.as_ref())?;
    writeln!(
        writer,
        "           {} sequences; {} total letters",
        format_with_commas(subjects.len() as u64),
        format_with_commas(total_subject_len.max(0) as u64),
    )?;
    writeln!(writer)?;
    writeln!(writer)?;
    writeln!(writer)?;
    Ok(())
}

fn write_pairwise_subject_query_header<W: Write>(
    writer: &mut W,
    query: &blast_rs::input::FastaRecord,
    subjects: &[blast_rs::input::FastaRecord],
    hits: &[&TabularHit],
    args: &BlastnArgs,
) -> io::Result<()> {
    let subject_deflines: std::collections::HashMap<&str, &str> = subjects
        .iter()
        .map(|rec| (rec.id.as_str(), rec.defline.as_str()))
        .collect();
    writeln!(writer, "Query= {}", query.defline.as_str())?;
    writeln!(writer)?;
    writeln!(writer, "Length={}", query.sequence.len())?;
    if hits.is_empty() {
        writeln!(writer)?;
        writeln!(writer)?;
        writeln!(writer, "***** No hits found *****")?;
        writeln!(writer)?;
        return Ok(());
    }
    write_pairwise_hit_summary_header(writer, args.sorthits)?;
    writeln!(writer)?;

    let mut seen = std::collections::HashSet::new();
    for hit in hits {
        if !seen.insert(hit.subject_id.as_str()) {
            continue;
        }
        let desc = truncate_description(
            subject_deflines
                .get(hit.subject_id.as_str())
                .copied()
                .unwrap_or(hit.subject_id.as_str()),
            68,
        );
        let subject_hits: Vec<&TabularHit> = hits
            .iter()
            .copied()
            .filter(|h| h.subject_id == hit.subject_id)
            .collect();
        write_pairwise_hit_summary_row(writer, &desc, &subject_hits, args.sorthits)?;
    }
    writeln!(writer)?;
    writeln!(writer)?;
    Ok(())
}

fn write_pairwise_subject_query_stats<W: Write>(
    writer: &mut W,
    args: &BlastnArgs,
    query: &blast_rs::input::FastaRecord,
    total_subject_len: i64,
    num_subjects: i32,
) -> io::Result<()> {
    writeln!(writer)?;
    writeln!(writer)?;
    writeln!(writer, "Lambda      K        H")?;
    writeln!(writer, "    1.37    0.711     1.31 ")?;
    writeln!(writer)?;
    writeln!(writer, "Gapped")?;
    writeln!(writer, "Lambda      K        H")?;
    writeln!(writer, "    1.37    0.711     1.31 ")?;
    writeln!(writer)?;
    writeln!(
        writer,
        "Effective search space used: {}",
        pairwise_subject_search_space(args, query, total_subject_len, num_subjects)
    )?;
    writeln!(writer)?;
    writeln!(writer)?;
    Ok(())
}

fn write_pairwise_subject_database_footer<W: Write>(
    writer: &mut W,
    subjects: &[blast_rs::input::FastaRecord],
    args: &BlastnArgs,
    total_subject_len: i64,
) -> io::Result<()> {
    write_pairwise_subject_database_line(writer, "  ", args.subject.as_ref())?;
    writeln!(writer, "    Posted date:  Unknown")?;
    writeln!(
        writer,
        "  Number of letters in database: {}",
        format_with_commas(total_subject_len.max(0) as u64)
    )?;
    writeln!(
        writer,
        "  Number of sequences in database:  {}",
        format_with_commas(subjects.len() as u64)
    )?;
    writeln!(writer)?;
    writeln!(writer)?;
    writeln!(writer)?;
    writeln!(
        writer,
        "Matrix: blastn matrix {} {}",
        args.reward(),
        args.penalty()
    )?;
    writeln!(
        writer,
        "Gap Penalties: Existence: {}, Extension: {}",
        args.gapopen(),
        args.gapextend()
    )
}

fn write_pairwise_subject_database_line<W: Write>(
    writer: &mut W,
    prefix: &str,
    subject: Option<&PathBuf>,
) -> io::Result<()> {
    let Some(subject) = subject else {
        return writeln!(writer, "{}Database: User specified sequence set.", prefix);
    };
    let path = subject.display().to_string();
    let inline = format!(
        "{}Database: User specified sequence set (Input: {}).",
        prefix, path
    );
    if inline.len() > 79 {
        writeln!(
            writer,
            "{}Database: User specified sequence set (Input:",
            prefix
        )?;
        writeln!(writer, "{}).", path)
    } else {
        writeln!(writer, "{}", inline)
    }
}

fn pairwise_subject_search_space(
    args: &BlastnArgs,
    query: &blast_rs::input::FastaRecord,
    total_subject_len: i64,
    num_subjects: i32,
) -> u64 {
    let Ok((loc_start, loc_end)) = query_loc_bounds(args, query.sequence.len()) else {
        return 0;
    };
    let query_plus_nomask: Vec<u8> = query.sequence[loc_start..loc_end]
        .iter()
        .map(|&b| iupacna_to_blastna(b))
        .collect();
    let (_, search_space, _) =
        blastn_subject_stats(args, &query_plus_nomask, total_subject_len, num_subjects);
    search_space.round() as u64
}

fn write_pairwise_db_report_preamble<W: Write>(writer: &mut W, db: &BlastDb) -> io::Result<()> {
    writeln!(writer, "BLASTN 2.12.0+")?;
    writeln!(writer)?;
    writeln!(writer)?;
    writeln!(
        writer,
        "Reference: Stephen F. Altschul, Thomas L. Madden, Alejandro A."
    )?;
    writeln!(
        writer,
        "Schaffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J."
    )?;
    writeln!(
        writer,
        "Lipman (1997), \"Gapped BLAST and PSI-BLAST: a new generation of"
    )?;
    writeln!(
        writer,
        "protein database search programs\", Nucleic Acids Res. 25:3389-3402."
    )?;
    writeln!(writer)?;
    writeln!(writer)?;
    writeln!(writer)?;
    writeln!(writer, "Database: {}", db.title)?;
    writeln!(
        writer,
        "           {} sequences; {} total letters",
        format_with_commas(db.num_oids as u64),
        format_with_commas(db.total_length),
    )?;
    writeln!(writer)?;
    writeln!(writer)?;
    writeln!(writer)?;
    Ok(())
}

fn write_pairwise_db_query_header<W: Write>(
    writer: &mut W,
    query: &blast_rs::input::FastaRecord,
    hits: &[&TabularHit],
    subject_deflines: &std::collections::HashMap<String, String>,
    args: &BlastnArgs,
) -> io::Result<()> {
    writeln!(writer, "Query= {}", query.defline.as_str())?;
    writeln!(writer)?;
    writeln!(writer, "Length={}", query.sequence.len())?;
    if hits.is_empty() {
        writeln!(writer)?;
        writeln!(writer)?;
        writeln!(writer, "***** No hits found *****")?;
        writeln!(writer)?;
        return Ok(());
    }
    write_pairwise_hit_summary_header(writer, args.sorthits)?;
    writeln!(writer)?;

    let mut seen = std::collections::HashSet::new();
    for hit in hits {
        if !seen.insert(hit.subject_id.as_str()) {
            continue;
        }
        let defline = subject_deflines
            .get(&hit.subject_id)
            .map(String::as_str)
            .unwrap_or(hit.subject_id.as_str());
        let desc = truncate_description(defline, 68);
        let subject_hits: Vec<&TabularHit> = hits
            .iter()
            .copied()
            .filter(|h| h.subject_id == hit.subject_id)
            .collect();
        write_pairwise_hit_summary_row(writer, &desc, &subject_hits, args.sorthits)?;
    }
    writeln!(writer)?;
    writeln!(writer)?;
    Ok(())
}

fn write_pairwise_db_query_stats<W: Write>(writer: &mut W, search_space: f64) -> io::Result<()> {
    writeln!(writer)?;
    writeln!(writer)?;
    writeln!(writer, "Lambda      K        H")?;
    writeln!(writer, "    1.37    0.711     1.31 ")?;
    writeln!(writer)?;
    writeln!(writer, "Gapped")?;
    writeln!(writer, "Lambda      K        H")?;
    writeln!(writer, "    1.37    0.711     1.31 ")?;
    writeln!(writer)?;
    writeln!(
        writer,
        "Effective search space used: {}",
        search_space.round() as u64
    )?;
    writeln!(writer)?;
    writeln!(writer)?;
    Ok(())
}

fn write_pairwise_db_database_footer<W: Write>(
    writer: &mut W,
    db: &BlastDb,
    args: &BlastnArgs,
) -> io::Result<()> {
    writeln!(writer, "  Database: {}", db.title)?;
    writeln!(writer, "    Posted date:  {}", db.date)?;
    writeln!(
        writer,
        "  Number of letters in database: {}",
        format_with_commas(db.total_length)
    )?;
    writeln!(
        writer,
        "  Number of sequences in database:  {}",
        format_with_commas(db.num_oids as u64)
    )?;
    writeln!(writer)?;
    writeln!(writer)?;
    writeln!(writer)?;
    writeln!(
        writer,
        "Matrix: blastn matrix {} {}",
        args.reward(),
        args.penalty()
    )?;
    writeln!(
        writer,
        "Gap Penalties: Existence: {}, Extension: {}",
        args.gapopen(),
        args.gapextend()
    )
}

fn truncate_description(s: &str, width: usize) -> String {
    if s.len() <= width {
        s.to_string()
    } else if width <= 3 {
        ".".repeat(width)
    } else {
        format!("{}...", &s[..width - 3])
    }
}

fn format_with_commas(n: u64) -> String {
    let s = n.to_string();
    let mut out = String::with_capacity(s.len() + s.len() / 3);
    for (i, b) in s.bytes().enumerate() {
        if i > 0 && (s.len() - i) % 3 == 0 {
            out.push(',');
        }
        out.push(b as char);
    }
    out
}

fn alignment_string_to_blastna(seq: &str) -> Vec<u8> {
    seq.bytes()
        .map(|b| match b {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' | b'U' | b'u' => 3,
            b'R' | b'r' => 4,
            b'Y' | b'y' => 5,
            b'M' | b'm' => 6,
            b'K' | b'k' => 7,
            b'W' | b'w' => 8,
            b'S' | b's' => 9,
            b'B' | b'b' => 10,
            b'D' | b'd' => 11,
            b'H' | b'h' => 12,
            b'V' | b'v' => 13,
            b'N' | b'n' => 14,
            _ => 15,
        })
        .collect()
}

fn write_pairwise_alignment<W: Write>(
    writer: &mut W,
    hit: &TabularHit,
    query_aln: &[u8],
    subject_aln: &[u8],
    show_subject_header: bool,
    subject_display: Option<&str>,
) -> io::Result<()> {
    let subject_display = subject_display.unwrap_or(hit.subject_id.as_str());
    let mut buf = Vec::new();
    blast_rs::format::format_pairwise_alignment(
        &mut buf,
        &hit.query_id,
        subject_display,
        query_aln,
        subject_aln,
        hit.query_start,
        hit.query_end,
        hit.subject_start,
        hit.subject_end,
        hit.raw_score,
        hit.bit_score,
        hit.evalue,
        hit.num_ident,
        hit.align_len,
        hit.gap_opens,
    )?;
    let start = buf.windows(2).position(|w| w == b"\n\n").map_or_else(
        || {
            buf.iter()
                .position(|&b| b == b'\n')
                .map_or(0, |idx| idx + 1)
        },
        |idx| idx + 1,
    );
    if show_subject_header {
        writer.write_all(&buf[..start])?;
        writeln!(writer, "Length={}", hit.subject_len)?;
    }
    writer.write_all(&buf[start..])
}

fn reverse_complement_alignment_blastna(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement_blastna(b)).collect()
}

/// DEPRECATED: No longer used. All KBP computation is now pure Rust.
/// Kept only as reference for the C FFI calling pattern.

/// Apply post-search filters (perc_identity, qcov_hsp_perc, max_hsps).
fn apply_filters(hits: &mut Vec<TabularHit>, args: &BlastnArgs, _query_len: i32) {
    let taxids = parse_taxid_filters(args.taxids.as_deref(), args.taxidlist.as_ref());
    if !taxids.is_empty() {
        hits.retain(|h| h.subject_taxids.iter().any(|taxid| taxids.contains(taxid)));
    }
    let negative_taxids = parse_taxid_filters(
        args.negative_taxids.as_deref(),
        args.negative_taxidlist.as_ref(),
    );
    if !negative_taxids.is_empty() {
        hits.retain(|h| {
            !h.subject_taxids
                .iter()
                .any(|taxid| negative_taxids.contains(taxid))
        });
    }
    // Filter by percent identity
    if args.perc_identity > 0.0 {
        hits.retain(|h| h.pct_identity >= args.perc_identity);
    }
    if let Some(min_score) = args.min_raw_gapped_score {
        hits.retain(|h| h.raw_score >= min_score);
    }
    // Filter by query coverage
    if args.qcov_hsp_perc > 0.0 {
        hits.retain(|h| {
            if h.query_len <= 0 {
                return false;
            }
            let query_span = (h.query_end - h.query_start).abs() + 1;
            let cov = 100.0 * query_span as f64 / h.query_len as f64;
            cov >= args.qcov_hsp_perc
        });
    }
    if args.culling_limit > 0 {
        apply_culling_limit(hits, args.culling_limit as usize);
    }
    if args.subject_besthit {
        apply_subject_besthit_filter(hits);
    }
    if args.best_hit_overhang.is_some() || args.best_hit_score_edge.is_some() {
        apply_best_hit_filter(
            hits,
            args.best_hit_overhang.unwrap_or(0.1),
            args.best_hit_score_edge.unwrap_or(0.1),
        );
    }
    // Limit HSPs per subject
    if let Some(max_hsps) = args.max_hsps {
        let max = max_hsps as usize;
        let mut counts: std::collections::HashMap<(String, String), usize> =
            std::collections::HashMap::new();
        hits.retain(|h| {
            let c = counts
                .entry((h.query_id.clone(), h.subject_id.clone()))
                .or_insert(0);
            *c += 1;
            *c <= max
        });
    }
}

fn apply_culling_limit(hits: &mut Vec<TabularHit>, culling_limit: usize) {
    if culling_limit == 0 || hits.len() <= 1 {
        return;
    }

    let mut kept: Vec<TabularHit> = Vec::with_capacity(hits.len());
    for hit in hits.drain(..) {
        let q_start = hit.query_start.min(hit.query_end);
        let q_end = hit.query_start.max(hit.query_end);
        let enveloping = kept
            .iter()
            .filter(|existing| {
                if existing.query_id != hit.query_id || existing.raw_score < hit.raw_score {
                    return false;
                }
                let existing_q_start = existing.query_start.min(existing.query_end);
                let existing_q_end = existing.query_start.max(existing.query_end);
                existing_q_start <= q_start && existing_q_end >= q_end
            })
            .take(culling_limit)
            .count();
        if enveloping < culling_limit {
            kept.push(hit);
        }
    }
    *hits = kept;
}

fn apply_subject_besthit_filter(hits: &mut Vec<TabularHit>) {
    if hits.len() <= 1 {
        return;
    }

    const SUBJECT_BESTHIT_NUCL_MAX_RANGE_DIFF: i32 = 3;
    let mut groups: std::collections::BTreeMap<(String, String), Vec<TabularHit>> =
        std::collections::BTreeMap::new();
    for hit in hits.drain(..) {
        groups
            .entry((hit.query_id.clone(), hit.subject_id.clone()))
            .or_default()
            .push(hit);
    }

    for group in groups.values_mut() {
        let mut keep = vec![true; group.len()];
        for i in 0..group.len() {
            if !keep[i] {
                continue;
            }
            let (offset, end) = hsp_context_query_range(&group[i]);
            let offset = (offset - SUBJECT_BESTHIT_NUCL_MAX_RANGE_DIFF).max(0);
            let end = end + SUBJECT_BESTHIT_NUCL_MAX_RANGE_DIFF;
            for j in (i + 1)..group.len() {
                if keep[j] && group[j].sframe == group[i].sframe {
                    let (other_offset, other_end) = hsp_context_query_range(&group[j]);
                    if other_offset >= offset && other_end <= end {
                        keep[j] = false;
                    }
                }
            }
        }

        for i in 0..group.len() {
            if !keep[i] {
                continue;
            }
            let (offset, end) = hsp_context_query_range(&group[i]);
            let offset = (offset - SUBJECT_BESTHIT_NUCL_MAX_RANGE_DIFF).max(0);
            let end = end + SUBJECT_BESTHIT_NUCL_MAX_RANGE_DIFF;
            let flipped_offset = group[i].query_len - end;
            let flipped_end = group[i].query_len - offset;
            for j in (i + 1)..group.len() {
                if keep[j] && group[j].sframe == -group[i].sframe {
                    let (other_offset, other_end) = hsp_context_query_range(&group[j]);
                    if other_offset >= flipped_offset && other_end <= flipped_end {
                        keep[j] = false;
                    }
                }
            }
        }

        let mut idx = 0;
        group.retain(|_| {
            let should_keep = keep[idx];
            idx += 1;
            should_keep
        });
    }

    hits.extend(groups.into_values().flatten());
}

fn hsp_context_query_range(hit: &TabularHit) -> (i32, i32) {
    if hit.sframe < 0 {
        (
            hit.query_len - hit.query_start.max(hit.query_end),
            hit.query_len - hit.query_start.min(hit.query_end) + 1,
        )
    } else {
        (
            hit.query_start.min(hit.query_end) - 1,
            hit.query_start.max(hit.query_end),
        )
    }
}

#[derive(Clone)]
struct BestHitNode {
    hit: TabularHit,
    begin: i32,
    end: i32,
    len: i32,
}

fn apply_best_hit_filter(hits: &mut Vec<TabularHit>, overhang: f64, score_edge: f64) {
    if hits.len() <= 1 {
        return;
    }

    let overhang = overhang.clamp(0.0, 0.499_999);
    let score_edge = score_edge.clamp(0.0, 0.5);
    let param_s = 1.0 - score_edge;
    let mut by_query: std::collections::BTreeMap<String, Vec<BestHitNode>> =
        std::collections::BTreeMap::new();

    for hit in hits.drain(..) {
        let query_id = hit.query_id.clone();
        let nodes = by_query.entry(query_id).or_default();
        let begin = hsp_query_order_start(&hit) - 1;
        let len = (hit.query_end - hit.query_start).abs() + 1;
        if len <= 0 {
            continue;
        }
        let end = begin + len;
        let score = hit.raw_score;
        let evalue = hit.evalue;
        let new_bad_density = score as f64 / len as f64 / param_s;

        let is_bad = nodes.iter().any(|node| {
            node.end >= end
                && node.begin <= begin
                && node.hit.evalue <= evalue
                && node.hit.raw_score as f64 / node.len as f64 > new_bad_density
        });
        if is_bad {
            continue;
        }

        let wide_overhang = (2.0 * len as f64 * overhang / (1.0 - 2.0 * overhang)) as i32;
        let allowed_begin = begin - wide_overhang;
        let allowed_end = end + wide_overhang;
        let stored_overhang = (len as f64 * overhang) as i32;
        let stored_begin = begin - stored_overhang;
        let stored_end = end + stored_overhang;
        let old_bad_density = score as f64 / len as f64 * param_s;

        nodes.retain(|node| {
            if node.begin < allowed_begin || node.begin >= allowed_end {
                return true;
            }
            let node_overhang = (node.end - node.begin - node.len) / 2;
            !(node.begin + node_overhang >= stored_begin
                && node.end - node_overhang <= stored_end
                && node.hit.evalue >= evalue
                && node.hit.raw_score as f64 / (node.len as f64) < old_bad_density)
        });

        let insert_at = nodes
            .iter()
            .position(|node| node.begin >= stored_begin)
            .unwrap_or(nodes.len());
        nodes.insert(
            insert_at,
            BestHitNode {
                hit,
                begin: stored_begin,
                end: stored_end,
                len,
            },
        );
    }

    hits.extend(
        by_query
            .into_values()
            .flat_map(|nodes| nodes.into_iter().map(|node| node.hit)),
    );
}

fn parse_taxid_filters(
    value: Option<&str>,
    list_path: Option<&PathBuf>,
) -> std::collections::HashSet<i32> {
    let mut taxids: std::collections::HashSet<i32> = value
        .into_iter()
        .flat_map(|s| s.split(','))
        .filter_map(|s| s.trim().parse::<i32>().ok())
        .collect();
    if let Some(path) = list_path {
        if let Ok(contents) = std::fs::read_to_string(path) {
            taxids.extend(
                contents
                    .split(|ch: char| ch == ',' || ch.is_ascii_whitespace())
                    .filter_map(|s| s.trim().parse::<i32>().ok()),
            );
        }
    }
    taxids
}

/// Complement a BLASTNA-encoded nucleotide.
fn complement_blastna(b: u8) -> u8 {
    match b {
        0 => 3, // A -> T
        1 => 2, // C -> G
        2 => 1, // G -> C
        3 => 0, // T -> A
        // Ambiguity codes: complement by swapping bits
        4 => 5,   // R(AG) -> Y(CT)
        5 => 4,   // Y(CT) -> R(AG)
        6 => 7,   // M(AC) -> K(GT)
        7 => 6,   // K(GT) -> M(AC)
        8 => 8,   // W(AT) -> W(AT)
        9 => 9,   // S(CG) -> S(CG)
        10 => 13, // B(CGT) -> V(ACG)
        11 => 12, // D(AGT) -> H(ACT)
        12 => 11, // H(ACT) -> D(AGT)
        13 => 10, // V(ACG) -> B(CGT)
        14 => 14, // N -> N
        _ => 15,  // gap
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_complement() {
        assert_eq!(complement_blastna(0), 3); // A -> T
        assert_eq!(complement_blastna(3), 0); // T -> A
        assert_eq!(complement_blastna(1), 2); // C -> G
        assert_eq!(complement_blastna(14), 14); // N -> N
    }

    #[test]
    fn test_blastn_omitted_task_defaults_to_megablast() {
        let cli = Cli::parse_from([
            "blast-cli",
            "blastn",
            "--query",
            "tests/fixtures/query_random_200.fa",
            "--db",
            "tests/fixtures/large_db/celegans",
        ]);
        let Commands::Blastn(mut args) = cli.command else {
            panic!("expected blastn command");
        };

        args.apply_task_defaults();

        assert_eq!(args.word_size(), 28);
        assert_eq!(args.reward(), 1);
        assert_eq!(args.penalty(), -2);
        assert_eq!(args.gapopen(), 0);
        assert_eq!(args.gapextend(), 0);
    }

    #[test]
    fn test_explicit_blastn_task_defaults_are_preserved() {
        let cli = Cli::parse_from([
            "blast-cli",
            "blastn",
            "--query",
            "tests/fixtures/query_random_200.fa",
            "--db",
            "tests/fixtures/large_db/celegans",
            "--task",
            "blastn",
        ]);
        let Commands::Blastn(mut args) = cli.command else {
            panic!("expected blastn command");
        };

        args.apply_task_defaults();

        assert_eq!(args.word_size(), 11);
        assert_eq!(args.reward(), 2);
        assert_eq!(args.penalty(), -3);
        assert_eq!(args.gapopen(), 5);
        assert_eq!(args.gapextend(), 2);
    }

    #[test]
    fn test_dc_megablast_task_defaults_match_ncbi() {
        let cli = Cli::parse_from([
            "blast-cli",
            "blastn",
            "--query",
            "tests/fixtures/query_random_200.fa",
            "--db",
            "tests/fixtures/large_db/celegans",
            "--task",
            "dc-megablast",
        ]);
        let Commands::Blastn(mut args) = cli.command else {
            panic!("expected blastn command");
        };

        args.apply_task_defaults();

        assert_eq!(args.word_size(), 28);
        assert_eq!(args.reward(), 2);
        assert_eq!(args.penalty(), -3);
        assert_eq!(args.gapopen(), 5);
        assert_eq!(args.gapextend(), 2);
    }

    #[test]
    fn test_ncbi_underscore_option_aliases_parse() {
        let cli = Cli::parse_from([
            "blast-cli",
            "blastn",
            "--query",
            "tests/fixtures/query_random_200.fa",
            "--db",
            "tests/fixtures/large_db/celegans",
            "--best_hit_overhang",
            "0.1",
            "--best_hit_score_edge",
            "0.1",
            "--subject_besthit",
            "--parse_deflines",
            "--negative_gilist",
            "tests/fixtures/gilist.txt",
            "--negative_seqidlist",
            "tests/fixtures/seqidlist.txt",
            "--template_type",
            "coding",
            "--template_length",
            "18",
            "--entrez_query",
            "txid9606[ORGN]",
            "--use_index",
            "--index_name",
            "nt",
            "--mt_mode",
            "1",
            "--off_diagonal_range",
            "2",
            "--sort_hits",
            "2",
            "--sort_hsps",
            "3",
            "--no_taxid_expansion",
        ]);
        let Commands::Blastn(args) = cli.command else {
            panic!("expected blastn command");
        };

        assert_eq!(args.best_hit_overhang, Some(0.1));
        assert_eq!(args.best_hit_score_edge, Some(0.1));
        assert!(args.subject_besthit);
        assert!(args.parse_deflines);
        assert!(args.negative_gilist.is_some());
        assert!(args.negative_seqidlist.is_some());
        assert_eq!(args.template_type.as_deref(), Some("coding"));
        assert_eq!(args.template_length, Some(18));
        assert_eq!(args.entrez_query.as_deref(), Some("txid9606[ORGN]"));
        assert!(args.use_index);
        assert_eq!(args.index_name.as_deref(), Some("nt"));
        assert_eq!(args.mt_mode, 1);
        assert_eq!(args.off_diagonal_range, 2);
        assert_eq!(args.sorthits, 2);
        assert_eq!(args.sorthsps, 3);
        assert!(args.no_taxid_expansion);
    }

    #[test]
    fn test_sort_option_ranges_match_ncbi() {
        assert!(Cli::try_parse_from([
            "blast-cli",
            "blastn",
            "--query",
            "tests/fixtures/query_random_200.fa",
            "--db",
            "tests/fixtures/large_db/celegans",
            "--sorthits",
            "5",
        ])
        .is_err());
        assert!(Cli::try_parse_from([
            "blast-cli",
            "blastn",
            "--query",
            "tests/fixtures/query_random_200.fa",
            "--db",
            "tests/fixtures/large_db/celegans",
            "--sorthsps",
            "-1",
        ])
        .is_err());
        assert!(Cli::try_parse_from([
            "blast-cli",
            "blastn",
            "--query",
            "tests/fixtures/query_random_200.fa",
            "--db",
            "tests/fixtures/large_db/celegans",
            "--max_target_seqs",
            "0",
        ])
        .is_err());
        assert!(Cli::try_parse_from([
            "blast-cli",
            "blastn",
            "--query",
            "tests/fixtures/query_random_200.fa",
            "--db",
            "tests/fixtures/large_db/celegans",
            "--max_target_seqs",
            "10",
            "--num_alignments",
            "10",
        ])
        .is_err());
        assert!(Cli::try_parse_from([
            "blast-cli",
            "blastn",
            "--query",
            "tests/fixtures/query_random_200.fa",
            "--db",
            "tests/fixtures/large_db/celegans",
            "--max_target_seqs",
            "10",
            "--num_descriptions",
            "10",
        ])
        .is_err());
        assert!(Cli::try_parse_from([
            "blast-cli",
            "blastn",
            "--query",
            "tests/fixtures/query_random_200.fa",
            "--db",
            "tests/fixtures/large_db/celegans",
            "--num_threads",
            "0",
        ])
        .is_err());
        assert!(Cli::try_parse_from([
            "blast-cli",
            "blastn",
            "--query",
            "tests/fixtures/query_random_200.fa",
            "--db",
            "tests/fixtures/large_db/celegans",
            "--max_hsps",
            "0",
        ])
        .is_err());
        assert!(Cli::try_parse_from([
            "blast-cli",
            "blastn",
            "--query",
            "tests/fixtures/query_random_200.fa",
            "--db",
            "tests/fixtures/large_db/celegans",
            "--perc_identity",
            "101",
        ])
        .is_err());
        assert!(Cli::try_parse_from([
            "blast-cli",
            "blastn",
            "--query",
            "tests/fixtures/query_random_200.fa",
            "--db",
            "tests/fixtures/large_db/celegans",
            "--qcov_hsp_perc",
            "-1",
        ])
        .is_err());
        assert!(Cli::try_parse_from([
            "blast-cli",
            "blastn",
            "--query",
            "tests/fixtures/query_random_200.fa",
            "--db",
            "tests/fixtures/large_db/celegans",
            "--culling_limit",
            "-1",
        ])
        .is_err());
        assert!(Cli::try_parse_from([
            "blast-cli",
            "blastn",
            "--query",
            "tests/fixtures/query_random_200.fa",
            "--db",
            "tests/fixtures/large_db/celegans",
            "--num_descriptions",
            "-1",
        ])
        .is_err());
        assert!(Cli::try_parse_from([
            "blast-cli",
            "blastn",
            "--query",
            "tests/fixtures/query_random_200.fa",
            "--db",
            "tests/fixtures/large_db/celegans",
            "--num_alignments",
            "-1",
        ])
        .is_err());
        assert!(Cli::try_parse_from([
            "blast-cli",
            "blastn",
            "--query",
            "tests/fixtures/query_random_200.fa",
            "--db",
            "tests/fixtures/large_db/celegans",
            "--line_length",
            "0",
        ])
        .is_err());
        assert!(Cli::try_parse_from([
            "blast-cli",
            "blastn",
            "--query",
            "tests/fixtures/query_random_200.fa",
            "--db",
            "tests/fixtures/large_db/celegans",
            "--window_size",
            "-1",
        ])
        .is_err());
        assert!(Cli::try_parse_from([
            "blast-cli",
            "blastn",
            "--query",
            "tests/fixtures/query_random_200.fa",
            "--db",
            "tests/fixtures/large_db/celegans",
            "--off_diagonal_range",
            "-1",
        ])
        .is_err());
        assert!(Cli::try_parse_from([
            "blast-cli",
            "blastn",
            "--query",
            "tests/fixtures/query_random_200.fa",
            "--db",
            "tests/fixtures/large_db/celegans",
            "--mt_mode",
            "2",
        ])
        .is_err());
        assert!(Cli::try_parse_from([
            "blast-cli",
            "blastn",
            "--query",
            "tests/fixtures/query_random_200.fa",
            "--db",
            "tests/fixtures/large_db/celegans",
            "--best_hit_overhang",
            "0.5",
        ])
        .is_err());
        assert!(Cli::try_parse_from([
            "blast-cli",
            "blastn",
            "--query",
            "tests/fixtures/query_random_200.fa",
            "--db",
            "tests/fixtures/large_db/celegans",
            "--best_hit_score_edge",
            "0",
        ])
        .is_err());
        assert!(Cli::try_parse_from([
            "blast-cli",
            "blastn",
            "--query",
            "tests/fixtures/query_random_200.fa",
            "--db",
            "tests/fixtures/large_db/celegans",
            "--task",
            "bogus",
        ])
        .is_err());
        assert!(Cli::try_parse_from([
            "blast-cli",
            "blastn",
            "--query",
            "tests/fixtures/query_random_200.fa",
            "--db",
            "tests/fixtures/large_db/celegans",
            "--strand",
            "bogus",
        ])
        .is_err());
        assert!(Cli::try_parse_from([
            "blast-cli",
            "blastn",
            "--query",
            "tests/fixtures/query_random_200.fa",
            "--db",
            "tests/fixtures/large_db/celegans",
            "--template_type",
            "bogus",
        ])
        .is_err());
        assert!(Cli::try_parse_from([
            "blast-cli",
            "blastn",
            "--query",
            "tests/fixtures/query_random_200.fa",
            "--db",
            "tests/fixtures/large_db/celegans",
            "--template_length",
            "17",
        ])
        .is_err());
        assert!(Cli::try_parse_from([
            "blast-cli",
            "blastn",
            "--query",
            "tests/fixtures/query_random_200.fa",
            "--db",
            "tests/fixtures/large_db/celegans",
            "--reward",
            "-1",
        ])
        .is_err());
        assert!(Cli::try_parse_from([
            "blast-cli",
            "blastn",
            "--query",
            "tests/fixtures/query_random_200.fa",
            "--db",
            "tests/fixtures/large_db/celegans",
            "--penalty",
            "1",
        ])
        .is_err());
        assert!(Cli::try_parse_from([
            "blast-cli",
            "blastn",
            "--query",
            "tests/fixtures/query_random_200.fa",
            "--db",
            "tests/fixtures/large_db/celegans",
            "--searchsp",
            "-1",
        ])
        .is_err());
    }

    fn tabular_hit_for_best_hit_filter(
        subject_id: &str,
        query_start: i32,
        query_end: i32,
        raw_score: i32,
        evalue: f64,
    ) -> TabularHit {
        let align_len = (query_end - query_start).abs() + 1;
        TabularHit {
            query_id: "q1".to_string(),
            subject_id: subject_id.to_string(),
            subject_title: String::new(),
            pct_identity: 100.0,
            align_len,
            mismatches: 0,
            gap_opens: 0,
            query_start,
            query_end,
            subject_start: 1,
            subject_end: align_len,
            evalue,
            bit_score: raw_score as f64,
            query_len: 40,
            subject_len: align_len,
            raw_score,
            qseq: None,
            sseq: None,
            qframe: 1,
            sframe: if query_start > query_end { -1 } else { 1 },
            subject_taxids: vec![],
            subject_sci_name: String::new(),
            subject_common_name: String::new(),
            subject_blast_name: String::new(),
            subject_kingdom: String::new(),
            num_ident: align_len,
        }
    }

    #[test]
    fn test_best_hit_filter_drops_lower_density_contained_hsp() {
        let mut hits = vec![
            tabular_hit_for_best_hit_filter("full", 1, 40, 80, 1.0e-30),
            tabular_hit_for_best_hit_filter("contained", 5, 32, 42, 1.0e-12),
        ];

        apply_best_hit_filter(&mut hits, 0.1, 0.1);

        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].subject_id, "full");
    }

    #[test]
    fn test_best_hit_filter_keeps_equal_density_contained_hsp() {
        let mut hits = vec![
            tabular_hit_for_best_hit_filter("full", 1, 40, 80, 1.0e-30),
            tabular_hit_for_best_hit_filter("contained", 5, 30, 52, 1.0e-20),
        ];

        apply_best_hit_filter(&mut hits, 0.1, 0.1);

        assert_eq!(hits.len(), 2);
    }

    #[test]
    fn test_subject_besthit_filter_removes_same_region_opposite_strand() {
        let mut reverse = tabular_hit_for_best_hit_filter("s1", 1, 40, 80, 1.0e-30);
        reverse.subject_start = 40;
        reverse.subject_end = 1;
        reverse.sframe = -1;
        let mut contained = tabular_hit_for_best_hit_filter("s1", 5, 32, 64, 1.0e-20);
        contained.subject_start = 1;
        contained.subject_end = 28;

        let mut hits = vec![
            tabular_hit_for_best_hit_filter("s1", 1, 40, 80, 1.0e-30),
            reverse,
            contained,
        ];

        apply_subject_besthit_filter(&mut hits);

        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].query_start, 1);
        assert_eq!(hits[0].query_end, 40);
        assert_eq!(hits[0].sframe, 1);
    }

    fn hsp_for_hitlist(score: i32, evalue: f64) -> blast_rs::search::SearchHsp {
        blast_rs::search::SearchHsp {
            query_start: 0,
            query_end: 28,
            subject_start: 0,
            subject_end: 28,
            score,
            bit_score: score as f64,
            evalue,
            num_ident: 28,
            align_length: 28,
            mismatches: 0,
            gap_opens: 0,
            context: 0,
            qseq: None,
            sseq: None,
        }
    }

    #[test]
    fn test_hitlist_prelim_pruning_prefers_high_oid_on_exact_ties() {
        let low = vec![hsp_for_hitlist(28, 1.0e-13)];
        let high = vec![hsp_for_hitlist(28, 1.0e-13)];

        assert_eq!(
            compare_oid_hsps_for_hitlist(5, &low, 19, &high),
            std::cmp::Ordering::Greater
        );
        assert_eq!(
            compare_oid_hsps_for_hitlist(19, &high, 5, &low),
            std::cmp::Ordering::Less
        );
    }

    #[test]
    fn test_hitlist_prelim_pruning_prefers_evalue_then_score() {
        let worse_evalue = vec![hsp_for_hitlist(100, 1.0e-10)];
        let better_evalue = vec![hsp_for_hitlist(50, 1.0e-20)];
        assert_eq!(
            compare_oid_hsps_for_hitlist(100, &worse_evalue, 1, &better_evalue),
            std::cmp::Ordering::Greater
        );

        let lower_score = vec![hsp_for_hitlist(40, 1.0e-20)];
        let higher_score = vec![hsp_for_hitlist(80, 1.0e-20)];
        assert_eq!(
            compare_oid_hsps_for_hitlist(100, &lower_score, 1, &higher_score),
            std::cmp::Ordering::Greater
        );
    }

    #[test]
    fn test_alignment_string_to_blastna_handles_gaps_and_ambiguity() {
        assert_eq!(
            alignment_string_to_blastna("ACGT-RYN"),
            vec![0, 1, 2, 3, 15, 4, 5, 14]
        );
    }

    #[test]
    fn test_reverse_complement_alignment_blastna_preserves_gaps() {
        assert_eq!(
            reverse_complement_alignment_blastna(&alignment_string_to_blastna("ACGT-RYN")),
            alignment_string_to_blastna("NRY-ACGT")
        );
    }
}
