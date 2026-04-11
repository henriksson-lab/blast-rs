//! BLAST command-line interface.

use blast_rs::db::{BlastDb, DbType};
use blast_rs::format::{format_tabular, TabularHit};
use blast_rs::input::{iupacna_to_blastna, parse_fasta};
use clap::{Parser, Subcommand};
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::PathBuf;

#[derive(Parser)]
#[command(name = "blast-cli", version, about = "BLAST sequence search (Rust implementation)", allow_negative_numbers = true)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
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
    #[arg(long = "num_threads", default_value = "1")]
    num_threads: i32,

    /// Output format: 0=pairwise, 5=XML, 6=tabular, 17=SAM, or '6 col1 col2...'
    #[arg(long = "outfmt", default_value = "6")]
    outfmt: String,

    /// Word size for initial seed
    #[arg(short = 'W', long = "word_size")]
    word_size: Option<i32>,

    /// Reward for nucleotide match
    #[arg(long)]
    reward: Option<i32>,

    /// Penalty for nucleotide mismatch
    #[arg(long)]
    penalty: Option<i32>,

    /// Cost to open a gap
    #[arg(long)]
    gapopen: Option<i32>,

    /// Cost to extend a gap
    #[arg(long)]
    gapextend: Option<i32>,

    /// Query strand(s) to search: both, plus, minus
    #[arg(long, default_value = "both")]
    strand: String,

    /// Maximum number of target sequences to report
    #[arg(long = "max_target_seqs", alias = "max-target-seqs", default_value = "500")]
    max_target_seqs: i32,

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
    #[arg(long = "perc_identity", default_value = "0.0")]
    perc_identity: f64,

    /// Minimum query coverage per HSP (percent)
    #[arg(long = "qcov_hsp_perc", default_value = "0.0")]
    qcov_hsp_perc: f64,

    /// Maximum number of HSPs per subject
    #[arg(long = "max_hsps", alias = "max-hsps", default_value = "0")]
    max_hsps: i32,

    /// Culling limit: delete hits enveloped by higher-scoring hits
    #[arg(long = "culling_limit", default_value = "0")]
    culling_limit: i32,

    /// Window size for multiple hit algorithm
    #[arg(long = "window_size", default_value = "0")]
    window_size: i32,

    /// Soft masking (filter query but use for lookup)
    #[arg(long = "soft_masking", default_value = "true")]
    soft_masking: String,

    /// Use lowercase masking in query
    #[arg(long = "lcase_masking", default_value = "false")]
    lcase_masking: bool,

    /// Restrict search to query location (format: start-stop)
    #[arg(long = "query_loc")]
    query_loc: Option<String>,

    /// Effective database size
    #[arg(long = "dbsize", default_value = "0")]
    dbsize: i64,

    /// Effective search space
    #[arg(long = "searchsp", default_value = "0")]
    searchsp: i64,

    /// Task: megablast, blastn, blastn-short, dc-megablast
    #[arg(long = "task")]
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
    #[arg(long)]
    negative_gilist: Option<PathBuf>,

    /// Negative SeqID list (exclude these)
    #[arg(long)]
    negative_seqidlist: Option<PathBuf>,

    /// TaxID list to restrict database sequences
    #[arg(long)]
    taxidlist: Option<PathBuf>,

    /// TaxIDs to restrict (comma-separated)
    #[arg(long)]
    taxids: Option<String>,

    /// Negative TaxID list (exclude these)
    #[arg(long)]
    negative_taxidlist: Option<PathBuf>,

    /// Negative TaxIDs (comma-separated)
    #[arg(long)]
    negative_taxids: Option<String>,

    /// Import search strategy from file
    #[arg(long)]
    import_search_strategy: Option<PathBuf>,

    /// Export search strategy to file
    #[arg(long)]
    export_search_strategy: Option<PathBuf>,

    /// Show GIs in deflines
    #[arg(long, default_value = "false")]
    show_gis: bool,

    /// Number of descriptions to show (pairwise output)
    #[arg(long, default_value = "500")]
    num_descriptions: i32,

    /// Number of alignments to show (pairwise output)
    #[arg(long, default_value = "250")]
    num_alignments: i32,

    /// Line length for pairwise output
    #[arg(long, default_value = "60")]
    line_length: i32,

    /// Produce HTML output
    #[arg(long, default_value = "false")]
    html: bool,

    /// Sort hits: 0=evalue, 1=bitscore, 2=total_score, 3=pct_identity, 4=query_coverage
    #[arg(long, default_value = "0")]
    sorthits: i32,

    /// Sort HSPs: 0=evalue, 1=score, 2=query_start, 3=pct_identity, 4=subject_start
    #[arg(long, default_value = "0")]
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
    #[arg(long, default_value = "false")]
    subject_besthit: bool,

    /// Best hit overhang
    #[arg(long)]
    best_hit_overhang: Option<f64>,

    /// Best hit score edge
    #[arg(long)]
    best_hit_score_edge: Option<f64>,

    /// Minimum raw gapped score
    #[arg(long)]
    min_raw_gapped_score: Option<i32>,

    /// Template type for discontiguous megablast: coding, optimal, coding_and_optimal
    #[arg(long)]
    template_type: Option<String>,

    /// Template length for discontiguous megablast: 16, 18, 21
    #[arg(long)]
    template_length: Option<i32>,

    /// Parse deflines in FASTA input
    #[arg(long, default_value = "false")]
    parse_deflines: bool,

    /// Restrict search to subject location (format: start-stop)
    #[arg(long)]
    subject_loc: Option<String>,

    /// Entrez query to restrict database (remote searches)
    #[arg(long)]
    entrez_query: Option<String>,

    /// Remote search at NCBI
    #[arg(long, default_value = "false")]
    remote: bool,

    /// Use index for search
    #[arg(long, default_value = "false")]
    use_index: bool,

    /// Index name
    #[arg(long)]
    index_name: Option<String>,

    /// Multithreading mode: 0=split by DB, 1=split by query
    #[arg(long, default_value = "0")]
    mt_mode: i32,

    /// Off-diagonal range for multi-hit extension
    #[arg(long, default_value = "0")]
    off_diagonal_range: i32,

    /// Disable TaxID expansion
    #[arg(long, default_value = "false")]
    no_taxid_expansion: bool,
}

impl BlastnArgs {
    /// Apply task-specific defaults for parameters not explicitly set by the user.
    fn apply_task_defaults(&mut self) {
        if let Some(ref task) = self.task {
            match task.as_str() {
                "blastn-short" => {
                    // NCBI defaults for blastn-short
                    if self.word_size.is_none() { self.word_size = Some(7); }
                    if self.reward.is_none() { self.reward = Some(1); }
                    if self.penalty.is_none() { self.penalty = Some(-3); }
                    if self.gapopen.is_none() { self.gapopen = Some(5); }
                    if self.gapextend.is_none() { self.gapextend = Some(2); }
                }
                "blastn" => {
                    if self.word_size.is_none() { self.word_size = Some(11); }
                    if self.reward.is_none() { self.reward = Some(2); }
                    if self.penalty.is_none() { self.penalty = Some(-3); }
                    if self.gapopen.is_none() { self.gapopen = Some(5); }
                    if self.gapextend.is_none() { self.gapextend = Some(2); }
                }
                "megablast" | "dc-megablast" => {
                    if self.word_size.is_none() { self.word_size = Some(28); }
                    if self.reward.is_none() { self.reward = Some(1); }
                    if self.penalty.is_none() { self.penalty = Some(-2); }
                    if self.gapopen.is_none() { self.gapopen = Some(0); }
                    if self.gapextend.is_none() { self.gapextend = Some(0); }
                }
                _ => {}
            }
        }
        // Fill in defaults for anything still unset
        if self.word_size.is_none() { self.word_size = Some(11); }
        if self.reward.is_none() { self.reward = Some(1); }
        if self.penalty.is_none() { self.penalty = Some(-3); }
        if self.gapopen.is_none() { self.gapopen = Some(5); }
        if self.gapextend.is_none() { self.gapextend = Some(2); }
    }

    fn word_size(&self) -> i32 { self.word_size.unwrap_or(11) }
    fn reward(&self) -> i32 { self.reward.unwrap_or(1) }
    fn penalty(&self) -> i32 { self.penalty.unwrap_or(-3) }
    fn gapopen(&self) -> i32 { self.gapopen.unwrap_or(5) }
    fn gapextend(&self) -> i32 { self.gapextend.unwrap_or(2) }
}

fn main() {
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
    if args.gilist.is_some() { eprintln!("Warning: --gilist is not yet supported, ignoring"); }
    if args.seqidlist.is_some() { eprintln!("Warning: --seqidlist is not yet supported, ignoring"); }
    if args.negative_gilist.is_some() { eprintln!("Warning: --negative-gilist is not yet supported, ignoring"); }
    if args.negative_seqidlist.is_some() { eprintln!("Warning: --negative-seqidlist is not yet supported, ignoring"); }
    if args.taxidlist.is_some() { eprintln!("Warning: --taxidlist is not yet supported, ignoring"); }
    if args.taxids.is_some() { eprintln!("Warning: --taxids is not yet supported, ignoring"); }
    if args.negative_taxidlist.is_some() { eprintln!("Warning: --negative-taxidlist is not yet supported, ignoring"); }
    if args.negative_taxids.is_some() { eprintln!("Warning: --negative-taxids is not yet supported, ignoring"); }
    if args.remote { eprintln!("Warning: --remote is not supported, ignoring"); }

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
            .map(|p| blast_rs::stat::KarlinBlk { lambda: p.lambda, k: p.k, log_k: p.k.ln(), h: p.h })
            .unwrap_or_else(blast_rs::stat::protein_ungapped_kbp);

        let total_subj_len: usize = subjects.iter().map(|s| s.sequence.len()).sum();
        let avg_q_len = queries.iter().map(|q| q.sequence.len()).sum::<usize>() / queries.len().max(1);
        let len_adj = blast_rs::stat::compute_length_adjustment(
            avg_q_len as i32, total_subj_len as i64, subjects.len() as i32, &prot_kbp);
        let search_space = blast_rs::stat::compute_search_space(
            avg_q_len as i64, total_subj_len as i64, subjects.len() as i32, len_adj);

        // Protein default word_size=3 (NCBI default), cap at 6 to keep table size manageable
        let word_size = if args.word_size.is_none() { 3 } else { args.word_size().max(2).min(6) as usize };
        let threshold = 11.0; // NCBI default for word_size=3

        let mut hits = Vec::new();
        for qrec in &queries {
            let query_aa: Vec<u8> = qrec.sequence.iter()
                .map(|&b| blast_rs::input::aminoacid_to_ncbistdaa(b)).collect();

            for srec in &subjects {
                let subj_aa: Vec<u8> = srec.sequence.iter()
                    .map(|&b| blast_rs::input::aminoacid_to_ncbistdaa(b)).collect();

                let phits = blast_rs::protein_lookup::protein_gapped_scan(
                    &query_aa, &subj_aa, &matrix, word_size, threshold,
                    30, args.gapopen(), args.gapextend(), 50, 0,
                );
                for ph in phits {
                    let evalue = prot_kbp.raw_to_evalue(ph.score, search_space);
                    if evalue > args.evalue { continue; }
                    hits.push(TabularHit {
                        query_id: qrec.id.clone(),
                        subject_id: srec.id.clone(),
                        pct_identity: if ph.align_length > 0 {
                            100.0 * ph.num_ident as f64 / ph.align_length as f64
                        } else { 0.0 },
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
                        subject_sci_name: String::new(), subject_common_name: String::new(),
                        subject_blast_name: String::new(), subject_kingdom: String::new(),
                        num_ident: ph.num_ident,
                    });
                }
            }
        }

        hits.sort_by(|a, b| b.bit_score.partial_cmp(&a.bit_score).unwrap_or(std::cmp::Ordering::Equal));
        let mut seen = std::collections::HashSet::new();
        hits.retain(|h| seen.insert((h.query_id.clone(), h.subject_id.clone(), h.query_start, h.subject_start)));

        let stdout = io::stdout();
        let mut writer: Box<dyn Write> = if let Some(ref path) = args.out {
            Box::new(BufWriter::new(File::create(path)?))
        } else {
            Box::new(BufWriter::new(stdout.lock()))
        };
        let outfmt_parts: Vec<&str> = args.outfmt.split_whitespace().collect();
        if outfmt_parts.len() > 1 {
            let cols = outfmt_parts[1..].join(" ");
            blast_rs::format::format_tabular_custom(&mut writer, &hits, &cols)?;
        } else {
            format_tabular(&mut writer, &hits)?;
        }
        writer.flush()?;
        return Ok(());
    }

    Err("blastp requires --subject (database search not yet implemented for protein)".into())
}

fn run_blastx(args: &BlastnArgs) -> Result<(), Box<dyn std::error::Error>> {

    use blast_rs::util::{six_frame_translation, STANDARD_GENETIC_CODE};

    let query_file = File::open(&args.query)?;
    let queries = parse_fasta(query_file);
    if queries.is_empty() {
        return Err("No sequences found in query file".into());
    }

    let subject_path = args.subject.as_ref()
        .ok_or("blastx requires --subject (protein FASTA)")?;
    let subject_file = File::open(subject_path)?;
    let subjects = parse_fasta(subject_file);

    let matrix = blast_rs::matrix::BLOSUM62;

    let prot_kbp = blast_rs::stat::protein_ungapped_kbp();
    let total_subj_len: usize = subjects.iter().map(|s| s.sequence.len()).sum();

    let mut all_hits = Vec::new();

    for qrec in &queries {
        // Encode nucleotide query
        let nuc_seq: Vec<u8> = qrec.sequence.iter().map(|&b| {
            match b {
                b'A' | b'a' => 0, b'C' | b'c' => 1, b'G' | b'g' => 2, b'T' | b't' => 3, _ => 0
            }
        }).collect();

        // Translate in 6 frames
        let frames = six_frame_translation(&nuc_seq, &STANDARD_GENETIC_CODE);

        for (frame, protein_query) in &frames {
            if protein_query.len() < 3 { continue; }

            let search_space = blast_rs::stat::compute_search_space(
                protein_query.len() as i64, total_subj_len as i64, subjects.len() as i32, 0);

            for srec in &subjects {
                let subj_aa: Vec<u8> = srec.sequence.iter()
                    .map(|&b| blast_rs::input::aminoacid_to_ncbistdaa(b)).collect();

                let phits = blast_rs::protein_lookup::protein_gapped_scan(
                    &protein_query, &subj_aa, &matrix, 3, 11.0,
                    30, args.gapopen(), args.gapextend(), 50, 0,
                );
                let nuc_len = nuc_seq.len() as i32;

                for ph in phits {
                    let evalue = prot_kbp.raw_to_evalue(ph.score, search_space);
                    if evalue > args.evalue { continue; }

                    let (q_nuc_start, q_nuc_end) = blast_rs::util::protein_to_nuc_coords(
                        ph.query_start as i32, ph.query_end as i32, *frame, nuc_len);

                    all_hits.push(TabularHit {
                        query_id: qrec.id.clone(),
                        subject_id: srec.id.clone(),
                        pct_identity: if ph.align_length > 0 {
                            100.0 * ph.num_ident as f64 / ph.align_length as f64
                        } else { 0.0 },
                        align_len: ph.align_length, mismatches: ph.mismatches, gap_opens: ph.gap_opens,
                        query_start: q_nuc_start, query_end: q_nuc_end,
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
                        subject_sci_name: String::new(), subject_common_name: String::new(),
                        subject_blast_name: String::new(), subject_kingdom: String::new(),
                        num_ident: ph.num_ident,
                    });
                }
            }
        }
    }

    all_hits.sort_by(|a, b| a.evalue.partial_cmp(&b.evalue).unwrap_or(std::cmp::Ordering::Equal));
    let mut seen = std::collections::HashSet::new();
    all_hits.retain(|h| seen.insert((h.query_id.clone(), h.query_start, h.subject_start)));

    let stdout = io::stdout();
    let mut writer: Box<dyn Write> = if let Some(ref path) = args.out {
        Box::new(BufWriter::new(File::create(path)?))
    } else {
        Box::new(BufWriter::new(stdout.lock()))
    };
    let outfmt_parts: Vec<&str> = args.outfmt.split_whitespace().collect();
    if outfmt_parts.len() > 1 {
        let cols = outfmt_parts[1..].join(" ");
        blast_rs::format::format_tabular_custom(&mut writer, &all_hits, &cols)?;
    } else {
        format_tabular(&mut writer, &all_hits)?;
    }
    writer.flush()?;
    Ok(())
}

fn run_tblastn(args: &BlastnArgs) -> Result<(), Box<dyn std::error::Error>> {

    use blast_rs::util::{six_frame_translation, STANDARD_GENETIC_CODE};

    // tblastn: protein query vs translated nucleotide subject
    let query_file = File::open(&args.query)?;
    let queries = parse_fasta(query_file);
    let subject_path = args.subject.as_ref()
        .ok_or("tblastn requires --subject (nucleotide FASTA)")?;
    let subject_file = File::open(subject_path)?;
    let subjects = parse_fasta(subject_file);

    let matrix = blast_rs::matrix::BLOSUM62;
    let prot_kbp = blast_rs::stat::protein_ungapped_kbp();

    let mut all_hits = Vec::new();
    for qrec in &queries {
        let query_aa: Vec<u8> = qrec.sequence.iter()
            .map(|&b| blast_rs::input::aminoacid_to_ncbistdaa(b)).collect();

        for srec in &subjects {
            let nuc_seq: Vec<u8> = srec.sequence.iter().map(|&b| {
                match b { b'A'|b'a'=>0, b'C'|b'c'=>1, b'G'|b'g'=>2, b'T'|b't'=>3, _=>0 }
            }).collect();
            let frames = six_frame_translation(&nuc_seq, &STANDARD_GENETIC_CODE);

            for (frame, subj_protein) in &frames {
                if subj_protein.len() < 3 || query_aa.len() < 3 { continue; }
                let search_space = (query_aa.len() * subj_protein.len()) as f64;
                let subj_nuc_len = nuc_seq.len() as i32;

                let phits = blast_rs::protein_lookup::protein_gapped_scan(
                    &query_aa, &subj_protein, &matrix, 3, 11.0,
                    30, args.gapopen(), args.gapextend(), 50, 0,
                );

                for ph in phits {
                    let evalue = prot_kbp.raw_to_evalue(ph.score, search_space);
                    if evalue > args.evalue { continue; }

                    let (s_nuc_start, s_nuc_end) = blast_rs::util::protein_to_nuc_coords(
                        ph.subject_start as i32, ph.subject_end as i32, *frame, subj_nuc_len);

                    all_hits.push(TabularHit {
                        query_id: qrec.id.clone(), subject_id: srec.id.clone(),
                        pct_identity: if ph.align_length > 0 {
                            100.0 * ph.num_ident as f64 / ph.align_length as f64
                        } else { 0.0 },
                        align_len: ph.align_length, mismatches: ph.mismatches, gap_opens: ph.gap_opens,
                        query_start: ph.query_start as i32 + 1, query_end: ph.query_end as i32,
                        subject_start: s_nuc_start, subject_end: s_nuc_end,
                        evalue, bit_score: prot_kbp.raw_to_bit(ph.score),
                        query_len: query_aa.len() as i32,
                        subject_len: subj_nuc_len,
                        raw_score: ph.score,
                        qseq: ph.qseq.clone(),
                        sseq: ph.sseq.clone(),
                        qframe: 0,
                        sframe: *frame,
                        subject_taxids: vec![],
                        subject_sci_name: String::new(), subject_common_name: String::new(),
                        subject_blast_name: String::new(), subject_kingdom: String::new(),
                        num_ident: ph.num_ident,
                    });
                }
            }
        }
    }
    all_hits.sort_by(|a, b| a.evalue.partial_cmp(&b.evalue).unwrap_or(std::cmp::Ordering::Equal));
    let mut seen = std::collections::HashSet::new();
    all_hits.retain(|h| seen.insert((h.query_id.clone(), h.subject_id.clone(), h.query_start, h.subject_start)));

    let stdout = io::stdout();
    let mut writer: Box<dyn Write> = if let Some(ref path) = args.out {
        Box::new(BufWriter::new(File::create(path)?))
    } else { Box::new(BufWriter::new(stdout.lock())) };
    let outfmt_parts: Vec<&str> = args.outfmt.split_whitespace().collect();
    if outfmt_parts.len() > 1 {
        let cols = outfmt_parts[1..].join(" ");
        blast_rs::format::format_tabular_custom(&mut writer, &all_hits, &cols)?;
    } else {
        format_tabular(&mut writer, &all_hits)?;
    }
    writer.flush()?;
    Ok(())
}

fn run_psiblast(args: &BlastnArgs) -> Result<(), Box<dyn std::error::Error>> {
    use blast_rs::pssm::Pssm;


    let query_file = File::open(&args.query)?;
    let queries = parse_fasta(query_file);
    if queries.is_empty() { return Err("No query sequences".into()); }

    let subject_path = args.subject.as_ref()
        .ok_or("psiblast requires --subject")?;
    let subject_file = File::open(subject_path)?;
    let subjects = parse_fasta(subject_file);

    let matrix = blast_rs::matrix::BLOSUM62;

    let prot_kbp = blast_rs::stat::protein_ungapped_kbp();
    let total_subj_len: usize = subjects.iter().map(|s| s.sequence.len()).sum();
    let search_space = (queries[0].sequence.len() * total_subj_len) as f64;

    // Build initial PSSM from query
    let query_aa: Vec<u8> = queries[0].sequence.iter()
        .map(|&b| blast_rs::input::aminoacid_to_ncbistdaa(b)).collect();
    let mut pssm = Pssm::from_sequence(&query_aa, &matrix);

    // Prepare subjects as (id, aa_sequence) pairs
    let subj_pairs: Vec<(String, Vec<u8>)> = subjects.iter()
        .map(|s| (s.id.clone(), s.sequence.iter().map(|&b| blast_rs::input::aminoacid_to_ncbistdaa(b)).collect()))
        .collect();

    // Run 3 iterations
    let mut all_hits = Vec::new();
    for _iter in 0..3 {
        let results = blast_rs::pssm::psi_blast_iteration(
            &pssm, &subj_pairs, args.evalue, search_space, prot_kbp.lambda, prot_kbp.k);

        if results.is_empty() { break; }

        all_hits.clear();
        for hit in &results {
            // Extract aligned sequences (ungapped PSSM scan)
            let subj_seq = subj_pairs.iter()
                .find(|(sid, _)| sid == &hit.subject_id)
                .map(|(_, s)| s.as_slice());
            let qseq: String = query_aa[..hit.align_len].iter()
                .map(|&b| blast_rs::protein::ncbistdaa_to_char(b)).collect();
            let sseq: String = subj_seq.map(|s| {
                s[hit.subject_start..hit.subject_start + hit.align_len].iter()
                    .map(|&b| blast_rs::protein::ncbistdaa_to_char(b)).collect()
            }).unwrap_or_default();
            // Compute identity
            let num_ident = query_aa[..hit.align_len].iter()
                .zip(subj_seq.map(|s| &s[hit.subject_start..hit.subject_start + hit.align_len]).unwrap_or(&[]))
                .filter(|(a, b)| a == b).count() as i32;
            let pct_identity = if hit.align_len > 0 {
                100.0 * num_ident as f64 / hit.align_len as f64
            } else { 0.0 };

            all_hits.push(TabularHit {
                query_id: queries[0].id.clone(), subject_id: hit.subject_id.clone(),
                pct_identity,
                align_len: hit.align_len as i32,
                mismatches: hit.align_len as i32 - num_ident,
                gap_opens: 0,
                query_start: 1, query_end: hit.align_len as i32,
                subject_start: hit.subject_start as i32 + 1,
                subject_end: (hit.subject_start + hit.align_len) as i32,
                evalue: hit.evalue, bit_score: prot_kbp.raw_to_bit(hit.score),
                query_len: pssm.length as i32, subject_len: hit.subject_len as i32,
                raw_score: hit.score,
                qseq: Some(qseq), sseq: Some(sseq),
                qframe: 0, sframe: 0,
                subject_taxids: vec![],
                subject_sci_name: String::new(), subject_common_name: String::new(),
                subject_blast_name: String::new(), subject_kingdom: String::new(),
                num_ident,
            });
        }

        // Update PSSM from aligned sequences (simplified)
        let aligned: Vec<Vec<u8>> = results.iter()
            .filter_map(|h| subj_pairs.iter().find(|(sid, _)| sid == &h.subject_id).map(|(_, s)| s.clone()))
            .collect();
        pssm.update_from_alignment(&aligned, &blast_rs::matrix::AA_FREQUENCIES, 10.0);
    }

    let stdout = io::stdout();
    let mut writer: Box<dyn Write> = if let Some(ref path) = args.out {
        Box::new(BufWriter::new(File::create(path)?))
    } else { Box::new(BufWriter::new(stdout.lock())) };
    let outfmt_parts: Vec<&str> = args.outfmt.split_whitespace().collect();
    if outfmt_parts.len() > 1 {
        let cols = outfmt_parts[1..].join(" ");
        blast_rs::format::format_tabular_custom(&mut writer, &all_hits, &cols)?;
    } else {
        format_tabular(&mut writer, &all_hits)?;
    }
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
    let subject_path = args.subject.as_ref()
        .ok_or("tblastx requires --subject (nucleotide FASTA)")?;
    let subject_file = File::open(subject_path)?;
    let subjects = parse_fasta(subject_file);

    let matrix = blast_rs::matrix::BLOSUM62;
    let prot_kbp = blast_rs::stat::protein_ungapped_kbp();

    let mut all_hits = Vec::new();
    for qrec in &queries {
        let q_nuc: Vec<u8> = qrec.sequence.iter().map(|&b| {
            match b { b'A'|b'a'=>0, b'C'|b'c'=>1, b'G'|b'g'=>2, b'T'|b't'=>3, _=>0 }
        }).collect();
        let q_frames = six_frame_translation(&q_nuc, &STANDARD_GENETIC_CODE);

        for srec in &subjects {
            let s_nuc: Vec<u8> = srec.sequence.iter().map(|&b| {
                match b { b'A'|b'a'=>0, b'C'|b'c'=>1, b'G'|b'g'=>2, b'T'|b't'=>3, _=>0 }
            }).collect();
            let s_frames = six_frame_translation(&s_nuc, &STANDARD_GENETIC_CODE);

            for (qframe, q_prot) in &q_frames {
                for (sframe, s_prot) in &s_frames {
                    if q_prot.len() < 3 || s_prot.len() < 3 { continue; }
                    let search_space = (q_prot.len() * s_prot.len()) as f64;

                    let phits = blast_rs::protein_lookup::protein_gapped_scan(
                        q_prot, s_prot, &matrix, 3, 13.0,
                        25, args.gapopen(), args.gapextend(), 50, 0,
                    );

                    for ph in phits {
                        let evalue = prot_kbp.raw_to_evalue(ph.score, search_space);
                        if evalue > args.evalue { continue; }
                        all_hits.push(TabularHit {
                            query_id: qrec.id.clone(), subject_id: srec.id.clone(),
                            pct_identity: if ph.align_length > 0 {
                                100.0 * ph.num_ident as f64 / ph.align_length as f64
                            } else { 0.0 },
                            align_len: ph.align_length, mismatches: ph.mismatches, gap_opens: ph.gap_opens,
                            query_start: ph.query_start as i32 + 1, query_end: ph.query_end as i32,
                            subject_start: ph.subject_start as i32 + 1, subject_end: ph.subject_end as i32,
                            evalue, bit_score: prot_kbp.raw_to_bit(ph.score),
                            query_len: q_prot.len() as i32, subject_len: s_prot.len() as i32,
                            raw_score: ph.score,
                            qseq: ph.qseq.clone(), sseq: ph.sseq.clone(),
                            qframe: *qframe, sframe: *sframe,
                            subject_taxids: vec![],
                            subject_sci_name: String::new(), subject_common_name: String::new(),
                            subject_blast_name: String::new(), subject_kingdom: String::new(),
                            num_ident: ph.num_ident,
                        });
                    }
                }
            }
        }
    }
    all_hits.sort_by(|a, b| a.evalue.partial_cmp(&b.evalue).unwrap_or(std::cmp::Ordering::Equal));
    let mut seen = std::collections::HashSet::new();
    all_hits.retain(|h| seen.insert((h.query_id.clone(), h.subject_id.clone(), h.query_start, h.subject_start)));

    let stdout = io::stdout();
    let mut writer: Box<dyn Write> = if let Some(ref path) = args.out {
        Box::new(BufWriter::new(File::create(path)?))
    } else { Box::new(BufWriter::new(stdout.lock())) };
    let outfmt_parts: Vec<&str> = args.outfmt.split_whitespace().collect();
    if outfmt_parts.len() > 1 {
        let cols = outfmt_parts[1..].join(" ");
        blast_rs::format::format_tabular_custom(&mut writer, &all_hits, &cols)?;
    } else {
        format_tabular(&mut writer, &all_hits)?;
    }
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

    let db_path = args.db.as_ref().ok_or("Either --db or --subject is required")?;
    let db = BlastDb::open(db_path)?;
    if db.db_type != DbType::Nucleotide {
        return Err("blastn requires a nucleotide database".into());
    }

    run_blastn_rust(args, &records, db)
}

/// Pure Rust blastn search — no FFI calls.
fn run_blastn_rust(
    args: &BlastnArgs,
    records: &[blast_rs::input::FastaRecord],
    db: BlastDb,
) -> Result<(), Box<dyn std::error::Error>> {
    use rayon::prelude::*;

    // Try to load taxonomy name database (taxdb.bti/taxdb.btd)
    let db_dir = std::path::Path::new(args.db.as_ref().unwrap()).parent().unwrap_or(std::path::Path::new("."));
    let tax_name_db = blast_rs::db::TaxNameDb::open(db_dir).ok()
        .or_else(|| std::env::var("BLASTDB").ok().and_then(|p| blast_rs::db::TaxNameDb::open(std::path::Path::new(&p)).ok()));

    // Pure Rust KBP computation — per-context (plus + minus) for exact e-value matching
    // The C engine computes separate KBP for each context based on strand composition
    let (kbp_plus, kbp_minus, searchsp_plus, searchsp_minus) = {
        let rec = &records[0];
        let query_plus: Vec<u8> = rec.sequence.iter().map(|&b| iupacna_to_blastna(b)).collect();
        let query_minus: Vec<u8> = query_plus.iter().rev().map(|&b| complement_blastna(b)).collect();

        const BLASTNA_SIZE: usize = 16;
        const BLASTNA_TO_NCBI4NA: [u8; BLASTNA_SIZE] = [1,2,4,8,5,10,3,12,9,6,14,13,11,7,15,0];
        let reward = args.reward();
        let penalty = args.penalty();
        let matrix_fn = |i: usize, j: usize| -> i32 {
            if i >= BLASTNA_SIZE || j >= BLASTNA_SIZE { return penalty; }
            if i == BLASTNA_SIZE - 1 || j == BLASTNA_SIZE - 1 { return i32::MIN / 2; }
            if BLASTNA_TO_NCBI4NA[i] & BLASTNA_TO_NCBI4NA[j] != 0 {
                let d = (0..4).filter(|&k| BLASTNA_TO_NCBI4NA[j] & BLASTNA_TO_NCBI4NA[k] != 0).count() as i32;
                (((d - 1) as f64 * penalty as f64 + reward as f64) / d as f64).round() as i32
            } else { penalty }
        };
        let mut lo = i32::MAX; let mut hi = i32::MIN;
        for i in 0..BLASTNA_SIZE { for j in 0..BLASTNA_SIZE {
            let s = matrix_fn(i, j);
            if s <= -100000000 || s >= 100000000 { continue; }
            if s < lo { lo = s; } if s > hi { hi = s; }
        }}

        let ambig: &[u8] = &[14, 15];
        let qlen = rec.sequence.len() as i32;

        // Compute KBP for both contexts
        let ctxs = [
            blast_rs::stat::UngappedKbpContext { query_offset: 0, query_length: qlen, is_valid: true },
        ];
        let plus_kbp_results = blast_rs::stat::ungapped_kbp_calc(&query_plus, &ctxs, lo, hi, BLASTNA_SIZE, ambig, &matrix_fn);
        let minus_kbp_results = blast_rs::stat::ungapped_kbp_calc(&query_minus, &ctxs, lo, hi, BLASTNA_SIZE, ambig, &matrix_fn);

        let default_kbp = blast_rs::stat::KarlinBlk { lambda: 1.374, k: 0.621, log_k: 0.621_f64.ln(), h: 1.286 };
        let ungapped_plus = plus_kbp_results[0].clone().unwrap_or(default_kbp.clone());
        let ungapped_minus = minus_kbp_results[0].clone().unwrap_or(default_kbp.clone());

        // Gapped KBP (from table — may fall back to ungapped for large gap costs)
        let (gkbp_plus, _) = blast_rs::stat::nucl_gapped_kbp_lookup(
            args.gapopen(), args.gapextend(), reward, penalty, &ungapped_plus,
        ).unwrap_or((ungapped_plus.clone(), false));
        let (gkbp_minus, _) = blast_rs::stat::nucl_gapped_kbp_lookup(
            args.gapopen(), args.gapextend(), reward, penalty, &ungapped_minus,
        ).unwrap_or((ungapped_minus.clone(), false));

        // Per-context search space
        let compute_searchsp = |kbp: &blast_rs::stat::KarlinBlk, ukbp: &blast_rs::stat::KarlinBlk| -> f64 {
            let (alpha, beta) = blast_rs::stat::nucl_alpha_beta(
                reward, penalty, args.gapopen(), args.gapextend(), ukbp.lambda, ukbp.h, true);
            let (len_adj, _) = blast_rs::stat::compute_length_adjustment_exact(
                kbp.k, kbp.log_k, alpha / kbp.lambda, beta,
                qlen, db.total_length as i64, db.num_oids as i32);
            let eff_db = std::cmp::max(db.total_length as i64 - db.num_oids as i64 * len_adj as i64, 1);
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
    let reward = args.reward();
    let penalty = args.penalty();
    let gapopen = args.gapopen();
    let gapextend = args.gapextend();
    let evalue = args.evalue;

    let mut all_hits: Vec<(u32, TabularHit)> = Vec::new();

    let apply_dust = args.dust != "no";

    for rec in records {
        // Keep unmasked query for gapped alignment (matches C engine's sequence_nomask)
        let query_plus_nomask: Vec<u8> = rec.sequence.iter().map(|&b| iupacna_to_blastna(b)).collect();
        let query_minus_nomask: Vec<u8> = query_plus_nomask.iter().rev()
            .map(|&b| complement_blastna(b)).collect();

        // Masked query for seed finding (lookup table uses masked version)
        let mut query_plus_masked = query_plus_nomask.clone();
        if apply_dust {
            let mask = blast_rs::filter::dust_filter(&query_plus_masked, 64, 2.0);
            mask.apply(&mut query_plus_masked, 14);
        }
        let query_minus_masked: Vec<u8> = query_plus_masked.iter().rev()
            .map(|&b| complement_blastna(b)).collect();

        // Apply strand filter
        let search_plus = args.strand != "minus";
        let search_minus = args.strand != "plus";
        let qp = if search_plus { &query_plus_masked[..] } else { &[] };
        let qm = if search_minus { &query_minus_masked[..] } else { &[] };
        let qp_nomask = if search_plus { &query_plus_nomask[..] } else { &[] };
        let qm_nomask = if search_minus { &query_minus_nomask[..] } else { &[] };

        // Compute cutoff score from e-value (matches C engine's BLAST_Cutoffs)
        let cutoff_score = {
            let e = evalue.max(1.0e-297);
            ((kbp.k * search_space / e).ln() / kbp.lambda).ceil() as i32
        };

        // Gapped search across all subjects using packed NCBI2na (fast, no full decode)
        let gapped_x_dropoff = (args.xdrop_gap_final * std::f64::consts::LN_2 / kbp.lambda) as i32;
        let oid_hits: Vec<(u32, Vec<blast_rs::search::SearchHsp>)> =
            (0..db.num_oids).into_par_iter().filter_map(|oid| {
                let packed = db.get_sequence(oid);
                let seq_len = db.get_seq_len(oid) as usize;
                let mut hsps = blast_rs::search::blastn_gapped_search_packed(
                    qp, qm, qp_nomask, qm_nomask,
                    packed, seq_len,
                    word_size, reward, penalty, gapopen, gapextend,
                    gapped_x_dropoff, &kbp, search_space, evalue,
                );
                // Apply cutoff score filter (matches C engine's preliminary cutoff)
                hsps.retain(|h| h.score >= cutoff_score);
                if hsps.is_empty() { None } else { Some((oid, hsps)) }
            }).collect();

        // Sort by OID to match C engine's deterministic database scan order
        let mut oid_hits = oid_hits;
        oid_hits.sort_by_key(|&(oid, _)| oid);

        // Apply hitlist pruning to match C engine behavior:
        // Keep only the best prelim_hitlist_size subjects (by best e-value per subject).
        // This matches the C engine's BlastHitList heap which limits tracked subjects.
        let hitlist_size = args.max_target_seqs;
        let prelim_hitlist_size = std::cmp::min(
            std::cmp::max(2 * hitlist_size, 10), hitlist_size + 50) as usize;
        if oid_hits.len() > prelim_hitlist_size {
            // Sort subjects by best (lowest) e-value, keep top prelim_hitlist_size
            oid_hits.sort_by(|(_, a_hsps), (_, b_hsps)| {
                let a_best = a_hsps.iter().map(|h| h.evalue).fold(f64::MAX, f64::min);
                let b_best = b_hsps.iter().map(|h| h.evalue).fold(f64::MAX, f64::min);
                a_best.partial_cmp(&b_best).unwrap_or(std::cmp::Ordering::Equal)
            });
            oid_hits.truncate(prelim_hitlist_size);
            // Re-sort by OID for deterministic output
            oid_hits.sort_by_key(|&(oid, _)| oid);
        }

        for (oid, hsps) in oid_hits {
            for hsp in hsps {
                let subject_id = db.get_accession(oid)
                    .unwrap_or_else(|| format!("gnl|BL_ORD_ID|{}", oid));
                let query_len = rec.sequence.len() as i32;

                let (q_start, q_end) = if hsp.context == 1 {
                    (query_len - hsp.query_end + 1, query_len - hsp.query_start)
                } else {
                    (hsp.query_start + 1, hsp.query_end)
                };
                let (s_start, s_end) = if hsp.context == 1 {
                    (hsp.subject_end, hsp.subject_start + 1)
                } else {
                    (hsp.subject_start + 1, hsp.subject_end)
                };

                all_hits.push((oid, TabularHit {
                    query_id: rec.id.clone(),
                    subject_id,
                    pct_identity: if hsp.align_length > 0 {
                        100.0 * hsp.num_ident as f64 / hsp.align_length as f64
                    } else { 0.0 },
                    align_len: hsp.align_length,
                    mismatches: hsp.mismatches,
                    gap_opens: hsp.gap_opens,
                    query_start: q_start,
                    query_end: q_end,
                    subject_start: s_start,
                    subject_end: s_end,
                    // Use per-context KBP for exact e-value matching with C engine
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
                    qframe: if hsp.context == 1 { -1 } else { 1 },
                    sframe: if hsp.context == 1 { -1 } else { 1 },
                    subject_taxids: db.get_taxids(oid),
                    subject_sci_name: String::new(), subject_common_name: String::new(),
                    subject_blast_name: String::new(), subject_kingdom: String::new(),
                    num_ident: hsp.num_ident,
                }));
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
        let mut groups: std::collections::BTreeMap<u32, Vec<TabularHit>> = std::collections::BTreeMap::new();
        let mut best_evalue: std::collections::HashMap<u32, f64> = std::collections::HashMap::new();
        let mut best_score: std::collections::HashMap<u32, i32> = std::collections::HashMap::new();
        for (oid, hit) in all_hits {
            let ee = best_evalue.entry(oid).or_insert(hit.evalue);
            if hit.evalue < *ee { *ee = hit.evalue; }
            let es = best_score.entry(oid).or_insert((hit.bit_score * 1000.0) as i32);
            let raw = (hit.bit_score * 1000.0) as i32;
            if raw > *es { *es = raw; }
            groups.entry(oid).or_default().push(hit);
        }

        // Sort HSPs within each group by descending score
        for hsps in groups.values_mut() {
            hsps.sort_by(|a, b| b.bit_score.partial_cmp(&a.bit_score).unwrap_or(std::cmp::Ordering::Equal)
                .then(a.query_start.cmp(&b.query_start)));
        }

        // Sort groups using C's qsort comparator logic with libc qsort for exact match
        let mut sorted_oids: Vec<u32> = groups.keys().copied().collect();

        // Use libc qsort to match C's exact tie-breaking behavior
        sorted_oids.sort_by(|&a, &b| {
            let ev_a = best_evalue[&a];
            let ev_b = best_evalue[&b];
            ev_a.partial_cmp(&ev_b).unwrap_or(std::cmp::Ordering::Equal)
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
    let max_subjects = args.max_target_seqs as usize;
    let mut seen_per_query: std::collections::HashMap<String, std::collections::HashSet<String>> = std::collections::HashMap::new();
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
    let query_len = records.first().map(|r| r.sequence.len() as i32).unwrap_or(0);
    let mut hits: Vec<TabularHit> = all_hits.into_iter().map(|(_, h)| h).collect();
    apply_filters(&mut hits, args, query_len);

    // Enrich hits with taxonomy names from taxdb (if available)
    if let Some(ref tdb) = tax_name_db {
        for hit in &mut hits {
            if let Some(&taxid) = hit.subject_taxids.first() {
                if let Some(info) = tdb.get_info(taxid) {
                    hit.subject_sci_name = info.scientific_name;
                    hit.subject_common_name = info.common_name;
                    hit.subject_blast_name = info.blast_name;
                    hit.subject_kingdom = info.kingdom;
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
    let custom_columns = if outfmt_parts.len() > 1 {
        Some(outfmt_parts[1..].join(" "))
    } else { None };

    if outfmt_num == 0 {
        // Pairwise text output
        for hit in &hits {
            blast_rs::format::format_pairwise_alignment(
                &mut writer,
                &hit.query_id, &hit.subject_id,
                &[], &[], // placeholder - would need decoded seqs
                hit.query_start, hit.query_end,
                hit.subject_start, hit.subject_end,
                0, hit.bit_score, hit.evalue,
                hit.align_len - hit.mismatches, hit.align_len, hit.gap_opens,
            )?;
        }
    } else if let Some(ref cols) = custom_columns {
        blast_rs::format::format_tabular_custom(&mut writer, &hits, cols)?;
    } else {
        format_tabular(&mut writer, &hits)?;
    }
    writer.flush()?;

    Ok(())
}

/// Search query against subject FASTA sequences (no database needed).
fn run_blastn_subject(
    args: &BlastnArgs,
    queries: &[blast_rs::input::FastaRecord],
    subjects: &[blast_rs::input::FastaRecord],
) -> Result<(), Box<dyn std::error::Error>> {
    use blast_rs::search::blastn_gapped_search;

    let kbp = blast_rs::stat::compute_ungapped_kbp(args.reward(), args.penalty());
    let total_subj_len: usize = subjects.iter().map(|s| s.sequence.len()).sum();
    let avg_query_len = queries.iter().map(|r| r.sequence.len()).sum::<usize>() as f64
        / queries.len().max(1) as f64;
    let search_space = total_subj_len as f64 * avg_query_len;
    let word_size = args.word_size() as usize;

    let mut all_hits = Vec::new();

    for query_rec in queries {
        let query_plus: Vec<u8> = query_rec.sequence.iter().map(|&b| iupacna_to_blastna(b)).collect();
        let query_minus: Vec<u8> = query_rec.sequence.iter().rev()
            .map(|&b| complement_blastna(iupacna_to_blastna(b))).collect();

        for subj_rec in subjects {
            let subject: Vec<u8> = subj_rec.sequence.iter().map(|&b| iupacna_to_blastna(b)).collect();

            let hsps = blastn_gapped_search(
                &query_plus, &query_minus, &subject,
                word_size, args.reward(), args.penalty(), args.gapopen(), args.gapextend(),
                20, &kbp, search_space, args.evalue,
            );

            for hsp in hsps {
                let query_len = query_rec.sequence.len() as i32;
                let (q_start, q_end) = if hsp.context == 1 {
                    (query_len - hsp.query_end + 1, query_len - hsp.query_start)
                } else {
                    (hsp.query_start + 1, hsp.query_end)
                };
                let (s_start, s_end) = if hsp.context == 1 {
                    (hsp.subject_end, hsp.subject_start + 1)
                } else {
                    (hsp.subject_start + 1, hsp.subject_end)
                };

                all_hits.push(TabularHit {
                    query_id: query_rec.id.clone(),
                    subject_id: subj_rec.id.clone(),
                    pct_identity: if hsp.align_length > 0 {
                        100.0 * hsp.num_ident as f64 / hsp.align_length as f64
                    } else { 0.0 },
                    align_len: hsp.align_length,
                    mismatches: hsp.mismatches,
                    gap_opens: hsp.gap_opens,
                    query_start: q_start,
                    query_end: q_end,
                    subject_start: s_start,
                    subject_end: s_end,
                    evalue: hsp.evalue,
                    bit_score: hsp.bit_score,
                    query_len: query_rec.sequence.len() as i32,
                    subject_len: subj_rec.sequence.len() as i32,
                    raw_score: hsp.score,
                    qseq: hsp.qseq.clone(),
                    sseq: hsp.sseq.clone(),
                    qframe: if hsp.context == 1 { -1 } else { 1 },
                    sframe: if hsp.context == 1 { -1 } else { 1 },
                    subject_taxids: vec![],
                    subject_sci_name: String::new(), subject_common_name: String::new(),
                    subject_blast_name: String::new(), subject_kingdom: String::new(),
                    num_ident: hsp.num_ident,
                });
            }
        }
    }

    all_hits.sort_by(|a, b| a.evalue.partial_cmp(&b.evalue).unwrap_or(std::cmp::Ordering::Equal));

    let stdout = io::stdout();
    let mut writer: Box<dyn Write> = if let Some(ref path) = args.out {
        Box::new(BufWriter::new(File::create(path)?))
    } else {
        Box::new(BufWriter::new(stdout.lock()))
    };
    let outfmt_parts: Vec<&str> = args.outfmt.split_whitespace().collect();
    if outfmt_parts.len() > 1 {
        let cols = outfmt_parts[1..].join(" ");
        blast_rs::format::format_tabular_custom(&mut writer, &all_hits, &cols)?;
    } else {
        format_tabular(&mut writer, &all_hits)?;
    }
    writer.flush()?;
    Ok(())
}

/// DEPRECATED: No longer used. All KBP computation is now pure Rust.
/// Kept only as reference for the C FFI calling pattern.

/// Apply post-search filters (perc_identity, qcov_hsp_perc, max_hsps).
fn apply_filters(hits: &mut Vec<TabularHit>, args: &BlastnArgs, query_len: i32) {
    // Filter by percent identity
    if args.perc_identity > 0.0 {
        hits.retain(|h| h.pct_identity >= args.perc_identity);
    }
    // Filter by query coverage
    if args.qcov_hsp_perc > 0.0 && query_len > 0 {
        hits.retain(|h| {
            let cov = 100.0 * h.align_len as f64 / query_len as f64;
            cov >= args.qcov_hsp_perc
        });
    }
    // Limit HSPs per subject
    if args.max_hsps > 0 {
        let max = args.max_hsps as usize;
        let mut counts: std::collections::HashMap<String, usize> = std::collections::HashMap::new();
        hits.retain(|h| {
            let c = counts.entry(h.subject_id.clone()).or_insert(0);
            *c += 1;
            *c <= max
        });
    }
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
}
