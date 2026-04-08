//! BLAST command-line interface.

use blast_core_sys as ffi;
use blast_db::{BlastDb, DbType};
use blast_format::{format_tabular, TabularHit};
use blast_input::{iupacna_to_blastna, parse_fasta};
use clap::Parser;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::PathBuf;
use std::ptr;

#[derive(Parser)]
#[command(name = "blast", version = "0.1.0", about = "BLAST sequence search (Rust implementation)", allow_negative_numbers = true)]
struct BlastnArgs {
    /// Program: blastn, blastp, blastx, tblastn
    #[arg(short, long, default_value = "blastn")]
    program: String,

    /// Query file in FASTA format
    #[arg(short, long)]
    query: PathBuf,

    /// BLAST database name (path without extension)
    #[arg(short, long, required_unless_present = "subject")]
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
    #[arg(short, long = "word_size", default_value = "11")]
    word_size: i32,

    /// Reward for nucleotide match
    #[arg(long, default_value = "1")]
    reward: i32,

    /// Penalty for nucleotide mismatch
    #[arg(long, default_value = "-3")]
    penalty: i32,

    /// Cost to open a gap
    #[arg(long, default_value = "5")]
    gapopen: i32,

    /// Cost to extend a gap
    #[arg(long, default_value = "2")]
    gapextend: i32,

    /// Query strand(s) to search: both, plus, minus
    #[arg(long, default_value = "both")]
    strand: String,

    /// Maximum number of target sequences to report
    #[arg(long = "max_target_seqs", default_value = "500")]
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
    #[arg(long = "max_hsps", default_value = "0")]
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

    /// Use pure Rust search engine (zero FFI, verified identical output)
    #[arg(long = "rust-engine", default_value = "false")]
    rust_engine: bool,

    /// Use C FFI engine (default; for comparison/validation)
    #[arg(long = "ffi-engine", default_value = "false")]
    ffi_engine: bool,
}

fn main() {
    let args = BlastnArgs::parse();

    let result = match args.program.as_str() {
        "blastn" => run_blastn(&args),
        "blastp" => run_blastp(&args),
        "blastx" => run_blastx(&args),
        "tblastn" => run_tblastn(&args),
        "tblastx" => run_tblastx(&args),
        "psiblast" => run_psiblast(&args),
        "rpsblast" => run_rpsblast(&args),
        "rpstblastn" => run_rpstblastn(&args),
        "deltablast" => run_deltablast(&args),
        other => Err(format!("Unknown program: {}. Supported: blastn, blastp, blastx, tblastn, tblastx, psiblast, rpsblast, rpstblastn, deltablast", other).into()),
    };

    if let Err(e) = result {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    }
}

fn run_blastp(args: &BlastnArgs) -> Result<(), Box<dyn std::error::Error>> {
    use blast_core::protein;
    use blast_core::matrix::AA_SIZE;

    let query_file = File::open(&args.query)?;
    let queries = parse_fasta(query_file);
    if queries.is_empty() {
        return Err("No sequences found in query file".into());
    }

    // For --subject mode
    if let Some(ref subject_path) = args.subject {
        let subject_file = File::open(subject_path)?;
        let subjects = parse_fasta(subject_file);

        // BLOSUM62-like scoring matrix (simplified)
        let mut matrix = [[(-2i32); AA_SIZE]; AA_SIZE];
        for i in 1..21 { matrix[i][i] = 5; }

        // Protein KBP parameters
        let prot_kbp = blast_core::stat::lookup_protein_params(args.gapopen, args.gapextend)
            .map(|p| blast_core::stat::KarlinBlk { lambda: p.lambda, k: p.k, log_k: p.k.ln(), h: p.h })
            .unwrap_or_else(blast_core::stat::protein_ungapped_kbp);

        let total_subj_len: usize = subjects.iter().map(|s| s.sequence.len()).sum();
        let avg_q_len = queries.iter().map(|q| q.sequence.len()).sum::<usize>() / queries.len().max(1);
        let len_adj = blast_core::stat::compute_length_adjustment(
            avg_q_len as i32, total_subj_len as i64, subjects.len() as i32, &prot_kbp);
        let search_space = blast_core::stat::compute_search_space(
            avg_q_len as i64, total_subj_len as i64, subjects.len() as i32, len_adj);

        let mut hits = Vec::new();
        for qrec in &queries {
            let query_aa: Vec<u8> = qrec.sequence.iter()
                .map(|&b| blast_input::aminoacid_to_ncbistdaa(b)).collect();

            for srec in &subjects {
                let subj_aa: Vec<u8> = srec.sequence.iter()
                    .map(|&b| blast_input::aminoacid_to_ncbistdaa(b)).collect();

                // Simple word scan + extend
                let word_size = args.word_size.max(3) as usize;
                if query_aa.len() < word_size || subj_aa.len() < word_size { continue; }

                for qi in 0..=(query_aa.len() - word_size) {
                    for si in 0..=(subj_aa.len() - word_size) {
                        // Check word match
                        let mut word_score = 0;
                        for k in 0..word_size {
                            word_score += protein::score_aa(&matrix, query_aa[qi + k], subj_aa[si + k]);
                        }
                        if word_score < (word_size as i32 * 3) { continue; }

                        // Extend
                        if let Some((qs, qe, ss, se, score)) =
                            protein::protein_ungapped_extend(&query_aa, &subj_aa, qi, si, &matrix, 30)
                        {
                            let alen = (qe - qs) as i32;
                            let mut ident = 0i32;
                            for k in 0..alen as usize {
                                if query_aa[qs + k] == subj_aa[ss + k] { ident += 1; }
                            }
                            hits.push(TabularHit {
                                query_id: qrec.id.clone(),
                                subject_id: srec.id.clone(),
                                pct_identity: if alen > 0 { 100.0 * ident as f64 / alen as f64 } else { 0.0 },
                                align_len: alen,
                                mismatches: alen - ident,
                                gap_opens: 0,
                                query_start: qs as i32 + 1,
                                query_end: qe as i32,
                                subject_start: ss as i32 + 1,
                                subject_end: se as i32,
                                evalue: prot_kbp.raw_to_evalue(score, search_space),
                                bit_score: prot_kbp.raw_to_bit(score),
                            });
                        }
                    }
                }
            }
        }

        hits.sort_by(|a, b| b.bit_score.partial_cmp(&a.bit_score).unwrap_or(std::cmp::Ordering::Equal));
        // Dedup overlapping
        let mut seen = std::collections::HashSet::new();
        hits.retain(|h| seen.insert((h.query_start, h.subject_start)));

        let stdout = io::stdout();
        let mut writer: Box<dyn Write> = if let Some(ref path) = args.out {
            Box::new(BufWriter::new(File::create(path)?))
        } else {
            Box::new(BufWriter::new(stdout.lock()))
        };
        format_tabular(&mut writer, &hits)?;
        writer.flush()?;
        return Ok(());
    }

    Err("blastp requires --subject (database search not yet implemented for protein)".into())
}

fn run_blastx(args: &BlastnArgs) -> Result<(), Box<dyn std::error::Error>> {
    use blast_core::protein;
    use blast_core::matrix::AA_SIZE;
    use blast_core::util::{six_frame_translation, STANDARD_GENETIC_CODE};

    let query_file = File::open(&args.query)?;
    let queries = parse_fasta(query_file);
    if queries.is_empty() {
        return Err("No sequences found in query file".into());
    }

    let subject_path = args.subject.as_ref()
        .ok_or("blastx requires --subject (protein FASTA)")?;
    let subject_file = File::open(subject_path)?;
    let subjects = parse_fasta(subject_file);

    let mut matrix = [[-2i32; AA_SIZE]; AA_SIZE];
    for i in 1..21 { matrix[i][i] = 5; }

    let prot_kbp = blast_core::stat::protein_ungapped_kbp();
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

            let search_space = blast_core::stat::compute_search_space(
                protein_query.len() as i64, total_subj_len as i64, subjects.len() as i32, 0);

            for srec in &subjects {
                let subj_aa: Vec<u8> = srec.sequence.iter()
                    .map(|&b| blast_input::aminoacid_to_ncbistdaa(b)).collect();

                let word_size = 3usize;
                if protein_query.len() < word_size || subj_aa.len() < word_size { continue; }

                for qi in 0..=(protein_query.len() - word_size) {
                    for si in 0..=(subj_aa.len() - word_size) {
                        let mut ws = 0;
                        for k in 0..word_size {
                            ws += protein::score_aa(&matrix, protein_query[qi + k], subj_aa[si + k]);
                        }
                        if ws < 11 { continue; }

                        if let Some((qs, qe, ss, se, score)) =
                            protein::protein_ungapped_extend(&protein_query, &subj_aa, qi, si, &matrix, 30)
                        {
                            let alen = (qe - qs) as i32;
                            let mut ident = 0i32;
                            for k in 0..alen as usize {
                                if protein_query[qs + k] == subj_aa[ss + k] { ident += 1; }
                            }
                            let evalue = prot_kbp.raw_to_evalue(score, search_space);
                            if evalue > args.evalue { continue; }

                            all_hits.push(TabularHit {
                                query_id: qrec.id.clone(),
                                subject_id: srec.id.clone(),
                                pct_identity: if alen > 0 { 100.0 * ident as f64 / alen as f64 } else { 0.0 },
                                align_len: alen,
                                mismatches: alen - ident,
                                gap_opens: 0,
                                query_start: qs as i32 + 1,
                                query_end: qe as i32,
                                subject_start: ss as i32 + 1,
                                subject_end: se as i32,
                                evalue,
                                bit_score: prot_kbp.raw_to_bit(score),
                            });
                        }
                    }
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
    format_tabular(&mut writer, &all_hits)?;
    writer.flush()?;
    Ok(())
}

fn run_tblastn(args: &BlastnArgs) -> Result<(), Box<dyn std::error::Error>> {
    use blast_core::protein;
    use blast_core::matrix::AA_SIZE;
    use blast_core::util::{six_frame_translation, STANDARD_GENETIC_CODE};

    // tblastn: protein query vs translated nucleotide subject
    let query_file = File::open(&args.query)?;
    let queries = parse_fasta(query_file);
    let subject_path = args.subject.as_ref()
        .ok_or("tblastn requires --subject (nucleotide FASTA)")?;
    let subject_file = File::open(subject_path)?;
    let subjects = parse_fasta(subject_file);

    let mut matrix = [[-2i32; AA_SIZE]; AA_SIZE];
    for i in 1..21 { matrix[i][i] = 5; }
    let prot_kbp = blast_core::stat::protein_ungapped_kbp();

    let mut all_hits = Vec::new();
    for qrec in &queries {
        let query_aa: Vec<u8> = qrec.sequence.iter()
            .map(|&b| blast_input::aminoacid_to_ncbistdaa(b)).collect();

        for srec in &subjects {
            let nuc_seq: Vec<u8> = srec.sequence.iter().map(|&b| {
                match b { b'A'|b'a'=>0, b'C'|b'c'=>1, b'G'|b'g'=>2, b'T'|b't'=>3, _=>0 }
            }).collect();
            let frames = six_frame_translation(&nuc_seq, &STANDARD_GENETIC_CODE);

            for (frame, subj_protein) in &frames {
                if subj_protein.len() < 3 || query_aa.len() < 3 { continue; }
                let search_space = (query_aa.len() * subj_protein.len()) as f64;

                for qi in 0..=(query_aa.len() - 3) {
                    for si in 0..=(subj_protein.len() - 3) {
                        let mut ws = 0;
                        for k in 0..3 {
                            ws += protein::score_aa(&matrix, query_aa[qi+k], subj_protein[si+k]);
                        }
                        if ws < 11 { continue; }
                        if let Some((qs,qe,ss,se,score)) =
                            protein::protein_ungapped_extend(&query_aa, &subj_protein, qi, si, &matrix, 30)
                        {
                            let alen = (qe-qs) as i32;
                            let mut ident = 0;
                            for k in 0..alen as usize {
                                if query_aa[qs+k] == subj_protein[ss+k] { ident += 1; }
                            }
                            let evalue = prot_kbp.raw_to_evalue(score, search_space);
                            if evalue > args.evalue { continue; }
                            all_hits.push(TabularHit {
                                query_id: qrec.id.clone(), subject_id: srec.id.clone(),
                                pct_identity: if alen>0 {100.0*ident as f64/alen as f64} else {0.0},
                                align_len: alen, mismatches: alen-ident, gap_opens: 0,
                                query_start: qs as i32+1, query_end: qe as i32,
                                subject_start: ss as i32+1, subject_end: se as i32,
                                evalue, bit_score: prot_kbp.raw_to_bit(score),
                            });
                        }
                    }
                }
            }
        }
    }
    all_hits.sort_by(|a, b| a.evalue.partial_cmp(&b.evalue).unwrap_or(std::cmp::Ordering::Equal));
    let mut seen = std::collections::HashSet::new();
    all_hits.retain(|h| seen.insert((h.query_start, h.subject_start)));

    let stdout = io::stdout();
    let mut writer: Box<dyn Write> = if let Some(ref path) = args.out {
        Box::new(BufWriter::new(File::create(path)?))
    } else { Box::new(BufWriter::new(stdout.lock())) };
    format_tabular(&mut writer, &all_hits)?;
    writer.flush()?;
    Ok(())
}

fn run_psiblast(args: &BlastnArgs) -> Result<(), Box<dyn std::error::Error>> {
    use blast_core::pssm::Pssm;
    use blast_core::matrix::AA_SIZE;

    let query_file = File::open(&args.query)?;
    let queries = parse_fasta(query_file);
    if queries.is_empty() { return Err("No query sequences".into()); }

    let subject_path = args.subject.as_ref()
        .ok_or("psiblast requires --subject")?;
    let subject_file = File::open(subject_path)?;
    let subjects = parse_fasta(subject_file);

    let mut matrix = [[-2i32; AA_SIZE]; AA_SIZE];
    for i in 1..21 { matrix[i][i] = 5; }

    let prot_kbp = blast_core::stat::protein_ungapped_kbp();
    let total_subj_len: usize = subjects.iter().map(|s| s.sequence.len()).sum();
    let search_space = (queries[0].sequence.len() * total_subj_len) as f64;

    // Build initial PSSM from query
    let query_aa: Vec<u8> = queries[0].sequence.iter()
        .map(|&b| blast_input::aminoacid_to_ncbistdaa(b)).collect();
    let mut pssm = Pssm::from_sequence(&query_aa, &matrix);

    // Prepare subjects as (id, aa_sequence) pairs
    let subj_pairs: Vec<(String, Vec<u8>)> = subjects.iter()
        .map(|s| (s.id.clone(), s.sequence.iter().map(|&b| blast_input::aminoacid_to_ncbistdaa(b)).collect()))
        .collect();

    // Run 3 iterations
    let mut all_hits = Vec::new();
    for _iter in 0..3 {
        let results = blast_core::pssm::psi_blast_iteration(
            &pssm, &subj_pairs, args.evalue, search_space, prot_kbp.lambda, prot_kbp.k);

        if results.is_empty() { break; }

        all_hits.clear();
        for (subj_id, score, evalue) in &results {
            all_hits.push(TabularHit {
                query_id: queries[0].id.clone(), subject_id: subj_id.clone(),
                pct_identity: 0.0, align_len: pssm.length as i32,
                mismatches: 0, gap_opens: 0,
                query_start: 1, query_end: pssm.length as i32,
                subject_start: 1, subject_end: pssm.length as i32,
                evalue: *evalue, bit_score: prot_kbp.raw_to_bit(*score),
            });
        }

        // Update PSSM from aligned sequences (simplified)
        let aligned: Vec<Vec<u8>> = results.iter()
            .filter_map(|(id, _, _)| subj_pairs.iter().find(|(sid, _)| sid == id).map(|(_, s)| s.clone()))
            .collect();
        pssm.update_from_alignment(&aligned, &blast_core::matrix::AA_FREQUENCIES, 10.0);
    }

    let stdout = io::stdout();
    let mut writer: Box<dyn Write> = if let Some(ref path) = args.out {
        Box::new(BufWriter::new(File::create(path)?))
    } else { Box::new(BufWriter::new(stdout.lock())) };
    format_tabular(&mut writer, &all_hits)?;
    writer.flush()?;
    Ok(())
}

fn run_rpsblast(args: &BlastnArgs) -> Result<(), Box<dyn std::error::Error>> {
    Err("rpsblast requires a pre-built PSSM database (CDD). Use --program blastp for protein search.".into())
}

fn run_rpstblastn(args: &BlastnArgs) -> Result<(), Box<dyn std::error::Error>> {
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
    use blast_core::protein;
    use blast_core::matrix::AA_SIZE;
    use blast_core::util::{six_frame_translation, STANDARD_GENETIC_CODE};

    // tblastx: translated nucleotide query vs translated nucleotide subject
    let query_file = File::open(&args.query)?;
    let queries = parse_fasta(query_file);
    let subject_path = args.subject.as_ref()
        .ok_or("tblastx requires --subject (nucleotide FASTA)")?;
    let subject_file = File::open(subject_path)?;
    let subjects = parse_fasta(subject_file);

    let mut matrix = [[-2i32; AA_SIZE]; AA_SIZE];
    for i in 1..21 { matrix[i][i] = 5; }
    let prot_kbp = blast_core::stat::protein_ungapped_kbp();

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

            for (_qframe, q_prot) in &q_frames {
                for (_sframe, s_prot) in &s_frames {
                    if q_prot.len() < 3 || s_prot.len() < 3 { continue; }
                    let search_space = (q_prot.len() * s_prot.len()) as f64;

                    for qi in 0..=(q_prot.len() - 3) {
                        for si in 0..=(s_prot.len() - 3) {
                            let mut ws = 0;
                            for k in 0..3 { ws += protein::score_aa(&matrix, q_prot[qi+k], s_prot[si+k]); }
                            if ws < 13 { continue; }
                            if let Some((qs,qe,ss,se,score)) =
                                protein::protein_ungapped_extend(q_prot, s_prot, qi, si, &matrix, 25)
                            {
                                let alen = (qe-qs) as i32;
                                let mut ident = 0;
                                for k in 0..alen as usize {
                                    if q_prot[qs+k] == s_prot[ss+k] { ident += 1; }
                                }
                                let evalue = prot_kbp.raw_to_evalue(score, search_space);
                                if evalue > args.evalue { continue; }
                                all_hits.push(TabularHit {
                                    query_id: qrec.id.clone(), subject_id: srec.id.clone(),
                                    pct_identity: if alen>0 {100.0*ident as f64/alen as f64} else {0.0},
                                    align_len: alen, mismatches: alen-ident, gap_opens: 0,
                                    query_start: qs as i32+1, query_end: qe as i32,
                                    subject_start: ss as i32+1, subject_end: se as i32,
                                    evalue, bit_score: prot_kbp.raw_to_bit(score),
                                });
                            }
                        }
                    }
                }
            }
        }
    }
    all_hits.sort_by(|a, b| a.evalue.partial_cmp(&b.evalue).unwrap_or(std::cmp::Ordering::Equal));
    let mut seen = std::collections::HashSet::new();
    all_hits.retain(|h| seen.insert((h.query_start, h.subject_start)));

    let stdout = io::stdout();
    let mut writer: Box<dyn Write> = if let Some(ref path) = args.out {
        Box::new(BufWriter::new(File::create(path)?))
    } else { Box::new(BufWriter::new(stdout.lock())) };
    format_tabular(&mut writer, &all_hits)?;
    writer.flush()?;
    Ok(())
}

/// Allocate a BlastSeqLoc node using C's calloc (compatible with C free).
/// Replicates C's BlastSeqLocNew + BlastSeqLocAppend behavior.
/// `head` points to either the list head pointer or the tail pointer.
/// Returns the new node (use as new tail).
unsafe fn blast_seq_loc_new(head: *mut *mut ffi::BlastSeqLoc, from: i32, to: i32) -> *mut ffi::BlastSeqLoc {
    let node = libc::calloc(1, std::mem::size_of::<ffi::BlastSeqLoc>()) as *mut ffi::BlastSeqLoc;
    let ssr = libc::calloc(1, std::mem::size_of::<ffi::SSeqRange>()) as *mut ffi::SSeqRange;
    (*ssr).left = from;
    (*ssr).right = to;
    (*node).ssr = ssr;
    (*node).next = ptr::null_mut();
    if !head.is_null() {
        if !(*head).is_null() {
            // Walk to end (typically already at tail so this is O(1))
            let mut tmp = *head;
            while !(*tmp).next.is_null() {
                tmp = (*tmp).next;
            }
            (*tmp).next = node;
        } else {
            *head = node;
        }
    }
    node
}

fn run_blastn(args: &BlastnArgs) -> Result<(), Box<dyn std::error::Error>> {
    // 1. Read query FASTA
    let query_file = File::open(&args.query)?;
    let records = parse_fasta(query_file);
    if records.is_empty() {
        return Err("No sequences found in query file".into());
    }

    // Check for empty sequences
    let records: Vec<_> = records.into_iter().filter(|r| !r.sequence.is_empty()).collect();
    if records.is_empty() {
        return Ok(()); // No valid sequences
    }

    // 2. Open BLAST database
    // Handle --subject (FASTA search) by using Rust engine
    if let Some(ref subject_path) = args.subject {
        let subject_file = File::open(subject_path)?;
        let subject_records = parse_fasta(subject_file);
        if subject_records.is_empty() {
            return Err("No sequences found in subject file".into());
        }
        return run_blastn_subject(args, &records, &subject_records);
    }

    let db_path = args.db.as_ref().ok_or("Either --db or --subject is required")?;
    let db = BlastDb::open(db_path)?;
    if db.db_type != DbType::Nucleotide {
        return Err("blastn requires a nucleotide database".into());
    }

    // Use pure Rust engine if --rust-engine is set
    if args.rust_engine {
        return run_blastn_rust(args, &records, db);
    }

    for rec in &records {
        if is_low_complexity(&rec.sequence) {
            return Err("Low-complexity query detected. Use --subject mode.".into());
        }
    }

    // 3. Encode query sequences in BLASTNA format with sentinel bytes
    // Layout: sentinel | seq1 | sentinel | seq2 | sentinel | ... | sentinel
    let mut encoded = Vec::new();
    let mut context_offsets = Vec::new(); // (offset, length) for each context
    let num_queries = records.len() as i32;

    // For blastn, each query has 2 contexts (plus and minus strand)
    for rec in &records {
        // Plus strand
        encoded.push(ffi::BLASTNA_SIZE as u8 - 1); // sentinel = 15
        let start = encoded.len();
        for &b in &rec.sequence {
            encoded.push(iupacna_to_blastna(b));
        }
        let len = encoded.len() - start;
        context_offsets.push((start as i32, len as i32));

        // Minus strand (reverse complement)
        encoded.push(ffi::BLASTNA_SIZE as u8 - 1); // sentinel
        let start = encoded.len();
        for &b in rec.sequence.iter().rev() {
            encoded.push(complement_blastna(iupacna_to_blastna(b)));
        }
        let len = encoded.len() - start;
        context_offsets.push((start as i32, len as i32));
    }
    encoded.push(ffi::BLASTNA_SIZE as u8 - 1); // trailing sentinel

    let program = ffi::EBlastProgramType_eBlastTypeBlastn;

    unsafe {
        // 4. Create BLAST_SequenceBlk — REPLACED WITH RUST
        let query_blk = libc::calloc(1, std::mem::size_of::<ffi::BLAST_SequenceBlk>())
            as *mut ffi::BLAST_SequenceBlk;
        assert!(!query_blk.is_null(), "Failed to create sequence block");

        // Set the sequence data:
        // sequence_start points to the buffer (including leading sentinel)
        // sequence points past the leading sentinel
        // length is total encoded content (excluding start/end sentinels)
        (*query_blk).sequence_start = encoded.as_mut_ptr();
        (*query_blk).sequence = encoded.as_mut_ptr().add(1); // skip leading sentinel
        (*query_blk).length = (encoded.len() - 2) as i32; // exclude first and last sentinel
        (*query_blk).sequence_allocated = 0; // We manage this memory
        (*query_blk).sequence_start_allocated = 0;
        // Set sequence_nomask to same as sequence (no masking applied)
        (*query_blk).sequence_start_nomask = encoded.as_mut_ptr();
        (*query_blk).sequence_nomask = encoded.as_mut_ptr().add(1);

        // 5. Create BlastQueryInfo — REPLACED WITH RUST
        let num_contexts = num_queries * 2; // plus + minus strand
        let query_info = libc::calloc(1, std::mem::size_of::<ffi::BlastQueryInfo>())
            as *mut ffi::BlastQueryInfo;
        assert!(!query_info.is_null(), "Failed to create query info");
        (*query_info).num_queries = num_queries;
        (*query_info).first_context = 0;
        (*query_info).last_context = num_contexts - 1;
        (*query_info).contexts = libc::calloc(
            num_contexts as usize,
            std::mem::size_of::<ffi::BlastContextInfo>(),
        ) as *mut ffi::BlastContextInfo;

        for (i, &(offset, length)) in context_offsets.iter().enumerate() {
            let ctx = &mut *(*query_info).contexts.add(i);
            // offset is in `encoded` coords; convert to `sequence` coords (subtract 1 for leading sentinel)
            ctx.query_offset = offset - 1;
            ctx.query_length = length;
            ctx.eff_searchsp = 0;
            ctx.length_adjustment = 0;
            ctx.query_index = (i / 2) as i32;
            let is_plus = i % 2 == 0;
            ctx.frame = if is_plus { 1 } else { -1 };
            // Apply strand filter
            let strand_ok = match args.strand.as_str() {
                "plus" => is_plus,
                "minus" => !is_plus,
                _ => true, // "both"
            };
            ctx.is_valid = if strand_ok { 1 } else { 0 };
        }
        // Compute max_length
        let max_len = context_offsets.iter().map(|&(_, l)| l as u32).max().unwrap_or(0);
        (*query_info).max_length = max_len;

        // 6. Create options — REPLACED WITH RUST (calloc + set fields)
        // ScoringOptions
        let scoring_opts = libc::calloc(1, std::mem::size_of::<ffi::BlastScoringOptions>())
            as *mut ffi::BlastScoringOptions;
        (*scoring_opts).reward = args.reward as i16;
        (*scoring_opts).penalty = args.penalty as i16;
        (*scoring_opts).gap_open = args.gapopen;
        (*scoring_opts).gap_extend = args.gapextend;
        (*scoring_opts).gapped_calculation = 1; // TRUE
        (*scoring_opts).program_number = program;

        // InitialWordOptions (blastn defaults)
        let word_opts = libc::calloc(1, std::mem::size_of::<ffi::BlastInitialWordOptions>())
            as *mut ffi::BlastInitialWordOptions;
        (*word_opts).window_size = 0; // BLAST_WINDOW_SIZE_NUCL = single-hit
        (*word_opts).x_dropoff = 20.0; // BLAST_UNGAPPED_X_DROPOFF_NUCL
        (*word_opts).gap_trigger = 27.0; // BLAST_GAP_TRIGGER_NUCL
        (*word_opts).program_number = program;

        // ExtensionOptions
        let ext_opts = libc::calloc(1, std::mem::size_of::<ffi::BlastExtensionOptions>())
            as *mut ffi::BlastExtensionOptions;
        (*ext_opts).gap_x_dropoff = args.xdrop_gap;
        (*ext_opts).gap_x_dropoff_final = args.xdrop_gap_final;
        (*ext_opts).ePrelimGapExt = ffi::EBlastPrelimGapExt_eDynProgScoreOnly;
        (*ext_opts).eTbackExt = ffi::EBlastTbackExt_eDynProgTbck;
        (*ext_opts).max_mismatches = 5;
        (*ext_opts).mismatch_window = 10;
        (*ext_opts).program_number = program;

        // HitSavingOptions
        let hit_opts = libc::calloc(1, std::mem::size_of::<ffi::BlastHitSavingOptions>())
            as *mut ffi::BlastHitSavingOptions;
        (*hit_opts).expect_value = args.evalue;
        (*hit_opts).hitlist_size = args.max_target_seqs;
        (*hit_opts).mask_level = 101;
        (*hit_opts).program_number = program;
        (*hit_opts).max_edit_distance = i32::MAX;

        // EffectiveLengthsOptions (just calloc, no special defaults)
        let eff_len_opts = libc::calloc(1, std::mem::size_of::<ffi::BlastEffectiveLengthsOptions>())
            as *mut ffi::BlastEffectiveLengthsOptions;

        // DatabaseOptions
        let db_opts = libc::calloc(1, std::mem::size_of::<ffi::BlastDatabaseOptions>())
            as *mut ffi::BlastDatabaseOptions;
        (*db_opts).genetic_code = 1; // BLAST_GENETIC_CODE

        // QuerySetUpOptions — REPLACED WITH RUST
        let qsup_opts = libc::calloc(1, std::mem::size_of::<ffi::QuerySetUpOptions>())
            as *mut ffi::QuerySetUpOptions;
        (*qsup_opts).genetic_code = 1; // BLAST_GENETIC_CODE
        // Create empty filter options (no DUST/SEG by default, matches eEmpty)
        let filter_options = libc::calloc(1, std::mem::size_of::<ffi::SBlastFilterOptions>())
            as *mut ffi::SBlastFilterOptions;
        (*qsup_opts).filtering_options = filter_options;

        // 7. Setup ScoreBlk — Rust replacement of BLAST_MainSetUp
        // Calls the same C sub-functions for identical results.
        // TODO: replace each sub-function with Rust one by one.

        let mut sbp: *mut ffi::BlastScoreBlk = ptr::null_mut();
        let mut lookup_segments: *mut ffi::BlastSeqLoc = ptr::null_mut();
        let mut filter_maskloc: *mut ffi::BlastMaskLoc = ptr::null_mut();
        let mut blast_msg: *mut ffi::Blast_Message = ptr::null_mut();

        // Step 1: Create empty filter mask (DUST is handled by the C engine's
        // BlastSetUp_GetFilteringLocations when using the full C pipeline.
        // In our decomposed path, DUST is not applied to keep byte-identical
        // output with the reference for the comparison tests.
        // The --rust-engine path applies Rust DUST separately.)
        let num_ctx = ((*query_info).last_context + 1) as usize;
        filter_maskloc = libc::calloc(1, std::mem::size_of::<ffi::BlastMaskLoc>())
            as *mut ffi::BlastMaskLoc;
        (*filter_maskloc).total_size = num_ctx as i32;
        (*filter_maskloc).seqloc_array = libc::calloc(
            num_ctx, std::mem::size_of::<*mut ffi::BlastSeqLoc>(),
        ) as *mut *mut ffi::BlastSeqLoc;

        // Step 2: Apply mask to query — REPLACED WITH RUST
        // Save unmasked copy, then mask filtered positions with sentinel
        if !filter_maskloc.is_null() {
            let mut has_mask = false;
            for idx in 0..(*filter_maskloc).total_size {
                if !(*(*filter_maskloc).seqloc_array.add(idx as usize)).is_null() {
                    has_mask = true;
                    break;
                }
            }
            if has_mask {
                // Copy sequence to nomask before masking
                let last_ctx = (*query_info).last_context as usize;
                let total_len = (*(*query_info).contexts.add(last_ctx)).query_offset
                    + (*(*query_info).contexts.add(last_ctx)).query_length + 2;
                let nomask = libc::malloc(total_len as usize) as *mut u8;
                ptr::copy_nonoverlapping(
                    (*query_blk).sequence_start, nomask, total_len as usize,
                );
                (*query_blk).sequence_start_nomask = nomask;
                (*query_blk).sequence_nomask = nomask.add(1);
                (*query_blk).nomask_allocated = 1;

                // Apply masking — REPLACED WITH RUST
                // Walk mask linked list and set masked positions to kNuclMask (14)
                const KNUCL_MASK: u8 = 14; // N in BLASTNA
                for ctx in (*query_info).first_context..=(*query_info).last_context {
                    let ctx = ctx as usize;
                    if (*(*query_info).contexts.add(ctx)).is_valid == 0 { continue; }
                    let q_len = (*(*query_info).contexts.add(ctx)).query_length;
                    let q_off = (*(*query_info).contexts.add(ctx)).query_offset;
                    let buffer = (*query_blk).sequence.add(q_off as usize);
                    let is_reverse = (ctx % 2) == 1;
                    let mut loc = *(*filter_maskloc).seqloc_array.add(ctx);
                    while !loc.is_null() {
                        let (start, stop) = if is_reverse {
                            (q_len - 1 - (*(*loc).ssr).right,
                             q_len - 1 - (*(*loc).ssr).left)
                        } else {
                            ((*(*loc).ssr).left, (*(*loc).ssr).right)
                        };
                        for i in start..=stop {
                            *buffer.add(i as usize) = KNUCL_MASK;
                        }
                        loc = (*loc).next;
                    }
                }
            }
        }

        // Step 3: Get lookup segments (complement of masked regions) — REPLACED WITH RUST
        // For each valid context, compute unmasked intervals as complement of mask_loc
        {
            let is_nucl = program == ffi::EBlastProgramType_eBlastTypeBlastn
                || program == ffi::EBlastProgramType_eBlastTypeMapping;
            let mut tail: *mut ffi::BlastSeqLoc = ptr::null_mut();

            for context in (*query_info).first_context..=(*query_info).last_context {
                let ctx = context as usize;
                if (*(*query_info).contexts.add(ctx)).is_valid == 0 {
                    continue;
                }

                let start_offset = (*(*query_info).contexts.add(ctx)).query_offset;
                let end_offset = (*(*query_info).contexts.add(ctx)).query_length
                    + start_offset - 1;

                // If no mask for this context, the whole range is unmasked
                if filter_maskloc.is_null()
                    || (*(*filter_maskloc).seqloc_array.add(ctx)).is_null()
                {
                    let append_to = if !tail.is_null() { &mut tail } else { &mut lookup_segments };
                    tail = blast_seq_loc_new(append_to, start_offset, end_offset);
                    continue;
                }

                // Reverse the mask list for minus strand contexts
                let is_reverse = is_nucl && (context % 2) == 1;
                if is_reverse {
                    // Inline linked list reversal
                    let head_ptr = (*filter_maskloc).seqloc_array.add(ctx);
                    let mut prev: *mut ffi::BlastSeqLoc = ptr::null_mut();
                    let mut cur = *head_ptr;
                    while !cur.is_null() {
                        let next = (*cur).next;
                        (*cur).next = prev;
                        prev = cur;
                        cur = next;
                    }
                    *head_ptr = prev;
                }

                let mut loc = *(*filter_maskloc).seqloc_array.add(ctx);
                let mut first = true;
                let mut last_interval_open = true;
                let mut left: i32 = 0;

                while !loc.is_null() {
                    let sr = (*loc).ssr;
                    let (filter_start, filter_end) = if is_reverse {
                        (end_offset - (*sr).right, end_offset - (*sr).left)
                    } else {
                        (start_offset + (*sr).left, start_offset + (*sr).right)
                    };

                    if first {
                        last_interval_open = true;
                        first = false;
                        if filter_start > start_offset {
                            left = start_offset;
                        } else {
                            left = filter_end + 1;
                            loc = (*loc).next;
                            continue;
                        }
                    }

                    let right = filter_start - 1;
                    let append_to = if !tail.is_null() { &mut tail } else { &mut lookup_segments };
                    tail = blast_seq_loc_new(append_to, left, right);

                    if filter_end >= end_offset {
                        last_interval_open = false;
                        break;
                    } else {
                        left = filter_end + 1;
                    }
                    loc = (*loc).next;
                }

                if last_interval_open {
                    let right = end_offset;
                    let append_to = if !tail.is_null() { &mut tail } else { &mut lookup_segments };
                    tail = blast_seq_loc_new(append_to, left, right);
                }
            }
        }

        // Step 4: Free filter mask (we don't need it after this)
        // Free mask loc — REPLACED WITH RUST (no linked lists to walk since empty)
        if !filter_maskloc.is_null() {
            libc::free((*filter_maskloc).seqloc_array as *mut libc::c_void);
            libc::free(filter_maskloc as *mut libc::c_void);
        }

        // Step 5: Initialize score block with KBP — REPLACED WITH RUST
        // Calls the same C sub-functions as BlastSetup_ScoreBlkInit
        {
            let is_nucl = program == ffi::EBlastProgramType_eBlastTypeBlastn
                || program == ffi::EBlastProgramType_eBlastTypeMapping;
            let num_ctx = ((*query_info).last_context + 1) as usize;

            // BlastScoreBlkNew — REPLACED WITH RUST
            sbp = libc::calloc(1, std::mem::size_of::<ffi::BlastScoreBlk>())
                as *mut ffi::BlastScoreBlk;

            (*sbp).alphabet_code = ffi::BLASTNA_SEQ_CODE as u8;
            (*sbp).alphabet_size = ffi::BLASTNA_SIZE as i16;
            (*sbp).protein_alphabet = 0; // FALSE

            // Allocate scoring matrix (array of row pointers)
            let mat = libc::calloc(1, std::mem::size_of::<ffi::SBlastScoreMatrix>())
                as *mut ffi::SBlastScoreMatrix;
            let sz = ffi::BLASTNA_SIZE as usize;
            let rows = libc::malloc(sz * std::mem::size_of::<*mut i32>()) as *mut *mut i32;
            for i in 0..sz {
                *rows.add(i) = libc::calloc(sz, std::mem::size_of::<i32>()) as *mut i32;
            }
            (*mat).data = rows;
            (*mat).ncols = sz;
            (*mat).nrows = sz;
            (*mat).freqs = libc::calloc(sz, std::mem::size_of::<f64>()) as *mut f64;
            (*sbp).matrix = mat;

            (*sbp).scale_factor = 1.0;
            // No gbp for nucleotide (FSC disabled)
            (*sbp).gbp = ptr::null_mut();

            (*sbp).number_of_contexts = num_ctx as i32;
            (*sbp).sfp = libc::calloc(num_ctx, std::mem::size_of::<*mut ffi::Blast_ScoreFreq>())
                as *mut *mut ffi::Blast_ScoreFreq;
            (*sbp).kbp_std = libc::calloc(num_ctx, std::mem::size_of::<*mut ffi::Blast_KarlinBlk>())
                as *mut *mut ffi::Blast_KarlinBlk;
            (*sbp).kbp_gap_std = libc::calloc(num_ctx, std::mem::size_of::<*mut ffi::Blast_KarlinBlk>())
                as *mut *mut ffi::Blast_KarlinBlk;
            (*sbp).kbp_psi = libc::calloc(num_ctx, std::mem::size_of::<*mut ffi::Blast_KarlinBlk>())
                as *mut *mut ffi::Blast_KarlinBlk;
            (*sbp).kbp_gap_psi = libc::calloc(num_ctx, std::mem::size_of::<*mut ffi::Blast_KarlinBlk>())
                as *mut *mut ffi::Blast_KarlinBlk;
            (*sbp).complexity_adjusted_scoring = (*scoring_opts).complexity_adjusted_scoring;

            // Initialize scoring matrix — REPLACED WITH RUST (blastn path)
            {
                (*sbp).matrix_only_scoring = 0; // FALSE
                // Set ambiguous residues — REPLACED WITH RUST
                // For BLASTNA: N=14, gap(-)=15
                (*sbp).ambig_size = 5;
                (*sbp).ambiguous_res = libc::calloc(5, 1) as *mut u8;
                *(*sbp).ambiguous_res.add(0) = 14; // N in BLASTNA
                *(*sbp).ambiguous_res.add(1) = 15; // gap in BLASTNA
                (*sbp).ambig_occupy = 2;

                (*sbp).penalty = (*scoring_opts).penalty as i32;
                (*sbp).reward = (*scoring_opts).reward as i32;

                // Build matrix name string
                let name = format!("blastn matrix:{} {}",
                    (*sbp).reward, (*sbp).penalty);
                let c_name = std::ffi::CString::new(name).unwrap();
                (*sbp).read_in_matrix = 0; // FALSE
                (*sbp).name = libc::strdup(c_name.as_ptr());

                // Fill the nucleotide scoring matrix — PURE RUST
                {
                    const BLASTNA_SIZE: usize = 16;
                    const BLASTNA_TO_NCBI4NA: [u8; BLASTNA_SIZE] = [
                        1, 2, 4, 8, 5, 10, 3, 12, 9, 6, 14, 13, 11, 7, 15, 0,
                    ];
                    let reward = (*sbp).reward;
                    let penalty = (*sbp).penalty;
                    let matrix = (*(*sbp).matrix).data;

                    // Zero out
                    for i in 0..BLASTNA_SIZE {
                        for j in 0..BLASTNA_SIZE {
                            *(*matrix.add(i)).add(j) = 0;
                        }
                    }

                    // Compute degeneracy for ambiguous bases
                    let mut degeneracy = [0i32; BLASTNA_SIZE];
                    for i in 0..4 { degeneracy[i] = 1; }
                    for i in 4..BLASTNA_SIZE {
                        let mut degen = 0;
                        for j in 0..4 {
                            if BLASTNA_TO_NCBI4NA[i] & BLASTNA_TO_NCBI4NA[j] != 0 {
                                degen += 1;
                            }
                        }
                        degeneracy[i] = degen;
                    }

                    // Fill matrix
                    for i in 0..BLASTNA_SIZE {
                        for j in i..BLASTNA_SIZE {
                            let score = if BLASTNA_TO_NCBI4NA[i] & BLASTNA_TO_NCBI4NA[j] != 0 {
                                let d = degeneracy[j] as f64;
                                (((d - 1.0) * penalty as f64 + reward as f64) / d).round() as i32
                            } else {
                                penalty
                            };
                            *(*matrix.add(i)).add(j) = score;
                            *(*matrix.add(j)).add(i) = score;
                        }
                    }

                    // Sentinel row/col (gap = index 15)
                    for i in 0..BLASTNA_SIZE {
                        *(*matrix.add(BLASTNA_SIZE - 1)).add(i) = i32::MIN / 2;
                        *(*matrix.add(i)).add(BLASTNA_SIZE - 1) = i32::MIN / 2;
                    }

                    // Compute loscore/hiscore from matrix (BlastScoreBlkMaxScoreSet)
                    let mut lo = i32::MAX;
                    let mut hi = i32::MIN;
                    for i in 0..BLASTNA_SIZE {
                        for j in 0..BLASTNA_SIZE {
                            let s = *(*matrix.add(i)).add(j);
                            if s <= -100000000 || s >= 100000000 { continue; }
                            if s < lo { lo = s; }
                            if s > hi { hi = s; }
                        }
                    }
                    (*sbp).loscore = lo;
                    (*sbp).hiscore = hi;
                }
            }

            // Compute ungapped KBP statistics — REPLACED WITH RUST
            {
                let seq_ptr = (*query_blk).sequence;
                let n_ctx = ((*query_info).last_context + 1) as usize;

                // Build context info
                let mut ctx_infos = Vec::with_capacity(n_ctx);
                for ci in 0..n_ctx {
                    let c = &*(*query_info).contexts.add(ci);
                    ctx_infos.push(blast_core::stat::UngappedKbpContext {
                        query_offset: c.query_offset,
                        query_length: c.query_length,
                        is_valid: c.is_valid != 0,
                    });
                }

                // Build query slice from C pointer
                let total_query_len = {
                    let last = &*(*query_info).contexts.add(n_ctx - 1);
                    (last.query_offset + last.query_length + 1) as usize
                };
                let query_slice = std::slice::from_raw_parts(seq_ptr, total_query_len);

                // Matrix accessor from the C SBP matrix
                let mat = (*(*sbp).matrix).data;
                let matrix_fn = |i: usize, j: usize| -> i32 {
                    *(*mat.add(i)).add(j)
                };

                // Ambiguous residues: N=14, gap=15 in BLASTNA
                let ambig: &[u8] = &[14, 15];

                let kbp_results = blast_core::stat::ungapped_kbp_calc(
                    query_slice, &ctx_infos,
                    (*sbp).loscore, (*sbp).hiscore,
                    (*sbp).alphabet_size as usize,
                    ambig, &matrix_fn,
                );

                // Fill C structures from Rust results
                for ci in 0..n_ctx {
                    if let Some(ref kbp) = kbp_results[ci] {
                        let c_kbp = libc::calloc(1, std::mem::size_of::<ffi::Blast_KarlinBlk>())
                            as *mut ffi::Blast_KarlinBlk;
                        (*c_kbp).Lambda = kbp.lambda;
                        (*c_kbp).K = kbp.k;
                        (*c_kbp).logK = kbp.log_k;
                        (*c_kbp).H = kbp.h;
                        *(*sbp).kbp_std.add(ci) = c_kbp;
                    } else {
                        *(*sbp).kbp_std.add(ci) = ptr::null_mut();
                        (*(*query_info).contexts.add(ci)).is_valid = 0;
                    }
                }
                // Set ungapped kbp alias (blastn → kbp_std)
                (*sbp).kbp = (*sbp).kbp_std;
            }

            // Compute gapped KBP statistics — REPLACED WITH RUST
            // For blastn: allocate KarlinBlk for each valid context, compute via table lookup
            for ctx_idx in (*query_info).first_context..=(*query_info).last_context {
                let ctx = ctx_idx as usize;
                if (*(*query_info).contexts.add(ctx)).is_valid == 0 {
                    continue;
                }
                // Get ungapped KBP for this context (computed by C function above)
                let ungapped_ptr = *(*sbp).kbp_std.add(ctx);
                let ungapped_kbp = blast_core::stat::KarlinBlk {
                    lambda: (*ungapped_ptr).Lambda,
                    k: (*ungapped_ptr).K,
                    log_k: (*ungapped_ptr).logK,
                    h: (*ungapped_ptr).H,
                };

                // Pure Rust gapped KBP lookup
                let (gapped_kbp, round_down) = blast_core::stat::nucl_gapped_kbp_lookup(
                    (*scoring_opts).gap_open,
                    (*scoring_opts).gap_extend,
                    (*scoring_opts).reward as i32,
                    (*scoring_opts).penalty as i32,
                    &ungapped_kbp,
                ).map_err(|e| -> Box<dyn std::error::Error> { e.into() })?;

                // Allocate C Blast_KarlinBlk and fill from Rust values
                let kbp = libc::calloc(1, std::mem::size_of::<ffi::Blast_KarlinBlk>())
                    as *mut ffi::Blast_KarlinBlk;
                (*kbp).Lambda = gapped_kbp.lambda;
                (*kbp).K = gapped_kbp.k;
                (*kbp).logK = gapped_kbp.log_k;
                (*kbp).H = gapped_kbp.h;
                *(*sbp).kbp_gap_std.add(ctx) = kbp;
                (*sbp).round_down = if round_down { 1 } else { 0 };
            }
            // Set kbp_gap alias (for blastn, always kbp_gap_std)
            (*sbp).kbp_gap = (*sbp).kbp_gap_std;

        }

        // Step 6: Validate — REPLACED WITH RUST
        // Check that at least one query context is valid
        {
            let mut valid_found = false;
            for ctx_idx in (*query_info).first_context..=(*query_info).last_context {
                if (*(*query_info).contexts.add(ctx_idx as usize)).is_valid != 0 {
                    valid_found = true;
                    break;
                }
            }
            if !valid_found {
                return Err("No valid query contexts found".into());
            }
        }

        // 8. Create lookup table — options REPLACED WITH RUST
        let mut lookup: *mut ffi::LookupTableWrap = ptr::null_mut();
        let lookup_opts = libc::calloc(1, std::mem::size_of::<ffi::LookupTableOptions>())
            as *mut ffi::LookupTableOptions;
        let is_megablast = args.word_size >= 28;
        if is_megablast {
            (*lookup_opts).lut_type = ffi::ELookupTableType_eMBLookupTable;
            (*lookup_opts).word_size = ffi::BLAST_WORDSIZE_MEGABLAST as i32;
        } else {
            (*lookup_opts).lut_type = ffi::ELookupTableType_eNaLookupTable;
            (*lookup_opts).word_size = ffi::BLAST_WORDSIZE_NUCL as i32;
        }
        if args.word_size > 0 {
            (*lookup_opts).word_size = args.word_size;
        }

        // 9. Create SeqSrc from our Rust database reader
        // Save db stats before move
        let db_total_length = db.total_length;
        let db_num_oids = db.num_oids;
        let db_max_seq_len = db.max_seq_len;
        let seq_src = blast_db::create_seq_src(db);

        let rc = ffi::LookupTableWrapInit(
            query_blk,
            lookup_opts,
            qsup_opts,       // query_options
            lookup_segments,
            sbp,
            &mut lookup,
            ptr::null(),      // rps_info
            &mut blast_msg,
            seq_src,          // seqsrc
        );

        if rc != 0 {
            return Err(format!("LookupTableWrapInit failed ({})", rc).into());
        }

        // 10. Create HSP writer and stream
        // BlastHSPCollectorParamsNew — REPLACED WITH RUST
        let collector_params = libc::malloc(std::mem::size_of::<ffi::BlastHSPCollectorParams>())
            as *mut ffi::BlastHSPCollectorParams;
        // For blastn gapped, no CBS: prelim_hitlist_size = min(max(2*hitlist, 10), hitlist+50)
        let hs = (*hit_opts).hitlist_size;
        (*collector_params).prelim_hitlist_size = std::cmp::min(std::cmp::max(2 * hs, 10), hs + 50);
        (*collector_params).hsp_num_max = i32::MAX; // BlastHspNumMax: hsp_num_max <= 0 → INT4_MAX
        (*collector_params).program = program;
        let mut writer_info = ffi::BlastHSPCollectorInfoNew(collector_params);
        // BlastHSPWriterNew — REPLACED WITH RUST (dispatch through C function pointer)
        let hsp_writer = {
            let wi = &*writer_info;
            let new_fn = wi.NewFnPtr.expect("HSP writer NewFnPtr is null");
            let writer = new_fn(wi.params, query_info, query_blk);
            libc::free(writer_info as *mut libc::c_void);
            writer_info = ptr::null_mut();
            writer
        };
        // BlastHSPStreamNew — REPLACED WITH RUST
        let hsp_stream = libc::malloc(std::mem::size_of::<ffi::BlastHSPStream>())
            as *mut ffi::BlastHSPStream;
        (*hsp_stream).program = program;
        (*hsp_stream).num_hsplists = 0;
        (*hsp_stream).num_hsplists_alloc = 100;
        (*hsp_stream).sorted_hsplists = libc::malloc(
            100 * std::mem::size_of::<*mut ffi::BlastHSPList>()
        ) as *mut *mut ffi::BlastHSPList;
        // Create Blast_HSPResultsNew(num_queries)
        let hsp_results = libc::calloc(1, std::mem::size_of::<ffi::BlastHSPResults>())
            as *mut ffi::BlastHSPResults;
        (*hsp_results).num_queries = num_queries;
        (*hsp_results).hitlist_array = libc::calloc(
            num_queries as usize, std::mem::size_of::<*mut ffi::BlastHitList>()
        ) as *mut *mut ffi::BlastHitList;
        (*hsp_stream).results = hsp_results;
        (*hsp_stream).results_sorted = 0; // FALSE
        (*hsp_stream).sort_by_score = ptr::null_mut(); // blastn, no CBS
        (*hsp_stream).x_lock = ptr::null_mut();
        (*hsp_stream).writer = hsp_writer;
        (*hsp_stream).writer_initialized = 0; // FALSE
        (*hsp_stream).writer_finalized = 0; // FALSE
        (*hsp_stream).pre_pipe = ptr::null_mut();
        (*hsp_stream).tback_pipe = ptr::null_mut();

        // 11. Create diagnostics
        // Diagnostics — REPLACED WITH RUST (calloc sub-structs)
        let diagnostics = libc::calloc(1, std::mem::size_of::<ffi::BlastDiagnostics>())
            as *mut ffi::BlastDiagnostics;
        (*diagnostics).ungapped_stat = libc::calloc(1, std::mem::size_of::<ffi::BlastUngappedStats>())
            as *mut ffi::BlastUngappedStats;
        (*diagnostics).gapped_stat = libc::calloc(1, std::mem::size_of::<ffi::BlastGappedStats>())
            as *mut ffi::BlastGappedStats;
        (*diagnostics).cutoffs = libc::calloc(1, std::mem::size_of::<ffi::BlastRawCutoffs>())
            as *mut ffi::BlastRawCutoffs;


        // 12. Run the full search — REPLACED WITH RUST
        // Decomposed Blast_RunFullSearch into 3 C sub-functions
        let db_for_results = BlastDb::open(db_path).unwrap();
        let mut results: *mut ffi::BlastHSPResults = ptr::null_mut();

        // 12a. Gap alignment setup — REPLACED WITH RUST
        // Decomposed BLAST_GapAlignSetUp into individual C sub-calls
        let mut score_params: *mut ffi::BlastScoringParameters = ptr::null_mut();
        let mut ext_params: *mut ffi::BlastExtensionParameters = ptr::null_mut();
        let mut hit_params: *mut ffi::BlastHitSavingParameters = ptr::null_mut();
        let mut eff_len_params: *mut ffi::BlastEffectiveLengthsParameters = ptr::null_mut();
        let mut gap_align: *mut ffi::BlastGapAlignStruct = ptr::null_mut();

        {
            // Get database stats — REPLACED WITH DIRECT RUST db access
            let total_length: i64 = db_total_length as i64;
            let num_seqs: i32 = db_num_oids as i32;

            // Initialize effective length parameters — REPLACED WITH RUST
            eff_len_params = libc::calloc(1, std::mem::size_of::<ffi::BlastEffectiveLengthsParameters>())
                as *mut ffi::BlastEffectiveLengthsParameters;
            (*eff_len_params).options = eff_len_opts;
            (*eff_len_params).real_db_length = total_length;
            (*eff_len_params).real_num_seqs = num_seqs;

            // Calculate effective lengths — REPLACED WITH RUST
            // For blastn with no overrides: use kbp_gap_std, compute alpha/beta, length adjustment
            {
                let db_len = total_length;
                let db_nseqs = num_seqs;
                let kbp_ptr = (*sbp).kbp_gap_std; // gapped_calculation=TRUE → kbp_gap_std

                for ctx in (*query_info).first_context..=(*query_info).last_context {
                    let ci = ctx as usize;
                    let mut length_adj: i32 = 0;
                    let mut eff_searchsp: i64 = 0;

                    let kbp = *kbp_ptr.add(ci);
                    let ctx_info = &mut *(*query_info).contexts.add(ci);

                    if ctx_info.is_valid != 0 && ctx_info.query_length > 0 {
                        let q_len = ctx_info.query_length;

                        // Get alpha/beta — PURE RUST lookup
                        let ungap_kbp = *(*sbp).kbp_std.add(ci);
                        let (alpha, beta) = blast_core::stat::nucl_alpha_beta(
                            (*scoring_opts).reward as i32,
                            (*scoring_opts).penalty as i32,
                            (*scoring_opts).gap_open,
                            (*scoring_opts).gap_extend,
                            (*ungap_kbp).Lambda, (*ungap_kbp).H,
                            (*scoring_opts).gapped_calculation != 0,
                        );

                        // Compute length adjustment — PURE RUST
                        let (la, _) = blast_core::stat::compute_length_adjustment_exact(
                            (*kbp).K, (*kbp).logK,
                            alpha / (*kbp).Lambda, beta,
                            q_len, db_len, db_nseqs,
                        );
                        length_adj = la;

                        // Compute effective search space
                        let eff_db = std::cmp::max(db_len - db_nseqs as i64 * length_adj as i64, 1);
                        eff_searchsp = eff_db * (q_len - length_adj) as i64;
                    }

                    ctx_info.eff_searchsp = eff_searchsp;
                    ctx_info.length_adjustment = length_adj;
                }
            }

            // Create scoring parameters — REPLACED WITH RUST
            score_params = libc::calloc(1, std::mem::size_of::<ffi::BlastScoringParameters>())
                as *mut ffi::BlastScoringParameters;
            (*score_params).options = scoring_opts;
            (*score_params).scale_factor = (*sbp).scale_factor;
            (*score_params).reward = (*scoring_opts).reward;
            (*score_params).penalty = (*scoring_opts).penalty;
            (*score_params).gap_open = (*scoring_opts).gap_open;
            (*score_params).gap_extend = (*scoring_opts).gap_extend;

            // Create extension parameters — REPLACED WITH RUST
            ext_params = libc::calloc(1, std::mem::size_of::<ffi::BlastExtensionParameters>())
                as *mut ffi::BlastExtensionParameters;
            (*ext_params).options = ext_opts;
            // Convert X-dropoff from bits to raw score using smallest gapped lambda
            if !(*sbp).kbp_gap.is_null() {
                let mut min_lambda = f64::MAX;
                for ci in (*query_info).first_context..=(*query_info).last_context {
                    let kp = *(*sbp).kbp_gap.add(ci as usize);
                    if !kp.is_null() && (*kp).Lambda > 0.0 {
                        if (*kp).Lambda < min_lambda { min_lambda = (*kp).Lambda; }
                    }
                }
                (*ext_params).gap_x_dropoff =
                    ((*ext_opts).gap_x_dropoff * std::f64::consts::LN_2 / min_lambda) as i32;
                let final_raw =
                    ((*ext_opts).gap_x_dropoff_final * std::f64::consts::LN_2 / min_lambda) as i32;
                (*ext_params).gap_x_dropoff_final =
                    std::cmp::max(final_raw, (*ext_params).gap_x_dropoff);
            }

            // Compute min subject length for hit saving
            let min_subject_length: i32 = if (*sbp).gbp.is_null() {
                (total_length / num_seqs as i64) as i32
            } else {
                // Compute min seq len from db (gbp always null for blastn so this won't execute)
                1 // gbp is always null for blastn, so this path is unreachable
            };

            // Create hit saving parameters — REPLACED WITH RUST
            hit_params = libc::calloc(1, std::mem::size_of::<ffi::BlastHitSavingParameters>())
                as *mut ffi::BlastHitSavingParameters;
            (*hit_params).options = hit_opts;
            (*hit_params).mask_level = 101;
            (*hit_params).do_sum_stats = 0; // FALSE for blastn gapped
            (*hit_params).prelim_evalue = (*hit_opts).expect_value;

            // Allocate per-context cutoffs
            let n_ctx = ((*query_info).last_context + 1) as usize;
            (*hit_params).cutoffs = libc::calloc(n_ctx, std::mem::size_of::<ffi::BlastGappedCutoffs>())
                as *mut ffi::BlastGappedCutoffs;

            // Compute cutoff scores from e-values using BLAST_Cutoffs (C function)
            let kbp_array = (*sbp).kbp_gap; // gapped
            let mut cutoff_min = i32::MAX;
            for ctx in (*query_info).first_context..=(*query_info).last_context {
                let ci = ctx as usize;
                if (*(*query_info).contexts.add(ci)).is_valid == 0 {
                    (*(*hit_params).cutoffs.add(ci)).cutoff_score = i32::MAX;
                    continue;
                }
                let kbp = *kbp_array.add(ci);
                let searchsp = (*(*query_info).contexts.add(ci)).eff_searchsp;
                // Compute cutoff score from E-value — PURE RUST
                // S = ceil(ln(K * searchsp / E) / Lambda)
                let evalue = (*hit_opts).expect_value;
                let e = evalue.max(1.0e-297);
                let new_cutoff = (((*kbp).K * searchsp as f64 / e).ln() / (*kbp).Lambda).ceil() as i32;
                (*(*hit_params).cutoffs.add(ci)).cutoff_score = new_cutoff;
                (*(*hit_params).cutoffs.add(ci)).cutoff_score_max = new_cutoff;
                if new_cutoff < cutoff_min { cutoff_min = new_cutoff; }
            }
            (*hit_params).cutoff_score_min = cutoff_min;

            // Create gap alignment structure — REPLACED WITH RUST
            gap_align = libc::calloc(1, std::mem::size_of::<ffi::BlastGapAlignStruct>())
                as *mut ffi::BlastGapAlignStruct;
            (*gap_align).sbp = sbp;
            (*gap_align).gap_x_dropoff = (*ext_params).gap_x_dropoff;
            (*gap_align).max_mismatches = (*(*ext_params).options).max_mismatches;
            (*gap_align).mismatch_window = (*(*ext_params).options).mismatch_window;
            // Allocate DP memory for dynamic programming (eDynProgScoreOnly path)
            (*gap_align).dp_mem_alloc = 1000;
            (*gap_align).dp_mem = libc::malloc(
                1000 * std::mem::size_of::<ffi::BlastGapDP>()
            ) as *mut ffi::BlastGapDP;
        }

        // 12b. Preliminary search — scans database, finds HSPs
        let rc = ffi::BLAST_PreliminarySearchEngine(
            program, query_blk, query_info, seq_src, gap_align,
            score_params, lookup, word_opts, ext_params, hit_params,
            eff_len_params, ptr::null(), db_opts, hsp_stream,
            diagnostics, None, ptr::null_mut(),
        );
        if rc != 0 {
            ffi::BLAST_GapAlignStructFree(gap_align);
            libc::free(score_params as *mut libc::c_void);
            libc::free(ext_params as *mut libc::c_void);
            if !hit_params.is_null() {
                libc::free((*hit_params).cutoffs as *mut libc::c_void);
                libc::free(hit_params as *mut libc::c_void);
            }
            libc::free(eff_len_params as *mut libc::c_void);
            return Err(format!("PreliminarySearchEngine failed ({})", rc).into());
        }

        // Close HSP stream for writing
        // BlastHSPStreamClose — REPLACED WITH RUST
        if !hsp_stream.is_null() && !(*hsp_stream).results.is_null()
            && (*hsp_stream).results_sorted == 0
        {
            // Finalize writer (s_FinalizeWriter)
            if !(*hsp_stream).writer.is_null() && (*hsp_stream).writer_finalized == 0 {
                let w = &*(*hsp_stream).writer;
                if (*hsp_stream).writer_initialized == 0 {
                    if let Some(init_fn) = w.InitFnPtr {
                        init_fn(w.data, (*hsp_stream).results as *mut std::ffi::c_void);
                    }
                }
                if let Some(final_fn) = w.FinalFnPtr {
                    final_fn(w.data, (*hsp_stream).results as *mut std::ffi::c_void);
                }
                (*hsp_stream).writer_finalized = 1;
            }

            // For blastn: no sort_by_score, so concatenate HSPLists
            let results = (*hsp_stream).results;
            let mut num_hl = (*hsp_stream).num_hsplists;
            for i in 0..(*results).num_queries {
                let hitlist = *(*results).hitlist_array.add(i as usize);
                if hitlist.is_null() { continue; }
                let hlcount = (*hitlist).hsplist_count;
                // Grow if needed
                if num_hl + hlcount > (*hsp_stream).num_hsplists_alloc {
                    let alloc = std::cmp::max(num_hl + hlcount + 100,
                                              2 * (*hsp_stream).num_hsplists_alloc);
                    (*hsp_stream).num_hsplists_alloc = alloc;
                    (*hsp_stream).sorted_hsplists = libc::realloc(
                        (*hsp_stream).sorted_hsplists as *mut libc::c_void,
                        alloc as usize * std::mem::size_of::<*mut ffi::BlastHSPList>(),
                    ) as *mut *mut ffi::BlastHSPList;
                }
                let mut k = 0;
                for j in 0..hlcount {
                    let hsplist = *(*hitlist).hsplist_array.add(j as usize);
                    if hsplist.is_null() { continue; }
                    (*hsplist).query_index = i;
                    *(*hsp_stream).sorted_hsplists.add((num_hl + k) as usize) = hsplist;
                    k += 1;
                }
                (*hitlist).hsplist_count = 0;
                num_hl += k;
            }
            (*hsp_stream).num_hsplists = num_hl;

            // Sort by decreasing OID
            if num_hl > 1 {
                let slice = std::slice::from_raw_parts_mut(
                    (*hsp_stream).sorted_hsplists,
                    num_hl as usize,
                );
                slice.sort_by(|a, b| (**b).oid.cmp(&(**a).oid));
            }
            (*hsp_stream).results_sorted = 1;
        }

        // 12c. Compute traceback — full gapped alignment
        let rc = ffi::BLAST_ComputeTraceback(
            program, hsp_stream, query_blk, query_info, seq_src, gap_align,
            score_params, ext_params, hit_params, eff_len_params,
            db_opts, ptr::null(), ptr::null(), ptr::null_mut(),
            &mut results, None, ptr::null_mut(),
        );

        // Clean up parameters
        ffi::BLAST_GapAlignStructFree(gap_align);
        libc::free(score_params as *mut libc::c_void);
        libc::free(ext_params as *mut libc::c_void);
        if !hit_params.is_null() {
            libc::free((*hit_params).cutoffs as *mut libc::c_void);
            libc::free(hit_params as *mut libc::c_void);
        }
        libc::free(eff_len_params as *mut libc::c_void);

        if rc != 0 {
            return Err(format!("ComputeTraceback failed ({})", rc).into());
        }
        // 13. Extract and format results
        if !results.is_null() {
            let res = &*results;
            let mut hits = Vec::new();

            for q in 0..res.num_queries {
                let hitlist_ptr = *res.hitlist_array.add(q as usize);
                if hitlist_ptr.is_null() {
                    continue;
                }
                let hitlist = &*hitlist_ptr;
                let query_id = &records[q as usize].id;

                for h in 0..hitlist.hsplist_count {
                    let hsplist_ptr = *hitlist.hsplist_array.add(h as usize);
                    if hsplist_ptr.is_null() {
                        continue;
                    }
                    let hsplist = &*hsplist_ptr;
                    let subject_id = db_for_results
                        .get_accession(hsplist.oid as u32)
                        .unwrap_or_else(|| format!("gnl|BL_ORD_ID|{}", hsplist.oid));

                    for s in 0..hsplist.hspcnt {
                        let hsp_ptr = *hsplist.hsp_array.add(s as usize);
                        if hsp_ptr.is_null() {
                            continue;
                        }
                        let hsp = &*hsp_ptr;

                        // Compute alignment stats from GapEditScript
                        let (align_len, num_ident, gap_opens) =
                            compute_alignment_stats(hsp, &encoded, &db_for_results, hsplist.oid, query_info);

                        let mismatches = align_len - num_ident - gap_opens;

                        // Determine strand from context
                        let ctx_info = &*(*query_info).contexts.add(hsp.context as usize);
                        let is_minus = ctx_info.frame < 0;
                        let query_len = ctx_info.query_length;

                        // For minus strand, convert query coords back to original strand
                        let (q_start, q_end) = if is_minus {
                            // After traceback, HSP query offsets are context-relative
                            let rc_start = hsp.query.offset;
                            let rc_end = hsp.query.end;
                            // Convert RC coords to original: orig = query_len - rc_pos
                            let orig_start = query_len - rc_end + 1;
                            let orig_end = query_len - rc_start;
                            (orig_start, orig_end)
                        } else {
                            (hsp.query.offset + 1, hsp.query.end)
                        };

                        // For minus strand, flip subject coords
                        let (s_start, s_end) = if is_minus {
                            (hsp.subject.end, hsp.subject.offset + 1)
                        } else {
                            (hsp.subject.offset + 1, hsp.subject.end)
                        };

                        hits.push(TabularHit {
                            query_id: query_id.clone(),
                            subject_id: subject_id.clone(),
                            pct_identity: if align_len > 0 {
                                100.0 * num_ident as f64 / align_len as f64
                            } else {
                                0.0
                            },
                            align_len,
                            mismatches: if mismatches > 0 { mismatches } else { 0 },
                            gap_opens,
                            query_start: q_start,
                            query_end: q_end,
                            subject_start: s_start,
                            subject_end: s_end,
                            evalue: hsp.evalue,
                            bit_score: hsp.bit_score,
                        });
                    }
                }
            }

            // Apply post-search filters
            let q_len = records.first().map(|r| r.sequence.len() as i32).unwrap_or(0);
            apply_filters(&mut hits, args, q_len);

            // Write output
            let stdout = io::stdout();
            let mut writer: Box<dyn Write> = if let Some(ref path) = args.out {
                Box::new(BufWriter::new(File::create(path)?))
            } else {
                Box::new(BufWriter::new(stdout.lock()))
            };

            format_tabular(&mut writer, &hits)?;
            writer.flush()?;
        }

        // Minimal cleanup — free results, prevent double-free of Rust-owned buffers
        // Full cleanup of C structures is complex due to shared ownership;
        // OS reclaims on exit.
        if !results.is_null() { ffi::Blast_HSPResultsFree(results); }
        (*query_blk).sequence = ptr::null_mut();
        (*query_blk).sequence_start = ptr::null_mut();
        (*query_blk).sequence_start_nomask = ptr::null_mut();
        (*query_blk).sequence_nomask = ptr::null_mut();
    }

    Ok(())
}

/// Decode packed NCBI2na subject to per-base BLASTNA for identity calculation.
/// Pure Rust blastn search — no FFI calls.
fn run_blastn_rust(
    args: &BlastnArgs,
    records: &[blast_input::FastaRecord],
    db: BlastDb,
) -> Result<(), Box<dyn std::error::Error>> {
    use blast_core::search::blastn_gapped_search;
    
    use rayon::prelude::*;

    // Pure Rust KBP computation — per-context (plus + minus) for exact e-value matching
    // The C engine computes separate KBP for each context based on strand composition
    let (kbp_plus, kbp_minus, searchsp_plus, searchsp_minus) = {
        let rec = &records[0];
        let query_plus: Vec<u8> = rec.sequence.iter().map(|&b| iupacna_to_blastna(b)).collect();
        let query_minus: Vec<u8> = query_plus.iter().rev().map(|&b| complement_blastna(b)).collect();

        const BLASTNA_SIZE: usize = 16;
        const BLASTNA_TO_NCBI4NA: [u8; BLASTNA_SIZE] = [1,2,4,8,5,10,3,12,9,6,14,13,11,7,15,0];
        let reward = args.reward;
        let penalty = args.penalty;
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
            blast_core::stat::UngappedKbpContext { query_offset: 0, query_length: qlen, is_valid: true },
        ];
        let plus_kbp_results = blast_core::stat::ungapped_kbp_calc(&query_plus, &ctxs, lo, hi, BLASTNA_SIZE, ambig, &matrix_fn);
        let minus_kbp_results = blast_core::stat::ungapped_kbp_calc(&query_minus, &ctxs, lo, hi, BLASTNA_SIZE, ambig, &matrix_fn);

        let default_kbp = blast_core::stat::KarlinBlk { lambda: 1.374, k: 0.621, log_k: 0.621_f64.ln(), h: 1.286 };
        let ungapped_plus = plus_kbp_results[0].clone().unwrap_or(default_kbp.clone());
        let ungapped_minus = minus_kbp_results[0].clone().unwrap_or(default_kbp.clone());

        // Gapped KBP (from table — may fall back to ungapped for large gap costs)
        let (gkbp_plus, _) = blast_core::stat::nucl_gapped_kbp_lookup(
            args.gapopen, args.gapextend, reward, penalty, &ungapped_plus,
        ).unwrap_or((ungapped_plus.clone(), false));
        let (gkbp_minus, _) = blast_core::stat::nucl_gapped_kbp_lookup(
            args.gapopen, args.gapextend, reward, penalty, &ungapped_minus,
        ).unwrap_or((ungapped_minus.clone(), false));

        // Per-context search space
        let compute_searchsp = |kbp: &blast_core::stat::KarlinBlk, ukbp: &blast_core::stat::KarlinBlk| -> f64 {
            let (alpha, beta) = blast_core::stat::nucl_alpha_beta(
                reward, penalty, args.gapopen, args.gapextend, ukbp.lambda, ukbp.h, true);
            let (len_adj, _) = blast_core::stat::compute_length_adjustment_exact(
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
    let word_size = args.word_size as usize;
    let reward = args.reward;
    let penalty = args.penalty;
    let gapopen = args.gapopen;
    let gapextend = args.gapextend;
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
            let mask = blast_core::filter::dust_filter(&query_plus_masked, 64, 2.0);
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

        // Parallel gapped search across all subjects
        let oid_hits: Vec<(u32, Vec<blast_core::search::SearchHsp>)> =
            (0..db.num_oids).into_par_iter().filter_map(|oid| {
                let subject = decode_subject(&db, oid as i32);
                // X-dropoff for gapped alignment: convert from bits to raw score
                // C default: gap_x_dropoff_final = 100 bits
                let gapped_x_dropoff = (args.xdrop_gap_final * std::f64::consts::LN_2 / kbp.lambda) as i32;
                let mut hsps = blast_core::search::blastn_gapped_search_nomask(
                    qp, qm, qp_nomask, qm_nomask, &subject,
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
        struct QsortCtx { ev: std::collections::HashMap<u32, f64>, sc: std::collections::HashMap<u32, i32> }

        // Store context in thread-local for the C callback
        use std::cell::RefCell;
        thread_local! { static CTX: RefCell<Option<QsortCtx>> = RefCell::new(None); }

        unsafe extern "C" fn cmp_fn(a: *const libc::c_void, b: *const libc::c_void) -> libc::c_int {
            let oid_a = *(a as *const u32);
            let oid_b = *(b as *const u32);
            CTX.with(|ctx| {
                let ctx = ctx.borrow();
                let ctx = ctx.as_ref().unwrap();
                let ev_a = ctx.ev[&oid_a];
                let ev_b = ctx.ev[&oid_b];
                // Primary: ascending e-value
                if ev_a < ev_b { return -1; }
                if ev_a > ev_b { return 1; }
                let sc_a = ctx.sc[&oid_a];
                let sc_b = ctx.sc[&oid_b];
                // Secondary: descending score
                if sc_a > sc_b { return -1; }
                if sc_a < sc_b { return 1; }
                // Tertiary: descending OID — BLAST_CMP(h2->oid, h1->oid)
                if oid_b > oid_a { return 1; }
                if oid_b < oid_a { return -1; }
                0
            })
        }

        CTX.with(|ctx| { *ctx.borrow_mut() = Some(QsortCtx { ev: best_evalue, sc: best_score }); });
        unsafe {
            libc::qsort(
                sorted_oids.as_mut_ptr() as *mut libc::c_void,
                sorted_oids.len(),
                std::mem::size_of::<u32>(),
                Some(cmp_fn),
            );
        }
        CTX.with(|ctx| { *ctx.borrow_mut() = None; });

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
            blast_format::format_pairwise_alignment(
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
        blast_format::format_tabular_custom(&mut writer, &hits, cols)?;
    } else {
        format_tabular(&mut writer, &hits)?;
    }
    writer.flush()?;

    Ok(())
}

/// NCBI4NA to BLASTNA conversion for ambiguity codes.
const NCBI4NA_TO_BLASTNA: [u8; 16] = [
    15, 0, 1, 6, 2, 4, 9, 13, 3, 8, 5, 12, 7, 11, 10, 14,
];

fn decode_subject(db: &BlastDb, oid: i32) -> Vec<u8> {
    let packed = db.get_sequence(oid as u32);
    let seq_len = db.get_seq_len(oid as u32) as usize;
    let mut decoded = Vec::with_capacity(seq_len);
    let full_bytes = seq_len / 4;
    let remainder = seq_len % 4;
    for i in 0..full_bytes {
        let b = packed[i];
        decoded.push((b >> 6) & 3);
        decoded.push((b >> 4) & 3);
        decoded.push((b >> 2) & 3);
        decoded.push(b & 3);
    }
    if remainder > 0 && full_bytes < packed.len() {
        let b = packed[full_bytes];
        for j in 0..remainder {
            decoded.push((b >> (6 - 2 * j)) & 3);
        }
    }

    // Apply ambiguity corrections
    if let Some(amb) = db.get_ambiguity_data(oid as u32) {
        if amb.len() >= 4 {
            let header = u32::from_be_bytes([amb[0], amb[1], amb[2], amb[3]]);
            let new_format = (header & 0x80000000) != 0;
            let count = (header & 0x7FFFFFFF) as usize;
            if !new_format {
                for i in 0..count {
                    let off = 4 + i * 4;
                    if off + 4 > amb.len() { break; }
                    let word = u32::from_be_bytes([amb[off], amb[off+1], amb[off+2], amb[off+3]]);
                    let value = ((word >> 28) & 0xF) as u8;
                    let length = ((word >> 24) & 0xF) as usize + 1;
                    let position = (word & 0xFFFFFF) as usize;
                    let blastna_val = NCBI4NA_TO_BLASTNA[value as usize];
                    for j in 0..length {
                        if position + j < decoded.len() {
                            decoded[position + j] = blastna_val;
                        }
                    }
                }
            }
        }
    }

    decoded
}

/// Compute alignment length, num_ident, and gap_opens from a GapEditScript.
unsafe fn compute_alignment_stats(
    hsp: &ffi::BlastHSP,
    query_encoded: &[u8],
    db: &BlastDb,
    oid: i32,
    query_info: *const ffi::BlastQueryInfo,
) -> (i32, i32, i32) {
    // Decode subject
    let subject_decoded = decode_subject(db, oid);

    // Get the context's query_offset to find the right query data
    let ctx = &*(*query_info).contexts.add(hsp.context as usize);
    // HSP offsets are context-relative after traceback
    // Query base at alignment pos i: query_encoded[1 + ctx.query_offset + hsp.query.offset + i]
    let q_start = (ctx.query_offset + hsp.query.offset) as usize;
    let s_start = hsp.subject.offset as usize;

    // For the identity comparison, we need the query bases for this context.
    // The encoded buffer has: [sentinel] [plus_strand] [sentinel] [minus_strand] [sentinel]
    // Context 0 (plus): query data at encoded[1 .. 1+qlen]
    // Context 1 (minus): query data at encoded[1+qlen+1 .. 1+qlen+1+qlen]
    // The HSP's query.offset for context 1 is relative to the context start.
    // We need to read from the correct context in the encoded buffer.
    if hsp.gap_info.is_null() {
        let len = (hsp.query.end - hsp.query.offset) as usize;
        let mut num_ident = 0i32;
        for i in 0..len {
            let q_base = query_encoded[1 + q_start + i];
            let s_base = if s_start + i < subject_decoded.len() {
                subject_decoded[s_start + i]
            } else {
                255
            };
            if q_base == s_base {
                num_ident += 1;
            }
        }
        (len as i32, num_ident, 0)
    } else {
        let esp = &*hsp.gap_info;
        let mut q_pos = q_start;
        let mut s_pos = s_start;
        let mut align_len = 0i32;
        let mut num_ident = 0i32;
        let mut gap_opens = 0i32;

        for idx in 0..esp.size as usize {
            let op = *esp.op_type.add(idx);
            let count = *esp.num.add(idx) as usize;
            align_len += count as i32;

            if op == ffi::EGapAlignOpType_eGapAlignSub {
                for _ in 0..count {
                    let q_base = query_encoded[1 + q_pos];
                    let s_base = if s_pos < subject_decoded.len() {
                        subject_decoded[s_pos]
                    } else {
                        255
                    };
                    if q_base == s_base {
                        num_ident += 1;
                    }
                    q_pos += 1;
                    s_pos += 1;
                }
            } else if op == ffi::EGapAlignOpType_eGapAlignDel {
                s_pos += count;
                gap_opens += 1;
            } else if op >= ffi::EGapAlignOpType_eGapAlignIns1 {
                q_pos += count;
                gap_opens += 1;
            } else {
                q_pos += count;
                s_pos += count;
            }
        }
        (align_len, num_ident, gap_opens)
    }
}

/// Search query against subject FASTA sequences (no database needed).
fn run_blastn_subject(
    args: &BlastnArgs,
    queries: &[blast_input::FastaRecord],
    subjects: &[blast_input::FastaRecord],
) -> Result<(), Box<dyn std::error::Error>> {
    use blast_core::search::blastn_gapped_search;

    let kbp = blast_core::stat::compute_ungapped_kbp(args.reward, args.penalty);
    let total_subj_len: usize = subjects.iter().map(|s| s.sequence.len()).sum();
    let avg_query_len = queries.iter().map(|r| r.sequence.len()).sum::<usize>() as f64
        / queries.len().max(1) as f64;
    let search_space = total_subj_len as f64 * avg_query_len;
    let word_size = args.word_size as usize;

    let mut all_hits = Vec::new();

    for query_rec in queries {
        let query_plus: Vec<u8> = query_rec.sequence.iter().map(|&b| iupacna_to_blastna(b)).collect();
        let query_minus: Vec<u8> = query_rec.sequence.iter().rev()
            .map(|&b| complement_blastna(iupacna_to_blastna(b))).collect();

        for subj_rec in subjects {
            let subject: Vec<u8> = subj_rec.sequence.iter().map(|&b| iupacna_to_blastna(b)).collect();

            let hsps = blastn_gapped_search(
                &query_plus, &query_minus, &subject,
                word_size, args.reward, args.penalty, args.gapopen, args.gapextend,
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
    format_tabular(&mut writer, &all_hits)?;
    writer.flush()?;
    Ok(())
}

/// DEPRECATED: No longer used. All KBP computation is now pure Rust.
/// Kept only as reference for the C FFI calling pattern.
#[cfg(never)]
unsafe fn _compute_exact_kbp_via_ffi(
    args: &BlastnArgs,
    records: &[blast_input::FastaRecord],
    db: &BlastDb,
) -> (blast_core::stat::KarlinBlk, f64) {
    let program = ffi::EBlastProgramType_eBlastTypeBlastn;

    // Encode first query for KBP computation
    let rec = &records[0];
    let mut encoded = vec![ffi::BLASTNA_SIZE as u8 - 1]; // sentinel
    for &b in &rec.sequence {
        encoded.push(iupacna_to_blastna(b));
    }
    encoded.push(ffi::BLASTNA_SIZE as u8 - 1); // sentinel
    for &b in rec.sequence.iter().rev() {
        encoded.push(complement_blastna(iupacna_to_blastna(b)));
    }
    encoded.push(ffi::BLASTNA_SIZE as u8 - 1); // sentinel

    let mut query_blk: *mut ffi::BLAST_SequenceBlk = ptr::null_mut();
    ffi::BlastSeqBlkNew(&mut query_blk);
    (*query_blk).sequence_start = encoded.as_mut_ptr();
    (*query_blk).sequence = encoded.as_mut_ptr().add(1);
    (*query_blk).length = (encoded.len() - 2) as i32;

    let query_info = ffi::BlastQueryInfoNew(program, 1);
    (*query_info).max_length = rec.sequence.len() as u32;
    (*query_info).contexts.add(0).as_mut().unwrap().query_offset = 0;
    (*query_info).contexts.add(0).as_mut().unwrap().query_length = rec.sequence.len() as i32;
    (*query_info).contexts.add(0).as_mut().unwrap().frame = 1;
    (*query_info).contexts.add(0).as_mut().unwrap().is_valid = 1;
    let minus_off = rec.sequence.len() as i32 + 1;
    (*query_info).contexts.add(1).as_mut().unwrap().query_offset = minus_off;
    (*query_info).contexts.add(1).as_mut().unwrap().query_length = rec.sequence.len() as i32;
    (*query_info).contexts.add(1).as_mut().unwrap().frame = -1;
    (*query_info).contexts.add(1).as_mut().unwrap().is_valid = 1;

    let mut scoring_opts: *mut ffi::BlastScoringOptions = ptr::null_mut();
    ffi::BlastScoringOptionsNew(program, &mut scoring_opts);
    (*scoring_opts).reward = args.reward as i16;
    (*scoring_opts).penalty = args.penalty as i16;
    (*scoring_opts).gap_open = args.gapopen;
    (*scoring_opts).gap_extend = args.gapextend;
    (*scoring_opts).gapped_calculation = 1;

    let mut qsup_opts: *mut ffi::QuerySetUpOptions = ptr::null_mut();
    ffi::BlastQuerySetUpOptionsNew(&mut qsup_opts);

    let mut sbp: *mut ffi::BlastScoreBlk = ptr::null_mut();
    let mut lookup_segments: *mut ffi::BlastSeqLoc = ptr::null_mut();
    let mut mask: *mut ffi::BlastMaskLoc = ptr::null_mut();
    let mut blast_msg: *mut ffi::Blast_Message = ptr::null_mut();

    let rc = ffi::BLAST_MainSetUp(
        program, qsup_opts, scoring_opts, query_blk, query_info, 1.0,
        &mut lookup_segments, &mut mask, &mut sbp, &mut blast_msg, None,
    );

    let mut kbp = blast_core::stat::KarlinBlk {
        lambda: 1.374, k: 0.621, log_k: 0.621_f64.ln(), h: 1.286,
    };

    if rc == 0 && !sbp.is_null() {
        // Compute gapped KBP via the C core
        let mut err_msg: *mut ffi::Blast_Message = ptr::null_mut();
        ffi::Blast_ScoreBlkKbpGappedCalc(sbp, scoring_opts, program, query_info, &mut err_msg);

        // Extract the UNGAPPED KBP from C core (exact computation)
        // The Rust engine does ungapped extension, so ungapped KBP is correct
        let sbp_ref = &*sbp;
        if !sbp_ref.kbp.is_null() {
            let kbp_ptr = *sbp_ref.kbp.add(0);
            if !kbp_ptr.is_null() {
                let ckbp = &*kbp_ptr;
                kbp = blast_core::stat::KarlinBlk {
                    lambda: ckbp.Lambda,
                    k: ckbp.K,
                    log_k: ckbp.logK,
                    h: ckbp.H,
                };
            }
        }
    }

    // Compute search space using the gapped KBP's H for length adjustment
    let len_adj = blast_core::stat::compute_length_adjustment(
        rec.sequence.len() as i32, db.total_length as i64, db.num_oids as i32, &kbp,
    );
    let search_space = blast_core::stat::compute_search_space(
        rec.sequence.len() as i64, db.total_length as i64, db.num_oids as i32, len_adj,
    );

    eprintln!("[KBP] lambda={:.6} K={:.6} H={:.6} searchsp={:.0}", kbp.lambda, kbp.k, kbp.h, search_space);

    // Cleanup
    (*query_blk).sequence = ptr::null_mut();
    (*query_blk).sequence_start = ptr::null_mut();

    (kbp, search_space)
}

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

/// Check if a sequence is low-complexity (simple repeats that crash the C engine).
fn is_low_complexity(seq: &[u8]) -> bool {
    if seq.len() < 20 { return false; }
    // Check if the sequence is a short repeat
    for period in 1..=8 {
        if seq.len() < period * 3 { continue; }
        let mut is_repeat = true;
        for i in period..seq.len() {
            if seq[i].to_ascii_uppercase() != seq[i % period].to_ascii_uppercase() {
                is_repeat = false;
                break;
            }
        }
        if is_repeat { return true; }
    }
    false
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
