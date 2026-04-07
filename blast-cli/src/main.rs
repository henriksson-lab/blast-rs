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
#[command(name = "blastn", version = "0.1.0", about = "Nucleotide-nucleotide BLAST (Rust implementation)", allow_negative_numbers = true)]
struct BlastnArgs {
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

    /// Output format (6 = tabular)
    #[arg(long = "outfmt", default_value = "6")]
    outfmt: i32,

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

    /// Use pure Rust search engine (no FFI)
    #[arg(long = "rust-engine", default_value = "false")]
    rust_engine: bool,
}

fn main() {
    let args = BlastnArgs::parse();

    if let Err(e) = run_blastn(&args) {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    }
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

    // Use pure Rust engine if requested, or fall back for safety on edge cases
    if args.rust_engine {
        return run_blastn_rust(args, &records, db);
    }

    // Check for low-complexity/repetitive queries that crash the C engine
    // The C core segfaults on queries like ACGTACGTACGT... due to excessive seed hits
    for rec in &records {
        if is_low_complexity(&rec.sequence) {
            eprintln!("Warning: low-complexity query detected, using Rust engine");
            return run_blastn_rust(args, &records, db);
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
        // 4. Create BLAST_SequenceBlk
        let mut query_blk: *mut ffi::BLAST_SequenceBlk = ptr::null_mut();
        let rc = ffi::BlastSeqBlkNew(&mut query_blk);
        assert!(rc == 0 && !query_blk.is_null(), "Failed to create sequence block");

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

        // 5. Create BlastQueryInfo
        let num_contexts = num_queries * 2; // plus + minus strand
        let query_info = ffi::BlastQueryInfoNew(program, num_queries);
        assert!(!query_info.is_null(), "Failed to create query info");

        // Fill in context information
        (*query_info).first_context = 0;
        (*query_info).last_context = num_contexts - 1;
        (*query_info).num_queries = num_queries;

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

        // 6. Create options manually (BLAST_InitDefaultOptions causes issues)
        let mut scoring_opts: *mut ffi::BlastScoringOptions = ptr::null_mut();
        ffi::BlastScoringOptionsNew(program, &mut scoring_opts);
        // Set scoring from command line args
        (*scoring_opts).reward = args.reward as i16;
        (*scoring_opts).penalty = args.penalty as i16;
        (*scoring_opts).gap_open = args.gapopen;
        (*scoring_opts).gap_extend = args.gapextend;

        let mut word_opts: *mut ffi::BlastInitialWordOptions = ptr::null_mut();
        ffi::BlastInitialWordOptionsNew(program, &mut word_opts);
        // BlastInitialWordOptionsNew already sets blastn defaults (window_size=0 = single-hit)

        let mut ext_opts: *mut ffi::BlastExtensionOptions = ptr::null_mut();
        ffi::BlastExtensionOptionsNew(program, &mut ext_opts, 1);
        (*ext_opts).gap_x_dropoff = 30.0;
        (*ext_opts).gap_x_dropoff_final = 100.0;

        let mut hit_opts: *mut ffi::BlastHitSavingOptions = ptr::null_mut();
        ffi::BlastHitSavingOptionsNew(program, &mut hit_opts, 1);
        (*hit_opts).expect_value = args.evalue;
        (*hit_opts).hitlist_size = args.max_target_seqs;

        let mut eff_len_opts: *mut ffi::BlastEffectiveLengthsOptions = ptr::null_mut();
        ffi::BlastEffectiveLengthsOptionsNew(&mut eff_len_opts);

        let mut db_opts: *mut ffi::BlastDatabaseOptions = ptr::null_mut();
        ffi::BlastDatabaseOptionsNew(&mut db_opts);

        let mut qsup_opts: *mut ffi::QuerySetUpOptions = ptr::null_mut();
        ffi::BlastQuerySetUpOptionsNew(&mut qsup_opts);

        // 7. Setup ScoreBlk via BLAST_MainSetUp

        let mut sbp: *mut ffi::BlastScoreBlk = ptr::null_mut();
        let mut lookup_segments: *mut ffi::BlastSeqLoc = ptr::null_mut();
        let mut mask: *mut ffi::BlastMaskLoc = ptr::null_mut();
        let mut blast_msg: *mut ffi::Blast_Message = ptr::null_mut();

        let rc = ffi::BLAST_MainSetUp(
            program,
            qsup_opts,
            scoring_opts,
            query_blk,
            query_info,
            1.0, // scale_factor
            &mut lookup_segments,
            &mut mask,
            &mut sbp,
            &mut blast_msg,
            None, // get_path
        );


        if rc != 0 {
            let msg = if !blast_msg.is_null() {
                let msg_str = (*blast_msg).message;
                if !msg_str.is_null() {
                    std::ffi::CStr::from_ptr(msg_str).to_string_lossy().to_string()
                } else {
                    "Unknown error".to_string()
                }
            } else {
                "Unknown error".to_string()
            };
            return Err(format!("BLAST_MainSetUp failed ({}): {}", rc, msg).into());
        }

        // 8. Create lookup table
        let mut lookup: *mut ffi::LookupTableWrap = ptr::null_mut();
        let mut lookup_opts: *mut ffi::LookupTableOptions = ptr::null_mut();
        ffi::LookupTableOptionsNew(program, &mut lookup_opts);
        // Use BLAST_FillLookupTableOptions for correct lut_type selection
        let is_megablast = if args.word_size >= 28 { 1u8 } else { 0u8 };
        ffi::BLAST_FillLookupTableOptions(
            lookup_opts, program, is_megablast, 0.0, args.word_size,
        );

        // 9. Create SeqSrc from our Rust database reader
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
        let collector_params = ffi::BlastHSPCollectorParamsNew(
            hit_opts,
            0, // compositionBasedStats
            1, // gapped_calculation = TRUE
        );
        let mut writer_info = ffi::BlastHSPCollectorInfoNew(collector_params);
        let hsp_writer = ffi::BlastHSPWriterNew(&mut writer_info, query_info, query_blk);
        let hsp_stream = ffi::BlastHSPStreamNew(
            program,
            ext_opts,
            0, // sort_on_read (FALSE)
            num_queries,
            hsp_writer,
        );

        // 11. Create diagnostics
        let diagnostics = ffi::Blast_DiagnosticsInit();


        // 12. Run the search manually (decomposed Blast_RunFullSearch)
        let mut score_params: *mut ffi::BlastScoringParameters = ptr::null_mut();
        let mut ext_params: *mut ffi::BlastExtensionParameters = ptr::null_mut();
        let mut hit_params: *mut ffi::BlastHitSavingParameters = ptr::null_mut();
        let mut eff_len_params: *mut ffi::BlastEffectiveLengthsParameters = ptr::null_mut();
        let mut gap_align: *mut ffi::BlastGapAlignStruct = ptr::null_mut();

        let rc = ffi::BLAST_GapAlignSetUp(
            program, seq_src, scoring_opts, eff_len_opts, ext_opts, hit_opts,
            query_info, sbp,
            &mut score_params, &mut ext_params, &mut hit_params,
            &mut eff_len_params, &mut gap_align,
        );

        if rc != 0 {
            return Err(format!("BLAST_GapAlignSetUp failed ({})", rc).into());
        }

        let rc = ffi::BLAST_PreliminarySearchEngine(
            program, query_blk, query_info, seq_src, gap_align, score_params,
            lookup, word_opts, ext_params, hit_params, eff_len_params,
            ptr::null(), db_opts, hsp_stream, diagnostics, None, ptr::null_mut(),
        );

        if rc != 0 {
            return Err(format!("PreliminarySearchEngine failed ({})", rc).into());
        }

        ffi::BlastHSPStreamClose(hsp_stream);

        // Re-open DB for result extraction
        let db_for_results = BlastDb::open(db_path).unwrap();

        let mut results: *mut ffi::BlastHSPResults = ptr::null_mut();
        let rc = ffi::BLAST_ComputeTraceback(
            program, hsp_stream, query_blk, query_info, seq_src,
            gap_align, score_params, ext_params, hit_params, eff_len_params,
            db_opts, ptr::null(), ptr::null(), ptr::null_mut(),
            &mut results, None, ptr::null_mut(),
        );
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

    let kbp = blast_core::stat::compute_ungapped_kbp(args.reward, args.penalty);
    let avg_query_len = records.iter().map(|r| r.sequence.len()).sum::<usize>() as f64
        / records.len().max(1) as f64;
    let search_space = (db.total_length as f64) * avg_query_len;
    let word_size = args.word_size as usize;
    let reward = args.reward;
    let penalty = args.penalty;
    let gapopen = args.gapopen;
    let gapextend = args.gapextend;
    let evalue = args.evalue;

    let mut all_hits = Vec::new();

    let apply_dust = args.dust != "no";

    for rec in records {
        let mut query_plus: Vec<u8> = rec.sequence.iter().map(|&b| iupacna_to_blastna(b)).collect();

        // Apply DUST masking: replace low-complexity regions with N (14)
        if apply_dust {
            let mask = blast_core::filter::dust_filter(&query_plus, 64, 10.0);
            mask.apply(&mut query_plus, 14);
        }

        let query_minus: Vec<u8> = query_plus.iter().rev()
            .map(|&b| complement_blastna(b)).collect();

        // Parallel gapped search across all subjects
        let oid_hits: Vec<(u32, Vec<blast_core::search::SearchHsp>)> =
            (0..db.num_oids).into_par_iter().filter_map(|oid| {
                let subject = decode_subject(&db, oid as i32);
                let hsps = blastn_gapped_search(
                    &query_plus, &query_minus, &subject,
                    word_size, reward, penalty, gapopen, gapextend,
                    20, &kbp, search_space, evalue,
                );
                if hsps.is_empty() { None } else { Some((oid, hsps)) }
            }).collect();

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

                all_hits.push(TabularHit {
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
                    evalue: hsp.evalue,
                    bit_score: hsp.bit_score,
                });
            }
        }
    }

    // Sort by e-value and limit to max_target_seqs unique subjects
    all_hits.sort_by(|a, b| a.evalue.partial_cmp(&b.evalue).unwrap_or(std::cmp::Ordering::Equal));
    let max_subjects = args.max_target_seqs as usize;
    let mut seen_subjects = std::collections::HashSet::new();
    all_hits.retain(|hit| {
        if seen_subjects.len() >= max_subjects && !seen_subjects.contains(&hit.subject_id) {
            false
        } else {
            seen_subjects.insert(hit.subject_id.clone());
            true
        }
    });

    // Output
    let stdout = io::stdout();
    let mut writer: Box<dyn Write> = if let Some(ref path) = args.out {
        Box::new(BufWriter::new(File::create(path)?))
    } else {
        Box::new(BufWriter::new(stdout.lock()))
    };
    if args.outfmt == 0 {
        // Pairwise text output
        for hit in &all_hits {
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
    } else {
        format_tabular(&mut writer, &all_hits)?;
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
