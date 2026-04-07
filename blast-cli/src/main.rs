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
#[command(name = "blastn", about = "Nucleotide-nucleotide BLAST")]
struct BlastnArgs {
    /// Query file in FASTA format
    #[arg(short, long)]
    query: PathBuf,

    /// BLAST database name (path without extension)
    #[arg(short, long)]
    db: PathBuf,

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
    #[arg(short, long = "word_size", default_value = "28")]
    word_size: i32,
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

    // 2. Open BLAST database
    let db = BlastDb::open(&args.db)?;
    if db.db_type != DbType::Nucleotide {
        return Err("blastn requires a nucleotide database".into());
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
            ctx.frame = if i % 2 == 0 { 1 } else { -1 };
            ctx.is_valid = 1;
        }
        // Compute max_length
        let max_len = context_offsets.iter().map(|&(_, l)| l as u32).max().unwrap_or(0);
        (*query_info).max_length = max_len;

        // 6. Create options manually (BLAST_InitDefaultOptions causes issues)
        let mut scoring_opts: *mut ffi::BlastScoringOptions = ptr::null_mut();
        ffi::BlastScoringOptionsNew(program, &mut scoring_opts);
        // BlastScoringOptionsNew already sets blastn defaults:
        // reward=1, penalty=-3, gap_open=5, gap_extend=2, gapped=TRUE

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

        // Get results from HSP stream BEFORE closing it
        // (skipping traceback for now — it segfaults with multiple subjects)
        let results = (*hsp_stream).results;
        (*hsp_stream).results = ptr::null_mut();
        ffi::BlastHSPStreamClose(hsp_stream);

        // Re-open DB for result extraction
        let db_for_results = BlastDb::open(&args.db).unwrap();
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
                    let subject_id = format!("gnl|BL_ORD_ID|{}", hsplist.oid);

                    for s in 0..hsplist.hspcnt {
                        let hsp_ptr = *hsplist.hsp_array.add(s as usize);
                        if hsp_ptr.is_null() {
                            continue;
                        }
                        let hsp = &*hsp_ptr;

                        // Compute alignment stats from GapEditScript
                        let (align_len, num_ident, gap_opens) =
                            compute_alignment_stats(hsp, &encoded, &db_for_results, hsplist.oid);

                        let mismatches = align_len - num_ident - gap_opens;

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
                            query_start: hsp.query.offset + 1,
                            query_end: hsp.query.end,
                            subject_start: hsp.subject.offset + 1,
                            subject_end: hsp.subject.end,
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

        // Cleanup — free the main results and let the rest be reclaimed on exit.
        // TODO: proper cleanup of all C structures (currently segfaults due to
        // ownership issues between manually-decomposed Blast_RunFullSearch components)
        ffi::Blast_HSPResultsFree(results);
        // Prevent double-free of Rust-owned query buffer
        (*query_blk).sequence = ptr::null_mut();
        (*query_blk).sequence_start = ptr::null_mut();
    }

    Ok(())
}

/// Decode packed NCBI2na subject to per-base BLASTNA for identity calculation.
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
    decoded
}

/// Compute alignment length, num_ident, and gap_opens from a GapEditScript.
unsafe fn compute_alignment_stats(
    hsp: &ffi::BlastHSP,
    query_encoded: &[u8],
    db: &BlastDb,
    oid: i32,
) -> (i32, i32, i32) {
    // Decode subject
    let subject_decoded = decode_subject(db, oid);

    // Query is in BLASTNA encoding; sequence starts at encoded[1]
    let q_start = hsp.query.offset as usize;
    let s_start = hsp.subject.offset as usize;

    if hsp.gap_info.is_null() {
        // Ungapped: simple comparison
        let len = (hsp.query.end - hsp.query.offset) as usize;
        let mut num_ident = 0i32;
        for i in 0..len {
            // query_encoded[0] is the leading sentinel, query data starts at [1]
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
        // Gapped: walk the GapEditScript
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
                // Deletion in query = gap in query, advance subject
                s_pos += count;
                gap_opens += 1;
            } else if op >= ffi::EGapAlignOpType_eGapAlignIns1 {
                // Insertion in query = gap in subject, advance query
                q_pos += count;
                gap_opens += 1;
            } else {
                // Other op types (Del1, Del2) - advance both
                q_pos += count;
                s_pos += count;
            }
        }
        (align_len, num_ident, gap_opens)
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
