//! High-level public API for running BLAST searches.
//!
//! Provides `blastp`, `blastn`, `blastx`, `tblastn`, and other search functions
//! with a simple builder-pattern interface, matching the API of the previous
//! blast-rs implementation.

use std::path::Path;
use std::io::{self, Write, BufWriter};
use std::fs::File;

use crate::db::{BlastDb, DbType};
use crate::matrix::AA_SIZE;
use crate::encoding::{AMINOACID_TO_NCBISTDAA, NCBISTDAA_TO_AMINOACID, IUPACNA_TO_BLASTNA};
use crate::search::{SearchHsp, blastn_gapped_search_nomask};
use crate::stat::{KarlinBlk, ungapped_kbp_calc, nucl_gapped_kbp_lookup,
                   nucl_alpha_beta, compute_length_adjustment_exact, UngappedKbpContext};
use crate::traceback::build_blastna_matrix;

// ── Result types ────────────────────────────────────────────────────────────

/// A single high-scoring segment pair from a BLAST search.
#[derive(Debug, Clone)]
pub struct Hsp {
    pub score: i32,
    pub bit_score: f64,
    pub evalue: f64,
    pub query_start: usize,
    pub query_end: usize,
    pub subject_start: usize,
    pub subject_end: usize,
    pub num_identities: usize,
    pub num_gaps: usize,
    pub alignment_length: usize,
    pub query_aln: Vec<u8>,
    pub midline: Vec<u8>,
    pub subject_aln: Vec<u8>,
    pub query_frame: i32,
    pub subject_frame: i32,
}

impl Hsp {
    /// Percent identity of this HSP.
    pub fn percent_identity(&self) -> f64 {
        if self.alignment_length == 0 { return 0.0; }
        100.0 * self.num_identities as f64 / self.alignment_length as f64
    }
}

/// Result for one subject sequence (may contain multiple HSPs).
#[derive(Debug, Clone)]
pub struct SearchResult {
    pub subject_oid: u32,
    pub subject_title: String,
    pub subject_accession: String,
    pub subject_len: usize,
    pub hsps: Vec<Hsp>,
    pub taxids: Vec<i32>,
}

impl SearchResult {
    pub fn best_evalue(&self) -> f64 {
        self.hsps.iter().map(|h| h.evalue).fold(f64::INFINITY, f64::min)
    }
}

// ── Search parameters ───────────────────────────────────────────────────────

/// Scoring matrix type for protein searches.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MatrixType {
    Blosum45, Blosum50, Blosum62, Blosum80, Blosum90,
    Pam30, Pam70, Pam250,
}

/// Scoring matrix wrapper (28x28 NCBIstdaa-indexed scores with metadata).
#[derive(Debug, Clone)]
pub struct ScoringMatrix {
    pub scores: [[i32; AA_SIZE]; AA_SIZE],
    pub min_score: i32,
    pub name: String,
}

impl ScoringMatrix {
    pub fn from_type(mt: MatrixType) -> Self {
        let scores = *get_matrix(mt);
        let min_score = scores.iter()
            .flat_map(|row| row.iter())
            .filter(|&&v| v > i32::MIN / 4)
            .copied()
            .min()
            .unwrap_or(-4);
        let name = format!("{:?}", mt).to_uppercase();
        ScoringMatrix { scores, min_score, name }
    }
    pub fn blosum62() -> Self { Self::from_type(MatrixType::Blosum62) }
    pub fn blosum45() -> Self { Self::from_type(MatrixType::Blosum45) }
    pub fn blosum50() -> Self { Self::from_type(MatrixType::Blosum50) }
    pub fn blosum80() -> Self { Self::from_type(MatrixType::Blosum80) }
    pub fn blosum90() -> Self { Self::from_type(MatrixType::Blosum90) }
    pub fn pam30() -> Self { Self::from_type(MatrixType::Pam30) }
    pub fn pam70() -> Self { Self::from_type(MatrixType::Pam70) }
    pub fn pam250() -> Self { Self::from_type(MatrixType::Pam250) }
    pub fn score(&self, a: u8, b: u8) -> i32 {
        self.scores[a as usize & 0x1F][b as usize & 0x1F]
    }
}

/// Search configuration with builder pattern.
#[derive(Debug, Clone)]
pub struct SearchParams {
    pub word_size: usize,
    pub matrix: MatrixType,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub evalue_threshold: f64,
    pub max_target_seqs: usize,
    pub x_drop_ungapped: i32,
    pub x_drop_gapped: i32,
    pub ungapped_cutoff: i32,
    pub min_score: i32,
    pub match_score: i32,
    pub mismatch: i32,
    pub num_threads: usize,
    pub filter_low_complexity: bool,
    pub comp_adjust: u8,
    pub strand: String,
    pub query_gencode: u8,
    pub db_gencode: u8,
    pub max_hsps: Option<usize>,
    pub culling_limit: Option<usize>,
    pub two_hit: bool,
    pub two_hit_window: usize,
    pub x_drop_final: i32,
    pub soft_masking: bool,
    pub lcase_masking: bool,
}

impl SearchParams {
    pub fn blastp() -> Self { Self::blastp_defaults() }
    pub fn blastn() -> Self { Self::blastn_defaults() }

    pub fn blastx() -> Self { Self::blastp_defaults() }
    pub fn tblastn() -> Self { Self::blastp_defaults() }
    pub fn tblastx() -> Self { Self::blastp_defaults() }

    pub fn blastp_defaults() -> Self {
        SearchParams {
            word_size: 3,
            matrix: MatrixType::Blosum62,
            gap_open: 11,
            gap_extend: 1,
            evalue_threshold: 10.0,
            max_target_seqs: 500,
            x_drop_ungapped: 7,   // BLAST_UNGAPPED_X_DROPOFF_PROT (bits)
            x_drop_gapped: 15,    // BLAST_GAP_X_DROPOFF_PROT (in bits, converted to raw by blastp())
            ungapped_cutoff: 0,
            min_score: 0,
            match_score: 1,
            mismatch: -2,
            num_threads: 1,
            filter_low_complexity: true,
            comp_adjust: 2,       // NCBI default: conditional compositional matrix adjust
            strand: "both".to_string(),
            query_gencode: 1,
            db_gencode: 1,
            max_hsps: None,
            culling_limit: None,
            two_hit: false,
            two_hit_window: 40,
            x_drop_final: 25,     // BLAST_GAP_X_DROPOFF_FINAL_PROT
            soft_masking: false,
            lcase_masking: false,
        }
    }

    pub fn blastn_defaults() -> Self {
        SearchParams {
            word_size: 11,
            matrix: MatrixType::Blosum62,
            gap_open: 5,
            gap_extend: 2,
            evalue_threshold: 10.0,
            max_target_seqs: 500,
            x_drop_ungapped: 20,  // BLAST_UNGAPPED_X_DROPOFF_NUCL
            x_drop_gapped: 30,    // BLAST_GAP_X_DROPOFF_NUCL
            ungapped_cutoff: 0,
            min_score: 0,
            match_score: 2,
            mismatch: -3,
            num_threads: 1,
            filter_low_complexity: true,
            comp_adjust: 0,
            strand: "both".to_string(),
            query_gencode: 1,
            db_gencode: 1,
            max_hsps: None,
            culling_limit: None,
            two_hit: false,
            two_hit_window: 40,
            x_drop_final: 100,    // BLAST_GAP_X_DROPOFF_FINAL_NUCL
            soft_masking: false,
            lcase_masking: false,
        }
    }

    // Builder methods
    pub fn evalue(mut self, v: f64) -> Self { self.evalue_threshold = v; self }
    pub fn max_target_seqs(mut self, v: usize) -> Self { self.max_target_seqs = v; self }
    pub fn num_threads(mut self, v: usize) -> Self { self.num_threads = v; self }
    pub fn filter_low_complexity(mut self, v: bool) -> Self { self.filter_low_complexity = v; self }
    pub fn comp_adjust(mut self, v: u8) -> Self { self.comp_adjust = v; self }
    pub fn strand(mut self, v: &str) -> Self { self.strand = v.to_string(); self }
    pub fn word_size(mut self, v: usize) -> Self { self.word_size = v; self }
    pub fn matrix(mut self, v: MatrixType) -> Self { self.matrix = v; self }
    pub fn gap_open(mut self, v: i32) -> Self { self.gap_open = v; self }
    pub fn gap_extend(mut self, v: i32) -> Self { self.gap_extend = v; self }
    pub fn match_score(mut self, v: i32) -> Self { self.match_score = v; self }
    pub fn mismatch(mut self, v: i32) -> Self { self.mismatch = v; self }
    pub fn query_gencode(mut self, v: u8) -> Self { self.query_gencode = v; self }
    pub fn db_gencode(mut self, v: u8) -> Self { self.db_gencode = v; self }
    pub fn max_hsps(mut self, v: Option<usize>) -> Self { self.max_hsps = v; self }
    pub fn culling_limit(mut self, v: Option<usize>) -> Self { self.culling_limit = v; self }
    pub fn two_hit(mut self, v: bool) -> Self { self.two_hit = v; self }
    pub fn two_hit_window(mut self, v: usize) -> Self { self.two_hit_window = v; self }
    pub fn x_drop_final(mut self, v: i32) -> Self { self.x_drop_final = v; self }
    pub fn soft_masking(mut self, v: bool) -> Self { self.soft_masking = v; self }
    pub fn lcase_masking(mut self, v: bool) -> Self { self.lcase_masking = v; self }
}

// ── Database builder ────────────────────────────────────────────────────────

/// Entry for building a database.
#[derive(Debug, Clone)]
pub struct SequenceEntry {
    pub title: String,
    pub accession: String,
    pub sequence: Vec<u8>,
    pub taxid: Option<u32>,
}

/// Builder for creating BLAST databases.
pub struct BlastDbBuilder {
    pub seq_type: DbType,
    pub db_title: String,
    pub entries: Vec<SequenceEntry>,
}

impl BlastDbBuilder {
    pub fn new(seq_type: DbType, db_title: impl Into<String>) -> Self {
        BlastDbBuilder {
            seq_type,
            db_title: db_title.into(),
            entries: Vec::new(),
        }
    }

    pub fn add(&mut self, entry: SequenceEntry) {
        self.entries.push(entry);
    }

    pub fn write(&self, base_path: &Path) -> io::Result<()> {
        match self.seq_type {
            DbType::Nucleotide => self.write_nucleotide(base_path),
            DbType::Protein => self.write_protein(base_path),
        }
    }

    fn write_nucleotide(&self, base_path: &Path) -> io::Result<()> {
        // Write .nsq
        let mut nsq = BufWriter::new(File::create(base_path.with_extension("nsq"))?);
        nsq.write_all(&[0u8])?; // sentinel

        let mut seq_offsets = vec![1u32];
        let mut amb_offsets = Vec::new();

        for entry in &self.entries {
            let seq = &entry.sequence;
            let iupac_to_2na = |b: u8| -> u8 {
                match b {
                    b'A' | b'a' => 0, b'C' | b'c' => 1,
                    b'G' | b'g' => 2, b'T' | b't' => 3, _ => 0,
                }
            };

            let mut packed = Vec::new();
            let full_bytes = seq.len() / 4;
            let remainder = seq.len() % 4;
            for i in 0..full_bytes {
                let b = (iupac_to_2na(seq[i*4]) << 6)
                    | (iupac_to_2na(seq[i*4+1]) << 4)
                    | (iupac_to_2na(seq[i*4+2]) << 2)
                    | iupac_to_2na(seq[i*4+3]);
                packed.push(b);
            }
            if remainder > 0 {
                let mut last = 0u8;
                for j in 0..remainder {
                    last |= iupac_to_2na(seq[full_bytes * 4 + j]) << (6 - 2 * j);
                }
                last |= remainder as u8;
                packed.push(last);
            } else {
                packed.push(0);
            }

            nsq.write_all(&packed)?;
            let seq_start = *seq_offsets.last().unwrap();
            let amb_offset = seq_start + packed.len() as u32;
            amb_offsets.push(amb_offset);
            seq_offsets.push(amb_offset);
        }
        amb_offsets.push(*seq_offsets.last().unwrap_or(&0));
        nsq.flush()?;

        // Write .nhr
        let mut nhr = BufWriter::new(File::create(base_path.with_extension("nhr"))?);
        let mut hdr_offsets = vec![0u32];
        for (oid, entry) in self.entries.iter().enumerate() {
            let hdr = format!("{} {}", entry.accession, entry.title);
            let asn1 = encode_defline_asn1(&hdr, oid as i32);
            nhr.write_all(&asn1)?;
            hdr_offsets.push(hdr_offsets.last().unwrap() + asn1.len() as u32);
        }
        nhr.flush()?;

        // Write .nin
        write_index_file(
            &base_path.with_extension("nin"),
            4, // format version
            DbType::Nucleotide,
            &self.db_title,
            self.entries.len() as u32,
            &hdr_offsets,
            &seq_offsets,
            Some(&amb_offsets),
        )
    }

    fn write_protein(&self, base_path: &Path) -> io::Result<()> {
        // Write .psq (protein sequences in NCBIstdaa)
        let mut psq = BufWriter::new(File::create(base_path.with_extension("psq"))?);
        psq.write_all(&[0u8])?; // sentinel

        let mut seq_offsets = vec![1u32];
        for entry in &self.entries {
            let encoded: Vec<u8> = entry.sequence.iter()
                .map(|&b| AMINOACID_TO_NCBISTDAA[b as usize & 0x7F])
                .collect();
            psq.write_all(&encoded)?;
            psq.write_all(&[0u8])?; // sentinel between sequences
            let prev = *seq_offsets.last().unwrap();
            seq_offsets.push(prev + encoded.len() as u32 + 1);
        }
        psq.flush()?;

        // Write .phr
        let mut phr = BufWriter::new(File::create(base_path.with_extension("phr"))?);
        let mut hdr_offsets = vec![0u32];
        for (oid, entry) in self.entries.iter().enumerate() {
            let hdr = format!("{} {}", entry.accession, entry.title);
            let asn1 = encode_defline_asn1(&hdr, oid as i32);
            phr.write_all(&asn1)?;
            hdr_offsets.push(hdr_offsets.last().unwrap() + asn1.len() as u32);
        }
        phr.flush()?;

        // Write .pin
        write_index_file(
            &base_path.with_extension("pin"),
            4,
            DbType::Protein,
            &self.db_title,
            self.entries.len() as u32,
            &hdr_offsets,
            &seq_offsets,
            None,
        )
    }
}

fn encode_defline_asn1(header: &str, _oid: i32) -> Vec<u8> {
    // Minimal ASN.1 BER encoding of Blast-def-line-set
    let title_bytes = header.as_bytes();
    let inner_len = 2 + title_bytes.len();
    let mut out = Vec::with_capacity(inner_len + 10);
    // SEQUENCE { VisibleString title }
    out.push(0x30); // SEQUENCE tag
    if inner_len < 128 {
        out.push(inner_len as u8);
    } else {
        out.push(0x81);
        out.push(inner_len as u8);
    }
    out.push(0x1A); // VisibleString tag
    out.push(title_bytes.len() as u8);
    out.extend_from_slice(title_bytes);
    out
}

fn write_index_file(
    path: &Path,
    format_version: u32,
    db_type: DbType,
    title: &str,
    num_oids: u32,
    hdr_offsets: &[u32],
    seq_offsets: &[u32],
    amb_offsets: Option<&[u32]>,
) -> io::Result<()> {
    use byteorder::{BigEndian, WriteBytesExt};
    let mut f = BufWriter::new(File::create(path)?);

    f.write_u32::<BigEndian>(format_version)?;
    let db_type_val: u32 = match db_type { DbType::Protein => 1, DbType::Nucleotide => 0 };
    f.write_u32::<BigEndian>(db_type_val)?;

    // Title (length-prefixed)
    f.write_u32::<BigEndian>(title.len() as u32)?;
    f.write_all(title.as_bytes())?;

    // Timestamp placeholder
    let ts = "2024-01-01T00:00:00";
    f.write_u32::<BigEndian>(ts.len() as u32)?;
    f.write_all(ts.as_bytes())?;

    f.write_u32::<BigEndian>(num_oids)?;

    // Total residue count (not critical for search)
    let total: u64 = 0;
    f.write_u64::<BigEndian>(total)?;
    // Max seq length
    f.write_u32::<BigEndian>(0)?;

    // Write offset arrays
    for &off in hdr_offsets { f.write_u32::<BigEndian>(off)?; }
    for &off in seq_offsets { f.write_u32::<BigEndian>(off)?; }
    if let Some(amb) = amb_offsets {
        for &off in amb { f.write_u32::<BigEndian>(off)?; }
    }

    f.flush()
}

// ── Search functions ────────────────────────────────────────────────────────

/// Run a blastp search (protein query vs protein database).
pub fn blastp(db: &BlastDb, query: &[u8], params: &SearchParams) -> Vec<SearchResult> {
    if query.is_empty() { return Vec::new(); }

    let query_aa: Vec<u8> = query.iter()
        .map(|&b| AMINOACID_TO_NCBISTDAA[b as usize & 0x7F])
        .collect();

    let matrix = *get_matrix(params.matrix);
    let word_size = params.word_size.max(2).min(6);
    let threshold = 11.0;

    let prot_kbp = crate::stat::lookup_protein_params(params.gap_open, params.gap_extend)
        .map(|p| crate::stat::KarlinBlk { lambda: p.lambda, k: p.k, log_k: p.k.ln(), h: p.h })
        .unwrap_or_else(crate::stat::protein_ungapped_kbp);

    // Convert bit-score x_dropoff values to raw scores.
    // NCBI uses UNGAPPED KBP for ungapped x_drop, GAPPED KBP for gapped x_drop.
    let ln2 = std::f64::consts::LN_2;
    let ungapped_kbp = crate::stat::protein_ungapped_kbp();
    let x_drop_ungapped = (params.x_drop_ungapped as f64 * ln2 / ungapped_kbp.lambda).ceil() as i32;
    let x_drop_gapped = (params.x_drop_gapped as f64 * ln2 / prot_kbp.lambda).ceil() as i32;
    let x_drop_final = (params.x_drop_final as f64 * ln2 / prot_kbp.lambda).ceil() as i32;

    // Gap trigger: BLAST_GAP_TRIGGER_PROT = 22.0 bits → raw score.
    // NCBI uses UNGAPPED KBP (kbp_std): (Int4)((bits * LN2 + logK) / Lambda)
    // = (Int4)((22 * 0.693 + ln(0.134)) / 0.3176) = (Int4)(41.7) = 41
    let gap_trigger_raw = ((22.0 * ln2 + ungapped_kbp.log_k) / ungapped_kbp.lambda) as i32;

    let total_subj_len: usize = (0..db.num_oids)
        .map(|oid| db.get_seq_len(oid) as usize)
        .sum();

    // Use exact length adjustment with alpha/beta from gapped params (matching NCBI C engine)
    let gapped_params = crate::stat::lookup_protein_params(params.gap_open, params.gap_extend);
    let (len_adj, search_space) = if let Some(ref gp) = gapped_params {
        let alpha_d_lambda = gp.alpha / prot_kbp.lambda;
        let (adj, _) = crate::stat::compute_length_adjustment_exact(
            prot_kbp.k, prot_kbp.log_k,
            alpha_d_lambda, gp.beta,
            query_aa.len() as i32, total_subj_len as i64, db.num_oids as i32,
        );
        let eff_q = (query_aa.len() as i64 - adj as i64).max(1);
        let eff_db = (total_subj_len as i64 - db.num_oids as i64 * adj as i64).max(1);
        (adj, eff_q as f64 * eff_db as f64)
    } else {
        let adj = crate::stat::compute_length_adjustment(
            query_aa.len() as i32, total_subj_len as i64, db.num_oids as i32, &prot_kbp);
        let ss = crate::stat::compute_search_space(
            query_aa.len() as i64, total_subj_len as i64, db.num_oids as i32, adj);
        (adj, ss)
    };

    // Build Gumbel block for Spouge FSC e-value (per-subject length correction)
    let gumbel_blk = crate::stat::protein_gumbel_blk(
        params.gap_open, params.gap_extend, total_subj_len as i64);

    // Build lookup table once per query (not per subject).
    let lookup_table = crate::protein_lookup::ProteinLookupTable::build(
        &query_aa, word_size, &matrix, threshold,
    );

    // Configure threading
    let num_threads = if params.num_threads == 0 { rayon::current_num_threads() } else { params.num_threads };

    let max_hsps = params.max_hsps;
    let evalue_threshold = params.evalue_threshold;
    let gap_open = params.gap_open;
    let gap_extend = params.gap_extend;

    // Process a single subject OID — extracted so it can be called
    // either sequentially or in parallel without per-call pool overhead.
    let search_oid = |oid: u32| -> Option<SearchResult> {
        let subj_raw = db.get_sequence(oid);
        let subj_len = db.get_seq_len(oid) as usize;
        if subj_len < word_size { return None; }

        // Use length-based slice — no allocation (matches NCBI C approach).
        // get_sequence() includes trailing sentinel; subj_len excludes it.
        let subj_aa = &subj_raw[..subj_len];

        let ungapped_hits = crate::protein_lookup::protein_scan_with_table(
            &query_aa, &subj_aa, &matrix, &lookup_table, x_drop_ungapped,
        );
        if ungapped_hits.is_empty() { return None; }

        // Use the best ungapped hit as seed for gapped alignment if it passes
        // the gap trigger threshold. The e-value on the final gapped score is the real filter.
        // Use 90% of gap_trigger_raw to account for minor scoring differences
        // in ungapped extension vs NCBI C.
        let adjusted_cutoff = (gap_trigger_raw * 9) / 10;
        let best_seed = ungapped_hits.iter()
            .filter(|uh| uh.score >= adjusted_cutoff)
            .max_by_key(|uh| uh.score);

        let mut phits = Vec::new();
        if let Some(uh) = best_seed {
            let seed_q = (uh.query_start + uh.query_end) / 2;
            let seed_s = (uh.subject_start + uh.subject_end) / 2;
            if let Some(gr) = crate::protein::protein_gapped_align(
                &query_aa, &subj_aa, seed_q, seed_s,
                &matrix, gap_open, gap_extend, x_drop_final,
            ) {
                let q_slice = &query_aa[gr.query_start..gr.query_end];
                let s_slice = &subj_aa[gr.subject_start..gr.subject_end];
                let (qs, ss) = gr.edit_script.render_alignment(
                    q_slice, s_slice, crate::protein::ncbistdaa_to_char,
                );
                phits.push(crate::protein_lookup::ProteinHit {
                    query_start: gr.query_start, query_end: gr.query_end,
                    subject_start: gr.subject_start, subject_end: gr.subject_end,
                    score: gr.score, num_ident: gr.num_ident,
                    align_length: gr.align_length, mismatches: gr.mismatches,
                    gap_opens: gr.gap_opens, qseq: Some(qs), sseq: Some(ss),
                });
            }
        }
        if phits.is_empty() { phits = ungapped_hits; }
        if phits.is_empty() { return None; }

        // Pre-filter: check e-value with Spouge FSC if available, else simple Karlin.
        let best_raw_ev = phits.iter()
            .map(|ph| {
                if let Some(ref gbp) = gumbel_blk {
                    crate::stat::spouge_evalue(ph.score, &prot_kbp, gbp,
                        query_aa.len() as i32, subj_len as i32)
                } else {
                    prot_kbp.raw_to_evalue(ph.score, search_space)
                }
            })
            .fold(f64::MAX, f64::min);
        if best_raw_ev > evalue_threshold { return None; }

        let accession = db.get_accession(oid).unwrap_or_else(|| format!("oid_{}", oid));
        let title = String::from_utf8_lossy(db.get_header(oid)).to_string();
        let sl = subj_aa.len();

        // Composition-based e-value adjustment (NCBI comp_based_stats).
        // Verbatim port of Blast_AdjustScores + s_AdjustEvaluesForComposition.
        let comp_mode = params.comp_adjust;

        // Determine adjustment rule and optionally build adjusted matrix.
        // adj_result: None = no adjustment, Some((adjusted_matrix_opt, lambda_ratio_opt))
        let adj_result: Option<(Option<[[i32; AA_SIZE]; AA_SIZE]>, Option<f64>)> = if comp_mode > 0 {
            let (qcomp28, qn) = crate::composition::read_composition(&query_aa, AA_SIZE);
            let (scomp28, sn) = crate::composition::read_composition(subj_aa, AA_SIZE);
            if qn == 0 || sn == 0 {
                None
            } else {
                let mut qp20 = [0.0f64; 20];
                let mut sp20 = [0.0f64; 20];
                crate::compo_mode_condition::gather_letter_probs(&qcomp28, &mut qp20);
                crate::compo_mode_condition::gather_letter_probs(&scomp28, &mut sp20);

                let rule = crate::compo_mode_condition::choose_matrix_adjust_rule(
                    query_aa.len(), subj_aa.len(), &qp20, &sp20, comp_mode,
                );

                use crate::compo_mode_condition::MatrixAdjustRule;
                match rule {
                    MatrixAdjustRule::DontAdjust => None,
                    MatrixAdjustRule::ScaleOldMatrix => {
                        // Port of NCBI Blast_CompositionBasedStats: rescale matrix
                        // using composition-specific lambda ratio, then re-align.
                        let ungapped_lambda = crate::stat::protein_ungapped_kbp().lambda;
                        // Build the frequency ratio matrix from the standard BLOSUM62
                        // joint probs (used as freq_ratios in s_ScaleSquareMatrix).
                        // For non-position-based, startFreqRatios is initialized from
                        // the standard matrix frequency ratios.
                        let freq_ratios = crate::matrix::get_blosum62_freq_ratios();
                        match crate::composition::composition_scale_matrix(
                            &matrix, &qcomp28, &scomp28, ungapped_lambda, &freq_ratios,
                        ) {
                            Some(adj_mat) => Some((Some(adj_mat), None)),
                            None => None,
                        }
                    }
                    MatrixAdjustRule::UserSpecifiedRelEntropy
                    | MatrixAdjustRule::UnconstrainedRelEntropy
                    | MatrixAdjustRule::RelEntropyOldMatrixNewContext
                    | MatrixAdjustRule::RelEntropyOldMatrixOldContext => {
                        // Full matrix optimization (Blast_CompositionMatrixAdj)
                        let (joint_probs, first_std, second_std) =
                            crate::composition::blosum62_workspace();
                        let mut adj_matrix = matrix;
                        // NCBI uses ungappedLambda (0.3176 for BLOSUM62) for matrix scaling,
                        // NOT the gapped lambda. See matrixInfo->ungappedLambda.
                        let ungapped_lambda = crate::stat::protein_ungapped_kbp().lambda;
                        let status = crate::composition::composition_matrix_adj(
                            &mut adj_matrix, AA_SIZE, rule,
                            qn, sn, &qcomp28, &scomp28,
                            20, // RE_pseudocounts
                            0.44, // kFixedReBlosum62
                            &joint_probs, &first_std, &second_std,
                            ungapped_lambda, &matrix,
                        );
                        if status == 0 {
                            // Optimization succeeded — use adjusted matrix
                            Some((Some(adj_matrix), None))
                        } else {
                            // Fall back to lambda ratio scaling
                            let lr = crate::composition::composition_lambda_ratio(
                                &matrix, &qcomp28, &scomp28, prot_kbp.lambda,
                            );
                            Some((None, lr))
                        }
                    }
                }
            }
        } else { None };

        // If we have an adjusted matrix, re-align and recompute scores
        let (final_phits, use_adj_matrix) = if let Some((Some(ref adj_mat), _)) = adj_result {
            // Re-do gapped alignment with adjusted matrix
            let mut new_phits = Vec::new();
            for ph in &phits {
                let seed_q = (ph.query_start + ph.query_end) / 2;
                let seed_s = (ph.subject_start + ph.subject_end) / 2;
                if let Some(gr) = crate::protein::protein_gapped_align(
                    &query_aa, subj_aa, seed_q, seed_s,
                    adj_mat, gap_open, gap_extend, x_drop_final,
                ) {
                    let q_slice = &query_aa[gr.query_start..gr.query_end];
                    let s_slice = &subj_aa[gr.subject_start..gr.subject_end];
                    let (qs, ss) = gr.edit_script.render_alignment(
                        q_slice, s_slice, crate::protein::ncbistdaa_to_char,
                    );
                    new_phits.push(crate::protein_lookup::ProteinHit {
                        query_start: gr.query_start, query_end: gr.query_end,
                        subject_start: gr.subject_start, subject_end: gr.subject_end,
                        score: gr.score, num_ident: gr.num_ident,
                        align_length: gr.align_length, mismatches: gr.mismatches,
                        gap_opens: gr.gap_opens, qseq: Some(qs), sseq: Some(ss),
                    });
                }
            }
            if new_phits.is_empty() {
                (phits.clone(), false)
            } else {
                (new_phits, true)
            }
        } else {
            (phits.clone(), false)
        };

        let lambda_ratio_opt = adj_result.as_ref().and_then(|(_, lr)| *lr);

        let hsps: Vec<Hsp> = final_phits.iter().filter_map(|ph| {
            // Use Spouge FSC for e-value when Gumbel params are available.
            // This gives per-subject-length corrected e-values matching NCBI.
            let evalue = if let Some(ref gbp) = gumbel_blk {
                let base_ev = crate::stat::spouge_evalue(
                    ph.score, &prot_kbp, gbp,
                    query_aa.len() as i32, sl as i32,
                );
                if use_adj_matrix {
                    base_ev
                } else if let Some(lr) = lambda_ratio_opt {
                    // Scale the e-value by the lambda ratio
                    let scaled_kbp = crate::stat::KarlinBlk {
                        lambda: prot_kbp.lambda / lr,
                        k: prot_kbp.k,
                        log_k: prot_kbp.log_k,
                        h: prot_kbp.h,
                    };
                    crate::stat::spouge_evalue(
                        ph.score, &scaled_kbp, gbp,
                        query_aa.len() as i32, sl as i32,
                    )
                } else {
                    base_ev
                }
            } else {
                let raw_evalue = prot_kbp.raw_to_evalue(ph.score, search_space);
                if use_adj_matrix {
                    raw_evalue
                } else if let Some(lr) = lambda_ratio_opt {
                    let scaled_lambda = prot_kbp.lambda / lr;
                    search_space * prot_kbp.k * (-scaled_lambda * ph.score as f64).exp()
                } else {
                    raw_evalue
                }
            };
            if evalue > evalue_threshold { return None; }
            let (q_aln, s_aln) = if let (Some(ref qs), Some(ref ss)) = (&ph.qseq, &ph.sseq) {
                (qs.as_bytes().to_vec(), ss.as_bytes().to_vec())
            } else {
                let q_aln: Vec<u8> = (0..ph.align_length as usize).map(|i| {
                    let idx = ph.query_start + i;
                    if idx < query_aa.len() { ncbistdaa_to_ascii(query_aa[idx]) } else { b'-' }
                }).collect();
                let s_aln: Vec<u8> = (0..ph.align_length as usize).map(|i| {
                    let idx = ph.subject_start + i;
                    if idx < sl { ncbistdaa_to_ascii(subj_aa[idx]) } else { b'-' }
                }).collect();
                (q_aln, s_aln)
            };
            let midline: Vec<u8> = q_aln.iter().zip(s_aln.iter()).map(|(&q, &s)| {
                if q == s { q } else if q != b'-' && s != b'-' { b'+' } else { b' ' }
            }).collect();
            Some(Hsp {
                score: ph.score, bit_score: prot_kbp.raw_to_bit(ph.score), evalue,
                query_start: ph.query_start, query_end: ph.query_end,
                subject_start: ph.subject_start, subject_end: ph.subject_end,
                num_identities: ph.num_ident as usize, num_gaps: ph.gap_opens as usize,
                alignment_length: ph.align_length as usize,
                query_aln: q_aln, midline, subject_aln: s_aln,
                query_frame: 0, subject_frame: 0,
            })
        }).collect();

        if hsps.is_empty() { return None; }
        let hsps = if let Some(max) = max_hsps { hsps.into_iter().take(max).collect() } else { hsps };
        Some(SearchResult {
            subject_oid: oid, subject_title: title, subject_accession: accession,
            subject_len: sl, hsps, taxids: vec![],
        })
    };

    // Run sequentially or in parallel depending on num_threads.
    let mut results: Vec<SearchResult> = if num_threads <= 1 {
        (0..db.num_oids).filter_map(|oid| search_oid(oid)).collect()
    } else {
        use rayon::prelude::*;
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .stack_size(64 * 1024 * 1024)
            .build()
            .unwrap_or_else(|_| rayon::ThreadPoolBuilder::new().build().unwrap());
        pool.install(|| {
            (0..db.num_oids).into_par_iter().filter_map(|oid| search_oid(oid)).collect()
        })
    };

    results.sort_by(|a, b| a.best_evalue().partial_cmp(&b.best_evalue()).unwrap_or(std::cmp::Ordering::Equal));
    if results.len() > params.max_target_seqs {
        results.truncate(params.max_target_seqs);
    }
    results
}

/// Run a batch blastp search: multiple protein queries vs one protein database.
///
/// This is much more efficient than calling `blastp()` per query because
/// subjects are scanned once and checked against all query lookup tables.
/// Each subject is loaded into cache once, then all queries are checked.
///
/// Returns one `Vec<SearchResult>` per query, in the same order as `queries`.
pub fn blastp_batch(
    db: &BlastDb,
    queries: &[&[u8]],
    params: &SearchParams,
) -> Vec<Vec<SearchResult>> {
    if queries.is_empty() { return Vec::new(); }

    let matrix = *get_matrix(params.matrix);
    let word_size = params.word_size.max(2).min(6);
    let threshold = 11.0;

    let prot_kbp = crate::stat::lookup_protein_params(params.gap_open, params.gap_extend)
        .map(|p| crate::stat::KarlinBlk { lambda: p.lambda, k: p.k, log_k: p.k.ln(), h: p.h })
        .unwrap_or_else(crate::stat::protein_ungapped_kbp);

    let total_subj_len: usize = (0..db.num_oids)
        .map(|oid| db.get_seq_len(oid) as usize)
        .sum();

    let max_hsps = params.max_hsps;
    let evalue_threshold = params.evalue_threshold;
    // Convert bit-score parameters to raw (same as blastp())
    let ln2_b = std::f64::consts::LN_2;
    let x_drop_ungapped = (params.x_drop_ungapped as f64 * ln2_b / prot_kbp.lambda).ceil() as i32;
    let x_drop_gapped = (params.x_drop_gapped as f64 * ln2_b / prot_kbp.lambda).ceil() as i32;
    let x_drop_final = (params.x_drop_final as f64 * ln2_b / prot_kbp.lambda).ceil() as i32;
    let gap_trigger_raw = (22.0 * ln2_b / prot_kbp.lambda).ceil() as i32;
    let gap_open = params.gap_open;
    let gap_extend = params.gap_extend;

    // Prepare all queries: encode, build lookup tables, compute search spaces
    struct PreparedQuery {
        aa: Vec<u8>,
        lookup: crate::protein_lookup::ProteinLookupTable,
        search_space: f64,
    }
    let prepared: Vec<PreparedQuery> = queries.iter().map(|q| {
        let aa: Vec<u8> = q.iter()
            .map(|&b| AMINOACID_TO_NCBISTDAA[b as usize & 0x7F])
            .collect();
        let len_adj = crate::stat::compute_length_adjustment(
            aa.len() as i32, total_subj_len as i64, db.num_oids as i32, &prot_kbp);
        let search_space = crate::stat::compute_search_space(
            aa.len() as i64, total_subj_len as i64, db.num_oids as i32, len_adj);
        let lookup = crate::protein_lookup::ProteinLookupTable::build(
            &aa, word_size, &matrix, threshold);
        PreparedQuery { aa, lookup, search_space }
    }).collect();

    let num_queries = queries.len();
    let num_threads = if params.num_threads == 0 { rayon::current_num_threads() } else { params.num_threads };

    // Build merged PV — bitwise OR of all query PVs.
    // Subject positions that don't match the merged PV can't match ANY query.
    let table_refs: Vec<&crate::protein_lookup::ProteinLookupTable> =
        prepared.iter().map(|pq| &pq.lookup).collect();
    let merged_pv = crate::protein_lookup::merge_pv(&table_refs);

    let query_refs: Vec<&[u8]> = prepared.iter().map(|pq| pq.aa.as_slice()).collect();

    // Process one subject: use merged PV scan, then gapped alignment on hits.
    let process_oid = |oid: u32| -> Vec<(usize, SearchResult)> {
        let subj_raw = db.get_sequence(oid);
        let subj_len = db.get_seq_len(oid) as usize;
        if subj_len < word_size { return Vec::new(); }

        let subj_aa: Vec<u8> = subj_raw.iter().filter(|&&b| b != 0).copied().collect();
        if subj_aa.is_empty() { return Vec::new(); }

        // Batch scan: one pass over subject, checks merged PV first
        let batch_ungapped = crate::protein_lookup::batch_scan_subject(
            &query_refs, &table_refs, &merged_pv,
            &subj_aa, &matrix, x_drop_ungapped,
        );

        let mut hits_for_queries: Vec<(usize, SearchResult)> = Vec::new();

        for (qi, ungapped) in batch_ungapped {
            let pq = &prepared[qi];

            // Phase 2: gapped extension on best seed
            let ungap_cutoff = gap_trigger_raw;
            let best_seed = ungapped.iter()
                .filter(|uh| uh.score >= ungap_cutoff)
                .max_by_key(|uh| uh.score);

            let mut phits = Vec::new();
            if let Some(uh) = best_seed.filter(|uh| {
                prot_kbp.raw_to_evalue(uh.score, pq.search_space) < 10.0
            }) {
                let seed_q = (uh.query_start + uh.query_end) / 2;
                let seed_s = (uh.subject_start + uh.subject_end) / 2;
                if let Some(gr) = crate::protein::protein_gapped_align(
                    &pq.aa, &subj_aa, seed_q, seed_s,
                    &matrix, gap_open, gap_extend, x_drop_final,
                ) {
                    let q_slice = &pq.aa[gr.query_start..gr.query_end];
                    let s_slice = &subj_aa[gr.subject_start..gr.subject_end];
                    let (qs, ss) = gr.edit_script.render_alignment(
                        q_slice, s_slice, crate::protein::ncbistdaa_to_char,
                    );
                    phits.push(crate::protein_lookup::ProteinHit {
                        query_start: gr.query_start, query_end: gr.query_end,
                        subject_start: gr.subject_start, subject_end: gr.subject_end,
                        score: gr.score, num_ident: gr.num_ident,
                        align_length: gr.align_length, mismatches: gr.mismatches,
                        gap_opens: gr.gap_opens, qseq: Some(qs), sseq: Some(ss),
                    });
                }
            }
            if phits.is_empty() { phits = ungapped; }
            if phits.is_empty() { continue; }

            let best_ev = phits.iter()
                .map(|ph| prot_kbp.raw_to_evalue(ph.score, pq.search_space))
                .fold(f64::MAX, f64::min);
            if best_ev > evalue_threshold { continue; }

            let accession = db.get_accession(oid).unwrap_or_else(|| format!("oid_{}", oid));
            let title = String::from_utf8_lossy(db.get_header(oid)).to_string();
            let sl = subj_aa.len();

            let hsps: Vec<Hsp> = phits.iter().filter_map(|ph| {
                let evalue = prot_kbp.raw_to_evalue(ph.score, pq.search_space);
                if evalue > evalue_threshold { return None; }
                let (q_aln, s_aln) = if let (Some(ref qs), Some(ref ss)) = (&ph.qseq, &ph.sseq) {
                    (qs.as_bytes().to_vec(), ss.as_bytes().to_vec())
                } else {
                    let qa: Vec<u8> = (0..ph.align_length as usize).map(|i| {
                        let idx = ph.query_start + i;
                        if idx < pq.aa.len() { ncbistdaa_to_ascii(pq.aa[idx]) } else { b'-' }
                    }).collect();
                    let sa: Vec<u8> = (0..ph.align_length as usize).map(|i| {
                        let idx = ph.subject_start + i;
                        if idx < sl { ncbistdaa_to_ascii(subj_aa[idx]) } else { b'-' }
                    }).collect();
                    (qa, sa)
                };
                let midline: Vec<u8> = q_aln.iter().zip(s_aln.iter()).map(|(&q, &s)| {
                    if q == s { q } else if q != b'-' && s != b'-' { b'+' } else { b' ' }
                }).collect();
                Some(Hsp {
                    score: ph.score, bit_score: prot_kbp.raw_to_bit(ph.score), evalue,
                    query_start: ph.query_start, query_end: ph.query_end,
                    subject_start: ph.subject_start, subject_end: ph.subject_end,
                    num_identities: ph.num_ident as usize, num_gaps: ph.gap_opens as usize,
                    alignment_length: ph.align_length as usize,
                    query_aln: q_aln, midline, subject_aln: s_aln,
                    query_frame: 0, subject_frame: 0,
                })
            }).collect();

            if hsps.is_empty() { continue; }
            let hsps = if let Some(max) = max_hsps { hsps.into_iter().take(max).collect() } else { hsps };
            hits_for_queries.push((qi, SearchResult {
                subject_oid: oid, subject_title: title, subject_accession: accession,
                subject_len: sl, hsps, taxids: vec![],
            }));
        }

        hits_for_queries
    };

    // Dispatch: sequential or parallel over subjects
    let all_hits: Vec<Vec<(usize, SearchResult)>> = if num_threads <= 1 {
        (0..db.num_oids).map(|oid| process_oid(oid)).collect()
    } else {
        use rayon::prelude::*;
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .stack_size(64 * 1024 * 1024)
            .build()
            .unwrap_or_else(|_| rayon::ThreadPoolBuilder::new().build().unwrap());
        pool.install(|| {
            (0..db.num_oids).into_par_iter().map(|oid| process_oid(oid)).collect()
        })
    };

    // Scatter results into per-query buckets
    let mut results: Vec<Vec<SearchResult>> = vec![Vec::new(); num_queries];
    for oid_hits in all_hits {
        for (qi, sr) in oid_hits {
            results[qi].push(sr);
        }
    }

    // Sort each query's results by e-value and truncate
    for r in &mut results {
        r.sort_by(|a, b| a.best_evalue().partial_cmp(&b.best_evalue()).unwrap_or(std::cmp::Ordering::Equal));
        if r.len() > params.max_target_seqs {
            r.truncate(params.max_target_seqs);
        }
    }

    results
}

/// Run a blastn search (nucleotide query vs nucleotide database).
pub fn blastn_search(db: &BlastDb, query: &[u8], params: &SearchParams) -> Vec<SearchResult> {
    if query.is_empty() { return Vec::new(); }

    // Encode query to BLASTNA
    let query_plus: Vec<u8> = query.iter()
        .map(|&b| IUPACNA_TO_BLASTNA[b as usize & 0x7F])
        .collect();
    let query_minus: Vec<u8> = crate::sequence::reverse_complement(&query_plus);

    let reward = params.match_score;
    let penalty = params.mismatch;

    let ungapped_kbp = crate::stat::KarlinBlk { lambda: 1.28, k: 0.46, log_k: 0.46_f64.ln(), h: 0.85 };
    let kbp = crate::stat::nucl_gapped_kbp_lookup(
        params.gap_open, params.gap_extend, reward, penalty, &ungapped_kbp
    ).map(|(k, _)| k).unwrap_or(ungapped_kbp);

    let total_subj_len = db.total_length;
    let len_adj = crate::stat::compute_length_adjustment(
        query.len() as i32, total_subj_len as i64, db.num_oids as i32, &kbp);
    let search_space = crate::stat::compute_search_space(
        query.len() as i64, total_subj_len as i64, db.num_oids as i32, len_adj);

    let x_dropoff = params.x_drop_ungapped;

    let (q_plus, q_minus) = match params.strand.as_str() {
        "plus" => (query_plus.as_slice(), &[] as &[u8]),
        "minus" => (&[] as &[u8], query_minus.as_slice()),
        _ => (query_plus.as_slice(), query_minus.as_slice()),
    };
    let prepared_query =
        crate::search::PreparedBlastnQuery::new(q_plus, q_minus, params.word_size);
    let mut last_hit_scratch = prepared_query.last_hit_scratch();

    let mut results = Vec::new();
    for oid in 0..db.num_oids {
        let subject_packed = db.get_sequence(oid);
        let subject_len = db.get_seq_len(oid) as usize;
        if subject_len < params.word_size { continue; }

        let hsps = crate::search::blastn_ungapped_search_packed_prepared_with_scratch(
            &prepared_query,
            subject_packed,
            subject_len,
            reward, penalty, x_dropoff,
            &kbp, search_space, params.evalue_threshold,
            &mut last_hit_scratch,
        );

        if hsps.is_empty() { continue; }

        let accession = db.get_accession(oid).unwrap_or_else(|| format!("oid_{}", oid));
        let title = String::from_utf8_lossy(db.get_header(oid)).to_string();

        let api_hsps: Vec<Hsp> = hsps.iter().map(|h| {
            Hsp {
                score: h.score,
                bit_score: h.bit_score,
                evalue: h.evalue,
                query_start: h.query_start as usize,
                query_end: h.query_end as usize,
                subject_start: h.subject_start as usize,
                subject_end: h.subject_end as usize,
                num_identities: h.num_ident as usize,
                num_gaps: h.gap_opens as usize,
                alignment_length: h.align_length as usize,
                query_aln: h.qseq.as_ref().map(|s| s.as_bytes().to_vec()).unwrap_or_default(),
                subject_aln: h.sseq.as_ref().map(|s| s.as_bytes().to_vec()).unwrap_or_default(),
                midline: build_midline(
                    h.qseq.as_deref().unwrap_or(""),
                    h.sseq.as_deref().unwrap_or(""),
                ),
                query_frame: if h.context == 1 { -1 } else { 1 },
                subject_frame: 0,
            }
        }).collect();

        results.push(SearchResult {
            subject_oid: oid,
            subject_title: title,
            subject_accession: accession,
            subject_len,
            hsps: api_hsps,
            taxids: db.get_taxids(oid),
        });
    }

    results.sort_by(|a, b| a.best_evalue().partial_cmp(&b.best_evalue()).unwrap_or(std::cmp::Ordering::Equal));
    if results.len() > params.max_target_seqs {
        results.truncate(params.max_target_seqs);
    }
    results
}

/// Run a blastx search (translated nucleotide query vs protein database).
pub fn blastx(db: &BlastDb, query: &[u8], params: &SearchParams) -> Vec<SearchResult> {
    if query.len() < 3 { return Vec::new(); }

    let query_2na = ascii_to_ncbi2na(query);
    let code = crate::util::lookup_genetic_code(params.query_gencode);
    let frames = crate::util::six_frame_translation(&query_2na, code);
    let matrix = *get_matrix(params.matrix);

    let prot_kbp = crate::stat::lookup_protein_params(params.gap_open, params.gap_extend)
        .map(|p| crate::stat::KarlinBlk { lambda: p.lambda, k: p.k, log_k: p.k.ln(), h: p.h })
        .unwrap_or_else(crate::stat::protein_ungapped_kbp);

    let total_subj_len: usize = (0..db.num_oids).map(|oid| db.get_seq_len(oid) as usize).sum();
    let search_space = (query.len() as f64 / 3.0) * total_subj_len as f64;

    let word_size = params.word_size.max(2).min(6);
    let threshold = 11.0;

    let mut results = Vec::new();
    for (frame, prot) in &frames {
        if prot.len() < word_size { continue; }

        // Build lookup table once per frame (not per subject).
        let lookup_table = crate::protein_lookup::ProteinLookupTable::build(
            prot, word_size, &matrix, threshold,
        );

        for oid in 0..db.num_oids {
            let subj_raw = db.get_sequence(oid);
            let subj_aa: Vec<u8> = subj_raw.iter().filter(|&&b| b != 0).copied().collect();
            if subj_aa.is_empty() { continue; }

            let phits = crate::protein_lookup::protein_scan_with_table(
                prot, &subj_aa, &matrix, &lookup_table, params.x_drop_gapped,
            );
            for ph in &phits {
                let evalue = prot_kbp.raw_to_evalue(ph.score, search_space);
                if evalue <= params.evalue_threshold {
                    let accession = db.get_accession(oid).unwrap_or_default();
                    let title = String::from_utf8_lossy(db.get_header(oid)).to_string();
                    results.push(SearchResult {
                        subject_oid: oid,
                        subject_title: title,
                        subject_accession: accession,
                        subject_len: subj_aa.len(),
                        hsps: vec![Hsp {
                            score: ph.score,
                            bit_score: prot_kbp.raw_to_bit(ph.score),
                            evalue,
                            query_start: ph.query_start,
                            query_end: ph.query_end,
                            subject_start: ph.subject_start,
                            subject_end: ph.subject_end,
                            num_identities: ph.num_ident as usize,
                            num_gaps: ph.gap_opens as usize,
                            alignment_length: ph.align_length as usize,
                            query_aln: Vec::new(),
                            midline: Vec::new(),
                            subject_aln: Vec::new(),
                            query_frame: *frame,
                            subject_frame: 0,
                        }],
                        taxids: vec![],
                    });
                }
            }
        }
    }

    results.sort_by(|a, b| a.best_evalue().partial_cmp(&b.best_evalue()).unwrap_or(std::cmp::Ordering::Equal));
    if results.len() > params.max_target_seqs { results.truncate(params.max_target_seqs); }
    results
}

/// Run a tblastn search (protein query vs translated nucleotide database).
pub fn tblastn(db: &BlastDb, query: &[u8], params: &SearchParams) -> Vec<SearchResult> {
    if query.is_empty() { return Vec::new(); }

    let query_aa: Vec<u8> = query.iter()
        .map(|&b| AMINOACID_TO_NCBISTDAA[b as usize & 0x7F])
        .collect();
    let matrix = *get_matrix(params.matrix);
    let word_size = params.word_size.max(2).min(6);
    let threshold = 11.0;

    let prot_kbp = crate::stat::lookup_protein_params(params.gap_open, params.gap_extend)
        .map(|p| crate::stat::KarlinBlk { lambda: p.lambda, k: p.k, log_k: p.k.ln(), h: p.h })
        .unwrap_or_else(crate::stat::protein_ungapped_kbp);

    // Build lookup table once per query.
    let lookup_table = crate::protein_lookup::ProteinLookupTable::build(
        &query_aa, word_size, &matrix, threshold,
    );

    let mut results = Vec::new();
    for oid in 0..db.num_oids {
        let subject_packed = db.get_sequence(oid);
        let subject_len = db.get_seq_len(oid) as usize;

        let subj_2na: Vec<u8> = (0..subject_len).map(|i| {
            let byte = subject_packed[i >> 2];
            (byte >> (6 - 2 * (i & 3))) & 3
        }).collect();

        let frames = crate::util::six_frame_translation(&subj_2na, crate::util::lookup_genetic_code(params.db_gencode));

        for (frame, prot) in &frames {
            if prot.len() < word_size { continue; }
            let phits = crate::protein_lookup::protein_scan_with_table(
                &query_aa, prot, &matrix, &lookup_table, params.x_drop_gapped,
            );
            for ph in &phits {
                let search_space = (query_aa.len() * subject_len) as f64;
                let evalue = prot_kbp.raw_to_evalue(ph.score, search_space);
                if evalue <= params.evalue_threshold {
                    let accession = db.get_accession(oid).unwrap_or_default();
                    let title = String::from_utf8_lossy(db.get_header(oid)).to_string();
                    results.push(SearchResult {
                        subject_oid: oid,
                        subject_title: title,
                        subject_accession: accession,
                        subject_len: prot.len(),
                        hsps: vec![Hsp {
                            score: ph.score,
                            bit_score: prot_kbp.raw_to_bit(ph.score),
                            evalue,
                            query_start: ph.query_start,
                            query_end: ph.query_end,
                            subject_start: ph.subject_start,
                            subject_end: ph.subject_end,
                            num_identities: ph.num_ident as usize,
                            num_gaps: ph.gap_opens as usize,
                            alignment_length: ph.align_length as usize,
                            query_aln: Vec::new(),
                            midline: Vec::new(),
                            subject_aln: Vec::new(),
                            query_frame: 0,
                            subject_frame: *frame,
                        }],
                        taxids: vec![],
                    });
                }
            }
        }
    }

    results.sort_by(|a, b| a.best_evalue().partial_cmp(&b.best_evalue()).unwrap_or(std::cmp::Ordering::Equal));
    if results.len() > params.max_target_seqs { results.truncate(params.max_target_seqs); }
    results
}

/// Run a tblastx search (translated nt query vs translated nt database).
pub fn tblastx(db: &BlastDb, query: &[u8], params: &SearchParams) -> Vec<SearchResult> {
    if query.len() < 3 { return Vec::new(); }

    let query_2na = ascii_to_ncbi2na(query);
    let q_code = crate::util::lookup_genetic_code(params.query_gencode);
    let q_frames = crate::util::six_frame_translation(&query_2na, q_code);
    let matrix = *get_matrix(params.matrix);
    let word_size = params.word_size.max(2).min(6);
    let threshold = 11.0;

    let prot_kbp = crate::stat::lookup_protein_params(params.gap_open, params.gap_extend)
        .map(|p| crate::stat::KarlinBlk { lambda: p.lambda, k: p.k, log_k: p.k.ln(), h: p.h })
        .unwrap_or_else(crate::stat::protein_ungapped_kbp);

    let mut results = Vec::new();
    for oid in 0..db.num_oids {
        let subject_packed = db.get_sequence(oid);
        let subject_len = db.get_seq_len(oid) as usize;

        let subj_2na: Vec<u8> = (0..subject_len).map(|i| {
            let byte = subject_packed[i >> 2];
            (byte >> (6 - 2 * (i & 3))) & 3
        }).collect();

        let s_frames = crate::util::six_frame_translation(&subj_2na, crate::util::lookup_genetic_code(params.db_gencode));

        for (qframe, q_prot) in &q_frames {
            if q_prot.len() < word_size { continue; }

            for (sframe, s_prot) in &s_frames {
                if s_prot.len() < word_size { continue; }
                let phits = crate::protein_lookup::protein_scan(
                    q_prot, s_prot, &matrix, word_size, threshold, params.x_drop_gapped,
                );
                for ph in &phits {
                    let search_space = (q_prot.len() * s_prot.len()) as f64;
                    let evalue = prot_kbp.raw_to_evalue(ph.score, search_space);
                    if evalue <= params.evalue_threshold {
                        let accession = db.get_accession(oid).unwrap_or_default();
                        let title = String::from_utf8_lossy(db.get_header(oid)).to_string();
                        results.push(SearchResult {
                            subject_oid: oid,
                            subject_title: title,
                            subject_accession: accession,
                            subject_len: s_prot.len(),
                            hsps: vec![Hsp {
                                score: ph.score,
                                bit_score: prot_kbp.raw_to_bit(ph.score),
                                evalue,
                                query_start: ph.query_start,
                                query_end: ph.query_end,
                                subject_start: ph.subject_start,
                                subject_end: ph.subject_end,
                                num_identities: ph.num_ident as usize,
                                num_gaps: ph.gap_opens as usize,
                                alignment_length: ph.align_length as usize,
                                query_aln: Vec::new(),
                                midline: Vec::new(),
                                subject_aln: Vec::new(),
                                query_frame: *qframe,
                                subject_frame: *sframe,
                            }],
                            taxids: vec![],
                        });
                    }
                }
            }
        }
    }

    results.sort_by(|a, b| a.best_evalue().partial_cmp(&b.best_evalue()).unwrap_or(std::cmp::Ordering::Equal));
    if results.len() > params.max_target_seqs { results.truncate(params.max_target_seqs); }
    results
}

// ── Utility functions ───────────────────────────────────────────────────────

/// Parse a multi-FASTA byte slice into (title, sequence) pairs.
pub fn parse_fasta(input: &[u8]) -> Vec<(String, Vec<u8>)> {
    let mut sequences = Vec::new();
    let mut current_title = String::new();
    let mut current_seq: Vec<u8> = Vec::new();

    for line in input.split(|&b| b == b'\n') {
        let line = line.strip_suffix(b"\r").unwrap_or(line);
        if line.is_empty() { continue; }
        if line.starts_with(b">") {
            if !current_title.is_empty() {
                sequences.push((current_title.clone(), current_seq.clone()));
                current_seq.clear();
            }
            current_title = String::from_utf8_lossy(&line[1..]).trim().to_string();
        } else {
            current_seq.extend_from_slice(line);
        }
    }
    if !current_title.is_empty() || !current_seq.is_empty() {
        sequences.push((current_title, current_seq));
    }
    sequences
}

/// Reverse complement an ASCII nucleotide sequence.
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| match b {
        b'A' | b'a' => b'T', b'T' | b't' => b'A',
        b'C' | b'c' => b'G', b'G' | b'g' => b'C',
        b'N' | b'n' => b'N',
        b'R' => b'Y', b'Y' => b'R',
        b'M' => b'K', b'K' => b'M',
        b'W' => b'W', b'S' => b'S',
        b'B' => b'V', b'V' => b'B',
        b'D' => b'H', b'H' => b'D',
        _ => b'N',
    }).collect()
}

/// Six-frame translation of an ASCII nucleotide sequence.
/// Returns array of 6 TranslatedFrame structs.
pub fn six_frame_translate(nt_seq: &[u8]) -> [TranslatedFrame; 6] {
    let nt_2na = ascii_to_ncbi2na(nt_seq);
    let frames = crate::util::six_frame_translation(&nt_2na, &crate::util::STANDARD_GENETIC_CODE);
    let mut result: Vec<TranslatedFrame> = frames.into_iter().map(|(frame, prot)| {
        // Convert NCBIstdaa to ASCII
        let ascii: Vec<u8> = prot.iter().map(|&b| ncbistdaa_to_ascii(b)).collect();
        TranslatedFrame { frame, protein: ascii, nt_len: nt_seq.len() }
    }).collect();
    // Ensure exactly 6 frames
    while result.len() < 6 {
        result.push(TranslatedFrame { frame: 0, protein: Vec::new(), nt_len: nt_seq.len() });
    }
    [
        result.remove(0), result.remove(0), result.remove(0),
        result.remove(0), result.remove(0), result.remove(0),
    ]
}

/// A translated reading frame.
pub struct TranslatedFrame {
    pub frame: i32,
    pub protein: Vec<u8>,
    pub nt_len: usize,
}

// ── Additional types ────────────────────────────────────────────────────────

/// Parsed BLAST database header/defline.
#[derive(Debug, Clone)]
pub struct BlastDefLine {
    pub title: String,
    pub accession: String,
    pub taxid: u32,
}

/// Apply SEG masking on NCBIstdaa-encoded sequence in place.
pub fn apply_seg_ncbistdaa(seq: &mut [u8]) {
    // SEG on raw NCBIstdaa: convert to ASCII, mask, convert back
    let mut ascii: Vec<u8> = seq.iter()
        .map(|&b| if (b as usize) < AA_SIZE { NCBISTDAA_TO_AMINOACID[b as usize] as u8 } else { b'X' })
        .collect();
    apply_seg(&mut ascii);
    for (i, &b) in ascii.iter().enumerate() {
        if b == b'X' && seq[i] != AMINOACID_TO_NCBISTDAA[b'X' as usize & 0x7F] {
            seq[i] = AMINOACID_TO_NCBISTDAA[b'X' as usize & 0x7F]; // masked
        }
    }
}

/// Compute a boolean mask indicating which positions are lowercase.
pub fn lowercase_mask(seq: &[u8]) -> Vec<bool> {
    seq.iter().map(|b| b.is_ascii_lowercase()).collect()
}

/// Apply repeat masking based on n-mer frequency.
pub fn apply_repeat_mask(seq: &mut [u8]) {
    let mask = repeat_mask(seq, 11, 2.0);
    for (i, masked) in mask.iter().enumerate() {
        if *masked { seq[i] = b'N'; }
    }
}

/// Compute repeat mask: positions where n-mer frequency exceeds threshold.
pub fn repeat_mask(seq: &[u8], nmer_size: usize, threshold: f64) -> Vec<bool> {
    let mut mask = vec![false; seq.len()];
    if seq.len() < nmer_size { return mask; }

    // Count n-mer occurrences
    let mut counts = std::collections::HashMap::new();
    for i in 0..=seq.len() - nmer_size {
        *counts.entry(&seq[i..i + nmer_size]).or_insert(0u32) += 1;
    }

    // Mark positions where n-mer count exceeds threshold * average
    let avg = seq.len() as f64 / 4f64.powi(nmer_size as i32).max(1.0);
    let cutoff = (threshold * avg).max(2.0) as u32;
    for i in 0..=seq.len() - nmer_size {
        if counts[&seq[i..i + nmer_size]] >= cutoff {
            for j in i..i + nmer_size {
                mask[j] = true;
            }
        }
    }
    mask
}

// ── BlastnSearch builder (moved from blastn.rs) ─────────────────────────────

/// Which query strand(s) to search.
#[derive(Clone, Copy, PartialEq)]
pub enum Strand {
    Both,
    Plus,
    Minus,
}

/// Builder for configuring and running a blastn search.
///
/// # Example
///
/// ```no_run
/// use blast_rs::BlastnSearch;
///
/// let results = BlastnSearch::new()
///     .query(b"ACGTACGTACGTACGTACGTACGTACGT")
///     .subject(b"NNNNNACGTACGTACGTACGTACGTACGTNNNNN")
///     .run();
///
/// for hit in &results {
///     println!("score={} evalue={:.2e}", hit.score, hit.evalue);
/// }
/// ```
pub struct BlastnSearch {
    pub word_size: usize,
    pub reward: i32,
    pub penalty: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub evalue: f64,
    pub dust: bool,
    pub strand: Strand,
    pub xdrop_gap_final: f64,
    query_raw: Vec<u8>,
    subject_raw: Vec<u8>,
}

impl Default for BlastnSearch {
    fn default() -> Self { Self::new() }
}

impl BlastnSearch {
    pub fn new() -> Self {
        BlastnSearch {
            word_size: 11, reward: 1, penalty: -3,
            gap_open: 5, gap_extend: 2, evalue: 10.0,
            dust: true, strand: Strand::Both, xdrop_gap_final: 100.0,
            query_raw: Vec::new(), subject_raw: Vec::new(),
        }
    }

    pub fn query(mut self, seq: &[u8]) -> Self { self.query_raw = seq.to_vec(); self }
    pub fn subject(mut self, seq: &[u8]) -> Self { self.subject_raw = seq.to_vec(); self }
    pub fn word_size(mut self, ws: usize) -> Self { self.word_size = ws; self }
    pub fn reward(mut self, r: i32) -> Self { self.reward = r; self }
    pub fn penalty(mut self, p: i32) -> Self { self.penalty = p; self }
    pub fn gap_open(mut self, g: i32) -> Self { self.gap_open = g; self }
    pub fn gap_extend(mut self, g: i32) -> Self { self.gap_extend = g; self }
    pub fn evalue(mut self, e: f64) -> Self { self.evalue = e; self }
    pub fn dust(mut self, d: bool) -> Self { self.dust = d; self }
    pub fn strand(mut self, s: Strand) -> Self { self.strand = s; self }

    pub fn run(&self) -> Vec<SearchHsp> {
        if self.query_raw.is_empty() || self.subject_raw.is_empty() {
            return Vec::new();
        }

        let mut query_plus: Vec<u8> = self.query_raw.iter()
            .map(|&b| blastn_iupac_to_blastna(b)).collect();

        if self.dust {
            let mask = crate::filter::dust_filter(&query_plus, 64, 2.0);
            mask.apply(&mut query_plus, 14);
        }

        let query_plus_nomask: Vec<u8> = self.query_raw.iter()
            .map(|&b| blastn_iupac_to_blastna(b)).collect();
        let query_minus: Vec<u8> = query_plus.iter().rev()
            .map(|&b| blastn_complement(b)).collect();
        let query_minus_nomask: Vec<u8> = query_plus_nomask.iter().rev()
            .map(|&b| blastn_complement(b)).collect();

        let qp = if self.strand != Strand::Minus { &query_plus[..] } else { &[] };
        let qm = if self.strand != Strand::Plus { &query_minus[..] } else { &[] };
        let qpn = if self.strand != Strand::Minus { &query_plus_nomask[..] } else { &[] };
        let qmn = if self.strand != Strand::Plus { &query_minus_nomask[..] } else { &[] };

        let subject: Vec<u8> = self.subject_raw.iter()
            .map(|&b| blastn_iupac_to_blastna(b)).collect();

        let m = build_blastna_matrix(self.reward, self.penalty);
        let matrix_fn = |i: usize, j: usize| -> i32 { m[i][j] };
        let mut lo = i32::MAX;
        let mut hi = i32::MIN;
        for i in 0..16 { for j in 0..16 {
            let s = m[i][j];
            if s <= -100000000 || s >= 100000000 { continue; }
            if s < lo { lo = s; }
            if s > hi { hi = s; }
        }}

        let ambig: &[u8] = &[14, 15];
        let ctx = UngappedKbpContext { query_offset: 0, query_length: query_plus.len() as i32, is_valid: true };
        let kbp_results = ungapped_kbp_calc(&query_plus, &[ctx], lo, hi, 16, ambig, &matrix_fn);
        let ungapped_kbp = kbp_results[0].clone().unwrap_or(KarlinBlk {
            lambda: 1.374, k: 0.621, log_k: 0.621_f64.ln(), h: 1.286,
        });
        let (gapped_kbp, _) = nucl_gapped_kbp_lookup(
            self.gap_open, self.gap_extend, self.reward, self.penalty, &ungapped_kbp,
        ).unwrap_or((ungapped_kbp.clone(), false));

        let (alpha, beta) = nucl_alpha_beta(
            self.reward, self.penalty, self.gap_open, self.gap_extend,
            ungapped_kbp.lambda, ungapped_kbp.h, true,
        );
        let db_length = self.subject_raw.len() as i64;
        let (len_adj, _) = compute_length_adjustment_exact(
            gapped_kbp.k, gapped_kbp.log_k,
            alpha / gapped_kbp.lambda, beta,
            self.query_raw.len() as i32, db_length, 1,
        );
        let eff_db = (db_length - len_adj as i64).max(1);
        let search_space = eff_db as f64 * (self.query_raw.len() as i32 - len_adj).max(1) as f64;
        let x_dropoff = (self.xdrop_gap_final * std::f64::consts::LN_2 / gapped_kbp.lambda) as i32;

        blastn_gapped_search_nomask(
            qp, qm, qpn, qmn,
            &subject, self.word_size, self.reward, self.penalty,
            self.gap_open, self.gap_extend, x_dropoff,
            &gapped_kbp, search_space, self.evalue,
        )
    }
}

fn blastn_iupac_to_blastna(b: u8) -> u8 {
    match b {
        b'A' | b'a' => 0, b'C' | b'c' => 1, b'G' | b'g' => 2, b'T' | b't' => 3,
        b'R' | b'r' => 4, b'Y' | b'y' => 5, b'M' | b'm' => 6, b'K' | b'k' => 7,
        b'W' | b'w' => 8, b'S' | b's' => 9, b'B' | b'b' => 10, b'D' | b'd' => 11,
        b'H' | b'h' => 12, b'V' | b'v' => 13, b'N' | b'n' => 14, _ => 15,
    }
}

fn blastn_complement(b: u8) -> u8 {
    match b {
        0 => 3, 1 => 2, 2 => 1, 3 => 0,
        4 => 5, 5 => 4, 6 => 7, 7 => 6,
        8 => 8, 9 => 9, 10 => 13, 11 => 12,
        12 => 11, 13 => 10, 14 => 14, _ => 15,
    }
}

// ── Internal helpers ────────────────────────────────────────────────────────

fn ascii_to_ncbi2na(seq: &[u8]) -> Vec<u8> {
    seq.iter().map(|&b| match b {
        b'A' | b'a' => 0, b'C' | b'c' => 1,
        b'G' | b'g' => 2, b'T' | b't' => 3, _ => 0,
    }).collect()
}

fn ncbistdaa_to_ascii(b: u8) -> u8 {
    if (b as usize) < AA_SIZE {
        NCBISTDAA_TO_AMINOACID[b as usize] as u8
    } else {
        b'X'
    }
}

fn build_midline(qseq: &str, sseq: &str) -> Vec<u8> {
    qseq.bytes().zip(sseq.bytes()).map(|(q, s)| {
        if q == s { b'|' } else { b' ' }
    }).collect()
}

// ── Aliases for API compatibility ───────────────────────────────────────────

/// Get the scoring matrix for a given MatrixType.
pub fn get_matrix(mt: MatrixType) -> &'static [[i32; AA_SIZE]; AA_SIZE] {
    match mt {
        MatrixType::Blosum45 => &crate::matrix::BLOSUM45,
        MatrixType::Blosum50 => &crate::matrix::BLOSUM50,
        MatrixType::Blosum62 => &crate::matrix::BLOSUM62,
        MatrixType::Blosum80 => &crate::matrix::BLOSUM80,
        MatrixType::Blosum90 => &crate::matrix::BLOSUM90,
        MatrixType::Pam30 => &crate::matrix::PAM30,
        MatrixType::Pam70 => &crate::matrix::PAM70,
        MatrixType::Pam250 => &crate::matrix::PAM250,
    }
}

/// Get a genetic code translation table by NCBI code number.
/// Returns a 64-byte table mapping codons to NCBIstdaa amino acid codes.
/// Supports all standard NCBI genetic codes (1-33).
pub fn get_codon_table(code: u8) -> &'static [u8; 64] {
    crate::util::lookup_genetic_code(code)
}

// ── Masking wrappers ────────────────────────────────────────────────────────

/// Apply DUST low-complexity masking in place (replaces masked positions with N).
pub fn apply_dust(seq: &mut [u8]) {
    let mask = crate::filter::dust_filter(seq, 64, 2.0);
    for r in &mask.regions {
        let start = r.start.max(0) as usize;
        let end = (r.end as usize).min(seq.len());
        for i in start..end {
            seq[i] = b'N';
        }
    }
}

/// Apply SEG low-complexity masking in place (replaces masked positions with X).
pub fn apply_seg(seq: &mut [u8]) {
    let mask = crate::filter::seg_filter(seq, 12, 2.2);
    for r in &mask.regions {
        let start = r.start.max(0) as usize;
        let end = (r.end as usize).min(seq.len());
        for i in start..end {
            seq[i] = b'X';
        }
    }
}

/// Replace lowercase characters with N (nucleotide masking).
pub fn apply_lowercase_mask_nucleotide(seq: &mut [u8]) {
    for b in seq.iter_mut() {
        if b.is_ascii_lowercase() { *b = b'N'; }
    }
}

/// Replace lowercase characters with X (protein masking).
pub fn apply_lowercase_mask_protein(seq: &mut [u8]) {
    for b in seq.iter_mut() {
        if b.is_ascii_lowercase() { *b = b'X'; }
    }
}

// ── Composition statistics ──────────────────────────────────────────────────

/// Compute amino acid composition (frequencies) from NCBIstdaa-encoded sequence.
pub fn composition_ncbistdaa(seq: &[u8]) -> [f64; 28] {
    let mut counts = [0u64; 28];
    for &b in seq {
        if (b as usize) < 28 { counts[b as usize] += 1; }
    }
    let total = seq.len().max(1) as f64;
    let mut freqs = [0.0; 28];
    for i in 0..28 { freqs[i] = counts[i] as f64 / total; }
    freqs
}

/// Adjust an E-value using per-sequence composition correction (mode 1).
///
/// q and r are amino acid frequency vectors (28 elements, NCBIstdaa indexed).
/// Returns the adjusted E-value, or the original if correction is inapplicable.
#[allow(clippy::too_many_arguments)]
pub fn adjust_evalue(
    raw_evalue: f64,
    score: i32,
    q: &[f64; 28],
    r: &[f64; 28],
    matrix: &ScoringMatrix,
    lambda_standard: f64,
    k: f64,
    eff_query_len: usize,
    eff_db_len: u64,
) -> f64 {
    match find_adjusted_lambda(q, r, matrix, lambda_standard) {
        None => raw_evalue,
        Some(lambda_prime) => {
            (eff_query_len as f64) * (eff_db_len as f64) * k
                * (-(lambda_prime * score as f64)).exp()
        }
    }
}

/// Adjust E-value with mode selection:
///   0 = no adjustment (returns raw_evalue)
///   1 = unconditional composition-based lambda adjustment
///   2 = conditional: only apply if composition diverges significantly
///   3 = unconditional (always applied, even if expected_score >= 0)
#[allow(clippy::too_many_arguments)]
pub fn adjust_evalue_with_mode(
    raw_evalue: f64,
    score: i32,
    q: &[f64; 28],
    r: &[f64; 28],
    matrix: &ScoringMatrix,
    lambda_standard: f64,
    k: f64,
    eff_query_len: usize,
    eff_db_len: u64,
    mode: u8,
) -> f64 {
    match mode {
        0 => raw_evalue,
        1 => adjust_evalue(raw_evalue, score, q, r, matrix, lambda_standard, k, eff_query_len, eff_db_len),
        2 => {
            let mu_bg = expected_score_with_bg(matrix);
            let mu_actual = compo_expected_score(q, r, matrix);
            let threshold = -0.2 * lambda_standard;
            if (mu_actual - mu_bg).abs() > threshold.abs() {
                adjust_evalue(raw_evalue, score, q, r, matrix, lambda_standard, k, eff_query_len, eff_db_len)
            } else {
                raw_evalue
            }
        }
        3 => {
            let eval_sum = |lam: f64| -> f64 {
                let mut sum = 0.0f64;
                for i in 1..23usize {
                    for j in 1..23usize {
                        let s = matrix.score(i as u8, j as u8) as f64;
                        sum += q[i] * r[j] * (lam * s).exp();
                    }
                }
                sum
            };
            let mut lo = 0.0f64;
            let mut hi = lambda_standard * 4.0;
            if eval_sum(hi) >= 1.0 || eval_sum(lo) < 1.0 { return raw_evalue; }
            for _ in 0..60 {
                let mid = (lo + hi) / 2.0;
                if eval_sum(mid) > 1.0 { lo = mid; } else { hi = mid; }
                if hi - lo < 1e-10 { break; }
            }
            let lambda_prime = (lo + hi) / 2.0;
            (eff_query_len as f64) * (eff_db_len as f64) * k
                * (-(lambda_prime * score as f64)).exp()
        }
        _ => raw_evalue,
    }
}

/// Find adjusted lambda via bisection for composition-based statistics.
fn find_adjusted_lambda(
    q: &[f64; 28], r: &[f64; 28], matrix: &ScoringMatrix, lambda_standard: f64,
) -> Option<f64> {
    if compo_expected_score(q, r, matrix) >= 0.0 { return None; }
    let eval_sum = |lam: f64| -> f64 {
        let mut sum = 0.0f64;
        for i in 1..23usize {
            for j in 1..23usize {
                sum += q[i] * r[j] * (lam * matrix.score(i as u8, j as u8) as f64).exp();
            }
        }
        sum
    };
    if eval_sum(0.0) < 1.0 { return None; }
    let mut lo = 0.0f64;
    let mut hi = lambda_standard * 4.0;
    if eval_sum(hi) >= 1.0 { return None; }
    for _ in 0..60 {
        let mid = (lo + hi) / 2.0;
        if eval_sum(mid) > 1.0 { lo = mid; } else { hi = mid; }
        if hi - lo < 1e-10 { break; }
    }
    Some((lo + hi) / 2.0)
}

fn compo_expected_score(q: &[f64; 28], r: &[f64; 28], matrix: &ScoringMatrix) -> f64 {
    let mut mu = 0.0f64;
    for i in 1..23usize {
        for j in 1..23usize {
            mu += q[i] * r[j] * matrix.score(i as u8, j as u8) as f64;
        }
    }
    mu
}

fn expected_score_with_bg(matrix: &ScoringMatrix) -> f64 {
    let bg = &crate::matrix::AA_FREQUENCIES;
    // AA_FREQUENCIES is [f64; 20] in ACDEFGHIKLMNPQRSTVWY order, need to map to NCBIstdaa
    let ncbi_order = [1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22]; // NCBIstdaa indices
    let mut freq28 = [0.0f64; 28];
    for (k, &idx) in ncbi_order.iter().enumerate() {
        if k < bg.len() { freq28[idx] = bg[k]; }
    }
    compo_expected_score(&freq28, &freq28, matrix)
}

// ── Low-level PSSM functions ────────────────────────────────────────────────

/// Build a PSSM from search results (for PSI-BLAST iteration).
pub fn build_pssm(
    query: &[u8],
    results: &[SearchResult],
    inclusion_evalue: f64,
    matrix_type: MatrixType,
    _lambda: f64,
) -> crate::pssm::Pssm {
    let query_aa: Vec<u8> = query.iter()
        .map(|&b| AMINOACID_TO_NCBISTDAA[b as usize & 0x7F])
        .collect();
    let matrix = get_matrix(matrix_type);
    let mut pssm = crate::pssm::Pssm::from_sequence(&query_aa, matrix);

    // Collect aligned subject sequences from results that pass inclusion threshold
    let aligned: Vec<Vec<u8>> = results.iter()
        .flat_map(|r| r.hsps.iter())
        .filter(|hsp| hsp.evalue <= inclusion_evalue)
        .filter_map(|hsp| {
            if !hsp.subject_aln.is_empty() {
                Some(hsp.subject_aln.iter()
                    .map(|&b| AMINOACID_TO_NCBISTDAA[b as usize & 0x7F])
                    .collect())
            } else { None }
        })
        .collect();

    if !aligned.is_empty() {
        pssm.update_from_alignment(&aligned, &crate::matrix::AA_FREQUENCIES, 0.5);
    }
    pssm
}

/// Search a database using a PSSM instead of a substitution matrix.
pub fn search_with_pssm(
    db: &BlastDb,
    _query: &[u8],
    pssm: &crate::pssm::Pssm,
    params: &SearchParams,
) -> Vec<SearchResult> {
    let prot_kbp = crate::stat::lookup_protein_params(params.gap_open, params.gap_extend)
        .map(|p| crate::stat::KarlinBlk { lambda: p.lambda, k: p.k, log_k: p.k.ln(), h: p.h })
        .unwrap_or_else(crate::stat::protein_ungapped_kbp);

    let total_subj_len: usize = (0..db.num_oids)
        .map(|oid| db.get_seq_len(oid) as usize)
        .sum();
    let search_space = (pssm.length as f64) * (total_subj_len as f64);

    let subj_pairs: Vec<(String, Vec<u8>)> = (0..db.num_oids).map(|oid| {
        let acc = db.get_accession(oid).unwrap_or_else(|| format!("oid_{}", oid));
        let raw = db.get_sequence(oid);
        let aa: Vec<u8> = raw.iter().filter(|&&b| b != 0).copied().collect();
        (acc, aa)
    }).collect();

    let hits = crate::pssm::psi_blast_iteration(
        pssm, &subj_pairs,
        params.evalue_threshold,
        search_space,
        prot_kbp.lambda,
        prot_kbp.k,
    );

    hits.iter().map(|h| {
        SearchResult {
            subject_oid: subj_pairs.iter().position(|(id, _)| id == &h.subject_id).unwrap_or(0) as u32,
            subject_title: h.subject_id.clone(),
            subject_accession: h.subject_id.clone(),
            subject_len: h.subject_len,
            hsps: vec![Hsp {
                score: h.score,
                bit_score: prot_kbp.raw_to_bit(h.score),
                evalue: h.evalue,
                query_start: 0,
                query_end: h.align_len,
                subject_start: h.subject_start,
                subject_end: h.subject_start + h.align_len,
                num_identities: 0,
                num_gaps: 0,
                alignment_length: h.align_len,
                query_aln: Vec::new(),
                midline: Vec::new(),
                subject_aln: Vec::new(),
                query_frame: 0,
                subject_frame: 0,
            }],
            taxids: vec![],
        }
    }).collect()
}

/// Alias for `blastn_search` matching the old API.
/// Note: at the crate root this is available as `blastn_search` since the `blastn` name
/// is used by the builder-pattern module. Use `blast_rs::api::blastn()` for the function form.
pub fn blastn(db: &BlastDb, query: &[u8], params: &SearchParams) -> Vec<SearchResult> {
    blastn_search(db, query, params)
}

/// Alias for `blastp` — low-level protein search backend.
pub fn blast_search(db: &BlastDb, query: &[u8], params: &SearchParams) -> Vec<SearchResult> {
    blastp(db, query, params)
}

/// Low-level blastx search (alias for `blastx`).
pub fn blastx_search(db: &BlastDb, query: &[u8], params: &SearchParams) -> Vec<SearchResult> {
    blastx(db, query, params)
}

/// Low-level tblastn search (alias for `tblastn`).
pub fn tblastn_search(db: &BlastDb, query: &[u8], params: &SearchParams) -> Vec<SearchResult> {
    tblastn(db, query, params)
}

/// Low-level tblastx search (alias for `tblastx`).
pub fn tblastx_search(db: &BlastDb, query: &[u8], params: &SearchParams) -> Vec<SearchResult> {
    tblastx(db, query, params)
}

/// Low-level PSI-BLAST iteration (alias for `psiblast`).
pub fn psiblast_search(db: &BlastDb, query: &[u8], params: &PsiblastParams) -> (Vec<SearchResult>, crate::pssm::Pssm) {
    psiblast(db, query, params)
}

/// Six-frame translation with a specific genetic code number (1=standard, 2=mito, etc.).
pub fn six_frame_translate_with_code(nt_seq: &[u8], genetic_code: u8) -> [TranslatedFrame; 6] {
    let nt_2na = ascii_to_ncbi2na(nt_seq);
    let table = crate::util::lookup_genetic_code(genetic_code);
    let frames = crate::util::six_frame_translation(&nt_2na, table);
    let mut result: Vec<TranslatedFrame> = frames.into_iter().map(|(frame, prot)| {
        let ascii: Vec<u8> = prot.iter().map(|&b| ncbistdaa_to_ascii(b)).collect();
        TranslatedFrame { frame, protein: ascii, nt_len: nt_seq.len() }
    }).collect();
    while result.len() < 6 {
        result.push(TranslatedFrame { frame: 0, protein: Vec::new(), nt_len: nt_seq.len() });
    }
    [
        result.remove(0), result.remove(0), result.remove(0),
        result.remove(0), result.remove(0), result.remove(0),
    ]
}

// ── PSI-BLAST ───────────────────────────────────────────────────────────────

/// PSI-BLAST specific parameters.
#[derive(Debug, Clone)]
pub struct PsiblastParams {
    pub search: SearchParams,
    /// Number of PSI-BLAST iterations (default: 3).
    pub num_iterations: u32,
    /// E-value threshold for including hits in PSSM construction (default: 0.001).
    pub inclusion_evalue: f64,
}

impl PsiblastParams {
    pub fn new(search: SearchParams) -> Self {
        PsiblastParams {
            search,
            num_iterations: 3,
            inclusion_evalue: 0.001,
        }
    }
    pub fn num_iterations(mut self, v: u32) -> Self { self.num_iterations = v; self }
    pub fn inclusion_evalue(mut self, v: f64) -> Self { self.inclusion_evalue = v; self }
}

/// Run iterative PSI-BLAST search.
/// Returns final-round hits and the resulting PSSM.
pub fn psiblast(db: &BlastDb, query: &[u8], params: &PsiblastParams) -> (Vec<SearchResult>, crate::pssm::Pssm) {
    let query_aa: Vec<u8> = query.iter()
        .map(|&b| AMINOACID_TO_NCBISTDAA[b as usize & 0x7F])
        .collect();
    let matrix = *get_matrix(params.search.matrix);

    // Initial PSSM from query
    let mut pssm = crate::pssm::Pssm::from_sequence(&query_aa, &matrix);

    let prot_kbp = crate::stat::lookup_protein_params(
        params.search.gap_open, params.search.gap_extend)
        .map(|p| crate::stat::KarlinBlk { lambda: p.lambda, k: p.k, log_k: p.k.ln(), h: p.h })
        .unwrap_or_else(crate::stat::protein_ungapped_kbp);

    let total_subj_len: usize = (0..db.num_oids)
        .map(|oid| db.get_seq_len(oid) as usize)
        .sum();
    let search_space = (query_aa.len() as f64) * (total_subj_len as f64);

    let mut final_results = Vec::new();

    // Build subject pairs for PSI-BLAST iteration
    let subj_pairs: Vec<(String, Vec<u8>)> = (0..db.num_oids).map(|oid| {
        let acc = db.get_accession(oid).unwrap_or_else(|| format!("oid_{}", oid));
        let raw = db.get_sequence(oid);
        let aa: Vec<u8> = raw.iter().filter(|&&b| b != 0).copied().collect();
        (acc, aa)
    }).collect();

    for _iter in 0..params.num_iterations {
        let hits = crate::pssm::psi_blast_iteration(
            &pssm, &subj_pairs,
            params.inclusion_evalue,
            search_space,
            prot_kbp.lambda,
            prot_kbp.k,
        );

        if hits.is_empty() { break; }

        final_results = hits.iter().map(|h| {
            SearchResult {
                subject_oid: subj_pairs.iter().position(|(id, _)| id == &h.subject_id).unwrap_or(0) as u32,
                subject_title: h.subject_id.clone(),
                subject_accession: h.subject_id.clone(),
                subject_len: h.subject_len,
                hsps: vec![Hsp {
                    score: h.score,
                    bit_score: prot_kbp.raw_to_bit(h.score),
                    evalue: h.evalue,
                    query_start: 0,
                    query_end: h.align_len,
                    subject_start: h.subject_start,
                    subject_end: h.subject_start + h.align_len,
                    num_identities: 0,
                    num_gaps: 0,
                    alignment_length: h.align_len,
                    query_aln: Vec::new(),
                    midline: Vec::new(),
                    subject_aln: Vec::new(),
                    query_frame: 0,
                    subject_frame: 0,
                }],
                taxids: vec![],
            }
        }).collect();

        // Update PSSM from aligned sequences
        let aligned: Vec<Vec<u8>> = hits.iter().filter_map(|h| {
            subj_pairs.iter()
                .find(|(id, _)| id == &h.subject_id)
                .and_then(|(_, seq)| {
                    if h.subject_start >= seq.len() { return None; }
                    let end = (h.subject_start + h.align_len).min(seq.len());
                    Some(seq[h.subject_start..end].to_vec())
                })
        }).collect();

        if !aligned.is_empty() {
            pssm.update_from_alignment(&aligned, &crate::matrix::AA_FREQUENCIES, 0.5);
        }
    }

    (final_results, pssm)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_builder_defaults() {
        let s = BlastnSearch::new();
        assert_eq!(s.word_size, 11);
        assert_eq!(s.reward, 1);
        assert_eq!(s.penalty, -3);
        assert_eq!(s.evalue, 10.0);
    }

    #[test]
    fn test_perfect_match() {
        let results = BlastnSearch::new()
            .word_size(7)
            .dust(false)
            .query(b"ACGTACGTACGTACGTACGTACGT")
            .subject(b"NNNNNACGTACGTACGTACGTACGTNNNNN")
            .run();
        assert!(!results.is_empty(), "Should find perfect match");
        assert!(results[0].score > 0);
    }

    #[test]
    fn test_no_match() {
        let results = BlastnSearch::new()
            .word_size(7)
            .dust(false)
            .strand(Strand::Plus)
            .query(b"AAAAAAAAAAAAAAAA")
            .subject(b"CCCCCCCCCCCCCCCC")
            .run();
        assert!(results.is_empty(), "Should find no match");
    }

    #[test]
    fn test_custom_scoring() {
        let results = BlastnSearch::new()
            .word_size(7)
            .reward(2)
            .penalty(-3)
            .gap_open(5)
            .gap_extend(2)
            .dust(false)
            .query(b"ACGTACGTACGTACGT")
            .subject(b"ACGTACGTACGTACGT")
            .run();
        assert!(!results.is_empty());
    }
}
