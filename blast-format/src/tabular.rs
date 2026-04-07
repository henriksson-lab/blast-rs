//! BLAST tabular output format (-outfmt 6/7).
//!
//! Default columns: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore

use std::io::Write;

/// A single alignment hit for tabular output.
pub struct TabularHit {
    pub query_id: String,
    pub subject_id: String,
    pub pct_identity: f64,
    pub align_len: i32,
    pub mismatches: i32,
    pub gap_opens: i32,
    pub query_start: i32,
    pub query_end: i32,
    pub subject_start: i32,
    pub subject_end: i32,
    pub evalue: f64,
    pub bit_score: f64,
}

/// Format an e-value matching BLAST reference output style.
/// Values >= 0.1 use fixed-point (e.g. "1.7"), smaller use scientific ("3.51e-23").
fn format_evalue(val: f64) -> String {
    if val == 0.0 {
        "0.0".to_string()
    } else if val >= 1e100 {
        format!("{:.2e}", val)
    } else if val >= 0.01 {
        // Use 2 significant digits like BLAST reference
        let digits = if val >= 100.0 { 0usize }
            else if val >= 10.0 { 1 }
            else if val >= 1.0 { 1 }
            else { 2 };
        let s = format!("{:.prec$}", val, prec = digits);
        s
    } else {
        format!("{:.2e}", val)
    }
}

/// Write tabular output (outfmt 6) for a set of hits.
pub fn format_tabular<W: Write>(writer: &mut W, hits: &[TabularHit]) -> std::io::Result<()> {
    for hit in hits {
        writeln!(
            writer,
            "{}\t{}\t{:.3}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.1}",
            hit.query_id,
            hit.subject_id,
            hit.pct_identity,
            hit.align_len,
            hit.mismatches,
            hit.gap_opens,
            hit.query_start,
            hit.query_end,
            hit.subject_start,
            hit.subject_end,
            format_evalue(hit.evalue),
            hit.bit_score,
        )?;
    }
    Ok(())
}
