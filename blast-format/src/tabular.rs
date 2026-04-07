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

/// Format an e-value matching BLAST tabular output (-outfmt 6).
/// Uses NStr::DoubleToString(evalue, 2, fDoubleScientific) which is
/// essentially %.2e for small values, transitioning to fixed for larger ones.
fn format_evalue(val: f64) -> String {
    if val < 1.0e-180 {
        "0.0".to_string()
    } else if val < 0.001 {
        // Scientific notation with 2 decimal places in mantissa
        format!("{:.2e}", val)
    } else if val < 0.1 {
        // Fixed with 3 decimal places
        format!("{:.3}", val)
    } else if val < 1.0 {
        // Fixed with 2 decimal places
        format!("{:.2}", val)
    } else if val < 10.0 {
        // Fixed with 1 decimal place
        format!("{:.1}", val)
    } else {
        // Integer
        format!("{:.0}", val)
    }
}

/// Format a bit score matching BLAST reference style.
/// >= 99.9: integer (e.g. "198"), < 99.9: one decimal (e.g. "56.0").
fn format_bitscore(val: f64) -> String {
    if val > 99999.0 {
        format!("{:.3e}", val)
    } else if val > 99.9 {
        format!("{:.0}", val as i64)
    } else {
        format!("{:.1}", val)
    }
}

/// Write tabular output (outfmt 6) for a set of hits.
pub fn format_tabular<W: Write>(writer: &mut W, hits: &[TabularHit]) -> std::io::Result<()> {
    for hit in hits {
        writeln!(
            writer,
            "{}\t{}\t{:.3}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
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
            format_bitscore(hit.bit_score),
        )?;
    }
    Ok(())
}
