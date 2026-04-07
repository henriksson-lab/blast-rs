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
        // Scientific notation with 2 decimal places, C-style zero-padded exponent
        let s = format!("{:.2e}", val);
        // Rust: "2.01e-4" → need "2.01e-04"
        // Fix: ensure exponent has at least 2 digits
        if let Some(e_pos) = s.find('e') {
            let (mantissa, exp_part) = s.split_at(e_pos);
            let exp_str = &exp_part[1..]; // skip 'e'
            let (sign, digits) = if exp_str.starts_with('-') {
                ("-", &exp_str[1..])
            } else if exp_str.starts_with('+') {
                ("", &exp_str[1..])
            } else {
                ("", exp_str)
            };
            if digits.len() < 2 {
                format!("{}e{}{:02}", mantissa, sign, digits.parse::<u32>().unwrap_or(0))
            } else {
                s
            }
        } else {
            s
        }
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

/// Get a field value by column name for custom tabular output.
pub fn get_field(hit: &TabularHit, column: &str) -> String {
    match column {
        "qseqid" | "qacc" => hit.query_id.clone(),
        "sseqid" | "sacc" => hit.subject_id.clone(),
        "pident" => format!("{:.3}", hit.pct_identity),
        "length" => hit.align_len.to_string(),
        "mismatch" => hit.mismatches.to_string(),
        "gapopen" => hit.gap_opens.to_string(),
        "qstart" => hit.query_start.to_string(),
        "qend" => hit.query_end.to_string(),
        "sstart" => hit.subject_start.to_string(),
        "send" => hit.subject_end.to_string(),
        "evalue" => format_evalue(hit.evalue),
        "bitscore" => format_bitscore(hit.bit_score),
        "score" => format!("{:.0}", hit.bit_score), // raw score approximation
        "nident" => (hit.align_len - hit.mismatches).to_string(),
        "gaps" => hit.gap_opens.to_string(),
        "qlen" => "0".to_string(), // not available in TabularHit
        "slen" => "0".to_string(),
        _ => "N/A".to_string(),
    }
}

/// Write tabular output with custom columns.
/// `columns` is a space-separated list of field names.
pub fn format_tabular_custom<W: Write>(
    writer: &mut W,
    hits: &[TabularHit],
    columns: &str,
) -> std::io::Result<()> {
    let cols: Vec<&str> = columns.split_whitespace().collect();
    for hit in hits {
        let fields: Vec<String> = cols.iter().map(|c| get_field(hit, c)).collect();
        writeln!(writer, "{}", fields.join("\t"))?;
    }
    Ok(())
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
