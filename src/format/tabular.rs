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
    pub query_len: i32,
    pub subject_len: i32,
    pub raw_score: i32,
    pub qseq: Option<String>,
    pub sseq: Option<String>,
    pub qframe: i32,
    pub sframe: i32,
    pub subject_taxids: Vec<i32>,
    pub subject_sci_name: String,
    pub subject_common_name: String,
    pub subject_blast_name: String,
    pub subject_kingdom: String,
}

/// Format an e-value matching BLAST tabular output (-outfmt 6).
/// Uses NStr::DoubleToString(evalue, 2, fDoubleScientific) which is
/// essentially %.2e for small values, transitioning to fixed for larger ones.
pub fn format_evalue(val: f64) -> String {
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
pub fn format_bitscore(val: f64) -> String {
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
        "qseqid" | "qacc" | "qaccver" => hit.query_id.clone(),
        "sseqid" | "sacc" | "saccver" => hit.subject_id.clone(),
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
        "score" => hit.raw_score.to_string(),
        "nident" => (hit.align_len - hit.mismatches).to_string(),
        "gaps" => hit.gap_opens.to_string(),
        "qlen" => hit.query_len.to_string(),
        "slen" => hit.subject_len.to_string(),
        "qcovs" => {
            if hit.query_len > 0 {
                let cov = 100.0 * (hit.query_end - hit.query_start + 1).abs() as f64 / hit.query_len as f64;
                format!("{:.0}", cov.min(100.0))
            } else { "0".to_string() }
        }
        "qcovhsp" => {
            if hit.query_len > 0 {
                let cov = 100.0 * hit.align_len as f64 / hit.query_len as f64;
                format!("{:.0}", cov.min(100.0))
            } else { "0".to_string() }
        }
        // Taxonomy ID fields — from .nto database file
        "staxid" => hit.subject_taxids.first().map(|t| t.to_string()).unwrap_or_else(|| "N/A".to_string()),
        "staxids" => if hit.subject_taxids.is_empty() { "N/A".to_string() }
            else { hit.subject_taxids.iter().map(|t| t.to_string()).collect::<Vec<_>>().join(";") },
        // Taxonomy name fields — from taxdb.bti/taxdb.btd
        "ssciname" | "sscinames" => if hit.subject_sci_name.is_empty() { "N/A".to_string() } else { hit.subject_sci_name.clone() },
        "scomname" | "scomnames" => if hit.subject_common_name.is_empty() { "N/A".to_string() } else { hit.subject_common_name.clone() },
        "sblastname" | "sblastnames" => if hit.subject_blast_name.is_empty() { "N/A".to_string() } else { hit.subject_blast_name.clone() },
        "sskingdom" | "sskingdoms" => if hit.subject_kingdom.is_empty() { "N/A".to_string() } else { hit.subject_kingdom.clone() },
        // Sequence hash
        "qseq" => hit.qseq.as_deref().unwrap_or("N/A").to_string(),
        "sseq" => hit.sseq.as_deref().unwrap_or("N/A").to_string(),
        // Frame fields
        "qframe" => hit.qframe.to_string(),
        "sframe" => hit.sframe.to_string(),
        "frames" => format!("{}/{}", hit.qframe, hit.sframe),
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

#[cfg(test)]
mod tests {
    use super::*;

    fn make_hit(qseq: Option<&str>, sseq: Option<&str>) -> TabularHit {
        TabularHit {
            query_id: "q1".to_string(),
            subject_id: "s1".to_string(),
            pct_identity: 95.0,
            align_len: 50,
            mismatches: 2,
            gap_opens: 1,
            query_start: 1,
            query_end: 50,
            subject_start: 10,
            subject_end: 59,
            evalue: 1e-20,
            bit_score: 80.0,
            query_len: 100,
            subject_len: 200,
            raw_score: 120,
            qseq: qseq.map(|s| s.to_string()),
            sseq: sseq.map(|s| s.to_string()),
            qframe: 1,
            sframe: 0,
            subject_taxids: vec![9606],
            subject_sci_name: "Homo sapiens".to_string(),
            subject_common_name: "human".to_string(),
            subject_blast_name: "primates".to_string(),
            subject_kingdom: "E".to_string(),
        }
    }

    #[test]
    fn test_qseq_sseq_returned_when_present() {
        let hit = make_hit(Some("ACGTACGT"), Some("ACGT-CGT"));
        assert_eq!(get_field(&hit, "qseq"), "ACGTACGT");
        assert_eq!(get_field(&hit, "sseq"), "ACGT-CGT");
    }

    #[test]
    fn test_qseq_sseq_na_when_none() {
        let hit = make_hit(None, None);
        assert_eq!(get_field(&hit, "qseq"), "N/A");
        assert_eq!(get_field(&hit, "sseq"), "N/A");
    }

    #[test]
    fn test_custom_format_includes_aligned_seqs() {
        let hit = make_hit(Some("AAACCC"), Some("AAA-CC"));
        let hits = vec![hit];
        let mut buf = Vec::new();
        format_tabular_custom(&mut buf, &hits, "qseqid qseq sseq").unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert_eq!(output.trim(), "q1\tAAACCC\tAAA-CC");
    }

    #[test]
    fn test_raw_score_field() {
        let hit = make_hit(None, None);
        assert_eq!(get_field(&hit, "score"), "120");
    }

    #[test]
    fn test_frame_fields() {
        let hit = make_hit(None, None);
        assert_eq!(get_field(&hit, "qframe"), "1");
        assert_eq!(get_field(&hit, "sframe"), "0");
        assert_eq!(get_field(&hit, "frames"), "1/0");
    }

    #[test]
    fn test_taxid_fields() {
        let hit = make_hit(None, None);
        assert_eq!(get_field(&hit, "staxid"), "9606");
        assert_eq!(get_field(&hit, "staxids"), "9606");
    }

    #[test]
    fn test_taxid_multiple() {
        let mut hit = make_hit(None, None);
        hit.subject_taxids = vec![9606, 9605];
        assert_eq!(get_field(&hit, "staxid"), "9606");
        assert_eq!(get_field(&hit, "staxids"), "9606;9605");
    }

    #[test]
    fn test_tax_name_fields() {
        let hit = make_hit(None, None);
        assert_eq!(get_field(&hit, "ssciname"), "Homo sapiens");
        assert_eq!(get_field(&hit, "scomname"), "human");
        assert_eq!(get_field(&hit, "sblastname"), "primates");
        assert_eq!(get_field(&hit, "sskingdom"), "E");
    }

    #[test]
    fn test_tax_name_empty_fallback() {
        let mut hit = make_hit(None, None);
        hit.subject_sci_name = String::new();
        assert_eq!(get_field(&hit, "ssciname"), "N/A");
    }

    #[test]
    fn test_taxid_empty() {
        let mut hit = make_hit(None, None);
        hit.subject_taxids = vec![];
        assert_eq!(get_field(&hit, "staxid"), "N/A");
        assert_eq!(get_field(&hit, "staxids"), "N/A");
    }
}
