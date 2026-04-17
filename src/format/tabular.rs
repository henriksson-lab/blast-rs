//! BLAST tabular output format (-outfmt 6/7).
//!
//! Default columns: qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore

use std::collections::HashMap;
use std::io::Write;

pub const DEFAULT_TABULAR_COLUMNS: &str =
    "qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore";

fn expand_column_token<'a>(token: &'a str, columns: &mut Vec<&'a str>) {
    if token == "std" {
        columns.extend(DEFAULT_TABULAR_COLUMNS.split_whitespace());
    } else {
        columns.push(token);
    }
}

fn normalize_column_tokens(cols: &mut Vec<&str>) {
    cols.retain(|col| field_display_name(col) != "unknown field");
    let mut seen = std::collections::HashSet::new();
    cols.retain(|col| seen.insert(*col));
    if cols.is_empty() {
        cols.extend(DEFAULT_TABULAR_COLUMNS.split_whitespace());
    }
}

pub fn expanded_column_tokens(columns: &str) -> Vec<&str> {
    let mut cols = Vec::new();
    for col in columns.split_whitespace() {
        if !col.starts_with("delim=") {
            expand_column_token(col, &mut cols);
        }
    }
    normalize_column_tokens(&mut cols);
    cols
}

/// A single alignment hit for tabular output.
#[derive(Clone)]
pub struct TabularHit {
    pub query_id: String,
    pub query_gi: Option<String>,
    pub query_acc: Option<String>,
    pub query_accver: Option<String>,
    pub subject_id: String,
    pub subject_gi: Option<String>,
    pub subject_acc: Option<String>,
    pub subject_accver: Option<String>,
    pub subject_title: String,
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
    pub num_ident: i32,
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
                format!(
                    "{}e{}{:02}",
                    mantissa,
                    sign,
                    digits.parse::<u32>().unwrap_or(0)
                )
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
    get_field_with_qcovs(hit, column, None)
}

pub fn field_display_name(column: &str) -> &'static str {
    match column {
        "qseqid" => "query id",
        "qgi" => "query gi",
        "qacc" => "query acc.",
        "qaccver" => "query acc.ver",
        "sseqid" => "subject id",
        "sgi" => "subject gi",
        "sallgi" => "all subject gis",
        "sacc" => "subject acc.",
        "saccver" => "subject acc.ver",
        "sallseqid" => "all subject ids",
        "sallacc" => "all subject acc.",
        "stitle" => "subject title",
        "salltitles" => "all subject title(s)",
        "pident" => "% identity",
        "length" => "alignment length",
        "mismatch" => "mismatches",
        "gapopen" => "gap opens",
        "qstart" => "q. start",
        "qend" => "q. end",
        "sstart" => "s. start",
        "send" => "s. end",
        "evalue" => "evalue",
        "bitscore" => "bit score",
        "score" => "score",
        "nident" => "identical",
        "positive" => "positives",
        "gaps" => "gaps",
        "ppos" => "% positives",
        "qlen" => "query length",
        "slen" => "subject length",
        "qcovs" => "% query coverage per subject",
        "qcovhsp" => "% query coverage per hsp",
        "qcovus" => "% query coverage per unique subject",
        "staxid" => "subject tax id",
        "staxids" => "unique subject tax ids",
        "ssciname" | "sscinames" => "subject sci name",
        "scomname" | "scomnames" => "subject common name",
        "sblastname" | "sblastnames" => "subject blast name",
        "sskingdom" | "sskingdoms" => "subject super kingdom",
        "qseq" => "aligned part of query sequence",
        "sseq" => "aligned part of subject sequence",
        "btop" => "BTOP",
        "sstrand" => "subject strand",
        "qframe" => "query frame",
        "sframe" => "subject frame",
        "frames" => "query/sbjct frames",
        _ => "unknown field",
    }
}

fn get_field_with_qcovs(hit: &TabularHit, column: &str, qcovs: Option<i32>) -> String {
    match column {
        "qseqid" => hit.query_id.clone(),
        "qacc" => hit.query_acc.as_ref().unwrap_or(&hit.query_id).clone(),
        "qaccver" => hit
            .query_accver
            .as_ref()
            .or(hit.query_acc.as_ref())
            .unwrap_or(&hit.query_id)
            .clone(),
        "qgi" => hit.query_gi.as_deref().unwrap_or("0").to_string(),
        "sseqid" | "sallseqid" => hit.subject_id.clone(),
        "sacc" | "sallacc" => hit.subject_acc.as_ref().unwrap_or(&hit.subject_id).clone(),
        "saccver" => hit
            .subject_accver
            .as_ref()
            .or(hit.subject_acc.as_ref())
            .unwrap_or(&hit.subject_id)
            .clone(),
        "sgi" | "sallgi" => hit.subject_gi.as_deref().unwrap_or("0").to_string(),
        "stitle" | "salltitles" => {
            if hit.subject_title.is_empty() {
                "N/A".to_string()
            } else {
                hit.subject_title.clone()
            }
        }
        "pident" => format!("{:.3}", hit.pct_identity),
        "length" => hit.align_len.to_string(),
        "mismatch" => {
            if let (Some(qseq), Some(sseq)) = (hit.qseq.as_deref(), hit.sseq.as_deref()) {
                qseq.bytes()
                    .zip(sseq.bytes())
                    .filter(|&(q, s)| q != b'-' && s != b'-' && q != s)
                    .count()
                    .to_string()
            } else {
                hit.mismatches.to_string()
            }
        }
        "gapopen" => hit.gap_opens.to_string(),
        "qstart" => hit.query_start.to_string(),
        "qend" => hit.query_end.to_string(),
        "sstart" => hit.subject_start.to_string(),
        "send" => hit.subject_end.to_string(),
        "evalue" => format_evalue(hit.evalue),
        "bitscore" => format_bitscore(hit.bit_score),
        "score" => hit.raw_score.to_string(),
        "nident" => hit.num_ident.to_string(),
        "positive" => hit.num_ident.to_string(),
        "gaps" => {
            if let (Some(qseq), Some(sseq)) = (hit.qseq.as_deref(), hit.sseq.as_deref()) {
                qseq.bytes()
                    .chain(sseq.bytes())
                    .filter(|&base| base == b'-')
                    .count()
                    .to_string()
            } else {
                (hit.align_len - hit.num_ident - hit.mismatches)
                    .max(0)
                    .to_string()
            }
        }
        "ppos" => {
            if hit.align_len > 0 {
                format!("{:.2}", 100.0 * hit.num_ident as f64 / hit.align_len as f64)
            } else {
                "0.00".to_string()
            }
        }
        "qlen" => hit.query_len.to_string(),
        "slen" => hit.subject_len.to_string(),
        "qcovs" | "qcovus" => {
            if let Some(qcovs) = qcovs {
                qcovs.to_string()
            } else if hit.query_len > 0 {
                let cov = 100.0 * (hit.query_end - hit.query_start + 1).abs() as f64
                    / hit.query_len as f64;
                format_query_coverage(cov)
            } else {
                "0".to_string()
            }
        }
        "qcovhsp" => {
            if hit.query_len > 0 {
                let query_span = (hit.query_end - hit.query_start).abs() + 1;
                let cov = 100.0 * query_span as f64 / hit.query_len as f64;
                format_query_coverage(cov)
            } else {
                "0".to_string()
            }
        }
        // Taxonomy ID fields — from .nto database file
        "staxid" => hit
            .subject_taxids
            .first()
            .map(|t| t.to_string())
            .unwrap_or_else(|| "N/A".to_string()),
        "staxids" => {
            if hit.subject_taxids.is_empty() {
                "N/A".to_string()
            } else {
                hit.subject_taxids
                    .iter()
                    .map(|t| t.to_string())
                    .collect::<Vec<_>>()
                    .join(";")
            }
        }
        // Taxonomy name fields — from taxdb.bti/taxdb.btd
        "ssciname" | "sscinames" => {
            if hit.subject_sci_name.is_empty() {
                "N/A".to_string()
            } else {
                hit.subject_sci_name.clone()
            }
        }
        "scomname" | "scomnames" => {
            if hit.subject_common_name.is_empty() {
                "N/A".to_string()
            } else {
                hit.subject_common_name.clone()
            }
        }
        "sblastname" | "sblastnames" => {
            if hit.subject_blast_name.is_empty() {
                "N/A".to_string()
            } else {
                hit.subject_blast_name.clone()
            }
        }
        "sskingdom" | "sskingdoms" => {
            if hit.subject_kingdom.is_empty() {
                "N/A".to_string()
            } else {
                hit.subject_kingdom.clone()
            }
        }
        // Sequence hash
        "qseq" => hit.qseq.as_deref().unwrap_or("N/A").to_string(),
        "sseq" => hit.sseq.as_deref().unwrap_or("N/A").to_string(),
        "btop" => format_btop(hit),
        "sstrand" => {
            if hit.subject_start <= hit.subject_end {
                "plus".to_string()
            } else {
                "minus".to_string()
            }
        }
        // Frame fields
        "qframe" => hit.qframe.to_string(),
        "sframe" => hit.sframe.to_string(),
        "frames" => format!("{}/{}", hit.qframe, hit.sframe),
        _ => "N/A".to_string(),
    }
}

fn format_btop(hit: &TabularHit) -> String {
    let (Some(qseq), Some(sseq)) = (hit.qseq.as_deref(), hit.sseq.as_deref()) else {
        return "N/A".to_string();
    };

    let mut out = String::new();
    let mut matches = 0;
    for (q, s) in qseq.bytes().zip(sseq.bytes()) {
        if q == s && q != b'-' {
            matches += 1;
        } else {
            if matches > 0 {
                out.push_str(&matches.to_string());
                matches = 0;
            }
            out.push(q as char);
            out.push(s as char);
        }
    }
    if matches > 0 {
        out.push_str(&matches.to_string());
    }
    out
}

fn compute_qcovs_by_query_subject(hits: &[TabularHit]) -> HashMap<(&str, &str), i32> {
    let mut intervals: HashMap<(&str, &str), (i32, Vec<(i32, i32)>)> = HashMap::new();
    for hit in hits {
        let start = hit.query_start.min(hit.query_end);
        let end = hit.query_start.max(hit.query_end);
        intervals
            .entry((hit.query_id.as_str(), hit.subject_id.as_str()))
            .or_insert_with(|| (hit.query_len, Vec::new()))
            .1
            .push((start, end));
    }

    let mut out = HashMap::new();
    for (key, (query_len, mut ranges)) in intervals {
        if query_len <= 0 || ranges.is_empty() {
            out.insert(key, 0);
            continue;
        }
        ranges.sort_unstable_by_key(|&(start, end)| (start, end));
        let mut covered = 0;
        let mut cur = ranges[0];
        for &(start, end) in &ranges[1..] {
            if start <= cur.1 + 1 {
                cur.1 = cur.1.max(end);
            } else {
                covered += cur.1 - cur.0 + 1;
                cur = (start, end);
            }
        }
        covered += cur.1 - cur.0 + 1;
        let cov = rounded_query_coverage(100.0 * covered as f64 / query_len as f64);
        out.insert(key, cov);
    }
    out
}

/// Write tabular output with custom columns.
/// `columns` is a space-separated list of field names.
pub fn format_tabular_custom<W: Write>(
    writer: &mut W,
    hits: &[TabularHit],
    columns: &str,
) -> std::io::Result<()> {
    format_tabular_custom_with_delimiter(writer, hits, columns, "\t")
}

pub fn format_tabular_custom_with_delimiter<W: Write>(
    writer: &mut W,
    hits: &[TabularHit],
    columns: &str,
    default_delimiter: &str,
) -> std::io::Result<()> {
    let mut delimiter = default_delimiter;
    let mut cols = Vec::new();
    for col in columns.split_whitespace() {
        if let Some(value) = col.strip_prefix("delim=") {
            delimiter = match value {
                r"\t" | "tab" => "\t",
                "space" => " ",
                "" => default_delimiter,
                value => value,
            };
        } else {
            expand_column_token(col, &mut cols);
        }
    }
    normalize_column_tokens(&mut cols);
    let qcovs_by_subject = if cols.iter().any(|&c| c == "qcovs" || c == "qcovus") {
        compute_qcovs_by_query_subject(hits)
    } else {
        HashMap::new()
    };
    for hit in hits {
        let qcovs = qcovs_by_subject
            .get(&(hit.query_id.as_str(), hit.subject_id.as_str()))
            .copied();
        let fields: Vec<String> = cols
            .iter()
            .map(|c| get_field_with_qcovs(hit, c, qcovs))
            .collect();
        writeln!(writer, "{}", fields.join(delimiter))?;
    }
    Ok(())
}

/// Write tabular output (outfmt 6) for a set of hits.
pub fn format_tabular<W: Write>(writer: &mut W, hits: &[TabularHit]) -> std::io::Result<()> {
    format_tabular_custom_with_delimiter(writer, hits, DEFAULT_TABULAR_COLUMNS, "\t")
}

fn format_query_coverage(cov: f64) -> String {
    rounded_query_coverage(cov).to_string()
}

fn rounded_query_coverage(cov: f64) -> i32 {
    ((cov.min(100.0) + 0.5).floor() as i32).min(100)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_hit(qseq: Option<&str>, sseq: Option<&str>) -> TabularHit {
        TabularHit {
            query_id: "q1".to_string(),
            query_gi: None,
            query_acc: None,
            query_accver: None,
            subject_id: "s1".to_string(),
            subject_gi: None,
            subject_acc: None,
            subject_accver: None,
            subject_title: "s1 synthetic title".to_string(),
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
            num_ident: 47,
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
    fn test_subject_title_fields() {
        let mut hit = make_hit(None, None);
        assert_eq!(get_field(&hit, "stitle"), "s1 synthetic title");
        assert_eq!(get_field(&hit, "salltitles"), "s1 synthetic title");
        hit.subject_title.clear();
        assert_eq!(get_field(&hit, "stitle"), "N/A");
        assert_eq!(get_field(&hit, "salltitles"), "N/A");
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

    #[test]
    fn test_tabular_all_standard_fields() {
        let hit = make_hit(None, None);
        let hits = vec![hit];
        let mut buf = Vec::new();
        format_tabular(&mut buf, &hits).unwrap();
        let output = String::from_utf8(buf).unwrap();
        let fields: Vec<&str> = output.trim().split('\t').collect();
        assert_eq!(
            fields.len(),
            12,
            "standard tabular format must have 12 fields"
        );
        assert_eq!(fields[0], "q1"); // qseqid
        assert_eq!(fields[1], "s1"); // sseqid
        assert_eq!(fields[2], "95.000"); // pident
        assert_eq!(fields[3], "50"); // length
        assert_eq!(fields[4], "2"); // mismatch
        assert_eq!(fields[5], "1"); // gapopen
        assert_eq!(fields[6], "1"); // qstart
        assert_eq!(fields[7], "50"); // qend
        assert_eq!(fields[8], "10"); // sstart
        assert_eq!(fields[9], "59"); // send
                                     // evalue and bitscore are formatted strings
        assert!(!fields[10].is_empty(), "evalue field must be present");
        assert!(!fields[11].is_empty(), "bitscore field must be present");
    }

    #[test]
    fn test_tabular_custom_columns() {
        let hit = make_hit(Some("ACGTACGT"), Some("ACGT-CGT"));
        let hits = vec![hit];
        let mut buf = Vec::new();
        format_tabular_custom(&mut buf, &hits, "qseq sseq qframe sframe score staxid").unwrap();
        let output = String::from_utf8(buf).unwrap();
        let fields: Vec<&str> = output.trim().split('\t').collect();
        assert_eq!(fields.len(), 6);
        assert_eq!(fields[0], "ACGTACGT"); // qseq
        assert_eq!(fields[1], "ACGT-CGT"); // sseq
        assert_eq!(fields[2], "1"); // qframe
        assert_eq!(fields[3], "0"); // sframe
        assert_eq!(fields[4], "120"); // score
        assert_eq!(fields[5], "9606"); // staxid
    }

    #[test]
    fn test_query_coverage_rounds_half_up() {
        let mut hit = make_hit(None, None);
        hit.query_len = 32;
        hit.query_start = 1;
        hit.query_end = 20;
        assert_eq!(get_field(&hit, "qcovhsp"), "63");
    }

    #[test]
    fn test_tabular_multiple_hits() {
        let hit1 = make_hit(None, None);
        let mut hit2 = make_hit(None, None);
        hit2.query_id = "q2".to_string();
        hit2.subject_id = "s2".to_string();
        let mut hit3 = make_hit(None, None);
        hit3.query_id = "q3".to_string();
        hit3.subject_id = "s3".to_string();
        let hits = vec![hit1, hit2, hit3];
        let mut buf = Vec::new();
        format_tabular(&mut buf, &hits).unwrap();
        let output = String::from_utf8(buf).unwrap();
        let lines: Vec<&str> = output.trim().split('\n').collect();
        assert_eq!(lines.len(), 3, "3 hits should produce 3 lines");
        assert!(lines[0].starts_with("q1\ts1\t"));
        assert!(lines[1].starts_with("q2\ts2\t"));
        assert!(lines[2].starts_with("q3\ts3\t"));
    }

    #[test]
    fn test_tabular_evalue_formatting() {
        // Very small evalue: scientific notation
        assert_eq!(format_evalue(1e-20), "1.00e-20");
        // Tiny evalue below threshold: "0.0"
        assert_eq!(format_evalue(0.0), "0.0");
        // Medium small: fixed with 3 decimals
        assert_eq!(format_evalue(0.005), "0.005");
        // Near 1: fixed with 2 decimals
        assert_eq!(format_evalue(0.5), "0.50");
        // Moderate: fixed with 1 decimal
        assert_eq!(format_evalue(5.0), "5.0");
        // Large: integer
        assert_eq!(format_evalue(100.0), "100");
        // Single-digit exponent should be zero-padded
        assert_eq!(format_evalue(1e-5), "1.00e-05");
    }

    #[test]
    fn test_tabular_identity_calculation() {
        let mut hit = make_hit(None, None);
        // 47 identities out of 50 alignment length = 94.0%
        hit.num_ident = 47;
        hit.align_len = 50;
        hit.pct_identity = 100.0 * 47.0 / 50.0; // 94.0
        assert_eq!(get_field(&hit, "pident"), "94.000");

        // Perfect identity
        hit.pct_identity = 100.0;
        assert_eq!(get_field(&hit, "pident"), "100.000");

        // Low identity
        hit.pct_identity = 33.333;
        assert_eq!(get_field(&hit, "pident"), "33.333");
    }
}
