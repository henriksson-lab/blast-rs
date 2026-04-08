//! SAM output format (-outfmt 17).

use std::io::Write;

/// Write SAM header.
pub fn write_sam_header<W: Write>(writer: &mut W, db_name: &str) -> std::io::Result<()> {
    writeln!(writer, "@HD\tVN:1.0\tSO:queryname")?;
    // Sanitize db_name: tabs/newlines would corrupt SAM format
    let safe_name: String = db_name.chars()
        .map(|c| if c == '\t' || c == '\n' || c == '\r' { ' ' } else { c })
        .collect();
    if !safe_name.is_empty() {
        writeln!(writer, "@PG\tID:blast\tPN:blastn\tCL:blast\tDS:{}", safe_name)?;
    } else {
        writeln!(writer, "@PG\tID:blast\tPN:blastn\tCL:blast")?;
    }
    Ok(())
}

/// Write a SAM alignment record.
///
/// Optional tags:
/// - `AS:i:` — alignment score (standard)
/// - `NM:i:` — edit distance (standard)
/// - `XQ:i:` — query start position (BLAST-specific, X-namespace)
/// - `XR:i:` — query end position (BLAST-specific, X-namespace)
pub fn write_sam_record<W: Write>(
    writer: &mut W,
    query_id: &str,
    subject_id: &str,
    query_start: i32,
    query_end: i32,
    subject_start: i32,
    subject_end: i32,
    score: i32,
    align_len: i32,
    num_ident: i32,
    is_minus: bool,
) -> std::io::Result<()> {
    let flag = if is_minus { 16 } else { 0 }; // 16 = reverse strand
    let mapq = 255; // mapping quality unavailable
    let cigar = format!("{}M", align_len);

    // SAM POS is the 1-based leftmost mapped position on the forward reference strand.
    // BLAST reports subject_start > subject_end for minus-strand hits, so min() gives
    // the leftmost position regardless of strand orientation.
    let pos = subject_start.min(subject_end);

    let edit_distance = align_len - num_ident;

    writeln!(writer, "{}\t{}\t{}\t{}\t{}\t{}\t*\t0\t0\t*\t*\tAS:i:{}\tNM:i:{}\tXQ:i:{}\tXR:i:{}",
        query_id,
        flag,
        subject_id,
        pos,
        mapq,
        cigar,
        score,
        edit_distance,
        query_start,
        query_end,
    )?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sam_record_forward() {
        let mut buf = Vec::new();
        write_sam_record(&mut buf, "query1", "chr1", 1, 50, 100, 149, 100, 50, 50, false).unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert!(output.contains("query1"));
        assert!(output.contains("chr1"));
        assert!(output.contains("50M"));
        assert!(output.contains("AS:i:100"));
        assert!(output.contains("NM:i:0"));
        assert!(output.contains("XQ:i:1"));
        assert!(output.contains("XR:i:50"));
        assert!(output.contains("\t100\t")); // POS = min(100, 149)
    }

    #[test]
    fn test_sam_record_reverse() {
        let mut buf = Vec::new();
        write_sam_record(&mut buf, "q1", "s1", 1, 20, 100, 81, 40, 20, 18, true).unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert!(output.contains("\t16\t")); // reverse strand flag
        assert!(output.contains("\t81\t")); // POS = min(100, 81)
    }

    #[test]
    fn test_sam_header_sanitize() {
        let mut buf = Vec::new();
        write_sam_header(&mut buf, "my\tdb\npath").unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert!(!output.contains("my\tdb")); // tab should be sanitized
        assert!(output.contains("my db path")); // replaced with spaces
    }
}
