//! SAM output format (-outfmt 17).

use std::io::Write;

/// Write SAM header.
pub fn write_sam_header<W: Write>(writer: &mut W, _db_name: &str) -> std::io::Result<()> {
    writeln!(writer, "@HD\tVN:1.0\tSO:queryname")?;
    writeln!(writer, "@PG\tID:blast\tPN:blastn\tCL:blast")?;
    Ok(())
}

/// Write a SAM alignment record.
pub fn write_sam_record<W: Write>(
    writer: &mut W,
    query_id: &str,
    subject_id: &str,
    _query_start: i32,
    _query_end: i32,
    subject_start: i32,
    _subject_end: i32,
    score: i32,
    align_len: i32,
    num_ident: i32,
    is_minus: bool,
) -> std::io::Result<()> {
    let flag = if is_minus { 16 } else { 0 }; // 16 = reverse strand
    let mapq = 255; // mapping quality unavailable
    let cigar = format!("{}M", align_len); // simplified: all matches

    writeln!(writer, "{}\t{}\t{}\t{}\t{}\t{}\t*\t0\t0\t*\t*\tAS:i:{}\tNM:i:{}",
        query_id,
        flag,
        subject_id,
        subject_start,
        mapq,
        cigar,
        score,
        align_len - num_ident, // NM = edit distance
    )?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sam_record() {
        let mut buf = Vec::new();
        write_sam_record(&mut buf, "query1", "chr1", 1, 50, 100, 149, 100, 50, 50, false).unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert!(output.contains("query1"));
        assert!(output.contains("chr1"));
        assert!(output.contains("50M")); // CIGAR
        assert!(output.contains("AS:i:100")); // alignment score
    }

    #[test]
    fn test_sam_minus_strand() {
        let mut buf = Vec::new();
        write_sam_record(&mut buf, "q1", "s1", 1, 20, 100, 81, 40, 20, 18, true).unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert!(output.contains("\t16\t")); // reverse strand flag
    }
}
