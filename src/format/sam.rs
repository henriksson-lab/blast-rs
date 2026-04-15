//! SAM output format (-outfmt 17).

use std::io::Write;

/// Write SAM header.
pub fn write_sam_header<W: Write>(writer: &mut W, db_name: &str) -> std::io::Result<()> {
    writeln!(writer, "@HD\tVN:1.0\tSO:queryname")?;
    // Sanitize db_name: tabs/newlines would corrupt SAM format
    let safe_name: String = db_name
        .chars()
        .map(|c| {
            if c == '\t' || c == '\n' || c == '\r' {
                ' '
            } else {
                c
            }
        })
        .collect();
    if !safe_name.is_empty() {
        writeln!(
            writer,
            "@PG\tID:blast\tPN:blastn\tCL:blast\tDS:{}",
            safe_name
        )?;
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
    let cigar = format!("{}M", align_len); // ungapped default

    // SAM POS is the 1-based leftmost mapped position on the forward reference strand.
    // BLAST reports subject_start > subject_end for minus-strand hits, so min() gives
    // the leftmost position regardless of strand orientation.
    let pos = subject_start.min(subject_end);

    let edit_distance = align_len - num_ident;

    writeln!(
        writer,
        "{}\t{}\t{}\t{}\t{}\t{}\t*\t0\t0\t*\t*\tAS:i:{}\tNM:i:{}\tXQ:i:{}\tXR:i:{}",
        query_id, flag, subject_id, pos, mapq, cigar, score, edit_distance, query_start, query_end,
    )?;
    Ok(())
}

/// Write a SAM alignment record with proper CIGAR from aligned sequences.
///
/// `query_aln` and `subject_aln` are gap-containing alignment strings (with '-' for gaps).
/// If empty, falls back to `{align_len}M`.
pub fn write_sam_record_gapped<W: Write>(
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
    query_aln: &[u8],
    subject_aln: &[u8],
) -> std::io::Result<()> {
    let flag = if is_minus { 16 } else { 0 };
    let mapq = 255;

    let cigar = if !query_aln.is_empty() && !subject_aln.is_empty() {
        build_cigar(query_aln, subject_aln)
    } else {
        format!("{}M", align_len)
    };

    let pos = subject_start.min(subject_end);
    let edit_distance = align_len - num_ident;

    writeln!(
        writer,
        "{}\t{}\t{}\t{}\t{}\t{}\t*\t0\t0\t*\t*\tAS:i:{}\tNM:i:{}\tXQ:i:{}\tXR:i:{}",
        query_id, flag, subject_id, pos, mapq, cigar, score, edit_distance, query_start, query_end,
    )?;
    Ok(())
}

/// Build a CIGAR string from aligned query and subject sequences.
/// M = match/mismatch, I = insertion in query (gap in subject), D = deletion in query (gap in subject).
fn build_cigar(query_aln: &[u8], subject_aln: &[u8]) -> String {
    let mut cigar = String::new();
    let mut current_op = ' ';
    let mut current_len = 0u32;

    for (&q, &s) in query_aln.iter().zip(subject_aln.iter()) {
        let op = if q == b'-' {
            'D' // deletion in query = gap in query = consuming subject
        } else if s == b'-' {
            'I' // insertion in query = gap in subject = consuming query
        } else {
            'M' // match or mismatch
        };

        if op == current_op {
            current_len += 1;
        } else {
            if current_len > 0 {
                cigar.push_str(&format!("{}{}", current_len, current_op));
            }
            current_op = op;
            current_len = 1;
        }
    }
    if current_len > 0 {
        cigar.push_str(&format!("{}{}", current_len, current_op));
    }

    if cigar.is_empty() {
        format!("{}M", query_aln.len())
    } else {
        cigar
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sam_record_forward() {
        let mut buf = Vec::new();
        write_sam_record(
            &mut buf, "query1", "chr1", 1, 50, 100, 149, 100, 50, 50, false,
        )
        .unwrap();
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
    fn test_sam_gapped_cigar() {
        let mut buf = Vec::new();
        // Query:   ACGT--ACGT
        // Subject: ACGTTTACGT
        let qaln = b"ACGT--ACGT";
        let saln = b"ACGTTTACGT";
        write_sam_record_gapped(
            &mut buf, "q1", "s1", 1, 8, 1, 10, 50, 10, 8, false, qaln, saln,
        )
        .unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert!(
            output.contains("4M2D4M"),
            "CIGAR should be 4M2D4M, got: {}",
            output
        );
    }

    #[test]
    fn test_sam_insertion_cigar() {
        let mut buf = Vec::new();
        // Query:   ACGTTTACGT
        // Subject: ACGT--ACGT
        let qaln = b"ACGTTTACGT";
        let saln = b"ACGT--ACGT";
        write_sam_record_gapped(
            &mut buf, "q1", "s1", 1, 10, 1, 8, 50, 10, 8, false, qaln, saln,
        )
        .unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert!(
            output.contains("4M2I4M"),
            "CIGAR should be 4M2I4M, got: {}",
            output
        );
    }

    #[test]
    fn test_sam_header_sanitize() {
        let mut buf = Vec::new();
        write_sam_header(&mut buf, "my\tdb\npath").unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert!(!output.contains("my\tdb")); // tab should be sanitized
        assert!(output.contains("my db path")); // replaced with spaces
    }

    #[test]
    fn test_sam_cigar_with_gaps() {
        // Query:   ACGT--TTACGT
        // Subject: ACGTGGTT--GT
        // Expected CIGAR: 4M2D2M2I2M
        let qaln = b"ACGT--TTACGT";
        let saln = b"ACGTGGTT--GT";
        let cigar = build_cigar(qaln, saln);
        assert_eq!(
            cigar, "4M2D2M2I2M",
            "CIGAR should encode both insertions and deletions"
        );

        // Verify it also works through write_sam_record_gapped
        let mut buf = Vec::new();
        write_sam_record_gapped(
            &mut buf, "q1", "s1", 1, 10, 1, 10, 50, 12, 8, false, qaln, saln,
        )
        .unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert!(
            output.contains("4M2D2M2I2M"),
            "SAM record should contain correct CIGAR"
        );
    }

    #[test]
    fn test_sam_multiple_records() {
        let mut buf = Vec::new();
        write_sam_header(&mut buf, "testdb").unwrap();
        // Write 3 records
        write_sam_record(
            &mut buf, "query1", "chr1", 1, 50, 100, 149, 100, 50, 48, false,
        )
        .unwrap();
        write_sam_record(
            &mut buf, "query1", "chr2", 1, 50, 200, 249, 90, 50, 45, false,
        )
        .unwrap();
        write_sam_record(
            &mut buf, "query1", "chr3", 1, 50, 300, 349, 80, 50, 42, false,
        )
        .unwrap();
        let output = String::from_utf8(buf).unwrap();
        // Count non-header lines (those not starting with @)
        let records: Vec<&str> = output.lines().filter(|l| !l.starts_with('@')).collect();
        assert_eq!(records.len(), 3, "should have 3 SAM alignment records");
        // Verify each record references the correct subject
        assert!(records[0].contains("chr1"));
        assert!(records[1].contains("chr2"));
        assert!(records[2].contains("chr3"));
    }

    #[test]
    fn test_sam_unmapped() {
        // When a query has no hits, no SAM alignment records should be produced.
        // Only the header is written.
        let mut buf = Vec::new();
        write_sam_header(&mut buf, "testdb").unwrap();
        // Do NOT write any records (no hits found)
        let output = String::from_utf8(buf).unwrap();
        let records: Vec<&str> = output.lines().filter(|l| !l.starts_with('@')).collect();
        assert_eq!(
            records.len(),
            0,
            "query with no hits should produce no alignment records"
        );
        // Header should still be present
        assert!(output.contains("@HD"));
        assert!(output.contains("@PG"));
    }
}
