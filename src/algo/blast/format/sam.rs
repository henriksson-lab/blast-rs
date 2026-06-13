//! SAM output format (-outfmt 17).

use std::io::Write;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SamReference<'a> {
    pub name: &'a str,
    pub len: i32,
}

#[derive(Debug, Clone, Copy)]
pub struct BlastnSamRecord<'a> {
    pub subject_id: &'a str,
    pub query_id: &'a str,
    pub query_start: i32,
    pub query_end: i32,
    pub subject_start: i32,
    pub subject_end: i32,
    pub subject_len: i32,
    pub raw_score: i32,
    pub bit_score: f64,
    pub evalue: f64,
    pub pct_identity: f64,
    pub align_len: i32,
    pub qseq: Option<&'a [u8]>,
    pub sseq: Option<&'a [u8]>,
}

/// Write SAM header.
/// blast-rs: Convenience wrapper for local outfmt 17 header emission.
pub fn write_sam_header<W: Write>(writer: &mut W, db_name: &str) -> std::io::Result<()> {
    write_sam_header_with_references(writer, db_name, &[])
}

/// Write SAM header with reference sequence lines.
///
/// NCBI blastn outfmt 17 emits `@SQ` records for query/reference sequences
/// represented in the alignment stream. This helper keeps the old db-name
/// `@PG` behavior while exposing those reference lines to library callers.
/// blast-rs: Library SAM header helper shaped around NCBI outfmt 17 output.
pub fn write_sam_header_with_references<W: Write>(
    writer: &mut W,
    db_name: &str,
    references: &[SamReference<'_>],
) -> std::io::Result<()> {
    writeln!(writer, "@HD\tVN:1.2\tSO:coordinate\tGO:reference")?;
    for reference in references {
        writeln!(
            writer,
            "@SQ\tSN:{}\tLN:{}",
            sanitize_sam_field(reference.name),
            reference.len.max(0)
        )?;
    }
    // Sanitize db_name: tabs/newlines would corrupt SAM format
    let safe_name = sanitize_sam_field(db_name);
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

/// blast-rs: SAM header-field sanitizer for local outfmt 17 helpers.
pub fn sanitize_sam_field(value: &str) -> String {
    value
        .chars()
        .map(|c| {
            if c == '\t' || c == '\n' || c == '\r' {
                ' '
            } else {
                c
            }
        })
        .collect()
}

/// Write a SAM alignment record.
///
/// Optional tags:
/// - `AS:i:` — alignment score (standard)
/// - `NM:i:` — edit distance (standard)
/// - `XQ:i:` — query start position (BLAST-specific, X-namespace)
/// - `XR:i:` — query end position (BLAST-specific, X-namespace)
/// blast-rs: Library SAM record helper; not a direct NCBI C port.
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
/// blast-rs: Library SAM gapped-record helper; not a direct NCBI C port.
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

/// Write a BLASTN outfmt 17-style SAM alignment record.
///
/// NCBI treats the database/subject sequence as the SAM read and the query as
/// the reference. This differs from the simpler helper above, which writes the
/// query as QNAME for generic SAM-like output.
/// blast-rs: Library BLASTN outfmt 17 record helper; not a direct NCBI C port.
pub fn write_blastn_sam_record<W: Write>(
    writer: &mut W,
    record: &BlastnSamRecord<'_>,
) -> std::io::Result<()> {
    let flag = if record.subject_start > record.subject_end {
        16
    } else {
        0
    };
    let pos = record.query_start.min(record.query_end);
    let cigar = build_subject_as_read_cigar(
        record.qseq,
        record.sseq,
        record.align_len,
        record.subject_start,
        record.subject_end,
        record.subject_len,
    );
    let nm = sam_gap_count(record.qseq, record.sseq);
    let pi = sam_pairwise_identity(record.qseq, record.sseq, record.pct_identity);

    writeln!(
        writer,
        "{}\t{}\t{}\t{}\t255\t{}\t*\t0\t0\t*\t*\tAS:i:{}\tEV:f:{}\tNM:i:{}\tPI:f:{:.2}\tBS:f:{}",
        record.subject_id,
        flag,
        record.query_id,
        pos,
        cigar,
        record.raw_score,
        format_sam_float(record.evalue),
        nm,
        pi,
        format_sam_float(record.bit_score),
    )
}

/// Build a CIGAR string from aligned query and subject sequences.
/// M = match/mismatch, I = gap in subject/consumes query, D = gap in query/consumes subject.
// blast-rs: SAM CIGAR synthesis from BLAST traceback strings; this is native
// outfmt 17 glue, not a port of an isolated NCBI helper.
pub fn build_cigar(query_aln: &[u8], subject_aln: &[u8]) -> String {
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

/// blast-rs: CIGAR helper for BLASTN outfmt 17's subject-as-read orientation.
pub fn build_subject_as_read_cigar(
    query_aln: Option<&[u8]>,
    subject_aln: Option<&[u8]>,
    align_len: i32,
    subject_start: i32,
    subject_end: i32,
    subject_len: i32,
) -> String {
    let s_lo = subject_start.min(subject_end);
    let s_hi = subject_start.max(subject_end);
    let pre_clip = (s_lo - 1).max(0);
    let post_clip = (subject_len - s_hi).max(0);

    let body = if let (Some(qseq), Some(sseq)) = (query_aln, subject_aln) {
        let mut cigar = String::new();
        let mut current_op = '\0';
        let mut current_len = 0usize;
        for (&q, &s) in qseq.iter().zip(sseq.iter()) {
            let op = if q == b'-' {
                'I'
            } else if s == b'-' {
                'D'
            } else {
                'M'
            };
            if op == current_op {
                current_len += 1;
            } else {
                if current_len > 0 {
                    cigar.push_str(&current_len.to_string());
                    cigar.push(current_op);
                }
                current_op = op;
                current_len = 1;
            }
        }
        if current_len > 0 {
            cigar.push_str(&current_len.to_string());
            cigar.push(current_op);
        }
        if cigar.is_empty() {
            format!("{}M", align_len)
        } else {
            cigar
        }
    } else {
        format!("{}M", align_len)
    };

    let mut result = String::new();
    if pre_clip > 0 {
        result.push_str(&format!("{}H", pre_clip));
    }
    result.push_str(&body);
    if post_clip > 0 {
        result.push_str(&format!("{}H", post_clip));
    }
    result
}

/// blast-rs: Gap-count helper matching BLASTN SAM `NM` tag behavior.
pub fn sam_gap_count(query_aln: Option<&[u8]>, subject_aln: Option<&[u8]>) -> i32 {
    match (query_aln, subject_aln) {
        (Some(qseq), Some(sseq)) => qseq
            .iter()
            .chain(sseq.iter())
            .filter(|&&base| base == b'-')
            .count() as i32,
        _ => 0,
    }
}

/// blast-rs: Pairwise identity helper for BLASTN SAM `PI` tag formatting.
pub fn sam_pairwise_identity(
    query_aln: Option<&[u8]>,
    subject_aln: Option<&[u8]>,
    fallback: f64,
) -> f64 {
    let (Some(qseq), Some(sseq)) = (query_aln, subject_aln) else {
        return fallback;
    };
    let mut compared = 0usize;
    let mut identical = 0usize;
    for (&q, &s) in qseq.iter().zip(sseq.iter()) {
        if q != b'-' && s != b'-' {
            compared += 1;
            if q == s {
                identical += 1;
            }
        }
    }
    if compared == 0 {
        0.0
    } else {
        100.0 * identical as f64 / compared as f64
    }
}

/// blast-rs: NCBI-like compact floating-point formatter for SAM optional tags.
pub fn format_sam_float(value: f64) -> String {
    if value != 0.0 && value.abs() < 1.0e-180 {
        "0".to_string()
    } else {
        format!("{value:.6e}")
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
        assert_eq!(
            output.lines().next(),
            Some("@HD\tVN:1.2\tSO:coordinate\tGO:reference")
        );
        assert!(!output.contains("my\tdb")); // tab should be sanitized
        assert!(output.contains("my db path")); // replaced with spaces
    }

    #[test]
    fn test_sam_header_with_references_writes_sq_lines() {
        let mut buf = Vec::new();
        write_sam_header_with_references(
            &mut buf,
            "testdb",
            &[
                SamReference {
                    name: "query1",
                    len: 12,
                },
                SamReference {
                    name: "bad\tname",
                    len: -5,
                },
            ],
        )
        .unwrap();
        let output = String::from_utf8(buf).unwrap();
        let lines: Vec<_> = output.lines().collect();
        assert_eq!(lines[0], "@HD\tVN:1.2\tSO:coordinate\tGO:reference");
        assert_eq!(lines[1], "@SQ\tSN:query1\tLN:12");
        assert_eq!(lines[2], "@SQ\tSN:bad name\tLN:0");
        assert!(lines[3].starts_with("@PG\t"));
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
    fn test_blastn_sam_record_uses_subject_as_read_with_ncbi_tags() {
        let mut buf = Vec::new();
        let record = BlastnSamRecord {
            subject_id: "subject1",
            query_id: "query1",
            query_start: 10,
            query_end: 19,
            subject_start: 5,
            subject_end: 14,
            subject_len: 20,
            raw_score: 42,
            bit_score: 23.5,
            evalue: 2.0e-12,
            pct_identity: 0.0,
            align_len: 12,
            qseq: Some(b"ACGT--TTACGT"),
            sseq: Some(b"ACGTGGTT--GT"),
        };

        write_blastn_sam_record(&mut buf, &record).unwrap();

        let output = String::from_utf8(buf).unwrap();
        let fields: Vec<&str> = output.trim_end().split('\t').collect();
        assert_eq!(fields[0], "subject1");
        assert_eq!(fields[1], "0");
        assert_eq!(fields[2], "query1");
        assert_eq!(fields[3], "10");
        assert_eq!(fields[5], "4H4M2I2M2D2M6H");
        assert!(output.contains("\tAS:i:42\t"));
        assert!(output.contains("\tEV:f:2.000000e-12\t"));
        assert!(output.contains("\tNM:i:4\t"));
        assert!(output.contains("\tPI:f:100.00\t"));
        assert!(output.contains("\tBS:f:2.350000e1"));
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
