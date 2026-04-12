//! Pairwise text output format (-outfmt 0).
//! The traditional BLAST output with alignment visualization.

use std::io::Write;

/// Format a single alignment in pairwise text format.
pub fn format_pairwise_alignment<W: Write>(
    writer: &mut W,
    _query_id: &str,
    subject_id: &str,
    query_seq: &[u8],   // BLASTNA encoded
    subject_seq: &[u8], // BLASTNA encoded
    q_start: i32,
    q_end: i32,
    s_start: i32,
    s_end: i32,
    score: i32,
    bit_score: f64,
    evalue: f64,
    num_ident: i32,
    align_len: i32,
    gap_opens: i32,
) -> std::io::Result<()> {
    let pident = if align_len > 0 {
        100.0 * num_ident as f64 / align_len as f64
    } else {
        0.0
    };

    writeln!(writer, "> {}", subject_id)?;
    writeln!(writer, "")?;
    writeln!(writer, " Score = {:.1} bits ({}),  Expect = {:.2e}",
        bit_score, score, evalue)?;
    writeln!(writer, " Identities = {}/{} ({}%),  Gaps = {}/{} ({}%)",
        num_ident, align_len,
        (pident + 0.5) as i32,
        gap_opens, align_len,
        if align_len > 0 { (100.0 * gap_opens as f64 / align_len as f64 + 0.5) as i32 } else { 0 })?;
    writeln!(writer, " Strand=Plus/{}", if s_start <= s_end { "Plus" } else { "Minus" })?;
    writeln!(writer, "")?;

    // Display alignment in 60-character lines
    let blastna_to_char = |b: u8| -> char {
        match b {
            0 => 'A', 1 => 'C', 2 => 'G', 3 => 'T',
            4 => 'R', 5 => 'Y', 6 => 'M', 7 => 'K',
            8 => 'W', 9 => 'S', 10 => 'B', 11 => 'D',
            12 => 'H', 13 => 'V', 14 => 'N', _ => '-',
        }
    };

    let line_width = 60;
    let q_len = (q_end - q_start) as usize;
    let mut pos = 0;
    let mut qi = q_start;
    let mut si = s_start;
    let s_dir: i32 = if s_start <= s_end { 1 } else { -1 };

    while pos < q_len {
        let chunk = line_width.min(q_len - pos);

        // Query line
        write!(writer, "Query  {:<6} ", qi)?;
        for i in 0..chunk {
            if pos + i < query_seq.len() {
                write!(writer, "{}", blastna_to_char(query_seq[pos + i]))?;
            }
        }
        writeln!(writer, "  {}", qi + chunk as i32 - 1)?;

        // Match line
        write!(writer, "              ")?;
        for i in 0..chunk {
            if pos + i < query_seq.len() && pos + i < subject_seq.len()
                && query_seq[pos + i] == subject_seq[pos + i] {
                write!(writer, "|")?;
            } else {
                write!(writer, " ")?;
            }
        }
        writeln!(writer)?;

        // Subject line
        write!(writer, "Sbjct  {:<6} ", si)?;
        for i in 0..chunk {
            if pos + i < subject_seq.len() {
                write!(writer, "{}", blastna_to_char(subject_seq[pos + i]))?;
            }
        }
        let s_end_pos = si + s_dir * (chunk as i32 - 1);
        writeln!(writer, "  {}", s_end_pos)?;
        writeln!(writer)?;

        pos += chunk;
        qi += chunk as i32;
        si += s_dir * chunk as i32;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pairwise_format() {
        let mut buf = Vec::new();
        let query = vec![0u8, 1, 2, 3, 0, 1, 2, 3]; // ACGTACGT
        let subject = vec![0u8, 1, 2, 3, 0, 1, 2, 3];
        format_pairwise_alignment(
            &mut buf, "query1", "subject1",
            &query, &subject,
            1, 8, 1, 8,
            16, 30.0, 1e-5, 8, 8, 0,
        ).unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert!(output.contains("> subject1"));
        assert!(output.contains("Score"));
        assert!(output.contains("Query"));
        assert!(output.contains("Sbjct"));
    }

    #[test]
    fn test_pairwise_with_gaps() {
        // Query:   ACGT-ACGT (gap at position 4 in query => '-' = code 15 or 0xFF)
        // Subject: ACGTGACGT
        // BLASTNA: A=0, C=1, G=2, T=3, gap=15
        let query   = vec![0u8, 1, 2, 3, 15, 0, 1, 2, 3];
        let subject = vec![0u8, 1, 2, 3,  2, 0, 1, 2, 3];
        let mut buf = Vec::new();
        format_pairwise_alignment(
            &mut buf, "query1", "subject1",
            &query, &subject,
            1, 9, 1, 9,
            14, 28.0, 1e-3, 8, 9, 1,
        ).unwrap();
        let output = String::from_utf8(buf).unwrap();
        // The gap character should appear as '-' in query line
        assert!(output.contains("-"), "query line should contain a gap character '-'");
        // Should show identities and gaps info
        assert!(output.contains("Gaps = 1/9"));
        // Match line should have pipe where bases match
        assert!(output.contains("|"), "match line should contain pipe characters");
    }

    #[test]
    fn test_pairwise_minus_strand() {
        let query   = vec![0u8, 1, 2, 3, 0, 1, 2, 3]; // ACGTACGT
        let subject = vec![0u8, 1, 2, 3, 0, 1, 2, 3];
        let mut buf = Vec::new();
        // Minus strand: s_start > s_end
        format_pairwise_alignment(
            &mut buf, "query1", "subject1",
            &query, &subject,
            1, 8, 100, 93,   // subject coords decrease for minus strand
            16, 30.0, 1e-5, 8, 8, 0,
        ).unwrap();
        let output = String::from_utf8(buf).unwrap();
        // Should indicate Minus strand
        assert!(output.contains("Strand=Plus/Minus"), "output should indicate minus strand");
        // Subject should start at 100 and decrease
        assert!(output.contains("Sbjct  100"), "subject line should start at position 100");
    }

    #[test]
    fn test_pairwise_long_alignment() {
        // Create an alignment longer than 60 characters to test wrapping
        let len = 130;
        let query: Vec<u8> = (0..len).map(|i| (i % 4) as u8).collect();   // repeating ACGT
        let subject: Vec<u8> = (0..len).map(|i| (i % 4) as u8).collect();
        let mut buf = Vec::new();
        format_pairwise_alignment(
            &mut buf, "query1", "subject1",
            &query, &subject,
            1, len as i32, 1, len as i32,
            260, 480.0, 1e-100, len as i32, len as i32, 0,
        ).unwrap();
        let output = String::from_utf8(buf).unwrap();
        // Count "Query" lines — should be ceil(130/60) = 3 blocks
        let query_lines: Vec<&str> = output.lines().filter(|l| l.starts_with("Query")).collect();
        assert_eq!(query_lines.len(), 3, "130-char alignment should wrap into 3 blocks (60+60+10)");
        let sbjct_lines: Vec<&str> = output.lines().filter(|l| l.starts_with("Sbjct")).collect();
        assert_eq!(sbjct_lines.len(), 3, "subject should also have 3 lines");
    }
}
