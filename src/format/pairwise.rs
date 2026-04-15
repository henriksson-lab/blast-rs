//! Pairwise text output format (-outfmt 0).
//! The traditional BLAST output with alignment visualization.

use std::io::Write;

pub fn format_pairwise_evalue(val: f64) -> String {
    if val < 1.0e-180 {
        "0.0".to_string()
    } else if val < 0.001 {
        let s = format!("{:.0e}", val);
        if let Some(e_pos) = s.find('e') {
            let (mantissa, exp_part) = s.split_at(e_pos);
            let exp_str = &exp_part[1..];
            let (sign, digits) = if let Some(rest) = exp_str.strip_prefix('-') {
                ("-", rest)
            } else if let Some(rest) = exp_str.strip_prefix('+') {
                ("", rest)
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
                format!("{}e{}{}", mantissa, sign, digits)
            }
        } else {
            s
        }
    } else if val < 0.1 {
        format!("{:.3}", val)
    } else if val < 1.0 {
        format!("{:.2}", val)
    } else if val < 10.0 {
        format!("{:.1}", val)
    } else {
        format!("{:.0}", val)
    }
}

/// Format a single alignment in pairwise text format.
pub fn format_pairwise_alignment<W: Write>(
    writer: &mut W,
    query_id: &str,
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
    format_pairwise_alignment_with_header(
        writer,
        query_id,
        subject_id,
        query_seq,
        subject_seq,
        q_start,
        q_end,
        s_start,
        s_end,
        score,
        bit_score,
        evalue,
        num_ident,
        align_len,
        gap_opens,
        true,
    )
}

/// Format one HSP in pairwise text format, optionally suppressing the subject
/// header when another HSP for the same subject has already introduced it.
pub fn format_pairwise_alignment_with_header<W: Write>(
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
    show_subject_header: bool,
) -> std::io::Result<()> {
    let pident = if align_len > 0 {
        100.0 * num_ident as f64 / align_len as f64
    } else {
        0.0
    };

    if show_subject_header {
        write_wrapped_subject_header(writer, subject_id)?;
    }
    writeln!(writer, "")?;
    writeln!(
        writer,
        " Score = {:.1} bits ({}),  Expect = {}",
        bit_score,
        score,
        format_pairwise_evalue(evalue)
    )?;
    writeln!(
        writer,
        " Identities = {}/{} ({}%), Gaps = {}/{} ({}%)",
        num_ident,
        align_len,
        (pident + 0.5) as i32,
        gap_opens,
        align_len,
        if align_len > 0 {
            (100.0 * gap_opens as f64 / align_len as f64 + 0.5) as i32
        } else {
            0
        }
    )?;
    writeln!(
        writer,
        " Strand=Plus/{}",
        if s_start <= s_end { "Plus" } else { "Minus" }
    )?;
    writeln!(writer, "")?;

    // Display alignment in 60-character lines
    let blastna_to_char = |b: u8| -> char {
        match b {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            3 => 'T',
            4 => 'R',
            5 => 'Y',
            6 => 'M',
            7 => 'K',
            8 => 'W',
            9 => 'S',
            10 => 'B',
            11 => 'D',
            12 => 'H',
            13 => 'V',
            14 => 'N',
            _ => '-',
        }
    };

    let line_width = 60;
    let coord_width = q_start
        .abs()
        .max(q_end.abs())
        .max(s_start.abs())
        .max(s_end.abs())
        .to_string()
        .len();
    let sequence_column = 5 + 2 + coord_width + 2;
    let q_len = if query_seq.is_empty() || subject_seq.is_empty() {
        (q_end - q_start).unsigned_abs() as usize + 1
    } else {
        query_seq.len().max(subject_seq.len())
    };
    let mut pos = 0;
    let mut qi = q_start;
    let mut si = s_start;
    let s_dir: i32 = if s_start <= s_end { 1 } else { -1 };

    while pos < q_len {
        let chunk = line_width.min(q_len - pos);

        // Query line
        write!(writer, "Query  {:<width$}  ", qi, width = coord_width)?;
        for i in 0..chunk {
            if pos + i < query_seq.len() {
                write!(writer, "{}", blastna_to_char(query_seq[pos + i]))?;
            }
        }
        writeln!(writer, "  {}", qi + chunk as i32 - 1)?;

        // Match line
        for _ in 0..sequence_column {
            write!(writer, " ")?;
        }
        for i in 0..chunk {
            if pos + i < query_seq.len()
                && pos + i < subject_seq.len()
                && query_seq[pos + i] == subject_seq[pos + i]
            {
                write!(writer, "|")?;
            } else {
                write!(writer, " ")?;
            }
        }
        writeln!(writer)?;

        // Subject line
        write!(writer, "Sbjct  {:<width$}  ", si, width = coord_width)?;
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

fn write_wrapped_subject_header<W: Write>(writer: &mut W, subject_id: &str) -> std::io::Result<()> {
    const MAX_HEADER_WIDTH: usize = 80;
    let mut line = format!(">{}", subject_id);
    while line.len() > MAX_HEADER_WIDTH {
        let split = line[..MAX_HEADER_WIDTH]
            .rfind(' ')
            .filter(|&idx| idx > 0)
            .unwrap_or(MAX_HEADER_WIDTH);
        writeln!(writer, "{}", &line[..split + 1])?;
        line = line[split + 1..].to_string();
    }
    writeln!(writer, "{}", line)
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
            &mut buf, "query1", "subject1", &query, &subject, 1, 8, 1, 8, 16, 30.0, 1e-5, 8, 8, 0,
        )
        .unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert!(output.contains(">subject1"));
        assert!(output.contains("Score"));
        assert!(output.contains("Query"));
        assert!(output.contains("Sbjct"));
    }

    #[test]
    fn test_pairwise_coordinate_width_matches_blast_columns() {
        let mut buf = Vec::new();
        let query = vec![0u8, 1, 2, 3];
        let subject = query.clone();
        format_pairwise_alignment(
            &mut buf, "query1", "subject1", &query, &subject, 1, 4, 10861758, 10861761, 8, 16.0,
            1.0, 4, 4, 0,
        )
        .unwrap();
        let output = String::from_utf8(buf).unwrap();
        let query_line = output
            .lines()
            .find(|line| line.starts_with("Query"))
            .unwrap();
        let match_line = output.lines().find(|line| line.contains("||||")).unwrap();
        let sbjct_line = output
            .lines()
            .find(|line| line.starts_with("Sbjct"))
            .unwrap();
        assert_eq!(query_line.find("ACGT"), Some(17));
        assert_eq!(match_line.find("||||"), Some(17));
        assert_eq!(sbjct_line.find("ACGT"), Some(17));
    }

    #[test]
    fn test_pairwise_wraps_long_subject_header() {
        let mut buf = Vec::new();
        write_wrapped_subject_header(
            &mut buf,
            "JQ990143.1 Caretta caretta glyceraldehyde-3-phosphate dehydrogenase mRNA, partial cds",
        )
        .unwrap();
        let output = String::from_utf8(buf).unwrap();
        let lines: Vec<_> = output.lines().collect();
        assert_eq!(
            lines,
            vec![
                ">JQ990143.1 Caretta caretta glyceraldehyde-3-phosphate dehydrogenase mRNA, ",
                "partial cds",
            ]
        );
    }

    #[test]
    fn test_pairwise_with_gaps() {
        // Query:   ACGT-ACGT (gap at position 4 in query => '-' = code 15 or 0xFF)
        // Subject: ACGTGACGT
        // BLASTNA: A=0, C=1, G=2, T=3, gap=15
        let query = vec![0u8, 1, 2, 3, 15, 0, 1, 2, 3];
        let subject = vec![0u8, 1, 2, 3, 2, 0, 1, 2, 3];
        let mut buf = Vec::new();
        format_pairwise_alignment(
            &mut buf, "query1", "subject1", &query, &subject, 1, 9, 1, 9, 14, 28.0, 1e-3, 8, 9, 1,
        )
        .unwrap();
        let output = String::from_utf8(buf).unwrap();
        // The gap character should appear as '-' in query line
        assert!(
            output.contains("-"),
            "query line should contain a gap character '-'"
        );
        // Should show identities and gaps info
        assert!(output.contains("Gaps = 1/9"));
        // Match line should have pipe where bases match
        assert!(
            output.contains("|"),
            "match line should contain pipe characters"
        );
    }

    #[test]
    fn test_pairwise_minus_strand() {
        let query = vec![0u8, 1, 2, 3, 0, 1, 2, 3]; // ACGTACGT
        let subject = vec![0u8, 1, 2, 3, 0, 1, 2, 3];
        let mut buf = Vec::new();
        // Minus strand: s_start > s_end
        format_pairwise_alignment(
            &mut buf, "query1", "subject1", &query, &subject, 1, 8, 100,
            93, // subject coords decrease for minus strand
            16, 30.0, 1e-5, 8, 8, 0,
        )
        .unwrap();
        let output = String::from_utf8(buf).unwrap();
        // Should indicate Minus strand
        assert!(
            output.contains("Strand=Plus/Minus"),
            "output should indicate minus strand"
        );
        // Subject should start at 100 and decrease
        assert!(
            output.contains("Sbjct  100"),
            "subject line should start at position 100"
        );
    }

    #[test]
    fn test_pairwise_long_alignment() {
        // Create an alignment longer than 60 characters to test wrapping
        let len = 130;
        let query: Vec<u8> = (0..len).map(|i| (i % 4) as u8).collect(); // repeating ACGT
        let subject: Vec<u8> = (0..len).map(|i| (i % 4) as u8).collect();
        let mut buf = Vec::new();
        format_pairwise_alignment(
            &mut buf, "query1", "subject1", &query, &subject, 1, len as i32, 1, len as i32, 260,
            480.0, 1e-100, len as i32, len as i32, 0,
        )
        .unwrap();
        let output = String::from_utf8(buf).unwrap();
        // Count "Query" lines — should be ceil(130/60) = 3 blocks
        let query_lines: Vec<&str> = output.lines().filter(|l| l.starts_with("Query")).collect();
        assert_eq!(
            query_lines.len(),
            3,
            "130-char alignment should wrap into 3 blocks (60+60+10)"
        );
        let sbjct_lines: Vec<&str> = output.lines().filter(|l| l.starts_with("Sbjct")).collect();
        assert_eq!(sbjct_lines.len(), 3, "subject should also have 3 lines");
    }
}
