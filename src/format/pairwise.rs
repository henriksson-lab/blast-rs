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
}
