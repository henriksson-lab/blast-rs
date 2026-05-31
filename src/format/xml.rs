//! BLAST XML output format (-outfmt 5).

use std::io::Write;

fn xml_escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for ch in s.chars() {
        match ch {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            '\'' => out.push_str("&apos;"),
            _ => out.push(ch),
        }
    }
    out
}

fn format_xml_double_g(value: f64) -> String {
    if value == 0.0 {
        return "0".to_string();
    }
    let abs = value.abs();
    if abs < 1e-4 || abs >= 1e6 {
        let s = format!("{value:.5e}");
        if let Some((mantissa, exponent)) = s.split_once('e') {
            let mut mantissa = mantissa.trim_end_matches('0').to_string();
            if mantissa.ends_with('.') {
                mantissa.pop();
            }
            let exponent = exponent
                .trim_start_matches('+')
                .trim_start_matches('0')
                .replace("-0", "-");
            format!("{mantissa}e{exponent}")
        } else {
            s
        }
    } else {
        let exp = abs.log10().floor() as i32;
        let decimals = (5 - exp).max(0) as usize;
        let mut s = format!("{value:.decimals$}");
        if s.contains('.') {
            while s.ends_with('0') {
                s.pop();
            }
            if s.ends_with('.') {
                s.pop();
            }
        }
        s
    }
}

fn format_xml_evalue(value: f64) -> String {
    if value == 0.0 || value < 1.0e-180 {
        "0".to_string()
    } else {
        format_xml_double_g(value)
    }
}

fn xml_hit_accession(subject_id: &str) -> &str {
    let parts: Vec<&str> = subject_id
        .split('|')
        .filter(|part| !part.is_empty())
        .collect();
    let accession = if parts.len() >= 3 {
        parts[parts.len() - 1]
    } else {
        subject_id
    };
    if let Some(dot) = accession.rfind('.') {
        let suffix = &accession[dot + 1..];
        if !suffix.is_empty() && suffix.bytes().all(|b| b.is_ascii_digit()) {
            return &accession[..dot];
        }
    }
    accession
}

const BLAST_XML_REFERENCE: &str = "Stephen F. Altschul, Thomas L. Madden, Alejandro A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), &quot;Gapped BLAST and PSI-BLAST: a new generation of protein database search programs&quot;, Nucleic Acids Res. 25:3389-3402.";

/// Write BLAST XML header.
pub fn write_xml_header<W: Write>(
    writer: &mut W,
    program: &str,
    version: &str,
    db: &str,
) -> std::io::Result<()> {
    writeln!(writer, "<?xml version=\"1.0\"?>")?;
    writeln!(writer, "<!DOCTYPE BlastOutput PUBLIC \"-//NCBI//NCBI BlastOutput/EN\" \"http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd\">")?;
    writeln!(writer, "<BlastOutput>")?;
    writeln!(
        writer,
        "  <BlastOutput_program>{}</BlastOutput_program>",
        xml_escape(program)
    )?;
    writeln!(
        writer,
        "  <BlastOutput_version>{} {}</BlastOutput_version>",
        xml_escape(program),
        xml_escape(version)
    )?;
    writeln!(
        writer,
        "  <BlastOutput_reference>{}</BlastOutput_reference>",
        BLAST_XML_REFERENCE
    )?;
    writeln!(
        writer,
        "  <BlastOutput_db>{}</BlastOutput_db>",
        xml_escape(db)
    )?;
    writeln!(writer, "  <BlastOutput_iterations>")?;
    Ok(())
}

/// Write XML footer.
pub fn write_xml_footer<W: Write>(writer: &mut W) -> std::io::Result<()> {
    writeln!(writer, "  </BlastOutput_iterations>")?;
    writeln!(writer, "</BlastOutput>")?;
    Ok(())
}

/// One HSP's XML fields, packed in NCBI's `<Hsp>` element order:
/// `(query_start, query_end, subject_start, subject_end, evalue,
/// bit_score, raw_score, num_ident, align_length, num_gaps)`.
pub type XmlHsp = (i32, i32, i32, i32, f64, f64, i32, i32, i32, i32);

/// Write one hit in XML format.
pub fn write_xml_hit<W: Write>(
    writer: &mut W,
    hit_num: i32,
    subject_id: &str,
    subject_def: &str,
    subject_len: i32,
    hsps: &[XmlHsp],
) -> std::io::Result<()> {
    writeln!(writer, "    <Hit>")?;
    writeln!(writer, "      <Hit_num>{}</Hit_num>", hit_num)?;
    writeln!(writer, "      <Hit_id>{}</Hit_id>", xml_escape(subject_id))?;
    writeln!(
        writer,
        "      <Hit_def>{}</Hit_def>",
        xml_escape(subject_def)
    )?;
    writeln!(
        writer,
        "      <Hit_accession>{}</Hit_accession>",
        xml_escape(xml_hit_accession(subject_id))
    )?;
    writeln!(writer, "      <Hit_len>{}</Hit_len>", subject_len)?;
    writeln!(writer, "      <Hit_hsps>")?;

    for (i, &(qs, qe, ss, se, eval, bit, raw_score, ident, alen, gaps)) in hsps.iter().enumerate() {
        writeln!(writer, "        <Hsp>")?;
        writeln!(writer, "          <Hsp_num>{}</Hsp_num>", i + 1)?;
        writeln!(
            writer,
            "          <Hsp_bit-score>{}</Hsp_bit-score>",
            format_xml_double_g(bit)
        )?;
        writeln!(writer, "          <Hsp_score>{}</Hsp_score>", raw_score)?;
        writeln!(
            writer,
            "          <Hsp_evalue>{}</Hsp_evalue>",
            format_xml_evalue(eval)
        )?;
        writeln!(writer, "          <Hsp_query-from>{}</Hsp_query-from>", qs)?;
        writeln!(writer, "          <Hsp_query-to>{}</Hsp_query-to>", qe)?;
        writeln!(writer, "          <Hsp_hit-from>{}</Hsp_hit-from>", ss)?;
        writeln!(writer, "          <Hsp_hit-to>{}</Hsp_hit-to>", se)?;
        writeln!(writer, "          <Hsp_identity>{}</Hsp_identity>", ident)?;
        writeln!(writer, "          <Hsp_gaps>{}</Hsp_gaps>", gaps)?;
        writeln!(writer, "          <Hsp_align-len>{}</Hsp_align-len>", alen)?;
        writeln!(writer, "        </Hsp>")?;
    }

    writeln!(writer, "      </Hit_hsps>")?;
    writeln!(writer, "    </Hit>")?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_xml_output() {
        let mut buf = Vec::new();
        write_xml_header(&mut buf, "blastn", "0.1.0", "testdb").unwrap();
        write_xml_hit(
            &mut buf,
            1,
            "subj1",
            "test subject",
            1000,
            &[(1, 50, 100, 149, 1e-10, 56.0, 104, 50, 50, 0)],
        )
        .unwrap();
        write_xml_footer(&mut buf).unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert!(output.contains("<BlastOutput>"));
        assert!(output.contains("subj1"));
        assert!(output.contains("</BlastOutput>"));
    }

    #[test]
    fn test_xml_well_formed() {
        let mut buf = Vec::new();
        write_xml_header(&mut buf, "blastn", "0.1.0", "mydb").unwrap();
        write_xml_hit(
            &mut buf,
            1,
            "hit1",
            "first hit",
            500,
            &[(10, 60, 200, 250, 1e-20, 100.0, 184, 48, 51, 0)],
        )
        .unwrap();
        write_xml_hit(
            &mut buf,
            2,
            "hit2",
            "second hit",
            800,
            &[(1, 30, 50, 80, 5e-5, 45.0, 82, 28, 31, 1)],
        )
        .unwrap();
        write_xml_footer(&mut buf).unwrap();
        let output = String::from_utf8(buf).unwrap();

        // Verify XML declaration
        assert!(
            output.starts_with("<?xml version=\"1.0\"?>"),
            "should start with XML declaration"
        );
        // Verify matching open/close tags for key elements
        assert!(output.contains("<BlastOutput>"));
        assert!(output.contains("</BlastOutput>"));
        assert!(output.contains("<BlastOutput_iterations>"));
        assert!(output.contains("</BlastOutput_iterations>"));
        assert!(output.contains("<Hit>"));
        assert!(output.contains("</Hit>"));
        assert!(output.contains("<Hit_hsps>"));
        assert!(output.contains("</Hit_hsps>"));
        assert!(output.contains("<Hsp>"));
        assert!(output.contains("</Hsp>"));
        // Count that opening and closing tags match
        let hit_opens = output.matches("<Hit>").count();
        let hit_closes = output.matches("</Hit>").count();
        assert_eq!(hit_opens, hit_closes, "Hit open/close tags must match");
        assert_eq!(hit_opens, 2, "should have 2 hits");
    }

    #[test]
    fn test_xml_contains_hit_data() {
        let mut buf = Vec::new();
        write_xml_header(&mut buf, "blastn", "0.1.0", "testdb").unwrap();
        write_xml_hit(
            &mut buf,
            1,
            "gi|12345|ref|NM_001.1|",
            "Homo sapiens gene",
            2500,
            &[(1, 100, 500, 599, 3.5e-40, 156.3, 289, 95, 100, 2)],
        )
        .unwrap();
        write_xml_footer(&mut buf).unwrap();
        let output = String::from_utf8(buf).unwrap();

        // Hit metadata
        assert!(output.contains("<Hit_num>1</Hit_num>"));
        assert!(output.contains("<Hit_id>gi|12345|ref|NM_001.1|</Hit_id>"));
        assert!(output.contains("<Hit_def>Homo sapiens gene</Hit_def>"));
        assert!(output.contains("<Hit_accession>NM_001</Hit_accession>"));
        assert!(output.contains("<Hit_len>2500</Hit_len>"));
        // HSP data fields
        assert!(output.contains("<Hsp_num>1</Hsp_num>"));
        assert!(output.contains("<Hsp_score>289</Hsp_score>"));
        assert!(output.contains("<Hsp_query-from>1</Hsp_query-from>"));
        assert!(output.contains("<Hsp_query-to>100</Hsp_query-to>"));
        assert!(output.contains("<Hsp_hit-from>500</Hsp_hit-from>"));
        assert!(output.contains("<Hsp_hit-to>599</Hsp_hit-to>"));
        assert!(output.contains("<Hsp_identity>95</Hsp_identity>"));
        assert!(output.contains("<Hsp_align-len>100</Hsp_align-len>"));
        assert!(output.contains("<Hsp_gaps>2</Hsp_gaps>"));
        // E-value and bit-score should be present
        assert!(output.contains("<Hsp_evalue>"));
        assert!(output.contains("<Hsp_bit-score>"));
        // Program info
        assert!(output.contains("<BlastOutput_program>blastn</BlastOutput_program>"));
        assert!(output.contains("<BlastOutput_reference>Stephen F. Altschul"));
        assert!(output.contains("<BlastOutput_db>testdb</BlastOutput_db>"));
    }

    #[test]
    fn test_xml_escapes_text_fields() {
        let mut buf = Vec::new();
        write_xml_header(&mut buf, "blast<xn", "v\"1", "db&a").unwrap();
        write_xml_hit(
            &mut buf,
            1,
            "hit<&\"'",
            "def>with&chars",
            10,
            &[(1, 2, 3, 4, 1e-10, 56.0, 104, 2, 2, 0)],
        )
        .unwrap();
        write_xml_footer(&mut buf).unwrap();
        let output = String::from_utf8(buf).unwrap();

        assert!(output.contains("<BlastOutput_program>blast&lt;xn</BlastOutput_program>"));
        assert!(output.contains("<BlastOutput_version>blast&lt;xn v&quot;1</BlastOutput_version>"));
        assert!(output.contains("<BlastOutput_db>db&amp;a</BlastOutput_db>"));
        assert!(output.contains("<Hit_id>hit&lt;&amp;&quot;&apos;</Hit_id>"));
        assert!(output.contains("<Hit_def>def&gt;with&amp;chars</Hit_def>"));
        assert!(output.contains("<Hit_accession>hit&lt;&amp;&quot;&apos;</Hit_accession>"));
    }

    #[test]
    fn test_xml_header_reference_matches_ncbi_order() {
        let mut buf = Vec::new();
        write_xml_header(&mut buf, "blastn", "2.17.0+", "testdb").unwrap();
        let output = String::from_utf8(buf).unwrap();

        assert!(output.contains("<BlastOutput_reference>Stephen F. Altschul"));
        assert!(
            output.find("<BlastOutput_version>").unwrap()
                < output.find("<BlastOutput_reference>").unwrap()
                && output.find("<BlastOutput_reference>").unwrap()
                    < output.find("<BlastOutput_db>").unwrap(),
            "BlastOutput_reference should follow version and precede db"
        );
    }

    #[test]
    fn test_xml_hit_accession_matches_ncbi_seqid_shape() {
        let mut buf = Vec::new();
        write_xml_hit(
            &mut buf,
            1,
            "gi|12345|ref|NM_001.1|",
            "def",
            10,
            &[(1, 2, 3, 4, 1e-10, 56.0, 104, 2, 2, 0)],
        )
        .unwrap();
        let output = String::from_utf8(buf).unwrap();

        assert!(output.contains("<Hit_id>gi|12345|ref|NM_001.1|</Hit_id>"));
        assert!(output.contains("<Hit_accession>NM_001</Hit_accession>"));
        assert!(
            output.find("<Hit_def>").unwrap() < output.find("<Hit_accession>").unwrap()
                && output.find("<Hit_accession>").unwrap() < output.find("<Hit_len>").unwrap(),
            "Hit_accession should follow BLAST XML Hit_def and precede Hit_len"
        );
    }

    #[test]
    fn test_xml_hsp_numbers_use_ncbi_general_format() {
        let mut buf = Vec::new();
        write_xml_hit(
            &mut buf,
            1,
            "hit",
            "def",
            10,
            &[
                (1, 2, 3, 4, 1e-10, 56.0, 104, 2, 2, 0),
                (1, 2, 3, 4, 0.000392118, 12.3456, 23, 2, 2, 0),
                (1, 2, 3, 4, 3.5e-40, 156.3, 289, 2, 2, 0),
                (1, 2, 3, 4, 1e-200, 1_234_567.0, 2_000_000, 2, 2, 0),
            ],
        )
        .unwrap();
        let output = String::from_utf8(buf).unwrap();

        assert!(output.contains("<Hsp_bit-score>56</Hsp_bit-score>"));
        assert!(output.contains("<Hsp_bit-score>156.3</Hsp_bit-score>"));
        assert!(output.contains("<Hsp_bit-score>1.23457e6</Hsp_bit-score>"));
        assert!(output.contains("<Hsp_evalue>1e-10</Hsp_evalue>"));
        assert!(output.contains("<Hsp_evalue>0.000392118</Hsp_evalue>"));
        assert!(output.contains("<Hsp_evalue>3.5e-40</Hsp_evalue>"));
        assert!(output.contains("<Hsp_evalue>0</Hsp_evalue>"));
    }

    #[test]
    fn test_xml_hsp_score_matches_ncbi_field_order() {
        let mut buf = Vec::new();
        write_xml_hit(
            &mut buf,
            1,
            "hit",
            "def",
            10,
            &[(1, 2, 3, 4, 1e-10, 56.0, 104, 2, 2, 0)],
        )
        .unwrap();
        let output = String::from_utf8(buf).unwrap();

        assert!(output.contains("<Hsp_score>104</Hsp_score>"));
        assert!(
            output.find("<Hsp_bit-score>").unwrap() < output.find("<Hsp_score>").unwrap()
                && output.find("<Hsp_score>").unwrap() < output.find("<Hsp_evalue>").unwrap(),
            "Hsp_score should follow Hsp_bit-score and precede Hsp_evalue"
        );
    }

    #[test]
    fn test_xml_hsp_gaps_precede_align_len() {
        let mut buf = Vec::new();
        write_xml_hit(
            &mut buf,
            1,
            "hit",
            "def",
            10,
            &[(1, 2, 3, 4, 1e-10, 56.0, 104, 2, 2, 1)],
        )
        .unwrap();
        let output = String::from_utf8(buf).unwrap();

        assert!(output.contains("<Hsp_gaps>1</Hsp_gaps>"));
        assert!(output.contains("<Hsp_align-len>2</Hsp_align-len>"));
        assert!(
            output.find("<Hsp_identity>").unwrap() < output.find("<Hsp_gaps>").unwrap()
                && output.find("<Hsp_gaps>").unwrap() < output.find("<Hsp_align-len>").unwrap(),
            "Hsp_gaps should follow Hsp_identity and precede Hsp_align-len"
        );
    }
}
