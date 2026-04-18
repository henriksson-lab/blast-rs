//! BLAST XML output format (-outfmt 5).

use std::io::Write;

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
        program
    )?;
    writeln!(
        writer,
        "  <BlastOutput_version>{} {}</BlastOutput_version>",
        program, version
    )?;
    writeln!(writer, "  <BlastOutput_db>{}</BlastOutput_db>", db)?;
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
/// bit_score, num_ident, align_length, num_gaps)`.
pub type XmlHsp = (i32, i32, i32, i32, f64, f64, i32, i32, i32);

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
    writeln!(writer, "      <Hit_id>{}</Hit_id>", subject_id)?;
    writeln!(writer, "      <Hit_def>{}</Hit_def>", subject_def)?;
    writeln!(writer, "      <Hit_len>{}</Hit_len>", subject_len)?;
    writeln!(writer, "      <Hit_hsps>")?;

    for (i, &(qs, qe, ss, se, eval, bit, ident, alen, gaps)) in hsps.iter().enumerate() {
        writeln!(writer, "        <Hsp>")?;
        writeln!(writer, "          <Hsp_num>{}</Hsp_num>", i + 1)?;
        writeln!(
            writer,
            "          <Hsp_bit-score>{:.4}</Hsp_bit-score>",
            bit
        )?;
        writeln!(writer, "          <Hsp_evalue>{:.2e}</Hsp_evalue>", eval)?;
        writeln!(writer, "          <Hsp_query-from>{}</Hsp_query-from>", qs)?;
        writeln!(writer, "          <Hsp_query-to>{}</Hsp_query-to>", qe)?;
        writeln!(writer, "          <Hsp_hit-from>{}</Hsp_hit-from>", ss)?;
        writeln!(writer, "          <Hsp_hit-to>{}</Hsp_hit-to>", se)?;
        writeln!(writer, "          <Hsp_identity>{}</Hsp_identity>", ident)?;
        writeln!(writer, "          <Hsp_align-len>{}</Hsp_align-len>", alen)?;
        writeln!(writer, "          <Hsp_gaps>{}</Hsp_gaps>", gaps)?;
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
            &[(1, 50, 100, 149, 1e-10, 56.0, 50, 50, 0)],
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
            &[(10, 60, 200, 250, 1e-20, 100.0, 48, 51, 0)],
        )
        .unwrap();
        write_xml_hit(
            &mut buf,
            2,
            "hit2",
            "second hit",
            800,
            &[(1, 30, 50, 80, 5e-5, 45.0, 28, 31, 1)],
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
            &[(1, 100, 500, 599, 3.5e-40, 156.3, 95, 100, 2)],
        )
        .unwrap();
        write_xml_footer(&mut buf).unwrap();
        let output = String::from_utf8(buf).unwrap();

        // Hit metadata
        assert!(output.contains("<Hit_num>1</Hit_num>"));
        assert!(output.contains("<Hit_id>gi|12345|ref|NM_001.1|</Hit_id>"));
        assert!(output.contains("<Hit_def>Homo sapiens gene</Hit_def>"));
        assert!(output.contains("<Hit_len>2500</Hit_len>"));
        // HSP data fields
        assert!(output.contains("<Hsp_num>1</Hsp_num>"));
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
        assert!(output.contains("<BlastOutput_db>testdb</BlastOutput_db>"));
    }
}
