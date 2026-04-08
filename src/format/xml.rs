//! BLAST XML output format (-outfmt 5).

use std::io::Write;

/// Write BLAST XML header.
pub fn write_xml_header<W: Write>(writer: &mut W, program: &str, version: &str, db: &str) -> std::io::Result<()> {
    writeln!(writer, "<?xml version=\"1.0\"?>")?;
    writeln!(writer, "<!DOCTYPE BlastOutput PUBLIC \"-//NCBI//NCBI BlastOutput/EN\" \"http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd\">")?;
    writeln!(writer, "<BlastOutput>")?;
    writeln!(writer, "  <BlastOutput_program>{}</BlastOutput_program>", program)?;
    writeln!(writer, "  <BlastOutput_version>{} {}</BlastOutput_version>", program, version)?;
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

/// Write one hit in XML format.
pub fn write_xml_hit<W: Write>(
    writer: &mut W,
    hit_num: i32,
    subject_id: &str,
    subject_def: &str,
    subject_len: i32,
    hsps: &[(i32, i32, i32, i32, f64, f64, i32, i32, i32)], // (qs,qe,ss,se,eval,bit,ident,alen,gaps)
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
        writeln!(writer, "          <Hsp_bit-score>{:.4}</Hsp_bit-score>", bit)?;
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
        write_xml_hit(&mut buf, 1, "subj1", "test subject", 1000,
            &[(1, 50, 100, 149, 1e-10, 56.0, 50, 50, 0)]).unwrap();
        write_xml_footer(&mut buf).unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert!(output.contains("<BlastOutput>"));
        assert!(output.contains("subj1"));
        assert!(output.contains("</BlastOutput>"));
    }
}
