use blast_rs::db::BlastDb;
use std::io::Write;
use std::path::PathBuf;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = std::env::args_os();
    let _prog = args.next();
    let db_path = PathBuf::from(args.next().ok_or("usage: dump_db_seq <db> <accession>")?);
    let accession = args.next().ok_or("usage: dump_db_seq <db> <accession>")?;
    let accession = accession.to_string_lossy();

    let db = BlastDb::open(&db_path)?;
    let oid = (0..db.num_oids)
        .find(|&oid| db.get_accession(oid).as_deref() == Some(accession.as_ref()))
        .ok_or("accession not found")?;

    let packed = db.get_sequence(oid);
    let seq_len = db.get_seq_len(oid) as usize;
    let defline = db
        .get_defline(oid)
        .unwrap_or_else(|| accession.to_string());
    let mut out = std::io::BufWriter::new(std::io::stdout().lock());
    writeln!(out, ">{}", defline)?;
    for chunk_start in (0..seq_len).step_by(60) {
        let chunk_end = (chunk_start + 60).min(seq_len);
        for pos in chunk_start..chunk_end {
            let base = packed[pos >> 2];
            let code = (base >> (6 - 2 * (pos & 3))) & 3;
            let ch = match code {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                3 => b'T',
                _ => b'N',
            };
            out.write_all(&[ch])?;
        }
        writeln!(out)?;
    }

    Ok(())
}
