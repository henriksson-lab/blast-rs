use blast_rs::db::{BlastDb, DbType};
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

    let sequence = db.get_sequence(oid);
    let seq_len = db.get_seq_len(oid) as usize;
    let defline = db.get_defline(oid).unwrap_or_else(|| accession.to_string());
    let mut out = std::io::BufWriter::new(std::io::stdout().lock());
    writeln!(out, ">{}", defline)?;
    match db.db_type {
        DbType::Nucleotide => {
            for chunk_start in (0..seq_len).step_by(60) {
                let chunk_end = (chunk_start + 60).min(seq_len);
                let codes: Vec<u8> = (chunk_start..chunk_end)
                    .map(|pos| blast_rs::encoding::ncbi2na_base_at(sequence, pos))
                    .collect();
                out.write_all(&blast_rs::encoding::blastna_to_iupacna_sequence(&codes))?;
                writeln!(out)?;
            }
        }
        DbType::Protein => {
            let decoded = blast_rs::encoding::ncbistdaa_to_aminoacid_sequence(sequence);
            for chunk in decoded.chunks(60) {
                out.write_all(chunk)?;
                writeln!(out)?;
            }
        }
    }

    Ok(())
}
