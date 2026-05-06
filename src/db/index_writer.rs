use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

use crate::db::{current_blastdb_date, DbType};

pub(crate) fn write_index_file(
    path: &Path,
    format_version: u32,
    db_type: DbType,
    title: &str,
    num_oids: u32,
    total_length: u64,
    max_seq_len: u32,
    hdr_offsets: &[u32],
    seq_offsets: &[u32],
    amb_offsets: Option<&[u32]>,
) -> io::Result<()> {
    let mut f = BufWriter::new(File::create(path)?);

    f.write_all(&format_version.to_be_bytes())?;
    let db_type_val: u32 = match db_type {
        DbType::Protein => 1,
        DbType::Nucleotide => 0,
    };
    f.write_all(&db_type_val.to_be_bytes())?;

    f.write_all(&(title.len() as u32).to_be_bytes())?;
    f.write_all(title.as_bytes())?;

    let date = current_blastdb_date();
    f.write_all(&(date.len() as u32).to_be_bytes())?;
    f.write_all(date.as_bytes())?;

    f.write_all(&num_oids.to_be_bytes())?;
    // NCBI v4 stores total length little-endian, followed by big-endian max length.
    f.write_all(&total_length.to_le_bytes())?;
    f.write_all(&max_seq_len.to_be_bytes())?;

    for &off in hdr_offsets {
        f.write_all(&off.to_be_bytes())?;
    }
    for &off in seq_offsets {
        f.write_all(&off.to_be_bytes())?;
    }
    if let Some(amb) = amb_offsets {
        for &off in amb {
            f.write_all(&off.to_be_bytes())?;
        }
    }

    f.flush()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn index_writer_uses_v4_endianness_for_lengths() {
        let dir = tempfile::TempDir::new().unwrap();
        let path = dir.path().join("test.nin");
        write_index_file(
            &path,
            4,
            DbType::Nucleotide,
            "db",
            1,
            0x0102_0304_0506_0708,
            99,
            &[0, 10],
            &[0, 5],
            Some(&[0, 5]),
        )
        .unwrap();
        let bytes = std::fs::read(path).unwrap();
        let total_offset = 4 + 4 + 4 + 2 + 4 + current_blastdb_date().len() + 4;
        assert_eq!(
            &bytes[total_offset..total_offset + 8],
            &0x0102_0304_0506_0708u64.to_le_bytes()
        );
        assert_eq!(
            &bytes[total_offset + 8..total_offset + 12],
            &99u32.to_be_bytes()
        );
    }
}
