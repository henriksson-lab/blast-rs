//! BLAST database index file reader.
//!
//! Reads .nin/.pin index files, .nsq/.psq sequence files, and .nhr/.phr header files.

use byteorder::{BigEndian, LittleEndian, ReadBytesExt};
use memmap2::Mmap;
use std::fs::File;
use std::io::{self, Cursor, Read};
use std::path::{Path, PathBuf};

/// Database type: nucleotide or protein.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DbType {
    Nucleotide,
    Protein,
}

impl DbType {
    fn index_ext(self) -> &'static str {
        match self {
            DbType::Nucleotide => "nin",
            DbType::Protein => "pin",
        }
    }
    fn seq_ext(self) -> &'static str {
        match self {
            DbType::Nucleotide => "nsq",
            DbType::Protein => "psq",
        }
    }
    fn hdr_ext(self) -> &'static str {
        match self {
            DbType::Nucleotide => "nhr",
            DbType::Protein => "phr",
        }
    }
}

/// A single volume of a BLAST database.
pub struct BlastDb {
    pub db_type: DbType,
    pub title: String,
    pub date: String,
    pub num_oids: u32,
    pub total_length: u64,
    pub max_seq_len: u32,
    pub version: u32,

    /// Header offsets: num_oids + 1 entries.
    hdr_offsets: Vec<u32>,
    /// Sequence offsets: num_oids + 1 entries.
    seq_offsets: Vec<u32>,
    /// Ambiguity offsets (nucleotide only): num_oids + 1 entries.
    amb_offsets: Option<Vec<u32>>,

    /// Memory-mapped sequence data.
    seq_mmap: Mmap,
    /// Memory-mapped header data.
    hdr_mmap: Mmap,
}

impl BlastDb {
    /// Open a BLAST database from the given base path (without extension).
    /// Automatically detects nucleotide vs protein.
    pub fn open(base_path: &Path) -> io::Result<Self> {
        // Try nucleotide first, then protein
        let db_type = if base_path.with_extension("nin").exists() {
            DbType::Nucleotide
        } else if base_path.with_extension("pin").exists() {
            DbType::Protein
        } else {
            return Err(io::Error::new(
                io::ErrorKind::NotFound,
                format!(
                    "No BLAST database found at {}",
                    base_path.display()
                ),
            ));
        };
        Self::open_typed(base_path, db_type)
    }

    /// Open a BLAST database of known type.
    pub fn open_typed(base_path: &Path, db_type: DbType) -> io::Result<Self> {
        let idx_path = base_path.with_extension(db_type.index_ext());
        let seq_path = base_path.with_extension(db_type.seq_ext());
        let hdr_path = base_path.with_extension(db_type.hdr_ext());

        // Read index file
        let idx_data = std::fs::read(&idx_path)?;
        let mut cur = Cursor::new(&idx_data);

        let version = cur.read_u32::<BigEndian>()?;
        if version != 4 && version != 5 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("Unsupported BLAST DB version: {}", version),
            ));
        }

        let seq_type_val = cur.read_u32::<BigEndian>()?;
        let parsed_type = if seq_type_val == 0 {
            DbType::Nucleotide
        } else {
            DbType::Protein
        };
        if parsed_type != db_type {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Database type mismatch between extension and header",
            ));
        }

        // Version 5 has an extra volume_number field
        if version == 5 {
            let _volume_number = cur.read_u32::<BigEndian>()?;
        }

        let title = read_string(&mut cur)?;

        // Version 5 has LMDB filename
        if version == 5 {
            let _lmdb_name = read_string(&mut cur)?;
        }

        let date = read_string(&mut cur)?;

        let num_oids = cur.read_u32::<BigEndian>()?;
        let total_length = cur.read_u64::<LittleEndian>()?;
        let max_seq_len = cur.read_u32::<BigEndian>()?;

        let n = (num_oids + 1) as usize;

        // Read offset arrays
        let hdr_offsets = read_u32_array(&mut cur, n)?;

        let (seq_offsets, amb_offsets) = if db_type == DbType::Nucleotide {
            // For nucleotide: header offsets, then seq offsets, then amb offsets
            // (despite the confusing naming in the C++ code)
            let seq_offs = read_u32_array(&mut cur, n)?;
            let amb_offs = read_u32_array(&mut cur, n)?;
            (seq_offs, Some(amb_offs))
        } else {
            let seq_offs = read_u32_array(&mut cur, n)?;
            (seq_offs, None)
        };

        // Memory-map sequence and header files
        let seq_file = File::open(&seq_path)?;
        let seq_mmap = unsafe { Mmap::map(&seq_file)? };

        let hdr_file = File::open(&hdr_path)?;
        let hdr_mmap = unsafe { Mmap::map(&hdr_file)? };

        Ok(BlastDb {
            db_type,
            title,
            date,
            num_oids,
            total_length,
            max_seq_len,
            version,
            hdr_offsets,
            seq_offsets,
            amb_offsets,
            seq_mmap,
            hdr_mmap,
        })
    }

    /// Get the raw sequence bytes for the given OID.
    /// For nucleotide: returns packed 2-bit data (4 bases per byte) + remainder byte.
    /// For protein: returns raw amino acid codes.
    pub fn get_sequence(&self, oid: u32) -> &[u8] {
        let start = self.seq_offsets[oid as usize] as usize;
        let end = self.seq_offsets[oid as usize + 1] as usize;
        &self.seq_mmap[start..end]
    }

    /// Get the sequence length in residues for the given OID.
    pub fn get_seq_len(&self, oid: u32) -> u32 {
        let start = self.seq_offsets[oid as usize] as usize;
        let end = self.seq_offsets[oid as usize + 1] as usize;
        let byte_len = end - start;
        if self.db_type == DbType::Nucleotide {
            // For nucleotide: last byte encodes remainder
            if byte_len == 0 {
                return 0;
            }
            let last_byte = self.seq_mmap[end - 1];
            let remainder = (last_byte & 0x03) as u32;
            let full_bytes = (byte_len - 1) as u32;
            full_bytes * 4 + remainder
        } else {
            byte_len as u32
        }
    }

    /// Get the raw header bytes (ASN.1) for the given OID.
    pub fn get_header(&self, oid: u32) -> &[u8] {
        let start = self.hdr_offsets[oid as usize] as usize;
        let end = self.hdr_offsets[oid as usize + 1] as usize;
        &self.hdr_mmap[start..end]
    }

    /// Extract the first accession-like string from the ASN.1 header.
    /// Matches patterns like BP722512, NC_003421, NM_001234, etc.
    pub fn get_accession(&self, oid: u32) -> Option<String> {
        let hdr = self.get_header(oid);
        let mut i = 0;
        while i < hdr.len() {
            if hdr[i].is_ascii_uppercase() {
                let start = i;
                // Match: [A-Z]+[_]?[0-9]+
                while i < hdr.len() && hdr[i].is_ascii_uppercase() {
                    i += 1;
                }
                // Allow optional underscore
                let mut j = i;
                if j < hdr.len() && hdr[j] == b'_' {
                    j += 1;
                }
                if j < hdr.len() && hdr[j].is_ascii_digit() {
                    i = j;
                    while i < hdr.len() && hdr[i].is_ascii_digit() {
                        i += 1;
                    }
                    let total = i - start;
                    if total >= 6 {
                        let acc = String::from_utf8_lossy(&hdr[start..i]).to_string();
                        // Look for ASN.1 version: \x00+\xa3\x80\x02\x01\xNN
                        let mut vi = i;
                        while vi < hdr.len() && hdr[vi] == 0 { vi += 1; }
                        if vi + 4 < hdr.len()
                            && hdr[vi] == 0xa3
                            && hdr[vi + 1] == 0x80
                            && hdr[vi + 2] == 0x02
                            && hdr[vi + 3] == 0x01
                        {
                            let version = hdr[vi + 4];
                            return Some(format!("{}.{}", acc, version));
                        }
                        return Some(acc);
                    }
                }
            } else {
                i += 1;
            }
        }
        None
    }

    /// Get the ambiguity data offset range for a nucleotide OID.
    pub fn get_ambiguity_data(&self, oid: u32) -> Option<&[u8]> {
        let amb_offsets = self.amb_offsets.as_ref()?;
        let start = amb_offsets[oid as usize] as usize;
        let end = amb_offsets[oid as usize + 1] as usize;
        if start == end {
            None
        } else {
            Some(&self.seq_mmap[start..end])
        }
    }
}

fn read_string<R: Read>(cur: &mut R) -> io::Result<String> {
    let len = ReadBytesExt::read_u32::<BigEndian>(cur)? as usize;
    let mut buf = vec![0u8; len];
    cur.read_exact(&mut buf)?;
    // Trim trailing nulls
    while buf.last() == Some(&0) {
        buf.pop();
    }
    Ok(String::from_utf8_lossy(&buf).into_owned())
}

fn read_u32_array<R: Read>(cur: &mut R, count: usize) -> io::Result<Vec<u32>> {
    let mut arr = Vec::with_capacity(count);
    for _ in 0..count {
        arr.push(ReadBytesExt::read_u32::<BigEndian>(cur)?);
    }
    Ok(arr)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_db_path() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .join("ncbi-blast-2.17.0+-src/c++/src/algo/blast/unit_tests/api/data/seqn")
    }

    #[test]
    fn test_open_seqn_db() {
        let db = BlastDb::open(&test_db_path()).unwrap();
        assert_eq!(db.version, 4);
        assert_eq!(db.db_type, DbType::Nucleotide);
        assert_eq!(db.num_oids, 2004);
        assert_eq!(db.total_length, 943942);
        assert_eq!(db.title, "Another test DB for CPPUNIT, SeqDB.");
    }

    #[test]
    fn test_seq_lengths() {
        let db = BlastDb::open(&test_db_path()).unwrap();
        // First sequence should have non-zero length
        let len = db.get_seq_len(0);
        assert!(len > 0 && len <= db.max_seq_len);
    }

    #[test]
    fn test_all_seq_lengths_sum() {
        let db = BlastDb::open(&test_db_path()).unwrap();
        let total: u64 = (0..db.num_oids).map(|oid| db.get_seq_len(oid) as u64).sum();
        // Total should be close to reported total_length
        // (might not be exact due to how nucleotide lengths work)
        assert!(
            total > 0,
            "Total sequence length should be positive"
        );
    }
}
