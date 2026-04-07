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
        } else if base_path.with_extension("nal").exists() || base_path.with_extension("pal").exists() {
            // Alias file — open first volume (multi-volume merge is TODO)
            let alias_ext = if base_path.with_extension("nal").exists() { "nal" } else { "pal" };
            let alias = crate::alias::parse_alias_file(&base_path.with_extension(alias_ext))?;
            if let Some(first_vol) = alias.dblist.first() {
                return Self::open(first_vol);
            }
            return Err(io::Error::new(io::ErrorKind::NotFound, "Empty alias file"));
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
            // For nucleotide, the index has 3 arrays after the header:
            //   1. Header offsets
            //   2. Sequence start offsets in .nsq (packed NCBI2na data)
            //   3. Ambiguity start offsets in .nsq (= end of sequence data for each OID)
            // Sequence data for OID i: nsq[seq_off[i] .. amb_off[i]]
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
    /// For nucleotide: returns packed 2-bit data from seq_start to amb_start.
    /// For protein: returns raw amino acid codes.
    pub fn get_sequence(&self, oid: u32) -> &[u8] {
        let start = self.seq_offsets[oid as usize] as usize;
        let end = if self.db_type == DbType::Nucleotide {
            // For nucleotide, sequence data ends at the ambiguity offset
            self.amb_offsets.as_ref().unwrap()[oid as usize] as usize
        } else {
            self.seq_offsets[oid as usize + 1] as usize
        };
        &self.seq_mmap[start..end]
    }

    /// Get the sequence length in residues for the given OID.
    pub fn get_seq_len(&self, oid: u32) -> u32 {
        let start = self.seq_offsets[oid as usize] as usize;
        if self.db_type == DbType::Nucleotide {
            // For nucleotide: sequence ends at amb_offset
            // Length = (whole_bytes * 4) + remainder
            // where whole_bytes = amb_off - seq_off - 1, remainder = last_byte & 3
            let end = self.amb_offsets.as_ref().unwrap()[oid as usize] as usize;
            let byte_len = end - start;
            if byte_len == 0 {
                return 0;
            }
            let whole_bytes = (byte_len - 1) as u32;
            let last_byte = self.seq_mmap[start + whole_bytes as usize];
            let remainder = (last_byte & 0x03) as u32;
            whole_bytes * 4 + remainder
        } else {
            let end = self.seq_offsets[oid as usize + 1] as usize;
            (end - start) as u32
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
                        // Check for text-format version: ".N" right after accession
                        if i + 1 < hdr.len() && hdr[i] == b'.' && hdr[i + 1].is_ascii_digit() {
                            let _dot_start = i;
                            i += 1;
                            while i < hdr.len() && hdr[i].is_ascii_digit() { i += 1; }
                            return Some(String::from_utf8_lossy(&hdr[start..i]).to_string());
                        }
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

    /// Get the ambiguity data for a nucleotide OID.
    /// Ambiguity data is stored between amb_offset and the next sequence's start.
    pub fn get_ambiguity_data(&self, oid: u32) -> Option<&[u8]> {
        let amb_offsets = self.amb_offsets.as_ref()?;
        let start = amb_offsets[oid as usize] as usize;
        // Ambiguity data ends at the next sequence's start offset
        let end = self.seq_offsets[oid as usize + 1] as usize;
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
        assert_eq!(total, db.total_length,
            "Sum of individual lengths must match reported total_length");
    }

    #[test]
    fn test_seq_len_oid0_matches_reference() {
        let db = BlastDb::open(&test_db_path()).unwrap();
        // Reference BLAST reports slen=386 for OID 0
        assert_eq!(db.get_seq_len(0), 386);
    }

    #[test]
    fn test_accession_seqn() {
        let db = BlastDb::open(&test_db_path()).unwrap();
        let acc = db.get_accession(0);
        assert_eq!(acc.as_deref(), Some("BP722512.1"));
    }

    #[test]
    fn test_accession_pombe() {
        let pombe = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent().unwrap()
            .join("ncbi-blast-2.17.0+-src/c++/src/algo/blast/unit_tests/api/data/pombe");
        if !pombe.with_extension("nin").exists() { return; }
        let db = BlastDb::open(&pombe).unwrap();
        assert_eq!(db.get_accession(0).as_deref(), Some("NC_003421.2"));
    }

    // ---- Validation / error handling tests ----

    fn corrupt_db_path(name: &str) -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent().unwrap()
            .join(format!("tests/fixtures/corrupt_db/{}", name))
    }

    #[test]
    fn test_missing_database() {
        let result = BlastDb::open(&PathBuf::from("/nonexistent/db"));
        assert!(result.is_err());
    }

    #[test]
    fn test_truncated_index() {
        let result = BlastDb::open(&corrupt_db_path("truncated"));
        assert!(result.is_err(), "Truncated index should fail to open");
    }

    #[test]
    fn test_bad_version() {
        let result = BlastDb::open(&corrupt_db_path("badversion"));
        assert!(result.is_err(), "Unsupported version should fail");
    }

    #[test]
    fn test_empty_database() {
        let result = BlastDb::open(&corrupt_db_path("empty"));
        // Empty DB should either open with 0 sequences or fail gracefully
        match result {
            Ok(db) => assert_eq!(db.num_oids, 0),
            Err(_) => {} // Also acceptable
        }
    }

    #[test]
    fn test_wrong_type_detection() {
        // Try opening a nucleotide DB as protein
        let path = test_db_path();
        let result = BlastDb::open_typed(&path, DbType::Protein);
        assert!(result.is_err(), "Opening nuc DB as protein should fail");
    }

    #[test]
    fn test_ambiguity_data() {
        let db = BlastDb::open(&test_db_path()).unwrap();
        // OID 0 should have ambiguity data (1 N at position 293)
        let amb = db.get_ambiguity_data(0);
        assert!(amb.is_some(), "OID 0 should have ambiguity data");
        let amb = amb.unwrap();
        assert!(amb.len() >= 8, "Ambiguity data should have header + at least 1 record");
    }

    #[test]
    fn test_header_data() {
        let db = BlastDb::open(&test_db_path()).unwrap();
        let hdr = db.get_header(0);
        assert!(!hdr.is_empty(), "Header should not be empty");
        // Should contain the title string somewhere
        assert!(hdr.windows(8).any(|w| w == b"Xenopus " || w == b"BP722512"),
            "Header should contain sequence info");
    }
}

    #[test]
    fn test_open_protein_db() {
        let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent().unwrap()
            .join("ncbi-blast-2.17.0+-src/c++/src/algo/blast/unit_tests/api/data/seqp");
        if !path.with_extension("pin").exists() { return; }
        let db = BlastDb::open(&path).unwrap();
        assert_eq!(db.db_type, DbType::Protein);
        assert_eq!(db.version, 4);
        assert_eq!(db.num_oids, 2005);
        assert!(db.total_length > 0);
    }

    #[test]
    fn test_protein_seq_length() {
        let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent().unwrap()
            .join("ncbi-blast-2.17.0+-src/c++/src/algo/blast/unit_tests/api/data/seqp");
        if !path.with_extension("pin").exists() { return; }
        let db = BlastDb::open(&path).unwrap();
        let len = db.get_seq_len(0);
        assert!(len > 0 && len <= db.max_seq_len, "Protein seq length should be valid");
    }
