//! BLAST database index file reader.
//!
//! Reads .nin/.pin index files, .nsq/.psq sequence files, and .nhr/.phr header files.

use byteorder::{BigEndian, LittleEndian, ReadBytesExt};
use memmap2::Mmap;
use std::fs::File;
use std::io::{self, Cursor, Read};
use std::path::Path;
#[cfg(test)]
use std::path::PathBuf;

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
    /// OID-to-taxid lookup (from .nto file, v5 databases).
    tax_lookup: Option<TaxIdLookup>,
}

/// OID-to-taxid lookup table parsed from a .nto file.
struct TaxIdLookup {
    /// Index: (num_oids+1) u64 entries — cumulative end offsets into `taxids`.
    index: Vec<u64>,
    /// Flat array of i32 taxids, indexed by the `index` array.
    taxids: Vec<i32>,
}

impl TaxIdLookup {
    /// Parse a .nto file. Supports two formats:
    ///
    /// **V4 format** (older databases):
    ///   u64         num_oids
    ///   u64[n+1]    cumulative end offsets into taxid data
    ///   i32[]       flat taxid array
    ///
    /// **V5 format** (BLAST+ v5 databases):
    ///   Flat per-OID records: for each OID: u32_LE count, then count × i32_LE taxids.
    ///   No header — the file is just a concatenation of per-OID records.
    fn from_file(path: &Path) -> io::Result<Self> {
        Self::from_file_with_hint(path, None)
    }

    fn from_file_with_hint(path: &Path, expected_oids: Option<u32>) -> io::Result<Self> {
        let data = std::fs::read(path)?;
        if data.len() < 8 {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "nto file too small"));
        }

        // Try V4 format first: starts with u64 num_oids.
        // Validate: num_oids must be reasonable (file must be large enough).
        let mut cur = Cursor::new(&data);
        let candidate_num_oids = cur.read_u64::<LittleEndian>()? as usize;
        let v4_index_bytes = (candidate_num_oids + 1) * 8;
        let v4_valid = candidate_num_oids > 0
            && candidate_num_oids < 1_000_000_000
            && data.len() >= 8 + v4_index_bytes
            && (expected_oids.is_none() || candidate_num_oids == expected_oids.unwrap() as usize);

        if v4_valid {
            // V4 format
            let num_oids = candidate_num_oids;
            let mut index = Vec::with_capacity(num_oids + 1);
            for _ in 0..=num_oids {
                index.push(cur.read_u64::<LittleEndian>()?);
            }
            let data_offset = 8 + v4_index_bytes;
            let remaining = data.len() - data_offset;
            let num_taxids = remaining / 4;
            let mut taxids = Vec::with_capacity(num_taxids);
            for _ in 0..num_taxids {
                taxids.push(cur.read_i32::<LittleEndian>()?);
            }
            return Ok(TaxIdLookup { index, taxids });
        }

        // V5 format: flat per-OID records (u32_LE count, then count × i32_LE taxids).
        Self::parse_v5_nto(&data)
    }

    /// Parse V5 .nto format: flat per-OID records.
    /// Each record: u32_LE count, then count × i32_LE taxids.
    fn parse_v5_nto(data: &[u8]) -> io::Result<Self> {
        let mut cur = Cursor::new(data);
        let mut index: Vec<u64> = Vec::new();
        let mut taxids: Vec<i32> = Vec::new();

        while (cur.position() as usize) + 4 <= data.len() {
            let count = cur.read_u32::<LittleEndian>()? as usize;
            // Sanity check: count should be small and data must have enough bytes
            if count > 1000 || (cur.position() as usize) + count * 4 > data.len() {
                break;
            }
            for _ in 0..count {
                taxids.push(cur.read_i32::<LittleEndian>()?);
            }
            index.push(taxids.len() as u64);
        }

        Ok(TaxIdLookup { index, taxids })
    }

    /// Get the taxid(s) for a given OID.
    fn get_taxids(&self, oid: u32) -> Vec<i32> {
        let oid = oid as usize;
        if oid >= self.index.len() { return Vec::new(); }
        let end = self.index[oid] as usize;
        let start = if oid == 0 { 0 } else { self.index[oid - 1] as usize };
        if start >= self.taxids.len() || end > self.taxids.len() || start >= end {
            return Vec::new();
        }
        self.taxids[start..end].to_vec()
    }
}

/// Taxonomy name information for a single taxid.
#[derive(Debug, Clone, Default)]
pub struct TaxInfo {
    pub scientific_name: String,
    pub common_name: String,
    pub blast_name: String,
    pub kingdom: String,
}

/// Taxonomy name database parsed from taxdb.bti/taxdb.btd files.
/// These are NCBI-distributed files, not per-database — typically found
/// alongside the BLAST databases or in $BLASTDB.
pub struct TaxNameDb {
    /// Sorted (taxid, offset) pairs from taxdb.bti.
    index: Vec<(i32, u32)>,
    /// Raw data from taxdb.btd.
    data: Vec<u8>,
}

impl TaxNameDb {
    /// Try to load taxdb from a directory. Looks for taxdb.bti + taxdb.btd.
    pub fn open(dir: &Path) -> io::Result<Self> {
        let bti_path = dir.join("taxdb.bti");
        let btd_path = dir.join("taxdb.btd");
        let bti_data = std::fs::read(&bti_path)?;
        let btd_data = std::fs::read(&btd_path)?;

        if bti_data.len() < 24 {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "taxdb.bti too small"));
        }
        let mut cur = Cursor::new(&bti_data);
        let magic = cur.read_u32::<BigEndian>()?;
        if magic != 0x8739 {
            return Err(io::Error::new(io::ErrorKind::InvalidData,
                format!("taxdb.bti bad magic: 0x{:04x}", magic)));
        }
        let count = cur.read_u32::<BigEndian>()? as usize;
        // Skip 4 reserved u32 fields
        for _ in 0..4 {
            let _ = cur.read_u32::<BigEndian>()?;
        }
        // Read index entries: (taxid, offset) pairs, each 8 bytes, big-endian
        let mut index = Vec::with_capacity(count);
        for _ in 0..count {
            let taxid = cur.read_i32::<BigEndian>()?;
            let offset = cur.read_u32::<BigEndian>()?;
            index.push((taxid, offset));
        }

        Ok(TaxNameDb { index, data: btd_data })
    }

    /// Look up taxonomy info by taxid. Uses binary search on sorted index.
    pub fn get_info(&self, taxid: i32) -> Option<TaxInfo> {
        let idx = self.index.binary_search_by_key(&taxid, |&(t, _)| t).ok()?;
        let offset = self.index[idx].1 as usize;
        if offset >= self.data.len() { return None; }

        // Find end of record (next entry's offset, or end of data)
        let end = if idx + 1 < self.index.len() {
            (self.index[idx + 1].1 as usize).min(self.data.len())
        } else {
            self.data.len()
        };

        let record = &self.data[offset..end];
        // Tab-delimited: scientific_name \t common_name \t blast_name \t kingdom
        let s = std::str::from_utf8(record).ok()?;
        let s = s.trim_end_matches('\0');
        let fields: Vec<&str> = s.split('\t').collect();
        Some(TaxInfo {
            scientific_name: fields.first().unwrap_or(&"").to_string(),
            common_name: fields.get(1).unwrap_or(&"").to_string(),
            blast_name: fields.get(2).unwrap_or(&"").to_string(),
            kingdom: fields.get(3).unwrap_or(&"").to_string(),
        })
    }
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
            let alias = crate::db::alias::parse_alias_file(&base_path.with_extension(alias_ext))?;
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

        // Try to load taxonomy data from .nto file (v5 databases)
        let nto_path = base_path.with_extension(
            if db_type == DbType::Nucleotide { "nto" } else { "pto" }
        );
        let tax_lookup = TaxIdLookup::from_file_with_hint(&nto_path, Some(num_oids)).ok();

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
            tax_lookup,
        })
    }

    /// Get the raw sequence bytes for the given OID.
    /// For nucleotide: returns packed 2-bit data from seq_start to amb_start.
    /// For protein: returns raw amino acid codes.
    pub fn get_sequence(&self, oid: u32) -> &[u8] {
        assert!((oid as usize) < self.seq_offsets.len().saturating_sub(1),
            "OID {} out of range (num_oids={})", oid, self.num_oids);
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
        assert!((oid as usize) < self.seq_offsets.len().saturating_sub(1),
            "OID {} out of range (num_oids={})", oid, self.num_oids);
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
        assert!((oid as usize) < self.hdr_offsets.len().saturating_sub(1),
            "OID {} out of range (num_oids={})", oid, self.num_oids);
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

    /// Get the taxid(s) for a given OID. Returns empty vec if no taxonomy data.
    pub fn get_taxids(&self, oid: u32) -> Vec<i32> {
        self.tax_lookup.as_ref().map(|t| t.get_taxids(oid)).unwrap_or_default()
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
            .join("tests/fixtures/seqn/seqn")
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
            .join("tests/fixtures/pombe/pombe");
        if !pombe.with_extension("nin").exists() { return; }
        let db = BlastDb::open(&pombe).unwrap();
        assert_eq!(db.get_accession(0).as_deref(), Some("NC_003421.2"));
    }

    // ---- Validation / error handling tests ----

    fn corrupt_db_path(name: &str) -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
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

    // ---- Tests ported from NCBI seqdb_unit_test.cpp ----

    fn pombe_db_path() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/pombe/pombe")
    }

    fn seqp_db_path() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/seqp/seqp")
    }

    /// Open the pombe nucleotide database and verify basic metadata.
    #[test]
    fn test_open_nucleotide_db() {
        let path = pombe_db_path();
        if !path.with_extension("nin").exists() { return; }
        let db = BlastDb::open(&path).unwrap();
        assert_eq!(db.db_type, DbType::Nucleotide);
        assert!(db.num_oids > 0, "pombe database should have sequences");
        // pombe has small number of chromosomes
        assert!(db.num_oids < 100, "pombe should have < 100 sequences");
    }

    /// Open the protein database and verify sequence count and type.
    #[test]
    fn test_open_protein_db_info() {
        let path = seqp_db_path();
        if !path.with_extension("pin").exists() { return; }
        let db = BlastDb::open(&path).unwrap();
        assert_eq!(db.db_type, DbType::Protein);
        assert_eq!(db.num_oids, 2005);
        assert!(db.total_length > 0);
    }

    /// For every OID, get_seq_len should return a positive value consistent
    /// with the raw bytes returned by get_sequence.
    #[test]
    fn test_sequence_length_consistency() {
        let path = seqp_db_path();
        if !path.with_extension("pin").exists() { return; }
        let db = BlastDb::open(&path).unwrap();
        for oid in 0..db.num_oids {
            let len = db.get_seq_len(oid);
            assert!(len > 0, "OID {} should have positive length", oid);
            let raw = db.get_sequence(oid);
            // For protein, raw length == residue count
            assert_eq!(
                raw.len() as u32, len,
                "OID {} raw bytes ({}) != reported length ({})", oid, raw.len(), len
            );
        }
    }

    /// Verify accessions can be retrieved for known OIDs.
    #[test]
    fn test_accession_retrieval() {
        let db = BlastDb::open(&test_db_path()).unwrap();
        // OID 0 should have a known accession
        let acc = db.get_accession(0);
        assert!(acc.is_some(), "OID 0 should have an accession");
        let acc = acc.unwrap();
        assert!(!acc.is_empty());
        // Accession format: letters + optional underscore + digits
        assert!(acc.chars().next().unwrap().is_ascii_uppercase());
    }

    /// Accessing OID >= num_oids should panic.
    #[test]
    #[should_panic]
    fn test_oid_out_of_bounds() {
        let db = BlastDb::open(&test_db_path()).unwrap();
        let _ = db.get_sequence(db.num_oids);
    }

    /// Iterate all OIDs in a database and verify each sequence is readable.
    #[test]
    fn test_all_sequences_readable() {
        let db = BlastDb::open(&test_db_path()).unwrap();
        for oid in 0..db.num_oids {
            let seq = db.get_sequence(oid);
            assert!(!seq.is_empty(), "OID {} should have non-empty sequence data", oid);
            let len = db.get_seq_len(oid);
            assert!(len > 0, "OID {} should have positive length", oid);
            let hdr = db.get_header(oid);
            assert!(!hdr.is_empty(), "OID {} should have non-empty header", oid);
        }
    }
}

    #[test]
    fn test_open_protein_db() {
        let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/seqp/seqp");
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
            .join("tests/fixtures/seqp/seqp");
        if !path.with_extension("pin").exists() { return; }
        let db = BlastDb::open(&path).unwrap();
        let len = db.get_seq_len(0);
        assert!(len > 0 && len <= db.max_seq_len, "Protein seq length should be valid");
    }

    #[test]
    fn test_taxid_lookup_parse() {
        // Build a .nto file in memory: 3 OIDs
        // OID 0: [6239], OID 1: [9606, 9605], OID 2: [7227]
        use std::io::Write;
        let dir = std::env::temp_dir().join("blast_tax_test");
        let _ = std::fs::create_dir_all(&dir);
        let nto_path = dir.join("test.nto");
        {
            let mut f = std::fs::File::create(&nto_path).unwrap();
            // u64 num_oids = 3
            f.write_all(&3u64.to_le_bytes()).unwrap();
            // (num_oids+1) u64 index entries: [1, 3, 4, 4]
            for v in &[1u64, 3, 4, 4] {
                f.write_all(&v.to_le_bytes()).unwrap();
            }
            // i32 taxid data: [6239, 9606, 9605, 7227]
            for t in &[6239i32, 9606, 9605, 7227] {
                f.write_all(&t.to_le_bytes()).unwrap();
            }
        }
        let lookup = TaxIdLookup::from_file(&nto_path).unwrap();
        assert_eq!(lookup.get_taxids(0), vec![6239]);
        assert_eq!(lookup.get_taxids(1), vec![9606, 9605]);
        assert_eq!(lookup.get_taxids(2), vec![7227]);
        assert_eq!(lookup.get_taxids(3), Vec::<i32>::new()); // out of range
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_taxid_lookup_empty_file() {
        // Truncated .nto file should gracefully return empty
        use std::io::Write;
        let dir = std::env::temp_dir().join("blast_tax_test2");
        let _ = std::fs::create_dir_all(&dir);
        let nto_path = dir.join("test.nto");
        {
            let mut f = std::fs::File::create(&nto_path).unwrap();
            f.write_all(&7u64.to_le_bytes()).unwrap(); // header only
        }
        let lookup = TaxIdLookup::from_file(&nto_path).unwrap();
        assert_eq!(lookup.get_taxids(0), Vec::<i32>::new());
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_tax_name_db() {
        let dir = std::env::temp_dir().join("blast_taxdb_test");
        let _ = std::fs::create_dir_all(&dir);

        // Build taxdb.btd (data file): tab-delimited records
        let record0 = "Caenorhabditis elegans\tnematode\tnematodes\tE";
        let record1 = "Homo sapiens\thuman\tprimates\tE";
        let mut btd = Vec::new();
        btd.extend_from_slice(record0.as_bytes());
        btd.push(0); // null terminator
        let offset1 = btd.len() as u32;
        btd.extend_from_slice(record1.as_bytes());
        btd.push(0);
        std::fs::write(dir.join("taxdb.btd"), &btd).unwrap();

        // Build taxdb.bti (index file): big-endian
        let mut bti = Vec::new();
        bti.extend_from_slice(&0x8739u32.to_be_bytes()); // magic
        bti.extend_from_slice(&2u32.to_be_bytes());      // count = 2
        for _ in 0..4 { bti.extend_from_slice(&0u32.to_be_bytes()); } // reserved
        // Entry 0: taxid=6239, offset=0
        bti.extend_from_slice(&6239i32.to_be_bytes());
        bti.extend_from_slice(&0u32.to_be_bytes());
        // Entry 1: taxid=9606, offset=offset1
        bti.extend_from_slice(&9606i32.to_be_bytes());
        bti.extend_from_slice(&offset1.to_be_bytes());
        std::fs::write(dir.join("taxdb.bti"), &bti).unwrap();

        let tdb = TaxNameDb::open(&dir).unwrap();
        let info = tdb.get_info(6239).unwrap();
        assert_eq!(info.scientific_name, "Caenorhabditis elegans");
        assert_eq!(info.common_name, "nematode");
        assert_eq!(info.blast_name, "nematodes");
        assert_eq!(info.kingdom, "E");

        let info2 = tdb.get_info(9606).unwrap();
        assert_eq!(info2.scientific_name, "Homo sapiens");
        assert_eq!(info2.common_name, "human");

        assert!(tdb.get_info(99999).is_none()); // not found

        let _ = std::fs::remove_dir_all(&dir);
    }
