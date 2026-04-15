//! BLAST database index file reader.
//!
//! Reads .nin/.pin index files, .nsq/.psq sequence files, and .nhr/.phr header files.

use byteorder::{BigEndian, LittleEndian, ReadBytesExt};
use memmap2::{Advice, Mmap, UncheckedAdvice};
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

fn db_path_with_ext(base_path: &Path, ext: &str) -> PathBuf {
    let mut path = base_path.as_os_str().to_os_string();
    path.push(".");
    path.push(ext);
    PathBuf::from(path)
}

/// A BLAST database, which may be a single physical volume or a multi-volume
/// alias database.
pub struct BlastDb {
    pub db_type: DbType,
    pub title: String,
    pub date: String,
    pub num_oids: u32,
    pub total_length: u64,
    pub max_seq_len: u32,
    pub version: u32,

    volumes: Vec<BlastDbVolume>,
    volume_starts: Vec<u32>,
    /// Optional global OID-to-taxid lookup (from alias-level .not/.pot files).
    tax_lookup: Option<TaxIdLookup>,
    /// Offset added to public OIDs before indexing `tax_lookup`.
    tax_lookup_oid_offset: u32,
}

/// A single physical volume of a BLAST database.
struct BlastDbVolume {
    db_type: DbType,
    title: String,
    date: String,
    num_oids: u32,
    total_length: u64,
    max_seq_len: u32,
    version: u32,

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
    /// OID-to-taxid lookup (from .not/.pot files, or legacy test .nto/.pto files).
    tax_lookup: Option<TaxIdLookup>,
}

/// OID-to-taxid lookup table.
enum TaxIdLookup {
    /// In-memory representation used by older/simple .nto files and tests.
    InMemory {
        /// Index: cumulative end offsets into `taxids`.
        index: Vec<u64>,
        /// Flat taxid array, indexed by the `index` array.
        taxids: Vec<i32>,
    },
    /// BLAST DB v5 OID-to-taxids lookup file (.not/.pot), mmap-backed.
    OidToTaxIdsMmap {
        mmap: Mmap,
        num_oids: u64,
        data_offset: usize,
    },
}

impl BlastDbVolume {
    fn advise(&self, advice: Advice) -> io::Result<()> {
        self.seq_mmap.advise(advice)?;
        self.hdr_mmap.advise(advice)?;
        Ok(())
    }

    fn advise_dontneed(&self) -> io::Result<()> {
        // SAFETY: these are read-only file-backed mappings, and callers only
        // use this after completing the scan for the current volume. Future
        // reads remain valid and will fault pages back from the database file.
        unsafe {
            self.seq_mmap.unchecked_advise(UncheckedAdvice::DontNeed)?;
            self.hdr_mmap.unchecked_advise(UncheckedAdvice::DontNeed)?;
        }
        Ok(())
    }

    fn get_sequence(&self, local_oid: u32) -> &[u8] {
        assert!(
            (local_oid as usize) < self.seq_offsets.len().saturating_sub(1),
            "local OID {} out of range (num_oids={})",
            local_oid,
            self.num_oids
        );
        let start = self.seq_offsets[local_oid as usize] as usize;
        let end = if self.db_type == DbType::Nucleotide {
            self.amb_offsets.as_ref().unwrap()[local_oid as usize] as usize
        } else {
            self.seq_offsets[local_oid as usize + 1] as usize
        };
        &self.seq_mmap[start..end]
    }

    fn get_seq_len(&self, local_oid: u32) -> u32 {
        assert!(
            (local_oid as usize) < self.seq_offsets.len().saturating_sub(1),
            "local OID {} out of range (num_oids={})",
            local_oid,
            self.num_oids
        );
        let start = self.seq_offsets[local_oid as usize] as usize;
        if self.db_type == DbType::Nucleotide {
            let end = self.amb_offsets.as_ref().unwrap()[local_oid as usize] as usize;
            let byte_len = end - start;
            if byte_len == 0 {
                return 0;
            }
            let whole_bytes = (byte_len - 1) as u32;
            let last_byte = self.seq_mmap[start + whole_bytes as usize];
            let remainder = (last_byte & 0x03) as u32;
            whole_bytes * 4 + remainder
        } else {
            let end = self.seq_offsets[local_oid as usize + 1] as usize;
            (end - start) as u32
        }
    }
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
        if matches!(
            path.extension().and_then(|s| s.to_str()),
            Some("not" | "pot")
        ) {
            return Self::parse_oid_to_taxids_mmap(path, expected_oids);
        }

        let data = std::fs::read(path)?;
        if data.len() < 8 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "nto file too small",
            ));
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
            return Ok(TaxIdLookup::InMemory { index, taxids });
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

        Ok(TaxIdLookup::InMemory { index, taxids })
    }

    /// Parse BLAST DB v5 OID-to-taxids lookup files (.not/.pot).
    ///
    /// Format from NCBI's blastdbv5_files.txt:
    /// u64 num_oids, then u64[num_oids] cumulative taxid-data ends, then
    /// i32 taxid data. This matches CLookupTaxIds in NCBI's seqdb_lmdb.cpp:
    /// data starts at 8 * (num_oids + 1) bytes from the file start.
    fn parse_oid_to_taxids_mmap(path: &Path, expected_oids: Option<u32>) -> io::Result<Self> {
        let file = File::open(path)?;
        let mmap = unsafe { Mmap::map(&file)? };
        if mmap.len() < 16 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "oid-to-taxids file too small",
            ));
        }

        let num_oids = read_u64_le_at(&mmap, 0)?;
        if num_oids == 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "oid-to-taxids file has zero OIDs",
            ));
        }
        if let Some(expected) = expected_oids {
            if num_oids != expected as u64 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "oid-to-taxids OID count does not match database",
                ));
            }
        }

        let index_bytes = (num_oids as usize).checked_mul(8).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "oid-to-taxids index too large")
        })?;
        let data_offset = 8usize.checked_add(index_bytes).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "oid-to-taxids offset overflow")
        })?;
        if mmap.len() < data_offset {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "oid-to-taxids truncated index",
            ));
        }
        let last_index_offset = 8 + (num_oids as usize - 1) * 8;
        let last_end = read_u64_le_at(&mmap, last_index_offset)? as usize;
        let data_bytes = last_end.checked_mul(4).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "oid-to-taxids data too large")
        })?;
        if mmap.len() < data_offset + data_bytes {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "oid-to-taxids truncated data",
            ));
        }

        Ok(TaxIdLookup::OidToTaxIdsMmap {
            mmap,
            num_oids,
            data_offset,
        })
    }

    /// Get the taxid(s) for a given OID.
    fn get_taxids(&self, oid: u32) -> Vec<i32> {
        match self {
            TaxIdLookup::InMemory { index, taxids } => {
                let oid = oid as usize;
                if oid >= index.len() {
                    return Vec::new();
                }
                let end = index[oid] as usize;
                let start = if oid == 0 { 0 } else { index[oid - 1] as usize };
                if start >= taxids.len() || end > taxids.len() || start >= end {
                    return Vec::new();
                }
                taxids[start..end].to_vec()
            }
            TaxIdLookup::OidToTaxIdsMmap {
                mmap,
                num_oids,
                data_offset,
            } => {
                let oid_u64 = oid as u64;
                if oid_u64 >= *num_oids {
                    return Vec::new();
                }
                let oid = oid as usize;
                let end = match read_u64_le_at(mmap, 8 + oid * 8) {
                    Ok(v) => v as usize,
                    Err(_) => return Vec::new(),
                };
                let start = if oid == 0 {
                    0
                } else {
                    match read_u64_le_at(mmap, 8 + (oid - 1) * 8) {
                        Ok(v) => v as usize,
                        Err(_) => return Vec::new(),
                    }
                };
                if start >= end {
                    return Vec::new();
                }
                let Some(start_byte) = data_offset.checked_add(start.saturating_mul(4)) else {
                    return Vec::new();
                };
                let Some(end_byte) = data_offset.checked_add(end.saturating_mul(4)) else {
                    return Vec::new();
                };
                if end_byte > mmap.len() || start_byte >= end_byte {
                    return Vec::new();
                }
                mmap[start_byte..end_byte]
                    .chunks_exact(4)
                    .map(|bytes| i32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]))
                    .collect()
            }
        }
    }
}

fn read_u64_le_at(buf: &[u8], offset: usize) -> io::Result<u64> {
    let end = offset
        .checked_add(8)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "u64 offset overflow"))?;
    let bytes = buf
        .get(offset..end)
        .ok_or_else(|| io::Error::new(io::ErrorKind::UnexpectedEof, "truncated u64"))?;
    Ok(u64::from_le_bytes([
        bytes[0], bytes[1], bytes[2], bytes[3], bytes[4], bytes[5], bytes[6], bytes[7],
    ]))
}

fn read_index_num_oids(base_path: &Path, db_type: DbType) -> io::Result<u32> {
    let idx_data = std::fs::read(db_path_with_ext(base_path, db_type.index_ext()))?;
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

    if version == 5 {
        let _volume_number = cur.read_u32::<BigEndian>()?;
    }
    let _title = read_string(&mut cur)?;
    if version == 5 {
        let _lmdb_name = read_string(&mut cur)?;
    }
    let _date = read_string(&mut cur)?;
    cur.read_u32::<BigEndian>()
}

fn alias_base_and_volume_start_for_member(
    base_path: &Path,
    db_type: DbType,
) -> io::Result<Option<(PathBuf, u64)>> {
    let Some(name) = base_path.file_name().and_then(|s| s.to_str()) else {
        return Ok(None);
    };
    if name.len() <= 3 || name.as_bytes()[name.len() - 3] != b'.' {
        return Ok(None);
    }
    let suffix = &name[name.len() - 2..];
    if !suffix.as_bytes().iter().all(u8::is_ascii_digit) {
        return Ok(None);
    }

    let alias_name = &name[..name.len() - 3];
    let alias_base = base_path.with_file_name(alias_name);
    let alias_path = if alias_base.with_extension("nal").exists() {
        alias_base.with_extension("nal")
    } else if alias_base.with_extension("pal").exists() {
        alias_base.with_extension("pal")
    } else {
        return Ok(None);
    };
    let alias = crate::db::alias::parse_alias_file(&alias_path)?;
    let canonical_target = std::fs::canonicalize(base_path).ok();
    let mut volume_start = 0u64;

    for volume_base in alias.dblist {
        let is_target = if let Some(target) = &canonical_target {
            std::fs::canonicalize(&volume_base)
                .map(|p| p == *target)
                .unwrap_or(false)
        } else {
            volume_base.file_name() == base_path.file_name()
        };
        if is_target {
            return Ok(Some((alias_base, volume_start)));
        }
        volume_start += read_index_num_oids(&volume_base, db_type)? as u64;
    }

    Ok(None)
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
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "taxdb.bti too small",
            ));
        }
        let mut cur = Cursor::new(&bti_data);
        let magic = cur.read_u32::<BigEndian>()?;
        if magic != 0x8739 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("taxdb.bti bad magic: 0x{:04x}", magic),
            ));
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

        Ok(TaxNameDb {
            index,
            data: btd_data,
        })
    }

    /// Look up taxonomy info by taxid. Uses binary search on sorted index.
    pub fn get_info(&self, taxid: i32) -> Option<TaxInfo> {
        let idx = self.index.binary_search_by_key(&taxid, |&(t, _)| t).ok()?;
        let offset = self.index[idx].1 as usize;
        if offset >= self.data.len() {
            return None;
        }

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
    /// Global OID ranges for each physical volume in database order.
    pub fn volume_oid_ranges(&self) -> Vec<(u32, u32)> {
        self.volumes
            .iter()
            .enumerate()
            .map(|(idx, vol)| {
                let start = self.volume_starts[idx];
                (start, start + vol.num_oids)
            })
            .collect()
    }

    /// Get sequence bytes and residue length by physical volume and local OID.
    ///
    /// This is intended for sequential volume scans that already know the
    /// active volume and need to avoid resolving every global OID in the hot
    /// loop.
    pub fn get_volume_sequence_and_len(&self, volume_idx: usize, local_oid: u32) -> (&[u8], u32) {
        let vol = &self.volumes[volume_idx];
        (vol.get_sequence(local_oid), vol.get_seq_len(local_oid))
    }

    /// Hint that a volume is about to be scanned sequentially.
    pub fn advise_volume_sequential(&self, volume_idx: usize) {
        let _ = self.volumes[volume_idx].advise(Advice::Sequential);
    }

    /// Hint that pages for a scanned volume can be dropped from resident memory.
    pub fn advise_volume_dontneed(&self, volume_idx: usize) {
        let _ = self.volumes[volume_idx].advise_dontneed();
    }

    /// Open a BLAST database from the given base path (without extension).
    /// Automatically detects nucleotide vs protein.
    pub fn open(base_path: &Path) -> io::Result<Self> {
        // Try nucleotide first, then protein
        let db_type = if db_path_with_ext(base_path, "nin").exists() {
            DbType::Nucleotide
        } else if db_path_with_ext(base_path, "pin").exists() {
            DbType::Protein
        } else if base_path.with_extension("nal").exists()
            || base_path.with_extension("pal").exists()
        {
            let alias_ext = if base_path.with_extension("nal").exists() {
                "nal"
            } else {
                "pal"
            };
            let alias = crate::db::alias::parse_alias_file(&base_path.with_extension(alias_ext))?;
            return Self::open_alias(base_path, alias);
        } else {
            return Err(io::Error::new(
                io::ErrorKind::NotFound,
                format!("No BLAST database found at {}", base_path.display()),
            ));
        };
        Self::open_typed(base_path, db_type)
    }

    fn open_alias(base_path: &Path, alias: crate::db::alias::AliasFile) -> io::Result<Self> {
        if alias.dblist.is_empty() {
            return Err(io::Error::new(io::ErrorKind::NotFound, "Empty alias file"));
        }

        let mut volumes = Vec::with_capacity(alias.dblist.len());
        for (idx, vol_path) in alias.dblist.iter().enumerate() {
            match Self::open_volume_auto(vol_path) {
                Ok(volume) => volumes.push(volume),
                Err(err) if err.kind() == io::ErrorKind::NotFound => {
                    let volume_name = alias
                        .raw_dblist
                        .get(idx)
                        .map(String::as_str)
                        .unwrap_or_else(|| vol_path.to_str().unwrap_or_default());
                    return Err(io::Error::new(
                        io::ErrorKind::NotFound,
                        format!(
                            "Could not find volume or alias file ({}) referenced in alias file ({}).",
                            volume_name,
                            base_path.display()
                        ),
                    ));
                }
                Err(err) => return Err(err),
            }
        }

        let db_type = volumes[0].db_type;
        for vol in &volumes {
            if vol.db_type != db_type {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Alias file references mixed nucleotide/protein volumes",
                ));
            }
        }

        Self::from_volumes(base_path, alias.title, volumes)
    }

    /// Open a BLAST database of known type.
    pub fn open_typed(base_path: &Path, db_type: DbType) -> io::Result<Self> {
        let volume = Self::open_volume_typed(base_path, db_type)?;
        Self::from_volumes(base_path, None, vec![volume])
    }

    fn open_volume_auto(base_path: &Path) -> io::Result<BlastDbVolume> {
        let db_type = if db_path_with_ext(base_path, "nin").exists() {
            DbType::Nucleotide
        } else if db_path_with_ext(base_path, "pin").exists() {
            DbType::Protein
        } else {
            return Err(io::Error::new(
                io::ErrorKind::NotFound,
                format!("No BLAST database volume found at {}", base_path.display()),
            ));
        };
        Self::open_volume_typed(base_path, db_type)
    }

    fn open_volume_typed(base_path: &Path, db_type: DbType) -> io::Result<BlastDbVolume> {
        let idx_path = db_path_with_ext(base_path, db_type.index_ext());
        let seq_path = db_path_with_ext(base_path, db_type.seq_ext());
        let hdr_path = db_path_with_ext(base_path, db_type.hdr_ext());

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

        // Try to load direct OID-to-taxid data. Modern v5 databases use
        // .not/.pot; older/simple test databases may still use .nto/.pto.
        let tax_lookup = {
            let oid_path = db_path_with_ext(
                base_path,
                if db_type == DbType::Nucleotide {
                    "not"
                } else {
                    "pot"
                },
            );
            TaxIdLookup::from_file_with_hint(&oid_path, Some(num_oids))
                .ok()
                .or_else(|| {
                    let legacy_path = db_path_with_ext(
                        base_path,
                        if db_type == DbType::Nucleotide {
                            "nto"
                        } else {
                            "pto"
                        },
                    );
                    TaxIdLookup::from_file_with_hint(&legacy_path, Some(num_oids)).ok()
                })
        };

        Ok(BlastDbVolume {
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

    fn from_volumes(
        base_path: &Path,
        alias_title: Option<String>,
        volumes: Vec<BlastDbVolume>,
    ) -> io::Result<Self> {
        if volumes.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::NotFound,
                "No database volumes",
            ));
        }

        let db_type = volumes[0].db_type;
        let version = volumes
            .iter()
            .map(|v| v.version)
            .max()
            .unwrap_or(volumes[0].version);
        let title = alias_title.unwrap_or_else(|| volumes[0].title.clone());
        let date = volumes[0].date.clone();
        let max_seq_len = volumes.iter().map(|v| v.max_seq_len).max().unwrap_or(0);
        let total_length = volumes.iter().map(|v| v.total_length).sum();

        let mut num_oids_u64 = 0u64;
        let mut volume_starts = Vec::with_capacity(volumes.len());
        for vol in &volumes {
            volume_starts.push(u32::try_from(num_oids_u64).map_err(|_| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Database has more than u32::MAX OIDs",
                )
            })?);
            num_oids_u64 += vol.num_oids as u64;
        }
        let num_oids = u32::try_from(num_oids_u64).map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "Database has more than u32::MAX OIDs",
            )
        })?;

        // Single-volume databases commonly have small per-volume taxid files.
        // Large alias databases such as core_nt have a global .not that can be
        // over a GB; loading it eagerly makes opening the DB expensive even
        // when taxonomy fields are not requested.
        let tax_lookup = if volumes.len() == 1 {
            let oid_path = db_path_with_ext(
                base_path,
                if db_type == DbType::Nucleotide {
                    "not"
                } else {
                    "pot"
                },
            );
            TaxIdLookup::from_file_with_hint(&oid_path, Some(num_oids))
                .ok()
                .or_else(|| {
                    let legacy_path = db_path_with_ext(
                        base_path,
                        if db_type == DbType::Nucleotide {
                            "nto"
                        } else {
                            "pto"
                        },
                    );
                    TaxIdLookup::from_file_with_hint(&legacy_path, Some(num_oids)).ok()
                })
        } else {
            None
        };

        Ok(BlastDb {
            db_type,
            title,
            date,
            num_oids,
            total_length,
            max_seq_len,
            version,
            volumes,
            volume_starts,
            tax_lookup,
            tax_lookup_oid_offset: 0,
        })
    }

    fn resolve_oid(&self, oid: u32) -> (&BlastDbVolume, u32) {
        assert!(
            oid < self.num_oids,
            "OID {} out of range (num_oids={})",
            oid,
            self.num_oids
        );
        let idx = match self.volume_starts.binary_search(&oid) {
            Ok(i) => i,
            Err(0) => 0,
            Err(i) => i - 1,
        };
        let local_oid = oid - self.volume_starts[idx];
        (&self.volumes[idx], local_oid)
    }

    /// Get the raw sequence bytes for the given OID.
    /// For nucleotide: returns packed 2-bit data from seq_start to amb_start.
    /// For protein: returns raw amino acid codes.
    pub fn get_sequence(&self, oid: u32) -> &[u8] {
        let (vol, local_oid) = self.resolve_oid(oid);
        vol.get_sequence(local_oid)
    }

    /// Get the sequence length in residues for the given OID.
    pub fn get_seq_len(&self, oid: u32) -> u32 {
        let (vol, local_oid) = self.resolve_oid(oid);
        vol.get_seq_len(local_oid)
    }

    /// Get the raw header bytes (ASN.1) for the given OID.
    pub fn get_header(&self, oid: u32) -> &[u8] {
        let (vol, local_oid) = self.resolve_oid(oid);
        assert!(
            (local_oid as usize) < vol.hdr_offsets.len().saturating_sub(1),
            "OID {} out of range (num_oids={})",
            oid,
            self.num_oids
        );
        let start = vol.hdr_offsets[local_oid as usize] as usize;
        let end = vol.hdr_offsets[local_oid as usize + 1] as usize;
        &vol.hdr_mmap[start..end]
    }

    /// Extract the best accession-like string from the ASN.1 header.
    /// Prefer versioned Seq-id candidates over title text, because deflines
    /// often contain gene/locus names before the real accession object.
    pub fn get_accession(&self, oid: u32) -> Option<String> {
        let hdr = self.get_header(oid);
        extract_accession_from_header(hdr)
    }

    /// Extract a BLAST-style subject defline from the ASN.1 header.
    /// This keeps the accession used by tabular output as the leading token,
    /// followed by the human-readable title when one is present.
    pub fn get_defline(&self, oid: u32) -> Option<String> {
        let hdr = self.get_header(oid);
        let accession = extract_accession_from_header(hdr);
        let title = extract_title_from_header(hdr);
        match (accession, title) {
            (Some(acc), Some(title)) if title == acc || title.starts_with(&format!("{} ", acc)) => {
                Some(title)
            }
            (Some(acc), Some(title)) => Some(format!("{} {}", acc, title)),
            (Some(acc), None) => Some(acc),
            (None, Some(title)) => Some(title),
            (None, None) => None,
        }
    }

    /// Get the taxid(s) for a given OID. Returns empty vec if no taxonomy data.
    pub fn get_taxids(&self, oid: u32) -> Vec<i32> {
        if let Some(tax_lookup) = &self.tax_lookup {
            let Some(global_oid) = oid.checked_add(self.tax_lookup_oid_offset) else {
                return Vec::new();
            };
            return tax_lookup.get_taxids(global_oid);
        }
        let (vol, local_oid) = self.resolve_oid(oid);
        vol.tax_lookup
            .as_ref()
            .map(|t| t.get_taxids(local_oid))
            .unwrap_or_default()
    }

    /// Load a database-level OID-to-taxid lookup from the matching .not/.pot
    /// file. Large alias databases such as core_nt keep taxonomy at the alias
    /// level, so this is called only when taxonomy output fields are requested.
    pub fn load_tax_lookup_from_base_path(&mut self, base_path: &Path) -> io::Result<bool> {
        if self.tax_lookup.is_some() {
            return Ok(true);
        }

        let oid_ext = if self.db_type == DbType::Nucleotide {
            "not"
        } else {
            "pot"
        };
        let legacy_ext = if self.db_type == DbType::Nucleotide {
            "nto"
        } else {
            "pto"
        };
        let mut candidates = vec![
            (db_path_with_ext(base_path, oid_ext), Some(self.num_oids), 0),
            (
                db_path_with_ext(base_path, legacy_ext),
                Some(self.num_oids),
                0,
            ),
        ];

        // NCBI can resolve taxonomy for an individual core_nt volume from the
        // alias-level core_nt.not. For a standalone volume path ending in
        // ".NN", compute that volume's global OID start from the alias DBLIST
        // and use it as an indexing offset into the alias lookup.
        if let Some((alias_base, volume_start)) =
            alias_base_and_volume_start_for_member(base_path, self.db_type)?
        {
            if volume_start <= u32::MAX as u64 {
                candidates.push((
                    db_path_with_ext(&alias_base, oid_ext),
                    None,
                    volume_start as u32,
                ));
            }
        }

        for (path, expected_oids, oid_offset) in candidates {
            match TaxIdLookup::from_file_with_hint(&path, expected_oids) {
                Ok(lookup) => {
                    self.tax_lookup = Some(lookup);
                    self.tax_lookup_oid_offset = oid_offset;
                    return Ok(true);
                }
                Err(err) if err.kind() == io::ErrorKind::NotFound => {}
                Err(err) => return Err(err),
            }
        }

        Ok(false)
    }

    /// Get the ambiguity data for a nucleotide OID.
    /// Ambiguity data is stored between amb_offset and the next sequence's start.
    pub fn get_ambiguity_data(&self, oid: u32) -> Option<&[u8]> {
        let (vol, local_oid) = self.resolve_oid(oid);
        let amb_offsets = vol.amb_offsets.as_ref()?;
        let start = amb_offsets[local_oid as usize] as usize;
        // Ambiguity data ends at the next sequence's start offset
        let end = vol.seq_offsets[local_oid as usize + 1] as usize;
        if start == end {
            None
        } else {
            Some(&vol.seq_mmap[start..end])
        }
    }
}

fn extract_accession_from_header(hdr: &[u8]) -> Option<String> {
    let mut i = 0;
    let mut first_text_versioned = None;
    let mut first_unversioned = None;
    let mut first_local = None;
    while i < hdr.len() {
        if hdr[i] == 0x1a && i + 1 < hdr.len() {
            let len = hdr[i + 1] as usize;
            let start = i + 2;
            let end = start.saturating_add(len);
            if len > 0
                && end <= hdr.len()
                && hdr[start..end]
                    .iter()
                    .all(|&b| b.is_ascii_alphanumeric() || matches!(b, b'_' | b'.' | b'-' | b'|'))
                && first_local.is_none()
            {
                first_local = Some(String::from_utf8_lossy(&hdr[start..end]).to_string());
            }
        }

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
                        let dot_start = i;
                        i += 1;
                        while i < hdr.len() && hdr[i].is_ascii_digit() {
                            i += 1;
                        }
                        if first_text_versioned.is_none() {
                            first_text_versioned =
                                Some(String::from_utf8_lossy(&hdr[start..i]).to_string());
                        }
                        let mut vi = i;
                        while vi < hdr.len() && hdr[vi] == 0 {
                            vi += 1;
                        }
                        if vi + 4 < hdr.len()
                            && hdr[vi] == 0xa3
                            && hdr[vi + 1] == 0x80
                            && hdr[vi + 2] == 0x02
                            && hdr[vi + 3] == 0x01
                        {
                            let version = hdr[vi + 4];
                            return Some(format!(
                                "{}.{}",
                                String::from_utf8_lossy(&hdr[start..dot_start]),
                                version
                            ));
                        }
                        continue;
                    }
                    let acc = String::from_utf8_lossy(&hdr[start..i]).to_string();
                    // Look for ASN.1 version: \x00+\xa3\x80\x02\x01\xNN
                    let mut vi = i;
                    while vi < hdr.len() && hdr[vi] == 0 {
                        vi += 1;
                    }
                    if vi + 4 < hdr.len()
                        && hdr[vi] == 0xa3
                        && hdr[vi + 1] == 0x80
                        && hdr[vi + 2] == 0x02
                        && hdr[vi + 3] == 0x01
                    {
                        let version = hdr[vi + 4];
                        return Some(format!("{}.{}", acc, version));
                    }
                    if first_unversioned.is_none() {
                        first_unversioned = Some(acc);
                    }
                }
            }
        } else {
            i += 1;
        }
    }
    first_text_versioned.or(first_unversioned).or(first_local)
}

fn extract_title_from_header(hdr: &[u8]) -> Option<String> {
    let mut i = 0;
    while i + 1 < hdr.len() {
        if matches!(hdr[i], 0x1a | 0x0c) {
            if let Some((len, len_len)) = read_ber_len(hdr, i + 1) {
                let start = i + 1 + len_len;
                let end = start.saturating_add(len);
                if len > 0 && end <= hdr.len() {
                    let bytes = &hdr[start..end];
                    if bytes
                        .iter()
                        .all(|&b| b == b'\t' || b == b' ' || b.is_ascii_graphic())
                    {
                        let s = String::from_utf8_lossy(bytes).trim().to_string();
                        if looks_like_title(&s) {
                            return Some(s);
                        }
                    }
                    i = end;
                    continue;
                }
            }
        }
        i += 1;
    }
    None
}

fn read_ber_len(buf: &[u8], pos: usize) -> Option<(usize, usize)> {
    let first = *buf.get(pos)?;
    if first & 0x80 == 0 {
        return Some((first as usize, 1));
    }
    let count = (first & 0x7f) as usize;
    if count == 0 || count > std::mem::size_of::<usize>() || pos + 1 + count > buf.len() {
        return None;
    }
    let mut len = 0usize;
    for &b in &buf[pos + 1..pos + 1 + count] {
        len = (len << 8) | b as usize;
    }
    Some((len, 1 + count))
}

fn looks_like_title(s: &str) -> bool {
    if s.is_empty() {
        return false;
    }
    let has_word_separator = s.bytes().any(|b| b == b' ' || b == b'\t');
    let has_lowercase = s.bytes().any(|b| b.is_ascii_lowercase());
    let not_seqid_text = s
        .bytes()
        .all(|b| b.is_ascii_alphanumeric() || matches!(b, b'_' | b'.' | b'-' | b'|'));
    (has_word_separator || has_lowercase) && !not_seqid_text
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
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/seqn/seqn")
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
        assert_eq!(
            total, db.total_length,
            "Sum of individual lengths must match reported total_length"
        );
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
        let pombe = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/pombe/pombe");
        if !pombe.with_extension("nin").exists() {
            return;
        }
        let db = BlastDb::open(&pombe).unwrap();
        assert_eq!(db.get_accession(0).as_deref(), Some("NC_003421.2"));
    }

    #[test]
    fn test_accession_local_id_fallback() {
        let hdr = [
            0x30, 0x80, 0x30, 0x80, 0xa0, 0x80, 0x1a, 0x00, 0xa1, 0x80, 0x30, 0x80, 0xa0, 0x80,
            0xa1, 0x80, 0x1a, 0x03, b's', b'1', b'9', 0x00,
        ];
        assert_eq!(extract_accession_from_header(&hdr).as_deref(), Some("s19"));
    }

    #[test]
    fn test_accession_prefers_versioned_refseq_over_local_id() {
        let hdr = [
            0x1a, 0x03, b's', b'1', b'9', b'X', b'R', b'_', b'1', b'5', b'4', b'0', b'5', b'2',
            0x00, 0x00, 0xa3, 0x80, 0x02, 0x01, 0x02,
        ];
        assert_eq!(
            extract_accession_from_header(&hdr).as_deref(),
            Some("XR_154052.2")
        );
    }

    #[test]
    fn test_accession_prefers_asn_version_over_title_token() {
        let hdr = [
            b'F', b'L', b'H', b'1', b'7', b'8', b'0', b'3', b'7', b'.', b'0', b'1', b'X', b';',
            b' ', b'c', b'l', b'o', b'n', b'e', b' ', b'D', b'Q', b'8', b'9', b'1', b'5', b'5',
            b'7', 0x00, 0x00, 0xa3, 0x80, 0x02, 0x01, 0x02,
        ];
        assert_eq!(
            extract_accession_from_header(&hdr).as_deref(),
            Some("DQ891557.2")
        );
    }

    #[test]
    fn test_extract_title_from_visible_string_header() {
        let title = b"Eukaryotic synthetic construct chromosome 13";
        let mut hdr = vec![
            0x30,
            0x80,
            0xa0,
            title.len() as u8 + 2,
            0x1a,
            title.len() as u8,
        ];
        hdr.extend_from_slice(title);
        hdr.extend_from_slice(b"CP034491");

        assert_eq!(
            extract_title_from_header(&hdr).as_deref(),
            Some("Eukaryotic synthetic construct chromosome 13")
        );
    }

    #[test]
    fn test_extract_title_skips_seqid_visible_string() {
        let hdr = [0x1a, 0x08, b'C', b'P', b'0', b'3', b'4', b'4', b'9', b'1'];
        assert_eq!(extract_title_from_header(&hdr), None);
    }

    #[test]
    fn test_open_multi_volume_alias() {
        let seqn = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/seqn/seqn");
        let pombe = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/pombe/pombe");
        if !seqn.with_extension("nin").exists() || !pombe.with_extension("nin").exists() {
            return;
        }

        let tmp = tempfile::tempdir().unwrap();
        let alias_path = tmp.path().join("multi.nal");
        std::fs::write(
            &alias_path,
            format!(
                "TITLE multi fixture\nDBLIST {} {}\n",
                seqn.display(),
                pombe.display(),
            ),
        )
        .unwrap();

        let seqn_db = BlastDb::open(&seqn).unwrap();
        let pombe_db = BlastDb::open(&pombe).unwrap();
        let db = BlastDb::open(&tmp.path().join("multi")).unwrap();

        assert_eq!(db.db_type, DbType::Nucleotide);
        assert_eq!(db.title, "multi fixture");
        assert_eq!(db.num_oids, seqn_db.num_oids + pombe_db.num_oids);
        assert_eq!(
            db.total_length,
            seqn_db.total_length + pombe_db.total_length
        );
        assert_eq!(db.get_seq_len(0), seqn_db.get_seq_len(0));

        let first_pombe_oid = seqn_db.num_oids;
        assert_eq!(db.get_seq_len(first_pombe_oid), pombe_db.get_seq_len(0));
        assert_eq!(db.get_accession(first_pombe_oid), pombe_db.get_accession(0));
        assert_eq!(
            db.volume_oid_ranges(),
            vec![(0, seqn_db.num_oids), (seqn_db.num_oids, db.num_oids)]
        );
    }

    #[test]
    #[ignore]
    fn test_core_nt_multivolume_metadata_if_available() {
        let core_nt = PathBuf::from("/husky/henriksson/for_claude/blast/core_nt/core_nt");
        if !core_nt.with_extension("nal").exists() {
            eprintln!("Skipping: core_nt alias not found at {:?}", core_nt);
            return;
        }

        let db = BlastDb::open(&core_nt).expect("core_nt multi-volume alias should open");
        assert_eq!(db.db_type, DbType::Nucleotide);
        assert_eq!(db.version, 5);
        assert_eq!(db.num_oids, 124_278_414);
        assert_eq!(db.total_length, 991_049_906_671);
        assert_eq!(db.max_seq_len, 99_993_016);
        assert_eq!(db.get_seq_len(0), 2440);
        assert!(!db.get_header(0).is_empty());

        let last_oid = db.num_oids - 1;
        assert!(db.get_seq_len(last_oid) > 0);
        assert!(!db.get_sequence(last_oid).is_empty());
    }

    #[test]
    #[ignore]
    fn test_core_nt_refseq_accession_preferred_over_locus_title_if_available() {
        let core_nt_00 = PathBuf::from("/husky/henriksson/for_claude/blast/core_nt/core_nt.00");
        if !db_path_with_ext(&core_nt_00, "nin").exists() {
            eprintln!("Skipping: core_nt.00 volume not found at {:?}", core_nt_00);
            return;
        }

        let db = BlastDb::open(&core_nt_00).expect("core_nt.00 should open");
        assert_eq!(db.get_accession(761_984).as_deref(), Some("XR_154052.2"));
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
        assert!(
            amb.len() >= 8,
            "Ambiguity data should have header + at least 1 record"
        );
    }

    #[test]
    fn test_header_data() {
        let db = BlastDb::open(&test_db_path()).unwrap();
        let hdr = db.get_header(0);
        assert!(!hdr.is_empty(), "Header should not be empty");
        // Should contain the title string somewhere
        assert!(
            hdr.windows(8).any(|w| w == b"Xenopus " || w == b"BP722512"),
            "Header should contain sequence info"
        );
    }

    // ---- Tests ported from NCBI seqdb_unit_test.cpp ----

    fn pombe_db_path() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/pombe/pombe")
    }

    fn seqp_db_path() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/seqp/seqp")
    }

    /// Open the pombe nucleotide database and verify basic metadata.
    #[test]
    fn test_open_nucleotide_db() {
        let path = pombe_db_path();
        if !path.with_extension("nin").exists() {
            return;
        }
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
        if !path.with_extension("pin").exists() {
            return;
        }
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
        if !path.with_extension("pin").exists() {
            return;
        }
        let db = BlastDb::open(&path).unwrap();
        for oid in 0..db.num_oids {
            let len = db.get_seq_len(oid);
            assert!(len > 0, "OID {} should have positive length", oid);
            let raw = db.get_sequence(oid);
            // For protein, raw length == residue count
            assert_eq!(
                raw.len() as u32,
                len,
                "OID {} raw bytes ({}) != reported length ({})",
                oid,
                raw.len(),
                len
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
            assert!(
                !seq.is_empty(),
                "OID {} should have non-empty sequence data",
                oid
            );
            let len = db.get_seq_len(oid);
            assert!(len > 0, "OID {} should have positive length", oid);
            let hdr = db.get_header(oid);
            assert!(!hdr.is_empty(), "OID {} should have non-empty header", oid);
        }
    }
}

#[test]
fn test_open_protein_db() {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/seqp/seqp");
    if !path.with_extension("pin").exists() {
        return;
    }
    let db = BlastDb::open(&path).unwrap();
    assert_eq!(db.db_type, DbType::Protein);
    assert_eq!(db.version, 4);
    assert_eq!(db.num_oids, 2005);
    assert!(db.total_length > 0);
}

#[test]
fn test_protein_seq_length() {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/seqp/seqp");
    if !path.with_extension("pin").exists() {
        return;
    }
    let db = BlastDb::open(&path).unwrap();
    let len = db.get_seq_len(0);
    assert!(
        len > 0 && len <= db.max_seq_len,
        "Protein seq length should be valid"
    );
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
fn test_oid_to_taxids_mmap_parse() {
    // Build a BLAST DB v5 .not file:
    // u64 num_oids, num_oids cumulative end offsets, then i32 taxids.
    use std::io::Write;
    let dir = std::env::temp_dir().join("blast_tax_not_test");
    let _ = std::fs::create_dir_all(&dir);
    let not_path = dir.join("test.not");
    {
        let mut f = std::fs::File::create(&not_path).unwrap();
        f.write_all(&3u64.to_le_bytes()).unwrap();
        for v in &[1u64, 3, 4] {
            f.write_all(&v.to_le_bytes()).unwrap();
        }
        for t in &[6239i32, 9606, 9605, 7227] {
            f.write_all(&t.to_le_bytes()).unwrap();
        }
    }
    let lookup = TaxIdLookup::from_file_with_hint(&not_path, Some(3)).unwrap();
    assert_eq!(lookup.get_taxids(0), vec![6239]);
    assert_eq!(lookup.get_taxids(1), vec![9606, 9605]);
    assert_eq!(lookup.get_taxids(2), vec![7227]);
    assert_eq!(lookup.get_taxids(3), Vec::<i32>::new());
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
    bti.extend_from_slice(&2u32.to_be_bytes()); // count = 2
    for _ in 0..4 {
        bti.extend_from_slice(&0u32.to_be_bytes());
    } // reserved
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
