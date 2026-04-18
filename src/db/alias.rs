//! BLAST database alias file (.nal/.pal) support.
//! Alias files list component database volumes for multi-volume databases.

use std::io::{self, BufRead, BufReader};
use std::path::{Path, PathBuf};

/// Parsed alias file.
#[derive(Debug)]
pub struct AliasFile {
    pub title: Option<String>,
    pub dblist: Vec<PathBuf>,
    pub raw_dblist: Vec<String>,
    pub nseq: Option<u64>,
    pub length: Option<u64>,
    pub stats_nseq: Option<u64>,
    pub stats_total_length: Option<u64>,
    /// One-based inclusive OID range used by BLAST alias files.
    pub first_oid: Option<u32>,
    pub last_oid: Option<u32>,
    pub oidlist: Option<PathBuf>,
    pub raw_oidlist: Option<String>,
    /// NCBI `MEMB_BIT` — which membership bit in the `.nin` v5 header to
    /// consult when filtering OIDs. Parsed for forward compatibility; the
    /// actual filtering is not yet applied at the DB-reader level.
    pub memb_bit: Option<u32>,
    /// NCBI `GILIST` — path to a text or binary GI list file.
    /// Parsed but not yet applied (no GI-based OID filtering).
    pub gilist: Option<PathBuf>,
    pub raw_gilist: Option<String>,
    /// NCBI `TILIST` — path to a Trace ID list file. Parsed only.
    pub tilist: Option<PathBuf>,
    pub raw_tilist: Option<String>,
    /// NCBI `SILIST` — path to a seqid list file. Parsed only.
    pub silist: Option<PathBuf>,
    pub raw_silist: Option<String>,
    /// NCBI `TAXIDLIST` — path to a taxonomy ID list file. Parsed only.
    pub taxidlist: Option<PathBuf>,
    pub raw_taxidlist: Option<String>,
}

impl AliasFile {
    /// List of unsupported-filter warnings for populated fields. Callers
    /// (e.g. CLI) should emit these to stderr when the user opens a
    /// database whose alias specifies filters we parse but don't yet apply.
    /// Matches NCBI's "supported but not applied" style where possible.
    pub fn unsupported_filter_warnings(&self) -> Vec<String> {
        let mut warnings = Vec::new();
        if let Some(n) = self.memb_bit {
            warnings.push(format!(
                "Warning: alias MEMB_BIT={} is parsed but membership-bit filtering is not yet applied",
                n
            ));
        }
        if let Some(p) = &self.raw_gilist {
            warnings.push(format!(
                "Warning: alias GILIST {} is parsed but GI-based filtering is not yet applied",
                p
            ));
        }
        if let Some(p) = &self.raw_tilist {
            warnings.push(format!(
                "Warning: alias TILIST {} is parsed but trace-id filtering is not yet applied",
                p
            ));
        }
        if let Some(p) = &self.raw_silist {
            warnings.push(format!(
                "Warning: alias SILIST {} is parsed but seqid-list filtering is not yet applied",
                p
            ));
        }
        if let Some(p) = &self.raw_taxidlist {
            warnings.push(format!(
                "Warning: alias TAXIDLIST {} is parsed but taxid-based filtering is not yet applied",
                p
            ));
        }
        warnings
    }
}

/// Parse a .nal or .pal alias file.
pub fn parse_alias_file(path: &Path) -> io::Result<AliasFile> {
    let file = std::fs::File::open(path)?;
    let reader = BufReader::new(file);
    let dir = path.parent().unwrap_or(Path::new("."));

    let mut alias = AliasFile {
        title: None,
        dblist: Vec::new(),
        raw_dblist: Vec::new(),
        nseq: None,
        length: None,
        stats_nseq: None,
        stats_total_length: None,
        first_oid: None,
        last_oid: None,
        oidlist: None,
        raw_oidlist: None,
        memb_bit: None,
        gilist: None,
        raw_gilist: None,
        tilist: None,
        raw_tilist: None,
        silist: None,
        raw_silist: None,
        taxidlist: None,
        raw_taxidlist: None,
    };

    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        if let Some(rest) = line.strip_prefix("TITLE ") {
            alias.title = Some(rest.to_string());
        } else if let Some(rest) = line.strip_prefix("DBLIST ") {
            alias.raw_dblist = rest.split_whitespace().map(str::to_string).collect();
            alias.dblist = alias.raw_dblist.iter().map(|name| dir.join(name)).collect();
        } else if let Some(rest) = line.strip_prefix("NSEQ ") {
            alias.nseq = rest.trim().parse().ok();
        } else if let Some(rest) = line.strip_prefix("LENGTH ") {
            alias.length = rest.trim().parse().ok();
        } else if let Some(rest) = line.strip_prefix("STATS_NSEQ ") {
            alias.stats_nseq = rest.trim().parse().ok();
        } else if let Some(rest) = line.strip_prefix("STATS_TOTLEN ") {
            alias.stats_total_length = rest.trim().parse().ok();
        } else if let Some(rest) = line.strip_prefix("FIRST_OID ") {
            alias.first_oid = rest.trim().parse().ok();
        } else if let Some(rest) = line.strip_prefix("LAST_OID ") {
            alias.last_oid = rest.trim().parse().ok();
        } else if let Some(rest) = line.strip_prefix("OIDLIST ") {
            let oidlist = rest.trim();
            if !oidlist.is_empty() {
                alias.raw_oidlist = Some(oidlist.to_string());
                alias.oidlist = Some(dir.join(oidlist));
            }
        } else if let Some(rest) = line.strip_prefix("MEMB_BIT ") {
            alias.memb_bit = rest.trim().parse().ok();
        } else if let Some(rest) = line.strip_prefix("GILIST ") {
            let v = rest.trim();
            if !v.is_empty() {
                alias.raw_gilist = Some(v.to_string());
                alias.gilist = Some(dir.join(v));
            }
        } else if let Some(rest) = line.strip_prefix("TILIST ") {
            let v = rest.trim();
            if !v.is_empty() {
                alias.raw_tilist = Some(v.to_string());
                alias.tilist = Some(dir.join(v));
            }
        } else if let Some(rest) = line.strip_prefix("SILIST ") {
            let v = rest.trim();
            if !v.is_empty() {
                alias.raw_silist = Some(v.to_string());
                alias.silist = Some(dir.join(v));
            }
        } else if let Some(rest) = line.strip_prefix("TAXIDLIST ") {
            let v = rest.trim();
            if !v.is_empty() {
                alias.raw_taxidlist = Some(v.to_string());
                alias.taxidlist = Some(dir.join(v));
            }
        }
    }

    Ok(alias)
}

/// Check if a database path has an alias file.
pub fn has_alias(base_path: &Path) -> bool {
    base_path.with_extension("nal").exists() || base_path.with_extension("pal").exists()
}

/// Get the alias file path if it exists.
pub fn alias_path(base_path: &Path) -> Option<PathBuf> {
    let nal = base_path.with_extension("nal");
    if nal.exists() {
        return Some(nal);
    }
    let pal = base_path.with_extension("pal");
    if pal.exists() {
        return Some(pal);
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn test_parse_alias() {
        let dir = std::env::temp_dir().join("blast_alias_test");
        std::fs::create_dir_all(&dir).ok();
        let alias_file = dir.join("test.nal");
        let mut f = std::fs::File::create(&alias_file).unwrap();
        writeln!(f, "# Comment").unwrap();
        writeln!(f, "TITLE Multi-volume test").unwrap();
        writeln!(f, "DBLIST vol1 vol2 vol3").unwrap();
        writeln!(f, "NSEQ 10000").unwrap();
        writeln!(f, "LENGTH 5000000").unwrap();
        writeln!(f, "STATS_NSEQ 9000").unwrap();
        writeln!(f, "STATS_TOTLEN 4500000").unwrap();
        writeln!(f, "FIRST_OID 2").unwrap();
        writeln!(f, "LAST_OID 4").unwrap();
        writeln!(f, "OIDLIST masks.msk").unwrap();
        writeln!(f, "MEMB_BIT 256").unwrap();
        writeln!(f, "GILIST taxa.gil").unwrap();
        writeln!(f, "TILIST traces.til").unwrap();
        writeln!(f, "SILIST ids.sil").unwrap();
        writeln!(f, "TAXIDLIST taxids.txt").unwrap();
        drop(f);

        let alias = parse_alias_file(&alias_file).unwrap();
        assert_eq!(alias.title.as_deref(), Some("Multi-volume test"));
        assert_eq!(alias.dblist.len(), 3);
        assert_eq!(alias.raw_dblist, vec!["vol1", "vol2", "vol3"]);
        assert_eq!(alias.nseq, Some(10000));
        assert_eq!(alias.length, Some(5000000));
        assert_eq!(alias.stats_nseq, Some(9000));
        assert_eq!(alias.stats_total_length, Some(4500000));
        assert_eq!(alias.first_oid, Some(2));
        assert_eq!(alias.last_oid, Some(4));
        assert_eq!(alias.raw_oidlist.as_deref(), Some("masks.msk"));
        assert_eq!(alias.oidlist, Some(dir.join("masks.msk")));
        assert_eq!(alias.memb_bit, Some(256));
        assert_eq!(alias.raw_gilist.as_deref(), Some("taxa.gil"));
        assert_eq!(alias.gilist, Some(dir.join("taxa.gil")));
        assert_eq!(alias.raw_tilist.as_deref(), Some("traces.til"));
        assert_eq!(alias.tilist, Some(dir.join("traces.til")));
        assert_eq!(alias.raw_silist.as_deref(), Some("ids.sil"));
        assert_eq!(alias.silist, Some(dir.join("ids.sil")));
        assert_eq!(alias.raw_taxidlist.as_deref(), Some("taxids.txt"));
        assert_eq!(alias.taxidlist, Some(dir.join("taxids.txt")));

        // `unsupported_filter_warnings` should surface all populated filters.
        let warnings = alias.unsupported_filter_warnings();
        assert_eq!(warnings.len(), 5);
        assert!(warnings[0].contains("MEMB_BIT=256"));
        assert!(warnings[1].contains("GILIST taxa.gil"));
        assert!(warnings[2].contains("TILIST traces.til"));
        assert!(warnings[3].contains("SILIST ids.sil"));
        assert!(warnings[4].contains("TAXIDLIST taxids.txt"));

        std::fs::remove_dir_all(&dir).ok();
    }
}
