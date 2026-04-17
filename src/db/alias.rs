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

        std::fs::remove_dir_all(&dir).ok();
    }
}
