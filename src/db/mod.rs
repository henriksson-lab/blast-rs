pub mod alias;
pub(crate) mod defline;
mod index;
pub(crate) mod index_writer;
pub mod makedb;

pub use index::{
    megablast_index_exists, megablast_index_volume_paths, BlastDb, DbType, TaxInfo, TaxNameDb,
};
pub use makedb::{make_db, make_nucleotide_db, make_protein_db};

pub(crate) fn current_blastdb_date() -> String {
    let days_since_epoch = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|duration| (duration.as_secs() / 86_400) as i64)
        .unwrap_or(0);
    let (year, month, day) = civil_from_days(days_since_epoch);
    format!("{year:04}-{month:02}-{day:02}")
}

fn civil_from_days(days_since_unix_epoch: i64) -> (i64, u32, u32) {
    let z = days_since_unix_epoch + 719_468;
    let era = if z >= 0 { z } else { z - 146_096 } / 146_097;
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1_460 + doe / 36_524 - doe / 146_096) / 365;
    let mut year = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let day = doy - (153 * mp + 2) / 5 + 1;
    let month = mp + if mp < 10 { 3 } else { -9 };
    if month <= 2 {
        year += 1;
    }
    (year, month as u32, day as u32)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn civil_from_days_matches_known_utc_dates() {
        assert_eq!(civil_from_days(0), (1970, 1, 1));
        assert_eq!(civil_from_days(10_957), (2000, 1, 1));
        assert_eq!(civil_from_days(20_454), (2026, 1, 1));
    }

    #[test]
    fn current_blastdb_date_is_yyyy_mm_dd() {
        let date = current_blastdb_date();
        assert_eq!(date.len(), 10);
        assert_eq!(date.as_bytes()[4], b'-');
        assert_eq!(date.as_bytes()[7], b'-');
        assert!(date
            .bytes()
            .enumerate()
            .all(|(i, b)| matches!(i, 4 | 7) || b.is_ascii_digit()));
    }

    #[test]
    fn makedb_helpers_are_reexported_from_db_module() {
        let _: fn(DbType, &std::path::Path, &std::path::Path, &str) -> std::io::Result<(u32, u64)> =
            make_db;
        let _: fn(&std::path::Path, &std::path::Path, &str) -> std::io::Result<(u32, u64)> =
            make_nucleotide_db;
        let _: fn(&std::path::Path, &std::path::Path, &str) -> std::io::Result<(u32, u64)> =
            make_protein_db;
    }
}
