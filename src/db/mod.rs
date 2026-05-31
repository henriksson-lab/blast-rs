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
    current_blastdb_date_from_system_time(std::time::SystemTime::now())
}

fn current_blastdb_date_from_system_time(time: std::time::SystemTime) -> String {
    let seconds_since_epoch = time
        .duration_since(std::time::UNIX_EPOCH)
        .map(|duration| duration.as_secs() as i64)
        .unwrap_or(0);
    let (days_since_epoch, seconds_of_day) = local_days_and_seconds(seconds_since_epoch);
    let (year, month, day) = civil_from_days(days_since_epoch);
    let hour24 = (seconds_of_day / 3_600) as u32;
    let minute = ((seconds_of_day % 3_600) / 60) as u32;
    let suffix = if hour24 < 12 { "AM" } else { "PM" };
    let hour12 = match hour24 % 12 {
        0 => 12,
        hour => hour,
    };
    format!(
        "{} {day:02}, {year:04}  {hour12:02}:{minute:02} {suffix}",
        month_abbrev(month)
    )
}

#[cfg(unix)]
fn local_days_and_seconds(seconds_since_epoch: i64) -> (i64, i64) {
    #[repr(C)]
    struct Tm {
        tm_sec: i32,
        tm_min: i32,
        tm_hour: i32,
        tm_mday: i32,
        tm_mon: i32,
        tm_year: i32,
        tm_wday: i32,
        tm_yday: i32,
        tm_isdst: i32,
        tm_gmtoff: isize,
        tm_zone: *const std::ffi::c_char,
    }

    unsafe extern "C" {
        fn localtime_r(timep: *const i64, result: *mut Tm) -> *mut Tm;
    }

    let mut tm = Tm {
        tm_sec: 0,
        tm_min: 0,
        tm_hour: 0,
        tm_mday: 1,
        tm_mon: 0,
        tm_year: 70,
        tm_wday: 0,
        tm_yday: 0,
        tm_isdst: 0,
        tm_gmtoff: 0,
        tm_zone: std::ptr::null(),
    };
    let mut seconds = seconds_since_epoch;
    let ok = unsafe { !localtime_r(&raw mut seconds, &raw mut tm).is_null() };
    if ok {
        let year = tm.tm_year as i64 + 1900;
        let month = (tm.tm_mon + 1) as u32;
        let day = tm.tm_mday as u32;
        let days = days_from_civil(year, month, day);
        let seconds_of_day =
            (tm.tm_hour as i64 * 3_600) + (tm.tm_min as i64 * 60) + tm.tm_sec as i64;
        (days, seconds_of_day)
    } else {
        utc_days_and_seconds(seconds_since_epoch)
    }
}

#[cfg(not(unix))]
fn local_days_and_seconds(seconds_since_epoch: i64) -> (i64, i64) {
    utc_days_and_seconds(seconds_since_epoch)
}

fn utc_days_and_seconds(seconds_since_epoch: i64) -> (i64, i64) {
    (
        seconds_since_epoch.div_euclid(86_400),
        seconds_since_epoch.rem_euclid(86_400),
    )
}

fn month_abbrev(month: u32) -> &'static str {
    match month {
        1 => "Jan",
        2 => "Feb",
        3 => "Mar",
        4 => "Apr",
        5 => "May",
        6 => "Jun",
        7 => "Jul",
        8 => "Aug",
        9 => "Sep",
        10 => "Oct",
        11 => "Nov",
        12 => "Dec",
        _ => "Jan",
    }
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

fn days_from_civil(year: i64, month: u32, day: u32) -> i64 {
    let year = year - i64::from(month <= 2);
    let era = if year >= 0 { year } else { year - 399 } / 400;
    let yoe = year - era * 400;
    let month = month as i64;
    let doy = (153 * (month + if month > 2 { -3 } else { 9 }) + 2) / 5 + day as i64 - 1;
    let doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
    era * 146_097 + doe - 719_468
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
    fn days_from_civil_round_trips_known_dates() {
        for date in [(1970, 1, 1), (2000, 1, 1), (2026, 1, 1), (2026, 5, 31)] {
            assert_eq!(
                civil_from_days(days_from_civil(date.0, date.1, date.2)),
                date
            );
        }
    }

    #[test]
    fn current_blastdb_date_is_ncbi_makeblastdb_style() {
        let date = current_blastdb_date();
        assert_eq!(date.len(), 22);
        assert!(matches!(
            &date[..3],
            "Jan"
                | "Feb"
                | "Mar"
                | "Apr"
                | "May"
                | "Jun"
                | "Jul"
                | "Aug"
                | "Sep"
                | "Oct"
                | "Nov"
                | "Dec"
        ));
        assert_eq!(&date[6..8], ", ");
        assert_eq!(&date[12..14], "  ");
        assert_eq!(&date[16..17], ":");
        assert_eq!(&date[19..20], " ");
        assert!(matches!(&date[20..], "AM" | "PM"));
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
