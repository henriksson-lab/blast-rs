use blast_rs::db::BlastDb;
use blast_rs::search::decode_packed_ncbi2na_with_ambiguity;
use std::path::PathBuf;

fn decode_raw_base(packed: &[u8], pos: usize) -> u8 {
    let byte = packed[pos >> 2];
    (byte >> (6 - 2 * (pos & 3))) & 3
}

fn code_to_char(code: u8) -> char {
    match code {
        0 => 'A',
        1 => 'C',
        2 => 'G',
        3 => 'T',
        4 => 'R',
        5 => 'Y',
        6 => 'M',
        7 => 'K',
        8 => 'W',
        9 => 'S',
        10 => 'B',
        11 => 'D',
        12 => 'H',
        13 => 'V',
        14 => 'N',
        15 => '-',
        _ => '?',
    }
}

fn render_raw_window(packed: &[u8], start: usize, end: usize) -> String {
    (start..end)
        .map(|pos| code_to_char(decode_raw_base(packed, pos)))
        .collect()
}

fn render_decoded_window(decoded: &[u8], start: usize, end: usize) -> String {
    decoded[start..end]
        .iter()
        .map(|&b| code_to_char(b))
        .collect()
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = std::env::args_os();
    let _prog = args.next();
    let db_path = PathBuf::from(
        args.next()
            .ok_or("usage: dump_db_region <db> <accession> <start_1based> <len>")?,
    );
    let accession = args.next().ok_or("missing accession")?;
    let accession = accession.to_string_lossy();
    let start_1based: usize = args
        .next()
        .ok_or("missing start_1based")?
        .to_string_lossy()
        .parse()?;
    let len: usize = args
        .next()
        .ok_or("missing len")?
        .to_string_lossy()
        .parse()?;

    let db = BlastDb::open(&db_path)?;
    let oid = (0..db.num_oids)
        .find(|&oid| db.get_accession(oid).as_deref() == Some(accession.as_ref()))
        .ok_or("accession not found")?;

    let packed = db.get_sequence(oid);
    let seq_len = db.get_seq_len(oid) as usize;
    let start = start_1based.saturating_sub(1);
    let end = (start + len).min(seq_len);

    println!("oid={oid} seq_len={seq_len} start0={start} end0={end}");
    println!("raw={}", render_raw_window(packed, start, end));

    match db.get_ambiguity_data(oid) {
        Some(amb) => {
            println!("ambiguity_len={}", amb.len());
            let decoded = decode_packed_ncbi2na_with_ambiguity(packed, seq_len, amb);
            println!("decoded={}", render_decoded_window(&decoded, start, end));
        }
        None => {
            println!("ambiguity_len=0");
        }
    }

    Ok(())
}
