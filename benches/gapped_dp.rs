use std::hint::black_box;
use std::time::{Duration, Instant};

use blast_rs::traceback::{blast_gapped_align, blast_gapped_score_only};

fn encode_blastna(seq: &str) -> Vec<u8> {
    seq.as_bytes()
        .iter()
        .map(|&b| match b {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' | b'U' | b'u' => 3,
            b'R' | b'r' => 4,
            b'Y' | b'y' => 5,
            b'M' | b'm' => 6,
            b'K' | b'k' => 7,
            b'W' | b'w' => 8,
            b'S' | b's' => 9,
            b'B' | b'b' => 10,
            b'D' | b'd' => 11,
            b'H' | b'h' => 12,
            b'V' | b'v' => 13,
            b'N' | b'n' => 14,
            _ => 15,
        })
        .collect()
}

fn timed(name: &str, iterations: usize, mut f: impl FnMut() -> i32) {
    let start = Instant::now();
    let mut checksum = 0i32;
    for _ in 0..iterations {
        checksum ^= black_box(f());
    }
    let elapsed = start.elapsed();
    let per_iter = elapsed / iterations as u32;
    println!(
        "{name}: {iterations} iterations in {}.{:03}s, {:?}/iter, checksum {}",
        elapsed.as_secs(),
        elapsed.subsec_millis(),
        per_iter,
        checksum
    );
}

fn main() {
    let iterations = std::env::var("BLAST_RS_GAPPED_DP_BENCH_ITERS")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or(20_000);

    let mut query = String::new();
    let mut subject = String::new();
    for i in 0..32 {
        query.push_str("ACGTACGTACGTACGT");
        if i % 7 == 0 {
            subject.push_str("ACGTACGTNNNNACGT");
        } else if i % 5 == 0 {
            subject.push_str("ACGTACGTTTACGT");
        } else {
            subject.push_str("ACGTACGTACGTACGT");
        }
    }

    let query = encode_blastna(&query);
    let subject = encode_blastna(&subject);
    let seed_q = query.len() / 2;
    let seed_s = subject.len() / 2;

    timed("score_only", iterations, || {
        blast_gapped_score_only(
            black_box(&query),
            black_box(&subject),
            seed_q,
            seed_s,
            1,
            -3,
            5,
            2,
            30,
        )
    });

    let traceback_iterations = (iterations / 10).max(1);
    timed("traceback", traceback_iterations, || {
        blast_gapped_align(
            black_box(&query),
            black_box(&subject),
            seed_q,
            seed_s,
            1,
            -3,
            5,
            2,
            30,
        )
        .map(|tb| tb.score)
        .unwrap_or_default()
    });

    let smoke_start = Instant::now();
    while smoke_start.elapsed() < Duration::from_millis(1) {
        black_box(blast_gapped_score_only(
            &query, &subject, seed_q, seed_s, 1, -3, 5, 2, 30,
        ));
    }
}
