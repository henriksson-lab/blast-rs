use std::hint::black_box;
use std::time::{Duration, Instant};

use blast_rs::encoding::encode_blastna_sequence;
use blast_rs::traceback::{blast_gapped_alignment_with_traceback, blast_semi_gapped_align};

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

    let query = encode_blastna_sequence(query.as_bytes());
    let subject = encode_blastna_sequence(subject.as_bytes());
    let seed_q = query.len() / 2;
    let seed_s = subject.len() / 2;

    timed("score_only", iterations, || {
        blast_semi_gapped_align(
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
        blast_gapped_alignment_with_traceback(
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
        black_box(blast_semi_gapped_align(
            &query, &subject, seed_q, seed_s, 1, -3, 5, 2, 30,
        ));
    }
}
