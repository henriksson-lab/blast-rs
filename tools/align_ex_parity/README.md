# align_ex parity harness

A minimal C translation of NCBI's `ALIGN_EX` (blast_gapalign.c:374), stripped
to the DP core (no `BlastGapAlignStruct`, no `GapPrelimEditBlock`, no fence
sentry). Run both Rust's `align_ex` and this C harness with tracehash-rs
instrumentation to prove (or disprove) DP-level parity row-by-row.

## Result (2026-04 run)

On the `adjacent_del_ins` fixture (`blastn-short`, reward=1, penalty=-3,
gap_open=5, gap_extend=2, x_dropoff=16, seed=(17,17)):

    left: M=18, N=18, score=18
    right: M=64, N=59, score=44, a_off=64, b_off=59
    total: 62

All **82 `align_ex_row` events** (18 left + 64 right) and all **82
`align_ex_tb` traceback-step events** produce **byte-identical FNV-64 hashes**
on both sides. See `/tmp/rust-trace.tsv` vs `/tmp/ncbi-trace.tsv`.

Conclusion: Rust's `src/traceback.rs::align_ex` is a character-exact port of
NCBI's `ALIGN_EX`. Any divergence from NCBI's actual blastn output on this
fixture (Rust produces BTOP `34A-G-G-G-T-43`, NCBI `35G-G-G-T-A-42`) lives
**outside `ALIGN_EX`** — likely in seed selection or downstream HSP
post-processing (`Blast_HSPReevaluateWithAmbiguitiesGapped`, endpoint
extension in `blast_hits.c:614-637`).

## Build

    TH_DIR=$(find ~/.cargo/registry/src -type d -name 'tracehash-rs-0.1.0' | head -1)/c
    cc -O0 -g -Wall -Wextra -I "$TH_DIR" \
       ncbi_align_ex.c "$TH_DIR/tracehash_c.c" -o ncbi_align_ex

## Run

    TRACEHASH_OUT=/tmp/ncbi-trace.tsv TRACEHASH_SIDE=ncbi \
      TRACEHASH_RUN_ID=adj_del_ins ./ncbi_align_ex

    TRACEHASH_OUT=/tmp/rust-trace.tsv TRACEHASH_SIDE=rust \
      TRACEHASH_RUN_ID=adj_del_ins \
      cargo test --release --lib trace_align_ex_adjacent_del_ins \
        -- --ignored --nocapture

    tracehash-compare /tmp/rust-trace.tsv /tmp/ncbi-trace.tsv

## Next step if the divergence returns

1. Confirm both sides still emit identical hashes. If so, divergence is
   downstream.
2. Instrument `blast_traceback.c:439` (`BlastGetOffsetsForGappedAlignment`)
   in NCBI and emit the computed seed as a `gapped_start` event.
3. Instrument `Blast_HSPReevaluateWithAmbiguitiesGapped` (`blast_hits.c:479`)
   and emit the edit-script state pre/post.
