#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RUST_BIN="${BLAST_RS_CLI_BIN:-$ROOT_DIR/target/debug/blast-cli}"

need_bin() {
    if [[ ! -x "$1" ]]; then
        echo "missing executable: $1" >&2
        exit 1
    fi
}

need_bin "$RUST_BIN"
need_bin /usr/bin/blastp
need_bin /usr/bin/blastx
need_bin /usr/bin/tblastn
need_bin /usr/bin/tblastx

tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT

printf ">q1 low complexity protein query\nAAAAAAAAAAAAAAAAAAAA\n" > "$tmpdir/protein_query.fa"
printf ">s1 low complexity protein subject\nAAAAAAAAAAAAAAAAAAAA\n>s2 mixed protein subject\nACDEFGHIKLMNPQRSTVWY\n" > "$tmpdir/protein_subject.fa"
printf ">q1 low complexity nt query\nGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCT\n" > "$tmpdir/nt_query.fa"
printf ">s1 low complexity nt subject\nGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCT\n>s2 mixed nt subject\nGCTTGTGATGAATTTGGTCATATAAAACTGATGAATCCTCAACGTTCTACTGTGGTA\n" > "$tmpdir/nt_subject.fa"

run_case() {
    local rust_prog="$1"
    local ncbi_prog="$2"
    local query="$3"
    local subject="$4"
    local label="$5"
    local rust_out="$tmpdir/${label}.rust.tsv"
    local ncbi_out="$tmpdir/${label}.ncbi.tsv"

    "$RUST_BIN" "$rust_prog" \
        --query "$query" \
        --subject "$subject" \
        --outfmt "6 qseqid sseqid length bitscore evalue" \
        --num_threads 1 \
        --out "$rust_out"

    "$ncbi_prog" \
        -query "$query" \
        -subject "$subject" \
        -outfmt "6 qseqid sseqid length bitscore evalue" \
        -num_threads 1 \
        -out "$ncbi_out"

    if ! diff -u "$ncbi_out" "$rust_out"; then
        echo "parity mismatch for $label" >&2
        exit 1
    fi
    echo "ok  $label"
}

run_case blastp /usr/bin/blastp "$tmpdir/protein_query.fa" "$tmpdir/protein_subject.fa" blastp_low_complexity
run_case blastx /usr/bin/blastx "$tmpdir/nt_query.fa" "$tmpdir/protein_subject.fa" blastx_low_complexity
run_case tblastn /usr/bin/tblastn "$tmpdir/protein_query.fa" "$tmpdir/nt_subject.fa" tblastn_low_complexity
run_case tblastx /usr/bin/tblastx "$tmpdir/nt_query.fa" "$tmpdir/nt_subject.fa" tblastx_low_complexity
