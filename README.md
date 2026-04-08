# blast-rs

Pure-Rust implementation of NCBI BLAST (Basic Local Alignment Search Tool) for nucleotide sequence similarity searching. Produces byte-identical output to the C reference implementation (NCBI BLAST+ 2.17.0) with no C dependencies.

## Features

- **Pure Rust** -- zero C/C++ FFI calls, no unsafe dependencies
- **Byte-identical output** to NCBI BLAST+ for blastn searches
- **Faster than C** -- 1.1-10x speedup on real genomes (see benchmarks)
- Reads standard BLAST database files (v4/v5 `.nin/.nsq/.nhr`)
- Supports all standard blastn options (word size, scoring, e-value, output formats)
- FASTA-vs-FASTA search (`--subject` mode) without pre-built database
- Multi-threaded search via rayon

## Installation

```bash
cargo install blast-cli
```

Or build from source:

```bash
git clone https://github.com/henriksson-lab/blast-rs
cd blast-rs
cargo build --release
# Binary at target/release/blast-cli
```

For maximum performance, compile for your native CPU:

```bash
RUSTFLAGS="-C target-cpu=native" cargo build --release
```

## CLI Usage

### Basic search against a BLAST database

```bash
blast-cli --query query.fa --db mydb --word_size 11 --outfmt 6
```

### Search against a FASTA file (no database needed)

```bash
blast-cli --query query.fa --subject reference.fa
```

### Common options

```bash
# Custom scoring
blast-cli --query query.fa --db mydb \
    --reward 2 --penalty -3 \
    --gapopen 5 --gapextend 2 \
    --evalue 0.001

# Plus strand only, no DUST masking
blast-cli --query query.fa --db mydb --strand plus --dust no

# Limit results
blast-cli --query query.fa --db mydb --max_target_seqs 100

# Pairwise text output
blast-cli --query query.fa --db mydb --outfmt 0

# Custom tabular columns
blast-cli --query query.fa --db mydb --outfmt "6 qseqid sseqid pident length evalue"
```

### All options

| Option | Default | Description |
|--------|---------|-------------|
| `--query` | required | Query FASTA file |
| `--db` | | BLAST database path (without extension) |
| `--subject` | | Subject FASTA file (alternative to `--db`) |
| `--out` | stdout | Output file |
| `--outfmt` | `6` | Output format: 0=pairwise, 5=XML, 6=tabular, 17=SAM |
| `--evalue` | `10.0` | E-value threshold |
| `--word_size` | `11` | Word size for initial seed |
| `--reward` | `1` | Match reward |
| `--penalty` | `-3` | Mismatch penalty |
| `--gapopen` | `5` | Gap opening cost |
| `--gapextend` | `2` | Gap extension cost |
| `--strand` | `both` | Query strand: `both`, `plus`, or `minus` |
| `--dust` | `yes` | DUST low-complexity filtering |
| `--max_target_seqs` | `500` | Maximum number of target sequences |
| `--num_threads` | `1` | Number of threads |
| `--xdrop_gap` | `30` | X-dropoff for gapped extension (bits) |
| `--xdrop_gap_final` | `100` | X-dropoff for final gapped extension (bits) |
| `--perc_identity` | `0` | Minimum percent identity |

## Library Usage

The project is split into reusable crates:

### blast-db -- Read BLAST databases

```rust
use blast_db::BlastDb;
use std::path::Path;

let db = BlastDb::open(Path::new("mydb")).unwrap();
println!("{} sequences, {} total bases", db.num_oids, db.total_length);

for oid in 0..db.num_oids {
    let packed_seq = db.get_sequence(oid);
    let seq_len = db.get_seq_len(oid);
    let accession = db.get_accession(oid);
    println!("{}: {} bases", accession.unwrap_or_default(), seq_len);
}
```

### blast-input -- Parse FASTA files

```rust
use blast_input::parse_fasta;
use std::fs::File;

let records = parse_fasta(File::open("query.fa").unwrap());
for rec in &records {
    println!(">{} ({} bases)", rec.id, rec.sequence.len());
}
```

### blast-core -- Search engine (high-level API)

```rust
use blast_core::blastn::BlastnSearch;

let results = BlastnSearch::new()
    .query(b"ACGTACGTACGTACGTACGTACGTACGT")
    .subject(b"NNNNNACGTACGTACGTACGTACGTACGTNNNNN")
    .run();

for hit in &results {
    println!("q={}-{} s={}-{} score={} evalue={:.2e}",
        hit.query_start, hit.query_end,
        hit.subject_start, hit.subject_end,
        hit.score, hit.evalue);
}
```

With custom parameters:

```rust
use blast_core::blastn::{BlastnSearch, Strand};

let results = BlastnSearch::new()
    .query(b"ACGTACGTACGT")
    .subject(b"ACGTACGTACGT")
    .word_size(7)
    .reward(2)
    .penalty(-3)
    .gap_open(5)
    .gap_extend(2)
    .evalue(0.001)
    .dust(false)
    .strand(Strand::Plus)
    .run();
```

### blast-format -- Output formatting

```rust
use blast_format::{format_tabular, TabularHit};
use std::io;

let hits: Vec<TabularHit> = /* ... */;
format_tabular(&mut io::stdout(), &hits).unwrap();
```

## Crate Structure

| Crate | Description |
|-------|-------------|
| `blast-cli` | Command-line binary |
| `blast-core` | Search engine, statistics, gapped alignment, traceback |
| `blast-db` | BLAST database reader (v4/v5 `.nin/.nsq/.nhr`) and writer |
| `blast-input` | FASTA parser with BLASTNA/NCBIstdaa encoding |
| `blast-format` | Output formatting (tabular, pairwise, XML, SAM) |

## Benchmarks

500bp query, word_size=11, compiled with `RUSTFLAGS="-C target-cpu=native"`.
Produces byte-identical output to the C reference.

### Single-threaded

| Database | Size | NCBI BLAST+ 2.17.0 | blast-rs | Speedup |
|----------|------|--------------------:|----------:|--------:|
| S. cerevisiae | 12 MB | 416 ms | **75 ms** | **5.5x** |
| C. elegans | 100 MB | 717 ms | **630 ms** | **1.1x** |

### Multi-threaded (4 threads)

| Database | Size | NCBI BLAST+ 2.17.0 | blast-rs | Speedup |
|----------|------|--------------------:|----------:|--------:|
| S. cerevisiae | 12 MB | 405 ms | **39 ms** | **10.4x** |
| C. elegans | 100 MB | 594 ms | **263 ms** | **2.3x** |

Both engines use BLAST database mode (`-db`) for fair comparison.
The speedup comes from eliminating C++ startup overhead (~400ms) and using efficient
Rust implementations of the same BLAST algorithms (packed byte scanning, PV array,
diagonal tracking, X-dropoff DP). Multi-threaded search parallelizes across database
sequences using rayon.

## Algorithms

This is a faithful port of the NCBI BLAST+ C engine algorithms:

- **Seed scanning**: Step-4 packed byte scanning with 8-mer lookup table and presence vector (PV) array, matching `s_BlastSmallNaScanSubject_8_4`
- **Word extension**: Left/right extension from 8-mer to full word size, matching `s_BlastSmallNaExtend`
- **Diagonal tracking**: Hash-based diagonal table to prevent redundant extensions, matching `BLAST_DiagTable`
- **Ungapped extension**: X-dropoff extension with scoring matrix, matching `s_NuclUngappedExtend`
- **Gapped alignment**: Bidirectional X-dropoff DP with traceback, matching `ALIGN_EX`
- **Score-only preliminary extension**: Fast DP without traceback for pre-filtering, matching `Blast_SemiGappedAlign`
- **Statistics**: Full Karlin-Altschul parameter computation (Newton-Raphson Lambda, DP K, Horner H), per-context KBP, length adjustment
- **Scoring matrix**: Full BLASTNA 16x16 matrix with ambiguity codes
- **DUST masking**: Low-complexity filtering with separate masked/unmasked query handling

## Compatibility

- Reads BLAST databases created by NCBI `makeblastdb` (v4 and v5 format)
- Output matches NCBI BLAST+ 2.17.0 byte-for-byte for blastn tabular format (`-outfmt 6`)
- Supports scoring parameters: reward/penalty 1/-1 through 5/-4, gap costs 0/0 through 12/8

## License

MIT OR Unlicense
