# blast-rs

Pure-Rust implementation of NCBI BLAST (Basic Local Alignment Search Tool). Produces byte-identical output to the C reference implementation (NCBI BLAST+ 2.17.0) for blastn, with no C dependencies. Supports blastn, blastp, blastx, tblastn, tblastx, and psiblast.

Based on [NCBI BLAST+ 2.17.0](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.17.0/) source distribution (`ncbi-blast-2.17.0+-src.tar.gz`, also on [GitHub](https://github.com/ncbi/blast)). The C core algorithms in `src/algo/blast/core/` were ported function-by-function to Rust.

This is a translation of the original code and not the authoritative implementation. This code should generate bitwise
equal output to the original. Please report any deviations.

The aim of this project is to increase performance, especially by providing this code through a type-safe library interface.
The code can also be compiled to be used for webassembly.


## Features

- **Pure Rust** -- zero C/C++ FFI calls, no unsafe dependencies
- **Byte-identical output** to NCBI BLAST+ for blastn searches
- **Faster than C** -- 1.1-10x speedup on real genomes (see benchmarks)
- **All major programs**: blastn, blastp, blastx, tblastn, tblastx, psiblast
- **High-level library API**: `blastp()`, `blastn_search()`, `blastx()`, `tblastn()`, `tblastx()` with builder-pattern `SearchParams`
- **Database creation**: `BlastDbBuilder` for creating protein and nucleotide databases programmatically
- **Subcommand CLI**: `blast-cli blastn`, `blast-cli blastp`, etc.
- **Task presets**: `--task blastn-short` for primer search with automatic parameter tuning
- Protein search with BLOSUM62 scoring and neighborhood-word lookup table (matching NCBI's algorithm)
- PSI-BLAST with Henikoff position-based weighting and proper PSSM computation
- Six-frame translation with correct nucleotide coordinate mapping
- Reads standard BLAST database files (v4/v5 `.nin/.nsq/.nhr`) with taxonomy support (`.nto`, `taxdb.bti`/`taxdb.btd`)
- **Full tabular field support**: `qseq`, `sseq`, `qframe`, `sframe`, `score`, `staxid`, `ssciname`, `scomname`, `sskingdom`, `sblastname`, and all standard columns
- FASTA-vs-FASTA search (`--subject` mode) without pre-built database
- Multi-threaded search via rayon
- **185 tests**: 140 unit tests + 43 integration tests + doc tests

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
blast-cli blastn --query query.fa --db mydb --word_size 11 --outfmt 6
```

### Search against a FASTA file (no database needed)

```bash
blast-cli blastn --query query.fa --subject reference.fa
```

### Short primer search

```bash
blast-cli blastn --query primers.fa --db mydb --task blastn-short \
    --max_hsps 1 --max_target_seqs 10000
```

### Common options

```bash
# Custom scoring
blast-cli blastn --query query.fa --db mydb \
    --reward 2 --penalty -3 \
    --gapopen 5 --gapextend 2 \
    --evalue 0.001

# Plus strand only, no DUST masking
blast-cli blastn --query query.fa --db mydb --strand plus --dust no

# Limit results
blast-cli blastn --query query.fa --db mydb --max_target_seqs 100

# Pairwise text output
blast-cli blastn --query query.fa --db mydb --outfmt 0

# Custom tabular columns
blast-cli blastn --query query.fa --db mydb \
    --outfmt "6 qseqid sseqid pident qseq sseq staxid ssciname"
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

### Protein search

```bash
# blastp: protein query vs protein subject
blast-cli blastp --query protein.fa --subject database.fa

# blastx: translated nucleotide query vs protein subject
blast-cli blastx --query nucleotide.fa --subject protein.fa

# tblastn: protein query vs translated nucleotide subject
blast-cli tblastn --query protein.fa --subject nucleotide.fa

# psiblast: iterative position-specific search
blast-cli psiblast --query protein.fa --subject database.fa
```

## Library Usage

### High-level search API

```rust
use blast_rs::api::{SearchParams, blastp, blastn_search, BlastDbBuilder, SequenceEntry};
use blast_rs::db::DbType;

// Build a protein database
let mut builder = BlastDbBuilder::new(DbType::Protein, "my db");
builder.add(SequenceEntry {
    title: "target protein".into(),
    accession: "P001".into(),
    sequence: b"MKFLILLFNILCLFPVLAADNHGVSMNAS".to_vec(),
    taxid: None,
});
builder.write(std::path::Path::new("mydb")).unwrap();

// Search
let db = blast_rs::db::BlastDb::open(std::path::Path::new("mydb")).unwrap();
let params = SearchParams::blastp()
    .evalue(10.0)
    .num_threads(1)
    .filter_low_complexity(false);

let results = blastp(&db, b"MKFLILLFNILCLFPVLAADNHGVSMNAS", &params);
for r in &results {
    for hsp in &r.hsps {
        println!("score={} evalue={:.2e} identity={:.1}%",
            hsp.score, hsp.evalue, hsp.percent_identity());
    }
}
```

### Blastn builder API

```rust
use blast_rs::BlastnSearch;

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

### Utility functions

```rust
use blast_rs::api::{reverse_complement, six_frame_translate, parse_fasta};

// Reverse complement
let rc = reverse_complement(b"ATGGCTAGCGATCG");
assert_eq!(&rc, b"CGATCGCTAGCCAT");

// Six-frame translation
let frames = six_frame_translate(b"ATGGCTAGCGATCGATCGATCGATCG");
assert_eq!(frames[0].frame, 1);
assert_eq!(frames[0].protein[0], b'M'); // ATG = Met

// FASTA parsing
let seqs = parse_fasta(b">seq1\nACGT\n>seq2\nTGCA\n");
assert_eq!(seqs.len(), 2);
```

### Read BLAST databases

```rust
use blast_rs::db::BlastDb;
use std::path::Path;

let db = BlastDb::open(Path::new("mydb")).unwrap();
println!("{} sequences", db.num_oids);

for oid in 0..db.num_oids {
    let packed_seq = db.get_sequence(oid);
    let seq_len = db.get_seq_len(oid);
    let accession = db.get_accession(oid);
    println!("{}: {} bases", accession.unwrap_or_default(), seq_len);
}
```

### Output formatting

```rust
use blast_rs::format::{format_tabular, TabularHit};
use std::io;

let hits: Vec<TabularHit> = /* ... */;
format_tabular(&mut io::stdout(), &hits).unwrap();
```

## Module Structure

Single crate `blast-rs` with modules:

| Module | Description |
|--------|-------------|
| `blast_rs::api` | High-level search API (`blastp`, `blastn`, `blastx`, `tblastn`, `tblastx`, `BlastnSearch` builder, DB builder) |
| `blast_rs::search` | Core nucleotide search engine (seed scanning, ungapped/gapped extension) |
| `blast_rs::protein_lookup` | Protein word lookup table with neighborhood word generation |
| `blast_rs::protein` | Protein ungapped extension, gapped alignment, and scoring |
| `blast_rs::pssm` | PSI-BLAST PSSM computation (Henikoff weighting, pseudocounts) |
| `blast_rs::stat` | Karlin-Altschul statistics, KBP computation |
| `blast_rs::traceback` | Gapped alignment with traceback |
| `blast_rs::matrix` | Scoring matrices (BLOSUM62, nucleotide match/mismatch) |
| `blast_rs::db` | BLAST database reader/writer (v4/v5), taxonomy ID/name lookups (`.nto`, `taxdb`) |
| `blast_rs::input` | FASTA parser (via noodles-fasta) with BLASTNA/NCBIstdaa encoding |
| `blast_rs::format` | Output formatting (tabular, pairwise, XML, SAM) |
| `blast_rs::filter` | DUST low-complexity masking |
| `blast_rs::util` | Six-frame translation, genetic code, coordinate mapping |

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

### Nucleotide (blastn)

- **Seed scanning**: Step-4 packed byte scanning with 8-mer lookup table and presence vector (PV) array, matching `s_BlastSmallNaScanSubject_8_4`
- **Word extension**: Left/right extension from 8-mer to full word size, matching `s_BlastSmallNaExtend`
- **Diagonal tracking**: Hash-based diagonal table to prevent redundant extensions, matching `BLAST_DiagTable`
- **Ungapped extension**: X-dropoff extension with scoring matrix, matching `s_NuclUngappedExtend`
- **Gapped alignment**: Bidirectional X-dropoff DP with traceback, matching `ALIGN_EX`
- **Score-only preliminary extension**: Fast DP without traceback for pre-filtering, matching `Blast_SemiGappedAlign`
- **Statistics**: Full Karlin-Altschul parameter computation (Newton-Raphson Lambda, DP K, Horner H), per-context KBP, length adjustment
- **Scoring matrix**: Full BLASTNA 16x16 matrix with ambiguity codes
- **DUST masking**: Low-complexity filtering with separate masked/unmasked query handling

### Protein (blastp, blastx, tblastn, tblastx)

- **BLOSUM62 scoring matrix**: Full 28x28 NCBIstdaa matrix from NCBI source
- **Neighborhood word lookup table**: For each query w-mer, all amino acid words scoring >= threshold are hashed into a table with PV array for fast rejection, matching `BlastAaLookupIndexQuery`
- **Two-phase extension**: Ungapped X-dropoff extension for filtering, then gapped X-dropoff DP with affine gap penalties for high-scoring seeds
- **Gapped alignment**: Bidirectional Smith-Waterman style DP with substitution matrix, gap open/extend penalties, and X-dropoff banding
- **Six-frame translation**: Forward and reverse frames with correct nucleotide coordinate mapping

### PSI-BLAST

- **Henikoff position-based sequence weighting**: Matches `_PSIComputeSequenceWeights`
- **Effective observations estimation**: Entropy-based calculation matching `s_effectiveObservations`
- **Column-specific pseudocounts**: Matches `s_columnSpecificPseudocounts` with NCBI constants
- **PSSM score conversion**: Frequency ratio blending with BLOSUM62 background, matching `_PSIConvertFreqRatiosToPSSM`

## Compatibility

- Reads BLAST databases created by NCBI `makeblastdb` (v4 and v5 format), including taxonomy (`.nto` taxid files, `taxdb.bti`/`taxdb.btd` name database)
- Output matches NCBI BLAST+ 2.17.0 byte-for-byte for blastn tabular format (`-outfmt 6`)
- Supports scoring parameters: reward/penalty 1/-1 through 5/-4, gap costs 0/0 through 12/8
- Protein programs (blastp/blastx/tblastn/tblastx) use BLOSUM62 with lookup-table seeding and gapped X-dropoff DP extension

## What's not ported

The NCBI BLAST+ 2.17.0 tarball includes a large subset of the NCBI C++ Toolkit (~140MB source). The following are **not** yet included in blast-rs:

- **Composition-based statistics** -- `composition_adjustment/` used by blastp/blastx for compositional score correction
- **IgBLAST** -- immunoglobulin/T-cell receptor analysis (`igblast/`)
- **Remote BLAST** -- network search against NCBI servers (`connect/`)
- **SRA/VDB integration** -- searching against Sequence Read Archive data (`vdb/`, `sra/`)
- **ASN.1 object model** -- the full `objects/` and `serial/` framework (~40MB) for Seq-entry, Bioseq, etc.
- **Object Manager** -- `objmgr/` for GenBank sequence retrieval and data loading
- **Database indexing** -- `dbindex/` for MegaBLAST indexed searches
- **Gumbel parameter estimation** -- `gumbel_params/` for domain-specific statistics
- **Protein k-mer search** -- `proteinkmer/` for fast protein pre-filtering
- **NCBI application framework** -- `corelib/` (logging, config, threading, diagnostics)
- **blast_formatter** -- standalone tool for reformatting archived BLAST results

## Changelog

### 0.7.0

- **High-level library API**: New `blast_rs::api` module with `blastp()`, `blastn_search()`, `blastx()`, `tblastn()`, `tblastx()` functions, `SearchParams` builder, `SearchResult`/`Hsp` types
- **Database creation**: `BlastDbBuilder` for programmatically creating protein and nucleotide BLAST databases
- **BLOSUM62 matrix**: Populated with correct values from NCBI source (was all zeros)
- **43 integration tests** ported from previous implementation covering all search programs, database operations, edge cases, and multi-subject indexing
- **17 bug fixes**: gapalign traceback, protein scoring matrix, i16 overflow, RPS search space, tabular output fields, log1p math, length adjustment formula, encoding tables, and more (see TODO.md)
- Utility functions: `reverse_complement()`, `six_frame_translate()`, `parse_fasta()` for ASCII sequences
- Lowercase support in encoding tables (IUPACNA, NCBI4NA, amino acid)
- OID bounds checking in database accessors

### 0.6.0

- Initial pure-Rust implementation with CLI and core search engine

## License

MIT OR Unlicense
