# blast-rs

Pure-Rust implementation of NCBI BLAST (Basic Local Alignment Search Tool). Produces byte-identical output to the C reference implementation (NCBI BLAST+ 2.17.0) for blastn, with no C dependencies. Supports blastn, blastp, blastx, tblastn, tblastx, and psiblast.

Based on [NCBI BLAST+ 2.17.0](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.17.0/) source distribution (`ncbi-blast-2.17.0+-src.tar.gz`, also on [GitHub](https://github.com/ncbi/blast)). The C core algorithms in `src/algo/blast/core/` were ported function-by-function to Rust.

## Features

- **Pure Rust** -- zero C/C++ FFI calls, no unsafe dependencies
- **Byte-identical output** to NCBI BLAST+ for blastn searches
- **Faster than C** -- 1.1-10x speedup on real genomes (see benchmarks)
- **All major programs**: blastn, blastp, blastx, tblastn, tblastx, psiblast
- **Subcommand CLI**: `blast-cli blastn`, `blast-cli blastp`, etc.
- **Task presets**: `--task blastn-short` for primer search with automatic parameter tuning
- Protein search with neighborhood-word lookup table (matching NCBI's algorithm)
- PSI-BLAST with Henikoff position-based weighting and proper PSSM computation
- Six-frame translation with correct nucleotide coordinate mapping
- Reads standard BLAST database files (v4/v5 `.nin/.nsq/.nhr`)
- Custom tabular columns: `qlen`, `slen`, `saccver`, `qcovs`, `qcovhsp`, and more
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
    --outfmt "6 qseqid qlen saccver slen sstart send bitscore pident"
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

### Read BLAST databases

```rust
use blast_rs::db::BlastDb;
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

### Parse FASTA files (powered by noodles-fasta)

```rust
use blast_rs::input::parse_fasta;
use std::fs::File;

let records = parse_fasta(File::open("query.fa").unwrap());
for rec in &records {
    println!(">{} ({} bases)", rec.id, rec.sequence.len());
}
```

### Search engine (high-level API)

```rust
use blast_rs::blastn::BlastnSearch;

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
use blast_rs::blastn::{BlastnSearch, Strand};

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
| `blast_rs::blastn` | High-level builder API for blastn searches |
| `blast_rs::search` | Core nucleotide search engine (seed scanning, ungapped/gapped extension) |
| `blast_rs::protein_lookup` | Protein word lookup table with neighborhood word generation |
| `blast_rs::protein` | Protein ungapped extension and scoring |
| `blast_rs::pssm` | PSI-BLAST PSSM computation (Henikoff weighting, pseudocounts) |
| `blast_rs::stat` | Karlin-Altschul statistics, KBP computation |
| `blast_rs::traceback` | Gapped alignment with traceback |
| `blast_rs::db` | BLAST database reader/writer (v4/v5 `.nin/.nsq/.nhr`) |
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

- **Neighborhood word lookup table**: For each query w-mer, all amino acid words scoring >= threshold are hashed into a table with PV array for fast rejection, matching `BlastAaLookupIndexQuery`
- **Two-phase extension**: Ungapped X-dropoff extension for filtering, then gapped X-dropoff DP with affine gap penalties for high-scoring seeds
- **Gapped alignment**: Bidirectional Smith-Waterman style DP with substitution matrix (e.g., BLOSUM62), gap open/extend penalties, and X-dropoff banding
- **Six-frame translation**: Forward and reverse frames with correct nucleotide coordinate mapping

### PSI-BLAST

- **Henikoff position-based sequence weighting**: Matches `_PSIComputeSequenceWeights`
- **Effective observations estimation**: Entropy-based calculation matching `s_effectiveObservations`
- **Column-specific pseudocounts**: Matches `s_columnSpecificPseudocounts` with NCBI constants
- **PSSM score conversion**: Frequency ratio blending with BLOSUM62 background, matching `_PSIConvertFreqRatiosToPSSM`

## Compatibility

- Reads BLAST databases created by NCBI `makeblastdb` (v4 and v5 format)
- Output matches NCBI BLAST+ 2.17.0 byte-for-byte for blastn tabular format (`-outfmt 6`)
- Supports scoring parameters: reward/penalty 1/-1 through 5/-4, gap costs 0/0 through 12/8
- Protein programs (blastp/blastx/tblastn/tblastx) use lookup-table seeding with gapped X-dropoff DP extension

## What's not ported

The NCBI BLAST+ 2.17.0 tarball includes a large subset of the NCBI C++ Toolkit (~140MB source). The following are **not** yet included in blast-rs:

- **Taxonomy lookups** -- `staxid`, `ssciname`, `sskingdom` output columns require `.ndb`/`.ntf`/`.nto` taxonomy files; not yet supported
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

## License

MIT OR Unlicense
