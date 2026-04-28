# blast-rs

Pure-Rust implementation of NCBI BLAST (Basic Local Alignment Search Tool). Produces byte-identical output to the C reference implementation (NCBI BLAST+ 2.17.0) for blastn, with no C dependencies. Supports blastn, blastp, blastx, tblastn, tblastx, and psiblast.

Based on [NCBI BLAST+ 2.17.0](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.17.0/) source distribution (`ncbi-blast-2.17.0+-src.tar.gz`, also on [GitHub](https://github.com/ncbi/blast)). The C core algorithms in `src/algo/blast/core/` were ported function-by-function to Rust.

* 2026-04-28: Errors in this library has been found. do NOT use it! (some things may work though)

## This is an LLM-mediated faithful (hopefully) translation, not the original code!

Most users should probably first see if the existing original code works for them, unless they have reason otherwise. The original source
may have newer features and it has had more love in terms of fixing bugs. In fact, we aim to replicate bugs if they are present, for the
sake of reproducibility! (but then we might have added a few more in the process)

There are however cases when you might prefer this Rust version. We generally agree with [this page](https://rewrites.bio/)
but more specifically:
* We have had many issues with ensuring that our software works using existing containers (Docker, PodMan, Singularity). One size does not fit all and it eats our resources trying to keep up with every way of delivering software
* Common package managers do not work well. It was great when we had a few Linux distributions with stable procedures, but now there are just too many ecosystems (Homebrew, Conda). Conda has an NP-complete resolver which does not scale. Homebrew is only so-stable. And our dependencies in Python still break. These can no longer be considered professional serious options. Meanwhile, Cargo enables multiple versions of packages to be available, even within the same program(!)
* The future is the web. We deploy software in the web browser, and until now that has meant Javascript. This is a language where even the == operator is broken. Typescript is one step up, but a game changer is the ability to compile Rust code into webassembly, enabling performance and sharing of code with the backend. Translating code to Rust enables new ways of deployment and running code in the browser has especial benefits for science - researchers do not have deep pockets to run servers, so pushing compute to the user enables deployment that otherwise would be impossible
* Old CLI-based utilities are bad for the environment(!). A large amount of compute resources are spent creating and communicating via small files, which we can bypass by using code as libraries. Even better, we can avoid frequent reloading of databases by hoisting this stage, with up to 100x speedups in some cases. Less compute means faster compute and less electricity wasted
* LLM-mediated translations may actually be safer to use than the original code. This article shows that [running the same code on different operating systems can give somewhat different answers](https://doi.org/10.1038/nbt.3820). This is a gap that Rust+Cargo can reduce. Typesafe interfaces also reduce coding mistakes and error handling, as opposed to typical command-line scripting

But:

* **This approach should still be considered experimental**. The LLM technology is immature and has sharp corners. But there are opportunities to reap, and the genie is not going back to the bottle. This translation is as much aimed to learn how to improve the technology and get feedback on the results.
* Translations are not endorsed by the original authors unless otherwise noted. **Do not send bug reports to the original developers**. Use our Github issues page instead.
* **Do not trust the benchmarks on this page**. They are used to help evaluate the translation. If you want improved performance, you generally have to use this code as a library, and use the additional tricks it offers. We generally accept performance losses in order to reduce our dependency issues
* **Check the original Github pages for information about the package**. This README is kept sparse on purpose. It is not meant to be the primary source of information



## Features

- **Pure Rust** -- zero C/C++ FFI calls, no unsafe dependencies
- **Byte-identical output** to NCBI BLAST+ for blastn searches
- **Low startup overhead** for CLI use, with a library API for avoiding repeated process and database setup
- **All major programs**: blastn, blastp, blastx, tblastn, tblastx, psiblast
- **High-level library API**: `blastp()`, `blastn_search()`, `blastx()`, `tblastn()`, `tblastx()` with builder-pattern `SearchParams`
- **Database creation**: `BlastDbBuilder` for creating protein and nucleotide databases programmatically
- **Subcommand CLI**: `blast-cli blastn`, `blast-cli blastp`, etc.
- **Task presets**: `--task blastn-short` for primer search with automatic parameter tuning
- Protein search with BLOSUM62 scoring and neighborhood-word lookup table (matching NCBI's algorithm)
- PSI-BLAST with Henikoff position-based weighting and proper PSSM computation
- Six-frame translation with correct nucleotide coordinate mapping
- Reads standard BLAST database files (v4/v5 `.nin/.nsq/.nhr`) with taxonomy support (`.nto`, `taxdb.bti`/`taxdb.btd`)
- **Full tabular field support**: `qseq`, `sseq`, `qframe`, `sframe`, `sstrand`, `score`, `staxid`, `ssciname`, `scomname`, `sskingdom`, `sblastname`, BTOP, commented tabular, CSV, and all standard columns
- FASTA-vs-FASTA search (`--subject` mode) without pre-built database
- Multi-threaded search via rayon
- **711 passing release tests**: 406 library unit tests + 13 CLI unit tests + 287 integration tests + 4 stress tests + 1 harness test, plus ignored parity tests

## Installation

```bash
cargo install blast-rs
```

Or build from source:

```bash
git clone https://github.com/henriksson-lab/blast-rs
cd blast-rs
cargo build --release --bin blast-cli
# Binary at target/release/blast-cli
```

You can also run the CLI directly from a checkout without installing it:

```bash
cargo run --release --bin blast-cli -- blastn --query query.fa --subject reference.fa
```

For maximum performance, compile for your native CPU:

```bash
RUSTFLAGS="-C target-cpu=native" cargo build --release --bin blast-cli
```

## CLI Usage

The installed command is `blast-cli`. Its first argument is the BLAST program
to run, such as `blastn`, `blastp`, `blastx`, `tblastn`, `tblastx`, or
`psiblast`.

Use `--db` for an existing BLAST database prefix. The prefix is the path without
the database extension; for example, pass `--db nt_subset` for files such as
`nt_subset.nhr`, `nt_subset.nin`, and `nt_subset.nsq`. Use `--subject` for a
plain FASTA file when you do not want to build a database first.

### Basic search against a BLAST database

```bash
blast-cli blastn --query query.fa --db mydb --word_size 11 --outfmt 6
```

### Search against a FASTA file (no database needed)

```bash
blast-cli blastn --query query.fa --subject reference.fa
```

### Write output to a file

```bash
blast-cli blastn --query query.fa --db mydb --out results.tsv --outfmt 6
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
| `--outfmt` | `6` | Output format. `blastn` supports 0=pairwise, 5=XML, 6=tabular, 7=commented tabular, 10=CSV, 17=SAM. Other CLI programs currently support tabular/CSV only. |
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

Protein and translated CLI programs currently support tabular output (`--outfmt`
`6` or `10`). Pairwise, XML, commented tabular, and SAM writers are only wired
for `blastn`; other programs fail explicitly for those formats until parity
writers are implemented.

## Library Usage

Add the crate to your `Cargo.toml`. The package name is `blast-rs`, while the
Rust import path is `blast_rs`.

```toml
[dependencies]
blast-rs = "0.12"
```

For local development against this checkout:

```toml
[dependencies]
blast-rs = { path = "/home/mahogny/github/claude/newblast" }
```

The library API returns structured Rust values instead of formatted BLAST text.
Open a database once, reuse it for many queries, and format results yourself or
with `blast_rs::format`.

### Search an existing nucleotide database

```rust
use blast_rs::{blastn, BlastDb, SearchParams};
use std::path::Path;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Pass the BLAST database prefix, not a specific .nin/.nsq/.nhr file.
    let db = BlastDb::open(Path::new("mydb"))?;

    let params = SearchParams::blastn()
        .evalue(1e-20)
        .strand("both")
        .filter_low_complexity(false);

    let query = b"ACGTACGTACGTACGTACGTACGTACGT";
    let results = blastn(&db, query, &params);

    for hit in &results {
        for hsp in &hit.hsps {
            println!(
                "{}\t{}..{}\t{}..{}\t{:.3}\t{:.2e}",
                hit.subject_accession,
                hsp.query_start,
                hsp.query_end,
                hsp.subject_start,
                hsp.subject_end,
                hsp.percent_identity(),
                hsp.evalue
            );
        }
    }

    Ok(())
}
```

### High-level search API

```rust
use blast_rs::{blastp, BlastDb, BlastDbBuilder, SearchParams, SequenceEntry};
use blast_rs::db::DbType;
use std::path::Path;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Build a protein database programmatically.
    let mut builder = BlastDbBuilder::new(DbType::Protein, "my db");
    builder.add(SequenceEntry {
        title: "target protein".into(),
        accession: "P001".into(),
        sequence: b"MKFLILLFNILCLFPVLAADNHGVSMNAS".to_vec(),
        taxid: None,
    });
    builder.write(Path::new("mydb"))?;

    let db = BlastDb::open(Path::new("mydb"))?;
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false);

    let results = blastp(&db, b"MKFLILLFNILCLFPVLAADNHGVSMNAS", &params);
    for r in &results {
        for hsp in &r.hsps {
            println!(
                "score={} evalue={:.2e} identity={:.1}%",
                hsp.score,
                hsp.evalue,
                hsp.percent_identity()
            );
        }
    }

    Ok(())
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
| `blast_rs::matrix` | Scoring matrices (BLOSUM45/50/62/80/90, PAM30/70/250, nucleotide) |
| `blast_rs::db` | BLAST database reader/writer (v4/v5), taxonomy ID/name lookups (`.nto`, `taxdb`) |
| `blast_rs::input` | FASTA parser (via noodles-fasta) with BLASTNA/NCBIstdaa encoding |
| `blast_rs::format` | Output formatting (tabular, pairwise, XML, SAM) |
| `blast_rs::filter` | DUST low-complexity masking |
| `blast_rs::util` | Six-frame translation, genetic code, coordinate mapping |

## Benchmarks

These numbers are only meant to track this translation against the local NCBI
BLAST+ 2.17.0 reference binary. They are not general performance claims.

Benchmark setup:

- Rust binary: `target/release/blast-cli`
- Reference binary: `ncbi-blast-2.17.0+-src/c++/ReleaseMT/bin/blastn`
- Build used for this run: release build, no LTO
- Host date: 2026-04-26
- Threads: 1
- Output: tabular file output, then `cmp` against the reference output
- Timing: median wall-clock seconds from 3 runs, using `/usr/bin/time`

Command shape:

```bash
blast-cli blastn --task blastn-short --dust no --evalue 10 \
    --query <query.fa> --db <db> \
    --outfmt '6 qseqid sseqid pident length qstart qend sstart send bitscore evalue' \
    --max_hsps 1 --num_threads 1 --out rust.tsv

blastn -task blastn-short -dust no -evalue 10 \
    -query <query.fa> -db <db> \
    -outfmt '6 qseqid sseqid pident length qstart qend sstart send bitscore evalue' \
    -max_hsps 1 -num_threads 1 -out ncbi.tsv
```

| Case | Query | Database | Output parity | Rows / bytes | Rust median | NCBI median | Speedup |
|------|-------|----------|---------------|--------------|-------------|-------------|---------|
| Short query vs yeast | `tests/fixtures/large_db/query_500.fa` | `tests/fixtures/large_db/yeast` | byte-identical (`cmp=0`) | 17 / 1260 | 0.33 s | 0.38 s | 1.15x |
| Short query vs C. elegans | `tests/fixtures/large_db/query_500.fa` | `tests/fixtures/large_db/celegans` | byte-identical (`cmp=0`) | 6 / 466 | 4.48 s | 5.55 s | 1.24x |
| Longer query vs C. elegans | `tests/fixtures/large_db/query_2000.fa` | `tests/fixtures/large_db/celegans` | byte-identical (`cmp=0`) | 7 / 566 | 11.61 s | 14.70 s | 1.27x |

Raw timings:

| Case | Rust runs | NCBI runs |
|------|-----------|-----------|
| Short query vs yeast | 0.34, 0.32, 0.33 s | 0.38, 0.39, 0.37 s |
| Short query vs C. elegans | 4.38, 4.55, 4.48 s | 5.25, 5.55, 5.61 s |
| Longer query vs C. elegans | 11.86, 11.61, 11.14 s | 14.79, 14.55, 14.70 s |


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
- Protein programs (blastp/blastx/tblastn/tblastx) use BLOSUM62 with lookup-table seeding and gapped X-dropoff DP extension; CLI output is currently limited to tabular/CSV formats.

See [docs/parity.md](docs/parity.md) for the current program, output-format,
database-feature, and option support matrix.

## What's not ported

The NCBI BLAST+ 2.17.0 tarball includes a large subset of the NCBI C++ Toolkit (~140MB source). The following are **not** yet included in blast-rs:

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
