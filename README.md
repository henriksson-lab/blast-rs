# blast-rs

Pure-Rust implementation of NCBI BLAST (Basic Local Alignment Search Tool). Produces byte-identical output to the C reference implementation (NCBI BLAST+ 2.17.0) for blastn, with no C dependencies. Supports blastn, blastp, blastx, tblastn, tblastx, and psiblast.

Based on [NCBI BLAST+ 2.17.0](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.17.0/) source distribution (`ncbi-blast-2.17.0+-src.tar.gz`, also on [GitHub](https://github.com/ncbi/blast)). The C core algorithms in `src/algo/blast/core/` were ported function-by-function to Rust.

This is a translation of the original code and not the authoritative implementation. This code should generate bitwise
equal output to the original. Please report any deviations.

The aim of this project is to increase performance, especially by providing this code through a type-safe library interface.
The code can also be compiled to be used for webassembly.

**Remaining open issues in ported code:**

  Remaining:
  1. Composition-based statistics performance — comp_adjust=2 (NCBI default) is extremely slow, had to disable it. Needs profiling/optimization
  2. Combined multi-query lookup table — scan subject truly once for all queries (would close the remaining multi-query gap: 0.68s vs NCBI's 0.41s)
  3. blastp finds extra hits — our ungapped extension scores differ by 1-2 points from NCBI, catching slightly more seeds 


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
- **Low startup overhead** -- near-zero startup vs ~400ms for NCBI BLAST+ (see benchmarks)
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
- **494 passing release tests**: 353 library unit tests + 13 CLI unit tests + 123 integration tests + 4 stress tests + doc tests, plus ignored real-database parity tests

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
| `blast_rs::matrix` | Scoring matrices (BLOSUM45/50/62/80/90, PAM30/70/250, nucleotide) |
| `blast_rs::db` | BLAST database reader/writer (v4/v5), taxonomy ID/name lookups (`.nto`, `taxdb`) |
| `blast_rs::input` | FASTA parser (via noodles-fasta) with BLASTNA/NCBIstdaa encoding |
| `blast_rs::format` | Output formatting (tabular, pairwise, XML, SAM) |
| `blast_rs::filter` | DUST low-complexity masking |
| `blast_rs::util` | Six-frame translation, genetic code, coordinate mapping |

## Benchmarks

Wall-clock time including process startup and database loading, measured on Linux x86-64.

### Current core_nt primer benchmark

Measured 2026-04-14/15 on Linux x86-64 with a native CPU release build:

```bash
RUSTFLAGS="-C target-cpu=native" cargo build --release
target/release/blast-cli blastn \
    --query /tmp/core_nt_realistic_primer.fa \
    --db /husky/henriksson/for_claude/blast/core_nt/core_nt.00 \
    --task blastn-short \
    --outfmt 6 \
    --max_target_seqs 500 \
    --num_threads 8
```

Taxonomy benchmark command:

```bash
target/release/blast-cli blastn \
    --task blastn-short \
    --evalue 5 \
    --max_target_seqs 10000 \
    --max_hsps 1 \
    --num_threads 8 \
    --perc_identity 90 \
    --db /husky/henriksson/for_claude/blast/core_nt/core_nt.00 \
    --outfmt '6 qseqid qlen qstart qend saccver slen sstart send bitscore staxid ssciname sskingdom length pident' \
    --query /tmp/core_nt_realistic_primer.fa \
    --out /tmp/primer.vsCore_nt.rsblastn
```

Query sequence:

```text
GTCTCCTCTGACTTCAACAGCG
```

| Test | Database | NCBI BLAST+ 2.17.0 | blast-rs native release | Output |
|------|----------|--------------------:|------------------------:|--------|
| 23bp primer, `blastn-short`, 8 threads, `outfmt 6` | `core_nt.00` | 2.04-2.05 s warm | 1.98-1.99 s warm | byte-identical |
| Same primer, taxonomy outfmt with `staxid ssciname sskingdom` | `core_nt.00` | 4.60 s with `BLASTDB` taxdb, ~5.59 s without | 2.01 s with `BLASTDB` taxdb, 2.00 s without | byte-identical |
| Same primer, taxonomy outfmt with `staxid ssciname sskingdom` | standalone `core_nt.12` and `core_nt.28` | matched in release ignored regression | matched in release ignored regression | byte-identical |
| Same primer, taxonomy outfmt with `staxid ssciname sskingdom` | full `core_nt` alias | 1217.89 s | 455.80 s | byte-identical |

On the first large real database volume, blast-rs is now roughly at parity to slightly faster than NCBI BLAST+ for the tested plain primer scan while producing byte-for-byte identical tabular output. This plain `outfmt 6` case is the better comparison for raw search-path speed.

The much larger taxonomy-output speedup should not be interpreted as Rust being intrinsically several times faster at BLAST alignment. That command also measures database metadata and taxonomy lookup behavior. For taxonomy-heavy output, blast-rs uses mmap-backed `.not/.pot` OID-to-taxid lookup, loads it only when taxonomy fields are requested, and touches only the OID index entries and taxid/name records needed by the returned hits. In the measured full `core_nt` alias run, NCBI BLAST+ used about 182 GB RSS versus about 4.4 GB for blast-rs, so the 455.80 s vs 1217.89 s gap is likely dominated by taxonomy/metadata working-set behavior rather than by the alignment scanner alone.

Standalone nonzero volumes use the alias DBLIST to translate local OIDs into the global OID namespace before indexing the alias-level `.not` file.

Ignored real-database parity tests for these cases:

```bash
RUSTFLAGS="-C target-cpu=native" cargo test --release -- --ignored test_core_nt00_primer_taxonomy_outfmt_matches_ncbi
RUSTFLAGS="-C target-cpu=native" cargo test --release -- --ignored test_core_nt_nonzero_volume_taxonomy_outfmt_matches_ncbi
RUSTFLAGS="-C target-cpu=native" cargo test --release -- --ignored test_core_nt00_primer_taxonomy_names_with_blastdb_match_ncbi
BLAST_RS_RUN_FULL_CORE_NT=1 RUSTFLAGS="-C target-cpu=native" cargo test --release -- --ignored test_core_nt_alias_primer_taxonomy_outfmt_matches_ncbi
```

The following older benchmark tables are single-threaded BLAST database mode (`--db`), median of 3 runs.

### blastn — single query

| Test | Database | NCBI BLAST+ 2.17.0 | blast-rs 0.10.0 |
|------|----------|--------------------:|----------------:|
| 200bp query | C. elegans (100 MB) | 411 ms | 786 ms |
| 200bp query | S. cerevisiae (12 MB) | 394 ms | 126 ms |
| 200bp query | S. pombe (2.5 MB) | 390 ms | 49 ms |
| 200bp query | seqn (1 MB, 2K seqs) | 392 ms | 68 ms |
| 1000bp query | C. elegans (100 MB) | 413 ms | 943 ms |

### blastn-short — primer search

| Test | Database | NCBI BLAST+ 2.17.0 | blast-rs 0.10.0 |
|------|----------|--------------------:|----------------:|
| 30bp primer | C. elegans (100 MB) | 493 ms | 723 ms |
| 19bp primer | 16S rRNA (10 MB, 27K seqs) | 1254 ms | 923 ms |

### blastn — multi-query

| Test | Database | NCBI BLAST+ 2.17.0 | blast-rs 0.10.0 |
|------|----------|--------------------:|----------------:|
| 10x 200bp queries | C. elegans (100 MB) | 420 ms | 8000 ms |
| 20 primers | 16S rRNA (10 MB, 27K seqs) | 598 ms | 7226 ms |

**Notes:**
- NCBI BLAST+ has ~390 ms fixed startup overhead (shared library loading, ASN.1 initialization). On small databases, NCBI wall-clock time is dominated by startup, not search.
- blast-rs has near-zero startup (<10 ms), so it is faster on wall-clock for small-to-medium databases with single queries.
- For large databases (C. elegans 100 MB) and multi-query workloads, NCBI BLAST+ is faster at raw search throughput. The inner scanning loop needs further optimization.
- Multi-threaded search (`--num_threads`) parallelizes across database sequences using rayon.

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

### 0.10.0

- **BLAST database v5 support**: Fixed `.nto` taxonomy file parser for v5 databases (used by all modern NCBI databases including core_nt, 16S_ribosomal_RNA, etc.)
- **Composition-based statistics**: Ported `composition_adjustment/` for blastp/blastx compositional score correction (matrix scaling, relative entropy optimization, lambda ratio adjustment)
- **Stack overflow fix**: Rayon thread pools now use 64 MB stack size and respect `--num_threads`, fixing crashes on large databases with many threads
- **Stress tests**: Added 4 new stress tests for blastn-short primer search scenarios (2000 subjects, long subjects, 15bp primers, real V5 database)
- Updated benchmarks against NCBI BLAST+ 2.17.0

### 0.9.1

- Maintenance release

### 0.9.0

- **204 new tests** (189 → 393): ported from NCBI BLAST+ unit test suite covering scoring/statistics, lookup tables, scanning, extension, HSP processing, traceback, gapped alignment, PSSM/PSI-BLAST, options, DUST filtering, database I/O, output formatting, encoding, and end-to-end integration
- Updated benchmarks with larger databases (C. elegans 100 MB, S. cerevisiae 12 MB) and wall-clock timing

### 0.8.0

- **33x blastp speedup**: Protein lookup table is now built once per query instead of per subject, eliminating redundant work when searching large databases
- **Gapped protein alignment**: blastp now performs two-phase search (ungapped seeds → Smith-Waterman gapped alignment), matching NCBI BLAST+ sensitivity for distant homologs
- **Smith-Waterman local alignment**: Protein gapped traceback uses local alignment so scores are seed-position-independent, producing optimal alignments matching NCBI BLAST+ scores
- **Corrected protein X-drop defaults**: `x_drop_ungapped` (7→40) and `x_drop_gapped` (38→260) now match NCBI BLAST+ defaults for protein searches
- **Pre-built lookup table API**: New `protein_scan_with_table()` and `protein_gapped_scan_with_table()` for amortizing lookup table construction across multiple subjects
- Inner-loop allocation optimizations in lookup table construction (`suffix_max`, `word_buf` reused across query positions)

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
