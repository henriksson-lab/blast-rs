//! blast-rs: Pure-Rust implementation of NCBI BLAST.
//!
//! # Public API
//!
//! - [`blastn`] — High-level builder API for blastn searches
//! - [`db`] — BLAST database reader/writer
//! - [`input`] — FASTA parser and sequence encoding
//! - [`format`] — Output formatting (tabular, pairwise, XML, SAM)
//! - [`search`] — Core search engine
//! - [`stat`] — Karlin-Altschul statistics
//! - [`filter`] — DUST low-complexity masking
//! - [`protein`] — Protein search utilities
//! - [`matrix`] — Scoring matrices (BLOSUM62, nucleotide)
//! - [`util`] — Translation utilities (six-frame, genetic code)

// Public API modules
pub mod blastn;
pub mod db;
pub mod filter;
pub mod format;
pub mod input;
pub mod matrix;
pub mod protein;
pub mod protein_lookup;
pub mod pssm;
pub mod search;
pub mod stat;
pub mod traceback;
pub mod util;

// Internal modules used by the search engine (partially used — allow dead items
// within them since they're ported from C and not all paths are active yet)
#[allow(dead_code)] pub(crate) mod gapinfo;
#[allow(dead_code)] pub(crate) mod hspstream;
#[allow(dead_code)] pub(crate) mod itree;
#[allow(dead_code)] pub(crate) mod options;

// Internal modules — ported from C engine, not yet fully wired up.
// Kept for completeness and future use (e.g. full protein search, engine orchestration).
#[allow(dead_code)] pub(crate) mod diagnostics;
#[allow(dead_code)] pub(crate) mod encoding;
#[allow(dead_code)] pub(crate) mod engine;
#[allow(dead_code)] pub(crate) mod extend;
#[allow(dead_code)] pub(crate) mod gapalign;
#[allow(dead_code)] pub(crate) mod greedy;
#[allow(dead_code)] pub(crate) mod hits;
#[allow(dead_code)] pub(crate) mod listnode;
#[allow(dead_code)] pub(crate) mod lookup;
#[allow(dead_code)] pub(crate) mod math;
#[allow(dead_code)] pub(crate) mod parameters;
#[allow(dead_code)] pub(crate) mod program;
#[allow(dead_code)] pub(crate) mod queryinfo;
#[allow(dead_code)] pub(crate) mod rps;
#[allow(dead_code)] pub(crate) mod sequence;
