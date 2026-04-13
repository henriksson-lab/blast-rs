//! blast-rs: Pure-Rust implementation of NCBI BLAST.
//!
//! # Public API
//!
//! The high-level API is in the [`api`] module. For convenience, the most
//! commonly used types and functions are re-exported at the crate root.
//!
//! - [`api`] — High-level search functions (`blastp`, `blastn`, `blastx`, etc.)
//! - [`blastn`] — Builder API for blastn searches
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
pub mod api;
pub mod composition;
pub mod compo_mode_condition;
pub mod nlm_linear_algebra;
pub mod optimize_target_freq;
pub mod semi_gapped_align;

// Re-export high-level API at crate root for convenience
pub use api::{
    // Search functions (high-level)
    blastp, blastp_batch, blastn, blastn_search, blastx, tblastn, tblastx, psiblast,
    // Search functions (low-level aliases)
    blast_search, blastx_search, tblastn_search, tblastx_search, psiblast_search,
    // PSSM functions
    build_pssm, search_with_pssm,
    // Types
    SearchParams, SearchResult, Hsp, MatrixType, ScoringMatrix,
    BlastDbBuilder, SequenceEntry, PsiblastParams, TranslatedFrame,
    BlastDefLine, BlastnSearch, Strand,
    // Utility functions
    parse_fasta, reverse_complement, six_frame_translate, six_frame_translate_with_code,
    get_codon_table, get_matrix,
    // Masking
    apply_dust, apply_seg, apply_seg_ncbistdaa,
    apply_lowercase_mask_nucleotide, apply_lowercase_mask_protein,
    lowercase_mask, apply_repeat_mask, repeat_mask,
    // Composition
    composition_ncbistdaa,
    adjust_evalue, adjust_evalue_with_mode,
};

// Re-export core types at crate root (matching old API)
pub use db::BlastDb;
pub use pssm::Pssm;
pub use stat::KarlinBlk as KarlinAltschul;
pub use matrix::AA_FREQUENCIES as BACKGROUND_FREQ;
// blastn builder (BlastnSearch) is now in the api module.
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
#[allow(dead_code)] pub mod gapinfo;
#[allow(dead_code)] pub(crate) mod hspstream;
#[allow(dead_code)] pub(crate) mod itree;
#[allow(dead_code)] pub(crate) mod options;

// Internal modules — ported from C engine, not yet fully wired up.
// Kept for completeness and future use (e.g. full protein search, engine orchestration).
#[allow(dead_code)] pub(crate) mod diagnostics;
#[allow(dead_code)] pub mod encoding;
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
#[allow(dead_code)] pub(crate) mod seqsrc;
