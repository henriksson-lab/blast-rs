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
//! - [`mod@format`] — Output formatting (tabular, pairwise, XML, SAM)
//! - [`search`] — Core search engine
//! - [`stat`] — Karlin-Altschul statistics
//! - [`filter`] — DUST low-complexity masking
//! - [`protein`] — Protein search utilities
//! - [`matrix`] — Scoring matrices (BLOSUM62, nucleotide)
//! - [`util`] — Translation utilities (six-frame, genetic code)

#![allow(clippy::doc_overindented_list_items)]
#![allow(clippy::duplicated_attributes)]
#![allow(clippy::field_reassign_with_default)]
#![allow(clippy::if_same_then_else)]
#![allow(clippy::identity_op)]
#![allow(clippy::iter_cloned_collect)]
#![allow(clippy::manual_div_ceil)]
#![allow(clippy::manual_is_multiple_of)]
#![allow(clippy::manual_range_contains)]
#![allow(clippy::mut_range_bound)]
#![allow(clippy::needless_borrow)]
#![allow(clippy::ptr_arg)]
#![allow(clippy::question_mark)]
#![allow(clippy::result_unit_err)]
#![allow(clippy::type_complexity)]
#![allow(clippy::unnecessary_cast)]
#![allow(clippy::unnecessary_map_or)]
#![allow(clippy::unusual_byte_groupings)]

// Public API modules
pub mod api;
pub mod blast_kappa;
pub mod blast_setup;
pub mod compo_mode_condition;
pub mod composition;
pub mod hspfilter_culling;
pub mod nlm_linear_algebra;
pub mod optimize_target_freq;
pub mod semi_gapped_align;
pub mod smith_waterman;

// Re-export high-level API at crate root for convenience
pub use api::{
    adjust_evalue,
    adjust_evalue_with_mode,
    // Masking
    apply_dust,
    apply_lowercase_mask_nucleotide,
    apply_lowercase_mask_protein,
    apply_repeat_mask,
    apply_seg,
    apply_seg_ncbistdaa,
    // Search functions (low-level aliases)
    blast_search,
    blastn,
    blastn_search,
    // Search functions (high-level)
    blastp,
    blastp_batch,
    blastx,
    blastx_search,
    // PSSM functions
    build_pssm,
    // Composition
    composition_ncbistdaa,
    get_codon_table,
    get_matrix,
    lowercase_mask,
    // Utility functions
    parse_fasta,
    psiblast,
    psiblast_search,
    repeat_mask,
    reverse_complement,
    search_with_pssm,
    six_frame_translate,
    six_frame_translate_with_code,
    tblastn,
    tblastn_search,
    tblastx,
    tblastx_search,
    BlastDbBuilder,
    BlastDefLine,
    BlastnSearch,
    Hsp,
    MatrixType,
    PsiblastParams,
    ScoringMatrix,
    // Types
    SearchParams,
    SearchResult,
    SequenceEntry,
    Strand,
    TranslatedFrame,
};

// Re-export core types at crate root (matching old API)
pub use db::BlastDb;
pub use link_hsps::{
    BLAST_LinkHsps, LinkBlastHsp, LinkBlastHspList, LinkBlastSeg, LinkHSPParameters, LinkScoreBlock,
};
pub use matrix::AA_FREQUENCIES as BACKGROUND_FREQ;
pub use program::BLASTN;
pub use pssm::Pssm;
pub use queryinfo::QueryInfo;
pub use stat::KarlinBlk as KarlinAltschul;
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
pub mod gapinfo;
#[allow(dead_code)]
pub(crate) mod hspstream;
#[allow(dead_code)]
pub(crate) mod itree;
#[allow(dead_code)]
pub(crate) mod options;

// Internal modules — ported from C engine, not yet fully wired up.
// Kept for completeness and future use (e.g. full protein search, engine orchestration).
#[allow(dead_code)]
pub(crate) mod diagnostics;
pub mod encoding;
#[allow(dead_code)]
pub(crate) mod engine;
#[allow(dead_code)]
pub(crate) mod extend;
#[allow(dead_code)]
pub(crate) mod gapalign;
#[allow(dead_code)]
pub(crate) mod greedy;
#[allow(dead_code)]
pub(crate) mod hits;
#[allow(dead_code)]
pub(crate) mod link_hsps;
#[allow(dead_code)]
pub(crate) mod listnode;
#[allow(dead_code)]
pub(crate) mod lookup;
pub mod math;
#[allow(dead_code)]
pub(crate) mod parameters;
#[allow(dead_code)]
pub(crate) mod program;
#[allow(dead_code)]
pub(crate) mod queryinfo;
#[allow(dead_code)]
pub(crate) mod rps;
#[allow(dead_code)]
pub(crate) mod seqsrc;
#[allow(dead_code)]
pub(crate) mod sequence;
