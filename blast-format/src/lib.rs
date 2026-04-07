//! BLAST output formatting.

pub use blast_core_sys as ffi;

mod pairwise;
mod tabular;

pub use pairwise::format_pairwise_alignment;
pub use tabular::{format_tabular, TabularHit};
