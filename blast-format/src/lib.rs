//! BLAST output formatting.

pub use blast_core_sys as ffi;

mod tabular;

pub use tabular::{format_tabular, TabularHit};
