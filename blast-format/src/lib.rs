//! BLAST output formatting.

pub use blast_core_sys as ffi;

mod pairwise;
mod sam;
mod tabular;
mod xml;

pub use pairwise::format_pairwise_alignment;
pub use sam::{write_sam_header, write_sam_record};
pub use tabular::{format_tabular, format_tabular_custom, TabularHit};
pub use xml::{write_xml_header, write_xml_hit, write_xml_footer};
