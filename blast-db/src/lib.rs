//! Pure Rust reader for BLAST database files (.nin/.nsq/.nhr, .pin/.psq/.phr).

pub use blast_core_sys as ffi;

pub mod alias;
mod index;
pub mod makedb;
mod seqsrc;

pub use index::{BlastDb, DbType};
pub use seqsrc::create_seq_src;
