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

pub mod algo;

#[path = "../src.orig/api.rs"]
pub mod api;
#[path = "../src.orig/blast_kappa.rs"]
pub mod blast_kappa;
#[path = "../src.orig/blast_setup.rs"]
pub mod blast_setup;
#[path = "../src.orig/compo_mode_condition.rs"]
pub mod compo_mode_condition;
#[path = "../src.orig/composition.rs"]
pub mod composition;
#[path = "../src.orig/db/mod.rs"]
pub mod db;
#[path = "../src.orig/diagnostics.rs"]
pub mod diagnostics;
#[path = "../src.orig/encoding.rs"]
pub mod encoding;
#[path = "../src.orig/engine.rs"]
pub mod engine;
#[path = "../src.orig/extend.rs"]
pub mod extend;
#[path = "../src.orig/filter.rs"]
pub mod filter;
#[path = "../src.orig/format/mod.rs"]
pub mod format;
#[path = "../src.orig/gapinfo.rs"]
pub mod gapinfo;
#[path = "../src.orig/greedy.rs"]
pub mod greedy;
#[path = "../src.orig/hits.rs"]
pub mod hits;
#[path = "../src.orig/hspfilter_culling.rs"]
pub mod hspfilter_culling;
#[path = "../src.orig/hspstream.rs"]
pub mod hspstream;
#[path = "../src.orig/index_ungapped.rs"]
pub mod index_ungapped;
#[path = "../src.orig/input/mod.rs"]
pub mod input;
#[path = "../src.orig/itree.rs"]
pub mod itree;
#[path = "../src.orig/link_hsps.rs"]
pub mod link_hsps;
#[path = "../src.orig/listnode.rs"]
pub mod listnode;
#[path = "../src.orig/lookup.rs"]
pub mod lookup;
#[path = "../src.orig/math.rs"]
pub mod math;
#[path = "../src.orig/matrix.rs"]
pub mod matrix;
#[path = "../src.orig/nlm_linear_algebra.rs"]
pub mod nlm_linear_algebra;
#[path = "../src.orig/optimize_target_freq.rs"]
pub mod optimize_target_freq;
#[path = "../src.orig/options.rs"]
pub mod options;
#[path = "../src.orig/parameters.rs"]
pub mod parameters;
#[path = "../src.orig/pattern.rs"]
pub mod pattern;
#[path = "../src.orig/program.rs"]
pub mod program;
#[path = "../src.orig/protein.rs"]
pub mod protein;
#[path = "../src.orig/protein_lookup.rs"]
pub mod protein_lookup;
#[path = "../src.orig/pssm.rs"]
pub mod pssm;
#[path = "../src.orig/queryinfo.rs"]
pub mod queryinfo;
#[path = "../src.orig/rps.rs"]
pub mod rps;
#[path = "../src.orig/search.rs"]
pub mod search;
#[path = "../src.orig/semi_gapped_align.rs"]
pub mod semi_gapped_align;
#[path = "../src.orig/seqsrc.rs"]
pub mod seqsrc;
#[path = "../src.orig/sequence.rs"]
pub mod sequence;
#[path = "../src.orig/smith_waterman.rs"]
pub mod smith_waterman;
#[path = "../src.orig/spliced_hits.rs"]
pub mod spliced_hits;
#[path = "../src.orig/split_query.rs"]
pub mod split_query;
#[path = "../src.orig/stat.rs"]
pub mod stat;
#[path = "../src.orig/traceback.rs"]
pub mod traceback;
#[path = "../src.orig/util.rs"]
pub mod util;

pub use api::*;
pub use db::{make_db, make_nucleotide_db, make_protein_db, BlastDb, DbType};
pub use hspstream::{BlastHSP, BlastHSPList, BlastSeg};
pub use link_hsps::{blast_link_hsp_list, LinkHSPParameters, LinkScoreBlock};
pub use matrix::AA_FREQUENCIES as BACKGROUND_FREQ;
pub use program::BLASTN;
pub use pssm::Pssm;
pub use queryinfo::QueryInfo;
pub use stat::KarlinBlk as KarlinAltschul;
