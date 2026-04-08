//! High-level blastn search API with builder pattern.

use crate::search::{SearchHsp, blastn_gapped_search};
use crate::stat::{KarlinBlk, ungapped_kbp_calc, nucl_gapped_kbp_lookup,
                   nucl_alpha_beta, compute_length_adjustment_exact, UngappedKbpContext};
use crate::filter;
use crate::traceback::build_blastna_matrix;

/// Builder for configuring and running a blastn search.
///
/// # Example
///
/// ```no_run
/// use blast_core::blastn::BlastnSearch;
///
/// let results = BlastnSearch::new()
///     .query(b"ACGTACGTACGT")
///     .subject(b"NNNNACGTACGTACGTNNNN")
///     .run();
///
/// for hit in &results {
///     println!("score={} evalue={:.2e}", hit.score, hit.evalue);
/// }
/// ```
pub struct BlastnSearch {
    pub word_size: usize,
    pub reward: i32,
    pub penalty: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub evalue: f64,
    pub dust: bool,
    pub strand: Strand,
    pub xdrop_gap_final: f64,
    query_raw: Vec<u8>,
    subject_raw: Vec<u8>,
}

/// Which query strand(s) to search.
#[derive(Clone, Copy, PartialEq)]
pub enum Strand {
    Both,
    Plus,
    Minus,
}

impl Default for BlastnSearch {
    fn default() -> Self {
        Self::new()
    }
}

impl BlastnSearch {
    /// Create a new search with default blastn parameters.
    pub fn new() -> Self {
        BlastnSearch {
            word_size: 11,
            reward: 1,
            penalty: -3,
            gap_open: 5,
            gap_extend: 2,
            evalue: 10.0,
            dust: true,
            strand: Strand::Both,
            xdrop_gap_final: 100.0,
            query_raw: Vec::new(),
            subject_raw: Vec::new(),
        }
    }

    /// Set the query sequence (IUPAC nucleotide, e.g. b"ACGTACGT").
    pub fn query(mut self, seq: &[u8]) -> Self {
        self.query_raw = seq.to_vec();
        self
    }

    /// Set the subject sequence (IUPAC nucleotide).
    pub fn subject(mut self, seq: &[u8]) -> Self {
        self.subject_raw = seq.to_vec();
        self
    }

    /// Set word size for initial seed (default: 11).
    pub fn word_size(mut self, ws: usize) -> Self {
        self.word_size = ws;
        self
    }

    /// Set match reward (default: 1).
    pub fn reward(mut self, r: i32) -> Self {
        self.reward = r;
        self
    }

    /// Set mismatch penalty (default: -3, should be negative).
    pub fn penalty(mut self, p: i32) -> Self {
        self.penalty = p;
        self
    }

    /// Set gap opening cost (default: 5).
    pub fn gap_open(mut self, g: i32) -> Self {
        self.gap_open = g;
        self
    }

    /// Set gap extension cost (default: 2).
    pub fn gap_extend(mut self, g: i32) -> Self {
        self.gap_extend = g;
        self
    }

    /// Set E-value threshold (default: 10.0).
    pub fn evalue(mut self, e: f64) -> Self {
        self.evalue = e;
        self
    }

    /// Enable/disable DUST low-complexity filtering (default: true).
    pub fn dust(mut self, d: bool) -> Self {
        self.dust = d;
        self
    }

    /// Set which strand(s) to search (default: Both).
    pub fn strand(mut self, s: Strand) -> Self {
        self.strand = s;
        self
    }

    /// Run the search and return hits.
    pub fn run(&self) -> Vec<SearchHsp> {
        if self.query_raw.is_empty() || self.subject_raw.is_empty() {
            return Vec::new();
        }

        // Encode query
        let mut query_plus: Vec<u8> = self.query_raw.iter()
            .map(|&b| iupac_to_blastna(b))
            .collect();

        // Apply DUST
        if self.dust {
            let mask = filter::dust_filter(&query_plus, 64, 2.0);
            mask.apply(&mut query_plus, 14);
        }

        let query_plus_nomask: Vec<u8> = self.query_raw.iter()
            .map(|&b| iupac_to_blastna(b))
            .collect();

        let query_minus: Vec<u8> = query_plus.iter().rev()
            .map(|&b| complement_blastna(b))
            .collect();
        let query_minus_nomask: Vec<u8> = query_plus_nomask.iter().rev()
            .map(|&b| complement_blastna(b))
            .collect();

        let qp = if self.strand != Strand::Minus { &query_plus[..] } else { &[] };
        let qm = if self.strand != Strand::Plus { &query_minus[..] } else { &[] };
        let _qpn = if self.strand != Strand::Minus { &query_plus_nomask[..] } else { &[] };
        let _qmn = if self.strand != Strand::Plus { &query_minus_nomask[..] } else { &[] };

        // Encode subject
        let subject: Vec<u8> = self.subject_raw.iter()
            .map(|&b| iupac_to_blastna(b))
            .collect();

        // Compute KBP
        let matrix_fn = |i: usize, j: usize| -> i32 {
            let m = build_blastna_matrix(self.reward, self.penalty);
            m[i][j]
        };
        let mut lo = i32::MAX;
        let mut hi = i32::MIN;
        let m = build_blastna_matrix(self.reward, self.penalty);
        for i in 0..16 { for j in 0..16 {
            let s = m[i][j];
            if s <= -100000000 || s >= 100000000 { continue; }
            if s < lo { lo = s; }
            if s > hi { hi = s; }
        }}

        let ambig: &[u8] = &[14, 15];
        let ctx = UngappedKbpContext { query_offset: 0, query_length: query_plus.len() as i32, is_valid: true };
        let kbp_results = ungapped_kbp_calc(&query_plus, &[ctx], lo, hi, 16, ambig, &matrix_fn);
        let ungapped_kbp = kbp_results[0].clone().unwrap_or(KarlinBlk {
            lambda: 1.374, k: 0.621, log_k: 0.621_f64.ln(), h: 1.286,
        });
        let (gapped_kbp, _) = nucl_gapped_kbp_lookup(
            self.gap_open, self.gap_extend, self.reward, self.penalty, &ungapped_kbp,
        ).unwrap_or((ungapped_kbp.clone(), false));

        // Compute search space
        let (alpha, beta) = nucl_alpha_beta(
            self.reward, self.penalty, self.gap_open, self.gap_extend,
            ungapped_kbp.lambda, ungapped_kbp.h, true,
        );
        let db_length = self.subject_raw.len() as i64;
        let (len_adj, _) = compute_length_adjustment_exact(
            gapped_kbp.k, gapped_kbp.log_k,
            alpha / gapped_kbp.lambda, beta,
            self.query_raw.len() as i32, db_length, 1,
        );
        let eff_db = (db_length - len_adj as i64).max(1);
        let search_space = eff_db as f64 * (self.query_raw.len() as i32 - len_adj).max(1) as f64;

        let x_dropoff = (self.xdrop_gap_final * std::f64::consts::LN_2 / gapped_kbp.lambda) as i32;

        blastn_gapped_search(
            qp, qm, &subject,
            self.word_size, self.reward, self.penalty,
            self.gap_open, self.gap_extend, x_dropoff,
            &gapped_kbp, search_space, self.evalue,
        )
    }
}

fn iupac_to_blastna(b: u8) -> u8 {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        b'R' | b'r' => 4,
        b'Y' | b'y' => 5,
        b'M' | b'm' => 6,
        b'K' | b'k' => 7,
        b'W' | b'w' => 8,
        b'S' | b's' => 9,
        b'B' | b'b' => 10,
        b'D' | b'd' => 11,
        b'H' | b'h' => 12,
        b'V' | b'v' => 13,
        b'N' | b'n' => 14,
        _ => 15,
    }
}

fn complement_blastna(b: u8) -> u8 {
    match b {
        0 => 3, 1 => 2, 2 => 1, 3 => 0,
        4 => 5, 5 => 4, 6 => 7, 7 => 6,
        8 => 8, 9 => 9, 10 => 13, 11 => 12,
        12 => 11, 13 => 10, 14 => 14, _ => 15,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_builder_defaults() {
        let s = BlastnSearch::new();
        assert_eq!(s.word_size, 11);
        assert_eq!(s.reward, 1);
        assert_eq!(s.penalty, -3);
        assert_eq!(s.evalue, 10.0);
    }

    #[test]
    fn test_perfect_match() {
        let results = BlastnSearch::new()
            .word_size(7)
            .dust(false)
            .query(b"ACGTACGTACGTACGTACGTACGT")
            .subject(b"NNNNNACGTACGTACGTACGTACGTNNNNN")
            .run();
        assert!(!results.is_empty(), "Should find perfect match");
        assert!(results[0].score > 0);
    }

    #[test]
    fn test_no_match() {
        let results = BlastnSearch::new()
            .word_size(7)
            .dust(false)
            .strand(Strand::Plus)
            .query(b"AAAAAAAAAAAAAAAA")
            .subject(b"CCCCCCCCCCCCCCCC")
            .run();
        assert!(results.is_empty(), "Should find no match");
    }

    #[test]
    fn test_custom_scoring() {
        let results = BlastnSearch::new()
            .word_size(7)
            .reward(2)
            .penalty(-3)
            .gap_open(5)
            .gap_extend(2)
            .dust(false)
            .query(b"ACGTACGTACGTACGT")
            .subject(b"ACGTACGTACGTACGT")
            .run();
        assert!(!results.is_empty());
    }
}
