//! Rust equivalent of blast_diagnostics.c — search diagnostics and statistics.
//! Tracks counts of lookup hits, extensions, and other search phases.

/// Diagnostics from the ungapped search phase.
#[derive(Debug, Clone, Default)]
pub struct UngappedStats {
    /// Number of successful lookup table hits
    pub lookup_hits: i64,
    /// Number of sequences with at least one lookup hit
    pub num_seqs_lookup_hits: i32,
    /// Number of initial words found and extended
    pub init_extends: i32,
    /// Number of successful initial extensions (HSPs saved)
    pub good_init_extends: i32,
    /// Number of sequences with at least one HSP after ungapped stage
    pub num_seqs_passed: i32,
}

/// Diagnostics from the gapped search phase.
#[derive(Debug, Clone, Default)]
pub struct GappedStats {
    /// Number of sequences with gapped extensions
    pub num_seqs_passed: i32,
    /// Number of gapped extensions performed
    pub extensions: i32,
    /// Number of gapped extensions producing significant HSPs
    pub good_extensions: i32,
}

/// Complete search diagnostics.
#[derive(Debug, Clone, Default)]
pub struct Diagnostics {
    pub ungapped: UngappedStats,
    pub gapped: GappedStats,
}

impl Diagnostics {
    pub fn new() -> Self {
        Diagnostics::default()
    }

    /// Record a lookup hit.
    pub fn add_lookup_hit(&mut self) {
        self.ungapped.lookup_hits += 1;
    }

    /// Record a successful ungapped extension.
    pub fn add_ungapped_hit(&mut self) {
        self.ungapped.good_init_extends += 1;
    }

    /// Record a successful gapped extension.
    pub fn add_gapped_hit(&mut self) {
        self.gapped.good_extensions += 1;
    }

    /// Summary string for logging.
    pub fn summary(&self) -> String {
        format!(
            "Lookup hits: {}, Ungapped extensions: {}/{}, Gapped extensions: {}/{}",
            self.ungapped.lookup_hits,
            self.ungapped.good_init_extends,
            self.ungapped.init_extends,
            self.gapped.good_extensions,
            self.gapped.extensions,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_diagnostics() {
        let mut diag = Diagnostics::new();
        diag.add_lookup_hit();
        diag.add_lookup_hit();
        diag.add_ungapped_hit();
        assert_eq!(diag.ungapped.lookup_hits, 2);
        assert_eq!(diag.ungapped.good_init_extends, 1);
        assert!(!diag.summary().is_empty());
    }

    #[test]
    fn test_diagnostics_default_zeros() {
        let diag = Diagnostics::new();
        assert_eq!(diag.ungapped.lookup_hits, 0);
        assert_eq!(diag.ungapped.num_seqs_lookup_hits, 0);
        assert_eq!(diag.ungapped.init_extends, 0);
        assert_eq!(diag.ungapped.good_init_extends, 0);
        assert_eq!(diag.ungapped.num_seqs_passed, 0);
        assert_eq!(diag.gapped.num_seqs_passed, 0);
        assert_eq!(diag.gapped.extensions, 0);
        assert_eq!(diag.gapped.good_extensions, 0);
    }

    #[test]
    fn test_diagnostics_gapped_tracking() {
        let mut diag = Diagnostics::new();
        diag.add_gapped_hit();
        diag.add_gapped_hit();
        diag.add_gapped_hit();
        assert_eq!(diag.gapped.good_extensions, 3);
        // Ungapped should remain zero
        assert_eq!(diag.ungapped.good_init_extends, 0);
    }

    #[test]
    fn test_diagnostics_mixed_tracking() {
        let mut diag = Diagnostics::new();
        for _ in 0..100 {
            diag.add_lookup_hit();
        }
        for _ in 0..10 {
            diag.add_ungapped_hit();
        }
        for _ in 0..3 {
            diag.add_gapped_hit();
        }
        assert_eq!(diag.ungapped.lookup_hits, 100);
        assert_eq!(diag.ungapped.good_init_extends, 10);
        assert_eq!(diag.gapped.good_extensions, 3);
    }

    #[test]
    fn test_diagnostics_summary_format() {
        let mut diag = Diagnostics::new();
        diag.add_lookup_hit();
        diag.ungapped.init_extends = 5;
        diag.add_ungapped_hit();
        diag.gapped.extensions = 2;
        diag.add_gapped_hit();
        let s = diag.summary();
        assert!(s.contains("Lookup hits: 1"));
        assert!(s.contains("Ungapped extensions: 1/5"));
        assert!(s.contains("Gapped extensions: 1/2"));
    }

    #[test]
    fn test_ungapped_stats_clone() {
        let mut stats = UngappedStats::default();
        stats.lookup_hits = 42;
        stats.num_seqs_passed = 7;
        let cloned = stats.clone();
        assert_eq!(cloned.lookup_hits, 42);
        assert_eq!(cloned.num_seqs_passed, 7);
    }

    #[test]
    fn test_gapped_stats_clone() {
        let mut stats = GappedStats::default();
        stats.extensions = 10;
        stats.good_extensions = 3;
        let cloned = stats.clone();
        assert_eq!(cloned.extensions, 10);
        assert_eq!(cloned.good_extensions, 3);
    }

    #[test]
    fn test_diagnostics_direct_field_mutation() {
        let mut diag = Diagnostics::new();
        diag.ungapped.num_seqs_lookup_hits = 5;
        diag.ungapped.num_seqs_passed = 3;
        diag.gapped.num_seqs_passed = 2;
        diag.gapped.extensions = 8;
        assert_eq!(diag.ungapped.num_seqs_lookup_hits, 5);
        assert_eq!(diag.ungapped.num_seqs_passed, 3);
        assert_eq!(diag.gapped.num_seqs_passed, 2);
        assert_eq!(diag.gapped.extensions, 8);
    }
}
