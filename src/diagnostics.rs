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
            self.ungapped.good_init_extends, self.ungapped.init_extends,
            self.gapped.good_extensions, self.gapped.extensions,
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
}
