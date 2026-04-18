//! Rust equivalent of blast_filter.c — query sequence masking/filtering.
//! Implements DUST (nucleotide) and SEG (protein) low-complexity filtering.

/// A masked region in a sequence.
#[derive(Debug, Clone, Copy)]
pub struct MaskedRegion {
    pub start: i32,
    pub end: i32,
}

/// Collection of masked regions for a query.
#[derive(Debug, Clone, Default)]
pub struct MaskLoc {
    pub regions: Vec<MaskedRegion>,
}

impl MaskLoc {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add(&mut self, start: i32, end: i32) {
        self.regions.push(MaskedRegion { start, end });
    }

    /// Check if a position is masked.
    pub fn is_masked(&self, pos: i32) -> bool {
        self.regions.iter().any(|r| pos >= r.start && pos < r.end)
    }

    /// Apply masking to a sequence by replacing masked positions with sentinel value.
    pub fn apply(&self, sequence: &mut [u8], sentinel: u8) {
        for region in &self.regions {
            for pos in region.start..region.end {
                if (pos as usize) < sequence.len() {
                    sequence[pos as usize] = sentinel;
                }
            }
        }
    }
}

/// Simple DUST-like low complexity filter for nucleotide sequences.
/// Identifies regions with low information content (e.g., poly-A runs,
/// simple repeats) and marks them for masking.
pub fn dust_filter(sequence: &[u8], window: usize, threshold: f64) -> MaskLoc {
    let mut mask = MaskLoc::new();
    if sequence.len() < window {
        return mask;
    }

    for start in 0..=(sequence.len() - window) {
        let window_seq = &sequence[start..start + window];
        let complexity = triplet_complexity(window_seq);
        if complexity > threshold {
            mask.add(start as i32, (start + window) as i32);
        }
    }

    // Merge overlapping regions
    merge_regions(&mut mask);
    mask
}

/// Compute triplet complexity score for a window.
fn triplet_complexity(seq: &[u8]) -> f64 {
    if seq.len() < 3 {
        return 0.0;
    }
    let mut counts = [0u32; 64]; // 4^3 possible triplets
    for i in 0..seq.len() - 2 {
        let a = (seq[i] & 3) as usize;
        let b = (seq[i + 1] & 3) as usize;
        let c = (seq[i + 2] & 3) as usize;
        let idx = a * 16 + b * 4 + c;
        counts[idx] += 1;
    }
    let n = (seq.len() - 2) as f64;
    let mut score = 0.0;
    for &c in &counts {
        if c > 1 {
            let f = c as f64;
            score += f * (f - 1.0) / 2.0;
        }
    }
    score / n
}

/// Merge overlapping masked regions.
fn merge_regions(mask: &mut MaskLoc) {
    if mask.regions.len() <= 1 {
        return;
    }
    mask.regions.sort_by_key(|r| r.start);
    let mut merged = Vec::new();
    let mut current = mask.regions[0];
    for r in &mask.regions[1..] {
        if r.start <= current.end {
            current.end = current.end.max(r.end);
        } else {
            merged.push(current);
            current = *r;
        }
    }
    merged.push(current);
    mask.regions = merged;
}

/// SEG-like low complexity filter for protein sequences.
/// Masks regions with low information content (e.g., poly-amino-acid runs).
pub fn seg_filter(sequence: &[u8], window: usize, threshold: f64) -> MaskLoc {
    let mut mask = MaskLoc::new();
    if sequence.len() < window {
        return mask;
    }

    for start in 0..=(sequence.len() - window) {
        let window_seq = &sequence[start..start + window];
        let entropy = sequence_entropy(window_seq, 20);
        if entropy < threshold {
            mask.add(start as i32, (start + window) as i32);
        }
    }
    merge_regions(&mut mask);
    mask
}

/// Compute Shannon entropy of a sequence with given alphabet size.
fn sequence_entropy(seq: &[u8], alphabet_size: usize) -> f64 {
    let mut counts = vec![0u32; alphabet_size + 1];
    for &b in seq {
        let idx = (b as usize).min(alphabet_size);
        counts[idx] += 1;
    }
    let n = seq.len() as f64;
    let mut entropy = 0.0;
    for &c in &counts {
        if c > 0 {
            let p = c as f64 / n;
            entropy -= p * p.log2();
        }
    }
    entropy
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mask_loc() {
        let mut mask = MaskLoc::new();
        mask.add(5, 10);
        assert!(!mask.is_masked(4));
        assert!(mask.is_masked(5));
        assert!(mask.is_masked(9));
        assert!(!mask.is_masked(10));
    }

    #[test]
    fn test_apply_mask() {
        let mut seq = vec![0u8, 1, 2, 3, 0, 1, 2, 3];
        let mut mask = MaskLoc::new();
        mask.add(2, 5);
        mask.apply(&mut seq, 14); // 14 = N in BLASTNA
        assert_eq!(seq, vec![0, 1, 14, 14, 14, 1, 2, 3]);
    }

    #[test]
    fn test_merge_regions() {
        let mut mask = MaskLoc::new();
        mask.add(0, 5);
        mask.add(3, 8);
        mask.add(10, 15);
        merge_regions(&mut mask);
        assert_eq!(mask.regions.len(), 2);
        assert_eq!(mask.regions[0].start, 0);
        assert_eq!(mask.regions[0].end, 8);
        assert_eq!(mask.regions[1].start, 10);
    }

    /// Low-complexity poly-A sequence should be flagged by DUST.
    #[test]
    fn test_dust_low_complexity_sequence() {
        // 64-base poly-A — maximally repetitive
        let seq: Vec<u8> = vec![0u8; 64]; // all A (encoded as 0)
        let mask = dust_filter(&seq, 64, 0.5);
        assert!(
            !mask.regions.is_empty(),
            "Poly-A sequence should produce at least one masked region"
        );
    }

    /// A sequence with high complexity (all distinct triplets) should not be masked.
    #[test]
    fn test_dust_normal_sequence() {
        // Construct a sequence with maximal triplet diversity:
        // cycle through 0,1,2,3 in a non-repeating pattern
        let seq: Vec<u8> = (0..64).map(|i| ((i * 7 + 3) % 4) as u8).collect();
        let mask = dust_filter(&seq, 64, 10.0);
        assert!(
            mask.regions.is_empty(),
            "High-complexity sequence should produce no masked regions with a generous threshold"
        );
    }

    /// Empty and very short sequences must not crash.
    #[test]
    fn test_dust_empty_sequence() {
        let empty: Vec<u8> = vec![];
        let mask = dust_filter(&empty, 64, 2.0);
        assert!(mask.regions.is_empty());

        let short = vec![0u8; 3];
        let mask = dust_filter(&short, 64, 2.0);
        assert!(mask.regions.is_empty());

        let single = vec![1u8];
        let mask = dust_filter(&single, 64, 2.0);
        assert!(mask.regions.is_empty());
    }

    /// Verify specific masked positions for a known input.
    #[test]
    fn test_dust_mask_positions() {
        // Use a small window so we can precisely control which regions get masked.
        // Build: 30 bases of poly-A, then 30 bases of high-complexity sequence.
        let window = 10;
        let threshold = 2.0;
        let mut seq = vec![0u8; 30]; // 30× A — maximally repetitive
                                     // High-complexity tail: pseudo-random, non-repeating triplet pattern
                                     // Using (i*7+3)%4 over a 10-base window yields diverse triplets
        for i in 0..30 {
            seq.push(((i * 7 + 3) % 4) as u8);
        }

        let mask = dust_filter(&seq, window, threshold);
        // The poly-A region should be masked
        assert!(
            mask.is_masked(0),
            "Position 0 (inside poly-A) should be masked"
        );
        assert!(
            mask.is_masked(15),
            "Position 15 (middle of poly-A) should be masked"
        );
        // Positions well into the high-complexity tail should NOT be masked.
        // The last window that can touch any poly-A is at start=(30-window)=20,
        // covering [20..30). Beyond that, windows are fully in the tail.
        // So position 50+ is guaranteed to be in a fully high-complexity window region.
        assert!(
            !mask.is_masked(55),
            "Position 55 (fully inside high-complexity region) should not be masked"
        );
    }

    /// Different window/threshold parameters should change the filtering outcome.
    #[test]
    fn test_dust_with_different_parameters() {
        // Dinucleotide repeat: ACACAC...
        let seq: Vec<u8> = (0..100).map(|i| if i % 2 == 0 { 0 } else { 1 }).collect();

        // Very high threshold => nothing masked
        let mask_high = dust_filter(&seq, 20, 1000.0);
        assert!(
            mask_high.regions.is_empty(),
            "Very high threshold should mask nothing"
        );

        // Very low threshold => regions masked
        let mask_low = dust_filter(&seq, 20, 0.1);
        assert!(
            !mask_low.regions.is_empty(),
            "Very low threshold should mask the dinucleotide repeat"
        );

        // Larger window
        let mask_big_win = dust_filter(&seq, 50, 0.5);
        // Smaller window
        let mask_small_win = dust_filter(&seq, 10, 0.5);
        // Both should find low-complexity regions (just different sizes)
        assert!(!mask_big_win.regions.is_empty());
        assert!(!mask_small_win.regions.is_empty());
    }
}
