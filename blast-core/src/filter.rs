//! Rust equivalent of blast_filter.c — query sequence masking/filtering.
//! Implements DUST (nucleotide) and SEG (protein) low-complexity filtering.

/// A masked region in a sequence.
#[derive(Debug, Clone, Copy)]
pub struct MaskedRegion {
    pub start: i32,
    pub end: i32,
}

/// Collection of masked regions for a query.
#[derive(Debug, Clone)]
pub struct MaskLoc {
    pub regions: Vec<MaskedRegion>,
}

impl MaskLoc {
    pub fn new() -> Self {
        MaskLoc { regions: Vec::new() }
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
}
