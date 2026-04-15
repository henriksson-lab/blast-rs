//! Rust equivalent of blast_itree.c — interval tree for HSP containment.
//! Used to efficiently check if a new HSP is contained within existing ones.

use std::cmp::{max, min};

/// An interval [start, end) on one axis.
#[derive(Debug, Clone, Copy)]
pub struct Interval {
    pub start: i32,
    pub end: i32,
}

impl Interval {
    pub fn new(start: i32, end: i32) -> Self {
        Interval { start, end }
    }

    pub fn contains(&self, other: &Interval) -> bool {
        self.start <= other.start && self.end >= other.end
    }

    pub fn overlaps(&self, other: &Interval) -> bool {
        self.start < other.end && other.start < self.end
    }

    pub fn length(&self) -> i32 {
        self.end - self.start
    }
}

/// A 2D interval (query range + subject range) representing an HSP.
#[derive(Debug, Clone)]
pub struct Interval2D {
    pub query: Interval,
    pub subject: Interval,
    pub score: i32,
}

/// Simple interval tree for 2D containment checking.
/// For now uses a flat list with linear scan — can be optimized later.
pub struct IntervalTree {
    intervals: Vec<Interval2D>,
}

impl IntervalTree {
    pub fn new(_q_max: i32, _s_max: i32) -> Self {
        IntervalTree {
            intervals: Vec::new(),
        }
    }

    /// Check if a new HSP is contained within any existing interval.
    /// Returns true if the new interval should be rejected (contained).
    pub fn is_contained(&self, query: Interval, subject: Interval) -> bool {
        self.is_contained_with_min_diag_separation(query, subject, 0)
    }

    /// Check if a new HSP is contained within any existing interval, using
    /// NCBI megablast's diagonal-separation rule when nonzero.
    pub fn is_contained_with_min_diag_separation(
        &self,
        query: Interval,
        subject: Interval,
        min_diag_separation: i32,
    ) -> bool {
        for existing in &self.intervals {
            if existing.query.contains(&query)
                && existing.subject.contains(&subject)
                && intervals_are_close_on_diagonal(
                    existing.query,
                    existing.subject,
                    query,
                    subject,
                    min_diag_separation,
                )
            {
                return true;
            }
        }
        false
    }

    /// Check if a new HSP overlaps with any existing interval by more than
    /// a given fraction.
    pub fn has_significant_overlap(
        &self,
        query: Interval,
        subject: Interval,
        max_overlap_fraction: f64,
    ) -> bool {
        self.has_significant_overlap_with_min_diag_separation(
            query,
            subject,
            max_overlap_fraction,
            0,
        )
    }

    pub fn has_significant_overlap_with_min_diag_separation(
        &self,
        query: Interval,
        subject: Interval,
        max_overlap_fraction: f64,
        min_diag_separation: i32,
    ) -> bool {
        for existing in &self.intervals {
            let q_overlap = max(
                0,
                min(existing.query.end, query.end) - max(existing.query.start, query.start),
            );
            let s_overlap = max(
                0,
                min(existing.subject.end, subject.end) - max(existing.subject.start, subject.start),
            );

            let q_frac = if query.length() > 0 {
                q_overlap as f64 / query.length() as f64
            } else {
                0.0
            };
            let s_frac = if subject.length() > 0 {
                s_overlap as f64 / subject.length() as f64
            } else {
                0.0
            };

            if q_frac > max_overlap_fraction
                && s_frac > max_overlap_fraction
                && intervals_are_close_on_diagonal(
                    existing.query,
                    existing.subject,
                    query,
                    subject,
                    min_diag_separation,
                )
            {
                return true;
            }
        }
        false
    }

    /// Add an interval to the tree.
    pub fn insert(&mut self, query: Interval, subject: Interval, score: i32) {
        self.intervals.push(Interval2D {
            query,
            subject,
            score,
        });
    }

    pub fn len(&self) -> usize {
        self.intervals.len()
    }
}

fn intervals_are_close_on_diagonal(
    existing_query: Interval,
    existing_subject: Interval,
    query: Interval,
    subject: Interval,
    min_diag_separation: i32,
) -> bool {
    if min_diag_separation == 0 {
        return true;
    }

    diagonal_distance(
        existing_query.start,
        existing_subject.start,
        query.start,
        subject.start,
    ) < min_diag_separation
        || diagonal_distance(
            existing_query.end,
            existing_subject.end,
            query.end,
            subject.end,
        ) < min_diag_separation
}

fn diagonal_distance(q1: i32, s1: i32, q2: i32, s2: i32) -> i32 {
    ((q1 - s1) - (q2 - s2)).abs()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_containment() {
        let mut tree = IntervalTree::new(100, 100);
        tree.insert(Interval::new(0, 50), Interval::new(0, 50), 100);

        assert!(tree.is_contained(Interval::new(10, 40), Interval::new(10, 40)));
        assert!(!tree.is_contained(Interval::new(0, 60), Interval::new(0, 50)));
    }

    #[test]
    fn test_overlap() {
        let mut tree = IntervalTree::new(100, 100);
        tree.insert(Interval::new(0, 50), Interval::new(0, 50), 100);

        assert!(tree.has_significant_overlap(Interval::new(5, 45), Interval::new(5, 45), 0.5));
        assert!(!tree.has_significant_overlap(Interval::new(45, 55), Interval::new(45, 55), 0.5));
    }

    #[test]
    fn test_megablast_diagonal_separation_keeps_distant_contained_hsp() {
        let mut tree = IntervalTree::new(100, 100);
        tree.insert(Interval::new(3, 40), Interval::new(3, 40), 37);

        assert!(!tree.is_contained_with_min_diag_separation(
            Interval::new(11, 40),
            Interval::new(3, 32),
            6
        ));
        assert!(tree.is_contained_with_min_diag_separation(
            Interval::new(5, 35),
            Interval::new(5, 35),
            6
        ));
    }
}
