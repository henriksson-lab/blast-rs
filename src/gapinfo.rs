//! Rust equivalent of gapinfo.c — gap edit script structures.

/// Type of gap alignment operation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u32)]
pub enum GapAlignOpType {
    Del = 0,   // Deletion in query (gap in query)
    Del2 = 1,
    Del1 = 2,
    Sub = 3,   // Substitution (aligned pair)
    Ins1 = 4,  // Insertion in query (gap in subject)
    Ins2 = 5,
    Ins = 6,
    Decline = 7,
}

/// A gap edit script representing a gapped alignment.
#[derive(Debug, Clone)]
pub struct GapEditScript {
    pub ops: Vec<(GapAlignOpType, i32)>, // (operation, count) pairs
}

impl GapEditScript {
    pub fn new() -> Self {
        GapEditScript { ops: Vec::new() }
    }

    pub fn with_capacity(size: usize) -> Self {
        GapEditScript {
            ops: Vec::with_capacity(size),
        }
    }

    pub fn push(&mut self, op: GapAlignOpType, count: i32) {
        self.ops.push((op, count));
    }

    /// Total alignment length (sum of all op counts).
    pub fn alignment_length(&self) -> i32 {
        self.ops.iter().map(|(_, n)| *n).sum()
    }

    /// Count identities given query and subject byte slices.
    pub fn count_identities(&self, query: &[u8], subject: &[u8]) -> (i32, i32, i32) {
        let mut q_pos = 0usize;
        let mut s_pos = 0usize;
        let mut align_len = 0i32;
        let mut num_ident = 0i32;
        let mut gap_opens = 0i32;

        for &(op, count) in &self.ops {
            let count = count as usize;
            align_len += count as i32;

            match op {
                GapAlignOpType::Sub => {
                    for _ in 0..count {
                        if q_pos < query.len() && s_pos < subject.len()
                            && query[q_pos] == subject[s_pos]
                        {
                            num_ident += 1;
                        }
                        q_pos += 1;
                        s_pos += 1;
                    }
                }
                GapAlignOpType::Del | GapAlignOpType::Del1 | GapAlignOpType::Del2 => {
                    s_pos += count;
                    gap_opens += 1;
                }
                GapAlignOpType::Ins | GapAlignOpType::Ins1 | GapAlignOpType::Ins2 => {
                    q_pos += count;
                    gap_opens += 1;
                }
                _ => {
                    q_pos += count;
                    s_pos += count;
                }
            }
        }
        (align_len, num_ident, gap_opens)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_perfect_match() {
        let mut esp = GapEditScript::new();
        esp.push(GapAlignOpType::Sub, 10);
        let query = &[0u8, 1, 2, 3, 0, 1, 2, 3, 0, 1];
        let subject = &[0u8, 1, 2, 3, 0, 1, 2, 3, 0, 1];
        let (len, ident, gaps) = esp.count_identities(query, subject);
        assert_eq!(len, 10);
        assert_eq!(ident, 10);
        assert_eq!(gaps, 0);
    }

    #[test]
    fn test_with_gap() {
        let mut esp = GapEditScript::new();
        esp.push(GapAlignOpType::Sub, 5);
        esp.push(GapAlignOpType::Del, 2);
        esp.push(GapAlignOpType::Sub, 3);
        let query = &[0u8, 1, 2, 3, 0, 1, 2, 3];
        let subject = &[0u8, 1, 2, 3, 0, 9, 9, 1, 2, 3];
        let (len, ident, gaps) = esp.count_identities(query, subject);
        assert_eq!(len, 10);
        assert_eq!(ident, 8);
        assert_eq!(gaps, 1);
    }
}
