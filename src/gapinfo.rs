//! Rust equivalent of gapinfo.c — gap edit script structures.

/// Type of gap alignment operation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u32)]
pub enum GapAlignOpType {
    Del = 0, // Deletion in query (gap in query)
    Del2 = 1,
    Del1 = 2,
    Sub = 3,  // Substitution (aligned pair)
    Ins1 = 4, // Insertion in query (gap in subject)
    Ins2 = 5,
    Ins = 6,
    Decline = 7,
}

/// A gap edit script representing a gapped alignment.
#[derive(Debug, Clone, Default)]
pub struct GapEditScript {
    pub ops: Vec<(GapAlignOpType, i32)>, // (operation, count) pairs
}

impl GapEditScript {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn with_capacity(size: usize) -> Self {
        GapEditScript {
            ops: Vec::with_capacity(size),
        }
    }

    /// Append an edit op, merging with the previous op if the type matches.
    /// Mirrors NCBI's `Blast_PrelimEditBlockToGapEditScript`
    /// (`blast_gapalign.c:2482`) which collapses consecutive same-type ops
    /// when concatenating the reverse + forward halves into a single
    /// edit script. Without this merge, our left-half-end and right-half-
    /// start SUB runs become two ops (`3:6 3:14`) instead of NCBI's
    /// single op (`3:20`), and downstream reevaluation walks them
    /// differently — affecting Kadane's-best-subarray decisions.
    pub fn push(&mut self, op: GapAlignOpType, count: i32) {
        if let Some(last) = self.ops.last_mut() {
            if last.0 == op {
                last.1 += count;
                return;
            }
        }
        self.ops.push((op, count));
    }

    /// Total alignment length (sum of all op counts).
    pub fn alignment_length(&self) -> i32 {
        self.ops.iter().map(|(_, n)| *n).sum()
    }

    /// Render aligned query and subject strings from edit script.
    /// `query` and `subject` are the byte slices covering the aligned region.
    /// `to_char` converts a single encoded byte to a display character.
    pub fn render_alignment(
        &self,
        query: &[u8],
        subject: &[u8],
        to_char: fn(u8) -> char,
    ) -> (String, String) {
        let mut q_str = String::new();
        let mut s_str = String::new();
        let mut q_pos = 0usize;
        let mut s_pos = 0usize;

        for &(op, count) in &self.ops {
            let count = count as usize;
            match op {
                GapAlignOpType::Sub => {
                    for _ in 0..count {
                        q_str.push(if q_pos < query.len() {
                            to_char(query[q_pos])
                        } else {
                            'N'
                        });
                        s_str.push(if s_pos < subject.len() {
                            to_char(subject[s_pos])
                        } else {
                            'N'
                        });
                        q_pos += 1;
                        s_pos += 1;
                    }
                }
                GapAlignOpType::Del | GapAlignOpType::Del1 | GapAlignOpType::Del2 => {
                    // Gap in query, bases in subject
                    for _ in 0..count {
                        q_str.push('-');
                        s_str.push(if s_pos < subject.len() {
                            to_char(subject[s_pos])
                        } else {
                            'N'
                        });
                        s_pos += 1;
                    }
                }
                GapAlignOpType::Ins | GapAlignOpType::Ins1 | GapAlignOpType::Ins2 => {
                    // Bases in query, gap in subject
                    for _ in 0..count {
                        q_str.push(if q_pos < query.len() {
                            to_char(query[q_pos])
                        } else {
                            'N'
                        });
                        s_str.push('-');
                        q_pos += 1;
                    }
                }
                _ => {
                    // Decline or unknown — treat as substitution
                    for _ in 0..count {
                        q_str.push(if q_pos < query.len() {
                            to_char(query[q_pos])
                        } else {
                            'N'
                        });
                        s_str.push(if s_pos < subject.len() {
                            to_char(subject[s_pos])
                        } else {
                            'N'
                        });
                        q_pos += 1;
                        s_pos += 1;
                    }
                }
            }
        }

        (q_str, s_str)
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
                        if q_pos < query.len()
                            && s_pos < subject.len()
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

    /// Helper: BLASTNA byte to IUPAC char (A=0,C=1,G=2,T=3).
    fn to_iupac(b: u8) -> char {
        match b {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            3 => 'T',
            _ => 'N',
        }
    }

    #[test]
    fn test_render_alignment_perfect_match() {
        let mut esp = GapEditScript::new();
        esp.push(GapAlignOpType::Sub, 5);
        let query = &[0u8, 1, 2, 3, 0]; // ACGTA
        let subject = &[0u8, 1, 2, 3, 0]; // ACGTA
        let (qseq, sseq) = esp.render_alignment(query, subject, to_iupac);
        assert_eq!(qseq, "ACGTA");
        assert_eq!(sseq, "ACGTA");
    }

    #[test]
    fn test_render_alignment_with_mismatch() {
        let mut esp = GapEditScript::new();
        esp.push(GapAlignOpType::Sub, 5);
        let query = &[0u8, 1, 2, 3, 0]; // ACGTA
        let subject = &[0u8, 3, 2, 1, 0]; // ATGCA
        let (qseq, sseq) = esp.render_alignment(query, subject, to_iupac);
        assert_eq!(qseq, "ACGTA");
        assert_eq!(sseq, "ATGCA");
    }

    #[test]
    fn test_render_alignment_with_gap_in_query() {
        // Del = gap in query
        let mut esp = GapEditScript::new();
        esp.push(GapAlignOpType::Sub, 3);
        esp.push(GapAlignOpType::Del, 2);
        esp.push(GapAlignOpType::Sub, 2);
        let query = &[0u8, 1, 2, 3, 0]; // ACGTA (5 bases, 2 gapped)
        let subject = &[0u8, 1, 2, 1, 1, 3, 0]; // ACGCCTA (7 bases)
        let (qseq, sseq) = esp.render_alignment(query, subject, to_iupac);
        assert_eq!(qseq, "ACG--TA");
        assert_eq!(sseq, "ACGCCTA");
    }

    #[test]
    fn test_render_alignment_with_gap_in_subject() {
        // Ins = gap in subject
        let mut esp = GapEditScript::new();
        esp.push(GapAlignOpType::Sub, 2);
        esp.push(GapAlignOpType::Ins, 2);
        esp.push(GapAlignOpType::Sub, 2);
        let query = &[0u8, 1, 2, 3, 0, 1]; // ACGTAC (6 bases)
        let subject = &[0u8, 1, 0, 1]; // ACAC (4 bases)
        let (qseq, sseq) = esp.render_alignment(query, subject, to_iupac);
        assert_eq!(qseq, "ACGTAC");
        assert_eq!(sseq, "AC--AC");
    }
}
