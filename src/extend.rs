//! Rust equivalent of blast_extend.c and na_ungapped.c
//! Ungapped and gapped extension structures.

/// Initial hit from word finding (before ungapped extension).
#[derive(Debug, Clone)]
pub struct InitHsp {
    pub query_offset: i32,
    pub subject_offset: i32,
    pub ungapped_data: Option<UngappedData>,
}

/// Data from an ungapped extension.
#[derive(Debug, Clone)]
pub struct UngappedData {
    pub q_start: i32,
    pub s_start: i32,
    pub length: i32,
    pub score: i32,
}

/// List of initial HSPs from the word-finding + ungapped extension phase.
#[derive(Debug)]
pub struct InitHitList {
    pub hits: Vec<InitHsp>,
}

impl InitHitList {
    pub fn new() -> Self {
        InitHitList { hits: Vec::new() }
    }

    pub fn add(&mut self, hsp: InitHsp) {
        self.hits.push(hsp);
    }

    pub fn total(&self) -> usize {
        self.hits.len()
    }

    pub fn reset(&mut self) {
        self.hits.clear();
    }
}

/// Perform ungapped extension on packed nucleotide sequences.
///
/// Extends a seed hit in both directions until the score drops below
/// the x-dropoff threshold. Returns the ungapped alignment data.
pub fn na_ungapped_extend(
    query: &[u8],   // BLASTNA encoded
    subject: &[u8], // NCBI2na packed (4 bases/byte)
    q_offset: i32,  // seed position in query
    s_offset: i32,  // seed position in subject
    reward: i32,
    penalty: i32,
    x_dropoff: i32,
) -> Option<UngappedData> {
    na_ungapped_extend_len(
        query,
        subject,
        subject.len() * 4,
        q_offset,
        s_offset,
        reward,
        penalty,
        x_dropoff,
    )
}

/// Ungapped extension with explicit subject length in bases.
pub fn na_ungapped_extend_len(
    query: &[u8],       // BLASTNA encoded
    subject: &[u8],     // NCBI2na packed (4 bases/byte)
    subject_len: usize, // actual number of bases
    q_offset: i32,      // seed position in query
    s_offset: i32,      // seed position in subject
    reward: i32,
    penalty: i32,
    x_dropoff: i32,
) -> Option<UngappedData> {
    // Extend right from the seed
    let q_len = query.len() as i32;
    let s_len_bases = subject_len as i32;

    let mut score = 0i32;
    let mut sum = 0i32;
    let mut best_len_right = 0i32;
    let mut qi = q_offset;
    let mut si = s_offset;
    let x_dropoff_neg = -x_dropoff;
    let mut x_current = x_dropoff_neg;

    while qi < q_len && si < s_len_bases {
        let q_base = query[qi as usize];
        // Decode subject base from packed format
        let s_byte = subject[(si / 4) as usize];
        let s_base = (s_byte >> (6 - 2 * (si % 4))) & 3;

        sum += blastna_score(q_base, s_base, reward, penalty);

        if sum > 0 {
            score += sum;
            best_len_right = qi - q_offset + 1;
            x_current = (-score).max(x_dropoff_neg);
            sum = 0;
        } else if sum < x_current {
            break;
        }
        qi += 1;
        si += 1;
    }

    // Extend left from the seed
    let mut score_left = 0i32;
    let mut sum_left = 0i32;
    let mut best_len_left = 0i32;
    qi = q_offset - 1;
    si = s_offset - 1;

    while qi >= 0 && si >= 0 {
        let q_base = query[qi as usize];
        let s_byte = subject[(si / 4) as usize];
        let s_base = (s_byte >> (6 - 2 * (si % 4))) & 3;

        sum_left += blastna_score(q_base, s_base, reward, penalty);

        if sum_left > 0 {
            score_left += sum_left;
            best_len_left = q_offset - qi;
            sum_left = 0;
        } else if sum_left < x_dropoff_neg {
            break;
        }
        qi -= 1;
        si -= 1;
    }

    let total_score = score + score_left;
    if total_score <= 0 {
        return None;
    }

    Some(UngappedData {
        q_start: q_offset - best_len_left,
        s_start: s_offset - best_len_left,
        length: best_len_left + best_len_right,
        score: total_score,
    })
}

#[inline(always)]
fn blastna_score(a: u8, b: u8, reward: i32, penalty: i32) -> i32 {
    const MASKS: [u8; 16] = [1, 2, 4, 8, 5, 10, 3, 12, 9, 6, 14, 13, 11, 7, 15, 0];
    const DEGEN: [i32; 16] = [1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 0];

    let a = a as usize;
    let b = b as usize;
    if a >= 15 || b >= 15 {
        return i32::MIN / 2;
    }
    if MASKS[a] & MASKS[b] == 0 {
        return penalty;
    }
    let degen = DEGEN[a.max(b)];
    (((degen - 1) * penalty + reward) as f64 / degen as f64).round() as i32
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_init_hitlist() {
        let mut list = InitHitList::new();
        list.add(InitHsp {
            query_offset: 0,
            subject_offset: 0,
            ungapped_data: None,
        });
        assert_eq!(list.total(), 1);
        list.reset();
        assert_eq!(list.total(), 0);
    }

    #[test]
    fn test_ungapped_extend_perfect_match() {
        // Query: ACGT (BLASTNA: 0,1,2,3)
        let query = vec![0u8, 1, 2, 3];
        // Subject: ACGT packed = 0b00_01_10_11 = 0x1B
        let subject = vec![0x1Bu8];

        let result = na_ungapped_extend(&query, &subject, 0, 0, 2, -3, 20);
        assert!(result.is_some());
        let data = result.unwrap();
        assert_eq!(data.score, 8); // 4 matches * reward=2
        assert_eq!(data.length, 4);
    }

    /// Helper: pack a slice of bases (0..3) into NCBI2na packed bytes.
    fn pack_ncbi2na(bases: &[u8]) -> Vec<u8> {
        let mut packed = vec![0u8; (bases.len() + 3) / 4];
        for (i, &b) in bases.iter().enumerate() {
            packed[i / 4] |= (b & 3) << (6 - 2 * (i % 4));
        }
        packed
    }

    #[test]
    fn test_ungapped_extend_with_xdrop() {
        // Build a sequence that starts matching then diverges.
        // Query:   A C G T  A A A A  (BLASTNA: 0,1,2,3,0,0,0,0)
        // Subject: A C G T  C C C C  (NCBI2na: 0,1,2,3,1,1,1,1)
        // First 4 positions match (+2 each = +8), then 4 mismatches (-3 each).
        // With x_dropoff=5 the extension should stop before consuming all mismatches.
        let query = vec![0u8, 1, 2, 3, 0, 0, 0, 0];
        let subject = pack_ncbi2na(&[0, 1, 2, 3, 1, 1, 1, 1]);

        let result = na_ungapped_extend(&query, &subject, 0, 0, 2, -3, 5);
        assert!(result.is_some());
        let data = result.unwrap();
        // Best score should be from the first 4 matching positions.
        assert_eq!(data.score, 8);
        // The alignment length should be 4 (only the matching region).
        assert_eq!(data.length, 4);
    }

    #[test]
    fn test_ungapped_extend_right_stops_when_total_score_goes_negative() {
        // NCBI's exact nucleotide extender uses a dynamic right-side x-drop:
        // after a +2 first base, a -3 mismatch makes the running total negative
        // and terminates the right extension even when x_dropoff is much larger.
        let query = vec![0u8, 0, 0, 0];
        let subject = pack_ncbi2na(&[0, 1, 0, 0]);

        let result = na_ungapped_extend_len(&query, &subject, 4, 0, 0, 2, -3, 20);
        assert!(result.is_some());
        let data = result.unwrap();
        assert_eq!(data.score, 2);
        assert_eq!(data.length, 1);
    }

    #[test]
    fn test_ungapped_extend_uses_blastna_ambiguity_score() {
        let query = vec![0u8, 0, 0, 0, 14, 0, 0, 0, 0, 0, 0, 0];
        let subject = pack_ncbi2na(&[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]);

        let result = na_ungapped_extend_len(&query, &subject, 12, 0, 0, 1, -5, 20);
        assert!(result.is_some());
        let data = result.unwrap();
        assert_eq!(data.score, 7);
        assert_eq!(data.length, 12);
    }

    #[test]
    fn test_ungapped_extend_at_boundary() {
        // Seed at position 0 — left extension has nothing to extend into.
        // Query and subject: ACGT (4 bases, all match)
        let query = vec![0u8, 1, 2, 3];
        let subject = pack_ncbi2na(&[0, 1, 2, 3]);

        let result = na_ungapped_extend(&query, &subject, 0, 0, 2, -3, 20);
        assert!(result.is_some());
        let data = result.unwrap();
        assert_eq!(data.q_start, 0);
        assert_eq!(data.s_start, 0);
        assert_eq!(data.score, 8);
        assert_eq!(data.length, 4);

        // Seed at the last position — right extension has nothing beyond it.
        let result2 = na_ungapped_extend(&query, &subject, 3, 3, 2, -3, 20);
        assert!(result2.is_some());
        let data2 = result2.unwrap();
        assert_eq!(data2.score, 8);
        assert_eq!(data2.length, 4);
        assert_eq!(data2.q_start, 0);
    }

    #[test]
    fn test_ungapped_extend_all_matches() {
        // 16 bases, all A, perfect match — extension should cover full length.
        let query = vec![0u8; 16];
        let subject = pack_ncbi2na(&vec![0u8; 16]);

        let result = na_ungapped_extend_len(&query, &subject, 16, 8, 8, 2, -3, 100);
        assert!(result.is_some());
        let data = result.unwrap();
        assert_eq!(data.score, 32); // 16 * 2
        assert_eq!(data.length, 16);
        assert_eq!(data.q_start, 0);
        assert_eq!(data.s_start, 0);
    }

    #[test]
    fn test_ungapped_extend_score_calculation() {
        // Manually verify score for a known pattern.
        // Query:   A C A C  G  (BLASTNA: 0,1,0,1,2)
        // Subject: A C G C  G  (NCBI2na: 0,1,2,1,2)
        // Position: 0:match(+2) 1:match(+2) 2:mismatch(-3) 3:match(+2) 4:match(+2)
        // Total = 2+2-3+2+2 = 5
        let query = vec![0u8, 1, 0, 1, 2];
        let subject = pack_ncbi2na(&[0, 1, 2, 1, 2]);

        let result = na_ungapped_extend_len(&query, &subject, 5, 0, 0, 2, -3, 20);
        assert!(result.is_some());
        let data = result.unwrap();
        assert_eq!(data.score, 5);
        assert_eq!(data.length, 5);
    }
}
