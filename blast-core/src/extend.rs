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
    query: &[u8],       // BLASTNA encoded
    subject: &[u8],     // NCBI2na packed (4 bases/byte)
    q_offset: i32,      // seed position in query
    s_offset: i32,      // seed position in subject
    reward: i32,
    penalty: i32,
    x_dropoff: i32,
) -> Option<UngappedData> {
    // Extend right from the seed
    let q_len = query.len() as i32;
    let s_len_bases = subject.len() as i32 * 4; // approximate

    let mut score = 0i32;
    let mut best_score = 0i32;
    let mut best_len_right = 0i32;
    let mut qi = q_offset;
    let mut si = s_offset;

    while qi < q_len && si < s_len_bases {
        let q_base = query[qi as usize];
        // Decode subject base from packed format
        let s_byte = subject[(si / 4) as usize];
        let s_base = (s_byte >> (6 - 2 * (si % 4))) & 3;

        score += if q_base == s_base { reward } else { penalty };

        if score > best_score {
            best_score = score;
            best_len_right = qi - q_offset + 1;
        }
        if best_score - score > x_dropoff {
            break;
        }
        qi += 1;
        si += 1;
    }

    // Extend left from the seed
    let mut score_left = 0i32;
    let mut best_score_left = 0i32;
    let mut best_len_left = 0i32;
    qi = q_offset - 1;
    si = s_offset - 1;

    while qi >= 0 && si >= 0 {
        let q_base = query[qi as usize];
        let s_byte = subject[(si / 4) as usize];
        let s_base = (s_byte >> (6 - 2 * (si % 4))) & 3;

        score_left += if q_base == s_base { reward } else { penalty };

        if score_left > best_score_left {
            best_score_left = score_left;
            best_len_left = q_offset - qi;
        }
        if best_score_left - score_left > x_dropoff {
            break;
        }
        qi -= 1;
        si -= 1;
    }

    let total_score = best_score + best_score_left;
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
}
