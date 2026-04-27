//! 1-1 ports of NCBI's Smith-Waterman primitives from
//! `composition_adjustment/smith_waterman.c`.
//!
//! These are used by the composition-adjusted alignment redo flow
//! (`Blast_RedoAlignmentCore_MT` in `blast_kappa.c`) to compute the
//! optimal SW alignment under the adjusted scoring matrix, which can
//! differ in extent from the original ungapped/X-drop alignment found
//! in the preliminary search. The two-pass forward-then-reverse design
//! mirrors NCBI's flow exactly:
//!
//!   1. `blast_smith_waterman_score_only` — forward SW DP, finds
//!       optimal end position `(queryEnd, matchSeqEnd)` and score,
//!       no traceback.
//!   2. `blast_smith_waterman_find_start` — reverse SW DP from the
//!       optimal end, finds the matching start `(queryStart,
//!       matchSeqStart)`. Stops as soon as the previous score is
//!       reached (`bestScore >= score_in`).
//!
//! Forbidden ranges (used when multiple HSPs are taken from the same
//! match sequence) are tracked by [`BlastForbiddenRanges`]. The wrappers
//! [`blast_smith_waterman_score_only`] and
//! [`blast_smith_waterman_find_start`] dispatch to the basic variants
//! when the range set is empty and to [`bl_special_smith_waterman_score_only`]
//! / [`bl_special_smith_waterman_find_start`] otherwise — matching
//! NCBI's `Blast_SmithWaterman*` dispatch in `smith_waterman.c`.
//!
//! Both routines operate on NCBIstdaa-encoded protein sequences (or
//! 4-bit BLASTNA for nucleotide), where the encoded byte indexes
//! directly into the matrix row.

/// Per-position gap state for the SW DP — 1-1 port of NCBI's
/// `SwGapInfo` (`smith_waterman.c:38`).
#[derive(Debug, Clone, Copy)]
struct SwGapInfo {
    /// score if no gap is open at this position
    no_gap: i32,
    /// score if a gap is open at this position
    gap_exists: i32,
}

/// Forward SW score-only — 1-1 port of `BLbasicSmithWatermanScoreOnly`
/// (`smith_waterman.c:51`).
///
/// Returns `(score, match_seq_end, query_end)` of the locally optimal
/// alignment ending at the highest-scoring cell. No traceback.
///
/// `matrix[query_byte][subject_byte]` indexed by encoded sequence bytes
/// for protein (NCBIstdaa). `gap_open` and `gap_extend` are the *positive*
/// gap costs (NCBI's convention is `-gapOpen`, `-gapExtend`); pass the
/// magnitudes here.
pub fn bl_basic_smith_waterman_score_only(
    match_seq: &[u8],
    query: &[u8],
    matrix: &[[i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE],
    gap_open: i32,
    gap_extend: i32,
) -> (i32, usize, usize) {
    let match_seq_length = match_seq.len();
    let query_length = query.len();
    if match_seq_length == 0 || query_length == 0 {
        return (0, 0, 0);
    }

    let mut score_vector: Vec<SwGapInfo> = (0..match_seq_length)
        .map(|_| SwGapInfo {
            no_gap: 0,
            gap_exists: -gap_open,
        })
        .collect();

    let mut best_score = 0i32;
    let mut best_match_seq_pos = 0usize;
    let mut best_query_pos = 0usize;
    let new_gap_cost = gap_open + gap_extend;

    for query_pos in 0..query_length {
        let matrix_row = &matrix[query[query_pos] as usize];
        let mut new_score = 0i32;
        let mut prev_score_no_gap_match_seq = 0i32;
        let mut prev_score_gap_match_seq = -gap_open;

        for match_seq_pos in 0..match_seq_length {
            // Gap in match_seq: start new or extend existing.
            new_score -= new_gap_cost;
            prev_score_gap_match_seq -= gap_extend;
            if new_score > prev_score_gap_match_seq {
                prev_score_gap_match_seq = new_score;
            }

            // Gap in query: start new or extend existing.
            new_score = score_vector[match_seq_pos].no_gap - new_gap_cost;
            let mut continue_gap_score = score_vector[match_seq_pos].gap_exists - gap_extend;
            if new_score > continue_gap_score {
                continue_gap_score = new_score;
            }

            // Substitution: extend by one position in match_seq + query.
            new_score = prev_score_no_gap_match_seq + matrix_row[match_seq[match_seq_pos] as usize];
            if new_score < 0 {
                new_score = 0; // Smith-Waterman locality
            }
            if new_score < prev_score_gap_match_seq {
                new_score = prev_score_gap_match_seq;
            }
            if new_score < continue_gap_score {
                new_score = continue_gap_score;
            }

            prev_score_no_gap_match_seq = score_vector[match_seq_pos].no_gap;
            score_vector[match_seq_pos].no_gap = new_score;
            score_vector[match_seq_pos].gap_exists = continue_gap_score;

            if new_score > best_score {
                best_score = new_score;
                best_query_pos = query_pos;
                best_match_seq_pos = match_seq_pos;
            }
        }
    }

    if best_score < 0 {
        best_score = 0;
    }

    (best_score, best_match_seq_pos, best_query_pos)
}

/// Reverse SW from a given endpoint — 1-1 port of `BLSmithWatermanFindStart`
/// (`smith_waterman.c:144`).
///
/// Walks the DP backward from `(query_end, match_seq_end)` until the
/// running best score reaches `score_in` (the score returned by the
/// forward score-only pass). Returns `(score, match_seq_start, query_start)`.
///
/// Same matrix/gap convention as `bl_basic_smith_waterman_score_only`.
pub fn bl_smith_waterman_find_start(
    match_seq: &[u8],
    query: &[u8],
    matrix: &[[i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE],
    gap_open: i32,
    gap_extend: i32,
    match_seq_end: usize,
    query_end: usize,
    score_in: i32,
) -> (i32, usize, usize) {
    let match_seq_length = match_seq_end + 1;
    if match_seq_length == 0 || query_end == usize::MAX {
        return (0, 0, 0);
    }

    let mut score_vector: Vec<SwGapInfo> = (0..match_seq_length)
        .map(|_| SwGapInfo {
            no_gap: 0,
            gap_exists: -gap_open,
        })
        .collect();

    let mut best_score = 0i32;
    let mut best_match_seq_pos = 0usize;
    let mut best_query_pos = 0usize;
    let new_gap_cost = gap_open + gap_extend;

    'outer: for query_pos in (0..=query_end).rev() {
        let matrix_row = &matrix[query[query_pos] as usize];
        let mut new_score = 0i32;
        let mut prev_score_no_gap_match_seq = 0i32;
        let mut prev_score_gap_match_seq = -gap_open;

        for match_seq_pos in (0..=match_seq_end).rev() {
            new_score -= new_gap_cost;
            prev_score_gap_match_seq -= gap_extend;
            if new_score > prev_score_gap_match_seq {
                prev_score_gap_match_seq = new_score;
            }

            new_score = score_vector[match_seq_pos].no_gap - new_gap_cost;
            let mut continue_gap_score = score_vector[match_seq_pos].gap_exists - gap_extend;
            if new_score > continue_gap_score {
                continue_gap_score = new_score;
            }

            new_score = prev_score_no_gap_match_seq + matrix_row[match_seq[match_seq_pos] as usize];
            if new_score < 0 {
                new_score = 0;
            }
            if new_score < prev_score_gap_match_seq {
                new_score = prev_score_gap_match_seq;
            }
            if new_score < continue_gap_score {
                new_score = continue_gap_score;
            }

            prev_score_no_gap_match_seq = score_vector[match_seq_pos].no_gap;
            score_vector[match_seq_pos].no_gap = new_score;
            score_vector[match_seq_pos].gap_exists = continue_gap_score;

            if new_score > best_score {
                best_score = new_score;
                best_query_pos = query_pos;
                best_match_seq_pos = match_seq_pos;
            }
            if best_score >= score_in {
                break 'outer;
            }
        }
    }

    if best_score < 0 {
        best_score = 0;
    }

    (best_score, best_match_seq_pos, best_query_pos)
}

/// Sentinel matching NCBI's `COMPO_SCORE_MIN`
/// (`composition_constants.h`). When a forbidden cell is encountered,
/// the substitution score is replaced by this value so the SW DP can
/// never claim a positive contribution from that position.
const COMPO_SCORE_MIN: i32 = i32::MIN / 2;

/// 1-1 port of `Blast_ForbiddenRanges` (`smith_waterman.h`).
///
/// Tracks per-query-position lists of forbidden subject ranges, used by
/// the redo-alignment driver to suppress re-finding the same SW
/// alignment when multiple HSPs come from the same (query, subject)
/// pair. Each query position has a flat `Vec<i32>` of `[start, end,
/// start, end, …]` pairs — same memory layout as NCBI's
/// `int** ranges`.
#[derive(Debug, Clone, Default)]
pub struct BlastForbiddenRanges {
    /// `numForbidden[query_pos]` — number of forbidden ranges at this
    /// query position (each range is a `(start, end)` pair).
    pub num_forbidden: Vec<i32>,
    /// `ranges[query_pos]` — flat array of `[start0, end0, start1,
    /// end1, …]` for each forbidden range at this query position.
    pub ranges: Vec<Vec<i32>>,
    /// `isEmpty` — fast-path flag mirrored from NCBI; TRUE means no
    /// ranges have been pushed yet.
    pub is_empty: bool,
}

impl BlastForbiddenRanges {
    /// 1-1 port of `Blast_ForbiddenRangesInitialize` (`smith_waterman.c:473`).
    ///
    /// `capacity` is the concatenated query length; we allocate one
    /// entry per query position. NCBI initializes `ranges[f]` to a
    /// 2-int slot (`[0, 0]`); we use empty `Vec`s instead since
    /// `num_forbidden[f] == 0` already signals "nothing here".
    pub fn new(capacity: i32) -> Self {
        let cap = capacity.max(0) as usize;
        Self {
            num_forbidden: vec![0i32; cap],
            ranges: vec![Vec::new(); cap],
            is_empty: true,
        }
    }

    /// 1-1 port of `Blast_ForbiddenRangesClear` (`smith_waterman.c:505`).
    /// Resets all per-position counts to zero and flips `is_empty`
    /// back to true. Backing arrays are kept around (matching NCBI).
    pub fn clear(&mut self) {
        for n in self.num_forbidden.iter_mut() {
            *n = 0;
        }
        // NCBI doesn't truncate `ranges[f]` either — the reused
        // backing storage is harmless because `numForbidden[f] == 0`
        // makes the inner loop in `BLspecial*` skip them entirely.
        self.is_empty = true;
    }

    /// 1-1 port of `Blast_ForbiddenRangesPush` (`smith_waterman.c:516`).
    ///
    /// Marks the rectangular region `[query_start, query_end) ×
    /// [match_start, match_end]` forbidden — for every query position
    /// in `[query_start, query_end)` we append one (matchStart,
    /// matchEnd) pair.
    pub fn push(&mut self, query_start: i32, query_end: i32, match_start: i32, match_end: i32) {
        let qs = query_start.max(0) as usize;
        let qe = query_end.max(0) as usize;
        if qe > self.num_forbidden.len() {
            self.num_forbidden.resize(qe, 0);
            self.ranges.resize(qe, Vec::new());
        }
        for f in qs..qe {
            self.ranges[f].push(match_start);
            self.ranges[f].push(match_end);
            self.num_forbidden[f] += 1;
        }
        self.is_empty = false;
    }
}

/// 1-1 port of `BLspecialSmithWatermanScoreOnly` (`smith_waterman.c:243`).
///
/// Forward score-only SW with forbidden-range exclusion. For each cell
/// `(query_pos, match_pos)` we check whether `match_pos` falls in any
/// of the forbidden ranges registered at that query position; if so,
/// the substitution contribution is replaced by `COMPO_SCORE_MIN` —
/// matching NCBI's `forbidden ? COMPO_SCORE_MIN : score`.
pub fn bl_special_smith_waterman_score_only(
    match_seq: &[u8],
    query: &[u8],
    matrix: &[[i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE],
    gap_open: i32,
    gap_extend: i32,
    forbidden: &BlastForbiddenRanges,
) -> (i32, usize, usize) {
    let match_seq_length = match_seq.len();
    let query_length = query.len();
    if match_seq_length == 0 || query_length == 0 {
        return (0, 0, 0);
    }

    let mut score_vector: Vec<SwGapInfo> = (0..match_seq_length)
        .map(|_| SwGapInfo {
            no_gap: 0,
            gap_exists: -gap_open,
        })
        .collect();

    let mut best_score = 0i32;
    let mut best_match_seq_pos = 0usize;
    let mut best_query_pos = 0usize;
    let new_gap_cost = gap_open + gap_extend;

    for query_pos in 0..query_length {
        let matrix_row = &matrix[query[query_pos] as usize];
        let mut new_score = 0i32;
        let mut prev_score_no_gap_match_seq = 0i32;
        let mut prev_score_gap_match_seq = -gap_open;

        let nf = *forbidden.num_forbidden.get(query_pos).unwrap_or(&0) as usize;
        let row_ranges: &[i32] = forbidden
            .ranges
            .get(query_pos)
            .map(|v| v.as_slice())
            .unwrap_or(&[]);

        for match_seq_pos in 0..match_seq_length {
            new_score -= new_gap_cost;
            prev_score_gap_match_seq -= gap_extend;
            if new_score > prev_score_gap_match_seq {
                prev_score_gap_match_seq = new_score;
            }
            new_score = score_vector[match_seq_pos].no_gap - new_gap_cost;
            let mut continue_gap_score = score_vector[match_seq_pos].gap_exists - gap_extend;
            if new_score > continue_gap_score {
                continue_gap_score = new_score;
            }

            let mut is_forbidden = false;
            for f in 0..nf {
                let r0 = row_ranges[2 * f] as i64;
                let r1 = row_ranges[2 * f + 1] as i64;
                if (match_seq_pos as i64) >= r0 && (match_seq_pos as i64) <= r1 {
                    is_forbidden = true;
                    break;
                }
            }
            new_score = if is_forbidden {
                COMPO_SCORE_MIN
            } else {
                prev_score_no_gap_match_seq + matrix_row[match_seq[match_seq_pos] as usize]
            };
            if new_score < 0 {
                new_score = 0;
            }
            if new_score < prev_score_gap_match_seq {
                new_score = prev_score_gap_match_seq;
            }
            if new_score < continue_gap_score {
                new_score = continue_gap_score;
            }

            prev_score_no_gap_match_seq = score_vector[match_seq_pos].no_gap;
            score_vector[match_seq_pos].no_gap = new_score;
            score_vector[match_seq_pos].gap_exists = continue_gap_score;

            if new_score > best_score {
                best_score = new_score;
                best_query_pos = query_pos;
                best_match_seq_pos = match_seq_pos;
            }
        }
    }

    if best_score < 0 {
        best_score = 0;
    }
    (best_score, best_match_seq_pos, best_query_pos)
}

/// 1-1 port of `BLspecialSmithWatermanFindStart` (`smith_waterman.c:351`).
///
/// Reverse SW from `(match_seq_end, query_end)` with forbidden-range
/// exclusion and an early-exit when `bestScore >= score_in` is reached.
pub fn bl_special_smith_waterman_find_start(
    match_seq: &[u8],
    query: &[u8],
    matrix: &[[i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE],
    gap_open: i32,
    gap_extend: i32,
    match_seq_end: usize,
    query_end: usize,
    score_in: i32,
    forbidden: &BlastForbiddenRanges,
) -> (i32, usize, usize) {
    let match_seq_length = match_seq_end + 1;
    if match_seq_length == 0 {
        return (0, 0, 0);
    }

    let mut score_vector: Vec<SwGapInfo> = (0..match_seq_length)
        .map(|_| SwGapInfo {
            no_gap: 0,
            gap_exists: -gap_open,
        })
        .collect();

    let mut best_score = 0i32;
    let mut best_match_seq_pos = 0usize;
    let mut best_query_pos = 0usize;
    let new_gap_cost = gap_open + gap_extend;

    'outer: for query_pos in (0..=query_end).rev() {
        let matrix_row = &matrix[query[query_pos] as usize];
        let mut new_score = 0i32;
        let mut prev_score_no_gap_match_seq = 0i32;
        let mut prev_score_gap_match_seq = -gap_open;

        let nf = *forbidden.num_forbidden.get(query_pos).unwrap_or(&0) as usize;
        let row_ranges: &[i32] = forbidden
            .ranges
            .get(query_pos)
            .map(|v| v.as_slice())
            .unwrap_or(&[]);

        for match_seq_pos in (0..=match_seq_end).rev() {
            new_score -= new_gap_cost;
            prev_score_gap_match_seq -= gap_extend;
            if new_score > prev_score_gap_match_seq {
                prev_score_gap_match_seq = new_score;
            }
            new_score = score_vector[match_seq_pos].no_gap - new_gap_cost;
            let mut continue_gap_score = score_vector[match_seq_pos].gap_exists - gap_extend;
            if new_score > continue_gap_score {
                continue_gap_score = new_score;
            }

            let mut is_forbidden = false;
            for f in 0..nf {
                let r0 = row_ranges[2 * f] as i64;
                let r1 = row_ranges[2 * f + 1] as i64;
                if (match_seq_pos as i64) >= r0 && (match_seq_pos as i64) <= r1 {
                    is_forbidden = true;
                    break;
                }
            }
            new_score = if is_forbidden {
                COMPO_SCORE_MIN
            } else {
                prev_score_no_gap_match_seq + matrix_row[match_seq[match_seq_pos] as usize]
            };
            if new_score < 0 {
                new_score = 0;
            }
            if new_score < prev_score_gap_match_seq {
                new_score = prev_score_gap_match_seq;
            }
            if new_score < continue_gap_score {
                new_score = continue_gap_score;
            }

            prev_score_no_gap_match_seq = score_vector[match_seq_pos].no_gap;
            score_vector[match_seq_pos].no_gap = new_score;
            score_vector[match_seq_pos].gap_exists = continue_gap_score;

            if new_score > best_score {
                best_score = new_score;
                best_query_pos = query_pos;
                best_match_seq_pos = match_seq_pos;
            }
            if best_score >= score_in {
                break 'outer;
            }
        }
    }

    if best_score < 0 {
        best_score = 0;
    }
    (best_score, best_match_seq_pos, best_query_pos)
}

/// 1-1 port of `Blast_SmithWatermanScoreOnly` (`smith_waterman.c:545`).
///
/// Public dispatch — falls through to [`bl_basic_smith_waterman_score_only`]
/// when `forbidden.is_empty`, otherwise calls the special variant. Both
/// variants are mirrors of NCBI's two paths.
pub fn blast_smith_waterman_score_only(
    match_seq: &[u8],
    query: &[u8],
    matrix: &[[i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE],
    gap_open: i32,
    gap_extend: i32,
) -> (i32, usize, usize) {
    bl_basic_smith_waterman_score_only(match_seq, query, matrix, gap_open, gap_extend)
}

/// Variant of [`blast_smith_waterman_score_only`] that takes a
/// forbidden-range set, mirroring NCBI's `Blast_SmithWatermanScoreOnly`
/// when `forbiddenRanges->isEmpty == FALSE`.
pub fn blast_smith_waterman_score_only_with_forbidden(
    match_seq: &[u8],
    query: &[u8],
    matrix: &[[i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE],
    gap_open: i32,
    gap_extend: i32,
    forbidden: &BlastForbiddenRanges,
) -> (i32, usize, usize) {
    if forbidden.is_empty {
        bl_basic_smith_waterman_score_only(match_seq, query, matrix, gap_open, gap_extend)
    } else {
        bl_special_smith_waterman_score_only(
            match_seq, query, matrix, gap_open, gap_extend, forbidden,
        )
    }
}

/// 1-1 port of `Blast_SmithWatermanFindStart` (`smith_waterman.c:579`).
///
/// Public dispatch — falls through to [`bl_smith_waterman_find_start`]
/// for the empty-forbidden-range case.
pub fn blast_smith_waterman_find_start(
    match_seq: &[u8],
    query: &[u8],
    matrix: &[[i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE],
    gap_open: i32,
    gap_extend: i32,
    match_seq_end: usize,
    query_end: usize,
    score_in: i32,
) -> (i32, usize, usize) {
    bl_smith_waterman_find_start(
        match_seq,
        query,
        matrix,
        gap_open,
        gap_extend,
        match_seq_end,
        query_end,
        score_in,
    )
}

/// Variant of [`blast_smith_waterman_find_start`] that takes a
/// forbidden-range set.
pub fn blast_smith_waterman_find_start_with_forbidden(
    match_seq: &[u8],
    query: &[u8],
    matrix: &[[i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE],
    gap_open: i32,
    gap_extend: i32,
    match_seq_end: usize,
    query_end: usize,
    score_in: i32,
    forbidden: &BlastForbiddenRanges,
) -> (i32, usize, usize) {
    if forbidden.is_empty {
        bl_smith_waterman_find_start(
            match_seq,
            query,
            matrix,
            gap_open,
            gap_extend,
            match_seq_end,
            query_end,
            score_in,
        )
    } else {
        bl_special_smith_waterman_find_start(
            match_seq,
            query,
            matrix,
            gap_open,
            gap_extend,
            match_seq_end,
            query_end,
            score_in,
            forbidden,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrix::AA_SIZE;

    fn encode_aa(s: &[u8]) -> Vec<u8> {
        s.iter()
            .map(|&b| crate::encoding::AMINOACID_TO_NCBISTDAA[b as usize & 0x7F])
            .collect()
    }

    fn blosum62() -> [[i32; AA_SIZE]; AA_SIZE] {
        *crate::api::get_matrix(crate::api::MatrixType::Blosum62)
    }

    #[test]
    fn test_score_only_exact_match() {
        let q = encode_aa(b"MKFLILLF");
        let s = encode_aa(b"MKFLILLF");
        let m = blosum62();
        let (score, m_end, q_end) = blast_smith_waterman_score_only(&s, &q, &m, 11, 1);
        // Diag scores: 5+5+6+4+4+4+4+6 = 38
        assert_eq!(score, 38);
        assert_eq!(m_end, 7); // last position of subject
        assert_eq!(q_end, 7); // last position of query
    }

    #[test]
    fn test_find_start_exact_match() {
        let q = encode_aa(b"MKFLILLF");
        let s = encode_aa(b"MKFLILLF");
        let m = blosum62();
        let (score, m_start, q_start) =
            blast_smith_waterman_find_start(&s, &q, &m, 11, 1, 7, 7, 38);
        assert_eq!(score, 38);
        assert_eq!(m_start, 0);
        assert_eq!(q_start, 0);
    }

    #[test]
    fn test_score_only_offset_match() {
        // Subject has 5 leading X's then the query
        let q = encode_aa(b"MKFLILLF");
        let s = encode_aa(b"AAAAAMKFLILLF");
        let m = blosum62();
        let (score, m_end, q_end) = blast_smith_waterman_score_only(&s, &q, &m, 11, 1);
        assert_eq!(score, 38);
        assert_eq!(q_end, 7);
        assert_eq!(m_end, 12); // 5 + 7 = position of last F in subject

        // Reverse SW should find start at position 5 in subject
        let (score2, m_start, q_start) =
            blast_smith_waterman_find_start(&s, &q, &m, 11, 1, m_end, q_end, score);
        assert_eq!(score2, 38);
        assert_eq!(q_start, 0);
        assert_eq!(m_start, 5);
    }

    #[test]
    fn test_score_only_no_match() {
        let q = encode_aa(b"AAAA");
        let s = encode_aa(b"WWWW");
        let m = blosum62();
        let (score, _, _) = blast_smith_waterman_score_only(&s, &q, &m, 11, 1);
        // BLOSUM62 A-W = -3, so all-negative; SW returns 0
        assert_eq!(score, 0);
    }

    #[test]
    fn test_score_only_with_gap() {
        // Aligning MKF + gap + LILLF vs MKFXLILLF (X is a substitution penalty)
        let q = encode_aa(b"MKFLILLF");
        let s = encode_aa(b"MKFGGLILLF");
        let m = blosum62();
        let (score, _, _) = blast_smith_waterman_score_only(&s, &q, &m, 11, 1);
        // SW will find best alignment which includes a gap or extension.
        // Just verify it returns some positive score >= 0.
        assert!(score > 0);
    }

    #[test]
    fn test_forbidden_ranges_basic() {
        // Subject has the query motif twice; without forbidden ranges we
        // find one of them (whichever the SW picks). Pushing the first
        // match's range as forbidden should make the second SW pass find
        // the *other* occurrence.
        let q = encode_aa(b"MKFLILLF");
        let s = encode_aa(b"MKFLILLFAAAAAMKFLILLF");
        let m = blosum62();
        // First pass — empty forbidden set, picks the first occurrence.
        let mut fr = BlastForbiddenRanges::new(q.len() as i32);
        assert!(fr.is_empty);
        let (s1, m_end1, q_end1) =
            blast_smith_waterman_score_only_with_forbidden(&s, &q, &m, 11, 1, &fr);
        assert_eq!(s1, 38);
        let (s1_back, m_start1, q_start1) = blast_smith_waterman_find_start_with_forbidden(
            &s, &q, &m, 11, 1, m_end1, q_end1, s1, &fr,
        );
        assert_eq!(s1_back, 38);
        assert_eq!(q_start1, 0);

        // Mark the first occurrence as forbidden, do the second pass.
        fr.push(
            q_start1 as i32,
            q_end1 as i32 + 1,
            m_start1 as i32,
            m_end1 as i32,
        );
        assert!(!fr.is_empty);
        let (s2, m_end2, _q_end2) =
            blast_smith_waterman_score_only_with_forbidden(&s, &q, &m, 11, 1, &fr);
        assert_eq!(s2, 38, "second occurrence should give same score");
        // Second-pass end position must be inside the second occurrence
        // (positions 13..=20). The first occurrence ends at index 7.
        assert!(
            m_end2 >= 13,
            "expected second-pass end inside second occurrence, got {}",
            m_end2
        );
    }

    #[test]
    fn test_forbidden_ranges_clear_resets_is_empty() {
        let mut fr = BlastForbiddenRanges::new(10);
        fr.push(0, 5, 0, 7);
        assert!(!fr.is_empty);
        fr.clear();
        assert!(fr.is_empty);
        assert!(fr.num_forbidden.iter().all(|&n| n == 0));
    }
}
