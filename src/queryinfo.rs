//! Rust equivalent of blast_query_info.c — query information management.

/// NCBI `ESequenceSegment` values used by `BlastContextInfo.segment_flags`.
pub const E_NO_SEGMENTS: i32 = 0;
pub const E_FIRST_SEGMENT: i32 = 1;
pub const E_LAST_SEGMENT: i32 = 2;

/// Information about a single query context (strand/frame).
#[derive(Debug, Clone)]
pub struct ContextInfo {
    pub query_offset: i32,
    pub query_length: i32,
    pub eff_searchsp: i64,
    pub length_adjustment: i32,
    pub query_index: i32,
    pub frame: i32,
    pub is_valid: bool,
    pub segment_flags: i32,
}

/// Information about all queries in a search.
#[derive(Debug, Clone)]
pub struct QueryInfo {
    pub num_queries: i32,
    pub contexts: Vec<ContextInfo>,
    pub max_length: u32,
}

impl QueryInfo {
    /// Create QueryInfo for protein-query programs with one context per query.
    pub fn new_blastp(query_lengths: &[usize]) -> Self {
        let num_queries = query_lengths.len() as i32;
        let mut contexts = Vec::new();
        let mut offset = 0i32;

        for (qi, &qlen) in query_lengths.iter().enumerate() {
            contexts.push(ContextInfo {
                query_offset: offset,
                query_length: qlen as i32,
                eff_searchsp: 0,
                length_adjustment: 0,
                query_index: qi as i32,
                frame: 0,
                is_valid: qlen > 0,
                segment_flags: E_NO_SEGMENTS,
            });
            offset += qlen as i32 + 1;
        }

        let max_length = query_lengths.iter().copied().max().unwrap_or(0) as u32;
        QueryInfo {
            num_queries,
            contexts,
            max_length,
        }
    }

    /// Create QueryInfo for blastn with the given query lengths.
    pub fn new_blastn(query_lengths: &[usize]) -> Self {
        let num_queries = query_lengths.len() as i32;
        let mut contexts = Vec::new();
        let mut offset = 0i32;

        for (qi, &qlen) in query_lengths.iter().enumerate() {
            // Plus strand
            contexts.push(ContextInfo {
                query_offset: offset,
                query_length: qlen as i32,
                eff_searchsp: 0,
                length_adjustment: 0,
                query_index: qi as i32,
                frame: 1,
                is_valid: true,
                segment_flags: E_NO_SEGMENTS,
            });
            offset += qlen as i32 + 1; // +1 for sentinel between strands

            // Minus strand
            contexts.push(ContextInfo {
                query_offset: offset,
                query_length: qlen as i32,
                eff_searchsp: 0,
                length_adjustment: 0,
                query_index: qi as i32,
                frame: -1,
                is_valid: true,
                segment_flags: E_NO_SEGMENTS,
            });
            offset += qlen as i32 + 1;
        }

        let max_length = query_lengths.iter().copied().max().unwrap_or(0) as u32;

        QueryInfo {
            num_queries,
            contexts,
            max_length,
        }
    }

    /// Create QueryInfo for a single translated nucleotide query from the
    /// offsets returned by `BLAST_GetAllTranslations`.
    ///
    /// The six contexts follow NCBI `BLAST_ContextToFrame` order:
    /// `+1,+2,+3,-1,-2,-3`. Context offsets point into the flat translation
    /// buffer after the leading NULLB sentinel for each frame.
    pub fn new_translated_query_from_offsets(
        frame_offsets: &[u32; crate::util::NUM_FRAMES + 1],
    ) -> Self {
        let contexts: Vec<ContextInfo> = (0..crate::util::NUM_FRAMES)
            .map(|ctx| {
                let begin = (frame_offsets[ctx] + 1) as i32;
                let end = frame_offsets[ctx + 1] as i32;
                let query_length = (end - begin).max(0);
                ContextInfo {
                    query_offset: begin,
                    query_length,
                    eff_searchsp: 0,
                    length_adjustment: 0,
                    query_index: 0,
                    frame: crate::util::blast_context_to_frame_blastx(ctx as u32),
                    is_valid: query_length > 0,
                    segment_flags: E_NO_SEGMENTS,
                }
            })
            .collect();
        let max_length = contexts
            .iter()
            .map(|ctx| ctx.query_length.max(0) as u32)
            .max()
            .unwrap_or(0);
        QueryInfo {
            num_queries: 1,
            contexts,
            max_length,
        }
    }

    /// Get context info for a given context index.
    pub fn get_context(&self, context: usize) -> Option<&ContextInfo> {
        self.contexts.get(context)
    }

    /// Get the query index for a given context.
    pub fn query_index(&self, context: usize) -> i32 {
        self.contexts.get(context).map_or(-1, |c| c.query_index)
    }

    pub fn query_segment_flags(&self, query_idx: usize) -> i32 {
        self.contexts
            .get(query_idx * crate::util::NUM_STRANDS)
            .map_or(E_NO_SEGMENTS, |context| context.segment_flags)
    }

    pub fn set_query_segment_flags(&mut self, query_idx: usize, segment_flags: i32) {
        let first_context = query_idx * crate::util::NUM_STRANDS;
        for context in self
            .contexts
            .iter_mut()
            .skip(first_context)
            .take(crate::util::NUM_STRANDS)
        {
            context.segment_flags = segment_flags;
        }
    }

    /// Number of contexts.
    pub fn num_contexts(&self) -> usize {
        self.contexts.len()
    }
}

/// 1-1 port of `Blast_GetQueryIndexFromContext`
/// (`blast_query_info.c:40`). Maps a context index to its parent
/// query index, accounting for per-program context multiplicity:
/// - protein/PSSM queries: `context` directly (1 context per query).
/// - translated queries (blastx/tblastx/RPS-tblastn): `context / 6`.
/// - nucleotide queries (blastn): `context / 2`.
pub fn blast_get_query_index_from_context(
    context: i32,
    program: crate::program::ProgramType,
) -> i32 {
    if crate::program::blast_query_is_protein(program) {
        context
    } else if crate::program::blast_query_is_translated(program) {
        context / crate::util::NUM_FRAMES as i32
    } else {
        context / crate::util::NUM_STRANDS as i32
    }
}

/// 1-1 port of `Blast_GetQueryIndexFromQueryOffset`
/// (`blast_query_info.c:52`).
pub fn blast_get_query_index_from_query_offset(
    query_offset: i32,
    program: crate::program::ProgramType,
    query_info: &QueryInfo,
) -> i32 {
    let context = bsearch_context_info(query_offset, query_info);
    blast_get_query_index_from_context(context, program)
}

/// 1-1 port of `BSearchContextInfo` (`blast_query_info.c:219`).
/// Returns the index of the context with the greatest query offset not
/// exceeding `n`.
pub fn bsearch_context_info(n: i32, query_info: &QueryInfo) -> i32 {
    let size = query_info.contexts.len();
    if size == 0 {
        return -1;
    }

    let mut b = 0usize;
    let mut e = size;

    while b < e - 1 {
        let m = (b + e) / 2;
        if query_info.contexts[m].query_offset > n {
            e = m;
        } else {
            b = m;
        }
    }

    b as i32
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new_blastn_single() {
        let qi = QueryInfo::new_blastn(&[50]);
        assert_eq!(qi.num_queries, 1);
        assert_eq!(qi.num_contexts(), 2);
        assert_eq!(qi.max_length, 50);
        assert_eq!(qi.contexts[0].query_offset, 0);
        assert_eq!(qi.contexts[0].query_length, 50);
        assert_eq!(qi.contexts[0].frame, 1);
        assert_eq!(qi.contexts[1].query_offset, 51);
        assert_eq!(qi.contexts[1].query_length, 50);
        assert_eq!(qi.contexts[1].frame, -1);
    }

    #[test]
    fn test_new_blastn_multi() {
        let qi = QueryInfo::new_blastn(&[100, 50]);
        assert_eq!(qi.num_queries, 2);
        assert_eq!(qi.num_contexts(), 4);
        assert_eq!(qi.max_length, 100);
        assert_eq!(qi.contexts[2].query_index, 1);
    }

    #[test]
    fn test_new_blastp_multi() {
        let qi = QueryInfo::new_blastp(&[150, 0, 40]);
        assert_eq!(qi.num_queries, 3);
        assert_eq!(qi.num_contexts(), 3);
        assert_eq!(qi.max_length, 150);
        assert_eq!(qi.contexts[0].query_offset, 0);
        assert_eq!(qi.contexts[0].frame, 0);
        assert_eq!(qi.contexts[1].query_offset, 151);
        assert!(!qi.contexts[1].is_valid);
        assert_eq!(qi.contexts[2].query_offset, 152);
        assert_eq!(qi.contexts[2].query_index, 2);
    }

    /// Port of NCBI queryinfo_unit_test: blastn with one query has 2 contexts (plus, minus).
    #[test]
    fn test_queryinfo_blastn_two_contexts() {
        let qi = QueryInfo::new_blastn(&[200]);
        assert_eq!(qi.num_queries, 1);
        assert_eq!(qi.num_contexts(), 2);

        // Context 0: plus strand
        let plus = qi.get_context(0).unwrap();
        assert_eq!(plus.frame, 1);
        assert_eq!(plus.query_length, 200);
        assert_eq!(plus.query_index, 0);
        assert!(plus.is_valid);

        // Context 1: minus strand
        let minus = qi.get_context(1).unwrap();
        assert_eq!(minus.frame, -1);
        assert_eq!(minus.query_length, 200);
        assert_eq!(minus.query_index, 0);
        assert!(minus.is_valid);
    }

    /// Port of NCBI queryinfo_unit_test: blastp should have 1 context per query.
    /// Since we only have new_blastn, we simulate blastp as a single-context case.
    #[test]
    fn test_queryinfo_blastp_one_context() {
        let query_len = 150;
        let qi = QueryInfo::new_blastp(&[query_len as usize]);
        assert_eq!(qi.num_queries, 1);
        assert_eq!(qi.num_contexts(), 1);
        assert_eq!(qi.contexts[0].frame, 0);
        assert_eq!(qi.contexts[0].query_length, query_len);
    }

    /// Port of NCBI queryinfo_unit_test: blastx should have 6 contexts (6 reading frames).
    #[test]
    fn test_queryinfo_blastx_six_contexts() {
        let nt = vec![1u8, 8, 4, 4, 2, 8]; // ATGGCT in NCBI4na.
        let (_buf, offsets) = crate::util::blast_get_all_translations_ncbi4na(
            &nt,
            nt.len(),
            &crate::util::STANDARD_GENETIC_CODE,
        );
        let qi = QueryInfo::new_translated_query_from_offsets(&offsets);
        assert_eq!(qi.num_queries, 1);
        assert_eq!(qi.num_contexts(), 6);
        assert_eq!(qi.contexts[0].frame, 1);
        assert_eq!(qi.contexts[0].query_length, 2);
        assert_eq!(qi.contexts[1].query_length, 1);
        assert_eq!(qi.contexts[2].query_length, 1);
        assert_eq!(qi.contexts[3].frame, -1);
        assert_eq!(qi.contexts[3].query_length, 2);
        assert_eq!(qi.contexts[5].frame, -3);
        assert_eq!(qi.max_length, 2);
    }

    #[test]
    fn test_queryinfo_blastx_short_sequence_contexts_are_invalid() {
        let nt = vec![1u8, 2]; // AC in NCBI4na; no complete codons.
        let (_buf, offsets) = crate::util::blast_get_all_translations_ncbi4na(
            &nt,
            nt.len(),
            &crate::util::STANDARD_GENETIC_CODE,
        );
        let qi = QueryInfo::new_translated_query_from_offsets(&offsets);
        assert_eq!(qi.num_contexts(), 6);
        assert_eq!(qi.max_length, 0);
        assert!(qi.contexts.iter().all(|ctx| ctx.query_length == 0));
        assert!(qi.contexts.iter().all(|ctx| !ctx.is_valid));
    }

    /// Port of NCBI queryinfo_unit_test: multiple queries multiply the context count.
    #[test]
    fn test_queryinfo_multi_query() {
        let qi = QueryInfo::new_blastn(&[100, 200, 50]);
        assert_eq!(qi.num_queries, 3);
        // blastn: 2 contexts per query (plus + minus)
        assert_eq!(qi.num_contexts(), 6);
        assert_eq!(qi.max_length, 200);

        // Verify query indices
        assert_eq!(qi.query_index(0), 0); // query 0 plus
        assert_eq!(qi.query_index(1), 0); // query 0 minus
        assert_eq!(qi.query_index(2), 1); // query 1 plus
        assert_eq!(qi.query_index(3), 1); // query 1 minus
        assert_eq!(qi.query_index(4), 2); // query 2 plus
        assert_eq!(qi.query_index(5), 2); // query 2 minus
    }

    /// Port of NCBI queryinfo_unit_test: context offsets are correctly sequenced.
    #[test]
    fn test_queryinfo_context_offsets() {
        let qi = QueryInfo::new_blastn(&[100, 50]);

        // Each context starts after the previous one ends (+1 for sentinel)
        // Query 0 plus: offset=0, length=100
        assert_eq!(qi.contexts[0].query_offset, 0);
        assert_eq!(qi.contexts[0].query_length, 100);

        // Query 0 minus: offset = 0 + 100 + 1 = 101
        assert_eq!(qi.contexts[1].query_offset, 101);
        assert_eq!(qi.contexts[1].query_length, 100);

        // Query 1 plus: offset = 101 + 100 + 1 = 202
        assert_eq!(qi.contexts[2].query_offset, 202);
        assert_eq!(qi.contexts[2].query_length, 50);

        // Query 1 minus: offset = 202 + 50 + 1 = 253
        assert_eq!(qi.contexts[3].query_offset, 253);
        assert_eq!(qi.contexts[3].query_length, 50);

        // Verify no overlapping: each context starts after previous ends
        for i in 1..qi.num_contexts() {
            let prev = &qi.contexts[i - 1];
            let curr = &qi.contexts[i];
            assert!(
                curr.query_offset >= prev.query_offset + prev.query_length,
                "Context {} overlaps with context {}: offset {} < {} + {}",
                i,
                i - 1,
                curr.query_offset,
                prev.query_offset,
                prev.query_length
            );
        }
    }

    #[test]
    fn test_bsearch_context_info_returns_preceding_context() {
        let qi = QueryInfo::new_blastn(&[10, 5]);

        assert_eq!(bsearch_context_info(-1, &qi), 0);
        assert_eq!(bsearch_context_info(0, &qi), 0);
        assert_eq!(bsearch_context_info(9, &qi), 0);
        assert_eq!(bsearch_context_info(10, &qi), 0);
        assert_eq!(bsearch_context_info(11, &qi), 1);
        assert_eq!(bsearch_context_info(21, &qi), 1);
        assert_eq!(bsearch_context_info(22, &qi), 2);
        assert_eq!(bsearch_context_info(27, &qi), 2);
        assert_eq!(bsearch_context_info(28, &qi), 3);
        assert_eq!(bsearch_context_info(1000, &qi), 3);
    }

    #[test]
    fn test_bsearch_context_info_empty_query_info() {
        let qi = QueryInfo {
            num_queries: 0,
            contexts: Vec::new(),
            max_length: 0,
        };

        assert_eq!(bsearch_context_info(0, &qi), -1);
    }

    #[test]
    fn test_blast_get_query_index_from_query_offset() {
        let blastn = QueryInfo::new_blastn(&[10, 5]);
        assert_eq!(
            blast_get_query_index_from_query_offset(0, crate::program::BLASTN, &blastn),
            0
        );
        assert_eq!(
            blast_get_query_index_from_query_offset(11, crate::program::BLASTN, &blastn),
            0
        );
        assert_eq!(
            blast_get_query_index_from_query_offset(22, crate::program::BLASTN, &blastn),
            1
        );

        let offsets = [0, 3, 7, 11, 14, 18, 22];
        let blastx = QueryInfo::new_translated_query_from_offsets(&offsets);
        assert_eq!(
            blast_get_query_index_from_query_offset(17, crate::program::BLASTX, &blastx),
            0
        );
    }
}
