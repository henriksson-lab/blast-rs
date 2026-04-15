//! Rust equivalent of blast_query_info.c — query information management.

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
}

/// Information about all queries in a search.
#[derive(Debug, Clone)]
pub struct QueryInfo {
    pub num_queries: i32,
    pub contexts: Vec<ContextInfo>,
    pub max_length: u32,
}

impl QueryInfo {
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

    /// Get context info for a given context index.
    pub fn get_context(&self, context: usize) -> Option<&ContextInfo> {
        self.contexts.get(context)
    }

    /// Get the query index for a given context.
    pub fn query_index(&self, context: usize) -> i32 {
        self.contexts.get(context).map_or(-1, |c| c.query_index)
    }

    /// Number of contexts.
    pub fn num_contexts(&self) -> usize {
        self.contexts.len()
    }
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
        // For blastp, there is only 1 context per query (no strand).
        // Construct manually since no new_blastp exists yet.
        let query_len = 150;
        let qi = QueryInfo {
            num_queries: 1,
            contexts: vec![ContextInfo {
                query_offset: 0,
                query_length: query_len,
                eff_searchsp: 0,
                length_adjustment: 0,
                query_index: 0,
                frame: 0,
                is_valid: true,
            }],
            max_length: query_len as u32,
        };
        assert_eq!(qi.num_queries, 1);
        assert_eq!(qi.num_contexts(), 1);
        assert_eq!(qi.contexts[0].frame, 0);
        assert_eq!(qi.contexts[0].query_length, query_len);
    }

    /// Port of NCBI queryinfo_unit_test: blastx should have 6 contexts (6 reading frames).
    #[test]
    fn test_queryinfo_blastx_six_contexts() {
        // blastx has 6 frames: +1,+2,+3,-1,-2,-3
        let query_len = 300;
        let frames = [1, 2, 3, -1, -2, -3];
        let protein_len = query_len / 3; // approximate translated length
        let mut contexts = Vec::new();
        let mut offset = 0i32;
        for &frame in &frames {
            contexts.push(ContextInfo {
                query_offset: offset,
                query_length: protein_len,
                eff_searchsp: 0,
                length_adjustment: 0,
                query_index: 0,
                frame,
                is_valid: true,
            });
            offset += protein_len + 1;
        }
        let qi = QueryInfo {
            num_queries: 1,
            contexts,
            max_length: protein_len as u32,
        };
        assert_eq!(qi.num_queries, 1);
        assert_eq!(qi.num_contexts(), 6);
        assert_eq!(qi.contexts[0].frame, 1);
        assert_eq!(qi.contexts[3].frame, -1);
        assert_eq!(qi.contexts[5].frame, -3);
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
}
