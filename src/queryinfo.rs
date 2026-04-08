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
}
