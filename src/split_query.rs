//! Rust port of NCBI `split_query.c`.

pub const K_BAD_PARAMETER: i16 = -1;
pub const K_OUT_OF_MEMORY: i16 = -2;

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SQueryChunkBoundary {
    pub left: u32,
    pub right: u32,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct SSplitQueryBlk {
    pub num_chunks: u32,
    pub chunk_query_map: Vec<Vec<u32>>,
    pub chunk_ctx_map: Vec<Vec<i32>>,
    pub chunk_offset_map: Vec<Vec<u32>>,
    pub chunk_bounds: Vec<SQueryChunkBoundary>,
    pub chunk_overlap_sz: usize,
    pub gapped_merge: bool,
}

fn valid_chunk(squery_blk: Option<&SSplitQueryBlk>, chunk_num: u32) -> bool {
    squery_blk
        .map(|blk| chunk_num < blk.num_chunks)
        .unwrap_or(false)
}

/// Port of NCBI `SplitQueryBlkNew` (`split_query.c:41`).
pub fn split_query_blk_new(num_chunks: u32, gapped_merge: bool) -> Option<SSplitQueryBlk> {
    if num_chunks == 0 {
        return None;
    }
    Some(SSplitQueryBlk {
        num_chunks,
        chunk_query_map: (0..num_chunks).map(|_| Vec::new()).collect(),
        chunk_ctx_map: (0..num_chunks).map(|_| Vec::new()).collect(),
        chunk_offset_map: (0..num_chunks).map(|_| Vec::new()).collect(),
        chunk_bounds: vec![SQueryChunkBoundary::default(); num_chunks as usize],
        chunk_overlap_sz: 0,
        gapped_merge,
    })
}

/// Rust ownership equivalent of NCBI `SplitQueryBlkFree`.
pub fn split_query_blk_free(mut squery_blk: Option<SSplitQueryBlk>) -> Option<SSplitQueryBlk> {
    if let Some(blk) = squery_blk.as_mut() {
        for queries in &mut blk.chunk_query_map {
            queries.clear();
        }
        for contexts in &mut blk.chunk_ctx_map {
            contexts.clear();
        }
        for offsets in &mut blk.chunk_offset_map {
            offsets.clear();
        }
        blk.chunk_bounds.clear();
        blk.num_chunks = 0;
    }
    None
}

/// Port of NCBI `SplitQueryBlk_AllowGap` (`split_query.c:154`).
pub fn split_query_blk_allow_gap(squery_blk: Option<&SSplitQueryBlk>) -> bool {
    squery_blk.map(|blk| blk.gapped_merge).unwrap_or(false)
}

/// Port of NCBI `SplitQueryBlk_SetChunkBounds` (`split_query.c:160`).
pub fn split_query_blk_set_chunk_bounds(
    squery_blk: Option<&mut SSplitQueryBlk>,
    chunk_num: u32,
    starting_offset: u32,
    ending_offset: u32,
) -> i16 {
    let Some(squery_blk) = squery_blk else {
        return K_BAD_PARAMETER;
    };
    if chunk_num >= squery_blk.num_chunks {
        return K_BAD_PARAMETER;
    }
    squery_blk.chunk_bounds[chunk_num as usize].left = starting_offset;
    squery_blk.chunk_bounds[chunk_num as usize].right = ending_offset;
    0
}

/// Port of NCBI `SplitQueryBlk_GetChunkBounds` (`split_query.c:174`).
pub fn split_query_blk_get_chunk_bounds(
    squery_blk: Option<&SSplitQueryBlk>,
    chunk_num: u32,
) -> (i16, Option<(usize, usize)>) {
    if !valid_chunk(squery_blk, chunk_num) {
        return (K_BAD_PARAMETER, None);
    }
    let squery_blk = squery_blk.expect("checked above");
    let bounds = squery_blk.chunk_bounds[chunk_num as usize];
    (0, Some((bounds.left as usize, bounds.right as usize)))
}

/// Port of NCBI `SplitQueryBlk_GetNumQueriesForChunk` (`split_query.c:189`).
pub fn split_query_blk_get_num_queries_for_chunk(
    squery_blk: Option<&SSplitQueryBlk>,
    chunk_num: u32,
) -> (i16, Option<usize>) {
    if !valid_chunk(squery_blk, chunk_num) {
        return (K_BAD_PARAMETER, None);
    }
    let squery_blk = squery_blk.expect("checked above");
    (
        0,
        Some(squery_blk.chunk_query_map[chunk_num as usize].len()),
    )
}

/// Port of NCBI `SplitQueryBlk_AddQueryToChunk` (`split_query.c:201`).
pub fn split_query_blk_add_query_to_chunk(
    squery_blk: Option<&mut SSplitQueryBlk>,
    query_index: u32,
    chunk_num: u32,
) -> i16 {
    let Some(squery_blk) = squery_blk else {
        return K_BAD_PARAMETER;
    };
    if chunk_num >= squery_blk.num_chunks {
        return K_BAD_PARAMETER;
    }
    squery_blk.chunk_query_map[chunk_num as usize].push(query_index);
    0
}

/// Port of NCBI `SplitQueryBlk_AddContextToChunk` (`split_query.c:211`).
pub fn split_query_blk_add_context_to_chunk(
    squery_blk: Option<&mut SSplitQueryBlk>,
    ctx_index: i32,
    chunk_num: u32,
) -> i16 {
    let Some(squery_blk) = squery_blk else {
        return K_BAD_PARAMETER;
    };
    if chunk_num >= squery_blk.num_chunks {
        return K_BAD_PARAMETER;
    }
    squery_blk.chunk_ctx_map[chunk_num as usize].push(ctx_index);
    0
}

/// Port of NCBI `SplitQueryBlk_AddContextOffsetToChunk` (`split_query.c:222`).
pub fn split_query_blk_add_context_offset_to_chunk(
    squery_blk: Option<&mut SSplitQueryBlk>,
    offset: u32,
    chunk_num: u32,
) -> i16 {
    let Some(squery_blk) = squery_blk else {
        return K_BAD_PARAMETER;
    };
    if chunk_num >= squery_blk.num_chunks {
        return K_BAD_PARAMETER;
    }
    squery_blk.chunk_offset_map[chunk_num as usize].push(offset);
    0
}

/// Port of NCBI `SplitQueryBlk_GetQueryIndicesForChunk` (`split_query.c:235`).
pub fn split_query_blk_get_query_indices_for_chunk(
    squery_blk: Option<&SSplitQueryBlk>,
    chunk_num: u32,
) -> (i16, Option<Vec<u32>>) {
    if !valid_chunk(squery_blk, chunk_num) {
        return (K_BAD_PARAMETER, None);
    }
    let squery_blk = squery_blk.expect("checked above");
    let mut retval = squery_blk.chunk_query_map[chunk_num as usize].clone();
    retval.push(u32::MAX);
    (0, Some(retval))
}

/// Port of NCBI `SplitQueryBlk_GetQueryContextsForChunk` (`split_query.c:262`).
pub fn split_query_blk_get_query_contexts_for_chunk(
    squery_blk: Option<&SSplitQueryBlk>,
    chunk_num: u32,
) -> (i16, Option<(Vec<i32>, u32)>) {
    if !valid_chunk(squery_blk, chunk_num) {
        return (K_BAD_PARAMETER, None);
    }
    let squery_blk = squery_blk.expect("checked above");
    let retval = squery_blk.chunk_ctx_map[chunk_num as usize].clone();
    let num_query_contexts = retval.len() as u32;
    (0, Some((retval, num_query_contexts)))
}

/// Port of NCBI `SplitQueryBlk_GetContextOffsetsForChunk` (`split_query.c:293`).
pub fn split_query_blk_get_context_offsets_for_chunk(
    squery_blk: Option<&SSplitQueryBlk>,
    chunk_num: u32,
) -> (i16, Option<Vec<u32>>) {
    if !valid_chunk(squery_blk, chunk_num) {
        return (K_BAD_PARAMETER, None);
    }
    let squery_blk = squery_blk.expect("checked above");
    let mut retval = squery_blk.chunk_offset_map[chunk_num as usize].clone();
    retval.push(u32::MAX);
    (0, Some(retval))
}

/// Port of NCBI `SplitQueryBlk_SetChunkOverlapSize` (`split_query.c:322`).
pub fn split_query_blk_set_chunk_overlap_size(
    squery_blk: Option<&mut SSplitQueryBlk>,
    size: usize,
) -> i16 {
    let Some(squery_blk) = squery_blk else {
        return K_BAD_PARAMETER;
    };
    squery_blk.chunk_overlap_sz = size;
    0
}

/// Port of NCBI `SplitQueryBlk_GetChunkOverlapSize` (`split_query.c:332`).
pub fn split_query_blk_get_chunk_overlap_size(squery_blk: Option<&SSplitQueryBlk>) -> usize {
    squery_blk
        .map(|blk| blk.chunk_overlap_sz)
        .unwrap_or(K_BAD_PARAMETER as usize)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn split_query_block_lifecycle_and_accessors_match_c_shape() {
        assert!(split_query_blk_new(0, false).is_none());
        let mut blk = split_query_blk_new(2, true).expect("split block");
        assert!(split_query_blk_allow_gap(Some(&blk)));
        assert_eq!(
            split_query_blk_set_chunk_bounds(Some(&mut blk), 1, 10, 25),
            0
        );
        assert_eq!(
            split_query_blk_get_chunk_bounds(Some(&blk), 1),
            (0, Some((10, 25)))
        );
        assert_eq!(
            split_query_blk_get_chunk_bounds(Some(&blk), 2),
            (K_BAD_PARAMETER, None)
        );

        assert_eq!(split_query_blk_add_query_to_chunk(Some(&mut blk), 7, 1), 0);
        assert_eq!(
            split_query_blk_add_context_to_chunk(Some(&mut blk), -3, 1),
            0
        );
        assert_eq!(
            split_query_blk_add_context_offset_to_chunk(Some(&mut blk), 99, 1),
            0
        );
        assert_eq!(
            split_query_blk_get_num_queries_for_chunk(Some(&blk), 1),
            (0, Some(1))
        );
        assert_eq!(
            split_query_blk_get_query_indices_for_chunk(Some(&blk), 1),
            (0, Some(vec![7, u32::MAX]))
        );
        assert_eq!(
            split_query_blk_get_query_contexts_for_chunk(Some(&blk), 1),
            (0, Some((vec![-3], 1)))
        );
        assert_eq!(
            split_query_blk_get_context_offsets_for_chunk(Some(&blk), 1),
            (0, Some(vec![99, u32::MAX]))
        );
        assert_eq!(
            split_query_blk_set_chunk_overlap_size(Some(&mut blk), 31),
            0
        );
        assert_eq!(split_query_blk_get_chunk_overlap_size(Some(&blk)), 31);
        assert_eq!(
            split_query_blk_set_chunk_overlap_size(None, 31),
            K_BAD_PARAMETER
        );
        assert!(split_query_blk_free(Some(blk)).is_none());
    }
}
