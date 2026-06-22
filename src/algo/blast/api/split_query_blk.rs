//! Rust port of NCBI `split_query_blk.cpp/.hpp`.

use std::fmt;
use std::ops::Range;

use super::split_query::{
    split_query_blk_add_context_offset_to_chunk, split_query_blk_add_context_to_chunk,
    split_query_blk_add_query_to_chunk, split_query_blk_free, split_query_blk_get_chunk_bounds,
    split_query_blk_get_chunk_overlap_size, split_query_blk_get_context_offsets_for_chunk,
    split_query_blk_get_num_queries_for_chunk, split_query_blk_get_query_contexts_for_chunk,
    split_query_blk_get_query_indices_for_chunk, split_query_blk_new,
    split_query_blk_set_chunk_bounds, split_query_blk_set_chunk_overlap_size, SSplitQueryBlk,
};

/// Range describing a query chunk.
///
/// NCBI C++ typedef: `COpenRange<TSeqPos> TChunkRange`.
pub type TChunkRange = Range<usize>;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SplitQueryBlkError {
    OutOfMemory,
    CoreCall(&'static str),
}

impl fmt::Display for SplitQueryBlkError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SplitQueryBlkError::OutOfMemory => f.write_str("SplitQueryBlkNew"),
            SplitQueryBlkError::CoreCall(name) => f.write_str(name),
        }
    }
}

impl std::error::Error for SplitQueryBlkError {}

/// Wrapper class around `SSplitQueryBlk`.
///
/// NCBI C++ class: `CSplitQueryBlk`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CSplitQueryBlk {
    split_query_blk: SSplitQueryBlk,
}

impl CSplitQueryBlk {
    pub fn new(num_chunks: u32, gapped_merge: bool) -> Result<Self, SplitQueryBlkError> {
        let Some(split_query_blk) = split_query_blk_new(num_chunks, gapped_merge) else {
            return Err(SplitQueryBlkError::OutOfMemory);
        };
        Ok(Self { split_query_blk })
    }

    pub fn get_num_chunks(&self) -> usize {
        self.split_query_blk.num_chunks as usize
    }

    pub fn get_num_queries_for_chunk(&self, chunk_num: usize) -> Result<usize, SplitQueryBlkError> {
        let (rv, retval) = split_query_blk_get_num_queries_for_chunk(
            Some(&self.split_query_blk),
            chunk_num as u32,
        );
        if rv != 0 {
            return Err(SplitQueryBlkError::CoreCall(
                "SplitQueryBlk_GetNumQueriesForChunk",
            ));
        }
        Ok(retval.expect("core call returned success without output"))
    }

    pub fn get_query_indices(&self, chunk_num: usize) -> Result<Vec<usize>, SplitQueryBlkError> {
        let (rv, query_indices) = split_query_blk_get_query_indices_for_chunk(
            Some(&self.split_query_blk),
            chunk_num as u32,
        );
        if rv != 0 {
            return Err(SplitQueryBlkError::CoreCall(
                "SplitQueryBlk_GetQueryIndicesForChunk",
            ));
        }
        let retval = query_indices
            .expect("core call returned success without output")
            .into_iter()
            .take_while(|&index| index != u32::MAX)
            .map(|index| index as usize)
            .collect();
        Ok(retval)
    }

    pub fn get_query_contexts(&self, chunk_num: usize) -> Result<Vec<i32>, SplitQueryBlkError> {
        let (rv, query_contexts) = split_query_blk_get_query_contexts_for_chunk(
            Some(&self.split_query_blk),
            chunk_num as u32,
        );
        if rv != 0 {
            return Err(SplitQueryBlkError::CoreCall(
                "SplitQueryBlk_GetQueryContextsForChunk",
            ));
        }
        let (retval, _num_contexts) =
            query_contexts.expect("core call returned success without output");
        Ok(retval)
    }

    pub fn get_context_offsets(&self, chunk_num: usize) -> Result<Vec<usize>, SplitQueryBlkError> {
        let (rv, context_offsets) = split_query_blk_get_context_offsets_for_chunk(
            Some(&self.split_query_blk),
            chunk_num as u32,
        );
        if rv != 0 {
            return Err(SplitQueryBlkError::CoreCall(
                "SplitQueryBlk_GetContextOffsetsForChunk",
            ));
        }
        let retval = context_offsets
            .expect("core call returned success without output")
            .into_iter()
            .take_while(|&offset| offset != u32::MAX)
            .map(|offset| offset as usize)
            .collect();
        Ok(retval)
    }

    pub fn get_chunk_bounds(&self, chunk_num: usize) -> Result<TChunkRange, SplitQueryBlkError> {
        let (rv, chunk_bounds) =
            split_query_blk_get_chunk_bounds(Some(&self.split_query_blk), chunk_num as u32);
        if rv != 0 {
            return Err(SplitQueryBlkError::CoreCall("SplitQueryBlk_GetChunkBounds"));
        }
        let (first, second) = chunk_bounds.expect("core call returned success without output");
        Ok(first..second)
    }

    pub fn set_chunk_bounds(
        &mut self,
        chunk_num: usize,
        chunk_range: TChunkRange,
    ) -> Result<(), SplitQueryBlkError> {
        let rv = split_query_blk_set_chunk_bounds(
            Some(&mut self.split_query_blk),
            chunk_num as u32,
            chunk_range.start as u32,
            chunk_range.end as u32,
        );
        if rv != 0 {
            return Err(SplitQueryBlkError::CoreCall("SplitQueryBlk_SetChunkBounds"));
        }
        Ok(())
    }

    pub fn add_query_to_chunk(
        &mut self,
        chunk_num: usize,
        query_index: i32,
    ) -> Result<(), SplitQueryBlkError> {
        let rv = split_query_blk_add_query_to_chunk(
            Some(&mut self.split_query_blk),
            query_index as u32,
            chunk_num as u32,
        );
        if rv != 0 {
            return Err(SplitQueryBlkError::CoreCall(
                "Failed to add query to SplitQueryBlk",
            ));
        }
        Ok(())
    }

    pub fn add_context_to_chunk(
        &mut self,
        chunk_num: usize,
        context_index: i32,
    ) -> Result<(), SplitQueryBlkError> {
        let rv = split_query_blk_add_context_to_chunk(
            Some(&mut self.split_query_blk),
            context_index,
            chunk_num as u32,
        );
        if rv != 0 {
            return Err(SplitQueryBlkError::CoreCall(
                "Failed to add context to SplitQueryBlk",
            ));
        }
        Ok(())
    }

    pub fn add_context_offset_to_chunk(
        &mut self,
        chunk_num: usize,
        context_offset: i32,
    ) -> Result<(), SplitQueryBlkError> {
        let rv = split_query_blk_add_context_offset_to_chunk(
            Some(&mut self.split_query_blk),
            context_offset as u32,
            chunk_num as u32,
        );
        if rv != 0 {
            return Err(SplitQueryBlkError::CoreCall(
                "Failed to add context offset to SplitQueryBlk",
            ));
        }
        Ok(())
    }

    pub fn get_c_struct(&self) -> &SSplitQueryBlk {
        &self.split_query_blk
    }

    pub fn get_c_struct_mut(&mut self) -> &mut SSplitQueryBlk {
        &mut self.split_query_blk
    }

    pub fn set_chunk_overlap_size(&mut self, size: usize) -> Result<(), SplitQueryBlkError> {
        let rv = split_query_blk_set_chunk_overlap_size(Some(&mut self.split_query_blk), size);
        if rv != 0 {
            return Err(SplitQueryBlkError::CoreCall(
                "Failed to set chunk overlap size in SplitQueryBlk",
            ));
        }
        Ok(())
    }

    pub fn get_chunk_overlap_size(&self) -> usize {
        let retval = split_query_blk_get_chunk_overlap_size(Some(&self.split_query_blk));
        if retval == 0 {
            eprintln!("Warning: Query-splitting Chunk overlap size was not set");
        }
        retval
    }
}

impl Drop for CSplitQueryBlk {
    fn drop(&mut self) {
        let split_query_blk = std::mem::take(&mut self.split_query_blk);
        let _ = split_query_blk_free(Some(split_query_blk));
    }
}

impl fmt::Display for CSplitQueryBlk {
    fn fmt(&self, out: &mut fmt::Formatter<'_>) -> fmt::Result {
        let num_chunks = self.get_num_chunks();

        writeln!(out)?;
        writeln!(out, "NumChunks = {num_chunks}")?;
        for chunk_num in 0..num_chunks {
            writeln!(
                out,
                "Chunk{chunk_num}Queries = {}",
                self.get_query_indices(chunk_num)
                    .map_err(|_| fmt::Error)?
                    .iter()
                    .map(|value| value.to_string())
                    .collect::<Vec<_>>()
                    .join(", ")
            )?;
        }
        writeln!(out)?;
        for chunk_num in 0..num_chunks {
            writeln!(
                out,
                "Chunk{chunk_num}Contexts = {}",
                self.get_query_contexts(chunk_num)
                    .map_err(|_| fmt::Error)?
                    .iter()
                    .map(|value| value.to_string())
                    .collect::<Vec<_>>()
                    .join(", ")
            )?;
        }
        writeln!(out)?;
        for chunk_num in 0..num_chunks {
            writeln!(
                out,
                "Chunk{chunk_num}ContextOffsets = {}",
                self.get_context_offsets(chunk_num)
                    .map_err(|_| fmt::Error)?
                    .iter()
                    .map(|value| value.to_string())
                    .collect::<Vec<_>>()
                    .join(", ")
            )?;
        }

        Ok(())
    }
}
