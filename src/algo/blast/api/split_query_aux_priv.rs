use std::fmt;
use std::sync::Arc;

use crate::algo::blast::api::query_data::{CBlastOptions, EBlastProgramType, ILocalQueryData};
use crate::algo::blast::api::setup_factory::SInternalData;
use crate::algo::blast::api::split_query_blk::CSplitQueryBlk;
use crate::algo::blast::core::blast_program::{
    blast_program_is_phi_blast, blast_query_is_translated, blast_subject_is_pssm,
};
use crate::util::CODON_LENGTH;

pub const K_INVALID_CONTEXT: i32 = -1;

#[derive(Debug, Clone)]
pub struct CContextTranslator {
    pub contexts_per_chunk: Vec<Vec<i32>>,
    pub starting_chunks: Vec<Vec<i32>>,
    pub absolute_contexts: Vec<Vec<i32>>,
    pub split_query_blk: Option<Arc<CSplitQueryBlk>>,
}

#[derive(Debug, Clone)]
pub struct CQueryDataPerChunk {
    pub program: EBlastProgramType,
    pub query_indices_per_chunk: Vec<Vec<usize>>,
    pub query_lengths: Vec<usize>,
    pub last_chunk_for_query_cache: Vec<i32>,
    pub local_query_data: Option<Arc<ILocalQueryData>>,
}

#[derive(Debug, Clone)]
pub struct SplitQueryAuxState {
    pub options: Option<Arc<CBlastOptions>>,
    pub full_data: Option<Arc<SInternalData>>,
    pub num_threaded: usize,
}

pub fn split_query_get_overlap_chunk_size(program: EBlastProgramType) -> usize {
    if let Ok(overlap_sz_str) = std::env::var("OVERLAP_CHUNK_SIZE") {
        if !overlap_sz_str.trim().is_empty() {
            return overlap_sz_str
                .parse::<usize>()
                .expect("OVERLAP_CHUNK_SIZE must be an integer");
        }
    }

    if blast_query_is_translated(program as u32) {
        297
    } else {
        100
    }
}

pub fn split_query_should_split(
    program: EBlastProgramType,
    _chunk_size: usize,
    _concatenated_query_length: usize,
    num_queries: usize,
) -> bool {
    if program == EBlastProgramType::Mapping {
        return false;
    }

    !(blast_subject_is_pssm(program as u32)
        || (program == EBlastProgramType::Blastx && num_queries > 1)
        || blast_program_is_phi_blast(program as u32))
}

pub fn split_query_calculate_num_chunks(
    program: EBlastProgramType,
    chunk_size: &mut usize,
    concatenated_query_length: usize,
    num_queries: usize,
) -> u32 {
    if !split_query_should_split(program, *chunk_size, concatenated_query_length, num_queries) {
        return 1;
    }

    let overlap_size = split_query_get_overlap_chunk_size(program);

    if blast_query_is_translated(program as u32) {
        let chunk_size_delta = *chunk_size % CODON_LENGTH;
        *chunk_size -= chunk_size_delta;
        debug_assert_eq!(*chunk_size % CODON_LENGTH, 0);
    }

    let mut num_chunks = 0usize;
    if *chunk_size > overlap_size {
        num_chunks = concatenated_query_length / (*chunk_size - overlap_size);
    }

    if num_chunks <= 1 {
        *chunk_size = concatenated_query_length;
        return 1;
    }

    if !blast_query_is_translated(program as u32) {
        *chunk_size = (concatenated_query_length + (num_chunks - 1) * overlap_size) / num_chunks;
        if num_chunks < *chunk_size - overlap_size {
            *chunk_size += 1;
        }
    }

    num_chunks as u32
}

impl CContextTranslator {
    pub fn get_absolute_context(&self, chunk_num: usize, context_in_chunk: i32) -> i32 {
        debug_assert!(chunk_num < self.contexts_per_chunk.len());
        debug_assert!(context_in_chunk >= 0);
        debug_assert!((context_in_chunk as usize) < self.contexts_per_chunk[chunk_num].len());
        self.contexts_per_chunk[chunk_num][context_in_chunk as usize]
    }

    pub fn get_context_in_chunk(&self, chunk_num: usize, absolute_context: i32) -> i32 {
        debug_assert!(chunk_num < self.contexts_per_chunk.len());
        let context_indices = &self.contexts_per_chunk[chunk_num];
        context_indices
            .iter()
            .position(|context| *context == absolute_context)
            .map(|index| index as i32)
            .unwrap_or(K_INVALID_CONTEXT)
    }

    pub fn get_starting_chunk(&self, curr_chunk: usize, context_in_chunk: i32) -> i32 {
        let absolute_context = self.get_absolute_context(curr_chunk, context_in_chunk);
        if absolute_context == K_INVALID_CONTEXT {
            return K_INVALID_CONTEXT;
        }

        let mut retval = curr_chunk;
        let mut chunk = curr_chunk;
        while chunk > 0 {
            chunk -= 1;
            if self.get_context_in_chunk(chunk, absolute_context) == K_INVALID_CONTEXT {
                break;
            }
            retval = chunk;
        }
        retval as i32
    }
}

impl fmt::Display for CContextTranslator {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self
            .starting_chunks
            .first()
            .is_none_or(|chunks| chunks.is_empty())
            || self
                .absolute_contexts
                .first()
                .is_none_or(|contexts| contexts.is_empty())
        {
            return Ok(());
        }

        let num_chunks = self.contexts_per_chunk.len();
        writeln!(f)?;
        writeln!(f, "NumChunks = {num_chunks}")?;

        for index in 0..num_chunks {
            write!(f, "Chunk{index}StartingChunks = ")?;
            if let Some(chunks) = self.starting_chunks.get(index) {
                for (offset, value) in chunks.iter().enumerate() {
                    if offset > 0 {
                        write!(f, ", ")?;
                    }
                    write!(f, "{value}")?;
                }
            }
            writeln!(f)?;
        }
        writeln!(f)?;

        for index in 0..num_chunks {
            write!(f, "Chunk{index}AbsoluteContexts = ")?;
            if let Some(contexts) = self.absolute_contexts.get(index) {
                for (offset, value) in contexts.iter().enumerate() {
                    if offset > 0 {
                        write!(f, ", ")?;
                    }
                    write!(f, "{value}")?;
                }
            }
            writeln!(f)?;
        }
        writeln!(f)
    }
}
