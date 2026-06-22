use std::any::Any;
use std::sync::Arc;

use crate::algo::blast::api::query_data::{
    EBlastProgramType, IBlastQuerySource, IQueryFactory, TSeqLocVector,
};
use crate::algo::blast::core::blast_encoding::EBlastEncoding;
use crate::algo::blast::core::blast_seqsrc::{
    BlastSeqSrc, BlastSeqSrcGetSeqArg, BlastSeqSrcIterator, _blast_seq_src_impl_set_data_structure,
    BLAST_SEQSRC_EOF, BLAST_SEQSRC_ERROR, BLAST_SEQSRC_SUCCESS,
};
use crate::algo::blast::core::blast_util::BLAST_SequenceBlk;

#[derive(Debug, Clone)]
pub struct CQueryFactoryInfo {
    pub is_prot: bool,
    pub seq_blk_vector: Vec<*mut BLAST_SequenceBlk>,
    pub max_length: u32,
    pub min_length: u32,
    pub avg_length: u32,
    pub query_source: Option<Arc<IBlastQuerySource>>,
    pub num_seqs: u32,
}

#[derive(Debug, Clone)]
pub struct SQueryFactorySrcNewArgs {
    pub query_factory: Option<Arc<IQueryFactory>>,
    pub subj_seqs: TSeqLocVector,
    pub program: EBlastProgramType,
}

impl CQueryFactoryInfo {
    pub fn get_max_length(&self) -> u32 {
        self.max_length
    }

    pub fn get_min_length(&self) -> u32 {
        self.min_length
    }

    pub fn get_avg_length(&self) -> u32 {
        self.avg_length
    }

    pub fn set_avg_length(&mut self, length: u32) {
        self.avg_length = length;
    }

    pub fn get_is_protein(&self) -> bool {
        self.is_prot
    }

    pub fn get_num_seqs(&self) -> u32 {
        self.num_seqs
    }

    pub fn get_seq_blk(&self, index: u32) -> *mut BLAST_SequenceBlk {
        if index >= self.get_num_seqs() {
            panic!("sequence index out of range");
        }
        self.seq_blk_vector[index as usize]
    }
}

pub fn s_query_factory_get_max_length(
    multiseq_handle: &dyn Any,
    _args: Option<&mut dyn Any>,
) -> i32 {
    multiseq_handle
        .downcast_ref::<CQueryFactoryInfo>()
        .map(|seq_info| seq_info.get_max_length() as i32)
        .unwrap_or(0)
}

pub fn s_query_factory_get_min_length(
    multiseq_handle: &dyn Any,
    _args: Option<&mut dyn Any>,
) -> i32 {
    multiseq_handle
        .downcast_ref::<CQueryFactoryInfo>()
        .map(|seq_info| seq_info.get_min_length() as i32)
        .unwrap_or(0)
}

pub fn s_query_factory_get_avg_length(
    multiseq_handle: &dyn Any,
    _args: Option<&mut dyn Any>,
) -> i32 {
    let Some(seq_info) = multiseq_handle.downcast_ref::<CQueryFactoryInfo>() else {
        return 0;
    };

    if seq_info.get_avg_length() > 0 {
        return seq_info.get_avg_length() as i32;
    }

    let num_seqs = seq_info.get_num_seqs();
    if num_seqs == 0 {
        return 0;
    }

    let total_length: i64 = seq_info
        .seq_blk_vector
        .iter()
        .filter_map(|&seq_blk| unsafe { seq_blk.as_ref() })
        .map(|seq_blk| seq_blk.length as i64)
        .sum();
    (total_length / num_seqs as i64) as i32
}

pub fn s_query_factory_get_num_seqs(multiseq_handle: &dyn Any, _args: Option<&mut dyn Any>) -> i32 {
    multiseq_handle
        .downcast_ref::<CQueryFactoryInfo>()
        .map(|seq_info| seq_info.get_num_seqs() as i32)
        .unwrap_or(0)
}

pub fn s_query_factory_get_num_seqs_stats(
    _multiseq_handle: &dyn Any,
    _args: Option<&mut dyn Any>,
) -> i32 {
    0
}

pub fn s_query_factory_get_tot_len(_multiseq_handle: &dyn Any, _args: Option<&mut dyn Any>) -> i64 {
    0
}

pub fn s_query_factory_get_tot_len_stats(
    _multiseq_handle: &dyn Any,
    _args: Option<&mut dyn Any>,
) -> i64 {
    0
}

pub fn s_query_factory_get_name(
    _multiseq_handle: &dyn Any,
    _args: Option<&mut dyn Any>,
) -> Option<String> {
    Some(String::new())
}

pub fn s_query_factory_get_is_prot(multiseq_handle: &dyn Any, _args: Option<&mut dyn Any>) -> bool {
    multiseq_handle
        .downcast_ref::<CQueryFactoryInfo>()
        .map(|seq_info| seq_info.get_is_protein())
        .unwrap_or(false)
}

pub fn s_query_factory_get_sequence(
    multiseq_handle: &dyn Any,
    args: &mut BlastSeqSrcGetSeqArg,
) -> i16 {
    let Some(seq_info) = multiseq_handle.downcast_ref::<CQueryFactoryInfo>() else {
        return BLAST_SEQSRC_ERROR;
    };
    if seq_info.num_seqs == 0 {
        return BLAST_SEQSRC_ERROR;
    }

    let index = args.oid;
    if index < 0 || index as u32 >= seq_info.num_seqs {
        return BLAST_SEQSRC_EOF as i16;
    }

    let seq_blk_ptr = seq_info.seq_blk_vector[index as usize];
    let Some(seq_blk) = (unsafe { seq_blk_ptr.as_ref() }) else {
        return BLAST_SEQSRC_ERROR;
    };
    let mut copied = seq_blk.clone();
    if args.encoding == EBlastEncoding::Nucleotide {
        copied.sequence = copied
            .sequence_start
            .as_ref()
            .and_then(|sequence_start| sequence_start.get(1..).map(|seq| seq.to_vec()));
    } else if args.encoding == EBlastEncoding::Ncbi4na {
        copied.sequence = copied.sequence_start.clone();
    }
    copied.lcase_mask = None;
    copied.lcase_mask_allocated = false;
    copied.oid = args.oid;

    unsafe {
        if args.seq.is_null() {
            args.seq = Box::into_raw(Box::new(copied));
        } else {
            *args.seq = copied;
        }
    }

    BLAST_SEQSRC_SUCCESS
}

pub fn s_query_factory_release_sequence(
    _multiseq_handle: &dyn Any,
    args: &mut BlastSeqSrcGetSeqArg,
) {
    if args.seq.is_null() {
        return;
    }
    unsafe {
        if (*args.seq).sequence_start_allocated {
            (*args.seq).sequence_start = None;
            (*args.seq).sequence_start_allocated = false;
        }
    }
}

pub fn s_query_factory_get_seq_len(multiseq_handle: &dyn Any, args: Option<&mut dyn Any>) -> i32 {
    let Some(seq_info) = multiseq_handle.downcast_ref::<CQueryFactoryInfo>() else {
        return BLAST_SEQSRC_ERROR as i32;
    };
    let Some(index) = args.and_then(|args| args.downcast_ref::<i32>()).copied() else {
        return BLAST_SEQSRC_ERROR as i32;
    };
    if index < 0 {
        return BLAST_SEQSRC_ERROR as i32;
    }
    unsafe {
        seq_info
            .seq_blk_vector
            .get(index as usize)
            .and_then(|&seq_blk| seq_blk.as_ref())
            .map(|seq_blk| seq_blk.length)
            .unwrap_or(BLAST_SEQSRC_ERROR as i32)
    }
}

pub fn s_query_factory_iterator_next(
    multiseq_handle: &dyn Any,
    itr: &mut BlastSeqSrcIterator,
) -> i32 {
    let status = s_query_factory_get_next_chunk(multiseq_handle, itr);
    if status == BLAST_SEQSRC_EOF {
        return status as i32;
    }

    let retval = itr.current_pos as i32;
    itr.current_pos += 1;
    retval
}

pub fn s_query_factory_get_next_chunk(
    multiseq_handle: &dyn Any,
    itr: &mut BlastSeqSrcIterator,
) -> i16 {
    let Some(seq_info) = multiseq_handle.downcast_ref::<CQueryFactoryInfo>() else {
        return BLAST_SEQSRC_EOF;
    };

    if itr.current_pos == u32::MAX {
        itr.current_pos = 0;
    }
    if itr.current_pos >= seq_info.num_seqs {
        return BLAST_SEQSRC_EOF;
    }

    BLAST_SEQSRC_SUCCESS
}

pub fn s_query_factory_reset_chunk_iter(_multiseq_handle: &dyn Any) {}

pub fn s_query_factory_src_free(seq_src: &mut BlastSeqSrc) -> Option<()> {
    _blast_seq_src_impl_set_data_structure(Some(seq_src), None);
    None
}

pub fn s_query_factory_src_copy(seq_src: &mut BlastSeqSrc) -> Option<()> {
    let data_structure = seq_src.data_structure.clone()?;
    _blast_seq_src_impl_set_data_structure(Some(seq_src), Some(data_structure));
    Some(())
}
