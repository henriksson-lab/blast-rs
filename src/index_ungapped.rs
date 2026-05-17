//! Rust port of NCBI `index_ungapped.c`.

pub const IR_HASH_SIZE: usize = 4 * 1024;
pub const FP_ENTRY_SIZE: usize = 1024 * 1024;

/// Port of NCBI `IR_DIAG` (`index_ungapped.c`).
pub fn ir_diag(qoff: u32, soff: u32) -> u32 {
    0x1000_0000u32.wrapping_add(soff).wrapping_sub(qoff)
}

/// Port of NCBI `IR_KEY` (`index_ungapped.c`).
pub fn ir_key(diag: u32) -> u32 {
    diag % IR_HASH_SIZE as u32
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct IrDiagData {
    pub diag: u32,
    pub qend: u32,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct IrHashEntry {
    pub diag_data: IrDiagData,
    pub next: Option<usize>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IrFpEntry {
    pub entries: Vec<IrHashEntry>,
    pub next: Option<usize>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IrDiagHash {
    pub entries: Vec<IrHashEntry>,
    pub free_pool: Vec<IrFpEntry>,
    pub free: Vec<usize>,
}

/// Rust ownership equivalent of NCBI `ir_fp_entry_destroy`
/// (`index_ungapped.c:45`).
pub fn ir_fp_entry_destroy(mut e: Option<IrFpEntry>) -> Option<IrFpEntry> {
    if let Some(entry) = e.as_mut() {
        entry.entries.clear();
        entry.next = None;
    }
    None
}

/// Port of NCBI `ir_fp_entry_create` (`index_ungapped.c:60`).
pub fn ir_fp_entry_create() -> Option<IrFpEntry> {
    let mut entries = vec![IrHashEntry::default(); FP_ENTRY_SIZE];
    for i in 0..FP_ENTRY_SIZE - 1 {
        entries[i].next = Some(i + 1);
    }
    Some(IrFpEntry {
        entries,
        next: None,
    })
}

/// Port of NCBI `ir_hash_create` (`index_ungapped.c:82`).
pub fn ir_hash_create() -> Option<IrDiagHash> {
    Some(IrDiagHash {
        entries: vec![IrHashEntry::default(); IR_HASH_SIZE],
        free_pool: Vec::new(),
        free: Vec::new(),
    })
}

/// Rust ownership equivalent of NCBI `ir_hash_destroy`.
pub fn ir_hash_destroy(mut hash: Option<IrDiagHash>) -> Option<IrDiagHash> {
    if let Some(hash) = hash.as_mut() {
        for entry in hash.free_pool.drain(..) {
            let _ = ir_fp_entry_destroy(Some(entry));
        }
        hash.entries.clear();
        hash.free.clear();
    }
    None
}

/// Port of NCBI `ir_locate` (`index_ungapped.c:119`).
///
/// Returns the index into `hash.entries` for the primary hash-table entry, or
/// `IR_HASH_SIZE + free_pool_index * FP_ENTRY_SIZE + entry_index` for overflow
/// entries allocated from the free-pool blocks.
pub fn ir_locate(hash: &mut IrDiagHash, diag: u32, key: u32) -> Option<usize> {
    let key = key as usize;
    if key >= hash.entries.len() {
        return None;
    }

    let mut current = hash.entries[key].next;
    while let Some(encoded_index) = current {
        let (pool_index, entry_index) = decode_pool_index(encoded_index)?;
        let entry = hash.free_pool.get(pool_index)?.entries.get(entry_index)?;
        if entry.diag_data.diag == diag {
            let tmp = entry.diag_data;
            hash.free_pool[pool_index].entries[entry_index].diag_data = hash.entries[key].diag_data;
            hash.entries[key].diag_data = tmp;
            return Some(key);
        }
        current = entry.next;
    }

    if hash.free.is_empty() {
        let fp = ir_fp_entry_create()?;
        let pool_index = hash.free_pool.len();
        hash.free_pool.push(fp);
        for entry_index in (0..FP_ENTRY_SIZE).rev() {
            hash.free.push(encode_pool_index(pool_index, entry_index));
        }
    }

    let encoded_index = hash.free.pop()?;
    let (pool_index, entry_index) = decode_pool_index(encoded_index)?;
    let old_next = hash.entries[key].next;
    let ce = &mut hash.free_pool[pool_index].entries[entry_index];
    ce.next = old_next;
    ce.diag_data.diag = diag;
    ce.diag_data.qend = 0;
    hash.entries[key].next = Some(encoded_index);
    Some(encoded_index)
}

/// Port of NCBI `IR_LOCATE` (`index_ungapped.c`).
/// naming: Rust includes `_macro` to distinguish this C macro-shaped helper
/// from the out-of-line `ir_locate` routine.
pub fn ir_locate_macro(hash: &mut IrDiagHash, diag: u32, key: u32) -> Option<usize> {
    let key = key as usize;
    let entry = hash.entries.get(key)?;
    if entry.diag_data.qend == 0 || diag == entry.diag_data.diag {
        Some(key)
    } else {
        ir_locate(hash, diag, key as u32)
    }
}

/// Mutable access to the entry index returned by [`ir_locate_macro`].
pub fn ir_hash_entry_mut(hash: &mut IrDiagHash, index: usize) -> Option<&mut IrHashEntry> {
    if index < IR_HASH_SIZE {
        return hash.entries.get_mut(index);
    }
    let (pool_index, entry_index) = decode_pool_index(index)?;
    hash.free_pool
        .get_mut(pool_index)?
        .entries
        .get_mut(entry_index)
}

fn encode_pool_index(pool_index: usize, entry_index: usize) -> usize {
    IR_HASH_SIZE + pool_index * FP_ENTRY_SIZE + entry_index
}

fn decode_pool_index(encoded_index: usize) -> Option<(usize, usize)> {
    let shifted = encoded_index.checked_sub(IR_HASH_SIZE)?;
    Some((shifted / FP_ENTRY_SIZE, shifted % FP_ENTRY_SIZE))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ir_hash_lifecycle_and_locate_match_c_shape() {
        assert_eq!(ir_diag(10, 25), 0x1000_000f);
        assert_eq!(ir_key(0x1000_000f), 15);

        let mut hash = ir_hash_create().expect("hash");
        assert_eq!(hash.entries.len(), IR_HASH_SIZE);
        let key = ir_key(123);
        hash.entries[key as usize].diag_data.qend = 1;
        hash.entries[key as usize].diag_data.diag = 99;

        let overflow = ir_locate(&mut hash, 123, key).expect("overflow");
        assert!(overflow >= IR_HASH_SIZE);
        assert_eq!(hash.free_pool.len(), 1);
        assert_eq!(hash.entries[key as usize].next, Some(overflow));

        hash.entries[key as usize].diag_data = IrDiagData { diag: 99, qend: 7 };
        let promoted = ir_locate(&mut hash, 123, key).expect("promoted");
        assert_eq!(promoted, key as usize);
        assert_eq!(hash.entries[key as usize].diag_data.diag, 123);
        assert_eq!(hash.entries[key as usize].diag_data.qend, 0);

        assert_eq!(ir_locate_macro(&mut hash, 123, key), Some(key as usize));
        assert_eq!(ir_locate_macro(&mut hash, 999, IR_HASH_SIZE as u32), None);
        hash.entries[key as usize].diag_data = IrDiagData { diag: 456, qend: 0 };
        assert_eq!(ir_locate_macro(&mut hash, 789, key), Some(key as usize));
        assert_eq!(hash.entries[key as usize].diag_data.diag, 456);

        assert!(ir_hash_destroy(Some(hash)).is_none());
        assert!(ir_fp_entry_destroy(ir_fp_entry_create()).is_none());
    }
}
