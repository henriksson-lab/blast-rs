# NCBI BLAST+ C → blast-rs Function Map

Systematic comparison of every function in the NCBI blastp path against its blast-rs equivalent.
Each function must match the C original or the difference must be documented.

## Status: IN PROGRESS

## NCBI C Functions (in call order for blastp)

### 1. Lookup Table Construction
| NCBI C Function | File | blast-rs Equivalent | Status |
|----------------|------|-------------------|--------|
| `BlastAaLookupTableNew` | blast_aalookup.c | `ProteinLookupTable::build()` | ⚠️ Check |
| `BlastAaLookupFinalize` | blast_aalookup.c | (inline in build) | ⚠️ Check |
| `_ComputeIndex` / `ComputeTableIndex` | blast_lookup.h | `word_hash()` | ✅ Matched (shift-based) |
| `ComputeTableIndexIncremental` | blast_lookup.h | Inline in scan loop | ✅ Matched |
| PV construction | blast_aalookup.c | Inline in build | ⚠️ Check PV size |

### 2. Subject Scanning
| NCBI C Function | File | blast-rs Equivalent | Status |
|----------------|------|-------------------|--------|
| `s_BlastAaScanSubject` | blast_aascan.c | (inline in scan loop) | ✅ Matched |
| `PV_TEST` macro | blast_extend.h | `pv[hash>>6] & (1<<(hash&63))` | ✅ Matched |

### 3. Two-Hit Word Finding
| NCBI C Function | File | blast-rs Equivalent | Status |
|----------------|------|-------------------|--------|
| `s_BlastAaWordFinder_TwoHit` | aa_ungapped.c:440 | `protein_scan_with_table_reuse()` | ⚠️ Partial |
| Diagonal array (`DiagStruct`) | blast_extend.h:57 | `diag_buf: Vec<i32>` | ⚠️ Different encoding |
| `diag_offset` trick | aa_ungapped.c:316 | NOT IMPLEMENTED | ❌ Missing |
| Two-hit window check | aa_ungapped.c:373 | `diff >= TWO_HIT_WINDOW` | ✅ Matched |
| Word overlap check | aa_ungapped.c:384 | `diff < ws` | ✅ Matched |

### 4. Ungapped Extension
| NCBI C Function | File | blast-rs Equivalent | Status |
|----------------|------|-------------------|--------|
| `s_BlastAaExtendTwoHit` | aa_ungapped.c:1089 | New: `extend_left`/`extend_right` | ⚠️ Just implemented |
| Word-best-start scan | aa_ungapped.c:1108-1118 | Inline in scan loop | ✅ Matched |
| Left must reach first hit | aa_ungapped.c:1139 | `reached_first` check | ✅ Matched |
| Right extend with init_score | aa_ungapped.c:1147 | `extend_right(..., left_score)` | ✅ Matched |
| `s_BlastAaExtendLeft` | aa_ungapped.c:903 | `extend_left()` | ⚠️ Check |
| `s_BlastAaExtendRight` | aa_ungapped.c:831 | `extend_right()` | ⚠️ Check |
| `s_BlastAaExtendOneHit` | aa_ungapped.c:1071 | `protein_ungapped_extend()` | ⚠️ Not used in two-hit |

### 5. Constants & Parameters
| NCBI Constant | Value | blast-rs Value | Status |
|--------------|-------|---------------|--------|
| `BLAST_UNGAPPED_X_DROPOFF_PROT` | 7 bits | 40 raw (should be 19 raw) | ❌ WRONG |
| `BLAST_GAP_X_DROPOFF_PROT` | 15 bits | 15 bits (converted to raw) | ✅ Matched |
| `BLAST_GAP_X_DROPOFF_FINAL_PROT` | 25 bits | 25 bits | ✅ Matched |
| `BLAST_GAP_TRIGGER_PROT` | 22 bits | 22 raw (should be 57 raw) | ❌ WRONG |
| `BLAST_WORD_THRESHOLD_BLASTP` | 11 | 11 | ✅ Matched |
| `AA_HITS_PER_CELL` | 3 | 3 | ✅ Matched |
| Two-hit window | 40 | 40 | ✅ Matched |
| charsize | 5 bits | 5 bits (CHARSIZE) | ✅ Matched |

### 6. Gapped Alignment
| NCBI C Function | File | blast-rs Equivalent | Status |
|----------------|------|-------------------|--------|
| `BlastGetStartForGappedAlignment` | blast_gapalign.c | Midpoint of ungapped hit | ⚠️ Different |
| `BLAST_GappedAlignmentWithTraceback` | blast_gapalign.c | `protein_gapped_align()` | ⚠️ Check |
| Gapped cutoff (trigger) | blast_parameters.c | `ungap_cutoff = 22` (raw) | ❌ Should be 57 raw |

### 7. E-value / Statistics
| NCBI C Function | File | blast-rs Equivalent | Status |
|----------------|------|-------------------|--------|
| `BLAST_ComputeSearchSpace` | blast_stat.c | `compute_search_space()` | ✅ Matched |
| `BLAST_KarlinBlkCalc` | blast_stat.c | `lookup_protein_params()` | ✅ Matched |
| `BlastRawScore_ToEvalue` | blast_stat.c | `raw_to_evalue()` | ✅ Matched |
| Bit-to-raw conversion | blast_parameters.c | Missing for ungapped | ❌ Need to add |

## TODO

1. Fix `x_drop_ungapped`: convert 7 bits to 19 raw correctly
2. Fix `ungap_cutoff`: convert 22 bits to 57 raw  
3. Fix `s_BlastAaExtendTwoHit` match — verify left-must-reach-first logic
4. Add `diag_offset` trick to avoid memset per subject
5. Verify PV table size matches NCBI (should be 32768 for charsize=5)
6. Re-add srtA to sensitivity test once parameters are correct
7. Add bit-to-raw conversion for all x_dropoff parameters
