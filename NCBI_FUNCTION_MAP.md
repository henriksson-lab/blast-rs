# NCBI BLAST+ C → blast-rs Function Map

Systematic comparison of every function in the NCBI blastp path against its blast-rs equivalent.
Each function must match the C original or the difference must be documented.

## Status: CORE PARITY MAP CURRENT

## NCBI C Functions (in call order for blastp)

### 1. Lookup Table Construction
| NCBI C Function | File | blast-rs Equivalent | Status |
|----------------|------|-------------------|--------|
| `BlastAaLookupTableNew` | blast_aalookup.c | `ProteinLookupTable::build()` | ✅ Matched (thick backbone + overflow) |
| `BlastAaLookupFinalize` | blast_aalookup.c | (inline in build) | ✅ Matched |
| `_ComputeIndex` / `ComputeTableIndex` | blast_lookup.h | `word_hash()` | ✅ Matched (shift-based) |
| `ComputeTableIndexIncremental` | blast_lookup.h | Inline in scan loop | ✅ Matched |
| PV construction | blast_aalookup.c | Inline in build / `merge_pv()` no-hit prefilter | ✅ Matched (32768 bits for word size 3; batch prefilter preserves two-hit scanner results) |

### 2. Subject Scanning
| NCBI C Function | File | blast-rs Equivalent | Status |
|----------------|------|-------------------|--------|
| `s_BlastAaScanSubject` | blast_aascan.c | (inline in scan loop) | ✅ Matched |
| `PV_TEST` macro | blast_extend.h | `pv[hash>>6] & (1<<(hash&63))` | ✅ Matched |

### 3. Two-Hit Word Finding
| NCBI C Function | File | blast-rs Equivalent | Status |
|----------------|------|-------------------|--------|
| `s_BlastAaWordFinder_TwoHit` | aa_ungapped.c:440 | `protein_scan_with_table_reuse()` | ✅ Matched |
| Diagonal array (`DiagStruct`) | blast_extend.h:57 | `diag_buf: Vec<(i32, bool)>` | ✅ Matched semantics |
| `diag_offset` trick | aa_ungapped.c:316 | `diag_offset = TWO_HIT_WINDOW` | ✅ Matched |
| Two-hit window check | aa_ungapped.c:373 | `diff >= TWO_HIT_WINDOW` | ✅ Matched |
| Word overlap check | aa_ungapped.c:384 | `diff < ws` | ✅ Matched |

### 4. Ungapped Extension
| NCBI C Function | File | blast-rs Equivalent | Status |
|----------------|------|-------------------|--------|
| `s_BlastAaExtendTwoHit` | aa_ungapped.c:1089 | `extend_left`/`extend_right` in scanner | ✅ Matched |
| Word-best-start scan | aa_ungapped.c:1108-1118 | Inline in scan loop | ✅ Matched |
| Left must reach first hit | aa_ungapped.c:1139 | `reached_first` check | ✅ Matched |
| Right extend with init_score | aa_ungapped.c:1147 | `extend_right(..., left_score)` | ✅ Matched |
| `s_BlastAaExtendLeft` | aa_ungapped.c:903 | `extend_left()` | ✅ Matched (`>=` xdrop) |
| `s_BlastAaExtendRight` | aa_ungapped.c:831 | `extend_right()` | ✅ Matched (`score <= 0`, `>=` xdrop, `s_last_off`) |
| `s_BlastAaExtendOneHit` | aa_ungapped.c:1071 | `protein_ungapped_extend()` | ✅ Matched for one-hit callers; two-hit path correctly uses `s_BlastAaExtendTwoHit` |

### 5. Constants & Parameters
| NCBI Constant | Value | blast-rs Value | Status |
|--------------|-------|---------------|--------|
| `BLAST_UNGAPPED_X_DROPOFF_PROT` | 7 bits | converted with ungapped lambda | ✅ Matched |
| `BLAST_GAP_X_DROPOFF_PROT` | 15 bits | 15 bits (converted to raw) | ✅ Matched |
| `BLAST_GAP_X_DROPOFF_FINAL_PROT` | 25 bits | 25 bits | ✅ Matched |
| `BLAST_GAP_TRIGGER_PROT` | 22 bits | converted with ungapped lambda/logK | ✅ Matched |
| `BLAST_WORD_THRESHOLD_BLASTP` | 11 | 11 | ✅ Matched |
| `AA_HITS_PER_CELL` | 3 | 3 | ✅ Matched |
| Two-hit window | 40 | 40 | ✅ Matched |
| charsize | 5 bits | 5 bits (CHARSIZE) | ✅ Matched |

### 6. Gapped Alignment
| NCBI C Function | File | blast-rs Equivalent | Status |
|----------------|------|-------------------|--------|
| `BlastGetStartForGappedAlignment` | blast_gapalign.c | `get_start_for_gapped_alignment()` | ✅ Matched |
| `BLAST_GappedAlignmentWithTraceback` | blast_gapalign.c | `protein_gapped_align()` | ✅ Matched for ALIGN_EX script concatenation and terminal-gap accounting |
| Gapped cutoff (trigger) | blast_parameters.c | per-query converted raw trigger | ✅ Matched |

### 7. E-value / Statistics
| NCBI C Function | File | blast-rs Equivalent | Status |
|----------------|------|-------------------|--------|
| `BLAST_ComputeSearchSpace` | blast_stat.c | `compute_search_space()` | ✅ Matched |
| `BLAST_KarlinBlkCalc` | blast_stat.c | `lookup_protein_params()` | ✅ Matched |
| `BlastRawScore_ToEvalue` | blast_stat.c | `raw_to_evalue()` | ✅ Matched |
| Bit-to-raw conversion | blast_parameters.c | implemented per program path | ✅ Matched |

### 8. Composition / PSI Kappa Helpers
| NCBI C Function | File | blast-rs Equivalent | Status |
|----------------|------|-------------------|--------|
| `s_GetStartFreqRatios` | blast_kappa.c | `get_start_freq_ratios()` | ✅ Matched |
| `s_GetPosBasedStartFreqRatios` | blast_kappa.c | `get_pos_based_start_freq_ratios()` | ✅ Matched |
| `s_ScalePosMatrix` | blast_kappa.c | `scale_pos_matrix()` / `psi_private_update_lambda_statistics()` | ✅ Matched boundary (position rows, private scaled PSI scores, and Lambda-ratio stat update hook) |
| `s_WithDistinctEnds` | redo_alignment.c | `with_distinct_ends()` | ✅ Matched |
| `s_WindowsFromProteinAligns` | redo_alignment.c | `windows_from_protein_aligns()` | ✅ Matched |
| `s_WindowsFromTranslatedAligns` | redo_alignment.c | `windows_from_translated_aligns()` | ✅ Matched |
| `s_WindowsFromAligns` | redo_alignment.c | `windows_from_aligns()` | ✅ Matched |
| `s_DistinctAlignmentsSort` | redo_alignment.c | `distinct_alignments_sort()` | ✅ Matched |
| `s_GetComposition` / `Blast_GetCompositionRange` | redo_alignment.c / composition_adjustment.c | `get_composition()` / `get_composition_range()` | ✅ Matched |
| `Blast_AdjustScores` scale/RE branches | composition_adjustment.c | `blast_adjust_scores_with_workspace()` / `blast_adjust_scores_scale_old_matrix()` / `blast_adjust_position_based_scores()` | ✅ Matched for square scale-old, RE optimization, p-value test, and PSSM scale-old |
| `s_EvalueFromScore` | redo_alignment.c | `evalue_from_score()` | ✅ Matched |
| `s_preliminaryTestNearIdentical` | redo_alignment.c | `preliminary_test_near_identical()` | ✅ Matched |
| `Blast_RedoAlignCallbacks` | redo_alignment.h | `BlastRedoAlignCallbacks` + callback type aliases | ✅ Matched API surface |
| `Blast_RedoOneMatch` | redo_alignment.c | `blast_redo_one_match_with_callbacks_and_adjustment()` / `blast_redo_one_match_with_callbacks_and_position_adjustment()` / `blast_redo_one_match_in_memory()` / `blast_redo_one_match_in_memory_with_adjustment()` / `blast_redo_one_match_in_memory_with_position_adjustment()` / wrappers | ✅ Matched for callback no-composition, callback non-position composition-adjusted, callback position-based/PSSM ordinary, in-memory no-composition BLASTN/BLASTP/BLASTX/TBLASTN/TBLASTX, in-memory non-position protein-space including BLASTX/TBLASTX translated protein-space, and position-based/PSSM ordinary redo |
| `Blast_RedoOneMatchSmithWaterman` significance gate | redo_alignment.c | `smith_waterman_alignment_is_significant()` | ✅ Matched |
| `Blast_RedoOneMatchSmithWaterman` | redo_alignment.c | `blast_redo_one_match_smith_waterman_with_callbacks_and_adjustment()` / `blast_redo_one_match_smith_waterman_with_callbacks_and_position_adjustment()` / `blast_redo_one_match_smith_waterman_with_callbacks()` / `blast_redo_one_match_smith_waterman_in_memory_nucl()` / `blast_redo_one_match_smith_waterman_in_memory_protein()` / `blast_redo_one_match_smith_waterman_in_memory_protein_with_adjustment()` / `blast_redo_one_match_smith_waterman_in_memory_protein_position_adjustment()` | ✅ Matched for callback no-composition, callback non-position composition-adjusted, callback position-based/PSSM, in-memory nucleotide/protein no-composition, materialized protein-space, BLASTX query-translated, TBLASTN translated-subject, TBLASTX query+subject-translated composition-adjusted, and position-based/PSSM score/start/final-X-drop paths; BLASTN/nucleotide composition-adjusted Smith-Waterman is intentionally rejected at the protein-space adjustment boundary |
| `Blast_RedoAlignmentCore_MT` materialized-subject path | blast_kappa.c | `blast_redo_alignment_core_mt_in_memory_subject()` / `blast_redo_alignment_core_mt_in_memory_subjects()` / `blast_redo_alignment_core_mt_seqsrc_subjects()` / `blast_redo_alignment_core_mt_with_callbacks()` / `merge_compo_thread_heaps()` | ✅ Matched for materialized, stream, SeqSrc, and callback-backed ordinary/SW redo branches currently represented in Rust, including per-context routing, heap retention/merge, early termination, program-specific SeqSrc encoding, supplied-context sum-stat/link-HSP evaluation, BLASTN composition-adjusted Smith-Waterman rejection, and PSI private Lambda-ratio stat update hook |
| `BlastCompo_EarlyTermination` | redo_alignment.c | `blast_compo_early_termination()` | ✅ Matched |
| `BlastCompo_AlignmentsFree` | redo_alignment.c | `alignments_free()` | ✅ Matched |
| `s_SequenceDataRelease` | redo_alignment.c | `sequence_data_release()` | ✅ Matched |
| `s_WindowInfoFree` | redo_alignment.c | `window_info_free()` | ✅ Matched |

## Status

No open parity gaps are tracked in this map; remaining work is broader fixture coverage and performance cleanup without changing translated behavior.
