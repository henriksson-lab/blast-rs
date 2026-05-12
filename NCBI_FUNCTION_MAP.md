# NCBI BLAST+ C → blast-rs Function Map

Systematic comparison of every function in the NCBI blastp path against its blast-rs equivalent.
Each function must match the C original or the difference must be documented.

## Upstream Source Layout

For ongoing audit, prefer the local BLAST+ 2.17.0 tarball-shaped source tree,
because its paths and file versions match the CCC reference:

`ncbi-blast-2.17.0+-src/c++/src/algo/blast/core/...`

For local navigation, this repo keeps an ignored symlink:

`.tmp/upstream-ncbi-core-tarball -> ../ncbi-blast-2.17.0+-src/c++/src/algo/blast/core`

The checked-out NCBI C++ Toolkit is also available as project-adjacent source:

`/home/mahogny/github/claude/ncbi-cxx-toolkit-rs/ncbi-cxx-toolkit-public`

That checkout comes from:

`https://github.com/ncbi/ncbi-cxx-toolkit-public.git`

The GitHub toolkit checkout has the same core source tree without the BLAST+
tarball packaging prefix:

`src/algo/blast/core/...`

`.tmp/upstream-ncbi-core -> /home/mahogny/github/claude/ncbi-cxx-toolkit-rs/ncbi-cxx-toolkit-public/src/algo/blast/core`

Use the GitHub checkout only when we intentionally want current toolkit source
rather than the exact BLAST+ 2.17.0 tarball snapshot.

## Status: CORE PARITY MAP IN PROGRESS

## NCBI C Functions (in call order for blastp)

### 1. Lookup Table Construction
| NCBI C Function | File | blast-rs Equivalent | Status |
|----------------|------|-------------------|--------|
| `BlastAaLookupTableNew` | blast_aalookup.c | `ProteinLookupTable::build()` | ✅ Matched (thick backbone + overflow) |
| `BlastAaLookupFinalize` | blast_aalookup.c | (inline in build) | ✅ Matched |
| `_ComputeIndex` / `ComputeTableIndex` | blast_lookup.h | `word_hash()` | ✅ Matched (shift-based) |
| `ComputeTableIndexIncremental` | blast_lookup.h | Inline in scan loop | ✅ Matched |
| PV construction | blast_aalookup.c | Inline in build / `merge_pv()` no-hit prefilter | ✅ Matched (32768 bits for word size 3; batch prefilter preserves two-hit scanner results) |
| `LookupTableWrapInit` / `LookupTableWrapInit_MT` | lookup_wrap.c | `lookup_table_wrap_init()` / `lookup_table_wrap_init_mt()` | ✅ Matched Rust typed lookup ownership transfer with C-shaped represented variant dispatch; Rust consumes most construction inputs before wrapper installation, and preserves `num_threads` at the boundary for a future represented NaHash construction path |
| `BlastChooseNaExtend` | na_ungapped.c | `blast_choose_na_extend()` | ✅ Matched represented nucleotide lookup and extension callback selection as typed enums instead of C raw function pointers |

### 2. Subject Scanning
| NCBI C Function | File | blast-rs Equivalent | Status |
|----------------|------|-------------------|--------|
| `s_BlastAaScanSubject` | blast_aascan.c | (inline in scan loop) | ✅ Matched |
| `s_BlastnDiagHashExtendInitialHit` | na_ungapped.c | `s_blastn_diag_hash_extend_initial_hit()` | ✅ Matched represented packed-nucleotide core: exact-word expansion, diagonal explored-area guard, ungapped extension, and HSP-save signal; mask, two-hit-window, and lookup-layout routing remain in scanner callers |
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
| `MBSpaceNew` / `s_GetMBSpace` / `s_RefreshMBSpace` | greedy_align.c | `mbspace_new()` / `s_get_mbspace()` / `s_refresh_mbspace()` | ✅ Matched pool shape with Rust-owned chunk vector replacing C linked nodes |
| `s_BlastGreedyAlignMemAlloc` | blast_gapalign.c | `s_blast_greedy_align_mem_alloc()` | ✅ Matched allocation sizing for non-affine and affine greedy work buffers |
| `BLAST_GappedAlignmentWithTraceback` | blast_gapalign.c | `protein_gapped_align()` | ✅ Matched for ALIGN_EX script concatenation and terminal-gap accounting |
| `s_ChainingAlignment` | blast_gapalign.c | `s_chaining_alignment()` | ✅ Matched chaining DP recurrence and seed-drop threshold over Rust `ProteinHit` seeds |
| `BLAST_GetGappedScore` | blast_gapalign.c | `blast_get_gapped_score()` | ✅ Matched ordinary protein seed loop: score ordering, preliminary interval-tree containment, start-point selection, preliminary gapped DP cutoff, traceback rerun, common-endpoint purge, and score sort |
| `Blast_TracebackFromHSPList` | blast_traceback.c | `blast_traceback_from_hsp_list()` | ✅ Matched ordinary protein traceback branch: score order, containment tree, final-xdrop traceback rerun, HSP field update, common-endpoint purge, and score sort |
| `BLAST_GetUngappedHSPList` | blast_gapalign.c | `blast_get_ungapped_hsp_list()` | ✅ Matched init-hit loop, query-context rebasing, HSP field mapping, list allocation/reuse, and score sort |
| `JumperExtendRightWithTraceback` | jumper.c | `jumper_extend_right_with_traceback()` | ✅ Matched direct-byte right extension traceback, raw jumper edit ops, recent-error trace window, trimming, identity count, and ungapped-extension length |
| `BlastNaExtendJumper` | jumper.c | `blast_na_extend_jumper()` | ✅ Matched represented mapper/short-read extension loop: diagonal/query/subject ordering, exact pre-extension, explored diagonal skip, compressed Jumper traceback, `JumperGoodAlign`, mapper `map_info` attachment on created HSPs, HSP save, and optional small-word subject-index rescue |
| `JumperNaWordFinder` | na_ungapped.c | `jumper_na_word_finder()` / `jumper_na_word_finder_with_word_hits()` / `jumper_na_word_finder_with_subject_ranges()` | ✅ Matched represented contiguous lookup dispatch, direct extension branch, `MapperWordHits` batching, adjacent same-diagonal suppression, flush-on-full, final bucket flush, masked subject `s_DetermineScanningOffsets` scan ranges, discontiguous megablast template scan routing, and C-shaped stats updates; real masked/discontig mapper fixtures remain in TODO |
| `MB_IndexedWordFinder` | na_ungapped.c | `mb_indexed_word_finder()` / `mb_indexed_word_finder_with_fallback()` | ✅ Matched indexed callback ingestion, non-indexed OID fallback hook, `ungapped_extension` gating, IR diagonal suppression, ungapped extension filtering, and score sort for represented path; real indexed database fixtures remain in TODO |
| `ShortRead_IndexedWordFinder` | na_ungapped.c | `short_read_indexed_word_finder()` / `short_read_indexed_word_finder_with_fallback()` | ✅ Matched indexed callback ingestion, non-indexed OID fallback hook, IR diagonal suppression, Jumper traceback, good-align HSP save, and gapped stats for represented path; real indexed database fixtures remain in TODO |
| `DoAnchoredScan` | jumper.c | `do_anchored_scan()` | ✅ Matched represented anchored flank scan: big-word selection, ambiguity/polyA skip, packed subject word scan, best Jumper traceback HSP selection, and list save |
| `DoAnchoredSearch` | jumper.c | `do_anchored_search()` / `do_anchored_search_to_stream()` | ✅ Matched represented partial-chain flank recovery, original chain HSP copy on recovered flank, C-shaped stream write/status wrapper, and focused stream-output test; broader mapper fixtures remain in TODO |
| `FindPartialyCoveredQueries` | hspfilter_mapper.c | `find_partialy_covered_queries()` | ✅ Matched saved-chain scan and clone rules for same OID, score floor, leading uncovered query, and trailing uncovered query |
| `s_FindSpliceJunctionsForGap` / `s_FindSpliceJunctionsForOverlaps` | hspfilter_mapper.c | `s_find_splice_junctions_for_gap()` / `s_find_splice_junctions_for_gap_from_map_info()` / `s_find_splice_junctions_for_overlaps()` | ✅ Matched represented mapper splice-signal search, overhang-backed gap extension, overlap split trimming, and mapper splice-edge bit updates; broader upstream mapper fixtures remain in TODO |
| `s_TrimChainStartToSubjPos` / `s_TrimChainEndToSubjPos` | hspfilter_mapper.c | `s_trim_chain_start_to_subj_pos()` / `s_trim_chain_end_to_subj_pos()` | ✅ Matched represented chain trimming, score refresh, containment drop, and mapper splice/exon edge-bit clearing |
| `s_SetAdapter` / `s_SetPolyATail` | hspfilter_mapper.c | `s_set_adapter()` / `s_set_poly_a_tail()` | ✅ Matched represented chain trimming/position propagation and C-shaped mapper adapter/poly-A/exon edge-bit side effects |
| `s_MergeHSPs` | hspfilter_mapper.c | `s_merge_hsps()` | ✅ Matched represented extent/gap-script merge plus mapper edit append, right-edge transfer, and right-overhang replacement |
| `s_TrimHSP` | hspfilter_mapper.c | `s_trim_hsp()` | ✅ Matched represented gap-script/coordinate trim plus mapper edit pruning and subject-overhang reconstruction |
| `s_BlastHSPMapperSplicedPairedRun` / `s_FindRearrangedPairs` | hspfilter_mapper.c | `s_blast_hspmapper_spliced_paired_run()` / `s_find_rearranged_pairs()` | ✅ Matched represented paired-read grouping via `BlastContextInfo.segment_flags == eFirstSegment`, including non-even query positions; Rust stores pair indexes/cached scores instead of raw `HSPChain*` peer pointers |
| `s_PruneChains` | hspfilter_mapper.c | `s_prune_chains()` | ✅ Matched represented two-pass pruning, score multiplicity assignment, pair bonus handling, and `HSPChainFree`-style mate-link clearing when a partner chain is pruned |
| Gapped cutoff (trigger) | blast_parameters.c | per-query converted raw trigger | ✅ Matched |

### 6a. PHI-BLAST Pattern Search / Gapping
| NCBI C Function | File | blast-rs Equivalent | Status |
|----------------|------|-------------------|--------|
| `FindPatternHits` | pattern.c | `find_pattern_hits()` | ✅ Matched low-level PHI pattern scan helper |
| `PHIBlastScanSubject` | phi_lookup.c | `phi_blast_scan_subject()` | ✅ Matched scan shape; Rust returns typed `PhiInitialHit` subject bounds |
| `s_PHISaveInitialHit` | phi_extend.c | `s_phi_save_initial_hit()` | ✅ Matched C save order (`s_start`, `s_end`) via `PhiInitialHit` view |
| `PHIBlastWordFinder` | phi_extend.c | `phi_blast_word_finder()` | ✅ Matched scan/save/stats loop over subject pattern hits |
| `s_PHIGappedAlignment` | phi_gapalign.c | `s_phi_gapped_alignment()` | ✅ Matched score-only PHI gapping around the pattern hit |
| `PHIGetGappedScore` | phi_gapalign.c | `phi_get_gapped_score()` | ✅ Matched driver shape over explicit Rust PHI initial-hit view |
| `s_PHIBlastAlignPatterns` | phi_gapalign.c | `s_phi_blast_align_patterns()` | ✅ Matched PHI pattern edit-script construction using translated pattern-bound helpers and `s_banded_align()` |
| `PHIGappedAlignmentWithTraceback` | phi_gapalign.c | `phi_gapped_alignment_with_traceback()` | ✅ Matched wrapper flow; PHI flanks now use translated `ALIGN_EX`-style affine traceback with generic score matrix |
| `Blast_HSPUpdateWithTraceback` | blast_traceback.c | `blast_hsp_update_with_traceback()` | ✅ Matched ownership transfer of `edit_script` from gap workspace to HSP |
| `s_PHITracebackFromHSPList` | blast_traceback.c | `s_phi_traceback_from_hsp_list()` | ✅ Matched PHI list loop, traceback update/drop, score sort, PHI e-values, e-value reap, and PHI bit scores |

### 7. E-value / Statistics
| NCBI C Function | File | blast-rs Equivalent | Status |
|----------------|------|-------------------|--------|
| `BLAST_ComputeSearchSpace` | blast_stat.c | `compute_search_space()` | ✅ Matched |
| `BLAST_KarlinBlkCalc` | blast_stat.c | `lookup_protein_params()` | ✅ Matched |
| `BlastRawScore_ToEvalue` | blast_stat.c | `raw_to_evalue()` | ✅ Matched |
| `Blast_HSPListGetEvalues` | blast_hits.c | `blast_hsp_list_get_evalues()` | ✅ Matched represented C path: context Karlin selection, RPS preliminary fallback, gap decay, scaling factor, Spouge/Karlin branch, and `round_down`; simplified callers use `blast_hsp_list_get_evalues_simple()` |
| Bit-to-raw conversion | blast_parameters.c | implemented per program path | ✅ Matched |

### 8. Composition / PSI Kappa Helpers
| NCBI C Function | File | blast-rs Equivalent | Status |
|----------------|------|-------------------|--------|
| `s_GetStartFreqRatios` | blast_kappa.c | `get_start_freq_ratios()` | ✅ Matched |
| `s_GetPosBasedStartFreqRatios` | blast_kappa.c | `get_pos_based_start_freq_ratios()` | ✅ Matched |
| `s_ScalePosMatrix` | blast_kappa.c | `scale_pos_matrix()` / `psi_private_update_lambda_statistics()` | ✅ Matched boundary (position rows, private scaled PSI scores, and Lambda-ratio stat update hook) |
| `Kappa_posSearchItemsNew` / `Kappa_posSearchItemsFree` | blast_posit.c | `kappa_pos_search_items_new()` / `kappa_pos_search_items_free()` | ✅ Matched Rust-owned lifecycle and matrix/frequency payload |
| `Kappa_compactSearchItemsNew` / `Kappa_compactSearchItemsFree` | blast_posit.c | `kappa_compact_search_items_new()` / `kappa_compact_search_items_free()` | ✅ Matched Rust-owned compact search payload over `BlastScoreBlk` |
| `fillSfp` | blast_posit.c | `fill_sfp()` | ✅ Matched true-residue score-frequency construction |
| `impalaScaleMatrix` / `Kappa_impalaScaling` | blast_posit.c | `impala_scale_matrix()` / `kappa_impala_scaling()` | ✅ Matched PSI private/public matrix scaling and lambda recomputation boundary |
| `_PSIComputeScoreProbabilities` / `_PSIUpdateLambdaK` | blast_psi_priv.c | `psi_compute_score_probabilities()` / `psi_update_lambda_k()` | ✅ Matched PSSM score-frequency, ungapped Karlin update, and gapped PSI K derivation |
| `s_PSIDiscardIfUnused` / `_PSIPurgeAlignedRegion` | blast_psi_priv.c | `s_psi_discard_if_unused()` / `psi_purge_aligned_region()` | ✅ Matched aligned-region clearing and sequence-use discard rule |
| `s_fillColumnProbabilities` | blast_psi_priv.c | `s_fill_column_probabilities()` | ✅ Matched NCBIstdaa-to-effective-alphabet column probability ordering |
| `_PSISaveDiagnostics` / `_PSISaveCDDiagnostics` | blast_psi_priv.c | `psi_save_diagnostics()` / `psi_save_cd_diagnostics()` | ✅ Matched requested diagnostics copy-out and information-content calculation |
| `_PSIStructureGroupCustomization` | blast_psi_priv.c | `psi_structure_group_customization()` | ✅ Matched query-row discard and position-count refresh |
| `_PSIPurgeBiasedSegments` / `s_PSIPurgeSelfHits` / `s_PSIPurgeSimilarAlignments` | blast_psi_priv.c | `psi_purge_biased_segments()` / `s_psi_purge_self_hits()` / `s_psi_purge_similar_alignments()` | ✅ Matched packed-MSA purge FSM, self-hit pass, and near-identical pass |
| `_PSIComputeAlignmentBlocks` / `_PSIGetLeftExtents` / `_PSIGetRightExtents` / `_PSIComputeAlignedRegionLengths` | blast_psi_priv.c | `psi_compute_alignment_blocks()` / `psi_get_left_extents()` / `psi_get_right_extents()` / `psi_compute_aligned_region_lengths()` | ✅ Matched extent propagation and X-adjusted aligned-region lengths |
| `_PSICalculateNormalizedSequenceWeights` / `_PSICalculateMatchWeights` / `_PSISpreadGapWeights` | blast_psi_priv.c | `psi_calculate_normalized_sequence_weights()` / `psi_calculate_match_weights()` / `psi_spread_gap_weights()` | ✅ Matched Henikoff-style normalized weights, match-weight accumulation, and gap dispersion |
| `_PSICheckSequenceWeights` | blast_psi_priv.c | `psi_check_sequence_weights()` | ✅ Matched column sum validation, matching-sequence skip rule, query-X skip rule, and NSG compatibility threshold |
| `_PSIComputeFrequenciesFromCDs` / `s_PSIComputeFrequenciesFromCDsCleanup` | blast_psi_priv.c | `psi_compute_frequencies_from_cds()` / `s_psi_compute_frequencies_from_cds_cleanup()` | ✅ Matched CD weighted-frequency accumulation, query-residue inclusion, normalization, and observation cap |
| `_PSIScaleMatrix` | blast_psi_priv.c | `psi_scale_matrix()` | ✅ Matched scaled-PSSM application, Lambda/K update, bracketing, and binary-search scaling loop |
| `_PSIComputeFreqRatiosFromCDs` | blast_psi_priv.c | `psi_compute_freq_ratios_from_cds()` | ✅ Matched CD frequency-ratio blending with matrix ratios, pseudocount handling, and X-position zeroing |
| `s_PSISavePssm` | blast_psi.c | `s_psi_save_pssm()` | ✅ Matched internal PSSM matrix copy plus gapped/ungapped PSI Karlin field copy into `PSIMatrix` |
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
| `s_SmithWatermanScoreOnly` | blast_sw.c | `s_smith_waterman_score_only()` | ✅ Matched traceback-stage protein/PSSM score-only DP, including C's non-PSSM sequence swap |
| `s_GetTraceback` / `SmithWatermanScoreWithTraceback` | blast_sw.c | `s_get_traceback()` / `smith_waterman_score_with_traceback()` | ✅ Matched traceback byte layout, path bookkeeping, gap-open flags, non-PSSM swap restoration, PSSM no-swap behavior, and subject offset adjustment; Rust returns recovered hits instead of mutating `BlastHSPList` |
| `BLAST_SmithWatermanGetGappedScore` | blast_sw.c | `blast_smith_waterman_get_gapped_score()` | ✅ Matched score-only SW driver shape: ignores initial hits, scans valid contexts, applies RPS/per-context cutoffs, saves placeholder HSPs for traceback, and handles protein/PSSM plus packed nucleotide subjects |
| `s_NuclSmithWaterman` | blast_sw.c | `s_nucl_smith_waterman()` | ✅ Matched packed NCBI2na subject score-only DP |
| `Blast_RedoOneMatchSmithWaterman` | redo_alignment.c | `blast_redo_one_match_smith_waterman_with_callbacks_and_adjustment()` / `blast_redo_one_match_smith_waterman_with_callbacks_and_position_adjustment()` / `blast_redo_one_match_smith_waterman_with_callbacks()` / `blast_redo_one_match_smith_waterman_in_memory_nucl()` / `blast_redo_one_match_smith_waterman_in_memory_protein()` / `blast_redo_one_match_smith_waterman_in_memory_protein_with_adjustment()` / `blast_redo_one_match_smith_waterman_in_memory_protein_position_adjustment()` | ✅ Matched for callback no-composition, callback non-position composition-adjusted, callback position-based/PSSM, in-memory nucleotide/protein no-composition, materialized protein-space, BLASTX query-translated, TBLASTN translated-subject, TBLASTX query+subject-translated composition-adjusted, and position-based/PSSM score/start/final-X-drop paths; BLASTN/nucleotide composition-adjusted Smith-Waterman is intentionally rejected at the protein-space adjustment boundary |
| `Blast_RedoAlignmentCore_MT` materialized-subject path | blast_kappa.c | `blast_redo_alignment_core_mt_in_memory_subject()` / `blast_redo_alignment_core_mt_in_memory_subjects()` / `blast_redo_alignment_core_mt_seqsrc_subjects()` / `blast_redo_alignment_core_mt_with_callbacks()` / `merge_compo_thread_heaps()` | ✅ Matched for materialized, stream, SeqSrc, and callback-backed ordinary/SW redo branches currently represented in Rust, including per-context routing, heap retention/merge, early termination, program-specific SeqSrc encoding, supplied-context sum-stat/link-HSP evaluation, BLASTN composition-adjusted Smith-Waterman rejection, and PSI private Lambda-ratio stat update hook |
| `BlastCompo_EarlyTermination` | redo_alignment.c | `blast_compo_early_termination()` | ✅ Matched |
| `BlastCompo_AlignmentsFree` | redo_alignment.c | `alignments_free()` | ✅ Matched |
| `s_SequenceDataRelease` | redo_alignment.c | `sequence_data_release()` | ✅ Matched |
| `s_WindowInfoFree` | redo_alignment.c | `window_info_free()` | ✅ Matched |

## Status

CCC currently reports no partial/stub pairs in the tracked core set.

Represented but still incomplete audit areas:

| Upstream Function | Source | Current Rust Coverage | Remaining Gap |
|------------------|--------|-----------------------|--------|
| `s_RPSGapAlignDataPrepare` | blast_traceback.c:944 | Represented by `s_rps_gap_align_data_prepare`; builds concatenated profile contexts and owned PSSM/frequency row maps with C alphabet-size selection | Uses an owned Rust profile view rather than raw mmapped row pointers attached to `BlastScoreBlk::psi_matrix` |
| `s_RPSComputeTraceback` | blast_traceback.c:1058 | Represented by `s_rps_compute_traceback`; drains HSP stream, selects RPS profile rows, applies RPS Karlin K multiplier, routes traceback callback, restores RPS orientation, inserts and sorts results | SeqSrc fetching, rescaled per-profile PSSM allocation, CBS frequency-ratio setup, and direct `Blast_RedoAlignmentCore`/`Blast_TracebackFromHSPList` parity still need fixture integration |
| `PSICreatePssmWithDiagnostics` | blast_psi.c:105 | Represented by `psi_create_pssm_with_diagnostics`; orchestrates packed MSA, purge, internal allocation, optional structure-group customization, validation, alignment blocks, sequence weights, frequency ratios, PSSM conversion/scaling, save-out, and optional diagnostics save-out | Needs broader upstream parity fixtures for diagnostics payloads, structure-group mode, and hard error paths |
| `s_OutOfFrameAlignWithTraceback` | blast_gapalign.c:1334 | Represented by `s_out_of_frame_align_with_traceback`; preserves C-shaped inputs, offset outputs, prelim traceback block, OOF edit-op space, the three-frame OOF score-state recurrence, and represented C gap-open/gap-extension traceback flags | Needs upstream fixture parity for complex frame-shift traceback rendering |
| `s_OutOfFrameGappedAlign` | blast_gapalign.c:1950 | Represented by `s_out_of_frame_gapped_align`; traceback dispatch and score-only three-frame OOF DP score-state recurrence live directly in the name-matched Rust function | Needs TBLASTX/BLASTX fixture parity for final traceback rendering |
| `s_OutOfFrameSemiGappedAlignWrap` | blast_gapalign.c:4223 | Represented by `s_out_of_frame_semi_gapped_align_wrap`; BLASTX switch-sequence routing is present | Needs TBLASTX/BLASTX fixture parity |
| `s_BlastOOFTracebackToGapEditScript` | blast_gapalign.c:4451 | Represented by `s_blast_oof_traceback_to_gap_edit_script`; shifted-substitution propagation, left/right merge, in-frame-gap flow-through at the left/right boundary, nucleotide-length truncation, frame-shift splitting, and substitution post-adjustment are present | Needs fixture parity against C for end-to-end complex frame-shift scripts |

CCC currently reports no unmapped upstream functions in the tracked core set.
