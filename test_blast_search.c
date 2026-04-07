/* Minimal test: run a BLAST search against a real database using our static library */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <algo/blast/core/blast_engine.h>
#include <algo/blast/core/blast_setup.h>
#include <algo/blast/core/blast_options.h>
#include <algo/blast/core/blast_util.h>
#include <algo/blast/core/blast_encoding.h>
#include <algo/blast/core/blast_seqsrc.h>
#include <algo/blast/core/blast_seqsrc_impl.h>
#include <algo/blast/core/blast_hspstream.h>
#include <algo/blast/core/blast_diagnostics.h>
#include <algo/blast/core/hspfilter_collector.h>

/* Simple SeqSrc that provides one hardcoded packed sequence */
static Uint1 s_subject_packed[] = {
    0x17, 0x51, 0xC9, 0xFE, 0xFE, 0xCB, 0x8E, 0x77, 0xB2,
    0x44, 0x9C, 0x01, 0x39, 0x09, 0x00  /* 60 bases from f555 */
};
static Int4 s_subject_len = 56; /* 14*4=56, last byte has 0 remainder */

static Int4 s_GetNumSeqs(void* d, void* a) { (void)d; (void)a; return 1; }
static Int4 s_GetMaxSeqLen(void* d, void* a) { (void)d; (void)a; return s_subject_len; }
static Int4 s_GetMinSeqLen(void* d, void* a) { (void)d; (void)a; return s_subject_len; }
static Int4 s_GetAvgSeqLen(void* d, void* a) { (void)d; (void)a; return s_subject_len; }
static Int8 s_GetTotLen(void* d, void* a) { (void)d; (void)a; return s_subject_len; }
static const char* s_GetName(void* d, void* a) { (void)d; (void)a; return "test"; }
static Boolean s_GetIsProt(void* d, void* a) { (void)d; (void)a; return FALSE; }

static Int4 s_GetSeqLen(void* d, void* a) {
    (void)d;
    if (!a) return -1;
    Int4 oid = *(Int4*)a;
    if (oid != 0) return -1;
    return s_subject_len;
}

static Int2 s_GetSequence(void* d, BlastSeqSrcGetSeqArg* arg) {
    (void)d;
    if (!arg || arg->oid != 0) return -1;

    fprintf(stderr, "[C] GetSequence called, encoding=%d\n", arg->encoding);

    /* Provide raw packed data, no sentinels */
    int buf_len = sizeof(s_subject_packed);
    Uint1* buf = (Uint1*)malloc(buf_len);
    memcpy(buf, s_subject_packed, buf_len);

    if (!arg->seq) {
        BlastSeqBlkNew(&arg->seq);
    }
    arg->seq->sequence = buf;
    arg->seq->sequence_start = NULL;
    arg->seq->length = s_subject_len;
    arg->seq->sequence_allocated = TRUE;
    arg->seq->oid = 0;
    return 0;
}

static void s_ReleaseSequence(void* d, BlastSeqSrcGetSeqArg* arg) {
    (void)d;
    if (arg && arg->seq) {
        BlastSequenceBlkClean(arg->seq);
    }
}

static Int4 s_next_oid = 0;
static Int4 s_IterNext(void* d, BlastSeqSrcIterator* itr) {
    (void)d;
    if (s_next_oid >= 1) return 0x7FFFFFFF; /* EOF */
    return s_next_oid++;
}

static void s_ResetChunkIter(void* d) { (void)d; s_next_oid = 0; }
static void s_SetNumThreads(void* d, int n) { (void)d; (void)n; }

static BlastSeqSrc* s_Constructor(BlastSeqSrc* ss, void* arg) {
    (void)arg;
    _BlastSeqSrcImpl_SetDeleteFnPtr(ss, NULL);
    _BlastSeqSrcImpl_SetDataStructure(ss, NULL);
    _BlastSeqSrcImpl_SetGetNumSeqs(ss, s_GetNumSeqs);
    _BlastSeqSrcImpl_SetGetNumSeqsStats(ss, s_GetNumSeqs);
    _BlastSeqSrcImpl_SetGetMaxSeqLen(ss, s_GetMaxSeqLen);
    _BlastSeqSrcImpl_SetGetMinSeqLen(ss, s_GetMinSeqLen);
    _BlastSeqSrcImpl_SetGetAvgSeqLen(ss, s_GetAvgSeqLen);
    _BlastSeqSrcImpl_SetGetTotLen(ss, s_GetTotLen);
    _BlastSeqSrcImpl_SetGetTotLenStats(ss, s_GetTotLen);
    _BlastSeqSrcImpl_SetGetName(ss, s_GetName);
    _BlastSeqSrcImpl_SetGetIsProt(ss, s_GetIsProt);
    _BlastSeqSrcImpl_SetGetSequence(ss, s_GetSequence);
    _BlastSeqSrcImpl_SetGetSeqLen(ss, s_GetSeqLen);
    _BlastSeqSrcImpl_SetReleaseSequence(ss, s_ReleaseSequence);
    _BlastSeqSrcImpl_SetIterNext(ss, s_IterNext);
    _BlastSeqSrcImpl_SetResetChunkIterator(ss, s_ResetChunkIter);
    _BlastSeqSrcImpl_SetSetNumberOfThreads(ss, s_SetNumThreads);
    return ss;
}

static Uint1 iupac_to_blastna(char c) {
    switch(c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 14;
    }
}
static Uint1 compl_blastna(Uint1 b) {
    Uint1 t[] = {3,2,1,0, 5,4,7,6, 8,9,13,12,11,10,14,15};
    return (b < 16) ? t[b] : 15;
}

int main() {
    const char* query_str = "TTAAGGAGGCTCATCTTTCAGAATCCATGCTGTGGGCCAGCAAGAGTTAA";
    int qlen = strlen(query_str);
    EBlastProgramType program = eBlastTypeBlastn;

    /* Encode query */
    int enc_len = 1 + qlen + 1 + qlen + 1;
    Uint1* encoded = (Uint1*)calloc(enc_len, 1);
    int pos = 0;
    encoded[pos++] = 15;
    for (int i = 0; i < qlen; i++) encoded[pos++] = iupac_to_blastna(query_str[i]);
    encoded[pos++] = 15;
    for (int i = qlen-1; i >= 0; i--) encoded[pos++] = compl_blastna(iupac_to_blastna(query_str[i]));
    encoded[pos++] = 15;

    /* Create sequence block */
    BLAST_SequenceBlk* query_blk = NULL;
    BlastSeqBlkNew(&query_blk);
    query_blk->sequence_start = encoded;
    query_blk->sequence = encoded + 1;
    query_blk->length = enc_len - 2;

    /* Create query info */
    BlastQueryInfo* query_info = BlastQueryInfoNew(program, 1);
    query_info->max_length = qlen;
    query_info->contexts[0].query_offset = 0;
    query_info->contexts[0].query_length = qlen;
    query_info->contexts[0].frame = 1;
    query_info->contexts[0].is_valid = TRUE;
    query_info->contexts[1].query_offset = qlen + 1;
    query_info->contexts[1].query_length = qlen;
    query_info->contexts[1].frame = -1;
    query_info->contexts[1].is_valid = TRUE;

    /* Create options */
    BlastScoringOptions* scoring_opts = NULL;
    BlastScoringOptionsNew(program, &scoring_opts);
    scoring_opts->reward = 2;
    scoring_opts->penalty = -3;
    scoring_opts->gap_open = 5;
    scoring_opts->gap_extend = 2;
    scoring_opts->gapped_calculation = TRUE;

    BlastInitialWordOptions* word_opts = NULL;
    BlastInitialWordOptionsNew(program, &word_opts);
    word_opts->window_size = 40;

    BlastExtensionOptions* ext_opts = NULL;
    BlastExtensionOptionsNew(program, &ext_opts, TRUE);
    ext_opts->gap_x_dropoff = 30.0;
    ext_opts->gap_x_dropoff_final = 100.0;

    BlastHitSavingOptions* hit_opts = NULL;
    BlastHitSavingOptionsNew(program, &hit_opts, TRUE);
    hit_opts->expect_value = 10.0;

    BlastEffectiveLengthsOptions* eff_len_opts = NULL;
    BlastEffectiveLengthsOptionsNew(&eff_len_opts);

    BlastDatabaseOptions* db_opts = NULL;
    BlastDatabaseOptionsNew(&db_opts);

    QuerySetUpOptions* qsup_opts = NULL;
    BlastQuerySetUpOptionsNew(&qsup_opts);

    /* Setup ScoreBlk */
    BlastScoreBlk* sbp = NULL;
    BlastSeqLoc* lookup_segments = NULL;
    BlastMaskLoc* mask = NULL;
    Blast_Message* blast_msg = NULL;

    fprintf(stderr, "[C] Calling BLAST_MainSetUp...\n");
    Int2 rc = BLAST_MainSetUp(program, qsup_opts, scoring_opts,
        query_blk, query_info, 1.0,
        &lookup_segments, &mask, &sbp, &blast_msg, NULL);
    fprintf(stderr, "[C] BLAST_MainSetUp returned: %d\n", rc);
    if (rc != 0) { fprintf(stderr, "Failed!\n"); return 1; }

    /* Create lookup table */
    LookupTableOptions* lookup_opts = NULL;
    LookupTableOptionsNew(program, &lookup_opts);
    lookup_opts->word_size = 28;

    /* Create SeqSrc */
    BlastSeqSrcNewInfo new_info;
    new_info.constructor = s_Constructor;
    new_info.ctor_argument = NULL;
    BlastSeqSrc* seq_src = BlastSeqSrcNew(&new_info);

    LookupTableWrap* lookup = NULL;
    rc = LookupTableWrapInit(query_blk, lookup_opts, qsup_opts,
        lookup_segments, sbp, &lookup, NULL, &blast_msg, seq_src);
    fprintf(stderr, "[C] LookupTableWrapInit returned: %d\n", rc);

    /* Create HSP stream */
    BlastHSPCollectorParams* coll_params = BlastHSPCollectorParamsNew(hit_opts, 0, TRUE);
    BlastHSPWriterInfo* writer_info = BlastHSPCollectorInfoNew(coll_params);
    BlastHSPWriter* writer = BlastHSPWriterNew(&writer_info, query_info, query_blk);
    BlastHSPStream* hsp_stream = BlastHSPStreamNew(program, ext_opts, FALSE, 1, writer);

    BlastDiagnostics* diagnostics = Blast_DiagnosticsInit();

    /* Run the full search */
    BlastHSPResults* results = NULL;
    fprintf(stderr, "[C] Calling Blast_RunFullSearch...\n");
    rc = Blast_RunFullSearch(program, query_blk, query_info, seq_src, sbp,
        scoring_opts, lookup, word_opts, ext_opts, hit_opts, eff_len_opts,
        NULL, db_opts, hsp_stream, NULL, diagnostics, &results, NULL, NULL);
    fprintf(stderr, "[C] Blast_RunFullSearch returned: %d\n", rc);

    if (results) {
        fprintf(stderr, "[C] num_queries=%d\n", results->num_queries);
        for (int q = 0; q < results->num_queries; q++) {
            if (results->hitlist_array[q]) {
                fprintf(stderr, "[C] query %d: %d hsp lists\n", q,
                    results->hitlist_array[q]->hsplist_count);
            }
        }
        Blast_HSPResultsFree(results);
    }

    /* Cleanup */
    Blast_DiagnosticsFree(diagnostics);
    BlastSeqSrcFree(seq_src);
    LookupTableWrapFree(lookup);
    LookupTableOptionsFree(lookup_opts);
    BlastScoreBlkFree(sbp);
    BlastQuerySetUpOptionsFree(qsup_opts);
    BlastDatabaseOptionsFree(db_opts);
    BlastEffectiveLengthsOptionsFree(eff_len_opts);
    BlastHitSavingOptionsFree(hit_opts);
    BlastExtensionOptionsFree(ext_opts);
    BlastInitialWordOptionsFree(word_opts);
    BlastScoringOptionsFree(scoring_opts);
    BlastQueryInfoFree(query_info);
    query_blk->sequence = NULL;
    query_blk->sequence_start = NULL;
    BlastSequenceBlkFree(query_blk);
    free(encoded);

    fprintf(stderr, "[C] Done!\n");
    return 0;
}
