/* Minimal test with a random query that won't match */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <algo/blast/core/blast_engine.h>
#include <algo/blast/core/blast_setup.h>
#include <algo/blast/core/blast_options.h>
#include <algo/blast/core/blast_util.h>
#include <algo/blast/core/blast_seqsrc.h>
#include <algo/blast/core/blast_seqsrc_impl.h>
#include <algo/blast/core/blast_hspstream.h>
#include <algo/blast/core/blast_diagnostics.h>
#include <algo/blast/core/hspfilter_collector.h>

/* Subject: 200 bases of packed ACGTACGT pattern (0x1B repeated) */
#define SUBJ_BASES 32
#define SUBJ_BYTES (SUBJ_BASES/4)
static Uint1 s_subj[SUBJ_BYTES];

static Int4 s_GetNumSeqs(void* d, void* a) { (void)d;(void)a; return 1; }
static Int4 s_GetMaxSeqLen(void* d, void* a) { (void)d;(void)a; return SUBJ_BASES; }
static Int4 s_GetMinSeqLen(void* d, void* a) { (void)d;(void)a; return SUBJ_BASES; }
static Int4 s_GetAvgSeqLen(void* d, void* a) { (void)d;(void)a; return SUBJ_BASES; }
static Int8 s_GetTotLen(void* d, void* a) { (void)d;(void)a; return SUBJ_BASES; }
static const char* s_GetName(void* d, void* a) { (void)d;(void)a; return "test"; }
static Boolean s_GetIsProt(void* d, void* a) { (void)d;(void)a; return FALSE; }
static Int4 s_GetSeqLen(void* d, void* a) { (void)d; if(!a) return -1; return (*(Int4*)a==0)?SUBJ_BASES:-1; }

static Int2 s_GetSequence(void* d, BlastSeqSrcGetSeqArg* arg) {
    (void)d;
    if (!arg || arg->oid != 0) return -1;
    fprintf(stderr, "[C2] GetSequence enc=%d\n", arg->encoding);
    if (!arg->seq) BlastSeqBlkNew(&arg->seq);
    /* Provide raw packed data directly */
    arg->seq->sequence = s_subj;
    arg->seq->sequence_start = NULL;
    arg->seq->length = SUBJ_BASES;
    arg->seq->sequence_allocated = FALSE; /* static data */
    arg->seq->oid = 0;
    return 0;
}
static void s_ReleaseSequence(void* d, BlastSeqSrcGetSeqArg* a) { (void)d;(void)a; }

static Int4 s_next_oid = 0;
static Int4 s_IterNext(void* d, BlastSeqSrcIterator* itr) {
    (void)d;(void)itr;
    if (s_next_oid >= 1) return BLAST_SEQSRC_EOF; /* -1 */
    return s_next_oid++;
}
static void s_ResetChunkIter(void* d) { (void)d; s_next_oid = 0; }
static void s_SetNumThreads(void* d, int n) { (void)d;(void)n; }

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
    switch(c) { case 'A': return 0; case 'C': return 1; case 'G': return 2; case 'T': return 3; default: return 14; }
}
static Uint1 compl_blastna(Uint1 b) {
    Uint1 t[]={3,2,1,0,5,4,7,6,8,9,13,12,11,10,14,15}; return (b<16)?t[b]:15;
}

int main(int argc, char** argv) {
    /* Subject: all 0x1B = ACGT repeated */
    memset(s_subj, 0x1B, SUBJ_BYTES);

    /* Query: random sequence that WON'T match ACGTACGT */
    const char* query_str = (argc > 1 && strcmp(argv[1], "match") == 0) ?
        "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT" :
        "AAGCCCAATAAACCACTCTGACTGGCCGAATAGGGATATAGGCAACGACAT";
    int qlen = strlen(query_str);

    int enc_len = 1 + qlen + 1 + qlen + 1;
    Uint1* encoded = (Uint1*)calloc(enc_len, 1);
    int pos = 0;
    encoded[pos++] = 15;
    for (int i=0; i<qlen; i++) encoded[pos++] = iupac_to_blastna(query_str[i]);
    encoded[pos++] = 15;
    for (int i=qlen-1; i>=0; i--) encoded[pos++] = compl_blastna(iupac_to_blastna(query_str[i]));
    encoded[pos++] = 15;

    BLAST_SequenceBlk* qblk = NULL;
    BlastSeqBlkNew(&qblk);
    qblk->sequence_start = encoded;
    qblk->sequence = encoded + 1;
    qblk->length = enc_len - 2;

    BlastQueryInfo* qi = BlastQueryInfoNew(eBlastTypeBlastn, 1);
    qi->max_length = qlen;
    qi->contexts[0].query_offset = 0;
    qi->contexts[0].query_length = qlen;
    qi->contexts[0].frame = 1;
    qi->contexts[0].is_valid = TRUE;
    qi->contexts[1].query_offset = qlen + 1;
    qi->contexts[1].query_length = qlen;
    qi->contexts[1].frame = -1;
    qi->contexts[1].is_valid = TRUE;

    BlastScoringOptions* so = NULL; BlastScoringOptionsNew(eBlastTypeBlastn, &so);
    so->reward = 2; so->penalty = -3; so->gap_open = 5; so->gap_extend = 2; so->gapped_calculation = TRUE;
    BlastInitialWordOptions* wo = NULL; BlastInitialWordOptionsNew(eBlastTypeBlastn, &wo);
    wo->window_size = 40;
    BlastExtensionOptions* eo = NULL; BlastExtensionOptionsNew(eBlastTypeBlastn, &eo, TRUE);
    eo->gap_x_dropoff = 30.0; eo->gap_x_dropoff_final = 100.0;
    eo->ePrelimGapExt = 1; /* eGreedyScoreOnly */
    eo->eTbackExt = 1; /* eGreedyTbck */
    BlastHitSavingOptions* ho = NULL; BlastHitSavingOptionsNew(eBlastTypeBlastn, &ho, TRUE);
    ho->expect_value = 10.0;
    BlastEffectiveLengthsOptions* elo = NULL; BlastEffectiveLengthsOptionsNew(&elo);
    BlastDatabaseOptions* dbo = NULL; BlastDatabaseOptionsNew(&dbo);
    QuerySetUpOptions* qo = NULL; BlastQuerySetUpOptionsNew(&qo);

    BlastScoreBlk* sbp = NULL;
    BlastSeqLoc* lseg = NULL; BlastMaskLoc* mask = NULL; Blast_Message* msg = NULL;
    BLAST_MainSetUp(eBlastTypeBlastn, qo, so, qblk, qi, 1.0, &lseg, &mask, &sbp, &msg, NULL);
    fprintf(stderr, "[C2] MainSetUp done, matrix[0][0]=%d\n", sbp->matrix->data[0][0]);

    LookupTableOptions* lo = NULL; LookupTableOptionsNew(eBlastTypeBlastn, &lo);
    BLAST_FillLookupTableOptions(lo, eBlastTypeBlastn, TRUE, 0.0, 28);
    fprintf(stderr, "[C2] lut_type=%d word_size=%d\n", lo->lut_type, lo->word_size);

    BlastSeqSrcNewInfo ni = { s_Constructor, NULL };
    BlastSeqSrc* ss = BlastSeqSrcNew(&ni);

    LookupTableWrap* lw = NULL;
    LookupTableWrapInit(qblk, lo, qo, lseg, sbp, &lw, NULL, &msg, ss);
    fprintf(stderr, "[C2] Lookup done, compressed=%p\n", qblk->compressed_nuc_seq);

    BlastHSPCollectorParams* cp = BlastHSPCollectorParamsNew(ho, 0, TRUE);
    BlastHSPWriterInfo* wi = BlastHSPCollectorInfoNew(cp);
    BlastHSPWriter* wr = BlastHSPWriterNew(&wi, qi, qblk);
    BlastHSPStream* hs = BlastHSPStreamNew(eBlastTypeBlastn, eo, FALSE, 1, wr);
    BlastDiagnostics* diag = Blast_DiagnosticsInit();

    BlastHSPResults* results = NULL;
    fprintf(stderr, "[C2] Calling Blast_RunFullSearch (query=%s)...\n",
        (argc > 1 && strcmp(argv[1], "match") == 0) ? "matching" : "random");

    Int4 rc = Blast_RunFullSearch(eBlastTypeBlastn, qblk, qi, ss, sbp, so, lw,
        wo, eo, ho, elo, NULL, dbo, hs, NULL, diag, &results, NULL, NULL);

    fprintf(stderr, "[C2] Done! rc=%d\n", rc);
    if (results) {
        fprintf(stderr, "[C2] queries=%d\n", results->num_queries);
        Blast_HSPResultsFree(results);
    }

    /* Cleanup (abbreviated) */
    Blast_DiagnosticsFree(diag);
    BlastSeqSrcFree(ss);
    free(encoded);
    return 0;
}
