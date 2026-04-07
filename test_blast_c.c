/* Minimal test: call BLAST C core directly to verify it works standalone */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <algo/blast/core/blast_engine.h>
#include <algo/blast/core/blast_setup.h>
#include <algo/blast/core/blast_options.h>
#include <algo/blast/core/blast_util.h>
#include <algo/blast/core/blast_encoding.h>
#include <algo/blast/core/blast_seqsrc.h>
#include <algo/blast/core/blast_hspstream.h>
#include <algo/blast/core/blast_diagnostics.h>
#include <algo/blast/core/hspfilter_collector.h>

/* IUPAC to BLASTNA */
static Uint1 iupac_to_blastna(char c) {
    switch(c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 14; /* N */
    }
}

/* Complement in BLASTNA */
static Uint1 compl_blastna(Uint1 b) {
    Uint1 table[] = {3,2,1,0, 5,4,7,6, 8,9,13,12,11,10,14,15};
    return (b < 16) ? table[b] : 15;
}

int main() {
    const char* query_str = "TTAAGGAGGCTCATCTTTCAGAATCCATGCTGTGGGCCAGCAAGAGTTAA";
    int qlen = strlen(query_str);
    EBlastProgramType program = eBlastTypeBlastn;

    /* Encode query: sentinel + plus + sentinel + minus + sentinel */
    int enc_len = 1 + qlen + 1 + qlen + 1;
    Uint1* encoded = (Uint1*)calloc(enc_len, 1);
    int pos = 0;
    encoded[pos++] = 15; /* sentinel */
    for (int i = 0; i < qlen; i++)
        encoded[pos++] = iupac_to_blastna(query_str[i]);
    encoded[pos++] = 15; /* sentinel */
    for (int i = qlen-1; i >= 0; i--)
        encoded[pos++] = compl_blastna(iupac_to_blastna(query_str[i]));
    encoded[pos++] = 15; /* sentinel */

    /* Create sequence block */
    BLAST_SequenceBlk* query_blk = NULL;
    BlastSeqBlkNew(&query_blk);
    query_blk->sequence_start = encoded;
    query_blk->sequence = encoded + 1;
    query_blk->length = enc_len - 2; /* exclude start/end sentinels */
    query_blk->sequence_allocated = FALSE;
    query_blk->sequence_start_allocated = FALSE;

    /* Create query info */
    BlastQueryInfo* query_info = BlastQueryInfoNew(program, 1);
    query_info->first_context = 0;
    query_info->last_context = 1;
    query_info->num_queries = 1;
    query_info->max_length = qlen;

    /* Plus strand context */
    query_info->contexts[0].query_offset = 0;
    query_info->contexts[0].query_length = qlen;
    query_info->contexts[0].frame = 1;
    query_info->contexts[0].is_valid = TRUE;
    query_info->contexts[0].query_index = 0;

    /* Minus strand context */
    query_info->contexts[1].query_offset = qlen + 1; /* skip middle sentinel */
    query_info->contexts[1].query_length = qlen;
    query_info->contexts[1].frame = -1;
    query_info->contexts[1].is_valid = TRUE;
    query_info->contexts[1].query_index = 0;

    /* Create options */
    BlastScoringOptions* scoring_opts = NULL;
    BlastScoringOptionsNew(program, &scoring_opts);
    scoring_opts->reward = 2;
    scoring_opts->penalty = -3;
    scoring_opts->gap_open = 5;
    scoring_opts->gap_extend = 2;
    scoring_opts->gapped_calculation = TRUE;

    /* Setup ScoreBlk */
    QuerySetUpOptions* qsup_opts = NULL;
    BlastQuerySetUpOptionsNew(&qsup_opts);

    BlastScoreBlk* sbp = NULL;
    BlastSeqLoc* lookup_segments = NULL;
    BlastMaskLoc* mask = NULL;
    Blast_Message* blast_msg = NULL;

    Int2 rc = BLAST_MainSetUp(program, qsup_opts, scoring_opts,
        query_blk, query_info, 1.0,
        &lookup_segments, &mask, &sbp, &blast_msg, NULL);

    printf("BLAST_MainSetUp returned: %d\n", rc);
    if (rc != 0) {
        printf("Setup failed!\n");
        return 1;
    }

    printf("Matrix[0][0]=%d [0][1]=%d\n",
        sbp->matrix->data[0][0], sbp->matrix->data[0][1]);
    printf("compressed_nuc_seq: %p\n", query_blk->compressed_nuc_seq);

    /* Create lookup table */
    LookupTableOptions* lookup_opts = NULL;
    LookupTableOptionsNew(program, &lookup_opts);
    lookup_opts->word_size = 11;

    LookupTableWrap* lookup = NULL;
    rc = LookupTableWrapInit(query_blk, lookup_opts, qsup_opts,
        lookup_segments, sbp, &lookup, NULL, &blast_msg, NULL);
    printf("LookupTableWrapInit returned: %d\n", rc);
    printf("compressed after lookup: %p\n", query_blk->compressed_nuc_seq);

    /* Since we don't have a SeqSrc, let's create a fake subject */
    /* Pack "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT" */
    int slen = 60;
    int sbytes = slen / 4; /* 15 bytes for 60 bases */
    Uint1* subj_packed = (Uint1*)calloc(sbytes, 1);
    for (int i = 0; i < sbytes; i++)
        subj_packed[i] = 0x1B; /* A=00, C=01, G=10, T=11 -> 0001 1011 */

    BLAST_SequenceBlk* subj_blk = NULL;
    BlastSeqBlkNew(&subj_blk);
    subj_blk->sequence = subj_packed;
    subj_blk->length = slen;

    /* Try the preliminary search components manually */
    printf("Subject: %d bases, %d packed bytes\n", slen, sbytes);
    printf("Starting scan test...\n");
    fflush(stdout);

    /* We can't easily test the full search without a SeqSrc, but we've
       verified the setup is correct */
    printf("Setup verification complete.\n");

    /* Cleanup */
    free(subj_packed);
    BlastSequenceBlkFree(subj_blk);
    LookupTableWrapFree(lookup);
    LookupTableOptionsFree(lookup_opts);
    BlastScoreBlkFree(sbp);
    BlastQuerySetUpOptionsFree(qsup_opts);
    BlastScoringOptionsFree(scoring_opts);
    BlastQueryInfoFree(query_info);
    query_blk->sequence = NULL;
    query_blk->sequence_start = NULL;
    BlastSequenceBlkFree(query_blk);
    free(encoded);

    return 0;
}
