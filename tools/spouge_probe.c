/*
 * spouge_probe — calls NCBI 2.17 BLAST_SpougeStoE / BLAST_SpougeEtoS with the
 * exact inputs blast-rs uses for the seqp DB e=10 boundary case (NP_982592).
 *
 * Build: see Makefile.spouge_probe (links against the bundled libblast.a tree).
 *
 * The Spouge function declarations and required structs are drawn straight
 * from `algo/blast/core/blast_stat.h` so the binary calls into the very same
 * symbols the production NCBI binary uses.
 */
#include <stdio.h>
#include <stdlib.h>

/* Reproduce only the bits of NCBI's headers we need to declare the symbols.
 * Layout MUST match `blast_stat.h` exactly. */
typedef int   Int4;
typedef long  Int8;
typedef unsigned char Uint1;
typedef Uint1 Boolean;

typedef struct Blast_KarlinBlk {
    double Lambda;
    double K;
    double logK;
    double H;
    double paramC;
} Blast_KarlinBlk;

typedef struct Blast_GumbelBlk {
    double Lambda;
    double C;
    double G;
    double a;
    double Alpha;
    double Sigma;
    double a_un;
    double Alpha_un;
    double b;
    double Beta;
    double Tau;
    Int8 db_length;
    Boolean filled;
} Blast_GumbelBlk;

extern double BLAST_SpougeStoE(Int4 S, Blast_KarlinBlk* kbp, Blast_GumbelBlk* gbp, Int4 qlen, Int4 slen);
extern Int4   BLAST_SpougeEtoS(double E, Blast_KarlinBlk* kbp, Blast_GumbelBlk* gbp, Int4 qlen, Int4 slen);

int main(int argc, char** argv) {
    /* Karlin block: BLOSUM62 gap_open=11 gap_extend=1 gapped values.
     * From NCBI blast_stat.c blosum62_values row {11,1,...}: lambda=0.267, K=0.041, H=0.14 */
    Blast_KarlinBlk kbp = {0};
    kbp.Lambda = 0.267;
    kbp.K = 0.041;
    kbp.logK = -1.0; /* placeholder, not used by SpougeStoE/EtoS */
    kbp.H = 0.14;
    kbp.paramC = 0.0;

    /* Gumbel block: from Blast_GumbelBlkLoadFromTables for BLOSUM62 11/1.
     * Row: {11, 1, MAX, 0.267, 0.041, 0.14, 1.9, -30, 0.669720, 42.6028, 43.6362}.
     * Ungapped row [0]: {MAX, MAX, MAX, 0.3176, 0.134, 0.4012, 0.7916, -3.2, 0.623757, 4.96466, 4.96466}. */
    Blast_GumbelBlk gbp = {0};
    double a_gapped = 1.9;
    double Alpha_gapped = 42.6028;
    double Sigma_gapped = 43.6362;
    double a_un = 0.7916;
    double Alpha_un = 4.96466;
    double G = 11.0 + 1.0;

    gbp.Lambda = 0.267;
    gbp.C = 0.669720;
    gbp.G = G;
    gbp.a = a_gapped;
    gbp.Alpha = Alpha_gapped;
    gbp.Sigma = Sigma_gapped;
    gbp.a_un = a_un;
    gbp.Alpha_un = Alpha_un;
    gbp.b = 2.0 * G * (a_un - a_gapped);
    gbp.Beta = 2.0 * G * (Alpha_un - Alpha_gapped);
    gbp.Tau = 2.0 * G * (Alpha_un - Sigma_gapped);
    gbp.db_length = 566396L;
    gbp.filled = 1;

    Int4 qlen = 61;
    Int4 slen = 7;
    double evalue = 10.0;

    /* default override */
    if (argc >= 5) {
        qlen = atoi(argv[1]);
        slen = atoi(argv[2]);
        evalue = atof(argv[3]);
        gbp.db_length = (Int8) atol(argv[4]);
    }

    Int4 cutoff = BLAST_SpougeEtoS(evalue, &kbp, &gbp, qlen, slen);

    printf("inputs: qlen=%d slen=%d evalue=%g db_length=%ld\n",
           qlen, slen, evalue, (long) gbp.db_length);
    printf("kbp: Lambda=%.10g K=%.10g\n", kbp.Lambda, kbp.K);
    printf("gbp: Lambda=%.10g a=%.10g b=%.10g Alpha=%.10g Beta=%.10g Sigma=%.10g Tau=%.10g\n",
           gbp.Lambda, gbp.a, gbp.b, gbp.Alpha, gbp.Beta, gbp.Sigma, gbp.Tau);
    printf("BLAST_SpougeEtoS(E=%g) = %d\n", evalue, cutoff);

    /* Also dump SpougeStoE for scores 34..40 to see the boundary */
    for (Int4 s = 34; s <= 40; ++s) {
        double e = BLAST_SpougeStoE(s, &kbp, &gbp, qlen, slen);
        printf("BLAST_SpougeStoE(S=%d) = %.10g\n", s, e);
    }
    return 0;
}
