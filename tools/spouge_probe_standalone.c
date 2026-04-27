/*
 * spouge_probe_standalone — re-implements BLAST_SpougeStoE / BLAST_SpougeEtoS
 * inline (verbatim from NCBI 2.17 blast_stat.c) and links against the
 * boost_erf.c ErfC implementation. No NCBI lib dependencies.
 *
 * Build:
 *   gcc -O2 -o spouge_probe_standalone \
 *     tools/spouge_probe_standalone.c \
 *     ncbi-blast-2.17.0+-src/c++/src/algo/blast/core/boost_erf.c \
 *     -I ncbi-blast-2.17.0+-src/c++/include \
 *     -I ncbi-blast-2.17.0+-src/c++/src/algo/blast/core \
 *     -lm
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef int    Int4;
typedef long   Int8;
typedef unsigned char Uint1;
typedef Uint1  Boolean;

/* From boost_erf.h */
extern double ErfC(double z);

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define ASSERT(x) ((void)0)

typedef struct Blast_KarlinBlk {
    double Lambda, K, logK, H, paramC;
} Blast_KarlinBlk;

typedef struct Blast_GumbelBlk {
    double Lambda, C, G;
    double a, Alpha, Sigma;
    double a_un, Alpha_un;
    double b, Beta, Tau;
    Int8 db_length;
    Boolean filled;
} Blast_GumbelBlk;

/* === Verbatim from NCBI 2.17 blast_stat.c:5176 === */
double
BLAST_SpougeStoE_local(Int4 y_,
                       Blast_KarlinBlk* kbp,
                       Blast_GumbelBlk* gbp,
                       Int4 m_, Int4 n_)
{
    double scale_factor = kbp->Lambda / gbp->Lambda;
    double db_scale_factor = (gbp->db_length) ?
            (double)gbp->db_length/(double)n_ : 1.0;

    double lambda_    = kbp->Lambda;
    double k_         = kbp->K;
    double ai_hat_    = gbp->a * scale_factor;
    double bi_hat_    = gbp->b;
    double alphai_hat_= gbp->Alpha * scale_factor;
    double betai_hat_ = gbp->Beta;
    double sigma_hat_ = gbp->Sigma * scale_factor;
    double tau_hat_   = gbp->Tau;

    double aj_hat_    = ai_hat_;
    double bj_hat_    = bi_hat_;
    double alphaj_hat_= alphai_hat_;
    double betaj_hat_ = betai_hat_;

    static double const_val = 0.39894228040143267793994605993438;

    double m_li_y, vi_y, sqrt_vi_y, m_F, P_m_F;
    double n_lj_y, vj_y, sqrt_vj_y, n_F, P_n_F;
    double c_y, p1, p2, area;
    double e_value;

    m_li_y = m_ - (ai_hat_*y_ + bi_hat_);
    vi_y = MAX(2.0*alphai_hat_/lambda_, alphai_hat_*y_+betai_hat_);
    sqrt_vi_y = sqrt(vi_y);
    m_F = m_li_y/sqrt_vi_y;
    P_m_F = ErfC(-m_F / sqrt(2.0)) / 2.0;
    p1 = m_li_y * P_m_F + sqrt_vi_y * const_val * exp(-0.5*m_F*m_F);

    n_lj_y = n_ - (aj_hat_*y_ + bj_hat_);
    vj_y = MAX(2.0*alphaj_hat_/lambda_, alphaj_hat_*y_+betaj_hat_);
    sqrt_vj_y = sqrt(vj_y);
    n_F = n_lj_y/sqrt_vj_y;
    P_n_F = ErfC(-n_F / sqrt(2.0)) / 2.0;
    p2 = n_lj_y * P_n_F + sqrt_vj_y * const_val * exp(-0.5*n_F*n_F);

    c_y = MAX(2.0*sigma_hat_/lambda_, sigma_hat_*y_+tau_hat_);
    area = p1 * p2 + c_y * P_m_F * P_n_F;

    e_value = area * k_ * exp(-lambda_ * y_) * db_scale_factor;
    ASSERT(e_value >= 0.0);

    return e_value;
}

/* === Verbatim from NCBI 2.17 blast_stat.c:5236 === */
Int4
BLAST_SpougeEtoS_local(double e0,
                       Blast_KarlinBlk* kbp,
                       Blast_GumbelBlk* gbp,
                       Int4 m, Int4 n)
{
    Int4 a=0, b, c;
    double e;
    double db_scale_factor = (gbp->db_length) ?
            (double)gbp->db_length : 1.0;

    b = MAX((int)(log(db_scale_factor/e0) / kbp->Lambda), 2);

    e = BLAST_SpougeStoE_local(b, kbp, gbp, m, n);

    if (e > e0) {
        while (e > e0) {
            a = b;
            b *= 2;
            e = BLAST_SpougeStoE_local(b, kbp, gbp, m, n);
        }
    } else {
        a = 0;
    }
    while (b-a > 1) {
        c = (a+b)/2;
        e = BLAST_SpougeStoE_local(c, kbp, gbp, m, n);
        if (e > e0) {
            a = c;
        } else {
            b = c;
        }
    }
    return a;
}

int main(int argc, char** argv) {
    Blast_KarlinBlk kbp = {0};
    kbp.Lambda = 0.267;
    kbp.K = 0.041;
    kbp.logK = -1.0;
    kbp.H = 0.14;

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

    if (argc >= 5) {
        qlen = atoi(argv[1]);
        slen = atoi(argv[2]);
        evalue = atof(argv[3]);
        gbp.db_length = (Int8) atol(argv[4]);
    }

    Int4 cutoff = BLAST_SpougeEtoS_local(evalue, &kbp, &gbp, qlen, slen);

    printf("inputs: qlen=%d slen=%d evalue=%g db_length=%ld\n",
           qlen, slen, evalue, (long) gbp.db_length);
    printf("kbp: Lambda=%.10g K=%.10g\n", kbp.Lambda, kbp.K);
    printf("gbp: Lambda=%.10g a=%.10g b=%.10g Alpha=%.10g Beta=%.10g Sigma=%.10g Tau=%.10g\n",
           gbp.Lambda, gbp.a, gbp.b, gbp.Alpha, gbp.Beta, gbp.Sigma, gbp.Tau);
    printf("BLAST_SpougeEtoS(E=%g) = %d\n", evalue, cutoff);

    for (Int4 s = 34; s <= 40; ++s) {
        double e = BLAST_SpougeStoE_local(s, &kbp, &gbp, qlen, slen);
        printf("BLAST_SpougeStoE(S=%d) = %.10g\n", s, e);
    }
    return 0;
}
