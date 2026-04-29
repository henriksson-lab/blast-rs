//! Rust equivalent of blast_stat.c — Karlin-Altschul statistics.
//! This is the mathematical core for computing E-values and bit scores.

/// Karlin-Altschul statistical parameters for one context.
#[derive(Debug, Clone, Default)]
pub struct KarlinBlk {
    pub lambda: f64, // Lambda parameter
    pub k: f64,      // K parameter
    pub log_k: f64,  // ln(K)
    pub h: f64,      // H (relative entropy)
    /// NCBI `BlastScoreBlk::round_down` (`blast_stat.c:1868`). When
    /// `true`, odd scores are rounded down to even before e-value and
    /// bit-score computation. Set by `nucl_gapped_kbp_lookup` for
    /// scoring systems whose table values are only valid at even
    /// scores (e.g. `(2,-3)`, `(2,-5)`, `(2,-7)`, `(3,-4)`).
    pub round_down: bool,
}

/// Apply NCBI's `score &= ~1` even-rounding when `round_down` is set.
#[inline]
fn round_down_score(score: i32, round_down: bool) -> i32 {
    if round_down {
        score & !1
    } else {
        score
    }
}

impl KarlinBlk {
    pub fn is_valid(&self) -> bool {
        self.lambda > 0.0 && self.k > 0.0 && self.h > 0.0
    }

    /// Convert a raw score to a bit score.
    /// Port of NCBI `Blast_HSPListGetBitScores` (`blast_hits.c:1927`):
    /// `(score * Lambda - logK) / NCBIMATH_LN2`. Note that unlike the
    /// E-value path, the bit-score formula does NOT apply `round_down`
    /// even-masking — `Blast_HSPListGetEvalues` (`blast_hits.c:1864-1869`)
    /// applies `score &= ~1`, but `Blast_HSPListGetBitScores` (`:1907`)
    /// only has a commented-out `#if 0` assertion for it.
    pub fn raw_to_bit(&self, raw_score: i32) -> f64 {
        (self.lambda * raw_score as f64 - self.log_k) / crate::math::NCBIMATH_LN2
    }

    /// Convert a raw score and search space to an E-value.
    /// Port of NCBI `BLAST_KarlinStoE_simple` (`blast_stat.c:4157`):
    /// returns `-1.0` when the Karlin block is degenerate, otherwise
    /// `searchsp * exp(-lambda * S + logK)` (NCBI `:4170` uses `logK`
    /// inside the exponent for numerical stability at large magnitudes).
    /// Applies `round_down` before the formula (`blast_hits.c:1864-1869`).
    pub fn raw_to_evalue(&self, raw_score: i32, search_space: f64) -> f64 {
        if self.lambda < 0.0 || self.k < 0.0 || self.h < 0.0 {
            return -1.0;
        }
        let score = round_down_score(raw_score, self.round_down);
        search_space * (-self.lambda * score as f64 + self.log_k).exp()
    }

    /// Convert an E-value to the minimum raw score needed.
    /// Port of NCBI `BlastKarlinEtoS_simple` (`blast_stat.c:4040`): the
    /// e-value is clamped to `K_SMALL_FLOAT` before conversion to avoid
    /// floating-point exceptions on extremely tight cutoffs. Returns
    /// `BLAST_SCORE_MIN` when lambda/K/H are invalid (NCBI `:4054-4057`).
    pub fn evalue_to_raw(&self, evalue: f64, search_space: f64) -> i32 {
        // NCBI checks `Lambda < 0. || K < 0. || H < 0.0`. We don't carry
        // `h` here at the call site, but zero/negative lambda or K are
        // definitive invalid markers; match NCBI's return value.
        if self.lambda < 0.0 || self.k < 0.0 || self.h < 0.0 {
            return BLAST_SCORE_MIN;
        }
        let denom = search_space * self.k;
        if denom <= 0.0 || self.lambda <= 0.0 {
            return BLAST_SCORE_MIN;
        }
        let evalue = evalue.max(K_SMALL_FLOAT);
        let score = -(evalue / denom).ln() / self.lambda;
        score.ceil() as i32
    }
}

/// `BLAST_GAP_PROB` (`blast_parameters.h:66`): gap probability for
/// ungapped sum-statistics (0.5).
pub const BLAST_GAP_PROB: f64 = 0.5;
/// `BLAST_GAP_PROB_GAPPED` (`blast_parameters.h:67`): gap probability
/// for gapped sum-statistics (1.0).
pub const BLAST_GAP_PROB_GAPPED: f64 = 1.0;
/// Port of NCBI `BLAST_GAP_DECAY_RATE` (`blast_parameters.h:68`):
/// default gap decay rate for ungapped search.
pub const BLAST_GAP_DECAY_RATE: f64 = 0.5;
/// Port of NCBI `BLAST_GAP_DECAY_RATE_GAPPED` (`blast_parameters.h:69`):
/// default gap decay rate for gapped search.
pub const BLAST_GAP_DECAY_RATE_GAPPED: f64 = 0.1;
/// `BLAST_GAP_SIZE` (`blast_parameters.h:70`): default gap size (nt
/// distance) used by HSP linking.
pub const BLAST_GAP_SIZE: i32 = 40;
/// `BLAST_OVERLAP_SIZE` (`blast_parameters.h:71`): default overlap
/// allowed between linked HSPs.
pub const BLAST_OVERLAP_SIZE: i32 = 9;
/// `RESTRICTED_ALIGNMENT_WORST_EVALUE` (`blast_parameters.h:193`):
/// composition-adjusted restricted alignment cutoff.
pub const RESTRICTED_ALIGNMENT_WORST_EVALUE: f64 = 10.0;

/// Port of NCBI `CUTOFF_E_BLASTN` (`blast_parameters.h:76`): default
/// e-value used by the preliminary ungapped cutoff for nucleotide
/// searches.
pub const CUTOFF_E_BLASTN: f64 = 0.05;
/// Port of NCBI `CUTOFF_E_BLASTP` (`blast_parameters.h:77`).
pub const CUTOFF_E_BLASTP: f64 = 1.0e-300;
/// Port of NCBI `CUTOFF_E_BLASTX` (`blast_parameters.h:78`).
pub const CUTOFF_E_BLASTX: f64 = 1.0;
/// Port of NCBI `CUTOFF_E_TBLASTN` (`blast_parameters.h:79`).
pub const CUTOFF_E_TBLASTN: f64 = 1.0;
/// Port of NCBI `CUTOFF_E_TBLASTX` (`blast_parameters.h:80`).
pub const CUTOFF_E_TBLASTX: f64 = 1.0e-300;

// X-dropoff defaults (in bits). NCBI `blast_options.h:122-148`.
/// `BLAST_UNGAPPED_X_DROPOFF_PROT` — ungapped dropoff for protein.
pub const BLAST_UNGAPPED_X_DROPOFF_PROT: i32 = 7;
/// `BLAST_UNGAPPED_X_DROPOFF_NUCL` — ungapped dropoff for nucleotide.
pub const BLAST_UNGAPPED_X_DROPOFF_NUCL: i32 = 20;
/// `BLAST_GAP_X_DROPOFF_PROT` — default preliminary-gapped dropoff for protein.
pub const BLAST_GAP_X_DROPOFF_PROT: i32 = 15;
/// `BLAST_GAP_X_DROPOFF_NUCL` — default preliminary-gapped dropoff for nucleotide.
pub const BLAST_GAP_X_DROPOFF_NUCL: i32 = 30;
/// `BLAST_GAP_X_DROPOFF_GREEDY` — default dropoff for greedy megablast.
pub const BLAST_GAP_X_DROPOFF_GREEDY: i32 = 25;
/// `BLAST_GAP_X_DROPOFF_FINAL_PROT` — final-gapped dropoff for protein.
pub const BLAST_GAP_X_DROPOFF_FINAL_PROT: i32 = 25;
/// `BLAST_GAP_X_DROPOFF_FINAL_NUCL` — final-gapped dropoff for nucleotide.
pub const BLAST_GAP_X_DROPOFF_FINAL_NUCL: i32 = 100;
/// `BLAST_GAP_TRIGGER_PROT` — protein gap-trigger bit threshold (`blast_options.h:137`).
pub const BLAST_GAP_TRIGGER_PROT: f64 = 22.0;
/// `BLAST_GAP_TRIGGER_NUCL` — nucleotide gap-trigger bit threshold (`blast_options.h:140`).
pub const BLAST_GAP_TRIGGER_NUCL: f64 = 27.0;

// Default word sizes (NCBI `blast_options.h:66-73`).
/// `BLAST_WORDSIZE_PROT` — default word size for protein searches (3).
pub const BLAST_WORDSIZE_PROT: i32 = 3;
/// `BLAST_WORDSIZE_NUCL` — default word size for blastn (11).
pub const BLAST_WORDSIZE_NUCL: i32 = 11;
/// `BLAST_WORDSIZE_MEGABLAST` — default word size for contiguous megablast (28).
pub const BLAST_WORDSIZE_MEGABLAST: i32 = 28;
/// `BLAST_WORDSIZE_MAPPER` — default word size for the magicblast mapper (18).
pub const BLAST_WORDSIZE_MAPPER: i32 = 18;

// Default gap costs (NCBI `blast_options.h:84-98`).
/// `BLAST_GAP_OPEN_PROT` — protein gap-open (11).
pub const BLAST_GAP_OPEN_PROT: i32 = 11;
/// `BLAST_GAP_OPEN_NUCL` — blastn gap-open (5).
pub const BLAST_GAP_OPEN_NUCL: i32 = 5;
/// `BLAST_GAP_OPEN_MEGABLAST` — megablast gap-open (0).
pub const BLAST_GAP_OPEN_MEGABLAST: i32 = 0;
/// `BLAST_GAP_EXTN_PROT` — protein gap-extend (1).
pub const BLAST_GAP_EXTN_PROT: i32 = 1;
/// `BLAST_GAP_EXTN_NUCL` — blastn gap-extend (2).
pub const BLAST_GAP_EXTN_NUCL: i32 = 2;
/// `BLAST_GAP_EXTN_MEGABLAST` — megablast gap-extend (0).
pub const BLAST_GAP_EXTN_MEGABLAST: i32 = 0;

// Default match/mismatch scores (NCBI `blast_options.h:151-152`).
/// `BLAST_PENALTY` — default nucleotide mismatch score (-3).
pub const BLAST_PENALTY: i32 = -3;
/// `BLAST_REWARD` — default nucleotide match score (1).
pub const BLAST_REWARD: i32 = 1;

// Default neighboring-word thresholds (NCBI `blast_options.h:104-116`).
/// `BLAST_WORD_THRESHOLD_BLASTP` (11).
pub const BLAST_WORD_THRESHOLD_BLASTP: f64 = 11.0;
/// `BLAST_WORD_THRESHOLD_BLASTP_FAST` — word-size 5 threshold for blastp/x-fast (20).
pub const BLAST_WORD_THRESHOLD_BLASTP_FAST: f64 = 20.0;
/// `BLAST_WORD_THRESHOLD_BLASTP_WD_SZ_6` (21).
pub const BLAST_WORD_THRESHOLD_BLASTP_WD_SZ_6: f64 = 21.0;
/// `BLAST_WORD_THRESHOLD_BLASTP_WD_SZ_7` (20.25).
pub const BLAST_WORD_THRESHOLD_BLASTP_WD_SZ_7: f64 = 20.25;
/// `BLAST_WORD_THRESHOLD_BLASTN` (0 = no threshold).
pub const BLAST_WORD_THRESHOLD_BLASTN: f64 = 0.0;
/// `BLAST_WORD_THRESHOLD_BLASTX` (12).
pub const BLAST_WORD_THRESHOLD_BLASTX: f64 = 12.0;
/// `BLAST_WORD_THRESHOLD_TBLASTN` (13).
pub const BLAST_WORD_THRESHOLD_TBLASTN: f64 = 13.0;
/// `BLAST_WORD_THRESHOLD_TBLASTX` (13).
pub const BLAST_WORD_THRESHOLD_TBLASTX: f64 = 13.0;
/// `BLAST_WORD_THRESHOLD_MEGABLAST` (0).
pub const BLAST_WORD_THRESHOLD_MEGABLAST: f64 = 0.0;

/// `BLAST_SCAN_RANGE_NUCL` — default scan range for blastn (0, `blast_options.h:63`).
pub const BLAST_SCAN_RANGE_NUCL: i32 = 0;
/// `BLAST_GAP_X_DROPOFF_TBLASTX` (0, `blast_options.h:134`).
pub const BLAST_GAP_X_DROPOFF_TBLASTX: i32 = 0;
/// `BLAST_GAP_X_DROPOFF_FINAL_TBLASTX` (0, `blast_options.h:148`).
pub const BLAST_GAP_X_DROPOFF_FINAL_TBLASTX: i32 = 0;

/// `BLAST_DEFAULT_MATRIX` — default scoring matrix for protein searches
/// (`blast_options.h:77`).
pub const BLAST_DEFAULT_MATRIX: &str = "BLOSUM62";

/// `BLAST_PENALTY_MAPPER` (`blast_options.h:154`): mapper nucleotide mismatch.
pub const BLAST_PENALTY_MAPPER: i32 = -4;
/// `BLAST_REWARD_MAPPER` (`blast_options.h:155`): mapper nucleotide match.
pub const BLAST_REWARD_MAPPER: i32 = 1;
/// `BLAST_GAP_OPEN_MAPPER` (`blast_options.h:89`).
pub const BLAST_GAP_OPEN_MAPPER: i32 = 0;
/// `BLAST_GAP_EXTN_MAPPER` (`blast_options.h:98`).
pub const BLAST_GAP_EXTN_MAPPER: i32 = 4;

/// `PSI_INCLUSION_ETHRESH` (`blast_options.h:163`): PSI-BLAST inclusion
/// threshold (e-value below which a hit is used for the next iteration).
pub const PSI_INCLUSION_ETHRESH: f64 = 0.002;
/// `PSI_PSEUDO_COUNT_CONST` (`blast_options.h:164`): PSI-BLAST
/// pseudo-count constant (0 disables).
pub const PSI_PSEUDO_COUNT_CONST: i32 = 0;
/// `DELTA_INCLUSION_ETHRESH` (`blast_options.h:165`): DELTA-BLAST inclusion.
pub const DELTA_INCLUSION_ETHRESH: f64 = 0.05;

/// `BLAST_GENETIC_CODE` (`blast_options.h:168`): the standard genetic
/// code (NCBI table #1).
pub const BLAST_GENETIC_CODE: i32 = 1;
/// `MAX_DB_WORD_COUNT_MAPPER` (`blast_options.h:174`): word cap for
/// the magicblast mapper.
pub const MAX_DB_WORD_COUNT_MAPPER: i32 = 30;

/// `HSP_MAX_WINDOW` — sliding-window size used by
/// `BlastGetOffsetsForGappedAlignment` / `BlastGetStartForGappedAlignment`
/// to pick a high-scoring seed offset near an ungapped HSP
/// (`blast_gapalign_priv.h:120`).
pub const HSP_MAX_WINDOW: usize = 11;

/// NCBI `MININT` (`blast_gapalign.c:58`): `INT4_MIN/2`. Used as the
/// sentinel for impossible DP cells; halving avoids underflow when
/// small mismatches are added to this value during the recurrence.
pub const MININT: i32 = i32::MIN / 2;

/// `BLAST_EXPECT_VALUE` — default e-value threshold
/// (`blast_options.h:158`).
pub const BLAST_EXPECT_VALUE: f64 = 10.0;
/// `BLAST_HITLIST_SIZE` — default number of database sequences to save hits
/// for (`blast_options.h:160`).
pub const BLAST_HITLIST_SIZE: usize = 500;

// Default two-hit window sizes (NCBI `blast_options.h:57-61`).
/// `BLAST_WINDOW_SIZE_PROT` — protein two-hit window (40).
pub const BLAST_WINDOW_SIZE_PROT: i32 = 40;
/// `BLAST_WINDOW_SIZE_NUCL` — blastn (non-MB) window (0 = disabled).
pub const BLAST_WINDOW_SIZE_NUCL: i32 = 0;
/// `BLAST_WINDOW_SIZE_MEGABLAST` — megablast window (0 = disabled).
pub const BLAST_WINDOW_SIZE_MEGABLAST: i32 = 0;
/// `BLAST_WINDOW_SIZE_DISC` — discontiguous-megablast window (40).
pub const BLAST_WINDOW_SIZE_DISC: i32 = 40;

/// Port of NCBI `BLAST_GapDecayDivisor` (`blast_stat.c:4079`).
/// Computes the divisor used by sum-statistics to compensate for the
/// effect of choosing the best among multiple alignments:
/// `(1 - decayrate) * decayrate^(nsegs - 1)`. Typical `decayrate` values
/// are [`BLAST_GAP_DECAY_RATE_GAPPED`] (0.1) and [`BLAST_GAP_DECAY_RATE`] (0.5).
pub fn gap_decay_divisor(decayrate: f64, nsegs: u32) -> f64 {
    // NCBI: `return (1. - decayrate) * BLAST_Powi(decayrate, nsegs - 1);`.
    (1.0 - decayrate) * crate::math::powi(decayrate, (nsegs as i32) - 1)
}

/// Port of NCBI `BLAST_Cutoffs` (`blast_stat.c:4089`). Given a desired
/// e-value `e_in`, a Karlin-Altschul block, a search space size, and an
/// optional decay-rate adjustment, returns a tuple
/// `(cutoff_score, effective_evalue)`.
///
/// * `dodecay=true, 0 < gap_decay_rate < 1` scales the input e-value by
///   `BLAST_GapDecayDivisor(gap_decay_rate, 1)` before converting to a
///   raw score (tightens the cutoff), then divides the recomputed
///   e-value back out on the return trip.
/// * The final cutoff is `max(s_floor, EtoS(e))`; callers pass a floor
///   (typically 1) to match NCBI's behavior where a user-specified
///   minimum score wins if it's larger than the statistics-derived one.
///
/// Returns `(1, e_in)` when the KarlinBlk is uncomputed (lambda/K/H == -1.0),
/// matching NCBI `BLAST_Cutoffs` (`blast_stat.c:4101-4102`) which uses the
/// exact -1.0 sentinel — not `< 0.0` — to detect "computation failed".
/// Note: this differs from `raw_to_evalue` and `evalue_to_raw`, which use
/// `< 0.0` to match their respective NCBI verbatim checks (NCBI itself is
/// inconsistent across these three functions; the Rust port preserves
/// each function's individual NCBI convention).
pub fn blast_cutoffs(
    s_floor: i32,
    e_in: f64,
    kbp: &KarlinBlk,
    searchsp: f64,
    dodecay: bool,
    gap_decay_rate: f64,
) -> (i32, f64) {
    if kbp.lambda == -1.0 || kbp.k == -1.0 || kbp.h == -1.0 {
        return (1, e_in);
    }
    let esave = e_in;
    let mut s = s_floor;
    let mut e = e_in;
    let mut s_changed = false;
    let mut es = 1i32;
    if e > 0.0 {
        if dodecay && gap_decay_rate > 0.0 && gap_decay_rate < 1.0 {
            e *= gap_decay_divisor(gap_decay_rate, 1);
        }
        es = kbp.evalue_to_raw(e, searchsp);
    }
    if es > s {
        s_changed = true;
        s = es;
    }
    // Recompute the e-value from the final cutoff when the input e was
    // non-positive or the cutoff didn't change. Mirrors NCBI's
    // `blast_stat.c:4134-4146`.
    let e_out = if esave <= 0.0 || !s_changed {
        let mut recomputed = kbp.raw_to_evalue(s, searchsp);
        if dodecay && gap_decay_rate > 0.0 && gap_decay_rate < 1.0 {
            recomputed /= gap_decay_divisor(gap_decay_rate, 1);
        }
        recomputed
    } else {
        esave
    };
    (s, e_out)
}

// NCBI s_BlastSumP interpolation tables (blast_stat.c:4359-4379).
// Retained verbatim so any drift is easy to spot against the C source.

#[rustfmt::skip]
const SUM_P_TAB2: [f64; 19] = [
    0.01669,  0.0249,   0.03683,  0.05390,  0.07794,  0.1111,   0.1559,   0.2146,
    0.2890,   0.3794,   0.4836,   0.5965,   0.7092,   0.8114,   0.8931,   0.9490,
    0.9806,   0.9944,   0.9989,
];

#[rustfmt::skip]
const SUM_P_TAB3: [f64; 38] = [
    0.9806,   0.9944,   0.9989,   0.0001682,0.0002542,0.0003829,0.0005745,0.0008587,
    0.001278, 0.001893, 0.002789, 0.004088, 0.005958, 0.008627, 0.01240,  0.01770,
    0.02505,  0.03514,  0.04880,  0.06704,  0.09103,  0.1220,   0.1612,   0.2097,
    0.2682,   0.3368,   0.4145,   0.4994,   0.5881,   0.6765,   0.7596,   0.8326,
    0.8922,   0.9367,   0.9667,   0.9846,   0.9939,   0.9980,
];

#[rustfmt::skip]
const SUM_P_TAB4: [f64; 55] = [
    2.658e-07,4.064e-07,6.203e-07,9.450e-07,1.437e-06,2.181e-06,3.302e-06,4.990e-06,
    7.524e-06,1.132e-05,1.698e-05,2.541e-05,3.791e-05,5.641e-05,8.368e-05,0.0001237,
    0.0001823,0.0002677,0.0003915,0.0005704,0.0008275,0.001195, 0.001718, 0.002457,
    0.003494, 0.004942, 0.006948, 0.009702, 0.01346,  0.01853,  0.02532,  0.03431,
    0.04607,  0.06128,  0.08068,  0.1051,   0.1352,   0.1719,   0.2157,   0.2669,
    0.3254,   0.3906,   0.4612,   0.5355,   0.6110,   0.6849,   0.7544,   0.8168,
    0.8699,   0.9127,   0.9451,   0.9679,   0.9827,   0.9915,   0.9963,
];

/// Port of NCBI `s_BlastSumP` (`blast_stat.c:4357`).
/// Estimates the Sum P-value by calculation or interpolation. Accuracy:
/// ~2.5 digits throughout the range of `r` (number of segments) and `s`
/// (total score in nats, adjusted by `-r*log(K*N)`).
///
/// For `r = 0` returns `0.0`; for `r = 1` uses the closed-form
/// `1 - exp(-exp(-s))`; for `r` in `2..=4` uses the table-interpolation
/// branches at `blast_stat.c:4394-4404` (either the analytic tail or
/// the `kTable` interpolation); for `r >= 5` delegates to `sum_p_calc`
/// (NCBI `s_BlastSumPCalc`, Romberg integration). Always returns `Some`.
pub fn sum_p(r: u32, s: f64) -> Option<f64> {
    if r == 1 {
        // NCBI `blast_stat.c:4384`: `return -BLAST_Expm1(-exp(-s));`.
        return Some(-crate::math::expm1(-(-s).exp()));
    }
    if r == 0 {
        return Some(0.0);
    }
    if (2..=4).contains(&r) {
        let r1 = (r - 1) as i32;
        let r_i = r as i32;
        if s >= (r * r + r - 1) as f64 {
            // NCBI `blast_stat.c:4394`: `a = BLAST_LnGammaInt(r+1)`.
            let a = crate::math::ln_gamma_int(r_i + 1);
            return Some(r_i as f64 * (r1 as f64 * s.ln() - s - a - a).exp());
        }
        if s > -2.0 * r as f64 {
            // Table interpolation — verbatim NCBI `blast_stat.c:4397-4403`:
            //   i = (Int4)(a = s+s+(4*r));
            //   a -= i;
            //   i = kTabsize[r2] - i;
            //   return a*kTable[r2][i-1] + (1.-a)*kTable[r2][i];
            let table: &[f64] = match r {
                2 => &SUM_P_TAB2,
                3 => &SUM_P_TAB3,
                4 => &SUM_P_TAB4,
                _ => unreachable!(),
            };
            let mut a = s + s + (4 * r) as f64;
            let mut i = a as i32;
            a -= i as f64;
            let tab_last = (table.len() - 1) as i32; // == NCBI `kTabsize[r2]`.
            i = tab_last - i;
            return Some(a * table[(i - 1) as usize] + (1.0 - a) * table[i as usize]);
        }
        return Some(1.0);
    }
    // r >= 5: delegate to `sum_p_calc` (Romberg integration).
    Some(sum_p_calc(r, s))
}

/// Port of NCBI `s_BlastSumPCalc` (`blast_stat.c:4269`).
/// Computes the Sum P-value via double Romberg integration for
/// `r ≥ 5` (callers with smaller `r` should go through `sum_p` instead).
/// Matches the Karlin-Altschul PNAS 1993 formula with the iteratively
/// tightened `itmin` that NCBI uses when the convergence is marginal.
pub fn sum_p_calc(r: u32, s: f64) -> f64 {
    let r_i = r as i32;
    if r == 1 {
        if s > 8.0 {
            return (-s).exp();
        }
        // NCBI `blast_stat.c:4271`: `return -BLAST_Expm1(-exp(-s));`.
        return -crate::math::expm1(-(-s).exp());
    }
    if r < 1 {
        return 0.0;
    }

    // Early-out bounds where the integral is essentially 1 ("no
    // significant distinction"). NCBI `blast_stat.c:4286-4305` uses
    // nested `if/else if { if (s <= …) return 1.0; }` — keep the nested
    // layout so the port is line-diffable against the C source.
    let rs = r as f64;
    #[allow(clippy::collapsible_if, clippy::collapsible_else_if)]
    if r < 8 {
        if s <= -2.3 * rs {
            return 1.0;
        }
    } else if r < 15 {
        if s <= -2.5 * rs {
            return 1.0;
        }
    } else if r < 27 {
        if s <= -3.0 * rs {
            return 1.0;
        }
    } else if r < 51 {
        if s <= -3.4 * rs {
            return 1.0;
        }
    } else if r < 101 {
        if s <= -4.0 * rs {
            return 1.0;
        }
    }

    let stddev = rs.sqrt();
    let stddev4 = 4.0 * stddev;
    let r1 = r - 1;

    if r > 100 {
        let est_mean = -(r as f64) * r1 as f64;
        if s <= est_mean - stddev4 {
            return 1.0;
        }
    }

    let logr = rs.ln();
    let mean = rs * (1.0 - logr) - 0.5;
    if s <= mean - stddev4 {
        return 1.0;
    }

    let (t0, itmin0) = if s >= mean {
        (s + 6.0 * stddev, 1)
    } else {
        (mean + 6.0 * stddev, 2)
    };

    let num_hsps = r_i;
    let num_hsps_minus_2 = r_i - 2;
    // NCBI `blast_stat.c:4338`: `adj1 = num_hsps_minus_2*logr
    //                                   - BLAST_LnGammaInt(r1)
    //                                   - BLAST_LnGammaInt(r)`.
    let adj1 = num_hsps_minus_2 as f64 * logr
        - crate::math::ln_gamma_int(r1 as i32)
        - crate::math::ln_gamma_int(r_i);
    /// NCBI `kSumpEpsilon` (`blast_stat.c:4276`): Romberg convergence
    /// epsilon for the sum-P double integral.
    const EPSILON: f64 = 0.002;

    // Inner integrand: the callback nested inside the outer Romberg.
    // Mirrors `s_OuterIntegralCback` + `s_InnerIntegralCback`.
    let integrand = |sv: f64| {
        let adj2 = adj1 - sv;
        let sdvir = sv / num_hsps as f64;
        let mx = if sv > 0.0 { sdvir + 3.0 } else { 3.0 };
        crate::math::romberg_integrate(
            |x| {
                let y = (x - sdvir).exp();
                if y == f64::INFINITY {
                    return 0.0;
                }
                if num_hsps_minus_2 == 0 {
                    return (adj2 - y).exp();
                }
                if x == 0.0 {
                    return 0.0;
                }
                (num_hsps_minus_2 as f64 * x.ln() + adj2 - y).exp()
            },
            0.0,
            mx,
            EPSILON,
            0,
            1,
        )
    };

    // NCBI iteratively tightens `itmin` if the outer integral returns a
    // marginal value (blast_stat.c:4345).
    let mut itmin = itmin0;
    let mut d;
    loop {
        d = crate::math::romberg_integrate(integrand, s, t0, EPSILON, 0, itmin);
        if d == f64::INFINITY {
            return d;
        }
        if !(s < mean && d < 0.4 && itmin < 4) {
            break;
        }
        itmin += 1;
    }
    d.min(1.0)
}

/// Port of NCBI `BLAST_KarlinPtoE` (`blast_stat.c:4175`).
/// Convert a P-value to an E-value. For `p = 1` returns `f64::INFINITY`
/// (NCBI returns `INT4_MAX`, which the sum-E callers cap below).
pub fn karlin_p_to_e(p: f64) -> f64 {
    if !(0.0..=1.0).contains(&p) {
        return i32::MIN as f64;
    }
    if p == 1.0 {
        return f64::INFINITY;
    }
    // NCBI: `return -BLAST_Log1p(-p)`.
    -crate::math::log1p(-p)
}

const SUM_E_CAP: f64 = i32::MAX as f64;

/// Port of NCBI `BLAST_SmallGapSumE` (`blast_stat.c:4418`).
/// Computes the e-value of a collection of distinct alignments
/// separated by small gaps. Matches NCBI's formula and cap at
/// `INT4_MAX`. Delegates the P-value step to `sum_p`, which handles
/// every `num` via closed-form / table interpolation / Romberg
/// (`sum_p_calc`). Returns `Some` on every finite input; the `Option`
/// return is retained for API symmetry with neighbouring helpers.
#[allow(clippy::too_many_arguments)]
pub fn small_gap_sum_e(
    starting_points: i32,
    num: u32,
    mut xsum: f64,
    query_length: i32,
    subject_length: i32,
    searchsp_eff: f64,
    weight_divisor: f64,
) -> Option<f64> {
    let mut sum_e = if num == 1 {
        searchsp_eff * (-xsum).exp()
    } else {
        let pair_search_space = subject_length as f64 * query_length as f64;
        xsum -= pair_search_space.ln() + 2.0 * (num - 1) as f64 * (starting_points as f64).ln();
        xsum -= crate::math::ln_factorial(num as i32);
        let p = sum_p(num, xsum)?;
        karlin_p_to_e(p) * (searchsp_eff / pair_search_space)
    };
    if weight_divisor == 0.0 {
        sum_e = SUM_E_CAP;
    } else {
        sum_e /= weight_divisor;
        if sum_e > SUM_E_CAP {
            sum_e = SUM_E_CAP;
        }
    }
    Some(sum_e)
}

/// Port of NCBI `BLAST_UnevenGapSumE` (`blast_stat.c:4491`).
/// Used for HSP collections with asymmetric gap widths — e.g. exons
/// separated by introns in translated searches.
#[allow(clippy::too_many_arguments)]
pub fn uneven_gap_sum_e(
    query_start_points: i32,
    subject_start_points: i32,
    num: u32,
    mut xsum: f64,
    query_length: i32,
    subject_length: i32,
    searchsp_eff: f64,
    weight_divisor: f64,
) -> Option<f64> {
    let mut sum_e = if num == 1 {
        searchsp_eff * (-xsum).exp()
    } else {
        let pair_search_space = subject_length as f64 * query_length as f64;
        xsum -= pair_search_space.ln()
            + (num - 1) as f64
                * ((query_start_points as f64).ln() + (subject_start_points as f64).ln());
        xsum -= crate::math::ln_factorial(num as i32);
        let p = sum_p(num, xsum)?;
        karlin_p_to_e(p) * (searchsp_eff / pair_search_space)
    };
    if weight_divisor == 0.0 {
        sum_e = SUM_E_CAP;
    } else {
        sum_e /= weight_divisor;
        if sum_e > SUM_E_CAP {
            sum_e = SUM_E_CAP;
        }
    }
    Some(sum_e)
}

/// Port of NCBI `BLAST_LargeGapSumE` (`blast_stat.c:4532`).
/// Computes the e-value of a collection of distinct alignments
/// separated by arbitrarily large gaps.
pub fn large_gap_sum_e(
    num: u32,
    mut xsum: f64,
    query_length: i32,
    subject_length: i32,
    searchsp_eff: f64,
    weight_divisor: f64,
) -> Option<f64> {
    let mut sum_e = if num == 1 {
        searchsp_eff * (-xsum).exp()
    } else {
        let q = query_length as f64;
        let s = subject_length as f64;
        xsum -= num as f64 * (s * q).ln() - crate::math::ln_factorial(num as i32);
        let p = sum_p(num, xsum)?;
        karlin_p_to_e(p) * (searchsp_eff / (q * s))
    };
    if weight_divisor == 0.0 {
        sum_e = SUM_E_CAP;
    } else {
        sum_e /= weight_divisor;
        if sum_e > SUM_E_CAP {
            sum_e = SUM_E_CAP;
        }
    }
    Some(sum_e)
}

/// Score frequency distribution.
#[derive(Debug, Clone)]
pub struct ScoreFreq {
    pub score_min: i32,
    pub score_max: i32,
    pub obs_min: i32,
    pub obs_max: i32,
    pub score_avg: f64,
    pub sprob: Vec<f64>, // probability for each score value
}

/// Gapped Karlin-Altschul parameters (precomputed tables).
#[derive(Debug, Clone)]
pub struct GappedParams {
    pub gap_open: i32,
    pub gap_extend: i32,
    pub lambda: f64,
    pub k: f64,
    pub h: f64,
    pub alpha: f64,
    pub beta: f64,
}

/// Gumbel block parameters for Spouge finite-size correction (FSC) e-value.
/// Port of NCBI Blast_GumbelBlk from blast_stat.h.
#[derive(Debug, Clone)]
pub struct GumbelBlk {
    pub lambda: f64,    // gbp->Lambda (gapped lambda from table)
    pub a: f64,         // gbp->a (alpha from array[6])
    pub b: f64,         // gbp->b = 2*G*(a_un - a)
    pub alpha: f64,     // gbp->Alpha (alpha_v from array[9])
    pub beta: f64,      // gbp->Beta = 2*G*(Alpha_un - Alpha)
    pub sigma: f64,     // gbp->Sigma (from array[10])
    pub tau: f64,       // gbp->Tau = 2*G*(Alpha_un - Sigma)
    pub db_length: i64, // total database length
}

#[derive(Debug, Clone, Copy)]
struct MatrixStatRow {
    gap_open: i32,
    gap_extend: i32,
    lambda: f64,
    k: f64,
    h: f64,
    alpha: f64,
    beta: f64,
    _theta: f64,
    alpha_v: f64,
    sigma: f64,
}

impl MatrixStatRow {
    const fn new(
        gap_open: i32,
        gap_extend: i32,
        lambda: f64,
        k: f64,
        h: f64,
        alpha: f64,
        beta: f64,
        _theta: f64,
        alpha_v: f64,
        sigma: f64,
    ) -> Self {
        Self {
            gap_open,
            gap_extend,
            lambda,
            k,
            h,
            alpha,
            beta,
            _theta,
            alpha_v,
            sigma,
        }
    }

    fn gapped_params(self) -> GappedParams {
        GappedParams {
            gap_open: self.gap_open,
            gap_extend: self.gap_extend,
            lambda: self.lambda,
            k: self.k,
            h: self.h,
            alpha: self.alpha,
            beta: self.beta,
        }
    }
}

/// Compute Spouge finite-size correction e-value.
/// Port of BLAST_SpougeStoE from blast_stat.c:5176.
/// Uses per-subject lengths for more accurate e-values than simple Karlin formula.
/// Spouge e-value with NCBI's HSP-linking gap-decay correction
/// applied. Used by translated searches (blastx/tblastn/tblastx) and
/// ungapped paths where `do_sum_stats=TRUE` (set in
/// `blast_options.c:1467`). NCBI applies this divisor in
/// `link_hsps.c:1791` via `link_hsp_params->gap_decay_rate`. The
/// default rate is `BLAST_GAP_DECAY_RATE_GAPPED = 0.1` for protein
/// gapped, giving an effective `/0.9 ≈ ×1.111` scaling vs. the raw
/// per-pair Spouge.
pub fn spouge_evalue_with_gap_decay(
    score: i32,
    kbp: &KarlinBlk,
    gbp: &GumbelBlk,
    query_length: i32,
    subject_length: i32,
) -> f64 {
    let raw = spouge_evalue(score, kbp, gbp, query_length, subject_length);
    raw / gap_decay_divisor(BLAST_GAP_DECAY_RATE_GAPPED, 1)
}

/// 1-1 port of `BLAST_SpougeStoE` (`blast_stat.c:5176`). Returns the raw
/// per-pair Spouge finite-size-corrected e-value. Some NCBI call sites
/// (e.g. translated searches with HSP linking via `link_hsps.c:1791`)
/// further divide by `BLAST_GapDecayDivisor(gap_decay_rate, 1)` —
/// callers in those paths should use [`spouge_evalue_with_gap_decay`]
/// instead. Standard blastp / direct traceback paths
/// (`blast_traceback.c:234`, `blast_kappa.c:419`) pass
/// `gap_decay_rate=0` and use the raw value unchanged.
pub fn spouge_evalue(
    score: i32,
    kbp: &KarlinBlk,
    gbp: &GumbelBlk,
    query_length: i32,
    subject_length: i32,
) -> f64 {
    let scale_factor = kbp.lambda / gbp.lambda;
    let db_scale_factor = if gbp.db_length > 0 {
        gbp.db_length as f64 / subject_length as f64
    } else {
        1.0
    };

    let ai_hat = gbp.a * scale_factor;
    let bi_hat = gbp.b;
    let alphai_hat = gbp.alpha * scale_factor;
    let betai_hat = gbp.beta;
    let sigma_hat = gbp.sigma * scale_factor;
    let tau_hat = gbp.tau;

    let aj_hat = ai_hat;
    let bj_hat = bi_hat;
    let alphaj_hat = alphai_hat;
    let betaj_hat = betai_hat;

    const CONST_VAL: f64 = 0.39894228040143267793994605993438; // 1/sqrt(2*PI)

    let y = score as f64;
    let m = query_length as f64;
    let n = subject_length as f64;

    let m_li_y = m - (ai_hat * y + bi_hat);
    let vi_y = (2.0 * alphai_hat / kbp.lambda).max(alphai_hat * y + betai_hat);
    let sqrt_vi_y = vi_y.sqrt();
    let m_f = m_li_y / sqrt_vi_y;
    let p_m_f = crate::math::erfc(-m_f / std::f64::consts::SQRT_2) / 2.0;
    let p1 = m_li_y * p_m_f + sqrt_vi_y * CONST_VAL * (-0.5 * m_f * m_f).exp();

    let n_lj_y = n - (aj_hat * y + bj_hat);
    let vj_y = (2.0 * alphaj_hat / kbp.lambda).max(alphaj_hat * y + betaj_hat);
    let sqrt_vj_y = vj_y.sqrt();
    let n_f = n_lj_y / sqrt_vj_y;
    let p_n_f = crate::math::erfc(-n_f / std::f64::consts::SQRT_2) / 2.0;
    let p2 = n_lj_y * p_n_f + sqrt_vj_y * CONST_VAL * (-0.5 * n_f * n_f).exp();

    let c_y = (2.0 * sigma_hat / kbp.lambda).max(sigma_hat * y + tau_hat);
    let area = p1 * p2 + c_y * p_m_f * p_n_f;

    (area * kbp.k * (-kbp.lambda * y).exp() * db_scale_factor).max(0.0)
}

/// Port of NCBI `BLAST_SpougeEtoS` (`blast_stat.c:5236`). Binary-search the
/// raw score `S` such that `BLAST_SpougeStoE(S, kbp, gbp, m, n) <= e0`.
///
/// Used for the per-context cutoff in `BlastHitSavingParametersUpdate`
/// (`blast_parameters.c:940`): the cutoff is computed at `m=query_length` and
/// `n=avg_subject_length`, with `gbp->db_length` set to the full DB. The
/// binary search uses Spouge's full FSC formula, so the cutoff is more
/// conservative than the simple Karlin `evalue_to_raw` formula.
pub fn spouge_etos(e0: f64, kbp: &KarlinBlk, gbp: &GumbelBlk, m: i32, n: i32) -> i32 {
    let db_scale_factor = if gbp.db_length > 0 {
        gbp.db_length as f64
    } else {
        1.0
    };

    // C: `b = MAX((int)(log(db_scale_factor/e0) / kbp->Lambda), 2);`
    let mut b = ((db_scale_factor / e0).ln() / kbp.lambda) as i32;
    if b < 2 {
        b = 2;
    }

    let mut a: i32 = 0;
    let mut e = spouge_evalue(b, kbp, gbp, m, n);

    if e > e0 {
        while e > e0 {
            a = b;
            // C: `b *= 2;` — note the integer overflow possibility on a
            // pathological input is mirrored from C verbatim
            b = b.saturating_mul(2);
            e = spouge_evalue(b, kbp, gbp, m, n);
            if b == i32::MAX {
                break;
            }
        }
    }
    while b - a > 1 {
        let c = (a + b) / 2;
        e = spouge_evalue(c, kbp, gbp, m, n);
        if e > e0 {
            a = c;
        } else {
            b = c;
        }
    }
    let _ = e; // suppress unused on the final iteration
    a
}

/// Build Gumbel block for protein BLOSUM62 with given gap costs.
/// Port of Blast_GumbelBlkLoadFromTables from blast_stat.c:3696.
pub fn protein_gumbel_blk(gap_open: i32, gap_extend: i32, db_length: i64) -> Option<GumbelBlk> {
    matrix_gumbel_blk("BLOSUM62", gap_open, gap_extend, db_length)
}

/// Build Gumbel block for the named protein scoring matrix.
/// Port of `Blast_GumbelBlkLoadFromTables` using NCBI `blast_stat.c`
/// matrix value rows.
pub fn matrix_gumbel_blk(
    matrix_name: &str,
    gap_open: i32,
    gap_extend: i32,
    db_length: i64,
) -> Option<GumbelBlk> {
    let rows = matrix_stat_rows(matrix_name)?;
    let ungapped = rows.first()?;
    let row = lookup_matrix_stat_row(matrix_name, gap_open, gap_extend)?;
    let g = (gap_open + gap_extend) as f64;

    Some(GumbelBlk {
        lambda: row.lambda,
        a: row.alpha,
        b: 2.0 * g * (ungapped.alpha - row.alpha),
        alpha: row.alpha_v,
        beta: 2.0 * g * (ungapped.alpha_v - row.alpha_v),
        sigma: row.sigma,
        tau: 2.0 * g * (ungapped.alpha_v - row.sigma),
        db_length,
    })
}

/// Build Gumbel block for NCBI's IDENTITY matrix.
/// Port of `Blast_GumbelBlkLoadFromTables` using `prot_idenity_values[]`.
pub fn identity_gumbel_blk(gap_open: i32, gap_extend: i32, db_length: i64) -> Option<GumbelBlk> {
    matrix_gumbel_blk("IDENTITY", gap_open, gap_extend, db_length)
}

/// Compute ungapped Karlin-Altschul lambda parameter for nucleotide scoring.
/// Solves: 0.25 * exp(lambda * reward) + 0.75 * exp(lambda * penalty) = 1
/// using Newton's method.
pub fn compute_ungapped_lambda(reward: i32, penalty: i32) -> f64 {
    let r = reward as f64;
    let p = penalty as f64;
    // Initial guess
    let mut lambda = 1.0;
    for _ in 0..100 {
        let er = (lambda * r).exp();
        let ep = (lambda * p).exp();
        let f = 0.25 * er + 0.75 * ep - 1.0;
        let fp = 0.25 * r * er + 0.75 * p * ep;
        if fp.abs() < 1e-30 {
            break;
        }
        let delta = f / fp;
        lambda -= delta;
        if delta.abs() < 1e-12 {
            break;
        }
    }
    lambda
}

/// Compute ungapped Karlin-Altschul K parameter for nucleotide scoring.
/// K = (1 - exp(-lambda)) * (1 - exp(-lambda)) / (sum of prob * exp(lambda*score) * score^2 / 2)
/// Simplified for uniform base frequencies.
pub fn compute_ungapped_k(lambda: f64, reward: i32, penalty: i32) -> f64 {
    let r = reward as f64;
    let p = penalty as f64;
    let er = (lambda * r).exp();
    let ep = (lambda * p).exp();
    // H = lambda * (0.25 * r * er + 0.75 * p * ep)
    let h = lambda * (0.25 * r * er + 0.75 * p * ep);
    // K approximation for simple scoring
    // K = H / (lambda * lambda * variance)
    // For nucleotide with uniform freqs, use the standard approximation
    let variance = 0.25 * r * r * er + 0.75 * p * p * ep;
    if variance.abs() < 1e-30 {
        return 0.1;
    }
    h / (lambda * variance)
}

/// Compute ungapped KarlinBlk for nucleotide scoring with given reward/penalty.
pub fn compute_ungapped_kbp(reward: i32, penalty: i32) -> KarlinBlk {
    let lambda = compute_ungapped_lambda(reward, penalty);
    let k = compute_ungapped_k(lambda, reward, penalty);
    let h = lambda
        * (0.25 * (reward as f64) * (lambda * reward as f64).exp()
            + 0.75 * (penalty as f64) * (lambda * penalty as f64).exp());
    KarlinBlk {
        lambda,
        k,
        log_k: k.ln(),
        h,
        round_down: false,
    }
}

/// Compute length adjustment using the Altschul-Gish formula.
/// This adjusts query and database lengths to account for edge effects.
/// length_adj satisfies: K * (q - length_adj) * (d - N * length_adj) * exp(-lambda * S) = threshold
/// Compute the length adjustment using NCBI `BLAST_ComputeLengthAdjustment`
/// (`blast_stat.c:5041`). For callers that have no explicit alpha/beta
/// (typically ungapped or missing-gapped-params paths), NCBI uses
/// `alpha/lambda = 1/H` and `beta = 0`; that convention is applied
/// here so the verbatim bracketing algorithm is exercised.
pub fn compute_length_adjustment(
    query_length: i32,
    db_length: i64,
    num_seqs: i32,
    kbp: &KarlinBlk,
) -> i32 {
    if !kbp.is_valid() || query_length <= 0 || db_length <= 0 {
        return 0;
    }
    if kbp.h <= 0.0 {
        return 0;
    }
    let (adj, _converged) = compute_length_adjustment_exact(
        kbp.k,
        kbp.log_k,
        1.0 / kbp.h,
        0.0,
        query_length,
        db_length,
        num_seqs,
    );
    adj
}

const I2MAX: i32 = i16::MAX as i32;

/// NCBI `blast_stat.c` protein matrix value rows. The first row is the
/// ungapped row; remaining rows are supported gapped costs.
const BLOSUM45_ROWS: &[MatrixStatRow] = &[
    MatrixStatRow::new(
        I2MAX, I2MAX, 0.2291, 0.0924, 0.2514, 0.9113, -5.7, 0.641318, 9.611060, 9.611060,
    ),
    MatrixStatRow::new(
        13, 3, 0.207, 0.049, 0.14, 1.5, -22.0, 0.671128, 35.855900, 35.963900,
    ),
    MatrixStatRow::new(
        12, 3, 0.199, 0.039, 0.11, 1.8, -34.0, 0.691530, 45.693600, 45.851700,
    ),
    MatrixStatRow::new(
        11, 3, 0.190, 0.031, 0.095, 2.0, -38.0, 0.691181, 62.874100, 63.103700,
    ),
    MatrixStatRow::new(
        10, 3, 0.179, 0.023, 0.075, 2.4, -51.0, 0.710529, 88.286800, 88.639100,
    ),
    MatrixStatRow::new(
        16, 2, 0.210, 0.051, 0.14, 1.5, -24.0, 0.666680, 36.279800, 36.452400,
    ),
    MatrixStatRow::new(
        15, 2, 0.203, 0.041, 0.12, 1.7, -31.0, 0.673871, 44.825700, 45.060400,
    ),
    MatrixStatRow::new(
        14, 2, 0.195, 0.032, 0.10, 1.9, -36.0, 0.685753, 60.736200, 61.102300,
    ),
    MatrixStatRow::new(
        13, 2, 0.185, 0.024, 0.084, 2.2, -45.0, 0.698480, 85.148100, 85.689400,
    ),
    MatrixStatRow::new(
        12, 2, 0.171, 0.016, 0.061, 2.8, -65.0, 0.713429, 127.758000, 128.582000,
    ),
    MatrixStatRow::new(
        19, 1, 0.205, 0.040, 0.11, 1.9, -43.0, 0.672302, 53.071400, 53.828200,
    ),
    MatrixStatRow::new(
        18, 1, 0.198, 0.032, 0.10, 2.0, -43.0, 0.682580, 72.342400, 73.403900,
    ),
    MatrixStatRow::new(
        17, 1, 0.189, 0.024, 0.079, 2.4, -57.0, 0.695035, 103.055000, 104.721000,
    ),
    MatrixStatRow::new(
        16, 1, 0.176, 0.016, 0.063, 2.8, -67.0, 0.712966, 170.100000, 173.003000,
    ),
];

const BLOSUM50_ROWS: &[MatrixStatRow] = &[
    MatrixStatRow::new(
        I2MAX, I2MAX, 0.2318, 0.112, 0.3362, 0.6895, -4.0, 0.609639, 5.388310, 5.388310,
    ),
    MatrixStatRow::new(
        13, 3, 0.212, 0.063, 0.19, 1.1, -16.0, 0.639287, 18.113800, 18.202800,
    ),
    MatrixStatRow::new(
        12, 3, 0.206, 0.055, 0.17, 1.2, -18.0, 0.644715, 22.654600, 22.777700,
    ),
    MatrixStatRow::new(
        11, 3, 0.197, 0.042, 0.14, 1.4, -25.0, 0.656327, 29.861100, 30.045700,
    ),
    MatrixStatRow::new(
        10, 3, 0.186, 0.031, 0.11, 1.7, -34.0, 0.671150, 42.393800, 42.674000,
    ),
    MatrixStatRow::new(
        9, 3, 0.172, 0.022, 0.082, 2.1, -48.0, 0.694326, 66.069600, 66.516400,
    ),
    MatrixStatRow::new(
        16, 2, 0.215, 0.066, 0.20, 1.05, -15.0, 0.633899, 17.951800, 18.092100,
    ),
    MatrixStatRow::new(
        15, 2, 0.210, 0.058, 0.17, 1.2, -20.0, 0.641985, 21.940100, 22.141800,
    ),
    MatrixStatRow::new(
        14, 2, 0.202, 0.045, 0.14, 1.4, -27.0, 0.650682, 28.681200, 28.961900,
    ),
    MatrixStatRow::new(
        13, 2, 0.193, 0.035, 0.12, 1.6, -32.0, 0.660984, 42.059500, 42.471600,
    ),
    MatrixStatRow::new(
        12, 2, 0.181, 0.025, 0.095, 1.9, -41.0, 0.678090, 63.747600, 64.397300,
    ),
    MatrixStatRow::new(
        19, 1, 0.212, 0.057, 0.18, 1.2, -21.0, 0.635714, 26.311200, 26.923300,
    ),
    MatrixStatRow::new(
        18, 1, 0.207, 0.050, 0.15, 1.4, -28.0, 0.643523, 34.903700, 35.734800,
    ),
    MatrixStatRow::new(
        17, 1, 0.198, 0.037, 0.12, 1.6, -33.0, 0.654504, 48.895800, 50.148600,
    ),
    MatrixStatRow::new(
        16, 1, 0.186, 0.025, 0.10, 1.9, -42.0, 0.667750, 76.469100, 78.443000,
    ),
    MatrixStatRow::new(
        15, 1, 0.171, 0.015, 0.063, 2.7, -76.0, 0.694575, 140.053000, 144.160000,
    ),
];

const BLOSUM62_ROWS: &[MatrixStatRow] = &[
    MatrixStatRow::new(
        I2MAX, I2MAX, 0.3176, 0.134, 0.4012, 0.7916, -3.2, 0.623757, 4.964660, 4.964660,
    ),
    MatrixStatRow::new(
        11, 2, 0.297, 0.082, 0.27, 1.1, -10.0, 0.641766, 12.673800, 12.757600,
    ),
    MatrixStatRow::new(
        10, 2, 0.291, 0.075, 0.23, 1.3, -15.0, 0.649362, 16.474000, 16.602600,
    ),
    MatrixStatRow::new(
        9, 2, 0.279, 0.058, 0.19, 1.5, -19.0, 0.659245, 22.751900, 22.950000,
    ),
    MatrixStatRow::new(
        8, 2, 0.264, 0.045, 0.15, 1.8, -26.0, 0.672692, 35.483800, 35.821300,
    ),
    MatrixStatRow::new(
        7, 2, 0.239, 0.027, 0.10, 2.5, -46.0, 0.702056, 61.238300, 61.886000,
    ),
    MatrixStatRow::new(
        6, 2, 0.201, 0.012, 0.061, 3.3, -58.0, 0.740802, 140.417000, 141.882000,
    ),
    MatrixStatRow::new(
        13, 1, 0.292, 0.071, 0.23, 1.2, -11.0, 0.647715, 19.506300, 19.893100,
    ),
    MatrixStatRow::new(
        12, 1, 0.283, 0.059, 0.19, 1.5, -19.0, 0.656391, 27.856200, 28.469900,
    ),
    MatrixStatRow::new(
        11, 1, 0.267, 0.041, 0.14, 1.9, -30.0, 0.669720, 42.602800, 43.636200,
    ),
    MatrixStatRow::new(
        10, 1, 0.243, 0.024, 0.10, 2.5, -44.0, 0.693267, 83.178700, 85.065600,
    ),
    MatrixStatRow::new(
        9, 1, 0.206, 0.010, 0.052, 4.0, -87.0, 0.731887, 210.333000, 214.842000,
    ),
];

const BLOSUM80_ROWS: &[MatrixStatRow] = &[
    MatrixStatRow::new(
        I2MAX, I2MAX, 0.3430, 0.177, 0.6568, 0.5222, -1.6, 0.564057, 1.918130, 1.918130,
    ),
    MatrixStatRow::new(
        25, 2, 0.342, 0.17, 0.66, 0.52, -1.6, 0.563956, 1.731000, 1.731300,
    ),
    MatrixStatRow::new(
        13, 2, 0.336, 0.15, 0.57, 0.59, -3.0, 0.570979, 2.673470, 2.692300,
    ),
    MatrixStatRow::new(
        9, 2, 0.319, 0.11, 0.42, 0.76, -6.0, 0.587837, 5.576090, 5.667860,
    ),
    MatrixStatRow::new(
        8, 2, 0.308, 0.090, 0.35, 0.89, -9.0, 0.597556, 7.536950, 7.686230,
    ),
    MatrixStatRow::new(
        7, 2, 0.293, 0.070, 0.27, 1.1, -14.0, 0.615254, 11.586600, 11.840400,
    ),
    MatrixStatRow::new(
        6, 2, 0.268, 0.045, 0.19, 1.4, -19.0, 0.644054, 19.958100, 20.441200,
    ),
    MatrixStatRow::new(
        11, 1, 0.314, 0.095, 0.35, 0.90, -9.0, 0.590702, 8.808610, 9.223320,
    ),
    MatrixStatRow::new(
        10, 1, 0.299, 0.071, 0.27, 1.1, -14.0, 0.609620, 13.833800, 14.533400,
    ),
    MatrixStatRow::new(
        9, 1, 0.279, 0.048, 0.20, 1.4, -19.0, 0.623800, 24.252000, 25.490400,
    ),
];

const BLOSUM90_ROWS: &[MatrixStatRow] = &[
    MatrixStatRow::new(
        I2MAX, I2MAX, 0.3346, 0.190, 0.7547, 0.4434, -1.4, 0.544178, 1.377760, 1.377760,
    ),
    MatrixStatRow::new(
        9, 2, 0.310, 0.12, 0.46, 0.67, -6.0, 0.570267, 4.232290, 4.334170,
    ),
    MatrixStatRow::new(
        8, 2, 0.300, 0.099, 0.39, 0.76, -7.0, 0.581580, 5.797020, 5.961420,
    ),
    MatrixStatRow::new(
        7, 2, 0.283, 0.072, 0.30, 0.93, -11.0, 0.600024, 9.040880, 9.321600,
    ),
    MatrixStatRow::new(
        6, 2, 0.259, 0.048, 0.22, 1.2, -16.0, 0.629344, 16.024400, 16.531600,
    ),
    MatrixStatRow::new(
        11, 1, 0.302, 0.093, 0.39, 0.78, -8.0, 0.576919, 7.143250, 7.619190,
    ),
    MatrixStatRow::new(
        10, 1, 0.290, 0.075, 0.28, 1.04, -15.0, 0.591366, 11.483900, 12.269800,
    ),
    MatrixStatRow::new(
        9, 1, 0.265, 0.044, 0.20, 1.3, -19.0, 0.613013, 21.408300, 22.840900,
    ),
];

const PAM250_ROWS: &[MatrixStatRow] = &[
    MatrixStatRow::new(
        I2MAX, I2MAX, 0.2252, 0.0868, 0.2223, 0.98, -5.0, 0.660059, 11.754300, 11.754300,
    ),
    MatrixStatRow::new(
        15, 3, 0.205, 0.049, 0.13, 1.6, -23.0, 0.687656, 34.578400, 34.928000,
    ),
    MatrixStatRow::new(
        14, 3, 0.200, 0.043, 0.12, 1.7, -26.0, 0.689768, 43.353000, 43.443800,
    ),
    MatrixStatRow::new(
        13, 3, 0.194, 0.036, 0.10, 1.9, -31.0, 0.697431, 50.948500, 51.081700,
    ),
    MatrixStatRow::new(
        12, 3, 0.186, 0.029, 0.085, 2.2, -41.0, 0.704565, 69.606500, 69.793600,
    ),
    MatrixStatRow::new(
        11, 3, 0.174, 0.020, 0.070, 2.5, -48.0, 0.722438, 98.653500, 98.927100,
    ),
    MatrixStatRow::new(
        17, 2, 0.204, 0.047, 0.12, 1.7, -28.0, 0.684799, 41.583800, 41.735800,
    ),
    MatrixStatRow::new(
        16, 2, 0.198, 0.038, 0.11, 1.8, -29.0, 0.691098, 51.635200, 51.843900,
    ),
    MatrixStatRow::new(
        15, 2, 0.191, 0.031, 0.087, 2.2, -44.0, 0.699051, 67.256700, 67.558500,
    ),
    MatrixStatRow::new(
        14, 2, 0.182, 0.024, 0.073, 2.5, -53.0, 0.714103, 96.315100, 96.756800,
    ),
    MatrixStatRow::new(
        13, 2, 0.171, 0.017, 0.059, 2.9, -64.0, 0.728738, 135.653000, 136.339000,
    ),
    MatrixStatRow::new(
        21, 1, 0.205, 0.045, 0.11, 1.8, -34.0, 0.683265, 48.728200, 49.218800,
    ),
    MatrixStatRow::new(
        20, 1, 0.199, 0.037, 0.10, 1.9, -35.0, 0.689380, 60.832000, 61.514100,
    ),
    MatrixStatRow::new(
        19, 1, 0.192, 0.029, 0.083, 2.3, -52.0, 0.696344, 84.019700, 84.985600,
    ),
    MatrixStatRow::new(
        18, 1, 0.183, 0.021, 0.070, 2.6, -60.0, 0.710525, 113.829000, 115.184000,
    ),
    MatrixStatRow::new(
        17, 1, 0.171, 0.014, 0.052, 3.3, -86.0, 0.727000, 175.071000, 177.196000,
    ),
];

const PAM30_ROWS: &[MatrixStatRow] = &[
    MatrixStatRow::new(
        I2MAX, I2MAX, 0.3400, 0.283, 1.754, 0.1938, -0.3, 0.436164, 0.161818, 0.161818,
    ),
    MatrixStatRow::new(
        7, 2, 0.305, 0.15, 0.87, 0.35, -3.0, 0.479087, 1.014010, 1.162730,
    ),
    MatrixStatRow::new(
        6, 2, 0.287, 0.11, 0.68, 0.42, -4.0, 0.499980, 1.688060, 1.951430,
    ),
    MatrixStatRow::new(
        5, 2, 0.264, 0.079, 0.45, 0.59, -7.0, 0.533009, 3.377010, 3.871950,
    ),
    MatrixStatRow::new(
        10, 1, 0.309, 0.15, 0.88, 0.35, -3.0, 0.474741, 1.372050, 1.788770,
    ),
    MatrixStatRow::new(
        9, 1, 0.294, 0.11, 0.61, 0.48, -6.0, 0.492716, 2.463920, 3.186150,
    ),
    MatrixStatRow::new(
        8, 1, 0.270, 0.072, 0.40, 0.68, -10.0, 0.521286, 5.368130, 6.763480,
    ),
    MatrixStatRow::new(
        15, 3, 0.339, 0.28, 1.70, 0.20, -0.5, 0.437688, 0.157089, 0.155299,
    ),
    MatrixStatRow::new(
        14, 2, 0.337, 0.27, 1.62, 0.21, -0.8, 0.440010, 0.206970, 0.198524,
    ),
    MatrixStatRow::new(
        14, 1, 0.333, 0.27, 1.43, 0.23, -1.4, 0.444817, 0.436301, 0.361947,
    ),
    MatrixStatRow::new(
        13, 3, 0.338, 0.27, 1.69, 0.20, -0.5, 0.439086, 0.178973, 0.175436,
    ),
];

const PAM70_ROWS: &[MatrixStatRow] = &[
    MatrixStatRow::new(
        I2MAX, I2MAX, 0.3345, 0.229, 1.029, 0.3250, -0.7, 0.511296, 0.633439, 0.633439,
    ),
    MatrixStatRow::new(
        8, 2, 0.301, 0.12, 0.54, 0.56, -5.0, 0.549019, 2.881650, 3.025710,
    ),
    MatrixStatRow::new(
        7, 2, 0.286, 0.093, 0.43, 0.67, -7.0, 0.565659, 4.534540, 4.785780,
    ),
    MatrixStatRow::new(
        6, 2, 0.264, 0.064, 0.29, 0.90, -12.0, 0.596330, 7.942630, 8.402720,
    ),
    MatrixStatRow::new(
        11, 1, 0.305, 0.12, 0.52, 0.59, -6.0, 0.543514, 3.681400, 4.108020,
    ),
    MatrixStatRow::new(
        10, 1, 0.291, 0.091, 0.41, 0.71, -9.0, 0.560723, 6.002970, 6.716570,
    ),
    MatrixStatRow::new(
        9, 1, 0.270, 0.060, 0.28, 0.97, -14.0, 0.585186, 11.360800, 12.636700,
    ),
    MatrixStatRow::new(
        11, 2, 0.323, 0.186, 0.80, 1.32, -27.0, 0.524062, 1.321301, 1.281671,
    ),
    MatrixStatRow::new(
        12, 3, 0.330, 0.219, 0.93, 0.82, -16.0, 0.516845, 0.818768, 0.811240,
    ),
];

const IDENTITY_ROWS: &[MatrixStatRow] = &[
    MatrixStatRow::new(
        I2MAX, I2MAX, 0.28768, 0.282, 1.69, 0.1703, -0.3, 0.43828, 0.16804, 0.16804,
    ),
    MatrixStatRow::new(
        15, 2, 0.2835, 0.255, 1.49, 0.19, -1.0, 0.44502, 0.24613, 0.22743,
    ),
];

fn matrix_stat_rows(matrix_name: &str) -> Option<&'static [MatrixStatRow]> {
    if matrix_name.eq_ignore_ascii_case("BLOSUM45") {
        Some(BLOSUM45_ROWS)
    } else if matrix_name.eq_ignore_ascii_case("BLOSUM50") {
        Some(BLOSUM50_ROWS)
    } else if matrix_name.eq_ignore_ascii_case("BLOSUM62") || matrix_name.is_empty() {
        Some(BLOSUM62_ROWS)
    } else if matrix_name.eq_ignore_ascii_case("BLOSUM80") {
        Some(BLOSUM80_ROWS)
    } else if matrix_name.eq_ignore_ascii_case("BLOSUM90") {
        Some(BLOSUM90_ROWS)
    } else if matrix_name.eq_ignore_ascii_case("PAM30") {
        Some(PAM30_ROWS)
    } else if matrix_name.eq_ignore_ascii_case("PAM70") {
        Some(PAM70_ROWS)
    } else if matrix_name.eq_ignore_ascii_case("PAM250") {
        Some(PAM250_ROWS)
    } else if matrix_name.eq_ignore_ascii_case("IDENTITY") {
        Some(IDENTITY_ROWS)
    } else {
        None
    }
}

fn lookup_matrix_stat_row(
    matrix_name: &str,
    gap_open: i32,
    gap_extend: i32,
) -> Option<MatrixStatRow> {
    matrix_stat_rows(matrix_name)?
        .iter()
        .copied()
        .find(|row| row.gap_open == gap_open && row.gap_extend == gap_extend)
}

/// Precomputed Karlin-Altschul parameters for BLOSUM62 with various gap costs.
/// Values from NCBI blast_stat.c blosum62_values[].
/// Format: (gap_open, gap_extend, lambda, K, H, alpha, beta)
pub const BLOSUM62_PARAMS: &[(i32, i32, f64, f64, f64, f64, f64)] = &[
    (11, 2, 0.297, 0.082, 0.27, 1.1, -10.0),
    (10, 2, 0.291, 0.075, 0.23, 1.3, -15.0),
    (9, 2, 0.279, 0.058, 0.19, 1.5, -19.0),
    (8, 2, 0.264, 0.045, 0.15, 1.8, -26.0),
    (7, 2, 0.239, 0.027, 0.10, 2.5, -46.0),
    (6, 2, 0.201, 0.012, 0.061, 3.3, -58.0),
    (13, 1, 0.292, 0.071, 0.23, 1.2, -11.0),
    (12, 1, 0.283, 0.059, 0.19, 1.5, -19.0),
    (11, 1, 0.267, 0.041, 0.14, 1.9, -30.0),
    (10, 1, 0.243, 0.024, 0.10, 2.5, -44.0),
    (9, 1, 0.206, 0.010, 0.052, 4.0, -87.0),
];

/// Precomputed Karlin-Altschul parameters for NCBI's IDENTITY matrix.
/// Values from `blast_stat.c` `prot_idenity_values[]`.
pub const IDENTITY_PARAMS: &[(i32, i32, f64, f64, f64, f64, f64)] = &[
    (
        i16::MAX as i32,
        i16::MAX as i32,
        0.28768,
        0.282,
        1.69,
        0.1703,
        -0.3,
    ),
    (15, 2, 0.2835, 0.255, 1.49, 0.19, -1.0),
];

/// Look up gapped KBP for protein scoring (BLOSUM62).
pub fn lookup_protein_params(gap_open: i32, gap_extend: i32) -> Option<GappedParams> {
    lookup_matrix_params("BLOSUM62", gap_open, gap_extend)
}

/// Look up gapped KBP for a named protein scoring matrix.
pub fn lookup_matrix_params(
    matrix_name: &str,
    gap_open: i32,
    gap_extend: i32,
) -> Option<GappedParams> {
    lookup_matrix_stat_row(matrix_name, gap_open, gap_extend).map(MatrixStatRow::gapped_params)
}

/// Look up gapped KBP for NCBI's IDENTITY matrix.
pub fn lookup_identity_params(gap_open: i32, gap_extend: i32) -> Option<GappedParams> {
    lookup_matrix_params("IDENTITY", gap_open, gap_extend)
}

/// Compute ungapped KBP for protein BLOSUM62. Values match NCBI
/// `blosum62_values[0]` (ungapped entry) from `blast_stat.c:259`:
/// Lambda = 0.3176, K = 0.134, H = 0.4012.
pub fn protein_ungapped_kbp() -> KarlinBlk {
    protein_ungapped_kbp_for_matrix("BLOSUM62")
}

/// Compute ungapped KBP for the named protein matrix from NCBI's standard
/// `blast_stat.c` row 0 table values.
pub fn protein_ungapped_kbp_for_matrix(matrix_name: &str) -> KarlinBlk {
    let row = matrix_stat_rows(matrix_name)
        .and_then(|rows| rows.first().copied())
        .unwrap_or(BLOSUM62_ROWS[0]);
    KarlinBlk {
        lambda: row.lambda,
        k: row.k,
        log_k: row.k.ln(),
        h: row.h,
        round_down: false,
    }
}

/// Compute the query-specific ungapped Karlin block for a protein query.
/// Mirrors NCBI's `Blast_ScoreBlkKbpUngappedCalc` (`blast_stat.c:2737`):
/// `Blast_ResFreqString` over the query → `BlastScoreFreqCalc` →
/// `Blast_KarlinBlkUngappedCalc`. The result drifts slightly from the ideal
/// kbp depending on query composition; the drift is what makes `gap_trigger`
/// match NCBI bit-for-bit at boundary cases (see iter 99).
pub fn query_specific_protein_ungapped_kbp(
    query_aa_ncbistdaa: &[u8],
    matrix: &[[i32; crate::matrix::AA_SIZE]; crate::matrix::AA_SIZE],
) -> KarlinBlk {
    let std_freq = protein_std_freq_ncbistdaa();
    // NCBIstdaa ambiguity letters that NCBI's Blast_ResFreqString skips.
    let ambiguous: [u8; 8] = [0, 2, 21, 23, 24, 25, 26, 27];
    let contexts = [UngappedKbpContext {
        query_offset: 0,
        query_length: query_aa_ncbistdaa.len() as i32,
        is_valid: true,
    }];
    let mat = |i: usize, j: usize| matrix[i][j];
    let r = ungapped_kbp_calc_with_std(
        query_aa_ncbistdaa,
        &contexts,
        -4,
        11,
        crate::matrix::AA_SIZE,
        &ambiguous,
        &std_freq,
        &mat,
    );
    r.into_iter()
        .next()
        .flatten()
        .unwrap_or_else(protein_ungapped_kbp)
}

/// NCBIstdaa-indexed standard amino-acid background frequencies.
/// Robinson & Robinson 1991 (matches NCBI's `Blast_ResFreqStdComp` table for
/// `BLASTAA_SEQ_CODE`). Entries for `*`, `-`, ambiguity letters (B/Z/X/U/O/J)
/// are zero — the C path leaves them out of the score-frequency calculation.
pub fn protein_std_freq_ncbistdaa() -> [f64; 28] {
    use crate::matrix::AA_FREQUENCIES;
    // AA_FREQUENCIES is alphabetical (A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y).
    // NCBIstdaa indices: -=0 A=1 B=2 C=3 D=4 E=5 F=6 G=7 H=8 I=9 K=10 L=11 M=12
    // N=13 P=14 Q=15 R=16 S=17 T=18 V=19 W=20 X=21 Y=22 Z=23 U=24 *=25 O=26 J=27
    let mut out = [0.0f64; 28];
    out[1] = AA_FREQUENCIES[0]; // A
    out[3] = AA_FREQUENCIES[1]; // C
    out[4] = AA_FREQUENCIES[2]; // D
    out[5] = AA_FREQUENCIES[3]; // E
    out[6] = AA_FREQUENCIES[4]; // F
    out[7] = AA_FREQUENCIES[5]; // G
    out[8] = AA_FREQUENCIES[6]; // H
    out[9] = AA_FREQUENCIES[7]; // I
    out[10] = AA_FREQUENCIES[8]; // K
    out[11] = AA_FREQUENCIES[9]; // L
    out[12] = AA_FREQUENCIES[10]; // M
    out[13] = AA_FREQUENCIES[11]; // N
    out[14] = AA_FREQUENCIES[12]; // P
    out[15] = AA_FREQUENCIES[13]; // Q
    out[16] = AA_FREQUENCIES[14]; // R
    out[17] = AA_FREQUENCIES[15]; // S
    out[18] = AA_FREQUENCIES[16]; // T
    out[19] = AA_FREQUENCIES[17]; // V
    out[20] = AA_FREQUENCIES[18]; // W
    out[22] = AA_FREQUENCIES[19]; // Y
    out
}

/// Compute effective search space.
/// eff_length = db_length - num_seqs * length_adjustment
/// eff_query_length = query_length - length_adjustment
/// search_space = eff_length * eff_query_length
pub fn compute_search_space(
    query_length: i64,
    db_length: i64,
    num_seqs: i32,
    length_adjustment: i32,
) -> f64 {
    let eff_query = (query_length - length_adjustment as i64).max(1);
    let eff_db = (db_length - num_seqs as i64 * length_adjustment as i64).max(1);
    eff_query as f64 * eff_db as f64
}

// ---------------------------------------------------------------------------
// Gapped Karlin-Altschul parameter lookup (exact C-compatible tables)
// Port of Blast_KarlinBlkNuclGappedCalc from blast_stat.c
// ---------------------------------------------------------------------------

/// Each row: [gap_open, gap_extend, lambda, K, H, alpha, beta, theta]
type KbpTableRow = [f64; 8];

const KBPT_1_5: &[KbpTableRow] = &[
    [0.0, 0.0, 1.39, 0.747, 1.38, 1.00, 0.0, 100.0],
    [3.0, 3.0, 1.39, 0.747, 1.38, 1.00, 0.0, 100.0],
];
const KBPT_1_4: &[KbpTableRow] = &[
    [0.0, 0.0, 1.383, 0.738, 1.36, 1.02, 0.0, 100.0],
    [1.0, 2.0, 1.36, 0.67, 1.2, 1.1, 0.0, 98.0],
    [0.0, 2.0, 1.26, 0.43, 0.90, 1.4, -1.0, 91.0],
    [2.0, 1.0, 1.35, 0.61, 1.1, 1.2, -1.0, 98.0],
    [1.0, 1.0, 1.22, 0.35, 0.72, 1.7, -3.0, 88.0],
];
const KBPT_2_7: &[KbpTableRow] = &[
    [0.0, 0.0, 0.69, 0.73, 1.34, 0.515, 0.0, 100.0],
    [2.0, 4.0, 0.68, 0.67, 1.2, 0.55, 0.0, 99.0],
    [0.0, 4.0, 0.63, 0.43, 0.90, 0.7, -1.0, 91.0],
    [4.0, 2.0, 0.675, 0.62, 1.1, 0.6, -1.0, 98.0],
    [2.0, 2.0, 0.61, 0.35, 0.72, 1.7, -3.0, 88.0],
];
const KBPT_1_3: &[KbpTableRow] = &[
    [0.0, 0.0, 1.374, 0.711, 1.31, 1.05, 0.0, 100.0],
    [2.0, 2.0, 1.37, 0.70, 1.2, 1.1, 0.0, 99.0],
    [1.0, 2.0, 1.35, 0.64, 1.1, 1.2, -1.0, 98.0],
    [0.0, 2.0, 1.25, 0.42, 0.83, 1.5, -2.0, 91.0],
    [2.0, 1.0, 1.34, 0.60, 1.1, 1.2, -1.0, 97.0],
    [1.0, 1.0, 1.21, 0.34, 0.71, 1.7, -2.0, 88.0],
];
const KBPT_2_5: &[KbpTableRow] = &[
    [0.0, 0.0, 0.675, 0.65, 1.1, 0.6, -1.0, 99.0],
    [2.0, 4.0, 0.67, 0.59, 1.1, 0.6, -1.0, 98.0],
    [0.0, 4.0, 0.62, 0.39, 0.78, 0.8, -2.0, 91.0],
    [4.0, 2.0, 0.67, 0.61, 1.0, 0.65, -2.0, 98.0],
    [2.0, 2.0, 0.56, 0.32, 0.59, 0.95, -4.0, 82.0],
];
const KBPT_1_2: &[KbpTableRow] = &[
    [0.0, 0.0, 1.28, 0.46, 0.85, 1.5, -2.0, 96.0],
    [2.0, 2.0, 1.33, 0.62, 1.1, 1.2, 0.0, 99.0],
    [1.0, 2.0, 1.30, 0.52, 0.93, 1.4, -2.0, 97.0],
    [0.0, 2.0, 1.19, 0.34, 0.66, 1.8, -3.0, 89.0],
    [3.0, 1.0, 1.32, 0.57, 1.0, 1.3, -1.0, 99.0],
    [2.0, 1.0, 1.29, 0.49, 0.92, 1.4, -1.0, 96.0],
    [1.0, 1.0, 1.14, 0.26, 0.52, 2.2, -5.0, 85.0],
];
const KBPT_2_3: &[KbpTableRow] = &[
    [0.0, 0.0, 0.55, 0.21, 0.46, 1.2, -5.0, 87.0],
    [4.0, 4.0, 0.63, 0.42, 0.84, 0.75, -2.0, 99.0],
    [2.0, 4.0, 0.615, 0.37, 0.72, 0.85, -3.0, 97.0],
    [0.0, 4.0, 0.55, 0.21, 0.46, 1.2, -5.0, 87.0],
    [3.0, 3.0, 0.615, 0.37, 0.68, 0.9, -3.0, 97.0],
    [6.0, 2.0, 0.63, 0.42, 0.84, 0.75, -2.0, 99.0],
    [5.0, 2.0, 0.625, 0.41, 0.78, 0.8, -2.0, 99.0],
    [4.0, 2.0, 0.61, 0.35, 0.68, 0.9, -3.0, 96.0],
    [2.0, 2.0, 0.515, 0.14, 0.33, 1.55, -9.0, 81.0],
];
const KBPT_3_4: &[KbpTableRow] = &[
    [6.0, 3.0, 0.389, 0.25, 0.56, 0.7, -5.0, 95.0],
    [5.0, 3.0, 0.375, 0.21, 0.47, 0.8, -6.0, 92.0],
    [4.0, 3.0, 0.351, 0.14, 0.35, 1.0, -9.0, 86.0],
    [6.0, 2.0, 0.362, 0.16, 0.45, 0.8, -4.0, 88.0],
    [5.0, 2.0, 0.330, 0.092, 0.28, 1.2, -13.0, 81.0],
    [4.0, 2.0, 0.281, 0.046, 0.16, 1.8, -23.0, 69.0],
];
const KBPT_1_1: &[KbpTableRow] = &[
    [3.0, 2.0, 1.09, 0.31, 0.55, 2.0, -2.0, 99.0],
    [2.0, 2.0, 1.07, 0.27, 0.49, 2.2, -3.0, 97.0],
    [1.0, 2.0, 1.02, 0.21, 0.36, 2.8, -6.0, 92.0],
    [0.0, 2.0, 0.80, 0.064, 0.17, 4.8, -16.0, 72.0],
    [4.0, 1.0, 1.08, 0.28, 0.54, 2.0, -2.0, 98.0],
    [3.0, 1.0, 1.06, 0.25, 0.46, 2.3, -4.0, 96.0],
    [2.0, 1.0, 0.99, 0.17, 0.30, 3.3, -10.0, 90.0],
];
const KBPT_4_5: &[KbpTableRow] = &[
    [0.0, 0.0, 0.22, 0.061, 0.22, 1.0, -15.0, 74.0],
    [6.0, 5.0, 0.28, 0.21, 0.47, 0.6, -7.0, 93.0],
    [5.0, 5.0, 0.27, 0.17, 0.39, 0.7, -9.0, 90.0],
    [4.0, 5.0, 0.25, 0.10, 0.31, 0.8, -10.0, 83.0],
    [3.0, 5.0, 0.23, 0.065, 0.25, 0.9, -11.0, 76.0],
];
const KBPT_3_2: &[KbpTableRow] = &[[5.0, 5.0, 0.208, 0.030, 0.072, 2.9, -47.0, 77.0]];
const KBPT_5_4: &[KbpTableRow] = &[
    [10.0, 6.0, 0.163, 0.068, 0.16, 1.0, -19.0, 85.0],
    [8.0, 6.0, 0.146, 0.039, 0.11, 1.3, -29.0, 76.0],
];

struct KbpTableMeta {
    table: &'static [KbpTableRow],
    gap_open_max: i32,
    gap_extend_max: i32,
    round_down: bool,
}

/// NCBI `kSmallFloat` (`blast_stat.c:4049`): smallest positive e-value
/// used as a clamp in `BlastKarlinEtoS_simple` to prevent FP exceptions.
pub const K_SMALL_FLOAT: f64 = 1.0e-297;

/// NCBI `NUM_FRAMES` (`blast_def.h:88`): number of translation frames
/// for protein searches (6 = ±3).
pub const NUM_FRAMES: usize = 6;

/// NCBI `CODON_LENGTH` (`blast_def.h:63`): nucleotides per amino acid.
pub const CODON_LENGTH: usize = 3;

/// NCBI `DEFAULT_LONGEST_INTRON` (`blast_def.h:78`): default longest
/// intron length (nt) for linking spliced HSPs in translated searches.
pub const DEFAULT_LONGEST_INTRON: usize = 122;

/// NCBI `COMPRESSION_RATIO` (`blast_def.h:83`): NCBI2na packs 4 bases
/// per byte.
pub const COMPRESSION_RATIO: usize = 4;

/// NCBI `NUM_STRANDS` (`blast_def.h:93`): plus / minus.
pub const NUM_STRANDS: usize = 2;

/// NCBI `GENCODE_STRLEN` (`blast_def.h:98`): fixed 64 codons in the
/// genetic-code table (all unique `NNN` triplets in NCBIstdaa).
pub const GENCODE_STRLEN: usize = 64;

/// Port of NCBI `BLAST_SCORE_MIN` (`blast_stat.h:121`): minimum allowed
/// score (one-letter comparison). NCBI defines it as `INT2_MIN`
/// (-32768); we keep the value literally so downstream "impossible
/// score" sentinels compare identically.
pub const BLAST_SCORE_MIN: i32 = i16::MIN as i32;

/// Port of NCBI `BLAST_SCORE_MAX` (`blast_stat.h:122`): maximum allowed
/// score (one-letter comparison) = `INT2_MAX` (32767).
pub const BLAST_SCORE_MAX: i32 = i16::MAX as i32;

/// NCBI `BLAST_MATRIX_NOMINAL` (`blast_stat.h:53`): matrix quality flag
/// for "acceptable values, not recommended".
pub const BLAST_MATRIX_NOMINAL: i32 = 0;
/// NCBI `BLAST_MATRIX_PREFERRED` (`blast_stat.h:54`): preferred values.
pub const BLAST_MATRIX_PREFERRED: i32 = 1;
/// NCBI `BLAST_MATRIX_BEST` (`blast_stat.h:55`): best value (one per matrix).
pub const BLAST_MATRIX_BEST: i32 = 2;

/// `DBSEQ_CHUNK_OVERLAP` (`blast_hits.h:192`): overlap between adjacent
/// chunks when a subject is split for parallel/restarted scanning.
pub const DBSEQ_CHUNK_OVERLAP: usize = 100;

/// `PV_ARRAY_BYTES` (`blast_lookup.h:42`): bytes in the PV-array native word.
pub const PV_ARRAY_BYTES: usize = 4;
/// `PV_ARRAY_BTS` (`blast_lookup.h:43`): bits-to-shift from lookup
/// index to pv-array index (word size 32 → shift 5).
pub const PV_ARRAY_BTS: usize = 5;
/// `PV_ARRAY_MASK` (`blast_lookup.h:44`): mask off low 5 bits to get
/// the bit-within-word offset.
pub const PV_ARRAY_MASK: u32 = 31;

/// `BLAST_SEQSRC_MINGAP` (`blast_seqsrc.h:203`): minimal gap allowed in
/// range list for a sequence source chunk reader.
pub const BLAST_SEQSRC_MINGAP: usize = 1024;
/// `BLAST_SEQSRC_OVERHANG` (`blast_seqsrc.h:204`): extension applied to
/// each new range added to the seq-src range list.
pub const BLAST_SEQSRC_OVERHANG: usize = 1024;
/// `BLAST_SEQSRC_MINLENGTH` (`blast_seqsrc.h:205`): minimum sequence
/// length a seq-src returns.
pub const BLAST_SEQSRC_MINLENGTH: usize = 10;

/// `BLAST_SEQSRC_EXCLUDED` (`blast_seqsrc.h:290`): seq-src status code —
/// sequence excluded due to filtering.
pub const BLAST_SEQSRC_EXCLUDED: i32 = -3;
/// `BLAST_SEQSRC_ERROR` (`blast_seqsrc.h:291`): generic error retrieving sequence.
pub const BLAST_SEQSRC_ERROR: i32 = -2;
/// `BLAST_SEQSRC_EOF` (`blast_seqsrc.h:292`): no more sequences available.
pub const BLAST_SEQSRC_EOF: i32 = -1;
/// `BLAST_SEQSRC_SUCCESS` (`blast_seqsrc.h:293`): successful retrieval.
pub const BLAST_SEQSRC_SUCCESS: i32 = 0;

// Program-type bitmask flags (`blast_program.h:48-64`).
/// `PROTEIN_QUERY_MASK`.
pub const PROTEIN_QUERY_MASK: u32 = 1 << 0;
/// `PROTEIN_SUBJECT_MASK`.
pub const PROTEIN_SUBJECT_MASK: u32 = 1 << 1;
/// `NUCLEOTIDE_QUERY_MASK`.
pub const NUCLEOTIDE_QUERY_MASK: u32 = 1 << 2;
/// `NUCLEOTIDE_SUBJECT_MASK`.
pub const NUCLEOTIDE_SUBJECT_MASK: u32 = 1 << 3;
/// `TRANSLATED_QUERY_MASK`.
pub const TRANSLATED_QUERY_MASK: u32 = 1 << 4;
/// `TRANSLATED_SUBJECT_MASK`.
pub const TRANSLATED_SUBJECT_MASK: u32 = 1 << 5;
/// `PSSM_QUERY_MASK`.
pub const PSSM_QUERY_MASK: u32 = 1 << 6;
/// `PSSM_SUBJECT_MASK`.
pub const PSSM_SUBJECT_MASK: u32 = 1 << 7;
/// `PATTERN_QUERY_MASK`.
pub const PATTERN_QUERY_MASK: u32 = 1 << 8;

fn get_kbp_table(reward: i32, penalty: i32) -> Option<KbpTableMeta> {
    match (reward, penalty) {
        (1, -5) => Some(KbpTableMeta {
            table: KBPT_1_5,
            gap_open_max: 3,
            gap_extend_max: 3,
            round_down: false,
        }),
        (1, -4) => Some(KbpTableMeta {
            table: KBPT_1_4,
            gap_open_max: 2,
            gap_extend_max: 2,
            round_down: false,
        }),
        (2, -7) => Some(KbpTableMeta {
            table: KBPT_2_7,
            gap_open_max: 4,
            gap_extend_max: 4,
            round_down: true,
        }),
        (1, -3) => Some(KbpTableMeta {
            table: KBPT_1_3,
            gap_open_max: 2,
            gap_extend_max: 2,
            round_down: false,
        }),
        (2, -5) => Some(KbpTableMeta {
            table: KBPT_2_5,
            gap_open_max: 4,
            gap_extend_max: 4,
            round_down: true,
        }),
        (1, -2) => Some(KbpTableMeta {
            table: KBPT_1_2,
            gap_open_max: 2,
            gap_extend_max: 2,
            round_down: false,
        }),
        (2, -3) => Some(KbpTableMeta {
            table: KBPT_2_3,
            gap_open_max: 6,
            gap_extend_max: 4,
            round_down: true,
        }),
        (3, -4) => Some(KbpTableMeta {
            table: KBPT_3_4,
            gap_open_max: 6,
            gap_extend_max: 3,
            round_down: true,
        }),
        (1, -1) => Some(KbpTableMeta {
            table: KBPT_1_1,
            gap_open_max: 4,
            gap_extend_max: 2,
            round_down: false,
        }),
        (3, -2) => Some(KbpTableMeta {
            table: KBPT_3_2,
            gap_open_max: 5,
            gap_extend_max: 5,
            round_down: false,
        }),
        (4, -5) => Some(KbpTableMeta {
            table: KBPT_4_5,
            gap_open_max: 12,
            gap_extend_max: 8,
            round_down: false,
        }),
        (5, -4) => Some(KbpTableMeta {
            table: KBPT_5_4,
            gap_open_max: 25,
            gap_extend_max: 10,
            round_down: false,
        }),
        _ => None,
    }
}

/// Look up gapped KBP params for nucleotide, matching C's Blast_KarlinBlkNuclGappedCalc exactly.
/// Returns Ok((lambda, k, log_k, h, round_down)) or Err message.
pub fn nucl_gapped_kbp_lookup(
    gap_open: i32,
    gap_extend: i32,
    reward: i32,
    penalty: i32,
    ungapped: &KarlinBlk,
) -> Result<(KarlinBlk, bool), String> {
    let divisor = crate::math::gcd(reward, penalty.abs());
    let (nr, np) = (reward / divisor, penalty / divisor);

    let meta = get_kbp_table(nr, np)
        .ok_or_else(|| format!("Unsupported scores {} {}", reward, penalty))?;

    let round_down = meta.round_down;

    // Split: first row with gap_open==0 && gap_extend==0 is the linear entry
    let (affine, linear) =
        if !meta.table.is_empty() && meta.table[0][0] == 0.0 && meta.table[0][1] == 0.0 {
            (&meta.table[1..], Some(&meta.table[0]))
        } else {
            (meta.table, None)
        };

    // GCD scaling (Rust extension beyond NCBI `Blast_KarlinBlkNuclGappedCalc`,
    // `blast_stat.c:3846`): NCBI's lookup reduces `(reward, penalty)` by
    // their GCD when selecting which table to use but then compares gap
    // costs against the reduced-system table values directly — users
    // must supply gap costs in the reduced units. Rust instead scales
    // the table entries by the divisor so users can pass gap costs in
    // the scaled system (e.g. `(reward=10, penalty=-20, gap_open=30,
    // gap_extend=10)` reduces to KBPT_1_2 with scaled lambda = row[2]/10).
    let (mut go_max, mut ge_max) = (meta.gap_open_max, meta.gap_extend_max);
    let scale = |row: &KbpTableRow, d: i32| -> (f64, f64, f64, f64, f64) {
        let go = row[0] * d as f64;
        let ge = row[1] * d as f64;
        let lam = row[2] / d as f64;
        (go, ge, lam, row[3], row[4])
    };
    if divisor != 1 {
        go_max *= divisor;
        ge_max *= divisor;
    }

    // Linear (non-affine) case
    if gap_open == 0 && gap_extend == 0 {
        if let Some(lin) = linear {
            let (_, _, lam, k, h) = scale(lin, divisor);
            return Ok((
                KarlinBlk {
                    lambda: lam,
                    k,
                    log_k: k.ln(),
                    h,
                    round_down,
                },
                round_down,
            ));
        }
    }

    // Search affine entries
    for row in affine {
        let (go, ge, lam, k, h) = scale(row, divisor);
        if go as i32 == gap_open && ge as i32 == gap_extend {
            return Ok((
                KarlinBlk {
                    lambda: lam,
                    k,
                    log_k: k.ln(),
                    h,
                    round_down,
                },
                round_down,
            ));
        }
    }

    // Fallback: gap costs exceed table maximum → use ungapped params
    if gap_open >= go_max && gap_extend >= ge_max {
        return Ok((ungapped.clone(), round_down));
    }

    Err(format!(
        "Unsupported gap costs {} {} for scores {} {}",
        gap_open, gap_extend, reward, penalty
    ))
}

/// Compute the length adjustment for effective search space (exact C-compatible).
/// Port of BLAST_ComputeLengthAdjustment from blast_stat.c.
/// Returns (length_adjustment, converged).
pub fn compute_length_adjustment_exact(
    k: f64,
    log_k: f64,
    alpha_d_lambda: f64,
    beta: f64,
    query_length: i32,
    db_length: i64,
    db_num_seqs: i32,
) -> (i32, bool) {
    let m = query_length as f64;
    let n = db_length as f64;
    let nn = db_num_seqs as f64;

    // Compute ell_max using quadratic formula
    let a = nn;
    let mb = m * nn + n;
    let c = n * m - m.max(n) / k;
    if c < 0.0 {
        return (0, false);
    }
    let ell_max_init = 2.0 * c / (mb + (mb * mb - 4.0 * a * c).sqrt());

    let mut ell_min = 0.0_f64;
    let mut ell_max = ell_max_init;
    let mut ell_next = 0.0_f64;
    let mut converged = false;

    for i in 1..=20 {
        let ell = ell_next;
        let ss = (m - ell) * (n - nn * ell);
        let ell_bar = alpha_d_lambda * (log_k + ss.ln()) + beta;

        if ell_bar >= ell {
            ell_min = ell;
            if ell_bar - ell_min <= 1.0 {
                converged = true;
                break;
            }
            if ell_min == ell_max {
                break;
            }
        } else {
            ell_max = ell;
        }

        if ell_min <= ell_bar && ell_bar <= ell_max {
            ell_next = ell_bar;
        } else {
            ell_next = if i == 1 {
                ell_max
            } else {
                (ell_min + ell_max) / 2.0
            };
        }
    }

    let mut adj = ell_min as i32;
    if converged {
        let ell = ell_min.ceil();
        if ell <= ell_max {
            let ss = (m - ell) * (n - nn * ell);
            if alpha_d_lambda * (log_k + ss.ln()) + beta >= ell {
                adj = ell as i32;
            }
        }
    }
    (adj, converged)
}

/// Look up alpha and beta for nucleotide gapped alignment.
/// Port of Blast_GetNuclAlphaBeta from blast_stat.c.
pub fn nucl_alpha_beta(
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    ungapped_lambda: f64,
    ungapped_h: f64,
    gapped: bool,
) -> (f64, f64) {
    if !gapped {
        return (
            ungapped_lambda / ungapped_h,
            get_ungapped_beta(reward, penalty),
        );
    }

    // Rust extension beyond NCBI `Blast_GetNuclAlphaBeta`
    // (`blast_stat.c:3965`): scale gap costs and alpha by the GCD
    // divisor so scaled scoring systems like `(10, -20)` match the
    // reduced `(1, -2)` table. NCBI compares against raw table values
    // and forces users to supply reduced-system gap costs.
    let divisor = crate::math::gcd(reward, penalty.abs());
    let (nr, np) = (reward / divisor, penalty / divisor);

    if let Some(meta) = get_kbp_table(nr, np) {
        let (affine, linear) =
            if !meta.table.is_empty() && meta.table[0][0] == 0.0 && meta.table[0][1] == 0.0 {
                (&meta.table[1..], Some(&meta.table[0]))
            } else {
                (meta.table, None)
            };

        if gap_open == 0 && gap_extend == 0 {
            if let Some(lin) = linear {
                let a = lin[5] / divisor as f64;
                let b = lin[6];
                return (a, b);
            }
        }

        for row in affine {
            let go = (row[0] * divisor as f64) as i32;
            let ge = (row[1] * divisor as f64) as i32;
            if go == gap_open && ge == gap_extend {
                let a = row[5] / divisor as f64;
                let b = row[6];
                return (a, b);
            }
        }
    }

    // Fallback: ungapped values. NCBI `Blast_GetNuclAlphaBeta`
    // (`blast_stat.c:4022-4026`): on lookup miss, returns
    // `alpha = Lambda/H` and `beta = s_GetUngappedBeta(reward, penalty)`.
    (
        ungapped_lambda / ungapped_h,
        get_ungapped_beta(reward, penalty),
    )
}

/// Port of NCBI `s_GetUngappedBeta` (`blast_stat.c:3955`): `(1,-1)` and
/// `(2,-3)` scoring systems use `beta = -2`; every other combination
/// uses `beta = 0`.
fn get_ungapped_beta(reward: i32, penalty: i32) -> f64 {
    if (reward == 1 && penalty == -1) || (reward == 2 && penalty == -3) {
        -2.0
    } else {
        0.0
    }
}

// ---------------------------------------------------------------------------
// Ungapped Karlin-Altschul parameter computation (exact C-compatible)
// Port of Blast_ScoreBlkKbpUngappedCalc and sub-functions from blast_stat.c
// ---------------------------------------------------------------------------

/// Score frequency distribution (internal to KBP computation). Rust
/// analog of NCBI `Blast_ScoreFreq` (`blast_stat.h:128`); `probs`
/// corresponds to NCBI's `sprob` array (shifted by `score_min`).
struct SfDist {
    score_min: i32,
    score_max: i32,
    obs_min: i32,
    obs_max: i32,
    score_avg: f64,
    probs: Vec<f64>, // indexed by [score - score_min]
}

impl SfDist {
    fn new(lo: i32, hi: i32) -> Self {
        let n = (hi - lo + 1) as usize;
        SfDist {
            score_min: lo,
            score_max: hi,
            obs_min: 0,
            obs_max: 0,
            score_avg: 0.0,
            probs: vec![0.0; n],
        }
    }
    #[inline]
    fn p(&self, s: i32) -> f64 {
        self.probs[(s - self.score_min) as usize]
    }
    #[inline]
    fn p_mut(&mut self, s: i32) -> &mut f64 {
        &mut self.probs[(s - self.score_min) as usize]
    }
}

/// Newton-Raphson solver for Lambda in x = exp(-lambda) space.
/// Port of NlmKarlinLambdaNR from blast_stat.c.
/// NCBI `BLAST_KARLIN_K_SUMLIMIT_DEFAULT` (`blast_stat.c:64`): inner-sum
/// convergence threshold in the K solver.
const BLAST_KARLIN_K_SUMLIMIT_DEFAULT: f64 = 0.0001;
/// NCBI `BLAST_KARLIN_LAMBDA_ACCURACY_DEFAULT` (`blast_stat.c:66`):
/// accuracy target for lambda solve.
const BLAST_KARLIN_LAMBDA_ACCURACY_DEFAULT: f64 = 1.0e-5;
/// NCBI `BLAST_KARLIN_LAMBDA_ITER_DEFAULT` (`blast_stat.c:68`): number
/// of bisection iterations after the Newton phase (20 + 17 = 37 total).
const BLAST_KARLIN_LAMBDA_ITER_DEFAULT: i32 = 17;
/// NCBI `BLAST_KARLIN_K_ITER_MAX` (`blast_stat.c:72`): upper limit on
/// iterations for the K solver.
const BLAST_KARLIN_K_ITER_MAX: usize = 100;

/// Port of NCBI `Blast_KarlinLambdaNR` (`blast_stat.c:2567`): combined
/// Newton/bisection solver for the Karlin-Altschul lambda parameter
/// over `x = exp(-lambda)`.
fn solve_lambda(sfp: &SfDist, d: i32, low: i32, high: i32, lambda0: f64) -> f64 {
    let x0 = (-lambda0).exp();
    let mut x = if x0 > 0.0 && x0 < 1.0 { x0 } else { 0.5 };
    let mut a = 0.0_f64;
    let mut b = 1.0_f64;
    let mut f = 4.0_f64;
    let mut is_newton = false;
    // NCBI's `Blast_KarlinLambdaNR` calls `NlmKarlinLambdaNR(... 20, 20 +
    // BLAST_KARLIN_LAMBDA_ITER_DEFAULT, ...)` — passing `itmax=20` and
    // `maxNewton=37`. In NCBI's implementation `itmax` is the loop bound and
    // `maxNewton` is the Newton-fallback threshold; with maxNewton > itmax
    // the threshold never fires within the loop. So the effective behavior
    // is "20 iterations, Newton always tried." Match NCBI exactly.
    let max_iter = 20;
    let max_newton = 20 + BLAST_KARLIN_LAMBDA_ITER_DEFAULT; // 37, threshold beyond max_iter
    let tolx = BLAST_KARLIN_LAMBDA_ACCURACY_DEFAULT;

    for k in 0..max_iter {
        let fold = f;
        let was_newton = is_newton;
        is_newton = false;

        // Horner evaluation of polynomial and derivative.
        let mut g = 0.0_f64;
        f = sfp.p(low);
        let mut i = low + d;
        while i < 0 {
            g = x * g + f;
            f = f * x + sfp.p(i);
            i += d;
        }
        g = x * g + f;
        f = f * x + sfp.p(0) - 1.0;
        i = d;
        while i <= high {
            g = x * g + f;
            f = f * x + sfp.p(i);
            i += d;
        }

        if f > 0.0 {
            a = x;
        } else if f < 0.0 {
            b = x;
        } else {
            break;
        }
        if b - a < 2.0 * a * (1.0 - b) * tolx {
            x = (a + b) / 2.0;
            break;
        }

        if k >= max_newton || (was_newton && f.abs() > 0.9 * fold.abs()) || g >= 0.0 {
            x = (a + b) / 2.0;
        } else {
            let p = -f / g;
            let y = x + p;
            if y <= a || y >= b {
                x = (a + b) / 2.0;
            } else {
                is_newton = true;
                x = y;
                if p.abs() < tolx * x * (1.0 - x) {
                    break;
                }
            }
        }
    }
    -x.ln() / d as f64
}

/// NCBI `BLAST_KARLIN_LAMBDA0_DEFAULT` (`blast_stat.c:70`): initial
/// guess for the Newton-Raphson lambda solver.
const BLAST_KARLIN_LAMBDA0_DEFAULT: f64 = 0.5;

/// Compute Lambda from score frequency distribution. Port of NCBI
/// `Blast_KarlinLambdaNR` (`blast_stat.c:2567`) — score-GCD reduction
/// followed by Newton-Raphson solve on `x = exp(-lambda)`.
fn compute_lambda(sfp: &SfDist) -> f64 {
    if sfp.score_avg >= 0.0 {
        return -1.0;
    }
    let low = sfp.obs_min;
    let high = sfp.obs_max;
    let mut d = -low;
    for i in 1..=(high - low) {
        if d <= 1 {
            break;
        }
        if sfp.p(low + i) != 0.0 {
            d = crate::math::gcd(d, i);
        }
    }
    solve_lambda(sfp, d, low, high, BLAST_KARLIN_LAMBDA0_DEFAULT)
}

/// Compute H (relative entropy) from score frequencies and Lambda.
fn compute_h(sfp: &SfDist, lambda: f64) -> f64 {
    if lambda < 0.0 {
        return -1.0;
    }
    let etonlam = (-lambda).exp();
    let mut sum = sfp.obs_min as f64 * sfp.p(sfp.obs_min);
    for s in (sfp.obs_min + 1)..=sfp.obs_max {
        sum = s as f64 * sfp.p(s) + etonlam * sum;
    }
    let scale = crate::math::powi(etonlam, sfp.obs_max);
    if scale > 0.0 {
        lambda * sum / scale
    } else {
        lambda * (lambda * sfp.obs_max as f64 + sum.ln()).exp()
    }
}

/// Compute K from Lambda and H using DP algorithm.
fn compute_k(sfp: &SfDist, lambda: f64, h: f64) -> f64 {
    if lambda <= 0.0 || h <= 0.0 || sfp.score_avg >= 0.0 {
        return -1.0;
    }
    let olow = sfp.obs_min;
    let ohigh = sfp.obs_max;
    let mut divisor = -olow;
    for i in 1..=(ohigh - olow) {
        if divisor <= 1 {
            break;
        }
        if sfp.p(olow + i) != 0.0 {
            divisor = crate::math::gcd(divisor, i);
        }
    }

    let high = ohigh / divisor;
    let low = olow / divisor;
    let lam_s = lambda * divisor as f64;
    let range = high - low;
    let ftcf = h / lam_s;
    let eml = (-lam_s).exp();

    // Special cases
    if low == -1 && high == 1 {
        let pl = sfp.p(olow);
        let ph = sfp.p(ohigh);
        return (pl - ph) * (pl - ph) / pl;
    }
    if low == -1 || high == 1 {
        let mut f = ftcf;
        if high != 1 {
            let sa = sfp.score_avg / divisor as f64;
            f = (sa * sa) / f;
        }
        return f * (1.0 - eml);
    }

    // DP
    let prob_at = |i: i32| -> f64 {
        let s = olow + i;
        if s >= sfp.score_min && s <= sfp.score_max {
            sfp.p(s)
        } else {
            0.0
        }
    };
    let ru = range as usize;
    let max_iter = BLAST_KARLIN_K_ITER_MAX;
    let mut asp = vec![0.0_f64; max_iter * ru + 1];
    asp[0] = 1.0;
    let mut outer_sum = 0.0_f64;
    let mut low_as = 0i32;
    let mut high_as = 0i32;
    let mut inner_sum = 1.0_f64;
    let mut oldsum;

    // NCBI's `BlastKarlinLHtoK` (`blast_stat.c:2346`) uses an idiomatic C
    // for-loop:
    //   for (iterCounter = 0;
    //        (iterCounter < iterlimit) && (innerSum > sumlimit);
    //        outerSum += innerSum /= ++iterCounter)
    // The increment expression DIVIDES innerSum in-place by the new
    // iterCounter BEFORE adding to outerSum. So the loop CONDITION at
    // iter K (K > 0) compares `innerSum_{K-1} / K > sumlimit` — i.e. the
    // divided value carries into the next iter's check. That makes the
    // effective sumlimit grow as K (sumlimit * K threshold). We mirror
    // this by mutating `inner_sum` after the body.
    let mut iter_counter = 0i32;
    while iter_counter < max_iter as i32 && inner_sum > BLAST_KARLIN_K_SUMLIMIT_DEFAULT {
        let mut first = range;
        let mut last = range;
        low_as += low;
        high_as += high;
        let span = (high_as - low_as) as isize;
        let mut pp = span;
        while pp >= 0 {
            let p1s = pp - first as isize;
            let p1e = pp - last as isize;
            let mut isum = 0.0;
            let mut p1 = p1s;
            let mut p2 = first;
            while p1 >= p1e {
                if p1 >= 0 {
                    isum += asp[p1 as usize] * prob_at(p2);
                }
                p1 -= 1;
                p2 += 1;
            }
            asp[pp as usize] = isum;
            if first > 0 {
                first -= 1;
            }
            if pp <= range as isize {
                last -= 1;
            }
            pp -= 1;
        }
        pp += 1;
        inner_sum = asp[pp as usize];
        let mut i = low_as + 1;
        while i < 0 {
            pp += 1;
            inner_sum = asp[pp as usize] + inner_sum * eml;
            i += 1;
        }
        inner_sum *= eml;
        while i <= high_as {
            pp += 1;
            inner_sum += asp[pp as usize];
            i += 1;
        }
        oldsum = inner_sum;
        let _ = oldsum;
        // Match NCBI's `outerSum += innerSum /= ++iterCounter`:
        // first increment iterCounter, then divide innerSum in place,
        // then add to outerSum.
        iter_counter += 1;
        inner_sum /= iter_counter as f64;
        outer_sum += inner_sum;
    }

    -(-2.0 * outer_sum).exp() / (ftcf * crate::math::expm1(-lam_s))
}

/// Context info for ungapped KBP computation.
pub struct UngappedKbpContext {
    pub query_offset: i32,
    pub query_length: i32,
    pub is_valid: bool,
}

/// Compute ungapped KBP for all contexts. Returns per-context `Option<KarlinBlk>`.
/// Port of Blast_ScoreBlkKbpUngappedCalc from blast_stat.c (nucleotide path —
/// std composition is fixed at 0.25 for A/C/G/T).
pub fn ungapped_kbp_calc(
    query: &[u8],
    contexts: &[UngappedKbpContext],
    loscore: i32,
    hiscore: i32,
    alphabet_size: usize,
    ambiguous: &[u8],
    matrix: &dyn Fn(usize, usize) -> i32,
) -> Vec<Option<KarlinBlk>> {
    // Standard composition: 0.25 each for A/C/G/T (indices 0-3)
    let mut std_freq = vec![0.0f64; alphabet_size];
    for i in 0..4.min(alphabet_size) {
        std_freq[i] = 0.25;
    }
    ungapped_kbp_calc_with_std(
        query,
        contexts,
        loscore,
        hiscore,
        alphabet_size,
        ambiguous,
        &std_freq,
        matrix,
    )
}

/// Like [`ungapped_kbp_calc`] but takes an explicit `std_freq` array.
/// Used by protein paths (blastx/tblastn/tblastx) where the standard
/// background is BLOSUM62's Robinson & Robinson amino acid frequencies
/// rather than the uniform 0.25 nucleotide composition.
///
/// Mirrors `Blast_ScoreBlkKbpUngappedCalc` (`blast_stat.c:2737`): per-context
/// `Blast_ResFreqString` → `BlastScoreFreqCalc` → `Blast_KarlinBlkUngappedCalc`.
pub fn ungapped_kbp_calc_with_std(
    query: &[u8],
    contexts: &[UngappedKbpContext],
    loscore: i32,
    hiscore: i32,
    alphabet_size: usize,
    ambiguous: &[u8],
    std_freq: &[f64],
    matrix: &dyn Fn(usize, usize) -> i32,
) -> Vec<Option<KarlinBlk>> {
    let std_freq = std_freq.to_vec();
    let mut results = Vec::with_capacity(contexts.len());
    for ctx in contexts {
        if !ctx.is_valid || ctx.query_length <= 0 {
            results.push(None);
            continue;
        }
        let off = ctx.query_offset as usize;
        let len = ctx.query_length as usize;
        let buf = &query[off..off + len];

        // Count residue composition. The 4-bit mask was needed only when the
        // caller was nucleotide-only (BLASTNA bytes are 0..15); now the
        // function is also used for protein (NCBIstdaa, bytes 0..27).
        // Removing the mask is a no-op for nucleotide (bytes already <16) and
        // unblocks protein. The `idx < alphabet_size` guard still bounds.
        let mut counts = vec![0i32; alphabet_size];
        for &b in buf {
            let idx = b as usize;
            if idx < alphabet_size {
                counts[idx] += 1;
            }
        }
        for &a in ambiguous {
            if (a as usize) < alphabet_size {
                counts[a as usize] = 0;
            }
        }
        let sum: f64 = counts.iter().map(|&c| c as f64).sum();
        let mut qfreq = vec![0.0f64; alphabet_size];
        if sum > 0.0 {
            for i in 0..alphabet_size {
                qfreq[i] = counts[i] as f64 / sum;
            }
        }

        // Compute score frequency
        let mut sfp = SfDist::new(loscore, hiscore);
        for i in 0..alphabet_size {
            for j in 0..alphabet_size {
                let s = matrix(i, j);
                if s >= loscore && s <= hiscore {
                    *sfp.p_mut(s) += qfreq[i] * std_freq[j];
                }
            }
        }
        // Find obs range and normalize
        let mut obs_min = i32::MAX;
        let mut obs_max = i32::MIN;
        let mut psum = 0.0;
        for s in loscore..=hiscore {
            if sfp.p(s) > 0.0 {
                psum += sfp.p(s);
                if s < obs_min {
                    obs_min = s;
                }
                obs_max = s;
            }
        }
        sfp.obs_min = obs_min;
        sfp.obs_max = obs_max;
        let mut savg = 0.0;
        if psum > 0.0 {
            for s in obs_min..=obs_max {
                *sfp.p_mut(s) /= psum;
                savg += s as f64 * sfp.p(s);
            }
        }
        sfp.score_avg = savg;

        // Compute Lambda, H, K
        let lambda = compute_lambda(&sfp);
        if lambda < 0.0 {
            results.push(None);
            continue;
        }
        let h = compute_h(&sfp, lambda);
        if h < 0.0 {
            results.push(None);
            continue;
        }
        let k = compute_k(&sfp, lambda, h);
        if k < 0.0 {
            results.push(None);
            continue;
        }
        results.push(Some(KarlinBlk {
            lambda,
            k,
            log_k: k.ln(),
            h,
            round_down: false,
        }));
    }
    results
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gap_decay_divisor_matches_ncbi_formula() {
        // NCBI `BLAST_GapDecayDivisor(decayrate, nsegs)` =
        // `(1 - decayrate) * decayrate^(nsegs - 1)`. Spot-check a few
        // typical values.
        let eps = 1e-12;
        // nsegs=1 → (1-r) * r^0 = 1-r.
        assert!((gap_decay_divisor(0.1, 1) - 0.9).abs() < eps);
        assert!((gap_decay_divisor(0.5, 1) - 0.5).abs() < eps);
        // nsegs=2 → (1-r) * r.
        assert!((gap_decay_divisor(0.1, 2) - 0.09).abs() < eps);
        assert!((gap_decay_divisor(0.5, 2) - 0.25).abs() < eps);
        // nsegs=3 → (1-r) * r^2.
        assert!((gap_decay_divisor(0.1, 3) - 0.009).abs() < eps);
    }

    #[test]
    fn test_blast_cutoffs_without_decay_matches_evalue_to_raw() {
        // With dodecay=false the adjustment is skipped and the cutoff is
        // just `EtoS(e)` (clamped to `s_floor`).
        let kbp = KarlinBlk {
            lambda: 0.625,
            k: 0.41,
            log_k: 0.41_f64.ln(),
            h: 0.78,
            round_down: false,
        };
        let searchsp = 1.0e9;
        let e_in = 10.0;
        let expected_s = kbp.evalue_to_raw(e_in, searchsp);
        let (s, _e_out) = blast_cutoffs(1, e_in, &kbp, searchsp, false, 0.1);
        assert_eq!(s, expected_s);
    }

    #[test]
    fn test_blast_cutoffs_with_decay_tightens_cutoff() {
        // With dodecay=true and gap_decay_rate=0.1, e is multiplied by
        // `gap_decay_divisor(0.1, 1) = 0.9` before conversion, giving a
        // slightly higher raw cutoff score than the no-decay path.
        let kbp = KarlinBlk {
            lambda: 0.625,
            k: 0.41,
            log_k: 0.41_f64.ln(),
            h: 0.78,
            round_down: false,
        };
        let searchsp = 1.0e9;
        let e_in = 10.0;
        let (s_no_decay, _) = blast_cutoffs(1, e_in, &kbp, searchsp, false, 0.1);
        let (s_decay, _) = blast_cutoffs(1, e_in, &kbp, searchsp, true, 0.1);
        assert!(
            s_decay >= s_no_decay,
            "decay={s_decay} vs no_decay={s_no_decay}"
        );
        // Independent calculation of the expected decayed cutoff.
        let expected = kbp.evalue_to_raw(0.9 * e_in, searchsp);
        assert_eq!(s_decay, expected);
    }

    #[test]
    fn test_blast_cutoffs_respects_floor() {
        // When the user supplies a floor larger than the stats-derived
        // cutoff, the floor wins and the recomputed e-value corresponds
        // to that floor (with the decay adjustment folded back in).
        let kbp = KarlinBlk {
            lambda: 0.625,
            k: 0.41,
            log_k: 0.41_f64.ln(),
            h: 0.78,
            round_down: false,
        };
        let searchsp = 1.0e9;
        let e_in = 10.0;
        let floor = 1000;
        let (s, e_out) = blast_cutoffs(floor, e_in, &kbp, searchsp, true, 0.1);
        assert_eq!(s, floor);
        // e_in did not win → recompute e from s=floor, then divide by
        // the decay divisor (0.9) to get the reported e-value.
        let expected_e = kbp.raw_to_evalue(floor, searchsp) / 0.9;
        assert!(
            (e_out - expected_e).abs() < 1e-14,
            "e_out={e_out} expected={expected_e}"
        );
    }

    #[test]
    fn test_blast_cutoffs_invalid_kbp_returns_one() {
        // Lambda/K/H < 0 → NCBI returns immediately with S=1.
        let kbp = KarlinBlk {
            lambda: -1.0,
            k: 0.41,
            log_k: 0.0,
            h: 0.78,
            round_down: false,
        };
        let (s, e_out) = blast_cutoffs(5, 10.0, &kbp, 1.0e9, true, 0.1);
        assert_eq!(s, 1);
        assert_eq!(e_out, 10.0);
    }

    #[test]
    fn test_sum_p_r1_formula() {
        // r=1 closed form: p = 1 - exp(-exp(-s)).
        // For s=0, exp(0)=1, 1-exp(-1) ≈ 0.632120.
        let p = sum_p(1, 0.0).unwrap();
        assert!((p - (1.0 - (-1.0_f64).exp())).abs() < 1e-12);
        // For s=5, small p.
        let p = sum_p(1, 5.0).unwrap();
        assert!(p > 0.0 && p < 0.01);
        // For s=-2, close to 1.
        let p = sum_p(1, -2.0).unwrap();
        assert!(p > 0.99);
    }

    #[test]
    fn test_sum_p_r2_tail_formula_against_hand_calc() {
        // r=2, s=6 is in the tail regime (s >= r^2 + r - 1 = 5).
        // NCBI formula: r * exp((r-1)*ln(s) - s - 2*ln(r!)).
        // = 2 * exp(ln(6) - 6 - 2*ln(2)) = 2 * (6/4) * exp(-6) = 3*exp(-6).
        let p = sum_p(2, 6.0).unwrap();
        let expected = 3.0 * (-6.0_f64).exp();
        assert!(
            (p - expected).abs() < 1e-12,
            "sum_p(2, 6) = {p}, expected {expected}"
        );
    }

    #[test]
    fn test_sum_p_r3_tail_formula_against_hand_calc() {
        // r=3, s=11 is in the tail regime (s >= r^2 + r - 1 = 11).
        // NCBI formula: r * exp((r-1)*ln(s) - s - 2*ln(r!)).
        // = 3 * exp(2*ln(11) - 11 - 2*ln(6)) = 3 * (121/36) * exp(-11).
        let p = sum_p(3, 11.0).unwrap();
        let expected = 3.0 * (121.0 / 36.0) * (-11.0_f64).exp();
        assert!(
            (p - expected).abs() < 1e-12,
            "sum_p(3, 11) = {p}, expected {expected}"
        );
    }

    #[test]
    fn test_sum_p_r2_table_interpolation() {
        // r=2 should produce a P-value in (0,1] for plausible s.
        for s in [-2.0, 0.0, 2.0, 5.0, 10.0, 15.0] {
            let p = sum_p(2, s).unwrap();
            assert!(
                (0.0..=1.0).contains(&p),
                "sum_p(2, {}) = {} not in [0,1]",
                s,
                p
            );
        }
        // At very negative s, should saturate to 1.
        assert_eq!(sum_p(2, -10.0).unwrap(), 1.0);
    }

    #[test]
    fn test_sum_p_r_gte_5_uses_romberg_branch() {
        // r >= 5 goes through `sum_p_calc` (Romberg integration). For
        // very negative s it saturates to 1 via the early-out; for
        // moderate s it produces a valid P-value in (0, 1].
        assert_eq!(sum_p(5, -20.0).unwrap(), 1.0);
        let p = sum_p(10, 5.0).unwrap();
        assert!(
            (0.0..=1.0).contains(&p),
            "sum_p(10, 5.0) = {p} not in [0,1]"
        );
        // r=1 via sum_p should still match sum_p_calc(r=1).
        let s = 3.0;
        assert!((sum_p(1, s).unwrap() - sum_p_calc(1, s)).abs() < 1e-14);
    }

    #[test]
    fn test_sum_p_calc_r1_matches_closed_form() {
        // sum_p_calc(1, s) has two branches (s>8 and s<=8); they are
        // both small-x approximations of the same Karlin-Altschul
        // P-value and agree to ~6 significant figures at the crossover.
        let a = sum_p_calc(1, 8.0);
        let b = sum_p_calc(1, 8.0 + 1e-12);
        assert!((a - b).abs() < 1e-6, "a={a} b={b}");
    }

    #[test]
    fn test_sum_e_num_gte_5_now_uses_romberg() {
        // Since sum_p handles r >= 5 via Romberg, the sum-E wrappers
        // now return Some for num >= 5 too.
        let e = small_gap_sum_e(40, 5, 30.0, 100, 1000, 1.0e9, 1.0);
        assert!(e.is_some(), "expected Some, got None");
    }

    #[test]
    fn test_karlin_p_to_e_basic() {
        // P = 0 → E = 0.
        assert_eq!(karlin_p_to_e(0.0), 0.0);
        // P = 1 → E = +infinity (NCBI caps at INT4_MAX downstream).
        assert_eq!(karlin_p_to_e(1.0), f64::INFINITY);
        // P = 0.5 → E = -ln(0.5) = ln(2) ≈ 0.6931.
        assert!((karlin_p_to_e(0.5) - std::f64::consts::LN_2).abs() < 1e-12);
        // Bad input returns sentinel.
        assert_eq!(karlin_p_to_e(-0.1), i32::MIN as f64);
        assert_eq!(karlin_p_to_e(1.5), i32::MIN as f64);
    }

    #[test]
    fn test_sum_e_num1_reduces_to_searchsp_exp() {
        // For num=1 all three sum-E variants collapse to
        // `searchsp_eff * exp(-xsum) / weight_divisor`.
        let xsum: f64 = 10.0;
        let searchsp: f64 = 1.0e9;
        let w: f64 = 0.9;
        let expected = searchsp * (-xsum).exp() / w;
        for got in [
            small_gap_sum_e(40, 1, xsum, 100, 1000, searchsp, w).unwrap(),
            uneven_gap_sum_e(40, 4000, 1, xsum, 100, 1000, searchsp, w).unwrap(),
            large_gap_sum_e(1, xsum, 100, 1000, searchsp, w).unwrap(),
        ] {
            assert!(
                (got - expected).abs() < 1e-6,
                "got={got} expected={expected}"
            );
        }
    }

    #[test]
    fn test_sum_e_caps_at_int4_max() {
        // Zero weight_divisor saturates to INT4_MAX.
        let capped = small_gap_sum_e(40, 1, 0.0, 100, 1000, 1.0e15, 0.0).unwrap();
        assert_eq!(capped, i32::MAX as f64);
    }

    #[test]
    fn test_large_gap_sum_e_num2_matches_hand_calc() {
        // Reconstruct the NCBI formula (`blast_stat.c:4557-4567`) for
        // num=2 step-by-step and verify `large_gap_sum_e` agrees.
        let q: i32 = 100;
        let s: i32 = 1000;
        let searchsp: f64 = 1.0e9;
        let w: f64 = 0.5;
        let raw_xsum: f64 = 120.0;
        let num: u32 = 2;

        // Manual: adjusted = xsum - num*ln(q*s) + ln_fact(num).
        let adjusted = raw_xsum - num as f64 * ((q as f64) * (s as f64)).ln()
            + crate::math::ln_factorial(num as i32);
        let p = sum_p(num, adjusted).unwrap();
        let expected = karlin_p_to_e(p) * (searchsp / (q as f64 * s as f64)) / w;

        let got = large_gap_sum_e(num, raw_xsum, q, s, searchsp, w).unwrap();
        assert!(
            (got - expected).abs() < expected.abs() * 1e-12 + 1e-12,
            "got={got} expected={expected}"
        );
    }

    #[test]
    fn test_small_gap_sum_e_num2_matches_hand_calc() {
        // Hand-verify NCBI `BLAST_SmallGapSumE` (`blast_stat.c:4440-4456`)
        // at num=2. Double-adjust xsum, fold in ln_factorial(2), then
        // `karlin_p_to_e(sum_p) * (searchsp / (q*s)) / w`.
        let q: i32 = 100;
        let s: i32 = 1000;
        let searchsp: f64 = 1.0e9;
        let w: f64 = 0.9;
        let raw_xsum: f64 = 120.0;
        let num: u32 = 2;
        let starting_points: i32 = 40;

        let pair_search_space = s as f64 * q as f64;
        let adjusted = raw_xsum
            - pair_search_space.ln()
            - 2.0 * (num - 1) as f64 * (starting_points as f64).ln()
            - crate::math::ln_factorial(num as i32);
        let p = sum_p(num, adjusted).unwrap();
        let expected = karlin_p_to_e(p) * (searchsp / pair_search_space) / w;

        let got = small_gap_sum_e(starting_points, num, raw_xsum, q, s, searchsp, w).unwrap();
        assert!(
            (got - expected).abs() < expected.abs() * 1e-12 + 1e-12,
            "got={got} expected={expected}"
        );
    }

    #[test]
    fn test_sum_p_r0_is_zero() {
        assert_eq!(sum_p(0, 0.0), Some(0.0));
        assert_eq!(sum_p(0, 5.0), Some(0.0));
    }

    #[test]
    fn test_evalue_to_raw_roundtrips_raw_to_evalue() {
        // For any positive raw score, `raw_to_evalue` then `evalue_to_raw`
        // should round-trip within ±1 (because `evalue_to_raw` rounds up
        // and `raw_to_evalue` is smooth).
        let kbp = KarlinBlk {
            lambda: 0.625,
            k: 0.41,
            log_k: 0.41_f64.ln(),
            h: 0.78,
            round_down: false,
        };
        let search_space = 1.0e9;
        for raw in [10, 25, 50, 100, 200] {
            let e = kbp.raw_to_evalue(raw, search_space);
            let back = kbp.evalue_to_raw(e, search_space);
            assert!(
                (back - raw).abs() <= 1,
                "raw={raw} e={e} back={back} (should be within 1)"
            );
        }
    }

    #[test]
    fn test_round_down_affects_evalue_not_bit_score() {
        // NCBI `blast_hits.c:1864-1869` applies `score &= ~1` before
        // E-value calculation when `sbp->round_down` is set, but
        // `Blast_HSPListGetBitScores` (`blast_hits.c:1907,1927`) does
        // NOT apply the same mask to bit scores — only a commented-out
        // `#if 0` assertion hints at it. Verify Rust matches that split.
        let mut kbp = KarlinBlk {
            lambda: 0.6,
            k: 0.4,
            log_k: 0.4_f64.ln(),
            h: 0.7,
            round_down: false,
        };
        // Without round_down, odd score 9 is used as-is.
        let bit_odd = kbp.raw_to_bit(9);
        let ev_odd = kbp.raw_to_evalue(9, 1.0e9);
        // Even score 8 is a reference point.
        let bit_even = kbp.raw_to_bit(8);
        let ev_even = kbp.raw_to_evalue(8, 1.0e9);
        assert_ne!(bit_odd, bit_even);
        // With round_down, odd score 9 rounds to 8 for E-value only.
        kbp.round_down = true;
        let bit_rd = kbp.raw_to_bit(9);
        let ev_rd = kbp.raw_to_evalue(9, 1.0e9);
        // Bit score IGNORES round_down (per NCBI):
        assert!((bit_rd - bit_odd).abs() < 1e-14);
        // E-value RESPECTS round_down (per NCBI):
        assert!((ev_rd - ev_even).abs() < ev_even.abs() * 1e-14);
        // round_down=true but even input → bit and evalue unchanged.
        let bit_even_rd = kbp.raw_to_bit(8);
        assert!((bit_even_rd - bit_even).abs() < 1e-14);
        let _ = ev_odd;
    }

    #[test]
    fn test_nucl_gapped_kbp_lookup_sets_round_down_for_2_minus_3() {
        // `(reward, penalty) = (2, -3)` scoring system has only
        // even-score entries, so NCBI sets `sbp->round_down = TRUE`.
        // Verify the returned KBP carries the flag.
        let ungapped = KarlinBlk {
            lambda: 0.55,
            k: 0.21,
            log_k: 0.21_f64.ln(),
            h: 0.46,
            round_down: false,
        };
        // Gap costs that exist in KBPT_2_3.
        let (kbp, rd) = nucl_gapped_kbp_lookup(4, 4, 2, -3, &ungapped).unwrap();
        assert!(rd, "expected round_down=true for (2,-3)");
        assert!(kbp.round_down, "returned KBP should carry round_down=true");
    }

    #[test]
    fn test_evalue_to_raw_pathological_inputs() {
        // Invalid lambda or searchsp returns `BLAST_SCORE_MIN`, matching
        // NCBI `BlastKarlinEtoS_simple` (`blast_stat.c:4054-4057`).
        let kbp = KarlinBlk {
            lambda: 0.0, // invalid: `denom <= 0.0 || lambda <= 0.0` triggers.
            k: 0.41,
            log_k: 0.41_f64.ln(),
            h: 0.78,
            round_down: false,
        };
        assert_eq!(kbp.evalue_to_raw(1.0, 1.0e9), BLAST_SCORE_MIN);
        let kbp_ok = KarlinBlk {
            lambda: 0.625,
            k: 0.41,
            log_k: 0.41_f64.ln(),
            h: 0.78,
            round_down: false,
        };
        // Zero search space → denom == 0 → return `BLAST_SCORE_MIN`.
        assert_eq!(kbp_ok.evalue_to_raw(1.0, 0.0), BLAST_SCORE_MIN);
        // Non-positive e-values are clamped to 1e-297 (NCBI
        // `BlastKarlinEtoS_simple`, `blast_stat.c:4059`) rather than
        // early-returning. The resulting raw score is finite and large.
        let clamped = kbp_ok.evalue_to_raw(-1.0, 1.0e9);
        let expected = kbp_ok.evalue_to_raw(1.0e-297, 1.0e9);
        assert_eq!(clamped, expected);
        assert_eq!(kbp_ok.evalue_to_raw(0.0, 1.0e9), expected);
    }

    #[test]
    fn test_karlin_blk_valid() {
        let kbp = KarlinBlk {
            lambda: 0.208,
            k: 0.049,
            log_k: 0.049_f64.ln(),
            h: 0.14,
            round_down: false,
        };
        assert!(kbp.is_valid());
    }

    #[test]
    fn test_raw_to_bit() {
        let kbp = KarlinBlk {
            lambda: 0.208,
            k: 0.049,
            log_k: 0.049_f64.ln(),
            h: 0.14,
            round_down: false,
        };
        let bit = kbp.raw_to_bit(50);
        assert!(bit > 0.0);
        // lambda * 50 / ln(2) - ln(K)/ln(2) ≈ 0.208*50/0.693 - ln(0.049)/0.693
        // ≈ 15.0 + 4.35 ≈ 19.35
        assert!((bit - 19.35).abs() < 0.5);
    }

    #[test]
    fn test_compute_ungapped_lambda() {
        // For reward=1, penalty=-3: known lambda ≈ 1.374
        let lambda = compute_ungapped_lambda(1, -3);
        assert!(
            (lambda - 1.374).abs() < 0.01,
            "lambda for 1/-3 should be ~1.374, got {}",
            lambda
        );

        // For reward=2, penalty=-3: different lambda
        let lambda2 = compute_ungapped_lambda(2, -3);
        assert!(
            lambda2 > 0.0 && lambda2 < lambda,
            "lambda for 2/-3 should be positive and < 1/-3 lambda"
        );
    }

    #[test]
    fn test_compute_ungapped_kbp() {
        let kbp = compute_ungapped_kbp(1, -3);
        assert!(kbp.lambda > 0.0);
        assert!(kbp.k > 0.0);
        assert!(kbp.h > 0.0);
        assert!(kbp.is_valid());
    }

    #[test]
    fn test_search_space() {
        let ss = compute_search_space(100, 1000000, 2000, 20);
        assert!(ss > 0.0);
        // (100-20) * (1000000 - 2000*20) = 80 * 960000 = 76800000
        assert_eq!(ss, 76800000.0);
    }

    // --- Ported from NCBI scoreblk_unit_test.cpp ---

    /// Port of GetScoreBlockNucl: verify nucleotide scoring properties.
    /// NCBI checks: alphabet_size=16, loscore=-3, hiscore=1, penalty=-3, reward=1
    #[test]
    fn test_scoreblk_nucl_properties() {
        let m = crate::matrix::nucleotide_matrix(1, -3);
        // Find min/max scores in the matrix
        let mut lo = i32::MAX;
        let mut hi = i32::MIN;
        for i in 0..16 {
            for j in 0..16 {
                if m[i][j] < lo {
                    lo = m[i][j];
                }
                if m[i][j] > hi {
                    hi = m[i][j];
                }
            }
        }
        assert_eq!(lo, -3, "nucleotide loscore should be -3");
        assert_eq!(hi, 1, "nucleotide hiscore should be 1");
    }

    /// Port of GetScoreBlockProtein: verify BLOSUM62 scoring properties.
    /// NCBI checks: alphabet_size=BLASTAA_SIZE(28), loscore=-4, hiscore=11
    #[test]
    fn test_scoreblk_protein_properties() {
        let m = &crate::matrix::BLOSUM62;
        let mut lo = i32::MAX;
        let mut hi = i32::MIN;
        for i in 0..crate::matrix::AA_SIZE {
            for j in 0..crate::matrix::AA_SIZE {
                if m[i][j] < lo {
                    lo = m[i][j];
                }
                if m[i][j] > hi {
                    hi = m[i][j];
                }
            }
        }
        assert_eq!(lo, -4, "BLOSUM62 loscore should be -4");
        assert_eq!(hi, 11, "BLOSUM62 hiscore should be 11");
    }

    /// Port of NuclGappedCalc: verify gapped KBP for reward=1, penalty=-2, gap_open=3, gap_extend=1.
    /// NCBI expects: Lambda≈1.32, K≈0.57, round_down=false
    #[test]
    fn test_nucl_gapped_calc_1_2_3_1() {
        let ungapped = compute_ungapped_kbp(1, -2);
        let (kbp, round_down) = nucl_gapped_kbp_lookup(3, 1, 1, -2, &ungapped).unwrap();
        assert!(!round_down);
        assert!(
            (kbp.lambda - 1.32).abs() < 0.01,
            "Lambda should be ~1.32, got {}",
            kbp.lambda
        );
        assert!(
            (kbp.k - 0.57).abs() < 0.01,
            "K should be ~0.57, got {}",
            kbp.k
        );
    }

    /// Port of NuclGappedCalc: verify alpha/beta for reward=1, penalty=-2, gap_open=3, gap_extend=1.
    /// NCBI expects: alpha≈1.3, beta≈-1.0
    #[test]
    fn test_nucl_alpha_beta_1_2_3_1() {
        let ungapped = compute_ungapped_kbp(1, -2);
        let (alpha, beta) = nucl_alpha_beta(1, -2, 3, 1, ungapped.lambda, ungapped.h, true);
        assert!(
            (alpha - 1.3).abs() < 0.01,
            "alpha should be ~1.3, got {}",
            alpha
        );
        assert!(
            (beta - (-1.0)).abs() < 0.01,
            "beta should be ~-1.0, got {}",
            beta
        );
    }

    /// Port of NuclGappedCalc: high gap costs fall back to ungapped params.
    /// reward=1, penalty=-2, gap_open=4, gap_extend=2 → copies ungapped Lambda/K.
    #[test]
    fn test_nucl_gapped_fallback_to_ungapped() {
        let ungapped = compute_ungapped_kbp(1, -2);
        let (kbp, round_down) = nucl_gapped_kbp_lookup(4, 2, 1, -2, &ungapped).unwrap();
        assert!(!round_down);
        assert_eq!(
            kbp.lambda, ungapped.lambda,
            "High gap costs should fall back to ungapped Lambda"
        );
        assert_eq!(
            kbp.k, ungapped.k,
            "High gap costs should fall back to ungapped K"
        );
    }

    /// Port of NuclGappedCalc: ungapped alpha/beta when gap costs exceed table maximum.
    /// Alpha = Lambda/H, beta = 0.0
    #[test]
    fn test_nucl_alpha_beta_ungapped_fallback() {
        let ungapped = compute_ungapped_kbp(1, -2);
        let (alpha, beta) = nucl_alpha_beta(1, -2, 4, 2, ungapped.lambda, ungapped.h, true);
        assert!(
            (alpha - ungapped.lambda / ungapped.h).abs() < 1e-6,
            "alpha should equal Lambda/H for unsupported gap costs"
        );
        assert_eq!(beta, 0.0);
    }

    /// Port of NuclGappedCalc: scaled-up values (reward=10, penalty=-20, gap_open=30, gap_extend=10).
    /// GCD=10, so internally maps to (1, -2, 3, 1). Lambda should be ~0.132 (1.32/10).
    #[test]
    fn test_nucl_gapped_scaled_values() {
        let ungapped = compute_ungapped_kbp(10, -20);
        let (kbp, round_down) = nucl_gapped_kbp_lookup(30, 10, 10, -20, &ungapped).unwrap();
        assert!(!round_down);
        assert!(
            (kbp.lambda - 0.132).abs() < 0.01,
            "Lambda for 10/-20/30/10 should be ~0.132 (scaled), got {}",
            kbp.lambda
        );
        assert!(
            (kbp.k - 0.57).abs() < 0.01,
            "K for 10/-20/30/10 should be ~0.57 (unscaled), got {}",
            kbp.k
        );
    }

    /// Port of NuclGappedCalc: reward=2, penalty=-7, gap_open=4, gap_extend=2 should set round_down=true.
    /// NCBI expects: Lambda≈0.675, K≈0.62, round_down=true
    #[test]
    fn test_nucl_gapped_round_down() {
        let ungapped = compute_ungapped_kbp(2, -7);
        let (kbp, round_down) = nucl_gapped_kbp_lookup(4, 2, 2, -7, &ungapped).unwrap();
        assert!(round_down, "2/-7 should set round_down=true");
        assert!(
            (kbp.lambda - 0.675).abs() < 0.01,
            "Lambda should be ~0.675, got {}",
            kbp.lambda
        );
        assert!(
            (kbp.k - 0.62).abs() < 0.01,
            "K should be ~0.62, got {}",
            kbp.k
        );
    }

    /// Port of NuclGappedCalc: invalid gap costs should return Err.
    /// reward=4, penalty=-5, gap_open=3, gap_extend=2 is not in the table.
    #[test]
    fn test_nucl_gapped_invalid_gap_costs() {
        let ungapped = compute_ungapped_kbp(4, -5);
        let result = nucl_gapped_kbp_lookup(3, 2, 4, -5, &ungapped);
        assert!(
            result.is_err(),
            "gap_open=3, gap_extend=2 for 4/-5 should fail"
        );
    }

    /// Port of NuclGappedCalc: invalid gap costs should return Err.
    /// reward=1, penalty=-2, gap_open=1, gap_extend=3 is not in the table.
    #[test]
    fn test_nucl_gapped_invalid_gap_costs_2() {
        let ungapped = compute_ungapped_kbp(1, -2);
        let result = nucl_gapped_kbp_lookup(1, 3, 1, -2, &ungapped);
        assert!(
            result.is_err(),
            "gap_open=1, gap_extend=3 for 1/-2 should fail"
        );
    }

    /// Port of NuclGappedCalc: unsupported reward/penalty pair returns Err.
    /// reward=2, penalty=-1 is not supported.
    #[test]
    fn test_nucl_gapped_unsupported_scores() {
        let ungapped = compute_ungapped_kbp(2, -1);
        let result = nucl_gapped_kbp_lookup(1, 3, 2, -1, &ungapped);
        assert!(result.is_err(), "2/-1 is not a supported scoring pair");
    }

    /// Port of EqualRewardPenaltyLHtoK: when reward==|penalty|, K should be 1/3.
    /// NCBI uses full Blast_ScoreBlkKbpIdealCalc with DP-based K computation.
    /// The full ungapped_kbp_calc function uses the same algorithm.
    #[test]
    fn test_equal_reward_penalty_k() {
        // Use the full KBP computation via ungapped_kbp_calc
        // Query: uniform ACGT, 1000 bases
        let query: Vec<u8> = (0..1000).map(|i| (i % 4) as u8).collect();
        let contexts = vec![UngappedKbpContext {
            query_offset: 0,
            query_length: 1000,
            is_valid: true,
        }];
        let m = crate::matrix::nucleotide_matrix(2, -2);
        let results = ungapped_kbp_calc(
            &query,
            &contexts,
            -2,
            2,
            16,
            &(4..16).collect::<Vec<u8>>(),
            &|i, j| m[i][j],
        );
        let kbp = results[0].as_ref().unwrap();
        assert!(
            (kbp.k - 1.0 / 3.0).abs() < 0.02,
            "K for equal reward/penalty should be ~1/3, got {}",
            kbp.k
        );
    }

    /// Port of BlastResFreqStdCompNucleotideTest: nucleotide base frequencies should be 0.25 each.
    #[test]
    fn test_nucleotide_standard_frequencies() {
        let freqs = &crate::matrix::NT_FREQUENCIES;
        // First 4 bases (A, C, G, T) should each be 0.25
        for i in 0..4 {
            assert!(
                (freqs[i] - 0.25).abs() < 0.001,
                "Base {} frequency should be 0.25, got {}",
                i,
                freqs[i]
            );
        }
        // Ambiguity codes should be 0
        for i in 4..freqs.len() {
            assert!(
                freqs[i].abs() < 0.001,
                "Ambiguity code {} frequency should be 0, got {}",
                i,
                freqs[i]
            );
        }
    }

    /// Port of BlastResFreqStdCompProteinTest: verify specific amino acid frequencies.
    /// Compact AA_FREQUENCIES indices (alphabetical): A=0,C=1,D=2,E=3,F=4,G=5,
    /// H=6,I=7,K=8,L=9,M=10,N=11,P=12,Q=13,R=14,S=15,T=16,V=17,W=18,Y=19.
    /// Expected x100000 from NCBI's `Robinson_prob` (`blast_stat.c:1795`)
    /// divided by 1000: F→3856, M→2243, Y→3216, C→1925.
    #[test]
    fn test_protein_standard_frequencies() {
        let freqs = &crate::matrix::AA_FREQUENCIES;
        let check = |idx: usize, expected_x100k: i32, name: &str| {
            let got = (freqs[idx] * 100000.0).round() as i32;
            assert!(
                (got - expected_x100k).abs() <= 5,
                "{} (idx {}) expected ~{}, got {} (raw: {})",
                name,
                idx,
                expected_x100k,
                got,
                freqs[idx]
            );
        };
        check(1, 1925, "C"); // Cys
        check(4, 3856, "F"); // Phe
        check(10, 2243, "M"); // Met
        check(19, 3216, "Y"); // Tyr
    }

    /// Port of EvalueForProteinFSC: E-values should never be negative.
    /// Uses BLOSUM62 with gap_open=11, gap_extend=1.
    #[test]
    fn test_evalue_never_negative() {
        let kbp = lookup_protein_params(11, 1).unwrap();
        let kbp = KarlinBlk {
            lambda: kbp.lambda,
            k: kbp.k,
            log_k: kbp.k.ln(),
            h: kbp.h,
            round_down: false,
        };
        // Test cases from NCBI stat_unit_test.cpp
        let cases = [
            (1201, 294, 422),
            (1204, 294, 416),
            (1179, 294, 418),
            (2332, 1801, 1671),
        ];
        for &(score, len1, len2) in &cases {
            let ss = len1 as f64 * len2 as f64;
            let evalue = kbp.raw_to_evalue(score, ss);
            assert!(
                evalue >= 0.0,
                "E-value should be >= 0 for score={}, lens=({},{}), got {}",
                score,
                len1,
                len2,
                evalue
            );
        }
    }

    /// Verify gapped KBP lookup for all supported nucleotide reward/penalty combos.
    /// Not all combos have a linear (0,0) entry; test the first affine entry from each table.
    #[test]
    fn test_all_nucl_reward_penalty_combos() {
        // (reward, penalty, gap_open, gap_extend) - first affine entry from each KBPT table
        let cases: &[(i32, i32, i32, i32)] = &[
            (1, -5, 3, 3),
            (1, -4, 1, 2),
            (1, -3, 2, 2),
            (1, -2, 2, 2),
            (1, -1, 3, 2),
            (2, -7, 2, 4),
            (2, -5, 2, 4),
            (2, -3, 4, 4),
            (3, -4, 6, 3),
            (3, -2, 5, 5),
            (4, -5, 6, 5),
            (5, -4, 10, 6),
        ];
        for &(reward, penalty, go, ge) in cases {
            let ungapped = compute_ungapped_kbp(reward, penalty);
            let result = nucl_gapped_kbp_lookup(go, ge, reward, penalty, &ungapped);
            assert!(
                result.is_ok(),
                "Gapped lookup should work for {}/{}/{}/{}",
                reward,
                penalty,
                go,
                ge
            );
            let (kbp, _) = result.unwrap();
            assert!(
                kbp.is_valid(),
                "KBP should be valid for {}/{}/{}/{}",
                reward,
                penalty,
                go,
                ge
            );
        }
    }

    /// Verify gapped KBP for all entries in KBPT_1_3 table via nucl_gapped_kbp_lookup.
    #[test]
    fn test_nucl_gapped_all_1_3_entries() {
        let ungapped = compute_ungapped_kbp(1, -3);
        // Entries from KBPT_1_3 (the C-compatible table)
        let entries: &[(i32, i32, f64, f64)] = &[
            (2, 2, 1.37, 0.70),
            (1, 2, 1.35, 0.64),
            (0, 2, 1.25, 0.42),
            (2, 1, 1.34, 0.60),
            (1, 1, 1.21, 0.34),
        ];
        for &(go, ge, expected_lambda, expected_k) in entries {
            let result = nucl_gapped_kbp_lookup(go, ge, 1, -3, &ungapped);
            assert!(
                result.is_ok(),
                "1/-3 gap_open={} gap_extend={} should work",
                go,
                ge
            );
            let (kbp, _) = result.unwrap();
            assert!(
                (kbp.lambda - expected_lambda).abs() < 0.01,
                "Lambda mismatch for 1/-3/{}/{}: expected {}, got {}",
                go,
                ge,
                expected_lambda,
                kbp.lambda
            );
            assert!(
                (kbp.k - expected_k).abs() < 0.01,
                "K mismatch for 1/-3/{}/{}: expected {}, got {}",
                go,
                ge,
                expected_k,
                kbp.k
            );
        }
    }

    /// Verify protein gapped KBP for all BLOSUM62 table entries.
    #[test]
    fn test_protein_gapped_all_blosum62_entries() {
        for &(go, ge, expected_lambda, expected_k, expected_h, _alpha, _beta) in BLOSUM62_PARAMS {
            let params = lookup_protein_params(go, ge);
            assert!(
                params.is_some(),
                "BLOSUM62 gap_open={} gap_extend={} should be in table",
                go,
                ge
            );
            let p = params.unwrap();
            assert!((p.lambda - expected_lambda).abs() < 1e-6);
            assert!((p.k - expected_k).abs() < 1e-6);
            assert!((p.h - expected_h).abs() < 1e-6);
        }
    }

    /// Verify length adjustment converges for typical search parameters.
    #[test]
    fn test_length_adjustment_exact_convergence() {
        let kbp = lookup_protein_params(11, 1).unwrap();
        let (adj, converged) = compute_length_adjustment_exact(
            kbp.k,
            kbp.k.ln(),
            kbp.alpha / kbp.lambda,
            kbp.beta,
            300,
            1_000_000,
            5000,
        );
        assert!(converged, "Length adjustment should converge");
        assert!(adj > 0, "Length adjustment should be positive");
        assert!(adj < 300, "Length adjustment should be < query length");
    }

    /// Verify E-value <-> raw score round-trip consistency.
    #[test]
    fn test_evalue_raw_score_roundtrip() {
        let kbp = KarlinBlk {
            lambda: 0.267,
            k: 0.041,
            log_k: 0.041_f64.ln(),
            h: 0.14,
            round_down: false,
        };
        let search_space = 1e9;
        // For a given raw score, convert to evalue, then back to raw score
        for raw in [20, 50, 100, 200] {
            let evalue = kbp.raw_to_evalue(raw, search_space);
            let recovered = kbp.evalue_to_raw(evalue, search_space);
            assert!(
                (recovered - raw).abs() <= 1,
                "Round-trip failed: raw={} -> evalue={:.2e} -> recovered={}",
                raw,
                evalue,
                recovered
            );
        }
    }

    /// Verify ungapped KBP computation for protein (BLOSUM62).
    #[test]
    fn test_protein_ungapped_kbp_values() {
        let kbp = protein_ungapped_kbp();
        assert!(kbp.is_valid());
        assert!((kbp.lambda - 0.3176).abs() < 0.01);
        assert!((kbp.k - 0.134).abs() < 0.01);
        assert!((kbp.h - 0.401).abs() < 0.01);
    }

    /// Verify ungapped KBP calc produces valid results for context-based computation.
    #[test]
    fn test_ungapped_kbp_calc_nucleotide() {
        // A simple BLASTNA-encoded query: ACGTACGT = [0,1,2,3,0,1,2,3]
        let query: Vec<u8> = vec![0, 1, 2, 3, 0, 1, 2, 3];
        let contexts = vec![UngappedKbpContext {
            query_offset: 0,
            query_length: 8,
            is_valid: true,
        }];
        let m = crate::matrix::nucleotide_matrix(1, -3);
        let results = ungapped_kbp_calc(
            &query,
            &contexts,
            -3,
            1,
            16,
            &(4..16).collect::<Vec<u8>>(), // ambiguous codes 4-15
            &|i, j| m[i][j],
        );
        assert_eq!(results.len(), 1);
        assert!(results[0].is_some(), "Should produce valid KBP");
        let kbp = results[0].as_ref().unwrap();
        assert!(kbp.is_valid());
        assert!(kbp.lambda > 0.0);
        assert!(kbp.k > 0.0);
        assert!(kbp.h > 0.0);
    }

    /// Verify invalid context (zero-length) produces None.
    #[test]
    fn test_ungapped_kbp_calc_invalid_context() {
        let query: Vec<u8> = vec![0, 1, 2, 3];
        let contexts = vec![UngappedKbpContext {
            query_offset: 0,
            query_length: 0,
            is_valid: false,
        }];
        let m = crate::matrix::nucleotide_matrix(1, -3);
        let results = ungapped_kbp_calc(
            &query,
            &contexts,
            -3,
            1,
            16,
            &(4..16).collect::<Vec<u8>>(),
            &|i, j| m[i][j],
        );
        assert!(results[0].is_none());
    }

    /// Verify search space with zero length adjustment.
    #[test]
    fn test_search_space_no_adjustment() {
        let ss = compute_search_space(500, 10_000_000, 1000, 0);
        assert_eq!(ss, 500.0 * 10_000_000.0);
    }

    /// Verify search space clamps to 1 when adjustment exceeds lengths.
    #[test]
    fn test_search_space_clamped() {
        let ss = compute_search_space(10, 100, 5, 50);
        // eff_query = max(10-50, 1) = 1
        // eff_db = max(100 - 5*50, 1) = max(-150, 1) = 1
        assert_eq!(ss, 1.0);
    }

    /// Verify length adjustment is 0 for degenerate inputs.
    #[test]
    fn test_length_adjustment_degenerate() {
        let kbp = KarlinBlk {
            lambda: 0.0,
            k: 0.0,
            log_k: f64::NEG_INFINITY,
            h: 0.0,
            round_down: false,
        };
        let adj = compute_length_adjustment(100, 1000000, 100, &kbp);
        assert_eq!(adj, 0);
    }
}
