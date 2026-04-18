//! Rust equivalent of ncbi_math.c — math utility functions for BLAST.

// NCBI fdlibm/gamma coefficients in this module are declared with >16
// significant digits. Keep them byte-identical to the C source; f64
// rounds the same as C so the extra digits are harmless.
#![allow(clippy::excessive_precision)]

/// Port of NCBI `NCBIMATH_LN2` (`ncbi_math.h:161`). Bit-identical to
/// `std::f64::consts::LN_2`; re-exported here so call sites mirror
/// NCBI naming.
pub const NCBIMATH_LN2: f64 = std::f64::consts::LN_2;
/// Port of NCBI `NCBIMATH_LNPI` (`ncbi_math.h:163`). Bit-identical to
/// `std::f64::consts::PI.ln()`.
pub const NCBIMATH_LNPI: f64 = 1.144_729_885_849_400_2;
/// `NCBIMATH_PI` = pi (`ncbi_math.h`). Bit-identical to `std::f64::consts::PI`.
pub const NCBIMATH_PI: f64 = std::f64::consts::PI;

/// exp(x) - 1 for small x (avoids catastrophic cancellation).
pub fn expm1(x: f64) -> f64 {
    let absx = x.abs();
    if absx > 0.33 {
        return x.exp() - 1.0;
    }
    if absx < 1e-16 {
        return x;
    }
    // Taylor series: x + x^2/2! + x^3/3! + ...
    x * (1.0
        + x * (1.0 / 2.0
            + x * (1.0 / 6.0
                + x * (1.0 / 24.0
                    + x * (1.0 / 120.0
                        + x * (1.0 / 720.0
                            + x * (1.0 / 5040.0
                                + x * (1.0 / 40320.0
                                    + x * (1.0 / 362880.0
                                        + x * (1.0 / 3628800.0
                                            + x * (1.0 / 39916800.0
                                                + x * (1.0 / 479001600.0
                                                    + x / 6227020800.0))))))))))))
}

/// ln(1+x) for small x (avoids catastrophic cancellation).
/// Verbatim port of NCBI `BLAST_Log1p` (`ncbi_math.c:64`): 500-term
/// cap, DBL_EPSILON convergence check inside a two-step loop that
/// processes one odd and one even Taylor term per iteration.
pub fn log1p(x: f64) -> f64 {
    if x.abs() >= 0.2 {
        return (x + 1.0).ln();
    }
    let mut sum = 0.0;
    let mut y = x;
    let mut i = 0i32;
    while i < 500 {
        i += 1;
        sum += y / i as f64;
        if y.abs() < f64::EPSILON {
            break;
        }
        y *= x;
        i += 1;
        sum -= y / i as f64;
        if y < f64::EPSILON {
            break;
        }
        y *= x;
    }
    sum
}

/// Port of NCBI `BLAST_Nint` (`ncbi_math.c:437`): round half-away-from-zero
/// and cast to integer. Equivalent to Rust's `f64::round() as i64`, but
/// named for parity-visibility at call sites that mirror NCBI code.
pub fn nint(x: f64) -> i64 {
    (x + if x >= 0.0 { 0.5 } else { -0.5 }) as i64
}

/// Port of NCBI `BLAST_Gcd` (`ncbi_math.c:405`). Only `b` is absolute-valued;
/// the swap ensures `a >= b` on entry to the Euclidean loop. In BLAST
/// practice both args are non-negative, so the weirdness where a negative
/// `a` could survive the swap is unreachable. Port is character-for-character.
pub fn gcd(a: i32, b: i32) -> i32 {
    let mut a = a;
    let mut b = b.abs();
    if b > a {
        std::mem::swap(&mut a, &mut b);
    }
    while b != 0 {
        let c = a % b;
        a = b;
        b = c;
    }
    a
}

/// Port of NCBI `BLAST_Powi` (`ncbi_math.c:444`): integer-power via
/// repeated squaring. Matches NCBI's special cases: `x^0 = 1` (even
/// when `x == 0`), `0^n = 0` for `n > 0`, `0^n = HUGE_VAL` for `n < 0`.
/// Semantically identical to Rust's `f64::powi`; named for parity
/// at sites that mirror NCBI.
pub fn powi(mut x: f64, mut n: i32) -> f64 {
    if n == 0 {
        return 1.0;
    }
    if x == 0.0 {
        if n < 0 {
            return f64::INFINITY;
        }
        return 0.0;
    }
    if n < 0 {
        x = 1.0 / x;
        n = -n;
    }
    let mut y = 1.0;
    while n > 0 {
        if n & 1 != 0 {
            y *= x;
        }
        n /= 2;
        x *= x;
    }
    y
}

/// Coefficients used by NCBI's Lanczos-style ln(gamma) approximation.
/// Verbatim from `ncbi_math.c:140-152` (`_default_gamma_coef`).
const GAMMA_COEF: [f64; 11] = [
    4.694_580_336_184_385e4,
    -1.560_605_207_784_446e5,
    2.065_049_568_014_106e5,
    -1.388_934_775_095_388e5,
    5.031_796_415_085_709e4,
    -9.601_592_329_182_778e3,
    8.785_855_930_895_250e2,
    -3.155_153_906_098_611e1,
    2.908_143_421_162_229e-1,
    -2.319_827_630_494_973e-4,
    1.251_639_670_050_933e-10,
];

/// Port of NCBI `s_GeneralLnGamma(x, 0)` (`ncbi_math.c:162`): computes
/// `ln(Gamma(x))` for `x >= 1` via a Lanczos-style series using
/// `GAMMA_COEF` and the identity `ln(Gamma(x)) = ln(y) + (ln(pi)+ln(2))/2
/// + (xx+0.5)*ln(tx+0.5) - (tx+0.5)`, where `xx = x-1`, `tx = xx+11`,
/// and `y = 1 + sum_i coef[i] / (tx - i)`. 10-digit accuracy per NCBI.
fn s_ln_gamma(x: f64) -> f64 {
    let xgamma_dim = GAMMA_COEF.len();
    let xx = x - 1.0;
    let tx = xx + xgamma_dim as f64;
    // Walk `coef` pointer backward as in the C code; equivalent to
    // summing coef[n-1]/tx + coef[n-2]/(tx-1) + ... + coef[0]/(tx-n+1).
    let mut tmp = tx;
    let mut value = GAMMA_COEF[xgamma_dim - 1] / tmp;
    for i in (0..(xgamma_dim - 1)).rev() {
        tmp -= 1.0;
        value += GAMMA_COEF[i] / tmp;
    }
    value += 1.0;
    // `s_LogDerivative(0, y)` == `log(y[0])`. NCBI `ncbi_math.c:197`:
    //   value += ((NCBIMATH_LNPI + NCBIMATH_LN2) / 2.) + (xx+0.5)*log(tmp) - tmp;
    let ln_y = value.ln();
    let tmp_half = tx + 0.5;
    ln_y + (NCBIMATH_LNPI + NCBIMATH_LN2) * 0.5 + (xx + 0.5) * tmp_half.ln() - tmp_half
}

/// Compute the natural log of the factorial of n.
/// For small `n` (< 35) uses an exact precomputed table matching
/// NCBI's `kPrecomputedFactorial` (35 entries, `ncbi_math.c:297`); for
/// larger `n` delegates to `s_ln_gamma` (NCBI `s_LnGamma`). Port of
/// NCBI `BLAST_LnGammaInt` (`ncbi_math.c:323`).
pub fn ln_factorial(n: i32) -> f64 {
    // Exact values of ln(n!) for n = 0..34. The table length matches
    // NCBI's `kPrecomputedFactorial` so the `s_ln_gamma` branch is only
    // reached for the same large-n inputs where NCBI falls back to
    // `s_LnGamma`.
    // ln(2!) = ln(2) happens to equal `std::f64::consts::LN_2`; the
    // literal is verbatim from NCBI `kPrecomputedFactorial` so we keep
    // it rather than rewriting to reference the stdlib constant.
    #[allow(clippy::approx_constant)]
    const LN_FACT_TABLE: [f64; 35] = [
        0.0,                             // ln(0!)
        0.0,                             // ln(1!)
        0.693_147_180_559_944_95,        // ln(2!) = LN_2
        1.791_759_469_228_055_4,         // ln(3!)
        3.178_053_830_347_944_86,        // ln(4!)
        4.787_491_742_782_046_7,         // ln(5!)
        6.579_251_212_010_102_1,         // ln(6!)
        8.525_161_361_065_414_67,        // ln(7!)
        10.604_602_902_745_249_1,        // ln(8!)
        12.801_827_480_081_467_3,        // ln(9!)
        15.104_412_573_075_514_1,        // ln(10!)
        17.502_307_845_873_886_5,        // ln(11!)
        19.987_214_495_661_884_7,        // ln(12!)
        22.552_163_853_123_421,          // ln(13!)
        25.191_221_182_738_683_4,        // ln(14!)
        27.899_271_383_840_890_3,        // ln(15!)
        30.671_860_106_080_671_9,        // ln(16!)
        33.505_073_450_136_890_8,        // ln(17!)
        36.395_445_208_033_052_6,        // ln(18!)
        39.339_884_187_199_494_6,        // ln(19!)
        42.335_616_460_753_485_1,        // ln(20!)
        45.380_138_898_476_907_6,        // ln(21!)
        48.471_181_351_835_220_1,        // ln(22!)
        51.606_675_567_764_376_9,        // ln(23!)
        54.784_729_398_112_318_7,        // ln(24!)
        58.003_605_222_980_517_9,        // ln(25!)
        61.261_701_761_002_008_5,        // ln(26!)
        64.557_538_627_006_337_6,        // ln(27!)
        67.889_743_137_181_525_9,        // ln(28!)
        71.257_038_967_168_000_5,        // ln(29!)
        74.658_236_348_830_172_3,        // ln(30!)
        78.092_223_553_315_307_1,        // ln(31!)
        81.557_959_456_115_028_7,        // ln(32!)
        85.054_467_017_581_515_6,        // ln(33!)
        88.580_827_542_197_681_6,        // ln(34!)
    ];
    if n < 0 {
        return 0.0;
    }
    let idx = n as usize;
    if idx < LN_FACT_TABLE.len() {
        return LN_FACT_TABLE[idx];
    }
    // `ln(n!) = ln(Gamma(n+1))`, matching NCBI `BLAST_LnFactorial` which
    // forwards to `s_LnGamma(x + 1.0)` (`ncbi_math.c:478`).
    s_ln_gamma(n as f64 + 1.0)
}

/// Port of NCBI `BLAST_LnGammaInt` (`ncbi_math.c:323`): returns
/// `ln(Gamma(n))` for positive integer `n`. Equivalent to
/// `ln_factorial(n - 1)`.
pub fn ln_gamma_int(n: i32) -> f64 {
    if n <= 0 {
        return 0.0;
    }
    // Match NCBI's fast path: small `n` uses the factorial table via
    // `log(kPrecomputedFactorial[n-1])`.
    ln_factorial(n - 1)
}

/// Port of NCBI `BLAST_RombergIntegrate` (`ncbi_math.c:346`).
/// Adaptive Romberg quadrature of `f` over `[p, q]`, stopping when
/// successive Romberg extrapolations agree within `eps`. Returns
/// `f64::INFINITY` when the integrand returns `HUGE_VAL` or the method
/// fails to converge within `MAX_DIAGS` iterations.
///
/// `epsit` — minimum number of consecutive iterations that must satisfy
/// `eps` before returning (clamped to `[1, 3]`).
/// `itmin` — minimum iterations (clamped to `[1, MAX_DIAGS-1]`).
pub fn romberg_integrate<F>(f: F, p: f64, q: f64, eps: f64, epsit: i32, itmin: i32) -> f64
where
    F: Fn(f64) -> f64,
{
    const MAX_DIAGS: usize = 20;
    let huge = f64::INFINITY;
    let mut romb = [0.0f64; MAX_DIAGS];

    // NCBI uses `itmin = MAX(1, itmin); itmin = MIN(itmin, MAX_DIAGS-1);`
    // (`ncbi_math.c:359-364`). `clamp` is semantically identical here.
    let itmin = itmin.clamp(1, (MAX_DIAGS - 1) as i32) as usize;
    let epsit = epsit.clamp(1, 3) as usize;
    let epsck = itmin as isize - epsit as isize;

    let mut npts: i64 = 1;
    let mut h = q - p;
    let x0 = f(p);
    if x0.abs() == huge {
        return x0;
    }
    let y0 = f(q);
    if y0.abs() == huge {
        return y0;
    }
    romb[0] = 0.5 * h * (x0 + y0);
    let mut epsit_cnt: usize = 0;
    let mut i = 1;
    while i < MAX_DIAGS {
        let mut sum = 0.0;
        let mut x = p + 0.5 * h;
        for _ in 0..npts {
            let y = f(x);
            if y.abs() == huge {
                return y;
            }
            sum += y;
            x += h;
        }
        romb[i] = 0.5 * (romb[i - 1] + h * sum);

        // Update Romberg array with new column.
        let mut n: i64 = 4;
        let mut j = (i as isize) - 1;
        while j >= 0 {
            romb[j as usize] =
                (n as f64 * romb[(j + 1) as usize] - romb[j as usize]) / (n - 1) as f64;
            n *= 4;
            j -= 1;
        }

        if (i as isize) > epsck {
            if (romb[1] - romb[0]).abs() > eps * romb[0].abs() {
                epsit_cnt = 0;
            } else {
                epsit_cnt += 1;
                if i >= itmin && epsit_cnt >= epsit {
                    return romb[0];
                }
            }
        }

        npts *= 2;
        h *= 0.5;
        i += 1;
    }
    huge
}

/// Port of NCBI `NCBI_ErfC` (`ncbi_erf.c:407`): Sun-style fdlibm erfc
/// with piecewise rational approximations and IEEE 754 bit manipulation.
/// Aliased by NCBI to `BLAST_ErfC` (`ncbi_erf.c:492`). Accurate to full
/// `f64` precision across all finite inputs.
pub fn erfc(x: f64) -> f64 {
    const TINY: f64 = 1e-300;
    const HALF: f64 = 0.5;
    const ONE: f64 = 1.0;
    const TWO: f64 = 2.0;
    // erx = (float)0.84506291151
    const ERX: f64 = 8.450_629_115_104_675e-1;

    // Coefficients for approximation to erf on [0, 0.84375]
    const PP0: f64 = 1.283_791_670_955_125_6e-1;
    const PP1: f64 = -3.250_421_072_470_015e-1;
    const PP2: f64 = -2.848_174_957_559_851e-2;
    const PP3: f64 = -5.770_270_296_489_442e-3;
    const PP4: f64 = -2.376_301_665_665_016_3e-5;
    const QQ1: f64 = 3.979_172_239_591_553_6e-1;
    const QQ2: f64 = 6.502_224_998_876_729e-2;
    const QQ3: f64 = 5.081_306_281_875_765_5e-3;
    const QQ4: f64 = 1.324_947_380_043_216_4e-4;
    const QQ5: f64 = -3.960_228_278_775_368e-6;

    // Coefficients for approximation to erf in [0.84375, 1.25]
    const PA0: f64 = -2.362_118_560_752_659_4e-3;
    const PA1: f64 = 4.148_561_186_837_483_3e-1;
    const PA2: f64 = -3.722_078_760_357_013e-1;
    const PA3: f64 = 3.183_466_199_011_617_5e-1;
    const PA4: f64 = -1.108_946_942_823_966_8e-1;
    const PA5: f64 = 3.547_830_432_561_823_5e-2;
    const PA6: f64 = -2.166_375_594_868_790_8e-3;
    const QA1: f64 = 1.064_208_804_008_442_2e-1;
    const QA2: f64 = 5.403_979_177_021_710_5e-1;
    const QA3: f64 = 7.182_865_441_419_627e-2;
    const QA4: f64 = 1.261_712_198_087_616_4e-1;
    const QA5: f64 = 1.363_708_391_202_905_1e-2;
    const QA6: f64 = 1.198_449_984_679_910_7e-2;

    // Coefficients for erfc in [1.25, 1/0.35]
    const RA0: f64 = -9.864_944_034_847_148e-3;
    const RA1: f64 = -6.938_585_727_071_817e-1;
    const RA2: f64 = -1.055_862_622_532_329_1e1;
    const RA3: f64 = -6.237_533_245_032_601e1;
    const RA4: f64 = -1.623_966_694_625_734_7e2;
    const RA5: f64 = -1.846_050_929_067_110_4e2;
    const RA6: f64 = -8.128_743_550_630_659e1;
    const RA7: f64 = -9.814_329_344_169_145e0;
    const SA1: f64 = 1.965_127_166_743_925_7e1;
    const SA2: f64 = 1.376_577_541_435_190_4e2;
    const SA3: f64 = 4.345_658_774_752_292_3e2;
    const SA4: f64 = 6.453_872_717_332_679e2;
    const SA5: f64 = 4.290_081_400_275_678_3e2;
    const SA6: f64 = 1.086_350_055_417_794_4e2;
    const SA7: f64 = 6.570_249_770_319_282e0;
    const SA8: f64 = -6.042_441_521_485_81e-2;

    // Coefficients for erfc in [1/0.35, 28]
    const RB0: f64 = -9.864_942_924_700_099e-3;
    const RB1: f64 = -7.992_832_376_805_23e-1;
    const RB2: f64 = -1.775_795_491_775_475_1e1;
    const RB3: f64 = -1.606_363_848_558_219e2;
    const RB4: f64 = -6.375_664_433_683_896e2;
    const RB5: f64 = -1.025_095_131_611_077_3e3;
    const RB6: f64 = -4.835_191_916_086_514e2;
    const SB1: f64 = 3.033_806_074_348_246e1;
    const SB2: f64 = 3.257_925_129_965_739e2;
    const SB3: f64 = 1.536_729_586_084_437e3;
    const SB4: f64 = 3.199_858_219_508_595_5e3;
    const SB5: f64 = 2.553_050_406_433_164_4e3;
    const SB6: f64 = 4.745_285_412_069_554e2;
    const SB7: f64 = -2.244_095_244_658_581_8e1;

    // NCBI uses __HI(x) to read the upper 32 bits of the IEEE754 double.
    // Rust equivalent: `(x.to_bits() >> 32) as u32`, then reinterpret as
    // signed for the `hx < 0` / `hx < 0x3fd00000` comparisons (NCBI's
    // `__HI` returns `int`).
    let bits = x.to_bits();
    let hx = (bits >> 32) as u32 as i32;
    let ix = hx & 0x7fff_ffff;
    if ix >= 0x7ff0_0000 {
        // erfc(nan) = nan; erfc(+inf) = 0, erfc(-inf) = 2
        return ((hx as u32 >> 31) << 1) as f64 + ONE / x;
    }

    if ix < 0x3feb_0000 {
        // |x| < 0.84375
        if ix < 0x3c70_0000 {
            // |x| < 2^-56
            return ONE - x;
        }
        let z = x * x;
        let r = PP0 + z * (PP1 + z * (PP2 + z * (PP3 + z * PP4)));
        let s = ONE + z * (QQ1 + z * (QQ2 + z * (QQ3 + z * (QQ4 + z * QQ5))));
        let y = r / s;
        if hx < 0x3fd0_0000 {
            // x < 1/4 (signed compare pulls in x < 0 via sign bit)
            return ONE - (x + x * y);
        }
        let mut r = x * y;
        r += x - HALF;
        return HALF - r;
    }

    if ix < 0x3ff4_0000 {
        // 0.84375 <= |x| < 1.25
        let s = x.abs() - ONE;
        let p = PA0 + s * (PA1 + s * (PA2 + s * (PA3 + s * (PA4 + s * (PA5 + s * PA6)))));
        let q = ONE + s * (QA1 + s * (QA2 + s * (QA3 + s * (QA4 + s * (QA5 + s * QA6)))));
        if hx >= 0 {
            return (ONE - ERX) - p / q;
        }
        return ONE + (ERX + p / q);
    }

    if ix < 0x403c_0000 {
        // |x| < 28
        let ax = x.abs();
        let s = ONE / (ax * ax);
        let (r, ss) = if ix < 0x4006_db6d {
            // |x| < 1/0.35
            let r = RA0
                + s * (RA1
                    + s * (RA2 + s * (RA3 + s * (RA4 + s * (RA5 + s * (RA6 + s * RA7))))));
            let ss = ONE
                + s * (SA1
                    + s * (SA2
                        + s * (SA3 + s * (SA4 + s * (SA5 + s * (SA6 + s * (SA7 + s * SA8)))))));
            (r, ss)
        } else {
            // |x| >= 1/0.35
            if hx < 0 && ix >= 0x4018_0000 {
                // x < -6
                return TWO - TINY;
            }
            let r = RB0
                + s * (RB1 + s * (RB2 + s * (RB3 + s * (RB4 + s * (RB5 + s * RB6)))));
            let ss = ONE
                + s * (SB1
                    + s * (SB2 + s * (SB3 + s * (SB4 + s * (SB5 + s * (SB6 + s * SB7))))));
            (r, ss)
        };
        // z has low 32 bits zeroed (single-precision version of |x|)
        // NCBI: `z = x; __LO(z) = 0;` after `x = fabs(x)`.
        let z = f64::from_bits(ax.to_bits() & 0xFFFF_FFFF_0000_0000);
        let rr = (-z * z - 0.5625).exp() * ((z - ax) * (z + ax) + r / ss).exp();
        if hx >= 0 {
            rr / ax
        } else {
            TWO - rr / ax
        }
    } else if hx >= 0 {
        TINY * TINY
    } else {
        TWO - TINY
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_expm1_small() {
        let x = 1e-10;
        let result = expm1(x);
        let expected = x.exp_m1(); // Use Rust's built-in as reference
        assert!(
            (result - expected).abs() < 1e-22,
            "expm1({}) = {}, expected ~{}",
            x,
            result,
            expected
        );
    }

    #[test]
    fn test_expm1_medium() {
        let result = expm1(0.5);
        let expected = 0.5_f64.exp() - 1.0;
        assert!((result - expected).abs() < 1e-14);
    }

    #[test]
    fn test_log1p_small() {
        let x = 1e-10;
        let result = log1p(x);
        assert!(
            (result - x).abs() < 1e-19,
            "log1p({}) = {}, expected ~{}",
            x,
            result,
            x
        );
    }

    #[test]
    fn test_log1p_medium() {
        let result = log1p(0.5);
        let expected = 1.5_f64.ln();
        assert!((result - expected).abs() < 1e-14);
    }

    #[test]
    fn test_erfc_matches_fdlibm() {
        // Reference values from glibc erfc (fdlibm-equivalent).
        // Match to within 1 ULP.
        let cases = [
            (0.0_f64, 1.0_f64),
            (1.0, 0.15729920705028513),
            (-1.0, 1.842700792949715),
            (0.5, 0.4795001221869535),
            // Piecewise boundaries:
            (0.84375, 0.23277433876765835),
            (1.25, 0.07709987174354177),
            // Small |x|:
            (1e-20, 1.0),
            (-1e-20, 1.0),
        ];
        for (x, expected) in cases {
            let got = erfc(x);
            let rel = if expected == 0.0 {
                got.abs()
            } else {
                (got - expected).abs() / expected.abs()
            };
            assert!(
                rel < 1e-14,
                "erfc({}): got {}, expected {}, rel {}",
                x, got, expected, rel
            );
        }
        // Beyond |x| = 6 the fdlibm erfc short-circuits to TINY / 2-TINY.
        assert_eq!(erfc(-6.5), 2.0);
        assert!(erfc(6.5) > 0.0 && erfc(6.5) < 1e-18);
        // erfc(+inf) = 0, erfc(-inf) = 2
        assert_eq!(erfc(f64::INFINITY), 0.0);
        assert_eq!(erfc(f64::NEG_INFINITY), 2.0);
    }

    #[test]
    fn test_gcd_matches_ncbi_semantics() {
        // Standard mathematical cases.
        assert_eq!(gcd(12, 8), 4);
        assert_eq!(gcd(8, 12), 4); // swap handled
        assert_eq!(gcd(7, 5), 1);
        assert_eq!(gcd(0, 5), 5);
        assert_eq!(gcd(100, 0), 100);
        // NCBI `BLAST_Gcd(reward, penalty)` with negative penalty
        // takes abs of second arg before the swap.
        assert_eq!(gcd(3, -12), 3);
        assert_eq!(gcd(5, -15), 5);
    }

    #[test]
    fn test_powi_matches_ncbi_semantics() {
        // Matches `BLAST_Powi`: `x^0 == 1`, `0^n == 0` for n > 0,
        // `0^n == INFINITY` for n < 0.
        assert_eq!(powi(2.0, 0), 1.0);
        assert_eq!(powi(2.0, 3), 8.0);
        assert_eq!(powi(2.0, -2), 0.25);
        assert_eq!(powi(0.0, 5), 0.0);
        assert_eq!(powi(0.0, 0), 1.0);
        assert_eq!(powi(0.0, -1), f64::INFINITY);
        assert_eq!(powi(-3.0, 3), -27.0);
    }

    #[test]
    fn test_nint_matches_ncbi_semantics() {
        // Half-away-from-zero rounding.
        assert_eq!(nint(0.0), 0);
        assert_eq!(nint(0.4), 0);
        assert_eq!(nint(0.5), 1);
        assert_eq!(nint(1.5), 2);
        assert_eq!(nint(2.5), 3);
        assert_eq!(nint(-0.4), 0);
        assert_eq!(nint(-0.5), -1);
        assert_eq!(nint(-1.5), -2);
        assert_eq!(nint(-2.5), -3);
        // Irrational-looking values round normally.
        assert_eq!(nint(3.7), 4);
        assert_eq!(nint(-3.7), -4);
    }

    #[test]
    fn test_romberg_integrates_polynomial_exactly() {
        // Integral of x^2 from 0 to 1 is 1/3. Romberg handles this
        // quickly with trapezoidal → 1/3 after a few refinements.
        let got = romberg_integrate(|x| x * x, 0.0, 1.0, 1e-10, 2, 3);
        assert!(
            (got - 1.0 / 3.0).abs() < 1e-10,
            "got {got}, expected ~0.333"
        );
    }

    #[test]
    fn test_romberg_integrates_exp_smooth() {
        // Integral of exp(-x) from 0 to 5: [-exp(-x)]_0^5 = 1 - exp(-5).
        let got = romberg_integrate(|x| (-x).exp(), 0.0, 5.0, 1e-12, 2, 3);
        let expected = 1.0 - (-5.0_f64).exp();
        assert!((got - expected).abs() < 1e-10, "got {got}, exp {expected}");
    }

    #[test]
    fn test_ln_factorial_small_values_exact() {
        // Every tabulated value must match the closed-form ln(n!) to
        // full f64 precision. The table is there specifically to avoid
        // Stirling's small-n error (was ~1e-4 for n=2).
        assert_eq!(ln_factorial(0), 0.0);
        assert_eq!(ln_factorial(1), 0.0);
        let eps = 1e-14;
        assert!((ln_factorial(2) - 2.0_f64.ln()).abs() < eps);
        assert!((ln_factorial(3) - 6.0_f64.ln()).abs() < eps);
        assert!((ln_factorial(4) - 24.0_f64.ln()).abs() < eps);
        assert!((ln_factorial(5) - 120.0_f64.ln()).abs() < eps);
        assert!((ln_factorial(10) - 3628800.0_f64.ln()).abs() < eps);
        // ln(19!) via product.
        let mut p: f64 = 1.0;
        for i in 2..=19 {
            p *= i as f64;
        }
        assert!((ln_factorial(19) - p.ln()).abs() < eps);
    }

    #[test]
    fn test_ln_factorial_negative_returns_zero() {
        assert_eq!(ln_factorial(-1), 0.0);
        assert_eq!(ln_factorial(-100), 0.0);
    }

    #[test]
    fn test_ln_factorial_large_n_uses_lanczos() {
        // n >= 35 goes through NCBI's Lanczos `s_LnGamma` via
        // `BLAST_LnFactorial`. Should agree with the boundary table
        // entry via ln(n!) = ln((n-1)!) + ln(n).
        let a = ln_factorial(35);
        let b = ln_factorial(34) + 35.0_f64.ln();
        assert!(
            (a - b).abs() < 1e-10,
            "lanczos(35) vs table(34)+ln(35) diverged: {a} vs {b}"
        );
    }

    #[test]
    fn test_ln_gamma_int_matches_ln_factorial_offset() {
        // NCBI `BLAST_LnGammaInt(n)` == ln(Gamma(n)) == ln((n-1)!).
        for n in 1..=50 {
            let got = ln_gamma_int(n);
            let expected = ln_factorial(n - 1);
            assert!(
                (got - expected).abs() < 1e-10,
                "ln_gamma_int({n}) = {got}, expected ln((n-1)!) = {expected}"
            );
        }
    }

    #[test]
    fn test_romberg_respects_itmin() {
        // Tiny eps but large itmin: converges trivially on constant fn.
        let got = romberg_integrate(|_| 1.0, 0.0, 1.0, 1e-15, 1, 5);
        assert!((got - 1.0).abs() < 1e-15);
    }
}
