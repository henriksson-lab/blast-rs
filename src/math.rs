//! Rust equivalent of ncbi_math.c — math utility functions for BLAST.

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
    x * (1.0 + x *
        (1.0/2.0 + x *
        (1.0/6.0 + x *
        (1.0/24.0 + x *
        (1.0/120.0 + x *
        (1.0/720.0 + x *
        (1.0/5040.0 + x *
        (1.0/40320.0 + x *
        (1.0/362880.0 + x *
        (1.0/3628800.0 + x *
        (1.0/39916800.0 + x *
        (1.0/479001600.0 +
         x/6227020800.0))))))))))))
}

/// ln(1+x) for small x (avoids catastrophic cancellation).
pub fn log1p(x: f64) -> f64 {
    if x.abs() >= 0.2 {
        return (x + 1.0).ln();
    }
    // Taylor series: ln(1+x) = x - x^2/2 + x^3/3 - x^4/4 + ...
    let mut sum = 0.0;
    let mut power = x; // x^i
    for i in 1..=500 {
        let term = power / i as f64;
        if i % 2 == 1 {
            sum += term;
        } else {
            sum -= term;
        }
        if term.abs() < f64::EPSILON {
            break;
        }
        power *= x;
    }
    sum
}

/// Log-sum of two values: ln(e^x + e^y) computed stably.
pub fn log_sum(x: f64, y: f64) -> f64 {
    let (big, small) = if x > y { (x, y) } else { (y, x) };
    if big - small > 100.0 {
        return big; // small contribution is negligible
    }
    big + log1p((small - big).exp())
}

/// Compute the natural log of the factorial of n using Stirling's approximation.
pub fn ln_factorial(n: i32) -> f64 {
    if n < 0 {
        return 0.0;
    }
    if n <= 1 {
        return 0.0;
    }
    // Use lgamma for larger values
    let n = n as f64;
    // Stirling's approximation + correction terms
    let ln2pi_half = 0.5 * (2.0 * std::f64::consts::PI).ln();
    ln2pi_half + (n + 0.5) * n.ln() - n
        + 1.0 / (12.0 * n)
        - 1.0 / (360.0 * n * n * n)
}

/// Compute n choose k (binomial coefficient) in log space.
pub fn ln_choose(n: i32, k: i32) -> f64 {
    if k < 0 || k > n {
        return 0.0;
    }
    ln_factorial(n) - ln_factorial(k) - ln_factorial(n - k)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_expm1_small() {
        let x = 1e-10;
        let result = expm1(x);
        let expected = x.exp_m1(); // Use Rust's built-in as reference
        assert!((result - expected).abs() < 1e-22, "expm1({}) = {}, expected ~{}", x, result, expected);
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
        assert!((result - x).abs() < 1e-19, "log1p({}) = {}, expected ~{}", x, result, x);
    }

    #[test]
    fn test_log1p_medium() {
        let result = log1p(0.5);
        let expected = 1.5_f64.ln();
        assert!((result - expected).abs() < 1e-14);
    }
}
