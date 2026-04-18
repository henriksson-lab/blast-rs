//! Verbatim port of NCBI nlm_linear_algebra.c.
//! Basic matrix and vector operations for composition adjustment.

/// Port of NCBI Nlm_FactorLtriangPosDef.
/// Cholesky factorization of a symmetric positive-definite matrix
/// stored as lower-triangular row pointers.
/// The matrix A is modified in-place; on exit `A[i][j]` for `j <= i`
/// holds the lower Cholesky factor L.
pub fn factor_ltriag_pos_def(a: &mut [Vec<f64>], n: usize) {
    for i in 0..n {
        for j in 0..i {
            let mut temp = a[i][j];
            for k in 0..j {
                temp -= a[i][k] * a[j][k];
            }
            a[i][j] = temp / a[j][j];
        }
        let mut temp = a[i][i];
        for k in 0..i {
            temp -= a[i][k] * a[i][k];
        }
        a[i][i] = temp.sqrt();
    }
}

/// Port of NCBI Nlm_SolveLtriangPosDef.
/// Solve L L^T x = b where L is the Cholesky factor.
/// On entry x = b; on exit x = solution.
pub fn solve_ltriag_pos_def(x: &mut [f64], n: usize, l: &[Vec<f64>]) {
    // Forward solve: L z = b
    for i in 0..n {
        let mut temp = x[i];
        for j in 0..i {
            temp -= l[i][j] * x[j];
        }
        x[i] = temp / l[i][i];
    }
    // Back solve: L^T y = z
    for j in (0..n).rev() {
        x[j] /= l[j][j];
        for i in 0..j {
            x[i] -= l[j][i] * x[j];
        }
    }
}

/// Port of NCBI Nlm_EuclideanNorm.
/// Scaled Euclidean norm to avoid overflow.
pub fn euclidean_norm(v: &[f64], n: usize) -> f64 {
    let mut sum = 1.0f64;
    let mut scale = 0.0f64;
    for i in 0..n {
        if v[i] != 0.0 {
            let absvi = v[i].abs();
            if scale < absvi {
                sum = 1.0 + sum * (scale / absvi) * (scale / absvi);
                scale = absvi;
            } else {
                sum += (absvi / scale) * (absvi / scale);
            }
        }
    }
    scale * sum.sqrt()
}

/// Port of NCBI Nlm_AddVectors.
/// `y[i] += alpha * x[i]` for `i` in `0..n`.
pub fn add_vectors(y: &mut [f64], n: usize, alpha: f64, x: &[f64]) {
    for i in 0..n {
        y[i] += alpha * x[i];
    }
}

/// Port of NCBI Nlm_StepBound.
/// Find the largest step alpha in [0, max] such that x + alpha*step_x >= 0.
pub fn step_bound(x: &[f64], n: usize, step_x: &[f64], max: f64) -> f64 {
    let mut alpha = max;
    for i in 0..n {
        let alpha_i = -x[i] / step_x[i];
        if alpha_i >= 0.0 && alpha_i < alpha {
            alpha = alpha_i;
        }
    }
    alpha
}
