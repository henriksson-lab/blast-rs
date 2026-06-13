//! Verbatim port of NCBI optimize_target_freq.c.
//! Finds optimal target frequencies for compositionally adjusted score matrices.
//!
//! The optimal target frequencies x minimize the Kullback-Liebler "distance"
//! `sum_k x[k] * ln(x[k]/q[k])` subject to marginal sum constraints and
//! optionally a relative entropy constraint.

use crate::algo::blast::composition_adjustment::nlm_linear_algebra::{
    add_vectors, euclidean_norm, factor_ltriag_pos_def, solve_ltriag_pos_def, step_bound,
};

/// Port of NCBI ScaledSymmetricProductA.
/// Compute A D A^T where A is the constraint matrix (used implicitly).
pub fn scaled_symmetric_product_a(w: &mut [Vec<f64>], diagonal: &[f64], alphsize: usize) {
    let m = 2 * alphsize - 1;
    for row_w in 0..m {
        for col_w in 0..=row_w {
            w[row_w][col_w] = 0.0;
        }
    }
    for i in 0..alphsize {
        for j in 0..alphsize {
            let dd = diagonal[i * alphsize + j];
            w[j][j] += dd;
            if i > 0 {
                w[i + alphsize - 1][j] += dd;
                w[i + alphsize - 1][i + alphsize - 1] += dd;
            }
        }
    }
}

/// Port of NCBI MultiplyByA.
/// y = beta * y + alpha * A * x.
pub fn multiply_by_a(beta: f64, y: &mut [f64], alphsize: usize, alpha: f64, x: &[f64]) {
    if beta == 0.0 {
        for i in 0..(2 * alphsize - 1) {
            y[i] = 0.0;
        }
    } else if beta != 1.0 {
        for i in 0..(2 * alphsize - 1) {
            y[i] *= beta;
        }
    }
    for i in 0..alphsize {
        for j in 0..alphsize {
            y[j] += alpha * x[i * alphsize + j];
        }
    }
    for i in 1..alphsize {
        for j in 0..alphsize {
            y[i + alphsize - 1] += alpha * x[i * alphsize + j];
        }
    }
}

/// Port of NCBI MultiplyByAtranspose.
/// y = beta * y + alpha * A^T * x.
pub fn multiply_by_a_transpose(beta: f64, y: &mut [f64], alphsize: usize, alpha: f64, x: &[f64]) {
    if beta == 0.0 {
        for k in 0..(alphsize * alphsize) {
            y[k] = 0.0;
        }
    } else if beta != 1.0 {
        for k in 0..(alphsize * alphsize) {
            y[k] *= beta;
        }
    }
    for i in 0..alphsize {
        for j in 0..alphsize {
            let k = i * alphsize + j;
            y[k] += alpha * x[j];
            if i > 0 {
                y[k] += alpha * x[i + alphsize - 1];
            }
        }
    }
}

/// Port of NCBI ResidualsLinearConstraints.
pub fn residuals_linear_constraints(
    r_a: &mut [f64],
    alphsize: usize,
    x: &[f64],
    row_sums: &[f64],
    col_sums: &[f64],
) {
    // NCBI `optimize_target_freq.c:218-226` uses explicit index loops
    // here; keep line-for-line parity rather than `copy_from_slice`.
    #[allow(clippy::manual_memcpy)]
    for i in 0..alphsize {
        r_a[i] = col_sums[i];
    }
    for i in 1..alphsize {
        r_a[i + alphsize - 1] = row_sums[i];
    }
    multiply_by_a(1.0, r_a, alphsize, -1.0, x);
}

/// Port of NCBI DualResiduals.
pub fn dual_residuals(
    resids_x: &mut [f64],
    alphsize: usize,
    grads: &[Vec<f64>],
    z: &[f64],
    constrain_rel_entropy: bool,
) {
    let n = alphsize * alphsize;
    if constrain_rel_entropy {
        let eta = z[2 * alphsize - 1];
        for i in 0..n {
            resids_x[i] = -grads[0][i] + eta * grads[1][i];
        }
    } else {
        for i in 0..n {
            resids_x[i] = -grads[0][i];
        }
    }
    multiply_by_a_transpose(1.0, resids_x, alphsize, 1.0, z);
}

/// Port of NCBI CalculateResiduals.
pub fn calculate_residuals(
    resids_x: &mut [f64],
    alphsize: usize,
    resids_z: &mut [f64],
    values: &[f64],
    grads: &[Vec<f64>],
    row_sums: &[f64],
    col_sums: &[f64],
    x: &[f64],
    z: &[f64],
    constrain_rel_entropy: bool,
    relative_entropy: f64,
) -> f64 {
    dual_residuals(resids_x, alphsize, grads, z, constrain_rel_entropy);
    let norm_resids_x = euclidean_norm(resids_x, alphsize * alphsize);

    residuals_linear_constraints(resids_z, alphsize, x, row_sums, col_sums);

    let norm_resids_z = if constrain_rel_entropy {
        resids_z[2 * alphsize - 1] = relative_entropy - values[1];
        euclidean_norm(resids_z, 2 * alphsize)
    } else {
        euclidean_norm(resids_z, 2 * alphsize - 1)
    };

    (norm_resids_x * norm_resids_x + norm_resids_z * norm_resids_z).sqrt()
}

/// Port of NCBI EvaluateReFunctions.
pub fn evaluate_re_functions(
    values: &mut [f64],
    grads: &mut [Vec<f64>],
    alphsize: usize,
    x: &[f64],
    q: &[f64],
    scores: &[f64],
    constrain_rel_entropy: bool,
) {
    values[0] = 0.0;
    values[1] = 0.0;
    for k in 0..(alphsize * alphsize) {
        let temp = (x[k] / q[k]).ln();
        values[0] += x[k] * temp;
        grads[0][k] = temp + 1.0;

        if constrain_rel_entropy {
            let temp2 = temp + scores[k];
            values[1] += x[k] * temp2;
            grads[1][k] = temp2 + 1.0;
        }
    }
}

/// Port of NCBI ComputeScoresFromProbs.
pub fn compute_scores_from_probs(
    scores: &mut [f64],
    alphsize: usize,
    target_freqs: &[f64],
    row_freqs: &[f64],
    col_freqs: &[f64],
) {
    for i in 0..alphsize {
        for j in 0..alphsize {
            let k = i * alphsize + j;
            scores[k] = (target_freqs[k] / (row_freqs[i] * col_freqs[j])).ln();
        }
    }
}

/// Newton system for the optimization.
/// Port of NCBI ReNewtonSystem.
pub struct ReNewtonSystem {
    alphsize: usize,
    constrain_rel_entropy: bool,
    w: Vec<Vec<f64>>,  // lower-triangular factor
    dinv: Vec<f64>,    // diagonal inverse
    grad_re: Vec<f64>, // gradient of RE constraint
}

impl ReNewtonSystem {
    /// Port of NCBI ReNewtonSystemNew.
    fn new(alphsize: usize) -> Self {
        let m = 2 * alphsize;
        // Lower-triangular storage: row i has (i+1) elements
        let mut w = Vec::with_capacity(m);
        for i in 0..m {
            w.push(vec![0.0f64; i + 1]);
        }
        let n = alphsize * alphsize;
        ReNewtonSystem {
            alphsize,
            constrain_rel_entropy: true,
            w,
            dinv: vec![0.0; n],
            grad_re: vec![0.0; n],
        }
    }

    /// Port of NCBI FactorReNewtonSystem.
    fn factor(
        &mut self,
        x: &[f64],
        z: &[f64],
        grads: &[Vec<f64>],
        constrain_rel_entropy: bool,
        workspace: &mut [f64],
    ) {
        let alphsize = self.alphsize;
        let n = alphsize * alphsize;
        let m = if constrain_rel_entropy {
            2 * alphsize
        } else {
            2 * alphsize - 1
        };
        self.constrain_rel_entropy = constrain_rel_entropy;

        // Find D^{-1}
        if constrain_rel_entropy {
            let eta = z[m - 1];
            for i in 0..n {
                self.dinv[i] = x[i] / (1.0 - eta);
            }
        } else {
            self.dinv[..n].copy_from_slice(&x[..n]);
        }

        // Compute J D^{-1} J^T using the constraint structure
        scaled_symmetric_product_a(&mut self.w, &self.dinv, alphsize);

        if constrain_rel_entropy {
            self.grad_re[..n].copy_from_slice(&grads[1][..n]);

            self.w[m - 1][m - 1] = 0.0;
            for i in 0..n {
                workspace[i] = self.dinv[i] * self.grad_re[i];
                self.w[m - 1][m - 1] += self.grad_re[i] * workspace[i];
            }
            // MultiplyByA(0.0, &W[m-1][0], alphsize, 1.0, workspace)
            // We need to write into w[m-1] but only the first (2*alphsize-1) elements
            let w_row = &mut self.w[m - 1];
            // Initialize the first (2*alphsize-1) elements to zero, keeping [m-1] intact
            for idx in 0..(2 * alphsize - 1) {
                w_row[idx] = 0.0;
            }
            // A * workspace: first alphsize elements are column sums
            for i in 0..alphsize {
                for j in 0..alphsize {
                    w_row[j] += workspace[i * alphsize + j];
                }
            }
            // Row sums for i >= 1
            for i in 1..alphsize {
                for j in 0..alphsize {
                    w_row[i + alphsize - 1] += workspace[i * alphsize + j];
                }
            }
        }
        factor_ltriag_pos_def(&mut self.w, m);
    }

    /// Port of NCBI SolveReNewtonSystem.
    fn solve(&self, x: &mut [f64], z: &mut [f64], workspace: &mut [f64]) {
        let alphsize = self.alphsize;
        let n = alphsize * alphsize;
        let m_a = 2 * alphsize - 1;
        let m = if self.constrain_rel_entropy {
            m_a + 1
        } else {
            m_a
        };

        // rzhat = rz - J D^{-1} rx
        for i in 0..n {
            workspace[i] = x[i] * self.dinv[i];
        }
        multiply_by_a(1.0, z, alphsize, -1.0, workspace);

        if self.constrain_rel_entropy {
            for i in 0..n {
                z[m - 1] -= self.grad_re[i] * workspace[i];
            }
        }

        // Solve using the factored matrix
        solve_ltriag_pos_def(z, m, &self.w);

        // Backsolve for x
        if self.constrain_rel_entropy {
            for i in 0..n {
                x[i] += self.grad_re[i] * z[m - 1];
            }
        }
        multiply_by_a_transpose(1.0, x, alphsize, 1.0, z);

        for i in 0..n {
            x[i] *= self.dinv[i];
        }
    }
}

/// Port of NCBI Blast_OptimizeTargetFrequencies.
/// Find optimal target frequencies that minimize KL-distance from q
/// subject to marginal sum constraints and optionally relative entropy.
///
/// Returns 0 on success (converged to minimizer), 1 on convergence failure.
pub fn optimize_target_frequencies(
    x: &mut [f64], // output: optimal target freqs (n = alphsize^2)
    alphsize: usize,
    q: &[f64],        // standard matrix target freqs
    row_sums: &[f64], // desired row marginals
    col_sums: &[f64], // desired column marginals
    constrain_rel_entropy: bool,
    relative_entropy: f64,
    tol: f64,
    maxits: usize,
) -> (i32, usize) {
    let n = alphsize * alphsize;
    let m_a = 2 * alphsize - 1;
    let m = if constrain_rel_entropy { m_a + 1 } else { m_a };

    let mut newton_system = ReNewtonSystem::new(alphsize);
    let mut resids_x = vec![0.0f64; n];
    let mut resids_z = vec![0.0f64; m_a + 1];
    let mut z = vec![0.0f64; m_a + 1]; // must be initialized to zero
    let mut old_scores = vec![0.0f64; n];
    let mut workspace = vec![0.0f64; n];
    let mut grads = vec![vec![0.0f64; n], vec![0.0f64; n]];
    let mut values = [0.0f64; 2];

    compute_scores_from_probs(&mut old_scores, alphsize, q, row_sums, col_sums);

    // Use q as initial value for x
    x[..n].copy_from_slice(&q[..n]);

    let mut its = 0usize;
    let mut rnorm = f64::MAX;

    while its <= maxits {
        evaluate_re_functions(
            &mut values,
            &mut grads,
            alphsize,
            x,
            q,
            &old_scores,
            constrain_rel_entropy,
        );
        rnorm = calculate_residuals(
            &mut resids_x,
            alphsize,
            &mut resids_z,
            &values,
            &grads,
            row_sums,
            col_sums,
            x,
            &z,
            constrain_rel_entropy,
            relative_entropy,
        );

        // NCBI `optimize_target_freq.c:527` uses `!(rnorm > tol)` rather
        // than `rnorm <= tol` so NaN breaks the loop (NaN > tol is false).
        #[allow(clippy::neg_cmp_op_on_partial_ord)]
        if !(rnorm > tol) {
            break;
        } else {
            its += 1;
            if its <= maxits {
                newton_system.factor(x, &z, &grads, constrain_rel_entropy, &mut workspace);
                newton_system.solve(&mut resids_x, &mut resids_z, &mut workspace);

                let mut alpha = step_bound(x, n, &resids_x, 1.0 / 0.95);
                alpha *= 0.95;

                add_vectors(x, n, alpha, &resids_x);
                add_vectors(&mut z, m, alpha, &resids_z);
            }
        }
    }

    let converged = its <= maxits && rnorm <= tol && (!constrain_rel_entropy || z[m - 1] < 1.0);

    let status = if converged { 0 } else { 1 };
    (status, its)
}
