# kernel/optimcon/lbfgs.m

- Signature: `direction=lbfgs(dx_hist,dg_hist,g)`

## Purpose

Computes an L-BFGS ascent direction from the current objective gradient and stored parameter-step and gradient-change columns. It does not form a Hessian matrix and does not perform a line search.

## Inputs and guards

All three arguments are required; the source defines no defaults. `g` must be a finite real column vector. `dx_hist` and `dg_hist` are column-stacked histories, with corresponding columns and the same row dimension as `g` when the history is non-empty; both must be finite and real. The source checks that the history arrays have equal column counts and that `g` has one column. Empty histories are allowed and return `g` directly.

History columns are processed newest to oldest. Before the two-loop recursion, each pair is retained only if its inner products are finite, both squared norms are positive, and

`dg_hist(:,n)'*dx_hist(:,n) < -0.01*sqrt(sum(dg_hist(:,n).^2)*sum(dx_hist(:,n).^2))`.

The 0.01 relative-curvature threshold is fixed in the function, not an input option. If every pair is rejected, the function returns the current gradient as the steepest-ascent direction. Otherwise it applies the L-BFGS two-loop recursion, scales its initial inverse-Hessian approximation using the newest retained pair, and negates the minimisation-form result to obtain an ascent direction. The number of history columns controls the stored history; there is no separate memory-size or line-search option here.

No control-freeze or phase-cycle mask is applied by this routine.

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/optimcon/lbfgs.m)
[Spinach Wiki page](https://spindynamics.org/wiki/index.php?title=lbfgs.m)
