# kernel/utilities/tikhol1n.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/tikhol1n.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/tikhol1n.m)

## Purpose

L1-norm Tikhonov regularised solver for `A*x=y` where `A` is an ill-conditioned matrix. The error functional `norm(A*x-y,2)^2 + lambda*norm(x,1)` is minimised using the FISTA algorithm. The user specifies the desired number of non-zeroes; the lambda parameter is then found by bracketing and bisection on the soft-thresholding value.

## Behaviour

- Syntax: `[x,err,reg]=tikhol1n(A,y,nnzt)`.
- Input consistency is enforced by an internal `grumble` function: `A` must be a numeric matrix, `y` a numeric column vector with as many elements as rows of `A`, and `nnzt` a positive integer scalar not exceeding the number of columns of `A`.
- The solver pre-computes `A'` (conjugate transpose), estimates the Lipschitz constant as `L = 2*(1+normest_tol)*normest(A,normest_tol)^2` with `normest_tol = 1e-3`, and sets the initial soft threshold `thr = 1/L`.
- Soft thresholding is applied as `sign(x).*max(abs(x)-thr,0)`, which handles real and complex inputs.
- FISTA iterations start from a zero initial guess with Nesterov momentum `t = 1`; the momentum is updated as `t_new = 0.5*(1+sqrt(1+4*t^2))` and the step is `x = x_prox + ((t-1)/t_new)*(x_prox-x_old)`.
- Convergence is declared when the relative step norm `norm(x_prox-x_old,2)/norm(x_prox,2)` falls below `step_norm_tol = 1e-6`.
- The zero threshold is bracketed between `thr_lower = 0` and `thr_upper = max(abs(A'*y))`. Every 1000 iterations, or on convergence, the solution's non-zero count is checked against the target `nnzt` with a tolerance of ±1. If the count is outside `[nnzt-1, nnzt, nnzt+1]` (or the solution is all zeros), the bracket is updated: `thr_lower = thr` when `nnz(x) > nnzt` (fewer non-zeroes needed), otherwise `thr_upper = thr`; the threshold is reset to the bracket midpoint and FISTA restarts with `t = 1`.
- Stagnation is detected when `(thr_upper-thr_lower)/(thr_upper+thr_lower) < 1e-6`, in which case the function errors with `nnz target unreachable`.
- Progress reports are printed at every 1000 iterations and at convergence, including iteration count, `nnz(x)`, relative squared error, 1-norm, relative step norm, and the current zero threshold.

## Inputs and outputs

**Inputs**

- `A` — a real or complex matrix.
- `y` — a real or complex column vector.
- `nnzt` — the target for the number of non-zeroes in the solution.

**Outputs**

- `x` — a real or complex vector.
- `err` — squared 2-norm of the fitting error divided by the squared 2-norm of the solution.
- `reg` — 1-norm of the solution.

## References

- Spin Dynamics Wiki page for this function: [https://spindynamics.org/wiki/index.php?title=tikhol1n.m](https://spindynamics.org/wiki/index.php?title=tikhol1n.m)
