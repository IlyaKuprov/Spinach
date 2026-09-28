# kernel/optimcon/lbfgs.m

- Signature: `direction=lbfgs(dx_hist,dg_hist,g)`

## Purpose

Computes a limited-memory BFGS approximation to the Newton–Raphson ascent direction without forming a Hessian matrix.

## Parameters / inputs

- `dx_hist` — columns of recent parameter-step vectors, newest first.
- `dg_hist` — matching gradient-change vectors, newest first.
- `g` — current real gradient column vector.

All vectors must be finite and real. The history arrays must have the same number of columns as each other and the same number of rows as `g`.

## Output

- `direction` — approximate maximisation step. The implementation filters out history pairs that fail its curvature test; if no pairs remain, it returns `g` (steepest ascent).

## Implementation

Uses the L-BFGS two-loop recursion over the retained history, scales its initial inverse-Hessian approximation from the newest retained pair, and negates the minimisation-form recursion to produce an ascent direction.

[Spinach Wiki page](https://spindynamics.org/wiki/index.php?title=lbfgs.m)
