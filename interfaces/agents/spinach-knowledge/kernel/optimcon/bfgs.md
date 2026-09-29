# kernel/optimcon/bfgs.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/optimcon/bfgs.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=bfgs.m)

## Purpose

Build a dense BFGS approximation to the negative Hessian for maximising an objective from a history of argument and gradient increments. The corresponding Newton-like ascent direction `p` is defined by `H*p=g`, where `g` is the current gradient. This routine uses gradient history only; it does not evaluate the objective or impose constraints.

## Syntax

`H=bfgs(dx_hist,dg_hist,g)`

## Inputs

- `dx_hist` — `n`-by-`N` array whose columns are argument increments, ordered newest to oldest.
- `dg_hist` — matching `n`-by-`N` array of gradient increments in the same order.
- `g` — current real finite gradient column vector. Its length sets the matrix dimension; the update uses the recorded history rather than the values of `g`.

The implementation checks that both histories have the same number of columns, `g` is a column, nonempty history arrays have the same row count as `g`, and all three inputs are real and finite. Empty histories are allowed.

## Output and update

- `H` — real symmetric `n`-by-`n` dense approximation to the negative objective Hessian. If no history pair passes the curvature filter, the result is the identity matrix.

The source filters history pairs using `dg_hist(:,i)' * dx_hist(:,i) < -0.01*norm(dg_hist(:,i))*norm(dx_hist(:,i))`, also requiring nonzero increments. For the BFGS update it reverses the gradient increment, using `y=-dg_hist(:,i)`. The initial scale is based on the newest retained pair; the remaining retained pairs are applied in chronological order. Unsafe update denominators are skipped, and each completed update is made real and symmetric.

No physical units are assigned by this routine; the increments and gradient use the caller's optimisation coordinates and objective.
