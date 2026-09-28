# kernel/utilities/tikhoind.m

- Signature: `[x,err,reg]=tikhoind(K,D,y,lam)`

## Purpose

Computes an unconstrained, sign-indefinite Tikhonov-regularised solution of `K*x=y`.

## Parameters / inputs

- `K` — kernel matrix; may be complex or non-square.
- `D` — regularisation matrix.
- `y` — column vector; may be complex.
- `lam` — nonnegative real scalar Tikhonov regularisation parameter.

## Outputs

- `x` — solution minimising `norm(K*x-y,2)^2 + lam*norm(D*x,2)^2` as computed by the normal equations. The source describes `x` as real, but the implementation does not enforce this; complex inputs may produce a complex solution.
- `err` — squared residual norm, `norm(K*x-y,2)^2`.
- `reg` — squared regularisation norm, `norm(D*x,2)^2`.

## Numerical / algorithmic content

The solution is computed as `x = (K'*K + lam*(D'*D)) \ (K'*y)`. For best numerical performance, scale `K` to have approximately unit 2-norm and `y` to have approximately unit 1-norm. See `tikhonov.m` for the positive-constrained solver.

The implementation checks that all inputs are numeric, that `size(K,1) == size(y,1)`, and that `lam` is a nonnegative real scalar. It computes `err` and `reg` only when those outputs are requested. The source does not explicitly check the dimensions of `D` or enforce that `y` is a column vector.

<https://spindynamics.org/wiki/index.php?title=tikhoind.m>