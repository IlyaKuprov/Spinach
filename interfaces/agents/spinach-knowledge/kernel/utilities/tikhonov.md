# kernel/utilities/tikhonov.m

- Signature: `[x,err,reg]=tikhonov(K,D,KtK,DtD,H,y,lambda)`

## Purpose

Find a nonnegative, real Tikhonov-regularised solution to `K*x=y` by minimising `norm(K*x-y,2)^2+lambda*norm(D*x,2)^2` subject to `x>=0`. The kernel and data may be complex.

## Parameters / inputs

- `K`: Kernel matrix; may be complex or non-square.
- `D`: Regularisation matrix. Leave empty to use `fdmat(size(K,2),5,2,'wall')`, a finite-difference second-derivative matrix.
- `KtK`: Precomputed `K'*K`; leave empty to compute it. Supplying it can speed up repeated calls.
- `DtD`: Precomputed `D'*D`; leave empty to compute it. Supplying it can speed up repeated calls.
- `H`: Precomputed Tikhonov Hessian `2*real(KtK+lambda*DtD)`; leave empty to compute it. Supplying it can speed up repeated calls.
- `y`: Column vector; may be complex.
- `lambda`: Nonnegative real scalar Tikhonov regularisation parameter.

## Outputs

- `x`: Real, nonnegative vector minimising the regularised objective subject to the positivity constraint.
- `err`: Squared residual norm `norm(K*x-y,2)^2`.
- `reg`: Squared regularisation norm `norm(D*x,2)^2`.

## Numerical / algorithmic content

The function uses `fmincon` with its interior-point algorithm, an initial vector of ones, lower bounds of zero, and upper bounds of infinity. It supplies the objective gradient `2*real(KtK*x-K'*y+lambda*DtD*x)` and the Hessian `H`. The optimisation is configured for at most 100 iterations, with an unlimited number of function evaluations and iteration-level display.

For best numerical performance, scale `K` to approximately unit 2-norm and `y` to approximately unit 1-norm. See `tikhoind.m` for the indeterminate solver.

Source reference: https://spindynamics.org/wiki/index.php?title=tikhonov.m