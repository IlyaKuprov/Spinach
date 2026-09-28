# kernel/utilities/tikhol1n.m

- Signature: `[x,err,reg]=tikhol1n(A,y,nnzt)`

## Purpose

Solves the ill-conditioned system `A*x=y` using L1-norm Tikhonov regularisation. It minimises `norm(A*x-y,2)^2+lambda*norm(x,1)` using FISTA. The user specifies the desired number of non-zero entries in the solution, and the regularisation parameter `lambda` is found by bracketing and bisection.

## Parameters / inputs

- `A` — a real or complex matrix.
- `y` — a real or complex column vector.
- `nnzt` — the target number of non-zero entries in the solution.

## Outputs

- `x` — a real or complex solution vector.
- `err` — the squared 2-norm of the fitting error divided by the squared 2-norm of the solution.
- `reg` — the 1-norm of the solution.