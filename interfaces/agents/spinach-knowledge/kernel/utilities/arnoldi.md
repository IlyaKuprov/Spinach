# kernel/utilities/arnoldi.m

- Signature: `[V,H]=arnoldi(Op,v0,niter)`

## Purpose

Construct an orthonormal Krylov basis by repeatedly applying an operator to a starting vector. The source cautions that this Arnoldi implementation is numerically unstable and should be used with care.

## Numerical content

Starting with the normalized `v0`, each iteration applies `Op` to the latest basis vector, orthogonalizes the result against all previously computed basis vectors with Gram-Schmidt, and records the coefficients in an extended Hessenberg matrix. If the residual norm is exactly zero, the function returns the completed invariant subspace with `V` and `H` truncated.

## Parameters / inputs

- `Op` - function handle that accepts a column vector and returns a column vector.
- `v0` - numeric column vector that starts the Arnoldi process.
- `niter` - non-negative integer number of iterations. Without exact breakdown, the Krylov basis has `niter+1` columns.

## Outputs

- `V` - matrix whose columns are the computed orthonormal Krylov basis vectors.
- `H` - extended Hessenberg matrix of Arnoldi coefficients.
