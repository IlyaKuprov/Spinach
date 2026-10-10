# kernel/utilities/arnoldi.m

Source: [kernel/utilities/arnoldi.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/arnoldi.m)

## Purpose

Generates an orthonormal Krylov basis from repeated action of an operator on a starting vector, using the Arnoldi procedure. The header comment notes that the procedure is numerically unstable and must be used with caution.

## Behaviour

- Syntax: `[V,H]=arnoldi(Op,v0,niter)`.
- The first basis vector is `v0` normalised by its 2-norm.
- Each iteration applies `Op` to the latest basis vector and orthogonalises the result against all previous basis vectors via a classical Gram–Schmidt loop, storing the projections in `H`.
- The subdiagonal entry `H(n+1,n)` is the 2-norm of the orthogonalised vector, which is then normalised to become the next basis vector.
- If an exact Krylov breakdown occurs (`H(n+1,n)==0`), `V` and `H` are truncated to the completed invariant subspace (`V=V(:,1:n)`, `H=H(1:n,1:n)`) and the function returns.
- `V` and `H` are preallocated as complex (`'like',1i`) with sizes `numel(v0)`-by-`niter+1` and `niter+1`-by-`niter` respectively.
- Input validation is performed by the internal `grumble` function, which errors when: `Op` is not a function handle; `v0` is not a numeric column vector; or `niter` is not a non-negative real integer scalar.

## Inputs and outputs

Inputs:

- `Op` — function handle taking a column vector and returning another column vector.
- `v0` — starting vector of the Arnoldi process.
- `niter` — number of iterations to take; the Krylov subspace will be `niter+1` dimensional.

Outputs:

- `V` — matrix containing the orthonormal basis vectors of the Krylov subspace in columns.
- `H` — extended Hessenberg matrix.

## References

- Spinach Wiki: [arnoldi.m](https://spindynamics.org/wiki/index.php?title=arnoldi.m)
