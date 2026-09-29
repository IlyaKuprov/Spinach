# kernel/utilities/svd_shrink.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/svd_shrink.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/svd_shrink.m)

## Purpose

Generates sets of vector-covector pairs for the parallel implementation of the time propagation algorithm described in [http://dx.doi.org/10.1063/1.3679656](http://dx.doi.org/10.1063/1.3679656) (Equation 9). The function decomposes a density matrix into a low-rank vector-covector representation by discarding singular values below a user-specified tolerance.

## Behaviour

- Syntax: `[vec,cov]=svd_shrink(spin_system,rho,tol)`.
- The function first validates its inputs via an internal consistency check (`grumble`):
  - `rho` must be a numeric square matrix, otherwise an error `'rho must be a square matrix.'` is raised.
  - `tol` must be a finite, non-negative, real scalar, otherwise an error `'tol must be a finite non-negative real scalar.'` is raised.
- The density matrix is converted to a full matrix and decomposed via the singular value decomposition `[vec,S,cov]=svd(full(rho))`, with singular values taken from `diag(S)`.
- A drop mask is formed from singular values satisfying `S<tol`.
- The user is informed through a report message stating how many insignificant vector-covector pairs were dropped, e.g. `'dropped N insignificant vector-covector pairs from the density matrix.'` where `N` is `nnz(drop_mask)`.
- Columns of `vec` and `cov` corresponding to dropped singular values are removed, and the corresponding singular values are deleted from `S`.
- The remaining singular values are spread across both factors: `vec=vec*diag(sqrt(S))` and `cov=cov*diag(sqrt(S))`, so that `vec*cov'` reconstructs the truncated density matrix.

## Inputs and outputs

**Inputs**

- `spin_system` — Spinach spin system object, used for user reporting.
- `rho` — density matrix (numeric square matrix).
- `tol` — singular value drop tolerance (finite non-negative real scalar).

**Outputs**

- `vec` — vectors as columns of a matrix.
- `cov` — covectors as columns of a matrix.

## References

1. Source code: [https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/svd_shrink.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/svd_shrink.m)
2. Spinach Wiki: [https://spindynamics.org/wiki/index.php?title=svd_shrink.m](https://spindynamics.org/wiki/index.php?title=svd_shrink.m)
3. Time propagation algorithm (Equation 9): [http://dx.doi.org/10.1063/1.3679656](http://dx.doi.org/10.1063/1.3679656)
