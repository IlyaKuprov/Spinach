# kernel/utilities/svd_shrink.m

- Signature: `[vec,cov]=svd_shrink(spin_system,rho,tol)`

## Purpose

Generates vector-covector pairs for the parallel implementation of the time propagation algorithm described in Equation 9 of http://dx.doi.org/10.1063/1.3679656.

## Parameters / inputs

- `spin_system` — spin system used to report the number of dropped pairs.
- `rho` — density matrix; must be a numeric square matrix.
- `tol` — singular value drop tolerance; must be a finite, non-negative real scalar.

## Outputs

- `vec` — vectors as columns of a matrix.
- `cov` — covectors as columns of a matrix.

## Numerical / algorithmic content

The function computes the singular value decomposition of `full(rho)`, then removes columns associated with singular values strictly below `tol`. It reports the number of dropped vector-covector pairs. Each retained vector and covector column is multiplied by the square root of its corresponding singular value, spreading that coefficient across the pair. With MATLAB's complex-conjugate transpose, the retained factors reconstruct the truncated density matrix as `vec*cov'`.

Source documentation: https://spindynamics.org/wiki/index.php?title=svd_shrink.m