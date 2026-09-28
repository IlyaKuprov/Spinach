# kernel/derivatives/fdkup.m

- Signature: `K=fdkup(npoints,extents,chi,nstenc)`

## Purpose

Returns a finite-difference representation of the Kuprov operator, acting on a three-dimensional array `rho` with axes ordered `[X Y Z]`:

`K[rho] = -(1/3) * Trace(Hessian[rho] * chi)`

The number of stencil points is specified by the user. For further information, see [the cited paper](http://dx.doi.org/10.1039/C4CP03106G).

## Parameters / inputs

- `npoints` — Three positive integer grid dimensions, ordered `[X Y Z]`. Each dimension must be at least `nstenc`.
- `extents` — Three positive axis extents in Angstroms, ordered `[X Y Z]`.
- `chi` — Real symmetric 3×3 electron magnetic susceptibility tensor in cubic Angstroms.
- `nstenc` — Odd integer number of finite-difference stencil points, at least 3. Periodic boundary conditions are used.

## Outputs

- `K` — Sparse matrix acting on the vectorisation of `rho`, whose dimensions are ordered `[X Y Z]`.

## Numerical / algorithmic content

The implementation constructs finite-difference second-derivative operators for all Cartesian axis pairs. Each operator is scaled by the grid-point counts divided by the corresponding axis extents, then weighted by the matching component of `chi` and summed with the factor `-1/3`.

[Spinach documentation](https://spindynamics.org/wiki/index.php?title=fdkup.m)