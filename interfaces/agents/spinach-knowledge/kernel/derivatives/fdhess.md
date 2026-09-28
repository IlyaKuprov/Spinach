# kernel/derivatives/fdhess.m

- Signature: `H=fdhess(A,nstenc)`

## Purpose

Computes the finite-difference Hessian of a numeric 3D array on a unit-spaced grid, with dimensions ordered `[X Y Z]`.

## Parameters / inputs

- `A`: numeric 3D array; each dimension must contain at least `nstenc` elements.
- `nstenc`: odd integer number of stencil points, at least 3. Periodic boundary conditions are used.

## Output

`H` is a 3×3 cell array of 3D derivative arrays, ordered as

`{d2A_dxdx d2A_dxdy d2A_dxdz; d2A_dydx d2A_dydy d2A_dydz; d2A_dzdx d2A_dzdy d2A_dzdz}`

## Implementation

Diagonal entries use second-derivative `fdmat` operators. Mixed entries apply first-derivative `fdmat` operators along both corresponding dimensions and are assembled with Kronecker products.
