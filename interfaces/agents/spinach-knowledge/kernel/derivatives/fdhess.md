# kernel/derivatives/fdhess.m

Source: [kernel/derivatives/fdhess.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/derivatives/fdhess.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=fdhess.m)

- Signature: `H=fdhess(A,nstenc)`

## Purpose and inputs

Computes the finite-difference Hessian of a numeric three-dimensional array `A`, with dimensions ordered `[X Y Z]`. The source specifies unit grid spacing, so the derivatives are per unit grid coordinate; this routine has no physical extents or unit-conversion input. The stencil uses periodic boundary conditions.

- `A`: numeric 3D array.
- `nstenc`: odd integer stencil-point count, at least 3. Each dimension of `A` must be at least this large.

## Output and assembly

`H` is a 3-by-3 cell array; every cell contains a 3D array with the same shape as `A`. Its rows and columns are the derivative axes `[X Y Z]`, in this order:

`{d2A_dxdx  d2A_dxdy  d2A_dxdz; d2A_dydx  d2A_dydy  d2A_dydz; d2A_dzdx  d2A_dzdy  d2A_dzdz}`

For each entry, the code applies second-derivative `fdmat` operators on a diagonal axis, or first-derivative operators on both axes for a mixed derivative; identity factors occupy untouched axes. Kronecker products act on `A(:)`, and each result is reshaped to `size(A)`. With MATLAB column-major vectorisation, the first (X) dimension is the fastest-varying factor. The source computes all nine ordered entries and places them in the displayed cell-array order.

## Guards

The source requires `A` to be numeric and three-dimensional, every dimension to meet the stencil size, and `nstenc` to be an odd integer of at least 3. The guard accepts 3 as the minimum stencil count. See related operators [fdkup.m](fdkup.md) and [fdlap.m](fdlap.md).
