# kernel/derivatives/fdkup.m

Source: [kernel/derivatives/fdkup.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/derivatives/fdkup.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=fdkup.m)

- Signature: `K=fdkup(npoints,extents,chi,nstenc)`

## Purpose

Builds a sparse finite-difference matrix for the Kuprov operator acting on the vectorisation of a 3D array `rho`, with array axes ordered `[X Y Z]`. The source states its action as `K[rho] = -(1/3) * Trace(Hessian[rho] * chi)` and cites [the original paper](http://dx.doi.org/10.1039/C4CP03106G).

## Inputs and output

- `npoints`: three positive integer grid counts in `[X Y Z]` order.
- `extents`: three positive real axis extents in the same order; the source documentation gives these in Angstroms.
- `chi`: real numeric symmetric 3-by-3 electron magnetic susceptibility tensor; the source documents cubic Angstrom units.
- `nstenc`: odd integer stencil size of at least 3; each grid count must be at least this large. The finite-difference scheme uses periodic boundary conditions.
- `K`: sparse square matrix acting on `rho(:)`, with one row and column per grid point (matrix dimension `prod(npoints)`). No `rho` values are passed to this constructor.

## Tensor weighting and grid scaling

The code builds second-derivative matrices with Kronecker products in the `[X Y Z]` array order. Diagonal derivatives use a second-derivative `fdmat` on one axis; mixed derivatives use first-derivative matrices on both axes, with identities on untouched axes. It aliases the reversed mixed-derivative matrices to the corresponding forward pair. Each derivative factor along axis `i` is scaled by `npoints(i)/extents(i)`; a diagonal derivative therefore receives the squared axis factor, while a mixed derivative receives the product of its two axis factors.

The final matrix is the explicit component sum `-(1/3) * sum(chi(i,j) * D_ij)` over all nine ordered axis pairs, where `D_ij` is the corresponding discretised second derivative. Thus `chi` is consumed directly in the supplied Cartesian component order; this routine has no separate tensor-frame rotation or coordinate-conversion input. The source comment describes `npoints` as dimensions “in Angstroms,” while the implementation uses these values as integer grid counts and applies physical extents through `extents`.

## Guards

The code checks that `npoints` has three real positive integer entries, `extents` has three positive real entries, each point count meets the stencil size, and `chi` is numeric, real, 3-by-3, and symmetric to an absolute 1-norm tolerance of `1e-10`. It requires an odd integer stencil count of at least 3. Related constructors: [fdhess.m](fdhess.md) and [fdlap.m](fdlap.md).
