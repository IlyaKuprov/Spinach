# kernel/derivatives/fdmat.m

- Signature: `D=fdmat(dim,nstenc,order,boundary)`

## Purpose

Returns a sparse, arbitrary-order finite-difference differentiation matrix for unit grid point spacing.

## Parameters / inputs

- `dim` — dimension of the column vector to be differentiated; must be an integer of at least 3.
- `nstenc` — number of points in the finite-difference stencil; must be a positive odd integer.
- `order` — derivative order; must be a positive integer smaller than `nstenc`.
- `boundary` — `'wall'` uses sided finite-difference schemes at the edges; `'pbc'` assumes periodic boundaries. Defaults to `'pbc'` and must be a character string.

## Outputs

- `D` — sparse finite-difference differentiation matrix of size `dim` by `dim`.

## Numerical / algorithmic content

The matrix is preallocated with space for `dim*nstenc` entries. Finite-difference coefficients are obtained from `fdweights`.

- For `'wall'`, the first `(nstenc-1)/2` rows use sided stencils spanning the first `nstenc` grid points. The corresponding rows at the opposite edge use the reversed coefficients, multiplied by `(-1)^order`. Interior rows use a centered stencil.
- For `'pbc'`, every row uses the same centered-stencil coefficients. Column indices wrap around the matrix using modulo indexing.

An unrecognized boundary type raises an error.

## Link

- <https://spindynamics.org/wiki/index.php?title=fdmat.m>
