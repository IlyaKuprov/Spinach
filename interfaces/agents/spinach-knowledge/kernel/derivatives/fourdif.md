# kernel/derivatives/fourdif.m

- Signature: `[x,DM]=fourdif(N,m)`

## Purpose

Computes the m-th derivative Fourier spectral differentiation matrix on the equispaced grid `x_j = 2*pi*j/N`, for `j=0,...,N-1`, in `[0,2*pi)`.

## Inputs

- `N` - positive real integer grid dimension.
- `m` - positive real integer derivative order.

## Outputs

- `x` - the grid points.
- `DM` - the m-th order differentiation matrix, constructed when two outputs are requested.

## Algorithm

For `m=1` and `m=2`, explicit formulae compute the first column; the implementation uses the flipping trick to improve accuracy. For `m>2`, it uses a discrete Fourier approach. The first row and column are assembled into `DM` with `toeplitz()`.

The code contains an `m=0` identity-matrix branch, but the input check requires `m` to be a positive integer, so validated calls cannot reach it.

## References

- S.C. Reddy and J.A.C. Weideman, [doi:10.1137/0916073](http://dx.doi.org/10.1137/0916073).
- [Spinach Wiki documentation](https://spindynamics.org/wiki/index.php?title=fourdif.m)
