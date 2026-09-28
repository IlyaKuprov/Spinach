# kernel/pulses/wave_basis.m

- Signature: `basis_waves=wave_basis(basis_type,n_func,n_points)`

## Purpose

Returns sampled sine, cosine, or Legendre basis functions for pulse-waveform expansion as columns of a matrix.

## Numerical / algorithmic content

The sine and cosine functions are sampled over `[-pi,pi]`; Legendre polynomials are sampled over `[-1,1]`. Rows are generated for the requested functions, then `orth(basis_waves')` orthogonalizes the sampled vectors and returns them as columns. Discretization means the raw functions may not be precisely orthogonal under the standard scalar product; orthogonalization can flip the sign of some functions. If the sampled basis is rank-deficient, the function errors and asks to reduce `n_func`.

## Parameters / inputs

- `basis_type` - character string: `'sine_waves'`, `'cosine_waves'`, or `'legendre'`
- `n_func` - positive integer number of functions; sine frequencies start at 1, cosine frequencies at 0, and Legendre polynomial ranks at 0
- `n_points` - positive integer number of discretization points

## Outputs

- `basis_waves` - matrix with the orthogonalized basis functions in columns

Source Wiki page: https://spindynamics.org/wiki/index.php?title=wave_basis.m
