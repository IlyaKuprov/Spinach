# kernel/derivatives/fourlap.m

- Signature: `L=fourlap(npoints,extents)`

## Purpose

Constructs the Fourier spectral Laplacian for data on a periodic grid. The documentation describes three axes ordered `[X Y Z]`; the implementation also has one- and two-dimensional branches.

## Inputs

- `npoints` - documented as positive integer grid-point counts ordered by axis; the implementation has branches for one, two, or three entries.
- `extents` - documented as positive real axis extents ordered `[X Y Z]`. The code checks positivity and reality but does not check that the input lengths match.

## Output

- `L` - Laplacian matrix acting on the vectorized data array. Each axis uses a second-derivative matrix from `fourdif()` scaled by `(2*pi/extents(i))^2`; the multidimensional cases combine these with Kronecker sums.

The documented boundary conditions are periodic.

## Documentation

- [Spinach Wiki documentation](https://spindynamics.org/wiki/index.php?title=fourlap.m)
