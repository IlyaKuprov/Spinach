# kernel/derivatives/fourlap.m

[Direct MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/derivatives/fourlap.m) · [Spinach Wiki documentation](https://spindynamics.org/wiki/index.php?title=fourlap.m)

## Purpose and signature

`L=fourlap(npoints,extents)` constructs a Fourier spectral Laplacian for a periodic rectangular grid. Axis order is `[X Y Z]`; the implementation has one-, two-, and three-axis branches. This is a matrix operator, not a routine that accepts or reshapes the data array itself.

## Inputs and output

- `npoints`: positive integer sample counts, one per axis, ordered `[X Y Z]` when three axes are supplied.
- `extents`: positive real periodic lengths for the corresponding axes, in the same order.
- `L`: square operator of size `prod(npoints)` by `prod(npoints)`, acting on the column-major vectorisation of the grid data. In one dimension the source assigns the second-derivative matrix directly; in two and three dimensions it builds Kronecker sums with sparse identity factors.

Each axis second-derivative matrix comes from `fourdif(spin_system,npoints(i),2)` and is scaled by `(2*pi/extents(i))^2`. The multidimensional Kronecker ordering makes X the fastest-varying index in the vectorised data. The resulting operator has inverse-length-squared scaling when `extents` are physical lengths; no physical unit is imposed by the function itself.

The boundary condition is periodic. The source validates positive real integer grid counts and positive real extents, selects branches by the number of entries in `npoints`, and errors for any other dimensionality. For the three-axis case the axis matrices are embedded and summed as the X, Y, and Z second derivatives. See also [`fourdif`](fourdif.md).
