# kernel/derivatives/fdlap.m

Source: [kernel/derivatives/fdlap.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/derivatives/fdlap.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=fdlap.m)

- Signature: `L=fdlap(dims,extents,nstenc)`

## Purpose and inputs

Constructs a sparse finite-difference Laplacian for a vectorised one-, two-, or three-dimensional array with axes ordered `[X Y Z]` (using only the leading axes for lower-dimensional inputs). Boundary conditions are periodic.

- `dims`: vector of one, two, or three positive integer grid counts.
- `extents`: corresponding vector of positive real axis extents.
- `nstenc`: odd integer stencil-point count of at least 3; every grid dimension must be at least this large.

## Output and assembly

`L` is a sparse square matrix with one row and column per grid point, acting on the column-major vectorisation of the array. The source obtains a second-derivative `fdmat` matrix for each axis and scales the axis-`i` term by `(dims(i)/extents(i))^2`. In two dimensions it adds the Y and X terms as Kronecker products; in three dimensions it adds Z, Y, and X terms with identity factors on the other axes. This ordering makes X (the first array dimension) the fastest-varying factor. For a one-dimensional input, the result is just the scaled X derivative matrix.

## Guards and naming clarification

The implementation accepts one, two, or three dimensions and rejects other dimension counts. It checks positive real integer grid counts, positive real extents, sufficient points for the stencil, and an odd integer stencil size of at least 3. The source comment's syntax line names the first argument `npoints`, while the MATLAB function signature and implementation call it `dims`.

For the related 3D Hessian and tensor-weighted derivative constructors, see [fdhess.m](fdhess.md) and [fdkup.m](fdkup.md).
