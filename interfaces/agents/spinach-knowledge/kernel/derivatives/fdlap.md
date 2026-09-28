# kernel/derivatives/fdlap.m

- Signature: `L=fdlap(dims,extents,nstenc)`

## Purpose

Constructs a sparse finite-difference Laplacian for a vectorized 1D, 2D, or 3D array whose dimensions are ordered as `[X Y Z]`. The finite-difference approximation uses periodic boundary conditions.

## Parameters / inputs

- `dims`: One-, two-, or three-element vector of positive integers giving the number of discretization points along each dimension, ordered as `[X Y Z]`.
- `extents`: Corresponding one-, two-, or three-element vector of positive real sizes, ordered as `[X Y Z]`.
- `nstenc`: Number of finite-difference stencil points; it must be an odd integer of at least 3, and every value in `dims` must be at least `nstenc`.

## Output

- `L`: Sparse Laplacian matrix acting on the vectorization of the array.

## Construction

For each dimension, `fdlap` obtains a second-derivative matrix using `fdmat(dims(i),nstenc,2)` and scales it by `(dims(i)/extents(i))^2`. In 1D, `L` is the scaled matrix `Dxx`. In 2D and 3D, `L` is the sum of the scaled second-derivative matrices expanded along the other dimensions with Kronecker products and sparse identity matrices (`speye`). The X dimension is the innermost factor, followed by Y and then Z.

The function checks that `dims` contains positive integers, `extents` contains positive real values, and the stencil satisfies the size and parity requirements. It rejects numbers of spatial dimensions other than one, two, or three.

## Reference

- [Spinach `fdlap.m` documentation](https://spindynamics.org/wiki/index.php?title=fdlap.m).