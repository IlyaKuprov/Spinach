# kernel/derivatives/fdmat.m

Direct source: [kernel/derivatives/fdmat.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/derivatives/fdmat.m)
Spin Dynamics Wiki: [fdmat.m](https://spindynamics.org/wiki/index.php?title=fdmat.m)

## Purpose and interface

D=fdmat(dim,nstenc,order,boundary) returns a sparse dim-by-dim matrix that applies a finite-difference derivative to a vector of dim samples. The stencil is built on unit-spaced integer sample positions. The optional boundary argument defaults to 'pbc'; the other supported value is 'wall'.

- dim is a positive integer of at least 3.
- nstenc is a positive odd integer, the number of grid points in the stencil.
- order is a positive integer smaller than nstenc.
- boundary must be a character string equal to 'pbc' or 'wall'; any other value errors. The source does not separately check that nstenc is no greater than dim.
- D is sparse and square. Applying D*x maps the sampled vector to derivative estimates at the same grid positions.

## How the matrix is assembled

The routine calls fdweights to obtain the coefficient row for the requested derivative, then inserts those coefficients into a sparse matrix preallocated for up to dim*nstenc entries. For interior points, the sample offsets are the centred integers from -(nstenc-1)/2 through (nstenc-1)/2, and the resulting row is placed at each interior grid point.

With 'wall', the first (nstenc-1)/2 rows use one-sided stencils drawn from the first nstenc samples. The matching rows at the far edge use the reversed coefficient order and the factor (-1)^order. The remaining rows use centred coefficients.

With 'pbc', every row uses the centred coefficient row. The column indices are wrapped into 1:dim by modulo indexing, which implements periodic indexing.

No grid-spacing argument is applied: coefficients correspond to unit sample spacing. If the sample spacing represents a physical interval other than one, the derivative must be scaled for that coordinate convention by the caller.

## Example

    D=fdmat(5,3,1,'wall');

This makes a 5-by-5 first-derivative matrix. The first row uses the three-point forward weights [-3/2, 2, -1/2], interior rows use the centred weights [-1/2, 0, 1/2], and the last row uses the corresponding backward weights. Multiplication by a sample vector returns one estimate at each of the five positions.

## Related routine

- [fdweights.m](fdweights.md) generates the coefficient rows used here.
