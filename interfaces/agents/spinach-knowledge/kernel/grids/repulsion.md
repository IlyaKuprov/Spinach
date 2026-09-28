# kernel/grids/repulsion.m

- Signature: `[alphas,betas,gammas,weights]=repulsion(npoints,ndims,niter)`

Generates repulsion grids on a unit hypersphere. For information on the algorithm, see Bak and Nielsen (http://dx.doi.org/10.1006/jmre.1996.1087).

## Inputs

- `npoints`: positive integer number of grid points.
- `ndims`: `2`, `3`, or `4`, producing a single-angle, two-angle, or three-angle grid, respectively.
- `niter`: positive integer number of repulsion iterations.

## Algorithm

The function starts with random points in `ndims` dimensions. At each iteration it computes pairwise normalized distance vectors, removes self-interactions, forms forces using point scalar products, moves the points by `ndims*F/npoints`, and reprojects them onto the unit sphere. It reports the maximum point displacement each iteration.

## Outputs

- `alphas`, `betas`, `gammas`: angular coordinates in radians. For `ndims=2`, `betas` contains polar angles and the other outputs are zero. For `ndims=3`, `betas` is `theta+pi/2`, `gammas` is `phi`, and `alphas` is zero. For `ndims=4`, the coordinates are obtained from `qter2euler`.
- `weights`: uniform point weights, each `1/npoints`. The source suggests using SHREWD to generate optimal weights.

With `ndims=3` and no requested outputs, the function also plots the grid.