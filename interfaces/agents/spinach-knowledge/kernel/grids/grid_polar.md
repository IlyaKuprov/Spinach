# kernel/grids/grid_polar.m

- Signature: `[phi,r,L]=grid_polar(ncircles,rmax)`

## Behaviour

Builds a balanced two-dimensional polar point set with a centre point and `ncircles-1` nonzero rings. The radii are equally spaced from zero to `rmax`; ring populations rise outward with circumference, so point density does not increase toward the centre. On ring `k` (counted outward from one), the source samples angles from zero through `2*pi` at `6*k` linspace points and removes the repeated endpoint, leaving `6*k-1` points. Thus the output vectors have `N=1+sum(6*k-1,k=1..ncircles-1)=3*(ncircles-1)^2+2*(ncircles-1)+1` entries, in matching order; `phi` is in radians and `r` has the same length unit as `rmax`.

`ncircles` must be a real integer at least two; `rmax` must be a positive real scalar. The centre has `r=0` and `phi=0`.

## Laplacian

When the third output is requested, the points are converted to Cartesian coordinates and Delaunay-triangulated. Each triangulation edge receives weight `1/d^2`, where `d` is its endpoint distance. The symmetric graph Laplacian uses negative off-diagonal edge weights and positive diagonal row sums, then is divided by its matrix 2-norm and converted to sparse form. `L` is `N` by `N`, follows the point ordering above, and is dimensionless after this normalisation.

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/grids/grid_polar.m) · [Wiki](https://spindynamics.org/wiki/index.php?title=grid_polar.m)