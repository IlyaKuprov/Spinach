# kernel/grids/get_hull.m

- Signature: `[hull,edges]=get_hull(theta_angles,phi_angles)`

## Purpose

Returns the triangulated convex-hull surface and its edge index pairs for a spherical angle grid, for example for 2D surface plotting.

## Mapping and outputs

Each paired angle is converted to a unit-sphere point using `x=sin(theta)*cos(phi)`, `y=sin(theta)*sin(phi)`, and `z=cos(theta)`; MATLAB `convhull(x,y,z)` supplies the triangular facet indices. `hull` is an `M-by-3` matrix, one row per returned triangular facet, with indices into the input angle vectors.

The routine forms edge pairs from the three sides of every facet, removes duplicate pairs, adds the reversed orientation of each pair, and removes self-edges. Consequently, each distinct undirected hull edge is represented in both directions in `edges` (an `N-by-2` index matrix). This is topology only: the function does not compute surface areas or quadrature weights.

## Parameters / inputs

- `theta_angles` — finite real column vector of polar angles in radians, using the ISO convention.
- `phi_angles` — finite real column vector of azimuthal angles in radians, using the ISO convention.

The vectors must have the same number of elements. The checker imposes no angle-range restriction; the trigonometric mapping places each pair on the unit sphere.

## Outputs

- `hull` — triangular facet indices, with three input-point indices per row.
- `edges` — two-column edge endpoint indices, including both orientations for every distinct hull edge.

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/grids/get_hull.m)
[Source reference](https://spindynamics.org/wiki/index.php?title=get_hull.m)
