# kernel/grids/get_hull.m

- Signature: `[hull,edges]=get_hull(theta_angles,phi_angles)`

## Purpose

Generates the convex hull of a two-angle grid for 2D surface plotting.

## Parameters / inputs

- `theta_angles`: column vector of theta angles in radians, using polar coordinates in the ISO convention.
- `phi_angles`: column vector of phi angles in radians, using polar coordinates in the ISO convention.

Both inputs must be real, finite numeric column vectors with the same number of elements.

## Outputs

- `hull`: N×3 matrix of point indices, where N is the number of triangular facets.
- `edges`: N×2 matrix of point indices, where N is the number of grid edges.

## Implementation

Converts the angles to Cartesian coordinates, computes the convex hull with `convhull(x,y,z)`, then extracts the facet edges. The returned edge list includes both directions of each edge and excludes self-edges.

[Source reference](https://spindynamics.org/wiki/index.php?title=get_hull.m)