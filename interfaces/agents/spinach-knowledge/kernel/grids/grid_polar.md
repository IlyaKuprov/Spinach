# kernel/grids/grid_polar.m

- Signature: `[phi,r,L]=grid_polar(ncircles,rmax)`

## Purpose

Generates a balanced polar grid whose point density does not increase towards the centre.

## Parameters / inputs

- `ncircles`: Number of radial circles; an integer greater than one.
- `rmax`: Maximum grid radius; a positive real number.

## Outputs

- `phi`: Column vector of polar angles in radians.
- `r`: Column vector of radii.
- `L`: Sparse Laplacian operator acting on values ordered as the grid; constructed when requested.

## Implementation

The radii are evenly spaced from zero to `rmax`. Each successive circle is assigned more angular grid points. When `L` is requested, the function triangulates the Cartesian coordinates, weights connected vertices by inverse squared distance, and normalizes the resulting Laplacian.

[Source documentation](https://spindynamics.org/wiki/index.php?title=grid_polar.m)