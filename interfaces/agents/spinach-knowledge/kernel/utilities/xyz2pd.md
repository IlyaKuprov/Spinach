# kernel/utilities/xyz2pd.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/xyz2pd.m`
- Signature: `density=xyz2pd(coords,x_range,y_range,z_range,...`
- Total lines: 110

## Purpose

Probability density estimation for a three-dimensional Cartesian point cloud on a user-specified regular grid. Syntax: density=xyz2pd(coords,x_range,y_range,z_range,... x_npts, y_npts, z_npts)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- coords -an N-by-3 array of Cartesian coordinates
- x_range -a two-element vector [xmin xmax] specifying the
- Cartesian grid extent along the x axis
- y_range -a two-element vector [ymin ymax] specifying the
- Cartesian grid extent along the y axis
- z_range -a two-element vector [zmin zmax] specifying the
- Cartesian grid extent along the z axis
- x_npts -the number of grid points along the x axis
- y_npts -the number of grid points along the y axis
- z_npts -the number of grid points along the z axis

## Outputs

- density -a three-dimensional array containing the num-
- ber of points falling into each grid cell

## Implementation structure

- Probability density estimation for a three-dimensional Cartesian
- point cloud on a user-specified regular grid. Syntax:
- density=xyz2pd(coords,x_range,y_range,z_range,...
- x_npts, y_npts, z_npts)
- coords -an N-by-3 array of Cartesian coordinates
- x_range -a two-element vector [xmin xmax] specifying the
- Cartesian grid extent along the x axis
- y_range -a two-element vector [ymin ymax] specifying the
- Cartesian grid extent along the y axis
- z_range -a two-element vector [zmin zmax] specifying the
- Cartesian grid extent along the z axis
- x_npts -the number of grid points along the x axis

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `x_range()`, `y_range()`, `z_range()`, `discretize()`, `coords()`, `isnan()`, `cell_idx_x()`, `cell_idx_y()`, `cell_idx_z()`, `sub2ind()`, `accumarray()`, `any()`, `isscalar()`.
