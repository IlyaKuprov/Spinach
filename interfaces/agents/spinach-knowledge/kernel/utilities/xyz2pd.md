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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 38-39: Check consistency; implemented by `grumble(coords,x_range,y_range,z_range,x_npts,y_npts,z_npts)`.
- Lines 41-42: Compute grid cell edges; implemented by `x_edges=linspace(x_range(1),x_range(2),x_npts+1)`.
- Lines 46-47: Assign each point to its grid cell; implemented by `cell_idx_x=discretize(coords(:,1),x_edges)`.
- Lines 51-54: Discard points outside the grid; implemented by `valid_mask=(~isnan(cell_idx_x))& (~isnan(cell_idx_y))& (~isnan(cell_idx_z))`.
- Lines 59-61: Assign a linear index to each cell; implemented by `linear_index=sub2ind([x_npts y_npts z_npts], cell_idx_x,cell_idx_y,cell_idx_z)`.
- Lines 63-64: Count the number of points falling into each cell; implemented by `density=accumarray(linear_index,1,[x_npts*y_npts*z_npts 1])`.
- Lines 66-67: Reshape to create a three-dimensional density; implemented by `density=reshape(density,[x_npts y_npts z_npts])`.

### Key state/data transformations

- Lines 42: computes `x_edges` using `x_edges=linspace(x_range(1),x_range(2),x_npts+1)`.
- Lines 43: computes `y_edges` using `y_edges=linspace(y_range(1),y_range(2),y_npts+1)`.
- Lines 44: computes `z_edges` using `z_edges=linspace(z_range(1),z_range(2),z_npts+1)`.
- Lines 47: computes `cell_idx_x` using `cell_idx_x=discretize(coords(:,1),x_edges)`.
- Lines 48: computes `cell_idx_y` using `cell_idx_y=discretize(coords(:,2),y_edges)`.
- Lines 49: computes `cell_idx_z` using `cell_idx_z=discretize(coords(:,3),z_edges)`.
- Lines 52-54: computes `valid_mask` using `valid_mask=(~isnan(cell_idx_x))& (~isnan(cell_idx_y))& (~isnan(cell_idx_z))`.
- Lines 60-61: computes `linear_index` using `linear_index=sub2ind([x_npts y_npts z_npts], cell_idx_x,cell_idx_y,cell_idx_z)`.
- Lines 64: computes `density` using `density=accumarray(linear_index,1,[x_npts*y_npts*z_npts 1])`.

### Local helper functions

- Line 72: `grumble()` — `function grumble(coords,x_range,y_range,z_range,`.
  - Representative operation: `x_npts, y_npts, z_npts)`.
  - Representative operation: `if (~isnumeric(coords))||(size(coords,2)~=3)|| (~isreal(coords))||any(~isfinite(coords),'all')`.

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
