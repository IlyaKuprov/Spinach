# kernel/utilities/xyz2pd.m

- Signature: `density=xyz2pd(coords,x_range,y_range,z_range,...`

## Purpose

Probability density estimation for a three-dimensional Cartesian point cloud on a user-specified regular grid. Syntax: density=xyz2pd(coords,x_range,y_range,z_range,... x_npts, y_npts, z_npts)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

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
