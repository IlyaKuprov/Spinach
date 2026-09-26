# kernel/grids/get_hull.m

- Signature: `[hull,edges]=get_hull(theta_angles,phi_angles)`

## Purpose

Generates a convex hull of a two-angle grid for 2D surface plotting. Syntax: [hull,edges]=get_hull(theta_angles,phi_angles)

## Physical / mathematical content

- Quadrature and geometry utilities. These files generate spherical/SO(3) grids, Voronoi weights, and adaptive integration tools for orientation averaging.

## Numerical / algorithmic content

## Parameters / inputs

- theta_angles -a column vector of theta angles,
- polar coordinates, ISO convention,
- radians
- phi_angles -a column vector of phi angles,
- polar coordinates, ISO convention,
- radians

## Outputs

- hull -a matrix of point indices of
- dimension Nx3, where N is the
- number of triangular facets
- edges -a matrix of point indices of
- dimension Nx2, where N is the
- number of grid edges

## Implementation structure

- Generates a convex hull of a two-angle grid for 2D
- surface plotting. Syntax:
- [hull,edges]=get_hull(theta_angles,phi_angles)
- theta_angles -a column vector of theta angles,
- polar coordinates, ISO convention,
- radians
- phi_angles -a column vector of phi angles,
- hull -a matrix of point indices of
- dimension Nx3, where N is the
- number of triangular facets
- edges -a matrix of point indices of
- dimension Nx2, where N is the
