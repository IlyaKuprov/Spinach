# kernel/grids/grid_polar.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/grids/grid_polar.m`
- Signature: `[phi,r,L]=grid_polar(ncircles,rmax)`
- Total lines: 98

## Purpose

Generates a balanced polar grid in which the density of points does not increase towards the centre. Syntax: [phi,r,L]=grid_polar(ncircles,rmax)

## Physical / mathematical content

- Quadrature and geometry utilities. These files generate spherical/SO(3) grids, Voronoi weights, and adaptive integration tools for orientation averaging.

## Numerical / algorithmic content

- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- ncircles -number of radial circles
- in the grid, an integer
- rmax -maximum radius that the
- grid must reach

## Outputs

- phi -a column vector of polar
- phi angles, radians
- rmax -a coluim vector of radii
- L -sparse Laplacian operator
- acting on functions defined
- as vector of values in the
- same order as the grid

## Implementation structure

- Generates a balanced polar grid in which the density of
- points does not increase towards the centre. Syntax:
- [phi,r,L]=grid_polar(ncircles,rmax)
- ncircles -number of radial circles
- in the grid, an integer
- rmax -maximum radius that the
- grid must reach
- phi -a column vector of polar
- phi angles, radians
- rmax -a coluim vector of radii
- L -sparse Laplacian operator
- acting on functions defined

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `phi_current()`, `radii()`, `r_current()`, `phi()`, `delaunay()`, `tri()`, `transpose()`, `isscalar()`.
