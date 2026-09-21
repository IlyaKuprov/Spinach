# kernel/grids/voronoisphere.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/grids/voronoisphere.m`
- Signature: `[vertices,indices,polygons,sangles]=voronoisphere(xyz)`
- Total lines: 118

## Purpose

Voronoi tessellation of the unit sphere around the specified po- ints, computed as the exact geometric dual of the Delaunay tri- angulation returned by the convex hull. Syntax: [vertices,indices,polygons,sangles]=voronoisphere(xyz)

## Physical / mathematical content

- Quadrature and geometry utilities. These files generate spherical/SO(3) grids, Voronoi weights, and adaptive integration tools for orientation averaging.

## Numerical / algorithmic content

- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- xyz -(3 x n) array, coordinates of n distinct vectors
- in R^3; these will be normalised

## Outputs

- vertices -(3 x m) array, coordinates of the vertices of the
- Voronoi tessellation
- indices -(n x 1) cell array, j-th element contains the in-
- dices of the Voronoi cell vertices that correspond
- to xyz(:,j). Vertices are oriented counterclockwise
- when looking from outside.
- polygons -(n x 1) cell array, j-th element contains the coor-
- dinates of the vertices of the j-th Voronoi cell,
- in the same counterclockwise order
- sangles -(n x 1) array, solid angles of each Voronoi cell
- Note: every Voronoi vertex is the circumcentre of a Delaunay tri-
- angle, placed on the side of the triangle plane that the
- convex hull property guarantees to be empty; every Voronoi
- cell is geodesically convex and therefore star-shaped around
- its generator point, which makes the angular sort used here
- an exact construction rather than a heuristic.

## Implementation structure

- Voronoi tessellation of the unit sphere around the specified po-
- ints, computed as the exact geometric dual of the Delaunay tri-
- angulation returned by the convex hull. Syntax:
- [vertices,indices,polygons,sangles]=voronoisphere(xyz)
- xyz -(3 x n) array, coordinates of n distinct vectors
- in R^3; these will be normalised
- vertices -(3 x m) array, coordinates of the vertices of the
- Voronoi tessellation
- indices -(n x 1) cell array, j-th element contains the in-
- dices of the Voronoi cell vertices that correspond
- to xyz(:,j). Vertices are oriented counterclockwise
- when looking from outside.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `convhull()`, `triangles()`, `all()`, `accumarray()`, `cross()`, `xyz()`, `sign()`, `pivot()`, `vertices()`, `atan2()`, `vcell_solidangle()`, `any()`, `uniquetol()`.
