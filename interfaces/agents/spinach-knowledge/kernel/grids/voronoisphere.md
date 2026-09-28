# kernel/grids/voronoisphere.m

- Signature: `[vertices,indices,polygons,sangles]=voronoisphere(xyz)`

## Purpose

Constructs the Voronoi tessellation on the unit sphere for the supplied vectors. It is the geometric dual of the Delaunay triangulation obtained from their convex hull.

## Mathematical content

Each Voronoi vertex is the circumcentre of a Delaunay triangle, placed on the side of the triangle plane that the convex hull guarantees to be empty. Each cell is geodesically convex and star-shaped around its generator, so sorting its vertices by angle gives their counterclockwise order when viewed from outside.

## Parameters

- `xyz` — a `(3 x n)` array of vectors in `R^3`; the vectors are normalised. At least four finite, real, nonzero vectors are required. Their directions must be distinct, must not all lie on a single great circle, and must form a closed convex-hull surface.

## Outputs

- `vertices` — a `(3 x m)` array of Voronoi vertex coordinates.
- `indices` — an `(n x 1)` cell array; element `j` contains the indices of the vertices of the Voronoi cell corresponding to `xyz(:,j)`, in counterclockwise order when viewed from outside.
- `polygons` — an `(n x 1)` cell array; element `j` contains the coordinates of the vertices of the `j`th Voronoi cell in the same counterclockwise order.
- `sangles` — an `(n x 1)` array of Voronoi-cell solid angles, returned when requested.
