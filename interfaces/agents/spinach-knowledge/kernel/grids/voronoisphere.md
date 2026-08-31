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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 41-42: Check consistency; implemented by `grumble(xyz)`.
- Lines 44-45: Normalise the points; implemented by `xyz=xyz./sqrt(sum(xyz.^2,1)); npts=size(xyz,2)`.
- Lines 47-48: Delaunay triangles of a spherical point set are convex hull facets; implemented by `triangles=convhull(xyz.'); ntrian=size(triangles,1)`.
- Lines 50-51: Refuse point sets whose hull is not a closed surface; implemented by `edges=sort([triangles(:,[1 2]); triangles(:,[2 3]); triangles(:,[3 1])],2)`.
- Lines 57-59: Voronoi vertices are triangle circumcentres on the empty side; implemented by `normals=cross(xyz(:,triangles(:,2))-xyz(:,triangles(:,1)), xyz(:,triangles(:,3))-xyz(:,triangles(:,1)),1)`.
- Lines 64-66: List the triangles incident at each point; implemented by `incidence=accumarray([triangles(:,1); triangles(:,2); triangles(:,3)], repmat((1:ntrian).',3,1),[npts 1],@(x){x})`.
- Lines 68-69: Build each cell by angular sorting around its generator; implemented by `indices=cell(npts,1); polygons=cell(npts,1)`.
- Lines 72-73: Right-handed tangent frame at the generator; implemented by `pivot=zeros(3,1); [~,pos]=min(abs(xyz(:,n))); pivot(pos)=1`.
- Lines 77-78: Sort the cell vertices counterclockwise from outside; implemented by `cell_verts=vertices(:,incidence{n})`.
- Lines 85-86: Return solid angles if needed; implemented by `if nargout>=4`.

### Control flow inferred from the code

- Line 53: conditional branch on `~all(accumarray(edge_ids,1)==2)`.
- Line 70: `for` loop over `n=1:npts`.
- Line 86: conditional branch on `nargout>=4`.

### Key state/data transformations

- Lines 45: computes `xyz` using `xyz=xyz./sqrt(sum(xyz.^2,1)); npts=size(xyz,2)`.
- Lines 48: computes `triangles` using `triangles=convhull(xyz.'); ntrian=size(triangles,1)`.
- Lines 51: computes `edges` using `edges=sort([triangles(:,[1 2]); triangles(:,[2 3]); triangles(:,[3 1])],2)`.
- Lines 52: computes `[~,~,edge_ids]` using `[~,~,edge_ids]=unique(edges,'rows')`.
- Lines 58-59: computes `normals` using `normals=cross(xyz(:,triangles(:,2))-xyz(:,triangles(:,1)), xyz(:,triangles(:,3))-xyz(:,triangles(:,1)),1)`.
- Lines 60: computes `vertices` using `vertices=normals./sqrt(sum(normals.^2,1))`.
- Lines 61: computes `empty_side` using `empty_side=sign(sum(vertices.*(xyz(:,triangles(:,1))-mean(xyz,2)),1))`.
- Lines 65-66: computes `incidence` using `incidence=accumarray([triangles(:,1); triangles(:,2); triangles(:,3)], repmat((1:ntrian).',3,1),[npts 1],@(x){x})`.
- Lines 69: computes `indices` using `indices=cell(npts,1); polygons=cell(npts,1)`.
- Lines 73: computes `pivot` using `pivot=zeros(3,1); [~,pos]=min(abs(xyz(:,n))); pivot(pos)=1`.
- Lines 74: computes `tang_a` using `tang_a=cross(xyz(:,n),pivot); tang_a=tang_a/norm(tang_a,2)`.
- Lines 75: computes `tang_b` using `tang_b=cross(xyz(:,n),tang_a)`.
- Lines 78: computes `cell_verts` using `cell_verts=vertices(:,incidence{n})`.
- Lines 79: computes `[~,ccw_order]` using `[~,ccw_order]=sort(atan2(tang_b.'*cell_verts,tang_a.'*cell_verts))`.
- Lines 80: computes `indices{n}` using `indices{n}=incidence{n}(ccw_order)`.
- Lines 81: computes `polygons{n}` using `polygons{n}=vertices(:,indices{n})`.
- Lines 87: computes `sangles` using `sangles=vcell_solidangle(vertices,indices,xyz)`.

### Local helper functions

- Line 93: `grumble()` — `function grumble(xyz)`.
  - Representative operation: `if (~isnumeric(xyz))||(~isreal(xyz))||(size(xyz,1)~=3)`.
  - Representative operation: `error('xyz must be a [3 x N] array of real numbers.')`.

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
