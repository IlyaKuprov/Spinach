# interfaces/comsol/mesh_crop.m

- Signature: `mesh=mesh_crop(mesh,ranges)`

## Purpose

Crops a 2D mesh to the rectangular coordinate window `[xmin,xmax] × [ymin,ymax]`. Vertices on the bounds are retained.

## Behavior

The routine removes cached Voronoi and plotting data when present, selects vertices within both coordinate ranges, keeps edges, triangles, and rectangles whose vertices all survive, and reindexes those elements. It crops coordinates and any present velocity or concentration arrays (`u`, `v`, and `c`). The active-vertex list is replaced by the vertices appearing in the retained triangles; an existing list triggers a warning.

## Parameters / inputs

- `mesh`: Spinach mesh object with vertex-index data.
- `ranges`: two-element cell array `{[xmin xmax],[ymin ymax]}`; each pair must contain two real, increasing bounds.

## Output

- `mesh`: cropped and reindexed mesh object.

## Source

[Spinach Wiki: mesh_crop.m](https://spindynamics.org/wiki/index.php?title=mesh_crop.m)
