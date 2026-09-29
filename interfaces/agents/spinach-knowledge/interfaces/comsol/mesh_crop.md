# interfaces/comsol/mesh_crop.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/interfaces/comsol/mesh_crop.m) · [Spinach Wiki: mesh_crop.m](https://spindynamics.org/wiki/index.php?title=mesh_crop.m)

## Purpose and call

`mesh = mesh_crop(mesh,ranges)` crops a two-dimensional mesh to a rectangular window in its x-y coordinates and returns the updated mesh structure. Coordinate values and any retained per-vertex data remain in their input units; no unit conversion is applied.

## Accepted data and transformation

`ranges` must be a two-element cell array, `{[xmin xmax],[ymin ymax]}`. Each bound pair must be numeric, real, contain two elements, and have its first value strictly less than its second. Vertices on either bound are retained. The routine keeps an edge, triangle, or rectangle only when all its vertex indices refer to retained vertices, then remaps those indices to the cropped coordinate arrays.

The routine crops `mesh.x` and `mesh.y`, and also crops `mesh.u`, `mesh.v`, and `mesh.c` along their vertex dimension when those fields exist. It removes cached `mesh.vor` and `mesh.plot` fields when present because they refer to the previous mesh. If `mesh.idx.active` already exists, it warns that the list is being overwritten; the new list is the unique vertex indices appearing in the retained triangles.

## Output and guardrails

The returned `mesh` contains the cropped coordinates and reindexed connectivity. The routine explicitly checks for `mesh.idx` and checks the range container, bound types, sizes, and ordering. It does not explicitly require finite bounds; the coordinates and connectivity arrays used by the crop are otherwise assumed to be present and compatible.
