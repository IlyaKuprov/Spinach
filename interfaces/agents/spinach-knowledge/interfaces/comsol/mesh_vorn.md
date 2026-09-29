# interfaces/comsol/mesh_vorn.m

[Canonical MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/interfaces/comsol/mesh_vorn.m) · [Spinach Wiki: mesh_vorn.m](https://spindynamics.org/wiki/index.php?title=mesh_vorn.m)

## Purpose and return value

`mesh=mesh_vorn(mesh)` computes a two-dimensional Voronoi tessellation from the mesh coordinates, keeps the cells belonging to active mesh vertices, adds area weights, and returns the updated mesh structure. This adapter operates on an existing MATLAB structure; it does not call COMSOL.

## Accepted data and dependencies

The input must provide coordinate columns `mesh.x` and `mesh.y` and an active-vertex index vector `mesh.idx.active`. The source checks for `mesh.idx.active`, then passes `[mesh.x mesh.y]` to MATLAB `voronoin`; compatible two-dimensional coordinates and valid active indices are therefore required. MATLAB `polyarea` is used for areas. The code does not rescale coordinates: the entries of `mesh.vor.weights` have the square of the coordinate unit (for example, if coordinates are in millimetres, the weights are in square millimetres).

## Tessellation transformation and guardrails

`mesh.vor.vertices` and `mesh.vor.cells` receive the `voronoin` result, after which cells are filtered by `mesh.idx.active`. A cell containing Voronoi vertex index `1` is treated as unbounded and causes an error naming the corresponding active mesh vertices; inactivate such vertices before this adapter is used. For each remaining cell, the polygon area is stored in `mesh.vor.weights`. The function also sets `mesh.vor.ncells` to the number of retained cells and `mesh.vor.max_cell_size` to the largest number of vertices in a retained cell. There is no separate status output: success is represented by the updated `mesh`.
