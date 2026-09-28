# interfaces/comsol/mesh_vorn.m

- Signature: `mesh=mesh_vorn(mesh)`

## Purpose

Computes a Voronoi tessellation of a two-dimensional COMSOL mesh and stores cells for the active vertices.

## Behavior

The routine calls MATLAB `voronoin` on the mesh coordinates, retains cells indexed by `mesh.idx.active`, and stops with an error if any retained cell is unbounded. It computes each retained cell's polygon area and stores the areas as Voronoi weights, along with the tessellation vertices and cells, cell count, and maximum cell size.

## Parameters / inputs

- `mesh`: Spinach mesh object containing coordinates and an active-vertex index list.

## Output

- `mesh`: updated mesh object with Voronoi tessellation and area-weight data in `mesh.vor`.

## Source

[Spinach Wiki: mesh_vorn.m](https://spindynamics.org/wiki/index.php?title=mesh_vorn.m)
