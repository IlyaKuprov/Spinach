# interfaces/comsol/mesh_preplot.m

- Signature: `mesh=mesh_preplot(mesh)`

## Purpose

Prepares mesh geometry arrays for plotting, including its Voronoi cells.

## Behavior

Using the mesh coordinates and connectivity, the routine builds coordinate arrays for edges, triangles, and rectangles, then constructs the boundary arrays for the Voronoi cells. NaN separators distinguish separate plotted segments, and each Voronoi boundary is closed by returning to its first vertex. The arrays are stored under `mesh.plot` for subsequent plotting.

## Parameters / inputs

- `mesh`: Spinach mesh object containing indexing data and a Voronoi tessellation in `mesh.vor`.

## Output

- `mesh`: updated mesh object with `plot.edg_a`, `plot.edg_b`, `plot.tri_a`, `plot.tri_b`, `plot.rec_a`, `plot.rec_b`, `plot.vor_a`, and `plot.vor_b` arrays.

## Source

[Spinach Wiki: mesh_preplot.m](https://spindynamics.org/wiki/index.php?title=mesh_preplot.m)
