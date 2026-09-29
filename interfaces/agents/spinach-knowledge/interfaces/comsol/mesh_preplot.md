# interfaces/comsol/mesh_preplot.m

[Canonical MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/interfaces/comsol/mesh_preplot.m) · [Spinach Wiki: mesh_preplot.m](https://spindynamics.org/wiki/index.php?title=mesh_preplot.m)

## Purpose and return value

`mesh=mesh_preplot(mesh)` converts mesh connectivity and Voronoi-cell data into coordinate arrays for later plotting, stores them in `mesh.plot`, and returns the updated mesh structure. It prepares data only; it does not draw a plot or call COMSOL.

## Accepted data and dependencies

The input is a mesh structure with `mesh.x`/`mesh.y` coordinates, `mesh.idx.edges`, `mesh.idx.triangles`, `mesh.idx.rectangles` connectivity, and `mesh.vor.cells` plus `mesh.vor.vertices` from a Voronoi tessellation (normally produced by `mesh_vorn`). Connectivity rows index the coordinate vectors; each cell list indexes rows of `mesh.vor.vertices`. MATLAB array and cell-array operations are used. The routine checks only that `mesh.idx` and `mesh.vor` exist, then assumes the listed fields and compatible indices are present.

## Transformations

For each edge it stores the two endpoint coordinates followed by `NaN` in `mesh.plot.edg_a` and `mesh.plot.edg_b`. Triangle and rectangle outlines are closed by repeating their first vertex and separated by `NaN`. Voronoi polygons are likewise closed; a `NaN` gap separates consecutive cells. The output fields are `edg_a`, `edg_b`, `tri_a`, `tri_b`, `rec_a`, `rec_b`, `vor_a`, and `vor_b` under `mesh.plot`. Coordinates are copied directly, without unit conversion, and no plotting scale is applied. The function returns the modified mesh, not a figure or status code.
