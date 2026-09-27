# interfaces/comsol/mesh_plot.m

- Signature: `mesh_plot(spin_system,qscale,nodelabels)`

## Purpose

Plots the prepared 2D mesh geometry and, optionally, its velocity field and vertex labels. The mesh must already contain plotting data, normally prepared by `mesh_preplot`.

## Behavior

The routine draws triangle, rectangle, and Voronoi-cell outlines, sets equal axis scaling, and labels the axes as X and Y position in millimetres. When `qscale` is positive it adds a velocity quiver plot; when `nodelabels` is 1 it labels vertices with their indices. The edge-array plotting block is commented out in the source.

## Parameters / inputs

- `spin_system`: Spinach spin-system structure containing mesh and prepared plotting data.
- `qscale`: non-negative real scalar controlling velocity-arrow scaling; zero disables the quiver plot.
- `nodelabels`: 1 to show vertex numbers, or 0 to suppress them.

## Output

The function draws the plot; it has no returned output argument.

## Source

[Spinach Wiki: mesh_plot.m](https://spindynamics.org/wiki/index.php?title=mesh_plot.m)
