# kernel/grids/grid_plot.m

- Signature: `grid_plot(x,y,z,vorn,c,options)`

## Purpose

Plots the Voronoi tessellation of a cloud of points on a sphere.

## Parameters / inputs

- `x`, `y`, `z`: Column vectors of Cartesian grid-point coordinates; they must be real, finite, numeric, and the same size.
- `vorn`: Cell array of Voronoi tessera. If omitted or empty, it is computed with `voronoisphere`.
- `c`: Tessellation-face colour. If omitted or empty, faces are white; a character value names a colour, while numeric values are applied by face index.
- `options.dots`: Whether to plot black dots at grid points; defaults to `true`.

## Outputs

- Plots a figure; returns no value.

## Implementation structure

The function plots the optional centre dots and draws each tessellation face with `patch`. It sets square axes, a camera position, and plot limits, and hides tick marks.

Source: <https://spindynamics.org/wiki/index.php?title=grid_plot.m>