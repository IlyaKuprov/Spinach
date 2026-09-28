# tests/kernel/test_dynamic_grid_plot.m

- Signature: `result=test_dynamic_grid_plot()`

## Purpose

Regression test for `grid_plot()` using invisible figures.

## Physical / mathematical content

- Uses four points at the vertices of a regular tetrahedron, normalized onto the unit sphere, and their spherical Voronoi tessellation.

## Numerical / algorithmic content

- Computes the tessellation with `voronoisphere(xyz)`.
- With supplied tessera, numeric colours `(1:4).'`, and `options.dots=false`, checks that `grid_plot()` creates one patch per Voronoi cell, passes the numeric colour data to the patches, and creates no centre-dot line object.
- With tessera and options omitted, checks that `grid_plot()` generates the same number of patches, draws one line object for the default centre dots, and leaves the axes with `PlotBoxAspectRatioMode` set to `manual` for square plotting.
- Sets the default figure visibility to `off` during the test, creates each figure with `Visible='off'`, closes the figures, and restores the previous default visibility through an `onCleanup` handler.

## Outputs

- `result` — regression test result with explanatory messages.

## Attribution

- ilya.kuprov@weizmann.ac.il