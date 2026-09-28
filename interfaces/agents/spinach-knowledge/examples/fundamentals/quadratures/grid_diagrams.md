# examples/fundamentals/quadratures/grid_diagrams.m

- Signature: `grid_diagrams()`

## Purpose

Create six spherical-grid diagrams for a book, with each panel labelled by grid type and point count.

## Physical / mathematical content

The function visualizes nodes on the unit sphere; it does not perform a quadrature-accuracy comparison. Cartesian plotting coordinates are formed from the spherical angles betas and gammas.

## Numerical / algorithmic content

The six displayed sets are a Fibonacci sequence grid (643 points), an icosahedral grid (642), an octahedral Stoll grid (678), a repulsion grid (620), an Igloo grid (616), and the rank-41 Lebedev grid (590).

## Implementation structure

- Create a compact 2-by-3 tiled figure.
- Generate or load each named grid, map angular coordinates to Cartesian sphere coordinates where needed, and plot the points with grid_plot.
- Add the grid labels and point counts shown in the panels.
