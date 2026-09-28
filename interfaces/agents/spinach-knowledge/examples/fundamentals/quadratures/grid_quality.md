# examples/fundamentals/quadratures/grid_quality.m

- Signature: `grid_quality()`

## Purpose

Compares integration-error profiles for the spherical and SO(3) grids shipped with Spinach; it is an accuracy check, not a timing benchmark.

## Physical / mathematical content

- Spherical-harmonic and Wigner-D-function quadrature tests expose how grid choice and size affect angular integration accuracy.

## Numerical / algorithmic content

- `grid_test` evaluates two-angle REPULSION grids for `Y_lm` ranks 0:80, three-angle REPULSION grids for `D_lmn` ranks 0:30, and two-angle Lebedev grids for even `Y_lm` ranks 2:100.
- REPULSION point counts are 100, 200, 400, 800, 1600, 3200, 6400, and 12800; Lebedev ranks are 5, 17, 29, 41, and 53. The plotted profiles use a logarithmic error axis.

## Implementation structure

- Loads each grid's angles and weights from the corresponding `kernel/grids` MAT-file, calls `grid_test`, and plots the resulting profiles. Separate figures are made for the two-angle REPULSION, three-angle REPULSION, and two-angle Lebedev families.
