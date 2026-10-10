# examples/fundamentals/quadratures/grid_quality.m

## Purpose and question

This example plots how the error profiles returned by `grid_test` vary with grid size and spherical rank for three shipped grid families. The first two panels assess two- and three-angle REPULSION grids; the third assesses two-angle Lebedev grids. It is an error-profile visualisation, not a timing benchmark.

## Inputs and numerical method

Each grid is loaded from its kernel MAT-file with Euler-angle arrays and weights. Grid sizes for both REPULSION families are 100, 200, 400, 800, 1600, 3200, 6400, and 12800 points. The two-angle grids are passed to `grid_test` for ranks 0 through 80 with `Y_lm`. The three-angle grids use ranks 0 through 30 with `D_lmn`. Their error profiles are plotted on logarithmic y axes.

The Lebedev panel loads two-angle grids at ranks 5, 17, 29, 41, and 53 and evaluates `Y_lm` at the even ranks 2 through 100. It explicitly plots the even-rank coordinates alongside the corresponding profile values. The REPULSION panels plot `grid_profile(2:end)` against MATLAB's default sample index; the source labels the horizontal axis spherical rank but does not pass explicit rank coordinates to those plot calls.

## Output and limitations

The function creates three figures with integration error on the y axis and adds a legend entry per grid size or Lebedev rank. No threshold, pass/fail test, or numerical results are encoded in the example; actual accuracy values must come from running it with the shipped grid data and helper routines. In particular, the REPULSION-panel horizontal positions should not be read as an explicit mapping to requested ranks without checking `grid_test`'s output convention. The source reports no citation.

Source: [grid_quality.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/quadratures/grid_quality.m).
