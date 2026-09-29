# examples/fundamentals/quadratures/grid_diagrams.m

## Purpose

This is a plotting example for a six-panel visual comparison of spherical-grid point patterns. It illustrates grid constructions and their displayed point counts; it does not measure integration accuracy.

## Grid inputs and construction

The function creates a compact 2-by-3 tiled figure. Its panels use these source-defined inputs and labels:

- Fibonacci sequence grid: `grid_fibon('fib',321)`, labelled 643 points.
- Icosahedral subdivision: angles loaded from `icos_2ang_642pts.mat`, labelled 642 points.
- Octahedral subdivision: `grid_trian('stoll',13)`, labelled 678 points.
- Repulsion grid: `repulsion(620,3,200)`, labelled 620 points.
- Igloo grid: `grid_igloo(23)`, labelled 616 points.
- Lebedev grid: angles loaded from `leb_2ang_rank_41.mat`, labelled 590 points.

For the loaded icosahedral and Lebedev angles, the displayed unit-sphere coordinates are formed from beta and gamma as `(sin(beta)*cos(gamma), sin(beta)*sin(gamma), cos(beta))`. The helper routines and plotting functions generate or display the other point sets.

## Output and limitations

The output is a MATLAB figure with panels marked SEQ, ICO, OCT, OPT, NAT, and LEB and with the point-count text shown above. Those counts are labels supplied by this example; it does not independently verify them, load quadrature weights, or test a spherical integral. The MAT-file paths are relative to the example's expected working directory. The source gives no numerical quality result or citation.

Source: [grid_diagrams.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/quadratures/grid_diagrams.m).
