# examples/fundamentals/quadratures/grid_leb_vs_etc.m

- Signature: `grid_leb_vs_etc()`

## Purpose

Compare spherical quadrature integration errors across Lebedev and several other grid constructions, as a function of spherical-harmonic rank. The source title calls this a heuristic-versus-Lebedev “bake-off”; it does not establish a general ranking beyond the plotted test.

## Physical / mathematical content

The benchmark tests spherical integration of Y_lm for ranks 4:2:60. Lebedev supplies its stored weights; the repulsion, ZCWn, Igloo, Stoll, ASG, and SOPHE grids are evaluated with the weights used in the source (including Voronoi weights for the heuristic grids).

## Numerical / algorithmic content

For each grid, grid_test returns an integration-error profile. The profiles are plotted against rank on a logarithmic error axis; no numerical result is asserted here without running the example.

## Implementation structure

- Load the rank-29 Lebedev grid and evaluate it over the rank sequence.
- Generate repulsion, ZCWn, Igloo, Stoll, ASG, and SOPHE grids; compute Voronoi weights where specified and normalize them.
- Plot and label the error profiles for comparison.
