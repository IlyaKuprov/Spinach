# kernel/grids/grid_trian.m

- Signature: `[alps,bets,gams,whts,vorn]=grid_trian(type,n)`

## Purpose

Generate triangular spherical quadrature grids as described in [Appendix A.6](http://dx.doi.org/10.1016/j.jmr.2014.05.009).

## Parameters

- `type` — character string selecting `'asg'`, `'sophe'`, or `'stoll'`.
- `n` — positive integer point-count parameter.

## Outputs

- `alps` — alpha Euler angles in radians; all zero because these are two-angle grids.
- `bets` — beta Euler angles in radians.
- `gams` — gamma Euler angles in radians.
- `whts` — Voronoi tessellation body-angle weights, divided by `4*pi`.
- `vorn` — cell array of matrices containing Voronoi-polyhedron vertex coordinates.

## Method

Each grid is built from points in the first octant and extended by reflection to the sphere. The `'stoll'` construction averages three SOPHE-derived octant point sets. The `'sophe'` and `'stoll'` constructions explicitly add pole points. Voronoi tessellation and weights are computed when more than three outputs are requested or when no outputs are requested. With no outputs, the function draws a grid schematic.