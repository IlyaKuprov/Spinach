# kernel/grids/grid_igloo.m

- Signature: `[alps,bets,gams,whts,vorn]=grid_igloo(n_long)`

## Purpose

Generate an igloo grid as described in Appendix A2 of [the cited paper](http://dx.doi.org/10.1016/j.jmr.2014.05.009).

## Syntax

```matlab
[alps,bets,gams,whts,vorn]=grid_igloo(n_long)
```

## Parameters / inputs

- `n_long` — number of longitudes in the grid; must be a positive integer.

## Outputs

- `alps` — alpha Euler angles in radians; zero because this is a two-angle grid.
- `bets` — beta Euler angles in radians.
- `gams` — gamma Euler angles in radians.
- `whts` — normalized Voronoi tessellation body-angle weights.
- `vorn` — cell array of matrices containing Voronoi polyhedron vertex coordinates.

If no outputs are requested, a schematic of the grid is drawn.

## Implementation

Beta angles are equally spaced from 0 to π. At each beta angle, the number of gamma points is `max(1,floor(2*(n_long-1)*sin(beta)+0.5))`; gamma angles are equally spaced around the circle without repeating the endpoint. When weights are requested or no outputs are requested, the points are converted to Cartesian coordinates and passed to `voronoisphere`. The resulting weights are divided by `4*pi`.

[Source documentation](https://spindynamics.org/wiki/index.php?title=grid_igloo.m)
