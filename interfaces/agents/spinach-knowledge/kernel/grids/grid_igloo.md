# kernel/grids/grid_igloo.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/grids/grid_igloo.m`
- Signature: `[alps,bets,gams,whts,vorn]=grid_igloo(n_long)`
- Total lines: 104

## Purpose

Igloo grid, as per Appendix A2 of

## Physical / mathematical content

- Quadrature and geometry utilities. These files generate spherical/SO(3) grids, Voronoi weights, and adaptive integration tools for orientation averaging.

## Numerical / algorithmic content

- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Syntax

```matlab
[alps,bets,gams,whts,vorn]=grid_igloo(n_long)
```

## Parameters / inputs

- n_long -number of longitudes in the grid

## Outputs

- alps -alpha Euler angles of the grid (radians),
- zeros because these are two-angle grids
- bets -beta Euler angles of the grid (radians)
- gams -gamma Euler angles of the grid (radians)
- whts -Voronoi tessellation body angle weights
- vorn -a cell array of matrices containing the
- coordinates of the vertices of the Voro-
- noi polyhedra
- If no outputs are requested, a schematic is drawn.

## Implementation structure

- Igloo grid, as per Appendix A2 of
- [alps,bets,gams,whts,vorn]=grid_igloo(n_long)
- n_long -number of longitudes in the grid
- alps -alpha Euler angles of the grid (radians),
- zeros because these are two-angle grids
- bets -beta Euler angles of the grid (radians)
- gams -gamma Euler angles of the grid (radians)
- whts -Voronoi tessellation body angle weights
- vorn -a cell array of matrices containing the
- coordinates of the vertices of the Voro-
- noi polyhedra
- If no outputs are requested, a schematic is drawn.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `bets()`, `gams()`, `cell2mat()`, `alps_bets_gams()`, `voronoisphere()`, `grid_plot()`, `isscalar()`.
