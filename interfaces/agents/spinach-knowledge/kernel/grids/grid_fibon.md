# kernel/grids/grid_fibon.m

- Signature: `[alps,bets,gams,whts,vorn]=grid_fibon(type,parm)`

## Purpose

Generate Fibonacci-type spherical quadrature grids, as described in Appendix A.5 of http://dx.doi.org/10.1016/j.jmr.2014.05.009.

## Parameters / inputs

- `type`: `'fib'`, `'zcw'`, or `'zcwn'`.
- `parm`: positive integer point-count parameter. The grids have `2*parm+1`, `fibonacci(parm+2)`, or `parm` points, respectively.

## Outputs

- `alps`: alpha Euler angles (radians); zero for these two-angle grids.
- `bets`: beta Euler angles (radians).
- `gams`: gamma Euler angles (radians).
- `whts`: Voronoi tessellation body-angle weights, normalized by `4*pi`.
- `vorn`: cell array of matrices containing Voronoi-polyhedron vertex coordinates.

## Implementation

The `'fib'` and `'zcwn'` grids use the golden ratio to place points; `'zcw'` uses Fibonacci numbers. Voronoi tessellation is computed when weights are requested or the function is called without outputs. With no outputs, the function plots a schematic of the grid.