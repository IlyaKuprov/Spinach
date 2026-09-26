# kernel/grids/grid_trian.m

- Signature: `[alps,bets,gams,whts,vorn]=grid_trian(type,n)`

## Purpose

Triangular spherical quadrature grids, as per Appendix A.6 of (http://dx.doi.org/10.1016/j.jmr.2014.05.009). Syntax: [alps,bets,gams,whts,vorn]=grid_trian(type,n)

## Physical / mathematical content

- Quadrature and geometry utilities. These files generate spherical/SO(3) grids, Voronoi weights, and adaptive integration tools for orientation averaging.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Parameters / inputs

- type -'asg', 'sophe', or 'stoll'
- n -point count parameter

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

- Triangular spherical quadrature grids, as per Appendix A.6 of
- (http://dx.doi.org/10.1016/j.jmr.2014.05.009). Syntax:
- [alps,bets,gams,whts,vorn]=grid_trian(type,n)
- type -'asg', 'sophe', or 'stoll'
- n -point count parameter
- alps -alpha Euler angles of the grid (radians),
- zeros because these are two-angle grids
- bets -beta Euler angles of the grid (radians)
- gams -gamma Euler angles of the grid (radians)
- whts -Voronoi tessellation body angle weights
- vorn -a cell array of matrices containing the
- coordinates of the vertices of the Voro-
