# kernel/grids/grid_test.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/grids/grid_test.m`
- Signature: `grid_profile=grid_test(alphas,betas,gammas,weights,ranks,sfun)`
- Total lines: 118

## Purpose

Plots grid integration quality as a function of spherical rank. The quality is defined as the norm of the residual of spherical harmon- ics or Wigner functions integrated using the grid provided. Syntax: grid_profile=grid_test(alphas,betas,gammas,weights,max_rank,sfun)

## Physical / mathematical content

- Quadrature and geometry utilities. These files generate spherical/SO(3) grids, Voronoi weights, and adaptive integration tools for orientation averaging.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- alphas -alpha Euler angles of the grid, in radians,
- zeros for single-angle grids
- betas -beta Euler angles of the grid, in radians
- gammas -gamma Euler angles of the grid, in radians,
- zeros for two-angle grids
- weights -point weights of the grid
- ranks -spherical ranks to consider
- sfun -spherical function type: for three-angle
- grids use 'D_lmn', for two-angle grids use
- 'Y_lm', for single-angle grids use 'Y_l0'.

## Outputs

- grid_profile -a vector of residual norms in each spherical rank

## Implementation structure

- Plots grid integration quality as a function of spherical rank. The
- quality is defined as the norm of the residual of spherical harmon-
- ics or Wigner functions integrated using the grid provided. Syntax:
- grid_profile=grid_test(alphas,betas,gammas,weights,max_rank,sfun)
- alphas -alpha Euler angles of the grid, in radians,
- zeros for single-angle grids
- betas -beta Euler angles of the grid, in radians
- gammas -gamma Euler angles of the grid, in radians,
- zeros for two-angle grids
- weights -point weights of the grid
- ranks -spherical ranks to consider
- sfun -spherical function type: for three-angle

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `ranks()`, `weights()`, `wigner()`, `alphas()`, `betas()`, `gammas()`, `strcmp()`, `grid_profile()`, `krondelta()`, `num2str()`, `kfigure()`, `kxlabel()`, `kylabel()`, `any()`, `isvector()`.
