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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 35-36: Check consistency; implemented by `grumble(alphas,betas,gammas,weights,ranks,sfun)`.
- Lines 38-39: Preallocate the answer; implemented by `grid_profile=zeros(size(ranks))`.
- Lines 41-42: Loop over spherical ranks; implemented by `for k=1:numel(ranks)`.
- Lines 44-45: Preallocate Wigner matrix; implemented by `D=zeros(2*ranks(k)+1,'like',1i)`.
- Lines 47-48: Loop over grid points; implemented by `parfor n=1:numel(alphas)`.
- Lines 53-54: Update grid profile; implemented by `if strcmp(sfun,'D_lmn')`.
- Lines 64-66: Update the user; implemented by `disp(['Spherical rank ' num2str(ranks(k)) ', residual ' sfun ' norm: ' num2str(grid_profile(k))])`.
- Lines 70-71: Do the plotting; implemented by `if nargout==0`.

### Control flow inferred from the code

- Line 42: `for` loop over `k=1:numel(ranks)`.
- Line 48: `parfor` loop over `n=1:numel(alphas)`.
- Line 54: conditional branch on `strcmp(sfun,'D_lmn')`.
- Line 71: conditional branch on `nargout==0`.

### Key state/data transformations

- Lines 39: computes `grid_profile` using `grid_profile=zeros(size(ranks))`.
- Lines 45: computes `D` using `D=zeros(2*ranks(k)+1,'like',1i)`.
- Lines 55: computes `grid_profile(k)` using `grid_profile(k)=norm(D,2)-krondelta(0,ranks(k))`.

### Local helper functions

- Line 80: `grumble()` — `function grumble(alphas,betas,gammas,weights,ranks,sfun)`.
  - Representative operation: `if (~isnumeric(alphas))||(~isreal(alphas))|| any(~isfinite(alphas))||(size(alphas,2)~=1)`.
  - Representative operation: `any(~isfinite(alphas))||(size(alphas,2)~=1)`.

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
