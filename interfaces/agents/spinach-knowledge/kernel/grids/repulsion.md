# kernel/grids/repulsion.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/grids/repulsion.m`
- Signature: `[alphas,betas,gammas,weights]=repulsion(npoints,ndims,niter)`
- Total lines: 140

## Purpose

Generates repulsion grids on a unit hypersphere. See the paper by Bak and Nielsen (http://dx.doi.org/10.1006/jmre.1996.1087) to get further information on the algorithm involved. Syntax: [alphas,betas,gammas,weights]=repulsion(npoints,ndims,niter)

## Physical / mathematical content

- Quadrature and geometry utilities. These files generate spherical/SO(3) grids, Voronoi weights, and adaptive integration tools for orientation averaging.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- npoints -number of points in the resulting spherical grid
- ndims -hypersphere dimension: 2 returns a single-angle
- (beta) grid, 3 returns a two-angle grid (alpha,
- beta), 4 returns a three-angle (alpha,beta,gam-
- ma) spherical grid
- niter -number of repulsion interations (simple clipped
- gradient descent at the moment)

## Outputs

- alphas -alpha Euler angles of the grid, in radians,
- zeros for two-angle grids
- betas -beta Euler angles of the grid, in radians
- gammas -gamma Euler angles of the grid, in radians,
- zeros for single-angle grids
- weights -point weights of the grid
- Note: uniform weights are assigned at the moment, use the supp-
- lied SHREWD function to generate optimal weights.

## Implementation structure

- Generates repulsion grids on a unit hypersphere. See the paper by
- Bak and Nielsen (http://dx.doi.org/10.1006/jmre.1996.1087) to get
- further information on the algorithm involved. Syntax:
- [alphas,betas,gammas,weights]=repulsion(npoints,ndims,niter)
- npoints -number of points in the resulting spherical grid
- ndims -hypersphere dimension: 2 returns a single-angle
- (beta) grid, 3 returns a two-angle grid (alpha,
- beta), 4 returns a three-angle (alpha,beta,gam-
- ma) spherical grid
- niter -number of repulsion interations (simple clipped
- gradient descent at the moment)
- alphas -alpha Euler angles of the grid, in radians,

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `dist_vecs()`, `num2str()`, `cart2pol()`, `cart2sph()`, `kfigure()`, `plot3()`, `qter2euler()`, `ismember()`.
