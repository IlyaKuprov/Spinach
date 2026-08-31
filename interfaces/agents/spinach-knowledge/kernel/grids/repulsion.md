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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 41-42: Check consistency; implemented by `grumble(npoints,ndims,niter)`.
- Lines 44-45: Generate guess points; implemented by `R=rand(ndims,npoints)-0.5`.
- Lines 47-48: Start the repulsion loop; implemented by `for m=1:niter`.
- Lines 50-51: Make space for distance vector array; implemented by `Rd=reshape(R,[ndims 1 npoints])`.
- Lines 53-55: Compute distance vector array; implemented by `dist_vecs=permute(Rd,[1 2 3])- permute(Rd,[1 3 2])`.
- Lines 57-58: Normalise distance vectors; implemented by `dist_vecs=dist_vecs./sqrt(sum(dist_vecs.^2,1))`.
- Lines 60-61: Eliminate self-interactions; implemented by `dist_vecs(~isfinite(dist_vecs))=0`.
- Lines 63-64: Get scalar products; implemented by `scalar_prods=reshape(R'*R,[1 npoints npoints])`.
- Lines 66-67: Get tangent forces; implemented by `F=sum(dist_vecs.*scalar_prods,3)`.
- Lines 69-70: Move under tangent forces; implemented by `R_new=R-ndims*F/npoints`.
- Lines 72-73: Reproject onto unit sphere; implemented by `R_new=R_new./sqrt(sum(R_new.^2,1))`.
- Lines 75-76: Report the difference; implemented by `max_diff=max(sqrt(sum((R-R_new).^2,2)))`.
- Lines 80-81: Close the loop; implemented by `R=R_new`.
- Lines 85-86: Get points; implemented by `switch ndims`.
- Lines 90-91: In 2D case return polar angles; implemented by `[phi,~]=cart2pol(R(1,:),R(2,:))`.
- Lines 96-97: In 3D case return spherical angles; implemented by `[phi,theta,~]=cart2sph(R(1,:),R(2,:),R(3,:))`.
- Lines 100-101: Display the grid; implemented by `if nargout==0`.
- Lines 110-111: In 4D case return Euler angles; implemented by `qter.u=R(1,:)'; qter.i=R(2,:)'`.

### Control flow inferred from the code

- Line 48: `for` loop over `m=1:niter`.
- Line 86: dispatches on `ndims`; cases `2`, `3`, `4`.
- Line 101: conditional branch on `nargout==0`.

### Key state/data transformations

- Lines 45: computes `R` using `R=rand(ndims,npoints)-0.5`.
- Lines 51: computes `Rd` using `Rd=reshape(R,[ndims 1 npoints])`.
- Lines 54-55: computes `dist_vecs` using `dist_vecs=permute(Rd,[1 2 3])- permute(Rd,[1 3 2])`.
- Lines 61: computes `dist_vecs(~isfinite(dist_vecs))` using `dist_vecs(~isfinite(dist_vecs))=0`.
- Lines 64: computes `scalar_prods` using `scalar_prods=reshape(R'*R,[1 npoints npoints])`.
- Lines 67: computes `F` using `F=sum(dist_vecs.*scalar_prods,3)`.
- Lines 70: computes `R_new` using `R_new=R-ndims*F/npoints`.
- Lines 76: computes `max_diff` using `max_diff=max(sqrt(sum((R-R_new).^2,2)))`.
- Lines 91: computes `[phi,~]` using `[phi,~]=cart2pol(R(1,:),R(2,:))`.
- Lines 92: computes `betas` using `betas=phi'; alphas=0*betas; gammas=0*betas`.
- Lines 97: computes `[phi,theta,~]` using `[phi,theta,~]=cart2sph(R(1,:),R(2,:),R(3,:))`.
- Lines 111: computes `qter.u` using `qter.u=R(1,:)'; qter.i=R(2,:)'`.
- Lines 112: computes `qter.j` using `qter.j=R(3,:)'; qter.k=R(4,:)'`.
- Lines 113: computes `[alphas,betas,gammas]` using `[alphas,betas,gammas]=qter2euler(qter)`.
- Lines 118: computes `weights` using `weights=ones(npoints,1)/npoints`.

### Local helper functions

- Line 123: `grumble()` — `function grumble(npoints,ndims,niter)`.
  - Representative operation: `if (~isnumeric(npoints))||(~isreal(npoints))||(numel(npoints)~=1)||(npoints<1)||(mod(npoints,1)~=0)`.
  - Representative operation: `error('npoints must be a positive integer.')`.

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
