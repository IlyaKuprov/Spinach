# kernel/utilities/rwalk.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/rwalk.m`
- Signature: `eulers=rwalk(npts,tau_c,dt)`
- Total lines: 83

## Purpose

Random walk on SO(3), isotropic rotational diffusion. Syntax: eulers=rwalk(npts,tau_c,dt)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 29-30: Check consistency; implemented by `grumble(npts,tau_c,dt)`.
- Lines 32-33: Random unit jump sequence; implemented by `jump_angles=randn(npts,3)/sqrt(3)`.
- Lines 35-36: Time step and diffusion coefficient; implemented by `jump_angles=sqrt(dt/tau_c)*jump_angles`.
- Lines 38-39: Require jump angles to be small; implemented by `if mean(abs(jump_angles))>(pi/32)`.
- Lines 43-44: Compute DCM trajectory; implemented by `DCM=zeros(3,3,npts); DCM(:,:,1)=eye(3)`.
- Lines 52-53: Preallocate the answer; implemented by `eulers=zeros(npts,3)`.
- Lines 55-56: Compute Euler angle trajectory; implemented by `parfor n=1:npts`.

### Control flow inferred from the code

- Line 39: conditional branch on `mean(abs(jump_angles))>(pi/32)`.
- Line 45: `for` loop over `n=2:npts`.
- Line 56: `parfor` loop over `n=1:npts`.

### Key state/data transformations

- Lines 33: computes `jump_angles` using `jump_angles=randn(npts,3)/sqrt(3)`.
- Lines 44: computes `DCM` using `DCM=zeros(3,3,npts); DCM(:,:,1)=eye(3)`.
- Lines 46-48: computes `R` using `R=[ 0 1 0; -1 0 0; 0 0 0]*jump_angles(n,1)+ [ 0 0 1; 0 0 0; -1 0 0]*jump_angles(n,2)+ [ 0 0 0; 0 0 1; 0 -1 0]*jump_angles(n,3)`.
- Lines 49: computes `DCM(:,:,n)` using `DCM(:,:,n)=expm(R)*DCM(:,:,n-1)`.
- Lines 53: computes `eulers` using `eulers=zeros(npts,3)`.
- Lines 57: computes `[alp,bet,gam]` using `[alp,bet,gam]=dcm2euler(DCM(:,:,n))`.
- Lines 58: computes `eulers(n,:)` using `eulers(n,:)=[alp bet gam]`.

### Local helper functions

- Line 64: `grumble()` — `function grumble(npts,tau_c,dt)`.
  - Representative operation: `if (~isnumeric(npts))||(~isreal(npts))|| (~isscalar(npts))||(npts<1)||(mod(npts,1)~=0)`.
  - Representative operation: `(~isscalar(npts))||(npts<1)||(mod(npts,1)~=0)`.

## Parameters / inputs

- npts -number of points in the trajectory
- tau_c -isotropic rotational correlation
- time, seconds
- dt -inter-point spacing, seconds

## Outputs

- eulers -npts x 3 array of Euler angles for
- for each trajectory point, radians
- Note: the angles are NOT increments relative to the previous
- points, they are angles relative to the starting point
- of the trajectory.

## Implementation structure

- Random walk on SO(3), isotropic rotational diffusion. Syntax:
- eulers=rwalk(npts,tau_c,dt)
- npts -number of points in the trajectory
- tau_c -isotropic rotational correlation
- time, seconds
- dt -inter-point spacing, seconds
- eulers -npts x 3 array of Euler angles for
- for each trajectory point, radians
- Note: the angles are NOT increments relative to the previous
- points, they are angles relative to the starting point
- of the trajectory.
- Check consistency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `DCM()`, `jump_angles()`, `dcm2euler()`, `eulers()`, `isscalar()`.
