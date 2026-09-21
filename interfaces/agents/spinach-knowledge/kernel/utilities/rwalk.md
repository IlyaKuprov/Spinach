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
