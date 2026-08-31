# kernel/conventions/transforms/xyz2sph.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/conventions/transforms/xyz2sph.m`
- Signature: `[r,theta,phi] = xyz2sph(x,y,z)`
- Total lines: 56

## Purpose

Converts Cartesian coordinates [x y z] into spherical coordinates according to the ISO convention. Syntax: [r,theta,phi] = xyz2sph(x,y,z)

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 25-26: Check consistency; implemented by `grumble(x,y,z)`.
- Lines 28-29: Radius 0 <= r < Inf; implemented by `r=sqrt(x.^2+y.^2+z.^2)`.
- Lines 31-32: Inclination 0 <= theta <= pi; implemented by `theta=acos(z./r)`.
- Lines 34-35: Azimuth 0 <= phi < 2*pi; implemented by `phi=mod(atan2(y,x),2*pi)`.

### Key state/data transformations

- Lines 29: computes `r` using `r=sqrt(x.^2+y.^2+z.^2)`.
- Lines 32: computes `theta` using `theta=acos(z./r)`.
- Lines 35: computes `phi` using `phi=mod(atan2(y,x),2*pi)`.

### Local helper functions

- Line 40: `grumble()` — `function grumble(x,y,z)`.
  - Representative operation: `if (~isnumeric(x))||(~isnumeric(y))||(~isnumeric(z))`.
  - Representative operation: `error('all inputs must be numeric.')`.

## Parameters / inputs

- x,y,z -arrays of X, Y and Z coordinates

## Outputs

- r -array of radii
- theta -array of inclinations
- phi -array of azimuth values

## Implementation structure

- Converts Cartesian coordinates [x y z] into spherical
- coordinates according to the ISO convention. Syntax:
- [r,theta,phi] = xyz2sph(x,y,z)
- x,y,z -arrays of X, Y and Z coordinates
- r -array of radii
- theta -array of inclinations
- phi -array of azimuth values
- Check consistency
- Radius 0 <= r < Inf
- Inclination 0 <= theta <= pi
- Azimuth 0 <= phi < 2*pi
- Consistency enforcement

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `acos()`, `atan2()`, `all()`.
