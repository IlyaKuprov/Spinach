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
