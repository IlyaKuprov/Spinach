# kernel/conventions/transforms/xyz2sph.m

- Signature: `[r,theta,phi] = xyz2sph(x,y,z)`

## Purpose

Converts Cartesian coordinates [x y z] into spherical coordinates according to the ISO convention. Syntax: [r,theta,phi] = xyz2sph(x,y,z)

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

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
