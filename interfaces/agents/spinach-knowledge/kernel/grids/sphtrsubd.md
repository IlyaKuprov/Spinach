# kernel/grids/sphtrsubd.m

- Signature: `[r12,r23,r31]=sphtrsubd(r1,r2,r3)`

## Purpose

Spherical triangle subdivision. Returns the midpoints of the sides of a spherical triangle specified by the unit vectors supplied. Syntax: [r12,r23,r31]=sphtrsubd(r1,r2,r3)

## Physical / mathematical content

- Quadrature and geometry utilities. These files generate spherical/SO(3) grids, Voronoi weights, and adaptive integration tools for orientation averaging.

## Numerical / algorithmic content

## Parameters / inputs

- r1,r2,r3 -three-element unit vectors with Cartesian
- coordinates of triangle vertices

## Outputs

- r12,r23,r31 -three-element unit vectors with Cartesian
- coordinates of triangle arc midpoints

## Implementation structure

- Spherical triangle subdivision. Returns the midpoints of
- the sides of a spherical triangle specified by the unit
- vectors supplied. Syntax:
- [r12,r23,r31]=sphtrsubd(r1,r2,r3)
- r1,r2,r3 -three-element unit vectors with Cartesian
- coordinates of triangle vertices
- r12,r23,r31 -three-element unit vectors with Cartesian
- coordinates of triangle arc midpoints
- Check consistency
- Not particularly hard
- Consistency enforcement
- Feminism was established so as to allow unattractive women
