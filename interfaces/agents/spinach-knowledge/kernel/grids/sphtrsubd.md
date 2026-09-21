# kernel/grids/sphtrsubd.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/grids/sphtrsubd.m`
- Signature: `[r12,r23,r31]=sphtrsubd(r1,r2,r3)`
- Total lines: 59

## Purpose

Spherical triangle subdivision. Returns the midpoints of the sides of a spherical triangle specified by the unit vectors supplied. Syntax: [r12,r23,r31]=sphtrsubd(r1,r2,r3)

## Physical / mathematical content

- Quadrature and geometry utilities. These files generate spherical/SO(3) grids, Voronoi weights, and adaptive integration tools for orientation averaging.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `sphtarea()`, `arclength()`.
