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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 23-24: Check consistency; implemented by `grumble(r1,r2,r3)`.
- Lines 26-27: Not particularly hard; implemented by `r12=r1+r2; r12=r12/norm(r12,2)`.

### Key state/data transformations

- Lines 27: computes `r12` using `r12=r1+r2; r12=r12/norm(r12,2)`.
- Lines 28: computes `r23` using `r23=r2+r3; r23=r23/norm(r23,2)`.
- Lines 29: computes `r31` using `r31=r3+r1; r31=r31/norm(r31,2)`.

### Local helper functions

- Line 34: `grumble()` — `function grumble(r1,r2,r3)`.
  - Representative operation: `if (~isnumeric(r1))||(~isreal(r1))||(numel(r1)~=3)|| (~isnumeric(r2))||(~isreal(r2))||(numel(r2)~=3)|| (~isnumeric(r3))||(~isreal(r3))||(numel(r3)~=3)`.
  - Representative operation: `(~isnumeric(r2))||(~isreal(r2))||(numel(r2)~=3)|| (~isnumeric(r3))||(~isreal(r3))||(numel(r3)~=3)`.

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
