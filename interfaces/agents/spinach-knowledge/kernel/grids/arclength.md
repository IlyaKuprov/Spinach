# kernel/grids/arclength.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/grids/arclength.m`
- Signature: `sig=arclength(r1,r2)`
- Total lines: 50

## Purpose

Arc length between two points on the unit sphere specified by the unit vectors supplied. Syntax: sig=arclength(r1,r2)

## Physical / mathematical content

- Quadrature and geometry utilities. These files generate spherical/SO(3) grids, Voronoi weights, and adaptive integration tools for orientation averaging.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 21-22: Check consistency; implemented by `grumble(r1,r2)`.
- Lines 24-25: Normalise the vectors; implemented by `r1=r1/norm(r1,2); r1=r1(:)`.
- Lines 28-29: Get the arc length; implemented by `sig=atan2(norm(cross(r1,r2),2),dot(r1,r2))`.

### Key state/data transformations

- Lines 25: computes `r1` using `r1=r1/norm(r1,2); r1=r1(:)`.
- Lines 26: computes `r2` using `r2=r2/norm(r2,2); r2=r2(:)`.
- Lines 29: computes `sig` using `sig=atan2(norm(cross(r1,r2),2),dot(r1,r2))`.

### Local helper functions

- Line 34: `grumble()` — `function grumble(r1,r2)`. Ah, there's nothing more exciting than science. You get all the
  - Representative operation: `if (~isnumeric(r1))||(~isreal(r1))||(numel(r1)~=3)|| (~isnumeric(r2))||(~isreal(r2))||(numel(r2)~=3)`.
  - Representative operation: `(~isnumeric(r2))||(~isreal(r2))||(numel(r2)~=3)`.

## Parameters / inputs

- r1,r2 -three-element unit vectors with Cartesian
- coordinates of the arc endpoints

## Outputs

- sig -arc length

## Implementation structure

- Arc length between two points on the unit sphere specified by
- the unit vectors supplied. Syntax:
- sig=arclength(r1,r2)
- r1,r2 -three-element unit vectors with Cartesian
- coordinates of the arc endpoints
- sig -arc length
- Check consistency
- Normalise the vectors
- Get the arc length
- Consistency enforcement
- Ah, there's nothing more exciting than science. You get all the
- fun of sitting still, being quiet, writing down numbers, paying

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `atan2()`, `cross()`, `dot()`.
