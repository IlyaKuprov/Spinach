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
