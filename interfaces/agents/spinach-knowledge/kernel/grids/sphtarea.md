# kernel/grids/sphtarea.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/grids/sphtarea.m`
- Signature: `S=sphtarea(r1,r2,r3,sflag)`
- Total lines: 68

## Purpose

Area of the curvilinear triangle on the unit sphere defined by the vertex coordinates supplied. Syntax: S=sphtarea(r1,r2,r3,sflag)

## Physical / mathematical content

- Quadrature and geometry utilities. These files generate spherical/SO(3) grids, Voronoi weights, and adaptive integration tools for orientation averaging.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- r1,r2,r3 -three-element unit vectors with Cartesian
- coordinates of triangle vertices
- sflag -'signed' would take into account surface
- normal direction, 'unsigned' (default)
- would always return a positive area

## Outputs

- S -spherical triangle surface area

## Implementation structure

- Area of the curvilinear triangle on the unit sphere defined
- by the vertex coordinates supplied. Syntax:
- S=sphtarea(r1,r2,r3,sflag)
- r1,r2,r3 -three-element unit vectors with Cartesian
- coordinates of triangle vertices
- sflag -'signed' would take into account surface
- normal direction, 'unsigned' (default)
- would always return a positive area
- S -spherical triangle surface area
- Default to unsigned area
- Check consistency
- Stretch the vectors

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `exist()`, `grumble()`, `atan2()`, `dot()`, `strcmp()`, `arclength()`, `ischar()`, `ismember()`.
