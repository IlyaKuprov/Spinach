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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 25-26: Default to unsigned area; implemented by `if ~exist('sflag','var'), sflag='unsigned'; end`.
- Lines 28-29: Check consistency; implemented by `grumble(r1,r2,r3,sflag)`.
- Lines 31-32: Stretch the vectors; implemented by `r1=r1(:); r2=r2(:); r3=r3(:)`.
- Lines 34-35: Get signed area; implemented by `S=2*atan2(det([r1 r2 r3]),dot(r1,r2)+dot(r2,r3)+dot(r3,r1)+1)`.
- Lines 37-38: Kill the sign if appropriate; implemented by `if strcmp(sflag,'unsigned'), S=abs(S); end`.

### Control flow inferred from the code

- Line 26: conditional branch on `~exist('sflag','var'), sflag='unsigned'; end`.
- Line 38: conditional branch on `strcmp(sflag,'unsigned'), S=abs(S); end`.

### Key state/data transformations

- Lines 32: computes `r1` using `r1=r1(:); r2=r2(:); r3=r3(:)`.
- Lines 35: computes `S` using `S=2*atan2(det([r1 r2 r3]),dot(r1,r2)+dot(r2,r3)+dot(r3,r1)+1)`.

### Local helper functions

- Line 43: `grumble()` — `function grumble(r1,r2,r3,sflag)`.
  - Representative operation: `if (~isnumeric(r1))||(~isreal(r1))||(numel(r1)~=3)|| (~isnumeric(r2))||(~isreal(r2))||(numel(r2)~=3)|| (~isnumeric(r3))||(~isreal(r3))||(numel(r3)~=3)`.
  - Representative operation: `(~isnumeric(r2))||(~isreal(r2))||(numel(r2)~=3)|| (~isnumeric(r3))||(~isreal(r3))||(numel(r3)~=3)`.

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
