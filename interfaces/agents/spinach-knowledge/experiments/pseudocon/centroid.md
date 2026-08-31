# experiments/pseudocon/centroid.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/pseudocon/centroid.m`
- Signature: `[x,y,z]=centroid(probden,ranges)`
- Total lines: 58

## Purpose

Finds the centre of mass point of a 3D probability density in a cube. Syntax: [x,y,z]=centroid(probden,ranges)

## Physical / mathematical content

- Paramagnetic-pseudocontact inference routines. The mathematics includes inverse problems, tensor parameterisation, interpolation, and regularisation.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 24-25: Check consistency; implemented by `grumble(probden,ranges)`.
- Lines 27-30: Get coordinate arrays; implemented by `[X,Y,Z]=ndgrid(linspace(ranges(1),ranges(2),size(probden,1)), linspace(ranges(3),ranges(4),size(probden,2)), linspace(ranges(5),ranges(6),size(probden,3)))`.
- Lines 32-33: Get the normalization; implemented by `n=trapz(trapz(trapz(probden)))`.
- Lines 35-36: Get centroid coordinates; implemented by `x=trapz(trapz(trapz(X.*probden)))/n`.

### Key state/data transformations

- Lines 28-30: computes `[X,Y,Z]` using `[X,Y,Z]=ndgrid(linspace(ranges(1),ranges(2),size(probden,1)), linspace(ranges(3),ranges(4),size(probden,2)), linspace(ranges(5),ranges(6),size(probden,3)))`.
- Lines 33: computes `n` using `n=trapz(trapz(trapz(probden)))`.
- Lines 36: computes `x` using `x=trapz(trapz(trapz(X.*probden)))/n`.
- Lines 37: computes `y` using `y=trapz(trapz(trapz(Y.*probden)))/n`.
- Lines 38: computes `z` using `z=trapz(trapz(trapz(Z.*probden)))/n`.

### Local helper functions

- Line 43: `grumble()` — `function grumble(probden,ranges)`.
  - Representative operation: `if (~isnumeric(ranges))||(~isreal(ranges))||(numel(ranges)~=6)`.
  - Representative operation: `error('ranges must be a real vector with six elements.')`.

## Parameters / inputs

- probden -probability density cube with dimensions
- ordered as [X Y Z]
- ranges -six-element vector giving axis extents
- as [xmin xmax ymin ymax zmin zmax]

## Outputs

- [x,y,z] -centre of mass coordinates

## Implementation structure

- Finds the centre of mass point of a 3D probability density
- in a cube. Syntax:
- [x,y,z]=centroid(probden,ranges)
- probden -probability density cube with dimensions
- ordered as [X Y Z]
- ranges -six-element vector giving axis extents
- as [xmin xmax ymin ymax zmin zmax]
- [x,y,z] -centre of mass coordinates
- Check consistency
- Get coordinate arrays
- Get the normalization
- Get centroid coordinates

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `ranges()`, `trapz()`, `ndims()`.
