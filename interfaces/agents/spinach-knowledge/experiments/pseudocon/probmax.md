# experiments/pseudocon/probmax.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/pseudocon/probmax.m`
- Signature: `[x,y,z]=probmax(probden,ranges)`
- Total lines: 57

## Purpose

Finds the maximum point of a 3D probability density in a cube. Syntax: [x,y,z]=probmax(probden,ranges)

## Physical / mathematical content

- Paramagnetic-pseudocontact inference routines. The mathematics includes inverse problems, tensor parameterisation, interpolation, and regularisation.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- probden -probability density cube with dimensions
- ordered as [X Y Z]
- ranges -six-element vector giving axis extents
- as [xmin xmax ymin ymax zmin zmax]

## Outputs

- [x,y,z] -maximum point coordinates

## Implementation structure

- Finds the maximum point of a 3D probability density in a
- cube. Syntax:
- [x,y,z]=probmax(probden,ranges)
- probden -probability density cube with dimensions
- ordered as [X Y Z]
- ranges -six-element vector giving axis extents
- as [xmin xmax ymin ymax zmin zmax]
- [x,y,z] -maximum point coordinates
- Check consistency
- Get coordinate arrays
- Get the max
- Get maximum coordinates

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `ranges()`, `probden()`, `ndims()`.
