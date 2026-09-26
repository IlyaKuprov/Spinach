# experiments/pseudocon/probmax.m

- Signature: `[x,y,z]=probmax(probden,ranges)`

## Purpose

Finds the maximum point of a 3D probability density in a cube. Syntax: [x,y,z]=probmax(probden,ranges)

## Physical / mathematical content

- Paramagnetic-pseudocontact inference routines. The mathematics includes inverse problems, tensor parameterisation, interpolation, and regularisation.

## Numerical / algorithmic content

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
