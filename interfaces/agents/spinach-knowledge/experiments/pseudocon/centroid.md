# experiments/pseudocon/centroid.m

- Signature: `[x,y,z]=centroid(probden,ranges)`

## Purpose

Finds the centre of mass point of a 3D probability density in a cube. Syntax: [x,y,z]=centroid(probden,ranges)

## Physical / mathematical content

- Paramagnetic-pseudocontact inference routines. The mathematics includes inverse problems, tensor parameterisation, interpolation, and regularisation.

## Numerical / algorithmic content

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
