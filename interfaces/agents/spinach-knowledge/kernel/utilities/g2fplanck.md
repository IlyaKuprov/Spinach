# kernel/utilities/g2fplanck.m

- Signature: `G=g2fplanck(spin_system,parameters)`

## Purpose

Returns gradient operators within the Fokker-Planck formalism used in the imaging module of Spinach. Syntax: G=g2fplanck(spin_system,parameters)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Parameters / inputs

- parameters.dims -a vector with one, two or three
- elements giving the dimensions
- of the box, metres
- parameters.npts -a vector with one, two or three
- elements giving number of points
- in each dimension of the box

## Outputs

- G -a cell array with the three gradient operators
- ordered as {Gx,Gy,Gz}, normalised to 1 T/m, empty
- matrices for non-exitent dimensions
- Note: gradients are assumed to be linear and centered on the
- middle of the sample.
- Note: the direct product order is Z(x)Y(x)X(x)Spin, this cor-
- responds to a column-wise vectorization of a 3D array
- with dimensions ordered as [X Y Z].
- Note: polyadic objects are returned, use inflate() to get the
- corresponding sparse matrix.

## Implementation structure

- Returns gradient operators within the Fokker-Planck formalism
- used in the imaging module of Spinach. Syntax:
- G=g2fplanck(spin_system,parameters)
- parameters.dims -a vector with one, two or three
- elements giving the dimensions
- of the box, metres
- parameters.npts -a vector with one, two or three
- elements giving number of points
- in each dimension of the box
- G -a cell array with the three gradient operators
- ordered as {Gx,Gy,Gz}, normalised to 1 T/m, empty
- matrices for non-exitent dimensions
