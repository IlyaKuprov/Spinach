# kernel/utilities/conmat.m

- Signature: `conmatrix=conmat(xyz,r0)`

## Purpose

Molecular connectivity matrix calculator with N*log(N) asymptotic complexity scaling with respect to the num- ber or atoms. Syntax: conmatrix=conmat(xyz,r0)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

## Parameters / inputs

- xyz -an array with N rows and three columns,
- giving the Cartesian coordinates of
- each particle
- r0 -the distance below which the particles
- are to be considered "connected"
- Output:
- conmatrix -a sparse logical matrix containing 1
- at the positions corresponding to the
- connected particles

## Implementation structure

- Molecular connectivity matrix calculator with N*log(N)
- asymptotic complexity scaling with respect to the num-
- ber or atoms. Syntax:
- conmatrix=conmat(xyz,r0)
- xyz -an array with N rows and three columns,
- giving the Cartesian coordinates of
- each particle
- r0 -the distance below which the particles
- are to be considered "connected"
- Output:
- conmatrix -a sparse logical matrix containing 1
- at the positions corresponding to the
