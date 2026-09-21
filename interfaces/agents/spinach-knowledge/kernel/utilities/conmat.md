# kernel/utilities/conmat.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/conmat.m`
- Signature: `conmatrix=conmat(xyz,r0)`
- Total lines: 106

## Purpose

Molecular connectivity matrix calculator with N*log(N) asymptotic complexity scaling with respect to the num- ber or atoms. Syntax: conmatrix=conmat(xyz,r0)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `xyz()`, `false()`, `x_sorted()`, `x_index()`, `y_sorted()`, `y_index()`, `z_sorted()`, `z_index()`, `row()`, `col()`, `conmatrix()`.
