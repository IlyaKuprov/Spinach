# kernel/utilities/conmat.m

- Signature: `conmatrix=conmat(xyz,r0)`

## Purpose

Calculates a molecular connectivity matrix from Cartesian coordinates and a distance threshold. The source describes `N*log(N)` asymptotic complexity scaling with the number of atoms.

## Physical / mathematical content

Particles are connected when their coordinate differences along all three axes are below `r0` and their Euclidean distance is not greater than `r0`.

## Numerical / algorithmic content

The function sorts each coordinate axis independently, scans nearby entries in sorted order to identify candidate pairs, and intersects the three candidate matrices. It then removes candidate entries whose Euclidean distance exceeds `r0` and converts the result to a sparse matrix.

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

Input validation requires `xyz` to be a real numeric matrix with three columns and `r0` to be a positive real number. The X, Y, and Z scans each build a symmetric logical candidate matrix; the function intersects these matrices, checks Euclidean distances, and returns the result as a sparse logical matrix.