# kernel/plotting/molplot.m

- Signature: `molplot(xyz,conmatrix)`

## Purpose

Plots a stick representation of a molecule from Cartesian coordinates supplied. Syntax: molplot(xyz,conmatrix)

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- `xyz` — Cartesian coordinates as an N-by-3 matrix, in Angstroms.
- `conmatrix` — N-by-N connectivity matrix indicating bonds to draw as sticks. If empty, connectivity is computed with `conmat(xyz,1.6)`, using a 1.6 Angstrom cutoff.

## Outputs

- this function creates a figure

## Implementation structure

- Validates that `xyz` is an N-by-3 coordinate array and that `conmatrix`, when supplied, is a logical square matrix with one row per atom.
- If `conmatrix` is empty, obtains it from `conmat(xyz,1.6)`.
- Builds NaN-separated coordinate arrays for each connected atom pair and draws the sticks with `plot3` in grey.

[Source reference](https://spindynamics.org/wiki/index.php?title=molplot.m)
