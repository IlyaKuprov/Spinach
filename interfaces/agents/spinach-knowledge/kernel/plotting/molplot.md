# kernel/plotting/molplot.m

- Signature: `molplot(xyz,conmatrix)`

## Purpose

Plots a stick representation of a molecule from Cartesian coordinates supplied. Syntax: molplot(xyz,conmatrix)

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- xyz -Cartesian coordinates, as Nx3 matrix, in
- Angstroms
- conmatrix -NxN connectivity matrix indicating chemical
- bonds that should be drawn as sticks. If an
- empty vector is supplied, 1.6 Angstrom cut-
- off distance is used

## Outputs

- this function creates a figure

## Implementation structure

- Plots a stick representation of a molecule from Cartesian coordinates
- supplied. Syntax:
- molplot(xyz,conmatrix)
- xyz -Cartesian coordinates, as Nx3 matrix, in
- Angstroms
- conmatrix -NxN connectivity matrix indicating chemical
- bonds that should be drawn as sticks. If an
- empty vector is supplied, 1.6 Angstrom cut-
- off distance is used
- this function creates a figure
- Check consistency
- Get the connectivity matrix
