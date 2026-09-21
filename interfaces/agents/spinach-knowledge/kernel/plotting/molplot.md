# kernel/plotting/molplot.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/plotting/molplot.m`
- Signature: `molplot(xyz,conmatrix)`
- Total lines: 69

## Purpose

Plots a stick representation of a molecule from Cartesian coordinates supplied. Syntax: molplot(xyz,conmatrix)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `conmat()`, `nnz()`, `xyz()`, `rows()`, `cols()`, `plot3()`.
