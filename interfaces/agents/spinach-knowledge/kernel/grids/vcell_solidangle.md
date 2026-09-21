# kernel/grids/vcell_solidangle.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/grids/vcell_solidangle.m`
- Signature: `S=vcell_solidangle(P,K,xyz)`
- Total lines: 85

## Purpose

Solid angle of a spherical Voronoi cell. Syntax: s=vcell_solidangle(P,K,xyz)

## Physical / mathematical content

- Quadrature and geometry utilities. These files generate spherical/SO(3) grids, Voronoi weights, and adaptive integration tools for orientation averaging.

## Numerical / algorithmic content

- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- P -(3 x m) array with coordinates of the
- vertices of the Voronoi cell
- K -(n x 1) cell, each K{j} contains the
- indices of the Voronoi cell
- xyz -optional (3 x n) knot points to guide
- vcell_solidangle to compute the solid
- angle of the "right" cell containing
- the node (and not the complement cell)

## Outputs

- S -the solid angle of the Voronoi cell

## Implementation structure

- Solid angle of a spherical Voronoi cell. Syntax:
- s=vcell_solidangle(P,K,xyz)
- P -(3 x m) array with coordinates of the
- vertices of the Voronoi cell
- K -(n x 1) cell, each K{j} contains the
- indices of the Voronoi cell
- xyz -optional (3 x n) knot points to guide
- vcell_solidangle to compute the solid
- angle of the "right" cell containing
- the node (and not the complement cell)
- S -the solid angle of the Voronoi cell
- Check consistency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `cellfun()`, `one_vcell_solidangle()`, `arrayfun()`, `xyz()`, `any()`, `iscell()`.
