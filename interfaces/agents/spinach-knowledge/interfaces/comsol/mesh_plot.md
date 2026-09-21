# interfaces/comsol/mesh_plot.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/interfaces/comsol/mesh_plot.m`
- Signature: `mesh_plot(spin_system,qscale,nodelabels)`
- Total lines: 101

## Purpose

2D microfluidic mesh plotting function. Draws the mesh, its Vo- ronoi tessellation, and a quiver plot of velocities. Syntax: mesh_plot(spin_system,qscale,nodelabels)

## Physical / mathematical content

- COMSOL interfaces. These files are mostly data-structure and numerical-geometry utilities for bringing concentration, velocity, and mesh data from finite-element simulations into Spinach transport calculations.

## Numerical / algorithmic content

- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- spin_system -Spinach spin system object containing
- mesh information
- qscale -scaling multiplier for the quiver plot
- of flow velocities, zero turns veloci-
- ty plotting off
- nodelabels -1 causes vertex numbers to be displa-
- yed, 0 turns that off

## Outputs

- the function creates a figure

## Implementation structure

- 2D microfluidic mesh plotting function. Draws the mesh, its Vo-
- ronoi tessellation, and a quiver plot of velocities. Syntax:
- mesh_plot(spin_system,qscale,nodelabels)
- spin_system -Spinach spin system object containing
- mesh information
- qscale -scaling multiplier for the quiver plot
- of flow velocities, zero turns veloci-
- ty plotting off
- nodelabels -1 causes vertex numbers to be displa-
- yed, 0 turns that off
- the function creates a figure
- Check consistency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `patch()`, `kxlabel()`, `kylabel()`, `quiver()`, `arrayfun()`, `text()`, `set()`, `isfield()`, `isscalar()`, `ismember()`.
