# interfaces/comsol/mesh_vorn.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/interfaces/comsol/mesh_vorn.m`
- Signature: `mesh=mesh_vorn(mesh)`
- Total lines: 67

## Purpose

Voronoi tessellation of a 2D COMSOL mesh. Syntax: mesh=mesh_vorn(mesh)

## Physical / mathematical content

- COMSOL interfaces. These files are mostly data-structure and numerical-geometry utilities for bringing concentration, velocity, and mesh data from finite-element simulations into Spinach transport calculations.

## Numerical / algorithmic content

- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- mesh -Spinach mesh object

## Outputs

- mesh -updated mesh object

## Implementation structure

- Voronoi tessellation of a 2D COMSOL mesh. Syntax:
- mesh=mesh_vorn(mesh)
- mesh -Spinach mesh object
- mesh -updated mesh object
- Check consistency
- Run Voronoi tessellation of the mesh
- Keep only active cells
- Refuse unbounded active cells
- Voronoi cell area calculation
- Add weights to mesh structure
- Find the maximum number of vertices making up the cell
- Consistency enforcement

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `voronoin()`, `cellfun()`, `any()`, `int2str()`, `vor_cell_areas()`, `polyarea()`, `isfield()`.
