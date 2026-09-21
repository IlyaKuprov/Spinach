# interfaces/comsol/mesh_inact.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/interfaces/comsol/mesh_inact.m`
- Signature: `mesh=mesh_inact(mesh,vertex_list)`
- Total lines: 63

## Purpose

Marks 2D microfluidic mesh vertices as inactive in hydrodyna- mic and diffusive transport processes. Syntax: mesh=mesh_inact(mesh,vertex_list)

## Physical / mathematical content

- COMSOL interfaces. These files are mostly data-structure and numerical-geometry utilities for bringing concentration, velocity, and mesh data from finite-element simulations into Spinach transport calculations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- mesh -Spinach mesh object
- vertex_list -row vector of integers specifying
- the vertices to be inactivated

## Outputs

- mesh -updated mesh object

## Implementation structure

- Marks 2D microfluidic mesh vertices as inactive in hydrodyna-
- mic and diffusive transport processes. Syntax:
- mesh=mesh_inact(mesh,vertex_list)
- mesh -Spinach mesh object
- vertex_list -row vector of integers specifying
- the vertices to be inactivated
- mesh -updated mesh object
- Check consistency
- Update the active vertex list
- Zero out velocities and concentrations, if present
- Consistency enforcement
- The basic principle of the new education is to be that dunces and

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `setdiff()`, `isfield()`, `isrow()`, `any()`.
