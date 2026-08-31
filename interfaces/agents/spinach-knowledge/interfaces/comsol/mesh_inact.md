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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 23-24: Check consistency; implemented by `grumble(mesh,vertex_list)`.
- Lines 26-27: Update the active vertex list; implemented by `mesh.idx.active=setdiff(mesh.idx.active,vertex_list)`.
- Lines 29-30: Zero out velocities and concentrations, if present; implemented by `vertex_list=setdiff(1:numel(mesh.x),mesh.idx.active)`.

### Control flow inferred from the code

- Line 31: conditional branch on `isfield(mesh,'u'), mesh.u(vertex_list)=0; end`.
- Line 32: conditional branch on `isfield(mesh,'v'), mesh.v(vertex_list)=0; end`.
- Line 33: conditional branch on `isfield(mesh,'c'), mesh.c(vertex_list,:)=0; end`.

### Key state/data transformations

- Lines 27: computes `mesh.idx.active` using `mesh.idx.active=setdiff(mesh.idx.active,vertex_list)`.
- Lines 30: computes `vertex_list` using `vertex_list=setdiff(1:numel(mesh.x),mesh.idx.active)`.

### Local helper functions

- Line 38: `grumble()` — `function grumble(mesh,vertex_list)`.
  - Representative operation: `if ~isfield(mesh,'idx')`.
  - Representative operation: `error('indexing information is missing from mesh structure.')`.

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
