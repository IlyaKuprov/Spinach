# interfaces/comsol/mesh_crop.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/interfaces/comsol/mesh_crop.m`
- Signature: `mesh=mesh_crop(mesh,ranges)`
- Total lines: 109

## Purpose

2D microfluidic mesh cropping. Updates the mesh object to remove anything outside the user-specified vertex coordi- nate ranges. Syntax: mesh=mesh_crop(mesh,ranges)

## Physical / mathematical content

- COMSOL interfaces. These files are mostly data-structure and numerical-geometry utilities for bringing concentration, velocity, and mesh data from finite-element simulations into Spinach transport calculations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- mesh -Spinach mesh object
- ranges -{[xmin xmax],[ymin ymax]}

## Outputs

- mesh -updated mesh object

## Implementation structure

- 2D microfluidic mesh cropping. Updates the mesh object to
- remove anything outside the user-specified vertex coordi-
- nate ranges. Syntax:
- mesh=mesh_crop(mesh,ranges)
- mesh -Spinach mesh object
- ranges -{[xmin xmax],[ymin ymax]}
- mesh -updated mesh object
- Check consistency
- Remove tessellation and preplot
- Find vertices in the user-specified range
- Find edges in the user-specified range
- Re-index edges with updated vertices

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `isfield()`, `rmfield()`, `ismember()`, `iscell()`.
