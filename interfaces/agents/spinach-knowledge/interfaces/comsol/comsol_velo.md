# interfaces/comsol/comsol_velo.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/interfaces/comsol/comsol_velo.m`
- Signature: `mesh=comsol_velo(mesh,file_name)`
- Total lines: 74

## Purpose

Imports ASCII 2D flow velocity files produced by COMSOL. Syntax: mesh=comsol_velo(mesh,file_name)

## Physical / mathematical content

- COMSOL interfaces. These files are mostly data-structure and numerical-geometry utilities for bringing concentration, velocity, and mesh data from finite-element simulations into Spinach transport calculations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- mesh -mesh object produced by comsol_mesh()
- file_name -a character string

## Outputs

- the following fields are added to the mesh object
- mesh.u, mesh.v -column vectors with velocities
- at each vertex of the mesh

## Implementation structure

- Imports ASCII 2D flow velocity files produced by COMSOL. Syntax:
- mesh=comsol_velo(mesh,file_name)
- mesh -mesh object produced by comsol_mesh()
- file_name -a character string
- the following fields are added to the mesh object
- mesh.u, mesh.v -column vectors with velocities
- at each vertex of the mesh
- Check consistency
- Open the file
- Velocity readout count
- Parse velocity readouts
- Close the file

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `fopen()`, `fgetl()`, `contains()`, `textscan()`, `num2str()`, `nan()`, `fclose()`, `ischar()`.
