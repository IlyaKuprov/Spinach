# interfaces/comsol/comsol_import.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/interfaces/comsol/comsol_import.m`
- Signature: `mesh=comsol_import(comsol)`
- Total lines: 100

## Purpose

COMSOL 2D mesh data import, cropping and preprocessing for Spinach. Syntax: mesh=comsol_import(comsol)

## Physical / mathematical content

- COMSOL interfaces. These files are mostly data-structure and numerical-geometry utilities for bringing concentration, velocity, and mesh data from finite-element simulations into Spinach transport calculations.

## Numerical / algorithmic content

- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 46-47: Check consistency; implemented by `grumble(comsol)`.
- Lines 49-50: Import the mesh; implemented by `mesh=comsol_mesh(comsol.mesh_file)`.
- Lines 52-53: Import the velocities; implemented by `mesh=comsol_velo(mesh,comsol.velo_file)`.
- Lines 55-56: Crop to the region of interest; implemented by `mesh=mesh_crop(mesh,comsol.crop)`.
- Lines 58-59: Inactivate user-specified vertices; implemented by `mesh=mesh_inact(mesh,comsol.inactivate)`.
- Lines 61-62: Run Voronoi tessellation; implemented by `mesh=mesh_vorn(mesh)`.
- Lines 64-65: Run graphical output preprocessing; implemented by `mesh=mesh_preplot(mesh)`.

### Key state/data transformations

- Lines 50: computes `mesh` using `mesh=comsol_mesh(comsol.mesh_file)`.

### Local helper functions

- Line 70: `grumble()` — `function grumble(comsol)`.
  - Representative operation: `if ~isstruct(comsol), error('comsol must be a structure.'); end`.
  - Representative operation: `if (~isfield(comsol,'mesh_file'))||(~ischar(comsol.mesh_file))`.

## Parameters / inputs

- comsol.mesh_file -name of an ASCII file with
- vertex coordinates and edge
- index produced by COMSOL
- comsol.velo_file -name of an ASCII file with
- vertex-centred flow veloci-
- ties produced by COMSOL
- comsol.crop -{[xmin xmax],[ymin ymax]}
- region of the mesh to retain
- comsol.inactivate -a row vector with mesh vertex
- indices to deactivate

## Outputs

- mesh – Spinach mesh object with
- ▸ geometry (.x, .y, .idx)
- ▸ flow velocities (.u, .v)
- ▸ Voronoi tessellation (.vor)
- ▸ fast-plot auxiliaries (.plot)
- Notes: internally this routine is just a convenience wrapper
- that calls the following functions
- ▸ comsol_mesh() – mesh import
- ▸ comsol_velo() – velocity import
- ▸ mesh_crop() – region trimming
- ▸ mesh_vorn() – Voronoi tessellation
- ▸ mesh_preplot() – plotting accelerators

## Implementation structure

- COMSOL 2D mesh data import, cropping and preprocessing for
- Spinach. Syntax:
- mesh=comsol_import(comsol)
- comsol.mesh_file -name of an ASCII file with
- vertex coordinates and edge
- index produced by COMSOL
- comsol.velo_file -name of an ASCII file with
- vertex-centred flow veloci-
- ties produced by COMSOL
- comsol.crop -{[xmin xmax],[ymin ymax]}
- region of the mesh to retain
- comsol.inactivate -a row vector with mesh vertex

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `comsol_mesh()`, `comsol_velo()`, `mesh_crop()`, `mesh_inact()`, `mesh_vorn()`, `mesh_preplot()`, `isstruct()`, `isfield()`, `ischar()`, `iscell()`, `isrow()`, `any()`.
