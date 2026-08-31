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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 24-25: Check consistency; implemented by `grumble(file_name)`.
- Lines 27-28: Open the file; implemented by `fid=fopen(file_name,'r')`.
- Lines 30-31: Velocity readout count; implemented by `while true`.
- Lines 41-42: Parse velocity readouts; implemented by `X=nan(nvel,1); Y=nan(nvel,1)`.
- Lines 50-51: Close the file; implemented by `fclose(fid)`.
- Lines 53-54: Check that vertex locations are the same; implemented by `if (norm(mesh.x-X,1)>1e-6)||(norm(mesh.y-Y,1)>1e-6)`.
- Lines 58-59: Store velocities; implemented by `mesh.u=U; mesh.v=V`.

### Control flow inferred from the code

- Line 31: `while` loop over `true`.
- Line 33: conditional branch on `contains(next_line,'% Nodes:')`.
- Line 44: `for` loop over `n=1:nvel`.
- Line 54: conditional branch on `(norm(mesh.x-X,1)>1e-6)||(norm(mesh.y-Y,1)>1e-6)`.

### Key state/data transformations

- Lines 28: computes `fid` using `fid=fopen(file_name,'r')`.
- Lines 32: computes `next_line` using `next_line=fgetl(fid)`.
- Lines 34: computes `nvel` using `nvel=textscan(next_line,'%s %s %f')`.
- Lines 42: computes `X` using `X=nan(nvel,1); Y=nan(nvel,1)`.
- Lines 43: computes `U` using `U=nan(nvel,1); V=nan(nvel,1)`.
- Lines 45: computes `VS` using `VS=textscan(fgetl(fid),'%f %f %f %f %f')`.
- Lines 46: computes `X(n)` using `X(n)=VS{1}; Y(n)=VS{2}`.
- Lines 47: computes `U(n)` using `U(n)=VS{4}; V(n)=VS{5}`.
- Lines 59: computes `mesh.u` using `mesh.u=U; mesh.v=V`.

### Local helper functions

- Line 64: `grumble()` — `function grumble(file_name)`. Words have no power to impress the mind without the exquisite horror of their reality.
  - Representative operation: `if ~ischar(file_name)`.
  - Representative operation: `error('file_name must be a character string.')`.

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
