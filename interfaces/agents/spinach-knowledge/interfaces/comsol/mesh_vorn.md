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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 20-21: Check consistency; implemented by `grumble(mesh)`.
- Lines 23-25: Run Voronoi tessellation of the mesh; implemented by `[V,C]=voronoin([mesh.x mesh.y])`.
- Lines 29-30: Keep only active cells; implemented by `mesh.vor.cells=mesh.vor.cells(mesh.idx.active)`.
- Lines 33-34: Refuse unbounded active cells; implemented by `unbounded=cellfun(@(x)any(x==1),mesh.vor.cells)`.
- Lines 40-41: Voronoi cell area calculation; implemented by `vor_cell_areas=zeros(mesh.vor.ncells,1)`.
- Lines 48-49: Add weights to mesh structure; implemented by `mesh.vor.weights=vor_cell_areas`.
- Lines 51-52: Find the maximum number of vertices making up the cell; implemented by `mesh.vor.max_cell_size=max(cellfun(@numel,mesh.vor.cells))`.

### Control flow inferred from the code

- Line 35: conditional branch on `any(unbounded)`.
- Line 42: `for` loop over `n=1:mesh.vor.ncells`.

### Key state/data transformations

- Lines 24-25: computes `[V,C]` using `[V,C]=voronoin([mesh.x mesh.y])`.
- Lines 26: computes `mesh.vor.vertices` using `mesh.vor.vertices=V`.
- Lines 27: computes `mesh.vor.cells` using `mesh.vor.cells=C`.
- Lines 31: computes `mesh.vor.ncells` using `mesh.vor.ncells=numel(mesh.vor.cells)`.
- Lines 41: computes `vor_cell_areas` using `vor_cell_areas=zeros(mesh.vor.ncells,1)`.
- Lines 43: computes `vor_xcoords` using `vor_xcoords=mesh.vor.vertices(mesh.vor.cells{n},1)`.
- Lines 44: computes `vor_ycoords` using `vor_ycoords=mesh.vor.vertices(mesh.vor.cells{n},2)`.
- Lines 45: computes `vor_cell_areas(n)` using `vor_cell_areas(n)=polyarea(vor_xcoords,vor_ycoords)`.
- Lines 49: computes `mesh.vor.weights` using `mesh.vor.weights=vor_cell_areas`.
- Lines 52: computes `mesh.vor.max_cell_size` using `mesh.vor.max_cell_size=max(cellfun(@numel,mesh.vor.cells))`.

### Local helper functions

- Line 57: `grumble()` — `function grumble(mesh)`. A man's rights rest in three boxes: the ballot box, the jury box, and the cartridge box.
  - Representative operation: `if (~isfield(mesh,'idx'))||(~isfield(mesh.idx,'active'))`.
  - Representative operation: `error('active vertex index is missing from mesh structure.')`.

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
