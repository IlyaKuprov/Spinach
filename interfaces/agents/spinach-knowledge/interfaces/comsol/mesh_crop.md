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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 24-25: Check consistency; implemented by `grumble(mesh,ranges)`.
- Lines 27-28: Remove tessellation and preplot; implemented by `if isfield(mesh,'vor')`.
- Lines 37-39: Find vertices in the user-specified range; implemented by `vertices_in_range=(mesh.x>=ranges{1}(1))&(mesh.x<=ranges{1}(2))& (mesh.y>=ranges{2}(1))&(mesh.y<=ranges{2}(2))`.
- Lines 42-44: Find edges in the user-specified range; implemented by `edges_in_range=ismember(mesh.idx.edges(:,1),vertices_in_range)& ismember(mesh.idx.edges(:,2),vertices_in_range)`.
- Lines 47-48: Re-index edges with updated vertices; implemented by `[~,mesh.idx.edges(:,1)]=ismember(mesh.idx.edges(:,1),vertices_in_range)`.
- Lines 51-54: Find triangles in the user-specified range; implemented by `tris_in_range=ismember(mesh.idx.triangles(:,1),vertices_in_range)& ismember(mesh.idx.triangles(:,2),vertices_in_range)& ismember(mesh.idx.triangles(:,3),vertices_in_rang…`.
- Lines 57-58: Re-index trianges with updated vertices; implemented by `[~,mesh.idx.triangles(:,1)]=ismember(mesh.idx.triangles(:,1),vertices_in_range)`.
- Lines 62-66: Find rectangles in the user-specified range; implemented by `rect_in_range=ismember(mesh.idx.rectangles(:,1),vertices_in_range)& ismember(mesh.idx.rectangles(:,2),vertices_in_range)& ismember(mesh.idx.rectangles(:,3),vertices_in_r…`.
- Lines 69-70: Re-index rectangles with updated vertices; implemented by `[~,mesh.idx.rectangles(:,1)]=ismember(mesh.idx.rectangles(:,1),vertices_in_range)`.
- Lines 75-76: Crop coordinates, velocities and concentrations; implemented by `mesh.x=mesh.x(vertices_in_range)`.
- Lines 82-83: First guess active vertices; implemented by `if isfield(mesh.idx,'active')`.

### Control flow inferred from the code

- Line 28: conditional branch on `isfield(mesh,'vor')`.
- Line 32: conditional branch on `isfield(mesh,'plot')`.
- Line 78: conditional branch on `isfield(mesh,'u'), mesh.u=mesh.u(vertices_in_range); end`.
- Line 79: conditional branch on `isfield(mesh,'v'), mesh.v=mesh.v(vertices_in_range); end`.
- Line 80: conditional branch on `isfield(mesh,'c'), mesh.c=mesh.c(vertices_in_range,:); end`.
- Line 83: conditional branch on `isfield(mesh.idx,'active')`.

### Key state/data transformations

- Lines 29: computes `mesh` using `mesh=rmfield(mesh,'vor')`.
- Lines 40: computes `vertices_in_range` using `vertices_in_range=find(vertices_in_range)`.
- Lines 43-44: computes `edges_in_range` using `edges_in_range=ismember(mesh.idx.edges(:,1),vertices_in_range)& ismember(mesh.idx.edges(:,2),vertices_in_range)`.
- Lines 45: computes `mesh.idx.edges` using `mesh.idx.edges=mesh.idx.edges(edges_in_range,:)`.
- Lines 48: computes `[~,mesh.idx.edges(:,1)]` using `[~,mesh.idx.edges(:,1)]=ismember(mesh.idx.edges(:,1),vertices_in_range)`.
- Lines 49: computes `[~,mesh.idx.edges(:,2)]` using `[~,mesh.idx.edges(:,2)]=ismember(mesh.idx.edges(:,2),vertices_in_range)`.
- Lines 52-54: computes `tris_in_range` using `tris_in_range=ismember(mesh.idx.triangles(:,1),vertices_in_range)& ismember(mesh.idx.triangles(:,2),vertices_in_range)& ismember(mesh.idx.triangles(:,3),vertices_in_rang…`.
- Lines 55: computes `mesh.idx.triangles` using `mesh.idx.triangles=mesh.idx.triangles(tris_in_range,:)`.
- Lines 58: computes `[~,mesh.idx.triangles(:,1)]` using `[~,mesh.idx.triangles(:,1)]=ismember(mesh.idx.triangles(:,1),vertices_in_range)`.
- Lines 59: computes `[~,mesh.idx.triangles(:,2)]` using `[~,mesh.idx.triangles(:,2)]=ismember(mesh.idx.triangles(:,2),vertices_in_range)`.
- Lines 60: computes `[~,mesh.idx.triangles(:,3)]` using `[~,mesh.idx.triangles(:,3)]=ismember(mesh.idx.triangles(:,3),vertices_in_range)`.
- Lines 63-66: computes `rect_in_range` using `rect_in_range=ismember(mesh.idx.rectangles(:,1),vertices_in_range)& ismember(mesh.idx.rectangles(:,2),vertices_in_range)& ismember(mesh.idx.rectangles(:,3),vertices_in_r…`.
- Lines 67: computes `mesh.idx.rectangles` using `mesh.idx.rectangles=mesh.idx.rectangles(rect_in_range,:)`.
- Lines 70: computes `[~,mesh.idx.rectangles(:,1)]` using `[~,mesh.idx.rectangles(:,1)]=ismember(mesh.idx.rectangles(:,1),vertices_in_range)`.
- Lines 71: computes `[~,mesh.idx.rectangles(:,2)]` using `[~,mesh.idx.rectangles(:,2)]=ismember(mesh.idx.rectangles(:,2),vertices_in_range)`.
- Lines 72: computes `[~,mesh.idx.rectangles(:,3)]` using `[~,mesh.idx.rectangles(:,3)]=ismember(mesh.idx.rectangles(:,3),vertices_in_range)`.
- Lines 73: computes `[~,mesh.idx.rectangles(:,4)]` using `[~,mesh.idx.rectangles(:,4)]=ismember(mesh.idx.rectangles(:,4),vertices_in_range)`.
- Lines 76: computes `mesh.x` using `mesh.x=mesh.x(vertices_in_range)`.

### Local helper functions

- Line 91: `grumble()` — `function grumble(mesh,ranges)`.
  - Representative operation: `if ~isfield(mesh,'idx')`.
  - Representative operation: `error('vertex index is missing from the mesh structure.')`.

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
