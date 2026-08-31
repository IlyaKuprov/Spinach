# interfaces/comsol/mesh_preplot.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/interfaces/comsol/mesh_preplot.m`
- Signature: `mesh=mesh_preplot(mesh)`
- Total lines: 86

## Purpose

Mesh preprocessing for drawing. Creates edge, triangle, and rectangle data structures needed for fast plotting later.

## Physical / mathematical content

- COMSOL interfaces. These files are mostly data-structure and numerical-geometry utilities for bringing concentration, velocity, and mesh data from finite-element simulations into Spinach transport calculations.

## Numerical / algorithmic content

- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 20-21: Check consistency; implemented by `grumble(mesh)`.
- Lines 23-24: Prepare edge array for plotting; implemented by `nlines=size(mesh.idx.edges,1); A=zeros(1,3*nlines); B=zeros(1,3*nlines)`.
- Lines 31-32: Prepare triangle array for plotting; implemented by `nlines=size(mesh.idx.triangles,1); A=zeros(1,5*nlines); B=zeros(1,5*nlines)`.
- Lines 41-42: Prepare rectangle array for plotting; implemented by `nlines=size(mesh.idx.rectangles,1); A=zeros(1,6*nlines); B=zeros(1,6*nlines)`.
- Lines 53-54: Prepare Voronoi cell array for plotting; implemented by `cell_sizes=cellfun(@numel,mesh.vor.cells)`.

### Control flow inferred from the code

- Line 25: `for` loop over `n=1:nlines`.
- Line 33: `for` loop over `n=1:nlines`.
- Line 43: `for` loop over `n=1:nlines`.
- Line 57: `for` loop over `n=1:numel(mesh.vor.cells)`.
- Line 66: conditional branch on `isempty(cell_sizes)`.

### Key state/data transformations

- Lines 24: computes `nlines` using `nlines=size(mesh.idx.edges,1); A=zeros(1,3*nlines); B=zeros(1,3*nlines)`.
- Lines 26: computes `A((3*(n-1)+1):(3*n))` using `A((3*(n-1)+1):(3*n))=[mesh.x(mesh.idx.edges(n,1)) mesh.x(mesh.idx.edges(n,2)) NaN]`.
- Lines 27: computes `B((3*(n-1)+1):(3*n))` using `B((3*(n-1)+1):(3*n))=[mesh.y(mesh.idx.edges(n,1)) mesh.y(mesh.idx.edges(n,2)) NaN]`.
- Lines 29: computes `mesh.plot.edg_a` using `mesh.plot.edg_a=A; mesh.plot.edg_b=B`.
- Lines 34-35: computes `A((5*(n-1)+1):(5*n))` using `A((5*(n-1)+1):(5*n))=[mesh.x(mesh.idx.triangles(n,1)) mesh.x(mesh.idx.triangles(n,2)) mesh.x(mesh.idx.triangles(n,3)) mesh.x(mesh.idx.triangles(n,1)) NaN]`.
- Lines 36-37: computes `B((5*(n-1)+1):(5*n))` using `B((5*(n-1)+1):(5*n))=[mesh.y(mesh.idx.triangles(n,1)) mesh.y(mesh.idx.triangles(n,2)) mesh.y(mesh.idx.triangles(n,3)) mesh.y(mesh.idx.triangles(n,1)) NaN]`.
- Lines 39: computes `mesh.plot.tri_a` using `mesh.plot.tri_a=A; mesh.plot.tri_b=B`.
- Lines 44-46: computes `A((6*(n-1)+1):(6*n))` using `A((6*(n-1)+1):(6*n))=[mesh.x(mesh.idx.rectangles(n,1)) mesh.x(mesh.idx.rectangles(n,2)) mesh.x(mesh.idx.rectangles(n,4)) mesh.x(mesh.idx.rectangles(n,3)) mesh.x(mesh.idx…`.
- Lines 47-49: computes `B((6*(n-1)+1):(6*n))` using `B((6*(n-1)+1):(6*n))=[mesh.y(mesh.idx.rectangles(n,1)) mesh.y(mesh.idx.rectangles(n,2)) mesh.y(mesh.idx.rectangles(n,4)) mesh.y(mesh.idx.rectangles(n,3)) mesh.y(mesh.idx…`.
- Lines 51: computes `mesh.plot.rec_a` using `mesh.plot.rec_a=A; mesh.plot.rec_b=B`.
- Lines 54: computes `cell_sizes` using `cell_sizes=cellfun(@numel,mesh.vor.cells)`.
- Lines 55: computes `array_size` using `array_size=sum(cell_sizes)+2*numel(mesh.vor.cells)`.
- Lines 56: computes `A` using `A=nan(array_size,1); B=nan(array_size,1); array_offset=0`.
- Lines 58: computes `cell_idx` using `cell_idx=mesh.vor.cells{n}; nvert=cell_sizes(n)`.
- Lines 59: computes `boundary_range` using `boundary_range=array_offset+(1:(nvert+1))`.
- Lines 60-61: computes `A(boundary_range)` using `A(boundary_range)=[mesh.vor.vertices(cell_idx,1); mesh.vor.vertices(cell_idx(1),1)]`.
- Lines 62-63: computes `B(boundary_range)` using `B(boundary_range)=[mesh.vor.vertices(cell_idx,2); mesh.vor.vertices(cell_idx(1),2)]`.
- Lines 64: computes `array_offset` using `array_offset=array_offset+nvert+2`.

### Local helper functions

- Line 74: `grumble()` — `function grumble(mesh)`. Life is a tragedy for those who feel, and a comedy for those who think.
  - Representative operation: `if ~isfield(mesh,'idx')`.
  - Representative operation: `error('vertex index is missing from spin_system.mesh structure.')`.

## Syntax

```matlab
mesh=mesh_preplot(mesh)
```

## Parameters / inputs

- mesh -Spinach mesh object

## Outputs

- mesh -updated mesh object

## Implementation structure

- Mesh preprocessing for drawing. Creates edge, triangle, and
- rectangle data structures needed for fast plotting later.
- mesh=mesh_preplot(mesh)
- mesh -Spinach mesh object
- mesh -updated mesh object
- Check consistency
- Prepare edge array for plotting
- Prepare triangle array for plotting
- Prepare rectangle array for plotting
- Prepare Voronoi cell array for plotting
- Consistency enforcement
- Life is a tragedy for those who feel,

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `cellfun()`, `nan()`, `cell_sizes()`, `cell_idx()`, `isfield()`.
