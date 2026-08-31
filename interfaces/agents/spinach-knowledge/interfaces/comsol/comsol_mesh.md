# interfaces/comsol/comsol_mesh.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/interfaces/comsol/comsol_mesh.m`
- Signature: `mesh=comsol_mesh(file_name)`
- Total lines: 146

## Purpose

Imports ASCII 2D mesh files produced by COMSOL. Syntax: mesh=comsol_mesh(file_name)

## Physical / mathematical content

- COMSOL interfaces. These files are mostly data-structure and numerical-geometry utilities for bringing concentration, velocity, and mesh data from finite-element simulations into Spinach transport calculations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 30-31: Check consistency; implemented by `grumble(file_name)`.
- Lines 33-34: Open the file; implemented by `fid=fopen(file_name,'r')`.
- Lines 36-37: Point count; implemented by `while true`.
- Lines 46-47: Scroll the file to point coordinates; implemented by `while ~contains(fgetl(fid),'Mesh point coordinates'), end`.
- Lines 49-50: Read vertex coordinates; implemented by `mesh.x=nan(npts,1)`.
- Lines 57-58: Scroll the file to edge specification; implemented by `while ~contains(fgetl(fid),'edg # type name'), end`.
- Lines 60-61: Edge count; implemented by `while true`.
- Lines 70-71: Read edges; implemented by `mesh.idx.edges=nan(nedg,2)`.
- Lines 78-79: Scroll the file to triangle specification; implemented by `while ~contains(fgetl(fid),'tri # type name'), end`.
- Lines 81-82: Triangle count; implemented by `while true`.
- Lines 91-92: Read triangles; implemented by `mesh.idx.triangles=nan(ntri,3)`.
- Lines 100-101: Scroll the file to rectangle specification; implemented by `while ~contains(fgetl(fid),'quad # type name'), end`.
- Lines 103-104: Rectangle count; implemented by `while true`.
- Lines 113-114: Read rectangles; implemented by `mesh.idx.rectangles=nan(nrec,4)`.
- Lines 123-124: Close the file; implemented by `fclose(fid)`.

### Control flow inferred from the code

- Line 37: `while` loop over `true`.
- Line 39: conditional branch on `contains(next_line,'number of mesh points')`.
- Line 47: `while` loop over `~contains(fgetl(fid),'Mesh point coordinates'), end`.
- Line 52: `for` loop over `n=1:npts`.
- Line 58: `while` loop over `~contains(fgetl(fid),'edg # type name'), end`.
- Line 61: `while` loop over `true`.
- Line 63: conditional branch on `contains(next_line,'number of elements')`.
- Line 72: `for` loop over `n=1:nedg`.
- Line 79: `while` loop over `~contains(fgetl(fid),'tri # type name'), end`.
- Line 82: `while` loop over `true`.
- Line 84: conditional branch on `contains(next_line,'number of elements')`.
- Line 93: `for` loop over `n=1:ntri`.
- Line 101: `while` loop over `~contains(fgetl(fid),'quad # type name'), end`.
- Line 104: `while` loop over `true`.

### Key state/data transformations

- Lines 34: computes `fid` using `fid=fopen(file_name,'r')`.
- Lines 38: computes `next_line` using `next_line=fgetl(fid)`.
- Lines 40: computes `npts` using `npts=textscan(next_line,'%f')`.
- Lines 50: computes `mesh.x` using `mesh.x=nan(npts,1)`.
- Lines 51: computes `mesh.y` using `mesh.y=nan(npts,1)`.
- Lines 53: computes `XY` using `XY=textscan(fgetl(fid),'%f %f')`.
- Lines 54: computes `mesh.x(n)` using `mesh.x(n)=XY{1}; mesh.y(n)=XY{2}`.
- Lines 64: computes `nedg` using `nedg=textscan(next_line,'%f')`.
- Lines 71: computes `mesh.idx.edges` using `mesh.idx.edges=nan(nedg,2)`.
- Lines 73: computes `ES` using `ES=textscan(fgetl(fid),'%f %f')`.
- Lines 74: computes `mesh.idx.edges(n,1)` using `mesh.idx.edges(n,1)=ES{1}+1`.
- Lines 75: computes `mesh.idx.edges(n,2)` using `mesh.idx.edges(n,2)=ES{2}+1`.
- Lines 85: computes `ntri` using `ntri=textscan(next_line,'%f')`.
- Lines 92: computes `mesh.idx.triangles` using `mesh.idx.triangles=nan(ntri,3)`.
- Lines 94: computes `TS` using `TS=textscan(fgetl(fid),'%f %f %f')`.
- Lines 95: computes `mesh.idx.triangles(n,1)` using `mesh.idx.triangles(n,1)=TS{1}+1`.
- Lines 96: computes `mesh.idx.triangles(n,2)` using `mesh.idx.triangles(n,2)=TS{2}+1`.
- Lines 97: computes `mesh.idx.triangles(n,3)` using `mesh.idx.triangles(n,3)=TS{3}+1`.

### Local helper functions

- Line 129: `grumble()` — `function grumble(file_name)`. Although Ray had a very much hands-off approach to men- toring, he occasionally would wander into the spectro-
  - Representative operation: `if ~ischar(file_name)`.
  - Representative operation: `error('file_name must be a character string.')`.

## Parameters / inputs

- file_name -a character string

## Outputs

- mesh.x, mesh.y -column vectors with vertex
- coordinates
- mesh.idx.edges -two-column array of integers
- containing edge index
- mesh.idx.triangles -three-column array of integers
- containing triangle index
- mesh.idx.rectangles -four-column array of integers
- containing rectangle index

## Implementation structure

- Imports ASCII 2D mesh files produced by COMSOL. Syntax:
- mesh=comsol_mesh(file_name)
- file_name -a character string
- mesh.x, mesh.y -column vectors with vertex
- coordinates
- mesh.idx.edges -two-column array of integers
- containing edge index
- mesh.idx.triangles -three-column array of integers
- containing triangle index
- mesh.idx.rectangles -four-column array of integers
- containing rectangle index
- Check consistency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `fopen()`, `fgetl()`, `contains()`, `textscan()`, `num2str()`, `nan()`, `fclose()`, `ischar()`.
