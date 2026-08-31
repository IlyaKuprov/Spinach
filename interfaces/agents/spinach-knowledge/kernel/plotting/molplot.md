# kernel/plotting/molplot.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/plotting/molplot.m`
- Signature: `molplot(xyz,conmatrix)`
- Total lines: 69

## Purpose

Plots a stick representation of a molecule from Cartesian coordinates supplied. Syntax: molplot(xyz,conmatrix)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 26-27: Check consistency; implemented by `grumble(xyz,conmatrix)`.
- Lines 29-30: Get the connectivity matrix; implemented by `if isempty(conmatrix)`.
- Lines 34-35: Prepare coordinate arrays; implemented by `nbonds=nnz(conmatrix)`.
- Lines 46-47: Draw the molecule; implemented by `plot3(X,Y,Z,'Color',[0.5 0.5 0.5],'LineWidth',1.5)`.

### Control flow inferred from the code

- Line 30: conditional branch on `isempty(conmatrix)`.
- Line 40: `for` loop over `n=1:nbonds`.

### Key state/data transformations

- Lines 31: computes `conmatrix` using `conmatrix=conmat(xyz,1.6)`.
- Lines 35: computes `nbonds` using `nbonds=nnz(conmatrix)`.
- Lines 36: computes `X` using `X=zeros(1,3*nbonds)`.
- Lines 37: computes `Y` using `Y=zeros(1,3*nbonds)`.
- Lines 38: computes `Z` using `Z=zeros(1,3*nbonds)`.
- Lines 39: computes `[rows,cols]` using `[rows,cols]=find(conmatrix)`.
- Lines 41: computes `X((3*(n-1)+1):(3*n))` using `X((3*(n-1)+1):(3*n))=[xyz(rows(n),1) xyz(cols(n),1) NaN]`.
- Lines 42: computes `Y((3*(n-1)+1):(3*n))` using `Y((3*(n-1)+1):(3*n))=[xyz(rows(n),2) xyz(cols(n),2) NaN]`.
- Lines 43: computes `Z((3*(n-1)+1):(3*n))` using `Z((3*(n-1)+1):(3*n))=[xyz(rows(n),3) xyz(cols(n),3) NaN]`.

### Local helper functions

- Line 52: `grumble()` — `function grumble(xyz,conmatrix)`.
  - Representative operation: `if (~isnumeric(xyz))||(~isreal(xyz))||(size(xyz,2)~=3)`.
  - Representative operation: `error('xyz must be an Nx3 real matrix.')`.

## Parameters / inputs

- xyz -Cartesian coordinates, as Nx3 matrix, in
- Angstroms
- conmatrix -NxN connectivity matrix indicating chemical
- bonds that should be drawn as sticks. If an
- empty vector is supplied, 1.6 Angstrom cut-
- off distance is used

## Outputs

- this function creates a figure

## Implementation structure

- Plots a stick representation of a molecule from Cartesian coordinates
- supplied. Syntax:
- molplot(xyz,conmatrix)
- xyz -Cartesian coordinates, as Nx3 matrix, in
- Angstroms
- conmatrix -NxN connectivity matrix indicating chemical
- bonds that should be drawn as sticks. If an
- empty vector is supplied, 1.6 Angstrom cut-
- off distance is used
- this function creates a figure
- Check consistency
- Get the connectivity matrix

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `conmat()`, `nnz()`, `xyz()`, `rows()`, `cols()`, `plot3()`.
