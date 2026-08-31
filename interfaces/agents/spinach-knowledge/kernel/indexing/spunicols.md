# kernel/indexing/spunicols.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/indexing/spunicols.m`
- Signature: `A=spunicols(A)`
- Total lines: 35

## Purpose

Sparse matrix unique-column utility. Syntax: A=spunicols(A)

## Physical / mathematical content

- Indexing utilities. These files build and transform compact index maps for basis states, matrix elements, trajectories, and tensor-product structures.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 20-21: Check consistency; implemented by `grumble(A)`.
- Lines 23-24: Return Matlab reference unique columns; implemented by `A=unique(A.','rows').'`.

### Key state/data transformations

- Lines 24: computes `A` using `A=unique(A.','rows').'`.

### Local helper functions

- Line 29: `grumble()` — `function grumble(A)`.
  - Representative operation: `if (~isnumeric(A))||(~issparse(A))||(~isreal(A))|| (~isa(A,'double'))||(~ismatrix(A))`.
  - Representative operation: `(~isa(A,'double'))||(~ismatrix(A))`.

## Parameters / inputs

- A -sparse real double matrix

## Outputs

- A -sparse real double matrix containing
- unique columns of the input matrix
- This file is a Matlab fallback for the compiled MEX function.

## Implementation structure

- Sparse matrix unique-column utility. Syntax:
- A=spunicols(A)
- A -sparse real double matrix
- A -sparse real double matrix containing
- unique columns of the input matrix
- This file is a Matlab fallback for the compiled MEX function.
- Check consistency
- Return Matlab reference unique columns
- Consistency enforcement

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `issparse()`, `ismatrix()`.
