# kernel/indexing/serpentine.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/indexing/serpentine.m`
- Signature: `S=serpentine(nlevels,idx_base)`
- Total lines: 65

## Purpose

Serpentine index matrix used in Spinach for single-index numbering of matrix elements. Syntax: S=serpentine(nlevels,idx_base)

## Physical / mathematical content

- Indexing utilities. These files build and transform compact index maps for basis states, matrix elements, trajectories, and tensor-product structures.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 35-36: Check consistency; implemented by `grumble(nlevels,idx_base)`.
- Lines 38-39: Build the serpentine matrix; implemented by `[rows,cols]=ndgrid(1:nlevels)`.
- Lines 43-44: Adjust the indexing base; implemented by `if idx_base==0, S=S-1; end`.

### Control flow inferred from the code

- Line 44: conditional branch on `idx_base==0, S=S-1; end`.

### Key state/data transformations

- Lines 39: computes `[rows,cols]` using `[rows,cols]=ndgrid(1:nlevels)`.
- Lines 40: computes `[~,idx]` using `[~,idx]=sortrows([rows(:)+cols(:), -rows(:)])`.
- Lines 41: computes `S` using `S=zeros(nlevels,nlevels); S(idx)=1:nlevels^2`.

### Local helper functions

- Line 49: `grumble()` — `function grumble(nlevels,idx_base)`. There are only two kinds of languages: the ones people
  - Representative operation: `if (~isnumeric(nlevels))||(~isreal(nlevels))|| (~isscalar(nlevels))||(nlevels<1)||(mod(nlevels,1)~=0)`.
  - Representative operation: `(~isscalar(nlevels))||(nlevels<1)||(mod(nlevels,1)~=0)`.

## Parameters / inputs

- nlevels -dimension of the matrix, a
- positive real integer
- idx_base -indexing base, 0 or 1

## Outputs

- S -serpentine matrix, for example (base 1):
- (1 )(3 )(6 )(10)
- (2 )(5 )(9 )(13)
- (4 )(8 )(12)(15)
- (7 )(11)(14)(16)
- or (with indexing set to base 0):
- (0 )(2 )(5 )(9 )
- (1 )(4 )(8 )(12)
- (3 )(7 )(11)(14)
- (6 )(10)(13)(15)

## Implementation structure

- Serpentine index matrix used in Spinach for single-index
- numbering of matrix elements. Syntax:
- S=serpentine(nlevels,idx_base)
- nlevels -dimension of the matrix, a
- positive real integer
- idx_base -indexing base, 0 or 1
- S -serpentine matrix, for example (base 1):
- (1 )(3 )(6 )(10)
- (2 )(5 )(9 )(13)
- (4 )(8 )(12)(15)
- (7 )(11)(14)(16)
- or (with indexing set to base 0):

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `sortrows()`, `rows()`, `cols()`, `isscalar()`, `ismember()`.
