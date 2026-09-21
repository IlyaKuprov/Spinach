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
