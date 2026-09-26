# kernel/indexing/serpentine.m

- Signature: `S=serpentine(nlevels,idx_base)`

## Purpose

Serpentine index matrix used in Spinach for single-index numbering of matrix elements. Syntax: S=serpentine(nlevels,idx_base)

## Physical / mathematical content

- Indexing utilities. These files build and transform compact index maps for basis states, matrix elements, trajectories, and tensor-product structures.

## Numerical / algorithmic content

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
