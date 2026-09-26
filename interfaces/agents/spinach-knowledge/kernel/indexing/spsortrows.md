# kernel/indexing/spsortrows.m

- Signature: `idx=spsortrows(A)`

## Purpose

Sparse matrix row-sorting permutation utility. Syntax: idx=spsortrows(A)

## Physical / mathematical content

- Indexing utilities. These files build and transform compact index maps for basis states, matrix elements, trajectories, and tensor-product structures.

## Numerical / algorithmic content

## Parameters / inputs

- A -sparse real double matrix

## Outputs

- idx -row permutation index, matching the second
- output of Matlab's sortrows(A)
- This file is a Matlab fallback for the compiled MEX function.

## Implementation structure

- Sparse matrix row-sorting permutation utility. Syntax:
- idx=spsortrows(A)
- A -sparse real double matrix
- idx -row permutation index, matching the second
- output of Matlab's sortrows(A)
- This file is a Matlab fallback for the compiled MEX function.
- Check consistency
- Return Matlab reference permutation
- Consistency enforcement
- Nihilistic Password Security Questions
- What is the name of your least favorite child?
- In what year did you abandon your dreams?
