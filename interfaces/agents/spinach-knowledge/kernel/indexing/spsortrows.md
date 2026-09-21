# kernel/indexing/spsortrows.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/indexing/spsortrows.m`
- Signature: `idx=spsortrows(A)`
- Total lines: 63

## Purpose

Sparse matrix row-sorting permutation utility. Syntax: idx=spsortrows(A)

## Physical / mathematical content

- Indexing utilities. These files build and transform compact index maps for basis states, matrix elements, trajectories, and tensor-product structures.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `sortrows()`, `issparse()`, `ismatrix()`.
