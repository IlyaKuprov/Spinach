# kernel/indexing/lin2lmn.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/indexing/lin2lmn.m`
- Signature: `[L,M,N]=lin2lmn(I)`
- Total lines: 66

## Purpose

Converts linear indices of Wigner D functions into L,M,N indices. In the linear indexing convention, Wigner D functions are listed in the order of increasing L rank. Within each L, the functions are listed in the order of decreasing left index M, and, for each M, in the or- der of decreasing N index. One base counting is used: I=1 -> (L=0,M=0,N=0) I=2 -> (L=1,M=1,N=1) I=3 -> (L=1,M=1,N=0), et cetera...

## Physical / mathematical content

- Indexing utilities. These files build and transform compact index maps for basis states, matrix elements, trajectories, and tensor-product structures.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Syntax

```matlab
[L,M,N]=lin2lmn(I)
```

## Parameters / inputs

- I -linear indices of Wigner D functions, with
- I=1 corresponding to L=0, M=0, N=0.

## Outputs

- L -ranks of Wigner D functions
- M -row indices of Wigner D functions
- N -column indices of Wigner D functions

## Implementation structure

- Converts linear indices of Wigner D functions into L,M,N indices. In
- the linear indexing convention, Wigner D functions are listed in the
- order of increasing L rank. Within each L, the functions are listed
- in the order of decreasing left index M, and, for each M, in the or-
- der of decreasing N index. One base counting is used:
- I=1 -> (L=0,M=0,N=0)
- I=2 -> (L=1,M=1,N=1)
- I=3 -> (L=1,M=1,N=0), et cetera...
- [L,M,N]=lin2lmn(I)
- I -linear indices of Wigner D functions, with
- I=1 corresponding to L=0, M=0, N=0.
- L -ranks of Wigner D functions

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `fix()`, `nnz()`, `lmn2lin()`, `any()`.
