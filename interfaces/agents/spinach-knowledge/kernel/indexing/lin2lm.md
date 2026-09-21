# kernel/indexing/lin2lm.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/indexing/lin2lm.m`
- Signature: `[L,M]=lin2lm(I)`
- Total lines: 56

## Purpose

Converts linear indexing of spin states into L,M indexing. In the linear indexing convention, spin states are listed in the order of increasing L rank, and, within ranks, in the order of decreasing M projection. Zero base counting is used: I=0 -> (L=0,M=0) I=1 -> (L=1,M=1) I=2 -> (L=1,M=0), et cetera...

## Physical / mathematical content

- Indexing utilities. These files build and transform compact index maps for basis states, matrix elements, trajectories, and tensor-product structures.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Syntax

```matlab
[L,M]=lin2lm(I)
```

## Parameters / inputs

- I -linear indices of spin states, with
- I=0 corresponding to L=0, M=0.

## Outputs

- L -ranks of the spin states
- M -projections of the spin states

## Implementation structure

- Converts linear indexing of spin states into L,M indexing. In
- the linear indexing convention, spin states are listed in the
- order of increasing L rank, and, within ranks, in the order of
- decreasing M projection. Zero base counting is used:
- I=0 -> (L=0,M=0)
- I=1 -> (L=1,M=1)
- I=2 -> (L=1,M=0), et cetera...
- [L,M]=lin2lm(I)
- I -linear indices of spin states, with
- I=0 corresponding to L=0, M=0.
- L -ranks of the spin states
- M -projections of the spin states

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `fix()`, `nnz()`, `lm2lin()`, `any()`.
