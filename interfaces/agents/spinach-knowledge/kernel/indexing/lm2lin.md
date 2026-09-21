# kernel/indexing/lm2lin.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/indexing/lm2lin.m`
- Signature: `I=lm2lin(L,M)`
- Total lines: 60

## Purpose

Converts L,M indexing of spin states into linear indexing. In the linear indexing convention, spin states are listed in the order of increasing L rank, and, within ranks, in the order of decreasing M projection. Zero base counting is used: (L=0,M=0) -> I=0 (L=1,M=1) -> I=1 (L=1,M=0) -> I=2, et cetera...

## Physical / mathematical content

- Indexing utilities. These files build and transform compact index maps for basis states, matrix elements, trajectories, and tensor-product structures.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Syntax

```matlab
I=lm2lin(L,M)
```

## Parameters / inputs

- L -ranks of the spin states
- M -projections of the spin states

## Outputs

- I -linear indices of spin states, with
- I=0 corresponding to L=0, M=0.

## Implementation structure

- Converts L,M indexing of spin states into linear indexing. In
- the linear indexing convention, spin states are listed in the
- order of increasing L rank, and, within ranks, in the order of
- decreasing M projection. Zero base counting is used:
- (L=0,M=0) -> I=0
- (L=1,M=1) -> I=1
- (L=1,M=0) -> I=2, et cetera...
- I=lm2lin(L,M)
- L -ranks of the spin states
- M -projections of the spin states
- I -linear indices of spin states, with
- I=0 corresponding to L=0, M=0.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `any()`.
