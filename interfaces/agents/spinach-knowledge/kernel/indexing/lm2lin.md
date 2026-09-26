# kernel/indexing/lm2lin.m

- Signature: `I=lm2lin(L,M)`

## Purpose

Converts L,M indexing of spin states into linear indexing. In the linear indexing convention, spin states are listed in the order of increasing L rank, and, within ranks, in the order of decreasing M projection. Zero base counting is used: (L=0,M=0) -> I=0 (L=1,M=1) -> I=1 (L=1,M=0) -> I=2, et cetera...

## Physical / mathematical content

- Indexing utilities. These files build and transform compact index maps for basis states, matrix elements, trajectories, and tensor-product structures.

## Numerical / algorithmic content

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
