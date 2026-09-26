# kernel/indexing/lin2lm.m

- Signature: `[L,M]=lin2lm(I)`

## Purpose

Converts linear indexing of spin states into L,M indexing. In the linear indexing convention, spin states are listed in the order of increasing L rank, and, within ranks, in the order of decreasing M projection. Zero base counting is used: I=0 -> (L=0,M=0) I=1 -> (L=1,M=1) I=2 -> (L=1,M=0), et cetera...

## Physical / mathematical content

- Indexing utilities. These files build and transform compact index maps for basis states, matrix elements, trajectories, and tensor-product structures.

## Numerical / algorithmic content

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
