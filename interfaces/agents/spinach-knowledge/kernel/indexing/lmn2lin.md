# kernel/indexing/lmn2lin.m

- Signature: `I=lmn2lin(L,M,N)`

## Purpose

Converts L,M,N indices of Wigner D functions into linear indices. In the linear indexing convention, Wigner D functions are listed in the order of increasing L rank. Within each L, the functions are listed in the order of decreasing left index M, and, for each M, in the or- der of decreasing N index. One base counting is used: (L=0,M=0,N=0) -> I=1 (L=1,M=1,N=1) -> I=2 (L=1,M=1,N=0) -> I=3, et cetera...

## Physical / mathematical content

- Indexing utilities. These files build and transform compact index maps for basis states, matrix elements, trajectories, and tensor-product structures.

## Numerical / algorithmic content

## Syntax

```matlab
I=lmn2lin(L,M,N)
```

## Parameters / inputs

- L -ranks of Wigner D functions
- M -row indices of Wigner D functions
- N -column indices of Wigner D functions

## Outputs

- I -linear indices of Wigner D functions, with
- I=1 corresponding to L=0, M=0, N=0.

## Implementation structure

- Converts L,M,N indices of Wigner D functions into linear indices. In
- the linear indexing convention, Wigner D functions are listed in the
- order of increasing L rank. Within each L, the functions are listed
- in the order of decreasing left index M, and, for each M, in the or-
- der of decreasing N index. One base counting is used:
- (L=0,M=0,N=0) -> I=1
- (L=1,M=1,N=1) -> I=2
- (L=1,M=1,N=0) -> I=3, et cetera...
- I=lmn2lin(L,M,N)
- L -ranks of Wigner D functions
- M -row indices of Wigner D functions
- N -column indices of Wigner D functions
