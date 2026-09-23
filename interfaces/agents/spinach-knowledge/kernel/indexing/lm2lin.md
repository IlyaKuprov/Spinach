# kernel/indexing/lm2lin.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/indexing/lm2lin.m`
- Signature: `I=lm2lin(L,M)`
- Total lines: 68

## Purpose

Converts L,M indexing of spin states into linear indexing. In the linear indexing convention, spin states are listed in the order of increasing L rank, and, within ranks, in the order of decreasing M projection. Zero base counting is used: (L=0,M=0) -> I=0, (L=1,M=1) -> I=1, (L=1,M=0) -> I=2, et cetera.

## Physical / mathematical content

- The linear index of a single-spin irreducible spherical tensor state is I=L^2+L-M; this is the inverse of `lin2lm` and the state numbering of the `sphten-liouv` basis descriptor.

## Numerical / algorithmic content

- The index is evaluated in double precision whatever the class of the inputs and cast back to the class of L, so integer rank and projection arrays return an integer index without the intermediate L^2+L saturating a narrow integer class.
- The grumbler requires real integer-valued arrays of the same size with non-negative ranks and projections within plus or minus the rank.

## Syntax

```matlab
I=lm2lin(L,M)
```

## Parameters / inputs

- L - ranks of the spin states; double, single, or a signed integer class
- M - projections of the spin states, same class as L

## Outputs

- I - linear indices of spin states, with I=0 corresponding to L=0, M=0; same class and sparsity as L
