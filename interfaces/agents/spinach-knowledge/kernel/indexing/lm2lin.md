# kernel/indexing/lm2lin.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/indexing/lm2lin.m`
- Signature: `I=lm2lin(L,M)`
- Total lines: 74

## Purpose

Converts L,M indexing of spin states into linear indexing. In the linear indexing convention, spin states are listed in the order of increasing L rank, and, within ranks, in the order of decreasing M projection. Zero base counting is used: (L=0,M=0) -> I=0, (L=1,M=1) -> I=1, (L=1,M=0) -> I=2, et cetera.

## Physical / mathematical content

- The linear index of a single-spin irreducible spherical tensor state is I=L^2+L-M; this is the inverse of `lin2lm` and the state numbering of the `sphten-liouv` basis descriptor.

## Numerical / algorithmic content

- The index is evaluated as `L^2-M+L` in the class of the inputs, an order in which the intermediate never exceeds the result, so no integer class that holds the result can overflow; double, single, and sparse inputs keep their class and sparsity, and every pair returned by `lin2lm` round-trips exactly.
- The grumbler requires real integer-valued arrays of the same size and class, non-negative ranks, projections within plus or minus the rank, and, for integer classes, ranks no higher than the last rank whose complete index range fits the class (10 for `int8`, 180 for `int16`, 46339 for `int32`); this is the bound by which `basis.m` chooses the descriptor class, and `lin2lm` admits the same domain. The bound is enforced before the arithmetic because a check after it cannot be trusted in saturating integer arithmetic: `int8` rank 12 with projection 12 evaluates to 127, and 127-127+12 is again 12.

## Syntax

```matlab
I=lm2lin(L,M)
```

## Parameters / inputs

- L - ranks of the spin states; double, single, or a signed integer class
- M - projections of the spin states, same class as L

## Outputs

- I - linear indices of spin states, with I=0 corresponding to L=0, M=0; same class and sparsity as L
