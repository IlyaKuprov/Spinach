# kernel/indexing/lin2lm.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/indexing/lin2lm.m`
- Signature: `[L,M]=lin2lm(I)`
- Total lines: 77

## Purpose

Converts linear indexing of spin states into L,M indexing. In the linear indexing convention, spin states are listed in the order of increasing L rank, and, within ranks, in the order of decreasing M projection. Zero base counting is used: I=0 -> (L=0,M=0), I=1 -> (L=1,M=1), I=2 -> (L=1,M=0), et cetera.

## Physical / mathematical content

- The linear index of a single-spin irreducible spherical tensor state is I=L^2+L-M, so that the rank is the integer part of the square root of I and the projection follows from the remainder; this is the state numbering used by the `sphten-liouv` basis descriptor, where 0 is the unit state of a spin.

## Numerical / algorithmic content

- The rank is the integer part of `sqrt(I)`, the one floating-point step, taken in double precision and returned in the class of the input, then stepped back by one wherever the root of an index just below a large perfect square rounded up; the projection is evaluated as `L^2-I+L` in that class, an order in which every intermediate lies between -2L and L, so no signed integer class that holds I can overflow. A final check that every projection lies within its rank guards the regime above `flintmax`, where `double(I)` is no longer exact.
- Integer inputs are accepted only in signed classes: projections are negative for half of the states of every rank, and an unsigned class would saturate them to zero, so the grumbler refuses `uint8`, `uint16`, `uint32`, and `uint64` inputs.
- Integer classes are admitted up to the last rank they hold completely (index 120 for `int8`, 32760 for `int16`), the rule by which `basis.m` picks the descriptor class, so every accepted input round-trips through `lm2lin` in its class.
- The `sphten-liouv` descriptor `spin_system.bas.basis` is sparse single, so its ranks and projections come back sparse single; the per-subgraph descriptor blocks inside `basis.m` are `int8` or `int16` and come back in that class.

## Syntax

```matlab
[L,M]=lin2lm(I)
```

## Parameters / inputs

- I - linear indices of spin states, with I=0 corresponding to L=0, M=0; double, single, or a signed integer class

## Outputs

- L - ranks of the spin states, same class and sparsity as the input
- M - projections of the spin states, same class and sparsity as the input
