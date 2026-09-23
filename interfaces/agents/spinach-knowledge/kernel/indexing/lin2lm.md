# kernel/indexing/lin2lm.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/indexing/lin2lm.m`
- Signature: `[L,M]=lin2lm(I)`
- Total lines: 67

## Purpose

Converts linear indexing of spin states into L,M indexing. In the linear indexing convention, spin states are listed in the order of increasing L rank, and, within ranks, in the order of decreasing M projection. Zero base counting is used: I=0 -> (L=0,M=0), I=1 -> (L=1,M=1), I=2 -> (L=1,M=0), et cetera.

## Physical / mathematical content

- The linear index of a single-spin irreducible spherical tensor state is I=L^2+L-M, so that the rank is the integer part of the square root of I and the projection follows from the remainder; this is the state numbering used by the `sphten-liouv` basis descriptor, where 0 is the unit state of a spin.

## Numerical / algorithmic content

- The rank is `fix(sqrt(I))` and the projection is `L^2+L-I`, both evaluated in double precision whatever the class of the input, and the pair is checked against `lm2lin` before it is returned; the outputs are cast back to the class of the input, so a sparse double descriptor gives sparse double ranks and projections, and an `int8` or `int16` descriptor block gives ranks and projections of the same class.
- Integer inputs are accepted only in signed classes: projections are negative for half of the states of every rank, and an unsigned class would saturate them to zero, so the grumbler refuses `uint8`, `uint16`, `uint32`, and `uint64` inputs.
- The double-precision evaluation is what makes integer inputs safe: the intermediate L^2+L exceeds the largest representable value of a narrow integer class long before the indices themselves do.

## Syntax

```matlab
[L,M]=lin2lm(I)
```

## Parameters / inputs

- I - linear indices of spin states, with I=0 corresponding to L=0, M=0; double, single, or a signed integer class

## Outputs

- L - ranks of the spin states, same class and sparsity as the input
- M - projections of the spin states, same class and sparsity as the input
