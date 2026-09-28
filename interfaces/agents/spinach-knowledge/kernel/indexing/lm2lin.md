# kernel/indexing/lm2lin.m

- Signature: `I=lm2lin(L,M)`

## Purpose

Converts total angular momentum rank `L` and projection `M` to the zero-based linear index of a spin state. Within each rank, projections are ordered from `L` down to `-L`; for example, (0,0) maps to 0, (1,1) to 1, and (1,0) to 2.

## Physical / mathematical content

The zero-based index is `I=L^2+L-M`. The allowed projection condition is `abs(M)<=L`, with `L>=0`.

## Numerical / algorithmic content

The implementation applies that element-wise formula after checking that `L` and `M` are real integer arrays of equal size, and that the rank and projection bounds are satisfied.

## Syntax

```matlab
I=lm2lin(L,M)
```

## Parameters / inputs

- `L` - non-negative integer ranks of the spin states.
- `M` - integer projections with `abs(M)<=L`; must have the same size as `L`.

## Outputs

- `I` - zero-based linear indices of spin states; `I=0` corresponds to `L=0, M=0`.
