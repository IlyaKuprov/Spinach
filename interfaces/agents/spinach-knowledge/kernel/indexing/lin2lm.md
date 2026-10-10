# kernel/indexing/lin2lm.m

- Signature: `[L,M]=lin2lm(I)`

## Purpose

Converts zero-based linear indices of Spinach spin states into rank `L` and projection `M`. Indices are grouped by increasing rank; within rank `L`, the order is `M=L,L-1,...,-L`. The preserved examples are `I=0 -> (L=0,M=0)`, `I=1 -> (L=1,M=1)`, and `I=2 -> (L=1,M=0)`.

## Physical / mathematical content

The outputs label a spin-state coordinate by angular-momentum rank and projection; this routine only converts its index and does not construct or modify the spin state, a Hamiltonian, or an operator. There is no interaction/formalism selection, physical sign convention, or unit conversion in this indexing routine; the indices are dimensionless.

## Numerical / algorithmic content

Elementwise, the code computes `L=floor(sqrt(I))` and `M=L^2+L-I`. Therefore the zero-based block for rank `L` runs from `I=L^2` through `I=(L+1)^2-1`, mapping to projections from `+L` down to `-L`. The input may be an array of any shape, and both outputs retain that shape. A consistency check calls `lm2lin(L,M)` and raises an error if the round trip does not reproduce `I`.

## Syntax

```matlab
[L,M]=lin2lm(I)
```

## Parameters / inputs

- `I` - array of non-negative real integers, in zero-based indexing; `I=0` corresponds to `L=0, M=0`.

## Outputs

- `L` - spin-state ranks.
- `M` - spin-state projections.

## Sources

- [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/indexing/lin2lm.m) (local path: `kernel/indexing/lin2lm.m`).
- [Existing Wiki page](https://spindynamics.org/wiki/index.php?title=lin2lm.m).
