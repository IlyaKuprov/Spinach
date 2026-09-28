# kernel/indexing/lin2lm.m

- Signature: `[L,M]=lin2lm(I)`

## Purpose

Converts zero-based linear indices of spin states to total angular momentum rank `L` and projection `M`. States are ordered by increasing `L`, and within each rank by decreasing `M`: `I=0` maps to (0,0), `I=1` to (1,1), and `I=2` to (1,0).

## Physical / mathematical content

For rank `L`, the allowed projections run from `L` down to `-L`; each rank therefore contributes `2*L+1` consecutive indices.

## Numerical / algorithmic content

For each input, the implementation sets `L=floor(sqrt(I))` and `M=L^2+L-I`. It then checks the conversion by applying `lm2lin` and comparing the result with `I`.

## Syntax

```matlab
[L,M]=lin2lm(I)
```

## Parameters / inputs

- `I` - array of non-negative real integers; `I=0` corresponds to `L=0, M=0`.

## Outputs

- `L` - ranks of the spin states.
- `M` - projections of the spin states.
