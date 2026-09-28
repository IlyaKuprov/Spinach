# kernel/indexing/lin2lmn.m

- Signature: `[L,M,N]=lin2lmn(I)`

## Purpose

Converts one-based linear indices of Wigner D functions to rank `L` and indices `M,N`. Indices are ordered by increasing `L`; within each rank, `M` decreases from `L` to `-L`, and for each `M`, `N` decreases from `L` to `-L`. Thus `I=1` maps to (0,0,0), `I=2` to (1,1,1), and `I=3` to (1,1,0).

## Physical / mathematical content

A rank-`L` block contains `(2*L+1)^2` Wigner D indices. The first linear index in that block is `(4*L^3-L)/3+1`.

## Numerical / algorithmic content

The implementation inverts the cumulative rank-block count to obtain `L`, computes the position within that block to obtain `M,N`, then checks the result by converting back with `lmn2lin`.

## Syntax

```matlab
[L,M,N]=lin2lmn(I)
```

## Parameters / inputs

- `I` - array of positive real integers; `I=1` corresponds to `L=0, M=0, N=0`.

## Outputs

- `L` - ranks of Wigner D functions.
- `M` - row indices of Wigner D functions.
- `N` - column indices of Wigner D functions.
