# kernel/indexing/lmn2lin.m

- Signature: `I=lmn2lin(L,M,N)`

## Purpose

Converts rank `L` and indices `M,N` of Wigner D functions to one-based linear indices. The order is increasing in `L`; within each rank, `M` decreases from `L` to `-L`, and within each `M`, `N` decreases from `L` to `-L`. For example, (0,0,0) maps to 1, (1,1,1) to 2, and (1,1,0) to 3.

## Physical / mathematical content

Each rank `L` contributes `(2*L+1)^2` indices. The within-rank ordering follows the listed descending `M` and `N` values.

## Numerical / algorithmic content

The element-wise formula used is `I=L*(4*L^2+6*(L-M)+5)/3-M-N+1`. The inputs must be real integers, with `L>=0`, `abs(M)<=L`, and `abs(N)<=L`; all three arrays must have the same size.

## Syntax

```matlab
I=lmn2lin(L,M,N)
```

## Parameters / inputs

- `L` - non-negative integer ranks of Wigner D functions.
- `M` - integer row indices satisfying `abs(M)<=L`.
- `N` - integer column indices satisfying `abs(N)<=L`.

## Outputs

- `I` - one-based linear indices of Wigner D functions; `I=1` corresponds to `L=0, M=0, N=0`.
