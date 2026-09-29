# kernel/indexing/lin2lmn.m

- Signature: `[L,M,N]=lin2lmn(I)`

## Purpose

Converts one-based linear indices of Wigner D functions into rank `L`, left index `M`, and right index `N`. Ranks increase first. Within each rank, `M` decreases from `L` to `-L`; for each `M`, `N` decreases from `L` to `-L`. The original examples are `I=1 -> (L=0,M=0,N=0)`, `I=2 -> (L=1,M=1,N=1)`, and `I=3 -> (L=1,M=1,N=0)`.

## Physical / mathematical content

This routine only decodes an array index into the three coordinates of a Wigner D function. It does not evaluate a Wigner D function, transform a spin state or operator, or specify a physical interaction. The coordinates are dimensionless, with no sign or unit convention introduced here.

## Numerical / algorithmic content

The code identifies the rank block by solving for `L` from the cumulative count `(4*L^3-L)/3`. It defines the zero-based position within that block as `p=I-(4*L^3-L)/3-1`, then computes `M=L-floor(p/(2*L+1))` and `N=L+(2*L+1)*(L-M)-p`. A rank-L block has `(2*L+1)^2` entries. The input must be positive real integers and may have any array shape; all three outputs preserve that shape. A consistency check calls `lmn2lin(L,M,N)` and raises an error if the round trip does not reproduce `I`.

## Syntax

```matlab
[L,M,N]=lin2lmn(I)
```

## Parameters / inputs

- `I` - array of positive real integers, using one-based indexing; `I=1` corresponds to `L=0, M=0, N=0`.

## Outputs

- `L` - Wigner D ranks.
- `M` - left (row) indices.
- `N` - right (column) indices.

## Sources

- [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/indexing/lin2lmn.m) (local path: `kernel/indexing/lin2lmn.m`).
- [Existing Wiki page](https://spindynamics.org/wiki/index.php?title=lin2lmn.m).
