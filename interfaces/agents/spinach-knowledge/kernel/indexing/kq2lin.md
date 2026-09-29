# kernel/indexing/kq2lin.m

- Signature: `I=kq2lin(N,K,Q,idx_base)`

## Purpose

Converts paired matrix row and column indices into Spinach's linear serpentine index. For an N-by-N matrix, coordinates are ordered by increasing row-plus-column sum; within one diagonal, rows are ordered from larger to smaller. The helper `serpentine` builds that index map, and this function looks up each (K,Q) pair in it.

The original 3-by-3 examples are:

- Base 1: `[1 3 6; 2 5 8; 4 7 9]`.
- Base 0: `[0 2 5; 1 4 7; 3 6 8]`.

Thus base 1 maps (1,2) to 3, while base 0 maps (0,1) to 2.

## Physical / mathematical content

This is a matrix-index conversion, not a Hamiltonian, interaction, or state/operator transformation. It changes neither matrix elements nor their physical interpretation. The indices are dimensionless; no sign or unit convention applies.

## Numerical / algorithmic content

The function constructs `S=serpentine(N,idx_base)`, then performs a direct lookup at `S(K,Q)` for base 1 or `S(K+1,Q+1)` for base 0. Each output is the index assigned to that matrix coordinate, and `I` has the same array shape as `K` and `Q`.

## Syntax

```matlab
I=kq2lin(N,K,Q,idx_base)
```

## Parameters / inputs

- `N` - positive integer matrix dimension, supplied as a scalar.
- `K` - numeric real integer row-index array.
- `Q` - numeric real integer column-index array, with the same size as `K`.
- `idx_base` - scalar numeric real indexing base, either 0 or 1. Each coordinate must be in the inclusive range `idx_base:N-1+idx_base`.

## Outputs

- `I` - linear serpentine indices, with the same size as the input index arrays. The range is 0 through N^2-1 for base 0, and 1 through N^2 for base 1.

## Sources

- [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/indexing/kq2lin.m) (local path: `kernel/indexing/kq2lin.m`).
- [Existing Wiki page](https://spindynamics.org/wiki/index.php?title=kq2lin.m).
