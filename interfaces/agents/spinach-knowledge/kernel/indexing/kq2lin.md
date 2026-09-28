# kernel/indexing/kq2lin.m

- Signature: `I=kq2lin(N,K,Q,idx_base)`

## Purpose

Maps paired row and column indices to the single linear index assigned by Spinach's serpentine ordering of an N-by-N matrix. The convention is selectable: base 1 indices range from 1 to N; base 0 indices range from 0 to N-1. For a 3-by-3 matrix, the base-1 map is `[1 3 6; 2 5 8; 4 7 9]`, and the base-0 map is that map minus 1.

## Physical / mathematical content

This is an indexing conversion; it does not alter matrix values or represent a physical operation.

## Numerical / algorithmic content

The function obtains the map from `serpentine(N,idx_base)` and looks up each (K,Q) pair. For base 0, it adds 1 to K and Q for MATLAB array indexing; the returned map values remain zero-based.

## Syntax

```matlab
I=kq2lin(N,K,Q,idx_base)
```

## Parameters / inputs

- `N` - scalar real integer matrix dimension; it must be at least `idx_base`.
- `K` - real integer row indices, in an array of any size.
- `Q` - real integer column indices, with the same size as `K`.
- `idx_base` - indexing base, either 0 or 1. Every index must lie between `idx_base` and `N-1+idx_base`.

## Outputs

- `I` - linear serpentine indices, with the same size as the input index arrays.
