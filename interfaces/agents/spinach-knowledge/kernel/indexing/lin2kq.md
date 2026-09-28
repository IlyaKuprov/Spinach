# kernel/indexing/lin2kq.m

- Signature: `[K,Q]=lin2kq(N,I,idx_base)`

## Purpose

Converts linear serpentine indices of an N-by-N matrix into row and column indices. It supports base-1 indices (1 through N) and base-0 indices (0 through N-1).

## Physical / mathematical content

This is the inverse indexing conversion to `kq2lin`; it does not alter matrix values or represent a physical operation.

## Numerical / algorithmic content

The function constructs `serpentine(N,idx_base)` and finds the matrix location whose value equals each element of `I`. For base 0, it subtracts 1 from the MATLAB row and column locations before returning them. The implementation validates the input types, integer values, base, and index range.

## Syntax

```matlab
[K,Q]=lin2kq(N,I,idx_base)
```

## Parameters / inputs

- `N` - positive real integer matrix dimension.
- `I` - real integer linear indices, in an array of any size, between `idx_base` and `N^2-1+idx_base`.
- `idx_base` - indexing base, either 0 or 1.

## Outputs

- `K` - row indices, with the same size as `I`.
- `Q` - column indices, with the same size as `I`.
