# kernel/indexing/lin2kq.m

- Signature: `[K,Q]=lin2kq(N,I,idx_base)`

## Purpose

Converts each linear serpentine matrix index back into its row and column coordinates. For an N-by-N matrix, the serpentine ordering runs over increasing row-plus-column diagonals, with row decreasing within each diagonal. The original 3-by-3 index maps are `[1 3 6; 2 5 8; 4 7 9]` in base 1 and `[0 2 5; 1 4 7; 3 6 8]` in base 0. For example, base 1 index 5 maps to (2,2); base 0 index 4 maps to (1,1).

## Physical / mathematical content

This is an indexing conversion only. It does not change matrix values, represent an interaction, or act on a physical state or operator. The inputs and outputs are dimensionless indices, so no sign or unit convention applies.

## Numerical / algorithmic content

The function builds `S=serpentine(N,idx_base)`. For every element of `I`, it finds the row and column of the matching entry in `S`. With base 0, it subtracts one from the row and column returned by MATLAB's one-based array lookup. Both outputs preserve the size of `I`.

## Syntax

```matlab
[K,Q]=lin2kq(N,I,idx_base)
```

## Parameters / inputs

- `N` - positive integer matrix dimension, supplied as a scalar.
- `I` - numeric real integer array of linear indices. The inclusive range is `idx_base:N^2-1+idx_base`.
- `idx_base` - scalar numeric real indexing base, either 0 or 1.

## Outputs

- `K` - row indices, with the same size as `I`.
- `Q` - column indices, with the same size as `I`.

## Sources

- [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/indexing/lin2kq.m) (local path: `kernel/indexing/lin2kq.m`).
- [Existing Wiki page](https://spindynamics.org/wiki/index.php?title=lin2kq.m).
