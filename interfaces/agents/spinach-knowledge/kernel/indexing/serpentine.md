# kernel/indexing/serpentine.m

- Signature: `S=serpentine(nlevels,idx_base)`

## Purpose

Builds the square serpentine index matrix used to number matrix elements with a single index. Entries increase by diagonals of constant row-plus-column, with ties ordered from the larger row index to the smaller one. For `nlevels=4`, the base-1 matrix is:

```text
1  3  6 10
2  5  9 13
4  8 12 15
7 11 14 16
```

The base-0 matrix has 1 subtracted from every entry.

## Physical / mathematical content

This matrix is an indexing map, not a physical matrix operation.

## Numerical / algorithmic content

The implementation creates all row and column coordinates, sorts them by increasing `row+column` and then decreasing row, and assigns the sequence `1:nlevels^2` in that order. For base 0 it subtracts 1 from the completed matrix.

## Parameters / inputs

- `nlevels` - positive real integer dimension of the square matrix.
- `idx_base` - indexing base, either 0 or 1.

## Outputs

- `S` - `nlevels`-by-`nlevels` serpentine index matrix, expressed in the selected base.
