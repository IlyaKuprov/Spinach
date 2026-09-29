# kernel/indexing/serpentine.m

- Signature: `S=serpentine(nlevels,idx_base)`
- Direct MATLAB source: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/indexing/serpentine.m>
- Spin Dynamics Wiki: <https://spindynamics.org/wiki/index.php?title=serpentine.m>

## Purpose

Builds the square index matrix used by Spinach to number matrix elements with a single index. Values are assigned in increasing order along diagonals of constant row-plus-column; ties are ordered from larger row number to smaller row number. This is a data-layout helper, not a physical interaction or dynamics operation: it does not change matrix values, states, operators, or Hamiltonians.

## Exact construction

The code creates 1-based row and column coordinate grids with `ndgrid(1:nlevels)`, sorts the coordinate pairs by `[row+column, -row]`, and assigns consecutive values `1:nlevels^2` to the sorted positions. When `idx_base=0`, it subtracts one from every completed entry; for `idx_base=1`, the assigned values are retained.

For `nlevels=4`, the base-1 matrix is:

```text
1  3  6 10
2  5  9 13
4  8 12 15
7 11 14 16
```

The corresponding base-0 matrix is:

```text
0  2  5  9
1  4  8 12
3  7 11 14
6 10 13 15
```

## Syntax and guards

`S=serpentine(nlevels,idx_base)`

- `nlevels` - positive real numeric integer scalar; determines the number of rows and columns.
- `idx_base` - real numeric scalar equal to 0 or 1.
- `S` - `nlevels`-by-`nlevels` index matrix in the requested base.
