# kernel/indexing/spsortrows.m

- Signature: `idx=spsortrows(A)`

## Purpose

Returns the row permutation that sorts a sparse matrix lexicographically by rows, matching the second output of MATLAB's `sortrows(A)`.

## Physical / mathematical content

This is a sparse-matrix ordering utility; it does not modify `A` or return the sorted matrix.

## Numerical / algorithmic content

This MATLAB fallback validates that `A` is a sparse, real, double matrix, then returns the permutation from `[~,idx]=sortrows(A)`. Spinach also provides a compiled MEX implementation.

## Syntax

```matlab
idx=spsortrows(A)
```

## Parameters / inputs

- `A` - sparse real double matrix.

## Outputs

- `idx` - row permutation index, matching the second output of MATLAB's `sortrows(A)`.
