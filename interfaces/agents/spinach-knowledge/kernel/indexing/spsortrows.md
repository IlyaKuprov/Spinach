# kernel/indexing/spsortrows.m

- Signature: `idx=spsortrows(A)`
- Direct MATLAB source: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/indexing/spsortrows.m>

## Purpose

Returns the row permutation for sorting a sparse matrix by rows, matching the second output of MATLAB's `sortrows(A)`. This is an indexing utility: it returns neither a sorted matrix nor a modified input, and has no physical interaction, state/operator effect, Hamiltonian, or unit convention.

## Execution and guards

This file is the MATLAB fallback for the compiled MEX function. It first requires `A` to be numeric, sparse, real, double precision, and two-dimensional (the source uses `ismatrix`). It then calls `[~,idx]=sortrows(A)`; sorting and tie behaviour therefore follow MATLAB's `sortrows` implementation. No additional row-count, column-count, or application-specific guard is present in this fallback.

## Syntax and arguments

`idx=spsortrows(A)`

- `A` - sparse real double matrix.
- `idx` - row permutation index, corresponding to the second output of `sortrows(A)`.
