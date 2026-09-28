# kernel/utilities/sparse2csr.m

- Signature: `[row_ptr,col_idx]=sparse2csr(A)`

## Purpose

Computes a partial compressed row storage (CSR) transformation for a given MATLAB sparse logical matrix. Adapted from code written by David Gleich. Returns the row-pointer and column-index arrays; it does not return the matrix values.

## Physical / mathematical content

- General mathematical and infrastructure utility for sparse-matrix storage.

## Numerical / algorithmic content

- Obtains the row and column coordinates of nonzero entries with `find(A)`.
- Counts entries per row, takes cumulative counts, and applies a final shift to produce a one-based CSR pointer array, `row_ptr`.
- Places column indices in `col_idx` according to the row pointers.

## Parameters / inputs

- `A` — a sparse logical matrix to be converted partially to CSR format.

## Outputs

- `row_ptr` — one-based CSR row-pointer array of length `size(A,1)+1`.
- `col_idx` — CSR column-index array of length `nnz(A)`.

## Implementation structure

- Checks that `A` is logical, two-dimensional, and sparse.
- Sets the row count with `size(A,1)` and the nonzero count with `nnz(A)`.
- Preallocates `row_ptr` and `col_idx` using those counts, then builds the row pointers and column indices.

Source: <https://spindynamics.org/wiki/index.php?title=sparse2csr.m>

dgleich@purdue.edu

ilya.kuprov@weizmann.ac.il