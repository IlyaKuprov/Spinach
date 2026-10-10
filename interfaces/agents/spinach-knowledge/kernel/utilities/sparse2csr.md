# kernel/utilities/sparse2csr.m

## Purpose

Computes a partial compressed row storage (CSR) transformation for a given MATLAB sparse matrix, adapted from code written by David Gleich. Only the index arrays are returned; the values are ignored.

## Behaviour

- Validates the input with a consistency check (`grumble`), which errors with `'A must be a sparse logical matrix.'` unless the input is simultaneously logical, a matrix, and sparse.
- Sets the problem dimensions from `size(A,1)` and `nnz(A)`.
- Obtains Cartesian indices of the nonzero entries with `find(A)`, returning row and column positions.
- Preallocates `col_idx` as a `n_nonzeros`-by-1 array and `row_ptr` as a `(matrix_dim+1)`-by-1 array.
- Counts elements per row by incrementing `row_ptr(rows(n)+1)` for each nonzero, then applies `cumsum` to obtain row offsets.
- Builds the column index array by placing each nonzero's column index at `row_ptr(rows(n))+1` and incrementing the corresponding row pointer.
- Rebuilds the row pointer array by shifting values one position (loop from `matrix_dim` down to 1 assigning `row_ptr(n+1)=row_ptr(n)`), sets `row_ptr(1)=0`, and finally adds 1 to all entries, yielding 1-based indices.

## Inputs and outputs

**Inputs**

- `A` — a MATLAB sparse matrix to be converted into the CSR format; must be a sparse logical matrix.

**Outputs**

- `row_ptr` — row pointer array of the CSR format.
- `col_idx` — column index array of the CSR format.

Syntax: `[row_ptr,col_idx]=sparse2csr(A)`.

## References

- Source: [kernel/utilities/sparse2csr.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/sparse2csr.m)
- Spin Dynamics Wiki: [sparse2csr.m](https://spindynamics.org/wiki/index.php?title=sparse2csr.m)
