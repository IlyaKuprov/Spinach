# kernel/utilities/sparse2csr.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/sparse2csr.m`
- Signature: `[row_ptr,col_idx]=sparse2csr(A)`
- Total lines: 71

## Purpose

Computes a partial compressed row storage (CSR) transformation for a given Matlab sparse matrix. Adapted from the code written by David Gleich. Only returns the index arrays and ignores the values. Syntax: [row_ptr,col_idx]=sparse2csr(A)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 26-27: Check consistency; implemented by `grumble(A)`.
- Lines 29-30: Set problem dimensions; implemented by `matrix_dim=size(A,1); n_nonzeros=nnz(A)`.
- Lines 32-33: Get Cartesian indices; implemented by `[rows,cols]=find(A)`.
- Lines 35-36: Preallocate the answer; implemented by `col_idx=zeros(n_nonzeros,1)`.
- Lines 39-40: Count row elements; implemented by `for n=1:n_nonzeros`.
- Lines 45-46: Build column index; implemented by `for n=1:n_nonzeros`.
- Lines 51-52: Build row index; implemented by `for n=matrix_dim:-1:1`.

### Control flow inferred from the code

- Line 40: `for` loop over `n=1:n_nonzeros`.
- Line 46: `for` loop over `n=1:n_nonzeros`.
- Line 52: `for` loop over `n=matrix_dim:-1:1`.

### Key state/data transformations

- Lines 30: computes `matrix_dim` using `matrix_dim=size(A,1); n_nonzeros=nnz(A)`.
- Lines 33: computes `[rows,cols]` using `[rows,cols]=find(A)`.
- Lines 36: computes `col_idx` using `col_idx=zeros(n_nonzeros,1)`.
- Lines 37: computes `row_ptr` using `row_ptr=zeros(matrix_dim+1,1)`.
- Lines 41: computes `row_ptr(rows(n)+1)` using `row_ptr(rows(n)+1)=row_ptr(rows(n)+1)+1`.
- Lines 47: computes `col_idx(row_ptr(rows(n))+1)` using `col_idx(row_ptr(rows(n))+1)=cols(n)`.
- Lines 48: computes `row_ptr(rows(n))` using `row_ptr(rows(n))=row_ptr(rows(n))+1`.
- Lines 53: computes `row_ptr(n+1)` using `row_ptr(n+1)=row_ptr(n)`.
- Lines 55: computes `row_ptr(1)` using `row_ptr(1)=0; row_ptr=row_ptr+1`.

### Local helper functions

- Line 60: `grumble()` — `function grumble(A)`. It had long since come to my attention that people of accomplishment rarely sat back and let things happen
  - Representative operation: `if (~islogical(A))||(~ismatrix(A))||(~issparse(A))`.
  - Representative operation: `error('A must be a sparse logical matrix.')`.

## Parameters / inputs

- A -a Matlab sparse matrix to be converted
- into the CSR format

## Outputs

- row_ptr -row pointer array of the CSR format
- col_idx -column index array of the CSR format

## Implementation structure

- Computes a partial compressed row storage (CSR) transformation
- for a given Matlab sparse matrix. Adapted from the code written
- by David Gleich. Only returns the index arrays and ignores the
- values. Syntax:
- [row_ptr,col_idx]=sparse2csr(A)
- A -a Matlab sparse matrix to be converted
- into the CSR format
- row_ptr -row pointer array of the CSR format
- col_idx -column index array of the CSR format
- Check consistency
- Set problem dimensions
- Get Cartesian indices

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `nnz()`, `row_ptr()`, `rows()`, `cumsum()`, `col_idx()`, `cols()`, `islogical()`, `ismatrix()`, `issparse()`.
