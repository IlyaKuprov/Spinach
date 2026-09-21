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
