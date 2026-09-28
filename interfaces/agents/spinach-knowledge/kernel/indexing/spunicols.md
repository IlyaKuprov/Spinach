# kernel/indexing/spunicols.m

- Signature: `A=spunicols(A)`

## Purpose

Returns a sparse matrix containing one copy of each distinct column of the input matrix.

## Physical / mathematical content

This utility operates on matrix columns; it does not model a physical system.

## Numerical / algorithmic content

The function transposes `A`, applies Matlab's `unique(...,'rows')` to the rows, then transposes the result back. Duplicate columns are therefore retained only once.

## Parameters / inputs

- `A` — a sparse, real, double matrix.

## Outputs

- `A` — a sparse real double matrix containing the unique columns of the input.

## Implementation structure

The input is checked for consistency, then the Matlab fallback computes `unique(A.','rows').'`. The file serves as the reference implementation for the compiled MEX function.
