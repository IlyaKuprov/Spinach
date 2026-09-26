# kernel/indexing/spunicols.m

- Signature: `A=spunicols(A)`

## Purpose

Sparse matrix unique-column utility. Syntax: A=spunicols(A)

## Physical / mathematical content

- Indexing utilities. These files build and transform compact index maps for basis states, matrix elements, trajectories, and tensor-product structures.

## Numerical / algorithmic content

## Parameters / inputs

- A -sparse real double matrix

## Outputs

- A -sparse real double matrix containing
- unique columns of the input matrix
- This file is a Matlab fallback for the compiled MEX function.

## Implementation structure

- Sparse matrix unique-column utility. Syntax:
- A=spunicols(A)
- A -sparse real double matrix
- A -sparse real double matrix containing
- unique columns of the input matrix
- This file is a Matlab fallback for the compiled MEX function.
- Check consistency
- Return Matlab reference unique columns
- Consistency enforcement
