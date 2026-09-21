# kernel/overloads/@rcv/rcv.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@rcv/rcv.m`
- Signature: `obj=rcv(varargin)`
- Total lines: 190

## Purpose

Creates an RCV (row-column-value storage) sparse matrix. Syntax: obj=rcv(M) obj=rcv(dim1,dim2) obj=rcv(R,C,V,dim1,dim2)

## Physical / mathematical content

- RCV sparse-matrix storage utilities. The focus is data structure design for sparse linear algebra and low-overhead composition of large matrices.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `isnumeric()`, `ismatrix()`, `isfloat()`, `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- M -a Matlab matrix
- dim1 -number of rows
- dim2 -number of columns
- R -row indices of non-zero entries
- C -column indices of non-zero entries
- V -values corresponding to entries in R and C

## Outputs

- obj -an RCV sparse matrix object

## Implementation structure

- Creates an RCV (row-column-value storage) sparse matrix. Syntax:
- obj=rcv(M)
- obj=rcv(dim1,dim2)
- obj=rcv(R,C,V,dim1,dim2)
- M -a Matlab matrix
- dim1 -number of rows
- dim2 -number of columns
- R -row indices of non-zero entries
- C -column indices of non-zero entries
- V -values corresponding to entries in R and C
- obj -an RCV sparse matrix object
- Check consistency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `classdef()`, `int64()`, `double()`, `grumble()`, `ismatrix()`, `row()`, `col()`, `val()`, `gpuArray()`, `true()`, `isfloat()`, `isscalar()`, `any()`.
