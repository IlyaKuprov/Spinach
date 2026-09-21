# kernel/overloads/@rcv/size.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@rcv/size.m`
- Signature: `[s,ncols]=size(A,dim)`
- Total lines: 61

## Purpose

Returns the size of an RCV sparse matrix. Syntax: s=size(A,dim) [s,ncols]=size(A)

## Physical / mathematical content

- RCV sparse-matrix storage utilities. The focus is data structure design for sparse linear algebra and low-overhead composition of large matrices.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- A -RCV sparse matrix
- dim -optional dimension index

## Outputs

- s -size vector, dimension length, or number
- of rows in the two-output form
- ncols -number of columns in the two-output form

## Implementation structure

- Returns the size of an RCV sparse matrix. Syntax:
- s=size(A,dim)
- [s,ncols]=size(A)
- A -RCV sparse matrix
- dim -optional dimension index
- s -size vector, dimension length, or number
- of rows in the two-output form
- ncols -number of columns in the two-output form
- Check consistency
- Refuse two outputs with a dimension query
- Mimic Matlab
- Consistency enforcement
- I did not succeed in life by intelligence. I succeeded
- because I have a long attention span.
- Charlie Munger

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `isscalar()`.
