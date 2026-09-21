# kernel/overloads/@rcv/mtimes.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@rcv/mtimes.m`
- Signature: `C=mtimes(A,B)`
- Total lines: 90

## Purpose

Multiplication for RCV sparse matrices. Syntax: C=mtimes(A,B)

## Physical / mathematical content

- RCV sparse-matrix storage utilities. The focus is data structure design for sparse linear algebra and low-overhead composition of large matrices.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- A -left operand
- B -right operand

## Outputs

- C -product A*B as a Matlab sparse matrix if
- both operands are RCV or Matlab matrices

## Implementation structure

- Multiplication for RCV sparse matrices. Syntax:
- C=mtimes(A,B)
- A -left operand
- B -right operand
- C -product A*B as a Matlab sparse matrix if
- both operands are RCV or Matlab matrices
- Check consistency
- RCV sparse by a scalar
- Scalar by RCV sparse
- RCV sparse by RCV sparse
- Check dimension consistency
- Result is Matlab sparse

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `isscalar()`, `issparse()`.
