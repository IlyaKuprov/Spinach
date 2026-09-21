# kernel/overloads/@rcv/times.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@rcv/times.m`
- Signature: `C=times(A,B)`
- Total lines: 57

## Purpose

Multiplies an RCV sparse matrix by a numeric scalar, in either operand order. Syntax: C=times(A,B)

## Physical / mathematical content

- RCV sparse-matrix storage utilities. The focus is data structure design for sparse linear algebra and low-overhead composition of large matrices.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- A,B -an RCV sparse matrix and a numeric
- scalar, in either order

## Outputs

- C -RCV sparse matrix

## Implementation structure

- Multiplies an RCV sparse matrix by a numeric scalar,
- in either operand order. Syntax:
- C=times(A,B)
- A,B -an RCV sparse matrix and a numeric
- scalar, in either order
- C -RCV sparse matrix
- Check consistency
- RCV sparse by a scalar
- Scalar by RCV sparse
- Consistency enforcement
- They say that the fish that gets away
- looks bigger than it really is.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `xor()`, `isscalar()`.
