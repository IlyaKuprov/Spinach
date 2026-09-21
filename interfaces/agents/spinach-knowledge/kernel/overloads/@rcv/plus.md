# kernel/overloads/@rcv/plus.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@rcv/plus.m`
- Signature: `C=plus(A,B)`
- Total lines: 98

## Purpose

Adds things to RCV sparse matrices. Syntax: C=plus(A,B)

## Physical / mathematical content

- RCV sparse-matrix storage utilities. The focus is data structure design for sparse linear algebra and low-overhead composition of large matrices.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- A -left operand
- B -right operand

## Outputs

- C -sum A+B, RCV sparse matrix

## Implementation structure

- Adds things to RCV sparse matrices. Syntax:
- C=plus(A,B)
- A -left operand
- B -right operand
- C -sum A+B, RCV sparse matrix
- Check consistency
- Process the special cases
- Explain the refusal to add a scalar to an RCV object
- Add a scalar to an RCV sparse matrix
- Add two RCV sparse matrices
- Check for dimension match
- Align locations

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `isscalar()`, `gpuArray()`, `issparse()`, `rcv()`.
