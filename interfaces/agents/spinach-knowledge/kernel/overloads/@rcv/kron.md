# kernel/overloads/@rcv/kron.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@rcv/kron.m`
- Signature: `C=kron(A,B)`
- Total lines: 58

## Purpose

Kronecker product between two RCV sparse matrices. Syntax: C=kron(A,B)

## Physical / mathematical content

- RCV sparse-matrix storage utilities. The focus is data structure design for sparse linear algebra and low-overhead composition of large matrices.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- A -left RCV sparse matrix
- B -right RCV sparse matrix

## Outputs

- C -RCV sparse matrix

## Implementation structure

- Kronecker product between two RCV sparse matrices. Syntax:
- C=kron(A,B)
- A -left RCV sparse matrix
- B -right RCV sparse matrix
- C -RCV sparse matrix
- Check consistency
- Compute the output dimensions
- Align locations
- Build the Cartesian product of indices and values
- Assemble the output RCV object
- Consistency enforcement
- Along the Yangzi River, apes moan ceaselessly.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `gpuArray()`, `rcv()`.
