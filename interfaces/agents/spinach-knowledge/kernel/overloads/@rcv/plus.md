# kernel/overloads/@rcv/plus.m

- Signature: `C=plus(A,B)`

## Purpose

Adds things to RCV sparse matrices. Syntax: C=plus(A,B)

## Physical / mathematical content

- RCV sparse-matrix storage utilities. The focus is data structure design for sparse linear algebra and low-overhead composition of large matrices.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

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
