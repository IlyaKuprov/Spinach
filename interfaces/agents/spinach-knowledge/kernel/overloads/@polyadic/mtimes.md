# kernel/overloads/@polyadic/mtimes.m

- Signature: `C=mtimes(A,B)`

## Purpose

Performs multiplications involving polyadics. Syntax: C=mtimes(A,B)

## Physical / mathematical content

Multiplication composes the polyadic terms with scalar or matrix operands; compatible small single-core polyadics can be multiplied core by core, while other matrix actions are retained in prefix or suffix buffers.

## Numerical / algorithmic content

## Parameters / inputs

- A, B: a polyadic object or a numerical array

## Outputs

- C: a polyadic object or a numerical array

## Implementation structure

- Performs multiplications involving polyadics. Syntax:
- C=mtimes(A,B)
- A,B -a polyadic or a numerical array
- C -a polyadic or a numerical array
- When A is a number
- Multiply smallest cores in the B buffer
- When A is a sparse matrix
- Attach as a prefix to B
- When A is a full matrix
- Issue a recursive call
- When B is a number
- Multiply smallest cores in the A buffer
