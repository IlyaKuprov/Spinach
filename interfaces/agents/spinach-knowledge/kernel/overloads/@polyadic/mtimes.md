# kernel/overloads/@polyadic/mtimes.m

- Signature: `C=mtimes(A,B)`

## Purpose

Performs multiplications involving polyadics. Syntax: C=mtimes(A,B)

## Physical / mathematical content

- Polyadic tensor-product linear algebra. The emphasis is compressed operator representation, deferred algebra, and efficient Kronecker-structured manipulations.

## Numerical / algorithmic content

## Parameters / inputs

- A,B -a polyadic or a numerical array

## Outputs

- C -a polyadic or a numerical array

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
