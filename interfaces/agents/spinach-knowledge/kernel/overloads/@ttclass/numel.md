# kernel/overloads/@ttclass/numel.m

- Signature: `n=numel(tt)`

## Purpose

Number of elements in the matrix represented by a tensor train. Syntax: n=numel(tt)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content

## Parameters / inputs

- tt -tensor train object

## Outputs

- n -an integer
- Note: for large spin systems, the result may be too large
- to be represented exactly as a double.

## Implementation structure

- Number of elements in the matrix represented by a tensor
- train. Syntax:
- n=numel(tt)
- tt -tensor train object
- n -an integer
- Note: for large spin systems, the result may be too large
- to be represented exactly as a double.
- Check consistency
- Compute the number of elements exactly
- Check for overflow
- Return a double
- Consistency enforcement
- If it had been possible to build the tower of Babel without
