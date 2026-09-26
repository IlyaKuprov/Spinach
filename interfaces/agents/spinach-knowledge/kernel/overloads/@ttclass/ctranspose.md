# kernel/overloads/@ttclass/ctranspose.m

- Signature: `ttrain=ctranspose(ttrain)`

## Purpose

Computes a Hermitian conjugate of a matrix in a tensor train representation. Syntax: ttrain=ctranspose(ttrain)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content

## Parameters / inputs

- ttrain -tensor train representation of a matrix

## Outputs

- ttrain -Hermitian conjugate of the input tensor train

## Implementation structure

- Computes a Hermitian conjugate of a matrix in a tensor train
- representation. Syntax:
- ttrain=ctranspose(ttrain)
- ttrain -tensor train representation of a matrix
- ttrain -Hermitian conjugate of the input tensor train
- Read tensor sizes and ranks
- Swap the middle dimensions of all cores
- Conjugate the result
- What gives the artist real prestige is his imitators.
- Igor Stravinsky
- #NGRUM
