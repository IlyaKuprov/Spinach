# kernel/overloads/@ttclass/diag.m

- Signature: `tt=diag(tt)`

## Purpose

Mimics the diag behaviour for tensor train matrix. Syntax: tt=diag(tt)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content

## Parameters / inputs

- tt -a tensor train representation of a matrix

## Outputs

- tt -if the input is a square matrix, returns a
- vector by computing diag of every core; if
- the input is a vector (one mode size is
- ones), returns a diagonal matrix

## Implementation structure

- Mimics the diag behaviour for tensor train matrix. Syntax:
- tt=diag(tt)
- tt -a tensor train representation of a matrix
- tt -if the input is a square matrix, returns a
- vector by computing diag of every core; if
- the input is a vector (one mode size is
- ones), returns a diagonal matrix
- Read tensor train sizes and ranks
- Decide the dimensions
- Vector on input, diagonal matrix on output
- Matrix on input, column vector on output
- The only mistake [the famous criminal finacier] Bernie Madoff
