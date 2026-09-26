# kernel/overloads/@ttclass/conj.m

- Signature: `tt=conj(tt)`

## Purpose

Conjugates all core elements and coefficients of a tensor train object. Syntax: tt=conj(tt)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content

## Parameters / inputs

- tt -tensor train object

## Outputs

- tt -tensor train object with complex-conjugated cores
- and coefficients

## Implementation structure

- Conjugates all core elements and coefficients of a tensor
- train object. Syntax:
- tt=conj(tt)
- tt -tensor train object
- tt -tensor train object with complex-conjugated cores
- and coefficients
- Read tensor train sizes and ranks
- Conjugate the cores
- Conjugate the coefficients
- I asked God for a bike, but I know God
- doesn't work that way. So I stole a bike
- and asked for forgiveness.
