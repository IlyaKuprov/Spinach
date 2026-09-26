# kernel/operators/oper2bm.m

- Signature: `[states,coeffs]=oper2bm(A)`

## Purpose

Bosonic monomial operator expansion of a user-specified square matrix. Syntax: [states,coeffs]=oper2bm(A)

## Physical / mathematical content

- Operator-construction utilities. They build bases and irreducible tensor representations for spin, bosonic, and transition operators.

## Numerical / algorithmic content

## Parameters / inputs

- A -a square matrix

## Outputs

- states -states, in the Spinach BM basis index-
- ing convention, that contribute to the
- operator in question; use lin2kq() fun-
- ction to convert to K,Q bosonic monomi-
- al indices
- coeffs -coefficients with which the BMs enter
- the linear combination

## Implementation structure

- Bosonic monomial operator expansion of a user-specified
- square matrix. Syntax:
- [states,coeffs]=oper2bm(A)
- A -a square matrix
- states -states, in the Spinach BM basis index-
- ing convention, that contribute to the
- operator in question; use lin2kq() fun-
- ction to convert to K,Q bosonic monomi-
- al indices
- coeffs -coefficients with which the BMs enter
- the linear combination
- Check consistency
