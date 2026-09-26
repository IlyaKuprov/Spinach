# kernel/operators/oper2ist.m

- Signature: `[states,coeffs]=oper2ist(A)`

## Purpose

Irreducible spherical tensor operator expansion of a user- specified square matrix. Syntax: [states,coeffs]=oper2ist(A)

## Physical / mathematical content

- Operator-construction utilities. They build bases and irreducible tensor representations for spin, bosonic, and transition operators.

## Numerical / algorithmic content

## Parameters / inputs

- A -a square matrix

## Outputs

- states -states, in the Spinach IST basis index-
- ing convention, that contribute to the
- operator in question; use lin2lm() fun-
- ction to convert to L,M spherical tens-
- or indices
- coeffs -coefficients with which the ISTs enter
- the linear combination

## Implementation structure

- Irreducible spherical tensor operator expansion of a user-
- specified square matrix. Syntax:
- [states,coeffs]=oper2ist(A)
- A -a square matrix
- states -states, in the Spinach IST basis index-
- ing convention, that contribute to the
- operator in question; use lin2lm() fun-
- ction to convert to L,M spherical tens-
- or indices
- coeffs -coefficients with which the ISTs enter
- the linear combination
- Check consistency
