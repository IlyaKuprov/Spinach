# kernel/operators/ct2ist.m

- Signature: `[states,coeffs]=ct2ist(mult,type)`

## Purpose

Irreducible spherical tensor expansion of central transition operators of half-integer spins. Syntax: [states,coeffs]=ct2ist(mult,type)

## Physical / mathematical content

- Operator-construction utilities. They build bases and irreducible tensor representations for spin, bosonic, and transition operators.

## Numerical / algorithmic content

## Parameters / inputs

- mult -multipicity of the spin in question, an
- even positive integer
- type -operator type: 'z' for polarisation, '+'
- for raising, '-' for lowering

## Outputs

- states -states, in the Spinach IST basis index-
- ing, that contribute to the operator in
- question; use lin2lm to convert to L,M
- spherical tensor indices
- coeffs -coefficients with which the ISTs enter
- the linear combination

## Implementation structure

- Irreducible spherical tensor expansion of central transition
- operators of half-integer spins. Syntax:
- [states,coeffs]=ct2ist(mult,type)
- mult -multipicity of the spin in question, an
- even positive integer
- type -operator type: 'z' for polarisation, '+'
- for raising, '-' for lowering
- states -states, in the Spinach IST basis index-
- ing, that contribute to the operator in
- question; use lin2lm to convert to L,M
- spherical tensor indices
- coeffs -coefficients with which the ISTs enter
