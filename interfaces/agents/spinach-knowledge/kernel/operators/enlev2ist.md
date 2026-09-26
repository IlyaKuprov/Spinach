# kernel/operators/enlev2ist.m

- Signature: `[states,coeffs]=enlev2ist(mult,lvl_num,particle)`

## Purpose

Irreducible spherical tensor expansion of specific Zeeman energy level projectors. Syntax: [states,coeffs]=enlev2ist(mult,lvl_num,particle)

## Physical / mathematical content

- Operator-construction utilities. They build bases and irreducible tensor representations for spin, bosonic, and transition operators.

## Numerical / algorithmic content

## Parameters / inputs

- mult -multipicity of the spin in question, a
- positive integer
- lvl_num -energy level number, counting from the
- bottom up for spins and from top down
- for bosons
- particle -particle type, 'S' for a spin and 'B'
- for a boson

## Outputs

- states -states, in the Spinach IST basis index-
- ing, that contribute to the operator in
- question; use lin2lm to convert to L,M
- spherical tensor indices
- coeffs -coefficients with which the ISTs enter
- the linear combination

## Implementation structure

- Irreducible spherical tensor expansion of specific Zeeman
- energy level projectors. Syntax:
- [states,coeffs]=enlev2ist(mult,lvl_num,particle)
- mult -multipicity of the spin in question, a
- positive integer
- lvl_num -energy level number, counting from the
- bottom up for spins and from top down
- for bosons
- particle -particle type, 'S' for a spin and 'B'
- for a boson
- states -states, in the Spinach IST basis index-
- ing, that contribute to the operator in
