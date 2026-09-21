# kernel/operators/ct2ist.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/operators/ct2ist.m`
- Signature: `[states,coeffs]=ct2ist(mult,type)`
- Total lines: 55

## Purpose

Irreducible spherical tensor expansion of central transition operators of half-integer spins. Syntax: [states,coeffs]=ct2ist(mult,type)

## Physical / mathematical content

- Operator-construction utilities. They build bases and irreducible tensor representations for spin, bosonic, and transition operators.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `centrans()`, `oper2ist()`, `isscalar()`, `ischar()`, `ismember()`.
