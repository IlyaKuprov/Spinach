# kernel/operators/centrans.m

- Signature: `A=centrans(mult,type)`

## Purpose

Central transition operators of half-integer spins in the Pauli basis. Syntax: A=centrans(mult,type)

## Physical / mathematical content

- Operator-construction utilities. They build bases and irreducible tensor representations for spin, bosonic, and transition operators.

## Numerical / algorithmic content

## Parameters / inputs

- mult -multipicity of the spin in question, an
- even positive integer
- type -operator type: 'z' for polarisation, '+'
- for raising, '-' for lowering

## Outputs

- A -central transition operator

## Implementation structure

- Central transition operators of half-integer spins in the
- Pauli basis. Syntax:
- A=centrans(mult,type)
- mult -multipicity of the spin in question, an
- even positive integer
- type -operator type: 'z' for polarisation, '+'
- for raising, '-' for lowering
- A -central transition operator
- Check consistency
- Build CT operator
- Sx on central transition
- Sy on central transition
