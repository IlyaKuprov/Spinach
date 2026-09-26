# kernel/operators/hilb2liouv.m

- Signature: `L=hilb2liouv(H,conv_type)`

## Purpose

Converts Hilbert space operators into Liouville space super- operators or state vectors. Syntax: L=hilb2liouv(H,conv_type)

## Physical / mathematical content

- Operator-construction utilities. They build bases and irreducible tensor representations for spin, bosonic, and transition operators.

## Numerical / algorithmic content

## Parameters / inputs

- H -a Hilbert space operator
- conv_type -the type of Liouville space superoperator
- to return:
- 'left' -left side product
- superoperator
- 'right' -right side product
- superoperator
- 'comm' -commutation superoperator
- 'acomm' -anticommutation
- superoperator
- 'statevec' -stretches the operator
- into a state vector

## Outputs

- L -the resulting superoperator or state vector

## Implementation structure

- Converts Hilbert space operators into Liouville space super-
- operators or state vectors. Syntax:
- L=hilb2liouv(H,conv_type)
- H -a Hilbert space operator
- conv_type -the type of Liouville space superoperator
- to return:
- 'left' -left side product
- superoperator
- 'right' -right side product
- 'comm' -commutation superoperator
- 'acomm' -anticommutation
- 'statevec' -stretches the operator
