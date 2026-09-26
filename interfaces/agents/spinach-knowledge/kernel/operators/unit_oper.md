# kernel/operators/unit_oper.m

- Signature: `A=unit_oper(spin_system)`

## Purpose

Returns a unit operator in the current formalism and basis. The operator has dimension equal to the basis size in sphten-liouv formalism, the dimension equal to the product of all spin multi- plicities in zeeman-hilb and zeeman-wavef formalisms, and the dimension of square of the product of all spin multiplicities in zeeman-liouv formalism. Syntax: A=unit_oper(spin_system)

## Physical / mathematical content

- Operator-construction utilities. They build bases and irreducible tensor representations for spin, bosonic, and transition operators.

## Numerical / algorithmic content

## Parameters / inputs

- spin_system -Spinach data object containing basis
- information (call basis.m first)

## Outputs

- A -a sparse unit matrix of appropriate
- dimension

## Implementation structure

- Returns a unit operator in the current formalism and basis. The
- operator has dimension equal to the basis size in sphten-liouv
- formalism, the dimension equal to the product of all spin multi-
- plicities in zeeman-hilb and zeeman-wavef formalisms, and the
- dimension of square of the product of all spin multiplicities in
- zeeman-liouv formalism. Syntax:
- A=unit_oper(spin_system)
- spin_system -Spinach data object containing basis
- information (call basis.m first)
- A -a sparse unit matrix of appropriate
- dimension
- Check consistency
