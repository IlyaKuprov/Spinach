# kernel/operators/twospinist.m

- Signature: `T=twospinist(spin_system,spin_a,spin_b,indices,type)`

## Purpose

Two-spin irreducible spherical tensor operators. Syntax: T=twospinist(spin_system,spin_a,spin_b,indices,type)

## Physical / mathematical content

- Operator-construction utilities. They build bases and irreducible tensor representations for spin, bosonic, and transition operators.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Parameters / inputs

- spin_a -number of the first spin
- spin_b -number of the second spin
- indices -two-element vector [L,M] containing
- rank and projection index; L=1,2 are
- available
- In Liouville space, type can be set to:
- 'left' -produces left side product superoperator
- 'right' -produces right side product superoperator
- 'comm' -produces commutation superoperator (default)
- 'acomm' -produces anticommutation superoperator
- In Hilbert space calculations, the type parameter is ig-
- nored, and the operator itself is always returned.

## Outputs

- T -irreducible spherical tensor operator

## Implementation structure

- Two-spin irreducible spherical tensor operators. Syntax:
- T=twospinist(spin_system,spin_a,spin_b,indices,type)
- spin_a -number of the first spin
- spin_b -number of the second spin
- indices -two-element vector [L,M] containing
- rank and projection index; L=1,2 are
- available
- In Liouville space, type can be set to:
- 'left' -produces left side product superoperator
- 'right' -produces right side product superoperator
- 'comm' -produces commutation superoperator (default)
- 'acomm' -produces anticommutation superoperator
