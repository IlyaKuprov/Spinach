# kernel/operators/hilb2liouv.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/operators/hilb2liouv.m`
- Signature: `L=hilb2liouv(H,conv_type)`
- Total lines: 94

## Purpose

Converts Hilbert space operators into Liouville space super- operators or state vectors. Syntax: L=hilb2liouv(H,conv_type)

## Physical / mathematical content

- Operator-construction utilities. They build bases and irreducible tensor representations for spin, bosonic, and transition operators.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `speye()`, `transpose()`, `ischar()`.
