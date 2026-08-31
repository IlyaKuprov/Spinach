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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 37-38: Check consistency; implemented by `grumble(H,conv_type)`.
- Lines 40-41: Prepare a unit matrix; implemented by `unit=speye(size(H))`.
- Lines 43-44: Decide how to proceed; implemented by `switch conv_type`.
- Lines 48-49: Compute a commutation superoperator; implemented by `L=kron(unit,H)-kron(transpose(H),unit)`.
- Lines 53-54: Compute an anticommutation superoperator; implemented by `L=kron(unit,H)+kron(transpose(H),unit)`.
- Lines 58-59: Compute a left side product superoperator; implemented by `L=kron(unit,H)`.
- Lines 63-64: Compute a right side product superoperator; implemented by `L=kron(transpose(H),unit)`.
- Lines 68-69: Stretch into a state vector; implemented by `L=H(:)`.
- Lines 73-74: Complain and bomb out; implemented by `error('unknown conversion type specification.')`.

### Control flow inferred from the code

- Line 44: dispatches on `conv_type`; cases `'comm'`, `'acomm'`, `'left'`, `'right'`, `'statevec'`.

### Key state/data transformations

- Lines 41: computes `unit` using `unit=speye(size(H))`.
- Lines 49: computes `L` using `L=kron(unit,H)-kron(transpose(H),unit)`.

### Local helper functions

- Line 81: `grumble()` — `function grumble(H,conv_type)`. If it wasn't for the fun and money, I really don't know why I'd bother.
  - Representative operation: `if ~isnumeric(H)`.
  - Representative operation: `error('H parameter must be numeric.')`.

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
