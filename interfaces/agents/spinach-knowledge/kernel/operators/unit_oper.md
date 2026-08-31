# kernel/operators/unit_oper.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/operators/unit_oper.m`
- Signature: `A=unit_oper(spin_system)`
- Total lines: 74

## Purpose

Returns a unit operator in the current formalism and basis. The operator has dimension equal to the basis size in sphten-liouv formalism, the dimension equal to the product of all spin multi- plicities in zeeman-hilb and zeeman-wavef formalisms, and the dimension of square of the product of all spin multiplicities in zeeman-liouv formalism. Syntax: A=unit_oper(spin_system)

## Physical / mathematical content

- Operator-construction utilities. They build bases and irreducible tensor representations for spin, bosonic, and transition operators.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 27-28: Check consistency; implemented by `grumble(spin_system)`.
- Lines 30-31: Decide how to proceed; implemented by `switch spin_system.bas.formalism`.
- Lines 35-36: Unit matrix; implemented by `A=speye(size(spin_system.bas.basis,1))`.
- Lines 40-41: Unit matrix; implemented by `A=speye(prod(spin_system.comp.mults.^2))`.
- Lines 45-46: Unit matrix; implemented by `A=speye(prod(spin_system.comp.mults))`.
- Lines 50-51: Complain and bomb out; implemented by `error('unknown formalism specification.')`.

### Control flow inferred from the code

- Line 31: dispatches on `spin_system.bas.formalism`; cases `'sphten-liouv'`, `'zeeman-liouv'`, `{'zeeman-hilb','zeeman-wavef'}`.

### Key state/data transformations

- Lines 36: computes `A` using `A=speye(size(spin_system.bas.basis,1))`.

### Local helper functions

- Line 58: `grumble()` — `function grumble(spin_system)`. The substance of this book, as it is expressed in the editor's preface, is that to measure "right" by the false philosophy of the Hebrew prophets and
  - Representative operation: `if (~isfield(spin_system,'bas'))||(~isfield(spin_system.bas,'formalism'))`.
  - Representative operation: `error('the spin_system object does not contain the required information.')`.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `speye()`, `isfield()`.
