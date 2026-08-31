# kernel/states/unit_state.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/states/unit_state.m`
- Signature: `rho=unit_state(spin_system)`
- Total lines: 67

## Purpose

Returns the unit state vector or matrix in the current formalism and basis. Syntax: rho=unit_state(spin_system)

## Physical / mathematical content

- State-construction utilities. These routines build equilibrium states, singlets, triplets, partner-state expansions, and physically meaningful density operators in the active basis.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 23-24: Check consistency; implemented by `grumble(spin_system)`.
- Lines 26-27: Decide how to proceed; implemented by `switch spin_system.bas.formalism`.
- Lines 31-32: Unit population of T(0,0) state; implemented by `rho=sparse(1,1,1,size(spin_system.bas.basis,1),1)`.
- Lines 36-37: Normalized stretched unit matrix; implemented by `rho=speye(prod(spin_system.comp.mults))`.
- Lines 42-43: Sparse unit matrix; implemented by `rho=speye(prod(spin_system.comp.mults))`.
- Lines 47-48: Complain and bomb out; implemented by `error('unknown formalism specification.')`.

### Control flow inferred from the code

- Line 27: dispatches on `spin_system.bas.formalism`; cases `'sphten-liouv'`, `'zeeman-liouv'`, `'zeeman-hilb'`.

### Key state/data transformations

- Lines 32: computes `rho` using `rho=sparse(1,1,1,size(spin_system.bas.basis,1),1)`.

### Local helper functions

- Line 55: `grumble()` — `function grumble(spin_system)`. There used to be a simple story about Russian literature, that we thought the good writers were the ones who opposed the regime. Once
  - Representative operation: `if (~isfield(spin_system,'bas'))||(~isfield(spin_system.bas,'formalism'))`.
  - Representative operation: `error('the spin_system object does not contain the required information.')`.

## Parameters / inputs

- spin_system -Spinach data object containing basis
- information (call basis.m first)

## Outputs

- rho -vector or matrix representation of
- the unit state

## Implementation structure

- Returns the unit state vector or matrix in the current formalism
- and basis. Syntax:
- rho=unit_state(spin_system)
- spin_system -Spinach data object containing basis
- information (call basis.m first)
- rho -vector or matrix representation of
- the unit state
- Check consistency
- Decide how to proceed
- Unit population of T(0,0) state
- Normalized stretched unit matrix
- Sparse unit matrix

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `speye()`, `rho()`, `isfield()`.
