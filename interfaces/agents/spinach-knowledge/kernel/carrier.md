# kernel/carrier.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/carrier.m`
- Signature: `H=carrier(spin_system,spins,operator_type)`
- Total lines: 83

## Purpose

Returns the "carrier" Hamiltonian -the part of the Zeeman interaction Hamiltonian that corresponds to all particles having the Zeeman frequ- ency prescribed by their isotropic free-particle magnetogyric ratio and Z axis magnet field specified by the user. This Hamiltonian is used in rotating frame transforms and average Hamiltonian theories. Syntax: H=carrier(spin_system,spins,operator_type)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 37-38: Set the default for the type; implemented by `if ~exist('operator_type','var'), operator_type='comm'; end`.
- Lines 40-41: Check consistency; implemented by `grumble(spin_system,spins,operator_type)`.
- Lines 43-44: Preallocate the answer; implemented by `H=mprealloc(spin_system,1)`.
- Lines 46-47: Find the spins; implemented by `if strcmp(spins,'all')`.
- Lines 53-54: Compute the answer; implemented by `for n=1:numel(spin_index)`.
- Lines 61-62: Clean up the answer; implemented by `H=clean_up(spin_system,(H+H')/2,spin_system.tols.liouv_zero)`.

### Control flow inferred from the code

- Line 38: conditional branch on `~exist('operator_type','var'), operator_type='comm'; end`.
- Line 47: conditional branch on `strcmp(spins,'all')`.
- Line 54: `for` loop over `n=1:numel(spin_index)`.
- Line 55: conditional branch on `abs(spin_system.inter.basefrqs(spin_index(n)))>0`.

### Key state/data transformations

- Lines 44: computes `H` using `H=mprealloc(spin_system,1)`.
- Lines 48: computes `spin_index` using `spin_index=1:spin_system.comp.nspins`.

### Local helper functions

- Line 67: `grumble()` — `function grumble(spin_system,spins,operator_type)`.
  - Representative operation: `if ~ischar(spins)`.
  - Representative operation: `error('spins argument must be a character array')`.

## Parameters / inputs

- spins -a string specifying the isotope, e.g. '1H';
- to select all spins, use 'all'.
- in Liouville space, operator_type can be set to
- 'left' -produces left side product superoperator
- 'right' -produces right side product superoperator
- 'comm' -produces commutation superoperator (default)
- 'acomm' -produces anticommutation superoperator
- in Hilbert space this parameter is ignored.

## Outputs

- H -a Hamiltonian (Hilbert space) or its superoperator
- of the specified type (Liouville space).

## Implementation structure

- Returns the "carrier" Hamiltonian -the part of the Zeeman interaction
- Hamiltonian that corresponds to all particles having the Zeeman frequ-
- ency prescribed by their isotropic free-particle magnetogyric ratio and
- Z axis magnet field specified by the user. This Hamiltonian is used in
- rotating frame transforms and average Hamiltonian theories. Syntax:
- H=carrier(spin_system,spins,operator_type)
- spins -a string specifying the isotope, e.g. '1H';
- to select all spins, use 'all'.
- in Liouville space, operator_type can be set to
- 'left' -produces left side product superoperator
- 'right' -produces right side product superoperator
- 'comm' -produces commutation superoperator (default)

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `exist()`, `grumble()`, `mprealloc()`, `strcmp()`, `spin_index()`, `operator()`, `clean_up()`, `ischar()`, `ismember()`.
