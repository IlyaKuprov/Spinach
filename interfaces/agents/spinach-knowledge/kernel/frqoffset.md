# kernel/frqoffset.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/frqoffset.m`
- Signature: `H=frqoffset(spin_system,H,parameters)`
- Total lines: 108

## Purpose

Adds omega*Lz Larmor frequency offsets to the Hamiltonian; this is useful in liquid state NMR experiments. Syntax: H=frqoffset(spin_system,H,parameters)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `numel()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 34-35: Check consistency; implemented by `grumble(spin_system,parameters)`.
- Lines 37-38: See if there are multiple channels on the same spin; implemented by `[unique_spins,forward_index,backward_index]=unique(parameters.spins)`.
- Lines 40-41: Decide how to proceed; implemented by `if numel(unique_spins)==numel(parameters.spins)`.
- Lines 43-44: Simply apply the offsets; implemented by `for n=find(parameters.offset~=0)`.
- Lines 52-53: Get unique offsets; implemented by `unique_offsets=parameters.offset(forward_index)`.
- Lines 55-56: Check offsets on duplicate channels; implemented by `if ~all(unique_offsets(backward_index)==parameters.offset)`.
- Lines 60-61: Apply the offsets; implemented by `for n=find(unique_offsets~=0)`.

### Control flow inferred from the code

- Line 41: conditional branch on `numel(unique_spins)==numel(parameters.spins)`.
- Line 44: `for` loop over `n=find(parameters.offset~=0)`.
- Line 56: conditional branch on `~all(unique_offsets(backward_index)==parameters.offset)`.
- Line 61: `for` loop over `n=find(unique_offsets~=0)`.

### Key state/data transformations

- Lines 38: computes `[unique_spins,forward_index,backward_index]` using `[unique_spins,forward_index,backward_index]=unique(parameters.spins)`.
- Lines 47: computes `H` using `H=H+2*pi*parameters.offset(n)*operator(spin_system,'Lz',parameters.spins{n})`.
- Lines 53: computes `unique_offsets` using `unique_offsets=parameters.offset(forward_index)`.

### Local helper functions

- Line 72: `grumble()` — `function grumble(spin_system,parameters)`.
  - Representative operation: `if ~isfield(parameters,'spins')`.
  - Representative operation: `error('parameters.spins variable must be present.')`.

## Parameters / inputs

- H -Hamiltonian operator or commutati-
- on superoperator
- parameters.spins -a cell array giving the spins that
- the offsets should be applied to,
- e.g. {'1H','13C'}
- parameters.offset -a vector of offsets (in Hz) on
- each of the spins listed in the
- parameters.spins array

## Outputs

- H -Hamiltonian operator or commutati-
- on superoperator
- Note: offset transformation of this kind is an approximati-
- on, use rotframe.m or intrep.m if a rigorous treat-
- ment of second order effects is required.

## Implementation structure

- Adds omega*Lz Larmor frequency offsets to the Hamiltonian;
- this is useful in liquid state NMR experiments. Syntax:
- H=frqoffset(spin_system,H,parameters)
- H -Hamiltonian operator or commutati-
- on superoperator
- parameters.spins -a cell array giving the spins that
- the offsets should be applied to,
- e.g. {'1H','13C'}
- parameters.offset -a vector of offsets (in Hz) on
- each of the spins listed in the
- parameters.spins array
- Note: offset transformation of this kind is an approximati-

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `report()`, `num2str()`, `operator()`, `all()`, `unique_offsets()`, `isfield()`, `iscell()`, `cellfun()`, `any()`, `ismember()`, `elseif()`, `isvector()`.
