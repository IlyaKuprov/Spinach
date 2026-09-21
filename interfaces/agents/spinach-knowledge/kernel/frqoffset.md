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
