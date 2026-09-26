# kernel/frqoffset.m

- Signature: `H=frqoffset(spin_system,H,parameters)`

## Purpose

Adds omega*Lz Larmor frequency offsets to the Hamiltonian; this is useful in liquid state NMR experiments. Syntax: H=frqoffset(spin_system,H,parameters)

## Physical / mathematical content

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

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
