# kernel/carrier.m

- Signature: `H=carrier(spin_system,spins,operator_type)`

## Purpose

Returns the "carrier" Hamiltonian -the part of the Zeeman interaction Hamiltonian that corresponds to all particles having the Zeeman frequ- ency prescribed by their isotropic free-particle magnetogyric ratio and Z axis magnet field specified by the user. This Hamiltonian is used in rotating frame transforms and average Hamiltonian theories. Syntax: H=carrier(spin_system,spins,operator_type)

## Physical / mathematical content

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

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
