# kernel/multiprop.m

- Signature: `rho=multiprop(spin_system,P,rho,N)`

## Purpose

Applies a propagator repeatedly by binary adaptive squaring. Syntax: rho=multiprop(spin_system,P,rho,N)

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- spin_system -Spinach spin system object
- P -propagator matrix
- rho -state vector or state-vector stack in Liouville space or
- wavefunction formalism, or a density matrix in Hilbert space
- formalism
- N -number of times to apply the propagator

## Outputs

- rho -state vector or density matrix after N applications of P
- Note: the algorithm expands N into binary powers, squares P successively,
- and applies only the active powers to rho. This avoids constructing
- P^N explicitly. Propagator squares are cleaned up using
- spin_system.tols.prop_chop.

## Implementation structure

- Applies a propagator repeatedly by binary adaptive squaring. Syntax:
- rho=multiprop(spin_system,P,rho,N)
- spin_system -Spinach spin system object
- P -propagator matrix
- rho -state vector or state-vector stack in Liouville space or
- wavefunction formalism, or a density matrix in Hilbert space
- formalism
- N -number of times to apply the propagator
- rho -state vector or density matrix after N applications of P
- Note: the algorithm expands N into binary powers, squares P successively,
- and applies only the active powers to rho. This avoids constructing
- P^N explicitly. Propagator squares are cleaned up using
