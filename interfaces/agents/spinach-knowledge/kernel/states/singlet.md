# kernel/states/singlet.m

- Signature: `S=singlet(spin_system,spin_a,spin_b)`

## Purpose

Returns a two-spin singlet state; both particles must be spin-1/2. Syntax: `S=singlet(spin_system,spin_a,spin_b)`

## Physical / mathematical content

- The function constructs the singlet operator `S=EE/4-(XX+YY+ZZ)` from the two spins' identity and Cartesian spin operators.

## Parameters / inputs

- spin_a -the number of the first spin in the
- singlet state
- spin_b -the number of the second spin in the
- singlet state

## Outputs

- S -a density matrix (Hilbert space) or
- a state vector (Liouville space)
