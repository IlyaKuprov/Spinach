# kernel/states/singlet.m

- Signature: `S=singlet(spin_system,spin_a,spin_b)`

## Purpose

Returns a two-spin singlet state; both particles must be spin-1/2. Syntax: rho=singlet(spin_system,spin_a,spin_b)

## Physical / mathematical content

- State-construction utilities. These routines build equilibrium states, singlets, triplets, partner-state expansions, and physically meaningful density operators in the active basis.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Parameters / inputs

- spin_a -the number of the first spin in the
- singlet state
- spin_b -the number of the second spin in the
- singlet state

## Outputs

- S -a density matrix (Hilbert space) or
- a state vector (Liouville space)

## Implementation structure

- Returns a two-spin singlet state; both particles must be
- spin-1/2. Syntax:
- rho=singlet(spin_system,spin_a,spin_b)
- spin_a -the number of the first spin in the
- singlet state
- spin_b -the number of the second spin in the
- S -a density matrix (Hilbert space) or
- a state vector (Liouville space)
- Check consistency
- Build the component operators
- Build the singlet state
- Consistency enforcement
