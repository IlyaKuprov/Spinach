# kernel/states/four_spin_states.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/states/four_spin_states.m`
- Signature: `rho=four_spin_states(spin_system,spins,spin_state)`
- Total lines: 194

## Purpose

Returns user-specified states for a system of four spin-1/2 particles; see also the enclosed Mathematica file. Syntax: rho=four_spin_states(spin_system,spins,spin_state)

## Physical / mathematical content

- State-construction utilities. These routines build equilibrium states, singlets, triplets, partner-state expansions, and physically meaningful density operators in the active basis.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- spins -a row vector of four spin numbers
- spin_state -one of the possible singlet-triplet
- product states, see below

## Outputs

- rho -a density matrix (Hilbert space) or
- a state vector (Liouville space)

## Implementation structure

- Returns user-specified states for a system of four spin-1/2
- particles; see also the enclosed Mathematica file. Syntax:
- rho=four_spin_states(spin_system,spins,spin_state)
- spins -a row vector of four spin numbers
- spin_state -one of the possible singlet-triplet
- product states, see below
- rho -a density matrix (Hilbert space) or
- a state vector (Liouville space)
- Check consistency
- Component operators: four-spin
- Component operators: three-spin
- Component operators: two-spin

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `state()`, `num2cell()`, `spins()`, `isvector()`, `any()`, `isrow()`, `ischar()`.
