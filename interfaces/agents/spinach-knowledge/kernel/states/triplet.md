# kernel/states/triplet.m

- Signature: `[TU,T0,TD]=triplet(spin_system,spin_a,spin_b)`

## Purpose

Returns the components of the two-spin triplet state; both particles must be spin-1/2. Syntax: [Tp,T0,Tm]=triplet(spin_system,spin_a,spin_b)

## Physical / mathematical content

- State-construction utilities. These routines build equilibrium states, singlets, triplets, partner-state expansions, and physically meaningful density operators in the active basis.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Parameters / inputs

- spin_a -the number of the first spin in the
- triplet state
- spin_b -the number of the second spin in the
- triplet state

## Outputs

- TU,T0,TD -density matrices (Hilbert space) or
- state vectors (Liouville space) of
- TU, T0, and TD projections

## Implementation structure

- Returns the components of the two-spin triplet state; both particles
- must be spin-1/2. Syntax:
- [Tp,T0,Tm]=triplet(spin_system,spin_a,spin_b)
- spin_a -the number of the first spin in the
- triplet state
- spin_b -the number of the second spin in the
- TU,T0,TD -density matrices (Hilbert space) or
- state vectors (Liouville space) of
- TU, T0, and TD projections
- Check consistency
- Build the component operators
- Build the triplet states
