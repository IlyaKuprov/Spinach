# kernel/states/triplet.m

- Signature: `[TU,T0,TD]=triplet(spin_system,spin_a,spin_b)`

## Purpose

Returns the three components of the two-spin triplet state for two distinct spin-1/2 particles. Syntax: `[TU,T0,TD]=triplet(spin_system,spin_a,spin_b)`.

## Physical / mathematical content

The routine constructs the triplet projectors from identity and spin operators on the selected pair: `TU=EE/4+(ZE+EZ)/2+ZZ`, `T0=EE/4+XX+YY-ZZ`, and `TD=EE/4-(ZE+EZ)/2+ZZ`. Here the component operators are formed using `state` with `E`, `Lx`, `Ly`, and `Lz`.

## Numerical / algorithmic content

Checks that the spin indices are distinct positive integers within the system and that both selected spins have multiplicity 2, then builds the three component operators and triplet states.

## Parameters / inputs

- `spin_a` - number of the first spin in the triplet state.
- `spin_b` - number of the second spin in the triplet state.

## Outputs

- `TU`, `T0`, `TD` - density matrices in Hilbert space or state vectors in Liouville space for the three triplet projections.

## Implementation structure

- Enforces input consistency, constructs the pair operators, then returns the up, middle, and down triplet components.
