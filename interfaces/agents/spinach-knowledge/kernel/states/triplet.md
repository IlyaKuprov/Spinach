# kernel/states/triplet.m

- Signature: `[TU,T0,TD]=triplet(spin_system,spin_a,spin_b)`

## Purpose

Constructs the three triplet projectors for two distinct spin-1/2 particles. The spin indices must be distinct positive integers within the spin system, and both selected spins must have multiplicity 2. Pair operators are built in the caller's `spin_a,spin_b` order.

## Construction and normalisation

The function obtains identity, single-spin z, and pairwise Cartesian operators from `state`, then forms

- `TU=EE/4+(ZE+EZ)/2+ZZ` (up projection)
- `T0=EE/4+XX+YY-ZZ` (middle projection)
- `TD=EE/4-(ZE+EZ)/2+ZZ` (down projection)

In Hilbert space each selected two-spin factor is a unit-trace triplet projector; the function applies no further rescaling. Any other spins retain their identity factors. In Liouville space each output is a state vector representing the corresponding operator in the configured basis, not a wavefunction.

## Outputs

- `TU`, `T0`, and `TD` are returned in that order: up, middle, and down. Each is a density matrix in Hilbert-space formalism or a state vector in Liouville-space formalism.

- Source: https://github.com/IlyaKuprov/Spinach/blob/main/kernel/states/triplet.m
- Wiki: https://spindynamics.org/wiki/index.php?title=triplet.m
- Related: [state](../state.md)
