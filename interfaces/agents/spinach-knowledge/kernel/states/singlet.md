# kernel/states/singlet.m

- Signature: `S=singlet(spin_system,spin_a,spin_b)`

## Purpose

Constructs the singlet projector for two distinct spin-1/2 particles. The spin indices must be distinct positive integers within the spin system, and both selected spins must have multiplicity 2.

## Construction

The function obtains pair operators from `state` in the caller's `spin_a,spin_b` order, then forms

`S=EE/4-(XX+YY+ZZ)`

where `EE`, `XX`, `YY`, and `ZZ` are the identity, x-, y-, and z-operator products on that pair. In Hilbert space the selected two-spin factor is the unit-trace singlet projector; the function applies no additional rescaling. Any other spins retain their identity factors. In Liouville space, `S` is the state vector representing the same operator, with coordinates in the configured basis rather than a wavefunction amplitude vector.

## Output

- `S` is a density matrix in Hilbert-space formalism or a state vector in Liouville-space formalism.

- Source: https://github.com/IlyaKuprov/Spinach/blob/main/kernel/states/singlet.m
- Wiki: https://spindynamics.org/wiki/index.php?title=singlet.m
- Related: [state](../state.md)
