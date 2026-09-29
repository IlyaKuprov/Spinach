# kernel/utilities/get_coupling.m

## Purpose

Extracts the 3x3 coupling tensor between a pair of spins from the `spin_system` data structure.

Source: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/get_coupling.m>

## Behaviour

- Syntax: `A=get_coupling(spin_system,n,k)`.
- Runs a consistency check (`grumble`) that errors if `spin_system` lacks the `inter` or `inter.coupling` fields, or if `n` or `k` is not a positive real integer.
- Retrieves the forward coupling `spin_system.inter.coupling.matrix{n,k}` and the backward coupling `spin_system.inter.coupling.matrix{k,n}`.
- Replaces empty entries with 3x3 zero matrices.
- Returns the sum of the forward and backward coupling tensors.

## Inputs and outputs

Inputs:

- `spin_system` — spin system data structure containing coupling information.
- `n`, `k` — indices of the two spins as they appear in `spin_system.comp.isotopes`; must be positive real integers.

Outputs:

- `A` — 3x3 coupling tensor in rad/s.

## References

- Spinach Wiki: <https://spindynamics.org/wiki/index.php?title=get_coupling.m>
