# kernel/utilities/nearest_spin.m

## Purpose

Returns the index of the spin nearest to a specified spin in a spin system. Only spins for which Cartesian coordinates are available are considered.

Source: [kernel/utilities/nearest_spin.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/nearest_spin.m)

## Behaviour

- Syntax: `[k,d]=nearest_spin(spin_system,n)`.
- The function first validates its inputs through an internal consistency check (`grumble`):
  - `n` must be a positive real integer scalar; otherwise it errors with `'n must be a positive real integer'`.
  - `n` must not exceed the number of isotopes in `spin_system.comp.isotopes`; otherwise it errors with `'the specified spin does not exist'`.
  - The specified spin must have coordinates assigned in `spin_system.inter.coordinates{n}`; otherwise it errors with `'the specified spin does not have coordinates'`.
- The search initialises `k` as empty and `d` as `inf`, then iterates over all spins `s` in `spin_system.inter.coordinates`.
- A candidate spin `s` is considered only if `s` differs from `n` and its coordinate cell is non-empty.
- For each candidate, the Euclidean distance (`norm(...,2)`) between its coordinates and those of spin `n` is computed; if it is smaller than the current best distance, `k` and `d` are updated.
- If no other spin has coordinates, the function errors with `'no other spin has coordinates'`.

## Inputs and outputs

Inputs:

- `spin_system` — spin system object containing inter-spin coordinates in `spin_system.inter.coordinates` (cell array) and isotope list in `spin_system.comp.isotopes`.
- `n` — index of the spin in question; positive real integer scalar.

Outputs:

- `k` — index of the nearest spin.
- `d` — distance to the nearest spin, in Angstrom.

## References

- Spinach Wiki: [nearest_spin.m](https://spindynamics.org/wiki/index.php?title=nearest_spin.m)
- Source file: [kernel/utilities/nearest_spin.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/nearest_spin.m)
