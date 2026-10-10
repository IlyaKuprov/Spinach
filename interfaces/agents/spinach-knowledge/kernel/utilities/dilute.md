# kernel/utilities/dilute.m

## Purpose

Splits a spin system into several independent subsystems, each containing only one instance (or one tuple) of a user-specified isotope that is deemed "dilute". All spin system data is updated accordingly. An existing basis is rebuilt for every isotopomer through `kill_spin`, including the local descriptors, offsets, and cache hash. Spin-free substance blocks are retained.

Source: [kernel/utilities/dilute.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/dilute.m)

## Behaviour

- Syntax: `subsystems=dilute(spin_system,isotope,tuples)`.
- If `tuples` is not supplied, it defaults to `1` (singles).
- Input consistency is enforced by an internal `grumble` check: `isotope` must be a character string, and `tuples` must be a real, scalar, positive integer.
- The user is informed via `report` that the isotope is being treated as dilute and which tuple size is being picked.
- The spins belonging to the dilute species are located by matching `spin_system.comp.isotopes` against the isotope string with `strcmp` via `cellfun`; the number of instances found is reported.
- If the number of dilute spins found is less than `tuples`, the function errors with `not enough <isotope> in the system to make the selection.`.
- All combinations of the dilute spins taken `tuples` at a time are generated with `nchoosek`, then sorted row-wise and made unique.
- For each combination, a new spin system is created by calling `kill_spin` on the original system, removing the dilute spins not in that combination (`setdiff(dilute_spins,combos(n,:))`).
- The number of spin subsystems returned is reported.

## Inputs and outputs

Inputs:

- `spin_system` — spin system object obtained from the `create.m` function.
- `isotope` — isotope specification string, for example `'13C'`.
- `tuples` — `1` returns all subsystems with a single instance of the isotope, `2` returns all subsystems where two instances of the isotope are present, etc. Optional; defaults to `1`.

Outputs:

- `subsystems` — a cell array of spin system objects, with each cell corresponding to one of the newly formed isotopomers.

## References

- Spinach Wiki page for dilute.m: https://spindynamics.org/wiki/index.php?title=dilute.m
- Source file: https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/dilute.m
