# tests/kernel/test_dynamic_chem_geometry_suite.m

- Signature: `result=test_dynamic_chem_geometry_suite()`

## Purpose

Checks deterministic chemistry and geometry helpers against small fixtures with known coordinates, tensor values, and metadata.

## Numerical / algorithmic content

- `cubic_lattice('13C',2,2)` produces eight isotope entries; the second coordinate is `[0 0 2]`, and the periodic vectors are `[4 0 0]`, `[0 4 0]`, and `[0 0 4]`.
- `dihedral()` returns `-90` degrees for the selected four-point geometry in Spinach convention.
- `xyz2pd()` bins two in-range points into the first two x bins of a `2x2x2` grid and discards the point outside the `[0 1]` coordinate ranges.
- In the three-spin fixture at `[0 0 0]`, `[2 0 0]`, and `[0.5 0 0]` Angstrom, `nearest_spin()` selects spin 3 at distance `0.5` Angstrom. `which_subst()` assigns spins 1 and 2 to part 1 and spin 3 to part 2.
- `get_coupling()` adds the forward and reverse cells to give `diag([5 7 9])`.
- `chemshifts()` returns `[1 2 0]` ppm and `[-100 -50 0]` Hz; `offsetof(spin_system,1)` returns `-100` Hz.
- With `ddscal=2*eye(3)`, `gtensorof()` is checked against `-ddscal*gamma*hbar/muB`.
- `shift_iso()` replaces the isotropic part of `diag([1 2 6])` by 10, giving `diag([8 9 13])`, and leaves an unselected identity tensor unchanged.

## Outputs

`result` records the regression assertions for these helper behaviors.
