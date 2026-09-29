# kernel/utilities/kill_spin.m

## Purpose

Removes specified particles (spins or bosonic modes) from the `spin_system` structure and updates all dependent data. Source: [kernel/utilities/kill_spin.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/kill_spin.m).

## Behaviour

- Validates the input via an internal `grumble` consistency check: a logical mask must have exactly `spin_system.comp.nspins` elements; a numeric list must contain positive integers not exceeding `spin_system.comp.nspins`.
- Converts a logical `hit_list` to indices with `find`.
- Reports the number of particles being removed via `report`.
- Deletes the corresponding entries from `spin_system.comp.isotopes`, `spin_system.comp.types`, `spin_system.comp.labels`, `spin_system.comp.mults`, `spin_system.inter.gammas`, `spin_system.inter.basefrqs`, `spin_system.inter.zeeman.matrix`, `spin_system.inter.zeeman.ddscal`, and `spin_system.inter.giant.coeff`.
- Recomputes `spin_system.comp.iso_hash` with `md5_hash` if that field exists.
- Decrements `spin_system.comp.nspins` by `numel(hit_list)`.
- Deletes both rows and columns indexed by `hit_list` from `spin_system.inter.coupling.matrix`, `spin_system.inter.proxmatrix`, and removes the matching rows from `spin_system.inter.coordinates`.
- If `spin_system.inter.modes` exists:
  - Removes particle-indexed entries from the scalar mode fields `frqs`, `carriers`, `anharms`, `damp`, `dephase`.
  - Removes rows and columns from the mode pair fields `exchange`, `kerr`, `longitudinal`, `dispersive`, `coupling_mod`, `zeeman_mod`.
  - For `coupling_mod` and `zeeman_mod`, reindexes spin leaves inside retained modulation derivative orders by deleting `hit_list` columns (and also `hit_list` rows for `coupling_mod`).
  - If no particles of type `C`, `V`, or `T` remain, removes the entire `modes` field; otherwise removes `spin_system.inter.modes.strength` if present.
- Removes `hit_list` entries from non-empty relaxation arrays: `spin_system.rlx.r1_rates`, `r2_rates`, `lind_r1_rates`, `lind_r2_rates`; removes rows and columns from `srfk_mdepth`, `weiz_r1d`, `weiz_r2d` when non-empty.
- Reindexes scalar relaxation source spins (`spin_system.rlx.srsk_sources`) by rebuilding the source mask over the pre-removal spin count and deleting `hit_list` positions.
- Reindexes each subsystem in `spin_system.chem.parts` the same way, and removes rows and columns from `spin_system.chem.flux_rate` when non-empty.
- Reindexes `spin_system.chem.rp_electrons`; if `spin_system.chem.rp_rates` is non-empty and fewer than two electrons remain, raises the error `cannot destroy an essential electron in a radical pair system system.`
- Removes `spin_system.bas` if present, with a warning that basis set information must be re-created.
- Removes `spin_system.inter.conmatrix` if present.
- Removes `spin_system.comp` fields `sym_group`, `sym_spins`, `sym_a1g_only` if `sym_group` is present.
- Removes `spin_system.inter.assumptions` if present, and removes `strength` subfields from `spin_system.inter.zeeman`, `spin_system.inter.giant`, and `spin_system.inter.coupling` when present, each with a warning that assumption information must be re-created.

## Inputs and outputs

**Syntax**

```matlab
spin_system = kill_spin(spin_system, hit_list)
```

**Inputs**

- `spin_system` — primary Spinach data structure.
- `hit_list` — vector of integers or logical vector giving particle numbers in the unified isotope list to be removed.

**Outputs**

- `spin_system` — the data structure with the indicated particles and dependent information (basis, assumptions) removed.

## References

- Spinach Wiki page: [kill_spin.m](https://spindynamics.org/wiki/index.php?title=kill_spin.m)
- Source file: [kernel/utilities/kill_spin.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/kill_spin.m)
