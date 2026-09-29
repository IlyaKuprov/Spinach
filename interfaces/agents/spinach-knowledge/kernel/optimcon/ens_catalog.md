# kernel/optimcon/ens_catalog.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/optimcon/ens_catalog.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=ens_catalog.m)

- Signature: `[catalog,ens_sizes]=ens_catalog(control)`

## Purpose and layout

Builds the case list used by ensemble optimal control. The Cartesian grid spans state-target pairs, drift generators, power levels, combinations of resonance offsets, phase-cycle rows, and distortion rows. `catalog` is an `n_cases x 6` numeric array; its columns index those six dimensions in that order. `ens_sizes` is a `1 x 6` vector of their sizes in the same order. The offset-combination count is the product of the number of values in each offset channel, or 1 when there are no offset channels.

`control` must be the optimal-control structure produced by `optimcon.m`. The function checks for the catalog inputs expected there, including `rho_init`, `rho_targ`, `ndrifts`, `pwr_levels`, `offsets`, `phase_cycle`, `distortion`, `ens_corrs`, and `budget`.

## Correlations and budget

The source applies these named correlation filters when present in `control.ens_corrs`:

- `rho_ens` removes the original state-target-pair column, deduplicates the remaining tuples, then rebuilds that column with a distinct sequential index.
- `rho_drift` retains rows whose state-target-pair and drift indices match.
- `power_drift` retains rows whose power-level and drift indices match.

Budgeting is applied after filtering. A finite budget at most 1 is treated as a fraction of the filtered case count, rounded to a sample count and clamped to at least 1. Otherwise the budget is used as a case count. If the resulting budget is smaller than the catalog, a subset is selected with a fixed Twister seed; the caller's prior MATLAB RNG state is restored afterward.
