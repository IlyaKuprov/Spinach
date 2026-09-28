# kernel/optimcon/ens_catalog.m

- Signature: `[catalog,ens_sizes]=ens_catalog(control)`

## Purpose

Builds an ensemble case catalog for optimal control. Enumerates the Cartesian product of state-target pairs, drift generators, control power levels, resonance-offset combinations, phase-cycle lines, and distortion functions, then applies ensemble correlation filters and the ensemble budget. Each catalog row identifies one case to simulate when evaluating control-sequence fidelity.

## Parameters / inputs

- `control` — control data structure produced by `optimcon.m`.

## Outputs

- `catalog` — `[n_cases x 6]` array of ensemble indices. Its columns index the state-target pair, drift generator, power level, offset combination, phase-cycle line, and distortion function, in that order.
- `ens_sizes` — `[1 x 6]` array of ensemble dimension sizes in the same column order, before correlation filters and budgeting.

## Numerical / algorithmic content

- The number of offset combinations is the product of the numbers of values in `control.offsets`, or 1 when `control.offsets` is empty.
- Correlation options in `control.ens_corrs` restrict cases: `rho_ens` assigns one state-target pair per remaining ensemble combination; `rho_drift` retains cases whose state-target and drift indices match; `power_drift` retains cases whose power and drift indices match.
- If `control.budget` is finite and at most 1, it is converted to a sample count by rounding that fraction of the filtered case count, with a minimum of 1. A budget smaller than the filtered case count selects a random subset using a fixed seed; the prior RNG state is restored afterward.

<https://spindynamics.org/wiki/index.php?title=ens_catalog.m>