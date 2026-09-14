# kernel/optimcon/ens_catalog.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/optimcon/ens_catalog.m`
- Signature: `[catalog,ens_sizes]=ens_catalog(control)`
- Total lines: 124

## Purpose

Ensemble case catalog for optimal control problems. Enumerates the Cartesian product of the state-target pairs, the drift generators, the control power levels, the resonance offsets, the phase cycle lines, and the distortion functions; then applies the ensemble correlation filters and the ensemble budget. Each row of the cata- log is one ensemble case to be simulated at each evaluation of the control sequence fidelit

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 33-34: Check consistency; implemented by `grumble(control)`.
- Lines 36-37: Get offset ensemble size; implemented by `off_ens_sizes=cellfun(@numel,control.offsets)`.
- Lines 44-45: Extract ensemble grid dimensions; implemented by `n_state_pairs=numel(control.rho_init)`.
- Lines 51-53: Record the grid dimensions; implemented by `ens_sizes=[n_state_pairs n_ens_systems n_power_levls n_offset_vals n_phase_specs n_distortions]`.
- Lines 55-56: Create a catalog of the ensemble; implemented by `catalog=(1:n_state_pairs)'`.
- Lines 63-64: Order the cases by drift generator so that a contiguous block of cases spans few drifts; implemented by `catalog=sortrows(catalog,2)`.
- Lines 66-67: Ensemble correlation: own state pair for each member; implemented by `if ismember('rho_ens',control.ens_corrs)`.
- Lines 73-74: Ensemble correlation: own state pair for each drift; implemented by `if ismember('rho_drift',control.ens_corrs)`.
- Lines 78-79: Ensemble correlation: own control power for each drift; implemented by `if ismember('power_drift',control.ens_corrs)`.
- Lines 83-84: Count the full ensemble size; implemented by `n_cases=size(catalog,1)`.
- Lines 86-87: Convert fractional budget into sample count; implemented by `ens_budget=control.budget`.
- Lines 93-94: Apply ensemble budget; implemented by `if ens_budget<n_cases`.
- Lines 96-97: Get RNG into a reproducible state; implemented by `rng_state=rng; rng(5318008,'twister')`.
- Lines 99-100: Draw a random subset of the ensemble, drift order restored; implemented by `catalog=sortrows(catalog(randperm(n_cases,ens_budget),:),2)`.
- Lines 102-103: Release RNG; implemented by `rng(rng_state)`.

### Control flow inferred from the code

- Line 38: conditional branch on `~isempty(off_ens_sizes)`.
- Line 67: conditional branch on `ismember('rho_ens',control.ens_corrs)`.
- Line 74: conditional branch on `ismember('rho_drift',control.ens_corrs)`.
- Line 79: conditional branch on `ismember('power_drift',control.ens_corrs)`.
- Line 88: conditional branch on `isfinite(ens_budget)&&(ens_budget<=1)`.
- Line 94: conditional branch on `ens_budget<n_cases`.

### Key state/data transformations

- Lines 37: computes `off_ens_sizes` using `off_ens_sizes=cellfun(@numel,control.offsets)`.
- Lines 39: computes `n_offset_vals` using `n_offset_vals=prod(off_ens_sizes)`.
- Lines 45: computes `n_state_pairs` using `n_state_pairs=numel(control.rho_init)`.
- Lines 46: computes `n_ens_systems` using `n_ens_systems=control.ndrifts`.
- Lines 47: computes `n_power_levls` using `n_power_levls=numel(control.pwr_levels)`.
- Lines 48: computes `n_phase_specs` using `n_phase_specs=size(control.phase_cycle,1)`.
- Lines 49: computes `n_distortions` using `n_distortions=size(control.distortion,1)`.
- Lines 52-53: computes `ens_sizes` using `ens_sizes=[n_state_pairs n_ens_systems n_power_levls n_offset_vals n_phase_specs n_distortions]`.
- Lines 56: computes `catalog` using `catalog=(1:n_state_pairs)'`.
- Lines 84: computes `n_cases` using `n_cases=size(catalog,1)`.
- Lines 87: computes `ens_budget` using `ens_budget=control.budget`.
- Lines 97: computes `rng_state` using `rng_state=rng; rng(5318008,'twister')`.

### Local helper functions

- Line 110: `grumble()` — `function grumble(control)`. The purpose of abstraction is not to be vague, but to create a new semantic level in which one can be absolutely precise.
  - Representative operation: `if ~isstruct(control)`.
  - Representative operation: `error('control must be a data structure produced by optimcon.m')`.

## Parameters / inputs

- control -control data structure produced by optimcon.m

## Outputs

- catalog -[n_cases x 6] array of ensemble indices; the col-
- umns index the state-target pair, the drift gene-
- rator, the power level, the offset combination,
- the phase cycle line, and the distortion function;
- rows are ordered by the drift generator index so
- that contiguous blocks of cases span few drifts
- ens_sizes -[1 x 6] array of the ensemble dimension sizes the
- catalog was built from, in the same column order

## Implementation structure

- Ensemble case catalog for optimal control problems. Enumerates the
- Cartesian product of the state-target pairs, the drift generators,
- the control power levels, the resonance offsets, the phase cycle
- lines, and the distortion functions; then applies the ensemble
- correlation filters and the ensemble budget. Each row of the cata-
- log is one ensemble case to be simulated at each evaluation of the
- control sequence fidelity. Syntax:
- [catalog,ens_sizes]=ens_catalog(control)
- control -control data structure produced by optimcon.m
- catalog -[n_cases x 6] array of ensemble indices; the col-
- umns index the state-target pair, the drift gene-
- rator, the power level, the offset combination,

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `cellfun()`, `numel()`, `prod()`, `kron()`, `ones()`, `sortrows()`, `ismember()`, `unique()`, `isfinite()`, `round()`, `max()`, `rng()`, `randperm()`.
