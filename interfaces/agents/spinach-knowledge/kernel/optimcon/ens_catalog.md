# kernel/optimcon/ens_catalog.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/optimcon/ens_catalog.m`
- Signature: `[catalog,ens_sizes]=ens_catalog(control)`
- Total lines: 119

## Purpose

Ensemble case catalog for optimal control problems. Enumerates the Cartesian product of the state-target pairs, the drift generators, the control power levels, the resonance offsets, the phase cycle lines, and the distortion functions; then applies the ensemble correlation filters and the ensemble budget. Each row of the cata- log is one ensemble case to be simulated at each evaluation of the control sequence fidelit

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- control -control data structure produced by optimcon.m

## Outputs

- catalog -[n_cases x 6] array of ensemble indices; the col-
- umns index the state-target pair, the drift gene-
- rator, the power level, the offset combination,
- the phase cycle line, and the distortion function
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

- Called routines detected from the main body: `grumble()`, `cellfun()`, `ismember()`, `catalog()`, `rng()`, `randperm()`, `isstruct()`, `all()`, `isfield()`.
