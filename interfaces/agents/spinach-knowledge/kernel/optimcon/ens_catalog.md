# kernel/optimcon/ens_catalog.m

- Signature: `[catalog,ens_sizes]=ens_catalog(control)`

## Purpose

Ensemble case catalog for optimal control problems. Enumerates the Cartesian product of the state-target pairs, the drift generators, the control power levels, the resonance offsets, the phase cycle lines, and the distortion functions; then applies the ensemble correlation filters and the ensemble budget. Each row of the cata- log is one ensemble case to be simulated at each evaluation of the control sequence fidelit

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.

## Numerical / algorithmic content

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
