# kernel/cce/cce_analyse.m

## Native role and invocation

Explicit decay-envelope analysis utility. This is not invoked automatically by clexpn and does not identify T1 from a generic signal.

Public simulations enter through `clexpn(spin_system,pulse_sequence,parameters,assumptions)` with a registered experiment.

Internal/expert signature: `analysis=cce_analyse(result,parameters)`. Engine records are not public experiment parameter structures.

[MATLAB source](../../../../../kernel/cce/cce_analyse.m)

## Purpose

Fit an explicitly declared decay envelope without an automatic T1 claim.

## Behaviour

- Preserves input values/masks and returns status not_requested when fitting is not enabled.
- Fits A*exp(-(t/T)^beta) in a user-declared interval; exponential fixes beta=1, stretched_exponential fits it. No offset is inferred.
- Uses the selected real or magnitude envelope and explicit inverse-variance weights; optional confidence is `none` or `linearized` at the requested level. Invalid data are diagnosed rather than silently fitted.

## Inputs and outputs

**Inputs**

- `result` — times, value, and valid.
- `parameters.analysis` — fit (default false); when enabled: observable column, envelope, model, interval in seconds, positive weights, confidence, and confidence_level.

**Outputs**

- `analysis` — retained data/masks, status, warnings, and fit parameters/diagnostics.

## References

- [cce_diagnostics.m](cce_diagnostics.md)

## Native context inputs

This standalone post-processing service validates only its retained-data and `parameters.analysis` contracts. Native `observed`, `times`, and nested `cce` algorithm options are permitted but do not alter fitting. Unknown analysis fields and nonscalar options are rejected locally.
