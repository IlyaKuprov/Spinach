# kernel/cce/cce_diagnostics.m

## Native role and invocation

Internal validity and uncertainty diagnostics for completed cluster reduction. Call this expert service on the native context’s raw and reduction records.

Public simulations enter through `clexpn(spin_system,pulse_sequence,parameters,assumptions)` with a registered experiment.

Internal/expert signature: `diagnostics=cce_diagnostics(spin_system,raw,result,parameters)`. Engine records are not public experiment parameter structures.

[MATLAB source](../../../../../kernel/cce/cce_diagnostics.m)

## Purpose

Report validity, conditioning, sampling, and observed-density diagnostics.

## Behaviour

- Constructs per-observable invalid reasons, valid-time masks/fraction, and numerical budget checks from result values, masks, and error bounds.
- Reports raw-record counts, invalid values, maximum spectral and orthogonality residuals, and supplied sampling standard errors; density outputs are checked for trace, Hermiticity, and minimum eigenvalue.
- Never clips observables or projects densities to hide failed validity or physicality checks.
- Density-specific diagnostic fields are always present; they are empty when no observed density was requested, and `physicality_checked=false` distinguishes an unavailable check from a passed one.

## Inputs and outputs

**Inputs**

- `spin_system` — compiled spin system; `raw` — cell of cluster-record arrays.
- `result` — CCE value, valid mask, abs_error_bound, and optional density/sampling_standard_error.
- `parameters` — numerical tolerances and observable declarations.

**Outputs**

- `diagnostics` — numerical validity and budgets, raw-record counts, spectral residuals, sampling uncertainty, and density diagnostic fields.

## References

- [cce_reduce.m](cce_reduce.md)

Requests pass the resolved `parameters`; tolerances are read from `spin_system.tols`. When `parameters.output` is `density`, the reduction’s time-by-vectorised-density values are inspected as square observed-spin density matrices using the compiled model dimensions. The routine never repairs these arrays; all coverage, error-budget, trace, Hermiticity, and minimum-eigenvalue decisions retain their original arithmetic.

Validation is local to this service and its consumed numerical records. Public experiment and approximation options are validated by `cce_options` before native context dispatch.

A supplied unavailable (`NaN`) spectral or orthogonality residual remains unavailable in the maximum across raw records; it is not discarded in favour of a finite residual from another record. Sampling standard errors, including `NaN`, are copied unchanged.

Native result shapes and numerical policy are validated before interpreting vectorised density values; the observed Hilbert dimension must match the supplied density columns.
