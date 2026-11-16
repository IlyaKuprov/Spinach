# kernel/cce/cce_diagnostics.m

## Native role and invocation

Internal validity and uncertainty diagnostics for completed cluster reduction. The native context calls this service on its result record.

Public simulations enter through `clexpn(spin_system,pulse_sequence,parameters,assumptions)` with a registered experiment.

Internal/expert signature: `diagnostics=cce_diagnostics(model,raw,result,parameters)`. Engine records are not public experiment parameter structures.

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_diagnostics.m)

## Purpose

Report validity, conditioning, sampling and observed-density diagnostics.

## Behaviour


- Constructs per-observable invalid reasons, valid-time masks/fraction and numerical budget checks from result values, masks and error bounds.
- Reports cluster cost/order information and sample statistics when available; density outputs are checked for trace, Hermiticity and minimum eigenvalue.
- Never clips observables or projects densities to hide failed validity or physicality checks.

- Density-specific diagnostic fields are always present; they are empty when no observed density was requested, and `physicality_checked=false` distinguishes an unavailable check from a passed one.

## Inputs and outputs

```matlab
diagnostics=cce_diagnostics(model,raw,result,parameters)
```

**Inputs**

- `model` — prepared physical model; `raw` — cell of cluster-record arrays.
- `result` — CCE value, valid mask, abs_error_bound, and optional density/sampling_standard_error.
- `parameters` — numerical tolerances and observable declarations.

**Outputs**

- `diagnostics` — numerical/conditioning, cluster, sample and density diagnostic fields.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_diagnostics.m)
- [cce_reduce.m](cce_reduce.md)
