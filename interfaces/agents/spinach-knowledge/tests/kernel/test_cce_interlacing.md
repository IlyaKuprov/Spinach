# tests/kernel/test_cce_interlacing.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_interlacing.m)

## Purpose

Check unique configurations, conditional weights, and shared denominator K.

## Behaviour

Checks the following named assertions in the source:

- `Overlapping cylinders yield six unique configurations, not eight occurrences`.
- `Projected probabilities are physical product weights`.
- `Projected raw average agrees with independent full-state enumeration`.
- `Duplicate-occurrence weighting gives a detectably different answer`.
- `Enumeration work exceeding the declared budget is refused`.
- `Included-supercluster family generates the same K automatically`.
- `Full-order W12 equals exact mixed bath`.

The test uses `new_test_result`, `test_true` and/or `test_close`; numerical tolerances and independent reference calculations are defined in the linked source.

## Inputs and outputs

```matlab
result=test_cce_interlacing()
```

**Inputs**

- None; run with Spinach and its test libraries on the MATLAB path.

**Outputs**

- `result` — standard test-harness result with checks, pass/fail status, and diagnostics. Printed numerical metrics are not additional result fields.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_interlacing.m)
- [Test manifest](https://github.com/IlyaKuprov/Spinach/blob/main/tests/lib/test_manifest.m)

Public optional-field contracts are explicit in reused fixtures: irrelevant method, event, geometry, and legacy numerical options are removed rather than passed silently.

Reference-spin_config tests preserve a non-default observed density and reject a bath density inconsistent with the reference basis assignment, before enumeration can overwrite it.

Native calls use clexpn with free_evolution or hahn_echo and the registered evaluator for direct conditional averages. Observed preparation belongs to experiment rho0, while sampled bath configurations own bath_rho and assignment. All original physical checks remain: conditioning unions and weights, memory refusal, reference consistency, worst residuals, shared denominators, full-bath bypass, Lindblad damping, driven responses, and convergence diagnostics. The former wrapper-schema check now checks numeric time-by-channel answers and separate diagnostic fields. The sampled-state memory check uses order zero so graph workspace does not pre-empt the intended configuration-count check.

Finite-control invocations use finite_drive with native physical descriptors, not a raw sequence field. Independent references retain their physical inputs and assertions; unused raw-event declarations are not part of the reference construction.
