# tests/kernel/test_gcce_density.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_gcce_density.m)

## Purpose

Baseline zeros mask truncated ratios and full order bypasses them.

## Native interface

Uses native `clexpn` density arrays (one column-major vectorised density per time row) and separate reduction diagnostics. Trace, Hermiticity, and minimum eigenvalue are computed directly without repairing the densities. A direct reducer probe disables the full-order bypass to retain the undefined-ratio masking assertion. A radius-selected singleton reproduces its raw density, including the exact zeros at zero duration. Reversing the excluded spin reverses the transverse coherence via its static mean field; admitting that site at the radius boundary restores the genuinely truncated mask at order one. Sampling checks inspect native physical configurations, including deterministic probabilities and unchanged global RNG.

## Behaviour

Checks the following named assertions in the source:

- `full-order bypass finite`.
- `initial density`.
- `unit trace`.
- `Hermitian`.
- `positive`.
- `observed baseline preserved`.
- `elementwise baseline zeros masked`.
- `raw density retained`.
- `negative input rejected`.
- `zero normalization rejected`.
- `global RNG unchanged`.
- `fixed probability state`.
- `reproducible assignments`.
- `static exterior exact`.
- `full-order sampling estimators`.
- `sample errors finite`.
- `explicit pulse timestamps`.
- `spin-one selective distinction`.

The test uses `new_test_result`, `test_true` and/or `test_close`; numerical tolerances and independent reference calculations are defined in the linked source.

## Inputs and outputs

```matlab
result=test_gcce_density()
```

**Inputs**

- None; run with Spinach and its test libraries on the MATLAB path.

**Outputs**

- `result` — standard test-harness result with checks, pass/fail status, and diagnostics. Printed numerical metrics are not additional result fields.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_gcce_density.m)
