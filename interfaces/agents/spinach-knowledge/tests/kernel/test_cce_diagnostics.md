# tests/kernel/test_cce_diagnostics.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_diagnostics.m)

## Purpose

Coverage, error budgets, unphysical density, and raw residuals.

## Behaviour

Checks the following named assertions in the source:

- `analytic reduced density`.
- `physical density accepted`.
- `negative density rejected`.
- `density untouched`.
- `trace and Hermiticity`.
- `invalid time coverage`.
- `nonfinite bound`.
- `numerical budget`.
- `raw diagnostics`.

The test uses `new_test_result`, `test_true` and/or `test_close`; numerical tolerances and independent reference calculations are defined in the linked source.

## Inputs and outputs

```matlab
result=test_cce_diagnostics()
```

**Inputs**

- None; run with Spinach and its test libraries on the MATLAB path.

**Outputs**

- `result` — standard test-harness result with checks, pass/fail status, and diagnostics. Printed numerical metrics are not additional result fields.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_diagnostics.m)
- [Test manifest](https://github.com/IlyaKuprov/Spinach/blob/main/tests/lib/test_manifest.m)

## Native contract

The test calls `clexpn` for native numeric density, signal, and noise answers, and audits its raw and reduction records with `cce_diagnostics`. Standard and partitioned projected/generalised Hilbert/Liouville calculations retain the same diagnostic schema. Negative populations, trace/Hermiticity violations, incomplete time coverage, nonfinite bounds, exceeded numerical budgets, and sampled order-prefix means are checked without modifying the audited matrices. Native density rows are column-vectorised observed-spin density matrices.
