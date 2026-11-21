# tests/kernel/test_cce_pycce_fixtures.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_pycce_fixtures.m)

## Purpose

Three fixtures, raw subsets, and complex final reductions.

## Behaviour

Checks the following named assertions in the source:

- `PyCCE raw subset`.
- `PyCCE final signal`.
- `PyCCE reduction only`.

The test uses `new_test_result`, `test_true` and/or `test_close`; numerical tolerances and independent reference calculations are defined in the linked source.

## Inputs and outputs

```matlab
result=test_cce_pycce_fixtures()
```

**Inputs**

- None; run with Spinach and its test libraries on the MATLAB path.

**Outputs**

- `result` — standard test-harness result with checks, pass/fail status, and diagnostics. Printed numerical metrics are not additional result fields.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_pycce_fixtures.m)

The native projected and generalised experiments both request normalised detected amplitude explicitly. All raw-trace and final-complex-signal comparisons remain unchanged.

## Native interface

Uses `clexpn` with `@hahn_echo`, explicit transition support and `normalized=true` for both projected and generalised physical signals. Approximation is selected with `parameters.cce.approximation`. Separate raw/subset/final and external-only reduction assertions remain unchanged.

Pinned external JSON traces remain unchanged. Receiver conversion is explicit: projected PyCCE Hahn traces already match the fixed physical receiver, whereas generalised pulse-rotated detection retains its conjugation. Raw-subset and final complex comparisons use the same convention.
