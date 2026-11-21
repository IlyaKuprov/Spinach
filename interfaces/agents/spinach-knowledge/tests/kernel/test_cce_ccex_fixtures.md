# tests/kernel/test_cce_ccex_fixtures.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_ccex_fixtures.m)

## Purpose

Analytic FID/Hahn and basic01 nearest-six CCE2.

## Behaviour

Checks the following named assertions in the source:

- `CCEX raw subset`.
- `CCEX printed curve`.

The test uses `new_test_result`, `test_true` and/or `test_close`; numerical tolerances and independent reference calculations are defined in the linked source.

## Inputs and outputs

```matlab
result=test_cce_ccex_fixtures()
```

**Inputs**

- None; run with Spinach and its test libraries on the MATLAB path.

**Outputs**

- `result` — standard test-harness result with checks, pass/fail status, and diagnostics. Printed numerical metrics are not additional result fields.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_ccex_fixtures.m)

## Native interface

Uses `clexpn` with `@free_evolution` or `@hahn_echo` and projected expansion policy. Physical support is looked up in `diagnostics.clusters`; raw responses and numerical bounds remain in `diagnostics.raw` and `diagnostics.reduction`. Every pinned CCEX complex-trace comparison is retained.

Pinned external JSON traces remain unchanged. Receiver conversion is explicit: CCEX Hahn spin-manifold overlaps are conjugated to the fixed physical receiver; FID references are unchanged. Raw-subset and final complex comparisons use the same convention.
