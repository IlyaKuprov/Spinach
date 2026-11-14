# tests/kernel/test_gcce_central_mediation.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_gcce_central_mediation.m)

## Purpose

Two bath spins with zero direct coupling share a flip-flopping centre.

## Behaviour

Checks the following named assertions in the source:

- `pair retained`.
- `exact pair`.
- `nonfactorization measured`.

The test uses `new_test_result`, `test_true` and/or `test_close`; numerical tolerances and independent reference calculations are defined in the linked source.

## Inputs and outputs

```matlab
result=test_gcce_central_mediation()
```

**Inputs**

- None; run with Spinach and its test libraries on the MATLAB path.

**Outputs**

- `result` — standard test-harness result with checks, pass/fail status and diagnostics. Printed numerical metrics are not additional result fields.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_gcce_central_mediation.m)
- [Test manifest](https://github.com/IlyaKuprov/Spinach/blob/main/tests/lib/test_manifest.m)
