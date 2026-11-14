# tests/kernel/test_gcce_exact.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_gcce_exact.m)

## Purpose

Four bath spins, complex pure/mixed states and physical Hahn pulses.

## Behaviour

Checks the following named assertions in the source:

- `raw full-register exact`.

The test uses `new_test_result`, `test_true` and/or `test_close`; numerical tolerances and independent reference calculations are defined in the linked source.

## Inputs and outputs

```matlab
result=test_gcce_exact()
```

**Inputs**

- None; run with Spinach and its test libraries on the MATLAB path.

**Outputs**

- `result` — standard test-harness result with checks, pass/fail status and diagnostics. Printed numerical metrics are not additional result fields.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_gcce_exact.m)
- [Test manifest](https://github.com/IlyaKuprov/Spinach/blob/main/tests/lib/test_manifest.m)

Public optional-field contracts are explicit in reused fixtures: irrelevant method, event, geometry, and legacy numerical options are removed rather than passed silently.
