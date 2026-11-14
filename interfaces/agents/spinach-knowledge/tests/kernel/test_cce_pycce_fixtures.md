# tests/kernel/test_cce_pycce_fixtures.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_pycce_fixtures.m)

## Purpose

Three fixtures, raw subsets and complex final reductions.

## Behaviour

Checks the following named assertions in the source:

- `PyCCE raw subset`.
- `PyCCE final coherence`.
- `PyCCE reduction only`.

The test uses `new_test_result`, `test_true` and/or `test_close`; numerical tolerances and independent reference calculations are defined in the linked source.

## Inputs and outputs

```matlab
result=test_cce_pycce_fixtures()
```

**Inputs**

- None; run with Spinach and its test libraries on the MATLAB path.

**Outputs**

- `result` — standard test-harness result with checks, pass/fail status and diagnostics. Printed numerical metrics are not additional result fields.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_pycce_fixtures.m)
- [Test manifest](https://github.com/IlyaKuprov/Spinach/blob/main/tests/lib/test_manifest.m)
