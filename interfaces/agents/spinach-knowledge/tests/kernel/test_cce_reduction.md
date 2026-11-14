# tests/kernel/test_cce_reduction.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_reduction.m)

## Purpose

Compare direct factors, signed powers and additive inclusion-exclusion.

## Behaviour

Checks the following named assertions in the source:

- `Resolved direct and powers masks`.
- `Signed powers versus known factors`.
- `Direct recursion versus signed powers`.
- `Unit initial value`.
- `Permutation-gather bitwise identity`.
- `All order prefixes`.
- `Additive polynomial full order`.
- `Additive quadratic saturation`.
- `Central-only exact zero`.

The test uses `new_test_result`, `test_true` and/or `test_close`; numerical tolerances and independent reference calculations are defined in the linked source.

## Inputs and outputs

```matlab
result=test_cce_reduction()
```

**Inputs**

- None; run with Spinach and its test libraries on the MATLAB path.

**Outputs**

- `result` — standard test-harness result with checks, pass/fail status and diagnostics. Printed numerical metrics are not additional result fields.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_reduction.m)
- [Test manifest](https://github.com/IlyaKuprov/Spinach/blob/main/tests/lib/test_manifest.m)
