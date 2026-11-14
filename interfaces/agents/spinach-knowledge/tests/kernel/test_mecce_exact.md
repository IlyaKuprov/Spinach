# tests/kernel/test_mecce_exact.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_mecce_exact.m)

## Purpose

Zero jumps, dense Liouville full bath and correlated support.

## Behaviour

Checks the following named assertions in the source:

- `C7 numerical check`.
- `ME-gCCE Hahn zero jumps`.
- `ME-CCE Hahn zero jumps`.
- `Full-bath conditional ME reference`.
- `Dissipative Hahn reference`.
- `C7 contract check`.
- `PyCCE reconstructed ME-CCE2`.
- `PyCCE reconstructed ME-gCCE2`.
- `ME tutorial validity`.

The test uses `new_test_result`, `test_true` and/or `test_close`; numerical tolerances and independent reference calculations are defined in the linked source.

## Inputs and outputs

```matlab
result=test_mecce_exact()
```

**Inputs**

- None; run with Spinach and its test libraries on the MATLAB path.

**Outputs**

- `result` — standard test-harness result with checks, pass/fail status and diagnostics. Printed numerical metrics are not additional result fields.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_mecce_exact.m)
- [Test manifest](https://github.com/IlyaKuprov/Spinach/blob/main/tests/lib/test_manifest.m)
