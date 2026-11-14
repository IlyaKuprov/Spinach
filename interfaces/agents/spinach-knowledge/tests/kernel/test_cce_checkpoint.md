# tests/kernel/test_cce_checkpoint.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_checkpoint.m)

## Purpose

Partial ledger resume, identity rejection and write integrity.

## Behaviour

Checks the following named assertions in the source:

- `partial raw restart`.
- `partial reduced restart`.
- `completed bitmap`.
- `time identity`.
- `model identity`.
- `seed identity`.
- `cluster identity`.
- `raw corruption`.
- `completed task immutability`.
- `failed replacement integrity`.
- `write failure integrity`.
- `actual save failure integrity`.
- `execution-neutral restart`.

The test uses `new_test_result`, `test_true` and/or `test_close`; numerical tolerances and independent reference calculations are defined in the linked source.

## Inputs and outputs

```matlab
result=test_cce_checkpoint()
```

**Inputs**

- None; run with Spinach and its test libraries on the MATLAB path.

**Outputs**

- `result` — standard test-harness result with checks, pass/fail status and diagnostics. Printed numerical metrics are not additional result fields.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_checkpoint.m)
- [Test manifest](https://github.com/IlyaKuprov/Spinach/blob/main/tests/lib/test_manifest.m)
