# tests/kernel/test_mecce_rates.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_mecce_rates.m)

## Purpose

Amplitude damping, dephasing, column convention and T1 limitation.

## Behaviour

Checks the following named assertions in the source:

- `C7 numerical check`.
- `C7 contract check`.
- `Zero-rate limit`.
- `Evolving exterior refusal`.

The test uses `new_test_result`, `test_true` and/or `test_close`; numerical tolerances and independent reference calculations are defined in the linked source.

## Inputs and outputs

```matlab
result=test_mecce_rates()
```

**Inputs**

- None; run with Spinach and its test libraries on the MATLAB path.

**Outputs**

- `result` — standard test-harness result with checks, pass/fail status and diagnostics. Printed numerical metrics are not additional result fields.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_mecce_rates.m)
- [Test manifest](https://github.com/IlyaKuprov/Spinach/blob/main/tests/lib/test_manifest.m)
