# tests/kernel/test_cce_analysis.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_analysis.m)

## Purpose

Synthetic windows, fit assumptions, and invalid/oscillatory refusal.

## Behaviour

Checks the following named assertions in the source:

- `C7 contract check`.
- `C7 numerical check`.

The test uses `new_test_result`, `test_true` and/or `test_close`; numerical tolerances and independent reference calculations are defined in the linked source.

## Inputs and outputs

```matlab
result=test_cce_analysis()
```

**Inputs**

- None; run with Spinach and its test libraries on the MATLAB path.

**Outputs**

- `result` — standard test-harness result with checks, pass/fail status, and diagnostics. Printed numerical metrics are not additional result fields.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_analysis.m)

Native parameter structures preserve fitted amplitude, time constant and exponent bitwise. A misspelled analysis option is rejected by the local analysis contract.

## Native contract

The standalone native analysis service accepts a numeric-record structure containing `times`, `value`, and `valid`, plus `parameters.analysis`. No legacy CCE driver, compiler, or method dispatch is used. The test retains all fitting, interval, confidence, invalid-window, nondecay, immutability, and unknown-field assertions.
