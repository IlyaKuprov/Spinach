# tests/kernel/test_cce_isotopologues.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_isotopologues.m)

## Purpose

H/D probabilities, actual Hilbert dimensions and local sampling.

## Behaviour

Checks the following named assertions in the source:

- `Unequal independent weights`.
- `H/D physical multiplicities`.
- `Complex configuration averaging occurs after each result`.
- `Weighted species-specific exact FIDs`.
- `Correlated joint isotope weights`.
- `Sampling agrees within a declared five-sigma bound`.
- `Vacancies remove spins entirely`.
- `Enumeration ceiling exposed`.

The test uses `new_test_result`, `test_true` and/or `test_close`; numerical tolerances and independent reference calculations are defined in the linked source.

## Inputs and outputs

```matlab
result=test_cce_isotopologues()
```

**Inputs**

- None; run with Spinach and its test libraries on the MATLAB path.

**Outputs**

- `result` — standard test-harness result with checks, pass/fail status and diagnostics. Printed numerical metrics are not additional result fields.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_isotopologues.m)
- [Test manifest](https://github.com/IlyaKuprov/Spinach/blob/main/tests/lib/test_manifest.m)
