# tests/kernel/test_cce_spectrum.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_spectrum.m)

## Purpose

Parseval, inverse normalisation, Lorentzian, sidedness, and units.

## Behaviour

Checks the following named assertions in the source:

- `C7 numerical check`.
- `C7 contract check`.

The test uses `new_test_result`, `test_true` and/or `test_close`; numerical tolerances and independent reference calculations are defined in the linked source.

## Inputs and outputs

```matlab
result=test_cce_spectrum()
```

**Inputs**

- None; run with Spinach and its test libraries on the MATLAB path.

**Outputs**

- `result` — standard test-harness result with checks, pass/fail status, and diagnostics. Printed numerical metrics are not additional result fields.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_spectrum.m)

## Native contract

The standalone native Fourier service uses `cce_spectrum(times,correlation,parameters)`. No legacy CCE driver, compiler, or method dispatch is used. All direct-transform, normalisation, Lorentzian, filter-unit, stationarity, one-sided ordering, and per-column scale assertions remain in this test.
