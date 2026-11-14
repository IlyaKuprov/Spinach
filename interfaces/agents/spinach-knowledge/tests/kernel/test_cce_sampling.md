# tests/kernel/test_cce_sampling.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_sampling.m)

## Purpose

Compare weighted samples, exact mixtures and clock-transition exterior fields.

## Behaviour

Checks the following named assertions in the source:

- `Exhaustive weighted samples equal exact product mixture`.
- `Average of products differs from product of averages`.
- `Clock-transition static eigenstate exterior matches exact retained bath`.
- `Omitting sampled exterior field changes clock-transition dynamics`.
- `Central gap has zero first-order longitudinal sensitivity`.

The test uses `new_test_result`, `test_true` and/or `test_close`; numerical tolerances and independent reference calculations are defined in the linked source.

## Inputs and outputs

```matlab
result=test_cce_sampling()
```

**Inputs**

- None; run with Spinach and its test libraries on the MATLAB path.

**Outputs**

- `result` — standard test-harness result with checks, pass/fail status and diagnostics. Printed numerical metrics are not additional result fields.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_sampling.m)
- [Test manifest](https://github.com/IlyaKuprov/Spinach/blob/main/tests/lib/test_manifest.m)
