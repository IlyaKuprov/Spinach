# tests/kernel/test_cce_interlacing.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_interlacing.m)

## Purpose

Check unique configurations, conditional weights and shared denominator K.

## Behaviour

Checks the following named assertions in the source:

- `Overlapping cylinders yield six unique configurations, not eight occurrences`.
- `Conditional probabilities are physical product weights`.
- `Conditional raw average agrees with independent full-state enumeration`.
- `Duplicate-occurrence weighting gives a detectably different answer`.
- `Enumeration work exceeding the declared budget is refused`.
- `Included-supercluster family generates the same K automatically`.
- `Full-order W12 equals exact mixed bath`.

The test uses `new_test_result`, `test_true` and/or `test_close`; numerical tolerances and independent reference calculations are defined in the linked source.

## Inputs and outputs

```matlab
result=test_cce_interlacing()
```

**Inputs**

- None; run with Spinach and its test libraries on the MATLAB path.

**Outputs**

- `result` — standard test-harness result with checks, pass/fail status and diagnostics. Printed numerical metrics are not additional result fields.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_interlacing.m)
- [Test manifest](https://github.com/IlyaKuprov/Spinach/blob/main/tests/lib/test_manifest.m)
