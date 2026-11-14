# tests/kernel/test_cce_parallel.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_parallel.m)

## Purpose

Serial, one/two workers, common samples and resource caps.

## Behaviour

Checks the following named assertions in the source:

- `one-worker raw identity`.
- `one-worker reduction`.
- `two-worker sample identity`.
- `two-worker raw identity`.
- `two-worker reduction`.
- `state scheduling`.
- `analytic sampled coherence`.
- `task cap`.
- `memory cap`.

The test uses `new_test_result`, `test_true` and/or `test_close`; numerical tolerances and independent reference calculations are defined in the linked source.

## Inputs and outputs

```matlab
result=test_cce_parallel()
```

**Inputs**

- None; run with Spinach and its test libraries on the MATLAB path.

**Outputs**

- `result` — standard test-harness result with checks, pass/fail status and diagnostics. Printed numerical metrics are not additional result fields.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_parallel.m)
- [Test manifest](https://github.com/IlyaKuprov/Spinach/blob/main/tests/lib/test_manifest.m)
