# tests/kernel/test_cce_projection.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_projection.m)

## Purpose

Check explicit states, conditional fields and resolved virtual levels.

## Behaviour

Checks the following named assertions in the source:

- `alpha expectation`.
- `beta expectation`.
- `second order off by default`.
- `virtual coefficient`.
- `second-order spectrum`.
- `near-degenerate denominator rejected`.
- `overlap guard`.
- `gCCE correction rejected`.
- `multiple central register`.

The test uses `new_test_result`, `test_true` and/or `test_close`; numerical tolerances and independent reference calculations are defined in the linked source.

## Inputs and outputs

```matlab
result=test_cce_projection()
```

**Inputs**

- None; run with Spinach and its test libraries on the MATLAB path.

**Outputs**

- `result` — standard test-harness result with checks, pass/fail status and diagnostics. Printed numerical metrics are not additional result fields.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_projection.m)
- [Test manifest](https://github.com/IlyaKuprov/Spinach/blob/main/tests/lib/test_manifest.m)
