# tests/kernel/test_cce_conditioning.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_conditioning.m)

## Purpose

Expose unresolved ratios rather than regularise or lose phase.

## Behaviour

Checks the following named assertions in the source:

- `Exact conditioning masks and reason codes`.
- `Explicit fail-fast zero policy`.
- `Signed phases across logarithm cuts`.
- `Caller error-budget refusal`.
- `Resolved logarithmic dynamic range`.
- `Overflow reason`.
- `Uncertain nonzero baseline`.
- `Exact signed-power cancellation`.
- `Direct singularity is reported honestly`.

The test uses `new_test_result`, `test_true` and/or `test_close`; numerical tolerances and independent reference calculations are defined in the linked source.

## Inputs and outputs

```matlab
result=test_cce_conditioning()
```

**Inputs**

- None; run with Spinach and its test libraries on the MATLAB path.

**Outputs**

- `result` — standard test-harness result with checks, pass/fail status and diagnostics. Printed numerical metrics are not additional result fields.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_conditioning.m)
- [Test manifest](https://github.com/IlyaKuprov/Spinach/blob/main/tests/lib/test_manifest.m)
