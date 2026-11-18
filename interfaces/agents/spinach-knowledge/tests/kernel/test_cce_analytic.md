# tests/kernel/test_cce_analytic.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_analytic.m)

## Purpose

Zero bath, isolated phase and one-spin Hahn modulation.

## Behaviour

Checks the following named assertions in the source:

- `complex FID`.
- `one-spin Hahn formula`.
- `zero time is unity`.
- `empty bath phase`.
- `uncoupled bath`.

The test uses `new_test_result`, `test_true` and/or `test_close`; numerical tolerances and independent reference calculations are defined in the linked source.

## Inputs and outputs

```matlab
result=test_cce_analytic()
```

**Inputs**

- None; run with Spinach and its test libraries on the MATLAB path.

**Outputs**

- `result` — standard test-harness result with checks, pass/fail status and diagnostics. Printed numerical metrics are not additional result fields.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_analytic.m)
- [Test manifest](https://github.com/IlyaKuprov/Spinach/blob/main/tests/lib/test_manifest.m)

Native calls use free_evolution and hahn_echo with the original complex phase, total-time mesh, analytic one-spin modulation, isolated observed-spin phase, and disconnected-bath assertions. The absolute pulse is a refocusing descriptor at 0.001 seconds; native raw pulse_applied metadata verifies exclusion before the timestamp and inclusion at equality. The independent manifold exponentials are unchanged.

Independent full-density exponentials also check absolute refocusing before, at, and after its timestamp for both eig and expm generalized propagation. The direct interval path applies its rotation only when that sample includes the pulse.

The absolute-time Hahn reference uses the fixed physical receiver: only samples containing the observed-spin inversion interchange the detected coherence. Both the analytic spin-manifold expression and independent full-density propagation verify the complex signal before, at, and after the pulse.
