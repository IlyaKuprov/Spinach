# tests/kernel/test_cce_exact.m

[MATLAB source](../../../../../tests/kernel/test_cce_exact.m)

## Purpose

Four bath spins, pure/mixed states, FID/Hahn, and raw subsets.

## Behaviour

Checks the following named assertions in the source:

- `raw subset exact`.
- `full order bypass`.

The test uses `new_test_result`, `test_true` and/or `test_close`; numerical tolerances and independent reference calculations are defined in the linked source.

## Inputs and outputs

```matlab
result=test_cce_exact()
```

**Inputs**

- None; run with Spinach and its test libraries on the MATLAB path.

**Outputs**

- `result` — standard test-harness result with checks, pass/fail status, and diagnostics. Printed numerical metrics are not additional result fields.

## References

- [Source file](../../../../../tests/kernel/test_cce_exact.m)

## Native interface

Uses `clexpn` with `@free_evolution` or `@hahn_echo`. `parameters.transition.spins` states physical transition support. Natural numeric output and separate raw-cluster diagnostics are checked against independent full-density calculations with an explicit fixed transition detector. Full-order bypass is checked by exact equality to the full-cluster raw response.

The dense reference embeds both independently assembled spin-manifold propagators into the full observed-spin/bath Hilbert space, applies the physical observed-spin inversion, propagates the initial density, and detects with an unchanged transition operator. Pure complex and maximally mixed bath preparations cover every retained subset; FID and zero-time controls remain included.
