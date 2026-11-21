# tests/kernel/test_cce_convergence.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_convergence.m)

## Purpose

Independent order/time ladders, resource exhaustion, and held-out checks.

## Behaviour

Checks the following named assertions in the source:

- `independent ladder success`.
- `independent work accounting`.
- `unchanged physics`.
- `resource cap`.
- `missing times refusal`.
- `held-out sampling requirement`.
- `held-out check execution`.
- `invalid coverage cannot converge`.

The test uses `new_test_result`, `test_true` and/or `test_close`; numerical tolerances and independent reference calculations are defined in the linked source.

## Inputs and outputs

```matlab
result=test_cce_convergence()
```

**Inputs**

- None; run with Spinach and its test libraries on the MATLAB path.

**Outputs**

- `result` — standard test-harness result with checks, pass/fail status, and diagnostics. Printed numerical metrics are not additional result fields.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_convergence.m)
- [Test manifest](https://github.com/IlyaKuprov/Spinach/blob/main/tests/lib/test_manifest.m)

## Native contract

Native convergence uses an explicit registered experiment and assumptions in `[answer,diagnostics]=cce_converge(spin_system,parameters,pulse_sequence,assumptions)`. Its independent native-path ladders preserve the original baseline, work cap, held-out seed requirement, shared-time coverage, and numerical validity. Agreement between later refinements never certifies an inaccurate baseline.

Four native-boundary checks reject NaN, infinity, fractional values, and text in an order ladder before any baseline evaluation. The original finite-bath order-three refinement remains accepted.
