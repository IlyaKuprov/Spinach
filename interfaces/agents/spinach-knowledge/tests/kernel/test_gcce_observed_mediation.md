# tests/kernel/test_gcce_observed_mediation.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_gcce_observed_mediation.m)

## Purpose

Two bath spins with zero direct coupling share a flip-flopping centre.

## Native interface

Uses numeric `clexpn(...,@free_evolution,parameters,'labframe')` outputs and the separate diagnostics structure. Physical observed IDs are in `parameters.observed`; graph, order, bath-state and averaging controls are in `parameters.cce`. The exact three-spin reference, nonfactorisation requirement, and joint-detection cases are retained.

## Behaviour

Checks the following named assertions in the source:

- `pair retained`.
- `exact pair`.
- `nonfactorization measured`.

The test uses `new_test_result`, `test_true` and/or `test_close`; numerical tolerances and independent reference calculations are defined in the linked source.

## Inputs and outputs

```matlab
result=test_gcce_observed_mediation()
```

**Inputs**

- None; run with Spinach and its test libraries on the MATLAB path.

**Outputs**

- `result` — standard test-harness result with checks, pass/fail status, and diagnostics. Printed numerical metrics are not additional result fields.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_gcce_observed_mediation.m)
- [Test manifest](https://github.com/IlyaKuprov/Spinach/blob/main/tests/lib/test_manifest.m)

The unitary and master-equation checks explicitly select `parameters.cce.basis.formalism` as `zeeman-hilb` and `zeeman-liouv`, respectively. An empty jump list does not itself select a representation. The zero-rate Liouville case retains the original joint-observable response and arithmetic-error bound, including internal bath averaging.
