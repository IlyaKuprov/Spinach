# tests/kernel/test_cce_partition.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_partition.m)

## Purpose

Dense exact unions, pair factors and bounded adaptation.

## Behaviour

Checks the following named assertions in the source:

- `Singleton pCCE equals ordinary CCE`.
- `One group equals dense exact echo`.
- `Nontrivial dense echo`.
- `pCCE group-order two equals dense exact`.
- `Disjoint-inter equals dense exact`.
- `Disjoint equals dense group product`.
- `Cross-group correction is nonzero`.
- `Three-group pair-factor identity`.
- `Canonical physical-ID groups`.
- `Empty groups mean singleton vertices`.
- `Deterministic strongest-edge grouping`.
- `Deterministic spatial grouping`.
- `Adaptive pilot and held-out gates`.
- `Probe ceiling is incomplete`.
- `Exact patch dimension ceiling is incomplete`.
- `Round ceiling is incomplete`.
- `Held-out probe cannot be skipped`.
- `Held-out numerical discrepancy is incomplete`.
- `Overlapping groups rejected`.
- `Missing physical IDs rejected`.
- `Nonfinite capacity rejected`.
- `Fractional probe ceiling rejected`.
- `Nonfinite patch tolerance rejected`.
- `Empty held-out contexts rejected`.
- `Fractional spatial ceiling rejected`.
- `Missing spatial coordinate rejected`.
- `One-group generalized density equals dense echo`.
- `Master-equation adaptive evaluator`.

The test uses `new_test_result`, `test_true` and/or `test_close`; numerical tolerances and independent reference calculations are defined in the linked source.

## Inputs and outputs

```matlab
result=test_cce_partition()
```

**Inputs**

- None; run with Spinach and its test libraries on the MATLAB path.

**Outputs**

- `result` — standard test-harness result with checks, pass/fail status and diagnostics. Printed numerical metrics are not additional result fields.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_partition.m)
- [Test manifest](https://github.com/IlyaKuprov/Spinach/blob/main/tests/lib/test_manifest.m)

Public optional-field contracts are explicit in reused fixtures: irrelevant method, event, geometry, and legacy numerical options are removed rather than passed silently.
