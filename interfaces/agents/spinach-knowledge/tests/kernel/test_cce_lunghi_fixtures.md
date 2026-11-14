# tests/kernel/test_cce_lunghi_fixtures.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_lunghi_fixtures.m)

## Purpose

Reduced VOTPP and Cu(mnt)2 source-input gCCE reconstruction.

## Behaviour

Requires `SPINACH_CCE_LUNGHI_FIXTURES` (external inputs/CSVs) and
`SPINACH_CCE_LUNGHI_STATES` (phase-preserving states JSON and diagnostic MAT
files). Without both, the test returns with an explanatory message and does
not run the external comparison. These unlicensed source structures are not
redistributed; reduced source-input validation is not published-figure validation.

Runs the assertions and tolerance checks declared in the source.

The test uses `new_test_result`, `test_true` and/or `test_close`; numerical tolerances and independent reference calculations are defined in the linked source.

## Inputs and outputs

```matlab
result=test_cce_lunghi_fixtures()
```

**Inputs**

- None; run with Spinach and its test libraries on the MATLAB path.

**Outputs**

- `result` — standard test-harness result with checks, pass/fail status and diagnostics. Printed numerical metrics are not additional result fields.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_lunghi_fixtures.m)
- [Test manifest](https://github.com/IlyaKuprov/Spinach/blob/main/tests/lib/test_manifest.m)
