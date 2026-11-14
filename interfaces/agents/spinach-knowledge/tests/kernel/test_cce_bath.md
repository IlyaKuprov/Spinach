# tests/kernel/test_cce_bath.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_bath.m)

## Purpose

Nonorthogonal cells, occupancy, density and covariance.

## Behaviour

Checks the following named assertions in the source:

- `Nonorthogonal row-vector cell convention`.
- `Unique physical site IDs`.
- `Seeded isotope occupancy is reproducible and local`.
- `Coordinate covariance`.
- `Tensor covariance`.
- `Surface density gives twenty unique sites`.
- `SI density and area conversion`.
- `Defect removal preserves physical group IDs`.
- `Duplicate physical sites are refused`.
- `Complete groups and isotope row convention`.
- `Scalar isotropic tensor`.
- `Native dipolar transform covariance`.
- `Unknown group site refused`.
- `Unknown pair site refused`.
- `Linear density and frozen coordinates`.
- `Packing ceiling exposed`.

The test uses `new_test_result`, `test_true` and/or `test_close`; numerical tolerances and independent reference calculations are defined in the linked source.

## Inputs and outputs

```matlab
result=test_cce_bath()
```

**Inputs**

- None; run with Spinach and its test libraries on the MATLAB path.

**Outputs**

- `result` — standard test-harness result with checks, pass/fail status and diagnostics. Printed numerical metrics are not additional result fields.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_bath.m)
- [Test manifest](https://github.com/IlyaKuprov/Spinach/blob/main/tests/lib/test_manifest.m)
