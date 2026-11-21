# tests/kernel/test_cce_io.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_io.m)

## Purpose

Checks the external CCE input boundary: coordinate and frequency unit conversion, active XYZ rotation, scalar and Cartesian tensor conventions, and cube voxel ordering. The spin-CSV test checks that a neighbouring JSON file does not override the table or its declared units. Missing frame and field gyromagnetic declarations, and unsupported MAT/JSON result formats, are rejected.

## Inputs and outputs

```matlab
result=test_cce_io()
```

**Inputs**

- None; run with Spinach and its test libraries on the MATLAB path.

**Outputs**

- `result` — standard test-harness result with checks, pass/fail status, and diagnostics. Printed numerical metrics are not additional result fields.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_io.m)
