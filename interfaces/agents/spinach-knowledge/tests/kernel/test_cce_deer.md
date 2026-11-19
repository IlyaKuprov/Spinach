# tests/kernel/test_cce_deer.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_deer.m)

## Purpose

Check strict frequency windows, isotope filters, and empty selections.

## Behaviour

Checks the following named assertions in the source:

- `Both exact boundaries excluded and species respected`.
- `Boundary rule is recorded`.
- `Explicit site selection intersects frequency window`.
- `Empty DEER window requires explicit opt-in`.
- `Explicit empty-selection opt-in accepted`.
- `Canonical frequencies are not converted twice`.

The test uses `new_test_result`, `test_true` and/or `test_close`; numerical tolerances and independent reference calculations are defined in the linked source.

## Inputs and outputs

```matlab
result=test_cce_deer()
```

**Inputs**

- None; run with Spinach and its test libraries on the MATLAB path.

**Outputs**

- `result` — standard test-harness result with checks, pass/fail status, and diagnostics. Printed numerical metrics are not additional result fields.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_deer.m)
- [Test manifest](https://github.com/IlyaKuprov/Spinach/blob/main/tests/lib/test_manifest.m)

Public optional-field contracts are explicit in reused fixtures: irrelevant method, event, geometry, and legacy numerical options are removed rather than passed silently.
The registered test also checks seven unsupported-controls entry points, ignored nested options, projected second-order applicability, unsupported Lindblad noise and noise weights, explicit raw suppression, and driver postprocessing parity.

Finite-radius compilation regressions cover missing/complete pre-spin_config geometry, named sequences on all four signal methods, enabled DEER compilation, generalised finite drives and bath-targeted pulses. Propagation still requires observed/exterior coordinates; malformed radii remain rejected before a spin configuration exists. Zero-field uncoupled controls retain analytically unit signal.

Native validation uses clexpn and experiment-owned pump selection, with numeric answers and separate diagnostics. Strict axial frequency filtering retains physical site and isotope selection, explicit empty-window consent, and canonical rad/s units. Single bath rotations use hahn_echo refocusing descriptors after the selector; bath_pumped_hahn itself retains the observed Hahn rotation. All original physical checks remain, including the finite-radius geometry prerequisite across direct, scheduled, and interlaced evaluation. Additional positive controls verify a nonzero exterior field crosses the radius boundary and that a Cartesian exterior inversion reverses the subsequent field. Obsolete controls are tested at native ownership boundaries; postprocessing uses explicit analysis and spectral services.
