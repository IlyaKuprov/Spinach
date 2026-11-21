# tests/kernel/test_cce_controls.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_controls.m)

## Purpose

Check finite drives, waveform refinement, and named pulse schedules.

## Behaviour

Checks the following named assertions in the source:

- `Finite detuned drive matches independent dense expm`.
- `Public finite-control driver`.
- `Finite pulse approaches ideal rotation`.
- `Piecewise midpoint refinement is second order`.
- `XY8 timestamps`.
- `XY8 phases`.
- `UDD timestamps`.
- `Exterior ideal pulse flips the sampled field`.
- `Exterior pulse changes subsequent Hamiltonian`.
- `Finite exterior drive is refused rather than frozen`.
- `Projected finite observed drive refused`.
- `Selective observed-spin system operator ordering`.
- `Physical spin-one pi pulse`.
- `Selective two-level pi pulse`.
- `Physical and selective pulses differ`.

The test uses `new_test_result`, `test_true` and/or `test_close`; numerical tolerances and independent reference calculations are defined in the linked source.

## Inputs and outputs

```matlab
result=test_cce_controls()
```

**Inputs**

- None; run with Spinach and its test libraries on the MATLAB path.

**Outputs**

- `result` — standard test-harness result with checks, pass/fail status, and diagnostics. Printed numerical metrics are not additional result fields.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_controls.m)

Public optional-field contracts are explicit in reused fixtures: irrelevant method, event, geometry, and legacy numerical options are removed rather than passed silently.

A supplied observed density and polarised bath state are checked against an independent dense finite-pulse propagator, including normalisation by the experiment-owned initial observable. The native subsystem rejects unsupported bath-configuration reference overrides.

The registered identifier is unchanged. All original physical checks run through
clexpn, finite_drive, hahn_echo, refocusing_train, or native compiled subsystem
services. Former omitted/empty operator-angle spellings are tested as equivalent
integrated-generator and Cartesian declarations; a physical angle is never
implicitly applied twice. Raw signal checks explicitly request normalized=false.
Named Hahn direct-segment arithmetic is selected with cce.numerics.propagator=expm.
The test retains dense detuned-drive, ideal-width, midpoint-waveform, tensor-order,
spin-one/selective-pulse, exterior-field, memory, and scale-invariant Hermiticity
assertions, without calling the retired public control/evaluator interface.
