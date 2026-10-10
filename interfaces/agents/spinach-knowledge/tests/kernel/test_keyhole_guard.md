# tests/kernel/test_keyhole_guard.m

[View source on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_keyhole_guard.m)

## Purpose

Regression test for the supported boundaries of keyhole optimal control in Spinach. It verifies that spin-state projectors acting between noncommuting spin-half pulse slices are handled correctly: density-vector and wavefunction bases must reject Newton/Goodwin keyhole methods at both setup and direct-engine entry while retaining first derivatives, and empty schedules and Hilbert-space keyhole Hessians must remain supported.

## Behaviour

- Declares a regression named `kernel/keyhole_guard` ("State-vector keyhole Hessian guard") with the requirement that unsupported Hessians are refused without changing supported methods.
- Builds a spin-half system with `spin_system.sys.output='hush'`, empty enable/disable lists, tolerances (`liouv_zero=1e-14`, `small_matrix=64`, `dense_matrix=0.5`, `prop_chop=1e-14`), and isotope `1H`.
- Uses Pauli operators `pauli(2)` to define initial and target states `rho_init=spin_ops.y+0.2*spin_ops.z` and `rho_targ=spin_ops.z-0.4*spin_ops.x`, a control waveform `[-4 2 1;1 -3 5]`, methods `{'lbfgs','rbfgs','newton','goodwin'}`, and formalisms `{'zeeman-liouv','sphten-liouv','zeeman-hilb','zeeman-wavef'}`.
- Transforms into normalised spherical tensors (identity, T11, T10, T1-1) using the explicit coordinate matrix `Q=[1/sqrt(2) 0 1/sqrt(2) 0; 0 0 0 1; 0 -1 0 0; 1/sqrt(2) 0 -1/sqrt(2) 0]` for the `sphten-liouv` formalism.
- For each formalism, constructs a minimal control problem with control operators, drift `0.7*spin_ops.z` (Liouville-transformed where applicable), source/target states, a keyhole projector (`diag(diag(rho))` for Hilbert, `P*rho` with `P=diag([1 0])` for wavefunction, `P=diag([1 0 0 1])` for density-vector), `pwr_levels=1`, `pulse_dt=[0.04 0.05 0.06]`, `max_iter=0`, no penalties, bounds `[-100,100]`, and no plotting.
- For Newton/Goodwin on non-Hilbert formalisms, requires `optimcon` setup to throw an error containing both `keyholes with Newton/Goodwin Hessians` and `not implemented`; then bypasses setup validation by setting `method='lbfgs'` before `optimcon` and restoring the method, requiring `grape_liouv` direct calls to throw the same error.
- Verifies that `optimcon` never substitutes another optimisation algorithm (`method retained` check for every method/formalism combination).
- For first-order methods (`lbfgs`, `rbfgs`) on non-Hilbert formalisms, requires both `grape_xy` and direct `grape_liouv` calls to refuse requested keyhole Hessians with errors mentioning `keyholes Hessians` / `not implemented`.
- Keeps empty keyhole schedules (`cell(1,3)`) available for both exact-Hessian methods and verifies their gradients and Hessians.
- Independently validates gradients (and Hessians for `method_idx>2`) by central finite differences with `step_size=1e-4` and tolerance `1e-9`, comparing against `grape_xy` outputs.
- For Liouville formalisms, additionally tests that dissipative Newton keyholes (drift modified by `-1i*0.3*(eye(4)-P)`) are refused at setup with the same unsupported-combination error.

## Inputs and outputs

```matlab
result = test_keyhole_guard()
```

- **Output**: `result` — regression test results and explanatory messages accumulated through `new_test_result`, `test_true`, and `test_close`.
- **Input**: none.

## References

- Spinach optimal control module functions used: `optimcon`, `grape_xy`, `grape_liouv`, `hilb2liouv`, `pauli`.
