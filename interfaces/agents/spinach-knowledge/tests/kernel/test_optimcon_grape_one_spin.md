# tests/kernel/test_optimcon_grape_one_spin.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_optimcon_grape_one_spin.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_optimcon_grape_one_spin.m)

## Purpose

Regression test for the one-spin optimal-control setup and Hilbert-space GRAPE in Spinach. The test checks that `optimcon()` accepts a minimal one-spin Hilbert-space control problem, and that `grape_hilb()` fidelity and gradient agree with independent matrix exponentiation and finite differences.

## Behaviour

- Announces the test target with `fprintf` and initialises a regression test result via `new_test_result('optimcon/grape_one_spin', ...)`.
- Ensures a parallel pool exists for the ensemble loop: if `gcp('nocreate')` returns empty, it calls `parpool('Processes',1)`.
- Builds a minimal Hilbert-space Spinach object through a local helper `local_spin_system()`, which sets `sys.output='hush'`, empty `sys.enable`/`sys.disable`, tolerances (`liouv_zero=1e-14`, `small_matrix=64`, `dense_matrix=0.5`, `prop_chop=1e-14`), `bas.formalism='zeeman-hilb'`, and `comp.isotopes={'E'}`.
- Defines one-spin Pauli operators via `pauli(2)`, a sparse 2-by-2 zero drift, and a two-step timing grid `pulse_dt=[0.02 0.03]`.
- Configures a minimal control problem with `control.isotopes={'E'}`, `control.channels=1`, `control.operators={S.y}`, `control.rho_init={S.x}`, `control.rho_targ={S.z}`, `control.pwr_levels=1`, `control.pulse_dt=pulse_dt`, `control.drifts={{drift}}`, `control.method='lbfgs'`, `control.max_iter=0`, `control.penalties={'none'}`, `control.p_weights=0`, `control.l_bound=-100`, `control.u_bound=100`, and empty plotting, then passes it to `optimcon()`.
- Verifies `optimcon()` absorbed the control metadata: one registered control (`spin_system.control.ncontrols==1`), the timing grid stored without modification (`test_close` with zero tolerances), and the rectangle grid using one waveform value per interval (`spin_system.control.pulse_ntpts==numel(pulse_dt)`).
- Re-runs `optimcon()` with `control.method` removed and `control.distortion={@no_dist}` supplied, checking that the method defaults to `'lbfgs'` before distortion processing and that the supplied distortion function is stored (`numel(spin_default.control.distortion)==1`).
- Checks that exact-Hessian methods remain incompatible with distortions: with `control.method='newton'` and `control.distortion={@no_dist}`, `optimcon()` is called inside a `try`/`catch`, and the test asserts the caught error message contains `'waveform distortions'`.
- Evaluates GRAPE fidelity and gradient for the non-trivial waveform `[7 11]` via `grape_hilb(spin_system,{drift},control.operators,waveform,S.x,S.z,'real')`.
- Builds an independent exact Hilbert-space trajectory by propagating `rho=S.x` with `P=expm(-1i*H*pulse_dt(n))` where `H=waveform(n)*S.y`, computing `fid_ref=real(hdot(rho,S.z))`, and compares the GRAPE fidelity with absolute and relative tolerances of `1e-13`.
- Computes a centred finite-difference gradient with step size `1e-6`, perturbing each waveform element by plus/minus the step and calling `grape_hilb()` for each perturbed waveform, then compares the analytic gradient with absolute and relative tolerances of `1e-8`.
- Checks that no heavy trajectory is returned when plotting is disabled: `isempty(traj_data.forward)` must hold.

## Inputs and outputs

**Syntax:**

```matlab
result = test_optimcon_grape_one_spin()
```

**Outputs:**

- `result` — regression test result with explanatory messages, accumulated through `test_true` and `test_close` assertions covering the `optimcon()` control metadata, distortion defaulting and Hessian rejection, `grape_hilb()` fidelity and gradient agreement, and trajectory suppression.

The function takes no inputs.

## References

- [Spinach GitHub repository — test_optimcon_grape_one_spin.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_optimcon_grape_one_spin.m)
- [Spinach GitHub repository (main page)](https://github.com/IlyaKuprov/Spinach)
