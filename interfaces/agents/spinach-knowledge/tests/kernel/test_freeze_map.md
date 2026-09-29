# tests/kernel/test_freeze_map.m

Regression test for frozen input derivatives in Spinach optimal-control waveform maps.

Source: [tests/kernel/test_freeze_map.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_freeze_map.m)

## Purpose

`test_freeze_map` verifies that freezing waveform coordinates during optimal-control optimization zeroes exactly the frozen derivatives while leaving the objective value and all unfrozen derivatives unchanged. The test exercises temporal filters, phase cycles, power scaling, and constrained Hessians in both the Hilbert and vectorised Liouville density formalisms, using noncommuting one-spin transfers.

## Behavior

- Creates a regression result via `new_test_result('optimcon/freeze_map', ...)` with the description "Frozen waveform-map derivatives" and the requirement that free derivatives must include every physical waveform coordinate.
- Builds a quiet one-spin `1H` problem with `spin_ops=pauli(2)` and iterates over the formalisms `{'zeeman-liouv','zeeman-hilb'}` and two non-collinear transfer fixtures:
  - Fixture 1: `rho_init=spin_ops.x`, `rho_targ=spin_ops.y+0.3*spin_ops.z`, `waveform=[3 -1 2;2 4 -3]`.
  - Fixture 2: `rho_init=spin_ops.y+0.2*spin_ops.z`, `rho_targ=spin_ops.z-0.4*spin_ops.x`, `waveform=[-4 2 1;1 -3 5]`.
- For the Liouville formalism, constructs left/right superoperators `lx`, `ly`, `lz` via Kronecker products and vectorises the initial and target states; for Hilbert, uses the Pauli operators directly.
- Common control settings: `control.channels=[1;1]`, `control.operators={lx,ly}`, `control.drifts={{0.7*lz}}`, `control.pwr_levels=1.3`, `control.pulse_dt=[0.04 0.05 0.06]`, `control.method='lbfgs'`, `control.max_iter=0`, `control.plotting={}`, `control.penalties={'none'}`, `control.p_weights=0`, `control.l_bound=-100`, `control.u_bound=100`.
- Runs six `grape_xy` variants per fixture/formalism covering identity, temporal (`firf(w,[0.8 0.4])` and `spf(w,0.3)` distortion maps), phase (`phase_cycle=[0 0.7 0]`), combined, empty-mask (freeze cleared to `[]`), and trapezium-integrator (`integrator='trapezium'`, `pulse_dt=[0.04 0.05]`) paths. For each variant it checks:
  - Frozen gradient entries are exactly zero.
  - Objective fidelity equals the unfrozen run (tolerances 0).
  - Unfrozen gradient entries match the free pullback to `1e-12`.
  - Active gradients agree with centred finite differences of the objective at increments `1e-3` and `1e-4` to `1e-8`.
- Runs nine `grape_curv` variants per fixture/formalism with curvilinear coordinate maps `u2x`/`dx_du` (identity, coupled linear, dimension-increasing, rank-reducing, and polar `u(1)*cos/sin(u(2))` maps), including empty masks, trapezium integration, and combined distortion plus phase cycling, with `penalties={'NS'}` and `p_weights=0.2`. It checks that frozen coordinates vanish in every objective channel, free coordinates retain all Cartesian contributions, values are unchanged, and full-objective derivatives match finite differences at increments `1e-4` and `1e-5` to `1e-8`; each finite-difference check prints a diagnostic line with the error and reference norms.
- Checks phase-only optimisation with `grape_phase`: `method='newton'`, `amplitudes=[2 3 4]`, `freeze=[true false false]`, `phase_cycle=[0 0.7 0]`, phases `[0.2 -0.4 0.7]`. Verifies that freezing preserves the fidelity, free gradient, and free Hessian (the frozen first phase's gradient and Hessian row/column are zeroed in the reference) to `1e-12`.
- Checks exact-Hessian methods `{'newton','goodwin','newton'}` (third skipped for the Hilbert formalism) with `phase_cycle=[0 0.7 0]` and mask entries `(1,1)` and `(2,3)` frozen; the third method adds a non-Hermitian drift `0.7*lz-1i*0.3*diag([0 1 1 0])`. Verifies frozen gradient entries, frozen Hessian rows and columns are exactly zero, the objective is unchanged, and active Hessians match finite differences of production gradients at increments `1e-3` and `1e-4` to `1e-8`.
- Checks the direct engines: for Liouville, `grape_liouv` must retain existing mask semantics (frozen entries compared against a free run with `free_grad(mask)=0`, exact match); for Hilbert, `grape_hilb` must retain existing unmasked derivatives (exact match of the full gradient).
- Final steady-state section switches to `sphten-liouv` with `stst_tol=1e-10`, operators `0.2*diag([0 1 -1 0])` and `0.3*diag([0 0 1 -1])`, a non-Hermitian drift with entries including `drift(4,1)=-0.02i`, `rho_init={[1;0;0;0]}`, `rho_targ={[0;0;0;1]}`, `pwr_levels=1`, `pulse_dt=[0.02 0.03 0.03 1e5]`, `method='rbfgs'`, `steady=true`, `budget=1`, and `distortion={@(w)firf(w,[0.8 0.4])}`:
  - With `freeze=[false(2,2) true(2,2)]` and waveform `[1 2 0 0;0 1 0 0]`, asserts that fidelity and gradient are all finite so a physically frozen long delay does not request its derivative.
  - With `freeze=[false(2,2) true(2,1) [true; false]]`, `distortion={@cancel_last}`, and `phase_cycle=[0 pi/4 0]`, asserts finiteness again, verifying that composed Jacobians identify cancelled physical derivatives.
- The local helper `cancel_last(waveform)` returns the waveform with the last column replaced by `[sum(waveform(:,end));0]` together with a sparse `[2*nsteps x 2*nsteps]` Jacobian built from `speye` with the two last-column rows zeroed and `J(rows(1),rows)=1`.

## Inputs and outputs

```matlab
result = test_freeze_map()
```

**Inputs:** none.

**Outputs:**

- `result` — regression test result structure with explanatory messages, accumulated through `test_close` and `test_true` assertions.

## References

- Spinach optimal control module functions used here: `optimcon`, `grape_xy`, `grape_curv`, `grape_phase`, `grape_liouv`, `grape_hilb`, `new_test_result`, `test_close`, `test_true`, `firf`, `spf`, `pauli`.
