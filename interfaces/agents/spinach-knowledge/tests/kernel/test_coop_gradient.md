# tests/kernel/test_coop_gradient.m

Source: [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_coop_gradient.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_coop_gradient.m)

## Purpose

Regression test for cooperative phase gradients and the initial gradient guard in Spinach's optimal control module. It verifies that the primary fidelity and the squared impurity penalty share one gradient, that the cooperative objective matches independent matrix propagation, and that the optimiser's assembled-initial-guess safeguards behave correctly.

## Behaviour

- Creates a regression result via `new_test_result('kernel/coop_gradient','Cooperative phase gradients','The primary fidelity and squared impurity must share one gradient.')`.
- Sets a small-system numerical environment: `spin_system.tols.liouv_zero=1e-14`, `spin_system.tols.small_matrix=64`, `spin_system.tols.dense_matrix=0.5`, `spin_system.tols.prop_chop=1e-14`, isotopes `{'1H'}`, and Pauli operators from `pauli(2)`.
- Sweeps four fixtures combining two formalisms (`'zeeman-liouv'`, `'zeeman-hilb'`) and two fidelity measures (`'real'`, `'square'`):
  - Fixture 1: unit target `spin_ops.x+0.3*spin_ops.z` normalised by Frobenius norm, single power level.
  - Fixture 2: nonunit complex detection operator `1.7*(spin_ops.z-0.4*spin_ops.x+1i*spin_ops.y)`, power ensemble `[0.8 1.1]`.
  - Fixture 3: rectangular initial state `[0 1;0 0]`, complex target `[1 3;1+1i 0]`, phases including `pi`, zero drift.
  - Fixture 4: identity initial and target states with zero amplitudes and zero drift.
- For each fixture, independently propagates both pulses in Hilbert space using `expm` and computes overlaps and squared impurity costs; the expected cooperative objective is `mean(primary,'all')-mean(dirt_cost)`.
- Calls `grape_coop(phase_pair,local_system)` and checks:
  - Objective value against the independent propagation with tolerance `1e-11`.
  - That both pulse trajectories are returned (`numel(traj_data)==2` with cell contents).
  - Phase gradient against centred finite differences at increments `steps=[1e-3 1e-4 1e-5]`, with tolerance `2e-9` plus `2*step_size^2` scaling.
- Fixture 3 additionally checks:
  - Purely imaginary auxiliary overlaps (`real(overlap)==0 && imag(overlap)~=0`).
  - Zero real auxiliary fidelity returning a nonzero derivative, verified against independent matrix exponentials with tolerance `2*step_size^2+2e-9`.
  - Zero impurity costate giving exact zero value and gradient.
  - Constant overlap (`identity` to `identity`) giving fidelity `2` with exact zero gradient.
  - For each optimiser method in `{'lbfgs','rbfgs','newton','goodwin'}`:
    - Zero-value admission: a nonzero derivative is admitted even when the value vanishes.
    - Zero-gradient guard: a stationary guess fails with `'gradient too small at iter 1, find a better guess.'`.
    - Objective-only evaluation with `max_iter=0` returns the unchanged point with `fx==1`, `gfx==0`, `hfx==0`.
    - A nonstationary physical transfer (`spin_ops.x` to `spin_ops.x+0.3*spin_ops.z`, guess `[0.2 -0.3;0.4 0.5]`, `max_iter=3`) improves fidelity; only `'newton'` and `'goodwin'` request Hessians, once per iteration.
    - Frozen coordinates (`freeze=true(size(guess))`) trigger the same gradient-too-small error.
- After the fixture loops, checks a nonstationary cooperative objective at a zero-score cancellation using `coop_guess=0.364110104613*ones(2,2)` with `grape_phase` primary scores; verifies `abs(coop_before(1))<1e-6`, `primary>1e-3`, and gradient norm `>1e-3`, then confirms `fmaxnewton` admits and improves it in one iteration.
- Repeats the cancellation admission through an anonymous forwarding adapter `@(wave,system)grape_coop(wave,system)`, requiring identical function-evaluation counts.
- Admits a nonstationary impurity direction with no primary transfer (target changed to `eye(2)`, guess perturbed by `pi/4`), verifying negative initial score, nonzero gradient, admission, and one-iteration improvement.
- Rejects a frozen cooperative gradient when all phases are frozen, expecting the same `'gradient too small at iter 1, find a better guess.'` error.
- Prints diagnostic lines prefixed `COOP`, `AUX`, and `OPTIMISER` with step sizes, errors, scales, and iteration counts.

## Inputs and outputs

**Signature:** `result=test_coop_gradient()`

- **Inputs:** none.
- **Outputs:**
  - `result` — regression result structure with explanatory messages, populated by `test_close` and `test_true` assertions.

## References

- Uses Spinach functions: `new_test_result`, `pauli`, `optimcon`, `grape_coop`, `grape_phase`, `grape_liouv`, `grape_hilb`, `grape_xy`, `fmaxnewton`, `hdot`, `test_close`, `test_true`.
