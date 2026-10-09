# tests/kernel/test_dynamic_optimcon_remaining.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_optimcon_remaining.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_optimcon_remaining.m)

## Purpose

Regression test for the remaining dynamic optimal-control helper functions in Spinach. It exercises waveform distortion models, FIR kernel estimation, quasi-Newton updates, Hessian handling, waveform utilities, GRAPE wrappers, Liouville-space GRAPE derivatives, TGRAPE duration gradients, `fmaxnewton` zero-iteration handling, and diagnostic plotting smoke paths, using small deterministic fixtures.

## Behaviour

The test function takes no arguments and returns a regression test result object with explanatory messages. It announces the target with `fprintf('TESTING: Remaining optimal-control dynamic helpers\n')`, creates a result via `new_test_result('kernel/dynamic_optimcon_remaining', ...)`, ensures a parallel pool exists (starting a one-worker `parpool('Processes',1)` if none is open), and then runs four independent check groups:

1. **Distortions** (`local_check_distortions`):
   - `no_dist`: identity distortion on a 2×4 waveform and its vectorised identity Jacobian (`speye(numel(waveform))`).
   - `non_orth`: non-orthogonal channel mixing at `xy_ang = 60` degrees, checked against the closed form `[waveform(1,:)+cosd(xy_ang)*waveform(2,:); sind(xy_ang)*waveform(2,:)]`, plus a Jacobian check.
   - `firf`: FIR filtering with kernel `[1.00; 0.25+0.10i; -0.15]`, compared to explicit truncated causal complex convolution per XY channel pair.
   - `spf`: single-pole filter with pole `0.35+0.15i`, compared to the causal recurrence `filtered(k) = (1-pole)*signal(k) + pole*filtered(k-1)`.
   - `szf`: single-zero filter with zero `0.20+0.10i`, compared to the causal recurrence `filtered(k) = (signal(k) - zero*signal(k-1))/(1-zero)`.
   - `amp_tanh`: hyperbolic-tangent amplifier compression at saturation level 2.0; phase is preserved, amplitudes stay below the saturation level.
   - `amp_root`: root-sigmoidal compression at saturation level 2.0 with shape 4; radial amplitudes are not increased.
   - `kernelest`: causal FIR kernel recovery of `h_ref = [0.60; -0.20; 0.10]` from `x = [1; 0; -1; 2; 1; -2; 3]` using the `backslash`, `pinv`, and `svd` methods (tolerance 1e-12), a Tikhonov-regularised call with regularisation parameter 1e-12 (tolerance 1e-8), and a `same`-alignment estimation using central Toeplitz convolution rows.
   - All Jacobians are verified by comparing `J*probe` against centred finite differences with a normalised sine probe and step size 1e-6 (tolerances 1e-9, or 1e-8 for the amplifier Jacobians).

2. **Quasi-Newton** (`local_check_quasi_newton`):
   - `bfgs_upd`: a good first curvature pair `([1;0], [-2;0])` initialises the Hessian to `2*eye(2)`; a bad first pair `([1;0], [2;0])` falls back to identity; a later bad pair leaves a `3*eye(2)` Hessian unchanged.
   - `bfgs`: orthogonal curvature histories reconstruct the diagonal Hessian `diag([2 3])`.
   - `lbfgs`: applies the inverse of the reconstructed Hessian, giving direction `[1; 1]`.
   - `hess_reorder`: swaps control-first and time-first ordering, verified against an explicit tensor permutation of `reshape(1:36,6,6)` with dimensions 2 and 3.
   - `hessreg`: RFO regularisation of the indefinite Hessian `[-1 0; 0 2]` with `reg_alpha=1`, `reg_phi=2`, `reg_max_iter=8`, `reg_max_cond=1e4`; the result must be positive definite (Cholesky flag 0), the RFO iteration counter must increase, and all elements must be finite.

3. **Wave utilities** (`local_check_wave_utils`):
   - `fapt2sfo`: frequency-amplitude-phase-time conversion of two events on the grid `0:0.25:1`, checked against explicit summation of masked rotating components; `dt` remains empty when a grid is supplied and the grid passes through unchanged.
   - `inst_freq`: instantaneous frequency of a quadratic-phase chirp (base frequency 3, chirp rate 4, sample dt 1e-3, stencil order 5, polynomial order 2) recovers `base_freq + chirp_rate*time_axis` to 1e-10; low-amplitude points (tolerance 0.5) and zero-magnitude points (tolerance 0) mask every stencil containing them (indices 1:4 for the fixture), leaving unaffected stencil estimates intact.
   - `drifts`: ensemble drift extraction through a context callback returning two members; the count is preserved, the first member equals `H + 3i*H` (Hamiltonian, relaxation, kinetics) and the second `2*H + 6i*H` (including hydrodynamics), with `H = sparse(diag(1:4))`; the classical subspace dimension is 2.
   - `aux_mat`: trapezium-product auxiliary matrices for a two-control Pauli fixture with time step 1e-3; the left block diagonal matches the explicit trapezium generator, the derivative blocks match explicit left/right directional derivatives, and the mixed-derivative call returns 6×6 blocks.
   - `ctrl_trajan`: an `xy_controls` plotting smoke path runs offscreen, and an `frq_controls` plot recovers the constant instantaneous frequency 123.4 of a linear-phase signal (65 points, time step 1e-4) to 1e-10.

4. **GRAPE family** (`local_check_grape_family`):
   - `ensemble` and `grape_xy`: the Cartesian wrapper reports the ensemble fidelity and gradient in its first channel and appends one gradient channel per penalty term (`size(grad_xy,3)==2`).
   - `grape_curv`: identity curvilinear coordinates (`local_u2x`, `local_dx_du`) preserve fidelity and gradient channels.
   - `grape_phase`: phase-parameterised waveform is equivalent to its Cartesian form, and its gradient matches centred finite differences (step 1e-6, tolerance 1e-6).
   - `fmaxnewton`: with `max_iter = 0`, returns the supplied point unchanged with zero iterations and zero gradient/Hessian function calls.
   - `grape_liouv`: Liouville-space GRAPE on a 2-dimensional vector-space fixture returns finite fidelity, gradient, and Hessian; the gradient matches centred finite differences (tolerance 1e-6) and the Hessian is symmetric to 1e-10 for a real fidelity.
   - `tgrape`: duration gradients for `dt_grid = [0.03; 0.04]` are finite and match centred finite differences (tolerance 1e-6).
   - `grape_coop`: cooperative phase wrapper on a spherical-tensor-like fixture returns finite fidelity and gradient channels with gradient size `[2 2 2]`, one phase-gradient row per cooperative pulse.

The test uses `test_close` for numerical comparisons and `test_true` for logical assertions, accumulating messages into the returned result. Local fixtures build minimal Spinach systems via `optimcon` with formalisms such as `zeeman-hilb`, `zeeman-liouv`, and `sphten-liouv`, using the electron isotope `'E'` and Pauli operators.

## Inputs and outputs

**Inputs:** none.

**Outputs:**

- `result` — regression test result object with explanatory messages, created by `new_test_result` and progressively updated by the check functions.

## References

- Tested functions: `no_dist`, `non_orth`, `firf`, `spf`, `szf`, `amp_tanh`, `amp_root`, `kernelest`, `bfgs_upd`, `bfgs`, `lbfgs`, `hess_reorder`, `hessreg`, `fapt2sfo`, `inst_freq`, `drifts`, `aux_mat`, `ctrl_trajan`, `ensemble`, `grape_xy`, `grape_curv`, `grape_phase`, `fmaxnewton`, `grape_liouv`, `tgrape`, `grape_coop`, `optimcon`, `pauli`.
- Test infrastructure: `new_test_result`, `test_close`, `test_true`.
- Source file: [tests/kernel/test_dynamic_optimcon_remaining.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_optimcon_remaining.m) in the Spinach repository.

The synthetic compiled fixtures use per-substance descriptor cells and offsets.
